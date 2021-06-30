#!/usr/local/bin/perl -w

################################################################################
# runMassChroQ.pl              1.9.15                                          #
# Authors: P. Poullet, V. Laigle, G. Arras, F. Yvon (Institut Curie)           #
# Contact: myproms@curie.fr                                                    #
# Interface between myProMS and MassChroQ software, developed by Inra.         #
# Website : http://pappso.inra.fr/bioinfo/masschroq                            #
# Reference : Valot B., Langella O., Nano E., Zivy M. MassChroQ: A versatile   #
# tool for mass spectrometry quantification. Proteomics (2011), 11, 3572-77    #
################################################################################
#----------------------------------CeCILL License-------------------------------
# This file is part of myProMS
#
# Copyright Institut Curie 2018
#
# This software is a computer program whose purpose is to process
# Mass Spectrometry-based proteomic data.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software. You can use,
# modify and/or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#-------------------------------------------------------------------------------

$| = 1;

use strict;
use POSIX qw(strftime :sys_wait_h);
use promsConfig;
use promsMod;
use File::Copy;
use XML::SAX;
# exit; # DEBUG!!!!!!!!!!!!!!!!!!!

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo('no_user');
my %clusterInfo=&promsConfig::getClusterInfo;
my $masschroqPath = $promsPath{'masschroq'}; # default
my $pathR = $promsPath{'R'}; # default
if ($clusterInfo{'on'}) {
	$masschroqPath = $clusterInfo{'path'}{'masschroq'};
	$pathR = $clusterInfo{'path'}{'R'};
}
my @aaList=('A','C'..'I','K'..'N','P'..'T','V','W','Y');

####################
####>Parameters<####
####################
my ($quantifID,$quantifDate)=@ARGV;
$quantifID=&promsMod::cleanNumericalParameters($quantifID);
my $quantifDir="$promsPath{tmp}/quantification/$quantifDate";
my $errorFile="$promsPath{tmp}/quantification/current/$quantifID\_$quantifDate\_error.txt";
my $maxStates=50;

############################
####> Connecting to DB <####
############################
my $dbh=&promsConfig::dbConnect('no_user');

####>Updating quantification status to 0 (running)<####
$dbh->do("UPDATE QUANTIFICATION SET STATUS=0,UPDATE_DATE=NOW() WHERE ID_QUANTIFICATION=$quantifID");
$dbh->commit;

#############################
###>Generating data files<###
#############################
my $fileStat="$quantifDir/status_$quantifID.out";
open(FILESTAT,">$fileStat");
print FILESTAT "Started ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";

####>Getting the parameters from quantif_info.txt<####
open (INFO,"$quantifDir/quantif_info.txt");
my $section='';
my (%params,%sampIDs,%pepInfo,%pep2Ana,%vmods,%spectralCount,%seqMissedCuts,%proteinScore,%channels,%isoLists,%isoInfo); #,%anaLabelModifs
while (<INFO>) {
	if (/^PARAMETERS:/) {
        $section='parameters';
		next;
	}
	last if /^(ANALYSIS|QUANTIFICATIONS)/;
	if ($section eq 'parameters' && /^\w+/) {
		chomp;
		my ($paramName,$paramValue)=split(/\t\t/,$_);
		if ($paramName eq 'MZXML'){
			my ($anaID,$mzXMLFile)=split(/:/,$paramValue);
			$params{$paramName}{$anaID}=$mzXMLFile;
			($params{'RAW'}{$anaID})=$dbh->selectrow_array("SELECT WIFF_FILE FROM ANALYSIS WHERE ID_ANALYSIS=$anaID");
		}
		else {
			$params{$paramName}=$paramValue;
		}
	}
}
close INFO;

my $numSteps=9;
print FILESTAT "1/$numSteps Generating data files\n";
close FILESTAT;

###> Label isotope information
my ($getLabelInfo,$lightChannelID)=(0,'light');
my $targetPos=0;
if ($params{'CHANNELS'}) {
	$getLabelInfo=1;
	$targetPos=-1;
	foreach my $channelString (split("::",$params{'CHANNELS'})) {
		my ($channelID,$name,$labels) = split(";",$channelString);
		my $hypTechDeltaMass=0.0000; # Hypothetical total delta mass = sum of all the label deltaMass
		foreach my $label (split('@',$labels)) {
			my ($labelModifAlias,$labelModifName,$modifRes,$searchModifName,$modID)=split("#",$label);
			next unless $modID; # Light label is not associated to a $modID
			my ($deltaMass)=($modID > 0 )? $dbh->selectrow_array("SELECT MONO_MASS FROM MODIFICATION WHERE ID_MODIFICATION=$modID") : 0.00000;
			@{$isoLists{$modID}}=($channelID,$deltaMass,$modifRes);
			$hypTechDeltaMass+=$deltaMass;
			#$anaLabelModifs{$modID}=1;
		}
		@{$channels{$channelID}}=($name,$labels,$hypTechDeltaMass);
		$lightChannelID=$channelID if $hypTechDeltaMass < 0.00001;
	}
}

####>Queries to fetch peptide data<####
my $sthGetInstrument=$dbh->prepare("SELECT INSTRUMENT FROM ANALYSIS WHERE ID_ANALYSIS=?");
#my $sthProtInfo=$dbh->prepare("SELECT PA.ID_PROTEIN,IDENTIFIER,SCORE,DB_RANK,VISIBILITY FROM PROTEIN P,PEPTIDE_PROTEIN_ATTRIB PPA,ANALYSIS_PROTEIN PA WHERE P.ID_PROTEIN=PPA.ID_PROTEIN AND P.ID_PROTEIN=PA.ID_PROTEIN AND PA.ID_ANALYSIS=? AND ID_PEPTIDE=?");
my $sthProtInfo=$dbh->prepare("SELECT PA.ID_PROTEIN,IDENTIFIER,SCORE,DB_RANK,VISIBILITY,PPA.ID_PEPTIDE
							  FROM PROTEIN P,PEPTIDE_PROTEIN_ATTRIB PPA,ANALYSIS_PROTEIN PA
							  WHERE P.ID_PROTEIN=PPA.ID_PROTEIN AND P.ID_PROTEIN=PA.ID_PROTEIN AND PA.ID_ANALYSIS=? ORDER BY VISIBILITY DESC"); # ORDER BY is important, Too many peptide matches for GROUP_CONCAT!!!!!
#my $sthGetNumMods=$dbh->prepare("SELECT COUNT(*) FROM PEPTIDE_MODIFICATION WHERE ID_MODIFICATION=? AND ID_PEPTIDE=?");
my $sthGetNumMods=$dbh->prepare("SELECT P.ID_PEPTIDE,COUNT(*) FROM PEPTIDE P,PEPTIDE_MODIFICATION PM WHERE P.ID_PEPTIDE=PM.ID_PEPTIDE AND ID_ANALYSIS=? AND ID_MODIFICATION=? GROUP BY P.ID_PEPTIDE");
my $skipLabelModStrg=(scalar keys %isoLists)? "AND PM.ID_MODIFICATION NOT IN (".join(',',keys %isoLists).")" : "";
#my $sthPepInfoTP=$dbh->prepare("SELECT ID_PEPTIDE,ELUTION_TIME,PEP_SEQ,MISS_CUT,MR_EXP,MR_CALC,MR_OBS,CHARGE,DATA,SPEC_COUNT FROM PEPTIDE WHERE ID_ANALYSIS=? AND VALID_STATUS>=1");# Avoid to take ghost peptides
my $sthPepInfoTP=$dbh->prepare("SELECT P.ID_PEPTIDE,PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),MISS_CUT,MR_EXP,MR_CALC,MR_OBS,CHARGE,DATA,SPEC_COUNT,ELUTION_TIME
								FROM PEPTIDE P
								LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
								WHERE ID_ANALYSIS=? AND VALID_STATUS>=1 $skipLabelModStrg GROUP BY P.ID_PEPTIDE");# Avoid to take ghost peptides
#my $sthPepInfoGP=$dbh->prepare("SELECT ID_PEPTIDE,PEP_SEQ,MISS_CUT,CHARGE,DATA FROM PEPTIDE WHERE ID_ANALYSIS=? AND VALID_STATUS=0");# Take ghost peptides only
##my $sthPepInfoGP=$dbh->prepare("SELECT P.ID_PEPTIDE,PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),MISS_CUT,CHARGE,DATA
##							   FROM PEPTIDE P
##							   LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
##							   WHERE ID_ANALYSIS=? AND VALID_STATUS=0 $skipLabelModStrg GROUP BY P.ID_PEPTIDE");# Take ghost peptides only
#my $sthCountObs=$dbh->prepare("SELECT COUNT(*) FROM OBSERVATION WHERE ID_ANALYSIS=?");
#my $sthInsObs=$dbh->prepare("INSERT INTO OBSERVATION (ID_ANALYSIS,TARGET_POS) VALUES (?,?)");

####>Generating masschroqML file content<####
my $mcqmlString = qq
|<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<masschroq>
|;

my ($sampID,$rawString,$pepFileString,$groupString,$anaRef)=(0,'','','','');
my $rootPepFileName='pepList_';
###> Get the right scan-number for QSTAR mass spectrometer...
my $deltaTime=5.0;
my $projectID;
#my (%mrObs,%mrCalc,%dbRank);
my %dbRank;
my %encodedPTMs; # to store aa-encoded PTM
my $numAna=scalar keys %{$params{'MZXML'}};
my $anaCount=0;
foreach my $anaID (keys %{$params{'MZXML'}}) {
	
	$anaCount++;
	open(FILESTAT,">>$fileStat");
	print FILESTAT "1/$numSteps Generating data files ($anaCount/$numAna)\n";
	close FILESTAT;

	###>Get list of proteins and matching peptides
	my %pep2Prot;
	$sthProtInfo->execute($anaID);
	while ( my ($protID,$identifier,$score,$dbR,$visibility,$pepID)=$sthProtInfo->fetchrow_array ) { #,$pepStrg
		$proteinScore{$anaID}{$identifier}=$score || 0;
		$dbRank{$identifier}=$dbR || 1;
		$pep2Prot{$pepID}=$identifier unless $pep2Prot{$pepID}; # visible prot 1st then others (in case previous extraction has modified MG)
	}
	
	###>Labeling modif
	my %numLabelMods;
	if ($getLabelInfo) {
		foreach my $modID (keys %isoLists) {
			$sthGetNumMods->execute($anaID,$modID);
			while (my ($pepID,$numMod)=$sthGetNumMods->fetchrow_array) {
				$numLabelMods{$modID}{$pepID}=$numMod;
			}
		}
	}
	
	###> Add Observation if not inserted already
	#$sthCountObs->execute($anaID);
	#my ($nbObs)=$sthCountObs->fetchrow_array;
	#if ($nbObs==0) {
	#	$sthInsObs->execute($anaID,$targetPos);
	#}
	$sthGetInstrument->execute($anaID);
	my ($instrument)=$sthGetInstrument->fetchrow_array;
	$instrument='' unless $instrument;

	my $dFormat =($params{'MZXML'}{$anaID} =~ /mzXML/) ?  'mzxml' : 'mzml';
	$rawString .="<data_file format=\"$dFormat\" id=\"samp$sampID\" path=\"$quantifDir/$params{'MZXML'}{$anaID}\" type=\"$params{'RAWDATA_ACQUISITION'}\"/>\n";
	$groupString .= "samp$sampID ";
	if ($params{REFERENCE} && $anaID == $params{REFERENCE}) {# In isotope case... no reference!
		$anaRef="samp$sampID";
	}

	####>Writing pepList file<####
	my $pepListFileName="$quantifDir/$rootPepFileName$anaID.txt";
	
	$projectID=&promsMod::getProjectID($dbh,$anaID,'analysis') unless $projectID;
	symlink("$promsPath{tmp}/upload/project_$projectID/$params{'MZXML'}{$anaID}","$promsPath{tmp}/quantification/$quantifDate/$params{'MZXML'}{$anaID}");
# symlink("$promsPath{tmp}/upload/project_25/$params{'MZXML'}{$anaID}","$promsPath{tmp}/quantification/$quantifDate/$params{'MZXML'}{$anaID}"); ###### TMP <-------------
	my ($handler,$xmlparser,$scan_allInfo);
	
	### Only for QSTAR files (WIFF)
	if ($instrument eq 'ESI-QUAD-TOF') {
		$handler=mzXMLHandler->new("$quantifDir/$params{'MZXML'}{$anaID}");
		$xmlparser = XML::SAX::ParserFactory->parser(Handler => $handler );
		$xmlparser->parse_uri("$quantifDir/$params{'MZXML'}{$anaID}");
		$scan_allInfo=$handler->get_scanAllInfo();
	}

	open(PEPFILE, ">$pepListFileName") || die "Cannot open $pepListFileName";
	print PEPFILE "scan\tsequence\tmh\tz\tproteins\tmods\n";
	my @sortedScanNumbers=sort{$scan_allInfo->{$a}{'RT'}<=>$scan_allInfo->{$b}{'RT'}} keys %{$scan_allInfo} if ($instrument eq 'ESI-QUAD-TOF');
	
	###> 1st loop: Get only validated peptides (no ghost ones with VALID_STATUS=0)
	$sthPepInfoTP->execute($anaID);
	my $pepDataString='';
	my $pepCount=0;
	while (my ($idPeptide,$pepSeq,$vmod,$missCut,$mrExp,$mrCalc,$mrObs,$charge,$data,$sc,$elutionTime) = $sthPepInfoTP->fetchrow_array) {
		#my $vmod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$idPeptide,$anaID,$pepSeq,\%anaLabelModifs);
		$seqMissedCuts{$pepSeq}=$missCut;
		$vmod = '' unless $vmod;
		$vmods{$idPeptide}=$vmod;
		$spectralCount{$idPeptide}=$sc;
		#$mrObs{"$pepSeq;$vmod;$charge"}=$mrObs;
		#$mrCalc{"$pepSeq;$vmod;$charge"}=$mrCalc;
		$mrCalc+=1.007825;
		my $scanNumber=-1;
		my ($timeMin,$timeMax);
		###> For Orbitrap, it is already there.
		if ($elutionTime=~/sc(\d+)/) { $scanNumber=$1; }
		if ($elutionTime =~ /to/) {
			($timeMin,$timeMax)=split(/ to /,$elutionTime);
		}
		elsif ($scanNumber == -1 && $elutionTime=~/et(\d+\.\d+)/) {
			($timeMin,$timeMax)=($1,$1);
		}
		else {
			($timeMin,$timeMax)=($elutionTime,$elutionTime);
		}

		### Only for QSTAR files (WIFF)
		#if ($instrument eq 'ESI-QUAD-TOF'){ # the scanNumber does not exist in DAT file... it has to be found by reading the mzXML file
		#	foreach my $scanEvent (@sortedScanNumbers) {
		#		next if $scan_allInfo->{$scanEvent}{'RT'} < $timeMin*60-$deltaTime;
		#		$scanNumber=$scanEvent;
		#		last if $scan_allInfo->{$scanEvent}{'RT'} > $timeMax*60+$deltaTime;
		#	}
		#}
		
		###> For labelled peptides...
		#if ($vmod =~ /Lys\+8/){# labelled peptide -> $mrExp has to be recomputed!
		#        my $nbK;
		#        foreach my $res (split(//,$pepSeq)){
		#                if ($res eq 'K'){
		#                    $nbK+=1;
		#                }
		#        }
		#        $mrExp-=8.0*$nbK;# Non labelled peptide mass
		#}
		$pep2Ana{$idPeptide}=$anaID; # Only for no-ghost peptides
		my $pepChannelID=-1;
		if ($getLabelInfo) {
			###> Check if the peptide contains a label modification and tag it !
			foreach my $modID (keys %isoLists) {
				##$sthGetNumMods->execute($modID,$idPeptide);
				##my ($numMod)=$sthGetNumMods->fetchrow_array;
				my ($channelID,$deltaMass,$modifRes)=@{$isoLists{$modID}};
				if (($numLabelMods{$modID} && $numLabelMods{$modID}{$idPeptide}) || ($deltaMass==0 && $pepChannelID==-1)) {
					$pepChannelID=$channelID;
				}
			}
		}
		$pepInfo{"$pepSeq;$vmod;$charge;$anaID;$pepChannelID"}{"SEQVMODCHARGE"}=$idPeptide;
		$pepInfo{"$pepSeq;$vmod;$charge;$anaID;$pepChannelID"}{"DATA"}=$data;
		
		##print PEPFILE "$scanNumber\t$pepSeq\t$mrCalc\t$charge\t$pep2Prot{$idPeptide}\t$idPeptide:$mrCalc:$pepChannelID\n";
		#$pepDataString.="$scanNumber\t$pepSeq\t$mrCalc\t$charge\t$pep2Prot{$idPeptide}\t$idPeptide:$mrCalc:$pepChannelID\n";
my $aaPtmCode=&aaEncodePTM($vmod,\%encodedPTMs,\@aaList);
		$pepDataString.="$scanNumber\t$pepSeq$aaPtmCode\t$mrCalc\t$charge\t$pep2Prot{$idPeptide}\t$idPeptide:$charge:$mrCalc:$pepChannelID\n"; # added charge in source pepID
		if (++$pepCount >= 1000) {
			print PEPFILE $pepDataString;
			$pepDataString='';
			$pepCount=0;
		}
	}
	print PEPFILE $pepDataString if $pepCount;
	close PEPFILE;
	$pepFileString .="<peptide_file data=\"samp$sampID\" path=\"$pepListFileName\"/>\n";

	###> 2nd loop: Get only ghost peptides (with VALID_STATUS=0)
	##$sthPepInfoGP->execute($anaID);
	##while (my ($idPeptide,$pepSeq,$vmod,$missCut,$charge,$data) = $sthPepInfoGP->fetchrow_array) {
	##	$seqMissedCuts{$pepSeq}=$missCut;
	##	#my $vmod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$idPeptide,$anaID,$pepSeq,\%anaLabelModifs);
	##	$vmod='' unless $vmod;
	##	$pep2Ana{$idPeptide}=$anaID;
	##	$spectralCount{$idPeptide}=0;
	##	my $pepChannelID=-1;
	##	if ($getLabelInfo) {
	##		###> Check if the peptide contains a label modification and tag it !
	##		foreach my $modID (keys %isoLists) {
	##			##$sthGetNumMods->execute($modID,$idPeptide);
	##			##my ($numMod)=$sthGetNumMods->fetchrow_array;
	##			my ($channelID,$deltaMass,$modifRes)=@{$isoLists{$modID}};
	##			if (($numLabelMods{$modID} && $numLabelMods{$modID}{$idPeptide}) || ($deltaMass==0 && $pepChannelID==-1)) {
	##				$pepChannelID=$channelID;
	##			}
	##		}
	##	}
	##	my $pepKey="$pepSeq;$vmod;$charge;$anaID;$pepChannelID";
	##	next if $pepInfo{$pepKey}; # {"SEQVMODCHARGE"}) # A no-ghost peptide was already referenced
	##	$pepInfo{$pepKey}{"SEQVMODCHARGE"}=$idPeptide; # Only for ghost peptides
	##	$pepInfo{$pepKey}{"DATA"}=$data;
	##	
	##	##$sthProtInfo->execute($anaID,$idPeptide); <===== Moved up for analysis-wide query (PP 22/05/20)
	##	#####> Add score of proteins or virtual proteins - 19/08/13
	##	##while ( my ($protID,$identifier,$score,$dbR,$visibility)=$sthProtInfo->fetchrow_array ){
	##	##	next unless defined($protID);
	##	##	$score=0 unless $score;
	##	##	$proteinScore{$anaID}{$protID}=$score unless defined($proteinScore{$anaID}{$protID});
	##	##}
	##}
	$sampIDs{"samp$sampID"}=$anaID;
	$sampID++;
}
$sthGetInstrument->finish;
$sthProtInfo->finish;
$sthPepInfoTP->finish;
##$sthPepInfoGP->finish;
$sthGetNumMods->finish;
#$sthInsObs->finish;
#$sthCountObs->finish;
#$dbh->commit;
$dbh->disconnect;

###>Check for errors
&checkForErrors;

###> File info
$mcqmlString.= qq
|<rawdata>
$rawString</rawdata>
<groups>
<group data_ids="$groupString" id="G1"/>
</groups>
<peptide_files_list>$pepFileString</peptide_files_list>
|;

####> Isotope info (get SILAC heavy and light info)
my @labelOrder;
if ($getLabelInfo) {
	$mcqmlString.="<isotope_label_list>\n";

	my (%channelToRes,%res,@channelOrder);
	foreach my $channelID (keys %channels) {
		my ($name,$labels,$hypTechDeltaMass)=@{$channels{$channelID}};
		foreach my $label (split('@',$labels)) {
			my ($labelModifAlias,$labelModifName,$modifRes,$searchModifName,$modID)=split("#",$label);
			my @isoInfos=@{$isoLists{$modID}};
			my $deltaMass=$isoInfos[1];
			next if $deltaMass==0;
			$channelToRes{$channelID}{$modifRes}=$deltaMass;
			$res{$modifRes}=1;
		}
		push @channelOrder, $channelID;
	}

	# Mandatory to start to write the ID of the isotope by some letters otherwise MassChroQ will crash because the ID is not valid from xs:ID definition
	for (my $i = 0 ; $i < $#channelOrder ; $i++) {
		for (my $j = $i+1 ; $j <= $#channelOrder ; $j++) {
			my ($chanID1,$chanID2)=($channelOrder[$i],$channelOrder[$j]);
			my ($isoPlusString,$isoNegString)=("<isotope_label id=\"ID${chanID1}_to_$chanID2\">\n","<isotope_label id=\"ID${chanID2}_to_$chanID1\">\n");
			foreach my $residu (keys %res) {
				next if !$channelToRes{$chanID1}{$residu} && !$channelToRes{$chanID2}{$residu};
				my $deltaM1=($channelToRes{$chanID1}{$residu})?$channelToRes{$chanID1}{$residu}:0.0000;
				my $deltaM2=($channelToRes{$chanID2}{$residu})?$channelToRes{$chanID2}{$residu}:0.0000;
				my $diff=$deltaM2-$deltaM1;
				$isoPlusString.="<mod at=\"$residu\" value=\"$diff\"/>\n";
				$diff=$deltaM1-$deltaM2;
				$isoNegString.="<mod at=\"$residu\" value=\"$diff\"/>\n";
			}
			$isoPlusString.="</isotope_label>\n";
			$isoNegString.="</isotope_label>\n";
			$mcqmlString.=$isoPlusString.$isoNegString;

			push @labelOrder, ("ID${chanID1}_to_$chanID2","ID${chanID2}_to_$chanID1");
		}
	}

	$mcqmlString.="</isotope_label_list>\n";
}

###> Alignment info
# Alignment is done only for non labelled extraction for the moment
$mcqmlString.="<alignments>\n<alignment_methods>\n";
my $alignmentID; # my_obiwarp or my_ms2
if ($params{'EXTRACTION_ALGO'} eq 'MS2') {
	$alignmentID="my_ms2";
	$mcqmlString.= qq 
|<alignment_method id="$alignmentID">
<ms2>
<ms2_tendency_halfwindow>$params{MS2_TENDENCY}</ms2_tendency_halfwindow>
|;
	$mcqmlString.="<ms2_smoothing_halfwindow>$params{MS2_SMOUTHING}</ms2_smoothing_halfwindow>\n" if $params{MS2_SMOUTHING};
	$mcqmlString.="<ms1_smoothing_halfwindow>$params{MS1_SMOUTHING}</ms1_smoothing_halfwindow>\n" if $params{MS1_SMOUTHING};
	$mcqmlString.="</ms2>\n";
}
else {
	$alignmentID="my_obiwarp";
	$mcqmlString.=qq
|<alignment_method id="$alignmentID">
<obiwarp>
<lmat_precision>1</lmat_precision>
<mz_start>$params{MZ_ALIGN_RANGE_MIN}</mz_start>
<mz_stop>$params{MZ_ALIGN_RANGE_MAX}</mz_stop>
</obiwarp>
|;
}
$mcqmlString.=qq
|</alignment_method>
</alignment_methods>
<align group_id="G1" method_id="$alignmentID" reference_data_id="$anaRef"/>
</alignments>
|;


###> Quanti info
$mcqmlString.=qq
|<quantification_methods>
<quantification_method id="quanti1">
<xic_extraction xic_type="$params{XIC_EXTRACTION_TYPE}">
<$params{XIC_RANGE}_range max="$params{MZTOL_MAX}" min="$params{MZTOL_MIN}"/>
</xic_extraction>
<xic_filters>
|;
$mcqmlString.="<anti_spike half=\"$params{ANTISPIKE}\"/>\n" if $params{ANTISPIKE};
$mcqmlString.="<background half_mediane=\"$params{MED_MIN}\" half_min_max=\"$params{MED_MAX}\"/>\n" if $params{MED_MIN};
$mcqmlString.="<smoothing half=\"$params{SMOOTH}\"/>\n" if $params{SMOOTH};
$mcqmlString.=qq
|</xic_filters>
<peak_detection>
<detection_zivy>
<mean_filter_half_edge>1</mean_filter_half_edge>
<minmax_half_edge>3</minmax_half_edge>
<maxmin_half_edge>2</maxmin_half_edge>
<detection_threshold_on_max>$params{DT_STOP}</detection_threshold_on_max>
<detection_threshold_on_min>$params{DT_START}</detection_threshold_on_min>
</detection_zivy>
</peak_detection>
</quantification_method>
</quantification_methods>
|;

###> Results format
$mcqmlString.=qq
|<quantification>
<quantification_results>
<quantification_result format="tsv" output_file="$quantifDir/result_quanti"/>
<quantification_result format="masschroqml" output_file="$quantifDir/result_quanti_mcqml"/>
</quantification_results>
|;
####> Test of quantification traces...
if ($params{TRACES}) {
	$mcqmlString.=qq
|<quantification_traces>
<all_xics_traces output_dir="$quantifDir/all_xics_traces" format="tsv"/>
</quantification_traces>
|;
}

my $niMinAbundStrg = (defined($params{NI_MIN_ABUND}) && $params{NI_MIN_ABUND}) ? " ni_min_abundance=\"$params{NI_MIN_ABUND}\"" : '';
$mcqmlString.=qq
|<quantify id="q1" quantification_method_id="quanti1" withingroup="G1">
<peptides_in_peptide_list mode="$params{XIC_VAL}"$niMinAbundStrg/>
</quantify>
|;

####> SILAC part-start
if ($getLabelInfo) {
	my $labelOrderStg=join(' ', @labelOrder);
	$mcqmlString.=qq
|<quantify id="q2" quantification_method_id="quanti1" withingroup="G1">
<peptides_in_peptide_list mode="$params{XIC_VAL}"$niMinAbundStrg isotope_label_refs="$labelOrderStg"/>
</quantify>
|;
}

$mcqmlString.=qq
|</quantification>
</masschroq>
|;

####>Writing masschroqML file<####
open (MASSCHROQML,">$quantifDir/quanti.masschroqML"); # item valueR valueDB
print MASSCHROQML $mcqmlString;
close MASSCHROQML;

###>Check for errors
&checkForErrors;

#####################################
####>1st MassChroQ (XML parsing)<####
#####################################
my $mcqTmpDir='/tmp/mcq_'.$quantifDate;
my ($massChroqError,$pbsError)=('','');
my $mcqOutFile1="$quantifDir/masschroq_status1.txt";
my $maxHours1=24;
if ($clusterInfo{'on'}) {
	
	open(FILESTAT,">>$fileStat");
	print FILESTAT "2/$numSteps Waiting for file parsing job to start...\n";
	close FILESTAT;

	my $scriptfile = "$quantifDir/script1.sh";
	open (BASH,"+>",$scriptfile);
	print BASH qq
|#!/bin/bash
export LC_ALL="C"
cd $quantifDir
mkdir -p $mcqTmpDir
$masschroqPath/masschroq --tmpdir $mcqTmpDir --parse-peptides $quantifDir/quanti.masschroqML > $mcqOutFile1 2> $quantifDir/masschroq_errors.txt
|;
	close BASH;
	my $modBash=0775;
	chmod $modBash, $scriptfile;
	my $timeStamp=strftime("%Y%m%d%H%M%S",localtime);
	# my $maxMem=int(3 + 0.1*$numAna).'Gb'; # 1.5 + 0.2*$numAna
	my $coeff=($numAna > 100)? 0.1 : ($numAna > 50)? 0.15 : ($numAna > 25)? 0.2 : ($numAna > 10)? 0.3 : 0.5;
	my $maxMem=int(3.5 + $coeff*$numAna).'Gb';
	my %jobParams=(
		maxMem=>$maxMem,
		numCPUs=>1,
		maxHours=>$maxHours1,
		jobName=>"myProMS_MCQ1_$quantifID",
		commandFile=>'PBScommand1.sh',
		jobIdTag=>'1',
		outFile=>'PBS1.txt',
		errorFile=>'PBSerror1.txt',
		jobEndFlag=>"_END_$timeStamp",
		noWatch=>1
	);
	($pbsError,my $pbsErrorFile)=$clusterInfo{'runJob'}->($quantifDir,$scriptfile,\%jobParams);
	sleep 15;
	my $hasStarted=0;
	my $nbWhile=0;
	my $numParsed=0;
	while ((!-e "$quantifDir/PBS1.txt" || !`tail -3 $quantifDir/PBS1.txt | grep _END_`) && !$massChroqError && !$pbsError) {
		if ($nbWhile > $maxHours1*4*60) { # sleep 15 => 4 loops/min * 60min/h * $maxHours1 h
			die "MassChroQ is taking too long or died before completion";
		}
		sleep 15;
		if (!$hasStarted && -e $mcqOutFile1) { # Job has started
			open(FILESTAT,">>$fileStat");
			print FILESTAT "2/$numSteps MassChroQ is parsing XML files\n";
			close FILESTAT;
			$hasStarted=1;
			next;
		}
		#>Check for progress
		$numParsed=&checkFileParseProgress($mcqOutFile1,$numParsed);
		
		#>Check for error
		$massChroqError=`head -5 $quantifDir/masschroq_errors.txt` if (-e "$quantifDir/masschroq_errors.txt" && -s _); # _ filehandle from previous -e
		$pbsError=$clusterInfo{'checkError'}->($pbsErrorFile);

		$nbWhile++;
	}
	#>Check for error AGAIN after waiting
	unless ($massChroqError) { # Check again
		sleep 10;
		$massChroqError=`head -5 $quantifDir/masschroq_errors.txt` if (-e "$quantifDir/masschroq_errors.txt" && -s _); # _ filehandle from previous -e
	}
}
else {
	open(FILESTAT,">>$fileStat");
	print FILESTAT "2/$numSteps MassChroQ is parsing XML files\n";
	close FILESTAT;

	my $childPid = fork;
	unless ($childPid) { # child here
		system "cd $quantifDir; mkdir -p $mcqTmpDir; $masschroqPath/masschroq --tmpdir $mcqTmpDir --parse-peptides $quantifDir/quanti.masschroqML > $mcqOutFile1 2> $quantifDir/masschroq_errors.txt"; # if no cd, this command does not work because no permission to use cgi-bin as temporary directory!
		exit;
	}
	my $nbWhile=0;
	my $numParsed=0;
	while (1) {
		if ($nbWhile > 4*60*$maxHours1) { # sleep 15 => 4 loops/min * 60min/h * $maxHours1 h
			die "MassChroQ is taking too long or died before completion";
		}
		sleep 15;
		
		#>Check for error
		$massChroqError=`head -5 $quantifDir/masschroq_errors.txt` if (-e "$quantifDir/masschroq_errors.txt" && -s _); # _ filehandle from previous -e
		last if $massChroqError;

		my $res = waitpid($childPid,WNOHANG); # WNOHANG (from POSIX ":sys_wait_h"): parent doesn't hang during waitpid
		$numParsed=&checkFileParseProgress($mcqOutFile1,$numParsed);
		last if $res; # child has ended
		
		$nbWhile++;
	}
}

####>Child job ERROR Management<####
my $errorStrg1=(!-d $quantifDir)? "The job has been manually interrupted"
			: ($pbsError)? "The computation cluster has generated the following Error:\n$pbsError" # && $pbsError !~ /WARNING/
			: ($massChroqError)? "MassChroQ has generated the following Error:\n$massChroqError"  # && $massChroqError !~ /WARNING/
			: '';
if ($errorStrg1) {
	$dbh=&promsConfig::dbConnect('no_user');
	$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2,UPDATE_DATE=NOW() WHERE ID_QUANTIFICATION=$quantifID");
	$dbh->commit;
	$dbh->disconnect;
	die $errorStrg1;
}

###> Checking & updating version
my ($version,$goodVersion)=('?',0);
open(MCQSTAT,$mcqOutFile1);
while (my $line=<MCQSTAT>) {
	if ($line =~ /Running MassChroQ (\S+)/) {
		$version=$1;
		$goodVersion=1 if $version=~/^2\.2/; #Running MassChroQ >=2.2.2 with XML file
		last;
	}
}
close MCQSTAT;
$dbh=&promsConfig::dbConnect('no_user');
$dbh->do("UPDATE QUANTIFICATION SET QUANTIF_ANNOT=REPLACE(QUANTIF_ANNOT,'::SOFTWARE=MCQ;?','::SOFTWARE=MCQ;$version') WHERE ID_QUANTIFICATION=$quantifID");
$dbh->commit;
$dbh->disconnect;

unless($goodVersion) {
	$dbh=&promsConfig::dbConnect('no_user');
	$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2,UPDATE_DATE=NOW() WHERE ID_QUANTIFICATION=$quantifID");
	$dbh->commit;
	$dbh->disconnect;
	die "It seems that you have a version of MassChroQ software not supported (2.2 or 2.2.12 optimized version)\n";
}

###> 2nd step: read and add allChargeStates peptides.
if ($params{'ALLCHARGESTATES'}) {# check-box was checked
	###> Read the parse-peptide XML file
	open(PEPTIDEPARSED,"$quantifDir/parsed-peptides_quanti.masschroqML"); # peptides parsed by masschroq in the 1st step
	open(ALLCHARGE,">$quantifDir/parsed-peptides2_quanti.masschroqML");
	my ($onPeptide,$mh)=0;
	my $chargeString='';
	my $lineCount=0;
	while (my $line=<PEPTIDEPARSED>) {
		$chargeString.=$line;
		$lineCount++;
		if ($line =~ /<peptide mods=\".* mh=\"(.*)\" id=/){
			$mh=$1;
			$onPeptide=1;
		}
		elsif ($line =~ /<\/peptide>/) {
			$onPeptide=0;
		}

		if ($onPeptide && $line =~ /<observed_in z=\"(\d+)\" data=\"(.*)\" scan=\"(\d+)\"\/>/) {
			my ($z,$sample,$scanNum)=($1,$2,$3);
			my $massPep=$z*($mh-1.007825);
			for (my $i=$z+1 ; $i < $maxStates ; $i++){
				my $chSt=($massPep+$i*1.007825)/$i;
				if ($chSt > $params{MZ_RANGE_MIN} && $chSt < $params{MZ_RANGE_MAX}) {
					$chargeString.="   <observed_in z=\"$i\" data=\"$sample\" scan=\"$scanNum\"/>\n";
					$lineCount++;
				}
			}
		}
		if ($lineCount >= 1000) {
			print ALLCHARGE $chargeString;
			$lineCount=0;
		}
	}
	print ALLCHARGE $chargeString if $lineCount;
	close ALLCHARGE;
	close PEPTIDEPARSED;
	system "mv $quantifDir/parsed-peptides2_quanti.masschroqML $quantifDir/parsed-peptides_quanti.masschroqML"; # overwrite original
}
###>Check for errors
&checkForErrors;

########################################
####>2nd MassChroQ (XIC extraction)<####
########################################
my $mcqOutFile2="$quantifDir/masschroq_status2.txt";
$massChroqError=$pbsError='';
my $maxHours2=10+12*$numAna;
####> Use cluster
if ($clusterInfo{'on'}) {

	open(FILESTAT,">>$fileStat");
	print FILESTAT "3/$numSteps Waiting for extraction job to start...\n";
	close FILESTAT;

	####> Script command <####
	my $nbProcs=($sampID > 19)? 8 : ($sampID > 9)? 6 : 4;
	my $scriptfile = "$quantifDir/script2.sh";
	open (BASH,"+>",$scriptfile);
	print BASH qq
|#!/bin/bash
export LC_ALL="C"
cd $quantifDir
mkdir -p $mcqTmpDir
$masschroqPath/masschroq --tmpdir $mcqTmpDir --cpus $nbProcs $quantifDir/parsed-peptides_quanti.masschroqML > $mcqOutFile2 2> $quantifDir/masschroq_errors.txt
|;
	close BASH;
	my $modBash=0775;
	chmod $modBash, $scriptfile;
	# my $maxHours=10+12*$numAna;
	#my $maxMem=(2*$numAna > 100) ? "100Gb" : 2*$numAna.'Gb';
	#my $maxMem=(2*$numAna > 100)? (100+($numAna-50)).'Gb' : 2*$numAna.'Gb';
	#my $maxMem=(1+$numAna).'Gb';
	#my $maxMem=int(3 + 0.25*$numAna).'Gb';
	my $coeff=($numAna > 100)? 0.3 : ($numAna > 50)? 0.4 : ($numAna > 25)? 0.6 : ($numAna > 10)? 0.8  : 1;
	my $maxMem=int(3.5 + $coeff*$numAna).'Gb';
	my $timeStamp=strftime("%Y%m%d%H%M%S",localtime);
	my %jobParams=(
		maxMem=>$maxMem,
		numCPUs=>$nbProcs,
		maxHours=>$maxHours2,
		jobName=>"myProMS_MCQ2_$quantifID",
		commandFile=>'PBScommand2.sh',
		jobIdTag=>'2',
		outFile=>'PBS2.txt',
		errorFile=>'PBSerror2.txt',
		jobEndFlag=>"_END_$timeStamp",
		noWatch=>1
	);
	($pbsError,my $pbsErrorFile)=$clusterInfo{'runJob'}->($quantifDir,$scriptfile,\%jobParams);
	sleep 15;

	###> When MassChroQ is finished, it is written DONE in the status_file -> need to check if mcq is finished every 30sec.
	my $hasStarted=0;
	my $nbWhile=0;
	my $progMessage='';
	while ((!-e $mcqOutFile2 || !`tail -3 $mcqOutFile2 | grep DONE`) && (!-e "$quantifDir/PBS2.txt" || !`tail -3 $quantifDir/PBS2.txt | grep _END_`) && !$massChroqError && !$pbsError) {
		if ($nbWhile > $maxHours2*60) { # sleep 60 => 1 loop/min * 60min/h * $maxHours2 h
			die "MassChroQ is taking too long or died before completion";
		}
		sleep 60;
		if (!$hasStarted && -e $mcqOutFile2) { # Job has started
			open(FILESTAT,">>$fileStat");
			print FILESTAT "3/$numSteps MassChroQ is performing XIC extraction\n";
			close FILESTAT;
			$hasStarted=1;
			next;
		}
		#>Check for progress
		$progMessage=&checkExtractionProgress($mcqOutFile2,$progMessage);
		
		#>Check for error
		$massChroqError=`head -5 $quantifDir/masschroq_errors.txt` if (-e "$quantifDir/masschroq_errors.txt" && -s _); # _ filehandle from previous -e
		$pbsError=$clusterInfo{'checkError'}->($pbsErrorFile);
		
		$nbWhile++;
	}
	#>Check for error AGAIN after waiting
	unless ($massChroqError) { # Check again
		sleep 10;
		$massChroqError=`head -5 $quantifDir/masschroq_errors.txt` if (-e "$quantifDir/masschroq_errors.txt" && -s _); # _ filehandle from previous -e
	}
}
else {
	open(FILESTAT,">>$fileStat");
	print FILESTAT "3/$numSteps MassChroQ is performing XIC extraction\n";
	close FILESTAT;

	my $childPid = fork;
	unless ($childPid) { # child here
		system "cd $quantifDir; mkdir -p $mcqTmpDir; $masschroqPath/masschroq --tmpdir $mcqTmpDir $quantifDir/parsed-peptides_quanti.masschroqML > $mcqOutFile2 2> $quantifDir/masschroq_errors.txt";
		exit;
	}
	my $nbWhile=0;
	my $progMessage='';
	while (1) {
		if ($nbWhile > $maxHours2*60) { # sleep 60 => 1 loop/min * 60min/h * $maxHours2 h
			die "MassChroQ is taking too long or died before completion";
		}
		sleep 60;
		
		#>Check for error
		$massChroqError=`head -5 $quantifDir/masschroq_errors.txt` if (-e "$quantifDir/masschroq_errors.txt" && -s _); # _ filehandle from previous -e
		last if $massChroqError;

		my $res = waitpid($childPid,WNOHANG); # WNOHANG (from POSIX ":sys_wait_h"): parent doesn't hang during waitpid
		$progMessage=&checkExtractionProgress($mcqOutFile2,$progMessage);
		last if $res; # child has ended
		
		$nbWhile++;
	}
}

####>ERROR Management<####
my $errorStrg2=(! -d $quantifDir)? "The job has been manually interrupted"
			: ($pbsError)? "The computation cluster has generated the following Error:\n$pbsError" #  && $pbsError !~ /WARNING/
			: ($massChroqError)? "MassChroQ has generated the following Error:\n$massChroqError" #  && $massChroqError !~ /WARNING/
			: '';
if ($errorStrg2) {
	$dbh=&promsConfig::dbConnect('no_user');
	$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2,UPDATE_DATE=NOW() WHERE ID_QUANTIFICATION=$quantifID");
	$dbh->commit;
	$dbh->disconnect;
	die $errorStrg2;
}

######################################################
####>Parsing MassChroQ quantification information<####
######################################################
open(FILESTAT,">>$fileStat");
print FILESTAT "4/$numSteps Checking for pre-existing MBWR\n";
close FILESTAT;

$dbh=&promsConfig::dbConnect('no_user');
my ($areaCode)=($getLabelInfo)? 'ISO_AREA' : 'XIC_AREA';
my %qparamsID;
my ($sthGetParams)=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,CODE FROM QUANTIFICATION_PARAMETER QP, QUANTIFICATION Q WHERE Q.ID_QUANTIFICATION_METHOD=QP.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=$quantifID");
$sthGetParams->execute;
while (my ($qpID,$code) = $sthGetParams->fetchrow_array) {
		$qparamsID{$code}=$qpID;
}
$sthGetParams->finish;


####>Load pre-existing ghost peptides data to prevent peptide duplication<####
# Moved down to optimize memery usage
$sthGetNumMods=$dbh->prepare("SELECT P.ID_PEPTIDE,COUNT(*) FROM PEPTIDE P,PEPTIDE_MODIFICATION PM WHERE P.ID_PEPTIDE=PM.ID_PEPTIDE AND ID_ANALYSIS=? AND ID_MODIFICATION=? GROUP BY P.ID_PEPTIDE"); # No "my" because already declared above
my $sthPepInfoGP=$dbh->prepare("SELECT P.ID_PEPTIDE,PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),MISS_CUT,CHARGE,DATA
							   FROM PEPTIDE P
							   LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
							   WHERE ID_ANALYSIS=? AND VALID_STATUS=0 $skipLabelModStrg GROUP BY P.ID_PEPTIDE");# Take ghost peptides only

foreach my $anaID (keys %{$params{'MZXML'}}) {
	
	###>Labeling modif
	my %numLabelMods;
	if ($getLabelInfo) {
		foreach my $modID (keys %isoLists) {
			$sthGetNumMods->execute($anaID,$modID);
			while (my ($pepID,$numMod)=$sthGetNumMods->fetchrow_array) {
				$numLabelMods{$modID}{$pepID}=$numMod;
			}
		}
	}
	
	$sthPepInfoGP->execute($anaID);
	while (my ($idPeptide,$pepSeq,$vmod,$missCut,$charge,$data) = $sthPepInfoGP->fetchrow_array) {
		$seqMissedCuts{$pepSeq}=$missCut;
		#my $vmod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$idPeptide,$anaID,$pepSeq,\%anaLabelModifs);
		$vmod='' unless $vmod;
		$pep2Ana{$idPeptide}=$anaID;
		$spectralCount{$idPeptide}=0;
		my $pepChannelID=-1;
		if ($getLabelInfo) {
			###> Check if the peptide contains a label modification and tag it !
			foreach my $modID (keys %isoLists) {
				##$sthGetNumMods->execute($modID,$idPeptide);
				##my ($numMod)=$sthGetNumMods->fetchrow_array;
				my ($channelID,$deltaMass,$modifRes)=@{$isoLists{$modID}};
				if (($numLabelMods{$modID} && $numLabelMods{$modID}{$idPeptide}) || ($deltaMass==0 && $pepChannelID==-1)) {
					$pepChannelID=$channelID;
				}
			}
		}
		my $pepKey="$pepSeq;$vmod;$charge;$anaID;$pepChannelID";
		next if $pepInfo{$pepKey}; # {"SEQVMODCHARGE"}) # A no-ghost peptide was already referenced
		$pepInfo{$pepKey}{"SEQVMODCHARGE"}=$idPeptide; # Only for ghost peptides
		$pepInfo{$pepKey}{"DATA"}=$data;
		$pepInfo{$pepKey}{"IS_GHOST"}=1;
	}
}
$sthGetNumMods->finish;
$sthPepInfoGP->finish;
$dbh->disconnect;

###>Check for errors
&checkForErrors;

my (%quantifPeptides,%retTime,%mergedPeptides); # Peptides inserted into PEPTIDE_QUANTIFICATION (identified in an analysis and valid)
my (%matchedGhostPeptides,%newPeptides,%referencePeptides); # Potential gost-peptides.

### Read XIC Traces files if this option is checked
my %traceFiles;
if ( $params{TRACES} ) {
	opendir (DIR, "$quantifDir/all_xics_traces");
	while (defined (my $traceFile = readdir (DIR))) {
		next if ($traceFile =~ m/^\./); # avoid . and .. filenames
		# ex of trace file: q1_G1_samp2_pep23_mz_741-904_rt_4680-88_l_741-604_h_742-204_z_2.tsv
		my @nameCut=split(/_|\.|-/,$traceFile);
		# $nameCut[0] <-> q1
		# $nameCut[1] <-> G1
		# $nameCut[2] <-> samp2
		# $nameCut[3] <-> pep23
		# $nameCut[5] <-> 741
		# $nameCut[-2] <-> 2
		$traceFiles{"$nameCut[0]:$nameCut[1]:$nameCut[2]:$nameCut[3]:$nameCut[5]:$nameCut[-2]"}=$traceFile;
	}
	close DIR;

	###>Check for errors
	&checkForErrors;
}

################################################
####>Remove merged PTM isoform from results<####
################################################
open(FILESTAT,">>$fileStat");
print FILESTAT "5/$numSteps Checking for co-eluted ions...\n";
close FILESTAT;
my %mergedXICs;
my ($totalIons,$totalMerged,$totalExcluded)=&cleanCoelutingIons("5/$numSteps Checking for co-eluted ions",'peptides_q1_G1.tsv',$sampID,\%pep2Ana,\%sampIDs,\%vmods,\%mergedXICs);

##################################################################################################################
####>Filter peptides which have a large dispersion of retention time across samples (considered not reliable)<####
##################################################################################################################
my $currStep=5;
my $mcqResults;
if ($getLabelInfo) {  # No filter on SILAC data
	$mcqResults = "peptides_q1_G1.tsv";
}
else {  # Filter peptides which have a large dispersion of retention time across samples (considered not reliable)
	$currStep++;
	$mcqResults = "peptides_q1_G1_filtered.tsv";
	
	my $filterOutFile="$quantifDir/result_quanti.d/filtering_log.out";
	my $filterErrorFile="$quantifDir/filtering_errors.txt";
	my $RoutFile = "$quantifDir/filterMassChroQ.Rout";
	my $filterError=$pbsError='';
	my $maxHoursR=10+12*$numAna;

	####> Use cluster
	if ($clusterInfo{'on'}) {

		open(FILESTAT,">>$fileStat");
		print FILESTAT "$currStep/$numSteps Waiting for filtering job to start...\n";
		close FILESTAT;
		
		####> Script command <####
		my $scriptfile = "$quantifDir/script3.sh";
		open (BASH,"+>",$scriptfile);
		print BASH qq
|#!/bin/bash
export LC_ALL="C"
cd $quantifDir
$pathR/R CMD BATCH --no-save --no-restore '--args $quantifDir $params{RT_SD_FILTER} $params{RT_SD_MAX} $params{RT_SD_MIN}' $promsPath{R_scripts}/filterMassChroQ.R 2> $filterErrorFile
|;
		close BASH;
		my $modBash=0775;
		chmod $modBash, $scriptfile;

		# my $maxMem=(int(1.5 + 0.75*$numAna)).'Gb';
		# my $maxMem=int(3 + 0.2*$numAna).'Gb';
		my $coeff=($numAna > 100)? 0.3 : ($numAna > 50)? 0.4 : ($numAna > 25)? 0.6 : ($numAna > 10)? 0.8  : 1;
		my $maxMem=int(3.5 + $coeff*$numAna).'Gb';
		my $timeStamp=strftime("%Y%m%d%H%M%S",localtime);
		my %jobParams=(
			maxMem=>$maxMem,
			numCPUs=>1,
			maxHours=>$maxHoursR,
			jobName=>"myProMS_filter.MCQ_$quantifID",
			commandFile=>'PBScommand3.sh',
			jobIdTag=>'3',
			outFile=>'PBS3.txt',
			errorFile=>'PBSerror3.txt',
			jobEndFlag=>"_END_$timeStamp",
			noWatch=>1
		);
		($pbsError,my $pbsErrorFile)=$clusterInfo{'runJob'}->($quantifDir,$scriptfile,\%jobParams);
		sleep 15;

		###> When filtering is finished, it is written "Ended..."" in the output file
		my $hasStarted=0;
		my $nbWhile=0;
		my $progMessage='';
		while ((!-e $filterOutFile || !`tail -3 $filterOutFile | grep ^Ended`) && (!-e "$quantifDir/PBS3.txt" || !`tail -3 $quantifDir/PBS3.txt | grep _END_`) && !$filterError && !$pbsError) {
			if ($nbWhile > $maxHoursR*60) { # sleep 60 => 1 loop/min * 60min/h * $maxHours2 h
				my $errorMessage=(-e "$quantifDir/PBS3.txt")? "Filtering job is taking too long or died silently before completion" : "Filtering job is taking too long to start";
				die $errorMessage;
			}
			sleep 60;
			if (!$hasStarted && -e $RoutFile) { # Job has started
				open(FILESTAT,">>$fileStat");
				print FILESTAT "$currStep/$numSteps Filtering retention time outliers\n";
				close FILESTAT;
				$hasStarted=1;
				next;
			}
			#>Check for progress
			$progMessage=&checkFilteringProgress($filterOutFile,$progMessage);
			
			#>Check for error
			$filterError=`head -5 $filterErrorFile` if (-e $filterErrorFile && -s _); # _ filehandle from previous -e
			$pbsError=$clusterInfo{'checkError'}->($pbsErrorFile);
			
			$nbWhile++;
		}
	}
	else { #### "local" launch
		open(FILESTAT,">>$fileStat");
		print FILESTAT "$currStep/$numSteps Filtering retention time outliers\n";
		close FILESTAT;
		my $childPid = fork;
		unless ($childPid) { # child here
			system "cd $quantifDir; $pathR/R CMD BATCH --no-save --no-restore '--args $quantifDir $params{RT_SD_FILTER} $params{RT_SD_MAX} $params{RT_SD_MIN}' $promsPath{R_scripts}/filterMassChroQ.R 2> $filterErrorFile";
			exit;
		}
		my $nbWhile=0;
		my $progMessage='';
		while (1) {
			if ($nbWhile > $maxHoursR*60) { # sleep 60 => 1 loop/min * 60min/h * $maxHoursR h
				die "Filtering is taking too long or died silently before completion";
			}
			sleep 60;

			#>Check for error
			$filterError=`head -5 $filterErrorFile` if (-e $filterErrorFile && -s _); # _ filehandle from previous -e
			last if $filterError;
			
			my $res = waitpid($childPid,WNOHANG); # WNOHANG (from POSIX ":sys_wait_h"): parent doesn't hang during waitpid
			$progMessage=&checkFilteringProgress($filterOutFile,$progMessage);
			last if $res; # child has ended
			
			$nbWhile++;
		}
	}
	
	# Check for errors after MassChroQ filter
	my $filterErrorStrg='';
	if (!-d $quantifDir) {$filterErrorStrg="The job has been manually interrupted";}
	elsif ($pbsError) {$filterErrorStrg="The computation cluster has generated the following Error:\n$pbsError";}
	else {
		if (-e $RoutFile) {
			my $RoutStrg = `tail -3 $RoutFile`;
			unless ($RoutStrg =~ /proc\.time\(\)/) {
				$RoutStrg = `tail -20 $RoutFile`;
				$filterErrorStrg = "Peptides filtering has generated the following error:\n";
				my $inError = 0;
				foreach my $line (split(/\n/, $RoutStrg)) {
					next if (!$inError && $line !~ /^Error in/ && $line !~ /^Error:/); # skip everything before "Error in..."
					$inError = 1;
					$filterErrorStrg .= "\n$line";
				}
			}
		}
		else {$filterErrorStrg = "No filtering output found! Check if job was started";}
	}
	if ($filterErrorStrg) {
		$dbh=&promsConfig::dbConnect('no_user');
		$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2,UPDATE_DATE=NOW() WHERE ID_QUANTIFICATION=$quantifID");
		$dbh->commit;
		$dbh->disconnect;
		die $filterErrorStrg;
	}
}

$currStep++;
open(FILESTAT,">>$fileStat");
print FILESTAT "$currStep/$numSteps Processing peptide XIC data\n";
close FILESTAT;
#open(PEPQUANTI,"$quantifDir/result_quanti_pep.tsv"); # former file of MassChroQ that was keeping q1 and q2 quantifications
#open(PEPQUANTI1,"$quantifDir/result_quanti.d/peptides_q1_G1.tsv"); # label-free area extraction data
open(PEPQUANTI1,"$quantifDir/result_quanti.d/$mcqResults"); # label-free area extraction data, filtered on peptides retention time dispersion
while (my $line=<PEPQUANTI1>) {
	# next if $line =~ /group/; #skip first line
	next if $.== 1; #skip first line
	chomp($line);
	# example: q1,G1,samp0,F5849FD,518.25061,1249.1633,409292.7,7127488.8,1201.6282,1264.3375,pep855,isotopeName,GWDVGVAGMK,2,613281
	my ($quantifMCQID,$group,$msrun,$msrunfile,$mz,$rt,$maxintensity,$area,$rtbegin,$rtend,$peptide,$isotope,$sequence,$z,$mods)=split("\t",$line);
	# Be careful for VMOD choice: some isobaric peptides can be grouped together in $mods string with different positions or from different samples <---- No longer true due to PTM encoded in sequence (PP 24/05/21)
	# Therefore, to avoid the creation of Ghost-Peptides with incorrect VMOD for this specific entry, 2 priority rules to follow :
	# 1st: choose the vmod from the same msrun peptide (and the highest number of spectral-counts if many possibilities)
	# 2nd: choose the vmod from the closest RT-time provided compared to the rt time available of other experiment
	my ($pepID,$mz2,$pepChannelID);
	my $sampXicKey="$msrun:$rt:$maxintensity";
	my $sameSample=0;
	foreach my $pepEntry (split(/\s/,$mods)) { # each pepEntry looks like this: $idPeptide:$z:$mz:$targetPos
		my ($pepIDloc,$zloc,$mz2loc,$pepChannelIDloc)=split(/:/,$pepEntry);
		next if $z != $zloc; # make sure we are comparing the same charge state
		if ($pep2Ana{$pepIDloc} == $sampIDs{$msrun}) { # make sure pepID comes from the current sample/Ana/XLS file
			if ($pepID) { # an ion was already matched to this XIC
				if ($sameSample) { # antoher ion also matches => keep ion with highest spectral count
					if ($spectralCount{$pepIDloc} > $spectralCount{$pepID}) {
						($pepID,$mz2,$pepChannelID)=($pepIDloc,$mz2loc,$pepChannelIDloc);
					}
				}
				else { # previously selected ion was not from current sample (kept to create ghost ion later) 
					($pepID,$mz2,$pepChannelID)=($pepIDloc,$mz2loc,$pepChannelIDloc); # select this ion
				}
			}
			else { # 1st ion to match & from good sample
				($pepID,$mz2,$pepChannelID)=($pepIDloc,$mz2loc,$pepChannelIDloc); # select this ion
			}
			$sameSample=1;
		}
		else { # ion is not from current sample => process anyway to create ghost ion later
			if ($pepID) { # antoher ion also matches => keep ion with highest spectral count
				if (!$sameSample) {
					if ($spectralCount{$pepIDloc} > $spectralCount{$pepID}) {
						($pepID,$mz2,$pepChannelID)=($pepIDloc,$mz2loc,$pepChannelIDloc);
					}
				}
			}
			else { # 1st ion to match & from another sample
				($pepID,$mz2,$pepChannelID)=($pepIDloc,$mz2loc,$pepChannelIDloc);
			}
		}
	}
	my $vmod=$vmods{$pepID};
	$vmod='' unless $vmod;
	my $mzInt= int $mz;
	if (!$isotope) {
		my $pepKey="$sequence;$vmod;$z;$sampIDs{$msrun};$pepChannelID";
		if ($pepInfo{$pepKey} && $pepInfo{$pepKey}{"SEQVMODCHARGE"}) { # peptide exists in current anaID: truely identified or ghost from a previous quanti
			if ($quantifPeptides{$pepKey}) { # already seen in previous line (should no longer happen with PTM encoded in sequence)
				if ($area > $quantifPeptides{$pepKey}) {
					$quantifPeptides{$pepKey}=$area;
					$retTime{$pepKey}{"BEGIN"}=$rtbegin;
					$retTime{$pepKey}{"END"}=$rtend;
					$retTime{$pepKey}{"APEX"}=$rt;
					$retTime{$pepKey}{"TRACE"}=$traceFiles{"q1:$group:$msrun:$peptide:$mzInt:$z"} if $params{TRACES};
					# no need to update $retTime{$pepKey}{"MERGED"} <- same
				}
			}
			else {
				$matchedGhostPeptides{ $sampIDs{$msrun} }{ $pepInfo{$pepKey}{"SEQVMODCHARGE"} }=1 if $pepInfo{$pepKey}{"IS_GHOST"}; # {anaID}{pepID}
				$quantifPeptides{$pepKey}=$area;
				$retTime{$pepKey}{"BEGIN"}=$rtbegin;
				$retTime{$pepKey}{"END"}=$rtend;
				$retTime{$pepKey}{"APEX"}=$rt;
				$retTime{$pepKey}{"TRACE"}=$traceFiles{"q1:$group:$msrun:$peptide:$mzInt:$z"} if $params{TRACES};
$retTime{$pepKey}{"MERGED"}=$mergedXICs{$sampXicKey} if $mergedXICs{$sampXicKey}; # list of representative pepID (1 per isoform accross whole dataset)
			}
		}
		else { # new ghost peptide... has to be created
			if ($newPeptides{$pepKey} && $newPeptides{$pepKey}{"AREA"}) { # already seen in a previous line
				if ($area > $newPeptides{$pepKey}{"AREA"}) {
					$newPeptides{$pepKey}{"AREA"}=$area;
					$newPeptides{$pepKey}{"BEGIN"}=$rtbegin;
					$newPeptides{$pepKey}{"END"}=$rtend;
					$newPeptides{$pepKey}{"APEX"}=$rt;
					$newPeptides{$pepKey}{"TRACE"}=$traceFiles{"q1:$group:$msrun:$peptide:$mzInt:$z"} if $params{TRACES};
				}
			}
			else {
				$referencePeptides{$pepID}=1;
				$newPeptides{$pepKey}{"AREA"}=$area;
				$newPeptides{$pepKey}{"RT"}=$rt/60;
				$newPeptides{$pepKey}{"MZ"}=$mz;
				$newPeptides{$pepKey}{"MREXP"}=$mz*$z-1.007825;
				$newPeptides{$pepKey}{"REFPEPID"}=$pepID;
				$newPeptides{$pepKey}{"BEGIN"}=$rtbegin;
				$newPeptides{$pepKey}{"END"}=$rtend;
				$newPeptides{$pepKey}{"APEX"}=$rt;
				$newPeptides{$pepKey}{"TRACE"}=$traceFiles{"q1:$group:$msrun:$peptide:$mzInt:$z"} if $params{TRACES};
				my $data=$pepInfo{"$sequence;$vmod;$z;$pep2Ana{$pepID};$pepChannelID"}{'DATA'};
				my ($PRS) = ($data =~ /(##PRS=[^#]+)/) if $data;
				$PRS='' unless $PRS;
				$newPeptides{$pepKey}{"DATA"}=$PRS; # in case there is PhosphoRS information
$newPeptides{$pepKey}{"MERGED"}=$mergedXICs{$sampXicKey} if $mergedXICs{$sampXicKey}; # list of representative pepID (1 per isoform accross whole dataset)
			}
		}
	}

	#if ( $getLabelInfo ) { # Analysis made to get labeling info
	#		my $traceFile=($params{TRACES})? $traceFiles{"q2:$group:$msrun:$peptide:$mzInt:$z"} : "";
	#		@{$isoInfo{$sampIDs{$msrun}}{"$sequence;$vmod;$z"}{$mz}}=($area,$rt/60,$mods,$isotope,$rtbegin,$rtend,$traceFile);
	#}
}
close PEPQUANTI1;

###>Check for errors
&checkForErrors;

if ($getLabelInfo) {
	$currStep++;
	open(FILESTAT,">>$fileStat");
	print FILESTAT "$currStep/$numSteps Processing isotope data\n";
	close FILESTAT;

	open(PEPQUANTI2,"$quantifDir/result_quanti.d/peptides_q2_G1.tsv"); # isotope area extraction data
	while (my $line=<PEPQUANTI2>) {
		next if $line =~ /group/; #skip first line
		chomp($line);
		# example: q2,G1,samp0,F5849FD,518.25061,1249.1633,409292.7,7127488.8,1201.6282,1264.3375,pep855,isotopeName,GWDVGVAGMK,2,613281
		my ($quantifMCQID,$group,$msrun,$msrunfile,$mz,$rt,$maxintensity,$area,$rtbegin,$rtend,$peptide,$isotope,$sequence,$z,$mods)=split("\t",$line);
# $sequence=~s/PPPPP[A-Z]+//; # remove aaPtmCode if any
		my ($pepID,$mz2,$pepChannelID);
		my $sameSample=0;

		foreach my $pepEntry (split(/\s/,$mods)) { # each pepEntry looks like this: $idPeptide:$mz:$targetPos
			my ($pepIDloc,$zloc,$mz2loc,$pepChannelIDloc)=split(/:/,$pepEntry);
			next if $z != $zloc; # make sure we are comparing the same charge state
			if ($pep2Ana{$pepIDloc} == $sampIDs{$msrun}) {
				if ($pepID) {
					if ($sameSample) {
						if ($spectralCount{$pepIDloc} > $spectralCount{$pepID}) {
							($pepID,$mz2,$pepChannelID)=($pepIDloc,$mz2loc,$pepChannelIDloc);
						}
					}
					else {
						($pepID,$mz2,$pepChannelID)=($pepIDloc,$mz2loc,$pepChannelIDloc);
					}
				}
				else {
					($pepID,$mz2,$pepChannelID)=($pepIDloc,$mz2loc,$pepChannelIDloc);
				}
				$sameSample=1;
			}
			elsif ($pepID) {
				if (!$sameSample) {
					if ($spectralCount{$pepIDloc} > $spectralCount{$pepID}) {
						($pepID,$mz2,$pepChannelID)=($pepIDloc,$mz2loc,$pepChannelIDloc);
					}
				}
			}
			else {
				($pepID,$mz2,$pepChannelID)=($pepIDloc,$mz2loc,$pepChannelIDloc);
			}
		}
		my $vmod=$vmods{$pepID};
		$vmod='' unless $vmod;
		my $mzInt= int $mz;

		my $traceFile=($params{TRACES})? $traceFiles{"q2:$group:$msrun:$peptide:$mzInt:$z"} : "";
		# @{$isoInfo{$sampIDs{$msrun}}{"$sequence;$vmod;$z"}{$mz}}=($area,$rt/60,$mods,$isotope,$rtbegin,$rtend,$traceFile);
		@{$isoInfo{$sampIDs{$msrun}}{"$sequence;$vmod;$z"}{$mz}}=($area,$rt/60,"$pepID:$z:$mz2:$pepChannelID",$isotope,$rtbegin,$rtend,$traceFile);
	}
	close PEPQUANTI2;

	######
	### Part for isotope extraction with MassChroQ
	######
	foreach my $anaID (keys %isoInfo) {
		my $qsetNum=0;
		foreach my $seqVmodCharge (keys %{$isoInfo{$anaID}}) {
			$qsetNum++;
			my ($seq,$vmod,$charge)=split(/;/,$seqVmodCharge);
			my %qsetValues=();
			###> 1st round: get isotope retrieved peptides.
			foreach my $mz (sort {$a <=> $b} keys %{$isoInfo{$anaID}{$seqVmodCharge}}) {
				my ($area,$rt,$pepIDmhchannelID,$isotope,$rtbegin,$rtend,$traceFile)=@{$isoInfo{$anaID}{$seqVmodCharge}{$mz}};
				#print "AREA=$area RT=$rt MODS=$pepIDmhchannelID ISOTOPE=$isotope\n";
				#my ($pepID,$mh,$channelIDOriginal,@pepIDs)=split(/ |:/,$pepIDmhchannelID);
				my ($pepID,$zloc,$mh,$channelIDOriginal)=split(/:/,$pepIDmhchannelID);
				if ($isotope) {
					$isotope=~s/ID//;
					my ($fromChanID,$toChanID)=split('_to_',$isotope);
					if ($fromChanID eq $channelIDOriginal && $channels{$toChanID}) {
						$referencePeptides{$pepID}=1;
						$qsetValues{$qsetNum}{$toChanID}{'AREA'}=$area;
						$qsetValues{$qsetNum}{$toChanID}{'REFPEPID'}=$pepID;
						$qsetValues{$qsetNum}{$toChanID}{'RT'}=$rt;
						$qsetValues{$qsetNum}{$toChanID}{'MZ'}=$mz;
						$qsetValues{$qsetNum}{$toChanID}{'Z'}=$charge;
						$qsetValues{$qsetNum}{$toChanID}{'BEGIN'}=$rtbegin;
						$qsetValues{$qsetNum}{$toChanID}{'END'}=$rtend;
						$qsetValues{$qsetNum}{$toChanID}{'TRACE'}=$traceFile if $params{TRACES};
						my $pepKey="$seq;$vmod;$charge;$anaID;$toChanID";
						if ($pepInfo{$pepKey}{"SEQVMODCHARGE"} && !$qsetValues{$qsetNum}{$toChanID}{'PEP_ID'}) { # This test prevent inflation of ghost peptides !
							$qsetValues{$qsetNum}{$toChanID}{'PEP_ID'}=$pepInfo{$pepKey}{"SEQVMODCHARGE"};
							if (!$quantifPeptides{$pepKey} ) { # Apply only for ghost peptides that were not quantified
								$quantifPeptides{$pepKey}=$area;
								$retTime{$pepKey}{"BEGIN"}=$rtbegin;
								$retTime{$pepKey}{"END"}=$rtend;
								$retTime{$pepKey}{"APEX"}=$rt*60;
								$retTime{$pepKey}{"TRACE"}=$traceFile if $params{TRACES};
							}
						}
					}
				}
			}
			###> 2nd round: get alignment retrieved peptides.
			foreach my $mz (sort {$a <=> $b} keys %{$isoInfo{$anaID}{$seqVmodCharge}}) {
				my ($area,$rt,$pepIDmhchannelID,$isotope,$rtbegin,$rtend,$traceFile)=@{$isoInfo{$anaID}{$seqVmodCharge}{$mz}};
				#print "AREA=$area RT=$rt MODS=$pepIDmhchannelID ISOTOPE=$isotope\n";
				#my ($pepID,$mh,$channelIDOriginal,@pepIDs)=split(/ |:/,$pepIDmhchannelID);
				my ($pepID,$zloc,$mh,$channelIDOriginal)=split(/:/,$pepIDmhchannelID);
				if (!$isotope) { # If isotope column empty -> original version of the isotope $channelIDOriginal
					#next if $pep2Ana{$pepID} && $pep2Ana{$pepID} != $anaID;
					$qsetValues{$qsetNum}{$channelIDOriginal}{'AREA'}=$area;
					if ($pepInfo{"$seq;$vmod;$charge;$anaID;$channelIDOriginal"}{"SEQVMODCHARGE"}) {
						$qsetValues{$qsetNum}{$channelIDOriginal}{'PEP_ID'}=$pepID if $pep2Ana{$pepID} == $anaID; # 1st priority: get the real area from the peptide of this analysis
						$qsetValues{$qsetNum}{$channelIDOriginal}{'PEP_ID'}=$pepID if !$qsetValues{$qsetNum}{$channelIDOriginal}{'PEP_ID'}; # 2nd: get the area retrieved by alignment ONLY if the peptide does not exist for this analysis!
					}
					else {
						$referencePeptides{$pepID}=1;
						$qsetValues{$qsetNum}{$channelIDOriginal}{'REFPEPID'}=$pepID;
						$qsetValues{$qsetNum}{$channelIDOriginal}{'RT'}=$rt;
						$qsetValues{$qsetNum}{$channelIDOriginal}{'MZ'}=$mz;
						$qsetValues{$qsetNum}{$channelIDOriginal}{'Z'}=$charge;
						$qsetValues{$qsetNum}{$channelIDOriginal}{'BEGIN'}=$rtbegin;
						$qsetValues{$qsetNum}{$channelIDOriginal}{'END'}=$rtbegin;
						$qsetValues{$qsetNum}{$channelIDOriginal}{'TRACE'}=$traceFile if $params{TRACES};
					}
				}
			}
			###> QSET VALUES
			my $data="##MCQSET_$quantifID=$qsetNum";
			#print "Into myProMS MCQSET_$quantiID=$qsetNum\n";
			foreach my $qChannelID (keys %{$qsetValues{$qsetNum}}) {
				next if !$qsetValues{$qsetNum}{$qChannelID}{'AREA'}; # Area for this Channel/Isotope was not found...
				my $pepKey="$seq;$vmod;$charge;$anaID;$qChannelID";
				if ( !$qsetValues{$qsetNum}{$qChannelID}{'PEP_ID'} ) {
					my $refPepID=$qsetValues{$qsetNum}{$qChannelID}{'REFPEPID'};
					$referencePeptides{$refPepID}=1;
					$newPeptides{$pepKey}{"AREA"}=$qsetValues{$qsetNum}{$qChannelID}{'AREA'};
					$newPeptides{$pepKey}{"RT"}=$qsetValues{$qsetNum}{$qChannelID}{'RT'};
					$newPeptides{$pepKey}{"MZ"}=$qsetValues{$qsetNum}{$qChannelID}{'MZ'};
					$newPeptides{$pepKey}{"MREXP"}=$qsetValues{$qsetNum}{$qChannelID}{'MZ'}*$qsetValues{$qsetNum}{$qChannelID}{'Z'}-1.007825;
					$newPeptides{$pepKey}{"REFPEPID"}=$refPepID;
					my $dataref=$pepInfo{"$seq;$vmod;$charge;$pep2Ana{$refPepID};$qChannelID"}{'DATA'};
					my ($PRS) = ($dataref =~ /(##PRS=[^#]+)/) if $dataref;
					$PRS='' unless $PRS;
					$newPeptides{$pepKey}{"DATA"}="$PRS$data";
					$newPeptides{$pepKey}{"BEGIN"}=$qsetValues{$qsetNum}{$qChannelID}{'BEGIN'};
					$newPeptides{$pepKey}{"END"}=$qsetValues{$qsetNum}{$qChannelID}{'END'};
					$newPeptides{$pepKey}{"APEX"}=$qsetValues{$qsetNum}{$qChannelID}{'RT'};
					$newPeptides{$pepKey}{"TRACE"}=$qsetValues{$qsetNum}{$qChannelID}{'TRACE'} if $params{TRACES};
				}
				else {
					$pepInfo{$pepKey}{"DATA"}.=$data;
				}
				#$pepInfo{"$seq;$vmod;$charge;$anaID"}{"CHANNEL"}=$qChannelID;
				#my $pepID=($qsetValues{$qsetNum}{$qChannelID}{'PEP_ID'})?$qsetValues{$qsetNum}{$qChannelID}{'PEP_ID'}: 'Ghost-Peptide';
				#print "CHANNEL=$qChannelID VALUE=$area PEP_ID=$pepID\n";
			}
		}
	}

	###>Check for errors
	&checkForErrors;
}

######################
####>Moving files<####
######################
$currStep++;
open(FILESTAT,">>$fileStat");
print FILESTAT "$currStep/$numSteps Storing peptide XIC data\n";
close FILESTAT;
my $finalDir="$promsPath{quantification}/project_$projectID/quanti_$quantifID";
# Copy everything but XML
system "mkdir -p $finalDir";
system "cp $quantifDir/*.txt $finalDir/."; # peplist & PBS files
system "cp $quantifDir/*.sh $finalDir/." if glob "$quantifDir/*.sh"; # cluster bash files
system "mv $quantifDir/result_quanti.d $finalDir/.";# tsv output files: result_quanti_pep.tsv, result_quanti_prot.tsv and result_quanti_compar.tsv
system "cp $quantifDir/*.time $finalDir/." if glob "$quantifDir/*.time"; # Alignment files -> important to perform again the search
system "cp $quantifDir/*.masschroqML $finalDir/.";
system "mv $quantifDir/all_xics_traces $finalDir/." if $params{TRACES};
system "cp $quantifDir/*.Rout $finalDir" if glob "$quantifDir/*.Rout";  # R output (RT filter)
system "cp $quantifDir/*.out $finalDir/.";

my %pepQuantHandle;
if ($getLabelInfo) { # labeled quanti
	foreach my $channelID (sort{$a<=>$b} keys %channels) {
		open($pepQuantHandle{$channelID},">$finalDir/peptide_quantification_$channelID.txt");
		print {$pepQuantHandle{$channelID}} "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n";
	}
}
else { # label-free
	open($pepQuantHandle{0},">$finalDir/peptide_quantification.txt");
	print {$pepQuantHandle{0}} "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\tMERGED\n";
}
if ($params{TRACES}) {
	open (TRACEFILE,">$finalDir/pepID_to_traces.txt");
	print TRACEFILE "PEPID\tTRACE\n";
}

$dbh=&promsConfig::dbConnect('no_user');
my $sthUpPep=$dbh->prepare("UPDATE PEPTIDE SET DATA=? WHERE ID_PEPTIDE=?");

my $numInsUp=0;
foreach my $seqVmodZanaIDChanID (keys %quantifPeptides ) {
	my ($seq,$vmod,$charge,$anaID,$channelID)=split(/;/,$seqVmodZanaIDChanID);
	my $targetPos=($getLabelInfo)? $channelID : 0;
	my $dataStrg="$qparamsID{$areaCode}\t$pepInfo{$seqVmodZanaIDChanID}{SEQVMODCHARGE}\t$quantifPeptides{$seqVmodZanaIDChanID}";
	$dataStrg.="\t$retTime{$seqVmodZanaIDChanID}{MERGED}" if $retTime{$seqVmodZanaIDChanID}{'MERGED'};
	$dataStrg.="\n";
	$dataStrg.="$qparamsID{RT_BEGIN}\t$pepInfo{$seqVmodZanaIDChanID}{SEQVMODCHARGE}\t$retTime{$seqVmodZanaIDChanID}{BEGIN}\n" if $qparamsID{'RT_BEGIN'};
	$dataStrg.="$qparamsID{RT_END}\t$pepInfo{$seqVmodZanaIDChanID}{SEQVMODCHARGE}\t$retTime{$seqVmodZanaIDChanID}{END}\n" if $qparamsID{'RT_END'};
	$dataStrg.="$qparamsID{RT_APEX}\t$pepInfo{$seqVmodZanaIDChanID}{SEQVMODCHARGE}\t$retTime{$seqVmodZanaIDChanID}{APEX}\n";
	print {$pepQuantHandle{$targetPos}} $dataStrg;

	$sthUpPep->execute($pepInfo{$seqVmodZanaIDChanID}{"DATA"},$pepInfo{$seqVmodZanaIDChanID}{"SEQVMODCHARGE"}) if $getLabelInfo;
	print TRACEFILE "$pepInfo{$seqVmodZanaIDChanID}{\"SEQVMODCHARGE\"}\t$retTime{$seqVmodZanaIDChanID}{\"TRACE\"}\n" if $params{TRACES};
	$numInsUp++;
	$dbh->commit if !$numInsUp % 5000;
}
$dbh->commit if $numInsUp;

#############################################################
###> Queries for insertion of ghost peptides information <###
#############################################################
my @allRefPeptides=sort{$a<=>$b} keys %referencePeptides; # <=> %newPeptides
my (%referencePepModifData,%referencePepPPA);
while (my @subPepList=splice(@allRefPeptides,0,1000)) { # chuncks of 1000 peptides
	my $pepIdStrg=join(',',@subPepList);
	my $queryGetPepM="SELECT ID_PEPTIDE,ID_MODIFICATION,MODIF_TYPE,POS_STRING,REF_POS_STRING FROM PEPTIDE_MODIFICATION WHERE ID_PEPTIDE IN ($pepIdStrg)";
	$queryGetPepM.=' AND ID_MODIFICATION NOT IN ('.join(',',keys %isoLists).')' if $getLabelInfo;
	my $sthGetPepM=$dbh->prepare($queryGetPepM);
	$sthGetPepM->execute;
	while (my ($pepID,@modifData)=$sthGetPepM->fetchrow_array) {
		push @{$referencePepModifData{$pepID}},\@modifData;
	}
	$sthGetPepM->finish;
	my $sthGetPPA=$dbh->prepare("SELECT PPA.ID_PEPTIDE,PPA.ID_PROTEIN,-ABS(PEP_BEG),-ABS(PEP_END),FLANKING_AA,IS_SPECIFIC,IDENTIFIER FROM PEPTIDE_PROTEIN_ATTRIB PPA,PROTEIN P WHERE P.ID_PROTEIN=PPA.ID_PROTEIN AND ID_PEPTIDE IN ($pepIdStrg)");
	$sthGetPPA->execute;
	while (my ($pepID,@ppaData)=$sthGetPPA->fetchrow_array) {
		push @{$referencePepPPA{$pepID}},\@ppaData;
	}
	$sthGetPPA->finish;
}
#my $sthInsPPA=$dbh->prepare("INSERT INTO PEPTIDE_PROTEIN_ATTRIB (ID_PROTEIN,ID_PEPTIDE,ID_ANALYSIS,PEP_BEG,PEP_END,FLANKING_AA) VALUES (?,?,?,?,?,?)");
my $sthAddPA=$dbh->prepare("INSERT INTO PEPTIDE_PROTEIN_ATTRIB (ID_PEPTIDE,ID_ANALYSIS,ID_PROTEIN,PEP_BEG,PEP_END,FLANKING_AA,IS_SPECIFIC) VALUES (?,?,?,?,?,?,?)");
#my $sthGetPPA=$dbh->prepare("SELECT PPA.ID_PROTEIN,PEP_BEG,PEP_END,FLANKING_AA,IS_SPECIFIC,IDENTIFIER FROM PEPTIDE_PROTEIN_ATTRIB PPA,PROTEIN P WHERE P.ID_PROTEIN=PPA.ID_PROTEIN AND ID_PEPTIDE=?");
my $sthInsGhostPep=$dbh->prepare("INSERT INTO PEPTIDE (ID_ANALYSIS,PEP_SEQ,PEP_LENGTH,MISS_CUT,MR_EXP,MR_OBS,CHARGE,ELUTION_TIME,VALID_STATUS,DATA) VALUES (?,?,?,?,?,?,?,?,0,?)"); # VALID_STATUS=0 and SCORE IS NULL characterize ghost peptides
#my $queryGetPepM="SELECT ID_MODIFICATION,MODIF_TYPE,POS_STRING,REF_POS_STRING FROM PEPTIDE_MODIFICATION WHERE ID_PEPTIDE=?";
#$queryGetPepM.=' AND ID_MODIFICATION NOT IN ('.join(',',keys %isoLists).')' if $getLabelInfo;
#my $sthGetPepM=$dbh->prepare($queryGetPepM);
my $sthInsPepM=$dbh->prepare("INSERT INTO PEPTIDE_MODIFICATION (ID_PEPTIDE,ID_MODIFICATION,MODIF_TYPE,POS_STRING,REF_POS_STRING) VALUES (?,?,?,?,?)");
#my $sthCheckProt=$dbh->prepare("SELECT 1 FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=? AND ID_PROTEIN=?");
my %virtualProteins;
foreach my $seqVmodZanaIDChanID (keys %newPeptides) {
	next if $quantifPeptides{$seqVmodZanaIDChanID}; # QUANTIF_INFORMATION already added
	my ($seq,$vmod,$charge,$anaID,$channelID)=split(/;/,$seqVmodZanaIDChanID);
	$vmod='' unless $vmod;
	###> 1st: New instance in PEPTIDE
	my $et=sprintf "rt%.2f",$newPeptides{$seqVmodZanaIDChanID}{"RT"};
	$sthInsGhostPep->execute($anaID,$seq,length($seq),$seqMissedCuts{$seq},$newPeptides{$seqVmodZanaIDChanID}{"MREXP"},$newPeptides{$seqVmodZanaIDChanID}{"MZ"},$charge,$et,$newPeptides{$seqVmodZanaIDChanID}{"DATA"});
	my $newpepID=$dbh->last_insert_id(undef,undef,'PEPTIDE','ID_PEPTIDE');
	$matchedGhostPeptides{$anaID}{$newpepID}=1;
	my $refPepID=$newPeptides{$seqVmodZanaIDChanID}{'REFPEPID'};
	#$sthUpPep->execute($newPeptides{$seqVmodZanaIDChanID}{"DATA"},$newpepID) if $getLabelInfo; # just for isotope-tagged peptides!!!
	###> 2nd : Add VMOD information
	#$sthGetPepM->execute($refPepID);
	#while (my ($modID,$modifType,$posString,$refPosString) = $sthGetPepM->fetchrow_array) {
	#	$sthInsPepM->execute($newpepID,$modID,$modifType,$posString,$refPosString);
	#}
	foreach my $refModData (@{$referencePepModifData{$refPepID}}) {
		$sthInsPepM->execute($newpepID,@{$refModData});
	}
	###> 2nd-bis : Add VMOD-labeled information
	if ($getLabelInfo) {
		foreach my $modID (keys %isoLists) {
			my ($isochannelID,$isodeltaMass,$isomodifRes)=@{$isoLists{$modID}};
			next if $isochannelID != $channelID;
			my %positions;
			foreach my $residueR (split(//,$isomodifRes)) { # Look for each position of the residue where the isotope can occur in pepSeq
				my $offset = length($seq);
				my $posMacth=rindex($seq, $residueR, $offset);
				while ($posMacth != -1) {
					$positions{$posMacth+1}=1;
					$offset=$posMacth-1;
					$posMacth=rindex($seq, $residueR, $offset)
				}
			}
			my $posString=join('.',sort{$a<=>$b} keys %positions);
			$sthInsPepM->execute($modID,$newpepID,'V',$posString,undef) if $posString;
		}
	}
	###> 3rd: New instance in PEPTIDE_PROTEIN_ATTRIB
	#$sthGetPPA->execute($refPepID);
	#while (my ($protID,$pepBeg,$pepEnd,$flankingAA,$isSpecific,$identifier) = $sthGetPPA->fetchrow_array) { #}
	foreach my $refPPA (@{$referencePepPPA{$refPepID}}) {
		#$sthAddPA->execute($newpepID,$anaID,$protID,-$pepBeg,-$pepEnd,$flankingAA,$isSpecific);
		my $identifier=$refPPA->[5];
		$sthAddPA->execute($newpepID,$anaID,@{$refPPA}[0..4]);
		#$dbh->commit;
		#$sthGetNumProt->execute($anaID,$protID);
		#my ($numID)=$sthGetNumProt->fetchrow_array;
		#$sthInsAP->execute($anaID,$protID) if $numID==0;
		#$dbh->commit;
		if (!defined($proteinScore{$anaID}) || !defined($proteinScore{$anaID}{$identifier})) {
			$proteinScore{$anaID}{$identifier}=0; # never seen in analysis => (new!) virtual prot
			#my $protID=$refPPA->[0];
			#$sthCheckProt->execute($anaID,$protID);
			#my ($ok)=$sthCheckProt->fetchrow_array;
			$virtualProteins{$anaID}{$identifier}=1; # unless $ok;
		}
	}
	###> 4th: Write PEPTIDE QUANTIFICATION to file
	my $targetPos=($getLabelInfo)? $channelID : 0;
	my $dataStrg="$qparamsID{$areaCode}\t$newpepID\t$newPeptides{$seqVmodZanaIDChanID}{AREA}";
	$dataStrg.="\t$newPeptides{$seqVmodZanaIDChanID}{MERGED}" if $newPeptides{$seqVmodZanaIDChanID}{'MERGED'};
	$dataStrg.="\n";
	$dataStrg.="$qparamsID{RT_BEGIN}\t$newpepID\t$newPeptides{$seqVmodZanaIDChanID}{BEGIN}\n" if $qparamsID{'RT_BEGIN'};
	$dataStrg.="$qparamsID{RT_END}\t$newpepID\t$newPeptides{$seqVmodZanaIDChanID}{END}\n" if $qparamsID{'RT_END'};
	$dataStrg.="$qparamsID{RT_APEX}\t$newpepID\t$newPeptides{$seqVmodZanaIDChanID}{APEX}\n";
	print {$pepQuantHandle{$targetPos}} $dataStrg;

	$sthUpPep->execute($newPeptides{$seqVmodZanaIDChanID}{"DATA"},$newpepID) if $getLabelInfo;
	print TRACEFILE "$newpepID\t$newPeptides{$seqVmodZanaIDChanID}{\"TRACE\"}\n" if $params{TRACES};
	#$dbh->commit;
	###> 5th: Update quantifPeptides hash information
	$quantifPeptides{"$seqVmodZanaIDChanID"}=1;

	$numInsUp++;
	$dbh->commit if !$numInsUp % 5000;
}
foreach my $targetPos (keys %pepQuantHandle) {close $pepQuantHandle{$targetPos};}
close TRACEFILE if $params{TRACES};
#$sthInsPPA->finish;
$sthAddPA->finish;
#$sthGetPPA->finish;
#$sthGetNumProt->finish;
#$sthInsAP->finish;
$sthInsGhostPep->finish;
#$sthGetPepM->finish;
$sthInsPepM->finish;
#$sthCheckProt->finish;
$sthUpPep->finish;
$dbh->commit;

###>Check for errors
&checkForErrors($dbh);

###########################################
###> Queries for update of matchGroups <###
###########################################
# IMPORTANT:
#	- Ghost peptides are not counted in ANALYSIS_PROTEIN.(NUM_PEP,NUM_MATCH) !!!
# 	- Ghost peptides are counted as normal peptides for Match Group (re)creation
$currStep++;
open(FILESTAT,">>$fileStat");
print FILESTAT "$currStep/$numSteps Updating match group data\n";
close FILESTAT;
my (%numProtTop,%bestProtVis);
my $sthProjVis=$dbh->prepare("SELECT P.ID_PROTEIN,IDENTIFIER,VISIBILITY,COUNT(ID_ANALYSIS) FROM ANALYSIS_PROTEIN AP,PROTEIN P WHERE AP.ID_PROTEIN=P.ID_PROTEIN AND ID_PROJECT=$projectID GROUP BY P.ID_PROTEIN,VISIBILITY");
$sthProjVis->execute;
while (my ($protID,$identifier,$vis,$countAna)=$sthProjVis->fetchrow_array) {
	$vis=0 unless $vis; # to be safe
	if ($vis==2) {$numProtTop{$identifier}=$countAna;}
	elsif (!defined $numProtTop{$identifier}) {$numProtTop{$identifier}=0;}
	$bestProtVis{$identifier}=$vis if (!$bestProtVis{$identifier} || $bestProtVis{$identifier} < $vis);
}
$sthProjVis->finish;

my $sthGetMatchInfo=$dbh->prepare("SELECT P.ID_PEPTIDE,PEP_SEQ,VALID_STATUS,PR.ID_PROTEIN,IDENTIFIER,PROT_LENGTH FROM PEPTIDE P,PEPTIDE_PROTEIN_ATTRIB PPA,PROTEIN PR WHERE P.ID_PEPTIDE=PPA.ID_PEPTIDE AND PPA.ID_PROTEIN=PR.ID_PROTEIN AND P.ID_ANALYSIS=?");
#my $sthUpAP=$dbh->prepare("INSERT INTO ANALYSIS_PROTEIN SET VISIBILITY=?,MATCH_GROUP=?,NUM_PEP=?,NUM_MATCH=?,PEP_SPECIFICITY=?,PEP_COVERAGE=? WHERE ID_PROTEIN=? AND ID_ANALYSIS=?");
#my $sthGetPepInfo=$dbh->prepare("SELECT ABS(PEP_BEG),ABS(PEP_END),PEP_SEQ,P.ID_PEPTIDE FROM PEPTIDE P,PEPTIDE_PROTEIN_ATTRIB PPA WHERE PPA.ID_PEPTIDE=P.ID_PEPTIDE AND PPA.ID_ANALYSIS=? AND ID_PROTEIN=? ORDER BY PEP_BEG ASC, PEP_END ASC"); # Get PEPTIDE and GHOST-PEPTIDE to compute a well coverage on the protein
#my $sthGetAPInfo=$dbh->prepare("SELECT IDENTIFIER,MATCH_GROUP,VISIBILITY FROM ANALYSIS_PROTEIN AP,PROTEIN P WHERE AP.ID_PROTEIN=P.ID_PROTEIN AND ID_ANALYSIS=?");
my ($protVisibility)=$dbh->selectrow_array("SELECT PROT_VISIBILITY FROM PROJECT WHERE ID_PROJECT=$projectID");
$protVisibility=0 unless $protVisibility;
#my $sthInsAP=$dbh->prepare("INSERT IGNORE INTO ANALYSIS_PROTEIN (ID_ANALYSIS,ID_PROTEIN,DB_RANK,CONF_LEVEL,SCORE,NUM_PEP,NUM_MATCH,PEP_COVERAGE,PEP_SPECIFICITY) VALUES (?,?,?,0,0,?,?,?,?)");
my $sthInsAP=$dbh->prepare("INSERT IGNORE INTO ANALYSIS_PROTEIN (ID_ANALYSIS,ID_PROTEIN,DB_RANK,CONF_LEVEL,SCORE,NUM_PEP,NUM_MATCH,PEP_COVERAGE,PEP_SPECIFICITY) VALUES (?,?,?,0,0,0,0,?,?)"); # only for adding virtual proteins
my $sthUpMG=$dbh->prepare("UPDATE ANALYSIS_PROTEIN SET VISIBILITY=?,MATCH_GROUP=? WHERE ID_ANALYSIS=? AND ID_PROTEIN=?");

foreach my $anaID (keys %{$params{'MZXML'}}) {
	my (%matchList,%matchGroup,%visibility,%identifiers,%proteinLength,%boundaryStatus,%pepSpecificity,%pepProteins,%maxCovLength,%numProtPeptides,%proteinPepDens); #%numProtTop,%bestProtVis
	$sthGetMatchInfo->execute($anaID);

	while (my ($pepID,$seq,$validStatus,$protID,$identifier,$protLen)=$sthGetMatchInfo->fetchrow_array) {
		next if (!$validStatus && (!$matchedGhostPeptides{$anaID} || !$matchedGhostPeptides{$anaID}{$pepID})); # skip previous ghosts not matched by this quantif
		#my $vmod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$anaID,$seq,\%anaLabelModifs);
		#$vmod='' unless $vmod;
		$matchList{$identifier}{$seq}=1; #"$seq$vmod"
		$pepProteins{$pepID}{$identifier}=1;
		$identifiers{$identifier}=$protID;
		$proteinLength{$identifier}=$protLen;
	}

	##my $sthTop=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=? AND VISIBILITY=2");
	##my $sthBestVis=$dbh->prepare("SELECT MAX(VISIBILITY) FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=?");
	foreach my $identifier (keys %matchList) {
		##my $protID=$identifiers{$identifier};
		##$sthTop->execute($protID);
		##($numProtTop{$protID})=$sthTop->fetchrow_array;
		##if ($protVisibility) {
		##	$sthBestVis->execute($protID);
		##	($bestProtVis{$identifier})=$sthBestVis->fetchrow_array;
		##}
		#>Peptide density
		$numProtPeptides{$identifier}=scalar (keys %{$matchList{$identifier}});
		$proteinPepDens{$identifier}=$proteinLength{$identifier}/$numProtPeptides{$identifier};
	}
	##$sthTop->finish;
	##$sthBestVis->finish;

	#####> Fill former match-group information so as to only add virtual proteins! <= COMMENTED: MG recreated from start
	##$sthGetAPInfo->execute($anaID);
	##my $maxTopM=0;
	##while (my ($identifier,$matchNumber,$vis) = $sthGetAPInfo->fetchrow_array) {
	##	$maxTopM=$matchNumber if $matchNumber>$maxTopM;
	##	$matchGroup{$identifier}=$matchNumber;
	##	$visibility{$identifier}=$vis;
	##}

	###> Computing the specificity of each peptide for a given protein
	foreach my $peptideID (keys %pepProteins) {
		my $specificity=sprintf "%.2f",100/(scalar keys %{$pepProteins{$peptideID}});
		$specificity=~s/\.00$//;
		foreach my $identifier (keys %{$pepProteins{$peptideID}}) {
			$pepSpecificity{$identifier}=$specificity if (!$pepSpecificity{$identifier} || $pepSpecificity{$identifier}<$specificity);
		}
	}

	###>Taking care of newly added virtual proteins<###
	if ($virtualProteins{$anaID}) {
		###> Pep boundaries
		my (@virtualProtIds,%id2identifier);
		foreach my $identifier (keys %{$virtualProteins{$anaID}}) {
			push @virtualProtIds,$identifiers{$identifier};
			$id2identifier{ $identifiers{$identifier} }=$identifier;
		}
		my $sthGetPepInfo=$dbh->prepare("SELECT ID_PROTEIN,ABS(PEP_BEG),ABS(PEP_END),PEP_SEQ,P.ID_PEPTIDE FROM PEPTIDE P,PEPTIDE_PROTEIN_ATTRIB PPA WHERE PPA.ID_PEPTIDE=P.ID_PEPTIDE AND PPA.ID_ANALYSIS=? AND ID_PROTEIN IN (".join(',',@virtualProtIds).") ORDER BY ID_PROTEIN ASC,PEP_BEG ASC, PEP_END ASC"); # Get PEPTIDE and GHOST-PEPTIDE to compute full coverage on the proteins
		$sthGetPepInfo->execute($anaID);
		while (my ($protID,$pepBeg,$pepEnd,$pSeq,$pepID)=$sthGetPepInfo->fetchrow_array) {
			my $identifier=$id2identifier{$protID};
			$boundaryStatus{$identifier}{$pepBeg}++;
			$boundaryStatus{$identifier}{$pepEnd}--;
			$maxCovLength{$identifier}=$pepEnd if (!$maxCovLength{$identifier} || $maxCovLength{$identifier} < $pepEnd);
		}
		$sthGetPepInfo->finish;
		
		foreach my $identifier (keys %{$virtualProteins{$anaID}}) {# Only virtual-proteins have to be taken care of.
	
			###> Pep Coverage
			##$sthGetPepInfo->execute($anaID,$identifiers{$identifier});
			###my $numMatch=0;
			###my $previousGP='';
			##while (my ($pepBeg,$pepEnd,$pSeq,$pepID)=$sthGetPepInfo->fetchrow_array) {
			##	#my $vmod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$anaID,$pSeq,\%anaLabelModifs);
			##	#$vmod='' unless $vmod;
			##	$boundaryStatus{$identifier}{$pepBeg}++;
			##	$boundaryStatus{$identifier}{$pepEnd}--;
			##	$maxCovLength{$identifier}=$pepEnd if (!$maxCovLength{$identifier} || $maxCovLength{$identifier} < $pepEnd);
			##	#$numMatch++ if ($previousGP ne "$pSeq:$vmod:$pepBeg");
			##	#$previousGP="$pSeq:$vmod:$pepBeg";
			##}
	
			my $pepCoverage=0;
			if ($proteinLength{$identifier} || $proteinLength{$identifier} > 0) {
				my $coverage=0;
				my $hasPeptide=0;
				my $boundaryNter=0;
				foreach my $boundary (sort{$a<=>$b} keys %{$boundaryStatus{$identifier}}) {
					if ($hasPeptide==0) { # start of peptide region (cannot become <0!)
						$boundaryNter=$boundary;
					}
					$hasPeptide+=$boundaryStatus{$identifier}{$boundary};
					if ($hasPeptide==0) { # end of peptide region (should be true for last boundary too)
						$coverage+=($boundary-$boundaryNter)+1;
					}
				}
	
				#$pepCoverage=sprintf "%.1f",(100*$coverage)/$proteinLength{$identifier};
				#$pepCoverage=~s/\.0//;
				$pepCoverage=($maxCovLength{$identifier} <= $proteinLength{$identifier})? sprintf "%.1f",(100*$coverage)/$proteinLength{$identifier} : sprintf "%.1f",(-100*$coverage)/$maxCovLength{$identifier}; # -: flag for protLength problem
				$pepCoverage*=1; # 25.0 -> 25
			}
			###> Pep Specificity
			#$sthInsAP->execute($anaID,$identifiers{$identifier},$dbRank{$identifier},scalar keys %{$matchList{$identifier}},$numMatch,$pepCoverage,$pepSpecificity{$identifier},$visibility{$identifier});
			#$sthInsAP->execute($anaID,$identifiers{$identifier},$dbRank{$identifier},scalar keys %{$matchList{$identifier}},$numMatch,$pepCoverage,$pepSpecificity{$identifier});
			$sthInsAP->execute($anaID,$identifiers{$identifier},$dbRank{$identifier},$pepCoverage,$pepSpecificity{$identifier});
			$numInsUp++;
			$dbh->commit if !$numInsUp % 5000;
		}
	}

####>DEBUG
# open(DEBUG,">>$quantifDir/debug.txt");
# print DEBUG ">ANA=$anaID-----------------------------------\n";
# foreach my $identifier (keys %matchList) {
# 	print DEBUG "numProtPeptides\t$identifier\n" unless defined $numProtPeptides{$identifier};
# 	print DEBUG "numProtTop\t$identifier\n" unless defined $numProtTop{$identifier};
# 	print DEBUG "proteinScore\t$identifier\n" unless defined $proteinScore{$anaID}{$identifier};
# 	print DEBUG "proteinLength\t$identifier\n" unless defined $proteinLength{$identifier};
# 	print DEBUG "proteinPepDens\t$identifier\n" unless defined $proteinPepDens{$identifier};
# 	print DEBUG "identifiers\t$identifier\n" unless defined $identifiers{$identifier};
# }
# close DEBUG;
####<DEBUG

	###>Recreate all Match Groups<###
    my @sortedIdentifiers=sort{$numProtPeptides{$b}<=>$numProtPeptides{$a} || $numProtTop{$b}<=>$numProtTop{$a} || $proteinScore{$anaID}{$b}<=>$proteinScore{$anaID}{$a} || &deltaLength($proteinLength{$a},$proteinLength{$b},$proteinPepDens{$a},$proteinPepDens{$b}) || $proteinLength{$a}<=>$proteinLength{$b} || $identifiers{$a}<=>$identifiers{$b}} keys %matchList;
	&promsMod::createMatchGroups(\%matchList,\%matchGroup,\@sortedIdentifiers,\%visibility,\%bestProtVis,$protVisibility,0);
	#&checkForErrors($dbh);
	
	foreach my $identifier (keys %matchList) {
		$sthUpMG->execute($visibility{$identifier},$matchGroup{$identifier},$anaID,$identifiers{$identifier});
		$numInsUp++;
		$dbh->commit if !$numInsUp % 5000;
	}
	
	###<Update project-wide data (%numProtTop & %bestProtVis) for next analysis
	foreach my $identifier (keys %visibility) {
		$numProtTop{$identifier}++ if $visibility{$identifier}==2;
		$bestProtVis{$identifier}=$visibility{$identifier} if $bestProtVis{$identifier} < $visibility{$identifier};
	}
	
}
#$sthGetPepInfo->finish;
$sthGetMatchInfo->finish;
#$sthGetAPInfo->finish;
$sthInsAP->finish;
$sthUpMG->finish;
$dbh->commit;

###>Check for errors
&checkForErrors($dbh);

#######################
#####>Moving files<####
#######################
#system "mkdir $promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
#system "mkdir $promsPath{quantification}/project_$projectID/quanti_$quantifID";
#system "mv $quantifDir/*.txt $promsPath{quantification}/project_$projectID/quanti_$quantifID";# peplist files
#system "mv $quantifDir/*.tsv $promsPath{quantification}/project_$projectID/quanti_$quantifID";# output files: result_quanti_pep.tsv, result_quanti_prot.tsv and result_quanti_compar.tsv
#system "mv $quantifDir/*.time $promsPath{quantification}/project_$projectID/quanti_$quantifID";# Alignment files -> important to perform again the search
#system "mv $quantifDir/*.masschroqML $promsPath{quantification}/project_$projectID/quanti_$quantifID";
#system "rm $quantifDir/*.mzXML";

####################################
####>Quantification is finished<####
####################################
$dbh->do("UPDATE QUANTIFICATION SET QUANTIF_ANNOT=CONCAT(QUANTIF_ANNOT,'::MERGED_IONS=$totalExcluded,$totalMerged,$totalIons'),STATUS=1,UPDATE_DATE=NOW() WHERE ID_QUANTIFICATION=$quantifID");
$dbh->commit;
$dbh->disconnect;

open(FILESTAT,">>$fileStat");
print FILESTAT "Ended ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
close FILESTAT;

sleep 10;
#unlink $fileStat;

############################
####>Deleting all files<####
############################
#system "rm -rf $quantifDir"; # Leave deleteion to monitoring system?


#####################################
####<Encoding PTMs into sequence>####
#####################################
# Needed for MCQ to distinguish position isoforms (eg. XXSXSR-Phospho:3 & XXSXSR-Phospho:5)
# Adds 5 P (to detect encoded a PTM) then encodes PTM by converting its occurrence rank in a base 20 with AA alphabet: 0=>A,1=>C,...,19=>Y,20=>CA,...
sub aaEncodePTM {
	my ($varMod,$refEncodedPTMs,$refAaList)=@_;
	return '' unless $varMod; # no PTM
	return $refEncodedPTMs->{$varMod} if $refEncodedPTMs->{$varMod}; # already encoded
	my $ptmRank=scalar keys %{$refEncodedPTMs};
	$ptmRank++;
	my $numAABase=20;
	my $aaCode='';
	while ($ptmRank >= $numAABase) {
		my $frac=$ptmRank % $numAABase;
#print "$ptmRank %> $frac\n";
		$aaCode=$refAaList->[$frac].$aaCode;
		$ptmRank-=$frac;
		$ptmRank/=$numAABase;
	}
	$aaCode=$refAaList->[$ptmRank].$aaCode if $ptmRank;
	$aaCode='PPPPP'.$aaCode;
	$refEncodedPTMs->{$varMod}=$aaCode;
	return $aaCode;
}
############################### OBSOLETE
####<Encoding phospho mods>#### Needed for MCQ to distinguish position isoforms (eg. XXSXSR-Phospho:3 & XXSXSR-Phospho:5)
############################### Encode PTM in AA sequence code 5 P followed by AA repeated x position (eg. Phospho:3 -> PPPPAAA, XXSXSR-Phospho:5 -> PPPPAAAAA)
# sub aaEncodePhosphoMod {
# 	my ($phosphoID,$vmod,$refAaList)=@_;
# 	return '' unless $vmod;
# 	my $codeStrg='';
# 	foreach my $modStrg (split(/&/,$vmod)) {
# 		my @mod=split(/[:\.]/,$modStrg);
# 		next unless $mod[0]==$phosphoID;
# 		$codeStrg='PPPPP';
# 		my $prevPos=0;
# 		for (my $i=1; $i<=$#mod; $i++) {
# 			$codeStrg.=$refAaList->[$i]x($mod[$i]-$prevPos);
# 			$prevPos=$mod[$i];
# 		}
# 		last;
# 	}
# 	return $codeStrg;
# }

sub cleanCoelutingIons { # globals: $quantifDir FILESTAT
	# Checks for coelution ions: same RT and max intensity
	# Reports merged isobaric PTM isoforms
	# Excludes merges ions with DIFFERENT sequences
	# Removes encoded PTM from sequence
	my ($rootMessage,$pepQuantiFile,$numSamples,$refPep2Ana,$refSampIDs,$refVmods,$refMergedXICs)=@_;	
	my $inFile="$quantifDir/result_quanti.d/$pepQuantiFile";
	my $outFile="$quantifDir/result_quanti.d/tmp_clean_ions.tsv";
	my $excludedFile="$quantifDir/result_quanti.d/excluded_ions.tsv";
	open(OUT,">$outFile");
	open(EXCL,">$excludedFile");
	my $headerLine=`head -1 $inFile`;
	print OUT $headerLine;
	print EXCL $headerLine;
	my ($totalIons,$totalMerged,$totalExcluded)=(0,0,0);
	for (my $msrunIndex=0; $msrunIndex < $numSamples; $msrunIndex++) {
		unless (($msrunIndex+1) % 5) {
			open(FILESTAT,">>$fileStat");
			print FILESTAT "$rootMessage (Processing sample ",($msrunIndex+1),"/$numSamples...)\n";
			close FILESTAT;
		}
		my $sampFile="$quantifDir/result_quanti.d/tmp_samp$msrunIndex.tsv";
		system "echo '$headerLine' > $sampFile";
		system "grep -w samp$msrunIndex $inFile >> $sampFile"; # "grep >>" adds "\n" between header and data!!!
		
		##<Scanning sampX file
		my %lineData;
		open(IN,$sampFile);
		while (my $line=<IN>) {
			next if $.== 1; # skip header
			# example: q1,G1,samp0,F5849FD,518.25061,1249.1633,409292.7,7127488.8,1201.6282,1264.3375,pep855,isotopeName,GWDVGVAGMK,2,613281
			my ($quantifMCQID,$group,$msrun,$msrunfile,$mz,$rt,$maxintensity,$area,$rtbegin,$rtend,$peptide,$isotope,$sequence)=split("\t",$line);
			next unless $group; # header is followed by empty line because of grep
			my $xicKey="$rt:$maxintensity";
			if (!$lineData{$xicKey} || !$lineData{$xicKey}{PEP}{$peptide}) { # duplicate lines exist with SAME matching ion!!!
				$sequence=~s/PPPPP[A-Z]+//; # remove aaPtmCode if any
				$line=~s/PPPPP[A-Z]+//; # remove aaPtmCode from sequence
				push @{$lineData{$xicKey}{DATA}},[$rt,$mz,$line];
				$lineData{$xicKey}{PEP}{$peptide}=1;
				$lineData{$xicKey}{SEQ}{$sequence}=1;
			}
			# print '.' unless $. % 5000;
		}
		close IN;
		# print " Done.\n";

		##<Looking for merged ions for sampX
		my $anaID=$refSampIDs->{"samp$msrunIndex"};
		foreach my $xicKey (sort {$lineData{$a}{DATA}[0][0]<=>$lineData{$b}{DATA}[0][0] || $lineData{$a}{DATA}[0][1]<=>$lineData{$b}{DATA}[0][1]} keys %lineData) {
			if (scalar keys %{$lineData{$xicKey}{SEQ}}==1) { # Single or merged isobaric ions with SAME sequence => keep and flag if merged
				if (scalar @{$lineData{$xicKey}{DATA}}==1) { # good ion
					print OUT $lineData{$xicKey}{DATA}[0][2];
					$totalIons++;
				}
				else { # merged isobaric ions
					my (%sameSample,%distinctPtms);
					foreach my $refLine (@{$lineData{$xicKey}{DATA}}) {
						print OUT $refLine->[2];
						my ($quantifMCQID,$group,$msrun,$msrunfile,$mz,$rt,$maxintensity,$area,$rtbegin,$rtend,$peptide,$isotope,$sequence,$z,$mods)=split("\t",$refLine->[2]);
						foreach my $pepEntry (split(/\s/,$mods)) {
							my ($pepIDloc,$zloc,$mz2loc,$pepChannelIDloc)=split(/:/,$pepEntry);
							next if $z != $zloc; # make sure we are comparing the same charge state
							if ($refPep2Ana->{$pepIDloc} == $anaID) {$sameSample{$pepIDloc}=1;} # counts ions only if found in current sample (otherwise -> ghost = 1)
							my $vMod=$refVmods->{$pepIDloc} || 'NONE'; # 'NONE' should never happen because sequence are the same
							push @{$distinctPtms{$vMod}},[$sameSample{$pepIDloc} || 0,$pepIDloc]; # [0/1,pepID]
						}
					}
					my @referencePepIDs;
					foreach my $vMod (keys %distinctPtms) { # record 1 (smallest) pepID per isoform
						my $refPepID=(sort{$b->[0]<=>$a->[0] || $a->[1]<=>$b->[1]} @{$distinctPtms{$vMod}})[0]; # lowest ID of current sample if possible
						push @referencePepIDs,$refPepID->[1];
					}
					$refMergedXICs->{"samp$msrunIndex:$xicKey"}=join(',',sort{$a<=>$b} @referencePepIDs);
					my $numMerged=scalar keys %sameSample || 1; # identified ions or 1 ghost to be created
					$totalMerged+=$numMerged;
					$totalIons+=$numMerged;
				}
			}
			else { # Merged with DIFFERENT sequence => exclude
				# print EXCL "#-------------------\n";
				my %sameSample;
				foreach my $refLine (@{$lineData{$xicKey}{DATA}}) {
					print EXCL $refLine->[2];
					my ($quantifMCQID,$group,$msrun,$msrunfile,$mz,$rt,$maxintensity,$area,$rtbegin,$rtend,$peptide,$isotope,$sequence,$z,$mods)=split("\t",$refLine->[2]);
					foreach my $pepEntry (split(/\s/,$mods)) {
						my ($pepIDloc,$zloc,$mz2loc,$pepChannelIDloc)=split(/:/,$pepEntry);
						next if $z != $zloc; # make sure we are comparing the same charge state
						$sameSample{$pepIDloc}=1 if $refPep2Ana->{$pepIDloc} == $anaID; # counts ions only if found in current sample (otherwise -> ghost = 1)
					}
				}
				my $numExcluded=scalar keys %sameSample; # identified ions or 0 ghost (no ghost created in this case)
				$totalExcluded+=$numExcluded;
				$totalIons+=$numExcluded;
			}
		}
		unlink $sampFile;
	}
	# print " Done: $totalMerged (",(sprintf "%.3f",100*$totalMerged / $totalIons),"%) merged, $totalExcluded (",(sprintf "%.3f",100*$totalExcluded / $totalIons),"%) excluded / $totalIons ions\n";
	close OUT;
	close EXCL;
# system "cp $inFile $inFile.back"; # temp DEBUG
	system "mv $outFile $inFile"; # replace original with clean file

	return ($totalIons,$totalMerged,$totalExcluded);
}


#############################
####<Checking for errors>####
#############################
sub checkForErrors {
	my ($mainDbh)=@_; # can be undef
	if ($errorFile && -s $errorFile) {
		#system "sed -i '/Permanently/ d' $errorFile"; # remove warning from ssh key addition
		my $dbh=$mainDbh || &promsConfig::dbConnect('no_user');
		$dbh->rollback;
		#<Check if error already handled (by parent/child/other) process)
		my ($handled)=$dbh->selectrow_array("SELECT 1 FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID AND STATUS=-2");
		if ($handled) { # already handled: end quietly
			$dbh->disconnect;
			exit;
		}
		#<Handle error
		$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID"); # failed
		$dbh->commit;
		$dbh->disconnect;
		die "Aborting quantification due to errors.";
	}
}


##############################################
####<Checking children processes progress>####
##############################################
sub checkFileParseProgress { # GLOBALS: $numAna, $fileStat, $numSteps, $rootPepFileName
	my ($mcqOutFile1,$prevNumParsed)=@_;
	my $numParsed=$prevNumParsed;
	if ($numAna > 5 && -e $mcqOutFile1) {
		my $numParsed=`grep -c $rootPepFileName $mcqOutFile1`;
		$numParsed=~s/\D//g;
		if ($numParsed > $prevNumParsed) {
			open(FILESTAT,">>$fileStat");
			print FILESTAT "2/$numSteps MassChroQ is parsing XML files ($numParsed/$numAna)\n";
			close FILESTAT;
		}
	}
	return $numParsed;
}
sub checkExtractionProgress { # GLOBALS: $numAna, $fileStat, $numSteps
	my ($mcqOutFile2,$prevMessage)=@_;
	my $message=$prevMessage;
	if (-e $mcqOutFile2) {
		my $progress=`grep -c '^post matching' $mcqOutFile2`;
		$progress=~s/\D//g;
		if ($progress) {$message='4/4: Post matching peaks';}
		else {
			$progress=`grep -c '^Quant.*preparing MS run' $mcqOutFile2`;
			$progress=~s/\D//g;
			if ($progress) {$message="3/4: Quantifying MS runs $progress/$numAna";}
			else {
				$progress=`grep -c 'Aligning MS run' $mcqOutFile2`;
				$progress=~s/\D//g;
				if ($progress) {$message="2/4: Aligning MS runs $progress/$numAna";}
				else {
					$progress=`grep -c '^MS run .*: added' $mcqOutFile2`;
					$progress=~s/\D//g;
					if ($progress) {$message="1/4: Adding XML files $progress/$numAna";}
				}
			}
		}
		if ($message ne $prevMessage) {
			open(FILESTAT,">>$fileStat");
			print FILESTAT "3/$numSteps MassChroQ is performing XIC extraction ($message)\n";
			close FILESTAT;
		}
	}
	return $message;
}
sub checkFilteringProgress { # GLOBALS: $numAna, $fileStat, $currStep, $numSteps
	my ($filterOutFile,$prevMessage)=@_;
	my $message=$prevMessage;
	if (-e $filterOutFile) {
		my $progress=`tail -1 $filterOutFile`;
		chomp($progress);
		$message=$progress if $progress=~/^\d/; # step number
	}
	if ($message ne $prevMessage) {
		open(FILESTAT,">>$fileStat");
		print FILESTAT "$currStep/$numSteps Filtering retention time outliers ($message)\n";
		close FILESTAT;
	}
	return $message;
}

############################ --> Copy of send2Biologist.cgi 06/02/13
####<Check delta length>#### Compares delta lengthes between 2 proteins with Peptide density of the smaller one
############################ (delta is considered not significant if smaller than pep density)
sub deltaLength {
	my ($l_a,$l_b,$d_a,$d_b)=@_;
	my $pepDensVal=($l_b > $l_a)? $d_a : $d_b;
	if (abs($l_b-$l_a) > $pepDensVal) {return $l_a<=>$l_b} else {return 0}
}

###############################################
####> Class that handles the mzXMLHandler <####
###############################################
package mzXMLHandler; {

    my %scan_type_values_array;#MS1 - MS2 (...)
    my %scan_acquisition_time_values_array;
    my %scanAllInfo;

    sub new {#Constructor of the mzXMLHandler
        my ($type)= @_;
        my $self = bless ({}, $type);
        return $self;
    }

    sub start_document {
        my ($self) = @_;
        #print "Starting reading mzXML file: $self->{'mzxmlfile'}<BR>\n";
    }

    sub end_document {
        my ($self) = @_;
        #print "Finishing reading mzXML file<BR>\n";
    }

    sub start_element {# Read an element of the mzXML <-> just scan and peaks event !
        my ($self, $element) = @_;

        if($element->{'Name'} eq "scan" ){
            $self->{'rt'} = "$element->{'Attributes'}{'{}retentionTime'}{'Value'}";
            $self->{'rt'} = substr($self->{'rt'},2,length($self->{'rt'})-3);
            $self->{'scannum'} = $element->{'Attributes'}{'{}num'}{'Value'};
            $self->{'scanlevel'} = $element->{'Attributes'}{'{}msLevel'}{'Value'};
            #print "Retention-Time=$self->{'rt'} ScanNumber=$self->{'scannum'}\n";
        }elsif ($element->{'Name'} eq "precursorMz"){
            $self->{'precursorMz'}=1;
            $self->{'precursorCharge'}=$element->{'Attributes'}{'{}precursorCharge'}{'Value'};
            $self->{'precursorMzData'}="";
        }
    }

    sub end_element {
        my ($self, $element) = @_;
        if($element->{'Name'} eq "precursorMz" ){
            $self->{'precursorMz'}=0;
            $scanAllInfo{$self->{'scannum'}}{'RT'}=$self->{'rt'};
            $scanAllInfo{$self->{'scannum'}}{'precursorMzData'}=$self->{'precursorMzData'};
            $scanAllInfo{$self->{'scannum'}}{'msLevel'}=$self->{'scanlevel'};
            if ($self->{'scanlevel'} == 2) {# Save the charge of the precursor
                $scanAllInfo{$self->{'scannum'}}{'precursorCharge'}=$self->{'precursorCharge'};
            }
        }
    }

    sub characters {
        my ($self, $element) = @_;
        if($self->{'precursorMz'} && $self->{'precursorMz'}==1){
            $self->{'precursorMzData'}.=$element->{'Data'};
        }
    }

    sub get_scan_type_values_array {
        my ($self)=@_;
        return \%scan_type_values_array;
    }

    sub get_acquisition_time_values_array {
        my ($self)=@_;
        return \%scan_acquisition_time_values_array;
    }

    sub get_scanAllInfo {
        my ($self)=@_;
        return \%scanAllInfo;
    }

}


####> Revision history
# 1.9.15 [CHANGE] Increased cluster memory in step-based computation (PP 28/06/21)
# 1.9.14 [CHANGE] Step-based cluster memory computation for both sub jobs (PP 24/06/21)
# 1.9.13 [CHANGE] Higher cluster memory for file parsing (step 1) (PP 21/06/21)
# 1.9.12 [CHANGE] Lower cluster memory estimation for all child jobs (PP 02/06/21)
# 1.9.11 [ENHANCEMENT] Add minimum threshold for peptide RT dispersion filter (VL 15/04/21)
# 1.9.10 [UPDATE] Distinguishable PTM-position isoforms by PTM-encoded peptide sequence and exclusion/flagging of co-eluting ions (PP 26/05/21)
# 1.9.5 [BUGFIX] Fix wrong reference peptide matching to ghost when multiple charge states preventing transfer of PRS probability (PP 12/02/21)
# 1.9.4 [ENHANCEMENT] Progress status for RT outliers filtering (PP 08/02/21)
# 1.9.3 [MERGE] Conflicts solve (PP 04/02/21)
# 1.9.2b [ENHANCEMENT] Added more error check points (PP 03/02/21)
# 1.9.2 [BUGFIX] Remove MassChroQ RT filter on SILAC data (VL 01/02/21)
# 1.9.1 [ENHANCEMENT] Use of fork in case MassChroQ local jobs to allow progress tracking (PP 01/02/21)
# 1.9.0 [ENHANCEMENT] Better cluster memory, error management and progress status display & disabled trace (PP 28/01/21)
# 1.8.6 [MINOR] Add regular commits to prevent timeout from other scripts (PP 07/01/21)
# 1.8.5 [MINOR] Add error management for MassChroQ filter (VL 09/10/20)
# 1.8.4 [BUGFIX] Change R path for MassChroQ filter, it must run on cluster (VL 27/08/20)
# 1.8.3 [CHANGE] Minor changes in cluster job handling (PP 29/07/20)
# 1.8.2 [ENHANCEMENT] Add natural isotopes minimum abundance to MassChroQ parameters (VL 22/07/20)
# 1.8.1 [ENHANCEMENT] Add filtering of MassChroQ output on peptides RT dispersion (VL 22/07/20)
# 1.8.0 [ENHANCEMENT] Optimized SQL queries and error management (PP 19/06/20)
# 1.7.7 [MODIF] Changed createMatchGroup to fit with promsMod.pm prototype (VS 03/07/19)
# 1.7.6 Change to use $cluster{'runJob'} for jobs (GA 12/11/18)
# 1.7.5 Add export LC_ALL="C" in all cluster calls (GA 09/11/18)
# 1.7.4 Also tests for 2.2.12 version when running MassChroQ locally (PP 11/10/18)
# 1.7.3 Customized version for demo/dev while using myproms_1.1.19-1.img (GA 10/10/18)<BR>TODO: remove code when new version of MassChroQ available !!!
# 1.7.2 Change vmod parsing in result file (GA 08/10/18)
# 1.7.1 Minor bug corrections (GA 11/06/18)
# 1.7.0 Peptide quantification data now written to file $promsPath{quantification}/project_$projectID/quanti_$quantifID/peptide_quantification(_$targetPos).txt (PP 11/05/18)
# 1.6.1 Modification of createMatchGroups arguments to include visibility computation (GA 27/03/18)
# 1.5.10 Removed OBSERVATION creation: Handled by manageDesignCondition.cgi (PP 13/02/18)
# 1.5.9 Add masschroqML in addition for output (GA 21/12/17)
# 1.5.8 Modify script so as to die more often (GA 20/12/17)
# 1.5.7 Add a limit on CentOS for max memory used (GA 11/12/17)
# 1.5.6 Move createMatchGroups to promsMod.pm and update number of CPUs (GA 01/12/17)
# 1.5.5 Modif in order to launch jobs on CentOS (GA 26/10/17)
# 1.5.4 Minor modification to avoid multiplication of infinite proteins (GA 11/09/17)
# 1.5.3 Update to MassChroQ 2.2.1 version (GA 18/08/17)
# 1.5.2 Minor modification to avoid warning (GA 02/05/17)
# 1.5.1 Modify $sthInsAP query to avoid insertion of duplicate entries (GA 15/03/17)
# 1.5.0 Distinguish label free and SILAC extraction like PD (GA 06/03/17)
# 1.4.9 Store PRS information for ghost-peptides (GA 26/01/17)
# 1.4.8 Minor modification (MLP 09/01/17)
# 1.4.7 Retrieve masschroqML groups definition removed between 1.4.4 and 1.4.5 script version (GA 22/11/2016)
# 1.4.6 Change to deal with autoincrement in PEPTIDE table (GA 19/10/2016)
# 1.4.5 Fix isotope bug on QSET (GA 11/10/16)
# 1.4.4 Drop DATA_SOURCE information, add OBSERVATION table and keep RT_APEX (GA 03/02/16)
# 1.4.3 Minor modification (GA 22/01/16)
# 1.4.2 Triple-tof updates so as to get scanNumber from mzXML reading (GA 17/11/14)
# 1.4.1 Now RT_BEGIN and RT_END stored and modification of mh computation (GA 02/10/14)
# 1.4.0 Modify quanti.masschroqML file (GA 05/09/14)
# 1.3.9 Add range extraction for all charge states (GA 09/04/14)
# 1.3.8 Uniform ANALYSIS in quantif_info.txt (GA 08/04/14)
# 1.3.7 Minor modification of XIC-Trace<BR>Include -e test for masschroq_status2.txt in wait loop (PP 28/03/14)
# 1.3.6 Add xic-trace extraction and modify mode in isotope extraction (GA 28/03/14)
# 1.3.5 Variable memory & duration time for qsub (PP 27/03/14)
# 1.3.4 Fix isotope alignment (GA 05/03/14)
# 1.3.3 Fix unlabeled channel annotation in QUANTIFICATION.QUANTIF_ANNOT<BR>Fix missing PEPTIDE_QUANTIFICATION.DATA_SOURCE for ghost peptides<BR>Fix extra/missing (or inverted!) label-modifs in PEPTIDE_MODIFICATION for ghost peptides<BR>Prevent ghost peptides inflation when multiple XIC-Isotope extractions are performed (PP 20/02/14)
# 1.3.2 Isotope extraction with MassChroQ (GA 29/01/14)
# 1.3.0 Error management like for runDesignQuantification.cgi (PP 09/09/13)
# 1.2.9 Modification in &createMatchGroups function so as to work with virtual proteins (GA 05/09/13)
# 1.2.8 Add option to use cluster-server for extractions (GA 23/08/13)
# 1.2.7 Add score=0 for proteins in 2nd loop that retrieve ghost peptides so as to avoid warnings in &createMatchGroups (GA 19/08/13)
# 1.2.6 Cancel error redirection to masschroq_statusN.txt (PP 25/07/13)
# 1.2.5 Redirect also error (12/07/13)
# 1.2.4 Added missing finish to some queries (PP 11/07/13)
# 1.2.3 Add ms2_method (GA 04/07/13)
# 1.2.2 Symbolic links to mzXML files instead of copy (PP 02/07/13)
# 1.2.1 Remove VAR_MOD and use PEPTIDE_MODIFICATION table to create VMODS for Ghost-Peptide (GA 27/05/13)
# 1.2.0 Refine match-group computing for "virtual proteins" (GA 08/03/13)
# 1.0.9 Improve createMatchGroups, change pepVisibility and insert mort information in PEPTIDE_PROTEIN_ATTRIB for ghosts (GA 07/02/13)
# 1.0.8 Added missed cut data for ghost peptides (PP 16/01/13)
# 1.0.7 Minor bug correction in SQL query for inserting ghost-peptides (GA 24/12/12)
# 1.0.6 Change in paring result_quanti_pep.tsv to keep the most intense XIC for repeated XICs<BR>Keep the m/z information for ghost peptides (GA 07/12/12)
# 1.0.5 Add multi charge state XIC extraction (GA 05/12/12)
# 1.0.4 Minor modification: take into accoun mannually validated peptides (GA 28/11/12)
# 1.0.3 Slight modification in mods column -> problem with MassChroQ suppressing spaces and creating new vmods for ghost peptides (GA 20/11/12)
# 1.0.2 add the matchGroup recomputing (GA 15/11/12)
# 1.0.1 New script lo launch XIC extraction with MassChroQ + Handle ghost peptides (GA 07/11/12)
# 1.0.0 New script lo launch XIC extraction with MassChroQ (GA 04/09/12)
