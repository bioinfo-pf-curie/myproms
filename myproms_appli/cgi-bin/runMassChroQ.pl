#!/usr/local/bin/perl -w

################################################################################
# runMassChroQ.pl              1.7.6                                           #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
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
use POSIX qw(strftime);
use promsConfig;
use promsMod;
use File::Copy;
use XML::SAX;
#exit; # DEBUG!!!!!!!!!!!!!!!!!!!

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo('no_user');
my %clusterInfo=&promsConfig::getClusterInfo;#('debian'); # default is 'centos'
my $masschroqPath=($clusterInfo{'on'})? $clusterInfo{'path'}{'masschroq'} : $promsPath{'masschroq'};

####################
####>Parameters<####
####################
my ($quantifID,$quantifDate)=@ARGV;
my $quantifDir="$promsPath{tmp}/quantification/$quantifDate";
my $maxStates=50;

############################
####> Connecting to DB <####
############################
my $dbh=&promsConfig::dbConnect('no_user');

####>Updating quantification status to 0 (running)<####
$dbh->do("UPDATE QUANTIFICATION SET STATUS=0 WHERE ID_QUANTIFICATION=$quantifID");
$dbh->commit;
my $sthPepInfoAll=$dbh->prepare("SELECT ID_PEPTIDE,ELUTION_TIME,PEP_SEQ,MISS_CUT,MR_EXP,MR_CALC,MR_OBS,CHARGE,DATA,SPEC_COUNT FROM PEPTIDE WHERE ID_ANALYSIS=? AND VALID_STATUS>=1");# Avoid to take ghost peptides
my $sthPepInfoGP=$dbh->prepare("SELECT ID_PEPTIDE,PEP_SEQ,MISS_CUT,CHARGE,DATA FROM PEPTIDE WHERE ID_ANALYSIS=? AND VALID_STATUS=0");# Take ghost peptides only
my $sthProtInfo=$dbh->prepare("SELECT PA.ID_PROTEIN,IDENTIFIER,SCORE,DB_RANK,VISIBILITY FROM PROTEIN P,PEPTIDE_PROTEIN_ATTRIB PPA,ANALYSIS_PROTEIN PA WHERE P.ID_PROTEIN=PPA.ID_PROTEIN AND P.ID_PROTEIN=PA.ID_PROTEIN AND PA.ID_ANALYSIS=? AND ID_PEPTIDE=?");
my $sthGetInstrument=$dbh->prepare("SELECT INSTRUMENT FROM ANALYSIS WHERE ID_ANALYSIS=?");
my $sthGetNumMods=$dbh->prepare("SELECT COUNT(*) FROM PEPTIDE_MODIFICATION WHERE ID_MODIFICATION=? AND ID_PEPTIDE=?");
#my $sthCountObs=$dbh->prepare("SELECT COUNT(*) FROM OBSERVATION WHERE ID_ANALYSIS=?");
#my $sthInsObs=$dbh->prepare("INSERT INTO OBSERVATION (ID_ANALYSIS,TARGET_POS) VALUES (?,?)");

#############################
###>Generating data files<###
#############################
my $fileStat="$quantifDir/status_$quantifID.out";
open(FILESTAT,">$fileStat");
print FILESTAT "Started ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
print FILESTAT "1/3 Generating data files\n";
close FILESTAT;

####>Getting the parameters from quantif_info.txt<####
open (INFO,"$quantifDir/quantif_info.txt");
my $section='';
my (%params,%sampIDs,%pepInfo,%vmods,%spectralCount,%seqMissedCuts,%proteinScore,%channels,%isoLists,%isoInfo,%anaLabelModifs);
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
		else{
			$params{$paramName}=$paramValue;
		}
	}
}
close INFO;

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
			$anaLabelModifs{$modID}=1;
		}
		@{$channels{$channelID}}=($name,$labels,$hypTechDeltaMass);
		$lightChannelID=$channelID if $hypTechDeltaMass < 0.00001;
	}
}

####>Writing masschroqML file<####
open (MASSCHROQML,">$quantifDir/quanti.masschroqML"); # item valueR valueDB
print MASSCHROQML "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
print MASSCHROQML "<masschroq>\n";

my ($sampID,$rawString,$pepString,$groupString,$anaRef,$instrument)=(0,'','','','','');
###> Get the right scan-number for QSTAR mass spectrometer...
my $deltaTime=5.0;
my $projectID;
#my (%mrObs,%mrCalc,%dbRank);
my (%dbRank);
foreach my $anaID (keys %{$params{'MZXML'}}) {
	###> Add Observation if not inserted already
	#$sthCountObs->execute($anaID);
	#my ($nbObs)=$sthCountObs->fetchrow_array;
	#if ($nbObs==0) {
	#	$sthInsObs->execute($anaID,$targetPos);
	#}
	$sthGetInstrument->execute($anaID);
	($instrument)=$sthGetInstrument->fetchrow_array;
	$instrument='' unless $instrument;

	my $dFormat =($params{'MZXML'}{$anaID} =~ /mzXML/) ?  'mzxml' : 'mzml';
	$rawString .="<data_file format=\"$dFormat\" id=\"samp$sampID\" path=\"$quantifDir/$params{'MZXML'}{$anaID}\" type=\"$params{'RAWDATA_ACQUISITION'}\"/>\n";
	$groupString .= "samp$sampID ";
	if ($params{REFERENCE} && $anaID == $params{REFERENCE}){# In isotope case... no reference!
		$anaRef="samp$sampID";
	}

	####>Writing pepList file<####
	my $pepListFileName="$quantifDir/pepList_$anaID.txt";
	###> 1st loop: Get only validated peptides (no ghost ones with VALID_STATUS=0)
	$sthPepInfoAll->execute($anaID);

	$projectID=&promsMod::getProjectID($dbh,$anaID,'analysis') unless $projectID;
	symlink("$promsPath{tmp}/upload/project_$projectID/$params{'MZXML'}{$anaID}","$promsPath{tmp}/quantification/$quantifDate/$params{'MZXML'}{$anaID}");
	my ($handler,$xmlparser,$scan_allInfo);
	### Only for QSTAR files (WIFF)
	if ($instrument eq 'ESI-QUAD-TOF'){
		$handler=mzXMLHandler->new("$quantifDir/$params{'MZXML'}{$anaID}");
		$xmlparser = XML::SAX::ParserFactory->parser(Handler => $handler );
		$xmlparser->parse_uri("$quantifDir/$params{'MZXML'}{$anaID}");
		$scan_allInfo=$handler->get_scanAllInfo();
	}

	open(PEPFILE, ">$pepListFileName") || die("Cannot open $pepListFileName");
	print PEPFILE "scan\tsequence\tmh\tz\tproteins\tmods\n";
	my @sortedScanNumbers;
	@sortedScanNumbers=sort{$scan_allInfo->{$a}{'RT'}<=>$scan_allInfo->{$b}{'RT'}} keys %{$scan_allInfo} if ($instrument eq 'ESI-QUAD-TOF');
	while (my ($idPeptide,$elutionTime,$pepSeq,$missCut,$mrExp,$mrCalc,$mrObs,$charge,$data,$sc) = $sthPepInfoAll->fetchrow_array) {
		my $vmod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$idPeptide,$anaID,$pepSeq,\%anaLabelModifs);
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
		if ($elutionTime=~/sc(\d+)/){ $scanNumber=$1; }
		if($elutionTime =~ /to/){
			($timeMin,$timeMax)=split(/ to /,$elutionTime);
		}
		elsif ($scanNumber == -1 && $elutionTime=~/et(\d+\.\d+)/) {
			($timeMin,$timeMax)=($1,$1);
		}
		else{
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

		$sthProtInfo->execute($anaID,$idPeptide);
		# Extract only peptides from top-match PROTEINS but keep the score of the others to know the ones that are not virtual
		my $protName=' ';
		while ( my ($protID,$identifier,$score,$dbR,$visibility)=$sthProtInfo->fetchrow_array ){
			next unless defined($protID);
			$proteinScore{$anaID}{$protID}=$score;
			$dbRank{$identifier}=$dbR;
			$protName = $identifier unless $visibility != 2;
		}
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
		$pepInfo{"ID"}{$idPeptide}=$anaID; # Only for no-ghost peptides
		my $pepChannelID=-1;
		if ($getLabelInfo) {
			###> Check if the peptide contains a label modification and tag it !
			foreach my $modID (keys %isoLists) {
				$sthGetNumMods->execute($modID,$idPeptide);
				my ($numMod)=$sthGetNumMods->fetchrow_array;
				my ($channelID,$deltaMass,$modifRes)=@{$isoLists{$modID}};
				if ($numMod>0 || ($deltaMass==0 && $pepChannelID==-1)) {
					$pepChannelID=$channelID;
				}
			}
		}
		$pepInfo{"$pepSeq;$vmod;$charge;$anaID;$pepChannelID"}{"SEQVMODCHARGE"}=$idPeptide;
		$pepInfo{"$pepSeq;$vmod;$charge;$anaID;$pepChannelID"}{"DATA"}=$data;
		#print PEPFILE "$scanNumber\t$pepSeq\t$mrExp\t$charge\t$protInfo\t$idPeptide;$vmod\n";
		print PEPFILE "$scanNumber\t$pepSeq\t$mrCalc\t$charge\t$protName\t$idPeptide:$mrCalc:$pepChannelID\n";
	}
	close PEPFILE;
	$pepString .="<peptide_file data=\"samp$sampID\" path=\"$pepListFileName\"/>\n";

	###> 2nd loop: Get only ghost peptides (with VALID_STATUS=0)
	$sthPepInfoGP->execute($anaID);
	while (my ($idPeptide,$pepSeq,$missCut,$charge,$data) = $sthPepInfoGP->fetchrow_array) {
		$seqMissedCuts{$pepSeq}=$missCut;
		my $vmod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$idPeptide,$anaID,$pepSeq,\%anaLabelModifs);
		$vmod='' unless $vmod;
		$pepInfo{"ID"}{$idPeptide}=$anaID;
		$spectralCount{$idPeptide}=0;
		my $pepChannelID=-1;
		if ($getLabelInfo) {
			###> Check if the peptide contains a label modification and tag it !
			foreach my $modID (keys %isoLists) {
				$sthGetNumMods->execute($modID,$idPeptide);
				my ($numMod)=$sthGetNumMods->fetchrow_array;
				my ($channelID,$deltaMass,$modifRes)=@{$isoLists{$modID}};
				if ($numMod>0 || ($deltaMass==0 && $pepChannelID==-1)) {
					$pepChannelID=$channelID;
				}
			}
		}
		next if $pepInfo{"$pepSeq;$vmod;$charge;$anaID;$pepChannelID"}{"SEQVMODCHARGE"}; # A no-ghost peptide was already referenced
		$pepInfo{"$pepSeq;$vmod;$charge;$anaID;$pepChannelID"}{"SEQVMODCHARGE"}=$idPeptide; # Only for ghost peptides
		$pepInfo{"$pepSeq;$vmod;$charge;$anaID;$pepChannelID"}{"DATA"}=$data;
		$sthProtInfo->execute($anaID,$idPeptide);
		###> Add score of proteins or virtual proteins - 19/08/13
		while ( my ($protID,$identifier,$score,$dbR,$visibility)=$sthProtInfo->fetchrow_array ){
			next unless defined($protID);
			$score=0 unless $score;
			$proteinScore{$anaID}{$protID}=$score unless defined($proteinScore{$anaID}{$protID});
		}
	}
	$sampIDs{"samp$sampID"}=$anaID;
	$sampID++;
}
$sthGetInstrument->finish;
$sthProtInfo->finish;
$sthPepInfoAll->finish;
$sthPepInfoGP->finish;
$sthGetNumMods->finish;
#$sthInsObs->finish;
#$sthCountObs->finish;
$dbh->commit;
$dbh->disconnect;

###> File info
print MASSCHROQML "<rawdata>\n$rawString</rawdata>\n";
print MASSCHROQML "<groups>\n<group data_ids=\"$groupString\" id=\"G1\"/>\n</groups>\n";
print MASSCHROQML "<peptide_files_list>$pepString</peptide_files_list>\n";

####> Isotope info (get SILAC heavy and light info)
my @labelOrder;
if ($getLabelInfo) {
	print MASSCHROQML "<isotope_label_list>\n";

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

			push @labelOrder, ("ID${chanID1}_to_$chanID2","ID${chanID2}_to_$chanID1");
			print MASSCHROQML $isoPlusString;
			print MASSCHROQML $isoNegString;
		}
	}

	print MASSCHROQML"</isotope_label_list>\n";
}

###> Alignment info
# Alignment is done only for non labelled extraction for the moment
print MASSCHROQML "<alignments>\n";
print MASSCHROQML "<alignment_methods>\n";
my $alignmentID; # my_obiwarp or my_ms2
if ($params{'EXTRACTION_ALGO'} eq 'MS2') {
	$alignmentID="my_ms2";
	print MASSCHROQML "<alignment_method id=\"$alignmentID\">\n";
	print MASSCHROQML "<ms2>\n";
	print MASSCHROQML "<ms2_tendency_halfwindow>$params{MS2_TENDENCY}</ms2_tendency_halfwindow>\n";
	print MASSCHROQML "<ms2_smoothing_halfwindow>$params{MS2_SMOUTHING}</ms2_smoothing_halfwindow>\n" if $params{MS2_SMOUTHING};
	print MASSCHROQML "<ms1_smoothing_halfwindow>$params{MS1_SMOUTHING}</ms1_smoothing_halfwindow>\n" if $params{MS1_SMOUTHING};
	print MASSCHROQML "</ms2>\n";
}
else {
	$alignmentID="my_obiwarp";
	print MASSCHROQML "<alignment_method id=\"$alignmentID\">\n";
	print MASSCHROQML "<obiwarp>\n";
	print MASSCHROQML "<lmat_precision>1</lmat_precision>\n";
	print MASSCHROQML "<mz_start>$params{MZ_ALIGN_RANGE_MIN}</mz_start>\n";
	print MASSCHROQML "<mz_stop>$params{MZ_ALIGN_RANGE_MAX}</mz_stop>\n";
	print MASSCHROQML "</obiwarp>\n";
}
print MASSCHROQML "</alignment_method>\n";
print MASSCHROQML "</alignment_methods>\n";
print MASSCHROQML "<align group_id=\"G1\" method_id=\"$alignmentID\" reference_data_id=\"$anaRef\"/>\n";
print MASSCHROQML "</alignments>\n";


###> Quanti info
print MASSCHROQML "<quantification_methods>\n";
print MASSCHROQML "<quantification_method id=\"quanti1\">\n";
print MASSCHROQML "<xic_extraction xic_type=\"$params{XIC_EXTRACTION_TYPE}\">\n";
print MASSCHROQML "<$params{XIC_RANGE}_range max=\"$params{MZTOL_MAX}\" min=\"$params{MZTOL_MIN}\"/>\n";
print MASSCHROQML "</xic_extraction>\n";
print MASSCHROQML "<xic_filters>\n";
print MASSCHROQML "	<anti_spike half=\"$params{ANTISPIKE}\"/>\n" if $params{ANTISPIKE};
print MASSCHROQML "	<background half_mediane=\"$params{MED_MIN}\" half_min_max=\"$params{MED_MAX}\"/>\n" if $params{MED_MIN};
print MASSCHROQML "	<smoothing half=\"$params{SMOOTH}\"/>\n" if $params{SMOOTH};
print MASSCHROQML "</xic_filters>\n";
print MASSCHROQML "<peak_detection>\n";
print MASSCHROQML "<detection_zivy>\n";
print MASSCHROQML "<mean_filter_half_edge>1</mean_filter_half_edge>\n";
print MASSCHROQML "<minmax_half_edge>3</minmax_half_edge>\n";
print MASSCHROQML "<maxmin_half_edge>2</maxmin_half_edge>\n";
print MASSCHROQML "<detection_threshold_on_max>$params{DT_STOP} </detection_threshold_on_max>\n";
print MASSCHROQML "<detection_threshold_on_min>$params{DT_START} </detection_threshold_on_min>\n";
print MASSCHROQML "</detection_zivy>\n";
print MASSCHROQML "</peak_detection>\n";
print MASSCHROQML "</quantification_method>\n";
print MASSCHROQML "</quantification_methods>\n";

###> Results format
print MASSCHROQML "<quantification>\n";
print MASSCHROQML "<quantification_results>\n";
print MASSCHROQML "<quantification_result format=\"tsv\" output_file=\"$quantifDir/result_quanti\"/>\n";
print MASSCHROQML "<quantification_result format=\"masschroqml\" output_file=\"$quantifDir/result_quanti_mcqml\"/>\n";
print MASSCHROQML "</quantification_results>\n";
####> Test of quantification traces...
if ($params{TRACES}) {
	print MASSCHROQML "<quantification_traces>\n";
	print MASSCHROQML "<all_xics_traces output_dir=\"$quantifDir/all_xics_traces\" format=\"tsv\"/>\n";
	print MASSCHROQML "</quantification_traces>\n";
}

####> End test
print MASSCHROQML "<quantify id=\"q1\" quantification_method_id=\"quanti1\" withingroup=\"G1\">\n";
print MASSCHROQML "<peptides_in_peptide_list mode=\"$params{XIC_VAL}\"/>\n";
print MASSCHROQML "</quantify>\n";

####> SILAC part-start
if ($getLabelInfo) {
	my $labelOrderStg=join(' ', @labelOrder);
	print MASSCHROQML "<quantify id=\"q2\" quantification_method_id=\"quanti1\" withingroup=\"G1\">\n";
	#print MASSCHROQML "<peptides_in_peptide_list mode=\"post_matching\"  isotope_label_refs=\"$labelOrderStg\"/>\n";
	print MASSCHROQML "<peptides_in_peptide_list mode=\"$params{XIC_VAL}\"  isotope_label_refs=\"$labelOrderStg\"/>\n";
	print MASSCHROQML "</quantify>\n";
}
####> SILAC part-end

####>
print MASSCHROQML "</quantification>\n";

print MASSCHROQML "</masschroq>\n";
close MASSCHROQML;

############################
####> Launch MassChroQ <####
############################
open(FILESTAT,">>$fileStat");
print FILESTAT "2/3 Waiting for MassChroQ results...\n";
close FILESTAT;

####>
####> PART TO REMOVE
####>
my $addParamStg='';
if ($clusterInfo{'on'}) {
		my $commandFile="$quantifDir/command.sh";
		my $timeStamp=strftime("%Y%m%d%H%M%S",localtime);
		open (COMM,">$commandFile");
		print COMM "#!/bin/bash\n";
		print COMM "export LC_ALL=\"C\"\n cd $quantifDir\n$masschroqPath/masschroq -v 2>&1\n";
		close COMM;
		my $modBash=0775;
		chmod $modBash,$commandFile;
		my %jobParams=(
maxMem=>'1Gb',
numCPUs=>1,
maxHours=>1,
jobName=>"myProMS_test_$timeStamp",
outFile=>'PBS.txt',
errorFile=>'PBSerror.txt',
jobEndFlag=>"_END_$timeStamp",
noWatch=>1
		);
		$clusterInfo{'runJob'}->($quantifDir,$commandFile,\%jobParams);
		sleep 60;
		while (!(-e "$quantifDir/PBS.txt" && `tail -3 $quantifDir/PBS.txt | grep _END_$timeStamp`)) {
				sleep 5;
		}
		open (PBS,"$quantifDir/PBS.txt");
		while (my $line=<PBS>) {
				chomp($line);
				if ($line =~ /MassChroQ [version ]*([\d\.]+)/) {
						if ($1=~ /2\.2\.12/) {
								$addParamStg="truc ";
						}
						last;
				}
		}
		close PBS;
}
else {
	my @response=`$masschroqPath/masschroq -v 2>&1`;
	$addParamStg="truc " if ($response[0] && $response[0]=~/2\.2\.12/);
}
####>
####> END PART TO REMOVE
####>

###> 1st step: Parse peptide option
if ($clusterInfo{'on'}) {
		#code
		my $scriptfile = "$quantifDir/script1";
		open (BASH,"+>",$scriptfile);
		print BASH qq
|#!/bin/bash
export LC_ALL="C"
cd $quantifDir
$masschroqPath/masschroq -p $addParamStg$quantifDir/quanti.masschroqML > $quantifDir/masschroq_status1.txt
|;
		close BASH;
		my $modBash=0775;
		chmod $modBash, $scriptfile;
		my $timeStamp=strftime("%Y%m%d%H%M%S",localtime);
		my %jobParams=(
maxMem=>'10Gb',
numCPUs=>1,
maxHours=>24,
jobName=>"mcq1_$quantifID",
outFile=>'PBS1.txt',
errorFile=>'PBSerror1.txt',
jobEndFlag=>"_END_$timeStamp",
noWatch=>1
		);
		$clusterInfo{'runJob'}->($quantifDir,$scriptfile,\%jobParams);
		sleep 60;
		while (!(-e "$quantifDir/PBS1.txt" && `tail -3 $quantifDir/PBS1.txt | grep _END_$timeStamp`)) {
				sleep 5;
		}
}
else{
		system "cd $quantifDir; $masschroqPath/masschroq -p $addParamStg$quantifDir/quanti.masschroqML > $quantifDir/masschroq_status1.txt"; # if no cd, this command does not work because no permission to use cgi-bin as temporary directory!
		sleep(5);
}

#print "cd $quantifDir; $promsPath{masschroq}/masschroq -p $quantifDir/quanti.masschroqML > $quantifDir/masschroq_status1.txt\n";
###> Check version
my $goodVersion=0;
open(MCQSTAT,"$quantifDir/masschroq_status1.txt");
while (my $line=<MCQSTAT>) {
		chomp($line);
		if ($line =~ /Running MassChroQ 2\.2/) { #Running MassChroQ 2.2.2 with XML file
				$goodVersion=1;
		}
}
close(MCQSTAT);
unless($goodVersion) {
	die "It seems that you have a version of MassChroQ software not supported (2.2 or 2.2.12 optimized version)\n";
}
###> 2nd step: read and add allChargeStates peptides.
if ($params{'ALLCHARGESTATES'}) {# check-box was checked
	###> Read the parse-peptide XML file
	open(PEPTIDEPARSED,"$quantifDir/parsed-peptides_quanti.masschroqML"); # peptides parsed by masschroq in the 1st step
	open(ALLCHARGE,">$quantifDir/parsed-peptides2_quanti.masschroqML");
	my ($onPeptide,$mh)=0;
	while (my $line=<PEPTIDEPARSED>) {
		print ALLCHARGE $line;
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
					print ALLCHARGE "   <observed_in z=\"$i\" data=\"$sample\" scan=\"$scanNum\"/>\n";
				}
			}
		}
	}
	close PEPTIDEPARSED;
	close ALLCHARGE;
	system ("mv $quantifDir/parsed-peptides2_quanti.masschroqML $quantifDir/parsed-peptides_quanti.masschroqML");
}
my ($massChroqError,$pbsError);
####> Use cluster
#$promsPath{'qsub'}=undef; # DEBUG
#if ($promsPath{'qsub'}) {
if ($clusterInfo{'on'}) {
    ##########################
	####> Script command <####
	##########################
	my ($nbProcs)=($sampID > 19)?8:($sampID > 9)?6:4;
	my $scriptfile = "$quantifDir/script2";
	open (BASH,"+>",$scriptfile);
	print BASH qq
|#!/bin/bash
export LC_ALL="C"
cd $quantifDir
$masschroqPath/masschroq --cpus $nbProcs -t $quantifDir $quantifDir/parsed-peptides_quanti.masschroqML > $quantifDir/masschroq_status2.txt 2> $quantifDir/masschroq_errors.txt
|;
	close BASH;
	my $modBash=0775;
	chmod $modBash, $scriptfile;
	my $numAna=scalar keys %{$params{'MZXML'}};
	my $maxHours=10+12*$numAna;
	my ($maxMem)=(2*$numAna > 100) ? "100Gb" : 2*$numAna.'Gb';
	my $timeStamp=strftime("%Y%m%d%H%M%S",localtime);
	my %jobParams=(
maxMem=>$maxMem,
numCPUs=>$nbProcs,
maxHours=>$maxHours,
jobName=>"mcq2_$quantifID",
outFile=>'PBS2.txt',
errorFile=>'PBSerror2.txt',
jobEndFlag=>"_END_$timeStamp",
		);
	$clusterInfo{'runJob'}->($quantifDir,$scriptfile,\%jobParams);
	sleep 60;

	###> When MassChroQ is finished, it is written DONE in the status_file -> need to check if mcq is finished every 30sec.
	my $nbWhile=0;

	while ((!-e "$quantifDir/masschroq_status2.txt" || !`tail -3 $quantifDir/masschroq_status2.txt | grep DONE`) && (!-e "$quantifDir/PBS2.txt" || !`tail -3 $quantifDir/PBS2.txt | grep _END_$quantifID`) && !$massChroqError && !$pbsError) {
		if ($nbWhile>$maxHours*60*2) {
			die "MassChroQ is taking too long or died before completion";
		}
		sleep 30;
		$massChroqError=`head -5 $quantifDir/masschroq_errors.txt` if -e "$quantifDir/masschroq_errors.txt";
		$pbsError=$clusterInfo{'checkError'}->("$quantifDir/PBSerror2.txt");
		$nbWhile++;
	}
}
else{
	system "cd $quantifDir; $masschroqPath/masschroq -t $quantifDir $quantifDir/parsed-peptides_quanti.masschroqML > $quantifDir/masschroq_status2.txt 2> $quantifDir/masschroq_errors.txt";
	#print "cd $quantifDir; $promsPath{masschroq}/masschroq -t $quantifDir $quantifDir/parsed-peptides_quanti.masschroqML > $quantifDir/masschroq_status2.txt";
	$massChroqError=`head -5 $quantifDir/masschroq_errors.txt` if -e "$quantifDir/masschroq_errors.txt";
}

####>ERROR Management<####
if ($massChroqError && $massChroqError !~ /WARNING/) {
	die "MassChroQ has generated the following Error:\n$massChroqError";
}
elsif ($pbsError && $pbsError !~ /WARNING/) {
	die "The computation cluster has generated the following Error:\n$pbsError";
}
elsif (! -d $quantifDir) {
	die "The job has been manually interrupted\n";
}

#########################################################
####> Fetching MassChroQ quantification information <####
#########################################################
open(FILESTAT,">>$fileStat");
print FILESTAT "3/3 Storing MassChroQ XIC data\n";
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

#my $sthInsPepQ=$dbh->prepare("INSERT INTO PEPTIDE_QUANTIFICATION (ID_QUANTIF_PARAMETER,ID_PEPTIDE,ID_QUANTIFICATION,QUANTIF_VALUE,TARGET_POS) VALUES (?,?,$quantifID,?,?)");
my $sthUpPep=$dbh->prepare("UPDATE PEPTIDE SET DATA=? WHERE ID_PEPTIDE=?");

my (%quantifPeptides,%rt); # Peptides inserted into PEPTIDE_QUANTIFICATION (identified in an analysis and valid)
my %newPeptides; # Potential gost-peptides.

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
}

#open(PEPQUANTI,"$quantifDir/result_quanti_pep.tsv"); # former file of MassChroQ that was keeping q1 and q2 quantifications
open(PEPQUANTI1,"$quantifDir/result_quanti.d/peptides_q1_G1.tsv"); # label-free area extraction data
while (my $line=<PEPQUANTI1>) {
        next if $line =~ /group/; #skip first line
        chomp($line);
        # example: q1,G1,samp0,F5849FD,518.25061,1249.1633,409292.7,7127488.8,1201.6282,1264.3375,pep855,isotopeName,GWDVGVAGMK,2,613281
        my ($quantifiMCQID,$group,$mrsrun,$msrunfile,$mz,$rt,$maxintensity,$area,$rtbegin,$rtend,$peptide,$isotope,$sequence,$z,$mods)=split("\t",$line);
		
		# Be careful for VMOD choice: some isobaric peptides can be grouped together in $mods string with different positions or from different samples
		# Therefore, to avoid the creation of Ghost-Peptides with incorrect VMOD for this specific entry, 2 priority rules to follow :
		# 1st: choose the vmod from the same msrun peptide (and the highest number of spectral-counts if many possibilities)
		# 2nd: choose the vmod from the closest RT-time provided compared to the rt time available of other experiment
		my ($pepID,$mz2,$pepChannelID);
		my $sameSample=0;
		foreach my $pepEntry (split(/\s/,$mods)) { # each pepEntry looks like this: $idPeptide:$mz:$targetPos
				my ($pepIDloc,$mz2loc,$pepChannelIDloc)=split(/:/,$pepEntry);
				if ($pepInfo{'ID'}{$pepIDloc} == $sampIDs{$mrsrun}) {
						if ($pepID) {
								if ($sameSample) {
										if ($spectralCount{$pepIDloc} > $spectralCount{$pepID}) {
												($pepID,$mz2,$pepChannelID)=($pepIDloc,$mz2loc,$pepChannelIDloc);
										}
								}
								else{
										($pepID,$mz2,$pepChannelID)=($pepIDloc,$mz2loc,$pepChannelIDloc);
								}
						}
						else{
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
				else{
						($pepID,$mz2,$pepChannelID)=($pepIDloc,$mz2loc,$pepChannelIDloc);
				}
		}
		my $vmod=$vmods{$pepID};
		$vmod='' unless $vmod;
		my $mzInt= int $mz;
		if (!$isotope) {
			if ( $pepInfo{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"SEQVMODCHARGE"} ){
				if ( $quantifPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}) {
					if ($area > $quantifPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}) {
						$quantifPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}=$area;
						$rt{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"BEGIN"}=$rtbegin;
						$rt{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"END"}=$rtend;
						$rt{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"APEX"}=$rt;
						$rt{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"TRACE"}=$traceFiles{"q1:$group:$mrsrun:$peptide:$mzInt:$z"} if $params{TRACES};
					}
				}
				else {
					$quantifPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}=$area;
					$rt{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"BEGIN"}=$rtbegin;
					$rt{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"END"}=$rtend;
					$rt{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"APEX"}=$rt;
					$rt{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"TRACE"}=$traceFiles{"q1:$group:$mrsrun:$peptide:$mzInt:$z"} if $params{TRACES};
				}
			}
			else{ # Ghost peptide... has to be added
				if ($newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"AREA"}) {
					if ($area > $newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"AREA"}) {
						$newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"AREA"}=$area;
						$newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"BEGIN"}=$rtbegin;
						$newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"END"}=$rtend;
						$newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"APEX"}=$rt;
						$newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"TRACE"}=$traceFiles{"q1:$group:$mrsrun:$peptide:$mzInt:$z"} if $params{TRACES};
					}
				}
				else {
					$newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"AREA"}=$area;
					$newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"RT"}=$rt/60;
					$newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"MZ"}=$mz;
					$newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"MREXP"}=$mz*$z-1.007825;
					$newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"REFPEPID"}=$pepID;
					$newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"BEGIN"}=$rtbegin;
					$newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"END"}=$rtend;
					$newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"APEX"}=$rt;
					$newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"TRACE"}=$traceFiles{"q1:$group:$mrsrun:$peptide:$mzInt:$z"} if $params{TRACES};
					my $data=$pepInfo{"$sequence;$vmod;$z;$pepInfo{'ID'}{$pepID};$pepChannelID"}{'DATA'};
				    my ($PRS) = ($data =~ /(##PRS=[^##]+)/) if $data;
					$PRS='' unless $PRS;
					$newPeptides{"$sequence;$vmod;$z;$sampIDs{$mrsrun};$pepChannelID"}{"DATA"}=$PRS; # in case there is PhosphoRS information
				}
			}
		}

		#if ( $getLabelInfo ) { # Analysis made to get labeling info
		#		my $traceFile=($params{TRACES})? $traceFiles{"q2:$group:$mrsrun:$peptide:$mzInt:$z"} : "";
		#		@{$isoInfo{$sampIDs{$mrsrun}}{"$sequence;$vmod;$z"}{$mz}}=($area,$rt/60,$mods,$isotope,$rtbegin,$rtend,$traceFile);
		#}
}
close PEPQUANTI1;

if ( $getLabelInfo ) {
		open(PEPQUANTI2,"$quantifDir/result_quanti.d/peptides_q2_G1.tsv"); # isotope area extraction data
		while (my $line=<PEPQUANTI2>) {
				next if $line =~ /group/; #skip first line
				chomp($line);
				# example: q2,G1,samp0,F5849FD,518.25061,1249.1633,409292.7,7127488.8,1201.6282,1264.3375,pep855,isotopeName,GWDVGVAGMK,2,613281
				my ($quantifiMCQID,$group,$mrsrun,$msrunfile,$mz,$rt,$maxintensity,$area,$rtbegin,$rtend,$peptide,$isotope,$sequence,$z,$mods)=split("\t",$line);
				my ($pepID,$mz2,$pepChannelID);
				my $sameSample=0;
				foreach my $pepEntry (split(/\s/,$mods)) { # each pepEntry looks like this: $idPeptide:$mz:$targetPos
						my ($pepIDloc,$mz2loc,$pepChannelIDloc)=split(/:/,$pepEntry);
						if ($pepInfo{'ID'}{$pepIDloc} == $sampIDs{$mrsrun}) {
								if ($pepID) {
										if ($sameSample) {
												if ($spectralCount{$pepIDloc} > $spectralCount{$pepID}) {
														($pepID,$mz2,$pepChannelID)=($pepIDloc,$mz2loc,$pepChannelIDloc);
												}
										}
										else{
												($pepID,$mz2,$pepChannelID)=($pepIDloc,$mz2loc,$pepChannelIDloc);
										}
								}
								else{
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
						else{
								($pepID,$mz2,$pepChannelID)=($pepIDloc,$mz2loc,$pepChannelIDloc);
						}
				}
				my $vmod=$vmods{$pepID};
				$vmod='' unless $vmod;
				my $mzInt= int $mz;

				my $traceFile=($params{TRACES})? $traceFiles{"q2:$group:$mrsrun:$peptide:$mzInt:$z"} : "";
				@{$isoInfo{$sampIDs{$mrsrun}}{"$sequence;$vmod;$z"}{$mz}}=($area,$rt/60,$mods,$isotope,$rtbegin,$rtend,$traceFile);
		}
		close(PEPQUANTI2);
}

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
			my ($pepID,$mh,$channelIDOriginal,@pepIDs)=split(/ |:/,$pepIDmhchannelID);
			if($isotope){
				$isotope=~s/ID//;
				my ($fromChanID,$toChanID)=split('_to_',$isotope);
				if ( $fromChanID eq $channelIDOriginal && $channels{$toChanID}) {
					$qsetValues{$qsetNum}{$toChanID}{'AREA'}=$area;
					$qsetValues{$qsetNum}{$toChanID}{'REFPEPID'}=$pepID;
					$qsetValues{$qsetNum}{$toChanID}{'RT'}=$rt;
					$qsetValues{$qsetNum}{$toChanID}{'MZ'}=$mz;
					$qsetValues{$qsetNum}{$toChanID}{'Z'}=$charge;
					$qsetValues{$qsetNum}{$toChanID}{'BEGIN'}=$rtbegin;
					$qsetValues{$qsetNum}{$toChanID}{'END'}=$rtend;
					$qsetValues{$qsetNum}{$toChanID}{'TRACE'}=$traceFile if $params{TRACES};
					if ($pepInfo{"$seq;$vmod;$charge;$anaID;$toChanID"}{"SEQVMODCHARGE"} && !$qsetValues{$qsetNum}{$toChanID}{'PEP_ID'}) { # This test prevent inflation of ghost peptides !
						$qsetValues{$qsetNum}{$toChanID}{'PEP_ID'}=$pepInfo{"$seq;$vmod;$charge;$anaID;$toChanID"}{"SEQVMODCHARGE"};
						if (!$quantifPeptides{"$seq;$vmod;$charge;$anaID;$toChanID"} ) { # Apply only for ghost peptides that were not quantified
							$quantifPeptides{"$seq;$vmod;$charge;$anaID;$toChanID"}=$area;
							$rt{"$seq;$vmod;$charge;$anaID;$toChanID"}{"BEGIN"}=$rtbegin;
							$rt{"$seq;$vmod;$charge;$anaID;$toChanID"}{"END"}=$rtend;
							$rt{"$seq;$vmod;$charge;$anaID;$toChanID"}{"APEX"}=$rt*60;
							$rt{"$seq;$vmod;$charge;$anaID;$toChanID"}{"TRACE"}=$traceFile if $params{TRACES};
						}
					}

				}
			}
		}
		###> 2nd round: get alignment retrieved peptides.
		foreach my $mz (sort {$a <=> $b} keys %{$isoInfo{$anaID}{$seqVmodCharge}}) {
			my ($area,$rt,$pepIDmhchannelID,$isotope,$rtbegin,$rtend,$traceFile)=@{$isoInfo{$anaID}{$seqVmodCharge}{$mz}};
			#print "AREA=$area RT=$rt MODS=$pepIDmhchannelID ISOTOPE=$isotope\n";
			my ($pepID,$mh,$channelIDOriginal,@pepIDs)=split(/ |:/,$pepIDmhchannelID);
			if(!$isotope) { # If isotope column empty -> original version of the isotope $channelIDOriginal
				#next if $pepInfo{"ID"}{$pepID} && $pepInfo{"ID"}{$pepID} != $anaID;
				$qsetValues{$qsetNum}{$channelIDOriginal}{'AREA'}=$area;
				if ($pepInfo{"$seq;$vmod;$charge;$anaID;$channelIDOriginal"}{"SEQVMODCHARGE"}) {
					$qsetValues{$qsetNum}{$channelIDOriginal}{'PEP_ID'}=$pepID if $pepInfo{"ID"}{$pepID} == $anaID; # 1st priority: get the real area from the peptide of this analysis
					$qsetValues{$qsetNum}{$channelIDOriginal}{'PEP_ID'}=$pepID if !$qsetValues{$qsetNum}{$channelIDOriginal}{'PEP_ID'}; # 2nd: get the area retrieved by alignment ONLY if the peptide does not exist for this analysis!
				}else{
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
			if ( !$qsetValues{$qsetNum}{$qChannelID}{'PEP_ID'} ) {
				$newPeptides{"$seq;$vmod;$charge;$anaID;$qChannelID"}{"AREA"}=$qsetValues{$qsetNum}{$qChannelID}{'AREA'};
				$newPeptides{"$seq;$vmod;$charge;$anaID;$qChannelID"}{"RT"}=$qsetValues{$qsetNum}{$qChannelID}{'RT'};
				$newPeptides{"$seq;$vmod;$charge;$anaID;$qChannelID"}{"MZ"}=$qsetValues{$qsetNum}{$qChannelID}{'MZ'};
				$newPeptides{"$seq;$vmod;$charge;$anaID;$qChannelID"}{"MREXP"}=$qsetValues{$qsetNum}{$qChannelID}{'MZ'}*$qsetValues{$qsetNum}{$qChannelID}{'Z'}-1.007825;
				$newPeptides{"$seq;$vmod;$charge;$anaID;$qChannelID"}{"REFPEPID"}=$qsetValues{$qsetNum}{$qChannelID}{'REFPEPID'};
				my $dataref=$pepInfo{"$seq;$vmod;$charge;$pepInfo{'ID'}{$qsetValues{$qsetNum}{$qChannelID}{'REFPEPID'}};$qChannelID"}{'DATA'};
				my ($PRS) = ($dataref =~ /(##PRS=[^##]+)/) if $dataref;
				$PRS='' unless $PRS;
				$newPeptides{"$seq;$vmod;$charge;$anaID;$qChannelID"}{"DATA"}="$PRS$data";
				$newPeptides{"$seq;$vmod;$charge;$anaID;$qChannelID"}{"BEGIN"}=$qsetValues{$qsetNum}{$qChannelID}{'BEGIN'};
				$newPeptides{"$seq;$vmod;$charge;$anaID;$qChannelID"}{"END"}=$qsetValues{$qsetNum}{$qChannelID}{'END'};
				$newPeptides{"$seq;$vmod;$charge;$anaID;$qChannelID"}{"APEX"}=$qsetValues{$qsetNum}{$qChannelID}{'RT'};
				$newPeptides{"$seq;$vmod;$charge;$anaID;$qChannelID"}{"TRACE"}=$qsetValues{$qsetNum}{$qChannelID}{'TRACE'} if $params{TRACES};
			}else{
				$pepInfo{"$seq;$vmod;$charge;$anaID;$qChannelID"}{"DATA"}.=$data;
			}
			#$pepInfo{"$seq;$vmod;$charge;$anaID"}{"CHANNEL"}=$qChannelID;
			#my $pepID=($qsetValues{$qsetNum}{$qChannelID}{'PEP_ID'})?$qsetValues{$qsetNum}{$qChannelID}{'PEP_ID'}: 'Ghost-Peptide';
			#print "CHANNEL=$qChannelID VALUE=$area PEP_ID=$pepID\n";
		}
	}
}




######################
####>Moving files<####
######################
system "mkdir $promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
system "mkdir $promsPath{quantification}/project_$projectID/quanti_$quantifID";
system "cp $quantifDir/*.txt $promsPath{quantification}/project_$projectID/quanti_$quantifID";# peplist files
system "mv $quantifDir/result_quanti.d $promsPath{quantification}/project_$projectID/quanti_$quantifID";# tsv output files: result_quanti_pep.tsv, result_quanti_prot.tsv and result_quanti_compar.tsv
system "cp $quantifDir/*.time $promsPath{quantification}/project_$projectID/quanti_$quantifID" if glob "$quantifDir/*.time"; # Alignment files -> important to perform again the search
system "cp $quantifDir/*.masschroqML $promsPath{quantification}/project_$projectID/quanti_$quantifID";
system "mv $quantifDir/all_xics_traces $promsPath{quantification}/project_$projectID/quanti_$quantifID" if $params{TRACES};

my %pepQuantHandle;
if ($getLabelInfo) { # labeled quanti
	foreach my $channelID (sort{$a<=>$b} keys %channels) {
		open($pepQuantHandle{$channelID},">$promsPath{quantification}/project_$projectID/quanti_$quantifID/peptide_quantification_$channelID.txt");
		print {$pepQuantHandle{$channelID}} "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n";
	}
}
else { # label-free
	open($pepQuantHandle{0},">$promsPath{quantification}/project_$projectID/quanti_$quantifID/peptide_quantification.txt");
	print {$pepQuantHandle{0}} "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n";
}

open (TRACEFILE,">$promsPath{quantification}/project_$projectID/quanti_$quantifID/pepID_to_traces.txt");
print TRACEFILE "PEPID\tTRACE\n";

#$sthInsPep->finish;
foreach my $seqVmodZanaIDChanID (keys %quantifPeptides ) {
	my ($seq,$vmod,$charge,$anaID,$channelID)=split(/;/,$seqVmodZanaIDChanID);
	#my $targetPos=($channelID > 0)? $channelID: undef;
	#$sthInsPepQ->execute($qparamsID{$areaCode},$pepInfo{$seqVmodZanaIDChanID}{"SEQVMODCHARGE"},$quantifPeptides{$seqVmodZanaIDChanID},$targetPos);
	#$sthInsPepQ->execute($qparamsID{'RT_BEGIN'},$pepInfo{$seqVmodZanaIDChanID}{"SEQVMODCHARGE"},$rt{$seqVmodZanaIDChanID}{"BEGIN"},$targetPos) if $qparamsID{'RT_BEGIN'};
	#$sthInsPepQ->execute($qparamsID{'RT_END'},$pepInfo{$seqVmodZanaIDChanID}{"SEQVMODCHARGE"},$rt{$seqVmodZanaIDChanID}{"END"},$targetPos) if $qparamsID{'RT_END'};
	#$sthInsPepQ->execute($qparamsID{'RT_APEX'},$pepInfo{$seqVmodZanaIDChanID}{"SEQVMODCHARGE"},$rt{$seqVmodZanaIDChanID}{"APEX"},$targetPos);
	my $targetPos=($getLabelInfo)? $channelID : 0;
	my $dataStrg="$qparamsID{$areaCode}\t$pepInfo{$seqVmodZanaIDChanID}{SEQVMODCHARGE}\t$quantifPeptides{$seqVmodZanaIDChanID}\n";
	$dataStrg.="$qparamsID{RT_BEGIN}\t$pepInfo{$seqVmodZanaIDChanID}{SEQVMODCHARGE}\t$rt{$seqVmodZanaIDChanID}{BEGIN}\n" if $qparamsID{'RT_BEGIN'};
	$dataStrg.="$qparamsID{RT_END}\t$pepInfo{$seqVmodZanaIDChanID}{SEQVMODCHARGE}\t$rt{$seqVmodZanaIDChanID}{END}\n" if $qparamsID{'RT_END'};
	$dataStrg.="$qparamsID{RT_APEX}\t$pepInfo{$seqVmodZanaIDChanID}{SEQVMODCHARGE}\t$rt{$seqVmodZanaIDChanID}{APEX}\n";
	print {$pepQuantHandle{$targetPos}} $dataStrg;

	$sthUpPep->execute($pepInfo{$seqVmodZanaIDChanID}{"DATA"},$pepInfo{$seqVmodZanaIDChanID}{"SEQVMODCHARGE"}) if $getLabelInfo;
	print TRACEFILE "$pepInfo{$seqVmodZanaIDChanID}{\"SEQVMODCHARGE\"}\t$rt{$seqVmodZanaIDChanID}{\"TRACE\"}\n" if $params{TRACES};
	#$dbh->commit;
}

#############################################################
###> Queries for insertion of ghost peptides information <###
#############################################################
#my $sthInsPPA=$dbh->prepare("INSERT INTO PEPTIDE_PROTEIN_ATTRIB (ID_PROTEIN,ID_PEPTIDE,ID_ANALYSIS,PEP_BEG,PEP_END,FLANKING_AA) VALUES (?,?,?,?,?,?)");
my $sthAddPA=$dbh->prepare("INSERT INTO PEPTIDE_PROTEIN_ATTRIB (ID_PEPTIDE,ID_PROTEIN,PEP_BEG,PEP_END,FLANKING_AA,IS_SPECIFIC,ID_ANALYSIS) VALUES (?,?,?,?,?,?,?)");
my $sthGetPPA=$dbh->prepare("SELECT PEPTIDE_PROTEIN_ATTRIB.ID_PROTEIN,PEP_BEG,PEP_END,FLANKING_AA,IS_SPECIFIC,IDENTIFIER FROM PEPTIDE_PROTEIN_ATTRIB,PROTEIN WHERE PROTEIN.ID_PROTEIN=PEPTIDE_PROTEIN_ATTRIB.ID_PROTEIN AND ID_PEPTIDE=?");
my $sthInsGhostPep=$dbh->prepare("INSERT INTO PEPTIDE (ID_ANALYSIS,PEP_SEQ,PEP_LENGTH,MISS_CUT,MR_EXP,MR_OBS,CHARGE,ELUTION_TIME,VALID_STATUS,DATA) VALUES (?,?,?,?,?,?,?,?,0,?)"); # VALID_STATUS=0 and SCORE IS NULL characterize ghost peptides
my $queryGetPepM="SELECT ID_MODIFICATION,MODIF_TYPE,POS_STRING,REF_POS_STRING FROM PEPTIDE_MODIFICATION WHERE ID_PEPTIDE=?";
$queryGetPepM.=' AND ID_MODIFICATION NOT IN ('.join(',',keys %anaLabelModifs).')' if $getLabelInfo;
my $sthGetPepM=$dbh->prepare($queryGetPepM);
my $sthInsPepM=$dbh->prepare("INSERT INTO PEPTIDE_MODIFICATION (ID_MODIFICATION,ID_PEPTIDE,MODIF_TYPE,POS_STRING,REF_POS_STRING) VALUES (?,?,?,?,?)");
my $sthCountPepSpec=$dbh->prepare("SELECT COUNT(PEP_SPECIFICITY) FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=? AND ID_PROTEIN=?");
my %virtualProteins;
foreach my $seqVmodZanaIDChanID ( keys %newPeptides ) {
	next if $quantifPeptides{$seqVmodZanaIDChanID}; # QUANTIF_INFORMATION already added
	my ($seq,$vmod,$charge,$anaID,$channelID)=split(/;/,$seqVmodZanaIDChanID);
	$vmod='' unless $vmod;
	###> 1st: New instance in PEPTIDE
	my $et=sprintf "rt%.2f",$newPeptides{$seqVmodZanaIDChanID}{"RT"};
	$sthInsGhostPep->execute($anaID,$seq,length($seq),$seqMissedCuts{$seq},$newPeptides{$seqVmodZanaIDChanID}{"MREXP"},$newPeptides{$seqVmodZanaIDChanID}{"MZ"},$charge,$et,$newPeptides{$seqVmodZanaIDChanID}{"DATA"});
	my $newpepID=$dbh->last_insert_id(undef,undef,'PEPTIDE','ID_PEPTIDE');
	#$sthUpPep->execute($newPeptides{$seqVmodZanaIDChanID}{"DATA"},$newpepID) if $getLabelInfo; # just for isotope-tagged peptides!!!
	###> 2nd : Add VMOD information
	$sthGetPepM->execute($newPeptides{$seqVmodZanaIDChanID}{'REFPEPID'});
	while (my ($modID,$modifType,$posString,$refPosString) = $sthGetPepM->fetchrow_array) {
		$sthInsPepM->execute($modID,$newpepID,$modifType,$posString,$refPosString);
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
	$sthGetPPA->execute($newPeptides{$seqVmodZanaIDChanID}{'REFPEPID'});
	while (my ($protID,$pepBeg,$pepEnd,$flankingAA,$isSpecific,$identifier) = $sthGetPPA->fetchrow_array) {
		#$sthInsPPA->execute($protID,$newpepID,$anaID,$pepBeg,$pepEnd,$flankAA);
		#$isSpecific='NULL' unless defined($isSpecific);
		$sthAddPA->execute($newpepID,$protID,-$pepBeg,-$pepEnd,$flankingAA,$isSpecific,$anaID);
		#$dbh->commit;
		#$sthGetNumProt->execute($anaID,$protID);
		#my ($numID)=$sthGetNumProt->fetchrow_array;
		#$sthInsAP->execute($anaID,$protID) if $numID==0;
		#$dbh->commit;
		if (!defined($proteinScore{$anaID}{$protID})) {
			$proteinScore{$anaID}{$protID}=0;
			$sthCountPepSpec->execute($anaID,$protID);
			my $nbSpec=$sthCountPepSpec->fetchrow_array;
			if ($nbSpec<1) {
				$virtualProteins{$anaID}{$identifier}=1;
			}
		}
	}
	###> 4th: Add PEPTIDE_QUANTIFICATION quantif-info
	#my $targetPos=($channelID >0 )? $channelID : undef;
	#$sthInsPepQ->execute($qparamsID{$areaCode},$newpepID,$newPeptides{$seqVmodZanaIDChanID}{"AREA"},$targetPos);
	#$sthInsPepQ->execute($qparamsID{"RT_BEGIN"},$newpepID,$newPeptides{$seqVmodZanaIDChanID}{"BEGIN"},$targetPos) if $qparamsID{"RT_BEGIN"};
	#$sthInsPepQ->execute($qparamsID{"RT_END"},$newpepID,$newPeptides{$seqVmodZanaIDChanID}{"END"},$targetPos) if $qparamsID{"RT_END"};
	#$sthInsPepQ->execute($qparamsID{"RT_APEX"},$newpepID,$newPeptides{$seqVmodZanaIDChanID}{"APEX"},$targetPos);
	my $targetPos=($getLabelInfo)? $channelID : 0;
	my $dataStrg="$qparamsID{$areaCode}\t$newpepID\t$newPeptides{$seqVmodZanaIDChanID}{AREA}\n";
	$dataStrg.="$qparamsID{RT_BEGIN}\t$newpepID\t$newPeptides{$seqVmodZanaIDChanID}{BEGIN}\n" if $qparamsID{'RT_BEGIN'};
	$dataStrg.="$qparamsID{RT_END}\t$newpepID\t$newPeptides{$seqVmodZanaIDChanID}{END}\n" if $qparamsID{'RT_END'};
	$dataStrg.="$qparamsID{RT_APEX}\t$newpepID\t$newPeptides{$seqVmodZanaIDChanID}{APEX}\n";
	print {$pepQuantHandle{$targetPos}} $dataStrg;

	$sthUpPep->execute($newPeptides{$seqVmodZanaIDChanID}{"DATA"},$newpepID) if $getLabelInfo;
	print TRACEFILE "$newpepID\t$newPeptides{$seqVmodZanaIDChanID}{\"TRACE\"}\n" if $params{TRACES};
	#$dbh->commit;
	###> 5th: Update quantifPeptides hash information
	$quantifPeptides{"$seqVmodZanaIDChanID"}=1;
}
foreach my $targetPos (keys %pepQuantHandle) {close $pepQuantHandle{$targetPos};}
close(TRACEFILE);
#$sthInsPPA->finish;
$sthAddPA->finish;
$sthGetPPA->finish;
#$sthGetNumProt->finish;
#$sthInsAP->finish;
#$sthInsPepQ->finish;
$sthInsGhostPep->finish;
$sthGetPepM->finish;
$sthInsPepM->finish;
$sthCountPepSpec->finish;
$sthUpPep->finish;
$dbh->commit;

my $sthgetMatchInfo=$dbh->prepare("SELECT PEPTIDE.ID_PEPTIDE,PEP_SEQ,PROTEIN.ID_PROTEIN,IDENTIFIER,PROT_LENGTH FROM PEPTIDE,PEPTIDE_PROTEIN_ATTRIB,PROTEIN WHERE PEPTIDE.ID_PEPTIDE=PEPTIDE_PROTEIN_ATTRIB.ID_PEPTIDE AND PEPTIDE_PROTEIN_ATTRIB.ID_PROTEIN=PROTEIN.ID_PROTEIN AND PEPTIDE.ID_ANALYSIS=?");
#my $sthUpAP=$dbh->prepare("INSERT INTO ANALYSIS_PROTEIN SET VISIBILITY=?,MATCH_GROUP=?,NUM_PEP=?,NUM_MATCH=?,PEP_SPECIFICITY=?,PEP_COVERAGE=? WHERE ID_PROTEIN=? AND ID_ANALYSIS=?");
#my $sthGetPepInfo=$dbh->prepare("SELECT PEP_BEG,PEP_END FROM PEPTIDE,PEPTIDE_PROTEIN_ATTRIB WHERE PEPTIDE_PROTEIN_ATTRIB.ID_PEPTIDE=PEPTIDE.ID_PEPTIDE AND VALID_STATUS>=1 AND PEPTIDE_PROTEIN_ATTRIB.ID_ANALYSIS=? AND ID_PROTEIN=? ORDER BY PEP_BEG ASC, PEP_END ASC");
my $sthGetPepInfo=$dbh->prepare("SELECT ABS(PEP_BEG),ABS(PEP_END),PEP_SEQ,P.ID_PEPTIDE FROM PEPTIDE P,PEPTIDE_PROTEIN_ATTRIB PPA WHERE PPA.ID_PEPTIDE=P.ID_PEPTIDE AND PPA.ID_ANALYSIS=? AND ID_PROTEIN=? ORDER BY PEP_BEG ASC, PEP_END ASC"); # Get PEPTIDE and GHOST-PEPTIDE to compute a well coverage on the protein
my $sthGetAPInfo=$dbh->prepare("SELECT IDENTIFIER,MATCH_GROUP,VISIBILITY FROM ANALYSIS_PROTEIN AP,PROTEIN P WHERE AP.ID_PROTEIN=P.ID_PROTEIN AND ID_ANALYSIS=?");
my ($protVisibility)=$dbh->selectrow_array("SELECT PROT_VISIBILITY FROM PROJECT WHERE ID_PROJECT=$projectID");
$protVisibility=0 unless $protVisibility;
#my $sthInsAP=$dbh->prepare("INSERT IGNORE INTO ANALYSIS_PROTEIN (ID_ANALYSIS,ID_PROTEIN,DB_RANK,CONF_LEVEL,SCORE,NUM_PEP,NUM_MATCH,PEP_COVERAGE,PEP_SPECIFICITY) VALUES (?,?,?,0,0,?,?,?,?)");
my $sthInsAP=$dbh->prepare("INSERT IGNORE INTO ANALYSIS_PROTEIN (ID_ANALYSIS,ID_PROTEIN,DB_RANK,CONF_LEVEL,SCORE,NUM_PEP,NUM_MATCH,PEP_COVERAGE,PEP_SPECIFICITY) VALUES (?,?,?,0,0,0,0,?,?)"); # only for adding virtual proteins
my $sthUpMG=$dbh->prepare("UPDATE ANALYSIS_PROTEIN SET VISIBILITY=?,MATCH_GROUP=? WHERE ID_ANALYSIS=? AND ID_PROTEIN=?");
##########################################
###> Queries for update of matchGroup <###
##########################################
foreach my $anaID (keys %{$params{'MZXML'}}) {
	my (%matchList,%matchGroup,%visibility,%identifiers,%proteinLength,%boundaryStatus,%pepSpecificity,%pepProteins,%maxCovLength,%numProtTop,%bestProtVis,%numProtPeptides,%proteinPepDens);
	$sthgetMatchInfo->execute($anaID);

	while(my ($pepID,$seq,$protID,$identifier,$protLen)=$sthgetMatchInfo->fetchrow_array) {
		my $vmod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$anaID,$seq,\%anaLabelModifs);
		$vmod='' unless $vmod;
		$matchList{$identifier}{"$seq$vmod"}=1;
		$pepProteins{$pepID}{$identifier}=1;
		$identifiers{$identifier}=$protID;
		$proteinLength{$protID}=$protLen;
	}

	my $sthTop=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=? AND VISIBILITY=2");
	my $sthBestVis=$dbh->prepare("SELECT MAX(VISIBILITY) FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=?");

	foreach my $identifier (keys %matchList) {
		my $protID=$identifiers{$identifier};
		$sthTop->execute($protID);
		($numProtTop{$protID})=$sthTop->fetchrow_array;
		if ($protVisibility) {
			$sthBestVis->execute($protID);
			($bestProtVis{$identifier})=$sthBestVis->fetchrow_array;
		}
		#>Peptide density
		$numProtPeptides{$protID}=scalar (keys %{$matchList{$identifier}});
		$proteinPepDens{$protID}=$proteinLength{$protID}/$numProtPeptides{$protID};
	}
	$sthTop->finish;
	$sthBestVis->finish;

	###> Fill former match-group information so as to only add virtual proteins!
	$sthGetAPInfo->execute($anaID);
	my $maxTopM=0;
	while (my ($identifier,$matchNumber,$vis) = $sthGetAPInfo->fetchrow_array) {
		$maxTopM=$matchNumber if $matchNumber>$maxTopM;
		$matchGroup{$identifier}=$matchNumber;
		$visibility{$identifier}=$vis;
	}

	###> Computing the specificity of each peptide for a given protein
	foreach my $peptideID (keys %pepProteins) {
		my $specificity=sprintf "%.2f",100/(scalar keys %{$pepProteins{$peptideID}});
		$specificity=~s/\.00//;
		foreach my $identifier (keys %{$pepProteins{$peptideID}}) {
			$pepSpecificity{$identifier}=$specificity if (!$pepSpecificity{$identifier} || $pepSpecificity{$identifier}<$specificity);
		}
	}

	foreach my $identifier (keys (%{$virtualProteins{$anaID}})) {# Only virtual-proteins have to be taking care of.

		###> Pep Coverage
		$sthGetPepInfo->execute($anaID,$identifiers{$identifier});
		#my $numMatch=0;
		#my $previousGP='';
		while (my ($pepBeg,$pepEnd,$pSeq,$pepID)=$sthGetPepInfo->fetchrow_array) {
			my $vmod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$anaID,$pSeq,\%anaLabelModifs);
			$vmod='' unless $vmod;
			$boundaryStatus{$identifiers{$identifier}}{$pepBeg}++;
			$boundaryStatus{$identifiers{$identifier}}{$pepEnd}--;
			$maxCovLength{$identifiers{$identifier}}=$pepEnd if (!$maxCovLength{$identifiers{$identifier}} || $maxCovLength{$identifiers{$identifier}} < $pepEnd);
			#$numMatch++ if ($previousGP ne "$pSeq:$vmod:$pepBeg");
			#$previousGP="$pSeq:$vmod:$pepBeg";
		}

		my $pepCoverage=0;
		if ($proteinLength{$identifiers{$identifier}} || $proteinLength{$identifiers{$identifier}}>0 ){
			my $coverage=0;
			my $hasPeptide=0;
			my $boundaryNter=0;
			foreach my $boundary (sort{$a<=>$b} keys %{$boundaryStatus{$identifiers{$identifier}}}) {
				if ($hasPeptide==0) { # start of peptide region (cannot become <0!)
					$boundaryNter=$boundary;
				}
				$hasPeptide+=$boundaryStatus{$identifiers{$identifier}}{$boundary};
				if ($hasPeptide==0) { # end of peptide region (should be true for last boundary too)
					$coverage+=($boundary-$boundaryNter)+1;
				}
			}

			#$pepCoverage=sprintf "%.1f",(100*$coverage)/$proteinLength{$identifiers{$identifier}};
			#$pepCoverage=~s/\.0//;
			$pepCoverage=($maxCovLength{$identifiers{$identifier}} <= $proteinLength{$identifiers{$identifier}})? sprintf "%.1f",(100*$coverage)/$proteinLength{$identifiers{$identifier}} : sprintf "%.1f",(-100*$coverage)/$maxCovLength{$identifiers{$identifier}}; # -: flag for protLength problem
			$pepCoverage*=1; # 25.0 -> 25
		}
		###> Pep Specificity
		#$sthInsAP->execute($anaID,$identifiers{$identifier},$dbRank{$identifier},scalar keys %{$matchList{$identifier}},$numMatch,$pepCoverage,$pepSpecificity{$identifier},$visibility{$identifier});
		#$sthInsAP->execute($anaID,$identifiers{$identifier},$dbRank{$identifier},scalar keys %{$matchList{$identifier}},$numMatch,$pepCoverage,$pepSpecificity{$identifier});
		$sthInsAP->execute($anaID,$identifiers{$identifier},$dbRank{$identifier},$pepCoverage,$pepSpecificity{$identifier});
	}

	###> Create Match Group.
    my @sortedIdentifiers=sort{$numProtPeptides{$identifiers{$b}}<=>$numProtPeptides{$identifiers{$a}} || $numProtTop{$identifiers{$b}}<=>$numProtTop{$identifiers{$a}} ||  $proteinScore{$anaID}{$identifiers{$b}}<=>$proteinScore{$anaID}{$identifiers{$a}} || &deltaLength($proteinLength{$identifiers{$a}},$proteinLength{$identifiers{$b}},$proteinPepDens{$identifiers{$a}},$proteinPepDens{$identifiers{$b}}) || $proteinLength{$identifiers{$a}}<=>$proteinLength{$identifiers{$b}} || $identifiers{$a}<=>$identifiers{$b}} keys %matchList;
	&promsMod::createMatchGroups(\%matchList,\%matchGroup,\%visibility,\%bestProtVis,$protVisibility,\@sortedIdentifiers,0);
	foreach my $identifier (keys %matchList){
		$sthUpMG->execute($visibility{$identifier},$matchGroup{$identifier},$anaID,$identifiers{$identifier});
	}

}
$sthGetPepInfo->finish;
$sthgetMatchInfo->finish;
$sthGetAPInfo->finish;
$sthInsAP->finish;
$sthUpMG->finish;
$dbh->commit;

############################
####>Deleting all files<####
############################
system "rm $quantifDir/*";


## Former code before modification on 05/03/2013 -> Update MatchGroup
##my $sthgetMatchInfo=$dbh->prepare("SELECT PEPTIDE.ID_PEPTIDE,PEP_SEQ,VAR_MOD,PROTEIN.ID_PROTEIN,IDENTIFIER,PROT_LENGTH FROM PEPTIDE,PEPTIDE_PROTEIN_ATTRIB,PROTEIN WHERE PEPTIDE.ID_PEPTIDE=PEPTIDE_PROTEIN_ATTRIB.ID_PEPTIDE AND PEPTIDE_PROTEIN_ATTRIB.ID_PROTEIN=PROTEIN.ID_PROTEIN AND PEPTIDE.ID_ANALYSIS=?");
##my $sthUpAP=$dbh->prepare("UPDATE ANALYSIS_PROTEIN SET VISIBILITY=?,MATCH_GROUP=?,NUM_PEP=?,NUM_MATCH=?,PEP_SPECIFICITY=?,PEP_COVERAGE=? WHERE ID_PROTEIN=? AND ID_ANALYSIS=?");
###my $sthGetPepInfo=$dbh->prepare("SELECT PEP_BEG,PEP_END FROM PEPTIDE,PEPTIDE_PROTEIN_ATTRIB WHERE PEPTIDE_PROTEIN_ATTRIB.ID_PEPTIDE=PEPTIDE.ID_PEPTIDE AND VALID_STATUS>=1 AND PEPTIDE_PROTEIN_ATTRIB.ID_ANALYSIS=? AND ID_PROTEIN=? ORDER BY PEP_BEG ASC, PEP_END ASC");
##my $sthGetPepInfo=$dbh->prepare("SELECT PEP_BEG,PEP_END FROM PEPTIDE,PEPTIDE_PROTEIN_ATTRIB WHERE PEPTIDE_PROTEIN_ATTRIB.ID_PEPTIDE=PEPTIDE.ID_PEPTIDE AND PEPTIDE_PROTEIN_ATTRIB.ID_ANALYSIS=? AND ID_PROTEIN=? ORDER BY PEP_BEG ASC, PEP_END ASC"); # Get PEPTIDE and GHOST-PEPTIDE to compute a well coverage on the protein
##
##my ($protVisibility)=$dbh->selectrow_array("SELECT PROT_VISIBILITY FROM PROJECT WHERE ID_PROJECT=$projectID");
##$protVisibility=0 unless $protVisibility;
##
##foreach my $anaID (keys %{$params{'MZXML'}}) {
##	my (%matchList,%matchGroup,%visibility,%identifiers,%proteinLength,%boundaryStatus,%pepSpecificity,%pepProteins,%maxCovLength,%numProtTop,%bestProtVis,%numProtPeptides,%proteinPepDens);
##	$sthgetMatchInfo->execute($anaID);
##
##	while(my ($pepID,$seq,$vmod,$protID,$identifier,$protLen)=$sthgetMatchInfo->fetchrow_array) {
##		$vmod='' unless $vmod;
##		$matchList{$identifier}{"$seq$vmod"}=1;
##		$pepProteins{$pepID}{$identifier}=1;
##		$identifiers{$identifier}=$protID;
##		$proteinLength{$protID}=$protLen;
##	}
##
##	my $sthTop=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=? AND VISIBILITY=2");
##	my $sthBestVis=$dbh->prepare("SELECT MAX(VISIBILITY) FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=?");
##	foreach my $identifier (keys %matchList) {
##		my $protID=$identifiers{$identifier};
##		$sthTop->execute($protID);
##		($numProtTop{$protID})=$sthTop->fetchrow_array;
##		if ($protVisibility) {
##			$sthBestVis->execute($protID);
##			$bestProtVis{$identifier}=$sthBestVis->fetchrow_array;
##		}
##		#>Peptide density
##		$numProtPeptides{$protID}=scalar (keys %{$matchList{$identifier}});
##		$proteinPepDens{$protID}=$proteinLength{$protID}/$numProtPeptides{$protID};
##	}
##	$sthTop->finish;
##	$sthBestVis->finish;
##
##	my @sortedIdentifiers=sort{$numProtPeptides{$identifiers{$b}}<=>$numProtPeptides{$identifiers{$a}} || $numProtTop{$identifiers{$b}}<=>$numProtTop{$identifiers{$a}} || $proteinScore{$identifiers{$b}}<=>$proteinScore{$identifiers{$a}} || &deltaLength($proteinLength{$identifiers{$a}},$proteinLength{$identifiers{$b}},$proteinPepDens{$identifiers{$a}},$proteinPepDens{$identifiers{$b}}) || $proteinLength{$identifiers{$a}}<=>$proteinLength{$identifiers{$b}} || $identifiers{$a}<=>$identifiers{$b}} keys %matchList;
##	&createMatchGroups(\%matchList,\%matchGroup,\%visibility,\%bestProtVis,$protVisibility,@sortedIdentifiers);
##
##	###> Computing the specificity of each peptide for a given protein
##	foreach my $peptideID (keys %pepProteins) {
##		my $specificity=sprintf "%.2f",100/(scalar keys %{$pepProteins{$peptideID}});
##		$specificity=~s/\.00//;
##		foreach my $identifier (keys %{$pepProteins{$peptideID}}) {
##			$pepSpecificity{$identifier}=$specificity if (!$pepSpecificity{$identifier} || $pepSpecificity{$identifier}<$specificity);
##		}
##	}
##
##	foreach my $identifier (keys (%matchList)) {
##
##		###> Pep Coverage
##		$sthGetPepInfo->execute($anaID,$identifiers{$identifier});
##		my $numMatch=0;
##		while (my ($pepBeg,$pepEnd)=$sthGetPepInfo->fetchrow_array) {
##			$pepBeg=abs($pepBeg);# Ghost-peptides could be with negatives pepBeg and pepEnd (a way to disinguish between regular one in PEPTIDE_PROTEIN_ATTRIB table)
##			$pepEnd=abs($pepEnd);
##			$boundaryStatus{$identifiers{$identifier}}{$pepBeg}++;
##			$boundaryStatus{$identifiers{$identifier}}{$pepEnd}--;
##			$maxCovLength{$identifiers{$identifier}}=$pepEnd if (!$maxCovLength{$identifiers{$identifier}} || $maxCovLength{$identifiers{$identifier}} < $pepEnd);
##			$numMatch++;
##		}
##
##		my $pepCoverage=0;
##		if ($proteinLength{$identifiers{$identifier}} || $proteinLength{$identifiers{$identifier}}>0 ){
##			my $coverage=0;
##			my $hasPeptide=0;
##			my $boundaryNter=0;
##			foreach my $boundary (sort{$a<=>$b} keys %{$boundaryStatus{$identifiers{$identifier}}}) {
##				if ($hasPeptide==0) { # start of peptide region (cannot become <0!)
##					$boundaryNter=$boundary;
##				}
##				$hasPeptide+=$boundaryStatus{$identifiers{$identifier}}{$boundary};
##				if ($hasPeptide==0) { # end of peptide region (should be true for last boundary too)
##					$coverage+=($boundary-$boundaryNter)+1;
##				}
##			}
##
##			#$pepCoverage=sprintf "%.1f",(100*$coverage)/$proteinLength{$identifiers{$identifier}};
##			#$pepCoverage=~s/\.0//;
##			$pepCoverage=($maxCovLength{$identifiers{$identifier}} <= $proteinLength{$identifiers{$identifier}})? sprintf "%.1f",(100*$coverage)/$proteinLength{$identifiers{$identifier}} : sprintf "%.1f",(-100*$coverage)/$maxCovLength{$identifiers{$identifier}}; # -: flag for protLength problem
##			$pepCoverage*=1; # 25.0 -> 25
##		}
##		###> Pep Specificity
##		$sthUpAP->execute($visibility{$identifier},$matchGroup{$identifier},scalar keys %{$matchList{$identifier}},$numMatch,$pepSpecificity{$identifier},$pepCoverage,$identifiers{$identifier},$anaID);
##	}
##}
##$dbh->commit;
##$sthgetMatchInfo->finish;
##$sthUpAP->finish;
#
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
$dbh->do("UPDATE QUANTIFICATION SET STATUS=1 WHERE ID_QUANTIFICATION=$quantifID");
$dbh->commit;
$dbh->disconnect;

open(FILESTAT,">>$fileStat");
print FILESTAT "Ended ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
close FILESTAT;

sleep 2;
system "rm $fileStat";


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
