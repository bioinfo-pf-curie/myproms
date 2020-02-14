#!/usr/local/bin/perl -w

################################################################################
# send2Biologist.cgi       3.3.4                                               #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Reports and/or terminates a validation (Analysis)                            #
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
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use POSIX qw(strftime); # to get the time
use strict;
use XML::Simple; # used by promsMod::extractData & quantif from MSF
use File::Path qw(rmtree); # remove_tree
use File::Copy qw(copy move); # also needed in &promsMod::removeFromMultiAna
use File::Copy::Recursive qw(dirmove);
use phosphoRS;
use HTML::Entities;

#############################################
####>Configuration & Fetching parameters<####
#############################################
my $call=($ARGV[0])? 'cmd' : (param('call'))? param('call') : 'std'; # called from user (std) or ajax from editBranch to end project (ajax) or command line from selectProject.cgi (cmd)
my %promsPath;
my (@analysisList,$userID,$endValidation,$jobID); # $jobID only defined for call=cmd
if ($call eq 'cmd') { # called from command line
	%promsPath=&promsConfig::getServerInfo('no_user');
	$jobID=$ARGV[0];
	@analysisList=split(':',$ARGV[1]);
	$endValidation=1;
	$userID='myproms'; # not used if only ending already reported validation(s)
}
else { # called from http protocol
	%promsPath=&promsConfig::getServerInfo; # sets $ENV{'REMOTE_USER'}!!!
	if (param('ID')) {$analysisList[0]=param('ID');} # obsolete since v3.0
	elsif (param('anaList')) {
		@analysisList=($call eq 'ajax')? split(':',param('anaList')) : param('anaList');
	}
	$endValidation=param('endVal'); # 0=send only, 1=end only, 2=send + end
	$userID=$ENV{'REMOTE_USER'};
}
#my $msType=param('MSTYPE');  #fp
####------>Test
#print header; warningsToBrowser(1);
#print "CALL=$call<BR>\n";
#print "USER=$userID<BR>\n";
#print "NUM_ANA=",scalar @analysisList," (@analysisList)<BR>\n";
#print "END_VAL=$endValidation<BR>\n";
#exit;
####<------End

##########################
####>Connecting to DB<####
##########################
my $dbh=($call eq 'cmd')? &promsConfig::dbConnect('no_user') : &promsConfig::dbConnect;

#############################
####>Fetching project ID<####
#############################
my $numAna=scalar @analysisList;
my @itemInfo=&promsMod::getItemInfo($dbh,'analysis',$analysisList[0]);
my $projectID=$itemInfo[0]{'ID'};
my ($protVisibility,$projectIdentMapID)=$dbh->selectrow_array("SELECT PROT_VISIBILITY,ID_IDENTIFIER FROM PROJECT WHERE ID_PROJECT=$projectID");
$protVisibility=0 unless $protVisibility;
$projectIdentMapID=0 unless $projectIdentMapID;
my ($giIdentID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GI'");

#######################
####>Starting HTML<####
#######################
if ($call eq 'cmd') {
	mkdir $promsPath{'logs'} unless -e $promsPath{'logs'};
	my $year=strftime("%Y",localtime);
	my $logFile="$promsPath{logs}/autoEndValidation_$year.log";
	open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
	open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
	open STDERR, ">>$logFile";
	open STDOUT, ">>$logFile";
	$| = 1;
	print "\n>JOB START $jobID, $numAna analyses for project #$projectID ($ARGV[1]):\n";
}
else {
	my $titleString=($endValidation==1)? 'Ending Validation...' : 'Reporting Validation...';
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	if ($call eq 'ajax') {
		print "<HTML><BODY>\n";
	}
	else {
		print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<IFRAME name="systemFrame" style="display:none"></IFRAME>
<CENTER>
<IMG src="$promsPath{images}/engrenage.gif">
</CENTER>
<BR>
<FONT class="title">$titleString</FONT>
<BR>
|;
	}
	#foreach my $i (1..200) {print "<!--BUFFER TEXT XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX-->\n";}
}

#####################################################################
####>Fetching list of existing Reference Retention Time peptides<####
#####################################################################
my %referenceRTdata;
unless ($endValidation==1) { # unless directly ending validation (No new data sent to biologists)
	my $sthAllRT=$dbh->prepare("SELECT ID_REFERENCE_RT,NAME,PROT_IDENTIFIER,DATA FROM REFERENCE_RT");
	$sthAllRT->execute;
	while (my ($refID,$refRTname,$protIdentifier,$pepData)=$sthAllRT->fetchrow_array) {
		$protIdentifier=~s/\s.*//; # remove any space & following characters
		$referenceRTdata{$refID}{'NAME'}=$refRTname;
		$referenceRTdata{$refID}{'IDENTIFIER'}=$protIdentifier;
		$protIdentifier=~s/^(gi|sp|tr|ipi)\|//i;
		my $identQueryStrg='';
		foreach my $subIdent (split(/\|/,$protIdentifier)) {
			$identQueryStrg.=" OR " if $identQueryStrg;
			$identQueryStrg.="IDENTIFIER REGEXP '(^|[[.|.]])$subIdent(-[[:digit:]]+)*(\$|[[.|.]])'"; # Now match single-peptide refRT protein identifier member (<identifier>-n)
		}
#print "<BR>'$identQueryStrg'<BR>\n";
		$referenceRTdata{$refID}{'QUERY'}=$dbh->prepare("SELECT ID_PROT_VALID,IDENTIFIER FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=? AND ($identQueryStrg) AND IDENTIFIER NOT LIKE 'DECOY_%'");
#print "<BR>$refRT_ID,$refRTname,$protIdentifier: SELECT ID_PROT_VALID,IDENTIFIER FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=? AND ($identQueryStrg) AND IDENTIFIER NOT LIKE 'DECOY_%'";
		my $xmlRT = new XML::Simple();
		my $xmlData = $xmlRT->XMLin($pepData);
		foreach my $pepInfo (@{$xmlData->{PEPTIDE}}) {
			next if $pepInfo->{excluded}; # peptide should not be used
			@{$referenceRTdata{$refID}{'PEPTIDES'}{$pepInfo->{sequence}}}=($pepInfo->{pepId},$pepInfo->{iRT},($pepInfo->{monoMass}-1.008)*$pepInfo->{charge}); # massExp
		}
	}
	$sthAllRT->finish;
}

######################################################################
####>Fetching list of existing proteins for corresponding project<####
######################################################################
####>Fetching all project's proteins (they will be substracted from the new list => no duplicates)
my (%projectProt,%proteinLength,%incompleteProtein,%noDescriptionProtein);
my @newAnaMapping; # if not empty => launches identifier mapping if set
unless ($endValidation==1) { # unless directly ending validation (No new data sent to biologists)
	print "<BR><FONT class='title3'>- Scanning project's proteins...";
	my $sthP=$dbh->prepare("SELECT ID_PROTEIN,IDENTIFIER,PROT_LENGTH,PROT_DES FROM PROTEIN WHERE ID_PROJECT=$projectID");
	#my $sthDesc=$dbh->prepare("SELECT PROT_DES FROM PROTEIN WHERE ID_PROTEIN=?");
	$sthP->execute || die $sthP->errstr;
	my $refProtData=$sthP->fetchall_arrayref; # reference to an array
	$sthP->finish;
	foreach my $refProt (@{$refProtData}) {
		next unless $refProt->[1]; # identifier is null (previous ID protection not removed)
		$projectProt{$refProt->[1]}=$refProt->[0]; # {identifier}=id
		$proteinLength{$refProt->[0]}=$refProt->[2]; #{id}=length
		if (!$refProt->[2]) { # no length
			$incompleteProtein{$refProt->[1]}=$refProt->[0]; # {identifier}=id
			if (!$refProt->[3] || $refProt->[3]=~/no\sdescription/i) { # no description
				$noDescriptionProtein{$refProt->[1]}=$refProt->[0]; # {identifier}=id
			}
		}
	}
	#foreach my $identifier (keys %incompleteProtein) {
	#	#print "$identifier,$incompleteProtein{$identifier}\n";
	#	$sthDesc->execute($incompleteProtein{$identifier});
	#	my $localDesc=$sthDesc->fetchrow_array;
	#	if (!$localDesc || $localDesc=~/no\sdescription/i) {
	#		$noDescriptionProtein{$identifier}=$incompleteProtein{$identifier};## {identifier}=id
	#	}
	#}
	#$sthDesc->finish;
	print " Done.</FONT><BR>\n";
}

##################################
####>Looping through Analyses<####
##################################
#my @querySummaryAttributes=('qmass','qexp','qmatch','qplughole'); # list of query attributes to be store in summary section of Pxxxxx.(dat/pdm) file
#my @queryPeptidesAttributes=('','_terms','_primary_nl','_subst'); # list of peptide attributes to be store in summary section of Pxxxxx.(dat/pdm) file
my $sthAnaInfo=$dbh->prepare("SELECT NAME,MS_TYPE,VALID_STATUS,DATA_FILE,FILE_FORMAT,INSTRUMENT,LABELING,MAX_RANK,DECOY,MIN_SCORE,ID_SAMPLE FROM ANALYSIS WHERE ID_ANALYSIS=?");
my $sthSpecies=$dbh->prepare("SELECT SCIENTIFIC_NAME FROM SPECIES SP,EXPERIMENT E,SAMPLE S WHERE S.ID_EXPERIMENT=E.ID_EXPERIMENT AND E.ID_SPECIES=SP.ID_SPECIES AND S.ID_SAMPLE=?");
my $sthLabMod=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_MODIFICATION AM,MODIFICATION M WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND PSI_MS_NAME LIKE 'Label:%' AND AM.ID_ANALYSIS=?");

my $sthPepRT=$dbh->prepare("SELECT QV.QUERY_NUM,EXT_SPECTRUMID,MASS_DATA,CHARGE,ELUTION_TIME,INFO_PEP1 FROM RANK_PROTEIN_MATCH RPM,QUERY_VALIDATION QV WHERE RPM.ID_ANALYSIS=? AND QV.ID_ANALYSIS=? AND RPM.QUERY_NUM=QV.QUERY_NUM AND IDENTIFIER=? ORDER BY MAX_SCORE DESC"); # AND VALID_STATUS>=1
my $sthInsART=$dbh->prepare("INSERT INTO ANALYSIS_REFRT (ID_REFERENCE_RT,ID_ANALYSIS,NUM_PEP,DATA) VALUES (?,?,?,?)");
my $sthInsQRT=$dbh->prepare("INSERT INTO QUANTIF_REFRT (ID_REFERENCE_RT,ID_QUANTIFICATION,SOURCE_RANK,NUM_PEP,SLOPE,Y_INTERCEPT,CORRELATION,DATA) VALUES (?,?,?,?,?,?,?,?)");

my $countAna=0;
my $prevValidation=0;
my $msfExtraction=0;
my %msfReportedList; # used by &endValidation for msf files only
foreach my $analysisID (@analysisList) {
	$countAna++;
	$sthAnaInfo->execute($analysisID);
	my ($anaName,$msType,$validStatus,$dataFileName,$fileFormat,$instrument,$labelMethod,$maxRank,$decoy,$minScorePep,$sampleID)=$sthAnaInfo->fetchrow_array;
	next if $validStatus==2; # in case call=cmd and multiple jobs called for same anaList
	$sthSpecies->execute($sampleID); # can be undef
	my ($preferredSpecies)=$sthSpecies->fetchrow_array;
	my %protSpeciesClass; # {idProt}=: 1 if preferred species is defined & matched, 0 otherwise (needed for match group decision)
	if (!$labelMethod || $labelMethod !~ /SILAC|iTRAQ/i) {
		$sthLabMod->execute($analysisID);
		my ($hasSilacMods)=$sthLabMod->fetchrow_array;
		$labelMethod=($hasSilacMods)? 'SILAC' : ($labelMethod)? $labelMethod : ''; # overwrites true method name
	}
	#$maxRank=10 if $msType eq 'PMF';
	$maxRank=&promsConfig::getMaxRank if (!$maxRank || $msType ne 'MIS'); # default to max
	##>Files
	mkdir "$promsPath{peptide}/proj_$projectID" unless -e "$promsPath{peptide}/proj_$projectID";
	my ($dataFile,$msfFile,$msfFileID); #,$pepFile,$oldPepFile; #,%sectionField
	if ($fileFormat eq 'MASCOT.DAT' || $fileFormat=~/\.PDM/) { # **$msType ne 'PMF' && ()** PMF also transferred since v2.7
		$dataFile="$promsPath{valid}/ana_$analysisID/$dataFileName";
		if ($fileFormat=~/\.PDM/) {
			($msfFile=$dataFileName)=~s/_\d+\.*\d*\.pdm/\.msf/;
			if ($dataFileName =~ /_\d+\.(\d+)\.pdm/) {
				$msfFileID=$1;
			}
		}
	}
	elsif ($fileFormat eq 'PHENYX.XML' || $fileFormat eq 'MASCOT.XML') {
		(my $realFileName=$dataFileName)=~s/\.xml/\.pgf/;
		$dataFile="$promsPath{valid}/ana_$analysisID/$realFileName";
	}
	elsif ($fileFormat eq 'PARAGON.XML') {
		$dataFile="$promsPath{valid}/ana_$analysisID/$dataFileName";
	}
	print "<BR><FONT class='title3'>- Processing Analysis $anaName ($countAna/$numAna):</FONT><BR>\n" if $call ne 'cmd';

	####################################
	####>Directly ending validation<#### No new data sent to biologists
	####################################
	if ($endValidation==1) {

		########################################
		####>Updating Analysis valid_status<####
		########################################
		#$validStatus=($endValidation)? 2 : 1;
		if ($validStatus==1) { # do not change original validator
			$dbh->do("UPDATE ANALYSIS SET VALID_STATUS=2 WHERE ID_ANALYSIS=$analysisID") || die $dbh->errstr;
		}
		else { # -1 or 0
			$dbh->do("UPDATE ANALYSIS SET VALID_STATUS=2,VALID_USER='$userID',FIRST_VALID_DATE=NOW(),VALID_DATE=NOW() WHERE ID_ANALYSIS=$analysisID") || die $dbh->errstr;
		}

		###########################
		####>Ending validation<####
		###########################
		&endValidation($analysisID,$fileFormat,$dataFileName); #,$dataFile,$pepFile,$oldPepFile,$decoyFile,$decoyPepFile,$oldDecoyPepFile

		$dbh->commit;

		#print " Done.</FONT><BR>\n";
		if ($call eq 'cmd') {
			print '.';
			print "\n" unless $countAna % 100;
		}
		next; # next Analysis
	}

	###>Set Analysis as partially validated (better if error after partial data transfer)
	$dbh->do("UPDATE ANALYSIS SET VALID_STATUS=1,VERIFIED_MG=0,VALID_USER='$userID',VALID_DATE=NOW(),UPDATE_USER='$userID',UPDATE_DATE=NOW() WHERE ID_ANALYSIS=$analysisID");

	##################################################
	####>Processing Quantification settings (MSF)<#### Done here because match between search varMod & SILAC labeling varMod is necessary early on
	##################################################
	my ($dbhLite,$pdVersion,$pdVersionComp,$searchNodeNum,$quantifNodeNum,$quantifName,$quantifAnnot,%labelingInfo,%reporterIonList,%qChannel2TargetPos,$noLabelChannel,%peptideExtraData,%termModifs,$labelResStrg,%mergedFiles,$guid);
	my $numMergedFiles=0;
	if ($msfFile) {
		my $dbFile="$promsPath{valid}/multi_ana/proj_$projectID/$msfFile";
		$dbhLite=DBI->connect("dbi:SQLite:$dbFile", "", "", {PrintError => 1,RaiseError => 1});
		($pdVersionComp)=$dbhLite->selectrow_array("SELECT SoftwareVersion FROM SchemaInfo ORDER BY rowid ASC LIMIT 1");
		$pdVersion=$pdVersionComp;
		$pdVersion=~s/^(\d\.\d+).*/$1/; # x.x.x.xx -> x.x (numerical value)

		##>Finding Search ProcessingNodeNumber & varMods
		if ($pdVersion >= 2.0) {
			my ($wfXML) = $dbhLite->selectrow_array("SELECT WorkflowXML FROM Workflows");
			my $xml = new XML::Simple();
			my $xmlData = $xml->XMLin($wfXML);
			my ($parenProcessingNumber,$processingNodeName,$processingNodeNumber);
			my $prsNodeID=0;
			my $softUsed=($fileFormat eq 'MASCOT.PDM')?'MASCOT':'SEQUEST';
			foreach my $processingNodeParameters (sort{$a->{ParentProcessingNodeNumber} cmp $b->{ParentProcessingNodeNumber}} @{$xmlData->{WorkflowTree}->{WorkflowNode}}) {
				next unless uc($processingNodeParameters->{FriendlyName}) =~ $softUsed;
				$searchNodeNum=$processingNodeParameters->{ProcessingNodeNumber};
				$guid=$processingNodeParameters->{Guid}; # use to get scoreID
				last;
			}
		}
		else {
			if ($fileFormat eq 'MASCOT.PDM') {
				($searchNodeNum)=$dbhLite->selectrow_array("SELECT ProcessingNodeNumber FROM ProcessingNodes WHERE NodeName='Mascot'");
			}
			else { # SEQUEST
				if ($pdVersion >= 1.4) {
					# For PD1.4 -> new name in the database in ProcessingNodes: NodeName='IseNode' & FriendlyName='Sequest HT' for SEQUEST searches
					# 2 different versions: SEQUEST (old version) or Sequest HT (high throughput for multi processing use)
					($searchNodeNum)=$dbhLite->selectrow_array("SELECT ProcessingNodeNumber FROM ProcessingNodes WHERE UPPER(FriendlyName) LIKE 'SEQUEST%'");
				}
				else { # 1.2/3
					($searchNodeNum)=$dbhLite->selectrow_array("SELECT ProcessingNodeNumber FROM ProcessingNodes WHERE NodeName='SequestNode'");
				}
			}
		}

		##>Fetching list of merged files
		my $fileTable=($pdVersion >= 2)? 'WorkflowInputFiles' : 'FileInfos';
		if ($msfFileID) { # sub import
			$numMergedFiles=1;
			my ($file)=$dbhLite->selectrow_array("SELECT FileName FROM $fileTable WHERE FileID=$msfFileID");
			my $fileName=(split(/[\/\\]/,$file))[-1];
			@{$mergedFiles{$msfFileID}}=($fileName,$numMergedFiles);
		}
		else { # must not be sub import
			my $sthMF=$dbhLite->prepare("SELECT FileID,FileName FROM $fileTable ORDER BY FileID");
			$sthMF->execute;
			while (my ($fileID,$file)=$sthMF->fetchrow_array) {
				$numMergedFiles++;
				my $fileName=(split(/[\/\\]/,$file))[-1];
				@{$mergedFiles{$fileID}}=($fileName,$numMergedFiles);
			}
			$sthMF->finish;

			if ($numMergedFiles > 1) {
				open(MERGE,">$promsPath{valid}/ana_$analysisID/mergedFiles.txt");
				foreach my $fileID (sort{$a<=>$b} keys %mergedFiles) { # sort by fileID = rank
					print MERGE "$mergedFiles{$fileID}[1]\t$mergedFiles{$fileID}[0]\t$fileID\n"; # rank,name,id
				}
				close MERGE;
			}
		}
	}

	if ($msfFile && $labelMethod=~/SILAC|TMT|iTRAQ/i) {
		#my $sthVarMods = $dbhLite->prepare("SELECT ValueDisplayString FROM ProcessingNodeParameters WHERE ProcessingNodeNumber=$searchNodeNum AND ParameterName LIKE ? ESCAPE ? AND ParameterName NOT LIKE ?");
		#my $sthDelMass = $dbhLite->prepare("SELECT DeltaMass FROM AminoAcidModifications WHERE ModificationName=? LIMIT 0,1"); # multiple rows are possible
		#$sthVarMods->execute('DynMod\_%','\\','Search%');
		my (%varModDeltaMass,%targetedResidue);
		#while (my ($varMod)=$sthVarMods->fetchrow_array) {
		#	(my $subVarMod=$varMod)=~s/ \([^\(]+\)\Z//; # Oxidation (M) -> Oxidation
		#	$sthDelMass->execute($subVarMod);
		#	($varModDeltaMass{$varMod})=$sthDelMass->fetchrow_array;
		#}
		#$sthVarMods->finish;
		#$sthDelMass->finish;
		my ($paramValue,$processMatched,$processCount)=('',0,0);
		#my %reporterIonList;
		if ($pdVersion >= 2.0 ) { # TODO : check TMT method for lower versions of PD
			my ($wfXML) = $dbhLite->selectrow_array("SELECT WorkflowXML FROM Workflows");
			my $xml = new XML::Simple();
			my $xmlData = $xml->XMLin($wfXML);
			###>Finding VMODs
			foreach my $processingNodeParameters (@{$xmlData->{WorkflowTree}->{WorkflowNode}}) {
				next unless $processingNodeParameters->{ProcessingNodeNumber} == $searchNodeNum;
				foreach my $processingNodeParameter (@{$processingNodeParameters->{ProcessingNodeParameters}->{ProcessingNodeParameter}}){
					if(($processingNodeParameter->{IntendedPurpose} =~ /Dynamic(Terminal)?Modification/) && $processingNodeParameter->{IsValueSet} eq "True"){
						my $dynmodXML=new XML::Simple();
						my $dynmodDesc=$dynmodXML->XMLin($processingNodeParameter->{content}); #<Modification Version="2" AminoAcids="C" Name="Carbamidomethyl" Abbreviation="Carbamidomethyl" ID="8" UnimodAccession="4" DeltaMass="57.02146" DeltaAverageMass="57.05130" IsSubstitution="False" LeavingGroup="" Substitution="H(3) C(2) N O" PositionType="Any" />
						my $res=($dynmodDesc->{AminoAcids})?$dynmodDesc->{AminoAcids} : $dynmodDesc->{Terminus};
						my $valueD="$dynmodDesc->{Name} ($res)";
						my $paramN=$processingNodeParameter->{Name};
						my $varModName=$valueD;
						$varModDeltaMass{$varModName}=$dynmodDesc->{DeltaMass};
						$targetedResidue{$varModName}=$dynmodDesc->{AminoAcids};
						$termModifs{$varModName}=($dynmodDesc->{AminoAcids}=~/N-term/i)? 'N' : ($dynmodDesc->{AminoAcids}=~/C-term/i)? 'C' : 0;
					}
				}
			}

			###>Finding Quantification NodeID
			foreach my $processingNodeParameters (@{$xmlData->{WorkflowTree}->{WorkflowNode}}) {
				### PrecursorIonsQuantifierNode -> SILAC
				### ReporterIonQuantifierNode -> TMT / ITRAQ
				next unless ($processingNodeParameters->{ProcessingNodeName} =~ /PrecursorIonsQuantifierNode|ReporterIonQuantifierNode|MinoraFeatureCreatorNode/);
				foreach my $processingNodeParameter (@{$processingNodeParameters->{ProcessingNodeParameters}->{ProcessingNodeParameter}}){
					if($processingNodeParameter->{IntendedPurpose} eq "QuantificationMethod"){
						$paramValue=$processingNodeParameter->{content};
						($quantifName)=($paramValue=~/^<ProcessingMethod name=\"(.*)\" version/);
						$quantifNodeNum=$processingNodeParameters->{ProcessingNodeNumber};
						if ($labelMethod =~ /TMT|iTRAQ/i) { # get the reporterIonList
								my $tmtData=$xml->XMLin($processingNodeParameter->{content});
								foreach my $repMass (keys %{$tmtData->{MethodPart}{MethodPart}}){
										#next unless $tmtData->{MethodPart}{MethodPart}{$repMass}{Parameter}{IsActive}{content}; # Commented on 2017/11/02 because of PD 2.0.0.802 version where IsActive do not exist
										my $channelID=$tmtData->{MethodPart}{MethodPart}{$repMass}{Parameter}{TagID}{content};
										my $monoMass=$tmtData->{MethodPart}{MethodPart}{$repMass}{Parameter}{MonoisotopicMZ}{content};
										next unless ($channelID && $monoMass);
										#my $tagName=$tmtData->{MethodPart}{MethodPart}{$tmtMass}{Parameter}{TagName}{content};
										@{$reporterIonList{$repMass}}=($channelID,$monoMass);
								}
						}

					}
				}
				if ($paramValue) {
					my $parentNodeStrg=$processingNodeParameters->{ParentProcessingNodeNumber};
					$processMatched=0;
					foreach my $parNodeNum (split(/;/,$parentNodeStrg)) {
						if ($parNodeNum==$searchNodeNum) {
							$processMatched=1;
							last;
						}
					}
				}
			}
		}
		else{
			my $sthVarMods = $dbhLite->prepare("SELECT ValueDisplayString FROM ProcessingNodeParameters WHERE ProcessingNodeNumber=$searchNodeNum AND IntendedPurpose IN (5,13)"); # !!! WARNING: 5,13 may not be stable nor complete !!!
			my $sthDelMass = $dbhLite->prepare("SELECT DeltaMass FROM AminoAcidModifications WHERE ModificationName=? LIMIT 0,1"); # multiple rows are possible
			$sthVarMods->execute;
			while (my ($varModData)=$sthVarMods->fetchrow_array) { # Carbamidomethyl / +57.021 Da (C)
				(my $varModName=$varModData)=~s/ \/ [\+\-].+//; # Carbamidomethyl / +57.021 Da (C) -> Carbamidomethyl
				my ($deltaMass)=($varModData=~/\/ ([\+\-\d\.]+)/); # Carbamidomethyl / +57.021 Da (C) -> +57.021 ONLY SEQUEST!!!
				unless ($deltaMass) { # MASCOT
					$varModName=~s/ \([^\(]+\)\Z//; # Oxidation (M) -> Oxidation
					$sthDelMass->execute($varModName);
					($deltaMass)=$sthDelMass->fetchrow_array;
				}
				$varModDeltaMass{$varModName}=$deltaMass;
				($targetedResidue{$varModName})=($varModData=~/\(([^\(]+)\)\Z/); # Carbamidomethyl / +57.021 Da (C) -> C
				$termModifs{$varModName}=($targetedResidue{$varModName}=~/N-term/i)? 'N' : ($targetedResidue{$varModName}=~/C-term/i)? 'C' : 0;
#print "VM DATA=$varModData NAME='$varModName' DELTA='$varModDeltaMass{$varModName}' RES='$targetedResidue{$varModName}'<BR>\n";
			}
			$sthVarMods->finish;
			$sthDelMass->finish;
#exit;
			###>Finding Quantification NodeID
			my $sthQN=$dbhLite->prepare("SELECT ProcessingNodes.ProcessingNodeNumber,ProcessingNodeParentNumber,ValueDisplayString,ParameterValue FROM ProcessingNodeParameters,ProcessingNodes WHERE ProcessingNodeParameters.ProcessingNodeID=ProcessingNodes.ProcessingNodeID AND ParameterName LIKE 'Quanti%ationMethod'"); # Quantification or Quantitation
			my $sthPN=$dbhLite->prepare("SELECT ProcessingNodeParentNumber FROM ProcessingNodes WHERE ProcessingNodeNumber=?");
			$sthQN->execute;
			#QNODE:while (($quantifNodeNum,my $parentNodeStrg,$quantifName,$paramValue)=$sthQN->fetchrow_array) {
			#	foreach my $parNodeNum (split(/;/,$parentNodeStrg)) {
			#		last QNODE if $parNodeNum==$searchNodeNum;
			#	}
			#}
			while (my ($qNodeNum,$parentNodeStrg,$qName,$pValue)=$sthQN->fetchrow_array) {
				$processCount++;
				$quantifNodeNum=$qNodeNum;
				$quantifName=$qName;
				$paramValue=$pValue;
				$processMatched=&checkSearchNode($sthPN,$searchNodeNum,$parentNodeStrg);
			}
			$sthQN->finish;
			$sthPN->finish;
		}

		if ($processMatched==0 && $processCount > 1) { # multiple quantifs that cannot be linked to the $processingNodeNumber used => skip quantif data
			$labelMethod=$quantifName=$paramValue='';
			print "[*Unable to match quantification process to search process: Skipping quantification data*]";
		}

		###>Processing parameter string
		if ($paramValue && $labelMethod=~/SILAC/i) {
			$noLabelChannel=0; # default in case no-unlabeled state
			$paramValue=~s/""/"/g;
			my $xml = new XML::Simple();
			my $xmlData = $xml->XMLin($paramValue);
			my ($labelXLM)=($pdVersion >= 2.0 ) ? $xmlData->{MethodPart}{MethodPart} : $xmlData->{MethodPart}{QuanChannels}{MethodPart}; # slight change for new PD 2.0 version
#		foreach my $labelName (keys %{$labelXLM}){
#			my $qChannelID=$labelXLM->{$labelName}{Parameter}{ChannelID}{content};
#			my $labelModifName='No label'; # not Labeled
#			my $modifRes='';
#			my $searchModifName='';
#			my $quotedSearchModif='';
#			if ($labelXLM->{$labelName}{MethodPart}{MethodPart}{MethodPart}{MethodPart}{SideChainModification}{Parameter}) {
#				my $labelParamStrg=$labelXLM->{$labelName}{MethodPart}{MethodPart}{MethodPart}{MethodPart}{SideChainModification}{Parameter}{content};
#				($labelModifName)=$labelParamStrg=~/ Name="([^"]+)/;
#				($modifRes)=$labelParamStrg=~/ AminoAcids="([^"]+)/;
#				my ($deltaMass)=$labelParamStrg=~/ DeltaMass="([^"]+)/;
#				##>Matching search varMod (could be the same label or a user-defined one! eg. 'Label:13C(8)' <=> 'Lys+8')
#				foreach my $varMod (keys %varModDeltaMass) {
#					if (abs($deltaMass-$varModDeltaMass{$varMod})<=0.2) { # match!
#						($searchModifName=$varMod)=~s/ \([A-Z]+\)\Z//; # remove targetted res
#						(my $trunckModif=$varMod)=~s/( \([A-Z]+)\)\Z/$1/; # 'Labelxxx (Res)' -> 'Labelxxx (Res' : To allow AA position-independent match
#						$quotedSearchModif=quotemeta($trunckModif);
#						last;
#					}
#				}
#			}
#			else {$noLabelChannel=$qChannelID;} # channelID for unlabeled peptides
#			@{$labelingInfo{$qChannelID}}=($labelName,$labelModifName,$modifRes,$searchModifName,$quotedSearchModif); # 'Heavy','Label:13C(6)','K','Lys+6'   quotemeta because of matches
#		}
#		$quantifAnnot='LABEL=SILAC';
#		foreach my $qChannelID (sort{$a<=>$b} keys %labelingInfo) {
#			$quantifAnnot.='::';
#			$quantifAnnot.=join(';',$qChannelID,@{$labelingInfo{$qChannelID}}[0..3]);
##print "$qChannelID => '@{$labelingInfo{$qChannelID}}'<BR>\n";
#		}


			foreach my $labelName (keys %{$labelXLM}){
#print ">$labelName:\n";
				my $qChannelID=$labelXLM->{$labelName}{Parameter}{ChannelID}{content};
				$labelingInfo{$qChannelID}{'NAME'}=$labelName;
#print Dumper($labelXLM->{$labelName}{MethodPart}{MethodPart}),"\n\n";
#foreach my $key (sort keys %{$labelXLM->{$labelName}{MethodPart}{MethodPart}}) {print "\t-$key\n";}
				if ($labelXLM->{$labelName}{MethodPart}{MethodPart}{MethodPart}) { # No label or single label
#print "*$labelXLM->{$labelName}{MethodPart}{MethodPart}{name}*\n";
					my $labelModifAlias=$labelXLM->{$labelName}{MethodPart}{MethodPart}{name};
#print "+$labelModifAlias:\n";
					if ($labelXLM->{$labelName}{MethodPart}{MethodPart}{MethodPart}{MethodPart}{SideChainModification}{Parameter}) {
						my $labelParamStrg=decode_entities(decode_entities($labelXLM->{$labelName}{MethodPart}{MethodPart}{MethodPart}{MethodPart}{SideChainModification}{Parameter}{content}));
#print "$labelParamStrg\n";
						my ($labelModifName)=$labelParamStrg=~/ Name="([^"]+)/;
						my ($modifRes)=$labelParamStrg=~/ AminoAcids="([^"]+)/;
						my ($deltaMass)=$labelParamStrg=~/ DeltaMass="([^"]+)/;
						my $searchModifName='';
						#my $quotedSearchModif='';
						##>Matching search varMod (could be the same label or a user-defined one! eg. 'Label:13C(8)' <=> 'Lys+8')
						foreach my $varModName (keys %targetedResidue) {
							#next if $varMod !~ /\($modifRes\)\Z/; # varMod does not target this residue
							next if $targetedResidue{$varModName} ne $modifRes; # varMod does not target this residue
#print "RES1: $varModName ($modifRes) DELTA=",abs($deltaMass-$varModDeltaMass{$varModName}),"<BR>\n";
							if (abs($deltaMass-$varModDeltaMass{$varModName})<=0.2) { # match!
								#($searchModifName=$varMod)=~s/ \([A-Z]+\)\Z//; # remove targeted res
								$searchModifName=$varModName;
								##(my $trunckModif=$varMod)=~s/( \([A-Z]+)\)\Z/$1/; # 'Labelxxx (Res)' -> 'Labelxxx (Res' : To allow AA position-independent match
								##$quotedSearchModif=quotemeta($trunckModif);
								last;
							}
						}
#print "OK1: $searchModifName<BR>\n";
						push @{$labelingInfo{$qChannelID}{'PARAMS'}},[$labelModifAlias,$labelModifName,$modifRes,$searchModifName]; #,$quotedSearchModif # 'Heavy','Label:13C(6)','K','Lys+6'   quotemeta because of matches
						$labelResStrg.=$modifRes;
					}
					else {
						push @{$labelingInfo{$qChannelID}{'PARAMS'}},[$labelModifAlias,'No label','','',''];
						$noLabelChannel=$qChannelID;
					} # channelID for unlabeled peptides
				}
				elsif (scalar keys %{$labelXLM->{$labelName}{MethodPart}{MethodPart}}) {
					foreach my $labelModifAlias (sort{lc($a) cmp lc($b)} keys %{$labelXLM->{$labelName}{MethodPart}{MethodPart}}) {
#print "+$labelModifAlias:\n";
#print Dumper($labelXLM->{$labelName}{MethodPart}{MethodPart}{$coLabelName}),"\n\n";
						my $labelParamStrg=decode_entities(decode_entities($labelXLM->{$labelName}{MethodPart}{MethodPart}{$labelModifAlias}{MethodPart}{MethodPart}{SideChainModification}{Parameter}{content}));
#print "$labelParamStrg\n";
						my ($labelModifName)=$labelParamStrg=~/ Name="([^"]+)/;
						my ($modifRes)=$labelParamStrg=~/ AminoAcids="([^"]+)/;
						my ($deltaMass)=$labelParamStrg=~/ DeltaMass="([^"]+)/;
						my $searchModifName='';
						#my $quotedSearchModif='';
						##>Matching search varMod (could be the same label or a user-defined one! eg. 'Label:13C(8) (K)' <=> 'Lys+8 (K)')
						foreach my $varModName (keys %targetedResidue) {
							#next if $varMod !~ /\($modifRes\)\Z/; # varMod does not target this residue
							next if $targetedResidue{$varModName} ne $modifRes; # varMod does not target this residue
#print "RES2: $varModName ($modifRes) DELTA=",abs($deltaMass-$varModDeltaMass{$varModName}),"<BR>\n";
							if (abs($deltaMass-$varModDeltaMass{$varModName})<=0.2) { # match!
								#($searchModifName=$varMod)=~s/ \([A-Z]+\)\Z//; # remove targetted res
								$searchModifName=$varModName;
								#(my $trunckModif=$varMod)=~s/( \([A-Z]+)\)\Z/$1/; # 'Labelxxx (Res)' -> 'Labelxxx (Res' : To allow AA position-independent match
								#$quotedSearchModif=quotemeta($trunckModif);
								last;
							}
						}
#print "OK2: $searchModifName<BR>\n";
						push @{$labelingInfo{$qChannelID}{'PARAMS'}},[$labelModifAlias,$labelModifName,$modifRes,$searchModifName]; #,$quotedSearchModif
						$labelResStrg.=$modifRes;
					}
				}
				else { # Assume no-label channel (PP 03/04/14)
					push @{$labelingInfo{$qChannelID}{'PARAMS'}},['None','No label','','','']; # No description for no-label channel
					$noLabelChannel=$qChannelID;
				}
			}
			$quantifAnnot="LABEL=SILAC::SOFTWARE=PD;$pdVersionComp";
			my $targetPos=0;
			foreach my $qChannelID (sort{$a<=>$b} keys %labelingInfo) {
				$targetPos++;
				$quantifAnnot.="::$targetPos;$labelingInfo{$qChannelID}{NAME};";
				my $count=0;
				foreach my $refParam (@{$labelingInfo{$qChannelID}{PARAMS}}) {
					$count++;
					$quantifAnnot.=join('#',@{$refParam}[0..3]);
					$quantifAnnot.='@' if $count < scalar @{$labelingInfo{$qChannelID}{PARAMS}};
				}
				$qChannel2TargetPos{$qChannelID}=$targetPos;
			}
#print "<BR>**$quantifAnnot**<BR>\n";
		}

		#### TMT / ITRAQ ###
		else {
			my ($qMethod)=($labelMethod =~/TMT/)? 'TMT':'ITRAQ';
			$quantifAnnot="LABEL=$qMethod"."::SOFTWARE=PD;$pdVersionComp";
			my $targetPos=0;
			foreach my $reporter (sort{$reporterIonList{$a}[0]<=>$reporterIonList{$b}[0]} keys %reporterIonList) {
				$targetPos++;
				$quantifAnnot.="::$reporterIonList{$reporter}[0];$reporter;$reporterIonList{$reporter}[1]";
				$qChannel2TargetPos{$reporterIonList{$reporter}[0]}=$targetPos;
			}
		}
		#$dbhLite->disconnect;
	}

	###>Looking for label-free extraction<###
	elsif ($msfFile) {
		my ($qNodeNum,$parentNodeStrg);
		if ($pdVersion >= 2.0) {
			# TODO : check this file!
			my ($wfXML) = $dbhLite->selectrow_array("SELECT WorkflowXML FROM Workflows");
			my $xml = new XML::Simple();
			my $xmlData = $xml->XMLin($wfXML);
			foreach my $processingNodeParameters (@{$xmlData->{WorkflowTree}->{WorkflowNode}}) {
				next unless $processingNodeParameters->{ProcessingNodeName} =~ /PrecursorIonsAreaDetectionNode|MinoraFeatureCreatorNode/; # for PD 2.2, the ProcessingNodeName has changed
				###> Label-free precursor Ion Area module of Proteome Discoverer
				($qNodeNum,$parentNodeStrg)=($processingNodeParameters->{ProcessingNodeNumber},$processingNodeParameters->{ParentProcessingNodeNumber});
				$quantifNodeNum=$qNodeNum;
				$quantifName='Precursor Ions Area Detector';
				$quantifAnnot="LABEL=FREE::SOFTWARE=PD;$pdVersionComp";
				$quantifAnnot.="::$processingNodeParameters->{ProcessingNodeParameters}{ProcessingNodeParameter}{FriendlyName}=$processingNodeParameters->{ProcessingNodeParameters}{ProcessingNodeParameter}{content}" if $pdVersion < 2.2;
			}
		}
		else {
			($qNodeNum,$parentNodeStrg)=$dbhLite->selectrow_array("SELECT ProcessingNodeNumber,ProcessingNodeParentNumber FROM ProcessingNodes WHERE NodeName='PrecursorIonsAreaDetectionNode'");
			if ($qNodeNum && $parentNodeStrg) {
				foreach my $parNodeNum (split(';',$parentNodeStrg)) {
					if ($parNodeNum==$searchNodeNum) {
						$quantifNodeNum=$qNodeNum;
						$quantifName='Precursor Ions Area Detector';
						$quantifAnnot="LABEL=FREE::SOFTWARE=PD;$pdVersionComp";
						my $sthParams=$dbhLite->prepare("SELECT FriendlyName,ValueDisplayString FROM ProcessingNodeParameters WHERE ProcessingNodeNumber=$quantifNodeNum");
						$sthParams->execute;
						while (my($paramName,$paramValue)=$sthParams->fetchrow_array) {
							$quantifAnnot.="::$paramName=$paramValue";
						}
						$sthParams->finish;
						last;
					}
				}
			}
		}
	}

	######################################
	####>Partially validated analysis<####
	######################################
	my %preValidatedProt;
	if ($validStatus==1) {
		print "&nbsp;&nbsp;&nbsp;- Looking for previously transferred data...";

		####>Deleting iRT data
		$dbh->do("DELETE FROM ANALYSIS_REFRT WHERE ID_ANALYSIS=$analysisID");
		foreach my $imgFile (glob("$promsPath{valid}/ana_$analysisID/regression_RT*.png")) {
			unlink $imgFile;
		}

		####>Deleting potential low-level quantification (2ndary quantification should prevent validation modification!!!)
		my $sthQ=$dbh->prepare("SELECT Q.ID_QUANTIFICATION,FOCUS,QUANTIF_ANNOT FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE Q.ID_QUANTIFICATION=AQ.ID_QUANTIFICATION AND ID_ANALYSIS=$analysisID ORDER BY Q.ID_QUANTIFICATION DESC");
		my $sthDelRT=$dbh->prepare("DELETE FROM QUANTIF_REFRT WHERE ID_QUANTIFICATION=?");
		#my $sthDelPep=$dbh->prepare("DELETE FROM PEPTIDE_QUANTIFICATION WHERE ID_QUANTIFICATION=?");
		my $sthDelProt=$dbh->prepare("DELETE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=?");
		my $sthDelAQ=$dbh->prepare("DELETE FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=?");
		my $sthDelQ=$dbh->prepare("DELETE FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
		$sthQ->execute;
		while (my ($quanID,$focus,$qAnnot)=$sthQ->fetchrow_array) {
			if ($focus eq 'protein') { # should not happen: Checked before allowing report
				$sthDelProt->execute($quanID);
				$sthDelAQ->execute($quanID);
				$sthDelQ->execute($quanID);
			}
			else { # peptide
				$sthDelRT->execute($quanID);
				#$sthDelPep->execute($quanID);
				rmtree("$promsPath{quantification}/proj_$projectID/quanti_$quanID") if -e "$promsPath{quantification}/project_$projectID/quanti_$quanID";
				#if ($qAnnot=~/EXTRACTION_ALGO/) { # MassChoQ
				#	$sthDelAQ->execute($quanID);
				#	$sthDelQ->execute($quanID);
				#}
				#else {
				#	$prevQuantifID=$quanID; # keep quantif entry
				#}
				$sthDelAQ->execute($quanID);
				$sthDelQ->execute($quanID);
			}
		}
		$sthQ->finish;
		$sthDelRT->finish;
		#$sthDelPep->finish;
		$sthDelProt->finish;
		$sthDelAQ->finish;
		$sthDelQ->finish;

		####>Fetching list already validated proteins from current analysis (if any)<####
		my $sthPV=$dbh->prepare("SELECT P.ID_PROTEIN,IDENTIFIER,ID_MASTER_PROTEIN FROM PROTEIN P,ANALYSIS_PROTEIN AP WHERE P.ID_PROTEIN=AP.ID_PROTEIN AND ID_ANALYSIS=$analysisID"); #,ID_CLUSTER
		$sthPV->execute;
		while (my ($protID,$identifier,$masterProtID)=$sthPV->fetchrow_array) { #,$clusID
			@{$preValidatedProt{$identifier}}=($protID,$masterProtID);
		}
		$sthPV->finish;
		if (scalar %preValidatedProt) {
			$prevValidation=1;

			####>Deleting peptides from both peptide tables and protein from ANALYSIS_PROTEIN(if any)<####
			$dbh->do("DELETE FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=$analysisID") || die $dbh->errstr;
			$dbh->do("DELETE FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_ANALYSIS=$analysisID") || die $dbh->errstr;
			$dbh->do("DELETE PEPTIDE_MODIFICATION FROM PEPTIDE_MODIFICATION INNER JOIN PEPTIDE ON PEPTIDE_MODIFICATION.ID_PEPTIDE=PEPTIDE.ID_PEPTIDE AND PEPTIDE.ID_ANALYSIS=$analysisID") || die $dbh->errstr;
			$dbh->do("DELETE FROM PEPTIDE WHERE ID_ANALYSIS=$analysisID") || die $dbh->errstr;
		}

		###>Deleting old peptide file (if any)<### transition to v2.7
		#unlink $oldPepFile if -e $oldPepFile;
		#if ($decoyFile) {
		#	unlink  $oldDecoyPepFile if -e $oldDecoyPepFile;
		#}
		#unlink $pepFile if $fileFormat=~/\.PDM/; # reaxtraction will be necessary

		print " Done.<BR>\n";
	}

	#####################################
	####>Fetching all validated data<#### from tables PROTEIN_VALIDATION, RANK_PROTEIN_MATCH and QUERY_VALIDATION
	##################################### Starting with proteins because 100% of sel_status=1 will be kept (!= queries:orphans are possible)
	print "&nbsp;&nbsp;&nbsp;- Fetching new validated data...";

	my %spectrum2fileID;
	if ($msfFile && $numMergedFiles > 1) {
		my $queryGetDataSrc=($pdVersion >=2.2)?"SELECT SpectrumID,SpectrumFileID FROM MSnSpectrumInfo" : "SELECT SpectrumID,FileID FROM SpectrumHeaders SH,MassPeaks MP WHERE SH.MassPeakID=MP.MassPeakID";
		my $sthDataSrc=$dbhLite->prepare($queryGetDataSrc);
		$sthDataSrc->execute;
		while (my ($spectrumID,$fileID)=$sthDataSrc->fetchrow_array) {
			$spectrum2fileID{$spectrumID}=$fileID;
		}
		$sthDataSrc->finish;
	}

	####>Reference RT<####
	my ($refRT_ID,%anaPepRefRTdata,%spectrum2seqRT);
	my ($skipRefRTidString,$skipRefRTidentString)=('','');
	my (%skipRefRTidStrg,%skipRefRTidentStrg,%numRefPepMatched);
	foreach my $refID (sort{$b<=>$a} keys %referenceRTdata) { # most recent first
		$referenceRTdata{$refID}{'QUERY'}->execute($analysisID);
		while (my ($refRTprotID,$refRTidentifier)=$referenceRTdata{$refID}{'QUERY'}->fetchrow_array) { # Match!
			if (!$skipRefRTidStrg{$refID}) {
				#$refRT_ID=$refID;
				$skipRefRTidStrg{$refID}="AND ID_PROT_VALID NOT IN ($refRTprotID";
				$skipRefRTidentStrg{$refID}="AND IDENTIFIER NOT IN ('$refRTidentifier'";
			}
			else {
				$skipRefRTidStrg{$refID}.=",$refRTprotID";
				$skipRefRTidentStrg{$refID}.=",'$refRTidentifier'";
			}

			##>Fetching reference peptides data
			my %usedDataSources; # PD
			$sthPepRT->execute($analysisID,$analysisID,$refRTidentifier);
			while (my ($qNum,$spectrumID,$massData,$charge,$elutionData,$pepInfo)=$sthPepRT->fetchrow_array) { # order by desc pep score (only best is recorded)
				next if ($msfFile && !$spectrumID); # PD
				$spectrumID=-$qNum unless $spectrumID; # non-PD
				my ($pepSeq)=($pepInfo=~/SEQ=(\w+)/);
				next unless defined $referenceRTdata{$refID}{'PEPTIDES'}{$pepSeq};
				my ($massExp)=($massData=~/EXP=([^,]+)/);
				next unless abs($massExp-$referenceRTdata{$refID}{'PEPTIDES'}{$pepSeq}[2]) < 1; # delta massExp (to exclude modified pep)
				my ($pepScore)=($pepInfo=~/SC=([^,]+)/);
				next if (defined($minScorePep) && $pepScore < $minScorePep);
				my $fileID;
				if ($msfFile) { # PD
					$fileID=($numMergedFiles > 1)? $spectrum2fileID{$spectrumID} : (keys %mergedFiles)[0];
					next if ($usedDataSources{$pepSeq} && $usedDataSources{$pepSeq}{$fileID}); # filtering duplicate identification for same source (no filtering on SEL=1 because of possible merge search)
				}
				my ($massCalc)=($pepInfo=~/CALC=([^,]+)/);
				my ($deltaMass)=($pepInfo=~/DELT=([^,]+)/);
				my ($massObs)=($massData=~/OBS=([^,]+)/);
				@{$anaPepRefRTdata{$pepSeq}{$spectrumID}}=($qNum,$charge,$massExp,$massObs,$massCalc,$deltaMass,$pepScore,$elutionData);
				$spectrum2seqRT{$spectrumID}=$pepSeq;
				$usedDataSources{$pepSeq}{$fileID}=1 if $msfFile; # PD
				$numRefPepMatched{$refID}++;
			}
		}
		if ($skipRefRTidStrg{$refID}) {
			$skipRefRTidStrg{$refID}.=')';
			$skipRefRTidentStrg{$refID}.=')';
			#print "[Retention time reference '$referenceRTdata{$refRT_ID}{NAME}' ($referenceRTdata{$refRT_ID}{IDENTIFIER}) detected]";
		}
		#last; # only 1 type of ref RT allowed?
	}
	if (scalar keys %numRefPepMatched) { # keep only best match (because a referenceRT identifier ambiguity) # only 1 type of ref RT allowed?
		$refRT_ID=(sort{$numRefPepMatched{$b}<=>$numRefPepMatched{$a}} keys %numRefPepMatched)[0];
		my $totRefPep=scalar keys %{$referenceRTdata{$refRT_ID}{'PEPTIDES'}};
		print "[Retention time reference '$referenceRTdata{$refRT_ID}{NAME}' ($referenceRTdata{$refRT_ID}{IDENTIFIER} with $numRefPepMatched{$refRT_ID}/$totRefPep matches) detected]";
		$skipRefRTidString=$skipRefRTidStrg{$refRT_ID};
		$skipRefRTidentString=$skipRefRTidentStrg{$refRT_ID};
	}

	####>Fetching all matched proteins from table PROTEIN_VALIDATION
	my $refScore=($msType eq 'PMF')? 'MAX_SCORE' : 'SCORE'; # take max score (score=0) if PMF ... Is it appropriate?...
	my $sthV=$dbh->prepare("SELECT IDENTIFIER,PROT_DES,ORGANISM,MW,$refScore,CONF_LEVEL,PROT_LENGTH,DB_RANK FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$analysisID AND SEL_STATUS>=1 AND IDENTIFIER NOT LIKE 'DECOY_%' $skipRefRTidString");
	$sthV->execute || die $sthV->errstr;
	my $refValidData=$sthV->fetchall_arrayref; # reference to an array
	$sthV->finish;
	print '.';

	####>Fetching all validated queries and associated peptides from table QUERY_VALIDATION
	my $selectString='ID_QUERY,QUERY_NUM,VALID_STATUS,MASS_DATA,CHARGE,ELUTION_TIME,EXT_SPECTRUMID';
	foreach my $r (1..$maxRank) {$selectString.=",INFO_PEP$r";}
	my $minQueryValStatus=($labelMethod=~/SILAC/i)? -3 : 1;
	my $sthQ=$dbh->prepare("SELECT $selectString FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID AND VALID_STATUS>=$minQueryValStatus AND QUERY_NUM>0"); # also contains refRT peptide data (no impact: filtered out by $queryList{$queryNum}{$rank})
	$sthQ->execute || die $sthQ->errstr;
	my $refQueryData=$sthQ->fetchall_arrayref; # (reference to an array)  Validated queries only!
	$sthQ->finish;
	print '.';

	####>Fetching starting ID and protecting ID space for table PROTEIN
	$dbh->do("DELETE FROM PROTEIN WHERE ID_PROJECT=$projectID AND IDENTIFIER is NULL") || $dbh->errstr; # in case of error in previous process (rare!)
	####>Auto-increment on 2018/03/20
	#my ($proteinID)=$dbh->selectrow_array("select max(ID_PROTEIN) from PROTEIN");
	#$proteinID=0 unless $proteinID;
	#my $maxProteinID=$proteinID + scalar @{$refValidData} + 1;
	#$dbh->do("INSERT INTO PROTEIN (ID_PROTEIN,ID_PROJECT) VALUES ($maxProteinID,$projectID)") || $dbh->errstr;
	#$dbh->commit;
    my $proteinID=-1;

	####>Fetching starting ID  and protecting ID space for table PEPTIDE
	$dbh->do("DELETE FROM PEPTIDE WHERE ID_ANALYSIS=$analysisID AND PEP_SEQ IS NULL") || $dbh->errstr; # in case of error in previous process (rare!)
	#print '.';

	####>Fetching query_num,pep_rank matching selected proteins from table RANK_PROTEIN_MATCH
	my %matchData; # All non-decoy data (validated or not)!
	my $sthAM=$dbh->prepare("SELECT IDENTIFIER,QUERY_NUM,PEP_RANK,MATCH_INFO FROM RANK_PROTEIN_MATCH WHERE ID_ANALYSIS=$analysisID AND QUERY_NUM>0 $skipRefRTidentString");
	$sthAM->execute;
	while (my ($identifier,@data) = $sthAM->fetchrow_array) {
		push @{$matchData{$identifier}},\@data;
	}
	$sthAM->finish;
	print '.';

	my (%proteinList,%protIdList, # validated protein lists
		%queryList # All queries (validated or not) matching validated proteins
		);
	foreach my $refProt (@{$refValidData}) { # validated proteins
		my $identifier=$$refProt[0];
		if (!$projectProt{$identifier}) { # protein is new to project=> create new id
			#$proteinList{$identifier}=++$proteinID;
			#$protIdList{$proteinID}=$identifier;
			$proteinList{$identifier}=--$proteinID;
		}
		else { # use existing id (no duplicates)
			$proteinList{$identifier}=$projectProt{$identifier};
			$protIdList{$projectProt{$identifier}}=$identifier;
		}
		#foreach my $refData (@{$matchData{$identifier}}) {
		#	my ($queryNum,$rank,$matchInfo)=@{$refData};
		#	push @{$queryList{$queryNum}{$rank}},"$proteinList{$identifier}:$matchInfo";
		#}
	}
	print " Done.<BR>\n";

	#################################################################
	####>Inserting data into tables PROTEIN and ANALYSIS_PROTEIN<####
	#################################################################

	####>Inserting data into PROTEIN and ANALYSIS_PROTEIN
	my %sequenceList;
	my %proteinScore;
	my %newValidatedProt;
	my %proteinAnnotQual;

	###########################################################
	####>Fetching protein sequences from sequence dataBank<####
	###########################################################
	print "&nbsp;&nbsp;&nbsp;- Extracting protein sequences from sequence databank...";
	####>Fetching protein Databank info...
	if (scalar keys %proteinList) {
		if (-e "$promsPath{valid}/ana_$analysisID/analysis.fasta") {
			open (FAS,"$promsPath{valid}/ana_$analysisID/analysis.fasta") || die "Fasta file not found!\n";
			my (@matchedProt,$matched);
			while (<FAS>) {
				if (/^>(.+)\n/) {
					@matchedProt=();
					$matched=0;
					foreach my $entryStrg (split(//,$1)) {
						my ($identifier)=($entryStrg=~/^(\S+)/);
						if ($proteinList{$identifier} && (!$projectProt{$identifier} || $incompleteProtein{$identifier})) {
							push @matchedProt,$identifier;
							$matched=1;
						}
					}
					next;
				}
				if ($matched) {
					next if /^#/; # skip ##END when fasta is generated by distant Mascot server (myproms4databanks.pl)
					chomp;
					foreach my $identifier (@matchedProt) {
						$sequenceList{$identifier}.=$_;
					}
				}
			}
			close FAS;
		}
		else {print 'Fasta file not found! Skipping...';}
	}
	print '.';

	foreach my $identifier (keys %proteinList) {
		next if ($projectProt{$identifier} && !defined($incompleteProtein{$identifier}));
		$sequenceList{$identifier}='-' unless $sequenceList{$identifier};# protein not found in fasta file
	}
	print " Done.<BR>\n";

	print "&nbsp;&nbsp;&nbsp;- Recording new proteins...";

	my $sthAnnot=$dbh->prepare("SELECT PROT_DES,ORGANISM,MW,PROT_LENGTH FROM PROTEIN WHERE IDENTIFIER=? AND PROT_LENGTH>0 ORDER BY ID_PROTEIN DESC LIMIT 0,1;");
	#my $sthProt=$dbh->prepare("INSERT INTO PROTEIN (ID_PROTEIN,ID_PROJECT,IDENTIFIER,ALIAS,PROT_DES,ORGANISM,MW,PROT_SEQ,PROT_LENGTH) VALUES (?,$projectID,?,?,?,?,?,?,?)");
	my $sthProt=$dbh->prepare("INSERT INTO PROTEIN (ID_PROJECT,IDENTIFIER,ALIAS,PROT_DES,ORGANISM,MW,PROT_SEQ,PROT_LENGTH) VALUES ($projectID,?,?,?,?,?,?,?)");
	my $sthAna=$dbh->prepare("INSERT INTO ANALYSIS_PROTEIN (ID_ANALYSIS,ID_PROTEIN,SCORE,CONF_LEVEL,DB_RANK) VALUES ($analysisID,?,?,?,?)");
	my $sthUpdProt=$dbh->prepare("UPDATE PROTEIN SET MW=?,PROT_SEQ=?,PROT_LENGTH=? WHERE ID_PROTEIN=?");
	my $sthUpdDesc=$dbh->prepare("UPDATE PROTEIN SET PROT_DES=?,ORGANISM=? WHERE ID_PROTEIN=?");
	my $newProteins=0;
	foreach my $refProt (@{$refValidData}) {
		my ($identifier,$des,$organism,$mw,$score,$confLevel,$protLength,$dbRank)=@{$refProt};
		$newValidatedProt{$identifier}=1;
		$confLevel=2 unless $confLevel;
		if (!$protLength) { # missing annotation but protein is in project
			$sthAnnot->execute($identifier);
			my ($dbDes,$dbOrg,$dbMW,$dbProtLength)=$sthAnnot->fetchrow_array;
			($des,$organism,$mw,$protLength)=($dbDes,$dbOrg,$dbMW,$dbProtLength) if $dbProtLength;
		}
		my $alias=($projectIdentMapID==$giIdentID && $identifier=~/(gi\|\d+)/)? $1 : $identifier; # "GI" restriction on alias
		unless ($projectProt{$identifier}) { # protein is new to project => insert
			$sthProt->execute($identifier,$alias,$des,$organism,$mw,$sequenceList{$identifier},$protLength) || die $sthProt->errstr;
			my ($proteinID)=$dbh->last_insert_id(undef, undef, 'PROTEIN', 'ID_PROTEIN');
			$proteinList{$identifier}=$proteinID;
			$protIdList{$proteinID}=$identifier;
			$newProteins=1; # new protein(s) added to project -> identifier mapping
		}
		foreach my $refData (@{$matchData{$identifier}}) {
			my ($queryNum,$rank,$matchInfo)=@{$refData};
			push @{$queryList{$queryNum}{$rank}},"$proteinList{$identifier}:$matchInfo";
		}
		#my $alias=($identifierTypes{$dbRank})? &promsMod::modifyIdentifier($identifier,$identModString{$identifierTypes{$dbRank}}) : $identifier;
		$proteinScore{$proteinList{$identifier}}=$score;
		$proteinLength{$proteinList{$identifier}}=$protLength;
		$protSpeciesClass{$proteinList{$identifier}}=($preferredSpecies && $organism eq $preferredSpecies)? 1 : 0;
		#$proteinAnnotQual{$proteinList{$identifier}}=($des=~/no\sdescription/i)? 3 : ($des=~/unnamed/i || $des=~/unknown/i)? 2 : ($des=~/hypothetical/i)? 1 : 0;  #fp
		# Annotation quality: the smaller the better:
		# 1: Swiss-Prot (EZRI_HUMAN) (<up to 5digits>_<species>)
		# 2: TrEMBL (Q12345_HUMAN) (<6 digits>_<species>)
		# 3: any other identifier with full desc
		# 4: any other identifier with 'hypothetical' in desc
		# 5: any other identifier with 'unknown'/'unnamed' in desc
		# 6: any other identifier with no desc
		$proteinAnnotQual{$proteinList{$identifier}}=($identifier=~/\b\w{1,5}_/)? 1 : ($identifier=~/\b\w{6}_/)? 2 : ($des=~/no\sdescription/i)? 6 : ($des=~/unnamed/i || $des=~/unknown/i)? 5 : ($des=~/hypothetical/i)? 4 : 3;  #PP

		if ($incompleteProtein{$identifier}) {# protein is in project, but incomplete, file MW, Sequence and length
			$sthUpdProt->execute($mw,$sequenceList{$identifier},$protLength,$incompleteProtein{$identifier});
			if ($noDescriptionProtein{$identifier} && $des && $des !~ /no\sdescription/i) {
				$sthUpdDesc->execute($des,$organism,$noDescriptionProtein{$identifier});
				delete $noDescriptionProtein{$identifier}; # remove from list
			}
			delete $incompleteProtein{$identifier} if $protLength; # remove from list
		}
		$sthAna->execute($proteinList{$identifier},$score,$confLevel,$dbRank) || die $sthAna->errstr;
	}
	$sthProt->finish;
	$sthAna->finish;
	$sthUpdDesc->finish;
	$sthUpdProt->finish;
	push @newAnaMapping,$analysisID if $newProteins;
	print " Done.<BR>\n";

	#######################################################################
	####>Inserting data into tables PEPTIDE and PEPTIDE_PROTEIN_ATTRIB<####
	#######################################################################
	print "&nbsp;&nbsp;&nbsp;- Recording new peptides...";
	#####>Fetching starting ID  and protecting ID space for table PEPTIDE
	#$dbh->do("DELETE FROM PEPTIDE WHERE ID_ANALYSIS=$analysisID AND PEP_SEQ IS NULL") || $dbh->errstr; # in case of error in previous process (rare!)
	#my ($peptideID)=$dbh->selectrow_array("SELECT MAX(ID_PEPTIDE) FROM PEPTIDE");
	#my $maxPeptideID=$peptideID + 2*$maxRank*(scalar @{$refQueryData}) + 1; # much larger than actual ID space used (but simpler)
	#$dbh->do("INSERT INTO PEPTIDE (ID_PEPTIDE,ID_ANALYSIS) VALUES ($maxPeptideID,$analysisID)") || $dbh->errstr;
	#$dbh->commit;

	#my %spectrum2fileID;
	#if ($msfFile && $numMergedFiles > 1) {
	#	my $sthDataSrc=$dbhLite->prepare("SELECT SpectrumID,FileID FROM SpectrumHeaders SH,MassPeaks MP WHERE SH.MassPeakID=MP.MassPeakID");
	#	$sthDataSrc->execute;
	#	while (my ($spectrumID,$fileID)=$sthDataSrc->fetchrow_array) {
	#		$spectrum2fileID{$spectrumID}=$fileID;
	#	}
	#	$sthDataSrc->finish;
	#}
	my %usedModifs;
	my $sthAnaMod=$dbh->prepare("SELECT ID_MODIFICATION,SPECIFICITY FROM ANALYSIS_MODIFICATION WHERE ID_ANALYSIS=$analysisID");
	$sthAnaMod->execute;
	while (my ($modID,$specificity)=$sthAnaMod->fetchrow_array) {
		$usedModifs{$modID}=$specificity;
	}
	$sthAnaMod->finish;

	my $insertPepString = 'ID_ANALYSIS,PEP_SEQ,PEP_LENGTH,QUERY_NUM,PEP_RANK,SEARCH_RANK,SCORE,MISS_CUT,';
	$insertPepString .= 'MR_EXP,MR_CALC,MR_OBS,MR_DELTA,SUBST,CHARGE,ELUTION_TIME,VALID_STATUS,SPEC_COUNT,COMMENTS,DATA';
	my $sthInsPep=$dbh->prepare("INSERT INTO PEPTIDE ($insertPepString) VALUES ($analysisID,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
	my $sthGetQueryMod=$dbh->prepare("SELECT ID_MODIFICATION,POS_STRING,REF_POS_STRING FROM QUERY_MODIFICATION WHERE ID_QUERY=? AND PEP_RANK=?");
	my $sthInsPepMod=$dbh->prepare("INSERT INTO PEPTIDE_MODIFICATION (ID_MODIFICATION,ID_PEPTIDE,MODIF_TYPE,POS_STRING,REF_POS_STRING) VALUES (?,?,?,?,?)");
	my $sthModName=$dbh->prepare("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES FROM MODIFICATION WHERE ID_MODIFICATION=?");

	my $insertAttString='ID_PEPTIDE,ID_PROTEIN,PEP_BEG,PEP_END,ID_ANALYSIS,FLANKING_AA';
	my $sthAtt=$dbh->prepare("INSERT INTO PEPTIDE_PROTEIN_ATTRIB ($insertAttString) VALUES (?,?,?,?,$analysisID,?)");
	my $sthUpAtt=$dbh->prepare("UPDATE PEPTIDE_PROTEIN_ATTRIB SET IS_SPECIFIC=1 WHERE ID_PEPTIDE=?");

	#my $sthGetData = $dbh->prepare("SELECT DATA FROM PEPTIDE WHERE ID_PEPTIDE=?");
	#my $sthUpData = $dbh->prepare("UPDATE PEPTIDE SET DATA=? WHERE ID_PEPTIDE=?");

	my (%validQueryData,%peptideInfo,%spectrum2PeptideIDs,%modifName2ID,%labelModifID,%qChannelID2ModifID); # needed for quanti  ###%pepSequence,%seqMissedCuts ,%validPeptideData,%unlabeledVarMods,%varModChannel
	# my %pepSeqList; # records all peptide sequences to prevent transmitting duplicates (value='peptideID:score')
	my %matchList; # records list of pepID matched by protID
	my %numMatch; # records total number of matching occurence / protein (>=num peptides)
	my %pepProteins; # records number of proteins matched by a peptide
	my %specificSequences; # records list of protein-specific peptide sequences (needed for SILAC ghost peptides)
	my %boundaryStatus; # records peptides boundaries
	my (%maxCovLength,%seqEndMatched); # used for coverage: Could be bigger than protLength due to isoforms with same identifier
	my $manualValidation=0; # Flag for manual validation
	my $isLabelPD=($labelMethod=~/SILAC/i)? 1 : 0; # Proteome Discoverer flag to extract lower-scoring as ghost/hidden peptides
	my $count=0;
	foreach my $refQuery (@{$refQueryData}) {
		my ($queryID,$queryNum,$qValStatus,$massData,$charge,$elutionData,$spectrumID,@pepInfo)=@{$refQuery};
		if ($fileFormat eq 'PARAGON.XML' && $labelMethod && $elutionData=~/sp([\d\.]+)/) {
			$spectrumID=$1;
		}
		$spectrumID=0 unless $spectrumID;
		$elutionData='' unless $elutionData;
		$charge='' unless $charge;
		my $rank=0;
		foreach my $info (@pepInfo) {
			$rank++;
			next unless $queryList{$queryNum}{$rank}; # no matched protein selected or orphan peptide (in case protein was excluded)
			my ($select)=($info=~/SEL=(-*\d)/);
			#next unless $select; # peptide is valid (selectable + selected)
			next if ($select <= -2 || (!$isLabelPD && $select <= 0));
			$manualValidation=1 if $select=~/(2|-3)/; # detection if at least 1 selection (validStatus>0)!!! (rejected queries will be scanned later)
			my $pepValStat=($select >= 0)? $select : 0; # ghost -> 0
			my ($massExp)=($massData=~/EXP=(\d+\.*\d*)/);
			my ($massObs)=($massData=~/OBS=(\d+\.*\d*)/);
			my ($searchRank)=($info=~/SRK=(\d+)/); $searchRank=$rank unless $searchRank; # mixed decoy DB
			my ($missCut)=($info=~/MIS=(\d)/);
			my ($seq)=($info=~/SEQ=(\w+)/);
			my ($score)=($info=~/SC=(-?\d+\.*\d*)/);
			next if (defined($minScorePep) && $score < $minScorePep); # in case ghost allowed ($select=0|-1)
			my ($delta)=($info=~/DELT=(-*\d+\.*\d*)/);
			my ($massCalc)=($info=~/CALC=(\d+\.*\d*)/);
			#my ($varMod)=&promsMod::toStringVariableModifications($dbh,'QUERY',$queryID,$analysisID,$seq,$rank);
			my ($sub)=($info=~/SUBST=([^,]+)/); $sub='' unless $sub;
			my ($specCount)=($info=~/SPC=(\d+)/);
			my ($comments)=($info=~/COM=(.+)/); $comments=~s/,PRS.+// if $comments; # COM comes at the end except in case of PhosphoRS data
			my $extraData;
			if ($rank==1) {
				($extraData)=($info=~/(PEP=.+,QVAL=[^,]+)/); # Qvality data
			}
			#my ($refMod)=($info=~/RMOD=([^,]+)/); # reference var_mod before user edition
			#if ($refMod) {
			#	$extraData.='##' if $extraData;
			#	$extraData.="REF_MOD=$refMod";
			#}
			my ($PRS)=($info=~/(PRS=[^,]+)/); # could be undef
			if ($PRS) {
				$extraData.='##' if $extraData;
				$extraData.=$PRS;
			}
			if ($spectrumID) {
				$extraData.='##' if $extraData;
				$extraData.="EXT_SPECTRUMID=$spectrumID";
				if ($numMergedFiles > 1) {# && !$labelMethod file rank (will be added later if $labelMethod)
					$extraData.="##SOURCE_RANK=$mergedFiles{$spectrum2fileID{$spectrumID}}[1]";
				}
			}
			#my $data = '';
			#if($PRS){
			#	$data = $PRS;
			#}

			$sthInsPep->execute($seq,length($seq),$queryNum,$rank,$searchRank,$score,$missCut,$massExp,$massCalc,$massObs,$delta,$sub,$charge,$elutionData,$pepValStat,$specCount,$comments,$extraData) || die $sthInsPep->errstr;
			my $peptideID=$dbh->last_insert_id(undef,undef,'PEPTIDE','ID_PEPTIDE');
			$sthGetQueryMod->execute($queryID,$rank);
			my @pepVmods=();
			while (my ($modificationID,$posStr,$refPosStr) = $sthGetQueryMod->fetchrow_array) {
				push @pepVmods,[$modificationID,$posStr,$refPosStr];
				$sthInsPepMod->execute($modificationID,$peptideID,'V',$posStr,$refPosStr) || die $sthInsPepMod->errstr;
				#$usedModifs{$modificationID}=1; # just to be safe: already generated from ANALYSIS_MODIFICATION table
			}
			if ($extraData && $labelMethod=~/SILAC/i) {$peptideExtraData{$peptideID}=$extraData;} # record for later update of DATA field in SILAC quanti
			my $ghostFlag=($pepValStat)? 1 : -1; # beg/end recorded as neg values for ghost peptides
			foreach my $match (@{$queryList{$queryNum}{$rank}}) {
				my ($prID,@matchInfo)=split(/:/,$match);
				$matchList{$prID}{$peptideID}=1 if $pepValStat; # not for ghost!
				$pepProteins{$peptideID}{$prID}=1;
				my %usedBeg; # Temporary: Prevents duplicate match from generating errors in DB (should be corrected during Ana import)
				foreach my $aaStr (@matchInfo) {
					my ($beg,$flAA_Nter,$flAA_Cter)=split(/,/,$aaStr);
					next if $usedBeg{$beg};
					$flAA_Nter='' unless $flAA_Nter;
					$flAA_Cter='' unless $flAA_Cter;
					my $end=$beg+length($seq)-1;
					if ($pepValStat) { # not for ghost!
						$boundaryStatus{$prID}{$beg}++;
						$boundaryStatus{$prID}{$end}--;
						$maxCovLength{$prID}=$end if (!$maxCovLength{$prID} || $maxCovLength{$prID} < $end);
						$seqEndMatched{$prID}=1 if $flAA_Cter eq '-'; # peptide matches end of sequence
						$numMatch{$prID}++; # not for ghost!
					}
					$sthAtt->execute($peptideID,$prID,$beg*$ghostFlag,$end*$ghostFlag,"$flAA_Nter$flAA_Cter") || die $sthAtt->errstr;
					$usedBeg{$beg}=1;
				}
			}

			##>Peptide is specific ?
			if (scalar keys %{$pepProteins{$peptideID}}==1) {
				$sthUpAtt->execute($peptideID);
				$specificSequences{$seq}=1 if ($labelMethod=~/SILAC/i); # needed later for ghost peptides
			}

			##>Quantif
			if ($labelMethod) {
				push @{$validQueryData{$queryNum}{'pepID'}},$peptideID; # needed for QUANTI (not for PARAGON)

				if ($msfFile && $labelMethod=~/SILAC|TMT|iTRAQ/i) { # SILAC and TMT
					@{$peptideInfo{$peptideID}}=($seq,length($seq),$missCut,$charge,\@pepVmods); # needed later for quantif data
					###$pepSequence{$peptideID}=$seq; # needed later for quantif data
					###$seqMissedCuts{$seq}=$missCut; # needed later for ghost peptides
#push @{$spectrum2PeptideIDs{$seq}{$spectrumID}},[$peptideID,$massObs,\@pepVmods]; # ($varMod) multiple pepID (with diff massObs are possible per SpectrumID!) <-- Updated PP 27/08/15: multiple MS/MS (SpectrumID) are linked to the same quanResultID but w/ diff mass & channel: Light(/Medium)/Heavy (PD!)
if ($rank==1) {
	@{$spectrum2PeptideIDs{$spectrumID}}=($peptideID,$score,$pepValStat,$massObs); #,\@pepVmods 1 spectrumID <=> 1 myProMS peptideID
}
				}
				elsif ($fileFormat eq 'PARAGON.XML' && $labelMethod=~/iTRAQ/i) { # PARAGON data structure != MSF
					push @{$spectrum2PeptideIDs{$spectrumID}{'PEP_ID'}},$peptideID;
				}
			}
			elsif ($msfFile && $quantifNodeNum) { # Label-free with PD ($labelMethod is undef)
				push @{$spectrum2PeptideIDs{$spectrumID}{'PEP_ID'}},$peptideID;
			}
		}
		$count++;
		if ($count>=1000) {
			print '.';
			$count=0;
		}
	}
	if ($labelMethod) {
		# Store modName <=> mod ID association
		foreach my $modifID (keys %usedModifs) {
			$sthModName->execute($modifID);
			my ($psiName,$interName,$synomymes)=$sthModName->fetchrow_array;
			my $modifName=($psiName)? $psiName : ($interName)? $interName : (split(/##/,$synomymes))[1];
			$modifName2ID{$modifName}=$modifID;
#print "# $modifName -> $modifID<BR>\n";
		}
		# Remember ID of labeling modifs
		foreach my $qChannelID (keys %labelingInfo) {
			%{$qChannelID2ModifID{$qChannelID}}=();
			foreach my $refParam (@{$labelingInfo{$qChannelID}{'PARAMS'}}) {
#print "## $qChannelID: '",join('::',@{$refParam}),"'<BR>\n";
				next if $refParam->[1] eq 'No label';
				#$labelModifID{$modifName2ID{$refParam->[3]}}=$qChannelID;
				if (!$modifName2ID{$refParam->[1]} || ($modifName2ID{$refParam->[1]} && $refParam->[2] !~ $usedModifs{$modifName2ID{$refParam->[1]}})) { # Somehow, this modification used for label-quantification was not used for identification
						print "<BR>\n<FONT color=\"#DD0000\">ERROR! Wrong label quantification!</FONT>";
						print "<BR>\nLabel modification <B>$refParam->[1]</B> on $refParam->[2] is used for quantification but has not been set as a dynamic modification for this residue.";
						exit;
				}
				$labelModifID{$modifName2ID{$refParam->[1]}}=$qChannelID; # [1] and no longer [3] because in storeAnalyses.cgi now records true unimod modif name
				$modifName2ID{$refParam->[3]}=$modifName2ID{$refParam->[1]}; # extend modID mapping to custom label modif name if any (because $sthVMIAll returns custom name!)
				push @{$qChannelID2ModifID{$qChannelID}{$modifName2ID{$refParam->[1]}}},$refParam->[2]; # modifRes
			}
		}
	}
	$sthGetQueryMod->finish;
	$sthInsPep->finish;
	$sthInsPepMod->finish;
	$sthAtt->finish;
	$sthUpAtt->finish;
	$sthModName->finish;

	#$sthGetData->finish;
	#$sthUpData->finish;

	print " Done.<BR>\n";

	###################################
	####>Updating ANALYSIS_PROTEIN<####
	###################################
	####>Finding number of times proteins are at top of match group hierarchy
	print "&nbsp;&nbsp;&nbsp;- Computing proteins visibility...";
	my (%numProtTop,%bestProtVis,%numProtPeptides,%proteinPepDens);
	my $sthTop=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=? AND VISIBILITY=2");
	my $sthBestVis=$dbh->prepare("SELECT MAX(VISIBILITY) FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=?");
	foreach my $protID (keys %matchList) {
		if ($projectProt{$protIdList{$protID}}) { # proteins are already in project
			$sthTop->execute($protID);
			($numProtTop{$protID})=$sthTop->fetchrow_array;
			if ($protVisibility) {
				$sthBestVis->execute($protID);
				$bestProtVis{$protID}=$sthBestVis->fetchrow_array;
			}
		}
		else { # proteins are new in project
			$numProtTop{$protID}=0;
			$projectProt{$protIdList{$protID}}=$protID; # add protein to project list
		}
		#>Peptide density
		$numProtPeptides{$protID}=scalar (keys %{$matchList{$protID}});
		$proteinPepDens{$protID}=$proteinLength{$protID}/$numProtPeptides{$protID};
	}
	$sthTop->finish;
	$sthBestVis->finish;
	print " Done.<BR>\n";

	####>Finding match groups and protein visibility
	print "&nbsp;&nbsp;&nbsp;- Building match groups...";
	# Rules of ordering: numPep > num@top > score > delta length > annotQuality > length > identifier
	#my @sortedProtIDs=sort{scalar (keys %{$matchList{$b}})<=>scalar (keys %{$matchList{$a}}) || $numProtTop{$b}<=>$numProtTop{$a} || $proteinScore{$b}<=>$proteinScore{$a} || $proteinAnnotQual{$a} <=> $proteinAnnotQual{$b} || $proteinLength{$a}<=>$proteinLength{$b} || $a<=>$b} keys %matchList;
	my @sortedProtIDs = sort{$numProtPeptides{$b}<=>$numProtPeptides{$a} || $protSpeciesClass{$b}<=>$protSpeciesClass{$a} || $numProtTop{$b}<=>$numProtTop{$a} || $proteinScore{$b}<=>$proteinScore{$a} || &deltaLength($proteinLength{$a},$proteinLength{$b},$proteinPepDens{$a},$proteinPepDens{$b}) || $proteinAnnotQual{$a}<=>$proteinAnnotQual{$b} || $proteinLength{$a}<=>$proteinLength{$b} || $a<=>$b} keys %matchList;
	my (%matchGroup, %visibility);
	&promsMod::createMatchGroups(\%matchList, \%matchGroup, \@sortedProtIDs, \%visibility, \%bestProtVis, $protVisibility, 1);

	
	####>Computing PEP_SPECIFICITY
	my %pepSpecificity;
	foreach my $peptideID (keys %pepProteins) {
		my $specificity=sprintf "%.2f",100/(scalar keys %{$pepProteins{$peptideID}});
		$specificity*=1; # 100.00 -> 100
		foreach my $protID (keys %{$pepProteins{$peptideID}}) {
			$pepSpecificity{$protID}=$specificity if (!$pepSpecificity{$protID} || $pepSpecificity{$protID}<$specificity);
		}
	}

	####>Computing PEP_COVERAGE
	my %pepCoverage;
	foreach my $protID (keys %boundaryStatus) {
		next unless $proteinLength{$protID};
		my $coverage=0;
		my $hasPeptide=0;
		my $boundaryNter=0;
		foreach my $boundary (sort{$a<=>$b} keys %{$boundaryStatus{$protID}}) {
			if ($hasPeptide==0) { # start of peptide region (cannot become <0!)
				$boundaryNter=$boundary;
			}
			$hasPeptide+=$boundaryStatus{$protID}{$boundary};
			if ($hasPeptide==0) { # end of peptide region (should be true for last boundary too)
				$coverage+=($boundary-$boundaryNter)+1;
			}
		}
		my $usedProtLength=($maxCovLength{$protID} <= $proteinLength{$protID})? $proteinLength{$protID} : ($seqEndMatched{$protID})? $maxCovLength{$protID} : -1*$maxCovLength{$protID}; # -: flag for protLength problem
		$pepCoverage{$protID}=sprintf "%.1f",(100*$coverage)/$usedProtLength;
		$pepCoverage{$protID}*=1; # 25.0 -> 25
	}

	####>Updating table
	print "&nbsp;&nbsp;&nbsp;- Recording proteins data...";
	my $sthUpAP=$dbh->prepare("UPDATE ANALYSIS_PROTEIN SET VISIBILITY=?,MATCH_GROUP=?,NUM_PEP=?,NUM_MATCH=?,PEP_SPECIFICITY=?,PEP_COVERAGE=? WHERE ID_PROTEIN=? AND ID_ANALYSIS=$analysisID");
	foreach my $protID (keys %matchList) {
		$sthUpAP->execute($visibility{$protID},$matchGroup{$protID},scalar keys %{$matchList{$protID}},$numMatch{$protID},$pepSpecificity{$protID},$pepCoverage{$protID},$protID);
	}
	$sthUpAP->finish;
	print " Done.<BR>\n";


	####################################################################
	####>Deleting old prevalidated proteins that were not confirmed<####
	####################################################################
	print "&nbsp;&nbsp;&nbsp;- Deleting obsolete protein data...";
	my %modifMasterProteins;
	my $sthNumAna=$dbh->prepare("SELECT 1 FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=? LIMIT 0,1");
	my $sthDelMS=$dbh->prepare("DELETE MS FROM MODIFICATION_SITE MS INNER JOIN CATEGORY_PROTEIN CP ON MS.ID_CATEGORY_PROTEIN=CP.ID_CATEGORY_PROTEIN AND CP.ID_PROTEIN=?");
	my $sthDelCP=$dbh->prepare("DELETE FROM CATEGORY_PROTEIN WHERE ID_PROTEIN=?");
	my $sthDelProt=$dbh->prepare("DELETE FROM PROTEIN WHERE ID_PROTEIN=?");
	foreach my $identifier (keys %preValidatedProt) {
		next if $newValidatedProt{$identifier}; # protein was confirmed
		$sthNumAna->execute($preValidatedProt{$identifier}[0]); # proteinID
		my ($exist)=$sthNumAna->fetchrow_array;
		next if $exist; # protein is not unique to current analysis
		$sthDelMS->execute($preValidatedProt{$identifier}[0]); # proteinID
		$sthDelCP->execute($preValidatedProt{$identifier}[0]); # proteinID
		$sthDelProt->execute($preValidatedProt{$identifier}[0]);
		delete $projectProt{$identifier}; # remove from project list
		$modifMasterProteins{$preValidatedProt{$identifier}[1]}=1 if $preValidatedProt{$identifier}[1]; # could be undef
	}
	$sthNumAna->finish;
	$sthDelMS->finish;
	$sthDelCP->finish;
	$sthDelProt->finish;
	print '.';

	#>Deleting unused master proteins (& species)
	&promsMod::deleteUnusedMasterProteins($dbh,\%modifMasterProteins);

	print " Done.<BR>\n";


	###########################################
	####>Transfering quantification values<####
	###########################################
	my $quantifID;
	#if ($prevQuantifID) {$quantifID=$prevQuantifID;}
	#else {
	#	($quantifID)=$dbh->selectrow_array("SELECT MAX(ID_QUANTIFICATION) FROM QUANTIFICATION");
	#	$quantifID++;
	#}
	if ($labelMethod) {
		my $numChannels=scalar keys %labelingInfo;

		####>SILAC/TMT<####
		if ($msfFile && $labelMethod=~/SILAC|TMT|iTRAQ/i && !$quantifNodeNum) {
			print "&nbsp;&nbsp;&nbsp;***WARNING: No quantification data found.***<BR>\n";
			$quantifID=0;
		}
		elsif ($msfFile && $labelMethod=~/SILAC/i && $quantifNodeNum) {
			print "&nbsp;&nbsp;&nbsp;- Extracting SILAC quantification data from $msfFile";
			my $quantiFocus='peptide';
			my ($quantifMethodID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='SILAC'");
			my ($areaParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE='ISO_AREA' AND ID_QUANTIFICATION_METHOD=$quantifMethodID");
			my ($apexRTParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE='RT_APEX' AND ID_QUANTIFICATION_METHOD=$quantifMethodID");
			$quantifName=$dbh->quote($quantifName);
			$quantifAnnot=$dbh->quote($quantifAnnot);
			#$dbh->do("INSERT IGNORE INTO QUANTIFICATION (ID_QUANTIFICATION,ID_QUANTIFICATION_METHOD,NAME,FOCUS,QUANTIF_ANNOT,STATUS,UPDATE_USER,UPDATE_DATE) VALUES ($quantifID,$quantifMethodID,$quantifName,'$quantiFocus',$quantifAnnot,-1,'$userID',NOW())");
			#$dbh->do("INSERT IGNORE INTO ANA_QUANTIFICATION (ID_ANALYSIS,ID_QUANTIFICATION) VALUES ($analysisID,$quantifID)");
			$dbh->do("INSERT INTO QUANTIFICATION (ID_QUANTIFICATION_METHOD,NAME,FOCUS,QUANTIF_ANNOT,STATUS,UPDATE_USER,UPDATE_DATE) VALUES ($quantifMethodID,$quantifName,'$quantiFocus',$quantifAnnot,-1,'$userID',NOW())");
			$quantifID=$dbh->last_insert_id(undef, undef, 'QUANTIFICATION', 'ID_QUANTIFICATION');
			$dbh->do("INSERT INTO ANA_QUANTIFICATION (ID_ANALYSIS,ID_QUANTIFICATION) VALUES ($analysisID,$quantifID)");

			$dbh->commit; # in case another user also generates quantification
			#print '.';

			###>Inserting quantification data<###

			#########################################>>>>PP 27/08/15>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
			# DATA STRUCTURE:
			# 1 MS/MS=1 SpectrumID (=1 query in PDM)
			# Heavy & Light have different SpectrumID but same QuanResultID (with different QuanChannelID)
			# 1 QuanResultID is linked to multiple SpectrumID (multiple successive MSMS matching the same pic/XIC)
			# In a QuanResultID, some QuanChannelID may NOT have a matching SpectrumID (no MSMS performed?)
			# !!!WARNING: There is no direct link between a SpectrumID and a QuanChannelID!!! Only mass(obs) comparison can allow this match!!!
			my (%spectrum2quantifSet,%quantifSetData,%LcmsFeatures);
			my $queryQDAll=($pdVersion >= 2.2)?
"SELECT MSnSpectrumInfoSpectrumID,QuanResultID,LF.ChannelID,MonoisotopicMassOverCharge,Area,ApexRT
	FROM TargetPsms TP, LcmsFeaturesTargetPsms LFTP, LcmsFeatures LF, TargetPsmsMSnSpectrumInfo TPSM
	WHERE LFTP.LcmsFeaturesWorkflowID=LF.WorkflowId
	AND LFTP.LcmsFeaturesID=LF.Id
	AND LFTP.TargetPsmsWorkflowID=TP.WorkflowID
	AND LFTP.TargetPsmsPeptideID=TP.PeptideID
	AND TPSM.TargetPsmsWorkflowID=TP.WorkflowID
	AND TPSM.TargetPsmsPeptideID=TP.PeptideID"
: ($msfFileID)?
"SELECT PIQRSS.SearchSpectrumID,PIQR.QuanResultID,QuanChannelID,PIQR.Mass,Area,RetentionTime
	FROM PrecursorIonQuanResultsSearchSpectra PIQRSS,PrecursorIonQuanResults PIQR,MassPeaks MP,SpectrumHeaders SH
	WHERE PIQRSS.QuanResultID=PIQR.QuanResultID AND PIQRSS.SearchSpectrumID=SH.SpectrumID AND SH.MassPeakID=MP.MassPeakID AND MP.FileID=$msfFileID AND ProcessingNodeNumber=$quantifNodeNum ORDER BY QuanChannelID ASC"
:
"SELECT SearchSpectrumID,PIQR.QuanResultID,QuanChannelID,PIQR.Mass,Area,RetentionTime
	FROM PrecursorIonQuanResultsSearchSpectra PIQRSS,PrecursorIonQuanResults PIQR
	WHERE PIQRSS.QuanResultID=PIQR.QuanResultID AND ProcessingNodeNumber=$quantifNodeNum";
			my $sthQDAll=$dbhLite->prepare($queryQDAll); # no label 1st
			#my $sthInsPQ=$dbh->prepare("INSERT INTO PEPTIDE_QUANTIFICATION (ID_QUANTIF_PARAMETER,ID_QUANTIFICATION,ID_PEPTIDE,QUANTIF_VALUE,TARGET_POS) VALUES (?,$quantifID,?,?,?)");
			mkdir "$promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
			my $quantifDir="$promsPath{quantification}/project_$projectID/quanti_$quantifID";
			mkdir $quantifDir;
			my %pepQuantHandle;
			foreach my $qChannelID (sort{$a<=>$b} keys %labelingInfo) {
				my $targetPos=$qChannel2TargetPos{$qChannelID} || $qChannelID;
				open ($pepQuantHandle{$targetPos},">$quantifDir/peptide_quantification_$targetPos.txt") || die $!;
				print {$pepQuantHandle{$targetPos}} "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n";
			}
			my $sthUpPep=$dbh->prepare("UPDATE PEPTIDE SET DATA=? WHERE ID_PEPTIDE=?");
			my $sthDelPPA=$dbh->prepare("DELETE FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_PEPTIDE=?");
			my $sthDelPM=$dbh->prepare("DELETE FROM PEPTIDE_MODIFICATION WHERE ID_PEPTIDE=?");
			my $sthDelGP=$dbh->prepare("DELETE FROM PEPTIDE WHERE ID_PEPTIDE=?");

			$sthQDAll->execute;
			while (my ($spectrumID,$quantifResID,$qChannelID,$mass,$area,$apexRT)=$sthQDAll->fetchrow_array) {
				next unless ($quantifResID && $qChannelID); # In PD 2_2, $quantifResID and $qChannelID could be null
				@{$LcmsFeatures{"$quantifResID:$qChannelID"}}=($spectrumID,$mass,$apexRT); # in PD 2_2, everything is retrieved with this query !
				if ($spectrum2seqRT{$spectrumID}) { # reference RT peptide
					$anaPepRefRTdata{$spectrum2seqRT{$spectrumID}}{$spectrumID}[8]=$apexRT;
					next;
				}
				$spectrum2quantifSet{$spectrumID}=$quantifResID;
				@{$quantifSetData{$quantifResID}{$qChannelID}}=($mass,$area,$apexRT); # repeated for multiple SpectrumID
			}
			$sthQDAll->finish;

			##>Looping through validated data
			my (%usedQuantifSet,%qSet2DataSource);
			my $count=0;
			foreach my $spectrumID (sort{$spectrum2PeptideIDs{$b}[1]<=>$spectrum2PeptideIDs{$a}[1]} keys %spectrum2PeptideIDs) { # sort by desc score
				$count++;
				if ($count>=1000) {
					print '.';
					$count=0;
				}
				my ($pepID,$score,$valStatus,$massVal)=@{$spectrum2PeptideIDs{$spectrumID}};
				my $quantifResID=$spectrum2quantifSet{$spectrumID};
				unless ($quantifResID) { # spectrum was not quantified
					#>Delete reported ghost peptide without quantif data (if any)
					if ($valStatus==0) { # Ghost
						$sthDelPPA->execute($pepID);
						$sthDelPM->execute($pepID);
						$sthDelGP->execute($pepID);
					}
					next;
				}
				#>Checking if peptide sequence is quantifiable
				if ($peptideInfo{$pepID}[0] !~ /[$labelResStrg]/) { # sequence does not contain labelable res (eg. protein C-term peptide)
#print "$peptideInfo{$pepID}[0]<BR>\n";
					#>Delete reported ghost peptide (if any)
					if ($valStatus==0) { # Ghost
						$sthDelPPA->execute($pepID);
						$sthDelPM->execute($pepID);
						$sthDelGP->execute($pepID);
					}
					foreach my $qChannelID (keys %{$quantifSetData{$quantifResID}}) { # set this QuanResultID as used for all QuanChannelIDs to prevent re-use during 2nd scan
						$usedQuantifSet{$quantifResID}{$qChannelID}=$pepID;
					}
					next;
				}
				#>Matching labeling channels based on mass
				foreach my $qChannelID (keys %{$quantifSetData{$quantifResID}}) {
					my ($massMsf,$area,$apexRT)=@{$quantifSetData{$quantifResID}{$qChannelID}};
					if (abs($massVal-$massMsf) < 0.5) { # match! TODO: Check seq coherence (rare abherent cases detected) & Phospho pos coherence if PhosphoRS
						if ($usedQuantifSet{$quantifResID} && $usedQuantifSet{$quantifResID}{$qChannelID}) { # quantifSet & qChannel already used by other (better scoring) peptide
							#>Delete reported ghost peptide (if any)
							if ($valStatus==0) { # Ghost
								$sthDelPPA->execute($pepID);
								$sthDelPM->execute($pepID);
								$sthDelGP->execute($pepID);
							}
							last;
						}
						$usedQuantifSet{$quantifResID}{$qChannelID}=$pepID;
						#>Store peptide quantif data & update PEPTIDE.DATA
						#$sthInsPQ->execute($areaParamID,$pepID,$area,$qChannelID); # Area
						#$sthInsPQ->execute($apexRTParamID,$pepID,$apexRT,$qChannelID); # RetentionTime
						my $targetPos=$qChannel2TargetPos{$qChannelID} || $qChannelID;
						print {$pepQuantHandle{$targetPos}} "$areaParamID\t$pepID\t$area\n$apexRTParamID\t$pepID\t$apexRT\n";
						my $qSetStrg=($peptideExtraData{$pepID})? $peptideExtraData{$pepID}.'##QSET='.$quantifResID : 'QSET='.$quantifResID;
						$sthUpPep->execute($qSetStrg,$pepID); # records QuanResultID
						$qSet2DataSource{$quantifResID}=$mergedFiles{$spectrum2fileID{$spectrumID}}[1] if $numMergedFiles > 1; # merged file rank (recorded in case missing qChannel with no SpectrumID at all-> 2nd Scan)
						last;
					}
				}
			}
			$sthDelPPA->finish;
			$sthDelPM->finish;
			$sthDelGP->finish;
			$sthUpPep->finish;

			##>Scanning for missing Channel peptide
			print '/';
			#my $sthVP=($pdVersion >= 2.2)?
			#$dbhLite->prepare("SELECT SpectrumID,LF.MonoisotopicMassOverCharge,LF.ApexRT
			#							    FROM MSnSpectrumInfo MS, TargetPsmsMSnSpectrumInfo TPMS, TargetPsms TP,  LcmsFeaturesTargetPsms LFTP, LcmsFeatures LF
			#							    WHERE TPMS.MSnSpectrumInfoWorkflowID=MS.WorkflowId
			#							    AND TPMS.MSnSpectrumInfoSpectrumID=MS.SpectrumID
			#							    AND TPMS.TargetPsmsWorkflowID=TP.WorkflowID
			#							    AND TPMS.TargetPsmsPeptideID=TP.PeptideID
			#							    AND LFTP.TargetPsmsWorkflowID=TP.WorkflowID
			#							    AND LFTP.TargetPsmsPeptideID=TP.PeptideID
			#							    AND LFTP.LcmsFeaturesWorkflowID=LF.WorkflowId
			#							    AND LFTP.LcmsFeaturesID=LF.Id
			#							    AND QuanResultID=? AND ChannelID=? AND MS.Charge=LF.ChargeState
			#							    ORDER BY ABS(MS.RetentionTime-ApexRT) ASC LIMIT 0,1"):
			my $sthVP;
			if($pdVersion < 2.2){
				$sthVP=$dbhLite->prepare("SELECT SpectrumID,MP.Mass,SH.RetentionTime
											FROM PrecursorIonQuanResultsSearchSpectra PIQRSS,PrecursorIonQuanResults PIQR,SpectrumHeaders SH,MassPeaks MP
											WHERE PIQRSS.QuanResultID=PIQR.QuanResultID AND PIQRSS.SearchSpectrumID=SH.SpectrumID AND SH.MassPeakID=MP.MassPeakID
											AND PIQRSS.ProcessingNodeNumber=$quantifNodeNum AND PIQR.QuanResultID=? AND QuanChannelID=? AND PIQR.Charge=MP.Charge AND ABS(PIQR.Mass-MP.Mass) <= 0.5
											ORDER BY ABS(PIQR.RetentionTime-SH.RetentionTime) ASC LIMIT 0,1");
			}
			my $sthInsVP=$dbh->prepare("INSERT INTO PEPTIDE (ID_ANALYSIS,VALID_STATUS,PEP_SEQ,PEP_LENGTH,MISS_CUT,CHARGE,MR_OBS,ELUTION_TIME,DATA) VALUE ($analysisID,0,?,?,?,?,?,?,?)");
			my $sthInsVPA=$dbh->prepare("INSERT INTO PEPTIDE_PROTEIN_ATTRIB (ID_ANALYSIS,ID_PEPTIDE,ID_PROTEIN,PEP_BEG,PEP_END,FLANKING_AA,IS_SPECIFIC) SELECT $analysisID,?,ID_PROTEIN,-ABS(PEP_BEG),-ABS(PEP_END),FLANKING_AA,IS_SPECIFIC FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_PEPTIDE=?");
			my $sthInsVPM=$dbh->prepare("INSERT INTO PEPTIDE_MODIFICATION (ID_PEPTIDE,ID_MODIFICATION,MODIF_TYPE,POS_STRING,REF_POS_STRING) VALUES (?,?,'V',?,?)");

			foreach my $quantifResID (keys %usedQuantifSet) {
				next if scalar keys %{$usedQuantifSet{$quantifResID}}==$numChannels;
				my $refPepID=$usedQuantifSet{$quantifResID}{(sort{$a<=>$b} keys %{$usedQuantifSet{$quantifResID}})[0]}; # qChannelID
				foreach my $qChannelID (keys %{$quantifSetData{$quantifResID}}) {
					next if $usedQuantifSet{$quantifResID}{$qChannelID}; # already used
					$count++;
					if ($count>=1000) {
						print '.';
						$count=0;
					}
					my ($spectrumID,$massObs,$ms2rt);
					if ($pdVersion < 2.2) {
						$sthVP->execute($quantifResID,$qChannelID);
						($spectrumID,$massObs,$ms2rt)=$sthVP->fetchrow_array; # can be null
					}
					else{
						($spectrumID,$massObs,$ms2rt)=$LcmsFeatures{"$quantifResID:$qChannelID"};
					}

					#Insert peptide: ***2 classes of virtual peptides depending on whether a spectrumID is matched or not***
					my $dataStrg;
					if ($spectrumID) {
						$dataStrg='SOURCE_RANK='.$mergedFiles{$spectrum2fileID{$spectrumID}}[1].'##' if $numMergedFiles > 1; # merged file rank
						$dataStrg.="EXT_SPECTRUMID=$spectrumID##QSET=$quantifResID";
					}
					else { # no matching SpectrumID
						$dataStrg='SOURCE_RANK='.$qSet2DataSource{$quantifResID}.'##' if $numMergedFiles > 1; # merged file rank
						$dataStrg.="QSET=$quantifResID";
						$massObs=$quantifSetData{$quantifResID}{$qChannelID}[0]; # use PIQR.Mass
					}
					$sthInsVP->execute(@{$peptideInfo{$refPepID}}[0..3],$massObs,$ms2rt,$dataStrg); #@{...}[0..3]=$seq,length($seq),$missCut,$charge
					my $peptideID=$dbh->last_insert_id(undef,undef,'PEPTIDE','ID_PEPTIDE');
					#Insert unlabeled varmods (if any)
					foreach my $refMod (@{$peptideInfo{$refPepID}[4]}) {
						next if $labelModifID{$refMod->[0]};
						$sthInsVPM->execute($peptideID,@{$refMod});
					}
					#Insert labeled modifs associated with QuanChannelID (if any)
					foreach my $modID (keys %{$qChannelID2ModifID{$qChannelID}}) {
						my $resStrg=join('|',@{$qChannelID2ModifID{$qChannelID}{$modID}}); # SILAC => normal res K/R (not =-+*)
						my @mPos;
						while ($peptideInfo{$refPepID}[0]=~/$resStrg/g) {
							push @mPos,pos($peptideInfo{$refPepID}[0]);
						}
						$sthInsVPM->execute($peptideID,$modID,join('.',@mPos),undef) if $mPos[0];
					}
					#Insert protein matches
					$sthInsVPA->execute($peptideID,$refPepID);
					#Store peptide quantif data
					my ($massMsf,$area,$apexRT)=@{$quantifSetData{$quantifResID}{$qChannelID}};
					#$sthInsPQ->execute($areaParamID,$peptideID,$area,$qChannelID); # Area
					#$sthInsPQ->execute($apexRTParamID,$peptideID,$apexRT,$qChannelID); # RetentionTime
					my $targetPos=$qChannel2TargetPos{$qChannelID} || $qChannelID;
					print {$pepQuantHandle{$targetPos}} "$areaParamID\t$peptideID\t$area\n$apexRTParamID\t$peptideID\t$apexRT\n";
					$usedQuantifSet{$quantifResID}{$qChannelID}=$peptideID;
				}
			}
			#$sthInsPQ->finish;
			foreach my $handle (values %pepQuantHandle) {close $handle;}
			$sthVP->finish if $pdVersion < 2.2;
			$sthInsVP->finish;
			$sthInsVPA->finish;
			$sthInsVPM->finish;

			$dbhLite->disconnect;

			print " Done<BR>\n";
		}
		elsif ($msfFile && $labelMethod=~/TMT|iTRAQ/i && $quantifNodeNum) {
			print "&nbsp;&nbsp;&nbsp;- Extracting $labelMethod quantification data from $msfFile";
			my $quantiFocus='peptide';
			my ($qMethod)=($labelMethod=~/TMT/)?'TMT':'ITRAQ';
			my ($quantifMethodID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='$qMethod'");
			my ($heightParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE='REP_INTENSITY' AND ID_QUANTIFICATION_METHOD=$quantifMethodID");
			my ($massParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE='REP_MASS' AND ID_QUANTIFICATION_METHOD=$quantifMethodID");
			$quantifName=$dbh->quote($quantifName);
			$quantifAnnot=$dbh->quote($quantifAnnot);
			$dbh->do("INSERT INTO QUANTIFICATION (ID_QUANTIFICATION_METHOD,NAME,FOCUS,QUANTIF_ANNOT,STATUS,UPDATE_USER,UPDATE_DATE) VALUES ($quantifMethodID,$quantifName,'$quantiFocus',$quantifAnnot,-1,'$userID',NOW())");
			$quantifID=$dbh->last_insert_id(undef, undef, 'QUANTIFICATION', 'ID_QUANTIFICATION');
			$dbh->do("INSERT INTO ANA_QUANTIFICATION (ID_ANALYSIS,ID_QUANTIFICATION) VALUES ($analysisID,$quantifID)");
			$dbh->commit; # in case another user also generates quantification
			#print '.';

			###>Inserting quantification data<###

			#########################################>>>>GA 30/03/17>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
			# DATA STRUCTURE:
			# 1 MS/MS=1 SpectrumID (=1 query in PDM)
			# Same SpectrumID for all TMT-reporter (with different QuanChannelID and Mass)
			# 1 SpectrumID in ReporterIonQuanResults table is linked to 1 SearchSpectrumID which is equal to SpectrumHeaders.SpectrumID
			my (%spectrum2quantifSet,%quantifSetData);
			my ($queryQDAll)=($pdVersion>=2.2)?
			"SELECT QuanChannelID,Mass,Height,QuanResultID,MSnSpectrumInfoSpectrumID
FROM ReporterQuanResults RQS, ReporterQuanResultsMSnSpectrumInfo RQRMS
WHERE RQS.WorkflowID=RQRMS.ReporterQuanResultsWorkflowID
AND RQS.QuanResultID=RQRMS.ReporterQuanResultsQuanResultID
AND RQS.QuanChannelID=RQRMS.ReporterQuanResultsQuanChannelID
AND Height>0 AND Mass>0"
			:"SELECT RQR.QuanChannelID,RQR.Mass,RQR.Height,RQR.SpectrumID,RQRSS.SearchSpectrumID
FROM ReporterIonQuanResults RQR, ReporterIonQuanResultsSearchSpectra RQRSS
WHERE RQR.SpectrumID=RQRSS.SpectrumID
AND RQR.Height>0 AND RQR.Mass>0";
			my $sthQDAll=$dbhLite->prepare($queryQDAll);
			#my $sthInsPQ=$dbh->prepare("INSERT INTO PEPTIDE_QUANTIFICATION (ID_QUANTIF_PARAMETER,ID_QUANTIFICATION,ID_PEPTIDE,QUANTIF_VALUE,TARGET_POS) VALUES (?,$quantifID,?,?,?)");
			mkdir "$promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
			my $quantifDir="$promsPath{quantification}/project_$projectID/quanti_$quantifID";
			mkdir $quantifDir;
			my %pepQuantHandle;
			foreach my $reporter (sort{$reporterIonList{$a}[0]<=>$reporterIonList{$b}[0]} keys %reporterIonList) {
				my $qChannelID=$reporterIonList{$reporter}[0];
				my $targetPos=$qChannel2TargetPos{$qChannelID} || $qChannelID;
				open($pepQuantHandle{$targetPos},">$quantifDir/peptide_quantification_$targetPos.txt");
				print {$pepQuantHandle{$targetPos}} "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n";
			}
			#my $sthUpPep=$dbh->prepare("UPDATE PEPTIDE SET DATA=? WHERE ID_PEPTIDE=?");

			$sthQDAll->execute;
			my $count=0;
			while (my ($qChannelID,$mass,$height,$quantifResID,$spectrumID)=$sthQDAll->fetchrow_array) {
				next unless $spectrum2PeptideIDs{$spectrumID};
				my ($pepID,$score,$valStatus,$massVal)=@{$spectrum2PeptideIDs{$spectrumID}};
				next unless $pepID;
				#$sthInsPQ->execute($heightParamID,$pepID,$height,$qChannelID);
				#$sthInsPQ->execute($massParamID,$pepID,$mass,$qChannelID);
				my $targetPos=$qChannel2TargetPos{$qChannelID} || $qChannelID;
				print {$pepQuantHandle{$targetPos}} "$heightParamID\t$pepID\t$height\n$massParamID\t$pepID\t$mass\n";
				if ($count>=1000) {
					print '.';
					$count=0;
				}
				$count++;
			}
			$sthQDAll->finish;
			$dbhLite->disconnect;
			foreach my $handle (values %pepQuantHandle) {close $handle;}

			print " Done<BR>\n";
		}

		####>iTRAQ (non-PD)<####
		elsif ($labelMethod=~/iTRAQ/i) {
			print "&nbsp;&nbsp;&nbsp;- Extracting iTRAQ quantification data from $dataFileName...";
			my $quantiFocus='peptide';
			my ($quantifMethodID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='ITRAQ'");
			my ($areaParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE='REP_AREA' AND ID_QUANTIFICATION_METHOD=$quantifMethodID");
			#my ($quantifID)=$dbh->selectrow_array("SELECT MAX(ID_QUANTIFICATION) FROM QUANTIFICATION");
			#$quantifID++;
			#my $sthInsPQ=$dbh->prepare("INSERT INTO PEPTIDE_QUANTIFICATION (ID_QUANTIF_PARAMETER,ID_QUANTIFICATION,ID_PEPTIDE,QUANTIF_VALUE,TARGET_POS) VALUES ($areaParamID,$quantifID,?,?,?)");
			mkdir "$promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
			my $quantifDir="$promsPath{quantification}/project_$projectID/quanti_$quantifID";
			mkdir $quantifDir;
			my (%pepQuantHandle);

			###>MASCOT<###
			if ($fileFormat eq 'MASCOT.DAT') {
				my $currentQuery=0;
				my $section='';
				my $quantificationXML='';
				my $reporterMaxMass=0;
				my $qCount=0;
				open (DAT,$dataFile) || die "Cannot open $dataFile";
				while (<DAT>) {
					if ($_=~/gc0p4Jq0M2Yt08jU534c0p/) {
						if ($section eq 'quantitation') { # end of quantitation section -> process XML string
							##>Extracting & storing quantification parameters
							my $xml = new XML::Simple();
							my $xmlData = $xml->XMLin($quantificationXML);
							my $quantifName=$xmlData->{method}{name};
							my $quantifAnnot='LABEL=ITRAQ::SOFTWARE=MAS';
							my $componentXML=$xmlData->{method}{component};
							my $reporterPos=0;
							foreach my $reporter (sort keys %{$componentXML}) {
								$reporterPos++;
								@{$reporterIonList{$reporter}}=($reporterPos,$componentXML->{$reporter}{moverz}{monoisotopic});
								$quantifAnnot.="::$reporterPos;$reporter;".$componentXML->{$reporter}{moverz}{monoisotopic};
								$reporterMaxMass=$componentXML->{$reporter}{moverz}{monoisotopic} if $reporterMaxMass < $componentXML->{$reporter}{moverz}{monoisotopic};
open($pepQuantHandle{$reporterPos},">$quantifDir/peptide_quantification_$reporterPos.txt");
print {$pepQuantHandle{$reporterPos}} "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n";
							}
							$reporterMaxMass+=0.2; # extending limit because of potential mass error
							##>Storing Quantification
							#($quantifID)=$dbh->selectrow_array("SELECT MAX(ID_QUANTIFICATION) FROM QUANTIFICATION");
							#$quantifID++;
							$quantifName=$dbh->quote($quantifName);
							$quantifAnnot=$dbh->quote($quantifAnnot);
							$dbh->do("INSERT INTO QUANTIFICATION (ID_QUANTIFICATION_METHOD,NAME,FOCUS,QUANTIF_ANNOT,UPDATE_USER,UPDATE_DATE) VALUES ($quantifMethodID,$quantifName,'$quantiFocus',$quantifAnnot,'$userID',NOW())");
							$quantifID=$dbh->last_insert_id(undef, undef, 'QUANTIFICATION', 'ID_QUANTIFICATION');
							$dbh->do("INSERT INTO ANA_QUANTIFICATION (ID_ANALYSIS,ID_QUANTIFICATION) VALUES ($analysisID,$quantifID)");
							$dbh->commit; # in case another user also generates quantification
							$quantificationXML=undef; # free memory
							print '.';
						}
						$section='';
						$currentQuery=0;
						next;
					}
					##>Quantification settings
					if ($section eq 'quantitation') {
						$quantificationXML.=$_ if $_=~/\w/;
						next;
					}
					##>Query sections
					if ($section eq '' && $_=~/name="(.+)"/) {
						$section=$1;
						if ($section =~ /query(\d+)/) {
							$currentQuery=$1;
							$section='query';
						}
						next;
					}
					##>Query spectrum
					if ($currentQuery && $validQueryData{$currentQuery} && $_=~/Ions1=(.*)/) {
						$qCount++;
						if ($qCount>=50) {
							print '.';
							$qCount=0;
						}
						my $spectrumStrg=$1;
						my %mass2intensities;
						foreach my $ionData (split(',',$spectrumStrg)) {
							$mass2intensities{$1}=$2 if ($ionData=~/(.+):(.+)/ && $1 < $reporterMaxMass);
						}
						foreach my $mass (sort{$a<=>$b} keys %mass2intensities) {
							foreach my $reporter (keys %reporterIonList) {
								if (abs($mass-$reporterIonList{$reporter}[1]) < 0.2) { # match!
									foreach my $pepID (@{$validQueryData{$currentQuery}{'pepID'}}) {
										#$sthInsPQ->execute($pepID,$mass2intensities{$mass},$reporterIonList{$reporter}[0]);
										print {$pepQuantHandle{$reporterIonList{$reporter}[0]}} "$areaParamID\t$pepID\t$mass2intensities{$mass}\n";
									}
									last;
								}
							}
						}
					}
					last if $section eq 'index';
				}
				close DAT;
			}

			###>PARAGON<###
			elsif ($fileFormat eq 'PARAGON.XML') {
				my $maxSpectra=scalar keys %spectrum2PeptideIDs;
				open (XML, $dataFile) || die "Cannot open $dataFile";
				my $section='';
				my ($reporterPos,$spectrumID,$dataSource,%sourceList);
				my $matchedSpectra=0;
				while (my $line=<XML>) {
					if ($line =~ /<QUANT_TYPE xml:id="$labelMethod"/){
						$section='LABELING';
						$reporterPos=0; # blocs are repeated
					}
					elsif ($section eq 'LABELING' && $line =~ /<\/QUANT_TYPE>/) {$section='';}
					elsif ($line =~ /<PEAKLIST .+ originalfilename="([^"]+)/) {
						my $fullSource=$1;
						($dataSource)=(split(/[\\\/]/,$fullSource))[-1];
						$sourceList{$dataSource}=1;
					}
					elsif ($line =~ /<SPECTRUM .+ xml:id="([^"]+)/ && $spectrum2PeptideIDs{$1}) {
						$spectrumID=$1;
						$section='SPECTRUM';
						$matchedSpectra++;
#print "<BR>>$spectrumID:<BR>\n";
					}
					elsif ($section eq 'SPECTRUM' && $line =~ /<ITRAQPEAKS /) {$section='ITRAQPEAKS';}
					elsif ($line =~ /<\/ITRAQPEAKS>/) {
						last if $matchedSpectra==$maxSpectra; # all validated spectra are matched: no need to continue scanning
						$section='';
					}
					#>iTRAQ parameters
					elsif ($section eq 'LABELING' && $line =~ /<QUANT_LABEL .+ display="([^"]+)/) {
						my $reporter=$1;
						$reporterPos++;
						@{$reporterIonList{$reporter}}=($reporterPos,''); # $reporterIonList{$reporter}[1] will not be updated if no corresponding quant data at all (bad quantif)
open($pepQuantHandle{$reporterPos},">$quantifDir/peptide_quantification_$reporterPos.txt");
print {$pepQuantHandle{$reporterPos}} "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n";
					}
					#>iTRAQ data
					elsif ($section eq 'ITRAQPEAKS') {
						chomp($line);
						$line =~ s/<\!\[CDATA\[//;
						$line =~ s/\]\]>//;
						my ($centroid,$area,$areaError)=split(/\t/,$line);
						my $reporter=int($centroid+0.5)*1; # 113.1067 -> 113
						$reporterIonList{$reporter}[1]=$centroid unless $reporterIonList{$reporter}[1];
						$spectrum2PeptideIDs{$spectrumID}{'QUANTIF'}{$reporter}=$area;
						$spectrum2PeptideIDs{$spectrumID}{'SOURCE'}=$dataSource;
#print "&nbsp;($spectrumID) $reporterIonList{$reporter}[0] $reporter: $area<BR>\n";
					}
				}
				close XML;

				##>Storing Quantification
				my $quantifAnnot='LABEL=ITRAQ::SOFTWARE=PAR';
				foreach my $reporter (sort{$reporterIonList{$a}[0]<=>$reporterIonList{$b}[0]} keys %reporterIonList) {
					$quantifAnnot.="::$reporterIonList{$reporter}[0];$reporter;$reporterIonList{$reporter}[1]";
				}
				my $quantifName=$dbh->quote($labelMethod);
				$quantifAnnot=$dbh->quote($quantifAnnot);
				$dbh->do("INSERT INTO QUANTIFICATION (ID_QUANTIFICATION_METHOD,NAME,FOCUS,QUANTIF_ANNOT,UPDATE_USER,UPDATE_DATE) VALUES ($quantifMethodID,$quantifName,'$quantiFocus',$quantifAnnot,'$userID',NOW())");
				$quantifID=$dbh->last_insert_id(undef, undef, 'QUANTIFICATION', 'ID_QUANTIFICATION');
				$dbh->do("INSERT INTO ANA_QUANTIFICATION (ID_ANALYSIS,ID_QUANTIFICATION) VALUES ($analysisID,$quantifID)");
				$dbh->commit; # in case another user also generates quantification
				print '.';

				##>Storing peptide quantification
				my $numSources=scalar keys %sourceList;
				foreach my $specID (sort{$spectrum2PeptideIDs{$a}{'PEP_ID'}[0]<=>$spectrum2PeptideIDs{$b}{'PEP_ID'}[0]} keys %spectrum2PeptideIDs) {
					foreach my $reporter (sort{$reporterIonList{$a}[0]<=>$reporterIonList{$b}[0]} keys %{$spectrum2PeptideIDs{$specID}{'QUANTIF'}}) {
						foreach my $pepID (@{$spectrum2PeptideIDs{$specID}{'PEP_ID'}}) {
							#my $usedSource=($numSources > 1)? $spectrum2PeptideIDs{$specID}{'SOURCE'} : undef;
							#$sthInsPQ->execute($pepID,$spectrum2PeptideIDs{$specID}{'QUANTIF'}{$reporter},$reporterIonList{$reporter}[0]); #,$usedSource
							print {$pepQuantHandle{$reporterIonList{$reporter}[0]}} "$areaParamID\t$pepID\t$spectrum2PeptideIDs{$specID}{'QUANTIF'}{$reporter}\n";
						}
					}
				}
			}
			#$sthInsPQ->finish;
			foreach my $handle (values %pepQuantHandle) {close $handle;}
			print " Done<BR>\n";
		}
		else {
			$quantifID=0;
			print "&nbsp;&nbsp;&nbsp;***WARNING: Unhandled quantification method detected: $labelMethod. Skipping quantification data import.***<BR>\n";
		}
	}
	else {
		###>LABEL-FREE with Proteome Discoverer<###
		if ($msfFile && $quantifNodeNum) {
			print "&nbsp;&nbsp;&nbsp;- Extracting Label-free quantification data from $msfFile:";
			my $quantiFocus='peptide';
			my ($quantifMethodID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='XIC'");
			my ($areaParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE='XIC_AREA' AND ID_QUANTIFICATION_METHOD=$quantifMethodID");
			my ($apexRTParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE='RT_APEX' AND ID_QUANTIFICATION_METHOD=$quantifMethodID");
			#my ($quantifID)=$dbh->selectrow_array("SELECT MAX(ID_QUANTIFICATION) FROM QUANTIFICATION");
			#$quantifID++;
			$quantifName=$dbh->quote($quantifName);
			$quantifAnnot=$dbh->quote($quantifAnnot);
			$dbh->do("INSERT IGNORE INTO QUANTIFICATION (ID_QUANTIFICATION_METHOD,NAME,FOCUS,QUANTIF_ANNOT,STATUS,UPDATE_USER,UPDATE_DATE) VALUES ($quantifMethodID,$quantifName,'$quantiFocus',$quantifAnnot,-1,'$userID',NOW())");
			$quantifID=$dbh->last_insert_id(undef, undef, 'QUANTIFICATION', 'ID_QUANTIFICATION');
			$dbh->do("INSERT IGNORE INTO ANA_QUANTIFICATION (ID_ANALYSIS,ID_QUANTIFICATION) VALUES ($analysisID,$quantifID)");
			$dbh->commit; # in case another user also generates quantification
			print '.';

			if ($pdVersion >= 2.2){
				###>Extracting quantification data<###
				my $sthArea=$dbhLite->prepare("SELECT MSnSpectrumInfoSpectrumID,Area,ApexRT
FROM LcmsFeatures LF, LcmsFeaturesTargetPsms LFTP, TargetPsms TP, TargetPsmsMSnSpectrumInfo TPMSI
WHERE LF.WorkflowID=LFTP.LcmsFeaturesWorkFlowId
AND LF.Id=LFTP.LcmsFeaturesId
AND LFTP.TargetPsmsWorkflowID=TP.WorkFlowID
AND LFTP.TargetPsmsPeptideID=TP.PeptideID
AND TP.WorkFlowID=TPMSI.TargetPsmsWorkflowID
AND TP.PeptideID=TPMSI.TargetPsmsPeptideID");
				$sthArea->execute;
				#my %sourceList;
				while (my ($spectrumID,$area,$retTime)=$sthArea->fetchrow_array) {
					if ($spectrum2PeptideIDs{$spectrumID}) { # in case Analysis is restricted to a specific msfFileID
						$spectrum2PeptideIDs{$spectrumID}{'QUANTIF'}=$area;
						$spectrum2PeptideIDs{$spectrumID}{'APEX_RT'}=$retTime;
						if ($spectrum2seqRT{$spectrumID}) { # iRT peptide
								$anaPepRefRTdata{$spectrum2seqRT{$spectrumID}}{$spectrumID}[8]=$retTime;
								next;
						}
						$spectrum2PeptideIDs{$spectrumID}{'APEX_RT'}=$retTime;
					}
				}
				$sthArea->finish;
				print '.';
			}
			else{
				my ($fieldID)=$dbhLite->selectrow_array("SELECT FieldID FROM CustomDataFields WHERE DisplayName='Area' AND SourceNodeNumber=$quantifNodeNum");
				my $sthArea=$dbhLite->prepare("SELECT CustomDataSpectra.SpectrumID,FieldValue,FileID FROM CustomDataSpectra,SpectrumHeaders,MassPeaks WHERE CustomDataSpectra.SpectrumID=SpectrumHeaders.SpectrumID AND SpectrumHeaders.MassPeakID=MassPeaks.MassPeakID AND FieldID=$fieldID");
				$sthArea->execute;
				#my %sourceList;
				while (my ($spectrumID,$area,$sourceID)=$sthArea->fetchrow_array) {
					if ($spectrum2PeptideIDs{$spectrumID}) { # in case Analysis is restricted to a specific msfFileID
						$spectrum2PeptideIDs{$spectrumID}{'QUANTIF'}=$area;
						#$spectrum2PeptideIDs{$spectrumID}{'SOURCE'}=$sourceID;
						#$sourceList{$sourceID}=1;
					}
				}
				$sthArea->finish;
				print '.';

				###>Extracting Apex retention times from MSF<###
				my $sthAllRT=$dbhLite->prepare("SELECT SearchSpectrumID,MAX(Area),RT FROM PrecursorIonAreaSearchSpectra PIASS,EventAreaAnnotations EAA,Events E WHERE PIASS.QuanResultID=EAA.QuanResultID AND EAA.EventID=E.EventID GROUP BY SearchSpectrumID");
				$sthAllRT->execute;
				while (my($spectrumID,$maxArea,$retTime)=$sthAllRT->fetchrow_array) {
					if ($spectrum2seqRT{$spectrumID}) { # iRT peptide
						$anaPepRefRTdata{$spectrum2seqRT{$spectrumID}}{$spectrumID}[8]=$retTime;
						next;
					}
					$spectrum2PeptideIDs{$spectrumID}{'APEX_RT'}=$retTime if $spectrum2PeptideIDs{$spectrumID}; # do not record if peptide not validated
				}
				$sthAllRT->finish;
				print '.';
			}

			###>Storing peptide quantification/RT<###
			#my $sthInsPQ=$dbh->prepare("INSERT INTO PEPTIDE_QUANTIFICATION (ID_QUANTIFICATION,ID_QUANTIF_PARAMETER,ID_PEPTIDE,QUANTIF_VALUE) VALUES ($quantifID,?,?,?)");
			mkdir "$promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
			my $quantifDir="$promsPath{quantification}/project_$projectID/quanti_$quantifID";
			mkdir $quantifDir;
			open(PEP_QUANTIF,">$quantifDir/peptide_quantification.txt"); # no TARGET_POS for LABEL-FREE
			print PEP_QUANTIF "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n";
			#my $numSources=scalar keys %sourceList;
			foreach my $specID (sort{$spectrum2PeptideIDs{$a}{'PEP_ID'}[0]<=>$spectrum2PeptideIDs{$b}{'PEP_ID'}[0]} keys %spectrum2PeptideIDs) {
				next unless $spectrum2PeptideIDs{$specID}{'QUANTIF'};
				foreach my $pepID (@{$spectrum2PeptideIDs{$specID}{'PEP_ID'}}) {
					#my $usedSource=($numSources > 1)? $mergedFiles{$spectrum2PeptideIDs{$specID}{'SOURCE'}} : undef;
					#my $usedSource=($numSources > 1)? $mergedFiles{$spectrum2PeptideIDs{$specID}{'SOURCE'}}[1] : undef;
					#$sthInsPQ->execute($pepID,$spectrum2PeptideIDs{$specID}{'QUANTIF'},$usedSource);
					#$sthInsPQ->execute($areaParamID,$pepID,$spectrum2PeptideIDs{$specID}{'QUANTIF'});
					print PEP_QUANTIF "$areaParamID\t$pepID\t$spectrum2PeptideIDs{$specID}{QUANTIF}\n";
					#$sthInsPQ->execute($apexRTParamID,$pepID,$spectrum2PeptideIDs{$specID}{'APEX_RT'}) if $spectrum2PeptideIDs{$specID}{'APEX_RT'};
					print PEP_QUANTIF "$apexRTParamID\t$pepID\t$spectrum2PeptideIDs{$specID}{APEX_RT}\n" if $spectrum2PeptideIDs{$specID}{'APEX_RT'};
				}
			}
			#$sthInsPQ->finish;
			close PEP_QUANTIF;
			$dbhLite->disconnect;

			print " Done<BR>\n";
		}
		else {
			$quantifID=0;
		}
	}
	$dbh->do("UPDATE QUANTIFICATION SET STATUS=1 WHERE ID_QUANTIFICATION=$quantifID") if $quantifID; # sets status as finished


	####>No need because of auto-increment
	#####>Removing protecting IDs
	#$dbh->do("DELETE FROM PROTEIN WHERE ID_PROTEIN=$maxProteinID") || $dbh->errstr;


	##################################################
	####>Recording & processing reference RT data<####
	##################################################
	if ($refRT_ID) {
		my (%rtVectors,%numPepUsed,%quantData); # PD
		my $anaNumpep=0;
		my $anaDataXml='<PEPTIDE_DATA version="1.0.2">';
		foreach my $pepSeq (sort{$referenceRTdata{$refRT_ID}{'PEPTIDES'}{$a}[1]<=>$referenceRTdata{$refRT_ID}{'PEPTIDES'}{$b}[1]} keys %anaPepRefRTdata) {
			my $pepID=$referenceRTdata{$refRT_ID}{'PEPTIDES'}{$pepSeq}[0];
			my %usedSource;
			foreach my $spectrumID (sort{$anaPepRefRTdata{$pepSeq}{$b}[6]<=>$anaPepRefRTdata{$pepSeq}{$a}[6]} keys %{$anaPepRefRTdata{$pepSeq}}) { # pepScore
				my ($qNum,$charge,$massExp,$massObs,$massCalc,$deltaMass,$pepScore,$elutionData,$rtApex)=@{$anaPepRefRTdata{$pepSeq}{$spectrumID}};
				my $spectrumStrg='';
				my $srcRk=0;
				my $srcStrg='';
				my $rtApexStrg='';
				if ($msfFile) { # PD
					#next if ($quantifNodeNum && !$rtApex);
					$spectrumStrg="";
					my $srcFileID;
					if ($numMergedFiles > 1) {
						$srcFileID=$spectrum2fileID{$spectrumID};
						$srcRk=$mergedFiles{$srcFileID}[1];
						$srcStrg=" sourceRank=\"$srcRk\"";
					}
					else {
						$srcFileID=(keys %mergedFiles)[0];
					}
					unless ($usedSource{$srcRk}) {
						$spectrumStrg=" extSpectrum=\"$spectrumID\"";
						if ($rtApex) { # => $quantifNodeNum
							push @{$rtVectors{$srcFileID}{'ANA'}},$rtApex;
							push @{$rtVectors{$srcFileID}{'REF'}},$referenceRTdata{$refRT_ID}{'PEPTIDES'}{$pepSeq}[1]; # ref RT
							$numPepUsed{$srcRk}++; # rank
							push @{$quantData{$srcRk}},[$pepID,$rtApex];
						}
					}
				}
				unless ($usedSource{$srcRk}) {
					$anaDataXml.=qq|\n<PEPTIDE pepId="$pepID"$srcStrg$spectrumStrg queryNumber="$qNum" massExp="$massExp" massObs="$massObs" massCalc="$massCalc" deltaMass="$deltaMass" charge="$charge" score="$pepScore" ms2RT="$elutionData"/>|;
					$anaNumpep++;
					$usedSource{$srcRk}=1;
				}
			}
		}
		$anaDataXml.="\n</PEPTIDE_DATA>";

		###>Insert in ANALYSIS_REFRT table<###
		$sthInsART->execute($refRT_ID,$analysisID,$anaNumpep,$anaDataXml);

		###>Compute RT correlation<###
		if ($msfFile && $quantifNodeNum && $anaNumpep >=2) { # PD
			print "&nbsp;&nbsp;&nbsp;- Computing retention time correction factors...";
			my %regressResults=&runRegressionR("$promsPath{valid}/ana_$analysisID",$refRT_ID,$quantifID,\%rtVectors,\%mergedFiles);

			###>Insert in QUANTIF_REFRT table<###
			foreach my $srcRk (sort{$a<=>$b} keys %regressResults) {
				my $usedRank=($numMergedFiles==1)? 0 : $srcRk;
				my $quantDataXml='<PEPTIDE_DATA version="1.0.1">';
				foreach my $refPep (@{$quantData{$usedRank}}) {
					$quantDataXml.=qq |<PEPTIDE pepId="$refPep->[0]" apexRT="$refPep->[1]"/>|;
				}
				$quantDataXml.='</PEPTIDE_DATA>';
				$sthInsQRT->execute($refRT_ID,$quantifID,$usedRank,$numPepUsed{$usedRank},$regressResults{$srcRk}{'SLOPE'},$regressResults{$srcRk}{'Y_INTERCEPT'},$regressResults{$srcRk}{'R2'},$quantDataXml);
#print "<BR>[SRC#$usedRank: S=$regressResults{$srcRk}{SLOPE}, Y=$regressResults{$srcRk}{Y_INTERCEPT}, R2=$regressResults{$srcRk}{R2}, $numPepUsed{$srcRk} peptides]\n";
			}
			print " Done.<BR>\n";
		}
	}


	##############################################
	####>Updating Analysis valid_status & FDR<####
	##############################################
	$validStatus=($endValidation)? 2 : 1;
	my $FDRstring='';
	if ($decoy) {
		my  $numTrueValidPeptides=0;
		my $numValidPeptides=0;
		my $numTrueValidDecoyPeptides=0;
		my $numValidDecoyPeptides=0;
		unless ($manualValidation) { # manual validation not always detected w/ selected peptides, checking rejected ones
			my $pepInfoStrg='';
			foreach my $rank (1..$maxRank) {$pepInfoStrg.=",INFO_PEP$rank";}
			#my $sthRQ=$dbh->prepare("SELECT INFO_PEP1$pepInfoStrg FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID AND QUERY_NUM>0 AND VALID_STATUS<=0 AND VALID_STATUS>=-2");
			my $sthRQ=$dbh->prepare("SELECT QUERY_NUM$pepInfoStrg FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID AND VALID_STATUS>=-2");
			$sthRQ->execute;
			QUERY:while (my ($queryNum,@pepInfo)=$sthRQ->fetchrow_array) {
				foreach my $info (@pepInfo) {
					last unless $info;
					if ($info=~/SEL=(-3|2)/) { # manual selection or rejection
						$manualValidation=1;
						last QUERY;
					}
					elsif ($info=~/SEL=1/) {
						my ($specCount)=($info=~/SPC=(\d+)/); # top + lower-scoring peptides
						$specCount=1 unless $specCount; # old data
						if ($queryNum > 0) {
							$numTrueValidPeptides++;
							$numValidPeptides+=$specCount;
						}
						else {
							$numTrueValidDecoyPeptides++;
							$numValidDecoyPeptides+=$specCount;
						}
					}
				}
			}
			$sthRQ->finish;
		}
		if ($manualValidation) {$FDRstring=',FALSE_DISCOVERY=-2';}
		else {
			$FDRstring=",FALSE_DISCOVERY='$numTrueValidDecoyPeptides:".($numValidDecoyPeptides-$numTrueValidDecoyPeptides).":$numTrueValidPeptides:".($numValidPeptides-$numTrueValidPeptides)."'";
		}
	}
	$dbh->do("UPDATE ANALYSIS SET VALID_STATUS=$validStatus$FDRstring WHERE ID_ANALYSIS=$analysisID") || die $dbh->errstr;
	$dbh->do("UPDATE ANALYSIS SET FIRST_VALID_DATE=VALID_DATE WHERE ID_ANALYSIS=$analysisID AND FIRST_VALID_DATE IS NULL") || die $dbh->errstr;

    &promsMod::updateAnalysisHistory($dbh,$analysisID,'','report');


	############################
	####>Ending validation?<####
	############################
    if ($endValidation) {
		&endValidation($analysisID,$fileFormat,$dataFileName); #,$dataFile,$pepFile,$oldPepFile,$decoyFile,$decoyPepFile,$oldDecoyPepFile
	}

	$dbh->commit;

	#print " Done.</FONT><BR>\n";
}

$sthAnaInfo->finish;
$sthSpecies->finish;
$sthLabMod->finish;
$sthPepRT->finish;
$sthInsART->finish;
$sthInsQRT->finish;
foreach my $refID (keys %referenceRTdata) {
	$referenceRTdata{$refID}{'QUERY'}->finish;
}


##>Update classifications<##
&promsMod::removeFromClassifications($dbh,$projectID) if $prevValidation; # in case prev visible proteins that are not anymore were included in a classification

$dbh->disconnect;

if ($call eq 'cmd') {
	unlink "$promsPath{tmp}/scratch/autoEndValidation/$jobID.flag";
	print " Done.\n";
	exit;
}

################################
####>New identifier Mapping<####
################################
if (scalar @newAnaMapping) { # true if new valid proteins in project
	my $anaStrg=join(',',@analysisList);
	#if ($call eq 'ajax') {
		system "./mapProteinIdentifiers.pl $userID $anaStrg";
	#}
	#else {
	#	print qq|<SCRIPT LANGUAGE="JavaScript">systemFrame.location="./send2Biologist.cgi?MAP=1&ANA_ID=$anaStrg";</SCRIPT>|;
	#}
}

##################
####>END HTML<####
##################
print "<BR><BR><FONT class='title3'>- Reporting is complete.</FONT><BR>\n";

if ($call eq 'ajax') {
	open (P_END,">$promsPath{tmp}/ajaxAnaReport_$userID.end");
	print P_END "__END__\n";
	close P_END;

	print "<!--__OK__--></BODY>\n</HTML>\n";
	exit;
}

sleep 3;
#exit; #debug

my $navUrlString.=($itemInfo[2]{'ITEM'} eq 'SAMPLE')? "&branchID=analysis:$analysisList[0]&ACT=open" : "&branchID=gel2d:$itemInfo[2]{ID}&itemBranchID=analysis:$analysisList[0]&ACT=open";
print qq
|<SCRIPT LANGUAGE="JavaScript">
top.promsFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID$navUrlString";
</SCRIPT>
</BODY>
</HTML>
|;

##################< SUBROUTINES >###################

#####################################
####<Launches identifier mapping>####
#####################################
sub mapIdentifiers {
	my $anaStrg=param('ANA_ID');

	print header(-'content-encoding'=>'no',-charset=>'utf-8'); # start_html,"\n"; # start_html required to force update
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Identifier Mapping</TITLE>
</HEAD>
<BODY>
|;
	system "./mapProteinIdentifiers.pl $userID $anaStrg";
	print qq
|</BODY>
</HTML>
|;
	exit;
}

############################
####<Check delta length>#### Compares delta lengthes between 2 proteins with Peptide density of the smaller one
############################ (delta is considered not significant if smaller than pep density)
sub deltaLength {
	my ($l_a,$l_b,$d_a,$d_b)=@_;
	my $pepDensVal=($l_b > $l_a)? $d_a : $d_b;
	if (abs($l_b-$l_a) > $pepDensVal) {return $l_a<=>$l_b} else {return 0}
}

###########################
####<Ending validation>####
###########################
sub endValidation {
	my ($analysisID,$fileFormat,$dataFileName)=@_; #,$dataFile,$pepFile,$oldPepFile,$decoyFile,$decoyPepFile,$oldDecoyPepFile
	my $validDir="$promsPath{valid}/ana_$analysisID";

	print "&nbsp;&nbsp;&nbsp;- Deleting temporary validation data..." if $call ne 'cmd';

	####<Cleaning-up RANK_PROTEIN_MATCH
	$dbh->do("DELETE FROM RANK_PROTEIN_MATCH WHERE ID_ANALYSIS=$analysisID") || die $dbh->errstr;
	print '.' if $call ne 'cmd';

	####<Cleaning-up QUERY_MODIFICATION,QUERY_VALIDATION
	$dbh->do("DELETE QUERY_MODIFICATION FROM QUERY_MODIFICATION INNER JOIN QUERY_VALIDATION ON QUERY_MODIFICATION.ID_QUERY=QUERY_VALIDATION.ID_QUERY WHERE QUERY_VALIDATION.ID_ANALYSIS=$analysisID") || die $dbh->errstr;
	$dbh->do("DELETE FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID") || die $dbh->errstr;
	print '.' if $call ne 'cmd';

	####<Cleaning-up PROTEIN_VALIDATION
	$dbh->do("DELETE FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$analysisID") || die $dbh->errstr;
	print '.' if $call ne 'cmd';

	####<Deleting useless files from valid dir
	unlink "$validDir/analysis.fasta" if -e "$validDir/analysis.fasta";
	(my $rootName=$dataFileName)=~s/\.[^\.]+\Z//;
	unlink "$validDir/$rootName.target" if -e "$validDir/$rootName.target";
	unlink "$validDir/$rootName.decoy" if -e "$validDir/$rootName.decoy";
	unlink "$validDir/$rootName.qvality" if -e "$validDir/$rootName.qvality";

	####<Decreasing PRS results size
	(my $prsFileName = $dataFileName)=~s/^(.+)\..+/PRS_$1\.xml/;
	if (-e "$validDir/$prsFileName") {
		my $sthPep = $dbh->prepare("SELECT QUERY_NUM FROM PEPTIDE WHERE DATA LIKE \"%PRS=%\" AND ID_ANALYSIS=? AND QUERY_NUM IS NOT NULL");
		$sthPep->execute($analysisID);
		my @queries;
		while (my ($queryNum) = $sthPep->fetchrow_array) {
			push @queries, $queryNum;
		}
		$sthPep->finish;
		phosphoRS::cleanResults("$validDir/$prsFileName", \@queries);
	}


	####<Moving data files to peptide directory
	my $pepDataDir="$promsPath{peptide}/proj_$projectID/ana_$analysisID";
	#remove_tree($pepDataDir) if -e $pepDataDir;
	rmtree($pepDataDir) if -e $pepDataDir;
	dirmove($validDir,$pepDataDir); # entire dir is moved!

	##<Multi Ana
	if ($fileFormat =~ /\.PDM/) {
		(my $msfFile = $dataFileName)=~s/_\d+\.*\d*\.pdm/\.msf/;

		#<Fetch list of already reported msf files
		unless (scalar keys %msfReportedList) { # global. Done only once for all msf files reported
			foreach my $directory (glob("$promsPath{peptide}/proj_$projectID/*")) { # /<full path>/ana_<anaID>
				opendir(DIR,$directory);
				while (defined (my $child = readdir(DIR))) {
					next if $child !~ /\.msf$/;
					if (-f "$directory/$child" && !-l "$directory/$child") { # msf files already stored in project and not as symbolic link
						$msfReportedList{$child}="$directory/$child";
					}
				}
				close DIR;
			}
		}

		(my $anaFile=$msfFile)=~s/\.msf/\.ana/;
		my $useCount;
		if (-e "$promsPath{valid}/multi_ana/proj_$projectID/$anaFile") { # just to be safe
			my $qAnaFile=quotemeta($anaFile);
			$useCount=`grep -c . $promsPath{valid}/multi_ana/proj_$projectID/$qAnaFile`;
			chomp $useCount;
#if (!$useCount && $call eq 'cmd') { # tmp
#	print "\n! Error in multi-ana file: ANA_ID=$analysisID, MSF=$msfFile, ANA_FILE=$anaFile QANA_FILE=$qAnaFile !\n";
#}
		}
		$useCount=1 unless $useCount; # assume msf is not shared (most likely. if shared, will be found in %msfReportedList)

		if ($useCount<=1) { # used only in 1 analysis => move
			###> Check if this MSF is not already in peptide directory
			my $moveOK=($msfReportedList{$msfFile})? symlink($msfReportedList{$msfFile},"$pepDataDir/$msfFile") : 0;
			if (-e "$promsPath{valid}/multi_ana/proj_$projectID/$msfFile") {
				if ($moveOK) {
					unlink("$promsPath{valid}/multi_ana/proj_$projectID/$msfFile");
				}
				else {
					$msfReportedList{$msfFile}="$pepDataDir/$msfFile";
					move("$promsPath{valid}/multi_ana/proj_$projectID/$msfFile",$msfReportedList{$msfFile}); # move msf to peptide dir!
				}
			}
			elsif (!$moveOK) {
				if ($call eq 'cmd') {print '!';}
				else {print "[Warning: msf file not found!]";}
			}
		}
		else { # used by other analys(ie)s => copy
			###> Check if this MSF is not already in peptide directory
			my $copyOK=($msfReportedList{$msfFile})? symlink($msfReportedList{$msfFile},"$pepDataDir/$msfFile") : 0;
			unless ($copyOK) {
				if (-e "$promsPath{valid}/multi_ana/proj_$projectID/$msfFile") {
					$msfReportedList{$msfFile}="$pepDataDir/$msfFile";
					copy("$promsPath{valid}/multi_ana/proj_$projectID/$msfFile",$msfReportedList{$msfFile}); # move msf to peptide dir!
				}
				else {
					if ($call eq 'cmd') {print '!';}
					else {print "[Warning: msf file not found!]";}
				}
			}
		}
		&promsMod::removeFromMultiAna($analysisID,$projectID,$anaFile); # projectID is global
	}
	print "\nDone.<BR>\n" if $call ne 'cmd';

	####<Updating validation history
	&promsMod::updateAnalysisHistory($dbh,$analysisID,'','endVal');
}

#############################################################################
####<Walk up the Processing Nodes tree to looking for parent search node>####
#############################################################################
sub checkSearchNode { # same in convertMsf2Pdm.cgi
	my ($sthPN,$searchNodeNum,$parentNodeStrg)=@_;
	my $match=0;
	foreach my $parNodeNum (split(/;/,$parentNodeStrg)) {
		if ($parNodeNum==$searchNodeNum) {
			$match=1;
			last;
		}
	}
	unless ($match) {
		QNODE:foreach my $parNodeNum (split(/;/,$parentNodeStrg)) {
			$sthPN->execute($parNodeNum);
			while (my ($grandParentNodeStrg)=$sthPN->fetchrow_array) {
				$match=&checkSearchNode($sthPN,$searchNodeNum,$grandParentNodeStrg);
				last QNODE if $match;
			}
		}
	}
	return $match;
}


####################################################
####<Computing RT regression/correlation with R>####
####################################################
sub runRegressionR { # Globals: %promsPath
	my ($anaDir,$refRT_ID,$pepQuantifID,$refVectors,$refMergedFiles)=@_;

	my $rScript="$anaDir/regression_RT$refRT_ID"."_Q$pepQuantifID.R";
	my $exportFile="$anaDir/regression_RT$refRT_ID"."_Q$pepQuantifID.txt";
	my $correlImg="$anaDir/regression_RT$refRT_ID"."_Q$pepQuantifID.png";

	my $numFiles=scalar keys %{$refMergedFiles};
	my $fontSize=($numFiles==1)? 5 : ($numFiles==2)? 4 : ($numFiles<=4)? 5 : ($numFiles<=9)? 4 : 3;
	open Rscript, ">$rScript";
	print Rscript qq
|require(ggplot2)

intercept <- c()
slope <- c()
eqX <- c()
eqY <- c()
eqL <- c()
r2 <- c()
|;
	my (@fileList,@refVectList,@anaVectList,%file2Rank); #@sourceList,%infoFile,
	my %results;
	my $okDraw=0;
	foreach my $fileID (sort{$a<=>$b} keys %{$refVectors}) {
		next unless scalar (@{$refVectors->{$fileID}{'ANA'}}) >= 2; # at least 2 peptides for regression
		$okDraw=1;
		my ($sourceFile,$sourceRank)=@{$refMergedFiles->{$fileID}};
		$file2Rank{$sourceFile}=$sourceRank;
		print Rscript "\n##>File $sourceRank: $sourceFile\n";
		#$infoFile{$sourceRank}="$anaDir/regression_RT$refRT_ID"."_Q$pepQuantifID"."_$sourceRank.txt";
		print Rscript "x_$sourceRank <- c(",join(',',@{$refVectors->{$fileID}{'ANA'}}),")\n";
		print Rscript "y_$sourceRank <- c(",join(',',@{$refVectors->{$fileID}{'REF'}}),")\n";
		print Rscript qq
|eq_$sourceRank <- lm(y_$sourceRank~x_$sourceRank)
eqstg_$sourceRank <- substitute(
	italic(Ref) == a_$sourceRank + b_$sourceRank %.% italic(Obs)*","~~italic(R)^2~"="~r2_$sourceRank,
	list(
		a_$sourceRank = format(coef(eq_$sourceRank)[1], digits = 4),
		b_$sourceRank = format(coef(eq_$sourceRank)[2], digits = 4),
		r2_$sourceRank = format(summary(eq_$sourceRank)\$r.squared, digits = 3)
	)
)
intercept <- c(intercept, coef(eq_$sourceRank)[1])
slope <- c(slope, coef(eq_$sourceRank)[2])
eqX <- c(eqX, min(x_$sourceRank))
eqY <- c(eqY, max(y_$sourceRank))
eqL <- c(eqL, as.character(as.expression(eqstg_$sourceRank)))
r2 <- c(r2, summary(eq_$sourceRank)\$r.squared)
|;
		#push @rankList,($sourceRank) x scalar @{$refVectors->{$fileID}{'ANA'}};
		push @fileList,('"'.$sourceFile.'"') x scalar @{$refVectors->{$fileID}{'ANA'}};
		push @refVectList,@{$refVectors->{$fileID}{'REF'}};
		push @anaVectList,@{$refVectors->{$fileID}{'ANA'}};
	}
	unless ($okDraw) {
		close Rscript;
		unlink $rScript;
		return %results;
	}
	print Rscript "\n##>Data frame #1\n";
	#print Rscript "rank <- c(",join(',',@rankList),")\n";
	print Rscript "file <- c(",join(',',@fileList),")\n";
	print Rscript "Observed_RT <- c(",join(',',@anaVectList),")\n";
	print Rscript "Reference_RT <- c(",join(',',@refVectList),")\n";
	print Rscript qq
|rtData <- data.frame(file, Observed_RT, Reference_RT)

##>Data frame #2
file <- unique(file)
lineInfo <- data.frame(file, intercept, slope, eqL, eqX, eqY)
exportData <- data.frame(file, intercept, slope, r2)

##>Plot
png(filename="$correlImg", width=1000, height=600, units = "px")
p <- ggplot(rtData, aes(Observed_RT, Reference_RT))
p + facet_wrap(~ file, scales="free") + stat_smooth(method="lm", se=TRUE) + geom_abline(data = lineInfo, aes(intercept=intercept, slope=slope), color="red", size=1) + geom_point() + geom_text(data=lineInfo, aes(label=eqL, x=eqX, y=eqY, hjust=0, size=$fontSize), parse=TRUE, show_guide=FALSE) #, vjust=0
dev.off()

##>Export data
write.table(exportData, file="$exportFile",sep=c('\\t'),quote=F,row.names=F,col.names=T)
|;
	close Rscript;

	system "cd $anaDir; $promsPath{R}/R CMD BATCH --no-save --no-restore $rScript";

	open(RES,$exportFile) || warn("**Error while opening $exportFile !**");
	while (my $line=<RES>) {
		next if $.==1;
		chomp($line);
		my ($sourceFile,$y_intercept,$slope,$rsquared)=split(/\t/,$line);
		%{$results{$file2Rank{$sourceFile}}}=('Y_INTERCEPT'=>$y_intercept,'SLOPE'=>$slope,'R2'=>$rsquared);
	}
	close RES;
	unlink($rScript,$rScript.'out',$exportFile);

	return %results;
}


####>Revision history<####
# 3.3.4 [MODIF] Restore lost CeCILL license text (PP 28/08/19)
# 3.3.3 [MODIF] Changed createMatchGroup to fit with promsMod.pm prototype (VS 03/07/19)
# 3.3.2 Removed all usages of PEPTIDE_QUANTIFICATION table (PP 27/06/18)
# 3.3.1 Minor bug correction for TMT (GA 18/06/18)
# 3.3.0 Peptide quantification data now written to file $promsPath{quantification}/project_$projectID/quanti_$quantifID/peptide_quantification(_$targetPos).txt (PP 24/05/18)
# 3.2.8 Launchable in command line mode by selectProject to auto-end validation of old analyses (PP 24/05/18)
# 3.2.7 Better matching for reference RT kit in case of identifier ambiguity (PP 30/04/18)
# 3.2.6 Update auto-increment PROTEIN and QUANTIFICATION insertions (GA 21/03/18)
# 3.2.5 Speeding up the queries for SILAC extraction with PD 2.2 (GA 12/03/18)
# 3.2.4 Minor modification for PD 2.2 and merged files (GA 08/12/17)
# 3.2.3 Minor modification for PD 2.0.0.802 TMT version bug (GA 02/11/17)
# 3.2.2 Add print '.' to large TMT files (GA 23/10/17)
# 3.2.1 Update to PD 2_2 version and works for SILAC, TMT, iTRAQ and label free. (GA 11/09/17)
# 3.2.0 Compatible with MODIFICATION_SITE (PP 02/08/17)
# 3.1.9 Minor modification for QUANTIF_ANNOT in TMT labeling<BR>TODO: make it compatible with Mascot/Paragon/Andromeda/X!Tandem (GA 03/04/2017)
# 3.1.8 TMT storing from PD 2.1 msf files (GA 29/03/2017)
# 3.1.7 Modify explicit error for label-modifications to check if the modification exist on the specific residue (GA 08/12/2016)
# 3.1.6 Add explicit error for label-modifications that are not defined in dynamic mofication in PD Workflows (GA 21/11/2016)
# 3.1.5 Change to deal with autoincrement in PEPTIDE table (GA 19/10/2016)
# 3.1.4 Update for matching and reassembling reference RT kit protein fragments <iRT_dentifier>-n (PP 24/06/16)
# 3.1.3 Fix bug duplicate iRT peptide selection & skip excluded iRT peptides (PP 07/04/16)
# 3.1.2 Extends iRT detection to non-quantif analysis & Experiment-level preferred species (PP 07/03/16)
# 3.1.1 Minor change in iRT detection to prevent matching normal identifiers (PP 03/02/16)
# 3.1.0 New management of PD SILAC quantification (PP 02/09/15)
# 3.0.0 Merged with 2.9.8 GA & iRT managment (PP 11/08/15)
# 2.9.8 Add label free for PD 2.0 (GA 24/06/15)
# 2.9.7 Add new SQLite queries for PD 2.0 version and add one package to handle HTML coding in XML workflow (GA 12/06/15)
# 2.9.6PP3 Fix bug for modif specificity with msf file with quantif when labeling specificity is declared on multiple lines (PP 16/07/15)
# 2.9.6PP2 New PD merged search files management & report RT_APEX for SILAC quantification (PP 13/07/15)
# 2.9.6PP Apply FDR filtering to ghost peptides && PARAGON iTRAQ quantif && Proteome Discoverer peptide XIC (PP ../04/15)
# 2.9.6 Change quantitation for msf-split version and update endValidation  (GA 30/03/15)
# 2.9.5 Minor modif due to split-mode that changes $msfFile (GA 27/02/15)
# 2.9.4 Minor modif to get Sequest HT or SEQUEST searches from PD1.4 (GA 28/11/14)
# 2.9.3 Check if MSF quantification was actually performed before extracting data (PP 17/09/14)
# 2.9.2 PD version explicitly retrieved from first entry in SchemaInfo table (PP 12/05/14)
# 2.9.1 Records number of validatable peptides including lower-scoring ones (PP 14/04/14)
# 2.9.0 Loads list of modificationIDs from ANALYSIS_MODIFICATION & improved match between quantification and search nodes (PP 01/04/14)
# 2.8.9 Handles no-label channel with poor XML description (PP 03/04/14)
# 2.8.8 Fix bug protein mapping call (PP 11/03/14)
# 2.8.7 Uses rmtree instead of remove_tree (PP 10/03/14)
# 2.8.6 Called by editBranch through ajax instead of LWP (PP 08/03/14)
# 2.8.5 ??
# 2.8.4 Bug fix in SEQUEST search query & system command removal (PP 12/11/13)
# 2.8.3 Bug fix in FIRST_VALID_DATE field naming (PP 23/10/13)
# 2.8.2 iTRAQ -> ITRAQ for LABEL value in QUANTIF_ANNOT field (PP 07/10/13)
# 2.8.1 No more check for complete quantif sets (PP 30/09/13)
# 2.8.0 Comment a part of the code to get SILAC information even though all isotopes are not available<BR>Minor modification in Processing Quantification settings (MSF) for SILAC to get $deltaMass and avoid warnings (GA 23/09/13)
# 2.7.9 Quantification deletion even if no pre-reported proteins (PP 19/09/13)
# 2.7.8 Fixed critical bug (FY 19/09/13)
# 2.7.7 PEPTIDE_SET tables no longer used (PP 16/09/13)
# 2.7.6 Minor change in multi-ana msf file management during end of validation (PP 13/09/13)
# 2.7.5 Merge of 2.7.3b & 2.7.4 (PP 25/07/13)
# 2.7.3b Bug fix for label detection in MASCOT.PDM files & update to change due to storeAnalyses.cgi 3.2.2f (PP 17/07/13)
# 2.7.3 Bug Fix for new modification management & better label detection & speed optimization (PP 19/06/13)
# 2.7.2 Remove VAR_MOD, VMOD and RMOD from script and insert PEPTIDE_MODIFICATION for SILAC virtual peptides (GA 31/05/13)
# 2.7.1 New PEPTIDE_MODIFICATION info (GA 22/04/13)
# 2.7.0 Reducing PRS result file size when ending validation (FY 11/04/13)
# 2.6.9 Add REF_MOD tag in PEPTIDE.DATA for recording unmodified VAR_MOD (PP 25/03/13)
# 2.6.8 Add Substitution for PARAGON (GA 21/03/13)
# 2.6.7 Handles new GI accession restriction on alias (PP 12/03/13)
# 2.6.6 No file transfert before end of validation, No spectrum extraction from msf & new peptide dir (PP 07/03/13)
# 2.6.5 Ghost peptides with missCut,-beg,-end,flankingAA (PP 06/02/13)
# 2.6.4 Flag on peptide coverage value if peptide position &gt protein length (PP 23/01/13)
# 2.6.3 Missed cut info for ghost peptides & additional multi-databank search compatibility (PP 16/01/12)
# 2.6.2 Multi-databank search compatibility (PP 12/12/12)
# 2.6.1 To make it work with PARAGON files (GA 11/12/12)
# 2.6.0 Handles lwp call from editBranch.cgi for ending validations<BR>& New identifier mapping procedure (PP 19/11/12)
# 2.5.3 Reports Qvality data & new attempt to prevent server time out at "Processing data file..." step (PP 07/08/12)
# 2.5.2 Prevents server disconnection at "Processing data file..." step with large PDM files (PP 06/06/12)
# 2.5.1 Symbolic link is preserved if validation file is already a link (PP 03/05/12)
# 2.5.0 Modification of queries for $sthVarMods and $sthAllVM to make it works for both PD 1.2 and 1.3 version (GA 03/05/12)
# 2.4.9 PEPTIDE.DATA field contains PRS & some SILAC data (PP 27/04/12)
# 2.4.8 Uses [ISO/REP]_AREA for SILAC & ITRAQ instead of XIC_AREA (PP 12/04/12)
# 2.4.7 Management of PhosphoRS data<BR>QUANTI_DATA column renamed to DATA (FY 29/03/12)
# 2.4.6 Update MSF SILAC management to skip incomplete quantif sets (PP 18/01/12)
# 2.4.5 Management of SILAC data from (merged-) MSF & iTRAQ data from DAT (PP 07/07/11)
# 2.4.4 Merge 2.4.2 & 2.4.3 (PP 27/04/11)<BR>& full datafile conservation (Mascot only) management (PP 01/06/11)
# 2.4.3 Correct the bug "Fasta file not found!" in report validation option (GA 02/05/11)
# 2.4.2 Improved protein ordering MG (PP 27/04/11)
# 2.4.1 Add "End validation" process to validation history (FY 18/04/11)
# 2.4.0 Adding updating analysis history (FY 02/2011)
# 2.3.9 Fix bug sort protLength for match group
# 2.3.8 Extract MSF w/o AJAX, Mixed decoy DB
