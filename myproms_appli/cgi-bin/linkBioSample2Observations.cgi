#!/usr/local/bin/perl -w

################################################################################
# linkBioSample2Observations.cgi     1.2.0                                     #
# Authors: P. Poullet, G. Arras, S.Liva (Institut Curie)              	       #
# Contact: myproms@curie.fr                                                    #
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
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use strict;

#print header; warningsToBrowser(1);#DEBUG

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};

#################################
####>Fetching Parameters... <####
#################################
my $action=param ('ACT') || 'edit';
my $projectID = param('projectID');
my $experimentID=param('expID') || 0;
my $biosampleID = param('biosampleID') || 0; # can be 0
($projectID,$experimentID,$biosampleID)=&promsMod::cleanNumericalParameters($projectID,$experimentID,$biosampleID);

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

##User Info and access
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};


#########################
####>Processing Form<####
#########################
if ($action eq 'store') {
	print header(-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Biological Sample Observation Link Management</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
|;

	my $sthUpSamp=$dbh->prepare("UPDATE OBSERVATION SET ID_BIOSAMPLE=? WHERE ID_OBSERVATION=?");
	my $sthAddObs=$dbh->prepare("INSERT INTO OBSERVATION (ID_ANALYSIS,ID_BIOSAMPLE,TARGET_POS) VALUES (?,?,?)");
	my $sthAddOMod=$dbh->prepare("INSERT INTO OBS_MODIFICATION (ID_OBSERVATION,ID_MODIFICATION) VALUES (?,?)");
	my %freeObservations;

	####>Called from a biological sample<####
	if ($biosampleID) {

		###>Cleaning Experiment's existing Observations
		foreach my $fullObsCode (param('obsList')) {
			if ($fullObsCode=~/:(\d+)/) { # true obs
				my $obsID=$1;
				$sthUpSamp->execute(undef,$obsID);
				$freeObservations{$obsID}=1;
			}
		}

		###>Updating/Creating Observations
		foreach my $fullObsCode (param('selObs')) {
			my ($obsCode,$obsID)=split(/:/,$fullObsCode);
			if ($obsID) { # Update Obs
				$sthUpSamp->execute($biosampleID,$obsID);
				delete $freeObservations{$obsID}; # remove from list of unused Obs
			}
			else { # Create Obs
				my ($anaID,$targetPos,@mods)=split(',',$obsCode);
				$sthAddObs->execute($anaID,$biosampleID,$targetPos);
				my $newObsID=$dbh->last_insert_id(undef,undef,'OBSERVATION','ID_OBSERVATION');
				foreach my $modID (@mods) { # can be empty
					$sthAddOMod->execute($newObsID,$modID);
				}
			}
		}
	}

	####>Called from the tree top branch<####
	else {
		my @selBioSamples=param('bioSamples');
		my $obsIndex=0;
		foreach my $fullObsData (param('obsList')) {
			my ($obsCode,$obsID,$prevBioSampID)=split(/[:=]/,$fullObsData);
			my $selBioSampID=$selBioSamples[$obsIndex];
			###>Add/Change biosample link
			if ($selBioSampID) { # a biosample is linked
				if ($obsID) { # true Obs
					if ($selBioSampID != $prevBioSampID) { # bioSample was changed or add (prev=0): update Obs
						$sthUpSamp->execute($selBioSampID,$obsID);
					} # else no change: do nothing
				}
				else { # Create Obs
					my ($anaID,$targetPos,@mods)=split(',',$obsCode);
					$sthAddObs->execute($anaID,$selBioSampID,$targetPos);
					my $newObsID=$dbh->last_insert_id(undef,undef,'OBSERVATION','ID_OBSERVATION');
					foreach my $modID (@mods) {
						$sthAddOMod->execute($newObsID,$modID);
					}
				}
			}
			###>Remove biosample link if any
			else { # no biosample link
				if ($obsID) { # true Obs
					if ($prevBioSampID) { # remove pre-existing bioSample
						$sthUpSamp->execute(undef,$obsID);
						$freeObservations{$obsID}=1; # add to list of unused Obs
					} # else no change: do nothing
				} # else fake Obs: do nothing
			}
			$obsIndex++;
		}
	}

	$sthUpSamp->finish;
	$sthAddObs->finish;
	$sthAddOMod->finish;


	####>Deleting any useless Observations<####
	if (scalar keys %freeObservations) {
		my $sthObsCond=$dbh->prepare("SELECT 1 FROM OBS_EXPCONDITION WHERE ID_OBSERVATION=? LIMIT 1");
		my $sthDelOMod=$dbh->prepare("DELETE FROM OBS_MODIFICATION WHERE ID_OBSERVATION=?");
		my $sthDelObs=$dbh->prepare("DELETE FROM OBSERVATION WHERE ID_OBSERVATION=?");
		foreach my $obsID (keys %freeObservations) {
			$sthObsCond->execute($obsID); # check if Obs is used by an ExpCondition
			my ($usedObs)=$sthObsCond->fetchrow_array;
			unless ($usedObs) {
				$sthDelOMod->execute($obsID);
				$sthDelObs->execute($obsID);
			}
		}
		$sthObsCond->finish;
		$sthDelOMod->finish;
		$sthDelObs->finish;
	}

	$dbh->commit;
	$dbh->disconnect;

	if ($biosampleID) {
		print qq
|<SCRIPT type="text/javascript">
top.promsFrame.selectedAction = 'summary';
parent.optionFrame.selectOption();
</SCRIPT>
</HEAD>
</HTML>
|;
	}
	else {
		print qq
|<CENTER>
<BR><BR>
<FONT class="title2">Links to Observations successfully updated!</FONT>&nbsp;&nbsp;<INPUT type="button" class="title3" value=" Continue " onclick="parent.itemFrame.location.reload(); parent.optionFrame.selectOption();"/>
</BR></BR>
</BODY>
</HTML>
|;
	}
    exit;
}

###################################################
####>Propagating links from another Experiment<####
###################################################
elsif ($action eq 'propagate') {
	
	my $targetExpID=&promsMod::cleanNumericalParameters(param('targetExp'));
	my $matchType=param('matchType') || 'analysis';
	
	print header(-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Biological Sample Observation Link Management</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
|;
	my ($expName)=$dbh->selectrow_array("SELECT NAME FROM EXPERIMENT WHERE ID_EXPERIMENT=$targetExpID");
	
	my $numNewLinks=0;
	my %linkedBioSamples;
	my $sthUpSamp=$dbh->prepare("UPDATE OBSERVATION SET ID_BIOSAMPLE=? WHERE ID_OBSERVATION=?");
	
	###>Match based on Ana name or WIFF file<###
	if ($matchType=~/analysis|datafile/) {
		
		my $keyField=($matchType eq 'analysis')? 'A.NAME' : 'A.WIFF_FILE';
		
		###>Source Experiment
		my $sthStartObs=$dbh->prepare("SELECT $keyField,O.ID_OBSERVATION,O.ID_BIOSAMPLE,TARGET_POS,GROUP_CONCAT(ID_MODIFICATION ORDER BY ID_MODIFICATION SEPARATOR ',')
										FROM OBSERVATION O
										LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION
										INNER JOIN ANALYSIS A ON O.ID_ANALYSIS=A.ID_ANALYSIS
										INNER JOIN SAMPLE S ON A.ID_SAMPLE=S.ID_SAMPLE
										WHERE S.ID_EXPERIMENT=? AND O.ID_BIOSAMPLE IS NOT NULL GROUP BY O.ID_OBSERVATION");
		my %obs2BioSample;
		$sthStartObs->execute($experimentID);
		while (my ($keyFieldValue,$obsID,$bioSampID,$targetPos,$modCode)=$sthStartObs->fetchrow_array) {
			my $obsCode="$keyFieldValue,$targetPos";
			$obsCode.=",$modCode" if $modCode;
			@{$obs2BioSample{$obsCode}}=($bioSampID,$targetPos);
			$obs2BioSample{$obsCode}[2]=$modCode if $modCode;
		}
		$sthStartObs->finish;
	
		###>Target Experiment
		my $sthGAI=$dbh->prepare("SELECT ID_ANALYSIS,$keyField FROM ANALYSIS A,SAMPLE S WHERE S.ID_SAMPLE=A.ID_SAMPLE AND VALID_STATUS > 0 AND S.ID_EXPERIMENT=?");
		my $sthLQ=$dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND FOCUS='peptide' AND ID_ANALYSIS=? LIMIT 1");
		my $sthLA=$dbh->prepare("SELECT DISTINCT M.ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES FROM MODIFICATION M,ANALYSIS_MODIFICATION AM WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND IS_LABEL=1 AND ID_ANALYSIS=? ORDER BY M.ID_MODIFICATION");
		my $sthObs=$dbh->prepare("SELECT O.ID_OBSERVATION,O.ID_BIOSAMPLE,O.TARGET_POS,GROUP_CONCAT(ID_MODIFICATION ORDER BY ID_MODIFICATION SEPARATOR ',')
									FROM OBSERVATION O
									LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION
									WHERE O.ID_ANALYSIS=? GROUP BY O.ID_OBSERVATION");
		my $sthAddObs=$dbh->prepare("INSERT INTO OBSERVATION (ID_ANALYSIS,ID_BIOSAMPLE,TARGET_POS) VALUES (?,?,?)");
		my $sthAddOMod=$dbh->prepare("INSERT INTO OBS_MODIFICATION (ID_OBSERVATION,ID_MODIFICATION) VALUES (?,?)");
		
		my %labelMods;
		$sthGAI->execute($targetExpID);
		while (my ($anaID,$keyFieldValue) = $sthGAI->fetchrow_array) {
			##>Fetch all possible quantif channels
			my %possibleObs;
			$sthLQ->execute($anaID);
			my ($quantifAnnot)=$sthLQ->fetchrow_array;
			if ($quantifAnnot) {
				$quantifAnnot=~s/::SOFTWARE=[^:]+::/::/; # remove software info for back compatibility
				my ($labelTypeStrg,@labelInfo)=split('::',$quantifAnnot);
				my ($labelType)=($labelTypeStrg =~ /^LABEL=(.+)/);
				if ($labelType =~ /FREE|NONE/) { # no labeling in analysis
					$possibleObs{"$keyFieldValue,0"}=1;
				}
				else {
					my $labelModID;
					foreach my $infoStrg (@labelInfo) {
						my ($targetPos,$chanName,$labelStrg)=split(';',$infoStrg);
						my $obsCode; #="$anaID,$targetPos";
						my $targetName;
						if ($labelType=~/SILAC/i) {
							$obsCode="$keyFieldValue,-1";
							#my @mods;
							my %mods; # hash because same mod can be used on different aa (mod is repeated)
							foreach my $modStrg (split ('@',$labelStrg)) { # multiple mods for same label channel
								my @labelData=split('#',$modStrg);
								next if $labelData[1] =~ /No label/i;
								my $modID=($labelData[4])? $labelData[4] : ($labelMods{$labelData[1]})? $labelMods{$labelData[1]} : &promsMod::getModificationIDfromString($dbh,$labelData[1]); # $labelData[4] if MassChroQ XIC
								$modID=0 unless $modID;
								$labelMods{$labelData[1]}=$modID;;
								$mods{$modID}=1;
							}
							if (scalar keys %mods) { # label channel
								$obsCode.=','.join(',',sort{$a<=>$b} keys %mods); # order is needed for later compare with %expCondition
							}
						}
						elsif ($labelType=~/itraq|TMT/i) {
							$obsCode="$keyFieldValue,$targetPos";
							unless ($labelModID) { # only once
								$sthLA->execute($anaID);
								($labelModID)=$sthLA->fetchrow_array;
							}
							$obsCode.=",$labelModID" if $labelModID;
						}
						$possibleObs{$obsCode}=1;
					}
				}
			}
			##>No peptide quantification associated with Analysis => guess from ANALYSIS_MODIFICATION table ****WARNING: NOT reliable
			else {
				$sthLA->execute($anaID);
				my $isLabeled=0;
				my $obsCode;
				my $channelPos=0;
				while (my ($modID,$psiName,$intName,$synName)=$sthLA->fetchrow_array) {
					$psiName='' unless $psiName;
					$intName='' unless $intName;
					$synName='' unless $synName;
					my $labelName=$psiName.'##'.$intName.'##'.$synName;
					if ($labelName=~/(itraq|tmt).+plex/i) { # Cannot work for SILAC
						$isLabeled=1;
						my %targetInfo;
						my $targetPos=0;
						my %infoSubmit=&promsMod::getSearchParam($dbh,$anaID);
						if ($infoSubmit{'q:Quantification'}) {
							my $first=1;
							foreach my $targetName (split(':',$infoSubmit{'q:Quantification'})) {
								if ($first) { # skip quantif Name
									$first=0;
									next;
								}
								$targetInfo{++$targetPos}=$targetName;
							}
						}
						else { # backup solution
							my ($plexNum)=($labelName=~/(\d+)\s*plex/i);
							foreach my $targetPos (1..$plexNum) {
								$targetInfo{$targetPos}='XXX';
							}
						}
						foreach my $targetPos (sort{$a<=>$b} keys %targetInfo) {
							$possibleObs{"$keyFieldValue,$targetPos,$modID"}=1;
						}
						next;
					}
					unless ($isLabeled) { # 1st loop: create a non-labeled channel
						#$channelPos++;
						$possibleObs{"$keyFieldValue,$channelPos"}=1;
						$isLabeled=1;
					}
					# !!!!!!!!!! TODO: Not compatible with single state defined by multiple modifications !!!!!!!
					$channelPos++;
					$possibleObs{"$keyFieldValue,$channelPos,$modID"}=1;
				}
				unless ($isLabeled) { # no labeling in analysis
					$possibleObs{"$keyFieldValue,0"}=1;
				}
			}
				
			##>Check if Ana already linked to obs
			$sthObs->execute($anaID);
			while (my ($obsID,$bioSampID,$targetPos,$modCode)=$sthObs->fetchrow_array) {
				#$targetPos=-1 if ($quantifAnnot && $quantifAnnot=~/SILAC/); # back compatibility for Design before 31/01/14
				my $obsCode="$keyFieldValue,$targetPos";
				$obsCode.=",$modCode" if $modCode;
				delete $possibleObs{$obsCode}; # remove from list to be created
				next if $bioSampID;
				if ($obs2BioSample{$obsCode}) {
					$sthUpSamp->execute($obs2BioSample{$obsCode}[0],$obsID); # bioSampID
					$numNewLinks++;
					$linkedBioSamples{$obs2BioSample{$obsCode}[0]}=1;
				}
			}
			##>Process remaining possible obs
			foreach my $obsCode (keys %possibleObs) {
				if ($obs2BioSample{$obsCode}) {
					$sthAddObs->execute($anaID,@{$obs2BioSample{$obsCode}}[0,1]); # bioSampID,targetPos
					if ($obs2BioSample{$obsCode}[2]) { # modID strg
						my $newObsID=$dbh->last_insert_id(undef,undef,'OBSERVATION','ID_OBSERVATION');
						foreach my $modID (split(',',$obs2BioSample{$obsCode}[2])) {
							$sthAddOMod->execute($newObsID,$modID);
						}
					}
					$numNewLinks++;
					$linkedBioSamples{$obs2BioSample{$obsCode}[0]}=1;
				}
			}
		}
		$sthLQ->finish;
		$sthLA->finish;
		$sthObs->finish;
		$sthAddObs->finish;
		$sthAddOMod->finish;
	}
	
	###>Match based on Design states name<###
	else {
		
		my $sthStateObs=$dbh->prepare("SELECT D.NAME,EC.NAME,O.ID_BIOSAMPLE,FRACTION_GROUP,TECH_REP_GROUP,O.ID_OBSERVATION
										FROM DESIGN D
										INNER JOIN EXPCONDITION EC ON D.ID_DESIGN=EC.ID_DESIGN
										INNER JOIN OBS_EXPCONDITION OEC ON EC.ID_EXPCONDITION=OEC.ID_EXPCONDITION
										INNER JOIN OBSERVATION O ON OEC.ID_OBSERVATION=O.ID_OBSERVATION
										WHERE D.ID_EXPERIMENT=?");
		###>Source Experiment
		my %sourceBioSamples;
		$sthStateObs->execute($experimentID);
		while (my ($desName,$condName,$bioSampID,$fractGroup,$techRepGroup,$obsID)=$sthStateObs->fetchrow_array) {
			next unless $bioSampID;
			$fractGroup=1 unless $fractGroup;
			$techRepGroup=1 unless $techRepGroup;
			$sourceBioSamples{"$condName,$fractGroup,$techRepGroup"}{$bioSampID}=$desName; # in case same state name linked to different biosamp depending on design (unlikely?)
		}
		
		###>Target Experiment
		$sthStateObs->execute($targetExpID);
		while (my ($desName,$condName,$tgtBioSampID,$fractGroup,$techRepGroup,$obsID)=$sthStateObs->fetchrow_array) {
			next if $tgtBioSampID; # keep only obs w/o bioSamp
			$fractGroup=1 unless $fractGroup;
			$techRepGroup=1 unless $techRepGroup;
			my $nameCode="$condName,$fractGroup,$techRepGroup";
			next unless $sourceBioSamples{$nameCode}; # not found in source Exp
			my $bioSampID;
			if (scalar keys %{$sourceBioSamples{$nameCode}} > 1) { # multiple bioSamples linked same state name+fracGr+techGr ==> check design name
				foreach my $bsID (keys %{$sourceBioSamples{$nameCode}}) {
					if ($sourceBioSamples{$nameCode}{$bsID} eq $desName) {
						$bioSampID=$bsID;
						last;
					}
				}
			}
			else { # only 1 bioSample 
				$bioSampID=(keys %{$sourceBioSamples{$nameCode}})[0];
			}
			next unless $bioSampID;
			$sthUpSamp->execute($bioSampID,$obsID);
			$numNewLinks++;
			$linkedBioSamples{$bioSampID}=1;
		}
		$sthStateObs->finish;
	}
	$sthUpSamp->finish;
	
	$dbh->commit;
	$dbh->disconnect;
	
	my $numBioSamp=scalar keys %linkedBioSamples;
	print qq
|<BR><BR>
<FONT class="title2">$numNewLinks links to Observations ($numBioSamp BioSamples) successfully propagated to Experiment <FONT color="#DD0000">$expName</FONT>!</FONT>
<BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="button" class="title2" value=" Continue " onclick="parent.optionFrame.selectOption();"/>
</BR></BR>
</BODY>
</HTML>
|;
	exit;
}

###################################
####>Importing links from file<####
###################################
elsif ($action eq 'import') {
	my $dataFile = tmpFileName(upload('dataFile'));
	
	my ($expName)=$dbh->selectrow_array("SELECT NAME FROM EXPERIMENT WHERE ID_EXPERIMENT=$experimentID");
	
	print header(-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Biological Sample Observation Link Management</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR><FONT class="title">Importing Biological Sample - Observation Link from file for Experiment <FONT color="#DD0000">$expName</FONT></FONT><BR><BR>
</CENTER>
|;
	
	####>Parsing file<####
	print "&bull; <FONT class=\"title3\">Reading data file...";
	my (%bioSampleLinks,%bioSamplePos,%anaTargets,%analysisLabel);
	my $fileHead=`head -2 $dataFile`;
	my $newLine='';
	$newLine.='CR' if $fileHead=~/\r/;
	$newLine.='LF' if $fileHead=~/\n/;
	local $/=($newLine eq 'CR')? "\r" : ($newLine eq 'CRLF')? "\r\n" : "\n"; # Perl new line separator 
	my $hasLabel=0;
	my $prevBsName='';
	open(DATA,$dataFile) || die "Unable to open $dataFile";
	while (my $line=<DATA>) {
		next if $.==1; # 1st line of the file
		$line=~s/\s+$//; # chomp is not enough <- Windows
		my ($bsName,$anaLink,$labelCode)=split(/ *\t */,$line); # remove starting/trailing spaces
		next unless $anaLink; # empty line or bug
		$bsName=$prevBsName unless $bsName; # allow empty field if same BS is usued in multiple adjacent rows
		$bsName=~s/^ +//; # remove any leading space
		$labelCode='None' unless $labelCode;
		push @{$bioSampleLinks{$bsName}},[$anaLink,$labelCode];
		$bioSamplePos{$bsName}=$. unless $bioSamplePos{$bsName};
		$anaTargets{$anaLink}=1;
		#if ($labelCode) {
		#	$hasLabel=1;
			$analysisLabel{$anaLink}{$labelCode}=1;
		#}
		$prevBsName=$bsName;
	}
	close DATA;

	my $numAllBioSamp=scalar keys %bioSampleLinks;
	my $numAllAna=scalar keys %anaTargets;
	if ($numAllBioSamp) {
		print " Done: $numAllBioSamp Biological samples and $numAllAna Analyses listed.</FONT><BR>\n";
	}
	else {
		print "</FONT><BR>";
		$dbh->disconnect;
		&exitOnError('No Biological samples found');
	}
	
	####>Mapping BioSample name to ID_BIOSAMPLE<####
	my %bioSampleIDs;
	my $sthAllBS=$dbh->prepare("SELECT BS.ID_BIOSAMPLE,BS.NAME FROM BIOSAMPLE BS INNER JOIN PROJECT_BIOSAMPLE PB ON PB.ID_BIOSAMPLE=BS.ID_BIOSAMPLE WHERE PB.ID_PROJECT=$projectID AND BS.NAME IN ('".(join("','",keys %bioSampleLinks))."')");
	$sthAllBS->execute;
	while (my ($bioSampID,$bsName)=$sthAllBS->fetchrow_array) {
		$bioSampleIDs{$bsName}=$bioSampID;
	}
	$sthAllBS->finish;
	
	###>Checking if everything was mapped<###
	my $numMapBioSamp=scalar keys %bioSampleIDs;
	if ($numMapBioSamp) {
		if ($numMapBioSamp < $numAllBioSamp) {
			my $bsStrg=($numAllBioSamp-$numMapBioSamp == 1)? 'Biological sample was' : ($numAllBioSamp-$numMapBioSamp)." Biological samples were";
			print "&bull; <FONT class=\"title3\">The following $bsStrg not found in current Project:<BR>\n";
			foreach my $bsName (sort{&promsMod::sortSmart($a,$b)} keys %bioSampleLinks) {
				next if $bioSampleIDs{$bsName};
				print "&nbsp;&nbsp;&nbsp;-$bsName<BR>\n";
			}
			print "</FONT><BR>\n";
		}
	}
	else {
		$dbh->disconnect;
		&exitOnError('Biological samples from file were not found in current Project');
	}
	
	####>Mapping Analyses<####
	#my %analysisIDs;
	#my $analysisTgtStrg="'".(join("','",keys %anaTargets))."'";
	#my $sthAna=$dbh->prepare("SELECT A.ID_ANALYSIS,A.NAME,DATA_FILE,WIFF_FILE,LAB_CODE FROM ANALYSIS A
	#							INNER JOIN SAMPLE S ON S.ID_SAMPLE=A.ID_SAMPLE
	#							WHERE S.ID_EXPERIMENT=$experimentID AND A.NAME IN ($analysisTgtStrg) OR DATA_FILE IN ($analysisTgtStrg) OR WIFF_FILE IN ($analysisTgtStrg) OR LAB_CODE IN ($analysisTgtStrg)");
	#$sthAna->execute;
	#while (my ($anaID,$name,$dFile,$wFile,$lCode)=$sthAna->fetchrow_array) {
	#	my $anaLink;
	#	foreach my $kw ($name,$dFile,$wFile,$lCode) {
	#		if ($anaTargets{$kw}) {
	#			$anaLink=$kw;
	#			last;
	#		}
	#	}
	#	$analysisIDs{$anaLink}=$anaID;
	#}
	#$sthAna->finish;

	###>Fetching all Analyses in Experiment<###
	my %expAnalyses;
	my $sthExpAna=$dbh->prepare("SELECT A.ID_ANALYSIS,A.NAME,DATA_FILE,WIFF_FILE,LAB_CODE FROM ANALYSIS A
								INNER JOIN SAMPLE S ON S.ID_SAMPLE=A.ID_SAMPLE WHERE S.ID_EXPERIMENT=?");
	$sthExpAna->execute($experimentID);
	while (my ($anaID,$name,$dFile,$wFile,$lCode)=$sthExpAna->fetchrow_array) {
		$expAnalyses{$name}=$anaID;
		$expAnalyses{$lCode}=$anaID if $lCode;
		foreach my $file ($dFile,$wFile) {
			$expAnalyses{$file}=$anaID;
			$file=~s/\.[^\.]+$//; # remove file extension
			$expAnalyses{$file}=$anaID;
		}
	}
	$sthExpAna->finish;
	
	####>Matching file links to Analyses in Experiment<####
	my %analysisIDs;
	foreach my $anaLink (keys %anaTargets) {
		my $modLink=$anaLink;
		$modLink=~s/\.[^\.]+$//; # remove file extension (if any)
		foreach my $link ($anaLink,$modLink) {
			if ($expAnalyses{$link}) {
				$analysisIDs{$anaLink}=$expAnalyses{$link};
				last;
			}
		}
	}
	
	###>Checking if everything was matched<###
	my $numMapAna=scalar keys %analysisIDs;
	my @unmappedLinks;
	if ($numMapAna) {
		if ($numMapAna < $numAllAna) {
			my $linkStrg=($numAllAna-$numMapAna == 1)? 'Analysis link was' : ($numAllAna-$numMapAna).' Analysis links were';
			print "&bull; <FONT class=\"title3\">The following $linkStrg not found in selected Experiment:<BR>\n";
			foreach my $anaLink (sort{&promsMod::sortSmart($a,$b)} keys %anaTargets) {
				next if $analysisIDs{$anaLink};
				print "&nbsp;&nbsp;&nbsp;-$anaLink<BR>\n";
				push @unmappedLinks,$anaLink;
			}
			print "</FONT><BR>\n";
			
			##>Checking impact on list of usable BioSamp<##
			my @lostBioSamps;
			foreach my $anaLink (@unmappedLinks) {
				foreach my $bsName (keys %bioSampleIDs) {
					my @unmapIdx;
					foreach my $l (0..$#{$bioSampleLinks{$bsName}}) {
						push @unmapIdx,$l if $bioSampleLinks{$bsName}[$l][0] eq $anaLink;
					}
					if (scalar @unmapIdx) { # at least 1
						if (scalar @unmapIdx == scalar @{$bioSampleLinks{$bsName}}) { # no more links for this BS
							delete $bioSampleIDs{$bsName};
							push @lostBioSamps,$bsName;
						}
						else {
							foreach my $l (@unmapIdx) {
								splice @{$bioSampleLinks{$bsName}},$l,1;
							}
						}
					}
				}
			}
			my $numBsLost=scalar @lostBioSamps;
			if ($numBsLost) {
				my $lostStrg=($numBsLost==1)? 'Biological sample' : "$numBsLost Biological samples";
				print "&bull; <FONT class=\"title3\">The following $lostStrg will not be linked due to unmatched Analysis link(s):<BR>\n";
				foreach my $bsName (sort{&promsMod::sortSmart($a,$b)} @lostBioSamps) {
					print "&nbsp;&nbsp;&nbsp;-$bsName<BR>\n";
				}
			}
			print "</FONT><BR>\n";
		}
	}
	else {
		$dbh->disconnect;
		&exitOnError('Analysis links from file were not found in selected Experiment');
	}
	
	####>Fetching (existing/potential observations linked to matched Analyses<####
	my (@anaList,%anaListIdx,%hierarchy,%observation,%position);
	foreach my $anaID (values %analysisIDs) {
		push @anaList,[$anaID];
		$anaListIdx{$anaID}=$#anaList;
	}
	&fetchObservations($dbh,\@anaList,\%hierarchy,\%observation,\%position); # @anaList is updated

	####>Recording links in DB<####
	my $sthBsName=$dbh->prepare("SELECT NAME FROM BIOSAMPLE WHERE ID_BIOSAMPLE=?");
	my $sthUpSamp=$dbh->prepare("UPDATE OBSERVATION SET ID_BIOSAMPLE=? WHERE ID_OBSERVATION=?");
	my $sthAddObs=$dbh->prepare("INSERT INTO OBSERVATION (ID_ANALYSIS,ID_BIOSAMPLE,TARGET_POS) VALUES (?,?,?)");
	my $sthAddOMod=$dbh->prepare("INSERT INTO OBS_MODIFICATION (ID_OBSERVATION,ID_MODIFICATION) VALUES (?,?)");
	my (%bioSampNames,%bsLinked);
	my $numNewLinks=0;
	
	foreach my $bsName (sort{$bioSamplePos{$a}<=>$bioSamplePos{$b}} keys %bioSampleIDs) {
		foreach my $refLink (@{$bioSampleLinks{$bsName}}) {
			my ($anaLink,$linkLabelCode)=@{$refLink};
			my $anaID=$analysisIDs{$anaLink};
			my ($labelType,$refTargets)=@{$anaList[$anaListIdx{$anaID}]}[-2,-1]; # <- from &fetchObservations
			
			###>Matching observation<###
			my $obsCode;
			if ($labelType=~/FREE|NONE/) {
				$obsCode="$anaID,0";
			}
			else {
				my @allObsCodes=sort{$position{$a}<=>$position{$b}} values %{$refTargets};
#print "==(",join('/',@allObsCodes),")<BR>\n";
				my $labelingPlex=@allObsCodes;
				if ($linkLabelCode=~/^\d+$/ && $linkLabelCode < 100) {$obsCode=$allObsCodes[$linkLabelCode-1];} # 100 < iTRAQ/TMT tags, SILAC: 1,2,3max
				elsif ($labelType=~/SILAC/i) {
#print "**$linkLabelCode<BR>\n";
					if ($linkLabelCode=~/^L(ight)*/i) {$obsCode=$allObsCodes[0];}
					elsif ($linkLabelCode=~/^H(eavy)*/i) {$obsCode=$allObsCodes[-1];}
					elsif ($linkLabelCode=~/M(edium)*/i) {$obsCode=($labelingPlex==3)? $allObsCodes[1] : ($allObsCodes[0] eq "$anaID,-1")? $allObsCodes[1] : $allObsCodes[0];} # Ligth/Medium/Heavy : Ligth/Medium : Medium/Heavy
				}
				else { # assume iTRAQ/TMT with tag names
					$obsCode=$refTargets->{$linkLabelCode}; # eg. 118 (tag name) => 1 (tgtPos)
				}
			}
			my $targetObsStrg=$anaLink;
			$targetObsStrg.=' / '.$linkLabelCode if $linkLabelCode;
			unless ($obsCode) {
				print "<BR>&nbsp;&nbsp;&nbsp;<B>-Warning: BioSample '$bsName' could not be linked to any Observation matching '$targetObsStrg'.<B>\n";
				next;
			}
			
			###>Recording link in DB<###
			##>Observation exists already
			if ($observation{$obsCode}) {
#print "-BS=$bsName, OBS=$obsCode<BR>\n";
				my ($obsID,$bsID)=@{$observation{$obsCode}};
				if ($bsID) {
					if ($bsID != $bioSampleIDs{$bsName}) {
						unless ($bioSampNames{$bsID}) {
							$sthBsName->execute($bsID);
							($bioSampNames{$bsID})=$sthBsName->fetchrow_array;
						}
						print "<BR>&nbsp;&nbsp;&nbsp;<B>-Warning: BioSample '$bsName' could not be linked to '$targetObsStrg'. This Observation is already linked to BioSample '$bioSampNames{$bsID}'";
					}
					else {
						# already linked: nothing to do
					}
				}
				else { # Update to link BS
					$sthUpSamp->execute($bioSampleIDs{$bsName},$obsID);
					$numNewLinks++;
					$bsLinked{$bsName}=1;
				}
			}
			##>Observation must be created
			else {
				my ($anaID,$targetPos,@mods)=split(',',$obsCode);
#print "+BS=$bsName, OBS=$obsCode<BR>\n";
				$sthAddObs->execute($anaID,$biosampleID,$targetPos);
				my $newObsID=$dbh->last_insert_id(undef,undef,'OBSERVATION','ID_OBSERVATION');
				foreach my $modID (@mods) { # can be empty
					$sthAddOMod->execute($newObsID,$modID);
				}
				$numNewLinks++;
				$bsLinked{$bsName}=1;
			}
		}
	}
	$sthBsName->finish;
	$sthUpSamp->finish;
	$sthAddObs->finish;
	$sthAddOMod->finish;

#$dbh->rollback;
	$dbh->commit;
	$dbh->disconnect;
	
	my $numBioSamp=scalar keys %bsLinked;
	my $bsStrg=($numBioSamp==1)? '' : 's';
	my $obsStrg=($numNewLinks==1)? '' : 's';
	print qq
|<BR><BR>
<FONT class="title2">$numBioSamp BioSample$bsStrg successfully linked to $numNewLinks Observation$obsStrg.</FONT>
<BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="button" class="title2" value=" Continue " onclick="parent.itemFrame.location.reload(); parent.optionFrame.selectOption();"/>
</BR></BR>
</BODY>
</HTML>
|;
	exit;
}

my (@experiments,%bioSampleNames,@sortedBioSamples,$title);
my $sthExp=$dbh->prepare("SELECT ID_EXPERIMENT,NAME
						  FROM EXPERIMENT
						  WHERE ID_PROJECT=$projectID
						  AND EXPERIMENT.ID_EXPERIMENT NOT IN (
							  SELECT E2.ID_EXPERIMENT
							  FROM EXPERIMENT E2
							  INNER JOIN USER_EXPERIMENT_LOCK EU ON EU.ID_EXPERIMENT=E2.ID_EXPERIMENT
							  WHERE E2.ID_PROJECT=$projectID AND EU.ID_USER='$userID'
						  )
						  ORDER BY DISPLAY_POS");
$sthExp->execute;
while (my ($expID,$expName)=$sthExp->fetchrow_array) {
	push @experiments,[$expID,$expName];
}
$sthExp->finish;

if ($biosampleID) { # called from a bioSample
	my ($bioSampleName)=$dbh->selectrow_array("SELECT NAME FROM BIOSAMPLE WHERE ID_BIOSAMPLE=$biosampleID");
	$title="Linking <FONT color=\"#DD0000\">$bioSampleName</FONT> to Observations";
}
else { # called to top branch
	#my $sthBS=$dbh->prepare("SELECT BS.ID_BIOSAMPLE,NAME,GROUP_CONCAT(ID_OBSERVATION SEPARATOR ':') FROM BIOSAMPLE BS,PROJECT_BIOSAMPLE PB,OBSERVATION O WHERE BS.ID_BIOSAMPLE=PB.ID_BIOSAMPLE AND BS.ID_BIOSAMPLE=O.ID_BIOSAMPLE AND PB.ID_PROJECT=$projectID GROUP BY BS.ID_BIOSAMPLE");
	my $sthBS=$dbh->prepare("SELECT BS.ID_BIOSAMPLE,NAME FROM BIOSAMPLE BS,PROJECT_BIOSAMPLE PB WHERE BS.ID_BIOSAMPLE=PB.ID_BIOSAMPLE AND PB.ID_PROJECT=$projectID");
	$sthBS->execute;
	while (my ($bsID,$bsName)=$sthBS->fetchrow_array) {
		$bioSampleNames{$bsID}=$bsName;
	}
	$sthBS->finish;
	$title="Linking Biological Samples to Observations";
	@sortedBioSamples=(sort{&promsMod::sortSmart(lc($bioSampleNames{$a}),lc($bioSampleNames{$b}))} keys %bioSampleNames);
}


#######################
####>Starting HTML<####
#######################
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Biological Sample Observation Link Management</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT type="text/javascript">
function extendSelection(obsIndex,chkStatus) {
	if (!document.obsForm.autoExtend.checked) return;
	var obsBoxList=document.obsForm.selObs;
	if (!obsBoxList.length) return;
	var selObsFullData=obsBoxList[obsIndex].value.split(':'); // "anaID,targetPos[,mod1,...,modN][:obsID]"
	var selObsData=selObsFullData[0].split(','); // "anaID,targetPos[,mod1,...,modN]"
	var selTargetPos=selObsData[1],selModStrg='';
	if (selTargetPos==-1) { // isotope labeling: SILAC
		var modData=selObsData.slice(2); // empty for No label channel "anaID,targetPos")
		selModStrg=modData.join(',');
	}
	for (var i=obsIndex+1; i<obsBoxList.length; i++) {
		if (selTargetPos==0) { // label free
			obsBoxList[i].checked=chkStatus;
		}
		else {
			var obsFullData=obsBoxList[i].value.split(':');
			var obsData=obsFullData[0].split(',');
			if (selTargetPos==-1) { // isotope labeling: SILAC
				if (obsData[1]==-1) { // isotope labeling
					var modData=obsData.slice(2);
					modStrg=modData.join(',');
					if (modStrg==selModStrg) {
						obsBoxList[i].checked=chkStatus;
					}
				}
			}
			else { // isobaric: iTRAQ (targetPos carries the channel info)
				if (obsData[1]==selTargetPos) {
					obsBoxList[i].checked=chkStatus;
				}
			}
		}
	}
}
function checkForm(myForm) {
	if (!myForm.targetExp.value) {
		alert('ERROR: Select a target Experiment.');
		return false;
	}
	return true;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">$title</FONT><BR><BR><BR><BR>
|;
my $selectStrg='';
if ($experimentID && scalar @experiments > 1 && !$biosampleID) {
	
}
elsif (!$experimentID) {
	$selectStrg='Select an';
}
print qq
|<FONT class="title2">$selectStrg Experiment:</FONT><SELECT name="expID" class="title2" onchange="window.location='$promsPath{cgi}/linkBioSample2Observations.cgi?biosampleID=$biosampleID&projectID=$projectID&expID='+this.value"><OPTION value="">-= Select =-</OPTION>
|;
foreach my $expData (@experiments) {
	print "<OPTION value=\"$expData->[0]\"";
	print ' selected' if $expData->[0]==$experimentID;
	print ">$expData->[1]</OPTION>\n";
}
print "</SELECT><BR>\n";

if ($experimentID) {
	my ($numObsLinked)=$dbh->selectrow_array("SELECT COUNT(O.ID_OBSERVATION) FROM OBSERVATION O
												INNER JOIN ANALYSIS A ON A.ID_ANALYSIS=O.ID_ANALYSIS
												INNER JOIN SAMPLE S ON S.ID_SAMPLE=A.ID_SAMPLE
												WHERE S.ID_EXPERIMENT=$experimentID AND O.ID_BIOSAMPLE IS NOT NULL");
	my ($numObsNotLinked)=$dbh->selectrow_array("SELECT COUNT(O.ID_OBSERVATION) FROM OBSERVATION O
												INNER JOIN ANALYSIS A ON A.ID_ANALYSIS=O.ID_ANALYSIS
												INNER JOIN SAMPLE S ON S.ID_SAMPLE=A.ID_SAMPLE
												WHERE S.ID_EXPERIMENT=$experimentID AND O.ID_BIOSAMPLE IS NULL");
	my $numAllObs=$numObsLinked+$numObsNotLinked;
	print qq|<FONT class="title3">($numObsNotLinked free/$numAllObs Observations)</FONT><BR><BR>|;
}
else {
	$dbh->disconnect;
	print "</BODY>\n</HTML>\n";
	exit;
}

unless ($biosampleID) {
	####>Import from file<####
	print qq
|<SCRIPT type="text/javascript">
function checkForm(myForm) {
	if (!myForm.dataFile.value) {
		alert('ERROR: No data file selected');
		return false;
	}
	return true;
}
</SCRIPT>
<BR>
<FORM name="importForm" method="post" enctype="multipart/form-data" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="ACT" value="import">
<INPUT type="hidden" name="projectID" value="$projectID">
<INPUT type="hidden" name="expID" value="$experimentID">
<TABLE cellspacing=0>
<TR bgcolor="$darkColor"><TH align="left" valign="top" nowrap colspan=2>&nbsp;<FONT class="title2">Import links from file<SUP>*</SUP>:</FONT></TH>
	<TD nowrap><INPUT type="file" name="dataFile"><BR>
<FONT class="font11"><B>Format:</B> tab-separated list (with <b>header</B>) of <B>Biological samples</B> and <B>Analysis names</B>/<B>Search result files</B>/<B>Lab. codes</B>.<BR>
For <B>labeled</B> samples, an <B>additional</B> column is required containing:<BR>
&nbsp;&nbsp;&nbsp;&nbsp;-<B>Chanel rank</B> (1, .., n) or<BR>
&nbsp;&nbsp;&nbsp;&nbsp;-<B>Tag name</B> (113, .., 121) for <B>iTRAQ/TMT</B> or<BR>
&nbsp;&nbsp;&nbsp;&nbsp;-<B>Label state</B> (None, Medium, Heavy) for <B>SILAC</B></FONT>
</TD><TH valign="middle">&nbsp;<INPUT type="submit" class="title3" value=" Proceed ">&nbsp;</TH></TR>
<TR><TD colspan=3>&nbsp;&nbsp;&nbsp;<SUP>*</SUP><I>Pre-existing links are not overwritten</I></TD></TR>
</TABLE>
</FORM>
|;

	####>Copy from another Experiment<####
	if (scalar @experiments > 1) {
		print qq
|<FONT class="title2">Or</FONT><BR><BR>
<FORM name="expForm" method="POST" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="ACT" value="propagate">
<INPUT type="hidden" name="projectID" value="$projectID">
<INPUT type="hidden" name="expID" value="$experimentID">
<TABLE cellspacing=0>
<TR bgcolor="$darkColor"><TD colspan=3></TD></TR>
<TR bgcolor="$darkColor"><TD align=right nowrap>&nbsp;<FONT class="title2">Propagate links from Experiment above to this one<SUP>*</SUP>:</FONT></TD><TD><SELECT name="targetExp" class="title2"><OPTION value="">-= Select =-</OPTION>
|;
		foreach my $expData (@experiments) {
			print "<OPTION value=\"$expData->[0]\"";
			print ' disabled' if $expData->[0]==$experimentID;
			print ">$expData->[1]</OPTION>\n";
		}
		print qq
|</SELECT>&nbsp;</TD><TH rowspan=2 valign="middle">&nbsp;&nbsp;<INPUT class="title3" type="submit" value=" Proceed ">&nbsp;&nbsp;</TH></TR>
<TR bgcolor="$darkColor"><TD align=right nowrap><FONT class="title2">based on matching:</FONT></TD><TD nowrap><SELECT name="matchType" class="title2"><OPTION value="analysis">Analyses name</OPTION><OPTION value="datafile">Analyses MS data file</OPTION><OPTION value="design">Design states</OPTION></SELECT></TD></TR>
<TR bgcolor="$darkColor"><TD colspan=3></TD></TR>
<TR><TD colspan=3>&nbsp;&nbsp;&nbsp;<SUP>*</SUP><I>Pre-existing links are not overwritten</I></TD></TR>
</TABLE>
</FORM>
<FONT class="title2">Or</FONT><BR><BR><FONT class="title2">Generate links manually:</FONT><BR>
|;
	}
}

my @sthGetAnaInfo=( # 2d-gel or sample
	$dbh->prepare("SELECT ID_ANALYSIS,'gel2d',G.NAME,'spot',SP.NAME,'analysis',A.NAME FROM ANALYSIS A,SAMPLE S,SPOT SP,GEL2D G WHERE G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT AND S.ID_SAMPLE=A.ID_SAMPLE AND VALID_STATUS > 0 AND G.ID_EXPERIMENT=? ORDER BY G.DISPLAY_POS,SP.NAME,S.DISPLAY_POS,A.DISPLAY_POS"),
	$dbh->prepare("SELECT ID_ANALYSIS,'sample',S.NAME,'analysis',A.NAME FROM ANALYSIS A,SAMPLE S WHERE S.ID_SAMPLE=A.ID_SAMPLE AND ID_SPOT IS NULL AND VALID_STATUS > 0 AND S.ID_EXPERIMENT=? ORDER BY S.DISPLAY_POS,A.DISPLAY_POS")
);

###> Get Analysis information in order to fill OBSERVATION, OBS_EXPCONDITION & OBS_MODIFICATION tables
my (%hierarchy,%position,%observation,%labelMods);
my $obsRank=0;
foreach my $sthGAI (@sthGetAnaInfo) {
my @anaList;
	$sthGAI->execute($experimentID);
	while (my ($anaID,@hierarchy) = $sthGAI->fetchrow_array) {
		push @anaList,[$anaID,@hierarchy];
	}
	&fetchObservations($dbh,\@anaList,\%hierarchy,\%observation,\%position,\$obsRank);
	$sthGAI->finish;
}

$dbh->disconnect;

my $maxHierarchyIndex=0;
foreach my $refH (values %hierarchy) {$maxHierarchyIndex=$#{$refH} if $maxHierarchyIndex < $#{$refH};}
my $obsColSpan=($maxHierarchyIndex+1)/2;
my $totColSpan=$obsColSpan+1;
my $headerText=($biosampleID)? 'Link' : 'Bio-Sample';
print qq
|<FORM name="obsForm" method="POST">
<INPUT type="hidden" name="ACT" value="store">
<INPUT type="hidden" name="biosampleID" value="$biosampleID">
<INPUT type="hidden" name="expID" value="$experimentID">
<INPUT type="hidden" name="projectID" value="$projectID">
|;
my @prevHierarchy=();
my $obsIndex=-1;
foreach my $obsCode (sort{$position{$a}<=>$position{$b}} keys %position){
	my $checkStrg;
	if ($biosampleID) { # called from a bioSample
		next if ($observation{$obsCode} && $observation{$obsCode}[1] && $observation{$obsCode}[1] != $biosampleID); # observation linked to another bioSample
		$checkStrg=($observation{$obsCode} && $observation{$obsCode}[1] && $observation{$obsCode}[1] == $biosampleID)? ' checked' : '';
	}
	$obsIndex++;
	if ($obsIndex==0) {
		if ($biosampleID) {
			print qq
|<INPUT type="checkbox" name="autoExtend" value="1" checked><FONT class="title3">Auto-extend selection</FONT>
|;
		}
		print qq
|<TABLE bgcolor="$darkColor" cellspacing=0>
<TR><TD class="title2 rbBorder" colspan="$obsColSpan" style="padding:2px">&nbsp;Available Observations</TD><TD class="title2 bBorder" style="padding:2px">&nbsp;$headerText&nbsp;</TD></TR>
|;
	}
	print "<TR class=\"list\" bgcolor=\"$lightColor\">";
	for (my $i=0;$i<$#{$hierarchy{$obsCode}}; $i+=2) {
		print "<TD valign='bottom' nowrap>";
		my $hideItem=1;
		for (my $j=0; $j<=$i; $j+=2) {
			if (!$prevHierarchy[$j] || $prevHierarchy[$j+1] ne $hierarchy{$obsCode}[$j+1]) {
				$hideItem=0;
				last;
			}
		}
		unless ($hideItem) {
			print "&nbsp;<B>".$hierarchy{$obsCode}[$i+1];
			print "&nbsp;>" if $i< $#{$hierarchy{$obsCode}}-1;
			print "</B>";
		}
		print "</TD>";
		if ($i==$#{$hierarchy{$obsCode}}-1) {
			for (my $j=$i+2; $j<$maxHierarchyIndex; $j+=2) {print "<TD></TD>";}
			my $fullObsCode=($observation{$obsCode})? "$obsCode:$observation{$obsCode}[0]" : $obsCode;
			print "<TH valign=\"bottom\" nowrap>&nbsp;";
			if ($biosampleID) { # called from a bioSample
				print "<INPUT type=\"checkbox\" name=\"selObs\" value=\"$fullObsCode\" onchange=\"extendSelection($obsIndex,this.checked)\" $checkStrg><INPUT type=\"hidden\" name=\"obsList\" value=\"$fullObsCode\">";
			}
			else { # called to top branch
				print "<SELECT name=\"bioSamples\"><OPTION value=\"0\">* Not linked *</OPTION>\n";
				my $selBsID=($observation{$obsCode} && $observation{$obsCode}[1])? $observation{$obsCode}[1] : 0;
				foreach my $bsID (@sortedBioSamples) {
					print "<OPTION value=\"$bsID\"";
					print ' selected' if $bsID==$selBsID;
					print ">$bioSampleNames{$bsID}</OPTION>\n";
				}
				print "</SELECT><INPUT type=\"hidden\" name=\"obsList\" value=\"$fullObsCode=$selBsID\">";
			}
			print "&nbsp;</TH>";
		}
	}
	@prevHierarchy=@{$hierarchy{$obsCode}};
	print "</TR>\n<TR bgcolor=\"$darkColor\"><TD colspan=$totColSpan></TD></TR>\n";
}
if ($obsIndex>=0) {
	print qq
|<TR bgcolor="$darkColor"><TD colspan=$totColSpan></TD></TR>
<TR bgcolor="$darkColor"><TH colspan=$totColSpan><INPUT type="submit" class="title3" value=" Save "></TH></TR>
<TR bgcolor="$darkColor"><TD colspan=$totColSpan></TD></TR>
</TABLE>
|;
}
else {
	print "<FONT class=\"title2\">No observations available.</FONT>\n";
}
print qq
|<FORM>
<BR><BR>
</BODY>
</HTML>
|;


sub exitOnError {
	my ($errorStrg)=@_;
	return unless $errorStrg;
	print qq
|<BR><BR><FONT class="title2" color="#DD0000">ERROR: $errorStrg!</FONT>
<BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="button" class="title2" value=" OK " onclick="parent.optionFrame.selectOption();"/>
</BR></BR>
</BODY>
</HTML>
|;
	exit;
}

sub fetchObservations {
	my ($dbh,$refAnaList,$refHierarchy,$refObservation,$refPosition,$obsRank)=@_;
	$refPosition={} unless $refPosition; # optional
	$$obsRank=0 unless $obsRank;
	
	my $sthLQ=$dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND FOCUS='peptide' AND ID_ANALYSIS=? LIMIT 1");
	my $sthLA=$dbh->prepare("SELECT DISTINCT M.ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES FROM MODIFICATION M,ANALYSIS_MODIFICATION AM WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND IS_LABEL=1 AND ID_ANALYSIS=? ORDER BY M.ID_MODIFICATION");
	my $sthObs=$dbh->prepare("SELECT O.ID_OBSERVATION,O.ID_BIOSAMPLE,O.TARGET_POS,GROUP_CONCAT(ID_MODIFICATION ORDER BY ID_MODIFICATION SEPARATOR ',')
								FROM OBSERVATION O
								LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION
								WHERE O.ID_ANALYSIS=? GROUP BY O.ID_OBSERVATION");
	my %labelMods;
	foreach my $refAnaData (@{$refAnaList}) {
		my ($anaID,@hierarchy)=@{$refAnaData};
		push @{$refAnaData},('FREE',{}); # default
					
		##>Fetch all possible quantif channels
		$sthLQ->execute($anaID);
		my ($quantifAnnot)=$sthLQ->fetchrow_array;

		##>A peptide quantification associated with Analysis
		# IMPORTANT NOTE (31/01/14 PP): TARGET_POS has meaning only for iTRAQ-like quantif
		# Otherwise: 0 for label-free and -1 for all channel pos of SILAC-like (true match given by OBS_MODIFCATION)
		if ($quantifAnnot) {
			$quantifAnnot=~s/::SOFTWARE=[^:]+::/::/; # remove software info for back compatibility
			my ($labelTypeStrg,@labelInfo)=split('::',$quantifAnnot);
			my ($labelType)=($labelTypeStrg =~ /^LABEL=(.+)/);
			if ($labelType =~ /FREE|NONE/) { # no labeling in analysis
				my $obsCode="$anaID,0";
				@{$refHierarchy->{$obsCode}}=(@hierarchy,'label','Label-free');
				$refPosition->{$obsCode}=++$$obsRank;
				$refAnaData->[-2]=$labelType;
				#$refAnaData->[-1]->{$chanName}=$obsCode;
			}
			else {
				my $labelModID;
				foreach my $infoStrg (@labelInfo) {
					my ($targetPos,$chanName,$labelStrg)=split(';',$infoStrg);
					my $obsCode; #="$anaID,$targetPos";
					my $targetName;
					if ($labelType=~/SILAC/i) {
						$obsCode="$anaID,-1";
						$targetName="$chanName ";
						$targetName.='[No label]' unless $labelStrg; # poorly annotated unlabeled channel
						my $first=1;
						#my @mods;
						my %mods; # hash because same mod can be used on different aa (mod is repeated)
						foreach my $modStrg (split ('@',$labelStrg)) { # multiple mods for same label channel
							my @labelData=split('#',$modStrg);
							if ($first) {$first=0;}
							else {$targetName.='+';}
							if ($labelData[1] =~ /No label/i) {
								$targetName.="[$labelData[1]]";
							}
							else {
								$targetName.="[$labelData[1]&raquo;$labelData[2]]"; #&raquo; &middot; &rarr;
								my $modID=($labelData[4])? $labelData[4] : ($labelMods{$labelData[1]})? $labelMods{$labelData[1]} : &promsMod::getModificationIDfromString($dbh,$labelData[1]); # $labelData[4] if MassChroQ XIC
								$modID=0 unless $modID;
								$labelMods{$labelData[1]}=$modID;
								#$obsCode.=",$modID";
								#push @mods,$modID;
								$mods{$modID}=1;
							}
						}
						if (scalar keys %mods) { # label channel
							$obsCode.=','.join(',',sort{$a<=>$b} keys %mods); # order is needed for later compare with %expCondition
						}
					}
					elsif ($labelType=~/itraq|TMT/i) {
						$obsCode="$anaID,$targetPos";
						$targetName=$chanName;
						unless ($labelModID) { # only once
							$sthLA->execute($anaID);
							($labelModID)=$sthLA->fetchrow_array;
						}
						$obsCode.=",$labelModID" if $labelModID;
					}
					@{$refHierarchy->{$obsCode}}=(@hierarchy,'label',$labelType,'target',$targetName);
					$refPosition->{$obsCode}=++$$obsRank;
					$refAnaData->[-2]=$labelType;
					$refAnaData->[-1]{$chanName}=$obsCode;
				}
			}
		}
		##>No peptide quantification associated with Analysis => guess from ANALYSIS_MODIFICATION table ****WARNING: NOT reliable
		else {
			$sthLA->execute($anaID);
			my $isLabeled=0;
			my $obsCode;
			my $channelPos=0;
			while (my ($modID,$psiName,$intName,$synName)=$sthLA->fetchrow_array) {
				$psiName='' unless $psiName;
				$intName='' unless $intName;
				$synName='' unless $synName;
				my $labelName=$psiName.'##'.$intName.'##'.$synName;
				if ($labelName=~/itraq.+plex/i) {
					$isLabeled=1;
					my %targetInfo;
					my $targetPos=0;
					my %infoSubmit=&promsMod::getSearchParam($dbh,$anaID);
					if ($infoSubmit{'q:Quantification'}) {
						my $first=1;
						foreach my $targetName (split(':',$infoSubmit{'q:Quantification'})) {
							if ($first) { # skip quantif Name
								$first=0;
								next;
							}
							$targetInfo{++$targetPos}=$targetName;
						}
					}
					else { # backup solution
						my ($plexNum)=($labelName=~/(\d+)\s*plex/i);
						foreach my $targetPos (1..$plexNum) {
							$targetInfo{$targetPos}=112+$targetPos;
							$targetInfo{$targetPos}=121 if $targetInfo{$targetPos}==120; # replaces 120 with 121 (120 does not exist...?)
						}
					}
					foreach my $targetPos (sort{$a<=>$b} keys %targetInfo) {
						$obsCode="$anaID,$targetPos,$modID";
						@{$refHierarchy->{$obsCode}}=(@hierarchy,'label',$labelName,'target',$targetInfo{$targetPos});
						$refPosition->{$obsCode}=++$$obsRank;
						$refAnaData->[-2]=$labelName;
						$refAnaData->[-1]->{$targetInfo{$targetPos}}=$obsCode;
					}
					next;
				}
				unless ($isLabeled) { # 1st loop: create a non-labeled channel
					#$channelPos++;
					$obsCode="$anaID,0";
					@{$refHierarchy->{$obsCode}}=(@hierarchy,'label','No label');
					$refPosition->{$obsCode}=++$$obsRank;
					$isLabeled=1;
					$refAnaData->[-2]='FREE';
					#$refAnaData->[-1]->{''}=$obsCode;
				}
				# !!!!!!!!!! TODO: Not compatible with single state defined by multiple modifications !!!!!!! <= SILAC
				$channelPos++;
				$obsCode="$anaID,$channelPos,$modID";
				@{$refHierarchy->{$obsCode}}=(@hierarchy,'label','SILAC','target',$labelName);
				$refPosition->{$obsCode}=++$$obsRank;
				$refAnaData->[-2]='SILAC';
				$refAnaData->[-1]->{$labelName}=$obsCode;

			}
			unless ($isLabeled) { # no labeling in analysis
				$obsCode="$anaID,0";
				@{$refHierarchy->{$obsCode}}=(@hierarchy,'label','Label-free');
				$refPosition->{$obsCode}=++$$obsRank;
				$refAnaData->[-2]='FREE';
				#$refAnaData->[-1]->{''}=$obsCode;
			}
		}

		# Get whether or not this analysis is already associated to an OBSERVATION & BIOSAMPLE!
		$sthObs->execute($anaID);
		while (my ($obsID,$bioSampID,$targetPos,$modCode)=$sthObs->fetchrow_array) {
			$targetPos=-1 if ($quantifAnnot && $quantifAnnot=~/SILAC/); # back compatibility for Design before 31/01/14
			my $obsCode="$anaID,$targetPos";
			$obsCode.=",$modCode" if $modCode;
			@{$refObservation->{$obsCode}}=($obsID,$bioSampID);
		}
	}
	$sthLA->finish;
	$sthLQ->finish;
	$sthObs->finish;
}

####>Revision history<####
# 1.2.0 [FEATURE] Added option to use link information from uploaded file (PP 16/03/20)
# 1.1.1 [FEATURE] Remove locked experiments from BioSamples searches (VS 08/08/19)
# 1.1.0 Added option to propagate links from one Exp to another based on analysis name/MS file or design states name (PP 25/04/19)
# 1.0.6 Also removes software version from QUANTIF_ANNOT (PP 09/12/18)
# 1.0.5 Minor modification for TMT (GA 03/04/17)
# 1.0.4 Compatible with SILAC cases where same labeling mod is used on different aa (PP 21/07/15)
# 1.0.3 Update for SOFTWARE tag in QUANTIF_ANNOT field (PP 28/04/15)
# 1.0.2 Auto-extend observation selection (PP 15/05/15)
# 1.0.1 Handles multi-sample Observation link (PP 28/08/17)
# 1.0.0 New script to link a biological sample to selected observerations (PP 17/07/14)
