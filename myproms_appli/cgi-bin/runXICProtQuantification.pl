#!/usr/local/bin/perl -w

################################################################################
# runXICProtQuantification.pl       2.7.0                                      #
# Component of site myProMS Web Server                                         #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
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
use promsConfig;
use promsMod;
use promsQuantif;
use strict;
#use XML::SAX;
use MIME::Base64;
use POSIX qw(strftime);
#use File::Copy qw(copy move);
use File::Copy::Recursive qw(dirmove dircopy);
#use File::Path qw(rmtree); # remove_tree

#exit;
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo('no_user');
#my $minProtNumber=30; # used to desactivate FDR option

###############################
####>Recovering parameters<####
###############################
my ($quantifID,$quantifDate,$quantItemID)=@ARGV; # only $quantifID & $quantifDate defined if Design-quantif
my $quantifDir="$promsPath{tmp}/quantification/$quantifDate";
my %quantifParameters=&promsQuantif::getQuantificationParameters("$quantifDir/quantif_info.txt");
my $Rdesign=$quantifParameters{'DB'}{'ALGO_TYPE'}[0]; # <=> {R}{design} for v3 only!
my ($call,$analysisID,$parentQuantifID); #,$ratioType
if ($quantItemID) { # no longer possible
	$call='ANALYSIS';
	($analysisID,$parentQuantifID)=split(/\./,$quantItemID);
	#$ratioType='Ratio';
}
else {
	$call='DESIGN';
	#$ratioType=$quantifParameters{'DB'}{'RATIO_TYPE'}[0]; # Ratio, SuperRatio or SimpleRatio;
}
my $ratioType=$quantifParameters{'DB'}{'RATIO_TYPE'}[0]; # Ratio, SuperRatio or SimpleRatio;
my $topN=($quantifParameters{'DB'}{'TOP_N'})? $quantifParameters{'DB'}{'TOP_N'}[0] : 0; # number peptide used for label-free (except old algo 'Ratio')
my $matchingPep=($quantifParameters{'DB'}{'MATCH_PEP'})? $quantifParameters{'DB'}{'MATCH_PEP'}[0] : 0; # label-free SimpleRatio: peptides must be the same across conditions

my $dbh=&promsConfig::dbConnect('no_user');
$dbh->do("UPDATE QUANTIFICATION SET STATUS=0 WHERE ID_QUANTIFICATION=$quantifID");
$dbh->commit;

my $fileStat="$quantifDir/status_$quantifID.out"; # : "$quantifDir/status_$quantItemID.out";
open(FILESTAT,">$fileStat") || die "Error while opening $fileStat";
print FILESTAT "Started ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
print FILESTAT "1/3 Generating data files\n";
close FILESTAT;

my $projectID=&promsMod::getProjectID($dbh,$quantifID,'QUANTIFICATION');

my ($quantiAnnotation,$quantifiedModifID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_MODIFICATION FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID");
my ($labeling)=($quantiAnnotation=~/LABEL=([^:]+)/);
$labeling=uc($labeling);

####> R part : statistical analysis of the proteins values
####> 1st: Create the directories for R analysis (via promsMod function)
#my $runDir=($call eq 'DESIGN')? "$quantifDir/quanti_$quantifID" : "$quantifDir/quanti_$quantItemID"."_$quantifID";
my $runDir="$quantifDir/quanti_$quantifID";
my $dataDir="$runDir/data";
my $resultDir="$runDir/results";
my $graphDir="$resultDir/graph";
#make_path($dataDir,$graphDir,{verbose=>0,mode=>0775}); # $runDir,$resultDir will be created automatically
mkdir $runDir unless -e $runDir; # already created if runXICProtQuantification.pl launched on cluster
mkdir $dataDir;
mkdir $resultDir;
mkdir $graphDir;

#open (DEBUG,">$promsPath{tmp}/quantification/debug.txt"); # DEBUG!!!!
open(DATA,">$dataDir/table.txt"); # Name of data table
if ($ratioType eq 'Ratio') { # Old algos (FC)
	print DATA "Protein_ID\tPep_Seq\tVar_Mod\tCharge\tPeptide_IDs\tProtein_Name\tProtein_Validity";
}
else { # Algos v2+
	print DATA "Protein_ID\tPeptide\tSample\tTreatment\tReplicate\tTechnical_Replicate\tExperiment\tQuantif_Set\tPeptide_IDs\tProtein_Name\tProtein_Validity\tValue";
}

################################
####>Get states composition<####
################################
my @quantiOrder;# In order to print the table in a good way for R
my %observationInfo; # for (Super/Simple)Ratio only
my (%poolGroup2Cond,%ana2Obs,%cond2Obs,%obs2Cond,%labelModifs,%anaLabelModifs,%xicEngine); #,%obsID2Cond,%ana2Cond
my (%ana2Experiment,%anaConds,%anaState1RatioPos); #<----- NEW EXPERIMENT assignment 31/10/17 -----
my $sthPQN=$dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?"); # peptide quantifs

my @ratioTarget;
my $numConditions=0;
my $poolObservations=0;
#my %experimentCode; # SILAC paired experiments
my @conditionList;
#my %designCheck; # Checks if experimental design matches algo type (label w. new Algo only) TODO: Implement actual test
my $sthALM=$dbh->prepare("SELECT M.ID_MODIFICATION FROM ANALYSIS_MODIFICATION AM,MODIFICATION M WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND AM.ID_ANALYSIS=? AND M.IS_LABEL=1");
if ($call eq 'DESIGN') {
	my $sthObsMod=$dbh->prepare("SELECT M.ID_MODIFICATION FROM OBSERVATION O LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION WHERE O.ID_OBSERVATION=?"); # returns NULL if no modif!!!

	###> Add on 26/10/12 -> To get new STATES info
	foreach my $annoStrg (split (/::/,$quantiAnnotation)) {
		#### STATES ####
		if ($annoStrg =~ /^STATES=(.*)/) { # 'state1;state2;---;stateN',  stateX='nbBioRep,bioRep1.bioRep2.---.bioRepN,#condID', bioRepX='techRep1=techRep2=---=techRepN', techRepX='frac1+frac2+---+fracN'
			##my (%replicAnaStrg,%anaExpCode,%commonRefCondAna); # ,$multiTechRep for (Super/Simple)Ratio only
			##my $lastReplicInRefIdx=-1; # for SuperRatio only: last index of all bio/tech replicates for reference condition in @quantiOrder
#----- NEW EXPERIMENT assignment 31/10/17 ----->
my (@anaOrder,%anaTargetPosCond); # to distinguish different Experiments (labeled ana only)
			foreach my $state (split(/;/,$1)) {
				$poolObservations=1 if $state=~/\+/;
				$numConditions++;
				### Tests fractions if Label-Free ###
				###if ($labeling eq 'FREE' && $poolObservations) {
				###	close DATA;
				###	$sthPQN->finish;
				###	$sthALM->finish;
				###	$sthObsMod->finish;
				###	$dbh->disconnect;
				###	die "ERROR in State #$numConditions: Fractions detected for Label-Free quantification. Correct your design.";
				###}
#----- EXPERIMENT assignment #OLD ----->
				##my $expCode=chr(63+$numConditions); # SuperRatio only 2 -> A, 3 -> B, ... (not used for ref cond=1)
				#$experimentCode{$expCode}=1 if $numConditions >= 2; # skip Ref (1st) Cond
				$state=~s/#//g; # remove all id tags
				my ($numBioRep,$quantiObsIDs,$condID)=split(/,/,$state);
				if ($ratioType eq 'Ratio') {
					push @{$quantifParameters{'R'}{'grp'}},$numBioRep; # Number of replicates of the condition
					push @{$quantifParameters{'R'}{'name.grp'}},"State$numConditions"; # $condID
				}
#elsif ($ratioType eq 'SuperRatio') { # to build primary ratio
	push @conditionList,$condID;
#}
				#elsif ($numConditions==1) {$lastReplicInRefIdx++;} # for SuperRatio only
				my $bioRepCount=0;
				my $allRepCount=0;
				my $multiRep=($quantiObsIDs=~/[.&]/)? 1 : 0; # for Ratio only
				#my %stateAna; # for SuperRatio only
				foreach my $bioReplicate (split(/\./,$quantiObsIDs)) {
					$bioRepCount++;
					my $techRepCount=0;
					foreach my $techReplicate (split(/&/,$bioReplicate)) {
						##$lastReplicInRefIdx++ if $numConditions==1; # for SuperRatio only
						$allRepCount++;
						$techRepCount++;
						###$multiTechRep=1 if $techRepCount>1; # records if exists multiple tech rep
						push @quantiOrder,$techReplicate;
						$obs2Cond{$techReplicate} = $condID;
#print DEBUG "\$techReplicate='$techReplicate' \$condID=$condID \$numConditions=$numConditions\n";
						if ($ratioType eq 'Ratio') {
							print DATA "\tState$numConditions"; # $condID
							print DATA ".rep$allRepCount" if $multiRep;
						}
						#>Required to detect +/- infinite ratios (except %labelModifs)
						$poolGroup2Cond{$techReplicate}=$condID;
						my %replicAna;
						foreach my $fraction (split(/\+/,$techReplicate)) {
							my @obsData=split(/:/,$fraction); # obsID,parQuantifID,anaID,targetPos
							my ($obsID,$parQuantifID,$anaID,$targetPos)=@obsData;
							unless ($xicEngine{$parQuantifID}) {
								$sthPQN->execute($parQuantifID);
								my ($parQuantiAnnot)=$sthPQN->fetchrow_array;
								$xicEngine{$parQuantifID}=($parQuantiAnnot=~/SOFTWARE=MCQ::|EXTRACTION_ALGO=/)? 'MCQ' : 'PD';
							}
							push @{$cond2Obs{$condID}},\@obsData;
							%{$labelModifs{$obsID}}=();
							if ($labeling eq 'FREE') { # LABEL-FREE
								$replicAna{0}=1; # no data linking
								#$stateAna{0}=1; # no data linking
								#$ana2Cond{$anaID}=$condID;
								push @{$ana2Obs{$anaID}},$obsID;
							}
							else {
$anaTargetPosCond{$anaID}{$targetPos}=$condID;
push @anaOrder,$anaID;
$anaConds{$anaID}{$condID}=1; # needed later for infinite ratio check
								$replicAna{$anaID}=1; # for data linking inside same analysis
								#$stateAna{$anaID}=1; # for data linking inside same analysis
								if ($labeling eq 'SILAC') { # SILAC QUANTI
									##if ($ratioType eq 'SuperRatio') {
									##	if ($numConditions==1) {push @{$commonRefCondAna{$techReplicate}},$anaID;} # Common ref condition
									##	else {$anaExpCode{$anaID}=$expCode;} # SuperRatio only
									##}
									unless ($anaLabelModifs{$anaID}) {
										$sthALM->execute($anaID);
										while (my ($modID)=$sthALM->fetchrow_array) {$anaLabelModifs{$anaID}{$modID}=1;}
									}
									$sthObsMod->execute($obsID);
									while (my ($modID)=$sthObsMod->fetchrow_array) {
										$modID=0 unless $modID; # undef for no-label channel
										#$obsID2Cond{$obsID}=$condID;
										#push @{$ana2Obs{$anaID}},$obsID;
										$labelModifs{$obsID}{$modID}=1; # if $modID;
										#$anaLabelModifs{$anaID}{$modID}=1 if $modID;
#print DEBUG "LABEL: $obsID ($targetPos)-> $modID\n";
									}
									#$obsID2Cond{$obsID}=$condID;
									push @{$ana2Obs{$anaID}},$obsID;
								}
								# iTRAQ case is handled when fetching peptide quanti data (<- ??? PP 20/03/15)
								elsif ($labeling=~/ITRAQ|TMT/) {
									my $modID;
									if ($anaLabelModifs{$anaID}) {$modID=(keys %{$anaLabelModifs{$anaID}})[0];}
									else {
										$sthALM->execute($anaID);
										($modID)=$sthALM->fetchrow_array;
						next unless $modID; # !!!TEMP!!! PARAGON quantif with no validated protein (PP 30/03/15)
										$anaLabelModifs{$anaID}{$modID}=1;
									}
									$labelModifs{$obsID}{$modID}=1; # if $modID;
									#$obsID2Cond{$obsID}=$condID;
									push @{$ana2Obs{$anaID}},$obsID;
								}
							}
						} # end of fraction
####						if ($ratioType=~/S\w+Ratio/) {
####							if ($Rdesign eq 'LABELED') { # SPECIAL CASE to handle iTRAQ: bioRep count carried by Exp (techRep?!!!)
####								# Symetrical design between bioRep is MANDATORY: num bioRep StateX = num bioRep StateY!!!!!!
####								#@{$observationInfo{$techReplicate}}=("State$numConditions",'rep1',"techRep$techRepCount",chr(64+$bioRepCount));
####@{$observationInfo{$techReplicate}}=("State$numConditions","rep$bioRepCount","techRep$techRepCount",'A'); #chr(64+$bioRepCount)
####								#$designCheck{chr(64+$bioRepCount)}{$numConditions}=1;
####							}
####							else {
####								@{$observationInfo{$techReplicate}}=("State$numConditions","rep$bioRepCount","techRep$techRepCount");
####								if ($ratioType eq 'SimpleRatio') {
####									push @{$observationInfo{$techReplicate}},'A'; ###'NA';
####									#$designCheck{'NA'}{$numConditions}=1;
####								} # Experiment is irrelevent in table.txt but needed later
####								elsif ($numConditions > 1) { # cond #1 has multiple exeriments
####									push @{$observationInfo{$techReplicate}},$expCode;
####									#$designCheck{$expCode}{$numConditions}=1;
####								}
####							}
####						}
#----- EXPERIMENT assignment #OLD ----->
if ($ratioType=~/S\w+Ratio/) {
	@{$observationInfo{$techReplicate}}=("State$numConditions","rep$bioRepCount","techRep$techRepCount");
	##if ($Rdesign eq 'SUPERRATIO') {
	##	push @{$observationInfo{$techReplicate}},$expCode if $numConditions > 1; # cond #1 has multiple experiments
	##} # LABELFREE & LABELED experiment assignment handled later
	####else { # LABELFREE & LABELED
	####	push @{$observationInfo{$techReplicate}},'A'; # temp for LABELED
	####}
}

					} # end of tech replicate (sample)
				} # end of bio-replicate
				##>Find Experiment letter code for test(s) condition(s)
				##if ($ratioType eq 'SuperRatio') {
				##	my $anaStrg=join('.',sort{$a<=>$b} keys %stateAna); # label-free => '0' (string of all anaID used for state)
				##	if ($numConditions > 1) { # not for common Ref ($numConditions=1)
				##		unless ($experimentCode{$anaStrg}) { # Should be the same between labeled channels to match same experiment
				##			$experimentCode{$anaStrg}=chr(65 + scalar keys %experimentCode); # SuperRatio: 0,1,... -> A,B,...
				##		}
				##		foreach my $bioReplicate (split (/\./,$quantiObsIDs)) { # update %observationInfo
				##			foreach my $techReplicate (split (/=/,$bioReplicate)) {
				##				push @{$observationInfo{$techReplicate}},$experimentCode{$anaStrg}; # index 3
				##				$replicAnaExp{$replicAnaStrg{$techReplicate}}=$experimentCode{$anaStrg};
				##			}
				##		}
				##	}
				##
				##}
			} # end of cond/state
			##>Finding experiment for reference tech replicates
			##if ($ratioType eq 'SuperRatio') {
			##	#my $lastReplicInRefIdx=$quantifParameters{'R'}{'grp'}[0]-1;
			##	foreach my $i (0..$lastReplicInRefIdx) {
			##		my $techReplicate=$quantiOrder[$i]; # update %observationInfo reference replicate only
			##		###push @{$observationInfo{$techReplicate}},$replicAnaExp{$replicAnaStrg{$techReplicate}}; # index 3
			##		my $expCode='NA';
			##		foreach my $anaID (@{$commonRefCondAna{$techReplicate}}) {
			##			if ($anaExpCode{$anaID}) { # 1st match is enough
			##				$expCode=$anaExpCode{$anaID};
			##				last;
			##			}
			##		}
			##		#push @{$observationInfo{$techReplicate}},($expCode,$expCode); # index 3
			##		push @{$observationInfo{$techReplicate}},$expCode; # index 3
			##		#$designCheck{$expCode}{$numConditions}=1;
			##	}
			##}
#----- NEW EXPERIMENT assignment 31/11/17 ----->
# Assumes same Experiment if same set of condIDs associated with same label channels (same targetPos)
# Creates a unique key of tgPos1:condID1[,tgPos2:condID2...] associated with a given Experiment
# 1 Ananlysis belongs to only 1 Experiment
# Should handle:
#    -fractions
#    -(bio/tech)Rep across different Experiments (eg. iTRAQ or TMT)
#    -SuperRatio: state1 found in multiple Experiments
if ($labeling ne 'FREE') {
	my $expCode='A';
	my %condSig2Experiment;
	foreach my $anaID (@anaOrder) { # Exp follows with State pos increase
		next if $ana2Experiment{$anaID};
		my $condSignatureStrg='';
		foreach my $targetPos (sort{$a<=>$b} keys %{$anaTargetPosCond{$anaID}}) {
				#foreach my $condID (sort{$a<=>$b} keys %{$anaTargetPosCond{$anaID}{$targetPos}}) {
				$condSignatureStrg.=',' if $condSignatureStrg;
				$condSignatureStrg.=$targetPos.':'.$anaTargetPosCond{$anaID}{$targetPos}; # targetPos:condID
				#}
		}
		unless ($condSig2Experiment{$condSignatureStrg}) {
			$condSig2Experiment{$condSignatureStrg}=$expCode;
			$expCode++; # 'A'++ -> 'B'!!!
		}
		$ana2Experiment{$anaID}=$condSig2Experiment{$condSignatureStrg};
	}
}
			###if ($ratioType=~/S.+Ratio/ && !$multiTechRep) {
			###	foreach my $techReplicate (keys %observationInfo) {
			###		$observationInfo{$techReplicate}[2]='NA'; # clear techRep info
			###	}
			###}
		} # end of STATES
		#### RATIOS ####
		elsif ($annoStrg =~ /^RATIOS=(.*)/) {
			@ratioTarget=split(/;/,$1);
		}
	}
	$sthObsMod->finish;


}
else { # $call eq 'ANALYSIS' !!!!!!!DEPRECATED!!!!!!!
###	$sthPQN->execute($parentQuantifID);
###	my ($parQuantiAnnot)=$sthPQN->fetchrow_array;
###	$xicEngine{$parentQuantifID}=($parQuantiAnnot=~/SOFTWARE=MCQ::|EXTRACTION_ALGO=/)? 'MCQ' : 'PD';
###	my ($labelStrg,@labelInfo)=split('::',$parQuantiAnnot);
###	if ($labeling eq 'SILAC') { # SILAC QUANTI
###		foreach my $infoStrg (@labelInfo) {
###			next unless $infoStrg=~/^\d+;/; # skip SOFTWARE=... & MCQ params or any future additional tags
###			#my ($chanNum,$chanName,$labelMod,$res,$searchMod)=split(';',$infoStrg); # '1;Light;No label;'   '2;Heavy;Label:13C(6);K'
###			my ($chanNum,$chanName,$labelDataStrg)=split(';',$infoStrg);
###			%{$labelModifs{$chanNum}}=();
###			if (!$labelDataStrg || $labelDataStrg=~/^None#/ || $labelDataStrg=~/No label/ || $labelDataStrg !~/#/) { # assume no-label if missing data
###				$labelModifs{$chanNum}{0}=1;
####print DEBUG "$labelDataStrg -> 0\n";
###			}
###			else {
###				foreach my $label (split(/\@/,$labelDataStrg)) {
###					my ($labelModifAlias,$labelModifName,$modifRes,$searchMod)=split(/#/,$label);
###					my $modID=&promsMod::getModificationIDfromString($dbh,$labelModifName);
####print DEBUG "$labelDataStrg -> $modID\n";
###					$labelModifs{$chanNum}{$modID}=1;
###					$anaLabelModifs{$analysisID}{$modID}=1;
###				}
###			}
###		}
###	}
###	elsif ($labeling=~/ITRAQ|TMT/) {
###		#my ($modID)=$dbh->selectrow_array("SELECT M.ID_MODIFICATION FROM ANALYSIS_MODIFICATION AM,MODIFICATION M WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND AM.ID_ANALYSIS=$analysisID AND M.IS_LABEL=1 LIMIT 0,1");
###		$sthALM->execute($analysisID);
###		my ($modID)=$sthALM->fetchrow_array;
###next unless $modID; # !!!TEMP!!! PARAGON quantif with no validated protein (PP 30/03/15)
###		$anaLabelModifs{$analysisID}{$modID}=1;
###		my $chanNum=0;
###		foreach my $infoStrg (@labelInfo) {
###			next unless $infoStrg=~/^\d+;/; # skip SOFTWARE=... & MCQ params or any future additional tags
###			my ($chanNum,$chanName,$labelMass)=split(';',$infoStrg);
###			$labelModifs{$chanNum}{$modID}=1;
###		}
###	}
###
###	foreach my $annoStrg (split (/::/,$quantiAnnotation)) {
###		if ($annoStrg =~ /^STATES=(.*)/) {
###			foreach my $state (split(/;/,$1)) {
###				$numConditions++;
###				my ($nbRep,$chanNameStrg,$chanNumStrg)=split(/,/,$state); # $quantifParameters{'R'}{'grp/name.grp'} already defined in quantif_info.txt
###				#push @{$quantifParameters{'R'}{'grp'}},$nbRep; # Number of replicates of the condition
###				#push @{$quantifParameters{'R'}{'name.grp'}},$chanNameStrg;
###				my @replicNames=split (/\./,$chanNameStrg);
###				my $repCount=0;
###				foreach my $observation (split(/\./,$chanNumStrg)) { # replicate=channel num
###					$repCount++;
###					push @quantiOrder,$observation;
###					#my $replic=$replicNames[$repCount];
###					#$replic=~s/^\s*-\s*/wo\./; # -XXX -> wo.XXX
###					#$replic=~s/^\s*\+\s*/w\./; # +XXX -> w.XXX
###					#$replic="L$replic" if $replic !~ /^[A-Z]/i; # 114 -> L114
###					#$replic=~s/\s+/_/g;
###					#$replic=&promsMod::shortenName($replic,28);
###					#print DATA "\t$replic";
###					if ($ratioType eq 'Ratio') {
###						print DATA "\tState$numConditions";
###						print DATA ".rep$repCount" if $nbRep > 1;
###					}
###					$obs2Cond{$observation}=$numConditions;
###
###					#>Required to detect +/- infinite ratios (except %labelModifs)
###					$poolGroup2Cond{$observation}=$numConditions;
###					push @{$cond2Obs{$numConditions}},[$observation,$parentQuantifID,$analysisID,$observation]; # to mimick design quantif (replicate=chanel num <=> condID)
###					#if ($labeling eq 'SILAC') {
###						#$obsID2Cond{$observation}=$numConditions;
###						push @{$ana2Obs{$analysisID}},$observation;
###					#}
###
###					if ($ratioType eq 'SimpleRatio') {
###						if ($Rdesign eq 'LABELED') { # SPECIAL CASE to handle iTRAQ: bioRep count carried by Exp (techRep?!!!)
###							# Symetrical design between bioRep is MANDATORY: num bioRep StateX = num bioRep StateY!!!!!!
###							@{$observationInfo{$observation}}=("State$numConditions",'rep1','techRep1',chr(64+$repCount)); # Only 1 (bio/tech)Rep possible
###						}
###						else {
###							@{$observationInfo{$observation}}=("State$numConditions","rep$repCount",'techRep1','A'); # assume bioRep if replicates!!!
###						}
###					}
###				}
###			}
###		}
###		elsif ($annoStrg =~ /^RATIOS=(.*)/) {
###			@ratioTarget=split(/;/,$1);
###		}
###	}
}
$sthPQN->finish;
$sthALM->finish;

#close DEBUG;
#$dbh->disconnect; die "Test is complete!";

###>Handling Modification Quantif phosphoRS threshold & PTM position ambiguity settings
#my $ambiguousModifPos=($quantifiedModifID && $quantifParameters{'DB'}{'AMBIGUOUS_QPTM_POS'})? $quantifParameters{'DB'}{'AMBIGUOUS_QPTM_POS'}[0] : 0;
my ($isPhosphoQuantif,$prsThreshold,$ambiguousModifPos,$quantifResMatchStrg,$qQuantifResMatchStrg)=(0,0,0,'',''); # default ($qQuantifResMatchStrg PP 2017/06/22)
if ($quantifiedModifID) {
	if ($quantifParameters{'DB'}{'PTM_POS'}[0]=~/PRS:(\d+)/) { # Phospho
		$prsThreshold=$1;
		$isPhosphoQuantif=1;
		$ambiguousModifPos=1 if $quantifParameters{'DB'}{'PTM_POS'}[1] eq 'ambiguous';
	}
	elsif ($quantifParameters{'DB'}{'PTM_POS'}[0] eq 'ambiguous') { # other PTMs
		$ambiguousModifPos=1;
	}

	#>Find list of targetted residues: needed even if !$ambiguousModifPos
	my $sthAMR=$dbh->prepare("SELECT SPECIFICITY FROM ANALYSIS_MODIFICATION WHERE ID_MODIFICATION=$quantifiedModifID AND ID_ANALYSIS=? AND MODIF_TYPE='V'");
	my %resList;
	foreach my $anaID (keys %ana2Obs) {
		$sthAMR->execute($anaID);
		my ($specifStrg)=$sthAMR->fetchrow_array;
		foreach my $res (split(',',$specifStrg)) {$resList{$res}=1;}
	}
	$quantifResMatchStrg=join('',sort keys %resList); # eg 'STY'
	$qQuantifResMatchStrg=quotemeta($quantifResMatchStrg); # +* -> \+\* (PP 2017/06/22)
	$sthAMR->finish;
}


#####################################################
####>Protein lists for filtering & normalization<####
#####################################################
my ($protSelectionType,%selectExcludeProteins);
my $sthList=$dbh->prepare("SELECT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=? ORDER BY ID_PROTEIN");

###>Selection/exclusion list
if ($quantifParameters{DB}{PROTEINS}) {
    ($protSelectionType,my $listID)=@{$quantifParameters{DB}{PROTEINS}};
    $listID=~s/#//; # remove id flag
    $sthList->execute($listID);
    while (my ($protID)=$sthList->fetchrow_array) {
        $selectExcludeProteins{$protID}=1;
    }
}

###>Reference proteins for bias correction
my $normProtUsage='';
my (%forNormalization,%notForNormalization);
if ($quantifParameters{DB}{BIAS_CORRECTION}[1]) {
	(my $listID,$normProtUsage)=@{$quantifParameters{DB}{BIAS_CORRECTION}}[1,2];
    $listID=~s/#//; # remove id flag
    $sthList->execute($listID);
	my $refList=($normProtUsage eq 'use')? \%forNormalization : \%notForNormalization;
	while (my ($protID)=$sthList->fetchrow_array) {
		$refList->{$protID}=1;
	}
}

$sthList->finish;


###> Write the table with all results
##################################################################
####>Query to get all the validated proteins and the peptides<####
##################################################################
my ($pepSpecificity,$pepMissedCleavage,$ptmFilter,$pepChargeState,$pepSource)=@{$quantifParameters{'DB'}{'PEPTIDES'}};
#my $pepQuantifMethodCode=($labeling eq 'FREE')? 'XIC' : uc($labeling);
#my $pepQuantifCode=($labeling eq 'FREE')? 'XIC_AREA' : ($labeling eq 'SILAC')? 'ISO_AREA' : 'REP_AREA'; # assumes ITRAQ
#my %pepParamCode2ID;
#my $sthPepQP=$dbh->prepare("SELECT QP.CODE,ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER QP,QUANTIFICATION_METHOD QM WHERE QP.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND QM.CODE='$pepQuantifMethodCode'");
#$sthPepQP->execute;
#while (my ($code,$paramID)=$sthPepQP->fetchrow_array) {
#	$pepParamCode2ID{$code}=$paramID;
#}
#$sthPepQP->finish;
my %peptideAreaParamID;
my $sthPepQP=($labeling eq 'TMT')? $dbh->prepare("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION Q,QUANTIFICATION_PARAMETER QP WHERE Q.ID_QUANTIFICATION_METHOD=QP.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=? AND CODE='REP_INTENSITY'")
								 : $dbh->prepare("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION Q,QUANTIFICATION_PARAMETER QP WHERE Q.ID_QUANTIFICATION_METHOD=QP.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=? AND CODE LIKE '%_AREA'"); # CRITICAL: '%_AREA' !!!!!!!!!!!!!!!!!!!!!!!!!!!
foreach my $pepQuantifID (keys %xicEngine) {
	$sthPepQP->execute($pepQuantifID);
	($peptideAreaParamID{$pepQuantifID})=$sthPepQP->fetchrow_array;
}
$sthPepQP->finish;


##>Normalization with all peptides or only those used for quantification<##
my $normalizeWithAll=0; # 0/1 (0<=>old version) <--- Will become a parameter

my $filterQueryStrg=($normalizeWithAll)? ',PPA.IS_SPECIFIC,AP.PEP_SPECIFICITY,MISS_CUT' : '';

#my $quantiQuery=qq
#|SELECT GROUP_CONCAT(DISTINCT PPA.ID_PROTEIN),P.ID_PEPTIDE,QUANTIF_VALUE,ABS(PEP_BEG),PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),CHARGE,PPA.IS_SPECIFIC,MISS_CUT,DATA,P.SCORE $filterQueryStrg
#	FROM PEPTIDE P
#	LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
#	INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
#	INNER JOIN PEPTIDE_QUANTIFICATION PQ ON P.ID_PEPTIDE=PQ.ID_PEPTIDE
#	INNER JOIN ANALYSIS_PROTEIN AP ON PPA.ID_PROTEIN=AP.ID_PROTEIN
#	WHERE PQ.ID_QUANTIFICATION=? AND P.ID_ANALYSIS=? AND AP.ID_ANALYSIS=? AND ID_QUANTIF_PARAMETER=? AND ABS(AP.VISIBILITY)=2|; # Critical to restrict DISTINCT on ID_PROTEIN for execution time !!! $pepParamCode2ID{$pepQuantifCode}
#unless ($normalizeWithAll) { # filter data
#	$quantiQuery.=($pepSpecificity eq 'unique')? ' AND PPA.IS_SPECIFIC=1' : ($pepSpecificity eq 'unique_shared')? ' AND AP.PEP_SPECIFICITY=100' : ''; # Filter at DB level
#	$quantiQuery.=' AND MISS_CUT=0' unless $pepMissedCleavage;
#}
#$quantiQuery.=($labeling eq 'FREE')? ' AND PQ.TARGET_POS IS NULL' : ' AND PQ.TARGET_POS=?';
#$quantiQuery.=' GROUP BY PPA.ID_PROTEIN,P.ID_PEPTIDE ORDER BY AP.NUM_PEP DESC,ABS(PEP_BEG) ASC';

my $quantiQuery=qq
|SELECT GROUP_CONCAT(DISTINCT PPA.ID_PROTEIN),P.ID_PEPTIDE,ABS(PEP_BEG),PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),CHARGE,PPA.IS_SPECIFIC,MISS_CUT,DATA,P.SCORE $filterQueryStrg
	FROM PEPTIDE P
	LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
	INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
	INNER JOIN ANALYSIS_PROTEIN AP ON PPA.ID_PROTEIN=AP.ID_PROTEIN
	WHERE P.ID_ANALYSIS=AP.ID_ANALYSIS AND AP.ID_ANALYSIS=? AND ABS(AP.VISIBILITY)=2|; # Critical to restrict DISTINCT on ID_PROTEIN for execution time !!! $pepParamCode2ID{$pepQuantifCode}
unless ($normalizeWithAll) { # filter data
	$quantiQuery.=($pepSpecificity eq 'unique')? ' AND PPA.IS_SPECIFIC=1' : ($pepSpecificity eq 'unique_shared')? ' AND AP.PEP_SPECIFICITY=100' : ''; # Filter at DB level
	$quantiQuery.=' AND MISS_CUT=0' unless $pepMissedCleavage;
}
$quantiQuery.=' GROUP BY PPA.ID_PROTEIN,P.ID_PEPTIDE ORDER BY AP.NUM_PEP DESC,ABS(PEP_BEG) ASC';


#print DEBUG ">QUERY---------------\n$quantiQuery\n----------------------\n";
my $sthGetPepData=$dbh->prepare($quantiQuery);

#my ($sthPepRT,$sthRefRT,%refRTcorrection);
#if ($ambiguousModifPos) {
#	#my $targetPosStrg=($labeling eq 'FREE'? 0 : 'TARGET_POS',
#	$sthPepRT=$dbh->prepare("SELECT ID_PEPTIDE,QUANTIF_VALUE FROM PEPTIDE_QUANTIFICATION WHERE ID_QUANTIFICATION=? AND ID_QUANTIF_PARAMETER=$pepParamCode2ID{RT_APEX}");
#	$sthRefRT=$dbh->prepare("SELECT SLOPE,Y_INTERCEPT FROM QUANTIF_REFRT WHERE ID_QUANTIFICATION=? LIMIT 0,1"); # TODO: Add real filtering if more than 1 reference RT
#}

#my $sthPTM=$dbh->prepare("SELECT ID_MODIFICATION FROM PEPTIDE_MODIFICATION WHERE ID_PEPTIDE=?");

$ptmFilter=~s/#//g; # remove id flag
my ($ptmAllowed,@selectedPTMs)=split(/[:,]/,$ptmFilter); # @selectedPTMs defined only if $ptmAllowed >= 2
my %allowedPtmID;
foreach my $modID (@selectedPTMs) {$allowedPtmID{$modID}=1;}
$allowedPtmID{$quantifiedModifID}=1 if $quantifiedModifID; # to allow modif quantification

#my $sthObs=$dbh->prepare("SELECT ID_ANALYSIS,TARGET_POS FROM OBSERVATION WHERE ID_OBSERVATION=?");
#my $sthObsMod=$dbh->prepare("SELECT ID_MODIFICATION FROM OBS_MODIFICATION WHERE ID_OBSERVATION=?");
#my $sthVis=$dbh->prepare("SELECT VISIBILITY FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=? AND ID_PROTEIN=?");
#my %replicatePeptides;

##########################!!!!!!! TEST 18/08/15 !!!!!!!!!!!!!!!!!!!
#my %preInfRatios;
#my $preQuantifID=1532; # MB42
#my $sthPreInf=$dbh->prepare("SELECT ID_PROTEIN FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=? AND ID_QUANTIF_PARAMETER=13 AND (QUANTIF_VALUE >= 1000 OR QUANTIF_VALUE <= 0.001) ORDER BY ID_PROTEIN");
#$sthPreInf->execute($preQuantifID);
#while (my($protID)=$sthPreInf->fetchrow_array) {
#	$preInfRatios{$protID}=1;
#}
#$sthPreInf->finish;
##########################

my (%proteins,%intensity,%proteinAlias,%excludedSeq,%proteinInConds,%peptideBeg,%qSetSequence); # %retentionTime,%bestVisibility,%visInReplicate,%visInAnalysis,%visibility
my (%qSetBestVarMod,%qSetIsBad); # SILAC only
my (%tempIntensity,%dataSource2ObsSet,%qSet2Protein); # for SILAC only
my %notUsed4Quantif;
my %pepID2QuantifSetSILAC; # PP 25/10/17: used to better compute %usedPeptideSets for SILAC when qSets are summed

##----- NEW EXPERIMENT assignment ----->
## RULES: 1 techRep <-> 1 experiment => 1 analysis <-> 1 experiment => 1 pepID <-> 1 experiment
## Should handle:
##    -fractions
##    -bioRep across different Experiments (LABELED with 3+ channels)
##    -SuperRatio: state1 found in multiple Experiments
##my (%analysis2Experiment,%techRep2Experiment); # %peptideID2Experiment,
##my $currentExperiment='A';


####>Retrieving peptides XIC from file (not DB anymore)<####
my (%peptideQuantifs,%pepQuantifValues);
my $sthPepAna=$dbh->prepare("SELECT MIN(ID_PEPTIDE),MAX(ID_PEPTIDE) FROM PEPTIDE WHERE ID_ANALYSIS=?");
foreach my $observationSet (@quantiOrder) {
	foreach my $fraction (split(/\+/,$observationSet)) {
		my ($obsID,$xicQuantif,$anaID,$targetPos)=split(/:/,$fraction);
		$targetPos=0 unless $targetPos; # label-free
		$sthPepAna->execute($anaID);
		my ($minPepID,$maxPepID)=$sthPepAna->fetchrow_array;
		@{$peptideQuantifs{$xicQuantif}{$targetPos}{$anaID}}=($minPepID,$maxPepID);
	}
}
$sthPepAna->finish;
foreach my $xicQuantif (keys %peptideQuantifs) {
	foreach my $targetPos (keys %{$peptideQuantifs{$xicQuantif}}) {
		my $pepQuantifFile=($labeling eq 'FREE')? 'peptide_quantification.txt' : "peptide_quantification_$targetPos.txt";
		open (PEP_QUANTIF,"$promsPath{quantification}/project_$projectID/quanti_$xicQuantif/$pepQuantifFile");
		while(<PEP_QUANTIF>) {
			next if $.==1;
			chomp;
			my ($paramID,$pepID,$quantifValue)=split(/\t/,$_);
			next if $paramID != $peptideAreaParamID{$xicQuantif};
			my $okPepID=0;
			foreach my $anaID (keys %{$peptideQuantifs{$xicQuantif}{$targetPos}}) { # filter for analyses used (multi-ana pep quantifs)
				if ($pepID >= $peptideQuantifs{$xicQuantif}{$targetPos}{$anaID}[0] && $pepID <= $peptideQuantifs{$xicQuantif}{$targetPos}{$anaID}[1]) {
					$pepQuantifValues{$pepID}{$targetPos}=$quantifValue;
					last;
				}
			}
		}
		close PEP_QUANTIF;
	}
}

foreach my $observationSet (@quantiOrder) {
	#my $labelObs=scalar keys %{$labelModifs{$observationSet}}; # observation belongs (or not) to a label channel
	my $anaInObsStrg=''; # list of anaIDs in observationSet
	if ($labeling eq 'FREE' && $ratioType ne 'Ratio') {
		my $condID=$obs2Cond{$observationSet};
		foreach my $refObsData (@{$cond2Obs{$condID}}) {
			my ($obs,$parQuantifID,$anaID,$targetPos)=@{$refObsData};
			$anaInObsStrg.=':' if $anaInObsStrg;
			$anaInObsStrg.=$anaID;
		}
	}
	foreach my $fraction (split(/\+/,$observationSet)) {
		my ($obsID,$xicQuantif,$anaID,$targetPos)=split(/:/,$fraction);
		#if ($labeling eq 'FREE') {
		#}
		#else { # LABELED
#----- NEW EXPERIMENT assignment ----->
#unless ($techRep2Experiment{$observationSet}) {
#	if ($analysis2Experiment{$anaID}) {
#		$techRep2Experiment{$observationSet}=$analysis2Experiment{$anaID};
#	}
#	else {
#		$analysis2Experiment{$anaID}=$techRep2Experiment{$observationSet}=$currentExperiment;
#		$currentExperiment++; # 'A'++ = 'B'!!!
#	}
#}
		#}

##----- NEW EXPERIMENT assignment ----->
#my $expCode='A'; # default
#if ($Rdesign eq 'LABELED' && $numConditions > 2) {
#	if ($analysis2Experiment{$anaID}) {$expCode=$analysis2Experiment{$anaID};}
#	else {
#		$expCode=$analysis2Experiment{$anaID}=$currentExperiment;
#		$currentExperiment++; # A++ = B!!!!
#	}
#}

		my $qSetStrg=($xicEngine{$xicQuantif} eq 'PD')? 'QSET' : "MCQSET_$xicQuantif";
		my %peptidesUsed; # only if non-proteotypic peptides are used => makes sure the same protein is selected (analysis-specific to respect Match group decision!!!)
		my $labelModStrg=join('|',keys %{$anaLabelModifs{$anaID}}) if $labeling eq 'SILAC'; # needed for SILAC
		my %mcqXICs; # SILAC only
		$sthGetPepData->execute($anaID);
		while (my ($protID,$pepID,$pepBeg,$pepSeq,$varModStrg,$charge,$specif,$misscut,$pepData,$score,$isSpecif,$bestSpecif,$missCut)=$sthGetPepData->fetchrow_array) { # executed multiple time for labeled quantif
			next if (!$pepQuantifValues{$pepID} || !$pepQuantifValues{$pepID}{$targetPos} || $pepQuantifValues{$pepID}{$targetPos} <= 0); # no quantif for pepID or Isobaric bug
			my $quantifValue=$pepQuantifValues{$pepID}{$targetPos};
			#next if $quantifValue<=0; # Isobaric bug

			#>Protein filtering
			if ($protSelectionType) {
				if ($protSelectionType eq 'exclude') {next if $selectExcludeProteins{$protID};}
				else {next unless $selectExcludeProteins{$protID};} # restrict
			}
			#>Normalization scope
			$forNormalization{$protID}=1 if ($normProtUsage eq 'not' && !$notForNormalization{$protID});

			$pepData='' unless $pepData;
			$score=0 unless $score;
########!!!!!! TEST 18/08/15
#next if $preInfRatios{$protID};
###################################
#next if $excludedProteins{$protID}; # beta
			if ($quantifiedModifID && (!$varModStrg || $varModStrg !~ /(^|&)$quantifiedModifID:/)) {
				if ($normalizeWithAll) {$notUsed4Quantif{$pepID}=1;}
				else {next;} # skip peptides not containing quantified modif
			}
			if ($pepSpecificity ne 'unique') {  # && $bestSpecif < 100 not proteotypic peptide
				next if ($peptidesUsed{$pepSeq} && $peptidesUsed{$pepSeq} != $protID); # Use peptide ONLY once even if shared by multiple proteins ($pepSpecificity= unique_shared or all) (Required by ML algos)
				$peptidesUsed{$pepSeq}=$protID;
			}
$proteinInConds{$protID}{$poolGroup2Cond{$observationSet}}=1; # protein "exists" in condition (detect before peptide filtering to prevent missleading info!!!)
#$visInReplicate{$protID}{$observationSet}=1; # protein is visible in replicate (VISIBILITY > 0 in query) needed for table.txt only
#$visInAnalysis{$protID}{$anaID}=1;
#$bestVisibility{$protID}=$vis if (!defined($bestVisibility{$protID}) || $bestVisibility{$protID} < $vis); # Keep the Best-Visibility of the protein (only needed for FREE).
			if ($normalizeWithAll) {
				$notUsed4Quantif{$pepID}=1 if (($pepSpecificity eq 'unique' && !$specif) || ($pepSpecificity eq 'unique_shared' && $bestSpecif<100)); # proteo-typic filter
				$notUsed4Quantif{$pepID}=1 if (!$pepMissedCleavage && (!defined $misscut || $misscut > 0)); # misscleavage filter (also skip if no info)
			}

			#if ($ptmAllowed <= 0) {
			#	my $numPTM=0;
			#	$sthPTM->execute($pepID);
			#	while (my ($modID)=$sthPTM->fetchrow_array) {
			#		next if ($allowedPtmID{$modID} || $labelModifs{$obsID}{$modID});
			#		$numPTM++;
			#		last;
			#	}
			#	if ($numPTM) { # modified peptides not allowed
			#		$excludedSeq{$pepSeq}=1 if $ptmAllowed==-1; # unmodified peptide not allowed if exists a modified version
			#		next;
			#	}
			#}
			#my $varMod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$anaID,$pepSeq,$labelModifs{$obsID});

			#>Processing var mods (Checking allowed & removing label if labeled quantif)
			my $varModCode='';
			if ($varModStrg && ($labeling ne 'FREE' || $ptmAllowed != 1)) { # There is filtering on PTM
				my $hasPTM=0;
				foreach my $vMod (split(/&/,$varModStrg)) {
					my ($modID)=($vMod=~/^(\d+)/);
					next if $labelModifs{$obsID}{$modID}; # skip labeling modifs
					if ($ptmAllowed != 1 && !$allowedPtmID{$modID}) {
						$hasPTM=1;
						last;
					}
					$varModCode.='&'.$vMod;
				}
				if ($hasPTM) { # this modification is not allowed
					$excludedSeq{$pepSeq}=1 if $ptmAllowed <= -1; # -1 or -2 unmodified peptide not allowed if exists a modified version
					next;
				}
			}
			elsif ($varModStrg) {$varModCode='&'.$varModStrg;}

			# Non-SILAC phospho quantif
			if ($isPhosphoQuantif && $labeling ne 'SILAC') {
				my ($prsProb)=($pepData=~/PRS=\d;([^;]+);/);
				next if (defined $prsProb && $prsProb < $prsThreshold);
			}

			#$intensity{$pepID}=$quantifValue;
			#unless ($dataSrc) {$dataSrc=($labeling eq 'FREE')? '-' : $anaID;}
			#$dataSrc=($labeling eq 'FREE')? '-' : ($dataSrc)? $dataSrc : '#'.$anaID;
#$dataSrc=$dataSrc || '#'.$anaID; # define dataSrc for label-free also
			my $dataSrc;
# mPP----->
			if ($labeling eq 'FREE' && $ratioType eq 'Ratio') {$dataSrc='-';}
			##if ($labeling eq 'FREE') {$dataSrc='-';} # <---- mPP
# <---- mPP
			else {
				$dataSrc='#'.$anaID; # TODO: if LABEL-FREE: check if still compatible with fraction merge ($anaInObsStrg in Qset) PP 27/11/17
				$dataSrc.='.'.$1 if $pepData=~/SOURCE_RANK=([^#]+)/;
			}
			my $qSet;
			if ($labeling eq 'FREE') { # No peptide pairing
# mPP----->
#$qSet=($ratioType eq 'Ratio')? 'A0_Q0' : "P$pepID"."_Q0"; # <---- mPP
#$qSet=($ratioType eq 'Ratio')? 'A0_Q0' : "A$anaInObsStrg"."_Q0"; # PP 10/08/17 Revert to "A$anaInObsStrg"."_Q0" to allow fractions merge !!!!!!!!!!!!!!! PP 18/04/16: anastrg incompatible with best Qset varMod coherence check
#$qSet=($ratioType eq 'Ratio')? 'A0_Q0' : "A$anaID"."_Q0"; #"A$anaID"."_Q0";!!!!!!!!!!!!!!! PP 10/08/17 WARNING: Should not be compatible with best Qset varMod coherence check
# <---- mPP
				if ($poolObservations) {
					$qSet="A$anaInObsStrg".'_Q0'; # += to cumulate fraction (<=> = if 1 fraction/replicate)
					#$qSet="A$anaID".'_Q0'; += to cumulate same identifications in analysis (<=> = if 1 peptideIon/analysis)
					#$intensity{"$pepSeq:$varModCode:$charge"}{$dataSrc}{$qSet}{$observationSet}[0]+=$quantifValue; #
				}
				else {
					$qSet=($ratioType eq 'Ratio')? 'A0_Q0' : "A$anaID"."_Q0"; # "P$pepID"."_Q0" in algo v2!
					#$intensity{"$pepSeq:$varModCode:$charge"}{$dataSrc}{$qSet}{$observationSet}[0]=$quantifValue;
				}
				$intensity{"$pepSeq:$varModCode:$charge"}{$dataSrc}{$qSet}{$observationSet}[0]=+$quantifValue; # pool just in case multiple occurences
			}
			else { # Labeling

#----- NEW EXPERIMENT assignment ----->
#$peptideID2Experiment{$pepID}=$expCode if $Rdesign eq 'LABELED'; # not for SUPERRATIO

				$qSet="A$anaID".'_Q';
				$qSet.=($pepData=~/$qSetStrg=(\d+)/)? $1 : ($labeling=~/ITRAQ|TMT/)? $pepID : 0;
				if ($labeling eq 'SILAC') {
					if ($varModCode) {
						next if $qSetIsBad{$qSet};
						if ($varModCode=~/&($labelModStrg):/) { # label modif found in wrong label channel => skip qSet
							$qSetIsBad{$qSet}=1;
							next;
						}
					}
					$qSetSequence{$qSet}{$pepSeq}=1; # to check if same pep sequence in qSet
	#push @{$intensity{"$pepSeq:$varModCode:$charge"}{$dataSrc}{$qSet}{$observationSet}},$quantifValue; # multiple QSET for same peptide (old PD SILAC reports before myProMS v3.5)
push @{$tempIntensity{"$pepSeq:#:$varModCode:#:$charge"}{$dataSrc}{$fraction}{$qSet}},$quantifValue; # $intensity computed later
$dataSource2ObsSet{$dataSrc}{$observationSet}=$fraction;
$qSet2Protein{$qSet}=$protID;

$pepID2QuantifSetSILAC{$pepID}=$qSet;
				}
				else {
					push @{$intensity{"$pepSeq:$varModCode:$charge"}{$dataSrc}{$qSet}{$observationSet}},$quantifValue;
				}
			}
			## <--- replaced by above code (PP 13/10/15)
			##if ($poolObservations) {
			##	if ($labeling eq 'FREE') { # All intensities found for a given peptide are summed into a single entry
			##		#$qSet=0;
			##		$qSet='A0_Q0'; # No peptide pairing
			##		$intensity{"$pepSeq:$varModCode:$charge"}{$dataSrc}{$qSet}{$observationSet}[0]+=$quantifValue; # += to cumulate fraction (<=> = if 1 fraction/replicate)
			##	}
			##	else { # SILAC (iTRAQ) intensities are not summed but pooled into the same State.Replicate
			##		$qSet=($pepData=~/$qSetStrg=(\d+)/)? $1 : ($labeling eq 'ITRAQ')? $pepID : 0;
			##		#$qSet.="_A$anaID";
			##		$qSet="A$anaID"."_Q$qSet";
			##		push @{$intensity{"$pepSeq:$varModCode:$charge"}{$dataSrc}{$qSet}{$observationSet}},$quantifValue; # multiple QSET for same peptide (SILAC)
			##	}
			##}
			##else {
			##	if ($labeling eq 'FREE' && $ratioType eq 'Ratio') {
			##		$qSet='A0_Q0'; # No peptide pairing
			##	}
			##	else {
			##		$qSet=($pepData=~/$qSetStrg=(\d+)/)? $1 : ($labeling eq 'ITRAQ')? $pepID : 0;
			##		$qSet="A$anaID"."_Q$qSet";
			##	}
			##	push @{$intensity{"$pepSeq:$varModCode:$charge"}{$dataSrc}{$qSet}{$observationSet}},$quantifValue; # multiple QSET for same peptide (SILAC)
			##}
			##$qSetSequence{$qSet}{$pepSeq}=1 if $labeling eq 'SILAC'; # to check if same pep sequence in qSet

			#my $usedObs=($ratioType eq 'Ratio')? 'ALL_OBS' : $observationSet;
			#push @{$proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet}{$usedObs}},$pepID;
push @{$proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet}{$observationSet}},$pepID;
			$peptideBeg{$protID}{$pepSeq}=$pepBeg; # for modif quantification only

##>Records duplicate XIC scans if MCQ (eg multiple PD XICs reunited in a single one by MCQ => keep only 1 qSet)
if ($xicEngine{$xicQuantif} eq 'MCQ' && $labeling eq 'SILAC') { # $labeling !~ /ITRAQ|TMT/
	$mcqXICs{"$pepSeq:$charge:$quantifValue"}{$qSet}=1; # duplicate if more than 1 qSet with the same pepSeq,charge (ideally mass) have exactly the same XIC value
} # assumes EXACTLY same XIC is computed by MCQ

			##>Record best varMod of qSet in case incoherence (assumes same pepSeq!)
			#if ($labeling ne 'FREE') { #} qSET has no meaning in label-free ***Also for FREE because of possible duplicate scan for same pepID with MassChroQ***
			if ($labeling eq 'SILAC') {
				if ($isPhosphoQuantif) { # PRS filtered earlier for non-SILAC
					my ($prsProb)=($pepData=~/PRS=\d;([^;]+);/);
					$prsProb=0 unless $prsProb;
					@{$qSetBestVarMod{$qSet}}=($prsProb,$varModCode) if (!$qSetBestVarMod{$qSet} || $qSetBestVarMod{$qSet}[0] < $prsProb); # record best phosphoRS score & vmod (both peptides in qSet will be filtered bas on best PRS)
				}
				else {
					@{$qSetBestVarMod{$qSet}}=($score,$varModCode) if (!$qSetBestVarMod{$qSet} || $qSetBestVarMod{$qSet}[0] < $score); # record best
				}
			}
			#if ($ambiguousModifPos) {
			#	$retentionTime{$pepID}=($elutionStrg && $elutionStrg=~/et([\d\.]+)/)? $1 : ($elutionStrg && $elutionStrg=~/^([\d\.]+)/)? $1 : 0;
			#}
		} # end of while

##>Remove duplicate XIC scans if MCQ
if ($labeling eq 'SILAC') {
	foreach my $pepXICkey (keys %mcqXICs) {
		next if scalar keys %{$mcqXICs{$pepXICkey}} == 1; # no duplicate
		my $bestQset;
		foreach my $qSet (sort{$qSetBestVarMod{$b}[0]<=>$qSetBestVarMod{$a}[0]} keys %{$mcqXICs{$pepXICkey}}) {
			if ($bestQset) {$qSetIsBad{$qSet}=1;} # to be excluded
			else {$bestQset=$qSet;} # keep only 1st in list
		}
	}
}

		###>Fetching Apex Retention time
		#if ($ambiguousModifPos) {
		#	$sthRefRT->execute($xicQuantif);
		#	my ($slope,$yIntercept)=$sthRefRT->fetchrow_array;
		#	if ($slope) {$refRTcorrection{'YES'}=1;}
		#	else {
		#		$slope=1;
		#		$yIntercept=0;
		#		$refRTcorrection{'NO'}=1;
		#	}
		#	$sthPepRT->execute($xicQuantif);
		#	while (my ($pepID,$retTime)=$sthPepRT->fetchrow_array) {
		#		$retentionTime{$pepID}=$slope*$retTime + $yIntercept if $peptidesUsed{$pepID}; # COMMENT: Some targetPos might not be used for Prot quantif...
		#	}
		#}
	}
}
$sthGetPepData->finish;
#$sthPTM->finish;
%pepQuantifValues=(); # hoping to free memory

if ($labeling eq 'SILAC') {

	##>Checking qSet sequence
	foreach my $qSet (keys %qSetSequence) { # It happens that identifications are different between Light and (Medium/)Heavy spectra!!!!
		$qSetIsBad{$qSet}=1 if scalar keys %{$qSetSequence{$qSet}} > 1;
	}

	##>Cumulate all qSets matching same peptide ion for a given source (& fraction) & update %proteins
	foreach my $pepIon (keys %tempIntensity) {
		my ($pepSeq,$varModCode,$charge)=split(/:#:/,$pepIon); $varModCode='' unless $varModCode;
		foreach my $dataSrc (keys %{$tempIntensity{$pepIon}}) {
			my (%qSetList,%sumQuantifValues,%qSetObsList);
			foreach my $observationSet (keys %{$dataSource2ObsSet{$dataSrc}}) {
				my $fraction=$dataSource2ObsSet{$dataSrc}{$observationSet};
				next unless $tempIntensity{$pepIon}{$dataSrc}{$fraction};
				foreach my $qSet (keys %{$tempIntensity{$pepIon}{$dataSrc}{$fraction}}) {
					next if $qSetIsBad{$qSet}; # bad qset (matches more than 1 sequence or labeling issues)
					next if ($isPhosphoQuantif && $qSetBestVarMod{$qSet}[0] < $prsThreshold && !$ambiguousModifPos);
#next if $qSetBestVarMod{$qSet}[1] ne $varModCode; # qSet should be associated with another PTM isoform
					my $usedVmod=$qSetBestVarMod{$qSet}[1]; # can be different from $varModCode (eg. incoherent Phospho pos in same qSet)
					$qSetList{$usedVmod}{$qSet}=1;
					push @{$qSetObsList{$observationSet}{$usedVmod}},$qSet;
					$sumQuantifValues{$observationSet}{$usedVmod}+=$tempIntensity{$pepIon}{$dataSrc}{$fraction}{$qSet}[0];
				}
			}
			foreach my $usedVmod (keys %qSetList) {
				my $qSetMerged=join('+',sort keys %{$qSetList{$usedVmod}}); # same qSet string for all observationSets used by dataSrc to make sure qSet pairing is maintained!!! (in case of singleton qSets)

				##>Merged qSet
				if ($qSetMerged=~/\+/) {
					$qSetMerged=~s/\+A\d+_/\+/g; # remove repeated A<anaID>_ except 1srt one: Axxx_Qyyy+Axxx_Qzzz -> Axxx_Qyyy+Qzzz
					foreach my $observationSet (keys %{$dataSource2ObsSet{$dataSrc}}) {
						next if (!$sumQuantifValues{$observationSet} || !$sumQuantifValues{$observationSet}{$usedVmod});

						$intensity{"$pepSeq:$usedVmod:$charge"}{$dataSrc}{$qSetMerged}{$observationSet}[0]=$sumQuantifValues{$observationSet}{$usedVmod};

						##>Correcting %proteins with merged qSets
						my $protID=$qSet2Protein{$qSetObsList{$observationSet}{$usedVmod}[0]};
						my @pepIdList;
						foreach my $qSet (@{$qSetObsList{$observationSet}{$usedVmod}}) {
							push @pepIdList,$proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet}{$observationSet}[0];
#print DEBUG "+PROT=$protID Q=$qSet PEP=$pepSeq VM=$varModCode Z=$charge S=$dataSrc OBS=$observationSet\n" unless $proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet}{$observationSet}[0];
						}
						$proteins{$protID}{$pepSeq}{$usedVmod}{$charge}{$dataSrc}{$qSetMerged}{$observationSet}[0]=join('+',@pepIdList);
					}

					foreach my $qSet (keys %{$qSetList{$usedVmod}}) {
						##>Deleting obsolete qSets
						delete $proteins{$qSet2Protein{$qSet}}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet};

						##>Updating %qSetBestVarMod
#print DEBUG "P=$pepIon S=$dataSrc Q=$qSet\n" unless $qSetBestVarMod{$qSet};
	#					next unless $qSetBestVarMod{$qSet};
						@{$qSetBestVarMod{$qSetMerged}}=@{$qSetBestVarMod{$qSet}} if (!$qSetBestVarMod{$qSetMerged} || $qSetBestVarMod{$qSetMerged}[0] < $qSetBestVarMod{$qSet}[0]); # record best
						# DO NOT DELETE $qSetBestVarMod{$qSet} -> Can be reused in another pepIon context!
					}
				}
				##>No qSet merge (qSetMerged = qSet) => simple update of %intensity
				else {
					my $qSet=$qSetMerged;
					my $protID=$qSet2Protein{$qSet};
					foreach my $observationSet (keys %{$dataSource2ObsSet{$dataSrc}}) {
						next if (!$sumQuantifValues{$observationSet} || !$sumQuantifValues{$observationSet}{$usedVmod}); # no data for this channel
						$intensity{"$pepSeq:$usedVmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0]=$sumQuantifValues{$observationSet}{$usedVmod}; # if ($sumQuantifValues{$observationSet} && $sumQuantifValues{$observationSet}{$usedVmod});
						if ($usedVmod ne $varModCode) {
#my $okString='';
#if ($proteins{$protID}) {
#	$okString.='protID_OK';
#	if ($proteins{$protID}{$pepSeq}) {
#		$okString.=' pepSeq_OK';
#		if ($proteins{$protID}{$pepSeq}{$varModCode}) {
#			$okString.=' varModCode_OK';
#			if ($proteins{$protID}{$pepSeq}{$varModCode}{$charge}) {
#				$okString.=' charge_OK';
#				if ($proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}) {
#					$okString.=' dataSrc_OK';
#					if ($proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet}) {
#						$okString.=' qSet_OK';
#						if ($proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet}{$observationSet}) {
#							$okString.=' observationSet_OK';
#						}
#					}
#				}
#			}
#		}
#	}
#}
#print DEBUG ">$protID: '$usedVmod' <=> '$varModCode' $pepSeq $charge $dataSrc $qSet ($observationSet) => $okString\n";
							@{$proteins{$protID}{$pepSeq}{$usedVmod}{$charge}{$dataSrc}{$qSet}{$observationSet}}=@{$proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet}{$observationSet}};
							#delete $proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet};
							delete $proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet}{$observationSet};
							delete $proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet} unless scalar keys %{$proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet}};
						}
					}
				}
			}
		} # end of dataSource
	}
} # end of SILAC

#close DEBUG;
#die "End of debug!";

#if ($ambiguousModifPos) {
#	$sthPepRT->finish;
#	$sthRefRT->finish;
#	if (scalar keys %refRTcorrection==2) {
#		$dbh->disconnect;
#		die "Incompatible RT scales used by peptide quantifications. Make sure all quantifications have an RT reference.";
#	}
#}

my $sthgetProtInfo=$dbh->prepare("SELECT ALIAS FROM PROTEIN WHERE ID_PROTEIN=?");
foreach my $protID (keys %proteins) {
	$sthgetProtInfo->execute($protID);
	($proteinAlias{$protID})=$sthgetProtInfo->fetchrow_array;
}
$sthgetProtInfo->finish;


######################################
####> Printing normalization file<####
######################################
if ($normProtUsage) {
	unless (scalar keys %forNormalization) {
		$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID"); # Failed
		$dbh->commit;
		$dbh->disconnect;
		die "ERROR: No proteins found for normalization!\n";
	}
	if ($ratioType eq 'Ratio') {
		foreach my $protID (keys %forNormalization) {
			push @{$quantifParameters{'R'}{'prot.ref'}},$protID;
		}
	}
	else {
		open(NORM,">$dataDir/normProtein.txt");
		print NORM "proteinId\n"; # header
		foreach my $protID (sort{$a<=>$b} keys %forNormalization) {
			print NORM "$protID\n";
		}
		close NORM;
	}
}

##############################
####> Printing data table<####
##############################
my (%infiniteRatioProteins,%noSuperRatioProteins); # populated by &checkInfiniteRatios as globals (infinite ratio proteins are allowed in data table but stats are ignored)
my %qSetUsage=('USED'=>{},'EXCLUDED'=>{}); # in case  qSet filtering (eg. PhosphoRS)
my $protNumber=0;
#my %prot2quantify;
my $countCompletePeptides = 0;
my $numPepValues=0;
PROT:foreach my $protID (sort{$a<=>$b} keys %proteins) {
#last PROT if $protNumber >=200; #!!!!!!!!!!!!!!!!
#print DEBUG ">$protID:\n";
	my (%usableData,%numPepInCond,%numPepInCond4ModProt);
	foreach my $pepSeq (keys %{$proteins{$protID}}) {
		next if $excludedSeq{$pepSeq}; # excluded by ptmAllowed <= -1 or shared by multiple top proteins
		#$prot2quantify{$protID}=1;
		foreach my $vmod (keys %{$proteins{$protID}{$pepSeq}}) {
			#>Finding usable charge state(s) & data sources
			#my %usableData;
			my ($bestCharge,$bestSource);
			my $bestValue=0;
			foreach my $charge (keys %{$proteins{$protID}{$pepSeq}{$vmod}}) {
				$bestValue=0 if $pepChargeState eq 'all'; # reset bestValue between charges
				foreach my $dataSrc (keys %{$proteins{$protID}{$pepSeq}{$vmod}{$charge}}) {
					if ($pepChargeState eq 'best' || $pepSource eq 'best') { # some data filtering applies
						my $setValue=0;
						foreach my $qSet (keys %{$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}}){
							next if $qSetIsBad{$qSet}; # only for SILAC
							foreach my $observationSet (keys %{$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}}){
								$setValue+=$intensity{"$pepSeq:$vmod:$charge"}{$observationSet}{$dataSrc}{$qSet}[0]; # take first! => TODO: find a way to deal with multiple values here!!!
							}
						}
						if ($setValue > $bestValue) {
							$bestValue=$setValue;
							$bestCharge=$charge;
							$bestSource=$dataSrc;
						}
					}
					else {$usableData{$pepSeq}{$vmod}{$charge}{$dataSrc}=1;} # all data
				}
				if ($pepChargeState eq 'all' && $pepSource eq 'best' && $bestSource) {$usableData{$pepSeq}{$vmod}{$charge}{$bestSource}=1;} # all charges but 1 source/charge
			}
			#if ($pepChargeState eq 'best' && $bestSource) {$usableData{$pepSeq}{$vmod}{$bestCharge}{$bestSource}=1;} # best charge => best source : 1 peptide set
			if ($pepChargeState eq 'best') { # PP 06/05/16 best charge used across all sources
				foreach my $dataSrc (keys %{$proteins{$protID}{$pepSeq}{$vmod}{$bestCharge}}) {
					$usableData{$pepSeq}{$vmod}{$bestCharge}{$dataSrc}=1;
				}
			}
		}
	}

	if ($labeling eq 'FREE') { # potential data filtering --->
		my %excludedData;

		##>Peptides must match across conditions (Never for PTM quantif)
		if ($matchingPep) { # && $ratioType eq 'SimpleRatio'
			my (%pepMeanValues,%checkExclusion);
			foreach my $pepSeq (keys %usableData) {
				foreach my $vmod (keys %{$usableData{$pepSeq}}) {
					foreach my $charge (sort{$a<=>$b} keys %{$usableData{$pepSeq}{$vmod}}) {
						next unless $intensity{"$pepSeq:$vmod:$charge"};
						my (%pepIonValues,%pepIonInCond);
						foreach my $dataSrc (keys %{$usableData{$pepSeq}{$vmod}{$charge}}) {
							next unless $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc};
							foreach my $qSet (keys %{$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}}){
								foreach my $observationSet (keys %{$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}}) {
									$pepIonValues{"$qSet:#:$observationSet:#:$pepSeq:#:$vmod:#:$charge:#:$dataSrc"}=$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0];
									$pepIonInCond{$obs2Cond{$observationSet}}++; # num pep in condition
								}
							}
						}
						my $pepIonKey="$pepSeq:$vmod:$charge";
						if (scalar keys %pepIonInCond < $numConditions) { # (< 2 keep any pairs? for 3+ conditions) peptide is missing in 1+ cond => flag to be deleted
							foreach my $dataKey (keys %pepIonValues) {
								#$excludedData{$dataKey}=1;
								$checkExclusion{$dataKey}=1; # do not exclude yet. could become +/-INF protein
							}
							next; # charge
						}
						foreach my $dataKey (keys %pepIonValues) {
							$pepMeanValues{$pepIonKey}{'MEAN'}+=$pepIonValues{$dataKey};
							push @{$pepMeanValues{$pepIonKey}{'DATA_KEY'}},$dataKey;
						}
						$pepMeanValues{$pepIonKey}{'NUM_VALUES'} = scalar (keys %pepIonValues);
						$pepMeanValues{$pepIonKey}{'MEAN'} /= $pepMeanValues{$pepIonKey}{'NUM_VALUES'}; # mean of peptide ions across all replicates of all conditions
$pepMeanValues{$pepIonKey}{'SEQ'}=$pepSeq;
					}
				}
			}

			##>TopN quantif
			if ($topN && !$quantifiedModifID) { # keep top<N>
my (%usedPepSeq,@backUpPepIons);
				my $numPep=0;
				foreach my $pepIonKey (sort{$pepMeanValues{$b}{'NUM_VALUES'}<=>$pepMeanValues{$a}{'NUM_VALUES'} || $pepMeanValues{$b}{'MEAN'}<=>$pepMeanValues{$a}{'MEAN'}} keys %pepMeanValues) { # most frequent then most intense
my $pepSeq=$pepMeanValues{$pepIonKey}{'SEQ'};
if ($usedPepSeq{$pepSeq}) { # Do not allow same pepSeq to be used twice even if charge & PTM are different (PP 14/06/18)
	push @backUpPepIons,$pepIonKey; # keep in case numPep stays < topN
	next;
}
$usedPepSeq{$pepSeq}=1;
					$numPep++;
					if ($numPep > $topN) { # flag to be deleted
						foreach my $dataKey (@{$pepMeanValues{$pepIonKey}{'DATA_KEY'}}) {$excludedData{$dataKey}=1;}
					}
				}
foreach my $pepIonKey (@backUpPepIons) { # Keep only backed up already-used pep seq for numPep <= topN
	$numPep++;
	if ($numPep > $topN) { # flag to be deleted
		foreach my $dataKey (@{$pepMeanValues{$pepIonKey}{'DATA_KEY'}}) {$excludedData{$dataKey}=1;}
	}
}
			}
			if (scalar keys %pepMeanValues) { # matching peptide(s) found => delete all non-matching peptides
				foreach my $dataKey (keys %checkExclusion) {$excludedData{$dataKey}=1;}
			}
#else { # "numPep=0" > infinite ratios
#print DEBUG "$protID\t",scalar keys %checkExclusion,"\n";
#	#foreach my $dataKey (keys %checkExclusion) {
#	#	$excludedData{$dataKey}=1;
#	#}
#}
		}

		##>TopN quantif (no peptide matching)
		elsif ($topN) { # && $ratioType eq 'SimpleRatio'	WARNING: +/-inf prot are also restricted to topN peptides
			my %topNvalues;
			foreach my $pepSeq (keys %usableData) {
				foreach my $vmod (keys %{$usableData{$pepSeq}}) {
					foreach my $charge (sort{$a<=>$b} keys %{$usableData{$pepSeq}{$vmod}}) {
						next unless $intensity{"$pepSeq:$vmod:$charge"};
						foreach my $dataSrc (keys %{$usableData{$pepSeq}{$vmod}{$charge}}) {
							next unless $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc};
							#foreach my $qSet (keys %{$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}}){ #}
							#	next unless $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet};
							foreach my $qSet (keys %{$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}}){
								foreach my $observationSet (keys %{$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}}) {
									my $newValue=$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0];
									my $qSetKey="$qSet:#:$observationSet";
									if ($topNvalues{$qSetKey} && scalar (keys %{$topNvalues{$qSetKey}}) == $topN) { # $qset = P<pepID>_Q0 for label free SimpleRatio (no fraction expected)
										foreach my $dataKey (sort{$topNvalues{$qSetKey}{$a} <=> $topNvalues{$qSetKey}{$b}} keys %{$topNvalues{$qSetKey}}) {
											if ($newValue > $topNvalues{$qSetKey}{$dataKey}) {
												$excludedData{"$qSetKey:#:$dataKey"}=$topNvalues{$qSetKey}{$dataKey}; # exclude previous value
												delete $topNvalues{$qSetKey}{$dataKey};
												$topNvalues{$qSetKey}{"$pepSeq:#:$vmod:#:$charge:#:$dataSrc"}=$newValue;
											}
											else {
												$excludedData{"$qSetKey:#:$pepSeq:#:$vmod:#:$charge:#:$dataSrc"}=$newValue; # exclude new value
											}
											last; # Only 1 test is necessary because of sort{$a<=>$b}
										}
									}
									else {
										$topNvalues{$qSetKey}{"$pepSeq:#:$vmod:#:$charge:#:$dataSrc"}=$newValue;
									}
								}
							}
						}
					}
				}
			}
		}

		##>Update usableData
		foreach my $dataKey (keys %excludedData) {
			my ($qSet,$observationSet,$pepSeq,$vmod,$charge,$dataSrc)=split(':#:',$dataKey);
			delete $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet};
			unless (scalar keys %{$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}}) {
				delete $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet};
				unless (scalar keys %{$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}}) {
					delete $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc};
					unless (scalar keys %{$intensity{"$pepSeq:$vmod:$charge"}}) {
						delete $intensity{"$pepSeq:$vmod:$charge"};
					}
delete $usableData{$pepSeq}{$vmod}{$charge}{$dataSrc};
unless (scalar keys %{$usableData{$pepSeq}{$vmod}{$charge}}) {
	delete $usableData{$pepSeq}{$vmod}{$charge};
	unless (scalar keys %{$usableData{$pepSeq}{$vmod}}) {
		delete $usableData{$pepSeq}{$vmod};
		unless (scalar keys %{$usableData{$pepSeq}}) {
			delete $usableData{$pepSeq};
		}
	}
}
				}
			}
		}
	}

	##>Printing data to file (table.txt)
	$protNumber++ if scalar keys %usableData;
	foreach my $pepSeq (sort{$peptideBeg{$protID}{$a}<=>$peptideBeg{$protID}{$b} || length($a)<=>length($b)} keys %usableData) {
		my @seq=split(//,$pepSeq);
		foreach my $vmod (keys %{$usableData{$pepSeq}}) {
			my $nonAmbigModProtID=$protID; # default
			my $ambigModProtID=$nonAmbigModProtID; # default: no ambiguity
			my $seqIsAmbig=0;
# (PP 2017/06/22) --->
			if ($quantifiedModifID) { # modif quantification
				my ($sep,$posStrg)=($vmod=~/(^|&)$quantifiedModifID:([\d+\.\-\=\+\*]+)/); # added N/C-term
				#>Non-ambiguous pos
				my @pos;
				foreach my $pepPos (split(/\./,$posStrg)) {
					if ($pepPos=~/\d+/) {push @pos,$seq[$pepPos-1].($peptideBeg{$protID}{$pepSeq}+$pepPos-1);}
					#elsif ($pepPos=~/[=\*]/) {next;} # skip peptide N/C-term: cannot be in protein TODO: add to table.txt for data normalisation but do not quantify
					elsif ($pepPos eq '=') {push @pos,'n'.$peptideBeg{$protID}{$pepSeq};} # peptide N-term: kept for normalization purpose
					elsif ($pepPos eq '-') {push @pos,'n0';} # protein N-term: kept for normalization purpose
					elsif ($pepPos eq '*') {push @pos,'c'.($peptideBeg{$protID}{$pepSeq}+$#seq);} # peptide C-term: kept for normalization purpose
					elsif ($pepPos eq '+') {push @pos,'c0';} # protein C-term: kept for normalization purpose
				}
				$nonAmbigModProtID=$protID.'-'.join('.',@pos);
				$ambigModProtID=$nonAmbigModProtID; # default: no ambiguity
				#if ($ambiguousModifPos) { # true for Phospho quantif only if <PRS allowed
					#my $numModifRes = () = $posStrg =~ /(\d+)/g;
					my $numModifRes=scalar @pos;
					my @targetResIdx;
					while ($pepSeq =~ /[$qQuantifResMatchStrg]/g) {push @targetResIdx,$-[0];} # position of match (1st=0)
					my $numTargetRes=scalar @targetResIdx;
					my ($NtermMatch,$CtermMatch)=('','');
					if ($quantifResMatchStrg=~/=/) {$NtermMatch='n'.$peptideBeg{$protID}{$pepSeq};} # peptide N-term
					elsif ($quantifResMatchStrg=~/\*/) {$CtermMatch='c'.($peptideBeg{$protID}{$pepSeq}+$#seq);} # peptide C-term
					elsif ($quantifResMatchStrg=~/-/) {$NtermMatch='n0' if $peptideBeg{$protID}{$pepSeq}==1;} # protein N-term
					#elsif ($quantifResMatchStrg=~/\+/) {$CtermMatch='c0' if $peptideBeg{$protID}{$pepSeq}+length($pepSeq)-1==<protein length>;} # TODO: fetch prot length
					my $numAlltargets=$numTargetRes;
					$numAlltargets++ if $NtermMatch;
					$numAlltargets++ if $CtermMatch;
					if ($numModifRes < $numAlltargets) { # Ambiguity! protID-<1stPos>~<lastPos>:<numModif>/<numSites>
						$seqIsAmbig=1;
						if ($ambiguousModifPos) {
							#$ambigModProtID=$protID.'-'.($peptideBeg{$protID}{$pepSeq}+$targetResIdx[0]).'~'.($peptideBeg{$protID}{$pepSeq}+$targetResIdx[-1]).':'.$numModifRes.'/'.$numTargetRes;
							$ambigModProtID=$protID.'-';
							if ($NtermMatch) {$ambigModProtID.=$NtermMatch;} # N-term
							else {$ambigModProtID.=($peptideBeg{$protID}{$pepSeq}+$targetResIdx[0]);} # internal pos
							$ambigModProtID.='~';
							if ($CtermMatch) {$ambigModProtID.=$CtermMatch;} # C-term
							else {$ambigModProtID.=($peptideBeg{$protID}{$pepSeq}+$targetResIdx[-1]);} # internal pos
							$ambigModProtID.=':'.$numModifRes.'/'.$numAlltargets;
						}
					}
				#}
# <--- (PP 2017/06/22)
			}
			#my %numPepInCond4ModProt;
			foreach my $charge (sort{$a<=>$b} keys %{$usableData{$pepSeq}{$vmod}}) {
				foreach my $dataSrc (sort keys %{$usableData{$pepSeq}{$vmod}{$charge}}) {
					my $infRatioDataSrc=($labeling eq 'FREE')? '-' : $dataSrc; # ignore dataSrc if label-free

					foreach my $qSet (sort{$a cmp $b} keys %{$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}}) { # sort cmp because Qset can be string
						next if $qSetIsBad{$qSet}; # only for SILAC or duplicate MCQ
						my $modProtID=$nonAmbigModProtID; # default
						my $usedVmod=($labeling eq 'SILAC')? $qSetBestVarMod{$qSet}[1] : $vmod; # can be different from $vmod in case qSet incoherence (eg. phospho pos after phosphoRS update on ambiguous pepSeq)
						if ($seqIsAmbig && $labeling eq 'SILAC') {
							if ($isPhosphoQuantif) {
								if ($qSetBestVarMod{$qSet}[0] >= $prsThreshold) { # check Phospho pos coherence inside qSet (=0 if no PRS data)
									if ($qSetBestVarMod{$qSet}[2]) { # best nonAmbigModProtID has already been generated
										###if ($qSetBestVarMod{$qSet}[1] ne $vmod) { # *** Position incoherence detected!!! ***
										###	$usedVmod=$qSetBestVarMod{$qSet}[1];
										###}
									}
									else { # best nonAmbigModProtID has not been generated yet => process now!
										if ($usedVmod eq $vmod) {
											$qSetBestVarMod{$qSet}[2]=$nonAmbigModProtID;
										}
										else { # *** Position incoherence detected!!! *** detected only if best vmod is 2nd to be scanned
											my ($sep,$posStrg)=($usedVmod=~/(^|&)$quantifiedModifID:([\d+\.\-\=\+\*]+)/); # (PP 2017/06/22)
											my @pos;
											#foreach my $pos (split(/\./,$posStrg)) {push @pos,$seq[$pos-1].($peptideBeg{$protID}{$pepSeq}+$pos-1);}
											foreach my $pos (split(/\./,$posStrg)) {
if (!$pos || ($pos=~/\d+/ && (!$seq[$pos-1])) || !$peptideBeg{$protID}{$pepSeq}) { # (PP 2017/06/22)
	die "Pos Error: PROT=$protID, PEP=$pepSeq, VMOD=$vmod, UVMOD=$usedVmod, QSET=$qSet, SEP=$sep, STRG=$posStrg, BEG=$peptideBeg{$protID}{$pepSeq}\n";
}

												push @pos,$seq[$pos-1].($peptideBeg{$protID}{$pepSeq}+$pos-1);
											}
											$qSetBestVarMod{$qSet}[2]=$protID.'-'.join('.',@pos);
											###$usedVmod=$qSetBestVarMod{$qSet}[1];
										}
									}
									$modProtID=$qSetBestVarMod{$qSet}[2];
								}
								else { # must be an ambiguous peptide (PRS < 100%)
									if ($ambiguousModifPos) {
										$modProtID=$ambigModProtID;
									} # allow quantif below PRS as ambigous
									else {
										$qSetUsage{'EXCLUDED'}{$qSet}=1;
										next;
									} # exclude pep below PRS
								}
							}
							else {$modProtID=$ambigModProtID;} # non-Phospho PTM quantified as ambiguous
						}

####---------------------------->
###$proteoForm{$modProtID}{$pepSeq}{$vmod}{$charge}{$dataSrc}{$qSet}=1 if (!$qSetUsage{'EXCLUDED'} || !$qSetUsage{'EXCLUDED'}{$qSet});
###
###					} # $qSet
###				} # $dataSrc
###			} # $vmod
###		} # $pepSeq
###
###foreach my $modProtID (sort keys %proteForm) {
###
###	if ($labeling eq 'FREE') { # potential data filtering --->
###		my %excludedData;
###
###		##>Peptides must match across conditions (Never for PTM quantif)
###		if ($matchingPep) { # && $ratioType eq 'SimpleRatio'
###			my (%pepMeanValues,%checkExclusion);
###			foreach my $pepSeq (keys %{$proteForm{$modProtID}}) {
###				foreach my $vmod (keys %{$proteForm{$modProtID}{$pepSeq}}) {
###					foreach my $charge (sort{$a<=>$b} keys %{$proteForm{$modProtID}{$pepSeq}{$vmod}}) {
###next unless $intensity{"$pepSeq:$vmod:$charge"};
###						my (%pepIonValues,%pepIonInCond);
###						foreach my $dataSrc (keys %{$proteForm{$modProtID}{$pepSeq}{$vmod}{$charge}}) {
###next unless $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc};
###							foreach my $qSet (keys %{$proteForm{$modProtID}{$pepSeq}{$vmod}{$charge}{$dataSrc}}){
###								foreach my $observationSet (keys %{$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}}) {
###									$pepIonValues{"$qSet:#:$observationSet:#:$pepSeq:#:$vmod:#:$charge:#:$dataSrc"}=$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0];
###									$pepIonInCond{$obs2Cond{$observationSet}}++; # num pep in condition
###								}
###							}
###						}
###						my $pepIonKey="$pepSeq:$vmod:$charge";
###						if (scalar keys %pepIonInCond < $numConditions) { # (< 2 keep any pairs? for 3+ conditions) peptide is missing in 1+ cond => flag to be deleted
###							foreach my $dataKey (keys %pepIonValues) {
###								#$excludedData{$dataKey}=1;
###								$checkExclusion{$dataKey}=1; # do not exclude yet. could become +/-INF protein
###							}
###							next; # charge
###						}
###						foreach my $dataKey (keys %pepIonValues) {
###							$pepMeanValues{$pepIonKey}{'MEAN'}+=$pepIonValues{$dataKey};
###							push @{$pepMeanValues{$pepIonKey}{'DATA_KEY'}},$dataKey;
###						}
###						$pepMeanValues{$pepIonKey}{'NUM_VALUES'} = scalar (keys %pepIonValues);
###						$pepMeanValues{$pepIonKey}{'MEAN'} /= $pepMeanValues{$pepIonKey}{'NUM_VALUES'}; # mean of peptide ions across all replicates of all conditions
###					}
###				}
###			}
###
###			##>TopN quantif
###			if ($topN && !$quantifiedModifID) { # keep top<N>
###				my $numPep=0;
###				foreach my $pepIonKey (sort{$pepMeanValues{$b}{'NUM_VALUES'}<=>$pepMeanValues{$a}{'NUM_VALUES'} || $pepMeanValues{$b}{'MEAN'}<=>$pepMeanValues{$a}{'MEAN'}} keys %pepMeanValues) { # most frequent then most intense
###					$numPep++;
###					if ($numPep > $topN) { # flag to be deleted
###						foreach my $dataKey (@{$pepMeanValues{$pepIonKey}{'DATA_KEY'}}) {$excludedData{$dataKey}=1;}
###					}
###				}
###			}
###			if (scalar keys %pepMeanValues) { # matching peptide(s) found => delete all non-matching peptides
###				foreach my $dataKey (keys %checkExclusion) {$excludedData{$dataKey}=1;}
###			}
####else { # "numPep=0" > infinite ratios
####print DEBUG "$protID\t",scalar keys %checkExclusion,"\n";
####	#foreach my $dataKey (keys %checkExclusion) {
####	#	$excludedData{$dataKey}=1;
####	#}
####}
###		}
###
###		##>TopN quantif (no peptide matching)
###		elsif ($topN) { # && $ratioType eq 'SimpleRatio'	WARNING: +/-inf prot are also restricted to topN peptides
###			my %topNvalues;
###			foreach my $pepSeq (keys %usableData) {
###				foreach my $vmod (keys %{$usableData{$pepSeq}}) {
###					foreach my $charge (sort{$a<=>$b} keys %{$usableData{$pepSeq}{$vmod}}) {
###						next unless $intensity{"$pepSeq:$vmod:$charge"};
###						foreach my $dataSrc (keys %{$usableData{$pepSeq}{$vmod}{$charge}}) {
###							next unless $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc};
###							#foreach my $qSet (keys %{$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}}){ #}
###							#	next unless $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet};
###							foreach my $qSet (keys %{$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}}){
###								foreach my $observationSet (keys %{$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}}) {
###									my $newValue=$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0];
###									my $qSetKey="$qSet:#:$observationSet";
###									if ($topNvalues{$qSetKey} && scalar (keys %{$topNvalues{$qSetKey}}) == $topN) { # $qset = P<pepID>_Q0 for label free SimpleRatio (no fraction expected)
###										foreach my $dataKey (sort{$topNvalues{$qSetKey}{$a} <=> $topNvalues{$qSetKey}{$b}} keys %{$topNvalues{$qSetKey}}) {
###											if ($newValue > $topNvalues{$qSetKey}{$dataKey}) {
###												$excludedData{"$qSetKey:#:$dataKey"}=$topNvalues{$qSetKey}{$dataKey}; # exclude previous value
###												delete $topNvalues{$qSetKey}{$dataKey};
###												$topNvalues{$qSetKey}{"$pepSeq:#:$vmod:#:$charge:#:$dataSrc"}=$newValue;
###											}
###											else {
###												$excludedData{"$qSetKey:#:$pepSeq:#:$vmod:#:$charge:#:$dataSrc"}=$newValue; # exclude new value
###											}
###											last; # Only 1 test is necessary because of sort{$a<=>$b}
###										}
###									}
###									else {
###										$topNvalues{$qSetKey}{"$pepSeq:#:$vmod:#:$charge:#:$dataSrc"}=$newValue;
###									}
###								}
###							}
###						}
###					}
###				}
###			}
###		}
###
###		##>Update usableData
###		foreach my $dataKey (keys %excludedData) {
###			my ($qSet,$observationSet,$pepSeq,$vmod,$charge,$dataSrc)=split(':#:',$dataKey);
###			delete $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet};
###			unless (scalar keys %{$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}}) {
###				delete $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet};
###				unless (scalar keys %{$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}}) {
###					delete $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc};
###					unless (scalar keys %{$intensity{"$pepSeq:$vmod:$charge"}}) {
###						delete $intensity{"$pepSeq:$vmod:$charge"};
###					}
###delete $usableData{$pepSeq}{$vmod}{$charge}{$dataSrc};
###unless (scalar keys %{$usableData{$pepSeq}{$vmod}{$charge}}) {
###	delete $usableData{$pepSeq}{$vmod}{$charge};
###	unless (scalar keys %{$usableData{$pepSeq}{$vmod}}) {
###		delete $usableData{$pepSeq}{$vmod};
###		unless (scalar keys %{$usableData{$pepSeq}}) {
###			delete $usableData{$pepSeq};
###		}
###	}
###}
###				}
###			}
###		}
###	}


						$qSetUsage{'USED'}{$qSet}=1;

						#my $pepIdStrg=join(';',@{$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}{$qSet}});
						my %pepIsComplete; # records in which cond peptide has value

						if ($ratioType eq 'Ratio') {
							#my $pepIdStrg=join(';',@{$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}{$qSet}{'ALL_OBS'}});
my $pepIdStrg='';
foreach my $observationSet (@quantiOrder) {
	if ($proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}{$qSet}{$observationSet}) {
		$pepIdStrg.=';' if $pepIdStrg;
		$pepIdStrg.=join(';',@{$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}{$qSet}{$observationSet}});
	}
}
							print DATA "\n$modProtID\t$pepSeq\t$usedVmod\t$charge\t$pepIdStrg\t$proteinAlias{$protID}\t1"; #$validProt
							foreach my $observationSet (@quantiOrder) { # order matters IMPORTANT: only 1 QSET per line for SILAC even with replicates to keep "linked" observations ratios computation (a/c)+(b/d) instead of (a+b)/(c+d) <- later good only for label-free
								#if ($visInReplicate{$protID}{$observationSet} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}) {#}
								if ($intensity{"$pepSeq:$vmod:$charge"}{$dataSrc} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}) {
									print DATA "\t",$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0];
									$pepIsComplete{$obs2Cond{$observationSet}}{$observationSet} = 1;
									#$numPepInCond{'NA'}{$obs2Cond{$observationSet}}++;
									#$numPepInCond4ModProt{$modProtID}{'NA'}{$obs2Cond{$observationSet}}++ if $quantifiedModifID;
									if ($quantifiedModifID) {$numPepInCond4ModProt{$modProtID}{'NA'}{$obs2Cond{$observationSet}}{"$pepSeq:$usedVmod:$charge:$infRatioDataSrc"}=$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet};}
									else {$numPepInCond{'NA'}{$obs2Cond{$observationSet}}{"$pepSeq:$usedVmod:$charge:$infRatioDataSrc"}=$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet};}
								}
								else {
									print DATA "\tNA";
									$pepIsComplete{$obs2Cond{$observationSet}}{$observationSet} = 0;
								}
								$numPepValues++;
							}
							# checking if peptide has a value in all conditions
							my $pepIsOk = 1;
							while (my ($condition, $completeChan) = each %pepIsComplete) {
								unless (scalar grep { $_ == 1 } values %{$completeChan}) {
									$pepIsOk = 0;
									last;
								}
							}
							$countCompletePeptides += $pepIsOk;
						}
						else { # (Super/Simple)Ratio
my ($anaID)=$qSet=~/^A(\d+)/;
							foreach my $observationSet (@quantiOrder) {
								my $itraqStrg=($labeling=~/ITRAQ|TMT/)? $observationInfo{$observationSet}[1] : ''; # Fake SILAC for ITRAQ: ..._repX (PP Tibor set1)
								#if ($visInReplicate{$protID}{$observationSet} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}) {#}
								if ($intensity{"$pepSeq:$vmod:$charge"}{$dataSrc} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}) {
#if (!$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}{$qSet}{$observationSet}) {
#	print DEBUG "OBSSET\t$modProtID\t$pepSeq${usedVmod}_${charge}_$dataSrc\t$observationInfo{$observationSet}[0]\tNA\t$observationInfo{$observationSet}[1]\t$observationInfo{$observationSet}[2]\t$observationInfo{$observationSet}[3]\t$qSet$itraqStrg\tXXXX\t$proteinAlias{$protID}\t1\t",$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0],"\n";
#	next;
#}
next unless $proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}{$qSet}{$observationSet}; # TODO: Sometimes is not defined for last key. why?
#----- NEW EXPERIMENT assignment (31/10/17) ----->
my ($anaID)=$qSet=~/^A(\d+)/;
my $expCode=($labeling eq 'FREE')? 'A' : $ana2Experiment{$anaID}; #$techRep2Experiment{$observationSet};
##my $expCode='A'; # default. for LABELFREE
##if ($Rdesign eq 'SUPERRATIO') {$expCode=$observationInfo{$observationSet}[3];}
##elsif ($Rdesign eq 'LABELED') {
##	my $pepID=$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}{$qSet}{$observationSet}[0];
##	$expCode=$peptideID2Experiment{$pepID} || 'A'; # just to be safe
##}

									#$pepIsComplete{$obs2Cond{$observationSet}}{$observationSet} = 1;
									my $pepIdStrg=join('=',@{$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}{$qSet}{$observationSet}}); # used to be ';'
									#print DATA "\n$protID\t$pepSeq${vmod}_${charge}_$dataSrc\t$observationInfo{$observationSet}[0]\tNA\t$observationInfo{$observationSet}[1]\tA$observationInfo{$observationSet}[3]_Q$qSet\t$pepIdStrg\t$proteinAlias{$protID}\t1\t",$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0];
									print DATA "\n$modProtID\t$pepSeq${usedVmod}_${charge}_$dataSrc\t$observationInfo{$observationSet}[0]\tNA\t$observationInfo{$observationSet}[1]\t$observationInfo{$observationSet}[2]\t$expCode\t$qSet$itraqStrg\t$pepIdStrg\t$proteinAlias{$protID}\t1\t",$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0]; #,"\t$usablePeptide";
									#if ($usablePeptide) {
									#$numPepInCond{$observationInfo{$observationSet}[3]}{$obs2Cond{$observationSet}}++; # numPepInCond{Exp}{condID} Exp needed because of SuperSILAC
									#$numPepInCond4ModProt{$modProtID}{$observationInfo{$observationSet}[3]}{$obs2Cond{$observationSet}}++ if $quantifiedModifID;

#----- INFINITE RATIOS
#my $usedExpCode=($ratioType eq 'SuperRatio')? $observationInfo{$observationSet}[3] : 'NA'; # relevant for SuperRatio only
#									if ($quantifiedModifID) {$numPepInCond4ModProt{$modProtID}{$usedExpCode}{$obs2Cond{$observationSet}}{"$pepSeq:$usedVmod:$charge:$infRatioDataSrc"}=$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0];}
#									else {
#										$numPepInCond{$usedExpCode}{$obs2Cond{$observationSet}}{"$pepSeq:$usedVmod:$charge:$infRatioDataSrc"}=$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0]; # numPepInCond{Exp}{condID} Exp needed because of SuperSILAC
#									}
# |
# |
# V
# Ratios are no longer Experiment-specific in SuperRatio => cannot use $usedExpCode as synonyme of ratio
my %contextRatios;
if ($ratioType eq 'SuperRatio' || $Rdesign eq 'PEP_RATIO') { # WARNING: Labeled w/ PEP_INTENSITY are not treated the same way (if more than 2 conds)
	if ($observationInfo{$observationSet}[0] eq 'State1') { # state is ref -> find matching state(s) given current analysis
		foreach my $i (1..$#conditionList) {
			if ($anaConds{$anaID}{$conditionList[$i]}) { # anaID covers this ratio (format <condIDy>/<condIDx>)
				$contextRatios{$conditionList[$i].'/'.$conditionList[0]}=1;
				push @{$anaState1RatioPos{$anaID}},$i; # $i<=> primary ratioPos
			}
		}
	}
	else { # State is test -> only 1 ref (State1). if SuperRatio 2ndary ratios are not included here
		my ($statePos)=$observationInfo{$observationSet}[0]=~/State(\d+)/;
		$contextRatios{$conditionList[$statePos-1].'/'.$conditionList[0]}=1; # anaID covers this ratio (format <condIDy>/<condIDx>)
	}
}
else {
	$contextRatios{'NA'}=1;
}
foreach my $context (keys %contextRatios) {
									if ($quantifiedModifID) {$numPepInCond4ModProt{$modProtID}{$context}{$obs2Cond{$observationSet}}{"$pepSeq:$usedVmod:$charge:$infRatioDataSrc"}=$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0];}
									else {
										$numPepInCond{$context}{$obs2Cond{$observationSet}}{"$pepSeq:$usedVmod:$charge:$infRatioDataSrc"}=$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0]; # numPepInCond{Exp}{condID} Exp needed because of SuperSILAC
									}
}

									$numPepValues++;
									#}
								}
								#else {
								#	$pepIsComplete{$obs2Cond{$observationSet}}{$observationSet} = 0;
								#}
							}
# TODO: code below should be more efficient (19/02/15)
###if ($intensity{"$pepSeq:$vmod:$charge"}{$dataSrc} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}) {
###	foreach my $observationSet (keys %{$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}}) {
###		#$pepIsComplete{$obs2Cond{$observationSet}}{$observationSet} = 1;
###		my $pepIdStrg=join(';',@{$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}{$qSet}{$observationSet}});
###		#print DATA "\n$protID\t$pepSeq${vmod}_${charge}_$dataSrc\t$observationInfo{$observationSet}[0]\tNA\t$observationInfo{$observationSet}[1]\tA$observationInfo{$observationSet}[3]_Q$qSet\t$pepIdStrg\t$proteinAlias{$protID}\t1\t",$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0];
###		print DATA "\n$modProtID\t$pepSeq${vmod}_${charge}_$dataSrc\t$observationInfo{$observationSet}[0]\tNA\t$observationInfo{$observationSet}[1]\t$observationInfo{$observationSet}[2]\t$observationInfo{$observationSet}[3]\t$qSet\t$pepIdStrg\t$proteinAlias{$protID}\t1\t",$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0]; #,"\t$usablePeptide";
###		#if ($usablePeptide) {
###	#$numPepInCond{$observationInfo{$observationSet}[3]}{$obs2Cond{$observationSet}}++; # numPepInCond{Exp}{condID} Exp needed because of SuperSILAC
###	#$numPepInCond4ModProt{$modProtID}{$observationInfo{$observationSet}[3]}{$obs2Cond{$observationSet}}++ if $quantifiedModifID;
###		$numPepInCond{$observationInfo{$observationSet}[3]}{$obs2Cond{$observationSet}}{"$pepSeq:$vmod:$charge:$dataSrc"}=1; # numPepInCond{Exp}{condID} Exp needed because of SuperSILAC
###		$numPepInCond4ModProt{$modProtID}{$observationInfo{$observationSet}[3]}{$obs2Cond{$observationSet}}{"$pepSeq:$vmod:$charge:$dataSrc"}=1 if $quantifiedModifID;
###		$numPepValues++;
###	}
###}
							#(PP 02/10/15) last if $labeling ne 'FREE'; # data from only 1 quantifSet allowed for each pep+seq+charge+src (true for SILAC before 10/2015)
						}

						## checking if peptide has a value in all conditions
						#my $pepIsOk = 1;
						#while (my ($condition, $completeChan) = each %pepIsComplete) {
						#	unless (scalar grep { $_ == 1 } values %{$completeChan}) {
						#		$pepIsOk = 0;
						#		last;
						#	}
						#}
						#$countCompletePeptides += $pepIsOk;
					}
				}
			}
			#&checkInfiniteRatios($modProtID,\%numPepInCond4ModProt) if $quantifiedModifID;
		}
	}

	##>New +/-inf ratio decision for each cond pairs ********!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	# infinite: if ratio >4 for 1,2 pep/protein or ratio >2 for 3+ pep/protein
	if ($quantifParameters{'DB'}{'MIN_INF_RATIOS'}[0]==0) { ##### $ratioType=~/S\w+Ratio/ &&
		if ($quantifiedModifID) {
			foreach my $modProtID (keys %numPepInCond4ModProt) {
				&checkInfiniteRatios($modProtID,$numPepInCond4ModProt{$modProtID});
			}
		}
		else {
			&checkInfiniteRatios($protID,\%numPepInCond);
		}
	}
}
close DATA;
#close DEBUG;
#die "Test table.txt is complete!!!";

####> Desactivate FDR option if the number of proteins sent to R is below a threshold : 1 (2 proteins needed)
#my $protNumber=scalar keys %prot2quantify; #%proteins;
#if ($quantifParameters{'R'}{'pAdj'}[0] eq 'TRUE' && $protNumber < $minProtNumber) {
if ($quantifParameters{'R'}{'pAdj'} && $quantifParameters{'R'}{'pAdj'}[0] eq 'TRUE' && $protNumber == 1) { # only for old algo
	@{$quantifParameters{'R'}{'pAdj'}}=('FALSE');
	delete($quantifParameters{'R'}{'pAdj.method'});
	my $newQuantiAnnotation='';
	foreach my $annoStrg (split (/::/,$quantiAnnotation)) {
		if ($annoStrg =~ /^FDR\_/ ) {
			$newQuantiAnnotation.='FDR_CONTROL=FALSE::' if ($annoStrg =~ /^FDR_CONTROL/);
		}
		else{
			$newQuantiAnnotation.=$annoStrg."::";
		}
	}
	chop($newQuantiAnnotation);
	chop($newQuantiAnnotation);
	$dbh->do("UPDATE QUANTIFICATION SET QUANTIF_ANNOT='$newQuantiAnnotation' WHERE ID_QUANTIFICATION=$quantifID");
	$dbh->commit;
}

# Checking if at least 1 complete peptide (no NA values)
#if ($ratioType eq 'Ratio' && $countCompletePeptides == 0) {
#	$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID"); # Failed
#	$dbh->commit;
#	my $hintStrg=($labeling eq 'SILAC' && $call eq 'DESIGN')? "\nMay be missmatch(es) between Design and Analysis labeling modification(s)?" : '';
#	die "There is no peptide with values in all conditions.$hintStrg";
#}

&promsQuantif::writeQuantifParameterFiles($dataDir,$quantifParameters{'R'});

&checkForErrors;

$dbh->disconnect;
#die "Data files are ready!!!"; # DEBUG


								############################################
			#############################> RUNNING R SCRIPTS <##########################
								############################################

####> 4th: Launch R script
open(FILESTAT,">>$fileStat");
print FILESTAT "2/3 Running quantification\n";
close FILESTAT;

my %cluster=&promsConfig::getClusterInfo; #('debian') # default is 'centos'
my $pathR=($cluster{'on'})? $cluster{'path'}{'R'} : $promsPath{'R'};

my $RcommandString='';
if ($ratioType eq 'Ratio') {
	$RcommandString="export LANG=en_US.UTF-8; cd $runDir; $pathR/R CMD BATCH --no-save --no-restore $promsPath{R_scripts}/LabeledAndUnlabeledQuanti.Function.R";
}
else { #  (Super/Simple)Ratio LabelFree
	#foreach my $scriptName ('AnalysisDiffLimma.R','FunctionLimma.R','AffectDefaultValuesLimma.R','AffectParametersLimma.R') {
	#	symlink("$promsPath{R_scripts}/$scriptName","$runDir/$scriptName");
	#}
	open(R_SCRIPT,">$runDir/AnalysisDiffLimma.R");
	print R_SCRIPT qq
|
###################################################################################
# Launcher for quantification R scripts (provides path for sourcing dependencies) #
###################################################################################

filepath="$promsPath{R_scripts}/"
source(paste(filepath,"AnalysisDiffLimma.R",sep=""))
|;
	close R_SCRIPT;
	$RcommandString="export LANG=en_US.UTF-8; cd $runDir; $pathR/R CMD BATCH --no-save --no-restore AnalysisDiffLimma.R";
}

### ALWAYS SKIP CLUSTER (PP 17/09/18)!!! ###
###if ($cluster{'on'}) {
###	my $bashFile = "$runDir/runXICProtQuantification.sh";
###	my $clusterRcommandString=$cluster{'buildCommand'}->($runDir,$RcommandString);
###	my $maxHours=int(0.5+($numPepValues/25000)); $maxHours=2 if $maxHours < 2; $maxHours=48 if $maxHours > 48; # 1 M lines -> 40 h
###	#my $maxMem=($numPepValues < 250000)? int(0.5+($numPepValues/25000)) : int(0.5+($numPepValues/8000)); $maxMem=1 if $maxMem < 1; $maxMem.='gb'; # 1 M lines -> 80 gb
###	#my $maxMem='2gb';
###	# RAM = 8E-7 * <num lines> + 0.1614
###	#my $maxMem=0.5 + 8E-7 * $numPepValues;
###	#my $maxMem=0.5 + 1E-6 * $numPepValues;
###	my $maxMem=int(1.5 + 1E-6 * $numPepValues);
###	$maxMem=50 if $maxMem > 50;
###	$maxMem.='Gb';
###	open (BASH,">$bashFile");
###	print BASH qq
###|#!/bin/bash
#####resources
####PBS -l mem=$maxMem
####PBS -l nodes=1:ppn=1
####PBS -l walltime=$maxHours:00:00
####PBS -q batch
#####Information
####PBS -N myProMS_protQuant_$quantifID
####PBS -m abe
####PBS -o $runDir/PBS.txt
####PBS -e $runDir/PBSerror.txt
####PBS -d $runDir
###
##### Command
###$clusterRcommandString
###echo _END_$quantifID
###|;
####close DEBUG;
###	close BASH;
###	my $modBash=0775;
###	chmod $modBash, $bashFile;
###	#system "$promsPath{qsub}/qsub $bashFile > $runDir/torqueID.txt";
###	#system "$cluster{command} $bashFile > $runDir/torqueID.txt 2> $runDir/connection_error.txt";
###	#system "ssh -q -o \"UserKnownHostsFile=/dev/null\" -o \"StrictHostKeyChecking no\" calcsub \"qsub $bashFile > $runDir/torqueID.txt\" 2>$runDir/ssh_error.txt ";
###	$cluster{'sendToCluster'}->($bashFile); # bash file, qsub output file[, ssh error file (only for centos cluster)]
###	sleep 30;
###
###	###>Waiting for R job to run
###	my $pbsError;
###	my $nbWhile=0;
###	my $maxNbWhile=$maxHours*60*2;
###	while ((!-e "$runDir/PBS.txt" || !`tail -3 $runDir/PBS.txt | grep _END_$quantifID`) && !$pbsError) {
###		if ($nbWhile > $maxNbWhile) {
###			$dbh=&promsConfig::dbConnect('no_user'); # reconnect
###			$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID"); # Failed
###			$dbh->commit;
###			$dbh->disconnect;
###			die "Aborting quantification: R is taking too long or died before completion";
###		}
###		sleep 30;
###		#$pbsError=`head -5 $runDir/PBSerror.txt` if -e "$runDir/PBSerror.txt";
###		$pbsError=$cluster{'checkError'}->("$runDir/PBSerror.txt");
###		$nbWhile++;
###	}
###}
###else { ###>Run job on Web server
	system $RcommandString;
###}
sleep 3;

$dbh=&promsConfig::dbConnect('no_user'); # reconnect

####>ERROR Management<####
my $RoutFile=($ratioType eq 'Ratio')? 'LabeledAndUnlabeledQuanti.Function.Rout' : 'AnalysisDiffLimma.Rout';
my $RoutStrg=`tail -3 $runDir/$RoutFile`;
unless ($RoutStrg=~/proc\.time\(\)/) {
	$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID");
	$dbh->commit;
	$dbh->disconnect;
	$RoutStrg=`tail -20 $runDir/$RoutFile`;
	my $RerrorStrg="R script has generated an error!";
	my $inError=0;
	foreach my $line (split(/\n/,$RoutStrg)) {
		next if (!$inError && $line !~ /^Error in/); # skip everything before "Error in..."
		$inError=1;
		$RerrorStrg.="\n$line";
	}
	die $RerrorStrg;
}

#die "R quantification is finished!!!"; # DEBUG


									############################################
			#############################> PARSING QUANTIFICATION RESULTS <##########################
									############################################

open(FILESTAT,">>$fileStat");
print FILESTAT "3/3 Parsing results\n";
close FILESTAT;
#$dbh->disconnect; exit; # DEBUG

my $numQuantifRatios=scalar @ratioTarget; # =numConditions*0.5*(numConditions-1)
my $adjStrg=($ratioType=~/S\w+Ratio/ || $quantifParameters{'DB'}{'FDR_CONTROL'}[0] eq 'TRUE')? '_ADJ' : '';

####>Fetching list of quantification parameters<####
my %quantifParamIDs;
my $sthQP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,P.CODE FROM QUANTIFICATION_PARAMETER P,QUANTIFICATION_METHOD M WHERE P.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND M.ID_QUANTIFICATION_METHOD=(SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID)");
$sthQP->execute;
while (my ($paramID,$code)=$sthQP->fetchrow_array) {
    $quantifParamIDs{$code}=$paramID;
}
$sthQP->finish;

my ($quantifMethod)=$dbh->selectrow_array("SELECT CODE FROM QUANTIFICATION_METHOD WHERE ID_QUANTIFICATION_METHOD=(SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID)");
my $sthInsProt=$dbh->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_QUANTIFICATION,ID_PROTEIN,ID_QUANTIF_PARAMETER,QUANTIF_VALUE,TARGET_POS) VALUES ($quantifID,?,?,?,?)"); # PK autoincrement
my ($sthInsModRes,$sthInsProtRes);
if ($quantifiedModifID) {
	$sthInsModRes=$dbh->prepare("INSERT INTO MODIFIED_RESIDUE (ID_QUANTIFICATION,RESIDUE,POSITION) VALUES ($quantifID,?,?)");
	$sthInsProtRes=$dbh->prepare("INSERT INTO PROTQUANTIF_MODRES (ID_MODIF_RES,ID_PROT_QUANTIF) VALUES (?,?)");
}



################################################################
####>Recording protein ratios (& p-values, Std. dev if any)<####
################################################################
my (%protQuantified,%modResiduesID);

####>Recording infinite ratios ratios<#### because prot might not be kept in quantif result files
if ($quantifiedModifID) { # modif quantification
	foreach my $modProtID (sort{&promsMod::sortSmart($a,$b)} keys %infiniteRatioProteins) {
		my ($protID,@modResidues)=split(/[-.]/,$modProtID);
		&promsQuantif::insertModifiedResidues($protID,\@modResidues,\%modResiduesID,$dbh,$sthInsModRes);
		foreach my $ratioPos (sort{$a<=>$b} keys %{$infiniteRatioProteins{$modProtID}}){
			$sthInsProt->execute($protID,$quantifParamIDs{'RATIO'},$infiniteRatioProteins{$modProtID}{$ratioPos},$ratioPos);
			&promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF'));
			$protQuantified{$modProtID}{$ratioPos}=1;
		}
	}
}
else { # whole protein quantif
	foreach my $protID (sort{$a<=>$b} keys %infiniteRatioProteins) {
		foreach my $ratioPos (sort{$a<=>$b} keys %{$infiniteRatioProteins{$protID}}){
			$sthInsProt->execute($protID,$quantifParamIDs{'RATIO'},$infiniteRatioProteins{$protID}{$ratioPos},$ratioPos);
			$protQuantified{$protID}{$ratioPos}=1;
		}
	}
}


###################################>(Super/Simple)Ratio (ML->AS scripts)<#############################################
my $numLostProteins=0;
if ($ratioType=~/S\w+Ratio/) {
	#unlink "$runDir/Rplots.pdf" if -e "$runDir/Rplots.pdf"; # temporary ?!!!!!!!!!!!!!!!!!!!!!!!!!!

	###>Generating ratio codes used by R<###
	my (%measurePos2Code,%measureCode2RatioPos,%state2RatioPos,%experiment2RatioPos); # a state can be involved in multiple ratios
	my $rPos=0;
####################### Old
	####foreach my $s1 (1..$numConditions-1) {
	####	foreach my $s2 ($s1+1..$numConditions) {
	####		$rPos++;
	####		$measurePos2Code{$rPos}="State$s2.State$s1";
	####		$measureCode2RatioPos{"State$s2.State$s1"}=$rPos;
	####		push @{$state2RatioPos{"State$s1"}},$rPos;
	####		push @{$state2RatioPos{"State$s2"}},$rPos;
	####	}
	####	last if ($quantifParameters{'DB'}{'SINGLE_REF'} || $Rdesign eq 'SUPERRATIO');
	####}
	####if ($Rdesign eq 'SUPERRATIO') {
	####	foreach my $s1 (2..$numConditions-1) {
	####		foreach my $s2 ($s1+1..$numConditions) {
	####			$rPos++;
	####			$measurePos2Code{$rPos}="State$s2.State1.State$s1.State1";
	####		}
	####	}
	####}
	####elsif ($Rdesign eq 'LABELFREE') { # also record state mean
	####	foreach my $s (1..$numConditions) {
	####		$rPos++;
	####		$measurePos2Code{$rPos}="State$s";
	####	}
	####}
##################### New
##if ($Rdesign eq 'SUPERRATIO') { # 2ndary ratios only (no primary ratios recorded!!!)
##	foreach my $s1 (2..$numConditions-1) { # skip State1 -> 2ndary ratios rely on Experiments (B/A,C/A,D/A,...,C/B,D/B,...,C/D,...)
##		foreach my $s2 ($s1+1..$numConditions) {
##			$rPos++;
##			$measurePos2Code{$rPos}="State$s2.State1.State$s1.State1"; # for ResultsDAProt.txt (used for SuperRatio only!!!)
##			push @{$experiment2RatioPos{chr(63+$s1)}},$rPos; # 63+2=65->A,... for resultsPep.txt & table.txt
##			push @{$experiment2RatioPos{chr(63+$s2)}},$rPos; # for resultsPep.txt & table.txt
##		}
##	}
##}
##else { # LABELED & LABELFREE
	foreach my $s1 (1..$numConditions-1) { # Primary ratios
		foreach my $s2 ($s1+1..$numConditions) {
			$rPos++;
			$measurePos2Code{$rPos}="State$s2.State$s1"; # for ResultsDAProt.txt
			if ($Rdesign eq 'PEP_INTENSITY') {
				push @{$measureCode2RatioPos{"State$s1"}},$rPos; # for resultsPep.txt (multiple rPos: when is ref & when is test state!)
				push @{$measureCode2RatioPos{"State$s2"}},$rPos;
			}
			else {@{$measureCode2RatioPos{"State$s2.State$s1"}}=($rPos);} # for resultsPep.txt
			push @{$state2RatioPos{"State$s1"}},$rPos if ($ratioType eq 'SimpleRatio' || $s1 > 1); # for inf ratios in table.txt. If SuperRatio: Skip primary ratios for State1 (recorded in %anaState1RatioPos because State1 ratios are analysis-dependant)
			push @{$state2RatioPos{"State$s2"}},$rPos;
		}
		last if ($quantifParameters{'DB'}{'SINGLE_REF'} || $ratioType eq 'SuperRatio');
	}
	if ($ratioType eq 'SuperRatio') { # 2ndary ratios
		foreach my $s1 (2..$numConditions-1) { # skip State1 -> 2ndary ratios rely on Experiments (B/A,C/A,D/A,...,C/B,D/B,...,C/D,...)
			foreach my $s2 ($s1+1..$numConditions) {
				$rPos++;
				$measurePos2Code{$rPos}="State$s2.State1.State$s1.State1"; # for ResultsDAProt.txt (used for SuperRatio only!!!)
				#push @{$experiment2RatioPos{chr(63+$s1)}},$rPos; # 63+2=65->A,... for resultsPep.txt & table.txt
				#push @{$experiment2RatioPos{chr(63+$s2)}},$rPos; # for resultsPep.txt & table.txt
				push @{$measureCode2RatioPos{"State$s1.State1"}},$rPos; # (norm/primary)ratio used in multiple secondary ratios
				push @{$measureCode2RatioPos{"State$s2.State1"}},$rPos; # (norm/primary)ratio used in multiple secondary ratios
				#push @{$state2RatioPos{"State1"}},$rPos; # Include all 2ndary ratio for State1
				push @{$state2RatioPos{"State$s1"}},$rPos; # for inf ratios in table.txt
				push @{$state2RatioPos{"State$s2"}},$rPos;
			}
		}
	}
	elsif ($Rdesign eq 'PEP_INTENSITY') { # also record state mean
		foreach my $s (1..$numConditions) {
			$rPos++;
			$measurePos2Code{$rPos}='State'.$s; # for ResultsDAProt.txt
		}
	}

#print DEBUG ">>measureCode2RatioPos\n";
#foreach my $state (keys %measureCode2RatioPos) {
#	print DEBUG ">$state: '@{$measureCode2RatioPos{$state}}'\n";
#}
##}

	#my $log2=log(2);

	####>Extracting ratio, (mean), p-values, std dev & conf intervals<####
	my @RDbCodeList=(['pVal','PVAL_ADJ'],['SD','SD_GEO'],['CI2.5','CONF_LOW'],['CI97.5','CONF_UP']);
	my %protHeader2Idx;
	open(PROT,"$resultDir/ResultsDAProt.txt");
	while (<PROT>) {
		chomp;
		s/"//g;
		my @lineData=split(/\t/,$_);
		####>Header------------------
		if ($.==1) { # 1st line of the file
			my $colIdx=0;
			foreach my $colName (@lineData) {
				$protHeader2Idx{$colName}=$colIdx;
				$colIdx++;
			}
			next;
		}
		####>Data------------------
		my $protID0=$lineData[0]; # $lineData[$protHeader2Idx{'proteinId'}];
		next if !$protID0; # || $protID=~/\D/ just to be safe
		next if ($quantifiedModifID && $protID0 !~/-/); # must be a mod prot
		my ($protID,@modResidues)=($quantifiedModifID)? split(/[-.]/,$protID0) : ($protID0,undef);
		&promsQuantif::insertModifiedResidues($protID,\@modResidues,\%modResiduesID,$dbh,$sthInsModRes) if $quantifiedModifID; # <-- %modResiduesID should be full at <PROT1>

		###>Conditions ratio data
		foreach my $targetPos (sort{$a<=>$b} keys %measurePos2Code) { #(1..$numQuantifRatios) {
			next if ($infiniteRatioProteins{$protID0} && $infiniteRatioProteins{$protID0}{$targetPos});
			#next if ($protQuantified{$protID0} && $protQuantified{$protID0}{$ratioPos});
			#>Fold change
			if ($targetPos <= $numQuantifRatios) { # foldChange
				my $log2FC=$lineData[$protHeader2Idx{'Log2FC_'.$measurePos2Code{$targetPos}}];
				next if $log2FC=~/NA/i; # NA or NaN
				my $fcValue=($log2FC eq '-Inf')? 0.001 : ($log2FC eq 'Inf')? 1000 : 2**$log2FC; #exp($log2FC*$log2);
				$sthInsProt->execute($protID,$quantifParamIDs{'RATIO'},$fcValue,$targetPos);
				if ($quantifiedModifID) {
					my $protQuantifID=$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF');
					#&promsQuantif::insertModifiedResidues($protID,\@modResidues,\%modResiduesID,$dbh,$sthInsModRes);
					&promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$protQuantifID);
				}
				$protQuantified{$protID0}{$targetPos}=1;
			}
			else { # mean state (LABELFREE)
next unless defined $protHeader2Idx{'Log2Mean_'.$measurePos2Code{$targetPos}}; # temp AS should had this again soon (PP 19/10/17)
				my $log2Mean=$lineData[$protHeader2Idx{'Log2Mean_'.$measurePos2Code{$targetPos}}];
				next if $log2Mean=~/NA/i; # NA or NaN
				$sthInsProt->execute($protID,$quantifParamIDs{'MEAN_STATE'},2**$log2Mean,$targetPos); #  exp($log2Mean*$log2)
				&promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $quantifiedModifID; # unlikely in LABELFREE...
			}
			#>Associated parameters
			foreach my $refParam (@RDbCodeList) {
				my ($paramR,$paramDB)=@{$refParam};
				my $paramValue=$lineData[$protHeader2Idx{$paramR.'_'.$measurePos2Code{$targetPos}}];
				next if $paramValue=~/NA/i; # NA or NaN
				$sthInsProt->execute($protID,$quantifParamIDs{$paramDB},$paramValue,$targetPos);
				&promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $quantifiedModifID;
			}
		} # next ratioPos

	} # next line
	close PROT;

	####>Computing peptide counts (NUM_PEP_TOTAL, NUM_PEP_USED & DIST_PEP_USED)<####
	my (%numPeptideSets,%usedPeptideSets,%distinctPeptideSets); # Set <=> num label channels (1 for Label-Free)

	##>Step 1: Read resultsPep.txt for NUM_PEP_USED & DIST_PEP_USED
	my %pepHeader2Idx;
	open(PEP,"$resultDir/resultsPep.txt");
	while(<PEP>) {
		#next if $.==1;
		#s/"//g; # remove all "
		chomp;
		#my ($condCode,$protID0,$peptideCharge,$experiment,$measPep,$usedPepIdStrg)=split(/\t/,$_);
		#my ($condCode,$experiment,$bioRep,$techRep,$protID0,$peptideCharge,$measPep,$pepIdStrg,$outlier)=split(/\t/,$_); # v3a
		#my ($protID0,$peptideCharge,$experiment,$bioRep,$techRep,$condCode,$protName,$pepIdStrg,$protValid,$measPep,$normProt,$outlier)=split(/\t/,$_); # v3b LF
		#ProteinID	Peptide	Experiment	Condition	proteinValidity	log2Measure	replicate	repTech	proteinName	PeptideId	normProtein	out # v3b LABELED
		my @lineData=split(/\t/,$_);
		####>Header------------------
		if ($.==1) { # 1st line of the file
			my $colIdx=0;
			foreach my $colName (@lineData) {
				$pepHeader2Idx{$colName}=$colIdx;
				$colIdx++;
			}
			next;
		}
		####>Data------------------
		next if $lineData[$pepHeader2Idx{'out'}] ne 'NA'; # skip outliers
		my $condCode=$lineData[$pepHeader2Idx{'Condition'}];
		next if ($Rdesign eq 'PEP_RATIO' && !$measureCode2RatioPos{$condCode}); # ratio not used (SINGLE_REF) StateY.StateX (not applicable for LABELFREE: StateX\nStateY)
		my ($seqVarMod,$charge)=split(/_/,$lineData[$pepHeader2Idx{'Peptide'}]);
		#my $numPepSets=scalar (split(';',$pepIdStrg));
		my $numPepSets=1; # default
		# TODO: True for SILAC & other Isotopes (not LF or isobar) --->
		if ($labeling eq 'SILAC' && $lineData[$pepHeader2Idx{'PeptideId'}]=~/\+/) { # XIC were summed
			my %distinctQsets;
			foreach my $pepID (split(/[\.\+\=]/,$lineData[$pepHeader2Idx{'PeptideId'}])) {
				$distinctQsets{$pepID2QuantifSetSILAC{$pepID}}=1;
			}
			$numPepSets=scalar keys %distinctQsets;
		}
#my $refRatioList;
#if ($Rdesign eq 'LABELFREE') {$refRatioList=$state2RatioPos{$condCode}} # StateX: multiple primary ratios
##elsif ($Rdesign eq 'LABELED') {$refRatioList=$measureCode2RatioPos{$condCode};} # StateY.StateX: only 1 primary ratio
##else {$refRatioList=$experiment2RatioPos{$lineData[$pepHeader2Idx{'Experiment'}]};} # SUPERRATIO, ##Old## Exp: multiple secondary ratios
		my $protID0=$lineData[$pepHeader2Idx{'ProteinID'}];
		next if !$protQuantified{$protID0};
		foreach my $ratioPos (@{$measureCode2RatioPos{$condCode}}) {
			next if !$protQuantified{$protID0}{$ratioPos};
			next if ($infiniteRatioProteins{$protID0} && $infiniteRatioProteins{$protID0}{$ratioPos}); # computed later for inf ratios
			$distinctPeptideSets{$protID0}{$ratioPos}{$seqVarMod}=1; # DIST_PEP_USED
			$usedPeptideSets{$protID0}{$ratioPos}+=$numPepSets;
			#$usedPeptideSets{$protID0}{$ratioPos}++;
		}
	}
	close PEP;

	##>Step 2: Read table.txt for NUM_PEP_TOTAL & inf ratios NUM_PEP_USED, DIST_PEP_USED
	####my %expCode2RatioPos; # used for SuperRatio only!!!
	####if ($Rdesign eq 'SUPERRATIO') {
	####	foreach my $primRatioPos (1..$numConditions-1) {
	####		$expCode2RatioPos{chr(64+$primRatioPos)}=$primRatioPos;
	####	}
	####}
	my %usedPepSets;
	my %lostProteins; # proteins in table.txt but not quantified at all!!!
	my $prevExp='-';
	my $ratioPos=0;
	open(DATA,"$dataDir/table.txt");
	while(<DATA>) {
		next if $.==1;
		#chomp;
		my ($protID0,$peptideData,$state,$treatment,$bioRep,$techRep,$experiment,$qSet,$pepIdStrg)=split(/\t/,$_);
		unless ($protQuantified{$protID0}) {
			$lostProteins{$protID0}=1;
			next;
		}
		$numPeptideSets{$protID0}{$peptideData}=1; # NUM_PEP_TOTAL (counted ONCE across all states) LABEL-FREE: requires that dataSrc tag != '-'!!!
#foreach my $pepID (split(/\+/),$pepIdStrg) { # v3: true number of peptides used (for all labeling channels)
#	$numPeptideSets{$protID0}{$pepID}=1; # NUM_PEP_TOTAL (counted ONCE across all states)
#}
		next unless $infiniteRatioProteins{$protID0}; # No source/replicate aggregation for NUM_PEP_USED in case of infinite ratio!!!
		#--- Only used for infinite ratios from here on ---->
		my ($seqVarMod,$charge,$src)=split(/_/,$peptideData); # DIST_PEP_USED: $peptideCharge

		####if ($Rdesign eq 'LABELFREE') {
		####	foreach my $rPos (keys %{$statePos2TargetPos{$state}}) { # for MEAN & RATIO
		####		if ($rPos <= $numQuantifRatios && $infiniteRatioProteins{$protID0}{$rPos}) { # ratios only
		####			$usedPepSets{$protID0}{$rPos}{$peptideData}=1; # $peptideCharge
		####			$distinctPeptideSets{$protID0}{$rPos}{$seqVarMod}=1;
		####		}
		####	}
		####}
		####else {
		####	if ($Rdesign eq 'SUPERRATIO') { # "Experiment" carries ratioPos info
		####		if ($experiment ne $prevExp) {
		####			$ratioPos=ord($experiment)-64; # A -> 1, B -> 2, ...
		####			$prevExp=$experiment;
		####		}
		####	}
		####	else { # LABELED. "State" carries ratioPos info (not experiment!): State1,State2,...
		####		if ($state eq 'State1') { # reference state => increment all infinite ratioPos
		####			foreach my $rPos (keys %{$infiniteRatioProteins{$protID0}}) {
		####				next if !$protQuantified{$protID0}{$rPos};
		####				$usedPepSets{$protID0}{$rPos}{$peptideData}=1; # $peptideCharge
		####				$distinctPeptideSets{$protID0}{$rPos}{$seqVarMod}=1;
		####			}
		####			next;
		####		}
		####		else {
		####			($ratioPos)=($state=~/State(\d+)/);
		####			$ratioPos--; # StateN -> ratioPos N-1
		####		}
		####	}
		####	if ($infiniteRatioProteins{$protID0}{$ratioPos}) {
		####		$usedPepSets{$protID0}{$ratioPos}{$peptideData}=1; # $peptideCharge
		####		$distinctPeptideSets{$protID0}{$ratioPos}{$seqVarMod}=1;
		####	}
		####}

		my $refRatioList; #=($Rdesign eq 'SUPERRATIO')? $experiment2RatioPos{$experiment} : $state2RatioPos{$state};
		# SUPERRATIO: State1-based primary ratios are analysis-specific. 2ndary ratio are not
		# LABELED, LABELFREE: "State" carries ratioPos value(s):
		if ($ratioType eq 'SuperRatio') {
			if ($state eq 'State1') {
				my ($anaID)=$qSet=~/^A(\d+)/;
				#my @ratioList=@{$anaState1RatioPos{$anaID}}; # primary ratios using anaID
				#push @ratioList,@{$state2RatioPos{'State1'}}; # all 2ndary ratios
				#>Add pos of ratios using ref State1 and associated test states in context of THIS anaID
				my %myRatioPos;
				foreach my $rPos1 (@{$anaState1RatioPos{$anaID}}) {
					$myRatioPos{$rPos1}=1; # primary ratioPos
					#next if $quantifParameters{'DB'}{'SINGLE_REF'};
					my $testState='State'.($rPos1+1); # in test StateX, X=ratioPos+1 (primary ratios only)
					foreach my $rPos2 (@{$state2RatioPos{$testState}}) { # 2ndary ratioPos
						$myRatioPos{$rPos2}=1;
					}
				}
				my @ratioList=(sort{$a<=>$b} keys %myRatioPos); # hash to revent duplicates
				$refRatioList=\@ratioList;
			}
			else {$refRatioList=$state2RatioPos{$state};}
		}
		else {$refRatioList=$state2RatioPos{$state};}

		foreach my $ratioPos (@{$refRatioList}) {
			next unless $infiniteRatioProteins{$protID0}{$ratioPos};
			$usedPepSets{$protID0}{$ratioPos}{$peptideData}=1; # $peptideCharge
			$distinctPeptideSets{$protID0}{$ratioPos}{$seqVarMod}=1;
		}
	}
	close DATA;
	foreach my $protID0 (keys %usedPepSets) { # inf ratios only
		foreach my $ratioPos (keys %{$usedPepSets{$protID0}}) {
			$usedPeptideSets{$protID0}{$ratioPos}=scalar keys %{$usedPepSets{$protID0}{$ratioPos}};
		}
	}

	##>Step 3 (SuperRatio only): Compute NUM_PEP_USED & DIST_PEP_USED for ratios of ratios (only defined for primary ratios so far!!!)
	if ($ratioType eq 'SuperRatio') {
		my $numPrimRatios=$ratioPos;
		foreach my $rPos1 (1..$numPrimRatios-1) {
			foreach my $rPos2 ($rPos1+1..$numPrimRatios) {
				$ratioPos++;
				#>NUM_PEP_USED
				foreach my $protID0 (keys %usedPeptideSets) {
					#next if ($noSuperRatioProteins{$protID0} && $noSuperRatioProteins{$protID0}{$ratioPos});
					next if (!$protQuantified{$protID0} || !$protQuantified{$protID0}{$ratioPos});
					if ($usedPeptideSets{$protID0}{$rPos1} && $usedPeptideSets{$protID0}{$rPos2}) {
						#$usedPeptideSets{$protID0}{$ratioPos}=$usedPeptideSets{$protID0}{$rPos1}+$usedPeptideSets{$protID0}{$rPos2};
						$usedPeptideSets{$protID0}{$ratioPos}=($usedPeptideSets{$protID0}{$rPos1} < $usedPeptideSets{$protID0}{$rPos2})? $usedPeptideSets{$protID0}{$rPos1} : $usedPeptideSets{$protID0}{$rPos2}; # keep smallest (18/02/15)
					}
				}
				#>DIST_PEP_USED
				foreach my $protID0 (keys %distinctPeptideSets) {
					#next if ($noSuperRatioProteins{$protID0} && $noSuperRatioProteins{$protID0}{$ratioPos});
					next if (!$protQuantified{$protID0} || !$protQuantified{$protID0}{$ratioPos});
					if ($distinctPeptideSets{$protID0}{$rPos1} && $distinctPeptideSets{$protID0}{$rPos2}) {
						#foreach my $rPos ($rPos1,$rPos2) {
						#	foreach my $seqVarMod (keys %{$distinctPeptideSets{$protID0}{$rPos}}) {
						#		$distinctPeptideSets{$protID0}{$ratioPos}{$seqVarMod}=1;
						#	}
						#}
						my ($numPos1,$numPos2)=(scalar keys %{$distinctPeptideSets{$protID0}{$rPos1}},scalar keys %{$distinctPeptideSets{$protID0}{$rPos1}});
						$distinctPeptideSets{$protID0}{$ratioPos}=($numPos1 < $numPos2)? $distinctPeptideSets{$protID0}{$rPos1} : $distinctPeptideSets{$protID0}{$rPos2}; # keep smallest (18/02/15)
					}
				}
			}
		}
	}

	##>Step 4: Insert peptide counts in DB
	#>DIST_PEP_USED
	foreach my $protID0 (keys %distinctPeptideSets) {
		my ($protID,@modResidues)=($quantifiedModifID)? split(/[-.]/,$protID0) : ($protID0,undef);
		foreach my $ratioPos (keys %{$distinctPeptideSets{$protID0}}) {
			next if (!$protQuantified{$protID0} || !$protQuantified{$protID0}{$ratioPos});
			$sthInsProt->execute($protID,$quantifParamIDs{'DIST_PEP_USED'},scalar keys %{$distinctPeptideSets{$protID0}{$ratioPos}},$ratioPos);
			&promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $quantifiedModifID;
		}
	}
	#>NUM_PEP_USED
	foreach my $protID0 (keys %usedPeptideSets) {
		my ($protID,@modResidues)=($quantifiedModifID)? split(/[-.]/,$protID0) : ($protID0,undef);
		foreach my $ratioPos (sort keys %{$usedPeptideSets{$protID0}}) {
			next if (!$protQuantified{$protID0} || !$protQuantified{$protID0}{$ratioPos});
			$sthInsProt->execute($protID,$quantifParamIDs{'NUM_PEP_USED'},$usedPeptideSets{$protID0}{$ratioPos},$ratioPos);
			&promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $quantifiedModifID;
		}
	}
	#>NUM_PEP_TOTAL
	foreach my $protID0 (keys %numPeptideSets) {
		next if !$protQuantified{$protID0};
		my ($protID,@modResidues)=($quantifiedModifID)? split(/[-.]/,$protID0) : ($protID0,undef);
		$sthInsProt->execute($protID,$quantifParamIDs{'NUM_PEP_TOTAL'},scalar keys %{$numPeptideSets{$protID0}},undef); # no ratioPos
		&promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $quantifiedModifID;
	}

	$numLostProteins=scalar keys %lostProteins;
	if ($numLostProteins) {
		open (LOST,">$resultDir/lostProteins.txt");
		print "Protein_ID\tProtein_Name";
		foreach my $modProtID (sort{&promsMod::sortSmart($a,$b)} keys %lostProteins) {
			my ($protID)=$modProtID=~/^(\d+)/;
			print LOST "\n$modProtID\t$proteinAlias{$protID}";
		}
		close LOST;
	}
}

else {
	###################################>Ratio or TnPQ (FC script)<#############################################
	# PP 09/12/14:
	# MAKE TnPQ obsolete
	# Ratio only if 2 conditions (1 ratioPos!)
	###> The same script launches TNPQ or Ratio method so according to what was chosen by the user, the
	###> parsing will not be the same because the R-scripts output will not be the same...

	####>Parsing prot_results.txt
	#my @codeProt1=($quantifMethod eq 'TNPQ')?("PVAL$adjStrg",'RATIO','CONF_LOW','CONF_UP'):("PVAL$adjStrg",'RATIO','SD_GEO','CONF_LOW','CONF_UP','NUM_PEP_USED'); # * numRatios + NUM_PEP_USED
	my (%numPeptideSets,%usedPeptideSets,%distinctPeptideSets);
	my @codeProt1=($quantifMethod eq 'TNPQ')?("PVAL$adjStrg",'RATIO','CONF_LOW','CONF_UP'):("PVAL$adjStrg",'RATIO','SD_GEO','CONF_LOW','CONF_UP'); # * numRatios + NUM_PEP_USED
	my $nbLine=0;
	if (-e "$resultDir/prot_results.txt") {
		open(PROT1, "$resultDir/prot_results.txt");
		while(<PROT1>) {
			$nbLine++;
			next if $.==1;
			chomp;
			my ($protID0,$protValid,@valueList)=split(/\t/,$_); # fdr ratio sd.geo confInf confSup numPep
			next if !$protID0; # just to be safe
			next if $infiniteRatioProteins{$protID0}; # should be only 1 ratio PP 09/12/14
			next if ($quantifiedModifID && $protID0 !~/-/); # must be a mod prot
			my ($protID,@modResidues)=($quantifiedModifID)? split(/[-.]/,$protID0) : ($protID0,undef);
			my $v=0;
			foreach my $ratioPos (1..$numQuantifRatios) { # should be only 1 PP 09/12/14
				if ($infiniteRatioProteins{$protID0} && $infiniteRatioProteins{$protID0}{$ratioPos}) {
					$v+=scalar @codeProt1;
				}
				else {
					foreach my $quantCode (@codeProt1) {
						if (defined($valueList[$v]) && $valueList[$v] ne 'NA' && $valueList[$v] ne 'Inf') {
							$sthInsProt->execute($protID,$quantifParamIDs{$quantCode},$valueList[$v],$ratioPos);
					$protQuantified{$protID0}{$ratioPos}=1;
							if ($quantifiedModifID) {
								my $protQuantifID=$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF');
								&promsQuantif::insertModifiedResidues($protID,\@modResidues,\%modResiduesID,$dbh,$sthInsModRes); # only 1 insertion/protein even if called in all loops
								&promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$protQuantifID);
							}
						}
						$v++;
					}
				}
			}
			if ($quantifMethod eq 'TNPQ') {
				#my $nbStates= scalar @{$quantifParameters{'R'}{'name.grp'}};
				## PROT_MEAN * num conditions
				#foreach my $state (1..$nbStates) {
				#	#$sthInsProt->execute($protID,$quantifParamIDs{'PROT_MEAN_COND'},$valueList[$v],$quantifParameters{'R'}{'name.grp'}[$state]) if (defined($valueList[$v]) && $valueList[$v] ne 'NA' && $valueList[$v] ne 'Inf');
				#	$v++;
				#}
				## PROT_MEAN * num conditions
				#foreach my $state (1..$nbStates) {
				#	#$sthInsProt->execute($protID,$quantifParamIDs{'COEF_VAR_COND'},$valueList[$v],$quantifParameters{'R'}{'name.grp'}[$state]) if (defined($valueList[$v]) && $valueList[$v] ne 'NA' && $valueList[$v] ne 'Inf');
				#	$v++;
				#}
				## Nb_Replicate * num conditions
				#foreach my $state (1..$nbStates) {
				#	#$sthInsProt->execute($protID,$quantifParamIDs{'NB_REP_COND'},$valueList[$v],$quantifParameters{'R'}{'name.grp'}[$state]) if (defined($valueList[$v]) && $valueList[$v] ne 'NA' && $valueList[$v] ne 'Inf');
				#	$v++;
				#}
				$v+=3 * scalar @{$quantifParameters{'R'}{'name.grp'}};
			}

			###> Add the number of pep used for the computation of the ratio
			if (defined($valueList[$v]) && $valueList[$v] ne 'NA' && $valueList[$v] ne 'Inf') {
				$usedPeptideSets{$protID0}=$valueList[$v]; # NUM_PEP_USED
				##$sthInsProt->execute($protID,$quantifParamIDs{'NUM_PEP_USED'},$valueList[$v],1);
				##if ($quantifiedModifID) {
				##	my $protQuantifID=$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF');
				##	&promsQuantif::insertModifiedResidues($protID,\@modResidues,\%modResiduesID,$dbh,$sthInsModRes);
				##	&promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$protQuantifID);
				##}
			#$protQuantified{$protID0}=1;
			}
		}
		close PROT1;
	}
	if ($nbLine < 2 ) {# File not read (do not exist) or with only headers inside
		$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID");
		$dbh->commit;
		mkdir "$promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
		dircopy($runDir,"$promsPath{quantification}/project_$projectID/quanti_$quantifID");
		$dbh->disconnect;
		die "No proteins were quantified!\n";
	}

	####>Parsing test_normal_power.txt (P_VAL_NORM * ratios + T_TEST_POWER if 2 cond);
	if ($quantifMethod ne 'TNPQ') {
		open(PROT2, "$resultDir/test_normal_power.txt");
		while(<PROT2>) {
			next if $.==1;
			chomp;
			my ($protID0,$protValid,@valueList)=split(/\t/,$_); # p.normFDR pwr
			next if !$protID0; # just to be safe
			next if ($quantifiedModifID && $protID0 !~/-/); # must be a mod prot
			my ($protID,@modResidues)=($quantifiedModifID)? split(/[-.]/,$protID0) : ($protID0,undef);
			my $v=0;
			foreach my $ratioPos (1..$numQuantifRatios) {
				if (defined($valueList[$v]) && $valueList[$v] ne 'NA' && (!$infiniteRatioProteins{$protID0} || !$infiniteRatioProteins{$protID0}{$ratioPos})) {
					$sthInsProt->execute($protID,$quantifParamIDs{"NORM_PVAL$adjStrg"},$valueList[$v],$ratioPos);
					if ($quantifiedModifID) {
						#&promsQuantif::insertModifiedResidues($protID,\@modResidues,\%modResiduesID,$dbh,$sthInsModRes);
						&promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF'));
					}
				}
				$v++;
			}
			if ($numQuantifRatios==1 && defined($valueList[$v]) && $valueList[$v] ne 'NA' && (!$infiniteRatioProteins{$protID0} || !$infiniteRatioProteins{$protID0}{1})) { # only if 2 conditions... so far
				$sthInsProt->execute($protID,$quantifParamIDs{"TTEST_POWER"},$valueList[$v],undef);
				if ($quantifiedModifID) {
					#&promsQuantif::insertModifiedResidues($protID,\@modResidues,\%modResiduesID,$dbh,$sthInsModRes);
					&promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF'));
				}
			}
		}
		close PROT2;
	}

	##>Step 1: Parse pep_mean_perm_ratio.txt for normal ratios
	#my @codePep1=('PEP_RATIO','PEP_MEAN_COND1','PEP_MEAN_COND2','COEF_VAR_COND1','COEF_VAR_COND2','COEF_VAR_OUT','P_VAL_OUT');
	open(PEP, "$resultDir/pep_mean_cv_ratio.txt");
	while(<PEP>) {
		next if $.==1;
		chomp;
		my ($protID0,$pepSeq,$varMod,$charge,$pepIdKey,$protName,$protValid,@valueList)=split(/\t/,$_); # ratio	mean.cond1	mean.cond2	PERM.Repl_Light	PERM.Repl_Heavy	PERM_Out	Pep_Out_0.05
next if $infiniteRatioProteins{$protID0}; # Make sure only 2 conds (1 ratio) for 'Ratio'!
next unless $protQuantified{$protID0};    # "" ""
		$numPeptideSets{$protID0}++; # NUM_PEP_TOTAL
		$distinctPeptideSets{$protID0}{"$pepSeq$varMod"}=1 if $valueList[-1]=~/NA|1/; # DIST_PEP_USED
	#next; # Skipping peptide DB recording for now!!!!!!!!!!!!!!!!!!!!!!!
	}
	close PEP;

	##>Step 2: Read table.txt for inf ratios NUM_PEP_USED, DIST_PEP_USED & NUM_PEP_TOTAL
	open(DATA,"$dataDir/table.txt");
	while(<DATA>) {
		next if $.==1;
		chomp;
		my ($protID0,$pepSeq,$varMod,$charge,$pepIdStrg,$protAlias,$valid,$refXIC,$testXIC)=split(/\t/,$_); # 2 states ONLY!!!
		next unless $infiniteRatioProteins{$protID0};
		$varMod='' unless $varMod;
		$numPeptideSets{$protID0}++; # NUM_PEP_TOTAL (counted ONCE across all states)
		$usedPeptideSets{$protID0}++; # NUM_PEP_USED (=NUM_PEP_TOTAL for inf ratios)
		$distinctPeptideSets{$protID0}{"$pepSeq$varMod"}=1; # DIST_PEP_USED
	}
	close DATA;

	##>Step 3: Save pep pair counts in DB
	foreach my $protID0 (keys %numPeptideSets) {
		my ($protID,@modResidues)=($quantifiedModifID)? split(/[-.]/,$protID0) : ($protID0,undef);
		# NUM_PEP_TOTAL
		$sthInsProt->execute($protID,$quantifParamIDs{'NUM_PEP_TOTAL'},$numPeptideSets{$protID0},undef); # no ratio pos
		if ($quantifiedModifID) {
			&promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF'));
		}
		# NUM_PEP_USED
		$sthInsProt->execute($protID,$quantifParamIDs{'NUM_PEP_USED'},$numPeptideSets{$protID0},undef); # same as NUM_PEP_TOTAL (for inf ratios)
		if ($quantifiedModifID) {
			&promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF'));
		}
		# DIST_PEP_USED
		if ($distinctPeptideSets{$protID0} && $quantifMethod ne 'TNPQ') {
			$sthInsProt->execute($protID,$quantifParamIDs{'DIST_PEP_USED'},scalar keys %{$distinctPeptideSets{$protID0}},undef);
			if ($quantifiedModifID) {
				&promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF'));
			}
		}
	}
}

#close DEBUG;

$sthInsProt->finish;
if ($quantifiedModifID) { # modif quantification
	$sthInsModRes->finish;
	$sthInsProtRes->finish;
}

my $extraQuantifAnnot=($labeling eq 'FREE')? '' : '::QSET_USED='.(scalar keys %{$qSetUsage{'USED'}}).'/'.((scalar keys %{$qSetUsage{'USED'}}) + (scalar keys %{$qSetUsage{'EXCLUDED'}}));

$dbh->do("UPDATE QUANTIFICATION SET QUANTIF_ANNOT=CONCAT(QUANTIF_ANNOT,'$extraQuantifAnnot\:\:LOST_PROTEINS=$numLostProteins') WHERE ID_QUANTIFICATION=$quantifID");

&checkForErrors;


###> Quantification is finished.
$dbh->do("UPDATE QUANTIFICATION SET STATUS=1 WHERE ID_QUANTIFICATION=$quantifID");
$dbh->commit;

$dbh->disconnect;
#exit; # DEBUG!!!
open(FILESTAT,">>$fileStat");
print FILESTAT "Ended ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
close(FILESTAT);

####> 7th: Move all the created files in a specific directory so as to clean the directory
if ($ratioType=~/S\w+Ratio/) {
	#foreach my $scriptName ('AnalysisDiffLimma.R','FunctionLimma.R','AffectDefaultValuesLimma.R','AffectParametersLimma.R') {
	#	unlink "$runDir/$scriptName";
	#}
	unlink "$runDir/AnalysisDiffLimma.R";
}

mkdir "$promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
dirmove($runDir,"$promsPath{quantification}/project_$projectID/quanti_$quantifID");


sleep 2;
unlink $fileStat;

#remove_tree($quantifDir);
#rmtree($quantifDir);


sub checkInfiniteRatios { # globals: $labeling,$ratioType, @ratioTarget, %infiniteRatioProteins, %noSuperRatioProteins
# RULES:
# +Case 1:
#		-FREE: Compute ratio if at least 1 peptide is shared
# +Case 2 (Other labelings + FREE if case 1 failed)
#		-Compute ratio if at least 3 peptides are shared
# +Case 3: (cases 1 & 2 failed)
#		-Compute total XIC for unique to test cond, unique to ref cond and shared by ref & test conds (sum of ref+test)
#		-The highest XIC assigns the ratio type: XICref:-inf, XICtest:+inf, XICshared:normal ratio
	my ($protID,$refNumPepInCond)=@_; # can be a modProtID
	my %noProteinStates; # needed for SuperRatio +/-inf
	my $ratioPos=0;
	foreach my $ratio (@ratioTarget) {
		$ratioPos++;
		$ratio=~s/#//g;
		my ($testCondID,$refCondID)=split(/\//,$ratio);
		if ($ratioType eq 'SuperRatio' && $ratio =~ /%/) { # ratios of ratios if super SILAC (tested AFTER primary ratios in ratios loop!)
			last unless scalar keys %noProteinStates; # inf primary ratios exist
			#$testCondID=~s/%\d+//;
			#$refCondID=~s/%\d+//;
			($testCondID,my $testRatioPos)=($testCondID=~/^(\d+)%(\d+)/);
			($refCondID,my $refRatioPos)=($refCondID=~/^(\d+)%(\d+)/);
			if ($noProteinStates{$testRatioPos}{$testCondID} && $noProteinStates{$refRatioPos}{$refCondID}) {$noSuperRatioProteins{$protID}{$ratioPos}=1;} # special case -> do not record any quantif data
			elsif ($noProteinStates{$testRatioPos}{$testCondID}==1 || $noProteinStates{$refRatioPos}{$refCondID}==-1) {$infiniteRatioProteins{$protID}{$ratioPos}=1000;}
			elsif ($noProteinStates{$testRatioPos}{$testCondID}==-1 || $noProteinStates{$refRatioPos}{$refCondID}==1) {$infiniteRatioProteins{$protID}{$ratioPos}=0.001;}
		}
		##>Primary ratios
		else {
			$noProteinStates{$ratioPos}{$testCondID}=0; # default
		#my $expCode=($ratioType eq 'SuperRatio')? chr(64+$ratioPos) : 'NA';
my $context=($Rdesign eq 'PEP_RATIO')? $ratio : 'NA'; # Experiment != context (ratio can use multiple experiments)
			my ($numShared,$xicShared,$xicMinusInf,$xicPlusInf)=(0,0,0,0);
			my (%usedPepCode,%sharedDistinct);
			foreach my $pepCode (keys %{$refNumPepInCond->{$context}{$testCondID}}) {
				if ($refNumPepInCond->{$context}{$refCondID} && $refNumPepInCond->{$context}{$refCondID}{$pepCode}) { # shared
					#$xicShared+=($refNumPepInCond->{$context}{$testCondID}{$pepCode} >= $refNumPepInCond->{$context}{$refCondID}{$pepCode})? $refNumPepInCond->{$context}{$testCondID}{$pepCode} : $refNumPepInCond->{$context}{$refCondID}{$pepCode};
					$xicShared+=($refNumPepInCond->{$context}{$testCondID}{$pepCode} + $refNumPepInCond->{$context}{$refCondID}{$pepCode}); # Sum of both channels!! (PP 27/10/15)
					$numShared++;
					my ($pepSeq,$vMod,$charge,$dataSrc)=split(':',$pepCode);
					$sharedDistinct{"$pepSeq:$vMod"}=1;
					$usedPepCode{$pepCode}=1;
				}
				else { # only in testCond
				$xicPlusInf+=$refNumPepInCond->{$context}{$testCondID}{$pepCode};
				}
			}
		foreach my $pepCode (keys %{$refNumPepInCond->{$context}{$refCondID}}) {
				next if $usedPepCode{$pepCode};
				$xicMinusInf+=$refNumPepInCond->{$context}{$refCondID}{$pepCode};
			}
			##<Assigning ratio type
			next if ($labeling eq 'FREE' && $numShared); # compute state mean with any available peptide
			next if scalar keys %sharedDistinct >= 3; # normal ratio <===================================================== OVERWRITE THRESHOLD
			next if ($xicShared >= $xicMinusInf && $xicShared >= $xicPlusInf); # normal ratio
			if ($xicPlusInf >= $xicMinusInf) { # +inf
				$infiniteRatioProteins{$protID}{$ratioPos}=1000;
				$noProteinStates{$ratioPos}{$refCondID}=1;
			}
			else { # -inf
				$infiniteRatioProteins{$protID}{$ratioPos}=0.001;
				$noProteinStates{$ratioPos}{$testCondID}=-1;
			}
		}
	}
}

#sub checkInfiniteRatios_bad { # globals: $labeling,$ratioType, @ratioTarget, %infiniteRatioProteins, %noSuperRatioProteins
## RULES:
## +Case 1:
##		-FREE: Compute ratio if at least 1 peptide is shared
## +Case 2 (Other labelings + FREE if case 1 failed)
##		-Compute ratio if at least 3 peptides are shared
## +Case 3: (cases 1 & 2 failed)
##		-Compute total XIC for unique to test cond, unique to ref cond and shared by ref & test conds (sum of ref+test)
##		-The highest XIC assigns the ratio type: XICref:-inf, XICtest:+inf, XICshared:normal ratio
#	my ($protID,$refNumPepInCond)=@_; # can be a modProtID
#	#<Build missing primary ratios for SuperRatio
#	my (@primaryRatios,%primInfRatioProt,%noProteinStates); # needed for SuperRatio +/-inf
#	if ($ratioType eq 'SuperRatio') {
#		my $primRatioPos=0;
#		foreach my $i (1..$#conditionList) {
#			push @primaryRatios,$conditionList[$i].'/'.$conditionList[0];
#		}
#	}
#
#	my ($refRatioList,$refInfRatio)=($ratioType eq 'SuperRatio')? (\@primaryRatios,\%primInfRatioProt) : (\@ratioTarget,\%infiniteRatioProteins);
#	my $ratioPos=0;
#	##>Normal/Primary ratios
#	foreach my $ratio (@{$refRatioList}) {
#		$ratioPos++;
#		$ratio=~s/#//g;
#		my ($testCondID,$refCondID)=split(/\//,$ratio);
#		$noProteinStates{$ratioPos}{$testCondID}=0; # default
#		#my $expCode=($ratioType eq 'SuperRatio')? chr(64+$ratioPos) : 'NA';
#my $context=($ratioType eq 'SuperRatio')? $ratio : 'NA'; # Experiment != context (ratio can use multiple experiments)
#		my ($numShared,$xicShared,$xicMinusInf,$xicPlusInf)=(0,0,0,0);
#		my (%usedPepCode,%sharedDistinct);
#		foreach my $pepCode (keys %{$refNumPepInCond->{$context}{$testCondID}}) {
#			if ($refNumPepInCond->{$context}{$refCondID} && $refNumPepInCond->{$context}{$refCondID}{$pepCode}) { # shared
#				#$xicShared+=($refNumPepInCond->{$context}{$testCondID}{$pepCode} >= $refNumPepInCond->{$context}{$refCondID}{$pepCode})? $refNumPepInCond->{$context}{$testCondID}{$pepCode} : $refNumPepInCond->{$context}{$refCondID}{$pepCode};
#				$xicShared+=($refNumPepInCond->{$context}{$testCondID}{$pepCode} + $refNumPepInCond->{$context}{$refCondID}{$pepCode}); # Sum of both channels!! (PP 27/10/15)
#				$numShared++;
#				my ($pepSeq,$vMod,$charge,$dataSrc)=split(':',$pepCode);
#				$sharedDistinct{"$pepSeq:$vMod"}=1;
#				$usedPepCode{$pepCode}=1;
#			}
#			else { # only in testCond
#				$xicPlusInf+=$refNumPepInCond->{$context}{$testCondID}{$pepCode};
#			}
#		}
#		foreach my $pepCode (keys %{$refNumPepInCond->{$context}{$refCondID}}) {
#			next if $usedPepCode{$pepCode};
#			$xicMinusInf+=$refNumPepInCond->{$context}{$refCondID}{$pepCode};
#		}
#		##<Assigning ratio type
#		next if ($labeling eq 'FREE' && $numShared); # compute state mean with any available peptide
#		next if scalar keys %sharedDistinct >= 3; # normal ratio <===================================================== OVERWRITE THRESHOLD
#		next if ($xicShared >= $xicMinusInf && $xicShared >= $xicPlusInf); # normal ratio
#		if ($xicPlusInf >= $xicMinusInf) { # +inf
#			$refInfRatio->{$protID}{$ratioPos}=1000;
#			$noProteinStates{$ratioPos}{$refCondID}=1;
#		}
#		else { # -inf
#			$refInfRatio->{$protID}{$ratioPos}=0.001;
#			$noProteinStates{$ratioPos}{$testCondID}=-1;
#		}
#	}
#	if ($ratioType eq 'SuperRatio' && scalar keys %noProteinStates) { # inf primary ratios exist => extend to 2ndary ratios
#		$ratioPos=0; # real ratioPos
#		foreach my $ratio (@ratioTarget) {
#			$ratioPos++;
#			$ratio=~s/#//g;
#			my ($testCondID,$refCondID)=split(/\//,$ratio);
#			($testCondID,my $testRatioPos)=($testCondID=~/(\d+)%(\d+)/);
#			($refCondID,my $refRatioPos)=($refCondID=~/(\d+)%(\d+)/);
#			if ($noProteinStates{$testRatioPos}{$testCondID} && $noProteinStates{$refRatioPos}{$refCondID}) {$noSuperRatioProteins{$protID}{$ratioPos}=1;} # special case -> do not record any quantif data
#			elsif ($noProteinStates{$testRatioPos}{$testCondID}==1 || $noProteinStates{$refRatioPos}{$refCondID}==-1) {$infiniteRatioProteins{$protID}{$ratioPos}=1000;}
#			elsif ($noProteinStates{$testRatioPos}{$testCondID}==-1 || $noProteinStates{$refRatioPos}{$refCondID}==1) {$infiniteRatioProteins{$protID}{$ratioPos}=0.001;}
#		}
#	}
#}
#
#sub checkInfiniteRatios_0 { # globals: $labeling,$ratioType, @ratioTarget, %infiniteRatioProteins, %noSuperRatioProteins
## RULES:
## +Case 1:
##		-FREE: Compute ratio if at least 1 peptide is shared
## +Case 2 (Other labelings + FREE if case 1 failed)
##		-Compute ratio if at least 3 peptides are shared
## +Case 3: (cases 1 & 2 failed)
##		-Compute total XIC for unique to test cond, unique to ref cond and shared by ref & test conds (sum of ref+test)
##		-The highest XIC assigns the ratio type: XICref:-inf, XICtest:+inf, XICshared:normal ratio
#	my ($protID,$refNumPepInCond)=@_; # can be a modProtID
#	my %noProteinStates; # needed for SuperRatio +/-inf
#	my $ratioPos=0;
#	foreach my $ratio (@ratioTarget) {
#		$ratioPos++;
#		$ratio=~s/#//g;
#		my ($testCondID,$refCondID)=split(/\//,$ratio);
#		if ($ratioType eq 'SuperRatio' && $ratio =~ /%/) { # ratios of ratios if super SILAC (tested AFTER primary ratios in ratios loop!)
#			last unless scalar keys %noProteinStates; # inf primary ratios exist
#			#$testCondID=~s/%\d+//;
#			#$refCondID=~s/%\d+//;
#			($testCondID,my $testRatioPos)=($testCondID=~/^(\d+)%(\d+)/);
#			($refCondID,my $refRatioPos)=($refCondID=~/^(\d+)%(\d+)/);
#			if ($noProteinStates{$testRatioPos}{$testCondID} && $noProteinStates{$refRatioPos}{$refCondID}) {$noSuperRatioProteins{$protID}{$ratioPos}=1;} # special case -> do not record any quantif data
#			elsif ($noProteinStates{$testRatioPos}{$testCondID}==1 || $noProteinStates{$refRatioPos}{$refCondID}==-1) {$infiniteRatioProteins{$protID}{$ratioPos}=1000;}
#			elsif ($noProteinStates{$testRatioPos}{$testCondID}==-1 || $noProteinStates{$refRatioPos}{$refCondID}==1) {$infiniteRatioProteins{$protID}{$ratioPos}=0.001;}
#		}
#		##>Primary ratios
#		else {
#			$noProteinStates{$ratioPos}{$testCondID}=0; # default
#			my $expCode=($ratioType eq 'SuperRatio')? chr(64+$ratioPos) : 'NA';
#			my ($numShared,$xicShared,$xicMinusInf,$xicPlusInf)=(0,0,0,0);
#			my (%usedPepCode,%sharedDistinct);
#			foreach my $pepCode (keys %{$refNumPepInCond->{$expCode}{$testCondID}}) {
#				if ($refNumPepInCond->{$expCode}{$refCondID} && $refNumPepInCond->{$expCode}{$refCondID}{$pepCode}) { # shared
#					#$xicShared+=($refNumPepInCond->{$expCode}{$testCondID}{$pepCode} >= $refNumPepInCond->{$expCode}{$refCondID}{$pepCode})? $refNumPepInCond->{$expCode}{$testCondID}{$pepCode} : $refNumPepInCond->{$expCode}{$refCondID}{$pepCode};
#					$xicShared+=($refNumPepInCond->{$expCode}{$testCondID}{$pepCode} + $refNumPepInCond->{$expCode}{$refCondID}{$pepCode}); # Sum of both channels!! (PP 27/10/15)
#					$numShared++;
#					my ($pepSeq,$vMod,$charge,$dataSrc)=split(':',$pepCode);
#					$sharedDistinct{"$pepSeq:$vMod"}=1;
#					$usedPepCode{$pepCode}=1;
#				}
#				else { # only in testCond
#					$xicPlusInf+=$refNumPepInCond->{$expCode}{$testCondID}{$pepCode};
#				}
#			}
#			foreach my $pepCode (keys %{$refNumPepInCond->{$expCode}{$refCondID}}) {
#				next if $usedPepCode{$pepCode};
#				$xicMinusInf+=$refNumPepInCond->{$expCode}{$refCondID}{$pepCode};
#			}
#			##<Assigning ratio type
#			next if ($labeling eq 'FREE' && $numShared); # compute state mean with any available peptide
#			next if scalar keys %sharedDistinct >= 3; # normal ratio <===================================================== OVERWRITE THRESHOLD
#			next if ($xicShared >= $xicMinusInf && $xicShared >= $xicPlusInf); # normal ratio
#			if ($xicPlusInf >= $xicMinusInf) { # +inf
#				$infiniteRatioProteins{$protID}{$ratioPos}=1000;
#				$noProteinStates{$ratioPos}{$refCondID}=1;
#			}
#			else { # -inf
#				$infiniteRatioProteins{$protID}{$ratioPos}=0.001;
#				$noProteinStates{$ratioPos}{$testCondID}=-1;
#			}
#		}
#	}
#}
#
#
#sub insertModifiedResidues {
#	my ($protID,$refModResList,$refModResiduesID,$dbh,$sthInsModRes)=@_;
#	my @newList;
#	if ($refModResList->[0]=~/~/) {&encodeModifPosAmbiguity($refModResList,\@newList)} else {@newList=@{$refModResList};} # $refModResList must not be overwritten!
##print DEBUG "\tMODRES_ID='@newList'\n"; # unless $refModResiduesID->{$modRes};
#	foreach my $modRes (@newList) {
#		if (!$refModResiduesID->{$protID} || !$refModResiduesID->{$protID}{$modRes}) {
#			my ($res,$pos)=($modRes=~/^(.)(.+)/); # Y180 -> (Y,180)
#			$sthInsModRes->execute($res,$pos);
#			$refModResiduesID->{$protID}{$modRes}=$dbh->last_insert_id(undef,undef,'MODIFIED_RESIDUE','ID_MODIF_RES');
#		}
#	}
#}
#sub insertProteinModResidues {
#	my ($refModResList,$sthInsProtRes,$refModResiduesID,$protQuantifID)=@_;
#	my @newList;
#	if ($refModResList->[0]=~/~/) {&encodeModifPosAmbiguity($refModResList,\@newList)} else {@newList=@{$refModResList};} # $refModResList must not be overwritten!
##print DEBUG "\tLIST='@newList'\n"; # unless $refModResiduesID->{$modRes};
#	foreach my $modRes (@newList) {
#		$sthInsProtRes->execute($refModResiduesID->{$modRes},$protQuantifID);
#	}
#}
#sub encodeModifPosAmbiguity { # (PP 2017/06/22)
#	# ex1 (res~res): 123~135:2/3@<rt> --> (-123,+135,@<rt>,23) WARNING: in x/y x,y must be between [1..9]
#	# ex2 (Nter~res): n0~7:2/3 --> (=0,+7,23)  "-" becomes "=" to indicate N-term (n0: 0 indicates Protein N-term)
#	# ex3 (res~Cter): 251~c0:2/3 --> (-251,*0,23)  "+" becomes "*" to indicate C-term (c0: 0 indicates Protein C-term)
#	# ex4 (Nter~Cter): n78~c85:2/3 --> (=78,*85,23)
#	my ($refModResList,$refCodedModResList)=@_;
#	my @modRes=split('[~:@/]',$refModResList->[0]);
#	#@{$refCodedModResList}=('-'.$modRes[0],'+'.$modRes[1]);
#	#push @{$refCodedModResList},'@'.$modRes[2] if $modRes[4]; # has RT data (optional)
#	#push @{$refCodedModResList},$modRes[-2].$modRes[-1];
#	# Starting pos
#	if ($modRes[0]=~/^n(\d+)/) {@{$refCodedModResList}=('='.$1);} # N-ter
#	else {@{$refCodedModResList}=('-'.$modRes[0]);} # res
#	# Ending pos
#	if ($modRes[1]=~/^c(\d+)/) {push @{$refCodedModResList},'*'.$1;} # C-ter
#	else {push @{$refCodedModResList},'+'.$modRes[1];} # res
#	# RT
#	push @{$refCodedModResList},'@'.$modRes[2] if $modRes[4]; # has RT data (optional)
#	# Matched & Available pos
#	push @{$refCodedModResList},$modRes[-2].$modRes[-1];
#}

sub checkForErrors {
	#my ($errorFile) = glob "$promsPath{tmp}/quantification/current/*_$quantifDate\_error.txt";
	my $errorFile = "$promsPath{tmp}/quantification/current/$quantifID\_$quantifDate\_error.txt";
	if ($errorFile && -s $errorFile) {
#system "sed -i '/Permanently/ d' $errorFile"; # remove warning from ssh key addition
#		if (-s $errorFile) {
			$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID"); # failed
			$dbh->commit;
			$dbh->disconnect;
			die "Aborting quantification due to errors.";
		#}
	}
}

# TODO: Make clear choice for labeled quantif done with PEP_INTENSITY algo: treat as 100% Label-Free or mixed LF/Label ?????
# TODO: Move label-free peptide matching check further downstream for compatibility with PTM quantif
####>Revision history<####
# 2.7.0 Now runs on cluster itself: R is launcher by system command (PP 17/09/18)
# 2.6.2 [Fix] bug in splitting peptide id string during results parsing (PP 13/09/18)
# 2.6.1 In TopX mode, prevents same peptide sequence to reused unless topX is not reached (PP 14/06/18)
# 2.6.0 Reads peptide XIC from file(s) not DB (PP 09/05/18)
# 2.5.7 [Fix] Bug no protein list-based normalization because wrong parameter was read (PP 20/04/18)
# 2.5.6 PBSerror is now handled by &clusterInfo{checkError}. Requires promsConfig.pm v2.9.0 or higher (PP 11/04/18)
# 2.5.5 [Fix] Bug lost of infinite ratios for PEP_RATIO (PP 06/04/18)
# 2.5.4 [Fix] Bug empty lostProteins.txt file (PP 13/02/18)
# 2.5.3 Quantif status updated to -2 for all die cases (PP 27/12/17)
# 2.5.2 Fix bug in LABEL-FREE with modif allowed due to usage of $qSetBestVarMod (PP 14/12/17)
# 2.5.1 Added list of non-quantified proteins & displays R error if any (PP 01/12/17)
# 2.5.0 Adapted to new quantif algo v3 & relies on &promsConfig::getClusterInfo (PP 27/11/17)
# 2.4.1 Bug fix for peptides being totally excluded when shared by 2+ proteins visible across different analyses even when non-unique option is selected (PP 14/11/17)
# 2.4.0 Compatible with Prot/Pep N/C-term quantification and position ambiguity (PP 02/08/17)
# 2.3.9 Updates quantification status to "failed" (-2) in case of error & minor bug fix (PP 21/03/17)
# 2.3.8b Minor update for TMT (GA 04/04/17)
# 2.3.8 Minor update in peptide main query to handle peptide multi-matches within same protein (PP 17/02/17)
# 2.3.7 Compatible with TMT labeling (PP 18/01/17)
# 2.3.6 Remove XIC duplication created by MassChroQ from PD data because some PD qSets can be associated with a fractionated single XIC peak (PP 31/12/16)
# 2.3.5 Bug fix for replicates in MultiSILAC & bug fix for PD/MCQ quantif selection & improved best charge selection (PP 06/10/16)
# 2.3.4 Bug fix in QSet merge for MultiSILAC & full exclusion of peptides shared by multiple top proteins (PP 12/05/16)
# 2.3.3 Bug fix modif quantif in label-free & with old algo (PP 28/04/16)
# 2.3.2 Peptide matching option for label-free topN (PP 25/02/16)
# 2.3.1 Bug fix for internal quantif due to SOFTWARE=xxx tag in peptide quantif annot (PP 18/02/16)
# 2.3.0 Modif-quantification,QSET merge, improve +/-inf decision, "Simple ratios"=multiSILAC algorithm for SILAC triplex (PP 29/01/16)
# 2.2.2 QSET unicity also when no pooling (PP 23/06/14)
# 2.2.1 Bug fix in 'experiment' declaration for SuperRatio when true replicates (PP 20/06/14)
# 2.2.0 Handles new quantification strategy ratio/ratio & non-pep ratio label-free (PP 19/06/14)
# 2.1.1 Bug fix of detection of SILAC unlabeled channel: impact of +/-inf ratios (PP 09/05/14)
# 2.1.0 Allows given PTMs to go through quantification when PTM not allowed (PP 11/04/14)
# 2.0.9 Fration intensities are no longer summed if labeled quantif (PP 08/04/14)
# 2.0.8 Assumes no-label channel if incomplete annot data (PP 03/04/14)
# 2.0.7 Bug fix for label-free when DATA_SOURCE is defined (PP 28/03/14)
# 2.0.6 Bug fix for forgotten mkdir (PP 19/03/14)
# 2.0.5 Uses rmtree instead of remove_tree (PP 10/03/14)
# 2.0.4 Change state/replicate label in R files to stateX(.repN) (PP 26/02/14)
# 2.0.3 Change PEPTIDE.DATA field parsing for compatibility with PD (QSET=) & MassChroQ (MCQSET_pepQuantID=) from isotope extraction (PP 07/02/14)
# 2.0.2 system command removal (PP 13/11/13)
# 2.0.1 Added infinite ratios switch (PP 21/10/13)
# 2.0.0 Merging runLabeledQuantification.pl with runDesignQuantification.pl (PP 09/10/13)
# 1.3.4 Corrections in 1.3.3 update & error check (PP 30/09/13)
# 1.3.3 Include partially hidden proteins (PP 27/09/13)
# 1.3.2 Remove forgotten debug.txt file deletion (PP 26/09/13)
# 1.3.1 Checking NA values before lauching R script (FY 09/09/13)
# 1.3.0 Generated from runTNPQQuantification.pl (PP 03/09/13)
# 1.2.1 Infinite ratio detection (PP 30/08/13)
# 1.2.0 Generic design-based quantification with fractions sum up (PP 11/07/13)
# 1.1.6 Minor changes & no error.log on directory manipulation (PP 08/07/13)
# 1.1.5 Move quantif run dir even when error occurs (PP 03/07/13)
# 1.1.4 Remove VAR_MOD from script (GA 28/05/13)
# 1.1.3 Use new R Script LabeledAndUnlabeledQuanti.Function.R -> require to add 'typeTest' in R parameters (GA 15/05/13)
# 1.1.2 Minor modif to inactivate FDR control where numValidProt=1 (GA 28/03/13)
# 1.1.1 Update validProt options to quantify top-match proteins only (GA 07/03/13)
# 1.1.0 Added more peptide filters: missed cleavage and PTM/unmodified peptide pair exclusion (PP 11/01/13)
# 1.0.9 Add bias correction with reference proteins (PP 09/01/13)
# 1.0.8 Add new columns for TNPQ in result-parsing (GA 08/11/12)
# 1.0.7 Modification to handle analyses extracted at the same time by MassChroQ + QUANTIF_ANNOT modification + multi-ratios handle (GA 07/11/12)
# 1.0.6 Minor modif for DIST_PEP_USED -> value not available for TnPQ  (GA 18/10/12)
# 1.0.5 Added distinct peptides used (DIST_PEP_USED) calculation (PP 10/05/12)
# 1.0.4 Add a column Peptide_IDs in the file sent to R for computation of ratio/tnpq (GA 11/04/12)
# 1.0.3 Clean version of the script (GA 10/04/12)
# 1.0.2 Update to PP parsing (GA 22/02/2012)
# 1.0.1 Update to new scripts of FC (a lot of new parameters have been added...) (GA 11/08/2011)
# 1.0.0 New script to handle TNPQ scripts of FC
