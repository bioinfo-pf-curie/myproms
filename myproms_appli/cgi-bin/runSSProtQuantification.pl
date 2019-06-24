#!/usr/local/bin/perl -w

################################################################################
# runSSProtQuantification.pl       1.2.3                                       #
# Component of site myProMS Web Server                                         #
# Authors: P. Poullet (Institut Curie)                                         #
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
use File::Copy::Recursive qw(dirmove);
use File::Path;
#exit;
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo('no_user');

###############################
####>Recovering parameters<####
###############################
my ($quantifID,$quantifDate)=@ARGV;
my $quantifDir="$promsPath{tmp}/quantification/$quantifDate";
unless (-e $quantifDir) {
	die "ERROR: Job directory '$quantifDate' not found!";
}
my %quantifParameters=&promsQuantif::getQuantificationParameters("$quantifDir/quantif_info.txt");

my $dbh=&promsConfig::dbConnect('no_user');
$dbh->do("UPDATE QUANTIFICATION SET STATUS=0 WHERE ID_QUANTIFICATION=$quantifID"); # running
$dbh->commit;

my $runDir="$quantifDir/quanti_$quantifID";
my $dataDir="$runDir/data";
my $resultDir="$runDir/results";
my $graphDir="$resultDir/graph";
#make_path($dataDir,$graphDir,{verbose=>0,mode=>0775}); # $runDir,$resultDir will be created automatically
mkdir $runDir;
mkdir $dataDir;
mkdir $resultDir;
mkdir $graphDir;

#open (DEBUG,">$promsPath{tmp}/quantification/debug.txt"); # DEBUG
#print DEBUG "PEPTIDES: '@{$quantifParameters{'DB'}{PEPTIDES}}'\n";
my $fileStat="$quantifDir/status_$quantifID.out";
open(FILESTAT,">$fileStat") || die "Error while opening $fileStat";
print FILESTAT "Started ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
print FILESTAT "1/4 Fetching data\n";
close FILESTAT;

my $projectID=&promsMod::getProjectID($dbh,$quantifID,'QUANTIFICATION');

my ($quantifAnnot,$quantifiedModifID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_MODIFICATION FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID");
$quantifAnnot=~s/::/\$/g;
my ($labeling)=($quantifAnnot=~/LABEL=([^\$]+)/);
$labeling=uc($labeling);
my ($stateStrg)=($quantifAnnot=~/STATES=([^\$]+)/);
$stateStrg=~s/#//g; # remove all id tags


################################
####>Get states composition<####
################################
my (%quantifDesign,%obs2Ana); # SSPA
my (%ana2Obs,%labelModifs,%anaLabelModifs);

my $sthALM=$dbh->prepare("SELECT M.ID_MODIFICATION FROM ANALYSIS_MODIFICATION AM,MODIFICATION M WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND AM.ID_ANALYSIS=? AND M.IS_LABEL=1");
my $sthObsMod=$dbh->prepare("SELECT M.ID_MODIFICATION FROM OBSERVATION O LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION WHERE O.ID_OBSERVATION=?"); # returns NULL if no modif!!!
#my $sthObsBioSamp=$dbh->prepare("SELECT ID_BIOSAMPLE FROM OBSERVATION WHERE ID_OBSERVATION=?");
#my $sthObsMsSamp=$dbh->prepare("SELECT ID_SAMPLE FROM ANALYSIS A WHERE ID_ANALYSIS=?");

my $numStates=0;
foreach my $state (split(/;/,$stateStrg)) { # 'state1;state2;---;stateN',  stateX='nbBioRep,bioRep1.bioRep2.---.bioRepN,#condID', bioRepX='techRep1=techRep2=---=techRepN', techRepX='frac1+frac2+---+fracN'
	my ($numBioRep,$quantiObsIDs,$condID)=split(/,/,$state);
#next if $condID==54; # SILAC ref
	$numStates++;
	@{$quantifDesign{$numStates}}=();
	#my $bioRepCount=0;
	foreach my $bioReplicate (split(/\./,$quantiObsIDs)) {
		#$bioRepCount++;
#last if $bioRepCount==2;
		my @bioRep;
		#my %nameList;
		my $numBioRepObs=0;
		foreach my $techReplicate (split(/&/,$bioReplicate)) {
			my @fractionObs;
#my $numObs=0;
			foreach my $fraction (split(/\+/,$techReplicate)) {
#last if ++$numObs==3;
				$numBioRepObs++;
				my ($obsID,$parQuantifID,$anaID,$targetPos)=split(/:/,$fraction); # $parQuantifID can be 0 for SSPA
				$obs2Ana{$obsID}=$anaID;
				push @fractionObs,$obsID;
				#$sthObsBioSamp->execute($obsID);
				#my ($bsID)=$sthObsBioSamp->fetchrow_array;
				#$nameList{BIOSAMPLE}{$bsID}++ if $bsID;
				#$sthObsMsSamp->execute($anaID);
				#my ($sampID)=$sthObsMsSamp->fetchrow_array;
				#$nameList{MSSAMPLE}{$sampID}++;
				%{$labelModifs{$obsID}}=();
				push @{$ana2Obs{$anaID}},$obsID;
				if ($labeling eq 'SILAC') { # SILAC QUANTI
					unless ($anaLabelModifs{$anaID}) { # get list of labeling modif for Analysis
						$sthALM->execute($anaID);
						while (my ($modID)=$sthALM->fetchrow_array) {$anaLabelModifs{$anaID}{$modID}=1;}
					}
					$sthObsMod->execute($obsID);
					while (my ($modID)=$sthObsMod->fetchrow_array) {
						$modID=0 unless $modID; # undef for no-label channel
						$labelModifs{$obsID}{$modID}=1; # {0} for no-label obs
						#$anaLabelModifs{$anaID}{$modID}=1 if $modID;
					}
				}
				# iTRAQ case is handled when fetching peptide quanti data (<- ??? PP 20/03/15)
				elsif ($labeling =~ /ITRAQ|TMT/) {
					my $modID;
					if ($anaLabelModifs{$anaID}) {$modID=(keys %{$anaLabelModifs{$anaID}})[0];}
					else {
						$sthALM->execute($anaID);
						($modID)=$sthALM->fetchrow_array;
next unless $modID; # !!!TEMP!!! PARAGON quantif with no validated protein (PP 30/03/15)
						$anaLabelModifs{$anaID}{$modID}=1;
					}
					$labelModifs{$obsID}{$modID}=1; # if $modID;
				}
				elsif ($labeling eq 'FREE') {$anaLabelModifs{$anaID}{0}=1; $labelModifs{$obsID}{0}=1;}
				else {die "ERROR: Unrecognized labeling ('$labeling')!";}
			} # end of fraction
			push @bioRep,\@fractionObs;

		} # end of tech replicate (sample)

		##>Guess most suitable name for bioRep
		#if ($nameList{BIOSAMPLE}) { # Bio-sample
		#	my $bsID=(sort{$nameList{BIOSAMPLE}{$b} <=> $nameList{BIOSAMPLE}{$a} || $a<=>$b} keys %{$nameList{BIOSAMPLE}})[0];
		#	$bioRep{ID}='b'.$bsID if $nameList{BIOSAMPLE}{$bsID}==$numBioRepObs;
		#}
		#unless ($bioRep{ID}) { # MS Sample
		#	my $sampID=(sort{$nameList{MSSAMPLE}{$b} <=> $nameList{MSSAMPLE}{$a} || $a <=> $b} keys %{$nameList{MSSAMPLE}})[0];
		#	$bioRep{ID}='s'.$sampID if $nameList{MSSAMPLE}{$sampID}==$numBioRepObs;
		#	$bioRep{ID}='c'.$bioRepCount unless $bioRep{ID};
		#}
		#$bioRep{LABEL}="State$statePos.BioRep$bioRepCount";
		push @{$quantifDesign{$numStates}},\@bioRep;

	} # end of bio-replicate

} # end of cond/state

$sthALM->finish;
$sthObsMod->finish;
#$sthObsBioSamp->finish;
#$sthObsMsSamp->finish;

#print Dumper %quantifDesign;
#print "",Dumper %ana2Obs
#$dbh->disconnect;
#exit;


#################################################################
####>Fetching and filtering peptides associated with samples<####
#################################################################
# Retrieve all visible & hidden proteins but keeps only those visible in at least 1 analysis
my ($pepSpecificity,$pepMissedCleavage,$ptmFilter,$pepFocus)=@{$quantifParameters{'DB'}{'PEPTIDES'}};
#my ($pepSpecificity,$pepMissedCleavage,$ptmFilter,$pepFocus)=('unique_shared',1,1,'psm');
my $chargeStrg=($pepFocus=~/sp_count|_ion/)? 'CHARGE' : '0'; # sp_count|all_ion|all_pep|dist_ion|dist_pep|dist_seq (before: peptide or pepSeq)
my $peptideQuery=qq|SELECT PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),ABS(PEP_BEG),CHARGE,SPEC_COUNT
		FROM PEPTIDE P
		LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
		INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
		WHERE P.ID_ANALYSIS=? AND PPA.ID_PROTEIN=?|; # No need for group_concat on PEP_BEG since only 1 position is usable
	$peptideQuery.=($pepSpecificity eq 'unique')? ' AND PPA.IS_SPECIFIC=1' : ($pepSpecificity eq 'unique_shared')? ' AND AP.PEP_SPECIFICITY=100' : ''; # Filter at DB level
	$peptideQuery.=' AND MISS_CUT=0' unless $pepMissedCleavage;
	$peptideQuery.=' GROUP BY P.ID_PEPTIDE'; # PPA.ID_PROTEIN,
#print DEBUG "\nQUERY=$peptideQuery\n";
my $sthGetPep=$dbh->prepare($peptideQuery);

$ptmFilter=~s/#//g; # remove id flag
my ($ptmAllowed,@selectedPTMs)=split(/[:,]/,$ptmFilter); # @selectedPTMs defined only if $ptmAllowed >= 2
my %allowedPtmID;
foreach my $modID (@selectedPTMs) {$allowedPtmID{$modID}=1;}
$allowedPtmID{$quantifiedModifID}=1 if $quantifiedModifID; # to allow modif quantification

##>Data file
my (@bioRepOrder,%anaData,%proteinData,%bestVisibility,%excludedSeq);
foreach my $statePos (sort{$a<=>$b} keys %quantifDesign) {
#print DEBUG ">STATE #$statePos ($stateLabel)\n";
	my $bioRepPos=0;
	foreach my $refBioRep (@{$quantifDesign{$statePos}}) {
		$bioRepPos++;
		my $bioRepLabel="State$statePos.BioRep$bioRepPos";
#print DEBUG "  -BioRep $bioRepLabel (x",scalar @{$refBioRep->{OBS}}," TechRep):\n";
		push @bioRepOrder,$bioRepLabel;

		####>Data
		my %techRepData;
		my $numTechRep=0;
		foreach my $refTechRep (@{$refBioRep}) {
			$numTechRep++;
#print DEBUG "#@{$refTechRep}#\n";
			my %fracSpCount;
			foreach my $obsID (@{$refTechRep}) { # 1 obs = 1 fraction
				%{$fracSpCount{$obsID}}=(); # defined so no need to checked later
				my $anaID=$obs2Ana{$obsID};
				my (%missingSpCount,%obsIonCount);
				if (!$anaData{$anaID}) {
#print DEBUG "($anaID,$obsID)"; #{@{$ana2Obs{$anaID}}}
					$sthGetPep->execute($anaID,$anaID);
					PEP:while (my ($protStrg,$pepSeq,$varModStrg,$charge,$spCount)=$sthGetPep->fetchrow_array) {
						next if $excludedSeq{$pepSeq}; # skip any future seq
						#$varModStrg='' unless $varModStrg;
						next if ($quantifiedModifID && (!$varModStrg || $varModStrg !~ /(^|&)$quantifiedModifID:/)); # skip peptides not containing quantified modif
						$spCount=0 unless $spCount; # in case only ghost (not defined), set to 1 if feature is not sp_count
						my @protList;
						foreach my $protData (split('&',$protStrg)) {
							my ($protID,$visibility)=split(':',$protData);
							last unless defined $visibility; # truncated string if peptide matches too many proteins (response max length=1024 chars)
							push @protList,$protID;
							$bestVisibility{$protID}=$visibility if (!$bestVisibility{$protID} || $bestVisibility{$protID} < $visibility);
						}
						#>Processing var mods & matching label channel (checking allowed & removing label if labeled quantif)
						my @varMods=($varModStrg)? split(/&/,$varModStrg) : ();
						my $varModCode='';
						##if ($labeling eq 'SILAC') {
						my $matchedChannel=0;
						foreach my $anaObsID (@{$ana2Obs{$anaID}}) {
							my $matchedObs=0;
							if ($labeling eq 'FREE') {
								$matchedObs=1;
							}
							else {
								my $isLabeled=0;
								foreach my $vMod (@varMods) {
									my ($modID)=($vMod=~/^(\d+)/);
									if ($labelModifs{$anaObsID}{$modID}) { # matched labeled channel
										$matchedObs=1;
										$isLabeled=1;
										##last;
									}
									elsif ($anaLabelModifs{$anaID}{$modID}) { # un-matched labeled channel
										$isLabeled=1;
									}
								}
								if (!$isLabeled && $labelModifs{$anaObsID}{0}) { # Non-labeled channel
									$matchedObs=1;
								}
							}
							if ($matchedObs) {
								my $hasExcludedPTM=0;
								foreach my $vMod (@varMods) {
									my ($modID)=($vMod=~/^(\d+)/);
									next if $labelModifs{$anaObsID}{$modID}; # skip labeling modifs
									if ($ptmAllowed != 1 && !$allowedPtmID{$modID}) {
										$hasExcludedPTM=1;
										last;
									}
									$varModCode.='&'.$vMod;
								}
								if ($hasExcludedPTM) { # this modification is not allowed
									if ($ptmAllowed <= -1) { # -1 or -2 unmodified peptide not allowed if exists a modified version (unmodif pep may have been recorded before detection => clean below)
										foreach my $protID (@protList) {$excludedSeq{$pepSeq}{$protID}=1;}
									}
									next PEP; # skip all other obs for this ana
								}
								if ($anaObsID == $obsID) {$matchedChannel=1;}
								else { # stored to be recalled when proper obsID is wanted
									push @{$anaData{$anaID}{$anaObsID}},[\@protList,$pepSeq,$varModCode,$charge,$spCount];
								}
								last; # obs in Ana: current anaObs was matched by peptide: no need to scan others
							}
						}

						###>Recording peptides
						if ($matchedChannel) {
							$varModCode='' if $pepFocus eq 'dist_seq';
							my $ionKey="$pepSeq:$varModCode:$charge";
							if ($pepFocus eq 'sp_count') { # record best spCount/ion in each fraction
								foreach my $protID (@protList) {
									#$techRepData{$protID}{$numTechRep}{"$pepSeq:$varModCode:$charge"}=$spCount if (!$techRepData{$protID} || !$techRepData{$protID}{$numTechRep} || !$techRepData{$protID}{$numTechRep}{"$pepSeq:$varModCode:$charge"} || $techRepData{$protID}{$numTechRep}{"$pepSeq:$varModCode:$charge"} < $spCount);
									$fracSpCount{$obsID}{$protID}{$ionKey}=$spCount if (!$fracSpCount{$obsID}{$protID} || !defined $fracSpCount{$obsID}{$protID}{$ionKey} || $fracSpCount{$obsID}{$protID}{$ionKey} < $spCount);
									$missingSpCount{$protID}{$ionKey}=1 if $spCount==0;
									$obsIonCount{$protID}{$ionKey}++;
								}
							}
							else { # all|distinct ion|pep|seq
								foreach my $protID (@protList) {$techRepData{$protID}{$numTechRep}{$ionKey}++;}
							}
#print DEBUG '+';
						}
					}
				}
				else { # analysis data already fetched & assigned to an Obs =>  $anaData{$anaID}{$obsID} defined
#print DEBUG "[$anaID,$obsID]=",scalar @{$anaData{$anaID}{$obsID}};
					foreach my $refPep (@{$anaData{$anaID}{$obsID}}) {
						my ($refProtList,$pepSeq,$varModCode,$charge,$spCount)=@{$refPep};
						$varModCode='' if $pepFocus eq 'dist_seq';
						my $ionKey="$pepSeq:$varModCode:$charge";
						if ($pepFocus eq 'sp_count') { # record best spCount for each ion in each fraction (cannot be directly summed because of multiple occurences of same ion with same spCount)
							foreach my $protID (@{$refProtList}) {
								#$techRepData{$protID}{$numTechRep}{"$pepSeq:$varModCode:$charge"}=$spCount if (!$techRepData{$protID} || !$techRepData{$protID}{$numTechRep} || !$techRepData{$protID}{$numTechRep}{"$pepSeq:$varModCode:$charge"} || $techRepData{$protID}{$numTechRep}{"$pepSeq:$varModCode:$charge"} < $spCount);
								$fracSpCount{$obsID}{$protID}{$ionKey}=$spCount if (!$fracSpCount{$obsID}{$protID} || !defined $fracSpCount{$obsID}{$protID}{$ionKey} || $fracSpCount{$obsID}{$protID}{$ionKey} < $spCount);
								$missingSpCount{$protID}{$ionKey}=1 if $spCount==0;
								$obsIonCount{$protID}{$ionKey}++;
							}
						}
						else { # all|distinct ion|pep|seq
							foreach my $protID (@{$refProtList}) {$techRepData{$protID}{$numTechRep}{$ionKey}++;}
						}
					}
				}
				##>Replace sp_count with all_ions where not enough sp_count
				if ($pepFocus eq 'sp_count') {
					foreach my $protID (keys %missingSpCount) {
						foreach my $ionKey (keys %{$missingSpCount{$protID}}) {
							$fracSpCount{$obsID}{$protID}{$ionKey}=$obsIonCount{$protID}{$ionKey} if $fracSpCount{$obsID}{$protID}{$ionKey} < $obsIonCount{$protID}{$ionKey};
						}
					}
				}
			} # End of obs/fraction

			###>Sum each ion's spCount across fractions
			if ($pepFocus eq 'sp_count') {
				foreach my $obsID (keys %fracSpCount) {
					foreach my $protID (keys %{$fracSpCount{$obsID}}) {
						foreach my $ionKey (keys %{$fracSpCount{$obsID}{$protID}}) {
							$techRepData{$protID}{$numTechRep}{$ionKey}+=$fracSpCount{$obsID}{$protID}{$ionKey};
						}
					}
				}
			}

			###>Clean excluded pepSeq ($ptmAllowed <= -1)
			foreach my $pepSeq (keys %excludedSeq) {
				foreach my $protID (keys %{$excludedSeq{$pepSeq}}) {
					my @peps=grep{/^$pepSeq:/} keys %{$techRepData{$protID}{$numTechRep}};
					foreach my $pep (@peps) {
						delete $techRepData{$protID}{$numTechRep}{$pep};
						delete $techRepData{$protID}{$numTechRep} unless scalar keys %{$techRepData{$protID}{$numTechRep}};
						delete $techRepData{$protID} unless scalar keys %{$techRepData{$protID}};
					}
				}
				%{$excludedSeq{$pepSeq}}=(); # reset list of matching prot to 0 before next techRep
			}

		} # End of techRep

		###>Counting peptides (aggregate techReps -> mean of pepCount)
		if ($pepFocus=~/^(sp|all)_/) { # sp_count | all_(ion|pep)
			foreach my $protID (keys %techRepData) {
				foreach my $techRepPos (keys %{$techRepData{$protID}}) {
					foreach my $count (values %{$techRepData{$protID}{$techRepPos}}) {
						$proteinData{$protID}{$bioRepLabel}+=$count;
					}
				}
				$proteinData{$protID}{$bioRepLabel}=int(0.5 + ($proteinData{$protID}{$bioRepLabel}/$numTechRep));
			}
		}
		else { # ion or peptide or pepSeq (distinct +/- charge)
			foreach my $protID (keys %techRepData) {
				foreach my $techRepPos (keys %{$techRepData{$protID}}) {
					$proteinData{$protID}{$bioRepLabel}+=scalar keys %{$techRepData{$protID}{$techRepPos}};
				}
				$proteinData{$protID}{$bioRepLabel}=int(0.5 + ($proteinData{$protID}{$bioRepLabel}/$numTechRep));
			}
		}
	} # End of bioRep

} # End of State
$sthGetPep->finish;


####################################
####>Printing data to table.txt<####
####################################
open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
print FILESTAT "2/4 Printing data to file\n";
close FILESTAT;

open(DATA,">$dataDir/table.txt"); # Name of data table
##>Header
print DATA 'Protein';
foreach my $bioRepLabel (@bioRepOrder) {print DATA "\t$bioRepLabel";}
##>Data
if (scalar keys %proteinData) {
	foreach my $protID (sort {$a<=>$b} keys %proteinData) {
		next if $bestVisibility{$protID}==0;
		print DATA "\n".$protID;
		foreach my $bioRepLabel (@bioRepOrder) {
			if ($proteinData{$protID}{$bioRepLabel}) {print DATA "\t".$proteinData{$protID}{$bioRepLabel};}
			else {print DATA "\t0";}
		}
	}
}
else {warn "ERROR: No proteins quantified!\n";}
close DATA;

&checkForErrors;

$dbh->disconnect;
#close DEBUG;
#die "Data file created!"; # DEBUG



								############################################
			#############################> RUNNING R SCRIPTS <##########################
								############################################
open(FILESTAT,">>$fileStat");
print FILESTAT "3/4 Running SSP Analysis\n";
close FILESTAT;

my %cluster=&promsConfig::getClusterInfo; #('debian') # default is 'centos'
my $pathR=($cluster{'on'})? $cluster{'path'}{'R'} : $promsPath{'R'};

open(R_SCRIPT,">$runDir/analysisCounting.R");
print R_SCRIPT qq
|
###################################################################################
# Launcher for quantification R scripts (provides path for sourcing dependencies) #
###################################################################################

filepath="$promsPath{R_scripts}/"
source(paste(filepath,"analysisCounting.R",sep=""))
|;
close R_SCRIPT;

my $RcommandString="export LANG=en_US.UTF-8; cd $runDir; $pathR/R CMD BATCH --no-save --no-restore analysisCounting.R";

if ($cluster{'on'}) { ###>Run job on cluster
#	my $bashFile = "$runDir/runSSProtQuantification.sh";
#	my $clusterRcommandString=$cluster{'buildCommand'}->($runDir,$RcommandString);
	my $numPepValues=(scalar keys %proteinData) * (scalar @bioRepOrder);
	my $maxHours=int(0.5+($numPepValues/25000)); $maxHours=2 if $maxHours < 2; $maxHours=48 if $maxHours > 48; # 1 M lines -> 40 h
	my $maxMem=int(1.5 + 1E-6 * $numPepValues);
	$maxMem.='Gb';
#	open (BASH,">$bashFile");
#	print BASH qq
#|#!/bin/bash
###resources
##PBS -l mem=$maxMem
##PBS -l nodes=1:ppn=1
##PBS -l walltime=$maxHours:00:00
##PBS -q batch
###Information
##PBS -N myProMS_SSPA_$quantifID
###PBS -M patrick.poullet\@curie.fr
##PBS -m abe
##PBS -o $runDir/PBS.txt
##PBS -e $runDir/PBSerror.txt
#
### Command
#$clusterRcommandString
#echo _END_$quantifID
#|;
#	close BASH;
#	my $modBash=0775;
#	chmod $modBash, $bashFile;
#	#system "$promsPath{qsub}/qsub $bashFile > $runDir/torqueID.txt";
#	$cluster{'sendToCluster'}->($bashFile);
#	sleep 30;
#
#	###>Waiting for R job to run
#	my $pbsError;
#	my $nbWhile=0;
#	my $maxNbWhile=$maxHours*60*2;
#	while ((!-e "$runDir/PBS.txt" || !`tail -3 $runDir/PBS.txt | grep _END_$quantifID`) && !$pbsError) {
#		if ($nbWhile > $maxNbWhile) {
#			$dbh=&promsConfig::dbConnect('no_user'); # reconnect
#			$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID"); # Failed
#			$dbh->commit;
#			$dbh->disconnect;
#			die "Aborting quantification: R is taking too long or died before completion";
#		}
#		sleep 30;
#		$pbsError=`head -5 $runDir/PBSerror.txt` if -e "$runDir/PBSerror.txt";
#		$nbWhile++;
#	}


	my %jobParameters=(
		commandFile=>"$runDir/runSSProtQuantification.sh",
		maxMem=>$maxMem,
		numCPUs=>1,
		maxHours=>$maxHours,
		jobName=>"myProMS_SSPA_$quantifID"
	);
	my ($pbsError,$pbsErrorFile)=$cluster{'runJob'}->($runDir,$RcommandString,\%jobParameters);
	if ($pbsError) { # move PBS error message to job error file
		system "cat $pbsErrorFile >> $promsPath{tmp}/quantification/current/$quantifID\_$quantifDate\_error.txt";
	}
}
else { ###>Run job on Web server
	system $RcommandString;
}
sleep 3;

&checkForErrors;

$dbh=&promsConfig::dbConnect('no_user'); # reconnect

####>ERROR Management<####
my $RoutStrg=`tail -3 $runDir/analysisCounting.Rout`;
unless ($RoutStrg=~/proc\.time\(\)/) {
	$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID");
	$dbh->commit;
	$dbh->disconnect;
	$RoutStrg=`tail -20 $runDir/analysisCounting.Rout`;
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



##################################################
####>Parsing R results & importing data in DB<####
##################################################
open(FILESTAT,">>$fileStat");
print FILESTAT "4/4 Parsing results\n";
close FILESTAT;

####>Fetching list of quantification parameters<####
my %quantifParamIDs;
my $sthQP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,P.CODE FROM QUANTIFICATION_PARAMETER P,QUANTIFICATION_METHOD M WHERE P.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND M.ID_QUANTIFICATION_METHOD=(SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID)");
$sthQP->execute;
while (my ($paramID,$code)=$sthQP->fetchrow_array) {
    $quantifParamIDs{$code}=$paramID;
}
$sthQP->finish;

####>Fetching list of possible sets & positions<####
my ($setStrg)=($quantifAnnot=~/SETS=([^\$]+)/);
my %setPosition;
my $tgtPos=0;
foreach my $set (split(/;/,$setStrg)) {
	my $setKey='';
	foreach my $pos (split('\+',$set)) {
		$setKey.='State'.$pos;
	}
	$setPosition{$setKey}=++$tgtPos;
}

my $sthInsProt=$dbh->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_QUANTIFICATION,ID_PROTEIN,ID_QUANTIF_PARAMETER,QUANTIF_VALUE,TARGET_POS) VALUES ($quantifID,?,?,?,?)"); # PK autoincrement

open(RES,"$resultDir/outTable.txt");
while(<RES>) {
	next if $.==1;
	chomp;
	s/"//g;
	my ($protID,$setStrg,$effect,$pVal,@stateMeans)=split(/\t/,$_); #"Protein"	"couple"	"effect"	"pvalue"	"State1.mean"	"Stata2.mean"	"State3.mean"	"State4.mean"	"pvalue_adjusted"
	my $pValAdj=pop(@stateMeans);
	my ($testStrg,$refStrg)=split(/ \/ /,$setStrg);
	my ($bestState,$testState,$refState);
	foreach my $state (split(/ \+ /,$testStrg)) {
		$bestState=$state if (!$bestState || $stateMeans[$setPosition{$state}-1] > $stateMeans[$setPosition{$bestState}-1]);
		$testState=$state if (!$testState || $stateMeans[$setPosition{$state}-1] < $stateMeans[$setPosition{$testState}-1]);
	}
	foreach my $state (split(/ \+ /,$refStrg)) {
		$refState=$state if (!$refState || $stateMeans[$setPosition{$state}-1] > $stateMeans[$setPosition{$refState}-1]);
	}
	(my $setKey=$testStrg)=~s/ \+ //g;

	#>States mean
	foreach my $i (0..$#stateMeans) {
		$sthInsProt->execute($protID,$quantifParamIDs{'MEAN'},$stateMeans[$i],$setPosition{'State'.($i+1)}); # states are ordered from table.txt
	}
	#>Effect
	$sthInsProt->execute($protID,$quantifParamIDs{'EFFECT'},$effect,$setPosition{$setKey});
	#>Best delta (%)
	my $bestMean=$stateMeans[$setPosition{$bestState}-1];
	my $testMean=$stateMeans[$setPosition{$testState}-1];
	my $refMean=$stateMeans[$setPosition{$refState}-1];
	my $deltaPc=100*($testMean-$refMean)/$bestMean;
	$sthInsProt->execute($protID,$quantifParamIDs{'DELTA_PC'},$deltaPc,$setPosition{$setKey});
	#>p-values
	$sthInsProt->execute($protID,$quantifParamIDs{'PVAL'},$pVal,$setPosition{$setKey}) unless $pVal==1;
	$sthInsProt->execute($protID,$quantifParamIDs{'PVAL_ADJ'},$pValAdj,$setPosition{$setKey}) unless $pValAdj==1;
}
close RES;

$sthInsProt->finish;


###> Quantification is finished.
$dbh->do("UPDATE QUANTIFICATION SET STATUS=1 WHERE ID_QUANTIFICATION=$quantifID");
$dbh->commit;

$dbh->disconnect;
#exit; # DEBUG!!!

open(FILESTAT,">>$fileStat");
print FILESTAT "Ended ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
close(FILESTAT);

####>Move all files outside tmp directory
unlink "$runDir/analysisCounting.R";

mkdir "$promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
dirmove($runDir,"$promsPath{quantification}/project_$projectID/quanti_$quantifID");

sleep 2;
unlink $fileStat;


sub checkForErrors {
	#my ($errorFile) = glob "$promsPath{tmp}/quantification/current/*_$quantifDate\_error.txt";
	my $errorFile = "$promsPath{tmp}/quantification/current/$quantifID\_$quantifDate\_error.txt";
	if ($errorFile && -s $errorFile) {
		$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID"); # Failed
		$dbh->commit;
		$dbh->disconnect;
		die "Aborting quantification due to errors.";
	}
}


####>Revision history<####
# 1.2.3 Uses single $cluster{runJob}->() command (PP 11/07/19)
# 1.2.2 Uses quanti_$quantifID instead of quantif_$quantifID (PP 18/05/18)
# 1.2.1 [FIX] Bug due to truncated protein:visibility SQL response string if &gt; 1024 characters (PP 18/01/18)
# 1.2.0 Now relies on &promsConfig::getClusterInfo (PP 30/11/17)
# 1.1.0 Handle multiple peptide counting methods including spectral count (PP 02/08/17)
# 1.0.0 Started (PP 04/08/16)
