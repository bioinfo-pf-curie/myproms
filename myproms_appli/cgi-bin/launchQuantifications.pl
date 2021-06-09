#!/usr/local/bin/perl -w

################################################################################
# launchQuantifications.pl       1.8.18                                        #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Launches multiple quantifications +/- paralellization                        #
# called by selAna4quantification.cgi & startDesignQuantification.cgi          #
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
use POSIX ":sys_wait_h"; # for WNOHANG
use File::Path qw(rmtree); # remove_tree
use promsConfig;
use promsMod;
use promsQuantif;

# exit;
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo('no_user');
my %cluster=&promsConfig::getClusterInfo;

###############################
####>Recovering parameters<####
###############################
my ($multiplicity,$userID,@extraArgs)=@ARGV;
my $MAX_PARALLEL_QUANTIFS=($multiplicity eq 'multi' && $cluster{'on'} && $extraArgs[0]=~/DESIGN/)? $cluster{'maxJobs'} : &promsConfig::getMaxParallelJobs; # multiCluster
#$MAX_PARALLEL_QUANTIFS=4; # tmp


                                                ########################
####################################################>Multi-job call<####################################################
                                                ########################
if ($multiplicity=~/multi/) {
	my @jobList=split(',',$extraArgs[0]);
	my $numJobs=scalar @jobList;
	my %runningJobs;
	$SIG{CHLD} = sub { # sub called when child communicates with parent (occur during parent sleep)
		local $!; # good practice. avoids changing errno.
		while (1) { # while because multiple children could finish at same time
			my $childPid = waitpid(-1,WNOHANG); # WNOHANG: parent doesn't hang during waitpid
			#my $childPid = waitpid(-1,POSIX->WNOHANG); # WNOHANG: parent doesn't hang during waitpid
			last unless ($childPid > 0); # No more to reap.
			delete $runningJobs{$childPid}; # untrack child.
		}
	};

	####>Looping through job list
	my $curJobIdx=-1;

	MAIN_LOOP:while (1) {

		###>Parallel launch of up to $MAX_PARALLEL_QUANTIFS jobs
        my $nJobs = scalar (keys %runningJobs);
		while (scalar (keys %runningJobs) < $MAX_PARALLEL_QUANTIFS) {
			my ($jobDir,$fullQuantifType,$anaIDstrg)=split('#',$jobList[$curJobIdx]);

			##>Also check other user/launches jobs
			my $numProc=`ps -edf | grep launchQuantifications.pl | grep $fullQuantifType | grep single | grep -v grep | wc -l`;
			chomp($numProc);
			last if $numProc >= $MAX_PARALLEL_QUANTIFS;

			##>Ok to go
			$curJobIdx++;
			$anaIDstrg='' unless $anaIDstrg; # for SIN|EMPAI only (PP A single anaID since 02/08/18)
			my $childPid = fork;
			unless ($childPid) { # child here
				system "./launchQuantifications.pl single $userID $jobDir $fullQuantifType $anaIDstrg";
	#system "./launchQuantifications.pl single $userID $jobDir $fullQuantifType $anaIDstrg 2> $promsPath{tmp}/quantification/$jobDir\_errorQuantif.txt";
				exit;
			}
			$runningJobs{$childPid}=$curJobIdx;
			last MAIN_LOOP if $curJobIdx==$#jobList; # no more jobs to launch => break MAIN_LOOP
			sleep 2;
		}

		###>Wait for potential ended job before next MAIN_LOOP
		sleep 60; # $SIG{CHLD} is active during sleep

	}

	###>Wait for last child processes to end
	while (scalar (keys %runningJobs) > 0) {
		sleep 60; # $SIG{CHLD} is active during sleep
	}
	exit;
}


                                                #########################
####################################################>Single-job call<########################################################
                                                #########################
my ($jobDir,$quantifType,$quantItemStrg)=@extraArgs; # $quantItemStrg defined only for SIN|EMPAI|XICCORR or if resuming DESIGN Abundance
my $quantifPath="$promsPath{tmp}/quantification/$jobDir";
my $currentQuantifPath="$promsPath{tmp}/quantification/current";
#my $quantifScript=($quantifType eq 'XIC')? 'runIonintensityQuantification.pl' : ($quantifType eq 'EMPAI')? 'runemPAIQuantification.pl' : ($quantifType eq 'SIN')? 'runSINQuantification.pl' : '';
my %quantifScripts=(
	'DESIGN'=>'runXICProtQuantification.pl',
	'DESIGN:TDA'=>'runXICProtQuantification.pl',
	'DESIGN:DIA'=>'runXICProtQuantification.pl',
	#'DESIGN:MSstats'=>'runSWATHProtQuantification.pl',
	'DESIGN:MSstats'=>'runXICProtQuantification.pl',
	'DESIGN:SSPA'=>'runSSProtQuantification.pl',
	'EMPAI'=>'runemPAIQuantification.pl',
	'HYDRO'=>'runHydrophobicity.pl',
	'SIN'=>'runSINQuantification.pl',
	'XICMCQ'=>'runMassChroQ.pl', # different XIC method
	'XICCORR'=>'correctPepQuantification.pl'
);

#$quantItemStrg='' unless $quantItemStrg;
#my @quantItemList=split(':',$quantItemStrg);
#my %quantifIdList;
my ($quantifID,$quantifItemID); my $numAllBioRep=0; # $quantifItemID != $quantifID for SIN|emPAI|XICCORR only (=)
my $resumeQuantif=$quantItemStrg || 0; # only for DESIGN Abundance
my $previousQuantifStrg=`ls -l | grep quanti_`;
my $jobHistoryType;  # Job type for monitoring
my $extraArgStrg='';
my (@pepQuantifFiles,%analysisList,$numStates,$numReplicates,$bestNumFractions); # needed for cluster run
my %quantifParams;
my $labelType='FREE'; # default
my (@quantifiedModifIDs,@trueQuantifiedModifIDs);
my $isDIA;
if ($quantifType !~ /XICMCQ|XICCORR|HYDRO/) {
	my $dbh=&promsConfig::dbConnect('no_user');
	#>Fetching parameters
	%quantifParams=&promsQuantif::getQuantificationParameters("$quantifPath/quantif_info.txt");
	my $qRootName=quotemeta($quantifParams{'DB'}{'QUANTIF_NAME'}[0]);
	$labelType=$quantifParams{'DB'}{'LABEL'}[0] if $quantifType =~ /^DESIGN/;
	my $quantifAnnotStrg="LABEL=$labelType";
	if ($quantifParams{'DB'}{'ALGO_TYPE'}) {
		if ($quantifParams{'DB'}{'ALGO_TYPE'}[0]=~/^(PEP_|TDA|DIA)/) {
			my $overwritten=0;
			if ($quantifParams{'R'}{'design'}[0] eq 'PEP_INTENSITY') { # Overwrite in case SuperRatio incompatibility
				if ($quantifParams{'DB'}{'ALGO_TYPE'}[0] eq 'PEP_RATIO') {
					$quantifParams{'DB'}{'ALGO_TYPE'}[0]=$quantifParams{'R'}{'design'}[0];
					$quantifParams{'DB'}{'RATIO_TYPE'}[0]='SimpleRatio';
					$overwritten=1;
				}
			}
			$quantifAnnotStrg.="::SOFTWARE=myProMS::ALGO_TYPE=$quantifParams{'DB'}{'ALGO_TYPE'}[0]"; # version added later by run script
			$quantifAnnotStrg.=';changed' if $overwritten;
		}
		elsif ($quantifParams{'DB'}{'ALGO_TYPE'}[0] eq 'MSstats') {
			$quantifAnnotStrg.='::SOFTWARE=DIA/MSstats'; # Version number inserted by runSWATHProtQuantification.pl
		}
		elsif ($quantifParams{'DB'}{'ALGO_TYPE'}[0] eq 'SSPA') {
			$quantifAnnotStrg.="::SOFTWARE=myProMS";
		}
	}
	$quantifAnnotStrg.='::TOP_N='.$quantifParams{'DB'}{'TOP_N'}[0] if $quantifParams{'DB'}{'TOP_N'};
	$quantifAnnotStrg.='::MATCH_PEP='.$quantifParams{'DB'}{'MATCH_PEP'}[0] if $quantifParams{'DB'}{'MATCH_PEP'};
	#my $ratioString='';
	my $conditionStrg='';
	my %anaCondStrg;
	#my $numerator;
	#my $denominator;
#open (OUT,">$promsPath{tmp}/quantification/debug.txt");
#print OUT "'$quantifType'\n";
	#my $quantifiedModifID=($quantifParams{'DB'}{'PTM_SCOPE'})? $quantifParams{'DB'}{'PTM_SCOPE'}[0] : undef;
	if ($quantifParams{'DB'}{'PTM_SCOPE'}) {
		@quantifiedModifIDs=@{$quantifParams{'DB'}{'PTM_SCOPE'}};
		@trueQuantifiedModifIDs=@quantifiedModifIDs;
		shift @trueQuantifiedModifIDs if $quantifiedModifIDs[0] <= 0;
	}
	my $referenceQuantifID=($quantifParams{'DB'}{'REF_QUANTIF'})? $quantifParams{'DB'}{'REF_QUANTIF'}[0] : 0;
	foreach my $param (sort{&promsMod::sortSmart(lc($a),lc($b))} keys %{$quantifParams{'DB'}}) { # sortSmart required for CONDITION_<n>
		next if $param=~/QUANTIF_NAME|LABEL|NUM_REP_|ALGO_TYPE|PTM_SCOPE|TOP_N|MATCH_PEP/;
		if ($param=~/CONDITION_(\d+)/) {
#print OUT ">$param:\n";
			my $condPos=$1;
			$conditionStrg.=';' if $conditionStrg;
			$conditionStrg.=join(',',($quantifParams{'DB'}{"NUM_REP_$condPos"}[0],@{$quantifParams{'DB'}{$param}})); # for each condition: numRep,nameRep1.nameRep2...nameRepN,posRep1.posRep2...posRepN
			next;
		}
		if ($quantifType !~ /DESIGN/ || $param !~ /QUANTITOCOND|DESIGN|STATE/){
			$quantifAnnotStrg.='::' if $quantifAnnotStrg;
			$quantifAnnotStrg.=($param=~/^PEPTIDES_(REF|NORM)$/)? "$param=$quantifParams{DB}{$param}[0]" : "$param=".join(';',@{$quantifParams{'DB'}{$param}}); # PEPTIDES_REF=manual/current or PEPTIDES_NORM=manual only (for now)
		}
		#if ($quantifType eq 'DESIGN' && $param=~/NUMERATOR/){
		#	my $numState=1;
		#	my ($condName)=$dbh->selectrow_array("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=$quantifParams{'DB'}{$param}[0]");
		#	#$ratioString.="$param,$condName,$quantifParams{'DB'}{$param}[0];";
		#	$numerator="$numState,$condName,$quantifParams{'DB'}{$param}[0]";
		#}
		#if ($quantifType eq 'DESIGN' && $param=~/DENOMINATOR/){
		#	my $numState=2;
		#	my ($condName)=$dbh->selectrow_array("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=$quantifParams{'DB'}{$param}[0]");
		#	$denominator="$numState,$condName,$quantifParams{'DB'}{$param}[0]";
		#}
	}
#close OUT;
	my $quantifMethod;
	if ($quantifType !~ /DESIGN/) { # emPAI,SIN,XICCOR
		$quantifAnnotStrg='';
		$quantifMethod=$quantifType;
		$jobHistoryType = uc($quantifType);
	}
	else { # DESIGN(:xxx)
		$quantifMethod=($quantifParams{'DB'}{'ALGO_TYPE'}[0] eq 'SSPA')? 'SSPA' : ($quantifParams{'DB'}{'ALGO_TYPE'}[0] eq 'TnPQ')? 'TNPQ' : ($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None')? 'PROT_ABUNDANCE' : 'PROT_RATIO_PEP';
		if ($quantifParams{'DB'}{'ALGO_TYPE'}[0] =~ /DIA|PEP_INTENSITY/) {
			if ($quantifParams{'DB'}{'RATIO_TYPE'}[0] && $quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None') {
				$jobHistoryType = 'DESIGN:' . $quantifParams{'DB'}{'ALGO_TYPE'}[0] . ':Abund';  # Char nb limitation
			} else {
				$jobHistoryType = 'DESIGN:' . $quantifParams{'DB'}{'ALGO_TYPE'}[0];
			}
		} elsif ($quantifParams{'DB'}{'ALGO_TYPE'}[0] =~ /SSPA|MSstats|TDA|PEP_RATIO/) {
			$jobHistoryType = 'DESIGN:' . $quantifParams{'DB'}{'ALGO_TYPE'}[0];
		} else {
			$jobHistoryType = 'DESIGN';
		}
		my @states=@{$quantifParams{'DB'}{'STATES'}};
		$numStates=scalar @states;

		$quantifAnnotStrg.='::NUM_TRUE_USED=1' if ($quantifParams{'DB'}{'NUM_TRUE_USED'} && $quantifParams{'DB'}{'NUM_TRUE_USED'}[0]); # PEP_INTENSITY	
		if ($quantifMethod ne 'SSPA' && $quantifParams{'DB'}{'RATIO_TYPE'}[0] ne 'None') {
			#$quantifAnnotStrg.='::MEAN_STATE=1' if ($quantifParams{'R'}{'design'} && $quantifParams{'R'}{'design'}[0] eq 'PEP_INTENSITY' && $quantifParams{'DB'}{'ALGO_TYPE'}[0] ne 'MSstats');
			$quantifAnnotStrg.='::MEAN_STATE=1' if $quantifParams{'DB'}{'ALGO_TYPE'}[0]=~/^(PEP_INTENSITY|TDA|DIA)/;
			$quantifAnnotStrg.="::RATIOS=";
			foreach my $x (0..$#states-1) {
				foreach my $y ($x+1..$#states) {
					$quantifAnnotStrg.="#$states[$y]/#$states[$x];";
				}
				last if ($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'SuperRatio' || $quantifParams{'DB'}{'SINGLE_REF'});
			}
			if ($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'SuperRatio') {
				foreach my $x (1..$#states-1) { # ratios of ratios
					foreach my $y ($x+1..$#states){
						$quantifAnnotStrg.="#$states[$y]%$y/#$states[$x]%$x;";
					}
				}
			}
			#chop($quantifAnnotStrg); # remove the last ';'
			$quantifAnnotStrg=~s/;\Z//; # remove trailing ';'
		}
		#my %labelName;
		#foreach my $anaToCondStrg (@{$quantifParams{'DB'}{'QUANTITOCOND'}}){
		#	my ($quantiID,$targetPos,$anaID,$condID)=split(/:/,$anaToCondStrg);
		#	push @{$labelName{$condID}} , "#$quantiID:#$anaID";
		#}

		my $sthAQ=$dbh->prepare("SELECT COUNT(*) FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=?");
		my %quantifNumAna;
		$numReplicates=1;
		$bestNumFractions=1;
		$quantifAnnotStrg.="::STATES="; #  numRep,#obsID1:#quantifID1(+#obsID2:#quantifID2)(.#obsID3:#quantifID3(+#obsID4:#quantifID4)),#condID;<data for state2>(;...)
		foreach my $condID (@states) {
			#my $curGroup=1000; # to make sure not overlap with user-set fracgroups (~ < 20)
			#my %fracGrpObs;
			#foreach my $obsData (@{$quantifParams{'DB'}{"QUANTITOCOND_$condID"}}) {
			#	my ($obsID,$parQuantifID,$anaID,$targetPos,$fracGroup)=split(/:/,$obsData);
			#	my $fgrp=($fracGroup)? $fracGroup : ++$curGroup;
			#	push @{$fracGrpObs{$fgrp}},"#$obsID:#$parQuantifID:#$anaID:$targetPos";
			#}
			#my $numRep=scalar keys %fracGrpObs;
			#$quantifAnnotStrg.="$numRep,";
			#my $firstGr=1;
			#foreach my $fgrp (sort{$a<=>$b} keys %fracGrpObs) {
			#	if ($firstGr) {$firstGr=0;}
			#	else {$quantifAnnotStrg.='.';}
			#	$quantifAnnotStrg.=join('+',@{$fracGrpObs{$fgrp}}); # sum up
			#}
			#$quantifAnnotStrg.=",#$condID;";
			my (%replicateHierarchy,%fracGrBioRep,%techRepPosition,%curTechRepPos);
			my $curTechRepGr=1000; # to make sure not overlap with user-set fracgroups (~ < 20)
			my $curFracGroup=1000;
			foreach my $obsData (@{$quantifParams{'DB'}{"QUANTITOCOND_$condID"}}) {
				my ($obsID,$parQuantifID,$anaID,$targetPos,$fracGroup,$bioRep)=split(/:/,$obsData);
				$fracGroup=++$curFracGroup unless $fracGroup;
				#$bioRep=++$curTechRepGr unless $bioRep;
				unless ($bioRep) {
					$bioRep=($fracGrBioRep{$fracGroup})? $fracGrBioRep{$fracGroup} : ++$curTechRepGr;
				}
				$fracGrBioRep{$fracGroup}=$bioRep;
				unless ($techRepPosition{$bioRep}) {
					%{$techRepPosition{$bioRep}}=();
					$curTechRepPos{$bioRep}=0;
				}
				unless ($techRepPosition{$bioRep}{$fracGroup}) { # tech rep pos in tech group (uses 1st one declared for fraction grp)
					$techRepPosition{$bioRep}{$fracGroup}=++$curTechRepPos{$bioRep};
				}
				push @{$replicateHierarchy{$bioRep}{ $techRepPosition{$bioRep}{$fracGroup} }},"#$obsID:#$parQuantifID:#$anaID:$targetPos";

				#>To estimate necessary cluster ressources
				if ($cluster{'on'}) {
					if ($quantifParams{'DB'}{'ALGO_TYPE'}[0] eq 'SSPA') {
						$analysisList{$anaID}=1;
					}
					else {
						unless ($quantifNumAna{$parQuantifID}) {
							$sthAQ->execute($parQuantifID);
							($quantifNumAna{$parQuantifID})=$sthAQ->fetchrow_array;
							unless (defined $isDIA) {
								($isDIA)=$dbh->selectrow_array("SELECT COUNT(*) FROM QUANTIFICATION Q,QUANTIFICATION_METHOD QM WHERE Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=$parQuantifID AND QM.CODE='DIA'");
							}
						}
						push @pepQuantifFiles,[$parQuantifID,$anaID,$targetPos,$quantifNumAna{$parQuantifID}];
					}
				}
			}
			my $numBioRep=scalar keys %replicateHierarchy;
			$quantifAnnotStrg.="$numBioRep,";
			$numAllBioRep+=$numBioRep;
			my $firstBioRep=1;
			foreach my $bioRep(sort{$a<=>$b} keys %replicateHierarchy) {
				if ($firstBioRep) {$firstBioRep=0;}
				else {$quantifAnnotStrg.='.';}
				my $firstTechRep=1;
				foreach my $techRep (sort{$a<=>$b} keys %{$replicateHierarchy{$bioRep}}) {
					if ($firstTechRep) {$firstTechRep=0;}
					else {$quantifAnnotStrg.='&';} # &: technical replicate separator
					$quantifAnnotStrg.=join('+',@{$replicateHierarchy{$bioRep}{$techRep}}); # sum up
					$numReplicates++; # counts all techReps
					my $numFrac=scalar @{$replicateHierarchy{$bioRep}{$techRep}}; 
					$bestNumFractions=$numFrac if $numFrac > $bestNumFractions;
				}
			}
			$quantifAnnotStrg.=",#$condID;";


			#if ($labelName{$condID}){
			#	$quantifAnnotStrg.=scalar(@{$labelName{$condID}}).',';
			#	$quantifAnnotStrg.=join('.', @{$labelName{$condID}} );
			#	$quantifAnnotStrg.=",#$condID;";
			#}
		}
		chop($quantifAnnotStrg); # remove the last ';'
		$sthAQ->finish;
	}
	#my $qQuantifAnnotStrg=$dbh->quote($quantifAnnotStrg);
	#my $quantifMethod=($quantifType=~/SILAC|ITRAQ|TMT/)? 'PROT_RATIO_PEP' : ($quantifType eq 'DESIGN')? $quantifParams{'DB'}{'QUANTIFICATION_METHOD'}[0] : $quantifType;
	my ($quantifMethodID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE=UCASE('$quantifMethod')");
	#my $sthN=$dbh->prepare("SELECT COUNT(*) FROM QUANTIFICATION Q,ANA_QUANTIFICATION A WHERE Q.ID_QUANTIFICATION=A.ID_QUANTIFICATION AND A.ID_ANALYSIS=? AND NAME LIKE '$qRootName%'");
	my $focus=($quantifType eq 'SIN' || $quantifType eq 'XIC')? 'peptide' : 'protein';
	#my $sthInsQ=$dbh->prepare("INSERT INTO QUANTIFICATION (ID_QUANTIFICATION_METHOD,ID_MODIFICATION,NAME,FOCUS,QUANTIF_ANNOT,STATUS,UPDATE_DATE,UPDATE_USER) VALUES ($quantifMethodID,?,?,'$focus',?,-1,NOW(),'$userID')");
	my $sthInsQ=$dbh->prepare("INSERT INTO QUANTIFICATION (ID_QUANTIFICATION_METHOD,NAME,FOCUS,QUANTIF_ANNOT,STATUS,UPDATE_DATE,UPDATE_USER) VALUES ($quantifMethodID,?,'$focus',?,-1,NOW(),'$userID')"); # NO LONGER USE ID_MODIFICATION
	my $sthInsModQ=$dbh->prepare("INSERT INTO MULTIMODIF_QUANTIFICATION (ID_QUANTIFICATION,ID_MODIFICATION,MODIF_RANK) VALUES (?,?,?)");
	my $sthInsAnaQ=$dbh->prepare("INSERT INTO ANA_QUANTIFICATION (ID_QUANTIFICATION,ID_ANALYSIS) VALUES (?,?)");
	my $sthInsParQ=$dbh->prepare("INSERT INTO PARENT_QUANTIFICATION (ID_QUANTIFICATION,ID_PARENT_QUANTIFICATION,PAR_FUNCTION) VALUES (?,?,?)");
	#@quantItemList=split(':',$quantItemStrg);

	if ($quantifType !~ /DESIGN/) { # SIN|EMPAI
		my $analysisID=$quantifItemID=$quantItemStrg;
		$extraArgStrg=$analysisID;
		$extraArgStrg.=" $userID" if $quantifType eq 'SIN';
		my ($anaName)=$dbh->selectrow_array("SELECT NAME FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
		my $quantifName=$quantifParams{'DB'}{'QUANTIF_NAME'}[0].'_'.$anaName; # Auto_quantifType_anaName
		my ($lastName)=$dbh->selectrow_array("SELECT NAME FROM QUANTIFICATION Q,ANA_QUANTIFICATION A WHERE Q.ID_QUANTIFICATION=A.ID_QUANTIFICATION AND A.ID_ANALYSIS=$analysisID AND FOCUS='$focus' AND NAME LIKE '$qRootName%' ORDER BY Q.ID_QUANTIFICATION DESC LIMIT 1");
		if ($lastName) {
			$quantifName.=($lastName=~/#(\d+)\Z/)? " #".($1+1) : " #2";
		}
		$sthInsQ->execute($quantifName,$quantifAnnotStrg);
		$quantifID=$dbh->last_insert_id(undef, undef,'QUANTIFICATION','ID_QUANTIFICATION');
		$sthInsAnaQ->execute($quantifID,$analysisID);

		###> Update quantif_info.txt in case needed later (PP 30/10/18)
		open (INFO,">>$quantifPath/quantif_info.txt");
		print INFO "QUANTIFICATIONS:\n$quantifID";
		close INFO;
	}
	else {# DESIGN updates
		my $quantifName=$quantifParams{'DB'}{'QUANTIF_NAME'}[0];
		if ($resumeQuantif) { # => Quantif already recorded in DB
			$quantifID=$resumeQuantif;
			#warn "Resuming quantification #$quantifID"; # DEBUG!!!!!!!!!!!!!
		}
		else {
			my $numQuantifiedModifs=scalar @quantifiedModifIDs;
			$sthInsQ->execute($quantifName,$quantifAnnotStrg);# Add the new created
			$quantifID=$dbh->last_insert_id(undef, undef, 'QUANTIFICATION', 'ID_QUANTIFICATION');
			if ($numQuantifiedModifs) {
				my $modifRank=0;
				foreach my $modID (@quantifiedModifIDs) { # already sorted in quantif_info.txt
					$modifRank++;
					my $usedModifRank=($numQuantifiedModifs==1)? undef : $modifRank;
					$sthInsModQ->execute($quantifID,$modID,$usedModifRank);
				}
			}

			my $designID=$quantifParams{'DB'}{'ID_DESIGN'}[0];
			$dbh->do("UPDATE QUANTIFICATION SET ID_DESIGN=$designID WHERE ID_QUANTIFICATION=$quantifID");

			my $sthInsEQ=$dbh->prepare("INSERT INTO EXPCONDITION_QUANTIF (ID_QUANTIFICATION,ID_EXPCONDITION) VALUES ($quantifID,?)");

			my (%usedParentQ,%usedAna);
			foreach my $condID (@{$quantifParams{'DB'}{'STATES'}}) {
				$sthInsEQ->execute($condID); # EXPCONDITION_QUANTIF
				foreach my $obsData (@{$quantifParams{'DB'}{"QUANTITOCOND_$condID"}}) {
					my ($obsID,$parQuantifID,$anaID,$targetPos,$fracGroup)=split(/:/,$obsData);
					$sthInsParQ->execute($quantifID,$parQuantifID,undef) unless $usedParentQ{$parQuantifID}; # PARENT_QUANTIFICATION
					$usedParentQ{$parQuantifID}=1;
					#$sthObs->execute($obsID);
					#my ($anaID)=$sthObs->fetchrow_array;
					$sthInsAnaQ->execute($quantifID,$anaID) unless $usedAna{$anaID}; # ANA_QUANTIFICATION
					$usedAna{$anaID}=1;
				}
			}

			#my @anaToCond=@{$quantifParams{'DB'}{'QUANTITOCOND'}};
			#my (%expAssoToQuanti,%assoToParent);
			#foreach my $anaToCondStrg (@anaToCond){
			#	my ($quantiID,$targetPos,$anaID,$condID)=split(/:/,$anaToCondStrg); # tragetPos=0 for label-free
			#	$sthInsEQ->execute($condID,$quantifID) unless $expAssoToQuanti{$condID};# Associate the new TNPQ quantification to the expcondition
			#	$expAssoToQuanti{$condID}=1;
			#	$sthInsParQ->execute($quantifID,$quantiID) unless $assoToParent{$quantiID};# Create new PARENT_QUANTIFICATION items and avoid duplicated entry
			#	$assoToParent{$quantiID}=1;
			#	$sthInsAQ->execute($anaID);
			#}

			$sthInsEQ->finish;

			# Insert reference quantif if any (in cas PTM quantif)
			$sthInsParQ->execute($quantifID,$referenceQuantifID,'reference') if $referenceQuantifID;

			# Update at the last minute the status of the desired value for R: whether Ratio or TnPQ
			#my $rType=$quantifParams{'R'}{'quantification.method'}[0];
			#if ($rType eq 'Ratio'){
			#	$dbh->do("UPDATE QUANTIFICATION SET ID_QUANTIFICATION_METHOD=(SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='PROT_RATIO_PEP') WHERE ID_QUANTIFICATION=$quantifID");
			#}

			#$quantifID=$maxQuantifID; # UPDATE $quantifID so as to launch properly the TNPQ script

			###> Update quantif_info.txt to make it work for watchQuantifications.cgi
			open (INFO,">>$quantifPath/quantif_info.txt");
			print INFO "QUANTIFICATIONS:\n$quantifID";
			close INFO;

		}
		$quantifItemID=$quantifID;

		#>Generating job list
		#@quantItemList=($quantifID);
		unlink "$currentQuantifPath/$jobDir\_request.flag" if -e "$currentQuantifPath/$jobDir\_request.flag"; # in case multi-job launch
		open(WAIT,">$currentQuantifPath/$quantifID\_$jobDir\_wait.flag"); # flag file (created by selAna4Quantification for other quantif types than DESIGN* & XICMCQ)
		print WAIT '#';
		close WAIT;
		#$quantifIdList{$quantifID}=$quantifID;
	}

	$sthInsQ->finish;
	$sthInsModQ->finish;
	$sthInsAnaQ->finish;
	$sthInsParQ->finish;

	$dbh->commit;
	$dbh->disconnect;
}
elsif ($quantifType eq 'XICMCQ') {
	$jobHistoryType = $quantifType;
	####>Getting the parameters from quantif_info.txt<####
	open (INFO,"$quantifPath/quantif_info.txt");
	my $section='';
	my (%params);
	while (<INFO>) {
		if (/^PARAMETERS:/) {
			$section='parameters';
			next;
		}
		###last if /^ANALYSES/; not listed anymore
		if ($section eq 'parameters' && /^\w+/) {
			chomp;
			my ($paramName,$paramValue)=split(/\t\t/,$_);
			if ($paramName eq 'MZXML'){
				my ($anaID,$mzXMLFile)=split(/:/,$paramValue);
				$params{$paramName}{$anaID}=$mzXMLFile;
				$analysisList{$anaID}=1;
			}
			else{
				$params{$paramName}=$paramValue;
			}
		}
	}
	close INFO;

	my $dbh=&promsConfig::dbConnect('no_user');
	###> Update Quantification
	my ($code)=($params{'CHANNELS'})? 'SILAC' :'XIC';
	my ($quantifMethodID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='$code'");

	my $sthInsQ=$dbh->prepare("INSERT INTO QUANTIFICATION (ID_QUANTIFICATION_METHOD,NAME,FOCUS,STATUS,QUANTIF_ANNOT,UPDATE_DATE,UPDATE_USER) VALUES ($quantifMethodID,?,'peptide',-1,?,NOW(),'$userID')");
	my $sthInsAnaQ=$dbh->prepare("INSERT INTO ANA_QUANTIFICATION (ID_QUANTIFICATION,ID_ANALYSIS,QUANTIF_FILE,IS_REFERENCE) VALUES (?,?,?,?)");

	my %quantifParams=&promsQuantif::getQuantificationParameters("$quantifPath/quantif_info.txt");
	my $qRootName=quotemeta($quantifParams{'DB'}{'QUANTIF_NAME'}[0]); # Name given by user
	my $anaIDS=join(',', keys (%{$params{'MZXML'}}) );
	my $quantifName=$quantifParams{'DB'}{'QUANTIF_NAME'}[0];
	my ($lastName)=$dbh->selectrow_array("SELECT NAME FROM QUANTIFICATION Q,ANA_QUANTIFICATION A WHERE Q.ID_QUANTIFICATION=A.ID_QUANTIFICATION AND A.ID_ANALYSIS IN ($anaIDS) AND FOCUS='peptide' AND NAME LIKE '$qRootName%' ORDER BY Q.ID_QUANTIFICATION DESC LIMIT 1");
	if ($lastName) {
		$quantifName.=($lastName=~/#(\d+)\Z/)? " #".($1+1) : " #2";
	}
	my $softwareStrg="::SOFTWARE=MCQ;?"; # Updated by runMassChroQ.pl with real version number
	my $quantifAnnotStrg=($params{'CHANNELS'})? "LABEL=SILAC$softwareStrg::$params{CHANNELS}" : "LABEL=FREE$softwareStrg";
	foreach my $param (sort{lc($a) cmp lc($b)} keys %{$quantifParams{'DB'}}) {
		$quantifAnnotStrg.="::$param=".join(';',@{$quantifParams{'DB'}{$param}}) unless ($param eq 'MZXML' || $param eq 'CHANNELS'|| $param eq 'QUANTIF_NAME');
	}
	$sthInsQ->execute($quantifName,$quantifAnnotStrg);
	$quantifID=$dbh->last_insert_id(undef, undef, 'QUANTIFICATION', 'ID_QUANTIFICATION');
	$sthInsQ->finish;
	$quantifItemID=$quantifID;

	foreach my $anaID (keys %{$params{'MZXML'}}) {
		my $isReference=($anaID==$params{'REFERENCE'})? 1 : 0;
		$sthInsAnaQ->execute($quantifID,$anaID,$params{MZXML}{$anaID},$isReference);
	}
	$sthInsAnaQ->finish;
	$dbh->commit;
	$dbh->disconnect;

	###> Update quantif_info.txt to make it work for watchQuantifications.cgi
	open (INFO,">>$quantifPath/quantif_info.txt");
	print INFO "QUANTIFICATIONS:\n$quantifID";
	close INFO;

#@quantItemList=($quantifID);
	open(WAIT,">$currentQuantifPath/$quantifID\_$jobDir\_wait.flag"); # flag file (created by selAna4Quantification for other quantif types than DESIGN*, XICMCQ & HYDRO)
	print WAIT '#';
	close WAIT;
#$quantifIdList{$quantifID}=$quantifID;
}
elsif ($quantifType eq 'XICCORR') {
	$jobHistoryType = $quantifType;
	$quantifItemID=$quantItemStrg; # "anaID.parentQuantifID"
	$quantifID=(split(/\./,$quantItemStrg))[1];

	####>Getting the parameters from quantif_info.txt<####
	my $productID=0;
	my $keepOld=0;
	open (INFO,"$quantifPath/quantif_info.txt");
	while (<INFO>) {
		if (/^ID_PRODUCT=(\d+)/) {
			$productID=$1;
		}
		elsif (/^KEEP_OLD=(\d)/) {
			$keepOld=$1;
		}
	}
	close INFO;
	$extraArgStrg="$quantifItemID $productID $keepOld $userID";
}
elsif ($quantifType eq 'HYDRO') {
	$jobHistoryType = $quantifType;

	# Get parameters from quantif_info.txt
	open(INFO, "$quantifPath/quantif_info.txt");
	my $section='';
	my %params;
	while (<INFO>) {
		if (/^PARAMETERS:/) {
			$section = 'parameters';
			next;
		} elsif (/^ANALYSES/) {
			$section = 'analyses';
			next;
		}
		if ($section eq 'parameters' && /^\w+/) {
			chomp;
			my ($paramName, $paramValue) = split(/\t\t/, $_);
			$params{$paramName} = $paramValue;
		} elsif ($section eq 'analyses' && /^\w+/) {
			chomp;
			my ($paramName, $paramValue) = split(/\t\t/, $_);
			if ($paramName eq 'ANA_PEPQUANTIF'){
				my ($anaID, $pepQuantifID) = split(/:/, $paramValue);
				$params{$paramName}{$anaID} = $pepQuantifID;
				$analysisList{$anaID} = 1;
			}
		}
	}
	close INFO;

	my $dbh=&promsConfig::dbConnect('no_user');
	###> Update Quantification
	my $code = "HYDRO";
	my ($quantifMethodID) = $dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='$code'");

	my $sthInsQ = $dbh->prepare("INSERT INTO QUANTIFICATION (ID_QUANTIFICATION_METHOD, NAME, FOCUS, STATUS, QUANTIF_ANNOT, UPDATE_DATE, UPDATE_USER) VALUES ($quantifMethodID, ?, 'peptide', -1, ?, NOW(), '$userID')");
	my $sthInsAnaQ = $dbh->prepare("INSERT INTO ANA_QUANTIFICATION (ID_QUANTIFICATION, ID_ANALYSIS) VALUES (?, ?)");
	my $sthInsParQ = $dbh->prepare("INSERT INTO PARENT_QUANTIFICATION (ID_QUANTIFICATION, ID_PARENT_QUANTIFICATION, PAR_FUNCTION) VALUES (?, ?, ?)");

	my %quantifParams = &promsQuantif::getQuantificationParameters("$quantifPath/quantif_info.txt");
	my $qRootName = quotemeta($quantifParams{'DB'}{'QUANTIF_NAME'}[0]); # Name given by user
	my $anaIDs = join(',', keys %{$params{'ANA_PEPQUANTIF'}});
	my ($lastName) = $dbh->selectrow_array("SELECT NAME FROM QUANTIFICATION Q, ANA_QUANTIFICATION A WHERE Q.ID_QUANTIFICATION = A.ID_QUANTIFICATION AND A.ID_ANALYSIS IN ($anaIDs) AND FOCUS = 'peptide' AND NAME LIKE '$qRootName%' ORDER BY Q.ID_QUANTIFICATION DESC LIMIT 1");
	my $quantifName = $quantifParams{'DB'}{'QUANTIF_NAME'}[0];
	if ($lastName) {
		$quantifName .= ($lastName =~ /#(\d+)\Z/)? " #".($1+1) : " #2";
	}

	my $softwareStrg="::SOFTWARE=myProMS";
	my $quantifAnnotStrg = "LABEL=FREE$softwareStrg";
	foreach my $param (sort{lc($a) cmp lc($b)} keys %{$quantifParams{'DB'}}) {
		$quantifAnnotStrg .= "::$param=" . join(';', @{$quantifParams{'DB'}{$param}}) unless ($param eq 'QUANTIF_NAME');
	}
	$quantifAnnotStrg .= "::ANA_PEPQUANTIF=";
	$quantifAnnotStrg .= join(';', map { if ($params{'ANA_PEPQUANTIF'}{$_} =~ /\d/) { "#$_,#$params{'ANA_PEPQUANTIF'}{$_}" } else { "#$_,$params{'ANA_PEPQUANTIF'}{$_}" } } (sort {$a <=> $b} keys %{$params{'ANA_PEPQUANTIF'}}));

	$sthInsQ->execute($quantifName, $quantifAnnotStrg);
	$quantifID = $dbh->last_insert_id(undef, undef, 'QUANTIFICATION', 'ID_QUANTIFICATION');
	$sthInsQ->finish;
	$quantifItemID = $quantifID;

	my (%usedParentQ,%usedAna);
	foreach my $anaID (keys %{$params{'ANA_PEPQUANTIF'}}) {
		my $pepQuantifID = $params{'ANA_PEPQUANTIF'}{$anaID};
		$sthInsAnaQ->execute($quantifID, $anaID) unless $usedAna{$anaID};
		$usedAna{$anaID} = 1;
		unless ($pepQuantifID eq "no_alignment" || $usedParentQ{$pepQuantifID}) {
			$sthInsParQ->execute($quantifID, $pepQuantifID, undef);
			$usedParentQ{$pepQuantifID} = 1;
		}
	}
	$sthInsAnaQ->finish;
	$sthInsParQ->finish;
	$dbh->commit;
	$dbh->disconnect;

	###> Update quantif_info.txt to make it work for watchQuantifications.cgi
	open(INFO, ">>$quantifPath/quantif_info.txt");
	print INFO "QUANTIFICATIONS:\n$quantifID";
	close INFO;

	open(WAIT,">$currentQuantifPath/$quantifID\_$jobDir\_wait.flag"); # flag file (created by selAna4Quantification for other quantif types than DESIGN*, XICMCQ & HYDRO)
	print WAIT '#';
	close WAIT;
}
#exit; # DEBUG

##########################################
####>Launching quantification process<#### only 1 child now (PP 03/08/18)
##########################################
my $dbh=&promsConfig::dbConnect('no_user');
my $projectID=&promsMod::getProjectID($dbh,$quantifID,'QUANTIFICATION');

my $quantiDBInfo = "ID_QUANTIFICATION=$quantifID;TYPE=$jobHistoryType";
$quantiDBInfo .= ";ID_ANALYSIS=$quantItemStrg" if ($quantItemStrg && $quantifType =~ /SIN|EMPAI/);
my $logFilePath = ($quantItemStrg) ? "$quantifPath/status\_$quantItemStrg.out" : "$quantifPath/status\_$quantifID.out";
my $errorFilePath = ($quantItemStrg) ? "$currentQuantifPath/$quantItemStrg\_$jobDir\_error.txt" : "$currentQuantifPath/$quantifID\_$jobDir\_error.txt";
$dbh->do("INSERT INTO JOB_HISTORY (ID_JOB, ID_USER, ID_PROJECT, TYPE, JOB_STATUS, FEATURES, SRC_PATH, LOG_PATH, ERROR_PATH, STARTED) VALUES('$jobDir', '$userID', $projectID, 'Quantification [$jobHistoryType]', 'Queued', '$quantiDBInfo', '$quantifPath', '$logFilePath', '$errorFilePath', NOW())");
$dbh->commit;

if ($cluster{'on'} && ($quantifType =~ /^DESIGN/ || $quantifType eq 'XICMCQ' || $quantifType eq 'HYDRO')) { #'DESIGN' or 'DESIGN:MSstats' or 'DESIGN:SSPA' or 'DESIGN:TDA'
	
	###>Estimating ressources required<###
	my ($maxHours,$maxMem,$numCPU,$jobName)=(24,1,1,"myProMS_protQuant_$quantifType\_$quantifID"); # default

	if ($quantifType =~ /DESIGN/ && $quantifType ne 'DESIGN:SSPA') {
		my $refSomeQuantFiles=\@pepQuantifFiles; # default
		my $numAllFiles=scalar @pepQuantifFiles;
		if ($numAllFiles > 15) { # scan only 15 representative files
			my $mid1=int(($numAllFiles-1)/2)-2;
			my $mid2=$mid1+4;
			@{$refSomeQuantFiles}=@pepQuantifFiles[0..4,$mid1..$mid2,-5..-1]; # select 5x3 files at begining, middle and end of list
		}
		my ($bestNumLines,$bestAnaID)=(0,0); # $numPepValues,
		# if ($quantifType eq 'DESIGN') {
		# 	my %usedQuantifFile;
		# 	foreach my $resQuantInfo (@{$refSomeQuantFiles}) {
		# 		my ($parQuantifID,$anaID,$targetPos,$numAna)=@{$resQuantInfo};
		# 		my $fileKey=$parQuantifID.':'.$targetPos;
		# 		unless ($usedQuantifFile{$fileKey}) {
		# 			my $pepQuantifFile=($targetPos)? "peptide_quantification_$targetPos.txt" : 'peptide_quantification.txt';
		# 			my $numLines=`wc -l $promsPath{quantification}/project_$projectID/quanti_$parQuantifID/$pepQuantifFile`;
		# 			#chomp($numLines);
		# 			$numLines=~s/^(\d+).*/$1/;
		# 			$usedQuantifFile{$fileKey}=int(0.5 + $numLines/$numAna);
		# 			if ($usedQuantifFile{$fileKey} > $bestNumLines) {
		# 				$bestAnaID=$anaID;
		# 				$bestNumLines=$usedQuantifFile{$fileKey};
		# 			}
		# 		}
		# 		$numPepValues+=$usedQuantifFile{$fileKey};
		# 	}
		# }
		# els
		if ($quantifType =~ /DESIGN:MSstats/) {
			foreach my $resQuantInfo (@{$refSomeQuantFiles}) {
				my ($parQuantifID,$anaID,$targetPos,$numAna)=@{$resQuantInfo};
				my $numLines=`wc -l $promsPath{quantification}/project_$projectID/quanti_$parQuantifID/swath_ana_$anaID.txt`;
				#chomp($numLines);
				$numLines=~s/^(\d+).*/$1/;
				$numLines*=2 if $quantifParams{'DB'}{'PEPTIDES'}[0] ne 'unique'; # more peptides => more memory
				# $numPepValues+=($numLines*4); #int($numLines/$numAna);
				if ($numLines > $bestNumLines) {
					$bestAnaID=$anaID;
					$bestNumLines=$numLines;
				}
			}
		}
		else { # myProMS algo   ($quantifType =~ /DESIGN:(TDA|DIA)/)
			my $sthNumPep=$dbh->prepare("SELECT COUNT(*) FROM PEPTIDE WHERE ID_ANALYSIS=?");
			foreach my $resQuantInfo (@{$refSomeQuantFiles}) {
				my ($parQuantifID,$anaID,$targetPos,$numAna)=@{$resQuantInfo};
				$sthNumPep->execute($anaID);
				my ($numPep)=$sthNumPep->fetchrow_array;
				# $numPepValues=$numPep if $numPep > $numPepValues; # only used for maxHour... not very useful
				if ($numPep > $bestNumLines) {
					$bestAnaID=$anaID;
					$bestNumLines=$numPep;
				}
			}
			$sthNumPep->finish;
			$bestNumLines*=5 if ($isDIA && $quantifParams{'DB'}{'RATIO_TYPE'}[0] ne 'None'); # only for myProMS ratio
		}
		# if ($numAllFiles > 15 && $quantifType =~ /^(DESIGN|DESIGN:MSstats)$/) {
		# 	$numPepValues*=($numAllFiles/15); # extend to all files <=> $numPepValues=$numAllFiles * $numPepValues/15;
		# }

		#>Adjust $bestNumLines if PTM quantif
		if (scalar @trueQuantifiedModifIDs) {
			my ($bestNumPep)=($quantifType =~ /DESIGN:MSstats/)? $dbh->selectrow_array("SELECT COUNT(*) FROM PEPTIDE WHERE ID_ANALYSIS=$bestAnaID") : $bestNumLines;
			my $ptmIdStrg=join(',',@trueQuantifiedModifIDs);
			my $numUsedPep=0;
			#my $dbh=&promsConfig::dbConnect('no_user');
			if ($quantifParams{'DB'}{'PTM_POS'} && $quantifParams{'DB'}{'PTM_POS'}[1] ne 'ambiguous' && $quantifParams{'DB'}{'PTM_POS'}[0] ne 'ambiguous') {
				my ($ptmProbSoftCode,$ptmProbThreshold)=$quantifParams{'DB'}{'PTM_POS'}[0]=~/(\w+):(\d+)/;
				if ($ptmProbSoftCode eq 'PRS') {
					my $sthModPep=$dbh->prepare("SELECT DATA FROM PEPTIDE WHERE ID_ANALYSIS=$bestAnaID AND DATA LIKE '%PRS=%'"); # AND VALID_STATUS > 0
					$sthModPep->execute;
					while (my ($pepData)=$sthModPep->fetchrow_array) {
						my ($ptmPosProb)=($pepData=~/PRS=\d;([^;]+);/);
						$numUsedPep++ if $ptmPosProb >= $ptmProbThreshold;
					}
				}
				else {
					$ptmProbThreshold/=100 if $ptmProbSoftCode =~ /^(MQ|SPC|PTMRS)$/; # 0-1
					my %modifQuantif;
					foreach my $modID (@trueQuantifiedModifIDs) {$modifQuantif{$modID}=1;}
					my $sthModPep=$dbh->prepare("SELECT P.ID_PEPTIDE,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',POS_STRING SEPARATOR '&'),GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,'%',COALESCE(PM.REF_POS_STRING,'') SEPARATOR '&') FROM PEPTIDE P INNER JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
													WHERE ID_MODIFICATION IN ($ptmIdStrg) AND ID_ANALYSIS=$bestAnaID AND VALID_STATUS > 0 GROUP BY P.ID_PEPTIDE"); #
					$sthModPep->execute;
					while (my ($pepID,$varModCode,$varProbStrg)=$sthModPep->fetchrow_array) {
						my $ptmPosProb=&promsMod::extractPtmPosProbability(\%modifQuantif,$varModCode,$varProbStrg);
						$numUsedPep++ if $ptmPosProb >= $ptmProbThreshold;
					}
					#>Account for proportion of ghosts (No proba)
					my ($numGhosts)=$dbh->selectrow_array("SELECT COUNT(*) FROM PEPTIDE WHERE ID_ANALYSIS=$bestAnaID AND VALID_STATUS=0");
					$numUsedPep*=$bestNumPep/($bestNumPep-$numGhosts);
				}
				# #>Add proportion of ghosts <=== WRONG here!!!! Ghosts should be counted because they also have PRS data (PP 16/02/21)
				# my ($numGhosts)=$dbh->selectrow_array("SELECT COUNT(*) FROM PEPTIDE WHERE ID_ANALYSIS=$bestAnaID AND VALID_STATUS=0");
				# $numUsedPep*=$bestNumLines/($bestNumLines-$numGhosts);
			}
			else { # just count modified ions <- ambiguous allowed
				($numUsedPep)=$dbh->selectrow_array("SELECT COUNT(DISTINCT P.ID_PEPTIDE) FROM PEPTIDE P INNER JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
														WHERE ID_MODIFICATION IN ($ptmIdStrg) AND ID_ANALYSIS=$bestAnaID");
			}
			#$dbh->disconnect;
			# $numPepValues*=($numUsedPep/$bestNumLines); # ratio of usable values
			# $bestNumLines*=($numUsedPep/$bestNumLines); # ratio of usable values
			$bestNumLines*=($numUsedPep/$bestNumPep); # ratio of usable values for best Ana
			$bestNumLines+=(0.05*$bestNumLines*$numStates); # accounts for ~5% state-specific sites
		}
		# if ($quantifType eq 'DESIGN') {
		# 	$maxHours=($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None')? 360 : int(0.5+($numPepValues/25000)); $maxHours=2 if $maxHours < 2; $maxHours=48 if $maxHours > 48; # 1 M lines -> 40 h
		# 	#$maxMem=($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None')? int(1.5 + 5E-7 * $numPepValues * $numReplicates) : int(1.5 + 1E-6 * $numPepValues * $numReplicates); #$numStates
		# 	#$maxMem=($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None')? int(1.5 + 5E-6 * $bestNumLines * $numReplicates) : int(1.5 + 2.5E-6 * $bestNumLines * $numReplicates * $numStates); #$numStates
		# 	$maxMem=($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None')?			int(1.5 + 5E-6 * $bestNumLines * $numReplicates) :							# ABUNDANCE
		# 			($quantifParams{'DB'}{'ALGO_TYPE'}[0] eq 'PEP_INTENSITY')?	int(1.5 + 2.5E-6 * $bestNumLines * $numReplicates * $numStates) :			# Intensity-based ratios
		# 																		int(1.5 + 2.5E-6 * $bestNumLines * $numReplicates * $bestNumFractions);		# Ratio-based ratios (PEP_RATIO)
		# 	$numCPU=2;
		# 	$jobName="myProMS_protQuant_DDA_$quantifID";
		# }
		# els
		if ($quantifType =~ /DESIGN:MSstats/) {
			$maxHours=96;
			#$maxMem=int(1.5 + 5E-6 * $numPepValues * $numStates);
			$maxMem=int(1.5 + 5E-6 * $bestNumLines * $numReplicates * $numStates); # 1/4 compared to above
			#$numCPU=($cluster{'maxCPUs'} < 10)? $cluster{'maxCPUs'} : 10;  done in $cluster{runJob}
			$numCPU=8; # adjust also in startDesignQuantification R 'clusters' parameter
			$jobName="myProMS_protQuant_DIA.MSstats\_$quantifID";
		}
		else { # myProMS algo ($quantifType eq 'DESIGN' || $quantifType =~ /DESIGN:(TDA|DIA)/)
			#$maxHours=($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None')? 360 : 96;
			$maxHours=($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None')? 4+$numReplicates : 1+$numStates**2;
			#$maxMem=($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None')? int(1.5 + 5E-6 * $numPepValues * $numReplicates) : int(1.5 + 5E-6 * $numPepValues * $numStates);
			#$maxMem=($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None')? 			int(1.5 + 5E-6 * $bestNumLines * $numReplicates) : int(1.5 + 5E-6 * $bestNumLines * $numReplicates * $numStates);
			$maxMem=($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None')?			int(1.5 + 1E-5 * $bestNumLines * $numReplicates) :							# ABUNDANCE
					($quantifParams{'DB'}{'ALGO_TYPE'}[0] eq 'PEP_INTENSITY')?	int(3 + 5E-6 * $bestNumLines * $numReplicates * $numStates) :			# Intensity-based ratios
																				int(3 + 5E-6 * $bestNumLines * $numReplicates * $bestNumFractions);		# Ratio-based ratios (PEP_RATIO)
			#$numCPU=($cluster{'maxCPUs'} < 10)? $cluster{'maxCPUs'} : 10;  done in $cluster{runJob}
			$numCPU=2; # adjust also in startDesignQuantification R 'clusters' parameter
			my $dataType=($quantifType =~ /DESIGN:(TDA|DIA)/)? $1 : 'DDA';
			my $jobAlgo=($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None')? "ABUND.$dataType" : "RATIO.$dataType";
			$jobName="myProMS_protQuant_$jobAlgo\_$quantifID";
		}
	}
	elsif ($quantifType eq 'DESIGN:SSPA') {
		my ($anaIdStrg,$numAnaUsed)=(scalar keys %analysisList > 10)? (join(',',(keys %analysisList)[0..9]),10) : (join(',',keys %analysisList),scalar keys %analysisList); # check only 1st 10 for performance 
		my $numProtValues;
		if (scalar @trueQuantifiedModifIDs) {
			my $ptmIdStrg=join(',',@trueQuantifiedModifIDs);
			($numProtValues)=$dbh->selectrow_array("SELECT COUNT(*) FROM
													(SELECT P.ID_PEPTIDE FROM PEPTIDE P INNER JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
														WHERE ID_MODIFICATION IN ($ptmIdStrg) AND ID_ANALYSIS IN ($anaIdStrg) GROUP BY P.ID_PEPTIDE
													) RES");
			$numProtValues/=$numAnaUsed; # because was sum across all numAnaUsed
		}
		else {
			($numProtValues)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS IN ($anaIdStrg) AND VISIBILITY >= 1");
		}
		$numProtValues*=$numAllBioRep;
		$maxHours=48;
		#$maxMem=int(1.5 + $numProtValues*2E-5);
		$maxMem=int(1.5 + $numProtValues*5E-6);
		$maxMem=1 if $maxMem < 1;
		$numCPU=2;
		$jobName="myProMS_protQuant_SSPA_$quantifID";
	}
	elsif ($quantifType eq 'XICMCQ') {
		my $numAna=scalar keys %analysisList;
		$maxHours=24+12*$numAna;
		# my ($numAllPep)=$dbh->selectrow_array("SELECT COUNT(*) FROM PEPTIDE WHERE ID_ANALYSIS IN (".join(',',keys %analysisList).")");
		# $maxMem=int(1.5 + $numAllPep*5E-5); # 0.5 Go/10000 peptides
		# $maxMem=int(1.5 + 0.75*$numAna); # Low mem because new MCQ jobs will be launched from this master job
		# my $coeff=1.5-0.03*$numAna; $coeff=0.75 if $coeff < 0.75; # coeff decreases with number of analyses: from 1.5 to 0.75 (for 25+ ana)
		my $coeff=1.5-0.03*$numAna; $coeff=0.15 if $coeff < 0.15; # coeff decreases with number of analyses: from 1.5 to 0.15 (for 45+ ana)
		$maxMem=int(1.5 + $coeff*$numAna);
		$numCPU=2;
		$jobName="myProMS_pepQuant_$quantifType\_$quantifID";
	}
	elsif ($quantifType eq 'HYDRO') {
		my $numAna = scalar keys %analysisList;
		$maxHours = ($numAna > 10)? 4 : 2;
		$maxMem = ($numAna > 50)? 25 : int(1 + 0.5 * $numAna);
		$maxMem .= 'Gb';
		$numCPU = 2;
		$jobName = "myProMS_$quantifType\_$quantifID";
	}
	
	$dbh->disconnect;
	
	###>Running job on cluster<###
	#(my $cgiUnixPath=$0)=~s/\/[^\/]+$//;
	my $cgiUnixPath=`pwd`;
	$cgiUnixPath=~s/\/*\s*$//;
	# cd is required for script to find myproms .pm files!!!!
	my $commandString="export LC_ALL=\"C\"; cd $cgiUnixPath; $cluster{path}{perl}/perl $quantifScripts{$quantifType} $quantifID $jobDir 2>> $errorFilePath";
	my %jobParameters=(
		maxMem=>$maxMem.'Gb',
		numCPUs=>$numCPU,
		maxHours=>$maxHours,
		jobName=>$jobName,
		pbsRunDir=>$cgiUnixPath,
		commandBefore=>"mv $currentQuantifPath/$quantifItemID\_$jobDir\_wait.flag $currentQuantifPath/$quantifItemID\_$jobDir\_run.flag", # run flag file
        noWatch=>'1' # WARNING: in multi-job context, all jobs will be launched on cluster at once because "noWatch set to 1"
	);
	my ($pbsError,$pbsErrorFile,$jobClusterID,$pbsFile)=$cluster{'runJob'}->($quantifPath,$commandString,\%jobParameters);

    # Add cluster job id to current job in DB
	$dbh=&promsConfig::dbConnect('no_user');
    $dbh->do("UPDATE JOB_HISTORY SET ID_JOB_CLUSTER='C$jobClusterID' WHERE ID_JOB='$jobDir'");
    $dbh->commit;

	if ($pbsError) { # move PBS error message to job error file
		system "cat $pbsErrorFile >> $errorFilePath";
		$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID"); # Failed
		$dbh->commit;
	}
	elsif (-e "$promsPath{quantification}/project_$projectID/quanti_$quantifID") { # Move cluster files to quantif final directory (except if error in quantif)
		system "mv $quantifPath/$pbsFile $promsPath{quantification}/project_$projectID/quanti_$quantifID/.";
		system "mv $quantifPath/*.sh $promsPath{quantification}/project_$projectID/quanti_$quantifID/.";
	}

} # end of cluster launch
else { # Local launch (web server): No cluster OR XICCORR|SIN|EMPAI
    open(WAIT,">>$currentQuantifPath/$quantifID\_$jobDir\_wait.flag"); # flag file (created by selAna4Quantification for other quantif types than DESIGN* & XICMCQ)
    print WAIT "./$quantifScripts{$quantifType} $quantifID $jobDir $extraArgStrg 2>> $errorFilePath\n";
    close WAIT;

    # Add process PID to current job in DB
    $dbh->do("UPDATE JOB_HISTORY SET ID_JOB_CLUSTER='L$$' WHERE ID_JOB='$jobDir'");
    $dbh->commit;
	$dbh->disconnect;
    
	rename("$currentQuantifPath/$quantifItemID\_$jobDir\_wait.flag","$currentQuantifPath/$quantifItemID\_$jobDir\_run.flag"); # run flag file
	#my $analysisStrg=($quantifType =~ /SIN|EMPAI/)? $quantifItemID : '';
	#my $extraArgStrg=($quantifType eq 'SIN')? " $userID" : '';

	system "./$quantifScripts{$quantifType} $quantifID $jobDir $extraArgStrg 2>> $errorFilePath";
	if (-e $errorFilePath && -s $errorFilePath) {
		#system "sed -i '/Permanently/ d' $errorFile"; # remove warning from ssh key addition
		my $dbh=&promsConfig::dbConnect('no_user');
		#>Check if error already handled (by parent/child/other) process)
		my ($handled)=$dbh->selectrow_array("SELECT 1 FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID AND STATUS=-2");
		if ($handled) { # already handled: end quietly
			$dbh->disconnect;
			exit;
		}
		#>Handle error
		$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID"); # failed
		$dbh->commit;
		$dbh->disconnect;
		die "Aborting quantification due to errors.";
	}

} # end of local launch
#$dbh->disconnect;

#####>Revision history<####
# 1.8.18 [CHANGE] Lower cluster memory estimation for MassChroQ job (PP 02/06/21)
# 1.8.17 [FEATURE] Add computation of hydrophobicity through quantification workflow (VL 09/02/21)
# 1.8.16 [CHANGE] Moved myProMS quantif/SSPA code version management to corresponding scripts (PP 23/04/21)
# 1.8.15 [UPDATE] Changed memory computation for MassChroQ and myProMS DIA ratios jobs (PP 25/03/21)
# 1.8.14 [BUGFIX] Set min memory to 10Go for MassChroQ jobs (VS 19/03/21)
# 1.8.13 [BUGFIX] Minor fix to prevent duplicate TOP_N in QUANTIF_ANNOT field (PP 15/03/21)
# 1.8.12 [UPDATE] Minor changes in error handling (PP 19/02/21)
# 1.8.11 [UPDATE] Double cluster memory for Abundance quantification (PP 18/02/21)
# 1.8.10 [UPDATE] Cluster memory adjustment for PTM quantification (PP 16/02/21)
# 1.8.9 [BUGFIX] Quantification name is no longer recorded in JOB_HISTORY.FEATURES (PP 10/02/21)
# 1.8.8 [UPDATE] Unified cluster ressources estimation for DESIGN DDA and DESIGN DIA (PP 03/02/21)
# 1.8.7 [UPDATE] Handle resuming quantification & new algo version 3.5 for multi-level bias correction (PP 27/01/21)
# 1.8.6 [UPDATE] Change cluster memory for PEP_RATIO (PP 08/12/20)
# 1.8.5 [ENHANCEMENT] Cluster memory adjusted for myProMS ratio-quantif and PTM-quantif & compatibility with Free-residue quantification (PP 23/11/20)
# 1.8.4 [ENHANCEMENT] Extend categories to abundance quantifications for job monitoring (VL 19/11/20)
# 1.8.3 [UPDATE] Changed JOB_HISTORY.STATUS to JOB_HISTORY.JOB_STATUS (PP 28/08/20)
# 1.8.2 [UPDATE] More memory allowed for DIA.MSstats jobs (PP 28/08/20)
# 1.8.1 [UPDATE] Changed myProMS algo version to 3.4 (PP 14/08/20)
# 1.8.0 [CHANGE] Algo MSstats now launches runXICProtQuantification.pl instead of runSWATHProtQuantification.pl (12/08/20)
# 1.7.17 [ENHANCEMENT] runMassChroQ.pl is launched directly on cluster (PP 29/07/20)
# 1.7.16 [BUGFIX] Changed hard-coded MassChroQ version to "?" to be updated during run (PP 08/06/20)
# 1.7.15 [CHANGE] Handles NUM_TRUE_USED and ABUND_MEASURES parameters (PP 20/05/20)
# 1.7.14 [CHANGES] SSPA: Added software version number and 4x less cluster memory requirement (PP 29/03/20)
# 1.7.13 [CHANGES] Changed memory requirement for Abundance & DESIGN (PP 07/02/20)
# 1.7.12 [BUGFIX] Added forgotten culster memory computation for DESIGN:DIA (PP 29/01/20)
# 1.7.11 [ENHANCEMENT] Improved cluster memory calculation for SSPA & skip SETS recording for SSPA (PP 28/01/20)
# 1.7.10 [ENHANCEMENT] Improved support for SSPA (PP 22/01/20)
# 1.7.9 [FEATURE] Compatible with DIA by myProMS and Protein Abundance (PP 10/01/20)
# 1.7.8 [FEATURE] Compatible with PEPTIDES_NORM parameter (PP 03/12/19)
# 1.7.7 [ENHANCEMENT] Remove call to singularity image 1.1.19 for DIA quantif & [BUGFIX] Set undef SSPA quantif cluster parameters (PP 21/11/19)
# 1.7.6 [CHANGES] Add quantification name in job information (VS 19/11/19)
# 1.7.5 [BUGFIX] Fix DIA quantification computation issue thanks to the new singularity image (VS 18/11/19)
# 1.7.4 [ENHANCEMENT] Improved cluster memory calculation (PP 31/10/19)
# 1.7.3 Compatible with the new job monitoring system (VS 09/10/19)
# 1.7.2 Compatible with multi-modif quantifications & TDA by myProMS (PP 09/07/19)
# 1.7.1 Cluster memory calculation now corrected for the number of analyses included in peptide quantification (PP 03/05/19)
# 1.7.0 Compatible with XICCORR quantif type for Isobaric correction (PP 19/02/19)
# 1.6.6 New cluster memory calculation for DIA (PP 22/01/19)
# 1.6.5 Compatible with protein-level normalization of site quantification (PP 07/11/18)
# 1.6.4 Update swath files path (VS 08/11/18)
# 1.6.3 Removed useless 'use File::Copy' package (PP 12/10/18)
# 1.6.2 Also checks for jobs launched by other parent launchQuantifications.pl (PP 01/10/18)
# 1.6.1 Switch from _wait to _run flag now done on cluster to account for job queueing (PP 21/09/18)
# 1.6.0 Uses cluster for design quantifications (PP 17/09/18)
# 1.5.2 Auto-increment QUANTIFICATION update (GA 21/03/18)
# 1.5.1 'SWATH' changed to 'DIA' & MSstats version number now handled by runSWATHProtQuantification.pl (PP 26/12/17)
# 1.5.0 Adapted to new quantif algo v3 parameters, removed internal quantif code & relies on &promsConfig::getClusterInfo (PP 27/11/17)
# 1.4.7 Add MassChroQ version 2.2.1 in QUANTIF_ANNOT (GA 14/09/17)
# 1.4.6 Reverse minor modification in $sthLM (GA 03/05/17)
# 1.4.5 Minor modification in $sthLM query because some label modifications do not contain PSI-MS Name (GA 04/04/17)
# 1.4.4 Different SILAC vs XIC QUANTIFICATION_METHOD for MassChroQ extractions (GA 06/03/17)
# 1.4.3 Compatible with TMT labeling (PP 18/01/17)
# 1.4.2 Bug fix to allow all ratio pairs for multiSILAC (PP 04/11/16)
# 1.4.1 Minor change (PP 19/08/16)
# 1.4.0 Adding Swath (MSstats) & SSPA (PP 02/08/16)
# 1.3.7 Added matching peptides option for label-free quantification (PP 25/02/16)
# 1.3.6 Bug fix for internal quantif due to SOFTWARE=xxx tag in peptide quantif annot (PP 18/02/16)
# 1.3.5 Records label-free topN selection in QUANTIF_ANNOT (PP 12/10/15)
# 1.3.4 Number of parallele jobs depends on $promsPath{qsub} (PP 03/08/15)
# 1.3.3 Minor change in SOFTWARE tag for SILAC XIC with MassChroQ (PP 29/04/15)
# 1.3.2 Add Software in QUANTIF_ANNOT (GA 24/04/15)
# 1.3.1 Handles "Simple ratios" algorithm for SILAC triplex (PP 23/04/15)
# 1.3.0 Handles "Super ratios"=SILAC algorithm (PP 19/06/14)
# 1.2.7 Update for emPAI and SIN (GA 16/04/14)
# 1.2.6 Assumes no-label channel if incomplete annot data (PP 03/04/14)
# 1.2.5 Uses rmtree instead of remove_tree (PP 10/03/14)
# 1.2.4 Updated for improved state selection for internal quantif & (PP 26/03/14)
# 1.2.3 Add distinction between FREE and SILAC for MassChroQ extractions (29/01/14)
# 1.2.2 system command removal (13/11/13)
# 1.2.1 Now calls runXICProtQuantification.pl for design & internal quantifications (PP 09/10/13)
# 1.2.0 Minor changes in error management (PP 30/09/13)
# 1.1.9 MCQ -> XICMCQ & launched like DESIGN (PP 09/09/13)
# 1.1.8 Term 'TNPQ'->'DESIGN' for quantif type (PP 29/08/13)
# 1.1.7 Waiting for child finishing before launching the next child when max parallel-jobs number is reached for multiple-children quantifications (FY 28/08/13)
# 1.1.6 Improved TNPQ management (PP 11/07/13)
# 1.1.5 Handles targetPos in observation parameters (PP 03/07/13)
# 1.1.4 Commit before runMassChroQ call (PP 02/07/13)
# 1.1.3 GPL license (PP 26/06/13)
# 1.1.2 Change for quantification names + Multi XIC extractions (GA 07/11/12)
# 1.1.1 Change for MassChroQ call (GA 04/09/12)
# 1.1.0 Minor changes & bug fixes (PP 30/04/12)
# 1.0.9 Update $quantifAnnotStrg for emPAI (should be empty)<BR> + XIC PARENT_QUANTIFICATION no longer associated to DESIGN (GA 12/04/12)
# 1.0.8 Merge 1.0.7GA & 1.0.6PP (PP 02/04/12)
# 1.0.7 Modification of XIC quantification names (GA 20/02/2012)
# 1.0.6 Quantification recorded with status=-1 (PP 29/11/11)
# 1.0.5 Added selAna4quantification.cgi as launcher (PP 22/11/11)
# 1.0.4 Adding SILAC & iTRAQ internal quantifications (PP 26/08/11)
# 1.0.3 Replacing t3pq by XIC (Ext. Ion Chromat.) in the tests (GA ...)
# 1.0.2 Max numb parallel processes includes other jobDirs
