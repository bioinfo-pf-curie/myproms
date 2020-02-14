#!/usr/local/bin/perl -w

################################################################################
# launchQuantifications.pl       1.7.13                                        #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Launches multiple quantifications +/- paralellization                        #
# called by selAna4quantification.cgi                                          #
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

### myProMS software versions for quantification (based on R scripts & data preprocessing in runXIC...)
# 1: FC scripts (Ratio/TnPQ)
# 2: ML scripts
# 3: AS scripts before site-ratio correction by reference protein ratio (AnaDiff <= 4.0.13)
# 3.1: AS site-ratio correction & IB (AnaDiff >= 4.2.0)
# 3.2: Multi-modif quantif (format of modProtID is now <protID>-<modifRank1>#<pos1>.<pos2>&...)
# 3.3: TDA, DIA and Protein Abundance

$| = 1;

use strict;
use POSIX ":sys_wait_h"; # for WNOHANG
use File::Path qw(rmtree); # remove_tree
use promsConfig;
use promsMod;
use promsQuantif;

#exit;
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo('no_user');
my %cluster=&promsConfig::getClusterInfo;

###############################
####>Recovering parameters<####
###############################
my ($jobType,$userID,@extraArgs)=@ARGV;
my $MAX_PARALLEL_QUANTIFS=($jobType eq 'multi' && $cluster{'on'} && $extraArgs[0]=~/DESIGN/)? $cluster{'maxJobs'} : &promsConfig::getMaxParallelJobs; # multiCluster
#$MAX_PARALLEL_QUANTIFS=4; # tmp


                                                ########################
####################################################>Multi-job call<####################################################
                                                ########################
if ($jobType=~/multi/) {
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
my ($jobDir,$quantifType,$quantItemStrg)=@extraArgs; # $quantItemStrg defined only for SIN|EMPAI|XICCORR
my $quantifDir="$promsPath{tmp}/quantification/$jobDir";
my $currentQuantifDir="$promsPath{tmp}/quantification/current";
#my $quantifScript=($quantifType eq 'XIC')? 'runIonintensityQuantification.pl' : ($quantifType eq 'EMPAI')? 'runemPAIQuantification.pl' : ($quantifType eq 'SIN')? 'runSINQuantification.pl' : '';
my %quantifScripts=(
	'DESIGN'=>'runXICProtQuantification.pl',
	'DESIGN:TDA'=>'runXICProtQuantification.pl',
	'DESIGN:DIA'=>'runXICProtQuantification.pl',
	'DESIGN:MSstats'=>'runSWATHProtQuantification.pl',
	'DESIGN:SSPA'=>'runSSProtQuantification.pl',
	'EMPAI'=>'runemPAIQuantification.pl',
	'SIN'=>'runSINQuantification.pl',
	#'SILAC'=>'runXICProtQuantification.pl',
	#'ITRAQ'=>'runXICProtQuantification.pl',
	#'TMT'=>'runXICProtQuantification.pl',
	'XICMCQ'=>'runMassChroQ.pl', # different XIC method
	'XICCORR'=>'correctPepQuantification.pl'
);

#$quantItemStrg='' unless $quantItemStrg;
#my @quantItemList=split(':',$quantItemStrg);
#my %quantifIdList;
my ($quantifID,$quantifItemID,$quantifName); my $numAllBioRep=0; # $quantifItemID != $quantifID for SIN|emPAI only (=)
my $extraArgStrg='';
my (@pepQuantifFiles,%analysisList,$numStates,$numReplicates); # needed for cluster run
my %quantifParams;
my $labelType='FREE'; # default
if ($quantifType !~ /XICMCQ|XICCORR/) {
	my $dbh=&promsConfig::dbConnect('no_user');
	#>Fetching parameters
	%quantifParams=&promsQuantif::getQuantificationParameters("$quantifDir/quantif_info.txt");
	my $qRootName=quotemeta($quantifParams{'DB'}{'QUANTIF_NAME'}[0]);
	$labelType=$quantifParams{'DB'}{'LABEL'}[0]; # if $quantifType =~ /^(DESIGN|DESIGN:TDA|DESIGN:DIA)$/; # eq 'DESIGN';
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
			$quantifAnnotStrg.="::SOFTWARE=myProMS;3.3::ALGO_TYPE=$quantifParams{'DB'}{'ALGO_TYPE'}[0]"; # name;version::PEP_(RATIO/INTENSITY)
			$quantifAnnotStrg.=';changed' if $overwritten;
		}
		elsif ($quantifParams{'DB'}{'ALGO_TYPE'}[0] eq 'MSstats') {
			$quantifAnnotStrg.='::SOFTWARE=DIA/MSstats'; # Version number inserted by runSWATHProtQuantification.pl
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
	my @quantifiedModifIDs=($quantifParams{'DB'}{'PTM_SCOPE'})? @{$quantifParams{'DB'}{'PTM_SCOPE'}} : ();
	my $referenceQuantifID=($quantifParams{'DB'}{'REF_QUANTIF'})? $quantifParams{'DB'}{'REF_QUANTIF'}[0] : 0;
	foreach my $param (sort{&promsMod::sortSmart(lc($a),lc($b))} keys %{$quantifParams{'DB'}}) { # sortSmart required for CONDITION_<n>
		next if $param=~/QUANTIF_NAME|LABEL|NUM_REP_|ALGO_TYPE|PTM_SCOPE|MATCH_PEP/;
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
	}
	else { # DESIGN(:xxx)
		$quantifMethod=($quantifParams{'DB'}{'ALGO_TYPE'}[0] eq 'SSPA')? 'SSPA' : ($quantifParams{'DB'}{'ALGO_TYPE'}[0] eq 'TnPQ')? 'TNPQ' : ($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None')? 'PROT_ABUNDANCE' : 'PROT_RATIO_PEP';
		my @states=@{$quantifParams{'DB'}{'STATES'}};
		$numStates=scalar @states;
		#if ($quantifMethod eq 'SSPA') {
		#	##>Compute and order all possible sets of size 1 to #states-1 (<=> abundance sets). Position=TARGET_POS of quantif value
		#	$quantifAnnotStrg.="::SETS=";
		#	## First algo --->
		#	#my (%sets,%key2StatePos);
		#	#foreach my $i (0..$#states) {
		#	#	my %stateKeys;
		#	#	$sets{"#$i"}=$i;
		#	#	@{$key2StatePos{"#$i"}}=($i+1); # state pos represented by this key (to prevent key split)
		#	#	@{$stateKeys{1}}=("#$i"); # {size}
		#	#	foreach my $s (2..$#states) {
		#	#		foreach my $j ($i+1..$#states) {
		#	#			foreach my $prevKey (@{$stateKeys{$s-1}}) {
		#	#				next if $j <= $sets{$prevKey}; # prevents "State3State2"
		#	#				my $key=$prevKey."#$j";
		#	#				@{$key2StatePos{$key}}=(@{$key2StatePos{$prevKey}},$j+1);
		#	#				push @{$stateKeys{$s}},$key;
		#	#				$sets{$key}=$j;
		#	#			}
		#	#		}
		#	#	}
		#	#}
		#	#foreach my $key (sort{scalar @{$key2StatePos{$a}} <=> @{$key2StatePos{$b}} || $key2StatePos{$a}[0]<=>$key2StatePos{$b}[0]} keys %key2StatePos) {
		#	#	$quantifAnnotStrg.=join('+',@{$key2StatePos{$key}}).';';
		#	#}
		#	## Second (better?) algo -->
		#	my @sets;
		#	my $n = scalar(@states);
		#	my (@flags) = (0) x $n;
		#	$flags[0] = 1; # skip empty subset
		#	my $stateIdx;
		#	while (1) {
		#		my @subset = ();
		#		for ($stateIdx=0; $stateIdx<$n; $stateIdx++) {
		#			push @subset, $stateIdx+1 if $flags[$stateIdx] == 1; # record state positions
		#		}
		#		last if scalar @subset == $n; # exclude complete subset (size=n) and exit while loop
		#		push @sets,\@subset;
		#		for ($stateIdx=0; $stateIdx<$n; $stateIdx++){
		#			if($flags[$stateIdx] == 0) {
		#				$flags[$stateIdx] = 1;
		#				last;
		#			}
		#			$flags[$stateIdx] = 0;
		#		}
		#	}
		#	foreach my $refSubset (sort{scalar @{$a} <=> scalar @{$b} || $a->[0]<=>$b->[0]} @sets) { # size then 1st state in subset (should be enough)
		#		$quantifAnnotStrg.= join('+',@{$refSubset}).';';
		#	}
		#	print scalar @sets,"\n";
		#
		#	$quantifAnnotStrg=~s/;\Z//; # remove trailing ';'
		#}
		#els
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
		$numReplicates=0;
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
		$quantifName=$quantifParams{'DB'}{'QUANTIF_NAME'}[0].'_'.$anaName; # Auto_quantifType_anaName
		my ($lastName)=$dbh->selectrow_array("SELECT NAME FROM QUANTIFICATION Q,ANA_QUANTIFICATION A WHERE Q.ID_QUANTIFICATION=A.ID_QUANTIFICATION AND A.ID_ANALYSIS=$analysisID AND FOCUS='$focus' AND NAME LIKE '$qRootName%' ORDER BY Q.ID_QUANTIFICATION DESC LIMIT 1");
		if ($lastName) {
			$quantifName.=($lastName=~/#(\d+)\Z/)? " #".($1+1) : " #2";
		}
		$sthInsQ->execute($quantifName,$quantifAnnotStrg);
		$quantifID=$dbh->last_insert_id(undef, undef,'QUANTIFICATION','ID_QUANTIFICATION');
		$sthInsAnaQ->execute($quantifID,$analysisID);
		
		###> Update quantif_info.txt in case needed later (PP 30/10/18)
		open (INFO,">>$quantifDir/quantif_info.txt");
		print INFO "QUANTIFICATIONS:\n$quantifID";
		close INFO;
	}
	else {# DESIGN updates
		$quantifName=$quantifParams{'DB'}{'QUANTIF_NAME'}[0];
		my $numQuantifiedModifs=scalar @quantifiedModifIDs;
		#my $singleQuantifModID=($numQuantifiedModifs==1)? abs($quantifiedModifIDs[0]) : undef;
		#$sthInsQ->execute($singleQuantifModID,$quantifName,$quantifAnnotStrg);# Add the new created
		$sthInsQ->execute($quantifName,$quantifAnnotStrg);# Add the new created
		$quantifID=$dbh->last_insert_id(undef, undef, 'QUANTIFICATION', 'ID_QUANTIFICATION');
		#if ($numQuantifiedModifs >= 2) { #}
		if ($numQuantifiedModifs) {
			my $modifRank=0;
			foreach my $qModID (@quantifiedModifIDs) { # already sorted in quantif_info.txt
				$modifRank++;
				my $usedModifRank=($numQuantifiedModifs==1)? undef : $modifRank;
				$sthInsModQ->execute($quantifID,abs($qModID),$usedModifRank);
			}
		}
		$quantifItemID=$quantifID;
		my $designID=$quantifParams{'DB'}{'ID_DESIGN'}[0];
		$dbh->do("UPDATE QUANTIFICATION SET ID_DESIGN=$designID WHERE ID_QUANTIFICATION=$quantifID");

		#my $sthUpQ=$dbh->prepare("UPDATE QUANTIFICATION SET ID_DESIGN=$designID WHERE ID_QUANTIFICATION=?"); --> PARENT_QUANTIFICATION like XIC should not be associated to a DESIGN <=> 12/04/12 Modification
		my $sthInsEQ=$dbh->prepare("INSERT INTO EXPCONDITION_QUANTIF (ID_QUANTIFICATION,ID_EXPCONDITION) VALUES ($quantifID,?)");
		#my $sthObs=$dbh->prepare("SELECT ID_ANALYSIS FROM OBSERVATION WHERE ID_OBSERVATION=?");

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
		#$sthObs->finish;

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
		open (INFO,">>$quantifDir/quantif_info.txt");
		print INFO "QUANTIFICATIONS:\n$quantifID";
		close INFO;

		#>Generating job list
		#@quantItemList=($quantifID);
		unlink "$currentQuantifDir/$jobDir\_request.flag" if -e "$currentQuantifDir/$jobDir\_request.flag"; # in case multi-job launch
		open(WAIT,">$currentQuantifDir/$quantifID\_$jobDir\_wait.flag"); # flag file (created by selAna4Quantification for other quantif types than DESIGN* & XICMCQ)
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
	####>Getting the parameters from quantif_info.txt<####
	open (INFO,"$quantifDir/quantif_info.txt");
	my $section='';
	my (%params,%sampIDs);
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

	my %quantifParams=&promsQuantif::getQuantificationParameters("$quantifDir/quantif_info.txt");
	my $qRootName=quotemeta($quantifParams{'DB'}{'QUANTIF_NAME'}[0]); # Name given by user
	my $anaIDS=join(',', keys (%{$params{'MZXML'}}) );
	$quantifName=$quantifParams{'DB'}{'QUANTIF_NAME'}[0];
	my ($lastName)=$dbh->selectrow_array("SELECT NAME FROM QUANTIFICATION Q,ANA_QUANTIFICATION A WHERE Q.ID_QUANTIFICATION=A.ID_QUANTIFICATION AND A.ID_ANALYSIS IN ($anaIDS) AND FOCUS='peptide' AND NAME LIKE '$qRootName%' ORDER BY Q.ID_QUANTIFICATION DESC LIMIT 1");
	if ($lastName) {
		$quantifName.=($lastName=~/#(\d+)\Z/)? " #".($1+1) : " #2";
	}

	my $quantifAnnotStrg=($params{'CHANNELS'})? "LABEL=SILAC::SOFTWARE=MCQ;2.2.1::$params{'CHANNELS'}":"LABEL=FREE::SOFTWARE=MCQ;2.2.1";
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
	open (INFO,">>$quantifDir/quantif_info.txt");
	print INFO "QUANTIFICATIONS:\n$quantifID";
	close INFO;

#@quantItemList=($quantifID);
	open(WAIT,">$currentQuantifDir/$quantifID\_$jobDir\_wait.flag"); # flag file (created by selAna4Quantification for other quantif types than DESIGN* & XICMCQ)
	print WAIT '#';
	close WAIT;
#$quantifIdList{$quantifID}=$quantifID;
}
elsif ($quantifType eq 'XICCORR') {
	$quantifItemID=$quantItemStrg; # "anaID.parentQuantifID"
	$quantifID=(split(/\./,$quantItemStrg))[1];
	
	####>Getting the parameters from quantif_info.txt<####
	my $productID=0;
	my $keepOld=0;
	open (INFO,"$quantifDir/quantif_info.txt");
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
#exit; # DEBUG

##########################################
####>Launching quantification process<#### only 1 child now (PP 03/08/18)
##########################################
my $dbh=&promsConfig::dbConnect('no_user');
my $projectID=&promsMod::getProjectID($dbh,$quantifID,'QUANTIFICATION');

$quantifName = '' if(!$quantifName);
my $quantiDBInfo = "ID_QUANTIFICATION=$quantifID;TYPE=$quantifType;QUANTIFICATION_NAME=$quantifName";
$quantiDBInfo .= ";ID_ANALYSIS=$quantItemStrg" if($quantItemStrg);
my $logFilePath = ($quantItemStrg) ? "$quantifDir/status\_$quantItemStrg.out" : "$quantifDir/status\_$quantifID.out";
my $errorFilePath = ($quantItemStrg) ? "$currentQuantifDir/$quantItemStrg\_$jobDir\_error.txt" : "$currentQuantifDir/$quantifID\_$jobDir\_error.txt";
$dbh->do("INSERT INTO JOB_HISTORY (ID_JOB, ID_USER, ID_PROJECT, TYPE, STATUS, FEATURES, SRC_PATH, LOG_PATH, ERROR_PATH, STARTED) VALUES('$jobDir', '$userID', $projectID, 'Quantification [$quantifType]', 'Queued', '$quantiDBInfo', '$quantifDir', '$logFilePath', '$errorFilePath', NOW())");
$dbh->commit;

if ($cluster{'on'} && $quantifType =~ /^DESIGN/) { #'DESIGN' or 'DESIGN:MSstats' or 'DESIGN:SSPA' or 'DESIGN:TDA'
	###>Estimating ressources required<###
	my ($maxHours,$maxMem,$numCPU,$jobName)=(24,'1Gb',1,"myProMS_protQuant_$quantifType\_$quantifID"); # default

	if ($quantifType eq 'DESIGN') {
		my $numPepValues=0;
		my %usedQuantifFile;
		foreach my $resQuantInfo (@pepQuantifFiles) {
			my ($parQuantifID,$anaID,$targetPos,$numAna)=@{$resQuantInfo};
			my $fileKey=$parQuantifID.':'.$targetPos;
			unless ($usedQuantifFile{$fileKey}) {
				my $pepQuantifFile=($targetPos)? "peptide_quantification_$targetPos.txt" : 'peptide_quantification.txt';
				my $numLines=`wc -l $promsPath{quantification}/project_$projectID/quanti_$parQuantifID/$pepQuantifFile`;
				#chomp($numLines);
				$numLines=~s/^(\d+).*/$1/;
				$usedQuantifFile{$fileKey}=int($numLines/$numAna);
			}
			$numPepValues+=$usedQuantifFile{$fileKey};
		}
		$maxHours=($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None')? 360 : int(0.5+($numPepValues/25000)); $maxHours=2 if $maxHours < 2; $maxHours=48 if $maxHours > 48; # 1 M lines -> 40 h
		$maxMem=($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None')? int(1.5 + 5E-7 * $numPepValues * $numReplicates) : int(1.5 + 1E-6 * $numPepValues * $numReplicates); #$numStates
##		$maxMem=$cluster{'maxMem'} if $maxMem > $cluster{'maxMem'};  done in $cluster{runJob}
		$maxMem.='Gb';
		$numCPU=2;
		$jobName="myProMS_protQuant_DDA_$quantifID";
	}
	elsif ($quantifType =~ /DESIGN:(MSstats|TDA|DIA)/) {
		my $algoType=(split(':',$quantifType))[1];
		my $numPepValues=0;
		my $sthNumPep=$dbh->prepare("SELECT COUNT(*) FROM PEPTIDE WHERE ID_ANALYSIS=?");
		foreach my $resQuantInfo (@pepQuantifFiles) {
			my ($parQuantifID,$anaID,$targetPos,$numAna)=@{$resQuantInfo};
			if ($algoType eq 'MSstats') {
				my $numLines=`wc -l $promsPath{quantification}/project_$projectID/quanti_$parQuantifID/swath_ana_$anaID.txt`;
				#chomp($numLines);
				$numLines=~s/^(\d+).*/$1/;
				$numPepValues+=$numLines; #int($numLines/$numAna);
			}
			else { # myProMS algo
				$sthNumPep->execute($anaID);
				my ($numPep)=$sthNumPep->fetchrow_array;
				$numPepValues=$numPep if $numPep > $numPepValues; # keep max
			}
		}
		$sthNumPep->finish;
		$maxHours=($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None')? 360 : 96;
		$maxMem=($quantifParams{'DB'}{'RATIO_TYPE'}[0] eq 'None')? int(1.5 + 5E-6 * $numPepValues * $numReplicates) : int(1.5 + 5E-6 * $numPepValues * $numStates);
		#$maxMem=int(2.5 + 2E-5 * $numPepValues * $numStates);
##		$maxMem=($maxMem < 10)? 10 : ($maxMem > $cluster{'maxMem'}) ? $cluster{'maxMem'} : $maxMem;  done in $cluster{runJob}
		$maxMem.='Gb';
##		$numCPU=($cluster{'maxCPUs'} < 10)? $cluster{'maxCPUs'} : 10;  done in $cluster{runJob}
		$numCPU=($algoType eq 'MSstats')? 8 : 2;
		my $jobAlgo=($algoType eq 'MSstats')? 'DIA.MSstats' : $algoType;
		$jobName="myProMS_protQuant_$jobAlgo\_$quantifID";
	}
	if ($quantifType eq 'DESIGN:SSPA') {
		my ($numProtValues)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS IN (".join(',',keys %analysisList).") AND VISIBILITY >= 1");
		$numProtValues*=$numAllBioRep;
		$maxHours=48;
		$maxMem=int(1.5 + $numProtValues*2E-5);
		$maxMem=1 if $maxMem < 1;
		$maxMem.='Gb';
		$numCPU=2;
		$jobName="myProMS_protQuant_SSPA_$quantifID";
	}
	
	###>Running job on cluster<###
##	my $pbsFile="$quantifDir/PBS.txt"; # !!!! Moved to parent directory !!!
##	my $pbsErrorFile="$quantifDir/PBSerror.txt"; # !!!! Moved to parent directory !!!
	###my $clusterQuantifDir="$quantifDir/quanti_$quantifID";
	###mkdir $clusterQuantifDir;
	#(my $cgiUnixDir=$0)=~s/\/[^\/]+$//;
	my $cgiUnixDir=`pwd`;
	$cgiUnixDir=~s/\/*\s*$//;
	# cd is required for script to find myproms .pm files!!!!
	my $commandString="export LC_ALL=\"C\"; cd $cgiUnixDir; $cluster{path}{perl}/perl $quantifScripts{$quantifType} $quantifID $jobDir 2> $currentQuantifDir/$quantifID\_$jobDir\_error.txt";
	my %jobParameters=(
		maxMem=>$maxMem,
		numCPUs=>$numCPU,
		maxHours=>$maxHours,
		jobName=>$jobName,
		pbsRunDir=>$cgiUnixDir,
		commandBefore=>"mv $currentQuantifDir/$quantifItemID\_$jobDir\_wait.flag $currentQuantifDir/$quantifItemID\_$jobDir\_run.flag", # run flag file
        noWatch=>'1',
	);
	my ($pbsError,$pbsErrorFile, $jobClusterID)=$cluster{'runJob'}->($quantifDir,$commandString,\%jobParameters);
    
    # Add cluster job id to current job in DB
    $dbh->do("UPDATE JOB_HISTORY SET ID_JOB_CLUSTER='C$jobClusterID' WHERE ID_JOB='$jobDir'");
    $dbh->commit;
    
	if ($pbsError) { # move PBS error message to job error file
		system "cat $pbsErrorFile >> $currentQuantifDir/$quantifID\_$jobDir\_error.txt";
		$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID"); # Failed
		$dbh->commit;
	}
	elsif (-e "$promsPath{quantification}/project_$projectID/quanti_$quantifID") { # Move cluster files to quantif final directory (except if error in quantif)
		system "mv $quantifDir/*.txt $promsPath{quantification}/project_$projectID/quanti_$quantifID/.";
		system "mv $quantifDir/*.sh $promsPath{quantification}/project_$projectID/quanti_$quantifID/.";
	}

} # end of cluster launch
else { # Local launch (web server): No cluster OR XICMCQ|SIN|EMPAI
    open(WAIT,">>$currentQuantifDir/$quantifID\_$jobDir\_wait.flag"); # flag file (created by selAna4Quantification for other quantif types than DESIGN* & XICMCQ)
    print WAIT "./$quantifScripts{$quantifType} $quantifID $jobDir $extraArgStrg 2> $currentQuantifDir/$quantifItemID\_$jobDir\_error.txt\n";
    close WAIT;
    
    # Add process PID to current job in DB
    $dbh->do("UPDATE JOB_HISTORY SET ID_JOB_CLUSTER='L$$' WHERE ID_JOB='$jobDir'");
    $dbh->commit;
    
	rename("$currentQuantifDir/$quantifItemID\_$jobDir\_wait.flag","$currentQuantifDir/$quantifItemID\_$jobDir\_run.flag"); # run flag file
	#my $analysisStrg=($quantifType =~ /SIN|EMPAI/)? $quantifItemID : '';
	#my $extraArgStrg=($quantifType eq 'SIN')? " $userID" : '';
    
	system "./$quantifScripts{$quantifType} $quantifID $jobDir $extraArgStrg 2> $currentQuantifDir/$quantifItemID\_$jobDir\_error.txt";
    
} # end of local launch
$dbh->disconnect;

#####>Revision history<####
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
