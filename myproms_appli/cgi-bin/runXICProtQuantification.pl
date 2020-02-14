#!/usr/local/bin/perl -w

################################################################################
# runXICProtQuantification.pl       2.14.0                                     #
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
use strict;
#use XML::SAX;
use MIME::Base64;
use POSIX qw(strftime :sys_wait_h);
use File::Copy qw(copy move);
use File::Copy::Recursive qw(dirmove dircopy);
use File::Path qw(rmtree); # remove_tree
use promsConfig;
use promsMod;
use promsQuantif;

#exit;
###############################
####>Recovering parameters<####
###############################
my ($quantifID,$quantifDate,$runMode)=@ARGV; # only $quantifID & $quantifDate defined if Design-quantif
$runMode='MASTER' unless $runMode; # MASTER, ABUND_NORM, ABUND_RSCRIPT, ABUND_SPLIT:<numJobs>, ABUND_RUN:<numJobs>, ABUND_JOB:<jobNum>, REF, NORM

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo('no_user');
my %cluster=&promsConfig::getClusterInfo; #('debian') # default is 'centos'
my $pathR=($cluster{'on'})? $cluster{'path'}{'R'} : $promsPath{'R'};
my $pathPython=($cluster{'on'})? $cluster{'path'}{'python'} : $promsPath{'python'};
#my $minProtNumber=30; # used to desactivate FDR option
#my $currentQuantifDir="$promsPath{tmp}/quantification/current";
my $quantifDir="$promsPath{tmp}/quantification/$quantifDate";
my $fileStat="$quantifDir/status_$quantifID.out"; # : "$quantifDir/status_$quantItemID.out";
my $runDir="$quantifDir/quanti_$quantifID";
my $dataDir="$runDir/data";
my $resultDir="$runDir/results";
my $graphDir="$resultDir/graph";

my %quantifParameters=&promsQuantif::getQuantificationParameters("$quantifDir/quantif_info.txt");
my $algoType=$quantifParameters{'DB'}{'ALGO_TYPE'}[0]; # <=> {R}{design} for v3+ only! [PEP_RATIO|PEP_INTENSITY|TDA|DIA]
my $ratioType=$quantifParameters{'DB'}{'RATIO_TYPE'}[0]; # SuperRatio or SimpleRatio or None (Abundance);
my $numSteps=($ratioType=~/S\w+Ratio/)? 3 : 4;


#>>>>>>>>> PROCESSING PROTEIN ABUNDANCE METHOD <<<<<<<<<<<<<
if ($runMode eq 'MASTER' && $ratioType eq 'None') { # 1st call for Abundance
	
	my $dbh=&promsConfig::dbConnect('no_user');
	$dbh->do("UPDATE QUANTIFICATION SET STATUS=0 WHERE ID_QUANTIFICATION=$quantifID");
	$dbh->commit;
	
	#>Normalizing peptide XICs
#if (-e "$promsPath{tmp}/quantification/abundance/data") {
#	my $dbh=&promsConfig::dbConnect('no_user');
#	$dbh->do("UPDATE QUANTIFICATION SET STATUS=0 WHERE ID_QUANTIFICATION=$quantifID");
#	$dbh->commit;
#	$dbh->disconnect;
#	open(FILESTAT,">>$fileStat");
#	print FILESTAT "Retrieving data from previous quantif\n";
#	close FILESTAT;
#	mkdir $runDir unless -e $runDir; # already created if runXICProtQuantification.pl launched on cluster
#	#mkdir $dataDir;
#	#mkdir $resultDir;
#	#mkdir $graphDir;
#	system "mv $promsPath{tmp}/quantification/abundance/data $runDir/.";
#	system "mv $promsPath{tmp}/quantification/abundance/results $runDir/.";
#	system "mv $promsPath{tmp}/quantification/abundance/AnalysisDiff* $runDir/.";
#}
#else {
	if (!-e $dataDir || !-e "$dataDir/table.txt") { # in case using data from previous incomplete quantif
		system "./runXICProtQuantification.pl $quantifID $quantifDate ABUND_NORM"; # 2>> $currentQuantifDir/$quantifID\_$quantifDate\_error.txt; # Default mode. Stops after R normalization results
	}
	
	#>Launch R normalization
	unless (-e "$resultDir/resultsPep.txt") { # in case using data from previous incomplete quantif
		system "./runXICProtQuantification.pl $quantifID $quantifDate ABUND_RSCRIPT"; # 2>> $currentQuantifDir/$quantifID\_$quantifDate\_error.txt; # Default mode. Stops after R normalization results
		&checkRsuccess;
	}
#}	
	#>Splitting resultsPep.txt for parallel jobs
	$dbh=&promsConfig::dbConnect('no_user'); # reconnect
	my ($quantifAnnot)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID");
	$quantifAnnot=~s/::/§/g;
	my ($stateStrg)=$quantifAnnot=~/STATES=([^§]+)/;
	my $numAllRep=scalar(split(/[;.&]/,$stateStrg));
	my $numLines=`wc -l $resultDir/resultsPep.txt`;
	$numLines=~s/^(\d+).*/$1/;
	my $numJobLn=int(1.5 + $numLines / 2.5E5);
	my $numJobSt=$numAllRep;
	my $numJobs=($numJobLn >= $numJobSt)? $numJobLn : $numJobSt;	
	$numJobs=5 if $numJobs < 5;
	
	unless (-e "$resultDir/ABUND$numJobs/resultsPep.txt") {
		system "./runXICProtQuantification.pl $quantifID $quantifDate ABUND_SPLIT:$numJobs"; # 2>> $currentQuantifDir/$quantifID\_$quantifDate\_error.txt; # Compute a Abundance & LFQ from R results
	}
	
	#>Launching parallel abundance computation jobs
	system "./runXICProtQuantification.pl $quantifID $quantifDate ABUND_RUN:$numJobs"; # 2>> $currentQuantifDir/$quantifID\_$quantifDate\_error.txt; # Compute a Abundance & LFQ from R results
	
	#>Parsing & importing results in DB
	&importProteinAbundance($numJobs);
	exit;
}
elsif ($runMode eq 'ABUND_RSCRIPT') {
	&launchRscript;
	exit;
}
elsif ($runMode=~/ABUND_SPLIT:(\d+)/) { # Split resultsPep.txt file for parallel LFQ jobs. Done in a different process to release memory after reading result file
	&splitResultsPepFile($1);
	exit;
}
elsif ($runMode=~/ABUND_RUN:(\d+)/) {
	&runProteinAbundance($1);
	exit;
}
elsif ($runMode=~/ABUND_JOB:(\d+)/) {
	&jobProteinAbundance($1);
	exit;
}

######>>>>>>>>>> $runMode: MASTER/REF/NORM (diffAna) or ABUND_NORM (normalization only) >>>>>>>>

my $topN=($quantifParameters{'DB'}{'TOP_N'})? $quantifParameters{'DB'}{'TOP_N'}[0] : 0; # number peptide used for label-free (except old algo 'Ratio')
my $matchingPep=($quantifParameters{'DB'}{'MATCH_PEP'})? $quantifParameters{'DB'}{'MATCH_PEP'}[0] : 0; # label-free SimpleRatio: peptides must be the same across conditions
my $referenceMode=($runMode =~/^(REF|NORM)$/)? $runMode : ''; # for table(Ref/Norm).txt
my $referenceStrg=($referenceMode eq 'REF')? ' reference' : ($referenceMode eq 'NORM')? ' normalization' : '';
#open (DEBUG,">$promsPath{tmp}/quantification/debug.txt") if !$referenceMode; # DEBUG!!!!
my $dbh=&promsConfig::dbConnect('no_user');

if (!$referenceMode) {
	if ($ratioType=~/S\w+Ratio/) {
		$dbh->do("UPDATE QUANTIFICATION SET STATUS=0 WHERE ID_QUANTIFICATION=$quantifID");
		$dbh->commit;
	}
	
	open(FILESTAT,">$fileStat") || die "Error while opening $fileStat";
	print FILESTAT "Started ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
	print FILESTAT "1/$numSteps Generating data files\n";
	close FILESTAT;
}
#elsif ($referenceMode eq 'REF') {
#	open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
#	print FILESTAT "1/$numSteps Generating data files (0/4 Writing$referenceStrg protein file)\n";
#	close FILESTAT;
#}

my $projectID=&promsMod::getProjectID($dbh,$quantifID,'QUANTIFICATION');

#my ($quantiAnnotation,$quantifiedModifID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_MODIFICATION FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID");
my ($quantiAnnotation,$quantifiedModifID,$quantifModStrg)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,Q.ID_MODIFICATION,GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ',')
																					 FROM QUANTIFICATION Q
																					 LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
																					 WHERE Q.ID_QUANTIFICATION=$quantifID
																					 GROUP BY Q.ID_QUANTIFICATION");
my $isModifQuantif=0;
my (%modifQuantifRank,%rankQuantifModif,@modifQuantifs);
my $matchQuantifModStrg='';
if ($quantifiedModifID || $quantifModStrg) {
	$isModifQuantif=1;
	$quantifModStrg=$quantifiedModifID if $quantifiedModifID;
	my $modifRank=0;
	foreach my $modID (split(',',$quantifModStrg)) {
		$modifRank++;
		$modifQuantifRank{$modID}=$modifRank;
		$rankQuantifModif{$modifRank}=$modID;
		push @modifQuantifs,$modID; # ordered list of modifs
	}
	$quantifiedModifID=$modifQuantifs[0] if $modifRank==1; # back compatibility in case QUANTIFICATION.ID_QUANTIFICATION no longer used
	$matchQuantifModStrg=join('|',@modifQuantifs); # Needed later for regexp on pep varmod
}
my ($labeling)=($quantiAnnotation=~/LABEL=([^:]+)/);
$labeling=uc($labeling);

my %skipFromRefModIDs;
if ($referenceMode) { # Switching to non-modif quantif
	%skipFromRefModIDs=%modifQuantifRank;
	$quantifiedModifID=0;
	$isModifQuantif=0;
}
my $quantifiedModifRemoved=0; # to handle site for which targeted PTM was removed by experimental protocole
if ($isModifQuantif && $quantifParameters{'DB'}{'CREATE_PTM'}) {
	$quantifiedModifRemoved=1;
}

####> R part : statistical analysis of the proteins values
####> 1st: Create the directories for R analysis
#my $runDir=($call eq 'DESIGN')? "$quantifDir/quanti_$quantifID" : "$quantifDir/quanti_$quantItemID"."_$quantifID";

unless ($referenceMode) {
	#make_path($dataDir,$graphDir,{verbose=>0,mode=>0775}); # $runDir,$resultDir will be created automatically
	mkdir $runDir unless -e $runDir; # already created if runXICProtQuantification.pl launched on cluster
	mkdir $dataDir;
	mkdir $resultDir;
	mkdir $graphDir;

	####>Protein-level correction for quantif sites (tableRef.txt)
	if ($quantifParameters{'DB'}{'INTRA_PROT_REF'}) {
		if ($quantifParameters{'DB'}{'INTRA_PROT_REF'}[0] eq '-1') { # current dataset
			system "./runXICProtQuantification.pl $quantifID $quantifDate REF"; # 2>> $currentQuantifDir/$quantifID\_$quantifDate\_error.txt; # referenceMode set to REF to create tableRef.txt
			$quantifParameters{'R'}{'normalization.ref.test'}[0]=$quantifParameters{'DB'}{'NORMALIZATION_METHOD'}[0]; # 'ref' part of ref.test
		}
		else {
			$quantifParameters{'R'}{'normalization.ref.test'}[0]=&generateReferenceProtFromQuantif; # 'ref' part of ref.test
		}
		$quantifParameters{'R'}{'normalization.ref.test'}[0].='.none' if $quantifParameters{'R'}{'normalization.ref.test'}[0] !~ /\./; # quantile # 'ref' part of ref.test
		$quantifParameters{'R'}{'normalization.ref.test'}[0].='.median.scale'; # 'test' part of ref.test (hard-coded for now!!) <== ************************
	}
	####>Peptide set experimental normalization for quantif sites (tableNorm.txt)
	if ($quantifParameters{'DB'}{'PEPTIDES_NORM'}) {
		system "./runXICProtQuantification.pl $quantifID $quantifDate NORM"; # 2>> $currentQuantifDir/$quantifID\_$quantifDate\_error.txt; # referenceMode set to NORM to create tableNorm.txt
	}
}

####>Main data file
my $dataFile=($referenceMode eq 'REF')? "$dataDir/tableRef.txt" : ($referenceMode eq 'NORM')? "$dataDir/tableNorm.txt" : "$dataDir/table.txt";
open(DATA,">$dataFile"); # Name of data table
#if ($ratioType eq 'Ratio') { # Old algos (FC)
#	print DATA "Protein_ID\tPep_Seq\tVar_Mod\tCharge\tPeptide_IDs\tProtein_Name\tProtein_Validity";
#}
#else { # Algos v2+
	print DATA "Protein_ID\tPeptide\tSample\tTreatment\tReplicate\tTechnical_Replicate\tExperiment\tQuantif_Set\tPeptide_IDs\tProtein_Name\tProtein_Validity\tValue" unless $referenceMode eq 'NORM';
#}
my %statesInFile;

################################
####>Get states composition<####
################################
my @quantiOrder;# In order to print the table in a good way for R
my %observationInfo; # for (Super/Simple)Ratio only
my (%poolGroup2Cond,%ana2Obs,%cond2Obs,%obs2Cond,%labelModifs,%anaLabelModifs,%xicEngine); #,%obsID2Cond,%ana2Cond
my (%ana2Experiment,%anaConds,%anaState1RatioPos); #<----- NEW EXPERIMENT assignment 31/10/17 -----
my $sthPQN=$dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?"); # peptide quantifs

my @ratioTarget;
my %replicateTargetPos; my $numRepTgtPos=0;
my $numConditions=0;
my $poolObservations=0;
#my %experimentCode; # SILAC paired experiments
my @conditionList;
#my %designCheck; # Checks if experimental design matches algo type (label w. new Algo only) TODO: Implement actual test
my $sthALM=$dbh->prepare("SELECT M.ID_MODIFICATION FROM ANALYSIS_MODIFICATION AM,MODIFICATION M WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND AM.ID_ANALYSIS=? AND M.IS_LABEL=1");
#if ($call eq 'DESIGN') {
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
				#if ($ratioType eq 'Ratio') {
				#	push @{$quantifParameters{'R'}{'grp'}},$numBioRep; # Number of replicates of the condition
				#	push @{$quantifParameters{'R'}{'name.grp'}},"State$numConditions"; # $condID
				#}
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
##>REPLICATE DEFINITION: BIOREP or TECHREP BASED on RESIDUAL_VAR!!!
$numRepTgtPos++ if ($techRepCount==1 || $quantifParameters{'DB'}{'RESIDUAL_VAR'}[0] eq 'technical'); # count either bioRep or techRep
$replicateTargetPos{"State$numConditions.rep$bioRepCount.techRep$techRepCount"}=$numRepTgtPos;
						###$multiTechRep=1 if $techRepCount>1; # records if exists multiple tech rep
						push @quantiOrder,$techReplicate;
						$obs2Cond{$techReplicate} = $condID;
#print DEBUG "\$techReplicate='$techReplicate' \$condID=$condID \$numConditions=$numConditions\n";
						#if ($ratioType eq 'Ratio') {
						#	print DATA "\tState$numConditions"; # $condID
						#	print DATA ".rep$allRepCount" if $multiRep;
						#}
						#>Required to detect +/- infinite ratios (except %labelModifs)
						$poolGroup2Cond{$techReplicate}=$condID;
						my %replicAna;
						foreach my $fraction (split(/\+/,$techReplicate)) {
							my @obsData=split(/:/,$fraction); # obsID,parQuantifID,anaID,targetPos
							my ($obsID,$parQuantifID,$anaID,$targetPos)=@obsData;
							unless ($xicEngine{$parQuantifID}) {
								$sthPQN->execute($parQuantifID);
								my ($parQuantiAnnot)=$sthPQN->fetchrow_array;
								$xicEngine{$parQuantifID}=($parQuantiAnnot=~/SOFTWARE=MCQ|EXTRACTION_ALGO=/)? 'MCQ' : ($parQuantiAnnot=~/SOFTWARE=MQ/)? 'MQ' : 'PD';
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
if ($ratioType=~/S\w+Ratio|None/) {
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

###>Computing number of targetPos used (needed later for recording of num pep per replicate)<###
my $usedTargetPos=scalar @ratioTarget; # ($quantifParameters{'DB'}{'SINGLE_REF'})? $numConditions-1 : ($numConditions * ($numConditions-1)) / 2;
$usedTargetPos+=$numConditions if $algoType=~/PEP_INTENSITY|TDA|DIA/; # num states (MEAN_INTENSITY)
foreach my $techReplicate (keys %replicateTargetPos) {
	$replicateTargetPos{$techReplicate}+=$usedTargetPos; # shift to free target pos
}

$sthPQN->finish;
$sthALM->finish;

#close DEBUG;
#$dbh->disconnect; die "Test is complete!";

###>Handling Modification Quantif phosphoRS/MaxQuant threshold & PTM position ambiguity settings
my ($ptmProbSoftCode,$ptmProbThreshold,$ambiguousModifPos,%quantifResMatchStrg,%qQuantifResMatchStrg)=('',0,0,'',''); # default ($qQuantifResMatchStrg PP 2017/06/22)
my (%targetableRes,$matchProtCterm); # for recreated PTMs only
if ($isModifQuantif) {
	if ($quantifParameters{'DB'}{'CREATE_PTM'}) { # assumes single-modif quantif!!!!!!!
		#$quantifResMatchStrg=join('',@{$quantifParameters{'DB'}{'CREATE_PTM'}});
		#$quantifResMatchStrg=~s/,.//g; # remove contexts
		foreach my $createModStrg (@{$quantifParameters{'DB'}{'CREATE_PTM'}}) {
			my ($modifRank,@targetRes)=split(/[:\.]/,$createModStrg);
			my $modID=$rankQuantifModif{$modifRank};
			foreach my $resCode (@targetRes) {
				my ($res,$context)=split(',',$resCode);
				push @{$targetableRes{$modID}},[$res,$context];
				$quantifResMatchStrg{$modID}.=$res;
				$matchProtCterm=1 if ($context && $context eq '+'); # for at least 1 modifID
			}
		}
	}
	else {
		if ($quantifParameters{'DB'}{'PTM_POS'}[0]=~/(\w+):(\d+)/) { # PRS Phospho or MQ any PTM
			$ptmProbSoftCode=$1; # PRS or MQ
			$ptmProbThreshold=$2; $ptmProbThreshold/=100 if $ptmProbSoftCode eq 'MQ'; # 0-1
			$ambiguousModifPos=1 if $quantifParameters{'DB'}{'PTM_POS'}[1] eq 'ambiguous';
		}
		elsif ($quantifParameters{'DB'}{'PTM_POS'}[0] eq 'ambiguous') { # other PTMs
			$ambiguousModifPos=1;
		}
	
		#>Find list of targeted residues: needed even if !$ambiguousModifPos
		my $sthAMR=$dbh->prepare("SELECT ID_MODIFICATION,SPECIFICITY FROM ANALYSIS_MODIFICATION WHERE ID_MODIFICATION IN ($quantifModStrg) AND ID_ANALYSIS=? AND MODIF_TYPE='V'");
		my %resList;
		foreach my $anaID (keys %ana2Obs) {
			$sthAMR->execute($anaID);
			my ($modID,$specifStrg)=$sthAMR->fetchrow_array;
			foreach my $res (split(',',$specifStrg)) {$resList{$modID}{$res}=1;}
		}
		$sthAMR->finish;
		my $sthMS=$dbh->prepare('SELECT SPECIFICITY FROM MODIFICATION WHERE ID_MODIFICATION=?');
		foreach my $modID (@modifQuantifs) {
			if ($resList{$modID}) {
				$quantifResMatchStrg{$modID}=join('',sort keys %{$resList{$modID}}); # eg 'STY'
			}
			else { # Fallback in case not specified during search(es)
				$sthMS->execute($modID);
				($quantifResMatchStrg{$modID})=$sthMS->fetchrow_array;
				$quantifResMatchStrg{$modID}=~s/,//g;
			}
		}
	}
	foreach my $modID (keys %quantifResMatchStrg) {
		$qQuantifResMatchStrg{$modID}=quotemeta($quantifResMatchStrg{$modID}); # +* -> \+\* (PP 2017/06/22)
	}
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
if ($quantifParameters{DB}{BIAS_CORRECTION}[1] && !$referenceMode) { # done only once in normal mode
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
my $manualPepSelection=0;
my (%manualPeptides,%manualProteins);
###>Protein-level fold-change correction
if ($referenceMode eq 'REF' && $quantifParameters{'DB'}{'PEPTIDES_REF'}[0] eq 'manual') { # if 'current', settings above are kept
	$manualPepSelection=1;
	($pepSpecificity,$pepMissedCleavage,$ptmFilter,$pepChargeState,$pepSource)=('all',1,1,'all','all'); # No filter at all
	foreach my $i (1..$#{$quantifParameters{'DB'}{'PEPTIDES_REF'}}) {
		my ($protID,@pepIDs)=split(/[:,]/,$quantifParameters{'DB'}{'PEPTIDES_REF'}[$i]);
		foreach my $pepID (@pepIDs) {$manualPeptides{$pepID}=1;}
#print DEBUG "MANUAL PEP=@pepIDs\n";
		$manualProteins{$protID}=1;
	}
}
###>Peptide set for experimental normalization (incompatible with Protein-level fold-change correction)
elsif ($referenceMode eq 'NORM' && $quantifParameters{'DB'}{'PEPTIDES_NORM'}[0] eq 'manual') {
	$manualPepSelection=1;
	$normProtUsage='use'; # needed to print normalization file
	foreach my $i (1..$#{$quantifParameters{'DB'}{'PEPTIDES_NORM'}}) {
		my ($protID,@pepIDs)=split(/[:,]/,$quantifParameters{'DB'}{'PEPTIDES_NORM'}[$i]);
		foreach my $pepID (@pepIDs) {$manualPeptides{$pepID}=1;}
#print DEBUG "MANUAL PEP=@pepIDs\n";
		$manualProteins{$protID}=1;
		$forNormalization{$protID}=1;
	}
}

#my $pepQuantifMethodCode=($labeling eq 'FREE')? 'XIC' : uc($labeling);
#my $pepQuantifCode=($labeling eq 'FREE')? 'XIC_AREA' : ($labeling eq 'SILAC')? 'ISO_AREA' : 'REP_AREA'; # assumes ITRAQ
#my %pepParamCode2ID;
#my $sthPepQP=$dbh->prepare("SELECT QP.CODE,ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER QP,QUANTIFICATION_METHOD QM WHERE QP.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND QM.CODE='$pepQuantifMethodCode'");
#$sthPepQP->execute;
#while (my ($code,$paramID)=$sthPepQP->fetchrow_array) {
#	$pepParamCode2ID{$code}=$paramID;
#}
#$sthPepQP->finish;
#my %peptideAreaParamID;
#my $sthPepQP=($labeling eq 'TMT')? $dbh->prepare("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION Q,QUANTIFICATION_PARAMETER QP WHERE Q.ID_QUANTIFICATION_METHOD=QP.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=? AND CODE='REP_INTENSITY'")
#								 : $dbh->prepare("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION Q,QUANTIFICATION_PARAMETER QP WHERE Q.ID_QUANTIFICATION_METHOD=QP.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=? AND CODE LIKE '%_AREA'"); # CRITICAL: '%_AREA' !!!!!!!!!!!!!!!!!!!!!!!!!!!
my %usedParamIDs;
my $sthPepQP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION Q,QUANTIFICATION_PARAMETER QP WHERE Q.ID_QUANTIFICATION_METHOD=QP.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=? AND (CODE LIKE '%_INTENSITY' OR CODE LIKE '%_AREA')");
foreach my $pepQuantifID (keys %xicEngine) {
	$sthPepQP->execute($pepQuantifID);
	#($peptideAreaParamID{$pepQuantifID})=$sthPepQP->fetchrow_array;
	while (my ($paramID)=$sthPepQP->fetchrow_array) {
		$usedParamIDs{$pepQuantifID}{$paramID}=1;
	}
}
$sthPepQP->finish;


##>Normalization with all peptides or only those used for quantification<##
my $normalizeWithAll=0; # 0/1 (0<=>old version) <--- Will become a parameter

#my $filterQueryStrg=($normalizeWithAll)? ',PPA.IS_SPECIFIC,AP.PEP_SPECIFICITY,MISS_CUT' : '';

my $ptmProbQueryStrg=($ptmProbSoftCode eq 'MQ')? "GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,'%',COALESCE(PM.REF_POS_STRING,'') ORDER BY PM.ID_MODIFICATION SEPARATOR '&')" : "'-'";
my $quantiQuery=qq
|SELECT GROUP_CONCAT(DISTINCT PPA.ID_PROTEIN),P.ID_PEPTIDE,ABS(PEP_BEG),PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),$ptmProbQueryStrg,CHARGE,DATA,P.SCORE,PPA.IS_SPECIFIC,AP.PEP_SPECIFICITY,MISS_CUT
	FROM PEPTIDE P
	LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
	INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
	INNER JOIN ANALYSIS_PROTEIN AP ON PPA.ID_PROTEIN=AP.ID_PROTEIN
	WHERE P.ID_ANALYSIS=AP.ID_ANALYSIS AND AP.ID_ANALYSIS=? AND ABS(AP.VISIBILITY) >= 1|; # VIS>=1 to allow multi-prot/matchGroup in manual case MG edition (PP 19/11/18) # Critical to restrict DISTINCT on ID_PROTEIN for execution time !!! $pepParamCode2ID{$pepQuantifCode}
if (!$normalizeWithAll && !$manualPepSelection) { # filter data
	$quantiQuery.=($pepSpecificity eq 'unique')? ' AND PPA.IS_SPECIFIC=1' : ($pepSpecificity eq 'unique_shared')? ' AND AP.PEP_SPECIFICITY=100' : ''; # Filter at DB level
	$quantiQuery.=' AND MISS_CUT=0' if $pepMissedCleavage==0; # 0: exclude only missed-cut peptides, 1: allow, -1: also exclude overlapping fully cut
}
if ($manualPepSelection) { # Do not apply user-defined peptide selection. Filter applied later on pepID
	$quantiQuery.=($referenceMode eq 'REF')? ' AND PPA.ID_PROTEIN IN ('.(join(',',keys %manualProteins)).')' : ' AND P.ID_PEPTIDE IN ('.(join(',',keys %manualPeptides)).')'; # NORM
}
$quantiQuery.=' GROUP BY PPA.ID_PROTEIN,P.ID_PEPTIDE ORDER BY ABS(AP.VISIBILITY) DESC,AP.NUM_PEP DESC,ABS(PEP_BEG) ASC';


#print DEBUG ">QUERY---------------\n$quantiQuery\n----------------------\n" if !$referenceMode;
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
foreach my $modID (@modifQuantifs) {
	$allowedPtmID{$modID}=1; # to allow modif quantification
}
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
my %notUsed4Quantif; # for future $normalizeWithAll management (Not yet used in table.txt)
my %pepID2QuantifSetSILAC; # PP 25/10/17: used to better compute %usedPeptideSets for SILAC when qSets are summed

##----- NEW EXPERIMENT assignment ----->
## RULES: 1 techRep <-> 1 experiment => 1 analysis <-> 1 experiment => 1 pepID <-> 1 experiment
## Should handle:
##    -fractions
##    -bioRep across different Experiments (LABELED with 3+ channels)
##    -SuperRatio: state1 found in multiple Experiments
##my (%analysis2Experiment,%techRep2Experiment); # %peptideID2Experiment,
##my $currentExperiment='A';

if (!$referenceMode || $referenceMode eq 'REF') {
	my $entityStrg=($algoType=~/^(TDA|DIA)$/)? 'fragment' : 'peptide';
	open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
	print FILESTAT "1/$numSteps Generating data files (1/4 Retrieving$referenceStrg $entityStrg intensity data)";
	close FILESTAT;
}
my %pepQuantifValues;
if ($algoType=~/^(TDA|DIA)$/) { # MS2 XIC
	my %fragQuantifs;
	foreach my $observationSet (@quantiOrder) {
		foreach my $fraction (split(/\+/,$observationSet)) {
			my ($obsID,$ms2Quantif,$anaID,$targetPos)=split(/:/,$fraction);
			$fragQuantifs{$ms2Quantif}{$anaID}=1;
		}
	}
	QUANTI:foreach my $ms2Quantif (keys %fragQuantifs) {
		foreach my $anaID (keys %{$fragQuantifs{$ms2Quantif}}) {
			unless (-e "$promsPath{quantification}/project_$projectID/quanti_$ms2Quantif/swath_ana_$anaID.txt") {
				warn "Missing fragment quantification file (#$ms2Quantif, Ana #$anaID)"; # warn to allow clean handling of exception by &checkForErrors
				last QUANTI;
			}
			open (FRAG_QUANTIF,"$promsPath{quantification}/project_$projectID/quanti_$ms2Quantif/swath_ana_$anaID.txt") || die $!;
#print DEBUG "XIC in $anaID:\n" if !$referenceMode;
			while (<FRAG_QUANTIF>) {
				next if /^PEP/; # $.=1
				s/\s\Z//g; # chomp is not enough. More hidden Windows character
				my ($pepID,$fragMZ,$fragCharge,$fragType,$fragRes,$fragArea,$fragRT)=split(/!/,$_);
				next if $fragArea=~/^(N\/A|None|0)$/;
				#next if (!$peptideList{$pepID});
				$pepQuantifValues{"$pepID:0"}{"$fragType$fragRes#$fragCharge"}=$fragArea;
#print DEBUG "$pepID\t$fragArea\n" if !$referenceMode;
			}
			close FRAG_QUANTIF;
		}
	}
}
else { # MS1 XIC
	my %peptideQuantifs;
	foreach my $observationSet (@quantiOrder) {
		foreach my $fraction (split(/\+/,$observationSet)) {
			my ($obsID,$xicQuantif,$anaID,$targetPos)=split(/:/,$fraction);
			$targetPos=0 unless $targetPos; # label-free
			$peptideQuantifs{$xicQuantif}{$targetPos}=1;
		}
	}
	QUANTI:foreach my $xicQuantif (keys %peptideQuantifs) {
		my $signalParamID;
		foreach my $targetPos (keys %{$peptideQuantifs{$xicQuantif}}) {
			my $pepQuantifFile=($labeling eq 'FREE')? 'peptide_quantification.txt' : "peptide_quantification_$targetPos.txt";
			unless (-e "$promsPath{quantification}/project_$projectID/quanti_$xicQuantif/$pepQuantifFile") {
				warn "Missing peptide quantification file (#$xicQuantif)"; # warn to allow clean handling of exception by &checkForErrors
				last QUANTI;
			}
			open (PEP_QUANTIF,"$promsPath{quantification}/project_$projectID/quanti_$xicQuantif/$pepQuantifFile") || die $!;
			while(<PEP_QUANTIF>) {
				next if $.==1;
				chomp;
				my ($paramID,$pepID,$quantifValue)=split(/\t/,$_);
				#next if $paramID != $peptideAreaParamID{$xicQuantif};
				next if $quantifValue <=0; # Isobaric bug
				unless ($signalParamID) {
					next unless $usedParamIDs{$xicQuantif}{$paramID};
					$signalParamID=$paramID;
				}
				next if $paramID != $signalParamID;
				$pepQuantifValues{"$pepID:$targetPos"}{PEP}=$quantifValue;
			}
			close PEP_QUANTIF;
		}
	}
}
&checkForErrors($dbh);

open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
print FILESTAT "1/$numSteps Generating data files (2/4 Retrieving$referenceStrg identification data)";
close FILESTAT;
my %recreatedVarMods; # to prevent re-creating same varMod multiple times
my %proteinLength; # Only if sites to be recreated at protein C-term
my (%missedCutPeptides,%beg2pepSeq); # only if $pepMissedCleavage=-1 (exclude missCut AND overlapping cut peptides)
my $sthProtL=$dbh->prepare('SELECT PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=?');
foreach my $observationSet (@quantiOrder) {
	#my $labelObs=scalar keys %{$labelModifs{$observationSet}}; # observation belongs (or not) to a label channel
	my $anaInObsStrg=''; # list of anaIDs in observationSet
	if ($labeling eq 'FREE') { # && $ratioType ne 'Ratio'
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

		my $qSetStrg=($xicEngine{$xicQuantif}=~/^(PD|MQ)$/)? 'QSET' : "MCQSET_$xicQuantif";
		my %peptidesUsed; # only if non-proteotypic peptides are used => makes sure the same protein is selected (analysis-specific to respect Match group decision!!!)
		my $labelModStrg=join('|',keys %{$anaLabelModifs{$anaID}}) if $labeling eq 'SILAC'; # needed for SILAC
		my %mcqXICs; # SILAC only
		$sthGetPepData->execute($anaID);
#print DEBUG "RUN_ANA $anaID ($targetPos)\n" if !$referenceMode;
		while (my ($protID,$pepID,$pepBeg,$pepSeq,$varModStrg,$varProbStrg,$charge,$pepData,$score,$specif,$bestSpecif,$misscut)=$sthGetPepData->fetchrow_array) { # executed multiple time for labeled quantif

			#>Protein filtering
			if ($protSelectionType) {
				if ($protSelectionType eq 'exclude') {next if $selectExcludeProteins{$protID};}
				else {next unless $selectExcludeProteins{$protID};} # restrict
			}
			
#$varModStrg='' unless $varModStrg;
#print DEBUG "$pepID,$pepBeg,$pepSeq,$varModStrg\n" if !$referenceMode;
			next if ($manualPepSelection && (!$manualProteins{$protID} || !$manualPeptides{$pepID}));
#print DEBUG "OK FILTER_1\n" if !$referenceMode;

if ($pepMissedCleavage==-1 && $misscut) {
	my $pepLength=length($pepSeq);
	$missedCutPeptides{$protID}{$pepBeg}=$pepLength if (!$missedCutPeptides{$protID} || !$missedCutPeptides{$protID}{$pepBeg} || $missedCutPeptides{$protID}{$pepBeg} < $pepLength);
	if ($normalizeWithAll) {$notUsed4Quantif{$pepID}=1;}
	else {next;}
}
											   
			next unless $pepQuantifValues{"$pepID:$targetPos"}; # no quantif for pepID
			my $quantifValue=$pepQuantifValues{"$pepID:$targetPos"}; # ref to last PEP/FRAGMENTS dimension
#print DEBUG "OK FILTER_2\n" if !$referenceMode;

			#>Normalization scope
			$forNormalization{$protID}=1 if ($normProtUsage eq 'not' && !$notForNormalization{$protID});
	
			### NEW (25/01/19): Handle quantif PTM removed from samples => unmodifed res ARE the sites!
			## PTM recreated in silico
			$varModStrg='' unless $varModStrg;
			if ($quantifParameters{'DB'}{'CREATE_PTM'}) {
				if ($matchProtCterm && !defined $proteinLength{$protID}) {
					$sthProtL->executed($protID);
					($proteinLength{$protID})=$sthProtL->fetchrow_array;
				}
				$varModStrg=&recreateModifSites($pepSeq,$varModStrg,$pepBeg,$proteinLength{$protID},\%targetableRes,\%recreatedVarMods);
			}

########!!!!!! TEST 18/08/15
#next if $preInfRatios{$protID};
###################################
#next if $excludedProteins{$protID}; # beta
			if (($isModifQuantif && (!$varModStrg || $varModStrg !~ /(^|&)($matchQuantifModStrg):/)) || ($varModStrg && $varModStrg=~/\d:(\.|&|$)/) ) { # missing mod pos data!
				if ($normalizeWithAll) {$notUsed4Quantif{$pepID}=1;}
				else {next;} # skip peptides not containing quantified modif
			}
#print DEBUG "OK FILTER_3\n" if !$referenceMode;
			if ($pepSpecificity ne 'unique' && !$manualPepSelection) {  # && $bestSpecif < 100 not proteotypic peptide && in manual mode (TDA ref): allow duplicate
				next if ($peptidesUsed{$pepSeq} && $peptidesUsed{$pepSeq} != $protID); # Use peptide ONLY once even if shared by multiple proteins ($pepSpecificity= unique_shared or all) (Required by ML algos)
				$peptidesUsed{$pepSeq}=$protID;
			}
$proteinInConds{$protID}{$poolGroup2Cond{$observationSet}}=1; # protein "exists" in condition (detect before peptide filtering to prevent missleading info!!!)
#$visInReplicate{$protID}{$observationSet}=1; # protein is visible in replicate (VISIBILITY > 0 in query) needed for table.txt only
#$visInAnalysis{$protID}{$anaID}=1;
#$bestVisibility{$protID}=$vis if (!defined($bestVisibility{$protID}) || $bestVisibility{$protID} < $vis); # Keep the Best-Visibility of the protein (only needed for FREE).
			if ($normalizeWithAll) {
				$notUsed4Quantif{$pepID}=1 if (($pepSpecificity eq 'unique' && !$specif) || ($pepSpecificity eq 'unique_shared' && $bestSpecif<100)); # proteo-typic filter
				$notUsed4Quantif{$pepID}=1 if ($pepMissedCleavage <= 0 && (!defined $misscut || $misscut > 0)); # misscleavage filter (also skip if no info)
			}

			$pepData='' unless $pepData;
			$score=0 unless $score;
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
			if ($varModStrg && ($labeling ne 'FREE' || $ptmAllowed != 1 || $referenceMode)) { # There is filtering on PTM
				my %hasPTMs;
				my $hasRefMod=0; # for referenceMode only
				foreach my $vMod (split('&',$varModStrg)) {
					my ($modID)=($vMod=~/^(\d+)/);
					next if $labelModifs{$obsID}{$modID}; # skip labeling modifs
					if ($ptmAllowed != 1 && !$allowedPtmID{$modID}) {
						$hasPTMs{$modID}=1;
						last;
					}
					$varModCode.='&'.$vMod;
					if ($skipFromRefModIDs{$modID}) {
						$hasRefMod=1;
						last;
					}
				}
				if (scalar %hasPTMs) { # this modification is not allowed
					$excludedSeq{$pepSeq}=1 if $ptmAllowed <= -1; # -1,-2: unmodified peptide not allowed if exists a modified version
					next;
				}
				if ($hasRefMod) { # referenceMode & target PTM matched: modified & unmodified peptides not allowed
					$excludedSeq{$pepSeq}=1;
					next;
				}
			}
			elsif ($varModStrg) {
				$varModCode='&'.$varModStrg;
			}
#print DEBUG "OK FILTER_4\n" if !$referenceMode;
			# Non-SILAC phospho quantif
			if ($ptmProbSoftCode && $labeling ne 'SILAC' && !$notUsed4Quantif{$pepID}) { # !$ptmProbSoftCode if 'CREATE_PTM'
				my $ptmPosProb;
				if ($ptmProbSoftCode eq 'PRS') { 
					($ptmPosProb)=($pepData=~/PRS=\d;([^;]+);/);
				}
				else { # assume MQ
					$ptmPosProb=&extractPtmPosProbability(\%modifQuantifRank,$varModCode,$varProbStrg);
				}
				if (!defined $ptmPosProb || $ptmPosProb < $ptmProbThreshold) {
					if ($normalizeWithAll) {$notUsed4Quantif{$pepID}=1;}
					else {next;} 
				}
			}
#print DEBUG "OK FILTER_5\n" if !$referenceMode;
			#$intensity{$pepID}=$quantifValue;
			#unless ($dataSrc) {$dataSrc=($labeling eq 'FREE')? '-' : $anaID;}
			#$dataSrc=($labeling eq 'FREE')? '-' : ($dataSrc)? $dataSrc : '#'.$anaID;
#$dataSrc=$dataSrc || '#'.$anaID; # define dataSrc for label-free also
			my $dataSrc;
# mPP----->
			#if ($labeling eq 'FREE' && $ratioType eq 'Ratio') {$dataSrc='-';}
			##if ($labeling eq 'FREE') {$dataSrc='-';} # <---- mPP
# <---- mPP
			#else {
				$dataSrc='#'.$anaID; # TODO: if LABEL-FREE: check if still compatible with fraction merge ($anaInObsStrg in Qset) PP 27/11/17
				$dataSrc.='.'.$1 if $pepData=~/SOURCE_RANK=([^#]+)/;
			#}
			my $qSet;
			if ($labeling eq 'FREE') { # No peptide pairing   || $algoType eq 'PEP_INTENSITY'
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
					$qSet="A$anaID"."_Q0"; # "P$pepID"."_Q0" in algo v2!
					#$qSet=($ratioType eq 'Ratio')? 'A0_Q0' : "A$anaID"."_Q0"; # "P$pepID"."_Q0" in algo v2!
					#$intensity{"$pepSeq:$varModCode:$charge"}{$dataSrc}{$qSet}{$observationSet}[0]=$quantifValue;
				}
				#$intensity{"$pepSeq:$varModCode:$charge"}{$dataSrc}{$qSet}{$observationSet}[0]+=$quantifValue; # pool just in case multiple occurences
				foreach my $feature (keys %{$quantifValue}) { # PEP of Fragments
					$intensity{"$pepSeq:$varModCode:$charge"}{$dataSrc}{$qSet}{$observationSet}{$feature}[0]+=$quantifValue->{$feature};
				}
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
					push @{$tempIntensity{"$pepSeq:#:$varModCode:#:$charge"}{$dataSrc}{$fraction}{$qSet}},$quantifValue->{PEP}; # $intensity computed later
					$dataSource2ObsSet{$dataSrc}{$observationSet}=$fraction;
					$qSet2Protein{$qSet}=$protID;

					$pepID2QuantifSetSILAC{$pepID}=$qSet;
				}
				else {
					push @{$intensity{"$pepSeq:$varModCode:$charge"}{$dataSrc}{$qSet}{$observationSet}{PEP}},$quantifValue->{PEP};
				}
			}

			#my $usedObs=($ratioType eq 'Ratio')? 'ALL_OBS' : $observationSet;
			#push @{$proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet}{$usedObs}},$pepID;
push @{$proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet}{$observationSet}},$pepID;
			$peptideBeg{$protID}{$pepSeq}=$pepBeg; # for modif quantification only
			
if ($pepMissedCleavage==-1 && !$misscut) {push @{$beg2pepSeq{$protID}{"$pepBeg:$pepSeq"}},$pepID;} # no need to include missed cleaved here
	
##>Records duplicate XIC scans if MCQ (eg multiple PD XICs reunited in a single one by MCQ => keep only 1 qSet)
if ($xicEngine{$xicQuantif} eq 'MCQ' && $labeling eq 'SILAC') { # $labeling !~ /ITRAQ|TMT/
	$mcqXICs{"$pepSeq:$charge:".$quantifValue->{PEP}}{$qSet}=1; # duplicate if more than 1 qSet with the same pepSeq,charge (ideally mass) have exactly the same XIC value
} # assumes EXACTLY same XIC is computed by MCQ

			##>Record best varMod of qSet in case incoherence (assumes same pepSeq!)
			#if ($labeling ne 'FREE') { #} qSET has no meaning in label-free ***Also for FREE because of possible duplicate scan for same pepID with MassChroQ***
			if ($labeling eq 'SILAC') {
				if ($ptmProbSoftCode) { # Position confidence filtered earlier for non-SILAC. !$ptmProbSoftCode if 'CREATE_PTM'
					my $ptmPosProb;
					if ($ptmProbSoftCode eq 'PRS') {
						($ptmPosProb)=($pepData=~/PRS=\d;([^;]+);/);
					}
					else { # assume MQ
						$ptmPosProb=&extractPtmPosProbability(\%modifQuantifRank,$varModCode,$varProbStrg);
					}
					$ptmPosProb=0 unless $ptmPosProb;
					@{$qSetBestVarMod{$qSet}}=($ptmPosProb,$varModCode) if (!$qSetBestVarMod{$qSet} || $qSetBestVarMod{$qSet}[0] < $ptmPosProb); # record best phosphoRS/MaxQuant prob & vmod (both peptides in qSet will be filtered based on best prob)
				} 
				else {
					@{$qSetBestVarMod{$qSet}}=($score,$varModCode) if (!$qSetBestVarMod{$qSet} || $qSetBestVarMod{$qSet}[0] < $score); # record best
				}
			}
			#if ($ambiguousModifPos) {
			#	$retentionTime{$pepID}=($elutionStrg && $elutionStrg=~/et([\d\.]+)/)? $1 : ($elutionStrg && $elutionStrg=~/^([\d\.]+)/)? $1 : 0;
			#}
#print DEBUG "OK\n" if !$referenceMode;
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
#%pepQuantifValues=(); # hoping to free memory

open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
print FILESTAT "1/$numSteps Generating data files (3/4 Processing$referenceStrg dataset)";
close FILESTAT;
	
##>Filter fully cut peptides included in missed-cut (if $pepMissedCleavage==-1). Match is CROSS-ANALYSIS!
if ($pepMissedCleavage==-1) {
	foreach my $protID (keys %missedCutPeptides) {
		next unless $beg2pepSeq{$protID}; # no fully cut peptides for protID
		foreach my $missBeg (keys %{$missedCutPeptides{$protID}}) {
			my $missEnd=$missBeg + $missedCutPeptides{$protID}{$missBeg} - 1;
			foreach my $begSeq (keys %{$beg2pepSeq{$protID}}) {
				my ($pepBeg,$pepSeq)=split(':',$begSeq);
				if ($pepBeg >= $missBeg && $pepBeg <= $missEnd) { # match! beg overlaps with missed cut seq
					if ($normalizeWithAll) {
						foreach my $pepID (@{$beg2pepSeq{$protID}{$begSeq}}) {$notUsed4Quantif{$pepID}=1;}
					}
					else {$excludedSeq{$pepSeq}=1;}
					delete $beg2pepSeq{$protID}{$begSeq};
					delete $beg2pepSeq{$protID} unless scalar keys %{$beg2pepSeq{$protID}};
				}
			}
		}
	}
}
undef %missedCutPeptides;
undef %beg2pepSeq;

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
					next if ($ptmProbSoftCode && $qSetBestVarMod{$qSet}[0] < $ptmProbThreshold && !$ambiguousModifPos);
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

						$intensity{"$pepSeq:$usedVmod:$charge"}{$dataSrc}{$qSetMerged}{$observationSet}{PEP}[0]=$sumQuantifValues{$observationSet}{$usedVmod};

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
						@{$qSetBestVarMod{$qSetMerged}}=@{$qSetBestVarMod{$qSet}} if (!$qSetBestVarMod{$qSetMerged} || $qSetBestVarMod{$qSetMerged}[0] > $qSetBestVarMod{$qSet}[0]); # record the worse (PP 04/01/19)
						# DO NOT DELETE $qSetBestVarMod{$qSet} -> Can be reused in another pepIon context!
					}
				}
				##>No qSet merge (qSetMerged = qSet) => simple update of %intensity
				else {
					my $qSet=$qSetMerged;
					my $protID=$qSet2Protein{$qSet};
					foreach my $observationSet (keys %{$dataSource2ObsSet{$dataSrc}}) {
						next if (!$sumQuantifValues{$observationSet} || !$sumQuantifValues{$observationSet}{$usedVmod}); # no data for this channel
						$intensity{"$pepSeq:$usedVmod:$charge"}{$dataSrc}{$qSet}{$observationSet}{PEP}[0]=$sumQuantifValues{$observationSet}{$usedVmod}; # if ($sumQuantifValues{$observationSet} && $sumQuantifValues{$observationSet}{$usedVmod});
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

if (scalar keys %proteins == 0) { # No data at all!!!
	my $dataType=($referenceMode eq 'REF')? ' protein reference' : ($referenceMode eq 'NORM')? ' normalizing peptide' : '';
	warn "ERROR: No$dataType data available for quantification!";
	&checkForErrors($dbh);
}
my $sthAlias=$dbh->prepare("SELECT ID_PROTEIN,ALIAS FROM PROTEIN WHERE ID_PROTEIN IN (".join(',',keys %proteins).")");
$sthAlias->execute;
while (my ($protID,$alias)=$sthAlias->fetchrow_array) {
	$proteinAlias{$protID}=$alias;
}
$sthAlias->finish;

#print DEBUG scalar keys %proteins," proteins\n";
#close DEBUG if !$referenceMode;
#die "END OF TEST";

######################################
####> Printing normalization file<####
######################################
if ($normProtUsage) {
	unless (scalar keys %forNormalization) {
		warn "ERROR: No proteins found for normalization!\n";
		&checkForErrors($dbh);
	}
	#if ($ratioType eq 'Ratio') {
	#	foreach my $protID (keys %forNormalization) {
	#		push @{$quantifParameters{'R'}{'prot.ref'}},$protID;
	#	}
	#}
	#else {
		open(NORM,">$dataDir/normProtein.txt");
		print NORM "proteinId\n"; # header
		foreach my $protID (sort{$a<=>$b} keys %forNormalization) {
			print NORM "$protID\n";
		}
		close NORM;
	#}
}

##############################
####> Printing data table<####
##############################
open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
print FILESTAT "1/$numSteps Generating data files (4/4 Printing$referenceStrg dataset to file)";
close FILESTAT;
my $protValidity=($referenceMode eq 'NORM')? 0 : 1;
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
				my $pepIonKey="$pepSeq:$vmod:$charge";
				$bestValue=0 if $pepChargeState eq 'all'; # reset bestValue between charges
				foreach my $dataSrc (keys %{$proteins{$protID}{$pepSeq}{$vmod}{$charge}}) {
					if ($pepChargeState eq 'best' || $pepSource eq 'best') { # some data filtering applies
						my $setValue=0;
						foreach my $qSet (keys %{$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}}){
							next if $qSetIsBad{$qSet}; # only for SILAC
							foreach my $observationSet (keys %{$intensity{$pepIonKey}{$dataSrc}{$qSet}}){
								foreach my $refIntens (values %{$intensity{$pepIonKey}{$dataSrc}{$qSet}{$observationSet}}) { # PEP or TDA=>Sum accross all fragments
									$setValue+=$refIntens->[0]; # take first! => TODO: find a way to deal with multiple values here!!!
								}
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

	if ($labeling eq 'FREE' || $algoType eq 'PEP_INTENSITY') { # potential data filtering --->
		my %excludedData;

		##>Peptides must match across conditions (Never for PTM quantif)
		if ($matchingPep) { # && $ratioType eq 'SimpleRatio'
			my (%pepMeanValues,%checkExclusion);
			foreach my $pepSeq (keys %usableData) {
				foreach my $vmod (keys %{$usableData{$pepSeq}}) {
					foreach my $charge (sort{$a<=>$b} keys %{$usableData{$pepSeq}{$vmod}}) {
						my $pepIonKey="$pepSeq:$vmod:$charge";
						next unless $intensity{$pepIonKey};
						my (%pepIonValues,%pepIonInCond);
						foreach my $dataSrc (keys %{$usableData{$pepSeq}{$vmod}{$charge}}) {
							next unless $intensity{$pepIonKey}{$dataSrc};
							foreach my $qSet (keys %{$intensity{$pepIonKey}{$dataSrc}}){
								foreach my $observationSet (keys %{$intensity{$pepIonKey}{$dataSrc}{$qSet}}) {
									#$pepIonValues{"$qSet:#:$observationSet:#:$pepSeq:#:$vmod:#:$charge:#:$dataSrc"}=$intensity{$pepIonKey}{$dataSrc}{$qSet}{$observationSet}[0];
									foreach my $refIntens (values %{$intensity{$pepIonKey}{$dataSrc}{$qSet}{$observationSet}}) { # PEP or TDA/DIA=>Sum accross all fragments
										$pepIonValues{"$qSet:#:$observationSet:#:$pepSeq:#:$vmod:#:$charge:#:$dataSrc"}+=$refIntens->[0];
									}
									$pepIonInCond{$obs2Cond{$observationSet}}++; # num pep in condition
								}
							}
						}
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
			if ($topN && !$isModifQuantif) { # keep top<N>
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
						my $pepIonKey="$pepSeq:$vmod:$charge";
						next unless $intensity{$pepIonKey};
						foreach my $dataSrc (keys %{$usableData{$pepSeq}{$vmod}{$charge}}) {
							next unless $intensity{$pepIonKey}{$dataSrc};
							#foreach my $qSet (keys %{$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}}){ #}
							#	next unless $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet};
							foreach my $qSet (keys %{$intensity{$pepIonKey}{$dataSrc}}){
								foreach my $observationSet (keys %{$intensity{$pepIonKey}{$dataSrc}{$qSet}}) {
									#my $newValue=$intensity{$pepIonKey}{$dataSrc}{$qSet}{$observationSet}[0];
									my $newValue=0;
									foreach my $refIntens (values %{$intensity{$pepIonKey}{$dataSrc}{$qSet}{$observationSet}}) { # PEP or TDA=>Sum accross all fragments
										$newValue+=$refIntens->[0];
									}
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
			my $pepIonKey="$pepSeq:$vmod:$charge";
			delete $intensity{$pepIonKey}{$dataSrc}{$qSet}{$observationSet};
			unless (scalar keys %{$intensity{$pepIonKey}{$dataSrc}{$qSet}}) {
				delete $intensity{$pepIonKey}{$dataSrc}{$qSet};
				unless (scalar keys %{$intensity{$pepIonKey}{$dataSrc}}) {
					delete $intensity{$pepIonKey}{$dataSrc};
					unless (scalar keys %{$intensity{$pepIonKey}}) {
						delete $intensity{$pepIonKey};
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

	##>Printing data to file
	$protNumber++ if scalar keys %usableData;
	foreach my $pepSeq (sort{$peptideBeg{$protID}{$a}<=>$peptideBeg{$protID}{$b} || length($a)<=>length($b)} keys %usableData) {
		my @seq=split(//,$pepSeq);
		foreach my $vmod (keys %{$usableData{$pepSeq}}) {
			my $nonAmbigModProtID=$protID; # default
			my $ambigModProtID=$protID; # default
			my $seqIsAmbig=0;
# (PP 2017/06/22) --->
			if ($isModifQuantif) { # modif quantification
my @pos;
				my $firstModQ=1;
				foreach my $modID (@modifQuantifs) {
					
					my ($sep,$posStrg)=($vmod=~/(^|&)$modID:([\d+\.\-\=\+\*]+)/); # added N/C-term
					next unless $posStrg; # in multi-modif quantif context: varmod does not necessarly contains all quantified modif
					#>Non-ambiguous pos
					my @pos=(); # reset
					foreach my $pepPos (split(/\./,$posStrg)) {
						if ($pepPos=~/\d+/) {push @pos,$seq[$pepPos-1].($peptideBeg{$protID}{$pepSeq}+$pepPos-1);}
						#elsif ($pepPos=~/[=\*]/) {next;} # skip peptide N/C-term: cannot be in protein TODO: add to table.txt for data normalisation but do not quantify
						elsif ($pepPos eq '=') {push @pos,'n'.$peptideBeg{$protID}{$pepSeq};} # peptide N-term: kept for normalization purpose
						elsif ($pepPos eq '-') {push @pos,'n0';} # protein N-term: kept for normalization purpose
						elsif ($pepPos eq '*') {push @pos,'c'.($peptideBeg{$protID}{$pepSeq}+$#seq);} # peptide C-term: kept for normalization purpose
						elsif ($pepPos eq '+') {push @pos,'c0';} # protein C-term: kept for normalization purpose
					}
					#>THIS WHERE THE PROTEIN & MODIFICATION SITE(S) ARE CODED!!!!!!
					# Single-modif quantif:
					#	-Non-ambiguous: "ProtID-Res1Pos1[.Res2Pos2...]" eg "1234-S56.T78"
					#	-Ambiguous: "ProtID-StartRes~EndRes:NumModifs/NumAcceptors" eg "1234-56~78:2/5"
					# Mutli-modif quantif:
					#	-Non-ambiguous: "ProtID-ModifXRank#Res1Pos1[.Res2Pos2...][&ModifYRank#Res3Pos3[.Res4Pos4...]]" eg "1234-1#S56.T78&2#C61"
					#	-Ambiguous: "ProtID-ModifXRank#StartRes~EndRes:NumModifs/NumAcceptors" eg "1234-1#56~78:2/5&2#C61"
					my $nonAmbigMod=($quantifiedModifID)? join('.',@pos) : $modifQuantifRank{$modID}.'#'.join('.',@pos); # Add '<modif rank>#' if multi-modif quantif
					$nonAmbigModProtID.=($firstModQ)? '-' : '&';
					$nonAmbigModProtID.=$nonAmbigMod;
					
					# ==> Look for ambiguity
					$ambigModProtID.=($firstModQ)? '-' : '&';
					$ambigModProtID.=$modifQuantifRank{$modID}.'#' unless $quantifiedModifID;
					#my $numModifRes = () = $posStrg =~ /(\d+)/g;
					my $numModifRes=scalar @pos;
					my @targetResIdx;
					while ($pepSeq =~ /[$qQuantifResMatchStrg{$modID}]/g) {push @targetResIdx,$-[0];} # position of match (1st=0)
					my $numTargetRes=scalar @targetResIdx;
					my ($NtermMatch,$CtermMatch)=('','');
					if ($quantifResMatchStrg{$modID}=~/=/) {$NtermMatch='n'.$peptideBeg{$protID}{$pepSeq};} # peptide N-term
					elsif ($quantifResMatchStrg{$modID}=~/\*/) {$CtermMatch='c'.($peptideBeg{$protID}{$pepSeq}+$#seq);} # peptide C-term
					elsif ($quantifResMatchStrg{$modID}=~/-/) {$NtermMatch='n0' if $peptideBeg{$protID}{$pepSeq}==1;} # protein N-term
					#elsif ($quantifResMatchStrg=~/\+/) {$CtermMatch='c0' if $peptideBeg{$protID}{$pepSeq}+length($pepSeq)-1==<protein length>;} # TODO: fetch prot length
					my $numAlltargets=$numTargetRes;
					$numAlltargets++ if $NtermMatch;
					$numAlltargets++ if $CtermMatch;
					if ($numModifRes < $numAlltargets) { # Ambiguity! protID-[modifPos=]<1stPos>~<lastPos>:<numModif>/<numSites>
						$seqIsAmbig=1;
						if ($ambiguousModifPos) {
							#$ambigModProtID=$protID.'-';
							if ($NtermMatch) {$ambigModProtID.=$NtermMatch;} # N-term
							else {$ambigModProtID.=($peptideBeg{$protID}{$pepSeq}+$targetResIdx[0]);} # internal pos
							$ambigModProtID.='~';
							if ($CtermMatch) {$ambigModProtID.=$CtermMatch;} # C-term
							else {$ambigModProtID.=($peptideBeg{$protID}{$pepSeq}+$targetResIdx[-1]);} # internal pos
							$ambigModProtID.=':'.$numModifRes.'/'.$numAlltargets;
						}
						else {$ambigModProtID.=$nonAmbigMod;}
					}
					else {$ambigModProtID.=$nonAmbigMod;}
					$firstModQ=0;
				}
# <--- (PP 2017/06/22)
			}
			#my %numPepInCond4ModProt;
			foreach my $charge (sort{$a<=>$b} keys %{$usableData{$pepSeq}{$vmod}}) {
				my $pepIonKey="$pepSeq:$vmod:$charge";
				foreach my $dataSrc (sort keys %{$usableData{$pepSeq}{$vmod}{$charge}}) {
					my $infRatioDataSrc=($labeling eq 'FREE')? '-' : $dataSrc; # ignore dataSrc if label-free

					foreach my $qSet (sort{$a cmp $b} keys %{$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}}) { # sort cmp because Qset can be string
						next if $qSetIsBad{$qSet}; # only for SILAC or duplicate MCQ
						my $modProtID=$nonAmbigModProtID; # default
						my $usedVmod=($labeling eq 'SILAC')? $qSetBestVarMod{$qSet}[1] : $vmod; # can be different from $vmod in case qSet incoherence (eg. phospho pos after phosphoRS update on ambiguous pepSeq)
						if ($seqIsAmbig && $labeling eq 'SILAC') {
							if ($ptmProbSoftCode) {
								if ($qSetBestVarMod{$qSet}[0] >= $ptmProbThreshold) { # check Phospho pos coherence inside qSet (=0 if no PRS data)
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
											$qSetBestVarMod{$qSet}[2]=$protID;
											my $firstModQ=1;
											foreach my $modID (@modifQuantifs) {
											
												my ($sep,$posStrg)=($usedVmod=~/(^|&)$modID:([\d+\.\-\=\+\*]+)/); # (PP 2017/06/22)
												next unless $posStrg; # in multi-modif quantif context: varmod does not necessarly contains all quantified modif
												my @pos;
												#foreach my $pos (split(/\./,$posStrg)) {push @pos,$seq[$pos-1].($peptideBeg{$protID}{$pepSeq}+$pos-1);}
												foreach my $pos (split(/\./,$posStrg)) {
if (!$pos || ($pos=~/\d+/ && (!$seq[$pos-1])) || !$peptideBeg{$protID}{$pepSeq}) { # (PP 2017/06/22)
	warn "PTM position Error: PROT=$protID, PEP=$pepSeq, VMOD=$vmod, UVMOD=$usedVmod, QSET=$qSet, SEP=$sep, STRG=$posStrg, BEG=$peptideBeg{$protID}{$pepSeq}\n";
	&checkForErrors($dbh);
}
													push @pos,$seq[$pos-1].($peptideBeg{$protID}{$pepSeq}+$pos-1);
												}
												#>ALSO HERE => MODIFICATION SITE(S) ARE CODED!!!!!!
												$qSetBestVarMod{$qSet}[2].=($firstModQ)? '-' : '&';
												$qSetBestVarMod{$qSet}[2].=($quantifiedModifID)? join('.',@pos) : $modifQuantifRank{$modID}.'#'.join('.',@pos); # Add '<modif rank>:' if multi-modif quantif
												###$usedVmod=$qSetBestVarMod{$qSet}[1];
												
												$firstModQ=0;
											}
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

#						if ($ratioType eq 'Ratio') {
#							#my $pepIdStrg=join(';',@{$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}{$qSet}{'ALL_OBS'}});
#my $pepIdStrg='';
#foreach my $observationSet (@quantiOrder) {
#	if ($proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}{$qSet}{$observationSet}) {
#		$pepIdStrg.=';' if $pepIdStrg;
#		$pepIdStrg.=join(';',@{$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}{$qSet}{$observationSet}});
#	}
#}
#							print DATA "\n$modProtID\t$pepSeq\t$usedVmod\t$charge\t$pepIdStrg\t$proteinAlias{$protID}\t1"; #$validProt
#							foreach my $observationSet (@quantiOrder) { # order matters IMPORTANT: only 1 QSET per line for SILAC even with replicates to keep "linked" observations ratios computation (a/c)+(b/d) instead of (a+b)/(c+d) <- later good only for label-free
#								#if ($visInReplicate{$protID}{$observationSet} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}) {#}
#								if ($intensity{"$pepSeq:$vmod:$charge"}{$dataSrc} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}) {
#									print DATA "\t",$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0];
#									$pepIsComplete{$obs2Cond{$observationSet}}{$observationSet} = 1;
#									#$numPepInCond{'NA'}{$obs2Cond{$observationSet}}++;
#									#$numPepInCond4ModProt{$modProtID}{'NA'}{$obs2Cond{$observationSet}}++ if $quantifiedModifID;
#									if ($quantifiedModifID) {$numPepInCond4ModProt{$modProtID}{'NA'}{$obs2Cond{$observationSet}}{"$pepSeq:$usedVmod:$charge:$infRatioDataSrc"}=$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet};}
#									else {$numPepInCond{'NA'}{$obs2Cond{$observationSet}}{"$pepSeq:$usedVmod:$charge:$infRatioDataSrc"}=$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet};}
#								}
#								else {
#									print DATA "\tNA";
#									$pepIsComplete{$obs2Cond{$observationSet}}{$observationSet} = 0;
#								}
#								$numPepValues++;
#							}
#							# checking if peptide has a value in all conditions
#							my $pepIsOk = 1;
#							while (my ($condition, $completeChan) = each %pepIsComplete) {
#								unless (scalar grep { $_ == 1 } values %{$completeChan}) {
#									$pepIsOk = 0;
#									last;
#								}
#							}
#							$countCompletePeptides += $pepIsOk;
#						}
#						else { # (Super/Simple)Ratio or None(Abundance)
#my ($anaID)=$qSet=~/^A(\d+)/;
							foreach my $observationSet (@quantiOrder) {
								my $itraqStrg=($labeling=~/ITRAQ|TMT/)? $observationInfo{$observationSet}[1] : ''; # Fake SILAC for ITRAQ: ..._repX (PP Tibor set1)
								#if ($visInReplicate{$protID}{$observationSet} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet} && $intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}) {#}
								if ($intensity{$pepIonKey}{$dataSrc} && $intensity{$pepIonKey}{$dataSrc}{$qSet} && $intensity{$pepIonKey}{$dataSrc}{$qSet}{$observationSet}) {
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
									my $sep=($labeling eq 'FREE')? '+' : '='; # XICs are summed for label FREE
									my $pepIdStrg=join($sep,@{$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}{$qSet}{$observationSet}}); # used to be ';'
									#print DATA "\n$protID\t$pepSeq${vmod}_${charge}_$dataSrc\t$observationInfo{$observationSet}[0]\tNA\t$observationInfo{$observationSet}[1]\tA$observationInfo{$observationSet}[3]_Q$qSet\t$pepIdStrg\t$proteinAlias{$protID}\t1\t",$intensity{"$pepSeq:$vmod:$charge"}{$dataSrc}{$qSet}{$observationSet}[0];
									
									#print DATA "\n$modProtID\t$pepSeq$usedVmod\_$charge\_$dataSrc\t$observationInfo{$observationSet}[0]\tNA\t$observationInfo{$observationSet}[1]\t$observationInfo{$observationSet}[2]\t$expCode\t$qSet$itraqStrg\t$pepIdStrg\t$proteinAlias{$protID}\t1\t",$intensity{$pepIonKey}{$dataSrc}{$qSet}{$observationSet}[0]; #,"\t$usablePeptide";
##foreach my $feature (keys %{$intensity{$pepIonKey}{$dataSrc}{$qSet}{$observationSet}}) {
##	my $fragStrg=($feature eq 'PEP')? '' : '&'.$feature;
##	print DATA "\n$modProtID\t$pepSeq$usedVmod\_$charge$fragStrg\_$dataSrc\t$observationInfo{$observationSet}[0]\tNA\t$observationInfo{$observationSet}[1]\t$observationInfo{$observationSet}[2]\t$expCode\t$qSet$itraqStrg\t$pepIdStrg\t$proteinAlias{$protID}\t1\t",$intensity{$pepIonKey}{$dataSrc}{$qSet}{$observationSet}{$feature}[0]; #,"\t$usablePeptide";
##}
my $pepIntensity=0;
foreach my $feature (keys %{$intensity{$pepIonKey}{$dataSrc}{$qSet}{$observationSet}}) {
	$pepIntensity+=$intensity{$pepIonKey}{$dataSrc}{$qSet}{$observationSet}{$feature}[0];
}
print DATA "\n$modProtID\t$pepSeq$usedVmod\_$charge\_$dataSrc\t$observationInfo{$observationSet}[0]\tNA\t$observationInfo{$observationSet}[1]\t$observationInfo{$observationSet}[2]\t$expCode\t$qSet$itraqStrg\t$pepIdStrg\t$proteinAlias{$protID}\t$protValidity\t$pepIntensity"; #,"\t$usablePeptide";
$statesInFile{$observationInfo{$observationSet}[0]}=1; # number of states actually written in table(Ref).txt
									
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
									if ($ratioType eq 'SuperRatio' || $algoType eq 'PEP_RATIO') { # WARNING: Labeled w/ PEP_INTENSITY are not treated the same way (if more than 2 conds)
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
									if ($ratioType ne 'None') { # S\w+Ratio
										foreach my $context (keys %contextRatios) {
											if ($isModifQuantif) {$numPepInCond4ModProt{$modProtID}{$context}{$obs2Cond{$observationSet}}{"$pepSeq:$usedVmod:$charge:$infRatioDataSrc"}=$intensity{$pepIonKey}{$dataSrc}{$qSet}{$observationSet};} # Perl ref
											else {
												$numPepInCond{$context}{$obs2Cond{$observationSet}}{"$pepSeq:$usedVmod:$charge:$infRatioDataSrc"}=$intensity{$pepIonKey}{$dataSrc}{$qSet}{$observationSet}; # Perl ref
											}
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
#						}

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
	if ($ratioType ne 'None') {
		if ($quantifParameters{'DB'}{'MIN_INF_RATIOS'}[0]==0 && !$referenceMode) { ##### $ratioType=~/S\w+Ratio/ &&
			if ($isModifQuantif) {
				foreach my $modProtID (keys %numPepInCond4ModProt) {
					&checkInfiniteRatios($modProtID,$numPepInCond4ModProt{$modProtID});
				}
			}
			else {
				&checkInfiniteRatios($protID,\%numPepInCond);
			}
		}
	}
}
close DATA;
if (scalar keys %statesInFile < 2) { # at least 2 states otherwise stop everything on error
	my $fileType=($referenceMode eq 'REF')? 'protein reference' : ($referenceMode eq 'NORM')? 'normalizing peptide' : 'data';
	warn "ERROR: Less than 2 state data found in $fileType file!";
	&checkForErrors($dbh);
}

if ($referenceMode) { ### <- End of Reference mode !
	$dbh->disconnect;
#close DEBUG;
	exit;
}
elsif ($quantifParameters{'DB'}{'PEPTIDES_NORM'}) { ### <- Merge tableNorm.txt with table.txt
	system "cat $dataDir/tableNorm.txt >> $dataFile";
	unlink "$dataDir/tableNorm.txt";
}
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

#&promsQuantif::writeQuantifParameterFiles($dataDir,$quantifParameters{'R'});

&checkForErrors($dbh);

$dbh->disconnect;
#die "Data files are ready!!!"; # DEBUG

####>Attempt to free-up memory
undef %observationInfo;
undef %poolGroup2Cond;
undef %ana2Obs;
undef %cond2Obs;
undef %obs2Cond;
undef %labelModifs;
undef %anaLabelModifs;
undef %xicEngine;
undef %ana2Experiment;
undef %anaConds;
undef @conditionList;
undef %forNormalization;
undef %notForNormalization;
undef %manualPeptides;
undef %usedParamIDs;
undef @selectedPTMs;
undef %allowedPtmID;
undef %proteins;
undef %intensity;
undef %excludedSeq;
undef %proteinInConds;
undef %peptideBeg;
undef %qSetSequence;
undef %qSetBestVarMod;
undef %qSetIsBad;
undef %tempIntensity;
undef %dataSource2ObsSet;
undef %qSet2Protein;
undef %notUsed4Quantif;
undef %pepQuantifValues;
undef %recreatedVarMods;
undef %noSuperRatioProteins; # used in sub


##################>ENDING SCRIPT HERE IF ABUNDANCE NORMALIZATION<#######################
if ($ratioType eq 'None') {
	$dbh->disconnect;
	exit;
}


								############################################
			#############################> RUNNING R SCRIPTS <##########################
								############################################

&launchRscript;
&checkRsuccess;

#####> 4th: Launch R script
#my $progressStrg=($ratioType=~/S\w+Ratio/)? "2/$numSteps Running quantification" : "2/$numSteps Running normalization";
#open(FILESTAT,">>$fileStat");
#print FILESTAT "$progressStrg\n";
#close FILESTAT;
#
#my $RcommandString='';
##if ($ratioType eq 'Ratio') {
##	$RcommandString="export LANG=en_US.UTF-8; cd $runDir; $pathR/R CMD BATCH --no-save --no-restore $promsPath{R_scripts}/LabeledAndUnlabeledQuanti.Function.R";
##}
##else { #  (Super/Simple)Ratio
#	#foreach my $scriptName ('AnalysisDiffLimma.R','FunctionLimma.R','AffectDefaultValuesLimma.R','AffectParametersLimma.R') {
#	#	symlink("$promsPath{R_scripts}/$scriptName","$runDir/$scriptName");
#	#}
#	open(R_SCRIPT,">$runDir/AnalysisDiffLimma.R");
#	print R_SCRIPT qq
#|
####################################################################################
## Launcher for quantification R scripts (provides path for sourcing dependencies) #
####################################################################################
#
#filepath <- "$promsPath{R_scripts}/"
#source(paste(filepath,"AnalysisDiffLimma.R",sep=""))
#|;
#	close R_SCRIPT;
#	$RcommandString="export LANG=en_US.UTF-8; cd $runDir; $pathR/R CMD BATCH --no-save --no-restore AnalysisDiffLimma.R";
##}
#
#### ALWAYS SKIP CLUSTER (PP 17/09/18)!!! ###
####if ($cluster{'on'}) {
####	my $bashFile = "$runDir/runXICProtQuantification.sh";
####	my $clusterRcommandString=$cluster{'buildCommand'}->($runDir,$RcommandString);
####	my $maxHours=int(0.5+($numPepValues/25000)); $maxHours=2 if $maxHours < 2; $maxHours=48 if $maxHours > 48; # 1 M lines -> 40 h
####	#my $maxMem=($numPepValues < 250000)? int(0.5+($numPepValues/25000)) : int(0.5+($numPepValues/8000)); $maxMem=1 if $maxMem < 1; $maxMem.='gb'; # 1 M lines -> 80 gb
####	#my $maxMem='2gb';
####	# RAM = 8E-7 * <num lines> + 0.1614
####	#my $maxMem=0.5 + 8E-7 * $numPepValues;
####	#my $maxMem=0.5 + 1E-6 * $numPepValues;
####	my $maxMem=int(1.5 + 1E-6 * $numPepValues);
####	$maxMem=50 if $maxMem > 50;
####	$maxMem.='Gb';
####	open (BASH,">$bashFile");
####	print BASH qq
####|#!/bin/bash
######resources
#####PBS -l mem=$maxMem
#####PBS -l nodes=1:ppn=1
#####PBS -l walltime=$maxHours:00:00
#####PBS -q batch
######Information
#####PBS -N myProMS_protQuant_$quantifID
#####PBS -m abe
#####PBS -o $runDir/PBS.txt
#####PBS -e $runDir/PBSerror.txt
#####PBS -d $runDir
####
###### Command
####$clusterRcommandString
####echo _END_$quantifID
####|;
#####close DEBUG;
####	close BASH;
####	my $modBash=0775;
####	chmod $modBash, $bashFile;
####	#system "$promsPath{qsub}/qsub $bashFile > $runDir/torqueID.txt";
####	#system "$cluster{command} $bashFile > $runDir/torqueID.txt 2> $runDir/connection_error.txt";
####	#system "ssh -q -o \"UserKnownHostsFile=/dev/null\" -o \"StrictHostKeyChecking no\" calcsub \"qsub $bashFile > $runDir/torqueID.txt\" 2>$runDir/ssh_error.txt ";
####	$cluster{'sendToCluster'}->($bashFile); # bash file, qsub output file[, ssh error file (only for centos cluster)]
####	sleep 30;
####
####	###>Waiting for R job to run
####	my $pbsError;
####	my $nbWhile=0;
####	my $maxNbWhile=$maxHours*60*2;
####	while ((!-e "$runDir/PBS.txt" || !`tail -3 $runDir/PBS.txt | grep _END_$quantifID`) && !$pbsError) {
####		if ($nbWhile > $maxNbWhile) {
####			$dbh=&promsConfig::dbConnect('no_user'); # reconnect
####			$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID"); # Failed
####			$dbh->commit;
####			$dbh->disconnect;
####			die "Aborting quantification: R is taking too long or died before completion";
####		}
####		sleep 30;
####		#$pbsError=`head -5 $runDir/PBSerror.txt` if -e "$runDir/PBSerror.txt";
####		$pbsError=$cluster{'checkError'}->("$runDir/PBSerror.txt");
####		$nbWhile++;
####	}
####}
####else { ###>Run job on Web server
##	system $RcommandString;
####}
#sleep 3;
#
#$dbh=&promsConfig::dbConnect('no_user'); # reconnect
#
#####>ERROR Management<####
##my $RoutFile=($ratioType eq 'Ratio')? 'LabeledAndUnlabeledQuanti.Function.Rout' : 'AnalysisDiffLimma.Rout';
#my $RoutFile='AnalysisDiffLimma.Rout';
#my $RoutStrg=`tail -3 $runDir/$RoutFile`;
#unless ($RoutStrg=~/proc\.time\(\)/) {
#	#$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID");
#	#$dbh->commit;
#	#$dbh->disconnect;
#	$RoutStrg=`tail -20 $runDir/$RoutFile`;
#	my $RerrorStrg="R script has generated an error!";
#	my $inError=0;
#	foreach my $line (split(/\n/,$RoutStrg)) {
#		next if (!$inError && $line !~ /^Error in/); # skip everything before "Error in..."
#		$inError=1;
#		$RerrorStrg.="\n$line";
#	}
#	warn $RerrorStrg;
#	&checkForErrors($dbh);
#}

#die "R quantification is finished!!!"; # DEBUG


									############################################
			#############################> PARSING QUANTIFICATION RESULTS <########################## ONLY Super/SimpleRatio from now on
									############################################

open(FILESTAT,">>$fileStat");
print FILESTAT "3/3 Parsing results\n";
close FILESTAT;
#$dbh->disconnect; exit; # DEBUG

my $numQuantifRatios=scalar @ratioTarget; # =numConditions*0.5*(numConditions-1)
my $adjStrg=($ratioType=~/S\w+Ratio/ || ($quantifParameters{'DB'}{'FDR_CONTROL'} && $quantifParameters{'DB'}{'FDR_CONTROL'}[0] eq 'TRUE'))? '_ADJ' : '';

####>Fetching list of quantification parameters<####
$dbh=&promsConfig::dbConnect('no_user'); # reconnect
my %quantifParamIDs;
my $sthQP=$dbh->prepare("SELECT QP.ID_QUANTIF_PARAMETER,QP.CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD QM,QUANTIFICATION_PARAMETER QP WHERE Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND QM.ID_QUANTIFICATION_METHOD=QP.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=$quantifID");
$sthQP->execute;
while (my ($paramID,$code)=$sthQP->fetchrow_array) {
    $quantifParamIDs{$code}=$paramID;
}
$sthQP->finish;

my $sthInsProt=$dbh->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_QUANTIFICATION,ID_PROTEIN,ID_QUANTIF_PARAMETER,QUANTIF_VALUE,TARGET_POS) VALUES ($quantifID,?,?,?,?)"); # PK autoincrement
my ($sthInsModRes,$sthInsProtRes);
if ($isModifQuantif) {
	$sthInsModRes=($quantifiedModifID)? $dbh->prepare("INSERT INTO MODIFIED_RESIDUE (ID_QUANTIFICATION,RESIDUE,POSITION) VALUES ($quantifID,?,?)")
									  : $dbh->prepare("INSERT INTO MODIFIED_RESIDUE (ID_QUANTIFICATION,RESIDUE,POSITION,MODIF_RANK) VALUES ($quantifID,?,?,?)");
	$sthInsProtRes=$dbh->prepare("INSERT INTO PROTQUANTIF_MODRES (ID_MODIF_RES,ID_PROT_QUANTIF) VALUES (?,?)");
}



my (%protQuantified,%modResiduesID,%lostProteins);

###################################>(Super/Simple)Ratio (ML->AS->IB scripts)<#############################################
if ($ratioType=~/S\w+Ratio/) {
	
	open(FILESTAT,">>$fileStat");
	print FILESTAT "3/$numSteps Parsing results\n";
	close FILESTAT;
	
	################################################################
	####>Recording protein ratios (& p-values, Std. dev if any)<####
	################################################################

	####>Recording infinite ratios ratios<#### because prot might not be kept in quantif result files
	if ($isModifQuantif) { # modif quantification
		foreach my $modProtID (sort{&promsMod::sortSmart($a,$b)} keys %infiniteRatioProteins) {
			my ($protID,@modQuantifs)=split(/[-&]/,$modProtID);
			my @allModResidues;
			foreach my $modQuantif (@modQuantifs) {
				my @modResidues;
				if ($quantifiedModifID) {@modResidues=split('\.',$modQuantif);}
				else {
					(my $modifRank,@modResidues)=split(/[#.]/,$modQuantif);
					foreach my $i (0..$#modResidues) {$modResidues[$i]=$modifRank.'#'.$modResidues[$i];} # Convert '1#V32.K44' into ['1#V32','1#K44']
				}
				&promsQuantif::insertModifiedResidues($protID,\@modResidues,\%modResiduesID,$dbh,$sthInsModRes); # $modifRank undef for single-modif quantif
				push @allModResidues,@modResidues;
			}
			foreach my $ratioPos (sort{$a<=>$b} keys %{$infiniteRatioProteins{$modProtID}}){
				$sthInsProt->execute($protID,$quantifParamIDs{'RATIO'},$infiniteRatioProteins{$modProtID}{$ratioPos},$ratioPos);
				&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF'));
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

	#unlink "$runDir/Rplots.pdf" if -e "$runDir/Rplots.pdf"; # temporary ?!!!!!!!!!!!!!!!!!!!!!!!!!!

	###>Generating ratio codes used by R<###
	my (%measurePos2Code,%measureCode2RatioPos,%state2RatioPos); # a state can be involved in multiple ratios #,%experiment2RatioPos
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
			if ($algoType =~ /PEP_INTENSITY|TDA|DIA/) {
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
	elsif ($algoType =~ /PEP_INTENSITY|TDA|DIA/) { # also record state mean
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
		my $protID=$protID0; # default
		my @allModResidues;
		if ($isModifQuantif) { # modif quantification
			next if $protID0 !~ /-/; # must be a mod prot
			($protID,my @modQuantifs)=split(/[-&]/,$protID0);
			foreach my $modQuantif (@modQuantifs) {
				my @modResidues;
				if ($quantifiedModifID) {@modResidues=split('\.',$modQuantif);}
				else {
					(my $modifRank,@modResidues)=split(/[#.]/,$modQuantif);
					foreach my $i (0..$#modResidues) {$modResidues[$i]=$modifRank.'#'.$modResidues[$i];} # Convert '1#V32.K44' into ['1#V32','1#K44']
				}
				&promsQuantif::insertModifiedResidues($protID,\@modResidues,\%modResiduesID,$dbh,$sthInsModRes); # <-- %modResiduesID should be full at <PROT1>
				push @allModResidues,@modResidues;
			}
		}
		
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
				&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
				$protQuantified{$protID0}{$targetPos}=1;
			}
			else { # mean state (LABELFREE)
next unless defined $protHeader2Idx{'Log2Mean_'.$measurePos2Code{$targetPos}}; # temp AS should had this again soon (PP 19/10/17)
				my $log2Mean=$lineData[$protHeader2Idx{'Log2Mean_'.$measurePos2Code{$targetPos}}];
				next if $log2Mean=~/NA/i; # NA or NaN
				$sthInsProt->execute($protID,$quantifParamIDs{'MEAN_STATE'},2**$log2Mean,$targetPos); #  exp($log2Mean*$log2)
				&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif; # unlikely in LABELFREE...
			}
			#>Associated parameters
			foreach my $refParam (@RDbCodeList) {
				my ($paramR,$paramDB)=@{$refParam};
				my $paramValue=$lineData[$protHeader2Idx{$paramR.'_'.$measurePos2Code{$targetPos}}];
				next if $paramValue=~/NA/i; # NA or NaN
				$sthInsProt->execute($protID,$quantifParamIDs{$paramDB},$paramValue,$targetPos);
				&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
			}
		} # next ratioPos

	} # next line
	close PROT;

	####>Computing peptide counts (NUM_PEP_TOTAL, NUM_PEP_USED & DIST_PEP_USED)<####
	my (%numPeptideSets,%usedPeptideSets,%distinctPeptideSets); # Set <=> num label channels (1 for Label-Free)
my (%usedPeptideReplicate,%distinctPeptideReplicate,%processedReplicates);
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
		next if ($algoType eq 'PEP_RATIO' && !$measureCode2RatioPos{$condCode}); # ratio not used (SINGLE_REF) StateY.StateX (not applicable for LABELFREE: StateX\nStateY)
		my $seqVarMod=(split(/_/,$lineData[$pepHeader2Idx{'Peptide'}]))[0];
		#my $numPepSets=scalar (split(';',$pepIdStrg));
		my @pepIdList=split(/[\.\+\=]/,$lineData[$pepHeader2Idx{'PeptideId'}]);
		my %numPepSets;
		##my $numPepSets=1; # default
		# TODO: True for SILAC & other Isotopes (not LF or isobar) --->
		if ($labeling eq 'SILAC') {
			#if ($lineData[$pepHeader2Idx{'PeptideId'}]=~/\+/) { # XIC were summed
				##my %distinctQsets;
				##foreach my $pepID (@pepIdList) {
				##	$distinctQsets{$pepID2QuantifSetSILAC{$pepID}}=1;
				##}
				##$numPepSets=scalar keys %distinctQsets;
				foreach my $pepID (@pepIdList) {
					$numPepSets{$pepID2QuantifSetSILAC{$pepID}}=1;
				}
			#}
			#else {$numPepSets{$lineData[$pepHeader2Idx{'PeptideId'}]}=1;} #???? TODO: Check again for SILAC (PP 03/07/19)
		}
		else {
			foreach my $pepID (@pepIdList) {$numPepSets{$pepID}=1;}
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
			##$usedPeptideSets{$protID0}{$ratioPos}+=$numPepSets;
foreach my $pepIdCode (keys %numPepSets) {
	$usedPeptideSets{$protID0}{$ratioPos}{$pepIdCode}=1;
}
		}
foreach my $state (split(/\./,$condCode)) {	# eg "StateY.StateX" for PEP_RATIO or "StateX" for PEP_INTENSITY/TDA
	my $replicateTgtPos=$replicateTargetPos{ $state.'.'.$lineData[$pepHeader2Idx{'replicate'}].'.'.$lineData[$pepHeader2Idx{'repTech'}] };
	$distinctPeptideReplicate{$protID0}{$replicateTgtPos}{NORM}{$seqVarMod}=1; # DIST_PEP_USED for replicates
	foreach my $pepIdCode (keys %numPepSets) {
		$usedPeptideReplicate{$protID0}{$replicateTgtPos}{NORM}{$pepIdCode}=1;
	}
	$processedReplicates{"$protID0:$replicateTgtPos"}=1;
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
	my $prevExp='-';
	my $ratioPos=0;
	open(DATA,"$dataDir/table.txt");
	while(<DATA>) {
		next if $.==1;
		#chomp;
		my ($protID0,$peptideData,$state,$treatment,$bioRep,$techRep,$experiment,$qSet,$pepIdStrg,$protName,$validity,$xicValue)=split(/\t/,$_);
		next unless $validity;
		unless ($protQuantified{$protID0}) {
			$lostProteins{$protID0}=1;
			next;
		}
		my @pepIdList=split(/[\.\+\=]/,$pepIdStrg);
		my %numPepSets;
		if ($labeling eq 'SILAC') {
			#if ($pepIdStrg=~/\+/) {
				foreach my $pepID (@pepIdList) {
					$numPepSets{$pepID2QuantifSetSILAC{$pepID}}=1;
					$numPeptideSets{$protID0}{$pepID2QuantifSetSILAC{$pepID}}=1; # used for NUM_PEP_TOTAL
				}
			#}
			#else { #???? TODO: Check again for SILAC (PP 03/07/19)
			#	$numPepSets{$pepIdStrg}=1;
			#	$numPeptideSets{$protID0}{$pepIdStrg}=1;
			#}
		}
		else {
			foreach my $pepID (@pepIdList) {
				$numPepSets{$pepID}=1;
				$numPeptideSets{$protID0}{$pepID}=1;
			}
		}	
##$peptideData=~s/^(.+_\d+)&.+(_.+)/$1$2/ if $algoType eq 'TDA'; # remove fragment data
##$numPeptideSets{$protID0}{$peptideData}=1; # NUM_PEP_TOTAL (counted ONCE across all states) LABEL-FREE: requires that dataSrc tag != '-'!!!
		
#foreach my $pepID (split(/\+/),$pepIdStrg) { # v3: true number of peptides used (for all labeling channels)
#	$numPeptideSets{$protID0}{$pepID}=1; # NUM_PEP_TOTAL (counted ONCE across all states)
#}
		next unless $infiniteRatioProteins{$protID0}; # No source/replicate aggregation for NUM_PEP_USED in case of infinite ratio!!!
		#--- Only used for infinite ratios from here on ---->
		my $seqVarMod=(split(/_/,$peptideData))[0]; # $src,$fragment DIST_PEP_USED: $peptideCharge

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

		my $stateInInfRatio=0;
		foreach my $ratioPos (@{$refRatioList}) {
			next unless $infiniteRatioProteins{$protID0}{$ratioPos};
			#$usedPepSets{$protID0}{$ratioPos}{$peptideData}=1; # $peptideCharge
foreach my $pepIdCode (keys %numPepSets) {
	$usedPeptideSets{$protID0}{$ratioPos}{$pepIdCode}=1; # PP 20/09/19
}
			$distinctPeptideSets{$protID0}{$ratioPos}{$seqVarMod}=1;
			$stateInInfRatio=1; # state is involved in +/-Inf ratio
		}
my $replicateTgtPos=$replicateTargetPos{"$state.$bioRep.$techRep"};
if (!$processedReplicates{"$protID0:$replicateTgtPos"} || $stateInInfRatio) { # WARNING: same state can be found in normal and inf ratios: num pep per replicates calculated will be based on normal ratio (from resultsPep.txt) => without excluded pep
	my $ratioTypeKey=($stateInInfRatio)? 'INF' : 'NORM';
	$distinctPeptideReplicate{$protID0}{$replicateTgtPos}{$ratioTypeKey}{$seqVarMod}=1; # DIST_PEP_USED for replicates
	foreach my $pepIdCode (keys %numPepSets) {
		$usedPeptideReplicate{$protID0}{$replicateTgtPos}{$ratioTypeKey}{$pepIdCode}=1; # $peptideData
	}
}
	}
	close DATA;
	#foreach my $protID0 (keys %usedPepSets) { # inf ratios only
	#	foreach my $ratioPos (keys %{$usedPepSets{$protID0}}) {
	#		##$usedPeptideSets{$protID0}{$ratioPos}=scalar keys %{$usedPepSets{$protID0}{$ratioPos}};
	#		$usedPeptideSets{$protID0}{$ratioPos}=$usedPepSets{$protID0}{$ratioPos};
	#	}
	#}

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
						##$usedPeptideSets{$protID0}{$ratioPos}=($usedPeptideSets{$protID0}{$rPos1} < $usedPeptideSets{$protID0}{$rPos2})? $usedPeptideSets{$protID0}{$rPos1} : $usedPeptideSets{$protID0}{$rPos2}; # keep smallest (18/02/15)
						$usedPeptideSets{$protID0}{$ratioPos}=(scalar keys %{$usedPeptideSets{$protID0}{$rPos1}} < scalar keys %{$usedPeptideSets{$protID0}{$rPos2}})? $usedPeptideSets{$protID0}{$rPos1} : $usedPeptideSets{$protID0}{$rPos2}; # keep smallest (18/02/15)
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
		#my ($protID,@modResidues)=($quantifiedModifID)? split(/[-.]/,$protID0) : ($protID0,undef);
		my $protID=$protID0; # default
		my @allModResidues;
		if ($isModifQuantif) {
			($protID,my @modQuantifs)=split(/[-&]/,$protID0);
			foreach my $modQuantif (@modQuantifs) {
				my @modResidues;
				if ($quantifiedModifID) {@modResidues=split('\.',$modQuantif);}
				else {
					(my $modifRank,@modResidues)=split(/[#.]/,$modQuantif);
					foreach my $i (0..$#modResidues) {$modResidues[$i]=$modifRank.'#'.$modResidues[$i];} # Convert '1#V32.K44' into ['1#V32','1#K44']
				}
				push @allModResidues,@modResidues;
			}
		}
		# Count for ratio
		foreach my $ratioPos (keys %{$distinctPeptideSets{$protID0}}) {
			next if (!$protQuantified{$protID0} || !$protQuantified{$protID0}{$ratioPos});
			$sthInsProt->execute($protID,$quantifParamIDs{'DIST_PEP_USED'},scalar keys %{$distinctPeptideSets{$protID0}{$ratioPos}},$ratioPos);
			&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
		}
		next if !$protQuantified{$protID0};
		# Count for replicate
		foreach my $replicateTgtPos (keys %{$distinctPeptideReplicate{$protID0}}) {
my $valueNorm=scalar keys %{$distinctPeptideReplicate{$protID0}{$replicateTgtPos}{NORM}}; # always defined
my $valueInf=($distinctPeptideReplicate{$protID0}{$replicateTgtPos}{INF})? scalar keys %{$distinctPeptideReplicate{$protID0}{$replicateTgtPos}{INF}} : 0;
if (!$valueInf || $valueInf==$valueNorm) {			
			$sthInsProt->execute($protID,$quantifParamIDs{'DIST_PEP_USED'},$valueNorm,$replicateTgtPos);
			&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
}
if ($valueInf && $valueInf != $valueNorm) {
	$sthInsProt->execute($protID,$quantifParamIDs{'DIST_PEP_USED'},$valueInf,-$replicateTgtPos);
	&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
}
		}
	}
	#>NUM_PEP_USED
	foreach my $protID0 (keys %usedPeptideSets) {
		#my ($protID,@modResidues)=($quantifiedModifID)? split(/[-.]/,$protID0) : ($protID0,undef);
		my $protID=$protID0; # default
		my @allModResidues;
		if ($isModifQuantif) {
			($protID,my @modQuantifs)=split(/[-&]/,$protID0);
			foreach my $modQuantif (@modQuantifs) {
				my @modResidues;
				if ($quantifiedModifID) {@modResidues=split('\.',$modQuantif);}
				else {
					(my $modifRank,@modResidues)=split(/[#.]/,$modQuantif);
					foreach my $i (0..$#modResidues) {$modResidues[$i]=$modifRank.'#'.$modResidues[$i];} # Convert '1#V32.K44' into ['1#V32','1#K44']
				}
				push @allModResidues,@modResidues;
			}
		}
		# Count for ratio
		foreach my $ratioPos (sort keys %{$usedPeptideSets{$protID0}}) {
			next if (!$protQuantified{$protID0} || !$protQuantified{$protID0}{$ratioPos});
			##$sthInsProt->execute($protID,$quantifParamIDs{'NUM_PEP_USED'},$usedPeptideSets{$protID0}{$ratioPos},$ratioPos);
			$sthInsProt->execute($protID,$quantifParamIDs{'NUM_PEP_USED'},scalar keys %{$usedPeptideSets{$protID0}{$ratioPos}},$ratioPos);
			&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
		}
		next if !$protQuantified{$protID0};
		# Count for replicate
		foreach my $replicateTgtPos (keys %{$usedPeptideReplicate{$protID0}}) {
my $valueNorm=scalar keys %{$usedPeptideReplicate{$protID0}{$replicateTgtPos}{NORM}}; # always defined
my $valueInf=($usedPeptideReplicate{$protID0}{$replicateTgtPos}{INF})? scalar keys %{$usedPeptideReplicate{$protID0}{$replicateTgtPos}{INF}} : 0;
if (!$valueInf || $valueInf==$valueNorm) {
			$sthInsProt->execute($protID,$quantifParamIDs{'NUM_PEP_USED'},$valueNorm,$replicateTgtPos);
			&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
}
if ($valueInf && $valueInf != $valueNorm) {
	$sthInsProt->execute($protID,$quantifParamIDs{'NUM_PEP_USED'},$valueInf,-$replicateTgtPos);
	&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
}
		}
	}
	#>NUM_PEP_TOTAL
	foreach my $protID0 (keys %numPeptideSets) {
		next if !$protQuantified{$protID0};
		#my ($protID,@modResidues)=($quantifiedModifID)? split(/[-.]/,$protID0) : ($protID0,undef);
		my $protID=$protID0; # default
		my @allModResidues;
		if ($isModifQuantif) {
			($protID,my @modQuantifs)=split(/[-&]/,$protID0);
			foreach my $modQuantif (@modQuantifs) {
				my @modResidues;
				if ($quantifiedModifID) {@modResidues=split('\.',$modQuantif);}
				else {
					(my $modifRank,@modResidues)=split(/[#.]/,$modQuantif);
					foreach my $i (0..$#modResidues) {$modResidues[$i]=$modifRank.'#'.$modResidues[$i];} # Convert '1#V32.K44' into ['1#V32','1#K44']
				}
				push @allModResidues,@modResidues;
			}
		}
		$sthInsProt->execute($protID,$quantifParamIDs{'NUM_PEP_TOTAL'},scalar keys %{$numPeptideSets{$protID0}},undef); # no ratioPos
		&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
	}

}


#close DEBUG;

$sthInsProt->finish;
if ($isModifQuantif) { # modif quantification
	$sthInsModRes->finish;
	$sthInsProtRes->finish;
}

my $extraQuantifAnnot=($labeling eq 'FREE')? '' : '::QSET_USED='.(scalar keys %{$qSetUsage{'USED'}}).'/'.((scalar keys %{$qSetUsage{'USED'}}) + (scalar keys %{$qSetUsage{'EXCLUDED'}}));

&endQuantification($dbh,$projectID,\%lostProteins,\%proteinAlias,$extraQuantifAnnot);



############################
####<Launching R script>####
############################
sub launchRscript {
	my $progressStrg=($ratioType=~/S\w+Ratio/)? "2/$numSteps Running quantification" : "2/$numSteps Running normalization";
	open(FILESTAT,">>$fileStat");
	print FILESTAT "$progressStrg\n";
	close FILESTAT;

	&promsQuantif::writeQuantifParameterFiles($dataDir,$quantifParameters{'R'});

	open(R_SCRIPT,">$runDir/AnalysisDiffLimma.R");
	print R_SCRIPT qq
|
###################################################################################
# Launcher for quantification R scripts (provides path for sourcing dependencies) #
###################################################################################

filepath <- "$promsPath{R_scripts}/"
source(paste(filepath,"AnalysisDiffLimma.R",sep=""))
|;
	close R_SCRIPT;
	system "export LANG=en_US.UTF-8; cd $runDir; $pathR/R CMD BATCH --no-save --no-restore AnalysisDiffLimma.R";
	
	sleep 5;
}
sub checkRsuccess {
	####<ERROR Management>####
	my $RoutFile='AnalysisDiffLimma.Rout';
	my $RoutStrg=`tail -3 $runDir/$RoutFile`;
	unless ($RoutStrg=~/proc\.time\(\)/) {
		$RoutStrg=`tail -20 $runDir/$RoutFile`;
		my $RerrorStrg="R script has generated an error!";
		my $inError=0;
		foreach my $line (split(/\n/,$RoutStrg)) {
			next if (!$inError && $line !~ /^Error in/); # skip everything before "Error in..."
			$inError=1;
			$RerrorStrg.="\n$line";
		}
		warn $RerrorStrg;
		my $dbh=&promsConfig::dbConnect('no_user');
		&checkForErrors($dbh);
	}
}

###############################
####<Endind Quantification>####
###############################
sub endQuantification { # Globals: $quantifID,%promsPath,%quantifParameters,$runDir,$resultDir,$fileStat
	my ($dbh,$projectID,$refLostProteins,$refProteinAlias,$extraQuantifAnnot)=@_;
	$extraQuantifAnnot='' unless $extraQuantifAnnot;
	
	my $numLostProteins=scalar keys %{$refLostProteins}; # proteins in table.txt but not quantified at all!!!
	if ($numLostProteins) {
		open (LOST,">$resultDir/lostProteins.txt");
		print LOST "Protein_ID\tProtein_Name";
		foreach my $modProtID (sort{&promsMod::sortSmart($a,$b)} keys %{$refLostProteins}) {
			my ($protID)=$modProtID=~/^(\d+)/;
			print LOST "\n$modProtID\t$refProteinAlias->{$protID}";
		}
		close LOST;
	}
	&checkForErrors($dbh);
	
	###> Quantification is finished.
	my $status=($quantifParameters{'DB'}{'VISIBILITY'})? $quantifParameters{'DB'}{'VISIBILITY'}[0] : 1;
	$dbh->do("UPDATE QUANTIFICATION SET STATUS=$status,QUANTIF_ANNOT=CONCAT(QUANTIF_ANNOT,'$extraQuantifAnnot\:\:LOST_PROTEINS=$numLostProteins') WHERE ID_QUANTIFICATION=$quantifID");
	
	$dbh->commit;
	
	$dbh->disconnect;
	#exit; # DEBUG!!!
	
	open(FILESTAT,">>$fileStat");
	print FILESTAT "Ended ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
	close(FILESTAT);
	
	####> 7th: Move all the created files in a specific directory so as to clean the directory
	#if ($ratioType=~/S\w+Ratio|None/) {
		#foreach my $scriptName ('AnalysisDiffLimma.R','FunctionLimma.R','AffectDefaultValuesLimma.R','AffectParametersLimma.R') {
		#	unlink "$runDir/$scriptName";
		#}
		unlink "$runDir/AnalysisDiffLimma.R";
	#}
	
	mkdir "$promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
	dirmove($runDir,"$promsPath{quantification}/project_$projectID/quanti_$quantifID");
	
	sleep 2;
	
}

################################################################
####<Computing protein abundance (mean, median, sum & LFQ?)>####
################################################################
sub splitResultsPepFile {
	my ($numJobs)=@_;
	
	open(FILESTAT,">>$fileStat");
	print FILESTAT "3/$numSteps Computing protein abundances\n";
	close FILESTAT;
	
	my ($header,%pepHeader2Idx,%proteinData);
	open(PEP,"$resultDir/resultsPep.txt") || die $!;
	while(<PEP>) {
		my @lineData=split(/\t/,$_);
		if ($.==1) { # 1st line of the file
			$header=$_;
			my $colIdx=0;
			foreach my $colName (@lineData) {
				$pepHeader2Idx{$colName}=$colIdx;
				$colIdx++;
			}
			next;
		}
		my $protID0=$lineData[$pepHeader2Idx{'ProteinID'}];
		push @{$proteinData{$protID0}},$_;
	}
	close PEP;
	
	my %fileHandle;
	foreach my $job (1..$numJobs) {
		rmtree "$resultDir/ABUND$job" if -e "$resultDir/ABUND$job"; # just to be safe
		mkdir "$resultDir/ABUND$job";
		open ($fileHandle{$job},">$resultDir/ABUND$job/resultsPep.txt") || die $!;
		print {$fileHandle{$job}} $header; # contains "\n"
	}
	my $currentJob=1;
	foreach my $protID0 (keys %proteinData) { # protein data are cyclicly spread into all files
		foreach my $lineData (@{$proteinData{$protID0}}) {
			print {$fileHandle{$currentJob}} $lineData; # contains "\n"
		}
		$currentJob++;
		$currentJob=1 if $currentJob > $numJobs;
	}
	foreach my $job (1..$numJobs) {
		close $fileHandle{$job};
	}
	
	exit;
}
sub runProteinAbundance {
	my ($numJobs)=@_;
	
	####<Launching parallel LFQ jobs>####
	my %cluster=&promsConfig::getClusterInfo;
	if ($cluster{'on'}) { #$numJobs > 1 && 
		my $MAX_PARALLEL_QUANTIFS=$cluster{'maxJobs'};
#$MAX_PARALLEL_QUANTIFS=25;
		my $cgiUnixDir=`pwd`;
		$cgiUnixDir=~s/\/*\s*$//; # cd is required for script to find myproms .pm files!!!!

		my $startMem=int(1.5 * $numJobs**2 / 4500); # numJobs <=> num bio/tech replicates
		#my $startMem=int(1.5 * $numJobs / 12.5);
		#my $startMem=32; # in Gb <-------------------------
		my %allJobsParameters=(
			#maxMem=>$startMem.'Gb',
			numCPUs=>1,
			maxHours=>48,
			pbsRunDir=>$cgiUnixDir,
			noWatch=>1 # Do not auto watch job
		);
		
		my (%launchedJobs,	# total num jobs launched (finished + running)
			%runningJobs,	# num jobs currently running
			%triedJobs		# num tries for each job (max. 3/job allowed before erroring quantif)
			);
		my $numJobToRun=$numJobs;
		MAIN_LOOP:while ($numJobToRun) {

			###<Parallel launch of up to $MAX_PARALLEL_QUANTIFS jobs
			while (scalar (keys %runningJobs) < $MAX_PARALLEL_QUANTIFS) {
				
				#<Scan all job lists for job to be run (in case a previously launched job has failed and must be relaunched)
				my $currentJob=0;
				foreach my $job (1..$numJobs) {
					if (!$launchedJobs{$job}) {
						if (-e "$resultDir/ABUND$job/end.flag") { # in case using data from previous incomplete quantif
							$launchedJobs{$job}=1;
							$numJobToRun--;
							next;
						}
						$currentJob=$job;
						last;
					}
				}
				last unless $currentJob; # no more jobs to launch (but some jobs could still be running)

				#<Launchjob
				my $jobDir="$resultDir/ABUND$currentJob";
				system "echo LAUNCH > $jobDir/launch.flag"; # launch flag file
				my %jobParameters=%allJobsParameters;
				$jobParameters{'commandBefore'}="echo RUN > $jobDir/run.flag", # run flag file
				$jobParameters{'commandAfter'}="echo END > $jobDir/end.flag", # end flag file
				
				my $commandString="export LC_ALL=\"C\"; cd $cgiUnixDir; $cluster{path}{perl}/perl runXICProtQuantification.pl $quantifID $quantifDate ABUND_JOB:$currentJob"; # 2> $jobDir/errorJob.txt
				
				while (1) { # try launching job (max 3 times)
					$triedJobs{$currentJob}++;
					$jobParameters{'jobName'}="myProMS_protQuant_ABUND_$quantifID.$currentJob.$triedJobs{$currentJob}";
					$jobParameters{'maxMem'}=($startMem * $triedJobs{$currentJob}).'Gb';
					my ($pbsError,$pbsErrorFile,$jobClusterID)=$cluster{'runJob'}->($jobDir,$commandString,\%jobParameters);
					if ($pbsError) {
						if ($triedJobs{$currentJob} == 3) {
							die "Job #$currentJob ($jobClusterID): $pbsError\n"; # error will be handled by parent process (MASTER)
						}
						else { # clean dir & retry
							opendir (DIR,$jobDir);
							while (my $file=readdir(DIR)) {
								next if ($file eq '.' || $file eq '..' || $file eq 'resultsPep.txt');
								unlink "$jobDir/$file";
							}
							close DIR;
							sleep 10;
						}
					}
					else {
						@{$runningJobs{$currentJob}}=($jobDir,$pbsErrorFile,$jobClusterID);
						$launchedJobs{$currentJob}=1;
						last;
					}
				}
				
				my $numJobsLaunched=scalar keys %launchedJobs;
				open(FILESTAT,">>$fileStat");
				if ($numJobToRun >= $MAX_PARALLEL_QUANTIFS) {
					print FILESTAT "3/$numSteps Computing abundances ($numJobsLaunched/$numJobs tasks launched)\n";
				}
				else { # less than $MAX_PARALLEL_QUANTIFS jobs left: display number of remaining jobs
					print FILESTAT "3/$numSteps Computing abundances ($numJobToRun/$numJobs tasks remaining)\n";
				}
				close FILESTAT;
				
				last if $numJobsLaunched==$numJobs; # not necessary
				sleep 5;
			}
			
			###<Watching running jobs
			sleep 30;
			foreach my $job (sort{$a<=>$b} keys %runningJobs) {
				my ($jobDir,$pbsErrorFile,$jobClusterID)=@{$runningJobs{$job}};
				if (-e $pbsErrorFile) {
					my $pbsError=$cluster{'checkError'}->($pbsErrorFile);
					if ($pbsError) {
						if ($triedJobs{$job} == 3) {
							die "Job #$job ($jobClusterID) has failed"; #error will be handled by parent process (MASTER)
							#$numJobToRun--;
						}
						else { # reset job as not run yet
							#>clean jobDir except resultsPep.txt
							opendir(DIR,$jobDir);
							while (my $file = readdir(DIR)) {
								next if ($file eq '.' || $file eq '..' || $file eq 'resultsPep.txt');
								unlink	"$jobDir/$file";
							}
							close DIR;
							delete $runningJobs{$job};
							delete $launchedJobs{$job};
							next; # job
						}
					}
				}
				if (-e "$jobDir/end.flag") { # job has ended
					delete $runningJobs{$job};
					$numJobToRun--;
				}
				sleep 2;
			}
		}
	}
	else { # local job
		my $MAX_PARALLEL_QUANTIFS=&promsConfig::getMaxParallelJobs;
#$MAX_PARALLEL_QUANTIFS=4;		
		my %runningJobs;
		my $numJobToRun=$numJobs;
		$SIG{CHLD} = sub { # sub called when child communicates with parent (occur during parent sleep)
			local $!; # good practice. avoids changing errno.
			while (1) { # while because multiple children could finish at same time
				my $childPid = waitpid(-1,WNOHANG); # WNOHANG: parent doesn't hang during waitpid
				#my $childPid = waitpid(-1,POSIX->WNOHANG); # WNOHANG: parent doesn't hang during waitpid
				last unless ($childPid > 0); # No more to reap.
				delete $runningJobs{$childPid}; # untrack child.
				$numJobToRun--;
			}
		};
	
		####>Looping through job list (No restart on job error!)
		my $currentJob=0;
		MAIN_LOOP:while ($numJobToRun) {
	
			###>Parallel launch of up to $MAX_PARALLEL_QUANTIFS jobs
			my $newJobs=0;
			while (scalar (keys %runningJobs) < $MAX_PARALLEL_QUANTIFS) {
					
				##>Also check other user/launches jobs
				my $numProc=`ps -edf | grep runXICProtQuantification.pl | grep ABUND_JOB | grep -v grep | wc -l`;
				chomp($numProc);
				last if $numProc >= $MAX_PARALLEL_QUANTIFS;
	
				##>Ok to go
				$currentJob++;
				my $jobDir="$resultDir/ABUND$currentJob";
				if (-e "$jobDir/end.flag") { # in case using data from previous incomplete quantif
					$numJobToRun--;
				}
				else {
					system "echo LAUNCH > $jobDir/launch.flag"; # launch flag file
					my $childPid = fork;
					unless ($childPid) { # child here	
						system "./runXICProtQuantification.pl $quantifID $quantifDate ABUND_JOB:$currentJob";					
						exit;
					}
					$runningJobs{$childPid}=$currentJob;
				}
				last MAIN_LOOP if $currentJob==$numJobs; # no more jobs to launch => break MAIN_LOOP
				$newJobs++;
				sleep 2;
			}
			
			if ($newJobs) {
				open(FILESTAT,">>$fileStat");
				print FILESTAT "3/$numSteps Computing abundances ($currentJob/$numJobs tasks launched)\n";
				close FILESTAT;
			}
	
			###>Wait for potential ended job before next MAIN_LOOP
			sleep 60; # $SIG{CHLD} is active during sleep
	
		}
		
		###>Wait for last child processes to end
		while (scalar (keys %runningJobs) > 0) {
			my $jobsLeft=$numJobToRun;
			sleep 60; # $SIG{CHLD} is active during sleep
			if ($jobsLeft > $numJobToRun) { # job(s) have ended
				open(FILESTAT,">>$fileStat");
				print FILESTAT "3/$numSteps Computing abundances ($numJobToRun/$numJobs tasks remaining)\n";
				close FILESTAT;
			}
		}
	}
	
	exit;
}


sub jobProteinAbundance {
	my ($jobNumber)=@_;
	my $jobDir="$resultDir/ABUND$jobNumber";
	
	###<Running LFQ>###
	my $minPepRatioLFQ=($quantifParameters{'DB'}{'NUM_LFQ_RATIOS'})? $quantifParameters{'DB'}{'NUM_LFQ_RATIOS'}[0] : 1;
	system "$pathPython/python3 $promsPath{python_scripts}/computeLFQ.py -i $jobDir/resultsPep.txt -o $jobDir/resultsLFQ.txt -m $minPepRatioLFQ -f 2> $jobDir/commentsLFQ.txt";
	
	####<Parsing peptide results & computing peptide counts (NUM_PEP_TOTAL, NUM_PEP_USED & DIST_PEP_USED)>####
	my (%allProteins,%peptideValues,%proteinData); # Set <=> num label channels (1 for Label-Free)

	open(PEP,"$jobDir/resultsPep.txt");
	my %pepHeader2Idx;
	while(<PEP>) {
		chomp;
		my @lineData=split(/\t/,$_);
		####<Header------------------
		if ($.==1) { # 1st line of the file
			my $colIdx=0;
			foreach my $colName (@lineData) {
				$pepHeader2Idx{$colName}=$colIdx;
				$colIdx++;
			}
			next;
		}
		####<Data------------------
		my $protID0=$lineData[$pepHeader2Idx{'ProteinID'}];
		my $numPeptides=0;
		foreach my $pepID (split(/[\.\+\=]/,$lineData[$pepHeader2Idx{'PeptideId'}])) {
			$allProteins{$protID0}{$pepID}=1; # pepID recorded in case Labeling
			$numPeptides++;
		}
		next if $lineData[$pepHeader2Idx{'out'}] ne 'NA'; # skip outliers
		my $state=$lineData[$pepHeader2Idx{'Condition'}];
		my ($statePos)=$state=~/State(\d+)/;
		my $bioRep=$lineData[$pepHeader2Idx{'replicate'}];
		my $techRep=$lineData[$pepHeader2Idx{'repTech'}];
		my $seqVarMod=(split(/_/,$lineData[$pepHeader2Idx{'Peptide'}]))[0];
		$proteinData{'NUM_PEP_USED'}{$protID0}{$statePos}+=$numPeptides;
		$proteinData{'DIST_PEP_USED'}{$protID0}{$statePos}{$seqVarMod}=1; # DIST_PEP_USED

		####<Normalized peptide values
		push @{$peptideValues{$protID0}{$statePos}{$bioRep}{$techRep}},2**$lineData[$pepHeader2Idx{'log2Measure'}];
	}
	close PEP;
	
	####<Saving all proteins & number of peptides>####
	open(PROT,">$jobDir/allProteins.txt") || die $!; # contains list of all proteins (before outlier filter) to compute lost proteins later
	print PROT "ProteinID\tnumPep";
	foreach my $protID0 (keys %allProteins) {
		print PROT "\n",$protID0,"\t",scalar keys %{$allProteins{$protID0}};
	}
	close PROT;
	
	####<Computing protein abundance (mean, median & sum)>####
	my @quantTypes=('SUM_INT','MEAN_INT','MEDIAN_INT');
	foreach my $protID0 (keys %peptideValues) {
		foreach my $statePos (keys %{$peptideValues{$protID0}}) {
			my %bioRepAbundance;
			foreach my $bioRep (keys %{$peptideValues{$protID0}{$statePos}}) {
				my %techAbundance;
				foreach my $techRep (keys %{$peptideValues{$protID0}{$statePos}{$bioRep}}) {
					my $numValues=scalar @{$peptideValues{$protID0}{$statePos}{$bioRep}{$techRep}};
					my $midPos=int($numValues/2) || 1; # cannot be 0
					my $midPos1=$midPos+1;
					my $mustAverage=($numValues % 2)? 0 : 1;
					my $sumValue=0;
					my $count=0;
					my $medianValue=0;
					foreach my $pepValue (sort{$a<=>$b} @{$peptideValues{$protID0}{$statePos}{$bioRep}{$techRep}}) {
						$sumValue+=$pepValue;
						$count++;
						if ($count > $midPos1) {next;}
						elsif ($count==$midPos) {$medianValue=$pepValue;}
						elsif ($mustAverage && $count==$midPos1) {
							$medianValue=($medianValue+$pepValue)/2;
						}
					}
					%{$techAbundance{$techRep}}=(
						SUM_INT		=>	$sumValue,
						MEAN_INT	=>	$sumValue/$numValues,
						MEDIAN_INT	=>	$medianValue
					);
				}
				##>Take mean of tech reps
				my $numTechReps=scalar keys %techAbundance;
				#foreach my $quantType (@quantTypes) {$bioRepAbundance{$bioRep}{$quantType}=0;}
				foreach my $techRep (keys %techAbundance) {
					foreach my $quantType (@quantTypes) {$bioRepAbundance{$bioRep}{$quantType}+=$techAbundance{$techRep}{$quantType};}
				}
				foreach my $quantType (@quantTypes) {$bioRepAbundance{$bioRep}{$quantType}/=$numTechReps;}
			}
			##>Take mean of bio reps
			my $numBioReps=scalar keys %bioRepAbundance;
			#foreach my $quantType (@quantTypes) {$bioRepAbundance{$bioRep}{$quantType}=0;}
			foreach my $bioRep (keys %bioRepAbundance) {
				foreach my $quantType (@quantTypes) {$proteinData{$quantType}{$protID0}{$statePos}+=$bioRepAbundance{$bioRep}{$quantType};}
			}
			foreach my $quantType (@quantTypes) {$proteinData{$quantType}{$protID0}{$statePos}=int(0.5 + ($proteinData{$quantType}{$protID0}{$statePos}/$numBioReps));}
		}
	}
	
	####<Printing Abundance to file>####
	my %dataFiles=(	SUM_INT			=> 'resultsSUM.txt',
					MEAN_INT		=> 'resultsMEAN.txt',
					MEDIAN_INT		=> 'resultsMEDIAN.txt',
					NUM_PEP_USED 	=> 'peptidesUSED.txt',
					DIST_PEP_USED	=> 'peptidesDIST.txt',
					NUM_PEP_TOTAL	=> 'peptidesALL.txt',
					);
	my $numStates=scalar @{$quantifParameters{'DB'}{'STATES'}};
	my $headerStrg="ProteinID\tState".join("\tState",1..$numStates)."\n";
	foreach my $measType (keys %dataFiles) {
		open (OUT,">$jobDir/$dataFiles{$measType}") || die $!;
		print OUT $headerStrg;
		foreach my $protID0 (keys %{$proteinData{$measType}}) {
			print OUT $protID0;
			foreach my $statePos (1..$numStates) {
				print OUT "\t";
				if ($proteinData{$measType}{$protID0}{$statePos}) {
					print OUT ($measType eq 'DIST_PEP_USED')? scalar keys %{$proteinData{$measType}{$protID0}{$statePos}} : $proteinData{$measType}{$protID0}{$statePos};
				}
				else {
					print OUT 'NA';
				}
			}
			print OUT "\n";
		}
		close OUT;
	}
	
	exit;
}

sub importProteinAbundance {
	my ($numJobs)=@_;
	
	####<Parsing results>####
	open(FILESTAT,">>$fileStat");
	print FILESTAT "4/$numSteps Parsing results\n";
	close FILESTAT;

	my %dataFiles=(	MY_LFQ			=> 'resultsLFQ.txt',
					SUM_INT			=> 'resultsSUM.txt',
					MEAN_INT		=> 'resultsMEAN.txt',
					MEDIAN_INT		=> 'resultsMEDIAN.txt',
					NUM_PEP_USED 	=> 'peptidesUSED.txt',
					DIST_PEP_USED	=> 'peptidesDIST.txt'
					);
	my (%proteinAbundance,%allPeptides,%lostProteins);
	JOB:foreach my $job (1..$numJobs) {
		my $jobDir="$resultDir/ABUND$job";
		
		###<Combining data after outlier filtering>###
		foreach my $measType (keys %dataFiles) {
			unless (-e "$jobDir/$dataFiles{$measType}") {
				warn "File $dataFiles{$measType} not found for job #$job";
				last JOB;
			}
			open(RES,"$jobDir/$dataFiles{$measType}") || die $! ;
			my @statePos;
			while(<RES>) {
				chomp;
				my @lineData=split(/\t/,$_);
				####>Header------------------
				if ($.==1) { # 1st line of the file
					$statePos[0]=0; # not used
					foreach my $idx (1..$#lineData) {
						my ($statePos)=$lineData[$idx]=~/State(\d+)/;
						$statePos[$idx]=$statePos; # 1..n
					}
					next;
				}
				####>Data------------------
				my $protID0=$lineData[0];
				foreach my $idx (1..$#lineData) {
					next if ($lineData[$idx] eq 'NA' || abs($lineData[$idx]) < 1);
					$proteinAbundance{$protID0}{$statePos[$idx]}{$measType}=$lineData[$idx];
				}
			}
			close RES;
		}
		
		###<Combining data before outlier filtering & lost proteins>###
		open (PROT,"$jobDir/allProteins.txt");
		while (<PROT>) {
			next if $.==1;
			chomp;
			my ($protID0,$numPep)=split(/\t/,$_);
			$allPeptides{$protID0}=$numPep;
			unless ($proteinAbundance{$protID0}) {
				$lostProteins{$protID0}=1;
			}
		}
		close PROT;
	}	
	
	####<Recording data in DB>####
	my $dbh=&promsConfig::dbConnect('no_user');

	my ($quantifiedModifID,$quantifModStrg)=$dbh->selectrow_array("SELECT Q.ID_MODIFICATION,GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ',')
																					 FROM QUANTIFICATION Q
																					 LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
																					 WHERE Q.ID_QUANTIFICATION=$quantifID
																					 GROUP BY Q.ID_QUANTIFICATION");
	my $isModifQuantif=($quantifiedModifID || $quantifModStrg)? 1 : 0;
	
	####<Fetching list of quantification parameters>####
	my %quantifParamIDs;
	my $sthQP=$dbh->prepare("SELECT QP.ID_QUANTIF_PARAMETER,QP.CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD QM,QUANTIFICATION_PARAMETER QP WHERE Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND QM.ID_QUANTIFICATION_METHOD=QP.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=$quantifID");
	$sthQP->execute;
	while (my ($paramID,$code)=$sthQP->fetchrow_array) {
		$quantifParamIDs{$code}=$paramID;
	}
	$sthQP->finish;
	
	
	my $sthInsProt=$dbh->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_QUANTIFICATION,ID_PROTEIN,ID_QUANTIF_PARAMETER,QUANTIF_VALUE,TARGET_POS) VALUES ($quantifID,?,?,?,?)"); # PK autoincrement
	my ($sthInsModRes,$sthInsProtRes);
	if ($isModifQuantif) {
		$sthInsModRes=($quantifiedModifID)? $dbh->prepare("INSERT INTO MODIFIED_RESIDUE (ID_QUANTIFICATION,RESIDUE,POSITION) VALUES ($quantifID,?,?)")
										  : $dbh->prepare("INSERT INTO MODIFIED_RESIDUE (ID_QUANTIFICATION,RESIDUE,POSITION,MODIF_RANK) VALUES ($quantifID,?,?,?)");
		$sthInsProtRes=$dbh->prepare("INSERT INTO PROTQUANTIF_MODRES (ID_MODIF_RES,ID_PROT_QUANTIF) VALUES (?,?)");
	}

	foreach my $protID0 (keys %proteinAbundance) {
		my $protID=$protID0; # default
		my @allModResidues;
		if ($isModifQuantif) { # modif quantification
			next if $protID0 !~ /-/; # must be a mod prot
			($protID,my @modQuantifs)=split(/[-&]/,$protID0);
			foreach my $modQuantif (@modQuantifs) {
				my @modResidues;
				if ($quantifiedModifID) {@modResidues=split('\.',$modQuantif);}
				else {
					(my $modifRank,@modResidues)=split(/[#.]/,$modQuantif);
					foreach my $i (0..$#modResidues) {$modResidues[$i]=$modifRank.'#'.$modResidues[$i];} # Convert '1#V32.K44' into ['1#V32','1#K44']
				}
				&promsQuantif::insertModifiedResidues($protID,\@modResidues,\%modResiduesID,$dbh,$sthInsModRes); # <-- %modResiduesID should be full at <PROT1>
				push @allModResidues,@modResidues;
			}
		}
		foreach my $statePos (sort{$a<=>$b} keys %{$proteinAbundance{$protID0}}) {
			foreach my $paramCode (keys %{$proteinAbundance{$protID0}{$statePos}}) { # includes (NUM/DIST)_PEP_USED
				$sthInsProt->execute($protID,$quantifParamIDs{$paramCode},$proteinAbundance{$protID0}{$statePos}{$paramCode},$statePos);
				&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
			}
		}
		#>All peptides
		$sthInsProt->execute($protID,$quantifParamIDs{'NUM_PEP_TOTAL'},$allPeptides{$protID0},undef);
		&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
	}
	
	###>Lost proteins
	my %proteinAlias;
	if (scalar keys %lostProteins) {
		foreach my $protID0 (keys %lostProteins) {
			my ($protID)=($isModifQuantif)? $protID0=~/^(\d+)/ : ($protID0);
			$proteinAlias{$protID}='-'; # default
		}
		
		my $sthAlias=$dbh->prepare("SELECT ID_PROTEIN,ALIAS FROM PROTEIN WHERE ID_PROTEIN IN (".join(',',keys %lostProteins).")");
		$sthAlias->execute;
		while (my ($protID,$alias)=$sthAlias->fetchrow_array) {
			$proteinAlias{$protID}=$alias;
		}
		$sthAlias->finish;
	}
	
	####<Remove temp job directories>####
	#foreach my $job (1..$numJobs) {
	#	rmtree "$resultDir/ABUND$job";
	#}
	
	####<End quantification process>####
	my $projectID=&promsMod::getProjectID($dbh,$quantifID,'QUANTIFICATION');
	&endQuantification($dbh,$projectID,\%lostProteins,\%proteinAlias);
	
	exit;
}



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
my $context=($algoType eq 'PEP_RATIO')? $ratio : 'NA'; # Experiment != context (ratio can use multiple experiments)
			my ($numShared,$xicShared,$xicMinusInf,$xicPlusInf)=(0,0,0,0);
			my (%usedPepCode,%sharedDistinct);
			foreach my $pepCode (keys %{$refNumPepInCond->{$context}{$testCondID}}) {
				my ($pepSeq,$vMod,$charge,$dataSrc)=split(':',$pepCode);
				
				foreach my $feature (keys %{$refNumPepInCond->{$context}{$testCondID}{$pepCode}}) {
					if ($refNumPepInCond->{$context}{$refCondID} && $refNumPepInCond->{$context}{$refCondID}{$pepCode} && $refNumPepInCond->{$context}{$refCondID}{$pepCode}{$feature}) { # shared
						#$xicShared+=($refNumPepInCond->{$context}{$testCondID}{$pepCode} >= $refNumPepInCond->{$context}{$refCondID}{$pepCode})? $refNumPepInCond->{$context}{$testCondID}{$pepCode} : $refNumPepInCond->{$context}{$refCondID}{$pepCode};

						$xicShared+=($refNumPepInCond->{$context}{$testCondID}{$pepCode}{$feature}[0] + $refNumPepInCond->{$context}{$refCondID}{$pepCode}{$feature}[0]); # Sum of both channels!! (PP 27/10/15)
##$numShared++;
$numShared=1;
						$sharedDistinct{"$pepSeq:$vMod:$feature"}=1;
						$usedPepCode{"$pepCode:$feature"}=1;
					}
					else { # only in testCond
						$xicPlusInf+=$refNumPepInCond->{$context}{$testCondID}{$pepCode}{$feature}[0];
					}
				}
			}
			foreach my $pepCode (keys %{$refNumPepInCond->{$context}{$refCondID}}) {
				foreach my $feature (keys %{$refNumPepInCond->{$context}{$refCondID}{$pepCode}}) {
					next if $usedPepCode{"$pepCode:$feature"};
					$xicMinusInf+=$refNumPepInCond->{$context}{$refCondID}{$pepCode}{$feature}[0];
				}
			}
			##<Assigning ratio type
			if ($algoType=~/^(TDA|DIA)$/) { # MS2
##next if $numShared >= 6; # compute state mean with at least 6 frag
next if $numShared;
				next if scalar keys %sharedDistinct >= 6; # normal ratio <===================================================== OVERWRITE THRESHOLD
			}
			else { # DDA MS1
				next if ($labeling eq 'FREE' && $numShared); # compute state mean with any available peptide
				next if scalar keys %sharedDistinct >= 3; # normal ratio <===================================================== OVERWRITE THRESHOLD
			}
			
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


sub recreateModifSites { ## NEW (25/01/19)
	my ($pepSeq,$varModStrg,$pepBeg,$protLength,$refQuantifiedModifRes,$refRecreatedVarMods)=@_;
	#$varModStrg='' unless $varModStrg;
	$protLength=0 unless $protLength;
	unless ($refRecreatedVarMods->{"$pepSeq:$varModStrg"}) { # never seen before => compute new var mod
		#<Look for already modified res (cannot be used as target)
		my %modifiedResPos;
		if ($varModStrg) {
			foreach my $modStrg (split('&',$varModStrg)) {
				$modStrg=~s/^\d+://;
				foreach my $pos (split(/\./,$modStrg)) {
					$modifiedResPos{$pos}=1;
				}
			}
		}
		#<Loop through all modif to be recreated
		foreach my $modID (sort{$a<=>$b} keys %{$refQuantifiedModifRes}) {
			#<Look for targetable res in pepSeq
			my (@numResPos,@termResPos);
			foreach my $refRes (@{$refQuantifiedModifRes->{$modID}}) {
				my ($res,$context)=@{$refRes};
				if ($res=~/\w/) { # QUESTION: Should we also restrict to non-modified termini (res='=-*+')???????
					if ($context) {
						if ($context eq '=') { # Any N-term
							if ($pepSeq=~/^$res/ && !$modifiedResPos{1}) {
								push @termResPos,1;
								$modifiedResPos{1}=1;
							}
						}
						elsif ($context eq '-') { # Protein N-term
							if ($pepSeq=~/^$res/ && ($pepBeg==1 || $pepBeg==2) && !$modifiedResPos{1}) { # 1 or 2: N-term Methione 
								push @termResPos,$pepBeg;
								$modifiedResPos{1}=1;
							}
						}
						elsif ($context eq '*') { # Any C-term
							my $lastPos=length($pepSeq);
							if ($pepSeq=~/$res$/ && !$modifiedResPos{$lastPos}) {
								push @termResPos,$lastPos;
								$modifiedResPos{$lastPos}=1;
							}
						}
						else { # '+' Protein C-term
							my $lastPos=length($pepSeq);
							if ($pepBeg+$lastPos-1==$protLength && !$modifiedResPos{$lastPos} && $pepSeq=~/$res$/) {
								push @termResPos,$lastPos;
								$modifiedResPos{$lastPos}=1;
							}
						}
					}
					else { # Anywhere
						while ($pepSeq=~/$res/g) {
							my $pos=$+[0];
							unless ($modifiedResPos{$pos}) {
								push @numResPos,$pos;
								$modifiedResPos{$pos}=1;
							}
						}
					}
				}
				else { # N/Cterm
					unless ($modifiedResPos{$res}) {
						push @termResPos,$res;
						$modifiedResPos{$res}=1;
					}
				}
			}
			if ($termResPos[0] || $numResPos[0]) {
				$varModStrg.='&' if $varModStrg;
				$varModStrg.=$modID.':'.join('.',@termResPos);
				if ($numResPos[0]) {
					$varModStrg.='.' if $termResPos[0];
					$varModStrg.=join('.',sort{$a<=>$b} @numResPos);
				}
			}
		}
		$refRecreatedVarMods->{"$pepSeq:$varModStrg"}=$varModStrg;
	}
	return $refRecreatedVarMods->{"$pepSeq:$varModStrg"};
}

sub extractPtmPosProbability {
	my ($refModList,$varModCode,$varProbStrg)=@_;
	
	my $worsePtmPosProb=1;
	foreach my $modID (keys %{$refModList}) {
		next unless $varModCode =~ /(^|&)$modID:/; # this mod is not in varMod
		my ($a,$posStrg)=($varModCode =~ /(^|&)$modID:([^&]+)/);
	return -1 unless $posStrg; # in case missing MQ data
		my ($b,$probStrg)=($varProbStrg =~ /(^|&)$modID.*##PRB_([^&]+)/); # (...&)9%##PRB_MQ=1:0.998,3:0.966,10:0.036(&...)
		my $ptmPosProb=1;
		if ($probStrg) {
			if ($probStrg=~/=-1$/) { # no prob data available during search import
				$ptmPosProb=-1;
			}
			else {
				$ptmPosProb=1; # max by default
				foreach my $pos (split(/\./,$posStrg)) {
					my ($prob)=($varProbStrg =~ /\D$pos:([\d\.]+)/);
					$ptmPosProb=$prob if ($prob && $ptmPosProb > $prob); # record worse prob if multiple PTMs
				}
			}
		}
		else {$ptmPosProb=1;} # no prob data => 100% (no recording of 100% prob at import)
		$worsePtmPosProb=$ptmPosProb if $worsePtmPosProb > $ptmPosProb;
	}
	return $worsePtmPosProb;
}

sub checkForErrors {
	my ($dbh)=@_;
	#my ($errorFile) = glob "$promsPath{tmp}/quantification/current/*_$quantifDate\_error.txt";
	my $errorFile = "$promsPath{tmp}/quantification/current/$quantifID\_$quantifDate\_error.txt";
	if ($errorFile && -s $errorFile) {
#system "sed -i '/Permanently/ d' $errorFile"; # remove warning from ssh key addition
#		if (-s $errorFile) {
			$dbh->rollback;
			$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID"); # failed
			$dbh->commit;
			$dbh->disconnect;
			die "Aborting quantification due to errors.";
		#}
	}
}

sub generateReferenceProtFromQuantif { # globals: %promsPath,%quantifParameters,$dbh,$runDir,$projectID
	#my ($refQuantifID,$ratioPos,$refCondID,$testCondID)=&promsMod::cleanNumericalParameters(split(/[_:]/,$quantifParameters{'DB'}{'INTRA_PROT_REF'}[0]));
	my $refQuantifID=$quantifParameters{'DB'}{'INTRA_PROT_REF'}[0];
	my $numStatesInTest=scalar @{$quantifParameters{'DB'}{'STATES'}};

	my ($quantifAnnot)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$refQuantifID");
	
	my ($numStatesInRef,$normMethod);
	#my @statePos;
	foreach my $quantInfo (split('::',$quantifAnnot)) {
		if ($quantInfo=~/^NORMALIZATION_METHOD=/) {
			$normMethod=(split(/[=;]/,$quantInfo))[1];
			#<Convert old normalization names (myProMS v2)
			$normMethod=~s/\.normalization//;
			$normMethod='median.scale' if $normMethod eq 'global.mad';
		}
		elsif ($quantInfo=~/^STATES=/) {
			$quantInfo=~s/^STATES=//;
			$numStatesInRef=scalar(split(';',$quantInfo));
			#my $pos=0;
			#foreach my $state (split(';',$quantInfo)) {
			#	$pos++;
			#	if ($state=~/,#($refCondID|$testCondID)$/) {
			#		push @statePos,$pos;
			#		last if scalar @statePos==2;
			#	}
			#}
			last;
		}
	}
	
	##<Converting reference quantif table.txt in current quantif tableRef.txt>##
	#system "head -1 $promsPath{quantification}/project_$projectID/quanti_$refQuantifID/data/table.txt > $runDir/data/tableRef.txt";
	#system "grep -E '[[:space:]]State($statePos[0]|$statePos[1])[[:space:]]' $promsPath{quantification}/project_$projectID/quanti_$refQuantifID/data/table.txt >> $runDir/data/tableRef.txt";
	#if ($statePos[0] > 1) {
	#	system "sed -i s/state$statePos[0]/state1/ $runDir/data/tableRef.txt";
	#}
	#if ($statePos[1] != 2) {
	#	system "sed -i s/state$statePos[1]/state2/ $runDir/data/tableRef.txt";
	#}
	if ($numStatesInRef > $numStatesInTest) { # exclude unused states
		open(OUT,">$runDir/data/tableRef.txt");
		open(IN,"$promsPath{quantification}/project_$projectID/quanti_$refQuantifID/data/table.txt");
		while(<IN>) {
			if ($.==1) {
				print OUT $_;
				next;
			}
			$_=~/\tstate(\d+)\t/;
			print OUT $_ if $1 <= $numStatesInTest;
		}
		close IN;
		close OUT;
	}
	else {
		copy("$promsPath{quantification}/project_$projectID/quanti_$refQuantifID/data/table.txt","$runDir/data/tableRef.txt");
	}
	
	return $normMethod;
}


# TODO: Make clear choice for labeled quantif done with PEP_INTENSITY algo: treat as 100% Label-Free or mixed LF/Label ?????
# TODO: Move label-free peptide matching check further downstream for compatibility with PTM quantif
####>Revision history<####
# 2.14.0 [FEATURE] Job split and parallel launch for Abundance (PP 07/02/20)
# 2.13.1 [BUGFIX] in python command for LFQ (PP 28/01/20)
# 2.13.0 [FEATURE] Compatible with DIA and Protein Abundance by myProMS (PP 27/12/19)
# 2.12.0 [FEATURE] Compatible with manual peptide selection for normalisation: PEPTIDES_NORM & [BUGFIX] in peptide count for replicates (PP 03/12/19)
# 2.11.1 [ENHANCEMENT] Error check on number of proteins available (PP 22/11/19)
# 2.11.0 [FEATURE] Added optional co-exclusion of fully cut peptides included in missed-cut ones (PP 18/11/19) 
# 2.10.5 [BUGFIX] Fixed forgotten syntaxe error (PP 12/11/19)
# 2.10.4 [BUGFIX] Fixed bug in TDA-based quantif when same MS2 quantif is shared by multiple MS Analyses (PP 12/11/19)
# 2.10.3 [FEATURE] 3+ States allowed with protein-level normalization & uses MULTI_MODIF table for single-PTM quantif (PP 19/09/19)
# 2.10.2 [MODIF] Do not remove status log file to match new monitoring system behavior (VS 10/10/19)
# 2.10.1 [FEATURE] Fragments are summed into a single peptide intensity (PP 15/07/19)
# 2.10.0 [FEATURE] Adapted for MS2 fragments data using fragments as peptides (PP 18/06/19)
# 2.9.5 [Fix] bug where '=+' was used instead of '+=' for label-free peptides intensities (PP 18/06/19)
# 2.9.4 Works whether XXX_AREA or XXX_INTENSITY is used as quantification parameter (PP 22/03/19)
# 2.9.3 Uses MODIFICATION.SPECIFITY in case ANALYSIS_MODIFICATION.SPECIFICITY is empty (PP 25/02/19) 
# 2.9.2 Removed problematic filtering on peptide_quantificationX.txt extraction (PP 23/02/19)
# 2.9.1 Handles dynamic creation of PTM-sites (PP 19/02/19)
# 2.9.0 Also handles MaxQuant any PTM probability (PP 24/01/19)
# 2.8.0 Added protein-level normalization for site quantification (PP 19/12/18)
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
