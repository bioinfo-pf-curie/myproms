#!/usr/local/bin/perl -w

################################################################################
# runXICProtQuantification.pl       2.19.7                                     #
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
### myProMS code version for quantification: $MYPROMS_XIC_VERSION (based on R scripts & data preprocessing)
# 1: FC scripts (Ratio/TnPQ)
# 2: ML scripts
# 3: AS scripts before site-ratio correction by reference protein ratio (AnaDiff <= 4.0.13)
# 3.1: AS site-ratio correction & IB (AnaDiff >= 4.2.0)
# 3.2: Multi-modif quantif (format of modProtID is now <protID>-<modifRank1>#<pos1>.<pos2>&...)
# 3.3: TDA, DIA and Protein Abundance
# 3.4: shared prot/pep for normalization & multi-file normalization boxplots if too many samples
# 3.5: Protein-level normalization (Abundance only) / rejected proteins file [PP 18/01/21]
# 3.6: New infinite ratio decision rules (PP 11/04/21)

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

# exit;
###############################
####>Recovering parameters<####
###############################
my ($quantifID,$quantifDate,$runMode)=@ARGV; # only $quantifID & $quantifDate defined if Design-quantif
die "ERROR: Missing quantification id or job directory" if (!$quantifID || !$quantifDate);
$runMode='MASTER' unless $runMode; # MASTER, ABUND_DATA, ABUND_RSCRIPT:<context>, ABUND_SPLIT:<numAbundJobs>, ABUND_RUN:<numAbundJobs>, ABUND_JOB:<jobNum>, REF, NORM
# REF: To generate tableRef.txt file for protein-level correction of sites
# NORM: To generate temporary tableNorm.txt file from manually-selected list of peptides for site (global) normalization (TDA only)

#######################
####>Configuration<####
#######################
my $MYPROMS_XIC_VERSION='3.6';
my %promsPath=&promsConfig::getServerInfo('no_user');
my %cluster=&promsConfig::getClusterInfo; #('debian') # default is 'centos'
my $pathR=($cluster{'on'})? $cluster{'path'}{'R'} : $promsPath{'R'};
my $pathPython3=($cluster{'on'})? $cluster{'path'}{'python3'} : $promsPath{'python3'};
#my $minProtNumber=30; # used to desactivate FDR option
my $currentQuantifDir="$promsPath{tmp}/quantification/current";
my $errorFile="$currentQuantifDir/$quantifID\_$quantifDate\_error.txt";
my $quantifDir="$promsPath{tmp}/quantification/$quantifDate";
my $fileStat="$quantifDir/status_$quantifID.out"; # : "$quantifDir/status_$quantItemID.out";
my $runDir="$quantifDir/quanti_$quantifID";
my $dataDir="$runDir/data";
my $resultDir="$runDir/results";
my $graphDir="$resultDir/graph";
my $absDir = "$quantifDir/absolute";

my %quantifParameters=&promsQuantif::getQuantificationParameters("$quantifDir/quantif_info.txt");
my $algoType=$quantifParameters{'DB'}{'ALGO_TYPE'}[0]; # <=> {R}{design} for v3+ only! [PEP_RATIO|PEP_INTENSITY|TDA|DIA] or 'MSstats
my $software=($algoType eq 'MSstats')? 'MSstats' : 'myProMS';
my $ratioType=$quantifParameters{'DB'}{'RATIO_TYPE'}[0]; # SuperRatio or SimpleRatio or None (Abundance);
my $numStates=scalar @{$quantifParameters{'DB'}{'STATES'}};
my $countIdentifiedPep=($quantifParameters{'DB'}{'NUM_TRUE_USED'})? $quantifParameters{'DB'}{'NUM_TRUE_USED'}[0] : 0;

#>>> For Abundance only >>>
#my $normalizeProteins=($ratioType eq 'None')? 1 : 0; # TEMP <- Flag for protein-level normalization
my $normalizeIons=($quantifParameters{'DB'}{'BIAS_CORRECTION'}[0]=~/^TRUE,ion/)? 1 : 0; # Only for Abundance (ratioType='None')
my $normalizeProteins=($ratioType eq 'None' && $quantifParameters{'DB'}{'BIAS_CORRECTION'}[0]=~/^TRUE,.*prot$/)? 1 : 0; # Only for Abundance (ratioType='None')
my %optDataFiles=(	SUM_INT			=> 'resultsSUM.txt',
					MEAN_INT		=> 'resultsMEAN.txt',
					GEO_MEAN_INT	=> 'resultsGEO.txt',
					MEDIAN_INT		=> 'resultsMEDIAN.txt',
					MY_LFQ			=> 'resultsLFQ.txt',
					NUM_TRUE_USED	=> 'peptidesTRUE.txt'
				);
my %dataFiles=(	NUM_PEP_USED 	=> 'peptidesUSED.txt',
				DIST_PEP_USED	=> 'peptidesDIST.txt'
				); # NUM_PEP_TOTAL in allProteins.txt
$dataFiles{'NUM_TRUE_USED'}=$optDataFiles{'NUM_TRUE_USED'} if $countIdentifiedPep;
my %absMeasures;  # To store absolute quantification measures
#<<<

my $numSteps=($ratioType=~/S\w+Ratio/)? 3 : ($normalizeProteins)? 5 : 4;


#############################################################
#>>>>>>>>> PROCESSING PROTEIN ABUNDANCE METHOD <<<<<<<<<<<<<#
#############################################################
if ($runMode eq 'MASTER' && $ratioType eq 'None') { # 1st call for Abundance
	
	my $dbh=&promsConfig::dbConnect('no_user');
	$dbh->do("UPDATE QUANTIFICATION SET STATUS=0,QUANTIF_ANNOT=REPLACE(QUANTIF_ANNOT,'SOFTWARE=myProMS','SOFTWARE=myProMS;$MYPROMS_XIC_VERSION') WHERE ID_QUANTIFICATION=$quantifID");
	$dbh->commit;
	my ($quantifAnnot)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID"); # needed later if restoring
	$dbh->disconnect;
	
	#>Resuming from previous failed quantification
	if (-e $dataDir && -e "$dataDir/endData.flag") { # table.txt
		open(FILESTAT,">>$fileStat");
		print FILESTAT "Resuming previous quantification job\n";
		close FILESTAT;
	}

	#>Generating data file from peptide/fragment XICs
	if (!-e $dataDir || !-e "$dataDir/endData.flag") { # skip if usable data from previous incomplete quantif
		system "./runXICProtQuantification.pl $quantifID $quantifDate ABUND_DATA 2>> $errorFile"; # Default mode. Stops after R normalization results
		
		#>Check resulting data file
		my $okData=(-e "$dataDir/table.txt")? `head -3 $dataDir/table.txt | wc -l` : 0;
		$okData=~s/\D//g;
		warn "ERROR: Not enough entries in starting data file!" if $okData < 3;
		sleep 3;
		&checkForErrors;
		system "echo END > $dataDir/endData.flag"; # Data table generated without error

		# my $dataError=(-e "$dataDir/endData.flag")? 0 : 1;
		# unless ($dataError) {
		# 	my $okData=(-e "$dataDir/table.txt")? `head -3 $dataDir/table.txt | wc -l` : 0;
		# 	$okData=~s/\D//g;
		# 	$dataError=1 if $okData < 3;
		# }
		# if ($dataError) {
		# 	warn "ERROR: Not enough data in resulting file!";
		# 	$dbh=&promsConfig::dbConnect('no_user');
		# 	&checkForErrors($dbh);
		# }
	}
	
	#>Launch R normalization
	unless (-e "$resultDir/resultsPep.txt" && -e "$runDir/endR.flag") { # in case using data from previous incomplete quantif
		system "./runXICProtQuantification.pl $quantifID $quantifDate ABUND_RSCRIPT:PRIMARY 2>> $errorFile"; # Default mode. Stops after R normalization results
		# or simply &launchRscript('PRIMARY'); # Runs in master process
	}

	#>Splitting resultsPep.txt for parallel jobs
	my $computeLFQ = 0;
	my $computeAbsolute = 0;
	my $computeOther = 0;  # SUM_INT, MEAN_INT, GEO_MEAN_INT, MEDIAN_INT
	my $computeNonAbsolute = 0;  # LFQ or other abundance measure that does not compute amounts
	my $numAbundJobs = 0;
	foreach my $meas (@{$quantifParameters{'DB'}{'ABUND_MEASURES'}}) {
		last if ($computeLFQ && $computeAbsolute && $computeOther);		
		if ($meas =~ /^MY_LFQ/) {
			$computeLFQ=1;
		} elsif ($meas =~ /^aLFQ/) {
			$computeAbsolute = 1;
		} else {
			$computeOther = 1;
		}
	}
	$computeNonAbsolute = ($computeLFQ || $computeOther);

	if ($computeAbsolute) {
		system "./runXICProtQuantification.pl $quantifID $quantifDate ABSOLUTE 2>> $errorFile"; # Compute absolute quantification with aLFQ R package
		&checkRsuccess("$absDir/aLFQ.Rout");  # Check that abs quantif computation has run correctly
	}

	if ($computeNonAbsolute) {  # No need to compute this if only absolute quantif was required

		$numAbundJobs = 1;
		if ($computeLFQ) {
			if (-e "$resultDir/endSplit.flag") { # split already done
				$numAbundJobs=`head -1 $resultDir/endSplit.flag`;
				$numAbundJobs=~s/\D//g;
			}
			else { # run split
				$quantifAnnot=~s/::/§/g;
				my ($stateStrg)=$quantifAnnot=~/STATES=([^§]+)/;
				my $numAllRep=scalar(split(/[;.&]/,$stateStrg));
				my $numLines=`wc -l $resultDir/resultsPep.txt`;
				$numLines=~s/^(\d+).*/$1/;
				my $numJobLn=int(1.5 + $numLines / 2.5E5);
				my $numJobSt=$numAllRep;
				$numAbundJobs=($numJobLn >= $numJobSt)? $numJobLn : $numJobSt;
				$numAbundJobs=5 if $numAbundJobs < 5;
			}
		}
		#>Launching data-splitting job
		if (!-e "$resultDir/endSplit.flag") { # skip if usable data from previous incomplete quantif
			#>Launching resultsPep splitting
			system "./runXICProtQuantification.pl $quantifID $quantifDate ABUND_SPLIT:$numAbundJobs 2>> $errorFile";
			&checkForErrors;
			system "echo $numAbundJobs > $resultDir/endSplit.flag"; # ends w/o error
		}
		#>Launching parallel abundance computation main job
		if (!-e "$resultDir/endAbund.flag") { # skip if usable data from previous incomplete quantif
			system "./runXICProtQuantification.pl $quantifID $quantifDate ABUND_RUN:$numAbundJobs 2>> $errorFile";
			&checkForErrors;
			system "echo END > $resultDir/endAbund.flag"; # ends w/o error
		}
	}
	else {  # Need to compute protein data for absolute if not done with other abundance measures
		&getAllProtPepNb($absDir);
	}

	my (%proteinAbundance,%allPeptides,%lostProteins); # records final data to import in DB
	my $stepNumber=4;
	
	#>Protein-level normalization
	if ($normalizeProteins) {
		#>Parse primary results
		my %primaryAbundance;
		&parsePrimaryAbundanceResults($numAbundJobs,\%primaryAbundance,\%allPeptides,\%lostProteins,\%proteinAbundance); # %proteinAbundance will only collect peptide counts at this stage
		
		unless (-e "$resultDir/end2ndNorm.flag") { # resumed Protein normalization generated without error
			&runSecondaryNormalizations($numAbundJobs,\%primaryAbundance,\%proteinAbundance,\%allPeptides,\%lostProteins);
		}
		$stepNumber++;
	}

	#>Processing final results (primary/2ndary normalization)
	&processFinalAbundanceResults($stepNumber,$numAbundJobs, $computeAbsolute, $computeNonAbsolute, \%proteinAbundance,\%allPeptides,\%lostProteins);

	#>Importing results in SQLite database
	&importAbundanceInDB($stepNumber,$numAbundJobs, $computeAbsolute, \%proteinAbundance,\%allPeptides,\%lostProteins); # DB import

	exit;
}
elsif ($runMode=~/ABUND_RSCRIPT:(.+)/) {
	&launchRscript($1);
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
elsif ($runMode=~/ABSOLUTE/) {
	&runAbsoluteQuantif();
	exit;
}


##################################################################################
#>>>>>>>>>> MASTER/REF/NORM (diffAna) or ABUND_DATA (normalization only) >>>>>>>>#
##################################################################################
my $topN=($quantifParameters{'DB'}{'TOP_N'})? $quantifParameters{'DB'}{'TOP_N'}[0] : 0; # number peptide used for PEP_INTENSITY (except MSstats)
my $matchingPep=($quantifParameters{'DB'}{'MATCH_PEP'})? $quantifParameters{'DB'}{'MATCH_PEP'}[0] : 0; # label-free SimpleRatio: peptides must be the same across conditions
my $referenceMode=($runMode =~/^(REF|NORM)$/)? $runMode : ''; # for table(Ref/Norm).txt (not for MSstats)
my $referenceStrg=($referenceMode eq 'REF')? ' reference' : ($referenceMode eq 'NORM')? ' normalization' : '';
# open (DEBUG,">$promsPath{tmp}/quantification/debug.txt") if !$referenceMode; # DEBUG!!!!
my $dbh=&promsConfig::dbConnect('no_user');

if (!$referenceMode) {
	if ($ratioType=~/S\w+Ratio/) {
		$dbh->do("UPDATE QUANTIFICATION SET STATUS=0,QUANTIF_ANNOT=REPLACE(QUANTIF_ANNOT,'SOFTWARE=myProMS','SOFTWARE=myProMS;$MYPROMS_XIC_VERSION') WHERE ID_QUANTIFICATION=$quantifID");
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
my ($quantiAnnotation,$quantifiedModifID,$quantifModStrg)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,Q.ID_MODIFICATION,GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ','),Q.ID_QUANTIFICATION
																					 FROM QUANTIFICATION Q
																					 LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
																					 WHERE Q.ID_QUANTIFICATION=$quantifID
																					 GROUP BY Q.ID_QUANTIFICATION");
my ($isModifQuantif,$trueModifQuantif,$isFreeResQuantif)=(0,0,0);
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
		if ($modID > 0) {$trueModifQuantif=1;} else {$isFreeResQuantif=1;} # -1 for fake "FreeResidue"
	}
	$quantifiedModifID=$modifQuantifs[0] if $modifRank==1; # single PTM: back compatibility in case QUANTIFICATION.ID_MODIFICATION no longer used
	$matchQuantifModStrg=join('|',@modifQuantifs); # Needed later for regexp on pep varmod
}
my $keepSeqContext=($isModifQuantif && $quantifParameters{'DB'}{'KEEP_SEQ_CONTEXT'} && $quantifParameters{'DB'}{'KEEP_SEQ_CONTEXT'}[0])? 1 : 0; #keep sequence context of PTM (=> Do not combine site from fully cleaved & missed cleaved peptides)

###>Protein enriched by PTMs?
my $isProteinByPTM=0;
my (%proteinPTMcontexts,%peptideContextMatch,%proteinSequence); # also used if FREE_RESIDUE
if ($quantifParameters{'DB'}{'PROTEIN_PTM'}) { # Restrict peptides to specific PTMs & context
	$isProteinByPTM=1;
	&promsQuantif::preparePTMcontexts(\%proteinPTMcontexts,$quantifParameters{'DB'}{'PROTEIN_PTM'});
}

my ($labeling)=($quantiAnnotation=~/LABEL=([^:]+)/);
$labeling=uc($labeling);

my %skipFromRefModIDs;
if ($referenceMode) { # Switching to non-modif quantif
	%skipFromRefModIDs=%modifQuantifRank;
	$quantifiedModifID=0;
	$isModifQuantif=0;
}
# my $quantifiedModifRemoved=0; # to handle site for which targeted PTM was removed by experimental protocole
# if ($isModifQuantif && $quantifParameters{'DB'}{'CREATE_PTM'}) {
# 	$quantifiedModifRemoved=1;
# }

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
			system "./runXICProtQuantification.pl $quantifID $quantifDate REF 2>> $errorFile"; # referenceMode set to REF to create tableRef.txt
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
		system "./runXICProtQuantification.pl $quantifID $quantifDate NORM 2>> $errorFile"; # referenceMode set to NORM to create tableNorm.txt
	}
}

####>Main data file
#my $dataFile=($referenceMode eq 'REF')? "$dataDir/tableRef.txt" : ($referenceMode eq 'NORM')? "$dataDir/tableNorm.txt" : "$dataDir/table.txt";
#open(DATA,">$dataFile"); # Name of data table
##if ($ratioType eq 'Ratio') { # Old algos (FC)
##	print DATA "Protein_ID\tPep_Seq\tVar_Mod\tCharge\tPeptide_IDs\tProtein_Name\tProtein_Validity";
##}
##else { # Algos v2+
#	print DATA "Protein_ID\tPeptide\tSample\tTreatment\tReplicate\tTechnical_Replicate\tExperiment\tQuantif_Set\tPeptide_IDs\tProtein_Name\tProtein_Validity\tValue" unless $referenceMode eq 'NORM';
##}
my %statesInFile;

################################
####>Get states composition<####
################################
my @quantiOrder;# In order to print the table in a good way for R
my %observationInfo; # for (Super/Simple)Ratio only
my (%poolGroup2Cond,%ana2Obs,%cond2Obs,%obs2Cond,%labelModifs,%anaLabelModifs,%xicEngine); #,%obsID2Cond,%ana2Cond
my (%ana2Experiment,%anaConds,%anaState1RatioPos); #<----- NEW EXPERIMENT assignment 31/10/17 -----
#my (%ana2Quant,%stateToCond,%conditionList,%condPeptideList); # for MSstats

my $sthPQN=$dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?"); # peptide quantifs

my @ratioTarget;
my %replicateTargetPos; my $numRepTgtPos=0;
#my ($hasMultiBioReps,$hasMultiTechReps)=(0,0);
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
		my $statePos=0;
		if ($annoStrg =~ /^STATES=(.*)/) { # 'state1;state2;---;stateN',  stateX='nbBioRep,bioRep1.bioRep2.---.bioRepN,#condID', bioRepX='techRep1=techRep2=---=techRepN', techRepX='frac1+frac2+---+fracN'
			##my (%replicAnaStrg,%anaExpCode,%commonRefCondAna); # ,$multiTechRep for (Super/Simple)Ratio only
			##my $lastReplicInRefIdx=-1; # for SuperRatio only: last index of all bio/tech replicates for reference condition in @quantiOrder
#----- NEW EXPERIMENT assignment 31/10/17 ----->
my (@anaOrder,%anaTargetPosCond); # to distinguish different Experiments (labeled ana only)
			my $MSstatsRun=0; # for MSstats
			foreach my $state (split(/;/,$1)) {
				$poolObservations=1 if $state=~/\+/;
				# if ($poolObservations && $algoType=~/^(TDA|DIA|MSstats)$/ ) {
				# 	warn "ERROR: Fractions are not compatible with TDA/DIA quantification. Correct your design.";
				# 	sleep 5;
				# 	&checkForErrors($dbh);
				# }
				$statePos++;
				
				### Tests fractions if Label-Free ###
				###if ($labeling eq 'FREE' && $poolObservations) {
				###	close DATA;
				###	$sthPQN->finish;
				###	$sthALM->finish;
				###	$sthObsMod->finish;
				###	$dbh->disconnect;
				###	die "ERROR in State #$statePos: Fractions detected for Label-Free quantification. Correct your design.";
				###}
#----- EXPERIMENT assignment #OLD ----->
				##my $expCode=chr(63+$statePos); # SuperRatio only 2 -> A, 3 -> B, ... (not used for ref cond=1)
				#$experimentCode{$expCode}=1 if $statePos >= 2; # skip Ref (1st) Cond
				$state=~s/#//g; # remove all id tags
				my ($numBioRep,$quantiObsIDs,$condID)=split(/,/,$state);
				#if ($ratioType eq 'Ratio') {
				#	push @{$quantifParameters{'R'}{'grp'}},$numBioRep; # Number of replicates of the condition
				#	push @{$quantifParameters{'R'}{'name.grp'}},"State$statePos"; # $condID
				#}
#elsif ($ratioType eq 'SuperRatio') { # to build primary ratio
	push @conditionList,$condID;
#}
				#elsif ($statePos==1) {$lastReplicInRefIdx++;} # for SuperRatio only
				my $bioRepCount=0;
				my $allRepCount=0;
				my $multiRep=($quantiObsIDs=~/[.&]/)? 1 : 0; # for Ratio only
				#my %stateAna; # for SuperRatio only
				foreach my $bioReplicate (split(/\./,$quantiObsIDs)) {
					$bioRepCount++;
#$hasMultiBioReps=1 if $bioRepCount > 1;
					my $techRepCount=0;
					foreach my $techReplicate (split(/&/,$bioReplicate)) {
						##$lastReplicInRefIdx++ if $statePos==1; # for SuperRatio only
						$allRepCount++;
						$techRepCount++;
						$MSstatsRun++;
#$hasMultiTechReps=1 if $techRepCount > 1;

##>REPLICATE DEFINITION: BIOREP or TECHREP BASED on RESIDUAL_VAR!!!
$numRepTgtPos++ if ($techRepCount==1 || ($quantifParameters{'DB'}{'RESIDUAL_VAR'} && $quantifParameters{'DB'}{'RESIDUAL_VAR'}[0] eq 'technical')); # count either bioRep or techRep
$replicateTargetPos{"State$statePos.rep$bioRepCount.techRep$techRepCount"}=$numRepTgtPos;
						###$multiTechRep=1 if $techRepCount>1; # records if exists multiple tech rep
						push @quantiOrder,$techReplicate;
						$obs2Cond{$techReplicate} = $condID;

						#>Required to detect +/- infinite ratios (except %labelModifs)
						$poolGroup2Cond{$techReplicate}=$condID;
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
							if ($labeling eq 'FREE') {  # Need labeling eq FREE and not algoType ne PEP_RATIO 17/09/20
								push @{$ana2Obs{$anaID}},$obsID;
							}
							else {
								$anaTargetPosCond{$anaID}{$targetPos}=$condID;
								push @anaOrder,$anaID;
								$anaConds{$anaID}{$condID}=1; # needed later for infinite ratio check
								if ($labeling eq 'SILAC') { # SILAC QUANTI
									##if ($ratioType eq 'SuperRatio') {
									##	if ($statePos==1) {push @{$commonRefCondAna{$techReplicate}},$anaID;} # Common ref condition
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
####								#@{$observationInfo{$techReplicate}}=("State$statePos",'rep1',"techRep$techRepCount",chr(64+$bioRepCount));
####@{$observationInfo{$techReplicate}}=("State$statePos","rep$bioRepCount","techRep$techRepCount",'A'); #chr(64+$bioRepCount)
####								#$designCheck{chr(64+$bioRepCount)}{$statePos}=1;
####							}
####							else {
####								@{$observationInfo{$techReplicate}}=("State$statePos","rep$bioRepCount","techRep$techRepCount");
####								if ($ratioType eq 'SimpleRatio') {
####									push @{$observationInfo{$techReplicate}},'A'; ###'NA';
####									#$designCheck{'NA'}{$statePos}=1;
####								} # Experiment is irrelevent in table.txt but needed later
####								elsif ($statePos > 1) { # cond #1 has multiple exeriments
####									push @{$observationInfo{$techReplicate}},$expCode;
####									#$designCheck{$expCode}{$statePos}=1;
####								}
####							}
####						}
						if ($ratioType=~/S\w+Ratio|None/) {
							@{$observationInfo{$techReplicate}}=($software eq 'MSstats')? ("State$statePos","State$statePos\_BioRep$bioRepCount",$MSstatsRun) : ("State$statePos","rep$bioRepCount","techRep$techRepCount");
#----- EXPERIMENT assignment #OLD ----->
							##if ($Rdesign eq 'SUPERRATIO') {
							##	push @{$observationInfo{$techReplicate}},$expCode if $statePos > 1; # cond #1 has multiple experiments
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
				##	if ($statePos > 1) { # not for common Ref ($statePos=1)
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
			##		#$designCheck{$expCode}{$statePos}=1;
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
if ($algoType eq 'PEP_RATIO') { # former $labeling ne 'FREE' 27/04/20
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
my $usedTargetPos=scalar @ratioTarget; # ($quantifParameters{'DB'}{'SINGLE_REF'})? $numStates-1 : ($numStates * ($numStates-1)) / 2;
$usedTargetPos+=$numStates if $algoType ne 'PEP_RATIO'; # former =~/PEP_INTENSITY|TDA|DIA/ 27/04/20 # num states (MEAN_INTENSITY)
foreach my $techReplicate (keys %replicateTargetPos) {
	$replicateTargetPos{$techReplicate}+=$usedTargetPos; # shift to free target pos
}

$sthPQN->finish;
$sthALM->finish;

###>Handling Modification Quantif phosphoRS/MaxQuant threshold & PTM position ambiguity settings
my ($ptmProbSoftCode,$ptmProbThreshold,$ambiguousModifPos,%quantifResMatchStrg,%qQuantifResMatchStrg)=('',0,0,'',''); # default ($qQuantifResMatchStrg PP 2017/06/22)
# my (%targetableRes,$matchProtCterm); # for free residues PTMs only
if ($isModifQuantif) {
	if ($isFreeResQuantif) { # assumes single-modif quantif!!!!!!!
		# my $context=''; # for now
		# my $modID=$modifQuantifs[0]; # -1 (fake "freeResidus")
		# foreach my $res (@{$quantifParameters{'DB'}{'FREE_RESIDUES'}}) {
		# 	push @{$targetableRes{$modID}},[$res,$context];
		# 	$quantifResMatchStrg{$modID}.=$res;
		# 	$matchProtCterm=1 if ($res eq '+' || $context eq '+'); # for at least 1 modifID
		# }
		&promsQuantif::preparePTMcontexts(\%proteinPTMcontexts,$quantifParameters{'DB'}{'SITE_CONTEXT'});
	}
	if ($trueModifQuantif) { # true PTMs to quantify
		if ($quantifParameters{'DB'}{'PTM_POS'}[0]=~/(\w+):(\d+)/) { # PRS Phospho or MQ/SPC/PTMRS any PTM
			$ptmProbSoftCode=$1; # PRS|MQ|SPC|PTMRS
			$ptmProbThreshold=$2; $ptmProbThreshold/=100 if $ptmProbSoftCode =~ /^(MQ|SPC|PTMRS)$/; # 0-1
			$ambiguousModifPos=1 if $quantifParameters{'DB'}{'PTM_POS'}[1] eq 'ambiguous';
		}
		elsif ($quantifParameters{'DB'}{'PTM_POS'}[0] eq 'ambiguous') { # other PTMs valid/ambiguous
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

###>Contaminants
my %contaminants;
if ($quantifParameters{DB}{EXCL_CONTAMINANTS} && $quantifParameters{DB}{EXCL_CONTAMINANTS}[0]) {
	my $sthCont=$dbh->prepare("SELECT DISTINCT ID_PROTEIN FROM ANALYSIS_PROTEIN AP
								INNER JOIN ANALYSIS_DATABANK AD ON AP.ID_ANALYSIS=AD.ID_ANALYSIS
								INNER JOIN DATABANK D ON AD.ID_DATABANK=D.ID_DATABANK
								WHERE AP.ID_ANALYSIS=? AND D.IS_CRAP=1 AND AP.DB_RANK LIKE CONCAT('%',AD.DB_RANK,'%')"); # LIKE because 2 can be in 2 or 12!!!
	foreach my $anaID (keys %ana2Obs) {
		$sthCont->execute($anaID);
		while (my ($protID)=$sthCont->fetchrow_array) {$contaminants{$protID}=1;}
	}
	$sthCont->finish;
}

###>Reference proteins for bias correction
#my ($sharedProtNorm,$normProtUsage)=(0,'use');
my $normProtUsage=''; # default: Normalization uses all proteins 
my (%forNormalization,%notForNormalization);
if ($quantifParameters{DB}{BIAS_CORRECTION}[1] && !$referenceMode) { # Correction with custom list. done only once in normal mode
	#(my $listID,$normProtUsage)=@{$quantifParameters{DB}{BIAS_CORRECTION}}[1,2];
	#$listID=~s/#//; # remove id flag
	#if ($listID  eq 'shared') {
	#	$sharedProtNorm=1;
	#	$normProtUsage='use';
	#}
	#elsif ($listID  eq 'sharedPep') {
	#	$sharedProtNorm=0;
	#	$normProtUsage=''; # shared item computation: done later by R
	#}
	if ($quantifParameters{DB}{BIAS_CORRECTION}[1]=~/^#*(\d+)$/) { # Custom list
		my $listID=$1;
		$normProtUsage=$quantifParameters{DB}{BIAS_CORRECTION}[2] || 'use';
		$sthList->execute($listID);
		my $refList=($normProtUsage eq 'use')? \%forNormalization : \%notForNormalization;
		while (my ($protID)=$sthList->fetchrow_array) {
			$refList->{$protID}=1;
		}
	}
	# Shared peptides handled by R
}

$sthList->finish;


###> Write the table with all results
##################################################################
####>Query to get all the validated proteins and the peptides<####
##################################################################
my ($pepSpecificity,$pepMissedCleavage,$ptmFilter,$pepChargeState,$pepSource)=@{$quantifParameters{'DB'}{'PEPTIDES'}};
my $MGSharedPep=($quantifParameters{'DB'}{'MG_SHARED_PEP'})? $quantifParameters{'DB'}{'MG_SHARED_PEP'}[0] : 'best'; # best|share|exclude. share only for Abundance (PP 20/01/21)
my $manualPepSelection=0; # used in reference mode (REF|NORM) only
my (%manualPeptides,%manualProteins); # used in reference mode (REF|NORM) only
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
elsif ($quantifParameters{'DB'}{'PEPTIDES_NORM'} && $quantifParameters{'DB'}{'PEPTIDES_NORM'}[0] eq 'manual') { # Bias correction with manually selected peptides
	if ($referenceMode eq 'NORM') {$manualPepSelection=1;}
	else {$normProtUsage='use';} # normal mode. needed to print normalization file
	foreach my $i (1..$#{$quantifParameters{'DB'}{'PEPTIDES_NORM'}}) {
		my ($protID,@pepIDs)=split(/[:,]/,$quantifParameters{'DB'}{'PEPTIDES_NORM'}[$i]);
		if ($referenceMode eq 'NORM') {
			foreach my $pepID (@pepIDs) {$manualPeptides{$pepID}=1;}
#print DEBUG "MANUAL PEP=@pepIDs\n";
			$manualProteins{$protID}=1;
		}
		else {$forNormalization{$protID}=1;} # normal mode
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

my $ptmProbQueryStrg=($ptmProbSoftCode && $ptmProbSoftCode ne 'PRS')? "GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,'%',COALESCE(PM.REF_POS_STRING,'') ORDER BY PM.ID_MODIFICATION SEPARATOR '&')" : "'-'";
my $quantiQuery=qq
|SELECT GROUP_CONCAT(DISTINCT PPA.ID_PROTEIN),P.ID_PEPTIDE,PEP_BEG,PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),$ptmProbQueryStrg,CHARGE,DATA,P.SCORE,PPA.IS_SPECIFIC,AP.PEP_SPECIFICITY,MISS_CUT
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
$quantiQuery.=' GROUP BY PPA.ID_PROTEIN,P.ID_PEPTIDE,PEP_BEG ORDER BY ABS(AP.VISIBILITY) DESC,AP.PEP_SPECIFICITY DESC,AP.NUM_PEP DESC,ABS(PEP_BEG) ASC';


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

my (%proteins,%intensity,%proteinAlias,%excludedSeq,%peptideBeg,%qSetSequence); #,%proteinInConds %retentionTime,%bestVisibility,%visInReplicate,%visInAnalysis,%visibility
my (%qSetBestVarMod,%qSetIsBad); # SILAC only
my (%tempIntensity,%dataSource2ObsSet,%qSet2Protein,%pepSeq2qSet,%qSet2AllProteins); # for SILAC only
my %notUsed4Quantif; # for future $normalizeWithAll management (Not yet used in table.txt)
my %pepID2QuantifSetSILAC; # PP 25/10/17: used to better compute %usedPeptideSets for SILAC when qSets are summed
my %isGhostPeptide; # list of ghosts for recording num true peptides/state. Needed at DB import step (global during table.txt generation & ratios import, local otherwise!)

##----- NEW EXPERIMENT assignment ----->
## RULES: 1 techRep <-> 1 experiment => 1 analysis <-> 1 experiment => 1 pepID <-> 1 experiment
## Should handle:
##    -fractions
##    -bioRep across different Experiments (LABELED with 3+ channels)
##    -SuperRatio: state1 found in multiple Experiments
##my (%analysis2Experiment,%techRep2Experiment); # %peptideID2Experiment,
##my $currentExperiment='A';

my $entityStrg=($algoType=~/^(TDA|DIA|MSstats)$/)? 'fragment' : 'peptide';
if (!$referenceMode || $referenceMode eq 'REF') {
	open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
	print FILESTAT "1/$numSteps Generating data files (1/4 Retrieving$referenceStrg $entityStrg intensity data)\n";
	close FILESTAT;
}
my %pepQuantifValues;
my $numAllObs=0;
if ($algoType=~/^(TDA|DIA|MSstats)$/) { # MS2 XIC
	my %fragQuantifs;
	foreach my $observationSet (@quantiOrder) {
		foreach my $fraction (split(/\+/,$observationSet)) {
			my ($obsID,$ms2Quantif,$anaID,$targetPos)=split(/:/,$fraction);
			$fragQuantifs{$ms2Quantif}{$anaID}=1;
			$numAllObs++;
		}
	}
	my $numObsData=0;
	QUANTI:foreach my $ms2Quantif (keys %fragQuantifs) {
		foreach my $anaID (keys %{$fragQuantifs{$ms2Quantif}}) {
			unless (-e "$promsPath{quantification}/project_$projectID/quanti_$ms2Quantif/swath_ana_$anaID.txt") {
				warn "Missing fragment quantification file (#$ms2Quantif, Ana #$anaID)"; # warn to allow clean handling of exception by &checkForErrors
				last QUANTI;
			}
			open(FRAG_QUANTIF,"$promsPath{quantification}/project_$projectID/quanti_$ms2Quantif/swath_ana_$anaID.txt") || die $!;
			while (<FRAG_QUANTIF>) {
				next if $.==1; # /^PEP/; #
				s/\s\Z//g; # chomp is not enough. More hidden Windows character
				my ($pepID,$fragMZ,$fragCharge,$fragType,$fragRes,$fragArea,$fragRT)=split(/!/,$_);
				#next if $fragArea=~/^(N\/*A|None)$/;
				next if $fragArea=~/[^\d\.]/;
				next if $fragArea <= 0; # can be '0.0'
				#next if (!$peptideList{$pepID});
				$pepQuantifValues{"$pepID:0"}{"$fragType$fragRes#$fragCharge"}=$fragArea;
			}
			close FRAG_QUANTIF;
			$numObsData++;
			unless ($numObsData % 10) {
				open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
				print FILESTAT "1/$numSteps Generating data files (1/4 Retrieving$referenceStrg $entityStrg intensity data: $numObsData/$numAllObs)\n";
				close FILESTAT;
			}
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
			$numAllObs++;
		}
	}
	my $numObsData=0;
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
				next if $quantifValue=~/[^\d\.]/;
				next if $quantifValue <=0; # Isobaric bug
				unless ($signalParamID) {
					next unless $usedParamIDs{$xicQuantif}{$paramID};
					$signalParamID=$paramID;
				}
				next if $paramID != $signalParamID;
				$pepQuantifValues{"$pepID:$targetPos"}{PEP}=$quantifValue;
			}
			close PEP_QUANTIF;
			$numObsData++;
			unless ($numObsData % 10) {
				open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
				print FILESTAT "1/$numSteps Generating data files (1/4 Retrieving$referenceStrg $entityStrg intensity data: $numObsData/$numAllObs)\n";
				close FILESTAT;
			}
		}
	}
}
&checkForErrors($dbh);

open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
print FILESTAT "1/$numSteps Generating data files (2/4 Retrieving$referenceStrg identification data)\n";
close FILESTAT;
my %recreatedVarMods; # to prevent re-creating same varMod multiple times
#my %proteinLength; # Only if sites to be recreated at protein C-term
my (%missedCutPeptides,%beg2pepSeq); # only if $pepMissedCleavage=-1 (exclude missCut AND overlapping cut peptides)
my (%fragmentsInPepIon);#,%peptidesUsedMSstats); # MSstats
my %peptideSpecificity; # to attribute a shared peptide to the "best" proteins
my %bestProteotypicity; # to attribute a shared peptide to the "best" proteins
my %ptmProbability; # for other than PEP_RATIO+SILAC: if $ambiguousModifPos to decide if ion should be ambiguous or not
# my $sthProtL=$dbh->prepare('SELECT PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=?');
my $sthProtSeq=$dbh->prepare('SELECT P.PROT_SEQ,MP.PROT_SEQ FROM PROTEIN P LEFT JOIN MASTER_PROTEIN MP ON P.ID_MASTER_PROTEIN=MP.ID_MASTER_PROTEIN WHERE P.ID_PROTEIN=?');
my $numObsProt=0;
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
		$numObsProt++;
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
#if ($Rdesign eq 'LABELED' && $numStates > 2) {
#	if ($analysis2Experiment{$anaID}) {$expCode=$analysis2Experiment{$anaID};}
#	else {
#		$expCode=$analysis2Experiment{$anaID}=$currentExperiment;
#		$currentExperiment++; # A++ = B!!!!
#	}
#}

		my $qSetStrg=($xicEngine{$xicQuantif}=~/^(PD|MQ)$/)? 'QSET' : "MCQSET_$xicQuantif"; # PEP_RATIO only
		##my %peptidesUsed; # only if non-proteotypic peptides are used => makes sure the same protein is selected (analysis-specific to respect Match group decision!!!)
		my $labelModStrg=join('|',keys %{$anaLabelModifs{$anaID}}) if $labeling eq 'SILAC'; # needed for SILAC
		my %mcqXICs; # SILAC only
		$sthGetPepData->execute($anaID);
#print DEBUG "RUN_ANA $anaID ($targetPos)\n" if !$referenceMode;
		my %usedPepQuantifValues; # prevent multiple insertion of same value due to peptides shared by multiple proteins if non-proteotypic allowed
		while (my ($protID,$pepID,$pepBeg0,$pepSeq,$varModStrg,$varProbStrg,$charge,$pepData,$score,$specif,$bestSpecif,$misscut)=$sthGetPepData->fetchrow_array) { # executed multiple time for labeled quantif
			$varModStrg='' unless $varModStrg;
			my $pepBeg=abs($pepBeg0); # $pepBeg0 required for ghost peptides detection

			#>Protein filtering
			next if $contaminants{$protID}; # only if option was checked
			if ($protSelectionType) {
				if ($protSelectionType eq 'exclude') {next if $selectExcludeProteins{$protID};}
				else {next unless $selectExcludeProteins{$protID};} # restrict
			}
			
#$varModStrg='' unless $varModStrg;
#print DEBUG "$pepID,$pepBeg,$pepSeq,$varModStrg\n" if $pepSeq eq 'FLYYEMGYK';
			next if ($manualPepSelection && (!$manualProteins{$protID} || !$manualPeptides{$pepID})); # reference mode

			if ($pepMissedCleavage==-1 && $misscut) {
				my $pepLength=length($pepSeq);
				$missedCutPeptides{$protID}{$pepBeg}=$pepLength if (!$missedCutPeptides{$protID} || !$missedCutPeptides{$protID}{$pepBeg} || $missedCutPeptides{$protID}{$pepBeg} < $pepLength);
				if ($normalizeWithAll) {$notUsed4Quantif{$pepID}=1;}
				else {next;}
			}

			next if ($labeling eq 'SILAC' && $pepSeq !~ /K|R/); # skip sequences not able to be labeled (C-term protein)

			next unless $pepQuantifValues{"$pepID:$targetPos"}; # no quantif for pepID
			my $quantifValue=$pepQuantifValues{"$pepID:$targetPos"}; # ref to last PEP/FRAGMENTS dimension

			### NEW (18/02/21) Protein enriched by specific PTM(s) & sequence motif(s) => restrict to matching peptides
			if ($isProteinByPTM) {
				my $okPeptide=&promsQuantif::checkPTMcontext($pepSeq,$varModStrg,$pepBeg,$protID,$sthProtSeq,\%proteinPTMcontexts,\%peptideContextMatch,\%proteinSequence);
				unless ($okPeptide) {
					if ($normalizeWithAll) {$notUsed4Quantif{$pepID}=1;}
					else {next;}
				}
			}

			#>Normalization scope
			$forNormalization{$protID}=1 if ($normProtUsage eq 'not' && !$notForNormalization{$protID});
	
			### NEW (25/01/19): Handle quantif PTM removed from samples => unmodifed res ARE the sites!
			## PTM recreated in silico
			if ($isFreeResQuantif) {
				# if ($matchProtCterm && !defined $proteinLength{$protID}) {
				# 	$sthProtL->executed($protID);
				# 	($proteinLength{$protID})=$sthProtL->fetchrow_array;
				# }
				# $varModStrg=&promsQuantif::modifyFreeResidues($pepSeq,$varModStrg,$pepBeg,$proteinLength{$protID},\%targetableRes,\%recreatedVarMods);
				$varModStrg=&promsQuantif::modifyFreeResidues($pepSeq,$varModStrg,$pepBeg,$protID,$sthProtSeq,\%proteinPTMcontexts,\%recreatedVarMods,\%proteinSequence);
			}

########!!!!!! TEST 18/08/15
#next if $preInfRatios{$protID};
###################################
#next if $excludedProteins{$protID}; # beta
			if (($isModifQuantif && (!$varModStrg || $varModStrg !~ /(^|&)($matchQuantifModStrg):/)) || ($varModStrg && $varModStrg=~/\d:(\.|&|$)/) ) { # missing mod pos data!
				if ($normalizeWithAll) {$notUsed4Quantif{$pepID}=1;}
				else {next;} # skip peptides not containing quantified modif
			}
#print DEBUG "OK FILTER_3\n" if $pepSeq eq 'FLYYEMGYK';			
			#if ($pepSpecificity ne 'unique' && !$manualPepSelection) {  # && $bestSpecif < 100 not proteotypic peptide && in manual mode (TDA ref): allow duplicate
				##if ($software eq 'MSstats') { # prot attribution must be the same across whole design (otherwise duplications occur)
				##	next if ($peptidesUsedMSstats{$pepSeq} && $peptidesUsedMSstats{$pepSeq} != $protID); # Use peptide ONLY once even if shared by multiple proteins ($pepSpecificity= unique_shared or all) (Required by ML algos)
				##	$peptidesUsedMSstats{$pepSeq}=$protID;
				##}
				##else {
				##	next if ($peptidesUsed{$pepSeq} && $peptidesUsed{$pepSeq} != $protID); # Use peptide ONLY once even if shared by multiple proteins ($pepSpecificity= unique_shared or all) (Required by ML algos)
				##	$peptidesUsed{$pepSeq}=$protID;
				##}
			#}			
#$proteinInConds{$protID}{$poolGroup2Cond{$observationSet}}=1; # protein "exists" in condition (detect before peptide filtering to prevent missleading info!!!)
#$visInReplicate{$protID}{$observationSet}=1; # protein is visible in replicate (VISIBILITY > 0 in query) needed for table.txt only
#$visInAnalysis{$protID}{$anaID}=1;
#$bestVisibility{$protID}=$vis if (!defined($bestVisibility{$protID}) || $bestVisibility{$protID} < $vis); # Keep the Best-Visibility of the protein (only needed for FREE).
			if ($normalizeWithAll) {
				$notUsed4Quantif{$pepID}=1 if (($pepSpecificity eq 'unique' && !$specif) || ($pepSpecificity eq 'unique_shared' && $bestSpecif<100)); # proteo-typic filter
				$notUsed4Quantif{$pepID}=1 if ($pepMissedCleavage <= 0 && (!defined $misscut || $misscut > 0)); # misscleavage filter (also skip if no info)
			}

			$pepData='' unless $pepData;
			$score=0 unless $score;
			
			$isGhostPeptide{$pepID}=1 if ($pepBeg0 < 0 && $countIdentifiedPep);
			
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
					my ($modID)=($vMod=~/^(-*\d+)/); # allow -1
					next if $labelModifs{$obsID}{$modID}; # skip labeling modifs
					next if $proteinPTMcontexts{$modID}; # skip context PTMs
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
#print DEBUG "OK FILTER_4\n" if $pepSeq eq 'FLYYEMGYK';
			# Non-SILAC phospho quantif
			if ($ptmProbSoftCode && $ptmProbThreshold > 0 && ($algoType ne 'PEP_RATIO' || $labeling ne 'SILAC') && !$notUsed4Quantif{$pepID}) { # SILAC case handled beloww
				my $ptmPosProb;
				if ($ptmProbSoftCode eq 'PRS') { 
					($ptmPosProb)=($pepData=~/PRS=\d;([^;]+);/);
				}
				else { # MQ|SPC|PTMRS
					$ptmPosProb=&promsMod::extractPtmPosProbability(\%modifQuantifRank,$varModCode,$varProbStrg); # returns 1 (100%) if no PTM in sites list
				}
				$ptmPosProb=-1 unless defined($ptmPosProb);
				if ($ambiguousModifPos) {
					my $pepIonKey="$pepSeq:$varModCode:$charge";
					$ptmProbability{$pepIonKey}=$ptmPosProb if (defined $ptmPosProb && (!defined $ptmProbability{$pepIonKey} || $ptmProbability{$pepIonKey} > $ptmPosProb)); # keep the worst across whole dataset
				}
				elsif (!defined $ptmPosProb || $ptmPosProb < $ptmProbThreshold) {
#print DEBUG "$pepID\t$varModCode\t$charge\t$varProbStrg\t$ptmPosProb\n";
					if ($normalizeWithAll) {$notUsed4Quantif{$pepID}=1;}
					else {next;} 
				}
			}
#print DEBUG "OK FILTER_5\n" if $pepSeq eq 'FLYYEMGYK';
			#$intensity{$pepID}=$quantifValue;
			#unless ($dataSrc) {$dataSrc=($labeling eq 'FREE')? '-' : $anaID;}
			#$dataSrc=($labeling eq 'FREE')? '-' : ($dataSrc)? $dataSrc : '#'.$anaID;
#$dataSrc=$dataSrc || '#'.$anaID; # define dataSrc for label-free also
			my $dataSrc;
# mPP----->
			#if ($labeling eq 'FREE' && $ratioType eq 'Ratio') {$dataSrc='-';}
			##if ($labeling eq 'FREE') {$dataSrc='-';} # <---- mPP
# <---- mPP
			if ($software eq 'MSstats') {
				$dataSrc='-';
			}
			else {
				$dataSrc='#'.$anaID; # TODO: if LABEL-FREE: check if still compatible with fraction merge ($anaInObsStrg in Qset) PP 27/11/17
				$dataSrc.='.'.$1 if $pepData=~/SOURCE_RANK=([^#]+)/;
			}
			my $qSet;
			if ($algoType eq 'PEP_RATIO') {

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
					push @{$tempIntensity{"$pepSeq:#:$varModCode:#:$charge"}{$dataSrc}{$fraction}{$qSet}},$quantifValue->{PEP} unless $usedPepQuantifValues{"$pepID:$targetPos"}; # $intensity computed later
					$dataSource2ObsSet{$dataSrc}{$observationSet}=$fraction;
if ($pepSpecificity ne 'unique' && !$manualPepSelection) {
	push @{$pepSeq2qSet{$pepSeq}},$qSet unless $usedPepQuantifValues{"$pepID:$targetPos"}; # For now multiple proteins can share the same qSet!!!
	$qSet2AllProteins{$qSet}{$protID}=1;
#print DEBUG "SHARED_PROT=$protID Q=$qSet SEQ=$pepSeq VM=$varModCode\n" if $pepSeq eq 'MSGFIYQGK';
}
else {
					$qSet2Protein{$qSet}=$protID;
}
					$pepID2QuantifSetSILAC{$pepID}=$qSet;
				}
				else {
					push @{$intensity{"$pepSeq:$varModCode:$charge"}{$dataSrc}{$qSet}{$observationSet}{PEP}},$quantifValue->{PEP} unless $usedPepQuantifValues{"$pepID:$targetPos"};
				}
			}
			else { # PEP_INTENSITY: No peptide pairing   former $labeling eq 'FREE' 27/04/20
# mPP----->
#$qSet=($ratioType eq 'Ratio')? 'A0_Q0' : "P$pepID"."_Q0"; # <---- mPP
#$qSet=($ratioType eq 'Ratio')? 'A0_Q0' : "A$anaInObsStrg"."_Q0"; # PP 10/08/17 Revert to "A$anaInObsStrg"."_Q0" to allow fractions merge !!!!!!!!!!!!!!! PP 18/04/16: anastrg incompatible with best Qset varMod coherence check
#$qSet=($ratioType eq 'Ratio')? 'A0_Q0' : "A$anaID"."_Q0"; #"A$anaID"."_Q0";!!!!!!!!!!!!!!! PP 10/08/17 WARNING: Should not be compatible with best Qset varMod coherence check
# <---- mPP
				if ($software eq 'MSstats') {
					$qSet='-';
				}
				elsif ($poolObservations) {
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

				unless ($usedPepQuantifValues{"$pepID:$targetPos"}) {
					foreach my $feature (keys %{$quantifValue}) { # PEP of Fragments
						$intensity{"$pepSeq:$varModCode:$charge"}{$dataSrc}{$qSet}{$observationSet}{$feature}[0]+=$quantifValue->{$feature};
						$fragmentsInPepIon{"$pepSeq:$varModCode:$charge"}{$feature}=1 if $software eq 'MSstats';
					}
				}
			}
		
			#my $usedObs=($ratioType eq 'Ratio')? 'ALL_OBS' : $observationSet;
			#push @{$proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet}{$usedObs}},$pepID;
			push @{$proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet}{$observationSet}},$pepID;
			if ($pepSpecificity ne 'unique' && !$manualPepSelection && $MGSharedPep ne 'share') { # not proteotypic peptide && in manual mode (TDA ref): allow duplicate
				$peptideSpecificity{$pepSeq}{$protID}++;
				$bestProteotypicity{$protID}=$bestSpecif if (!$bestProteotypicity{$protID} || $bestProteotypicity{$protID} < $bestSpecif);
			}
			$peptideBeg{$protID}{$pepSeq}=$pepBeg; # for modif quantification only
			if ($pepMissedCleavage==-1 && !$misscut) {push @{$beg2pepSeq{$protID}{"$pepBeg:$pepSeq"}},$pepID;} # no need to include missed cleaved here
			
			next if $usedPepQuantifValues{"$pepID:$targetPos"}; # Do not add it again
			$usedPepQuantifValues{"$pepID:$targetPos"}=1;
			
			##>Record best varMod of qSet in case incoherence (assumes same pepSeq!)
			#if ($labeling ne 'FREE') { #} qSET has no meaning in label-free ***Also for FREE because of possible duplicate scan for same pepID with MassChroQ***
			if ($algoType eq 'PEP_RATIO' && $labeling eq 'SILAC') { # Specific to PEP_RATIO with SILAC (isotope labeling: qSet links multiple pepIDs)
				
##>Records duplicate XIC scans if MCQ (eg multiple PD XICs reunited in a single one by MCQ => keep only 1 qSet)
if ($xicEngine{$xicQuantif} eq 'MCQ') {
	$mcqXICs{"$pepSeq:$varModCode:$charge:".$quantifValue->{PEP}}{$qSet}=1; # duplicate if more than 1 qSet with the same pepSeq,charge (ideally mass) have exactly the same XIC value
} # assumes EXACTLY same XIC is computed by MCQ		
				
				if ($ptmProbSoftCode) { # Position confidence filtered earlier for non-SILAC. !$ptmProbSoftCode if $isFreeResQuantif
					my $ptmPosProb;
					if ($ptmProbSoftCode eq 'PRS') {
						($ptmPosProb)=($pepData=~/PRS=\d;([^;]+);/);
					}
					else { # MQ|PTMRS (SPC unlikely for SILAC)
						$ptmPosProb=&promsMod::extractPtmPosProbability(\%modifQuantifRank,$varModCode,$varProbStrg);
					}
					$ptmPosProb=-1 unless defined($ptmPosProb);
					@{$qSetBestVarMod{$qSet}}=($ptmPosProb,$varModCode) if (!$qSetBestVarMod{$qSet} || $qSetBestVarMod{$qSet}[0] < $ptmPosProb); # record best ptm probability & vmod (both peptides in qSet will be filtered based on best prob)
#print DEBUG "$pepID\t$varModCode\t$charge\t$varProbStrg\t$ptmPosProb\n" if (!defined $ptmPosProb || $ptmPosProb < $ptmProbThreshold);
				} 
				else {
					@{$qSetBestVarMod{$qSet}}=($score,$varModCode) if (!$qSetBestVarMod{$qSet} || $qSetBestVarMod{$qSet}[0] < $score); # record best
				}
			}
			#if ($ambiguousModifPos) {
			#	$retentionTime{$pepID}=($elutionStrg && $elutionStrg=~/et([\d\.]+)/)? $1 : ($elutionStrg && $elutionStrg=~/^([\d\.]+)/)? $1 : 0;
			#}
#print DEBUG "OK\n" if $pepSeq eq 'FLYYEMGYK';
		} # end of while

		##>Remove duplicate XIC scans if MCQ
		if ($algoType eq 'PEP_RATIO' && $labeling eq 'SILAC') {
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

		unless ($numObsProt % 10) { # every 10 obs 
			open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
			print FILESTAT "1/$numSteps Generating data files (2/4 Retrieving$referenceStrg identification data: $numObsProt/$numAllObs)\n";
			close FILESTAT;
		}
	}
}
$sthGetPepData->finish;
# $sthProtL->finish;
$sthProtSeq->finish;
# $sthPTM->finish;
# %pepQuantifValues=(); # hoping to free memory

open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
print FILESTAT "1/$numSteps Generating data files (3/4 Processing$referenceStrg dataset)\n";
close FILESTAT;

# Count initial number of potential peptides identified for each prot (before attribution)
my %protPepNb;
foreach my $protID (keys %proteins) {
	$protPepNb{$protID} = scalar keys %{$proteins{$protID}};
}

##>Attributing shared peptides to their UNIQUE "best" protein (higher num pepSeq > higher num data points in design > smaller id)
my %sharedPepProteins;
if ($pepSpecificity ne 'unique' && !$manualPepSelection && !$MGSharedPep ne 'share') { # MGSharedPep can be best|exclude
	foreach my $pepSeq (keys %peptideSpecificity) {
		next unless scalar keys %{$peptideSpecificity{$pepSeq}} > 1;
		my @excludedProt=(sort{$bestProteotypicity{$a} <=> $bestProteotypicity{$b} || $protPepNb{$b} <=> $protPepNb{$a} || $peptideSpecificity{$pepSeq}{$b} <=> $peptideSpecificity{$pepSeq}{$a} || $a <=> $b} keys %{$peptideSpecificity{$pepSeq}});
#print DEBUG "ALL_SHARED_PROT=",join(',',@excludedProt)," SEQ=$pepSeq\n" if scalar @excludedProt > 1; #$pepSeq eq 'MSGFIYQGK';		
		my $keptProtID=($MGSharedPep eq 'best')? shift @excludedProt : 0; # for THIS pepSeq: 'best': remove "best" protein from rejection list OR 'exclude': keep all proteins in rejection list
		foreach my $protID (@excludedProt) {
			delete $proteins{$protID}{$pepSeq}; # remove pepSeq from list of matching sequences
			if ($pepMissedCleavage==-1) {
				my $pepBeg=$peptideBeg{$protID}{$pepSeq};
				delete $beg2pepSeq{$protID}{"$pepBeg:$pepSeq"} if ($beg2pepSeq{$protID} && $beg2pepSeq{$protID}{"$pepBeg:$pepSeq"});
				delete $missedCutPeptides{$protID}{$pepBeg} if ($missedCutPeptides{$protID} && $missedCutPeptides{$protID}{$pepBeg});
			}
			delete $peptideBeg{$protID}{$pepSeq};
			$sharedPepProteins{$protID}{$keptProtID}=1;
		}
	}

	my %rejectedProteins;
	foreach my $protID (keys %sharedPepProteins) {
		if (scalar keys %{$proteins{$protID}}==0) { # no pepSeq left
			delete $proteins{$protID};
			delete $beg2pepSeq{$protID};
			delete $peptideBeg{$protID};
			delete $missedCutPeptides{$protID};
			$rejectedProteins{$protID}=1;
		}
		if (scalar keys %rejectedProteins) {
			open(REJECT,">$dataDir/rejectedProteins.txt"); # list of exluded proteins due to de-attribution of peptides shared with other proteins
			print REJECT "Protein_ID\tOther\n";
			foreach my $protID (keys %rejectedProteins) {
				print REJECT "$protID\t",join(',',keys %{$sharedPepProteins{$protID}}),"\n"; # "protID 0" if $MGSharedPep is 'exclude'
			}
			close REJECT;
		}
	}
	
	#>Build final qSet2Protein (SILAC only)
	if ($labeling eq 'SILAC') {
		foreach my $pepSeq (keys %pepSeq2qSet) {
			my $protID=(keys %{$peptideSpecificity{$pepSeq}})[0]; # there should be only 1 protein left!!!
			foreach my $qSet (@{$pepSeq2qSet{$pepSeq}}) {
				$qSet2Protein{$qSet}=$protID if $qSet2AllProteins{$qSet}{$protID}; # <= Protein may be missing (hidden) in qSet's observation => qSet is orphan => not used!!!
			}
		}
		undef %pepSeq2qSet;
	}
	
	undef %peptideSpecificity;
	undef %qSet2AllProteins;
}

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

if ($algoType eq 'PEP_RATIO' && $labeling eq 'SILAC') { # former $labeling eq 'SILAC' 27/04/20

	##>Checking qSet sequence
	foreach my $qSet (keys %qSetSequence) { # It happens that identifications are different between Light and (Medium/)Heavy spectra!!!!
		$qSetIsBad{$qSet}=1 if scalar keys %{$qSetSequence{$qSet}} > 1;
#print DEBUG "BAD_QSET=$qSet\n" if scalar keys %{$qSetSequence{$qSet}} > 1;
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

#print DEBUG "QSET_OBS: S=$dataSrc VM=$varModCode UVM=$usedVmod QSETs=$qSetMerged ($pepSeq)\n" if $pepSeq eq 'MSGFIYQGK'; 
				##>Merged qSet
				if ($qSetMerged=~/\+/) {
					$qSetMerged=~s/\+A\d+_/\+/g; # remove repeated A<anaID>_ except 1srt one: Axxx_Qyyy+Axxx_Qzzz -> Axxx_Qyyy+Qzzz
					foreach my $observationSet (keys %{$dataSource2ObsSet{$dataSrc}}) {
						next if (!$sumQuantifValues{$observationSet} || !$sumQuantifValues{$observationSet}{$usedVmod});

						$intensity{"$pepSeq:$usedVmod:$charge"}{$dataSrc}{$qSetMerged}{$observationSet}{PEP}[0]=$sumQuantifValues{$observationSet}{$usedVmod};

						##>Correcting %proteins with merged qSets
						my $protID=$qSet2Protein{$qSetObsList{$observationSet}{$usedVmod}[0]};
#print DEBUG "*ORPHAN QSETM=$qSetMerged SEQ=$pepSeq VM=$varModCode UVM=$usedVmod Z=$charge S=$dataSrc\n" unless $protID;
						next unless $protID; # orphan qSet (best protein is missing: only for "shared peptides")
						my @pepIdList;
						my $okObs=0;
						foreach my $qSet (@{$qSetObsList{$observationSet}{$usedVmod}}) {
							next if ($sharedPepProteins{$protID} && (!$proteins{$protID} || !$proteins{$protID}{$pepSeq}));	# shares peptides with other proteins & has no unshared peptides
							push @pepIdList,$proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet}{$observationSet}[0];
							$okObs=1;
#print DEBUG "+PROT=$protID Q=$qSet PEP=$pepSeq VM=$varModCode Z=$charge S=$dataSrc OBS=$observationSet\n" unless $proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet}{$observationSet}[0];
						}
						
						$proteins{$protID}{$pepSeq}{$usedVmod}{$charge}{$dataSrc}{$qSetMerged}{$observationSet}[0]=join('+',@pepIdList) if $okObs;

#print DEBUG "+NO: PROT=$protID Q=$qSetMerged PEP=$pepSeq VM=$varModCode Z=$charge S=$dataSrc OBS=$observationSet\n" unless $okObs;

					}

					foreach my $qSet (keys %{$qSetList{$usedVmod}}) {
						##>Deleting obsolete qSets
						my $protID=$qSet2Protein{$qSet};
						delete $proteins{$protID}{$pepSeq}{$varModCode}{$charge}{$dataSrc}{$qSet} if ($protID && $proteins{$protID} && $proteins{$protID}{$pepSeq});

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
#print DEBUG "*ORPHAN QSET=$qSet SEQ=$pepSeq VM=$varModCode UVM=$usedVmod Z=$charge S=$dataSrc\n" unless $protID;
					next unless $protID; # orphan qSet (best protein is missing: only for "shared peptides")
					foreach my $observationSet (keys %{$dataSource2ObsSet{$dataSrc}}) {
						next if (!$sumQuantifValues{$observationSet} || !$sumQuantifValues{$observationSet}{$usedVmod}); # no data for this channel
						next if ($sharedPepProteins{$protID} && (!$proteins{$protID} || !$proteins{$protID}{$pepSeq}));	# shares peptides with other proteins & has no unshared peptides
						$intensity{"$pepSeq:$usedVmod:$charge"}{$dataSrc}{$qSet}{$observationSet}{PEP}[0]=$sumQuantifValues{$observationSet}{$usedVmod}; # if ($sumQuantifValues{$observationSet} && $sumQuantifValues{$observationSet}{$usedVmod});
						if ($usedVmod ne $varModCode) {
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
my $numProteins=scalar keys %proteins;
if ($numProteins == 0) { # No data at all!!!
	my $dataType=($referenceMode eq 'REF')? ' protein reference' : ($referenceMode eq 'NORM')? ' normalizing peptide' : '';
	warn "ERROR: No$dataType data available for quantification!";
	&checkForErrors($dbh);
}
if ($software eq 'myProMS') {
	my $computeiBAQ = 0;
	my $computeAbsolute = 0;
	foreach my $meas (@{$quantifParameters{'DB'}{'ABUND_MEASURES'}}) {
		if ($meas =~ /^aLFQ_iBAQ/) {
			$computeiBAQ = 1;
			$computeAbsolute = 1;
			last;
		} elsif ($meas =~ /^aLFQ_TOP/) {
			$computeAbsolute = 1;
			last;
		}
	}

	if ($computeAbsolute) {  # Absolute quantification requires some more processing
		if (-d $absDir && -f "$absDir/aLFQ_info.txt") {
			# Must convert user quantities file to custom format
			my $aLFQInfoFile  = "$absDir/aLFQ_info.txt";
			my $userQtyFile   = "$absDir/user_prot_quantities.csv";
			my $customQtyFile = "$dataDir/prot_quantities.csv";
			my $protMWFile 	  = "$dataDir/prot_molecular_weights.tsv";

			open(aLFQ_INFO, ">>$aLFQInfoFile");  # The first part of the file was written in startDesignQuantification
			# Only the molar/mass percent will be computed if user didn't provide quantities
			print aLFQ_INFO "quantities\t$customQtyFile\n" if (-f $userQtyFile);
			# If absolute quantif is from DDA experiment, we start from resultsPep.txt, not when from DIA
			print aLFQ_INFO "resultsPep\t$resultDir/resultsPep.txt\n" if ($algoType eq 'PEP_INTENSITY');
			print aLFQ_INFO "fasta\t$dataDir/protein_iBAQ.fasta\n" if ($computeiBAQ);
			print aLFQ_INFO "molecular_weights\t$protMWFile\n";
			close aLFQ_INFO;

			my %protAlias2IDs;
			my %protMolWeights;
			my %protSequences;
			my $sthAlias = $dbh->prepare("SELECT P.ID_PROTEIN, P.ALIAS, P.PROT_SEQ, MP.PROT_SEQ, P.MW, MP.MW FROM PROTEIN P LEFT JOIN MASTER_PROTEIN MP ON P.ID_MASTER_PROTEIN = MP.ID_MASTER_PROTEIN WHERE ID_PROTEIN IN (".join(',', keys %proteins).")");
			$sthAlias->execute;
			while (my ($protID, $alias, $protSeq, $masterProtSeq, $protMW, $masterProtMW) = $sthAlias->fetchrow_array) {
				$proteinAlias{$protID} = $alias;
				$protAlias2IDs{$alias} = $protID;

				if ($protMW) {
					$protMolWeights{$protID} = $protMW;
				} elsif ($masterProtMW) {
					$protMolWeights{$protID} = $masterProtMW;
				} else {
					next;  # Discard protein if no molecular weight (mass computation becomes meaningless)
				}

				if ($computeiBAQ) {  # Get sequences for fasta file if quantification uses aLFQ_iBAQ
					if ($protSeq && $masterProtSeq) {
						$protSequences{$protID} = ($protSeq eq '+')? $masterProtSeq : $protSeq;
					} elsif ($masterProtSeq) {
						$protSequences{$protID} = $masterProtSeq;
					} elsif ($protSeq && $protSeq ne '+' && $protSeq ne '-') {
						$protSequences{$protID} = $protSeq;
					} else {
						delete $protMolWeights{$protID} if $protMolWeights{$protID};
						next;  # Discard protein if no sequence (aLFQ computation is meaningless if included in fasta)
					}
				}
			}
			$sthAlias->finish;

			open(MOL_WEIGHTS,">$protMWFile");  # Record molecular weight of each prot to compute mass%
			print MOL_WEIGHTS "ID_PROTEIN\tMW\n";
			foreach my $protID (sort{$a<=>$b} keys %protMolWeights) {
				print MOL_WEIGHTS "$protID\t$protMolWeights{$protID}\n";
			}
			close MOL_WEIGHTS;

			if ($computeiBAQ) {  # Creation of the fasta file if quantification uses aLFQ_iBAQ
				open(IBAQ_FASTA,">$dataDir/protein_iBAQ.fasta");
				foreach my $protID (sort{$a<=>$b} keys %protSequences) {
					print IBAQ_FASTA ">$protID\n$protSequences{$protID}\n";
				}
				close IBAQ_FASTA;
			}

			&formatQtyFile($quantifID, $userQtyFile, $customQtyFile, \%protAlias2IDs, \%protMolWeights) if (-f $userQtyFile);

		} else {
			warn "$absDir/aLFQ_info.txt was not written, impossible to launch absolute quantification without the required information !";
		}

	} else {
		my $sthAlias=$dbh->prepare("SELECT ID_PROTEIN,ALIAS FROM PROTEIN WHERE ID_PROTEIN IN (".join(',',keys %proteins).")");
		$sthAlias->execute;
		while (my ($protID,$alias)=$sthAlias->fetchrow_array) {
			$proteinAlias{$protID}=$alias;
		}
		$sthAlias->finish;
	}
}

#print DEBUG scalar keys %proteins," proteins\n";
#close DEBUG if !$referenceMode;
#die "END OF TEST";

##############################
####> Printing data table<####
##############################
open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
print FILESTAT "1/$numSteps Generating data files (4/4 Printing$referenceStrg dataset to file)\n";
close FILESTAT;
my $protValidity=($referenceMode eq 'NORM')? 0 : 1;
my (%infiniteRatioProteins,%noSuperRatioProteins); # populated by &checkInfiniteRatios as globals (infinite ratio proteins are allowed in data table but stats are ignored)
my %qSetUsage=('USED'=>{},'EXCLUDED'=>{}); # in case  qSet filtering (eg. PhosphoRS)
my $protNumber=0;
#my %prot2quantify;
my $countCompletePeptides = 0;
my $numPepValues=0;

####>Main data file
my $dataFile=($referenceMode eq 'REF')? "$dataDir/tableRef.txt" : ($referenceMode eq 'NORM')? "$dataDir/tableNorm.txt" : "$dataDir/table.txt";
open(DATA,">$dataFile"); # Name of data table
if ($software eq 'MSstats') {
	print DATA "ProteinName\tPeptideSequence\tPrecursorCharge\tFragmentIon\tProductCharge\tIsotopeLabelType\tCondition\tBioReplicate\tRun\tIntensity\n";
}
else {
	print DATA "Protein_ID\tPeptide\tSample\tTreatment\tReplicate\tTechnical_Replicate\tExperiment\tQuantif_Set\tPeptide_IDs\tProtein_Name\tProtein_Validity\tValue" unless $referenceMode eq 'NORM';
}
my (%protein2Sites); #,%proteinSharing);
#my %pepList; # MSstats
my $pcStep=1;
my $pcCurrent=0;
my $countStep=($numProteins/100)*$pcStep;
my $nextStep=$countStep;
my $protCount=0;
foreach my $protID (sort{$a<=>$b} keys %proteins) {
#last if $protNumber >=200; #!!!!!!!!!!!!!!!!
#print DEBUG ">$protID:\n";
	$protCount++;
	if ($protCount >= $nextStep) {
		$pcCurrent+=$pcStep;
		if ($numProteins >= 1000) {
			open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
			print FILESTAT "1/$numSteps Generating data files (4/4 Printing$referenceStrg dataset to file: $pcCurrent%)\n";
			close FILESTAT;
		}
		$nextStep+=$countStep;
	}
	my (%usableData,%numPepInCond,%numPepInCond4ModProt);
	foreach my $pepSeq (keys %{$proteins{$protID}}) { # %{$peptideBeg{$protID}} safer than %{$proteins{$protID}} in case shared peptides are used
		next if $excludedSeq{$pepSeq}; # excluded by ptmAllowed <= -1 or shared by multiple top proteins
		#$prot2quantify{$protID}=1;
#unless ($peptideBeg{$protID}{$pepSeq}) { # To prevent use of theorically deleted $proteins{$protID}{$pepSeq}
#	print DEBUG "ERR: PROT=$protID SEQ=$pepSeq\n";
#	next;
#}
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
								foreach my $refIntens (values %{$intensity{$pepIonKey}{$dataSrc}{$qSet}{$observationSet}}) { # PEP or TDA/DIA=>Sum accross all fragments
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

	if ($algoType ne 'PEP_RATIO') { # potential data filtering ---> 27/04/20
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
						if (scalar keys %pepIonInCond < $numStates) { # (< 2 keep any pairs? for 3+ conditions) peptide is missing in 1+ cond => flag to be deleted
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
#print DEBUG "OK FILTER_6\n" if $pepSeq eq 'FLYYEMGYK';				
			}
		}

		##>Update usableData
		foreach my $dataKey (keys %excludedData) {
			my ($qSet,$observationSet,$pepSeq,$vmod,$charge,$dataSrc)=split(':#:',$dataKey);
			my $pepIonKey="$pepSeq:$vmod:$charge";
#print DEBUG "OUT FILTER_7.0\n" if $pepSeq eq 'FLYYEMGYK';
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
#print DEBUG "OUT FILTER_7.1\n" if $pepSeq eq 'FLYYEMGYK';
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
#print DEBUG "OK DATA_1\n" if $pepSeq eq 'FLYYEMGYK';
		my $pepSeqBeg=$peptideBeg{$protID}{$pepSeq};
		my $pepSeqContext=($keepSeqContext && !$referenceMode)? '['.$pepSeqBeg.'.'.length($pepSeq).']' : ''; # only if $isModifQuantif
		my @seq=split(//,$pepSeq);
		foreach my $vmod (keys %{$usableData{$pepSeq}}) {
			my $nonAmbigModProtID=$protID; # default
			my $ambigModProtID=$protID; # default
			my $seqIsAmbig=0;
			if ($isModifQuantif) { # modif quantification
#>Required for conversion to function: Externals: $protID, $vmod, @modifQuantifs, @seq, $pepSeqBeg, $quantifiedModifID, %modifQuantifRank, %qQuantifResMatchStrg, %quantifResMatchStrg, $ambiguousModifPos
#>Returns: $nonAmbigModProtID, $ambigModProtID, $seqIsAmbig
				my $firstModQ=1;
				foreach my $modID (@modifQuantifs) {
					my ($sep,$posStrg)=($vmod=~/(^|&)$modID:([\d+\.\-\=\+\*]+)/); # added N/C-term
					next unless $posStrg; # in multi-modif quantif context: varmod does not necessarly contains all quantified modif
					#>Non-ambiguous pos
					my @pos=(); # reset
					foreach my $pepPos (split(/\./,$posStrg)) {
						if ($pepPos=~/\d+/) {push @pos,$seq[$pepPos-1].($pepSeqBeg+$pepPos-1);}
						#elsif ($pepPos=~/[=\*]/) {next;} # skip peptide N/C-term: cannot be in protein TODO: add to table.txt for data normalisation but do not quantify
						elsif ($pepPos eq '=') {push @pos,'n'.$pepSeqBeg;} # peptide N-term: kept for normalization purpose
						elsif ($pepPos eq '-') {push @pos,'n0';} # protein N-term: kept for normalization purpose
						elsif ($pepPos eq '*') {push @pos,'c'.($pepSeqBeg+$#seq);} # peptide C-term: kept for normalization purpose
						elsif ($pepPos eq '+') {push @pos,'c0';} # protein C-term: kept for normalization purpose
					}
					#>THIS IS WHERE THE PROTEIN & MODIFICATION SITE(S) ARE CODED!!!!!!
					# Single-modif quantif:
					#	-Non-ambiguous: "ProtID-Res1Pos1[.Res2Pos2...]" eg "1234-S56.T78"
					#	-Ambiguous: "ProtID-StartRes~EndRes:NumModifs/NumAcceptors" eg "1234-56~78:2/5"
					# Mutli-modif quantif:
					#	-Non-ambiguous: "ProtID-ModifXRank#Res1Pos1[.Res2Pos2...][&ModifYRank#Res3Pos3[.Res4Pos4...]]" eg "1234-1#S56.T78&2#C61" sorted by increasing rank!!!
					#	-Ambiguous: "ProtID-ModifXRank#StartRes~EndRes:NumModifs/NumAcceptors" eg "1234-1#56~78:2/5&2#C61" (PP: Not allowed from start form as of 30/06/20)
					my $nonAmbigMod=($quantifiedModifID)? join('.',@pos) : $modifQuantifRank{$modID}.'#'.join('.',@pos); # Add '<modif rank>#' if multi-modif quantif
					$nonAmbigModProtID.=($firstModQ)? '-' : '&';
					$nonAmbigModProtID.=$nonAmbigMod;
					
					# ==> Look for ambiguity
			if ($ambiguousModifPos) {
					$ambigModProtID.=($firstModQ)? '-' : '&';
					$ambigModProtID.=$modifQuantifRank{$modID}.'#' unless $quantifiedModifID;
					#my $numModifRes = () = $posStrg =~ /(\d+)/g;
					my $numModifRes=scalar @pos;
					my @targetResIdx;
					while ($pepSeq =~ /[$qQuantifResMatchStrg{$modID}]/g) {push @targetResIdx,$-[0];} # position of match (1st=0)
					my $numTargetRes=scalar @targetResIdx;
					my ($NtermMatch,$CtermMatch)=('','');
					if ($quantifResMatchStrg{$modID}=~/=/) {$NtermMatch='n'.$pepSeqBeg;} # peptide N-term
					elsif ($quantifResMatchStrg{$modID}=~/\*/) {$CtermMatch='c'.($pepSeqBeg+$#seq);} # peptide C-term
					elsif ($quantifResMatchStrg{$modID}=~/-/) {$NtermMatch='n0' if $pepSeqBeg==1;} # protein N-term
					#elsif ($quantifResMatchStrg=~/\+/) {$CtermMatch='c0' if $pepSeqBeg+length($pepSeq)-1==<protein length>;} # TODO: fetch prot length
					my $numAlltargets=$numTargetRes;
					$numAlltargets++ if $NtermMatch;
					$numAlltargets++ if $CtermMatch;
					if ($numModifRes < $numAlltargets) { # Ambiguity! protID-[modifPos=]<1stPos>~<lastPos>:<numModif>/<numSites>
						$seqIsAmbig=1;
						if ($ambiguousModifPos) {
							#$ambigModProtID=$protID.'-';
							if ($NtermMatch) {$ambigModProtID.=$NtermMatch;} # N-term
							else {$ambigModProtID.=($pepSeqBeg+$targetResIdx[0]);} # internal pos
							$ambigModProtID.='~';
							if ($CtermMatch) {$ambigModProtID.=$CtermMatch;} # C-term
							else {$ambigModProtID.=($pepSeqBeg+$targetResIdx[-1]);} # internal pos
							$ambigModProtID.=':'.$numModifRes.'/'.$numAlltargets;
						}
						else {$ambigModProtID.=$nonAmbigMod;}
					}
					else {$ambigModProtID.=$nonAmbigMod;}
			}
					$firstModQ=0;
				}
			}
			foreach my $charge (sort{$a<=>$b} keys %{$usableData{$pepSeq}{$vmod}}) {
				my $pepIonKey="$pepSeq:$vmod:$charge";
				foreach my $dataSrc (sort keys %{$usableData{$pepSeq}{$vmod}{$charge}}) {
					my $infRatioDataSrc=($algoType ne 'PEP_RATIO')? '-' : $dataSrc; # ignore dataSrc if label-free former $labeling eq 'FREE' 27/04/20

					foreach my $qSet (sort{$a cmp $b} keys %{$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}}) { # sort cmp because Qset can be string
						next if $qSetIsBad{$qSet}; # only for PEP_RATIO SILAC or duplicate MCQ
						my $modProtID=$nonAmbigModProtID; # default
						my $usedVmod=($algoType eq 'PEP_RATIO' && $labeling eq 'SILAC')? $qSetBestVarMod{$qSet}[1] : $vmod; # can be different from $vmod in case qSet incoherence (eg. phospho pos after phosphoRS update on ambiguous pepSeq)
						if ($seqIsAmbig) {
							if ($ptmProbSoftCode) {
								if ($algoType eq 'PEP_RATIO' && $labeling eq 'SILAC') {
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
													#foreach my $pos (split(/\./,$posStrg)) {push @pos,$seq[$pos-1].($pepSeqBeg+$pos-1);}
													foreach my $pos (split(/\./,$posStrg)) {
														# if (!$pos || ($pos=~/\d+/ && (!$seq[$pos-1])) || !$pepSeqBeg) { # (PP 2017/06/22)
														# 	warn "PTM position Error: PROT=$protID, PEP=$pepSeq, VMOD=$vmod, UVMOD=$usedVmod, QSET=$qSet, SEP=$sep, STRG=$posStrg, BEG=$pepSeqBeg\n";
														# 	&checkForErrors($dbh);
														# }
														push @pos,$seq[$pos-1].($pepSeqBeg+$pos-1);
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
								elsif ($ambiguousModifPos &&  $ptmProbability{$pepIonKey} < $ptmProbThreshold) { # for everything except 'PEP_RATIO' with 'SILAC'
									$modProtID=$ambigModProtID; # WARNING: ion is ambiguous across whole dataset if its ptm prob < $ptmProbThreshold in at least 1 obs!!!!
								}
							}
							else {$modProtID=$ambigModProtID;} # non-Phospho PTM quantified as ambiguous
						}
						
						###>Injecting sequence context to PTM
						$modProtID.=$pepSeqContext;
						
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
###						if (scalar keys %pepIonInCond < $numStates) { # (< 2 keep any pairs? for 3+ conditions) peptide is missing in 1+ cond => flag to be deleted
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
						#my %pepIsComplete; # records in which cond peptide has value

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
					
					
							#>MSstats
							if ($software eq 'MSstats') {
								foreach my $feature (sort{&promsMod::sortSmart($a,$b)} keys %{$fragmentsInPepIon{$pepIonKey}}) { # sort: type/res/charge
									my ($fragType,$fragCharge)=split('#',$feature);
									
									foreach my $observationSet (@quantiOrder) {
										# No longer need to fill missing values with "NA" (PP 12/08/20)
										if ($intensity{$pepIonKey} && $intensity{$pepIonKey}{'-'}{'-'}{$observationSet} && $intensity{$pepIonKey}{'-'}{'-'}{$observationSet}{$feature}) {
											my $fragArea=$intensity{$pepIonKey}{'-'}{'-'}{$observationSet}{$feature}[0]*1 || 'NA'; # extra check
											next if $fragArea eq 'NA';
											#$pepList{$modProtID}{$observationInfo{$observationSet}[0]}{"$pepSeq$usedVmod\_$charge"}{"$fragType\_$fragCharge"}=1; # To count peptides, can be created after R by reading table.txt
											
											$statesInFile{$observationInfo{$observationSet}[0]}=1; # number of states with valid data
if ($normProtUsage) { # Normalization uses a subset of proteins/sites
	#if ($sharedProtNorm) {$proteinSharing{$modProtID}{$observationSet}=1;} # Normalization is restricted to shared proteins/sites across all techReps ($observationSet)
	#els
	if ($isModifQuantif && !$quantifParameters{'DB'}{'PEPTIDES_NORM'}) {$protein2Sites{$protID}{$modProtID}=1;} # to convert protID from normalization list to modProtID in quantif
}																						
											#           ProteinName PeptideSequence PrecursorCharge FragmentIon ProductCharge IsotopeLabelType Condition(StateX) BioReplicate(StateX_BioRepY)          Run                                   Intensity
											print DATA "$modProtID\t$pepSeq$usedVmod\t$charge\t$fragType\t$fragCharge\tL\t$observationInfo{$observationSet}[0]\t$observationInfo{$observationSet}[1]\t$observationInfo{$observationSet}[2]\t$fragArea\n";
										}
									}
								}	
								next; # qSet only '-' for MSstats => all myProMS code below is skipped
							}
									
							#>myProMS
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
my $expCode=($algoType ne 'PEP_RATIO')? 'A' : $ana2Experiment{$anaID}; #$techRep2Experiment{$observationSet}; # former  $labeling eq 'FREE' 27/04/40
##my $expCode='A'; # default. for LABELFREE
##if ($Rdesign eq 'SUPERRATIO') {$expCode=$observationInfo{$observationSet}[3];}
##elsif ($Rdesign eq 'LABELED') {
##	my $pepID=$proteins{$protID}{$pepSeq}{$vmod}{$charge}{$dataSrc}{$qSet}{$observationSet}[0];
##	$expCode=$peptideID2Experiment{$pepID} || 'A'; # just to be safe
##}

									#$pepIsComplete{$obs2Cond{$observationSet}}{$observationSet} = 1;
									my $sep=($algoType ne 'PEP_RATIO')? '+' : '='; # XICs are summed for label FREE, combined for SILAC former $labeling eq 'FREE' 27/04/20
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

if ($normProtUsage) { # Normalization uses a subset of proteins/sites
	#if ($sharedProtNorm) {$proteinSharing{$modProtID}{$observationSet}=1;} # Normalization is restricted to shared proteins/sites across all techReps ($observationSet)
	#els
	if ($isModifQuantif && !$quantifParameters{'DB'}{'PEPTIDES_NORM'}) {$protein2Sites{$protID}{$modProtID}=1;} # to convert protID from normalization list to modProtID in quantif
}
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
									# Ratios are no longer Experiment-specific in SuperRatio => cannot use $usedExpCode as synonym of ratio
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
							} # end of observationSet
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
	if ($ratioType ne 'None' && $software eq 'myProMS') {
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
	if ($referenceMode) {
		$dbh->disconnect;
		my $fileType=($referenceMode eq 'REF')? 'protein reference' : 'normalizing peptide';
		die "ERROR: Less than 2 state data found in $fileType file!<BR>Maybe not enough peptides passed your selection criteria."; # error will be handled by parent (MASTER) process
	}
	# MASTER here
	warn "ERROR: Less than 2 state data found in data file!<BR>Maybe not enough peptides passed your selection criteria.";
	sleep 2;
	&checkForErrors($dbh);
}

if ($referenceMode) { ### <- End of Reference mode (REF|NORM)!
	$dbh->disconnect;
#close DEBUG;
	exit;
}

######################################################
####> Merge tmp tableNorm.txt file into table.txt<#### protein_validity=0 => not included in diff ana process
######################################################
if ($quantifParameters{'DB'}{'PEPTIDES_NORM'}) { ### <- Merge tableNorm.txt with table.txt
	system "cat $dataDir/tableNorm.txt >> $dataFile";
	unlink "$dataDir/tableNorm.txt";
}

######################################
####> Printing normalization file<####
######################################
if ($normProtUsage) {
	
	#if ($sharedProtNorm) {
	#	my $numTechReps=scalar @quantiOrder; # total num techRep
	#	foreach my $modProtID (keys %proteinSharing) {
	#		$forNormalization{$modProtID}=1 if scalar keys %{$proteinSharing{$modProtID}} == $numTechReps;
	#	}
	#}
	
	if ($software eq 'MSstats') {
		open(NORM,">$dataDir/globalStandards.txt"); # no header
	}
	else {
		open(NORM,">$dataDir/normProtein.txt");
		print NORM "proteinId\n"; # header
	}
	if (!$isModifQuantif || $quantifParameters{'DB'}{'PEPTIDES_NORM'}) { # || $sharedProtNorm { # protID or modProtID in $sharedProtNorm
		foreach my $protID (keys %forNormalization) {
			print NORM "$protID\n";
		}
	}
	else { # sites: convert protID from custom List to all matching modProtIDs
		foreach my $protID (keys %forNormalization) {
			if ($protein2Sites{$protID}) { # prot has quantified sites
				foreach my $modProtID (keys %{$protein2Sites{$protID}}) {print NORM "$modProtID\n";}
			}
		}
	}
	close NORM;
	
	# Error handling
	my $numLines=`wc -l $dataDir/normProtein.txt`;
	chomp $numLines;
	$numLines--;
	if ($numLines < 1) {
		my $item=($isModifQuantif)? 'sites' : 'proteins';
		warn "ERROR: No $item found for normalization!\n";
		&checkForErrors($dbh);
	}
	#elsif ($sharedProtNorm) { # store num shared prot in QUANTIF_ANNOT
	#	$dbh->do("UPDATE QUANTIFICATION SET QUANTIF_ANNOT=REPLACE(QUANTIF_ANNOT,';shared',';shared;$numLines') WHERE ID_QUANTIFICATION=$quantifID");
	#	$dbh->commit;
	#}
}
# sleep 5;
# close DEBUG;
# die "Test table.txt is complete!!!";

####> Desactivate FDR option if the number of proteins sent to R is below a threshold : 1 (2 proteins needed)
#my $protNumber=scalar keys %prot2quantify; #%proteins;
#if ($quantifParameters{'R'}{'pAdj'}[0] eq 'TRUE' && $protNumber < $minProtNumber) {
#if ($quantifParameters{'R'}{'pAdj'} && $quantifParameters{'R'}{'pAdj'}[0] eq 'TRUE' && $protNumber == 1) { # only for old algo
#	@{$quantifParameters{'R'}{'pAdj'}}=('FALSE');
#	delete($quantifParameters{'R'}{'pAdj.method'});
#	my $newQuantiAnnotation='';
#	foreach my $annoStrg (split (/::/,$quantiAnnotation)) {
#		if ($annoStrg =~ /^FDR\_/ ) {
#			$newQuantiAnnotation.='FDR_CONTROL=FALSE::' if ($annoStrg =~ /^FDR_CONTROL/);
#		}
#		else{
#			$newQuantiAnnotation.=$annoStrg."::";
#		}
#	}
#	chop($newQuantiAnnotation);
#	chop($newQuantiAnnotation);
#	$dbh->do("UPDATE QUANTIFICATION SET QUANTIF_ANNOT='$newQuantiAnnotation' WHERE ID_QUANTIFICATION=$quantifID");
#	$dbh->commit;
#}

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
#undef %proteinInConds;
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
undef %protein2Sites;
#undef %proteinSharing;

# system "echo END > $dataDir/endData.flag"; # Data table generated without error

##################>ENDING SCRIPT HERE IF ABUNDANCE<#######################
if ($ratioType eq 'None') {
	if ($countIdentifiedPep && scalar keys %isGhostPeptide) { # Write ghost peptides to file for later
		open (GHOST,">$dataDir/ghosts.txt");
		foreach my $pepIdCode (keys %isGhostPeptide) {
			print GHOST $pepIdCode,"\n";
		}
		close GHOST;
	}
	exit;
}


								############################################
			#############################> RUNNING R SCRIPTS <##########################
								############################################

if ($software eq 'MSstats') { # MSstats
	&launchMSstatsRscript;
}
else { # myProMS
	&launchRscript('PRIMARY');
}
#&checkRsuccess;

#die "TEST complete"; # TEST !!!

									############################################
			#############################> PARSING QUANTIFICATION RESULTS <########################## ONLY Super/SimpleRatio from now on
									############################################

open(FILESTAT,">>$fileStat");
print FILESTAT "3/3 Processing results\n";
close FILESTAT;
#$dbh->disconnect; exit; # DEBUG


if ($software eq 'MSstats') { # MSstats
	&processMSstatsResults;
	exit;
}

### myProMS >>>
my $numQuantifRatios=scalar @ratioTarget; # =numStates*0.5*(numStates-1)
my $adjStrg=($ratioType=~/S\w+Ratio/ || ($quantifParameters{'DB'}{'FDR_CONTROL'} && $quantifParameters{'DB'}{'FDR_CONTROL'}[0] eq 'TRUE'))? '_ADJ' : '';

####>Fetching list of quantification parameters<####
$dbh=&promsConfig::dbConnect('no_user'); # reconnect
my $dbhLite=&promsQuantif::dbCreateProteinQuantification($quantifID,$runDir); # SQLite

my %quantifParamIDs;
my $sthQP=$dbh->prepare("SELECT QP.ID_QUANTIF_PARAMETER,QP.CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD QM,QUANTIFICATION_PARAMETER QP WHERE Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND QM.ID_QUANTIFICATION_METHOD=QP.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=$quantifID");
$sthQP->execute;
while (my ($paramID,$code)=$sthQP->fetchrow_array) {
    $quantifParamIDs{$code}=$paramID;
}
$sthQP->finish;

##my $sthInsProt=$dbh->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_QUANTIFICATION,ID_PROTEIN,ID_QUANTIF_PARAMETER,QUANTIF_VALUE,TARGET_POS) VALUES ($quantifID,?,?,?,?)"); # PK autoincrement
##my ($sthInsModRes,$sthInsProtRes);
##if ($isModifQuantif) {
##	$sthInsModRes=($quantifiedModifID)? $dbh->prepare("INSERT INTO MODIFIED_RESIDUE (ID_QUANTIFICATION,RESIDUE,POSITION) VALUES ($quantifID,?,?)")
##									  : $dbh->prepare("INSERT INTO MODIFIED_RESIDUE (ID_QUANTIFICATION,RESIDUE,POSITION,MODIF_RANK) VALUES ($quantifID,?,?,?)");
##	$sthInsProtRes=$dbh->prepare("INSERT INTO PROTQUANTIF_MODRES (ID_MODIF_RES,ID_PROT_QUANTIF) VALUES (?,?)");
##}
my $sthInsProt=$dbhLite->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_PROTEIN,SITE_CODE,ID_QUANTIF_PARAMETER,QUANTIF_VALUE,TARGET_POS) VALUES (?,?,?,?,?)"); # PK autoincrement

my (%protQuantified,%lostProteins); #%modResiduesID,

###################################>(Super/Simple)Ratio (ML->AS->IB scripts)<#############################################
if ($ratioType=~/S\w+Ratio/) {
	
	open(FILESTAT,">>$fileStat");
	print FILESTAT "3/$numSteps Processing results\n";
	close FILESTAT;
	
	################################################################
	####>Recording protein ratios (& p-values, Std. dev if any)<####
	################################################################

	####>Recording infinite ratios<#### because prot might not be kept in quantif result files
	##if ($isModifQuantif) { # modif quantification
	##	foreach my $modProtID (sort{&promsMod::sortSmart($a,$b)} keys %infiniteRatioProteins) {
	##		my ($protID,@modQuantifs)=split(/[-&]/,$modProtID);
	##		my @allModResidues;
	##		foreach my $modQuantif (@modQuantifs) {
	##			my @modResidues;
	##			if ($quantifiedModifID) {@modResidues=split('\.',$modQuantif);}
	##			else {
	##				(my $modifRank,@modResidues)=split(/[#.]/,$modQuantif);
	##				foreach my $i (0..$#modResidues) {$modResidues[$i]=$modifRank.'#'.$modResidues[$i];} # Convert '1#V32.K44' into ['1#V32','1#K44']
	##			}
	##			&promsQuantif::insertModifiedResidues($protID,\@modResidues,\%modResiduesID,$dbh,$sthInsModRes); # $modifRank undef for single-modif quantif
	##			push @allModResidues,@modResidues;
	##		}
	##		foreach my $ratioPos (sort{$a<=>$b} keys %{$infiniteRatioProteins{$modProtID}}){
	##			$sthInsProt->execute($protID,$quantifParamIDs{'RATIO'},$infiniteRatioProteins{$modProtID}{$ratioPos},$ratioPos);
	##			&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF'));
	##			$protQuantified{$modProtID}{$ratioPos}=1;
	##		}
	##	}
	##}
	##else { # whole protein quantif
	##	foreach my $protID (sort{$a<=>$b} keys %infiniteRatioProteins) {
	##		foreach my $ratioPos (sort{$a<=>$b} keys %{$infiniteRatioProteins{$protID}}){
	##			$sthInsProt->execute($protID,$quantifParamIDs{'RATIO'},$infiniteRatioProteins{$protID}{$ratioPos},$ratioPos);
	##			$protQuantified{$protID}{$ratioPos}=1;
	##		}
	##	}
	##}
	foreach my $modProtID (sort{&promsMod::sortSmart($a,$b)} keys %infiniteRatioProteins) {
		my ($protID,$modStrg)=split(/-/,$modProtID);
		foreach my $ratioPos (sort{$a<=>$b} keys %{$infiniteRatioProteins{$modProtID}}){
			$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'RATIO'},$infiniteRatioProteins{$modProtID}{$ratioPos},$ratioPos);
			$protQuantified{$modProtID}{$ratioPos}=1;
		}
	}

	#unlink "$runDir/Rplots.pdf" if -e "$runDir/Rplots.pdf"; # temporary ?!!!!!!!!!!!!!!!!!!!!!!!!!!

	###>Generating ratio codes used by R<###
	my (%measurePos2Code,%measureCode2RatioPos,%state2RatioPos); # a state can be involved in multiple ratios #,%experiment2RatioPos
	my %state2StatePos;
	my $rPos=0;
	foreach my $s1 (1..$numStates-1) { # Primary ratios
		foreach my $s2 ($s1+1..$numStates) {
			$rPos++;
			$measurePos2Code{$rPos}="State$s2.State$s1"; # for ResultsDAProt.txt
			if ($algoType =~ /PEP_INTENSITY|TDA|DIA/) {
				push @{$measureCode2RatioPos{"State$s1"}},$rPos; # for resultsPep.txt (multiple rPos: when is ref & when is test state!)
				push @{$measureCode2RatioPos{"State$s2"}},$rPos;
			}
			else {@{$measureCode2RatioPos{"State$s2.State$s1"}}=($rPos);} # for resultsPep.txt
			push @{$state2RatioPos{"State$s1"}},$rPos if ($ratioType eq 'SimpleRatio' || $s1 > 1); # for inf ratios in table.txt. If SuperRatio/PEP_RATIO: Skip primary ratios for State1 (recorded in %anaState1RatioPos because State1 ratios are analysis-dependant)
			push @{$state2RatioPos{"State$s2"}},$rPos;
		}
		last if ($quantifParameters{'DB'}{'SINGLE_REF'} || $ratioType eq 'SuperRatio');
	}
	if ($ratioType eq 'SuperRatio') { # 2ndary ratios
		foreach my $s1 (2..$numStates-1) { # skip State1 -> 2ndary ratios rely on Experiments (B/A,C/A,D/A,...,C/B,D/B,...,C/D,...)
			foreach my $s2 ($s1+1..$numStates) {
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
	#elsif ($algoType ne 'PEP_RATIO') { # also record state mean former =~ /PEP_INTENSITY|TDA|DIA/ 27/04/20
		foreach my $s (1..$numStates) {
			$rPos++;
			my $state='State'.$s;
			$measurePos2Code{$rPos}=$state if $algoType ne 'PEP_RATIO'; # for ResultsDAProt.txt
$state2StatePos{$state}=$rPos; # for num true pep in State
		}
	#}

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
		my ($protID,$modStrg)=split(/-/,$protID0);
		
		# my $protID=$protID0; # default
		# my @allModResidues;
		# if ($isModifQuantif) { # modif quantification
		# 	next if $protID0 !~ /-/; # must be a mod prot
		# 	($protID,my @modQuantifs)=split(/[-&]/,$protID0);
		# 	foreach my $modQuantif (@modQuantifs) {
		# 		my @modResidues;
		# 		if ($quantifiedModifID) {@modResidues=split('\.',$modQuantif);}
		# 		else {
		# 			(my $modifRank,@modResidues)=split(/[#.]/,$modQuantif);
		# 			foreach my $i (0..$#modResidues) {$modResidues[$i]=$modifRank.'#'.$modResidues[$i];} # Convert '1#V32.K44' into ['1#V32','1#K44']
		# 		}
		# 		&promsQuantif::insertModifiedResidues($protID,\@modResidues,\%modResiduesID,$dbh,$sthInsModRes); # <-- %modResiduesID should be full at <PROT1>
		# 		push @allModResidues,@modResidues;
		# 	}
		# }
		
		###>Conditions ratio data
		foreach my $targetPos (sort{$a<=>$b} keys %measurePos2Code) { #(1..$numQuantifRatios) {
			next if ($infiniteRatioProteins{$protID0} && $infiniteRatioProteins{$protID0}{$targetPos});
			#next if ($protQuantified{$protID0} && $protQuantified{$protID0}{$ratioPos});
			#>Fold change
			if ($targetPos <= $numQuantifRatios) { # foldChange
				my $log2FC=$lineData[$protHeader2Idx{'Log2FC_'.$measurePos2Code{$targetPos}}];
				next if $log2FC=~/NA/i; # NA or NaN
				my $fcValue=($log2FC eq '-Inf')? 0.001 : ($log2FC eq 'Inf')? 1000 : 2**$log2FC; #exp($log2FC*$log2);
				$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'RATIO'},$fcValue,$targetPos);
				##&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
				$protQuantified{$protID0}{$targetPos}=1;
			}
			else { # mean state (PEP_INTENSITY)
next unless defined $protHeader2Idx{'Log2Mean_'.$measurePos2Code{$targetPos}}; # temp AS should had this again soon (PP 19/10/17)
				my $log2Mean=$lineData[$protHeader2Idx{'Log2Mean_'.$measurePos2Code{$targetPos}}];
				next if $log2Mean=~/NA/i; # NA or NaN
				$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'MEAN_STATE'},2**$log2Mean,$targetPos); #  exp($log2Mean*$log2)
				##&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif; # unlikely in LABELFREE...
			}
			#>Associated parameters
			foreach my $refParam (@RDbCodeList) {
				my ($paramR,$paramDB)=@{$refParam};
				my $paramValue=$lineData[$protHeader2Idx{$paramR.'_'.$measurePos2Code{$targetPos}}];
				next if $paramValue=~/NA/i; # NA or NaN
				$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{$paramDB},$paramValue,$targetPos);
				##&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
			}
		} # next ratioPos

	} # next line
	close PROT;

	####>Computing peptide counts (NUM_PEP_TOTAL, NUM_PEP_USED & DIST_PEP_USED)<####
	my (%numPeptideSets,%usedPeptideSets,%distinctPeptideSets); # Set <=> num label channels (1 for PEP_INTENSITY)
my (%usedPeptideReplicate,%distinctPeptideReplicate,%processedReplicates);
my %usedTruePeptideState;
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
		my $condCode=$lineData[$pepHeader2Idx{'Condition'}]; # StateX
		next if ($algoType eq 'PEP_RATIO' && !$measureCode2RatioPos{$condCode}); # ratio not used (SINGLE_REF) StateY.StateX (not applicable for LABELFREE: StateX\nStateY)
		my $seqVarMod=(split(/_/,$lineData[$pepHeader2Idx{'Peptide'}]))[0];
		#my $numPepSets=scalar (split(';',$pepIdStrg));
		my @pepIdList=split(/[\.\+\=]/,$lineData[$pepHeader2Idx{'PeptideId'}]);
		my %numPepSets;
		##my $numPepSets=1; # default
		# TODO: True for SILAC & other Isotopes (not LF or isobar) --->
		if ($labeling eq 'SILAC' && $algoType eq 'PEP_RATIO') { # former $labeling eq 'SILAC'
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
next unless $usedPeptideSets{$protID0}; # only +/-INF ratio(s) for protein => Skip 'NORM' (PP 18/09/20)
		##>Num true peptides used in State
		if ($countIdentifiedPep) { # only for 'PEP_INTENSITY' for now (not for PEP_RATIO: Value would be specific to ratio, not to State!)
			my $stateTgtPos=$state2StatePos{$condCode};
			%{$usedTruePeptideState{$protID0}{$stateTgtPos}}=() if (!$usedTruePeptideState{$protID0} || !$usedTruePeptideState{$protID0}{$stateTgtPos}); # defined because scalar=0 is recorded
			foreach my $pepIdCode (keys %numPepSets) {
				$usedTruePeptideState{$protID0}{$stateTgtPos}{$pepIdCode}=1 unless $isGhostPeptide{$pepIdCode};
			}
		}
		##>Peptides in Replicates
		foreach my $state (split(/\./,$condCode)) {	# eg "StateY.StateX" for PEP_RATIO or "StateX" for PEP_INTENSITY/TDA
			my $replicateTgtPos=$replicateTargetPos{ $state.'.'.$lineData[$pepHeader2Idx{'replicate'}].'.'.$lineData[$pepHeader2Idx{'repTech'}] };
			$distinctPeptideReplicate{$protID0}{$replicateTgtPos}{$seqVarMod}=1; # DIST_PEP_USED for replicates
			foreach my $pepIdCode (keys %numPepSets) {
				$usedPeptideReplicate{$protID0}{$replicateTgtPos}{$pepIdCode}=1;
			}
			$processedReplicates{"$protID0:$replicateTgtPos"}=1;
		}		
	}
	close PEP;

	##>Step 2: Read table.txt for NUM_PEP_TOTAL & inf ratios NUM_PEP_USED, DIST_PEP_USED
	####my %expCode2RatioPos; # used for SuperRatio only!!!
	####if ($Rdesign eq 'SUPERRATIO') {
	####	foreach my $primRatioPos (1..$numStates-1) {
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
		if ($algoType eq 'PEP_RATIO' && $labeling eq 'SILAC') {
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
		if ($ratioType eq 'SuperRatio' || $algoType eq 'PEP_RATIO') {
			if ($state eq 'State1') {
				my ($anaID)=$qSet=~/^A(\d+)/;
				#my @ratioList=@{$anaState1RatioPos{$anaID}}; # primary ratios using anaID
				#push @ratioList,@{$state2RatioPos{'State1'}}; # all 2ndary ratios
				#>Add pos of ratios using ref State1 and associated test states in context of THIS anaID
				my %myRatioPos;
				foreach my $rPos1 (@{$anaState1RatioPos{$anaID}}) {
					$myRatioPos{$rPos1}=1; # primary ratioPos
					#next if $quantifParameters{'DB'}{'SINGLE_REF'};
					next if $ratioType ne 'SuperRatio';
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
next unless $stateInInfRatio; # Computed from ResultsPep.txt for normal ratios (PP 18/09/20)

		##>Num true peptides used in State
		if ($countIdentifiedPep) { # only for 'PEP_INTENSITY' for now (not for PEP_RATIO: Value would be specific to ratio, not to State!)
			my $stateTgtPos=$state2StatePos{$state}; # Use -$stateTgtPos for data from INF ratio due to possible NORM & INF contexts
			%{$usedTruePeptideState{$protID0}{-$stateTgtPos}}=() if (!$usedTruePeptideState{$protID0} || !$usedTruePeptideState{$protID0}{-$stateTgtPos}); # defined because scalar=0 is recorded
			foreach my $pepIdCode (keys %numPepSets) {
				$usedTruePeptideState{$protID0}{-$stateTgtPos}{$pepIdCode}=1 unless $isGhostPeptide{$pepIdCode};
			}
		}

		##>Peptides in Replicates
		my $replicateTgtPos=$replicateTargetPos{"$state.$bioRep.$techRep"};
		#if (!$processedReplicates{"$protID0:$replicateTgtPos"} || $stateInInfRatio) { # WARNING: same state can be found in normal and inf ratios: num pep per replicates calculated will be based on normal ratio (from resultsPep.txt) => without excluded pep
			#my $ratioTypeKey=($stateInInfRatio)? 'INF' : 'NORM';
			$distinctPeptideReplicate{$protID0}{-$replicateTgtPos}{$seqVarMod}=1; # DIST_PEP_USED for replicates
			foreach my $pepIdCode (keys %numPepSets) {
				$usedPeptideReplicate{$protID0}{-$replicateTgtPos}{$pepIdCode}=1; # $peptideData
			}
		#}
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
		my ($protID,$modStrg)=split(/-/,$protID0);
		# my $protID=$protID0; # default
		# my @allModResidues;
		# if ($isModifQuantif) {
		# 	($protID,my @modQuantifs)=split(/[-&]/,$protID0);
		# 	foreach my $modQuantif (@modQuantifs) {
		# 		my @modResidues;
		# 		if ($quantifiedModifID) {@modResidues=split('\.',$modQuantif);}
		# 		else {
		# 			(my $modifRank,@modResidues)=split(/[#.]/,$modQuantif);
		# 			foreach my $i (0..$#modResidues) {$modResidues[$i]=$modifRank.'#'.$modResidues[$i];} # Convert '1#V32.K44' into ['1#V32','1#K44']
		# 		}
		# 		push @allModResidues,@modResidues;
		# 	}
		# }
		# Count for ratio
		foreach my $ratioPos (keys %{$distinctPeptideSets{$protID0}}) {
			next if (!$protQuantified{$protID0} || !$protQuantified{$protID0}{$ratioPos});
			$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'DIST_PEP_USED'},scalar keys %{$distinctPeptideSets{$protID0}{$ratioPos}},$ratioPos);
			##&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
		}
		next if !$protQuantified{$protID0};
		# Count for replicate
		foreach my $replicateTgtPos (sort{$a<=>$b} keys %{$distinctPeptideReplicate{$protID0}}) {
			$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'DIST_PEP_USED'},scalar keys %{$distinctPeptideReplicate{$protID0}{$replicateTgtPos}},$replicateTgtPos);
			# my $valueNorm=($distinctPeptideReplicate{$protID0}{$replicateTgtPos}{NORM})? scalar keys %{$distinctPeptideReplicate{$protID0}{$replicateTgtPos}{NORM}} : 0; # defined unless all outliers
			# my $valueInf=($distinctPeptideReplicate{$protID0}{$replicateTgtPos}{INF})? scalar keys %{$distinctPeptideReplicate{$protID0}{$replicateTgtPos}{INF}} : 0;
			# if ($valueNorm) {		
			# 	$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'DIST_PEP_USED'},$valueNorm,$replicateTgtPos);
			# 	##&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
			# }
			# if ($valueInf) { # tgtPos always < 0 for INF (PP 18/09/20)
			# 	$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'DIST_PEP_USED'},$valueInf,-$replicateTgtPos);
			# 	##&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
			# }
		}
	}
	#>NUM_PEP_USED
	foreach my $protID0 (keys %usedPeptideSets) {
		my ($protID,$modStrg)=split(/-/,$protID0);
		# my $protID=$protID0; # default
		# my @allModResidues;
		# if ($isModifQuantif) {
		# 	($protID,my @modQuantifs)=split(/[-&]/,$protID0);
		# 	foreach my $modQuantif (@modQuantifs) {
		# 		my @modResidues;
		# 		if ($quantifiedModifID) {@modResidues=split('\.',$modQuantif);}
		# 		else {
		# 			(my $modifRank,@modResidues)=split(/[#.]/,$modQuantif);
		# 			foreach my $i (0..$#modResidues) {$modResidues[$i]=$modifRank.'#'.$modResidues[$i];} # Convert '1#V32.K44' into ['1#V32','1#K44']
		# 		}
		# 		push @allModResidues,@modResidues;
		# 	}
		# }
		# Count for ratio
		foreach my $ratioPos (sort{$a<=>$b} keys %{$usedPeptideSets{$protID0}}) {
			next if (!$protQuantified{$protID0} || !$protQuantified{$protID0}{$ratioPos});
			##$sthInsProt->execute($protID,$quantifParamIDs{'NUM_PEP_USED'},$usedPeptideSets{$protID0}{$ratioPos},$ratioPos);
			$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'NUM_PEP_USED'},scalar keys %{$usedPeptideSets{$protID0}{$ratioPos}},$ratioPos);
			##&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
		}
		next if !$protQuantified{$protID0};
		
		# Count for replicate
		foreach my $replicateTgtPos (sort{$a<=>$b} keys %{$usedPeptideReplicate{$protID0}}) {
			$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'NUM_PEP_USED'}, scalar keys %{$usedPeptideReplicate{$protID0}{$replicateTgtPos}},$replicateTgtPos);
			# my $valueNorm=($usedPeptideReplicate{$protID0}{$replicateTgtPos}{NORM})? scalar keys %{$usedPeptideReplicate{$protID0}{$replicateTgtPos}{NORM}} : 0; # defined unless all outliers
			# my $valueInf=($usedPeptideReplicate{$protID0}{$replicateTgtPos}{INF})? scalar keys %{$usedPeptideReplicate{$protID0}{$replicateTgtPos}{INF}} : 0;
			# if ($valueNorm) {
			# 	$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'NUM_PEP_USED'},$valueNorm,$replicateTgtPos);
			# 	##&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
			# }
			# if ($valueInf) { # tgtPos always < 0 for INF (PP 18/09/20)
			# 	$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'NUM_PEP_USED'},$valueInf,-$replicateTgtPos);
			# 	##&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
			# }
		}
	}
	#>NUM_TRUE_USED Count for non-ghost for each State
	if ($countIdentifiedPep) {
		foreach my $protID0 (keys %usedTruePeptideState) {
			next if !$protQuantified{$protID0};
			my ($protID,$modStrg)=split(/-/,$protID0);
			# my $protID=$protID0; # default
			# my @allModResidues;
			# if ($isModifQuantif) {
			# 	($protID,my @modQuantifs)=split(/[-&]/,$protID0);
			# 	foreach my $modQuantif (@modQuantifs) {
			# 		my @modResidues;
			# 		if ($quantifiedModifID) {@modResidues=split('\.',$modQuantif);}
			# 		else {
			# 			(my $modifRank,@modResidues)=split(/[#.]/,$modQuantif);
			# 			foreach my $i (0..$#modResidues) {$modResidues[$i]=$modifRank.'#'.$modResidues[$i];} # Convert '1#V32.K44' into ['1#V32','1#K44']
			# 		}
			# 		push @allModResidues,@modResidues;
			# 	}
			# }
			foreach my $stateTgtPos (sort{$a<=>$b} keys %{$usedTruePeptideState{$protID0}}) {
				my $value=scalar keys %{$usedTruePeptideState{$protID0}{$stateTgtPos}};
				#next if $value > 5; # record only for if <= 5 (to save DB storage) => MUST store 0 values!!!
				$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'NUM_TRUE_USED'},$value,$stateTgtPos);
				##&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
			}
		}
	}
	#>NUM_PEP_TOTAL
	foreach my $protID0 (keys %numPeptideSets) {
		next if !$protQuantified{$protID0};
		my ($protID,$modStrg)=split(/-/,$protID0);
		# my $protID=$protID0; # default
		# my @allModResidues;
		# if ($isModifQuantif) {
		# 	($protID,my @modQuantifs)=split(/[-&]/,$protID0);
		# 	foreach my $modQuantif (@modQuantifs) {
		# 		my @modResidues;
		# 		if ($quantifiedModifID) {@modResidues=split('\.',$modQuantif);}
		# 		else {
		# 			(my $modifRank,@modResidues)=split(/[#.]/,$modQuantif);
		# 			foreach my $i (0..$#modResidues) {$modResidues[$i]=$modifRank.'#'.$modResidues[$i];} # Convert '1#V32.K44' into ['1#V32','1#K44']
		# 		}
		# 		push @allModResidues,@modResidues;
		# 	}
		# }
		$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'NUM_PEP_TOTAL'},scalar keys %{$numPeptideSets{$protID0}},undef); # no ratioPos
		##&promsQuantif::insertProteinModResidues(\@allModResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $isModifQuantif;
	}

}


#close DEBUG;

$sthInsProt->finish;
# if ($isModifQuantif) { # modif quantification
# 	$sthInsModRes->finish;
# 	$sthInsProtRes->finish;
# }
$dbhLite->commit;
$dbhLite->disconnect;

my $extraQuantifAnnot=($algoType ne 'PEP_RATIO')? '' : '::QSET_USED='.(scalar keys %{$qSetUsage{'USED'}}).'/'.((scalar keys %{$qSetUsage{'USED'}}) + (scalar keys %{$qSetUsage{'EXCLUDED'}})); # former  $labeling eq 'FREE' 27/04/20


&endQuantification($dbh,$projectID,\%lostProteins,\%proteinAlias,$extraQuantifAnnot);

exit;

#####################################<<<<< SUBROUTINES >>>>>#################################################

############################
####<Launching R script>####
############################
sub launchRscript { # launched in $runDir by default (peptide-level normalization +/- quantification) or in directory of protein measure to be normalized (protein-level normalization)
	my $context=$_[0] || 'PRIMARY';

	if ($context eq 'PRIMARY') { # primary (ion-level) normalization
		my $progressStrg=($ratioType eq 'None' && $normalizeIons)? "2/$numSteps Running peptide-level normalization" : ($ratioType eq 'None')? "2/$numSteps Removing peptide-level outliers" : "2/$numSteps Running quantification";
		open(FILESTAT,">>$fileStat");
		print FILESTAT "$progressStrg\n";
		close FILESTAT;

		#<Abundance with protein-only bias correction
		if ($ratioType eq 'None' && !$normalizeIons) { # overwrite R options
			@{$quantifParameters{'R'}{'normalization.method'}}=('none.none');
			delete $quantifParameters{'R'}{'common.peptides'};
		}
	}

	my ($localDataDir,$localRunDir)=($context eq 'PRIMARY')? ($dataDir,$runDir) : ("$resultDir/resultsNorm/$context/data","$resultDir/resultsNorm/$context");
	
	&promsQuantif::writeQuantifParameterFiles($localDataDir,$quantifParameters{'R'});

	open(R_SCRIPT,">$localRunDir/AnalysisDiffLimma.R");
	print R_SCRIPT qq
|
###################################################################################
# Launcher for quantification R scripts (provides path for sourcing dependencies) #
###################################################################################

filepath <- "$promsPath{R_scripts}/"
source(paste(filepath,"AnalysisDiffLimma.R",sep=""))
|;
	close R_SCRIPT;

	system "echo RUN > $localRunDir/runR.flag"; # run flag file
	system "export LANG=en_US.UTF-8; cd $localRunDir; $pathR/R CMD BATCH --no-save --no-restore AnalysisDiffLimma.R";
	sleep 5;
	&checkRsuccess("$localRunDir/AnalysisDiffLimma.Rout");
	system "rm -f $localRunDir/runR.flag; echo END > $localRunDir/endR.flag"; # end flag file (witten only if no R error, otherwise killed by &checkRsuccess)
}

sub launchMSstatsRscript { # MSstats
	open(FILESTAT,">>$fileStat");
	print FILESTAT "2/$numSteps Running quantification\n";
	close FILESTAT;

	&promsQuantif::writeQuantifParameterFiles($dataDir,$quantifParameters{'R'});
	
	open(R_SCRIPT,">$runDir/quantiSwath.R");
	print R_SCRIPT qq
|
#########################################
# Launcher for quantification R scripts #
#########################################

filepath="$promsPath{R_scripts}/"
source(paste(filepath,"quantiSwath.R",sep=""))
|;
	close R_SCRIPT;
	my $RcommandString="export LANG=en_US.UTF-8; cd $runDir; $pathR/R CMD BATCH --no-save --no-restore quantiSwath.R";
	system $RcommandString;
	sleep 5;
	
	if (-e "$resultDir/sessionInfo.txt") {
		my $versionLine=`grep MSstats_ $resultDir/sessionInfo.txt`;
		my ($version)=($versionLine=~/MSstats_(\S+)/);
		if ($version) {
			$dbh=&promsConfig::dbConnect('no_user'); # reconnect
			my $MSstatsStrg='MSstats;'.$version.'::';
			$dbh->do("UPDATE QUANTIFICATION SET QUANTIF_ANNOT=REPLACE(QUANTIF_ANNOT,'MSstats::','$MSstatsStrg') WHERE ID_QUANTIFICATION=$quantifID");
			$dbh->commit;
			$dbh->disconnect;
		}
	}
	
	&checkRsuccess("$runDir/quantiSwath.Rout");
}

sub checkRsuccess {
	####<ERROR Management>####
	my $hasError=0;
	my ($RoutFile) = @_;
	unless ($RoutFile && -e $RoutFile) {
		warn "The .Rout file was not even written ! Check that the R script is actually run.\n";
	}
	my $RoutStrg = `tail -3 $RoutFile`;
	unless ($RoutStrg=~/proc\.time\(\)/) {
		$RoutStrg=`tail -20 $RoutFile`;
		my $RerrorStrg="R script has generated an error!";
		my $inError=0;
		foreach my $line (split(/\n/,$RoutStrg)) {
			next if (!$inError && $line !~ /^Error in/ && $line !~ /^Error:/); # skip everything before "Error in..."
			$inError=1;
			$RerrorStrg.="\n$line";
		}
		warn $RerrorStrg;
		my $dbh=&promsConfig::dbConnect('no_user');
		&checkForErrors($dbh);
	}
}

#################################
####<Parsing MSstats results>####
#################################
sub processMSstatsResults {
	
	$dbh=&promsConfig::dbConnect('no_user'); # reconnect

	unless (-e "$resultDir/ResultsDAProt.txt") { # File does not exist
		warn "Missing result file !\n";
		&checkForErrors($dbh);
	}
	
	my %quantifParamIDs;
	my $sthQP=$dbh->prepare("SELECT QP.ID_QUANTIF_PARAMETER,QP.CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD QM,QUANTIFICATION_PARAMETER QP WHERE Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND QM.ID_QUANTIFICATION_METHOD=QP.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=$quantifID");
	$sthQP->execute;
	while (my ($paramID,$paramCode)=$sthQP->fetchrow_array) {
		$quantifParamIDs{$paramCode}=$paramID;
	}
	$sthQP->finish;
	
	###<Generating ratio codes used by R>###
	my (%measurePos2Code,%measureCode2RatioPos,%ratioPos2states); # a state can be involved in multiple ratios
	my $rPos=0;
	foreach my $s1 (1..$numStates-1) {
		foreach my $s2 ($s1+1..$numStates) {
			$rPos++;
			$measurePos2Code{$rPos}="State$s2-State$s1"; # for ResultsDAProt.txt
			@{$ratioPos2states{$rPos}}=("State$s1","State$s2");
		}
		last if $quantifParameters{'DB'}{'SINGLE_REF'};
	}
	
	####<Inserting data into DB>####
	my $dbhLite=&promsQuantif::dbCreateProteinQuantification($quantifID,$runDir); # SQLite
	my $sthInsProt=$dbhLite->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_PROTEIN,SITE_CODE,ID_QUANTIF_PARAMETER,QUANTIF_VALUE,TARGET_POS) VALUES (?,?,?,?,?)"); # PK autoincrement
	
	my (%protQuantified,%lostProteins,%pepList);
	my @RDbCodeList=(['pval','PVAL_ADJ'],['SD','SD_GEO']);
	my %protHeader2Idx;
	open(PROT,"$resultDir/ResultsDAProt.txt") or die "[$resultDir/ResultsDAProt.txt] $!";
	while (<PROT>) {
		$_=~s/\s*\Z//;
		my @lineData=split(/\t/,$_);
		#<Header------------------
		if ($.==1) { # 1st line of the file
			my $colIdx=0;
			foreach my $colName (@lineData) {
				$protHeader2Idx{$colName}=$colIdx;
				$colIdx++;
			}
			next;
		}
		#<Data------------------
		my $protID0=$lineData[0]; # $lineData[$protHeader2Idx{'Protein'}];
		next if !$protID0; # just to be safe
		my ($protID,$modStrg)=split(/-/,$protID0);
		next if ($quantifiedModifID && !$modStrg); # must be a mod prot

		###<Conditions ratio data
		foreach my $targetPos (sort{$a<=>$b} keys %measurePos2Code) {
			#<Fold change
			my $log2FC=$lineData[$protHeader2Idx{'log2FC_'.$measurePos2Code{$targetPos}}];
			next if $log2FC=~/NA/i; # NA or NaN
			my $fcValue=($log2FC eq '-Inf')? 0.001 : ($log2FC eq 'Inf')? 1000 : 2**$log2FC; #exp($log2FC*$log2);
			$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'RATIO'},$fcValue,$targetPos);
			$protQuantified{$protID0}{$targetPos}=1;
			
			#<Associated parameters
			foreach my $refParam (@RDbCodeList) {
				my ($paramR,$paramDB)=@{$refParam};
				my $paramValue=$lineData[$protHeader2Idx{$paramR.'_'.$measurePos2Code{$targetPos}}];
				next if $paramValue=~/NA/i; # NA or NaN
				$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{$paramDB},$paramValue,$targetPos);
			}
		} # next ratioPos

	} # next line
	close PROT;
	
	####<Computing peptide counts>####
	###<Read table.txt
	open(DATA,"$dataDir/table.txt");
	while(<DATA>) {
		next if $.==1;
		chomp;
		#ProteinName     PeptideSequence PrecursorCharge FragmentIon     ProductCharge
		my ($protID0,$peptide,$charge,$fragType,$fragCharge,$labelType,$state,$bioRep,$run,$fragArea)=split(/\t/,$_);
		if ($protQuantified{$protID0}) {
			$pepList{$protID0}{$state}{"$peptide\_$charge"}{"$fragType\_$fragCharge"}=1 if $fragArea ne 'NA';
		}
		else {$lostProteins{$protID0}=1;}
	}
	close DATA;
	
	###<Read CVoutPep.txt to recover outliers
	if (-e "$resultDir/CVoutPep.txt") {
		open(OUTLIERS,"$resultDir/CVoutPep.txt") or die $!;
		while (<OUTLIERS>) {
			next if $.==1;
			my ($protID0,$pep,$state,$value)=split(/\s/,$_);
			next unless $protQuantified{$protID0};
			my ($peptide,$pepCharge,$pepIon,$ionCharge)=split(/_/,$pep);
			delete $pepList{$protID0}{$state}{"$peptide\_$pepCharge"}{"$pepIon\_$ionCharge"};
			delete $pepList{$protID0}{$state}{"$peptide\_$pepCharge"} if scalar keys %{$pepList{$protID0}{$state}{"$peptide\_$pepCharge"}} == 0;
		}
		close OUTLIERS;
	}
	
	###<Read excluded.txt
	if (-e "$dataDir/excluded.txt") {
		open(EXCLUDED,"$dataDir/excluded.txt") or die $!;
		while (<EXCLUDED>) {
			next if $.==1;
			my @lineInfo=split(/\s/,$_);
			next unless $protQuantified{$lineInfo[0]};
			delete $pepList{$lineInfo[0]}{$lineInfo[6]}{"$lineInfo[1]\_$lineInfo[2]"}{"$lineInfo[3]\_$lineInfo[4]"};    ##$lineInfo[0] -> protID ; $lineInfo[1] -> peptide seq ; $lineInfo[2] -> peptide charge ; $lineInfo[3] -> fragment ; $lineInfo[4] -> fragment charge ; $lineInfo[6] -> state
			delete $pepList{$lineInfo[0]}{$lineInfo[6]}{"$lineInfo[1]\_$lineInfo[2]"} if scalar keys %{$pepList{$lineInfo[0]}{$lineInfo[6]}{"$lineInfo[1]\_$lineInfo[2]"}} == 0;
		}
		close EXCLUDED;
	}
	#else { # File not read (does not exist)
	#	$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID");
	#	$dbh->commit;
	#	$dbh->disconnect;
	#	die "Missing file 'excluded.txt'!\n";
	#}
	
	###<Compute peptide counts
	foreach my $protID0 (keys %pepList) {
		next if !$protQuantified{$protID0};
		my ($protID,$modStrg)=split(/-/,$protID0);
		next if ($quantifiedModifID && !$modStrg); # must be a mod prot
		my %pepUsed;
		my $numPepTotal=0;
		foreach my $ratioPos (sort {$a <=> $b} keys %ratioPos2states) {
			next if !$protQuantified{$protID0}{$ratioPos};
			my %pepDistinct;
			my $numPepUsed=0;
			foreach my $state (@{$ratioPos2states{$ratioPos}}) {
				unless ($pepUsed{$state}) { # same state can be used in different ratios
					$pepUsed{$state}=scalar keys %{$pepList{$protID0}{$state}};
					$numPepTotal+=$pepUsed{$state};
				}
				foreach my $pep (keys %{$pepList{$protID0}{$state}}) {
					my ($pepVmod,$pepCharge)=split(/_/,$pep);
					$pepDistinct{$pepVmod}=1;
				}
				$numPepUsed+=$pepUsed{$state}; # WARNING: Changed from runSWATHProtQuantification.pl to match counting done in runXICProtQuantification.pl
			}
			$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'NUM_PEP_USED'},$numPepUsed,$ratioPos);
			$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'DIST_PEP_USED'},scalar keys %pepDistinct,$ratioPos);
		}
		$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'NUM_PEP_TOTAL'},$numPepTotal,undef);
	}
	$sthInsProt->finish;
	$dbhLite->commit;

	###<Lost proteins
	my %proteinAlias;
	&getProteinAlias($dbh,\%lostProteins,\%proteinAlias);
	
	####<End quantification process>####
	&endQuantification($dbh,$projectID,\%lostProteins,\%proteinAlias);
	#my $status=($quantifParameters{'DB'}{'VISIBILITY'})? $quantifParameters{'DB'}{'VISIBILITY'}[0] : 1;
	#$dbh->do("UPDATE QUANTIFICATION SET STATUS=$status WHERE ID_QUANTIFICATION=$quantifID");
	#$dbh->commit;
	#
	#$dbh->disconnect;
	#
	#open(FILESTAT,">>$fileStat");
	#print FILESTAT "Ended ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
	#close FILESTAT;
	
	####> Moving final files
	#mkdir "$promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
	#dirmove($runDir,"$promsPath{quantification}/project_$projectID/quanti_$quantifID");
	#
	#sleep 2;

	exit;
}

sub getProteinAlias {
	my ($dbh,$refLostProteins,$refProtAlias)=@_;
	if (scalar keys %{$refLostProteins}) {
		foreach my $protID0 (keys %{$refLostProteins}) {
			my ($protID)=$protID0=~/^(\d+)/;
			$refProtAlias->{$protID}='-'; # default
		}
		
		my $sthAlias=$dbh->prepare("SELECT ID_PROTEIN,ALIAS FROM PROTEIN WHERE ID_PROTEIN IN (".join(',',keys %{$refProtAlias}).")");
		$sthAlias->execute;
		while (my ($protID,$alias)=$sthAlias->fetchrow_array) {
			$refProtAlias->{$protID}=$alias;
		}
		$sthAlias->finish;
	}
}



###############################
####<Endind Quantification>####
###############################
sub endQuantification { # Globals: $quantifID,%promsPath,%quantifParameters,$runDir,$resultDir,$fileStat
	my ($dbh,$projectID,$refLostProteins,$refProteinAlias,$extraQuantifAnnot)=@_;
	$extraQuantifAnnot='' unless $extraQuantifAnnot;
	
	my $numLostProteins=scalar keys %{$refLostProteins}; # proteins in table.txt but not quantified at all!!! Sites if site-quantif!
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

	##<LFQ lost proteins
	my $lostLFQstrg='';
	if ($ratioType eq 'None') {
		my $hasLFQ=0;
		foreach my $measType (@{$quantifParameters{'DB'}{'ABUND_MEASURES'}}) {
			if ($measType =~ /MY_LFQ/) {$hasLFQ=1; last;}
		}
		if ($hasLFQ) { # Only if Abundance with LFQ
			$lostLFQstrg='::LOST_PROTEINS_LFQ=';
			my $numLost=0;
			if (-e "$resultDir/lostProteinsLFQ.txt") {
				$numLost=`wc -l $resultDir/lostProteinsLFQ.txt`;
				$numLost=~s/^(\d+).*/$1/;
				$numLost--; # do not count header line
			}
			$lostLFQstrg.=$numLost;
		}
	}

	##<Proteins (NOT sites) lost due to shared peptides
	my $numRejectedProteins=0;
	if (-e "$dataDir/rejectedProteins.txt") {
		$numRejectedProteins=`wc -l $dataDir/rejectedProteins.txt`;
		$numRejectedProteins=~s/^(\d+).*/$1/;
		$numRejectedProteins--; # do not count header line
	}

	###> Quantification is finished.
	my $status=($quantifParameters{'DB'}{'VISIBILITY'})? $quantifParameters{'DB'}{'VISIBILITY'}[0] : 1;
	$dbh->do("UPDATE QUANTIFICATION SET STATUS=$status,QUANTIF_ANNOT=CONCAT(QUANTIF_ANNOT,'$extraQuantifAnnot\:\:REJECTED_PROTEINS=$numRejectedProteins\:\:LOST_PROTEINS=$numLostProteins$lostLFQstrg') WHERE ID_QUANTIFICATION=$quantifID");
	
	$dbh->commit;
	
	$dbh->disconnect;
	#exit; # DEBUG!!!
	
	open(FILESTAT,">>$fileStat");
	print FILESTAT "Ended ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
	close FILESTAT;
	
	####> 7th: Move all the created files in a specific directory so as to clean the directory
	#if ($ratioType=~/S\w+Ratio|None/) {
		#foreach my $scriptName ('AnalysisDiffLimma.R','FunctionLimma.R','AffectDefaultValuesLimma.R','AffectParametersLimma.R') {
		#	unlink "$runDir/$scriptName";
		#}
		#unlink "$runDir/AnalysisDiffLimma.R";
	#}

	mkdir "$promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
	dirmove($runDir,"$promsPath{quantification}/project_$projectID/quanti_$quantifID");

	sleep 2;
	
}

################################################################
####<Computing protein abundance (mean, median, sum & LFQ?)>####
################################################################
sub splitResultsPepFile {
	my ($numAbundJobs)=@_;
	
	open(FILESTAT,">>$fileStat");
	print FILESTAT "3/$numSteps Computing protein abundances\n";
	close FILESTAT;
				
	system "rm -rf $resultDir/ABUND*" if glob("$resultDir/ABUND*"); # clean up, just to be safe
	if ($numAbundJobs==1) {
		mkdir "$resultDir/ABUND1";
		system "ln -s $resultDir/resultsPep.txt $resultDir/ABUND1/resultsPep.txt";
		exit;
	}
	
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
	foreach my $job (1..$numAbundJobs) {
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
		$currentJob=1 if $currentJob > $numAbundJobs;
	}
	foreach my $job (1..$numAbundJobs) {
		close $fileHandle{$job};
	}
	
	exit;
}

sub runProteinAbundance { # all non-absolute measures: SUM,MEDIAN,...,LFQ
	my ($numAbundJobs)=@_;

	system "rm -f $resultDir/results_* errorLFQ.txt"; # in case restored from quantif w. LFQ bug (job file is copied into $resultDir)

	####<Launching parallel LFQ jobs>####
	open (ABUND_LOG,">$resultDir/abundance.log");
	print ABUND_LOG "$numAbundJobs job(s) to run\n";
	my %cluster=&promsConfig::getClusterInfo;

	if ($numAbundJobs > 1 && $cluster{'on'}) {

		#<Detect PTM quantif
		my $dataLine=`head $dataDir/table.txt | tail -1`;
		my $isModifQuantif=($dataLine=~/^\d+\-\S/)? 1 : 0;

		my $MAX_PARALLEL_QUANTIFS=$cluster{'maxJobs'};
		#$MAX_PARALLEL_QUANTIFS=25;
		my $cgiUnixDir=`pwd`;
		$cgiUnixDir=~s/\/*\s*$//; # cd is required for script to find myproms .pm files!!!!

		#my $startMem=int(1.5 + $numAbundJobs**2 / 4500); # numAbundJobs <=> num bio/tech replicates
		my $startMem=($isModifQuantif)? int(1.5 + $numAbundJobs**2 / 1500) : int(1.5 + $numAbundJobs**2 / 4500); # numAbundJobs <=> num bio/tech replicates, x3 for PTM quantif
		#my $startMem=int(1.5 + $numAbundJobs / 12.5);
		#my $startMem=32; # in Gb <-------------------------
		print ABUND_LOG $startMem."Gb RAM/job\n";
		close ABUND_LOG;
		my %allJobsParameters=(
			numCPUs=>1,
			maxHours=>48,
			pbsRunDir=>$cgiUnixDir,
			noWatch=>1 # Do not auto watch job
		);
		
		my (%launchedJobs,	# total num jobs launched (finished + running)
			%runningJobs,	# num jobs currently running
			%durationJobs,	# records job duration statistics
			%triedJobs,		# num tries for each job (max. 3/job allowed before erroring quantif)
			%restartJobs	# jobs restarted because of error or bad cluster nodes
			);
		my $numJobToRun=$numAbundJobs;
		my $numJobsFinished=0;
		MAIN_LOOP:while ($numJobToRun) {

			###<Parallel launch of up to $MAX_PARALLEL_QUANTIFS jobs
			open (ABUND_LOG,">>$resultDir/abundance.log");
			while (scalar (keys %runningJobs) < $MAX_PARALLEL_QUANTIFS) {
				
				#<Scan all job lists for job to be run (in case a previously launched job has failed and must be relaunched)
				my $currentJob=0;
				foreach my $job (1..$numAbundJobs) {
					my $jobDir="$resultDir/ABUND$job";
					if (!$launchedJobs{$job}) {
						if (-e "$jobDir/end.flag") { # in case using data from previous incomplete quantif
							$launchedJobs{$job}=1;
							$numJobToRun--;
							my @statStart=stat("$jobDir/run.flag");
							$durationJobs{START}{$job}=$statStart[9];
							my @statEnd=stat("$jobDir/end.flag");
							push @{$durationJobs{ALL}},$statEnd[9]-$durationJobs{START}{$job} if $statEnd[9]-$durationJobs{START}{$job} > 10; # job duration in seconds/ check in case restored w/o timestamp
							print ABUND_LOG "Job #$job already run. Skipping\n";
							next;
						}
						elsif (-e "$resultDir/ABUND$job/launch.flag") { # in case incomplete job launched by previously failed master quantification => clean dir
							opendir (DIR,"$resultDir/ABUND$job");
							while (my $file=readdir(DIR)) {
								next if ($file eq '.' || $file eq '..' || $file eq 'resultsPep.txt');
								unlink "$resultDir/ABUND$job/$file";
							}
							close DIR;
						}
						$currentJob=$job;
						last;
					}
				}
				last unless $currentJob; # no more jobs to launch (but some jobs could still be running)

				#<Launching a new job
				my $jobDir="$resultDir/ABUND$currentJob";
				system "echo LAUNCH > $jobDir/launch.flag"; # launch flag file
				my %jobParameters=%allJobsParameters;
				# $jobParameters{'commandBefore'}="echo RUN > $jobDir/run.flag"; # run flag file
				# $jobParameters{'commandAfter'}="echo END > $jobDir/end.flag"; # end flag file
				
				my $commandString="export LC_ALL=\"C\"; cd $cgiUnixDir; $cluster{path}{perl}/perl runXICProtQuantification.pl $quantifID $quantifDate ABUND_JOB:$currentJob"; # 2> $jobDir/errorJob.txt
				
				while (1) { # try launching job (max 3 times)
					$triedJobs{$currentJob}++;
					$restartJobs{$currentJob}++ if $triedJobs{$currentJob} > 1;
					$jobParameters{'jobName'}="myProMS_protQuant_ABUND_$quantifID.$currentJob.$triedJobs{$currentJob}";
					$jobParameters{'maxMem'}=int($startMem + 0.5 * $startMem * ($triedJobs{$currentJob}-1)).'Gb';
					
					print ABUND_LOG "Launching job #$currentJob ($jobParameters{maxMem} RAM)... ";
					my ($pbsError,$pbsErrorFile,$jobClusterID)=$cluster{'runJob'}->($jobDir,$commandString,\%jobParameters);
					if ($pbsError) {
						if ($triedJobs{$currentJob} == 5) {
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
							print ABUND_LOG "failed (x$triedJobs{$currentJob})!\n";
						}
					}
					else {
						@{$runningJobs{$currentJob}}=($jobDir,$pbsErrorFile,$jobClusterID);
						$launchedJobs{$currentJob}=1;
						print ABUND_LOG "running on $jobClusterID\n";
						last;
					}
				}
				
				my $numJobsLaunched=scalar keys %launchedJobs;
				my $numJobsRestarted=scalar keys %restartJobs;
				my $restartStrg=($numJobsRestarted)? " [$numJobsRestarted restarted]" : '';
				open(FILESTAT,">>$fileStat");
				if ($numJobToRun >= $MAX_PARALLEL_QUANTIFS) {
					print FILESTAT "3/$numSteps Computing abundances ($numJobsLaunched/$numAbundJobs tasks launched$restartStrg)\n";
				}
				else { # less than $MAX_PARALLEL_QUANTIFS jobs left: display number of remaining jobs
					print FILESTAT "3/$numSteps Computing abundances ($numJobToRun/$numAbundJobs tasks remaining$restartStrg)\n";
				}
				close FILESTAT;
				
				last if $numJobsLaunched==$numAbundJobs; # not necessary
				sleep 5;
			}
			
			###<Watching running jobs
			sleep 30;
			foreach my $job (sort{$a<=>$b} keys %runningJobs) {
				my ($jobDir,$pbsErrorFile,$jobClusterID)=@{$runningJobs{$job}};

				#<Check for error
				if (-e $pbsErrorFile) {
					my $pbsError=$cluster{'checkError'}->($pbsErrorFile);
					if ($pbsError) {
						print ABUND_LOG "**Job #$job on $jobClusterID has failed $triedJobs{$job} time(s)!**\n";
						if ($triedJobs{$job} == 5) {
							close ABUND_LOG;
							system "echo FAILED > $jobDir/failed.flag";
							unlink "$jobDir/end.flag" if -e "$jobDir/end.flag";
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
							delete $durationJobs{START}{$job} if $durationJobs{START};
							next; # job
						}
					}
				}

				if (-e "$jobDir/end.flag") { #  new job has ended without (cluster launch) error
					delete $runningJobs{$job};
					if (!-e "$jobDir/resultsLFQ.txt") { # unexpected script error
						print ABUND_LOG "**Job #$job on $jobClusterID has ended prematurely!**\n";
						#<Clean jobDir except resultsPep.txt
						opendir(DIR,$jobDir);
						while (my $file = readdir(DIR)) {
							next if ($file eq '.' || $file eq '..' || $file eq 'resultsPep.txt');
							unlink	"$jobDir/$file";
						}
						close DIR;
						delete $launchedJobs{$job};
						delete $durationJobs{START}{$job} if $durationJobs{START};
						next; # job
					}

					#<End with no error
					$numJobToRun--;
					$numJobsFinished++;
					print ABUND_LOG "Job #$job is finished\n";

					#<Computing job duration statistics
					if (!$durationJobs{START} || !$durationJobs{START}{$job}) {
						my @statStart=stat("$jobDir/run.flag");
						$durationJobs{START}{$job}=$statStart[9];
					}
					my @statEnd=stat("$jobDir/end.flag");
					push @{$durationJobs{ALL}},$statEnd[9]-$durationJobs{START}{$job}; # job duration in seconds
					my $numDurations=scalar @{$durationJobs{ALL}};
					if (scalar  keys %launchedJobs > $numJobToRun && $numDurations >= 3) { # more than half jobs finished
						my @sortedDurations=sort{$a<=>$b} @{$durationJobs{ALL}};
						my $midRank=int($numDurations/2);
						my $median=($midRank % 2)? $sortedDurations[$midRank] : ($sortedDurations[$midRank-1]+$sortedDurations[$midRank])/2;
						my $lowQuatRank=int($numDurations/4);
						my $lowerQuatLimit=($lowQuatRank % 2)? $sortedDurations[$lowQuatRank] : ($sortedDurations[$lowQuatRank-1]+$sortedDurations[$lowQuatRank])/2;
						my $upQuatRank=int(3*$numDurations/4);
						my $upperQuatLimit=($upQuatRank % 2)? $sortedDurations[$upQuatRank] : ($sortedDurations[$upQuatRank-1]+$sortedDurations[$upQuatRank])/2;
						my $outLimit=$upperQuatLimit + 1.5*($upperQuatLimit-$lowerQuatLimit); # Q3 + 1.5 * IQR = upper outlier limit
						$durationJobs{LIMIT}=($outLimit < 900)? 900: ($outLimit < 2*$median)? 2*$median : int(0.5+$outLimit); # allow at least 15 min or 2 medians
						#$durationJobs{MEDIAN_STRG}=int($median/60).'m'.($median % 60).'s';
						print ABUND_LOG "[Job duration statistics for $numDurations recorded jobs: median=$median sec, max. allowed=$durationJobs{LIMIT} sec\n";
					}

					if ($numJobToRun < $MAX_PARALLEL_QUANTIFS) { # less than $MAX_PARALLEL_QUANTIFS jobs left: display number of remaining jobs
						open(FILESTAT,">>$fileStat");
						print FILESTAT "3/$numSteps Computing abundances ($numJobToRun/$numAbundJobs tasks remaining)\n";
						close FILESTAT;
					}
				}

				elsif (-e "$jobDir/run.flag") { # job is running
					if (!$durationJobs{START} || !$durationJobs{START}{$job}) {
						my @statStart=stat("$jobDir/run.flag");
						$durationJobs{START}{$job}=$statStart[9];
					}
					if ($durationJobs{LIMIT}) { # Check if job is not taking too long (Some cluster nodes are deficient!!!)
						my $now=time; # time in sec
						if ($now-$durationJobs{START}{$job} > $durationJobs{LIMIT} * $triedJobs{$job}) { # job has been running for too long => Kill & prepare for relaunch
							my $jobClusterID=$runningJobs{$job}[2];
							print ABUND_LOG "Job #$job on $jobClusterID is taking too long (attempt #$triedJobs{$job}). Killing... ";
							$cluster{'killJob'}->([$jobClusterID]);
							if ($triedJobs{$job} == 3) {
								print ABUND_LOG "No more retry. Quantification has failed!\n";
								close ABUND_LOG;
								system "echo FAILED > $jobDir/failed.flag";
								unlink "$jobDir/end.flag" if -e "$jobDir/end.flag";
								die "Job #$job ($jobClusterID) has failed"; #error will be handled by parent process (MASTER)
							}
							opendir (DIR,$jobDir);
							while (my $file=readdir(DIR)) {
								next if ($file eq '.' || $file eq '..' || $file eq 'resultsPep.txt');
								unlink "$jobDir/$file";
							}
							close DIR;
							delete $runningJobs{$job};
							delete $launchedJobs{$job};
							delete $durationJobs{START}{$job};
							print ABUND_LOG "Successfully killed\n";
							sleep 10;
						}
					}
				}

				sleep 2;
			}
			close ABUND_LOG;
		}
	}
	else { # single job or no cluster ("local" run)
		close ABUND_LOG;
		my $MAX_PARALLEL_QUANTIFS=&promsConfig::getMaxParallelJobs;
		#$MAX_PARALLEL_QUANTIFS=4;		
		my %runningJobs;
		my $numJobToRun=$numAbundJobs;
		$SIG{CHLD} = sub { # sub called when child communicates with parent (occur during parent sleep)
			local $!; # good practice. avoids changing errno.
			while (1) { # while because multiple children could finish at same time
				my $childPid = waitpid(-1,WNOHANG); # WNOHANG: parent doesn't hang during waitpid
				#my $childPid = waitpid(-1,POSIX->WNOHANG); # WNOHANG: parent doesn't hang during waitpid
				last unless ($childPid > 0); # No more to reap.
				open (ABUND_LOG,">>$resultDir/abundance.log");
				print ABUND_LOG "Job #$runningJobs{$childPid} is finished\n";
				close ABUND_LOG;
				delete $runningJobs{$childPid}; # untrack child.
				$numJobToRun--;
			}
		};
	
		####>Looping through job list (No restart on job error!)
		my $currentJob=0;
		MAIN_LOOP:while ($currentJob < $numAbundJobs) {
	
			###>Parallel launch of up to $MAX_PARALLEL_QUANTIFS jobs
			open (ABUND_LOG,">>$resultDir/abundance.log");
			my $newJobs=0;
			while (scalar (keys %runningJobs) < $MAX_PARALLEL_QUANTIFS) {
					
				##>Also check other user/launches jobs
				my $numProc=`ps -edf | grep ABUND_JOB | grep -v grep | wc -l`;
				chomp($numProc);
				last if $numProc >= $MAX_PARALLEL_QUANTIFS;
	
				##>Ok to go
				$currentJob++;
				my $jobDir="$resultDir/ABUND$currentJob";
				if (-e "$jobDir/end.flag") { # in case using data from previous incomplete quantif
					$numJobToRun--;
					print ABUND_LOG "Job #$currentJob already run. Skipping\n";
				}
				else {
					print ABUND_LOG "Launching job #$currentJob... ";
					my $childPid = fork;
					unless ($childPid) { # child here
						close ABUND_LOG;
						#system "echo RUN > $jobDir/run.flag", # run flag file
						system "./runXICProtQuantification.pl $quantifID $quantifDate ABUND_JOB:$currentJob";					
						#system "echo END > $jobDir/end.flag", # end flag file
						exit;
					}
					$runningJobs{$childPid}=$currentJob;
					print ABUND_LOG " OK\n";
				}
				$newJobs++;
				last if $currentJob==$numAbundJobs; # no more jobs to launch
				sleep 2;
			}
			close ABUND_LOG;

			if ($newJobs) {
				open(FILESTAT,">>$fileStat");
				print FILESTAT "3/$numSteps Computing abundances ($currentJob/$numAbundJobs tasks launched)\n";
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
				print FILESTAT "3/$numSteps Computing abundances ($numJobToRun/$numAbundJobs tasks remaining)\n";
				close FILESTAT;
			}
		}
	}

	exit;
}


sub jobProteinAbundance {
	my ($job)=@_;
	my $jobDir="$resultDir/ABUND$job";
	system "echo RUN > $jobDir/run.flag"; # run flag file

	my %measureUsed;
	my ($minPepRatioLFQ,$stabilzeLargRatios)=(1,1); # default
	foreach my $measType (@{$quantifParameters{'DB'}{'ABUND_MEASURES'}}) {
		next if ($measType =~ /^aLFQ/);
		if ($measType =~ /MY_LFQ,(\d),(\d)/) {
			$measType='MY_LFQ';
			($minPepRatioLFQ,$stabilzeLargRatios)=($1,$2);
		}
		$measureUsed{$measType}=1;
	}
	
	###<Running LFQ>###
	if ($measureUsed{'MY_LFQ'}) {
		#my $minPepRatioLFQ=($quantifParameters{'DB'}{'NUM_LFQ_RATIOS'})? $quantifParameters{'DB'}{'NUM_LFQ_RATIOS'}[0] : 1;
		my $SLRflag=($stabilzeLargRatios)? '-s' : '';
		my $aggregRepFlag=($normalizeProteins)? '' : '-a';
		my $LFQcommand="$pathPython3/python3 $promsPath{python_scripts}/computeLFQ.py -i $jobDir/resultsPep.txt -o $jobDir/resultsLFQ.txt -m $minPepRatioLFQ $SLRflag $aggregRepFlag -f 2> $jobDir/commentsLFQ.txt";
		system $LFQcommand;

		##<Check LFQ results (1 failure allowed)
		my $killedLFQ;
		if (-e "$jobDir/commentsLFQ.txt") {
			$killedLFQ=`grep Killed $jobDir/commentsLFQ.txt`;
			$killedLFQ=~s/\s//g;
		}
		if (!-e "$jobDir/resultsLFQ.txt" || -z "$jobDir/resultsLFQ.txt" || $killedLFQ) { # no/empty file or killed
			unlink "$jobDir/resultsLFQ.txt" if -e "$jobDir/resultsLFQ.txt";
			unlink "$jobDir/commentsLFQ.txt" if -e "$jobDir/commentsLFQ.txt"; 
			system $LFQcommand;
		}
	}
	
	####<Loading list of ghost peptides>####
	my %isGhostPeptide;
	if ($countIdentifiedPep && -e "$dataDir/ghosts.txt") {
		open(GHOST,"$dataDir/ghosts.txt");
		while(<GHOST>) {
			chomp;
			$isGhostPeptide{$_}=1;
		}
		close GHOST;
	}
	
	####<Parsing peptide results & computing peptide counts (NUM_PEP_TOTAL, NUM_PEP_USED & DIST_PEP_USED)>####
	my (%allProteins,%peptideValues,%proteinData); # Set <=> num label channels (1 for Label-Free)
my ($multiBioRep,$multiTechRep)=(0,0);
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
		my $numTruePeptides=0;
		foreach my $pepID (split(/[\.\+\=]/,$lineData[$pepHeader2Idx{'PeptideId'}])) {
			$allProteins{$protID0}{$pepID}=1; # for NUM_PEP_TOTAL & in case Labeling
			$numPeptides++;
			$numTruePeptides++ if ($countIdentifiedPep && !$isGhostPeptide{$pepID});
		}
		next if ($lineData[$pepHeader2Idx{'out'}] ne 'NA' && $lineData[$pepHeader2Idx{'out'}] ne 'outBalanceStates'); # skip outliers but not those considered outliers for ratios
		my $state=$lineData[$pepHeader2Idx{'Condition'}];
		my ($statePos)=$state=~/State(\d+)/;
		my $bioRep=$lineData[$pepHeader2Idx{'replicate'}];
		my $techRep=$lineData[$pepHeader2Idx{'repTech'}];
		my $seqVarMod=(split(/_/,$lineData[$pepHeader2Idx{'Peptide'}]))[0];
		$proteinData{'NUM_PEP_USED'}{$protID0}{$statePos}+=$numPeptides;
		$proteinData{'DIST_PEP_USED'}{$protID0}{$statePos}{$seqVarMod}=1; # DIST_PEP_USED
		$proteinData{'NUM_TRUE_USED'}{$protID0}{$statePos}+=$numTruePeptides if $countIdentifiedPep;

		####<Normalized peptide values
		push @{$peptideValues{$protID0}{$statePos}{$bioRep}{$techRep}},2**$lineData[$pepHeader2Idx{'log2Measure'}];

		$multiTechRep=1 if scalar keys %{$peptideValues{$protID0}{$statePos}{$bioRep}} > 1;
		$multiBioRep=1 if scalar keys %{$peptideValues{$protID0}{$statePos}} > 1;
	}
	close PEP;
	
	####<Saving all proteins & number of peptides>####
	open(PROT,">$jobDir/allProteins.txt") || die $!; # contains list of all proteins (before outlier filter) to compute lost proteins later
	print PROT "ProteinID\tnumPep";
	foreach my $protID0 (keys %allProteins) {
		print PROT "\n",$protID0,"\t",scalar keys %{$allProteins{$protID0}};
	}
	close PROT;

	####<Computing protein abundance (mean, geo mean, median & sum)>####
	my $replicType=(!$normalizeProteins)? '' : ($multiBioRep || !$multiTechRep)? 'bio' : 'tech'; # in last case: only 1 bio rep anyway
	my @quantTypes=('SUM_INT','MEAN_INT','GEO_MEAN_INT','MEDIAN_INT');
	foreach my $protID0 (keys %peptideValues) {
		foreach my $statePos (keys %{$peptideValues{$protID0}}) {
			my %replicAbundance;
			foreach my $bioRep (keys %{$peptideValues{$protID0}{$statePos}}) {
				my %techAbundance;
				foreach my $techRep (keys %{$peptideValues{$protID0}{$statePos}{$bioRep}}) {
					my $numValues=scalar @{$peptideValues{$protID0}{$statePos}{$bioRep}{$techRep}};
					my $midPos=int($numValues/2) || 1; # cannot be 0
					my $midPos1=$midPos+1;
					my $mustAverage=($numValues % 2)? 0 : 1;
					my $sumValue=0;
					my $sumLogValue=0;
					my $medianValue=0;
					my $count=0;
					foreach my $pepValue (sort{$a<=>$b} @{$peptideValues{$protID0}{$statePos}{$bioRep}{$techRep}}) { # aggregate all ions in this tech rep
						$sumValue+=$pepValue;
						$sumLogValue+=log($pepValue);
						$count++;
						if ($count > $midPos1) {next;}
						elsif ($count==$midPos) {$medianValue=$pepValue;}
						elsif ($mustAverage && $count==$midPos1) {
							$medianValue=($medianValue+$pepValue)/2;
						}
					}
					%{$techAbundance{$techRep}}=(
						SUM_INT			=>	$sumValue,
						MEAN_INT		=>	$sumValue/$numValues,
						GEO_MEAN_INT	=>	exp($sumLogValue/$numValues),
						MEDIAN_INT		=>	$medianValue
					);
				}
				if ($replicType eq 'tech') { # $normalizeProteins & no bio rep & multiple tech reps => replicates are tech reps
					%replicAbundance=%techAbundance;
				}
				else { # No protein norm OR multiple bio reps OR no tech rep => replicate(s) are bio rep(s)
					##>Take MEAN of tech reps
					my $numTechReps=scalar keys %techAbundance;
					foreach my $techRep (keys %techAbundance) {
						foreach my $quantType (@quantTypes) {$replicAbundance{$bioRep}{$quantType}+=$techAbundance{$techRep}{$quantType};}
					}
					foreach my $quantType (@quantTypes) {$replicAbundance{$bioRep}{$quantType}/=$numTechReps;}
				}	
			}
			if ($replicType eq 'tech') { # Normalize prot & no bio rep & multiple tech reps
				foreach my $techRep (keys %replicAbundance) {
					foreach my $quantType (@quantTypes) {$proteinData{$quantType}{$protID0}{$statePos}{$techRep}=$replicAbundance{$techRep}{$quantType};}
				}
			}
			elsif ($replicType eq 'bio') { # Normalize prot & multiple bio reps (tech reps already averaged if multiple)
				foreach my $bioRep (keys %replicAbundance) {
					foreach my $quantType (@quantTypes) {$proteinData{$quantType}{$protID0}{$statePos}{$bioRep}=int(0.5 + $replicAbundance{$bioRep}{$quantType});} # round potential mean or tech reps
				}
			}
			else { # agragate both
				##>Take mean of bio reps
				my $numBioReps=scalar keys %replicAbundance;
				foreach my $bioRep (keys %replicAbundance) {
					foreach my $quantType (@quantTypes) {$proteinData{$quantType}{$protID0}{$statePos}+=$replicAbundance{$bioRep}{$quantType};}
				}
				foreach my $quantType (@quantTypes) {$proteinData{$quantType}{$protID0}{$statePos}=int(0.5 + ($proteinData{$quantType}{$protID0}{$statePos}/$numBioReps));}
			}
		}
	}
	
	####<Printing Abundance to file>####
	# my %optDataFiles=(	SUM_INT			=> 'resultsSUM.txt',
	# 					MEAN_INT		=> 'resultsMEAN.txt',
	# 					GEO_MEAN_INT	=> 'resultsGEO.txt',
	# 					MEDIAN_INT		=> 'resultsMEDIAN.txt',
	# 					NUM_TRUE_USED	=> 'peptidesTRUE.txt'
	# 				);
	# my %dataFiles=(	NUM_PEP_USED 	=> 'peptidesUSED.txt',
	# 				DIST_PEP_USED	=> 'peptidesDIST.txt'
	# 				); # NUM_PEP_TOTAL in allProteins.txt
	foreach my $measType (@quantTypes) {
		$dataFiles{$measType}=$optDataFiles{$measType} if $measureUsed{$measType};
	}
	# $dataFiles{'NUM_TRUE_USED'}=$optDataFiles{'NUM_TRUE_USED'} if $countIdentifiedPep;
	
	##>>>>> Old (large) file format
	# my $headerStrg="ProteinID\tState".join("\tState",1..$numStates)."\n";
	# foreach my $measType (keys %dataFiles) {
	# 	open (OUT,">$jobDir/$dataFiles{$measType}") || die $!;
	# 	print OUT $headerStrg;
	# 	foreach my $protID0 (keys %{$proteinData{$measType}}) {
	# 		print OUT $protID0;
	# 		foreach my $statePos (1..$numStates) {
	# 			print OUT "\t";
	# 			if ($proteinData{$measType}{$protID0}{$statePos}) {
	# 				print OUT ($measType eq 'DIST_PEP_USED')? scalar keys %{$proteinData{$measType}{$protID0}{$statePos}} : $proteinData{$measType}{$protID0}{$statePos};
	# 			}
	# 			else {
	# 				print OUT 'NA';
	# 			}
	# 		}
	# 		print OUT "\n";
	# 	}
	# 	close OUT;
	# }
	##>>>>> New (long) file format
	my $headerStrg=($replicType)? "ProteinID\tCondition\tReplicate\tValue\n" : "ProteinID\tCondition\tValue\n";
	foreach my $measType (keys %dataFiles) {
		open (OUT,">$jobDir/$dataFiles{$measType}") || die $!;
		my ($isPeptide,$header)=($measType=~/_PEP_|_TRUE_/)? (1,"ProteinID\tCondition\tValue\n") : (0,$headerStrg);
		print OUT $header;
		foreach my $protID0 (keys %{$proteinData{$measType}}) {
			foreach my $statePos (sort {$a<=>$b} keys %{$proteinData{$measType}{$protID0}}) {
				if ($replicType && !$isPeptide) { # 4 column format
					foreach my $replic (keys %{$proteinData{$measType}{$protID0}{$statePos}}) {
						print OUT "$protID0\tState$statePos\t$replic\t$proteinData{$measType}{$protID0}{$statePos}{$replic}\n";
					}
				}
				else { # 3 columns format (peptide counts are fully aggregated for final import)
					my $value=($measType eq 'DIST_PEP_USED')? scalar keys %{$proteinData{$measType}{$protID0}{$statePos}} : $proteinData{$measType}{$protID0}{$statePos};
					print OUT "$protID0\tState$statePos\t$value\n";
				}
			}
		}
		close OUT;
	}

	system "echo END > $jobDir/end.flag"; # end flag file
	exit;
}

sub runSecondaryNormalizations {
	my ($numAbundJobs,$refPrimaryAbundance,$refProteinAbundance,$refAllPeptides,$refLostProteins)=@_;

	####<Printing data table.txt for normalization>####
	my %fileHandle;
	mkdir "$resultDir/resultsNorm" unless -e "$resultDir/resultsNorm";
	foreach my $measType (keys %dataFiles) {
		next if $measType=~/_PEP_|_TRUE_/; # skip peptide count measures
		my $measDir="$resultDir/resultsNorm/$measType";
		next if (-e $measDir && -e "$measDir/endR.flag");  # already computed
		rmtree $measDir if -e $measDir; # clean up in case resuming
		mkdir $measDir;
		mkdir "$measDir/results";
		mkdir "$measDir/results/graph";
		my $measDataDir="$measDir/data";
		mkdir $measDataDir;
		#copy("$dataDir/param_char.txt","$measDataDir/param_char.txt"); <- done by ABUND_RSCRIPT
		#system "ln -s $dataDir/normProtein.txt $measDataDir/normProtein.txt" if -e "$dataDir/normProtein.txt"; # protein normalization file
		system "cd $measDataDir && ln -s ../../../../data/normProtein.txt $measDataDir/normProtein.txt" if -e "$dataDir/normProtein.txt"; # protein normalization file
		open($fileHandle{$measType},">$measDataDir/table.txt");
		print {$fileHandle{$measType}} "Protein_ID\tPeptide\tSample\tTreatment\tReplicate\tTechnical_Replicate\tExperiment\tQuantif_Set\tPeptide_IDs\tProtein_Name\tProtein_Validity\tValue";
	}

	my $allRowCount=0;
	my %measRowCount;	
	foreach my $protID0 (keys %{$refPrimaryAbundance}) {
		foreach my $statePos (sort{$a<=>$b} keys %{$refPrimaryAbundance->{$protID0}}) {
			foreach my $replic (keys %{$refPrimaryAbundance->{$protID0}{$statePos}}) { # repX or techRepX
				my ($bioRep,$techRep)=($replic=~/^tech/)? ('rep1',$replic) : ($replic,'techRep1');
				foreach my $measType (keys %{$refPrimaryAbundance->{$protID0}{$statePos}{$replic}}) {
					$allRowCount++;
					$measRowCount{$measType}++;
					next unless $fileHandle{$measType}; # in case resuming
					print {$fileHandle{$measType}} "\n$protID0\tPEP#$protID0\tState$statePos\tNA\t$bioRep\t$techRep\tA\tA0_Q0\t$allRowCount\tProt_$protID0\t1\t$refPrimaryAbundance->{$protID0}{$statePos}{$replic}{$measType}";
				}
			}
		}
	}

	foreach my $measType (keys %fileHandle) {
		close $fileHandle{$measType};
		if ($measRowCount{$measType} && $measRowCount{$measType} > 3) {
			system "echo END > $resultDir/resultsNorm/$measType/data/endData.flag"; # Data table generated without error
		}
		else {
			warn "ERROR: Not enough entries in $measType starting data file!";
			sleep 3;
			&checkForErrors;
		}
	}

	####<Launching R normalization for each measure>####
	open(FILESTAT,">>$fileStat");
	print FILESTAT "4/$numSteps Running protein-level normalization (2/2 Correcting biases)\n";
	close FILESTAT;

	if ($cluster{'on'}) { # Cluster
		my $cgiUnixDir=`pwd`;
		$cgiUnixDir=~s/\/*\s*$//; # cd is required for script to find myproms .pm files!!!!
		my $normMaxMem=int(1.5 + 5E-6 * (scalar keys %{$refPrimaryAbundance}) * (scalar @{$quantifParameters{'DB'}{'STATES'}}));
		my %allJobsParameters=(
			numCPUs=>1,
			maxMem=>$normMaxMem.'Gb',
			maxHours=>48,
			pbsRunDir=>$cgiUnixDir,
			noWatch=>1 # Do not auto watch job
		);
		#<Launching all norm jobs
		my %runningJobs;
		foreach my $measType (keys %fileHandle) {
			my $measDir="$resultDir/resultsNorm/$measType";
			system "echo LAUNCH > $measDir/launch.flag"; # launch flag file
			my %jobParameters=%allJobsParameters;
			# $jobParameters{'commandBefore'}="echo RUN > $measDir/runR.flag", # run flag file
			# $jobParameters{'commandAfter'}="rm -f $measDir/runR.flag; echo END > $measDir/endR.flag", # end flag file
			$jobParameters{'jobName'}="myProMS_protQuant_ABUND_$quantifID.NORM.$measType";
			my $commandString="export LC_ALL=\"C\"; cd $cgiUnixDir; $cluster{path}{perl}/perl runXICProtQuantification.pl $quantifID $quantifDate ABUND_RSCRIPT:$measType 2>> $errorFile"; # 2> $measDir/errorJob.txt
			
			my ($pbsError,$pbsErrorFile,$jobClusterID)=$cluster{'runJob'}->($measDir,$commandString,\%jobParameters);
			@{$runningJobs{$measType}}=($measDir,$pbsErrorFile,$jobClusterID);
			sleep 10;
		}

		###<Watching running jobs
		sleep 30;
		my $numJobToRun=scalar keys %runningJobs;
		while ($numJobToRun) {
			foreach my $measType (keys %runningJobs) {
				my ($measDir,$pbsErrorFile,$jobClusterID)=@{$runningJobs{$measType}};

				#<Check for cluster error
				if (-e $pbsErrorFile) {
					my $pbsError=$cluster{'checkError'}->($pbsErrorFile);
					if ($pbsError) {
						system "echo FAILED > $measDir/failed.flag";
						unlink "$measDir/endR.flag" if -e "$measDir/endR.flag";
						die "Normalization job for $measType ($jobClusterID) has failed"; #error will be handled by parent process (MASTER)
					}
				}
				if (-e "$measDir/endR.flag") { #  new job has ended without error (ABUND_RSCRIPT already checks for R success)
					delete $runningJobs{$measType};
					$numJobToRun--;
				}
			}
			sleep 30;
		}
	}
	else { # Local
		my $MAX_PARALLEL_QUANTIFS=&promsConfig::getMaxParallelJobs;		
		my %runningJobs;
		my $numJobToRun = my $numNorms = scalar keys %fileHandle;
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

		####<Looping through norm jobs
		my %measLaunched;
		my $currentJob=0;
		MAIN_LOOP:while ($currentJob < $numNorms) {
	
			##<Launching all norm jobs up to $MAX_PARALLEL_QUANTIFS jobs
			foreach my $measType (keys %fileHandle) {
				next if $measLaunched{$measType}; # already launched

				##>Also check other user/launches jobs
				my $numProc=`ps -edf | grep ABUND_RSCRIPT | grep -v grep | wc -l`;
				chomp($numProc);
				last if $numProc >= $MAX_PARALLEL_QUANTIFS;
	
				##>Ok to go
				$currentJob++;
				# my $measDir="$resultDir/resultsNorm/$measType";
				# if (-e "$measDir/endR.flag") { # in case using data from previous incomplete quantif
				# 	$numJobToRun--;
				# }
				# else {
					my $childPid = fork;
					unless ($childPid) { # child here
						system "./runXICProtQuantification.pl $quantifID $quantifDate ABUND_RSCRIPT:$measType";					
						exit;
					}
					$runningJobs{$childPid}=$currentJob;
					$measLaunched{$measType}=1;
				# }
				sleep 2;
			}
	
			###>Wait for potential ended job before next MAIN_LOOP
			sleep 60; # $SIG{CHLD} is active during sleep
	
		}
		
		###>Wait for all child processes to end
		while (scalar (keys %runningJobs) > 0) {
			sleep 60; # $SIG{CHLD} is active during sleep
		}

	}

	&checkForErrors();

	system "echo END > $resultDir/end2ndNorm.flag"; # Protein normalization generated without error
}

sub processFinalAbundanceResults {
	my ($stepNumber,$numAbundJobs,$computeAbsolute,$computeNonAbsolute,$refProteinAbundance,$refAllPeptides,$refLostProteins)=@_;
	#my %absMeasures;  # To store the absolute measures if needed

	####<Parsing abundance results>####
	open(FILESTAT,">>$fileStat");
	print FILESTAT "$stepNumber/$numSteps Processing results (1/2 Parsing files)\n";
	close FILESTAT;

	#####################>>>
	# if ($computeAbsolute) {
	# 	&parseAbsoluteQuantif($refProteinAbundance) unless $normalizeProteins; # otherwise will be done after prot normalization
	# }
	if ($computeNonAbsolute) {  # Abundance measures other than absolute amounts were computed
		if ($normalizeProteins) {
			&parseSecondaryAbundanceResults($refProteinAbundance);
		}
		else {
			&parsePrimaryAbundanceResults($numAbundJobs,$refProteinAbundance,$refAllPeptides,$refLostProteins);
		}
	}
	if ($computeAbsolute) {  # Peptides data and lost proteins not yet computed in that case, so do it
		&parseAbsoluteQuantif($refProteinAbundance);

		unless ($computeNonAbsolute) { # Tasks below not yet done if only absolute quantif
			###<Parsing peptide counts for each State>###
			foreach my $measType (keys %dataFiles) {
				unless (-e "$absDir/$dataFiles{$measType}") {
					warn "File $dataFiles{$measType} not found in $absDir";
					last;
				}
				&parseMeasureFile("$absDir/$dataFiles{$measType}", $measType, $refProteinAbundance);
			}
			###<Parsing peptide counts for whole dataset before outlier filtering & lost proteins>###
			&getAllProtInfo("$absDir/allProteins.txt", $refProteinAbundance, $refAllPeptides, $refLostProteins);
		}
	}
	if (!$computeNonAbsolute && !$computeAbsolute) {
		warn "No abundance measure selected. How did it get this far ?";
		exit;
	}
	#####################<<<
	# # my %optDataFiles=(	SUM_INT			=> 'resultsSUM.txt',
	# # 					MEAN_INT		=> 'resultsMEAN.txt',
	# # 					GEO_MEAN_INT	=> 'resultsGEO.txt',
	# # 					MEDIAN_INT		=> 'resultsMEDIAN.txt',
	# # 					NUM_TRUE_USED	=> 'peptidesTRUE.txt',
	# # 					MY_LFQ			=> 'resultsLFQ.txt'
	# # 				 );
	# # my %dataFiles=(	NUM_PEP_USED 	=> 'peptidesUSED.txt',
	# # 				DIST_PEP_USED	=> 'peptidesDIST.txt'
	# # 				);
	# foreach my $measType (@{$quantifParameters{'DB'}{'ABUND_MEASURES'}}) {
	# 	if ($measType =~ /^aLFQ_/) {  # aLFQ_iBAQ or aLFQ_TOP
	# 		# foreach my $absMeas (split(';', $measType)) {  # aLFQ_iBAQ,abs_model,normalization;MOL_PERCENT;MOL_CONC;...
	# 		# 	if ($absMeas =~  /^(aLFQ_[^,]+),/) {
	# 		# 		$absMeasures{$1} = 1;
	# 		# 	} else {
	# 		# 		$absMeasures{$absMeas} = 1;
	# 		# 	}
	# 		# }
	# 	} else {
	# 		$measType =~ s/,.+//; # MY_LFQ,x,x
	# 		$dataFiles{$measType} = $optDataFiles{$measType};
	# 	}
	# }
	# # $dataFiles{'NUM_TRUE_USED'}=$optDataFiles{'NUM_TRUE_USED'} if $countIdentifiedPep;

	# my (%proteinAbundance,%allPeptides,%lostProteins);

	# if ($computeAbsolute) {
	# 	&parseAbsoluteQuantif(\%proteinAbundance);
	# }

	# if ($computeNonAbsolute) {  # Abundance measures other than absolute amounts are computed
	# 	JOB:foreach my $job (1..$numAbundJobs) {
	# 		my $jobDir="$resultDir/ABUND$job";

	# 		###<Combining data after outlier filtering>###
	# 		foreach my $measType (keys %dataFiles) {
	# 			unless (-e "$jobDir/$dataFiles{$measType}") {
	# 				if ($measType eq 'MY_LFQ') { # Do not fail if no LFQ data
	# 					open(BAD_LFQ,">>$resultDir/errorLFQ.txt");
	# 					print BAD_LFQ "File $dataFiles{$measType} not found for job #$job\n";
	# 					close BAD_LFQ;
	# 					copy("$jobDir/resultsPep.txt","$resultDir/resultsPep_$job.txt");
	# 					next;
	# 				}
	# 				warn "File $dataFiles{$measType} not found for job #$job";
	# 				last JOB;
	# 			}
	# 			&parseMeasureFile("$jobDir/$dataFiles{$measType}", $measType, \%proteinAbundance);
	# 		}
	# 		###<Combining data before outlier filtering & lost proteins>###
	# 		&getAllProtInfo("$jobDir/allProteins.txt", \%proteinAbundance, \%allPeptides, \%lostProteins);
	# 	}
	# } elsif ($computeAbsolute) {  # Peptides data and lost proteins not yet computed in that case, so do it
	# 	###<Combining data after outlier filtering>###
	# 		foreach my $measType (keys %dataFiles) {
	# 			unless (-e "$absDir/$dataFiles{$measType}") {
	# 				warn "File $dataFiles{$measType} not found in $absDir";
	# 				last;
	# 			}
	# 			&parseMeasureFile("$absDir/$dataFiles{$measType}", $measType, \%proteinAbundance);
	# 		}
	# 		###<Combining data before outlier filtering & lost proteins>###
	# 		&getAllProtInfo("$absDir/allProteins.txt", \%proteinAbundance, \%allPeptides, \%lostProteins);
	# } else {
	# 	warn "No abundance measure selected. How did it get this far ?";
	# 	exit;
	# }
	#################################<<<

}

sub importAbundanceInDB {
	my ($stepNumber,$numAbundJobs, $computeAbsolute, $refProteinAbundance,$refAllPeptides,$refLostProteins)=@_;

	my @allMeasures;  # Store all measures from dataFiles AND absolute measures
	push @allMeasures, (keys %dataFiles);
	push @allMeasures, (keys %absMeasures);

	####<Recording data in DB>####
	my $dbh=&promsConfig::dbConnect('no_user');
	&checkForErrors($dbh);
	open(FILESTAT,">>$fileStat");
	print FILESTAT "$stepNumber/$numSteps Processing results (2/2 Storing results)\n";
	close FILESTAT;
	
	my ($quantifiedModifID,$quantifModStrg)=$dbh->selectrow_array("SELECT Q.ID_MODIFICATION,GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ',')
																					 FROM QUANTIFICATION Q
																					 LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
																					 WHERE Q.ID_QUANTIFICATION=$quantifID
																					 GROUP BY Q.ID_QUANTIFICATION");
	my $isModifQuantif=0;
	if ($quantifiedModifID || $quantifModStrg) {
		$isModifQuantif=1;
		$quantifiedModifID=$quantifModStrg if (!$quantifiedModifID && $quantifModStrg !~ /,/); # single PTM: back compatibility in case QUANTIFICATION.ID_MODIFICATION no longer used
	}
	
	####<Fetching list of quantification parameters>####
	my %quantifParamIDs;
	my $sthQP=$dbh->prepare("SELECT QP.ID_QUANTIF_PARAMETER,QP.CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD QM,QUANTIFICATION_PARAMETER QP WHERE Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND QM.ID_QUANTIFICATION_METHOD=QP.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=$quantifID");
	$sthQP->execute;
	while (my ($paramID,$code)=$sthQP->fetchrow_array) {
		$quantifParamIDs{$code}=$paramID;
	}
	$sthQP->finish;
	
	my $dbhLite=&promsQuantif::dbCreateProteinQuantification($quantifID,$runDir); # SQLite
	my $sthInsProt=$dbhLite->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_PROTEIN,SITE_CODE,ID_QUANTIF_PARAMETER,QUANTIF_VALUE,TARGET_POS) VALUES (?,?,?,?,?)"); # PK autoincrement

	my $numEntities=scalar keys %{$refProteinAbundance};
	my (%numMissingValues,%protLostLFQ);
	foreach my $paramCode (@allMeasures) {
		next if ($paramCode=~/_(USED|PEP)/  || $paramCode=~/^CI_(INF|SUP)/ || $paramCode=~/^GEO_(SD|CV)/); # skip peptide and variability measures
		foreach my $statePos (1..$numStates) {
			$numMissingValues{$paramCode}{$statePos}=0;
		}
	}
	my $pcCount=int($numEntities/20); # 5%
	my $pcEnt=0;
	my $countEnt=0;
	foreach my $protID0 (keys %{$refProteinAbundance}) {
		$countEnt++;
		if ($numEntities >= 200 && $countEnt==$pcCount) {
			$pcEnt+=5; # 5%
			$countEnt=0;
			open(FILESTAT,">>$fileStat");
			print FILESTAT "4/$numSteps Processing results (2/2 Storing results [$pcEnt%])\n";
			close FILESTAT;
		}
		my $numMissingLFQ=0;
		my ($protID,$modStrg)=split(/-/,$protID0);
		#foreach my $statePos (sort{$a<=>$b} keys %{$refProteinAbundance->{$protID0}}) { # }
		foreach my $statePos (1..$numStates) {
			#foreach my $paramCode (keys %{$refProteinAbundance->{$protID0}{$statePos}}) { # also includes (NUM/DIST)_(PEP/TRUE)_USED
			foreach my $paramCode (@allMeasures) {
				if (defined $refProteinAbundance->{$protID0}{$statePos}{$paramCode}) {
					$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{$paramCode},$refProteinAbundance->{$protID0}{$statePos}{$paramCode},$statePos);
				}
				elsif ($paramCode !~ /_(USED|PEP)/  && $paramCode !~ /^CI_(INF|SUP)/ && $paramCode !~ /^GEO_(SD|CV)/) {
					$numMissingValues{$paramCode}{$statePos}++;
					$numMissingLFQ++ if $paramCode eq 'MY_LFQ'; # For LFQ, prot can be missing in all states -> do not count as missing values but as missing prot
				}
			}
		}
		if ($numMissingLFQ==$numStates) { # fully missing in LFQ
			$protLostLFQ{$protID0}=1;
			foreach my $statePos (keys %{$numMissingValues{'MY_LFQ'}}) {$numMissingValues{'MY_LFQ'}{$statePos}--;}
		}
		#>All peptides
		$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'NUM_PEP_TOTAL'},$refAllPeptides->{$protID0},undef);
	}
	$dbhLite->commit;
	$dbhLite->disconnect;
	
	###<% NA protein/site-level
	my $numProteins=scalar keys %{$refProteinAbundance};
	foreach my $paramCode (keys %numMissingValues) {
		my ($numTotProteins,$fileName)=($paramCode eq 'MY_LFQ')? ($numProteins - (scalar keys %protLostLFQ),'percentageNAProtLFQ.txt') : (defined($dataFiles{$paramCode}))? ($numProteins,'percentageNAProt.txt') : ($paramCode =~ /^MOL/)? ($numProteins,'percentageNAProtMol.txt') : ($numProteins,'percentageNAProtMass.txt');
		unless (-s "$resultDir/$fileName") {
			open (MISS,">$resultDir/$fileName");
			print MISS "State\tPC_NA\n";
			foreach my $statePos (1..$numStates) {
				print MISS "State$statePos\t",100*$numMissingValues{$paramCode}{$statePos}/$numTotProteins,"\n";
			}
			close MISS;
		}
	}

	###<LFQ lost proteins
	if (scalar keys %protLostLFQ) {
		my %proteinAliasLFQ;
		&getProteinAlias($dbh,\%protLostLFQ,\%proteinAliasLFQ);
		open (LOST,">$resultDir/lostProteinsLFQ.txt");
		print LOST "Protein_ID\tProtein_Name";
		foreach my $modProtID (sort{&promsMod::sortSmart($a,$b)} keys %protLostLFQ) {
			my ($protID)=$modProtID=~/^(\d+)/;
			print LOST "\n$modProtID\t$proteinAliasLFQ{$protID}";
		}
		close LOST;
	}
	
	###<Lost proteins (All)
	my %proteinAlias;
	&getProteinAlias($dbh,$refLostProteins,\%proteinAlias);



	
	####<Remove temp data>####
	foreach my $job (1..$numAbundJobs) {
		rmtree "$resultDir/ABUND$job";
	}
	dirmove($absDir, "$resultDir/absolute") if ($computeAbsolute && -e $absDir);
	unlink "$dataDir/ghosts.txt" if -e "$dataDir/ghosts.txt";
	
	####<End quantification process>####
	my $projectID=&promsMod::getProjectID($dbh,$quantifID,'QUANTIFICATION');
	&endQuantification($dbh,$projectID,$refLostProteins,\%proteinAlias);
	
	exit;
}



sub checkInfiniteRatios { # globals: $labeling,$ratioType, @ratioTarget, %infiniteRatioProteins, %noSuperRatioProteins, $isModifQuantif
	# RULES:
	# 1. Compute normal ratio if num distinct peptides (seq+vmod) shared >= min. count (1 for PTM quantif, 3 otherwise).
	# 2. Compute normal ratio if num distinct peptides (seq+vmod) shared >= num. distinct in Test or in Ref
	# 3. If not:
	#	-Compute sum of xics ponderated by num. distinct peptides for the 3 options: Unique Test / Shared (sum of 2 states) / Unique Ref
	# 	-Assign ratio type to option with highest xic: +Inf / ratio / -Inf
	# OLD RULES below are obsolete (PP 09/04/21)
	## +Case 1:
	##		-FREE: Compute ratio if at least 1 peptide is shared
	## +Case 2 (Other labelings + FREE if case 1 failed)
	##		-Compute ratio if at least 3 peptides are shared
	## +Case 3: (cases 1 & 2 failed)
	##		-Compute total XIC for unique to test cond, unique to ref cond and shared by ref & test conds (sum of ref+test)
	##		-The highest XIC assigns the ratio type: XICref:-inf, XICtest:+inf, XICshared:normal ratio
	my ($protID,$refNumPepInCond)=@_; # can be a modProtID
	my $minNumDistinctShared=($isModifQuantif)? 1 : 3; # min num peptides (seq:vmod) to compute normal ratio without considering state-specific ones
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
			my ($xicShared,$xicMinusInf,$xicPlusInf)=(0,0,0,0); #$numShared,
			my (%usedPepCode,%sharedDistinct,%plusInfDistinct,%minusInfDistinct);
			foreach my $pepCode (keys %{$refNumPepInCond->{$context}{$testCondID}}) {
				my ($pepSeq,$vMod,$charge,$infRatioDataSrc)=split(':',$pepCode); # $infRatioDataSrc='-' for label-free
				
				foreach my $feature (keys %{$refNumPepInCond->{$context}{$testCondID}{$pepCode}}) {
					if ($refNumPepInCond->{$context}{$refCondID} && $refNumPepInCond->{$context}{$refCondID}{$pepCode} && $refNumPepInCond->{$context}{$refCondID}{$pepCode}{$feature}) { # shared
						#$xicShared+=($refNumPepInCond->{$context}{$testCondID}{$pepCode} >= $refNumPepInCond->{$context}{$refCondID}{$pepCode})? $refNumPepInCond->{$context}{$testCondID}{$pepCode} : $refNumPepInCond->{$context}{$refCondID}{$pepCode};

						$xicShared+=($refNumPepInCond->{$context}{$testCondID}{$pepCode}{$feature}[0] + $refNumPepInCond->{$context}{$refCondID}{$pepCode}{$feature}[0]); # Sum of both channels!! (PP 27/10/15)
	##$numShared++;
	#$numShared=1;
						# $sharedDistinct{"$pepSeq:$vMod:$feature"}=1;
						$sharedDistinct{"$pepSeq:$vMod"}=1; # feature is no longer considered
						$usedPepCode{"$pepCode:$feature"}=1;
					}
					else { # only in testCond
						$xicPlusInf+=$refNumPepInCond->{$context}{$testCondID}{$pepCode}{$feature}[0];
						#$plusInfDistinct{"$pepSeq:$vMod:$feature"}=1;
						$plusInfDistinct{"$pepSeq:$vMod"}=1; # feature is no longer considered
					}
				}
			}
			foreach my $pepCode (keys %{$refNumPepInCond->{$context}{$refCondID}}) {
				my ($pepSeq,$vMod,$charge,$infRatioDataSrc)=split(':',$pepCode); # $infRatioDataSrc='-' for label-free
				foreach my $feature (keys %{$refNumPepInCond->{$context}{$refCondID}{$pepCode}}) {
					next if $usedPepCode{"$pepCode:$feature"};
					$xicMinusInf+=$refNumPepInCond->{$context}{$refCondID}{$pepCode}{$feature}[0];
					#$minusInfDistinct{"$pepSeq:$vMod:$feature"}=1;
					$minusInfDistinct{"$pepSeq:$vMod"}=1; # feature is no longer considered
				}
			}
			##<Assigning ratio type
			# if ($algoType=~/^(TDA|DIA|MStats)$/) { # MS2 counts
			# 	##next if $numShared >= 6; # compute state mean with at least 6 frag
			# 	next if $numShared;
			# 	next if scalar keys %sharedDistinct >= 6; # normal ratio <===================================================== OVERWRITE THRESHOLD
			# }
			# else { # DDA MS1
			# 	next if ($algoType ne 'PEP_RATIO' && $numShared); # compute state mean with any available peptide former $labeling eq 'FREE' 27/04/20
			# 	next if scalar keys %sharedDistinct >= 3; # normal ratio <===================================================== OVERWRITE THRESHOLD
			# }
			my $numDistinctShared=scalar keys %sharedDistinct;
			next if $numDistinctShared >= $minNumDistinctShared; # Rule 1
			
			my $numDistinctPlusInf=scalar keys %plusInfDistinct;
			my $numDistinctMinusInf=scalar keys %minusInfDistinct;
			next if ($numDistinctShared >= $numDistinctPlusInf && $numDistinctShared >= $numDistinctMinusInf); # Rule 2
			
			#<Ponderate xic by num. distincts peptides
			$xicShared *= $numDistinctShared;
			$xicPlusInf *= $numDistinctPlusInf;
			$xicMinusInf *= $numDistinctMinusInf;

			next if ($xicShared >= $xicMinusInf && $xicShared >= $xicPlusInf); # Rule 3
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

sub checkForErrors {
	my ($mainDbh)=@_; # can be undef
	if (-e $errorFile && -s $errorFile) {
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


sub formatQtyFile {
	my ($quantifID, $userFile, $aLFQFile, $refProtIDs, $refProtMW) = @_;  # Uses also global $dbh

	my %mapAnaChan2StateRep;  # To match pairs analysis/channel to the state/bioRep/techRep combination
	my %obs2StateRep;  # We need to match observations to the state/rep combinations
	my %obs2AnaChannel;  # And the observations to the analysis/channel pairs
	my $hasFractions = 0;
	my $isMassCalib = 0;  # Calibration model with mass given as anchor protein quantity
	my $qtyMultiplier = 1;  # To accomodate the qty units order of magnitude
	my $sthQuantifAnnot = $dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION = ?");
	$sthQuantifAnnot->execute($quantifID);
	my ($quantifAnnot) = $sthQuantifAnnot->fetchrow_array;

	foreach my $annotStrg (split('::', $quantifAnnot)) {
		if ($annotStrg =~ /^STATES=(.*)/) {
		# 'state1;state2;---;stateN'
		# stateX = 'nbBioRep,bioRep1.bioRep2.---.bioRepN,#condID'
		# bioRepX = 'techRep1=techRep2=---=techRepN'
		# techRepX = 'frac1+frac2+---+fracN'
			my $statePos = 0;
			foreach my $state (split(/;/, $1)) {
				$statePos++;
				$hasFractions = 1 if ($state =~ /\+/);
				$state =~ s/#//g;  # remove all id tags
				my ($numBioRep, $quantiObsIDs, $condID) = split(/,/, $state);
				my $bioRepCount = 0;
				foreach my $bioReplicate (split(/\./, $quantiObsIDs)) {
					$bioRepCount++;
					my $techRepCount = 0;
					foreach my $techReplicate (split(/&/, $bioReplicate)) {
						$techRepCount++;
						my $fracCount = 0;
						foreach my $fraction (split(/\+/, $techReplicate)) {
							$fracCount++;
							my $channel;
							my ($obsID, $parQuantifID, $anaID, $targetPos) = split(/:/, $fraction);

							$sthQuantifAnnot->execute($parQuantifID);
							my ($parQuantifAnnot) = $sthQuantifAnnot->fetchrow_array;
							my ($xicLabeling) = ($parQuantifAnnot=~/LABEL=([^:]+)/);
							$xicLabeling = uc($xicLabeling);

							if ($xicLabeling =~ /SILAC|TMT|ITRAQ/) {
								foreach my $feature (split('::', $parQuantifAnnot)) {
									if ($feature !~ /=/) {
										# Matching between targetPos and channel
										my ($annotTargetPos, $channelName, $modif) = split(';', $feature);
										if ($annotTargetPos == $targetPos) {
											$channel = $channelName;
											last;
										} else {
											next;
										}
									}
								}
							}
							$obs2AnaChannel{$obsID}{'Analysis'} = $anaID;
							$obs2AnaChannel{$obsID}{'Channel'} = ($channel) ? $channel : 'NA';
							$obs2AnaChannel{$obsID}{'ChannelPos'} = ($targetPos) ? $targetPos : 0;

							$obs2StateRep{$obsID}{'StateID'}  = $condID;
							$obs2StateRep{$obsID}{'StatePos'} = $statePos;
							$obs2StateRep{$obsID}{'bioRep'}   = $bioRepCount;
							$obs2StateRep{$obsID}{'techRep'}  = $techRepCount;
							$obs2StateRep{$obsID}{'fraction'} = $fracCount;
						} # end of fraction
					} # end of tech replicate (sample)
				} # end of bio-replicate
			} # end of cond/state
		} elsif ($annotStrg =~ /^ABUND_MEASURES=(.*)/) {
			foreach my $meas (split(/;/, $1)) {
				if ($meas =~ /^aLFQ_(iBAQ|TOP)/) {
					my $pqi = $1;
					my  ($absMeas, $absModel, $pqiNorm, $outType, $topX, $topStrictness, $qtyUnit);
					if ($pqi eq 'TOP') {
						($absMeas, $absModel, $pqiNorm, $outType, $topX, $topStrictness, $qtyUnit) = split(',', $meas);
					} else {
						($absMeas, $absModel, $pqiNorm, $outType, $qtyUnit) = split(',', $meas);
					}
					if ($qtyUnit) {
						$isMassCalib = 1 if ($absModel eq 'calibration' && $qtyUnit =~ /mass/);
						$qtyMultiplier = ($qtyUnit =~ /^mmol|^conc_mmol_l|^mass_mg|^mass_conc_mg_l/)? 1e-03 :
										 ($qtyUnit =~ /^mass_ug|^mass_conc_ug_l/)? 1e-06 :
										 ($qtyUnit =~ /^copy_number/) ? 1 / 6.02e+23 :
										 1;  # default is 1 for /^mol|^conc_mol_l|^mass_g|^mass_conc_g_l/
					}
					last;
				}
			}
		}
	}
	$sthQuantifAnnot->finish;

	# Get analysis name from DB and map to state/rep/techRep
	my $sthGetAnaName = $dbh->prepare("SELECT NAME, DATA_FILE, WIFF_FILE FROM ANALYSIS WHERE ID_ANALYSIS = ?");
	foreach my $obsID (keys %obs2AnaChannel) {
		$sthGetAnaName->execute($obs2AnaChannel{$obsID}{'Analysis'});
		my ($anaName, $anaDataFile, $wiffFile) = $sthGetAnaName->fetchrow_array;
		$anaDataFile =~ s/\..*$//;
		$wiffFile =~ s/\..*$//;
		my $channelName = $obs2AnaChannel{$obsID}{'Channel'};
		my $channelPos = $obs2AnaChannel{$obsID}{'ChannelPos'};

		my $stateRepStrg = "State".$obs2StateRep{$obsID}{'StatePos'};
		$stateRepStrg   .= "_rep".$obs2StateRep{$obsID}{'bioRep'};
		$stateRepStrg   .= "_techRep".$obs2StateRep{$obsID}{'techRep'};

		# Multiple possibilities of analysis names and channel name/pos
		$mapAnaChan2StateRep{$anaName}{$channelName} = $stateRepStrg;
		$mapAnaChan2StateRep{$anaDataFile}{$channelName} = $stateRepStrg;
		$mapAnaChan2StateRep{$wiffFile}{$channelName} = $stateRepStrg;
		$mapAnaChan2StateRep{$anaName}{$channelPos} = $stateRepStrg;
		$mapAnaChan2StateRep{$anaDataFile}{$channelPos} = $stateRepStrg;
		$mapAnaChan2StateRep{$wiffFile}{$channelPos} = $stateRepStrg;
	}
	$sthGetAnaName->finish;

	# Now do the actual matching on the user-provided file
	my %colName2Idx;
	my $hasProteins = 0;
	my $hasChannels = 0;
	my %outLines;  # What will be written to build aLFQ input
	my %protNotFound;  # In case some proteins in the user-provided file are not found in the experiment

	open(USER_FILE, $userFile);
	while (<USER_FILE>) {
		chomp($_);
		if ($. == 1) {  # Check headers
			my $anaPattern 		= qr/[Aa]nalys[ie]s|[Ss]amples?/;
			my $channelPattern 	= qr/[Cc]hannels?/;
			my $protPattern 	= qr/[Pp]roteins?|[Gg]enes?/;
			my $qtyPattern 		= qr/[Qq]uantit(?:y|ies)|[Vv]alues?|[Aa]mouts?|[Mm]ass(?:es)?|[Cc]oncentrations?/;

			$hasProteins = 1 if ($_ =~ /$protPattern/);  # Quantity given by protein ? -> calibration model
			$hasChannels = 1 if ($_ =~ /$channelPattern/);  # SILAC, TMT or ITRAQ quantifs
			my @columns = split(',', $_);

			# Get index of each column
			$colName2Idx{'Analysis'} = (grep { $columns[$_] =~ /$anaPattern/ } 0..$#columns)[0];
			$colName2Idx{'Channel'}  = (grep { $columns[$_] =~ /$channelPattern/ } 0..$#columns)[0] if ($hasChannels);
			$colName2Idx{'Protein'}  = (grep { $columns[$_] =~ /$protPattern/ } 0..$#columns)[0] if ($hasProteins);
			$colName2Idx{'Quantity'} = (grep { $columns[$_] =~ /$qtyPattern/ } 0..$#columns)[0];
		} else {
			my @line = split(',', $_);
			my $analysis = $line[$colName2Idx{'Analysis'}];
			my $quantity = $line[$colName2Idx{'Quantity'}] * $qtyMultiplier;  # Convert quantity to mol, mol/L, g or g/L
			my $channel  = ($hasChannels) ? $line[$colName2Idx{'Channel'}] : 'NA';
			my $protein  = ($hasProteins) ? $line[$colName2Idx{'Protein'}] : '';

			if ($hasProteins) {  # Calibration model
				if (!$mapAnaChan2StateRep{$analysis}{$channel}) {  # Check if we can find the analysis and channel in the design
					die "Impossible to match the names provided in the user file to the experimental design. Be sure to follow the guidelines for the file containing the protein quantities.";
				} elsif (!$refProtIDs->{$protein}) {
					$protNotFound{$protein} = 1;  # TODO : report calibration proteins not quantified
				} else {
					my $stateRepStrg = $mapAnaChan2StateRep{$analysis}{$channel};
					my $proteinID = $refProtIDs->{$protein};
					if ($refProtMW->{$proteinID}) {
						if ($isMassCalib) {  # Convert mass to mol for calibration proteins
							$quantity = $quantity / $refProtMW->{$proteinID};  # MW in Da = g/mol in database
						}
						if ($hasFractions && $outLines{$stateRepStrg}{$proteinID}) {
							$outLines{$stateRepStrg}{$proteinID} += $quantity;  # In case quantity was given per fraction
						} else {
							$outLines{$stateRepStrg}{$proteinID} = $quantity;
						}
					} else {
						$protNotFound{$protein} = 1;  # Calibration protein does not have MW in DB
					}
				}
			} else {
				if (!$mapAnaChan2StateRep{$analysis}{$channel}) {  # Check if we can find the analysis and channel in the design
					die "Impossible to match the names provided in the user file to the experimental design. Be sure to follow the guidelines for the file containing the protein quantities.";
				} else {
					my $stateRepStrg = $mapAnaChan2StateRep{$analysis}{$channel};
					if ($outLines{$stateRepStrg}) {
						$outLines{$stateRepStrg} += $quantity;  # In case quantity was given per fraction
					} else {
						$outLines{$stateRepStrg} = $quantity;
					}
				}
			}
		}
	}
	close USER_FILE;

	open(aLFQ_FILE, ">$aLFQFile");
	if ($hasProteins) {  # Calibration model with anchor proteins
		print aLFQ_FILE "run_id,protein_id,concentration\n";
		foreach my $sample (keys %outLines) {
			foreach my $protID (keys %{$outLines{$sample}}) {
				print aLFQ_FILE "$sample,\"$protID\",$outLines{$sample}{$protID}\n";
			}
		}
	} else {  # Proportionality model
		print aLFQ_FILE "run_id,concentration\n";
		foreach my $sample (keys %outLines) {
			print aLFQ_FILE "$sample,$outLines{$sample}\n";
		}
	}
	close aLFQ_FILE;
}


sub runAbsoluteQuantif {  # Uses global $promsPath, $pathR, $fileStat, $numSteps, $absDir

	open(FILESTAT,">>$fileStat");
	print FILESTAT "3/$numSteps Computing protein quantification indices and absolute quantities\n";
	close FILESTAT;

	# Launching absolute quantification with aLFQ
	if (-d $absDir && -f "$absDir/aLFQ_info.txt") {
		# R package "aLFQ" needed
		system "cd $absDir; $pathR/R CMD BATCH --no-save --no-restore '--args $absDir $promsPath{R_scripts}' $promsPath{R_scripts}/aLFQ.R 2> $absDir/errors_aLFQ.txt";  # R CMD BATCH does not redirect errors, it stays in .Rout file. Should probably swith to Rscript call
	} else {
		die "$absDir/aLFQ_info.txt was not written, impossible to launch absolute quantification without the required information !";
	}
	exit;
}


sub parseAbsoluteQuantif {  # Globals : $absDir, %absMeasures

	my ($refProteinAbundance) = @_;

	my %absFiles=(	INFO		=> "$absDir/aLFQ_info.txt",
					PQI_VALUES	=> "$absDir/aLFQ_pqi.txt",
					ABS_VALUES	=> "$absDir/aLFQ_abs.txt"
	);

	my %measPatterns = (
		'MOL_PERCENT'		=> qr/\bMol_percent\b/,
		'CI_INF_MOL_PCT'	=> qr/\bMol_percent_CI_inf\b/,
		'CI_SUP_MOL_PCT'	=> qr/\bMol_percent_CI_sup\b/,
		'MOL'				=> qr/\bMol_amount\b/,
		'MOL_CONC'			=> qr/\bMol_concentration\b/,
		'CI_INF_MOL'		=> qr/\bMol_CI_inf\b/,
		'CI_SUP_MOL'		=> qr/\bMol_CI_sup\b/,
		'MASS_PERCENT'		=> qr/\bMass_percent\b/,
		'CI_INF_MASS_PCT'	=> qr/\bMass_percent_CI_inf\b/,
		'CI_SUP_MASS_PCT'	=> qr/\bMass_percent_CI_sup\b/,
		'MASS'				=> qr/\bMass_amount\b/,
		'MASS_CONC'			=> qr/\bMass_concentration\b/,
		'CI_INF_MASS'		=> qr/\bMass_CI_inf\b/,
		'CI_SUP_MASS'		=> qr/\bMass_CI_sup\b/,
		'GEO_SD'			=> qr/\bM_geo_sd\b/,
		'GEO_CV'			=> qr/\bM_geo_cv\b/
	);

	# Parse results of absolute values
	my %colName2Idx;
	open(ABS, "$absFiles{ABS_VALUES}") || die $!;
	while(<ABS>) {
		my $line = $_;
		chomp $line;
		my @lineData = split(/\t/, $line);
		if ($. == 1) { # Header : e.g. Condition ProteinID Mol_percent Quantity Mass_percent Mass
			# Check which measure has been computed and get index of each column for these measures
			foreach my $absMeas (keys %measPatterns) {
				if ($line =~ $measPatterns{$absMeas}) {
					$colName2Idx{$absMeas} = (grep { $lineData[$_] =~ /$measPatterns{$absMeas}/ } 0..$#lineData)[0];
					$absMeasures{$absMeas} = 1;  # Fill hash with abs measures recorded
				}
			}
		} else {
			my $state      = $lineData[0];
			my ($statePos) = $state =~ /State(\d+)/;
			my $protID0    = $lineData[1];
			foreach my $absMeas (keys %colName2Idx) {
				unless ($lineData[$colName2Idx{$absMeas}] eq 'NA') {
					$refProteinAbundance->{$protID0}{$statePos}{$absMeas} = $lineData[$colName2Idx{$absMeas}];
				}
			}
		}
	}
	close ABS;
}


sub getAllProtPepNb {  # Globals : $resultDir, %quantifParameters, %dataFiles, $countIdentifiedPep, %isGhostPeptide ?
	# In case only absolute quantif (aLFQ) was computed, otherwise, already done in jobProteinAbundance
	# Parse peptide results and
	# Compute peptide counts (NUM_PEP_TOTAL, NUM_PEP_USED & DIST_PEP_USED)

	# TODO : remove repetition of this part of code in &jobProteinAbundance
	# Something like &writeDataFiles($folder) should be ok

	my ($outDir) = @_;
	my (%allProteins, %proteinData);

	####<Loading list of ghost peptides>####
	my %isGhostPeptide;
	if ($countIdentifiedPep && -e "$dataDir/ghosts.txt") {
		open(GHOST,"$dataDir/ghosts.txt");
		while(<GHOST>) {
			chomp;
			$isGhostPeptide{$_}=1;
		}
		close GHOST;
	}

	open(PEP, "$resultDir/resultsPep.txt");
	my %pepHeader2Idx;
	while(<PEP>) {
		chomp;
		my @lineData = split(/\t/, $_);

		if ($.==1) { # 1st line of the file
			# Header
			my $colIdx=0;
			foreach my $colName (@lineData) {
				$pepHeader2Idx{$colName} = $colIdx;
				$colIdx++;
			}
		}
		else {
			# Data
			my $protID0 = $lineData[$pepHeader2Idx{'ProteinID'}];
			my $numPeptides = 0;
			my $numTruePeptides = 0;
			foreach my $pepID (split(/[\.\+\=]/, $lineData[$pepHeader2Idx{'PeptideId'}])) {
				$allProteins{$protID0}{$pepID} = 1; # for NUM_PEP_TOTAL & in case Labeling
				$numPeptides++;
				$numTruePeptides++ if ($countIdentifiedPep && !$isGhostPeptide{$pepID});
			}
			next if ($lineData[$pepHeader2Idx{'out'}] ne 'NA' && $lineData[$pepHeader2Idx{'out'}] ne 'outBalanceStates'); # skip outliers but not those considered outliers for ratios
			my $state 		= $lineData[$pepHeader2Idx{'Condition'}];
			my ($statePos) 	= $state =~ /State(\d+)/;
			my $bioRep 		= $lineData[$pepHeader2Idx{'replicate'}];
			my $techRep 	= $lineData[$pepHeader2Idx{'repTech'}];
			my $seqVarMod 	= (split(/_/, $lineData[$pepHeader2Idx{'Peptide'}]))[0];
			$proteinData{'NUM_PEP_USED'}{$protID0}{$statePos}  += $numPeptides;
			$proteinData{'NUM_TRUE_USED'}{$protID0}{$statePos} += $numTruePeptides if $countIdentifiedPep;
			$proteinData{'DIST_PEP_USED'}{$protID0}{$statePos}{$seqVarMod} = 1;
		}
	}
	close PEP;

	####<Saving all proteins & number of peptides>####
	open(PROT, ">$outDir/allProteins.txt") || die $!; # contains list of all proteins (before outlier filter) to compute lost proteins later
	print PROT "ProteinID\tnumPep";
	foreach my $protID0 (keys %allProteins) {
		print PROT "\n",$protID0,"\t",scalar keys %{$allProteins{$protID0}};
	}
	close PROT;

	my $headerStrg="ProteinID\tCondition\tValue\n";
	foreach my $measType (keys %dataFiles) {
		open (OUT,">$outDir/$dataFiles{$measType}") || die $!;
		print OUT $headerStrg;
		foreach my $protID0 (keys %{$proteinData{$measType}}) {
			foreach my $statePos (1..$numStates) {
				if ($proteinData{$measType}{$protID0}{$statePos}) {
					my $value = ($measType eq 'DIST_PEP_USED')? scalar keys %{$proteinData{$measType}{$protID0}{$statePos}} : $proteinData{$measType}{$protID0}{$statePos};
					print OUT "$protID0\tState$statePos\t$value\n";
				}
			}
		}
		close OUT;
	}
}


sub parsePrimaryAbundanceResults { # Globals: %dataFiles, %optDataFiles, %quantifParameters, $resultDir
	my ($numAbundJobs,$refProteinAbundance,$refAllPeptides,$refLostProteins,$refPeptideCount)=@_; # $refPeptideCount (true %proteinAbundance) is defined only in protein normalization context. Collects only peptide counts 

	if ($refPeptideCount) { # Protein normalization context
		open(FILESTAT,">>$fileStat");
		print FILESTAT "4/$numSteps Running protein-level normalization (1/2 Preparing data)\n";
		close FILESTAT;
	}

	foreach my $measType (@{$quantifParameters{'DB'}{'ABUND_MEASURES'}}) {
		if ($measType =~ /^aLFQ_/) {  # aLFQ_iBAQ or aLFQ_TOP
			# foreach my $absMeas (split(';', $measType)) {  # aLFQ_iBAQ,abs_model,normalization;MOL_PERCENT;MOL_CONC;...
			# 	if ($absMeas =~  /^(aLFQ_[^,]+),/) {
			# 		$absMeasures{$1} = 1;
			# 	} else {
			# 		$absMeasures{$absMeas} = 1;
			# 	}
			# }
		} else {
			$measType =~ s/,.+//; # MY_LFQ,x,x
			$dataFiles{$measType} = $optDataFiles{$measType};
		}
	}

	JOB:foreach my $job (1..$numAbundJobs) {
		my $jobDir="$resultDir/ABUND$job";

		###<Combining data after outlier filtering>###
		foreach my $measType (keys %dataFiles) {
			unless (-e "$jobDir/$dataFiles{$measType}") {
				if ($measType eq 'MY_LFQ') { # Do not fail if no LFQ data
					open(BAD_LFQ,">>$resultDir/errorLFQ.txt");
					print BAD_LFQ "File $dataFiles{$measType} not found for job #$job\n";
					close BAD_LFQ;
					copy("$jobDir/resultsPep.txt","$resultDir/resultsPep_$job.txt");
					next;
				}
				warn "File $dataFiles{$measType} not found for job #$job";
				last JOB;
			}
			my $usedRefAbundance=($refPeptideCount && $measType=~/_PEP_|_TRUE_/)? $refPeptideCount : $refProteinAbundance; # if defined, $refPeptideCount (true %proteinAbundance) collects only peptide counts
			&parseMeasureFile("$jobDir/$dataFiles{$measType}",$measType,$usedRefAbundance);
		}

		###<Combining data before outlier filtering & lost proteins>###
		&getAllProtInfo("$jobDir/allProteins.txt", $refProteinAbundance, $refAllPeptides, $refLostProteins);
	}
	
}

sub parseSecondaryAbundanceResults {
	my ($refProteinAbundance)=@_;

	####<Parsing R resultsPep.txt for each measure>####
	open(FILESTAT,">>$fileStat");
	print FILESTAT "6/$numSteps Processing results (1/2 Parsing files)\n";
	close FILESTAT;

	foreach my $measType (keys %dataFiles) {
		next if $measType=~/_PEP_|_TRUE_/;
		open(NORM,"$resultDir/resultsNorm/$measType/results/resultsPep.txt");
		my %pepHeader2Idx;
		my %normalizedValues;
		while(<NORM>) {
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
			my $state=$lineData[$pepHeader2Idx{'Condition'}];
			my ($statePos)=$state=~/(\d+)$/;
			my $bioRep=$lineData[$pepHeader2Idx{'replicate'}];
			my $techRep=$lineData[$pepHeader2Idx{'repTech'}];
			#<Normalized values
			$normalizedValues{$protID0}{$statePos}{$bioRep}{$techRep}=2**$lineData[$pepHeader2Idx{'log2Measure'}];
		}
		close NORM;

		##<Combined replicates (mean of replicates)
		foreach my $protID0 (keys %normalizedValues) {
			foreach my $statePos (keys %{$normalizedValues{$protID0}}) {
				my %bioRepAbundance;
				foreach my $bioRep (keys %{$normalizedValues{$protID0}{$statePos}}) {
					##<Take mean of tech reps
					my $numTechRep=0;
					foreach my $techRep (keys %{$normalizedValues{$protID0}{$statePos}{$bioRep}}) {
						$bioRepAbundance{$bioRep}+=$normalizedValues{$protID0}{$statePos}{$bioRep}{$techRep};
						$numTechRep++;
					}
					$bioRepAbundance{$bioRep}/=$numTechRep;
				}
				##<Take mean of bio reps
				my $numBioReps=scalar keys %bioRepAbundance;
				foreach my $bioRep (keys %bioRepAbundance) {
					$refProteinAbundance->{$protID0}{$statePos}{$measType}+=$bioRepAbundance{$bioRep};
				}
				$refProteinAbundance->{$protID0}{$statePos}{$measType}=int(0.5 + ($refProteinAbundance->{$protID0}{$statePos}{$measType}/$numBioReps));
			}
		}
	}
}

# sub parseMeasureFile { # Old (large) file format
# 	my ($dataFile, $measType, $refProteinAbundance) = @_;

# 	open(RES, "$dataFile") || die $!;
# 	my @statePos;
# 	while(<RES>) {
# 		chomp;
# 		my @lineData = split(/\t/, $_);
# 		if ($. == 1) {  # Header
# 			$statePos[0] = 0;  # not used
# 			foreach my $idx (1..$#lineData) {
# 				my ($statePos) = $lineData[$idx] =~ /State(\d+)/;
# 				$statePos[$idx]=$statePos;  # 1..n
# 			}
# 		}
# 		else {  # Data
# 			my $protID0 = $lineData[0];
# 			foreach my $idx (1..$#lineData) {
# 				next if ($lineData[$idx] eq 'NA' || abs($lineData[$idx]) < 1);
# 				$refProteinAbundance->{$protID0}{$statePos[$idx]}{$measType} = $lineData[$idx];
# 			}
# 		}
# 	}
# 	close RES;
# }
sub parseMeasureFile { # New (long) file format
	my ($dataFile,$measType,$refProteinAbundance) = @_;

	open(RES, "$dataFile") || die $!;
	my %header2Idx;
	my $hasReplicates=0;
	while(<RES>) {
		chomp;
		my @lineData = split(/\t/, $_);
		if ($.==1) {  # Header
			my $colIdx=0;
			foreach my $colName (@lineData) {
				$header2Idx{$colName}=$colIdx;
				$colIdx++;
			}
			$hasReplicates=1 if defined $header2Idx{'Replicate'};
			next;
		}
		# Data
		next if ($lineData[$header2Idx{'Value'}] eq 'NA' || abs($lineData[$header2Idx{'Value'}]) < 1);
		my ($statePos)=$lineData[$header2Idx{'Condition'}]=~/(\d+)/; # keep only pos
		if ($hasReplicates) {
			$refProteinAbundance->{$lineData[$header2Idx{'ProteinID'}]}{$statePos}{$lineData[$header2Idx{'Replicate'}]}{$measType} = $lineData[$header2Idx{'Value'}];
		}
		else {
			$refProteinAbundance->{$lineData[$header2Idx{'ProteinID'}]}{$statePos}{$measType} = $lineData[$header2Idx{'Value'}];
		}
	}
	close RES;
}

sub getAllProtInfo {
	my ($allProtFile, $refProteinAbundance, $refAllPeptides, $refLostProteins) = @_;
	open (PROT,"$allProtFile");
	while (<PROT>) {
		next if ($. == 1);
		chomp;
		my ($protID0, $numPep)=split(/\t/, $_);
		if ($refProteinAbundance->{$protID0}) {
			$refAllPeptides->{$protID0} = $numPep;
		}
		else {
			$refLostProteins->{$protID0} = 1;
		}
	}
	close PROT;
}


# TODO: Make clear choice for labeled quantif done with PEP_INTENSITY algo: treat as 100% Label-Free or mixed LF/Label ?????
# TODO: Move label-free peptide matching check further downstream for compatibility with PTM quantif
####>Revision history<####
# 2.19.7 [BUGFIX] Fix no retry limit for abundance sub job taking too long (PP 03/06/21)
# 2.19.6 [BUGFIX] Fix +/-INF peptide count for SILAC with multi-replicate reference state in non-Super Ratio context (PP 26/04/21)
# 2.19.5 [UPDATE] Modified +/-INF decision rules & quantif code version management moved here (PP 13/04/21) 
# 2.19.4 [BUGFIX] Fix wrong data structure for Abundance with tech replicates and protein-level normalization (PP 13/04/21)
# 2.19.3 [BUGFIX] Fix misspelled quantiSwath.R script name (PP 25/03/21)
# 2.19.2 [BUGFIX] Fix no data in PTM-enriched proteins when all PTMs are excluded (PP 15/03/21)
# 2.19.1 [BUGFIX] Fix import and pep nb files for absolute measures (VL 02/03/21)
# 2.19.0 [FEATURE] Context/motif-based PTM filtering for PTM-enriched proteins (PP 09/03/21)
# 2.18.3 [FEATURE] Record % missing values at protein/site-level for abundance (PP 12/02/21)
# 2.18.2 [BUGFIX] Add forgotten aggregation option flag for LFQ command (PP 02/02/21) 
# 2.18.1 [MERGE] master into ppoullet (PP 02/02/21)
# 2.18.0 [UPDATE] Multi-level normalization for Abundance & MG shared peptides rule & contaminant exlusion (PP 01/02/21)
# 2.17.2b [MINOR] Minor update in recovery of previously failed abundance quantification (PP 08/01/21)
# 2.17.2 [ENHANCEMENT] Add variability and uncertainty parameters for absolute quantification (VL 07/12/20)
# 2.17.1 [UPDATE] Improved run based on data from another quantification that has failed (PP 11/12/20)
# 2.17.0 [CHANGE] Compatible with Free residues quantification & &modifyFreeResidues added to promsQuantif.pm (PP 26/11/20)
# 2.16.12 [CHANGE] Ghost peptides detection now based on negative PEP_BEG (PP 26/11/20)
# 2.16.11 [BUGFIX] Broken bias correction with manually selected peptides works again (PP 17/11/20)
# 2.16.10 [BUGFIX] Ignore sequence context in reference mode (PP 14/11/20)
# 2.16.9b [ENHANCEMENT] Handle LFQ jobs on deficient cluster nodes & x3 more cluster memory for Site LFQ (PP 04/11/20)
# 2.16.9 [ENHANCEMENT] Add absolute quantification with aLFQ package (VL 13/10/20)
# 2.16.8b [CHANGE] Moved &recreateModifSites to promsQuantif.pm (PP 12/10/20)
# 2.16.8 [BUGFIX] Fix count of peptides for each prot before attribution of protein for shared peptides (VL 30/09/20)
# 2.16.7 [ENHANCEMENT] Use of systematic negative targetPos for state/replicate-peptide count of +/-INF ratio (PP 21/09/20)
# 2.16.6 [BUGFIX] Use of sortSmart to order DIA fragments in table.txt to avoid warnings on classical sort (PP 02/09/20)
# 2.16.5 [FEATURE] Compatible with PTMRS & Spectronaut site probability (PP 26/08/20)
# 2.16.4 [BUGFIX] Major fix in data processing due to dataset-wide attribution of "best" protein for shared peptides (PP 19/08/20)
# 2.16.3 [UPDATE] Change in "shared peptides" parameter handling & remove "shared proteins" (PP 17/08/20)
# 2.16.2 [UPDATE] Compatiblie with 'sharedPep' normalization filter (PP 14/08/20)
# 2.16.1 [BUGFIX] Fix bug in peptide count for MSstats (PP 12/08/20)
# 2.16.0 [FEATURE] Uses SQLite database to store quantification data & compatibility with site sequence context & python3 path (PP 30/07/20)
# 2.15.3 [BUGFIX] Restricts use of %pepID2QuantifSetSILAC to SILAC in PEP_RATIO mode (PP 30/07/20)
# 2.15.2 [FEATURE] Compatibility with Bias correction using shared proteins and LFQ stabilize large ratios options (PP 18/06/20)
# 2.15.1 [BUGFIX] Peptides with no PTM position score are no longer excluded if score threshold=0 (PP 05/06/20)
# 2.15.0 [FEATURE] Records number of non-ghost (true) peptides used (0 up to 5) in each state & abundance measures are optional (PP 20/05/20)
# 2.14.7 [BUGFIX] Force Experiment as 'A' for Labeled quantif if non-ratio algo is used (PP 06/04/20)
# 2.14.6 [ENHANCEMENT] Added % progression for data import in DB (PP 01/04/20)
# 2.14.5 [ENHANCEMENT] Checks for valid data.txt in Abundance context (PP 22/03/20)
# 2.14.4 [BUGFIX] Fix forgotten $dbh->disconnect in Abundance (PP 17/03/20)
# 2.14.3 [BUGFIX] in Abundance with technical replicates (PP 09/03/20)
# 2.14.2 [BUGFIX] Abundance now works properly with PTMs (PP 03/03/20)
# 2.14.1 [BUGFIX] Added forgotten newlines in progress status output (PP 16/02/20) 
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
