#!/usr/local/bin/perl -w

################################################################################
# launchProtRulerQuantif.pl        1.0.4                                       #
# Authors: V. Laigle (Institut Curie)                                          #
# Contact: myproms@curie.fr                                                    #
# called by startDesignProtRulerQuantif.cgi                                    #
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

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo('no_user');
my %cluster=&promsConfig::getClusterInfo;

###############################
####>Recovering parameters<####
###############################
my ($jobType,$userID,@extraArgs)=@ARGV;
my ($jobDir,$quantifType)=@extraArgs;
my $quantifDir="$promsPath{tmp}/quantification/$jobDir";
my $currentQuantifDir="$promsPath{tmp}/quantification/current";

my %quantifScripts=(
	'DESIGN:PROT_RULER'=>'runProtRulerQuantification.pl'  # Consistency with launchQuantifications.pl
);

my $quantifID;
my $extraArgStrg;
my $labelType='FREE';

my $dbh=&promsConfig::dbConnect('no_user');
#>Fetching parameters
my %quantifParams=&promsQuantif::getQuantificationParameters("$quantifDir/quantif_info.txt");
my $quantifName=$quantifParams{'DB'}{'QUANTIF_NAME'}[0];
my $qRootName=quotemeta($quantifName);
my $quantifAnnotStrg="LABEL=$labelType";
if ($quantifParams{'DB'}{'ALGO_TYPE'}) {
	my $quantifMethod=$quantifParams{'DB'}{'ALGO_TYPE'}[0];  # PROT_RULER
	if ($quantifMethod=~/PROT_RULER/) {
		$quantifAnnotStrg.="::SOFTWARE=myProMS;3.1::ALGO_TYPE=$quantifMethod"; # name;version::PROT_RULER
	} else {
		print STDERR "Error, the ALGO_TYPE was set to something else than PROT_RULER";
		$dbh->disconnect;
		exit;
	}
	if ($quantifParams{'DB'}{'PROTEINS'}) {
		$quantifAnnotStrg.="::PROTEINS=$quantifParams{'DB'}{'PROTEINS'}[0];$quantifParams{'DB'}{'PROTEINS'}[1]"; # Selected/Excluded proteins list
	}
	my ($quantifMethodID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE=UCASE('$quantifMethod')");
	my $focus='protein';  # Only protein quantif with Prot. Ruler (no peptide quantif)
	
	my $sthInsQ=$dbh->prepare("INSERT INTO QUANTIFICATION (ID_QUANTIFICATION_METHOD,NAME,FOCUS,STATUS,QUANTIF_ANNOT,UPDATE_DATE,UPDATE_USER) VALUES ($quantifMethodID,?,'$focus',-1,?,NOW(),'$userID')");
	$sthInsQ->execute($quantifName,$quantifAnnotStrg);
	$quantifID=$dbh->last_insert_id(undef, undef, 'QUANTIFICATION', 'ID_QUANTIFICATION');
	my $designID=$quantifParams{'DB'}{'ID_DESIGN'}[0];
	$dbh->do("UPDATE QUANTIFICATION SET ID_DESIGN=$designID WHERE ID_QUANTIFICATION=$quantifID");
	
	###> Update quantif_info.txt to make it work for watchQuantifications.cgi
	open (INFO,">>$quantifDir/quantif_info.txt");
	print INFO "QUANTIFICATIONS:\n$quantifID";
	close INFO;
	
	#>Generating job list
	unlink "$currentQuantifDir/$jobDir\_request.flag" if -e "$currentQuantifDir/$jobDir\_request.flag"; # in case multi-job launch
	open(WAIT,">$currentQuantifDir/$quantifID\_$jobDir\_wait.flag"); # flag file
	print WAIT '#';
	close WAIT;
	
	$sthInsQ->finish;

	$dbh->commit;
	
} else {
	print STDERR "Error, the ALGO_TYPE was not set to PROT_RULER";
	$dbh->disconnect;
	exit;		
}

##########################################
####>Launching quantification process<####
##########################################
my $projectID=&promsMod::getProjectID($dbh,$quantifID,'QUANTIFICATION');
$dbh->do("INSERT INTO JOB_HISTORY (ID_JOB, ID_USER, ID_PROJECT, TYPE, STATUS, FEATURES, SRC_PATH, LOG_PATH, ERROR_PATH, STARTED) VALUES('$jobDir', '$userID', $projectID, 'Quantification [DESIGN:PROT_RULER]', 'Queued', 'ID_QUANTIFICATION=$quantifID;TYPE=DESIGN:PROT_RULER;QUANTIFICATION_NAME=$quantifName', '$quantifDir', '$quantifDir/status\_$quantifID.out', '$currentQuantifDir/$quantifID\_$jobDir\_error.txt', NOW())");
$dbh->commit;

if ($cluster{'on'}) {
	
	###>Estimating ressources required<###
	my $dbh=&promsConfig::dbConnect('no_user');
	my $projectID=&promsMod::getProjectID($dbh,$quantifID,'QUANTIFICATION');
	
	# Counting number of proteins quantified for memory usage
	my @selectedQuantifs=split(';', $quantifParams{'DB'}{'QUANTIF_IDS'}[0]);
	my @selectedQuantifIDs;
	foreach my $quantif (@selectedQuantifs) {
		my ($selQuantifID, $tarPos)=split('_', $quantif);
		push @selectedQuantifIDs, $selQuantifID;
	}
	my $selectedQuantifStrg=join(',', @selectedQuantifIDs);
	my $nbProts=$dbh->selectrow_array("SELECT COUNT(1) FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION IN ($selectedQuantifStrg) AND ID_QUANTIF_PARAMETER=$quantifParams{'R'}{'int_metric'}[0]");

	my ($maxHours,$maxMem,$numCPU,$jobName);
	if ($quantifType eq 'DESIGN:PROT_RULER') {
		$maxHours=2;  # Work with proteins, usually less than 10k lines and less than 20 samples/columns
		$maxMem=int(1.3 + 1E-6 * $nbProts);  # Usually 1Gb, 2Gb if more than 200 000 proteins in total (proteins are counted multiple times if they are in different samples/states). Around 100Mb used with 50 000 proteins 
		$maxMem.='Gb';
		$numCPU=1;  # A priori, no need for more 
		$jobName="myProMS_protRulerQuant_$quantifID";
	}

	###>Running job on cluster<###
	my $cgiUnixDir=`pwd`;
	$cgiUnixDir=~s/\/*\s*$//;
	# cd is required for script to find myproms .pm files
	my $commandString="export LC_ALL=\"C\"; cd $cgiUnixDir; $cluster{path}{perl}/perl $quantifScripts{$quantifType} $quantifID $jobDir 2> $currentQuantifDir/$quantifID\_$jobDir\_error.txt";
	my %jobParameters=(
		maxMem=>$maxMem,
		numCPUs=>$numCPU,
		maxHours=>$maxHours,
		jobName=>$jobName,
		pbsRunDir=>$cgiUnixDir,
		commandBefore=>"mv $currentQuantifDir/$quantifID\_$jobDir\_wait.flag $currentQuantifDir/$quantifID\_$jobDir\_run.flag" # run flag file
	);
	my ($pbsError,$pbsErrorFile,$jobClusterID)=$cluster{'runJob'}->($quantifDir,$commandString,\%jobParameters);
	
	# Add cluster job id to current job in DB
    $dbh->do("UPDATE JOB_HISTORY SET ID_JOB_CLUSTER='C$jobClusterID' WHERE ID_JOB='$jobDir'");
    $dbh->commit;
	
	if ($pbsError) { # move PBS error message to job error file
		system "cat $pbsErrorFile >> $currentQuantifDir/$quantifID\_$jobDir\_error.txt";
		my $dbh=&promsConfig::dbConnect('no_user'); # reconnect
		$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID"); # Failed
		$dbh->commit;
	}
	elsif (-e "$promsPath{quantification}/project_$projectID/quanti_$quantifID") { # Move cluster files to quantif final directory (except if error in quantif)
		system "mv $quantifDir/*.txt $promsPath{quantification}/project_$projectID/quanti_$quantifID/.";
		system "mv $quantifDir/*.sh $promsPath{quantification}/project_$projectID/quanti_$quantifID/.";
	}

} # end of cluster launch
else { # Local launch (web server): No cluster
	# Add process PID to current job in DB
    $dbh->do("UPDATE JOB_HISTORY SET ID_JOB_CLUSTER='L$$' WHERE ID_JOB='$jobDir'");
    $dbh->commit;
	
	rename("$currentQuantifDir/$quantifID\_$jobDir\_wait.flag","$currentQuantifDir/$quantifID\_$jobDir\_run.flag"); # run flag file
	system "./$quantifScripts{$quantifType} $quantifID $jobDir $extraArgStrg 2> $currentQuantifDir/$quantifID\_$jobDir\_error.txt";
} # end of local launch

$dbh->disconnect;

unlink "$currentQuantifDir/$quantifID\_$jobDir\_run.flag";
sleep 15;

####>Check for error & clean root quantif directory<####
unless (-s "$currentQuantifDir/$quantifID\_$jobDir\_error.txt") {
	if (-e $quantifDir) {
		my $numFiles=`ls -l $quantifDir | wc -l`;
		chomp($numFiles);
		if ($numFiles*1 <= 2) {  # total + only quantif_info.txt
			rmtree($quantifDir);
			unlink "$currentQuantifDir/$quantifID\_$jobDir\_error.txt";
		}
	}
	else {unlink "$currentQuantifDir/$quantifID\_$jobDir\_error.txt";}
}

#####>Revision history<#####
# 1.0.4 [CHANGES] Add quantification name in job information (VS 19/11/19)
# 1.0.3 Compatibility with the new job monitoring system (VS 09/10/19)
# 1.0.2 Adapt code to new intensity metric + minor bug correction (VL 25/09/19)
# 1.0.1 Add a simple estimation of resources usage (VL 23/08/19)
# 1.0.0 Created (VL 26/07/2019)

