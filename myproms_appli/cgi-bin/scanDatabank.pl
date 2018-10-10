#!/usr/local/bin/perl -w

################################################################################
# scanDatabank.pl               3.0.7                                          #
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
#-------------------------------------------------------------------------------# Performes background import of multiple analysis data into Database.
# Replace masterImport.pl & importAnalysis.pl
# Launched by storeAnalyses.cgi & by itself
# Tables: ANALYSIS and PROTEIN_VALIDATION
#############################################################################

$| = 1;
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;


########################
####>Immediate fork<####
########################
my $childPid = fork;
if ($childPid) { # parent here
	#sleep 2;
	exit; # returns control to storeAnalysis
}

##############################
####>Arguments/Parameters<####
##############################
my ($userID,$analysisStrg)=@ARGV;
$ENV{'REMOTE_USER'}=$userID;
die "ERROR: No Analyses specified!\n" unless $analysisStrg;
my @analysisList=split(',',$analysisStrg);


#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo('no_user'); # required for proper $dbh initialization
my %massAAave=&promsConfig::getMassAAave; # Average mass, needed for protein mass calculation
my %massATave=&promsConfig::getMassATave; # Average mass, needed for protein mass calculation
mkdir $promsPath{'logs'} unless -e $promsPath{'logs'};
my $year=strftime("%Y",localtime);

####>Disconnecting from parent's STDs<####
open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
open STDERR, ">>$promsPath{logs}/databankScan_$year.log";
open STDOUT, ">>$promsPath{logs}/databankScan_$year.log";
$| = 1;

#########################################
####>Waiting for previous job to end<####
#########################################
my $jobFlagFile="$promsPath{logs}/databankScan.job";
while (-e $jobFlagFile) {
	open(OLD_JOB,$jobFlagFile);
	my @jobTime=<OLD_JOB>;
	close OLD_JOB;
	unless ($jobTime[0]) {
		unlink($jobFlagFile);
		last;
	}
	chomp($jobTime[0]);
	if (time-$jobTime[0] > 900) { # more than 15 min ago => problem!
		unlink($jobFlagFile);
		last;
	}
	sleep 15;
}

###>New job flag file<###
open(JOB,">$jobFlagFile");
print JOB time;
close JOB;

my $jobStartTime=strftime("%Y%m%d%H%M%S",localtime);
print "\n>JOB START $jobStartTime user=$userID ana=$analysisStrg\n";


##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect('no_user');


####################################################
####>Fetching list of proteins for all analyses<####
####################################################
my (%protList,%protOrg,%protDes,%protMW,%protLength);
my $sthProt=$dbh->prepare("SELECT IDENTIFIER,ID_PROT_VALID,DB_RANK FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=? AND IDENTIFIER NOT LIKE 'DECOY_%'");
my $numProt=0;
foreach my $analysisID (@analysisList) {
	$sthProt->execute($analysisID);
	while (my ($identifier,$protValidID,$dbRank)=$sthProt->fetchrow_array) {
		$protList{$dbRank}{$identifier}{$analysisID}=$protValidID;
		$numProt++;
	}
}
$sthProt->finish;
print "+$numProt proteins to annotate\n";

###########################################################
####>Extracting protein annotations from databank file<####
###########################################################
#my ($databankID)=$dbh->selectrow_array("SELECT ID_DATABANK FROM ANALYSIS WHERE ID_ANALYSIS=$analysisList[0]"); # same DB set for all analyses
my ($numDatabanks,$protInfoErrors)=(0,0);
my $sthDB=$dbh->prepare("SELECT AD.ID_DATABANK,NAME,DB_RANK FROM ANALYSIS_DATABANK AD,DATABANK D WHERE ID_ANALYSIS=$analysisList[0] AND AD.ID_DATABANK=D.ID_DATABANK ORDER BY DB_RANK"); # same DB set for all analyses
$sthDB->execute;
while (my ($dbID,$dbName,$dbRank)=$sthDB->fetchrow_array) {
	$numDatabanks++;
	print "+ Scanning databank #$numDatabanks ($dbName)...";
	my ($error)=&promsMod::getProtInfo('silent',$dbh,$dbID,\@analysisList,\%protDes,\%protMW,\%protOrg,\%protLength,$protList{$dbRank});
	if ($error) {
		print $error;
		$protInfoErrors++;
	}
	else {print " Done.\n";}
}
$sthDB->finish;
if ($protInfoErrors == $numDatabanks) {
	$dbh->disconnect;
	unlink $jobFlagFile;
	die "\n**Too many errors. Databank scan failed**\n";
}

###########################################
####>Updating table PROTEIN_VALIDATION<####
###########################################
my $sthUpProt=$dbh->prepare("UPDATE PROTEIN_VALIDATION SET MW=?,PROT_LENGTH=?,PROT_DES=?,ORGANISM=? WHERE ID_PROT_VALID=?");
foreach my $dbRank (keys %protList) {
	foreach my $identifier (keys %{$protList{$dbRank}}) {
		#my $des=&promsMod::resize($protDes{$identifier},250); # max length allowed in table
		#my $organism=&promsMod::resize($protOrg{$identifier},100); # max length allowed in table
		foreach my $analysisID (keys %{$protList{$dbRank}{$identifier}}) {
			$sthUpProt->execute($protMW{$identifier},$protLength{$identifier},$protDes{$identifier},$protOrg{$identifier},$protList{$dbRank}{$identifier}{$analysisID}) || die $sthUpProt->errstr();
		}
	}
}
$sthUpProt->finish;

#################################
####>Updating table ANALYSIS<####
#################################
my $sthUpAna=$dbh->prepare("UPDATE ANALYSIS SET VALID_STATUS=0 WHERE ID_ANALYSIS=?");
foreach my $analysisID (@analysisList) {
	$sthUpAna->execute($analysisID) || die $sthUpAna->errstr();
}
$sthUpAna->finish;

$dbh->commit;
$dbh->disconnect;
unlink $jobFlagFile;

my $jobEndTime=strftime("%Y%m%d%H%M%S",localtime);
print "<JOB END $jobEndTime\n";


####>Revision history<####
# 3.0.7 Bug fix by better separation of multi-databanks (PP 25/03/14)
# 3.0.6 Uses mkdir instead of make_path (PP 10/03/14)
# 3.0.5 Minor typo correction (PP 30/10/13)
# 3.0.4 GPL license (PP 23/09/13)
# 3.0.3 Added _&lt;year&gt; tag to log file (PP 22/02/13)
# 3.0.2 Multi-databank search & databankScan log (PP 13/12/12)
# 3.0.1 Error handling after &getProtInfo (PP 13/09/12)
# 3.0.0 Merge of masterImport.pl & importAnalysis.pl (PP 09/08/12)
# 2.1.4 Comment on %promsPath declaration (PP 04/05/12)
# 2.1.3 script moved to cgi directory on myproms.curie.fr: no need for use lib path (PP 05/01/12)
