#!/usr/local/bin/perl -w

################################################################################
# launchExploratoryAnalyses.pl       1.0.5                                     #
# Authors: P. Poullet, S.Liva (Institut Curie)                                 #
# Contact: myproms@curie.fr                                                    #
# Launches PCA and Clustering analyses                                         #
# called by startExploratoryAnalyses.cgi                                       #
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
use File::Path qw(rmtree); # remove_tree
use File::Copy qw(move);
use promsConfig;
use promsMod;
my %promsPath=&promsConfig::getServerInfo;
my $dbh = &promsConfig::dbConnect;
my ($explorID, $anaType, $projectID, $metric, $method, $itemMetric) = @ARGV;

my $sthUpdateStatus = $dbh -> prepare("UPDATE EXPLORANALYSIS set STATUS=? where ID_EXPLORANALYSIS = ?");

my $explorDIR = "$promsPath{tmp}/exploratory_analysis/$explorID";
my $explorScript_Rout;
if($anaType eq "PCA"){
	system "cd $explorDIR; $promsPath{R}/R CMD BATCH --no-save --no-restore $promsPath{R_scripts}/PCA.R &";
	$explorScript_Rout="$explorDIR/PCA.Rout";
}
else {
	system "cd $promsPath{tmp}/exploratory_analysis/$explorID; $promsPath{R}/R CMD BATCH --no-save --no-restore '--args $metric $method $itemMetric' $promsPath{R_scripts}/cluster.R &";
	$explorScript_Rout="$explorDIR/cluster.Rout";
}

my $tempError = "$explorDIR/error.txt";

my $wait = 1;
my $error = 0;
my $count = 0;
my ($Rprocess, $RprocessFr);

while ($wait == 1) {
	sleep 2;
	if (!-e "$explorScript_Rout") {
		if ($count > 1200) {
			open(ERROR, ">$tempError");
			print ERROR "\n***No server Response or No R.out file***\n";
			close ERROR;
			$sthUpdateStatus -> execute(-2, $explorID) or die "Cannot execute: " . $sthUpdateStatus -> errstr();
			$sthUpdateStatus -> finish;
			$dbh -> commit;
			exit;
		}
		$count++;
	}
	elsif (! -z  $tempError) {
		$sthUpdateStatus -> execute(-2, $explorID) or die "Cannot execute: " . $sthUpdateStatus -> errstr();
		$sthUpdateStatus -> finish;
		$dbh -> commit;
		exit;
	}
	else {
		$Rprocess = `grep -c '> proc.time()' $explorScript_Rout`;
		chomp $Rprocess;
		if ($Rprocess) {
			$wait = 0;
		}
		else {
			$Rprocess = `grep -c '^Execution halted' $explorScript_Rout`;
			#$RprocessFr = `grep -c '^Exécution' $explorScript_Rout`;
			chomp $Rprocess;
			#chomp $RprocessFr;
			if($Rprocess) { # || $RprocessFr
				open(ERROR, ">$tempError");
				print  ERROR "\n**execution halted from R.**\n";
				close ERROR;
				$sthUpdateStatus -> execute(-2, $explorID) or die "Cannot execute: " . $sthUpdateStatus -> errstr();
				$sthUpdateStatus -> finish;
				$dbh -> commit;
				exit;
			}
		}
	}
}

if ($wait == 0) {
	$sthUpdateStatus -> execute(1, $explorID) or die "Cannot execute: " . $sthUpdateStatus -> errstr();
	$sthUpdateStatus -> finish;
	$dbh -> commit;
	move($explorDIR,"$promsPath{data}/exploratory_data/project_$projectID/$explorID");
}



#####>Revision history<####
# 1.0.5 Remove unsued variable (PP 18/07/15)
# 1.0.4 add argument $itemMetric for clustering  (cor or dist) (SL 16/12/14)
# 1.0.3 Minor changes (PP 09/07/14)
# 1.0.2 Move R excecution from startExploratoryAnalysis.cgi (GA 08/07/14)
# 1.0.1 Minor changes (PP 07/07/14)
# 1.0.0 first version (SL 10/03/14)
