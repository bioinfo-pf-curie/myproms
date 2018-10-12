#!/usr/local/bin/perl -w
################################################################################
# launchMotifEnrichmentAnalyses.pl       1.0.1                                 #
# Authors: P. Poullet, S.Liva (Institut Curie)                                 #
# Contact: myproms@curie.fr                                                    #
# called by startMotifAnalyses.cgi											   #
#launch R script for maotif enrichment analysis                                #
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
use File::Copy::Recursive qw(dirmove);
use promsConfig;
use promsMod;
use POSIX qw(strftime); # to get the time
my %promsPath=&promsConfig::getServerInfo;
my $dbh = &promsConfig::dbConnect;
my ($motifID, $projectID,$centralRes, $occurence,$pvalCutoff,$bgValue,$width,$randomNbSeq) = @ARGV;
my $sthUpdateStatus = $dbh -> prepare("UPDATE EXPLORANALYSIS set STATUS=? where ID_EXPLORANALYSIS = ?");
my $motifDIR = "$promsPath{tmp}/motifX/$motifID";
my $logFile="$motifDIR/motif.log";
open (LOG,">$logFile");
print LOG ">PARAMETERS: @ARGV\n\n";
print LOG ">Running Motif Enrichment Analysis:\n";

#############################
####<Launching R scripts>####
#############################
print LOG ">COMMAND: $promsPath{R}/R CMD BATCH --no-save --no-restore '--args foreground_sequence.txt  $centralRes $occurence $pvalCutoff $bgValue $width $randomNbSeq' \n\n";
print LOG " -Starting (",strftime("%Y-%m-%d %H:%M:%S",localtime),")...\n"; #<BR>

##script bash pour PBS
#my $bashFile="$motifDIR/makeMotifAnalysis_${userID}.sh";
#open (BASH,">$bashFile");
##print BASH qq|
###!/bin/bash
##cd $motifDIR
##$promsPath{R}/R CMD BATCH --no-save --no-restore '--args foreground_sequence.txt  $centralRes $occurence $pvalCutoff $bgValue $width $randomNbSeq' $promsPath{R_scripts}/motifEnrichment.R $RoutFile &"
##|;
#
#print BASH qq
#|#!/bin/bash
###resources
##PBS -l walltime=02:00:00
##PBS -l mem=900mb
##PBS -l nodes=1:ppn=1
##PBS -q batch
###Information
##PBS -N makeMotifAnalysis_${userID}
##PBS -M stephane.liva\@curie.fr,patrick.poullet\@curie.fr
##PBS -m ae
##PBS -o $motifDIR/PBS.txt
##PBS -e $motifDIR/PBSerror.txt
#
###sur les serveurs de calcul
#cd $motifDIR
#$promsPath{R}/R CMD BATCH --no-save --no-restore '--args foreground_sequence.txt $centralRes $occurence $pvalCutoff $bgValue $width $randomNbSeq' $promsPath{R_scripts}/motifEnrichment.R $RoutFile
#|;

#close BASH;
##system "chmod 775 $bashFile";
#my $modBash=0775;
#chmod $modBash, $bashFile;

###sur les serveurs de calcul
#system "qsub $bashFile > $motifDIR/torqueID.txt";
my $Rscript="$promsPath{R_scripts}/motifx-functions.R";
system "cd $motifDIR;$promsPath{R}/R CMD BATCH --no-save --no-restore '--args foreground_sequence.txt $centralRes $occurence $pvalCutoff $bgValue $width $randomNbSeq $Rscript' $promsPath{R_scripts}/motifEnrichment.R&";
my $RoutFile="$motifDIR/motifEnrichment.Rout";
my $tempError = "$motifDIR/error.txt";
my $wait = 1;
my $error = 0;
my $count = 0;
my ($Rprocess, $RprocessFr);

while ($wait == 1) {
	sleep 2;
	if (!-e "$RoutFile") {
		if ($count > 1200) {
			open(ERROR, ">$tempError");
			print ERROR "\n***No server Response or No R.out file***\n";
			close ERROR;
			$sthUpdateStatus -> execute(-2, $motifID) or die "Cannot execute: " . $sthUpdateStatus -> errstr();
			$sthUpdateStatus -> finish;
			$dbh -> commit;
			print LOG " -End [ERROR] (",strftime("%Y-%m-%d %H:%M:%S",localtime),")...\n";
			close LOG;
			exit;
		}
		$count++;
	}
	elsif (! -z  $tempError) {
		$sthUpdateStatus -> execute(-2, $motifID) or die "Cannot execute: " . $sthUpdateStatus -> errstr();
		$sthUpdateStatus -> finish;
		$dbh -> commit;
		print LOG " -End [ERROR] (",strftime("%Y-%m-%d %H:%M:%S",localtime),")...\n";
		close LOG;
		exit;
	}
	else {
		$Rprocess = `grep -c '> proc.time()' $RoutFile`;
		chomp $Rprocess;
		if ($Rprocess) {
			$wait = 0;
		}
		else {
			$Rprocess = `grep -c '^Execution halted' $RoutFile`;
			chomp $Rprocess;
			if($Rprocess) { # || $RprocessFr
				open(ERROR, ">$tempError");
				print  ERROR "\n**execution halted from R.**\n";
				close ERROR;
				$sthUpdateStatus -> execute(-2, $motifID) or die "Cannot execute: " . $sthUpdateStatus -> errstr();
				$sthUpdateStatus -> finish;
				$dbh -> commit;
				print LOG " -End [ERROR] (",strftime("%Y-%m-%d %H:%M:%S",localtime),")...\n";
				close LOG;
				exit;
			}
		}
	}
}
print LOG " -End (",strftime("%Y-%m-%d %H:%M:%S",localtime),")...\n";
close(LOG);

if ($wait == 0) {
	$sthUpdateStatus -> execute(1, $motifID) or die "Cannot execute: " . $sthUpdateStatus -> errstr();
	$sthUpdateStatus -> finish;
	$dbh -> commit;
	dirmove($motifDIR,"$promsPath{data}/exploratory_data/project_$projectID/$motifID");
}

#####>Revision history<####
# 1.0.1 Uses dirmove instead of move (PP 12/10/18)
# 1.0.0 first version (SL 07/06/17)
