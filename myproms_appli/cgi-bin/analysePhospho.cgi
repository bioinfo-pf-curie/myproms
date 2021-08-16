#!/usr/local/bin/perl -w

#############################################################################
# analysePhospho.cgi         2.1.0                                          #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                   #
# Contact: myproms@curie.fr                                                 #
# Script processing PhosphoRS analyses started by selectAnalyses.cgi        #
#############################################################################
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

use strict;
use POSIX ":sys_wait_h"; # for WNOHANG
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use POSIX qw(strftime); # to get the time
use File::Path qw(rmtree);
use promsConfig;
use promsMod;
use phosphoRS;
#use IO::Uncompress::Unzip qw(unzip $UnzipError); # needed for &promsMod::extractSpectrumMSF
use XML::Simple; # needed for &promsMod::extractSpectrumMSF

####>Configuration<####
my %promsPath=&promsConfig::getServerInfo('no_user');
my $phosphoRSPath="$promsPath{tmp}/phosphoRS";
my $currentPRSPath="$phosphoRSPath/current";

                      ######################
##########################>>>CGI mode<<<##########################
                      ######################
if ($#ARGV < 0) {
	
	####################
	####>Delete PRS<####
	####################
	if (param('deleteID')) {
		my $branchID = param('branchID');
		print header(-'content-encoding' => 'no', -'charset' => 'UTF-8');
		warningsToBrowser(1);
		print qq
|<HEAD>
<TITLE>Deleting PhosphoRS Analysis</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER><FONT class="title">Deleting Previous PhosphoRS Analyses</FONT></CENTER>
<BR><BR>
|;
		my $dbh=&promsConfig::dbConnect;
		my $sthAN = $dbh->prepare("SELECT NAME FROM ANALYSIS WHERE ID_ANALYSIS=?");
		####>Looping through analyses<####
		foreach my $anaID (split(',',param('deleteID'))) {
			$anaID=~s/\D+//g; # clean parameter
			next unless $anaID;
			$sthAN->execute($anaID);
			my ($anaName)=$sthAN->fetchrow_array;
			print "<BR><FONT class=\"title3\">&nbsp;&nbsp;- Analysis $anaName...";
			&deletePRS($dbh,$anaID);
			print " Done.</FONT><BR>\n";
		}
		$sthAN->finish;
		$dbh->disconnect;
		sleep 3;
		print qq
|<SCRIPT type="text/javascript">
window.location="$promsPath{cgi}/selectAnalyses.cgi?ID=$branchID&callType=phosphoRS";
</SCRIPT>
</BODY>
|;
		exit;
	}
	######################
	####>Clear errors<####
	######################
	elsif (param('clearErrors')) { # Called from selectAnalyses.cgi
		my @analysisList=split(',',param('anaList'));
		my $branchID = param('branchID');
		
		my $anaKeyword=(scalar @analysisList > 1)? 'analyses' : 'analysis';
		print header(-'content-encoding' => 'no', -'charset' => 'UTF-8');
		warningsToBrowser(1);
		print qq
|<HEAD>
<TITLE>Clearing PhosphoRS Errors</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<BR><BR>
<FONT class="title2">Clearing failed PhosphoRS $anaKeyword:...|;

		my $dbh = &promsConfig::dbConnect;
		my $sthDelJob=$dbh->prepare("DELETE FROM JOB_HISTORY WHERE ID_JOB=?");
		my $sthPRS=$dbh->prepare("SELECT 1 FROM QUERY_VALIDATION WHERE ID_ANALYSIS=? AND INFO_PEP1 LIKE '%,PRS=%' LIMIT 1");
		foreach my $anaID (@analysisList) {
			my $jobDir;
			foreach my $currFile (glob("$currentPRSPath/$anaID\_*")) {
				unless ($jobDir) {
					my $fileName=(split(/\//,$currFile))[-1];
					($jobDir)=$fileName=~/_([^_]+)_/;
				}
				unlink $currFile;
			}
			rmtree("$phosphoRSPath/$jobDir");
			my $prsParamFile = "$promsPath{valid}/ana_$anaID/PRSparam_ana_$anaID.txt";
			if (-e $prsParamFile) {
				#>Check for partially recorded data
				$sthPRS->execute($anaID);
				my ($hasPRS)=$sthPRS->fetchrow_array;
				if ($hasPRS) {
					&deletePRS($dbh,$anaID);
				}
				else {
					unlink $prsParamFile;
				}
			}
			$sthDelJob->execute($jobDir);
			print '.';
		}
		$dbh->commit;
		$dbh->disconnect;
		print " Done.</FONT>\n";
		sleep 3;

		print qq
|<SCRIPT type="text/javascript">
window.location="$promsPath{cgi}/selectAnalyses.cgi?ID=$branchID&callType=phosphoRS";
</SCRIPT>
</BODY>
|;
		exit;
	}

	#------------------------------------#
	# Preparing Phosphorylation Analyses #
	#------------------------------------#

	####>parameters<####
	my $probThreshold = param('probThreshold');
	my $massTolerance = param('massTolerance');
	my $activationType = param('activationType');
	my $overWrite=param('overwrite') || 0;
	my @analysisList = param('anaList');

	my $dbh=&promsConfig::dbConnect;
	
    # my $projectID = $dbh->selectrow_array("SELECT DISTINCT E.ID_PROJECT FROM EXPERIMENT E INNER JOIN SAMPLE S ON S.ID_EXPERIMENT=E.ID_EXPERIMENT INNER JOIN ANALYSIS A ON A.ID_SAMPLE=S.ID_SAMPLE WHERE ID_ANALYSIS IN (".join(',', @analysisList).") LIMIT 1");
	my $projectID = &promsMod::getProjectID($dbh,$analysisList[0],'analysis');

	print header(-'content-encoding' => 'no', -'charset' => 'UTF-8');
	warningsToBrowser(1);
	print qq
|<HEAD>
<TITLE>Preparing Phosphorylation Analyses</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT type="text/javascript">
function newPhosphoRS() {
	var branchID=top.promsFrame.navFrame.getSelectedBranchID();
	window.location="./selectAnalyses.cgi?ID="+branchID+"&callType=phosphoRS";
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER><FONT class="title">Launching Phosphorylation Site Analyses</FONT></CENTER>
<BR>
<BR>
|;
	####>job directories<####
	mkdir $phosphoRSPath unless -e $phosphoRSPath;
	mkdir $currentPRSPath unless -e $currentPRSPath;
	
	my $masterJobCode=strftime("%Y%m%d%H%M%S",localtime);
	while (glob "$phosphoRSPath/$masterJobCode*") { # to prevent multi-user launch collision
		sleep 2;
		$masterJobCode=strftime("%Y%m%d%H%M%S",localtime);
	}
	my $numJobs=scalar @analysisList;
	print '<FONT class="title3">Launching PhosphoRS process';
	print "es" if $numJobs > 1;
	print " [Master job #$masterJobCode]:\n";
	
	####>Looping through analyses<####
	my $sthAN = $dbh->prepare("SELECT NAME FROM ANALYSIS WHERE ID_ANALYSIS=?");
	my $jobPos=0; 
	foreach my $anaID (@analysisList) {
		$jobPos++;
		sleep 1 if $jobPos > 1; # wait 1 sec between jobs
		
		$sthAN->execute($anaID);
		my ($anaName)=$sthAN->fetchrow_array;
		print "<BR>&nbsp;&nbsp;-$jobPos/$numJobs: Analysis '$anaName'";
		if (-e "$promsPath{valid}/ana_$anaID/PRSparam_ana_$anaID.txt") { # $paramFile
			print ' [';
			&deletePRS($dbh,$anaID,1);
			print ']';
		}
		my $jobDir=$masterJobCode.'.'.$jobPos;
		open(FLAG,">$currentPRSPath/$anaID\_$jobDir\_wait.flag"); # flag file used by watchQuantif
		print FLAG "#";
		close FLAG;
		my $jobPath="$phosphoRSPath/$jobDir";
		mkdir $jobPath || die "ERROR detected: $!";
		open (INFO,">$jobPath/prs_info.txt"); # item valueR valueDB
		print INFO qq
|USER=$ENV{REMOTE_USER}
ANAID=$anaID
probThreshold=$probThreshold
massTolerance=$massTolerance
activationType=$activationType
overWrite=$overWrite
|;
		close INFO;	
	}
	$sthAN->finish;
    $dbh->disconnect;
	print "<BR><BR> Done.</FONT><BR>\n";

	####>Forking to launch master job in background<####
	my $childPid = fork;
	unless ($childPid) { # child here
		#>Disconnecting STDOUT & STDIN from Apache (STDERR still goes to Apache error log!)
		open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
		open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
		&launchAllJobs($masterJobCode,\@analysisList,$probThreshold,$massTolerance,$activationType,$overWrite);
		exit;
	}
	
	print qq
|<CENTER>
<BR><FONT class="title2"><BR>PhosphoRS is running in background mode. You can continue using myProMS.</FONT>
<BR><BR><INPUT type="button" value="New PhosphoRS Analysis" onclick="newPhosphoRS();">
</CENTER>
|;
	####>Calling monitoring popup window<####
	sleep 3;
	print qq
|<SCRIPT type="text/javascript">
var monitorJobsWin=window.open("$promsPath{cgi}/monitorJobs.cgi?filterType=Phospho&filterDateNumber=1&filterDateType=DAY&filterStatus=Queued&filterStatus=Running&filterProject=$projectID",'monitorJobsWindow','width=1200,height=500,scrollbars=yes,resizable=yes');
monitorJobsWin.focus();
</SCRIPT>
</BODY>
</HTML>
|;
	exit;

} # <--- End of CGI mode
	
	
                      #####################################
##########################>>>BASH mode (Single job)<<<##########################
                      #####################################

####>Recovering parameters<####
my ($userID,$jobDir,$anaID,$probThreshold,$massTolerance,$activationType,$overWrite)=@ARGV;
my $jobPath="$phosphoRSPath/$jobDir";
my $fileStat="$jobPath/status.out";

rename("$currentPRSPath/$anaID\_$jobDir\_wait.flag","$currentPRSPath/$anaID\_$jobDir\_run.flag"); # run flag file

open(FILESTAT,">$fileStat") || die "Error while opening $fileStat";
print FILESTAT "Started ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
print FILESTAT "1/3 Preparing data\n";
close FILESTAT;

my ($infoPepStrg,@infoPepList);
@infoPepList = map { "INFO_PEP$_" } (1..10);
$infoPepStrg = join(', ',@infoPepList);

####>Database connection<####
my $dbh=&promsConfig::dbConnect('no_user');

my ($phosphoModID)=$dbh->selectrow_array("SELECT ID_MODIFICATION FROM MODIFICATION WHERE UNIMOD_ACC=21");

my ($specificity)=$dbh->selectrow_array("SELECT SPECIFICITY FROM ANALYSIS_MODIFICATION WHERE ID_MODIFICATION=$phosphoModID AND ID_ANALYSIS=$anaID");
$specificity=~s/,//g; # S,T,Y -> STY

my ($countUnchanged,$countFlagged,$countChanged,$countConfirmed,$countTotal) = (0,0,0,0,0);
my ($anaName,$dataFile,$fileFormat,$maxRank) = $dbh->selectrow_array("SELECT NAME,DATA_FILE,FILE_FORMAT,MAX_RANK FROM ANALYSIS WHERE ID_ANALYSIS=$anaID");
$dbh->disconnect;
my $anadir = "$promsPath{valid}/ana_$anaID";
my $fullDataFileName = "$anadir/$dataFile";

####>Preparing data for PhosphoRS<####
my $phosphoRS = new phosphoRS(AnaID => $anaID, fullJobDir =>$jobPath, File => $fullDataFileName, FileFormat => $fileFormat, MassTolerance => $massTolerance, ActivationType => $activationType);

my %noAmbiguityQueries;
if ($fullDataFileName =~ /([^\/]+)\.pdm$/) {
	
	open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
	print FILESTAT "1.1/3 Fetching analysis information from MSF file\n";
	close FILESTAT;
	my $cutFileName = $1;
	$cutFileName =~ s/_\d_?$//; # use s/_\d_?.*$// to be compatible with current splitted file format
	
	$dbh=&promsConfig::dbConnect('no_user');
	
	my $projectID = &promsMod::getProjectID($dbh, $anaID, 'ANALYSIS');
	
	my $dir = "$promsPath{valid}/multi_ana/proj_$projectID";
	my $msfFileName = "$dir/$cutFileName.msf";
	die "File: $dir/$cutFileName.msf not found!" unless -e $msfFileName;

	#my @infoPepList = map { "INFO_PEP$_ LIKE \"%Phospho%\"" } (1..10);
	#my $pepInfoQuery = join ' OR ', @infoPepList;
	#my $sthSpecId = $dbh->prepare("SELECT EXT_SPECTRUMID, QUERY_NUM FROM QUERY_VALIDATION WHERE VALID_STATUS>=? AND ID_ANALYSIS=? AND QUERY_NUM>0 AND ($pepInfoQuery)");
	my @spectrumList;
	my %sthPep;
	foreach my $i (1..$maxRank){
		$sthPep{$i} = $dbh->prepare("SELECT INFO_PEP$i FROM QUERY_VALIDATION WHERE ID_QUERY=?");
	}
	#my $sthSpecId = $dbh->prepare("SELECT EXT_SPECTRUMID,QUERY_NUM,V.ID_QUERY,GROUP_CONCAT(CONCAT(PEP_RANK,':',POS_STRING) SEPARATOR ',') FROM QUERY_VALIDATION V,QUERY_MODIFICATION M WHERE V.ID_QUERY=M.ID_QUERY AND VALID_STATUS>=-3 AND QUERY_NUM>0 AND ID_MODIFICATION=$phosphoModID AND ID_ANALYSIS=$anaID GROUP BY QUERY_NUM");
	my ($addQuery)=($overWrite)? '' : ' M.REF_POS_STRING IS NULL AND';###> Avoid to take into account modified queries
	my $sthSpecId = $dbh->prepare("SELECT EXT_SPECTRUMID,QUERY_NUM,V.ID_QUERY,GROUP_CONCAT(CONCAT(PEP_RANK,':',POS_STRING) SEPARATOR ',') FROM QUERY_VALIDATION V,QUERY_MODIFICATION M WHERE V.ID_QUERY=M.ID_QUERY AND$addQuery VALID_STATUS>=-3 AND QUERY_NUM>0 AND ID_MODIFICATION=$phosphoModID AND ID_ANALYSIS=$anaID GROUP BY V.ID_QUERY,QUERY_NUM,EXT_SPECTRUMID");
	$sthSpecId->execute;
	while(my ($extSpectrumID,$queryNum,$queryID,$phosPosStrg) = $sthSpecId->fetchrow_array){
		#>Filetring for obvious unambiguous positions
		my $ambiguity=0;
		foreach my $pepPosData (split(',',$phosPosStrg)) {
			my ($rank,$posStrg)=split(':',$pepPosData);
			$sthPep{$rank}->execute($queryID);
			my ($infoPep)=$sthPep{$rank}->fetchrow_array;
			my ($seq) = ($infoPep=~ /SEQ=([^,]+),/);
			my @acceptors=$seq=~/[$specificity]/g;
			my @phosphoPos=split(/\./,$posStrg);
			if (scalar @acceptors > scalar @phosphoPos) { # more acceptors than modifs -> phosphoRS
				$ambiguity=1;
				last;
			}
		}
		if ($ambiguity) {
			push @spectrumList,[$queryNum,$extSpectrumID];
		}
		else { # no need for phosphoRS
			$noAmbiguityQueries{$queryID}=$phosPosStrg;
		}
	}
	$sthSpecId->finish;
	foreach my $sth (values %sthPep){
		$sth->finish;
	}
	$dbh->disconnect;

	open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
	print FILESTAT "1.2/3 Fetching spectra from MSF file...\n";
	close FILESTAT;
	my $dbsqlite = DBI->connect("dbi:SQLite:$msfFileName","","",{PrintError=>1,RaiseError=>1});
	
	my $numSpectra=scalar @spectrumList;
	my @fracList;
	foreach my $i (1..10) {
		push @fracList,int($i*$numSpectra/10);
	}
	my $spCount=0; my $fracIdx=0;
    open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
	print FILESTAT "1.3/3 Extracting $numSpectra spectra from MSF file: 0% ($spCount/$numSpectra)\n";
	close FILESTAT;

	foreach my $refSpectrum (@spectrumList) {
		my ($queryNum,$extSpectrumID)=@{$refSpectrum};
		my ($parentFile,@spectrumInfos) = &promsMod::extractSpectrumMSF($dbsqlite, $extSpectrumID); #, $dir, 'no_user');
		my $refIons = $spectrumInfos[0];
		my $charge = $spectrumInfos[4];
		my $peaks = '';

		foreach my $value (@{$refIons}){
			$peaks .= $value->{X} . ':' . $value->{Y} . ',';
		}
		$peaks =~ s/,$//;

		$phosphoRS->addSpectrumInfo(Peaks => $peaks, Charge => $charge, Activation => $activationType, Query => $queryNum) or die "Cannot add spectrum info for query $queryNum";

		$spCount++;
		if ($spCount==$fracList[$fracIdx]) {
			$fracIdx++;
		}
        
        if ($spCount%50 == 0) {
			open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
			print FILESTAT "1.3/3 Extracting $numSpectra spectra from MSF file: ",(($fracIdx+1)*10),"% ($spCount/$numSpectra)\n";
			close FILESTAT;
		}
	}
	$dbsqlite->disconnect;
	sleep 3;
}

####>Running PhosphoRS<####
open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
print FILESTAT "2/3 Running PhosphoRS analysis\n";
close FILESTAT;

my $code = $phosphoRS->startAnalysis;
#if($code != 0){
#	print "<BR><FONT class=\"title2\" color=\"red\">Failed, code:$code</FONT><BR><BR>";
#	$phosphoRS->cleanFiles;
#	next;
#}


####>Updating sites position in DB<####
open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
print FILESTAT "3/3 Updating sites position\n";
close FILESTAT;

$dbh = &promsConfig::dbConnect('no_user');
my %sthUpQ;
foreach my $i (1..$maxRank){
	$sthUpQ{$i} = $dbh->prepare("UPDATE QUERY_VALIDATION SET INFO_PEP$i=? WHERE ID_QUERY=?");
}
my $sthUpQM1 = $dbh->prepare("UPDATE QUERY_MODIFICATION SET REF_POS_STRING=POS_STRING WHERE ID_MODIFICATION=$phosphoModID AND ID_QUERY=? AND PEP_RANK=? AND REF_POS_STRING IS NULL");
my $sthUpQM2 = $dbh->prepare("UPDATE QUERY_MODIFICATION SET POS_STRING=? WHERE ID_MODIFICATION=$phosphoModID AND ID_QUERY=? AND PEP_RANK=?");
my $sthQuery = $dbh->prepare("SELECT ID_QUERY,QUERY_NUM,$infoPepStrg FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$anaID AND QUERY_NUM>0 AND VALID_STATUS >=-3");
my $sthPhosphoPosStg = $dbh->prepare("SELECT POS_STRING,REF_POS_STRING FROM QUERY_MODIFICATION WHERE ID_QUERY=? AND PEP_RANK=? AND ID_MODIFICATION=$phosphoModID");
$sthQuery->execute;
while (my ($idQuery,$queryNum,@infoPep) = $sthQuery->fetchrow_array) {
	if ($noAmbiguityQueries{$idQuery}) {
		foreach my $pepPosData (split(',',$noAmbiguityQueries{$idQuery})) {
			my ($rank,$posStrg)=split(':',$pepPosData);
			my $infoPepStrg = $infoPep[$rank-1];
			$infoPepStrg .= "PRS=3;100;,";
			$sthUpQ{$rank}->execute($infoPepStrg,$idQuery);
			$countConfirmed++;
			$countTotal++;
		}
		next;
	}
	# Browsing each info pep #
	my $i=0;
	while( my $infoPepStrg = shift @infoPep){
		$i++;
		$sthPhosphoPosStg->execute($idQuery,$i);
		my ($positionStrg,$refpos)=$sthPhosphoPosStg->fetchrow_array;
		$positionStrg=$refpos if ($refpos && $positionStrg =~ /-\d/);# one site is ambigous
		if($infoPepStrg =~ /VMOD=([^,]+),/){# keep VMOD for RMOD
			my $vmodStrg = $1;
			my $rmodStrg;
			if ($infoPepStrg !~ /RMOD=/) {
				$rmodStrg = $vmodStrg;
			}
			#my $positionStrg;
			# Fetching all phosphorylation positions #
			my %positions;
			#while($vmodStrg =~ /Phospho \(([A-Z]+):([^\)]+)\)/g){
			#	$positionStrg .= ($positionStrg)? ".$2" : $2;
			#	$positions{$1} = $2;
			#}
			$positions{'STY'}=$positionStrg;
			next unless $positionStrg;
			my ($seq) = ($infoPepStrg=~ /SEQ=([^,]+),/);
			#next unless ($seq =~ /[STY].*[STY]/);makePhosphoRS.sh

			$countTotal++;

			# Fetching best isoform from PhosphoRS #
			#print "Query: $queryNum , i: $i , Positions: $positionStrg<BR>\n";
			if(my $isoformRef = $phosphoRS->getIsoforms($queryNum,$i)){
				my @isoforms = @{$isoformRef};
				my $bestIsoformRef = (sort {$b->[0] <=> $a->[0]} @isoforms )[0];
				my ($proba,$positionsRef) = @{$bestIsoformRef};
				$proba *= 100;
				my $bestPositions = join('.',@{$positionsRef});
				$infoPepStrg =~ s/PRS=[^,]+,//g;

				# Updating info pep data strings #
				if($bestPositions ne $positionStrg){
					if($proba >= $probThreshold){
						# Change Mascot positions to phosphoRS positions
						my ($sequence) = ($infoPepStrg =~ /SEQ=([^,]+)/);
						my @sequence = split(//,$sequence);
						my %newPositions;
						foreach my $position (@{$positionsRef}){
							my $aa = $sequence[$position-1];
							#my $phosphoType = ($aa eq 'S' or $aa eq 'T')? 'ST' : $aa; # assuming that S & T phosphorylations are always regrouped
							my $phosphoType = ($aa=~/[$specificity]/)? $specificity : $aa; # uses search specif (PP 19/12/14)
							push @{$newPositions{$phosphoType}}, $position;
						}
						$vmodStrg =~ s/ \+ Phospho[^\)]+\)//g;
						foreach my $phosphoType (keys %newPositions){
							my $posStrg = join '.', @{$newPositions{$phosphoType}};
							$vmodStrg .= " + Phospho ($phosphoType:$posStrg)";
						}
						$infoPepStrg =~ s/VMOD=[^,]+/VMOD=$vmodStrg/;
						# Store old positions
						$infoPepStrg .= "RMOD=$rmodStrg," if $rmodStrg;
						$sthUpQM1->execute($idQuery,$i); # Record original position
						$sthUpQM2->execute($bestPositions,$idQuery,$i);
						my $aaModPos;
						foreach my $aa ( keys %positions ){ # to keep phospho aa type (lost by PhosphoRS)
							$aaModPos .= "[$aa]";
							$aaModPos .= $positions{$aa};
						}
						$infoPepStrg .= "PRS=2;$proba;$aaModPos,"; # (status ; PRS score ; previous Mascot positions)
						$countChanged++;
					}
					else {
						# PRS != Mascot, but unchanged because PRS proba < threshold
						$infoPepStrg .= "PRS=1;$proba;$bestPositions,"; # (status ; PRS score ; PRS best positions)
						$countFlagged++;
					}
				}
				else {
					# PRS = Mascot
					if($proba >= $probThreshold){
						$infoPepStrg .= "PRS=3;$proba;,"; # same isoform for PRS and Mascot
						$countConfirmed++;
					}
					else {
						$infoPepStrg .= "PRS=0;$proba;,"; # same isoform but PRS score < threshold
						$countUnchanged++;
					}
				}
			}
			else {
			   $infoPepStrg .= "PRS=4;;,";
			}
			$sthUpQ{$i}->execute($infoPepStrg,$idQuery);
		}
	}
}
foreach my $sth (values %sthUpQ) {
	$sth->finish;
}
$sthUpQM1->finish;
$sthUpQM2->finish;
$sthPhosphoPosStg->finish;

my $paramFile = "$promsPath{valid}/ana_$anaID/PRSparam_ana_$anaID.txt";
open (PRSPARAM, ">$paramFile");
print PRSPARAM qq
|Threshold:$probThreshold%
Mass Tolerance:$massTolerance Da
Activation Type:$activationType|;
close PRSPARAM;

####>Updating validation history<####
my $paramStrg = "threshold:$probThreshold;massTolerance:$massTolerance;activationType:$activationType;";
&promsMod::updateAnalysisHistory($dbh,$anaID,$paramStrg,'prs',$userID);
$dbh->commit;
$dbh->disconnect;
		
rename("$currentPRSPath/$anaID\_$jobDir\_run.flag","$currentPRSPath/$anaID\_$jobDir\_end.flag"); # end flag file

open(FILESTAT,">>$fileStat");
print FILESTAT "PhosphoRS Ended.",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
close(FILESTAT);


sleep 30;

##>Cleaning tmp job directory
# unlink "$currentPRSPath/$anaID\_$jobDir\_error.txt";
# unlink "$currentPRSPath/$anaID\_$jobDir\_end.flag"; Do NOT delete




######################################> SUBROUTINES <########################################

sub launchAllJobs { # Globals: $phosphoRSPath,$currentPRSPath
	
	my ($masterJobCode,$refAnalysisList,$probThreshold,$massTolerance,$activationType,$overWrite)=@_;
	my $userID=$ENV{REMOTE_USER};
	my $numJobToRun=scalar @{$refAnalysisList};
	my %runningJobs;
	my $cgiUnixPath=`pwd`;
	$cgiUnixPath=~s/\/*\s*$//; # trim any trailing '/' & spaces
	my $dbh=&promsConfig::dbConnect;
	my $projectID = &promsMod::getProjectID($dbh,$refAnalysisList->[0],'ANALYSIS');
	my %cluster=&promsConfig::getClusterInfo;
    
	####>Cluster jobs<####
	if ($cluster{'on'}) {
		my $checkForEndedJobs=sub { # This sub can access all variables global to parent bloc code!!!
			foreach my $anaID (keys %runningJobs) {
				my $jobDir=$masterJobCode.'.'.$runningJobs{$anaID};
				my $jobPath="$phosphoRSPath/$jobDir";
				if (-e "$currentPRSPath/$anaID\_$jobDir\_end.flag" || $cluster{'checkError'}->("$currentPRSPath/$anaID\_$jobDir\_error.txt",['Cache disabled'])) {
					delete $runningJobs{$anaID}; # untrack job.
					$numJobToRun--;
				}
			}
		};

		####>Fetching amount of data to estimate maxMem
		my %peptidesInAna;
		my $sthNumPep=$dbh->prepare("SELECT ID_ANALYSIS,COUNT(*) FROM QUERY_VALIDATION WHERE QUERY_NUM>0 AND VALID_STATUS >=-3 AND ID_ANALYSIS IN (".join(',',@{$refAnalysisList}).") GROUP BY ID_ANALYSIS");
		$sthNumPep->execute;
		while (my ($anaID,$numPep)=$sthNumPep->fetchrow_array) {
			$peptidesInAna{$anaID}=$numPep;
		}
		$sthNumPep->finish;

		my %baseJobParameters=(
			numCPUs=>2,
			maxHours=>168, # 7 days
			pbsRunDir=>$cgiUnixPath,
			noWatch=>1 # do not wait for job to end
		);
		my $MAX_PARALLEL_JOBS=$cluster{'maxJobs'};
        
		####>Looping through job list
		my $curJobIdx=-1;
		MAIN_LOOP:while ($numJobToRun) {
	
			###>Parallel launch of up to $MAX_PARALLEL_JOBS jobs
			while (scalar (keys %runningJobs) < $MAX_PARALLEL_JOBS) {
	
				##>Also check other user/launches jobs
				my @runs=glob "$currentPRSPath/*_run.flag";
				last if scalar @runs >= $MAX_PARALLEL_JOBS;
				
				##>Ok to go
				$curJobIdx++;
				my $jobPos=$curJobIdx+1;
				my $anaID=$refAnalysisList->[$curJobIdx];
				my $jobDir=$masterJobCode.'.'.$jobPos;
				my $jobPath="$phosphoRSPath/$jobDir";
				
				my %jobParameters=%baseJobParameters;
				$jobParameters{jobName}="myProMS_phosphoRS_$anaID";
				$jobParameters{maxMem}=int(1 + $peptidesInAna{$anaID}/750).'Gb';
				my $commandString="$cgiUnixPath/analysePhospho.cgi $userID $jobDir $anaID $probThreshold $massTolerance $activationType $overWrite 2> $currentPRSPath/$anaID\_$jobDir\_error.txt";
				my ($pbsError, $pbsErrorFile, $jobClusterID) = $cluster{'runJob'}->($jobPath, $commandString, \%jobParameters);
				$runningJobs{$anaID}=$jobPos;
				
				unless ($dbh) { # reconnect in case disconnected
					$dbh=&promsConfig::dbConnect('no_user');
				}
                $dbh->do("INSERT INTO JOB_HISTORY (ID_JOB, ID_USER, ID_PROJECT, ID_JOB_CLUSTER, TYPE, JOB_STATUS, FEATURES, SRC_PATH, LOG_PATH, ERROR_PATH, STARTED) VALUES('$jobDir', '$userID', $projectID, 'C$jobClusterID', 'Phospho', 'Queued', 'ID_ANALYSIS=$anaID;PROB_THRESHOLD=$probThreshold;MASS_TOLERANCE=$massTolerance;ACTIVATION_TYPE=$activationType;OVERWRITE=$overWrite', '$jobPath', '$jobPath/status.out', '$currentPRSPath/$anaID\_$jobDir\_error.txt', NOW())");
                $dbh->commit;
                
				last MAIN_LOOP if $curJobIdx==$#{$refAnalysisList}; # no more jobs to launch => break MAIN_LOOP
				sleep 5;
			}
			if ($dbh) { # disconnect because wait for next launch is too long
				$dbh->disconnect;
				undef $dbh;
			}
			
			###>Wait for potential ended job before next MAIN_LOOP
			sleep 60; # $SIG{CHLD} is active during sleep
			&{$checkForEndedJobs};
	
		}
		###>Wait for last jobs to end
		while (scalar (keys %runningJobs) > 0) {
			sleep 60;
			&{$checkForEndedJobs};
		}
	}
	####>Local jobs<####
	else {
		$dbh->disconnect; # disconnect before fork

		$SIG{CHLD} = sub { # sub called when child communicates with parent (occur during parent sleep)
			local $!; # good practice. avoids changing errno.
			while (1) { # while because multiple children could finish at same time
				my $childPid = waitpid(-1,WNOHANG); # WNOHANG: parent doesn't hang during waitpid
				#my $childPid = waitpid(-1,POSIX->WNOHANG); # WNOHANG: parent doesn't hang during waitpid
				last unless ($childPid > 0); # No more to reap.
				delete $runningJobs{$childPid}; # untrack child.
			}
		};
		
		my $MAX_PARALLEL_JOBS=&promsConfig::getMaxParallelJobs;

		####>Looping through job list
		my $curJobIdx=-1;
		MAIN_LOOP:while (1) {
	
			###>Parallel launch of up to $MAX_PARALLEL_JOBS jobs
			while (scalar (keys %runningJobs) < $MAX_PARALLEL_JOBS) {
	
				##>Also check other user/launches jobs
				my $numProc=`ps -edf | grep analysePhospho.cgi | grep ' 2> ' | grep -v grep | wc -l`;
				chomp($numProc);
				last if $numProc >= $MAX_PARALLEL_JOBS;
	
				##>Ok to go
				$curJobIdx++;
	
				my $childPid = fork;
				unless ($childPid) { # child here
					my $jobPos=$curJobIdx+1;
					my $jobDir=$masterJobCode.'.'.$jobPos;
					my $anaID=$refAnalysisList->[$curJobIdx];
                    my $jobPath="$phosphoRSPath/$jobDir";
                    
					my $dbh=&promsConfig::dbConnect('no_user');
                    $dbh->do("INSERT INTO JOB_HISTORY (ID_JOB, ID_USER, ID_PROJECT, ID_JOB_CLUSTER, TYPE, JOB_STATUS, FEATURES, SRC_PATH, LOG_PATH, ERROR_PATH, STARTED) VALUES('$jobDir', '$userID', $projectID, 'L$$', 'Phospho', 'Queued', 'ID_ANALYSIS=$anaID;PROB_THRESHOLD=$probThreshold;MASS_TOLERANCE=$massTolerance;ACTIVATION_TYPE=$activationType;OVERWRITE=$overWrite', '$jobPath', '$jobPath/status.out', '$currentPRSPath/$anaID\_$jobDir\_error.txt', NOW())");
                    $dbh->commit;
					$dbh->disconnect;
					system "$cgiUnixPath/analysePhospho.cgi $userID $jobDir $anaID $probThreshold $massTolerance $activationType $overWrite 2> $currentPRSPath/$anaID\_$jobDir\_error.txt";
					exit;
				}
				$runningJobs{$childPid}=$curJobIdx;
                
				last MAIN_LOOP if $curJobIdx==$#{$refAnalysisList}; # no more jobs to launch => break MAIN_LOOP
				sleep 2;
			}
	
			###>Wait for potential ended job before next MAIN_LOOP
			sleep 60; # $SIG{CHLD} is active during sleep
	
		}
		###>Wait for last child processes to end
		while (scalar (keys %runningJobs) > 0) {
			sleep 60; # $SIG{CHLD} is active during sleep
		}
	}
}

sub deletePRS{
    my ($dbh,$anaID,$verbose) = @_;

	print "Deleting previous PhosphoRS analysis..." if $verbose;
    ## Restore DB entries ##
	my ($maxRank) = $dbh->selectrow_array("SELECT MAX_RANK FROM ANALYSIS WHERE ID_ANALYSIS=$anaID");
	my @infoPepList = map { "INFO_PEP$_" } (1..$maxRank);
	my $infoPepStrg = join(', ',@infoPepList);

	# table QUERY_MODIFICATION #
	my ($phosphoModID)=$dbh->selectrow_array("SELECT ID_MODIFICATION FROM MODIFICATION WHERE UNIMOD_ACC=21");
	$dbh->do("UPDATE QUERY_MODIFICATION M INNER JOIN QUERY_VALIDATION V ON M.ID_QUERY=V.ID_QUERY AND V.ID_ANALYSIS=$anaID AND M.ID_MODIFICATION=$phosphoModID AND M.REF_POS_STRING IS NOT NULL SET M.POS_STRING=M.REF_POS_STRING");
	$dbh->do("UPDATE QUERY_MODIFICATION M INNER JOIN QUERY_VALIDATION V ON M.ID_QUERY=V.ID_QUERY SET M.REF_POS_STRING=NULL WHERE V.ID_ANALYSIS=$anaID AND M.ID_MODIFICATION=$phosphoModID");

	# table QUERY_VALIDATION #
    my $sthInfoPep = $dbh->prepare("SELECT ID_QUERY,$infoPepStrg FROM QUERY_VALIDATION WHERE ID_ANALYSIS=? AND QUERY_NUM>0");
    my %sthUpInfoPep;
    for(my $i=1; $i<=$maxRank; $i++){
		$sthUpInfoPep{$i} = $dbh->prepare("UPDATE QUERY_VALIDATION SET INFO_PEP$i=? WHERE ID_QUERY=?");
    }
    $sthInfoPep->execute($anaID);
    while(my ($idQuery,@infoPep) = $sthInfoPep->fetchrow_array){
		my $i = 0;
		while(my $infoPep = shift @infoPep){
			$i++;
			if($infoPep =~ /VMOD=([^,]+Phospho[^,]+),/){
				my $vmodStrg = $1;
				if($infoPep =~ /PRS=(\d);/){
					if($1 == 2){
						# Restore old positions
						my ($aaModPos) = ($infoPep =~ /PRS=\d;[^;]+;([^,]+)/);
						my $restoredPhosphoStrg = '';
						while($aaModPos =~ /\[(\w+)\]([^\[,]+)/g){
							$restoredPhosphoStrg .= " + Phospho ($1:$2)";
						}
						$vmodStrg =~ s/ \+ Phospho \([^\)]+\)//g;
						$vmodStrg .= "$restoredPhosphoStrg";

						$infoPep =~ s/VMOD=[^,]+/VMOD=$vmodStrg/;
						$infoPep =~ s/RMOD=[^,]+,//; # Deleting PRS data
					}
					$infoPep =~ s/PRS=[^,]+,//; # Deleting PRS data
					$sthUpInfoPep{$i}->execute($infoPep,$idQuery);
				}
			}
		}
    }
    $sthInfoPep->finish;
    foreach my $i (keys %sthUpInfoPep){
		$sthUpInfoPep{$i}->finish;
    }

    # Update history #
    $dbh->do("DELETE FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$anaID AND VAL_TYPE='prs' AND STATUS<1");
    $dbh->do("UPDATE VALIDATION_HISTORY SET STATUS=-1 WHERE ID_ANALYSIS=$anaID AND VAL_TYPE='prs' AND STATUS=1");
    $dbh->commit;

    # Delete PRS files #
    #my ($dataFileName) = $dbh->selectrow_array("SELECT DATA_FILE FROM ANALYSIS WHERE ID_ANALYSIS=$anaID");
    #$dataFileName =~ s/\.\w{3}//;
    my $dir = "$promsPath{valid}/ana_$anaID";
    #unlink "$dir/PRS_$dataFileName.xml" or warn "Cannot delete PRS file: $!";
    #unlink "$dir/PRSparam_$dataFileName.txt" or warn "Cannot delete PRS file: $!";
	#unlink "$dir/PRS_ana_$anaID.xml" or warn "Cannot delete PRS output file: $!";
	#unlink "$dir/PRSparam_ana_$anaID.txt" or warn "Cannot delete PRS parameter file: $!";
	unlink glob "$dir/PRSparam_*.txt"; # old & new param naming
	unlink glob "$dir/PRS_*"; # old & new result.xml naming + status.txt
	#unlink "$dir/PRSerrors.txt" if -e "$dir/PRSerrors.txt";
	print " Done." if $verbose;
}

####>Revision history<####
# 2.1.0 [FEATURE] Handle reversion of multiple Analyses at once (PP 02/08/21)
# 2.0.9 [BUGFIX] Multiple fixes including ended jobs detection in cluster context (PP 21/04/21) 
# 2.0.8 [FEATURE] Clears PhosphoRS error & MySQL8 compatibility (PP 12/02/21)
# 2.0.7 [UPDATE] Changed JOB_HISTORY.STATUS to JOB_HISTORY.JOB_STATUS (PP 28/08/20)
# 2.0.6 [MINOR] Added project selection when opening monitor jobs windows (VS 02/09/20)
# 2.0.5 [BUGFIX] Handles analysis imported from splitted MSF files (VS 20/04/20)
# 2.0.4 [CHANGES] Use new job monitoring window opening parameters (VS 18/11/19)
# 2.0.3 [MODIF] Switch from watchPhosphoAnalyses to monitorJobs script (VS 21/10/19)
# 2.0.2 [ENHANCEMENT] Add insertion into JOB_HISTORY table for job monitoring (VS 08/10/19)
# 2.0.1 Added 10% step progression for msf spectra extraction (PP 07/01/19)
# 2.0.0 Major code rewrite for full background run support with call of watchPhosphoAnalyses.cgi (PP 08/11/18)
# 1.1.3 Minor bug fix in call for a new PhosphoRS analysis (PP 15/06/18)
# 1.1.2 Modification to avoid positionStg stored in INFO_PEP (GA 31/01/18)
# 1.1.1 Uses xxx_ana_&lt;anaID&gt;.xxx instead of data file name for PRS parameter and output files (PP 21/03/17)
# 1.1.0 Add a fork for big xml files that take too long with XMLin in phosphoRS.pm (GA 11/03/16)
# 1.0.9 Minor change (PP 26/08/15)
# 1.0.8 Compatibility update du to new output of &promsMod::extractSpectrumMSF v3.3.9 (PP 18/02/15)
# 1.0.7 Full compatibily with MSF files, updates table QUERY_MODIFICATION & speed improvement (PP 23/12/14)
# 1.0.6 Added 'use XML::Simple' (PP 13/11/13)
# 1.0.5 Added file format to phosphoRS parameters (FY 30/09/13)
# 1.0.4 Using RMOD tag in INFO_PEP string to store previous positions (FY 11/04/13)
# 1.0.3 New status (4) if no phospho-isoform for a given peptide (FY 06/11/12)
# 1.0.2 Management of any kind of phosphorylation (FY 26/07/12)
# 1.0.1 Management of Proteome Discoverer files +<BR>Add phosphopeptides filtering (FY 24/05/12)
# 1.0.0 New script processing PhosphoRS analyses (FY 28/03/12)