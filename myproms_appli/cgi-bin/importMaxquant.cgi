#!/usr/local/bin/perl -w

################################################################################
# importMaxquant.cgi       2.1.12                                              #
# Component of site myProMS Web Server                                         #
# Authors: P. Poullet, G. Arras, S. Liva (Institut Curie)                      #
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
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use POSIX qw(strftime :sys_wait_h);  # Core module
#use IO::Uncompress::Gunzip qw(gunzip);
use Cwd; #?
#use XML::Simple;
use File::Path qw(rmtree); # remove_tree
use File::Copy qw(copy move); # Core module
use File::Spec::Functions qw(splitpath); # Core module
use promsConfig;
use promsMod;
#exit;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my ($color1,$color2)=&promsConfig::getRowColors;
my $MAX_DB=3;

my $experimentID=&promsMod::cleanNumericalParameters(param('ID'));
my $projectID=&promsMod::cleanNumericalParameters(param('id_project'));
my $action=param('ACT') || 'form';

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);

if ($action eq 'form') {

	my (%DBlist,%isCrapDB);

	my $sthDB=$dbh->prepare("SELECT D.ID_DATABANK,D.NAME,VERSION_NAME,FASTA_FILE,DT.NAME,IS_CRAP FROM DATABANK D,DATABANK_TYPE DT WHERE D.ID_DBTYPE=DT.ID_DBTYPE AND USE_STATUS='yes'");
	$sthDB->execute;
	while (my ($dbID,$name,$version,$fastaFile,$dbankType,$isCrap)=$sthDB->fetchrow_array) {
		my $dbSource;
		if ($fastaFile=~/:/) {
			my ($mascotServer,$dbankDir,$fileName)=split(':',$fastaFile);
			$dbSource=$mascotServer;
			$version=$dbankDir;
		}
		else {
			if (!-e "$promsPath{banks}/db_$dbID/$fastaFile") {next;}
			$dbSource='Local';
		}
		$DBlist{$dbSource}{$dbID}=$name;
		$DBlist{$dbSource}{$dbID}.=" ($version)" if $version;
		$DBlist{$dbSource}{$dbID}.=" [$dbankType]";
		$isCrapDB{$dbID}=$isCrap;
	}
	$sthDB->finish;
	my $databaseString="<OPTION selected value=\"\">-= Select =-</OPTION>\n";
	my $databaseContString="<OPTION selected value=\"\">*Contaminants not searched*</OPTION>\n";

	foreach my $dbSource (sort{lc($a) cmp lc($b)} keys %DBlist) {
		$databaseString.="<OPTGROUP label=\"$dbSource:\">\n";
		foreach my $dbID (sort{lc($DBlist{$dbSource}{$a}) cmp lc($DBlist{$dbSource}{$b})} keys %{$DBlist{$dbSource}}) {
			$databaseString.="<OPTION value=\"$dbID\">$DBlist{$dbSource}{$dbID}</OPTION>\n";
			$databaseContString.="<OPTION value=\"$dbID\">$DBlist{$dbSource}{$dbID}</OPTION>\n" if $isCrapDB{$dbID};
		}
		$databaseString.="</OPTGROUP>\n";
		$databaseContString.="</OPTGROUP>\n";
	}

	$dbh->disconnect;

	print qq
|<HTML>
<HEAD>
<TITLE>Select MaxQuant files</TITLE>
<LINK rel=stylesheet href="$promsPath{html}/promsStyle.css" type=text/css>
<SCRIPT type="text/javascript">
|;
	&promsMod::browseDirectory_JavaScript if $promsPath{'shared'};
	print qq
|function updateDatabankSelection(dbNum,optIdx) {
	if (dbNum==$MAX_DB) return;
	var dbNext=dbNum+1;
	if (optIdx) { // selection
		document.getElementById('dbDIV_'+dbNext).style.display='';
		document.getElementById('db_'+dbNext).disabled=false;
	}
	else {
		for (var d=dbNext; d<=$MAX_DB; d++) {
			document.getElementById('dbDIV_'+d).style.display='none';
			document.getElementById('db_'+d).disabled=true;
			document.getElementById('db_'+d).options.selectedIndex=0;
		}
	}
}
function cancelAction() {
	window.location="./processAnalyses.cgi?ID=$experimentID";
}
function checkFileForm(importForm) {
	if (document.getElementById('formSubmit').disabled==true) return false; // in case "Enter" key is pressed after 1st submission
	var missingFiles=false;
	if (importForm.uplArch.files.length == 0) { // nothing in upload option
		missingFiles=true;
		if (importForm.sharedDirFiles) { // shared dir available
			if (importForm.sharedDirFiles.length) { // multiple files found
				var numFiles=0;
				for (let i=0; i<importForm.sharedDirFiles.length; i++) {
					if (importForm.sharedDirFiles[i].checked) {
						if (importForm.sharedDirFiles[i].value.match(/\\.(zip\|gz)\$/)) {
							missingFiles=false;
							break;
						}
						else if (importForm.sharedDirFiles[i].value.match(/(evidence\|peptides)\\.txt\$/)) {
							if (++numFiles==2) {
								missingFiles=false;
								break;
							}
						}
					}
				}
			}
			else if (importForm.sharedDirFiles.checked && importForm.sharedDirFiles.value.match(/\\.(zip\|gz)\$/)) {missingFiles=false;} // single file found
		}
	}
	if (missingFiles) {
		alert('ERROR: Missing MaxQuant data file(s).');
		return false;
	}
	if (!importForm.databank1.value) {
		alert('ERROR: Select at least 1 sequence databank to be used !');
		return false;
	}
	if (!importForm.mgType.value) {
		alert('ERROR: Select a rule for Match groups generation !');
		return false;
	}
	if (importForm.uplArch.files.length) {
		document.getElementById('formSubmit').disabled=true;
		document.getElementById('waitDIV').style.display='block';
	}
	if (importForm.sharedDirFiles.length) {
		document.getElementById('formSubmit').disabled=true;
	}
	return true;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Select MaxQuant file to be imported</FONT>
<BR><BR>
<FORM name="mqListForm" action="./importMaxquant.cgi" method="post" onsubmit="return checkFileForm(this);" enctype="multipart/form-data">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="ID" value="$experimentID">
<INPUT type="hidden" name="ACT" value="upload">
<TABLE>
<TR><TD colspan=2>
	<TABLE bgcolor=$color2 style="min-width:800px">
	<TR>
		<TH align=right valign=top>MaxQuant files* :</TH>
		<TD bgcolor=$color1 nowrap>
		<TABLE cellpadding=0>
			<TR><TH align="right">Upload archive file**:</TH><TD><INPUT type="file" name="uplArch"></TD></TR>
|;
	if ($promsPath{'shared'}) {
		print qq
|			<TR><TH align="right">or</TH><TH></TH></TR>
			<TR><TH align="right" valign=top nowrap>Use shared directory:</TH><TD><DIV id="sharedDirDIV" style="width:600px;max-height:300px;overflow:auto;border:2px solid $color2">
|;
		&promsMod::browseDirectory_getFiles($promsPath{'shared'},'sharedDirFiles',{fileMatch=>qr/\.(txt|xml|gz|zip)$/i}); # List of handled file extensions
		print qq
|		</DIV></TD>
	</TR>
|;
	}
	print qq
|		</TABLE>
		</TD>
	</TR>
	<TR>
		<TH align=right valign="top" nowrap>&nbsp;Databank(s) used for search :</TH>
		<TD bgcolor=$color1 nowrap>|;
	foreach my $d (1..$MAX_DB) {
		my ($dispStrg,$disabStrg)=($d==1)? ('','') : ('style="display:none"',' disabled');
		print qq
|<DIV id="dbDIV_$d" $dispStrg><B>#$d:</B><SELECT id="db_$d" name="databank$d" style="width:550px" onchange="updateDatabankSelection($d,this.options.selectedIndex)">$databaseString</SELECT></DIV>
|;
	}
	print qq
|</TD>
	</TR>
	<TR>
		<TH align=right valign="top">Contaminants :</TH>
		<TD bgcolor=$color1 nowrap><INPUT type="checkbox" name="excludeCON" value="1" onclick="document.getElementById('dbDIV_con').style.display=(this.checked)? 'none' : ''"><B>Exclude from protein list<BR><DIV id="dbDIV_con">&nbsp;<B>Databank:</B><SELECT name="databankCON" style="max-width:600px">$databaseContString</SELECT></DIV></TD>
	</TR>
	<TR>
		<TH align=right valign="top">Quantification :</TH>
		<TD bgcolor=$color1 nowrap><INPUT type="checkbox" name="protQuantif"  value="1" checked><B>Import protein quantification data</B></TD>
	</TR>
	<TR>
		<TH align=right>Match groups rule :</TH>
		<TD bgcolor=$color1 nowrap><SELECT name="mgType"><OPTION value="">-= Select =-</OPTION><OPTION value="MaxQuant">MaxQuant</OPTION><OPTION value="myProMS">myProMS</OPTION></SELECT></TD>
	</TR>
	<TR><TH colspan=2>
		<INPUT type="submit" name="save" id="formSubmit" value=" Proceed ">
		&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="cancelAction();"></TH>
	</TR>
</TABLE></TD></TR>
<TR>
	<TD valign=top nowrap>
	<DIV style="float:left">*</DIV><DIV style="float:left;height:100px"><B>mqpar.xml</B> <I>(optional)</I><BR>and files from MaxQuant 'combined/txt' directory:
	<UL style="margin:0">
		<LI><B>evidence.txt</B></LI>
		<LI><B>msms.txt</B> <I>(optional. For MS/MS spectrum display)</I></LI>
		<LI><B>peptides.txt</B></LI>
		<LI><B>proteinGroups.txt</B> <I>(required only if protein quantification is needed)</I></LI>
		<LI><B>parameters.txt</B> <I>(required only if <B>mqpar.xml</B> is not provided)</I></LI>
		<LI><B>summary.txt</B> <I>(required only if <B>mqpar.xml</B> is not provided)</I></LI>
	</UL>
	</DIV>
	<DIV>&nbsp;**zip or gz archive&nbsp;</DIV>
	</TD>
</TR>
</TABLE>
</FORM>
<DIV id="waitDIV" style="display:none">
<BR><BR>
<FONT class="title2">Uploading file. Please wait...</FONT>
</DIV>
</CENTER>
</BODY>
</HTML>
|;
	exit;
}

### Import Archive file
print qq
|<HTML>
<HEAD>
<TITLE>Preparing MaxQuant Import</TITLE>
<LINK rel=stylesheet href="$promsPath{html}/promsStyle.css" type=text/css>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER><FONT class="title1">Preparing MaxQuant Import</FONT></CENTER>
<BR>
|;

my $evidenceFile='evidence.txt'; # Contains all peptide identified and peptide quantitation
my $peptideFile='peptides.txt'; # Contains all peptide identified and peptide quantitation
my $proteinGroupsFile='proteinGroups.txt'; # Contains all quantitation information at protein level (LFQ, iBAQ, Intensity,...)
my $msmsFile='msms.txt'; # Contains m/z and intensities information
my $summaryFile='summary.txt'; # Contains some redundant information found in mqpar.xml (varmod, label type, experimental design)
my $parametersFile='parameters.txt'; # Contains some redundant information found in mqpar.xml (fixmod, peptide parameters, software version)

#mkdir "$promsPath{tmp}/MaxQuant" unless -e "$promsPath{tmp}/MaxQuant";
#&promsMod::cleanDirectory("$promsPath{tmp}/MaxQuant",'1d'); # [1 day]  also created if does not exists yet
#my $tmpFilesDir="$promsPath{tmp}/MaxQuant/".(strftime("%Y%m%d%H%M%S",localtime)).'_'.$userID;
#rmtree $tmpFilesDir if -e $tmpFilesDir; # clean directory
#mkdir $tmpFilesDir;

mkdir "$promsPath{tmp}/quantification" unless -e "$promsPath{tmp}/quantification";
my $currentQuantifDir="$promsPath{tmp}/quantification/current";
mkdir $currentQuantifDir unless -e $currentQuantifDir;
my $jobID=strftime("%Y%m%d%H%M%S",localtime);
while (-e "$promsPath{tmp}/quantification/$jobID") { # to prevent multi-user launch collision
	sleep 2;
	$jobID=strftime("%Y%m%d%H%M%S",localtime);
}
my $tmpFilesDir="$promsPath{tmp}/quantification/$jobID";
mkdir $tmpFilesDir;

$dbh->disconnect;  # Must disconnect before copying files and reconnect after, otherwise loose connection
my $importPepQuantif=1;
my $importProtQuantif=param('protQuantif');

print "<BR><FONT class=\"title3\">Retrieving File(s)...";
print "\0" x 1024, "<BR>";  # Fill buffer for the first print
my $newFile;
my $numFiles=0;
my $archiveFile;
if (param('uplArch')) {
	my (undef,$path,$uplFile)=splitpath(param('uplArch'));
	$newFile="$tmpFilesDir/$uplFile";
	move(tmpFileName(upload("uplArch")),$newFile); # name of temp file being uploaded
	$numFiles=1;
	$archiveFile=$uplFile;
}
else { # from shared directory
	foreach my $sharedFile (param('sharedDirFiles')) {
		my (undef,$path,$fileName)=splitpath($sharedFile);
		$newFile="$tmpFilesDir/$fileName";
		&copyAndPrint("$promsPath{shared}/$sharedFile", "$newFile");
		if ($fileName =~/\.(gz|zip)\Z/) { # keep only archive if extra files
			$numFiles=1;
			$archiveFile=$fileName;
			last;
		}
		else {$numFiles++;}
	}
}
if ($numFiles==1 && $newFile !~ /\.(gz|zip)\Z/) {
	print "[<FONT color=\"#DD0000\">ERROR: The archive type is not recognized (Use zip or gz only).</FONT>]\n";
	rmtree $tmpFilesDir;
	exit;
}
print "<BR> Done.</FONT><BR>";

###>Writing parameter file
my @databankIDs;
foreach my $d (1..$MAX_DB) {
	next unless param("databank$d");
	push @databankIDs,&promsMod::cleanNumericalParameters(param("databank$d"));
}
my $excludeCON=param('excludeCON') || 0;
my $contaminantDB;
if (!$excludeCON && param('databankCON')) {
	$contaminantDB=&promsMod::cleanNumericalParameters(param('databankCON'));
	#push @databankIDs,$contaminantDB;
}
my $matchGroupType=param('mgType') || 'myProMS';
open (INFO,">$tmpFilesDir/quantif_info.txt");
print INFO qq
|USER=$userID
TYPE=IMPORT_MQ
ID_EXPERIMENT=$experimentID
ID_PROJECT=$projectID
DATABANKS=@databankIDs
EXCLUDE_CON=$excludeCON
MG_TYPE=$matchGroupType
PROT_QUANTIF=$importProtQuantif
|;
print INFO "ID_CONT_DB=$contaminantDB\n" if $contaminantDB;
if ($archiveFile) {
	my (undef,$path,$archFile)=splitpath($newFile);
	print INFO "ARCHIVE=$archFile\n";
}
close INFO;

open(FLAG,">$currentQuantifDir/$experimentID\_$jobID\_wait.flag"); # flag file used by watchQuantif
print FLAG "#";
close FLAG;

$dbh=&promsConfig::dbConnect;
$dbh->do("INSERT INTO JOB_HISTORY (ID_JOB, ID_USER, ID_PROJECT, TYPE, JOB_STATUS, FEATURES, SRC_PATH, LOG_PATH, ERROR_PATH, STARTED) VALUES('$jobID', '$userID', $projectID, 'Import [MaxQuant]', 'Queued', 'ID_EXPERIMENT=$experimentID', '$tmpFilesDir', '$tmpFilesDir/status_$experimentID.out', '$currentQuantifDir/$experimentID\_$jobID\_error.txt', NOW())");
$dbh->commit;
$dbh->disconnect;

###>Forking
my $childPid = fork;
unless ($childPid) { # child here
	#>Disconnecting from server
	open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
	open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
	#open STDERR, ">>$promsPath{logs}/launchQuantification.log";
	#system "./launchQuantifications.pl single $userID $jobID MAXQUANT $experimentID";
	
	my $dbh=&promsConfig::dbConnect;
	my $error;
	my %cluster=&promsConfig::getClusterInfo;
#$cluster{'on'}=0; # TEMP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	my $eviSize=("$tmpFilesDir/$evidenceFile")? -s "$tmpFilesDir/$evidenceFile" : 5 * (-s $archiveFile);
	my $onCluster=0;
	if ($cluster{'on'} && $eviSize > 50000000) { # 50 Mb
		$onCluster=1;
		my $cgiUnixDir=`pwd`;
		$cgiUnixDir=~s/\/*\s*$//;
		# cd is required for script to find myproms .pm files!!!!
		my $commandString="export LC_ALL=\"C\"; cd $cgiUnixDir; $cluster{path}{perl}/perl runMaxquantImport.pl $experimentID $jobID $onCluster 2> $currentQuantifDir/$experimentID\_$jobID\_error.txt";
		my $maxMem=int(2 + (15 * $eviSize / 1024**3)); # 2 + 15 * size in Gb
		my %jobParameters=(
			maxMem=>$maxMem.'Gb',
			numCPUs=>1,
			maxHours=>24,
			jobName=>"myProMS_MxQtImport_$experimentID\_$jobID",
			pbsRunDir=>$cgiUnixDir,
			commandBefore=>"mv $currentQuantifDir/$experimentID\_$jobID\_wait.flag $currentQuantifDir/$experimentID\_$jobID\_run.flag" # run flag file
		);
		my ($pbsError,$pbsErrorFile,$jobClusterID)=$cluster{'runJob'}->($tmpFilesDir,$commandString,\%jobParameters);
		
		# Add cluster job id to current job in DB
		$dbh->do("UPDATE JOB_HISTORY SET ID_JOB_CLUSTER='C$jobClusterID' WHERE ID_JOB='$jobID'");
		$dbh->commit;
		
		if ($pbsError) { # move PBS error message to job error file
			system "cat $pbsErrorFile >> $currentQuantifDir/$experimentID\_$jobID\_error.txt";
			$error=1;
		}
	}
	else { # Local launch (web server): No cluster OR XICMCQ|SIN|EMPAI
		rename("$currentQuantifDir/$experimentID\_$jobID\_wait.flag","$currentQuantifDir/$experimentID\_$jobID\_run.flag"); # run flag file

		# Add process PID to current job in DB
		$dbh->do("UPDATE JOB_HISTORY SET ID_JOB_CLUSTER='L$$' WHERE ID_JOB='$jobID'");
		$dbh->commit;
		
		system "./runMaxquantImport.pl $experimentID $jobID $onCluster 2> $currentQuantifDir/$experimentID\_$jobID\_error.txt";

		if (-s "$currentQuantifDir/$experimentID\_$jobID\_error.txt") {
			$error=1;
		}
		
	} # end of local launch
	if ($error) { # Update Quantifications status
		my $sthQ=$dbh->prepare("SELECT DISTINCT Q.ID_QUANTIFICATION FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ,ANALYSIS A,SAMPLE S
								WHERE Q.ID_QUANTIFICATION=AQ.ID_QUANTIFICATION AND AQ.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE
								AND S.ID_EXPERIMENT=? AND Q.STATUS=0");
		my $sthUpQ=$dbh->prepare("UPDATE QUANTIFICATION SET STATUS=-1,UPDATE_DATE=NOW() WHERE ID_QUANTIFICATION=?");
		$sthQ->execute($experimentID);
		while (my ($quantifID)=$sthQ->fetchrow_array) {
			$sthUpQ->execute($quantifID);
		}
		$sthQ->finish;
		$sthUpQ->finish;
		$dbh->commit;
	}
	
	unlink "$currentQuantifDir/$experimentID\_$jobID\_run.flag";
	sleep 5;
	unless (-s "$currentQuantifDir/$experimentID\_$jobID\_error.txt") {
		if (-e "$tmpFilesDir/ana2map.txt") { # New identifiers to map
			open(MAPPING,"$tmpFilesDir/ana2map.txt") || die "Error while opening ana2map.txt: $!";
			my $anaStrg=<MAPPING>;
			system "./mapProteinIdentifiers.pl $userID $anaStrg";
		}

		#unlink "$currentQuantifDir/$quantifItemID\_$jobID\_end.flag";
		unlink "$currentQuantifDir/$experimentID\_$jobID\_error.txt";
	}
	
	$dbh->disconnect;
	exit;
} # end of child


print qq
|<BR><FONT class="title3">Launching MaxQuant import in background.<BR>
Progress can be tracked in "Monitor Quantifications" window.<BR><BR>
This page will refresh itself in a few seconds.</FONT>
<SCRIPT type="text/javascript">
var monitorJobsWin=window.open("$promsPath{cgi}/monitorJobs.cgi?filterType=Import [MaxQuant]&filterDateNumber=1&filterDateType=DAY&filterStatus=Queued&filterStatus=Running&filterProject=$projectID",'monitorJobsWindow','width=1200,height=500,scrollbars=yes,resizable=yes');
monitorJobsWin.focus();
</SCRIPT>
|;
sleep 5;
print qq
|<SCRIPT type="text/javascript">
// top.promsFrame.selectedAction='summary';
parent.optionFrame.selectOption(parent.optionFrame.document.getElementById('summary')); // refresh optionFrame with summary option
</SCRIPT>
</BODY>
</HTML>
|;


sub copyAndPrint {
	my ($sourceFile, $targetFile) = @_;
	my (undef, $path, $fileName) = splitpath($sourceFile);

	my $sourceFileSize = (stat $sourceFile)[7];
	my $sourceSizeGo = $sourceFileSize / (1024 ** 3);
	
	if ($sourceSizeGo < 0.1) {
		copy($sourceFile, $targetFile);
	}
	else { # fork and wait
		my $maxCopyTime = ($sourceSizeGo > 1)? $sourceSizeGo * 20 : 20;  # max 20min or 20min per Go
		my $targetFileSize = 0;
		my $loopNb = 0;  # To avoid infinite loop
	
		# Fork to copy files (makes it possible to keep the page active and avoid timeout for big files)
		my $childPid = fork;
		unless ($childPid) { # child here
			#>Disconnecting from server
			open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
			open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
			open STDERR, ">>$promsPath{logs}/importMaxQuant.log";
			copy($sourceFile, $targetFile);
			exit;
		}
		print "<BR>&nbsp;&nbsp;&nbsp;&nbsp;- Storing $fileName.";
		sleep 1;  # To avoid sleeping 10 sec after in most cases
	
		# Start tracking target file and check if the copy is actually being done
		while ($loopNb < 6 * $maxCopyTime) {  # 10sec loops -> loopNb = 6 * nb of minutes passed
			my $res = waitpid($childPid, WNOHANG); # WNOHANG (from POSIX ":sys_wait_h"): parent doesn't hang during waitpid
			last if $res; # child has ended
	
			sleep 10;
			print ".";
			print "<BR>" if ($loopNb > 30);  # End line every 5 minutes
			$loopNb++;
		}
		if ($loopNb >= 6 * $maxCopyTime) {
			print "Killing child process $childPid" if ($childPid);
			kill "SIGKILL", $childPid if ($childPid);
			die "Retrieving file is taking too long or process died silently before completion";
		} else {
			$targetFileSize = (-e "$targetFile") ? (stat $targetFile)[7] : 0;
			if ($targetFileSize == $sourceFileSize) {
				print("<BR>&nbsp;&nbsp;&nbsp;&nbsp;Done.");
			} else {
				die "Copy of file was not done properly. Exiting...";
			}
		}
	}
}


####>Revision history<####
# 2.1.12 [ENHANCEMENT] Minor change to &amp;copyAndPrint to prevent fork on small files (PP 07/06/21)
# 2.1.11 [BUGFIX] Print on web page during copy from scratch to avoid timeout for large size files (VL 19/05/21)
# 2.1.10 [BUGFIX] Pass cluster parameter to runMaxquantImport for proper management in &getProtInfo (VL 19/05/21)
# 2.1.9 [UPDATE] Changed JOB_HISTORY.STATUS to JOB_HISTORY.JOB_STATUS (PP 28/08/20)
# 2.1.8 [MINOR] Added project selection when opening monitor jobs windows (VS 02/09/20)
# 2.1.7 [CHANGES] Use new job monitoring window opening parameters (VS 18/11/19)
# 2.1.6 [ENHANCEMENT] Handles new job monitoring system (VS 16/10/19)
# 2.1.5 [FEATURE] Improved cluster memory calculation & removed 'PARAM_GR=0' in quantif_info.txt (PP 23/08/19)
# 2.1.4 [BUGFIX] missing experimentID argument and archive parameter for local job launch (PP 10/07/19)
# 2.1.3 Better cluster memory estimation for small job (PP 08/07/19)
# 2.1.2 [Fix] mispelled runMaxQuantImport.pl to runMaxquantImport.pl (PP 02/07/19)
# 2.1.1 [Fix] bug full path for archive file (PP 20/06/19)
# 2.1.0 Improved error management & cluster resources calculation (PP 11/06/19)
# 2.0.0 Restricted to form submission & data file transfert. Import is performed in background by runMaxquantImport.cgi (PP 31/05/19)
# 1.3.11 Improved contaminents exclusion & memory management (PP 17/05/19)
# 1.3.10 Improved handling of TMT/iTRAQ and fix $massErrorDa bug recomputing for isobaric data (PP 10/04/19)
# 1.3.9 Minor modification for TMT as Terminal label has to be parsed as well (GA 25/02/19)
# 1.3.8 Minor modification due to mqpar.xml IsobaricLabelInfo/TMT changes in 1.6.3.4 maxQuantVersion (GA 20/02/19)
# 1.3.7 parameters.txt & summary.txt files are copied to analyses directory if mqpar.xml is not provided (PP 01/02/19)
# 1.3.6 mqpar.xml file now managed like the other files & bug fix linked to missing position probability in evidence.txt (PP 10/01/19)
# 1.3.5 Adapted to MaxQuant versions 1.6 (PP 19/12/18)
# 1.3.4 Minor change on the getProtInfo call (VS 16/11/2018)
# 1.3.3 Minor change to revert code block for modifications allowed for quantification to v1.3.1 (PP 01/10/18)
# 1.3.2 Add extra files for maxquant import other than TMT or iTRAQ: summary.txt and parameters.txt (GA 09/08/18)
# 1.3.1 [Fix] bug peptide specificity computation now excludes decoy and optionally contaminent matches (PP 15/06/18)
# 1.3.0 Peptide quantification data now written to file $promsPath{quantification}/project_$projectID/quanti_$quantifID/peptide_quantification(_$targetPos).txt (PP 09/05/18)
# 1.2.1 Minor modification in myProMS import (GA 10/04/18)
# 1.2.0 Updated for auto-increment in tables ANALYSIS,PROTEIN,QUANTIFICATION & added compatibility with shared directory (PP 23/03/18)
# 1.1.5 Modification in <EVIDENCE> parsing because of GlyGly modification (GA 06/12/17)
# 1.1.4 Added MaxQuant version string in QUANTIF_ANNOT field (PP 13/09/17)
# 1.1.3 Matches Analysis display pos with fraction number if used (PP 09/08/17)
# 1.1.2 Bug fix in peptide start position for non-razor proteins & records identification evidence id in PEPTIDE.DATA (PP 08/03/17)
# 1.1.1 Records modification position probability & bug fix in SILAC ghost peptides & optional import of protein quantifications (PP 17/02/17)
# 1.1.0 Major revision in import mechanism and data storage design. Not compatible with older imports (PP 20/01/17)
# 1.0.6 Add UPDATE_DATE and UPDATE_USER in INSERT queries for ANALYSIS and QUANTIFICATION tables (GA 19/10/2016)
# 1.0.5 Change to deal with autoincrement in PEPTIDE table (GA 19/10/2016)
# 1.0.4 Minor modification (GA 05/10/2016)
# 1.0.3 No import for quantitative data if values are articially put to zero (GA 12/09/2016)
# 1.0.2 Minor change for MaxQuant query. MAXQUANT code becomes MQ (GA 27/07/2016)
# 1.0.1 Add labeled import (SILAC) & update design information (GA 09/05/2016)
# 1.0.0 Import for MaxQuant quantitated folder<BR> TODO: check for replicate<->fraction depedencies and update Observation according to technical replicates (GA 24/09/15)<BR>TODO Import as well Light and Heavy information
