#!/usr/local/bin/perl -w

################################################################################
# watchPhosphoAnalyses.cgi      1.0.2                                          #
# Authors: P. Poullet (Institut Curie)                                         #
# Contact: myproms@curie.fr                                                    #
# Monitors all on-going PhosphoRS analyses                                     #
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
use POSIX qw(strftime); # to get the time
use File::stat;
use File::Path qw(rmtree); # remove_tree
use promsConfig;
use promsMod;

#print header(-charset=>'UTF-8'); warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my $prsHomeDir="$promsPath{tmp}/phosphoRS";

####################
####>Parameters<####
####################
if (param('AJAX')) {
	if (param('AJAX') eq 'update') {&ajaxUpdateAnalysisStatus;}
	elsif (param('AJAX') eq 'delete') {&ajaxDeleteJob;}
	exit;
}
my $projectID=(param('projectID'))? param('projectID') : 0;

#######################
####>Starting HTML<####
#######################
my ($lightColor,$darkColor)=&promsConfig::getRowColors;

print header(-'content-encoding'=>'no',-charset=>'UTF-8'); #(-cache_control=>"no-cache, no-store, must-revalidate"); for AJAX in IE but works better in script called by AJAX function
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Monitoring PhosphoRS Analyses</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT language="JavaScript">
function doNothing() {}
// AJAX --->
function deleteJob(jobDir,anaID) {
	autoUpdate=0; // stop job status autoupdate
	//Wait for autoupdate XHR to finish if exists
	var count=0; // just to be safe
	while (XHR && XHR.readyState != 4 && count<5) {count++; setTimeout('doNothing()',1000);}
	//Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if (!XHR) {return false;}
	XHR.open("GET","./watchPhosphoAnalyses.cgi?AJAX=delete&jobDir="+jobDir+"&anaID="+anaID,true); //+"&rand="+Math.random()  (for AJAX in IE or use cache-control in child-script header)
    XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			var responseStrg=(XHR.responseText.match('##OK'))? 'Job deleted' : XHR.responseText;
			document.getElementById('status:'+jobDir+':'+anaID).innerHTML='<B>'+responseStrg+'</B>';
			//Update watchPhospho window
			if (XHR.responseText.match('##RELOAD')) { // an entire job dir was deleted => reload window
				//window.location.reload();
				setTimeout('window.location.reload()',5000);
				return;
			}
			else {
				autoUpdate=1;
				setTimeout('ajaxWatchAnalysis()',5000);
			}
		}
	}
	XHR.send(null);
}
var XHR=null;
function getXMLHTTP(){
	var xhr=null;
	if(window.XMLHttpRequest) {// Firefox & others
		xhr = new XMLHttpRequest();
	}
	else if(window.ActiveXObject){ // Internet Explorer
		try {
		  xhr = new ActiveXObject("Msxml2.XMLHTTP");
		} catch (e) {
			try {
				xhr = new ActiveXObject("Microsoft.XMLHTTP");
			} catch (e1) {
				xhr = null;
			}
		}
	}
	else { // XMLHttpRequest not supported by browser
		alert("Your browser does not support XMLHTTPRequest objects...");
	}
	return xhr;
}
// <--- AJAX
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" >
<CENTER>
<FONT class="title">Monitoring PhosphoRS Analyses</FONT>
<BR><BR><INPUT type="button" value="Refresh window" onclick="window.location.reload()"/>&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Close window" onclick="window.close()"/><BR><BR>
|;

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;
my $sthUsr=$dbh->prepare("SELECT USER_NAME FROM USER_LIST WHERE ID_USER=?");
my $sthPj=$dbh->prepare("SELECT NAME FROM PROJECT WHERE ID_PROJECT=?");

############################################
####>Scanning PhosphoRS job directories<####
############################################
my @jobDirList=&getJobList;
my (%projectJobs,%userList);
foreach my $jobDir (sort{&promsMod::sortSmart($a,$b)} @jobDirList) {
	next if (!-e "$prsHomeDir/$jobDir" || !-e "$prsHomeDir/$jobDir/prs_info.txt"); # job finished between 2 scans 
	my ($jobUserID,$jobUserName,$anaID,%prsParameters);
	my $section='';
	open (INFO,"$prsHomeDir/$jobDir/prs_info.txt");
	while (<INFO>) {
		if (/^USER=(\S+)/) {
			$jobUserID=$1;
			if ($userList{$jobUserID}) {$jobUserName=$userList{$jobUserID};}
			else {
				$sthUsr->execute($jobUserID);
				($jobUserName)=$sthUsr->fetchrow_array;
				$userList{$jobUserID}=$jobUserName;
			}
		}
		elsif (/^ANAID=(\d+)/) {$anaID=$1;}
		else {
			chomp;
			my ($param,$value)=split('=',$_);
			$prsParameters{$param}=$value;
		}
	}
	close INFO;

	next unless $anaID; # happens if job has finished
	#my $projID;
	#if ($anaID) {
	my $projID=&promsMod::getProjectID($dbh,$anaID,'analysis');

	##>Checking access rights
	my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projID);
	my $projectAccess=${$userInfo[2]}{$projID};
	next if ($projectAccess !~ /bioinfo|mass|manag/ && $userID ne $jobUserID);

	unless ($projectJobs{$projID}) {
		$sthPj->execute($projID);
		($projectJobs{$projID}{'NAME'})=$sthPj->fetchrow_array;
	}
	#}
	#else { # No matching project => PROBLEM!!!
	#	$projID=0;
	#	$projectJobs{0}{'NAME'}='**Unknown**';
	#}
	$jobUserName='Unknown' unless $jobUserName;
	#push @{$projectJobs{$projID}{'JOBS'}},[$jobDir,$jobUserName,$anaID,\%prsParameters];
	my ($masterJobCode,$jobRank)=split('\.',$jobDir);
	push @{$projectJobs{$projID}{'JOBS'}{$masterJobCode}},[$jobRank,$anaID];
	@{$projectJobs{$projID}{'PARAMS'}{$masterJobCode}}=($jobUserName,\%prsParameters) unless $projectJobs{$projID}{'PARAMS'}{$masterJobCode}; # only once per $masterJobCode
}
$sthUsr->finish;
$sthPj->finish;

####>No jobs<####
unless (scalar keys %projectJobs) {
	$dbh->disconnect;
	print qq
|<FONT class=\"title2\">No PhosphoRS Analyses found.</FONT>
<SCRIPT language="JavaScript">
setTimeout('window.location.reload()',10000);
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

####################################
####>Displaying Analyses Status<####
####################################
my $sthQN=$dbh->prepare("SELECT NAME FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
my @watchedAnalyses;
print "<TABLE>";
foreach my $projectID (sort{lc($projectJobs{$a}{'NAME'}) cmp lc($projectJobs{$b}{'NAME'})} keys %projectJobs) {
	print "<TR><TD><TABLE cellspacing=0><TR><TD class=\"title2\" colspan=3>Project '$projectJobs{$projectID}{NAME}':</TD></TR>\n";
	unless (scalar keys %{$projectJobs{$projectID}}) { # some data are missing in prs_info.txt
		print "<TR><TD class=\"title3\" colspan=3 nowrap>&nbsp;&bull;&nbsp;Not enough information on job</TD></TR></TABLE></TD></TR>\n";
		next;
	}
	foreach my $masterJobCode (sort{$a<=>$b} keys %{$projectJobs{$projectID}{'JOBS'}}) {
		my ($jobUserName,$refParam)=@{$projectJobs{$projectID}{'PARAMS'}{$masterJobCode}};
		print qq
|<TR><TH align=left colspan=3 nowrap>&nbsp;&bull;&nbsp;<FONT class="title3">Master job: $masterJobCode&nbsp;&nbsp;-&nbsp;&nbsp;Parameters:</FONT> Threshold=$refParam->{probThreshold}%, Activation type=$refParam->{activationType}, Mass tolerance=$refParam->{massTolerance} Da<FONT class="title3">&nbsp;&nbsp;-&nbsp;&nbsp;Launched by $jobUserName</FONT></TH></TR>
<TR bgcolor="$darkColor">
	<TH class="rbBorder" nowrap witdh=20>Job #</TH>
	<TH class="rbBorder" nowrap>Analysis Name</TH>
	<TH class="bBorder" nowrap width=600>Status</TH>
</TR>
|;	
		my $rowColor=$lightColor;
		foreach my $refJob (@{$projectJobs{$projectID}{'JOBS'}{$masterJobCode}}) {
			my ($jobRank,$anaID)=@{$refJob};
			my $jobDir=$masterJobCode.'.'.$jobRank;
			my ($scanAgain,$anaStatus,$anaError)=&getPhosphoAnaInfo($jobDir,$anaID);
			if ($scanAgain==3) { # jobDir was deleted => reload everything
				print qq
|</TABLE></TD></TR></TABLE>
<SCRIPT language="JavaScript">
window.location.reload();
</SCRIPT>
</BODY>
</HTML>
|;
				exit;
			}
			
			my @anaInfo=&promsMod::getItemInfo($dbh,'ANALYSIS',$anaID);
			my @anaHierarchy;
			foreach my $i (1..$#anaInfo) {
				push @anaHierarchy,$anaInfo[$i]{'NAME'};
			}
			my $anaName=join(' > ',@anaHierarchy);
	
			$anaError = "ok" unless $anaError;
			$anaError=~s/\n/<BR>/g;
			my $anaStatusStrg;
			if ($scanAgain) {
				$anaStatusStrg="<SPAN id='status:$jobDir:$anaID'>$anaStatus</SPAN>\n<DIV id=\"error:$jobDir:$anaID\" style=\"display:none\"><FIELDSET><LEGEND><B>Error message:</B></LEGEND>$anaError</FIELDSET></DIV>";
				push @watchedAnalyses,$jobDir.':'.$anaID;
			}
			else {
				$anaStatusStrg=$anaStatus;
			}
			print "<TR bgcolor=\"$rowColor\" valign=top><TH>$jobRank</TH><TH align=left valign=top nowrap>&nbsp;$anaName&nbsp;</TH><TD valign=top nowrap>&nbsp;$anaStatusStrg</TD></TR>\n";
			$rowColor=($rowColor eq $lightColor)? $darkColor : $lightColor;
		}
	}
	print "</TABLE></TD></TR>\n";
}
print "</TABLE>\n";
$sthQN->finish;

$dbh->disconnect;

if (scalar @watchedAnalyses) {
	my $watchedAnaStrg=join(',',@watchedAnalyses);
	print qq
|<SCRIPT LANGUAGE="JavaScript">
var watchedAnaStrg='$watchedAnaStrg';
var autoUpdate=1;
var noJobErrorUpdate={}; // list of jobs not to be updated because user is reading error message
function viewError(errorDivID) {
	var errorDiv = document.getElementById(errorDivID);
	if (errorDiv.style.display == 'none') {
		errorDiv.style.display = '';
		noJobErrorUpdate[errorDivID]=1;
	}
	else {
		errorDiv.style.display = 'none';
		delete noJobErrorUpdate[errorDivID];
	}
}
function updateAnalysisStatus(anaStatusTxt) {
	if (anaStatusTxt.match('##RELOAD')) { // new job(s) detected => reload window
		if (Object.keys(noJobErrorUpdate).length) { // Prevent reload to maintain error message display
			setTimeout(function(){updateAnalysisStatus(anaStatusTxt)},10000); // re-call function
		}
		else {
			//window.location.reload();
			setTimeout('window.location.reload()',10000);
			return;
		}
		return;
	}
	var statusData=anaStatusTxt.split('\\n');
	watchedAnaStrg=''; // empty list of analyses to be watched
	for (let i=0; i<statusData.length; i++) { // last line is empty
		if (!statusData[i].match('#')) continue;
		var anaData=statusData[i].split('##');
		var statusSpan=document.getElementById('status:'+anaData[1]+':'+anaData[2]);
		if (!statusSpan) {
//alert(statusData[i]+' => id="status:'+anaData[1]+':'+anaData[2]+'" not found => Reload');
			//window.location.reload();
		}
		statusSpan.innerHTML=anaData[3];
		var errorDivID='error:'+anaData[1]+':'+anaData[2];
		if (!noJobErrorUpdate[errorDivID]) {
			document.getElementById(errorDivID).innerHTML='<FIELDSET><LEGEND><B>Error message:</B></LEGEND>'+anaData[4]+'</FIELDSET>';
		}
		if (anaData[0]==1) { // unfinished ana
			if (!watchedAnaStrg.match(anaData[1])) { // job is not listed yet
				if (watchedAnaStrg) {watchedAnaStrg+=',';}
				watchedAnaStrg+=anaData[1]; // add jobDir
			}
			watchedAnaStrg+=':'+anaData[2]; // add anaID
		}
	}
	if (watchedAnaStrg && autoUpdate) {
		setTimeout('ajaxWatchAnalysis()',5000);
	}
}
function ajaxWatchAnalysis() {
	//If XHR object already exists, the request is canceled & the object is deleted
	if(XHR && XHR.readyState != 0){
		XHR.abort();
		delete XHR;
	}
	//Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if (!XHR) {return false;}
	XHR.open("GET","./watchPhosphoAnalyses.cgi?AJAX=update&anaStrg="+watchedAnaStrg,true); //+"&rand="+Math.random()  (for AJAX in IE or use cache-control in child-script header)
    XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			updateAnalysisStatus(XHR.responseText);
		}
	}
	XHR.send(null);
}
window.onload=setTimeout('ajaxWatchAnalysis()',5000); // First ajax call!
</SCRIPT>
|;
}
print "</BODY>\n</HTML>\n";


##############################<<< SUBROUTINE >>>################################
sub ajaxUpdateAnalysisStatus {
	print header(-charset=>'UTF-8',-type=>'text/plain',-cache_control=>"no-cache, no-store, must-revalidate"); # for AJAX in IE
	warningsToBrowser(1);

	###<Fetching job list
	my @jobDirList=&getJobList;

	###<Scanning job directories
	my %jobAnalyses;
	foreach my $jobAnaStrg (split(',',param('anaStrg'))) {
		my ($jobDir,$anaID)=split(':',$jobAnaStrg);
		$jobAnalyses{$jobDir}=$anaID;
	}

	###<Checking for new job(s)
	foreach my $jobDir (@jobDirList) {
		unless ($jobAnalyses{$jobDir}) { #  new job Detected!!!
			print "##RELOAD (New job $jobDir)##\n";
			exit;
		}
	}

	###<Checking on-going jobs
	my $responseStrg='';
	foreach my $jobDir (sort{&promsMod::sortSmart($a,$b)} keys %jobAnalyses) {
		my $anaID=$jobAnalyses{$jobDir};
		my ($scanAgain,$anaStatus,$anaError)=&getPhosphoAnaInfo($jobDir,$anaID);
		if ($scanAgain==3) {
			$responseStrg="##RELOAD (job recorded $jobDir)##\n";
			last;
		}
		$anaError = '' unless $anaError;
		$anaError=~s/\n/<BR>/g;
		$responseStrg.="$scanAgain##$jobDir##$anaID##$anaStatus##$anaError\n";
	}
	print $responseStrg;
	exit;
}

sub getJobList {
	my @jobDirList;
	opendir (DIR, $prsHomeDir);
	while (defined (my $jobDir = readdir (DIR))) {
		if ($jobDir =~ /^(\d{14})/ && -d "$prsHomeDir/$jobDir") { # directory with numbers only 20130808100527(.123)
			my $jobNum=$1;
			unless (-e "$prsHomeDir/$jobDir/prs_info.txt") {
				#my $now=strftime("%Y%m%d%H%M%S",localtime);
				rmtree("$prsHomeDir/$jobDir") if -e "$prsHomeDir/$jobDir"; # if $now-$jobNum > 604800; # 1 week
				rmdir "$prsHomeDir/$jobDir" if -e "$prsHomeDir/$jobDir";
				next;
			}
			push @jobDirList,$jobDir;
		}
	}
	close DIR;
	return @jobDirList;
}

sub getPhosphoAnaInfo {
	my ($jobDir,$anaID)=@_;
	my ($anaStatus,$scanAgain,$anaError)=('',1,'');
	my $delButtonStrg="&nbsp;<INPUT type=\"button\" value=\"Delete job\" onclick=\"deleteJob('$jobDir',$anaID)\">";

	unless (-e "$prsHomeDir/$jobDir") {return(3,0,'');} # dir has just been deleted!
	if (-e "$prsHomeDir/current/$anaID\_$jobDir\_wait.flag") {
		$anaStatus='<B>Pending...</B> [Waiting for job to start]';
		$scanAgain=1; # tagged as unfinished so always scanned
	}
	elsif (-e "$prsHomeDir/$jobDir/status.out") {
		open(OUT, "$prsHomeDir/$jobDir/status.out");
		my ($firstLine,$lastLine)=('','');
		while (<OUT>) {
			chomp;
			if ($.==1) {$firstLine=$_;}
			else {$lastLine=$_;}
		}
		close OUT;
		if ($lastLine=~/^Ended /) { # ended
			$anaStatus="<B>Finished!</B> [$firstLine -> $lastLine]";
			$scanAgain=0;
		}
		else { # started
			$anaStatus=($firstLine)? "<B>$firstLine</B>" : '<B>Running...</B>'; # in case pb with $|
			my $progressStrg=" [$lastLine" if ($lastLine && $lastLine ne $firstLine);
			$scanAgain=1;
			my $refFile=(-e "$prsHomeDir/$jobDir/output_data.xml")? "$prsHomeDir/$jobDir/output_data.xml"
					  : (-e "$prsHomeDir/$jobDir/input_data.xml")? "$prsHomeDir/$jobDir/input_data.xml"
					  : "$prsHomeDir/$jobDir/status.out";
			open(my $fh, "<", $refFile);
			my $lastChangeInSec=stat($fh)->[9]; # in sec
			close $fh;
			my $now=strftime("%s",localtime); # in sec
			my $waitTime=strftime("%Hh %Mm %Ss",localtime($now-$lastChangeInSec-3600));
			$anaStatus.="$progressStrg (Updated $waitTime ago)]&nbsp;";
			my $maxWaitTime=14400; # 4 h in sec
			$anaStatus.="<BR><B><FONT color='#DD0000'>Wait &gt ".($maxWaitTime/3600)." h! Possible job failure. </FONT>$delButtonStrg</B>" if $now-$lastChangeInSec > $maxWaitTime; # in sec
		}

		my $error=0;
		$error= -s "$prsHomeDir/current/$anaID\_$jobDir\_error.txt" if -e "$prsHomeDir/current/$anaID\_$jobDir\_error.txt";
		chomp $error if $error;
		if ($error) {
			my $errorString;
			open(ERRORFILE, "$prsHomeDir/current/$anaID\_$jobDir\_error.txt") or die $!;
			while (<ERRORFILE>) {
				$errorString .= $_;
			}
			close ERRORFILE;
			$anaStatus.="<B><FONT color=\"#DD0000\">**ERROR**</FONT>";
			$anaStatus.=$delButtonStrg unless $anaStatus=~/deleteJob/; # delete button already added
			$anaStatus.=qq|&nbsp;<INPUT type="button" value="Show/Hide error" onclick="viewError('error:$jobDir:$anaID');"></B>|;
			$anaError = $errorString;
			$scanAgain=1;
		}
	}
	return ($scanAgain,$anaStatus,$anaError);
}

sub ajaxDeleteJob {
	print header(-charset=>'UTF-8',-type=>'text/plain',-cache_control=>"no-cache, no-store, must-revalidate"); # for AJAX in IE
	warningsToBrowser(1);

	my $jobDir=param('jobDir');
	my $anaID=param('anaID');

	unlink glob "$prsHomeDir/current/*_$jobDir\_*" if glob "$prsHomeDir/current/*_$jobDir\_*";

	##<Delete job directory
	if (-e "$prsHomeDir/$jobDir") {
		rmtree("$prsHomeDir/$jobDir");
		rmdir "$prsHomeDir/$jobDir" if -e "$prsHomeDir/$jobDir";
		print "##RELOAD\n##OK";
	}
	else {print 'Job directory not found.';}
	exit;
}

####>Revision history<####
# 1.0.2 [Fix] JS bug in deleteJob function (PP 07/01/19)
# 1.0.1 [Fix] Bug undef $anaID when job finishes during prs_info.txt file scan (PP 05/12/18)
# 1.0.0 Forked from watchQuantifications.cgi 1.3.6 (PP 08/11/18)
