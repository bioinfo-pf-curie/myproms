#!/usr/local/bin/perl -w

################################################################################
# watchQuantifications.cgi       1.5.0                                         #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Monitors all on-going quantifications                                        #
# called from processAnalyses.cgi                                              #
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
use promsQuantif;
#exit;

#print header(-charset=>'UTF-8'); warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my %quantifProcesses=('XIC'=>'Ext. ion chrom.','EMPAI'=>'emPAI','SIN'=>'SI<SUB>N</SUB>','XICCORR'=>'XIC correction',
					  'IMPORT_MQ'=>'MaxQuant Import',
					  'DESIGN'=>'Protein ratio','PROT_RATIO_PEP'=>'Protein ratio',
					  'DESIGN:SimpleRatio'=>'Protein ratio','DESIGN:SuperRatio'=>'Super ratio',
					  'DESIGN:myProMS'=>'Protein ratio',
					  'DESIGN:PEP_INTENSITY'=>'Protein ratio (Peptide intensity)',
					  'DESIGN:PEP_RATIO'=>'Protein ratio (Peptide ratio)',
					  'DESIGN:MSstats'=>'DIA-based Protein ratio (MSstats)',
					  'DESIGN:SSPA'=>'SSP Analysis'
					  ); #'SILAC'=>'SILAC','ITRAQ'=>'iTRAQ', 'DESIGN:SWATH/MSstats'=>'SWATH-based Protein ratio (MSstats)',
my $quantifHomeDir="$promsPath{tmp}/quantification";

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
<TITLE>Monitoring Quantifications</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT type="text/javascript">
function doNothing() {}
// AJAX --->
function deleteJob(jobDir,quantItem) {
	autoUpdate=0; // stop quantif status autoupdate
	//Wait for autoupdate XHR to finish if exists
	var count=0; // just to be safe
	while (XHR && XHR.readyState != 4 && count<5) {count++; setTimeout('doNothing()',1000);}
	//Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if (!XHR) {return false;}
	XHR.open("GET","./watchQuantifications.cgi?AJAX=delete&jobDir="+jobDir+"&quantItem="+quantItem,true); //+"&rand="+Math.random()  (for AJAX in IE or use cache-control in child-script header)
    XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			//Refresh Quantifications frame if displayed
			if (opener) {
				var selView=opener.parent.itemFrame.document.getElementById('view');
				if (selView && selView.value=='quanti') opener.parent.itemFrame.location.reload();
			}
			//Update watchQuanti window
			if (XHR.responseText.match('##DELETE_EXP')) {
				alert('Experiment data were not deleted.\\nUse "Edit" for Experiment then "Delete Experiment" to force-delete all data');
				window.location.reload();
			}
			if (XHR.responseText.match('##RELOAD')) { // an entire job dir was deleted => reload window
				window.location.reload();
				//setTimeout('window.location.reload()',5000);
				//return;
			}
			else {
				if (document.getElementById('status:'+jobDir+':'+quantItem)) { // can be undef after 1st delete
					//var responseStrg=(XHR.responseText.match('##OK'))? 'Job deleted' : XHR.responseText;
					//document.getElementById('status:'+jobDir+':'+quantItem).innerHTML='<B>'+responseStrg+'</B>';
					document.getElementById('status:'+jobDir+':'+quantItem).innerHTML='<B>'+XHR.responseText+'</B>';
				}
				else {alert('Deletion of job '+jobDir+': '+XHR.responseText);}
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
<FONT class="title">Monitoring Quantification Jobs</FONT>
<BR><BR><INPUT type="button" value="Refresh window" onclick="window.location.reload()"/>&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Close window" onclick="window.close()"/><BR><BR>
|;

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;
my $sthUsr=$dbh->prepare("SELECT USER_NAME FROM USER_LIST WHERE ID_USER=?");
my $sthPj=$dbh->prepare("SELECT NAME FROM PROJECT WHERE ID_PROJECT=?");

#############################################
####>Scanning Quantification directories<####
#############################################
my @jobDirList=&getJobList;
my (%projectJobs,%userList);
foreach my $jobDir (sort{&promsMod::sortSmart($a,$b)} @jobDirList) {
	my ($jobUserID,$jobUserName,$quantifType,$designID,$projID,@quantItemList,@quantificationList,@quantifNameList,@designList,@experimentList);
	my $section='';
	open (INFO,"$quantifHomeDir/$jobDir/quantif_info.txt");
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
		elsif (/^TYPE=(\S+)/) {$quantifType=$1 unless $quantifType;} #$quantifProcesses{$1}
		#elsif (/^ALGO_TYPE\s+(\S+)/) {$quantifType='DESIGN:'.$1;} # Modified on 2017/05/04
		elsif (/^ALGO_TYPE\s+(\S+)/) {
			$quantifType.=($quantifType)? ':'.$1 : $1; # for internal quantifications, no design
		}
		elsif (/^PARAMETERS:/) {$section='parameters'; next;}
		elsif (/^ANALYSES/) {$section='analyses'; next;}
		elsif (/^QUANTIFICATIONS:/) {$section='quantifications'; next;}
		elsif ($section eq 'analyses' && /^(\S+)/) {
			push @quantItemList,$1;
		}
		elsif ($section eq 'parameters') {
			if (/^QUANTIF_NAME\t\t(.+)\n/) {push @quantifNameList,$1;} # in case not quantifID yet
			elsif (/^ID_DESIGN\s+(\d+)/) {$designID=$1; push @designList,-$designID;}
		}
		elsif ($section eq 'quantifications' && /^(\d+)/) {push @quantificationList,$1;}
		elsif ($quantifType eq 'IMPORT_MQ') {
			if (/^ID_EXPERIMENT=(\d+)/) {
				my $expID=$1;
				my ($expName)=$dbh->selectrow_array("SELECT NAME FROM EXPERIMENT WHERE ID_EXPERIMENT=$expID");
				push @experimentList,[$expID,$expName];
			}
			elsif (/^ID_PROJECT=(\d+)/) {$projID=$1;}
		}
	}
	close INFO;

	if ($quantItemList[0]) {
		my ($anaID,$parentQuantif)=split(/\./,$quantItemList[0]); # anaID or anaID.parentQuantifID
		$projID=&promsMod::getProjectID($dbh,$anaID,'analysis');
	}
	elsif ($quantificationList[0]) {
		$projID=&promsMod::getProjectID($dbh,$quantificationList[0],'quantification');
	}
	elsif ($designID) {
		$projID=&promsMod::getProjectID($dbh,$designID,'design');
	}
	if ($projID) {
		##>Checking access rights
		my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projID);
		my $projectAccess=${$userInfo[2]}{$projID};
		next if ($projectAccess !~ /bioinfo|mass|manag/ && $userID ne $jobUserID);

		unless ($projectJobs{$projID}) {
			$sthPj->execute($projID);
			($projectJobs{$projID}{'NAME'})=$sthPj->fetchrow_array;
		}
	}
	else { # No matching project => PROBLEM!!!
		$projID=0;
		$projectJobs{0}{'NAME'}='**Unknown**';
	}
	#if ($quantItemList[0]) {push @{$projectJobs{$projID}{'JOBS'}},[$jobDir,$quantifType,$jobUserName,\@quantItemList];}
	#if ($quantificationList[0]) {push @{$projectJobs{$projID}{'JOBS'}},[$jobDir,$quantifType,$jobUserName,\@quantificationList];}
	$quantifProcesses{$quantifType}="Unknown ($quantifType)" unless $quantifProcesses{$quantifType};
	$jobUserName='Unknown' unless $jobUserName;
	my @jobData=($jobDir,$quantifType,$jobUserName);
	if ($quantItemList[0]) {push @jobData,\@quantItemList;}
	elsif ($quantificationList[0]) {push @jobData,\@quantificationList;}
	elsif ($quantifNameList[0]) {push @jobData,(\@designList,\@quantifNameList);} # for design-based multi-launch with no quantifID yet
	elsif ($experimentList[0]) {push @jobData,\@experimentList;} # MaxQuant import
	push @{$projectJobs{$projID}{'JOBS'}},\@jobData;
}
$sthUsr->finish;
$sthPj->finish;

####>No jobs<####
unless (scalar keys %projectJobs) {
	$dbh->disconnect;
	print qq
|<FONT class=\"title2\">No quantifications found.</FONT>
<SCRIPT type="text/javascript">
setTimeout('window.location.reload()',30000);
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

##########################################
####>Displaying Quantification Status<####
##########################################
my $sthQN=$dbh->prepare("SELECT NAME FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
my %watchedAnalyses;
print "<TABLE>";
foreach my $projectID (sort{lc($projectJobs{$a}{'NAME'}) cmp lc($projectJobs{$b}{'NAME'})} keys %projectJobs) {
	print "<TR><TD><TABLE cellspacing=0><TR><TD class=\"title2\" colspan=3>Project $projectJobs{$projectID}{NAME}:</TD></TR>\n";
	unless ($projectJobs{$projectID}{'JOBS'}) { # some data are missing in quantif_info.txt (probably not yet added by launchQuantifications.pl)
		print "<TR><TD class=\"title3\" colspan=3 nowrap>&nbsp;&bull;&nbsp;Not enough information on job</TD></TR></TABLE></TD></TR>\n";
		next;
	}
	my %designHierarchyStrg;
	foreach my $refJob (@{$projectJobs{$projectID}{'JOBS'}}) {
		my $jobDir=$refJob->[0];
		print "<TR><TD class=\"title3\" colspan=3 nowrap>&nbsp;&bull;&nbsp;Job: $jobDir&nbsp;&nbsp;-&nbsp;&nbsp;Quantification method: ",$quantifProcesses{$refJob->[1]},"&nbsp;&nbsp;-&nbsp;&nbsp;Launched by user: ",$refJob->[2];
		if ($projectID==0) {
			print "&nbsp;<INPUT type=\"button\" value=\"Delete job\" onclick=\"deleteJob($jobDir,'')\">&nbsp;</TD></TR>\n";
			next;
		}
		else {print "&nbsp;</TD></TR>\n";}
		
		if ($refJob->[1] eq 'IMPORT_MQ') {
			print qq
|<TR bgcolor="$darkColor">
	<TH class="rbBorder" nowrap width=20>#</TH>
	<TH class="rbBorder" nowrap>Experiment Name</TH>
	<TH class="bBorder" nowrap width=600>Status</TH>
</TR>
|;
			my $rowColor=$lightColor;
			my $countImport=0;
			foreach my $refExp (@{$refJob->[3]}) {
				my ($expID,$expName)=@{$refExp};
				my ($scanAgain,$importStatus,$importError,$detailsFile)=&getAnaQuantifInfo($jobDir,$expID);
				if ($scanAgain==3) { # jobDir was deleted => reload everything
					print qq
|</TABLE></TD></TR></TABLE>
<SCRIPT type="text/javascript">
window.location.reload();
</SCRIPT>
</BODY>
</HTML>
|;
					exit;
				}
				$countImport++;

				$importError = "ok" unless $importError;
				$importError=~s/\n/<BR>/g;
				my $importStatusStrg;
				if ($scanAgain) {
					$importStatusStrg="<SPAN id='status:$jobDir:$expID'>$importStatus</SPAN>\n<DIV id=\"error:$jobDir:$expID\" style=\"display:none\"><FIELDSET><LEGEND><B>Error message:</B></LEGEND>$importError</FIELDSET></DIV>";
					push @{$watchedAnalyses{$jobDir}},$expID; # quantItems to be watched
				}
				else {
					$importStatusStrg=$importStatus;
				}
				print "<TR bgcolor=\"$rowColor\" valign=top><TH align=right valign=top>$countImport&nbsp;</TH><TH align=left valign=top nowrap>&nbsp;$expName&nbsp;</TH><TD valign=top nowrap>&nbsp;$importStatusStrg</TD></TR>\n";

				$rowColor=($rowColor eq $lightColor)? $darkColor : $lightColor;
			}
		}
		elsif (uc($refJob->[1]) !~ /DESIGN/ && uc($refJob->[1]) ne 'XIC') { # MassChroQ=XIC
			print qq
|<TR bgcolor="$darkColor">
	<TH class="rbBorder" nowrap width=20>#</TH>
	<TH class="rbBorder" nowrap>Analysis Name</TH>
	<TH class="bBorder" nowrap width=600>Status</TH>
</TR>
|;
			#>Analyses quantif
			my $rowColor=$lightColor;
			my $countItem=0;
			foreach my $quantItemID (@{$refJob->[3]}) {
				my ($scanAgain,$anaStatus,$anaError)=&getAnaQuantifInfo($jobDir,$quantItemID);
				if ($scanAgain==3) { # jobDir was deleted => reload everything
					print qq
|</TABLE></TD></TR></TABLE>
<SCRIPT type="text/javascript">
window.location.reload();
</SCRIPT>
</BODY>
</HTML>
|;
					exit;
				}
				$countItem++;
				my ($anaID,$parentQuantifID)=split(/\./,$quantItemID); # anaID or anaID.parentQuantifID
				my @anaInfo=&promsMod::getItemInfo($dbh,'ANALYSIS',$anaID);
				my @anaHierarchy;
				foreach my $i (1..$#anaInfo) {
					push @anaHierarchy,$anaInfo[$i]{'NAME'};
				}
				my $anaName=join(' > ',@anaHierarchy);

				#>Parent Quantification
				if ($parentQuantifID && uc($refJob->[1]) eq 'XICCORR') { # =~ /SILAC|ITRAQ|TMT/
					$sthQN->execute($parentQuantifID);
					my ($parQuantifName)=$sthQN->fetchrow_array;
					$anaName.="&nbsp;<BR>&nbsp;&nbsp;&nbsp;&nbsp;using Quantification: $parQuantifName";
				}

				$anaError = "ok" unless $anaError;
				$anaError=~s/\n/<BR>/g;
				my $anaStatusStrg;
				if ($scanAgain) {
					$anaStatusStrg="<SPAN id='status:$jobDir:$quantItemID'>$anaStatus</SPAN>\n<DIV id=\"error:$jobDir:$quantItemID\" style=\"display:none\"><FIELDSET><LEGEND><B>Error message:</B></LEGEND>$anaError</FIELDSET></DIV>";
					push @{$watchedAnalyses{$jobDir}},$quantItemID; # quantItems to be watched
				}
				else {
					$anaStatusStrg=$anaStatus;
				}
				print "<TR bgcolor=\"$rowColor\" valign=top><TH align=right valign=top>$countItem&nbsp;</TH><TH align=left valign=top nowrap>&nbsp;$anaName&nbsp;</TH><TD valign=top nowrap>&nbsp;$anaStatusStrg</TD></TR>\n";

				$rowColor=($rowColor eq $lightColor)? $darkColor : $lightColor;
			}
			#print "</TABLE></TD></TR>\n";
		}
		else {
			print qq
|<TR bgcolor="$darkColor">
	<TH class="rbBorder" nowrap width=10>#</TH>
	<TH class="rbBorder" nowrap width=350>Quantification Name</TH>
	<TH class="bBorder" nowrap width=600>Status</TH>
</TR>
|;
			my $rowColor=$lightColor;
			my $countQuant=0;
			foreach my $refenceID (@{$refJob->[3]}) {
				my ($quantID,$quantName);
				if ($refenceID > 0) { # normal case
					$quantID=$refenceID;
					my @quantInfo=&promsMod::getItemInfo($dbh,'QUANTIFICATION',$quantID);
					my @quantHierarchy;
					foreach my $i (1..$#quantInfo) {
						push @quantHierarchy,$quantInfo[$i]{'NAME'};
					}
					$quantName=join(' > ',@quantHierarchy);
				}
				else { # < 0 => -designID: quantif not yet created in DB (multi-job launch))
					$quantID=0;
					my $designID=abs($refenceID);
					unless ($designHierarchyStrg{$designID}) {
						my @quantInfo=&promsMod::getItemInfo($dbh,'DESIGN',$designID);
						my @quantHierarchy;
						foreach my $i (1..$#quantInfo) {
							push @quantHierarchy,$quantInfo[$i]{'NAME'};
						}
						$designHierarchyStrg{$designID}=join(' > ',@quantHierarchy);
					}
					$quantName=$designHierarchyStrg{$designID}.' > '.$refJob->[4][$countQuant-1]; # quantifName always [0] in fact
				}
				my ($scanAgain,$quantStatus,$anaError,$RoutFile)=&getAnaQuantifInfo($jobDir,$quantID);
				if ($scanAgain==3) { # jobDir was deleted => reload everything
					print qq
|</TABLE></TD></TR></TABLE>
<SCRIPT type="text/javascript">
window.location.reload();
</SCRIPT>
</BODY>
</HTML>
|;
					exit;
				}
				$countQuant++;
				$anaError = "ok" unless $anaError;
				$anaError=~s/\n/<BR>/g;
				my $quantStatusStrg;
				if ($scanAgain) {
					if ($quantID) {
						$quantStatusStrg="<SPAN id='status:$jobDir:$quantID'>$quantStatus</SPAN><DIV id=\"error:$jobDir:$quantID\" style=\"display:none\"><FIELDSET><LEGEND><B>Error message:</B></LEGEND>$anaError</FIELDSET></DIV>";
						#if ($RoutFile && -e "$quantifHomeDir/$jobDir/quanti_$quantID/$RoutFile") {
						#	open (ROUT,"$quantifHomeDir/$jobDir/quanti_$quantID/$RoutFile");
						#	$quantStatusStrg.="\n<DIV id=\"rout:$jobDir:$quantID\" style=\"display:none\"><FIELDSET><LEGEND><B>R execution out file:</B></LEGEND><CODE>";
						#	while (my $line=<ROUT>) {
						#		$quantStatusStrg.="$line<BR>";
						#	}
						#	close INFO;
						#	$quantStatusStrg.="\n</CODE></FIELDSET></DIV>";
						#}
					}
					else {
						$quantStatusStrg=$quantStatus;
					}
					push @{$watchedAnalyses{$jobDir}},$quantID; # analyses to be watched ($quantID=0 for not yet created)
				}
				else {
					$quantStatusStrg=$quantStatus;
				}
				print "<TR bgcolor=\"$rowColor\"><TH align=right valign=top>$countQuant&nbsp;</TH><TH align=left valign=top nowrap>&nbsp;$quantName</TH><TD valign=top nowrap>&nbsp;$quantStatusStrg</TD></TR>\n";

				$rowColor=($rowColor eq $lightColor)? $darkColor : $lightColor;
			}

		}
	}
	print "</TABLE></TD></TR>\n";
}
print "</TABLE>\n";
$sthQN->finish;

$dbh->disconnect;

if (scalar keys %watchedAnalyses) {
	my $watchedAnaStrg='';
	foreach my $jobDir (sort{&promsMod::sortSmart($a,$b)} keys %watchedAnalyses) {
		$watchedAnaStrg.=',' if $watchedAnaStrg;
		$watchedAnaStrg.="$jobDir:".join(':',@{$watchedAnalyses{$jobDir}});
	}
	print qq
|<SCRIPT type="text/javascript">
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
			setTimeout('window.location.reload()',30000);
			return;
		}
		return;
	}
	var statusData=anaStatusTxt.split('\\n');
	watchedAnaStrg=''; // empty list of analyses to be watched
	for (var i=0; i<statusData.length; i++) { // last line is empty
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
			watchedAnaStrg+=':'+anaData[2]; // add anaID or anaID.parentQuantifID
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
	XHR.open("GET","./watchQuantifications.cgi?AJAX=update&anaStrg="+watchedAnaStrg,true); //+"&rand="+Math.random()  (for AJAX in IE or use cache-control in child-script header)
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

	###<Scanning quantification directories
	my %jobAnalyses;
	foreach my $jobAnaStrg (split(',',param('anaStrg'))) {
		my ($jobDir,@quantItemList)=split(':',$jobAnaStrg);
		@{$jobAnalyses{$jobDir}}=@quantItemList;
	}

	###<Checking for new job(s)
	foreach my $jobDir (@jobDirList) {
		unless ($jobAnalyses{$jobDir}) { #  new job Detected!!!
			print "##RELOAD (New job $jobDir)##\n";
			exit;
		}
	}

	###<Checking on-going quantifications
	my $responseStrg='';
	JOB:foreach my $jobDir (sort{&promsMod::sortSmart($a,$b)} keys %jobAnalyses) {
		foreach my $quantItemID (@{$jobAnalyses{$jobDir}}) {
			my ($scanAgain,$anaStatus,$anaError)=&getAnaQuantifInfo($jobDir,$quantItemID);
			if ($scanAgain==3) {
				$responseStrg="##RELOAD (quantif recorded $jobDir)##\n";
				last JOB;
			}
			elsif ($scanAgain==2) {next JOB;} # Do not add to responseStrg (no matching SPAN to update because fix text '<B>Pending...</B> [Request to be processed]')
			$anaError = '' unless $anaError;
			$anaError=~s/\n/<BR>/g;
			$responseStrg.="$scanAgain##$jobDir##$quantItemID##$anaStatus##$anaError\n";
		}
	}
	print $responseStrg;
	exit;
}

sub getJobList {
	my @jobDirList;
	opendir (DIR, $quantifHomeDir);
	while (defined (my $jobDir = readdir (DIR))) {
		if ($jobDir =~ /^(\d{14})/ && -d "$quantifHomeDir/$jobDir") { # directory with numbers only 20130808100527(.123)
			my $jobNum=$1;
			unless (-e "$quantifHomeDir/$jobDir/quantif_info.txt") {
				my $now=strftime("%Y%m%d%H%M%S",localtime);
				#remove_tree("$quantifHomeDir/$jobDir") if $now-$jobNum > 604800; # 1 week
				rmtree("$quantifHomeDir/$jobDir") if $now-$jobNum > 604800; # 1 week
				next;
			}
			push @jobDirList,$jobDir;
		}
	}
	close DIR;
	return @jobDirList;
}

sub getAnaQuantifInfo {
	my ($jobDir,$quantItemID)=@_;
	$quantItemID=0 unless $quantItemID; # for job with quantif not yet recorded in DB (multi-job launch)
	my ($anaStatus,$scanAgain,$anaError)=('',0,'');
	my $delButtonStrg="&nbsp;<INPUT type=\"button\" value=\"Delete job\" onclick=\"deleteJob($jobDir,'$quantItemID')\">";
	my $RoutFile=''; # needed for displaying in DIV
	my $OKquantifID=0;
	unless (-e "$quantifHomeDir/$jobDir") {return(3,0,'','');} # dir has just been deleted!
	if (-e "$quantifHomeDir/$jobDir/quantif_info.txt") {
		my ($quantifType,$ratioType,$algoType);
		open (INFO,"$quantifHomeDir/$jobDir/quantif_info.txt");
		while (<INFO>) {
			if (/^TYPE=(\S+)/) {$quantifType=$1;}
			elsif (/^RATIO_TYPE\t.*\t(\S+)/) {$ratioType=$1;}
			elsif (/^ALGO_TYPE\t.*\t(\S+)/) {$algoType=$1;}
			elsif (/^QUANTIFICATIONS:/) {$OKquantifID=1;}
		}
		close INFO;
		$ratioType='' unless $ratioType;
		$algoType='myProMS' if $ratioType eq 'Ratio';
		#$RoutFile=($ratioType eq 'Ratio')? 'LabeledAndUnlabeledQuanti.Function.Rout' : ($ratioType=~/^S*Ratio/)? 'AnalysisDiffLimma.Rout' : '';
		if ($quantItemID && $quantifType =~ /PROT_PEP_RATIO|DESIGN/ && -d "$quantifHomeDir/$jobDir/quanti_$quantItemID") {
			my $RoutFileInfo=`ls $quantifHomeDir/$jobDir/quanti_$quantItemID | grep Rout\$`; # grep better than `ls .../*.Rout` because no error if no Rout file
			chomp $RoutFileInfo;
			$RoutFile=(split(/\s/,$RoutFileInfo))[-1] if $RoutFileInfo;
		}
		#my $quantiTime=($quantifType eq 'XIC')? 300 : ($quantifType eq 'SIN')? 40 : 5; # min
		my $quantiTime=($quantifType eq 'XICCORR')? 300 : ($quantifType eq 'XIC')? 600 : ($quantifType eq 'SIN')? 1800 : ($algoType && $algoType=~/MSstats/)? 14400 : ($quantifType eq 'IMPORT_MQ')? 900 : 3600; # sec

		my $error=0;
		my $errorDisplayStrg='';
		if ($quantItemID && -e "$quantifHomeDir/current/$quantItemID\_$jobDir\_error.txt") {
			$error= -s "$quantifHomeDir/current/$quantItemID\_$jobDir\_error.txt" if -e "$quantifHomeDir/current/$quantItemID\_$jobDir\_error.txt";
			chomp $error if $error;
			if ($error) {
				my $errorString;
				open(ERRORFILE, "$quantifHomeDir/current/$quantItemID\_$jobDir\_error.txt") or die $!;
				while (<ERRORFILE>) {
					$errorString .= $_;
					if ($. > 25) {
						$errorString.='[Truncated]...';
						last;
					}
				}
				close ERRORFILE;
				$errorDisplayStrg="<B><FONT color=\"#DD0000\">**ERROR**</FONT>";
				$errorDisplayStrg.=$delButtonStrg; #### unless $anaStatus=~/deleteJob/; # delete button already added
				$errorDisplayStrg.=qq|&nbsp;<INPUT type="button" value="Show/Hide error" onclick="viewError('error:$jobDir:$quantItemID');"></B>|; # id="errorButton:$jobDir:$quantItemID"
				#$anaStatus.=qq|&nbsp;<INPUT type="button" value="Show/Hide Rout" onclick="viewError('rout:$jobDir:$quantItemID');"></B>| if -e "$quantifHomeDir/$jobDir/quanti_$quantItemID/$RoutFile"; # id="routButton:$jobDir:$quantItemID"
				$anaError = $errorString;
				$scanAgain=1;
			}
		}
		
		if ($quantItemID==0 || -e "$quantifHomeDir/current/$jobDir\_request.flag") { # quantif not yet recorded in DB: design-based multi job launch
			$anaStatus='<B>Pending...</B> [Request to be processed]';
			$scanAgain=($OKquantifID)? 3 : 2; # 2: no change, 3: quantification was created in DB => RELOAD
		}
		elsif (-e "$quantifHomeDir/current/$quantItemID\_$jobDir\_wait.flag") {
			$anaStatus='<B>Pending...</B> [Waiting for job to start]';
			$anaStatus.=$errorDisplayStrg if $error;
			$scanAgain=1; # tagged as unfinished so always scanned
		}
		else {
			if (-e "$quantifHomeDir/$jobDir/status_$quantItemID.out") {
				open(OUT, "$quantifHomeDir/$jobDir/status_$quantItemID.out");
				my ($firstLine,$lastLine)=('','');
				while (<OUT>) {
					chomp;
					if ($.==1) {$firstLine=$_;}
					else {$lastLine=$_;}
				}
				#my $startInSec=stat($out)->[9]; # File::stat
				close OUT;
				#my $firstLine=`head -1 $quantifHomeDir/$jobDir/status_$quantItemID.out`;
				#chomp $firstLine;
				#my $lastLine=`tail -1 $quantifHomeDir/$jobDir/status_$quantItemID.out`;
				#chomp $lastLine;

				if ($lastLine=~/^Ended /) { # ended
					$anaStatus="<B>Finished!</B> [$firstLine -> $lastLine]";
					$scanAgain=0;
				}
				else { # started
					$anaStatus=($firstLine)? "<B>$firstLine</B>" : '<B>Running...</B>'; # in case pb with $|
					my $progressStrg=" [$lastLine" if ($lastLine && $lastLine ne $firstLine);
					$scanAgain=1;
					#<Check for run time exceeded
					#if ($firstLine) {
					#	my ($h0,$m0,$s0,$d0,$M0,$y0)=($firstLine=~/Started (\d+):(\d+):(\d+) (\d+)\/(\d+)\/(\d+)/);
					#	my $start=((((($M0*31)+$d0)*24)+$h0)*60)+$m0; # approx. minutes in current year
					#	my $date=strftime("%H:%M:%S %d/%m/%Y",localtime);
					#	my ($h1,$m1,$s1,$d1,$M1,$y1)=($date=~/(\d+):(\d+):(\d+) (\d+)\/(\d+)\/(\d+)/);
					#	my $now=((((($M1*31)+$d1)*24)+$h1)*60)+$m1; # approx. minutes in current year
					#	$anaStatus="<B><FONT color='#DD0000'>**RUN TIME EXCEEDED**</FONT>$delButtonStrg</B>" if $now-$start > $quantiTime; # minutes
					#	#$scanAgain=0;
					#}
					#my $now=strftime("%s",localtime); # in sec
					#my $totalDuration=strftime("%H:%M:%S",localtime($now-$startInSec-3600));
					#$anaStatus.=" $totalDuration]";
					my $refFile;
					if ($quantifType eq 'XIC' && $lastLine =~ /^2\/3 /) { # very long run not yet ended
						$refFile=(-e "$quantifHomeDir/$jobDir/masschroq_status2.txt")? "$quantifHomeDir/$jobDir/masschroq_status2.txt" : "$quantifHomeDir/$jobDir/masschroq_status1.txt";
						if ($refFile=~/masschroq_status2/) {
							open(MCQ, $refFile);
							my $status='';
							my $sampNum=0;
							while (<MCQ>) {
								if (/^Parsing XML/) {
									$status='parsing';
									$progressStrg=' [2.1/3 Parsing XML input file...';
									$sampNum=0;
								}
								elsif (/^Alignment method/) {
									$status='aligning';
									$progressStrg=' [2.2/3 Aligning MS run';
									$sampNum=0;
								}
								elsif (/^Quantification method/) {
									$status='quantification';
									$progressStrg=' [2.3/3 Quantifying in MS run';
									$sampNum=0;
								}
								elsif (/ quantification finished/) {
									$status='post matching';
									$progressStrg=' [2.4/3 Post matching in MS run';
									$sampNum=0;
								}
								elsif (($status eq 'parsing' && /MS run 'samp\d+'/) || / MS run 'samp\d+'/) {
									$sampNum++;
								}
							}
							close MCQ;
							$progressStrg.=" #$sampNum" if $sampNum;
						}
					}
					elsif ($quantifType =~ /PROT_PEP_RATIO|DESIGN/) {
						$refFile=($RoutFile && -e "$quantifHomeDir/$jobDir/quanti_$quantItemID/$RoutFile")? "$quantifHomeDir/$jobDir/quanti_$quantItemID/$RoutFile" : "$quantifHomeDir/$jobDir/quanti_$quantItemID/data/table.txt";
					}
					else {
						$refFile="$quantifHomeDir/$jobDir/status_$quantItemID.out";
					}
					unless ($error) {
						$refFile="$quantifHomeDir/$jobDir/status_$quantItemID.out" unless -e $refFile;
						open(my $fh, "<", $refFile);
						my $lastChangeInSec=stat($fh)->[9]; # in sec
						close $fh;
						my $now=strftime("%s",localtime); # in sec
						my $waitTime=strftime("%Hh %Mm %Ss",localtime($now-$lastChangeInSec-3600));
						$anaStatus.="$progressStrg (Updated $waitTime ago)]&nbsp;";
						$anaStatus.="<BR><B><FONT color='#DD0000'>Wait &gt ".($quantiTime/60)." min! Possible job failure. </FONT>$delButtonStrg</B>" if $now-$lastChangeInSec > $quantiTime; # in sec
					}
				}
			}
			
			#my $error=0;
			#$error= -s "$quantifHomeDir/current/$quantItemID\_$jobDir\_error.txt" if -e "$quantifHomeDir/current/$quantItemID\_$jobDir\_error.txt";
			#chomp $error if $error;
			#if ($error) {
			#	my $errorString;
			#	open(ERRORFILE, "$quantifHomeDir/current/$quantItemID\_$jobDir\_error.txt") or die $!;
			#	while (<ERRORFILE>) {
			#		$errorString .= $_;
			#	}
			#	close ERRORFILE;
			#	$anaStatus.="<B><FONT color=\"#DD0000\">**ERROR**</FONT>";
			#	$anaStatus.=$delButtonStrg unless $anaStatus=~/deleteJob/; # delete button already added
			#	$anaStatus.=qq|&nbsp;<INPUT type="button" value="Show/Hide error" onclick="viewError('error:$jobDir:$quantItemID');"></B>|; # id="errorButton:$jobDir:$quantItemID"
			#	#$anaStatus.=qq|&nbsp;<INPUT type="button" value="Show/Hide Rout" onclick="viewError('rout:$jobDir:$quantItemID');"></B>| if -e "$quantifHomeDir/$jobDir/quanti_$quantItemID/$RoutFile"; # id="routButton:$jobDir:$quantItemID"
			#	$anaError = $errorString;
			#	$scanAgain=1;
			#}
			$anaStatus.=$errorDisplayStrg if $error;

			if (!$error && !-e "$quantifHomeDir/$jobDir/status_$quantItemID.out") {
				if (-e "$quantifHomeDir/current/$quantItemID\_$jobDir\_run.flag") { # delay between run flag and run dir creation OR during files move at the end
					$anaStatus='<B>Moving files...</B>';
					$scanAgain=1;
				}
				else {
					$anaStatus='<B>Finished?</B> [Run file deleted]';
					unless (glob "$quantifHomeDir/$jobDir/status_*.out") { # in case multiple quantifs in same jobs
						my $now=strftime("%Y%m%d%H%M%S",localtime);
						$anaStatus.="&nbsp;<INPUT type=\"button\" value=\"Clean job\" onclick=\"deleteJob($jobDir,'')\">" if $now-$jobDir > 3600; # 1h / Deletes jobDir ONLY, not quanti in DB!
					}
					$scanAgain=0;
				}
			}
		}
	}
	else {
		$anaStatus='<B>Finished!</B> [Run file deleted]';
		$scanAgain=0;
		my $now=strftime("%Y%m%d%H%M%S",localtime);
		$anaStatus.="&nbsp;<INPUT type=\"button\" value=\"Clean Job\" onclick=\"deleteJob($jobDir,'')\">" if $now-$jobDir > 300; # 5 min / Deletes jobDir ONLY, not quanti in DB!
	}
	return ($scanAgain,$anaStatus,$anaError,$RoutFile);
}

sub ajaxDeleteJob {
	print header(-charset=>'UTF-8',-type=>'text/plain',-cache_control=>"no-cache, no-store, must-revalidate"); # for AJAX in IE
	warningsToBrowser(1);

	my $jobDir=param('jobDir');
	my $quantItemID=param('quantItem');

	if ($quantItemID) {
		if (-e "$quantifHomeDir/$jobDir") {
			my ($runDir,$quantifID);
			my $quantifType='?';
			if (-d "$quantifHomeDir/$jobDir/quanti_$quantItemID") { # Design
				$quantifID=$quantItemID;
				#$runDir="$quantifHomeDir/$jobDir/quanti_$quantItemID";
				$runDir="quanti_$quantItemID";
			}
			else { # other
				opendir (DIR,"$quantifHomeDir/$jobDir");
				while (defined (my $content = readdir (DIR))) {
					if ($content=~/quanti_${quantItemID}_(\d+)/) {
						$quantifID=$1;
						$runDir=$content;
						last;
					}
				}
				close DIR;
			}
			unless ($quantifID) { # IMPORT_MQ or error at in launchQuantification.cgi ?
				my $section='';
				open (INFO,"$quantifHomeDir/$jobDir/quantif_info.txt");
				while (<INFO>) {
					if (/^TYPE=(\S+)/) {$quantifType=$1; next;}
					elsif (/^QUANTIFICATIONS:/) {$section='quantifications'; next;}
					elsif ($section eq 'quantifications' && /^(\d+)/) {
						$quantifID=$1;
						$runDir=''; # bug occured before creation
						last;
					}
				}
				close INFO;
				unless ($quantifID) { # for emPAI & SIN
					(my $anaID,$quantifID)=split(/\./,$quantItemID);
					$runDir=0;
				}
			}
#print "**jobDir=$jobDir\n**quantItemID=$quantItemID\n**runDir=$runDir\n**quantifID=$quantifID\n";

			##<Delete from DB
			if ($quantifType eq 'IMPORT_MQ') {
				rmtree("$quantifHomeDir/$jobDir"); unlink "$quantifHomeDir/$jobDir" if -e "$quantifHomeDir/$jobDir"; # empty dir can remain
				my $dbh=&promsConfig::dbConnect;
				my ($hasChild)=$dbh->selectrow_array("SELECT 1 FROM SAMPLE WHERE ID_EXPERIMENT=$quantItemID");
				$dbh->disconnect;
				if ($hasChild) {print "##DELETE_EXP\n";}
				else {print "##RELOAD\n";}
			}
			elsif ($quantifID && $quantifType ne 'XICCORR') { # quantif is recorded in DB BUT all data still in temp dir
#die "BAD: $quantifType!";
				my $dbh=&promsConfig::dbConnect;
				my $projID=&promsMod::getProjectID($dbh,$quantifID,'quantification');
				&promsQuantif::deleteQuantification($dbh,$projID,$quantifID);
				$dbh->commit;
				$dbh->disconnect;
			}
			##<Delete files
			unlink "$quantifHomeDir/$jobDir/status_$quantItemID.out";
		}
		#unlink glob "$quantifHomeDir/current/$quantItemID\_$jobDir\_*" if glob "$quantifHomeDir/current/$quantItemID\_$jobDir\_*";
		#unlink "$quantifHomeDir/current/$jobDir\_request.flag" if -e "$quantifHomeDir/current/$jobDir\_request.flag";
	}
	unlink glob "$quantifHomeDir/current/*_$jobDir\_*" if glob "$quantifHomeDir/current/*_$jobDir\_*";
	unlink "$quantifHomeDir/current/$jobDir\_request.flag" if -e "$quantifHomeDir/current/$jobDir\_request.flag";

	##<delete quantif directory if no job left
	if (-e "$quantifHomeDir/$jobDir") {
		#my $numFiles=`ls -l $quantifHomeDir/$jobDir | wc -l`;
		#chomp($numFiles);
		#if ($numFiles*1 <= 2 || !$quantItemID) { # total + only quantif_info.txt
			#remove_tree("$quantifHomeDir/$jobDir");
			rmtree("$quantifHomeDir/$jobDir"); unlink "$quantifHomeDir/$jobDir" if -e "$quantifHomeDir/$jobDir"; # empty dir can remain
			print "##RELOAD\n";
		#}
		##print '##OK';
		#else {
		#	print 'Could not delete all temporary files.';
		#}
	}
	else {
		#print 'Data directory not found.';
		print "##RELOAD (Data directory not found.')\n";
	}
	exit;
}

# TODO: Check also for error at cluster connection
####>Revision history<####
# 1.5.0 Adapted for MaxQuant import (PP 20/05/19) 
# 1.4.1 [Fix] minor bug fix in quantifType detection (PP 22/02/19)
# 1.4.0 Adapted for XICCORR quantif type (PP 19/02/19) 
# 1.3.9 Checks for error earlier in pipeline (PP 11/02/19)
# 1.3.8 [Fix] occasional JS bug following deletion request (PP 04/12/18)
# 1.3.7 Various improvements (PP 14/11/18)
# 1.3.6 Minor fix of undef variable (PP 30/10/18)
# 1.3.5 Slower window.reload frequency in case multi-job & fix bug modif-quantif deletion (PP 21/09/18)
# 1.3.4 Scans also for a new wait.flag file used for design-based multi-launch (PP 30/08/18)
# 1.3.3 Removes peptide data file(s) not from PEPTIDE_QUANTIFICATION & uses quanti_$quantifID instead of quantif_$quantifID (PP 18/05/18)
# 1.3.2 [Fix] Bug in directory scanning for MaassChroQ quantif (PP 13/02/18)
# 1.3.1 Also prevents window reload caused by New/Orphan jobs when user is reading an error message (PP 08/02/18)
# 1.3.0 Prevents AJAX error message refresh when user is reading it (PP 17/01/18)
# 1.2.9 Added refresh window button (PP 01/12/17)
# 1.2.8 Handles algo v3 & bug fix for full design-based quantification/job deletion at 1st click (PP 27/11/17)
# 1.2.7 Minor bug fix in quantification method detection and display (PP 24/10/17)
# 1.2.6 Minor change (SWATH -> DIA) (MLP 26/09/17)
# 1.2.5 Minor modification for internal quantification (GA 04/05/17)
# 1.2.4 Handles SWATH/MSstats & SSPA (PP 01/08/16)
# 1.2.3 Minor change (GA 04/02/16)
# 1.2.2 Minor changes (PP 25/03/15)
# 1.2.1 Add R.out DIV (GA 27/06/14)
# 1.2.0 Handles SuperRatio quantification (PP 10/06/14)
# 1.1.9 Minor update in &ajaxDeleteJob for emPAI and SIN (GA 16/04/14)
# 1.1.8 Better run duration estimation/display & non-fatal/time out error management (PP 04/04/14)
# 1.1.7 Uses rmtree instead of remove_tree (PP 10/03/14)
# 1.1.6 system command removal (PP 12/11/13)
# 1.1.5 Minor changes (PP 02/10/13)
# 1.1.4 Better handling of remaining job directories and files (PP 26/09/13)
# 1.1.3 Update for multi-line errors display (PP 09/09/13)
# 1.1.2 Refresh quantification frame after job deletion if displayed (PP 03/09/13)
# 1.1.1 Fix bug job deletion for label-free quantification (PP 30/08/13)
# 1.1.0 Displaying error message for failed quantifications (FY 27/08/13)
# 1.0.9 Added project detection on designID for TNPQ & bug fix in &promsMod::getItemInfo call (PP 11/07/13)
# 1.0.8 Add $quantiTime in &getAnaQuantifInfo for quanti time processing (GA 02/07/13)
# 1.0.7 Minor bug fix: uninitialized $error when no error (PP 19/06/13)
# 1.0.6 Checks if error file exists before grep (PP 30/04/12)
# 1.0.5 Job deletion option on error (PP 22/11/11)
# 1.0.4 cache-control on AJAX calls for IE (PP 26/10/11)
# 1.0.3 Error & labeled quantification management (PP 26/08/11)
# 1.0.2 Modified (GA ...)
