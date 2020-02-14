#!/usr/local/bin/perl -w

################################################################################
# editProjectItem.cgi    3.2.10                                                 #
# Authors: P. Poullet, G. Arras, F. Yvon, V. Sabatet (Institut Curie)          #
# Contact: myproms@curie.fr                                                    #
# Generates the interface allowing                                             #
# the user to create new project items in myProMS                              #
# (projects, experiments, gels, spots, samples or analysis).                   #
# Submits data to storeProjectItem.cgi.                                        #
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
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use POSIX qw(strftime); # to get the time
use Time::Piece;
use Archive::Tar;
use promsConfig;
use promsMod;
use promsQuantif;
use strict;
use File::Basename;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my %massSpecType=&promsConfig::getMsType;
my $maxRank = &promsConfig::getMaxRank;
my $userID=$ENV{'REMOTE_USER'};
my $date = strftime("%Y-%m-%d %H:%M:%S",localtime);

#print header; warningsToBrowser(1); # DEBUG
#############################
####>Fetching parameters<####
#############################
&reImportData if param('reval'); # a validated Analysis is being reValidated
my $action=param('ACT');
my $rowID=param('ID'); # if action=add => rowID=parentID (undef if item is a project)
$rowID=0 unless $rowID;
($rowID)=&promsMod::cleanParameters($rowID);
my $item=uc(param('ITEM')); # item being processed (child if 'add')

if ($action eq 'ajaxMissedCleavage') {&computeMissedCleavages;}
elsif ($action eq 'ajaxProteinSummary') {&proteinSummary;}
elsif ($action eq 'ajaxComputeFDR') {&computeFDR;}
elsif ($action eq 'downloadSearchFiles') {&downloadSearchFiles;}
elsif ($action eq 'updateMetadata') { &updateMetadata; }
my $parentItem=param('PARENT'); # only for add (except PROJECT ($item eq 'PROJECT')? '' : ($item eq 'EXPERIMENT')? 'PROJECT' : ($item eq 'GEL2D' || $item eq 'SAMPLE')? 'EXPERIMENT' : 'SAMPLE'; # for add only
#my ($item,$parentItem)=($action eq 'add' && $rowID)? (uc(param('CHILD')),param('ITEM')) : (uc(param('ITEM')),undef); # add project: rowID is null;
my $showParameters=(param('showParam'))? 1 : 0;


####################################
####>(Re)run identifier mapping<####
####################################
if ($action eq 'force') {
	print header(-charset=>'UTF-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<HEAD>
<TITLE>Re-Mapping Project</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<BODY>
<CENTER>
<BR><BR><FONT class="title2">Annotations mapping has been launched as a background task. It may take several minutes.
<BR><BR>You can continue using myProMS in the mean time. </FONT><INPUT type="button" class="title3" value=" Continue " onclick="parent.optionFrame.selectOption()"/>
</CENTER>
</BODY>
</HTML>
|;
	($userID,$rowID)=&promsMod::cleanParameters($userID,$rowID);
	system "./mapProteinIdentifiers.pl $userID $rowID 0 force";
	#window.location="./editProjectItem.cgi?ACT=summary&ITEM=PROJECT&ID=$rowID";
	exit;
}
#if ($action eq 'mapMissing' || $action eq 'mapAll') { # project only
#	$dbh->disconnect;
#	my $mappingDir="$promsPath{tmp}/mapping";
#	system "mkdir $promsPath{tmp}" unless -e "$promsPath{tmp}";
#	system "mkdir $mappingDir" unless -e $mappingDir;
#	my $jobListFile="$mappingDir/job_list.txt";
#	my $jobTimeFile="$mappingDir/job_time.txt";
#
#	#>Check time stamp
#	my $runningJobTime;
#	if (-e $jobTimeFile) {
#		$runningJobTime=`head -1 $jobTimeFile`;
#		chomp($runningJobTime);
#	}
#	my $safetyStrg='';
#	if (!$runningJobTime || time-$runningJobTime > 180) { # more than 3 min since last job was started => problem!
#		system "rm $jobListFile $jobTimeFile";
#		$safetyStrg='top.logoFrame.mappingJobRunning=0;'; # JS. just to be safe
#	}
#
#	#>Record new job
#	my $jobStartTime=strftime("%Y-%m-%d.%H:%M:%S",localtime); # just for info
#	open (JOBS,">>$jobListFile");
#	print JOBS "$rowID\tproject\t$rowID\t$userID\t$jobStartTime\n";
#	close JOBS;
#	print header(-charset=>'UTF-8');
#	warningsToBrowser(1);
#	print qq
#|<HTML><HEAD><SCRIPT type="text/javascript">
#$safetyStrg
#if (top.logoFrame.mappingJobRunning==0) {
#	top.logoFrame.wgetFrame.location="$promsPath{cgi}/wgetMapIdentifiers.cgi?ACT=$action";
#	top.logoFrame.mappingJobRunning=1;
#}
#window.location="./editProjectItem.cgi?ACT=summary&ITEM=PROJECT&ID=$rowID";
#</SCRIPT></HEAD></HTML>
#|;
#	exit;
#}


##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;
my %allPostTransModifs=&promsMod::getVariableModifications($dbh);

####>Checking Project mapping status<####
#if ($item eq 'PROJECT' && $action eq 'summary') {
#	my ($identMappingID)=$dbh->selectrow_array("SELECT IDENT_MAPPING FROM PROJECT WHERE ID_PROJECT=$rowID");
#	unless ($identMappingID) { # pre-existing project: IDENT_MAPPING=NULL
#		$dbh->do("UPDATE PROJECT SET IDENT_MAPPING='NONE' WHERE ID_PROJECT=$rowID");
#		$dbh->commit;
#		$action='mapMissing';
#	}
#}

#################################
####>Stop identifier mapping<#### Deletes all jobs programmed for current project
#################################
#if ($action eq 'stopMapping') { # project only
#	my $jobListFile="$promsPath{tmp}/mapping/job_list.txt";
#	my $tempFile="$promsPath{tmp}/mapping/job_list.tmp";
#	open (JOBS,$jobListFile);
#	open (TMP,">$tempFile") || die "PROBLEM!!";
#	my $otherJobs=0;
#	my $localUserOtherJobs=0;
#	while (<JOBS>) {
#		next if /^$rowID\s/; #skip selected project jobs
#		$localUserOtherJobs=1 if /\s$userID\s/; # local user has other on-going jobs
#		print TMP $_;
#		$otherJobs++;
#	}
#	close TMP;
#	close JOBS;
#	if ($otherJobs) {system "mv $tempFile $jobListFile";}
#	#else {system "rm $jobListFile $tempFile";}
#	else {system "rm -r $promsPath{tmp}/mapping";}
#	#$action='summary'; # Change action and continue!
#	print header(-charset=>'UTF-8');
#	warningsToBrowser(1);
#	print qq
#|<HTML><HEAD><SCRIPT type="text/javascript">
#top.logoFrame.mappingJobRunning=$localUserOtherJobs;
#window.location="./editProjectItem.cgi?ACT=summary&ITEM=PROJECT&ID=$rowID";
#</SCRIPT></HEAD></HTML>
#|;
#	exit;
#}

###################################
####>Fetching user information<####
###################################
# my ($userName,$userStatus,$refProfile,$userEmail,$userInfo)=&promsMod::getUserInfo($dbh,$userID);
# Temp -->
my ($userStatus,$userWorkgroup)=$dbh->selectrow_array("SELECT USER_STATUS,WORK_GROUP FROM USER_LIST WHERE ID_USER='$userID'");
$userWorkgroup='' unless $userWorkgroup;
# <--

############################
####>Fetching item info<####
############################
my $itemID;
my $itemType=&promsMod::getItemType($item);
my $itemSubType = &promsMod::cleanParameters(param("TYPE"));
my @itemInfo;
if ($action eq 'add') {
	###>Fetching ParentItem info
	@itemInfo=&promsMod::getItemInfo($dbh,$parentItem,$rowID) unless($item eq 'PROJECT');

	###>Completing refItemInfo with item data
	push @itemInfo,{'ITEM'=>$item,'TYPE'=>$itemType,'ID'=>undef,'NAME'=>''}; # new item id and name are undefined (SPOT->ANALYSIS!)
}
else { #if ($action=~/^(edit|summary)\Z/) { #} Edit or Summary
	$itemID=$rowID;
	print(($itemSubType) ? $itemSubType : $item);
	@itemInfo=&promsMod::getItemInfo($dbh, ($itemSubType) ? $itemSubType : $item, $itemID);
}
my ($projectID,$projectAccess); # not defined for Project addition
if ($itemInfo[0]{ID}) { # projectID is defined
	my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
	$projectID=$itemInfo[0]{'ID'};
	$projectAccess=$userInfo[2]->{$projectID};
}
my (@nameText,@valueText); # Globals
my $metaAna; # item is meta-analysis
my ($fileFormat,$decoy); # for analysis only     #,$fileName
my %infoSubmit ;
my $numParents=0;
my $hasValidAna=0;
my ($validStatus,$okScoreDistrib,@regressionPlots); # global, initialized in &analysis      ,$validDate
my ($missingProteins,$selIdentifierType,@identifierTypes); # global, initialized in &analysis
my $fdrDivString="<DIV id=\"fdrDIV\">&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Compute\" onclick=\"ajaxComputeFDR()\"/></DIV>";
if ($item eq 'PROJECT') {&project;}
elsif ($item eq 'ANALYSIS') {
	#$projectID=$itemInfo[0]{'ID'};
	#my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
	#$projectAccess=$userInfo[2]->{$projectID}; #${$userInfo[2]}{$projectID};
	&analysis;
	##>Checking for missing annotations
	if ($action eq 'summary' && $validStatus>=1 && $projectAccess ne 'guest') {
		#my ($vYear,$vMonth,$vDay,$vHour,$vMin,$vSec)=split(/-| |:/,$validDate);
		#my ($sec,$min,$hour,$mday,$month,$year,$wday,$yday,$isdst) = localtime(time); $year+=1900; $month++;
		#my $deltaValidTime=($year-$vYear)*525600 + ($month-$vMonth)*45360 + ($mday-$vDay)*1440 + ($hour-$vHour)*60 + ($min-$vMin); # delta time in minutes
		($missingProteins)=$dbh->selectrow_array("SELECT COUNT(PROTEIN.ID_PROTEIN) FROM PROTEIN,ANALYSIS_PROTEIN WHERE PROTEIN.ID_PROTEIN=ANALYSIS_PROTEIN.ID_PROTEIN AND (MW=0 OR ORGANISM='unknown organism' OR ORGANISM='' OR (PROT_SEQ='-' AND ID_MASTER_PROTEIN IS NOT NULL)) AND ID_ANALYSIS=$itemID");
		if ($missingProteins) {
			($selIdentifierType)=$dbh->selectrow_array("SELECT IDENTIFIER_TYPE FROM ANALYSIS_DATABANK AD,DATABANK D WHERE AD.ID_DATABANK=D.ID_DATABANK AND ID_ANALYSIS=$itemID LIMIT 1"); # first one only
			my $sthIT=$dbh->prepare("SELECT DISTINCT DEF_IDENT_TYPE FROM DATABANK_TYPE WHERE DEF_IDENT_TYPE IS NOT NULL");
			$sthIT->execute;
			while (my ($idenType)=$sthIT->fetchrow_array) {push @identifierTypes,$idenType;}
			$sthIT->finish;
		}
	}

	###>Showing Search parameters
	if ($showParameters==1) {
		%infoSubmit=promsMod::getSearchParam($dbh,$rowID);
	}
}
elsif ($item eq 'GEL2D') {&gel2D;}
elsif ($item eq 'SPOT') {&spot;}
else {&experiment_sample;} # experiment or sample

$dbh->disconnect;


################
####> HTML <####
################
my $titleString;
if ($action eq 'add') {
	$itemType.='(s)' if ($item eq 'EXPERIMENT' || $item eq 'SAMPLE');
	my $actionText=($item eq 'PROJECT')? 'Creating' : 'Adding';
	$titleString=$actionText;
	$titleString.=" a" if ($item ne 'EXPERIMENT' && $item ne 'SAMPLE');
	$titleString.=" new $itemType";
}
#elsif ($metaAna) {$titleString="Combined Analysis <FONT color=#DD0000>$itemInfo[-1]{NAME}</FONT>";}
else {$titleString="$itemType <FONT color=#DD0000>$itemInfo[-1]{NAME}</FONT>";}
print header(-'content-encoding'=>'no',-charset=>'UTF-8');
warningsToBrowser(1);
my ($light,$dark)=&promsConfig::getRowColors;
print qq
|<HTML>
<HEAD>
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.4.1/jquery.min.js"></script>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.8.2/css/all.css" integrity="sha384-oS3vJWv+0UjzBfQzYUhtDYW+Pj2yciDJxpsK1OYPAYjqT085Qq/1cq5FLXAZQ7Ay" crossorigin="anonymous">
<STYLE>
|;
foreach my $varMod (sort{$a cmp $b} keys %allPostTransModifs) {
	print ".$allPostTransModifs{$varMod}[2] \{font-weight:bold;$allPostTransModifs{$varMod}[3]\}\n";
}
print qq |
	.annot .value, .annot .desc, .meta .name, .meta .desc {
		display: inline-block;
	}
	
	.annot {
		padding: 6px;
	}
	
	.meta, .annotMetaLike {
		margin-top: 6px;
	}
	
	.meta > div {
		margin-bottom: 5px;
		display: inline-block;
	}
	
	.annot > div {
		display: inline-block;
	}
	
	.meta .name {
		font-weight: bold;
		text-decoration: underline;
		font-size: 12.5px;
	}
	
	.meta .desc {
		font-style: italic;
	}
	
	.annot .desc {
		font-style: normal;
	}
	
	.annot .value {
		//min-width: 120px; 
		padding-right: 15px;
	}

	.annot .desc {
		width: 100%;
	}
	
	button.action {
		width: 25px;
		height: 15px;
		font-size: 10px;
	}
	
	.actionsDiv {
		margin-right: 15px;
		min-width: 113px;
	}
	
	.imgFile {
		position: relative;
		top: 2px;
	}
	
	li.draggable-line {
		font-size: 11px;
	}
	
	li.annotation-item-row {
		padding: 1px 0px 3.5px 7px;
	}
	
	.annotMetaLike {
		padding-left: 0px !important;
	}
	
	#addMetaButton {
		margin: 5px 0 10px 6px;
	}
</STYLE>

<script type="text/javascript" src="$promsPath{html}/js/other/html5sortable.min.js"></script>
<SCRIPT type="text/javascript">
var popupWin;
function closeWindow() {
	if (popupWin != null) {
		if (!popupWin.closed) {
			popupWin.close();
		}
	}
}

// AJAX --->
function getXMLHTTP() {
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
|;
if ($action eq 'summary') {
	print qq
|
function ajaxMissedCleavage() {
	var mcSpan=document.getElementById('missCleavSPAN');
	mcSpan.innerHTML="<IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\" width=200>&nbsp;";
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/editProjectItem.cgi?ACT=ajaxMissedCleavage&ID=$rowID&ITEM=$item",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			mcSpan.innerHTML=XHR.responseText;
		}
	}
	XHR.send(null);
}
var proteinSummary=false; // flag for displayed or not protein summary
function ajaxProteinSummary() {
	if (proteinSummary) return; // in case called by setTimeout but already called by user
	proteinSummary=true;
	var protSpan2=document.getElementById('protSPAN2');
	protSpan2.innerHTML="<IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\" width=300>&nbsp;";

	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/editProjectItem.cgi?ACT=ajaxProteinSummary&ID=$rowID&ITEM=$item",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			var [totProt,visProt,pcMapped,numDbTypes]=XHR.responseText.split(',');
			var [warnStrg1,warnStrg2]=(numDbTypes*1 < 2)? ['',''] : ['<FONT style="color:#DD0000"><SUP>*</SUP></FONT>','<FONT style="font-size:11px;color:#DD0000"><BR>&nbsp;'+numDbTypes+' identifier types used: Different instances of the same protein may occur.</FONT>'];
			document.getElementById('protSPAN1').innerHTML='&nbsp;'+visProt+warnStrg1;
			var sStrg=(visProt*1 > 1)? 's' : '';
			protSpan2.innerHTML='&nbsp;Protein'+sStrg+' ('+totProt+' max) <INPUT type="button" class="font11" value="List ambiguities" onclick="window.location=\\'$promsPath{cgi}/listConflicts.cgi?TYPE=$item&ID=$rowID&ACT=ambiguity\\'"/> <INPUT type="button" class="font11" value="Visible & Hidden" onclick="window.location=\\'$promsPath{cgi}/listConflicts.cgi?TYPE=$item&ID=$rowID&ACT=switch\\'"/>'+warnStrg2;
			document.getElementById('protSPAN3').innerHTML='&nbsp;'+pcMapped;
			document.getElementById('protSPAN4').innerHTML='&nbsp;% Proteins mapped to biological resources.';
			document.getElementById('protMapTR').style.display='';
		}
	}
	XHR.send(null);
}
function ajaxComputeFDR() {
	var resDiv=document.getElementById('fdrDIV');
	resDiv.innerHTML="<IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\">";

	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/editProjectItem.cgi?ACT=ajaxComputeFDR&ID=$rowID&ITEM=$item",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			resDiv.innerHTML=XHR.responseText;
		}
	}
	XHR.send(null);
}

// <--- AJAX
|;
}
else {
	if ($item eq 'PROJECT') {
		print qq
|function updateWorkgroup(wkGr) {
	var visStatus=(wkGr=='#new#')? 'visible' : 'hidden';
	document.getElementById('newWorkgroup').style.visibility=visStatus;
}
|;
	}
	print qq
|function checkForm(myForm) {
	if ((myForm.name && !myForm.name.value) && (myForm.metaName && !myForm.metaName.value) && (myForm.ANNOT_VALUE && !myForm.ANNOT_VALUE.value)) {
		alert('Type a name for this $itemType.');
		return false;
	}
|;
	if ($item eq 'PROJECT') {
		print qq
|	var okVis=false;
	for (var i=0; i<myForm.protVis.length; i++) {
		if (myForm.protVis[i].checked==true) {
			okVis=true;
			break;
		}
	}
	if (!okVis) {
		alert('Select a Protein visibility rule');
		return false;
	}
	//workgroup
	if (myForm.workgroup.value=='#new#' && !myForm.newWorkgroup.value) {
		alert('Enter a valid Workgroup name.');
		return false;
	}
|;
	}
	elsif ($item eq 'SPOT') { # edit only
		print qq
|	if (myForm.assoSampleID.value==-1) {
		if (!myForm.newSampleName.value) {
			alert('Type a name for Associated Sample.');
			return false;
		}
		if (myForm.newSampleName.value>45) {
			alert('Name for Associated Sample is too long.');
			return false;
		}
	}
|;
	}
	elsif ($item eq 'GEL2D') {
		print qq
|	if (!myForm.image_file.value) {
		alert('Gel image file is missing!');
		return false;
	}
	if (!myForm.image_file.value.toUpperCase().match(/\\.JPE?G/)) {
		alert('Gel image file is not in jpeg format!');
		return false;
	}
|;
	}
#	elsif ($item eq 'ANALYSIS') {
#		print qq
#|	if (myForm.databankID.value==-1) {
#		alert('Select the sequence databank used for this analysis.');
#		return false;
#	}
#	if (!myForm.data_file.value) {
#		alert('Select a search data file (DAT or XML) for this analysis.');
#		return false;
#	}
#	var dirNames=new Array();
#	if (myForm.data_file.value.match(/\\\\/)) { // from Windows
#		dirNames=myForm.data_file.value.split(/\\\\/); // array
#	}
#	else {dirNames=myForm.data_file.value.split(/\\\//);} // from other OS
#	var fileName=dirNames[dirNames.length-1];
#	if (fileName.length > 50) {
#		if (!confirm('WARNING: Data file name will be truncated to 50 characters.\\n                                      Proceed?')) {
#			return false;
#		}
#	}
#	if ('$action'=='add') {
#		popupWindow();
#	}
#|;
#	}
	print qq
|	return true;
}
|;
#function popupWindow() {
#	popupWin=window.open('','atWorkWindow','width=500,height=200');
#	var doc=popupWin.document;
#	doc.open();
#	doc.writeln('<HTML><HEAD><TITLE>Data File Upload</TITLE>');
#	doc.writeln('<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">');
#	doc.writeln('</HEAD><BODY>');
#	doc.writeln('<CENTER>\\n<BR><BR><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">');
#	doc.writeln('<FONT style="font-size:18px;font-weight:bold;">Uploading data file...</FONT>');
#	doc.writeln('</CENTER></BODY></HTML>');
#	doc.close;
#}
	print qq
|function cancelAction(fallback) {
	top.promsFrame.selectedAction=fallback; // set default action to 'summary'
	top.promsFrame.optionFrame.selectOption();
}
|;
} # end of action ne 'summary'

if ($item eq 'ANALYSIS') { # && $action ne 'add'
	print qq
|var regPlotHidden=true;
function showHideRegressionPlots(numPlots) {
	var plotStrg=(numPlots > 1)? 'plots' : 'plot';
	var displayValue,buttonValue;
	if (regPlotHidden) {
		displayValue='';
		regPlotHidden=false;
		buttonValue='Hide regression '+plotStrg;
	}
	else {
		displayValue='none';
		regPlotHidden=true;
		buttonValue='Show regression '+plotStrg;
	}
	document.getElementById('regressionDIV').style.display=displayValue;
	document.getElementById('regressionBUTTON').value=buttonValue;
}
var scoreImgHidden=true;
function showHideScoresDistribution() {
	var displayValue,buttonValue;
	if (scoreImgHidden) {
		displayValue='';
		scoreImgHidden=false;
		buttonValue='Hide scores distribution';
	}
	else {
		displayValue='none';
		scoreImgHidden=true;
		buttonValue='Show scores distribution';
	}
	document.getElementById('scoresDIV').style.display=displayValue;
	document.getElementById('scoresBUTTON').value=buttonValue;
}
function getMissingProteins(identifierType) {
	if (!identifierType) return;
	mappingFrame.location="$promsPath{cgi}/editProjectItem.cgi?ITEM=ANALYSIS&ID=$itemID&ACT=mapMissing&IDENT="+identifierType;
	document.getElementById('mappingDIV').innerHTML='<FONT class="title3" color="#DD0000">Annotation retrieval is running in background.<BR>You can continue using myProMS.</FONT>';
}
|;
	if ($action eq 'summary') {
		print qq
|var historyWin;
function displayValidationHistory() {
	var winLocation='$promsPath{cgi}/displayValidationHistory.cgi?ID=$itemID';
	if (!historyWin \|\| historyWin.closed) {
		historyWin=window.open(winLocation,'historyWindow','width=900,height=800,location=no,resizable=yes,scrollbars=yes');
	}
	else {
		historyWin.location=winLocation;
	}
	historyWin.focus();
}
|;
	}
}
if ($action eq 'edit') {
		print qq
|function displayNewSample(sampID) { // spot only
	var newSampVis=(sampID==-1)? 'visible' : 'hidden';
	document.getElementById('newSample').style.visibility=newSampVis;
}

function delMetadataRow(row) {
	var parent = row.parentNode;
	
	parent.removeChild(row);
	
	var nbChilds = parent.children.length;
	if(nbChilds == 0 && parent.classList.contains("meta-annotation-list")) {
		var metaRow = parent.parentNode;
		metaRow.parentNode.removeChild(metaRow);
	}
}

function activateSortableMetadata() {
	sortable('.js-sortable', {
		"items": ":not(.disabled)"
	});
	
	sortable('.js-sortable-inner', {
		"items": ":not(.disabled)"
	});
	
	var annotationsSubLists = document.getElementsByClassName('js-sortable-inner');
	setUpdateMetadataAjax(annotationsSubLists, "sub");
	
	var annotationsGlobalLists = document.getElementsByClassName('js-sortable');
	setUpdateMetadataAjax(annotationsGlobalLists, "global");
}

function setUpdateMetadataAjax(items, type) {
	for(var i = 0; i < items.length; i++) {
		items[i].addEventListener('sortstop', function (e) {
			var element = e.detail.item
			var parent = element.parentNode;
			var currentIndex = 1;
			var positions = "";
			
			for(var y=0; y < parent.children.length; y++) {
				if(y != 0) {
					positions += ",";
				}
				
				var currentElementId = (parent.children[y].dataset.idMetadata) ? parent.children[y].dataset.idMetadata : parent.children[y].dataset.idAnnotation;
				
				positions += currentElementId + ":" + currentIndex;
				currentIndex++;
			}

			var XHR = getXMLHTTP();
			XHR.onreadystatechange=function () {
				if(XHR.readyState==4 && XHR.responseText ) {
					console.log("Updated all annotations for metadata " + parent.dataset.id);
				}
			}
			XHR.open("GET", "$promsPath{cgi}/editProjectItem.cgi?ITEM=$item&ID=$itemID&ACT=updateMetadata&metaType=" + type + "&positions=" + positions, true);
			XHR.send(null);
		});
	}
}

window.onload = activateSortableMetadata;
|;
}

print qq |
	function editBranch(action) {
		if (action=='end' && !confirm('Ending Project will also end all on-going Analysis validations. Proceed?')) {return;}
		else if (action=='archive' && !confirm('Archiving will disable access to Project. Proceed?')) {return;}
		else if (action=='continue' && !confirm('Continue Project?')) {return;}
		window.location="./editBranch.cgi?ITEM=$item&ID=$rowID&ACT="+action;
	}
	
	function switchFieldType() {
		var type = document.getElementById("ANNOT_TYPE").value;
		var valueField = document.getElementById("ANNOT_VALUE");
		var styleField = '';
		
		if(type == 'link' \|\| type == 'text') {
			styleField = 'width: 541px;';
		} else if(type == 'image') {
			type = 'file';
		} else {
			styleField = 'width: 250px;';
		}
		
		var html = "<input style='" + styleField + "' id='ANNOT_VALUE' name='ANNOT_VALUE' type='" + type + "' value='' />";
		
		valueField.insertAdjacentHTML('afterend', html);
		valueField.parentNode.removeChild(valueField);
	}
	
	function resetForm() {
		document.getElementById("addItemForm").reset();
		switchFieldType();
		
		var value = document.getElementById("currentValue").innerHTML;
		if(value !== null && value !== '') {
			value = value.substring(11, value.length-1);
		}
		
		var type = document.getElementById("ANNOT_TYPE").value;
		if(type == "date") {
			var date = value.split("/");
			value = date[2] + "-" + date[1] + "-" + date[0];
		}
		document.getElementById("ANNOT_VALUE").value = value;
	}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" onunload="closeWindow()">
<CENTER>
<FONT class="title">$titleString</FONT><BR><BR>
|;

######################
####>Writing form<####
######################
# my $userCode=(split(//,$projectAccess))[0];
my $actionString=($item eq 'ANALYSIS' && $action eq 'add')? 'storeAnalyses.cgi' : 'storeProjectItem.cgi';

if ($action ne 'summary') {
	print qq
|<FORM name="addItemForm" id="addItemForm" onsubmit="return(checkForm(this));" method="post" action="$actionString" enctype="multipart/form-data\">
<INPUT type="hidden" name="ACT" value="$action">
<INPUT type="hidden" name="ITEM" value="$item">
|;
	print "<INPUT type=\"hidden\" name=\"ITEM_ID\" value=\"$itemID\">\n" if $action eq 'edit'; 	# transmitting hidden fields with submit
	if ($action eq 'add' && $item ne 'PROJECT') {
		print "<INPUT type=\"hidden\" name=\"PARENT_ID\" value=\"$itemInfo[-2]{ID}\">\n"; # unless ($item eq 'PROJECT' || $action eq 'edit');#
		print "<INPUT type=\"hidden\" name=\"PARENT\" value=\"$itemInfo[-2]{ITEM}\">\n";
	}
	if ($item ne 'PROJECT') {
		print "<INPUT type=\"hidden\" name=\"PROJECT_ID\" value=\"$itemInfo[0]{ID}\">\n";
	}
}

####>Printing data in table<####
print "<TABLE border=0 width=850>\n";
print "<TR><TD bgcolor=$dark>";
print "<TABLE border=0 cellpadding=2 width=100%>\n";
for my $i (0..$#nameText) {
	my $width=($item eq 'ANALYSIS')? 220 : 180;
	print "<TR>\n\t<TH align=right valign=top bgcolor=$dark nowrap width=$width>&nbsp;$nameText[$i] :</TH>\n";
	print "\t<TD bgcolor=$light nowrap>$valueText[$i]</TD>\n</TR>\n";
}
if ($action ne 'summary') {
	print "<TR><TD colspan=2 align=center>\n";
	###>SUBMIT button
	print '<INPUT type="submit" name="save" value=" Save " >',"\n"; # style="color: #A0A0A0; background-color: #FF9999"
	###>CLEAR button
	my $resetString=($action eq 'add')? '  Clear  ' : ' Cancel changes ';
	my $onclickEvent = ($item eq 'METADATA') ? 'resetForm();' : 'reset();';
	print "&nbsp;&nbsp;&nbsp;<INPUT type='button' value='$resetString' onclick='$onclickEvent'>";
	###>CANCEL button
	my $urlString;
	if ($item eq 'PROJECT' && $action eq 'add') {
		$urlString='history.back()'; # 'top.promsFrame.location=\"./selectProject.cgi\"';
	}
	else {
		my $selectedItem;
		if ($action eq 'add') {
			$selectedItem=$parentItem;
		}
		else { # action=edit
			$selectedItem=$item;
		}
		my $cancelFallback = ($item eq 'METADATA') ? 'edit' : 'summary';
		$urlString="cancelAction('$cancelFallback');";
	}

	print "&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" value=\"Cancel\" onclick=\"$urlString\" >\n";
	print "</TD></TR>\n";
	print "</TABLE>\n";
	print "</TD></TR>\n";
	#if ($action eq 'add' && $item eq 'ANALYSIS') {
	#	print "<TR><TD><FONT style=\"font-size:11px;font-weight:bold;font-style:italic;\"><SUP>*</SUP>File must be in DAT or XML format.";
	#	print "&nbsp&nbsp&nbsp;°&nbsp;Only for Mascot MS/MS analysis.\n" ;
	#	print "</FONT></TD></TR>\n";
	#}
	print "</TABLE>\n";
	print "</FORM>\n";
}
else {print "</TABLE></TD></TR>\n</TABLE>\n";}


if ($action eq 'summary' && $item eq 'ANALYSIS' && $showParameters==1) {
	print qq
|<BR><BR>
<TABLE border=0 width=800>
<TR><TD bgcolor=$dark><TABLE border=0 cellpadding=2 width=100%>
<TR><TH colspan=2 bgcolor=$dark><FONT class="title2">Search Parameters</FONT>
</TH></TR>
|;
	foreach my $param (sort {$a cmp $b} keys %infoSubmit) {
		(my $trueParam=$param)=~s/^\w+://; # remove sort tag
		if (ref($infoSubmit{$param}) eq 'HASH') { # => multi-db search
			#>Fetching number of databanks searched
			my $lastIndex=0;
			foreach my $subParam (keys %{$infoSubmit{$param}}) {
				$lastIndex=$#{$infoSubmit{$param}{$subParam}};
				last;
			}
			foreach my $i (0..$lastIndex) {
				my $bdRankStrg=($lastIndex>=1)? ' #'.($i+1) : '';
				print "<TR valign=top><TH align=right bgcolor=$dark nowrap>&nbsp;$trueParam$bdRankStrg :</TH><TD nowrap bgcolor=$light>";
				foreach my $subParam (sort {$a cmp $b} keys %{$infoSubmit{$param}}) {
					(my $trueSubParam=$subParam)=~s/^\w://; # remove sort tag
					print "&nbsp;&bull;<B>$trueSubParam:</B>&nbsp;$infoSubmit{$param}{$subParam}[$i]<BR>";
				}
				print "</TD></TR>\n";
			}
		} else {
			print "<TR valign=top><TH align=right bgcolor=$dark nowrap>&nbsp;$trueParam :</TH><TD nowrap bgcolor=$light>&nbsp;$infoSubmit{$param}</TD></TR>\n";
		}
	}
	unless (scalar keys %infoSubmit) {
		print "<TR><TH align=left bgcolor=$light colspan=2>No search parameters found.</TH></TR>\n";
	}
	print "</TABLE>\n</TD></TR></TABLE>\n";
}

if ($missingProteins) { # validated analysis in summary mode
	my $protString=($missingProteins>1)? 's were' : ' was';
	print qq
|<BR><DIV id="mappingDIV">
<FONT class="title3">$missingProteins protein$protString found with incomplete annotation.<BR>
|;
	if (-e "$promsPath{tmp}/mapping.job") {
		print "An annotation retrieval is already on-going. Wait for job to end and check again.</FONT></DIV>\n";
	}
	else {
		print qq
|myProMS can attempt to retrieve the missing annotations from the internet<BR>
using the following identifier type :</FONT><SELECT id="identType" class="title3"><OPTION value="">-= Select =-</OPTION>
|;
		foreach my $identType (sort @identifierTypes) {
			print "<OPTION value=\"$identType\"";
			print ' selected' if ($selIdentifierType && $identType eq $selIdentifierType);
			print ">$identType</OPTION>\n";
		}
		print qq
|</SELECT><INPUT type="button" value=" Proceed " onclick="getMissingProteins(document.getElementById('identType').value)">
</DIV><BR>
<IFRAME name="mappingFrame" src="" width=150 height=40 frameborder=0 style="display:none"></IFRAME>
|;
	}
}

if ($action eq 'summary' && $item eq 'PROJECT' && $itemInfo[0]{'STATUS'} > 0 && $projectAccess=~/bioinfo|mass|manag/) {
	print "<TABLE cellpadding=4><TR>";
	print "<TD bgcolor=\"#00DD00\"><INPUT type=\"button\" class=\"title3\" value=\"Continue Project\" style=\"width:170px;\" onclick=\"editBranch('continue')\"></TD>";
	print "<TD>&nbsp;&nbsp;</TD><TD bgcolor=\"#DD0000\"><INPUT type=\"button\" class=\"title3\" value=\"Archive Project\" style=\"width:170px;\" onclick=\"editBranch('archive')\"></TD>" if $itemInfo[0]{STATUS}==1; # already ended
	print "</TR></TABLE>\n";
}
if ($action eq 'edit') { # edit => user is not a guest & project status=0
	print "<TABLE cellpadding=4><TR>\n";
	if ($item eq 'PROJECT') { # status<=0
		if ($projectAccess=~/bioinfo|mass|manag/) {
			print "<TD bgcolor=\"#0000DD\"><INPUT type=\"button\" class=\"title3\" value=\"End Project\" style=\"width:170px;\" onclick=\"editBranch('end')\"></TD>";
			print "<TD>&nbsp;&nbsp;</TD><TD bgcolor=\"#DD0000\"><INPUT type=\"button\" class=\"title3\" value=\"Archive Project\" style=\"width:170px;\" onclick=\"editBranch('archive')\"></TD>";
		}
	}
	else {
		if ($item ne 'SPOT' && $numParents>1) { # item can be moved to a different parent
			print "<TD><INPUT type=\"button\" value=\"Move $itemType\" style=\"width:170px;\" onclick=\"editBranch('move')\"></TD>\n";
			print "<TD width=30></TD>\n";
		}
		if ($item ne 'ANALYSIS' && $item ne 'METADATA') {
			print "<TD><INPUT type=\"button\" value=\"Delete $itemType\" style=\"width:170px;\" onclick=\"editBranch('delete')\"></TD>";
		}
		#elsif ($validStatus==2 && !$metaAna && $fileFormat eq 'MASCOT.DAT' && $decoy !~ /EXT:/ && $projectAccess=~/bioinfo|mass|manag|super/) {
		#	print "<TD><INPUT type=\"button\" value=\"Restart Validation\" style=\"width:170px;\" onclick=\"window.location='./editProjectItem.cgi?ID=$itemID&reval=1'\"></TD>";
		#}
	}
	print "</TR></TABLE>\n";
}

##>Score distribution<##
if ($item eq 'ANALYSIS') {
	my $anaDirHtml=($validStatus <= 1)? "$promsPath{data_html}/validation/ana_$itemID" : "$promsPath{data_html}/peptide_data/proj_$projectID/ana_$itemID";
	if (scalar @regressionPlots) {
		print "<DIV id=\"regressionDIV\" style=\"display:none\">";
		foreach my $regData (@regressionPlots) {
			my ($refRTname,$pepQuantName,$refRT_ID,$pepQuantifID)=@{$regData};
			my $regImg="regression_RT$refRT_ID"."_Q$pepQuantifID.png";
			print "<BR><FIELDSET style=\"width:50%\"><LEGEND class=\"title3\">Regression Plot for $refRTname w/ $pepQuantName:</LEGEND><IMG style=\"border:solid 3px $dark\" src=\"$anaDirHtml/$regImg\"/></FIELDSET>\n";
		}
		print "</DIV>\n";
	}
	if ($okScoreDistrib) {
		#my $scoreImage=($validStatus <= 1)? "$promsPath{data_html}/validation/ana_$itemID/scores.png" : "$promsPath{data_html}/peptide_data/proj_$projectID/ana_$itemID/scores.png";
		print "<DIV id=\"scoresDIV\" style=\"display:none\"><BR><IMG style=\"border:solid 3px $dark\" src=\"$anaDirHtml/scores.png\"/><BR></DIV>\n";
	}
}
if ($action eq 'summary' && $item ne 'ANALYSIS' && $hasValidAna) { # auto display Protein summary
	print qq
|<SCRIPT type="text/javascript">
setTimeout(ajaxProteinSummary,10000);
</SCRIPT>
|;
}
print "</CENTER>\n</BODY>\n</HTML>\n";


##############################<<< SUBROUTINES >>>##############################

#################
####<Project>####
#################
sub project {
	@nameText=('Name','Description','Protein visibility','Identifier conversion','Relevant PTMs','Project owner','Workgroup','Start date','Status');
	push @nameText,('Summary','False discovery') if $action eq 'summary';
	push @nameText,'Comments';
	push @nameText,'Metadata' if $action ne 'add';
	my @protVisText=('A protein is <B>Visible</B> only when <B>Alias</B> of a Match Group.',
					 'A protein is <B>Visible</B> everywhere if <B>Alias</B> of at least 1 Match Group.',
					 'A protein is <B>Visible</B> everywhere if <B>Alias</B> or made <B>Visible</B> in at least 1 Match Group.',
					 'A protein is <B>Visible</B> if it has the maximum number of peptides in a Match Group.'
					);
	#my ($name,$description,$protVisibility,$identMappingID,$ptmString,$owner,$workgroup,$startDate,$status,$comments);
	my ($name,$description,$protVisibility,$identMappingID,$owner,$workgroup,$startDate,$status,$comments);
	my $statusStrg;
	my %relevantPTMs;

	#<Identifier conversion
	my %identifierMapping;
	my $sthIC=$dbh->prepare("SELECT ID_IDENTIFIER,NAME FROM IDENTIFIER");
	$sthIC->execute;
	while (my ($identID,$identName)=$sthIC->fetchrow_array) {
		$identifierMapping{$identID}=$identName;
	}
	$sthIC->finish;

	if ($action ne 'add') { # summary or edit
		#($name,$description,$protVisibility,$identMappingID,$ptmString,$owner,$workgroup,$startDate,$status,$comments)=$dbh->selectrow_array("SELECT NAME,DES,PROT_VISIBILITY,ID_IDENTIFIER,RELEVANT_PTMS,OWNER,WORK_GROUP,START_DATE,STATUS,COMMENTS FROM PROJECT WHERE ID_PROJECT=$itemID"); # pep_thres,
		($name,$description,$protVisibility,$identMappingID,$owner,$workgroup,$startDate,$status,$comments)=$dbh->selectrow_array("SELECT NAME,DES,PROT_VISIBILITY,ID_IDENTIFIER,OWNER,WORK_GROUP,START_DATE,STATUS,COMMENTS FROM PROJECT WHERE ID_PROJECT=$itemID"); # pep_thres,
		$status=0 unless $status; # same as $itemInfo[0]{STATUS}
		##> Deprecated since 06/05/2013
		##<Check list of PTMs used
		#if ($ptmString) {
		#	foreach my $ptm (split(';',$ptmString)) {$relevantPTMs{$ptm}=$allPostTransModifs{$ptm}[0];}
		#}
		#<Check list of PTMs used
		my $sthGetPM=$dbh->prepare("SELECT ID_MODIFICATION FROM PROJECT_MODIFICATION WHERE ID_PROJECT=$itemID");
		$sthGetPM->execute;
		while (my ($modID)=$sthGetPM->fetchrow_array) {
			$relevantPTMs{$modID}=$allPostTransModifs{$modID}[0];
		}
		$sthGetPM->finish;
		if ($status <= 0) {
			$statusStrg='Ongoing';
			if ($userStatus eq 'bioinfo') {
				if ($action eq 'summary') {
					$statusStrg.='&nbsp;&nbsp;(Auto-end validations not allowed)' if $status==-1;
				}
				elsif ($action eq 'edit') {
					my ($act,$butStrg)=($status==0)? ('noAutoEndVal','No auto-end validations') : ('autoEndVal','Auto-end validations');
					$statusStrg.="&nbsp;&nbsp;<INPUT type=\"button\" value=\"$butStrg\" onclick=\"editBranch('$act')\">";
				}
			}
		}
		elsif ($status==1) {$statusStrg='Ended';}
		else {$statusStrg='Archived';}

	}
	else { # add
		($owner)=$dbh->selectrow_array("SELECT USER_NAME FROM USER_LIST WHERE ID_USER='$userID'");
		#$workgroup=$userWorkgroup;
		$startDate=$date;
		$statusStrg='Starting';
	}
	$identMappingID=0 unless $identMappingID;

	($name,$description,$owner,$workgroup,$comments)=&promsMod::chkDef($name,$description,$owner,$workgroup,$comments);
	$protVisibility=0 unless $protVisibility;
	$startDate=&promsMod::formatDate($startDate);
	if ($action eq 'summary') {
		$name=&promsMod::HTMLcompatible($name);
		push @valueText,"&nbsp;$name";
		$description=&promsMod::HTMLcompatible($description);
		push @valueText,"&nbsp;$description"; # description [1]
		push @valueText,"&nbsp;$protVisText[$protVisibility]"; # [2]
		my $identMapStrg=($identMappingID)? $identifierMapping{$identMappingID} : 'None';
		push @valueText,"&nbsp;$identMapStrg"; # identifier mapping [3]
		my @ptmList;
		foreach my $ptm (sort{lc($relevantPTMs{$a}) cmp lc($relevantPTMs{$b})} keys %relevantPTMs) {
			push @ptmList,"$allPostTransModifs{$ptm}[0] (<FONT class=\"$allPostTransModifs{$ptm}[2]\">$allPostTransModifs{$ptm}[1]</FONT>)";
		}
		my $ptmDisplayString=join(', ',@ptmList);
		$ptmDisplayString='None' unless $ptmDisplayString;
		push @valueText,"&nbsp;$ptmDisplayString"; # relevant PTMs
		push @valueText,"&nbsp;$owner"; # owner [4]
		$workgroup='None' unless $workgroup;
		push @valueText,"&nbsp;<B>$workgroup</B>"; # workgroup [5]
		push @valueText,"&nbsp;$startDate"; # start date [6]
		push @valueText,"&nbsp;<B>$statusStrg</B>"; # [7]
		# Summary
		my ($numGel2D,$numSpot,$numSamp,$numAna)=(0,0,0,0); #,$totProt,$visProt,$pcMapped)=(0,0,0);
		my ($numExp)=$dbh->selectrow_array("SELECT COUNT(ID_EXPERIMENT) FROM EXPERIMENT WHERE ID_PROJECT=$itemID");
		if ($numExp) {
			($numGel2D)=$dbh->selectrow_array("SELECT COUNT(ID_GEL2D) FROM EXPERIMENT,GEL2D WHERE EXPERIMENT.ID_EXPERIMENT=GEL2D.ID_EXPERIMENT AND ID_PROJECT=$itemID");
			if ($numGel2D) {
				($numSpot)=$dbh->selectrow_array("SELECT COUNT(ID_SPOT) FROM EXPERIMENT,GEL2D,SPOT WHERE EXPERIMENT.ID_EXPERIMENT=GEL2D.ID_EXPERIMENT AND GEL2D.ID_GEL2D=SPOT.ID_GEL2D AND ID_PROJECT=$itemID");
			}
			($numSamp)=$dbh->selectrow_array("SELECT COUNT(ID_SAMPLE) FROM EXPERIMENT,SAMPLE WHERE EXPERIMENT.ID_EXPERIMENT=SAMPLE.ID_EXPERIMENT AND ID_PROJECT=$itemID");
			if ($numSamp) {
				($numAna,$hasValidAna)=$dbh->selectrow_array("SELECT COUNT(ID_ANALYSIS),MAX(VALID_STATUS) FROM EXPERIMENT,SAMPLE,ANALYSIS WHERE EXPERIMENT.ID_EXPERIMENT=SAMPLE.ID_EXPERIMENT AND SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND ID_PROJECT=$itemID");
			}
		}
		#($totProt)=$dbh->selectrow_array("SELECT COUNT(*) FROM PROTEIN WHERE ID_PROJECT=$itemID"); # check even if no items in project in case of bug (only at project level)
		#if ($totProt) {
		#	($visProt)=$dbh->selectrow_array("SELECT COUNT(DISTINCT PROTEIN.ID_PROTEIN) FROM PROTEIN,ANALYSIS_PROTEIN WHERE PROTEIN.ID_PROTEIN=ANALYSIS_PROTEIN.ID_PROTEIN AND VISIBILITY>0 AND ID_PROJECT=$itemID");
		#	my ($numMapped)=$dbh->selectrow_array("SELECT COUNT(*) FROM PROTEIN WHERE ID_PROJECT=$itemID AND ID_MASTER_PROTEIN IS NOT NULL");
		#	$pcMapped=1*(sprintf "%.1f",($numMapped/$totProt)*100);
		#}
		my $pepModifStrg = ($hasValidAna)? '&nbsp;<SPAN id="missCleavSPAN"><INPUT type="button" class="font11" value="% Missed-cleavages" onclick="ajaxMissedCleavage()"/></SPAN>&nbsp;' : '';
		my $expString=($numExp>1)? "<TH align=right>&nbsp;$numExp</TH><TH align=left>&nbsp;Experiments</TH>" : "<TH align=right>&nbsp;$numExp</TH><TH align=left>&nbsp;Experiment</TH>";
		my $gel2dString=($numGel2D>1)? "<TH align=right>&nbsp;$numGel2D</TH><TH align=left>&nbsp;2D Gels</TH>" : "<TH align=right>&nbsp;$numGel2D</TH><TH align=left>&nbsp;2D Gel</TH>";
		my $spotString=($numSpot>1)? "<TH align=right>&nbsp;$numSpot</TH><TH align=left>&nbsp;Spots</TH>" : "<TH align=right>&nbsp;$numSpot</TH><TH align=left>&nbsp;Spot</TH>";
		my $sampString=($numSamp>1)? "<TH align=right>&nbsp;$numSamp</TH><TH align=left>&nbsp;Samples</TH>" : "<TH align=right>&nbsp;$numSamp</TH><TH align=left>&nbsp;Sample</TH>";
		my $anaString=($numAna>1)? "<TH align=right valign=top>&nbsp;$numAna</TH><TH align=left>&nbsp;Analyses $pepModifStrg</TH>" : "<TH align=right valign=top>&nbsp;$numAna</TH><TH align=left>&nbsp;Analysis $pepModifStrg</TH>";
		#my $sStrg=($visProt>1)? 's' : '';
		#my $protString="<TR><TH align=right>&nbsp;$visProt</TH><TH align=left>&nbsp;Protein$sStrg ($totProt max)&nbsp;".&getAmbiguityButton."&nbsp;".&getConflictButton."</TH></TR>";
		#my $mapString="<TR><TH align=right>&nbsp;$pcMapped</TH><TH align=left>&nbsp;% Proteins mapped to biological resources.".&checkOnGoingMapping."</TH></TR>";
		my ($protString,$mapString)=('','');
		if ($hasValidAna) {
			$protString="<TR><TH align=right valign=top><SPAN id=\"protSPAN1\"></SPAN></TH><TH align=left><SPAN id=\"protSPAN2\"><INPUT type=\"button\" class=\"font11\" value=\"Protein summary\" onclick=\"ajaxProteinSummary()\"/></SPAN></TH></TR>";
			$mapString="<TR id=\"protMapTR\" style=\"display:none\"><TH align=right><SPAN id=\"protSPAN3\"></SPAN></TH><TH align=left><SPAN id=\"protSPAN4\"></SPAN>".&checkOnGoingMapping."</TH></TR>";
		}
		my ($numBioSamples)=$dbh->selectrow_array("SELECT COUNT(ID_BIOSAMPLE) FROM PROJECT_BIOSAMPLE WHERE ID_PROJECT=$itemID");
		my $bioSampWord=($numBioSamples >= 2)? 'samples' : 'sample';
		my $bioSampString="<TH align=right>&nbsp;$numBioSamples</TH><TH align=left>&nbsp;Biological $bioSampWord referenced.</TH>";
		push @valueText,"<TABLE border=0 cellpadding=0 cellspacing=0><TR>$expString</TR><TR>$gel2dString</TR><TR>$spotString</TR><TR>$sampString</TR><TR>$anaString</TR>$protString$mapString<TR>$bioSampString</TR></TABLE>"; # [8]
		push @valueText,$fdrDivString; # [9]
		$comments=&promsMod::HTMLcompatible($comments);
		push @valueText,"&nbsp;$comments"; # [10]
		push @valueText, getMetadata($item, $itemID, $action) if($item ne 'METADATA');
	}
	else { # add or edit
		push @valueText,"<INPUT type='text' name='name' value='$name' size='50' maxlength='50'>";
		push @valueText,"<TEXTAREA name='des' rows='2' cols='65' maxlength='100'>$description</TEXTAREA>";
		my $protVisStrg='';
		foreach my $i (0..3) {
			my $chkStrg=($i==$protVisibility)? ' checked' : '';
			$protVisStrg.="<INPUT type='radio' name='protVis' value='$i'$chkStrg> $protVisText[$i]<BR>\n";
		}
		if ($action eq 'edit') {
			$protVisStrg.=qq
|&nbsp;<B>Or</B><BR><INPUT type='checkbox' name='upVis' value='1'> Update all proteins' visibility based on current rule.<BR>
&nbsp;<I><B>Warning:</B><BR>
&nbsp;&nbsp;&nbsp;&nbsp;-Changing current selection may modify <B>Custom lists</B> contents.<BR>
&nbsp;&nbsp;&nbsp;&nbsp;-Analyses with protein quantification or GO data will not be updated.</I>
<INPUT type='hidden' name='oldProtVis' value='$protVisibility'/>|;
		}
		push @valueText,$protVisStrg;
		my $identMapStrg="<SELECT name=\"mapping\"><OPTION value=\"\">None</OPTION>";
		foreach my $identID (2,1,3,10,50,60,110,120,130) { # Warning: IDs are hard coded here!!!
			$identMapStrg.="<OPTION value=\"$identID\"";
			$identMapStrg.=' selected' if $identID==$identMappingID;
			$identMapStrg.=">$identifierMapping{$identID}";
			$identMapStrg.=' (restrict to)' if $identifierMapping{$identID} eq 'GI Accession';
			$identMapStrg.="</OPTION>";
		}
		$identMapStrg.="</SELECT>";
		push @valueText,$identMapStrg;
		my @ptmList=(sort{lc($allPostTransModifs{$a}[1]) cmp lc($allPostTransModifs{$b}[1])} keys %allPostTransModifs);
		my $numRows=(scalar @ptmList)/3;
		$numRows++ if (int($numRows)<$numRows);
		$numRows=int($numRows);
		my $ptmIndex=0;
		my $ptmDisplayString="<TABLE cellspacing=0>";
		foreach my $row (1..$numRows) {
			$ptmDisplayString.="<TR>";
			foreach my $col (1..3) {
				if ($ptmIndex<=$#ptmList) {
					my $chkStrg=($relevantPTMs{$ptmList[$ptmIndex]})? ' checked' : '';
					$ptmDisplayString.="<TD nowrap><INPUT type=\"checkbox\" name=\"relevantPTMs\" value=\"$ptmList[$ptmIndex]\"$chkStrg>$allPostTransModifs{$ptmList[$ptmIndex]}[0] (<FONT class=\"$allPostTransModifs{$ptmList[$ptmIndex]}[2]\">$allPostTransModifs{$ptmList[$ptmIndex]}[1]</FONT>)&nbsp;</TD>";
				}
				else {$ptmDisplayString.="<TD></TD>";}
				$ptmIndex++;
			}
			$ptmDisplayString.="</TR>\n";
		}
		$ptmDisplayString.="</TABLE>\n";
		push @valueText,$ptmDisplayString;
		push @valueText,"<INPUT type='text' name='owner' value='$owner' size='30'>";
		my $workgroupStrg;
		if ($userStatus eq 'manag') {
			$workgroupStrg=($userWorkgroup)? "<B>$userWorkgroup</B>" : '<B>None</B>';
			$workgroupStrg.="<INPUT type=\"hidden\" name=\"workgroup\" value=\"$userWorkgroup\"/>";
		}
		else { # bioinfo|mass
			my %workgroupList;
			#Users
			my $sthWG1=$dbh->prepare('SELECT DISTINCT(WORK_GROUP) FROM USER_LIST WHERE WORK_GROUP IS NOT NULL');
			$sthWG1->execute;
			while (my ($wkGr)=$sthWG1->fetchrow_array) {
				next unless $wkGr;
				$workgroupList{$wkGr}=1;
			}
			$sthWG1->finish;
			#Projects
			my $sthWG2=$dbh->prepare('SELECT DISTINCT(WORK_GROUP) FROM PROJECT WHERE WORK_GROUP IS NOT NULL');
			$sthWG2->execute;
			while (my ($wkGr)=$sthWG2->fetchrow_array) {
				next unless $wkGr;
				$workgroupList{$wkGr}=1;
			}
			$sthWG2->finish;
			$workgroupStrg="<SELECT name=\"workgroup\" class=\"title3\" onchange=\"updateWorkgroup(this.value)\"><OPTION value=''>None</OPTION>\n";
			foreach my $wkGr (sort{lc($a) cmp lc($b)} keys %workgroupList) {
				$workgroupStrg.="<OPTION value=\"$wkGr\"";
				$workgroupStrg.=' selected' if $wkGr eq $workgroup;
				$workgroupStrg.=">$wkGr</OPTION>\n";
			}
			$workgroupStrg.=qq
|<OPTION value="#new#">-= New Workgroup =-</OPTION>
</SELECT>&nbsp;<INPUT type="text" name="newWorkgroup" id="newWorkgroup" value="" style="visibility:hidden"/>&nbsp;
|;
		}
		push @valueText,$workgroupStrg;
		push @valueText,"&nbsp;$startDate";
		push @valueText,"&nbsp;<B>$statusStrg</B>";
		push @valueText,"<TEXTAREA name='comments' rows='5' cols='65' maxlength='250'>$comments</TEXTAREA>";
		push @valueText, getMetadata($item, $itemID, $action) if ($action eq 'edit' && $item ne 'METADATA');
	}
}

##############################
####<Experiment or Sample>####
##############################
sub experiment_sample {
	my ($name, $description, $displayPos, $startDate, $comments, $prefSpeciesID, $prefSpeciesStrg);
	my ($metaName, $metaDes, $annotType, $annotValue, $accessibility, $metadataOwner) = ('', '', '', '', '', '');
	if ($action ne 'add') { # summary or edit
		@nameText=('Name','Description','Position','Start date');
		push @nameText,'Preferred species' if $item eq 'EXPERIMENT';
		push @nameText,('Summary','False discovery') if($action eq 'summary');
		push @nameText,'Comments';
		push @nameText,'List' if $action eq 'summary';
		push @nameText, 'Experiment lock' if($item eq 'EXPERIMENT' && $action eq 'edit' && $projectAccess =~ /bioinfo|mass/);
		push @nameText,'Metadata';
		if ($item eq 'SAMPLE') {($name,$description,$displayPos,$startDate,$comments) = $dbh->selectrow_array("SELECT NAME,DES,DISPLAY_POS,START_DATE,COMMENTS FROM SAMPLE WHERE ID_SAMPLE=$itemID");}
		elsif($item eq 'METADATA') {
			($metaName, $metaDes, $annotType, $annotValue, $comments, $accessibility, $metadataOwner) = $dbh->selectrow_array("SELECT NAME, M.DES, ANNOT_TYPE, ANNOT_VALUE, A.DES, ACCESSIBILITY, RECORD_USER FROM META_ANNOTATION M INNER JOIN ANNOTATION_ITEM A ON A.ID_META_ANNOTATION = M.ID_META_ANNOTATION WHERE ID_ANNOTATION_ITEM=$itemID");
			@nameText=('Category name', 'Category description', 'Annot. accessibility', 'Annot. type', 'Annot. value', 'Annot. comments', ) if($metadataOwner eq $userID);
		} else {
			($name,$description,$displayPos,$startDate,$comments,$prefSpeciesID,my $speciesSc,my $speciesCom)=$dbh->selectrow_array("SELECT NAME,DES,DISPLAY_POS,START_DATE,COMMENTS,E.ID_SPECIES,SCIENTIFIC_NAME,COMMON_NAME FROM EXPERIMENT E LEFT JOIN SPECIES S ON E.ID_SPECIES=S.ID_SPECIES WHERE ID_EXPERIMENT=$itemID");
			$prefSpeciesStrg=($prefSpeciesID)? "$speciesSc ($speciesCom)" : 'None';
		}
	}
	else { # add
		if($item eq 'METADATA') {
			if($itemSubType && $itemSubType eq 'ANNOTATION_ITEM') {
				@nameText=('Annot. type', 'Annot. value', 'Annot. comments');
			} else {
				@nameText=('Category name', 'Category description', 'Annot. accessibility', 'Annot. type', 'Annot. value', 'Annot. comments');
			}
		} elsif ($item eq 'EXPERIMENT' || $itemInfo[-2]{'ITEM'} eq 'EXPERIMENT') {
			@nameText=('Name','Multiple entries labels','Description','Start date');
			push @nameText,'Preferred species' if $item eq 'EXPERIMENT';
			push @nameText,'Comments';
		}
		else { # SAMPLE with not EXPERIMENT as parent
			@nameText=('Name','Description','Start date','Comments');
		}
		$startDate=$date;
	}
	
	($name,$description,$comments)=&promsMod::chkDef($name,$description,$comments);
	$startDate=&promsMod::formatDate($startDate);
	if ($action eq 'summary') {
		my $summaryString;
		if ($item eq 'SAMPLE') {
			#my ($totProt,$visProt,$pcMapped)=(0,0,0);
			(my $numAna,$hasValidAna)=$dbh->selectrow_array("SELECT COUNT(ID_ANALYSIS),MAX(VALID_STATUS) FROM ANALYSIS WHERE ID_SAMPLE=$itemID");
			#if ($numAna) {
			#	($totProt)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM ANALYSIS_PROTEIN,ANALYSIS WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ID_SAMPLE=$itemID");
			#	if ($totProt) {
			#		($visProt)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM ANALYSIS_PROTEIN,ANALYSIS WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND VISIBILITY>0 AND ID_SAMPLE=$itemID");
			#		my ($numMapped)=$dbh->selectrow_array("SELECT COUNT(DISTINCT P.ID_PROTEIN) FROM PROTEIN P,ANALYSIS_PROTEIN AP,ANALYSIS A WHERE ID_SAMPLE=$itemID AND P.ID_PROTEIN=AP.ID_PROTEIN AND AP.ID_ANALYSIS=A.ID_ANALYSIS AND ID_MASTER_PROTEIN IS NOT NULL");
			#		$pcMapped=1*(sprintf "%.1f",($numMapped/$totProt)*100);
			#	}
			#}
			my $pepModifStrg = ($hasValidAna)? '&nbsp;<SPAN id="missCleavSPAN"><INPUT type="button" class="font11" value="% Missed-cleavages" onclick="ajaxMissedCleavage()"/></SPAN>&nbsp;'.&getPtmDistributionButton() : '';
			my $exportAnalysisStrg = ($hasValidAna)? &getExportAnalysisButton : '';
			my $anaString=($numAna>1)? "<TH align=right valign=top>&nbsp;$numAna</TH><TH align=left>&nbsp;Analyses $pepModifStrg $exportAnalysisStrg</TH>" : "<TH align=right valign=top>&nbsp;$numAna</TH><TH align=left>&nbsp;Analysis $pepModifStrg $exportAnalysisStrg</TH>";
			#my $sStrg=($visProt>1)? 's' : '';
			#my $protString="<TH align=right>&nbsp;$visProt</TH><TH align=left>&nbsp;Protein$sStrg ($totProt max)";
			#$protString.="&nbsp;".&getAmbiguityButton."&nbsp;".&getConflictButton if $totProt;
			#$protString.="&nbsp;</TH>";
			##my $geneString=($numGenes>1)? "<TH align=right>&nbsp;$numGenes</TH><TH align=left>&nbsp;Genes mapped.</TH>" : "<TH align=right>&nbsp;$numGenes</TH><TH align=left>&nbsp;Gene mapped.</TH>";
			#my $mapString="<TH align=right>&nbsp;$pcMapped</TH><TH align=left>&nbsp;% Proteins mapped to biological resources.".&checkOnGoingMapping."</TH>";
			my ($protString,$mapString)=('','');
			if ($hasValidAna) {
				$protString="<TR><TH align=right><SPAN id=\"protSPAN1\"></SPAN></TH><TH align=left><SPAN id=\"protSPAN2\"><INPUT type=\"button\" class=\"font11\" value=\"Protein summary\" onclick=\"ajaxProteinSummary()\"/></SPAN></TH></TR>";
				$mapString="<TR id=\"protMapTR\" style=\"display:none\"><TH align=right><SPAN id=\"protSPAN3\"></SPAN></TH><TH align=left><SPAN id=\"protSPAN4\"></SPAN>".&checkOnGoingMapping."</TH></TR>";
			}
			$summaryString="<TABLE border=0 cellpadding=0 cellspacing=0><TR>$anaString</TR>$protString$mapString</TABLE>";
		}
		else { # EXPERIMENT
			my ($numSpot,$numSamp,$numAna)=(0,0,0); #$totProt,$visProt,$pcMapped)=(0,0,0);
			my ($numGel2D)=$dbh->selectrow_array("SELECT COUNT(ID_GEL2D) FROM GEL2D WHERE ID_EXPERIMENT=$itemID");
			if ($numGel2D) {
				($numSpot)=$dbh->selectrow_array("SELECT COUNT(ID_SPOT) FROM GEL2D,SPOT WHERE GEL2D.ID_GEL2D=SPOT.ID_GEL2D AND ID_EXPERIMENT=$itemID");
			}
			($numSamp)=$dbh->selectrow_array("SELECT COUNT(ID_SAMPLE) FROM SAMPLE WHERE ID_EXPERIMENT=$itemID");
			if ($numSamp) {
				($numAna,$hasValidAna)=$dbh->selectrow_array("SELECT COUNT(ID_ANALYSIS),MAX(VALID_STATUS) FROM SAMPLE,ANALYSIS WHERE SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND ID_EXPERIMENT=$itemID");
				#if ($numAna) {
				#	($totProt)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM ANALYSIS_PROTEIN,ANALYSIS,SAMPLE WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_EXPERIMENT=$itemID");
				#	if ($totProt) {
				#		($visProt)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM ANALYSIS_PROTEIN,ANALYSIS,SAMPLE WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND VISIBILITY>0 AND ID_EXPERIMENT=$itemID");
				#		my ($numMapped)=$dbh->selectrow_array("SELECT COUNT(DISTINCT P.ID_PROTEIN) FROM PROTEIN P,ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S WHERE ID_EXPERIMENT=$itemID AND P.ID_PROTEIN=AP.ID_PROTEIN AND AP.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND ID_MASTER_PROTEIN IS NOT NULL");
				#		$pcMapped=1*(sprintf "%.1f",($numMapped/$totProt)*100);
				#	}
				#}
			}
			my $gel2dString=($numGel2D>1)? "<TH align=right>&nbsp;$numGel2D</TH><TH align=left>&nbsp;2D Gels</TH>" : "<TH align=right>&nbsp;$numGel2D</TH><TH align=left>&nbsp;2D Gel</TH>";
			my $spotString=($numSpot>1)? "<TH align=right>&nbsp;$numSpot</TH><TH align=left>&nbsp;Spots</TH>" : "<TH align=right>&nbsp;$numSpot</TH><TH align=left>&nbsp;Spot</TH>";
			my $sampString=($numSamp>1)? "<TH align=right>&nbsp;$numSamp</TH><TH align=left>&nbsp;Samples</TH>" : "<TH align=right>&nbsp;$numSamp</TH><TH align=left>&nbsp;Sample</TH>";
			my $pepModifStrg = ($hasValidAna)? '&nbsp;<SPAN id="missCleavSPAN"><INPUT type="button" class="font11" value="% Missed-cleavages" onclick="ajaxMissedCleavage()"/></SPAN>'.&getPtmDistributionButton() : '';
			my $exportAnalysisStrg = ($hasValidAna)? &getExportAnalysisButton() : '';
			my $anaString=($numAna>1)? "<TH align=right valign=top>&nbsp;$numAna</TH><TH align=left>&nbsp;Analyses $pepModifStrg $exportAnalysisStrg</TH>" : "<TH align=right valign=top>&nbsp;$numAna</TH><TH align=left>&nbsp;Analysis $pepModifStrg $exportAnalysisStrg</TH>";
			#my $sStrg=($visProt>1)? 's' : '';
			#my $protString="<TH align=right>&nbsp;$visProt</TH><TH align=left>&nbsp;Protein$sStrg ($totProt max)";
			#$protString.="&nbsp;".&getAmbiguityButton."&nbsp;".&getConflictButton if $totProt;
			#$protString.="&nbsp;</TH>";
			##my $geneString=($numGenes>1)? "<TH align=right>&nbsp;$numGenes</TH><TH align=left>&nbsp;Genes mapped.</TH>" : "<TH align=right>&nbsp;$numGenes</TH><TH align=left>&nbsp;Gene mapped.</TH>";
			#my $mapString="<TH align=right>&nbsp;$pcMapped</TH><TH align=left>&nbsp;% Proteins mapped to biological resources.".&checkOnGoingMapping."</TH>";
			my ($protString,$mapString,$quantiString,$explorString,$funcString)=('','','','','');
			if ($hasValidAna) {
				$protString="<TR><TH align=right><SPAN id=\"protSPAN1\"></SPAN></TH><TH align=left><SPAN id=\"protSPAN2\"><INPUT type=\"button\" class=\"font11\" value=\"Protein summary\" onclick=\"ajaxProteinSummary()\"/></SPAN></TH></TR>";
				$mapString="<TR id=\"protMapTR\" style=\"display:none\"><TH align=right><SPAN id=\"protSPAN3\"></SPAN></TH><TH align=left><SPAN id=\"protSPAN4\"></SPAN>".&checkOnGoingMapping."</TH></TR>";
				##>Child items
				my ($numProtQuant)=$dbh->selectrow_array("SELECT COUNT(ID_QUANTIFICATION) FROM QUANTIFICATION Q,DESIGN D WHERE D.ID_DESIGN=Q.ID_DESIGN AND Q.FOCUS='protein' AND D.ID_EXPERIMENT=$itemID");
				my ($numProtQuant2)=$dbh->selectrow_array("SELECT COUNT(DISTINCT Q.ID_QUANTIFICATION) FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ,ANALYSIS A,SAMPLE S WHERE S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_DESIGN IS NULL AND Q.FOCUS='protein' AND S.ID_EXPERIMENT=$itemID");
				$numProtQuant+=$numProtQuant2; # +=non-design quantif (old)
				$quantiString=($numProtQuant > 1)? "<TR><TH align=right>&nbsp;$numProtQuant</TH><TH align=left>&nbsp;Protein/Isoform quantifications</TH></TR>" : "<TR><TH align=right>&nbsp;$numProtQuant</TH><TH align=left>&nbsp;Protein/Isoform quantification</TH></TR>";
				my ($numExplorAna)=$dbh->selectrow_array("SELECT COUNT(ID_EXPLORANALYSIS) FROM EXPLORANALYSIS WHERE ID_EXPERIMENT=$itemID AND ANA_TYPE NOT LIKE '%MOTIF'");
				$explorString=($numExplorAna > 1)? "<TR><TH align=right>&nbsp;$numExplorAna</TH><TH align=left>&nbsp;Exploratory analyses</TH></TR>" : "<TR><TH align=right>&nbsp;$numExplorAna</TH><TH align=left>&nbsp;Exploratory analysis</TH></TR>";
				my ($numMotifAna)=$dbh->selectrow_array("SELECT COUNT(ID_EXPLORANALYSIS) FROM EXPLORANALYSIS WHERE ID_EXPERIMENT=$itemID AND ANA_TYPE LIKE '%MOTIF'");
				my ($numGOAna)=$dbh->selectrow_array("SELECT COUNT(ID_GOANALYSIS) FROM GO_ANALYSIS WHERE ID_EXPERIMENT=$itemID AND ID_PARENT_GOANA IS NULL");
				my ($numPathAna)=$dbh->selectrow_array("SELECT COUNT(ID_PATHWAY_ANALYSIS) FROM PATHWAY_ANALYSIS WHERE ID_EXPERIMENT=$itemID");
				my $numFuncAna=$numMotifAna+$numGOAna+$numPathAna;
				$funcString=($numFuncAna > 1)? "<TR><TH align=right>&nbsp;$numFuncAna</TH><TH align=left>&nbsp;Functional analyses</TH></TR>" : "<TR><TH align=right>&nbsp;$numFuncAna</TH><TH align=left>&nbsp;Functional analysis</TH></TR>";
			}
			$summaryString="<TABLE border=0 cellpadding=0 cellspacing=0><TR>$gel2dString</TR><TR>$spotString</TR><TR>$sampString</TR><TR>$anaString</TR>$protString$mapString$quantiString$explorString$funcString</TABLE>";
		}
		$name=&promsMod::HTMLcompatible($name);
		push @valueText,"&nbsp;$name";
		$description=&promsMod::HTMLcompatible($description);
		push @valueText,"&nbsp;$description";
		push @valueText,"&nbsp;$displayPos";
		push @valueText,"&nbsp;$startDate";
		push @valueText,"&nbsp;$prefSpeciesStrg" if $item eq 'EXPERIMENT';
		push @valueText,$summaryString;
		push @valueText,$fdrDivString;
		$comments=&promsMod::HTMLcompatible($comments);
		push @valueText,"&nbsp;$comments";
		push @valueText,"&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Show item list\" onclick=\"window.location='./selectAnalyses.cgi?callType=list&ID=$item:$itemID'\"/>";
	}
	else { # add or edit
		if($item eq 'METADATA') {
			my $hasEditRights = ($action eq 'add' || $projectAccess =~ /bioinfo|mass/ || ($accessibility =~ /bio|private/ && $userID eq $metadataOwner)) ? 1 : 0;
			if($hasEditRights) {
				if($action eq 'edit' || !$itemSubType || $itemSubType ne 'ANNOTATION_ITEM') {
					push @valueText,"<input type='text' name='name' value=\"$metaName\" size='50' maxlength='50'>";
					push @valueText, "<textarea name='M.DES' rows='5' cols='65' maxlength='250'>$metaDes</textarea>";
	
						my @accessibilityValues = ('private', 'bio');
						push(@accessibilityValues, ('mass', 'bioinfo')) if($projectAccess eq 'bioinfo');
						push(@accessibilityValues, 'mass') if($projectAccess eq 'mass');
						my $strAccessibilitySelect = "";
						$strAccessibilitySelect .= "<select name='ACCESSIBILITY'>";
						foreach my $accessibilityValue (@accessibilityValues) {
							my $selected = ($accessibility eq $accessibilityValue) ? 'selected' : '';
							$strAccessibilitySelect .= "<option value='$accessibilityValue' $selected>".ucfirst($accessibilityValue)."</option>";
						}
						$strAccessibilitySelect .= "</select>";
						push @valueText, $strAccessibilitySelect;
				}
				
				my @annotTypeValues = ('text', 'file', 'image', 'date', 'link', 'number');
				my $strAnnotTypeSelect = "";
				$strAnnotTypeSelect .= "<select id='ANNOT_TYPE' name='ANNOT_TYPE' onchange='switchFieldType()'>";
				$annotType = 'text' if(!$annotType);
				foreach my $annotTypeValue (@annotTypeValues) {
					my $selected = ($annotType eq $annotTypeValue) ? 'selected' : '';
					$strAnnotTypeSelect .= "<option value='$annotTypeValue' $selected>".ucfirst($annotTypeValue)."</option>";
				}
				$strAnnotTypeSelect .= "</select>";
				
				$strAnnotTypeSelect .=  "<input type='hidden' name='PREVIOUS_ID' value=\'".param("PREVIOUS_ID")."'/>" if(param("PREVIOUS_ID"));
				$strAnnotTypeSelect .=  "<input type='hidden' name='ID' value=\'".param("ID")."'/>" if(param("ID"));
				$strAnnotTypeSelect .=  "<input type='hidden' name='PARENT_ID' value=\'".param("PARENT_ID")."'/>" if(param("PARENT_ID"));
				$strAnnotTypeSelect .=  "<input type='hidden' name='PARENT' value=\'".param("PARENT")."'/>" if(param("PARENT"));
				push @valueText, $strAnnotTypeSelect;
				
				my $displayAnnotValue = $annotValue;
				my $fieldValue = ($annotValue) ? $annotValue : '';
				my $fieldStyle = "";
				
				if($annotType eq 'date') {
					my $dateFormater = Time::Piece->strptime($annotValue, "%Y-%m-%d");
					$displayAnnotValue = $dateFormater->strftime("%d/%m/%Y");
				} elsif($annotType eq 'file' || $annotType eq 'image') {
					$annotType = 'file';
					$displayAnnotValue = basename($annotValue);
					$fieldValue = $displayAnnotValue;
				} elsif($annotType eq 'link' || $annotType eq 'text') {
					$annotType = 'text';
					$fieldStyle = 'width: 541px;';
					$displayAnnotValue = substr($annotValue, 0, ($annotType eq 'link') ? 30 : 50)."...";
				}
				
				my $currentValue = ($action eq 'edit') ? "(current : $displayAnnotValue)" : "";
				push @valueText, "<input type='$annotType' style='$fieldStyle' name='ANNOT_VALUE' id='ANNOT_VALUE' value='$fieldValue' /> <span id='currentValue'>$currentValue</span>";
				push @valueText,"<TEXTAREA name='comments' rows='5' cols='65' maxlength='250'>$comments</TEXTAREA>";
			}
		} else {
			push @valueText,"<INPUT type='text' name='name' value=\"$name\" size='50' maxlength='50'>";
			if ($action eq 'add' && ($item eq 'EXPERIMENT' || $itemInfo[-2]{'ITEM'} eq 'EXPERIMENT')) {
				push @valueText,"<INPUT type='text' name='labels' value='' size='20'/>&nbsp<FONT style=\"font-size:11px\">Use ',' between single values and '-' for range (eg. <B>1,3,5-10</B>).</FONT>"; # labels
			}
			push @valueText,"<TEXTAREA name='des' rows='2' cols='65' maxlength='100'>$description</TEXTAREA>";
			if ($action eq 'edit') {
				my $parentItem=($item eq 'EXPERIMENT')? 'PROJECT' : 'EXPERIMENT';
				my ($maxDisplayPos)=$dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM $item WHERE ID_$parentItem=(SELECT ID_$parentItem FROM $item WHERE ID_$item=$itemID)");
				my $positionStrg='<SELECT name="displayPos">';
				foreach my $pos (1..$maxDisplayPos) {
					$positionStrg.="<OPTION value=\"$pos\"";
					$positionStrg.=' selected' if $pos==$displayPos;
					$positionStrg.=">$pos</OPTION>";
				}
				$positionStrg.="</SELECT><INPUT type='hidden' name='oldDisplayPos' value='$displayPos'/>";
				push @valueText,$positionStrg;
			}
			push @valueText,"&nbsp;$startDate";
			if ($item eq 'EXPERIMENT') {
				my $sthSpecies=$dbh->prepare("SELECT ID_SPECIES,CONCAT(SCIENTIFIC_NAME,' (',COMMON_NAME,')') FROM SPECIES WHERE IS_REFERENCE=1 ORDER BY SCIENTIFIC_NAME");
				$sthSpecies->execute;
				my $speciesSelStrg="<SELECT name=\"prefSpecies\"><OPTION value=\"\">-= Select =-</OPTION>\n";
				while (my ($spID,$spNames)=$sthSpecies->fetchrow_array) {
					$speciesSelStrg.="<OPTION value=\"$spID\"";
					$speciesSelStrg.=' selected' if ($prefSpeciesID && $prefSpeciesID==$spID);
					$speciesSelStrg.=">$spNames</OPTION>\n";
				}
				$sthSpecies->finish;
				$speciesSelStrg.="</SELECT>\n";
				push @valueText,$speciesSelStrg;
			}
			push @valueText,"<TEXTAREA name='comments' rows='5' cols='65' maxlength='250'>$comments</TEXTAREA>";
			
			if ($action eq 'edit' && $item eq 'EXPERIMENT' && $projectAccess =~ /bioinfo|mass/) {
				my %usersLockStatus = ();
				my $lockStatusStr = '';
				my $lockMsg = '';
				
				my $sthUsers=$dbh->prepare("SELECT PA.ID_USER, USER_NAME, LOCK_MSG FROM EXPERIMENT E INNER JOIN PROJECT_ACCESS PA ON PA.ID_PROJECT=E.ID_PROJECT INNER JOIN USER_LIST UL ON UL.ID_USER=PA.ID_USER LEFT JOIN USER_EXPERIMENT_LOCK UEL ON UEL.ID_EXPERIMENT=E.ID_EXPERIMENT AND UEL.ID_USER=UL.ID_USER WHERE E.ID_EXPERIMENT=$itemID");
				$sthUsers->execute;
				while (my ($userID, $userName, $lockMessage)=$sthUsers->fetchrow_array) {
					@{$usersLockStatus{$userID}} = ($userName, $lockMessage);
					$lockMsg = $lockMessage if(!$lockMsg && $lockMessage);
				}
				
				if(scalar keys %usersLockStatus) {
					$lockStatusStr .= qq |
						<table>
							<tr>
								<td>
					|;
					
					foreach $userID (keys %usersLockStatus) {
						my $selected = ($usersLockStatus{$userID}[1]) ? 'checked' : '';
						$lockStatusStr .= "<label><input type='checkbox' name='lockUsers' value='$userID' $selected>".$usersLockStatus{$userID}[0]."</label><br/>"; 
					}
					
					$lockStatusStr .= qq |
								</td>
								<td><TEXTAREA name='lockMsg' rows='4' cols='58' style='max-width:401px' maxlength='250' placeholder='Message to display to locked user(s)'>$lockMsg</TEXTAREA></td>
							</tr>
						</table>
					|;
				} else {
					$lockStatusStr .= "No project users."
				}
				
				push @valueText, $lockStatusStr;
			}
		}
	}
	
	push @valueText, getMetadata($item, $itemID, $action) if ($action ne 'add' && $item ne 'METADATA');
	
	if ($action eq 'edit' && $item ne 'METADATA') {
		if ($item eq 'EXPERIMENT') {
			my @userInfoAll=&promsMod::getUserInfo($dbh,$userID); # $userInfoAll[2] is ref to access info for all projects
			foreach my $projID (keys %{$userInfoAll[2]}) {
				$numParents++ if ${$userInfoAll[2]}{$projID} ne 'guest';
			}
		}
		else { # SAMPLE
			($numParents)=$dbh->selectrow_array("SELECT COUNT(*) FROM EXPERIMENT WHERE ID_PROJECT=$itemInfo[0]{ID}");
		}
	}
}


################
####<2D Gel>####
################
sub gel2D {
	my ($name,$description,$displayPos,$imageFile,$startDate,$comments);
	if ($action ne 'add') { # summary or edit
		@nameText=('Name','Gel image file','Description','Position','Start date');
		push @nameText,('Summary','False discovery') if $action eq 'summary';
		push @nameText,'Comments';
		push @nameText,'List' if $action eq 'summary';
		($name,$description,$imageFile,$displayPos,$startDate,$comments)=$dbh->selectrow_array("SELECT NAME,DES,IMAGE_FILE,DISPLAY_POS,START_DATE,COMMENTS FROM GEL2D WHERE ID_GEL2D=$itemID");
	}
	else { # add
		@nameText=('Name','Gel image file','Description','Start date','Comments'); #,'Spot mapping file'
		$startDate=$date;
	}
	($name,$description,$comments)=&promsMod::chkDef($name,$description,$comments);
	$startDate=&promsMod::formatDate($startDate);
	if ($action eq 'summary') {
		my ($numSamp,$numAna)=(0,0); # $totProt,$visProt,$pcMapped)=(0,0,0);
		my ($numSpot)=$dbh->selectrow_array("SELECT COUNT(*) FROM SPOT WHERE ID_GEL2D=$itemID");
		if ($numSpot) {
			($numSamp)=$dbh->selectrow_array("SELECT COUNT(ID_SAMPLE) FROM SPOT,SAMPLE WHERE SPOT.ID_SPOT=SAMPLE.ID_SPOT AND ID_GEL2D=$itemID");
			if ($numSamp) {
				($numAna,$hasValidAna)=$dbh->selectrow_array("SELECT COUNT(ID_ANALYSIS),MAX(VALID_STATUS) FROM SPOT,SAMPLE,ANALYSIS WHERE SPOT.ID_SPOT=SAMPLE.ID_SPOT AND SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND ID_GEL2D=$itemID");
				#if ($numAna) {
				#	($totProt)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM ANALYSIS_PROTEIN,ANALYSIS,SAMPLE,SPOT WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SPOT.ID_SPOT=SAMPLE.ID_SPOT AND ID_GEL2D=$itemID");
				#	if ($totProt) {
				#		($visProt)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM ANALYSIS_PROTEIN,ANALYSIS,SAMPLE,SPOT WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SPOT.ID_SPOT=SAMPLE.ID_SPOT AND VISIBILITY>0 AND ID_GEL2D=$itemID");
				#		my ($numMapped)=$dbh->selectrow_array("SELECT COUNT(DISTINCT P.ID_PROTEIN) FROM PROTEIN P,ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S,SPOT SP WHERE ID_GEL2D=$itemID AND P.ID_PROTEIN=AP.ID_PROTEIN AND AP.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_SPOT=SP.ID_SPOT AND ID_MASTER_PROTEIN IS NOT NULL");
				#		$pcMapped=1*(sprintf "%.1f",($numMapped/$totProt)*100);
				#	}
				#}
			}
		}
		my $spotString=($numSpot>1)? "<TH align=right>&nbsp;$numSpot</TH><TH align=left>&nbsp;Spots</TH>" : "<TH align=right>&nbsp;$numSpot</TH><TH align=left>&nbsp;Spot</TH>";
		my $sampString=($numSamp>1)? "<TH align=right>&nbsp;$numSamp</TH><TH align=left>&nbsp;Samples</TH>" : "<TH align=right>&nbsp;$numSamp</TH><TH align=left>&nbsp;Sample</TH>";
		my $pepModifStrg = ($hasValidAna)? &getPtmDistributionButton : '';
		my $exportAnalysisStrg = ($hasValidAna)? &getExportAnalysisButton : '';
		my $anaString=($numAna>1)? "<TH align=right valign=top>&nbsp;$numAna</TH><TH align=left>&nbsp;Analyses $pepModifStrg $exportAnalysisStrg</TH>" : "<TH align=right valign=top>&nbsp;$numAna</TH><TH align=left>&nbsp;Analysis $pepModifStrg $exportAnalysisStrg</TH>";
		#my $sStrg=($visProt>1)? 's' : '';
		#my $protString="<TH align=right>&nbsp;$visProt</TH><TH align=left>&nbsp;Protein$sStrg ($totProt max)";
		#$protString.="&nbsp;".&getAmbiguityButton."&nbsp;".&getConflictButton if $totProt;
		#$protString.="</TH>";
		##my $geneString=($numGenes>1)? "<TH align=right>&nbsp;$numGenes</TH><TH align=left>&nbsp;Genes mapped.</TH>" : "<TH align=right>&nbsp;$numGenes</TH><TH align=left>&nbsp;Gene mapped.</TH>";
		#my $mapString="<TH align=right>&nbsp;$pcMapped</TH><TH align=left>&nbsp;% Proteins mapped to biological resources.".&checkOnGoingMapping."</TH>";
		my ($protString,$mapString)=('','');
		if ($hasValidAna) {
			$protString="<TR><TH align=right><SPAN id=\"protSPAN1\"></SPAN></TH><TH align=left><SPAN id=\"protSPAN2\"><INPUT type=\"button\" class=\"font11\" value=\"Protein summary\" onclick=\"ajaxProteinSummary()\"/></SPAN></TH></TR>";
			$mapString="<TR id=\"protMapTR\" style=\"display:none\"><TH align=right><SPAN id=\"protSPAN3\"></SPAN></TH><TH align=left><SPAN id=\"protSPAN4\"></SPAN>".&checkOnGoingMapping."</TH></TR>";
		}
		my $summaryString="<TABLE border=0 cellpadding=0 cellspacing=0><TR>$spotString</TR><TR>$sampString</TR><TR>$anaString</TR>$protString$mapString</TABLE>";

		$name=&promsMod::HTMLcompatible($name);
		push @valueText,"&nbsp;$name";
		$imageFile=&promsMod::HTMLcompatible($imageFile);
		push @valueText,"&nbsp;$imageFile";
		$description=&promsMod::HTMLcompatible($description);
		push @valueText,"&nbsp;$description";
		push @valueText,"&nbsp;$displayPos";
		push @valueText,"&nbsp;$startDate";
		push @valueText,$summaryString;
		push @valueText,$fdrDivString;
		$comments=&promsMod::HTMLcompatible($comments);
		push @valueText,"&nbsp;$comments";
		push @valueText,"&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Show item list\" onclick=\"window.location='./selectAnalyses.cgi?callType=list&ID=$item:$itemID'\"/>";
	}
	else { # add or edit
		push @valueText,"<INPUT type='text' name='name' value='$name' size='50' maxlength='50'>";
		if ($action eq 'add') {
			push @valueText,'<INPUT type="file" name="image_file" size="65"> <I><B>(Jpeg only)</B></I>';
			#push @valueText,'<INPUT type="file" name="spot_file" size="65">';
		}
		else { # edit
			push @valueText,'&nbsp;'.&promsMod::HTMLcompatible($imageFile);
		}
		push @valueText,"<TEXTAREA name='des' rows='2' cols='65' maxlength='100'>$description</textarea>";
		if ($action eq 'edit') {
			my ($maxDisplayPos)=$dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM GEL2D WHERE ID_EXPERIMENT=(SELECT ID_EXPERIMENT FROM GEL2D WHERE ID_GEL2D=$itemID)");
			my $positionStrg='<SELECT name="displayPos">';
			foreach my $pos (1..$maxDisplayPos) {
				$positionStrg.="<OPTION value=\"$pos\"";
				$positionStrg.=' selected' if $pos==$displayPos;
				$positionStrg.=">$pos</OPTION>";
			}
			$positionStrg.="</SELECT><INPUT type='hidden' name='oldDisplayPos' value='$displayPos'/>";
			push @valueText,$positionStrg;
		}
		push @valueText,"&nbsp;$startDate";
		push @valueText,"<TEXTAREA name='comments' rows='5' cols='65' maxlength='250'>$comments</textarea>";
	}
	if ($action eq 'edit') {
		($numParents)=$dbh->selectrow_array("SELECT COUNT(*) FROM EXPERIMENT WHERE ID_PROJECT=$itemInfo[0]{ID}");
	}
}


##############
## METADATA ##
##############
sub getMetadata {
	my ($item, $itemID, $action) = @_;
	my	$strMetaData .= "<ul style='padding: 0px 0px 0px 6px; position: relative; bottom: 12px; margin: 15px 0px -6px 0px;' class='js-sortable' aria-grabbed='false' aria-dropeffect='move'>\n";
	my $idMeta = -1;
	my %itemIcones=&promsConfig::getItemIcones;
	my $draggableOptions = ($action eq 'edit') ? "draggable='true' aria-grabbed='false' aria-dropeffect='move'" : '';
	
	# Filter by item type
	my ($projectID, $experimentID, $sampleID) = $dbh->selectrow_array("SELECT PROJECT.ID_PROJECT, EXPERIMENT.ID_EXPERIMENT, SAMPLE.ID_SAMPLE FROM PROJECT LEFT JOIN EXPERIMENT ON EXPERIMENT.ID_PROJECT = PROJECT.ID_PROJECT LEFT JOIN SAMPLE ON SAMPLE.ID_EXPERIMENT = EXPERIMENT.ID_EXPERIMENT WHERE $item.ID_$item=$itemID LIMIT 1");
	my $sampleIDFilter = ($item ne 'SAMPLE') ? "IS NULL" : "=$sampleID";
	my $experimentIDFilter = ($item eq 'PROJECT') ? "IS NULL" : "=$experimentID";
	my $itemFilter = "ID_PROJECT=$projectID AND ID_EXPERIMENT $experimentIDFilter AND ID_SAMPLE $sampleIDFilter";
	
	# Filter by rights
	my $rightsFilter = "";
	if($projectAccess !~ /bioinfo|manag/) {
		$rightsFilter .= " AND (";
		$rightsFilter .= "(ACCESSIBILITY='private' AND RECORD_USER='$userID') OR ACCESSIBILITY='bio'";
		$rightsFilter .= " OR ACCESSIBILITY='mass'" if($projectAccess eq 'mass');
		$rightsFilter .= ")";
	}

	my $sthMetaData=$dbh->prepare("SELECT M.ID_META_ANNOTATION, ID_ANNOTATION_ITEM, ID_PROJECT, ID_EXPERIMENT, ID_SAMPLE, NAME, A.DES as DES_ANNOT, M.DES as META_DES, M.DISPLAY_POS as META_DISPLAY_POS, A.DISPLAY_POS as ANNOT_DISPLAY_POS, ANNOT_TYPE, ANNOT_VALUE, ACCESSIBILITY, RECORD_USER, RECORD_DATE FROM META_ANNOTATION M INNER JOIN ANNOTATION_ITEM A ON A.ID_META_ANNOTATION = M.ID_META_ANNOTATION WHERE $itemFilter $rightsFilter ORDER BY M.DISPLAY_POS, A.DISPLAY_POS");
	$sthMetaData->execute;
	my $previousMetaHasName = 0;
	while (my ($idMetaItem, $idAnnotationItem, $idProject, $idExperiment, $idSample, $metaName, $annotDes, $metaDes, $annotDisplayPos, $metaDisplayPost, $type, $value, $accessibility, $recordUser, $recordDate) = $sthMetaData->fetchrow_array) {
		my $hasEditRights = ($projectAccess =~ /bioinfo|mass/ || ($accessibility =~ /bio|private/ && $userID eq $recordUser)) ? 1 : 0;
		my $disabled = ($hasEditRights) ? '' : 'disabled';
		my $classMetaName = '';
		my $mayPrintMetaId = '';
		
		if($idMeta != $idMetaItem) {
			$strMetaData .= "</ul>\n" if($previousMetaHasName);
			
			if($metaName) {
				$strMetaData .= "<li data-id-metadata='$idMetaItem' class='$disabled meta draggable-line meta-annotation-row' style='list-style: none;' $draggableOptions>\n";
				if($action eq 'edit') {
					$strMetaData .= qq |
							<input type='hidden' name='metaAnnot' value='$idMetaItem' />\n
							<div class='actionsDiv'>\n
					|;
					
					if($hasEditRights) {
						$strMetaData .= qq |
									<button type="button" class='action' style='cursor: grab'><i class="fas fa-arrows-alt"></i></button>\n
									<!-- <button type="button" class='action' onclick='window.location="./editProjectItem.cgi?ACT=edit&ITEM=METADATA&ID=$idMetaItem&TYPE=METADATA&PARENT=$item&PARENT_ID=$itemID";'><i class="fas fa-edit"></i></button></i>\n -->
									<button type="button" class='action' onclick='window.location="./editProjectItem.cgi?ACT=add&ITEM=METADATA&ID=$itemID&TYPE=ANNOTATION_ITEM&PARENT=$item&PREVIOUS_ID=$idMetaItem";'><i class="fas fa-plus"></i></button>\n
									<!-- <button type="button" class='action' onclick='if(confirm("Do you really want to delete this metadata ? (Save is also required after confirmation)")) { delMetadataRow(this.parentNode.parentNode) }'><i class="fas fa-minus"></i></button>\n -->
						|;
					}
					
					$strMetaData .= "</div>\n";
				}
				$strMetaData .= qq |
						<div>\n
							<span class='name'>$metaName</span> &nbsp;<span class='desc'>$metaDes</span>\n
						</div>\n
						<ul style='padding: 0px 6px 8px 6px;' class='js-sortable-inner meta-annotation-list' $draggableOptions>\n
				|;
				
				$previousMetaHasName = 1;
			} else {
				$classMetaName = 'annotMetaLike';
				$mayPrintMetaId = "data-id-metadata='$idMetaItem'";
				$previousMetaHasName = 0;
			}
		}
		$strMetaData .= "<li $mayPrintMetaId data-id-annotation='$idAnnotationItem' class='$disabled $classMetaName annot draggable-line annotation-item-row' style='list-style: none;' $draggableOptions>\n";
		
		if($action eq 'edit') {
			$strMetaData .= "<input type='hidden' name='metaAnnot' value='$idMetaItem' />\n" if(!$metaName);
			$strMetaData .= qq |
								<input type="hidden" name="annotationItem" value="$idAnnotationItem" />\n
								<div class='actionsDiv'>\n
			|;
			
			if($hasEditRights) {
				$strMetaData .= qq |
									<button type="button" class='action' onclick='return false;' style='cursor: grab' formnovalidate><i class="fas fa-arrows-alt"></i></button>\n
									<button type="button" class='action' onclick='window.location="./editProjectItem.cgi?ACT=edit&ITEM=METADATA&ID=$idAnnotationItem&TYPE=ANNOTATION_ITEM&PARENT=$item&PARENT_ID=$itemID";'><i class="fas fa-edit"></i></button></i>\n
									<!-- <button type="button" class='action' onclick='window.location="./editProjectItem.cgi?ACT=add&ITEM=METADATA&ID=$itemID&TYPE=ANNOTATION_ITEM&PARENT=$item&PREVIOUS_ID=$idAnnotationItem";'><i class="fas fa-plus"></i></button>\n -->
									<button type="button" class='action' onclick='if(confirm("Do you really want to delete this metadata ? (Save is also required after confirmation)")) { delMetadataRow(this.parentNode.parentNode) }'><i class="fas fa-minus"></i></button>\n
				|;
			}
			
			$strMetaData .= "</div>\n";
		}
		
		$strMetaData .= "<div>\n";
		
		if($type eq 'file' || $type eq 'image') {
			my $linkFeature = ($type eq 'file') ? "download" : "target='_blank'";
			my $filePath = $promsPath{"metadata_html"}."/proj_$projectID/";
			$filePath .= "exp_$experimentID/" if($experimentID && $item ne 'PROJECT');
			$filePath .= "samp_$sampleID/" if($sampleID && $item eq 'SAMPLE');
			$filePath .= "$value";
			
			$strMetaData .= "<span class='value' style='text-align:left'>\n";
			$strMetaData .= "<a href='$filePath' $linkFeature>".basename($value)." <img class='imgFile' src='$promsPath{images}/$itemIcones{file}'/></a>\n";
			$strMetaData .= "</span><span class='desc'>$annotDes</span>\n";
		} elsif($type eq 'date') {
			my $dateFormater = Time::Piece->strptime($value, "%Y-%m-%d");
			$value = $dateFormater->strftime("%d/%m/%Y");
			$strMetaData .= "<span class='value' style='text-align:left'>\n";
			$strMetaData .= "<b>$value</b>\n";
			$strMetaData .= "</span><span class='desc'>$annotDes</span>\n";
		} elsif($type eq 'number') {
			$strMetaData .= "<span class='value'>\n";
			$strMetaData .= "<b>$annotDes : $value</b>\n";
		} elsif($type eq 'link') {
			my $linkValue = substr($value, 0, 20).'...';
			$strMetaData .= "<span class='value' style='text-align:left'>\n";
			$strMetaData .= "<b><a href ='$value' target='_blank'>$linkValue</a></b>\n";
			$strMetaData .= "</span><span class='desc'>$annotDes</span>\n";
		} else {
			$strMetaData .= "<span class='value' style='text-align:left'>\n";
			$strMetaData .= "<b>$value</b>\n";
			$strMetaData .= "</span><span class='desc'>$annotDes</span>\n";
		}
		
		$strMetaData .= "	</div>\n";
		$strMetaData .= "</li>\n";
		
		$idMeta = $idMetaItem;
	}
	$sthMetaData->finish;
	
	$strMetaData .= "</ul></ul>";
	$strMetaData .= qq | <button type="button" class='action' id='addMetaButton' onclick='window.location="./editProjectItem.cgi?ACT=add&ITEM=METADATA&ID=$itemID&TYPE=META_ANNOTATION&PARENT=$item";'><i class="fas fa-plus"></i></button> | if($action eq 'edit');
	
	return $strMetaData;
}

sub updateMetadata {
	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
	
	my $metaType = &promsMod::cleanParameters(param("metaType"));
	my $metaID = &promsMod::cleanParameters(param("metaID"));
	my $dbh = &promsConfig::dbConnect;
	my @positions = split(/,/, param("positions")) if(param("positions"));
	
	if(@positions) {
		if($metaType eq "sub") {
			foreach my $position (@positions) {
				my ($annotID, $annotPos) = split(/:/, $position);
				my $affected = $dbh->do("UPDATE ANNOTATION_ITEM SET DISPLAY_POS=$annotPos WHERE ID_ANNOTATION_ITEM=$annotID");
				print("Affected rows : $affected");
			}
		} elsif($metaType eq "global") {
			foreach my $position (@positions) {
				my ($annotID, $annotPos) = split(/:/, $position);
				my $affected = $dbh->do("UPDATE META_ANNOTATION SET DISPLAY_POS=$annotPos WHERE ID_META_ANNOTATION=$annotID");
				print("Affected rows : $affected");
			}
		}
	}
	
	$dbh->commit;
	$dbh->disconnect;
	exit;
}

##############
####<Spot>####
##############
sub spot {
	my ($name,$description,$coordString,$isoPoint,$molWeight,$intensity,$externalID,$comments);
	@nameText=('Name','Description','Coordinates','Isoelectric point','Molecular weight','Intensity','External identifier','Associated Sample');
	push @nameText,('Summary','False discovery') if $action eq 'summary';
	push @nameText,'Comments';
	push @nameText,'List' if $action eq 'summary';
	($name,$description,my $xPos,my $yPos,$isoPoint,$molWeight,$intensity,$externalID,$comments)=$dbh->selectrow_array("SELECT NAME,DES,X_POS,Y_POS,ISOELECTRIC_POINT,MOLECULAR_WEIGHT,INTENSITY,EXTERNAL_ID,COMMENTS FROM SPOT WHERE ID_SPOT=$itemID");
	$coordString="<B>X:</B>$xPos&nbsp;&nbsp;&nbsp;<B>Y:</B>$yPos";

	##>Associated sample
	my ($sampID,$sampName)=$dbh->selectrow_array("SELECT ID_SAMPLE,NAME FROM SAMPLE WHERE ID_SPOT=$itemID");
	($sampID,$sampName)=(0,'None') unless $sampID;

	($name,$description,$comments)=&promsMod::chkDef($name,$description,$comments);
	if ($action eq 'summary') {
		$isoPoint='Unknown' unless $isoPoint;
		$molWeight='Unknown' unless $molWeight;
		$intensity='Unknown' unless defined $intensity;
		$externalID='Unknown' unless $externalID;
		my $numAna=0; #$totProt,$visProt,$pcMapped)=(0,0,0);
		#my ($numSamp)=$dbh->selectrow_array("SELECT COUNT(ID_SAMPLE) FROM SAMPLE WHERE ID_SPOT=$itemID");
		if ($sampID) {
			($numAna,$hasValidAna)=$dbh->selectrow_array("SELECT COUNT(ID_ANALYSIS),MAX(VALID_STATUS) FROM SAMPLE,ANALYSIS WHERE SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND ID_SPOT=$itemID");
			#if ($numAna) {
			#	($totProt)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM ANALYSIS_PROTEIN,ANALYSIS,SAMPLE WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_SPOT=$itemID");
			#	if ($totProt) {
			#		($visProt)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM ANALYSIS_PROTEIN,ANALYSIS,SAMPLE WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND VISIBILITY>0 AND ID_SPOT=$itemID");
			#		my ($numMapped)=$dbh->selectrow_array("SELECT COUNT(DISTINCT P.ID_PROTEIN) FROM PROTEIN P,ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S WHERE ID_SPOT=$itemID AND P.ID_PROTEIN=AP.ID_PROTEIN AND AP.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND ID_MASTER_PROTEIN IS NOT NULL");
			#		$pcMapped=1*(sprintf "%.1f",($numMapped/$totProt)*100);
			#	}
			#}
		}
		#my $sampString=($numSamp>1)? "<TH align=right>&nbsp;$numSamp</TH><TH align=left>&nbsp;Samples</TH>" : "<TH align=right>&nbsp;$numSamp</TH><TH align=left>&nbsp;Sample</TH>";
		my $pepModifStrg = ($hasValidAna)? &getPtmDistributionButton : '';
		my $exportAnalysisStrg = ($hasValidAna)? &getExportAnalysisButton : '';
		my $anaString = ($numAna>1)? "<TH align=right valign=top>&nbsp;$numAna</TH><TH align=left>&nbsp;Analyses $pepModifStrg $exportAnalysisStrg</TH>" : "<TH align=right valign=top>&nbsp;$numAna</TH><TH align=left>&nbsp;Analysis $pepModifStrg $exportAnalysisStrg</TH>";
		#my $sStrg=($visProt>1)? 's' : '';
		#my $protString="<TH align=right>&nbsp;$visProt</TH><TH align=left>&nbsp;Protein$sStrg ($totProt max)";
		#$protString.="&nbsp;".&getAmbiguityButton."&nbsp;".&getConflictButton if $totProt;
		#$protString.="&nbsp;</TH>";
		##my $geneString=($numGenes>1)? "<TH align=right>&nbsp;$numGenes</TH><TH align=left>&nbsp;Genes mapped.</TH>" : "<TH align=right>&nbsp;$numGenes</TH><TH align=left>&nbsp;Gene mapped.</TH>";
		#my $mapString="<TH align=right>&nbsp;$pcMapped</TH><TH align=left>&nbsp;% Proteins mapped to biological resources.".&checkOnGoingMapping."</TH>";
		my ($protString,$mapString)=('','');
		if ($hasValidAna) {
			$protString="<TR><TH align=right><SPAN id=\"protSPAN1\"></SPAN></TH><TH align=left><SPAN id=\"protSPAN2\"><INPUT type=\"button\" class=\"font11\" value=\"Protein summary\" onclick=\"ajaxProteinSummary()\"/></SPAN></TH></TR>";
			$mapString="<TR id=\"protMapTR\" style=\"display:none\"><TH align=right><SPAN id=\"protSPAN3\"></SPAN></TH><TH align=left><SPAN id=\"protSPAN4\"></SPAN>".&checkOnGoingMapping."</TH></TR>";
		}
		my $summaryString="<TABLE border=0 cellpadding=0 cellspacing=0><TR>$anaString</TR>$protString$mapString</TABLE>"; #<TR>$sampString</TR>
		$name=&promsMod::HTMLcompatible($name);
		push @valueText,"&nbsp;$name";
		$description=&promsMod::HTMLcompatible($description);
		push @valueText,"&nbsp;$description";
		push @valueText,"&nbsp;$coordString";
		push @valueText,"&nbsp;$isoPoint";
		push @valueText,"&nbsp;$molWeight kDa";
		push @valueText,"&nbsp;$intensity";
		push @valueText,"&nbsp;$externalID";
		push @valueText,"&nbsp;$sampName";
		push @valueText,$summaryString;
		push @valueText,$fdrDivString;
		$comments=&promsMod::HTMLcompatible($comments);
		push @valueText,"&nbsp;$comments";
		push @valueText,"&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Show item list\" onclick=\"window.location='./selectAnalyses.cgi?callType=list&ID=$item:$itemID'\"/>";
	}
	else { # edit (cannot be add)
		($isoPoint,$molWeight,$intensity,$externalID)=&promsMod::chkDef($isoPoint,$molWeight,$intensity,$externalID);
		push @valueText,"<INPUT type='text' name='name' value='$name' size='50' maxlength='50'>";
		push @valueText,"<TEXTAREA name='des' rows='2' cols='65' maxlength='100'>$description</textarea>";
		push @valueText,"&nbsp;$coordString";
		push @valueText,"<INPUT type='text' name='isoPoint' value='$isoPoint' size='5'>";
		push @valueText,"<INPUT type='text' name='molWeight' value='$molWeight' size='5'> kDa";
		push @valueText,"<INPUT type='text' name='intensity' value='$intensity' size='5'>";
		push @valueText,"<INPUT type='text' name='externalID' value='$externalID' size='25'>";
		##>List of candidate samples<##
		my $sampleString="<TABLE cellpadding=0 cellspacing=0><TR><TD><SELECT name=\"assoSampleID\" onchange=\"displayNewSample(this.value)\"><OPTION value=\"0\">-= None =-</OPTION>\n";
		my $sthCS=$dbh->prepare("SELECT ID_SAMPLE,NAME FROM SAMPLE WHERE ID_EXPERIMENT=$itemInfo[1]{ID} AND (ID_SPOT IS NULL OR ID_SPOT=$itemID) ORDER BY DISPLAY_POS ASC");
		$sthCS->execute;
		while (my ($sID,$sName)=$sthCS->fetchrow_array) {
			my $selStrg=($sID==$sampID)? ' selected' : '';
			$sampleString.="<OPTION value=\"$sID\"$selStrg>$sName</OPTION>\n";
		}
		$sampleString.='<OPTION value="-1">-= New =-</OPTION></SELECT></TD><TD id="newSample" style="visibility:hidden">&nbsp;<B>Name:</B><INPUT type="text" name="newSampleName" size=20 value=""/></TD></TR></TABLE>';
		$sampleString.="\n<INPUT type='hidden' name='oldSampleID' value='$sampID'/>";
		push @valueText,$sampleString;
		push @valueText,"<TEXTAREA name='comments' rows='5' cols='65' maxlength='250'>$comments</textarea>";
	}
	#($numParents)=$dbh->selectrow_array("SELECT COUNT(*) FROM GEL2D WHERE ID_EXPERIMENT=$itemInfo[1]{ID}");
}


##################
####<Analysis>####
##################
sub analysis {
	my ($name,$description,$displayPos,$startDate,$msType,$instrument,$dataFile,$wiffFile,$labCode,$FDRdata,$taxonomy,$validUserID,$validDate,$comments,$verifMG,$labelMethod,$maxRankImport,$minScore,$lowerScores);
	my ($dataBankString,$statusString,$pepQuantifStrg,$bioSampleStrg,$searchEngine,$minScoreString,$fdrString,$validUser,$refRetTimeString);
	my $stringMG='';
	my $scanDBString='<SELECT name="scanDB"><OPTION value="auto" selected>auto (* recommanded *)</OPTION><OPTION value="now">now</OPTION></SELECT>';
	#if ($action ne 'add') { # summary or edit (no longer add)
	@nameText=('Name','Laboratory code','Description','Position','Start date','MS type','Search engine','Search file','MS data file','Databank(s)','Taxonomy','Reference RT','Labeling','Max. rank','Min. score','False discovery');
	my $anaQuery=qq |SELECT NAME,DES,DISPLAY_POS,START_DATE,MS_TYPE,INSTRUMENT,DATA_FILE,FILE_FORMAT,WIFF_FILE,LAB_CODE,DECOY,FALSE_DISCOVERY,TAXONOMY,VALID_STATUS,VALID_USER,VALID_DATE,COMMENTS,VERIFIED_MG,MIN_SCORE,MAX_RANK,LABELING,LOWER_SCORES
					FROM ANALYSIS WHERE ID_ANALYSIS=$itemID|;
	($name,$description,$displayPos,$startDate,$msType,$instrument,$dataFile,$fileFormat,$wiffFile,$labCode,$decoy,$FDRdata,$taxonomy,$validStatus,$validUserID,$validDate,$comments,$verifMG,$minScore,$maxRankImport,$labelMethod,$lowerScores)=$dbh->selectrow_array($anaQuery);
	$lowerScores=0 unless $lowerScores;
	$taxonomy='Unknown' unless $taxonomy;
	if ($validStatus>=1) {
		#push @nameText,'Match groups';
		#push @nameText,'Validation date' if $validStatus>=2;
		push @nameText,('Match groups','Status','Peptide quantification(s)','Biological sample(s)','Validation date');
		my (%referenceRT,%pepQuantifName);
		my $maxSrcRank=0;
		#my $sthQRT=$dbh->prepare("SELECT RT.ID_REFERENCE_RT,RT.NAME,Q.ID_QUANTIFICATION,Q.NAME,QR.SOURCE_RANK,QR.NUM_PEP,CORRELATION FROM REFERENCE_RT RT,ANALYSIS_REFRT A,QUANTIF_REFRT QR,QUANTIFICATION Q,ANA_QUANTIFICATION AQ WHERE RT.ID_REFERENCE_RT=A.ID_REFERENCE_RT AND RT.ID_REFERENCE_RT=QR.ID_REFERENCE_RT AND QR.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND A.ID_ANALYSIS=$itemID AND AQ.ID_ANALYSIS=$itemID");
		my $sthART=$dbh->prepare("SELECT RT.ID_REFERENCE_RT,RT.NAME,RT.NUM_PEP,AR.NUM_PEP FROM REFERENCE_RT RT,ANALYSIS_REFRT AR WHERE RT.ID_REFERENCE_RT=AR.ID_REFERENCE_RT AND AR.ID_ANALYSIS=$itemID");
		my $sthQRT=$dbh->prepare("SELECT Q.ID_QUANTIFICATION,Q.NAME,QR.SOURCE_RANK,QR.NUM_PEP,CORRELATION FROM QUANTIF_REFRT QR,QUANTIFICATION Q,ANA_QUANTIFICATION AQ WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION=QR.ID_QUANTIFICATION AND QR.ID_REFERENCE_RT=? AND AQ.ID_ANALYSIS=$itemID");
		$sthART->execute;
		while (my ($refRT_ID,$refRTname,$refNumPep,$AnaNumPep)=$sthART->fetchrow_array) {
			$referenceRT{$refRT_ID}{NAME}=$refRTname;
			$referenceRT{$refRT_ID}{REF_NUM_PEP}=$refNumPep;
			$referenceRT{$refRT_ID}{ANA_NUM_PEP}=$AnaNumPep;

			$sthQRT->execute($refRT_ID);
			while (my ($quantifID,$quantifName,$srcRank,$numPep,$r2)=$sthQRT->fetchrow_array) {
				$maxSrcRank=$srcRank if $srcRank > $maxSrcRank;
				if ($referenceRT{$refRT_ID}{'DATA'} && $referenceRT{$refRT_ID}{'DATA'}{$quantifID}) {
					$referenceRT{$refRT_ID}{'DATA'}{$quantifID}[0]=$numPep if $numPep < $referenceRT{$refRT_ID}{'DATA'}{$quantifID}[0];
					$referenceRT{$refRT_ID}{'DATA'}{$quantifID}[1]=$numPep if $numPep > $referenceRT{$refRT_ID}{'DATA'}{$quantifID}[1];
					$referenceRT{$refRT_ID}{'DATA'}{$quantifID}[2]=$r2 if $r2 < $referenceRT{$refRT_ID}{'DATA'}{$quantifID}[2];
					$referenceRT{$refRT_ID}{'DATA'}{$quantifID}[3]=$r2 if $r2 > $referenceRT{$refRT_ID}{'DATA'}{$quantifID}[3];
				}
				else {
					@{$referenceRT{$refRT_ID}{'DATA'}{$quantifID}}=($numPep,$numPep,$r2,$r2); #min/max
				}
				$pepQuantifName{$quantifID}=$quantifName;
			}
		}
		$sthART->finish;
		$sthQRT->finish;

		my ($numRef,$numAlignProc)=(0,0);
		foreach my $refRT_ID (sort{lc($referenceRT{$a}{'NAME'}) cmp lc($referenceRT{$b}{'NAME'})} keys %referenceRT) {
			$numRef++;
			$refRetTimeString.='<BR>&nbsp;' if $refRetTimeString;
			$refRetTimeString.="<B>$referenceRT{$refRT_ID}{NAME}</B> ($referenceRT{$refRT_ID}{ANA_NUM_PEP}/$referenceRT{$refRT_ID}{REF_NUM_PEP} peptides validated)";
			my $numQuantifs=scalar keys %{$referenceRT{$refRT_ID}{'DATA'}};
			if ($numQuantifs) {
				$numAlignProc++;
				foreach my $quantifID (sort{$a<=>$b} keys %{$referenceRT{$refRT_ID}{'DATA'}}) {
					$refRetTimeString.=":<BR> w/ <B>".&promsMod::resize($pepQuantifName{$quantifID},15)."</B>" if $numQuantifs > 1;
					$refRetTimeString.=": R<SUP>2</SUP>=";
					my $r2min=sprintf "%.3f",$referenceRT{$refRT_ID}{DATA}{$quantifID}[2];
					if ($maxSrcRank > 0) { #merged search
						my $r2max=sprintf "%.3f",$referenceRT{$refRT_ID}{DATA}{$quantifID}[3];
						$refRetTimeString.="$r2min~$r2max based on $referenceRT{$refRT_ID}{DATA}{$quantifID}[0]~$referenceRT{$refRT_ID}{DATA}{$quantifID}[1] peptides ($maxSrcRank regressions).";
					}
					else {
						$refRetTimeString.="$r2min based on $referenceRT{$refRT_ID}{DATA}{$quantifID}[0] peptides.";
					}
					push @regressionPlots,[$referenceRT{$refRT_ID}{NAME},$pepQuantifName{$quantifID},$refRT_ID,$quantifID];
				}
			}
		}
		if ($numRef) {
			if ($numAlignProc) {
				my $plotStrg=($numAlignProc > 1)? 'plots' : 'plot';
				$refRetTimeString.="&nbsp;<INPUT type=\"button\" id=\"regressionBUTTON\" class=\"font11\" value=\"Show regression $plotStrg\" onclick=\"showHideRegressionPlots($numAlignProc)\"/>";
			}
			else {$refRetTimeString.=': No alignment performed.';}
		}
		else {$refRetTimeString='None detected';}
	}
	else {
		push @nameText,'Status';
		$refRetTimeString='N/A (Data not reported)';
	}
	push @nameText,'Comments';

	#$fileName = $dataFile;
	#if ($QuantifID) {$quantiMethod = $dbh->selectrow_array("SELECT NAME FROM QUANTIF_METHOD WHERE ID_QUANTIF_METHOD=$QuantifID");}
	#else {
	#	$QuantifID=0;
	#	$quantiMethod='None';
	#}
	$labelMethod='None' unless $labelMethod;
	$minScoreString=($decoy && $decoy=~/precomputed/)? 'N/A' : (!defined $minScore)? '-' : $minScore; # could be =0
	$maxRankImport= "-" unless $maxRankImport ;
	$metaAna=1 if ($fileFormat=~/MASCOT/ && $dataFile=~/^C/);
	if ($validUserID) {
		($validUser)=$dbh->selectrow_array("SELECT USER_NAME FROM USER_LIST WHERE ID_USER='$validUserID'");
		$validUser=$validUserID unless $validUser;
	}
	$searchEngine=($fileFormat =~ /MASCOT/)? 'Mascot' : ($fileFormat =~ /PHENYX/)? 'Phenyx' : ($fileFormat =~ /SEQUEST/)? 'Sequest' : ($fileFormat =~ /MAXQUANT/)? 'Andromeda from MaxQuant' : 'Unknown';
	$searchEngine.=' from Proteome Discoverer' if $fileFormat =~ /\.PDM/;

	###>Fetching databank info
	#my ($databankName,$versionName,$databankType)=$dbh->selectrow_array("SELECT DATABANK.NAME,VERSION_NAME,DATABANK_TYPE.NAME FROM DATABANK,DATABANK_TYPE WHERE DATABANK.ID_DBTYPE=DATABANK_TYPE.ID_DBTYPE AND ID_DATABANK=$dataBankID");
	#$dataBankString="$databankName";
	#$dataBankString.=" ($versionName)" if $versionName;
	#$dataBankString.="&nbsp;&nbsp;&nbsp<B>Type:</B> $databankType";
	my $sthDB=$dbh->prepare("SELECT D.NAME,D.VERSION_NAME,DT.NAME FROM ANALYSIS_DATABANK AD,DATABANK D,DATABANK_TYPE DT WHERE AD.ID_ANALYSIS=$itemID AND AD.ID_DATABANK=D.ID_DATABANK AND D.ID_DBTYPE=DT.ID_DBTYPE ORDER BY DB_RANK");
	$sthDB->execute;
	while (my ($dbName,$dbVersion,$dbType)=$sthDB->fetchrow_array) {
		$dataBankString.='<BR>&nbsp;' if $dataBankString;
		$dataBankString.="&bull;$dbName";
		$dataBankString.=" ($dbVersion)" if $dbVersion;
		$dataBankString.="&nbsp;&nbsp;&nbsp<B>Type:</B> $dbType";
	}
	$sthDB->finish;

	###>Validation history button
	my $historyButtonStrg = '';
	my ($valHistory) = $dbh->selectrow_array("SELECT COUNT(*) FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$itemID");
	if ($action eq 'summary'){
		my $existValHistory;
		if ($projectAccess =~ /bioinfo|mass|manag|super/) {
			$existValHistory=$valHistory;
		}
		else {
			($existValHistory) = $dbh->selectrow_array("SELECT COUNT(*) FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$itemID AND (STATUS=1 OR STATUS=-1)");
		}
		#$historyButtonStrg = "<INPUT type=button  class=\"font11\" value=\"Show validation history\" onclick=\"window.location='./displayValidationHistory.cgi?ID=$rowID'\"";
		$historyButtonStrg = "<INPUT type=button  class=\"font11\" value=\"Show validation history\" onclick=\"displayValidationHistory()\"";
		if(!$existValHistory){ $historyButtonStrg .= ' disabled' };
		$historyButtonStrg .= '>';
	}

	###>Lower-scoring peptides
	my $lowerScoresString='';
	if ($lowerScores) {
		$lowerScoresString=($lowerScores==2)? '<BR>&nbsp;All' : '&nbsp;<BR>Some';
		$lowerScoresString.=' lower-scoring peptides activated.';
	}

	###>Validation status
	#my $projectID; # only if (partially) validated
	if ($validStatus>=1) { # (partially) validated
		my $sthcountProt=$dbh->prepare("SELECT VISIBILITY,CONF_LEVEL,COUNT(ID_PROTEIN) FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=$itemID GROUP BY VISIBILITY,CONF_LEVEL");
		$sthcountProt->execute;
		my ($totProt,$visProt,$beforeMCQtotProt,$beforeMCQvisProt)=(0,0,0,0);
		while (my ($visibility,$confLevel,$count)=$sthcountProt->fetchrow_array){
			$totProt+=$count;
			if ($visibility){
				$beforeMCQvisProt=$count if $confLevel==2;
				$visProt+=$count;
			}
			$beforeMCQtotProt+=$count if $confLevel==2;
		}

		#my ($beforeMCQProt)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM ANALYSIS_PROTEIN WHERE VISIBILITY>0 AND CONF_LEVEL=2 AND ID_ANALYSIS=$itemID");
		#my ($visProt)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM ANALYSIS_PROTEIN WHERE VISIBILITY>0 AND ID_ANALYSIS=$itemID");
		#my ($totProt)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=$itemID");
		#$projectID=&promsMod::getProjectID($dbh,$itemID,'ANALYSIS');
		my ($projStatus)=$dbh->selectrow_array("SELECT STATUS FROM PROJECT WHERE ID_PROJECT=$projectID");
		$projStatus=0 unless $projStatus;
		#push @nameText,'Validation date' if $validStatus>=2;
		#push @nameText,'Match groups';
		if ($visProt) {
			my ($numMapped)=$dbh->selectrow_array("SELECT COUNT(DISTINCT P.ID_PROTEIN) FROM PROTEIN P,ANALYSIS_PROTEIN AP WHERE ID_ANALYSIS=$itemID AND P.ID_PROTEIN=AP.ID_PROTEIN AND ID_MASTER_PROTEIN IS NOT NULL");
			my $pcMapped=1*(sprintf "%.1f",($numMapped/$totProt)*100);
			my $protString=($visProt>1)? 'Proteins' : 'Protein';
			#my $geneString=($numGenes>1)? "$numGenes Genes mapped" : "$numGenes Gene mapped";
			my $mapString="<TH align=right valign=top>&nbsp;$pcMapped</TH><TH align=left>&nbsp;% Proteins mapped to biological resources.".&checkOnGoingMapping."</TH>";
			$statusString=($validStatus==1)? '<B>Partially v' : '<B>V';
			$statusString.="alidated. $historyButtonStrg$lowerScoresString";
			$statusString.='&nbsp;<SPAN id="missCleavSPAN"><INPUT type="button" class="font11" value="% Missed-cleavages" onclick="ajaxMissedCleavage()"/></SPAN>&nbsp;'.&getPtmDistributionButton."&nbsp;" if $totProt ;
			my $visProtString=($visProt != $beforeMCQvisProt)? "$beforeMCQvisProt~$visProt": $visProt;
			my $totProtString=($totProt != $beforeMCQtotProt)? "$beforeMCQtotProt~$totProt": $totProt;
			$statusString.="<BR><TABLE border=0 cellpadding=0 cellspacing=0><TR><TH align=right>$visProtString</TH><TH align=left>&nbsp;$protString ($totProtString max)&nbsp;".&getAmbiguityButton."</TH><TR><TR>$mapString</TR></TABLE>";
		}
		else {
			$statusString=($validStatus==1)? '<B>Partially v' : '<B>V';
			$statusString.="alidated. $historyButtonStrg$lowerScoresString<BR>&nbsp;No proteins selected.</B>";
		}
		##>Match groups status
		if ($verifMG) {
			$stringMG="<B>Verified</B>";
			$stringMG.="&nbsp;&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Set as 'Not verified'\" onclick=\"window.location='./storeProjectItem.cgi?MG=0&ITEM=ANALYSIS&ITEM_ID=$itemID'\">" if ($action eq 'summary' && $projStatus<=0 && $projectAccess ne 'guest');
		}
		else {
			$stringMG="<FONT style=\"font-weight:bold;color:#DD0000\">Not verified</FONT>";
			$stringMG.="&nbsp;&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Set as 'Verified'\" onclick=\"window.location='./storeProjectItem.cgi?MG=1&ITEM=ANALYSIS&ITEM_ID=$itemID'\">" if ($action eq 'summary' && $projStatus<=0 && $projectAccess ne 'guest');
		}
		##>Peptide quantification(s)
		my $sthPQ=$dbh->prepare("SELECT Q.ID_QUANTIFICATION,ID_QUANTIFICATION_METHOD,NAME,STATUS,QUANTIF_ANNOT,IS_REFERENCE FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND ID_ANALYSIS=$itemID AND FOCUS='peptide' ORDER BY Q.ID_QUANTIFICATION");
		my $sthQM=$dbh->prepare("SELECT NAME FROM QUANTIFICATION_METHOD WHERE ID_QUANTIFICATION_METHOD=?");
		my $sthOA=$dbh->prepare("SELECT COUNT(*) FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=?");
		my %quantifMehtods;
		$sthPQ->execute;
		while (my ($quantifID,$quantifMethID,$qName,$status,$quantifAnnot)=$sthPQ->fetchrow_array) {
			my $statusImage=($status==-2)? 'lightRed1.gif' : ($status==-1)? 'lightGray1.gif' : ($status==0)? 'lightYellow1.gif' : 'lightGreen1.gif';
			unless ($quantifMehtods{$quantifMethID}) {
				$sthQM->execute($quantifMethID);
				($quantifMehtods{$quantifMethID})=$sthQM->fetchrow_array;
			}
			$sthOA->execute($quantifID);
			my ($numAna)=$sthOA->fetchrow_array;
			my ($xicSoftCode,$xicSoftware,$version)=&promsQuantif::getXicSoftware($dbh,$quantifID,$quantifAnnot);
			my $versionStrg=($version)? " v$version" :'';
			$pepQuantifStrg.="&nbsp;&bull;<IMG src=\"$promsPath{images}/$statusImage\"/>$qName";
			$pepQuantifStrg.=" x$numAna Analyses" if $numAna > 1;
			$pepQuantifStrg.=" [<B>$quantifMehtods{$quantifMethID}</B> - $xicSoftware$versionStrg]<BR>\n";
		}
		$sthOA->finish;
		$sthQM->finish;
		$sthPQ->finish;
		$pepQuantifStrg='&nbsp;None' unless $pepQuantifStrg;

		##>BioSamples
		my $sthBS=$dbh->prepare("SELECT BS.NAME FROM BIOSAMPLE BS INNER JOIN OBSERVATION O ON BS.ID_BIOSAMPLE=O.ID_BIOSAMPLE WHERE O.ID_ANALYSIS=$itemID ORDER BY BS.NAME");
		$sthBS->execute;
		while (my ($bsName)=$sthBS->fetchrow_array) {
			$bioSampleStrg.="&nbsp;&bull;$bsName<BR>\n";
		}
		$bioSampleStrg='&nbsp;None' unless $bioSampleStrg;
		$sthBS->finish;
	}
	elsif ($validStatus==0) { # not validated
		##>Fetching number of proteins already selected (if on-going validation)
		my ($selProt)=$dbh->selectrow_array("SELECT COUNT(*) FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$itemID AND SEL_STATUS>=1 AND IDENTIFIER NOT LIKE 'DECOY_%'");
		my ($totProt)=$dbh->selectrow_array("SELECT COUNT(*) FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$itemID AND IDENTIFIER NOT LIKE 'DECOY_%'");
		my ($filtProt)=$dbh->selectrow_array("SELECT COUNT(*) FROM PROTEIN_VALIDATION WHERE SEL_STATUS=-3 AND ID_ANALYSIS=$itemID AND IDENTIFIER NOT LIKE 'DECOY_%'");
		my $protString;
		if ($selProt || $valHistory) {
			#my $protString=($selProt>1)? "$selProt proteins" : "$selProt protein";
			#$statusString="<B>Validation in progress ($protString selected so far). $historyButtonStrg";
			$statusString="<B>Validation in progress. $historyButtonStrg";
			my $leftProt=$totProt-$filtProt;
			$protString=($selProt>1)? "<BR>&nbsp;$selProt/$leftProt proteins selected so far" : "&nbsp;$selProt/$leftProt protein selected so far";
			$protString.=($filtProt>1)? "&nbsp;($filtProt excluded proteins)." : ($filtProt==1)? "&nbsp;($filtProt excluded protein)." : '.';
		}
		else {
			$statusString="<B>Not validated. $historyButtonStrg";
			$protString=($filtProt>1)? "<BR>&nbsp;$filtProt excluded proteins." : ($filtProt==1)? "<BR>&nbsp;$filtProt excluded protein." : '';
		}
		$statusString.="$lowerScoresString$protString";
	}
	else { # valid_status=-1 databank data not imported
		$statusString="<B>Protein annotations not yet imported. $historyButtonStrg</B>$lowerScoresString<BR>";
		if ($action eq 'edit') {
			$statusString.="&nbsp;<B>Databank scan:</B>&nbsp;$scanDBString";
		}
	}
	###>FDR
	if ($decoy) {
		my ($decoyMethod,$fdrData)=split(',',$decoy);
		my ($maxFDR,$FDRalgo)=split(/:/,$fdrData) if $fdrData;
		if ($maxFDR) { # if Qvality, DTcount or precomputed
			$maxFDR=~s/FDR=//;
			$FDRalgo='qvality' unless $FDRalgo;
			$fdrString="<B>Targeted FDR:</B> $maxFDR% - ";
			$minScoreString.=($FDRalgo eq 'precomputed')? ' [precomputed FDR]' : ' [FDR-based]';
			#if ($minScore || $FDRalgo eq 'precomputed') {#}
			if (-e "$promsPath{peptide}/proj_$projectID/ana_$itemID/scores.png") {
				$minScoreString.="&nbsp;&nbsp;<INPUT type=\"button\" id=\"scoresBUTTON\" class=\"font11\" value=\"Show scores distribution\" onclick=\"showHideScoresDistribution()\"/>";
				$okScoreDistrib=1;
			}
			if (!$minScore && $FDRalgo ne 'precomputed') {$minScoreString.=" - No decoy matches found.";}
		}
		else {$minScoreString.=' [User-defined]';}
		$fdrString.='<B>Observed FDR:</B> ';
		#if ($validStatus<=0 || !defined $numDecoyPeptides) {$fdrString.='...';}
		if ($validStatus<=0) {
			if ($valHistory) {
				$fdrString.='[Data not reported: ';
				my $numValidPeptides=0;
				my $numValidDecoyPeptides=0;
				###>Checking for manual validation
				my $pepInfoStrg='';
				my $maxRank=($maxRankImport eq '-')? 10 : $maxRankImport; # '-' old data
				foreach my $rank (1..$maxRank) {$pepInfoStrg.=",INFO_PEP$rank";}
				my $sthQ=$dbh->prepare("SELECT QUERY_NUM$pepInfoStrg FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$itemID AND VALID_STATUS >= -2"); #AND VALID_STATUS != -1
				$sthQ->execute || die $sthQ->errstr;
				QUERY:while (my ($queryNum,@pepData)=$sthQ->fetchrow_array) {
					foreach my $pepInfo (@pepData) {
						last unless $pepInfo; # empty field
						if ($pepInfo=~/SEL=(2|-3)/) { # manual selection or rejection
							$fdrString.='N/A (Manual validation).';
							last QUERY;
						}
						if ($pepInfo=~/SEL=1/) { # includes lower-scoring (approx. because best scoring might not be validated)
							#my ($specCount)=($pepInfo=~/SPC=(\d+)/); # top + lower-scoring peptides
							#$specCount=1 unless $specCount; # old data
							#if ($queryNum > 0) {$numValidPeptides+=$specCount;}
							#else {$numValidDecoyPeptides+=$specCount;}
							if ($queryNum > 0) {$numValidPeptides++;} # No longer uses spectral count for FDR calculation #2 (PP 06/08/15)
							else {$numValidDecoyPeptides++;} # No longer uses spectral count for FDR calculation #2 (PP 06/08/15)
						}
					}
				}
				$sthQ->finish;
				if ($numValidPeptides==0) {$fdrString.='No peptides validated.';}
				else {
					if ($numValidDecoyPeptides==0) {$fdrString.='0% (No decoy peptide validated).';}
					else {
						my $fdr=sprintf "%.2f",100*$numValidDecoyPeptides/$numValidPeptides;
						$fdr*=1; # 1.00 -> 1
						$fdrString.="[$numValidPeptides] ~$fdr%"; # ($numValidDecoyPeptides decoy $pepStrg validated).";
					}
				}
				$fdrString.=']';
			}
			else {$fdrString.=' N/A (Validation not started).';}
		}
		elsif (!defined $FDRdata) {$fdrString.='...';}
		elsif ($FDRdata eq '-2') {$fdrString.='N/A (Manual validation).';}
		#elsif ($FDRdata eq '0') {$fdrString.='0% (No decoy peptide validated).';}
		else { # num decoy >=1
			my ($numGoodPeptides,$numDecoyPeptides,$approxStrg,$validStrg);
			if ($FDRdata=~/:/) { # new structure
				my ($d1,$d2,$p1,$p2)=split(':',$FDRdata);
				$numDecoyPeptides=$d1; #+$d2; +N2=spectral count! why? (PP 07/05/15)
				$numGoodPeptides=$p1; #+$p2;
				$approxStrg='';
				$validStrg='selected';
			}
			else { # old structure
				$numDecoyPeptides=$FDRdata;
				($numGoodPeptides)=$dbh->selectrow_array("SELECT COUNT(*) FROM PEPTIDE WHERE ID_ANALYSIS=$itemID AND VALID_STATUS=1"); # valid status cannot =2 since autoval
				$approxStrg='~';
				$validStrg='validated';
			}
			if ($numGoodPeptides) {
				#my $fdr=sprintf "%.2f",100*$numDecoyPeptides/($numGoodPeptides+$numDecoyPeptides);
				my $fdr=sprintf "%.2f",100*$numDecoyPeptides/$numGoodPeptides; # !!!************** Change in decoy calculation (PP 07/08/12) ******************** !!!
				$fdr*=1;
				my $pepStrg=($numDecoyPeptides==1)? 'peptide' : 'peptides';
				$fdrString.="$approxStrg$fdr% ($approxStrg$numDecoyPeptides/$numGoodPeptides $pepStrg $validStrg).";
			}
			else {$fdrString.='N/A (No peptides validated).';}
		}
		$fdrString.='<BR>&nbsp;<B>Decoy method:</B> ';
		if ($decoyMethod eq 'INT:SEARCH') {$fdrString.='Automated decoy search.';}
		elsif ($decoyMethod eq 'INT:BANK') {$fdrString.='Mixed databank with valid and decoy sequences.';}
		else { # External decoy search
			my ($decoyFile)=($decoyMethod=~/EXT:(.+)/);
			$fdrString.="External search file ($decoyFile).";
		}
		if ($maxFDR && $minScore) {
			$fdrString.='<BR>&nbsp;<B>FDR algorithm:</B> ';
			$fdrString.=(!$FDRalgo || $FDRalgo=~/qval/i)? 'Qvality.' : 'DT count.';
		}
	}
	else {$fdrString='N/A (No decoy search).';}
	#}
	#else { # add
	#	@nameText=('Name','Description','Start date','Local data file<SUP>*</SUP>','Databank','Scan databank','Minimal score°','Maximal rank°');
	#	$startDate=$date;
	#}
	#push @nameText,'Comments';
	($name,$description,$comments,$decoy,$labCode)=&promsMod::chkDef($name,$description,$comments,$decoy,$labCode); # also if action=add
	$startDate=&promsMod::formatDate($startDate);
	my $elutionString='';
	if (!$metaAna && $msType eq 'MIS' && ($fileFormat =~ /\.DAT/ || $fileFormat =~ /\.PDM/ || $fileFormat eq 'PHENYX.XML' ) && $projectAccess =~ /bioinfo|mass/ ) { # ($action ne 'add' && ). not manag|super
		if ($instrument && $instrument eq 'MALDI-TOF-TOF') {
			$elutionString="&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Import Elution Data\" onclick=\"window.location='./importElution.cgi?anaList=$itemID&id_project=$itemInfo[0]{ID}'\">";
		}
		my $okExport;
		if ($validStatus<=0) {
			($okExport)=$dbh->selectrow_array("SELECT 1 FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$itemID AND ELUTION_TIME IS NOT NULL LIMIT 1");
		}
		else {
			($okExport)=$dbh->selectrow_array("SELECT 1 FROM PEPTIDE WHERE ID_ANALYSIS=$itemID AND ELUTION_TIME IS NOT NULL LIMIT 1");
		}
		my $disabExp=($okExport)? '' : 'disabled';
		$elutionString.="&nbsp;&nbsp;" unless $elutionString;
		$elutionString.="&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Export Elution Data\" onclick=\"window.location='./exportElution.cgi?ana_ID=$itemID&proj_ID=$itemInfo[0]{ID}'\" $disabExp>";
	}
	if ($action eq 'summary') {
# Temp --> && $userStatus ne 'bio'
		#my $quantiButtonString = ($quantiMethod ne 'None' && $userStatus ne 'bio')? "<INPUT type=\"button\" class=\"font11\" value=\"Quantification Settings\" onclick=\"window.location='./editProtQuantification.cgi?start=1&ACT=listAna&ana=$rowID'\">" : '';
		$name=&promsMod::HTMLcompatible($name);
		push @valueText,"&nbsp;$name&nbsp;";
		push @valueText,"&nbsp;$labCode";
		$description=&promsMod::HTMLcompatible($description);
		push @valueText,"&nbsp;$description";
		push @valueText,"&nbsp;$displayPos";
		push @valueText,"&nbsp;$startDate";
		push @valueText,"&nbsp;$massSpecType{$msType}$elutionString";
		push @valueText,"&nbsp;$searchEngine";

		#my $dataFileStrg=($remoteFile)? "<B>local:</B> $dataFile &nbsp&nbsp&nbsp<B>remote:</B> $remoteFile" : "$dataFile";
		my $dataFileStrg=($fileFormat=~/MAXQUANT/)? "Multiple files" : $dataFile; #no remote file now !
		#if ($msType ne 'SWATH') {
			my $jsShowParam=($showParameters==1)? 0 : 1;
			$dataFileStrg.="&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Show search parameters\" onclick=\"window.location='./editProjectItem.cgi?ACT=summary&ITEM=ANALYSIS&ID=$rowID&showParam=$jsShowParam'\">";
			$dataFileStrg .= " ".getExportAnalysisButton();
		#}
		push @valueText,"&nbsp;$dataFileStrg";
		push @valueText,"&nbsp;$wiffFile";
		push @valueText,"&nbsp;$dataBankString";
		push @valueText,"&nbsp;$taxonomy";
		push @valueText,"&nbsp;$refRetTimeString";
		#push @valueText,"&nbsp;$quantiMethod $quantiButtonString";
		push @valueText,"&nbsp;$labelMethod";
		push @valueText,"&nbsp;$maxRankImport";
		push @valueText,"&nbsp;$minScoreString";
		push @valueText,"&nbsp;$fdrString";
		push @valueText,"&nbsp;$stringMG" if $validStatus>=1;
		push @valueText,"&nbsp;$statusString";
		#if ($statusString !~ /not |in progress/i) { #} fully or partially validated (Not validated | ... annotations not imported | Validation in progress)
		if ($validStatus >= 1) {
			push @valueText,$pepQuantifStrg;
			push @valueText,$bioSampleStrg;
			my $validDateStrg='&nbsp;'.&promsMod::formatDate($validDate);
			$validDateStrg.=" by <B>$validUser</B>" if $projectAccess=~/bioinfo|mass|manag/;
			push @valueText,$validDateStrg;
		}
		$comments=&promsMod::HTMLcompatible($comments);
		push @valueText,"&nbsp;$comments";
	}
	else { # edit (no longer 'add': no single analysis addition)
		push @valueText,"<INPUT type='text' name='name' value='$name' size='50' maxlength='50'>";
		push @valueText,"<INPUT type='text' name='labCode' value='$labCode' size='50' maxlength='50'>";
		push @valueText,"<TEXTAREA name='des' rows='2' cols='65' maxlength='100'>$description</TEXTAREA>";
		if ($action eq 'edit') {
			my ($maxDisplayPos)=$dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM ANALYSIS WHERE ID_SAMPLE=(SELECT ID_SAMPLE FROM ANALYSIS WHERE ID_ANALYSIS=$itemID)");
			my $positionStrg='<SELECT name="displayPos">';
			foreach my $pos (1..$maxDisplayPos) {
				$positionStrg.="<OPTION value=\"$pos\"";
				$positionStrg.=' selected' if $pos==$displayPos;
				$positionStrg.=">$pos</OPTION>";
			}
			$positionStrg.="</SELECT><INPUT type='hidden' name='oldDisplayPos' value='$displayPos'/>";
			push @valueText,$positionStrg;
		}
		push @valueText,$startDate;
		#if ($action eq 'add') {
		#	###>MASCOT file
		#	push @valueText,"<INPUT type='file' name='data_file' size='60'>";
		#	###>Databank
		#	my $sthDB=$dbh->prepare("SELECT D.ID_DATABANK,D.NAME,VERSION_NAME,FASTA_FILE,DT.NAME FROM DATABANK D,DATABANK_TYPE DT WHERE D.ID_DBTYPE=DT.ID_DBTYPE AND USE_STATUS='yes' ORDER BY ID_DATABANK ASC");
		#	$sthDB->execute;
		#	my $refDbData=$sthDB->fetchall_arrayref(); # reference to an array
		#	$sthDB->finish;
		#	my $dbString="<SELECT name=\"databankID\">\n<OPTION selected value=\"-1\">Choose from list</OPTION>\n";
		#	foreach my $refDbRow (@{$refDbData}) {
		#		my ($dbID,$name,$version,$fastaFile,$dbankType)=@{$refDbRow};
		#		if ($fastaFile=~/:/) {
		#			my ($mascotServer,$dbankDir,$fileName)=split(':',$fastaFile);
		#			$version="$mascotServer > $dbankDir";
		#		}
		#		elsif (!-e "$promsPath{banks}/db_$dbID/$fastaFile") {next;}
		#		$dbString.="<OPTION value=\"$dbID\">$name";
		#		$dbString.=" ($version)" if $version;
		#		$dbString.=" [Type: $dbankType]</OPTION>\n";
		#	}
		#	$dbString.="</SELECT>";
		#	push @valueText,$dbString;
		#	push @valueText,$scanDBString;
		#	push @valueText,"<INPUT type='text' name='minScore' value='default' size='10'>&nbsp;&nbsp;(default is set to: <B>".&promsConfig::getMinScore('MASCOT')."</B> for Mascot, <B>".&promsConfig::getMinScore('PHENYX')."</B> for Phenyx)";
		#	my $maxRankOptionString ="";
		#	foreach my $rank (1..10) {
		#		my $IsSelected=($rank==$maxRank)?"selected":"";
		#		$maxRankOptionString .="<OPTION $IsSelected value='$rank'>$rank</OPTION>" ;
		#	}
		#	push @valueText,"<SELECT name='maxRank'>$maxRankOptionString</SELECT>";
		#}
		#else { # edit
			push @valueText,"&nbsp;$massSpecType{$msType}$elutionString";
			push @valueText,"&nbsp;$searchEngine";
			#my $dataFileStrg=($remoteFile)? "<B>local:</B> $dataFile &nbsp&nbsp&nbsp<B>remote:</B> $remoteFile" : "$dataFile";
			my $dataFileStrg="$dataFile";  #no remote file
			#my $changeQuantiMethod;
			#if ($validStatus<=1) { # validation is not complete
			#	my $sthQuanti=$dbh->prepare("SELECT ID_QUANTIF_METHOD,NAME FROM QUANTIF_METHOD");
			#	$sthQuanti->execute;
			#	my %listQuantiMethod;
			#	while (my ($IdQuanti,$QuantiName)=$sthQuanti->fetchrow_array) {
			#		$listQuantiMethod{$IdQuanti}=$QuantiName;
			#	}
			#	$changeQuantiMethod = "<SELECT name=\"changeQuanti\" onchange=\"changeQuantification()\"><OPTION value=\"\">None</OPTION>";
			#	foreach my $localIDQuanti (sort {lc($listQuantiMethod{$a}) cmp lc($listQuantiMethod{$b})} keys %listQuantiMethod) {
			#		my $selected = ($localIDQuanti==$QuantifID)? 'selected' : '';
			#		$changeQuantiMethod .= "<OPTION value=$localIDQuanti $selected>$listQuantiMethod{$localIDQuanti}</OPTION>" ;
			#	}
			#	$changeQuantiMethod .= "</SELECT>\n";
			#}
			#else {$changeQuantiMethod=$quantiMethod;}
			push @valueText,"&nbsp;$dataFileStrg<INPUT type=\"hidden\" name=\"data_file\" value=\"$dataFile\" ><INPUT type=\"hidden\" name=\"file_format\" value=\"$fileFormat\" >";
			push @valueText,"&nbsp;$wiffFile";
			push @valueText,"&nbsp;$dataBankString"; #"<INPUT type=\"hidden\" name=\"databankID\" value=\"$dataBankID\" >";
			push @valueText,"&nbsp;$taxonomy";
			push @valueText,"&nbsp;$refRetTimeString";
			#push @valueText,"&nbsp;$changeQuantiMethod";
			push @valueText,"&nbsp;$labelMethod";
			push @valueText,"&nbsp;$maxRankImport";
			push @valueText,"&nbsp;$minScoreString";
			push @valueText,"&nbsp;$fdrString";
			push @valueText,"&nbsp;$stringMG" if $validStatus>=1;
			push @valueText,"&nbsp;$statusString";
			#if ($statusString !~ /not /i) { #} fully or partially validated (Not validated | ... annotations not imported)
			if ($validStatus >= 1) {
				push @valueText,$pepQuantifStrg;
				push @valueText,$bioSampleStrg;
				my $validDateStrg='&nbsp;'.&promsMod::formatDate($validDate);
				$validDateStrg.=" by <B>$validUser</B>" if $projectAccess=~/bioinfo|mass|manag/;
				push @valueText,$validDateStrg;
			}
		#}
		push @valueText,"<TEXTAREA name='comments' rows='5' cols='65' maxlength='500'>$comments</TEXTAREA>";
	}

	if ($action eq 'edit') { # all samples in project
		($numParents)=$dbh->selectrow_array("SELECT COUNT(ID_SAMPLE) FROM SAMPLE,EXPERIMENT WHERE SAMPLE.ID_EXPERIMENT=EXPERIMENT.ID_EXPERIMENT AND ID_PROJECT=$itemInfo[0]{ID}");
	}
}

################################################
####<ReImport form for a validated Analysis>####
################################################
sub reImportData {
	my $analysisID=param('ID');

	####<Connecting to DB>####
	my $dbh=&promsConfig::dbConnect;
	my ($anaName,$dataFile)=$dbh->selectrow_array("SELECT NAME,DATA_FILE FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
	my $sthDB=$dbh->prepare("SELECT D.ID_DATABANK,DB_RANK,NAME,VERSION_NAME,USE_STATUS FROM ANALYSIS_DATABANK AD,DATABANK D WHERE AD.ID_DATABANK=D.ID_DATABANK AND ID_ANALYSIS=$analysisID");
	$sthDB->execute;
	my (%dbUsed,$dbEnded);
	while (my($dbID,$dbRank,$dbName,$dbVersion,$dbUse)=$sthDB->fetchrow_array) {
		@{$dbUsed{$dbRank}}=($dbID,$dbName,$dbVersion,$dbUse);
		$dbEnded=1 if $dbUse eq 'no';
	}
	#my $dbString=($dbVersion)? "$dbName ($dbVersion)" : $dbName;
	my $refDbData;
	if ($dbEnded) {
		my $sthDB=$dbh->prepare("SELECT ID_DATABANK,NAME,VERSION_NAME FROM DATABANK WHERE USE_STATUS='yes' ORDER BY ID_DATABANK ASC");
		$sthDB->execute;
		$refDbData=$sthDB->fetchall_arrayref; # reference to an array
		$sthDB->finish;
	}

	#	$dbString.=": <FONT color=\"#DD0000\">Not available</FONT>.<BR>Please select a new one :\n";
	#	$dbString.="<SELECT name=\"databankID\">\n<OPTION selected value=\"-1\">Choose from list</OPTION>\n";
	#	foreach my $refDbRow (@{$refDbData}) {
	#		my ($dbID,$name,$version)=@{$refDbRow};
	#		$dbString.="<OPTION value=\"$dbID\">$name";
	#		$dbString.=" ($version)" if $version;
	#		$dbString.="</OPTION>\n";
	#	}
	#	$dbString.="</SELECT>";
	#}
	#else {$dbString.=".<INPUT type=\"hidden\" name=\"databankID\" value=\"$databankID\" >";}

	my $projectID=&promsMod::getProjectID($dbh,$analysisID,'ANALYSIS');
	$dbh->disconnect;

	####<HTML>####
	my ($light,$dark)=&promsConfig::getRowColors;
	(my $JSfile=$dataFile)=~s/\./\\\./;
	print header(-'content-encoding'=>'no',-charset=>'UTF-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT type="text/javascript">
function cancelReImport(fallback) {
	top.promsFrame.selectedAction=fallback;
	top.promsFrame.optionFrame.location='$promsPath{cgi}/selectOption.cgi?ID='+top.promsFrame.selectedBranchID;
}
function checkForm(myForm) {
	for (var i=1; i<=3; i++) {
		if (myForm['databankID'+i] && myForm['databankID'+i].value==-1) {
			alert('You must select a new sequence databank#'+i+'.');
			return false;
		}
	}
	if (!myForm.reval_file_ok) {
		if (!myForm.data_file.value) {
			alert('You must provide a path for search results file $dataFile.');
			return false;
		}
		var okFile=myForm.data_file.value.search(/$JSfile/);
		if (okFile==-1) {
			alert('The file specified is not $dataFile.');
			return false;
		}
	}
	return true;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Reimporting Search Results File for Analysis <FONT color="#DD0000">$anaName</FONT>.</FONT><BR><BR>

<FORM method="post" enctype="multipart/form-data" action="./storeAnalyses.cgi" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="ACT" value="reval" >
<INPUT type="hidden" name="ITEM" value="ANALYSIS" >
<INPUT type="hidden" name="ITEM_ID" value="$analysisID" >
<INPUT type="hidden" name="PROJECT_ID" value="$projectID" >
<!-- <INPUT type="hidden" name="scanDB" value="now"> OBSOLETE -->
<TABLE border=0>
<TR><TD bgcolor=$dark>
 <TABLE border=0 cellpadding=2 width=100%>
|;
	my $numDb=scalar keys %dbUsed;
	foreach my $dbRank (sort{$a<=>$b} keys %dbUsed) {
		print "<TR><TH align=right valign=top bgcolor=$dark width=170>Sequence databank";
		print " #$dbRank" if $numDb > 1;
		my $dbName=$dbUsed{$dbRank}[1];
		$dbName.= " ($dbUsed{$dbRank}[2])" if $dbUsed{$dbRank}[2]; # version
		print " :</TH><TH align=left bgcolor=$light nowrap>$dbName";
		if ($dbUsed{$dbRank}[3] eq 'no') {
			print qq
|: <FONT color="#DD0000">Not available</FONT>.<BR>Please select a new one :
<SELECT name="databankID$dbRank"><OPTION selected value="-1">Choose from list</OPTION>
|;
			foreach my $refDbRow (@{$refDbData}) {
				my ($dbID,$name,$version)=@{$refDbRow};
				print "<OPTION value=\"$dbID\">$name";
				print " ($version)" if $version;
				print "</OPTION>\n";
			}
			print "</SELECT>";
		}
		else {print "<INPUT type=\"hidden\" name=\"databankID$dbRank\" value=\"$dbUsed{$dbRank}[0]\">";}
		print "</TH></TR>\n";
	}
	print qq
|  <TR>
   <TH align=right valign=top bgcolor=$dark nowrap>Search file (<FONT color="#DD0000">$dataFile</FONT>) :</TH>
   <TD bgcolor=$light nowrap>|;
	my $fullDataFile="$promsPath{peptide}/proj_$projectID/ana_$analysisID/$dataFile";
	my $okDataFile=0;
	if (-e $fullDataFile) {
		if (-l $fullDataFile) {
			my $realFile=readlink($fullDataFile);
			$okDataFile=1 if -e $realFile;
		}
		else {$okDataFile=1;}
	}
	if ($okDataFile) {print " <B>Found.</B><INPUT type=\"hidden\" name=\"reval_file_ok\" value=\"1\"><INPUT type=\"hidden\" name=\"data_file\" value=\"$fullDataFile\">";}
	else {print " <B>Not found!<BR>Upload file:</B><INPUT type=\"file\" name=\"data_file\" size=\"70\">";}
	my $fallBack = ($item eq 'METADATA') ? 'edit' : 'summary';
	print qq
|</TD>
  </TR>
  <TR><TH colspan=2>
	<INPUT type="submit" name="save" value=" Proceed " >&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="cancelReImport('$fallBack')" >
  </TH></TR>
 </TABLE>
</TD></TR>
</TABLE>
</FORM>
</CENTER>
</BODY>
</HTML>
|;
	exit;
}


sub downloadSearchFiles {
	print header(-'content-encoding'=>'no',-charset=>'UTF-8');
	warningsToBrowser(1);
	my ($light,$dark)=&promsConfig::getRowColors;
	print qq
|<HTML>
<HEAD>
<TITLE>Downloading search files</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT type="text/javascript">
</SCRIPT>
<BODY>
|;
	my $jobID=strftime("%Y%m%d%H%M%S",localtime);
	if (-e "$promsPath{tmp}/export") {
		##>Clean previous downloads
		opendir(EXPORT,"$promsPath{tmp}/export");
		while (defined (my $child = readdir(EXPORT))) {
			next unless $child=~/^(\d+)_/;
			my $jID=$1;
			if ($jobID-$jID > 86400) { # ~24h (not normal numbers)
				if ($child=~/\.tgz/) {unlink "$promsPath{tmp}/export/$child";} # old archive
				else { # old dir
					unlink glob "$promsPath{tmp}/export/$child/*.txt";
					rmdir "$promsPath{tmp}/export/$child";
				}
			}
		}
		close EXPORT;
	}
	else {mkdir "$promsPath{tmp}/export";}
	my $exportDir=$jobID."_$userID";
	my $exportPath="$promsPath{tmp}/export/$exportDir";
	mkdir $exportPath;
	chdir $exportPath;
#my $tar=Archive::Tar->new();
#$Archive::Tar::FOLLOW_SYMLINK=1; # replace symlinks with real files

	my $projectID=()? $rowID : &promsMod::getProjectID($rowID,$item);
	my %dataFileQueries=(
		PROJECT=>"SELECT A.ID_ANALYSIS,A.FILE_FORMAT,A.DATA_FILE,A.VALID_STATUS,E.NAME,S.NAME FROM EXPERIMENT E,SAMPLE S,ANALYSIS A WHERE E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SAMPLE=A.ID_SAMPLE AND E.ID_PROJECT=?",
		EXPERIMENT=>"SELECT A.ID_ANALYSIS,A.FILE_FORMAT,A.DATA_FILE,A.VALID_STATUS,S.NAME FROM SAMPLE S,ANALYSIS A WHERE S.ID_SAMPLE=A.ID_SAMPLE AND S.ID_EXPERIMENT=?",
		GEL2D=>"SELECT A.ID_ANALYSIS,A.FILE_FORMAT,A.DATA_FILE,A.VALID_STATUS,SP.NAME FROM SPOT SP,SAMPLE S,ANALYSIS A WHERE SP.ID_SPOT=S.ID_SPOT AND S.ID_SAMPLE=A.ID_SAMPLE AND SP.GEL2D=?",
		SPOT=>"SELECT A.ID_ANALYSIS,A.FILE_FORMAT,A.DATA_FILE,A.VALID_STATUS FROM SAMPLE S,ANALYSIS A WHERE S.ID_SAMPLE=A.ID_SAMPLE AND S.ID_SPOT=?",
		SAMPLE=>"SELECT A.ID_ANALYSIS,A.FILE_FORMAT,A.DATA_FILE,A.VALID_STATUS FROM ANALYSIS A WHERE A.ID_SAMPLE=?",
		ANALYSIS=>"SELECT A.ID_ANALYSIS,A.FILE_FORMAT,A.DATA_FILE,A.VALID_STATUS FROM ANALYSIS A WHERE A.ID_ANALYSIS=?"
	);
	my $dbh=&promsConfig::dbConnect;
	my $sthDF=$dbh->prepare($dataFileQueries{$item});
	$sthDF->execute($rowID);
	while (my ($anaID,$fileFormat,$dataFileName,$validStatus,@parentInfo)=$sthDF->fetchrow_array) {
		my $localPath='';
		if ($parentInfo[0]) { # not Ana nor Samp
			foreach my $i (0..$#parentInfo) {
				$localPath.='/' if $localPath;
				$localPath.=quotemeta($parentInfo[$i]);
				mkdir "$exportPath/$localPath" unless -e "$exportPath/$localPath";
			}
		}
		my $searchFile;
		if ($fileFormat=~/\.PDM/) {
			($searchFile=$dataFileName)=~s/_\d+\.*\d*\.pdm/\.msf/;
		}
		elsif ($fileFormat eq 'PHENYX.XML' || $fileFormat eq 'MASCOT.XML') {
			($searchFile=$dataFileName)=~s/\.xml/\.pgf/;
		}
		else {$searchFile=$dataFileName;}
		my $filePath;
		if ($validStatus <= 1) {
			if ($fileFormat=~/\.PDM/) {$filePath="$promsPath{validation}/multi_ana/project_$projectID/$searchFile";}
			else {$filePath="$promsPath{validation}/ana_$anaID/$searchFile";}
		}
		else {$filePath="$promsPath{peptide}/proj_$projectID/ana_$anaID/$searchFile";}
		if ($localPath) {symlink($filePath,"$exportPath/$localPath/$searchFile");}
		else {symlink($filePath,"$exportPath/$searchFile");}
	}
	$sthDF->finish;

	my $archiveFile="dataset_$jobID";

	exit;

}

############################################
####<Compute % missed cleavage peptides>####
############################################
sub computeMissedCleavages {
	
	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);

	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;

	my $anaStrg=($item eq 'ANALYSIS')? $rowID : &getItemAnalyses($dbh,$item,$rowID);
	my $sthMC=$dbh->prepare("SELECT MISS_CUT,VALID_STATUS,COUNT(ID_PEPTIDE) FROM PEPTIDE WHERE ID_ANALYSIS IN ($anaStrg) GROUP BY MISS_CUT,VALID_STATUS");
	$sthMC->execute;
	my %pepCount=(VALID_CUT=>0,VALID_MISS=>0,GHOST_CUT=>0,GHOST_MISS=>0);
	while (my ($missCut,$validStatus,$count)=$sthMC->fetchrow_array) {
		my $isMissKey=($missCut)? 'MISS' : 'CUT';
		my $validKey=($validStatus)? 'VALID' : 'GHOST';
		$pepCount{$validKey.'_'.$isMissKey}+=$count;
	}
	$sthMC->finish;
	$dbh->disconnect;
	
	my $pcValid=1 * sprintf "%.1f",100*$pepCount{VALID_MISS}/($pepCount{VALID_CUT}+$pepCount{VALID_MISS});
	print "<FONT class=\"font11\">Missed cleav.:$pcValid%";
	if ($pepCount{GHOST_CUT} || $pepCount{GHOST_MISS}) {
		my $pcGhost=1 * sprintf "%.1f",100*($pepCount{VALID_MISS}+$pepCount{GHOST_MISS})/($pepCount{VALID_CUT}+$pepCount{VALID_MISS}+$pepCount{GHOST_CUT}+$pepCount{GHOST_MISS});
		print "($pcGhost% with MBR peptides)";
	}
	print '</FONT>';
	exit;
}

###############################################################
####<Compute protein summary for an item (except Analysis)>####
###############################################################
sub proteinSummary {

	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);

	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;

	my $anaStrg=&getItemAnalyses($dbh,$item,$rowID);

	my (%allProteins,%visProteins);
	my $sthProt=$dbh->prepare("SELECT ID_PROTEIN,VISIBILITY FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS IN ($anaStrg)");
	$sthProt->execute;
	while (my ($protID,$vis)=$sthProt->fetchrow_array) {
		$allProteins{$protID}=1;
		$visProteins{$protID}=1 if $vis;
	}
	$sthProt->finish;
	my $totProt=scalar keys %allProteins;
	my $visProt=scalar keys %visProteins;
	my ($numMapped,$numDbTypes)=(0,0);
	if ($totProt) {
		($numMapped)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM PROTEIN WHERE ID_MASTER_PROTEIN IS NOT NULL AND ID_PROTEIN IN (".(join(',',keys %allProteins)).")");			
		($numDbTypes)=$dbh->selectrow_array("SELECT COUNT(DISTINCT DB.ID_DBTYPE) FROM DATABANK DB,ANALYSIS_DATABANK AD WHERE DB.ID_DATABANK=AD.ID_DATABANK AND DB.IS_CRAP=0 AND ID_ANALYSIS IN ($anaStrg)");
	}
	$dbh->disconnect;
	
	my $pcMapped=($totProt && $numMapped)? 1*(sprintf "%.1f",($numMapped/$totProt)*100) : 0;

	print "$totProt,$visProt,$pcMapped,$numDbTypes";
	exit;
}

###################################################
####<Compute FDR for an item (except Analysis)>####
###################################################
sub computeFDR {

	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);

	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;
	my @queryList;
	if ($item eq 'PROJECT') {
		push @queryList,"SELECT A.ID_ANALYSIS,A.FALSE_DISCOVERY FROM EXPERIMENT E,SAMPLE S,ANALYSIS A WHERE E.ID_PROJECT=$rowID AND E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SAMPLE=A.ID_SAMPLE AND S.ID_SPOT IS NULL AND A.VALID_STATUS>=1 AND DECOY IS NOT NULL AND FALSE_DISCOVERY != '-2'";
		push @queryList,"SELECT A.ID_ANALYSIS,A.FALSE_DISCOVERY FROM EXPERIMENT E,GEL2D G,SPOT SP,SAMPLE S,ANALYSIS A WHERE E.ID_PROJECT=$rowID AND E.ID_EXPERIMENT=G.ID_EXPERIMENT AND G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT AND S.ID_SAMPLE=A.ID_SAMPLE AND A.VALID_STATUS>=1 AND DECOY IS NOT NULL AND FALSE_DISCOVERY != '-2'";
	}
	elsif ($item eq 'EXPERIMENT') {
		push @queryList,"SELECT A.ID_ANALYSIS,A.FALSE_DISCOVERY FROM SAMPLE S,ANALYSIS A WHERE S.ID_EXPERIMENT=$rowID AND S.ID_SAMPLE=A.ID_SAMPLE AND S.ID_SPOT IS NULL AND A.VALID_STATUS>=1 AND DECOY IS NOT NULL AND FALSE_DISCOVERY != '-2'";
		push @queryList,"SELECT A.ID_ANALYSIS,A.FALSE_DISCOVERY FROM GEL2D G,SPOT SP,SAMPLE S,ANALYSIS A WHERE G.ID_EXPERIMENT=$rowID AND G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT AND S.ID_SAMPLE=A.ID_SAMPLE AND A.VALID_STATUS>=1 AND DECOY IS NOT NULL AND FALSE_DISCOVERY != '-2'";
	}
	elsif ($item eq 'GEL2D') {
		push @queryList,"SELECT A.ID_ANALYSIS,A.FALSE_DISCOVERY FROM SPOT SP,SAMPLE S,ANALYSIS A WHERE SP.ID_GEL2D=$rowID AND SP.ID_SPOT=S.ID_SPOT AND S.ID_SAMPLE=A.ID_SAMPLE AND A.VALID_STATUS>=1 AND DECOY IS NOT NULL AND FALSE_DISCOVERY != '-2'";
	}
	elsif ($item eq 'SPOT') {
		push @queryList,"SELECT A.ID_ANALYSIS,A.FALSE_DISCOVERY FROM SAMPLE S,ANALYSIS A WHERE S.ID_SPOT=$rowID AND S.ID_SAMPLE=A.ID_SAMPLE AND A.VALID_STATUS>=1 AND DECOY IS NOT NULL AND FALSE_DISCOVERY != '-2'";
	}
	elsif ($item eq 'SAMPLE') {
		push @queryList,"SELECT ID_ANALYSIS,FALSE_DISCOVERY FROM ANALYSIS WHERE ID_SAMPLE=$rowID AND VALID_STATUS>=1 AND DECOY IS NOT NULL AND FALSE_DISCOVERY != '-2'"; # -2: manual validation
	}
	my $sthPep=$dbh->prepare("SELECT COUNT(*) FROM PEPTIDE WHERE ID_ANALYSIS=? AND VALID_STATUS=1");
	my ($numAnalyses,$numGoodPeptides,$numDecoyPeptides,$approxFlag)=(0,0,0,0);
	foreach my $query (@queryList) {
		my $sthAna=$dbh->prepare($query);
		$sthAna->execute;
		while (my ($anaID,$FDRdata)=$sthAna->fetchrow_array) {
			$numAnalyses++;
			my ($decoyPep,$goodPep)=(0,0);
			if ($FDRdata=~/:/) { # new data structure
				my ($d1,$d2,$p1,$p2)=split(':',$FDRdata);
				$decoyPep=$d1; #+$d2;
				$goodPep=$p1; #+$p2;
			}
			else { # old data structure
				$decoyPep=$FDRdata;
				$sthPep->execute($anaID);
				($goodPep)=$sthPep->fetchrow_array;
				$approxFlag=1;
			}
			$numDecoyPeptides+=$decoyPep;
			$numGoodPeptides+=$goodPep;
		}
		$sthAna->finish;
	}
	$sthPep->finish;

	$dbh->disconnect;

	#my $numAllPeptides=$numGoodPeptides+$numDecoyPeptides;
	#if ($numAllPeptides) {
	#	my $fdr=sprintf "%.2f",100*$numDecoyPeptides/$numAllPeptides;
	#	$fdr=~s/\.00//;
	#	my $pepStrg=($numAllPeptides==1)? 'peptide' : 'peptides';
	#	my $anaStrg=($numAnalyses==1)? 'Analysis' : 'Analyses';
	#	print "&nbsp;<B>Observed FDR:</B> $fdr% ($numDecoyPeptides decoy in $numAllPeptides validated $pepStrg based on $numAnalyses $anaStrg).";
	#}
	#else {print "No decoy search.";}

	if ($numGoodPeptides) {
		my $fdr=sprintf "%.2f",100*$numDecoyPeptides/$numGoodPeptides;
		$fdr=~s/\.00//;
		my $pepStrg=($numGoodPeptides==1)? 'peptide' : 'peptides';
		my $anaStrg=($numAnalyses==1)? 'Analysis' : 'Analyses';
		my $validStrg=($approxFlag)? 'approx. validated' : 'selected';
		print "&nbsp;<B>Observed FDR:</B> $fdr% ($numDecoyPeptides/$numGoodPeptides $validStrg $pepStrg based on $numAnalyses $anaStrg).";
	}
	else {print "No decoy search.";}

	exit;
}


####################################
####<Check for on-going mapping>#### .&getRemapButton.
####################################
sub checkOnGoingMapping {
	my $returnStrg='';
	my $jobFlagFile="$promsPath{logs}/mapping.job";
	if (-e $jobFlagFile) {
		open(JOB,$jobFlagFile);
		my @jobInfo=<JOB>;
		close JOB;
		my ($timeStamp,$job,$userID,$projID)=split(/\s/,$jobInfo[0]);
		$returnStrg=($projID==$projectID)? "&nbsp;(mapping is on-going)" : "&nbsp;".&getRemapButton;
	}
	elsif ($item eq 'PROJECT') {$returnStrg="&nbsp;".&getRemapButton;}
}



sub getItemAnalyses {
	my ($dbh,$item,$rowID)=@_;
	my $sthAna;
	if ($item eq 'PROJECT') {
		$sthAna=$dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS A,SAMPLE S,EXPERIMENT E WHERE A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_EXPERIMENT=E.ID_EXPERIMENT AND E.ID_PROJECT=?");
	}
	elsif ($item eq 'EXPERIMENT') {
		$sthAna=$dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS A,SAMPLE S WHERE A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_EXPERIMENT=?");
	}
	elsif ($item eq 'GEL2D') {
		$sthAna=$dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS A,SAMPLE S,SPOT SP WHERE A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_SPOT=SP.ID_SPOT AND SP.ID_GEL2D=?");
	}
	elsif ($item eq 'SPOT') {
		$sthAna=$dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS A,SAMPLE S WHERE A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_SPOT=?");
	}
	elsif ($item eq 'SAMPLE') {
		$sthAna=$dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS WHERE ID_SAMPLE=?");
	}
	my @anaList;
	$sthAna->execute($rowID);
	while (my ($anaID)=$sthAna->fetchrow_array) {push @anaList,$anaID;}
	$sthAna->finish;
	return join(',',@anaList);
}


#################
####<Buttons>####
#################
sub getConflictButton {
	return qq|<INPUT type="button" class="font11" value="Visible & Hidden" onclick="window.location='$promsPath{cgi}/listConflicts.cgi?TYPE=$item&ID=$rowID&ACT=switch'"/>|;
}
sub getAmbiguityButton {
	return qq|<INPUT type="button" class="font11" value="List ambiguities" onclick="window.location='$promsPath{cgi}/listConflicts.cgi?TYPE=$item&ID=$rowID&ACT=ambiguity'"/>|;
}
sub getPtmDistributionButton {
	my $lcItem = lc($item);
	return qq|<INPUT type="button" class="font11" value="Modifications distribution" onclick="window.location='$promsPath{cgi}/selectAnalyses.cgi?ID=$lcItem:$rowID&callType=peptidePTM'"/>|;
}
sub getExportAnalysisButton {
	my $lcItem = lc($item);
	return qq|<INPUT type="button" class="font11" value="Export search file(s)" onclick="window.location='$promsPath{cgi}/selectAnalyses.cgi?ID=$lcItem:$rowID&id_project=$projectID&callType=export'"/>|;
}
sub getRemapButton {
	my $disabButton = ($projectAccess=~/bioinfo|mass/)? " " : " disabled";
	return qq|<INPUT type="button"$disabButton class="font11" value="Re-map identifiers" onclick="if (confirm('Clear and map all identifiers again?')) {window.location='$promsPath{cgi}/editProjectItem.cgi?ACT=force&ID=$rowID&ITEM=PROJECT';}"/>|; ##rowID == projectID
}

####>Revision history<####
# 3.2.10 [BUGFIX] Changed default value of metadata ID to numeric instead of string to avoid NaN error (VS 08/01/20)
# 3.2.9 [FEATURE] Option to compute % of missed-cleavage peptides (PP 15/11/19)
# 3.2.8 [ENHANCEMENT] Optimized SQL queries for Protein summary display (PP 07/11/19)
# 3.2.7 [FEATURE] Add a lock system for experiments (VS 08/08/19)
# 3.2.6 [FEATURE] Allow metadata to have empty categories (VS 28/06/2019)
# 3.2.5 [FIX] Change project data path for metadata path (VS 11/06/19)
# 3.2.4 [FEATURE] Add image type to metadata (VS 06/06/19)
# 3.2.3 [FIX] Bug due to call of &getMetadata in 'add' context (PP/VS 06/06/19)
# 3.2.2 Add project items' metadata handling (VS 01/06/19)
# 3.2.1 Multiple sequence databanks are now listed in order in Analysis summary (PP 31/05/19) 
# 3.2.0 Add export search file button in summary (VS 13/05/19)
# 3.1.9 [FIX] minor display bugs (PP 10/01/19)
# 3.1.8 Restored classical $promsPath{cgi} instead of ${promsPath{cgi} (PP 21/11/18)
# 3.1.7 Fix handling of complex search parameters (hash), used in DIA and TDA (VS 20/11/18)
# 3.1.6 Minor bug fix on paths building (VS 18/10/18)
# 3.1.5 Added check on number of identifier types used in ajaxProteinSummary call (PP 19/09/18)
# 3.1.4 Added count of prot quantifications and exploratory/functional analyses to Experiment summary (PP 01/08/18)
# 3.1.3 [FIX] Minor typo in &getRemapButton (PP 29/06/18)
# 3.1.2 Handles project status=-1 [no auto-end validation] (PP 07/06/18)
# 3.1.1 Minor bug fix in mapped protein value when no validated proteins (PP 13/04/18)
# 3.1.0 Faster display by delaying protein summary through ajax (PP 11/04/2018)
# 3.0.2 Minor modif to display the number of protein present before running massChroQ (MLP 06/04/2018)
# 3.0.1 Minor modif on editing sample's name : to take into account if there are a quote on that name (MLP 06/04/18)
# 3.0.0 Added test on on-going mapping & confirmation on identifier re-mapping (PP 19/03/18)
# 2.9.8 Minor modif to display library export option (MLP 06/12/17)
# 2.9.7 Added peptide quantifications to Analyses<BR>TODO: complete downloadSearchFiles option (PP 11/10/17)
# 2.9.6 Minor modif on $projectID initialization in &analysis (MLP 02/05/17)
# 2.9.5 Minor change in Search file info for MaxQuant (PP 29/12/16)
# 2.9.4 Display search param button for SWATH analysis to show PeakView parameters (MLP 20/12/2016)
# 2.9.3 Compatible with MaxQuant (PP 26/11/16)
# 2.9.2 Uses &cleanParameters function (PP 22/09/16)
# 2.9.1 Displays information for user following Remapping (PP 23/06/16)
# 2.9.0 Added BioSample count in project summary & no search param button for SWATH Analysis (PP 02/05/16)
# 2.8.9 Minor change to print multiple tolerances for MS/MS for MaxQuant (GA 29/04/16)
# 2.8.8 Minor display change (PP 07/04/16)
# 2.8.7 Reference RT detection no longer relies solely on quantification (PP 21/03/16)
# 2.8.6 Added Preferred species selection for Experiment (PP 01/03/16)
# 2.8.5 No longer uses spectral count for FDR calculation also for non-reported data (PP 06/08/15)
# 2.8.4 Reference RT detection and regression plot display [Requires iRT tables] (PP 16/06/15)
# 2.8.3 No longer uses spectral count for FDR calculation (PP 07/05/15)
# 2.8.2 add force update identifier and remap Button (SL 26/01/15)
# 2.8.1 Extends "Peptide Modif Distribution" button to all project items (PP 02/12/14)
# 2.8.0 Added "Peptide Modif Distribution" button for Analysis (PP 25/11/14)
# 2.7.9 Conflict button splitted into Ambiguity and V&H (PP 22/09/14)
# 2.7.8 Lab Code for Analysis & linked biological sample listing (PP 25/07/14)
# 2.7.7 Fix FDR calculation when 0 peptide are validated after reporting (PP 05/06/14)
# 2.7.6 Handles PD percolator precomputed FDR & corrected decoy/target count (PP 07/05/14)
# 2.7.5 Minor text change "errors" to "ambiguities" (PP 23/04/14)
# 2.7.4 Handles new ANALYSIS.FALSE_DISCOVERY data structure (PP 14/04/14)
# 2.7.3 Fix uninitialized variable in FDR data parsing (PP 18/03/14)
# 2.7.2 Minor display changes (PP 05/03/14)
# 2.7.1 Added FDR algorithm for Analysis (PP 04/03/14)
# 2.7.0 Fix minor display bug in conflict button (PP 06/09/13)
# 2.6.9 Button to access conflict-list page (FY 03/09/13)
# 2.6.8 Change Classifications by Custom lists (PP 19/08/13)
# 2.6.7 Minor bug fix due to uninitialized $identMappingID when adding new project (PP 12/07/13)
# 2.6.6 Change getVariableModifications call (promsConfig->promsMod) (GA 07/05/13)
# 2.6.5 Displays lower-scoring peptides activation status (PP 29/04/13)
# 2.6.4 Adds GI Accession to project mapping option (PP 12/03/13)
# 2.6.3 Allows identifier mapping at project creation (PP 28/02/13)
# 2.6.2 Bug correction for multi-databank searches (PP 04/02/13)
# 2.6.1 Multi-databank searches (PP 12/12/12)
# 2.6.0 Added Project status management & new mapping procedure (PP 22/11/12)
# 2.5.9 Minor bug fix for comments display (PP 12/09/12)
# 2.5.8 Added Targeted FDR (qvality) to Analysis display & change in Observed FDR calculation (PP 06/08/12)
# 2.5.7 Fix bug analysis FDR when ghost peptides (PP 04/06/12)
# 2.5.6 Analysis import: Max length of data file name extended to 50 characters (PP 09/01/12)
# 2.5.5 Manager & workgroup management (PP 20/09/11)
# 2.5.4 Labeling Method replaces Quantification for Analyses & bug FDR calculation (PP 07/07/11)
# 2.5.3 Global FDR for any item (PP 14/06/11)
# 2.5.2 Merge FY validation history button code + changes for popup window (PP 12/04/11)
# 2.5.1 Modification of the javascript test in the gel image name so as to consider <BR>with the UpperCase() function more extensions: jpg, jpeg, JPG and JPEG. (GA 05/04/2011)
