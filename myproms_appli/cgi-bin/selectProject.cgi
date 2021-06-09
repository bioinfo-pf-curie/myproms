#!/usr/local/bin/perl -w

################################################################################
# selectProject.cgi         2.5.0                                              #
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
$| = 1;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use strict;
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;

#print header(-'content-encoding'=>'no',-charset=>'utf-8'),"DEBUG<BR>\n"; warningsToBrowser(1); # DEBUG
####################
####>Parameters<####
####################
my $action=param('ACT') || ''; # can be undef
my $display=param('DISPLAY') || '';
my $lastOpen=param('LAST_OPEN') || 30;
my ($color1,$color2)=&promsConfig::getRowColors;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my %massSpecType=&promsConfig::getMsType;
my %dataFiles=('DATA_FILE'=>'Search result file','WIFF_FILE'=>'MS data file');

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

###############################
####>Fetching user details<####
###############################
my $userID=$ENV{'REMOTE_USER'};
my ($userName,$userStatus,$refProfile,$userLab,$userTel,$userEmail,$userInfo,$refMascotIDs,$userWorkgroup)=&promsMod::getUserInfo($dbh,$userID);
#my $superBio=($userStatus eq 'bio' && scalar(@{$refMascotIDs}))? 1 : 0;
if ($userStatus eq 'bio') {
	$action='list';
	$display='project';
	$lastOpen=-1; # no limit for biologists
}

#####################
####>Main Window<####
#####################
if (!$action) { # never for biologist

	####>Scan Analyses for validations to be auto-ended<####
	my $dayDiff=&promsConfig::getMaxPartialValidationDuration; #365; # 1 year # 92 <=> 3 months
	my ($analysisID)=$dbh->selectrow_array("SELECT ID_ANALYSIS FROM ANALYSIS WHERE VALID_STATUS=1 AND DATEDIFF(NOW(),VALID_DATE) >= $dayDiff ORDER BY VALID_DATE LIMIT 1");
	if ($analysisID) {
		my $flagDir="$promsPath{tmp}/scratch/autoEndValidation";
		if (-e $flagDir) {
			&promsMod::cleanDirectory($flagDir,'24h'); # remove every file older than 24h
		}
		else {
			mkdir "$promsPath{tmp}/scratch" unless -e "$promsPath{tmp}/scratch";
			mkdir $flagDir;
		}
		unless (glob "$flagDir/*.flag") { # no existing flag file => proceed
			my ($projectID)=&promsMod::getProjectID($dbh,$analysisID,'analysis');
			my ($projStatus)=$dbh->selectrow_array("SELECT STATUS FROM PROJECT WHERE ID_PROJECT=$projectID");
			$projStatus=0 unless $projStatus;
			my $sthAna=$dbh->prepare("SELECT A.ID_ANALYSIS FROM EXPERIMENT E,SAMPLE S,ANALYSIS A
									 WHERE E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SAMPLE=A.ID_SAMPLE AND E.ID_PROJECT=$projectID
									 AND VALID_STATUS=1 AND DATEDIFF(NOW(),VALID_DATE) >= $dayDiff
									 ORDER BY A.ID_ANALYSIS");
			$sthAna->execute;
			my @anaList;
			while (my ($anaID) = $sthAna ->fetchrow_array) {
				push @anaList,$anaID;
			}
			$sthAna->finish;

			if ($projStatus==0) { # usual case
				##>Launch process in background<##
				my $jobID=strftime("%Y%m%d%H%M%S",localtime).'_'.$projectID;
				my $anaListStrg=join(':',@anaList);
				system "echo $anaListStrg > $flagDir/$jobID.flag";
				system "./send2Biologist.cgi $jobID $anaListStrg &";
			}
			elsif ($projStatus==-1) { # set for no auto-update => reset valid date
				my $sthUpAna=$dbh->prepare("UPDATE ANALYSIS SET VALID_DATE=NOW() WHERE ID_ANALYSIS=?");
				foreach my $anaID (@anaList) {
					$sthUpAna->execute($anaID);
				}
				$sthUpAna->finish;
				$dbh->commit;
			}
		}
	}

	$dbh->disconnect;

	####>HTML<####
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Select Project - Menu</TITLE>
<SCRIPT LANGUAGE="JavaScript">
//Clear all tree status
top.projectTreeStatus='';
top.itemTreeStatus='';
</SCRIPT>
</HEAD>
<FRAMESET cols="40%,*" border=0>
	<FRAME name="ownerFrame" src="./selectProject.cgi?ACT=wg">
	<FRAME name="projectFrame" src="$promsPath{html}/nothing.html">
</FRAMESET>
</HTML>
|;
	exit;
}

####>SQL query for last time open filter (changed for display=analysis)<####
my $sqlTimeStrg=($lastOpen <= 0 || $action eq 'all' || $action eq 'archive')? '' : "AND DATEDIFF(NOW(),LAST_OPEN) <= $lastOpen";


#######################
####>Action Window<####
#######################
if ($action ne 'list') {
	my @actionList=(['wg','Browse workgroups from projects'],['owner','Browse owners of projects'],['all','List all projects'],['archived','List all archived projects'],['analysis','List on-going analyses'],['search','Search']);
	my @timeList=([1,'day'],[3,'3 days'],[7,'week'],[30,'month'],[90,'3 months'],[180,'6 months'],[365,'year'],[1095,'3 years'],[-1,'**any**']);

	####>Fetching list of owners<####
	my %ownerList;
	if ($action eq 'wg') {
		my $sthOw=$dbh->prepare("SELECT 1,WORK_GROUP,OWNER,STATUS FROM PROJECT WHERE ID_PROJECT=? $sqlTimeStrg");
		foreach my $projID (keys %{$refProfile}) {
			$sthOw->execute($projID);
			my ($match,$workgroup,$pjOwner,$status)=$sthOw->fetchrow_array;
			next if (!$match || ($status && $status > 1)); # archived project
			$workgroup='' unless $workgroup;
			$pjOwner='*No owner*' unless $pjOwner;
			$ownerList{$workgroup}{$pjOwner}=1;
		}
		$sthOw->finish;
	}
	elsif ($action eq 'owner') {
		my $sthOw=$dbh->prepare("SELECT 1,OWNER,STATUS FROM PROJECT WHERE ID_PROJECT=? $sqlTimeStrg");
		foreach my $projID (keys %{$refProfile}) {
			$sthOw->execute($projID);
			my ($match,$pjOwner,$status)=$sthOw->fetchrow_array;
			next if (!$match || ($status && $status > 1)); # archived project
			$pjOwner='*No owner*' unless $pjOwner;
			$ownerList{'all'}{$pjOwner}=1;
		}
		$sthOw->finish;
	}
	$dbh->disconnect;

	####>HTML<####
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Projects List</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
A {font-weight:bold;}
A:hover {text-decoration:none;}
A.wg,A.selWg {font-size:20px;}
A.owner,A.selOwner {font-size:18px;}
A.wg,A.owner {color:black;}
A.wg:hover,A.owner:hover {color:#DD8800;}
A.selWg,A.selOwner {color:#DD0000;}
IMG.buttonBig {border-width:3px;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
function showHideOwners(tableID) {
	var ownerTable=document.getElementById(tableID);
	ownerTable.style.display=(ownerTable.style.display=='none')? 'block' : 'none';
}
function selectOwner(ownerID,owner,workgroup) {
	updateSelection(ownerID);
	if (!workgroup) workgroup='None';
	//parent.projectFrame.location="$promsPath{cgi}/selectProject.cgi?ACT=list&LAST_OPEN=$lastOpen&DISPLAY=project&OWNER="+owner+"&WG="+workgroup;
	var myForm=document.selectForm;
	myForm.OWNER.value=owner;
	myForm.WG.value=workgroup;
	myForm.submit();
}
function selectWorkgroup(wgID,workgroup) {
	updateSelection(wgID);
	if (!workgroup) workgroup='None';
	//parent.projectFrame.location="$promsPath{cgi}/selectProject.cgi?ACT=list&LAST_OPEN=$lastOpen&DISPLAY=project&WG="+workgroup;
	var myForm=document.selectForm;
	myForm.OWNER.value='';
	myForm.WG.value=workgroup;
	myForm.submit();
}
function updateSelection(newID) {
	if (selectedID == newID) return;
	var oldSelected = document.getElementById(selectedID);
	var newSelected = document.getElementById(newID);
	if (oldSelected) {
		oldSelected.className=(selectedID.match('wg:'))? 'wg' : 'owner';
	}
	newSelected.className=(newID.match('wg:'))? 'selWg' : 'selOwner';
	selectedID = newID;
}
function selectDisplay(action,lastOpen) {
	if (!action) {action='$action';} // change in LAST_OPEN
	if (!lastOpen) {lastOpen='$lastOpen';} // change in ACT
	if (action=='search') { //action=='owner' \|\| action=='wg' \|\|
		parent.projectFrame.location="$promsPath{html}/nothing.html";
	}
	else {
		if (action=='all') {
			parent.projectFrame.location="$promsPath{cgi}/selectProject.cgi?ACT=list&DISPLAY=project&LAST_OPEN=-1";
		}
		else if (action=='archived') {
			parent.projectFrame.location="$promsPath{cgi}/selectProject.cgi?ACT=list&DISPLAY=archived&LAST_OPEN=-1";
		}
		else if (action=='analysis') {
			parent.projectFrame.location="$promsPath{cgi}/selectProject.cgi?ACT=list&DISPLAY=analysis&LAST_OPEN="+lastOpen;
		}
	}
	window.location="./selectProject.cgi?ACT="+action+"&LAST_OPEN="+lastOpen;
}
function createProject() {
	top.promsFrame.location='$promsPath{cgi}/editProjectItem.cgi?ACT=add&ITEM=project';
}
function checkForm(myForm) {
	if (myForm.searchString.value=='') {
		alert('Type a search string.');
		return false;
	}
	if (myForm.name.checked \|\| myForm.file.checked \|\| myForm.protId.checked \|\| myForm.protDes.checked) {return true;}
	alert('Select a search field.');
	return false;
}
var selectedID='';
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR>
<DIV padding=4px>
<INPUT type="button" style="font-weight=bold;font-size:20px;width:185px;" value="Main Window" onclick="parent.location='./promsMain.cgi'">
<INPUT type="button" style="font-weight=bold;font-size:20px;width:185px;" value="Close Session" onclick="top.closeSession()">
<FORM name="selectForm" method="post" target="projectFrame">
<INPUT type="hidden" name="ACT" value="list">
<INPUT type="hidden" name="DISPLAY" value="project">
<INPUT type="hidden" name="LAST_OPEN" value="$lastOpen">
<INPUT type="hidden" name="WG" value="">
<INPUT type="hidden" name="OWNER" value="">
</FORM>
<BR><INPUT type="button" value="Create a new Project" style="font-weight=bold;font-size:20px;" onclick="createProject()"><BR><FONT class="title2">or</FONT><BR>
<TABLE border=0 cellspacing=0 cellpadding=0>
<TR bgcolor="$color2"><TD colspan=2>
	<TABLE cellpadding=4>
|;
	print qq
|<TR><TD nowrap><SELECT class="title2 action" onchange="selectDisplay(this.value)">
|;
	foreach my $refAction (@actionList) {
		my $selStrg=($refAction->[0] eq $action)? ' selected' : '';
		print "\t<OPTION value=\"$refAction->[0]\"$selStrg>$refAction->[1]</OPTION>\n";
	}
	print "</SELECT>\n";

	if ($action !~ /^(all|archived|search)$/) {
		my $timeTypeStrg=($action eq 'analysis')? 'imported' : 'open';
		print qq
|&nbsp;<FONT class="title2">$timeTypeStrg within last</FONT>&nbsp;<SELECT class="title2 action" onchange="selectDisplay(null,this.value)">
|;
		foreach my $refLimit (@timeList) {
			my $selStrg=($refLimit->[0]==$lastOpen)? ' selected' : '';
			print "\t<OPTION value=\"$refLimit->[0]\"$selStrg>$refLimit->[1]</OPTION>\n";
		}
		print "</SELECT>\n";
	}
	print qq
|</TD></TR>
	</TABLE>
</TD></TR>
<TR><TD colspan=2><BR></TD></TR>
|;
	###>List of owners<###
	my $matchedWorkgroup=''; # wgName
	my $matchedWgID=''; # pseudoID: wg:wgRank
	my $matchedOwnerID=''; # pseudoID: owner:ownerRank
	my $firstWgID='';
	my $firstWgName='';
	my $firstOwnerID='';
	my $firstUserName='';
	if ($action=~/wg|owner/) {
		my $countWg=0;
		my $countOw=0;
		my $defaultDisplay;
		foreach my $workgroup (sort{lc($a) cmp lc($b)} keys %ownerList) {
			$countWg++;
			if ($action eq 'wg') {
				$defaultDisplay='none';
				if ($workgroup eq $userWorkgroup || $countWg==1) { # select user's wg (manager only) or 1st in list
					$matchedWorkgroup=$workgroup;
					$matchedWgID="wg:$countWg";
				}
				my $wkGrStrg=($workgroup)? "Wg: $workgroup" : 'No Workgroup';
				#print "<TR><TD colspan=2 valign=middle nowrap><IMG src='$promsPath{images}/blank.gif' width=20 height=40><BUTTON style=\"width:49px\" onclick=\"showHideOwners('table_wg:$countWg')\"><IMG src='$promsPath{images}/owners.gif' width=39 height=40></BUTTON>&nbsp;<A id=\"wg:$countWg\" class=\"wg\" href=\"javascript:selectWorkgroup('wg:$countWg','$workgroup');\">$wkGrStrg</A></TD></TR>\n";
				print "<TR><TD colspan=2 valign=middle nowrap><IMG src='$promsPath{images}/blank.gif' width=20 height=40><IMG class=\"button buttonBig\" src='$promsPath{images}/owners.gif' width=39 height=40 onclick=\"showHideOwners('table_wg:$countWg')\">&nbsp;<A id=\"wg:$countWg\" class=\"wg\" href=\"javascript:selectWorkgroup('wg:$countWg','$workgroup');\">$wkGrStrg</A></TD></TR>\n";
			}
			else { # owner
				$defaultDisplay='block';
				print "<TR><TD colspan=2 valign=middle nowrap><IMG src='$promsPath{images}/blank.gif' width=20 height=40><IMG src='$promsPath{images}/owners.gif'>&nbsp;<FONT class=\"title1\">All owners</FONT></TD></TR>\n";
			}
			my $numOwners=scalar keys %{$ownerList{$workgroup}};
			my $count=0;
			print "<TR><TD colspan=2 nowrap><TABLE id='table_wg:$countWg' border=0 cellspacing=0 cellpadding=0 style='display:$defaultDisplay'>\n";
			foreach my $pjOwner (sort{lc($a) cmp lc($b)} keys %{$ownerList{$workgroup}}) {
				$count++;
				$countOw++;
				my $nodeImg=($count==$numOwners)? 'lastnode.gif' : 'node.gif';
				#my $usedWg=($action eq 'wg')? $workgroup : '';
				print "<TR><TD nowrap><IMG src='$promsPath{images}/blank.gif' width=26 height=22><IMG src='$promsPath{images}/$nodeImg' width=16 height=22><A id=\"owner:$countOw\" class=\"owner\" style=\"vertical-align:top\" href=\"javascript:selectOwner('owner:$countOw','$pjOwner','$workgroup');\">&nbsp;$pjOwner&nbsp;</A></TD></TR>\n";
				if ($countOw==1) { # && $action eq 'owner'
					$firstWgID="wg:$countWg";
					$firstWgName=$workgroup;
					$firstOwnerID="owner:$countOw";
					$firstUserName=$pjOwner;
				}
				if (!$matchedOwnerID && $pjOwner eq $userName) { # owner name can appear more than once
					$matchedOwnerID="owner:$countOw";
					$matchedWorkgroup=$workgroup;
					$matchedWgID="wg:$countWg";
				}
			}
			print "<TR><TD>&nbsp;</TD></TR>\n</TABLE></TD></TR>\n";
		}
		unless (scalar keys %ownerList) {
			print "<TR><TD width=45 valign=middle nowrap><IMG src='$promsPath{images}/blank.gif' width=26 height=22><IMG src='$promsPath{images}/lastnode.gif' width=16 height=22></TD><TD valign=middle nowrap width=100%><FONT class=\"title2\">No projects found</FONT></TD></TR>\n";
		}
	}

	###>Search form<###
	elsif ($action eq 'search') {
		print qq
|<FORM name="searchForm" method="post" onsubmit="return checkForm(this);" target="projectFrame">
<INPUT type="hidden" name="ACT" value="list">
<INPUT type="hidden" name="DISPLAY" value="search">
<TR bgcolor="$color2"><TD>
<TABLE border=0 cellpadding=0 cellspacing=5>
	<TR>
	<TH align=right valign=top width="120px"><FONT style="font-size:18px;">Search for :</FONT></TH>
	<TD align=left><INPUT type='text' name="searchString" value="" size="45"><BR>
	<FONT style="font-size:11px;">Match is <B>not</B> case sensitive.<BR>
	<B>Match :</B><SELECT name="allAnyEx" style="font-size:11px;font-weight:bold;"><OPTION value="AND">all words</OPTION><OPTION value="OR">any word</OPTION><OPTION value="exact">exactly</OPTION></SELECT>
	</FONT></TD>
	</TR>
	<TR>
	<TH align=right valign=top><FONT style="font-size:18px;">in :</FONT></TH>
	<TH align=left>
	<INPUT type="checkbox" name="name" value="1" id="name">&nbsp;<LABEL for="name">Item name (Project, Experiment, ...)</LABEL><BR>
	<INPUT type="checkbox" name="file" value="1" id="file">&nbsp;<LABEL for="file">Data files (raw data or search results)</LABEL><BR>
	<INPUT type="checkbox" name="protId" value="1" id="protId">&nbsp;<LABEL for="protId">Protein identifier</LABEL><BR>
	<INPUT type="checkbox" name="protDes" value="1" id="protDes">&nbsp;<LABEL for="protDes">Protein description</LABEL>
	</TH>
	</TR>
	<TR>
	<TH align=right valign=top><FONT style="font-size:18px;">Restrict to :</FONT></TH>
	<TH align=left><SELECT name="restrict" style="font-weight:bold;"><OPTION value="0">all</OPTION><OPTION value="1">owned</OPTION></SELECT> projects.
	<INPUT type="checkbox" name="archived" value="1" id="archived"/><LABEL for="archived">Include archived projects</LABEL></TH>
	</TR>
	<TR><TD colspan=2 align=center><INPUT type="submit" name="search" value="Search" style="width:100px"/></TD></TR>
</TABLE></TD><TD></TD></TR>
</FORM>
|;
	}
	print "</TABLE>\n";
	if ($matchedOwnerID) {
		print qq
|<SCRIPT LANGUAGE="JavaScript">
document.getElementById('table_$matchedWgID').style.display='block';
selectOwner('$matchedOwnerID','$userName','$matchedWorkgroup');
</SCRIPT>
|;
	}
	elsif ($matchedWorkgroup) { # $action = 'wg'
		print qq
|<SCRIPT LANGUAGE="JavaScript">
selectWorkgroup('$matchedWgID','$matchedWorkgroup');
</SCRIPT>
|;
	}
	elsif ($firstOwnerID) {
		if ($action eq 'wg') {
			print qq
|<SCRIPT LANGUAGE="JavaScript">
document.getElementById('table_$firstWgID').style.display='block';
selectOwner('$firstOwnerID','$firstUserName','$firstWgName');
</SCRIPT>
|;
		}
		else { # owner
			print qq
|<SCRIPT LANGUAGE="JavaScript">
selectOwner('$firstOwnerID','$firstUserName','all');
</SCRIPT>
|;
		}
	}
	elsif ($action=~/wg|owner/) {
		print qq
|<SCRIPT LANGUAGE="JavaScript">
//parent.projectFrame.location="$promsPath{html}/nothing.html";
parent.projectFrame.location="$promsPath{cgi}/selectProject.cgi?ACT=list&DISPLAY=project&LAST_OPEN=$lastOpen";
</SCRIPT>
|;
	} # other actions are triggered in parallele by calling script
	print qq
|</CENTER>
</BODY>
</HTML>
|;
	exit;

}


#######################
####>Search Window<####
#######################
if ($display eq 'search') {

#print header(-'content-encoding'=>'no',-charset=>'utf-8'),"DEBUG<BR>\n"; warningsToBrowser(1); # DEBUG

	my $searchString=param('searchString');
	my $allAnyEx=param('allAnyEx');
	my $search_name=param('name');
	my $search_file=param('file');
	my $search_protId=param('protId');
	my $search_protDes=param('protDes');
	my $restrict=param('restrict');
	my $includeArchived=param('archived');

	####>Processing search parameters<####
	$searchString=uc($searchString);

	my @searchWords=($allAnyEx eq 'exact')? ($searchString) : split(/\s+/,$searchString);
	for (my $i=0;$i<=$#searchWords;$i++) {$searchWords[$i]="LIKE '%$searchWords[$i]%'";}

	my %matchingItems;

	####>Fetching info on item names<####
	if ($search_name) {
		foreach my $item ('PROJECT','EXPERIMENT','GEL2D','SPOT','SAMPLE','ANALYSIS') {
			my $restrictString=($restrict && $item eq 'PROJECT')? "AND OWNER='$userName'" : '';
			my $qString=&createQuery('NAME',$allAnyEx,@searchWords);
			my $sthN=$dbh->prepare("SELECT ID_$item FROM $item WHERE $qString $restrictString");
			$sthN->execute;
			while (my ($itemID)=$sthN->fetchrow_array) {
				%{$matchingItems{$item}{$itemID}}=();
			}
			$sthN->finish;
		}
	}
	####>Fetching info on files<####
	if ($search_file) {
		foreach my $file (keys %dataFiles) {
			my $qString=&createQuery($file,$allAnyEx,@searchWords);
			my $sthF=$dbh->prepare("SELECT ID_ANALYSIS,$file FROM ANALYSIS WHERE $qString");
			$sthF->execute;
			while (my ($anaID,$fileValue)=$sthF->fetchrow_array) {
				push @{$matchingItems{'ANALYSIS'}{$anaID}{'FILE'}},[$file,$fileValue];
			}
			$sthF->finish;
		}
	}

	####>Fetching info on proteins<####
	if ($search_protId) {
		my $qString1=&createQuery('IDENTIFIER',$allAnyEx,@searchWords);
		my $qString2=&createQuery('ALIAS',$allAnyEx,@searchWords);
		my $sthId=$dbh->prepare("SELECT ID_ANALYSIS,PROTEIN.ID_PROTEIN,IDENTIFIER,ALIAS,PROT_DES,PROT_LENGTH,ORGANISM FROM PROTEIN,ANALYSIS_PROTEIN WHERE PROTEIN.ID_PROTEIN=ANALYSIS_PROTEIN.ID_PROTEIN  AND ($qString1 OR $qString2)");
		$sthId->execute;
		while (my ($anaID,$protID,@protInfo)=$sthId->fetchrow_array) {
			@{$matchingItems{'ANALYSIS'}{$anaID}{'PROTEIN'}{$protID}}=@protInfo;
		}
		$sthId->finish;
	}
	if ($search_protDes) {
		my $qString=&createQuery('PROT_DES',$allAnyEx,@searchWords);
		my $sthDes=$dbh->prepare("SELECT ID_ANALYSIS,PROTEIN.ID_PROTEIN,IDENTIFIER,ALIAS,PROT_DES,PROT_LENGTH,ORGANISM FROM PROTEIN,ANALYSIS_PROTEIN WHERE PROTEIN.ID_PROTEIN=ANALYSIS_PROTEIN.ID_PROTEIN AND $qString");
		$sthDes->execute;
		while (my ($anaID,$protID,@protInfo)=$sthDes->fetchrow_array) {
			@{$matchingItems{'ANALYSIS'}{$anaID}{'PROTEIN'}{$protID}}=@protInfo;
		}
		$sthDes->finish;
	}

	####>Processing search results<####
	my $sthPj=$dbh->prepare("SELECT DES,OWNER,WORK_GROUP,LAST_OPEN,LAST_USER FROM PROJECT WHERE ID_PROJECT=?");
	my $sthUser=$dbh->prepare("SELECT USER_NAME FROM USER_LIST WHERE ID_USER=?");
	my (%projectList,%itemInfo,%okProjects,%lastUserNames);
	foreach my $item (keys %matchingItems) {
# print ">$item:<BR>\n";
		foreach my $itemID (keys %{$matchingItems{$item}}) {
			@{$itemInfo{"$item:$itemID"}}=&promsMod::getItemInfo($dbh,$item,$itemID);
			my $projStatus=$itemInfo{"$item:$itemID"}[0]{'STATUS'};
			next if (!$includeArchived && $projStatus==2); # project is archived
			my $projID=$itemInfo{"$item:$itemID"}[0]{'ID'};
			next if (defined($okProjects{$projID}) && $okProjects{$projID}==0); # manag or super bio
			unless ($projectList{$projID}) { # done only once
				$sthPj->execute($projID);
				my ($des,$pjOwner,$wkgr,$lastOpenDate,$lastUserID)=$sthPj->fetchrow_array;
				$des='None' unless $des;
				$wkgr='' unless $wkgr;
				my $lastOpenStrg='Never'; # default
				if ($lastOpenDate) {
					unless ($lastUserNames{$lastUserID}) {
						$sthUser->execute($lastUserID);
						($lastUserNames{$lastUserID})=$sthUser->fetchrow_array;
					}
					$lastOpenStrg=&promsMod::formatDate($lastOpenDate).' by '.$lastUserNames{$lastUserID};
				}
				@{$projectList{$projID}{'INFO'}}=($itemInfo{"$item:$itemID"}[0]{'NAME'},$des,$pjOwner,$wkgr,$projStatus,$lastOpenStrg); # name
			}
			unless (defined($okProjects{$projID})) {
				$okProjects{$projID}=0;
				if ($userStatus=~/(manag|bio\Z)/) {
					if ($userStatus eq 'manag') {
						$okProjects{$projID}=1 if $userWorkgroup eq $projectList{$projID}{'INFO'}[3];
					}
					unless ($okProjects{$projID}) {
						$okProjects{$projID}=1 if ($refProfile->{$projID} && $refProfile->{$projID}=~/super_/);
					}
				}
				else {$okProjects{$projID}=1;} # bioinfo & mass
			}
			push @{$projectList{$projID}{'ITEM'}{$item}},[$itemInfo{"$item:$itemID"}[-1]{'ID'},$itemInfo{"$item:$itemID"}[-1]{'NAME'}] if ($okProjects{$projID} && $item ne 'PROJECT');
		}
	}
	$sthPj->finish;
	$sthUser->finish;
	$dbh->disconnect;

	####>HTML<####
	my $bgColor=$color1;
	my ($popBgColor,$popBdColor)=&promsConfig::getPopupColors;
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Analyses List</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TD {font-size:13px;}
A.protein:hover{color:#DD0000;text-decoration:none;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
|;
&promsMod::popupInfo(400,'rel',-200,'rel',20);
print qq
|function openProject(projID,branchID,itemBranchID) {
	var URLstrg="$promsPath{cgi}/openProject.cgi?ACT=open&ID="+projID+"&branchID="+branchID;
	if (itemBranchID) {URLstrg+="&itemBranchID="+itemBranchID;}
	top.promsFrame.location=URLstrg;
}
function restoreProject(projID) {
	if (confirm('Restore Project?')) {
		window.location="./editBranch.cgi?ITEM=PROJECT&ID="+projID+"&ACT=continue";
	}
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR><BR>
<FONT class="title1">List of matching projects</FONT>
<TABLE width=95% cellspacing=0 border=0>
|;
	###>Looping through all projects<###
	my $numMatches=0;
	foreach my $projID (sort{lc($projectList{$a}{'INFO'}[0]) cmp lc($projectList{$b}{'INFO'}[0]) || $a<=>$b} keys %projectList) { # alphabetical sort on project_name
		next unless $okProjects{$projID};
		$numMatches++;
		my ($name,$des,$pjOwner,$wkgr,$status,$lastOpen)=@{$projectList{$projID}{'INFO'}};
		$wkgr='None' unless $wkgr;
		my $image=($status==1)? 'lightGreen1.gif' : ($status==2)? 'lightRed1.gif' : ($status==-1)? 'lightOrange1.gif' : 'lightYellow1.gif';
		print qq
|<TR bgcolor="$bgColor">
<TD width=10>&nbsp;</TD>
<TD><IMG src="$promsPath{images}/$image"/>&nbsp;<FONT class="title2">$name</FONT>
<BR><B>&nbsp;&nbsp;Description:</B> $des
<BR><B>&nbsp;&nbsp;Workgroup:</B> $wkgr
<BR><B>&nbsp;&nbsp;Owner:</B> $pjOwner
<BR><B>&nbsp;&nbsp;Last time open:</B> $lastOpen</TD>
<TH width=110>
|;
		if ($status==2) {
			my $disabStrg=($refProfile->{$projID}=~/bioinfo|mass|manag/)? '' : 'disabled'; # manager can be simple bio on some projects
			print "<INPUT type=\"button\" class=\"title3\" style=\"width:100px\" value=\"Restore\" onclick=\"restoreProject($projID)\" $disabStrg >";
		}
		else {print "<INPUT type=\"button\" class=\"title3\" style=\"width:60px\" value=\"Open\" onclick=\"openProject($projID,'project:$projID')\">";}
		print "</TH>\n</TR>\n";

		if ($projectList{$projID}{'ITEM'}) {
			foreach my $item ('EXPERIMENT','GEL2D','SPOT','SAMPLE','ANALYSIS') { # loop through matched project items if any
				next unless $projectList{$projID}{'ITEM'}{$item};
				foreach my $refItemInfo (@{$projectList{$projID}{'ITEM'}{$item}}) {
					my ($itemID,$itemName)=@{$refItemInfo};
					#my $iconeImg;
					#if ($item eq 'ANALYSIS') {
					#	my $validStatus=$itemInfo{"$item:$itemID"}[-1]{'VALID'};
					#	$iconeImg=$iconeImages{"analysis:$analysisClass{$validStatus}"}
					#}
					#else {$iconeImg=$iconeImages{lc($item)};}
					my $itemType=$itemInfo{"$item:$itemID"}[-1]{'TYPE'};
					print "<TR bgcolor=\"$bgColor\"><TD>&nbsp</TD>\n";
					#print "<TD><HR width=80%>&nbsp;&nbsp;<IMG src=\"$promsPath{images}/$iconeImg\"><FONT class=\"title3\">$itemType:</FONT> $itemName";
					print "<TD><HR width=80%>&nbsp;&nbsp;<FONT class=\"title3\">$itemType:</FONT> $itemName";
					my $alignStrg;
					if ($item eq 'ANALYSIS' && scalar keys %{$matchingItems{'ANALYSIS'}{$itemID}}) { # loop through match analysis items if any
						print "<BR><TABLE>";
						if ($matchingItems{'ANALYSIS'}{$itemID}{'FILE'}) {
							print "<TR><TD colspan=2>";
							foreach my $refFile (@{$matchingItems{'ANALYSIS'}{$itemID}{'FILE'}}) {
								my ($fileType,$fileName)=($dataFiles{$refFile->[0]},$refFile->[1]);
								print "<B>&nbsp;&nbsp;&nbsp;&nbsp;$fileType:</B> $fileName<BR>\n";
							}
							print "</TD></TR>\n";
						}
						if ($matchingItems{'ANALYSIS'}{$itemID}{'PROTEIN'}) {
							print "<TR><TD nowrap valign=top>&nbsp;&nbsp;&nbsp;&nbsp;<B>Matching proteins:</B></TD><TD>";
							my $protString='';
							foreach my $protID (keys %{$matchingItems{'ANALYSIS'}{$itemID}{'PROTEIN'}}) {
								my ($identifier,$alias,$protDes,$protSize,$species)=@{$matchingItems{'ANALYSIS'}{$itemID}{'PROTEIN'}{$protID}};
								$protString.=', ' if $protString;
								$protString.="<A class=\"protein\" href=\"javascript:void(null);\" onmouseover=\"popup('<B>Protein:</B> $alias<BR><B>Identifier:</B> $identifier<BR><B>Description:</B> $protDes.<BR><B>Size:</B> $protSize aa.<BR><B>Organism:</B> <I><U>$species</U></I>.')\" onmouseout=\"popout()\">$alias</A>";
							}
							print $protString;
							print "</TD></TR>\n";
						}
						print "</TABLE>\n";
						$alignStrg='valign=middle';
					} else {$alignStrg='valign=bottom';}
					# (parent,child)=(no gel)? (item,'') : (gel sample)? (gel,spot) : (gel,gel item);
					my ($branchID,$itemBranchID)=(($item ne 'SPOT' && $item ne 'SAMPLE' && $item ne 'ANALYSIS') || $itemInfo{"$item:$itemID"}[2]{'ITEM'} eq 'SAMPLE')? (lc($item).":$itemID",'') : ($item eq 'SAMPLE')? (lc($itemInfo{"$item:$itemID"}[2]{ITEM}).":".$itemInfo{"$item:$itemID"}[2]{ID},lc($itemInfo{"$item:$itemID"}[3]{ITEM}).":".$itemInfo{"$item:$itemID"}[3]{ID}) : (lc($itemInfo{"$item:$itemID"}[2]{ITEM}).":".$itemInfo{"$item:$itemID"}[2]{ID},lc($item).":$itemID"); #Gel Sample replaced by Spot
					print "</TD><TH width=80 $alignStrg><INPUT type=\"button\" style=\"width:60px\" value=\"Open\" onclick=\"openProject($projID,'$branchID','$itemBranchID')\"></TH>\n";
				}
			}
		}
		$bgColor=($bgColor eq $color1)? $color2 : $color1;
	}
	####>No match<####
	unless ($numMatches) {
		print "<TR><TH nowrap><FONT class=\"title2\" style=\"color:#DD0000\">No match found.</FONT></TH></TR>\n";
	}

	print qq
|</TABLE><BR><BR><BR><BR>
</CENTER>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup()</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

######################################################
####>Fetching info for projects available to user<####
######################################################
my $owner=(param('OWNER'))? param('OWNER') : '';
my $selWorkgroup=(param('WG'))? param('WG') : '';

my (%projectInfo,%lastUserNames);
my $sthPj=$dbh->prepare("SELECT 1,NAME,DES,OWNER,WORK_GROUP,STATUS,LAST_OPEN,LAST_USER FROM PROJECT WHERE ID_PROJECT=?$sqlTimeStrg");
my $sthUser=$dbh->prepare("SELECT USER_NAME FROM USER_LIST WHERE ID_USER=?");
foreach my $projID (keys %{$refProfile}) {
	$sthPj->execute($projID);
	my ($match,$name,$des,$pjOwner,$workgroup,$status,$lastOpenDate,$lastUserID)=$sthPj->fetchrow_array;
	next unless $match;
	$status=0 unless $status;
	if ($display eq 'archived') {next if $status != 2;}
	else {next if $status==2;}
	next if ($owner && (($pjOwner && $pjOwner ne $owner) || (!$pjOwner && $owner !~ /no owner/i)));
	$workgroup='None' unless $workgroup;
	next if ($selWorkgroup && $selWorkgroup ne 'all' && $workgroup ne $selWorkgroup);
	$des='None' unless $des;
	my $lastOpenStrg='Never'; # default
	if ($lastOpenDate) {
		unless ($lastUserNames{$lastUserID}) {
			$sthUser->execute($lastUserID);
			($lastUserNames{$lastUserID})=$sthUser->fetchrow_array;
		}
		$lastOpenStrg=&promsMod::formatDate($lastOpenDate).' by '.$lastUserNames{$lastUserID};
	}
	@{$projectInfo{$projID}}=&promsMod::chkDef($name,$des,$pjOwner,$workgroup,$status,$lastOpenStrg);
}
$sthPj->finish;
$sthUser->finish;


########################
####>Project Window<####
########################
if ($display=~/project|archived/) {

	$dbh->disconnect;

	####>HTML<####
	my $bgColor=$color1;
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Projects List</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TD {font-size:13px;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
function openProject(projID,branchID) {
	top.promsFrame.location="$promsPath{cgi}/openProject.cgi?&ACT=open&ID="+projID+"&branchID="+branchID;
}
function restoreProject(projID) {
	if (confirm('Restore Project?')) {
		window.location="./editBranch.cgi?ITEM=PROJECT&ID="+projID+"&ACT=continue";
	}
}
</SCRIPT>
</HEAD>
<BODY topmargin=0 background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR><BR>
|;
	if ($userStatus eq 'bio') {
		print qq
|<INPUT type="button" style="font-weight=bold;font-size:20px;width:185px;" value="Main Window" onclick="window.location='./promsMain.cgi'">
<INPUT type="button" style="font-weight=bold;font-size:20px;width:185px;" value="Close Session" onclick="top.closeSession()">
<BR><BR>
|;
	}
	my $projectString='';
	if (scalar keys %projectInfo) {
		$projectString.=($owner)? 'Projects' : ($display eq 'project')? 'All active projects' : 'Archived projects';
		$projectString.=(!$selWorkgroup)? '' : ($selWorkgroup eq 'None')? ' not in a workgroup' : ($selWorkgroup eq 'all')? '' : " from Wg '$selWorkgroup'";
		$projectString.=($owner && $owner=~/no owner/i)? " with no declared owners:" : ($owner)? " owned by $owner:" : ':';
	}
	else {$projectString='No Projects available.';}
	print "<FONT class=\"title1\">$projectString</FONT>\n";
	my $width=($userStatus eq 'bio')? 600 : '95%';
	if (scalar keys %projectInfo) {
		print "<TABLE width=$width cellspacing=0>\n";
		###>Looping through all projects<###
		foreach my $projID (sort{lc($projectInfo{$a}[0]) cmp lc($projectInfo{$b}[0]) || $a<=>$b} keys %projectInfo) { # alphabetical sort on project_name
			my $access=&promsMod::getProfileAlias($$refProfile{$projID});
			my ($name,$des,$pjOwner,$workgroup,$status,$lastOpenStrg)=@{$projectInfo{$projID}};
			my $image=($status==1)? 'lightGreen1.gif' : ($status==2)? 'lightRed1.gif' : ($status==-1)? 'lightOrange1.gif' : 'lightYellow1.gif';
			print "<TR bgcolor=\"$bgColor\"><TD width=10>&nbsp</TD><TD><IMG src=\"$promsPath{images}/$image\"/>&nbsp;<FONT class=\"title2\">$name</FONT><BR><B>&nbsp;&nbsp;Description:</B> $des";
			print "<BR><B>&nbsp;&nbsp;Owner:</B> $pjOwner" unless $owner;
			print "<BR><B>&nbsp;&nbsp;Workgroup:</B> $workgroup" if (!$selWorkgroup || $selWorkgroup eq 'all');
			print "<BR><B>&nbsp;&nbsp;Last time open:</B> $lastOpenStrg";
			print "<BR><B>&nbsp;&nbsp;My access right: $access</B></TD>\n";
			if ($status==2) { # archived
				print "<TH width=110><INPUT type=\"button\" class=\"title3\" style=\"width:100px\" value=\"Restore\" onclick=\"restoreProject($projID)\"></TH>\n";
			}
			else {
				print "<TH width=80><INPUT type=\"button\" class=\"title3\" style=\"width:60px\" value=\"Open\" onclick=\"openProject($projID,'project:$projID')\"></TH>\n";
			}
			print "</TR>\n";
			$bgColor=($bgColor eq $color1)? $color2 : $color1;
		}
		print "</TABLE>\n";
	}
	print qq
|</CENTER>
</BODY>
</HTML>
|;
}

#########################
####>Analysis Window<####
#########################
elsif ($display eq 'analysis') {
	my $sort=(param('sort'))? param('sort') : 'date';
	my $sortString;
	if ($sort eq 'name') {$sortString='NAME ASC,START_DATE ASC,ID_ANALYSIS ASC';}
	elsif ($sort eq 'msType') {$sortString='MS_TYPE ASC,START_DATE ASC,NAME ASC,ID_ANALYSIS ASC';}
	elsif ($sort eq 'file') {$sortString='DATA_FILE ASC,START_DATE ASC,NAME ASC,ID_ANALYSIS ASC';}
	elsif ($sort eq 'status') {$sortString='VALID_STATUS ASC,NAME ASC,START_DATE ASC,ID_ANALYSIS ASC';}
	else {$sortString='START_DATE ASC,NAME ASC,ID_ANALYSIS ASC';}

	####>Fetching unvalidated analyses<####
	my $sthPj=$dbh->prepare("SELECT WORK_GROUP FROM PROJECT WHERE ID_PROJECT=?");
	$sqlTimeStrg=($lastOpen <= 0)? '' : "AND DATEDIFF(NOW(),START_DATE) <= $lastOpen";
	my $sthA=$dbh->prepare("SELECT ID_ANALYSIS,NAME,DES,MS_TYPE,DATA_FILE,START_DATE,VALID_STATUS FROM ANALYSIS WHERE VALID_STATUS<=1 $sqlTimeStrg ORDER BY $sortString");
	$sthA->execute;
	my @analysesData;
	my %projectWG;
	while (my @anaData=$sthA->fetchrow_array) {
		my $okAna=0;
		if ($userStatus=~/(manag|bio\Z)/) {
			my ($projID)=&promsMod::getProjectID($dbh,$anaData[0],'ANALYSIS');
			if ($userStatus eq 'manag') {
				if (!defined $projectWG{$projID}) {
					$sthPj->execute($projID);
					($projectWG{$projID})=$sthPj->fetchrow_array;
					$projectWG{$projID}='' unless $projectWG{$projID};
				}
				$okAna=1 if $userWorkgroup eq $projectWG{$projID};
			}
			unless ($okAna) {
				$okAna=1 if ($refProfile->{$projID} && $refProfile->{$projID}=~/super_/);
			}
		}
		else {$okAna=1;} # bioinfo & mass
		push @analysesData,\@anaData if $okAna;
	}
	$sthPj->finish;
	$sthA->finish;
	my $numAna=scalar @analysesData;

	####>HTML<####
	my $bgColor=$color1;
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Analyses List</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TD {font-size:13px;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
function go2Validation(anaID) {
	top.promsFrame.location="$promsPath{cgi}/startValidation.cgi?alertTax=1&ID="+anaID;
}
function openProject(projID,branchID,itemBranchID) {
	var URLstrg="$promsPath{cgi}/openProject.cgi?ACT=open&ID="+projID+"&branchID="+branchID;
	if (itemBranchID) {URLstrg+="&itemBranchID="+itemBranchID;}
	top.promsFrame.location=URLstrg;
}
</SCRIPT>
</HEAD>
<BODY topmargin="0" background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR>
<FORM name="sortForm">
<INPUT type="hidden" name="ACT" value="list">
<INPUT type="hidden" name="DISPLAY" value="analysis">
<INPUT type="hidden" name="LAST_OPEN" value="$lastOpen">
|;
	if ($numAna) {
		print qq
|<FONT class="title1">$numAna non-validated Analyses:</FONT><BR>
<FONT class="title2">Sort by :</FONT><SELECT name="sort" style="font-weight:bold;" onchange="document.sortForm.submit()">
<OPTION value="date"
|;
		print " selected" if $sort eq 'data';
		print ">Import Date</OPTION><OPTION value=\"name\"";
		print " selected" if $sort eq 'name';
		print ">Name</OPTION><OPTION value=\"msType\"";
		print " selected" if $sort eq 'msType';
		print ">MS Type</OPTION><OPTION value=\"status\"";
		print " selected" if $sort eq 'status';
		print ">Validation Status</OPTION><OPTION value=\"file\"";
		print " selected" if $sort eq 'file';
		print ">Data file</OPTION></SELECT><BR>\n";

		###>Looping through analyses<###
		print "<TABLE width=95% cellspacing=0>\n";
		my %analysisClass=(-1=>'no_scan',0=>'no_val',1=>'part_val',2=>'val');
		my %iconeImages=&promsConfig::getItemIcones;
		foreach my $refRow (@analysesData) { #@{$refAnalysisData}
			my ($analysisID,$name,$des,$msType,$file,$date,$validStatus)=@{$refRow};
			$des='None' unless $des;
			my $anaClass="analysis:$analysisClass{$validStatus}";
			$date=&promsMod::formatDate($date);
			my @anaInfo=&promsMod::getItemInfo($dbh,'analysis',$analysisID);
			my $onclickStrg;
			if ($validStatus<=0) { # go to validation
				$onclickStrg="go2Validation($analysisID)";
			}
			else { # open project
				my ($branchID,$itemBranchID)=($anaInfo[2]{'ITEM'} ne 'SAMPLE')? (lc($anaInfo[2]{ITEM}).":$anaInfo[2]{ID}","analysis:$analysisID") : ("analysis:$analysisID",'');
				$onclickStrg="openProject($anaInfo[0]{ID},'$branchID','$itemBranchID')";
			}
			print qq
|<TR bgcolor="$bgColor">
 <TD width=3></TD>
 <TD>
  <IMG src="$promsPath{images}/$iconeImages{$anaClass}">&nbsp<FONT style="font-size:18px;font-weight:bold;">$name</FONT> <B>($massSpecType{$msType})</B>
  <BR><B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Description:</B> $des<BR><B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Data file:</B> $file<BR><B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Creation date:</B> $date
  <BR><B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Project:</B> $anaInfo[0]{NAME}
 </TD>
 <TH width=80><INPUT type="button" style="width:60px" value="Open" onclick="$onclickStrg"></TH>
</TR>
|;
			$bgColor=($bgColor eq $color1)? $color2 : $color1;
		}
		print "</TABLE>\n";
		#$sthProj->finish;
	}
	else {# no analysis in list
		print "<FONT class=\"title1\"><BR><BR>No on-going Analyses found.</FONT>\n";
	}

	$dbh->disconnect;

	print qq
|</FORM>
</CENTER>
</BODY>
</HTML>
|;

}

sub createQuery {
	my ($field,$join,@words)=@_;
	for (my $i=0;$i<=$#words;$i++) {$words[$i]="UPPER($field) $words[$i]";}
	my $queryStrg=join(" $join ",@words);
	$queryStrg="($queryStrg)";
	return $queryStrg;
}

####>Revision history<####
# 2.5.0 Auto-ends analysis validation if too old & handles project status=-1 [no auto-end validation] (PP 07/06/18)
# 2.4.1 Changed default last-open time limit to 30 days (PP 17/11/17)
# 2.4.0 Open/import time filtering & simplified display (PP 11/10/17)
# 2.3.8 Minor display change (PP 22/03/16)
# 2.3.7 GPL license (PP 23/09/13)
# 2.3.6 Project status management & frames reorganization (PP 13/12/12)
# 2.3.5 Improved item auto-selection & display for manager (PP 31/07/12)
# 2.3.4 Fixes bug sequence databank conditional access (PP 25/07/12)
# 2.3.3 Conditional access to annotation for mass & manag (PP 07/05/12)
# 2.3.2 Improved workgroups/owners/projects listing (PP 17/03/12)
# 2.3.1 Add GO files managing button & species managing button (FY 17/11/11)
# 2.3.0 Bug correction in on-going analyses & search modes (PP 17/10/11)
# 2.2.9 Manager & Workgroup management (PP 18/09/11)
# 2.2.8 Super User/Admin full access to project creation/list (PP 22/08/2011)
# 2.2.7 Templates managing button activation if existing templates for a biologist user (FY 20/06/2011)
# 2.2.6 Add templates managing button (FY 02/2011)
