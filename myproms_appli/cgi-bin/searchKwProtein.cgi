#!/usr/local/bin/perl -w

################################################################################
# searchKwProtein.cgi      2.2.10                                              #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Searches DB for proteins using keywords                                      #
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

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};

#####################################
####>Arguments / Processing form<####
#####################################
my $itemID=param('ID');
my $item=param('ITEM');
my $projectID=param('id_project');
my $refItem=(param('searchItem'))? param('searchItem') : $item;
my ($refItemID,$itemChkStrg,$listChkStrg,$projChkStrg)=($refItem eq 'PROJECT')? ($projectID,'','','checked') : ($refItem eq 'LIST')? (param('list'),'','checked','') : ($itemID,'checked','','');
$itemChkStrg='checked' if $item eq 'PROJECT';
my $searchString=(param('searchString'))? param('searchString') : '';
$searchString=~s/['",\*;\(\)\[\]]//g; # ignoring these characters
my $selOrAnd='or';
my $orString='selected';
my $andString='';
if (param('selOrAnd') && param('selOrAnd') eq 'and') {
	$orString='';
	$andString='selected';
	$selOrAnd='and';
}
my $chk_partial=(!$searchString || param('partialMatch'))? 'checked=true' : ''; # checked by default
my $chk_UID=(!$searchString || param('UID'))? 'checked=true' : ''; # checked by default
my $chk_name=(!$searchString || param('name'))? 'checked=true' : ''; # checked by default
my $chk_des=(!$searchString || param('des'))? 'checked=true' : ''; # checked by default
my $chk_org=(param('org'))? 'checked=true' : '';

&storeResult if param('saveRes');


##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;
my ($itemName)=$dbh->selectrow_array("SELECT NAME FROM $item WHERE ID_$item=$itemID");
my $itemType=&promsMod::getItemType($item);
my $listName; # if a list has been selected
my ($projectStatus)=$dbh->selectrow_array("SELECT STATUS FROM PROJECT WHERE ID_PROJECT=$projectID");
$projectStatus=0 unless $projectStatus;
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
my %classificationList=&promsMod::getListClass($dbh,$projectID);
my %categoryList;
my $sthCat=$dbh->prepare("SELECT ID_CATEGORY,NAME FROM CATEGORY WHERE ID_CLASSIFICATION=? AND LIST_TYPE='PROT' ORDER BY DISPLAY_POS ASC");
foreach my $classID (keys %classificationList) {
	$sthCat->execute($classID);
	while (my ($catID,$catName)=$sthCat->fetchrow_array) {
		push @{$categoryList{$classID}},[$catID,$catName];
	}
}
$sthCat->finish;

#####################################
####>Starting HTML / Search form<####
#####################################
my ($color1,$color2)=&promsConfig::getRowColors;
print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
FONT.matched {background-color:#F6C2F9;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
|;
&promsMod::popupInfo();
print qq
|function changeSearch(item) {
	//top.promsFrame.selectedSearch='sequence';
	//window.location="./searchSqProtein.cgi?ITEM=$item&ID=$itemID&id_project=$projectID";
	if (item == 'inter') {
		window.location="./searchInteractors.cgi?ITEM=$item&ID=$itemID&id_project=$projectID";
	}
}
function checkForm(myForm) {
	if (myForm.searchItem[1].checked && !myForm.list.value) {
		alert('Select a List to be searched.');
		return false;
	}
	if (myForm.searchString.value=='') {
		alert('Type a search string.');
		return false;
	}
	if (myForm.UID.checked \|\| myForm.name.checked \|\| myForm.des.checked \|\| myForm.org.checked) {return true;} //.checked=true/false (no '' or "")
	alert('Select a search field.');
	return false;
}
function sequenceView(proteinId,analysisId){
	var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+analysisId+"&id_prot="+proteinId+"&msdata="+top.showMSdata;
	top.openProtWindow(winLocation);
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<!--<FONT class="title">Search for Proteins</FONT>-->
<FONT class="title">Search for </FONT>
<SELECT name="SearchType" onchange="changeSearch(this.value)" style="font-weight:bold;font-size:18px">
<OPTION value="prot" selected> Proteins</OPTION>
<OPTION value="inter"> Interactors</OPTION>
</SELECT>
<!-- commented sequence search option
<FONT class="title">Search for Proteins using</FONT>
<SELECT name="SearchType" onchange="changeSearch()" style="font-weight:bold;font-size:18px">
<OPTION value="kw" selected>Keywords</OPTION>
<OPTION value="sq">a Sequence</OPTION>
</SELECT>
-->
<BR><BR>
<FORM method="post" onsubmit="return checkForm(this);">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="ID" value="$itemID">
<INPUT type="hidden" name="ITEM" value="$item">
<TABLE border=0 bgcolor="$color2" cellpadding=0 cellspacing=5>
<TR>
	<TH align=right valign=top width="120px"><FONT style="font-size:18px;">Search in :</FONT></TH>
	<TH align=left><INPUT type="radio" name="searchItem" value="$item" $itemChkStrg><FONT style="font-size:18px;">$itemType <FONT style="color:#DD0000">$itemName</FONT></FONT><BR>
	<INPUT type="radio" name="searchItem" value="LIST" $listChkStrg><FONT style="font-size:18px;">List :</FONT><SELECT name="list" style="font-size:18px;"><OPTION value="">-= Select =-</OPTION>
|;
foreach my $classID (sort{lc($classificationList{$a}[0]) cmp lc($classificationList{$b}[0])} keys %classificationList) {
	print "<OPTGROUP label=\"-Theme: $classificationList{$classID}[0]\">\n";
	foreach my $refCat (@{$categoryList{$classID}}) {
		print '<OPTION value="',$refCat->[0],'"';
		if ($refItem eq 'LIST' && $refCat->[0]==$refItemID) {
			print ' selected';
			$listName=$refCat->[1];
		}
		print '>',$refCat->[1],"</OPTION>\n";
	}
	print "</OPTGROUP>\n";
}
print "</SELECT><BR>\n";
if ($item ne 'PROJECT') {
	print "<INPUT type=\"radio\" name=\"searchItem\" value=\"PROJECT\" $projChkStrg><FONT style=\"font-size:18px;\">Whole Project</FONT>";
}
print qq
|</TH>
</TR>
<TR>
	<TH align=right valign=top><FONT style="font-size:18px;">for :</FONT></TH>
	<TH align=left><INPUT type='text' name="searchString" value="$searchString" size="50"/><BR>
	<FONT style="font-size:11px;font-weight:normal;">Match is <B>not</B> case sensitive. Non-[a-z 0-9] are ignored.</FONT><BR>
	<INPUT type="checkbox" name="partialMatch" $chk_partial/>&nbsp;Partial word match.</TH>
</TR>
<TR>
	<TH align=right valign=top><FONT style="font-size:18px;">use :</FONT></TH>
	<TH align=left><SELECT name="selOrAnd" style="font-weight:bold;"><OPTION $orString value="or">or</OPTION><OPTION $andString value="and">and</OPTION></SELECT> between words.</TH>
</TR>
<TR>
	<TH align=right valign=top><FONT style="font-size:18px;">in fields :</FONT></TH>
	<TH align=left>
	<INPUT type="checkbox" name="UID" $chk_UID/>&nbsp;Protein/gene identifiers<BR>
	<INPUT type="checkbox" name="name" $chk_name/>&nbsp;Name in myProMS<BR>
	<INPUT type="checkbox" name="des" $chk_des/>&nbsp;Description<BR>
	<INPUT type="checkbox" name="org" $chk_org/>&nbsp;Organism
	</TH>
</TR>
<TR><TD colspan=2 align=center><INPUT type="submit" name="search" value="Search" style="width:100px"/></TD></TR>
</TABLE>
</FORM>
</CENTER>

<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup()</SCRIPT>
|;
# print "Search String=$searchString<BR>\n";
unless ($searchString) {
	$dbh->disconnect;
	print "</BODY>\n</HTML>\n";
	exit;
}


#----------- A search is performed ----------#
my $refItemName=($refItem eq 'PROJECT')? 'Project' : ($refItem eq 'LIST')? "List <FONT style=\"color:#DD0000\">$listName</FONT>" : "$itemType <FONT style=\"color:#DD0000\">$itemName</FONT>";
print qq
|<CENTER>
<FONT style="font-size:22px;font-weight:bold">Search results for $refItemName <SPAN id="resultCountSPAN"></SPAN></FONT>
<BR>
<DIV id="waitDIV">
<BR><BR><FONT class="title3">Searching...</FONT><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">
</DIV>
</CENTER>
|;

##################################
####>Processing search string<####
##################################
#($searchString)=($searchString=~/\A\s*(.*\S+)\s*\Z/); # removing starting and trailing spaces
$searchString=~s/^\s+//; # removing starting spaces
$searchString=~s/\s+\Z//; # removing trailing spaces
$searchString=~s/[^\w\s\.\-]+//g; # remove non-word/space except '-' or '.'
#$searchString=~s/\|/\\\|/g; # converting | into \|
# $searchString=~s/["',\.]//g; # ignoring these characters

my @wordList=split(/\s+/,lc($searchString));
my @trueWordList=@wordList; # unprocessed words
my $orPattern=join('|',@wordList); # used even for $selOrAnd='and' because each word can match different item. AND filtering is performed at result display
$orPattern='[[:<:]]('.$orPattern.')[[:>:]]' unless $chk_partial;
$orPattern="REGEXP '$orPattern'";

my @sthSrch;
if ($chk_UID) {
	push @sthSrch,$dbh->prepare("SELECT DISTINCT ID_PROTEIN FROM PROTEIN P,MASTER_PROTEIN M,MASTERPROT_IDENTIFIER MI WHERE P.ID_MASTER_PROTEIN=M.ID_MASTER_PROTEIN AND M.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND ID_PROJECT=$projectID AND (LOWER(IDENTIFIER) $orPattern OR LOWER(VALUE) $orPattern)");
}
if ($chk_name) {
	push @sthSrch,$dbh->prepare("SELECT ID_PROTEIN FROM PROTEIN P WHERE ID_PROJECT=$projectID AND (LOWER(ALIAS) $orPattern)");
}
if ($chk_des) {
	push @sthSrch,$dbh->prepare("SELECT ID_PROTEIN FROM PROTEIN P WHERE ID_PROJECT=$projectID AND (LOWER(PROT_DES) $orPattern)");
}
if ($chk_org) {
	push @sthSrch,$dbh->prepare("SELECT ID_PROTEIN FROM PROTEIN P WHERE ID_PROJECT=$projectID AND (LOWER(ORGANISM) $orPattern)");
}
my %matchedProjectProtIDs;
foreach my $sth (@sthSrch) {
	$sth->execute;
	while (my ($protID)=$sth->fetchrow_array) {
		$matchedProjectProtIDs{$protID}=0;
	}
	$sth->finish;
}

my $existMatches=0;
my @sthItem;
if ($refItem eq 'PROJECT') {
	$sthItem[0]=$dbh->prepare("SELECT P.ID_PROTEIN,ID_ANALYSIS FROM PROTEIN P,ANALYSIS_PROTEIN AP WHERE P.ID_PROTEIN=AP.ID_PROTEIN AND ID_PROJECT=$refItemID ORDER BY NUM_PEP DESC,VISIBILITY DESC,CONF_LEVEL DESC"); # best ana 1st
}
elsif ($refItem eq 'EXPERIMENT') {
	$sthItem[0]=$dbh->prepare("SELECT ID_PROTEIN,AP.ID_ANALYSIS FROM ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S WHERE AP.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_SPOT IS NULL AND ID_EXPERIMENT=$refItemID ORDER BY NUM_PEP DESC,VISIBILITY DESC,CONF_LEVEL DESC");
	$sthItem[1]=$dbh->prepare("SELECT ID_PROTEIN,AP.ID_ANALYSIS FROM ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S,SPOT SP,GEL2D G WHERE AP.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_SPOT=SP.ID_SPOT AND SP.ID_GEL2D=G.ID_GEL2D AND G.ID_EXPERIMENT=$refItemID ORDER BY NUM_PEP DESC,VISIBILITY DESC,CONF_LEVEL DESC");
}
elsif ($refItem eq 'GEL2D') {
	$sthItem[0]=$dbh->prepare("SELECT ID_PROTEIN,AP.ID_ANALYSIS FROM ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S,SPOT SP WHERE AP.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_SPOT=SP.ID_SPOT AND ID_GEL2D=$refItemID ORDER BY NUM_PEP DESC,VISIBILITY DESC,CONF_LEVEL DESC");
}
elsif ($refItem eq 'SPOT') {
	$sthItem[0]=$dbh->prepare("SELECT ID_PROTEIN,AP.ID_ANALYSIS FROM ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S WHERE AP.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND ID_SPOT=$refItemID ORDER BY NUM_PEP DESC,VISIBILITY DESC,CONF_LEVEL DESC");
}
elsif ($refItem eq 'SAMPLE') {
	$sthItem[0]=$dbh->prepare("SELECT ID_PROTEIN,AP.ID_ANALYSIS FROM ANALYSIS_PROTEIN AP,ANALYSIS A WHERE AP.ID_ANALYSIS=A.ID_ANALYSIS AND ID_SAMPLE=$refItemID ORDER BY NUM_PEP DESC,VISIBILITY DESC,CONF_LEVEL DESC");
}
elsif ($refItem eq 'ANALYSIS') {
	$sthItem[0]=$dbh->prepare("SELECT ID_PROTEIN,ID_ANALYSIS FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=$refItemID");
}
else { #LIST
	$sthItem[0]=$dbh->prepare("SELECT ID_PROTEIN,-1 FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$refItemID");
}
foreach my $sth (@sthItem) {
	$sth->execute;
	while (my ($protID,$anaID)=$sth->fetchrow_array) {
		if (defined($matchedProjectProtIDs{$protID}) && $matchedProjectProtIDs{$protID}==0) {
			$matchedProjectProtIDs{$protID}=$anaID;
			$existMatches=1;
		}
	}
	$sth->finish;
}

####>No Match<####
unless ($existMatches) {
	$dbh->disconnect;
	print qq
|<SCRIPT type="text/javascript">document.getElementById('waitDIV').style.display='none';</SCRIPT>
<FONT class="title2" style="text-align:center;color:#DD0000">No match found.</FONT>
</BODY>
</HTML>
|;
	exit;
}

####>Match<####
my ($geneIdentID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
my $sthPI=$dbh->prepare("SELECT ALIAS,ID_MASTER_PROTEIN,MW,PROT_DES,ORGANISM FROM PROTEIN WHERE ID_PROTEIN=?");
my $sthAP=$dbh->prepare("SELECT CONF_LEVEL,VISIBILITY FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=? AND ID_ANALYSIS=?");
my $sthGN=$dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER MI WHERE ID_MASTER_PROTEIN=? AND ID_IDENTIFIER=$geneIdentID ORDER BY RANK");

my (%proteinInfo,%master2genes);
foreach my $protID (keys %matchedProjectProtIDs) {
	$sthPI->execute($protID);
	@{$proteinInfo{$protID}}=$sthPI->fetchrow_array;
	@{$master2genes{$proteinInfo{$protID}[1]}}=() if $proteinInfo{$protID}[1]; # master protein ID;
	if ($refItem eq 'LIST') {
		push @{$proteinInfo{$protID}},(2,2);
	}
	else {
		$sthAP->execute($protID,$matchedProjectProtIDs{$protID}); # anaID
		push @{$proteinInfo{$protID}},$sthAP->fetchrow_array;
	}
}
foreach my $masterID (keys %master2genes) {
	$sthGN->execute($masterID);
	while (my ($geneName)=$sthGN->fetchrow_array) {
		push @{$master2genes{$masterID}},$geneName;
	}
}
$sthPI->finish;
$sthAP->finish;
$sthGN->finish;

if ($projectAccess ne 'guest' && $projectStatus <= 0) {
	print qq
|<SCRIPT LANGUAGE="JavaScript">
function changeCatList(classIdx){
	var myForm=document.storeForm;
	myForm.catList.options.length=0;
	for (i=0; i<menu[classIdx].length; i++){
		myForm.catList.options[i]=new Option(menu[classIdx][i][0],menu[classIdx][i][1]);
		if (menu[classIdx][i][2]==1) {myForm.catList.options[i].disabled=true;}
	}
	myForm.catList.selectedIndex=0;
	if (menu[classIdx].length==1) { // no category (only create option)
		document.getElementById('newCat').style.display='block';
	}
}
function selectClassification(myForm){
	var classValue=myForm.classList.value;
	if (classValue == -1) { // create new classification
		document.getElementById('newClass').style.display='block';
		document.getElementById('newCat').style.display='block';
		myForm.catList.disabled=true;
	}
	else {
		document.getElementById('newClass').style.display='none';
		document.getElementById('newCat').style.display='none';
		if (classValue == "null"){
			myForm.catList.disabled=true;
		}
		else{
			myForm.catList.disabled=false;
			window.location.href = classValue;
		}
	}
}
function checkCreateCat(optionValue){
	var IDs=optionValue.split(":");
	if (IDs[1] && IDs[1]==-1) {
		document.getElementById('newCat').style.display='block';
	}
	else {
		document.getElementById('newCat').style.display='none';
	}
}
function checkall(field,checkStatus) {
	if (field.length) { // more than 1 checkboxes
		for (var i=0; i < field.length; i++) {field[i].checked=checkStatus;}
	}
	else {
		field.checked=checkStatus;
	}
	document.getElementById('storeList').style.display=(checkStatus)? 'block' : 'none';
}
function showClassifications(boxStatus) {
	var classDisplay;
	if (boxStatus==true) {classDisplay='block';}
	else {
		if (document.storeForm.chkProt.length) {
			for (var i = 0; i < document.storeForm.chkProt.length; i++){
				if (document.storeForm.chkProt[i].checked==true) {return;}
			}
		}
		classDisplay='none';
	}
	document.getElementById('storeList').style.display=classDisplay;
}
function checkStoreForm(myForm){
	if (myForm.classList.value == 'null') { // no classification selected
		alert('ERROR: You must select a Theme and a List.');
		return false;
	}
	if (myForm.classList.value == -1) { // create new classification
		if (!myForm.newClassName.value \|\| myForm.newClassName.value=='New Theme name') {
			alert('ERROR: You must give a name to the new Theme.');
			return false;
		}
		if (!myForm.newCatName.value \|\| myForm.newCatName.value=='New List name') {
			alert('ERROR: You must give a name to the new List.');
			return false;
		}
	}
	else {
		// Existing classification selected
		var IDs=myForm.catList.options[myForm.catList.selectedIndex].value.split(":");
		if (IDs[1]=="null"){ // no category selected
			alert('ERROR: You must select a List.');
			return false;
		}
		if (IDs[1]==-1 && (!myForm.newCatName.value \|\| myForm.newCatName.value=='New List name')) { // create new category
			alert('ERROR: You must give a name to the new List.');
			return false;
		}
	}
	// Proteins checked
	var protChecked=0;
	if (myForm.chkProt) {
		if (myForm.chkProt.length) {
			for (var i=0;i<myForm.chkProt.length;i++){
				if (myForm.chkProt[i].checked) {
					protChecked=1;
					break;
				}
			}
		}
		else if (myForm.chkProt.checked) {protChecked=1;}
	}
	if (protChecked==0) {
		alert('ERROR: No proteins checked !');
		return false;
	}
	return true;
}
|;

	##################################################################
	####>Fetching list of available classification and categories<####
	##################################################################
	my (%classifications,%categories);
	my $sthClass=$dbh->prepare("SELECT ID_CLASSIFICATION,NAME FROM CLASSIFICATION WHERE ID_PROJECT=$projectID");
	my $sthCat=$dbh->prepare("SELECT ID_CATEGORY,NAME,DISPLAY_POS FROM CATEGORY WHERE ID_CLASSIFICATION=?");
	my @sthNoMod=(
		$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,COMPARISON CO WHERE CO.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=?"), # COMPARISON
		$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,CAT_COMPARISON CC WHERE CC.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=?"), # CAT_COMPARISON
		$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,GO_ANALYSIS GA WHERE GA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=?") # GO_ANALYSIS
		#$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,EXPLORANALYSIS EA WHERE EA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=?") # EXPLORANALYSIS
	);

	$sthClass->execute;
	while (my ($classID,$className)=$sthClass->fetchrow_array) {
		$classifications{$classID}=$className;
		$sthCat->execute($classID);
		while (my ($catID,$catName,$displayPos)=$sthCat->fetchrow_array) {
			@{$categories{$classID}{$catID}}=($catName,$displayPos,0);
		}
		foreach my $sth (@sthNoMod) {
			$sth->execute($classID);
			while (my ($catID)=$sth->fetchrow_array) {
				$categories{$classID}{$catID}[2]=1; # update modif flag
			}
		}
	}
	$sthClass->finish;
	$sthCat->finish;
	foreach my $sth (@sthNoMod) {$sth->finish;}

	print "var menu=new Array();\n";
	my $classIdx=0;
	print "menu[$classIdx]=new Array();\n";
	foreach my $classID (sort{lc($classifications{$a}) cmp lc($classifications{$b})} keys %classifications) {
		$classIdx++;
		print "menu[$classIdx]=new Array();\n";
		my $catIdx=0;
		if ($categories{$classID}) { # at least 1 category
			print "menu[$classIdx][$catIdx]=['-=Select a custom List=-','$classID:null',0];\n";
			foreach my $catID (sort{$categories{$classID}{$a}[1] cmp $categories{$classID}{$b}[1]} keys %{$categories{$classID}}) {
				$catIdx++;
				print "menu[$classIdx][$catIdx]=['$categories{$classID}{$catID}[0]','$classID:$catID',$categories{$classID}{$catID}[2]];\n";
			}
			$catIdx++;
		}
		print "menu[$classIdx][$catIdx]=['-=Create a new List=-','$classID:-1',0];\n";
	}
	print qq
|</SCRIPT>

<FORM name="storeForm" method="post" onsubmit="return checkStoreForm(this);">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="ID" value="$itemID">
<INPUT type="hidden" name="ITEM" value="$item">
<INPUT type="hidden" name="searchItem" value="$refItem">
<INPUT type="hidden" name="searchString" value="$searchString">
<INPUT type="hidden" name="selOrAnd" value="$selOrAnd">
<INPUT type="hidden" name="partialMatch" value="$chk_partial">
<INPUT type="hidden" name="UID" value="$chk_UID">
<INPUT type="hidden" name="name" value="$chk_name">
<INPUT type="hidden" name="des" value="$chk_des">
<INPUT type="hidden" name="org" value="$chk_org">

<TABLE border=0 width=900 bgcolor=$color2 align=center id="storeList" style="display:none"><TR>
<TH valign=top nowrap>&nbsp<FONT style="font-size:15px;">Save selected proteins in</FONT></TH>
<TD valign=top>
<SELECT name="classList" style="font-weight:bold" onchange="selectClassification(this.form)">
|;
	my $displayStatus;
	if (scalar keys %classifications) { # exists at least 1 classification
		print "<OPTION selected value=\"null\">-=Select a Theme=-</OPTION>\n";
		my $classIdx=0;
		foreach my $classID (sort{lc($classifications{$a}) cmp lc($classifications{$b})} keys %classifications) {
			$classIdx++;
			print "<OPTION value=\"Javascript:changeCatList($classIdx)\">$classifications{$classID}</OPTION>\n";
		}
		$displayStatus='none';
	}
	else {$displayStatus='block';}
	print qq
|<OPTION value="-1">-=Create a new Theme=-</OPTION>
</SELECT><BR>
<TABLE cellpadding=0 cellspacing=0 id="newClass" style="display:$displayStatus"><TR><TD>
<INPUT type="text" name="newClassName" value="New Theme name" style="width:250"/>
</TD></TR></TABLE>
</TD>
<TD valign=top>
<SELECT name="catList" style="width:250; font-weight:bold" onchange="checkCreateCat(this.value)" disabled>
<OPTION value="null"></OPTION>
</SELECT><BR>
<TABLE cellpadding=0 cellspacing=0 id="newCat" style="display:$displayStatus"><TR><TD>
<INPUT type="text" name="newCatName" value="New List name" style="width:250"/>
</TD></TR></TABLE>
</TD>
<TD><INPUT type="submit" name="saveRes" value=" Save "></TD>
<TR><TD></TD>
<TH align=left colspan=3><INPUT type="radio" name="replace" value="1" checked>Replace List contents with new proteins<BR>
<INPUT type="radio" name="replace" value="0">Add new proteins to List</TH>
</TR></TABLE>
|;
}

$dbh->disconnect;
print qq
|<SCRIPT type="text/javascript">
document.getElementById('waitDIV').style.display='none';
</SCRIPT>
|;


################################
####>Preparing result table<####
################################
print qq
|<BR>
<TABLE border=0 cellspacing=0>
<TR bgcolor="$color2">
|;
if ($projectAccess ne 'guest' && $projectStatus <= 0) {
	print qq
|<TH width=20 class="rbBorder" onmouseover="popup('Check protein boxes to display custom Lists form.')" onmouseout="popout()"><INPUT type="checkbox" onclick="checkall(document.storeForm.chkProt,this.checked)"></TH>
|;
}
else {print '<TH width=20 class="bBorder">&nbsp;</TH>';}
print qq
|<TH width=250 class="rbBorder">Protein</TH>
<TH width=150 class="rbBorder">Gene name</TH>
<TH width=75 class="rbBorder">MW <SUP>(kDa)</SUP></TH>
<TH width=900 class="bBorder">Description - Organism</TH>
</TR>
|;

my $matchStrg=join('|',@trueWordList); # used later for text highlight. '|' even if $selOrAnd='and' because results will be filtered accordingly
my $bgColor=$color1;
my $numProtFound=0;
foreach my $protID (sort{lc($proteinInfo{$a}[0]) cmp lc($proteinInfo{$b}[0])} keys %proteinInfo) { # %matchedProjectProtIDs doesn't work due to sort on %proteinInfo
	next unless $matchedProjectProtIDs{$protID};
	my ($alias,$masterID,$mw,$protDes,$organism,$bestConf,$bestVis)=@{$proteinInfo{$protID}};
	#>Furher filtering required if $selOrAnd='and'
	if ($selOrAnd eq 'and') {
		my %matchedWords;
		WORD:foreach my $word (@trueWordList) {
			if ($chk_name && $alias=~/$word/i) {$matchedWords{$word}=1; next WORD;}
			if ($chk_des && $protDes=~/$word/i) {$matchedWords{$word}=1; next WORD;}
			if ($chk_org && $organism=~/$word/i) {$matchedWords{$word}=1; next WORD;}
			if ($chk_UID && $masterID && $master2genes{$masterID}[0]) {
				foreach my $ident (@{$master2genes{$masterID}}) {
					if ($ident=~/$word/i) {$matchedWords{$word}=1; next WORD;}
				}
			}
		}
		next if scalar keys %matchedWords < scalar @trueWordList; # all words were not matched
	}
	$numProtFound++;
	#>Higlighting matched pattern
	$alias=~s/($matchStrg)/<FONT class="matched">$1<\/FONT>/ig if $chk_name;
	$protDes=~s/($matchStrg)/<FONT class="matched">$1<\/FONT>/ig if $chk_des;
	$organism=~s/($matchStrg)/<FONT class="matched">$1<\/FONT>/ig if $chk_org;

	my $bestAnaID=($matchedProjectProtIDs{$protID} < 0)? 0 : $matchedProjectProtIDs{$protID};
	print "<TR class=\"list\" bgcolor=\"$bgColor\">\n";
	my ($class,$titleString)=&promsMod::getProteinClass($bestConf,$bestVis);
	$titleString=&promsMod::HTMLcompatible($titleString);
	my $chkBoxStrg=($projectAccess eq 'guest' || $projectStatus>0)? '&nbsp;' : ($bestVis>=1)? "<INPUT type=\"checkbox\" name=\"chkProt\" value=\"$protID\" onclick=\"showClassifications(this.checked)\">" : '<INPUT type="checkbox" name="noChkProt" style="visibility:hidden">';
	print "<TH valign=top>$chkBoxStrg</TH><TD valign=top nowrap><A class=\"$class\" href=\"javascript:sequenceView($protID,$bestAnaID)\" onmouseover=\"popup('$titleString')\" onmouseout=\"popout()\">$alias</A></TD>\n";
	#>Genes name
	if ($masterID && $master2genes{$masterID}[0]) {
		if ($chk_UID) {
			$master2genes{$masterID}[0]=~s/($matchStrg)/<FONT class="matched">$1<\/FONT>/ig;
			foreach my $i (1..$#{$master2genes{$masterID}}) {
				$master2genes{$masterID}[$i]=~s/($matchStrg)/<FONT class=\\'matched\\'>$1<\/FONT>/ig; # " & ' already used in popup
			}
		}
		if (scalar @{$master2genes{$masterID}} > 1) {
			my $geneName1=shift @{$master2genes{$masterID}};
			print "<TD valign=top>&nbsp;<A class=\"$class\" href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U></B><FONT class=\\'$class\\'><BR>&nbsp;&nbsp;-",join('<BR>&nbsp;&nbsp;-',@{$master2genes{$masterID}}),"</FONT>')\" onmouseout=\"popout()\">$geneName1</A></TD>\n";
		}
		else {print "<TD valign=top class=\"$class\">&nbsp;$master2genes{$masterID}[0]</TD>\n";}
	}
	else {print "<TD class=\"$class\">&nbsp;-</TD>\n";}
	##>Mass
	print "<TD class=\"$class\" valign=top align=right>";
	if ($mw) {printf "%.1f&nbsp;&nbsp;</TH>\n",$mw/1000;}
	else {print "-&nbsp;&nbsp;</TH>\n";}
	##>Description- organism
	print "<TD>$protDes <FONT class=\"org\">$organism</FONT></TD>\n";
	print "</TR>\n";
	$bgColor=($bgColor eq $color1)? $color2 : $color1;
}

print qq
|<TR><TH colspan=5 align=left><I>&nbsp;&nbsp;End of list</I></TH></TR>
</TABLE>
|;
print "</FORM>\n" if ($projectAccess ne 'guest' && $projectStatus <= 0);
print qq
|<SCRIPT type="text/javascript">
var matchStr=($numProtFound > 1)? 'matches' : 'match';
document.getElementById('resultCountSPAN').innerHTML=' ($numProtFound '+matchStr+' found)';
</SCRIPT>
</BODY>
</HTML>
|;


##########################################
####<Storing result in Classification>####
##########################################
sub storeResult {

	####<Loading list of checked proteins>####
	my %listProteins;
	foreach my $protID (param('chkProt')) {
		$listProteins{$protID}=1;
	}

	####<Connecting to DB>####
	my $dbh=&promsConfig::dbConnect;

	my ($classID,$catID,$className,$catName);
	if (param('classList') eq '-1') {
		###>Create Classification
		$className=param('newClassName');
		($classID)=$dbh->selectrow_array("SELECT MAX(ID_CLASSIFICATION) FROM CLASSIFICATION");
		$classID++;
		my $sthInsCL=$dbh->prepare("INSERT INTO CLASSIFICATION (ID_CLASSIFICATION,NAME,ID_PROJECT,UPDATE_USER,UPDATE_DATE) VALUES (?,?,?,?,NOW())");
		$sthInsCL->execute($classID,$className,$projectID,$userID);
		$sthInsCL->finish;
		$catID=-1;
	}
	else {
		($classID,$catID)=split(/:/,param('catList'));
		($className)=$dbh->selectrow_array("SELECT NAME FROM CLASSIFICATION WHERE ID_CLASSIFICATION=$classID");
		($catName)=$dbh->selectrow_array("SELECT NAME FROM CATEGORY WHERE ID_CATEGORY=$catID");
	}
	if ($catID==-1) {
		###>Create Category
		$catName=param('newCatName');
		($catID)=$dbh->selectrow_array("SELECT MAX(ID_CATEGORY) FROM CATEGORY");
		$catID++;
		my ($displayPos)=$dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM CATEGORY WHERE ID_CLASSIFICATION=$classID");
		$displayPos++;
		my $sthInsCat=$dbh->prepare("INSERT INTO CATEGORY (ID_CATEGORY,NAME,DISPLAY_POS,ID_CLASSIFICATION,LIST_TYPE,UPDATE_USER,UPDATE_DATE) VALUES (?,?,?,?,'PROT',?,NOW())");
		$sthInsCat->execute($catID,$catName,$displayPos,$classID,$userID);
		$sthInsCat->finish;
	}
	my %oldList;
	if (param('replace')) {
		###>Delete previous category contents
		$dbh->do("DELETE FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$catID");
	}
	else { # Add to old contents
		###>Fetching existing protein contents
		my $sthOld=$dbh->prepare("SELECT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$catID");
		$sthOld->execute();
		while (my ($protID)=$sthOld->fetchrow_array) {
			$oldList{$protID}=1;
		}
		$sthOld->finish;
	}

	###>Storing protein list
	my $sthInsCP=$dbh->prepare("INSERT INTO CATEGORY_PROTEIN (ID_PROTEIN,ID_CATEGORY) VALUES (?,?)");
	foreach my $protID (keys %listProteins) {
		next if $oldList{$protID}; # protein is already in category
		$sthInsCP->execute($protID,$catID);
	}
	$sthInsCP->finish;

	###>Updating classification
	$dbh->do("UPDATE CLASSIFICATION SET UPDATE_USER='$userID',UPDATE_DATE=NOW() WHERE ID_CLASSIFICATION=$classID");

	$dbh->commit;
	$dbh->disconnect;

	####<Starting HTML>####
	print header; warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
function displayClassification() {
	top.promsFrame.selectedMode='classification:$classID';
	top.promsFrame.selectedAction='proteins';
	top.promsFrame.optionFrame.selectOption();
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR>
<FONT class="title2">Selected proteins were saved in custom List <FONT color="#DD0000">$catName</FONT> (Theme <FONT color="#DD0000">$className</FONT>).</FONT>
<BR><BR><BR><BR>
<FORM>
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="ID" value="$itemID">
<INPUT type="hidden" name="ITEM" value="$item">
<INPUT type="hidden" name="searchItem" value="$refItem">
<INPUT type="hidden" name="searchString" value="$searchString">
<INPUT type="hidden" name="selOrAnd" value="$selOrAnd">
<INPUT type="hidden" name="partialMatch" value="$chk_partial">
<INPUT type="hidden" name="UID" value="$chk_UID">
<INPUT type="hidden" name="name" value="$chk_name">
<INPUT type="hidden" name="des" value="$chk_des">
<INPUT type="hidden" name="org" value="$chk_org">

<INPUT type="button" value="Display Custom Lists" style="width:180px" onclick="displayClassification()">
&nbsp;&nbsp;&nbsp;<INPUT type="submit" value="Back to Search" style="width:180px">
</FORM>
</CENTER>
</BODY>
</HTML>
|;
	exit;
}

####>Revision history<####
# 2.2.10 Added number of matches found (PP 14/06/18)
# 2.2.9 Handles project status=-1 [no auto-end validation] (PP 07/06/18)
# 2.2.8 [Fix] bug: forgotten LOWER(ALIAS) in SQL query (PP 31/05/18)
# 2.2.7 Improved 'and' logic and partial word match handling, default search item selection & matched pattern highlighting (PP 18/10/17)
# 2.2.6 Compatible with MODIFICATION_SITE.<BR>TODO: Use promsMod/showProtQuantification functions to store proteins in Lists (PP 02/08/17)
# 2.2.5 $searchString can contain . and is not replaced anymore (GA 24/08/16)
# 2.2.4 Fix problem with REGEXP not being case insentivity (PP 17/08/15)
# 2.2.3 add select menu (protein/interactors) (SL 29/09/14)
# 2.2.2 Check if Custom lists can be modified (PP 15/05/14)
# 2.2.1 Search on List is now possible (PP 27/05/13)
# 2.2.0 Minor modification of $sth to remove ambiguity on select ID_ANALYSIS for experiment, gel2d, spot and sample (GA 16/05/13)
# 2.1.9 Uses new protein mapping & &promsMod::getProteinClass instead of getLinkClass (PP 12/02/13)
# 2.1.8 Project status management & user lists renaming (PP 19/09/12)
# 2.1.7 new protein list management options (PP 23/02/2011)
