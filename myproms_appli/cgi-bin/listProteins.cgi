#!/usr/local/bin/perl -w

################################################################################
# listProteins.cgi           2.4.7											   #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Lists the validated proteins in a project's item                             #
# List Modes:                                                                  #
# -Project hierarchy: raw,child,experiment,sample_gel,sample_spot,analysis     #
# -Classifications: #classificationID                                          #
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
$|=1;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use strict;

#print header(-'content-encoding'=>'no'); warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my %itemIcones=&promsConfig::getItemIcones;
my ($color1,$color2)=&promsConfig::getRowColors;

##############
####>Main<####
##############
# Connect to the database
my $dbh=&promsConfig::dbConnect;

####################
####>Parameters<####
####################
my $ajax=param('AJAX') || '';
if ($ajax eq 'customLists') {
	&ajaxUpdateCustomLists;
	exit;
}

my ($geneNameID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'"); # Global
my %allPostTransModifs=&promsMod::getVariableModifications($dbh);
my %itemHierarchy=('project'=>[0],'experiment'=>[1],'gel2d'=>[2],'spot'=>[3],'sample'=>[2,3],'analysis'=>[4]); #'sample_gel'=>2,'sample_spot'=>3,
my $itemID=param('ID');
my $item=lc(param('TYPE'));
my $projectID=&promsMod::getProjectID($dbh,$itemID,$item);
my $listMode=param('listMode'); $listMode='child' unless $listMode;
#my $classificationID=param('id_classification');
my $classificationID=($listMode=~/classification:(\d+)/)? $1 : 0;
my $showUnClass=param('unClass');
if ($classificationID==0 && $listMode ne 'raw' && $listMode ne 'child' && $listMode !~ /$item/) {
	#>Correcting listMode for hierarchy
	my $OKhierarchy=0;
	foreach my $it (split/_/,$listMode) {
		foreach my $pos (@{$itemHierarchy{$item}}) {
			if ($itemHierarchy{$it}>=$pos) {
				$OKhierarchy=1;
				last;
			}
		}
	}
	$listMode='child' unless $OKhierarchy;
}
##my $addRemStrg=((param('what') && param('what') eq 'addRemove'))? &addRemoveSelProt : ''; # can modifies classificationID <= commented
my $view=(param('view'))? param('view') : 'peptide';
my $expandMode=param('expMode');
my $listFilter=(param('listFilter'))? param('listFilter') : '';
my ($listFilterType,$listFilterValue)=($listFilter)? split(':',$listFilter) : ('','');
#my $ptmFilterPep=param('ptmFilterPep');
#my $extPtmFilterPep=($listFilterType eq 'ptm' && $ptmFilterPep)? 1 : 0;
my $selPepType=(param('pepType'))? param('pepType') : 'ba';
#my ($ptmString,$projectStatus)=$dbh->selectrow_array("SELECT RELEVANT_PTMS,STATUS FROM PROJECT WHERE ID_PROJECT=$projectID");
#$projectStatus=0 unless $projectStatus;
#my %projectVarMods;
#if ($ptmString) {
#	foreach my $varMod (split(';',$ptmString)) {$projectVarMods{$varMod}=$allPostTransModifs{$varMod};}
#}
my ($projectStatus)=$dbh->selectrow_array("SELECT STATUS FROM PROJECT WHERE ID_PROJECT=$projectID");
$projectStatus=0 unless $projectStatus;
my $sthGetPM=$dbh->prepare("SELECT ID_MODIFICATION FROM PROJECT_MODIFICATION WHERE ID_PROJECT=$projectID");
$sthGetPM->execute;
my %projectVarMods;
while (my ($modID)=$sthGetPM->fetchrow_array) {
	$projectVarMods{$modID}=$allPostTransModifs{$modID};
}
$sthGetPM->finish;
my ($geneNameWidth,$tableWidth)=(220,1296); #(120,1196);

#######################
####>Starting HTML<####
#######################
print header(-'content-encoding'=>'no',-charset=>'utf-8');
warningsToBrowser(1); # writes PERL warnings as HTML comments <!-- --> must be activated after header

my $ITEM=uc($item);
my ($itemName)=$dbh->selectrow_array("SELECT NAME FROM $ITEM WHERE ID_$ITEM=$itemID");
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};

my (%classificationList, %categoryList, %dbList);
&getProjectItemDBList($dbh, $projectID, \%dbList);
&getCustomList($dbh,$projectID,\%classificationList,\%categoryList);

$listMode='child' if ($classificationID && !$classificationList{$classificationID}); # in unlikely case of deletion of selected classification

print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.fakeOptGr {font-weight:bold; color:#000;}
.popPTM {font-weight:bold; background-color:$color1;}
.missing1 {text-decoration:line-through;}
.missing2 {}
|;
foreach my $varMod (sort{$a cmp $b} keys %projectVarMods) {
	print ".$projectVarMods{$varMod}[2] \{$projectVarMods{$varMod}[3]\}\n";
}
print qq
|</STYLE>
<TITLE>Proteins</TITLE>
<SCRIPT LANGUAGE="JavaScript">
|;
&promsMod::popupInfo();

&promsMod::printAjaxManageSaveProteins($projectID,\%promsPath,'document.protView.chkProt','ajaxUpdateCustomLists');

print qq
|function ajaxUpdateCustomLists(themeSel) {
	//Fetch updated Themes & List with AJAX
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/listProteins.cgi?AJAX=customLists&projectID=$projectID",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			var optList=XHR.responseText.split('\\n');

			/*** Restriction Lists ***/
			var filterSel=document.getElementById('listFilterType');
			var selectedValue=null;
			//Remove all Theme & List options
			for (var i=0; i<filterSel.options.length; i++) {
				if (filterSel.options[i].value.match('classification:')) {
					if (i <= filterSel.selectedIndex) {selectedValue=filterSel.value;} // remember selected value
					filterSel.length=i;
					break;
				}
			}
			//Append updated Themes & Lists to filterSel Classification OPTGROUP
			var clOptGr=document.getElementById('restClassifOPTGR');
			for (var i=0; i<optList.length-1; i++) { // last is null because of \\n
				var optData=optList[i].split(':=:');
				var opt;
				if (optData[0]=='T') {
					opt=new Option('-Theme '+optData[1]+':','classification:'+optData[2]);
					opt.disabled=true;
					//opt.style='font-weight:bold;color:black';
					opt.className='fakeOptGr';
				}
				else if (optData[0]=='L') {
					opt=new Option(unescape('%A0%A0%A0')+optData[1],'category:'+optData[2]); // %A0 works but not %20!
					if (optData[2]==1) opt.disabled=true;
				}
				filterSel.options.add(opt);
				clOptGr.appendChild(opt);
			}
			//Re-select preselected value
			if (selectedValue) {
				for (var i=0; i<filterSel.options.length; i++) {
					if (filterSel.options[i].value==selectedValue) {
						filterSel.selectedIndex=i;
						break;
					}
				}
			}

			/*** Display Lists ***/
			if (themeSel.value=='getLists:-1') { // new Theme
				var displaySel=document.getElementById('listMode');
				selectedValue=null;
				//Remove all Themes
				for (var i=0; i<displaySel.options.length; i++) {
					if (displaySel.options[i].value.match('classification:')) {
						if (i <= displaySel.selectedIndex) {selectedValue=displaySel.value;} // remember selected value
						displaySel.length=i;
						break;
					}
				}
				//Append updated Themes to displaySel Classification OPTGROUP
				var clOptGr=document.getElementById('dispClassifOPTGR');
				for (var i=0; i<optList.length-1; i++) { // last is null because of \\n
					var optData=optList[i].split(':=:');
					if (optData[0]=='T') {
						var opt=new Option(optData[1],'classification:'+optData[2]);
						displaySel.options.add(opt);
						clOptGr.appendChild(opt);
					}
				}
				//Re-select preselected value
				if (selectedValue) {
					for (var i=0; i<displaySel.options.length; i++) {
						if (displaySel.options[i].value==selectedValue) {
							displaySel.selectedIndex=i;
							break;
						}
					}
				}
			}
		}
	}
	XHR.send(null);
}
function extendProteinCheck(chkStatus,startPos) {
	if (!document.getElementById('serialCheck').checked) return; // No serial check
	var checkBoxList=document.protView.chkProt;
	for (var i=startPos; i<checkBoxList.length; i++) {
		checkBoxList[i].checked=chkStatus;
	}
}
function sequenceView(id_protein,id_analysis){
	var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+id_analysis+"&id_prot="+id_protein+"&msdata="+top.showMSdata;
	top.openProtWindow(winLocation);
}
function graphicalView(what,list){
	var myForm=document.protView;
	if (what=='checkbox') {
		var okChecked=false;
		if (myForm.chkProt.length) {
			for (var i=0; i<myForm.chkProt.length; i++) {
				if (myForm.chkProt[i].checked) {
					okChecked=true;
					break;
				}
			}
		}
		else {okChecked=myForm.chkProt.checked;}
		if (!okChecked) {
			alert('No protein selected!');
			return;
		}
	}
	else { // what=group
		myForm.groupInfo.value=list;
	}
	top.openProtWindow('');
	myForm.what.value=what;
	myForm.graphView.value=top.graphView;
	myForm.largWindow.value=(document.body)? top.protWindow.document.body.clientWidth : top.protWindow.innerWidth;
	myForm.action="./graphicalView.cgi";
	myForm.target='ProteinWindow';
	myForm.submit();
}
function editMatchGroup(anaId,group) {
	window.location="./editMatchGroup.cgi?id_item=$itemID&item=$item&id_analysis="+anaId+"&group="+group;
}
function changeListMode(listMode){
	top.promsFrame.selectedMode=listMode; // update promsFrame variable
	//update pepType if necessary
	if (listMode.match(/classification/) && (listMode=='analysis' \|\| '$item'=='analysis' \|\| (('$item'=='sample' \|\| '$item'=='spot') && listMode=='child'))) { // final item is 'analysis'
		top.promsFrame.selectedPepType=('$selPepType'=='all')? 'ba' : ('$selPepType'=='anr')? 'banr' : '$selPepType';
	}
	parent.optionFrame.selectOption();
}
function changeListFilter(filType) {
	var spanValue;
	var textValue;
	var textSize;
	if (!filType \|\| filType.match(':')) {
		applyListFilter(); // auto apply filter if no filter
		return;
	}
	document.getElementById('filterDiv').style.visibility='visible';
	var filterSpan=document.getElementById('filterSpan');
	var textField=document.getElementById('listFilterValue');
	if (filType=='species') {
		spanValue='&#126;';
		textSize='100px';
		filValue='';
	}
	else {
		if (!textField.value.match(/\\d/)) {filValue='';}
		else {filValue=textField.value;}
		textSize='40px';
		if (filType=='top') {spanValue=':';}
		else {spanValue='&#8805;';}
	}
	filterSpan.innerHTML=spanValue;
	textField.style.width=textSize;
	textField.value=filValue;
}
function applyListFilter() {
	var listFilType=document.getElementById('listFilterType').value;
	var listFilValue=document.getElementById('listFilterValue').value;
	if (listFilType && !listFilType.match(':') && !listFilValue) {
		alert('ERROR: Missing value for filter!');
		return;
	}
	var listFilter=(listFilType.match(':'))? listFilType : (listFilType && listFilValue)? listFilType+':'+listFilValue : '';
	top.promsFrame.selectedFilter=listFilter; // update promsFrame variable
	parent.optionFrame.selectOption();
}
function selectView(view){
	top.promsFrame.selectedView=view; // update promsFrame variable
	parent.optionFrame.selectOption();
}
function selectPeptideType(selPepType) {
	top.promsFrame.selectedPepType=selPepType; // update promsFrame variable
	parent.optionFrame.selectOption();
}
function changeMode(modeBox) {
	var expMode;
	if (modeBox.checked) {
		expMode=1;
		//update pepType if necessary
		top.promsFrame.selectedPepType=('$selPepType'=='all')? 'ba' : ('$selPepType'=='anr')? 'banr' : '$selPepType';
	}
	else {expMode=0;}
	top.expandMode=expMode; // update top variable
	parent.optionFrame.selectOption();
}
function showHideUnClass(showBox) {
	var unClass;
	if (showBox.checked) {unClass=1;}
	else {unClass=0;}
	top.showUnClass=unClass; // update top variable
	parent.optionFrame.selectOption();
}
function checkall(field,checkStatus) {
	if (field.length) { // more than 1 checkboxes
		for (var i=0; i < field.length; i++) {field[i].checked=checkStatus;}
	}
	else {
		field.checked=checkStatus;
	}
}
/*
function compareItems(){
	top.promsFrame.optionFrame.autoSelectButton('compare'); //change button selection in optionFrame
	var item1=document.protView.comparison1.value;
	var item2=document.protView.comparison2.value;
	window.location="./compareItems.cgi?FRAMES=1&id_project=$projectID&refITEM1="+item1+"&refITEM2="+item2+"&view=$view";
}
*/
function toggleMatchGroupVisibility(img,rowId) {
	var otherProtRows=document.getElementsByName('other_'+rowId);
	var newVis;
	if (img.src.match('plus')) { //otherProtRows[0].style.display=='none'
		newVis='';
		img.src='$promsPath{images}/minus1.gif';
	}
	else { //block is shown => hide block
		newVis='none';
		img.src='$promsPath{images}/plus.gif';
	}
	for (let i=0; i<otherProtRows.length; i++) {
		otherProtRows[i].style.display=newVis;
	}
}
function showMatchGroup(rowId) {
	var img=document.getElementById('img_'+rowId);
	if (img.src.match('plus')) {
		toggleMatchGroupVisibility(img,rowId);
	}
	var topProtRow=document.getElementById('top_'+rowId);
	topProtRow.scrollIntoView({block:"start",inline:"nearest",behavior:"smooth"});
	var trueBg=topProtRow.style.backgroundColor;
	topProtRow.style.backgroundColor='#FF5555'; // '#EDCA7D';
	setTimeout(function() {topProtRow.style.backgroundColor=trueBg;},2000);
}
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
// <--- AJAX
</SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<FORM name="protView" method="post">
<DIV align=center style="font-weight:bold;">
<FONT class="title">List of proteins in <FONT color=#DD0000>$itemName</FONT></FONT><BR>
<INPUT type="hidden" name="id_item" value="$itemID">
<INPUT type="hidden" name="ID" value="$itemID">
<INPUT type="hidden" name="item" value="$item">
<INPUT type="hidden" name="TYPE" value="$item">
<INPUT type="hidden" name="id_classification" value="$classificationID">
<INPUT type="hidden" name="largWindow">
<INPUT type="hidden" name="what">
<INPUT type="hidden" name="groupInfo">
<INPUT type="hidden" name="graphView">
<INPUT type="hidden" name="view" value="$view">
<INPUT type="hidden" name="expMode" value="$expandMode">
<INPUT type="hidden" name="unClass" value="$showUnClass">
<INPUT type="hidden" name="listFilter" value="$listFilter">
<INPUT type="hidden" name="pepType" value="$selPepType">
|;

my $bgColor=$color1;

print "<TABLE><TR>\n";

###>List mode
my $classPopupStrg=($classificationID && $classificationList{$classificationID}[1])? "onmouseover=\"popup('<B><U>Description:</U></B><BR>$classificationList{$classificationID}[1]')\" onmouseout=\"popout()\"" : "";
print qq
|<TH nowrap><FONT class="title2" $classPopupStrg>Display: </FONT>
<SELECT name="listMode" id="listMode" onchange="changeListMode(this.value)" style="font-weight:bold;font-size:16px;width:300px" $classPopupStrg>
<OPTGROUP label="Project hierarchy:">
|;
my @listModeOptions=(['raw','Raw list']);
push @listModeOptions,['child','Child Items'] if $item ne 'analysis';
my $OKselOption=0;
foreach my $refMode (['experiment','Experiments'],['sample_gel2d','Samples/Gels'],['sample_spot','Samples/Spots'],['analysis','Analyses']) {
	next if ($item eq 'sample' && $refMode->[0]=~/_/); # skip sample_gel2d & sample_spot
	my $it=(split/_/,$refMode->[0])[-1]; # unique or last one
	push @listModeOptions,$refMode if $itemHierarchy{$item}[0] < $itemHierarchy{$it}[0];
}
foreach my $refMode (@listModeOptions) {
	print "<OPTION value=\"$refMode->[0]\"";
	if ($refMode->[0] eq $listMode) {
		print ' selected';
		$OKselOption=1;
	}
	print ">$refMode->[1]</OPTION>\n";
}
print "</OPTGROUP>\n";
if (scalar %classificationList) {
	print "<OPTGROUP id=\"dispClassifOPTGR\" label=\"Custom lists (by Themes):\">\n";
	foreach my $classID (sort{lc($classificationList{$a}[0]) cmp lc($classificationList{$b}[0])} keys %classificationList){
		print "<OPTION value=\"classification:$classID\"";
		if ($classID==$classificationID) {
			print ' selected';
			$OKselOption=1;
		}
		print ">$classificationList{$classID}[0]</OPTION>\n";
	}
	print "</OPTGROUP>\n";
}
print "</SELECT>\n</TH>\n";
$listMode='raw' unless $OKselOption;
print "<TD width=10></TD>\n";

###>Filter
print qq
|<TH nowrap><FONT class="title2">Restrict to: </FONT>
<SELECT name="listFilterType" id="listFilterType" class="preSpace" style="font-weight:bold;font-size:16px" onchange="changeListFilter(this.value)">
<OPTION value="">All proteins</OPTION>
|;
my @filterOptions=(['pep_all','All peptides...'],
				   ['pep_anr','All distinct peptides...'],
				   ['pep_ba','All peptides in best Analysis...'],
				   ['pep_banr','Distinct pept. in best Analysis...'],
				   ['specificity','Peptide specificity (%)...'],
				   ['coverage','Peptide coverage (%)...'],
				   ['top','First n proteins...'],
				   ['species','Species...']
				   ); #,['mw','Molecular weight']
if (scalar keys %projectVarMods) {
	push @filterOptions,(['OPTGROUP','PTMs'],
						 ['ptm:none','&nbsp;&nbsp;&nbsp;No PTMs'],
						 ['ptm:any','&nbsp;&nbsp;&nbsp;Any PTMs']
						);
}
foreach my $refFilterOpt (@filterOptions) {
	if ($refFilterOpt->[0] eq 'OPTGROUP') {
		print "<OPTGROUP label=\"+$refFilterOpt->[1]:\">\n";
		next;
	}
	print "<OPTION value=\"$refFilterOpt->[0]\"";
	print ' selected' if ($refFilterOpt->[0] eq $listFilterType || $refFilterOpt->[0] eq $listFilter);
	print ">$refFilterOpt->[1]</OPTION>\n";
}
if (scalar keys %projectVarMods) {
	foreach my $varMod (sort{lc($projectVarMods{$a}[0]) cmp lc($projectVarMods{$b}[0])} keys %projectVarMods) {
		print "<OPTION value=\"ptm:$varMod\"";
		print ' selected' if ($listFilterType eq 'ptm' && $listFilterValue eq $varMod);
		print ">&nbsp;&nbsp;&nbsp;$projectVarMods{$varMod}[0]</OPTION>\n";
	}
	print "</OPTGROUP>\n";
}
if (scalar %classificationList) {
	print "<OPTGROUP id=\"restClassifOPTGR\" label=\"+Custom lists:\">\n";
	foreach my $classID (sort{lc($classificationList{$a}[0]) cmp lc($classificationList{$b}[0])} keys %classificationList){
		print "<OPTION value=\"classification:$classID\" class=\"fakeOptGr\" disabled>-Theme $classificationList{$classID}[0]:</OPTION>\n"; # style=\"font-weight:bold;color:black\"
		foreach my $catID (sort{$categoryList{$classID}{$a}[1]<=>$categoryList{$classID}{$b}[1]} keys %{$categoryList{$classID}}) {
			print "<OPTION value=\"category:$catID\"";
			print ' selected' if ($listFilter eq "category:$catID");
			print ">&nbsp;&nbsp;&nbsp;$categoryList{$classID}{$catID}[0]</OPTION>\n";
		}
	}
	print "</OPTGROUP>\n";
}

if (scalar %dbList) {
	my $foundDB = 0;
	print "<OPTGROUP id=\"restDBOPTGR\" label=\"+Databanks:\">\n";
	foreach my $dbID (keys %dbList){
		print "<OPTION value='db:$dbID'";
		if ($listFilterType eq 'db' && $listFilter eq "db:$dbID") {
			print ' selected' if ($listFilterType eq 'db' && $listFilter eq "db:$dbID");
			$foundDB = 1;
		}
		print ">&nbsp;&nbsp;&nbsp;$dbList{$dbID}</OPTION>\n";
	}
	print "</OPTGROUP>\n";
	
	if(!$foundDB && $listFilterType eq 'db') {
		$listFilterType = 'pep_all';
		$listFilter = '';
	}
}

my ($divVisibility,$usedValue)=(!$listFilter || $listFilterType =~ /ptm|category|db/)? ('hidden','') : ('visible',$listFilterValue);
my ($spanStrg,$textWidth)=($listFilterType eq 'species')? ('&#126;','100px') : ($listFilterType eq 'top')? (':','40px') : ('&#8805;','40px');
print qq
|</SELECT>
<INPUT type="button" value="Test" onclick="ajaxUpdateCustomLists()"/>
</TH>
<TD nowrap><DIV id="filterDiv" style="visibility:$divVisibility">
<SPAN id="filterSpan" class="title2">$spanStrg</SPAN>
<INPUT type="text" name="listFilterValue" id="listFilterValue" class="title3" style="width:$textWidth" value="$usedValue"/>
<INPUT type="button" class="title3" value="Apply" onclick="applyListFilter()"/></DIV></TD>
</TR></TABLE>
|;

###>Expand option
if ($listMode eq 'analysis' || ($classificationID==0 && $item eq 'analysis')) { # expand Mode
	my $expModeString=($expandMode)? 'checked' : '';
	print "<INPUT type=\"checkbox\" name=\"expMode\" $expModeString onclick=\"changeMode(this)\">&nbsp<FONT class=\"title3\">Show Match Groups</FONT>\n";
	print "<BR><FONT style=\"font-size:11px;font-style:italic\">Click on <IMG border=0 align=top src=\"$promsPath{images}/plus.gif\"> to expand Match Groups</FONT>\n" if $expandMode;
}
elsif ($classificationID > 0) { # show unclassified proteins
	my $showUnClassString=($showUnClass)? 'checked' : '';
	print "<FONT class=\"title3\"><INPUT type=\"checkbox\" name=\"unClass\" $showUnClassString onclick=\"showHideUnClass(this)\">&nbsp;Also show proteins not in Theme '$classificationList{$classificationID}[0]'</FONT>\n";
}
else {$expandMode=0;} # expand only possible in analysis
print "<BR>\n";

my $disabSave=($projectAccess eq 'guest')? ' disabled' : '';
print qq
|<TABLE><TR>
	<TD nowrap><INPUT type="checkbox" name="serialCheck" id="serialCheck" value="1"><FONT class="title3">Serial (un)check</FONT></TD>
	<TD>&nbsp;&nbsp;&nbsp;</TD>
	<TD><INPUT type="button" id="saveFormBUTTON" class="title3" value="Save Proteins..." style="display:block;" onclick="ajaxManageSaveProteins('getThemes')"$disabSave/><DIV id="saveProtDIV"></DIV></TD>
	<TD>&nbsp;&nbsp;&nbsp;</TD>
	<TD valign="middle"><INPUT type="button" class="title3" value="Peptide Distribution" onclick="graphicalView('checkbox')"></TD>
</TR></TABLE>
</DIV>
|;

##>Global variables<##
my @listChildren;
my (%listGroups,%timesFound,%listProteins,%bestAnalysis); # global because of subroutine ####---> %listClusters: obsolete!
my (%allPeptides,%nonRedundPeptides,%nonRedundPepBestAna,%anaPepVmods);
my %classProteins; # proteins in user-defined selected classification
my %catFilterProteins; # proteins in category (only if filter is a category)
my %dbFilterProteins; # proteins in a databank (only if filter is a databank)
my $itemPTMtext;
my $totNumChkboxes=0;

# %listGroups:
# @{$listGroups{$protID}{$pID}}=($numPep,$numMatch,$score°,$conf°,$vis°,$cov°,$specif°);
# %listProteins:
# $listProteins{protID}[0]: alias
# $listProteins{protID}[1]: description
# $listProteins{protID}[2]: organism
# $listProteins{protID}[3]: clusterID (obsolete!!!)
# $listProteins{protID}[4]: MW
# $listProteins{protID}[5]*: num peptides
# $listProteins{protID}[6]*: num matches
# $listProteins{protID}[7]*°: score
# $listProteins{protID}[8]*°: confidence
# $listProteins{protID}[9]*°: coverage
# $listProteins{protID}[10]*°: specificity
# * only for classification
# ° if expandMode -> properties for best analysis; if not -> best properties for item

#my $CHILD_ITEM;
if ($listFilterType eq 'category') {
	my $sthCP=$dbh->prepare("SELECT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$listFilterValue");
	$sthCP->execute;
	while (my ($protID)=$sthCP->fetchrow_array) {
		$catFilterProteins{$protID}=1;
	}
	$sthCP->finish;
}

if ($listFilterType eq 'db') {
	my @idAnalysis = ();
	getProjectItemAnaIDs(\@idAnalysis);
	
	if(@idAnalysis) {
		my $sthCP=$dbh->prepare("SELECT DISTINCT(AP.ID_PROTEIN) FROM ANALYSIS_PROTEIN AP INNER JOIN ANALYSIS_DATABANK AD ON AD.ID_ANALYSIS = AP.ID_ANALYSIS AND AD.DB_RANK = AP.DB_RANK WHERE AD.ID_DATABANK=$listFilterValue AND AP.ID_ANALYSIS IN (".join(',', @idAnalysis).")");
		$sthCP->execute;
		while (my ($protID)=$sthCP->fetchrow_array) {
			$dbFilterProteins{$protID}=1;
		}
		$sthCP->finish;
	}
}

###############################
#>User-defined classification<#
###############################
if ($classificationID) {
	##>Proteins in categories
	@listChildren=&promsMod::getChildrenList($dbh,$classificationID,'CLASSIFICATION');
	##&selectedProtOptions;
	foreach my $categoryID (sort{$listChildren[0]{'LIST'}{$a}[2] <=> $listChildren[0]{'LIST'}{$b}[2]} keys %{$listChildren[0]{'LIST'}}) {
		##print "<INPUT type=\"hidden\" name=\"childList\" value=\"CATEGORY_$categoryID\"/>\n";
		&listClassProteins($categoryID,0);
	}
	##>Unclassified proteins
	if ($showUnClass) {
		##print "<INPUT type=\"hidden\" name=\"childList\" value=\"CATEGORY_0\"/>\n";
		&listClassProteins(0,1);
	}
}
#####################
# Project hierarchy #
#####################
else {
	if ($listMode eq 'raw' || $item eq 'analysis') { # || $listMode=~/$item/
		##&selectedProtOptions;
		&listItemProteins($ITEM,$itemID);
	}
	else {
		#@listChildren=&promsMod::getChildrenList($dbh,$itemID,$item);

		####>Compare two items<####
		#my $numChildren=0;
		#my @childrenTypes;
		#foreach my $refChild (@listChildren) {
		#	$numChildren+=scalar keys %{$refChild->{'LIST'}};
		#	push @childrenTypes,&promsMod::getItemPlurial($refChild->{'ITEM'});
		#}
		#if ($numChildren > 1) { # at least 2 children
		#	print "<TABLE align=center bgcolor=$color2><TR><TD nowrap>\n";
		#	print "<FONT class=\"title2\">&nbsp;Compare two ",join(' or ',@childrenTypes)," :</FONT>\n";
		#	foreach my $compName ('comparison1','comparison2') {
		#		print '&nbsp;' if $compName eq 'comparison2';
		#		print "<SELECT name=\"$compName\" style=\"font-weight:bold;font-size:16px\">\n";
		#		my $count=0;
		#		foreach my $i (0..$#listChildren) {
		#			print "<OPTGROUP label=\"$childrenTypes[$i]\">\n" if $#listChildren >= 1;
		#			foreach my $childID (sort{$listChildren[$i]{'LIST'}{$a}[2]<=>$listChildren[$i]{'LIST'}{$b}[2]} keys %{$listChildren[$i]{'LIST'}}){
		#				$count++;
		#				#last if $count==$numChildren; # skip last child
		#				my $selStrg=($compName eq 'comparison2' && $count==2)? ' selected' : '';
		#				print "<OPTION value=\"$listChildren[$i]{ITEM}:$childID\"$selStrg>$listChildren[$i]{LIST}{$childID}[0]</OPTION>\n";
		#			}
		#			print "</OPTGROUP>\n" if $#listChildren >= 1;
		#		}
		#		print "</SELECT>\n";
		#	}
		#	print "<INPUT type=\"button\" style=\"width:50;font-size:14px\" value=\"OK\" onclick=\"javascript:compareItems()\">&nbsp\n";
		#	print "</TD></TR></TABLE>\n<BR>\n";
		#}

		####>Display protein lists<####
		##&selectedProtOptions;
		if ($listMode eq 'child') {
			@listChildren=&promsMod::getChildrenList($dbh,$itemID,$item);
			if (scalar @listChildren) {
				foreach my $i (0..$#listChildren) {
					foreach my $childID (sort{$listChildren[$i]{'LIST'}{$a}[2]<=>$listChildren[$i]{'LIST'}{$b}[2]} keys %{$listChildren[$i]{'LIST'}}) {
						##print "<INPUT type=\"hidden\" name=\"childList\" value=\"$listChildren[$i]{ITEM}_$childID\"/>\n";
						&listItemProteins($listChildren[$i]{'ITEM'},$childID,$i);
					}
				}
			}
			#else {&listItemProteins($ITEM,$itemID,-1,'<I>No child items</I>');}
			else {&listItemProteins($ITEM,$itemID,-1,'<I>No '.&promsMod::getItemType(join(' ',&promsMod::getItemChild($ITEM))).'</I>');}
		}
		else {&displayItemChildren($item,$itemID);}
	}
}
$dbh->disconnect;

print qq
|<TABLE align=center border=0 cellpadding=0 width=1020><TR><TH align=left><I>End of list.</I></TH></TR></TABLE>
</FORM>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">
top.promsFrame.selectedMode='$listMode'; // in case changed during display
setPopup();
</SCRIPT>
</BODY>
</HTML>
|;


#########################################################################
####<<<Walking down hierarchy from selected item to listMode item >>>####
#########################################################################
sub displayItemChildren {
	my ($curItem,$curItemID,$parentStrg)=@_;
	$parentStrg='' unless $parentStrg;
	my @tmpListChildren=&promsMod::getChildrenList($dbh,$curItemID,$curItem);
	if (scalar @tmpListChildren) {
		foreach my $i (0..$#tmpListChildren) {
			my $childItem=lc($tmpListChildren[$i]{'ITEM'});
			if ($listMode =~ /$childItem/) {
				@listChildren=($tmpListChildren[$i]); # only 1 index [0]
				foreach my $i (0..$#listChildren) {
					foreach my $childID (sort{$listChildren[$i]{'LIST'}{$a}[2]<=>$listChildren[$i]{'LIST'}{$b}[2]} keys %{$listChildren[$i]{'LIST'}}) {
						##print "<INPUT type=\"hidden\" name=\"childList\" value=\"$listChildren[$i]{ITEM}_$childID\"/>\n";
						&listItemProteins($listChildren[$i]{'ITEM'},$childID,$i,$parentStrg);
					}
				}
			}
			else {
				my $CHILD_ITEM=$tmpListChildren[$i]{'ITEM'}; #uc($curItem);
				my $childItem=lc($CHILD_ITEM);
				foreach my $childID (sort{$tmpListChildren[$i]{'LIST'}{$a}[2]<=>$tmpListChildren[$i]{'LIST'}{$b}[2]} keys %{$tmpListChildren[$i]{'LIST'}}) {
					my ($childName)=$dbh->selectrow_array("SELECT NAME FROM $CHILD_ITEM WHERE ID_$CHILD_ITEM=$childID");
					my $childStrg=$parentStrg;
					$childStrg.=' > ' if $childStrg;
					$childStrg.="<IMG src=\"$promsPath{images}/$itemIcones{$childItem}\"> $childName";
					&displayItemChildren($childItem,$childID,$childStrg);
				}
			}
		}
	}
	else {
		$parentStrg.=' > ' if $parentStrg;
		#&listItemProteins(uc($curItem),$curItemID,-1,"$parentStrg<I>No child items</I>");
		&listItemProteins(uc($curItem),$curItemID,-1,"$parentStrg<I>No ".&promsMod::getItemType(join(' ',&promsMod::getItemChild($curItem))).'</I>');
	}
}

###########################################################################
####<<<Fetching & listing proteins for project hierarchy or raw list>>>####
###########################################################################
sub listItemProteins { # $subjectID is for childItem except if item is ANALYSIS
	my ($subject, $subjectID, $childIdx, $parentStrg) = @_;
	my (%masterProteins, %bestNumPep, %bestScore, %bestProtMG, %protPeptides, %anaProtPeptides, %descStatus, %protVarMods, %protVarModsMG, %labelingInfo, %anaProtConfKey, %anaHasQuantifOrGO); #%anaProtVarMods,
	%listGroups = %listProteins = %timesFound = %bestAnalysis = (); # reset to () %listClusters not reset (obsolete)
	%allPeptides = %nonRedundPeptides = %nonRedundPepBestAna = %anaPepVmods = ();

	## Increase GROUP_CONCAT max length because it will contain the peptides ID list, for each protein
	$dbh->do("SET SESSION group_concat_max_len = 1000000");
	
	## Build visible proteins query
	my $visibleProtQuery = "SELECT P.ID_PROTEIN, A.ID_ANALYSIS, MATCH_GROUP, NUM_PEP, NUM_MATCH, AP.SCORE, CONF_LEVEL, VISIBILITY, PEP_COVERAGE, PEP_SPECIFICITY,
						ALIAS, PROT_DES, ORGANISM, ID_MASTER_PROTEIN, MW, 
						GROUP_CONCAT(PPA.ID_PEPTIDE) AS ID_PEPTIDE, 
						AQ.ID_QUANTIFICATION, GOA.ID_GOANALYSIS
						FROM PROTEIN P
						INNER JOIN ANALYSIS_PROTEIN AP ON P.ID_PROTEIN=AP.ID_PROTEIN
						INNER JOIN ANALYSIS A ON A.ID_ANALYSIS=AP.ID_ANALYSIS
						INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON PPA.ID_PROTEIN=P.ID_PROTEIN
						INNER JOIN PEPTIDE PEP ON PPA.ID_PEPTIDE=PEP.ID_PEPTIDE AND PPA.ID_ANALYSIS=A.ID_ANALYSIS
						LEFT JOIN ANA_QUANTIFICATION AQ ON AQ.ID_ANALYSIS=A.ID_ANALYSIS LEFT JOIN QUANTIFICATION Q ON AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND (Q.FOCUS IS NULL OR Q.FOCUS='protein')
						LEFT JOIN GOANA_ANALYSIS GOA ON GOA.ID_ANALYSIS=A.ID_ANALYSIS";
						
	my $visibleProtQueryFilter = " WHERE VISIBILITY>0 AND (CONF_LEVEL=0 || ((CONF_LEVEL=1 || CONF_LEVEL=2) && PEP.VALID_STATUS > 0))"; 
	
	if ($subject !~ /SAMPLE|ANALYSIS/) {
		$visibleProtQuery .= " INNER JOIN SAMPLE S ON A.ID_SAMPLE=S.ID_SAMPLE ";
		
		if ($subject eq 'GEL2D') {
			$visibleProtQuery .= " INNER JOIN SPOT SP ON S.ID_SPOT=SP.ID_SPOT ";
		} elsif ($subject eq 'PROJECT') {
			$visibleProtQuery .= " INNER JOIN EXPERIMENT E ON E.ID_PROJECT=P.ID_PROJECT AND S.ID_EXPERIMENT=E.ID_EXPERIMENT ";
			$visibleProtQueryFilter  .= " AND E.ID_EXPERIMENT IN (
								SELECT E2.ID_EXPERIMENT
								FROM EXPERIMENT E2
								LEFT JOIN USER_EXPERIMENT_LOCK EU ON EU.ID_EXPERIMENT=E2.ID_EXPERIMENT
								WHERE E2.ID_PROJECT=$subjectID AND NOT (EU.ID_USER IS NOT NULL AND EU.ID_USER='$userID')
							)";
		}
	}
	
	$visibleProtQueryFilter .= ($subject =~ /ANALYSIS|PROJECT/) ? " AND ".substr($subject, 0, 1).".ID_$subject=$subjectID " : "AND ID_$subject=$subjectID";
	$visibleProtQuery .= $visibleProtQueryFilter;
	$visibleProtQuery .= " GROUP BY P.ID_PROTEIN, A.ID_ANALYSIS";
	
	my $sthVis=$dbh->prepare($visibleProtQuery);
	my $sthPepMod=$dbh->prepare("SELECT PM.ID_MODIFICATION, P.ID_PEPTIDE, POS_STRING
								 FROM PEPTIDE_MODIFICATION PM 
								 INNER JOIN MODIFICATION M ON PM.ID_MODIFICATION=M.ID_MODIFICATION 
								 INNER JOIN ANALYSIS_MODIFICATION AM ON AM.ID_MODIFICATION=M.ID_MODIFICATION 
								 INNER JOIN PEPTIDE P ON PM.ID_PEPTIDE=P.ID_PEPTIDE 
								 WHERE AM.ID_ANALYSIS=? AND P.ID_ANALYSIS=? AND P.VALID_STATUS > 0");


	## Find best analysis for each visible protein 
	$sthVis->execute;
	while (my ($protID,$anaID,$matchGroup,$numPep,$numMatch,$score,$conf,$vis,$cov,$specif,$alias,$desc,$org,$masterProtID,$mw,$pepIDs,$quantiID,$goID)=$sthVis->fetchrow_array) {
		
		# Retrieve allowed PTMs (only once)
		if(!defined($anaPepVmods{$anaID})) {
			$sthPepMod->execute($anaID, $anaID);
			while (my ($modID,$pepID,$posString) = $sthPepMod->fetchrow_array) {
				$anaPepVmods{$anaID}{$pepID}{$modID}=$posString;
			}
		}
		
		next if (($listFilterType eq 'category' && !$catFilterProteins{$protID}) || ($listFilterType eq 'db' && !$dbFilterProteins{$protID})); # category filter
		$timesFound{$protID}++;
		next if ($listFilterType eq 'pep_ba' && $numPep<$listFilterValue); # type #1 (pep_ba)
		
		# Store protein's information
		$desc = ($desc) ? &promsMod::HTMLcompatible($desc) : '';
		$descStatus{$protID} = ($desc =~ /no\sdescription/) ? 2 : ($desc =~ /unnamed/) ? 1 : 0;
		$cov = 0 unless($cov);
		$org = '' unless($org);
		$mw = sprintf("%.1f", ($mw) ? $mw/1000 : 0); # MW
		@{$masterProteins{$masterProtID}} = () if($masterProtID);
		@{$listProteins{$protID}} = ($alias, $desc, $org, $masterProtID, $mw);
		
		# Check if it corresponds to the best protein
		$bestProtMG{$anaID}{$matchGroup} = $protID if($vis==2);
		if (!$bestNumPep{$protID} || $bestNumPep{$protID}<$numPep || ($bestNumPep{$protID}==$numPep && $bestScore{$protID}<$score)) {
			$bestNumPep{$protID} = $numPep;
			$bestScore{$protID} = $score;
			@{$bestAnalysis{$protID}} = ($anaID, $matchGroup);

			if ($selPepType =~ /^ba/) { # ba or banr
				@{$listGroups{$protID}{$protID}} = ($numPep, $numMatch, $score, $conf, $vis, $cov, $specif);
				undef $protVarMods{$protID} if($protVarMods{$protID}); # reset to match best ana only
			} else {
				@{$listGroups{$protID}{$protID}}[0,1] = ($numPep, $numMatch);
			}
			
			# Check if bestAna has protein quantif or GO analysis
			if(!$anaHasQuantifOrGO{$anaID} && ($goID || $quantiID)) {
				$anaHasQuantifOrGO{$anaID} = 1;
			}
		}
		$allPeptides{$protID} += $numPep;

		
		# Add PTMs to protein annotations
		my @pepIDs = split(/,/, $pepIDs);
		
		foreach my $pepID (@pepIDs) {
			my $varModStrg='';
			if ($anaPepVmods{$anaID} && $anaPepVmods{$anaID}{$pepID}) {
				foreach my $modID (sort{$a<=>$b} keys %{$anaPepVmods{$anaID}{$pepID}}) {
					$varModStrg.="+$modID:$anaPepVmods{$anaID}{$pepID}{$modID}";
					$protVarMods{$protID}{$modID}=1 if($projectVarMods{$modID});
				}
			}
			$protPeptides{$protID}{"$pepID$varModStrg"}++;
			$anaProtPeptides{$protID}{$anaID}{"$pepID$varModStrg"}++;
		}
		
		# Store best value for each property
		if ($selPepType !~ /^ba/) {
			$listGroups{$protID}{$protID}[2]=$score if (!$listGroups{$protID}{$protID}[2] || $listGroups{$protID}{$protID}[2] < $score);
			$listGroups{$protID}{$protID}[3]=$conf if (!$listGroups{$protID}{$protID}[3] || $listGroups{$protID}{$protID}[3] < $conf);
			$listGroups{$protID}{$protID}[4]=$vis if (!$listGroups{$protID}{$protID}[4] || $listGroups{$protID}{$protID}[4] < $vis);
			$listGroups{$protID}{$protID}[5]=$cov if (!$listGroups{$protID}{$protID}[5] || $listGroups{$protID}{$protID}[5] < $cov);
			$listGroups{$protID}{$protID}[6]=$specif if (!$listGroups{$protID}{$protID}[6] || $listGroups{$protID}{$protID}[6] < $specif);
		}
	}
	$sthVis->finish;
	
	## Non-redundant peptides
	foreach my $protID (keys %protPeptides) {
		$nonRedundPeptides{$protID}=scalar keys %{$protPeptides{$protID}};
		$nonRedundPepBestAna{$protID}=scalar keys %{$anaProtPeptides{$protID}{$bestAnalysis{$protID}[0]}};
	}

	## Apply filter (type #2)
	if ($listFilter && $listFilterType ne 'pep_ba' && $listFilterType !~ /species|top|category|db/) {
		foreach my $protID (keys %listGroups) {
			next if ($listFilterType eq 'pep_all' && $allPeptides{$protID}>=$listFilterValue);
			next if ($listFilterType eq 'pep_anr' && $nonRedundPeptides{$protID}>=$listFilterValue);
			next if ($listFilterType eq 'pep_banr' && $nonRedundPepBestAna{$protID}>=$listFilterValue);
			next if ($listFilterType eq 'specificity' && $listGroups{$protID}{$protID}[6]>=$listFilterValue);
			next if ($listFilterType eq 'coverage' && $listGroups{$protID}{$protID}[5]>=$listFilterValue);
			if ($listFilterType eq 'ptm') {
				if (scalar keys %{$protVarMods{$protID}}) { # prot is modified
					next if $listFilterValue eq 'any'; # any PTM
				}
				else { # prot is not modified
					next if $listFilterValue eq 'none'; # any PTM
					delete $listGroups{$protID};
					next;
				}
				next if $protVarMods{$protID}{$listFilterValue};
			}
			delete $listGroups{$protID};
		}
	}

	## Fetching list of proteins in best groups
	my %proteinGroups; # record in how many groups a protein is found (needed for filter 'top' & 'species')
	if ($expandMode && %listGroups) {
		my (%matchGroups, %idAnalysis);
		foreach my $protID (keys %listGroups) {
			#@{$listProteins{$protID}}=();
			next if($listGroups{$protID}{$protID}[4]==1); # only for best of group (vis=2) cannot be 0
			my ($idAna, $matchGroup) = @{$bestAnalysis{$protID}};
			$matchGroups{$matchGroup} = 1;
			$idAnalysis{$idAna} = 1;
		}
		
		my $sthGroup=$dbh->prepare("SELECT AP.ID_PROTEIN,
									GROUP_CONCAT(PPA.ID_PEPTIDE) AS ID_PEPTIDE,
									ALIAS, PROT_DES, ORGANISM, ID_MASTER_PROTEIN, MW, 
									NUM_PEP, NUM_MATCH, AP.SCORE, CONF_LEVEL, VISIBILITY, PEP_COVERAGE, PEP_SPECIFICITY 
									FROM ANALYSIS_PROTEIN AP
									INNER JOIN PROTEIN P ON P.ID_PROTEIN=AP.ID_PROTEIN
									INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON PPA.ID_PROTEIN=AP.ID_PROTEIN
									INNER JOIN PEPTIDE PEP ON PPA.ID_PEPTIDE=PEP.ID_PEPTIDE AND PPA.ID_ANALYSIS=AP.ID_ANALYSIS
									WHERE AP.ID_ANALYSIS IN (".join(', ', keys %idAnalysis).") AND MATCH_GROUP IN (".join(', ', keys %matchGroups).")
										AND (CONF_LEVEL=0 || ((CONF_LEVEL=1 || CONF_LEVEL=2) && PEP.VALID_STATUS > 0))
									GROUP BY AP.ID_PROTEIN, AP.ID_ANALYSIS
									ORDER BY AP.ID_ANALYSIS, MATCH_GROUP, VISIBILITY DESC");
		
		$sthGroup->execute();
		my $matchGroupHeadID;
		while (my ($pID, $pepIDs, $alias, $desc, $org, $masterProtID, $mw, @protInfo) = $sthGroup->fetchrow_array) {
			$matchGroupHeadID = $pID if($protInfo[4] == 2);
			$proteinGroups{$matchGroupHeadID}{$pID} = 1;
			next if($pID == $matchGroupHeadID); # already recorded
			@{$listGroups{$matchGroupHeadID}{$pID}} = @protInfo;
			
			# Fetching non-redondant peptides for all proteins in MG
			my %protPeptidesMG;
			my @pepIDs = split(/,/, $pepIDs);
			
			foreach my $pepID (@pepIDs) {
				my $varModStrg='';
				if ($anaPepVmods{$bestAnalysis{$matchGroupHeadID}[0]} && $anaPepVmods{$bestAnalysis{$matchGroupHeadID}[0]}{$pepID}) {
					foreach my $modID (sort{$a<=>$b} keys %{$anaPepVmods{$bestAnalysis{$matchGroupHeadID}[0]}{$pepID}}) {
						$varModStrg.="+$modID:$anaPepVmods{$bestAnalysis{$matchGroupHeadID}[0]}{$pepID}{$modID}";
						$protVarModsMG{$pID}{$modID}=1 if($projectVarMods{$modID});
					}
				}
				$protPeptidesMG{"$pepID$varModStrg"}++;
			}
			
			$desc = ($desc) ? &promsMod::HTMLcompatible($desc) : '';
			$org = '' unless($org);
			$mw = sprintf("%.1f", ($mw) ? $mw/1000 : 0); # MW
			$descStatus{$pID} = ($desc =~ /no\sdescription/) ? 2 : ($desc =~ /unnamed/) ? 1 : 0;
			@{$masterProteins{$masterProtID}} = () if($masterProtID);
			@{$listProteins{$pID}} = ($alias, $desc, $org, $masterProtID, $mw);
			
			$nonRedundPepBestAna{$pID}=scalar keys %protPeptidesMG if($selPepType eq 'banr'); # Non-redondant peptides for all proteins in MG
		}
			
		$sthGroup->finish;
	}
	$sthPepMod->finish;

	## Fetching master proteins info
	if(%masterProteins) {
		my $sthMP=$dbh->prepare("SELECT ID_MASTER_PROTEIN,VALUE FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN IN (".join(', ', keys %masterProteins).") AND ID_IDENTIFIER=$geneNameID ORDER BY RANK");
		$sthMP->execute();
		while (my ($masterProtID, $gene)=$sthMP->fetchrow_array) {
			push @{$masterProteins{$masterProtID}}, $gene;
		}
		$sthMP->finish;
	}

	## Applying filter (type #3: top or species)
	if ($listFilter && ($listFilterType eq 'top' || $listFilterType eq 'species')) {
		foreach my $protID (keys %proteinGroups) {$proteinGroups{$protID}=scalar keys %{$proteinGroups{$protID}};}
		my $count=0;
		foreach my $protID (sort {&itemSort($view)} keys %listGroups) {
			$count++;
			next if ($listFilterType eq 'top' && $count<=$listFilterValue);
			next if ($listFilterType eq 'species' && (!$listProteins{$protID}[2] || $listProteins{$protID}[2]=~/unknown/i || $listProteins{$protID}[2]=~/$listFilterValue/i));
			if ($expandMode) {
				foreach my $pID (keys %{$listGroups{$protID}}) {
					if ($proteinGroups{$pID}==1) {delete $listProteins{$pID};}
					else {$proteinGroups{$pID}--;}
				}
			}
			delete $listGroups{$protID};
		}
	}

	## Generating PTMs string
	my (%proteinPTMs,%itemPTMs);
	foreach my $protID (keys %listGroups) {
		##>Visible protein
		$proteinPTMs{$protID}{$protID}=&convertVarMods(\%{$protVarMods{$protID}},'protein');
		foreach my $varMod (keys %{$protVarMods{$protID}}) {$itemPTMs{$varMod}=1;}

		next if $listGroups{$protID}{$protID}[4]==1; # only for best of group (vis=2) cannot be 0
		##>Rest of match group
		foreach my $pID (keys %{$listGroups{$protID}}) {
			next if $pID==$protID;
			$proteinPTMs{$protID}{$pID}=&convertVarMods(\%{$protVarModsMG{$pID}},'protein');
		}
	}
	$itemPTMtext=&convertVarMods(\%itemPTMs,'item');

	## Printing data in table
	my $itemType=&promsMod::getItemType($subject);
	&prepareTable($subject,$subjectID,$itemType,$childIdx,$parentStrg,scalar keys %listGroups,scalar keys %listProteins);
	$bgColor=$color1;
	foreach my $protID (sort {&itemSort($view)} keys %listGroups) {
		my $protPosInMG=0;
		my $rowID="$subjectID"."_$protID";
		foreach my $pID (sort {$listGroups{$protID}{$b}[4]<=>$listGroups{$protID}{$a}[4] || $listGroups{$protID}{$b}[0]<=>$listGroups{$protID}{$a}[0] || $listGroups{$protID}{$b}[2]<=>$listGroups{$protID}{$a}[2] || $descStatus{$a} <=> $descStatus{$b} || $listProteins{$a}[0] cmp $listProteins{$b}[0]} keys %{$listGroups{$protID}}) { # sort by visibility>peptide>score>desc_val>alias
			my ($foundString,$expandString);
			my ($protClass,$protPopup)=&promsMod::getProteinClass($listGroups{$protID}{$pID}[3],$listGroups{$protID}{$pID}[4]); # conf,vis
			$protPosInMG++;
			if ($pID==$protID) { # top of the list visible protein
				$foundString=($timesFound{$protID}>1)? "<ACRONYM  class=\"$protClass\" onmouseover=\"popup('Protein is found in $timesFound{$protID} Analyses.')\" onmouseout=\"popout()\"><SUP>x$timesFound{$protID}</SUP></ACRONYM>" : '' ;
				if ($expandMode) {
					my $numProtMGr=scalar keys %{$listGroups{$protID}};
					if ($numProtMGr > 1){ # there are hidden proteins
						#$rowID="$subjectID"."_$protID";
						$expandString="&nbsp;<IMG border=0 align=top id=\"img_$rowID\" src=\"$promsPath{images}/plus.gif\" onclick=\"toggleMatchGroupVisibility(this,'$rowID')\" onmouseover=\"popup('<B>Best Analysis:<BR>$numProtMGr</B> proteins in <B>Match Group</B>.<BR>Click to show/hide proteins.')\" onmouseout=\"popout()\">";
					}
					else {
						if ($listGroups{$protID}{$pID}[4]==1) { # visiblity=1: protein is listed outside it's match group
# 							$expandString="&nbsp<ACRONYM class=\"$protClass\" onmouseover=\"popup('Protein was made <B>visible</B> by a user.<BR>It is listed outside its <B>Match Group</B>.')\" onmouseout=\"popout()\">";
							my $refRowID="$subjectID"."_$bestProtMG{$bestAnalysis{$protID}[0]}{$bestAnalysis{$protID}[1]}"; # {anaID}{matchGroup}
							$expandString.="&nbsp;<IMG width=16 height=22 align=top src=\"$promsPath{images}/out_group.gif\" onclick=\"showMatchGroup('$refRowID')\" onmouseover=\"popup('Protein was made <B>visible</B>.<BR>It is listed outside its <B>Match Group</B>.<BR>Click to show/hide proteins in <B>Match Group</B>.')\" onmouseout=\"popout()\">"; # </ACRONYM>";
						}
						else {
							$expandString="&nbsp;<IMG width=16 height=22 align=top src=\"$promsPath{images}/space.gif\">";
						}
					}
				}
				else {$expandString="&nbsp;";}
				# checkbox for visible protein
				$totNumChkboxes++;
				$expandString.="<INPUT type=\"checkbox\" name=\"chkProt\" value=\"$bestAnalysis{$protID}[0]:$pID\" onclick=\"extendProteinCheck(this.checked,$totNumChkboxes)\">";
			}
			else { # 2ndary proteins (must be expandMode)
				$foundString='';
				$expandString="&nbsp;<IMG width=38 height=22 align=middle src=\"$promsPath{images}/space.gif\">";
			}

			##<Begining of visible/hidden row(s)
			##if ($protPosInMG==2) { # create a new table for the 2ndary proteins
			##	print "</TABLE>\n"; # end table
			##	print "<TABLE align=center width=$tableWidth border=0 cellspacing=0 cellpadding=0 id=\"$rowID\" style=\"display:none\">\n";
			##}
			print "<TR class=\"list\" bgcolor=\"$bgColor\" valign=top";
			if ($protPosInMG==1) { # top prot
				print " id=\"top_$rowID\"";
			}
			else { # lower prot in MG (hidden)
				print " name=\"other_$rowID\" style=\"display:none\"";
			}
			print ">\n";
			#<Protein
			my $displayedName=$listProteins{$pID}[0];
			my $nameIsShorten=0;
			if (length($listProteins{$pID}[0])>22) {
				$displayedName=&promsMod::shortenName($listProteins{$pID}[0],22);
				$nameIsShorten=1;
			}
			print qq
|<TD align=left nowrap>$expandString<A class="$protClass" href="javascript:sequenceView($pID,$bestAnalysis{$protID}[0])" onmouseover="popup('<FONT class=\\'$protClass\\'>$listProteins{$pID}[0]</FONT><BR>$protPopup')" onmouseout="popout()">$displayedName$proteinPTMs{$protID}{$pID}</A>$foundString</TD>
|;
			#<Gene name
			if ($listProteins{$pID}[3] && $masterProteins{$listProteins{$pID}[3]}[0]) { # gene is not always defined even if master prot is
				if (scalar @{$masterProteins{$listProteins{$pID}[3]}} > 1) {
					my $geneName1=$masterProteins{$listProteins{$pID}[3]}[0];
					print "<TD width=$geneNameWidth align=left>&nbsp;<A class=\"$protClass\" href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U></B><FONT class=\\'$protClass\\'><BR>&nbsp;&nbsp;-",join('<BR>&nbsp;&nbsp;-',@{$masterProteins{$listProteins{$pID}[3]}}[1..$#{$masterProteins{$listProteins{$pID}[3]}}]),"</FONT>')\" onmouseout=\"popout()\">$geneName1</A></TD>";
				}
				else {print "<TD class=\"$protClass\" width=$geneNameWidth align=left>&nbsp;$masterProteins{$listProteins{$pID}[3]}[0]</TD>\n";}
			}
			else {print "<TD class=\"$protClass\" width=$geneNameWidth align=left>&nbsp;-</TD>\n";}
			#<Mass
			my $mass=($listProteins{$pID}[4]==0)? '-' : $listProteins{$pID}[4];
			print "<TD class=\"$protClass\" align=right>$mass&nbsp;&nbsp;</TD>\n";
			#<Peptides
			print "<TD class=\"$protClass\" align=center>";
			my $pepValue=($selPepType eq 'all')? $allPeptides{$pID} : ($selPepType eq 'anr')? $nonRedundPeptides{$pID} : ($selPepType eq 'banr')? $nonRedundPepBestAna{$pID} : $listGroups{$protID}{$pID}[0];
			if ($pID==$protID && ($timesFound{$protID}>1 || $listGroups{$protID}{$protID}[0] != $nonRedundPepBestAna{$protID})) {
				print "<ACRONYM onmouseover=\"popup('<B><U>Peptides:</U><BR>&nbsp;-All in $itemType:</B> $allPeptides{$protID}";
				if ($subject eq 'ANALYSIS') {
					print "<BR><B>&nbsp;-Distinct in $itemType:</B> $nonRedundPeptides{$protID}";
				}
				else {
					print "<BR><B>&nbsp;-Distinct in $itemType:</B> $nonRedundPeptides{$protID}<BR><B>&nbsp;-All in best Analysis:</B> $listGroups{$protID}{$pID}[0]<BR><B>&nbsp;-Distinct in best Analysis :</B> $nonRedundPepBestAna{$protID}";
				}
				print "')\" onmouseout=\"popout()\"><U>$pepValue</U></ACRONYM>";
			}
			else {print $pepValue;}
			print "<ACRONYM class=\"$protClass\" onmouseover=\"popup('$listGroups{$protID}{$pID}[1] matches due to sequence repeats.')\" onmouseout=\"popout()\"><SUP>($listGroups{$protID}{$pID}[1])</SUP></ACRONYM>" if $listGroups{$protID}{$pID}[1]>$listGroups{$protID}{$pID}[0];
			print "</TD>\n";
# 			#<Score
# 			print "<TD width=70 align=center>$listGroups{$protID}{$pID}[2]</TD>\n";

			#<Specificity
			print "<TD class=\"$protClass\" align=center>$listGroups{$protID}{$pID}[6]</TD>\n";
			#<Coverage
			my $covText=(!$listGroups{$protID}{$pID}[5])? '?' : ($listGroups{$protID}{$pID}[5] > 0)? $listGroups{$protID}{$pID}[5] : '&lt;'.abs($listGroups{$protID}{$pID}[5]);
			print "<TD class=\"$protClass\" align=center>$covText</TD>\n";

			#<Description - Species
			print "<TD>$listProteins{$pID}[1] <FONT class=\"org\">$listProteins{$pID}[2]</FONT></TD>\n";
	    	print "</TR>\n";
    	} # end of pID loop

    	##<End of visible/hidden row(s)
    	if ($protPosInMG>1 && $protPosInMG == scalar keys %{$listGroups{$protID}}) {
	    	# Edit & Peptide Distribution buttons
	    	print "<TR bgcolor=\"$bgColor\" name=\"other_$rowID\" style=\"display:none\"><TD colspan=7>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\n";
	    	print "<INPUT type=\"button\" class=\"font11\" value=\"Peptide Distribution\" onclick=\"graphicalView('group','$bestAnalysis{$protID}[0]:$bestAnalysis{$protID}[1]')\">\n";
	    	if ($projectAccess ne 'guest' && $projectStatus <= 0) {
		    	print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;";
				my $disabStrg=($anaHasQuantifOrGO{$bestAnalysis{$protID}[0]})? ' disabled' : '';
		    	print "<INPUT type=\"button\" class=\"font11\"  value=\"Edit Match Group\" onclick=\"editMatchGroup($bestAnalysis{$protID}[0],$bestAnalysis{$protID}[1])\"$disabStrg>\n";
	    	}
			print "</TD></TR>\n";
	    	##print "</TABLE>\n"; # end visible/hidden table
	    	##print "<TABLE align=center border=0 cellpadding=0 cellspacing=0 width=$tableWidth>\n"; # start next table
    	}
    	$bgColor=($bgColor eq $color1)? $color2 : $color1;
	} # end of protID loop
	print "</TABLE><BR>\n";
}

##########################################################################
####<<<Fetching & listing proteins for user-defined classifications>>>####
##########################################################################
sub listClassProteins {
	my ($categoryID, $filter)=@_;
	%listProteins = %timesFound = %bestAnalysis = (); # reset to () %listClusters not reset (obsolete)
	%allPeptides = %nonRedundPeptides = %nonRedundPepBestAna = %anaPepVmods = ();
	my $subject = uc($item);
	
	## Increase GROUP_CONCAT max length because it will contain the peptides ID list, for each protein
	$dbh->do("SET SESSION group_concat_max_len = 1000000");
	
	## Build visible proteins query
	my $visibleProtQuery = "SELECT P.ID_PROTEIN, A.ID_ANALYSIS, NUM_PEP, NUM_MATCH, AP.SCORE, CONF_LEVEL, PEP_COVERAGE, PEP_SPECIFICITY,
						ALIAS, PROT_DES, ORGANISM, ID_MASTER_PROTEIN, MW, 
						GROUP_CONCAT(PPA.ID_PEPTIDE) AS ID_PEPTIDE 
						FROM PROTEIN P 
						INNER JOIN ANALYSIS_PROTEIN AP ON P.ID_PROTEIN=AP.ID_PROTEIN 
						INNER JOIN ANALYSIS A ON A.ID_ANALYSIS=AP.ID_ANALYSIS 
						INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON PPA.ID_PROTEIN=P.ID_PROTEIN 
						INNER JOIN PEPTIDE PEP ON PPA.ID_PEPTIDE=PEP.ID_PEPTIDE AND PPA.ID_ANALYSIS=A.ID_ANALYSIS";
						
	my $visibleProtQueryFilter = " WHERE VISIBILITY>0 AND (CONF_LEVEL=0 || ((CONF_LEVEL=1 || CONF_LEVEL=2) && PEP.VALID_STATUS > 0))"; 
	
	if ($subject !~ /SAMPLE|ANALYSIS/) {
		$visibleProtQuery .= " INNER JOIN SAMPLE S ON A.ID_SAMPLE=S.ID_SAMPLE ";
		
		if ($subject eq 'GEL2D') {
			$visibleProtQuery .= " INNER JOIN SPOT SP ON S.ID_SPOT=SP.ID_SPOT ";
		} elsif ($subject eq 'PROJECT') {
			$visibleProtQuery .= " INNER JOIN EXPERIMENT E ON E.ID_PROJECT=P.ID_PROJECT AND S.ID_EXPERIMENT=E.ID_EXPERIMENT ";
			$visibleProtQueryFilter  .= " AND E.ID_EXPERIMENT IN (
								SELECT E2.ID_EXPERIMENT
								FROM EXPERIMENT E2
								LEFT JOIN USER_EXPERIMENT_LOCK EU ON EU.ID_EXPERIMENT=E2.ID_EXPERIMENT
								WHERE E2.ID_PROJECT=$itemID AND NOT (EU.ID_USER IS NOT NULL AND EU.ID_USER='$userID')
							)";
		}
	}
	
	$visibleProtQueryFilter .= ($subject =~ /ANALYSIS|PROJECT/) ? " AND ".substr($subject, 0, 1).".ID_$subject=$itemID " : "AND ID_$subject=$itemID";
	$visibleProtQuery .= $visibleProtQueryFilter;
	$visibleProtQuery .= " GROUP BY P.ID_PROTEIN, A.ID_ANALYSIS";
	
	my $sthVis=$dbh->prepare($visibleProtQuery);
	my $sthPepMod=$dbh->prepare("SELECT PM.ID_MODIFICATION,P.ID_PEPTIDE,POS_STRING FROM PEPTIDE_MODIFICATION PM, MODIFICATION M, ANALYSIS_MODIFICATION AM, PEPTIDE P WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND PM.ID_MODIFICATION=M.ID_MODIFICATION AND PM.ID_PEPTIDE=P.ID_PEPTIDE AND AM.ID_ANALYSIS=? AND P.ID_ANALYSIS=? AND P.VALID_STATUS > 0");

	####<List of proteins in current Category>####
	my %catProteins;
	if (!$filter) { # proteins in classification
		my $sthCP=$dbh->prepare("SELECT DISTINCT CP.ID_PROTEIN, P.ALIAS, NUM_PEP, PEP_COVERAGE, PEP_SPECIFICITY, PROT_DES, MW, ORGANISM
								FROM CATEGORY_PROTEIN CP
								INNER JOIN PROTEIN P ON CP.ID_PROTEIN=P.ID_PROTEIN
								INNER JOIN ANALYSIS_PROTEIN AP ON P.ID_PROTEIN=AP.ID_PROTEIN
								INNER JOIN ANALYSIS A ON A.ID_ANALYSIS=AP.ID_ANALYSIS
								INNER JOIN SAMPLE S ON S.ID_SAMPLE=A.ID_SAMPLE
								INNER JOIN EXPERIMENT E ON E.ID_PROJECT=P.ID_PROJECT AND S.ID_EXPERIMENT=E.ID_EXPERIMENT
								WHERE ID_CATEGORY=$categoryID AND P.ID_PROJECT=$projectID AND E.ID_EXPERIMENT IN (
									SELECT E2.ID_EXPERIMENT
									FROM EXPERIMENT E2
									LEFT JOIN USER_EXPERIMENT_LOCK EU ON EU.ID_EXPERIMENT=E2.ID_EXPERIMENT
									WHERE E2.ID_PROJECT=$projectID AND NOT (EU.ID_USER IS NOT NULL AND EU.ID_USER='$userID')
								)");
		$sthCP->execute;
		while (my ($protID, $alias, $numPep, $pepCoverage, $pepSpe, $desc, $mw, $org)=$sthCP->fetchrow_array) {
			next if ($listFilterType eq 'category' && !$catFilterProteins{$protID});
			$catProteins{$protID}=1;
			$classProteins{$protID}=1; # global
			@{$listProteins{$protID}}=($alias,$desc,$org,'','','','','','','',''); # default
		}
		$sthCP->finish;
	}
		
	####<List of visible proteins in selected Project item>####
	my (%bestNumPep,%bestScore,%protPeptides,%anaProtPeptides,%protVarMods,%masterProteins); #,%anaProtVarMods
	$sthVis->execute;
	while (my ($protID,$anaID,$numPep,$numMatch,$score,$conf,$cov,$specif,$alias,$desc,$org,$masterProtID,$mw,$pepIDs)=$sthVis->fetchrow_array) {
		next if ((!$filter && !$catProteins{$protID}) || ($filter && $classProteins{$protID})); # next if (classProt wanted && prot not a cat) || (classProt not wanted && prot in any cat))
		next if ($listFilterType eq 'category' && !$catFilterProteins{$protID});
		$timesFound{$protID}++;
		next if ($listFilterType eq 'pep_ba' && $numPep<$listFilterValue);
		
		# Retrieve allowed PTMs (only once)
		if(!defined($anaPepVmods{$anaID})) {
			$sthPepMod->execute($anaID, $anaID);
			while (my ($modID,$pepID,$posString) = $sthPepMod->fetchrow_array) {
				$anaPepVmods{$anaID}{$pepID}{$modID}=$posString;
			}
		}
		
		# Store protein's information
		$desc = ($desc) ? &promsMod::HTMLcompatible($desc) : '';
		$cov = 0 unless($cov);
		$org = '' unless($org);
		$mw = sprintf("%.1f", ($mw) ? $mw/1000 : 0); # MW
		@{$masterProteins{$masterProtID}} = () if($masterProtID);
		$allPeptides{$protID} += $numPep;
		@{$listProteins{$protID}}[0,1,2,3,4] = ($alias, $desc, $org, $masterProtID, $mw);
		
		# Check if it corresponds to the best protein
		if (!$bestNumPep{$protID} || $bestNumPep{$protID}<$numPep || ($bestNumPep{$protID}==$numPep && $bestScore{$protID}<$score)) {
			$bestNumPep{$protID} = $numPep;
			$bestScore{$protID} = $score;
			$bestAnalysis{$protID} = $anaID;

			if ($selPepType =~ /^ba/) { # ba or banr
				@{$listProteins{$protID}}[5,6,7,8,9,10] = ($numPep, $numMatch, $score, $conf, $cov, $specif);
				undef $protVarMods{$protID} if($protVarMods{$protID}); # reset to match best ana only
			} else {
				@{$listProteins{$protID}}[5,6] = ($numPep, $numMatch);
			}
		}
		
		# Add PTMs to protein annotations
		my @pepIDs = split(/,/, $pepIDs);
		
		foreach my $pepID (@pepIDs) {
			my $varModStrg='';
			if ($anaPepVmods{$anaID} && $anaPepVmods{$anaID}{$pepID}) {
				foreach my $modID (sort{$a<=>$b} keys %{$anaPepVmods{$anaID}{$pepID}}) {
					$varModStrg.="+$modID:$anaPepVmods{$anaID}{$pepID}{$modID}";
					$protVarMods{$protID}{$modID}=1 if($projectVarMods{$modID});
				}
			}
			$protPeptides{$protID}{"$pepID$varModStrg"}++;
			$anaProtPeptides{$protID}{$anaID}{"$pepID$varModStrg"}++;
		}
		
		# Store best value for each property
		if ($selPepType !~ /^ba/) {
			$listProteins{$protID}[7] = $score if(!$listProteins{$protID}[7] || $listProteins{$protID}[7] < $score);
			$listProteins{$protID}[8] = $conf if(!$listProteins{$protID}[8] || $listProteins{$protID}[8] < $conf);
			$listProteins{$protID}[9] = $cov if(!$listProteins{$protID}[9] || $listProteins{$protID}[9] < $cov);
			$listProteins{$protID}[10] = $specif if(!$listProteins{$protID}[10] || $listProteins{$protID}[10] < $specif);
		}
	}
	$sthVis->finish;
	$sthPepMod->finish;

	## Non-redundant peptides
	foreach my $protID (keys %protPeptides) {
		$nonRedundPeptides{$protID} = scalar keys %{$protPeptides{$protID}};
		$nonRedundPepBestAna{$protID} = scalar keys %{$anaProtPeptides{$protID}{$bestAnalysis{$protID}}};
	}

	## Fetching master proteins info
	if(%masterProteins) {
		my $sthMP=$dbh->prepare("SELECT ID_MASTER_PROTEIN,VALUE FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN IN (".join(', ', keys %masterProteins).") AND ID_IDENTIFIER=$geneNameID ORDER BY RANK");
		$sthMP->execute();
		while (my ($masterProtID, $gene)=$sthMP->fetchrow_array) {
			push @{$masterProteins{$masterProtID}}, $gene;
		}
		$sthMP->finish;
	}
	
	## Applying filter (type 2 & 3)
	if ($listFilter && $listFilterType !~ /pep_ba|category|db/) {
		my $count=0;
		foreach my $protID (sort{&classSort($view)} keys %listProteins) {
			$count++;
			next if ($listFilterType eq 'pep_all' && $allPeptides{$protID}>=$listFilterValue);
			next if ($listFilterType eq 'pep_anr' && $nonRedundPeptides{$protID}>=$listFilterValue);
			next if ($listFilterType eq 'pep_banr' && $nonRedundPepBestAna{$protID}>=$listFilterValue);
			next if ($listFilterType eq 'specificity' && $listProteins{$protID}[10]>=$listFilterValue);
			next if ($listFilterType eq 'coverage' && $listProteins{$protID}[9]>=$listFilterValue);
			next if ($listFilterType eq 'top' && $count<=$listFilterValue);
			next if ($listFilterType eq 'species' && (!$listProteins{$protID}[2] || $listProteins{$protID}[2]=~/unknown/i || $listProteins{$protID}[2]=~/$listFilterValue/i));
			if ($listFilterType eq 'ptm') {
				if (scalar keys %{$protVarMods{$protID}}) { # prot is modified
					next if $listFilterValue eq 'any'; # any PTM
				}
				else { # prot is not modified
					next if $listFilterValue eq 'none'; # any PTM
					delete $listProteins{$protID};
					next;
				}
				next if $protVarMods{$protID}{$listFilterValue};
			}
			delete $listProteins{$protID};
		}
	}

	
	## Generating PTMs string
	my (%proteinPTMs, %itemPTMs);
	foreach my $protID (keys %listProteins) {
		$proteinPTMs{$protID} = &convertVarMods(\%{$protVarMods{$protID}},'protein');
		foreach my $varMod (keys %{$protVarMods{$protID}}) { $itemPTMs{$varMod}=1; }
	}
	$itemPTMtext=&convertVarMods(\%itemPTMs,'item');
	
	
	## Printing data in table
	my $itemType=&promsMod::getItemType($item);
	&prepareTable('CATEGORY',$categoryID,$itemType,0,'',(scalar keys %allPeptides),(scalar keys %listProteins),$filter);
	$bgColor=$color1;
	foreach my $protID (sort{&classSort($view)} keys %listProteins) {
		print "<TR $geneNameWidth class=\"list\" bgcolor=\"$bgColor\" valign=top>\n";
		#<Protein
		print "<TD width=250>&nbsp;";
		my $displayedName=$listProteins{$protID}[0];
		my $nameIsShorten=0;
		if (length($listProteins{$protID}[0])>22) {
			$displayedName=&promsMod::shortenName($listProteins{$protID}[0],22);
			$nameIsShorten=1;
		}
		my $protClass;
		if ($listProteins{$protID}[6]) { # Found in project Item
			($protClass,my $protPopup)=&promsMod::getProteinClass($listProteins{$protID}[8],1); # conf,vis >= 1
			my $foundString=($timesFound{$protID}>1)? "<ACRONYM class=\"$protClass\" onmouseover=\"popup('Protein is found in $timesFound{$protID} Analyses.')\" onmouseout=\"popout()\"><SUP>x$timesFound{$protID}</SUP></ACRONYM>" : '' ;
			$totNumChkboxes++;
			print qq
|<INPUT type="checkbox" name="chkProt" value="$bestAnalysis{$protID}:$protID" onclick="extendProteinCheck(this.checked,$totNumChkboxes)"><A class="$protClass" href="javascript:sequenceView($protID,$bestAnalysis{$protID})" onmouseover="popup('<FONT class=\\'$protClass\\'>$listProteins{$protID}[0]</FONT><BR>$protPopup')" onmouseout="popout()">$displayedName$proteinPTMs{$protID}</A>$foundString|;
		}
		else {
			$protClass='missing2';
			print qq
|<INPUT type="checkbox" name="chkProt_bad" value="0:$protID" disabled><A class="missing1" href="javascript:sequenceView($protID,0)" onmouseover="popup('<FONT class=\\'$protClass\\'>$listProteins{$protID}[0]</FONT><BR>Not found/visible in selected Project item')" onmouseout="popout()">$displayedName</A>|;
		}
		print "</TD>\n";
		#<Gene name
		if ($listProteins{$protID}[3] && $masterProteins{$listProteins{$protID}[3]}[0]) { # gene is not always defined even if master prot is
			if (scalar @{$masterProteins{$listProteins{$protID}[3]}} > 1) {
				my $geneName1=shift(@{$masterProteins{$listProteins{$protID}[3]}});
				print "<TD width=$geneNameWidth align=left>&nbsp;<A class=\"$protClass\" href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U></B><FONT class=\\'$protClass\\'><BR>&nbsp;&nbsp;-",join('<BR>&nbsp;&nbsp;-',@{$masterProteins{$listProteins{$protID}[3]}}),"</FONT>')\" onmouseout=\"popout()\">$geneName1</A></TD>";
			}
			else {print "<TD class=\"$protClass\" width=$geneNameWidth align=left>&nbsp;$masterProteins{$listProteins{$protID}[3]}[0]</TD>\n";}
		}
		else {print "<TD class=\"$protClass\" width=$geneNameWidth align=left>&nbsp;-</TD>\n";}
		#<Mass
		my $mass=($listProteins{$protID}[4]==0)? '-' : $listProteins{$protID}[4];
		print "<TD class=\"$protClass\" width=72 align=right>$mass&nbsp;&nbsp;</TD>\n";
		#<Peptides (index 5 = num peptides, index 6 = num matches)
		if ($listProteins{$protID}[6]) { # found in project item
			print "<TD class=\"$protClass\" width=60 align=center>";
			my $pepValue=($selPepType eq 'all')? $allPeptides{$protID} : ($selPepType eq 'anr')? $nonRedundPeptides{$protID} : ($selPepType eq 'banr')? $nonRedundPepBestAna{$protID} : $listProteins{$protID}[5];
			if ($timesFound{$protID}>1 ||  $listProteins{$protID}[5] != $nonRedundPepBestAna{$protID}) {
				print "<ACRONYM onmouseover=\"popup('<B><U>Peptides in $itemType:</U><BR>&nbsp;-All:</B> $allPeptides{$protID}";
				if ($item eq 'analysis') {
					print "<BR><B>&nbsp;-Distinct:</B> $nonRedundPeptides{$protID}";
				}
				else {
					print "<BR><B>&nbsp;-All distinct:</B> $nonRedundPeptides{$protID}<BR><B>&nbsp;-Best Analysis:</B> $listProteins{$protID}[5]<BR><B>&nbsp;-Best Analysis distinct:</B> $nonRedundPepBestAna{$protID}";
				}
				print "')\" onmouseout=\"popout()\"><U>$pepValue</U></ACRONYM>";
			}
			else {print $pepValue;}
			print "<ACRONYM onmouseover=\"popup('$listProteins{$protID}[6] matches due to sequence repeats.')\" onmouseout=\"popout()\"><SUP>($listProteins{$protID}[6])</SUP></ACRONYM>" if $listProteins{$protID}[6]>$listProteins{$protID}[5];
			print "</TD>\n";
# 			#<Score
# 			print "<TH width=70 align=center>$listProteins{$protID}[7]</TH>\n";

			#<Specificity
			print "<TD class=\"$protClass\" width=65 align=center>$listProteins{$protID}[10]</TD>\n";
			#<Coverage
			my $covText=(!$listProteins{$protID}[9])? '?' : ($listProteins{$protID}[9] > 0)? $listProteins{$protID}[9] : '&lt;'.abs($listProteins{$protID}[9]);
			print "<TD class=\"$protClass\" width=65 align=center>$covText</TD>\n";
		}
		else { # no in project item
			print qq
|<TD class="visible" width=60 align=center>-</TD><TD class="$protClass" width=60 align=center>-</TD><TD class="$protClass" width=60 align=center>-</TD>
|;
		}

		#<Description - Species
		print "<TD width=575>$listProteins{$protID}[1] <FONT class=\"org\">$listProteins{$protID}[2]</FONT></TD>\n";
		print "</TR>\n";
		$bgColor=($bgColor eq $color1)? $color2 : $color1;
	}
	print "</TABLE><BR>\n";
}

###########################
####<<<Prepare table>>>####
###########################
sub prepareTable {
	my ($subject,$subjectID,$itemType,$childIdx,$parentStrg,$minListSize,$maxListSize,$filter)=@_; #listSize

	print "<TABLE align=center border=0 cellpadding=0 cellspacing=0 width=$tableWidth>\n<TR><TD colspan=7 class=\"bBorder\">\n";
	if (!$filter) {
		#if ($classificationID>0 || ($item ne 'analysis' && $classificationID>-1)) { #}
		if (defined($childIdx)) {
			if ($childIdx>=0) {
				if ($parentStrg) {$parentStrg.=' > ';} else {$parentStrg='';}
				my $itemCode='';
				if ($listChildren[$childIdx]{'ITEM'} eq 'ANALYSIS') {
					$itemCode=($listChildren[$childIdx]{'LIST'}{$subjectID}[3]==-1)? ':no_scan' : ($listChildren[$childIdx]{'LIST'}{$subjectID}[3]==0)? ':no_val' : ($listChildren[$childIdx]{'LIST'}{$subjectID}[3]==1)? ':part_val' : ':val';
				}
				#my $imgIcon=$itemIcones{lc($listChildren[$childIdx]{'ITEM'})."$itemCode"};
				print "<FONT class=\"title3\">$parentStrg<IMG src=\"$promsPath{images}/",$itemIcones{lc($listChildren[$childIdx]{'ITEM'})."$itemCode"},"\"> $listChildren[$childIdx]{LIST}{$subjectID}[0] :<FONT style=\"font-weight:normal\"> $listChildren[$childIdx]{LIST}{$subjectID}[1]</FONT></FONT><BR>\n";
			}
			elsif ($parentStrg) { # childIdx=-1 (no child items)
				print "<FONT class=\"title3\">$parentStrg</FONT><BR>";
			}
		}
	}
	else {print "<FONT class=\"title2\">Proteins not in Theme '$classificationList{$classificationID}[0]':</FONT><BR>\n";}

	print "</TD></TR>\n<TR bgcolor=$color2>\n";

	# Printing according to selected sort
	my $sStrg=($maxListSize>1)? 's' : '';
	my $proteinValue=($minListSize==$maxListSize)? "$minListSize Protein$sStrg" : "$minListSize/$maxListSize Protein$sStrg";
	my $massValue='MW<FONT style="font-size:9px;"> kDa</FONT>';
	#my $peptideValue='Pept.';
	my $peptideValue='P';
	my ($selAP,$selAD,$selBA,$selBD)=($selPepType eq 'all')? (' selected','','','') : ($selPepType eq 'anr')? ('',' selected','','') : ($selPepType eq 'banr')? ('','','',' selected') : ('','',' selected',''); # default is best Ana
	my $analysisStrg;
	my $pepSelectStrg='<SELECT name="selPepType" style="width:42px" onchange="selectPeptideType(this.value)">';
	if ($item ne 'analysis' && $subject ne 'ANALYSIS') {
		$pepSelectStrg.="<OPTION value=\"all\"$selAP>AP: All Peptides in $itemType</OPTION><OPTION value=\"anr\"$selAD>DP: Distinct pept. in $itemType</OPTION>"  if (!$expandMode || $classificationID>0);
		$pepSelectStrg.="<OPTION value=\"ba\"$selBA>AB: All peptides in Best analysis</OPTION><OPTION value=\"banr\"$selBD>DB: Distinct pept. in Best analysis</OPTION>";
		$analysisStrg=($selPepType=~/^ba/)? 'From best Analysis' : 'Best value';
	}
	else {
		$pepSelectStrg.="<OPTION value=\"ba\"$selBA>AA: All peptides in Analysis</OPTION><OPTION value=\"banr\"$selBD>DA: Distinct pept. in Analysis</OPTION>";
		$analysisStrg='From current Analysis';
	}
	$pepSelectStrg.='</SELECT>';
# 	my $scoreValue='Score';
	my $specificityValue='Spec.<FONT style="font-size:9px;"> %</FONT>';
	my $coverageValue='Cov.<FONT style="font-size:9px;"> %</FONT>';
	my $organismValue='Species';
	if ($view eq 'protein'){
		$proteinValue="<FONT color=#DD0000>$proteinValue</FONT>";
	}
	elsif ($view eq 'mass') {
		$massValue="<FONT color=#DD0000>$massValue</FONT>";
	}
# 	elsif ($view eq 'score'){
# 		$scoreValue='<FONT color=#DD0000>Score</FONT>';
# 	}
	elsif ($view eq 'peptide') {
		$peptideValue="<FONT color=#DD0000>$peptideValue</FONT>";
	}
	elsif ($view eq 'specificity') {
		$specificityValue="<FONT color=#DD0000>$specificityValue</FONT>";
	}
	elsif ($view eq 'coverage') {
		$coverageValue="<FONT color=#DD0000>$coverageValue</FONT>";
	}
	elsif ($view eq 'organism') {
		$organismValue="<FONT color=#DD0000>$organismValue</FONT>";
	}
	#print "<TH width=250 align=left>&nbsp;$checkAllStrg&nbsp;<A href=\"javascript:selectView('protein')\" onmouseover=\"popup('Click to sort proteins by <B>ascending name</B>.')\" onmouseout=\"popout()\">$proteinValue</A>";
	print "<TH width=250 align=left class=\"rbBorder\">&nbsp;<A href=\"javascript:selectView('protein')\" onmouseover=\"popup('Click to sort proteins by <B>ascending name</B>.')\" onmouseout=\"popout()\">$proteinValue</A>";
	print "+<FONT class=\"help\" onmouseover=\"popup('<B>Post-Trans. Modifications:</B>$itemPTMtext')\" onmouseout=\"popout()\">PTMs</FONT>" if $itemPTMtext;
	print "</TH>\n";
	print qq
|	<TH width=$geneNameWidth class="rbBorder">Gene name</TH>
	<TH width=70 class="rbBorder"><A href="javascript:selectView('mass')" onmouseover="popup('Click to sort proteins by <B>decreasing mass</B>.')" onmouseout="popout()">$massValue</A></TH>
	<TH width=55 class="rbBorder"><A href="javascript:selectView('peptide')" onmouseover="popup('<B><U>Peptides:</U></B> All or distinct in item or best Analysis.<BR>Click to sort proteins by <B>decreasing number of peptides</B>.')" onmouseout="popout()">$peptideValue</A>$pepSelectStrg</TH>
	<TH width=60 class="rbBorder"><A href="javascript:selectView('specificity')" onmouseover="popup('<B>$analysisStrg.</B><BR>Click to sort proteins by <B>decreasing peptide specificity</B>.')" onmouseout="popout()">$specificityValue</A></TH>
	<TH width=60 class="rbBorder"><A href="javascript:selectView('coverage')" onmouseover="popup('<B>$analysisStrg.</B><BR>Click to sort proteins by <B>decreasing peptide coverage</B>.')" onmouseout="popout()">$coverageValue</A></TH>
	<TH width=570 class="bBorder">Description - <A href="javascript:selectView('organism')" onmouseover="popup('Click to sort proteins by <B>ascending organism</B>.')" onmouseout="popout()">$organismValue</A></TH>
</TR>
|;
##	print qq|</TABLE>
##<TABLE align=center border=0 cellspacing=0 cellpadding=0 width=$tableWidth>
##|;
#	<TH width=65><A href="javascript:selectView('score')" onmouseover="popup('Click to sort proteins by <B>decreasing score</B>.')" onmouseout="popout()">$scoreValue</A></TH>

	if ($maxListSize == 0) {
		print "<TR bgcolor=\"$color1\">";
		print "<TH align=left colspan=7>&nbsp;No protein</TH>\n";
		print "</TR>\n";
	}
}

######################################
####<<<Protein sort subroutines>>>####
######################################

####<Project hierarchy & raw list>####
sub itemSort {
	my ($sortItem)=@_;
	# Alias asc. (alias > peptide > score)
	if ($sortItem eq 'protein') {return lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0]) || $listGroups{$b}{$b}[0]<=>$listGroups{$a}{$a}[0] || $listGroups{$b}{$b}[2]<=>$listGroups{$a}{$a}[2]}
	# Protein mass desc (mass > peptide > score > alias)
	elsif ($sortItem eq 'mass') {return $listProteins{$b}[4]<=>$listProteins{$a}[4] || $listGroups{$b}{$b}[0]<=>$listGroups{$a}{$a}[0] || $listGroups{$b}{$b}[2]<=>$listGroups{$a}{$a}[2] || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	# Protein score desc (score > peptide > alias)
	# elsif ($sortItem eq 'score') {return $listGroups{$b}{$b}[2]<=>$listGroups{$a}{$a}[2] || $listGroups{$b}{$b}[0]<=>$listGroups{$a}{$a}[0] || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	# Peptide desc. (peptide > score > alias)
	elsif ($sortItem eq 'peptide') {return $listGroups{$b}{$b}[0]<=>$listGroups{$a}{$a}[0] || $listGroups{$b}{$b}[2]<=>$listGroups{$a}{$a}[2]  || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}

	# Specificity desc. (specificity > peptide > score > alias)
	elsif ($sortItem eq 'specificity') {return $listGroups{$b}{$b}[6]<=>$listGroups{$a}{$a}[6] || $listGroups{$b}{$b}[0]<=>$listGroups{$a}{$a}[0] || $listGroups{$b}{$b}[2]<=>$listGroups{$a}{$a}[2] || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	# Coverage desc. (coverage > peptide > score > alias)
	elsif ($sortItem eq 'coverage') {return $listGroups{$b}{$b}[5]<=>$listGroups{$a}{$a}[5] || $listGroups{$b}{$b}[0]<=>$listGroups{$a}{$a}[0] || $listGroups{$b}{$b}[2]<=>$listGroups{$a}{$a}[2] || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}

	# Organism asc. (organism > peptide > score > alias)
	elsif ($sortItem eq 'organism') {return lc($listProteins{$a}[2]) cmp lc($listProteins{$b}[2]) || $listGroups{$b}{$b}[0]<=>$listGroups{$a}{$a}[0] || $listGroups{$b}{$b}[2]<=>$listGroups{$a}{$a}[2] || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	else {return 0}
}

####<User-defined classification>####
sub classSort {
	my ($sortItem)=@_;
	# Alias asc. (alias > peptide > score)
	if ($sortItem eq 'protein') {return lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0]) || $listProteins{$b}[5]<=>$listProteins{$a}[5] || $listProteins{$b}[7]<=>$listProteins{$a}[7]}
	# Protein mass desc (mass > peptide > score > alias)
	elsif ($sortItem eq 'mass') {return $listProteins{$b}[4]<=>$listProteins{$a}[4] || $listProteins{$b}[5]<=>$listProteins{$a}[5] || $listProteins{$b}[7]<=>$listProteins{$a}[7] || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	# Protein score desc (score > peptide > alias)
	# elsif ($sortItem eq 'score') {return $listProteins{$b}[7]<=>$listProteins{$a}[7] || $listProteins{$b}[5]<=>$listProteins{$a}[5] || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	# Peptide desc. (peptide > score > alias)
	elsif ($sortItem eq 'peptide') {return $listProteins{$b}[5]<=>$listProteins{$a}[5] || $listProteins{$b}[7]<=>$listProteins{$a}[7] || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}

	# Specificity desc. (specificity > peptide > score > alias)
	elsif ($sortItem eq 'specificity') {return $listProteins{$b}[10]<=>$listProteins{$a}[10] || $listProteins{$b}[5]<=>$listProteins{$a}[5] || $listProteins{$b}[7]<=>$listProteins{$a}[7] || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	# Coverage desc. (coverage > peptide > score > alias)
	elsif ($sortItem eq 'coverage') {return $listProteins{$b}[9]<=>$listProteins{$a}[9] || $listProteins{$b}[5]<=>$listProteins{$a}[5] || $listProteins{$b}[7]<=>$listProteins{$a}[7] || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}

	# Organism asc. (organism > peptide > score > alias)
	elsif ($sortItem eq 'organism') {return lc($listProteins{$a}[2]) cmp lc($listProteins{$b}[2]) || $listProteins{$b}[5]<=>$listProteins{$a}[5] || $listProteins{$b}[7]<=>$listProteins{$a}[7] || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	else {return 0}
}


########################################
####<Convert VarMods hash to string>####
########################################
sub convertVarMods {
	my ($refVarMods,$type)=@_;
	my $varModText='';
	foreach my $varMod (sort{$a cmp $b} keys %{$refVarMods}) {
		if ($type eq 'protein') {
			$varModText.='+' unless $varModText;
			$varModText.="<FONT class=\"$projectVarMods{$varMod}[2]\">$projectVarMods{$varMod}[1]</FONT>";
		}
		else {
			$varModText.="<BR><FONT class=\\'popPTM $projectVarMods{$varMod}[2]\\'>$projectVarMods{$varMod}[1]</FONT>: $projectVarMods{$varMod}[0]";
		}
	}
	#return '+<FONT class="Acetyl">A</FONT><FONT class="Carbamidomethyl">C</FONT><FONT class="Methyl">M</FONT><FONT class="Oxidation">O</FONT><FONT class="Phospho">P</FONT>';
	return $varModText;
}

##############################################
#####<Generate VarMod code from DB string>####
##############################################
#sub getVarModCodes {
#	my ($refProtMods,$varModStrg)=@_;
#	$varModStrg=~s/^ *\+ //;
#	foreach my $varModPos (split(/ \+ /,$varModStrg)) {
#		#my ($varMod)=($varModPos=~/^(.+) \(/);
#		##(my $varMod=$varModPos)=~s/\s*\(.+//;
#		my ($varModCode,$resStrg,$posStrg)=&promsMod::convertVarModString($varModPos,1);
#		next unless $varModCode;
#		$refProtMods->{$varModCode}=1 if $projectVarMods{$varModCode}; # must be in project list
#	}
#}

###################################################
####<Correct peptide count based on PTM filter>####
###################################################
#sub correctPeptideCount {
#	my ($refListGroup,$refPepCount)=@_;
#
#
#
#}


sub ajaxUpdateCustomLists { # $dbh is global here
	my $projectID=param('projectID');

	print header(-type=>'text/plain',-charset=>'utf-8');
	warningsToBrowser(1);

	my (%classificationList,%categoryList);
	&getCustomList($dbh,$projectID,\%classificationList,\%categoryList);
	$dbh->disconnect;

	foreach my $classID (sort{lc($classificationList{$a}[0]) cmp lc($classificationList{$b}[0])} keys %classificationList) {
		if ($categoryList{$classID}) { # at least 1 category
			print "T:=:$classificationList{$classID}[0]:=:$classID\n";
			foreach my $catID (sort{$categoryList{$classID}{$a}[1]<=>$categoryList{$classID}{$b}[1]} keys %{$categoryList{$classID}}) {
				print "L:=:$categoryList{$classID}{$catID}[0]:=:$catID:=:$categoryList{$classID}{$catID}[2]\n";
			}
		}
	}
	exit;
}

sub getCustomList {
	my ($dbh,$projectID,$refClass,$refCat)=@_;

	%{$refClass}=&promsMod::getListClass($dbh,$projectID);
	my $sthCat=$dbh->prepare("SELECT ID_CATEGORY,NAME,DISPLAY_POS FROM CATEGORY WHERE ID_CLASSIFICATION=?");
	##>Check for modifiable Lists<##
	my @sthNoMod=(
		$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,COMPARISON CO WHERE CO.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=?"), # COMPARISON
		$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,CAT_COMPARISON CC WHERE CC.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=?"), # CAT_COMPARISON
		$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,GO_ANALYSIS GA WHERE GA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=?"), # GO_ANALYSIS
		$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,EXPLORANALYSIS EA WHERE EA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=?") # EXPLORANALYSIS
	);
	foreach my $classID (keys %{$refClass}) {
		$sthCat->execute($classID);
		while (my ($catID,$catName,$dispPos)=$sthCat->fetchrow_array) {
			@{$refCat->{$classID}{$catID}}=($catName,$dispPos,0); # 0 modif flag
		}
		foreach my $sth (@sthNoMod) {
			$sth->execute($classID);
			while (my ($catID)=$sth->fetchrow_array) {
				$refCat->{$classID}{$catID}[2]=1; # update modif flag
			}
		}
	}
	$sthCat->finish;
	foreach my $sth (@sthNoMod) {$sth->finish;}
}

sub getProjectItemDBList {
	my ($dbh, $projectID, $refDBList) = @_;
	my @idAnalysis = ();
	getProjectItemAnaIDs(\@idAnalysis);
	
	if(@idAnalysis) {
		my $sthDB = $dbh->prepare("SELECT D.ID_DATABANK, D.NAME FROM ANALYSIS A INNER JOIN ANALYSIS_DATABANK AD ON AD.ID_ANALYSIS = A.ID_ANALYSIS INNER JOIN DATABANK D ON D.ID_DATABANK = AD.ID_DATABANK WHERE A.ID_ANALYSIS IN (".join(',', @idAnalysis).")");
		$sthDB->execute();
		while (my ($dbID, $dbName)=$sthDB->fetchrow_array) {
			$refDBList->{$dbID}=$dbName;
		}
	}
}

sub getProjectItemAnaIDs {
	my ($refIDAnalysis) = @_;
	my $queryIDAnalysis = "SELECT ID_ANALYSIS FROM ANALYSIS A";
	my $projectItem = uc($item);
	
	if($projectItem eq 'ANALYSIS') {
		push @$refIDAnalysis, $itemID;
	} else {
		$queryIDAnalysis .= " INNER JOIN SAMPLE S ON S.ID_SAMPLE = A.ID_SAMPLE ";
		if($projectItem ne 'SAMPLE') {
			$queryIDAnalysis .= " INNER JOIN EXPERIMENT E ON E.ID_EXPERIMENT = S.ID_EXPERIMENT ";
		}
		
		if($projectItem eq 'PROJECT') {
			$queryIDAnalysis .= " WHERE E.ID_PROJECT = $projectID ";
		} else {
			$queryIDAnalysis .= " WHERE ".substr($projectItem, 0, 1).".ID_$projectItem = $itemID";
		}
		
		my $sthIdAna = $dbh->prepare($queryIDAnalysis);
		$sthIdAna->execute();
		while (my ($idAna)=$sthIdAna->fetchrow_array) {
			push(@$refIDAnalysis, $idAna); 
		}
	}
}

####>Revision history<####
# 2.4.7 [ENHANCEMENT][UX] Simplifies DOM elements management for Match group display (PP 13/11/19)
# 2.4.6 [BUGFIX] Proteins information was not displayed properly on class listing (VS 11/10/19)
# 2.4.5 [ENHANCEMENT] Speed up proteins listing (VS 06/08/19)
# 2.4.4 [BUGFIX] Error thrown when listing proteins of an empty experiment (VS 24/06/19)
# 2.4.3 [FEATURE] Add Databank selection to List filtering (VS 19/06/19)
# 2.4.2 Handles project status=-1 [no auto-end validation] (PP 07/06/18)
# 2.4.1 Fix minor bug in gene display of hidden proteins (PP 08/02/16)
# 2.4.0 Custom lists management function moved to promsMod.pm + ajax (PP 26/04/14)
# 2.3.9 Check if Custom lists can be modified (PP 15/05/14)
# 2.3.8 -charset=>'utf-8' (PP 22/04/14)
# 2.3.7 "help" font class now in promsStyle.css (FY 08/10/13)
# 2.3.6 Optimisation in var mod queries (PP 27/09/13)
# 2.3.5 Extends 2.3.4 to GO Analysis & fix width for Gene name column (PP 19/08/13)
# 2.3.4 No match group editing if analysis has protein quantification (PP 24/07/13)
# 2.3.3 Minor modification of column alignment -> Gene name width was not big enough... Add $geneNameWidth and $tableWidth to fix this adaptive column size (GA 26/06/13)
# 2.3.2 Minor modification of $sthPep query for ambiguous ID_PEPTIDE selection and virtual proteins (GA 17/05/13)
# 2.3.1 Proteins in Category but not in Project item are also displayed (PP 15/05/13)
# 2.3.0 Minor modification in $sthPepMod query to gain efficiency (GA 15/05/13)
# 2.2.9 New PTMs handling: suppression of &getVarModCodes, add/update of SQL queries and &getVariableModifications is now in promsMod (GA 13/05/13)
# 2.2.8 Handles virtual proteins (PP 06/03/13)
# 2.2.7 Handles '-' flag for peptide coverage of problematic sequence (PP 25/01/13)
# 2.2.6 Handles new mapping procedure. DAVID mapping is obsolete (PP 21/11/12)
# 2.2.5 Project status management & user lists renaming (PP 18/09/12)
# 2.2.4 Handles better long protein names (PP 29/05/12)
# 2.2.3 More ghost peptide skipping (PP 29/11/11)
# 2.2.2 Correction in varMods Class & ghost peptide skipping (PP 13/07/11)
# 2.2.1 header(-'content-encoding'=>'no') (PP 01/07/11)
# 2.2.0 varMod parsing done by promsMod (PP 29/04/11)
# 2.1.9 Adds Classification/Category as filter (PP 28/02/11)
# 2.1.8 New protein list management options: raw, to child/item, classification (PP 21/02/2011)
# 2.1.7 Adds PTMs to list of filter options (PP 15/12/2010)
