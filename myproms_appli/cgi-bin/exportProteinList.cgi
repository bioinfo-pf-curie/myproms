#!/usr/local/bin/perl -w

################################################################################
# exportProteinList.cgi        2.6.0	                                         #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Exports a list of proteins in MS Excel or HTML format                        #
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
# data to be ensured and, more generally, to use and operate it ina the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#-------------------------------------------------------------------------------
$|=1;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use Encode 'decode_utf8';
#use utf8; # needed only if utf8 char from script itself must be printed to Excel
use POSIX qw(strftime);
use Archive::Zip;
use Spreadsheet::WriteExcel;
use promsConfig;
use promsMod;
use promsQuantif;
use phosphoRS;
use XML::Simple;
#use Storable 'dclone';
use strict;

#print header; warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my ($color1,$color2)=&promsConfig::getRowColors;
my $userID=$ENV{'REMOTE_USER'};

#################################
####>Connect to the database<####
#################################
my $dbh=&promsConfig::dbConnect;

if (param('ACT')) {
	if (param('ACT') eq 'export') {&exportAnalysisFiles; exit;}
	elsif (param('ACT') eq 'download') {&downloadFile; exit;}
}

####################
####>Parameters<####
####################
my $itemID=param('itemID');
my $item=lc(param('item'));
my $ITEM=uc($item);
if (param('AJAX')) {
	if (param('AJAX') eq 'ajaxCheckMapping') {&checkMapping;}
	exit;
}
my $projectID=param('projectID');
my $listMode=param('listMode'); # not defined after submit
my $classificationID=(!$listMode)? param('classID') : ($listMode=~/classification:(\d+)/)? $1 : ($listMode eq 'raw')? -1 : 0; # converts new list prot option to old classification
my $listID=(param('catList'))? param('catList') : ''; # selected list
my $restrictListID=param('restrictList') || 0; # exclusion/selection List ID
my $listAction=param('listAction') || 'restrict';
my $view=(param('view'))? param('view') : 'pepAll'; # depth=analysis by default

##>Global variables<##
my (%listProteins,%classProteins,%restrictProteins); # global because of listItemProteins subroutine
my (%bestScore,%numPeptides,%pepCoverage,%timesFound,%proteinGeneSymbol,%protSpecificity,%protPSM); # global because of sort subroutine
#my %fontSize=('PROJECT'=>'24px','EXPERIMENT'=>'20px','GEL2D'=>'17px','SPOT'=>'15px','SAMPLE'=>'15px','ANALYSIS'=>'13px','TITLE'=>'28px','CATEGORY'=>'20px','COMMENTS'=>'11px');
my (%fontSize,@selColName,@selPepColName,@protColumns,@numPepColumns,@checkPosDef,@mappingColumns,%mappingColURLs,%protIdentifierLinks,%mappingLinks);
my ($exportFormat,$colSpan,$rowSpan,$colSpanNumPep,$colSpanPep,$badConfidence,$ghostProteins,$geneSymbolPos,$geneNamePos,$mwColPos,$newLine,$restrictThemeID);
my ($workbook,$worksheet,$worksheet2,%format,%itemFormat,@columnSizes,$xlsRow,$rowColor1,$rowColor2);
my %relevantPTMs;

##################
####>>>Main<<<####
##################
my ($itemName,$itemDes)=$dbh->selectrow_array("SELECT NAME,DES FROM $ITEM WHERE ID_$ITEM=$itemID");
($itemDes)=&promsMod::chkDef($itemDes);

####>Fetching list of classifications<####
#my %classificationList=&promsMod::getListClass($dbh,$projectID);
#$classificationList{0}='Project hierarchy';
my $sthL=$dbh->prepare("SELECT T.ID_CLASSIFICATION,T.NAME,L.ID_CATEGORY,L.NAME,L.DISPLAY_POS,LIST_TYPE FROM CLASSIFICATION T,CATEGORY L WHERE T.ID_CLASSIFICATION=L.ID_CLASSIFICATION AND T.ID_PROJECT=$projectID");
$sthL->execute;
my (%customList,%classificationList);
while (my ($themeID,$themeName,$listID,$listName,$listPos,$type)=$sthL->fetchrow_array) {
	$type='PROT' unless $type;
	$classificationList{$themeID}=$themeName;
	@{$customList{$themeID}{$listID}}=($listName,$listPos,$type);
	$restrictThemeID=$themeID if ($restrictListID && $restrictListID==$listID);
}
$sthL->finish;

####>Fetching list of identifier mappings<####
my %identifierCodes;
my $sthI=$dbh->prepare('SELECT CODE,NAME,RESRC_NAME,RESRC_URL,ID_IDENTIFIER FROM IDENTIFIER');
$sthI->execute;
while (my ($code,$identName,$resrcName,$resrcURL,$identID)=$sthI->fetchrow_array) {
	$resrcName=$identName if !$resrcName;
	@{$identifierCodes{$code}}=($identName,$resrcName,$resrcURL,$identID);
}
$sthI->finish;

####>Processing Form<####
&exportProteinList if param('export');

$dbh->disconnect;

############################
####>Starting HTML Form<####
############################
print header(-charset=>'utf-8');
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Option Menu</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.{font-weight:bold;}
.notMapped{color:#555555; text-decoration:line-through;}
#varMod {
	display: inline-block;
    position: relative;
    right: 3px;
    margin-top: 4px;
    margin-bottom: 2px;
}
</STYLE>
<SCRIPT type="text/javascript">
|;
&promsMod::popupInfo();
print qq
|// AJAX --->
function ajaxCheckMapping() {
	var mapDiv=document.getElementById('mappingDIV');
	mapDiv.innerHTML="<FONT class=\\"font11\\">Checking...</FONT><BR><IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\" width=200>";
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/exportProteinList.cgi?AJAX=ajaxCheckMapping&itemID=$itemID&item=$item",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			var usedIdentCodes=XHR.responseText.split(',');
			var selMapping=document.exportForm.sel_mapping;
			var numMapped=0;
			for (let i=0; i<selMapping.length; i++) {
				var matched=false;
				for (let c=0; c<usedIdentCodes.length; c++) {
					if (usedIdentCodes[c]==selMapping[i].value) {
						matched=true;
						break;
					}
				}
				if (matched) {numMapped++;}
				else {document.getElementById('mapText_'+selMapping[i].value).className='notMapped';}
			}
			mapDiv.style.display='none';
			var alertStrg='Proteins were mapped to '+numMapped+'/'+selMapping.length+' identifier';
			if (numMapped > 1) {alertStrg+='s';}
			alert(alertStrg);
		}
	}
	XHR.send(null);
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
function chooseClass(classID) {
	var myForm=document.exportForm;
	if (classID > 0) {
		document.getElementById('itemSpan').style.display='none';
		document.getElementById('catSpan').style.display='';
		document.getElementById('subFocusTitle').innerHTML = 'Selected List(s):';
		document.getElementById('catList').disabled = false;
		
		var elThemeSelect = document.getElementById("classID");
		var themeListsStr = elThemeSelect.options[elThemeSelect.selectedIndex].dataset.lists;
		var themeLists = themeListsStr.split('\$\$');
		var elListSelectHTML = "<option value='0'>All</option>";
		themeLists.forEach(function(list) {
			var listInfos = list.split('::');
			elListSelectHTML += "<option value='" + listInfos[0] + "'>" + listInfos[1] + "</option>";
		});
		document.getElementById('catList').innerHTML= elListSelectHTML;
		
		myForm.depth.disabled=true;
		myForm.listAction.disabled=true;
		myForm.restrictList.disabled=true;
		myForm.unclass.disabled=false;
		if ('$item'=='analysis') {
			myForm.sel_matchGr.disabled=false;
		}
		else {
			myForm.sel_times.disabled=false;
			myForm.sel_matchGr.disabled=true;
			myForm.sel_cumPepAll.disabled=false;
			myForm.sel_cumPepNr.disabled=false;
			myForm.sel_cumCov.disabled=false;
			myForm.sel_bestSpec.disabled=false;
			myForm.sel_cumSeqPep.disabled=false;
		}
		top.promsFrame.selectedMode='classification:'+classID; // update promsFrame variable
	}
	else { // project hierarchy
		document.getElementById('catSpan').style.display='none';
		document.getElementById('itemSpan').style.display='';
		document.getElementById('subFocusTitle').innerHTML = 'Depth:';
		document.getElementById('catList').disabled=true;
		myForm.depth.disabled=false;
		myForm.listAction.disabled=false;
		myForm.restrictList.disabled=false;
		myForm.unclass.disabled=true;
		if (myForm.depth.value == 'analysis') {
			myForm.sel_times.disabled=true;
			myForm.sel_matchGr.disabled=false;
			myForm.sel_cumPepAll.disabled=true;
			myForm.sel_cumPepNr.disabled=true;
			myForm.sel_cumCov.disabled=true;
			myForm.sel_bestSpec.disabled=true;
			myForm.sel_cumSeqPep.disabled=true;
			if (myForm.view.value.match('cum') \|\| myForm.view.value=='times_found') {
				myForm.view.selectedIndex=0;
			}
		}
		else {
			myForm.sel_matchGr.disabled=true;
		}
		top.promsFrame.selectedMode=('$item'=='analysis')? 'raw' : 'child'; // no need to update promsFrame variable...?
	}
}
function chooseDepth(depth) {
	var myForm=document.exportForm;
	if (depth=='analysis') {
		myForm.sel_times.disabled=true;
		myForm.sel_matchGr.disabled=false;
		myForm.sel_cumPepAll.disabled=true;
		myForm.sel_cumPepNr.disabled=true;
		myForm.sel_cumCov.disabled=true;
		myForm.sel_bestSpec.disabled=true;
		myForm.sel_cumSeqPep.disabled=true;
		// just in case...
		if (myForm.view.value.match('cum') \|\| myForm.view.value=='times_found') {
			myForm.view.selectedIndex=0;
		}
		if (myForm.sel_cumSeqPep.checked==true) {
			myForm.sel_cumSeqPep.checked=false;
			myForm.sel_seqPep.disabled=false;
			myForm.sel_modPep.disabled=true;
			myForm.sel_modPepProt.disabled=true;
			myForm.sel_etPep.disabled=true;
			myForm.sel_posPep.disabled=true;
			myForm.sel_flkPep.disabled=true;
			myForm.sel_occPep.disabled=true;
			myForm.sel_isSpec.disabled=true;
			myForm.sel_scorePep.disabled=true;
			myForm.sel_qvalPep.disabled=true;
			myForm.sel_chargePep.disabled=true;
			myForm.sel_massPep.disabled=true;
			myForm.sel_titlePep.disabled=true;
		}
	}
	else {
		myForm.sel_times.disabled=false;
		myForm.sel_matchGr.disabled=true;
		myForm.sel_cumPepAll.disabled=false;
		myForm.sel_cumPepNr.disabled=false;
		myForm.sel_cumCov.disabled=false;
		myForm.sel_bestSpec.disabled=false;
		//if (myForm.sel_cumSeqPep.checked==true \|\| myForm.sel_seqPep.checked==false) {
			myForm.sel_cumSeqPep.disabled=false;
		//}
	}
}
function chooseView(view) {
	if (checkedView && checkedView != 'identifier') {document.getElementById(checkedView).checked=false;} // uncheck previous view (sortItem)
	if (!view) {return;}
	if (document.exportForm.depth.disabled==false && document.exportForm.depth.value=='analysis') {
		if (view.match('cum') \|\| view=='timesFound') {
			alert('Not applicable when listing proteins by Analyses.');
			document.exportForm.view.selectedIndex=0;
			return;
		}
	}
	if (document.getElementById(view)) {document.getElementById(view).checked=true;}
	checkedView=view;
}
function checkSelection(box) {
	if (box.value==checkedView \|\| box.id==checkedView) { // box.id for geneName
		box.checked=true;
		alert('Data used to sort proteins must be exported.');
	}
}
var peptideOptions=['sel_modPep','sel_modPepProt','sel_etPep','sel_posPep','sel_flkPep','sel_occPep','sel_isSpec','sel_scorePep','sel_qvalPep','sel_chargePep','sel_massPep','sel_titlePep','sel_commentPep','sel_phosphoRS','sel_pepPSM']
function checkPepSequence(box) {
	var disStatus;
	if (box.checked) {
		disStatus=false;
		if (box.name=='sel_seqPep') {
			document.exportForm.sel_cumSeqPep.checked=false;
		}
		else { //sel_cumSeqPep
			document.exportForm.sel_seqPep.checked=false;
		}
	}
	else { // unchecked
		disStatus=true;
		//else if (document.exportForm.depth.disabled==false && document.exportForm.depth.value != 'analysis') { //sel_seqPep
		//	document.exportForm.sel_cumSeqPep.disabled=false;
		//}
	}
	for (var i=0; i<peptideOptions.length; i++) {document.exportForm[peptideOptions[i]].disabled=disStatus;}
}

function checkForm(myForm) {
	if (!myForm.view.value) {
		alert('Select a sort option');
		return false;
	}
	return true;
}
var checkedView;
var reportWin;
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" onload="chooseView('$view')">
<CENTER>
<FONT class="title">Export Proteins Contained in <FONT color="#DD0000">$itemName</FONT></FONT>
<BR><BR>
<FORM name="exportForm" method="post" onsubmit="return checkForm(this)">
<INPUT type="hidden" name="itemID" value="$itemID">
<INPUT type="hidden" name="item" value="$item">
<INPUT type="hidden" name="projectID" value="$projectID">
<TABLE align=center border=0 bgcolor="$color2" cellpadding=0 cellspacing=10>
<TR><TD colspan=4><TABLE width=100% border=0>
<TR><TD align=right><FONT class="title2">Focus:</FONT></TD>
<TD colspan=3 nowrap><SELECT id='classID' name="classID" onchange="chooseClass(this.value)" style="font-size:18px;">
|;
# my $selString=(==-1)? 'selected' : '';
# print "<OPTION value=-1 $selString>Raw list</OPTION>\n";
my $selString=($classificationID<=0)? 'selected' : '';
print qq
|<OPTION value=0 $selString>Project hierarchy</OPTION>
<OPTGROUP label="Themes:">
|;
foreach my $themeID (sort{lc($classificationList{$a}) cmp lc($classificationList{$b})} keys %classificationList){
	next if $themeID<=0; # raw list, proj hierarchy
	$selString=($classificationID==$themeID)? 'selected' : '';
	my @themeLists = ();
	foreach my $listID (sort{lc($customList{$themeID}{$a}[1]) cmp lc($customList{$themeID}{$b}[1])} keys %{$customList{$themeID}}) {
		push(@themeLists, "$listID\::".$customList{$themeID}{$listID}[0]);
	}
	print "<OPTION data-lists='".join('$$', @themeLists)."' value=\"$themeID\" $selString>$classificationList{$themeID}</OPTION>\n";
}
print "</OPTGROUP></SELECT>\n";

##>Expand tree options
my ($visItemSpan,$visCatSpan,$disabledString1, $disabledString5)=($classificationID>0)? ('none','','disabled','') : ('','none','', 'none');
my $displayFieldName = ($classificationID>0) ? 'Selected List(s)' : 'Depth';
print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span id='subFocusTitle' class=\"title2\">$displayFieldName:</span>&nbsp;";
my $itemOK=0;
print "<SPAN id=\"itemSpan\" style=\"display:$visItemSpan\"><SELECT name=\"depth\" style=\"font-size:18px;\" onchange=\"chooseDepth(this.value)\" $disabledString1>\n";
foreach my $it ('project','experiment','gel2d sample','spot','analysis') {
	next if ($item eq 'sample' && $it eq 'spot');
	my $rawString='';
	if ($it=~/$item/) {$itemOK=1; $rawString=' (Raw List)';}
	next unless $itemOK; # skip parent
	#my ($extraIt,$extraTxt)=($it eq 'gel2d' || $it eq 'spot')? (' sample',' or free Sample') : ('','');
	my $extraIt=($it eq 'spot' && $item ne 'spot')? ' sample' : '';
	my $extraTxt=(($it=~/gel2d/ || $it eq 'spot') && $item ne 'gel2d' && $item ne 'spot' && $item ne 'sample')? ' or free Samples' : '';
	print "<OPTION value='$it$extraIt'";
	print ' selected' if $it eq 'analysis';
	my $refIt=($item eq 'sample' && $it=~/gel2d/)? 'sample' : ($it=~/gel2d/)? 'gel2d' : $it;
	my $itemStrg=($refIt eq 'analysis')? 'Analyses' : &promsMod::getItemType($refIt).'s';
	print ">$itemStrg$extraTxt$rawString</OPTION>\n";
}
print qq |
		</SELECT><SUP>*</SUP>
	</SPAN>
	<SPAN id=\"catSpan\" style=\"display:$visCatSpan\">
		<SELECT id='catList' name='catList' style='font-size:18px;' $disabledString5>
			<OPTION value='0'>All</OPTION>
|;

if($classificationID>0) {
	foreach my $listID (sort{lc($customList{$classificationID}{$a}[1]) cmp lc($customList{$classificationID}{$b}[1])} keys %{$customList{$classificationID}}) {
		print("<OPTION value='$listID'>".$customList{$classificationID}{$listID}."</OPTION>");
	}
}

print qq |
		</SELECT>
	</SPAN></TD></TR>\n
|;

##>Restrict/exclude List option
print qq
|<TR><TD align=right nowrap><SPAN id="listActSpan">
	<SELECT name="listAction" class="title3"><OPTION value="restrict">Restrict to</OPTION><OPTION value="exclude">Exclude</OPTION></SELECT><FONT class="title3">:</FONT></SPAN>
	</TD>
	<TD><SELECT name="restrictList" class="title3"><OPTION value="0">-= No List filtering =-</OPTION>
|;
foreach my $themeID (sort{lc($classificationList{$a}) cmp lc($classificationList{$b})} keys %classificationList){
	print "<OPTGROUP label=\"$classificationList{$themeID}\">\n";
	foreach my $listID (sort{lc($customList{$themeID}{$a}[1]) cmp lc($customList{$themeID}{$b}[1])} keys %{$customList{$themeID}}) {
		my $typeStrg=($customList{$themeID}{$listID}[2] eq 'SITE')? ' [Sites]' : '';
		print "<OPTION value=\"$listID\">$customList{$themeID}{$listID}[0]$typeStrg</OPTION>\n";
	}
	print "</OPTGROUP>\n";
}
print "</SELECT></TD>\n";

##>Sort options
my %viewOptions=(
	1=>['identifier','Identifier'],
	2=>['alias','Name in myProMS'],
	3=>['geneName','Gene Name'],
	4=>['mw','Molecular weight'],
	5=>['bestScore','Protein score in best Analysis'],
	6=>['bestCov','Coverage in best Analysis'],
	9=>['pepAll','All peptides in best Analysis'],
	10=>['pepNr','Distinct peptides in best Analysis']
);
if ($item ne 'analysis') {
	@{$viewOptions{7}}=('cumCov','Coverage (cumulated in displayed Item)');
	@{$viewOptions{8}}=('timesFound','Times found (visible)');
	@{$viewOptions{11}}=('cumPepAll','All peptides in displayed Item');
	@{$viewOptions{12}}=('cumPepNr','Distinct peptides in displayed Item');
}
print "<TD align=right nowrap><FONT class=\"title3\">&nbsp;Sort proteins by:</FONT></TD>\n";
print "<TD><SELECT name=\"view\" class=\"title3\" onchange=\"chooseView(this.value)\">\n";
print "<OPTION value=''>-=Select a sort option=-</OPTION>\n"; # used if currently selected option is no longer compatible with form
foreach my $sortRank (sort{$a<=>$b} keys %viewOptions){
	my $selString=($viewOptions{$sortRank}[0] eq $view)? 'selected' : '';
	print "<OPTGROUP label=\"Matching peptides:\">\n" if $sortRank==9;
	print "<OPTION value=\"$viewOptions{$sortRank}[0]\" $selString>$viewOptions{$sortRank}[1]</OPTION>\n";
}
my $disabledString2=($item eq 'analysis' || $classificationID==0)? 'disabled' : ''; # because default depth is 'analysis'
my $disabledString3=($classificationID>0)? '' : 'disabled';
my $disabledString4=($item ne 'analysis' && $classificationID>0)? 'disabled' : '';
print qq
|</OPTGROUP></SELECT></TD></TR></TABLE></TD></TR>
<TR><TH colspan=4><FONT style="font-size:18px;">Data to be exported:</FONT></TH></TR>
<TR><TH align=left valign=top rowspan=2 nowrap>
	&nbsp;&nbsp;&bull;&nbsp;<FONT style="font-size:15px;">Proteins:</FONT><BR>
	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_prot" id="identifier" value="IDENTIFIER" checked disabled><B>&nbsp;Identifier<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_prot" id="alias" value="ALIAS" onclick="checkSelection(this)">&nbsp;Name in myProMS<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_prot" value="PROT_DES">&nbsp;Description<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_prot" value="ORGANISM">&nbsp;Species<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_prot" value="PROT_LENGTH">&nbsp;Protein length<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_prot" id="mw" value="MW">&nbsp;Molecular weight<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_prot" value="PROT_SEQ">&nbsp;Sequence<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_protVmod">&nbsp;Relevant PTMs<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_prot" value="COMMENTS">&nbsp;Comments<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_matchGr" value="1" $disabledString4>&nbsp;Match group<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_times" id="timesFound" value="1" $disabledString2 onclick="checkSelection(this)">&nbsp;Times found (visible)<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_score" id="bestScore" value="1" onclick="checkSelection(this)">&nbsp;Score in best Analysis<BR>
	<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="unclass" value="1" $disabledString3>&nbsp;Also export proteins not in Theme
</TH>
<TH align=left valign=top rowspan=2 nowrap>
	&nbsp;&nbsp;&bull;&nbsp;<FONT style="font-size:15px;">Resource mapping:&nbsp;&nbsp;</FONT><BR>
	<DIV id="mappingDIV">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="button" class="font11" value="Check availability" onclick="ajaxCheckMapping()"><BR></DIV>
|;
foreach my $code (sort{lc($identifierCodes{$a}[0]) cmp lc($identifierCodes{$b}[0])} keys %identifierCodes) {
	print "&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type=\"checkbox\" name=\"sel_mapping\" value=\"$code\"";
	print " id=\"geneName\" onclick=\"checkSelection(this)\"" if $code eq 'GN'; # only for gene symbol (used in sort option)
	print ">&nbsp;<FONT id=\"mapText_$code\" onmouseover=\"popup('<B>$identifierCodes{$code}[1]</B>')\" onmouseout=\"popout()\">$identifierCodes{$code}[0]</FONT><BR>\n";
}
print qq
|</TH>
<TH align=left valign=top nowrap>
	&nbsp;&nbsp;&bull;&nbsp;<FONT style="font-size:15px;">Matching peptides:</FONT><BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-&nbsp;Peptide type<SELECT name="sel_pepType" id="pepType"><OPTION value="ALL">All</OPTION><OPTION value="TYPIC">Proteotypic</OPTION></SELECT><BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_pepPSM" value="1" disabled>&nbsp;PSM<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-&nbsp;Number in best Analysis<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_pepAll" id="pepAll" value="1" onclick="checkSelection(this)">&nbsp;All peptides<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_pepNr" id="pepNr" value="1" onclick="checkSelection(this)">&nbsp;Distinct peptides<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-&nbsp;Number in displayed Item<SUP>*</SUP><BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_cumPepAll" id="cumPepAll" value="1" onclick="checkSelection(this)" disabled>&nbsp;All peptides<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_cumPepNr" id="cumPepNr" value="1" onclick="checkSelection(this)" disabled>&nbsp;Distinct peptides<BR>

	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_cov" id="bestCov" value="1" onclick="checkSelection(this)">&nbsp;Coverage in best Analysis<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_cumCov" id="cumCov" value="1" $disabledString2 onclick="checkSelection(this)">&nbsp;Cumulated coverage in displayed Item<SUP>*</SUP><BR>
	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_spec" value="1">&nbsp;Specificity in best Analysis<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_bestSpec" value="1" disabled>&nbsp;Best specificity in displayed Item<SUP>*</SUP><BR>
</TH>
<TH align=left valign=top nowrap>
	&bull;&nbsp;<FONT style="font-size:15px;">Peptide data:</FONT><BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_noDupPep" value="1" checked>&nbsp;Distinct peptides only<BR>
	&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_seqPep" value="1" onclick="checkPepSequence(this)">&nbsp;Sequences in best analysis&nbsp&nbsp&nbsp<BR>
	&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_cumSeqPep" value="1" onclick="checkPepSequence(this)" $disabledString2>&nbsp;Sequences in displayed Item<SUP>*</SUP><BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_flkPep" value="1" disabled>&nbsp;Flanking residues<BR>
	<span id='varMod'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Variable modifications<BR/>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_modPep" style="margin-top: 5px" value="1" disabled>&nbsp;Position in peptide<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_modPepProt" value="1" disabled>&nbsp;Position in protein</span><br/>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_etPep" value="1" disabled>&nbsp;Retention time<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_phosphoRS" value="1" disabled>&nbsp;PhosphoRS results<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_posPep" value="1" disabled>&nbsp;Position(s)<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_occPep" value="1" disabled>&nbsp;Occurence in protein<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_isSpec" value="1" disabled>&nbsp;Specificity for protein<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_scorePep" value="1" disabled>&nbsp;Score<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_qvalPep" value="1" disabled>&nbsp;Q-value<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_chargePep" value="1" disabled>&nbsp;Charge<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_massPep" value="1" disabled>&nbsp;Mass (observed)<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_commentPep" value="1" disabled>&nbsp;Validation comments<BR>
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sel_titlePep" value="1" disabled>&nbsp;Spectra title<BR>
</TH></TR>
<TR><TH valign=bottom colspan=2>
	<FONT class="title3">Export format:</FONT><SELECT name="expFormat" class="title3"><OPTION value="HTML">HTML</OPTION><OPTION value="XLS">Excel</OPTION></SELECT>
	&nbsp;&nbsp;&nbsp;<INPUT type="submit" name="export" value="   Export List   " style="font-weight:bold;"/>
</TH></TR>
</TABLE>
<I>*Only visible proteins will be exported.</I>
</FORM>
</CENTER>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">
	setPopup()
	chooseClass(document.getElementById('classID').value);
</SCRIPT>
</BODY>
</HTML>
|;

########################################################################################################

###################################################
####<<<Checking for valid mapped identifiers>>>####
###################################################
sub checkMapping {

	####<Starting HTML>####
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);

	####<Connect to the database>####
	my $dbh=&promsConfig::dbConnect;

	####<Fetching list of mappings available for selected item>####
	my %identCodes;
	my $sthI=$dbh->prepare('SELECT ID_IDENTIFIER,CODE FROM IDENTIFIER');
	$sthI->execute;
	while (my ($identID,$code,)=$sthI->fetchrow_array) {
		$identCodes{$identID}=$code;
	}
	$sthI->finish;

	my @sthAI;
	if ($ITEM eq 'PROJECT') {
		$sthAI[0]=$dbh->prepare("SELECT DISTINCT ID_IDENTIFIER FROM MASTERPROT_IDENTIFIER MI
									INNER JOIN PROTEIN P ON P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN
									INNER JOIN ANALYSIS_PROTEIN AP ON AP.ID_PROTEIN=P.ID_PROTEIN
									WHERE P.ID_PROJECT=$itemID AND VISIBILITY>=1");
		#AND P.ID_PROJECT=$itemID AND E.ID_EXPERIMENT NOT IN (SELECT E2.ID_EXPERIMENT FROM EXPERIMENT E2 INNER JOIN USER_EXPERIMENT_LOCK EU ON EU.ID_EXPERIMENT=E2.ID_EXPERIMENT WHERE E2.ID_PROJECT=$itemID AND EU.ID_USER='$userID')");
	}
	elsif ($ITEM eq 'EXPERIMENT') {
		$sthAI[0]=$dbh->prepare("SELECT DISTINCT ID_IDENTIFIER FROM MASTERPROT_IDENTIFIER MI,PROTEIN P,ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S WHERE AP.ID_PROTEIN=P.ID_PROTEIN AND P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND VISIBILITY>=1 AND AP.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_SPOT IS NULL AND ID_EXPERIMENT=$itemID");
		$sthAI[1]=$dbh->prepare("SELECT DISTINCT ID_IDENTIFIER FROM MASTERPROT_IDENTIFIER MI,PROTEIN P,ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S,SPOT SP,GEL2D G WHERE AP.ID_PROTEIN=P.ID_PROTEIN AND P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND VISIBILITY>=1 AND AP.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_SPOT=SP.ID_SPOT AND SP.ID_GEL2D=G.ID_GEL2D AND G.ID_EXPERIMENT=$itemID");
	}
	elsif ($ITEM eq 'GEL2D') {
		$sthAI[0]=$dbh->prepare("SELECT DISTINCT ID_IDENTIFIER FROM MASTERPROT_IDENTIFIER MI,PROTEIN P,ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S,SPOT SP WHERE AP.ID_PROTEIN=P.ID_PROTEIN AND P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND VISIBILITY>=1 AND AP.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_SPOT=SP.ID_SPOT AND ID_GEL2D=$itemID");
	}
	elsif ($ITEM eq 'SPOT') {
		$sthAI[0]=$dbh->prepare("SELECT DISTINCT ID_IDENTIFIER FROM MASTERPROT_IDENTIFIER MI,PROTEIN P,ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S WHERE AP.ID_PROTEIN=P.ID_PROTEIN AND P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND VISIBILITY>=1 AND AP.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND ID_SPOT=$itemID");
	}
	elsif ($ITEM eq 'SAMPLE') {
		$sthAI[0]=$dbh->prepare("SELECT DISTINCT ID_IDENTIFIER FROM MASTERPROT_IDENTIFIER MI,PROTEIN P,ANALYSIS_PROTEIN AP,ANALYSIS A WHERE AP.ID_PROTEIN=P.ID_PROTEIN AND P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND VISIBILITY>=1 AND AP.ID_ANALYSIS=A.ID_ANALYSIS AND ID_SAMPLE=$itemID ORDER BY NUM_PEP DESC,VISIBILITY DESC,CONF_LEVEL DESC");
	}
	else { # ANALYSIS
		$sthAI[0]=$dbh->prepare("SELECT DISTINCT ID_IDENTIFIER FROM MASTERPROT_IDENTIFIER MI,PROTEIN P,ANALYSIS_PROTEIN AP WHERE AP.ID_PROTEIN=P.ID_PROTEIN AND P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND VISIBILITY>=1 AND ID_ANALYSIS=$itemID");
	}
	my %identAvalaible;
	foreach my $sth (@sthAI) {
		$sth->execute;
		while (my ($identID)=$sth->fetchrow_array) {
			$identAvalaible{$identCodes{$identID}}=1;
		}
		$sth->finish;
	}

	$dbh->disconnect;

	print join(',',sort keys %identAvalaible);

	exit;
}

######################################################
####<<<Processing Form: Exporting List to Excel>>>####
######################################################
sub exportProteinList {
	
	#######################
	####<Starting HTML>####
	#######################
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Exporting data files</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Exporting protein data from <FONT color="#DD0000">$itemName</FONT></FONT>
<BR><BR>
<DIV id="waitDIV">
<BR><BR><BR><BR><BR><FONT class="title3">Fetching data...</FONT><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">
<BR><BR><FONT class="title3">Status:&nbsp;<SPAN id="waitSPAN"><SPAN>...</FONT>
</DIV>
|;
	
	###>Config<####
	mkdir "$promsPath{tmp}/scratch" unless -d "$promsPath{tmp}/scratch";
	my $exportDir="$promsPath{tmp}/scratch/export";
	&promsMod::cleanDirectory($exportDir,'15m');
	my $jobID=strftime("%Y%m%d%H%M%S",localtime).'_'.$userID;
	
	if (param('sel_protVmod')) {
		my %allPostTransModifs=&promsMod::getVariableModifications($dbh);
		my $sthGetPM=$dbh->prepare("SELECT ID_MODIFICATION FROM PROJECT_MODIFICATION WHERE ID_PROJECT=$projectID");
		$sthGetPM->execute;
		while (my ($modID)=$sthGetPM->fetchrow_array) {
			$relevantPTMs{$allPostTransModifs{$modID}[0]} = 1;
		}
	}

	@mappingColumns=param('sel_mapping');
	my ($addColMatchPeptides)=(param('sel_pepType') eq 'TYPIC')? ' Proteotypic ' : ' ';


	##################################
	####<Writing to EXCEL or HTML>####
	##################################
	$exportFormat=param('expFormat');
	my $exportFile;
	$xlsRow=0;
	($rowColor1,$rowColor2)=('#FFFFFF','#DDDDDD');
	my $focusStrg= ($classificationID && $listID) ? $customList{$classificationID}{$listID}[0] : ($classificationID) ? $classificationList{$classificationID} : 'Project hierarchy';
	if ($exportFormat eq 'XLS') { # EXCEL
		$exportFile="export_protein_list_$jobID.xls";
		$newLine="\n";
		%fontSize=('PROJECT'=>24,'EXPERIMENT'=>20,'GEL2D'=>17,'SPOT'=>15,'SAMPLE'=>15,'ANALYSIS'=>13,'TITLE'=>28,'CATEGORY'=>20,'COMMENTS'=>11);
		$workbook=Spreadsheet::WriteExcel->new("$exportDir/$exportFile");
		eval { # Perl 5.8 compatibility
			$workbook->set_properties(title=>'MS-based protein identification data',
									  author=>'myProMS server',
									  comments=>'Automatically generated with Perl and Spreadsheet::WriteExcel'
									  );
		};
		#$workbook->set_custom_color(40,224,224,255); # light light color #E0E0FF
		#$workbook->set_custom_color(41,176,176,255); # dark blue color  #B0B0FF
		$workbook->set_custom_color(40,189,215,255); # light light color  #BDD7FF (new myProMS colors V3.5+)
		$workbook->set_custom_color(41,128,179,255); # dark blue color  #80B3FF
		$worksheet=$workbook->add_worksheet('Data');
		$worksheet2=$workbook->add_worksheet('Links') if scalar @mappingColumns;
		$itemFormat{'TITLE'}=$workbook->add_format(align=>'center',size=>$fontSize{'TITLE'},bold=>1,border=>1);
		foreach my $item (keys %fontSize) {
			next if $item eq 'TITLE';
			$itemFormat{$item}=$workbook->add_format(align=>'left',size=>$fontSize{$item},bold=>1,bg_color=>41,border=>1);
		}
		$itemFormat{'red'}=$workbook->add_format(bg_color=>10);
		$itemFormat{'yellow'}=$workbook->add_format(bg_color=>13);
		$itemFormat{'green'}=$workbook->add_format(bg_color=>17);
		$itemFormat{'header'}=$workbook->add_format(align=>'center',size=>10,text_wrap=>1,bold=>1,bg_color=>40,border=>1);
		$itemFormat{'mergeHeader'}=$workbook->add_format(align=>'center',valign=>'vcenter',size=>10,text_wrap=>1,bold=>1,bg_color=>40,border=>1);
		$itemFormat{'mergeURL'}=$workbook->add_format(align=>'left',size=>10,bold=>1,italic=>1);
		$format{'ident'}=$workbook->add_format(align=>'left',valign=>'top',size=>10,border=>1);
		$format{'mergeIdent'}=$workbook->add_format(align=>'left',valign=>'top',size=>10,border=>1);
		$format{'identBC'}=$workbook->add_format(align=>'left',size=>10,color=>'grey',border=>1);
		$format{'mergeIdentBC'}=$workbook->add_format(align=>'left',valign=>'top',size=>10,color=>'grey',border=>1);
		$format{'identVP'}=$workbook->add_format(align=>'left',size=>10,italic=>1,border=>1);
		$format{'mergeIdentVP'}=$workbook->add_format(align=>'left',valign=>'top',size=>10,italic=>1,border=>1);
		$format{'text'}=$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>1,border=>1);
		$format{'mergeText'}=$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>1,border=>1);
		$format{'mergeColText'}=$workbook->add_format(align=>'left',size=>10,text_wrap=>1,border=>1);
		$format{'number'}=$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1);
		$format{'mergeNumber'}=$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1);
		$format{'number1d'}=$workbook->add_format(size=>10,align=>'center',valign=>'top',num_format=>'0.0',border=>1);
		$format{'mergeNumber1d'}=$workbook->add_format(size=>10,align=>'center',valign=>'top',num_format=>'0.0',border=>1);
		$format{'seq'}=$workbook->add_format(align=>'left',size=>10,font=>'Monospace',valign=>'top',border=>1);

		#print header(-type=>"application/vnd.ms-excel",-attachment=>"$itemName ($classificationList{$classificationID}).xls");
	}
	else { # HTML
		$newLine='<BR>';
		%fontSize=('PROJECT'=>'24px','EXPERIMENT'=>'20px','GEL2D'=>'17px','SPOT'=>'15px','SAMPLE'=>'15px','ANALYSIS'=>'13px','TITLE'=>'28px','CATEGORY'=>'20px','COMMENTS'=>'11px');
		#print header(-attachment=>"$itemName ($classificationList{$classificationID}).html");
		#warningsToBrowser(1);
		$exportFile="export_protein_list_$jobID.html";
		open(HTML,">$exportDir/$exportFile");
		print HTML qq
|<HTML>
<HEAD>
<TITLE>$itemName ($focusStrg)</TITLE>
<STYLE type="text/css">
.lowConf{color:DimGray;}
.noConf{font-style:italic;}
</STYLE>
</HEAD>
<BODY>
<TABLE border=1 cellspacing=1>
|;
	}

	my $depthITEM= (param('depth')) ? uc(param('depth')) : '';
	@protColumns=param('sel_prot');
	unshift @protColumns,'IDENTIFIER';
	push @selColName,('#','Identifier'); # add numb prot & identifier (mandatory)
	my $colPos=0;
	foreach my $col (@protColumns) {
		if ($col eq 'ALIAS') {push @selColName,'Name in myProMS'; push @checkPosDef,$colPos;}
		elsif ($col eq 'PROT_DES') {push @selColName,'Description'; push @checkPosDef,$colPos;}
		elsif ($col eq 'ORGANISM') {push @selColName,'Species'; push @checkPosDef,$colPos;}
		elsif ($col eq 'PROT_LENGTH') {push @selColName,'Size (aa)';}
		elsif ($col eq 'MW') {push @selColName,'Mol. weight (Da)'; $mwColPos=$colPos;}
		elsif ($col eq 'PROT_SEQ') {push @selColName,'Sequence'; push @checkPosDef,$colPos;}
		elsif ($col eq 'COMMENTS') {push @selColName,'Comments'; push @checkPosDef,$colPos;}
		$colPos++;
	}
	push @selColName,'Relevant PTMs' if param('sel_protVmod');
	#$mappingColURLs{0}=1; # identifier
	foreach my $i (0..$#mappingColumns) {
		push @selColName,$identifierCodes{$mappingColumns[$i]}[0];
		$mappingColURLs{$i}=1 if $identifierCodes{$mappingColumns[$i]}[2]; # link to resource
	}

	push @selColName,'Match group' if param('sel_matchGr');
	push @selColName,"Times found$newLine(visible)" if param('sel_times');
	if (param('sel_score')) {
		my $colName='Score';
		$colName.="$newLine(best Analysis)" if ($ITEM ne 'ANALYSIS' && $depthITEM ne 'ANALYSIS');
		push @selColName,$colName;
	}
	foreach my $colName (@selColName) {($colName eq '#')? push @columnSizes,4 : ($colName eq 'Ensembl')? push @columnSizes,18 : push @columnSizes,15;}
	foreach my $pepParam ('sel_pep','sel_cumPep') {
		foreach my $countType ('PSM','All','Nr') {
			next unless param("$pepParam$countType");
			my $colName='';
			if ($countType eq 'Nr') {$colName.='Distinct';}
			else {$colName.=$countType;}
			if ($pepParam eq 'sel_pep') {
				if ($ITEM ne 'ANALYSIS' && $depthITEM ne 'ANALYSIS') {$colName.=' in best Analysis';}
			}
			elsif ($pepParam eq 'sel_cumPep') {$colName.=' in #ITEM#';}
			push @numPepColumns,$colName;
		}
	}
	$colSpanNumPep=scalar @numPepColumns;
	foreach my $colName (@numPepColumns) {($colName eq 'All' && $colSpanNumPep==1)? push @columnSizes,25 : push @columnSizes,length($colName)*1.2;}
	$rowSpan=1;
	if ($colSpanNumPep) {
		push @selColName,"Matching${addColMatchPeptides}peptides";
		$rowSpan=2;
	}
	if (param('sel_cov')) {
		my $colName='Coverage (%)';
		$colName.="$newLine(best Analysis)" if ($ITEM ne 'ANALYSIS' && $depthITEM ne 'ANALYSIS');
		push @selColName,$colName;
		push @columnSizes,12;
	}
	if (param('sel_cumCov')) {
		push @selColName,"Coverage (%)$newLine(cum. in #ITEM#)";
		push @columnSizes,12; # anticipates #ITEM# replacement
	}
	if (param('sel_spec')) {
		my $colName='Specificity (%)';
		$colName.="$newLine(best Analysis)" if $ITEM ne 'ANALYSIS';
		push @selColName,$colName;
		push @columnSizes,15; # anticipates #ITEM# replacement
	}
	if (param('sel_bestSpec')) {
		my $colName="Specificity (%)$newLine(best in #ITEM#)";
		push @selColName,$colName;
		push @columnSizes,15; # anticipates #ITEM# replacement
	}
	$colSpan=scalar @selColName;
	$colSpan+=($colSpanNumPep-1) if $colSpanNumPep;
	if (param('sel_seqPep') || param('sel_cumSeqPep')) {
		my $pepLabel='Peptide data (';
		$pepLabel.=(param('sel_noDupPep'))? 'distinct' : 'all';
		if (param('sel_seqPep')) {
			$pepLabel.=($ITEM ne 'ANALYSIS' && $classificationID<=0)? ' in best Analysis' : ($ITEM ne 'ANALYSIS')? ' in #ITEM#' : '';
		}
		else {$pepLabel.=' in #ITEM#';}
		#$pepLabel.=($ITEM ne 'ANALYSIS' && param('sel_seqPep'))? ' in best Analysis' : ($ITEM ne 'ANALYSIS')? ' in #ITEM#' : '';
		$pepLabel.=')';
		push @selColName,$pepLabel;
		push @selPepColName,('#','Sequence');
		if (param('sel_modPep') || param('sel_modPepProt') || param('sel_etPep') || param('sel_phosphoRS')|| param('sel_posPep')|| param('sel_occPep') ||  param('sel_isSpec') || param('sel_scorePep') || param('sel_qvalPep') || param('sel_chargePep') || param('sel_massPep') || param('sel_commentPep') || param('sel_titlePep')) {
			$rowSpan=2;
			push @selPepColName,'Variable modifications (in peptide)' if param('sel_modPep');
			push @selPepColName,'Variable modifications (in protein)' if param('sel_modPepProt');
			push @selPepColName,'Variable modifications summary' if param('sel_modPepProt') || param('sel_modPep');
			push @selPepColName,'Retention time' if param('sel_etPep');
			push @selPepColName,'PhosphoRS results' if param('sel_phosphoRS');
			push @selPepColName,'PhosphoRS Probabilities' if param('sel_phosphoRS');
			push @selPepColName,'Start..End' if param('sel_posPep');
			push @selPepColName,'Occurence' if param('sel_occPep');
			push @selPepColName,'Is specific' if param('sel_isSpec');
			push @selPepColName,'Score' if param('sel_scorePep');
			push @selPepColName,'Q-value' if param('sel_qvalPep');
			push @selPepColName,'Charge' if param('sel_chargePep');
			push @selPepColName,'Mass Obs.' if param('sel_massPep');
			push @selPepColName,'Validation comments' if param('sel_commentPep');
			push @selPepColName,'Title' if param('sel_titlePep');
		}
	}
	foreach my $colName (@selPepColName) {
		my $colSize = ($colName eq '#') ? 4 : ($colName eq 'Sequence') ? 35 : length($colName)*1.2;
		push @columnSizes,$colSize;
	}

	$colSpanPep=scalar @selPepColName;
	$colSpan+=$colSpanPep;

	if ($exportFormat eq 'XLS') {
		foreach my $i (0..$#columnSizes) {$worksheet->set_column($i,$i,$columnSizes[$i]);}
	}

	########################
	####<Classification>####
	########################

	####<Project hierarchy + Raw list>####
	if ($classificationID == 0) {
		my $listRestricInfo='';
		if ($restrictListID) {
			$listRestricInfo=($listAction eq 'restrict')? 'Restricted to' : 'Excluding';
			$listRestricInfo.=' proteins in List '.$classificationList{$restrictThemeID}.' > '.$customList{$restrictThemeID}{$restrictListID}[0];
			&promsQuantif::fetchCustomList($dbh,$restrictListID,\%restrictProteins,1);
		}
		if ($exportFormat eq 'XLS') {
			$listRestricInfo=' ('.$listRestricInfo.')' if $listRestricInfo;
			$worksheet->merge_range(0,0,$xlsRow,$colSpan-1,"$focusStrg$listRestricInfo",$itemFormat{'TITLE'});
		}
		else {
			$listRestricInfo="<FONT style=\"font-size:$fontSize{CATEGORY}\"> ($listRestricInfo)</FONT>" if $listRestricInfo;
			print HTML "<TR><TH colspan=$colSpan><FONT style=\"font-size:$fontSize{TITLE}\">$focusStrg</FONT>$listRestricInfo</TH></TR>\n";
		}
		&displayItemContents(uc($item),$itemID,$depthITEM,'');
	}
	####<User-defined classification>####
	else {
		my $focusType = ($listID) ? "List" : "Theme";
		#my ($className,$classDes)=$dbh->selectrow_array("SELECT NAME,DES FROM CLASSIFICATION WHERE ID_CLASSIFICATION=$classificationID");
		if ($exportFormat eq 'XLS') {
			$worksheet->merge_range(0,0,0,$colSpan-1,&promsMod::getItemType($item).": ".decode_utf8($itemName).", $focusType: ".decode_utf8($focusStrg),$itemFormat{'TITLE'});
			#$worksheet->write_comment(0,0,$classDes) if $classDes;
		}
		else {
			print HTML "<TR><TH colspan=$colSpan><FONT style=\"font-size:$fontSize{TITLE}\">",&promsMod::getItemType($item),": $itemName, $focusType: $focusStrg";
			#print HTML "<BR><FONT style=\"font-size:$fontSize{CATEGORY}\">($classDes)</FONT>" if $classDes;
			print HTML "</FONT></TH></TR>\n";
		}
		
		##>Proteins in categories
		if($listID) {
			my ($catName) = $dbh->selectrow_array("SELECT NAME FROM CATEGORY WHERE ID_CATEGORY=$listID");
			&printWaitingMsg("Processing Custom List $catName");
			&listItemProteins(uc($item),$listID,0);
		} else {
			my @listChildren=&promsMod::getChildrenList($dbh,$classificationID,'classification');
			foreach my $catID (sort{$listChildren[0]{'LIST'}{$a}[2]<=>$listChildren[0]{'LIST'}{$b}[2]} keys %{$listChildren[0]{'LIST'}}){
				my ($catDes)=&promsMod::chkDef($listChildren[0]{'LIST'}{$catID}[1]);
				$catDes=" ($catDes)" if $catDes;
				
				if ($exportFormat eq 'XLS') {
					$worksheet->merge_range(++$xlsRow,0,$xlsRow,$colSpan-1,decode_utf8($listChildren[0]{'LIST'}{$catID}[0].$catDes),$itemFormat{'CATEGORY'});
				}
				else {
					print HTML "<TR><TH align=left colspan=$colSpan bgcolor=\"$color2\"><FONT style=\"font-size:$fontSize{CATEGORY}\">$listChildren[0]{LIST}{$catID}[0]$catDes</FONT></TH></TR>\n";
				}
				&printWaitingMsg("Processing Custom List $listChildren[0]{LIST}{$catID}[0]");
				&listItemProteins(uc($item),$catID,0);
			}
		}
		##>Unclassified proteins
		if (param('unclass')) {
			if ($exportFormat eq 'XLS') {
				$worksheet->merge_range(++$xlsRow,0,$xlsRow,$colSpan-1,'Unclassified Proteins',$itemFormat{'CATEGORY'});
			}
			else {
				print HTML "<TR><TH align=left colspan=$colSpan bgcolor=\"$color2\"><FONT style=\"font-size:$fontSize{CATEGORY}\">Proteins not in List(s)</FONT></TH></TR>\n";
			}
			&printWaitingMsg("Processing Proteins not in List(s)");
			&listItemProteins(uc($item),0,1);
		}
	}

	$dbh->disconnect;

	if ($exportFormat eq 'XLS') {
		$format{'mergeIdentBC'}->set_bold();
		$format{'mergeIdentVP'}->set_bold();
		$worksheet->merge_range(++$xlsRow,0,$xlsRow,5,'Proteins in grey indicate a bad confidence.',$format{'mergeIdentBC'}) if $badConfidence;
		$worksheet->merge_range(++$xlsRow,0,$xlsRow,5,'Proteins in italics are virtual.',$format{'mergeIdentVP'}) if $ghostProteins;

		####################
		####<Link Sheet>####
		####################
		if (scalar @mappingColumns) {
			$worksheet2->set_column(0,scalar @mappingColumns,20); # $#... + identifier
			$worksheet2->write_string(0,0,'Identifier',$format{'header'});
			my $col2=0;
			foreach my $i (0..$#mappingColumns) {
				next unless $mappingColURLs{$i};
				$worksheet2->write_string(0,++$col2,$identifierCodes{$mappingColumns[$i]}[0],$format{'header'});
			}
			my $xlsRow2=1;
			foreach my $protIdent (sort{lc($a) cmp lc($b)} keys %mappingLinks) {
				my ($linkURL,$maxLinkRowSpan)=@{$protIdentifierLinks{$protIdent}};
				next if (!$linkURL && !$maxLinkRowSpan);
				$maxLinkRowSpan=1 unless $maxLinkRowSpan; # can be 0
				if ($linkURL) {$worksheet2->write_url($xlsRow2,0,$linkURL,$protIdent);} #,$format{'ident2'}
				else {$worksheet2->write_string($xlsRow2,0,$protIdent);} #,$format{'ident2'}
				$col2=0;
				foreach my $i (0..$#mappingColumns) {
					next unless $mappingColURLs{$i}; #
					$col2++;
					next unless $mappingLinks{$protIdent}{$mappingColumns[$i]};
					my $row=$xlsRow2-1;
					foreach my $refLink (@{$mappingLinks{$protIdent}{$mappingColumns[$i]}}) {
						$worksheet2->write_url(++$row,$col2,$refLink->[0],$refLink->[1]);
					}
				}
				$xlsRow2+=$maxLinkRowSpan;
			}

			#<Link to myProMS publication
			$worksheet->write_url_range(++$xlsRow,0,$xlsRow,8,'http://www.ncbi.nlm.nih.gov/pubmed/17610305','Data exported from myProMS Server (Poullet et al. Proteomics. 2007 (15):2553-6)',$itemFormat{'mergeURL'});
			$worksheet2->write_url_range($xlsRow2,0,$xlsRow2,8,'http://www.ncbi.nlm.nih.gov/pubmed/17610305','Data exported from myProMS Server (Poullet et al. Proteomics. 2007 (15):2553-6)',$itemFormat{'mergeURL'});
		}
		$workbook->close();
	}
	else {
		print HTML "<TR><TD colspan=$colSpan><FONT style=\"font-size:$fontSize{COMMENTS}\">Proteins in italics indicate a bad confidence.</FONT></TD></TR>\n" if $badConfidence;
		print HTML "<TR><TD colspan=$colSpan><FONT style=\"font-size:$fontSize{COMMENTS}\">Proteins in italics are virtual.</FONT></TD></TR>\n" if $ghostProteins;
		print HTML qq
|</TABLE>
<FONT style="font-size:$fontSize{COMMENTS}">Data exported from myProMS Server <I>(<A href="http://www.ncbi.nlm.nih.gov/pubmed/17610305" target="_blank">Poullet et al. Proteomics. 2007 (15):2553-6</A>)</I></FONT>
</BODY>
</HTML>
|;
		close HTML;
	}
	
	##<Link to download
		print qq
|<BR><BR><BR>
<INPUT type="button" class="title2" value="Download file" onclick="window.location='$promsPath{cgi}/exportProteinList.cgi?ACT=download&FILE=$exportFile'"/>
</CENTER>
<SCRIPT type="text/javascript">
document.getElementById('waitDIV').style.display='none';
top.promsFrame.selectedAction = 'summary';
window.location='$promsPath{cgi}/exportProteinList.cgi?ACT=download&FILE=$exportFile'; // auto-start download
</SCRIPT>
</BODY>
</HTML>
|;	
	exit;
}

################################################################################
####   		Triggers download of exported protein file (HTML/XLS)   	    ####
################################################################################
sub downloadFile {
	my $fileName=param('FILE');
	print header(-type=>"application/octet-stream", -attachment=>$fileName);
	
	if(-e "$promsPath{tmp}/scratch/export/$fileName") {
	open(FILE,"$promsPath{tmp}/scratch/export/$fileName");
	while (<FILE>) {print $_;}
	close FILE;
	}
	exit;
}

################################################################################
####   		Recover and export search files from selected analyses   		####
################################################################################
sub exportAnalysisFiles {
	my $projectID = param('id_project');
	my @anaIDs = param('anaList');
	my $item = (param('item')) ? lc param('item') : '';
	my $genDesign = param('gen_design');
	
	#######################
	####<Starting HTML>####
	#######################
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq|<HTML>
<HEAD>
<TITLE>Exporting data files</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Exporting search files from selected <FONT color="#DD0000">Analyses</FONT></FONT>
<BR><BR>

|;

	# Creating archive that contain all msf files
	if ($projectID && scalar @anaIDs) {
		print qq|<DIV id="waitDIV">
<BR><BR><BR><BR><BR><FONT class="title3">Fetching data...</FONT><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">
<BR><BR><FONT class="title3">Status:&nbsp;<SPAN id="waitSPAN"><SPAN>...</FONT>
</DIV>
|;	
		###>Config<####
		mkdir "$promsPath{tmp}/scratch" unless -d "$promsPath{tmp}/scratch";
		&promsMod::cleanDirectory("$promsPath{tmp}/scratch/export",'15m');
		my $jobID=strftime("%Y%m%d%H%M%S",localtime);
		my $exportDir="$promsPath{tmp}/scratch/export/$jobID\_$userID";
		mkdir "$exportDir";

		my (@msfFiles, $currentAnaFile, $currentSampName);
		my $designFile = "$exportDir/export_project$projectID\_analyses_from_$item.xls"; # $promsPath{'data'}."/
		my ($expName, $sampName, $projName) = $dbh->selectrow_array("SELECT E.NAME, S.NAME, P.NAME FROM SAMPLE S INNER JOIN ANALYSIS A ON A.ID_SAMPLE=S.ID_SAMPLE INNER JOIN EXPERIMENT E ON E.ID_EXPERIMENT=S.ID_EXPERIMENT JOIN PROJECT P ON P.ID_PROJECT=E.ID_PROJECT WHERE P.ID_PROJECT=$projectID AND A.ID_ANALYSIS=$anaIDs[0]");
		my ($workbook, $worksheet, $bodyFormat, $headerFormat);
		my $zip = Archive::Zip->new();

		if($genDesign) {
			# Initialize design file parameters
			$workbook=Spreadsheet::WriteExcel->new($designFile);
			eval { # Perl 5.8 compatibility
				$workbook->set_properties(title=>'Design of ',
											author=>'myProMS server',
											comments=>"Describes analyses carried out for project $projName"
											);
			};
			
			# Set design formats
			$headerFormat = $workbook->add_format(size=>10, bold=>1, underline=>1, font=>'Arial', align=>'center');
			$bodyFormat = $workbook->add_format(size=>11, align=>'center', font=>'Calibri');
			$worksheet = $workbook->add_worksheet('Design');
			
			# Write design header
			$worksheet->write(0, 0, "Experiment Name", $headerFormat);
			$worksheet->write(0, 1, "Sample Name", $headerFormat);
			$worksheet->write(0, 2, "Identification", $headerFormat);
			$worksheet->write(0, 3, "Raw file", $headerFormat);
		}
		
		# Parse relevant analyses
		my $totAna=scalar @anaIDs;
		my $iAna = 1;
		my %usedDataFile;
		my $sthAna=$dbh->prepare("SELECT E.ID_EXPERIMENT, A.ID_ANALYSIS, A.NAME, DATA_FILE, A.FILE_FORMAT, S.NAME, WIFF_FILE, VALID_STATUS FROM SAMPLE S INNER JOIN ANALYSIS A ON A.ID_SAMPLE=S.ID_SAMPLE INNER JOIN EXPERIMENT E ON E.ID_EXPERIMENT=S.ID_EXPERIMENT JOIN PROJECT P ON P.ID_PROJECT=E.ID_PROJECT WHERE P.ID_PROJECT=$projectID AND A.ID_ANALYSIS IN (".join(', ', @anaIDs).")");
		my $sthFirstQ=$dbh->prepare("SELECT MIN(ID_QUANTIFICATION) FROM ANA_QUANTIFICATION WHERE ID_ANALYSIS=?");
		$sthAna->execute;
		while (my ($expID,$anaID,$anaName,$dataFile,$fileFormat,$currentSampName,$wiffFile,$validStatus) = $sthAna->fetchrow_array) {
			&printWaitingMsg("Processing Analysis $iAna/$totAna ($anaName)");
			#$currentSampName =~ s/\//-/;
			my ($fileExt, $fullFileName);
			
			if ($fileFormat =~ /MAXQUANT/ && $dataFile=~/mqpar.*\.xml|parameters\.txt/) { # MaxQuant
				$fileExt='';
				$fullFileName=$wiffFile;
				$currentAnaFile='proteinGroups.txt';
				unless ($usedDataFile{$currentAnaFile}) {
					$sthFirstQ->execute($anaID);
					my ($firstQuantID)=$sthFirstQ->fetchrow_array;
					$zip->addFile("$promsPath{quantification}/project_$projectID/quanti_$firstQuantID/proteinGroups.txt","proteinGroups.txt");
					$zip->addFile("$promsPath{quantification}/project_$projectID/quanti_$firstQuantID/peptides.txt","peptides.txt");
				}
			} elsif($fileFormat =~ /SWATH|SPECTRONAUT|SKYLINE/) { # OpenSwath, Spectronaut or Skyline
				$fullFileName = $dataFile;
				$dataFile = "$promsPath{peptide}/proj_$projectID/exp_$expID/$fullFileName";
				if (-e $dataFile) { # Actual storage path
					$zip->addFile($dataFile, $fullFileName) unless $usedDataFile{$dataFile};
				} else { # Old storage path
					$dataFile = "$promsPath{peptide}/proj_$projectID/ana_$anaID/$fullFileName";
					if (-e $dataFile) {
						$zip->addFile($dataFile, $fullFileName) unless $usedDataFile{$dataFile};
					}
				}
			} else { # Mascot or PD
				$fileExt = ($fileFormat eq 'SEQUEST.PDM') ? 'msf' : 'dat';
				(my $fileName = $dataFile ) =~ s/(?:_[0-9]+)?\..*$//;
				$fullFileName = "$fileName.$fileExt";
				my ($dataFile, $dataFileValid);
				$dataFile = ($fileExt eq 'msf') ? "$promsPath{valid}/multi_ana/proj_$projectID/$fullFileName" : "$promsPath{valid}/ana_$anaID/$fullFileName";
				$dataFileValid = "$promsPath{peptide}/proj_$projectID/ana_$anaID/$fullFileName";

				if (-e $dataFile || -e $dataFileValid) {
					my $zipOutFolder = ($item eq 'experiment') ? "$currentSampName/" : "";
					$zip->addFile((-e $dataFile)? $dataFile : $dataFileValid, $zipOutFolder.$dataFile) unless $usedDataFile{$dataFile};
				}
			}

			if($genDesign) {
				$worksheet->write($iAna, 0, $expName, $bodyFormat);
				$worksheet->write($iAna, 1, $currentSampName, $bodyFormat);
				$worksheet->write($iAna, 2, $fullFileName, $bodyFormat);
				$worksheet->write($iAna, 3, $wiffFile, $bodyFormat);
			}
			$usedDataFile{$dataFile}=1;
			$iAna++;
		}
		
		if($genDesign) {
			$worksheet->set_column(0, 0, 26);
			$worksheet->set_column(1, 1, 26);
			$worksheet->set_column(2, 2, 13);
			$worksheet->set_column(3, 3, 13);
			
			$workbook->close();
		
			$zip->addFile($designFile,"design.xls");
		}
		
		my $archiveFile = $projName;
		$archiveFile .= ($item eq 'experiment') ? " - $expName.zip" : ($item eq 'sample') ? " - $sampName.zip" : " - Analyses.zip";
		$archiveFile =~ s/[\/\\]+/-/g;
		
		# Export built archive
		&printWaitingMsg("Building archive");
		$zip->writeToFileNamed("$promsPath{tmp}/scratch/export/$archiveFile"); # my $status = 

		unlink $designFile;
		rmdir $exportDir;
		
		##<Link to download
		print qq |
				<BR><BR><BR>
				<INPUT type="button" class="title2" value="Download archive" onclick="window.location='$promsPath{cgi}/exportProteinList.cgi?ACT=download&FILE=$archiveFile'"/></CENTER>
				<SCRIPT type="text/javascript">
					document.getElementById('waitDIV').style.display='none';
					top.promsFrame.selectedAction = 'summary';
					window.location='$promsPath{cgi}/exportProteinList.cgi?ACT=download&FILE=$archiveFile'; // auto-start download
				</SCRIPT>
			</BODY>
		</HTML>|;	
	}
	else {
		print qq
|<BR><BR><FONT class="title2">Error: No analysis selected.</FONT>
</BODY>
</HTML>
|;
	}
	exit;
}

################################################################################
####   Generate modifications summary text based on a given varMod string   ####
################################################################################
sub generateVarModSummary {
	my ($varMod, $pepSeq, $shift) = @_;
	my %varMods; # All parsed modifications
	my $comments = "";
	my @allMods = split(/\s\+\s/, $varMod);
	
	if(@allMods) {
		foreach my $mod (@allMods) {
			# Phospho (ST:995.?.?[974,979,989]) => $modType ($modAA:$modPos[$modAmbigousPos])
			my ($modType, $modAA, $modPos, $modAmbigousPos) = ( $mod =~ /(.+)\s\(([A-Za-z-\s]+):?([0-9.?]+)?\[?([0-9,]*)\]?\)$/ );
			my $posInPeptide;
			
			next unless($modType);
			
			# Has positions
			if($modPos) {
				my @varModPos = split('\.', $modPos);
				foreach my $pos (@varModPos) {
					# Add nb of sites to total for the modType
					$varMods{$modType}{"amount"}++;
					
					next if ($pos eq '?');
					
					my $posInPeptide = ($shift) ? $pos-$shift : $pos;
					$modAA = substr($pepSeq, $posInPeptide-1, 1);
					push(@{$varMods{$modType}{"regular"}}, "$modAA$pos");
				}
				
				# Has ambigous positions
				if($modAmbigousPos) {
					@varModPos = split(',', $modAmbigousPos);
					foreach my $pos (@varModPos) {
						my $posInPeptide = ($shift) ? $pos-$shift : $pos;
						$modAA = substr($pepSeq, $posInPeptide-1, 1);
						push(@{$varMods{$modType}{"ambigous"}}, "$modAA$pos");
					}
				}
			# Modification is located on Protein N-term / N-term / C-term / Protein C-term
			} else {
				push(@{$varMods{$modType}{"regular"}}, "$modAA");
				$varMods{$modType}{"amount"}++;
			}
		}
		
		foreach my $modType (keys %varMods) {
			my $nbModType = $varMods{$modType}{"amount"};
			my $nbModAmbigous = ($varMods{$modType}{"regular"}) ? $nbModType - scalar @{$varMods{$modType}{"regular"}} : 0;
			my $ambigousTxt = ($nbModAmbigous > 1) ? "and/or" : "or";
			
			if ($comments) {
				$comments .= ($exportFormat eq 'XLS') ? "\n" : "<br/>";
			}
			$comments .= "$nbModType $modType";
			$comments .= " : ".join(', ', @{$varMods{$modType}{"regular"}}) if($varMods{$modType}{"regular"});
			$comments .= " and $nbModAmbigous" if($varMods{$modType}{"regular"} && $varMods{$modType}{"ambigous"});
			$comments .= " ambigous : ".join(" $ambigousTxt ", @{$varMods{$modType}{"ambigous"}}) if ($varMods{$modType}{"ambigous"});
		}
	}
	
	return ($comments) ? $comments : '-';
}

###############################################################
####<<<Fetching & displaying item's contents (recursive)>>>####
###############################################################
sub displayItemContents { # item is upper cased
	my ($curITEM,$curItemID,$depthITEM,$shiftSize)=@_;
	my ($itemColSpan,$validStrg)=($curITEM eq 'ANALYSIS')? ($colSpan-1,',VALID_STATUS') : ($colSpan,'');
	my ($name,$des,$comments,$validStatus)=$dbh->selectrow_array("SELECT NAME,DES,COMMENTS$validStrg FROM $curITEM WHERE ID_$curITEM=$curItemID");
	$comments=~s/<Data[^>]+>// if $comments; # remove auto-comment on analysis annotation import

	if ($exportFormat eq 'XLS') {
		$des=($des)? " ($des)" : '';
		$xlsRow++;
		my $xlsCol=0;
		if ($curITEM eq 'ANALYSIS') {
			my $validColor=($validStatus<=0)? 'red' : ($validStatus==1)? 'yellow' : 'green';
			$worksheet->write_blank($xlsRow,$xlsCol,$itemFormat{$validColor});
			$xlsCol++;
		}
		$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow,$xlsCol+$itemColSpan-1,$shiftSize.&promsMod::getItemType($curITEM).": ".decode_utf8($name.$des),$itemFormat{$curITEM});
		$worksheet->write_comment($xlsRow,$xlsCol,decode_utf8($comments)) if $comments;
	}
	else { # HTML
		$des=&promsMod::HTMLcompatible($des);
		$comments=&promsMod::HTMLcompatible($comments);
		$des=" ($des)" if $des;
		$comments="<BR><FONT style=\"font-size:$fontSize{ANALYSIS};font-style:italic\">Comments: </FONT><FONT style=\"font-size:$fontSize{ANALYSIS};font-weight:normal;font-style:italic\">$comments</FONT>" if $comments;
		print HTML "<TR>";
		if ($curITEM eq 'ANALYSIS') {
			my $validColor=($validStatus<=0)? 'red' : ($validStatus==1)? 'yellow' : 'green';
			print HTML "<TH bgcolor=\"$validColor\">&nbsp;</TH>";
		}
		print HTML "<TH align=left colspan=$itemColSpan bgcolor=\"$color2\">$shiftSize<FONT style=\"font-size:$fontSize{$curITEM}\">",&promsMod::getItemType($curITEM),": $name$des</FONT>$comments</TH></TR>\n";
	}
	if ($depthITEM =~ /$curITEM/) {
		&printWaitingMsg("Processing ".&promsMod::getItemType($curITEM)." $name");
		&listItemProteins($curITEM,$curItemID,0);
	}
	else {
		my @listChildren=&promsMod::getChildrenList($dbh,$curItemID,$curITEM);
		if (scalar keys %{$listChildren[0]{'LIST'}} || ($listChildren[1] && scalar keys %{$listChildren[1]{'LIST'}})) {
			foreach my $i (0..$#listChildren) {
				my $childITEM=$listChildren[$i]{'ITEM'};
				foreach my $childID (sort{$listChildren[$i]{'LIST'}{$a}[2]<=>$listChildren[$i]{'LIST'}{$b}[2]} keys %{$listChildren[$i]{'LIST'}}) {
					my $newShift=($exportFormat eq 'XLS')?  $shiftSize.'     ' : $shiftSize.'&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;';
					&displayItemContents($childITEM,$childID,$depthITEM,$newShift);
				}
			}
		}
		else {
			if ($exportFormat eq 'XLS') {$worksheet->merge_range(++$xlsRow,0,$xlsRow,$colSpan-1,'No proteins',$format{'mergeText'});}
			else {print HTML "<TR><TD colspan=$colSpan>&nbsp;No proteins</TD></TR>\n";}
			#print HTML "<TR><TD colspan=$colSpan>&nbsp</TD></TR>\n";
		}
	}
}


########################################################
####<<<Fetching & listing item's protein contents>>>####
########################################################
sub listItemProteins {
	my ($subjectITEM,$subjectID,$filter)=@_;
	%listProteins=%bestScore=%protSpecificity=%numPeptides=%pepCoverage=%protPSM=(); # reset to () %listClusters not reset
	my (%bestAnalysis,%protLength,%confLevel,%proteinMGrHead); #%listAnalyses,
	my %analysisIdentType;

	#################################
	####<Retrieving data from DB>####
	#################################

	####<List of queries>####
	my $fieldStrg='ID_ANALYSIS,SCORE,NUM_PEP,CONF_LEVEL,MATCH_GROUP,VISIBILITY,PEP_COVERAGE,PEP_SPECIFICITY,DB_RANK';
	my $visibleQuery;
	if ($subjectITEM eq 'ANALYSIS') {
		if ($filter) { # unclassified proteins
	 		$visibleQuery="SELECT ID_PROTEIN,$fieldStrg FROM ANALYSIS_PROTEIN WHERE VISIBILITY>0 AND ID_ANALYSIS=$itemID";
		}
		elsif ($classificationID > 0) { # user-defined
			$visibleQuery="SELECT ANALYSIS_PROTEIN.ID_PROTEIN,$fieldStrg FROM ANALYSIS_PROTEIN,CATEGORY_PROTEIN WHERE ANALYSIS_PROTEIN.ID_PROTEIN=CATEGORY_PROTEIN.ID_PROTEIN AND VISIBILITY>0 AND ID_ANALYSIS=$itemID AND ID_CATEGORY=$subjectID";
		}
		else { # Project hierarchy
			$visibleQuery="SELECT ID_PROTEIN,$fieldStrg FROM ANALYSIS_PROTEIN WHERE VISIBILITY>0 AND ID_ANALYSIS=$subjectID";
		}
	}
	elsif ($subjectITEM eq 'SAMPLE') {
		if ($filter) { # unclassified proteins
			$visibleQuery="SELECT ID_PROTEIN,ANALYSIS.$fieldStrg FROM ANALYSIS_PROTEIN,ANALYSIS WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND VISIBILITY>0 AND ID_SAMPLE=$itemID";
		}
		elsif ($classificationID > 0) { # user-defined
			$visibleQuery="SELECT ANALYSIS_PROTEIN.ID_PROTEIN,ANALYSIS.$fieldStrg FROM ANALYSIS_PROTEIN,CATEGORY_PROTEIN,ANALYSIS WHERE ANALYSIS.ID_ANALYSIS=ANALYSIS_PROTEIN.ID_ANALYSIS AND ANALYSIS_PROTEIN.ID_PROTEIN=CATEGORY_PROTEIN.ID_PROTEIN AND VISIBILITY>0 AND ANALYSIS.ID_SAMPLE=$itemID AND ID_CATEGORY=$subjectID";
		}
		else { # project hierarchy
			$visibleQuery="SELECT ID_PROTEIN,ANALYSIS.$fieldStrg FROM ANALYSIS_PROTEIN,ANALYSIS WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND VISIBILITY>0 AND ID_SAMPLE=$subjectID";
		}
	}
	elsif ($subjectITEM eq 'SPOT') {
		if ($filter) { # unclassified proteins
			$visibleQuery=qq|SELECT ID_PROTEIN,ANALYSIS.$fieldStrg
								FROM ANALYSIS_PROTEIN,ANALYSIS,SAMPLE
								WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE
									AND VISIBILITY>0 AND ID_SPOT=$itemID|;
		}
		elsif ($classificationID > 0) { # user-defined
			$visibleQuery=qq|SELECT ANALYSIS_PROTEIN.ID_PROTEIN,ANALYSIS.$fieldStrg
								FROM ANALYSIS_PROTEIN,SAMPLE,CATEGORY_PROTEIN,ANALYSIS
								WHERE ANALYSIS.ID_ANALYSIS=ANALYSIS_PROTEIN.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ANALYSIS_PROTEIN.ID_PROTEIN=CATEGORY_PROTEIN.ID_PROTEIN
									AND VISIBILITY>0 AND ID_SPOT=$itemID AND ID_CATEGORY=$subjectID|;
		}
		else { # project hierarchy
			$visibleQuery=qq|SELECT ID_PROTEIN,ANALYSIS.$fieldStrg
								FROM ANALYSIS_PROTEIN,ANALYSIS,SAMPLE
								WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE
									AND VISIBILITY>0 AND ID_SPOT=$subjectID|;
		}
	}
	elsif ($subjectITEM eq 'GEL2D') {
		if ($filter) { # unclassified proteins
			$visibleQuery=qq|SELECT ID_PROTEIN,ANALYSIS.$fieldStrg
								FROM ANALYSIS_PROTEIN,ANALYSIS,SAMPLE,SPOT
								WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT=SPOT.ID_SPOT
									AND VISIBILITY>0 AND ID_GEL2D=$itemID|;
		}
		elsif ($classificationID > 0) { # user-defined
			$visibleQuery=qq|SELECT ANALYSIS_PROTEIN.ID_PROTEIN,ANALYSIS.$fieldStrg
								FROM ANALYSIS_PROTEIN,SAMPLE,SPOT,CATEGORY_PROTEIN,ANALYSIS
								WHERE ANALYSIS.ID_ANALYSIS=ANALYSIS_PROTEIN.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT=SPOT.ID_SPOT AND ANALYSIS_PROTEIN.ID_PROTEIN=CATEGORY_PROTEIN.ID_PROTEIN
									AND VISIBILITY>0 AND ID_SPOT=$itemID AND ID_CATEGORY=$subjectID|;
		}
		else { # project hierarchy
			$visibleQuery=qq|SELECT ID_PROTEIN,ANALYSIS.$fieldStrg
								FROM ANALYSIS_PROTEIN,ANALYSIS,SAMPLE,SPOT
								WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT=SPOT.ID_SPOT
									AND VISIBILITY>0 AND ID_GEL2D=$subjectID|;
		}
	}
	elsif ($subjectITEM eq 'EXPERIMENT') {
		if ($filter) { # unclassified proteins
			$visibleQuery="SELECT ID_PROTEIN,ANALYSIS.$fieldStrg FROM ANALYSIS_PROTEIN,ANALYSIS,SAMPLE WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND VISIBILITY>0 AND ID_EXPERIMENT=$itemID";
		}
		elsif ($classificationID > 0) { # user-defined
			$visibleQuery=qq|SELECT ANALYSIS_PROTEIN.ID_PROTEIN,ANALYSIS.$fieldStrg
								FROM ANALYSIS,ANALYSIS_PROTEIN,SAMPLE,CATEGORY_PROTEIN
								WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ANALYSIS_PROTEIN.ID_PROTEIN=CATEGORY_PROTEIN.ID_PROTEIN
									AND VISIBILITY>0 AND ID_EXPERIMENT=$itemID AND ID_CATEGORY=$subjectID|;
 		}
		else { # project hierarchy
			$visibleQuery="SELECT ID_PROTEIN,ANALYSIS.$fieldStrg FROM ANALYSIS_PROTEIN,ANALYSIS,SAMPLE WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND VISIBILITY>0 AND ID_EXPERIMENT=$subjectID";
		}
	}
	elsif ($subjectITEM eq 'PROJECT') {
		if ($filter) { # unclassified proteins
			$visibleQuery="SELECT P.ID_PROTEIN,$fieldStrg
						   FROM PROTEIN P INNER JOIN ANALYSIS_PROTEIN AP ON P.ID_PROTEIN=AP.ID_PROTEIN INNER JOIN PROJECT PJ ON PJ.ID_PROJECT=P.ID_PROJECT INNER JOIN EXPERIMENT E ON E.ID_PROJECT=PJ.ID_PROJECT
						   WHERE VISIBILITY>0 AND P.ID_PROJECT=$itemID AND E.ID_EXPERIMENT NOT IN (
							 SELECT E2.ID_EXPERIMENT
							 FROM EXPERIMENT E2
                             INNER JOIN USER_EXPERIMENT_LOCK EU ON EU.ID_EXPERIMENT=E2.ID_EXPERIMENT
                             WHERE E2.ID_PROJECT=$itemID AND EU.ID_USER='$userID'
	   					   )";
		}
		elsif ($classificationID > 0) { # user-defined
			$visibleQuery=qq|SELECT P.ID_PROTEIN,$fieldStrg
							 FROM PROTEIN P INNER JOIN ANALYSIS_PROTEIN AP ON P.ID_PROTEIN=AP.ID_PROTEIN INNER JOIN PROJECT PJ ON PJ.ID_PROJECT=P.ID_PROJECT INNER JOIN EXPERIMENT E ON E.ID_PROJECT=PJ.ID_PROJECT INNER JOIN CATEGORY_PROTEIN CP ON CP.ID_PROTEIN=P.ID_PROTEIN
							 WHERE VISIBILITY>0 AND PJ.ID_PROJECT=$itemID AND ID_CATEGORY=$subjectID
							 AND E.ID_EXPERIMENT NOT IN (
							   SELECT E2.ID_EXPERIMENT
							   FROM EXPERIMENT E2
							   INNER JOIN USER_EXPERIMENT_LOCK EU ON EU.ID_EXPERIMENT=E2.ID_EXPERIMENT
							   WHERE E2.ID_PROJECT=$itemID AND EU.ID_USER='$userID'
							 )|;
		}
		else { # project hierarchy
			$visibleQuery="SELECT PROTEIN.ID_PROTEIN,$fieldStrg FROM PROTEIN,ANALYSIS_PROTEIN WHERE PROTEIN.ID_PROTEIN=ANALYSIS_PROTEIN.ID_PROTEIN AND VISIBILITY>0 AND ID_PROJECT=$subjectID";
		}
	}

	my $sthVis=$dbh->prepare($visibleQuery);
	my $sthAIT=$dbh->prepare("SELECT IDENTIFIER_TYPE,DB_RANK FROM ANALYSIS_DATABANK AD,DATABANK D WHERE ID_ANALYSIS=? AND AD.ID_DATABANK=D.ID_DATABANK");
	my ($addSpecificStg)=(param('sel_pepType') eq 'TYPIC')?' AND IS_SPECIFIC=1':'';
	#my $sthnbPep=$dbh->prepare("SELECT COUNT(ID_PEPTIDE) FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_PROTEIN=? AND ID_ANALYSIS=? AND PEP_BEG>0 AND PEP_END>0$addSpecificStg");


	####<Finding best analysis for each protein>####
	%timesFound=();
	my (%topMatchGroup,%protMatchGroup);
	my %virtualProt;
	my %analysisList;

	$sthVis->execute;
	my $protCount=0;
	if (!$filter) { # proteins in classification
		while (my ($protID,$anaID,$score,$numPep,$conf,$matchGr,$visibility,$coverage,$specificity,$dbRank)=$sthVis->fetchrow_array) { # $matchGr,$visibility defined if depthITEM=ANALYSIS
			if ($restrictListID) {
				next if (($listAction eq 'restrict' && !$restrictProteins{$protID}) || ($listAction eq 'exclude' && $restrictProteins{$protID}));
			}
			$dbRank=1 unless $dbRank;
			#my $nbPep; # nbPep <-> retrieve the number of proteotypic peptides or not, depending on user choice.
			#if ($addSpecificStg) {
			#	$sthnbPep->execute($protID,$anaID);
			#	($nbPep)=$sthnbPep->fetchrow_array;
			#}
			#else {$nbPep=$numPep;}
$analysisList{$anaID}=1;
			$timesFound{$protID}++;
			$virtualProt{$protID}{$anaID}=1 if $conf==0;
			@{$listProteins{$protID}}=();
			if(!$numPeptides{$protID}{'BA_ALL'} || $numPep > $numPeptides{$protID}{'BA_ALL'} || ($numPep == $numPeptides{$protID}{'BA_ALL'} && (!$bestScore{$protID} || $score > $bestScore{$protID}))) {
				$bestScore{$protID}=$score;
				$numPeptides{$protID}{'BA_ALL'} = $numPep; # always computed (needed for sort if no pepSort)
				$confLevel{$protID}=$conf;
				@{$bestAnalysis{$protID}}=($anaID,$dbRank);
				$pepCoverage{$protID}{'BA'}=$coverage if param('sel_cov');
				$protSpecificity{$protID}{'BA'}=$specificity if param('sel_spec');
				#$numPeptides{$protID}{'BA_ALL_TYPIC'}=$nbPep if param('sel_pepType') ne 'TYPIC';
				$numPeptides{$protID}{'BA_ALL_TYPIC'}=$numPep if param('sel_pepType') ne 'TYPIC';
			}
			if (param('sel_bestSpec')) {
				$protSpecificity{$protID}{'BI'}=$specificity if (!$protSpecificity{$protID}{'BI'} || $protSpecificity{$protID}{'BI'}<$specificity);
			}
			if (param('sel_matchGr')) { # depthITEM=ANALYSIS
				$topMatchGroup{$matchGr}=$protID if $visibility==2;
				$protMatchGroup{$protID}=$matchGr;
			}
			$numPeptides{$protID}{'CUM_ALL'}+=$numPep if param('sel_cumPepAll');
			#$numPeptides{$protID}{'CUM_ALL_TYPIC'}+=$nbPep if param('sel_cumPepAll');
$numPeptides{$protID}{'CUM_ALL_TYPIC'}+=$numPep if (param('sel_cumPepAll') && param('sel_pepType') ne 'TYPIC');
			#push @{$listAnalyses{$protID}},$anaID; # if param('sel_cumCov');
##########################################
			unless ($analysisIdentType{$anaID}) {
				$sthAIT->execute($anaID);
				while (my($identType,$dbRank)=$sthAIT->fetchrow_array) {
					$analysisIdentType{$anaID}{$dbRank}=$identType;
				}
				$analysisIdentType{$anaID}{1}='unknown' unless $analysisIdentType{$anaID};
			}
##########################################
			$classProteins{$protID}=1 if $classificationID>0;
			$protCount++;
			&printWaitingMsg('.','+=') unless $protCount % 1000;
		}
	}
	else { # unclassified proteins
		while (my ($protID,$anaID,$score,$numPep,$conf,$matchGr,$visibility,$coverage,$specificity,$dbRank)=$sthVis->fetchrow_array) {
			next if $classProteins{$protID}; # filter proteins in classification
			if ($restrictListID) {
				next if (($listAction eq 'restrict' && !$restrictProteins{$protID}) || ($listAction eq 'exclude' && $restrictProteins{$protID}));
			}
			$dbRank=1 unless $dbRank;
			#my $nbPep; # nbPep <-> retrieve the number of proteotypic peptides or not, depending on user choice.
			#if ($addSpecificStg) {
			#	$sthnbPep->execute($protID,$anaID);
			#	($nbPep)=$sthnbPep->fetchrow_array;
			#}
			#else {$bnPep=$numPep;}
$analysisList{$anaID}=1;
			$timesFound{$protID}++;
			$virtualProt{$protID}{$anaID}=1 if $conf==0;
			@{$listProteins{$protID}}=();
			if(!$numPeptides{$protID}{'BA_ALL'} || $numPep > $numPeptides{$protID}{'BA_ALL'} || ($numPep == $numPeptides{$protID}{'BA_ALL'} && (!$bestScore{$protID} || $score > $bestScore{$protID}))) {
				$bestScore{$protID}=$score;
				$numPeptides{$protID}{'BA_ALL'}= $numPep; # always computed (needed for sort if no pepSort)
				$confLevel{$protID}=$conf;
				@{$bestAnalysis{$protID}}=($anaID,$dbRank); # if param('sel_cov');
				$pepCoverage{$protID}{'BA'}=$coverage if param('sel_cov');
				$protSpecificity{$protID}{'BA'}=$specificity if param('sel_spec');
				#$numPeptides{$protID}{'BA_ALL_TYPIC'}=$nbPep;
				$numPeptides{$protID}{'BA_ALL_TYPIC'}=$numPep if param('sel_pepType') ne 'TYPIC';
			}
			if (param('sel_bestSpec')) {
				$protSpecificity{$protID}{'BI'}=$specificity if (!$protSpecificity{$protID}{'BI'} || $protSpecificity{$protID}{'BI'}<$specificity);
			}
			if (param('sel_matchGr')) { # depthITEM=ANALYSIS
				$topMatchGroup{$matchGr}=$protID if $visibility==2;
				$protMatchGroup{$protID}=$matchGr;
			}
			$numPeptides{$protID}{'CUM_ALL'}+=$numPep if param('sel_cumPepAll');
			#$numPeptides{$protID}{'CUM_ALL_TYPIC'}+=$nbPep if param('sel_cumPepAll');
$numPeptides{$protID}{'CUM_ALL_TYPIC'}+=$numPep if (param('sel_cumPepAll') && param('sel_pepType') ne 'TYPIC');
			#push @{$listAnalyses{$protID}},$anaID; # if param('sel_cumCov');
##########################################
			unless ($analysisIdentType{$anaID}) {
				$sthAIT->execute($anaID);
				while (my($identType,$dbRank)=$sthAIT->fetchrow_array) {
					$analysisIdentType{$anaID}{$dbRank}=$identType;
				}
				$analysisIdentType{$anaID}{1}='unknown' unless $analysisIdentType{$anaID};
			}
##########################################
			$protCount++;
			&printWaitingMsg('.','+=') unless $protCount % 1000;
		}
	}
	$sthVis->finish;
	$sthAIT->finish;
	#$sthnbPep->finish;

	####<Fetching remaining peptides statistics if 'TYPIC'>####
	if (param('sel_pepType') eq 'TYPIC') {
		my $sthAnaPep=$dbh->prepare("SELECT ID_PROTEIN,COUNT(ID_PEPTIDE) FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_ANALYSIS=? AND PEP_BEG>0 AND PEP_END>0 AND IS_SPECIFIC=1 GROUP BY ID_PROTEIN");
		foreach my $anaID (keys %analysisList) {
			$sthAnaPep->execute($anaID);
			while (my ($protID,$nbPep)=$sthAnaPep->fetchrow_array) {
				next unless $listProteins{$protID};
				$numPeptides{$protID}{'BA_ALL_TYPIC'}=$nbPep if $bestAnalysis{$protID}[0]==$anaID;
				$numPeptides{$protID}{'CUM_ALL_TYPIC'}+=$nbPep if param('sel_cumPepAll');
			}
		}
		$sthAnaPep->finish;
		&printWaitingMsg('.','+=');
	}

	####<Fetching info on selected proteins>####
	my $selProtString=(join ',',@protColumns);
	my $protSeqIdx;
	if ($selProtString=~/PROT_SEQ/) {
		foreach my $i (0..$#protColumns) {
			#if ($protColumns[$i] eq 'PROT_SEQ') {$protSeqIdx=$i+1; last;} # +1 because PROT_LENGTH is added to query
			if ($protColumns[$i] eq 'PROT_SEQ') {$protSeqIdx=$i; last;}
		}
	}
	#my $sthProtInfo=$dbh->prepare("SELECT PROT_LENGTH,$selProtString FROM PROTEIN WHERE ID_PROTEIN=?");
	#my $sthMPseq=$dbh->prepare("SELECT MP.PROT_SEQ FROM PROTEIN P,MASTER_PROTEIN MP WHERE P.ID_MASTER_PROTEIN=MP.ID_MASTER_PROTEIN AND ID_PROTEIN=?");
	#foreach my $protID (keys %listProteins) {
	#	$sthProtInfo->execute($protID);
	#	($protLength{$protID},@{$listProteins{$protID}})=$sthProtInfo->fetchrow_array;
	#	foreach my $pos (@checkPosDef) {$listProteins{$protID}[$pos]='-' unless $listProteins{$protID}[$pos];}
	#	if ($protSeqIdx && $listProteins{$protID}[$protSeqIdx] eq '+') { # '+' means replace with masterProt seq
	#		$sthMPseq->execute($protID);
	#		($listProteins{$protID}[$protSeqIdx])=$sthMPseq->fetchrow_array;
	#	}
	#}
	#$sthProtInfo->finish;
	#$sthMPseq->finish;
	my @proteins=keys %listProteins;
	my @noSeqProteins;
	while (my @subProtIDs=splice(@proteins,0,2000)) {
		my $sthProtInfo=$dbh->prepare("SELECT ID_PROTEIN,PROT_LENGTH,$selProtString FROM PROTEIN WHERE ID_PROTEIN IN (".join(',',@subProtIDs).")");
		$sthProtInfo->execute;
		while (my ($protID,$length,@info)=$sthProtInfo->fetchrow_array) {
			$protLength{$protID}=$length;
			@{$listProteins{$protID}}=@info;
			foreach my $pos (@checkPosDef) {$listProteins{$protID}[$pos]='-' unless $listProteins{$protID}[$pos];}
			push @noSeqProteins,$protID if ($protSeqIdx && $listProteins{$protID}[$protSeqIdx] eq '+');
		}
		$sthProtInfo->finish;
		&printWaitingMsg('.','+=');
	}
	if (scalar @noSeqProteins) {
		while (my @subProtIDs=splice(@noSeqProteins,0,2000)) {
			my $sthMPseq=$dbh->prepare("SELECT ID_PROTEIN,MP.PROT_SEQ FROM PROTEIN P,MASTER_PROTEIN MP WHERE P.ID_MASTER_PROTEIN=MP.ID_MASTER_PROTEIN AND ID_PROTEIN IN (".join(',',@subProtIDs).")");
			$sthMPseq->execute;
			while (my ($protID,$protSeq)=$sthMPseq->fetchrow_array) {
				$listProteins{$protID}[$protSeqIdx]=$protSeq;
			}
			$sthMPseq->finish;
		}
		&printWaitingMsg('.','+=');
	}

	if (param('sel_matchGr')) { # depthITEM=ANALYSIS
		my $sthMGrH=$dbh->prepare("SELECT P.ID_PROTEIN,IDENTIFIER,MATCH_GROUP FROM PROTEIN P,ANALYSIS_PROTEIN AP WHERE P.ID_PROTEIN=AP.ID_PROTEIN AND ID_ANALYSIS=$itemID AND VISIBILITY=2");
		my %matchGroups;
		$sthMGrH->execute;
		while (my ($protID,$ident,$matchGr)=$sthMGrH->fetchrow_array) {@{$matchGroups{$matchGr}}=($protID,$ident);}
		$sthMGrH->finish;
		
		my %externalHeadMGr;
		foreach my $protID (keys %listProteins) {
			if ($topMatchGroup{$protMatchGroup{$protID}}) {
				$proteinMGrHead{$protID}=$listProteins{$topMatchGroup{$protMatchGroup{$protID}}}[0]; # identifier
			}
			else { # head of MGr is listed in another block (user classification or unclassified proteins)
				unless ($externalHeadMGr{$protMatchGroup{$protID}}) {
					#$sthMGrH->execute($protMatchGroup{$protID}); # -> fetch it again
					#($externalHeadMGr{$protMatchGroup{$protID}}{'ID'},$externalHeadMGr{$protMatchGroup{$protID}}{'IDENT'})=$sthMGrH->fetchrow_array;
					($externalHeadMGr{$protMatchGroup{$protID}}{'ID'},$externalHeadMGr{$protMatchGroup{$protID}}{'IDENT'})=@{$matchGroups{ $protMatchGroup{$protID} }};
				}
				$proteinMGrHead{$protID}=$externalHeadMGr{$protMatchGroup{$protID}}{'IDENT'};
			}
		}
		&printWaitingMsg('.','+=');	
	}

	####<Fetching identifier mapping data>####
	my %proteinMappings;
	if (scalar @mappingColumns) {
		my %identID2code;
		my $identQueryStrg='';
		foreach my $code (@mappingColumns) {
			$identQueryStrg.=' OR ' if $identQueryStrg;
			$identQueryStrg.=" ID_IDENTIFIER=$identifierCodes{$code}[3]";
			$identID2code{$identifierCodes{$code}[3]}=$code;
		}
		#my $sthIV=$dbh->prepare("SELECT ID_IDENTIFIER,VALUE FROM PROTEIN P,MASTERPROT_IDENTIFIER MI WHERE P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND ID_PROTEIN=? AND ($identQueryStrg) ORDER BY IDENT_RANK");
		#foreach my $protID (keys %listProteins) {
		#	$sthIV->execute($protID);
		#	while (my ($identID,$value)=$sthIV->fetchrow_array) {
		#		push @{$proteinMappings{$protID}{$identID2code{$identID}}},$value; # ordered by rank
		#	}
		#	if ($proteinMappings{$protID} && $proteinMappings{$protID}{'GN'}) { # geneName/Symbol
		#		$proteinGeneSymbol{$protID}=$proteinMappings{$protID}{'GN'}[0];
		#	}
		#	else {$proteinGeneSymbol{$protID}='zzz';}
		#}
		#$sthIV->finish;
		my @proteins=keys %listProteins;
		while (my @subProtIDs=splice(@proteins,0,2000)) {
			my $sthIV=$dbh->prepare("SELECT ID_PROTEIN,ID_IDENTIFIER,GROUP_CONCAT(VALUE ORDER BY IDENT_RANK SEPARATOR ',') FROM PROTEIN P,MASTERPROT_IDENTIFIER MI
									WHERE P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND ID_PROTEIN IN (".join(',',@subProtIDs).") AND ($identQueryStrg)
									GROUP BY ID_PROTEIN,ID_IDENTIFIER");
			$sthIV->execute;
			while (my ($protID,$identID,$valueStrg)=$sthIV->fetchrow_array) {
				@{$proteinMappings{$protID}{$identID2code{$identID}}}=split(',',$valueStrg);
			}
			$sthIV->finish;
			&printWaitingMsg('.','+=');
		}		
		foreach my $protID (keys %listProteins) {
			if ($proteinMappings{$protID} && $proteinMappings{$protID}{'GN'}) { # geneName/Symbol
				$proteinGeneSymbol{$protID}=$proteinMappings{$protID}{'GN'}[0];
			}
			else {$proteinGeneSymbol{$protID}='zzz';}
		}
	}

	####<Cumulated coverage>####
	if (param('sel_cumCov')) {
		my $sthPep=$dbh->prepare("SELECT ID_PROTEIN,GROUP_CONCAT(PEP_BEG SEPARATOR ','),GROUP_CONCAT(PEP_END SEPARATOR ',') FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_ANALYSIS=? $addSpecificStg GROUP BY ID_PROTEIN");
		my %boundaryStatus;
		foreach my $anaID (keys %analysisList) {
			$sthPep->execute($anaID);
			while (my($protID,$begStrg,$endStrg)=$sthPep->fetchrow_array) {
				next unless $listProteins{$protID};
				my @begs=split(',',$begStrg);
				my @ends=split(',',$endStrg);
				foreach my $i (0..$#begs) {
					if ($virtualProt{$protID} && $virtualProt{$protID}{$anaID}) {
						$boundaryStatus{$protID}{abs($begs[$i])}++; # boundary are negative for ghost peptides
						$boundaryStatus{$protID}{abs($ends[$i])}--;
					}
					else {
						next if $begs[$i]<0; # ghost peptides
						$boundaryStatus{$protID}{$begs[$i]}++;
						$boundaryStatus{$protID}{$ends[$i]}--;
					}
				}
			}
			&printWaitingMsg('.','+=');	
		}
		$sthPep->finish;
		
		foreach my $protID (keys %boundaryStatus) {
			$pepCoverage{$protID}{'CUM'}=($protLength{$protID})? &promsMod::getCoverage($protLength{$protID},\%{$boundaryStatus{$protID}}) : 0;
		}
	}

	###<Fetching PhosphoRS results if needed>###
	#my %phosphoRS;
	#if (param('pepPhosphoRS')) {
	#	my %analyses;
	#	foreach my $protID (keys %listAnalyses){
	#		foreach my $anaID (@{$listAnalyses{$protID}}){
	#			$analyses{$anaID} = 1;
	#		}
	#	}
	#
	#	my $sthFile = $dbh->("SELECT DATA_FILE, VALID_STATUS FROM ANALYSIS WHERE ID_ANALYSIS=?");
	#	foreach my $anaID (keys %analyses){
	#		$sthFile->execute($anaID);
	#		my ($dataFile, $validStatus) = $sthFile->fetchrow_array;
	#		(my $prsFile = $dataFile) =~ s/^(.+)\..+/PRS_$1\.xml/;
	#		my $prsPath = ($validStatus == 2)? "$promsPath{peptide}/proj_$projectID/$anaID" : "$promsPath{valid}/ana_$anaID";
	#		if (-e "$prsPath/$prsFile") {
	#			$phosphoRS{$anaID} = phosphoRS->newFromResults(File => "$prsPath/$prsFile");
	#		}
	#
	#	}
	#	$sthFile->finish;
	#}


	###<Peptide data>###
	my (%pepSequences,%pepBoundaries,%flankingAA,%pepScores,%pepVarMods,%pepVarModsProt,%pepVarModsComment,%pepElutionTime,%firstPepStart,%pepTitle,%pepQValue,%pepCharge,%pepMass,%isSpecific,%pepComments,%pepPhosphoRS,%probPhosphoRS,%varModName);
	if (param('sel_pepNr') || param('sel_cumPepNr') || scalar @selPepColName || param('sel_protVmod')) {
		#my $sthPep2=$dbh->prepare("SELECT ID_PROTEIN,P.ID_PEPTIDE,PEP_SEQ,PEP_BEG,PEP_END,FLANKING_AA,SCORE,CHARGE,MR_OBS,IS_SPECIFIC,COMMENTS,DATA,SPEC_COUNT FROM PEPTIDE_PROTEIN_ATTRIB PPA,PEPTIDE P WHERE P.ID_PEPTIDE=PPA.ID_PEPTIDE AND P.ID_ANALYSIS=? $addSpecificStg ORDER BY ID_PROTEIN");
		my $sthPep2=$dbh->prepare("SELECT ID_PROTEIN,P.ID_PEPTIDE,QUERY_NUM,PEP_SEQ,ELUTION_TIME,GROUP_CONCAT(PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),PEP_BEG,PEP_END,FLANKING_AA,SCORE,CHARGE,MR_OBS,IS_SPECIFIC,COMMENTS,DATA,SPEC_COUNT
											FROM PEPTIDE P
											LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
											INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
											WHERE P.ID_ANALYSIS=? $addSpecificStg GROUP BY P.ID_PEPTIDE,ID_PROTEIN ORDER BY ID_PROTEIN");
		if (param('sel_seqPep') || param('sel_pepNr')) { # best analysis
			foreach my $anaID (keys %analysisList) {
				$sthPep2->execute($anaID);
				my ($prsOutputFile,%dataPhosphoRS);
				if (param('sel_phosphoRS')) {
					($prsOutputFile,my $validStatus)=$dbh->selectrow_array("SELECT DATA_FILE,VALID_STATUS FROM ANALYSIS WHERE ID_ANALYSIS=$anaID");
					$prsOutputFile =~ s/\.\w{3}$/\.xml/;
					if ($validStatus == 2) {
						$prsOutputFile = (-e "$promsPath{peptide}/proj_$projectID/ana_$anaID/PRS_$prsOutputFile")? "$promsPath{peptide}/proj_$projectID/ana_$anaID/PRS_$prsOutputFile" : "$promsPath{peptide}/proj_$projectID/ana_$anaID/PRS_ana_$anaID.xml";
					}
					else {
						$prsOutputFile = (-e "$promsPath{valid}/ana_$anaID/PRS_$prsOutputFile")?"$promsPath{valid}/ana_$anaID/PRS_$prsOutputFile" : "$promsPath{valid}/ana_$anaID/PRS_ana_$anaID.xml";
					}
				}
				my $prevProtID=0;
				my (%bestPepScore,%refPeptide,%distinctPep); # in case no sel_modPep => do not duplicate peptides with same sequence
				while (my ($protID,$pepID,$qNum,$pepSeq,$elutionTime,$modCode,$beg,$end,$flankingAA,$score,$charge,$massObs,$isSpec,$comments,$data,$specCount)=$sthPep2->fetchrow_array) {
					next if (!$listProteins{$protID} || $bestAnalysis{$protID}[0] != $anaID);
					if ($protID != $prevProtID) { # reset for each prot
						$numPeptides{$prevProtID}{'BA_NR'}=scalar keys %distinctPep if (param('sel_pepNr') && $listProteins{$prevProtID});
						%bestPepScore=%refPeptide=%distinctPep=();
						$prevProtID=$protID;
					}
					if ($virtualProt{$protID} && $virtualProt{$protID}{$bestAnalysis{$protID}[0]}) { # virtual protein
						$beg=abs($beg);
						$end=abs($end);
					}
					elsif ($beg<0) {next;} # ghost peptides
					#my $vMod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$bestAnalysis{$protID}[0],$pepSeq); #!!! TODO: improve this. TOO SLOW!!!!!!!!!!!!!
					my $vMod=&promsMod::decodeVarMod($dbh,$pepSeq,$modCode,\%varModName);
					$vMod='-' unless $vMod;
					
					my $vModProt=&promsMod::decodeVarMod($dbh,$pepSeq,$modCode,\%varModName, $beg-1);
					$vModProt = '-' unless $vModProt;
					
					#$vMod=~s/([:.])-1/$1?/g; # -1 -> ? (position ambiguity)
					$comments='-' unless $comments;
					if (param('sel_pepNr')) {$distinctPep{"$pepSeq$vMod"}=1;}
					if (param('sel_seqPep')) {
						$vMod=~s/^ \+ //;
						#my $pepStrg=(param('sel_modPep'))? "$pepSeq:$vMod" : $pepSeq;
						my $pepStrg="$pepSeq:$vMod";
						if (param('sel_noDupPep') && $bestPepScore{$pepStrg}) { # || !param('sel_modPep')
							if ($score>$bestPepScore{$pepStrg}) {
								$bestPepScore{$pepStrg}=$score;
								$pepScores{$protID}{$refPeptide{$pepStrg}}=$score; # always keep highest score
							}
							next;
						}
						$pepTitle{$pepID}{'analysis'}=$bestAnalysis{$protID}[0];
						$refPeptide{$pepStrg}=$pepID;
						$bestPepScore{$pepStrg}=$score;
						my ($flkNter,$flkCter)=('','');
						if (param('sel_flkPep')) {
							($flkNter,$flkCter)=($flankingAA)? split(//,$flankingAA) : ('?','?');
							$flkNter.='.';
							$flkCter=".$flkCter";
						}
						$pepSequences{$protID}{$pepID}=$flkNter.$pepSeq.$flkCter;
						$pepBoundaries{$protID}{$pepID}{$beg}=$end; # always computed
						$pepScores{$protID}{$pepID}=$score if param('sel_scorePep');

						if (param('sel_modPep') || param('sel_modPepProt') || param('sel_protVmod')) {
							($pepVarMods{$protID}{$pepID}= $vMod) =~ s/\[[0-9,]+\]//g;
							($pepVarModsProt{$protID}{$pepID} = $vModProt) =~ s/\[[0-9,]+\]//g;
							$pepVarModsComment{$protID}{$pepID} = &generateVarModSummary($vModProt, $pepSeq, $beg-1);
						}
						
						if (param('sel_etPep')) {
							my @elutionTimeParsed = ( $elutionTime =~ /et([0-9.]+);?/);
							$elutionTime = (scalar @elutionTimeParsed == 1)? $elutionTimeParsed[0] : '-';
							$pepElutionTime{$pepID} = $elutionTime;
						}
						
						$pepCharge{$pepID}=$charge if param('sel_chargePep');
						if (param('sel_qvalPep') ){
							if ($data && $data =~ /QVAL=([^,]+)/) {
								$pepQValue{$protID}{$pepID}=$1;
							}
							else {
								$pepQValue{$protID}{$pepID}='-';
							}
						}
						$pepMass{$pepID}=$massObs if param('sel_massPep');
						$pepComments{$pepID}=$comments if param('sel_commentPep');
						if (param('sel_isSpec')) {$isSpecific{$protID}{$pepID}=($isSpec)? 'Yes' : 'No';}
					}
					if (param('sel_phosphoRS')) {
						if ($data && $data =~ /PRS=([^#]+)/) {
							$pepPhosphoRS{$pepID} = &getPhosphoRsString($1);
							#$probPhosphoRS{$pepID} = &getPhosphoRsProbString($pepID,$pepSeq,$dbh,$prsOutputFile);
							$dataPhosphoRS{$qNum}=$pepID;
							
						} else {
							$pepPhosphoRS{$pepID} = '-';
						}
						$probPhosphoRS{$pepID} = '-'; # default
					}
					if (param('sel_pepPSM')) {
						$specCount=0 unless defined($specCount);
						$protPSM{$protID}+=$specCount;
					}
				}
				$numPeptides{$prevProtID}{'BA_NR'}=scalar keys %distinctPep if (param('sel_pepNr') && $listProteins{$prevProtID}); # for last protID in while loop
				
				#<Fetch all PRS probabilities
				if (param('sel_phosphoRS') && scalar keys %pepPhosphoRS) {
					&getPhosphoRsProbString($prsOutputFile,\%pepPhosphoRS,\%dataPhosphoRS,\%probPhosphoRS);
				}

				&printWaitingMsg('.','+=');
			}
		}
		if (param('sel_cumSeqPep') || param('sel_cumPepNr') || param('sel_protVmod')) { # cumulated
			foreach my $anaID (keys %analysisList) {
				$sthPep2->execute($anaID);
				my $tot='';
				my ($prsOutputFile,%dataPhosphoRS);
				if (param('sel_phosphoRS')) {
					($prsOutputFile,my $validStatus)=$dbh->selectrow_array("SELECT DATA_FILE,VALID_STATUS FROM ANALYSIS WHERE ID_ANALYSIS=$anaID");
					$prsOutputFile =~ s/\.\w{3}$/\.xml/;
					if ($validStatus == 2) {
						$prsOutputFile = (-e "$promsPath{peptide}/proj_$projectID/ana_$anaID/PRS_$prsOutputFile")? "$promsPath{peptide}/proj_$projectID/ana_$anaID/PRS_$prsOutputFile" : "$promsPath{peptide}/proj_$projectID/ana_$anaID/PRS_ana_$anaID.xml";
					}
					else {
						$prsOutputFile = (-e "$promsPath{valid}/ana_$anaID/PRS_$prsOutputFile")?"$promsPath{valid}/ana_$anaID/PRS_$prsOutputFile" : "$promsPath{valid}/ana_$anaID/PRS_ana_$anaID.xml";
					}
				}
				my $prevProtID=0;
				my (%bestPepScore,%refPeptide,%distinctPep);
				while (my ($protID,$pepID,$qNum,$pepSeq,$elutionTime,$modCode,$beg,$end,$flankingAA,$score,$charge,$massObs,$isSpec,$comments,$data,$specCount)=$sthPep2->fetchrow_array) {
					next unless $listProteins{$protID};
					if ($protID != $prevProtID) { # reset for each prot
						$numPeptides{$prevProtID}{'CUM_NR'}=scalar keys %distinctPep if (param('sel_cumPepNr') && $listProteins{$prevProtID});
						%bestPepScore=%refPeptide=%distinctPep=();
						$prevProtID=$protID;
					}
					#my $vMod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$anaID,$pepSeq); #!!! TODO: improve this. TOO SLOW!!!!!!!!!!!!!
					my $vMod=&promsMod::decodeVarMod($dbh,$pepSeq,$modCode,\%varModName);
					$vMod='-' unless $vMod;
					
					my $vModProt=&promsMod::decodeVarMod($dbh,$pepSeq,$modCode,\%varModName, $beg-1);
					$vModProt = '-' unless $vModProt;
					
					$comments='-' unless $comments;
					if (param('sel_cumPepNr')) {$distinctPep{"$pepSeq$vMod"}=1;}
					if (param('sel_cumSeqPep')) {
						$vMod=~s/^ \+ //;
						#my $pepStrg=(param('sel_modPep'))? "$pepSeq:$vMod" : $pepSeq;
						my $pepStrg="$pepSeq:$vMod";
						if (param('sel_noDupPep') && $bestPepScore{$pepStrg}) { # || !param('sel_modPep')
							if ($score>$bestPepScore{$pepStrg}) {
								$bestPepScore{$pepStrg}=$score;
								$pepScores{$protID}{$refPeptide{$pepStrg}}=$score; # always keep highest score
							}
							next;
						}
						$pepTitle{$pepID}{'analysis'}=$anaID;
						$refPeptide{$pepStrg}=$pepID;
						$bestPepScore{$pepStrg}=$score;
						my ($flkNter,$flkCter)=('','');
						if (param('sel_flkPep')) {
							($flkNter,$flkCter)=($flankingAA)? split(//,$flankingAA) : ('?','?');
							$flkNter.='.';
							$flkCter=".$flkCter";
						}
						$pepSequences{$protID}{$pepID}=$flkNter.$pepSeq.$flkCter;
						$pepBoundaries{$protID}{$pepID}{$beg}=$end ;#if param('sel_posPep'); #order by pos !
						$pepScores{$protID}{$pepID}=$score if param('sel_scorePep');
						
						if (param('sel_modPep') || param('sel_modPepProt') || param('sel_protVmod')) {
							($pepVarMods{$protID}{$pepID}= $vMod) =~ s/\[[0-9,]+\]//g;
							($pepVarModsProt{$protID}{$pepID} = $vModProt) =~ s/\[[0-9,]+\]//g;
							$pepVarModsComment{$protID}{$pepID} = &generateVarModSummary($vModProt, $pepSeq, $beg-1);
						}
						
						if (param('sel_etPep')) {
							my @elutionTimeParsed = ( $elutionTime =~ /et([0-9.]+);?/);
							$elutionTime = (scalar @elutionTimeParsed == 1) ? $elutionTimeParsed[0] : '-';
							$pepElutionTime{$pepID} = $elutionTime;
						}
						
						$pepCharge{$pepID}=$charge if param('sel_chargePep');
						if (param('sel_qvalPep')) {
							if ($data && $data =~ /QVAL=([^,]+)/) {
								$pepQValue{$protID}{$pepID}=$1;
							}
							else {
								$pepQValue{$protID}{$pepID}='-';
							}
						}
						$pepMass{$pepID}=$massObs if param('sel_massPep');
						$pepComments{$pepID}=$comments if param('sel_commentPep');
						if (param('sel_isSpec')) {$isSpecific{$protID}{$pepID}=($isSpec)? 'Yes' : 'No';}
						if (param('sel_phosphoRS')) {
							if ($data && $data =~ /PRS=([^#]+)/) {
								$pepPhosphoRS{$pepID} = &getPhosphoRsString($1);
								#$probPhosphoRS{$pepID} = &getPhosphoRsProbString($pepID,$pepSeq,$dbh,$prsOutputFile);
								$dataPhosphoRS{$qNum}=$pepID;
							}
							else {
								$pepPhosphoRS{$pepID} = '-';
							}
							$probPhosphoRS{$pepID} = '-'; # default
						}
					}
					if (param('sel_protVmod')) {
						$pepVarMods{$protID}{$pepID}=$vMod;
						$pepBoundaries{$protID}{$pepID}{$beg}=$end;
					}
					if (param('sel_pepPSM')) {
						$specCount=0 unless defined($specCount);
						$protPSM{$protID}+=$specCount;
					}
					$numPeptides{$protID}{'CUM_ALL_TYPIC'}+=1;
				}
				$numPeptides{$prevProtID}{'CUM_NR'}=scalar keys %distinctPep if (param('sel_cumPepNr') && $listProteins{$prevProtID}); # for last protID in while loop
				
				#<Fetch all PRS probabilities
				if (param('sel_phosphoRS') && scalar keys %pepPhosphoRS) {
					&getPhosphoRsProbString($prsOutputFile,\%pepPhosphoRS,\%dataPhosphoRS,\%probPhosphoRS);
				}
				
				&printWaitingMsg('.','+=');	
			}
		}
		$sthPep2->finish;
		
		##<finding 1st start for all pep/prot>## (used to sort peptides in list)
		foreach my $protID (keys %pepBoundaries) {
			foreach my $pepID (keys %{$pepBoundaries{$protID}}) {
				$firstPepStart{$protID}{$pepID}=(sort{$a<=>$b} keys %{$pepBoundaries{$protID}{$pepID}})[0];
			}
		}

		if (param('sel_titlePep')) {
			my %tempAnaPep;
			foreach my $pepID (keys (%pepTitle)) {
				$tempAnaPep{$pepTitle{$pepID}{'analysis'}}{$pepID}=1;
			}

			#my $sthquery = $dbh->prepare("SELECT PEPTIDE.QUERY_NUM, ANALYSIS.FILE_FORMAT FROM PEPTIDE, ANALYSIS WHERE PEPTIDE.ID_PEPTIDE = ? AND ANALYSIS.ID_ANALYSIS IN (SELECT PEPTIDE.ID_ANALYSIS FROM PEPTIDE WHERE ID_PEPTIDE = ?)");
			my $sthAna=$dbh->prepare("SELECT DATA_FILE,FILE_FORMAT FROM ANALYSIS WHERE ID_ANALYSIS=?");
			#my $sthPep=$dbh->prepare("SELECT QUERY_NUM FROM PEPTIDE WHERE ID_PEPTIDE=?");

			foreach my $anaID (keys (%tempAnaPep)) {
				$sthAna->execute($anaID);
				my ($fileName,$fileFormat)=$sthAna->fetchrow_array;
				my %queryNumList; # queryNumList {query_Number} = pep_ID
				#foreach my $pepID (keys (%{$tempAnaPep{$anaID}})) {
				#	$sthPep->execute($pepID);
				#	my ($queryNum)=$sthPep->fetchrow_array;
				#	$queryNumList{$queryNum}=$pepID;
				#	$pepTitle{$pepID}{'title'}='Not found'; # default
				#}
				my @pepList=(keys (%{$tempAnaPep{$anaID}}));
				while (my @subPepIDs=splice(@pepList,0,2000)) {
					my $sthPep=$dbh->prepare("SELECT ID_PEPTIDE,QUERY_NUM FROM PEPTIDE WHERE ID_PEPTIDE IN (".join(',',@subPepIDs).")");
					$sthPep->execute;
					while (my ($pepID,$queryNum)=$sthPep->fetchrow_array) {
						$queryNumList{$queryNum}=$pepID;
						$pepTitle{$pepID}{'title'}='Not found'; # default
					}
					$sthPep->finish;
				}

				my $onQuery=0;
				my $localQueryNumber;
				my $dataFile;
				#my $fileName=sprintf "P%06d",$anaID;
				#if ($fileFormat eq 'MASCOT.DAT') {$dataFile="$promsPath{peptide}/proj_$projectID/$fileName.dat";}
				#elsif ($fileFormat eq 'PHENYX.XML' || $fileFormat eq 'MASCOT.XML') {$dataFile="$promsPath{peptide}/proj_$projectID/$fileName.pgf";}
				#$dataFile="$promsPath{peptide}/proj_$projectID/ana_$anaID/$fileName";
				if ($fileFormat eq 'MASCOT.DAT') {$dataFile=~s/\.dat\Z/_min\.dat/ unless -e $dataFile;} # assumes Fxxxxxx_min.dat file
				elsif ($fileFormat=~/(PHENYX|MASCOT)\.XML/) {$dataFile=~s/\.xml\Z/\.pgf/;}

				open (DESC, $dataFile) || next; # skip if problem

				while (my $line=<DESC>) {
					if ($line=~/name="query(\d+)"/) {
						$localQueryNumber = $1;
						if (exists $queryNumList{$localQueryNumber}) {$onQuery=1;}
						else {$onQuery=0;}
					}
					elsif ($onQuery==1 && $line=~/^--gc0p4Jq0M2Yt08jU534c0p/) {$onQuery=0;}
					elsif ($onQuery==1 && $line=~/title=(.+)/) {
						my $title=$1;
						$title=~s/%([\dA-Fa-f][\dA-Fa-f])/pack ("C", hex($1))/eg;
						$pepTitle{$queryNumList{$localQueryNumber}}{'title'}=$title;
					}
				}
				close DESC;
				&printWaitingMsg('.','+=');	
			}
			#$sthPep->finish;
			$sthAna->finish;
		}
	}

	# Fetching peptide var mods for each protein
	my %protVarMods;
	if (param('sel_protVmod')) {
		foreach my $protID (keys %listProteins) {
			my %varMod;
			foreach my $pepID (keys %{$pepVarMods{$protID}}){
				my $vmodString = $pepVarMods{$protID}{$pepID};
				next if $vmodString eq '-';
				$vmodString =~ s/^ \+ //;
				my @begs = keys %{$pepBoundaries{$protID}{$pepID}};
				my $beg = $begs[0];
				foreach my $vmod (split / \+ /, $vmodString){
					next if (($vmod =~ /N-term/ || $vmod =~ /C-term/) && $vmod !~ /Protein/);
					if ($vmod =~ /Protein/) {
						my ($name, $aa) = ($vmod =~ /^(\w+) \((.+)\)/);
						$varMod{$name}{$aa}{0} = 1;
					}
					else {
						my ($name, $aa, $posString) = ($vmod =~ /^(\w+) \((.+):(.+)\)/);
						foreach my $pos (split /\./,$posString){
							my $protPos=($pos=~/\D/)? $pos : $pos+$beg-1; # eg '?'
							$varMod{$name}{$aa}{$protPos} = 1;
						}
					}
				}
			}
			foreach my $name (sort {$a cmp $b} keys %varMod){
				next unless $relevantPTMs{$name};
				foreach my $aa (sort {$a cmp $b} keys %{$varMod{$name}}){
					if ($aa =~ /Protein/) {
						$aa =~ s/Protein //;
						$protVarMods{$protID} .= "$name ($aa) + ";
					} else {
						$protVarMods{$protID} .= "$name ($aa:" . join('.', sort {$a <=> $b} keys %{$varMod{$name}{$aa}}) . ') + ';
					}
				}
			}
			if ($protVarMods{$protID}) {
				$protVarMods{$protID} =~ s/ \+ $//;
			} else {
				$protVarMods{$protID} = '-';
			}
			push @{$listProteins{$protID}}, $protVarMods{$protID};
		}
	}

	################################
	####<Printing data in table>####
	################################
	&printWaitingMsg('/','+=');
	
	####<No proteins>####
	if (scalar (keys %listProteins)==0) {
		if ($exportFormat eq 'XLS') {
			$worksheet->merge_range(++$xlsRow,0,$xlsRow,$colSpan-1,'No proteins',$format{'mergeText'});
		}
		else {print HTML "<TR><TD colspan=$colSpan>&nbsp;No proteins</TD></TR>\n";}
		#print HTML "<TR><TD colspan=$colSpan>&nbsp</TD></TR>\n";
		return;
	}

	####<Columns name>####
	my $itemType=&promsMod::getItemType($subjectITEM);
	my ($addColMatchPeptides)=(param('sel_pepType') eq 'TYPIC')?' Proteotypic ':' ';
	if ($exportFormat eq 'XLS') {
		$xlsRow++;
		my $xlsCol=0;
		my $pepDataCol = 0;
		my $pepMatchCol;
		foreach my $colName (@selColName) {
			$colName=~s/#ITEM#/$itemType/;
			if ($colName eq "Matching${addColMatchPeptides}peptides") {
				if ($colSpanNumPep > 1) {$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow,$xlsCol+$colSpanNumPep-1,$colName,$itemFormat{'mergeHeader'});}
				else {$worksheet->write_string($xlsRow,$xlsCol,$colName,$itemFormat{'header'});}
				$pepMatchCol=$xlsCol-1;
				$xlsCol+=$colSpanNumPep;
			}
			elsif ($colName=~/^Peptide data/) {
				if ($colSpanPep > 1) {$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow,$xlsCol+$colSpanPep-1,$colName,$itemFormat{'mergeHeader'});}
				else {$worksheet->write_string($xlsRow,$xlsCol,$colName,$itemFormat{'header'});}
				$pepDataCol=$xlsCol-1;
				$xlsCol+=$colSpanPep;
			}
			else { # no col merge
				if ($rowSpan==2) {$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow+$rowSpan-1,$xlsCol,$colName,$itemFormat{'mergeHeader'});} # row merge
				else {$worksheet->write_string($xlsRow,$xlsCol,$colName,$itemFormat{'header'});}
				$xlsCol++;
			}
		}
		if ($rowSpan==2) {
			$xlsRow++;
			foreach my $colName (@numPepColumns) {$colName=~s/#ITEM#/$itemType/; $worksheet->write_string($xlsRow,++$pepMatchCol,$colName,$itemFormat{'header'});}
			foreach my $colName (@selPepColName) {$worksheet->write_string($xlsRow,++$pepDataCol,$colName,$itemFormat{'header'});}
		}
		$xlsRow++;
	}
	else { # HTML
		print HTML "<TR>\n";
		foreach my $colName (@selColName) {
			$colName=~s/#ITEM#/$itemType/;
			if ($colName eq "Matching${addColMatchPeptides}peptides") {print HTML "<TH colspan=$colSpanNumPep bgcolor=\"$color1\">$colName</TH>\n";}
			elsif ($colName=~/^Peptide data/) {print HTML "<TH colspan=$colSpanPep bgcolor=\"$color1\">$colName</TH>\n";}
			else {print HTML "<TH rowspan=$rowSpan bgcolor=\"$color1\">$colName</TH>\n";}
		}
		print HTML "</TR>\n";
		if ($rowSpan==2) {
			print HTML "<TR>\n";
			foreach my $colName (@numPepColumns) {$colName=~s/#ITEM#/$itemType/; print HTML "<TH bgcolor=\"$color1\">$colName</TH>\n";}
			foreach my $colName (@selPepColName) {print HTML "<TH bgcolor=\"$color1\">$colName</TH>\n";}
			print HTML "</TR>\n";
		}
	}

	####<Listing proteins>####
	my $sortPep='';
	if ($view=~/pep/i) { # sort by peptides
		$sortPep=($view eq 'pepAll')? 'BA_ALL' : ($view eq 'pepNr')? 'BA_NR' : ($view eq 'cumPepAll')? 'CUM_ALL' : 'CUM_NR';
	}
	else {
		$sortPep=(param('sel_cumPepAll'))? 'CUM_ALL' : (param('sel_cumPepNr'))? 'CUM_NR' : (param('sel_pepNr'))? 'BA_NR' : 'BA_ALL';
	}
	my $rowColor=$rowColor1;
	$protCount=0;
	foreach my $protID (sort{&sortRule($view,$sortPep)} keys %listProteins) {
		$protCount++;
		my $numListedPep=scalar keys %{$pepSequences{$protID}};
		my $xlsEndRow=$xlsRow+$numListedPep-1;
		my $rowSpanStrg=($pepSequences{$protID})? "rowspan=$numListedPep" : '';
		print HTML "<TR valign=top>\n" if $exportFormat ne 'XLS';

		##<Protein data
		my $xlsCol=0;
		my $linkURL; # for NCBI & SWISSPROT only
		my $identType=$analysisIdentType{$bestAnalysis{$protID}[0]}{$bestAnalysis{$protID}[1]};
		if ($identType=~/NCBI_ALL|GI_ACCESSION/) {
			my ($modExtIdent)=($listProteins{$protID}[0]=~/gi\|(\d+)/);
			$linkURL="http://www.ncbi.nlm.nih.gov/protein/$modExtIdent?report=genpept" if $modExtIdent;
		}
		elsif ($identType=~/UNIPROT_ALL|UNIPROT_ID/) {
			my ($modExtIdent)=($listProteins{$protID}[0]=~/(\w+_\w+)/);
			($linkURL=$identifierCodes{'ID'}[2])=~s/XXXX/$modExtIdent/ if $modExtIdent;
		}
		elsif ($identType eq 'AC') {($linkURL=$identifierCodes{'AC'}[2])=~s/XXXX/$listProteins{$protID}[0]/;}
		@{$protIdentifierLinks{$listProteins{$protID}[0]}}=($linkURL,0) unless $protIdentifierLinks{$listProteins{$protID}[0]}; # linkURL can be undef [1]= max numb mapped ident per code

		if ($exportFormat eq 'XLS') {
			my $statusStrg=($confLevel{$protID}==0)? 'VP' : ($confLevel{$protID}==1)? 'BC' : '';
			if ($numListedPep > 1) {
				$worksheet->merge_range($xlsRow,$xlsCol,$xlsEndRow,$xlsCol,$protCount,$format{'mergeNumber'});
				$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$listProteins{$protID}[0],$format{'mergeText'});
				foreach my $colData (@{$listProteins{$protID}}[1..$#{$listProteins{$protID}}]) {$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$colData,$format{"mergeIdent$statusStrg"});}
			}
			else {
				$worksheet->write_string($xlsRow,$xlsCol,$protCount,$format{'number'});
				$worksheet->write_string($xlsRow,++$xlsCol,$listProteins{$protID}[0],$format{'text'});
				foreach my $colData (@{$listProteins{$protID}}[1..$#{$listProteins{$protID}}]) {$worksheet->write_string($xlsRow,++$xlsCol,$colData,$format{"ident$statusStrg"});}
			}
		}
		else { # HTML
			print HTML "<TD $rowSpanStrg bgcolor=\"$rowColor\" align=\"right\">$protCount</TD>";
			my $linkedIdentifier=($linkURL)? "<A href=\"$linkURL\" target=\"_blank\">$listProteins{$protID}[0]</A>" : $listProteins{$protID}[0];
			if ($confLevel{$protID}==0) { # virtual prot
				print HTML "<TD $rowSpanStrg bgcolor=\"$rowColor\"><I>$linkedIdentifier</I></TD>";
				$badConfidence=1;
			}
			elsif ($confLevel{$protID}==1) { # bad confidence
				print HTML "<TD $rowSpanStrg bgcolor=\"$rowColor\"><FONT color=\"DimGray\">$linkedIdentifier</FONT></TD>";
				$badConfidence=1;
			}
			else {print HTML "<TD $rowSpanStrg bgcolor=\"$rowColor\">$linkedIdentifier</TD>";}
			my $colPos=1;
			foreach my $colData (@{$listProteins{$protID}}[1..$#{$listProteins{$protID}}]) {
				print HTML "<TD $rowSpanStrg bgcolor=\"$rowColor\">$colData</TD>";
				$colPos++;
			}
			print HTML "\n";
		}

		##<Mapping data
		if (scalar @mappingColumns) {
			unless ($mappingLinks{$listProteins{$protID}[0]}) { # compute URL for each mapping type only once per protID
				if ($proteinMappings{$protID}) { # mappings exist for prot
					foreach my $code (keys %{$proteinMappings{$protID}}) {
						my $resrcURL=$identifierCodes{$code}[2];
						next unless $resrcURL; # no link for code (eg. GN)
						my $numIdent=scalar (@{$proteinMappings{$protID}{$code}});
						$protIdentifierLinks{$listProteins{$protID}[0]}[1]=$numIdent if $numIdent > $protIdentifierLinks{$listProteins{$protID}[0]}[1];
						foreach my $identValue (@{$proteinMappings{$protID}{$code}}) {
							(my $identURL=$resrcURL)=~s/XXXX/$identValue/;
							if ($code eq 'RefSeq') {
								my $extraURLstrg=($identValue=~/^NP_/)? 'entrez/viewer.fcgi?db=protein&id=' : 'nuccore/';
								$identURL=~s/VVVV/$extraURLstrg/;
							}
							elsif ($code eq 'UniGene') {
								my ($org,$val)=split('\.',$identValue);
								$identURL=~s/YYYY/$org/;
								$identURL=~s/ZZZZ/$val/;
							}
							push @{$mappingLinks{$listProteins{$protID}[0]}{$code}},[$identURL,$identValue]; # @{___{protIdent}{mapIdentCode}}=([$url,$Xident],[...])
						}
					}
				}
			}

			foreach my $i (0..$#mappingColumns) {
				my $code=$mappingColumns[$i];
				if ($proteinMappings{$protID}{$code}) { # mapped ident exists
					if ($exportFormat eq 'XLS') { # links not handled on this worksheet
						if ($numListedPep > 1) {
							$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,join("\n",@{$proteinMappings{$protID}{$code}}),$format{'mergeText'});
						}
						else {$worksheet->write_string($xlsRow,++$xlsCol,join("\n",@{$proteinMappings{$protID}{$code}}),$format{'text'});}
					}
					else { # HTML
						if ($mappingLinks{$listProteins{$protID}[0]}{$code}) { # link for mapped ident
							my $listURL='';
							foreach my $refLink (@{$mappingLinks{$listProteins{$protID}[0]}{$code}}) {
								$listURL.='<BR>' if $listURL;
								$listURL.='<A href="'.$refLink->[0].'" target="_blank">'.$refLink->[1].'</A>';
							}
							print HTML "<TD $rowSpanStrg bgcolor=\"$rowColor\">$listURL</TD>";
						}
						else { # no link for this code
							print HTML "<TD $rowSpanStrg bgcolor=\"$rowColor\">",join('<BR>',@{$proteinMappings{$protID}{$code}}),"</TD>";
						}
					}
				}
				else { # no mapping
					if ($exportFormat eq 'XLS') {
						if ($numListedPep > 1) {$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,'-',$format{'mergeText'});}
						else {$worksheet->write_string($xlsRow,++$xlsCol,'-',$format{'text'});}
					}
					else {print HTML "<TD $rowSpanStrg bgcolor=\"$rowColor\">&nbsp;</TD>";}
				}
			}
			print HTML "\n" if $exportFormat ne 'XLS';
		}

		##<Match data
		if ($exportFormat eq 'XLS') {
			if ($numListedPep > 1) {
				$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$proteinMGrHead{$protID},$format{'mergeNumber'}) if param('sel_matchGr');
				$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$timesFound{$protID},$format{'mergeNumber'}) if param('sel_times');
				$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$bestScore{$protID},$format{'mergeNumber'}) if param('sel_score');
				$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$protPSM{$protID},$format{'mergeNumber'}) if param('sel_pepPSM');
				$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$numPeptides{$protID}{'BA_ALL'},$format{'mergeNumber'}) if param('sel_pepAll') && param('sel_pepType') eq 'ALL';
				$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$numPeptides{$protID}{'BA_ALL_TYPIC'},$format{'mergeNumber'}) if param('sel_pepAll') && param('sel_pepType') eq 'TYPIC';
				$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$numPeptides{$protID}{'BA_NR'},$format{'mergeNumber'}) if param('sel_pepNr');
				$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$numPeptides{$protID}{'CUM_ALL'},$format{'mergeNumber'}) if param('sel_cumPepAll') && param('sel_pepType') eq 'ALL';
				$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$numPeptides{$protID}{'CUM_ALL_TYPIC'},$format{'mergeNumber'}) if param('sel_cumPepAll') && param('sel_pepType') eq 'TYPIC';
				$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$numPeptides{$protID}{'CUM_NR'},$format{'mergeNumber'}) if  param('sel_cumPepNr');
				if (param('sel_cov')) {
					if ($pepCoverage{$protID} && $pepCoverage{$protID}{'BA'}) {$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$pepCoverage{$protID}{'BA'},$format{'mergeNumber1d'});}
					else {$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,'N/A',$format{'mergeText'});}
				}
				if (param('sel_cumCov')) {
					if ($pepCoverage{$protID} && $pepCoverage{$protID}{'CUM'}) {$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$pepCoverage{$protID}{'CUM'},$format{'mergeNumber1d'});}
					else {$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,'N/A',$format{'mergeText'});}
				}
				$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$protSpecificity{$protID}{'BA'},$format{'mergeNumber'}) if  param('sel_spec');
				$worksheet->merge_range($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$protSpecificity{$protID}{'BI'},$format{'mergeNumber'}) if  param('sel_bestSpec');
			}
			else {
				$worksheet->write_string($xlsRow,++$xlsCol,$proteinMGrHead{$protID},$format{'number'}) if param('sel_matchGr');
				$worksheet->write_number($xlsRow,++$xlsCol,$timesFound{$protID},$format{'number'}) if param('sel_times');
				$worksheet->write_number($xlsRow,++$xlsCol,$xlsEndRow,$xlsCol,$protPSM{$protID},$format{'Number'}) if param('sel_pepPSM');
				$worksheet->write_number($xlsRow,++$xlsCol,$bestScore{$protID},$format{'number'}) if param('sel_score');
				$worksheet->write_number($xlsRow,++$xlsCol,$numPeptides{$protID}{'BA_ALL'},$format{'number'}) if param('sel_pepAll') && param('sel_pepType') eq 'ALL';
				$worksheet->write_number($xlsRow,++$xlsCol,$numPeptides{$protID}{'BA_ALL_TYPIC'},$format{'number'}) if param('sel_pepAll') && param('sel_pepType') eq 'TYPIC';
				$worksheet->write_number($xlsRow,++$xlsCol,$numPeptides{$protID}{'BA_NR'},$format{'number'}) if param('sel_pepNr');
				$worksheet->write_number($xlsRow,++$xlsCol,$numPeptides{$protID}{'CUM_ALL'},$format{'number'}) if param('sel_cumPepAll') && param('sel_pepType') eq 'ALL';
				$worksheet->write_number($xlsRow,++$xlsCol,$numPeptides{$protID}{'CUM_ALL_TYPIC'},$format{'number'}) if param('sel_cumPepAll') && param('sel_pepType') eq 'TYPIC';
				$worksheet->write_number($xlsRow,++$xlsCol,$numPeptides{$protID}{'CUM_NR'},$format{'number'}) if  param('sel_cumPepNr');
				if (param('sel_cov')) {
					if ($pepCoverage{$protID} && $pepCoverage{$protID}{'BA'}) {$worksheet->write_number($xlsRow,++$xlsCol,$pepCoverage{$protID}{'BA'},$format{'number1d'});}
					else {$worksheet->write_string($xlsRow,++$xlsCol,'N/A',$format{'text'});}
				}
				if (param('sel_cumCov')) {
					if ($pepCoverage{$protID} && $pepCoverage{$protID}{'CUM'}) {$worksheet->write_number($xlsRow,++$xlsCol,$pepCoverage{$protID}{'CUM'},$format{'number1d'});}
					else {$worksheet->write_string($xlsRow,++$xlsCol,'N/A',$format{'text'});}
				}
				$worksheet->write_number($xlsRow,++$xlsCol,$protSpecificity{$protID}{'BA'},$format{'number'}) if  param('sel_spec');
				$worksheet->write_number($xlsRow,++$xlsCol,$protSpecificity{$protID}{'BI'},$format{'number'}) if  param('sel_bestSpec');
			}
		}
		else { # HTML
			print HTML "<TD $rowSpanStrg bgcolor=\"$rowColor\">$proteinMGrHead{$protID}</TD>" if param('sel_matchGr');
			print HTML "<TD $rowSpanStrg align=center bgcolor=\"$rowColor\">$timesFound{$protID}</TD>" if param('sel_times');
			print HTML "<TD $rowSpanStrg align=center bgcolor=\"$rowColor\">$bestScore{$protID}</TD>" if param('sel_score');
			print HTML "<TD $rowSpanStrg align=center bgcolor=\"$rowColor\">$protPSM{$protID}</TD>" if param('sel_pepPSM');
			print HTML "<TD $rowSpanStrg align=center bgcolor=\"$rowColor\">$numPeptides{$protID}{BA_ALL}</TD>" if param('sel_pepAll') && param('sel_pepType') eq 'ALL';
			print HTML "<TD $rowSpanStrg align=center bgcolor=\"$rowColor\">$numPeptides{$protID}{BA_ALL_TYPIC}</TD>" if param('sel_pepAll') && param('sel_pepType') eq 'TYPIC';
			print HTML "<TD $rowSpanStrg align=center bgcolor=\"$rowColor\">$numPeptides{$protID}{BA_NR}</TD>" if param('sel_pepNr');
			print HTML "<TD $rowSpanStrg align=center bgcolor=\"$rowColor\">$numPeptides{$protID}{CUM_ALL}</TD>" if param('sel_cumPepAll') && param('sel_pepType') eq 'ALL';
			print HTML "<TD $rowSpanStrg align=center bgcolor=\"$rowColor\">$numPeptides{$protID}{CUM_ALL_TYPIC}</TD>" if param('sel_cumPepAll') && param('sel_pepType') eq 'TYPIC';
			print HTML "<TD $rowSpanStrg align=center bgcolor=\"$rowColor\">$numPeptides{$protID}{CUM_NR}</TD>" if param('sel_cumPepNr');
			if (param('sel_cov')) {
				if ($pepCoverage{$protID} && $pepCoverage{$protID}{'BA'}) {printf HTML "<TD $rowSpanStrg align=center bgcolor=\"$rowColor\">%.1f</TD>",$pepCoverage{$protID}{'BA'};}
				else {print HTML "<TD $rowSpanStrg align=center bgcolor=\"$rowColor\">N/A</TD>";}
			}
			if (param('sel_cumCov')) {
				if ($pepCoverage{$protID} && $pepCoverage{$protID}{'CUM'}) {printf HTML "<TD $rowSpanStrg align=center bgcolor=\"$rowColor\">%.1f</TD>",$pepCoverage{$protID}{'CUM'};}
				else {print HTML "<TD $rowSpanStrg align=center bgcolor=\"$rowColor\">N/A</TD>";}
			}
			print HTML "<TD $rowSpanStrg align=center bgcolor=\"$rowColor\">$protSpecificity{$protID}{BA}</TD>" if param('sel_spec');
			print HTML "<TD $rowSpanStrg align=center bgcolor=\"$rowColor\">$protSpecificity{$protID}{BI}</TD>" if param('sel_bestSpec');
			print HTML "\n";
		}

		###<Listing peptides>###
		if ($exportFormat eq 'XLS') {
			if (scalar @selPepColName) {
				my $pepCount=0;
				my $pepRow=$xlsRow;
				foreach my $pepID (sort{$firstPepStart{$protID}{$a} <=> $firstPepStart{$protID}{$b} || length($pepSequences{$protID}{$a}) <=> length($pepSequences{$protID}{$b})} keys %{$pepSequences{$protID}}) {
					my $pepCol=$xlsCol;
					$pepCount++;
					$worksheet->write_number($pepRow,++$pepCol,$pepCount,$format{'number'});
					$worksheet->write_string($pepRow,++$pepCol,$pepSequences{$protID}{$pepID},$format{'seq'});
					$worksheet->write_string($pepRow,++$pepCol,$pepVarMods{$protID}{$pepID},$format{'text'}) if (param('sel_modPep'));
					$worksheet->write_string($pepRow,++$pepCol,$pepVarModsProt{$protID}{$pepID},$format{'text'}) if (param('sel_modPepProt'));
					$worksheet->write_string($pepRow,++$pepCol,$pepVarModsComment{$protID}{$pepID},$format{'text'}) if (param('sel_modPepProt') || param('sel_modPep'));
					$worksheet->write_string($pepRow,++$pepCol,$pepElutionTime{$pepID},$format{'text'}) if (param('sel_etPep'));
					$worksheet->write_string($pepRow,++$pepCol,$pepPhosphoRS{$pepID},$format{'text'}) if param('sel_phosphoRS');
					$worksheet->write_string($pepRow,++$pepCol,$probPhosphoRS{$pepID},$format{'text'}) if param('sel_phosphoRS');
					if (param('sel_posPep')) {
						my $posString='';
						foreach my $beg (sort{$a<=>$b} keys %{$pepBoundaries{$protID}{$pepID}}) {
							$posString.=',' if $posString;
							$posString.="$beg..$pepBoundaries{$protID}{$pepID}{$beg}";
						}
						$worksheet->write_string($pepRow,++$pepCol,$posString,$format{'text'});
					}
					$worksheet->write_number($pepRow,++$pepCol,scalar keys %{$pepBoundaries{$protID}{$pepID}},$format{'number'}) if param('sel_occPep');
					$worksheet->write_string($pepRow,++$pepCol,$isSpecific{$protID}{$pepID},$format{'text'}) if param('sel_isSpec');
					$worksheet->write_number($pepRow,++$pepCol,$pepScores{$protID}{$pepID},$format{'number'}) if param('sel_scorePep');
					$worksheet->write_number($pepRow,++$pepCol,$pepQValue{$protID}{$pepID},$format{'text'}) if param('sel_qvalPep');
					$worksheet->write_number($pepRow,++$pepCol,$pepCharge{$pepID},$format{'number'}) if param('sel_chargePep');
					$worksheet->write_number($pepRow,++$pepCol,$pepMass{$pepID},$format{'number3d'}) if param('sel_massPep');
					$worksheet->write_string($pepRow,++$pepCol,$pepComments{$pepID},$format{'text'}) if param('sel_commentPep');
					$worksheet->write_string($pepRow,++$pepCol,$pepTitle{$pepID}{title},$format{'text'}) if param('sel_titlePep');
					$pepRow++;
				}
			}
			$xlsRow+=$numListedPep;
		}
		else { # HTML
			if (scalar @selPepColName) {
				#my $firstPep=1;
				my $pepCount=0;
				foreach my $pepID (sort{$firstPepStart{$protID}{$a} <=> $firstPepStart{$protID}{$b} || length($pepSequences{$protID}{$a}) <=> length($pepSequences{$protID}{$b})} keys %{$pepSequences{$protID}}) {
					$pepCount++;
					#print HTML "<TR>" unless $firstPep;
					print HTML "<TR>" if $pepCount>1;
					#$firstPep=0;
					print HTML "<TD align=right bgcolor=\"$rowColor\">$pepCount</TD><TD bgcolor=\"$rowColor\">$pepSequences{$protID}{$pepID}</TD>";
					print HTML "<TD bgcolor=\"$rowColor\">$pepVarMods{$protID}{$pepID}</TD>" if (param('sel_modPep')); 
					print HTML "<TD bgcolor=\"$rowColor\">$pepVarModsProt{$protID}{$pepID}</TD>" if (param('sel_modPepProt'));
					print HTML "<TD bgcolor=\"$rowColor\">$pepVarModsComment{$protID}{$pepID}</TD>" if (param('sel_modPepProt') || param('sel_modPep'));
					print HTML "<TD bgcolor=\"$rowColor\">$pepElutionTime{$pepID}</TD>" if (param('sel_etPep')); 
					print HTML "<TD bgcolor=\"$rowColor\">$pepPhosphoRS{$pepID}</TD>" if param('sel_phosphoRS');
					print HTML "<TD bgcolor=\"$rowColor\">$probPhosphoRS{$pepID}</TD>" if param('sel_phosphoRS');

					if (param('sel_posPep')) {
						my $posString='';
						foreach my $beg (sort{$a<=>$b} keys %{$pepBoundaries{$protID}{$pepID}}) {
							$posString.=',' if $posString;
							$posString.="$beg..$pepBoundaries{$protID}{$pepID}{$beg}";
						}
						print HTML "<TD bgcolor=\"$rowColor\">$posString</TD>";
					}
					print HTML "<TD align=center bgcolor=\"$rowColor\">",scalar keys %{$pepBoundaries{$protID}{$pepID}},"</TD>" if param('sel_occPep');
					print HTML "<TD align=center bgcolor=\"$rowColor\">$isSpecific{$protID}{$pepID}</TD>" if param('sel_isSpec');
					print HTML "<TD align=center bgcolor=\"$rowColor\">$pepScores{$protID}{$pepID}</TD>" if param('sel_scorePep');
					print HTML "<TD align=center bgcolor=\"$rowColor\">$pepQValue{$protID}{$pepID}</TD>" if param('sel_qvalPep');
					print HTML "<TD align=center bgcolor=\"$rowColor\">$pepCharge{$pepID}</TD>" if param('sel_chargePep');
					printf HTML "<TD align=center bgcolor=\"$rowColor\">%.3f</TD>",$pepMass{$pepID} if param('sel_massPep');
					print HTML "<TD bgcolor=\"$rowColor\">$pepComments{$pepID}</TD>" if param('sel_commentPep');
					print HTML "<TD bgcolor=\"$rowColor\">$pepTitle{$pepID}{title}</TD>" if param('sel_titlePep');
					print HTML "</TR>\n";
				}
			}
			else {print HTML "\n</TR>\n";}
		}
		$xlsRow++ unless $numListedPep;
		# set row color
		$rowColor=($rowColor eq $rowColor1)? $rowColor2: $rowColor1;
		&printWaitingMsg('.','+=') unless $protCount % 1000;	
	}
	$xlsRow-- if $exportFormat eq 'XLS';
	#print HTML "<TR><TD colspan=$colSpan>&nbsp</TD></TR>\n";
}

############## Protein sort subroutines #########
sub sortRule {
	my ($sortItem,$sortPep)=@_;
	#<Identifier asc.
	if ($sortItem eq 'identifier') {return lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	#<Alias asc. (alias > identifier)
	elsif ($sortItem eq 'alias') {return lc($listProteins{$a}[1]) cmp lc($listProteins{$b}[1]) || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	#<Gene symbol asc. (gene symbol > num pep > identifier)
	elsif ($sortItem eq 'geneName') {return lc($proteinGeneSymbol{$a}) cmp lc($proteinGeneSymbol{$b}) || $numPeptides{$b}{$sortPep}<=>$numPeptides{$a}{$sortPep} || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	#<MW desc. (mw > num pep > identifier)
	elsif ($sortItem eq 'mw') {return $listProteins{$b}[$mwColPos]<=>$listProteins{$a}[$mwColPos] || $numPeptides{$b}{$sortPep}<=>$numPeptides{$a}{$sortPep} || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	#<Peptide desc. (num pep > identifier)
	elsif ($sortItem=~/pep/i) {return $numPeptides{$b}{$sortPep}<=>$numPeptides{$a}{$sortPep} || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	#<Score desc. (score > identifier)
	elsif ($sortItem eq 'bestScore') {return $bestScore{$b}<=>$bestScore{$a} || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	#<Coverage desc. (coverage > identifier)
	elsif ($sortItem eq 'bestCov') {return $pepCoverage{$b}{'BA'}<=>$pepCoverage{$a}{'BA'} || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	#<Cumulated coverage desc. (coverage > identifier)
	elsif ($sortItem eq 'cumCov') {return $pepCoverage{$b}{'CUM'}<=>$pepCoverage{$a}{'CUM'} || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	#<Times found desc. (times found > num pep > identifier)
	elsif ($sortItem eq 'times_found') {return $timesFound{$b}<=>$timesFound{$a} || $numPeptides{$b}{$sortPep}<=>$numPeptides{$a}{$sortPep} || lc($listProteins{$a}[0]) cmp lc($listProteins{$b}[0])}
	else {return 0}
}

sub getPhosphoRsString{
	my ($prsDBString) = @_;

	my ($status, $score, $positions) = split /;/, $prsDBString;

	my $prsString;
	if ($status == 3 || $status == 0) {
		$prsString = "Confirmed ($score%)";
	} elsif ($status == 2){
		$prsString = "Changed ($score%)";
	} elsif ($status == 1) {
		$prsString = "Insufficient probability";
	} else {
		$prsString = "No phosphorylation evidence (or manually validated site)";
	}

	return $prsString;
}

sub getPhosphoRsProbString {
	my ($prsOutputFile,$refPepPhosphoRS,$refDataPhosphoRS,$refProbPhosphoRS) = @_;
	return unless -e $prsOutputFile;
	
	my ($spectrumOK,$sitePredictionOK,$currentPepID) = (0,0,0);
	my @positions=();
	open (XML,$prsOutputFile);
    while(<XML>){
		if (/<Spectrum ID="(\d+)">/ && $refDataPhosphoRS->{$1}) {
			$currentPepID=$refDataPhosphoRS->{$1};
			$spectrumOK = 1;
		}
		elsif ($spectrumOK) {
			if (/<SitePrediction>/) {$sitePredictionOK = 1;}
			elsif (/<\/SitePrediction>/) {
				$refProbPhosphoRS->{$currentPepID}=join(', ',@positions);
				($spectrumOK,$sitePredictionOK,$currentPepID) = (0,0,0);
				@positions=();
			}
			elsif ($sitePredictionOK) {
				my ($info1,$pos,$info2,$prob,@infos)=split(/\"/,$_);
				$prob = sprintf("%.2f", 100*$prob) . "%";
				push @positions,"$pos($prob)";
			}
		}
    }
    close XML;
}

sub printWaitingMsg {
	my ($msg, $append) = @_;
	my $appendHTML = ($append) ? "+=" : "=";
	print qq
|<SCRIPT type="text/javascript">document.getElementById('waitSPAN').innerHTML $appendHTML '$msg';</SCRIPT>
|;
}

####>Revision history<####
# 2.6.0 [CHANGE] Make the generation of the experiment design file optional (VS 04/02/21)
# 2.5.9 [ENHANCEMENT] Added export of skyline data file format (VS 02/02/21)
# 2.5.8 [BUGFIX] Fix undef depth when exporting a list of proteins (VS 22/07/20)
# 2.5.7 [BUGFIX] Fix analysis file name displaying (VS 08/06/20)
# 2.5.6 [BUGFIX] Fix path computation for search files export: old analysis (e.g. 2012) with ended validation status did not have their files stored at the same place (VS 17/04/2020)
# 2.5.5 [ENHANCEMENT] Add possibility to export a specific list of the selected theme (VS 09/04/20)
# 2.5.4 [ENHANCEMENT] Handles DIA/Spectronaut search files export  (VS 06/04/20)
# 2.5.4a [BUGFIX] in mapping resource check for Project (PP 03/04/20)
# 2.5.3 [UPDATE] Changed RANK field to IDENT_RANK for compatibility with MySQL 8 (PP 04/03/20) 
# 2.5.2 [ENHANCEMENT] Optimized PhosphoRS file parsing for positions probablility (PP 29/01/20)
# 2.5.1 [CHANGES] Add more robustness to search files export (VS 08/01/20)
# 2.5.0 [ENHANCEMENT] File(s) stored server-side to prevents server timeout & [FEATURE] List exclusion/restriction for project hierarchy (PP 07/11/19) 
# 2.4.6 [ENHANCEMENT] Handles search files runned multiple times : searchFile-(X).msf/dat (VS 30/10/19)
# 2.4.5 [BUGFIX] Fix amount of peptides in best analysis and cumulative coverage value computing (VS 28/10/19)
# 2.4.4 [FEATURE] Remove locked experiments from exportable proteins queries (VS 08/08/19)
# 2.4.3 Fix issue with search files export having a different name than sample/analysis, depending on context (VS 13/03/19)
# 2.4.2 Added possibility to export whole analyses files: .msf (VS 26/02/19)
# 2.4.1 Added new export options: variable modifications based on protein sequence + retention time (VS 14/12/18)
# 2.4.0 Faster form display by moving mapped identifiers detection to user-dependent ajax call (PP 11/04/18)
# 2.3.8 Export PhosphoRS probabilities (GA 31/01/18)
# 2.3.7 fix little bug with index L1029 (SL 17/06/16)
# 2.3.6 Speed optimization (PP 20/01/16)
# 2.3.5 Add charge and Q-value in peptide data checkbox option so as to fit with http://www.microvesicles.org/ submission format (GA 04/01/16)
# 2.3.4 Minor changes in displayed text (PP 10/06/15)
# 2.3.3 Add proteotypic peptides extraction and PSM (GA 19/01/15)
# 2.3.2 Fix bug in relevant PTM management (PP ../06/14)
# 2.3.1 Support for multi-databank Analysis & better check on external links & UTF-8 issues (PP 18/04/14)
# 2.3.0 Minor bug fix (PP 22/01/13)
# 2.2.9 Remove VAR_MOD from script (GA 22/05/13)
# 2.2.8 Changed Protein varMod to Relevant PTMs (PP 22/04/13)
# 2.2.7 Add variable-modification option for proteins (FY 11/04/13)
# 2.2.6 Add phosphoRS status support and minor fixes (FY 10/04/13)
# 2.2.5 Converts varMod position ambiguity from -1 to ? (PP 21/03/13)
# 2.2.4 New GPL license (PP 12/03/13)
# 2.2.3 True Excel export & compatible with new identifier mapping strategy (PP 19/02/13)
# 2.2.2 Add subItem comments if any (PP 06/04/12)
# 2.2.1 fix bug unwanted Specificity column + html as true attachment (PP 13/03/12)
# 2.2.0 adapting to new protein list management options (PP 23/02/2011)
