#!/usr/local/bin/perl -w

################################################################################
# autoSelect.cgi              1.6.7                                            #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Displays myProMS main entry page with links to different sections            #
# Called after login to server                                                 #
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
use promsConfig;
use promsMod;

#print header(-'content-encoding'=>'no'),"\n"; warningsToBrowser(1); # DEBUG
##########################################
####>Configuration & Connection to DB<####
##########################################
my %promsPath=&promsConfig::getServerInfo;
#my %massSpecType=&promsConfig::getMsType;
my %msTypeName=&promsConfig::getMsType;
my $absMaxRank=&promsConfig::getMaxRank;
my $dbh=&promsConfig::dbConnect;

####################
####>Parameters<####
####################
my $userID = $ENV{'REMOTE_USER'}; #for templates
my $ITEM=(param('ITEM'))? uc(param('ITEM')) : ''; # not defined for &clearValidations for !ANALYSIS
my $itemID=param('ID');
my $projectID=param('id_project');
my $action=(param('ACT'))? param('ACT') : (param('clearMode'))? param('clearMode') : 'select';
#my $selMsType=(param('MSTYPE'))? param('MSTYPE') : 'MIS'; # only for auto-select within validation mode (1 analysis)
my $selMsType=param('MSTYPE'); # not defined for 1st call of autoselect outside validation mode
my @analysisList=($ITEM eq 'ANALYSIS')? ($itemID) : (param('anaList'))? param('anaList') : (); #Global!!!!!
if (scalar @analysisList) { # Only 1 ana or after form submission (except auto-select)
	if ($action eq 'restore') {&restoreValidation;} # no longer maintained
	elsif ($action eq 'clearSelection' && $ITEM eq 'ANALYSIS') {&selectClearLevel;} # only for in validation mode (1 analysis)
	elsif ($action eq 'clearAuto') {&clearValidations(1);} # clears auto validation only
	elsif ($action eq 'clearAll') {&clearValidations(2);} # clears manual + auto validation
	elsif ($action eq 'lowScores') {&activateLowScores;}
}

####>$action=select from here on ----------------->
my ($selFileFormat,$maxRank,$absMinScore,$disableMS1,$disableMS2);
if (param('SELFORMAT')) {
	$selFileFormat=param('SELFORMAT');
}
elsif ($ITEM eq 'ANALYSIS') {
	(my $anaFormat,$maxRank,$absMinScore)=$dbh->selectrow_array("SELECT FILE_FORMAT,MAX_RANK,MIN_SCORE FROM ANALYSIS WHERE ID_ANALYSIS=$itemID");
	($selFileFormat)=($anaFormat=~/(\w+)\./);
}
$selFileFormat='SEQUEST' unless $selFileFormat;

if (param('submit')) { # select or clearSelection
	&processSelectionForm;
	exit;
}

my ($itemName)=$dbh->selectrow_array("SELECT NAME FROM $ITEM WHERE ID_$ITEM=$itemID");
my $targetString='';
my $titleString='<FONT class="title">';
if ($ITEM eq 'ANALYSIS') {
	$titleString.="<SELECT class=\"title\" style=\"text-align:left\" onchange=\"goTocomparative();\"><OPTION selected>Qualitative</OPTION><OPTION>Comparative</OPTION></SELECT>";
}
else {$titleString.="Qualitative";}
$titleString.=" Peptide/Protein Selection</FONT><BR>\n";
#my $subTitleString='';
#my $fileFormat;
my $checkChkBox='';
my (%listDataBank,@itemAnalyses);
my (%lowerScPepAnalyses,%existMsType,%sthLS);
my $fileFormatString='';
my $lowScoreString='';
my $selectMsTypeString='';

if ($ITEM eq 'ANALYSIS') {
	foreach my $rank (1..$absMaxRank) {
		$sthLS{$rank}=$dbh->prepare("SELECT 1 FROM QUERY_VALIDATION WHERE INFO_PEP$rank LIKE '%SEL=-1%' AND ID_ANALYSIS=? LIMIT 0,1");
	}
	$targetString='target="spectrumFrame"';
	#my ($fileFormat)=$dbh->selectrow_array("SELECT FILE_FORMAT FROM ANALYSIS WHERE ID_ANALYSIS=$itemID");
	#$absMinScore=&promsConfig::getMinScore($fileFormat);

	#>Checking for inactivated lower score peptides (SEL=-1)
	$lowScoreString='<BR>';
	foreach my $rank (sort{$a<=>$b} keys %sthLS) {
		$sthLS{$rank}->execute($itemID);
		my ($inactiveLS)=$sthLS{$rank}->fetchrow_array;
		if ($inactiveLS) {
			$lowScoreString.='<FONT class="title3">Make <U>all</U> lower scoring interpretations selectable <INPUT type="button" value="Proceed" onclick="activateLowScores()"></FONT>';
			last;
		}
	}
	$existMsType{$selMsType}=1; # selMsType passed as param
	$selectMsTypeString="<INPUT type=\"hidden\" name=\"MSTYPE\" value=\"$selMsType\">";
}
else {
	$titleString="<BR>$titleString<BR>";
	#@analysesMIS=&getAnalysesList('MIS');
	# if (scalar @analysesMIS) {
		# $subTitleString.="<BR><FONT class=\"title3\">(Apply to all <SELECT class='title3' onchange='changeMsType(this.value);'><OPTION value='MIS'";
		# $subTitleString.=' selected' if $msType eq 'MIS';
		# $subTitleString.=">$massSpecType{'MIS'}</OPTION><OPTION value='PMF'";
		# $subTitleString.=' selected' if $msType eq 'PMF';
		# $subTitleString.=">$massSpecType{'PMF'}</OPTION></SELECT> Analyses in ";
		# $subTitleString.=&promsMod::getItemType($ITEM)." <FONT color=\"#DD0000\">$itemName</FONT>)</FONT><BR>\n";
	# }
	
	my ($selMascot,$selPhenyx,$selSequest,$selParagon,$selTandem);
	if($selFileFormat eq 'MASCOT'){
		($selMascot,$selPhenyx,$selSequest,$selParagon,$selTandem)=(' selected','','','','');
	}elsif($selFileFormat eq 'PHENYX'){
		($selMascot,$selPhenyx,$selSequest,$selParagon,$selTandem)=('',' selected','','','');
	}elsif($selFileFormat eq 'SEQUEST'){
		($selMascot,$selPhenyx,$selSequest,$selParagon,$selTandem)=('','',' selected','','');
	}elsif($selFileFormat eq 'PARAGON'){
		($selMascot,$selPhenyx,$selSequest,$selParagon,$selTandem)=('','','',' selected','');
	}elsif($selFileFormat eq 'TDM'){
		($selMascot,$selPhenyx,$selSequest,$selParagon,$selTandem)=('','','','',' selected');
	}
	$fileFormatString="<FONT class=\"title2\">&nbsp;Search engine:</FONT><SELECT name=\"SELFORMAT\" class=\"title3\" onchange=\"changeSearch()\"><OPTION value=\"MASCOT\"$selMascot>Mascot</OPTION><OPTION value=\"PHENYX\"$selPhenyx>Phenyx</OPTION><OPTION value=\"SEQUEST\"$selSequest>Sequest</OPTION><OPTION value=\"PARAGON\"$selParagon>Paragon</OPTION><OPTION value=\"TDM\"$selTandem>Tandem</OPTION></SELECT>";

	my @sthItem;
	my $baseFieldString='MAX_RANK,MIN_SCORE,ID_ANALYSIS,VALID_STATUS,MS_TYPE,DATA_FILE,FILE_FORMAT,WIFF_FILE,TAXONOMY,0'; # 0 will be replaced by list of DB used
	if ($ITEM eq 'SAMPLE') {
		$sthItem[0]=$dbh->prepare("SELECT $baseFieldString,NAME FROM ANALYSIS WHERE ID_SAMPLE=$itemID ORDER BY DISPLAY_POS ASC");
	}
	elsif ($ITEM eq 'SPOT') {
		$sthItem[0]=$dbh->prepare("SELECT $baseFieldString,ANALYSIS.NAME FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_SPOT=$itemID ORDER BY ANALYSIS.DISPLAY_POS ASC");
	}
	elsif ($ITEM eq 'GEL2D') {
		$sthItem[0]=$dbh->prepare("SELECT $baseFieldString,'SPOT',SPOT.NAME,ANALYSIS.NAME FROM ANALYSIS,SAMPLE,SPOT WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT=SPOT.ID_SPOT AND ID_GEL2D=$itemID ORDER BY SPOT.NAME ASC, ANALYSIS.DISPLAY_POS ASC");
	}
	elsif ($ITEM eq 'EXPERIMENT') {
		$sthItem[0]=$dbh->prepare("SELECT $baseFieldString,'GEL2D',GEL2D.NAME,'SPOT',SPOT.NAME,ANALYSIS.NAME FROM ANALYSIS,SAMPLE,SPOT,GEL2D WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT=SPOT.ID_SPOT AND SPOT.ID_GEL2D=GEL2D.ID_GEL2D AND GEL2D.ID_EXPERIMENT=$itemID ORDER BY GEL2D.DISPLAY_POS ASC, SPOT.NAME ASC, ANALYSIS.DISPLAY_POS ASC");
		$sthItem[1]=$dbh->prepare("SELECT $baseFieldString,'SAMPLE',SAMPLE.NAME,ANALYSIS.NAME FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT IS NULL AND ID_EXPERIMENT=$itemID ORDER BY SAMPLE.DISPLAY_POS ASC, ANALYSIS.DISPLAY_POS ASC");
	}
	my $sthAD = $dbh->prepare("SELECT D.ID_DATABANK,NAME,VERSION_NAME FROM ANALYSIS_DATABANK AD,DATABANK D WHERE AD.ID_DATABANK=D.ID_DATABANK AND AD.ID_ANALYSIS=?");
	$maxRank=0;
	$absMinScore=100;
	foreach my $sth (@sthItem) {
		$sth->execute;
		while (my ($maxRk,$minSc,@anaData)=$sth->fetchrow_array) {
			$sthAD->execute($anaData[0]);
			my @dbUsed;
			while (my ($dbID,$dbName,$version)=$sthAD->fetchrow_array) {
				push @dbUsed,$dbID;
				$listDataBank{$dbID}=$dbName;
				$listDataBank{$dbID}.=" ($version)" if $version;
			}
			$anaData[7]=\@dbUsed;
			push @itemAnalyses,\@anaData;
			$maxRank=$maxRk if ($maxRk && $maxRank<$maxRk);
			$absMinScore=$minSc if (defined $minSc && $absMinScore>$minSc);
			$existMsType{$anaData[2]}++; # MS type
		}
		$sth->finish;
	}
	$sthAD->finish;
	$checkChkBox=qq
|	if (testCheckbox(myForm)<1) {
		alert('No Analyses selected.');
		return false;
	}
|;
	##>Guessing selected MS Type
	unless ($selMsType) {
		$selMsType=($existMsType{'SQ'} || scalar keys %existMsType > 1)? 'SQ' : ($existMsType{'PMF'})? 'PMF' : 'MIS';
	}
	$selectMsTypeString="&nbsp;<FONT class=\"title2\">type:</FONT><SELECT name=\"MSTYPE\" onchange=\"changeSearch()\" class=\"title3\" style=\"width:220px\">";
	foreach my $msType (sort{$a cmp $b} keys %existMsType) {
		$selectMsTypeString.="<OPTION value=\"$msType\"";
		$selectMsTypeString.=' selected' if $msType eq $selMsType;
		$selectMsTypeString.=">$msTypeName{$msType}</OPTION>";
	}
	$selectMsTypeString.="</SELECT>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;";
}
foreach my $rank (keys %sthLS) {$sthLS{$rank}->finish;}

if ($selMsType eq 'MIS') {
	$maxRank=$absMaxRank unless $maxRank;
	$absMinScore=&promsConfig::getMinScore($selFileFormat) unless defined($absMinScore);
	$disableMS1='disabled';
	$disableMS2='';
}
elsif ($selMsType eq 'SQ') {
	$maxRank=$absMaxRank unless $maxRank;
	$absMinScore=&promsConfig::getMinScore($selFileFormat) unless defined($absMinScore);
	$disableMS1='';
	$disableMS2='';
}
else { # PMF
	$maxRank=$absMaxRank; # not definable during import
	$absMinScore=0;
	$disableMS1='';
	$disableMS2='disabled';
}

###>Template management<###
my $selectTemplateString='';
my $defaultIndex = 0;
#if ($action eq 'select') {
	my $sthGetTemplate = $dbh->prepare("SELECT NAME, PARAM_STRG, DES, IS_DEFAULT FROM VALIDATION_TEMPLATE WHERE ID_USER='$userID' AND MS_TYPE='$selMsType' AND SEARCH_ENGINE='$selFileFormat' AND VAL_TYPE='quali'");
	$sthGetTemplate->execute;
	my $refTemplateList = $sthGetTemplate->fetchall_arrayref;
	$sthGetTemplate->finish;
	$selectTemplateString = "<FONT class=\"title2\">Template:</FONT><SELECT name=\"templateList\" onchange=\"useTemplate(this.selectedIndex)\" class=\"title3\" style=\"width:250px\">\n<OPTION value=''>None</OPTION>\n";
	my $i=0;
	foreach my $template (@{$refTemplateList}){
		$i++;
		my $templateName = $template->[0];
		my $parameters = $template->[1];
		my $des = ($template->[2])? "<B>Description : </B>".&promsMod::HTMLcompatible($template->[2]) : "No description" ;
		$des =~ s/\/\//\\\/\\\//g; $des =~ s/&apos;*/\\\'/g; ### javascript compatible
		$selectTemplateString .= "<OPTION value=\"$parameters\" onmouseover=\"popup('$des')\" onmouseout=\"popout()\">$templateName</OPTION>\n";
		if($template->[3]){
			$defaultIndex = $i;
		}
	}
	$selectTemplateString .= "</SELECT>";
#}

$dbh->disconnect;


#######################
####>Starting HTML<####
#######################
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
print header(-'content-encoding'=>'no',-charset=>'utf-8');
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Auto Validation Form</TITLE>
<LINK rel=stylesheet href="$promsPath{html}/promsStyle.css" type=text/css>
<STYLE type="text/css">
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
|;
if ($ITEM eq 'ANALYSIS') {
	print qq
|function activateLowScores() {
	window.location="$promsPath{cgi}/autoSelect.cgi?id_project=$projectID&ITEM=$ITEM&ID=$itemID&ACT=lowScores";
}
function goTocomparative() {
	window.location="$promsPath{cgi}/filterPeptideAnalysis.cgi?ITEM=$ITEM&ID=$itemID&id_project=$projectID&SELFORMAT=$selFileFormat&MSTYPE=$selMsType";
}
|;
}
else {
	#&promsMod::popupInfo();
	print qq
|function changeSearch() {
	var myForm=document.autoValForm;
	window.location="$promsPath{cgi}/autoSelect.cgi?id_project=$projectID&ITEM=$ITEM&ID=$itemID&ACT=select&SELFORMAT="+myForm.SELFORMAT.value+"&MSTYPE="+myForm.MSTYPE.value;
}
|;
}
print qq
|function updateScores() {
	var myForm=document.autoValForm;
	var scProt=parseFloat(myForm.minProtSc.value);
|;
if ($selFileFormat=~/SEQUEST/) {
	for(my $i=2 ; $i<=4; $i++){
		for(my $j=1 ; $j<=3 ; $j++){
			print "\tvar sc$i$j=parseFloat(myForm.minInt$i"."Sc$j.value);\n";
		}
	}
	for(my $i=2 ; $i<=4; $i++){
		for(my $j=1 ; $j<=3 ; $j++){
			print "\tif (!sc$i$j \|\| sc$i"."1 < $absMinScore) {sc$i$j=$absMinScore;}\n";
		}
	}
	for(my $i=2 ; $i<=4; $i++){
		for(my $j=2 ; $j<=3 ; $j++){
			my $previous=$j-1;
			print "\tif (sc$i$j > sc$i$previous) {sc$i$j=sc$i$previous;}\n";
		}
	}
	for(my $i=2 ; $i<=4; $i++){
		print "\tif (scProt < sc$i"."3) {scProt=sc$i"."3;}\n";
	}
	for(my $i=2 ; $i<=4; $i++){
		for(my $j=1 ; $j<=3 ; $j++){
			print "\tmyForm.minInt$i"."Sc$j.value=sc$i$j;\n";
		}
	}
}
else {
	print qq
|	var sc1=parseFloat(myForm.minIntSc1.value);
	var sc2=parseFloat(myForm.minIntSc2.value);
	var sc3=parseFloat(myForm.minIntSc3.value);
	if (!sc1 \|\| sc1 < $absMinScore) {sc1=$absMinScore;}
	if (!sc2 \|\| sc2 < $absMinScore) {sc2=$absMinScore;}
	if (!sc3 \|\| sc3 < $absMinScore) {sc3=$absMinScore;}
	if (sc2 > sc1) {sc2=sc1;}
	if (sc3 > sc2) {sc3=sc2;}
	if (scProt < sc3) {scProt=sc3;}
	myForm.minIntSc1.value=sc1;
	myForm.minIntSc2.value=sc2;
	myForm.minIntSc3.value=sc3;
	myForm.minProtSc.value=scProt;
|;
}
print qq
|}
function activatePepSel() {
	var myForm=document.autoValForm;
	var disab;
	if (myForm.selGoodInt.checked==false && myForm.rejBadInt.checked==false){disab=true;}
	else {disab=false;}
|;
if ($selMsType ne 'PMF') {
	if ($selFileFormat=~/SEQUEST/) {
		for(my $i=2 ; $i<=4; $i++){
			for(my $j=1 ; $j<=3 ; $j++){
				print "\tmyForm.minInt$i"."Sc$j.disabled=disab;\n";
			}
		}
	}
	else {
		print qq
|	myForm.minIntSc1.disabled=disab;
	myForm.minIntSc2.disabled=disab;
	myForm.minIntSc3.disabled=disab;
|;
	}
}
print qq
|	myForm.noReject.disabled=disab;
	myForm.minSize.disabled=disab;
	myForm.maxSize.disabled=disab;
	myForm.minDelta.disabled=disab;
	//myForm.newMaxRank.disabled=disab;
|;
print "\tmyForm.newMaxRankMS1.disabled=disab;\n" if $selMsType ne 'MIS';
print "\tmyForm.newMaxRankMS2.disabled=disab;\n" if $selMsType ne 'PMF';
print qq
|	myForm.flaggedUp.disabled=disab;
	myForm.flaggedDown.disabled=disab;
	myForm.oneInt.disabled=disab;
	myForm.overPep.disabled=disab;
	if (myForm.selGoodInt.checked==false) {
		myForm.oneInt.disabled=true;
	}
}
function activateProtSel() {
	var myForm=document.autoValForm;
	var disab;
	if (myForm.selProt.checked==false){disab=true;} else {disab=false;}
	myForm.minPep.disabled=disab;
	myForm.minCov.disabled=disab;
|;
print "\tmyForm.minProtSc.disabled=disab;\n" if $selMsType ne 'PMF';
print qq
|}
function showTemplateBox() {
	var myForm=document.autoValForm;
	myForm.template_name.style.visibility=(myForm.actTemplate.checked==false)?"hidden" : "visible";
}
function checkExistingTemplates(myForm){
	if(myForm.template_name.value){
	    for (var i=0;i<myForm.templateList.length;i++){
			if (myForm.template_name.value == myForm.templateList[i].value){
				return confirm("Are you sure to want to overwrite "+myForm.templateList[i].value+" ?");
				break;
			}
		}
    } else {
	    return true;
    }
}
function checkForm(myForm) {
    if (myForm.selGoodInt.checked==false && myForm.rejBadInt.checked==false && myForm.selProt.checked==false && myForm.bestMatch.checked==false){
		alert('Select a validation option.');
		return false;
	}
$checkChkBox
    if(myForm.actTemplate.checked){
	    if(myForm.template_name.value){ //checking template fields
			for (var i=0;i<myForm.templateList.length;i++){
		        if (myForm.template_name.value.toLowerCase() == myForm.templateList[i].text.toLowerCase()){
					return confirm("Are you sure to want to overwrite "+myForm.templateList[i].text+" ?");
		        }
			}
	    } else {
			alert("Enter a template name.");
			return false;
	    }
    }
    return true;
}
|;
if ($ITEM eq 'ANALYSIS') {
	print qq
|function cancelAction() {
	if (parent.selectedView=='myList') {parent.writeMenu(parent.selectedProtID,2);}
	else {parent.setView(parent.selectedView);}
}
top.promsFrame.spectrumFrame.location="$promsPath{html}/nothing.html";
|;
}
else {
	print qq
|function cancelAction() {
	//top.promsFrame.selectedAction='view'; // set default action to 'view'
	top.promsFrame.optionFrame.selectOption();

}
function testCheckbox(myForm) {
	var selected=0;
	var anaBox=myForm.anaList;
	if (anaBox.length) { // more than 1 checkboxes
		for (i=0; i < anaBox.length; i++){
			if (anaBox[i].checked==true) {selected=1; break;}
		}
	}
	else if (anaBox.checked==true){selected=1;}
	return selected;
}
function checkall(checkStatus){
	var anaBox=document.autoValForm.anaList;
	if (!anaBox) return; // no matching analyses
	if (anaBox.length) { // more than 1 checkboxes
		for (i=0; i < anaBox.length; i++){
			if (anaBox[i].disabled) continue;
			anaBox[i].checked=checkStatus;
		}
	}
	else {anaBox.checked=checkStatus;} // Only 1 checkboxes
}
|;
}
&promsMod::popupInfo;
print qq
|function useTemplate(index){
	var myForm=document.autoValForm;
	if (index==0) { // No template
		document.autoValForm.reset();
		return;
	}
	var paramString = myForm.templateList.options[index].value;
	if (!paramString){
		myForm.reset();
	}
	var paramList = paramString.split(';');
	for(var i=0; i<paramList.length; i++){
		var parameters = paramList[i].split(':');
		if (!myForm[parameters[0]]) {continue;} // just to be safe
		if (parameters[1] == 'check'){
			myForm[parameters[0]].checked = (parameters[2] == 1)? true : false;
		} else if (parameters[1] == 'text'){
			myForm[parameters[0]].value = parameters[2];
		} else if (parameters[1] == 'select'){
			myForm[parameters[0]].selectedIndex = parameters[2] - 1; //only for maxRank
		} else if (parameters[1] == 'radio'){
			for(var j=0;j<myForm[parameters[0]].length;j++){
				if(myForm[parameters[0]][j].value == parameters[2]){
					myForm[parameters[0]][j].checked = true;
				} else {
					myForm[parameters[0]][j].checked = false;
				}
			}
		}
	}
	myForm.template_name.value = (index > 0)? myForm.templateList.options[index].text : '';
}
</SCRIPT>
</HEAD>
<BODY topmargin=0 background="$promsPath{images}/bgProMS.gif" onload="document.autoValForm.templateList.selectedIndex = $defaultIndex; useTemplate($defaultIndex);">
<CENTER>
<FORM name="autoValForm" method="post" onsubmit="return(checkForm(this));" $targetString>
$titleString
$lowScoreString
<INPUT type="hidden" name="ITEM" value="$ITEM">
<INPUT type="hidden" name="ID" value="$itemID">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="ACT" value="$action">
<TABLE align=center border=0>
<TR><TD colspan=3>$fileFormatString$selectMsTypeString$selectTemplateString</TD></TR>
<TR><TD valign=top><TABLE cellspacing=0>
<TR bgcolor=$darkColor><TH colspan=2><FONT class="title2">&nbsp&nbsp;Interpretation Selection Rules</FONT></TH></TR>
<TR bgcolor=$lightColor><TD width=30></TD><TH nowrap align=left>
<INPUT type="checkbox" name="selGoodInt" value=1 onclick="activatePepSel()" checked > Select interpretations meeting the following criteria.<BR>
<INPUT type="checkbox" name="rejBadInt" value=1 onclick="activatePepSel()" checked > Reject interpretations that <U>do not</U> meet these criteria.<BR>
<FONT style="font-size:4px"><BR></FONT>+ Criteria:<BR>
|;
if ($selFileFormat=~/SEQUEST/) {
	print qq
|<TABLE cellspacing=0 cellpadding=0>
	<TR><TH align=left>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- minimum score<SUP>*</SUP> for:&nbsp;&nbsp;</TH>
		<TH align=left nowrap>2<SUP>+</SUP> peptides
1:<INPUT type="text" name="minInt2Sc1" id="minInt2Sc1" value="$absMinScore" size="4" onchange="updateScores()" $disableMS2>,
2:<INPUT type="text" name="minInt2Sc2" id="minInt2Sc2" value="$absMinScore" size="4" onchange="updateScores()" $disableMS2>,
<FONT class="title3">&#8805</FONT>3:<INPUT type="text" name="minInt2Sc3" id="minInt2Sc3" value="$absMinScore" size="4" onchange="updateScores()" $disableMS2> peptides/protein.<TH>
	</TR>
	<TR><TH></TH><TH align=left nowrap>3<SUP>+</SUP> peptides
1:<INPUT type="text" name="minInt3Sc1" id="minInt3Sc1" value="$absMinScore" size="4" onchange="updateScores()" $disableMS2>,
2:<INPUT type="text" name="minInt3Sc2" id="minInt3Sc2" value="$absMinScore" size="4" onchange="updateScores()" $disableMS2>,
<FONT class="title3">&#8805</FONT>3:<INPUT type="text" name="minInt3Sc3" id="minInt3Sc3" value="$absMinScore" size="4" onchange="updateScores()" $disableMS2> peptides/protein.<TH>
	</TR>
	<TR><TH></TH><TH align=left nowrap>4<SUP>+</SUP> peptides
1:<INPUT type="text" name="minInt4Sc1" id="minInt4Sc1" value="$absMinScore" size="4" onchange="updateScores()" $disableMS2>,
2:<INPUT type="text" name="minInt4Sc2" id="minInt4Sc2" value="$absMinScore" size="4" onchange="updateScores()" $disableMS2>,
<FONT class="title3">&#8805</FONT>3:<INPUT type="text" name="minInt4Sc3" id="minInt4Sc3" value="$absMinScore" size="4" onchange="updateScores()" $disableMS2> peptides/protein.<TH>
	<TR></TABLE>
|;
}
else {
	print qq
|&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- minimum score<SUP>*</SUP> if
1:<INPUT type="text" name="minIntSc1" value="$absMinScore" size="4" onchange="updateScores()" $disableMS2>,
2:<INPUT type="text" name="minIntSc2" value="$absMinScore" size="4" onchange="updateScores()" $disableMS2>,
<FONT class="title3">&#8805</FONT>3:<INPUT type="text" name="minIntSc3" value="$absMinScore" size="4" onchange="updateScores()" $disableMS2> peptides/protein.<BR>
|;
}
print qq
|&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;
<INPUT type="checkbox" name="noReject" value=1 checked /> Exclude already rejected interpretations (if any) from count.</FONT><BR>
&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- and size from <INPUT type="text" name="minSize" value="" size="1" placeholder="1">aa to <INPUT type="text" name="maxSize" value="" size="1" placeholder="999">aa.<BR>
&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- and mass error <FONT class="title3">&#8804</FONT><INPUT type="text" name="minDelta" value="" size="4" placeholder="10"> Da.<BR>
&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- and Rank for MS <FONT class=\"title3\">&#8804</FONT><SELECT name=\"newMaxRankMS1\" $disableMS1>
|;
foreach my $rank (1..$absMaxRank) {
	print "<OPTION value=\"$rank\"";
	print ' selected' if $rank==$absMaxRank;
	print ">$rank</OPTION>";
}
print "</SELECT> and/or for MS/MS<SUP>*</SUP> <FONT class=\"title3\">&#8804</FONT><SELECT name=\"newMaxRankMS2\" $disableMS2>\n";
foreach my $rank (1..$absMaxRank) {
	print "<OPTION value=\"$rank\"";
	print ' selected' if $rank==$maxRank;
	print ">$rank</OPTION>";
}
print qq
|</SELECT><BR>
&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;-<INPUT type="checkbox" name="flaggedUp" value="1"> Interpretation was <IMG src="$promsPath{images}/lightYellow1.gif" hspace=0 border=0 height=11 width=11> flagged.<BR>
&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;-<INPUT type="checkbox" name="flaggedDown" value="1"> Interpretation was not <IMG src="$promsPath{images}/lightOrange1.gif" hspace=0 border=0 height=11 width=11> flagged.<BR>


<FONT style="font-size:4px"><BR></FONT>
<INPUT type="checkbox" name="oneInt" value=1 checked /> Select only 1 (best) interpretation/query.<BR>
<INPUT type="checkbox" name="overPep" value=1 /> Overwrite previous selections/rejections.<BR>
<FONT style="font-size:4px"><BR></FONT>
</TH></TR>
<TR bgcolor=$lightColor><TH colspan=2>
<FONT style="font-style:italic;">All proteins matched by selected peptides will be selected unless<BR>Protein Selection Rules are applied.</FONT>
</TH></TR>

<TR bgcolor=$darkColor><TH colspan=2><FONT class="title2">&nbsp&nbsp;Protein Selection Rules</FONT></TH></TR>
<TR bgcolor=$lightColor><TD width=30></TD><TH width=670 align=left>
<INPUT type="checkbox" name="selProt" value=1 onclick="activateProtSel()" checked /> Select only proteins meeting the following criteria:<BR>
&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- containing at least <INPUT type="text" name="minPep" value="1" size="1"> peptide(s).<BR>
&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- and with score<SUP>*</SUP> <FONT class="title3">&#8805</FONT><INPUT type="text" name="minProtSc" value="$absMinScore" size="3" $disableMS2><BR>
&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- and with peptide coverage <FONT class="title3">&#8805</FONT><INPUT type="text" name="minCov" value="" size="2" placeholder="0">%.<BR>
<INPUT type="checkbox" name="bestMatch" value=1 /> Keep only best matching protein(s) for each match group.<BR>
<INPUT type="checkbox" name="overProt" value=1 /> Overwrite previous exclusions.
<FONT style="font-size=7px;"><BR></FONT>
</TH></TR>
<TR bgcolor=$darkColor><TD width=30></TD><TH align=left><INPUT type="checkbox" name="actTemplate" value=1 onclick="showTemplateBox()"/>Save parameters as template <INPUT type="text" name="template_name" style='visibility:hidden'>
</TH>
</TR>
<TR bgcolor=$darkColor><TH colspan=2><INPUT type="submit" name="submit" value="Proceed"/>&nbsp &nbsp &nbsp<INPUT type="button" value="Cancel" onclick="cancelAction();">
</TH>
</TR>
</TABLE><B><SUP>*</SUP>Applies to MS/MS interpretations only.
</TD>
|;
if ($ITEM eq 'ANALYSIS') {
	print "<INPUT type=\"hidden\" name=\"anaList\" value=\"$itemID\">";
}
else {
	print qq
|<TD>&nbsp</TD>
<TD valign=top><TABLE border=0 cellspacing=0 cellpadding=1>
<TR bgcolor="$darkColor"><TH class="rbBorder"><INPUT type="checkbox" onclick="checkall(this.checked)"></TH><TH class="bBorder" colspan=2>&nbsp;<FONT class="title2">Select Analyses</FONT>&nbsp;</TH></TR>
|;
	my %itemIcones=&promsConfig::getItemIcones;
	my $bgColor=($ITEM eq 'SAMPLE' || $ITEM eq 'SPOT')? $lightColor : $darkColor;
	my %prevItemName;
	foreach my $refAnaData (@itemAnalyses) {
		my ($anaID,$valStat,$msType,$dataFile,$fileFormat,$wiffFile,$taxonomy,$refDbUsed,@projHierarchy)=@{$refAnaData}; #$dbID
		$wiffFile=~s/\n//; # brakes popup otherwise
		$taxonomy='Unknown' unless $taxonomy;
		$taxonomy=~s/\(.*\)//;
		$fileFormat=~s/\..*//;
		##>Row color
		my $fatherIt=$projHierarchy[-3];
		if ($fatherIt && (!$prevItemName{$fatherIt} || $prevItemName{$fatherIt} ne $projHierarchy[-2])) { # keep color if same analysis parent item (SAMPLE or SPOT)
			$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
		}
		##>Checkbox (or -)
		#my $disabStrg=($valStat>1 || $fileFormat !~ /$selFileFormat/ || $msType ne $selMsType)? ' disabled' : '';
		my $boxStr = ($valStat > 1 || $fileFormat !~ /$selFileFormat/ || $msType ne $selMsType)? '-' : "<INPUT type=\"checkbox\" name=\"anaList\" value=\"$anaID\">"; #$disabStrg
		##>Parents
		my $parentStrg='';
		for (my $i=0;$i<=$#projHierarchy-2;$i+=2) { # stops before ana name
			my $IT=$projHierarchy[$i];
			my $itName=$projHierarchy[$i+1];
			if ($prevItemName{$IT} && $prevItemName{$IT} eq $itName) {
				$parentStrg.="<FONT color=\"$bgColor\">$itName&nbsp;&nbsp;&gt;&nbsp;&nbsp;</FONT>";
			}
			else {
				$parentStrg.="$itName&nbsp;&nbsp;&gt;&nbsp;&nbsp;";
				$prevItemName{$projHierarchy[$i]}=$itName;
				for (my $j=$i+2;$j<$#projHierarchy-1;$j+=2) {$prevItemName{$projHierarchy[$j]}='';}
			}
		}
		##>Analysis
		my $anaCode=($valStat==-1)? 'analysis:no_scan' : ($valStat==0)? 'analysis:no_val' : ($valStat==1)? 'analysis:part_val' : 'analysis:val';
		my $dbString;
		if (scalar @{$refDbUsed}==1) {$dbString="<BR><B>Databank:</B> $listDataBank{$refDbUsed->[0]}";}
		else {
			foreach my $i (0..$#{$refDbUsed}) {
				$dbString.="<BR><B>Databank #".($i+1).":</B> $listDataBank{$refDbUsed->[$i]}";
			}
		}
		print qq
|<TR valign=top bgcolor=$bgColor><TH valign=middle>$boxStr</TH><TH nowrap align=left valign=middle>&nbsp;&nbsp;$parentStrg</TH>
<TH nowrap align=left valign=middle><IMG src="$promsPath{images}/$itemIcones{$anaCode}">&nbsp;
<A href="javascript:void(null)" onmouseover="popup('<B>Type:</B> $msTypeName{$msType}<BR><B>Search file:</B> $dataFile<BR><B>MS file:</B> $wiffFile<BR><B>Search engine:</B> $fileFormat<BR><B>Taxonomy:</B> $taxonomy$dbString')" onmouseout="popout()">$projHierarchy[-1]</A>&nbsp;
</TH></TR>
|;
	}
}
print qq
|</TABLE>
</TD></TR>
</TABLE>
</FORM>
</CENTER>
|;
if ($ITEM ne 'ANALYSIS') {
	print qq
|<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup()</SCRIPT>
|;
}
print "</BODY>\n<HTML>\n";

##################################<<< SUBROUTINES >>>###########################################

#########################
####<<<ProcessSelectionForm>>>#### Already connected to DB
#########################
sub processSelectionForm {

	#######################
	####<Starting HTML>####
	#######################
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Auto Selection Process</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<BR><FONT class="title2">Performing Automated Selection:<BR><BR></FONT>
|;
	####################
	####<Parameters>####
	####################
	my @analysesList=param('anaList');
	my $MStype=param('MSTYPE');
	my $selGoodInt=(param('selGoodInt'))? 1 : 0;
	my $rejBadInt=(param('rejBadInt'))? 1 : 0;
	my $oneInt=(param('oneInt'))? 1 : 0;
	my $overPep=(param('overPep'))? 1 : 0;
	$absMinScore=(param('minScore'))? param('minScore') : 0; # => no peptide selection at all (Re-validation pass on protein selection)
	my (%minIntScore,%paramMaxRank);
	if ($selFileFormat=~/SEQUEST/) { #Take into account the score variation with the charge of the peptide
		foreach my $i (2..4) { # charges
			foreach my $j (1..3 ) { # num pep in prot
				#Definition of $minIntScore{charge:peptide} -> For Sequest research, charge influence the score (2+, 3+ and 4+)
				$minIntScore{"$i:$j"}=(param("minInt$i"."Sc$j"))? param("minInt$i"."Sc$j") : $absMinScore; $minIntScore{"$i:$j"}=$absMinScore if $minIntScore{"$i:$j"} < $absMinScore;
			}
		}
	}
	else {
		my ($paramSc1,$paramSc2,$paramSc3)=&promsMod::cleanParameters(param('minIntSc1'),param('minIntSc2'),param('minIntSc3'));
		$minIntScore{1}=($paramSc1)? $paramSc1 : $absMinScore; $minIntScore{1}=$absMinScore if $minIntScore{1} < $absMinScore;
		$minIntScore{2}=($paramSc2)? $paramSc2 : $absMinScore; $minIntScore{2}=$absMinScore if $minIntScore{2} < $absMinScore;
		$minIntScore{3}=($paramSc3)? $paramSc3 : $absMinScore; $minIntScore{3}=$absMinScore if $minIntScore{3} < $absMinScore;
	}

	$paramMaxRank{'MS1'}=(param('newMaxRankMS1'))? param('newMaxRankMS1') : 0;
	$paramMaxRank{'MS2'}=(param('newMaxRankMS2'))? param('newMaxRankMS2') : 0;
	my $noRejected=(param('noReject'))? 1 : 0;
	my ($paramMinSize,$paramMaxSize,$paramDelta,$paramMinPep,$paramMinProtSc,$paramMinCov)=&promsMod::cleanParameters(param('minSize'),param('maxSize'),param('minDelta'),param('minPep'),param('minProtSc'),param('minCov'));
	my $minSize=($paramMinSize)? $paramMinSize : 1; $minSize=1 if $minSize < 1;
	my $maxSize=($paramMaxSize)? $paramMaxSize : 1000; $maxSize=$minSize if $maxSize < $minSize;
	my $minDelta=$paramDelta? abs($paramDelta) : 10; # 0 is ignored!!!

	my $selProt=(param('selProt'))? 1 : 0;
	my $bestMatch=(param('bestMatch'))? 1 : 0;
	my $overProt=(param('overProt'))? 1 : 0;
	my $minMatchPep=($paramMinPep)? $paramMinPep : 1; $minMatchPep=1 if $minMatchPep < 1;
	my $minProtScore=($paramMinProtSc)? $paramMinProtSc : $absMinScore;  $minProtScore=$absMinScore if $minProtScore < $absMinScore; # only for MIS
	my $minPepCov=($paramMinCov)? $paramMinCov : 0;

	my $flaggedUp=(param('flaggedUp'))? 1 : 0;
	my $flaggedDown=(param('flaggedDown'))? 1 : 0;

	##################################################
	####<Validation History: Building paramString>####
	##################################################
	my $paramStrg = "selGoodInt:check:$selGoodInt;rejBadInt:check:$rejBadInt;";
	if ($selFileFormat=~/SEQUEST/){
		foreach my $i (2..4) { # charges
			foreach my $j (1..3) { # num pep in prot
				$paramStrg .= "minInt$i"."Sc$j".":text:".$minIntScore{"$i:$j"}.";";
			}
		}
	}
	else {
		$paramStrg .= "minIntSc1:text:$minIntScore{1};minIntSc2:text:$minIntScore{2};minIntSc3:text:$minIntScore{3};";
	}
	$paramStrg .= "noReject:check:$noRejected;" if $noRejected;
	$paramStrg .= 'minSize:text:'.param('minSize').';' if param('minSize');
	$paramStrg .= 'maxSize:text:'.param('maxSize').';' if param('maxSize');
	$paramStrg .= 'minDelta:text:'.param('minDelta').';' if param('minDelta');
	$paramStrg .= "newMaxRankMS1:select:$paramMaxRank{MS1};" if $paramMaxRank{'MS1'};
	$paramStrg .= "newMaxRankMS2:select:$paramMaxRank{MS2};" if $paramMaxRank{'MS2'};
	$paramStrg .= "flaggedUp:check:$flaggedUp;" if $flaggedUp;
	$paramStrg .= "flaggedDown:check:$flaggedDown;" if $flaggedDown;
	$paramStrg .= "oneInt:check:$oneInt;" if $oneInt;
	$paramStrg .= "overPep:check:$overPep;" if $overPep;
	$paramStrg .= "selProt:check:$selProt;" if $selProt;
	$paramStrg .= "minPep:text:$minMatchPep;";
	$paramStrg .= "minProtSc:text:$minProtScore;" if $MStype eq 'MIS';
	$paramStrg .= "minCov:text:$minPepCov;" if $minPepCov;
	$paramStrg .= "bestMatch:check:$bestMatch;" if $bestMatch;
	$paramStrg .= "overProt:check:$overProt;" if $overProt;

	####<Validation template>####
	if (param('actTemplate')) {
		my $quotedTemplateName=$dbh->quote(param('template_name'));
	    my ($existTemplateID) = $dbh->selectrow_array("SELECT ID_VAL_TEMPLATE FROM VALIDATION_TEMPLATE WHERE NAME=$quotedTemplateName AND ID_USER='$userID'");
	    if ($existTemplateID) { # modify template if exists
			$dbh->do("UPDATE VALIDATION_TEMPLATE SET VAL_TYPE='quali',SEARCH_ENGINE='$selFileFormat',MS_TYPE='$MStype',PARAM_STRG='$paramStrg' WHERE ID_VAL_TEMPLATE=$existTemplateID");
	    }
		else { # else inserts a new template
			$dbh->do("INSERT INTO VALIDATION_TEMPLATE(ID_USER,NAME,VAL_TYPE,SEARCH_ENGINE,MS_TYPE,PARAM_STRG) VALUES ('$userID',$quotedTemplateName,'quali','$selFileFormat','$MStype','$paramStrg')");
	    }
	    $dbh->commit;
	}

	######################################################
	####<Looping through all Analyses to be validated>####
	######################################################
#	my @analysesList=&getAnalysesList;
	my $maxAnalyses=scalar (@analysesList);
	my $countAna=0;
	foreach my $analysisID (@analysesList) {
		$countAna++;
		print "<FONT class=\"title3\">Processing Analysis $countAna/$maxAnalyses...</FONT><BR>\n";

		#my ($newMaxRank) = $dbh->selectrow_array("SELECT MAX_RANK FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
		#$newMaxRank = $paramMaxRank unless ($newMaxRank && $paramMaxRank > $newMaxRank);
		my ($anaMStype,$fileFormat)=$dbh->selectrow_array("SELECT MS_TYPE,FILE_FORMAT FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
		my $isPDM=($fileFormat=~/PDM\Z/)? 1 : 0;

		#############################################
		####<Fetching Query & Protein match data>####
		#############################################
		print "&nbsp&nbsp&nbsp<B>Fetching Peptide/Protein Match Data...";
		my (%queryStatus,%goodPeptides,%startingSelect,%noSelect,%allRankInfo,%trueRank,%queryID2Num,%queryNum2ID,%queryCharge,%pepSize,%pepScore,%pepFlag,%queryMStype);

		my $infoString=''; #'INFO_PEP0';
		foreach my $rank (1..$absMaxRank) {$infoString.=",INFO_PEP$rank";}
		my $sthQD=$dbh->prepare("SELECT ID_QUERY,QUERY_NUM,MAX_SCORE,CHARGE$infoString FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID AND VALID_STATUS >= -3");
		$sthQD->execute;

		while (my($queryID,$queryNum,$isMS2,$charge,@pepInfo)=$sthQD->fetchrow_array) { # MAX_SCORE=0 if MS1
			$queryID2Num{$queryID}=$queryNum;
			$queryNum2ID{$queryNum}=$queryID;
			$queryMStype{$queryID}=($isMS2)? 'MS2' : 'MS1';
			$queryCharge{$queryID}=($charge)? $charge : 2; # just in case: defaults to 2+ <=> lower score
			my $trueRk=1;
			my $prevScore;
			my $rank=0;
			my $refNumPep=($isPDM && $queryNum<0)? 1 : 3;
			foreach my $info (@pepInfo) {
				$rank++;
				last unless $info;
				($pepScore{"$queryNum:$rank"})=($info=~/SC=(\d*\.*\d*)/);
				if ($isMS2) { # MS/MS
					if ($prevScore) { # rank>1; peptides with same score are considered to be of same true rank !!!!!!!!!!
						$trueRk++ if $pepScore{"$queryNum:$rank"}<$prevScore;
					}
				}
				else { # MS (score=0)
					$trueRk=$rank;
				}
				$trueRank{"$queryNum:$rank"}=$trueRk;
				$prevScore=$pepScore{"$queryNum:$rank"};
				$pepFlag{"$queryNum:$rank"}= ($info=~/FLT=(-?\d)/)? $1 : 0;
				my ($sequence)=($info=~/SEQ=(\w+)/);
				$pepSize{"$queryNum:$rank"}=length($sequence);
				my ($select)=($info=~/SEL=(-*\d)/);
				if ($select==-3 || $select==-1 || ($select==-2 && $noRejected)) {	# skip if manual rejection  || better score exists ...
					$allRankInfo{$queryID}{$rank}=$info;							# || exclude already rejected interpretations
					next;
				}
				if ($selGoodInt || $rejBadInt) { # normal case
					$noSelect{$queryNum}=$rank if (($select==2 || ($select==1 && !$overPep)) && $oneInt);
					if ($selFileFormat=~/SEQUEST/) {
						my $usedCharge=($charge<=4)? $charge : 4;
						if ($isMS2 && $select<2 && ($overPep || $select==0) && $pepScore{"$queryNum:$rank"} < $minIntScore{"$usedCharge:$refNumPep"}) { # quick check for very bad MS2 peptides
							if ($rejBadInt) {$info=~s/SEL=-*\d/SEL=-2/;}
							#!!! bug fix 03/06 ->
							elsif ($overPep) { # in case pep is already selected
								$goodPeptides{"$queryNum:$rank"}=$startingSelect{"$queryNum:$rank"}=$select;
							}
							#<-!!!
						}
						else {$goodPeptides{"$queryNum:$rank"}=$startingSelect{"$queryNum:$rank"}=$select;}
					}
					else {
						if ($isMS2 && $select<2 && ($overPep || $select==0) && $pepScore{"$queryNum:$rank"} < $minIntScore{$refNumPep}) { # quick check for very bad MS2 peptides
							if ($rejBadInt) {$info=~s/SEL=-*\d/SEL=-2/;}
							#!!! bug fix 03/06 ->
							elsif ($overPep) { # in case pep is already selected
								$goodPeptides{"$queryNum:$rank"}=$startingSelect{"$queryNum:$rank"}=$select;
							}
							#<-!!!
						}
						else {$goodPeptides{"$queryNum:$rank"}=$startingSelect{"$queryNum:$rank"}=$select;}
					}
				}
				else { # no peptide (re)selection
					$goodPeptides{"$queryNum:$rank"}=$startingSelect{"$queryNum:$rank"}=$select if $select >=0;
				}
				$allRankInfo{$queryID}{$rank}=$info;
			}
		}
		$sthQD->finish;

		####<Cleaning list of auto-selectable peptides>#### rejecting already selected 2ndary peptides if $oneInt
		foreach my $queryNum (keys %noSelect) {
			foreach my $rank (1..$absMaxRank) {
				next if $rank==$noSelect{$queryNum}; # preselected rank
				if (defined($goodPeptides{"$queryNum:$rank"})) {
					delete $goodPeptides{"$queryNum:$rank"};
					$allRankInfo{$queryNum2ID{$queryNum}}{$rank}=~s/SEL=-*\d/SEL=-2/ if $rejBadInt; # \d cannot be -3
				}
			}
		}
		
		
		####<Fetching all protein hits>####
		my (%numPepProt,%pepProtHit);
		my %matchedProt;
		my $sthMP=$dbh->prepare("SELECT IDENTIFIER,QUERY_NUM,PEP_RANK,MATCH_INFO FROM RANK_PROTEIN_MATCH WHERE ID_ANALYSIS=$analysisID");
		$sthMP->execute;
		my $count=0;
		while (my ($identifier,$queryNum,$rank,$begInfo)=$sthMP->fetchrow_array) {
			$count++;
			if ($count==1000) {$count=0; print '.';}
			next if $allRankInfo{$queryNum2ID{$queryNum}}{$rank}=~/SEL=-1/; # better score exists
			push @{$matchedProt{$identifier}},"$queryNum:$rank:$begInfo"; # for protein validation
			next unless (defined($goodPeptides{"$queryNum:$rank"}));
			$numPepProt{$identifier}++;
			push @{$pepProtHit{"$queryNum:$rank"}},$identifier;
		}
		$sthMP->finish;

		####<Fetching all proteins info>####
		my %protInfo;
		my $sthPV=$dbh->prepare("SELECT IDENTIFIER,ID_PROT_VALID,SEL_STATUS,PROT_LENGTH,MAX_MATCH,MATCH_GROUP,CONF_LEVEL,NUM_MATCH FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$analysisID");
		$sthPV->execute;
		$count=0;
		while (my ($identifier,@info)=$sthPV->fetchrow_array) {
			$count++;
			if ($count==1000) {$count=0; print '.';}
			next unless $matchedProt{$identifier};
			@{$protInfo{$identifier}}=@info;
		}
		$sthPV->finish;
		print " Done.</B><BR>\n";

		#######################################
		####<Applying Interpretation Rules>####
		#######################################
		print "&nbsp&nbsp&nbsp<B>Applying Interpretation Selection rules (delta=$minDelta Da)...";
		my $modified=1;
		my $firstLoop=1;
		while ($modified==1 && ($selGoodInt || $rejBadInt)) {
			$modified=0;
			foreach my $queryID (keys %allRankInfo) { # important for locNoSelect (oneInt)
				my $queryNum=$queryID2Num{$queryID};
				my ($totRank,$selRank,$rejRank)=(0,0,0);
				my $locNoSelect; # smaller range
				if ($firstLoop) {
					$firstLoop=0;
					if ($oneInt) {
						##<Quick scan to detect already selected rank
						foreach my $rank (sort{$a<=>$b} keys %{$allRankInfo{$queryID}}) {
							if ($allRankInfo{$queryID}{$rank}=~/SEL=2/ || ($allRankInfo{$queryID}{$rank}=~/SEL=1/ && !$overPep)) {$locNoSelect=$rank; last;}
						}
					}
				}
				foreach my $rank (sort{$a<=>$b} keys %{$allRankInfo{$queryID}}) {
					next unless defined($goodPeptides{"$queryNum:$rank"});

					##<Fetching max number of co-matching peptides
					my $maxCoPep=0;
					foreach my $identifier (@{$pepProtHit{"$queryNum:$rank"}}) {
						$maxCoPep=$numPepProt{$identifier} if $maxCoPep<$numPepProt{$identifier};
					}
					my $refNumPep=(($isPDM && $queryNum<0) || $maxCoPep<=1)? 1 : ($maxCoPep==2)? 2 : 3;

					##<Applying validation rules
					my $select=$goodPeptides{"$queryNum:$rank"};
					next if $select==2;
					if ($overPep || $startingSelect{"$queryNum:$rank"}==0) {
						my ($absDelta)=($allRankInfo{$queryID}{$rank}=~/DELT=-*(\d+\.*\d*)/); # takes abs(value)!!!
						my $rejected=0;
						if ($selFileFormat=~/SEQUEST/) {
							my $usedCharge=($queryCharge{$queryID}<=4)? $queryCharge{$queryID} : 4;
							if (($queryMStype{$queryID} eq 'MS1' || $pepScore{"$queryNum:$rank"}>=$minIntScore{"$usedCharge:$refNumPep"}) && $pepSize{"$queryNum:$rank"}>=$minSize && $pepSize{"$queryNum:$rank"}<=$maxSize && $absDelta<=$minDelta && $trueRank{"$queryNum:$rank"}<=$paramMaxRank{$queryMStype{$queryID}} && ($pepFlag{"$queryNum:$rank"}==1 || $flaggedUp!=1) && ($pepFlag{"$queryNum:$rank"} !=-1 || $flaggedDown!=1 )) {
								if (defined($locNoSelect) && $rank != $locNoSelect) {$rejected=1;} # another inter already selected
								elsif ($selGoodInt) { # important to go from rank0 to last
									$select=1;
									$locNoSelect=$rank if $oneInt;
								}
							}
							else {$rejected=1;} # Pep did not pass selection
						}
						else {
							if (($queryMStype{$queryID} eq 'MS1' || $pepScore{"$queryNum:$rank"}>=$minIntScore{$refNumPep}) && $pepSize{"$queryNum:$rank"}>=$minSize && $pepSize{"$queryNum:$rank"}<=$maxSize && $absDelta<=$minDelta && $trueRank{"$queryNum:$rank"}<=$paramMaxRank{$queryMStype{$queryID}} && ($pepFlag{"$queryNum:$rank"}==1 || $flaggedUp!=1) && ($pepFlag{"$queryNum:$rank"} !=-1 || $flaggedDown!=1 )) {
								if (defined($locNoSelect) && $rank != $locNoSelect) {$rejected=1;} # another inter already selected
								elsif ($selGoodInt) { # important to go from rank0 to last
									$select=1;
									$locNoSelect=$rank if $oneInt;
								}
							}
							else {$rejected=1;} # Pep did not pass selection
						}

						my $tmpSelect=$select; #!!! bug fix 03/06
						if ($rejected) {
							$select=-2 if $rejBadInt; # failed => reject
							$tmpSelect=-2; #!!! bug fix 03/06
							foreach my $identifier (@{$pepProtHit{"$queryNum:$rank"}}) {
								$numPepProt{$identifier}--;
							}
						}
						#!!! bug fix 03/06
# 						if ($select != $goodPeptides{"$queryNum:$rank"}) { # change compare to previous loop?
# 							$allRankInfo{$queryID}{$rank}=~s/SEL=-*\d/SEL=$select/;
# 							$modified=1;
# 							if ($rejected) {delete $goodPeptides{"$queryNum:$rank"};}
# 							else {$goodPeptides{"$queryNum:$rank"}=$select;}
# 						}
						if ($tmpSelect != $goodPeptides{"$queryNum:$rank"}) { # change compare to previous loop?
							$modified=1;
							$allRankInfo{$queryID}{$rank}=~s/SEL=-*\d/SEL=$select/;
							if ($rejected) { # $tmpSelect=-2
								delete $goodPeptides{"$queryNum:$rank"};
								$allRankInfo{$queryID}{$rank}=~s/SEL=1/SEL=0/ unless $rejBadInt; # deselect if selected but do not reject
							}
							else {$goodPeptides{"$queryNum:$rank"}=$select;}
						}
						#!!!
					}
				}
			}
		}
		####<Updating query selection status>####
		foreach my $queryID (keys %allRankInfo) {
			my ($totRank,$selRank,$rejRank,$badRank)=(0,0,0,0);
			foreach my $info (values %{$allRankInfo{$queryID}}) {
				next unless $info; # no interpretation for rank
				$totRank++;
				if ($info=~/SEL=-1/) {$badRank++;} # better score exists
				elsif ($info=~/SEL=[12]/) {$selRank++;}
				elsif ($info=~/SEL=-[23]/) {$rejRank++;}
			}
			###<Computing VALID_STATUS
			if ($badRank==$totRank) {$queryStatus{$queryID}=-3;}
			elsif ($selRank) {$queryStatus{$queryID}=$selRank;}
			elsif ($rejRank==$totRank) {$queryStatus{$queryID}=0;}
			elsif ($rejRank) {$queryStatus{$queryID}=-2;} # rejected + not verified
			else {$queryStatus{$queryID}=-1;}
		}

		################################
		####<Updating queries in DB>####
		################################
		print '.';
		my $updateString='VALID_STATUS=?';
		foreach my $rank (1..$absMaxRank) {$updateString.=",INFO_PEP$rank=?";}
		my $sthUpQ=$dbh->prepare("UPDATE QUERY_VALIDATION SET $updateString WHERE ID_QUERY=?");
		foreach my $queryID (sort{$a<=>$b} keys %queryStatus) {
			my @rankData;
			foreach my $rank (1..$absMaxRank) {
				if ($allRankInfo{$queryID}{$rank}) {push @rankData,$allRankInfo{$queryID}{$rank};}
				else {
					push @rankData,undef;
					delete $allRankInfo{$queryID}{$rank}; # just to be safe
				}
			}
# print "><B>qID=</B>$queryID <B>VS=</B>$queryStatus{$queryID} <B>PEP: </B>",join (' * ',@rankData),' (',scalar keys %{$allRankInfo{$queryID}},")<BR>\n";
			$sthUpQ->execute($queryStatus{$queryID},@rankData,$queryID);
		}
		$sthUpQ->finish;
		print " Done.</B><BR>\n";

		###########################################
		####<Updating protein selection status>####
		###########################################
		if ($selProt || $bestMatch) {print "&nbsp&nbsp&nbsp<B>Applying Protein Selection rules...";}
		else {print "&nbsp&nbsp&nbsp<B>Updating Protein Selection Status...";}
		my %groupMaxMatch;
		foreach my $identifier (keys %matchedProt) {
			# -3 and no overProt -2 already removed
			next unless $protInfo{$identifier};
			my ($protID,$selStatus,$protLength,$maxMatch,$matchGroup,$oldConfLevel,$oldNumMatch)=@{$protInfo{$identifier}};
			my ($numMatch,$protScore,$rejectedMatch)=(0,0,0);
			my %posPeptide;
			foreach my $matchInfo (@{$matchedProt{$identifier}}) {
				my ($queryNum,$rank,@aaData)=split(/:/,$matchInfo);
				my $queryID=$queryNum2ID{$queryNum};
				if ($allRankInfo{$queryID}{$rank}=~/SEL=[12]/) {
					$numMatch++;
					#$protScore+=(scalar (@aaData) * $pepScore{"$queryNum:$rank"}); # matchFreq * pepScore
					$protScore+=$pepScore{"$queryNum:$rank"}; # pepScore
				}
				elsif ($allRankInfo{$queryID}{$rank}=~/SEL=-[23]/) {$rejectedMatch++;}
				foreach my $aaStr (@aaData) { # aaStr=begAa,flankNterAA,flankCterAA
					my $beg=(split(/,/,$aaStr))[0];
					$posPeptide{$beg}=$pepSize{"$queryNum:$rank"};
				}
			}
			$protInfo{$identifier}[6]=$numMatch; # numMatch
			push @{$protInfo{$identifier}},$protScore; # score (index 7)
			if ($numMatch>$oldNumMatch) {$protInfo{$identifier}[5]=2;} # confLevel
			else {$protInfo{$identifier}[5]=($oldConfLevel)? $oldConfLevel : 2;} # else unchanged

			####<Applying exclusion rules>####
			if ($protInfo{$identifier}[1]>-3 && ($overProt || $protInfo{$identifier}[1]>-2)) { # !filtered or overProt+excluded
# print "**$identifier => $protInfo{$identifier}[1]<BR>\n";
				if ($numMatch) { # at least 1 peptide selected
					my $excludeProt;
					if ($selProt) {
						##<Min numMatch & min score
						if ($numMatch<$minMatchPep || ($anaMStype eq 'MIS' && $protScore<$minProtScore)) {$excludeProt=1;}
						##<Min peptide coverage
						elsif ($minPepCov) {
							my $coverage=0;
							my $curBeg=0;
							my $curEnd=0;
							unless ($protLength) { # set protLength=last covering peptide
								foreach my $aa (keys %posPeptide) {
									$protLength=$aa+$posPeptide{$aa}-1 if $protLength<$aa+$posPeptide{$aa}-1;
								}
							}
							foreach my $aa (sort{$a<=>$b} keys %posPeptide) {
								if ($aa>$curEnd) {
									$coverage+=($curEnd-$curBeg+1) unless $curBeg==0;
									$curBeg=$aa;
									$curEnd=$aa+$posPeptide{$aa}-1;
								}
								elsif ($aa+$posPeptide{$aa}-1>$curEnd) {
									$curEnd=$aa+$posPeptide{$aa}-1;
								}
							}
							$coverage+=($curEnd-$curBeg+1);
							$coverage=100*($coverage/$protLength);
							$excludeProt=1 if $coverage<$minPepCov;
						}
					}
					if ($excludeProt) {$protInfo{$identifier}[1]=-2;}
					elsif ($selStatus>=-1 || $overProt) {
						$protInfo{$identifier}[1]=($numMatch+$rejectedMatch==$maxMatch)? 2 : 1;
					}
				}
				else { # no peptides selected
					$protInfo{$identifier}[1]=($rejectedMatch==$maxMatch)? 0 : -1;
				}

				###<Best Group Matches (Only validated proteins are considered)>###
				$groupMaxMatch{$matchGroup}=$numMatch if ($protInfo{$identifier}[1]>=1 && (!$groupMaxMatch{$matchGroup} || $groupMaxMatch{$matchGroup}<$numMatch));
			}
		}

		####<Exclude if not best matching proteins in group>####
		if ($bestMatch) { # 2nd pass to exclude prot with less peptides than best prot in MGr
			foreach my $identifier (keys %matchedProt) {
				$protInfo{$identifier}[1]=-2 if ($protInfo{$identifier}[1]>=1 && $protInfo{$identifier}[6]<$groupMaxMatch{$protInfo{$identifier}[4]});
			}
		}

		#################################
		####<Updating Proteins in DB>####
		#################################
		print '.';
		my $sthUpProt=$dbh->prepare("UPDATE PROTEIN_VALIDATION SET SEL_STATUS=?,CONF_LEVEL=?,NUM_MATCH=?,SCORE=? WHERE ID_PROT_VALID=?");
		foreach my $identifier (keys %matchedProt) {
# print "($protInfo{$identifier}[1], @{$protInfo{$identifier}}[5..7], $protInfo{$identifier}[0])<BR>\n";
			$sthUpProt->execute($protInfo{$identifier}[1],@{$protInfo{$identifier}}[5..7],$protInfo{$identifier}[0]);
		}
		$sthUpProt->finish;

		print " Done.</B><BR><BR>\n";
		$dbh->commit;
		#my @paramList = ($selGoodInt,$rejBadInt,$minIntScore{1},$minIntScore{2},$minIntScore{3},$noRejected,param('minSize'),param('maxSize'),param('minDelta'),$paramMaxRank,$flaggedUp,$flaggedDown,$oneInt,$overPep,$selProt,$minMatchPep,$minProtScore,$minPepCov,$bestMatch,$overProt);
		&promsMod::updateAnalysisHistory($dbh,$analysisID,$paramStrg,'quali');
	}

	print "<FONT class=\"title2\">Automated Selection is Complete.</FONT><BR>\n";

	$dbh->disconnect;

	####<Reload all frames>####
	&reloadFrames;

	exit;
}


#####################################################
####<<<Restoring last 'saved' validation level>>>####
#####################################################
sub restoreValidation { # Only for 1 analysis so far...

	####<Starting HTML>####
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Restore Selection</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
</STYLE>
|;
	if ($ITEM eq 'ANALYSIS') {
		print qq
|<SCRIPT LANGUAGE="JavaScript">
top.promsFrame.spectrumFrame.location="$promsPath{html}/nothing.html";
</SCRIPT>
|;
	}
	print qq
|</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<BR><FONT class="title2">Restoring last saved selection:<BR><BR></FONT>
|;

	#############################
	####<Clearing validation>####
	#############################
	&clearValidations(2);

	############################################
	####<Restoring peptide/query validation>####
	############################################
	my $sthSelPep=$dbh->prepare("SELECT ID_PEPTIDE,QUERY_NUM,PEP_RANK,VALID_STATUS,COMMENTS FROM PEPTIDE WHERE ID_ANALYSIS=$itemID");
	#my $sthSelR0=$dbh->prepare("SELECT SCORE,MISS_CUT,MR_CALC,MR_DELTA,PEP_SEQ FROM PEPTIDE WHERE ID_PEPTIDE=?");
	#my $sthSelRM=$dbh->prepare("SELECT ID_PROTEIN,PEP_BEG FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_PEPTIDE=?");
	my %sthSelRank;
	my %sthUpRank;
	foreach my $rank (1..$absMaxRank) {
		$sthSelRank{$rank}=$dbh->prepare("SELECT ID_QUERY,VALID_STATUS,INFO_PEP$rank FROM QUERY_VALIDATION WHERE QUERY_NUM=? AND ID_ANALYSIS=$itemID");
		$sthUpRank{$rank}=$dbh->prepare("UPDATE QUERY_VALIDATION SET INFO_PEP$rank=?,VALID_STATUS=? WHERE ID_QUERY=?");
	}
	my %rankList; #,%rank0Proteins,%rank0Score);
	$sthSelPep->execute;
	while (my ($peptideID,$queryNum,$rank,$select,$pepComments)=$sthSelPep->fetchrow_array) {
		$rankList{"$queryNum:$rank"}=1;

		###<Fetching validated data>###
		$sthSelRank{$rank}->execute($queryNum);
		my ($queryID,$validStatus,$rankInfo)=$sthSelRank{$rank}->fetchrow_array;
		###<Rank 0: Recreate rank info>##
		#if ($rank==0) {
		#	$sthSelR0->execute($peptideID);
		#	my ($score,$missCut,$MrCalc,$delta,$sequence)=$sthSelR0->fetchrow_array;
		#	$sthSelRM->execute($peptideID);
		#	#<Listing matched proteins>#
		#	my %numProteins;
		#	while (my ($proteinID,$pepBeg,$pepEnd)=$sthSelRM->fetchrow_array) {
		#		push @{$rank0Proteins{$proteinID}{$queryNum}},$pepBeg; # rank=0
		#		$numProteins{$proteinID}=1;
		#	}
		#	my $match=scalar keys %numProteins;
		#	$rankInfo="SEL=1,MIS=$missCut,CALC=$MrCalc,DELT=$delta,SEQ=$sequence,SC=$score,MATCH=$match,";
		#	$rank0Score{$queryNum}=$score; # needed later for (MAX_)SCORE
		#}
		###<Normal rank>##
		#else {
		$select=1 unless $select; # default=auto
		$rankInfo=~s/SEL=0/SEL=$select/;
		#}

		###<Peptide comments>###
		$rankInfo.="COM=$pepComments" if $pepComments;

		###<Updating QUERY_VALIDATION>###
		$validStatus=($validStatus<=-1)? 1 : $validStatus+1;
		$sthUpRank{$rank}->execute($rankInfo,$validStatus,$queryID);
	}
	$sthSelPep->finish;
	#$sthSelR0->finish;
	#$sthSelRM->finish;
	foreach my $rank (keys %sthSelRank) {
		$sthSelRank{$rank}->finish;
		$sthUpRank{$rank}->finish;
	}

	################################
	####<Fetching protein lists>####
	################################
	####<Fetching list of effectively validated proteins>####
	my (%validProteins,%validProtIdent);
	my $sthSelVP=$dbh->prepare("SELECT PROTEIN.ID_PROTEIN,IDENTIFIER,NUM_PEP,SCORE,CONF_LEVEL FROM ANALYSIS_PROTEIN,PROTEIN WHERE ANALYSIS_PROTEIN.ID_PROTEIN=PROTEIN.ID_PROTEIN AND ID_ANALYSIS=$itemID");
	$sthSelVP->execute;
	while (my ($proteinID,$identifier,@protInfo)=$sthSelVP->fetchrow_array) {
		@{$validProteins{$identifier}}=($proteinID,@protInfo);
		$validProtIdent{$proteinID}=$identifier;
	}
	$sthSelVP->finish;

	####<Fetching list of theorically validated proteins>#### (including excluded ones if any)
	my %theoValidProteins;
	#my @missingMatches;
	my $sthSelTP=$dbh->prepare("SELECT IDENTIFIER FROM RANK_PROTEIN_MATCH WHERE QUERY_NUM=? AND PEP_RANK=? AND ID_ANALYSIS=$itemID");
	foreach my $qrKey (keys %rankList) {
		my ($queryNum,$rank)=split(/:/,$qrKey);
		$sthSelTP->execute($queryNum,$rank);
		#my $match=0;
		while (my ($identifier)=$sthSelTP->fetchrow_array) {
			@{$theoValidProteins{$identifier}{$qrKey}}=1;
			#$match=1;
		}
		#push @missingMatches,$qrKey unless $match; # happens only for rank0
	}
	$sthSelVP->finish;

	####################################################
	####<Adding rank0 matches to RANK_PROTEIN_MATCH>#### can be updated before PROTEIN_VALIDATION because requires 'IDENTIFIER' not ID_PROT_VALID
	####################################################
	#my $sthInsM0=$dbh->prepare("INSERT INTO RANK_PROTEIN_MATCH (IDENTIFIER,QUERY_NUM,MATCH_MULTI,MATCH_INFO,PEP_RANK,ANALYSIS_ID) VALUES (?,?,?,?,0,$itemID)");
	#foreach my $proteinID (keys %rank0Proteins) {
	#	foreach my $queryNum (keys %{$rank0Proteins{$proteinID}}) {
	#		my $matchMulti=scalar @{$rank0Proteins{$proteinID}{$queryNum}};
	#		my $matchInfo=join(':',@{$rank0Proteins{$proteinID}{$queryNum}});
	#		$sthInsM0->execute($validProtIdent{$proteinID},$queryNum,$matchMulti,$matchInfo);
	#	}
	#}
	#$sthInsM0->finish;

	######################################
	####<Restoring protein validation>#### &clearValidations does not remove filters!
	######################################

	####<Updating existing proteins in PROTEIN_VALIDATION>####
	my $sthSelMM=$dbh->prepare("SELECT ID_PROT_VALID,MAX_MATCH,MAX_SCORE FROM PROTEIN_VALIDATION WHERE IDENTIFIER=? AND ID_ANALYSIS=$itemID");
	my $sthUpVP=$dbh->prepare("UPDATE PROTEIN_VALIDATION SET SEL_STATUS=?,NUM_MATCH=?,SCORE=?,CONF_LEVEL=? WHERE ID_PROT_VALID=?");
	foreach my $identifier (keys %theoValidProteins) {
		$sthSelMM->execute($identifier);
		my ($protValidID,$maxMatch,$maxScore)=$sthSelMM->fetchrow_array;
		##<Effectively validated protein
		if ($validProteins{$identifier}) {
			#<Also match by rank0 : Correct MAX_MATCH & MAX_SCORE
			#my $proteinID=$validProteins{$identifier}[0];
			#if ($rank0Proteins{$proteinID}) {
			#	foreach my $queryNum (keys %{$rank0Proteins{$proteinID}}) {
			#		$maxMatch++;
			#		$maxScore+=$rank0Score{$queryNum};
			#	}
			#}
			my $selStatus=($validProteins{$identifier}[1]==$maxMatch)? 2 : 1;
			$sthUpVP->execute($selStatus,$validProteins{$identifier}[1],$validProteins{$identifier}[2],$validProteins{$identifier}[3],$protValidID);
		}
		##>Rejected protein
		else {
			# WARNING: Sequence must be scanned to check if matched by rank0 (MAX_MATCH & MAX_SCORE will be affected)   (<--- rank0 obsolete)!!!!!!
			$sthUpVP->execute(-2,0,0,2,$protValidID);
		}
	}
	$sthSelMM->finish;
	$sthUpVP->finish;

	#####<Adding proteins ONLY matched by missing rank0 to PROTEIN_VALIDATION>#### !!!!RANK0 OBSOLETE!!!!!
	#my (%proteinSeq,%queryMatchGroup);
	#my ($maxTempProtID)=$dbh->selectrow_array("SELECT MAX(ID_PROT_VALID) FROM PROTEIN_VALIDATION");
	#my ($dataFile)=$dbh->selectrow_array("SELECT DATA_FILE FROM ANALYSIS WHERE ID_ANALYSIS=$itemID");
	#my ($matchGroup)=$dbh->selectrow_array("SELECT max(MATCH_GROUP) FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$itemID");
	#my $sthSelM0=$dbh->prepare("SELECT COUNT(*) FROM RANK_PROTEIN_MATCH WHERE IDENTIFIER=? AND ID_ANALYSIS=$itemID");
	#my $sthSelPI=$dbh->prepare("SELECT PROT_SEQ,PROT_DES,PROT_LENGTH,MW,ORGANISM FROM PROTEIN WHERE ID_PROTEIN=?");
	#my $sthInsMP=$dbh->prepare("INSERT INTO PROTEIN_VALIDATION (ID_PROT_VALID,IDENTIFIER,PROT_DES,PROT_LENGTH,MW,ORGANISM,NUM_MATCH,MAX_MATCH,SCORE,MAX_SCORE,MATCH_GROUP,SEL_STATUS,DATA_FILE,ID_ANALYSIS) VALUES (?,?,?,?,?,?,?,?,?,?,2,$dataFile,$itemID)");
	#foreach my $identifier (keys %validProteins) {
	#	my $proteinID=$validProteins{$identifier}[0];
	#	next unless $rank0Proteins{$proteinID}; # prot is matched by other ranks
	#	###<Checking if prot is only matched by rank 0>###
	#	my ($numMatch)=$sthSelM0->fetchrow_array;
	#	if ($numMatch==scalar keys %{$rank0Proteins{$proteinID}}) {
	#		my $score=0;
	#		foreach my $queryNum (keys %{$rank0Proteins{$proteinID}}) {
	#			$score+=$rank0Score{$queryNum};
	#			$queryMatchGroup{$queryNum}=++$matchGroup unless $queryMatchGroup{$queryNum};
	#		}
	#		##<Fetching other info>##
	#		$sthSelPI->execute($proteinID);
	#		($proteinSeq{$identifier},my @protInfo)=$sthSelPI->fetchrow_array;
	#		##<Inserting protein>##
	#		$sthInsMP->execute(++$maxTempProtID,$identifier,@protInfo,$numMatch,$numMatch,$score,$score,$matchGroup); # NUM_MATCH=MAX_MATCH idem for scores
	#	}
	#}
	#$sthSelM0->finish;
	#$sthSelPI->finish;
	#$sthInsMP->finish;
	$dbh->commit;
	$dbh->disconnect;

	##################################################
	####<Adding rank0-only proteins to fasta file>#### !!!!RANK0 OBSOLETE!!!!!
	##################################################

	#####<checking if rank0-only proteins have not already beed added to fasta file (multiple Restores)>####
	#if (scalar keys %proteinSeq) {
	#	open (FAS,">>$promsPath{valid}/ana_$itemID/analysis.fasta");
	#	while (<FAS>) {
	#		if (/^>(\S+)/) {
	#			my @line=split(//,$1);
	#			foreach my $identifier (@line) {
	#				delete $proteinSeq{$identifier} if $proteinSeq{$identifier};
	#			}
	#		}
	#	}
	#	close FAS;
	#}

	#####<Adding remaining sequences from rank0-only proteins to fasta file>####
	#if (scalar keys %proteinSeq) {
	#	open (FAS,">>$promsPath{valid}/ana_$itemID/analysis.fasta");
	#	foreach my $identifier (keys %proteinSeq) {
	#		print FAS ">$identifier\n$proteinSeq{identifier}\n";
	#	}
	#	close FAS;
	#}

	###########################
	####<Reload all frames>####
	###########################
	&reloadFrames;
}

#########################################################
####<<<Select auto- or all for clearing selections>>>####
#########################################################
sub selectClearLevel {
	####<Starting HTML>####
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Select Clear Level</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
function clearSelection() {
	parent.spectrumFrame.location="$promsPath{cgi}/autoSelect.cgi?id_project=$projectID&ITEM=$ITEM&ID=$itemID&MSTYPE=$selMsType&ACT="+document.getElementById('ACT').value;
}
function cancelAction() {
	if (parent.selectedView=='myList') {parent.writeMenu(parent.selectedProtID,2);}
	else {parent.setView(parent.selectedView);}
}
top.promsFrame.spectrumFrame.location="$promsPath{html}/nothing.html";
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Clear Peptide/Protein Selection</FONT>
<BR><BR>
<FONT class="title3">Clear <SELECT id="ACT" style="font-weight:bold;font-size:14px"><OPTION value="clearAuto">Auto-</OPTION><OPTION value="clearAll">All</OPTION></SELECT>
selection in current analysis</FONT>
<BR><BR><INPUT type="button" value="Proceed" onclick="clearSelection()"/>&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="cancelAction();\">
</CENTER>
</BODY>
</HTML>
|;
	exit;
}
###############################################
####<<<Clearing all previous validations>>>#### Already connected to DB
###############################################
sub clearValidations {
	my $clearLevel=$_[0];
	my $clearStrg=($clearLevel==1)? 'auto-selection' : 'all selection';

	if ($action=~/clear/) { # also called by &restoreValidation
		####<Starting HTML>####
		print header(-'content-encoding'=>'no',-charset=>'utf-8');
		warningsToBrowser(1);
		print qq
|<HTML>
<HEAD>
<TITLE>Clear Selection</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<BR><FONT class="title2">Deleting $clearStrg data:<BR><BR></FONT>
|;
	}

	####################################################
	####<Looping through all Analyses to be cleared>####
	####################################################
	my $sthMsType=$dbh->prepare("SELECT MS_TYPE FROM ANALYSIS WHERE ID_ANALYSIS=?");
	#my @analysesList=&getAnalysesList;
	my $maxAnalyses=scalar(@analysisList);
	my $countAna=0;
	foreach my $analysisID (@analysisList) {
		$countAna++;
		print "<FONT class=\"title3\">Processing Analysis $countAna/$maxAnalyses...</FONT><BR>\n";

		$sthMsType->execute($analysisID);
		my ($msType)=$sthMsType->fetchrow_array; # overwrites global

		####<Setting all validatable queries to 'not verified' (VALID_STATUS=-1)>#### (RANK 0 set to NULL)
		print "&nbsp&nbsp&nbsp<B>Clearing peptide selections...";
		my $infoString='ID_QUERY,QUERY_NUM';
		my %sthRank;
		foreach my $rank (1..$absMaxRank) {
			$infoString.=",INFO_PEP$rank";
			#if ($rank==0) {
			#	$sthRank{0}=$dbh->prepare("UPDATE QUERY_VALIDATION SET INFO_PEP0=NULL WHERE ID_QUERY=?");
			#}
			#else {
				$sthRank{$rank}=$dbh->prepare("UPDATE QUERY_VALIDATION SET INFO_PEP$rank=? WHERE ID_QUERY=?");
			#}
		}
		print '.';
		my $sthSelQ=$dbh->prepare("SELECT $infoString FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID AND VALID_STATUS>=-3");
		my $sthUpQ=$dbh->prepare("UPDATE QUERY_VALIDATION SET VALID_STATUS=? WHERE ID_QUERY=?");

		$sthSelQ->execute;
		print '.';
		#my %rank0Score; # records the scores of rank 0 peptides
		my %selectedRanks;
		while (my ($queryID,$queryNum,@pepInfoList)=$sthSelQ->fetchrow_array) {
			my %pepInfo;
			my $i=0;
			foreach my $rank (1..$absMaxRank) {
				$pepInfo{$rank}=$pepInfoList[$i];
				$i++;
			}
			my $validStatus;
			my ($totRank,$badRank,$selRank,$rejRank)=(0,0,0,0);
			#if ($pepInfo{0}) { # rank0 defined
			#	#($rank0Score{$queryNum})=($pepInfo{0}=~/SC=(\d+\.*\d*)/);
			#	$sthRank{0}->execute($queryID);
			#}
			foreach my $rank (1..$absMaxRank) {
				next unless $pepInfo{$rank};
				$totRank++;
				my ($select)=($pepInfo{$rank}=~/SEL=(-*\d)/);
				$pepInfo{$rank}=~s/COM=.*// if $clearLevel==2;
				if ($clearLevel==2 && ($select>=1 || $select<=-2)) { # selected or rejected (not -1 and 0)
					$pepInfo{$rank}=~s/SEL=$select/SEL=0/; # 0=never verified
					$sthRank{$rank}->execute($pepInfo{$rank},$queryID);
				}
				elsif ($clearLevel==1) {
					if ($select==1 || $select==-2) { # selected or rejected
						$pepInfo{$rank}=~s/SEL=-*\d/SEL=0/; # 0=never verified
						$sthRank{$rank}->execute($pepInfo{$rank},$queryID);
					}
					elsif ($select==2) {
						$selRank++;
						($selectedRanks{"$queryNum:$rank"})=($pepInfo{$rank}=~/SC=(\d+\.*\d*)/);
					}
					elsif ($select==-3) {$rejRank++;}
				}
				#if ($pepInfo{$rank}=~/SEL=1/ || $pepInfo{$rank}=~/SEL=-2/) { # selected or rejected
				#	$pepInfo{$rank}=~s/SEL=-*\d/SEL=0/; # 0=never verified
				#	$sthRank{$rank}->execute($pepInfo{$rank},$queryID);
				#}
				if ($select==-1) {$badRank++;}
			}
			$validStatus=($totRank==0)? -4 : ($badRank==$totRank)? -3 : ($selRank)? $selRank : ($rejRank)? -2 : -1;
			$sthUpQ->execute($validStatus,$queryID);
		}
		foreach my $rank (keys %sthRank) {$sthRank{$rank}->finish;}
		$sthUpQ->finish;
		#print '.';
		#print " Done.</B><BR>\n";

		#####<Updating protein matched by RANK 0>#### !!!!RANK0 OBSOLETE!!!!!
		#print "&nbsp&nbsp&nbsp<B>Clearing protein selections...";
		#my $sthSelM=$dbh->prepare("SELECT IDENTIFIER,QUERY_NUM,MATCH_MULTI FROM RANK_PROTEIN_MATCH WHERE PEP_RANK=0 AND ID_ANALYSIS=$analysisID");
		#my $sthSelP=$dbh->prepare("SELECT ID_PROT_VALID,MAX_MATCH,MAX_SCORE FROM PROTEIN_VALIDATION WHERE IDENTIFIER=? AND ID_ANALYSIS=$analysisID");
		#my $sthUpM=$dbh->prepare("UPDATE PROTEIN_VALIDATION SET MAX_MATCH=?,MAX_SCORE=? WHERE ID_PROT_VALID=?");
		#$sthSelM->execute();
		#print '.';
		#while (my ($identifier,$queryNum,$matchMulti)=$sthSelM->fetchrow_array) {
		#	$sthSelP->execute($identifier);
		#	my ($protID,$maxMatch,$maxScore)=$sthSelP->fetchrow_array;
		#	my $newMaxMatch=$maxMatch-1;
		#	my $newMaxScore=$maxScore-($rank0Score{$queryNum}*$matchMulti);
		#	$sthUpM->execute($newMaxMatch,$newMaxScore,$protID);
		#}
		###<Deleting protein matches with RANK 0>##
		#print '.';
		#$dbh->do("DELETE FROM RANK_PROTEIN_MATCH WHERE PEP_RANK=0 AND ID_ANALYSIS=$analysisID");
		###<Deleting proteins only matched by RANK 0>##
		#print '.';
		#$dbh->do("DELETE FROM PROTEIN_VALIDATION WHERE MAX_MATCH=0 AND ID_ANALYSIS=$analysisID");

		####<Setting all non filtered proteins to 'not verified' (SEL_STATUS=-1)>####
		print '.';
		my $scoreString=($msType eq 'PMF')? '' : ',SCORE=0';
		$dbh->do("UPDATE PROTEIN_VALIDATION SET SEL_STATUS=-1,NUM_MATCH=0$scoreString,CONF_LEVEL=NULL WHERE SEL_STATUS>-3 AND ID_ANALYSIS=$analysisID");

		###<Computing remaining protein/peptide matches (if any)
		if ($clearLevel==1 && scalar keys %selectedRanks) {
			my (%selectedProteins,%proteinInfo);
			my $sthSelM=$dbh->prepare("SELECT IDENTIFIER,MATCH_MULTI FROM RANK_PROTEIN_MATCH WHERE ID_ANALYSIS=$analysisID AND QUERY_NUM=? AND PEP_RANK=?");
			foreach my $peptide (keys %selectedRanks) {
				my ($queryNum,$rank)=split(':',$peptide);
				$sthSelM->execute($queryNum,$rank);
				while (my ($identifier,$matchMulti)=$sthSelM->fetchrow_array) {
					$selectedProteins{$identifier}[0]++;
					$selectedProteins{$identifier}[1]+=$selectedRanks{$peptide}*$matchMulti;
				}
			}
			$sthSelM->finish;

			my $sthSelP=$dbh->prepare("SELECT ID_PROT_VALID,IDENTIFIER,MAX_MATCH FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$analysisID");
			$sthSelP->execute;
			while (my ($protValidID,$identifier,$maxMatch)=$sthSelP->fetchrow_array) {
				next unless $selectedProteins{$identifier};
				@{$proteinInfo{$identifier}}=($protValidID,$maxMatch);
			}
			$sthSelP->finish;

			my $sthUpP=$dbh->prepare("UPDATE PROTEIN_VALIDATION SET SEL_STATUS=?,NUM_MATCH=?,SCORE=?,CONF_LEVEL=2 WHERE ID_PROT_VALID=?");
			foreach my $identifier (keys %selectedProteins) {
				my $selStatus=($selectedProteins{$identifier}[0]==$proteinInfo{$identifier}[1])? 2 : 1;
				$sthUpP->execute($selStatus,@{$selectedProteins{$identifier}},$proteinInfo{$identifier}[0]);
			}
			$sthUpP->finish;
		}

		print " Done.</B><BR><BR>\n";
		$dbh->commit;
		if($clearLevel == 1){
			&promsMod::updateAnalysisHistory($dbh,$analysisID,'','clear_auto');
		} elsif ($clearLevel == 2){
			&promsMod::updateAnalysisHistory($dbh,$analysisID,'','clear_all');
		}
	}
	$sthMsType->finish;

	if ($action=~/clear/) {
		print "<FONT class=\"title2\">Selection data have been deleted.</FONT><BR>\n";

		$dbh->disconnect;

		####<Reload all frames>####
		&reloadFrames;
	}
}


####################################################
####<<<Making low scoring peptides selectable>>>####
####################################################
sub activateLowScores {

	####<Starting HTML>####
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Activate Peptides</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<BR><FONT class="title2">Activating all lower-scoring interpretations:<BR><BR></FONT>
|;

	####<Fetching list of MIS analyses>####
	#my @analysesMIS;
	#if ($ITEM eq 'ANALYSIS') {push @analysesMIS,$itemID if $msType eq 'MIS';}
	#else {@analysesMIS=&getAnalysesList('MIS');}
	#my $numAnalyses=scalar (@analysesMIS);
	#my $countAna=0;

	####<Looping through analyses>####
	my (%sthLS,%sthUpLS);
	foreach my $rank (1..$absMaxRank) {
		$sthLS{$rank}=$dbh->prepare("SELECT INFO_PEP$rank,ID_QUERY,QUERY_NUM,VALID_STATUS FROM QUERY_VALIDATION WHERE INFO_PEP$rank LIKE '%SEL=-1%' AND ID_ANALYSIS=?");
		$sthUpLS{$rank}=$dbh->prepare("UPDATE QUERY_VALIDATION SET VALID_STATUS=?,INFO_PEP$rank=? WHERE ID_QUERY=?");
	}
	my $sthMax=$dbh->prepare("SELECT QUERY_NUM,PEP_RANK,ID_PROT_VALID,MAX_MATCH,MAX_SCORE FROM PROTEIN_VALIDATION PV,RANK_PROTEIN_MATCH RPM WHERE PV.IDENTIFIER=RPM.IDENTIFIER AND PV.ID_ANALYSIS=? AND RPM.ID_ANALYSIS=?");
	my $sthUpMax=$dbh->prepare("UPDATE PROTEIN_VALIDATION SET MAX_MATCH=?,MAX_SCORE=? WHERE ID_PROT_VALID=?");
	my $sthUpAna=$dbh->prepare("UPDATE ANALYSIS SET LOWER_SCORES=2 WHERE ID_ANALYSIS=?");

	print "<FONT class=\"title3\">\n";
	my $numAnalyses=scalar @analysisList;
	my $countAna=0;
	foreach my $anaID (@analysisList) { #@analysesMIS
		$countAna++;
		print "+Processing Analysis $countAna/$numAnalyses...";
		##<Fetching and updating lowser-scoring ranks
		my %rankScores;
		foreach my $rank (1..$absMaxRank) {
			print '.';
			$sthLS{$rank}->execute($anaID);
			while (my ($pepInfo,$queryID,$queryNum,$valStatus)=$sthLS{$rank}->fetchrow_array) {
				$pepInfo=~s/SEL=-1/SEL=0/;
				$valStatus=-1 if $valStatus==-3;
				$sthUpLS{$rank}->execute($valStatus,$pepInfo,$queryID);
				(@{$rankScores{"$queryNum:$rank"}})=($pepInfo=~/SC=(-?\d+\.?\d*)/);
			}
		}
		##<Fetching matching proteins (match frequency is not considered)
		print '+';
		my (%proteinMaxMatch,%proteinMaxScore);
		my $pCount=0;
		$sthMax->execute($anaID,$anaID);
		while (my ($queryNum,$rank,$protID,$maxMatch,$maxScore)=$sthMax->fetchrow_array) {
			$pCount++;
			if ($pCount==500) {
				print '.';
				$pCount=1;
			}
			next unless $rankScores{"$queryNum:$rank"};
			unless ($proteinMaxScore{$protID}) {
				$proteinMaxMatch{$protID}=$maxMatch;
				$proteinMaxScore{$protID}=$maxScore;
			}
			$proteinMaxMatch{$protID}++;
			$proteinMaxScore{$protID}+=$rankScores{"$queryNum:$rank"};
		}
		##<Updating matching proteins
		print '+';
		$pCount=0;
		foreach my $protID (keys %proteinMaxMatch) {
			$pCount++;
			if ($pCount==500) {
				print '.';
				$pCount=1;
			}
			$sthUpMax->execute($proteinMaxMatch{$protID},$proteinMaxScore{$protID},$protID);
		}
		$sthUpAna->execute($anaID);

		promsMod::updateAnalysisHistory($dbh, $anaID, undef, 'lowP_a');

		$dbh->commit;
		print " Done.<BR><BR>\n";
	}
	print "</FONT>\n";
	print "<FONT class=\"title2\">All peptides are now selectable.</FONT><BR>\n";

	foreach my $rank (1..$absMaxRank) {$sthLS{$rank}->finish; $sthUpLS{$rank}->finish;}
	$sthMax->finish;
	$sthUpMax->finish;
	$sthUpAna->finish;

	$dbh->disconnect;

	####<Reload all frames>####
	&reloadFrames;
}


################################
####<<Reloading all frames>>####
################################
sub reloadFrames {
	sleep 3;

	print qq
|<SCRIPT LANGUAGE="JavaScript">
if (top.protWindow && !top.protWindow.closed) {top.protWindow.location.reload(true);} // in case changes in validation affect protein being displayed in protWindow
|;
	if ($ITEM eq 'ANALYSIS') {
		print "top.promsFrame.location=\"$promsPath{cgi}/startValidation.cgi?ID=$itemID\"\n";
	}
	else {
		print qq
|top.promsFrame.selectedAction='summary';
parent.optionFrame.selectOption();
|;
	}
	print qq
|</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

####>Revision history<####
# 1.6.7 Modif to select 'SEQUEST' as default search engine for Qualitative Peptide/Protein selection (MLP 16/04/18)
# 1.6.6 Minor modif on protein SCORE calculation (PROTEIN_VALIDATION) (MLP 04/05/17) 
# 1.6.5 Minor modif for TANDEM format (MLP 05/04/17)
# 1.6.4 Added modifications for TANDEM format & added &promsMod::cleanParameters verification (MLP 16/12/16)
# 1.6.3 Check on undefined $taxonomy (PP 21/01/15)
# 1.6.2 Minor Changes in default values display (PP 10/10/14)
# 1.6.1 Handles charge states > 4 for SEQUEST (PP 23/10/13)
# 1.6.0 Updating analysis history during lower-scorer peptide activation (FY 17/09/13)
# 1.5.9 Fix undef $ITEM in case of &clearValidations (PP 18/06/13)
# 1.5.8 Fix absolute-min-score miscalculation if analysis min score = 0 (FY 17/05/13)
# 1.5.7 Better server time-out control & utf-8 & LOWER_SCORES (PP 24/04/13)
# 1.5.6 Multi-databank search (PP 17/12/12)
# 1.5.5 Minor modifications to make it work with Paragon converted XML (GA 11/12/12)
# 1.5.4 Use 1 pep score for all decoys if PDM format (no matching decoy proteins in MSF files!) (PP 09/01/12)
# 1.5.3 Correction bug template not saved & form reset if back to no template (PP 23/08/11)
# 1.5.2 header(-'content-encoding'=>'no') & bug undef values after submission (PP 01/08/11)
# 1.5.1 Added search type selection & code cleaning (PP 09/06/11)
# 1.5.0 Compatibility with mixed (SQ) search (PP 27/05/11)
# 1.5.0 Minor changes on $fileFormatString & $selectTemplateString (PP 11/04/2011)
# 1.4.9 Adding template management (FY 02/2011)