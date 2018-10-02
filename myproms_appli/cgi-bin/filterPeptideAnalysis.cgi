#!/usr/local/bin/perl -w

################################################################################
# filterPeptideAnalysis.cgi          1.2.6                                     #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Select/reject or Flag 1 or more Analyses based on                            #
# peptides (not)shared by 1 or more other analyses                             #
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
use promsConfig;
use promsMod;
use strict;


#print header;warningsToBrowser(1);
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my %msTypeName=&promsConfig::getMsType;

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

###################
####>Arguments<####
###################
my $ITEM=uc(param('ITEM'));
my $item=lc($ITEM);
my $itemID=param('ID');
my $projectID=param('id_project');
my $maxRank= &promsConfig::getMaxRank;
my $allrank=10; # !!!** Only top ranking peptides are compared **!!!
my $userID = $ENV{'REMOTE_USER'};


#########################
####>Displaying form<####
#########################
if (!param('apply')) {

	############################
	####>Fetching item info<####
	############################
	my ($itemName)=$dbh->selectrow_array("SELECT NAME FROM $ITEM WHERE ID_$ITEM=$itemID");
	my $targetString='';
	my $specFrameComment='';
	my $brString='';
	my $analysisString='';
	my $titleString='';
	#my $subTitleString='';
	my (%listDataBank,@itemAnalyses);
	my $sthAnaDb=$dbh->prepare("SELECT D.ID_DATABANK,NAME,VERSION_NAME FROM ANALYSIS_DATABANK AD,DATABANK D WHERE ID_ANALYSIS=? AND AD.ID_DATABANK=D.ID_DATABANK ORDER BY DB_RANK");
	if ($ITEM eq 'ANALYSIS') {
		$targetString='target="spectrumFrame"';
		$titleString.='<SELECT class="title" style="text-align:left" onchange="goToQualitative();"><OPTION>Qualitative</OPTION><OPTION selected>Comparative</OPTION></SELECT>';

		my $sthSamp=$dbh->prepare("SELECT ID_ANALYSIS,VALID_STATUS,MS_TYPE,DATA_FILE,FILE_FORMAT,WIFF_FILE,TAXONOMY,0,NAME FROM ANALYSIS WHERE ID_SAMPLE=(SELECT ID_SAMPLE FROM ANALYSIS WHERE ID_ANALYSIS=$itemID) ORDER BY DISPLAY_POS ASC");
		$sthSamp->execute;
		while (my @anaData=$sthSamp->fetchrow_array) {
			my @dbUsed;
			$sthAnaDb->execute($anaData[0]);
			while (my($dbID,$dbName,$version)=$sthAnaDb->fetchrow_array) {
				push @dbUsed,$dbID;
				$listDataBank{$dbID}=$dbName; # if ($anaData[0]==$itemID);
				$listDataBank{$dbID}.=" ($version)" if $version;
			}
			$anaData[7]=\@dbUsed; # ref to list

			push @itemAnalyses,\@anaData;
		}
		$sthSamp->finish;
	}
	else {
		$titleString.='Comparative';
		#$subTitleString.='<BR><BR><FONT class="title3">(Apply to all MS/MS Ion Search Analyses in';
		#$subTitleString.=&promsMod::getItemType($ITEM)." <FONT color=\"#DD0000\">$itemName</FONT>)</FONT>\n";
		$brString='<BR>';
		$specFrameComment='//';
		$analysisString=' selected Analyses of';

		#####>Info from all dataBanks
		#my $sthDB = $dbh->prepare("SELECT ID_DATABANK,NAME FROM DATABANK");
		#$sthDB->execute;
		#while (my ($idDbk,$DbkName)= $sthDB->fetchrow_array) {
		#	$listDataBank{$idDbk}=$DbkName;
		#	#$listDataBank{$idDbk} .=" ($version)" if $version;
		#}
		#$sthDB->finish;

		my @sthItem;
		my $baseFieldString='ID_ANALYSIS,VALID_STATUS,MS_TYPE,DATA_FILE,FILE_FORMAT,WIFF_FILE,TAXONOMY,0';
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
		foreach my $sth (@sthItem) {
			$sth->execute;
			while (my @anaData=$sth->fetchrow_array) {
				my @dbUsed;
				$sthAnaDb->execute($anaData[0]);
				while (my($dbID,$dbName,$version)=$sthAnaDb->fetchrow_array) {
					push @dbUsed,$dbID;
					$listDataBank{$dbID}=$dbName;
					$listDataBank{$dbID}.=" ($version)" if $version;
				}
				$anaData[7]=\@dbUsed; # ref to list
				push @itemAnalyses,\@anaData;
			}
			$sth->finish;
		}
	}
	$sthAnaDb->finish;
	$titleString.=" Peptide/Protein Selection\n";

	my $selectTemplateString = "<SELECT name='templateList' onchange='useTemplate(this.selectedIndex)' class=\"title2\"><OPTION value=''>No template</OPTION>";
	my $sthGetTemplate = $dbh->prepare("SELECT NAME, PARAM_STRG, DES, IS_DEFAULT FROM VALIDATION_TEMPLATE WHERE ID_USER='$userID' AND VAL_TYPE='compa'");
	$sthGetTemplate->execute;
	my $refTemplateList = $sthGetTemplate->fetchall_arrayref;
	$sthGetTemplate->finish;
	$dbh->disconnect;

	my $defaultIndex = 0;
	my $i = 0;
	foreach my $template (@{$refTemplateList}){
		$i++;
		my $templateName = $template->[0];
		my $parameters = $template->[1];
		my $des = ($template->[2])? "<B>Description : </B>".&promsMod::HTMLcompatible($template->[2]) : "No description" ;
		$des =~ s/\/\//\\\/\\\//g; $des =~ s/&apos;*/\\\'/g; ### javascript compatible
		$selectTemplateString .= "<OPTION value=$parameters onmouseover=\"popup('$des')\" onmouseout=\"popout()\">$templateName</OPTION>\n";
		if($template->[3]){ $defaultIndex = $i}
	}
	$selectTemplateString .= "</SELECT>";

	##############
	####>HTML<####
	##############
	my ($lightColor,$darkColor)=&promsConfig::getRowColors;
	print header(-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
	TD {font-weight:bold;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
|;
	&promsMod::popupInfo();
	if ($ITEM eq 'ANALYSIS') {
		my $selMsType = param('MSTYPE');
		print qq
|function goToQualitative() {
	window.location='./autoSelect.cgi?ITEM=$ITEM&ID=$itemID&id_project=$projectID&ACT=select&MSTYPE=$selMsType';
}
|;
	}
	print qq
|function checkall() {
	var anaBox=document.filterForm.anaList;
	if (anaBox.length) { // more than 1 checkboxes
		for (i=0; i < anaBox.length; i++){
			anaBox[i].checked=checkStatus;
			updateSelect(i);
		}
	}
	else {anaBox.checked=checkStatus;} // Only 1 checkbox
	checkStatus=(checkStatus==true)? false : true;
}
function updateSelect(boxIndex) {
	var myForm=document.filterForm;
	var newStatus=(myForm.anaList[boxIndex].checked)? false : true;
	myForm.anaComp[boxIndex].disabled=newStatus;
}
function checkForm(action) {
	var myForm=document.filterForm;

	//Selected items
	var anaBox=myForm.anaList;
	var compSelect=myForm.anaComp;
	var itemChecked=0;
	if (anaBox.length) { // more than 1 checkboxes
		for (i=0; i < anaBox.length; i++){
			if (anaBox[i].checked==true) {
				itemChecked++;
				if (action == ' Proceed ') {
					if (!compSelect[i].value) {
						alert('Select Analysis to be used for comparison.');
						return;
					}
				}
				else {break;}
			}
		}
	}
	else {
		if (anaBox.checked==true) {
			itemChecked++;
			if (action == ' Proceed ' && !compSelect.value) {
				alert('Select Analysis to be used for comparison.');
				return;
			}
		}
	}
	if (itemChecked==0) {
		alert('Select at least 1 Analysis to be processed.');
		return;
	}

	if (action == ' Proceed ') {

		// Main action choice
		if (!myForm.actionPep[0].checked && !myForm.actionPep[1].checked) {
			alert("Choose either 'Select/Reject' or 'Flag'.");
			return;
		}
		// Select/Reject
		if (myForm.actionPep[0].checked) {
			if (!myForm.selectPep.checked && !myForm.rejectPep.checked) {
				alert("Select at least 1 action: 'Select' and/or 'Reject'.");
				return;
			}
			else if (!myForm.anaAction[0].checked && !myForm.anaAction[1].checked) {
				alert("Select which peptides to 'Select' and/or 'Reject'.");
				return;
			}
		}
		//Flag
		if (myForm.actionPep[1].checked) {
			if (!myForm.FlagUp.checked && !myForm.FlagDown.checked) {
				alert('Select at least 1 flag type.');
				return;
			}
			else if (!myForm.FltAction[0].checked && !myForm.FltAction[1].checked) {
				alert('Select which peptides to flag.');
				return;
			}
		}
		if(myForm.actTemplate.checked){
			if(myForm.templateName.value){ //checking template fields
				for (var i=0;i<myForm.templateList.length;i++){
					if (myForm.templateName.value.toLowerCase() == myForm.templateList[i].text.toLowerCase()){
						return confirm("Are you sure to want to overwrite "+myForm.templateList[i].text+" ?");
					}
				}
			}
			else {
				alert("Enter a template name.");
				return false;
			}
        }
	}
	myForm.apply.value=action;
	myForm.submit();
}
function cancelAction() {
	if ('$ITEM'=='ANALYSIS') {
		if (parent.selectedView=='myList') {parent.writeMenu(parent.selectedProtID,2);}
		else {parent.setView(parent.selectedView);}
	}
	else {
		//top.promsFrame.selectedAction='summary';
		top.promsFrame.optionFrame.selectOption(); // back to 'Process Analyses'
	}
}
function showTemplateBox(){
	var myForm=document.filterForm;
	myForm.templateName.style.visibility=(myForm.actTemplate.checked==false)?"hidden" : "visible";
}
function useTemplate(index){
	var myForm=document.filterForm;
	var paramString = myForm.templateList.options[index].value;
	if (!paramString){
		myForm.reset();
	}
	var paramList = paramString.split(';');
	for(var i=0; i<paramList.length; i++){
		var parameters = paramList[i].split(':');
		if (parameters[1] == 'check'){
			myForm[parameters[0]].checked = (parameters[2] == 1)? true : false;
		}
		else if (parameters[1] == 'text'){
			myForm[parameters[0]].value = parameters[2];
		}
		else if (parameters[1] == 'select'){
			myForm[parameters[0]].selectedIndex = parameters[2] - 1; //only for maxRank
		}
		else if (parameters[1] == 'radio'){
			for(var j=0;j<myForm[parameters[0]].length;j++){
				if(myForm[parameters[0]][j].value == parameters[2]){
					myForm[parameters[0]][j].checked = true;
				}
				else {
					myForm[parameters[0]][j].checked = false;
				}
			}
		}
	}
	myForm.templateName.value = (index > 0)? myForm.templateList.options[index].text : '';
}
var checkStatus=true;
$specFrameComment top.promsFrame.spectrumFrame.location="$promsPath{html}/nothing.html";
</SCRIPT>
</HEAD>
<BODY topmargin=0 background="$promsPath{images}/bgProMS.gif" onload="document.filterForm.templateList.selectedIndex = $defaultIndex; useTemplate($defaultIndex);">
<CENTER>
<FORM name="filterForm" method="post" $targetString>
$brString<FONT class="title">$titleString</FONT>$brString$brString
<INPUT type="hidden" name="ID" value="$itemID" />
<INPUT type="hidden" name="ITEM" value="$ITEM" />
<INPUT type="hidden" name="id_project" value="$projectID" />
<INPUT type="hidden" name="apply" value="" />
<BR><BR>

<TABLE border=0 cellpadding=0 cellspacing=0>
<TR>
<TD align='left'>
<FONT class="title2">Parameters template : $selectTemplateString</FONT>
</TD>
</TR>
<TR>
<TD nowrap valign=top>
<TABLE bgcolor=$darkColor border=0>
	<TR>
		<TH nowrap><FONT class="title2"><INPUT type="radio" name="actionPep" value="selRej">Select/Reject Peptides</FONT></TH>
	</TR>
	<TR bgcolor=$lightColor>
		<TD><FONT style="line-height:5px;"><BR></FONT>
			<INPUT type="checkbox" name="selectPep">Select peptides matching the option below:<BR>
			<INPUT type="checkbox" name="rejectPep">Reject peptides that <u>do not </u> match this option:<BR>
			<FONT style="line-height:5px;"><BR></FONT>
			&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;<INPUT type="radio" name="anaAction" value="inter">&nbsp;common to Reference and$analysisString <FONT color='#DD0000'>$itemName</FONT>.<BR>
			&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;<INPUT type="radio" name="anaAction" value="diff">&nbsp;unique to$analysisString <FONT color='#DD0000'>$itemName</FONT>.<BR>
			&nbsp;(Already selected/rejected peptides will not be affected)<BR>
			<FONT style="line-height:5px;"><BR></FONT>
		</TD>
	</TR>

	<TR>
		<TH nowrap>
			<FONT class="title2">or&nbsp<INPUT type="radio" name="actionPep" value="flag">Flag Peptides</FONT>
			&nbsp &nbsp &nbsp<INPUT type="button" value="Remove Flags" onclick="checkForm(this.value);"/>
		</TH>
	</TR>
	<TR bgcolor=$lightColor>
		<TD><FONT style="line-height:5px;"><BR></FONT>
			<INPUT type="checkbox" name="FlagUp">&nbsp;Flag (<IMG src="$promsPath{images}/lightYellow1.gif" hspace=0 border=0 height=11 width=11>) peptides matching the option below:<BR>
			<INPUT type="checkbox" name="FlagDown">&nbsp;Flag (<IMG src="$promsPath{images}/lightOrange1.gif" hspace=0 border=0 height=11 width=11>) peptides that <u>do not </u> match this option:<BR>
			<FONT style="line-height:5px;"><BR></FONT>
			&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;<INPUT type="radio" name="FltAction" value="inter">&nbsp;common to Reference and$analysisString <FONT color='#DD0000'>$itemName</FONT>.<BR>
			&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;<INPUT type="radio" name="FltAction" value="diff">&nbsp;unique to$analysisString <FONT color='#DD0000'>$itemName</FONT>.<BR>
			<FONT style="line-height:5px;"><BR></FONT>
		</TD>
	</TR>

	<TR>
		<TH nowrap><FONT class="title2">Options</FONT></TH>
	</TR>
	<TR bgcolor=$lightColor>
		<TD nowrap><FONT style="line-height:5px;"><BR></FONT>
			&nbsp;Compare only peptides with rank &le;<SELECT name="selMaxRank" style="font-weight:bold;">
|;
	foreach my $rank (1..$allrank) {print "<OPTION value=\"$rank\">$rank</OPTION>";}
	print qq
|</SELECT><BR>
			<INPUT type="checkbox" name="validPept">Use only peptides already selected in Reference.<BR>
			<!--<INPUT type="checkbox" name="rejPept" checked>Ignore already rejected peptides in$analysisString <FONT color='#DD0000'>$itemName</FONT>.<BR>-->
			<INPUT type="checkbox" name="rejMinValPept" checked>Ignore peptides in$analysisString <FONT color='#DD0000'>$itemName</FONT> with score below <INPUT type="text" name="smPepSc" value="0" size="3" >. <BR>
			<INPUT type="checkbox" name="ignRefMinScore" >Ignore peptides in Reference with score below <INPUT type="text" name="refMinScore" value="0" size="3" >. <BR>
			<INPUT type="checkbox" name="discardEltTime">Ignore peptides with elution time differing by more than
			<INPUT type="text" name="maxEltTime" value="2" size="2" > mn.<BR>
			<FONT style="line-height:5px;"><BR></FONT>
		</TD>
	</TR>
	<TR>
		<TH nowrap align='left'><INPUT type="checkbox" name="actTemplate" onclick="showTemplateBox();">Save parameters as template. &nbsp
		<INPUT type="text" name="templateName" style='visibility:hidden'>
		</TH>
	</TR>
	<TR>
		<TH nowrap><INPUT type="button" value=" Proceed " onclick="checkForm(this.value);"/>
			&nbsp &nbsp &nbsp<INPUT type="reset" value="  Clear  " />
			&nbsp &nbsp &nbsp<INPUT type="button" value="Cancel" onClick="cancelAction();">
		</TH>
	</TR>
</TABLE>
</TD>
<TD>&nbsp</TD>
<TD valign=top><TABLE border=0 cellspacing=0 cellpadding=1>
<TR bgcolor="$darkColor"><TH class="rbBorder"><IMG src="$promsPath{images}/good.gif" onclick="checkall()"></TH><TH class="rbBorder" colspan=2>&nbsp;<FONT class="title2">Select Analyses</FONT>&nbsp;</TH><TH class="bBorder">&nbsp;<FONT class="title2">Reference</FONT>&nbsp;</TH></TR>
|;
	my %itemIcones=&promsConfig::getItemIcones;
	my $bgColor=($ITEM eq 'SAMPLE' || $ITEM eq 'SPOT' || $ITEM eq 'ANALYSIS')? $lightColor : $darkColor;
	my %prevItemName;
	my $boxIndex=-1;
	foreach my $refAnaData (@itemAnalyses) {
		my ($anaID,$valStat,$msType,$dataFile,$fileFormat,$wiffFile,$taxonomy,$refDbUsed,@projHierarchy)=@{$refAnaData};
		next if ($ITEM eq 'ANALYSIS' && $anaID != $itemID);
		$wiffFile=~s/\n//; # brakes popup otherwise
		$taxonomy=~s/\(.*\)//;
		$msType=$msTypeName{$msType}; # global
		$fileFormat=~s/\..*//;
		##>Row color
		my $fatherIt=$projHierarchy[-3];
		if ($fatherIt && (!$prevItemName{$fatherIt} || $prevItemName{$fatherIt} ne $projHierarchy[-2])) { # keep color if same analysis parent item (SAMPLE or SPOT)
			$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
		}
		##>Checkbox
		my ($boxStr,$disabSelect) = ($ITEM eq 'ANALYSIS')? ("<IMG src=\"$promsPath{images}/good.gif\"><INPUT type=\"checkbox\" name=\"anaList\" value=\"$itemID\" style=\"display:none\" checked>",'') : ($valStat>1)? ('-',undef) : ("<INPUT type=\"checkbox\" name=\"anaList\" value=\"$anaID\" onclick=\"updateSelect(".++$boxIndex.")\">",' disabled');
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
|<TR valign=top bgcolor=$bgColor><TH valign=middle>$boxStr</TH><TH nowrap align=left valign=middle>$parentStrg</TH>
<TH nowrap align=left valign=middle><IMG src="$promsPath{images}/$itemIcones{$anaCode}">&nbsp;
<A href="javascript:void(null)" onmouseover="popup('<B>Type:</B> $msType<BR><B>Search file:</B> $dataFile<BR><B>MS file:</B> $wiffFile<BR><B>Search engine:</B> $fileFormat<BR><B>Taxonomy:</B> $taxonomy$dbString')" onmouseout="popout()">$projHierarchy[-1]</A>&nbsp;
</TH><TD>
|;
		if ($valStat<2) {
			print "<SELECT name=\"anaComp\" style=\"width:150px\"$disabSelect><OPTION value=\"\">-= Select=-</OPTION>\n";
			my $prevParNames='';
			foreach my $refCompAna (@itemAnalyses) {
				my ($compID,$vStat,$type,$dFile,$fFormat,$wFile,$tax,$dID,@hierarchy)=@{$refCompAna};
				next if $compID==$anaID;
				my $parNames='';
				for (my $i=0;$i<=$#hierarchy-2;$i+=2) { # stops before ana name
					$parNames.=' > ' if $parNames;
					$parNames.=$hierarchy[$i+1];
				}
				if ($parNames ne $prevParNames) {
					print "</OPTGROUP>\n" if $prevParNames;
					print "<OPTGROUP label=\"$parNames\">\n";
					$prevParNames=$parNames;
				}
				print "<OPTION value=\"$compID\">$hierarchy[-1]</OPTION>\n";
			}
			print "</OPTGROUP>\n" if $prevParNames;
			print "</SELECT>\n";
		}
		print "</TR>\n";
	}
	print qq
|</TABLE></TD></TR>
</TABLE>
</FORM>
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


#############################################################
###################>Processing form<#########################
#############################################################
print header(-charset=>'utf-8'); warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<BR>
<FONT class="title3">
|;
my $apply = param('apply');

####>Recovering list of Analysis<####
my @listAnalyses = param('anaList'); #list of analysis to be validated

	# if ($item eq 'analysis') {
		# @listAnalyses=($itemID);
	# }
	# else {
	# my $selAnaQuery;
	# if ($item eq 'sample') {
		# $selAnaQuery="SELECT ID_ANALYSIS FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND VALID_STATUS<=1 AND SAMPLE.ID_SAMPLE=$itemID";
	# }
	# elsif ($item eq 'experiment') {
		# $selAnaQuery="SELECT ID_ANALYSIS FROM ANALYSIS,SAMPLE,EXPERIMENT WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_EXPERIMENT=EXPERIMENT.ID_EXPERIMENT AND VALID_STATUS<=1 AND EXPERIMENT.ID_EXPERIMENT=$itemID";
	# }
	# my $sthSelAna=$dbh->prepare($selAnaQuery);
	# $sthSelAna->execute();
	# while(my($anaID)=$sthSelAna->fetchrow_array) {
		# push @listAnalyses,$anaID;
	# }
	# $sthSelAna->finish;
# }


#################################################
                ####>APPLY<####
#################################################
if ($apply ne "Remove Flags") {

	####>Recovering other parameters<####
	my @anaComp=param('anaComp');
#print 'CHK: ',join(',',@listAnalyses),"<BR>\n";
#print 'SEL: ',join(',',@anaComp),"<BR>\n";
#exit;
	my $selMaxRank=param('selMaxRank');
	my $actionPep = param ('actionPep');

	#Validate/Reject
	my $anaAction = param('anaAction');
	my $selectPep = (param('selectPep'))? 1 : 0;
	my $rejectPep = (param('rejectPep'))? 1 : 0;

	#Select/Hide
	my $fltAction = param('FltAction');
	my $FlagUp = (param('FlagUp'))? 1 : 0;
	my $FlagDown = (param('FlagDown'))? 1 : 0;

	#Options
	my $validPept = (param('validPept'))? 1 : 0; #ok
	my $minValidStatus = ($validPept==1)? 1 : -1;
	#my $rejPept = (param('rejPept'))? 1 : 0;
	my $rejMinValPept = (param('rejMinValPept'))? 1 : 0;
	my $smPepSc = param('smPepSc');
	$smPepSc=0 unless $smPepSc;
	my $discardEltTime = (param('discardEltTime'))? 1 : 0;
	my $maxEltTime = abs(param('maxEltTime'));
	$maxEltTime = 2 unless $maxEltTime ;
	my $ignRefMinScore= (param('ignRefMinScore'))? 1 : 0;
	my $refMinScore = param('refMinScore');
	$refMinScore=0 unless $refMinScore;

	#one2one List
	my %one2oneItems;
	for (my $i=0; $i<=$#listAnalyses; $i++) {
		$one2oneItems{$listAnalyses[$i]}=$anaComp[$i];
	}

	my $dbQuery;
	my %sthUpRk;
	foreach my $rank (1..$allrank) {
		$dbQuery.=",INFO_PEP$rank";
		$sthUpRk{$rank}=$dbh->prepare("UPDATE QUERY_VALIDATION SET INFO_PEP$rank=?,VALID_STATUS=? WHERE ID_QUERY=?");
	}
	my $sthQueryAna=$dbh->prepare("SELECT ID_QUERY,QUERY_NUM,VALID_STATUS,MASS_DATA,ELUTION_TIME$dbQuery FROM QUERY_VALIDATION WHERE ID_ANALYSIS=? AND VALID_STATUS>=?");
	my $sthPepAna=$dbh->prepare("SELECT ID_PEPTIDE,PEP_SEQ,MR_EXP,MR_OBS,ELUTION_TIME,SCORE FROM PEPTIDE WHERE ID_ANALYSIS=?");
	my $sthAnaNV=$dbh->prepare("SELECT NAME,VALID_STATUS FROM ANALYSIS WHERE ID_ANALYSIS=?");


	foreach my $analysisID (@listAnalyses) {

		my ($anaName,$validStatus);

		#############################################################################
		####>Recovering information from Reference analysis (compare one To one)<####
		#############################################################################
		my %peptideList;
		my  $anaID =  $one2oneItems{$analysisID};
		$sthAnaNV->execute($anaID);
		($anaName,$validStatus) = $sthAnaNV->fetchrow_array;
		print "Recovering information from Reference analysis $anaName...";

		if ($validStatus<1) { # not-validated analysis
			my $count=0;
			$sthQueryAna->execute($anaID,$minValidStatus);
			while (my ($queryID,$queryNum,$qValStatus,$massData,$elutionTime,@infoQuery) = $sthQueryAna->fetchrow_array) {
				if ($elutionTime && $elutionTime =~ /(\d+\.?\d*)\sto\s(\d+\.?\d*)/) {
					$elutionTime = ($1+$2)/2;
				}
				my $trueRank=1;
				my $prevScore;
				my $rank=0;
				foreach my $pepString (@infoQuery) {
					last unless $pepString;
					$rank++;
					my ($score) = ($pepString =~/SC=(-?\d+\.?\d*)/);
					if ($prevScore) { # peptides with same score are considered to be of same true rank !!!!!!!!!!
						$trueRank++ if $score<$prevScore;
						last if $trueRank>$selMaxRank;
					}
					$prevScore=$score;
					last if ($ignRefMinScore==1 && $score<$refMinScore);
					my ($sequence) = ($pepString =~/SEQ=(\w+)/);
					next unless $sequence;
					my ($sel) = ($pepString =~/SEL=(-?\d)/);
					next if $sel<=-1; # skip already rejected & lower scoring
					next if ($validPept==1 && $sel<1);
					#my ($varMod)=($pepString=~/VMOD=([^,]+)/);
					my $varMod=&promsMod::toStringVariableModifications($dbh,'QUERY',$queryID,$anaID,$sequence,$rank);
					$varMod='' unless $varMod;
					$varMod=~ s/Phospho \(\w+:/Phospho \(/ ; #compare Mascot/phenyx
					#$peptideList{"$sequence:$varMod"}{$queryID}{"selStatus"} = $sel;
					$peptideList{"$sequence:$varMod"}{$queryID}{"ElutionTime"} = ($elutionTime =~ /et(\d+\.\d+);/) ? $1 : $elutionTime;
					($peptideList{"$sequence:$varMod"}{$queryID}{"massExp"}) = ($massData=~/EXP=(\d+\.?\d*)/);
					($peptideList{"$sequence:$varMod"}{$queryID}{"massObs"}) = ($massData=~/OBS=(\d+\.?\d*)/);
				}
				$count++;
				if ($count == 100) {$count=0; print '.';}  #avoid loosing connection
			}
			#$sthQueryAna->finish;
		}
		elsif ($validStatus>=1) { #(partially) validated analysis
			$sthPepAna->execute($anaID);
			while (my($idPeptide,$pepSeq,$massExp,$massObs,$elutionTime,$score)=$sthPepAna->fetchrow_array) {
				my $varMod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$idPeptide,$anaID,$pepSeq);
				$varMod='' unless $varMod;
				$varMod=~ s/Phospho \(\w+:/Phospho \(/ ; #compare Mascot/phenyx
				next if ($ignRefMinScore==1 && $score<$refMinScore);
				if ($elutionTime && $elutionTime =~ /(\d+\.?\d*)\sto\s(\d+\.?\d*)/) {
						$elutionTime = ($1+$2)/2;
				}
				#$peptideList{"$pepSeq:$varMod"}{$idPeptide}{"selStatus"} = 3;
				$peptideList{"$pepSeq:$varMod"}{$idPeptide}{"ElutionTime"} = ($elutionTime =~ /et(\d+\.\d+);/) ? $1 : $elutionTime;
				$peptideList{"$pepSeq:$varMod"}{$idPeptide}{"massExp"} = $massExp;
				$peptideList{"$pepSeq:$varMod"}{$idPeptide}{"massObs"} = $massObs;
			}
			#$sthPepAna->finish;
		}
		print " Done.<BR>\n";

		######################################################
		####>Recovering information from current analysis<####
		######################################################
		my (%infoPepList,%queryList,%protList);
		$sthAnaNV->execute($analysisID);
		($anaName,$validStatus) = $sthAnaNV->fetchrow_array;
		print "Recovering data from analysis $anaName...";
		my $count;
		$sthQueryAna->execute($analysisID,-1);
		while (my ($queryID,$queryNum,$qValStatus,$massData,$elutionTime,@infoQuery) = $sthQueryAna->fetchrow_array) {
			$queryList{$queryID}{"validStatus"} = $qValStatus;
			if ($elutionTime && $elutionTime =~ /(\d+\.?\d*)\sto\s(\d+\.?\d*)/) {
				$elutionTime = ($1+$2)/2;
			}
			my $trueRank=1;
			my $prevScore;
			my $rank=0;
			foreach my $pepString (@infoQuery) {
				$rank++;
				#my $pepString = $infoQuery[$rank-1];
				last unless $pepString;
				my ($score) = ($pepString =~/SC=(-?\d+\.?\d*)/);
				if ($prevScore) { # rank>1; peptides with same score are considered to be of same true rank !!!!!!!!!!
					$trueRank++ if $score<$prevScore;
					last if $trueRank>$selMaxRank;
				}
				$prevScore=$score;
				last if ($rejMinValPept==1 && $score<$smPepSc);
				my ($sequence) = ($pepString=~/SEQ=(\w+)/);
				next unless $sequence;
				my ($sel) = ($pepString=~/SEL=(-?\d)/);
				next unless $sel==0; # ignore already selected/rejected peptides!!!!!!!!
				#my ($varMod)=($pepString=~/VMOD=([^,]+)/);
				my $varMod=&promsMod::toStringVariableModifications($dbh,'QUERY',$queryID,$analysisID,$sequence,$rank);
				$varMod='' unless $varMod;
				$varMod=~ s/Phospho \(\w+:/Phospho \(/ ; #compare Mascot/phenyx
				$infoPepList{"$sequence:$varMod"}{$queryID}{"SEL"} = $sel;
				($infoPepList{"$sequence:$varMod"}{$queryID}{"massExp"}) = ($massData=~/EXP=(\d+\.?\d*)/);
				($infoPepList{"$sequence:$varMod"}{$queryID}{"massObs"}) = ($massData=~/OBS=(\d+\.?\d*)/);
				#$infoPepList{"$sequence:$varMod"}{"query"} = $queryID;
				$infoPepList{"$sequence:$varMod"}{$queryID}{"rank"} = $rank;
				$infoPepList{"$sequence:$varMod"}{$queryID}{"pepString"} = $pepString;
				$infoPepList{"$sequence:$varMod"}{$queryID}{"queryNum"} = $queryNum;
				$infoPepList{"$sequence:$varMod"}{$queryID}{"ElutionTime"} = ($elutionTime =~ /[er]t(\d+\.\d+)/) ? $1 : $elutionTime;
				$infoPepList{"$sequence:$varMod"}{$queryID}{"ElutionTime"} = '' unless($infoPepList{"$sequence:$varMod"}{$queryID}{"ElutionTime"} !~ /s[cp]/);
				$infoPepList{"$sequence:$varMod"}{$queryID}{"score"} = $score;
			}
			$count++;
			if ($count == 100) {$count=0; print '.';}  #avoid loosing connection
		}
		#$sthAna->finish;
		print " Done.<BR>\n";

		############################
		####>Peptide comparison<####
		############################
		print "Comparing peptides...";
		foreach my $peptide (keys %infoPepList) {
			(my $truePeptide=$peptide)=~s/^_//;

			foreach my $queryID (keys %{$infoPepList{$peptide}}) { # selected analysis
				my $pepMatched=0;

				PEP:foreach my $pepVersion ($truePeptide,"_$truePeptide") { # normal & decoy versions

					if (defined (%{$peptideList{$pepVersion}})) { # reference analysis
						foreach my $ID (keys %{$peptideList{$pepVersion}}) {
							#$pepMatched = 1 ;
							my $thisMatch=1;
							$thisMatch=0 unless abs($peptideList{$pepVersion}{$ID}{"massExp"}-$infoPepList{$peptide}{$queryID}{"massExp"})<=1;
							$thisMatch=0 unless abs($peptideList{$pepVersion}{$ID}{"massObs"}-$infoPepList{$peptide}{$queryID}{"massObs"})<=1;
							$thisMatch=0 unless ($discardEltTime==0 || ($infoPepList{$peptide}{$queryID}{"ElutionTime"} && $peptideList{$pepVersion}{$ID}{"ElutionTime"} && abs($peptideList{$pepVersion}{$ID}{"ElutionTime"}-$infoPepList{$peptide}{$queryID}{"ElutionTime"})<=$maxEltTime)); #ElutionTime filter
							if ($thisMatch==1) { #1 match is found !
								$pepMatched=1;
								last PEP;
							}
						}
					}
					else {
						#$pepMatched=0;
						delete $peptideList{$pepVersion}; # just to be safe
					}

				}

				if ($actionPep eq "selRej") { # Select or reject peptides
					if ($anaAction eq "diff") { # Process peptides unique to selected analysis
						if ($selectPep == 1 && $pepMatched == 0) {
							$infoPepList{$peptide}{$queryID}{"updatedSEL"} = 1;
							$infoPepList{$peptide}{$queryID}{"pepString"}=~s/SEL=0/SEL=1/;
							$queryList{$queryID}{"validStatus"} = ($queryList{$queryID}{"validStatus"} < 1)? 1 : $queryList{$queryID}{"validStatus"}++;
						}
						if ($rejectPep == 1 && $pepMatched == 1) {
							$infoPepList{$peptide}{$queryID}{"updatedSEL"} = -2;
							$infoPepList{$peptide}{$queryID}{"pepString"} =~s/SEL=0/SEL=-2/;
							$queryList{$queryID}{"validStatus"} = ($queryList{$queryID}{"validStatus"} < 1)? -2 : $queryList{$queryID}{"validStatus"}--; # 0,-1,-2
						}
					}
					elsif ($anaAction eq "inter") { # Process peptides common to selected analysis & ref Set
						if ($selectPep == 1 && $pepMatched == 1) {
							$infoPepList{$peptide}{$queryID}{"updatedSEL"} = 1;
							$infoPepList{$peptide}{$queryID}{"pepString"} =~s/SEL=0/SEL=1/;
							$queryList{$queryID}{"validStatus"} = ($queryList{$queryID}{"validStatus"} < 1)? 1 : $queryList{$queryID}{"validStatus"}++;
						}
						if ($rejectPep == 1 && $pepMatched == 0) {
							$infoPepList{$peptide}{$queryID}{"updatedSEL"} = -2;
							$infoPepList{$peptide}{$queryID}{"pepString"} =~s/SEL=0/SEL=-2/;
							$queryList{$queryID}{"validStatus"} = ($queryList{$queryID}{"validStatus"} < 1)? -2 : $queryList{$queryID}{"validStatus"}--; # 0,-1,-2
						}
					}
				}
				elsif ($actionPep eq "flag") {
					if ($fltAction eq "diff") {
						if ($FlagUp == 1 && $pepMatched == 0) {
							$infoPepList{$peptide}{$queryID}{"flt"} = 1;
							if ($infoPepList{$peptide}{$queryID}{"pepString"} =~/FLT=/) {
								$infoPepList{$peptide}{$queryID}{"pepString"} =~s/FLT=-?\d/FLT=1/;
							}
							else {
								$infoPepList{$peptide}{$queryID}{"pepString"} =~s/SEL=/FLT=1,SEL=/;
							}
						}
						elsif ($FlagDown == 1 && $pepMatched == 1) {
							$infoPepList{$peptide}{$queryID}{"flt"} = -1;
							if ($infoPepList{$peptide}{$queryID}{"pepString"} =~/FLT=/) {
								$infoPepList{$peptide}{$queryID}{"pepString"} =~s/FLT=-?\d/FLT=-1/;
							}
							else {
								$infoPepList{$peptide}{$queryID}{"pepString"} =~s/SEL=/FLT=-1,SEL=/;
							}
						}
					}

					elsif ($fltAction eq "inter") {
						if ($FlagUp == 1 && $pepMatched == 1) {
							$infoPepList{$peptide}{$queryID}{"flt"} = 1;
							if ($infoPepList{$peptide}{$queryID}{"pepString"} =~/FLT=/) {
								$infoPepList{$peptide}{$queryID}{"pepString"} =~s/FLT=-?\d/FLT=1/;
							}
							else {
								$infoPepList{$peptide}{$queryID}{"pepString"} =~s/SEL=/FLT=1,SEL=/;
							}
						}
						elsif ($FlagDown == 1 && $pepMatched == 0) {
							$infoPepList{$peptide}{$queryID}{"flt"} = -1;
							if ($infoPepList{$peptide}{$queryID}{"pepString"} =~/FLT=/) {
								$infoPepList{$peptide}{$queryID}{"pepString"} =~s/FLT=-?\d/FLT=-1/;
							}
							else {
								$infoPepList{$peptide}{$queryID}{"pepString"} =~s/SEL=/FLT=-1,SEL=/;
							}
						}
					}
				}
			}
		}
		print " Done.<BR>\n";

		##################################
		####>Computing protein status<####
		##################################
		if ($actionPep eq "selRej") {
			print "Updating protein data...";
			my $sthAllProt=$dbh->prepare("SELECT ID_PROT_VALID,IDENTIFIER,SEL_STATUS,NUM_MATCH,MAX_MATCH,SCORE FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$analysisID");
			my $sthAllMatches=$dbh->prepare("SELECT IDENTIFIER,QUERY_NUM,PEP_RANK,MATCH_MULTI FROM RANK_PROTEIN_MATCH WHERE ID_ANALYSIS=$analysisID");
			my (%allPepProt,%ident2ProtID,%protPeptides);
			$sthAllProt->execute;
			while (my ($protID,$identifier,$selStatus,$numMatch,$maxMatch,$score) = $sthAllProt->fetchrow_array) {
				$protList{$protID}{"selStatus"} = $selStatus;
				$protList{$protID}{"numMatch"} = $numMatch;
				$protList{$protID}{"maxMatch"} = $maxMatch;
				$protList{$protID}{"score"} = $score;
				$protList{$protID}{"badMatch"} = 0;
				$ident2ProtID{$identifier}=$protID;
			}
			$sthAllProt->finish;
			$sthAllMatches->execute;
			while (my ($identifier,$queryNum,$rank,$matchMulti) = $sthAllMatches->fetchrow_array) {
				$protList{$ident2ProtID{$identifier}}{"matchMulti"}{"$queryNum:$rank"} = $matchMulti;
				$protPeptides{$ident2ProtID{$identifier}}{"$queryNum:$rank"}=1;
				push @{$allPepProt{"$queryNum:$rank"}},$ident2ProtID{$identifier};
			}
			$sthAllMatches->finish;

			my %checkProtStatus;
			foreach my $peptide (keys %infoPepList) {
				foreach my $queryID (keys %{$infoPepList{$peptide}}) {
					next unless defined($infoPepList{$peptide}{$queryID}{"updatedSEL"});
					#$sthProtVal->execute($infoPepList{$peptide}{$queryID}{"queryNum"},$infoPepList{$peptide}{$queryID}{"rank"});
					#while (my ($protID,$selStatus,$numMatch,$maxMatch,$matchMulti,$score) = $sthProtVal->fetchrow_array) { #}
						#if (!defined(%{$protList{$protID}})) {
						#	$protList{$protID}{"selStatus"} = $selStatus;
						#	$protList{$protID}{"numMatch"} = $numMatch;
						#	$protList{$protID}{"maxMatch"} = $maxMatch;
						#	$protList{$protID}{"score"} = $score;
						#	$protList{$protID}{"badMatch"} = 0;
						#}
					foreach my $protID (@{$allPepProt{"$infoPepList{$peptide}{$queryID}{queryNum}:$infoPepList{$peptide}{$queryID}{rank}"}}) { # @{allPepProt{queryNum:rank}}
						$protPeptides{$protID}{"$infoPepList{$peptide}{$queryID}{queryNum}:$infoPepList{$peptide}{$queryID}{rank}"}=0;
						if ($infoPepList{$peptide}{$queryID}{"updatedSEL"}==1) { # && $infoPepList{$peptide}{$queryID}{"SEL"}==0 <-- always true
							$protList{$protID}{"score"} += ($infoPepList{$peptide}{$queryID}{"score"}*$protList{$protID}{"matchMulti"}{"$infoPepList{$peptide}{$queryID}{queryNum}:$infoPepList{$peptide}{$queryID}{rank}"});
							$protList{$protID}{"numMatch"}++;
							$protList{$protID}{"selStatus"}=(($protList{$protID}{"numMatch"}+$protList{$protID}{"badMatch"})==$protList{$protID}{"maxMatch"})? 2 : 1;
						}
						elsif ($infoPepList{$peptide}{$queryID}{"updatedSEL"}==-2) { # && $infoPepList{$peptide}{$queryID}{"SEL"}==0 <-- always true
							$protList{$protID}{"badMatch"}++;
							$protList{$protID}{"selStatus"}=($protList{$protID}{"numMatch"}>=1 && (($protList{$protID}{"numMatch"}+$protList{$protID}{"badMatch"})==$protList{$protID}{"maxMatch"}))? 2 : ($protList{$protID}{"numMatch"}>0)? 1 : ($protList{$protID}{"maxMatch"}==$protList{$protID}{"badMatch"})? 0 : -1;
						}
						if ($protList{$protID}{"selStatus"}==1) {$checkProtStatus{$protID}=1;} elsif ($checkProtStatus{$protID}) {delete $checkProtStatus{$protID};}
						#elsif ($infoPepList{$peptide}{$queryID}{"updatedSEL"}==-2 && $infoPepList{$peptide}{$queryID}{"SEL"}==1) {
						#	$protList{$protID}{"score"} -= ($infoPepList{$peptide}{"score"}*$matchMulti);
						#	$protList{$protID}{"numMatch"}--;
						#	$protList{$protID}{"badMatch"}++;
						#	$protList{$protID}{"selStatus"}=($protList{$protID}{"numMatch"}>0)? 1 : ($protList{$protID}{"maxMatch"}==$protList{$protID}{"badMatch"})? 0 : -1;
						#}
						#elsif ($infoPepList{$peptide}{$queryID}{"updatedSEL"}==-2 && $infoPepList{$peptide}{$queryID}{"SEL"}==0) {
						#	$sthProtVal->execute($analysisID,$infoPepList{$peptide}{$queryID}{"queryNum"},$infoPepList{$peptide}{$queryID}{"rank"});
						#	while (my ($protID,$selStatus,$numMatch,$maxMatch,$score) = $sthProtVal->fetchrow_array) {
						#		if (!defined(%{$protList{$protID}})) {
						#				$protList{$protID}{"selStatus"} = $selStatus;
						#				$protList{$protID}{"numMatch"} = $numMatch;
						#				$protList{$protID}{"maxMatch"} = $maxMatch;
						#				$protList{$protID}{"score"} = $score;
						#				$protList{$protID}{"badMatch"} = 0;
						#		}
						#		$protList{$protID}{"badMatch"}++;
						#		$protList{$protID}{"selStatus"}=($protList{$protID}{"numMatch"}>=1 && (($protList{$protID}{"numMatch"}+$protList{$protID}{"badMatch"})==$protList{$protID}{"maxMatch"}))? 2 : ($protList{$protID}{"numMatch"}>0)? 1 : ($protList{$protID}{"maxMatch"}==$protList{$protID}{"badMatch"})? 0 : -1;
						#	}
						#}
						#elsif ($infoPepList{$peptide}{"updatedSEL"}==-2 && $infoPepList{$peptide}{"SEL"}==1) {
						#	$sthProtVal->execute($analysisID,$infoPepList{$peptide}{"queryNum"},$infoPepList{$peptide}{"rank"});
						#	while (my ($protID,$selStatus,$numMatch,$maxMatch,$score) = $sthProtVal->fetchrow_array) {
						#		if (!defined( %{$protList{$protID}})) {
						#			$protList{$protID}{"selStatus"} = $selStatus ;
						#			$protList{$protID}{"numMatch"} = $numMatch ;
						#			$protList{$protID}{"maxMatch"} = $maxMatch ;
						#			$protList{$protID}{"score"} = $score ;
						#			$protList{$protID}{"badMatch"} = 0;
						#		}
						#		$protList{$protID}{"score"} -= $infoPepList{$peptide}{"score"};
						#		$protList{$protID}{"numMatch"}--;
						#		$protList{$protID}{"badMatch"}++;
						#		$protList{$protID}{"selStatus"}=($protList{$protID}{"numMatch"}>0)? 1 : ($protList{$protID}{"maxMatch"}==$protList{$protID}{"badMatch"})? 0 : -1;
						#	}
						#}
					}
				}
			}

			###>Checking for previously rejected matches (for prot status=1)<### (not counted in test above)
			if (scalar keys %checkProtStatus) {
				my %sthExtPep;
				foreach my $rank (1..$allrank) {
					$sthExtPep{$rank}=$dbh->prepare("SELECT INFO_PEP$rank FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID AND QUERY_NUM=?");
				}
				foreach my $protID (keys %checkProtStatus) {
					foreach my $queryRank (keys %{$protPeptides{$protID}}) {
						next if $protPeptides{$protID}{$queryRank}==0; # pep data already retrieved
						my ($queryNum,$rank)=split(':',$queryRank);
						$sthExtPep{$rank}->execute($queryNum);
						my ($pepString)=$sthExtPep{$rank}->fetchrow_array;
						$protList{$protID}{"badMatch"}++ if $pepString=~/SEL=-[23]/;
					}
					$protList{$protID}{"selStatus"}=2 if $protList{$protID}{"numMatch"}+$protList{$protID}{"badMatch"}==$protList{$protID}{"maxMatch"};
				}
				foreach my $rank (keys %sthExtPep) {$sthExtPep{$rank}->finish;}
			}
			print " Done.<BR>\n";
		}


		####>Storing data in Database<####
		print "Updating Database...";

		###>Updating query_Validation<###
		#my $sthUpQuery1=$dbh->prepare("UPDATE QUERY_VALIDATION SET INFO_PEP1=?,VALID_STATUS=? WHERE ID_QUERY=?");
		#my $sthUpQuery2=$dbh->prepare("UPDATE QUERY_VALIDATION SET INFO_PEP2=?,VALID_STATUS=? WHERE ID_QUERY=?");
		#my $sthUpQuery3=$dbh->prepare("UPDATE QUERY_VALIDATION SET INFO_PEP3=?,VALID_STATUS=? WHERE ID_QUERY=?");
		foreach my $peptide (keys %infoPepList) {
			foreach my $queryID (keys %{$infoPepList{$peptide}}) {
				next unless (defined($infoPepList{$peptide}{$queryID}{"updatedSEL"}) || defined($infoPepList{$peptide}{$queryID}{"flt"}));
				#my $rank = $infoPepList{$peptide}{$queryID}{"rank"};
				#my $pepString =$infoPepList{$peptide}{$queryID}{"pepString"};
				#my $queryID = $infoPepList{$peptide}{"query"} ;
				#my $validStatus = $queryList{$queryID}{"validStatus"};
				$sthUpRk{$infoPepList{$peptide}{$queryID}{"rank"}}->execute($infoPepList{$peptide}{$queryID}{"pepString"},$queryList{$queryID}{"validStatus"},$queryID) || die $dbh->errstr;
			}
		}
		print "...";

		if ($actionPep eq "selRej") {
			#updating protein_Validation
			#my $sthUpProt = $dbh->prepare("UPDATE PROTEIN_VALIDATION SET SEL_STATUS=?,NUM_MATCH=?,SCORE=?,CONF_LEVEL=2 WHERE ID_PROT_VALID=?");
			my $sthUpProt=$dbh->prepare("UPDATE PROTEIN_VALIDATION SET SEL_STATUS=?,NUM_MATCH=?,SCORE=?,CONF_LEVEL=2 WHERE ID_PROT_VALID=?");
			foreach my $protID (keys %protList) {
				$sthUpProt->execute($protList{$protID}{"selStatus"},$protList{$protID}{"numMatch"},$protList{$protID}{"score"},$protID) || die $dbh->errstr;
			}
			$sthUpProt->finish;
		}
#$dbh->rollback;
		$dbh->commit; # for each analysis
		print " Done.<BR><BR>\n";
	}
	$sthQueryAna->finish;
	$sthPepAna->finish;
	$sthAnaNV->finish;
	foreach my $rank (keys %sthUpRk) {$sthUpRk{$rank}->finish;}
	print "All selected analyses have been validated!<BR>\n";

	### Saving Parameters ###
	my $paramStrg = "selMaxRank:select:$selMaxRank;actionPep:radio:$actionPep;anaAction:radio:$anaAction;selectPep:check:$selectPep;rejectPep:check:$rejectPep;FltAction:radio:$fltAction;FlagUp:check:$FlagUp;FlagDown:check:$FlagDown;validPept:check:$validPept;rejMinValPept:check:$rejMinValPept;smPepSc:text:$smPepSc;discardEltTime:check:$discardEltTime;maxEltTime:text:$maxEltTime;ignRefMinScore:check:$ignRefMinScore;refMinScore:text:$refMinScore;" ;
		## Saving template ##
	if(param('actTemplate')){
		print "<BR> Saving template...\n";
		my $name=param('templateName');


		my $existIDTemplate = $dbh->selectrow_array("SELECT ID_VAL_TEMPLATE FROM VALIDATION_TEMPLATE WHERE NAME='$name' AND ID_USER='$userID'");
		my $sthAddTemplate;
		if($existIDTemplate){
			$sthAddTemplate = $dbh->prepare("UPDATE VALIDATION_TEMPLATE SET PARAM_STRG=? WHERE ID_VAL_TEMPLATE=$existIDTemplate");
			$sthAddTemplate->execute($paramStrg);
		} else {
			$sthAddTemplate = $dbh->prepare("INSERT INTO VALIDATION_TEMPLATE(NAME,PARAM_STRG,VAL_TYPE,ID_USER) VALUES (?,?,?,?)");
			$sthAddTemplate->execute($name,$paramStrg,'compa',$userID);
		}
		$sthAddTemplate->finish;
		$dbh->commit;
		print "<BR> Template $name saved !<BR>\n";
	}
		## Saving history ##
	my $valType = ($actionPep eq 'selRej')? "comp_s" : ($actionPep eq 'flag')? "comp_f" : '' ;
	foreach my $analysisID (@listAnalyses){
		my $IDcomp = $one2oneItems{$analysisID};
		my $fileComp = $dbh->selectrow_array("SELECT DATA_FILE FROM ANALYSIS WHERE ID_ANALYSIS=$IDcomp");
		$paramStrg .= "COMP:$fileComp;";
		&promsMod::updateAnalysisHistory($dbh,$analysisID,$paramStrg,$valType);
	}
}


##########################################################
                ####>Removing Flag<####
##########################################################
elsif ($apply eq "Remove Flags") {
	my $dbQuery;
	my %sthUpRk;
	foreach my $rank (1..$allrank) {
		$dbQuery.=",INFO_PEP$rank";
		$sthUpRk{$rank}=$dbh->prepare("UPDATE QUERY_VALIDATION SET INFO_PEP$rank=? WHERE ID_QUERY=?");
	}
	my $sthQuery=$dbh->prepare("SELECT ID_QUERY$dbQuery FROM QUERY_VALIDATION WHERE ID_ANALYSIS=?");
	my $sthAnaN=$dbh->prepare("SELECT NAME FROM ANALYSIS WHERE ID_ANALYSIS=?");

	foreach my $analysisID (@listAnalyses) {
		my %localPepList;

		###>Recovering information from current analysis<###
		$sthAnaN->execute($analysisID);
		my ($anaName) = $sthAnaN->fetchrow_array;
		print "Recovering data from analysis $anaName...";
		$sthQuery->execute($analysisID);
		my $i;
		while (my ($queryID,@infoQuery) = $sthQuery->fetchrow_array) {
			my $rank=0;
			foreach my $pepString (@infoQuery) {
				$rank++;
				last unless $pepString;
				next unless $pepString =~/FLT=/;
				$pepString =~s/FLT=-?\d,?//;
				$localPepList{$queryID}{$rank}=$pepString;
			}
			$i++;
			if ($i == 100) {$i=0; print ".";}  #avoid loosing connection
		}
		print " Done.<BR>\n";

		###>Updating QUERY_VALIDATION<###
		print "Updating Database...";
		foreach my $queryID (keys %localPepList) {
			foreach my $rank (keys (%{$localPepList{$queryID}})) {
				$sthUpRk{$rank}->execute($localPepList{$queryID}{$rank},$queryID) || die $dbh->errstr;
			}
		}
		print " Done.<BR><BR>\n";
		$dbh->commit;
		&promsMod::updateAnalysisHistory($dbh,$analysisID,'','r_flag');
	}

	$sthQuery->finish;
	$sthAnaN->finish;
	foreach my $rank (keys %sthUpRk) {$sthUpRk{$rank}->finish;}
	print "All flags have been removed!<BR>\n";
}

print "</FONT>\n";
$dbh->disconnect;
#exit ; #debug
sleep 3;

print qq
|<SCRIPT LANGUAGE="JavaScript">
if (top.protWindow && !top.protWindow.closed) {top.protWindow.location.reload(true);} // in case changes in validation affect protein being displayed in protWindow
|;
if ($ITEM eq 'ANALYSIS') {
	print "top.promsFrame.location=\"$promsPath{cgi}/startValidation.cgi?ID=$itemID\";\n";
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

####>Revision history<####
# 1.2.6 Remove VAR_MOD and VMOD from script (GA 23/05/13)
# 1.2.5 Multi-databank search (PP 18/12/12)
# 1.2.4 Minor code formatting (PP 11/04/11)
# 1.2.3 Add template management (FY 04/03/11)
# 1.2.2 Minor update to filter well on the elution time (garras 06/01/11)<BR>See 2.6.5 modification in storeAnalyses.cgi