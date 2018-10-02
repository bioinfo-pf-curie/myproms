#!/usr/local/bin/perl -w

################################################################################
# filterAnalysis.cgi           2.0.1                                           #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Applies and removes filters on 1 or multiple Analyses                        #
# Filter is a list of proteins from another Analysis and/or a file             #
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
use File::Copy::Recursive qw(dircopy);
use File::Copy;


#print header;warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my $minPepScore=&promsConfig::getMinScore;

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

#####################################
####>Arguments / Processing form<####
#####################################
my @listAnalyses=param('anaList');
my $projectID=param('id_project');
my $removeFilter=param('remFilter');
my $analysisStrg=(scalar @listAnalyses > 1)? 'Analyses' : 'Analysis';

#########################
####>Processing form<####
#########################
&processForm if param('apply');

#########################
####>Removing filter<####
#########################
&removeFilter if $removeFilter;

############################################
####>Fetching selected analyses parents<####
############################################
my (@anaHierarchy,@anaValidStatus);
foreach my $anaID (@listAnalyses) {
	my @anaInfo=&promsMod::getItemInfo($dbh,'ANALYSIS',$anaID);
	my @anaParents;
	for (my $i=1;$i<=$#anaInfo;$i++) { # $i=1 (skip PROJECT)
		next if $anaInfo[$i-1]{'ITEM'} eq 'SPOT'; # skip sample if parent is SPOT
		push @anaParents,$anaInfo[$i]{'ITEM'} if $anaInfo[$i]{'ITEM'} ne 'ANALYSIS';
		push @anaParents,$anaInfo[$i]{'NAME'};
	}
	push @anaHierarchy,\@anaParents;
	push @anaValidStatus,$anaInfo[-1]{'VALID'};
}

###############################
####>Building project tree<####
###############################
my %analysisClass=(-1=>'no_scan',0=>'no_val',1=>'part_val',2=>'val');

####>Queries<####
my ($projectName)=$dbh->selectrow_array("SELECT NAME FROM PROJECT WHERE ID_PROJECT=$projectID");
my $sthExp=$dbh->prepare("SELECT ID_EXPERIMENT,NAME FROM EXPERIMENT WHERE ID_PROJECT=$projectID ORDER BY DISPLAY_POS ASC");
my $sthG2D=$dbh->prepare("SELECT ID_GEL2D,NAME FROM GEL2D WHERE ID_EXPERIMENT=? ORDER BY DISPLAY_POS ASC");
my $sthSpot=$dbh->prepare("SELECT SPOT.ID_SPOT,SPOT.NAME,ID_SAMPLE,SAMPLE.NAME FROM SPOT,SAMPLE WHERE SPOT.ID_SPOT=SAMPLE.ID_SPOT AND ID_GEL2D=? ORDER BY SPOT.NAME ASC");
my $sthSamp=$dbh->prepare("SELECT ID_SAMPLE,NAME FROM SAMPLE WHERE ID_EXPERIMENT=? AND ID_SPOT IS NULL ORDER BY DISPLAY_POS ASC");
my $sthAna=$dbh->prepare("SELECT ID_ANALYSIS,NAME,VALID_STATUS,VERIFIED_MG FROM ANALYSIS WHERE ID_SAMPLE=? ORDER BY DISPLAY_POS ASC");

####>Project<####
my %treeOptions;
foreach my $anaID (@listAnalyses) {$treeOptions{'DISCHECK'}{"analysis:$anaID"}=1;}
@{$treeOptions{'SELECTION'}}=('analysis',$listAnalyses[0]);
$treeOptions{'BUTTON'}='font-size:9px;width:55px';
%{$treeOptions{'CHECKBOX'}}=('experiment'=>1,'gel2d'=>1,'spot'=>1,'sample'=>1,'analysis'=>1); # all items are potentially checkable
my @projectTree=(0,'project',$projectID,'','',1,$projectName,''); #(depth,type,ID,labelClass,imageClass,selectable,name,popup,refChild1,...2,...3)

###>Experiments<###
$sthExp->execute;
while (my ($expID,$expName)=$sthExp->fetchrow_array) {
	my @experimentTree=(1,'experiment',$expID,'','',1,$expName,'');
	my $chkBoxExp=0;

	##>Gels<##
	$sthG2D->execute($expID);
	while (my ($gelID,$gelName)=$sthG2D->fetchrow_array) {
		my @gelTree=(2,'gel2d',$gelID,'','',1,$gelName,'');
		my $chkBoxGel=0;

		##>Spots<##
		$sthSpot->execute($gelID);
		while (my ($spotID,$spotName,$sampID,$sampName)=$sthSpot->fetchrow_array) {
			#my $sampPopup=($sampID)? "&nbsp;<IMG src=\\\'$promsPath{images}/$itemIcones{sample}\\\'/><B>$sampName</B>" : '<B>No sample</B>';
			my @spotTree=(3,'spot',$spotID,'','',1,$spotName,'');
			my $chkBoxSpot=0;

			#>Analyses<#
			if ($sampID) {
				$sthAna->execute($sampID);
				while (my ($anaID,$anaName,$validStatus,$verifMG)=$sthAna->fetchrow_array) {
					push @spotTree,[4,'analysis',$anaID,'',$analysisClass{$validStatus},1,$anaName,''];
					$chkBoxSpot++ unless $treeOptions{'DISCHECK'}{"analysis:$anaID"};
				}
			}
			push @gelTree,\@spotTree;
			if ($chkBoxSpot > 0) {$chkBoxGel++;}
			else {$treeOptions{'DISCHECK'}{"spot:$spotID"}=1;}
		}
		push @experimentTree,\@gelTree;
		if ($chkBoxGel > 0) {$chkBoxExp++;}
		else {$treeOptions{'DISCHECK'}{"gel2d:$gelID"}=1;}
	}

	##>Free Samples<##
	$sthSamp->execute($expID);
	while (my ($sampID,$sampName)=$sthSamp->fetchrow_array) {
		my @sampleTree=(2,'sample',$sampID,'','',1,$sampName,'');
		my $chkBoxSamp=0;

		#>Analyses<#
		$sthAna->execute($sampID);
		while (my ($anaID,$anaName,$validStatus,$verifMG)=$sthAna->fetchrow_array) {
			push @sampleTree,[3,'analysis',$anaID,'',$analysisClass{$validStatus},1,$anaName,''];
			$chkBoxSamp++ unless $treeOptions{'DISCHECK'}{"analysis:$anaID"};
		}
		push @experimentTree,\@sampleTree;
		if ($chkBoxSamp > 0) {$chkBoxExp++;}
		else {$treeOptions{'DISCHECK'}{"sample:$sampID"}=1;}
	}
	push @projectTree,\@experimentTree;
	$treeOptions{'DISCHECK'}{"experiment:$expID"}=1 if $chkBoxExp<=0;
}
$sthExp->finish;
$sthG2D->finish;
$sthSpot->finish;
$sthSamp->finish;
$sthAna->finish;

#######################################################
####>Fetching list of species from item's analyses<####
#######################################################
my %organismList;
#my $sthOrg=$dbh->prepare("SELECT DISTINCT(ORGANISM) FROM PROTEIN_VALIDATION WHERE SEL_STATUS>-2 AND ID_ANALYSIS=?");
#my $sthNum=$dbh->prepare("SELECT COUNT(*) FROM PROTEIN_VALIDATION WHERE SEL_STATUS>-2 AND ID_ANALYSIS=? AND ORGANISM=?"); # ignore already excluded/filtered proteins
my $sthOrgNum=$dbh->prepare("SELECT ORGANISM,COUNT(ORGANISM) FROM PROTEIN_VALIDATION WHERE SEL_STATUS>-2 AND ID_ANALYSIS=? AND IDENTIFIER NOT LIKE 'DECOY_%' GROUP BY ORGANISM"); # ignore already excluded/filtered proteins
foreach my $anaID (@listAnalyses) {
	$sthOrgNum->execute($anaID);
	while (my ($org,$numOrg)=$sthOrgNum->fetchrow_array) {
		#$sthNum->execute($anaID,$org);
		#(my $numOrg)=$sthNum->fetchrow_array;
		$organismList{$org}+=$numOrg;
	}
}
#$sthOrg->finish;
#$sthNum->finish;
$sthOrgNum->finish;

$dbh->disconnect;


##############
####>HTML<####
##############
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
print header(-'content-encoding'=>'no',-charset=>'UTF-8');
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
	TD {font-size:13px;font-weight:bold;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
|;
	#&promsMod::popupInfo();
	&promsMod::writeJsTreeFunctions(\@projectTree,\%treeOptions);
	print qq
|//--->Local functions<---//
function actionOnSelect(tabIndex) {}
function actionOnCheck(tabIndex) {}

function enableSmProt(smPbox) {
	if (smPbox.checked==true) {
		document.filterForm.smPepSc.disabled=false;
	}
	else {
		document.filterForm.smPepSc.disabled=true;
	}
}
|;

###>Species list<####
print qq
|function writeOrganisms() {
	document.filterForm.apply.disabled=true; // apply safety
	var doc=orgFrame.window.document;
	doc.open();
	doc.writeln('<HTML><HEAD>');
	doc.writeln('<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">');
	doc.writeln('<STYLE type="text/css">');
	doc.writeln('ACRONYM {cursor:help;}');
	doc.writeln('</STYLE></HEAD>');
	doc.writeln('<BODY bgcolor="white" topmargin="0" leftmargin="0" rigthmargin="0">');
	doc.writeln('<TABLE border=0 cellpadding=0 width=470>');
|;
my ($pos,$count)=(0,0);
my $disabledStrg=(scalar (keys %organismList) > 1)? '' : 'disabled'; # not selectable if only 1 organism
foreach my $org (sort{lc($a) cmp lc($b)} keys %organismList) { # sort order is critical!!!
	print "\tdoc.writeln('<TR>')\n" if $count==0;
	print "\tdoc.writeln('<TD nowrap width=235><INPUT type=checkbox value=\"$pos\" onClick=\"javascript:parent.updateOrganisms($pos,this.checked)\" $disabledStrg>";
	my $orgShort=&promsMod::resize($org,25);
	$orgShort=&promsMod::HTMLcompatible($orgShort);
	(my $orgSafe=$org)=~s/['"]/\?/g;
	if ($org eq 'unknown organism') {print "<ACRONYM title=\"Species information could not be extracted from sequence file.\">$org</ACRONYM>";}
	elsif ($orgShort =~/\.{3}/) {print "<ACRONYM title=\"$orgSafe\">$orgShort</ACRONYM>";} # => was truncated
	else {print $orgShort;}
	print " ($organismList{$org})</TD>');\n";
	$count++;
	if ($count==2) {
		print "\tdoc.writeln('</TR>')\n";
		$count=0;
	}
	$pos++;
}
print qq
|	doc.writeln('</TABLE>');
	doc.writeln('</BODY>');
	doc.writeln('</HTML>');
	doc.close();
	document.filterForm.apply.disabled=false; // release safety
}
function updateOrganisms(boxPos,checkedStatus) {
	document.filterForm.organisms[boxPos].checked=checkedStatus;
}
function checkForm(myForm) {
	//File
	if (myForm.protFile.value && (myForm.fileAction[0].checked==false && myForm.fileAction[1].checked==false)) {
		alert('Specify what should be done with proteins from file.');
		return false;
	}
	//Organism
	var orgChecked=0;
	if (myForm.organisms.length) { // if displayed => there is more than 1 organism
		for (i = 0; i < myForm.organisms.length; i++){
			if (myForm.organisms[i].checked == true) {
				orgChecked=1;
				break;
			}
		}
	}
	if (orgChecked && (myForm.orgAction[0].checked==false && myForm.orgAction[1].checked==false)) {
		alert('Specify what should be done with proteins from selected organism(s).');
		return false;
	}
	//Project item
	myForm.checkedItems.value=getCheckedItems();
	if (myForm.checkedItems.value && (myForm.anaAction[0].checked==false && myForm.anaAction[1].checked==false)) {
		alert('Specify what should be done with proteins from selected $analysisStrg.');
		return false;
	}
	if (myForm.protFile.value \|\| orgChecked \|\| myForm.checkedItems.value) { // filter(s) selected
		return true;
	}
	//Duplicate only?
	if (myForm.duplicate.checked==true) { // no filter selected
		if (confirm('No filter is selected. Duplicate $analysisStrg anyway ?')) {
			return true;
		}
		return false;
	}
	alert('You must select a valid filter.');
	return false;
}
function cancelAction() {
	//top.promsFrame.selectedAction='summary'; // set default action to 'summary'
	top.promsFrame.optionFrame.selectOption();
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" onload="writeOrganisms()">
<CENTER>
<FONT class="title">Apply Filter(s) on Selected $analysisStrg</FONT><BR>
<FORM name="filterForm" method="post" onsubmit="return checkForm(this);" enctype="multipart/form-data">
<INPUT type="hidden" name="id_project" value="$projectID" />
<INPUT type="hidden" name="checkedItems" value="" />
|;
my %itemIcones=&promsConfig::getItemIcones;
foreach my $anaID (@listAnalyses) {print "<INPUT type=\"hidden\" name=\"anaList\" value=\"$anaID\"/>\n";}
print "<BR><TABLE border=0 cellspacing=0 cellpadding=1>\n<TR bgcolor=\"$darkColor\"><TH align=left colspan=2>&nbsp;<FONT class=\"title2\">Selected Analyses</FONT>&nbsp;</TH></TR>\n";
my $bgColor=$darkColor;
my %prevItemName;
for (my $a=0;$a<=$#anaHierarchy;$a++) {
	my $refParents=$anaHierarchy[$a];
	##>Row color
	my $fatherIt=$refParents->[-3];
	if ($fatherIt && (!$prevItemName{$fatherIt} || $prevItemName{$fatherIt} ne $refParents->[-2])) { # keep color if same analysis parent item (SAMPLE or SPOT)
		$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
	}
	my $parentStrg='';
	for (my $i=0;$i<=$#{$refParents}-2;$i+=2) { # stops before ana name
		my $IT=$refParents->[$i];
		my $itName=$refParents->[$i+1];
		if ($prevItemName{$IT} && $prevItemName{$IT} eq $itName) {
			$parentStrg.="<FONT color=\"$bgColor\">$itName&nbsp;&nbsp;&gt;&nbsp;&nbsp;</FONT>";
		}
		else {
			$parentStrg.="$itName&nbsp;&nbsp;&gt;&nbsp;&nbsp;";
			$prevItemName{$refParents->[$i]}=$itName;
			for (my $j=$i+2;$j<$#{$refParents}-1;$j+=2) {$prevItemName{$refParents->[$j]}='';}
		}
	}
	##>Analysis
	my $anaCode=($anaValidStatus[$a]==-1)? 'analysis:no_scan' : ($anaValidStatus[$a]==0)? 'analysis:no_val' : 'analysis:part_val'; # cannot be 2
	print qq
|<TR valign=top bgcolor=$bgColor><TH nowrap align=left valign=middle>&nbsp;&nbsp;$parentStrg</TH>
<TH nowrap align=left valign=middle><IMG src="$promsPath{images}/$itemIcones{$anaCode}">&nbsp;$refParents->[-1]&nbsp;</TH>
</TR>
|;
}
print qq |</TABLE>
<FONT class="title3"><FONT color='#DD0000'>Warning:</FONT> Filters apply to data in <FONT color='#DD0000'>Validation</FONT> section only.
<BR>You must <FONT color='#DD0000'>Report Validation</FONT> in order to transmit changes to validated data.</FONT>
<BR>
</CENTER>\n<DIV style="display:none">
|;
foreach my $org (sort{lc($a) cmp lc($b)} keys %organismList) { # sort order is critical!!!
	$org=~s/'/\\q/g;
	$org=~s/"/\\Q/g;
	print "<INPUT type=\"checkbox\" name=\"organisms\" value=\"$org\">\n";
}
print qq
|</DIV>

<TABLE border=0 cellpadding=0 cellspacing=0>
<TR>
	<TD width=40px></TD>
	<TD nowrap valign=top><TABLE bgcolor=$darkColor border=0>
	<TR><TH nowrap><FONT style="font-size:18px;">Use proteins from items selected on the right &nbsp&nbsp&nbsp>>></FONT></TH></TR>
	<TR bgcolor=$lightColor>
		<TD><FONT style="line-height:5px;"><BR></FONT><FONT style="line-height:5px;"><BR></FONT>&nbsp&nbsp&nbsp;If partially-validated Analyses are used as filter: Only validated proteins will be used.</FONT><BR>
			<FONT style="line-height:5px;"><BR></FONT>
			&nbsp&nbsp&nbsp;If non-validated Analyses are used as filter:<BR>
			&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp<INPUT type="checkbox" name="rejProt" checked>Ignore already rejected proteins.<BR>
			&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp<INPUT type="checkbox" name="smPep" onclick="enableSmProt(this)">
			Ignore proteins matched by a single peptide with score &lt
			<INPUT type="text" name="smPepSc" value="$minPepScore" size="4" disabled/> (MS/MS Analyses only).
			<BR>
			<FONT style="line-height:5px;"><BR></FONT>
			<TABLE border=0 cellspacing=0 cellpadding=0>
			<TR>
				<TD align=right valign=top>&nbsp&nbsp&nbsp;Keep proteins</TD>
				<TD><INPUT type="radio" name="anaAction" value="diff">unique to <FONT color='#DD0000'>selected</FONT> $analysisStrg.<BR>
					<INPUT type="radio" name="anaAction" value="inter">common to filter and <FONT color='#DD0000'>selected</FONT> $analysisStrg.</TD>
				<TD align=left valign=top>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp<INPUT type="checkbox" name="groupAnaAction" value="byGroup">Extend filter to Match Group</TD>
			</TR></TABLE>
		</TD>
	</TR>

	<TR><TH nowrap><FONT style="font-size:18px;line-height:17px;">and/or<BR>Use proteins from selected species</FONT></TH></TR>
	<TR bgcolor=$lightColor>
		<TD nowrap><FONT style="line-height:5px;"><BR></FONT>&nbsp&nbsp&nbsp;Choose from List : <IFRAME name="orgFrame" src="" width=500 height=70 align=top></IFRAME><BR>
			<TABLE border=0 cellspacing=0 cellpadding=0>
			<TR>
				<TD align=right valign=top>&nbsp&nbsp&nbsp;Proteins from selected species must be</TD>
				<TD><INPUT type="radio" name="orgAction" value="keep">kept.<BR>
				<INPUT type="radio" name="orgAction" value="exclude">excluded.</TD>
			</TR>
			</TABLE>
		</TD>
	</TR>

	<TR><TH nowrap><FONT style="font-size:18px;line-height:17px;">and/or<BR>Use a list of proteins from a local file</FONT></TH></TR>
	<TR bgcolor=$lightColor>
		<TD nowrap><FONT style="line-height:5px;"><BR></FONT>&nbsp&nbsp&nbsp;Choose file : <INPUT type='file' name='protFile' size='70'/>&nbsp&nbsp&nbsp<BR>
			<TABLE border=0 cellspacing=0 cellpadding=0>
			<TR>
				<TD align=right valign=top>&nbsp&nbsp&nbsp;Keep proteins</TD>
				<TD><INPUT type="radio" name="fileAction" value="diff">unique to <FONT color='#DD0000'>selected</FONT> $analysisStrg.<BR>
					<INPUT type="radio" name="fileAction" value="inter">common to file and <FONT color='#DD0000'>selected</FONT> $analysisStrg.</TD>
				<TD align=left valign=top>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp<INPUT type="checkbox" name="groupFileAction" value="byGroup">Extend filter to Match Group</TD>
			</TR>
			</TABLE>
		</TD>
	</TR>

	<TR><TH nowrap><FONT style="font-size:18px;">Duplication and/or Filtering</FONT></TH></TR>
	<TR bgcolor=$lightColor>
		<TD nowrap><FONT style="line-height:5px;"><BR></FONT>&nbsp&nbsp&nbsp<INPUT type="checkbox" name="duplicate">Duplicate <FONT color='#DD0000'>selected</FONT> $analysisStrg before applying filter.
			<FONT style="line-height:5px;"><BR></FONT>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;(Validated proteins from partially-validated Analyses will not be duplicated)<BR><BR>
			<TABLE border=0>
			<TR>
				<TD align=right valign=top>&nbsp&nbsp&nbsp;Filter Name :</TD>
				<TD><INPUT type="text" name="filtName" value="" size="15"/><FONT style="font-size:9px;">&nbsp;(Will be added to filtered analyses name)</FONT></TD>
			</TR>
			<TR>
				<TD align=right valign=top>Comments :</TD>
				<TD><TEXTAREA name="filtComments" value="" rows="2" cols="60"></TEXTAREA><BR><FONT style="font-size:9px;">(Will be added to filtered analyses comments)</FONT></TD>
			</TR>
			</TABLE>
		</TD>
	</TR>
	<TR>
		<TH nowrap><INPUT type="submit" name="apply" value="Filter / Duplicate" />&nbsp &nbsp &nbsp<INPUT type="reset" value="  Clear  " />&nbsp &nbsp &nbsp<INPUT type="button" value="Cancel" onclick="cancelAction();"></TH>
	</TR>
	</FORM>
	</TABLE></TD>
|;

####>Project Tree<####
print qq
|	<TD valign=top><TABLE bgcolor=$darkColor border=0 width=100%><TR><TH nowrap><FONT class="title2">Select filter</FONT></TH></TR></TABLE>
|;
my $tableIndex=0;
&promsMod::printTree(\@projectTree,\%treeOptions,\$tableIndex,1); # (refArray, is last child)

print qq
|	</TD>
</TR></TABLE>
</BODY>
</HTML>
|;
exit;


#########################################################################################
################################<<<< Processing Form >>>>################################
#########################################################################################
sub processForm {

	####<Parameters>####
	my $anaAction=param('anaAction');
	my $totItemStrg=param('checkedItems');
	my $rejProt=param('rejProt');
	my $smPep=param('smPep');
	my $smPepScore=param('smPepSc');
	my @organisms=param('organisms');
	my $orgAction=param('orgAction');
	my $protFile=upload('protFile');
	my $fileAction=param('fileAction');
	my $duplicate=param('duplicate');
	my $filterName=(param('filtName'))? '-'.param('filtName') : ($duplicate)? '-flt' : '';
	my $filterComments=param('filtComments');
	my $applyFilter=($totItemStrg || @organisms || $protFile)? 1 : 0;
	my $maxRank=10; # &promsConfig::getMaxRank; # Maximum number of interpretations allowed per query
	my @percent=(10,20,30,40,50,60,70,80,90,100); # needed to monitor process progression
	my $groupFileAction=param('groupFileAction') ;  #fp
	my $groupAnaAction=param('groupAnaAction') ;

	###Building param string for history###
	my $paramStrg = '';
	$paramStrg.="anaAction:radio:$anaAction;" if $anaAction;
	$paramStrg.="groupAnaAction:check:$groupAnaAction;" if($groupAnaAction);
	$paramStrg.="rejProt:check:$rejProt;" if($rejProt);
	$paramStrg.="smPepScore:text:$smPepScore;" if($smPep);
	if($orgAction && scalar(@organisms)){
		$paramStrg.="orgAction:radio:$orgAction;orgList:list:";
		foreach my $orgName (@organisms){
			$paramStrg.="$orgName,";
		}
		$paramStrg=~ s/,\Z//;
		$paramStrg.=";";
	}
	if($protFile){
		my $protFileName = param('protFile');
		$paramStrg.="protFile:file:$protFileName;";
		$paramStrg.="fileAction:radio:$fileAction;";
	}
	$paramStrg.="groupFileAction:check:$groupFileAction;" if($groupFileAction);
	#$paramStrg.="totItemStrg:strg:$totItemStrg;" if($totItemStrg);

	##############
	####>HTML<####
	##############
	print header(-charset=>'utf-8'),"\n";
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
//BODY {font-weight:bold;}
</STYLE>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Applying Filter(s) on Selected Analyses</FONT><BR><BR>
<IMG src="$promsPath{images}/engrenage.gif">
</CENTER>
|;

	##############################
	####<Generating filter(s)>####
	##############################
	my @filteredAnalyses;
	my %keepProteins; # list of proteins to be kept
	my %removeProteins; # list of proteins to be removed

	if ($applyFilter) {
		print "<FONT class=\"title2\">Generating filter..." if $applyFilter;

		####<Filter from file>####
		if ($protFile) {
			###<Fetching protein list from file>###
			my %tempListProt;  # FP
			my $refList=($fileAction eq 'diff')? \%removeProteins : \%keepProteins;
			while(<$protFile>) {
				# $$refList{$1}=1 if /^\s*>*(\S+)/;  # FP
				if (/^\s*>*(\S+)/) {
					$$refList{$1}=1 ;
					$tempListProt{$1}=1;
				}
			}
			print '.';

			####Match group filter####  #FP
			if ($groupFileAction) {
				my $sthSelProtMatchGroup=$dbh->prepare("SELECT IDENTIFIER FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=? AND MATCH_GROUP=(SELECT MATCH_GROUP FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=? AND IDENTIFIER=?)");
				foreach my $anaID (@listAnalyses) {
					#$sthSelProtMatchGroup->execute($anaID);
					foreach my $protID (keys (%tempListProt)) {
						$sthSelProtMatchGroup->execute($anaID,$anaID,$protID);
						while (my ($protIDsup)=$sthSelProtMatchGroup->fetchrow_array) {
							$$refList{$protIDsup}=1;
						}
					}
				}
				$sthSelProtMatchGroup->finish;
			}
		}

		####<Filter from organisms>####
		if (@organisms) {
			my $sthSelOrg=$dbh->prepare("SELECT IDENTIFIER FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=? AND ORGANISM=?");
			my $refList=($orgAction eq 'exclude')? \%removeProteins : \%keepProteins;
			foreach my $anaID (@listAnalyses) {
				foreach my $org (@organisms) {
		# 			print "*$org=>";
					$org=~s/\\q/'/g; # reverse conversion
					$org=~s/\\Q/"/g;
		# 			print "$org*";
					$sthSelOrg->execute($anaID,$org);
					while(my ($identifier)=$sthSelOrg->fetchrow_array) {
						$$refList{$identifier}=1;
					}
					print '.';
				}
			}
			$sthSelOrg->finish;
		}

		####<Filter from Project items>####
		if ($totItemStrg) {

			my %sthAnaList=(
				'experiment'=>$dbh->prepare('SELECT ID_ANALYSIS,VALID_STATUS,MS_TYPE,DATA_FILE FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_EXPERIMENT=?'),
				'gel2d'=>$dbh->prepare('SELECT ID_ANALYSIS,VALID_STATUS,MS_TYPE,DATA_FILE FROM ANALYSIS,SAMPLE,SPOT WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT=SPOT.ID_SPOT AND ID_GEL2D=?'),
				'spot'=>$dbh->prepare('SELECT ID_ANALYSIS,VALID_STATUS,MS_TYPE,DATA_FILE FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_SPOT=?'),
				'sample'=>$dbh->prepare('SELECT ID_ANALYSIS,VALID_STATUS,MS_TYPE,DATA_FILE FROM ANALYSIS WHERE ID_SAMPLE=?'),
				'analysis'=>$dbh->prepare('SELECT ID_ANALYSIS,VALID_STATUS,MS_TYPE,DATA_FILE FROM ANALYSIS WHERE ID_ANALYSIS=?') # after submit only
			);

			my $sthSelVal=$dbh->prepare("SELECT IDENTIFIER FROM PROTEIN,ANALYSIS_PROTEIN WHERE ANALYSIS_PROTEIN.ID_PROTEIN=PROTEIN.ID_PROTEIN AND ID_ANALYSIS=?");
			my (%sthSelTmp,%tempListProt);
			my $selString=($rejProt)? '(SEL_STATUS=-1 OR SEL_STATUS>=1)' : 'SEL_STATUS>=-3';
			if ($smPep) { # single match exclusion
				$sthSelTmp{'MIS'}=$dbh->prepare("SELECT IDENTIFIER FROM PROTEIN_VALIDATION WHERE $selString AND (MAX_MATCH>1 OR MAX_SCORE>=$smPepScore) AND ID_ANALYSIS=?");
				$sthSelTmp{'PMF'}=$dbh->prepare("SELECT IDENTIFIER FROM PROTEIN_VALIDATION WHERE $selString AND ID_ANALYSIS=?");
			}
			else {
				$sthSelTmp{'MIS'}=$sthSelTmp{'PMF'}=$dbh->prepare("SELECT IDENTIFIER FROM PROTEIN_VALIDATION WHERE $selString AND ID_ANALYSIS=?");
			}
			my $refList=($anaAction ne 'inter')? \%removeProteins : \%keepProteins;

			$paramStrg.="anaFileList:";
			foreach my $itemInfo (split(/\+/,$totItemStrg)) { # ex: experiment:1,2,3+sample:1,3+analysis:2
				my ($chkItem,@listIDs)=split(/[:,]/,$itemInfo);
				foreach my $itID (@listIDs) {
					$sthAnaList{$chkItem}->execute($itID);
					while (my ($anaID,$validStatus,$msType,$dataFile)=$sthAnaList{$chkItem}->fetchrow_array) {
						$paramStrg.="$dataFile,";
						my $sthSelProt=($validStatus<=0)? $sthSelTmp{$msType} : $sthSelVal;
						$sthSelProt->execute($anaID);
						while (my ($identifier)=$sthSelProt->fetchrow_array) {
							$$refList{$identifier}=1;
							$tempListProt{$identifier}=1;
						}
						print '.';
					}
				}
			}
			$paramStrg=~ s/,\Z/;/;
			foreach my $sth (values %sthAnaList) {$sth->finish;}
			$sthSelVal->finish;
			$sthSelTmp{'MIS'}->finish;
			$sthSelTmp{'PMF'}->finish;

			####Match group filter####  #FP
			if ($groupAnaAction) {
				my $sthSelProtMatchGroup=$dbh->prepare("SELECT IDENTIFIER FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=? AND MATCH_GROUP=(SELECT MATCH_GROUP FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=? AND IDENTIFIER=?)");
				foreach my $anaID (@listAnalyses) {
					#$sthSelProtMatchGroup->execute($anaID);
					foreach my $protID (keys (%tempListProt)) {
						$sthSelProtMatchGroup->execute($anaID,$anaID,$protID);
						while (my ($protIDsup)=$sthSelProtMatchGroup->fetchrow_array) {
							$$refList{$protIDsup}=1;
						}
					}
				}
				$sthSelProtMatchGroup->finish;
			}
		}

		print " Done.";
		my $numFiltProt=scalar (keys %removeProteins) + scalar (keys %keepProteins);
		if ($numFiltProt) {print "&nbsp;($numFiltProt proteins in filter)";}
		else {
			$dbh->disconnect;
			print "<BR><BR>WARNING: There are no proteins in filter.<BR>Filtering was terminated. Try again with other settings.</FONT><BR>\n";
			print "</BODY>\n</HTML>\n";
			exit;
		}
		print "</FONT><BR><BR>\n";
	}
	my $keepProt=(scalar keys %keepProteins)? 1 : 0;

	my $numAnaID=scalar @listAnalyses;
	my $anaString=($numAnaID==1)? 'Analysis' : 'Analyses';

	##############################
	####<Duplicating analyses>####
	##############################
	if ($duplicate) {
		print "<FONT class=\"title2\">Duplicating $numAnaID $anaString:</FONT><BR>\n";

		####<Fetching number of ANALYSIS rows to be duplicated>####
		my ($numQueryID,$numProtID);
		my %numProteins; # needed later for % progression
		my $sthNumQ=$dbh->prepare("SELECT COUNT(*) FROM QUERY_VALIDATION WHERE ID_ANALYSIS=?");
		my $sthNumP=$dbh->prepare("SELECT COUNT(*) FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=?");

		foreach my $anaID (@listAnalyses) {
			##<Queries
			$sthNumQ->execute($anaID);
			my ($numQ)=$sthNumQ->fetchrow_array;
			$numQueryID+=$numQ;
			##<Proteins
			$sthNumP->execute($anaID);
			($numProteins{$anaID})=$sthNumP->fetchrow_array;
			$numProtID+=$numProteins{$anaID};
		}
		$sthNumQ->finish;
		$sthNumP->finish;

		####<Reserving rows>####
		###<ANALYSIS
		my ($maxAnaID)=$dbh->selectrow_array("SELECT MAX(ID_ANALYSIS) FROM ANALYSIS");
		my ($minSampID)=$dbh->selectrow_array("SELECT MIN(ID_SAMPLE) FROM SAMPLE"); # required for topAnaID row (min better than max)
		my $topAnaID=$maxAnaID+$numAnaID+1;
		$dbh->do("INSERT INTO ANALYSIS (ID_ANALYSIS,ID_SAMPLE) VALUES ($topAnaID,$minSampID)");

		###<QUERY_VALIDATION
		my ($maxQueryID)=$dbh->selectrow_array("SELECT MAX(ID_QUERY) FROM QUERY_VALIDATION");
		my $topQueryID=$maxQueryID+$numQueryID+1;
		$dbh->do("INSERT INTO QUERY_VALIDATION (ID_QUERY,ID_ANALYSIS) VALUES ($topQueryID,$topAnaID)");

		###<PROTEIN_VALIDATION
		my ($maxProtID)=$dbh->selectrow_array("SELECT MAX(ID_PROT_VALID) FROM PROTEIN_VALIDATION");
		my $topProtID=$maxProtID+$numProtID+1;
		$dbh->do("INSERT INTO PROTEIN_VALIDATION (ID_PROT_VALID,ID_ANALYSIS) VALUES ($topProtID,$topAnaID)");
#		$dbh->commit;

		####<Coping data>####
		my $sthCpAna=$dbh->prepare("SELECT NAME,VALID_STATUS,ID_SAMPLE,DES,MS_TYPE,INSTRUMENT,DATA_FILE,FILE_FORMAT,WIFF_FILE,LABELING,DECOY,FALSE_DISCOVERY,MIN_SCORE,LOWER_SCORES,MAX_RANK,TAXONOMY,COMMENTS,LAB_CODE FROM ANALYSIS WHERE ID_ANALYSIS=?");
		my $sthCpAnaDb=$dbh->prepare("SELECT ID_DATABANK,DB_RANK FROM ANALYSIS_DATABANK WHERE ID_ANALYSIS=?");
		my $sthMaxDP=$dbh->prepare("SELECT MAX(DISPLAY_POS) FROM ANALYSIS WHERE ID_SAMPLE=?");
		my $sthInAna=$dbh->prepare("INSERT INTO ANALYSIS (ID_ANALYSIS,NAME,DISPLAY_POS,VALID_STATUS,ID_SAMPLE,DES,MS_TYPE,INSTRUMENT,DATA_FILE,FILE_FORMAT,WIFF_FILE,LABELING,DECOY,FALSE_DISCOVERY,MIN_SCORE,LOWER_SCORES,MAX_RANK,TAXONOMY,COMMENTS,LAB_CODE,START_DATE,UPDATE_DATE,UPDATE_USER) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,'$userID')");
		my $sthInAnaDb=$dbh->prepare("INSERT INTO ANALYSIS_DATABANK (ID_ANALYSIS,ID_DATABANK,DB_RANK) VALUES (?,?,?)");

		my $pepString='INFO_PEP1';
		my $valueString='?';
		for (my $r=2;$r<=$maxRank;$r++) {
			$pepString.=",INFO_PEP$r";
			$valueString.=',?';
		}
		my $sthCpQuery=$dbh->prepare("SELECT ID_QUERY,QUERY_NUM,EXT_SPECTRUMID,CHARGE,VALID_STATUS,MASS_DATA,MAX_SCORE,ELUTION_TIME,$pepString FROM QUERY_VALIDATION WHERE ID_ANALYSIS=?");
		my $sthInQuery=$dbh->prepare("INSERT INTO QUERY_VALIDATION (ID_ANALYSIS,ID_QUERY,QUERY_NUM,EXT_SPECTRUMID,CHARGE,VALID_STATUS,MASS_DATA,MAX_SCORE,ELUTION_TIME,$pepString) VALUES (?,?,?,?,?,?,?,?,?,$valueString)");

		my $sthCpProt=$dbh->prepare("SELECT IDENTIFIER,DB_RANK,MW,PROT_DES,SEL_STATUS,NUM_MATCH,MAX_MATCH,SCORE,MAX_SCORE,ORGANISM,CONF_LEVEL,PROT_LENGTH,MATCH_GROUP FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=? AND IDENTIFIER NOT LIKE 'DECOY_%'");
		my $sthInProt=$dbh->prepare("INSERT INTO PROTEIN_VALIDATION (ID_ANALYSIS,ID_PROT_VALID,IDENTIFIER,DB_RANK,MW,PROT_DES,SEL_STATUS,NUM_MATCH,MAX_MATCH,SCORE,MAX_SCORE,ORGANISM,CONF_LEVEL,PROT_LENGTH,MATCH_GROUP) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");

		my $sthCpMatch=$dbh->prepare("SELECT QUERY_NUM,PEP_RANK,IDENTIFIER,MATCH_INFO,MATCH_MULTI FROM RANK_PROTEIN_MATCH WHERE ID_ANALYSIS=?");
		my $sthInMatch=$dbh->prepare("INSERT INTO RANK_PROTEIN_MATCH (ID_ANALYSIS,QUERY_NUM,PEP_RANK,IDENTIFIER,MATCH_INFO,MATCH_MULTI) VALUES (?,?,?,?,?,?)");

		my $sthCpAnaMods=$dbh->prepare("SELECT ID_MODIFICATION,MODIF_TYPE,SPECIFICITY FROM ANALYSIS_MODIFICATION WHERE ID_ANALYSIS=?");
		my $sthInAnaMods=$dbh->prepare("INSERT INTO ANALYSIS_MODIFICATION (ID_ANALYSIS,ID_MODIFICATION,MODIF_TYPE,SPECIFICITY) VALUES (?,?,?,?)");
		my $sthCpMods=$dbh->prepare("SELECT ID_MODIFICATION,PEP_RANK,POS_STRING,REF_POS_STRING FROM QUERY_MODIFICATION WHERE ID_QUERY=?");
		my $sthInMods=$dbh->prepare("INSERT INTO QUERY_MODIFICATION (ID_MODIFICATION,ID_QUERY,PEP_RANK,POS_STRING,REF_POS_STRING) VALUES (?,?,?,?,?)");

		my $sthCpHis=$dbh->prepare("SELECT VAL_TYPE,QUERY_VAL_STRG,PEP_VAL_STRG,PROT_VAL_STRG,STEP,PARAM_STRG,START_DATE FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=? AND STATUS !=-1 ORDER BY ID_VAL_HISTORY"); # -1: action was performed again on validation data after report
		my $sthInHis=$dbh->prepare("INSERT INTO VALIDATION_HISTORY (ID_ANALYSIS,VAL_TYPE,QUERY_VAL_STRG,PEP_VAL_STRG,PROT_VAL_STRG,STATUS,STEP,PARAM_STRG,START_DATE) VALUES (?,?,?,?,?,?,?,?,?)");

		#Observations & reference RTs are not copied: must be reported first!

		####<Looping through Analyses>####
		my %anaMaxDispPos;
		my $numProcAnaID=0;
		foreach my $anaID (@listAnalyses) {
			$numProcAnaID++;
			print "<FONT class=\"title3\">&nbsp;&nbsp;&nbsp;Duplicating Analysis $numProcAnaID/$numAnaID...";

			###<ANALYSIS
			my $newAnaID=++$maxAnaID;
			push @filteredAnalyses,$newAnaID;
			$sthCpAna->execute($anaID);
			my ($name,$validStatus,$sampID,@anaInfo)=$sthCpAna->fetchrow_array;
			my %dbUsed;
			$sthCpAnaDb->execute($anaID);
			while (my($dbID,$dbRank)=$sthCpAnaDb->fetchrow_array) {
				$dbUsed{$dbRank}=$dbID;
			}

			##<Updating fields
			$name.='-copy';
			#$name=&promsMod::resize($name,50);

			##<Changing displayPos
			unless ($anaMaxDispPos{$sampID}) {
				$sthMaxDP->execute($sampID);
				($anaMaxDispPos{$sampID})=$sthMaxDP->fetchrow_array;
			}

			my $date=strftime("%Y-%m-%d %H:%M:%S",localtime);
			$validStatus=0 if $validStatus==1; # valid status 1 -> 0 (partially to not)-validated
			$sthInAna->execute($newAnaID,$name,++$anaMaxDispPos{$sampID},$validStatus,$sampID,@anaInfo,$date,$date);
			foreach my $dbRank (keys %dbUsed) {
				$sthInAnaDb->execute($newAnaID,$dbUsed{$dbRank},$dbRank);
			}
			$sthCpAnaMods->execute($anaID);
			while (my @anaMods=$sthCpAnaMods->fetchrow_array) {
				$sthInAnaMods->execute($newAnaID,@anaMods);
			}
			print '.';

			###<QUERY_VALIDATION
			$sthCpQuery->execute($anaID);
			while (my ($queryID,@queryInfo)=$sthCpQuery->fetchrow_array) {
				$sthInQuery->execute($newAnaID,++$maxQueryID,@queryInfo);
				$sthCpMods->execute($queryID);
				while (my ($modID,$pepR,$posStr,$refPosStr)=$sthCpMods->fetchrow_array) {
					$sthInMods->execute($modID,$maxQueryID,$pepR,$posStr,$refPosStr);
				}
			}
			print '.';

			###<PROTEIN_VALIDATION
			$sthCpProt->execute($anaID);
			while (my @protInfo=$sthCpProt->fetchrow_array) {
				$sthInProt->execute($newAnaID,++$maxProtID,@protInfo);
			}
			print '.';

			###<RANK_PROTEIN_MATCH
			$sthCpMatch->execute($anaID);
			while (my @matchInfo=$sthCpMatch->fetchrow_array) {
				$sthInMatch->execute($newAnaID,@matchInfo);
			}

			###<VALIDATION HISTORY
			$sthCpHis->execute($anaID);
			while (my ($valType,$queryVal,$pepVal,$protVal,$step,$paramStrg,$startData)=$sthCpHis->fetchrow_array) {
				$queryVal=~s/;R_[^;]+//g; # remove reported info
				$pepVal=~s/;R_[^;]+//g;
				$protVal=~s/;R_[^;]+//g;
				$sthInHis->execute($newAnaID,$valType,$queryVal,$pepVal,$protVal,0,$step,$paramStrg,$startData); # all status set to 0 (not reported)
			}

			###<Copying files
			dircopy("$promsPath{valid}/ana_$anaID","$promsPath{valid}/ana_$newAnaID");

			###<ProteomeDiscoverer file? => Update the ana file to take into this account this new analysis
			my ($dataFile,$fileFormat)=@anaInfo[3,4];
			if ($fileFormat=~/\.PDM/) {
				$dataFile =~ s/_\d+\.pdm/\.ana/; # reconstruct the name of the MSF file and change the extension
				open (ANAFILE, ">>$promsPath{valid}/multi_ana/proj_$projectID/$dataFile");
				print ANAFILE "$newAnaID\n";
				close ANAFILE;
			}
			###<Reference retention time image(s)
			foreach my $imgFile (glob("$promsPath{valid}/ana_$newAnaID/regression_RT*.png")) {
				unlink $imgFile;
			}
			
			###>PRS files? => Must renamed them !
			if (-e "$promsPath{valid}/ana_$newAnaID/PRS_ana_$anaID.xml") {
				move("$promsPath{valid}/ana_$newAnaID/PRS_ana_$anaID.xml","$promsPath{valid}/ana_$newAnaID/PRS_ana_$newAnaID.xml");
				move("$promsPath{valid}/ana_$newAnaID/PRSparam_ana_$anaID.txt","$promsPath{valid}/ana_$newAnaID/PRSparam_ana_$newAnaID.txt");
			}
			

			print " Done.</FONT><BR><BR>\n";
		}
		$sthCpAna->finish;
		$sthCpAnaDb->finish;
		$sthMaxDP->finish;
		$sthInAna->finish;
		$sthInAnaDb->finish;
		$sthCpQuery->finish;
		$sthInQuery->finish;
		$sthCpMatch->finish;
		$sthInMatch->finish;
		$sthCpAnaMods->finish;
		$sthInAnaMods->finish;
		$sthCpMods->finish;
		$sthInMods->finish;
		$sthCpHis->finish;
		$sthInHis->finish;

		####<Removing row protection>####
		$dbh->do("DELETE FROM QUERY_VALIDATION WHERE ID_QUERY=$topQueryID");
		$dbh->do("DELETE FROM PROTEIN_VALIDATION WHERE ID_PROT_VALID=$topProtID");
		$dbh->do("DELETE FROM ANALYSIS WHERE ID_ANALYSIS=$topAnaID");
	}
	else {@filteredAnalyses=@listAnalyses;} # no duplication

#	$dbh->commit;

	#########################
	####<Applying filter>####
	#########################
	if ($applyFilter) {
		print "<FONT class=\"title2\">Applying filter to $numAnaID $anaString:</FONT><BR><BR>\n";

		####<Applying filter>####
		my $sthSelAna=$dbh->prepare("SELECT NAME,COMMENTS FROM ANALYSIS WHERE ID_ANALYSIS=?");
		my $sthUpAna=$dbh->prepare("UPDATE ANALYSIS SET NAME=?,COMMENTS=?,UPDATE_DATE=?,UPDATE_USER='$userID' WHERE ID_ANALYSIS=?");
		my $sthNumProt=$dbh->prepare("SELECT COUNT(*) FROM PROTEIN_VALIDATION WHERE SEL_STATUS>=-1 AND ID_ANALYSIS=?");
		my $sthSelProt=$dbh->prepare("SELECT ID_PROT_VALID,IDENTIFIER FROM PROTEIN_VALIDATION WHERE SEL_STATUS>=-1 AND ID_ANALYSIS=?");
		my $sthUpProt=$dbh->prepare("UPDATE PROTEIN_VALIDATION SET SEL_STATUS=-3 WHERE ID_PROT_VALID=?");

		my $numProcAnaID=0;
		foreach my $anaID (@filteredAnalyses) {
			$numProcAnaID++;
			print "<FONT class=\"title3\">&nbsp&nbsp&nbsp;Filtering Analysis $numProcAnaID/$numAnaID:</FONT> <B>0%";

			###<Updating name and comments>###
			$sthSelAna->execute($anaID);
			my ($name,$comments)=$sthSelAna->fetchrow_array;
			$name.=$filterName if $filterName;
			$name=&promsMod::resize($name,50);
			if ($filterComments) {
				$comments.="\n" if $comments;
				$comments.=$filterComments;
				$comments=&promsMod::resize($comments,500);
			}
			my $date=strftime("%Y-%m-%d %H:%M:%S", localtime);
			$sthUpAna->execute($name,$comments,$date,$anaID);

			###<Filtering proteins>###
			$sthNumProt->execute($anaID);
			my ($numProteins)=$sthNumProt->fetchrow_array;
			my @limitValue;
			foreach my $pc (@percent) {push @limitValue,int(0.5+($pc*$numProteins/100));}
			my $index=0;
			my $counter1=0;
			my $counter2=0;
			my $maxCounter2=int(0.5+($numProteins/100));

			$sthSelProt->execute($anaID);
			while (my ($protID,$identifier)=$sthSelProt->fetchrow_array) {

				##<Counter>##
				$counter1++;
				if ($counter1>=$limitValue[$index]) {
					print '.....' if $maxCounter2 < 10;
					print "$percent[$index]%"; # keeping connection alive
					$index++;
					$counter2=0;
				}
				if ($maxCounter2>=10) {
					$counter2++;
					if ($counter2==$maxCounter2) {print '.'; $counter2=0;}
				}

				##<Filter>##
				if (($keepProt && !defined($keepProteins{$identifier})) || defined($removeProteins{$identifier})) {
					$sthUpProt->execute($protID);
				}
			}
			print " Done.</B><BR><BR>\n";

			&promsMod::updateAnalysisHistory($dbh,$anaID,$paramStrg,'filter');

		}
		$sthSelProt->finish;
		$sthUpProt->finish;
	}
	$dbh->commit;

	print "<BR><FONT class=\"title2\">Filtering is complete.</FONT><BR>\n";

	####<Updating navigation frames>####
	&reloadFrames($filteredAnalyses[0]);

	exit;
}

#########################################################################################
################################<<<< Removing Filter >>>>################################
#########################################################################################
sub removeFilter {

	####<Parameters>####
	my @percent=(10,20,30,40,50,60,70,80,90,100); # needed to monitor process progression

	##############
	####>HTML<####
	##############
	print header,"\n";
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
//BODY {font-weight:bold;}
</STYLE>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Removing Filter(s) from Selected Analyses</FONT><BR><BR>
<IMG src="$promsPath{images}/engrenage.gif">
</CENTER>
<BR>
|;

	#########################
	####<Removing filter>####
	#########################

	####<DB queries>####
	my $sthMsType=$dbh->prepare("SELECT MS_TYPE FROM ANALYSIS WHERE ID_ANALYSIS=?");
	my $sthSelFP=$dbh->prepare("SELECT ID_PROT_VALID,IDENTIFIER,MAX_MATCH,MAX_SCORE FROM PROTEIN_VALIDATION WHERE SEL_STATUS=-3 AND ID_ANALYSIS=?");
	my $sthSelQR=$dbh->prepare("SELECT QUERY_NUM,PEP_RANK,MATCH_MULTI FROM RANK_PROTEIN_MATCH WHERE IDENTIFIER=? AND ID_ANALYSIS=?");
	my $sthUpFP=$dbh->prepare("UPDATE PROTEIN_VALIDATION SET NUM_MATCH=?,SEL_STATUS=?,SCORE=?,CONF_LEVEL=2 WHERE ID_PROT_VALID=?");

	my $numAnaID=scalar @listAnalyses;
	my $numProcAnaID=0;
	foreach my $anaID (@listAnalyses) {
		$numProcAnaID++;

		$sthMsType->execute($anaID);
		my ($msType)=$sthMsType->fetchrow_array;

		####<Processing filtered proteins>####
		my %allRankInfo;
		$sthSelFP->execute($anaID);
		my $refFilProt=$sthSelFP->fetchall_arrayref;

		my $numEntry=scalar @{$refFilProt};
		my @limitValue=();
		foreach my $pc (@percent) {push @limitValue,int(0.5+($pc*$numEntry/100));}
		my $index=0;
		my $counter1=0;
		my $counter2=0;
		my $maxCounter2=int(0.5+($numEntry/100));
		print "<FONT class=\"title3\">&nbsp&nbsp&nbsp;Processing Analysis $numProcAnaID/$numAnaID:</FONT><BR>\n";
		print '<B>&nbsp&nbsp&nbsp&nbsp&nbsp; 0%';
		foreach my $refData (@{$refFilProt}) {
			my ($protID,$identifier,$maxMatch,$maxScore)=@{$refData};

			###<Fetching matching peptides and their selection status>###
			my ($numMatch,$protScore,$rejectedMatch)=(0,0,0);
			my $selStatus;
			$sthSelQR->execute($identifier,$anaID);
			while (my ($queryNum,$rank,$matchFreq) = $sthSelQR->fetchrow_array) {
				my $qrKey="$queryNum:$rank";
				unless ($allRankInfo{$qrKey}) { # rankInfo not already extracted from DB
					($allRankInfo{$qrKey})=$dbh->selectrow_array("SELECT INFO_PEP$rank FROM QUERY_VALIDATION WHERE QUERY_NUM=$queryNum AND ID_ANALYSIS=$anaID");
				}
				if ($allRankInfo{$qrKey} =~ /SEL=1/) {
					$numMatch++;
					if ($msType ne 'PMF') {
						my ($pepScore)=($allRankInfo{$qrKey}=~/SC=(\d+\.*\d*)/); # Score corresponding to peptide
						$protScore+=($matchFreq*$pepScore);
					}
				}
				elsif ($allRankInfo{$qrKey} =~ /SEL=-2/) {$rejectedMatch++;}
			}

			if ($numMatch) { # at least 1 peptide selected
				$selStatus=($numMatch+$rejectedMatch==$maxMatch)? 2 : 1;
			}
			else { # no peptides selected
				$selStatus=($rejectedMatch==$maxMatch)? 0 : -1;
			}

			###<Updating table with new status>###
			$protScore=$maxScore if $msType eq 'PMF';
			$sthUpFP->execute($numMatch,$selStatus,$protScore,$protID) || die $sthUpFP->errstr;

			###<Counter
			$counter1++;
			if ($counter1>=$limitValue[$index]) {
				print '.....' if $maxCounter2<10;
				print "$percent[$index]%"; # keeping connection alive
				$index++;
				$counter2=0;
			}
			if ($maxCounter2>=10) {
				$counter2++;
				if ($counter2>=$maxCounter2) {print '.'; $counter2=0;}
			}
		}
		print " Done.</B><BR><BR>\n";

		&promsMod::updateAnalysisHistory($dbh,$anaID,'','r_filt');
	}
	$sthMsType->finish;
	$sthSelFP->finish;
	$sthSelQR->finish;
	$sthUpFP->finish;
	$dbh->commit;

	print "<BR><FONT class=\"title2\">Filter removal is complete.</FONT><BR>\n";

	####<Updating navigation frames>####
	&reloadFrames($listAnalyses[0]);

	exit;
}

#################################
####<<<< Updating Frames >>>>####
#################################
sub reloadFrames {
	my $selAnaID=$_[0];
	my @anaInfo=&promsMod::getItemInfo($dbh,'ANALYSIS',$selAnaID);

	$dbh->disconnect;

	sleep 3;

	my $anaBranchID="analysis:$selAnaID";
	my $updateURL="parent.navFrame.location=\"$promsPath{cgi}/openProject.cgi?ID=$projectID&ACT=nav";
	my ($treeFrame,$targetTree);
	if ($anaInfo[2]{'ITEM'} eq 'SAMPLE') {
		$updateURL.="&branchID=$anaBranchID\"";
		$treeFrame='navFrame';
		$targetTree='projectTreeStatus';
	}
	else {
		$updateURL.="&branchID=".lc($anaInfo[2]{ITEM}).":$anaInfo[2]{ID}&itemBranchID=$anaBranchID\"";
		$treeFrame='itemFrame';
		$targetTree='itemTreeStatus';
	}

	print qq
|<SCRIPT LANGUAGE="JavaScript">
	top.promsFrame.selectedAction='summary';
	top.$targetTree=parent.$treeFrame.addItemToTree('$anaBranchID');
	$updateURL;
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

####>Revision history<####
# 2.0.1 Rename PhosphoRS .xml and .txt files with good anaID (GA 31/01/18)
# 2.0.0 Update for duplication compatibility with features up to iRT (PP 10/08/15)
# 1.9.9 Compatible with new modifications handling (GA 18/11/13)
# 1.9.8 system command removal (PP 08/11/13)
# 1.9.7 Minor bug correction in duplicate option -> remove DECOY proteins selection while filtering analysis (GA 26/06/13)
# 1.9.6 Minor bug correction in duplicate option -> organism has to be selected on selected analysis only (GA 17/06/13)
# 1.9.5 Adds LOWER_SCORES to analysis duplication (PP 29/04/13)
# 1.9.4 GPL license (PP 11/04/13)
# 1.9.3 Compatible with new Quantification & multi-databank search (PP 04/02/13)
# 1.9.2 no longer uses DATA_FILE from QUERY_VALIDATION & PROTEIN_VALIDATION tables (PP 14/06/11)
# 1.9.1 Minor bug (FY 18/04/11)
# 1.9.0 Validation history managing (FY 29/03/11)
# 1.8.9 Warning added for protein filter (PP 25/02/11)
