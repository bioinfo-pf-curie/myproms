#!/usr/local/bin/perl -w

################################################################################
# selectOption.cgi    2.7.3													   #
# Authors: P. Poullet, G. Arras & F. Yvon (Institut Curie)                     #
# Contact: myproms@curie.fr                                                    #
# Generates list of options available to user                                  #
# Displayed in optionFrame                                                     #
# Called by openProject JavaScript tree navigation function                    #
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
use promsConfig;
use promsMod;
use strict;

#print header; warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};

###################
####>Arguments<####
###################
my $branchID=param('branchID');
my ($item,$itemID)=split(':',$branchID);
$item=uc($item);

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

####>Fetching item info<####
my @itemInfo=&promsMod::getItemInfo($dbh,$item,$itemID);
my $itemType=$itemInfo[-1]{'TYPE'};
my @childItems=&promsMod::getItemChild($item); # '' if $item is ANALYSIS

####>Fetching user information<####
my $projectID=$itemInfo[0]{'ID'};
my $experimentID=($item eq 'EXPERIMENT')? $itemID : ($itemInfo[1])? $itemInfo[1]{'ID'} : 0;
my $notEditable=($itemInfo[0]{'STATUS'} <= 0)? 0 : 1; # project status <= 0: editable
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
my $projectFullAccess=($projectAccess =~ /bioinfo|mass|manag|super/)? 1 : 0;
my $anaDeleteAccess=($projectAccess =~ /bioinfo|mass|manag/)? 1 : 0;
my $nbExpCond=0;

####>Detecting validation status<####
my ($validQuery,$filterQuery,$numUnValidDat, $numDesign);
if ($item eq 'PROJECT') {
	$validQuery="SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS,SAMPLE,EXPERIMENT WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_EXPERIMENT=EXPERIMENT.ID_EXPERIMENT AND ID_PROJECT=$itemID AND VALID_STATUS=?";
	$numUnValidDat=0;
	$numDesign=0;
}
elsif ($item eq 'EXPERIMENT') {
	$validQuery="SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_EXPERIMENT=$itemID AND VALID_STATUS=?";
	$filterQuery="SELECT COUNT(*) FROM PROTEIN_VALIDATION,ANALYSIS,SAMPLE WHERE PROTEIN_VALIDATION.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SEL_STATUS=-3 AND ID_EXPERIMENT=$itemID";
	($numDesign) = $dbh -> selectrow_array("SELECT COUNT(D.ID_DESIGN) from EXPERIMENT E INNER JOIN DESIGN D on E.ID_EXPERIMENT = D.ID_EXPERIMENT where E.ID_EXPERIMENT = $itemID");
	($numUnValidDat)=$dbh->selectrow_array("SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_EXPERIMENT=$itemID AND (VALID_STATUS=0 OR VALID_STATUS=1) AND FILE_FORMAT='MASCOT.DAT'");

}
#elsif ($item eq 'GEL2D') {
#	$validQuery="SELECT COUNT(ID_ANALYSIS) FROM SPOT,SAMPLE,ANALYSIS WHERE SPOT.ID_SPOT=SAMPLE.ID_SPOT AND SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND ID_GEL2D=$itemID AND VALID_STATUS=?";
#	$filterQuery="SELECT COUNT(*) FROM PROTEIN_VALIDATION,ANALYSIS,SAMPLE,SPOT WHERE PROTEIN_VALIDATION.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SPOT.ID_SPOT=SAMPLE.ID_SPOT AND SEL_STATUS=-3 AND ID_GEL2D=$itemID";
#	($numUnValidDat)=$dbh->selectrow_array("SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS,SAMPLE,SPOT WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SPOT.ID_SPOT=SAMPLE.ID_SPOT AND ID_GEL2D=$itemID AND (VALID_STATUS=0 OR VALID_STATUS=1) AND FILE_FORMAT='MASCOT.DAT'");
#}
#elsif ($item eq 'SPOT') {
#	$validQuery="SELECT COUNT(ID_ANALYSIS) FROM SAMPLE,ANALYSIS WHERE SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND ID_SPOT=$itemID AND VALID_STATUS=?";
#	$filterQuery="SELECT COUNT(*) FROM PROTEIN_VALIDATION,ANALYSIS,SAMPLE WHERE PROTEIN_VALIDATION.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SEL_STATUS=-3 AND ID_SPOT=$itemID";
#	($numUnValidDat)=$dbh->selectrow_array("SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_SPOT=$itemID AND (VALID_STATUS=0 OR VALID_STATUS=1) AND FILE_FORMAT='MASCOT.DAT'");
#}
elsif ($item eq 'SAMPLE') {
	$validQuery="SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS WHERE ID_SAMPLE=$itemID AND VALID_STATUS=?";
	$filterQuery="SELECT COUNT(*) FROM PROTEIN_VALIDATION,ANALYSIS WHERE PROTEIN_VALIDATION.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND SEL_STATUS=-3 AND ID_SAMPLE=$itemID";
	($numUnValidDat)=$dbh->selectrow_array("SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS WHERE ID_SAMPLE=$itemID AND (VALID_STATUS=0 OR VALID_STATUS=1) AND FILE_FORMAT='MASCOT.DAT'");
}
else { # analysis
	$validQuery="SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS WHERE ID_ANALYSIS=$itemID AND VALID_STATUS=?";
	$filterQuery="SELECT COUNT(*) FROM PROTEIN_VALIDATION WHERE SEL_STATUS=-3 AND ID_ANALYSIS=$itemID";
	$numUnValidDat=0;
# 	($numUnValidDat)=$dbh->selectrow_array("SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS WHERE ID_ANALYSIS=$itemID AND VALID_STATUS=0 AND FILE_FORMAT='MASCOT.DAT'"); # only non-validated!
}
my $sthValid=$dbh->prepare($validQuery);
##>Unvalidated + not-imported Analyses
$sthValid->execute(-1);
my ($numNotImport)=$sthValid->fetchrow_array;
##>Unvalidated Analyses
$sthValid->execute(0);
my ($numUnValid)=$sthValid->fetchrow_array;
$numUnValid+=$numNotImport;
##>Partially validated Analyses
$sthValid->execute(1);
my ($numPartValid)=$sthValid->fetchrow_array;
##>Validated Analyses
$sthValid->execute(2);
my ($numValid)=$sthValid->fetchrow_array;
$sthValid->finish;
#my $numAnalyses=$numUnValid+$numPartValid+$numValid;
$numUnValid+=$numPartValid;
$numValid+=$numPartValid;
$numValid=2 if ($item eq 'ANALYSIS' && $numValid && !$numPartValid);

##>Filtered Analyses
my $filter;
if ($numUnValid && $item ne 'PROJECT') {
	($filter)=$dbh->selectrow_array($filterQuery);
}

####>Checking if item has children<####
my $hasChildren=0;
my $hasAnaQuantifs=0;
my $okReport=1;
if ($item eq 'ANALYSIS') {
	#>Quantification
	if ($projectFullAccess && $numNotImport==0 && $numValid <= 2) { # validStatus==0, 1 or 2
		my ($userQuantif)=$dbh->selectrow_array("SELECT 1 FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q,QUANTIFICATION_METHOD QM WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND QM.CODE !='SILAC' AND QM.CODE !='ITRAQ' AND QUANTIF_ANNOT NOT LIKE '%::SOFTWARE=PD%' AND AQ.ID_ANALYSIS=$itemID LIMIT 1");
		if ($userQuantif) {
			$hasChildren++;
			#$okReport=0;
		}
	}
	$okReport=0 if ($numValid==2 || $hasChildren);
	($hasAnaQuantifs)=($hasChildren)? 1 : $dbh->selectrow_array("SELECT 1 FROM ANA_QUANTIFICATION A,QUANTIFICATION Q WHERE ID_ANALYSIS=$itemID AND A.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_DESIGN IS NULL"); # internal quantifs only

	if ($anaDeleteAccess && !$hasChildren) { # Scanning ANALYSIS n/n tables
		foreach my $table ('GOANA_ANALYSIS','PATHWAYANA_ANALYSIS') { # ok to delete if ANA_COMPARISON, OBSERVATION
			($hasChildren)=$dbh->selectrow_array("SELECT 1 FROM $table WHERE ID_ANALYSIS=$itemID LIMIT 1");
			last if $hasChildren;
		}
	}
}
else {
	if ($item eq 'SPOT') { # spot virtual child is analysis (not sample)!
		($hasChildren)=$dbh->selectrow_array("SELECT 1 FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_SPOT=$itemID LIMIT 1");
	}
	else {
		foreach my $childItem (@childItems) {
			($hasChildren)=$dbh->selectrow_array("SELECT 1 FROM $childItem WHERE ID_$item=$itemID LIMIT 1");
			last if $hasChildren;
		}
	}
}

####>Checking if project contains validated proteins<####
my ($validProt)=$dbh->selectrow_array("SELECT 1 FROM PROTEIN WHERE ID_PROJECT=$projectID LIMIT 1"); # 1 is enough

####>Checking if experiment contains protein quanifications<####
my $hasProtQuantifs=0;
if ($experimentID && $validProt) {
	($hasProtQuantifs)=$dbh->selectrow_array("SELECT 1 FROM SAMPLE S,ANALYSIS A,ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND S.ID_EXPERIMENT=$experimentID AND FOCUS='protein' LIMIT 1");
}

####>Checking if project contains classifications<####
#my $classif;
#if ($projectAccess ne 'user') {
#	($classif)=$dbh->selectrow_array("SELECT COUNT(*) FROM CLASSIFICATION WHERE ID_PROJECT=$projectID");
#}

$dbh->disconnect;

####>Checking if on-going clustering<#### OBSOLETE
#my $clusteringStatus=0;
#if (-e "$promsPath{tmp}/clustering/project_$projectID/item.fsa") {
#	if (-e "$promsPath{tmp}/clustering/project_$projectID/blast_end.txt") {$clusteringStatus=2;} # blast is finished
#	else {$clusteringStatus=1;} # blast is on-going
#}

####>Checking if files in private batch directory<####
#my $batchString='';
#my $batchFilesDir="$promsPath{tmp}/batch/$userID";
#if (-e $batchFilesDir && ($projectAccess eq 'bioinfo' || $projectAccess eq 'mass') && ($item eq 'EXPERIMENT' || $item eq 'SAMPLE')) {
#	opendir (DIR, $batchFilesDir);
#	while (defined (my $file = readdir (DIR)))	{
#		next if -d "$batchFilesDir/$file"; # directory
#		$batchString='*';
#		last;
#	}
#	close DIR;
#}

#######################
####>Starting HTML<####
#######################
my ($popBgColor,$popBdColor)=&promsConfig::getPopupColors;
my $navSubClass=($itemInfo[-1]->{ITEM} eq 'SPOT' || ($itemInfo[-3] && $itemInfo[-3]->{ITEM} eq 'SPOT'))? 'navBorder2' : 'navBorder1';
print header(-'content-encoding'=>'no',-charset=>'UTF-8');
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Option Menu</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
function selectOption(selectedButton) {
	if (!selectedButton) { // 1st time this page is loaded
		selectedButton=document.getElementById(top.promsFrame.selectedAction); // default value is stored in promsFrame;
		if (!selectedButton \|\| selectedButton.id == 'monitor') { // button is not defined for current project item => switch to 'summary'
			selectedButton=document.getElementById('summary');
		}
	}
	if (selectedButton.id=='delete') {
		var confString;
		if ('$item'=='ANALYSIS' && $numValid) {
			confString='WARNING: This Analysis has already been ';
			if ($numValid==1) {confString+='partially ';}
			confString+='validated.\\nDeletion will remove all associated proteins.\\nProceed ?';
		}
		else {confString='Delete $itemType ?';}
		if (confirm(confString)) {
			top.promsFrame.resultFrame.location="$promsPath{cgi}/deleteProjectItem.cgi?ITEM=$item&ID=$itemID&ACT=delete&PROJECT_ID=$projectID";
		}
		return;
	}
	if (currentButton) {
		//currentButton.style.color='#000000';
		currentButton.className='';
	}
	currentButton=selectedButton;
	//currentButton.style.color='#DD0000';
	currentButton.className='selectedButton';
	var action=currentButton.id;

	//***Update promsFrame variable
	if (action != 'access' && action != 'remFilter') {
		top.promsFrame.selectedAction=action;
	}

	//***List of possible actions
	if (action=='summary' \|\| action=='edit') { //\|\| action=='launch'
		top.promsFrame.resultFrame.location="$promsPath{cgi}/editProjectItem.cgi?ACT="+action+"&ITEM=$item&ID=$itemID&PROJECT_ID=$projectID";
	}
	else if (action=='add') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/editProjectItem.cgi?ACT=add&ITEM=$childItems[0]&ID=$itemID&PARENT=$item&PROJECT_ID=$projectID";
	}
	else if (action=='addDesign') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/manageDesignCondition.cgi?ACT=add&ITEM=DESIGN&ID=$itemID&PARENT=$item&PROJECT_ID=$projectID";
	}
	else if (action=='addGel2D') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/editProjectItem.cgi?ACT=add&ITEM=GEL2D&ID=$itemID&PARENT=$item";
	}
	else if (action=='show2dGel') {
		var gelURL="$promsPath{cgi}/view2dGel.cgi?id_project=$projectID&id_gel=$itemID&ACT=main";
		if (!top.gel2dWindow \|\| top.gel2dWindow.closed) {
			top.gel2dWindow=window.open(gelURL,'2D_Gel','width=1050,height=500');
		}
		else if (top.gel2dWindow.id_gel != $itemID) {
			top.gel2dWindow.location=gelURL;
		}
		top.gel2dWindow.focus();
		selectOption(document.getElementById('summary')); // switch main window to summary
	}
	else if (action=='processAnalyses') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/processAnalyses.cgi?ID=$branchID";
	}
	else if (action=='reportValid') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/selectAnalyses.cgi?ID=$branchID&callType=report";
	}
	else if (action=='proteins') {
		if (top.promsFrame.selectedView.match(/ana:/) \|\| top.promsFrame.selectedView=='relDelta') top.promsFrame.selectedView='peptide';
		if (top.expandMode==null) {
			if (confirm('Show Match Groups?\\n(Match Group display is not compatible with all peptide display options)')) {top.expandMode=1;}
			else {top.expandMode=0;}
		}
		top.promsFrame.resultFrame.location="$promsPath{cgi}/listProteins.cgi?TYPE=$item&ID=$itemID&listMode="+top.promsFrame.selectedMode+"&listFilter="+top.promsFrame.selectedFilter+"&view="+top.promsFrame.selectedView+"&expMode="+top.expandMode+"&pepType="+top.promsFrame.selectedPepType+"&unClass="+top.showUnClass;
	}
	else if (action=='export') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/exportProteinList.cgi?projectID=$projectID&item=$item&itemID=$itemID&listMode="+top.promsFrame.selectedMode;
	}
	else if (action=='classification') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/listClassifications.cgi?id_project=$itemID";
	}
	else if (action=='access') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/projectAccess.cgi?ID=$itemID";
	}
	else if (action=='search') {
		if (top.promsFrame.selectedSearch=='sequence') {
			top.promsFrame.resultFrame.location="$promsPath{cgi}/searchSqProtein.cgi?ITEM=$item&ID=$itemID&id_project=$projectID";
		}
		else { // keyword
			top.promsFrame.resultFrame.location="$promsPath{cgi}/searchKwProtein.cgi?ITEM=$item&ID=$itemID&id_project=$projectID";
		}
	}
	else if (action=='compare') {
		if (top.promsFrame.selectedComparison=='items') {
			if (top.promsFrame.selectedView.match(/ana:/) \|\| top.promsFrame.selectedView=='relDelta') top.promsFrame.selectedView='peptide';
			top.promsFrame.resultFrame.location="$promsPath{cgi}/compareItems.cgi?FRAMES=1&id_project=$projectID&refITEM1=$item:$itemID&view="+top.promsFrame.selectedView;
		}
		else {
			top.promsFrame.resultFrame.location="$promsPath{cgi}/compareAnalyses.cgi?ACT=frames&id_project=$projectID&parentItem=$item:$itemID&sort="+top.promsFrame.selectedView;
		}
	}
	else if (action=='appFilter') { // only if item=ANALYSIS
		//top.promsFrame.resultFrame.location="$promsPath{cgi}/filterAnalysis.cgi?ITEM=$item&ID=$itemID&id_project=$projectID";
		top.promsFrame.resultFrame.location="$promsPath{cgi}/filterAnalysis.cgi?anaList=$itemID&id_project=$projectID";
	}
	else if (action=='remFilter') { // only if item=ANALYSIS
		if (confirm('Remove filter from all unvalidated Analyses found in selected item?')) {
			//top.promsFrame.resultFrame.location="$promsPath{cgi}/filterAnalysis.cgi?ITEM=$item&ID=$itemID&id_project=$projectID&remFilter=1";
			top.promsFrame.resultFrame.location="$promsPath{cgi}/filterAnalysis.cgi?anaList=$itemID&id_project=$projectID&remFilter=1";
		} else {autoSelectButton('summary');}
	}
	/*else if (action=='goAnalysis') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/startGOAnalysis.cgi?id_exp=$itemID";
	}
	else if (action=='goQuantiAnalysis') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/selectAnalyses.cgi?ID=$itemType:$itemID&id_project=$projectID&callType=goQuantiAnalysis";
	}*/
	else if (action=='anaQuantifs') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/showProtQuantification.cgi?id_ana=$itemID&CALL=ana";
	}
	else if (action=='compareQuantifs') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/compareQuantifications.cgi?&id_project=$projectID&parentItem=$item:$itemID";
	}
	else if (action=='exportQuantifs') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/startExploratoryAnalysis.cgi?ID=$experimentID&ACT=export"; // also used for export since v2.0.0
	}
	else if (action=='goVisualisation') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/listGO.cgi?ITEM=$item&ID=$itemID";
	}
	else if (action == 'startExploratoryAnalyses') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/startExploratoryAnalysis.cgi?ID=$experimentID";
	}
	else if (action == 'monitor'){
		var monitorWindow=window.open("$promsPath{cgi}/monitorJobs.cgi?filterDateNumber=1&filterDateType=DAY&filterStatus=Queued&filterStatus=Running",'monitorJobsWindow','width=1000,height=500,scrollbars=yes,resizable=yes');
	}
	//else if (action == 'mergeUniprotIDs') {
	//	top.promsFrame.resultFrame.location="$promsPath{cgi}/mergeUniprotIDs.cgi?ID=$itemID";
	//}
	else {alert('Unrecognized option');}
|;
#	else if (action=='seeProteins') {
#		top.promsFrame.resultFrame.location="$promsPath{cgi}/$editScript?ACT="+action+"&ITEM=$item&ID=$itemID&PROJECT_ID=$projectID";
#	}
#	else if (action=='combAna') {
#		top.promsFrame.resultFrame.location="./combineAnalyses.cgi?projID=$projectID&ITEM=$item&ID=$itemID";
#	}
#	else if (action=='autoSelect') {
#		top.promsFrame.resultFrame.location="./autoSelect.cgi?id_project=$projectID&ITEM=$item&ID=$itemID&ACT=select";
#	}
#	else if (action=='clearSelect') {
#		if (confirm('Clear peptide selection in all Analyses found in selected item?')) {
#			top.promsFrame.resultFrame.location="./autoSelect.cgi?id_project=$projectID&ITEM=$item&ID=$itemID&ACT=clear";
#		} else {autoSelectButton('summary');}
#	}
#//	else if (action=='deleteBatch') {
#//		top.promsFrame.resultFrame.location="./selectAnalyses.cgi?ID=$branchID&callType=delete";
#//	}
#	else if (action=='batchAnalyses') {
#		top.promsFrame.resultFrame.location="./importBatchAnalyses.cgi?ID=$branchID&action=start&numAna=$numAnalyses";
#	}
print qq
|}
function autoSelectButton(buttonId) { // an action was cancelled or selected from resultFrame
	currentButton.style.color='#000000';
	currentButton=document.getElementById(buttonId);
	currentButton.style.color='#DD0000';
	top.promsFrame.selectedAction=buttonId;
}
function startValidation() {
	if (top.gel2dWindow && !top.gel2dWindow.closed) {top.gel2dWindow.close();}
	top.promsFrame.location="$promsPath{cgi}/startValidation.cgi?ID=$itemID&alertTax=1";
}
var currentButton; // currently selected button
//var clusStatus=\$clusteringStatus; OBSOLETE
</SCRIPT>
</HEAD>
<BODY style="margin:0px;" onload="selectOption()">
<DIV class="navigation $navSubClass">
|;

#################################
####>Printing Item hierarchy<####
#################################
my @titleString;
my $titleChar='';
foreach my $i (0..$#itemInfo) {
	if ($i==$#itemInfo) {
		push @titleString,"<FONT color=#DD0000>$itemInfo[$i]{NAME}</FONT>\n";
		$titleChar.=$itemInfo[$i]{NAME};
	}
	else {
		push @titleString,$itemInfo[$i]{NAME};
		$titleChar.="$itemInfo[$i]{NAME} > ";
	}
}
print '<FONT class="title2">';
if (length($titleChar)<=85) {print join (' > ',@titleString);}
else {
	my @skipItems=($#titleString==3)? ('...','...') : ('...'); # item must be analysis or sample
	print join (' > ',$titleString[0],@skipItems,$titleString[-1]);
}
print "</FONT><BR>\n";

#########################
####>List of options<####
#########################
print "<TABLE cellpadding=0 cellspacing=1 border=0>\n";
print "<TR align=center valign=top>\n";


if ($notEditable || $projectAccess eq 'guest') {
	print "<TD nowrap>",&displayItemButton,"</TD>\n";
	#print "<TD nowrap>",&show2dGelButton,"</TD>\n" if $item eq 'GEL2D';
	print "</TD>\n";
	if ($validProt) {
		print "<TD nowrap>\n",&listProtButton,"<BR>\n",&exportListButton,"</TD>\n";
		print "<TD nowrap>",&compareItemsButton;
		if ($hasProtQuantifs) {
			print "<BR>\n",&compQuantifsButton,"</TD>\n";
			print "<TD nowrap>",&exportQuantifsButton,"</TD>\n" if ($item ne 'ANALYSIS' || !$hasAnaQuantifs);
		}
		else {print "</TD>\n";}
	}
	#print "<TD nowrap>",&anaQuantifsButton,"</TD>\n" if ($item eq 'ANALYSIS' && $hasAnaQuantifs);
	if ($item eq 'ANALYSIS' && $hasAnaQuantifs) {
		print "<TD nowrap>",&anaQuantifsButton;
		print "<BR>",&exportQuantifsButton if $hasProtQuantifs;
		print "</TD>\n";
	}
	print "<TD nowrap>\n",&goVisualisation,"</TD>\n" if $validProt;
	print "<TD nowrap>",&searchButton,"</TD>\n" if $validProt;
	print "<TD nowrap>",&projectAccessButton,"</TD>\n" if ($item eq 'PROJECT' && $projectAccess ne 'guest');
}
####>Bioinformatician, Massist, Manager or Super<####
elsif ($projectFullAccess) {
	print "<TD nowrap>",&displayItemButton,"<BR>\n",&editItemButton,"</TD>\n";
	print "<TD nowrap>";
	if ($item eq 'EXPERIMENT') {
		print &delItemButton,"</TD>\n<TD nowrap>" unless $hasChildren;
		print &addChildButton('SAMPLE',1),"</TD>\n"; #,"<BR>\n",&addChildButton('GEL2D',0)
		print "<TD nowrap>",&addChildButton('DESIGN',0),"</TD>\n"; # Mulitple samples can be added
	}
	else {
		if ($item eq 'PROJECT') {print &addChildButton('EXPERIMENT',1),"<BR>\n";} # Multiple experiments can be added
		#elsif ($item eq 'GEL2D') {print &show2dGelButton,"<BR>\n";}
		if ($item eq 'ANALYSIS') {
			print &validButton,"<BR>\n" if $numUnValid;
			print &delItemButton if ($anaDeleteAccess && !$hasChildren);
		}
		else {print &delItemButton unless $hasChildren;}
		#print &combAnaButton if ($numUnValidDat>=2 && ($item eq 'EXPERIMENT' || $item eq 'SAMPLE'));
		print "</TD>\n";
	}
	#if ($numUnValid && ($item eq 'EXPERIMENT' || $item eq 'SAMPLE')) {
	#	print "<TD nowrap>",&autoSelectButton,"<BR>\n",&clearSelectButton,"</TD>\n";
	#}
	if ($item ne 'PROJECT' && $item ne 'ANALYSIS') {print "<TD nowrap>",&processAnalysesButton,"</TD>\n";} # && $item ne 'QUANTIFICATION' && $item ne 'DESIGN'
	elsif ($item eq 'ANALYSIS' && $numUnValid) {
		print "<TD nowrap>";
		print &appFilterButton,"<BR>\n";
		print &remFilterButton if $filter;
		print "</TD>\n";
		print "<TD nowrap>",&reportValidationButton,"</TD>\n" if (!$numNotImport && $okReport); #unless $numNotImport;
	}
	#print "<TD nowrap>",&searchButton,"<BR>\n",&compareItemsButton,"</TD>\n";
	if ($validProt) {
		print "<TD nowrap>\n",&listProtButton,"<BR>\n",&exportListButton,"</TD>\n";
		print "<TD nowrap>",&compareItemsButton;
		if ($hasProtQuantifs) {
			print "<BR>\n",&compQuantifsButton,"</TD>\n";
			print "<TD nowrap>",&exportQuantifsButton,"</TD>\n" if ($item ne 'ANALYSIS' || !$hasAnaQuantifs);
		}
		else {print "</TD>\n";}
	}
	#print "<TD nowrap>",&anaQuantifsButton,"</TD>\n" if ($item eq 'ANALYSIS' && $hasAnaQuantifs);
	if ($item eq 'ANALYSIS' && $hasAnaQuantifs) {
		print "<TD nowrap>",&anaQuantifsButton;
		print "<BR>",&exportQuantifsButton if $hasProtQuantifs;
		print "</TD>\n";
	}
	print "<TD nowrap>\n",&goVisualisation if $validProt;
	print "<BR>",&startExploratoryAnalyses if ($item eq 'EXPERIMENT' && $hasProtQuantifs); # $validProt && $numDesign
	print "</TD>";
	#print "<TD nowrap>";
	#print &classificationButton,"<BR>\n" if $item eq 'PROJECT';
	#print &goAnalysisButton,"<BR>\n",&goQuantiAnalysisButton,"</TD>\n<TD nowrap>" if ($item eq 'EXPERIMENT' && $validProt);
	#print &searchButton,"</TD>\n";
	my $optionStrg='';
	$optionStrg.=&classificationButton."<BR>\n" if $item eq 'PROJECT';
	#$optionStrg.=&goAnalysisButton."<BR>\n".&goQuantiAnalysisButton."</TD>\n<TD nowrap>" if ($item eq 'EXPERIMENT' && $validProt);
	$optionStrg.=&searchButton if $validProt;
	print "<TD nowrap>$optionStrg</TD>\n" if $optionStrg;
	print "<TD nowrap>",&projectAccessButton,"</TD>\n" if $item eq 'PROJECT'; # if $projectAccess !~ /user/;	# admin, mass & bioinfo only
	print "<TD nowrap>",&monitorButton,"</TD>\n";
	#print "<TD nowrap>",&projectAccessButton,"<BR>",&mergeUniprotIDs,"</TD>\n" if $item eq 'PROJECT'; # if $projectAccess !~ /user/;	# admin, mass & bioinfo only
}

####>(Power) User or Administrator<####
else {
	print "<TD nowrap>",&displayItemButton,"<BR>\n";
	print &editItemButton unless $item eq 'ANALYSIS';
	print "</TD>\n";
	print "<TD nowrap>";
	if ($item eq 'EXPERIMENT') {
		print "<TD nowrap>",&delItemButton,"</TD>\n" unless $hasChildren; #,&addChildButton('GEL2D',0),"<BR>\n";
		print "<TD nowrap>",&addChildButton('SAMPLE',1),"<BR>\n",&addChildButton('DESIGN',0),"</TD>\n"; # Mulitple samples can be added
	}
	else {
		if ($item eq 'PROJECT') {print &addChildButton('EXPERIMENT',1),"<BR>\n";} # Mulitple experiments can be added
		#elsif ($item eq 'GEL2D') {print &show2dGelButton,"<BR>\n";}
		elsif ($item eq 'SAMPLE') {print &addChildButton('ANALYSIS',0),"<BR>\n";} # only 1 analysis can be added at once    $item eq 'SPOT' || 
		elsif ($item eq 'ANALYSIS' && $numUnValid && $projectAccess =~ /power/) {print &validButton,"<BR>\n";}
		print &delItemButton unless $hasChildren;
		#print &combAnaButton if ($numUnValidDat>=2 && ($item eq 'EXPERIMENT' || $item eq 'SAMPLE'));
		print "</TD>\n";
	}
	#print &combAnaButton if ($numUnValidDat>=2 && ($item eq 'EXPERIMENT' || $item eq 'SAMPLE') && $projectAccess =~ /super/);
	#if ($numUnValid && ($projectAccess =~ /power/ || $projectAccess =~ /super/)) {
	if ($item ne 'PROJECT' && $item ne 'ANALYSIS') {
		#print "<TD nowrap>",&autoSelectButton,"<BR>\n",&clearSelectButton,"</TD>\n";
		print "<TD nowrap>",&processAnalysesButton,"</TD>\n" if $numUnValid;
	}
	if ($item eq 'ANALYSIS' && $numUnValid) {
		print "<TD nowrap>";
		print &appFilterButton,"<BR>\n";
		print &remFilterButton if $filter;
		print "</TD>\n";
	}
	if ($validProt) {
		print "<TD nowrap>\n",&listProtButton,"<BR>\n",&exportListButton,"</TD>\n";
		print "<TD nowrap>",&compareItemsButton;
		if ($hasProtQuantifs) {
			print "<BR>\n",&compQuantifsButton,"</TD>\n";
			print "<TD nowrap>",&exportQuantifsButton,"</TD>\n" if ($item ne 'ANALYSIS' || !$hasAnaQuantifs);
		}
		else {print "</TD>\n";}
	}
	#print "<TD nowrap>",&anaQuantifsButton,"</TD>\n" if ($item eq 'ANALYSIS' && $hasAnaQuantifs);
	if ($item eq 'ANALYSIS' && $hasAnaQuantifs) {
		print "<TD nowrap>",&anaQuantifsButton;
		print "<BR>",&exportQuantifsButton if $hasProtQuantifs;
		print "</TD>\n";
	}
	print "<TD nowrap>\n",&goVisualisation if $validProt;
	print "<BR>",&startExploratoryAnalyses if ($item eq 'EXPERIMENT' && $hasProtQuantifs); # $validProt && $numDesign
	print "</TD>" if $validProt;
	#print "<TD nowrap>";
	#print &classificationButton,"<BR>\n" if $item eq 'PROJECT';
	#print &goAnalysisButton,"<BR>\n",&goQuantiAnalysisButton,"</TD>\n<TD nowrap>" if ($item eq 'EXPERIMENT' && $validProt);
	#print &searchButton,"</TD>\n";
	my $optionStrg='';
	$optionStrg.=&classificationButton."<BR>\n" if $item eq 'PROJECT';
	#$optionStrg.=&goAnalysisButton."<BR>\n".&goQuantiAnalysisButton."</TD>\n<TD nowrap>" if ($item eq 'EXPERIMENT' && $validProt);
	$optionStrg.=&searchButton if $validProt;
	print "<TD nowrap>$optionStrg</TD>\n" if $optionStrg;
	print "<TD nowrap>",&projectAccessButton,"</TD>\n" if $item eq 'PROJECT';
	#print "<TD nowrap>",&projectAccessButton,"<BR>",&mergeUniprotIDs,"</TD>\n" if $item eq 'PROJECT';
}

print qq
|</TR></TABLE>
</DIV>
</BODY>
</HTML>
|;

##############################<<< Buttons Subroutines >>>############################

####<Display Item>####
sub displayItemButton {
	return "<INPUT type=\"button\" id=\"summary\" style=\"width:80px\" value=\"Summary\" onclick=\"selectOption(this)\">";
}

####<Edit Item>####
sub editItemButton {
	my ($disabled)=@_;
	$disabled='' unless $disabled;
	return "<INPUT type=\"button\" id=\"edit\" style=\"width:80px\" value=\"Edit\" onclick=\"selectOption(this)\" $disabled>";
}
####<Add Child Item>####
sub addChildButton {
	#my ($addedChild,$buttonID)=($_[0])? ($_[0],'addGel2D') : ($childItems[0],'add');
	my ($addedChild,$multi)=@_;
	my $childType=&promsMod::getItemType($addedChild);
	$childType.='(s)' if $multi;
	my $buttonID=($addedChild eq 'GEL2D')? 'addGel2D' : ($addedChild eq 'DESIGN')? 'addDesign' : 'add';
	return "<INPUT type=\"button\" id=\"$buttonID\" style=\"width:140px\" value=\"Add $childType\" onclick=\"selectOption(this)\">";
}
####<Add 2D gel>####
sub show2dGelButton {
	return "<INPUT type=\"button\" id=\"show2dGel\" style=\"width:140px;font-weight:bold\" value=\"Display 2D Gel\" onclick=\"selectOption(this)\">";
}
####<Delete Item>####
sub delItemButton {
	return "<INPUT type=\"button\" id=\"delete\" style=\"width:140px\" value=\"Delete\" onclick=\"selectOption(this)\">";
}
####<Start Validation>####
sub validButton {
	return "<INPUT type=\"button\" style=\"width:140px;font-weight:bold\" value=\"Validation\" onclick=\"startValidation()\">";
}
####<Clear Peptide Selection>####
sub clearSelectButton {
	return "<INPUT type=\"button\" id=\"clearSelect\" style=\"width:130px\" value=\"Clear Selection\" onclick=\"selectOption(this)\">";
}
####<Apply Filter>####
sub appFilterButton {
	return "<INPUT type=\"button\" id=\"appFilter\" style=\"width:130px\" value=\"Filter / Duplicate\" onclick=\"selectOption(this)\">";
}
####<Remove Filter>####
sub remFilterButton {
	return "<INPUT type=\"button\" id=\"remFilter\" style=\"width:130px\" value=\"Remove Filter\" onclick=\"selectOption(this)\">";
}
####<Display Protein List>####
sub listProtButton {
	return "<INPUT type=\"button\" id=\"proteins\" style=\"width:120px\" value=\"List Proteins\" onclick=\"selectOption(this)\">";
}
####<Export Protein List>####
sub exportListButton {
	return "<INPUT type=\"button\" id=\"export\" style=\"width:120px\" value=\"Export Proteins\" onclick=\"selectOption(this)\">";
}
####<Manage Classification>####
sub classificationButton {
	return "<INPUT type=\"button\" id=\"classification\" style=\"width:160px\" value=\"Manage Custom Lists\" onclick=\"selectOption(this)\">";
}
####<Project Access>####
sub projectAccessButton {
	return "<INPUT type=\"button\" id=\"access\" style=\"width:160px\" value=\"Project Accessibility\" onclick=\"selectOption(this)\">";
}
####<Search for a Protein>####
sub searchButton {
	return "<INPUT type=\"button\" id=\"search\" style=\"width:160px\" value=\"Search for Proteins\" onclick=\"selectOption(this)\">";
}
####<Compare Items>####
sub compareItemsButton {
	return "<INPUT type=\"button\" id=\"compare\" style=\"width:175px\" value=\"Compare Project Items\" onclick=\"selectOption(this)\">";
}
####<Send2Biologists>####
sub reportValidationButton {
	return "<INPUT type=\"button\" id=\"reportValid\" style=\"width:140px\" value=\"Report Validation\" onclick=\"selectOption(this)\">";
}
####<Process Analyses>####
sub processAnalysesButton {
	return "<INPUT type=\"button\" id=\"processAnalyses\" style=\"width:140px;font-weight:bold\" value=\"Process Analyses\" onclick=\"selectOption(this)\">";
}
####<GO analysis>####
#sub goAnalysisButton {
#	return "<INPUT type=\"button\" id=\"goAnalysis\" style=\"width:160px\" value=\"Start GO Analysis\" onclick=\"selectOption(this)\">";
#}
####<GO analysis>####
#sub goQuantiAnalysisButton {
#	return "<INPUT type=\"button\" id=\"goQuantiAnalysis\" style=\"width:160px\" value=\"Start Q. GO Analysis\" onclick=\"selectOption(this)\">";
#}
####<Ana-level quantification>####
sub anaQuantifsButton {
	return "<INPUT type=\"button\" id=\"anaQuantifs\" style=\"width:175px\" value=\"Internal Quantifications\" onclick=\"selectOption(this)\">";
}
####<Jobs monitoring>####
sub monitorButton {
	return "<INPUT type=\"button\" id=\"monitor\" style=\"width:125px\" value=\"Monitor Jobs\" onclick=\"selectOption(this)\">";
}
####<Compare Quantifications>####
sub compQuantifsButton {
	return "<INPUT type=\"button\" id=\"compareQuantifs\" style=\"width:175px\" value=\"Compare Quantifications\" onclick=\"selectOption(this)\">";
}
####<Export Quantifications>####
sub exportQuantifsButton {
	return "<INPUT type=\"button\" id=\"exportQuantifs\" style=\"width:175px\" value=\"Export Quantifications\" onclick=\"selectOption(this)\">";
}
####<GO visualisation>####
sub goVisualisation {
	return "<INPUT type=\"button\" id=\"goVisualisation\" style=\"width:190px\" value=\"Gene Ontology Summary\" onclick=\"selectOption(this)\">";
}

sub startExploratoryAnalyses {
	return "<INPUT type=\"button\" id=\"startExploratoryAnalyses\" style=\"width:190px\" value=\"Start Exploratory Analyses\" onclick=\"selectOption(this)\">";
}
#############  OBSOLETE ###############

####<Combine Analyses>####
#sub combAnaButton {
#	return "&nbsp<INPUT type=\"button\" id=\"combAna\" style=\"width:140px\" value=\"Combine Analyses\" onclick=\"selectOption(this)\">";
#}
#####<Batch Analyses>####
#sub batchAnalysesButton {
#	return "&nbsp<INPUT type=\"button\" id=\"batchAnalyses\" style=\"width:140px\" value=\"Batch Processes$batchString\" onclick=\"selectOption(this)\">";
#}
####< Clustering >#####sub clusteringButton {
#	my $title='Clustering';
#	$title.='...' if $clusteringStatus==1;
#	$title.='*' if $clusteringStatus==2;
#	return "&nbsp<INPUT type=\"button\" id=\"clustering\" style=\"width:120px;font-weight:bold\" value=\"$title\" onclick=\"selectOption(this)\">";
#}
####< Display Cluster List >#### OBSOLETE
#sub displayClusterButton {
#	return "&nbsp<INPUT type=\"button\" id=\"clusters\" style=\"width:120px\" value=\"Cluster List\" onclick=\"selectOption(this)\">";
#}
####<Process Quantification>####
#sub addQuantificationButton {
#	my ($disabled)=@_;
#	return "&nbsp<INPUT type=\"button\" id=\"addQuantification\" style=\"width:160px\" value=\"Add Quantification\" onclick=\"selectOption(this)\" $disabled>";
#}
####<Process Quantification>####
#sub launchQuantificationButton {
#	my ($disabled)=@_;
#	return "&nbsp<INPUT type=\"button\" id=\"launch\" style=\"width:160px\" value=\"Launch Quantification\" onclick=\"selectOption(this)\" $disabled>";
#}
####<Monitor Quantification>####
#sub monitorQuantificationButton {
#	return "&nbsp<INPUT type=\"button\" id=\"monitor\" style=\"width:180px\" value=\"Monitor Quantification(s)\" onclick=\"watchQuantifications()\">";
#}
####<Process Quantification>####
#sub seeRatiosButton {
#	my ($disabled)=@_;
#	return "&nbsp<INPUT type=\"button\" id=\"seeProteins\" style=\"width:180px\" value=\"See Ratios\" onclick=\"selectOption(this)\" $disabled>";
#}

####>Revision history<####
# 2.7.3 [CHANGES] Removed support for 2D-gel and spot (PP 21/01/20)
# 2.7.2 [CHANGES] Set same name for all monitor window so it does not open multiple instances of it (VS 08/01/20)
# 2.7.1 [CHANGES] Go to summary if monitor button is focused on page changing (VS 08/01/20)
# 2.7.0 [CHANGES] Use new job monitoring window opening parameters (VS 18/11/19)
# 2.6.10 [BUGFIX] Removed useless print of undef $hasChildren variable in case of Analysis (PP 13/11/19)
# 2.6.9 [MODIF] Switch monitoring to monitorJobs script (VS 21/10/19)
# 2.6.8 [BUGFIX] Removed Delete button for fully validated Analysis with linked Quantifications (PP 04/09/19)
# 2.6.7 Handles project status=-1 [no auto-end validation] (PP 07/06/18)
# 2.6.6 Add "Monitor DIA import" button (MLP 09/02/18)
# 2.6.5 Commented merge identifiers for now & restored 'Project Accessibility' for closed projects (PP 02/08/17)
# 2.6.4 Add option merge distinct identifiers (GA 31/08/16)
# 2.6.3 Removed commented code referencing manageQuantification.cgi (PP 25/07/16)
# 2.6.2 CSS style for active button (PP 15/03/16)
# 2.6.1 Add "Export quantifications" button. Requires startExploratoryAnalyses.cgi v2.0.0 (PP 08/01/16)
# 2.6.0 Left border color code to help navigation (PP 19/10/15)
# 2.5.9 Allow new report on Proteome Discoverer label-free quantification (PP 22/04/15)
# 2.5.8 No longer passes okReport param to startValidation (PP 06/03/15)
# 2.5.7 Restrict ANALYSIS deletion to bioinfo/mass/manag users only (PP 13/11/14)
# 2.5.6 move start GO and start Q GO to selectGOOptions (SL 16/10/14)
# 2.5.5 add new button exploratory analyses (SL 14/05/14)
# 2.5.4 Better control on GO option buttons (PP 04/11/13)
# 2.5.3 Now calls showProtQuantification.cgi for single analysis quantification data (PP 03/09/13)
# 2.5.2 Minor change in button display (PP 28/08/13)
# 2.5.1 Authorizing biologists to start GO analyses (FY 25/07/13)
# 2.5.0 Authorizing guests to visualize GO bar plots (FY 23/04/13)
# 2.4.9 Add Gene Ontology visualisation and Start Q. GO analysis buttons (FY 27/03/13)
# 2.4.8 Project status management (PP 03/10/12)
# 2.4.7 Single analysis can no longer be added: Use Process Analyses (PP 09/08/12)
# 2.4.6 Added "Compare Quantifications" (PP 30/07/12)
# 2.4.5 Better check for deletability due to quantification (PP 07/05/12)
# 2.4.4 Restored "Add Design" option to Experiment & code cleaning (PP 19/04/12)
# 2.4.3 Added CALL parameter for listAnaQuantifications (PP 11/04/12)
# 2.4.2 Removed check for no wget iframe in logoFrame (PP 23/03/12)
# 2.4.1 Updating the visual for quantification case (GA 15/12/11)
# 2.4.0 Added Quantifiation checks on analysis deletion (PP 24/11/11)
# 2.3.9 Adding quantification buttons (monitor, design,...) (GA 25/10/11)
# 2.3.8 Manager has full access validation data (PP 01/10/11)
# 2.3.7 Adding Ana quantification button (PP 25/07/2011)
# 2.3.6 Adding GO analysis button on experiments (FY 17/06/2011)
# 2.3.5 Super User/Admin full access to project validation data (PP 14/03/2011)
# 2.3.4 Add Comparison & Condition checks for analysis deletion (PP 14/03/2011)
# 2.3.3 new protein list management options (PP 23/02/2011)
