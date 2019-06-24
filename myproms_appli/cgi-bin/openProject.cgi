#!/usr/local/bin/perl -w

################################################################################
# openProject.cgi        1.6.9                                                 #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Generates the project's main navigation frames                               #
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

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my %analysisClass=(-1=>'no_scan',0=>'no_val',1=>'part_val',2=>'val');
my $userID=$ENV{'REMOTE_USER'};

# print header(-charset=>'UTF-8'); warningsToBrowser(1); # DEBUG
####################
####>Parameters<####
####################
my $action=(param('ACT'))? param('ACT') : 'open';
my $projectID=param('ID'); # open & nav
my $selBranchID=lc(param('branchID'));
my $selItemBranchID=(param('itemBranchID'))? lc(param('itemBranchID')) : '';
#my $experimentView=(param('VIEW'))? param('VIEW') : ''; # itemFrame
#my $isnavFrame=param('ISNAVFRAME');

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

#######################################> ACT=open <###############################################
if ($action eq 'open') {

	###########################
	####>Opening a project<####
	###########################

	###>Record date & user<###
	#my ($userName,$userStatus,$refProfile,$userEmail,$userInfo)=&promsMod::getUserInfo($dbh,$userID,$projectID);
	#if ($userStatus ne 'bioinfo') {
		my $sthUp=$dbh->prepare("UPDATE PROJECT SET LAST_OPEN=NOW(),LAST_USER=? WHERE ID_PROJECT=?");
		$sthUp->execute($userID,$projectID);
		$sthUp->finish;
		$dbh->commit;
	#}

	#####>Checking for designs in project<####
	#my ($existDesigns)=$dbh->selectrow_array("SELECT COUNT(*) FROM DESIGN WHERE ID_EXPERIMENT IN (SELECT ID_EXPERIMENT FROM EXPERIMENT WHERE ID_PROJECT=$projectID)");

	#####>Checking for gels in project<####
	#my ($existGels)=$dbh->selectrow_array("SELECT COUNT(*) FROM GEL2D WHERE ID_EXPERIMENT IN (SELECT ID_EXPERIMENT FROM EXPERIMENT WHERE ID_PROJECT=$projectID)");
	$dbh->disconnect;

	####>Starting HTML<####
	print header(-charset=>'UTF-8');
	warningsToBrowser(1);
	my $itemFrameString=($selItemBranchID)? "&itemBranchID=$selItemBranchID" : '';
	print qq
|<HTML>
<HEAD>
<TITLE>Open Project</TITLE>
<SCRIPT language="JavaScript" type="text/javascript">
//***** Default values for multiple script variables (values updated by user) *******
var selectedBranchID; // keeps track of selected item in navigation trees
var selectedAction='summary'; // selected action (used by selectOption.cgi,...)
var selectedMode='child'; // display mode for protein list (used by selectOption.cgi, listProteins.cgi, exportProteinList.cgi & others)
var selectedFilter=''; // protein list filter (used by selectOption.cgi and listProteins.cgi)
var selectedComparison='analyses'; // comparison mode ('analyses' or 'items'; used by selectOption.cgi, compareItems.cgi and compareAnalyses.cgi)
var selectedView='peptide'; // sort mode for protein list (used by selectOption.cgi, listProteins.cgi, compareItems.cgi and compareAnalyses.cgi)
var selectedPepType='ba'; // peptide count method (used by selectOption.cgi, listProteins.cgi, compareItems.cgi and compareAnalyses.cgi)
var selectedSearch='keyword'; // used by selectOption.cgi and search(Kw/Sq)Protein.cgi
var bioSampleView='name'; // default view for bioSample tree in itemFrame
//var selectedClus; // OBSOLETE record a selected cluster if any
//var expandMode=0; //*** used by listProteins.cgi & calling scripts
//var showUnClass=0; //*** used by listProteins.cgi & all calling scripts
//var graphView='StdExp'; //*** used by graphicalView.cgi & all calling scripts
//var showMSdata=0; //*** used by sequenceView.cgi & calling scripts
//var gel2dWindow; //*** window used to display 2D Gels
</SCRIPT>
</HEAD>
<FRAMESET cols="250,*" bordercolor="$darkColor">
	<FRAMESET rows="60%,*" bordercolor="$darkColor">
		<FRAME name="navFrame" scrolling=no src="./openProject.cgi?ACT=nav&ID=$projectID&branchID=$selBranchID$itemFrameString">
		<FRAME name="itemFrame">
	</FRAMESET>
	<FRAMESET rows="70,*" bordercolor="$darkColor">
		<FRAME name="optionFrame">
		<FRAME name="resultFrame">
	</FRAMESET>
</FRAMESET>
</HTML>
|;
	exit;
}

#######################################> ACT=nav (Draw navFrame)<###############################################
elsif ($action eq 'nav') {
#print header(-charset=>'UTF-8'); warningsToBrowser(1); # DEBUG

	####>Fetching user info<####
	#my ($userName,$userStatus,$refProfile,$userEmail,$userInfo)=&promsMod::getUserInfo($dbh,$userID,$projectID);
	#my $projectAccess=$$refProfile{$projectID};
	my $expandMode=0;

	###############################
	####>Building project tree<####
	###############################

	my %treeOptions;
	$treeOptions{'AUTOSAVE'}='top.projectTreeStatus'; # name of JS variable storing tree status

	####>Queries<####
	my ($projectName)=$dbh->selectrow_array("SELECT NAME FROM PROJECT WHERE ID_PROJECT=$projectID");
	my $sthExp=$dbh->prepare("SELECT ID_EXPERIMENT,NAME FROM EXPERIMENT WHERE ID_PROJECT=$projectID ORDER BY DISPLAY_POS ASC");
	my $sthG2D=$dbh->prepare("SELECT ID_GEL2D,NAME FROM GEL2D WHERE ID_EXPERIMENT=? ORDER BY DISPLAY_POS ASC");
	my $sthSamp=$dbh->prepare("SELECT ID_SAMPLE,NAME FROM SAMPLE WHERE ID_EXPERIMENT=? AND ID_SPOT IS NULL ORDER BY DISPLAY_POS ASC");
	my $sthAna=$dbh->prepare("SELECT ID_ANALYSIS,NAME,VALID_STATUS,VERIFIED_MG FROM ANALYSIS WHERE ID_SAMPLE=? ORDER BY DISPLAY_POS ASC");

	####>Project<####
	my @projectTree=(0,'project',$projectID,'','',1,$projectName,''); #(depth,type,ID,labelClass,imageClass,selectable,name,popup,refChild1,...2,...3)

	###>Experiments<###
	$sthExp->execute;
	while (my ($expID,$expName)=$sthExp->fetchrow_array) {
		my @experimentTree=(1,'experiment',$expID,'','',1,$expName,'');

		##>Gels<##
		$sthG2D->execute($expID);
		while (my ($gelID,$gelName)=$sthG2D->fetchrow_array) {
			push @experimentTree,[2,'gel2d',$gelID,'','',1,$gelName,''];
		}

		##>Free Samples<##
		$sthSamp->execute($expID);
		while (my ($sampID,$sampName)=$sthSamp->fetchrow_array) {
			my @sampleTree=(2,'sample',$sampID,'','',1,$sampName,'');

			#>Analyses<#
			$sthAna->execute($sampID);
			while (my ($anaID,$anaName,$validStatus,$verifMG)=$sthAna->fetchrow_array) {
				push @sampleTree,[3,'analysis',$anaID,'',$analysisClass{$validStatus},1,$anaName,''];
				$expandMode=1 if ($validStatus>=1 && !$verifMG); # set expandMode to 1 if at least 1 valid analysis is not verified for MG
			}
			push @experimentTree,\@sampleTree;
		}
		push @projectTree,\@experimentTree;
	}
	$sthExp->finish;
	$sthG2D->finish;
	$sthSamp->finish;
	$sthAna->finish;
	$dbh->disconnect;

	@{$treeOptions{'SELECTION'}}=($selBranchID)? split(':',$selBranchID) : ('project',$projectID);

	#######################
	####>Starting HTML<####
	#######################
	print header(-charset=>'UTF-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Navigate Project Tree</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
|;
	&promsMod::popupInfo();
	&promsMod::writeJsTreeFunctions(\@projectTree,\%treeOptions);
	print qq
|//--->Local functions<---//
var selItemBranchID='$selItemBranchID'; // modified/added spot
parent.expandMode=$expandMode;
//var view='\$experimentView'; // for itemFrame: 'quanti' or 'go'
//var isnavframe=0;
function actionOnSelect(tabIndex) {
	var selBranchID=getSelectedBranchID(); //tableArray[tabIndex][4];
	if (selBranchID.match(/project/)) {
		parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=project&ID=$projectID&VIEW="+parent.bioSampleView+"&branchID="+selBranchID; // IMPORTANT: branchID="project:projID"!
		selItemBranchID=null;
		top.itemTreeStatus='';
	}
	else if (selBranchID.match(/experiment/)) {
		parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&ID=$projectID&EXPERIMENT="+selBranchID+"&branchID="+selBranchID+"&VIEW="+top.experimentView; //+"&ISNAVFRAME=1";
		selItemBranchID=null;
		top.itemTreeStatus='';
	}
	else if (selBranchID.match(/gel2d/)) {
		parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=gel&GEL="+selBranchID+"&branchID="+selItemBranchID; //+"&ISNAVFRAME=1";
		selItemBranchID=null;
	}
	else {
		/*
		if (selBranchID.match(/design/)) {
			parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=design&ID=$projectID&DESIGN="+selBranchID+"&branchID="+selItemBranchID+"&VIEW="+top.experimentView; //+"&ISNAVFRAME=1";
			selItemBranchID=null;
		}else{
			parent.selectedBranchID=selBranchID;
			// parent.optionFrame.location="$promsPath{cgi}/selectOption.cgi?branchID="+selBranchID;
			if (parent.itemFrame) {clearItemFrame();}
		}
		*/
		parent.selectedBranchID=selBranchID;
		clearItemFrame();

	}
	parent.optionFrame.location="$promsPath{cgi}/selectOption.cgi?branchID="+selBranchID;
}
function closeProject() {
	if (top.gel2dWindow && !top.gel2dWindow.closed) {top.gel2dWindow.close();}
	parent.location='$promsPath{cgi}/selectProject.cgi'; // will clear all tree status
}
function clearItemFrame() {
	var doc=parent.itemFrame.document;
	doc.clear();
	doc.writeln('<HTML><BODY>');
	doc.writeln("</BODY></HTML>");
	doc.close();
	top.itemTreeStatus='';
}
</SCRIPT>
</HEAD>
<BODY style="margin:0px;" onload="applyTreeStatus(top.projectTreeStatus);document.getElementById(getSelectedBranchID()).scrollIntoView(0);">
<DIV class="navBorder1" style="background-color:#FFFFFF;position:absolute;top:0px;left:0px;z-index:2">
<TABLE cellspacing=0 cellpadding=2><TR><TD width=35 rowspan=2></TD>
<TH><INPUT type="button" value="Close Project" style="width:120px;font-weight:bold;" onclick="closeProject()"/></TH></TR>
<TR><TH><INPUT type="button" value="Expand" style="font-size:9px;width:60px" onclick="expandCollapseBranch('item','block')"/><INPUT type="button" value="Collapse" style="font-size:9px;width:60px" onclick="expandCollapseBranch('item','none')"/></TH></TR>
</TABLE>
</DIV>
<DIV class="navigation navBorder1">
<BR><BR><BR>
|;
	####>Printing Project Tree<####
	my $tableIndex=0;
	&promsMod::printTree(\@projectTree,\%treeOptions,\$tableIndex,1); # (refArray, is last child)

	print qq
|</DIV>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup();</SCRIPT>
</BODY>
</HTML>
|;
}

#######################################> ACT=gel (Draw gelContent)<###############################################
elsif ($action eq 'gel') {
#print header(-charset=>'UTF-8'); warningsToBrowser(1); # DEBUG

	my %itemIcones=&promsConfig::getItemIcones;
	my ($item,$itemID)=($selBranchID)? split(':',$selBranchID) : split(':',param('GEL'));
	my $gelID;
	if (param('GEL')) {(my $gelType,$gelID)=split(':',param('GEL'));}
	elsif ($item eq 'spot') {
		($gelID)=$dbh->selectrow_array("SELECT ID_GEL2D FROM SPOT WHERE ID_SPOT=$itemID");
	}
	#elsif ($item eq 'sample') {
	#	($gelID)=$dbh->selectrow_array("SELECT ID_GEL2D FROM SPOT,SAMPLE WHERE SPOT.ID_SPOT=SAMPLE.ID_SPOT AND ID_SAMPLE=$itemID");
	#}
	elsif ($item eq 'analysis') {
		($gelID)=$dbh->selectrow_array("SELECT ID_GEL2D FROM SPOT,SAMPLE,ANALYSIS WHERE SPOT.ID_SPOT=SAMPLE.ID_SPOT AND SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND ID_ANALYSIS=$itemID");
	}

	###########################
	####>Building gel tree<####
	###########################

	my %treeOptions;
	$treeOptions{'BUTTON'}='font-size:9px;width:60px';
	$treeOptions{'AUTOSAVE'}='top.itemTreeStatus'; # name of JS variable storing tree status

	####>Queries<####
	my ($gelName)=$dbh->selectrow_array("SELECT NAME FROM GEL2D WHERE ID_GEL2D=$gelID");
	my $sthSpot=$dbh->prepare("SELECT ID_SPOT,NAME FROM SPOT WHERE ID_GEL2D=? ORDER BY NAME ASC");
	my $sthSamp=$dbh->prepare("SELECT ID_SAMPLE,NAME FROM SAMPLE WHERE ID_SPOT=? ORDER BY NAME ASC");
	my $sthAna=$dbh->prepare("SELECT ID_ANALYSIS,NAME,VALID_STATUS,VERIFIED_MG FROM ANALYSIS WHERE ID_SAMPLE=? ORDER BY DISPLAY_POS ASC");

	####>Gel<####
	my @gelTree=(0,'gel2d',$gelID,'','',1,$gelName,''); #(depth,type,ID,labelClass,imageClass,selectable,name,popup,refChild1,...2,...3)

	###>Spots<###
	$sthSpot->execute($gelID);
	while (my ($spotID,$spotName,$sampID,$sampName)=$sthSpot->fetchrow_array) {
		##>Spot/Sample fusion<##
		$sthSamp->execute($spotID);
		my ($sampID,$sampName)=$sthSamp->fetchrow_array;
		my $sampPopup=($sampID)? "&nbsp;<IMG src=\\\'$promsPath{images}/$itemIcones{sample}\\\'/><B>$sampName</B>" : '<B>No sample</B>';
		my @spotTree=(1,'spot',$spotID,'','',1,$spotName,$sampPopup);

		#>Analyses<#
		if ($sampID) {
			$sthAna->execute($sampID);
			while (my ($anaID,$anaName,$validStatus,$verifMG)=$sthAna->fetchrow_array) {
				push @spotTree,[2,'analysis',$anaID,'',$analysisClass{$validStatus},1,$anaName,''];
			}
		}
		push @gelTree,\@spotTree;
	}
	$sthSpot->finish;
	$sthSamp->finish;
	$sthAna->finish;
	$dbh->disconnect;

	@{$treeOptions{'SELECTION'}}=($item,$itemID);

	#######################
	####>Starting HTML<####
	#######################
	print header(-charset=>'UTF-8');
	warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Navigate Project Tree</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
|;
	&promsMod::popupInfo();
	&promsMod::writeJsTreeFunctions(\@gelTree,\%treeOptions);
	print qq
|//--->Local functions<---//
function actionOnSelect(tabIndex) {
	var selBranchID=getSelectedBranchID();
	parent.selectedBranchID=selBranchID;
	parent.optionFrame.location="$promsPath{cgi}/selectOption.cgi?branchID="+selBranchID;
}
</SCRIPT>
</HEAD>
<BODY style="margin:0px;" onload="if (top.itemTreeStatus) {applyTreeStatus(top.itemTreeStatus)} else {click2OpenClose(0,'open',1)};document.getElementById(getSelectedBranchID()).scrollIntoView(0);">
<DIV class="navigation navBorder2">
|;

	####>Printing Gel Tree<####
	my $tableIndex=0;
	&promsMod::printTree(\@gelTree,\%treeOptions,\$tableIndex,1); # (refArray, is last child)

	print qq
|</DIV>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup()</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

#######################################>ACT=experiment<###############################################
###> Print the button that allows the choice between GO-mode and QUANTI-mode
elsif ($action eq 'experiment') {# If experiment -> 2 possibilities: GO ANALYSIS or Design creation
	my $experimentView=(param('VIEW'))? param('VIEW') : ''; # itemFrame
	my $expID=(param('EXPERIMENT'))? (split(':',param('EXPERIMENT')))[1] : undef;
	my ($item,$itemID)=split(':',$selBranchID);
	$experimentView = ($experimentView eq 'go')? 'functAna' : $experimentView;
#	my $itemFrameString=($selItemBranchID)? "&itemBranchID=$selItemBranchID" : '';
#
#	#######################
#	####>Starting HTML<####
#	#######################
#	print header(-charset=>'UTF-8');
#	warningsToBrowser(1);
#	print qq
#|<HTML>
#<HEAD>
#<TITLE>Navigate Project Tree</TITLE>
#<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
#<SCRIPT LANGUAGE="JavaScript">
#function updateView(newView) {
#	parent.navFrame.view=newView;
#/*
#	var selBranchID=parent.navFrame.getSelectedBranchID();
#	if (parent.navFrame.view.match(/GO/)) {
#		parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=go&ID=\$projectID&EXPERIMENT="+selBranchID+"&branchID=$selBranchID$itemFrameString&VIEW="+newView; //+"&ISNAVFRAME=0";
#	}else if (parent.navFrame.view.match(/QUANTI/)) {
#		parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=design&ID=\$projectID&EXPERIMENT="+selBranchID+"&branchID=$selBranchID$itemFrameString&VIEW="+newView; //+"&ISNAVFRAME=0";
#	}else {
#		parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&ID=\$projectID&EXPERIMENT="+selBranchID+"&branchID=$selBranchID$itemFrameString&VIEW="+newView; //+"&ISNAVFRAME=0";
#	}
#*/
#	window.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=$selBranchID&branchID=$selBranchID&VIEW="+newView;
#}
#</SCRIPT>
#</HEAD>
#<BODY>
#<CENTER>
#<FONT style="font-size:9px;font-weight:bold;">Display :</FONT><SELECT name="type" style="font-size:9px;font-weight:bold;" onchange="updateView(this.value)">
#|;
#	###> Print the select options according to what was chosen by the user
#	print "<OPTION value=\"\">-= Select =-</OPTION>\n" unless $experimentView;
#	my ($selGo,$selQuanti)=(!$experimentView)? ('','') : ($experimentView eq 'go')? (' selected','') : ('',' selected');
#	print qq
#|<OPTION value="go"$selGo>GO Analyses</OPTION>
#<OPTION value="quanti"$selQuanti>Quantifications</OPTION>
#</SELECT>
#</CENTER>
#|;

	########################> No view <#############################
	if (!$experimentView) {
		$dbh->disconnect;
		print header(-charset=>'UTF-8');
		warningsToBrowser(1);
		print qq
|<HTML>
<HEAD>
<TITLE>Navigate Experiment Tree</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY style="margin:0px;">
<DIV class="navigation navBorder2">
|;
		####>View selection menu<####
		&printExperimentViewSelection($experimentView);
		print qq
|<BR>
<FONT class="title3">&nbsp;No view selected.</FONT>
</DIV>
</BODY>
</HTML>
|;
		exit;
	}

	########################> VIEW=quanti <#############################
	#elsif ($action eq 'design' || ($experimentView eq 'QUANTI'  && $action eq 'experiment')) {#}
	elsif ($experimentView eq 'quanti') {
		#my ($item,$itemID)=split(':',param('EXPERIMENT'));
		unless ($expID) {
			if ($item eq 'experiment') {
				$expID = $itemID;
			}
			elsif ($item eq 'design') {
				($expID)=$dbh->selectrow_array("SELECT ID_EXPERIMENT FROM DESIGN WHERE ID_DESIGN=$itemID");
			}
			elsif ($item eq 'quantification') {
				my ($anaID)=$dbh->selectrow_array("SELECT ID_ANALYSIS FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$itemID LIMIT 0,1");
				($expID)=$dbh->selectrow_array("SELECT S.ID_EXPERIMENT FROM ANALYSIS A, SAMPLE S WHERE A.ID_SAMPLE=S.ID_SAMPLE AND A.ID_ANALYSIS=$anaID"); # as long as analyses are directly linked to sample!
			}
		}
		#my ($nbDesign)=$dbh->selectrow_array("SELECT COUNT(ID_DESIGN) FROM DESIGN WHERE ID_EXPERIMENT=$expID");
		my ($expName)=$dbh->selectrow_array("SELECT NAME FROM EXPERIMENT WHERE ID_EXPERIMENT=$expID");


		###########################
		####>Building gel tree<####
		###########################

		my %treeOptions;
		$treeOptions{'BUTTON'}='font-size:9px;width:60px';
		$treeOptions{'AUTOSAVE'}='top.itemTreeStatus'; # name of JS variable storing tree status

		my @experimentTree=(0,'experiment',$expID,'','',1,$expName,'');
		#if($nbDesign>0){# No design has been created yet for this experiment
		my $existDesign=0;
		my $sthQM=$dbh->prepare("SELECT ID_QUANTIFICATION_METHOD,CODE FROM QUANTIFICATION_METHOD");
		my $sthD=$dbh->prepare("SELECT ID_DESIGN,NAME FROM DESIGN WHERE ID_EXPERIMENT=$expID ORDER BY NAME ASC");
		#my $sthQuanti=$dbh->prepare("SELECT ID_QUANTIFICATION,NAME FROM QUANTIFICATION WHERE ID_DESIGN=? AND ID_QUANTIFICATION NOT IN (SELECT DISTINCT(ID_PARENT_QUANTIFICATION) FROM PARENT_QUANTIFICATION) ORDER BY NAME ASC");
		my $sthDQ=$dbh->prepare('SELECT ID_QUANTIFICATION,NAME,FOCUS,QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_DESIGN=?'); # quanti with design
		my $sthQP=$dbh->prepare('SELECT P.ID_PARENT_QUANTIFICATION,NAME,FOCUS,QUANTIF_ANNOT,ID_QUANTIFICATION_METHOD FROM PARENT_QUANTIFICATION P,QUANTIFICATION Q WHERE P.ID_PARENT_QUANTIFICATION=Q.ID_QUANTIFICATION AND P.ID_QUANTIFICATION=?'); # quanti parents
		my $sthSA=$dbh->prepare("SELECT S.ID_SAMPLE,S.NAME,S.DISPLAY_POS,A.ID_ANALYSIS,A.NAME,A.DISPLAY_POS FROM ANA_QUANTIFICATION AQ,ANALYSIS A,SAMPLE S WHERE AQ.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND AQ.ID_QUANTIFICATION=?");
		my $sthNDPQ=$dbh->prepare("SELECT DISTINCT Q.ID_QUANTIFICATION,Q.NAME,FOCUS,QUANTIF_ANNOT,ID_QUANTIFICATION_METHOD FROM QUANTIFICATION Q, ANA_QUANTIFICATION AQ, SAMPLE S, ANALYSIS A WHERE S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS=AQ.ID_ANALYSIS AND Q.ID_QUANTIFICATION=AQ.ID_QUANTIFICATION AND S.ID_EXPERIMENT=$expID AND FOCUS='protein' AND ID_DESIGN IS NULL");
		my $sthPQ=$dbh->prepare("SELECT DISTINCT Q.ID_QUANTIFICATION,Q.NAME,FOCUS,QUANTIF_ANNOT,ID_QUANTIFICATION_METHOD FROM QUANTIFICATION Q, ANA_QUANTIFICATION AQ, SAMPLE S, ANALYSIS A WHERE S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS=AQ.ID_ANALYSIS AND Q.ID_QUANTIFICATION=AQ.ID_QUANTIFICATION AND S.ID_EXPERIMENT=$expID AND FOCUS='peptide'"); # peptide quantif AND QUANTIF_ANNOT LIKE 'ALGO_TYPE=%'"); # MassChroq only 'LABEL=FREE%'
		#my $sthQuantiParentID=$dbh->prepare("SELECT ID_PARENT_QUANTIFICATION FROM PARENT_QUANTIFICATION WHERE ID_QUANTIFICATION=?");
		#my $sthQuantiParentName=$dbh->prepare("SELECT NAME FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");

		my (%quantifMethods,%fetchedQuantifications);
		$sthQM->execute;
		while (my ($qMethID,$code)=$sthQM->fetchrow_array) {$quantifMethods{$qMethID}=$code;}
		$sthQM->finish;

		####>Design<####
		$sthD->execute;
		while (my ($designID,$designName) = $sthD->fetchrow_array) {
			$existDesign=1;
			my @designTree=(1,'design',$designID,'','',1,$designName,''); #(depth,type,ID,labelClass,imageClass,selectable,name,popup,refChild1,...2,...3)
			###>Quantifications with design<###
			$sthDQ->execute($designID);
			my %quantifs;
			while (my ($quantiID,$quantiName,$focus,$annot)=$sthDQ->fetchrow_array) {
				@{$quantifs{$quantiID}}=($quantiName,$focus,$annot); # pre-storing because sortSmart required
			}
			foreach my $quantiID (sort{&promsMod::sortSmart(lc($quantifs{$a}[0]),lc($quantifs{$b}[0]))} keys %quantifs) {
				my ($label)=($quantifs{$quantiID}[2] && $quantifs{$quantiID}[2]=~/^LABEL=([^:]+)/)? $1 : 'FREE';
				my $isLabeled=($label eq 'FREE')? '_no_label' : '_label';
				##>Quantification<##
				my @quantiTree=(2,'quantification',$quantiID,'',"$quantifs{$quantiID}[1]$isLabeled",1,$quantifs{$quantiID}[0],'');
				$sthQP->execute($quantiID);
				my $refQuantiDesign=$sthQP->fetchall_arrayref; # reference to an array
				&createPeptideQuantificationBranches($refQuantiDesign,$sthSA,$quantiID,\@quantiTree,'design',\%fetchedQuantifications,\%quantifMethods);
				push @designTree,\@quantiTree;
			}
			push @experimentTree,\@designTree;
		}
		push @experimentTree,[1,'none',0,'','',0,'No Design',''] unless $existDesign;

		####>Add non-design protein quantif<####
		my @noDesignTree=(1,'none',0,'','',0,'Other protein quantifications','');
		$sthNDPQ->execute;
		my $refQuantiNoDesign=$sthNDPQ->fetchall_arrayref;
		&createPeptideQuantificationBranches($refQuantiNoDesign,$sthSA,undef,\@noDesignTree,'no-design',\%fetchedQuantifications,\%quantifMethods);
		push @experimentTree,\@noDesignTree if $noDesignTree[8]; # has children

		####> Add Quantifications multiple 30/09/2014
		my @xicTree=(1,'none',0,'','',0,'Unused peptide quantifications',''); #(depth,type,ID,labelClass,imageClass,selectable,name,popup,refChild1,...2,...3)
		$sthPQ->execute;
		my $refPepQuantif=$sthPQ->fetchall_arrayref;
		&createPeptideQuantificationBranches($refPepQuantif,$sthSA,undef,\@xicTree,'peptide',\%fetchedQuantifications,\%quantifMethods);
		push @experimentTree,\@xicTree if $xicTree[8];

		$sthD->finish;
		$sthDQ->finish;
		$sthQP->finish;
		$sthSA->finish;
		$sthNDPQ->finish;
		$sthPQ->finish;
		#}
		$dbh->disconnect;
		#$item=lc($item);
		###> To select the good item in the navigation Tree
		#@{$treeOptions{'SELECTION'}}=(param('branchID') && lc(param('branchID')) !~ /expcondition/ )? split(':',lc(param('branchID'))): (lc(param('branchID')) !~ /expcondition/ )? split(':',lc(param('branchID'))) : ('experiment',$itemID);
		@{$treeOptions{'SELECTION'}}=($selBranchID)? split(':',$selBranchID) : undef;

		#######################
		####>Starting HTML<####
		#######################
		print header(-charset=>'UTF-8');
		warningsToBrowser(1);
		print qq
|<HTML>
<HEAD>
<TITLE>Navigate Project Tree</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
|;
		&promsMod::popupInfo();
		&promsMod::writeJsTreeFunctions(\@experimentTree,\%treeOptions);
		print qq
|//--->Local functions<---//
//parent.navFrame.isnavframe=\$isnavFrame;
function actionOnSelect(tabIndex) {
	var selBranchID=getSelectedBranchID();
	/*
	if (parent.navFrame.isnavframe==0) {
		parent.optionFrame.location="$promsPath{cgi}/selectOptionQuanti.cgi?branchID="+selBranchID;
	}else{
		parent.navFrame.isnavframe=0;
	}
	*/
	if (selBranchID.match(/experiment/)) {
		if (parent.selectedBranchID != selBranchID) {parent.optionFrame.location="$promsPath{cgi}/selectOption.cgi?branchID="+selBranchID;} // do not load Exp options twice
	}
	else {
		parent.optionFrame.location="$promsPath{cgi}/selectOptionQuanti.cgi?branchID="+selBranchID;
	}
	parent.selectedBranchID=selBranchID;
}
</SCRIPT>
</HEAD>
<BODY style="margin:0px;" onload="if (top.itemTreeStatus) {applyTreeStatus(top.itemTreeStatus)} else {click2OpenClose(0,'open',1)};document.getElementById(getSelectedBranchID()).scrollIntoView(0);">
<DIV class="navigation navBorder2">
|;
		####>View selection menu<####
		&printExperimentViewSelection($experimentView);

		####>Printing Design Tree<####
		my $tableIndex=0;
		&promsMod::printTree(\@experimentTree,\%treeOptions,\$tableIndex,1); # (refArray, is last child)

		print qq
|</DIV>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">
	setPopup();
</SCRIPT>
</BODY>
</HTML>
|;
		exit;
	}


	########################> VIEW=go <#############################
	#elsif ($experimentView eq 'go') {
	elsif ($experimentView eq 'functAna') {
		unless ($expID) {
			if ($item eq 'experiment') {
				$expID = $itemID;
			}
			elsif ($item eq 'go_analysis'){
				($expID) = $dbh->selectrow_array("SELECT ID_EXPERIMENT FROM GO_ANALYSIS WHERE ID_GOANALYSIS=$itemID");
			}
			elsif ($item eq 'pathwayanalysis') {
				($expID) = $dbh->selectrow_array("SELECT ID_EXPERIMENT FROM PATHWAY_ANALYSIS WHERE ID_PATHWAY_ANALYSIS=$itemID");
			}
			elsif ($item eq 'motifanalysis') {
				($expID) = $dbh->selectrow_array("SELECT ID_EXPERIMENT FROM EXPLORANALYSIS WHERE ID_EXPLORANALYSIS=$itemID");
			}
		}
		my ($expName) = $dbh->selectrow_array("SELECT NAME FROM EXPERIMENT WHERE ID_EXPERIMENT=$expID");
		my $sthGOAna=$dbh->prepare("SELECT ID_GOANALYSIS, NAME FROM GO_ANALYSIS WHERE ID_EXPERIMENT=$expID AND ID_PARENT_GOANA IS NULL ORDER BY NAME ASC");
		my $sthPathwayAna=$dbh->prepare("SELECT ID_PATHWAY_ANALYSIS, NAME FROM PATHWAY_ANALYSIS WHERE ID_EXPERIMENT=$expID ORDER BY NAME ASC");
		my $sthMotifAna=$dbh->prepare("SELECT ID_EXPLORANALYSIS, NAME, ANA_TYPE FROM EXPLORANALYSIS WHERE ID_EXPERIMENT=$expID AND ANA_TYPE like '%MOTIF' ORDER BY NAME ASC");
		#my $sthHMmotifAna=$dbh->prepare("SELECT ID_EXPLORANALYSIS, NAME FROM EXPLORANALYSIS WHERE ID_EXPERIMENT=$expID AND ANA_TYPE='HM_MOTIF' ORDER BY NAME ASC");

		################################
		####>Building func ana tree<####
		################################
		my %treeOptions;
		$treeOptions{'BUTTON'}='font-size:9px;width:60px';
		$treeOptions{'AUTOSAVE'}='top.itemTreeStatus'; # name of JS variable storing tree status

		my @expTree=(0,'experiment',$expID,'','',1,$expName,''); #(depth,type,ID,labelClass,imageClass,selectable,name,popup,refChild1,...2,...3)
		my @goTree=(1,'go',$expID,'','',1,'GO Analyses','');
		#my $existGoAna=0;
		$sthGOAna->execute;
		while(my ($goAnaID,$goAnaName) = $sthGOAna->fetchrow_array){
			#$existGoAna=1;
			my @goAnaTree=(2,'go_analysis',$goAnaID,'','',1,$goAnaName,'');
			push @goTree,\@goAnaTree;
		}
		#push @goTree,[2,'none',0,'','',0,'None',''] unless $existGoAna;
		push @expTree,\@goTree;

		my @pathTree=(1,'pathway',$expID,'','',1,'Pathway Analyses','');
		#my $existPathAna=0;
		$sthPathwayAna->execute;
		while (my($pathwayID,$pathwayName)=$sthPathwayAna->fetchrow_array) {
			#$existPathAna=1;
			my @pathAnaTree=(2,'pathwayanalysis',$pathwayID,'','',1,$pathwayName,'');
			push @pathTree,\@pathAnaTree;
		}
		#push @pathTree,[2,'none',0,'','',0,'None',''] unless $existPathAna;
		push @expTree,\@pathTree;

		my @motifTree=(1,'motif',$expID,'','',1,'Motif Analyses','');
		$sthMotifAna->execute;
		while(my($explorAnaID, $explorAnaName, $anaType)=$sthMotifAna->fetchrow_array) {
			my @motifAnaTree=(2,'motifanalysis',$explorAnaID,'',lc($anaType),1,$explorAnaName);
			push @motifTree,\@motifAnaTree;
		}

		push @expTree, \@motifTree;

		$sthGOAna->finish;
		$sthPathwayAna->finish;
		$sthMotifAna->finish;
		$dbh->disconnect;
		@{$treeOptions{'SELECTION'}}=($selBranchID)? split(':',$selBranchID) : undef;

		#######################
		####>Starting HTML<####
		#######################
		print header(-charset=>'UTF-8');
		warningsToBrowser(1);
		print qq
|<HTML>
<HEAD>
<TITLE>Navigate GO Analysis Tree</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
|;
		&promsMod::popupInfo();
		&promsMod::writeJsTreeFunctions(\@expTree,\%treeOptions);
		print qq
|//--->Local functions<---//
function actionOnSelect(tabIndex) {
	var selBranchID=getSelectedBranchID();
	if (selBranchID.match(/experiment/)) {
		if (parent.selectedBranchID != selBranchID) {parent.optionFrame.location="$promsPath{cgi}/selectOption.cgi?branchID="+selBranchID;} // do not load Exp options twice
	}
	else if (selBranchID.match(/go/)) {
		parent.optionFrame.location="$promsPath{cgi}/selectGOOption.cgi?ID="+selBranchID;
	}
	else if (selBranchID.match(/pathway/)){
		parent.optionFrame.location="$promsPath{cgi}/selectOptionPathway.cgi?ID="+selBranchID;
	}
	else {
		parent.optionFrame.location="$promsPath{cgi}/selectOptionMotif.cgi?ID="+selBranchID;
	}
	parent.selectedBranchID=selBranchID;
}
</SCRIPT>
</HEAD>
<BODY style="margin:0px;" onload="if (top.itemTreeStatus) {applyTreeStatus(top.itemTreeStatus)} else {click2OpenClose(0,'open',1)};document.getElementById(getSelectedBranchID()).scrollIntoView(0);">
<DIV class="navigation navBorder2">
|;
		####>View selection menu<####
		&printExperimentViewSelection($experimentView);

		####>Printing GO Tree<####
		my $tableIndex=0;
		&promsMod::printTree(\@expTree,\%treeOptions,\$tableIndex,1); # (refArray, is last child)

		print qq
|</DIV>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">
setPopup();
</SCRIPT>
</BODY>
</HTML>
|;
		exit;
	}
	elsif ($experimentView eq 'explorAna') {

		unless ($expID) {
			if ($item eq 'experiment') {
				$expID = $itemID;
			}
		}
		my ($expName) = $dbh->selectrow_array("SELECT NAME FROM EXPERIMENT WHERE ID_EXPERIMENT=$expID");

		my %treeOptions;
		$treeOptions{'BUTTON'}='font-size:9px;width:60px';
		$treeOptions{'AUTOSAVE'}='top.itemTreeStatus'; # name of JS variable storing tree status

		my @expTree=(0,'experiment',$expID,'','',1,$expName,'');
		#my $existDesign=0;
		#my $sthD=$dbh->prepare("SELECT ID_DESIGN,NAME FROM DESIGN WHERE ID_EXPERIMENT=$expID ORDER BY NAME ASC");
		#$sthD->execute;
		#while (my ($designID,$designName) = $sthD->fetchrow_array) {
		#	$existDesign=1;
		#	my @designTree=(1,'design',$designID,'','',1,$designName,'');
		#	push @expTree,\@designTree;
		#}
		#
		#push @expTree,[1,'none',0,'','',0,'None',''] unless $existDesign;
		#$sthD->finish;

		my $existExplorAna = 0;
		my $sthEA = $dbh -> prepare("SELECT ID_EXPLORANALYSIS, NAME, ANA_TYPE from EXPLORANALYSIS where ID_EXPERIMENT = $expID and ANA_TYPE != 'MOTIF' AND ANA_TYPE != 'HM_MOTIF' ORDER BY ANA_TYPE ASC,NAME ASC");
		$sthEA -> execute;
		while (my($explorAnaID, $explorAnaName, $anaType) = $sthEA -> fetchrow_array) {
			$existExplorAna = 1;
			my @explorAnaTree = (1,'exploranalysis',$explorAnaID,'',lc($anaType),1,$explorAnaName);
			push @expTree,\@explorAnaTree;
		}
		$sthEA -> finish;
		$dbh->disconnect;

		@{$treeOptions{'SELECTION'}}=($selBranchID)? split(':',$selBranchID) : undef;

		#######################
		####>Starting HTML<####
		#######################
		print header(-charset=>'UTF-8');
		warningsToBrowser(1);
		print qq
|<HTML>
<HEAD>
<TITLE>Navigate Project Tree</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
|;
		&promsMod::popupInfo();
		&promsMod::writeJsTreeFunctions(\@expTree,\%treeOptions);
		print qq
|//--->Local functions<---//
//parent.navFrame.isnavframe=\$isnavFrame;
function actionOnSelect(tabIndex) {
	var selBranchID=getSelectedBranchID();
	if (selBranchID.match(/experiment/)) {
		if (parent.selectedBranchID != selBranchID) {parent.optionFrame.location="$promsPath{cgi}/selectOption.cgi?branchID="+selBranchID;} // do not load Exp options twice
	}
	else if (selBranchID.match(/design/)) {
		parent.optionFrame.location="$promsPath{cgi}/selectOptionQuanti.cgi?fromExplorAna=1&branchID="+selBranchID;
	}
	else {
		parent.optionFrame.location="$promsPath{cgi}/selectOptionExplorAna.cgi?branchID="+selBranchID;
	}
	parent.selectedBranchID=selBranchID;
}
</SCRIPT>
</HEAD>
<BODY style="margin:0px;" onload="if (top.itemTreeStatus) {applyTreeStatus(top.itemTreeStatus)} else {click2OpenClose(0,'open',1)};document.getElementById(getSelectedBranchID()).scrollIntoView(0);">
<DIV class="navigation navBorder2">
|;
		####>View selection menu<####
		&printExperimentViewSelection($experimentView);

		####>Printing Design Tree<####
		my $tableIndex=0;
		&promsMod::printTree(\@expTree,\%treeOptions,\$tableIndex,1); # (refArray, is last child)

		print qq
|</DIV>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
</BODY>
</HTML>
|;
		exit;
	}
}

#######################################>ACT=project<###############################################
elsif ($action eq 'project') {
	my $sampleView = param('VIEW')? param('VIEW') : "name";
	my %treeOptions;
	$treeOptions{'BUTTON'}='font-size:9px;width:60px';
	$treeOptions{'AUTOSAVE'}='top.itemTreeStatus';

	@{$treeOptions{'SELECTION'}}=($selBranchID)? split(':',$selBranchID) : undef;
	
	my $sthUsedBioSample=$dbh->prepare("SELECT 'used' FROM OBSERVATION WHERE ID_BIOSAMPLE=? LIMIT 1");

	my $existBiosample = 0;
 	my @bioSampTree=(0,'all',$projectID,'','',1,'All BioSamples','');
	if ($sampleView eq 'name') {
		my $sthSelBioSample = $dbh->prepare("SELECT BS.ID_BIOSAMPLE, BS.NAME FROM BIOSAMPLE BS
												INNER JOIN PROJECT_BIOSAMPLE PB ON BS.ID_BIOSAMPLE = PB.ID_BIOSAMPLE
												WHERE PB.ID_PROJECT = $projectID");
		$sthSelBioSample->execute;
		my %bioSamples;
		while (my ($bioSampID, $bioSampName) = $sthSelBioSample->fetchrow_array) {
			$existBiosample=1;
			$sthUsedBioSample->execute($bioSampID);
			my ($usedStatus)=$sthUsedBioSample->fetchrow_array;
			$usedStatus='not_used' unless $usedStatus;
			@{$bioSamples{$bioSampID}}=($bioSampName,$usedStatus);
		}
		$sthSelBioSample->finish;

		foreach my $bioSampID (sort{&promsMod::sortSmart(lc($bioSamples{$a}[0]),lc($bioSamples{$b}[0]))} keys %bioSamples) {
			push @bioSampTree,[1,'biosample',$bioSampID,'',$bioSamples{$bioSampID}[1],1,$bioSamples{$bioSampID}[0],''];
		}
	}
	elsif ($sampleView eq 'species') {
		my $sthSelBioSample = $dbh -> prepare("SELECT BS.ID_BIOSAMPLE,BS.NAME,S.ID_SPECIES,S.SCIENTIFIC_NAME,S.COMMON_NAME FROM BIOSAMPLE BS
												INNER JOIN SPECIES S ON BS.ID_SPECIES = S.ID_SPECIES
												INNER JOIN PROJECT_BIOSAMPLE PB ON BS.ID_BIOSAMPLE = PB.ID_BIOSAMPLE
												WHERE PB.ID_PROJECT = $projectID");
		$sthSelBioSample -> execute;
		my (%bioSamples,%species);
		while (my ($bioSampID,$bioSampName,$speciesID,$scientificName,$commonName) = $sthSelBioSample -> fetchrow_array) {
			$existBiosample=1;
			$sthUsedBioSample->execute($bioSampID);
			my ($usedStatus)=$sthUsedBioSample->fetchrow_array;
			$usedStatus='not_used' unless $usedStatus;
			@{$bioSamples{$speciesID}{$bioSampID}}=($bioSampName,$usedStatus);
			@{$species{$speciesID}}=($scientificName,$commonName) unless $species{$speciesID};
		}
		$sthSelBioSample->finish;
		foreach my $speciesID (sort{lc($species{$a}[0]) cmp lc($species{$b}[0])} keys %species) {
			my @speciesTree=(1,'species',$speciesID,'','',0,"$species{$speciesID}[0] ($species{$speciesID}[1])",'');
			foreach my $bioSampID (sort{&promsMod::sortSmart(lc($bioSamples{$speciesID}{$a}[0]),lc($bioSamples{$speciesID}{$b}[0]))} keys %{$bioSamples{$speciesID}}) {
				push @speciesTree,[2,'biosample',$bioSampID,'',$bioSamples{$speciesID}{$bioSampID}[1],1,$bioSamples{$speciesID}{$bioSampID}[0],''];
			}
			push @bioSampTree,\@speciesTree;
		}
	}
	elsif ($sampleView eq 'treatments') {
		my $sthSelBioSample = $dbh -> prepare("SELECT BS.ID_BIOSAMPLE,BS.NAME,BP.ID_PROPERTY
												FROM BIOSAMPLE BS
												INNER JOIN PROJECT_BIOSAMPLE PB ON BS.ID_BIOSAMPLE = PB.ID_BIOSAMPLE
												LEFT JOIN BIOSAMPLE_PROPERTY BP ON BS.ID_BIOSAMPLE=BP.ID_BIOSAMPLE
												WHERE PB.ID_PROJECT = $projectID");
		my $sthSelProp = $dbh -> prepare("SELECT NAME,PROPERTY_TYPE FROM PROPERTY WHERE ID_PROPERTY=?");
		my (%properties,%bioSampleName,%treatBioSamples);
		$properties{0}=['','T'];
		$treatBioSamples{0}={};
		$sthSelBioSample -> execute;
		while (my ($bioSampID,$bioSampName,$propID) = $sthSelBioSample -> fetchrow_array) {
			$existBiosample=1;
			$sthUsedBioSample->execute($bioSampID);
			my ($usedStatus)=$sthUsedBioSample->fetchrow_array;
			$usedStatus='not_used' unless $usedStatus;
			@{$bioSampleName{$bioSampID}}=($bioSampName,$usedStatus);
			if ($propID) {
				unless ($properties{$propID}) {
					$sthSelProp->execute($propID);
					@{$properties{$propID}}=$sthSelProp->fetchrow_array;
				}
				if ($properties{$propID}[1] eq 'T') {$treatBioSamples{$propID}{$bioSampID}=1;}
				else {$treatBioSamples{0}{$bioSampID}=1;}
			}
			else {$treatBioSamples{0}{$bioSampID}=1;}
		}
		$sthSelBioSample->finish;
		$sthSelProp->finish;

		my $selBioSampID=($selBranchID=~/biosample:(\d+)/)? $1 : 0; # in case $selBranchID = biosample:<bioSampID>:<treatID>
		foreach my $treatID (sort{lc($properties{$a}[0]) cmp lc($properties{$b}[0])} keys %treatBioSamples) {
			next if $treatID==0;
			my @treatmentTree=(1,'treatment',$treatID,'','',0,$properties{$treatID}[0],'');
			foreach my $bioSampID (sort{&promsMod::sortSmart(lc($bioSampleName{$a}[0]),lc($bioSampleName{$b}[0]))} keys %{$treatBioSamples{$treatID}}) {
				push @treatmentTree,[2,'biosample',"$bioSampID:$treatID",'',$bioSampleName{$bioSampID}[1],1,$bioSampleName{$bioSampID}[0],''];
				delete $treatBioSamples{0}{$bioSampID} if $treatBioSamples{0}{$bioSampID}; # IMPORTANT:sample can be "untreated" for another treatment!
				if ($bioSampID==$selBioSampID) { # matches 1st occurence of $bioSampID in tree
					@{$treeOptions{'SELECTION'}}=('biosample',"$bioSampID:$treatID"); # overwrites previous init
					$selBioSampID=0;
				}
			}
			push @bioSampTree,\@treatmentTree;
		}
		if (scalar keys %{$treatBioSamples{0}}) {
			my @treatmentTree=(1,'treatment',0,'','',0,'*Untreated*','');
			foreach my $bioSampID (sort{lc($bioSampleName{$a}[0]) cmp lc($bioSampleName{$b}[0])} keys %{$treatBioSamples{0}}) {
				push @treatmentTree,[2,'biosample',"$bioSampID:0",'',$bioSampleName{$bioSampID}[1],1,$bioSampleName{$bioSampID}[0],''];
				if ($bioSampID==$selBioSampID) { # matches 1st occurence of $bioSampID in tree
					@{$treeOptions{'SELECTION'}}=('biosample',"$bioSampID:0"); # overwrites previous init
					$selBioSampID=0;
				}
			}
			push @bioSampTree,\@treatmentTree;
		}
	}
	elsif ($sampleView=~/^property:(.+)/) {
		my $propertyID=$1;
		my $sthSelBioSample = $dbh -> prepare("SELECT BS.ID_BIOSAMPLE,BS.NAME,BP.PROPERTY_VALUE FROM BIOSAMPLE BS
												LEFT JOIN BIOSAMPLE_PROPERTY BP ON BS.ID_BIOSAMPLE=BP.ID_BIOSAMPLE AND BP.ID_PROPERTY=$propertyID
												INNER JOIN PROJECT_BIOSAMPLE PB ON BS.ID_BIOSAMPLE = PB.ID_BIOSAMPLE
												WHERE PB.ID_PROJECT = $projectID");
		my %bioSamples;
		$sthSelBioSample -> execute;
		while (my ($bioSampID, $bioSampName, $propValue) = $sthSelBioSample -> fetchrow_array) {
			$existBiosample=1;
			$propValue=' ' unless $propValue; # too make sure it comme first in sort
			$sthUsedBioSample->execute($bioSampID);
			my ($usedStatus)=$sthUsedBioSample->fetchrow_array;
			$usedStatus='not_used' unless $usedStatus;
			@{$bioSamples{$propValue}{$bioSampID}}=($bioSampName,$usedStatus);
		}
		$sthSelBioSample->finish;
		my $valueRank=0;
		foreach my $propValue (sort{&promsMod::sortSmart(lc($a),lc($b))} keys %bioSamples) {
			$valueRank++;
			my $dispValue=($propValue eq ' ')? '*No value*' : $propValue;
			my @propertyTree=(1,'property',$valueRank,'','',0,$dispValue,'');
			foreach my $bioSampID (sort{&promsMod::sortSmart(lc($bioSamples{$propValue}{$a}[0]),lc($bioSamples{$propValue}{$b}[0]))} keys %{$bioSamples{$propValue}}) {
				push @propertyTree,[2,'biosample',$bioSampID,'',$bioSamples{$propValue}{$bioSampID}[1],1,$bioSamples{$propValue}{$bioSampID}[0],''];
			}
			push @bioSampTree,\@propertyTree;
		}
	}

	$sthUsedBioSample->finish;
	
	push @bioSampTree,[1,'none',$projectID,'','',0,'None',''] unless $existBiosample;

	#######################
	####>Starting HTML<####
	#######################

	print header(-charset=>'UTF-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Navigate Biosample Tree</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
|;
	&promsMod::popupInfo();
	&promsMod::writeJsTreeFunctions(\@bioSampTree,\%treeOptions);
	print qq
|var firstLoad = true;
function actionOnSelect(tabIndex) {
	var selBranchID=getSelectedBranchID();
	if (firstLoad == true && '$selBranchID'.match('project:') ) { // just called by project
		firstLoad = false;
		return;
	}
	else if (selBranchID.match('all')) {
		parent.optionFrame.location="$promsPath{cgi}/selectOptionBioSample.cgi?branchID="+selBranchID;
	}
	else if (selBranchID.match('biosample') ) {
		if (top.promsFrame.selectedAction == 'properties' \|\| top.promsFrame.selectedAction == 'treatments') {
			top.promsFrame.selectedAction = 'summary';
		}
		parent.optionFrame.location="$promsPath{cgi}/selectOptionBioSample.cgi?projectID=$projectID&branchID="+selBranchID;
	}

	parent.selectedBranchID=selBranchID;
}
</SCRIPT>
</HEAD>
<BODY style="margin:0px;" onload="if (top.itemTreeStatus) {applyTreeStatus(top.itemTreeStatus)} else {click2OpenClose(0,'open',1)};document.getElementById(getSelectedBranchID()).scrollIntoView(0);">
<DIV class="navigation navBorder2">
|;
		####>View selection menu<####
		&printSampleViewSelection($dbh,$sampleView);

		$dbh -> disconnect;

		####>Printing Design Tree<####
		my $tableIndex=0;
		&promsMod::printTree(\@bioSampTree,\%treeOptions,\$tableIndex,1); # (refArray, is last child)
print qq
|</DIV>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
</BODY>
</HTML>
|;

}

sub printExperimentViewSelection {
	my ($expView)=@_;
	print qq
|<SCRIPT LANGUAGE="JavaScript">
function updateView(newView) {
	top.itemTreeStatus=''; // reset item tree
	top.experimentView=newView;
	window.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT="+parent.navFrame.getSelectedBranchID()+"&branchID=$selBranchID&VIEW="+newView;
}
</SCRIPT>
<TABLE align="center" bgcolor="$darkColor"><TR>
	<TH align=right nowrap>&nbsp;View :</TH>
	<TD><SELECT name="view" id="view" style="font-weight:bold;" onchange="updateView(this.value)">
|; # id is used by watchQuantification.cgi
	###> Print the select options according to what was chosen by the user
	print "<OPTION value=\"\">-= Select =-</OPTION>\n" unless $expView;
	my ($selFunctAna,$selQuanti,$selExplorAna,)=(!$expView)? ('','','') : ($expView eq 'functAna')? (' selected','','') : ($expView eq 'quanti')? ('',' selected','') : ('','',' selected');
	print qq
|	<OPTION value="quanti"$selQuanti>Quantifications</OPTION>
	<OPTION value="explorAna"$selExplorAna>Exploratory Analyses</OPTION>
	<OPTION value="functAna"$selFunctAna>Functional Analyses</OPTION>
	</SELECT></TD>
</TR></TABLE>
|;
}

sub printSampleViewSelection {
	my ($dbh,$sampleView) = @_;
	print qq
|<SCRIPT LANGUAGE="JavaScript">
function updateBioSampleView(newView) {
	top.itemTreeStatus='';
	parent.bioSampleView=newView;
	window.location="$promsPath{cgi}/openProject.cgi?ACT=project&ID=$projectID&branchID="+getSelectedBranchID()+"&VIEW="+newView;
}
</SCRIPT>
<TABLE align="center" bgcolor="$darkColor">
	<TR><TH class="title3 bBorder" nowrap>&nbsp;Biological Samples&nbsp;</TH>
	<TR>
		<TH align="left" nowrap>&nbsp;Sort by:<SELECT name="view" id="view" style="font-weight:bold;" onchange="updateBioSampleView(this.value)">
|;
	my $sthSelProperty = $dbh -> prepare("SELECT P.ID_PROPERTY,NAME
										FROM PROPERTY P
										LEFT OUTER JOIN PROJECT_PROPERTY PP ON P.ID_PROPERTY=PP.ID_PROPERTY
										WHERE P.PROPERTY_TYPE='O' AND (PP.ID_PROJECT=$projectID OR P.IS_VERIFIED=1) ORDER BY NAME"); # linked to current project or verified
	my ($selName,$selSpecies,$selTreatments) = ($sampleView eq 'name')? (' selected','','') : ($sampleView eq 'species')? ('',' selected','') : ($sampleView eq 'treatments')? ('','',' selected') : ('','','');
	print qq
|		<OPTION value="name"$selName>Name</OPTION>
		<OPTION value="species"$selSpecies>Species</OPTION>
		<OPTION value="treatments"$selTreatments>Treatments</OPTION>
		<OPTGROUP label="Properties:">
|;
	$sthSelProperty->execute;
	my $existProp=0;
	while (my($propertyID,$propertyName) = $sthSelProperty -> fetchrow_array) {
		$existProp=1;
		print "<OPTION value=\"property:$propertyID\"";
		print ' selected' if $sampleView eq "property:$propertyID";
		print ">$propertyName</OPTION>\n";
	}
	$sthSelProperty->finish;
	print "<OPTION value=\"\" disabled>None</OPTION>\n" unless $existProp;
	print qq
|	</OPTGROUP>
	</SELECT></TH>
	</TR>
</TABLE>
|;
}

sub createPeptideQuantificationBranches {
	my ($refQuantiData,$sthSA,$quantiID,$refQuantiTree,$type,$refFetched,$refMethodCodes)=@_;
	#$sthQP->execute($quantiID);
	my (%hierarchy,%itemInfo,%distQuantifNames);
	#while (my ($parentQuantiID,$parentQuantiName,$parentFocus,$parentAnnot) = $sthQP->fetchrow_array) { #}
	foreach my $refData (@{$refQuantiData}) {
		my ($parentQuantiID,$parentQuantiName,$parentFocus,$parentAnnot,$parentMehtodID) = @{$refData};
		$parentQuantiName='?#'.$parentQuantiID.'?' unless $parentQuantiName;
		$distQuantifNames{$parentQuantiName}=1;
		my ($sampID,$sampName,$sampPos,$anaID,$anaName,$anaPos);
		if ($refFetched->{$parentQuantiID}) { # info already retrieved: Do not query DB again
			next if $type eq 'peptide'; # keep unused pep quantif
			($sampID,$sampName,$sampPos,$anaID,$anaName,$anaPos)=@{$refFetched->{$parentQuantiID}};
		}
		else {
			$sthSA->execute($parentQuantiID);
			($sampID,$sampName,$sampPos,$anaID,$anaName,$anaPos)=$sthSA->fetchrow_array;
			@{$refFetched->{$parentQuantiID}}=($sampID,$sampName,$sampPos,$anaID,$anaName,$anaPos);
		}
		@{$hierarchy{$sampID}{$anaID}{$parentQuantiID}}=($parentQuantiName,$parentFocus,$parentAnnot,$parentMehtodID);
		@{$itemInfo{SAMPLE}{$sampID}}=($sampPos,$sampName);
		@{$itemInfo{ANALYSIS}{$anaID}}=($anaPos,$anaName);
	}
	my $numDistQuantName=scalar keys %distQuantifNames;
	foreach my $sampID (sort{$itemInfo{SAMPLE}{$a}[0]<=>$itemInfo{SAMPLE}{$b}[0]} keys %hierarchy) {
		my $sampName=$itemInfo{SAMPLE}{$sampID}[1];
		foreach my $anaID (sort{$itemInfo{ANALYSIS}{$a}[0]<=>$itemInfo{ANALYSIS}{$b}[0]} keys %{$hierarchy{$sampID}}) {
			my $anaName=$itemInfo{ANALYSIS}{$anaID}[1];
			foreach my $parentQuantiID (sort{$a<=>$b} keys %{$hierarchy{$sampID}{$anaID}}) {
				my ($parentQuantiName,$parentFocus,$parentAnnot,$parentMehtodID)=@{$hierarchy{$sampID}{$anaID}{$parentQuantiID}};
				my ($parentLabel)=($parentAnnot && $parentAnnot=~/^LABEL=([^:]+)/)? $1 : 'FREE';
				my $parentIsLabeledFlag=($parentLabel eq 'FREE')? '_no_label' : '_label';
				my $parentMethodFlag=($refMethodCodes->{$parentMehtodID} eq 'SWATH')? '_swath' : '';
				my $fullQuantName="$sampName > $anaName > $parentQuantiName";
				my $dispName=($numDistQuantName > 1)? &promsMod::shortenName($fullQuantName,25) : &promsMod::shortenName("$sampName > $anaName > ~",25);
				if ($type eq 'design') { # parent of design quantif
					push @{$refQuantiTree},[3,'quantification',"$parentQuantiID:$quantiID",'',"$parentFocus$parentIsLabeledFlag$parentMethodFlag",1,$dispName,$fullQuantName]; # double id since quanti can be found multiple times
				}
				else {
					push @{$refQuantiTree},[2,'quantification',$parentQuantiID,'',"$parentFocus$parentIsLabeledFlag$parentMethodFlag",1,$dispName,$fullQuantName];
				}
			}
		}
	}

}

####>Revision history<####
# 1.6.9 Display different icons for used and unused BioSample. Requires promsConfig &ge; 2.9.10 (PP 11/04/19)
# 1.6.8 Minor update to handle unexpected undefined quantification name (PP 13/09/18)
# 1.6.7 order cluster and PCA by ana_type (SL 31/01/18)
# 1.6.6 Last open date also updated by bioinfo (PP 11/10/17)
# 1.6.5 Add motif analyses and heatMap motif analyses (SL 31/05/17)
# 1.6.4 Records last open date and user (PP 17/02/17)
# 1.6.3 Icon for SWATH peptide quantification (PP 02/08/16)
# 1.6.2 Change border color (PP 14/03/16)
# 1.6.1 Left border color code to help navigation (PP 19/10/15)
# 1.6.0 Bug fix in quantif nav tree for no-design protein quantifs (PP 13/05/15)
# 1.5.9 Improved bioSample sorting & parent quantif display in nav trees (PP 05/05/15)
# 1.5.8 Order Exp/Func Analyses by name (PP 28/01/15)
# 1.5.7 Minor display changes (PP 07/11/14)
# 1.5.6 add functional analyses item (GO analyses + Pathway analyses) (SL 16/10/14)
# 1.5.5 Add Peptide quantification in order to show multiple-quantifications aligned (GA 01/10/14)
# 1.5.4 change selectedAction button to view summary when comes from properties/treatment (SL 26/09/14)
# 1.5.3 Added Properties & Treatments views (PP 27/08/14)
# 1.5.2 Minor changes in Explor Ana view (PP 17/07/14)
# 1.5.1 display biological samples by project, new function printSampleViewSelection and sort by view (SL 24/06/14)
# 1.5.0 add new branch exploratory analyses (SL 22/04/14)
# 1.4.9 New icons for label(-free) quantification (PP 03/09/13)
# 1.4.8 Added forgotten disconnect if !$experimentView (PP 17/04/13)
# 1.4.7 Minor bug fix for undef selectOption in gel2d frame (PP 10/04/13)
# 1.4.6 Disable display of child GO/quantification analyses (FY 06/12/12)
# 1.4.5 User status is no longer displayed in navFrame (PP 23/11/12)
# 1.4.4 Added -charset=>'UTF-8' to header (PP 24/09/12)
# 1.4.3 Code re-organization & optimization (PP 19/04/12)
# 1.4.2 Merge 1.3.9bPP and 1.4.1GA (03/04/12)
# 1.4.1 Fix the bug with 'setPopup is not defined' for GO/QUANTI frames + Group design for XIC (GA 07/02/2012)
# 1.4.0 Remove the design from the hierarchy (GA 16/12/2011) + Fusion with GO ANALYSIS
# 1.3.9b Minor fix (FY 27/01/12)
# 1.3.9 Add GO analysis items to item frame (fy 05/06/11)
# 1.3.7 New protein list management (ppoullet 11/01/2011)
