#!/usr/local/bin/perl -w

################################################################################
# selectOptionQuanti.cgi    1.2.7                                              #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
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

#print header(-charset=>'UTF-8'); warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};

###################
####>Arguments<####
###################
my $branchID=param('branchID');
my ($item,$itemID,$childID)=split(':',$branchID);# For XIC extractions linked to a design (TnPQ or Pep-Ratio) quantification
$item=uc($item);
my $fromExplorAna = (param('fromExplorAna'))? param('fromExplorAna') : '';
#print "branchID$branchID selectOptionQuanti.cgi<BR>\n";
##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

####>Fetching item info<####
my @itemInfo=&promsMod::getItemInfo($dbh,$item,$itemID,$childID);
my $itemType=$itemInfo[-1]{'TYPE'};
my @childItems=&promsMod::getItemChild($item);
my @itemParents=&promsMod::getItemParent($item);

####>Fetching user information<####
my $projectID=$itemInfo[0]{'ID'};
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=$userInfo[2]->{$projectID};
my $projectFullAccess=($projectAccess =~ /bioinfo|mass|manag|super/)? 1 : 0;
my $nbExpCond=0;
my $okAddQuantif=0;
my $quantifType='';
my $status;# for quantification only

####>Checking if parent Experiment has quantifs<####
my $experimentID=$itemInfo[1]{'ID'};
my ($hasProtQuantifs)=$dbh->selectrow_array("SELECT 1 FROM SAMPLE S,ANALYSIS A,ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND ID_EXPERIMENT=$experimentID AND FOCUS='protein' AND STATUS > 0 LIMIT 0,1");

####>Checking if item has children<####
my $hasChildren=0;
my $focus='';
if ($item eq 'DESIGN') {
	($nbExpCond)=$dbh->selectrow_array("SELECT COUNT(*) FROM EXPCONDITION WHERE ID_DESIGN=$itemID");
	if ($nbExpCond) {
		$hasChildren=$nbExpCond;
		#($okAddQuantif)=$dbh->selectrow_array("SELECT COUNT(*) FROM EXPCONDITION EC,OBSERVATION O WHERE EC.ID_EXPCONDITION=O.ID_EXPCONDITION AND EC.ID_DESIGN=$itemID");
		($okAddQuantif)=$dbh->selectrow_array("SELECT 1 FROM EXPCONDITION EC,OBS_EXPCONDITION OE WHERE EC.ID_EXPCONDITION=OE.ID_EXPCONDITION AND EC.ID_DESIGN=$itemID LIMIT 0,1");
	}
	else {
		foreach my $childItem (@childItems) {
			next if $childItem eq 'EXPCONDITION';
			($hasChildren)=$dbh->selectrow_array("SELECT 1 FROM $childItem WHERE ID_DESIGN=$itemID LIMIT 0,1");
			last if $hasChildren;
		}
	}
}
elsif ($item eq 'QUANTIFICATION') {
	($quantifType,$status,$focus)=$dbh->selectrow_array("SELECT QUANTIFICATION_METHOD.CODE,STATUS,FOCUS FROM QUANTIFICATION,QUANTIFICATION_METHOD WHERE QUANTIFICATION.ID_QUANTIFICATION_METHOD=QUANTIFICATION_METHOD.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$itemID");
	#$anaString='&ACT=display' if $quantifType eq 'XIC';
	#if ($quantifType eq 'XIC') {
	#	my ($analysisID)=$dbh->selectrow_array("SELECT ID_ANALYSIS FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$itemID");
	#	$anaString="id_ana=$analysisID";
	#}
	if ($quantifType eq 'SWATH') {
        $hasChildren=1;
    }
    else{
		($hasChildren)=$dbh->selectrow_array("SELECT 1 FROM PARENT_QUANTIFICATION WHERE ID_PARENT_QUANTIFICATION=$itemID LIMIT 0,1");
	}
	unless ($hasChildren) {
		($hasChildren) = $dbh->selectrow_array("SELECT 1 FROM GOANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$itemID LIMIT 0,1");
	}
	unless ($hasChildren) {
		($hasChildren) = $dbh->selectrow_array("SELECT 1 FROM EXPLORANA_QUANTIF WHERE ID_QUANTIFICATION=$itemID LIMIT 0,1");
	}
}

$dbh->disconnect;

#######################
####>Starting HTML<####
#######################
my ($popBgColor,$popBdColor)=&promsConfig::getPopupColors;
print header(-charset=>'UTF-8');
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
		if (!selectedButton) { // button is not defined for current project item => switch to 'summary'
			selectedButton=document.getElementById('summary');
		}
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
	if (action != 'monitor') {
		top.promsFrame.selectedAction=action;
	}

	if (action=='summary' \|\| action=='delete' \|\| action=='edit') {
		if (action=='delete' && !confirm('Delete selected item ?')) {
			selectOption(document.getElementById('summary'));
		}
		else if ('$item'=='QUANTIFICATION') {
			if (action=='summary' && '$focus'=='protein') {
				if ('$quantifType'=='SSPA') {
					top.promsFrame.resultFrame.location="$promsPath{cgi}/showSSProtQuantification.cgi?CALL=quanti&id_quantif=$itemID&ACT="+action;
				}
				else {
					top.promsFrame.resultFrame.location="$promsPath{cgi}/showProtQuantification.cgi?CALL=quanti&id_quantif=$itemID&ACT="+action;
				}
			}
			else {
				top.promsFrame.resultFrame.location="$promsPath{cgi}/manageQuantification.cgi?ID=$itemID&PROJECT_ID=$projectID&branchID=$branchID&ACT="+action;
			}
		}
		else { // DESIGN
			top.promsFrame.resultFrame.location="$promsPath{cgi}/manageDesignCondition.cgi?ACT="+action+"&ITEM=$item&ID=$itemID&PARENT=$itemParents[0]&PROJECT_ID=$projectID";
		}
	}
	else if (action=='quantiData') {
		if ('$focus'=='peptide') {
			top.promsFrame.resultFrame.location="$promsPath{cgi}/showPepQuantification.cgi?CALL=quanti&id_quantif=$itemID";
		}
		else { // protein
			if ('$quantifType'=='SSPA') {
				top.promsFrame.resultFrame.location="$promsPath{cgi}/showSSProtQuantification.cgi?CALL=quanti&id_quantif=$itemID&ACT=select";
			}
			else {
				top.promsFrame.resultFrame.location="$promsPath{cgi}/showProtQuantification.cgi?CALL=quanti&id_quantif=$itemID&ACT=select";
			}
		}
	}
	else if (action=='compareQuantifs') {
		if ('$item'=='QUANTIFICATION') {
			top.promsFrame.resultFrame.location="$promsPath{cgi}/compareQuantifications.cgi?&id_project=$projectID&parentItem=$itemInfo[-2]{TYPE}:$itemInfo[-2]{ID}"; // Design
		}
		else { // DESIGN
			top.promsFrame.resultFrame.location="$promsPath{cgi}/compareQuantifications.cgi?&id_project=$projectID&parentItem=$item:$itemID";
		}
	}
	else if (action=='exportQuantifs') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/startExploratoryAnalysis.cgi?ID=$experimentID&ACT=export";
	}
	else if (action == 'goAnalysis') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/startGOQuantiAnalysis.cgi?quanti=$itemID";
	}
	//else if ('$item'=='QUANTIFICATION' && action != 'launch'){
	//	top.promsFrame.resultFrame.location="$promsPath{cgi}/manageQuantification.cgi?ACT="+action+"&ID=$itemID&branchID=$branchID&PROJECT_ID=$projectID";
	//}
	//else if (action=='add') {
	//	top.promsFrame.resultFrame.location="$promsPath{cgi}/manageDesignCondition.cgi?ACT=add&ITEM=$childItems[0]&ID=$itemID&PARENT=$item&PROJECT_ID=$projectID";
	//}
	else if (action=='addDesign') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/manageDesignCondition.cgi?ACT=add&ITEM=DESIGN&ID=$itemID&PARENT=$item&PROJECT_ID=$projectID";
	}
	else if (action=='addQuantification') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/selAna4Quantification.cgi?ID=$branchID&quantifType=DESIGN";
	}
	else {alert('Unrecognized option');}
}
function autoSelectButton(buttonId) { // an action was cancelled or selected from resultFrame
	//currentButton.style.color='#000000';
	currentButton.className='';
	currentButton=document.getElementById(buttonId);
	//currentButton.style.color='#DD0000';
	currentButton.className='selectedButton';
	top.promsFrame.selectedAction=buttonId;
}
function watchQuantifications() {
	var watchQuantifWin=window.open("$promsPath{cgi}/watchQuantifications.cgi",'WatchQuantifWindow','width=1000,height=500,scrollbars=yes,resizable=yes');
	watchQuantifWin.focus();
	parent.optionFrame.selectOption(parent.optionFrame.document.getElementById('summary')); // refresh optionFrame with summary option
}
var currentButton; // currently selected button
</SCRIPT>
</HEAD>
<BODY style="margin:0px;" onload="selectOption()">
<DIV class="navigation navBorder2">
|;

#################################
####>Printing Item hierarchy<####
#################################
my @titleString;
my $titleChar='';
foreach my $i (1..$#itemInfo) { # 1 not 0 (skip project name)
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
print "<TABLE cellpadding=0 cellspacing=0>\n";
print "<TR align=center valign=top>\n";

####>Bioinformatician or Massist<####
if ($projectAccess ne 'guest') { #$projectFullAccess
	if ($item eq 'DESIGN') {
		print "<TD nowrap>",&displayItemButton;
		print "<BR>\n",&editItemButton,"</TD>\n";
		print "<TD>",&delItemButton,"</TD>\n";
		#my $isSelectable=($nbExpCond <= 1)? ' disabled':'';
		print "<TD nowrap>"; #,&addQuantificationButton($isSelectable);
		print &addChildButton('QUANTIFICATION',0),"<BR>\n" if $nbExpCond >=2;
		print "</TD>\n";
		print "<TD nowrap>",&monitorQuantificationButton,"<BR>\n",&compQuantifsButton,"</TD>\n";
		print "<TD nowrap>",&exportQuantifsButton,"</TD>\n" if $hasProtQuantifs;
	}
	elsif ($item eq 'QUANTIFICATION') {
		#my $dbh=&promsConfig::dbConnect;
		#my $disabled=($status!=-1)? ' disabled': '';
		print "<TD nowrap>",&displayItemButton;
		#print "<BR>\n",&editItemButton if $status < 0;
		print "<BR>\n",&editItemButton;# if $status < 0;
		print "</TD>\n";
		print "<TD nowrap>",&delItemButton,"</TD>\n";
		print "<TD nowrap>",&quantiDataButton,"</TD>\n" if $status > 0;
		print "<TD nowrap>",&monitorQuantificationButton,"<BR>\n",&compQuantifsButton,"</TD>\n";
		print "<TD nowrap>",&exportQuantifsButton,"</TD>\n" if $hasProtQuantifs;
		print "<TD nowrap>",&goAnalysisButton,"</TD>\n" if ($quantifType=~/PROT_RATIO_PEP|TNPQ/ && $status > 0);
	}
}

####>Guest<####
else { # ($projectAccess eq 'guest')
	print "<TD nowrap>",&displayItemButton,"</TD>\n";
	print "<TD nowrap>",&quantiDataButton,"</TD>\n" if ($item eq 'QUANTIFICATION' && $status > 0);
	print "<TD nowrap>",&exportQuantifsButton,"</TD>\n" if $hasProtQuantifs;
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
	return "<INPUT type=\"button\" id=\"edit\" style=\"width:80px\" value=\"Edit\" onclick=\"selectOption(this)\">";
}
####<Add Child Item>####
sub addChildButton { # only for adding quantification so far
	#my ($addedChild,$buttonID)=($_[0])? ($_[0],'addGel2D') : ($childItems[0],'add');
	my ($addedChild,$multi)=@_;
	my $childType=&promsMod::getItemType($addedChild);
	$childType.='(s)' if $multi;
	my ($buttonID,$disabStrg);
	if ($addedChild eq 'QUANTIFICATION') {
		$buttonID='addQuantification';
		$disabStrg=($okAddQuantif && (!$fromExplorAna))? '' : 'disabled';
	}
	else {
		$buttonID='add';
		$disabStrg= '';
	}
	return "<INPUT type=\"button\" id=\"$buttonID\" style=\"width:160px\" value=\"Add $childType\" onclick=\"selectOption(this)\" $disabStrg>";
}
####<Delete Item>####
sub delItemButton {
	my $disabledString = ($hasChildren)?' disabled':'';
	return "<INPUT type=\"button\" id=\"delete\" style=\"width:80px\" value=\"Delete\" onclick=\"selectOption(this)\"$disabledString>";
}


####<Process Quantification>####
#sub launchQuantificationButton {
#	my ($disabled)=@_;
#	return "<INPUT type=\"button\" id=\"launch\" style=\"width:180px\" value=\"Launch Quantification\" onclick=\"selectOption(this)\" $disabled>";
#}
####<Monitor Quantification>####
sub monitorQuantificationButton {
	my $disabledString = ($fromExplorAna)? ' disabled' : '';
	return "<INPUT type=\"button\" id=\"monitor\" style=\"width:180px\" value=\"Monitor Quantification(s)\" onclick=\"watchQuantifications()\"$disabledString>";
}
####<Call listAnaQuantification>####
sub quantiDataButton {
	#my ($disabled)=@_;
	return "<INPUT type=\"button\" id=\"quantiData\" style=\"width:180px\" value=\"Display Results\" onclick=\"selectOption(this)\">";
}
####<Compare Quantification>####
sub compQuantifsButton {
	return "<INPUT type=\"button\" id=\"compareQuantifs\" style=\"width:180px\" value=\"Compare Quantifications\" onclick=\"selectOption(this)\">";
}
####<Export Quantifications>####
sub exportQuantifsButton {
	return "<INPUT type=\"button\" id=\"exportQuantifs\" style=\"width:175px\" value=\"Export Quantifications\" onclick=\"selectOption(this)\">";
}
####<GO analysis>####
sub goAnalysisButton{
	return "<INPUT type=\"button\" id=\"goAnalysis\" style=\"width:180px\" value=\"Start Q. GO Analysis\" onclick=\"selectOption(this)\">";
}
####<Add Analysis to EXPCONDITION>####
#sub addAnalysisButton {
#	return "<INPUT type=\"button\" id=\"addAnalysis\" style=\"width:140px\" value=\"Manage content\" onclick=\"selectOption(this)\">";
#}


####>Revision history<####
# 1.2.7 Added Export Quantifications options (PP 07/11/17)
# 1.2.6 Calls showSSProtQuantification.cgi for SSPA quantification summary (PP 17/08/16)
# 1.2.5 Calls showProtQuantification.cgi for design-based protein quantification summary (PP 25/07/16)
# 1.2.4 Disable 'delete' button if param quantifType eq 'SWATH' (MLP 14/06/16)
# 1.2.3 CSS style for active button (PP 15/03/16)
# 1.2.2 Left border color code to help navigation (PP 19/10/15)
# 1.2.1 Add Compare Quantification button to Design and Protein Quantification levels (PP 22/09/15)
# 1.2.0 Add Exploratory Analysis in list of children for quantif deletability (PP 28/01/15)
# 1.1.9 New Observation management: OBS_EXPCONDITION table (PP 15/07/14)
# 1.1.8 Add scrollbar option for WatchQuantifWindow (GA 27/06/14)
# 1.1.7 disable action on quantifications if design come from exploratory analyses (SL 28/04/14)
# 1.1.6 Disable 'delete' button if item has children, including GO analyses (FY 10/10/13)
# 1.1.5 charset=> utf-8 (PP 09/09/13)
# 1.1.4 Changed value of param quantifType from TNPQ to DESIGN for selAna4Quantification.cgi call<BR>and action button reorganization (PP 29/08/13)
# 1.1.3 Add edition of name for label-free quantifications (GA 13/08/13)
# 1.1.2 Updated for OBSERVATION table (PP 02/07/13)
# 1.1.1 Disable Add Quantif button when no analysis in exp_conditions (PP 03/06/13)
# 1.1.0 Adding GO analysis starting button (FY 02/04/13)
# 1.0.9 Minor modif to take into account error handling, i.e. $status=-2 (GA 07/03/13)
# 1.0.8 Code cleaning (PP 04/01/13)
# 1.0.7 Modification of quantification name in HTML to print the DESIGN hierarchy (GA 07/11/12)
# 1.0.6 Fixed delete check on quantification (PP 26/04/12)
# 1.0.5 Modifications for optimal call of listAnaQuantification (PP 19/04/12)
# 1.0.4 Clean version of the script (GA 11/04/12)
# 1.0.3 Make the call to fillCondition.cgi disappear (GA 03/04/12)
# 1.0.2 Update to current quantification methods: TNPQ and Protein Ratio (GA 07/03/12)
# 1.0.1 Update to current quantification methods: only TNPQ (GA 06/01/12)
# 1.0.0 New sccript to handle quantification navigation