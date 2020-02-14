#!/usr/local/bin/perl -w

################################################################################
# selectOptionExplorAna.cgi    1.0.6                                          #
# Authors: P. Poullet, S.Liva (Institut Curie)                                 #
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

#print "branchID$branchID selectOptionQuanti.cgi<BR>\n";
##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

#####>Fetching item info<####
my @itemInfo=&promsMod::getItemInfo($dbh,$item,$itemID,$childID);
my $itemType=$itemInfo[-1]{'TYPE'};
my @childItems=&promsMod::getItemChild($item); # '' if $item is ANALYSIS
my @itemParents=&promsMod::getItemParent($item);

#####>Fetching user information<####
my $projectID=$itemInfo[0]{'ID'};
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
my $projectFullAccess=($projectAccess =~ /bioinfo|mass|manag|super/)? 1 : 0;

####>Checking if item has children<####
my $numChildren=0;
my $focus='';
my ($anaType, $expID ,$status) = $dbh -> selectrow_array("SELECT ANA_TYPE, ID_EXPERIMENT, STATUS from EXPLORANALYSIS where ID_EXPLORANALYSIS = $itemID");


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
		else {
			top.promsFrame.resultFrame.location="$promsPath{cgi}/manageExploratoryAnalyses.cgi?ACT="+action+"&ITEM=$item&explorID=$itemID&experimentID=$expID&PROJECT_ID=$projectID";
		}
	}
	else if (action == 'displayPCA' \|\| action == 'displayClustering') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/"+action+".cgi?&explorID=$itemID&experimentID=$expID&PROJECT_ID=$projectID";
	}
	else if (action == 'displayPCAPeptide' \|\| action == 'displayClusteringPeptide') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/"+action+".cgi?&explorID=$itemID&experimentID=$expID&PROJECT_ID=$projectID";
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
print "<TABLE cellpadding=0 cellspacing=0>\n";
print "<TR align=center valign=top>\n";

####>Bioinformatician or Massist<####
if ($projectAccess ne 'guest') { #$projectFullAccess
	if ($item eq 'EXPLORANALYSIS') {
		print "<TD nowrap>",&displayItemButton;
		print "<BR>\n",&editItemButton,"</TD>\n";
		print "<TD>",&delItemButton,"</TD>\n";
	}
}

####>Guest<####
else { # ($projectAccess eq 'guest')
	print "<TD nowrap>",&displayItemButton,"</TD>\n";
}

####>All Users<####
if ($anaType eq ('PCA') || $anaType eq ('PCAPEP')) {print "<TD nowrap>",&displayPCAButton,"</TD>\n";}
else {print "<TD nowrap>",&displayClusterButton,"</TD>\n";} # if $anaType eq 'cluster';

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
####<Delete Item>####
sub delItemButton {
	my $disabledString = ($numChildren)?' disabled':'';
	return "<INPUT type=\"button\" id=\"delete\" style=\"width:80px\" value=\"Delete\" onclick=\"selectOption(this)\"$disabledString>";
}

sub displayPCAButton {
	my $disabledString = ($status != 1)? ' disabled' : '';
	my $strgPcaID=($anaType eq 'PCA')? "displayPCA" : "displayPCAPeptide";
	#return "<INPUT type=\"button\" id=\"displayPCA\" style=\"width:120px\" value=\"Display PCA\" onclick=\"selectOption(this)\"$disabledString>";
	return "<INPUT type=\"button\" id=\"$strgPcaID\" style=\"width:120px\" value=\"Display PCA\" onclick=\"selectOption(this)\"$disabledString>";

}

sub displayClusterButton {
	my $disabledString = ($status != 1)? ' disabled' : '';
	my $strgClusterID=($anaType eq 'cluster')? "displayClustering" : "displayClusteringPeptide";
	#return "<INPUT type=\"button\" id=\"displayClustering\" style=\"width:140px\" value=\"Display Clustering\" onclick=\"selectOption(this)\"$disabledString>";
	return "<INPUT type=\"button\" id=\"$strgClusterID\" style=\"width:140px\" value=\"Display Clustering\" onclick=\"selectOption(this)\"$disabledString>";
}


####>Revision history<####
# 1.0.6 Compatible with Petide exploratory analyses (SL ../11/17)
# 1.0.5 CSS style for active button (PP 15/03/16)
# 1.0.4 Left border color code to help navigation (PP 19/10/15)
# 1.0.3 Allow guest to display Analysis (PP 16/09/14)
# 1.0.2 Clean useless parameters for display(PCA|CLustering) (PP 08/08/14)
# 1.0.1 Minor bug fixes (PP 27/06/14)
# 1.0.0 New script to handle exploratory analyses navigation
