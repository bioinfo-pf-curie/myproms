#!/usr/local/bin/perl -w

################################################################################
# selectOptionBioSample.cgi    1.0.6                                           #
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
my ($item,$itemID)=split(':',$branchID);# For XIC extractions linked to a design (TnPQ or Pep-Ratio) quantification
$item=uc($item);
my $projectID;
my $hasObs = 0;
my $titleStrg='';

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

if ($item eq 'ALL') {
	$projectID=$itemID;
	$itemID=0;
	$titleStrg='All Biological Samples';
}
elsif ($item eq 'BIOSAMPLE') {
	$projectID = param('projectID');
	($titleStrg)=$dbh->selectrow_array("SELECT NAME FROM BIOSAMPLE WHERE ID_BIOSAMPLE=$itemID");
	#my $sthSelObs = $dbh -> prepare("SELECT ID_OBSERVATION from OBSERVATION where ID_BIOSAMPLE = ?");
	#$sthSelObs -> execute($itemID);
	#while (my($obsID) = $sthSelObs -> fetchrow_array) {
	#	$hasObs = 1;#disable biosample if has observations
	#	last;
	#}
	#$sthSelObs -> finish;
	($hasObs)=$dbh->selectrow_array("SELECT COUNT(O.ID_OBSERVATION) FROM OBSERVATION O
										INNER JOIN ANALYSIS A ON O.ID_ANALYSIS=A.ID_ANALYSIS
										INNER JOIN SAMPLE S ON A.ID_SAMPLE=S.ID_SAMPLE
										INNER JOIN EXPERIMENT E ON S.ID_EXPERIMENT=E.ID_EXPERIMENT
										WHERE O.ID_BIOSAMPLE=$itemID AND E.ID_PROJECT=$projectID"); # BioSample is used in an observation of the current project!
}

######>Fetching user information<####
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
my $projectFullAccess=($projectAccess =~ /guest/)? 0 : 1;
$dbh -> disconnect;

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
			if ('$item' == 'ALL') {
				selectedButton=document.getElementById('properties');
			}
			else {
				selectedButton=document.getElementById('summary');
			}
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
	top.promsFrame.selectedAction=action;

	if (action == 'summary' \|\| action == 'edit' \|\| action == 'delete') {
		parent.resultFrame.location="$promsPath{cgi}/manageBioSample.cgi?ACT="+action+"&projectID=$projectID&biosampleID=$itemID";
	}
	else if (action == 'createBS') { //create a new biological sample
		parent.resultFrame.location="$promsPath{cgi}/manageBioSample.cgi?ACT=add&projectID=$projectID";
	}
	else if (action == 'importBS') { //import biological samples from file
		parent.resultFrame.location="$promsPath{cgi}/importBioSampleData.cgi?projectID=$projectID";
	}
	else if (action == 'createSP') { //create a new species
		parent.resultFrame.location="$promsPath{cgi}/manageSpecies.cgi?ACT=add";
	}
	else if (action == 'properties') {
		parent.resultFrame.location="$promsPath{cgi}/manageProperties.cgi?ACT=summary&projectID=$projectID&type=O";
	}
	else if (action == 'treatments') {
		parent.resultFrame.location="$promsPath{cgi}/manageProperties.cgi?ACT=summary&projectID=$projectID&type=T";
	}
	else if (action == 'link') {
		parent.resultFrame.location="$promsPath{cgi}/linkBioSample2Observations.cgi?biosampleID=$itemID&projectID=$projectID";
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

##########################
#####>List of options<####
##########################
print qq
|<FONT class="title2" color="#DD0000">$titleStrg</FONT>
<BR><TABLE cellpadding=0 cellspacing=0>
<TR align=center valign=top>
|;
if ($item eq 'ALL') {
	print "<TD nowrap>",&createBioSample,"<BR>",&importBioSample,"</TD>\n";
	print "<TD nowrap>",&sampleProperties,"<BR>",&sampleTreatments,"</TD>\n";
	print "<TD nowrap>",&createSpecies,"</TD>\n";
}
elsif ($item eq 'BIOSAMPLE') {
	print "<TD nowrap>",&summaryBiosample;
	print "<BR>",&editBiosample,"</TD>\n";
	print "<TD nowrap>",&delBiosample,"</TD>\n";
	print "<TD nowrap>",&sampleProperties,"<BR>",&sampleTreatments,"</TD>\n";
}
print "<TD nowrap>",&linkToObservations,"</TD>\n";
print qq
|</TR></TABLE>
</DIV>
</BODY>
</HTML>
|;

##############################<<< Buttons Subroutines >>>############################
sub summaryBiosample {
	return "<INPUT type=\"button\" id=\"summary\" style=\"width:80px\" value=\"Summary\" onclick=\"selectOption(this)\">";
}
sub editBiosample {
	my $disabEdit = ($projectFullAccess)? '' : ' disabled';
	return "<INPUT type=\"button\" id=\"edit\" style=\"width:80px\" value=\"Edit\" onclick=\"selectOption(this)\"$disabEdit>";
}
sub delBiosample {
	my $disabDelete = ($hasObs || !$projectFullAccess)? ' disabled' : '';
	return "<INPUT type=\"button\" id=\"delete\" style=\"width:80px\" value=\"Delete\" onclick=\"selectOption(this)\"$disabDelete>";
}
sub createBioSample {
	my $disabSamp = ($projectFullAccess)? '' : ' disabled' ;
	return "<INPUT type=\"button\" id=\"createBS\" style=\"width:170px\" value=\"New Biological Sample\" onclick=\"selectOption(this)\"$disabSamp>";
}
sub importBioSample {
	my $disabSamp = ($projectFullAccess)? '' : ' disabled' ;
	return "<INPUT type=\"button\" id=\"importBS\" style=\"width:170px\" value=\"Import from File\" onclick=\"selectOption(this)\"$disabSamp>";
}
sub createSpecies {
	my $disabSpecies = ($projectFullAccess)? '' : ' disabled' ;
	return "<INPUT type=\"button\" id=\"createSP\" style=\"width:110px\" value=\"New Species\" onclick=\"selectOption(this)\"$disabSpecies>";
}
sub sampleProperties {
	return "<INPUT type=\"button\" id=\"properties\" style=\"width:150px\" value=\"Manage Properties\" onclick=\"selectOption(this)\">";
}
sub sampleTreatments {
	return "<INPUT type=\"button\" id=\"treatments\" style=\"width:150px\" value=\"Manage Treatments\" onclick=\"selectOption(this)\">";
}
sub linkToObservations {
	my $disabLinkObs = ($projectFullAccess)? '' : ' disabled' ;
	return "<INPUT type=\"button\" id=\"link\" style=\"width:150px\" value=\"Link to Observations\" onclick=\"selectOption(this)\"$disabLinkObs>";
}

####>Revision history<####
# 1.0.6 [FEATURE] Added option button for biological sample import from file (PP 10/03/20)
# 1.0.5 CSS style for active button (PP 15/03/16)
# 1.0.4 Left border color code to help navigation (PP 19/10/15)
# 1.0.3 Disable Observation Link option for guest (PP 16/09/14)
# 1.0.2 Added multi-sample Observation link option (PP 27/08/17)
# 1.0.1 Added item title and other minor changes (PP 16/07/14)
# 1.0.0 New script to handle biological sample navigation (SL 02/07/14)
