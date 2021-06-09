#!/usr/local/bin/perl -w

#################################################################################
# selectOptionPathway.cgi    1.0.5                                              #
# Authors: P. Poullet, G. Arras & S. Liva (Institut Curie)                      #
# Contact: myproms@curie.fr                                                     #
# Generates list of options available to user, in case of GO analysis item      #
# Displayed in optionFrame                                                      #
# Called by openProject JavaScript tree navigation for Pathway analyses         #
#################################################################################
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
my $userID=$ENV{'REMOTE_USER'};

###################
####>Arguments<####
###################
my ($itemType,$itemID)=split(":",param('ID'));

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

####>Fetching item informations<####
my ($pathwayName, $status) = ($itemType eq 'pathwayanalysis')? $dbh->selectrow_array("SELECT NAME, STATUS FROM PATHWAY_ANALYSIS WHERE ID_PATHWAY_ANALYSIS=$itemID") : ('New Pathway analysis', '');
my ($expName,$projectID) = ($itemType eq 'pathwayanalysis')? $dbh->selectrow_array("SELECT EXPERIMENT.NAME,EXPERIMENT.ID_PROJECT FROM EXPERIMENT,PATHWAY_ANALYSIS WHERE EXPERIMENT.ID_EXPERIMENT=PATHWAY_ANALYSIS.ID_EXPERIMENT AND PATHWAY_ANALYSIS.ID_PATHWAY_ANALYSIS=$itemID") : $dbh->selectrow_array("SELECT NAME, ID_PROJECT from EXPERIMENT where ID_EXPERIMENT=$itemID");

####>Fetching user information<####
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
my $projectGuestAccess = ($projectAccess eq 'guest')? 1 : 0;
$dbh->disconnect;

####>HTML<####
print header(-'content-encoding'=>'no',-charset=>'UTF-8');
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<TITLE>Option Menu</TITLE>
<SCRIPT LANGUAGE="JavaScript">
function selectOption(selectedButton) {
	if (!selectedButton) { // 1st time this page is loaded
		selectedButton=document.getElementById(top.promsFrame.selectedAction); // default value is stored in promsFrame;
		if (!selectedButton) { // button is not defined for current project item => switch to 'summary'
			if ('$itemType' == 'pathwayanalysis') {
				selectedButton=document.getElementById('summary');
			}
			else {
				selectedButton=document.getElementById('startPathAnalysis');
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
|;
if ($itemType eq 'pathwayanalysis') {
	print qq
|		if (action=='summary' \|\| action=='edit') {
			top.promsFrame.resultFrame.location="$promsPath{cgi}/managePathwayAnalysis.cgi?ACT="+action+"&ID=$itemID";
		} else if (action=='delete'){
			if(confirm('Delete pathway analysis $pathwayName ?')){
				top.promsFrame.resultFrame.location="$promsPath{cgi}/managePathwayAnalysis.cgi?ACT="+action+"&ID=$itemID";
			}
		} else if (action=='displayPathAnalysis'){
			top.promsFrame.resultFrame.location="$promsPath{cgi}/displayPathwayAnalysis.cgi?ID=$itemID";
		}
|;
}
else {
	print qq
|		if (action=='startPathAnalysis') {
			top.promsFrame.resultFrame.location="$promsPath{cgi}/startPathwayAnalysis.cgi?id_exp=$itemID";
		}
|;
	}
print qq
|}
var currentButton; // currently selected button
</SCRIPT>
</HEAD>
<BODY style="margin:0px;" onload="selectOption()">
<DIV class="navigation navBorder2">
<FONT class="title2">$expName > Pathway Analyses > <FONT color=#DD0000>$pathwayName</FONT></FONT><BR>
<TABLE cellpadding=0 cellspacing=0>
<TR align=center valign=top>
|;
if ($itemType eq 'pathwayanalysis') {
	print "<TD nowrap>",&summaryOption;
	unless ($projectGuestAccess) {
		print "<BR>",&editOption,"</TD>\n";
		print "<TD nowrap>",&deleteOption;
	}
	print "</TD>\n";
	print "<TD nowrap>",&displayPathwayAnalysisButton,"</TD>\n";
}
else {
	print "<TD nowrap>",&pathwayAnalysisButton,"</TD>\n";
}
print qq
|</TR></TABLE>
</DIV>
</BODY>
</HTML>
|;

sub summaryOption {
	return "<INPUT type=\"button\" style=\"width:80px\" id=\"summary\" value=\"Summary\" onclick=\"selectOption(this)\">";
}
sub editOption {
	return "<INPUT type=\"button\" style=\"width:80px\" id=\"edit\" value=\"Edit\" onclick=\"selectOption(this)\">";
}
sub deleteOption {
	return "<INPUT type=\"button\" style=\"width:80px\" id=\"delete\" value=\"Delete\" onclick=\"selectOption(this)\">";
}

###<display pathway analysis>####
sub displayPathwayAnalysisButton {
	my $disabled = ($status == 1)? " " : " disabled";
	return "<INPUT type=\"button\" id=\"displayPathAnalysis\" value=\"Display Pathway Analysis\" onclick=\"selectOption(this)\"$disabled>";
}
###<start pathway analysis>####
sub pathwayAnalysisButton {
	return "<INPUT type=\"button\" id=\"startPathAnalysis\" value=\"Start Pathway Analysis\" onclick=\"selectOption(this)\">";
}

####>Revision history<####
# 1.0.5 [BUGFIX] Disable edition and deletion by guest users (VL 28/01/21)
# 1.0.4 CSS style for active button (PP 15/03/16)
# 1.0.3 Left border color code to help navigation (PP 19/10/15)
# 1.0.2 change runAndDisplayPathwayAnalysis to displayPathwayAnalysis.cgi, add disabled displayPathwayAnalysis (SL 03/09/15)
# 1.0.1 header with UTF-8 (PP 13/01/14)
# 1.0.0 New script for action buttons in option frame for pathway analyses (SL 16/10/14)
