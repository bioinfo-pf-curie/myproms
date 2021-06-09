#!/usr/local/bin/perl -w
#################################################################################
# selectOptionMotif.cgi    1.0.2                                                #
# Authors: P. Poullet, S. Liva (Institut Curie)                   				#
# Contact: myproms@curie.fr                                                     #
# Generates list of options available to user, in case of GO analysis item      #
# Displayed in optionFrame                                                      #
# Called by openProject JavaScript tree navigation for Motif/HeatMap analyses   #
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
my ($anaType) = $dbh->selectrow_array("SELECT ANA_TYPE FROM EXPLORANALYSIS WHERE ID_EXPLORANALYSIS=$itemID");
my ($isParent)=$dbh->selectrow_array("SELECT 1 FROM PARENT_EXPLORANALYSIS WHERE ID_PARENT_EXPLORANALYSIS=$itemID");
my ($countAna)=$dbh->selectrow_array("SELECT COUNT(ID_EXPLORANALYSIS) FROM EXPLORANALYSIS WHERE ANA_TYPE ='MOTIF'");

my ($motifName, $status) = ($itemType eq 'motifanalysis')? $dbh->selectrow_array("SELECT NAME, STATUS FROM EXPLORANALYSIS WHERE ID_EXPLORANALYSIS=$itemID and ANA_TYPE='$anaType'") : ('New Motif analysis', '');
my ($expName,$projectID) = ($itemType eq 'motifanalysis')? $dbh->selectrow_array("SELECT EXPERIMENT.NAME,EXPERIMENT.ID_PROJECT FROM EXPERIMENT,EXPLORANALYSIS WHERE EXPERIMENT.ID_EXPERIMENT=EXPLORANALYSIS.ID_EXPERIMENT AND EXPLORANALYSIS.ID_EXPLORANALYSIS=$itemID AND EXPLORANALYSIS.ANA_TYPE='$anaType'") : $dbh->selectrow_array("SELECT NAME, ID_PROJECT from EXPERIMENT where ID_EXPERIMENT=$itemID");
my ($strgMotif, $strgTitle) = ($itemType eq 'motifanalysis' && $anaType eq "MOTIF")? ("Motif Analysis", "MotifAnalysis") : ("Heatmap Motif Analysis", "HeatmapMotifAnalysis");

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
			if ('$itemType' == 'motifanalysis') {
				selectedButton=document.getElementById('summary');
			}
			else {
				selectedButton=document.getElementById('startMotifAnalysis');
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
if ($itemType eq 'motifanalysis') {
	print qq
|
		if (action=='summary' \|\| action=='edit') {
			top.promsFrame.resultFrame.location="$promsPath{cgi}/manageMotifAnalysis.cgi?ACT="+action+"&ID=$itemID&TYPE=$anaType";
		} else if (action=='delete'){

			if(confirm('Delete motif analysis $motifName?')){
				top.promsFrame.resultFrame.location="$promsPath{cgi}/manageMotifAnalysis.cgi?ACT="+action+"&ID=$itemID&TYPE=$anaType";
			}
		} else if (action=='display$strgTitle'){
		console.log(action);
			top.promsFrame.resultFrame.location="$promsPath{cgi}/display$strgTitle.cgi?ID=$itemID&PROJECT_ID=$projectID&TYPE=$anaType";
		}
|;
}
else {
	print qq
|		if (action=='startMotifAnalysis') {
			top.promsFrame.resultFrame.location="$promsPath{cgi}/startMotifAnalysis.cgi?id_exp=$itemID";
		}else if (action == 'heatMapMotifAnalysis') {
			top.promsFrame.resultFrame.location="$promsPath{cgi}/heatMapMotifAnalysis.cgi?id_exp=$itemID";
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
<FONT class="title2">$expName > Motif Analyses > <FONT color=#DD0000>$motifName</FONT></FONT><BR>
<TABLE cellpadding=0 cellspacing=0>
<TR align=center valign=top>
|;
if ($itemType eq 'motifanalysis') {
	print "<TD nowrap>",&summaryOption;
	unless ($projectGuestAccess) {
		print "<BR>",&editOption,"</TD>\n";
		print "<TD nowrap>",&deleteOption;
	}
	print "</TD>\n";
	print "<TD nowrap>",&displayAnalysisButton,"</TD>\n";
}
else {
	print "<TD nowrap>",&motifAnalysisButton,"</TD>\n";
	print "<TD nowrap>",&heatMapMotifAnalysisButton,"</TD>\n";
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
	my $disab = (!$isParent)? " " : " disabled";
	return "<INPUT type=\"button\" style=\"width:80px\" id=\"delete\" value=\"Delete\" onclick=\"selectOption(this)\"$disab>";
}

###<display motif analysis>####
sub displayAnalysisButton {
	my $disabled = ($anaType eq "MOTIF" && $status == 1 || $anaType eq "HM_MOTIF")? " " : " disabled";
	my $strgID="display$strgTitle";
	return "<INPUT type=\"button\" id=\"$strgID\" value=\"Display $strgMotif\" onclick=\"selectOption(this)\"$disabled>";
}
###<start motif analysis>####
sub motifAnalysisButton {
	return "<INPUT type=\"button\" id=\"startMotifAnalysis\" value=\"Start Motif Analysis\" onclick=\"selectOption(this)\">";
}
###<heatMap motif analysis>###
sub heatMapMotifAnalysisButton {
	my $disabHM = ($countAna<2 || $projectGuestAccess)? " disabled" : " ";
	return "<INPUT type=\"button\" id=\"heatMapMotifAnalysis\" value=\"HeatMap Motif Analysis\" onclick=\"selectOption(this)\"$disabHM>";
}

####>Revision history<####
# 1.0.2 [BUGFIX] Disable edition and deletion by guest users (VL 28/01/21)
# 1.0.1 add revision history tag (SL 14/03/18)
# 1.0.0 New script for action buttons in option frame for motif analyses (SL 31/05/17)