#!/usr/local/bin/perl -w

#############################################################################
# selectGOOption.cgi    1.0.6                                               #
# Authors: P. Poullet, G. Arras & F. Yvon (Institut Curie)                  #
# Contact: myproms@curie.fr                                                 #
# Generates list of options available to user, in case of GO analysis item  #
# Displayed in optionFrame                                                  #
# Called by openProject JavaScript tree navigation function for GO analyses #
#############################################################################
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
#my ($anaType, $goAnaID)=split(":",param('ID'));
my ($itemType,$itemID)=split(":",param('ID'));

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

####>Fetching item informations<####
my ($goAnaName,$aspectStrg,$type) = ($itemType eq 'go_analysis')? $dbh->selectrow_array("SELECT NAME,ASPECT,GOA_TYPE FROM GO_ANALYSIS WHERE ID_GOANALYSIS=$itemID") : ('New GO analysis','','');
my ($expName,$projectID) = ($itemType eq 'go_analysis')? $dbh->selectrow_array("SELECT EXPERIMENT.NAME,EXPERIMENT.ID_PROJECT FROM EXPERIMENT,GO_ANALYSIS WHERE EXPERIMENT.ID_EXPERIMENT=GO_ANALYSIS.ID_EXPERIMENT AND GO_ANALYSIS.ID_GOANALYSIS=$itemID") : $dbh->selectrow_array("SELECT NAME, ID_PROJECT from EXPERIMENT where ID_EXPERIMENT=$itemID");
my ($validProt)=$dbh->selectrow_array("SELECT ID_PROTEIN FROM PROTEIN WHERE ID_PROJECT=$projectID LIMIT 0,1"); # 1 is enough
$validProt=1 if $validProt;

####>Fetching user information<####
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
my $projectFullAccess=($projectAccess eq 'bioinfo' ||$projectAccess eq 'mass' || $projectAccess =~ /super/)? 1 : 0;

$dbh->disconnect;

my %aspectName = ( 'P' => 'Biological Process',
                  'C' => 'Cellular Component',
                  'F' => 'Molecular Function'
                  );

####>HTML<####
print header(-charset=>'utf-8');
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
			if ('$itemType' == 'go_analysis') {
				selectedButton=document.getElementById('summary');
			}
			else {
				selectedButton=document.getElementById('goAnalysis');
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
if ($itemType eq 'go_analysis') {
	print qq
|	if (action=='summary' \|\| action=='edit') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/editGOAnalysis.cgi?ACT="+action+"&ID=$itemID";
	} else if (action=='delete'){
		if(confirm('Delete GO analysis $goAnaName ?')){
			top.promsFrame.resultFrame.location="$promsPath{cgi}/editGOAnalysis.cgi?ACT="+action+"&ID=$itemID";
		}
	} else if (action=='C' \|\| action=='F' \|\| action=='P'){
		top.promsFrame.resultFrame.location="$promsPath{cgi}/displayGOAnalysis.cgi?ASPECT="+action+"&ID=$itemID";
	} else if (action == 'heatmap') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/displayGOQuantiAnalysis.cgi?ID=$itemID";
	}
|;
}
else {
	print qq
|	if (action=='goAnalysis') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/startGOAnalysis.cgi?id_exp=$itemID";
	}
	else if (action=='goQuantiAnalysis') {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/selectAnalyses.cgi?ID=experiment:$itemID&id_project=$projectID&callType=goQuantiAnalysis";
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
<FONT class="title2">$expName > GO Analyses > <FONT color=#DD0000>$goAnaName</FONT></FONT><BR>
<TABLE cellpadding=0 cellspacing=0>
<TR align=center valign=top>
|;
if ($itemType eq 'go_analysis') {
	print "<TD nowrap><INPUT type=\"button\" id=\"summary\" style=\"width:80px\" value=\"Summary\" onclick=\"selectOption(this)\">";
	#if($projectFullAccess){
		print "<BR>\n<INPUT type=\"button\" id=\"edit\" style=\"width:80px\" value=\"Edit\" onclick=\"selectOption(this)\">";
		print "</TD><TD nowrap><INPUT type=\"button\" id=\"delete\" style=\"width:80px\" value=\"Delete\" onclick=\"selectOption(this)\">";
	#}
	print "</TD>";
	if(lc($type) eq 'heatmap'){
		print qq
		|<TD nowrap><INPUT type="button" id="heatmap" style="font-weight:bold" value="Heatmap" onclick="selectOption(this);"></TD>
		|;
	}
	else {
		foreach my $aspect (sort { $aspectName{$a} cmp $aspectName{$b} } keys(%aspectName)){
			my $actionString;
			if($aspectStrg =~ /$aspect/){
				$actionString = "onclick=\"selectOption(this)\"";
			} else {
				$actionString = "disabled";
			}
			print "<TD nowrap><INPUT type=\"button\" id=\"$aspect\" style=\"font-weight:bold\" value=\"$aspectName{$aspect}\" $actionString></TD>"
		}
	}
}
else {
	print "<TD nowrap>",&goAnalysisButton,"</TD>";
	print "<TD nowrap>",&goQuantiAnalysisButton,"</TD>";
}
print qq
|</TR></TABLE>
</DIV>
</BODY>
</HTML>
|;


###<GO analysis>####
sub goAnalysisButton {
	return "<INPUT type=\"button\" id=\"goAnalysis\" style=\"width:160px\" value=\"Start GO Analysis\" onclick=\"selectOption(this)\">";
}
###<GO analysis>####
sub goQuantiAnalysisButton {
	return "<INPUT type=\"button\" id=\"goQuantiAnalysis\" style=\"width:160px\" value=\"Start Q. GO Analysis\" onclick=\"selectOption(this)\">";
}

####>Revision history<####
# 1.0.6 Minor modification to make biologists able to remove GO analysis (GA 20/06/17)
# 1.0.5 CSS style for active button (PP 15/03/16)
# 1.0.4 Left border color code to help navigation (PP 19/10/15)
# 1.0.3 UTF-8 (PP 18/12/14)
# 1.0.2 Add start GO and start Q GO , now call from openProject functional analyses (SL )
# 1.0.1 Add specific options for GO analyses applied to quantification data (FY 07/12/12)
# 1.0.0 New script for action buttons in option frame for go analyses (FY 06/07/11)