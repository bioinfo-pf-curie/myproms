#!/usr/local/bin/perl -w

#################################################################################
# selectOptionGSEA.cgi    1.0.0                                                 #
# Authors: V. Laigle (Institut Curie)                                           #
# Contact: myproms@curie.fr                                                     #
# Generates list of options available to user for GSEA item                     #
# Displayed in optionFrame                                                      #
# Called by openProject JavaScript tree navigation for GSEA analyses            #
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
use warnings;

#####################
### Configuration ###
#####################
my %promsPath = &promsConfig::getServerInfo;
my $userID = $ENV{'REMOTE_USER'};
my $dbh=&promsConfig::dbConnect;

my ($itemType, $itemID) = split(":", param('ID'));

# Fetching item informations
my ($gseaName, $status) = ($itemType eq 'gseanalysis')? $dbh->selectrow_array("SELECT NAME, STATUS FROM PATHWAY_ANALYSIS WHERE ID_PATHWAY_ANALYSIS = $itemID") : ('New Gene Set Enrichment Analysis', '');
my ($expName, $projectID) = ($itemType eq 'gseanalysis')? $dbh->selectrow_array("SELECT E.NAME, E.ID_PROJECT FROM EXPERIMENT E, PATHWAY_ANALYSIS PA WHERE E.ID_EXPERIMENT = PA.ID_EXPERIMENT AND PA.ID_PATHWAY_ANALYSIS = $itemID") : $dbh->selectrow_array("SELECT NAME, ID_PROJECT FROM EXPERIMENT WHERE ID_EXPERIMENT = $itemID");

# Fetching user information
my @userInfo = &promsMod::getUserInfo($dbh, $userID, $projectID);
my $projectAccess = ${$userInfo[2]}{$projectID};
my $projectGuestAccess = ($projectAccess eq 'guest')? 1 : 0;
$dbh->disconnect;

############
### HTML ###
############
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
        selectedButton = document.getElementById(top.promsFrame.selectedAction); // default value stored in promsFrame
        if (!selectedButton) { // button is not defined for current project item => switch to 'summary'
            if ('$itemType' == 'gseanalysis') {
                selectedButton = document.getElementById('summary');
            } else {
                selectedButton = document.getElementById('startGSEA');
            }
        }
    }
    if (currentButton) {
        currentButton.className = '';
    }
    currentButton = selectedButton;
    currentButton.className = 'selectedButton';
    var action = currentButton.id;
|;

if ($itemType eq 'gseanalysis') {
    print qq
|        if (action == 'summary' \|\| action == 'edit') {
            top.promsFrame.resultFrame.location = "$promsPath{cgi}/manageGSEA.cgi?ACT="+action+"&ID=$itemID";
        } else if (action == 'delete') {
            if(confirm('Delete Gene Set Enrichment Analysis $gseaName ?')){
                top.promsFrame.resultFrame.location = "$promsPath{cgi}/manageGSEA.cgi?ACT="+action+"&ID=$itemID";
            }
        } else if (action == 'displayGSEA') {
            top.promsFrame.resultFrame.location = "$promsPath{cgi}/displayGSEA.cgi?ID=$itemID";
        } else if (action == 'monitorJobs') {
            var monitorWindow = window.open("$promsPath{cgi}/monitorJobs.cgi?filterType=Functional Analysis [GSEA]&filterDateNumber=1&filterDateType=DAY&filterStatus=Queued&filterStatus=Running&filterStatus=Error&filterProject=$projectID",'monitorJobsWindow','width=1200,height=500,scrollbars=yes,resizable=yes');
        }
|;

} else {
    print qq
|        if (action == 'startGSEA') {
            top.promsFrame.resultFrame.location = "$promsPath{cgi}/chooseQuantiGSEA.cgi?action=choose&exp=$itemID";
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
<FONT class="title2">$expName > Gene Set Enrichment Analyses > <FONT color=#DD0000>$gseaName</FONT></FONT><BR>
<TABLE cellpadding=0 cellspacing=0>
<TR align=center valign=top>
|;

if ($itemType eq 'gseanalysis') {
    print "<TD nowrap>", &summaryOption;
    unless ($projectGuestAccess) {
        print "<BR>", &editOption, "</TD>\n";
        print "<TD nowrap>", &deleteOption;
    }
    print "</TD>\n";
    print "<TD nowrap>", &displayGSEAButton, "</TD>\n";
    print "<TD nowrap>", &monitorJobsButton, "</TD>\n";
}
else {
    print "<TD nowrap>",&startGSEAButton,"</TD>\n";
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

sub displayGSEAButton {
    my $disabled = ($status == 1)? " " : " disabled";
    return "<INPUT type=\"button\" id=\"displayGSEA\" value=\"Display GSEA\" onclick=\"selectOption(this)\"$disabled>";
}

sub monitorJobsButton {
    return "<INPUT type=\"button\" id=\"monitorJobs\" value=\"Monitor Jobs\" onclick=\"selectOption(this)\">";
}

sub startGSEAButton {
    return "<INPUT type=\"button\" id=\"startGSEA\" value=\"Start GSEA\" onclick=\"selectOption(this)\">";
}

####>Revision history<####
# 1.0.0 Created, forked from selectOptionPathway.cgi (VL 05/11/20)
