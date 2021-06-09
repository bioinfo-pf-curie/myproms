#!/usr/local/bin/perl -w

#################################################################################
# selectOptionHydro.cgi    1.0.0                                                #
# Authors: V. Laigle (Institut Curie)                                           #
# Contact: myproms@curie.fr                                                     #
# Generates list of options available to user for Hydrophobicity computations   #
# Displayed in optionFrame                                                      #
# Called from editProjectItem experiment, sample or analysis sections           #
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

my $quantifID = param('ID');
my ($fromItem, $fromID) = split(":", param('FROM'));

# Fetching item informations
my ($quantifName, $status) = $dbh->selectrow_array("SELECT NAME, STATUS FROM QUANTIFICATION WHERE ID_QUANTIFICATION = $quantifID");
my ($projectID, $expName, $sampName, $anaName);
my $itemStrg;
if ($fromItem eq 'analysis') {
    ($projectID, $expName, $sampName, $anaName) = $dbh->selectrow_array("SELECT E.ID_PROJECT, E.NAME, S.NAME, A.NAME FROM ANALYSIS A INNER JOIN SAMPLE S ON A.ID_SAMPLE = S.ID_SAMPLE INNER JOIN EXPERIMENT E ON S.ID_EXPERIMENT = E.ID_EXPERIMENT WHERE ID_ANALYSIS = $fromID");
    $itemStrg = "<FONT class=\"title2\">$expName > $sampName > $anaName > Hydrophobicity > <FONT color=#DD0000>$quantifName</FONT></FONT><BR>";
} elsif ($fromItem eq 'sample') {
    ($projectID, $expName, $sampName) = $dbh->selectrow_array("SELECT E.ID_PROJECT, E.NAME, S.NAME FROM SAMPLE S INNER JOIN EXPERIMENT E ON S.ID_EXPERIMENT = E.ID_EXPERIMENT WHERE ID_SAMPLE = $fromID");
    $itemStrg = "<FONT class=\"title2\">$expName > $sampName > Hydrophobicity > <FONT color=#DD0000>$quantifName</FONT></FONT><BR>";
} elsif ($fromItem eq 'experiment') {
    ($projectID, $expName) = $dbh->selectrow_array("SELECT ID_PROJECT, NAME FROM EXPERIMENT WHERE ID_EXPERIMENT = $fromID");
    $itemStrg = "<FONT class=\"title2\">$expName > Hydrophobicity > <FONT color=#DD0000>$quantifName</FONT></FONT><BR>";
} else {
    $projectID = &promsMod::getProjectID($dbh, $quantifID, 'QUANTIFICATION');
    $itemStrg = "<FONT class=\"title2\">Hydrophobicity > <FONT color=#DD0000>$quantifName</FONT></FONT><BR>";
}

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
            selectedButton = document.getElementById('summary');
        }
    }
    if (currentButton) {
        currentButton.className = '';
    }
    currentButton = selectedButton;
    currentButton.className = 'selectedButton';
    var action = currentButton.id;
    if (action == 'summary' \|\| action == 'edit') {
        top.promsFrame.resultFrame.location = "$promsPath{cgi}/displayHydro.cgi?ACT="+action+"&ID=$quantifID&FROM=$fromItem:$fromID";
    } else if (action == 'delete') {
        if(confirm('Delete Hydrophobicity computation $quantifName ?')){
            top.promsFrame.resultFrame.location = "$promsPath{cgi}/displayHydro.cgi?ACT="+action+"&ID=$quantifID&FROM=$fromItem:$fromID";
        }
    } else if (action == 'displayHydro') {
        top.promsFrame.resultFrame.location = "$promsPath{cgi}/displayHydro.cgi?ACT=display&ID=$quantifID&FROM=$fromItem:$fromID";
    } else if (action == 'monitorJobs') {
        var monitorWindow = window.open("$promsPath{cgi}/monitorJobs.cgi?filterType=Quantification [HYDRO]&filterDateNumber=1&filterDateType=DAY&filterStatus=Queued&filterStatus=Running&filterStatus=Error&filterProject=$projectID",'monitorJobsWindow','width=1200,height=500,scrollbars=yes,resizable=yes');
    }
}
var currentButton; // currently selected button
</SCRIPT>
</HEAD>
<BODY style="margin:0px;" onload="selectOption()">
<DIV class="navigation navBorder2">
$itemStrg
<TABLE cellpadding=0 cellspacing=0>
<TR align=center valign=top>
|;

print "<TD nowrap>", &summaryOption;
unless ($projectGuestAccess) {
    print "<BR>", &editOption, "</TD>\n";
    print "<TD nowrap>", &deleteOption;
}
print "</TD>\n";
print "<TD nowrap>", &displayHydroButton, "</TD>\n";
print "<TD nowrap>", &monitorJobsButton, "</TD>\n";
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

sub displayHydroButton {
    my $disabled = ($status == 1)? " " : " disabled";
    return "<INPUT type=\"button\" id=\"displayHydro\" value=\"Display Hydrophobicity\" onclick=\"selectOption(this)\"$disabled>";
}

sub monitorJobsButton {
    return "<INPUT type=\"button\" id=\"monitorJobs\" value=\"Monitor Jobs\" onclick=\"selectOption(this)\">";
}


####>Revision history<####
# 1.0.0 Created, forked from selectOptionGSEA.cgi (VL 16/02/21)
