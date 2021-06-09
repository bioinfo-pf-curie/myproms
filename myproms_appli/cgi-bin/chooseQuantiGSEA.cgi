#!/usr/local/bin/perl -w

################################################################################
# chooseQuantiGSEA.cgi       1.0.0                                             #
# Authors: V. Laigle (Institut Curie)                                          #
# Contact: myproms@curie.fr                                                    #
# Initialisation, launch and data import of Gene Set Enrichment Analysis       #
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
use strict;
use warnings;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use promsQuantif;
use POSIX qw(strftime); # to get the time
use File::Path qw(rmtree); # remove_tree
use File::Copy;


my %promsPath = &promsConfig::getServerInfo;
my %clusterInfo = &promsConfig::getClusterInfo;
my $dbh = &promsConfig::dbConnect;
my $userID = $ENV{'REMOTE_USER'};
my $projectGuestAccess;  # Must be global
my ($lightColor, $darkColor) = &promsConfig::getRowColors;

my $action = param('start')? 'launch' : param('action')? param('action') : 'choose';
my ($quantifID, $projectID);

if ($action eq 'choose') {  # GSEA launched from functional analysis pannel -> need to choose quantif
    my $expID = &promsMod::cleanNumericalParameters(param('exp'));
    $projectID = &promsMod::getProjectID($dbh, $expID, 'EXPERIMENT');
    my @userInfo = &promsMod::getUserInfo($dbh, $userID, $projectID);
    my $projectAccess = ${$userInfo[2]}{$projectID};
    $projectGuestAccess = ($projectAccess eq 'guest')? 1 : 0;

    &printQuantiChoice($expID);
    $dbh->disconnect;
    exit;
} elsif ($action eq 'launch') {
    $quantifID = &promsMod::cleanNumericalParameters(param('quantifChoice'));
    $projectID = &promsMod::getProjectID($dbh, $quantifID, 'QUANTIFICATION');
    &callStartGSEA($quantifID);
    $dbh->disconnect;
    exit;
}


sub printQuantiChoice {  # Globals: %promsPath, $dbh, $projectGuestAccess, $lightColor, $darkColor
    my ($expID) = @_;
    my %quantifTree;

    print header(-'content-encoding' => 'no',-charset=>'utf-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<TITLE>Choose quantification for GSEA</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT type="text/javascript">
function checkForm(myForm){
    if (!myForm.quantifChoice.value) {
        alert('You must choose a quantification to perform GSEA.');
        return false;
    }
    return true;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Starting a Gene Set Enrichment Analysis</FONT>
<BR><BR>
|;

    if ($projectGuestAccess) {  # Guests cannot start GSEA
        print qq
|<BR>
<FONT class="title2">Your guest status does not allow you to start a GSEA process.<BR>You may only check the results of existing GSEAs.</FONT>
|;
    } else {
        my $sthGetTree = $dbh->prepare(
            "SELECT D.ID_DESIGN, D.NAME, Q.ID_QUANTIFICATION, Q.NAME
            FROM DESIGN D
            INNER JOIN QUANTIFICATION Q ON D.ID_DESIGN = Q.ID_DESIGN
            INNER JOIN QUANTIFICATION_METHOD QM ON Q.ID_QUANTIFICATION_METHOD = QM.ID_QUANTIFICATION_METHOD
            WHERE D.ID_EXPERIMENT = ? AND QM.CODE = 'PROT_RATIO_PEP' AND Q.STATUS >= 1"
        );
        $sthGetTree->execute($expID);
        while (my ($designID, $designName, $quantifID, $quantifName) = $sthGetTree->fetchrow_array) {
            $quantifTree{$designID}{'DESIGN_NAME'} = $designName;
            $quantifTree{$designID}{$quantifID} = $quantifName;
        }
        $sthGetTree->finish;

        my $selStrg = "<OPTION value=\"\">-= Select =-</OPTION>\n";
        foreach my $designID (sort{&promsMod::sortSmart($quantifTree{$a}{'DESIGN_NAME'}, $quantifTree{$b}{'DESIGN_NAME'})} keys %quantifTree) {
            $selStrg .= "<OPTGROUP label=\"$quantifTree{$designID}{'DESIGN_NAME'}\">\n";
            foreach my $quantifID (sort{&promsMod::sortSmart($quantifTree{$designID}{$a}, $quantifTree{$designID}{$b})} keys %{$quantifTree{$designID}}) {
                next if ($quantifID eq 'DESIGN_NAME');
                $selStrg .= "<OPTION value=\"$quantifID\">$quantifTree{$designID}{$quantifID}</OPTION>\n";
            }
            $selStrg .= "</OPTGROUP>\n";
        }

        print qq
|<FORM name="gseaChoiceForm" id="gseaChoiceForm" method="post" onsubmit="return(checkForm(this));">
<BR>
<FONT class="title2">First, choose the quantification on which GSEA should be performed</FONT>
<TABLE bgcolor=$darkColor border=0 width=500>
    <TR>
        <TH valign='top' align='center'>&nbsp;&nbsp;Quantification to perform GSEA on</TH>
    </TR>
    <TR>
        <TD bgcolor=$lightColor align=center>
            <SELECT name="quantifChoice" id="quantifChoice">
                $selStrg
            </SELECT>
        </TD>
    </TR>
    <TR>
        <TD align="center">
            <INPUT type="submit" name="start" value="Start" />
            &nbsp;&nbsp;&nbsp;&nbsp;
            <INPUT type="reset" value="Clear" />
        </TD>
    </TR>
</TABLE>
</FORM>
|;
    }

    print qq
|</CENTER>
</BODY>
</HTML>
|;
}


sub callStartGSEA {
    my ($quantifID) = @_;
    
    print header(-'content-encoding' => 'no',-charset=>'utf-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<TITLE>GSEA start on choosen Quantification</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">GSEA start on choosen Quantification</FONT>

</CENTER>
<SCRIPT type="text/javascript">
    top.promsFrame.resultFrame.location="$promsPath{cgi}/startGSEA.cgi?quanti=$quantifID";
</SCRIPT>
</BODY>
</HTML>
|;
}

####>Revision history<####
# 1.0.0 Created (VL 09/11/20)