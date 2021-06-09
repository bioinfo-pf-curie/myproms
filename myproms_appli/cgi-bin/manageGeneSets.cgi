#!/usr/local/bin/perl -w

################################################################################
# manageGeneSets.cgi         1.0.0                                             #
# Authors: V. Laigle (Institut Curie)                                          #
# Contact: myproms@curie.fr                                                    #
# List Gene Sets available in myProMS and allow to add/edit/delete them        #
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
use strict;
use warnings;
use promsConfig;
use promsMod;
use File::Basename;

#####################
### Configuration ###
#####################
my %promsPath = &promsConfig::getServerInfo;
my $userID = $ENV{'REMOTE_USER'};
my $dbh = &promsConfig::dbConnect;
my ($lightColor, $darkColor) = &promsConfig::getRowColors;
my $action = param("ACT") || 'display';

# Check which gene sets the user has access to. If not massist or bioinfo, grant
# access to public gene sets + those from projects where he/she is at least administrator
my ($userStatus, $adminProjects) = $dbh->selectrow_array(
    "SELECT UL.USER_STATUS, GROUP_CONCAT(PA.ID_PROJECT SEPARATOR ',')
    FROM USER_LIST UL
    LEFT JOIN PROJECT_ACCESS PA ON UL.ID_USER = PA.ID_USER
    LEFT JOIN USER_PROFILE UP ON PA.ID_PROFILE = UP.ID_PROFILE
    WHERE UL.ID_USER = '$userID' AND (UP.NAME IS NULL OR UP.NAME LIKE '\%administrator\%')
    GROUP BY UL.ID_USER"
);

if ($action eq 'display') {
    &displayGeneSets;
}

sub displayGeneSets {  # Globals: $dbh, $userStatus, $adminProjects, $lightColor, $darkColor
    my (%geneSetsInfo, %geneSetsNotUsed, %geneSetsError);
    my $bgColor = $lightColor;
    my $projectStrg = ($userStatus =~ /mass|bioinfo/)? " " : " WHERE GS.ID_PROJECT IN ($adminProjects)";  # Exclude gene sets from projects where user is not at least administrator

    my $sthGetGSInfo = $dbh->prepare(
        "SELECT GS.ID_GENESETS, GS.NAME, GS.DES, GS.ID_PROJECT, GS.GENESETS_FILE, GS.GENESETS_TYPE, GS.NUM_SETS, GS.GENESETS_STATUS, GS.VERSION_DATE,
        GS.OTHER_VERSIONS, S.COMMON_NAME, S.SCIENTIFIC_NAME, I.NAME, P.NAME
        FROM GENESETS GS
        LEFT JOIN SPECIES S ON GS.ID_SPECIES = S.ID_SPECIES
        LEFT JOIN IDENTIFIER I ON GS.ID_IDENTIFIER = I.ID_IDENTIFIER
        LEFT JOIN PROJECT P ON GS.ID_PROJECT = P.ID_PROJECT
        $projectStrg"
    );
    $sthGetGSInfo->execute;

    while (my ($geneSetsID, $gsName, $gsDes, $gsProjectID, $gmtFile, $gsFileType, $numSets, $gsStatus, $gsVersion, $gsOtherVersions, $spCName, $spSName, $identifier, $gsProjectName) = $sthGetGSInfo->fetchrow_array) {
        $gsOtherVersions =~ s/;/, /g if ($gsOtherVersions);
        my %gsInfo = (
            'NAME'           => $gsName,
            'DESCRIPTION'    => ($gsDes)? $gsDes : '-',
            'ID_PROJECT'     => ($gsProjectID)? $gsProjectID : 0,
            'PROJECT_NAME'   => ($gsProjectName)? $gsProjectName : '',
            'GENESETS_FILE'  => basename($gmtFile),
            'GENESETS_TYPE'  => lc($gsFileType),
            'NUM_SETS'       => $numSets,
            'STATUS'         => ($gsStatus == 1)? "In use" : ($gsStatus == 0)? "Not in use" : "Error",
            'VERSION_DATE'   => $gsVersion,
            'OTHER_VERSIONS' => ($gsOtherVersions)? ($gsOtherVersions) : '',
            'SPECIES'        => ($spSName && $spCName)? "$spSName ($spCName)" :
                                ($spSName)? $spSName :
                                ($spCName)? $spCName :
                                'No species',
            'IDENTIFIER'     => ($identifier)? $identifier : 'Custom identifier'
        );

        if ($gsStatus == 1) {
            $geneSetsInfo{$geneSetsID} = \%gsInfo;
        } elsif ($gsStatus == 0) {
            $geneSetsNotUsed{$geneSetsID} = \%gsInfo;
        } else {
            $geneSetsError{$geneSetsID} = \%gsInfo;
        }
    }
    $sthGetGSInfo->finish;

    ####################
    ### Display HTML ###
    ####################
    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<TITLE>List of Gene Sets Databanks</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
    TD {font-size:13px;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
function deleteGeneSetsDB(geneSetsID, gsName) {
    if (confirm('Delete Gene Sets databank: ' + gsName + ' ?')) {
        window.location = './editGeneSets.cgi?ACT=delete&ID=' + geneSetsID;
    }
}

function updateGeneSetsDB(geneSetsID, gsName) {
    if (confirm('Update Gene Sets Databank: ' + gsName + ' ? The current databank will be archived.')) {
        window.location = './editGeneSets.cgi?ACT=update&ID=' + geneSetsID;
    }
}

function restoreGeneSetsDB(geneSetsID, gsName) {
    if (confirm('Restore another version of Gene Sets Databank: ' + gsName + ' ? The current databank will be archived.')) {
        window.location = './editGeneSets.cgi?ACT=restore&ID=' + geneSetsID;
    }
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<DIV id="backMenuDiv" style="float:top">
<BR>
<TABLE>
    <TR>
        <TH bgcolor="$darkColor">
            <FONT class="title2">&nbsp;Go to:</FONT>
            <SELECT class="title action" onchange="if (this.value) {window.location='$promsPath{cgi}/'+this.value;}">
                <OPTION value="">-= Select =-</OPTION>
                <OPTION value="promsMain.cgi">Main Window</OPTION>
                <OPTION value="selectProject.cgi">Project Selection</OPTION>
            </SELECT>
        </TH>
    </TR>
</TABLE>
</DIV>
<CENTER>
<FONT class="title">List of Gene Sets Databanks <FONT color='#DD0000'>in Use</FONT></FONT>
&nbsp;&nbsp;
<INPUT type="button" class="title2" value="Add new Gene Sets" onclick="window.location='./editGeneSets.cgi?ACT=add'">
<BR><BR>
<TABLE border=0 cellspacing=0>
|;
    if (scalar keys %geneSetsInfo == 0) {
        print "<TR><TH colspan=3>(No Gene Sets Databank in use)</TH></TR>\n";
    } else {
        foreach my $geneSetsID (sort{&promsMod::sortSmart($geneSetsInfo{$a}{'NAME'}, $geneSetsInfo{$b}{'NAME'})} keys %geneSetsInfo) {
            &printGeneSetsInfo($geneSetsID, $geneSetsInfo{$geneSetsID}, $bgColor);
            $bgColor = ($bgColor eq $lightColor)? $darkColor : $lightColor;
        }
    }

    print qq
|</TABLE>
<BR><BR>
<FONT class="title">List of Gene Sets Databanks <FONT color='#DD0000'>Not</FONT> in Use</FONT>
<BR><BR>
<TABLE border=0 cellspacing=0>
|;
    if (scalar keys %geneSetsNotUsed == 0) {
        print "<TR><TH colspan=3>(No Gene Sets Databank Not in use)</TH></TR>\n";
    } else {
        foreach my $geneSetsID (sort{&promsMod::sortSmart($geneSetsNotUsed{$a}{'NAME'}, $geneSetsNotUsed{$b}{'NAME'})} keys %geneSetsNotUsed) {
            &printGeneSetsInfo($geneSetsID, $geneSetsNotUsed{$geneSetsID}, $bgColor);
            $bgColor = ($bgColor eq $lightColor)? $darkColor : $lightColor;
        }
    }

    print qq
|</TABLE>
<BR><BR>
<FONT class="title">List of Gene Sets Databanks with <FONT color='#DD0000'>Errors</FONT></FONT>
<BR><BR>
<TABLE border=0 cellspacing=0>
|;
    if (scalar keys %geneSetsError == 0) {
        print "<TR><TH colspan=3>(No Gene Sets Databank with errors)</TH></TR>\n";
    } else {
        foreach my $geneSetsID (sort{&promsMod::sortSmart($geneSetsError{$a}{'NAME'}, $geneSetsError{$b}{'NAME'})} keys %geneSetsError) {
            &printGeneSetsInfo($geneSetsID, $geneSetsError{$geneSetsID}, $bgColor);
            $bgColor = ($bgColor eq $lightColor)? $darkColor : $lightColor;
        }
    }

    print qq
|</TABLE>
</CENTER>
</BODY>
</HTML>
|;

}


sub printGeneSetsInfo {
    my ($geneSetsID, $refGSInfo, $bgColor) = @_;

    my $projectStrg;
    if ($refGSInfo->{'PROJECT_NAME'}) {
        $projectStrg = "&nbsp;&nbsp;<FONT style=\"font-size:14px;\">(Specific to Project<FONT color='#DD0000'>&nbsp;$refGSInfo->{'PROJECT_NAME'}</FONT>)</FONT>";
    } else {
        $projectStrg = "&nbsp;&nbsp;<FONT style=\"font-size:14px;\">(Public)</FONT>";
    }

    print qq
|<TR id='gs$geneSetsID' style='margin-bottom:10px'>
    <TD>
        <TABLE bgcolor='$bgColor'>
            <TR>
                <TD id='gsInfos$geneSetsID' style='width:600px'>
                    <P style='font-size:18px;font-weight:bold;margin:3px 0 0 12px;'>
                        $refGSInfo->{'NAME'}$projectStrg
                    </P>
                    <P style='margin-left:12px;'>
                        <B>&nbsp;&nbsp;Version date:</B> $refGSInfo->{'VERSION_DATE'}
|;
    print "&nbsp;&nbsp;(Other:</B> $refGSInfo->{'OTHER_VERSIONS'})" if ($refGSInfo->{'OTHER_VERSIONS'});
    print "<BR><B>&nbsp;&nbsp;Description:</B> ", &promsMod::HTMLcompatible($refGSInfo->{'DESCRIPTION'});
    print qq
|                       <BR>
                        <B>&nbsp;&nbsp;Status:</B> $refGSInfo->{'STATUS'}
                        <BR>
                        <B>&nbsp;&nbsp;Identifier type:</B> $refGSInfo->{'IDENTIFIER'}
                        <BR>
                        <B>&nbsp;&nbsp;Species:</B> $refGSInfo->{'SPECIES'}
                        <BR>
                        <B>&nbsp;&nbsp;Gene Sets file:</B> $refGSInfo->{'GENESETS_FILE'}
                        <BR>
                        <B>&nbsp;&nbsp;File type:</B> $refGSInfo->{'GENESETS_TYPE'}
                        <BR>
                        <B>&nbsp;&nbsp;Number of Sets:</B> $refGSInfo->{'NUM_SETS'}
                    </P>
                </TD>
                <TD id="gsActions$geneSetsID" style='width:105px;'>
                    <DIV id="edit$geneSetsID" style="margin:2px 0 0 0;">
                        <INPUT type="button" value="Edit" style="width:100px" onclick="window.location='./editGeneSets.cgi?ACT=edit&ID=$geneSetsID'">
                    </DIV>
                    <DIV id="update$geneSetsID" style="margin:2px 0 0 0;">
                        <INPUT type="button" value="Update" style="width:100px" onclick="updateGeneSetsDB($geneSetsID, '$refGSInfo->{NAME}')">
                    </DIV>
|;
    if ($refGSInfo->{'OTHER_VERSIONS'}) {
        print qq
|                   <DIV id="restore$geneSetsID" style="margin:2px 0 0 0;">
                        <INPUT type="button" value="Restore" style="width:100px" onclick="restoreGeneSetsDB($geneSetsID, '$refGSInfo->{NAME}')">
                    </DIV>
|;
    }

    print qq
|                    <DIV id="delete$geneSetsID" style="margin:2px 0 0 0;">
                        <INPUT type="button" value="Delete" style="width:100px" onclick="deleteGeneSetsDB($geneSetsID, '$refGSInfo->{NAME}')">
                    </DIV>
                </TD>
            </TR>
        </TABLE>
    </TD>
</TR>
|;
}


####>Revision history<####
# 1.0.0 Created (VL 10/11/20)
