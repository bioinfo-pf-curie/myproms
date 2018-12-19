#!/usr/local/bin/perl -w

################################################################################
# managePathwayAnalysis.cgi       1.0.4                                        #
# Authors: P. Poullet, S.Liva (Institut Curie)                                 #
# Contact: myproms@curie.fr                                                    #
# summary, edit, delete a pathway analysis                                     #
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
$|=1;
use strict;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use POSIX qw(strftime); # to get the time
use File::Path qw(rmtree);

#print header(-charset=>'utf-8');warningsToBrowser(1);#DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $dbh=&promsConfig::dbConnect;
my $action=param('ACT')? param('ACT') : 'summary';
my $pathwayID=param('ID');
my $projectID=&promsMod::getProjectID($dbh,$pathwayID,'PATHWAY_ANALYSIS');
my ($expID, $catID, $name, $desc, $param, $status, $recDate, $upDate, $user)=$dbh->selectrow_array("SELECT ID_EXPERIMENT, ID_CATEGORY, NAME, DES, PARAM_STRG, STATUS, RECORD_DATE, UPDATE_DATE, UPDATE_USER from PATHWAY_ANALYSIS where ID_PATHWAY_ANALYSIS=$pathwayID");
$desc='' unless $desc;

if ($action eq "delete") {
    print header(-'content-encoding'=>'no',-charset=>'utf-8');
    warningsToBrowser(1); # DEBUG

    if (!$catID) {
	my $sthDeletePathAnaAna=$dbh->do("DELETE from PATHWAYANA_ANALYSIS where ID_PATHWAY_ANALYSIS=$pathwayID");
	$dbh->commit;
    }
    my $sthDeletePathway=$dbh->do("DELETE from PATHWAY_ANALYSIS where ID_PATHWAY_ANALYSIS=$pathwayID");
    $dbh->commit;
    $dbh->disconnect;

    my $pathToFile = "$promsPath{pathAna}/project_$projectID/$pathwayID";
    rmtree($pathToFile);
    print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
    parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=experiment:$expID&branchID=pathway:$expID&VIEW=functAna";
</SCRIPT>
</HEAD>
</HTML>
|;
    exit;
}

if (param('submit')) {

    my $description = param('description')? param('description') : "";
    my $pathName = param('pathName');
    my $sthUpdatePathwayAna=$dbh->do("UPDATE PATHWAY_ANALYSIS set NAME = '$pathName', DES = '$description' where ID_PATHWAY_ANALYSIS = $pathwayID") or die "Cannot prepare: " . $dbh->errstr();
    $dbh -> commit;
    $dbh -> disconnect;

    ####>Updating all frames<###
    print header(-'content-encoding'=>'no',-charset=>'utf-8');
    warningsToBrowser(1); # DEBUG
    print qq
|<HTML>
<HEAD>
<TITLE>Update All Frames</TITLE>
<SCRIPT LANGUAGE="JavaScript">
	top.promsFrame.selectedAction = 'summary';
	parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=experiment:$expID&branchID=pathwayanalysis:$pathwayID&VIEW=functAna";
</SCRIPT>
</HEAD>
</HTML>
|;
    exit;
}

my $statusStrg = ($status == -1)? "<IMG src=$promsPath{images}/bad.gif>&nbsp;&nbsp;<FONT color=\"red\">***ERROR***</FONT>&nbsp;&nbsp;<INPUT type=\"button\" value=\"open Error File\" onclick=\"openError()\">" : ($status == 0)? "<FONT color=\"orange\"><B>0n-going...</B></FONT>" : "<IMG src=$promsPath{images}/good.gif>&nbsp;<FONT color=\"green\">Finished</FONT>";

#### START HTML
print header(-charset=>'utf-8');
warningsToBrowser(1);
my $title = ($action eq 'summary')? "Pathway Analysis : <FONT color=\"red\">$name</FONT>" : "Editing Pathway Analysis : <FONT color=\"red\">$name</FONT>" ;
print qq
|<HTML>
<HEAD>
<TITLE>Pathway Analysis</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT langage="Javascript">
|;
if ($action eq 'edit') {
    print qq
|function checkForm(myForm) {
	if (!myForm.pathName.value) {
		alert('A name is expected for the analysis');
		return false;
	}
	return true;
}
|;
}
else { # summary
	if ($status == 0) {
		print qq
|//setTimeout(function(){window.location.reload(true);},5000);// Reload page every 5 seconds
setTimeout(function(){parent.optionFrame.location.reload(1);},5000);// Reload page every 5 seconds
|;
	}
}
print qq
|</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">$title</FONT><BR><BR>
|;
if ($action eq 'edit') {
    print qq
|<FORM NAME="displayMenu" method="POST" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="ID" value="$pathwayID" />
|;
}
my ($light,$dark)=&promsConfig::getRowColors;
print qq
|<TABLE border="0" width="500" cellspacing="2" cellpadding="2" bgcolor="$dark">
|;
if ($action eq 'edit') {
	print qq
|<TR><TH align="right" nowrap>Analysis name :</TH><TD bgcolor="$light"><INPUT type="text" name="pathName" size="50" maxlength="100" value="$name"/></TD></TR>
<TR><TH align="right" valign="top" nowrap>&nbsp;Description :</TH><TD bgcolor="$light"><TEXTAREA name="description" rows="2" cols="50">$desc</TEXTAREA></TD></TR>
<TR><TH colspan="2"><INPUT type="submit" name="submit" value="Update"></TD></TR>
|;
}
else {#summary
    $desc=&promsMod::HTMLcompatible($desc);
    print qq
|<TR><TH align="right" width="150" nowrap>Analysis name :</TH><TD bgcolor="$light">$name</TD></TR>
<TR><TH align="right" valign="top" nowrap>&nbsp;Description :</TH><TD bgcolor="$light">$desc</TD></TR>
<TR><TH align="right" valign="top" nowrap>&nbsp;Feature selection :</TH>
    <TD nowrap bgcolor=$light>
|;
    my ($strgFDR, $strgPVAL) = split(";", $param);
    my $FDRthreshold = (split("=",$strgFDR))[1];
    my $PVALthreshold = (split("=",$strgPVAL))[1];
    print qq|&nbsp;&bull;&nbsp;FDR=$FDRthreshold%<br>&nbsp;&bull;&nbsp;p-value=$PVALthreshold</TD>
</TR>
<TR><TH align="right" bgcolor="$dark" nowrap>Status :</TH><TD nowrap bgcolor="$light">&nbsp;&nbsp;$statusStrg</TD></TR>
<TR><TH align="right" valign="top" bgcolor="$dark" nowrap>Proteins used :</TH><TD nowrap bgcolor="$light">|;
    if ($catID) {
		my ($catName)=$dbh->selectrow_array("SELECT CONCAT(CL.NAME,' > ',CA.NAME) FROM CATEGORY CA,CLASSIFICATION CL WHERE CL.ID_CLASSIFICATION=CA.ID_CLASSIFICATION AND CA.ID_CATEGORY=$catID");
		print "&nbsp;List: $catName";
    }
    else {
		print qq|<INPUT type="button" class="font11" value="Show item list" onclick="window.location='./selectAnalyses.cgi?callType=list&ID=PATHWAY_ANALYSIS:$pathwayID'">|;
    }
	print "</TD></TR>\n";
}
print "</TABLE><BR>\n";
if ($action eq 'edit'){
    print "</FORM>\n";
}
else {
    print qq
|<FONT class="font11" style="font-weight:bold">Analysis performed by <A href="http://www.reactome.org" target="_blank">Reactome</A> (<A href="http://www.mdpi.com/2072-6694/4/4/1180" target="_blank">Milacic M. et al. Cancers 4, 2012</A>, <A href="http://nar.oxfordjournals.org/content/42/D1/D472" target="_blank">Croft D. et al. Nucleic Acid Research 42, 2013</A>.)
|;
}
print qq
|</CENTER>
</BODY>
</HTML>
|;
$dbh -> disconnect;

####>Revision history<####
# 1.0.4 Added classification name to list name (PP 07/12/18)
# 1.0.3 Change so that "Pathway Analyses" branch is selected after deletion instead of Experiment (PP 26/10/15)
# 1.0.2 Change reload window by reload optionFrame (SL 03/09/15)
# 1.0.1 Auto reload page if analysis is on-going (PP 13/11/14)
# 1.0.0 new script to manage, edit, delete pathway analysis (SL 21/10/14)
