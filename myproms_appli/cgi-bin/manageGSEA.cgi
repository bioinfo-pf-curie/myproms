#!/usr/local/bin/perl -w

################################################################################
# manageGSEA.cgi       1.0.0                                                   #
# Authors: V. Laigle (Institut Curie)                                          #
# Contact: myproms@curie.fr                                                    #
# summary, edit, delete a Gene Set Enrichment Analysis                         #
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
use warnings;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use File::Path qw(rmtree);

#####################
### Configuration ###
#####################
my %promsPath = &promsConfig::getServerInfo;
my $dbh = &promsConfig::dbConnect;
my $action = param('ACT')? param('ACT') : 'summary';
my $gseaID = param('ID');
my $projectID = &promsMod::getProjectID($dbh, $gseaID, 'PATHWAY_ANALYSIS');  # GSEA in table PATHWAY_ANALYSIS
my ($expID, $name, $desc, $param, $status, $recDate, $upDate, $user) = $dbh->selectrow_array("SELECT ID_EXPERIMENT, NAME, DES, PARAM_STRG, STATUS, RECORD_DATE, UPDATE_DATE, UPDATE_USER from PATHWAY_ANALYSIS where ID_PATHWAY_ANALYSIS = $gseaID");
$desc='' unless $desc;

if ($action eq "delete") {
    print header(-'content-encoding'=>'no',-charset=>'utf-8');
    warningsToBrowser(1);

    $dbh->do("DELETE from PATHWAYANA_QUANTIFICATION where ID_PATHWAY_ANALYSIS = $gseaID");
    $dbh->do("DELETE from PATHWAY_ANALYSIS where ID_PATHWAY_ANALYSIS = $gseaID");
    $dbh->commit;
    $dbh->disconnect;

    my $pathToFile = "$promsPath{'gsea'}/project_$projectID/gsea_$gseaID";
    rmtree($pathToFile);

    print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
    parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=experiment:$expID&branchID=gsea:$expID&VIEW=functAna";
</SCRIPT>
</HEAD>
</HTML>
|;
    exit;
}

if (param('submit')) {

    my $description = param('description')? param('description') : "";
    my $gseaName = param('gseaName');
    $dbh->do("UPDATE PATHWAY_ANALYSIS set NAME = '$gseaName', DES = '$description' where ID_PATHWAY_ANALYSIS = $gseaID") or die "Cannot prepare: " . $dbh->errstr();
    $dbh -> commit;
    $dbh -> disconnect;

    ####>Updating all frames<###
    print header(-'content-encoding'=>'no',-charset=>'utf-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<TITLE>Update All Frames</TITLE>
<SCRIPT LANGUAGE="JavaScript">
    top.promsFrame.selectedAction = 'summary';
    parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=experiment:$expID&branchID=gseanalysis:$gseaID&VIEW=functAna";
</SCRIPT>
</HEAD>
</HTML>
|;
    exit;
}

my $statusStrg;
if ($status == -1) {
    $statusStrg = "<IMG src=$promsPath{images}/bad.gif>";
    $statusStrg .= "&nbsp;&nbsp;<FONT color=\"red\">***ERROR***</FONT>";
    $statusStrg .= "&nbsp;&nbsp;(Check what happened with the \"Monitor Jobs\" button at the top)"
} elsif ($status == 0) {
    $statusStrg = "<FONT color=\"orange\"><B>0n-going...</B></FONT>";
} else {
    $statusStrg = "<IMG src=$promsPath{images}/good.gif>";
    $statusStrg .= "&nbsp;<FONT color=\"green\">Finished</FONT>";
}

# START HTML
print header(-charset=>'utf-8');
warningsToBrowser(1);

my $title = ($action eq 'summary')? "Gene Set Enrichment Analysis : <FONT color=\"red\">$name</FONT>" : "Editing Gene Set Enrichment Analysis : <FONT color=\"red\">$name</FONT>" ;
print qq
|<HTML>
<HEAD>
<TITLE>Gene Set Enrichment Analysis</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT langage="Javascript">
|;

if ($action eq 'edit') {
    print qq
|function checkForm(myForm) {
    if (!myForm.gseaName.value) {
        alert('A name is expected for the analysis');
        return false;
    }
    return true;
}
function cancelAction() {
    top.promsFrame.selectedAction='summary';
    top.promsFrame.optionFrame.selectOption();
}
|;

} else {  # summary
    if ($status == 0) {
        print qq
|setTimeout(function(){parent.optionFrame.location.reload(1);},5000);  // Reload page every 5 seconds
|;
    }
}

print qq
|</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">$title</FONT>
<BR><BR>
|;

if ($action eq 'edit') {
    print qq
|<FORM NAME="displayMenu" method="POST" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="ID" value="$gseaID" />
|;
}

my ($light, $dark)=&promsConfig::getRowColors;

print qq
|<TABLE border="0" width="500" cellspacing="2" cellpadding="2" bgcolor="$dark">
|;

if ($action eq 'edit') {
    print qq
|<TR>
    <TH align="right" nowrap>GSEA name :</TH>
    <TD bgcolor="$light">
        <INPUT type="text" name="gseaName" size="50" maxlength="100" value="$name"/>
    </TD>
</TR>
<TR>
    <TH align="right" valign="top" nowrap>&nbsp;Description :</TH>
    <TD bgcolor="$light">
        <TEXTAREA name="description" rows="2" cols="50">$desc</TEXTAREA>
    </TD>
</TR>
<TR>
    <TH colspan="2">
        <INPUT type="submit" name="submit" value="Update">
        &nbsp;&nbsp;
        <INPUT type="button" value="Cancel" onclick="cancelAction()">
    </TH>
</TR>
|;

} else {  # summary
    $desc = &promsMod::HTMLcompatible($desc);
    my ($quantifID, $targetPos, $isRatio, $geneSetsID, $gsVersion);
    my ($speciesID, $pepWeights, $pvalCutoff, $pvalCorrection, $maxRelevantGS);
    my ($quantifStrg, $geneSetsStrg, $speciesStrg, $pvalCorrStrg);
    foreach my $feature (split(';', $param)) {
        $feature =~ s/#//g;  # Clean ID tags
        if ($feature =~ /^quantification=(\d+)$/) {
            $quantifID = $1;
        } elsif ($feature =~ /^targetPos=(.*)$/) {
            $targetPos = $1;
        } elsif ($feature =~ /^isRatio=(\d)$/) {
            $isRatio = $1;
        } elsif ($feature =~ /^species=(.+)$/) {  # ID_SPECIES but may also equal "Custom"
            $speciesID = $1;
        } elsif ($feature =~ /^geneSetsDB=(\d+)$/) {
            $geneSetsID = $1;
        } elsif ($feature =~ /^geneSetsVersion=(.+)$/) {
            $gsVersion = $1;
        } elsif ($feature =~ /^weightByPep=(.+)$/) {
            $pepWeights = $1;
        } elsif ($feature =~ /^pvalCutoff=(.+)$/) {
            $pvalCutoff = $1;
        } elsif ($feature =~ /^pvalCorr=(.+)$/) {
            $pvalCorrection = $1;
        } elsif ($feature =~ /^maxGS=(\d+)$/) {
            $maxRelevantGS = $1;
        }
    }

    # Info about parent quantification
    my ($designName, $quantifName, $quantifAnnot) = $dbh->selectrow_array("SELECT D.NAME, Q.NAME, Q.QUANTIF_ANNOT FROM QUANTIFICATION Q INNER JOIN DESIGN D ON Q.ID_DESIGN = D.ID_DESIGN WHERE ID_QUANTIFICATION = $quantifID");
    my $sthExpCond = $dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION = ?");

    $quantifStrg = "$designName > $quantifName";

    if ($isRatio) {  # Only case possible for now, prevision for GSEA on one state (after abundance quantif)
        my ($testCondID, $refCondID);
        my ($testCondName, $refCondName);
        my $supRatioTag = '';  # In case of SuperRatio 
        foreach my $features (split('::', $quantifAnnot)) {
            if ($features =~ /^RATIOS=(.*)$/) {
                my $ratioStrg = $1;
                $ratioStrg =~ s/#//g;  # Clean ID tags
                my @allRatios = split(';', $ratioStrg);
                my $interestRatio = $allRatios[$targetPos - 1];  # Retrieve state IDs of the wanted ratio
                ($testCondID, $refCondID) = split('/', $interestRatio);
                if ($testCondID =~ /%/) {
                    $testCondID =~ s/%.+//;
                    $refCondID =~ s/%.+//;
                    $supRatioTag = '°';
                }
                last;
            }
        }
        $sthExpCond->execute($testCondID);
        ($testCondName) = $sthExpCond->fetchrow_array;
        $sthExpCond->execute($refCondID);
        ($refCondName) = $sthExpCond->fetchrow_array;

        $quantifStrg .= " : ${testCondName}${supRatioTag}/${refCondName}${supRatioTag}";
    }
    $sthExpCond->finish;

    # Weights on protein scores/ratios ?
    my $pepWeightsStrg = "Protein scores are ";
    if ($pepWeights eq "none") {
        $pepWeightsStrg .= "not weighted with peptides number";
    } else {
        $pepWeightsStrg .= "weighted with number of <FONT color=\"red\"><B>";
        if ($pepWeights eq "msms") {
            $pepWeightsStrg .= "Truly identified";
        } else {
            $pepWeightsStrg .= ucfirst($pepWeights);
        }
        $pepWeightsStrg .= "</B></FONT> peptides";
    }

    # Info about Gene Sets databank
    my ($gsName, $gsDes, $gsProject, $gsSpecies, $gsIdentifier) = $dbh->selectrow_array(
        "SELECT GS.NAME, GS.DES, GS.ID_PROJECT, S.SCIENTIFIC_NAME, I.NAME
        FROM GENESETS GS 
        LEFT JOIN SPECIES S ON GS.ID_SPECIES = S.ID_SPECIES
        LEFT JOIN IDENTIFIER I ON GS.ID_IDENTIFIER = I.ID_IDENTIFIER
        WHERE ID_GENESETS = $geneSetsID"
    );
    $gsSpecies    = "No specific organism" unless ($gsSpecies);
    $gsIdentifier = "No specific identifier" unless ($gsIdentifier);

    # Make p-value correction readable
    my %fullPvalCorr = (
        "none"       => "No p-value correction",
        "BH"         => "Benjamini-Hochberg",
        "BY"         => "Benjamini-Yekutieli",
        "hochberg"   => "Hochberg",
        "bonferroni" => "Bonferroni",
        "holm"       => "Holm",
        "hommel"     => "Hommel"
    );

    # Creation date and user info
    $upDate=&promsMod::formatDate($upDate);
    if ($user) {
        my ($userName) = $dbh->selectrow_array("SELECT USER_NAME FROM USER_LIST WHERE ID_USER = '$user'");
        $user = $userName || $user;
    } else {
        $user = '?';
    }

    print qq
|<TR>
    <TH align="right" width="150" nowrap>GSEA name :</TH>
    <TD bgcolor="$light">&nbsp;&nbsp;$name</TD>
</TR>
<TR>
    <TH align="right" valign="top" nowrap>&nbsp;Description :</TH>
    <TD bgcolor="$light">&nbsp;&nbsp;$desc</TD>
</TR>
<TR>
    <TH align="right" valign="top" nowrap>&nbsp;Parent quantification :</TH>
    <TD nowrap bgcolor=$light>&nbsp;&nbsp;$quantifStrg</TD>
</TR>
<TR>
    <TH align="right" valign="top" nowrap>&nbsp;Gene sets databank :</TH>
    <TD bgcolor=$light>&nbsp;&nbsp;$gsName
        <BR>
        &nbsp;&nbsp;&bull;&nbsp;Version date : $gsVersion
        <BR>
        &nbsp;&nbsp;&bull;&nbsp;Organism : $gsSpecies
        <BR>
        &nbsp;&nbsp;&bull;&nbsp;Identifier type : $gsIdentifier
|;
    if ($gsDes) {
        print qq
|       <BR>
        &nbsp;&nbsp;&bull;&nbsp;Gene sets description : $gsDes
|;    
    }
    print qq
|   </TD>
</TR>
<TR>
    <TH align="right" valign="top" nowrap>&nbsp;Parameters :</TH>
    <TD nowrap bgcolor=$light>
        &nbsp;&nbsp;&bull;&nbsp;$pepWeightsStrg
        <BR>
        &nbsp;&nbsp;&bull;&nbsp;Adjusted p-value cutoff : $pvalCutoff
        <BR>
        &nbsp;&nbsp;&bull;&nbsp;p-value correction : $fullPvalCorr{$pvalCorrection}
        <BR>
        &nbsp;&nbsp;&bull;&nbsp;Max enriched Gene Sets reported : $maxRelevantGS
    </TD>
</TR>
<TR>
    <TH align="right" bgcolor="$dark" nowrap>Status :</TH>
    <TD nowrap bgcolor="$light">&nbsp;&nbsp;$statusStrg</TD>
</TR>
<TR>
    <TH align="right" bgcolor="$dark" nowrap>Creation date :</TH>
    <TD nowrap bgcolor="$light">&nbsp;&nbsp;$upDate by <B>$user</B></TD></TR>

|;
}

print qq
|</TABLE><BR>
|;

if ($action eq 'edit'){
    print "</FORM>\n";
} else {
    print qq
|<FONT class="font11" style="font-weight:bold">Analysis performed with <A href="https://guangchuangyu.github.io/software/clusterProfiler/" target="_blank">clusterProfiler</A> (<A href="https://doi.org/10.1089/omi.2011.0118" target="_blank">Yu G. et al. OMICS, 2012</A>), based on the <A href="https://www.gsea-msigdb.org/gsea/" target="_blank">GSEA method</A> (<A href="https://doi.org/10.1073/pnas.0506580102" target="_blank">Subramanian A. et al. PNAS, 2005</A>).</FONT>
|;
}
print qq
|</CENTER>
</BODY>
</HTML>
|;
$dbh -> disconnect;

####>Revision history<####
# 1.0.0 Created, forked from managePathwayAnalysis.cgi (VL 05/11/20)
