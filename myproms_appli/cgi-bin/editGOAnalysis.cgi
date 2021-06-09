#!/usr/local/bin/perl -w

################################################################################
# editGOAnalysis.cgi    1.1.7                                                  #
# Authors: P. Poullet, G. Arras & F. Yvon (Institut Curie)                     #
# Contact: myproms@curie.fr                                                    #
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
use promsQuantif;
use strict;
use goAnalysis;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};

###################
####>Arguments<####
###################
my $goAnaID=param('ID');
my $action=(param('ACT'))?param('ACT'):'summary';

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

my ($projectID,$expID)=$dbh->selectrow_array("SELECT ID_PROJECT,GO_ANALYSIS.ID_EXPERIMENT FROM EXPERIMENT,GO_ANALYSIS WHERE GO_ANALYSIS.ID_EXPERIMENT=EXPERIMENT.ID_EXPERIMENT AND ID_GOANALYSIS=$goAnaID");
my $resultDirUnix="$promsPath{go_unix}/project_$projectID/$goAnaID";

my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);

if(param('save')){
    &processForm; exit;
}
if($action eq 'delete'){
    &deleteAna; exit;
}


#####################################
####>Fetching goAna informations<####
#####################################

my ($name,$des,$aspectStrg,$paramStrg,$type,$updateDate,$updateUser,$goaName,$goaVersionDate,$ontoName,$ontoDataVersion,$ontoVersionDate) = $dbh->selectrow_array(
                "SELECT GO_ANALYSIS.NAME,GO_ANALYSIS.DES,ASPECT,PARAM_STRG,GOA_TYPE,GO_ANALYSIS.UPDATE_DATE,GO_ANALYSIS.UPDATE_USER,
                 GOANNOTATION.NAME,GOANNOTATION.VERSION_DATE,
                 ONTOLOGY.NAME,ONTOLOGY.DATA_VERSION,ONTOLOGY.VERSION_DATE
                 FROM GO_ANALYSIS, GOANNOTATION, ONTOLOGY
                 WHERE ID_GOANALYSIS=$goAnaID
                 AND GO_ANALYSIS.ID_GOANNOTATION=GOANNOTATION.ID_GOANNOTATION
                 AND GO_ANALYSIS.ID_ONTOLOGY=ONTOLOGY.ID_ONTOLOGY"
                 ) or die "Unfound GO analysis";

if(param('AJAX')){
    print header(-'content-encoding'=>'no',-charset=>'utf-8');
    warningsToBrowser(1);
    if($type eq 'heatmap'){
        my %binToId;
        my $sthChildren = $dbh->prepare("SELECT ID_GOANALYSIS, PARAM_STRG FROM GO_ANALYSIS WHERE ID_PARENT_GOANA=?");
        my $sthProt = $dbh->prepare("SELECT ALIAS FROM PROTEIN WHERE ID_PROTEIN=?");
        $sthChildren->execute($goAnaID);
        while(my ($childID, $paramStrg) = $sthChildren->fetchrow_array){
            my ($bin) = ($paramStrg =~ /bin=(\d)/);
            $binToId{$bin} = $childID;
        }
        $sthChildren->finish;
        foreach my $bin (sort {$a <=> $b } keys %binToId){
            print '<B>Bin '.($bin+1).":</B><BR>\n";
            getWarnings($projectID, $binToId{$bin});
            print "<BR>\n";
            # Checking unannotated proteins #
            my @unannotatedProteins;
            my $resultDirUnix="$promsPath{go_unix}/project_$projectID/$binToId{$bin}";
            open RESULTS, "$resultDirUnix/results_$aspectStrg.txt" or do {print "Cannot read result file<BR>\n"; next};
            while (<RESULTS>) {
                next if /^###/;
                my @infos = split /\t/;

                if ($infos[0] eq 'unannotated') {
                    foreach my $protID (split /;/, $infos[5]){
                        $sthProt->execute($protID);
                        my ($protName) = $sthProt->fetchrow_array;
                        push @unannotatedProteins, $protName;
                    }
                }
            }
            close RESULTS;
            if (scalar @unannotatedProteins) {
                print "Unannotated proteins: ";
                print join ', ', @unannotatedProteins;
                print "<BR><BR>\n";
            }

        }
        $sthProt->finish;
    } else {
        getWarnings($projectID, $goAnaID)
    }
    $dbh->disconnect;
    exit;
}

my %parameters;
foreach my $param (split(';',$paramStrg)){
    my ($key,$value) = ($param =~ /([^=]+)=([^=]+)/);
    $parameters{$key} = $value;
}

my ($light,$dark)=&promsConfig::getRowColors;
my %aspectName = ( 'P' => 'Biological Process',
                  'C' => 'Cellular Component',
                  'F' => 'Molecular Function'
                  );

if($action eq 'summary'){
    ####>Format value list<####
    my @entries = ( ['Name', $name], ['Description', $des]);
    # Aspects
    my $aspectNameStrg = '';
    foreach my $aspect (split(//,$aspectStrg)){
        $aspectNameStrg.=$aspectName{$aspect}.", ";
    }
    $aspectNameStrg =~ s/\, \Z//;
    push @entries, ['Aspect', $aspectNameStrg];
    $entries[-1]->[0] .= 's' if length($aspectStrg) > 1;
    # Parameters
    my $paramValStrg = '';
    if (defined $parameters{criteria}) {
        if ($parameters{criteria} eq 'FDR'){
            my $FDRmethod = ($parameters{method} eq 'BH') ? 'Benjamini & Hochberg' : ($parameters{method} eq 'FDRsimulation') ? 'simulation' : $parameters{method};
            $paramValStrg .= "&bull;FDR &le; $parameters{threshold}% ($FDRmethod)<BR>\n";
        }
        elsif ($parameters{criteria} eq 'pval'){
            $paramValStrg .= "&bull;p-value threshold: $parameters{threshold}<BR>\n";
            $paramValStrg .= "&bull;Multitest correction: $parameters{method}<BR>\n";
        }
    }
    if ($parameters{showAllNodes}){$parameters{showAllNodes} = 'Yes'}else{$parameters{showAllNodes} = 'No'}
    $paramValStrg .= "&bull;Show non-significant terms: $parameters{showAllNodes}<BR>\n";
    if (defined $parameters{minPep}){ $paramValStrg .= "&bull;Include only proteins containing at least $parameters{minPep} peptide(s)<BR>\n"}
    if (defined $parameters{logRatios}) {
        $paramValStrg .= "&bull;log2-ratios thresholds: " . join(', ',split(',',$parameters{logRatios})) . "<BR>\n";
    }
    if (defined $parameters{ratioPvalue}) {
        $paramValStrg .= "&bull;Only keep ratios with a (adj.) p-value &le; $parameters{ratioPvalue}<BR>\n";
    }
    $paramValStrg =~ s/<BR>\n*\Z//;
    push @entries, ['Parameters', $paramValStrg];
    push @entries, ['Ontology File', "$ontoName (v. $ontoDataVersion, $ontoVersionDate)"];
    $entries[-1]->[1] =~ s/v\. N\/A/unknown version/;
    if($parameters{'depth'}){
        push @entries, ['Depth', $parameters{'depth'}];
    }
    push @entries, ['Annotation File', "$goaName (v. $goaVersionDate)"];
    # Protein set #
    if ($type eq 'heatmap') {
        my ($quantifID,$ratioPos,$quantifModID,$multiModifStrg) = $dbh->selectrow_array("SELECT Q.ID_QUANTIFICATION,RATIO,Q.ID_MODIFICATION,GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ',')
                                                                        FROM QUANTIFICATION Q
                                                                        LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
                                                                        INNER JOIN GOANA_QUANTIFICATION G ON G.ID_QUANTIFICATION=Q.ID_QUANTIFICATION
                                                                        WHERE ID_GOANALYSIS=$goAnaID GROUP BY Q.ID_QUANTIFICATION");
        my %quantifAllNames=&promsQuantif::getDistinctQuantifNames($dbh,$quantifID.'_'.$ratioPos);
        my $modifStrg='';
        if ($quantifModID) {
            my ($modifName)=$dbh->selectrow_array("SELECT PSI_MS_NAME FROM MODIFICATION WHERE ID_MODIFICATION=$quantifModID");
            $modifStrg=' ['.$modifName.'-sites]';
        }
        elsif ($multiModifStrg) {$modifStrg=' [Multi-sites]';}
        push @entries, ['Analyzed set', 'Quantification '.$quantifAllNames{'FULL'}{$quantifID.'_'.$ratioPos}.$modifStrg];
    }
    else {
        if (defined $parameters{protSetCategory}){
            my $categoryID = $parameters{protSetCategory};
            my ($categoryName) = $dbh->selectrow_array("SELECT NAME FROM CATEGORY WHERE ID_CATEGORY=$categoryID");
            push @entries, ['List', "$categoryName&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Show protein list\" onclick=\"window.location='./listProteins.cgi?TYPE=EXPERIMENT&ID=$expID&listMode=child&listFilter=category:$categoryID&view=peptide&expMode=1&pepType=ba&unClass=0'\">"];
        }
        else {
            push @entries, ['List', "<INPUT type=\"button\" class=\"font11\" value=\"Show item list\" onclick=\"window.location='./selectAnalyses.cgi?callType=list&ID=GO_ANALYSIS:$goAnaID'\">"];
        }
    }

    # Population #
    if (defined $parameters{popCategory}){
        my $categoryID = $parameters{popCategory};
        if ($categoryID > 0){
            my ($categoryName) = $dbh->selectrow_array("SELECT NAME FROM CATEGORY WHERE ID_CATEGORY=$categoryID");
            push @entries,  ['Population', "$categoryName&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Show protein list\" onclick=\"window.location='./listProteins.cgi?TYPE=EXPERIMENT&ID=$expID&listMode=child&listFilter=category:$categoryID&view=peptide&expMode=1&pepType=ba&unClass=0'\">"];
        } else { # =0
            push @entries, ['Population', 'Unspecified'];
        }
    } elsif (defined $parameters{popFile}) {
        push @entries, ['Population', "Local file: $parameters{popFile}"];
    } elsif (defined $parameters{pop} && $parameters{pop} eq 'allannotatedgenes') {
        push @entries, ['Population', "all annotated genes"];
    } elsif ($type eq 'heatmap') {
        my $numProtStrg='';
        if (defined $parameters{numProtSel}) {
            $numProtStrg="$parameters{numSiteSel} sites, " if defined $parameters{numSiteSel};
            $numProtStrg.="$parameters{numProtSel} proteins"; # Added PP 09/02/18
            $numProtStrg=' ('.$numProtStrg.')';
        }
        push @entries, ['Population', "All merged bins$numProtStrg"];
    }

    # Unannotated proteins #
    if (defined $parameters{unannotated}) {
        $parameters{unannotated} =~ /^(\d+.\d+),(\d+.\d+),(\d+.\d+),(\d+.\d+),(\d+.\d+)/;
        push @entries, ['Unannotated proteins',"Bin1=$1% Bin2=$2% Bin3=$3% Bin4=$4% Bin5=$5%"];
    }

    # Checking for warnings in error file
    my $errorFileName = "$resultDirUnix/errorLog.txt";
    if (! -z $errorFileName) { # file contains warnings
        my $warningString = "&nbsp;<INPUT id='showWarnsButton' type='button' class=\"font11\" value=\"Show\" onclick=\"getWarnings(this);\"><DIV id='warnings' style='max-height:400px;overflow:auto;display:none'></DIV>";
        push @entries, ['Warnings', $warningString];
    }
    my ($userName)=$dbh->selectrow_array("SELECT USER_NAME FROM USER_LIST WHERE ID_USER='$updateUser'") || $updateUser;
    push @entries, ['Last changed',"On ".&promsMod::formatDate($updateDate)." by $userName"];

    $dbh->disconnect;

    ####>HTML<####
    print header(-'content-encoding'=>'no',-charset=>'utf-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT language="Javascript">
// AJAX --->
var XHR=null;
function ajaxGetWarnings(warningDiv) {
        //wait image
        warningDiv.innerHTML='<IMG src="$promsPath{images}/scrollbarGreen.gif">';
        warningDiv.style.display='block';

        //If XHR object already exists, the request is canceled & the object is deleted
        if(XHR && XHR.readyState != 0){
                XHR.abort();
                delete XHR;
        }

        //Creation of the XMLHTTPRequest object
        XHR = getXMLHTTP();
        if (!XHR) {
                return false;
        }
        XHR.open("GET","$promsPath{cgi}/editGOAnalysis.cgi?ID=$goAnaID&AJAX=1",true);
        XHR.onreadystatechange=function() {
                if (XHR.readyState==4 && XHR.responseText){
                        warningDiv.innerHTML = XHR.responseText;
                        warningDiv.scrollIntoView();
                }
        }
        XHR.send(null);
}
function getXMLHTTP(){
        var xhr=null;
        if(window.XMLHttpRequest) {// Firefox & others
                xhr = new XMLHttpRequest();
        }
        else if(window.ActiveXObject){ // Internet Explorer
                try {
                  xhr = new ActiveXObject("Msxml2.XMLHTTP");
                } catch (e) {
                        try {
                                xhr = new ActiveXObject("Microsoft.XMLHTTP");
                        } catch (e1) {
                                xhr = null;
                        }
                }
        }
        else { // XMLHttpRequest not supported by browser
                alert("Your browser does not support XMLHTTPRequest objects...");
        }
        return xhr;
}
// <--- AJAX
function getWarnings(warningButton){
    var warningDiv=document.getElementById("warnings");
    if(warningButton.value == 'Show'){
        ajaxGetWarnings(warningDiv);
        warningButton.value = 'Hide';
    } else if (warningButton.value == 'Hide'){
        warningDiv.style.display = 'none';
        warningButton.value = 'Show';
    }
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER><FONT class="title1">GO Analysis <FONT color=red>$name</FONT></FONT><BR><BR>
<TABLE border=0 width=600>
<TR><TD bgcolor=$dark>
<TABLE width=100% cellpadding=2 border=0>
|;
    foreach my $entry (@entries){
        print "<TR><TH align=right valign=top nowrap>&nbsp;",$entry->[0]," : </TH><TD bgcolor=$light>",$entry->[1]," &nbsp;</TD></TR>\n";
    }
    print qq
|</TABLE>
</TD></TR></TABLE>
</CENTER>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
|;
}
elsif ($action eq 'edit') {
    print header(-'content-encoding'=>'no',-charset=>'utf-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
function checkForm(myForm){
    if(!myForm.name.value){
        alert('Type a name for this GO analysis');
        return false;
    } else {
        return true;
    }
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER><FONT class="title1">GO Analysis <FONT color=red>$name</FONT></FONT><BR><BR>
<FORM name="editGO" method="post" onsubmit="return(checkForm(this))">
<INPUT type="hidden" name="ID" value=$goAnaID>
<TABLE border=0>
<TR><TD bgcolor=$dark>
<TABLE width=100% cellpadding=2 border=0>
<TR><TH align=right valign=top width=130>Name: </TH><TD bgcolor=$light><INPUT type="text" name="name" value="$name" size=50></TD></TR>
<TR><TH align=right valign=top width=130>Description: </TH><TD bgcolor=$light><TEXTAREA cols=50 name="des">$des</TEXTAREA></TD></TR>
<TR><TD colspan=2 align=center><INPUT type="submit" name="save" value="Save"> &nbsp <INPUT type="reset" value="Clear"> &nbsp <INPUT type="button" value="Cancel" onclick="top.promsFrame.optionFrame.location='$promsPath{cgi}/selectGOOption.cgi?ID=go_analysis:$goAnaID';"></TD></TR>
</TABLE></TD></TR></TABLE>
</FORM>
|;
}

sub processForm{
    my $name = param('name');
    my $des = param('des');
    my $sthUpGO = $dbh->prepare("UPDATE GO_ANALYSIS SET NAME=?, DES=?, UPDATE_DATE=NOW(), UPDATE_USER=? WHERE ID_GOANALYSIS=$goAnaID");
    $sthUpGO->execute($name,$des,$userInfo[0]);
    $sthUpGO->finish;
    $dbh->commit;
    $dbh->disconnect;
    print header;
    print qq
|<HTML><HEAD>
<SCRIPT LANGUAGE="JavaScript">
parent.itemFrame.location='$promsPath{cgi}/openProject.cgi?ACT=experiment&VIEW=go&branchID=go_analysis:$goAnaID';
</SCRIPT>
</HTML>
|;
}

sub deleteAna{
    &goAnalysis::deleteGOAnalysis($dbh, $projectID, $goAnaID);
    $dbh->commit;
    $dbh->disconnect;
    print header;
    print qq
|<HTML><HEAD>
<SCRIPT LANGUAGE="JavaScript">
    //parent.selectedAction='summary';
    //parent.navFrame.selectItem(parent.navFrame.selectedTab); // must be current experiment
    parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=experiment:$expID&branchID=go:$expID&VIEW=functAna";
 </SCRIPT>
</BODY></HTML>
|;
}

sub getWarnings{
    my ($projectID, $goAnaID) = @_;

    my $resultDirUnix="$promsPath{go_unix}/project_$projectID/$goAnaID";
    my $errorFileName = "$resultDirUnix/errorLog.txt";
    my $warnings = 0;

    open(ERRLOG,$errorFileName) || do { print "<FONT color=\"red\">Cannot open $errorFileName: $!</FONT>"; return };
        while(<ERRLOG>){
            chomp;
            $warnings = 1;
            print "$_<BR>\n";
        }
    close(ERRLOG);
    print 'No warnings' unless $warnings;
}

####>Revision history<####
# 1.1.7 [BUGFIX] Uses proper update user name in 'Last changed' field (PP 25/03/20)
# 1.1.6 [ENHANCEMENT] Made display compatible with multi-PTM GO quantifications (PP 26/02/20)
# 1.1.5 Explicit display of PTM if a PTM-quantif was used for a quanti-GO (PP 31/08/18)
# 1.1.4 Minor bug fixes and changes to handle parameters from PTM quantifications (PP 26/07/18)
# 1.1.3 Added number proteins used in Quanti-GO (PP 09/02/18)
# 1.1.2 Change so that "GO Analyses" branch is selected after deletion instead of Experiment (PP 26/10/15)
# 1.1.1 Uses &promsQuantif::getDistinctQuantifNames for description of quantif used (PP 18/12/14)
# 1.1.0 Change population Summary and add unannotated genes in Q.GO Analysis (GA 14/05/14)
# 1.0.9 Remove File::Path (PP 10/03/14)
# 1.0.8 -charset=>'utf-8' (PP 27/01/14)
# 1.0.7 Displaying ratio P-value threshold if exists and quantifications in Q. GO analysis case (03/04/13)
# 1.0.6 Displaying all error log content without filtering (FY 21/03/13)
# 1.0.5 Using goAnalysis::deleteGOAnalysis to delete current go analysis (FY 26/11/12)+<BR> Managing GO Quanti Analyses (FY 07/12/12)
# 1.0.4 Compatible with new data path .../results/project_ID/GoID (PP 13/09/12)
# 1.0.3 Add name and version of files in GO analysis summary<BR>& Minor code improvements (FY 23/04/12)
# 1.0.2 Change in itemFrame refresh parameters (PP 16/04/12)
# 1.0.1 Matching with openProject.cgi GO changes (FY 15/02/12)
# 1.0.0 New script to display GO analysis summary or edit it in resultFrame (FY 07/07/11)