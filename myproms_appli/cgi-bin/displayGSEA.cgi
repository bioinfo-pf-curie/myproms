#!/usr/local/bin/perl -w

################################################################################
# displayGSEA.cgi       1.1.1                                                  #
# Authors: V. Laigle (Institut Curie)                                          #
# Contact: myproms@curie.fr                                                    #
# Display the results of Gene Set Enrichment Analysis                          #
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
use POSIX qw(floor ceil);
use promsConfig;
use promsMod;
use Spreadsheet::WriteExcel;

#####################
### Configuration ###
#####################
my %promsPath = &promsConfig::getServerInfo;
my $dbh = &promsConfig::dbConnect;
my $userID = $ENV{'REMOTE_USER'};
my ($lightColor, $darkColor) = &promsConfig::getRowColors;
my ($red, $green, $black) = ("#AC0000", "#147E0D", "#000000");  # Colors for enrichment in test and reference conditions respectively

my $gseaID = &promsMod::cleanNumericalParameters(param('ID'));
my $action = (param('ACTION'))? param('ACTION') : 'display';
my $view = (param('VIEW'))? lc(param('VIEW')) : 'graph';
my $sortKey = (param('SORT'))? param('SORT') : 'padj';
my $searchText = param("SEARCH") || '';
my $searchIdType = param("SEARCH_ID_TYPE") || '';
my $projectID = &promsMod::getProjectID($dbh, $gseaID,'PATHWAY_ANALYSIS');
my @userInfo = &promsMod::getUserInfo($dbh, $userID, $projectID);
my $projectAccess = ${$userInfo[2]}{$projectID};
my $disabSave = ($projectAccess eq 'guest')? ' disabled' : '';

my $gseaDir = "$promsPath{'gsea'}/project_$projectID/gsea_$gseaID";
my $dataFile = "$gseaDir/gsea_data.txt";
my $resultsFile = "$gseaDir/gsea_results.txt";
my $genesUsedFile = "$gseaDir/gsea_genes_used.txt";
my (%gseaInfo, %gseaResults, %geneIsPresent, %gseaData, %protResults);
my (%terms2Genes, %genes2Terms);

my ($expID, $name, $desc, $param, $status, $recDate, $upDate, $user) = $dbh->selectrow_array(
    "SELECT ID_EXPERIMENT, NAME, DES, PARAM_STRG, STATUS, RECORD_DATE, UPDATE_DATE, UPDATE_USER
    FROM PATHWAY_ANALYSIS WHERE ID_PATHWAY_ANALYSIS = $gseaID"
);
$desc='' unless $desc;

&getGSEAInfo(\%gseaInfo);
&parseGSEAResults(\%gseaResults, \%geneIsPresent, \%protResults, $resultsFile, $genesUsedFile);
&parseGSEAData(\%gseaData, $dataFile);
&linkTermsAndGenes(\%terms2Genes, \%genes2Terms);
$dbh->disconnect;

#####################
### Start Display ###
#####################
if ($action eq 'display') {
    &displayVisualization;
    exit;
} elsif ($action eq 'details') {  # AJAX
    my $geneSet = param('geneSet');
    &printGeneSetDetails($geneSet, \%gseaInfo, \%gseaResults);
    exit;
} elsif ($action eq 'graph') {  # AJAX
    my $geneSet = param('geneSet') || '';
    &showGraph($geneSet);
    exit;
} elsif ($action eq 'listProts') {  # AJAX
    my $protIDsArray = (param('protID')) ? [param('protID')] : '';
    &displayProtGeneSets($protIDsArray);
    exit;
} elsif ($action eq 'export') {
    if (param('XLS') eq 'results') {
        &exportGSEAResults;
        exit;
    } elsif (param('XLS') eq 'proteins'){
        my $protIDsArray = (param('protID')) ? [param('protID')] : '';
        &exportProtList($protIDsArray);
        exit;
    }
}


sub displayVisualization {  # Globals: %promsPath, $gseaID, $projectID, $name, $view

    print header(-'content-encoding' => 'no',-charset=>'utf-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
    <TITLE>Display Gene Set Enrichment Analysis</TITLE>
    <LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
    <STYLE type="text/css">
        .popup {
            z-index:999;
            background-color:#FFFFFF;
            border:solid 3px #999999;
            padding:5px;
            position:absolute;
            display:none;
            box-shadow:10px 10px 20px #808080;
        }
    </STYLE>
    <SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
    <SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
    <SCRIPT src="$promsPath{html}/js/local/genericPlot.js"></SCRIPT> <!-- Needed ? -->
    <SCRIPT src="$promsPath{html}/js/local/vennDiagram.js"></SCRIPT> <!-- Needed ? -->
    <SCRIPT src="$promsPath{html}/js/local/peptidePlot.js"></SCRIPT> <!-- Needed ? -->
    <SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
    <SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
    <SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
    <SCRIPT src="$promsPath{html}/js/idvis/idvis.js"></SCRIPT>
    <SCRIPT src="$promsPath{html}/js/idvis/idvis.categoryPlot.js"></SCRIPT>
    <SCRIPT src="$promsPath{html}/js/idvis/idvis.heatMap.js"></SCRIPT>
    <SCRIPT LANGUAGE="JavaScript">
// Dictionnary of help popup texts (avoid using '"')
const helpText = {
    rank: 'The number of genes/proteins until we reach the Enrichment Score.',
    tags: 'The percentage of genes contributing to the enrichment score.',
    list: 'At which point is the enrichment score attained in the list.',
    signal: 'The enrichment signal strength.',
    genesUsed: 'Genes in bold are the genes/proteins found in this set that are present in the quantification and used during the GSEA.',
    leadingEdge: 'These are the genes/proteins that contribute the most to the enrichment of the current gene set in the specified condition. Specifically, they are the proteins which score contribute to the maximum Enrichment Score. They are sorted from the most enriched to the least enriched in the quantification.'
};
|;
    &promsMod::popupInfo();
    print qq
|function getElementPosition(e) {
    var left=0;
    var top=0;
    while (e.offsetParent != undefined && e.offsetParent != null) { // Adding parent item position
        left += e.offsetLeft + (e.clientLeft != null ? e.clientLeft : 0);
        top += e.offsetTop + (e.clientTop != null ? e.clientTop : 0);
        e = e.offsetParent;
    }
    return [top,left];
}

function getXMLHTTP(){
    var xhr=null;
    if(window.XMLHttpRequest) {  // Firefox & others
        xhr = new XMLHttpRequest();
    } else if(window.ActiveXObject){  // Internet Explorer
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

var XHR=null;
function ajaxGSDetails(geneSet, searchedProtID) {
    saveListGlobals.themesFetched=false; // used for themes management, need to reset at each call
    var geneSetDiv = document.getElementById('divGeneSets');
    geneSetDiv.innerHTML = '<IMG src="$promsPath{images}/scrollbarGreen.gif">';
    geneSetDiv.style.display = 'block';

    // If XHR object already exists, the request is canceled & the object is deleted
    if(XHR && XHR.readyState != 0){
        XHR.abort();
        delete XHR;
    }

    // Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if (!XHR) {
        return false;
    }
    var GETstring = 'ID=' + $gseaID + '&ACTION=details&geneSet=' + geneSet + '&SEARCH_PROT=' + searchedProtID;
    var URL = "$promsPath{cgi}/displayGSEA.cgi?" + GETstring;
    XHR.open("GET", URL, true);
    XHR.onreadystatechange = function() {
        if (XHR.readyState == 4 && XHR.responseText){
            geneSetDiv.innerHTML = XHR.responseText;
            geneSetDiv.scrollIntoView();
        }
    }
    XHR.send(null);
}

function showGraph(button, geneSet) {
    var displayDiv = document.getElementById('displayDIV');
    var infoDiv = document.getElementById('infoDIV');

    var [top,left] = getElementPosition(button);
    top += 30;
    left -= 230;
    if (left < 0) {left = 0;}
    displayDiv.style.left = left + 'px';
    displayDiv.style.top = top + 'px';
    displayDiv.style.display = 'block';

    // If XHR object already exists, the request is canceled & the object is deleted
    if(XHR && XHR.readyState != 0){
        XHR.abort();
        delete XHR;
    }

    // Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if (!XHR) {
        return false;
    }

    var GETstring = 'ID=' + $gseaID + '&ACTION=graph';
    if (geneSet) {
        GETstring += '&geneSet=' + geneSet;
    }
    var URL = "$promsPath{cgi}/displayGSEA.cgi?" + GETstring;
    XHR.open("GET", URL, true);
    XHR.onreadystatechange = function() {
        if (XHR.readyState == 4 && XHR.responseText){
            infoDiv.innerHTML = XHR.responseText;
        }
    }
    XHR.send(null);
}
|;
    &promsMod::printAjaxManageSaveProteins($projectID, \%promsPath, 'document.protForm.chkProt');
    print qq
|function sequenceView(id_protein, anaIdStrg) {
    var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+anaIdStrg+"&id_prot="+id_protein+"&msdata="+top.showMSdata;
    top.openProtWindow(winLocation);
}
function sortTable(sortItem) {
    window.location="./displayGSEA.cgi?ID=$gseaID&ACTION=$action&VIEW=$view&SORT="+sortItem+"&SEARCH=$searchText&SEARCH_ID_TYPE=$searchIdType";
}
function checkAllProteins(chkStatus) {
    var checkBoxList = document.protForm.chkProt;
    if (checkBoxList.length) {
        for (var i=0; i < checkBoxList.length; i++) {
            checkBoxList[i].checked = chkStatus;
        }
    } else {
        checkBoxList.checked = chkStatus;
    }
}
function checkLEdgeProts(chkStatus, protList) {
    var checkBoxList = protList.map(
        function(protID) {
            return document.getElementById('chk_' + protID);
        }
    );
    if (checkBoxList.length) {
        for (var i=0; i < checkBoxList.length; i++) {
            checkBoxList[i].checked = chkStatus;
        }
    } else {
        checkBoxList.checked = chkStatus;
    }
}
function filterGeneSets(filter) {
    var searchText = document.getElementById('SEARCH').value;
    var searchIdType = document.getElementById('SEARCH_ID_TYPE').value;
    if (filter) {
        window.location="$promsPath{cgi}/displayGSEA.cgi?ID=$gseaID&ACTION=$action&VIEW=$view&SEARCH="+searchText+"&SEARCH_ID_TYPE="+searchIdType;
    } else {
        window.location="$promsPath{cgi}/displayGSEA.cgi?ID=$gseaID&ACTION=$action&VIEW=$view";
    }

}
   </SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<DIV id="displayDIV" class="popup">
    <DIV id="infoDIV"></DIV>
</DIV>
<CENTER>
<FONT class="title">Display Gene Set Enrichment Analysis for <FONT color="red">$name</FONT></FONT>
&nbsp;&nbsp;
<INPUT type="button" class="title3" value="Export" onclick="window.location='$promsPath{cgi}/displayGSEA.cgi?ACTION=export&ID=$gseaID&XLS=results';">
<BR><BR>
Test condition: <FONT color="$red" style="font-weight:bold">$gseaInfo{'testCondName'}</FONT>
&nbsp;&nbsp;&nbsp;&nbsp;
Reference condition: <FONT color="$green" style="font-weight:bold">$gseaInfo{'refCondName'}</FONT>
&nbsp;&nbsp;&nbsp;&nbsp;
<LABEL for='dataVizType' style='font-weight:bold'>Visualization type: </LABEL>
<SELECT id='dataVizType' onchange="window.location='$promsPath{cgi}/displayGSEA.cgi?ACTION=display&VIEW=' + this.value + '&ID=$gseaID';">
|;

    my @vizTypes = ('Graph', 'HeatMap', 'Table');
    foreach my $vizType (@vizTypes) {
        my $disabStrg = (lc($vizType) eq $view)? " selected disabled" : " ";
        print "<OPTION value=\"$vizType\"$disabStrg>$vizType</OPTION>\n";
    }
    print qq
|</SELECT>
<BR><BR>
|;

    if ($view eq 'graph') {
        &displayDotplotGSEA;
    } elsif ($view eq 'heatmap') {
        &displayHeatmapGSEA;
    } else {
        &displayTableGSEA;
        print "<BR><BR><BR>";
    }

    print qq
|<BR><BR><BR>
<DIV id="divGeneSets" style="float:center;display:none">
<!--Empty div-->
</DIV>
</CENTER>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup();</SCRIPT>
</BODY>
</HTML>
|;

}

sub displayTableGSEA {  # Globals: %gseaResults, $lightColor, $darkColor, $sortKey, $searchText, $searchIdType

    my @filteredGeneSets;
    my $searchedProtID = "null";  # To avoid printing an empty string. It is passed to JS
    if ($searchText) {
        ($searchedProtID, my $searchedGene) = &searchConvertIdentifier();
        if ($searchedGene) {
            @filteredGeneSets = @{$genes2Terms{$searchedGene}};
        } else {
            @filteredGeneSets = ();
        }
    }

    my %allIdentifiers = &getAllIdentifierTypes();
    my $identSelStrg;
    foreach my $identID (sort {$a <=> $b} keys %allIdentifiers) {
        my $selectedStrg = ($searchIdType && $identID == $searchIdType)? " selected" : " ";
        $identSelStrg .= "<OPTION value=\"$identID\"$selectedStrg>$allIdentifiers{$identID}</OPTION>\n";
    }

    my ($gsTermColor, $gsNbColor, $esColor, $nesColor, $pAdjColor, $qValColor) = 
        ($sortKey eq "setSize")? ("black", "red", "black", "black", "black", "black") :
        ($sortKey eq "ES")?      ("black", "black", "red", "black", "black", "black") :
        ($sortKey eq "NES")?     ("black", "black", "black", "red", "black", "black") :
        ($sortKey eq "padj")?    ("black", "black", "black", "black", "red", "black") :
        ($sortKey eq "qval")?    ("black", "black", "black", "black", "black", "red") :
        ("red", "black", "black", "black", "black", "black");  # $sortKey eq "gsTerm"

    sub sortTable{
        if ($sortKey eq 'padj') {
            ($gseaResults{$a}{'padj'} <=> $gseaResults{$b}{'padj'}) ||
            (abs($gseaResults{$b}{'NES'}) <=> abs($gseaResults{$a}{'NES'}));
        } elsif ($sortKey eq 'setSize') {
            ($gseaResults{$b}{'setSize'} <=> $gseaResults{$a}{'setSize'}) ||
            ($gseaResults{$a}{'padj'} <=> $gseaResults{$b}{'padj'});
        } elsif ($sortKey eq 'ES') {
            ($gseaResults{$b}{'ES'} <=> $gseaResults{$a}{'ES'}) ||
            ($gseaResults{$b}{'NES'} <=> $gseaResults{$a}{'NES'}) ||
            ($gseaResults{$a}{'padj'} <=> $gseaResults{$b}{'padj'});
        } elsif ($sortKey eq 'NES') {
            ($gseaResults{$b}{'NES'} <=> $gseaResults{$a}{'NES'}) ||
            ($gseaResults{$b}{'ES'} <=> $gseaResults{$a}{'ES'}) ||
            ($gseaResults{$a}{'padj'} <=> $gseaResults{$b}{'padj'});
        } elsif ($sortKey eq 'qval') {
            ($gseaResults{$a}{'qval'} <=> $gseaResults{$b}{'qval'}) ||
            (abs($gseaResults{$b}{'NES'}) <=> abs($gseaResults{$a}{'NES'}));
        } else {  # $sortKey eq "gsTerm"
            ($a cmp $b);
        }
    }

    print qq
|<BR>
<DIV id="searchDIV">
    Filter gene sets related to a specific protein or gene:
    &nbsp;&nbsp;
    <INPUT type="text" id="SEARCH" name="SEARCH" maxlength="30" placeholder="MYC" value="$searchText">
    &nbsp;&nbsp;
    Search by:
    &nbsp;&nbsp;
    <SELECT id="SEARCH_ID_TYPE" name="SEARCH_ID_TYPE">
        $identSelStrg
    </SELECT>
    &nbsp;&nbsp;
    <INPUT type="button" id="SRCH_BUTTON" name="SRCH_BUTTON" value="Filter" onclick="filterGeneSets(true);">
    &nbsp;&nbsp;
    <INPUT type="button" id="FULL_BUTTON" name="FULL_BUTTON" value="Full Table" onclick="filterGeneSets(false);">
</DIV>
<BR>
<TABLE id="resultsTable" border=0 cellspacing=0 align=center>
<TR bgcolor=$darkColor>
    <TH class="rbBorder" rowspan="2" nowrap>
        <A href="javascript:sortTable('gsTerm')"><FONT color="$gsTermColor">
            &nbsp;&nbsp;Gene Set Term&nbsp;&nbsp;
        </A>
    </TH>
    <TH class="rbBorder" rowspan="2" nowrap>
        <A href="javascript:sortTable('setSize')"><FONT color="$gsNbColor">
            &nbsp;&nbsp;Nb of Genes&nbsp;&nbsp;<br>&nbsp;&nbsp;found in Set&nbsp;&nbsp;
        </A>
    </TH>
    <TH class="rbBorder" colspan="2" nowrap>&nbsp;&nbsp;Enrichment Score&nbsp;&nbsp;</TH>
    <TH class="rbBorder" rowspan="2" nowrap>
        <A href="javascript:sortTable('padj')"><FONT color="$pAdjColor">
            &nbsp;&nbsp;Adjusted&nbsp;&nbsp;<br>&nbsp;&nbsp;p-value&nbsp;&nbsp;
        </A>
    </TH>
    <TH class="rbBorder" rowspan="2" nowrap>
        <A href="javascript:sortTable('qval')"><FONT color="$qValColor">
            &nbsp;&nbsp;Q-value&nbsp;&nbsp;
        </A>
    </TH>
    <TH class="rbBorder" colspan="2" nowrap>&nbsp;&nbsp;Display&nbsp;&nbsp;</TH>
</TR>
<TR bgcolor=$darkColor>
    <TH class="rbBorder" nowrap>
        <A href="javascript:sortTable('ES')"><FONT color="$esColor">
            &nbsp;&nbsp;&nbsp;&nbsp;Raw&nbsp;&nbsp;&nbsp;&nbsp;
        </A>
    </TH>
    <TH class="rbBorder" nowrap>
        <A href="javascript:sortTable('NES')"><FONT color="$nesColor">
            &nbsp;&nbsp;Normalized&nbsp;&nbsp;
        </A>
    </TH>
    <TH class="rbBorder" nowrap>&nbsp;&nbsp;&nbsp;&nbsp;Gene Set&nbsp;&nbsp;&nbsp;&nbsp;</TH>
    <TH class="rbBorder" nowrap>&nbsp;&nbsp;Enrichment graph&nbsp;&nbsp;</TH>
</TR>
|;

    my $lineColor = $lightColor;
    my $geneSetCount = 0;
    foreach my $geneSet (sort { sortTable } keys %gseaResults) {
        next if ($searchText && !grep(/$geneSet/, @filteredGeneSets));
        $geneSetCount++;
        (my $geneSetStrg = $geneSet) =~ s/_/ /g;
        $geneSetStrg =~ s/([\w']+)/\u\L$1/g;  # capitalize first letter of each word
        my $graphPath = "$promsPath{'gsea'}/project_$projectID/gsea_$gseaID/gsea_graph_$geneSet.png";
        my $graphStrg = '';
        if (-e $graphPath) {
            $graphStrg = "<INPUT type=\"button\" class=\"font11\" id=\"$geneSet\" name=\"$geneSet\" value=\"Show Graph\" onclick=\"showGraph(this, '$geneSet')\">";
        }
        my $detailsStrg = "<INPUT type=\"button\" class=\"font11\" id=\"$geneSet\" name=\"$geneSet\" value=\"Details\" onclick=\"ajaxGSDetails('$geneSet', $searchedProtID)\">";
        my $gsEnrichment = ($gseaResults{$geneSet}{'NES'} > 0)? "<IMG src=\"$promsPath{images}/up_red.png\"/>" : "<IMG src=\"$promsPath{images}/down_green.png\"/>";
        print qq
|<TR bgcolor="$lineColor" valign="top">
    <TD class="title3 rBorder" align="left" id="td_$geneSet" nowrap>
        &nbsp;&nbsp;$gsEnrichment&nbsp;&nbsp;
        <A href="javascript:ajaxGSDetails('$geneSet', $searchedProtID)">$geneSetStrg</A>
        &nbsp;&nbsp;
    </TD>
    <TD class="rBorder" align="center" nowrap>&nbsp;&nbsp;$gseaResults{$geneSet}{'setSize'}&nbsp;&nbsp;</TD>
    <TD class="rBorder" align="center" nowrap>&nbsp;&nbsp;${\sprintf("%.3g", $gseaResults{$geneSet}{'ES'})}&nbsp;&nbsp;</TD>
    <TD class="rBorder" align="center" nowrap>&nbsp;&nbsp;${\sprintf("%.3g", $gseaResults{$geneSet}{'NES'})}&nbsp;&nbsp;</TD>
    <TD class="rBorder" align="center" nowrap>&nbsp;&nbsp;${\sprintf("%.3e", $gseaResults{$geneSet}{'padj'})}&nbsp;&nbsp;</TD>
    <TD class="rBorder" align="center" nowrap>&nbsp;&nbsp;${\sprintf("%.3e", $gseaResults{$geneSet}{'qval'})}&nbsp;&nbsp;</TD>
    <TD class="rBorder" align="center" nowrap>&nbsp;&nbsp;$detailsStrg&nbsp;&nbsp;</TD>
    <TD class="rBorder" align="center" nowrap>&nbsp;&nbsp;$graphStrg&nbsp;&nbsp;</TD>
</TR>
|;
        $lineColor = ($lineColor eq $lightColor)? $darkColor : $lightColor;
    }
    print "</TABLE>\n";
    if ($searchText && $geneSetCount == 0) {
        print "</CENTER><BR>These filters did not match any gene set.<BR>Check that the identifier type corresponds to your search or that you did not mispelled your protein of interest.<BR>Otherwise try searching for another protein or reset the full table.<CENTER>\n";
    } else {
        print "<BR>End of Gene Sets.\n";
    }
}


sub displayHeatmapGSEA {
    print qq
|<BR>
<DIV id="heatmapDIV" style="float:center;padding:2px;"></DIV>  <!-- To display the actual heatmap -->
<DIV id="protTable" style="display:none"></DIV>  <!-- For lists of proteins -->
<BR>
|;
    my ($jsColumnsStr, $jsAddRowStr, $jsDendroStr, $protDataStrg, $dataStrg, $min, $max) = computeHeatMap();
    printHeatMap($jsColumnsStr, $jsAddRowStr, $jsDendroStr, $protDataStrg, $dataStrg, $min, $max);
}


sub computeHeatMap {  # Globals: %gseaResults, %geneIsPresent, %gseaData, %protResults, $gseaDir
    my ($jsColumnsStr, $jsAddRowStr, $jsDendroStr, $protDataStrg, $dataStrg, $min, $max) = ('', '', '', '', '', 10000000000, -10000000000);
    my ($geneSetsNb, $protNb);
    my (%allCoreGenes, %genesToSets);
    my (%protInfo, %masterProts);
    my %protValues;

    foreach my $geneSet (sort {$gseaResults{$a}{'NES'} <=> $gseaResults{$b}{'NES'}} keys %gseaResults) {
        my $setSize = $gseaResults{$geneSet}{'setSize'};
        my $NES     = sprintf("%.2f", $gseaResults{$geneSet}{'NES'});
        my $qvalue  = sprintf("%.2e", $gseaResults{$geneSet}{'qval'});
        my $padj    = sprintf("%.2e", $gseaResults{$geneSet}{'padj'});

        $jsColumnsStr .= ',' if ($jsColumnsStr);
        $jsColumnsStr .= "[\"$geneSet\", \"$geneSet\", \"NES: $NES\\nAdj. p-value: $padj\\nq-value: $qvalue\\nSet size: $setSize\"]";

        foreach my $gene (@{$gseaResults{$geneSet}{'coreGenes'}}) {
            $allCoreGenes{$gene} = $geneIsPresent{$gene} unless ($allCoreGenes{$gene});
            $genesToSets{$gene}{$geneSet} = 1;
        }
    }
    $jsColumnsStr = "HM_GSEA.setColumns([$jsColumnsStr]);";
    
    $geneSetsNb = scalar keys %gseaResults;
    $protNb     = scalar keys %allCoreGenes;

    my @allCoreProtIDs = values %allCoreGenes;
    my ($protInfoRef, $masterProtsRef) = &getProtInfo(\@allCoreProtIDs);
    %protInfo = %{$protInfoRef};
    %masterProts = %{$masterProtsRef};

    # Build prots Dendogram (if it does not exist yet)
    my $hmMatrixFile = "matrix_no_cluster.tsv";
    my $heatmapDir = "$gseaDir/Heatmap";
    mkdir $heatmapDir unless (-d $heatmapDir);
    # if (!-e "$heatmapDir/$hmMatrixFile" && $geneSetsNb > 0 && $protNb > 0) {
        open(MATRIX, ">$heatmapDir/$hmMatrixFile") || die "Could not open matrix file";        
        print MATRIX "Gene Sets";
        foreach my $geneSet (sort {$gseaResults{$a}{'NES'} <=> $gseaResults{$b}{'NES'}} keys %gseaResults) {
            print MATRIX "\t$geneSet";
        }
        print MATRIX "\n";

        foreach my $gene (keys %allCoreGenes) {
            # Double check that data corresponds to gene/prot
            next unless ($gene eq $gseaData{$allCoreGenes{$gene}}{'geneSetId'});
            my $protID = $allCoreGenes{$gene};
            my $valueStrg = $gseaData{$protID}{'value'};
            my $value = $valueStrg;
            if ($value =~ /1\/(.*)/) {
                my $denom = $1;
                $value = ($denom eq 'Inf')? 0.001 : 1 / $denom;
            } elsif ($value =~ /Inf/) {
                $value = 1000;
            }
            $value = sprintf("%.2f", log($value) / log(2));
            
            # Value used for backcompatibility when score was not recorded. Shoud not actually be used
            my $score = ($protResults{$protID}{'score'})? $protResults{$protID}{'score'} : ($value)? $value : 0.00;
            $max = $score if ($score > $max);
            $min = $score if ($score < $min);

            my @valuesLine;
            foreach my $geneSet (sort {$gseaResults{$a}{'NES'} <=> $gseaResults{$b}{'NES'}} keys %gseaResults) {
                my $geneExistsInSet = (defined $genesToSets{$gene}{$geneSet})? $score : 0.00;
                push @valuesLine, $geneExistsInSet;
            }
            $protValues{$protID} = \@valuesLine;
            
            print MATRIX "$protID\t" . join("\t", @valuesLine) . "\n";
        }
        close MATRIX;

        # Reuse exiting script for GO
        system "cd $heatmapDir; $promsPath{R}/R CMD BATCH --no-save --no-restore '--args $hmMatrixFile $geneSetsNb $protNb' $promsPath{R_scripts}/goHierarchCluster.R";
    # }
    
    # Load dendogram from results file (names are go... because we are reusing a script from GO analysis)
    my @protOrder;
    if(-s "$heatmapDir/goDendro_$geneSetsNb\_$protNb.txt" && -s "$heatmapDir/goOrder_$geneSetsNb\_$protNb.txt") {
        my @dendroElements;
        my %orderedProtIds;
        open(ORDER, "$heatmapDir/goOrder_$geneSetsNb\_$protNb.txt") or die $!;
        while (my $lineOrder = <ORDER>){
            next if ($. == 1);
            chomp($lineOrder);
            my ($index, $pos, $protID) = split("\t", $lineOrder);
            push @protOrder, $protID;
            $orderedProtIds{$pos} = $protID;
        }
        close ORDER;
        
        if (@protOrder) {
            my $maxIndex = scalar @protOrder;
            open(DENDRO, "$heatmapDir/goDendro_$geneSetsNb\_$protNb.txt") or die $!;
            while (my $line = <DENDRO>) {
                next if ($. == 1);
                chomp($line);
                my @line = split(/\t/, $line);
                $line[1] =~ s/-(\d+)/-$orderedProtIds{$1}/ if ($line[1] =~ /-\d+/);
                $line[2] =~ s/-(\d+)/-$orderedProtIds{$1}/ if ($line[2] =~ /-\d+/);
                push @dendroElements, join(',', @line);
            }
            close DENDRO;
            $jsDendroStr = "HM_GSEA.addTree('row', '" . join(';', @dendroElements) . "');";
        }
    }
    
    # Build protein rows
    foreach my $protID (@protOrder) {
        # Fill heatmap row cells
        my $protValuesStr = '';
        my $protName = ($protInfo{$protID}{'Name'}) ? $protInfo{$protID}{'Name'} : $protID;
        my $protDes = ($protInfo{$protID}{'Des'}) ? $protInfo{$protID}{'Des'} : '';
        $protDes =~ s/"/'/;
        $jsAddRowStr .= "HM_GSEA.addRow([\"$protName\", \"$protID\", \"$protDes\"], [";
        my $geneSetIdx = 0;
        foreach my $geneSet (sort {$gseaResults{$a}{'NES'} <=> $gseaResults{$b}{'NES'}} keys %gseaResults) {    
            $protValuesStr .= ',' if ($geneSetIdx > 0);
            $protValuesStr .= ($protValues{$protID}->[$geneSetIdx])? $protValues{$protID}->[$geneSetIdx] : 'null';
            $geneSetIdx++;
        }
        $jsAddRowStr .= "$protValuesStr]);\n";

        # Fill heatmap cells popups
        my $quantifValue = $gseaData{$protID}{'value'};
        my $quantifPval  = $gseaData{$protID}{'pvalue'};
        my $pepUsed      = $gseaData{$protID}{'pepNb'};
        $protDataStrg .= ',' if ($protDataStrg);
        $protDataStrg .= "\"$protID\": ['$quantifValue', '$quantifPval', '$pepUsed']";
        
        # Fill data on prot clicked
        $dataStrg .= ',' if ($dataStrg);
        my $anaID = ($protInfo{$protID}{'AnaID'}) ? $protInfo{$protID}{'AnaID'} : -1;
        $dataStrg .= "\"$protID\":$anaID";
    }
    $protDataStrg = "{$protDataStrg}";
    $dataStrg = "{$dataStrg}";
    return ($jsColumnsStr, $jsAddRowStr, $jsDendroStr, $protDataStrg, $dataStrg, $min, $max);
}


sub printHeatMap {
    my ($jsColumnsStr, $jsAddRowStr, $jsDendroStr, $protDataStrg, $dataStrg, $min, $max) = @_;
    $protDataStrg = "{}" if (!$protDataStrg);
    $dataStrg     = "{}" if (!$dataStrg);
    $min = floor($min);
    $max = ceil($max);
    my $limits = ($min < 0)? "min: $min, ref: 0, max: $max" : "min: 0, max: $max";
    my $pepStrg = ($gseaInfo{'pepWeights'} eq "msms")? "MS/MS peptides" : ($gseaInfo{'pepWeights'} eq "distinct")? "Dist. peptides" : "All peptides";

    # Draw HeatMap
    print qq 
|<SCRIPT>
var data = $dataStrg;
var protData = $protDataStrg;
var HM_GSEA;
var rowList_GSEA;
var columnList_GSEA;

function drawHeatMap() {
    HM_GSEA = new idvis.heatMap({
        div: 'heatmapDIV',
        exportAsImage:['Export as image', 'GSEA_heatmap', '$promsPath{cgi}/exportSVG.cgi'],
        cellHeight: 8,
        moveLabel: true,
        editable: true,
        editableItems: {scope :false, type: false, color: true},
        entities: {row: 'Proteins', column: 'Gene Sets'},
        normalization: {scope: 'row', reference: 'user', colors: 'BbR', limitValues: {$limits}},
        showTrace: false,
        labelsOnSelect: {row: {text: 'List proteins', action: listProteins, ranking: true}},
        cellOnMouseOver: showValue,
        singleCellSelection: true,
        cellOnClick: displayProtGeneSets,
        columnOnClick: displayGeneSet,
        rowOnClick: displayProtInfo
    });

$jsColumnsStr
$jsAddRowStr
$jsDendroStr

    rowList_GSEA = HM_GSEA.getRowList();
    columnList_GSEA = HM_GSEA.getColumnList();
    
    HM_GSEA.draw();
}

function showValue(value, hCell) {
    var valueStrg = '';
    if (value == null) {
        valueStrg = 'Not Leading Edge or no Score';
    } else {
        valueStrg = 'Protein Score: ' + value.toLocaleString('en-US', {maximumFractionDigits: 2});

        if (hCell) {
            var rowIdData = rowList_GSEA\[hCell.rowIdx].id.split(':')
            var colId = columnList_GSEA\[hCell.columnIdx].id;

            if (protData[rowIdData[0]]) {
                valueStrg += '\\nQuantif. Value: ' + protData[rowIdData[0]][0];
                valueStrg += '\\nQuantif. p-value: ' + protData[rowIdData[0]][1];
                valueStrg += '\\n$pepStrg: ' + protData[rowIdData[0]][2];
            }
        }
    }
    return valueStrg;
}

function displayProtGeneSets(rowId, colId, status) {
    if (!status \|\| !rowId \|\| rowId == '') {
        return;
    }
    var geneSetsDiv = document.getElementById("divGeneSets");
    geneSetsDiv.style.display = 'none';

    var protTable = document.getElementById('protTable');
    protTable.style.display = '';
    protTable.innerHTML = '<BR><BR><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">';

    // Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if (!XHR) {
        return false;
    }

    var URL = "$promsPath{cgi}/displayGSEA.cgi?ID=$gseaID&ACTION=listProts&protID=" + rowId;
    XHR.open("GET", URL, true);
    XHR.onreadystatechange = function() {
        if (XHR.readyState == 4 && XHR.responseText){
            protTable.innerHTML = XHR.responseText;
            protTable.scrollIntoView();
        }
    }
    XHR.send(null);
}

function displayGeneSet(colId) {
    if (!colId) return;
    var protTable = document.getElementById("protTable");
    protTable.style.display = 'none';

    ajaxGSDetails(colId);
}

function displayProtInfo(rowId) {
    if (!rowId) return;
    if (data === undefined \|\| data[rowId] === undefined) return;
    sequenceView(rowId, data[rowId]);
}

function listProteins(selectedProteins) {
    var geneSetsDiv = document.getElementById("divGeneSets");
    geneSetsDiv.style.display = 'none';

    var protTable = document.getElementById("protTable");
    protTable.style.display = '';
    protTable.innerHTML = '<BR><BR><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">';

    // Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if (!XHR) {
        return false;
    }

    var URL = "$promsPath{cgi}/displayGSEA.cgi?ID=$gseaID&ACTION=listProts&protID=" + selectedProteins.join("&protID=");
    XHR.open("GET", URL, true);
    XHR.onreadystatechange = function() {
        if (XHR.readyState == 4 && XHR.responseText){
            protTable.innerHTML = XHR.responseText;
            protTable.scrollIntoView();
        }
    }
    XHR.send(null);
}

drawHeatMap()
</SCRIPT>
|;
}


sub getProtInfo {  # Globals: $gseaID
    my ($protListRef) = @_;
    my (%protInfo, %masterProteins);
    my $dbh = &promsConfig::dbConnect;

    my ($uniAccID) = $dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='AC'");
    my ($geneNameID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
    
    my $sthProt = $dbh->prepare("SELECT ALIAS, ID_MASTER_PROTEIN, PROT_DES, ORGANISM FROM PROTEIN WHERE ID_PROTEIN = ?");
    my $sthUniprotACC = $dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN = ? AND ID_IDENTIFIER = $uniAccID ORDER BY IDENT_RANK");
    my $sthGeneNames = $dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN = ? AND ID_IDENTIFIER = $geneNameID ORDER BY IDENT_RANK");
    
    my $sthAna = $dbh->prepare("SELECT AP.ID_ANALYSIS, AP.NUM_PEP
        FROM PATHWAYANA_QUANTIFICATION PAQ
        INNER JOIN ANA_QUANTIFICATION AQ ON PAQ.ID_QUANTIFICATION = AQ.ID_QUANTIFICATION
        INNER JOIN ANALYSIS_PROTEIN AP ON AQ.ID_ANALYSIS = AP.ID_ANALYSIS
        WHERE PAQ.ID_PATHWAY_ANALYSIS = $gseaID AND AP.ID_PROTEIN = ? AND AP.VISIBILITY > 0
        ORDER BY NUM_PEP DESC LIMIT 0,1"
    );

    foreach my $protID (@{$protListRef}) {
        $sthProt->execute($protID);
        my ($protName, $masterProtID, $des, $org) = $sthProt->fetchrow_array;
        $sthAna->execute($protID);
        my ($anaID, $numPep) = $sthAna->fetchrow_array;

        $protInfo{$protID}{'Name'}       = $protName;
        $protInfo{$protID}{'Des'}        = $des;
        $protInfo{$protID}{'Org'}        = $org;
        $protInfo{$protID}{'AnaID'}      = $anaID;
        $protInfo{$protID}{'NumPep'}     = $numPep;
        $protInfo{$protID}{'MasterProt'} = $masterProtID if ($masterProtID);

        if($masterProtID) {
            # Retrieve UniProt ACC
            if(!defined $masterProteins{$masterProtID}{'uniACC'}) {
                @{$masterProteins{$masterProtID}{'uniACC'}} = ();
                $sthUniprotACC->execute($masterProtID);
                while (my ($uniprotACC) = $sthUniprotACC->fetchrow_array) {
                    push @{$masterProteins{$masterProtID}{'uniACC'}}, $uniprotACC;
                }
            }
            # Retrieve Gene Name
            if(!defined $masterProteins{$masterProtID}{'gene'}) {
                @{$masterProteins{$masterProtID}{'gene'}} = ();
                $sthGeneNames->execute($masterProtID);
                while (my ($gene) = $sthGeneNames->fetchrow_array) {
                    push @{$masterProteins{$masterProtID}{'gene'}}, $gene;
                }
            }
        }
    }
    $sthProt->finish;
    $sthAna->finish;
    $sthUniprotACC->finish;
    $sthGeneNames->finish;
    $dbh->disconnect;

    return (\%protInfo, \%masterProteins);
}


sub displayDotplotGSEA {  # Globals: %promsPath, %gseaResults, %gseaInfo, $project_ID, $gseaID
    my $geneSetsNb = scalar keys %gseaResults;
    my $width = 350;
    my $height = 45 * $geneSetsNb + 20;
    my $maxGeneSetSize = 0;
    my $maxQvalue = 0;
    foreach my $geneSet (keys %gseaResults) {
        $maxGeneSetSize = ($gseaResults{$geneSet}{'setSize'} > $maxGeneSetSize)? $gseaResults{$geneSet}{'setSize'} : $maxGeneSetSize;
        $maxQvalue = ($gseaResults{$geneSet}{'qval'} > $maxQvalue)? $gseaResults{$geneSet}{'qval'} : $maxQvalue;
    }
    $maxGeneSetSize = sprintf("%.0f", $maxGeneSetSize / 40) * 40;  # Round to nearest multiple of 40 (LCM of 10 and 8)
    $maxQvalue = sprintf("%.2f", $maxQvalue) + 0.01;  # Nearest hundredth up

    my $catListStrg = '';
    my $gsCount = 0;
    foreach my $geneSet (sort {$gseaResults{$b}{'NES'} <=> $gseaResults{$a}{'NES'}} keys %gseaResults) {
        $gsCount++;
        my $popupStrg = ($gseaResults{$geneSet}{'NES'} > 0)? "Enriched in $gseaInfo{'testCondName'}" : "Enriched in $gseaInfo{'refCondName'}";
        $catListStrg .= "             ['$geneSet', '$geneSet', null, '$popupStrg']";
        if ($gsCount == $geneSetsNb) {
            $catListStrg .=  "\n";
        } else {
            $catListStrg .=  ",\n";
        }
    }

    my $datasetValuesStrg = '';
    $gsCount = 0;
    foreach my $geneSet (sort {$gseaResults{$b}{'NES'} <=> $gseaResults{$a}{'NES'}} keys %gseaResults) {
        $gsCount++;
        $datasetValuesStrg .= qq
|               {value: $gseaResults{$geneSet}{'NES'}, 
                label: '$geneSet',
                id: '$geneSet',
                propValues: [$gseaResults{$geneSet}{'setSize'}, $gseaResults{$geneSet}{'qval'}, $gseaResults{$geneSet}{'ES'}]
               }|;
        if ($gsCount == $geneSetsNb) {
            $datasetValuesStrg .= "\n";
        } else {
            $datasetValuesStrg .= ",\n";
        }
    }


    print qq
|<DIV id="enrichmentDIV" style="float:center;padding:2px;"></DIV>
<SCRIPT type="text/javascript">
var EP;
function drawEnrichmentPlot() {
    EP = new idvis.categoryPlot({
        div:'enrichmentDIV',
        width:$width,height:$height,
        orientation:'horizontal', //vertical=default
        noColorEdition:true,
        categories:{
            title: 'Gene sets',
            list: [
                $catListStrg
            ]
        },
        valueAxes: {title: 'Normalised Enrichment Score'},
        categoryOnClick: categoryClicked,
        pointOnClick: pointClicked,
        scales: {
            color: [{map: 'q-value', gradient: 'RB', range: [0, $maxQvalue]}],
            pointSize: [{map: 'Set size', range: [0, $maxGeneSetSize]}]
        },
        datasets: [{   
            id: 'ds1',
            name: 'Enrichment',
            type: 'point',
            properties: ['Set size', 'q-value', 'Raw Enrichment Score'],
            appliedScales: ['q-value', 'Set size'],
            values: [
                $datasetValuesStrg
            ]
        }],
        exportAsImage: ['Export as image', 'EnrichmentPlot', './exportSVG.cgi']
    });

    EP.addFeature({type: 'line', axis: 'N', label: 'No enrichment', value: 0, color: '#ff8c00', width: 2, pattern: '.'});
    EP.draw();
}

function categoryClicked(catID) {
    ajaxGSDetails(catID, null);
}

function pointClicked(dsIdx, pLabel, pId) {
    ajaxGSDetails(pLabel, null);
}

drawEnrichmentPlot();
</SCRIPT>
<BR>
|;
}


sub printGeneSetDetails {  # Globals: $param, $gseaDir, %geneIsPresent, $projectID, $gseaID, %terms2Genes
    my ($geneSet, $refGseaInfo, $refGseaResults) = @_;
    my $searchedProtID = (param('SEARCH_PROT') && param('SEARCH_PROT') ne "null")? &promsMod::cleanNumericalParameters(param('SEARCH_PROT')) : 0;
    my %gseaInfo = %{$refGseaInfo};
    my %gseaResults = %{$refGseaResults};
    my $gsIsEnriched = (exists($gseaResults{$geneSet}))? 1 : 0;
    my $enrichedCond = ($gsIsEnriched)? ($gseaResults{$geneSet}{'ES'} >= 0)? $gseaInfo{'testCondName'} : $gseaInfo{'refCondName'} : "None";
    my $enrichedColor = ($gsIsEnriched)? ($gseaResults{$geneSet}{'ES'} >= 0)? $red : $green : $black;
    (my $geneSetStrg = $geneSet) =~ s/_/ /g;
    $geneSetStrg =~ s/([\w']+)/\u\L$1/g;  # capitalize first leter of each word

    # Retrieve analyses corresponding to GSEA quantification and all genes in the gene set from gmt file
    my $dbh = &promsConfig::dbConnect;
    my ($anaList) = $dbh->selectrow_array("SELECT GROUP_CONCAT(ID_ANALYSIS) FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION = $gseaInfo{'quantifID'} GROUP BY ID_QUANTIFICATION");
    $dbh->disconnect;

    my $gsDesc = ($gsIsEnriched)? $gseaResults{$geneSet}{'description'} : '';
    my @geneSetAll = @{$terms2Genes{$geneSet}};
    my $totalGSSize = scalar @geneSetAll;
    my @geneSetPresent = grep {exists($geneIsPresent{$_})} @geneSetAll;
    my $analysisGSSize = ($gsIsEnriched)? $gseaResults{$geneSet}{'setSize'} : scalar @geneSetPresent;
    my @geneSetAbsent = grep {!exists($geneIsPresent{$_})} @geneSetAll;
    my $leadEdgeProts = ($gsIsEnriched)? join(',', map($geneIsPresent{$_}, @{$gseaResults{$geneSet}{'coreGenes'}})) : '';

    # Popup quantification values for each gene
    sub makePopupStrg {
        my ($gene) = @_;
        my $popStrg = "<B>Protein" . "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" . ": " .
            "<FONT color=red>$gseaData{$geneIsPresent{$gene}}{'alias'}</FONT></B>" . "&nbsp;&nbsp;<BR>" .
            "<B>Gene Name" . "&nbsp;". ":</B> $gseaData{$geneIsPresent{$gene}}{'geneName'}" . "&nbsp;&nbsp;<BR>" .
            "<B>Pep. nb." . "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;". ":</B> " . 
            "$gseaData{$geneIsPresent{$gene}}{'pepNb'}" . "&nbsp;&nbsp;<BR>" .
            "<B>Ratio" . "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" . ":</B> " .
            "$gseaData{$geneIsPresent{$gene}}{'value'}" . "&nbsp;&nbsp;<BR>";
        if (exists($gseaData{$geneIsPresent{$gene}}{'pvalue'})) {
            $popStrg .= "<B>P-value" . "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" . ":</B> " .
                "$gseaData{$geneIsPresent{$gene}}{'pvalue'}" . "&nbsp;&nbsp;";
        }
    }
    my %popupDataStrgs = map { $_ => makePopupStrg($_) } @geneSetPresent;

    # Parameters for gene names display
    my $cellWidth = "80px";
    my $colNbGeneDisplay = 12;  # Nb of columns on which gene names are displayed (even nb for columns of same color)
    my $colSpanFullSet = (scalar @geneSetAll > $colNbGeneDisplay)? $colNbGeneDisplay : scalar @geneSetAll;
    my $colSpanLeadEdge = ($gsIsEnriched)? (scalar @{$gseaResults{$geneSet}{'coreGenes'}} > $colNbGeneDisplay)? $colNbGeneDisplay : scalar @{$gseaResults{$geneSet}{'coreGenes'}} : 0;

    # Parameters for graph display (if exists)
    my $size = 540;
    my $graphPath = "$promsPath{gsea_html}/project_$projectID/gsea_$gseaID/gsea_graph_$geneSet.png";
    my $hasGraph = (-e "$promsPath{'gsea'}/project_$projectID/gsea_$gseaID/gsea_graph_$geneSet.png") ? 1 : 0;

    print header(-'content-encoding' => 'no',-charset=>'utf-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
    <TITLE>Display Gene Set Details</TITLE>
    <LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY>
<CENTER>
<FONT class="title">GSEA Details for <FONT color="red">$geneSetStrg</FONT></FONT>
&nbsp;&nbsp;
<INPUT type="button" class="font11" value=" Hide " onclick="document.getElementById('divGeneSets').innerHTML='';document.getElementById('divGeneSets').style.display='none'">
<BR><BR>
<FORM name="protForm" method="POST">
<TABLE>
    <TR>
        <TD valign="top" align="center">
            <TABLE>
                <TR valign="top">
                    <TD style="padding-right:1cm">
                        <FONT class="title2">About the Gene Set</FONT>
                        <BR>
                        <TABLE id="detailsAbout" bgColor="$darkColor" border=0 cellspacing=1>
                            <TR>
                                <TH align="right" nowrap>&nbsp;&nbsp;Gene Set Term :&nbsp;</TH>
                                <TD class="lightBg" align="left" nowrap class="title3">&nbsp;&nbsp;$geneSet&nbsp;&nbsp;</TD>
                            </TR>
                            <TR>
                                <TH align="right" nowrap>&nbsp;&nbsp;Gene Set Databank :&nbsp;</TH>
                                <TD class="lightBg" align="left" nowrap>&nbsp;&nbsp;$gseaInfo{'gsDbName'}&nbsp;&nbsp;</TD>
                            </TR>
                            <TR>
                                <TH align="right" nowrap>&nbsp;&nbsp;Gene Set Description :&nbsp;</TH>
                                <TD class="lightBg" align="left" nowrap>&nbsp;&nbsp;$gsDesc&nbsp;&nbsp;</TD>
                            </TR>
                            <TR>
                                <TH align="right" nowrap>&nbsp;&nbsp;Total number of genes in Set :&nbsp;</TH>
                                <TD class="lightBg" align="left" nowrap>&nbsp;&nbsp;$totalGSSize&nbsp;&nbsp;</TD>
                            </TR>
                            <TR>
                                <TH align="right" nowrap>&nbsp;&nbsp;Genes from Set in Analysis :&nbsp;</TH>
                                <TD class="lightBg" align="left" nowrap>&nbsp;&nbsp;$analysisGSSize&nbsp;&nbsp;</TD>
                            </TR>
                        </TABLE>
                    </TD>
|;

    if ($gsIsEnriched) {
        print qq
|                    <TD>
                        <FONT class="title2">GSEA results</FONT>
                        <BR>
                        <TABLE id="detailsResults" bgcolor="$darkColor" border=0 cellspacing=1>
                            <TR>
                                <TH align="right" nowrap>&nbsp;&nbsp;Enrichment Score :&nbsp;</TH>
                                <TD class="lightBg" align="center" nowrap>&nbsp;&nbsp;${\sprintf("%.3g", $gseaResults{$geneSet}{'ES'})}&nbsp;&nbsp;</TD>
                            </TR>
                            <TR>
                                <TH align="right" nowrap>&nbsp;&nbsp;Normalized Enrichment Score :&nbsp;</TH>
                                <TD class="lightBg" align="center" nowrap>&nbsp;&nbsp;${\sprintf("%.3g", $gseaResults{$geneSet}{'NES'})}&nbsp;&nbsp;</TD>
                            </TR>
                            <TR>
                                <TH align="right" nowrap>&nbsp;&nbsp;P-value :&nbsp;</TH>
                                <TD class="lightBg" align="center" nowrap>&nbsp;&nbsp;${\sprintf("%.3e", $gseaResults{$geneSet}{'pval'})}&nbsp;&nbsp;</TD>
                            </TR>
                            <TR>
                                <TH align="right" nowrap>&nbsp;&nbsp;Adjusted p-value :&nbsp;</TH>
                                <TD class="lightBg" align="center" nowrap>&nbsp;&nbsp;${\sprintf("%.3e", $gseaResults{$geneSet}{'padj'})}&nbsp;&nbsp;</TD>
                            </TR>
                            <TR>
                                <TH align="right" nowrap>&nbsp;&nbsp;Q-value :&nbsp;</TH>
                                <TD class="lightBg" align="center" nowrap>&nbsp;&nbsp;${\sprintf("%.3e", $gseaResults{$geneSet}{'qval'})}&nbsp;&nbsp;</TD>
                            </TR>
                            <TR>
                                <TH align="right" nowrap>&nbsp;&nbsp;Rank <SUP onmouseover="popup(helpText.rank)" onmouseout="popout()">?</SUP>:&nbsp;</TH>
                                <TD class="lightBg" align="center" nowrap>&nbsp;&nbsp;$gseaResults{$geneSet}{'rank'}&nbsp;&nbsp;</TD>
                            </TR>
                            <TR>
                                <TH align="right" nowrap>&nbsp;&nbsp;Tags <SUP onmouseover="popup(helpText.tags)" onmouseout="popout()">?</SUP>:&nbsp;</TH>
                                <TD class="lightBg" align="center" nowrap>&nbsp;&nbsp;$gseaResults{$geneSet}{'tags'} %&nbsp;&nbsp;</TD>
                            </TR>
                            <TR>
                                <TH align="right" nowrap>&nbsp;&nbsp;List <SUP onmouseover="popup(helpText.list)" onmouseout="popout()">?</SUP>:&nbsp;</TH>
                                <TD class="lightBg" align="center" nowrap>&nbsp;&nbsp;$gseaResults{$geneSet}{'list'} %&nbsp;&nbsp;</TD>
                            </TR>
                            <TR>
                                <TH align="right" nowrap>&nbsp;&nbsp;Signal <SUP onmouseover="popup(helpText.signal)" onmouseout="popout()">?</SUP>:&nbsp;</TH>
                                <TD class="lightBg" align="center" nowrap>&nbsp;&nbsp;$gseaResults{$geneSet}{'signal'} %&nbsp;&nbsp;</TD>
                            </TR>
                        </TABLE>
                    </TD>
|;
    }
    print qq
|                </TR>
            </TABLE>
        </TD>
    </TR>
    <TR>
        <TD valign="top" align="left">
            <BR>
            <FONT class="title2">Gene Set composition</FONT>
            &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
            <INPUT type="button" valign="top" value="Check all" onclick="checkAllProteins(true)">
            <INPUT type="button" valign="top" value="Uncheck all" onclick="checkAllProteins(false)">
            &nbsp;&nbsp;&nbsp;&nbsp;
|;
    if ($gsIsEnriched) {
        print qq
|
            <INPUT type="button" valign="top" value="Check Lead. Edge" onclick="checkLEdgeProts(true, [$leadEdgeProts])">
            <INPUT type="button" valign="top" value="Uncheck Lead. Edge" onclick="checkLEdgeProts(false, [$leadEdgeProts])">
            &nbsp;&nbsp;&nbsp;&nbsp;
|;
    }
    print qq
|            <DIV id="saveProtDIV" style="display:none;"></DIV>
            <INPUT type="button" valign="top" id="saveFormBUTTON" value="Save proteins" onclick="ajaxManageSaveProteins('getThemes','PROT')"$disabSave>
            <BR>
            <TABLE border=0 cellspacing=1>
                <TR>
                    <TH class="darkBg" nowrap colspan=$colSpanFullSet>&nbsp;&nbsp;All Genes/Proteins in the Gene Set<SUP onmouseover="popup(helpText.genesUsed)" onmouseout="popout()">?</SUP>&nbsp;&nbsp;</TH>
                </TR>
                <TR>
|;
    my $cellCount = 0;
    my $cellColor = $lightColor;
    # First display genes that were in analysis
    foreach my $gene (sort{&promsMod::sortSmart($a, $b)} @geneSetPresent) {
        my $geneColor = ($searchedProtID && $geneIsPresent{$gene} == $searchedProtID)? "red" : "black";
        $cellColor = ($cellCount % 2)? $darkColor : $lightColor;
        my $geneStrg = "<INPUT type=\"checkbox\" name=\"chkProt\" id=\"chk_$geneIsPresent{$gene}\" value=\"$geneIsPresent{$gene}\">";
        $geneStrg .= "<A href=\"javascript:sequenceView($geneIsPresent{$gene}, '$anaList')\" onmouseover=\"popup('$popupDataStrgs{$gene}')\" onmouseout=\"popout()\"><B><FONT color=$geneColor>$gene</FONT></B></A>";
        unless ($cellCount == 0 || $cellCount % $colNbGeneDisplay) {
            print "</TR>\n<TR>\n";
        }
        print "<TD bgcolor=\"$cellColor\" valign=\"top\" align=\"left\" style=\"min-width:$cellWidth\" nowrap>$geneStrg</TD>\n";
        $cellCount++;
    }    
    unless ($colSpanFullSet < $colNbGeneDisplay) {
        while ($cellCount % $colNbGeneDisplay) {  # Finish the table line with alternating colors
            $cellColor = ($cellCount % 2)? $darkColor : $lightColor;
            print "<TD bgcolor=\"$cellColor\" valign=\"top\" align=\"left\" style=\"min-width:$cellWidth\" nowrap>&nbsp;</TD>\n";
            $cellCount++;
        }
    }

    # Then display genes from gene set that were not in analysis
    foreach my $gene (sort{&promsMod::sortSmart($a, $b)} @geneSetAbsent) {
        $cellColor = ($cellCount % 2)? $darkColor : $lightColor;
        unless ($cellCount % $colNbGeneDisplay) {
            print "</TR>\n<TR>\n";
        }
        print "<TD bgcolor=\"$cellColor\" valign=\"top\" align=\"center\" style=\"min-width:$cellWidth\" nowrap>$gene</TD>\n";
        $cellCount++;
    }
        
    unless ($colSpanFullSet < $colNbGeneDisplay) {
        while ($cellCount % $colNbGeneDisplay) {  # Finish the table line with alternating colors
            $cellColor = ($cellCount % 2)? $darkColor : $lightColor;
            print "<TD bgcolor=\"$cellColor\" valign=\"top\" align=\"center\" style=\"min-width:$cellWidth\" nowrap>&nbsp;</TD>\n";
            $cellCount++;
        }
    }
    
    print qq
|               </TR>
            </TABLE>
        </TD>
    </TR>
|;

    if ($gsIsEnriched) {
        print qq
|    <TR>
        <TD valign="top" align="left">
            <BR>
            <FONT class="title2">Leading Edge Genes/Proteins (most enriched in <FONT color="$enrichedColor">$enrichedCond</FONT>)</FONT>
            <BR>
            <TABLE border=0 cellspacing=1>
                <TR>
                    <TH class="darkBg" nowrap colspan=$colSpanLeadEdge>&nbsp;&nbsp;Core enrichment Genes/Proteins<SUP onmouseover="popup(helpText.leadingEdge)" onmouseout="popout()">?</SUP>&nbsp;&nbsp;</TH>
                </TR>
                <TR>
|;
        $cellCount = 0;
        $cellColor = $lightColor;
        foreach my $gene (@{$gseaResults{$geneSet}{'coreGenes'}}) {
            my $geneColor = ($searchedProtID && $geneIsPresent{$gene} == $searchedProtID)? "red" : "black";
            $cellColor = ($cellCount % 2)? $darkColor : $lightColor;
            unless ($cellCount == 0 || $cellCount % $colNbGeneDisplay) {
                print "</TR>\n<TR>\n";
            }
            print qq
|                   <TD bgcolor="$cellColor" valign="top" align="center" style="min-width:$cellWidth" nowrap>
                        <A href="javascript:sequenceView($geneIsPresent{$gene}, '$anaList')" onmouseover="popup('$popupDataStrgs{$gene}')" onmouseout="popout()">
                            <B><FONT color="$geneColor">$gene</FONT></B>
                        </A>
                    </TD>
|;
            $cellCount++;
        }
        unless ($colSpanLeadEdge < $colNbGeneDisplay) {
            while ($cellCount % $colNbGeneDisplay) {  # Finish the table line with alternating colors
                $cellColor = ($cellCount % 2)? $darkColor : $lightColor;
                print qq
|                   <TD bgcolor="$cellColor" valign="top" align="center" style="min-width:$cellWidth" nowrap>&nbsp;</TD>
|;
                $cellCount++;
            }
        }
        print qq
|               </TR>
            </TABLE>
        </TD>
    </TR>
|;
    } else {
        print qq
|   <TR>
        <TD valign="top" align="center">
            <BR>
            <FONT class="title3">Gene Set <FONT color="red"><B>NOT</B></FONT> enriched in analysis. Limited information available.</FONT>
            <BR>
        </TD>
    </TR>
|;
    }
    
    if ($hasGraph) {
        print qq
|   <TR>
        <TD valign="top" align="center">
            <BR>
            <FONT class="title2">Gene Set GSEA Graph</FONT>&nbsp;&nbsp;&nbsp;&nbsp;
            <A href="$graphPath" download="gsea_graph_$geneSet.png">
                <INPUT type="button" class="font11" value=" Export as PNG ">
            </A>
            <BR>
            <IMG src="$graphPath" width=$size heigth=$size>;
        </TD>
    </TR>
|;
    }
    print qq
|</TABLE>
</FORM>
</CENTER>
</BODY>
</HTML>
|;

}


sub showGraph {  # Globals: %promsPath, $projectID, $gseaID
    my ($geneSet) = @_;
    my $graphPath;
    if ($geneSet) {
        $graphPath = "$promsPath{gsea_html}/project_$projectID/gsea_$gseaID/gsea_graph_$geneSet.png";
    } else {
        $graphPath = "$promsPath{gsea_html}/project_$projectID/gsea_$gseaID/gsea_graph_dotplot.png";
    }
    my $size = ($geneSet)? 540 : 650;
    print header(-'content-encoding' => 'no',-charset=>'utf-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
    <TITLE>Display GSEA graph</TITLE>
    <LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY>
<CENTER>
    <FONT class="title">GSEA graph</FONT>
    &nbsp;&nbsp;
    <A href="$graphPath" download="gsea_graph_$geneSet.png">
        <INPUT type="button" class="font11" value=" Export as PNG ">
    </A>
    &nbsp;&nbsp;
    <INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('infoDIV').innerHTML='';document.getElementById('displayDIV').style.display='none'">
    <BR><BR>
    <IMG src="$graphPath" width=$size heigth=$size>
</CENTER>
</BODY>
</HTML>
|;

}


sub displayProtGeneSets {  # Globals: %protResults, %genes2Terms, %gseaResults, %gseaData, %gseaInfo, $lightColor, $darkColor
    my ($refProtIDs) = @_;
    my @protIDs = @{$refProtIDs};
    
    my ($protInfoRef, $masterProtsRef) = &getProtInfo(\@protIDs);
    my %protInfo = %{$protInfoRef};
    my %masterProts = %{$masterProtsRef};

    my $exportFilterParam = "";
    $exportFilterParam .= "&protID=" . join("&protID=", @protIDs) if (@protIDs);

    my ($pepStrg, $pepPopup);
    if ($gseaInfo{'pepWeights'} eq 'msms') {
        $pepStrg = "&nbsp;MS/MS&nbsp;<BR>&nbsp;Peptides&nbsp;";
        $pepPopup = "Truly identified peptides found in quantification and used to weight protein scores in GSEA";
    } elsif ($gseaInfo{'pepWeights'} eq 'distinct') {
        $pepStrg = "&nbsp;Dist.&nbsp;<BR>&nbsp;Peptides&nbsp;";
        $pepPopup = "Distinct Peptides found in quantification and used to weight protein scores in GSEA";
    } else {
        $pepStrg = "&nbsp;Peptides&nbsp;";
        $pepPopup = "All Peptides used in quantification";
    }

    my $titleStrg = "Protein Information";
    print header(-'content-encoding' => 'no',-charset=>'utf-8');
    warningsToBrowser(1);

    print qq
|<HTML>
<HEAD>
    <TITLE>Protein Information</TITLE>
    <LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY>
<CENTER>
<BR><BR><BR>
<FONT class="title1">$titleStrg</FONT>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<INPUT type="button" class="font11" value="Export List" onclick="window.location='$promsPath{cgi}/displayGSEA.cgi?ID=$gseaID&ACTION=export&XLS=proteins&$exportFilterParam&VIEW=$view';">
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<INPUT type="button" class="font11" value=" Hide " onclick="document.getElementById('protTable').innerHTML='';document.getElementById('protTable').style.display='none'">
<BR><BR>
<DIV id="tableDiv">
    <TABLE border=0 cellspacing=0>
    <TR bgcolor=$darkColor>
        <TH class="rbBorder" style="min-width:400px;" align="center">&nbsp;Gene Set(s)&nbsp;</TH>
        <TH class="rbBorder" style="min-width:125px;" align="center">&nbsp;Protein&nbsp;</TH>
        <TH class="rbBorder" align="center">&nbsp;Uniprot&nbsp;ACC&nbsp;</TH>
        <TH class="rbBorder" align="center">&nbsp;Gene&nbsp;name&nbsp;</TH>
        <TH class="rbBorder" style="min-width:130px;" align="center">&nbsp;Quantif. value&nbsp;</TH>
        <TH class="rbBorder" style="min-width:130px;" align="center">&nbsp;Quantif. p-value&nbsp;</TH>
        <TH class="rbBorder" align="center">
            <FONT style="cursor:help;" onmouseover="popup('$pepPopup');" onmouseout="popout();">
                $pepStrg
            </FONT>
        </TH>
        <TH class="bBorder" width=700 align="center">Description - Species</TH>
    </TR>
|;

    my $bgColor = $lightColor;
    foreach my $protID (sort { $gseaData{$b}{'pepNb'} <=> $gseaData{$a}{'pepNb'} } @protIDs) {
        my $gene = $protResults{$protID}{'gene'};

        my $geneSetsStr = '';
        if (defined $genes2Terms{$gene}) {
            foreach my $geneSet (@{$genes2Terms{$gene}}) {
                my $gsEnrichment = ($gseaResults{$geneSet})? ($gseaResults{$geneSet}{'NES'} > 0)? 
                    "<IMG src=\"$promsPath{images}/up_red.png\" width=\"11\" height=\"11\"/>" : 
                    "<IMG src=\"$promsPath{images}/down_green.png\" width=\"11\" height=\"11\"/>" : 
                    "<span style=\"padding-left:11px\"></span>";
                $geneSetsStr .= "<BR>" if ($geneSetsStr);
                $geneSetsStr .= "&nbsp;&nbsp;$gsEnrichment&nbsp;&nbsp;<A href=\"javascript:ajaxGSDetails('$geneSet', $protID);\">$geneSet</A>";
            }
        }
        my $geneStrg = '-';
        if ($protInfo{$protID}{'MasterProt'} && $masterProts{$protInfo{$protID}{'MasterProt'}}{'gene'}[0]) {
            my $geneName = @{$masterProts{$protInfo{$protID}{'MasterProt'}}{'gene'}}[0];
            if (scalar @{$masterProts{$protInfo{$protID}{'MasterProt'}}{'gene'}} > 1) {
                my $genesList = "<BR>&nbsp;&nbsp;-" . join('<BR>&nbsp;&nbsp;-', @{$masterProts{$protInfo{$protID}{'MasterProt'}}{'gene'}});
                $geneStrg = "<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U>" . $genesList . "</FONT>')\" onmouseout=\"popout()\">$geneName</A>";
            } else {
                $geneStrg = $geneName;
            }
        }
        my $uniAccStrg = '-';
        if ($protInfo{$protID}{'MasterProt'} && $masterProts{$protInfo{$protID}{'MasterProt'}}{'uniACC'}[0]) {
            my $uniACC = @{$masterProts{$protInfo{$protID}{'MasterProt'}}{'uniACC'}}[0];
            $uniAccStrg = $uniACC;
        }
        my $quantifVal = $gseaData{$protID}{'value'};
        my $pval = $gseaData{$protID}{'pvalue'};
        $pval = '-' unless ($pval && $pval ne "NA");
        my $pepNb = $gseaData{$protID}{'pepNb'};

        print qq 
|   <TR bgcolor=$bgColor>
        <TH align="left" style="min-width:400px; padding: 0 20px 0 5px;" valign="top">$geneSetsStr</TH>
        <TH align="center" valign="center">&nbsp;<A href="javascript:sequenceView($protID,$protInfo{$protID}{AnaID})">$protInfo{$protID}{'Name'}</A>&nbsp;</TH>
        <TH align="center" valign="center">&nbsp;$uniAccStrg&nbsp;</TH>
        <TH align="center" valign="center">&nbsp;$geneStrg&nbsp;</TH>
        <TH align="center" valign="center">&nbsp;$quantifVal&nbsp;</TH>
        <TH align="center" valign="center">&nbsp;$pval&nbsp;</TH>
        <TH align="center" valign="center">&nbsp;$pepNb&nbsp;</TH>
        <TD>&nbsp;$protInfo{$protID}{'Des'} <FONT class="org">$protInfo{$protID}{'Org'}</FONT></TD>
    </TR>
|;
        $bgColor = ($bgColor eq $lightColor)? $darkColor : $lightColor;
    }

    print qq
|    </TABLE>
<BR><BR>
</CENTER>
</DIV>
</BODY>
</HTML>
|;
}


sub getGSEAInfo {  # Globals: $dbh, $expID, $param, $upDate, $user
    my ($refGseaInfo) = @_;

    # Experiment info
    ($refGseaInfo->{'expName'}) = $dbh->selectrow_array("SELECT NAME FROM EXPERIMENT WHERE ID_EXPERIMENT = $expID");

    # Parameter string info
    ($refGseaInfo->{'quantifID'})      = ($param =~ /quantification=#([^;]+)/);
    ($refGseaInfo->{'targetPos'})      = ($param =~ /targetPos=([^;]+)/);
    ($refGseaInfo->{'isRatio'})        = ($param =~ /isRatio=([^;]+)/);
    ($refGseaInfo->{'geneSetsID'})     = ($param =~ /geneSetsDB=#([^;]+)/);
    ($refGseaInfo->{'gsDbVersion'})    = ($param =~ /geneSetsVersion=([^;]+)/);
    ($refGseaInfo->{'pepWeights'})     = ($param =~ /weightByPep=([^;]+)/);
    ($refGseaInfo->{'pvalCutoff'})     = ($param =~ /pvalCutoff=([^;]+)/);
    ($refGseaInfo->{'pvalCorrection'}) = ($param =~ /pvalCorr=([^;]+)/);

    # Info about parent quantification
    my ($designName, $quantifName, $quantifAnnot) = $dbh->selectrow_array("SELECT D.NAME, Q.NAME, Q.QUANTIF_ANNOT FROM QUANTIFICATION Q INNER JOIN DESIGN D ON Q.ID_DESIGN = D.ID_DESIGN WHERE ID_QUANTIFICATION = $refGseaInfo->{'quantifID'}");
    my $sthExpCond = $dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION = ?");

    $refGseaInfo->{'quantifStrg'} = "$designName > $quantifName";

    if ($refGseaInfo->{'isRatio'}) {  # Only case possible for now, prevision for GSEA on one state (after abundance quantif)
        my ($testCondID, $refCondID);
        my ($testCondName, $refCondName);
        my $supRatioTag = '';  # In case of SuperRatio 
        foreach my $features (split('::', $quantifAnnot)) {
            if ($features =~ /^RATIOS=(.*)$/) {
                my $ratioStrg = $1;
                $ratioStrg =~ s/#//g;  # Clean ID tags
                my @allRatios = split(';', $ratioStrg);
                my $interestRatio = $allRatios[$refGseaInfo->{'targetPos'} - 1];  # Retrieve state IDs of the wanted ratio
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

        $refGseaInfo->{'quantifStrg'} .= " : ${testCondName}${supRatioTag}/${refCondName}${supRatioTag}";
        $refGseaInfo->{'testCondID'}   = $testCondID;
        $refGseaInfo->{'testCondName'} = $testCondName;
        $refGseaInfo->{'refCondID'}    = $refCondID;
        $refGseaInfo->{'refCondName'}  = $refCondName;
        $refGseaInfo->{'supRatioTag'}  = $supRatioTag;
    }
    $sthExpCond->finish;

    # Info about Gene Sets databank
    my ($gsName, $gsDes, $gsProject, $gsSpeciesID, $gsSpecies, $gsIdentID, $gsIdentifier, $gsFile, $gsFileType) = $dbh->selectrow_array(
        "SELECT GS.NAME, GS.DES, GS.ID_PROJECT, GS.ID_SPECIES, S.SCIENTIFIC_NAME, I.ID_IDENTIFIER, I.NAME, GS.GENESETS_FILE, GS.GENESETS_TYPE
        FROM GENESETS GS 
        LEFT JOIN SPECIES S ON GS.ID_SPECIES = S.ID_SPECIES
        LEFT JOIN IDENTIFIER I ON GS.ID_IDENTIFIER = I.ID_IDENTIFIER
        WHERE ID_GENESETS = $refGseaInfo->{'geneSetsID'}"
    );
    $refGseaInfo->{'gsDbName'}       = $gsName;
    $refGseaInfo->{'gsDbDes'}        = $gsDes;
    $refGseaInfo->{'gsDbProject'}    = ($gsProject)? $gsProject : '';
    $refGseaInfo->{'gsDbSpeciesID'}  = ($gsSpeciesID)? $gsSpeciesID : undef;
    $refGseaInfo->{'gsDbSpecies'}    = ($gsSpecies)? $gsSpecies : "No specific organism";
    $refGseaInfo->{'gsDbIdentID'}    = ($gsIdentID)? $gsIdentID : undef;
    $refGseaInfo->{'gsDbIdentifier'} = ($gsIdentifier)? $gsIdentifier : "No specific identifier";
    $refGseaInfo->{'gsDbFile'}       = "$promsPath{'genesets'}/$gsFile";
    $refGseaInfo->{'gsDbFileType'}   = $gsFileType;

    # Creation date and user info
    $refGseaInfo->{'updateDate'} = &promsMod::formatDate($upDate);
    if ($user) {
        my ($userName) = $dbh->selectrow_array("SELECT USER_NAME FROM USER_LIST WHERE ID_USER = '$user'");
        $refGseaInfo->{'user'} = $userName || $user;
    } else {
        $refGseaInfo->{'user'} = '?';
    }

}


sub parseGSEAResults {
    my ($refGseaResults, $refGeneIsPresent, $refProtResults, $results, $genesUsed) = @_;

    open(RESULTS, $results);
    # ID Description setSize enrichmentScore NES pvalue p.adjust qvalues rank leading_edge core_enrichment
    while (my $line = <RESULTS>) {
        next if ($. == 1);  # Skip headers
        chomp $line;
        my @values = split(/\t/, $line);
        my $geneSet = $values[0];
        $refGseaResults->{$geneSet}{'description'} = ($values[1] ne $geneSet)? $values[1] : '-';
        $refGseaResults->{$geneSet}{'setSize'} = $values[2];
        $refGseaResults->{$geneSet}{'ES'} = $values[3];
        $refGseaResults->{$geneSet}{'NES'} = $values[4];
        $refGseaResults->{$geneSet}{'pval'} = $values[5];
        $refGseaResults->{$geneSet}{'padj'} = $values[6];
        $refGseaResults->{$geneSet}{'qval'} = $values[7];
        $refGseaResults->{$geneSet}{'rank'} = $values[8];
        my @leadingEdge = split(/,\s/, $values[9]);
        ($refGseaResults->{$geneSet}{'tags'}) = ($leadingEdge[0] =~ /tags=(\d+)%/);
        ($refGseaResults->{$geneSet}{'list'}) = ($leadingEdge[1] =~ /list=(\d+)%/);
        ($refGseaResults->{$geneSet}{'signal'}) = ($leadingEdge[2] =~ /signal=(\d+)%/);
        my @coreGenes = split('/', $values[10]);
        $refGseaResults->{$geneSet}{'coreGenes'} = \@coreGenes;
    }
    close RESULTS;

    # Retrieve all genes that were used for GSEA
    open(DATA, $genesUsed);
    while (my $line = <DATA>) {
        next if ($. == 1);  # Header line
        chomp $line;
        my @lineData = split("\t", $line);  # Protein_id Gene_name
        my $lineLength = scalar @lineData;
        $refGeneIsPresent->{$lineData[1]} = $lineData[0];
        $refProtResults->{$lineData[0]}{'gene'} = $lineData[1];
        if ($lineLength > 2) {
            $refProtResults->{$lineData[0]}{'score'} = $lineData[2];
        }
    }
    close DATA;
}


sub parseGSEAData {  # Globals: %gseaInfo
    my ($refGseaData, $originalDataFile) = @_;

    my $hasPvalue;
    open(DATA, $originalDataFile);
    # Protein_ID Gene_name Value Peptides_nb [Pvalue]
    while (my $line = <DATA>) {
        chomp $line;
        my @lineData = split("\t", $line);
        if ($. == 1) {  # Header line
            $hasPvalue = (scalar @lineData > 4)? 1 : 0;
            next;
        }
        my $protID       = $lineData[0];
        my $gsIdentifier = $lineData[1];
        my $log2value    = $lineData[2];
        my $pepNb        = $lineData[3];
        my $pvalue       = $lineData[4] if ($hasPvalue);

        $refGseaData->{$protID}{'geneSetId'} = $gsIdentifier;  # Identifier in gene set (not ID of the gene set)
        $refGseaData->{$protID}{'pepNb'} = $pepNb;

        if ($gseaInfo{'isRatio'}) {
            my $log2InfRatio = log(1000) / log(2);
            my $tol = 1e-8;  # Tolerance to check if ratio is infinite
            my $isInf = (abs(abs($log2value) - $log2InfRatio) < $tol) ? 1 : 0;  # Ratio is infinite ?
            if ($log2value < 0) {
                if ($isInf) {
                    $refGseaData->{$protID}{'value'} = "1/Inf";
                } else {
                    $refGseaData->{$protID}{'value'} = "1/" . sprintf("%.2f", 2 ** (-$log2value));
                }
            } else {
                if ($isInf) {
                    $refGseaData->{$protID}{'value'} = "Inf";
                } else {
                    $refGseaData->{$protID}{'value'} = sprintf("%.2f", 2 ** ($log2value));
                }
            }
            $refGseaData->{$protID}{'pvalue'} = ($pvalue)? ($pvalue eq "NA")? "NA" : sprintf("%.2e", $pvalue) : '';
        } else {
            $refGseaData->{$protID}{'value'} = sprintf("%.2e", 2 ** $log2value);
        }
    }
    close DATA;

    # Get protein and gene name for each prot in data file
    my $protIDStrg = join(',', keys %{$refGseaData});
    my $sthGetIdentifiers = $dbh->prepare("SELECT P.ID_PROTEIN, P.ALIAS, MI.VALUE 
        FROM PROTEIN P 
        LEFT JOIN MASTERPROT_IDENTIFIER MI ON P.ID_MASTER_PROTEIN = MI.ID_MASTER_PROTEIN
        LEFT JOIN IDENTIFIER I ON MI.ID_IDENTIFIER = I.ID_IDENTIFIER
        WHERE P.ID_PROTEIN IN ($protIDStrg) AND 
        I.CODE = 'GN' AND MI.IDENT_RANK = 1");
    $sthGetIdentifiers->execute;
    while (my ($protID, $alias, $geneName) = $sthGetIdentifiers->fetchrow_array) {
        $refGseaData->{$protID}{'alias'} = $alias;
        $refGseaData->{$protID}{'geneName'} = $geneName;
    }
    $sthGetIdentifiers->finish; 
}


sub linkTermsAndGenes {  # Globals: %promsPath, %gseaInfo
    my ($refTerms2Genes, $refGenes2Terms) = @_;
    if ($gseaInfo{'gsDbFileType'} eq "GMT") {
        open(GMTFILE, "$gseaInfo{'gsDbFile'}");
        while (my $line = <GMTFILE>) {
            chomp $line;
            my @geneSetAll = split(/\t/, $line);
            my $term = shift @geneSetAll;
            shift @geneSetAll;  # shift twice because .gmt files have 2 columns before genes
            $refTerms2Genes->{$term} = \@geneSetAll;
        }
        close GMTFILE;
    }

    foreach my $term (keys %{$refTerms2Genes}) {
        foreach my $gene (@{$refTerms2Genes->{$term}}) {
            push @{$refGenes2Terms->{$gene}}, $term;
        }
    }

}


sub getAllIdentifierTypes {
    my %allIdentifiers;
    my $dbh = &promsConfig::dbConnect;
    my $sthGetIdent = $dbh->prepare("SELECT ID_IDENTIFIER, NAME FROM IDENTIFIER");
    $sthGetIdent->execute;
    while (my ($identID, $identName) = $sthGetIdent->fetchrow_array) {
        $allIdentifiers{$identID} = $identName;
    }
    $sthGetIdent->finish;
    $dbh->disconnect;
    return %allIdentifiers;
}


sub searchConvertIdentifier {  # Globals: %gseaInfo, $projectID
    my $searchText = lc(param('SEARCH'));
    my $searchIdType = param('SEARCH_ID_TYPE');
    return (undef, undef) unless ($searchText && $searchIdType);

    my $dbh = &promsConfig::dbConnect;

    # Check if needed to restrict to specific species according to Gene Set Databank
    my @restrictSpecies = ('', '');
    if ($gseaInfo{'gsDbSpeciesID'}) {
        @restrictSpecies = ("LEFT JOIN MASTER_PROTEIN MP ON P.ID_MASTER_PROTEIN = MP.ID_MASTER_PROTEIN",
                            "(MP.ID_SPECIES IS NULL OR MP.ID_SPECIES = $gseaInfo{'gsDbSpeciesID'}) AND ");
    }
    
    # Search protein in DB based on projectID, speciesID and search text (exact match)
    my ($searchProtID) = $dbh->selectrow_array("SELECT ID_PROTEIN 
        FROM PROTEIN P $restrictSpecies[0]
        LEFT JOIN MASTERPROT_IDENTIFIER MI ON P.ID_MASTER_PROTEIN = MI.ID_MASTER_PROTEIN
        LEFT JOIN IDENTIFIER I ON MI.ID_IDENTIFIER = I.ID_IDENTIFIER
        WHERE P.ID_PROJECT = $projectID AND $restrictSpecies[1] (
            (LOWER(P.IDENTIFIER) = '$searchText')
            OR
            (LOWER(P.ALIAS) = '$searchText')
            OR
            (I.ID_IDENTIFIER = $searchIdType AND LOWER(MI.VALUE) = '$searchText')
        )"
    );
    
    unless ($searchProtID) {
        $dbh->disconnect;
        return (undef, undef);
    }

    # Find Identifier of the searched protein in Gene Set Databank format (not always the same as the search text)
    my ($searchProtGsId) = $dbh->selectrow_array("SELECT MI.VALUE 
        FROM PROTEIN P 
        LEFT JOIN MASTERPROT_IDENTIFIER MI ON P.ID_MASTER_PROTEIN = MI.ID_MASTER_PROTEIN
        WHERE P.ID_PROTEIN = $searchProtID AND MI.ID_IDENTIFIER = $gseaInfo{'gsDbIdentID'} AND MI.IDENT_RANK = 1"
    );

    $dbh->disconnect;
    return ($searchProtID, $searchProtGsId);
}


sub exportGSEAResults {  # Globals: %gseaResults, %gseaInfo, $name, $desc
    my $workbook = Spreadsheet::WriteExcel->new("-");
    my $worksheet1 = $workbook->add_worksheet('Results');
    my $worksheet2 = $workbook->add_worksheet('Parameters');

    print header(-type=>"application/vnd.ms-excel",-attachment=>"$name.xls");

    # Cell formatting
    my $headerMergeFormat = $workbook->add_format(align => 'center', valign => 'vcenter');
    $headerMergeFormat->set_bold();
    $headerMergeFormat->set_text_wrap();
    my $headerFormat = $workbook->add_format(align => 'center', valign => 'vcenter');
    $headerFormat->set_bold();
    $headerFormat->set_text_wrap();
    my $formatExp = $workbook->add_format();
    $formatExp->set_num_format("0.00E+00");
    
    ### Writing Results sheet
    # Format columns
    $worksheet1->set_column("A:A", 40);
    $worksheet1->set_column("B:B", 20);
    $worksheet1->set_column("C:C", 15);
    $worksheet1->set_column("E:E", 15);
    $worksheet1->set_column("L:L", 15);

    ## Write column headers
    $worksheet1->merge_range("A1:A2", "Gene Set Term", $headerMergeFormat);
    $worksheet1->merge_range("B1:B2", "Gene Set Description", $headerMergeFormat);
    $worksheet1->merge_range("C1:C2", "Nb of Genes\nfound in Set", $headerMergeFormat);
    $worksheet1->merge_range("D1:E1", "Enrichment Score", $headerMergeFormat);
    $worksheet1->write_string("D2", "Raw", $headerFormat);
    $worksheet1->write_string("E2", "Normalized", $headerFormat);
    $worksheet1->merge_range("F1:G1", "P-values", $headerMergeFormat);
    $worksheet1->write_string("F2", "Raw", $headerFormat);
    $worksheet1->write_string("G2", "Adjusted", $headerFormat);
    $worksheet1->merge_range("H1:H2", "Q-values", $headerMergeFormat);
    $worksheet1->merge_range("I1:I2", "Rank", $headerMergeFormat);
    $worksheet1->merge_range("J1:J2", "Tags (%)", $headerMergeFormat);
    $worksheet1->merge_range("K1:K2", "List (%)", $headerMergeFormat);
    $worksheet1->merge_range("L1:L2", "Signal (%)", $headerMergeFormat);

    ## Write Gene Sets results
    my $i = 1;
    my $maxLEdgeNb = 0;
    foreach my $geneSet (sort {$gseaResults{$a}{'padj'} <=> $gseaResults{$b}{'padj'} || abs($gseaResults{$b}{'NES'}) <=> abs($gseaResults{$a}{'NES'})} keys %gseaResults) {
        $i++;
        (my $geneSetStrg = $geneSet) =~ s/_/ /g;        
        my $gsDesc        = $gseaResults{$geneSet}{'description'};
        my $setSize     = $gseaResults{$geneSet}{'setSize'};
        my $ES          = sprintf("%.3f", $gseaResults{$geneSet}{'ES'});
        my $NES         = sprintf("%.3f", $gseaResults{$geneSet}{'NES'});
        my $pval        = sprintf("%.3e", $gseaResults{$geneSet}{'pval'});
        my $padj        = sprintf("%.3e", $gseaResults{$geneSet}{'padj'});
        my $qval        = sprintf("%.3e", $gseaResults{$geneSet}{'qval'});
        my @coreGenes   = @{$gseaResults{$geneSet}{'coreGenes'}};
        
        $worksheet1->write_string($i, 0, $geneSetStrg);
        $worksheet1->write_string($i, 1, $gsDesc);
        $worksheet1->write_number($i, 2, $setSize);
        $worksheet1->write_number($i, 3, $ES);
        $worksheet1->write_number($i, 4, $NES);
        $worksheet1->write_number($i, 5, $pval, $formatExp);
        $worksheet1->write_number($i, 6, $padj, $formatExp);
        $worksheet1->write_number($i, 7, $qval, $formatExp);
        $worksheet1->write_number($i, 8, $gseaResults{$geneSet}{'rank'});
        $worksheet1->write_number($i, 9, $gseaResults{$geneSet}{'tags'});
        $worksheet1->write_number($i, 10, $gseaResults{$geneSet}{'list'});
        $worksheet1->write_number($i, 11, $gseaResults{$geneSet}{'signal'});
        my $j = 11;  # Leading Edge genes start at column M
        foreach my $gene (@coreGenes) {
            $j++;
            $worksheet1->write_string($i, $j, $gene);
        }
        $maxLEdgeNb = ($maxLEdgeNb < $j)? $j : $maxLEdgeNb;
    }

    # Write this after loops because need to catch the max number of leading edge genes
    $worksheet1->merge_range(0, 12, 1, $maxLEdgeNb, "Leading Edge Genes", $headerMergeFormat);

    ### End of Results sheet

    # Weights on protein scores/ratios ?
    my $pepWeightsStrg = "Protein scores are ";
    if ($gseaInfo{'pepWeights'} eq "none") {
        $pepWeightsStrg = "Not weighted with peptides number";
    } else {
        $pepWeightsStrg = "weighted with number of ";
        if ($gseaInfo{'pepWeights'} eq "msms") {
            $pepWeightsStrg .= "Truly identified";
        } else {
            $pepWeightsStrg .= ucfirst($gseaInfo{'pepWeights'});
        }
        $pepWeightsStrg .= " peptides";
    }

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

    ### Write Parameters sheet (corresponds to GSEA summary)
    # Format columns
    $worksheet2->set_column("A:A", 30);
    $worksheet2->set_column("B:B", 30);
    $worksheet2->set_column("C:C", 90);
    
    # Write row headers
    $worksheet2->merge_range("A1:B1", "GSEA name", $headerMergeFormat);
    $worksheet2->merge_range("A2:B2", "Description", $headerMergeFormat);
    $worksheet2->merge_range("A3:B3", "Experiment name", $headerMergeFormat);
    $worksheet2->merge_range("A4:B4", "Parent quantification", $headerMergeFormat);
    $worksheet2->merge_range("A5:A9", "Gene sets databank", $headerMergeFormat);
    $worksheet2->write_string("B5", "Name", $headerFormat);
    $worksheet2->write_string("B6", "Version Date", $headerFormat);
    $worksheet2->write_string("B7", "Organism", $headerFormat);
    $worksheet2->write_string("B8", "Identifier type", $headerFormat);
    $worksheet2->write_string("B9", "Gene sets description", $headerFormat);
    $worksheet2->merge_range("A10:A12", "Parameters", $headerMergeFormat);
    $worksheet2->write_string("B10", "Protein scores are :", $headerFormat);
    $worksheet2->write_string("B11", "Adjusted p-value cutoff", $headerFormat);
    $worksheet2->write_string("B12", "P-value correction", $headerFormat);
    $worksheet2->merge_range("A13:B13", "Creation date", $headerMergeFormat);
    
    # Write parameters values
    $worksheet2->write_string("C1", $name);
    $worksheet2->write_string("C2", $desc);
    $worksheet2->write_string("C3", $gseaInfo{'expName'});
    $worksheet2->write_string("C4", $gseaInfo{'quantifStrg'});
    $worksheet2->write_string("C5", $gseaInfo{'gsDbName'});
    $worksheet2->write_string("C6", $gseaInfo{'gsDbVersion'});
    $worksheet2->write_string("C7", $gseaInfo{'gsDbSpecies'});
    $worksheet2->write_string("C8", $gseaInfo{'gsDbIdentifier'});
    $worksheet2->write_string("C9", $gseaInfo{'gsDbDes'});
    $worksheet2->write_string("C10", $pepWeightsStrg);
    $worksheet2->write_string("C11", $gseaInfo{'pvalCutoff'});
    $worksheet2->write_string("C12", $fullPvalCorr{$gseaInfo{'pvalCorrection'}});
    $worksheet2->write_string("C13", "$gseaInfo{'updateDate'} by $gseaInfo{'user'}");

    $worksheet2->write_string("A15", "Analysis performed with clusterProfiler (Yu G. et al. OMICS, 2012), based on the GSEA method (Subramanian A. et al. PNAS, 2005)");

    $workbook->close();
    exit;
}


sub exportProtList {  # Globals: %protResults, %genes2Terms, %gseaData, %gseaInfo
    my ($refProtIDs) = @_;
    my @protIDs = @{$refProtIDs};

    my ($protInfoRef, $masterProtsRef) = &getProtInfo(\@protIDs);
    my %protInfo = %{$protInfoRef};
    my %masterProts = %{$masterProtsRef};

    my $pepStrg = ($gseaInfo{'pepWeights'} eq "msms")? "MS/MS peptides" : ($gseaInfo{'pepWeights'} eq "distinct")? "Dist. peptides" : "All peptides";

    print header(-type => "application/vnd.ms-excel", -attachment => "Protein_list.xls");
    warningsToBrowser(1);
    my $workbook = Spreadsheet::WriteExcel->new("-");
    my $worksheet = $workbook->add_worksheet();

    # Cell formatting
    my $headerMergeFormat = $workbook->add_format(align => 'center', valign => 'vcenter');
    $headerMergeFormat->set_bold();
    my $headerFormat = $workbook->add_format(align => 'center');
    $headerFormat->set_bold();
    $headerFormat->set_text_wrap();
    my $contentFormat = $workbook->add_format(align => 'left',  valign => 'vcenter');
    $contentFormat->set_text_wrap();

    # Writting column headers
    $worksheet->write_string("A1", "Gene Set(s)", $headerFormat);
    $worksheet->write_string("B1", "Proteins", $headerFormat);
    $worksheet->write_string("C1", "Uniprot ACC", $headerFormat);
    $worksheet->write_string("D1", "Gene Name", $headerFormat);
    $worksheet->write_string("E1", "Quantif. value", $headerFormat);
    $worksheet->write_string("F1", "Quantif. p-value", $headerFormat);
    $worksheet->write_string("G1", $pepStrg, $headerFormat);
    $worksheet->write_string("H1", "Description - Species", $headerFormat);

    my $row=1;
    my $maxCharGS = 0;
    my $maxCharDesc = 0;

    foreach my $protID (sort { $gseaData{$b}{'pepNb'} <=> $gseaData{$a}{'pepNb'} } @protIDs) {
        my $gene = $protResults{$protID}{'gene'};

        my $geneSetCount = 0;
        my $geneSetsStr = '';
        if (defined $genes2Terms{$gene}) {
            $geneSetCount = scalar @{$genes2Terms{$gene}};
            foreach my $geneSet (@{$genes2Terms{$gene}}) {
                $geneSetsStr .= ($geneSetsStr)? "\n$geneSet" : $geneSet;
                $maxCharGS = length($geneSet) if (length($geneSet) > $maxCharGS);
            }
        }
        $geneSetsStr = '-' unless ($geneSetsStr);

        my $geneStrg = '-';
        if ($protInfo{$protID}{'MasterProt'} && $masterProts{$protInfo{$protID}{'MasterProt'}}{'gene'}[0]) {
            $geneStrg = shift @{$masterProts{$protInfo{$protID}{'MasterProt'}}{'gene'}};
        }
        my $uniAccStrg = '-';
        if ($protInfo{$protID}{'MasterProt'} && $masterProts{$protInfo{$protID}{'MasterProt'}}{'uniACC'}[0]) {
            $uniAccStrg = shift @{$masterProts{$protInfo{$protID}{'MasterProt'}}{'uniACC'}};
        }
        my $pval = $gseaData{$protID}{'pvalue'};
        $pval = '-' unless ($pval && $pval ne "NA");

        $worksheet->write_string($row, 0, $geneSetsStr, $contentFormat);
        $worksheet->write_string($row, 1, $protInfo{$protID}{'Name'}, $contentFormat);
        $worksheet->write_string($row, 2, $uniAccStrg, $contentFormat);
        $worksheet->write_string($row, 3, $geneStrg, $contentFormat);
        $worksheet->write_string($row, 4, $gseaData{$protID}{'value'}, $contentFormat);
        if ($pval eq "-") {
            $worksheet->write_string($row, 5, $pval, $contentFormat);
        } else {
            $worksheet->write_number($row, 5, $pval, $contentFormat);
        }
        $worksheet->write_number($row, 6, $gseaData{$protID}{'pepNb'}, $contentFormat);
        
        my $protDesc = $protInfo{$protID}{'Des'} . " - " . $protInfo{$protID}{'Org'};
        $maxCharDesc = length($protDesc) if (length($protDesc) > $maxCharDesc);
        $worksheet->write_string($row, 7, $protDesc, $contentFormat);

        my $rowHeight = ($geneSetCount > 1)? $geneSetCount * 12 : 14;
        $worksheet->set_row($row, $rowHeight);

        $row++;
    }
    
    $worksheet->set_column(0, 0, $maxCharGS * 1.2);
    $worksheet->set_column(1, 1, 15);
    $worksheet->set_column(2, 2, 15);
    $worksheet->set_column(3, 3, 15);
    $worksheet->set_column(4, 4, 15);
    $worksheet->set_column(5, 5, 15);
    $worksheet->set_column(6, 6, 15);
    $worksheet->set_column(7, 7, $maxCharDesc);
    $workbook->close();
}

####>Revision history<####
# 1.1.1 [ENHANCEMENT] Small improvement to protein lists display (VL 10/06/21)
# 1.1.0 [ENHENCEMENT] Add interactive dotplot and heatmap (VL 27/05/21)
# 1.0.0 Created (VL 09/11/20)
