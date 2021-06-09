#!/usr/local/bin/perl -w

################################################################################
# startGSEA.cgi       1.0.0                                                    #
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
use POSIX qw(strftime ceil); # to get the time and ceil function
use File::Path qw(rmtree); # remove_tree
use File::Copy;


my %promsPath = &promsConfig::getServerInfo;
my %clusterInfo = &promsConfig::getClusterInfo;
my $dbh = &promsConfig::dbConnect;
my $userID = $ENV{'REMOTE_USER'};
my ($lightColor, $darkColor) = &promsConfig::getRowColors;

my $quantifID = &promsMod::cleanNumericalParameters(param('quanti'));
my $action = param('action') || "form";
my $tmpDir = param('tmpDir') || strftime("%Y%m%d%H%M%S", localtime);

my $projectID = &promsMod::getProjectID($dbh, $quantifID,'QUANTIFICATION');

if ($action eq 'fetchGeneSets') {
    &ajaxFetchGeneSets;
    $dbh->disconnect;
    exit;
}

###############################################
### Get basic info about the quantification ###
###############################################
my ($experimentID, $designID, $quantiName, $quantifAnnot) = $dbh->selectrow_array("SELECT D.ID_EXPERIMENT, D.ID_DESIGN, Q.NAME, Q.QUANTIF_ANNOT FROM QUANTIFICATION Q INNER JOIN DESIGN D ON Q.ID_DESIGN = D.ID_DESIGN WHERE ID_QUANTIFICATION = $quantifID");

# Check what type of quantification we are dealing with
my $isRatio = 0;
my $ratioPos;
if ($quantifAnnot =~ /RATIO_TYPE=([^:]+)/) {
    my $ratioType = $1;
    if ($ratioType ne "None") {
        $isRatio = 1;
        $ratioPos = param('subQuantif') // 1;
    }
}

##########################
### Form was submitted ###
##########################
if (param('start')) {
    &startGSEA;
    exit;
}

#########################
### Display HTML Form ###
#########################
if ($action eq "form") {
    ##########################
    ### Fetch info from DB ###
    ##########################
    my (%parameters, %quantifInfo, %quantifValues);
    my $quantif;
    if ($isRatio) {
        $quantif = $quantifID.'_'.$ratioPos;
        %parameters = (
            QUANTIF_FAMILY  => 'RATIO',
            VIEW            => 'log2',
            NUM_PEP_CODE    => 'NUM_PEP_USED',
            QUANTIF_LIST    => [$quantif]
        );

        &promsQuantif::fetchQuantificationData($dbh, \%parameters, \%quantifInfo, \%quantifValues);
    }

    my $maxProtNb = 0;  # To assess cluster memory needs
    foreach my $quantifPos (keys %quantifValues) {
        my $protNb = scalar keys %{$quantifValues{$quantifPos}};
        $maxProtNb = $protNb if ($protNb > $maxProtNb);
    }

    # HTML string for selection of sub quantification on which GSEA will be performed
    my $rank = 0;
    my $subQuantifStrg = "<SELECT name=\"subQuantif\" id=\"subQuantif\">\n";
    if ($isRatio) {
        foreach my $ratio (@{$quantifInfo{$quantifID}[1]{RATIOS}}) {
            $rank++;
            my $selected = ($rank == $ratioPos)? 'selected' : '';
            my ($testCond, $refCond) = split(/\//, $ratio);
            my $supRatioTag = '';
            if ($testCond =~ /%/) {
                $testCond =~ s/%.+//;
                $refCond =~ s/%.+//;
                $supRatioTag = '°';
            }
            $subQuantifStrg .= "<OPTION value=\"$rank\" $selected>$quantifInfo{$quantifID}[2]{$testCond}{NAME}$supRatioTag/$quantifInfo{$quantifID}[2]{$refCond}{NAME}$supRatioTag</OPTION>\n";
        }
    }
    $subQuantifStrg .= "</SELECT>\n";

    # HTML string for species (deduced from the proteins in project)
    my $speciesStrg = "<OPTION value=\"\" selected>-= Select =-</OPTION>\n";
    $speciesStrg .= "<OPTION value=\"Custom\">Not relevant (only for custom gene sets)</OPTION>\n";

    my $sthSpeciesSN = $dbh->prepare("SELECT GROUP_CONCAT(DISTINCT ORGANISM SEPARATOR \"','\") FROM PROTEIN WHERE ID_PROJECT = $projectID");
    $sthSpeciesSN->execute;
    my ($speciesNameStrg) = $sthSpeciesSN->fetchrow_array;
    $sthSpeciesSN->finish;
    $speciesNameStrg = "'" . $speciesNameStrg . "'";
    my $sthSpecies = $dbh->prepare("SELECT ID_SPECIES, SCIENTIFIC_NAME FROM SPECIES WHERE SCIENTIFIC_NAME IN ($speciesNameStrg)");
    $sthSpecies->execute;
    while (my ($speciesID, $speciesSN) = $sthSpecies->fetchrow_array) {
        $speciesStrg .= "<OPTION value=\"$speciesID\">$speciesSN</OPTION>\n";
    }
    $sthSpecies->finish;

    # HTML string for protein lists (exclude/restrict)
    my $protLists = "<OPTION value=\"\">-= Select =-</OPTION>\n";
    my %categoryList;
    my $sthCat = $dbh->prepare(
        "SELECT CL.ID_CLASSIFICATION, CL.NAME, ID_CATEGORY, CA.NAME, CA.LIST_TYPE
        FROM CLASSIFICATION CL
        INNER JOIN CATEGORY CA ON CL.ID_CLASSIFICATION = CA.ID_CLASSIFICATION
        WHERE ID_PROJECT = $projectID
        ORDER BY CL.NAME ASC, CA.DISPLAY_POS ASC"
    );
    $sthCat->execute;
    while (my ($classID, $className, $catID, $catName, $listType) = $sthCat->fetchrow_array) {
        $listType = 'PROT' unless $listType;
        next unless $listType eq 'PROT'; # restrict to protein lists
        $categoryList{$classID}{'CL_NAME'} = $className unless $categoryList{$classID};
        push @{$categoryList{$classID}{'CAT'}}, [$catID, $catName, $listType];
    }
    $sthCat->finish;

    foreach my $classID (sort{lc($categoryList{$a}{'CL_NAME'}) cmp lc($categoryList{$b}{'CL_NAME'})} keys %categoryList) {
        $protLists .= "<OPTGROUP label=\"$categoryList{$classID}{CL_NAME}:\">\n";
        foreach my $refCat (@{$categoryList{$classID}{'CAT'}}) {
            my $typeStrg = ($refCat->[2] eq 'SITE')? ' [Sites]' : '';  # Not used for GSEA yet
            $protLists .= "<OPTION value=\"$refCat->[0]\">$refCat->[1]$typeStrg</OPTION>\n";
        }
        $protLists .= "</OPTGROUP>\n";
    }

    # Set parameters for protein filter (on peptides)
    my $hasTruePepCount=($quantifInfo{$quantifID}[1]{NUM_TRUE_USED})? $quantifInfo{$quantifID}[1]{NUM_TRUE_USED}[0] : 0;
    my $disabTruePep = ($hasTruePepCount)? ' ' : ' disabled';


    ##################
    ### Start HTML ###
    ##################
    print header(-'content-encoding' => 'no',-charset=>'utf-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<TITLE>Gene Set Enrichment Analysis On Quantification Data</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT type="text/javascript">
// Dictionnary of help popup texts (avoid using '"')
const helpText = {
    gsea: 'GSEA searches for gene or protein sets (e.g. proteins belonging to the same pathway or associated with the same GO term) that are all enriched at the same time in a given condition.<br>The idea is that a small but consistent increase in the expression of proteins from the same gene set will most probably have an impact on phenotype.',
    geneSets: 'From which database should we take the gene sets ?<br>Numbers in parentheses are the number of gene sets in each databank and the date is the version of the file used.<br>You can get in contact with the lab if there is not any Gene Sets databank available for your species of interest yet or if you want to add your own.',
    protSelection: 'Restrict GSEA to a subset of proteins.',
    weights: 'Whether to weight the quantification scores (e.g. ratios) with the number of peptides of each protein (used as a measure of confidence).<br>Proteins with infinite ratios are discarded if this is not used (impossible to sort the protein scores).',
    pval: 'The p-value cutoff for the Gene Set Enrichment. Only the gene sets enriched with an adjusted p-value under this threshold will be reported. It does not always mean that these gene sets are the most enriched, only that we are more confident that they are actually enriched.',
    padj: 'The method to adjust p-values for multiple testing. Change this only if you know what you are doing.<br>* Recommended option',
    maxGS: 'The maximum number of enriched Gene Sets that you want to analyse. The GSEA plots will not be drawn for more gene sets than this number.<br>Cannot equal 0, it is reset to 25 if you set it to this value.<br>Set it to -1 to consider all possible gene sets (that pass the p-value cutoff).'
    // simplify: 'Gene Sets are often redondant (see for example GO terms that are embeded in one another). This function reduces this redundancy using semantic similarity among terms and remove those highly similar terms by keeping one representative term.'
};
|;
    &promsMod::popupInfo();
    print qq
|function checkForm(myForm){
    if (!myForm.gseaName.value) {
        alert('Enter a name for the analysis.');
        return false;
    } else if (!myForm.species.value) {
        alert('Select a species.');
        return false;
    } else if (!myForm.geneSets.value) {
        alert('Select a Gene Set database.');
        return false;
    }
    return true;
}

function getXMLHTTP(){
    var xhr = null;
    if(window.XMLHttpRequest) {  // Firefox & others
        xhr = new XMLHttpRequest();
    } else if (window.ActiveXObject) {  // Internet Explorer
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

var XHR = null;
function updateSpecies(species) {
    var selGeneSets = document.getElementById('geneSets');

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
    var GETstring = 'action=fetchGeneSets&ID=' + $projectID + '&species=' + species;
    var URL = "$promsPath{cgi}/startGSEA.cgi?" + GETstring;
    XHR.open("GET", URL, true);
    XHR.onreadystatechange = function() {
        if (XHR.readyState == 4 && XHR.responseText){
            selGeneSets.innerHTML = XHR.responseText;
        }
    }
    XHR.send(null);
}
    </SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Gene Set Enrichment Analysis On Quantification Data</FONT>
&nbsp;
<FONT class="help" onmouseover="popup(helpText.gsea);" onmouseout="popout();">[?]</FONT>
<BR>

<FORM name="gseaForm" id="gseaForm" method="post" enctype="multipart/form-data" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="quanti" value="$quantifID"/>
<INPUT type="hidden" name="protNb" value="$maxProtNb"/>
<BR>
<TABLE bgcolor=$darkColor border=0>
    <TR>
        <TH valign='top' align='right'>&nbsp;&nbsp;Name : </TH>
        <TD bgcolor=$lightColor><INPUT type='text' name='gseaName' size="30"></TD>
    </TR>
    <TR>
        <TH valign='top' align='right'>&nbsp;&nbsp;Description : </TH>
        <TD bgcolor=$lightColor><TEXTAREA rows="2" cols="60" name="description"></TEXTAREA></TD>
    </TR>
    <TR>
        <TH valign='top' align='right'>&nbsp;&nbsp;Quantification : </TH>
        <TD bgcolor=$lightColor><B>$quantiName</B></TD>
    </TR>
    <TR>
        <TH valign='top' align='right'>&nbsp;&nbsp;Test : </TH>
        <TD bgcolor=$lightColor>$subQuantifStrg</TD>
    </TR>
    <TR>
        <TH align=right nowrap>
            Protein lists<SUP onmouseover="popup(helpText.protSelection)" onmouseout="popout()">?</SUP> :
        </TH>
        <TD bgcolor=$lightColor nowrap>
            <SELECT name="protSelType">
                <OPTION value="exclude">Exclude</OPTION>
                <OPTION value="restrict">Restrict to</OPTION>
            </SELECT>
            &nbsp;proteins from List:
            <SELECT name="protSelection">
                $protLists
            </SELECT>
        </TD>
    </TR>
    <TR>
        <TH></TH>
        <TD valign='top' align='center'><FONT class="title3">Gene Sets parameters</FONT></TD>
    </TR>
    <TR>
        <TH valign='top' align='right'>&nbsp;&nbsp;Species : </TH>
        <TD bgcolor=$lightColor>
            <SELECT name="species" id="species"  style="width:400px" onchange="updateSpecies(this.value);">
                $speciesStrg
            </SELECT>
        </TD>
    </TR>
    <TR>
        <TH valign='top' align='right'>&nbsp;&nbsp;Gene set<SUP onmouseover="popup(helpText.geneSets);" onmouseout="popout();">?</SUP> : </TH>
        <TD bgcolor=$lightColor>
            <SELECT name="geneSets" id="geneSets" style="width:400px">
                <OPTION value="" selected>-= Select a species first =-</OPTION>
            </SELECT>
        </TD>
    </TR>
    <TR>
        <TH></TH>
        <TD valign='top' align='center'><FONT class="title3">GSEA parameters</FONT></TD>
    </TR>
    <TR>
        <TH valign='top' align='right'>&nbsp;&nbsp;Peptides to weight protein scores<SUP onmouseover="popup(helpText.weights);" onmouseout="popout();">?</SUP> : </TH>
        <TD bgcolor=$lightColor>
            <SELECT name="pepWeights" id="pepWeights" style="width:180px">
                <OPTION value="none" selected>None</OPTION>
                <OPTION value="all">All</OPTION>
                <OPTION value="distinct">Distinct</OPTION>
                <OPTION value="msms"$disabTruePep>Truly identified</OPTION>
            </SELECT>
        </TD>
    </TR>
    <TR>
        <TH valign='top' align='right'>&nbsp;&nbsp;P-value cutoff<SUP onmouseover="popup(helpText.pval);" onmouseout="popout();">?</SUP> : </TH>
        <TD bgcolor=$lightColor>
            <INPUT type="number" name="pvalCutoff" id="pvalCutoff" min="0" max="1" step="any" value="0.05" style="width:180px">
        </TD>
    </TR>
    <TR>
        <TH valign='top' align='right'>&nbsp;&nbsp;P-value adjustment<SUP onmouseover="popup(helpText.padj);" onmouseout="popout();">?</SUP> : </TH>
        <TD bgcolor=$lightColor>
            <SELECT name="pvalCorrection" id="pvalCorrection" style="width:180px">
                <OPTION value="none">None (Not recommended)</OPTION>
                <OPTION value="BH" selected>Benjamini-Hochberg*</OPTION>
                <OPTION value="BY">Benjamini-Yekutieli</OPTION>
                <OPTION value="hochberg">Hochberg</OPTION>
                <OPTION value="bonferroni">Bonferroni</OPTION>
                <OPTION value="holm">Holm</OPTION>
                <OPTION value="hommel">Hommel</OPTION>
            </SELECT>
        </TD>
    </TR>
    <TR>
        <TH valign='top' align='right'>&nbsp;&nbsp;Max enriched GS to report<SUP onmouseover="popup(helpText.maxGS);" onmouseout="popout();">?</SUP> : </TH>
        <TD bgcolor=$lightColor>
            <INPUT type="number" name="maxRelevantGS" id="maxRelevantGS" min="-1" step="1" value="25" style="width:180px">
        </TD>
    </TR>
    <TR>
        <TD align="center" colspan=2>
            <INPUT type="submit" name="start" value="Start" />
            &nbsp;&nbsp;&nbsp;&nbsp;
            <INPUT type="reset" value="Clear" />
            <INPUT type="hidden" name="tmpDir" value="$tmpDir" />
        </TD>
    </TR>
</TABLE>
</FORM>
</CENTER>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">
    setPopup();
</SCRIPT>
</BODY>
</HTML>
|;

}
$dbh->disconnect;


sub startGSEA {  # Globals: $promsPath, $quantifID, $isRatio
    ##################
    ### Start HTML ###
    ##################
    print header(-'content-encoding' => 'no', -charset => 'utf-8');
    warningsToBrowser(1);
    print qq
|<HTML>
    <HEAD>
        <TITLE>Launching GSEA</TITLE>
        <LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
    </HEAD>
    <BODY background="$promsPath{images}/bgProMS.gif">
        <CENTER>
        <FONT class="title1">Launching Gene Set Enrichment Analysis</FONT>
        </CENTER>
        <BR><BR><BR>
        <FONT class="title3">Launching GSEA process
|;

    # GSEA DB parameters
    my ($gseaName, $gseaDes, $subQuantif, $tmpDir);
    $gseaName   = param('gseaName');
    $gseaDes    = param('description') || '';
    $subQuantif = param('subQuantif');
    $tmpDir     = param('tmpDir');
    # Protein selection parameters
    my ($protSelType, $protSelection);
    my $protSelType   = param('protSelType');
    my $protSelection = param('protSelection') || 0;
    # Gene Sets parameters
    my ($speciesID, $geneSetsID);
    $geneSetsID = param('geneSets');
    $speciesID  = param('species');
    # GSEA analysis parameters
    my ($pepWeights, $pvalCutoff, $pvalCorr, $maxRelGS);
    $pepWeights = param('pepWeights');
    $pvalCutoff = param('pvalCutoff');
    $pvalCorr   = param('pvalCorrection');
    $maxRelGS   = param('maxRelevantGS') || 25;

    my ($gseaID, $gseaDir, $geneSetsVersion, $gseaParamStrg);
    $geneSetsVersion = $dbh->selectrow_array("SELECT VERSION_DATE FROM GENESETS WHERE ID_GENESETS = $geneSetsID");
    $gseaDir = "$promsPath{'tmp'}/gsea/$tmpDir";

    unless (-d "$gseaDir") {
        # Clean older dirs if any
        &promsMod::cleanDirectory("$promsPath{'tmp'}/gsea",'7d');  # 7 days to avoid removing errors too quickly
        # Create tmp directory to store results of GSEA
        mkdir "$promsPath{'tmp'}/gsea" unless (-d "$promsPath{'tmp'}/gsea");
        mkdir "$gseaDir";
    }

    $gseaParamStrg  = "quantification=#$quantifID";
    $gseaParamStrg .= ";targetPos=$subQuantif";
    $gseaParamStrg .= ";isRatio=$isRatio";
    $gseaParamStrg .= ($speciesID eq "Custom")? ";species=$speciesID" : ";species=#$speciesID";
    $gseaParamStrg .= ";geneSetsDB=#$geneSetsID";
    $gseaParamStrg .= ";geneSetsVersion=$geneSetsVersion";
    $gseaParamStrg .= ";protSel=$protSelType,#$protSelection" if ($protSelection);
    $gseaParamStrg .= ";weightByPep=$pepWeights";
    $gseaParamStrg .= ";pvalCutoff=$pvalCutoff";
    $gseaParamStrg .= ";pvalCorr=$pvalCorr";
    $gseaParamStrg .= ";maxGS=$maxRelGS";
    # $gseaParamStrg .= ";simplify=$simplify";

    ##############################
    ### Insert GSEA info in DB ###
    ##############################
    my $sthInsGSEA = $dbh->prepare("INSERT INTO PATHWAY_ANALYSIS (ID_PARENT_PATHWAYANA, ID_EXPERIMENT, ID_CATEGORY, ANALYSIS_TYPE, NAME, DES, PARAM_STRG, STATUS, RECORD_DATE, UPDATE_DATE, UPDATE_USER) VALUES (NULL, $experimentID, NULL, \"GSEA\", \"$gseaName\", \"$gseaDes\", \"$gseaParamStrg\", 0, NOW(), NOW(), \"$userID\")");
    $sthInsGSEA->execute;
    $sthInsGSEA->finish;
    $gseaID = $dbh->last_insert_id(undef, undef,'PATHWAY_ANALYSIS','ID_PATHWAY_ANALYSIS');

    #############################################################
    ### Insert link between GSEA and quantif/subquantif in DB ###
    #############################################################
    my $sthInsGseaQuanti = $dbh->prepare("INSERT INTO PATHWAYANA_QUANTIFICATION (ID_QUANTIFICATION, ID_PATHWAY_ANALYSIS, RATIO) VALUES ($quantifID, $gseaID, \"$subQuantif\")");  # RATIO used as target_pos
    $sthInsGseaQuanti->execute;
    $sthInsGseaQuanti->finish;

    #############################
    ### Insert GSEA job in DB ###
    #############################
    my $gseaDBInfo = "ID_GSEA=$gseaID;TYPE=GSEA;GSEA_NAME=$gseaName;ID_PARENT_QUANTIF=$quantifID";
    my $logFilePath = "$gseaDir/status_$gseaID.out";
    my $errorFilePath = "$gseaDir/error_$gseaID.txt";
    $dbh->do("INSERT INTO JOB_HISTORY (ID_JOB, ID_USER, ID_PROJECT, TYPE, JOB_STATUS, FEATURES, SRC_PATH, LOG_PATH, ERROR_PATH, STARTED) VALUES('$tmpDir', '$userID', $projectID, 'Functional Analysis [GSEA]', 'Queued', '$gseaDBInfo', '$gseaDir', '$logFilePath', '$errorFilePath', NOW())");
    $dbh->commit;
    $dbh->disconnect;

    ##################################################
    ### Fork to launch GSEA computation on cluster ###
    ##################################################
    my $childPid = fork;
    unless ($childPid) {  # child here
        # Disconnecting from server
        open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
        open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
        
        my $cgiUnixDir = `pwd`;
        $cgiUnixDir =~ s/\/*\s*$//;

        if ($clusterInfo{'on'}) {
            # cd required for script to find myproms .pm files
            my $commandString = "cd $cgiUnixDir; $clusterInfo{path}{perl}/perl runGSEA.pl $gseaID $userID $tmpDir 2> $errorFilePath";
            my $maxMem = ceil(log(param('protNb'))/log(10));  # Adjust max cluster memory to ceiling of log10(protNb)
            my %jobParameters = (
                maxMem          => $maxMem . "Gb",
                numCPUs         => 1,
                maxHours        => 1,
                jobName         => "myProMS_GSEA_$tmpDir",
                pbsRunDir       => $gseaDir,
                commandBefore   => "echo \"Launched GSEA $gseaID on cluster\" >> $logFilePath",
                noWatch         => '1',
            );
            my ($pbsError, $pbsErrorFile, $jobClusterID) = $clusterInfo{'runJob'}->($gseaDir, $commandString, \%jobParameters);

            # Add cluster job id to current job in DB
            $dbh = &promsConfig::dbConnect;
            $dbh->do("UPDATE JOB_HISTORY SET ID_JOB_CLUSTER = 'C$jobClusterID' WHERE ID_JOB = '$tmpDir'");
            $dbh->commit;
            $dbh->disconnect;

            if ($pbsError) {  # move PBS error message to job error file
                system "cat $pbsErrorFile >> $errorFilePath";
                $dbh = &promsConfig::dbConnect;
                $dbh->do("UPDATE PATHWAY_ANALYSIS SET STATUS = -1, UPDATE_DATE = NOW() WHERE ID_PATHWAY_ANALYSIS = $gseaID"); # Failed
                $dbh->commit;
                $dbh->disconnect;
                die "An error occured during cluster computation, check $errorFilePath for more information.";
            } elsif (-d "$promsPath{'gsea'}/project_$projectID/gsea_$gseaID") {  # Move cluster files to final directory
                system "mv $gseaDir/*.txt $promsPath{'gsea'}/project_$projectID/gsea_$gseaID/.";
                system "mv $gseaDir/*.sh $promsPath{'gsea'}/project_$projectID/gsea_$gseaID/.";
            }
        } else {  # Local launch
            system "echo \"Launched GSEA $gseaID on local server\" >> $logFilePath; cd $cgiUnixDir; ./runGSEA.pl $gseaID $userID $tmpDir 2> $errorFilePath; echo \"Computation of GSEA $gseaID on local server ended\" >> $logFilePath";
        }
        exit;
    }

    print qq
|       <BR>
        Done.</FONT>
        <BR>
        <FONT class="title3">This page will refresh itself in a few seconds.</FONT>
|;
    # Calling watch popup window
    sleep 3;
    print qq
|<SCRIPT type="text/javascript">
    var monitorJobsWin = window.open("$promsPath{cgi}/monitorJobs.cgi?filterType=Functional Analysis [GSEA]&filterDateNumber=1&filterDateType=DAY&filterStatus=Queued&filterStatus=Running&filterProject=$projectID",'monitorJobsWindow','width=1200,height=500,scrollbars=yes,resizable=yes');
monitorJobsWin.focus();
</SCRIPT>
|;
    sleep 5;
    print qq
|<SCRIPT type="text/javascript">
    parent.itemFrame.location = "$promsPath{cgi}/openProject.cgi?ACT=experiment&VIEW=functAna&branchID=gseanalysis:$gseaID";
</SCRIPT>
</CENTER>
</BODY>
</HTML>
|;
}


sub ajaxFetchGeneSets {  # Globals: $dbh, $projectID
    my $projectID = param('ID');
    my $speciesID = param('species');
    my $geneSetsStrg = "<OPTION value=\"\" selected>-= Select =-</OPTION>\n";
    my %geneSets;

    if (!$speciesID) {
        $geneSetsStrg = "<OPTION value=\"\" selected>-= Select a species first =-</OPTION>\n";
    } else {
        my $speciesSthStrg = ($speciesID eq "Custom")? "IS NULL" : "= $speciesID";
        my $sthGeneSets = $dbh->prepare("SELECT ID_GENESETS, NAME, NUM_SETS, VERSION_DATE FROM GENESETS WHERE ID_SPECIES $speciesSthStrg AND (ID_PROJECT IS NULL OR ID_PROJECT = $projectID) AND GENESETS_STATUS = 1");
        $sthGeneSets->execute;
        while (my ($geneSetsID, $geneSetsName, $geneSetsNb, $geneSetsVersion) = $sthGeneSets->fetchrow_array) {
            $geneSets{$geneSetsID} = "$geneSetsName ($geneSetsNb) - $geneSetsVersion";
        }
        $sthGeneSets->finish;

        if (scalar %geneSets) {
            foreach my $geneSetsID (sort{&promsMod::sortSmart($geneSets{$a}, $geneSets{$b})} keys %geneSets) {
                $geneSetsStrg .= "<OPTION value=\"$geneSetsID\">$geneSets{$geneSetsID}</OPTION>\n";
            }
        } else {
            $geneSetsStrg = "<OPTION value=\"\" selected>-= No Gene Sets available for this species (yet) =-</OPTION>\n";
        }
    }

    print header(-type=>'text/plain',-charset=>'UTF-8');
    warningsToBrowser(1);
    print "$geneSetsStrg";
}

####>Revision history<####
# 1.0.0 Created (VL 21/10/20)
