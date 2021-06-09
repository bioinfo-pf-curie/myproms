#!/usr/local/bin/perl -w

################################################################################
# displayHydro.cgi       1.0.0                                                 #
# Authors: V. Laigle (Institut Curie)                                          #
# Contact: myproms@curie.fr                                                    #
# Display the results of Hydrophobicity computations                           #
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
use File::Path qw(rmtree);
use Spreadsheet::WriteExcel;

#####################
### Configuration ###
#####################
my %promsPath = &promsConfig::getServerInfo;
my $dbh = &promsConfig::dbConnect;
my $userID = $ENV{'REMOTE_USER'};
my ($lightColor, $darkColor) = &promsConfig::getRowColors;

my $action = (param('ACT'))? param('ACT') : 'summary';
my $quantifID = &promsMod::cleanNumericalParameters(param('ID'));
my ($fromItem, $fromID) = split(':', param('FROM'));
my $view = (param('VIEW'))? lc(param('VIEW')) : 'graph';
my $sortKey = param('SORT') || 'HI';
my $sortOrder = param('SORT_ORDER') || 'asc';
my $pct_cutoff = (defined(param('PCT_CUTOFF')))? param('PCT_CUTOFF') : 99;
my $measure = param('MEAS') || 'shift';

my $projectID = &promsMod::getProjectID($dbh, $quantifID, 'QUANTIFICATION');
my @userInfo = &promsMod::getUserInfo($dbh, $userID, $projectID);
my $projectAccess = ${$userInfo[2]}{$projectID};
my $disabSave = ($projectAccess eq 'guest')? ' disabled' : '';

my $quantifDir        = "$promsPath{'quantification'}/project_$projectID/quanti_$quantifID";
my $hydroInfoFile     = "$quantifDir/hydro_info.txt";
my $hydroResultsFile  = "$quantifDir/hydro_results.txt";
my $linModFile        = "$quantifDir/lin_mod_results.txt";
my $unmodifPropsFile  = "$quantifDir/unmodified_proportions.txt";
my $hydroGraph        = "$quantifDir/hydro_graph";
my $hydroGraphDisplay = "$promsPath{quanti_html}/project_$projectID/quanti_$quantifID/hydro_graph";
my (%hydroInfo, %hydroResults, %linMods, %unmodifProps);

my ($name, $quantifAnnot, $status, $upDate, $user) = $dbh->selectrow_array(
    "SELECT NAME, QUANTIF_ANNOT, STATUS, UPDATE_DATE, UPDATE_USER
    FROM QUANTIFICATION WHERE ID_QUANTIFICATION = $quantifID"
);

#####################
### Start Display ###
#####################
if ($action eq 'summary') {
    &getHydroInfo(\%hydroInfo);
    $dbh->disconnect;
    &displaySummary;
} elsif ($action eq 'edit') {
    $dbh->disconnect;
    if (param('submit')) {
        &editHydro;
    } else {
        &displayFormEdit;
    }
} elsif ($action eq 'delete') {
    $dbh->disconnect;
    &deleteHydro;
} elsif ($action eq 'display') {
    &getHydroInfo(\%hydroInfo);
    &parseHydroResults(\%hydroResults, $hydroResultsFile) if (-e $hydroResultsFile && $view eq 'table');
    &parseLinMods(\%linMods, $linModFile) if (-e $linModFile && $view eq 'linmod');
    # &parseUnmodifProps(\%unmodifProps, $unmodifPropsFile) if (-e $unmodifPropsFile);
    $dbh->disconnect;
    &displayVisualization;
} elsif ($action eq 'export') {
    &getHydroInfo(\%hydroInfo);
    &parseHydroResults(\%hydroResults, $hydroResultsFile) if (-e $hydroResultsFile);
    &parseLinMods(\%linMods, $linModFile) if (-e $linModFile);
    &parseUnmodifProps(\%unmodifProps, $unmodifPropsFile) if (-e $unmodifPropsFile);
    $dbh->disconnect;
    &exportHydroResults;
}
exit;


sub displayVisualization {  # Globals: %promsPath, $quantifID, $projectID, $name, $view, $fromItem, $fromID
    print header(-'content-encoding' => 'no',-charset=>'utf-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
    <TITLE>Display Hydrophobicity vs. Retention Time</TITLE>
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
    <SCRIPT src="$promsPath{html}/js/local/genericPlot.js"></SCRIPT>
    <SCRIPT src="$promsPath{html}/js/local/vennDiagram.js"></SCRIPT>
    <SCRIPT src="$promsPath{html}/js/local/peptidePlot.js"></SCRIPT>
    <SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
    <SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
    <SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
    <SCRIPT LANGUAGE="JavaScript">
// Dictionnary of help popup texts (avoid using '"')
const helpText = {
    HI: 'The Hydrophobicity Index, expressed in % of the Acetonitrile gradient. It is the % of ACN at which the peptide should exit the chromatography. If the ACN gradient is linear, HI is supposed to be correlated with the retention time.',
    rt_shift: 'Difference between observed retention time and retention time predicted through SSRCalc algorithm. Parentheses indicate that in the corresponding analysis, the peptide shift in RT is not over the specified percentile cutoff (but that is still true for at least one analysis). Values in grey indicate peptides reextracted with MBR.',
    rt_obs: 'Observed retention time. Parentheses indicate that in the corresponding analysis, the peptide shift in RT is not over the specified percentile cutoff (but that is still true for at least one analysis). Values in grey indicate peptides reextracted with MBR.',
    rt_pred: 'Retention time predicted through SSRCalc algorithm. Parentheses indicate that in the corresponding analysis, the peptide shift in RT is not over the specified percentile cutoff (but that is still true for at least one analysis). Values in grey indicate peptides reextracted with MBR.',
    pct_cutoff: 'This cutoff applies on the shift between observed and predicted retention times. All peptides for which the shift is above the specified percentile in at least one analysis will be displayed (and only these ones).<br>CAUTION: Make sure not to set a value too low if it concerns a lot of peptides, otherwise it might make the application quite laggy.'
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
function sortTable(sortItem) {
    var sortOrder = (sortItem != "$sortKey")? "$sortOrder" : ("$sortOrder" == "asc")? "desc" : "asc";
    window.location="./displayHydro.cgi?ID=$quantifID&ACT=$action&VIEW=$view&FROM=$fromItem:$fromID&PCT_CUTOFF=$pct_cutoff&MEAS=$measure&SORT=" + sortItem + "&SORT_ORDER=" + sortOrder;
}
   </SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<DIV id="displayDIV" class="popup">
    <DIV id="infoDIV"></DIV>
</DIV>
|;
    print &displayBackButton;
    print qq
|<CENTER>
<FONT class="title">Display Hydrophobicity vs. Retention Time for <FONT color="red">$name</FONT></FONT>
&nbsp;&nbsp;
<INPUT type="button" class="title3" value="Export" onclick="window.location='$promsPath{cgi}/displayHydro.cgi?ACT=export&ID=$quantifID&FROM=$fromItem:$fromID';">
<BR><BR>
<LABEL for='dataVizType' style='font-weight:bold'>Visualization type: </LABEL>
<SELECT id='dataVizType' onchange="window.location='$promsPath{cgi}/displayHydro.cgi?ACT=display&VIEW=' + this.value + '&ID=$quantifID&FROM=$fromItem:$fromID';">
|;

    my @vizTypes = ('Graph', 'Table', 'Linear Models');
    foreach my $vizType (@vizTypes) {
        my $vizTypeCode = ($vizType eq 'Linear Models')? 'linmod' : $vizType;
        my $disabStrg = (lc($vizTypeCode) eq $view)? " selected disabled" : " ";
        print "<OPTION value=\"$vizTypeCode\"$disabStrg>$vizType</OPTION>\n";
    }
    print "</SELECT>\n";

    if ($view eq 'graph') {
        &displayGraphHydroRT;
    } elsif ($view eq 'linmod') {
        &displayLinearModels;
    } else {
        &displayTableHydro;
    }

    print qq
|<BR><BR><BR><BR><BR><BR>
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

sub displayTableHydro {  # Globals: %hydroInfo, %hydroResults, $lightColor, $darkColor, $sortKey, $sortOrder
    
    my $numAna = scalar keys %{$hydroInfo{'ana'}};
    my $numSamp = scalar keys %{$hydroInfo{'samp'}};
    my $numRowHead = ($numSamp > 1)? 3 : 2;
    my (%sampOrder, %anaOrder);
    my @measures = ('shift', 'obsRT', 'predRT');
    my @measNames = ('Shift obs. <-> pred.', 'Observed RT', 'Predicted RT');
    my ($tableTitle, $helpKey) = ($measure eq 'shift')? ('Shift between observed and predicted Retention Time (min)', 'rt_shift') : ($measure eq 'obsRT')? ('Observed Retention Time (min)', 'rt_obs') : ('Predicted Retention Time (min)', 'rt_pred');
    $tableTitle = "<SUP onmouseover=\"popup(helpText.$helpKey)\" onmouseout=\"popout()\">?</SUP>" . $tableTitle;

    my ($seqColor, $hiColor, %rtColor);
    $seqColor = "black";
    $hiColor = "black";
    foreach my $anaID (keys %{$hydroInfo{'ana'}}) {
        $rtColor{$anaID}  = "black";
    }
    if ($sortKey eq 'HI') {
        $hiColor = "red";
    } elsif ($sortKey =~ /\d+/) {  # Sort by a given analysis (ID)
        $rtColor{$sortKey} = "red";
    } else {  # $sortKey eq "seq"
        $seqColor = "red";
    }

    sub sortTable {
        if ($sortKey eq 'HI') {  # Take HI value of any analysis for this sequence
            return ($sortOrder eq 'desc')? 
                ($hydroResults{$b}{(keys %{$hydroResults{$b}})[0]}->[2] <=> $hydroResults{$a}{(keys %{$hydroResults{$a}})[0]}->[2]) || ($a cmp $b) :
                ($hydroResults{$a}{(keys %{$hydroResults{$a}})[0]}->[2] <=> $hydroResults{$b}{(keys %{$hydroResults{$b}})[0]}->[2]) || ($a cmp $b);
        } elsif ($sortKey =~ /\d+/) {  # Check if value is defined for this analysis before sorting (undefs go last)
            my $anaID = $sortKey;
            my $idx = ($measure =~ /obsRT/)? 1 : ($measure =~ /predRT/)? 4 : 5;  # Index of each value in array
            if (!defined($hydroResults{$a}{$anaID}) && !defined($hydroResults{$b}{$anaID})) {
                return ($a cmp $b);
            } elsif (!defined($hydroResults{$a}{$anaID})) {
                return +1;
            } elsif (!defined($hydroResults{$b}{$anaID})) {
                return -1;
            } else {
                return ($sortOrder eq 'desc')? 
                ($hydroResults{$b}{$anaID}->[$idx] <=> $hydroResults{$a}{$anaID}->[$idx]) || ($a cmp $b) :
                ($hydroResults{$a}{$anaID}->[$idx] <=> $hydroResults{$b}{$anaID}->[$idx]) || ($a cmp $b);
            }
        } else {  # $sortKey eq "seq"
            return ($sortOrder eq 'desc')? ($b cmp $a) : ($a cmp $b);
        }
    }

    print qq
|   &nbsp;&nbsp;&nbsp;&nbsp;
    <LABEL for="measure" style='font-weight:bold'>Measure: </LABEL>
    <SELECT id="measure" onchange="window.location='$promsPath{cgi}/displayHydro.cgi?ACT=$action&VIEW=$view&ID=$quantifID&FROM=$fromItem:$fromID&PCT_CUTOFF=$pct_cutoff&SORT=$sortKey&SORT_ORDER=$sortOrder&MEAS=' + this.value;">
|;
    for my $idx (0..$#measures) {
        my $selected = ($measures[$idx] eq $measure)? ' selected' : '';
        print "<OPTION value=\"$measures[$idx]\"$selected>$measNames[$idx]</OPTION>"
    }
    print qq
|   </SELECT>
    &nbsp;&nbsp;&nbsp;&nbsp;
    <LABEL for="pct_cutoff" style='font-weight:bold'><SUP onmouseover="popup(helpText.pct_cutoff)" onmouseout="popout()">?</SUP>Percentile cutoff :&nbsp;</LABEL>
    <INPUT type="number" id="pct_cutoff" name="pct_cutoff" value="$pct_cutoff" min="0" max="100" step="1">
    &nbsp;
    <INPUT type="button" name="Update" value="Update" onclick="window.location='$promsPath{cgi}/displayHydro.cgi?ACT=$action&VIEW=$view&ID=$quantifID&FROM=$fromItem:$fromID&PCT_CUTOFF=' + document.getElementById('pct_cutoff').value + '&SORT=$sortKey&SORT_ORDER=$sortOrder&MEAS='+document.getElementById('measure').value">
    <BR><BR>
|;
    print qq
|<BR><BR><BR>
<TABLE id="resultsTable" border=0 cellspacing=0 align=center>
<TR bgcolor=$darkColor>
    <TH class="rbBorder" id="pepSeqTH" rowspan="$numRowHead" nowrap>
        <A href="javascript:sortTable('seq')"><FONT color="$seqColor">
            &nbsp;&nbsp;Peptide sequence&nbsp;&nbsp;
        </A>
    </TH>
    <TH class="rbBorder" rowspan="$numRowHead" nowrap>
        &nbsp;&nbsp;Charge&nbsp;&nbsp;
    </TH>
    <TH class="rbBorder" rowspan="$numRowHead" nowrap>
        <A href="javascript:sortTable('HI')"><FONT color="$hiColor">
            &nbsp;&nbsp;Hydrophobicity&nbsp;&nbsp;<br>&nbsp;&nbsp;Index&nbsp;&nbsp;<br>&nbsp;&nbsp;(ACN%)&nbsp;&nbsp;
        </A>
    </TH>
    <TH class="rbBorder" colspan="$numAna" nowrap>
        &nbsp;&nbsp;$tableTitle&nbsp;&nbsp;
    </TH>
</TR>
<TR bgcolor=$darkColor>
|;
    my $sampCol = 0;
    if ($numSamp > 1) {  # Print sample names in displayPos order
        foreach my $sampID (sort {$hydroInfo{'samp'}{$a}->[1] <=> $hydroInfo{'samp'}{$b}->[1]} keys %{$hydroInfo{'samp'}}) {
            $sampCol++;
            $sampOrder{$sampID} = $sampCol;
            my $sampColNb = scalar @{$hydroInfo{'samp'}{$sampID}->[2]};
            print qq
|   <TH class="rbBorder" colspan=$sampColNb nowrap>
        &nbsp;&nbsp;$hydroInfo{'samp'}{$sampID}->[0]&nbsp;&nbsp;
    </TH>
|;
        }
        print "</TR>\n<TR bgcolor=$darkColor>\n";
    } else {
        $sampOrder{(keys %{$hydroInfo{'samp'}})[0]} = 1;  # Get the only sample ID in hydroInfo{'samp'}
    }

    my $anaCol = 0;
    foreach my $sampID (sort {$sampOrder{$a} <=> $sampOrder{$b}} keys %sampOrder) {
        foreach my $anaID (sort {$hydroInfo{'ana'}{$a}->[6] <=> $hydroInfo{'ana'}{$b}->[6]} @{$hydroInfo{'samp'}{$sampID}->[2]}) {
            $anaCol++;
            $anaOrder{$anaID} = $anaCol;
            print qq
|   <TH class="rbBorder" nowrap>
        <A href="javascript:sortTable('$anaID')"><FONT color="$rtColor{$anaID}">
            &nbsp;&nbsp;$hydroInfo{'ana'}{$anaID}->[5]&nbsp;&nbsp;
        </A>
    </TH>
|;
        }
    }
    print "</TR>\n";

    my $lineColor = $lightColor;
    my $dispIonCount = 0;
    my $totalIonCount = scalar keys %hydroResults;
    foreach my $seqZ (sort { sortTable } keys %hydroResults) {
        my $displayPeptide = 0;  # Check if peptide is above percent cutoff for delta between predicted and observed RT
        my $pepLine = '';  # line to print if ok
        my $anaKey = (keys %{$hydroResults{$seqZ}})[0];
        my $HI = $hydroResults{$seqZ}{$anaKey}->[2];
        my ($seq, $charge) = split('_', $seqZ);

        $pepLine = qq
|<TR bgcolor="$lineColor">
<TD class="rBorder" align="left" nowrap>&nbsp;&nbsp;$seq&nbsp;&nbsp;</TD>
<TD class="rBorder" align="center" nowrap>&nbsp;&nbsp;$charge+&nbsp;&nbsp;</TD>
<TD class="rBorder" align="center" nowrap>&nbsp;&nbsp;${\sprintf("%.2f", $HI)}&nbsp;&nbsp;</TD>
|;

        foreach my $anaID (sort { $anaOrder{$a} <=> $anaOrder{$b} } keys %anaOrder) {
            my ($pepID, $RT, $hi, $validStatus, $predRT, $timeShift, $percentile, $rtClass);
            my $displayValue;
            my $popupStrg;
            if ($hydroResults{$seqZ}{$anaID}) {
                ($pepID, $RT, $hi, $validStatus, $predRT, $timeShift, $percentile) = @{$hydroResults{$seqZ}{$anaID}};
                $displayValue = ($measure eq 'shift')? $timeShift : ($measure eq 'obsRT')? $RT : $predRT;
                if (100 * $percentile >= $pct_cutoff) {
                    $displayPeptide = 1;
                    $displayValue = ${\sprintf("%.2f", $displayValue)};
                } else {
                    $displayValue = "(" . ${\sprintf("%.2f", $displayValue)} . ")";
                }
                $popupStrg = "<B>$seq $charge+</B>&nbsp;&nbsp;&nbsp;&nbsp;";
                $popupStrg .= "<br>Analysis : $hydroInfo{'ana'}{$anaID}->[5]";
                $popupStrg .= "<br><U>Retention Time (min)</U>";
                $popupStrg .= "<br>Obs. : " . ${\sprintf("%.2f", $RT)};
                $popupStrg .= "<br>Pred. : " . ${\sprintf("%.2f", $predRT)};
                $popupStrg .= "<br>Shift : " . ${\sprintf("%.2f", $timeShift)};
                $popupStrg .= "<br>Percentile : " . ${\sprintf("%.2f", 100 * $percentile)};
                $popupStrg = "onmouseover=\"popup('$popupStrg')\" onmouseout=\"popout()\"";
            } else {
                $displayValue = "-";
                $popupStrg = "";
                $validStatus = 1;
            }
            $rtClass = ($validStatus == 1)? "" : "class=\"virtualProt\"";
            $pepLine .= "<TD class=\"rBorder\" align=\"center\" nowrap $popupStrg><FONT $rtClass>&nbsp;&nbsp;$displayValue&nbsp;&nbsp;</FONT></TD>\n";
        }
        $pepLine .= "</TR>\n";
        if ($displayPeptide) {
            print "$pepLine";
            $dispIonCount++;
            $lineColor = ($lineColor eq $lightColor)? $darkColor : $lightColor;
        }
    }
    my $pepNbStrg = "<A href=\"javascript:sortTable(\\'seq\\')\"><FONT color=\"$seqColor\">&nbsp;&nbsp;Peptide ions ($dispIonCount / $totalIonCount)&nbsp;&nbsp;</FONT></A>";
    print qq
|   </TABLE>
    End of peptides list.
<SCRIPT type="text/javascript">
    document.getElementById('pepSeqTH').innerHTML='$pepNbStrg';
</SCRIPT>
|;
}


sub displayLinearModels {  # Globals: %hydroInfo, %linMods
    my $numAna = scalar keys %{$hydroInfo{'ana'}};
    my $numSamp = scalar keys %{$hydroInfo{'samp'}};

    print qq
|<BR><BR><BR>
Linear models corresponding to the equation : <B>HI = slope x RT + intercept</B>
<BR>
The slope should correspond to the Acetonitrile gradient of the HPLC.
<BR><BR>
<TABLE id="linModTable" border=0 cellspacing=0 align=center>
<TR bgcolor=$darkColor>
    <TH class="rbBorder" valign="center" nowrap>
        &nbsp;&nbsp;Sample&nbsp;&nbsp;
    </TH>
    <TH class="rbBorder" valign="center" nowrap>
        &nbsp;&nbsp;Analysis&nbsp;&nbsp;
    </TH>
    <TH class="rbBorder" valign="center" nowrap>
        &nbsp;&nbsp;Slope&nbsp;&nbsp;
    </TH>
    <TH class="rbBorder" valign="center" nowrap>
        &nbsp;&nbsp;Intercept&nbsp;&nbsp;
    </TH>
    <TH class="rbBorder" valign="top" nowrap>
        &nbsp;&nbsp;R<SUP>2</SUP>&nbsp;&nbsp;
    </TH>
</TR>
|;

    my $lineColor = $lightColor;
    foreach my $sampID (sort {$hydroInfo{'samp'}{$a}->[1] <=> $hydroInfo{'samp'}{$b}->[1]} keys %{$hydroInfo{'samp'}}) {
        print "<TR bgcolor=$lineColor>\n";
        my $sampAnaNb = scalar @{$hydroInfo{'samp'}{$sampID}->[2]};
        print qq
|   <TH class="rbBorder" bgColor=$darkColor valign="center" rowspan=$sampAnaNb nowrap>
        &nbsp;$hydroInfo{'samp'}{$sampID}->[0]&nbsp;
    </TH>
|;

        my $anaCount = 1;
        foreach my $anaID (sort {$hydroInfo{'ana'}{$a}->[6] <=> $hydroInfo{'ana'}{$b}->[6]} @{$hydroInfo{'samp'}{$sampID}->[2]}) {
            print "<TR bgcolor=$lineColor>\n" unless ($anaCount == 1);
            print qq
|   <TH class="rbBorder" valign="center" nowrap>
        &nbsp;$hydroInfo{'ana'}{$anaID}->[5]&nbsp;
    </TH>
    <TD class="rbBorder" valign="center" nowrap>
        &nbsp;&nbsp;$linMods{$anaID}->[0]&nbsp;&nbsp;
    </TD>
    <TD class="rbBorder" valign="center" nowrap>
        &nbsp;&nbsp;$linMods{$anaID}->[1]&nbsp;&nbsp;
    </TD>
    <TD class="rbBorder" valign="center" nowrap>
        &nbsp;&nbsp;$linMods{$anaID}->[3]&nbsp;&nbsp;
    </TD>
</TR>
|;
            $anaCount++;
            $lineColor = ($lineColor eq $lightColor)? $darkColor : $lightColor;
        }
    }

    print "</TABLE>\n";
}


sub displayGraphHydroRT {  # Globals: $fromItem, $fromID, $hydroGraph, $hydroGraphDisplay
    
    my @scopeItems = ($fromItem eq "analysis")? ("experiment", "sample", "analysis") : 
                    ($fromItem eq "sample")? ("experiment", "sample") : ("experiment");
    my $showItem = param('SCOPE') || $fromItem;
    my $size = 1000;

    print qq
|   &nbsp;&nbsp;&nbsp;&nbsp;
    <LABEL for='scopeType' style='font-weight:bold'>Visualization Scope: </LABEL>
    <SELECT id='scopeType' onchange="window.location='$promsPath{cgi}/displayHydro.cgi?ACT=display&VIEW=$view&ID=$quantifID&FROM=$fromItem:$fromID&SCOPE=' + this.value;">
|;
    foreach my $scope (@scopeItems) {
        my $disabStrg = (lc($scope) eq $showItem)? " selected disabled" : " ";
        my $scopeStrg = ucfirst($scope);
        print "<OPTION value=\"$scope\"$disabStrg>$scopeStrg</OPTION>\n";
    }
    print "</SELECT>\n<BR><BR>\n";

    my @hydroGraphs;
    my @hydroGraphsDisplay;
    if ($showItem eq "experiment") {
        push @hydroGraphs, "$hydroGraph\_exp.png";
        push @hydroGraphs, "$hydroGraph\_exp_shift.png";
        push @hydroGraphs, "$hydroGraph\_exp_shift_RT.png";
        push @hydroGraphsDisplay, "$hydroGraphDisplay\_exp.png";
        push @hydroGraphsDisplay, "$hydroGraphDisplay\_exp_shift.png";
        push @hydroGraphsDisplay, "$hydroGraphDisplay\_exp_shift_RT.png";
        $size = 1200;
    } elsif ($fromItem eq "analysis" && $showItem eq "analysis") {
        my $showID = $fromID;
        push @hydroGraphs, "$hydroGraph\_ana_$showID.png";
        push @hydroGraphs, "$hydroGraph\_ana_shift_$showID.png";
        push @hydroGraphs, "$hydroGraph\_ana_shift_RT_$showID.png";
        push @hydroGraphsDisplay, "$hydroGraphDisplay\_ana_$showID.png";
        push @hydroGraphsDisplay, "$hydroGraphDisplay\_ana_shift_$showID.png";
        push @hydroGraphsDisplay, "$hydroGraphDisplay\_ana_shift_RT_$showID.png";
        $size = 600;
    } else {  # showItem must be "sample"
        my $showID;
        if ($fromItem eq "sample") {  # $showItem must be "sample" too then
            $showID = $fromID;
        } elsif ($fromItem eq "analysis") {
            my $dbh = &promsConfig::dbConnect;
            ($showID) = $dbh->selectrow_array("SELECT ID_SAMPLE FROM ANALYSIS WHERE ID_ANALYSIS = $fromID");
            $dbh->disconnect;
        }
        push @hydroGraphs, "$hydroGraph\_samp_$showID.png";
        push @hydroGraphs, "$hydroGraph\_samp_shift_$showID.png";
        push @hydroGraphs, "$hydroGraph\_samp_shift_RT_$showID.png";
        push @hydroGraphsDisplay, "$hydroGraphDisplay\_samp_$showID.png";
        push @hydroGraphsDisplay, "$hydroGraphDisplay\_samp_shift_$showID.png";
        push @hydroGraphsDisplay, "$hydroGraphDisplay\_samp_shift_RT_$showID.png";
        $size = 900;
    }

    my $graphsStrg = '';
    print "Export as PNG : ";

    for my $idx (0..$#hydroGraphs) {
        my $graphType = ($idx == 0)? "HI vs observed RT" : ($idx == 1)? "Shift distribution" : "Shift vs observed RT";
        $graphsStrg .= "<BR><BR><FONT class=\"title3\">$graphType</FONT><BR>";
        if (-e $hydroGraphs[$idx]) {
            my $graph = $hydroGraphsDisplay[$idx];            
            print qq
|   &nbsp;&nbsp;
    <A href="$graph" download="hydro_graph.png">
        <INPUT type="button" class="font11" value=" $graphType ">
    </A>
|;
            $graphsStrg .= "<IMG src=\"$graph\" width=$size heigth=$size><BR><BR>";
        } else {
            $graphsStrg .= "<FONT class=\"title2\" color=\"red\">Graph not available</FONT><BR><BR>\n";
        }
    }
    print $graphsStrg;
}


sub getHydroInfo {  # Globals: $dbh, $quantifID, $quantifAnnot, $upDate, $user
    my ($refHydroInfo) = @_;

    # Experiment, Sample & Analysis Info
    my $sthGetInfo = $dbh->prepare("SELECT A.ID_ANALYSIS, A.NAME, A.DISPLAY_POS, 
        S.ID_SAMPLE, S.NAME, S.DISPLAY_POS, E.ID_EXPERIMENT, E.NAME
        FROM QUANTIFICATION Q
        INNER JOIN ANA_QUANTIFICATION AQ ON Q.ID_QUANTIFICATION = AQ.ID_QUANTIFICATION
        INNER JOIN ANALYSIS A ON AQ.ID_ANALYSIS = A.ID_ANALYSIS
        INNER JOIN SAMPLE S ON A.ID_SAMPLE = S.ID_SAMPLE
        INNER JOIN EXPERIMENT E ON S.ID_EXPERIMENT = E.ID_EXPERIMENT
        WHERE Q.ID_QUANTIFICATION = $quantifID"
    );
    $sthGetInfo->execute;
    while (my ($anaID, $anaName, $anaPos, $sampID, $sampName, $sampPos, $expID, $expName) = $sthGetInfo->fetchrow_array) {
        $refHydroInfo->{'ana'}{$anaID} = [$expID, $expName, $sampID, $sampName, $sampPos, $anaName, $anaPos];
        $refHydroInfo->{'samp'}{$sampID} = [$sampName, $sampPos, []] unless $refHydroInfo->{'samp'}{$sampID};
        push @{$refHydroInfo->{'samp'}{$sampID}[2]}, $anaID;
    }
    $sthGetInfo->finish;
    
    # Quantif Annot Info
    ($refHydroInfo->{'rtOffset'}) = ($quantifAnnot =~ /RT_OFFSET=([^:]+)/);
    my ($anaPepQuantifStrg) = ($quantifAnnot =~ /ANA_PEPQUANTIF=([^:]+)/);
    $anaPepQuantifStrg =~ s/#//g;  # Remove ID tags
    my @anaPepQuantifs = split(';', $anaPepQuantifStrg);
    my $sthPepQuantif = $dbh->prepare("SELECT NAME FROM QUANTIFICATION WHERE ID_QUANTIFICATION = ?");
    foreach my $anaPepQuantif (@anaPepQuantifs) {
        my ($anaID, $pepQuantifID) = split(',', $anaPepQuantif);
        if ($pepQuantifID eq "no_alignment") {
            push @{$refHydroInfo->{'ana'}{$anaID}}, ($pepQuantifID, "No peptide quantification");
        } else {
            $sthPepQuantif->execute($pepQuantifID);
            my ($pepQuantifName) = $sthPepQuantif->fetchrow_array;
            push @{$refHydroInfo->{'ana'}{$anaID}}, ($pepQuantifID, $pepQuantifName);
        }
    }
    $sthPepQuantif->finish;

    # Creation date and user info
    $refHydroInfo->{'updateDate'} = &promsMod::formatDate($upDate);
    if ($user) {
        my ($userName) = $dbh->selectrow_array("SELECT USER_NAME FROM USER_LIST WHERE ID_USER = '$user'");
        $refHydroInfo->{'user'} = $userName || $user;
    } else {
        $refHydroInfo->{'user'} = '?';
    }
}


sub parseHydroResults {
    my ($refHydroResults, $results) = @_;

    open(RESULTS, $results);
    # Analysis_ID Sample_ID Peptide_ID Sequence RT Valid_status HI time_shift percentile
    while (my $line = <RESULTS>) {
        next if ($. == 1);  # Skip headers
        chomp $line;
        my @values = split(/\t/, $line);
        my $anaID       = $values[0];
        my $pepID       = $values[2];
        my $seq         = $values[3];
        my $RT          = $values[4];
        my $charge      = $values[5];
        my $validStatus = $values[6];
        my $HI          = $values[7];
        my $predRT      = $values[8];
        my $timeShift   = $values[9];
        my $percentile  = $values[10];
        my $seqZ        = "$seq\_$charge";

        $refHydroResults->{$seqZ}{$anaID} = [$pepID, $RT, $HI, $validStatus, $predRT, $timeShift, $percentile] unless (exists($refHydroResults->{$seqZ}) && exists($refHydroResults->{$seqZ}{$anaID}));
    }
    close RESULTS;
}


sub parseLinMods {
    my ($refLinMods, $lmResults) = @_;
    open(RESULTS, $lmResults);
    # Analysis_ID Sample_ID slope intercept R R2
    while (my $line = <RESULTS>) {
        next if ($. == 1);  # Skip headers
        chomp $line;
        my @values = split(/\t/, $line);
        my $anaID     = $values[0];
        my $slope     = $values[2];
        my $intercept = $values[3];
        my $Rcorr     = $values[4];
        my $R2        = $values[5];
        
        $refLinMods->{$anaID} = [$slope, $intercept, $Rcorr, $R2];
    }
    close RESULTS;
}


sub parseUnmodifProps {
    my ($refUnmodifProps, $unmodifPropsFile) = @_;
    open(RESULTS, $unmodifPropsFile);
    # Analysis_ID Analysis Sample_ID Sample Unmodified_pep_nb Total_pep_nb Unmodified_proportion
    while (my $line = <RESULTS>) {
        next if ($. == 1);  # Skip headers
        chomp $line;
        my @values = split(/\t/, $line);
        my $anaID          = $values[0];
        my $unmodifiedNb   = $values[2];
        my $totalPepNb     = $values[3];
        my $unmodifiedProp = $values[4];

        $refUnmodifProps->{$anaID} = [$unmodifiedNb, $totalPepNb, $unmodifiedProp];
    }
    close RESULTS;
}


sub displaySummary {  # Globals: $quantifID, %hydroInfo, $name, $status, $lightColor, $darkColor
    print header(-charset=>'utf-8');
    warningsToBrowser(1);
    my $title = "Hydrophobicity analysis : <FONT color=\"red\">$name</FONT>";

    print qq
|<HTML>
<HEAD>
<TITLE>Hydrophobicity Analysis</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT langage="Javascript">
|;

    if ($status == 0) {
        print qq
|setTimeout(function(){parent.optionFrame.location.reload(1);},5000);  // Reload page every 5 seconds
|;
    }

    print qq
|</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
|;
    print &displayBackButton;
    print qq
|<CENTER>
<FONT class="title">$title</FONT>
<BR><BR>
<TABLE border="0" cellspacing="2" cellpadding="2" bgcolor="$darkColor">
|;

    # Format parameters and analyses
    my $numAna = scalar keys %{$hydroInfo{'ana'}};
    my $anaStrg = ($numAna > 1)? "Samples, Analyses \&<BR>Peptide Quantifications" : "Sample, Analysis \&<BR>Peptide Quantification";

    my $statusStrg;
    if ($status == -2) {
        $statusStrg = "<IMG src=$promsPath{images}/bad.gif>";
        $statusStrg .= "&nbsp;&nbsp;<FONT color=\"red\">***ERROR***</FONT>";
        $statusStrg .= "&nbsp;&nbsp;(Check what happened with the \"Monitor Jobs\" button at the top)"
    } elsif ($status == -1) {
        $statusStrg = "<B>Not launched yet...</B>";
    } elsif ($status == 0) {
        $statusStrg = "<FONT color=\"orange\"><B>0n-going...</B></FONT>";
    } else {
        $statusStrg = "<IMG src=$promsPath{images}/good.gif>";
        $statusStrg .= "&nbsp;<FONT color=\"green\">Finished</FONT>";
    }

    print qq
|<TR>
    <TH align="right" width="150" nowrap>Name :</TH>
    <TD bgcolor="$lightColor">&nbsp;&nbsp;$name</TD>
</TR>
<TR>
    <TH align="right" valign="top" nowrap>&nbsp;$anaStrg</TH>
    <TD nowrap bgcolor=$lightColor style="min-width:400px">
        <DIV style="max-height:200px;overflow:auto">
            <TABLE width=100% style="white-space:nowrap;">
|;
    foreach my $sampID (sort {$hydroInfo{'samp'}{$a}->[1] <=> $hydroInfo{'samp'}{$b}->[1]} keys %{$hydroInfo{'samp'}}) {  # Sort by sample display pos
        my $sampAnaNb = scalar @{$hydroInfo{'samp'}{$sampID}->[2]};
        print "<TR>\n";
        print "<TH valign=\"top\" align=\"right\" rowspan=\"$sampAnaNb\">&nbsp;$hydroInfo{'samp'}{$sampID}->[0] > </TH>\n";

        my $anaCount = 0;
        foreach my $anaID (sort {$hydroInfo{'ana'}{$a}->[6] <=> $hydroInfo{'ana'}{$b}->[6]} @{$hydroInfo{'samp'}{$sampID}->[2]}) {  # sort by analysis display pos
            $anaCount++;
            my $anaPepQuantifStrg = "&nbsp;&nbsp;&bull;&nbsp;$hydroInfo{'ana'}{$anaID}->[5] : $hydroInfo{'ana'}{$anaID}->[8]";
            if ($sampAnaNb == 1) {
                print "<TD valign=\"top\">$anaPepQuantifStrg</TD>\n";
            } else {
                print "<TR>\n" unless ($anaCount == 1);
                print "<TD valign=\"top\">$anaPepQuantifStrg</TD>\n";
                print "</TR>\n" unless ($anaCount == $sampAnaNb);
            }
        }
        print "</TR>\n";
    }
    print qq
|           </TABLE>
        </DIV>
    </TD>
</TR>
<TR>
    <TH align="right" valign="top" nowrap>&nbsp;Parameters :</TH>
    <TD nowrap bgcolor=$lightColor>    
        &nbsp;&nbsp;&bull;&nbsp;RT offset : $hydroInfo{'rtOffset'} min
    </TD>
</TR>
<TR>
    <TH align="right" bgcolor="$darkColor" nowrap>Status :</TH>
    <TD nowrap bgcolor="$lightColor">&nbsp;&nbsp;$statusStrg</TD>
</TR>
<TR>
    <TH align="right" bgcolor="$darkColor" nowrap>Creation date :</TH>
    <TD nowrap bgcolor="$lightColor">&nbsp;&nbsp;$hydroInfo{'updateDate'} by <B>$hydroInfo{'user'}</B></TD>
</TR>
</TABLE>
<BR>
<FONT class="font11" style="font-weight:bold">Hydrophobicity computation based on <A href="http://hs2.proteome.ca/SSRCalc/SSRCalcQ.html" target="_blank">SSRCalc</A> algorithm (<A href="https://doi.org/10.1074/mcp.M400031-MCP200" target="_blank">Krokhin O.V. et al. Mol. Cell. Proteomics, 2004</A>).
</CENTER>
</BODY>
</HTML>
|;

}


sub displayFormEdit {  # Globals: $quantifID, $name, $lightColor, $darkColor
    print header(-charset=>'utf-8');
    warningsToBrowser(1);
    my $title = "Editing Hydrophobicity analysis : <FONT color=\"red\">$name</FONT>";
    
    print qq
|<HTML>
<HEAD>
<TITLE>Hydrophobicity analysis</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT langage="Javascript">
function checkForm(myForm) {
    if (!myForm.quantifName.value) {
        alert('A name is expected for the analysis');
        return false;
    }
    return true;
}
function cancelAction() {
    top.promsFrame.selectedAction = 'summary';
    top.promsFrame.optionFrame.selectOption();
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">$title</FONT>
<BR><BR>
<FORM NAME="displayMenu" method="POST" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="ACT" value="$action">
<INPUT type="hidden" name="ID" value="$quantifID">
<INPUT type="hidden" name="FROM" value="$fromItem:$fromID">
<TABLE border="0" cellspacing="2" cellpadding="2" bgcolor="$darkColor">
<TR>
    <TH align="right" nowrap>Name :</TH>
    <TD bgcolor="$lightColor">
        <INPUT type="text" name="quantifName" size="50" maxlength="100" value="$name"/>
    </TD>
</TR>
<TR>
    <TH colspan="2">
        <INPUT type="submit" name="submit" value="Update">
        &nbsp;&nbsp;
        <INPUT type="button" value="Cancel" onclick="cancelAction()">
    </TH>
</TR>
</TABLE><BR>
</FORM>
</CENTER>
</BODY>
</HTML>
|;
}


sub editHydro {  # Globals: $quantifID, $fromItem, $fromID
    my $newQuantifName = param('quantifName');
    
    my $dbh = &promsConfig::dbConnect;
    $dbh->do("UPDATE QUANTIFICATION SET NAME = '$newQuantifName' WHERE ID_QUANTIFICATION = $quantifID") or die "Cannot prepare: " . $dbh->errstr();
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
    top.promsFrame.resultFrame.location = "$promsPath{cgi}/displayHydro.cgi?ACT=summary&ID=$quantifID&FROM=$fromItem:$fromID";
</SCRIPT>
</HEAD>
</HTML>
|;
}


sub deleteHydro {  # Globals: $quantifID, %promsPath, $projectID, $fromItem, $fromID
    print header(-'content-encoding'=>'no',-charset=>'utf-8');
    warningsToBrowser(1);

    my $dbh = &promsConfig::dbConnect;
    $dbh->do("DELETE from PARENT_QUANTIFICATION where ID_QUANTIFICATION = $quantifID");
    $dbh->do("DELETE from ANA_QUANTIFICATION where ID_QUANTIFICATION = $quantifID");
    $dbh->do("DELETE from QUANTIFICATION where ID_QUANTIFICATION = $quantifID");
    $dbh->commit;
    $dbh->disconnect;

    my $pathToFile = "$promsPath{'quantification'}/project_$projectID/quanti_$quantifID";
    rmtree($pathToFile);

    print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
    parent.navFrame.location="$promsPath{cgi}/openProject.cgi?ACT=nav&ID=$projectID&branchID=$fromItem:$fromID";
</SCRIPT>
</HEAD>
</HTML>
|;
}


sub displayBackButton {  # Globals: $projectID, $fromItem, $fromID
    my $item = ucfirst(lc($fromItem));
    return "<INPUT type=\"button\" id=\"backButton\" value=\"<< Back to $item\" onclick=\"parent.navFrame.location='$promsPath{cgi}/openProject.cgi?ACT=nav&ID=$projectID&branchID=$fromItem:$fromID'\">";
}


sub exportHydroResults {  # Globals: %hydroResults, %hydroInfo, %linMods, %unmodifProps $name, $desc
    my $workbook = Spreadsheet::WriteExcel->new("-");
    my $worksheet1 = $workbook->add_worksheet('Parameters');
    my $worksheet2 = $workbook->add_worksheet('Results');
    my $worksheet3 = $workbook->add_worksheet('Linear_models');
    my $worksheet4 = $workbook->add_worksheet('Pep_proportions');
    my $worksheet5 = $workbook->add_worksheet('Info');

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
    my $mergeFormat = $workbook->add_format(align => 'left', valign => 'vcenter');
    $mergeFormat->set_text_wrap();
    my $multilineFormat = $workbook->add_format(align => 'left', valign => 'vcenter');
    $multilineFormat->set_text_wrap();
    my $boldFormat = $workbook->add_format(align => 'left', valign => 'vcenter');
    $boldFormat->set_bold();
    
    # General variables
    my $numAna = scalar keys %{$hydroInfo{'ana'}};
    my $numSamp = scalar keys %{$hydroInfo{'samp'}};
    my $totalIonCount = scalar keys %hydroResults;
    my $pepNbStrg = "Peptide Ions ($totalIonCount)";
    my (%sampCols, %anaCols);
    my @measures = ('Shift', 'Obs. RT', 'Pred. RT');
    
    my $row = 0;
    my $colNb = 3;
    foreach my $sampID (sort {$hydroInfo{'samp'}{$a}->[1] <=> $hydroInfo{'samp'}{$b}->[1]} keys %{$hydroInfo{'samp'}}) {
        # my $sampAnaNb = scalar @{$hydroInfo{'samp'}{$sampID}->[2]};
        $sampCols{$sampID}{'start'} = $colNb;

        foreach my $anaID (sort {$hydroInfo{'ana'}{$a}->[6] <=> $hydroInfo{'ana'}{$b}->[6]} @{$hydroInfo{'samp'}{$sampID}->[2]}) {
            $anaCols{$anaID}{'start'} = $colNb;
            foreach my $meas (@measures) {
                $anaCols{$anaID}{$meas} = $colNb;
                $colNb++;
            }
            $anaCols{$anaID}{'end'} = $colNb - 1;
        }
        $sampCols{$sampID}{'end'} = $colNb - 1;
    }


    ### Writing Parameters sheet
    # Format columns
    $worksheet1->set_row(0, 30);
    $worksheet1->set_column("A:A", 30);
    $worksheet1->set_column("B:B", 20);
    $worksheet1->set_column("C:C", 10);
    $worksheet1->set_column("D:D", 40);
    
    # Write column headers
    $worksheet1->write_string("A1", "Hydrophobicity analysis name", $headerFormat);
    $worksheet1->merge_range("B1:D1", $name, $headerMergeFormat);
    $worksheet1->merge_range(2, 0, 2 + $numAna, 0, "Data used", $headerMergeFormat);
    $worksheet1->write_string("B3", "Sample", $headerFormat);
    $worksheet1->write_string("C3", "Analysis", $headerFormat);
    $worksheet1->write_string("D3", "Peptide Quantification", $headerFormat);

    $row = 3;
    foreach my $sampID (sort {$hydroInfo{'samp'}{$a}->[1] <=> $hydroInfo{'samp'}{$b}->[1]} keys %{$hydroInfo{'samp'}}) {  # Sort by sample display pos
        my $sampAnaNb = scalar @{$hydroInfo{'samp'}{$sampID}->[2]};
        if ($sampAnaNb > 1) {
            $worksheet1->merge_range($row, 1, $row + $sampAnaNb - 1, 1, $hydroInfo{'samp'}{$sampID}->[0], $headerMergeFormat);
        } else {
            $worksheet1->write_string($row, 1, $hydroInfo{'samp'}{$sampID}->[0], $headerFormat);
        }

        foreach my $anaID (sort {$hydroInfo{'ana'}{$a}->[6] <=> $hydroInfo{'ana'}{$b}->[6]} @{$hydroInfo{'samp'}{$sampID}->[2]}) {  # sort by analysis display pos
            $worksheet1->write_string($row, 2, $hydroInfo{'ana'}{$anaID}->[5]);
            $worksheet1->write_string($row, 3, $hydroInfo{'ana'}{$anaID}->[8]);
            $row++;
        }
    }
    
    $row++;  # Keep a blank line
    $worksheet1->write_string($row, 0, "RT offset", $headerFormat);
    $worksheet1->write_string($row++, 1, "$hydroInfo{'rtOffset'} min");
    $worksheet1->write_string($row, 0, "Creation date", $headerFormat);
    $worksheet1->write_string($row++, 1, "$hydroInfo{'updateDate'} by $hydroInfo{'user'}");
    $row++;  # Another blank line
    $worksheet1->set_row($row, 20);
    $worksheet1->merge_range($row, 0, $row, 3, "Hydrophobicity computation based on SSRCalc algorithm (http://hs2.proteome.ca/SSRCalc/SSRCalcQ.html)", $mergeFormat);
    $row++;
    $worksheet1->set_row($row, 20);
    $worksheet1->merge_range($row, 0, $row, 3, "Krokhin O.V. et al. Mol. Cell. Proteomics, 2004 (https://doi.org/10.1074/mcp.M400031-MCP200)", $mergeFormat);
    ### End of Parameters sheet

    ### Writing Results sheet
    $row = 0;
    # Format columns
    $worksheet2->set_column("A:A", 40);
    $worksheet2->set_column("B:B", 10);
    $worksheet2->set_column("C:C", 20);
    
    # Write column headers
    $worksheet2->merge_range("A1:B2", $pepNbStrg, $headerMergeFormat);
    $worksheet2->merge_range("A3:A4", "Peptide Sequence", $headerMergeFormat);
    $worksheet2->merge_range("B3:B4", "Charge", $headerMergeFormat);
    $worksheet2->merge_range("C1:C4", "Hydrophobicity\nIndex\n(ACN%)", $headerMergeFormat);
    $worksheet2->merge_range(0, 3, 0, 2 + $numAna * 3, "Retention Time (min)", $headerMergeFormat);
    foreach my $sampID (keys %sampCols) {
        $worksheet2->merge_range(1, $sampCols{$sampID}{'start'}, 1, $sampCols{$sampID}{'end'}, $hydroInfo{'samp'}{$sampID}->[0], $headerMergeFormat);
    }
    foreach my $anaID (keys %anaCols) {
        $worksheet2->merge_range(2, $anaCols{$anaID}{'start'}, 2, $anaCols{$anaID}{'end'}, $hydroInfo{'ana'}{$anaID}->[5], $headerMergeFormat);
        foreach my $meas (@measures) {
            $worksheet2->write_string(3, $anaCols{$anaID}{$meas}, $meas, $headerFormat);
        }
    }

    $row = 4;
    foreach my $seqZ (sort {$a cmp $b} keys %hydroResults) {
        my $anaKey = (keys %{$hydroResults{$seqZ}})[0];
        my $HI = ${\sprintf("%.2f", $hydroResults{$seqZ}{$anaKey}->[2])};
        my ($seq, $charge) = split('_', $seqZ);
        $charge = "$charge+";

        $worksheet2->write_string($row, 0, $seq);
        $worksheet2->write_string($row, 1, $charge);
        $worksheet2->write_number($row, 2, $HI);

        foreach my $anaID (keys %anaCols) {
            my ($pepID, $RT, $hi, $validStatus, $predRT, $timeShift, $percentile, $rtClass);
            
            if ($hydroResults{$seqZ}{$anaID}) {
                ($pepID, $RT, $hi, $validStatus, $predRT, $timeShift, $percentile) = @{$hydroResults{$seqZ}{$anaID}};
                $RT = ${\sprintf("%.2f", $RT)};
                $predRT = ${\sprintf("%.2f", $predRT)};
                $timeShift = ${\sprintf("%.2f", $timeShift)};
            } else {
                $RT = "NA";
                $predRT = "NA";
                $timeShift = "NA";
                $validStatus = 1;
            }

            ($timeShift eq "NA")? $worksheet2->write_string($row, $anaCols{$anaID}{'Shift'}, $timeShift) :
                                  $worksheet2->write_number($row, $anaCols{$anaID}{'Shift'}, $timeShift);
            ($RT eq "NA")? $worksheet2->write_string($row, $anaCols{$anaID}{'Obs. RT'}, $RT) :
                           $worksheet2->write_number($row, $anaCols{$anaID}{'Obs. RT'}, $RT);
            ($predRT eq "NA")? $worksheet2->write_string($row, $anaCols{$anaID}{'Pred. RT'}, $predRT) :
                               $worksheet2->write_number($row, $anaCols{$anaID}{'Pred. RT'}, $predRT);
        }
        $row++;
    }
    ### End of results sheet
    
    ### Writing Linear models sheet
    $row = 0;
    # Format columns
    $worksheet3->set_column("A:A", 30);
    $worksheet3->set_column("B:B", 10);
    
    $worksheet3->write_string($row, 0, "Sample", $headerFormat);
    $worksheet3->write_string($row, 1, "Analysis", $headerFormat);
    $worksheet3->write_string($row, 2, "Slope", $headerFormat);
    $worksheet3->write_string($row, 3, "Intercept", $headerFormat);
    $worksheet3->write_string($row, 4, "R2", $headerFormat);

    $row++;
    foreach my $sampID (sort {$hydroInfo{'samp'}{$a}->[1] <=> $hydroInfo{'samp'}{$b}->[1]} keys %{$hydroInfo{'samp'}}) {
        my $sampAnaNb = scalar @{$hydroInfo{'samp'}{$sampID}->[2]};
        if ($sampAnaNb > 1) {
            $worksheet3->merge_range($row, 0, $row + $sampAnaNb - 1, 0, $hydroInfo{'samp'}{$sampID}->[0], $headerMergeFormat);
        } else {
            $worksheet3->write_string($row, 0, $hydroInfo{'samp'}{$sampID}->[0], $headerFormat);
        }
        
        foreach my $anaID (sort {$hydroInfo{'ana'}{$a}->[6] <=> $hydroInfo{'ana'}{$b}->[6]} @{$hydroInfo{'samp'}{$sampID}->[2]}) {
            $worksheet3->write_string($row, 1, $hydroInfo{'ana'}{$anaID}->[5]);
            $worksheet3->write_number($row, 2, $linMods{$anaID}->[0]);
            $worksheet3->write_number($row, 3, $linMods{$anaID}->[1]);
            $worksheet3->write_number($row, 4, $linMods{$anaID}->[3]);
            $row++;
        }
    }
    $row++;  # Blank line
    $worksheet3->set_row($row, 20);
    $worksheet3->merge_range($row, 0, $row, 4, "Equation to fit linear models : HI = slope x RT + intercept", $mergeFormat);
    $row++;
    $worksheet3->set_row($row, 20);
    $worksheet3->merge_range($row, 0, $row, 4, "The slope should correspond to the Acetonitrile gradient of the HPLC.", $mergeFormat);
    ### End of Linear models sheet


    ### Writing Peptide Proportions sheet
    $row = 0;
    # Format columns
    $worksheet4->set_row($row, 30);
    $worksheet4->set_column("A:A", 30);
    $worksheet4->set_column("B:B", 10);
    $worksheet4->set_column("C:C", 25);
    $worksheet4->set_column("D:D", 25);
    $worksheet4->set_column("E:E", 30);

    $worksheet4->write_string($row, 0, "Sample", $headerFormat);
    $worksheet4->write_string($row, 1, "Analysis", $headerFormat);
    $worksheet4->write_string($row, 2, "Nb of unmodified peptides (or carbamidomethyl)", $headerFormat);
    $worksheet4->write_string($row, 3, "Total number of peptides", $headerFormat);
    $worksheet4->write_string($row, 4, "Proportion of unmodified peptides (or carbamidomethyl)", $headerFormat);

    $row++;
    foreach my $sampID (sort {$hydroInfo{'samp'}{$a}->[1] <=> $hydroInfo{'samp'}{$b}->[1]} keys %{$hydroInfo{'samp'}}) {
        my $sampAnaNb = scalar @{$hydroInfo{'samp'}{$sampID}->[2]};
        if ($sampAnaNb > 1) {
            $worksheet4->merge_range($row, 0, $row + $sampAnaNb - 1, 0, $hydroInfo{'samp'}{$sampID}->[0], $headerMergeFormat);
        } else {
            $worksheet4->write_string($row, 0, $hydroInfo{'samp'}{$sampID}->[0], $headerFormat);
        }
        
        foreach my $anaID (sort {$hydroInfo{'ana'}{$a}->[6] <=> $hydroInfo{'ana'}{$b}->[6]} @{$hydroInfo{'samp'}{$sampID}->[2]}) {
            $worksheet4->write_string($row, 1, $hydroInfo{'ana'}{$anaID}->[5]);
            $worksheet4->write_number($row, 2, $unmodifProps{$anaID}->[0]);
            $worksheet4->write_number($row, 3, $unmodifProps{$anaID}->[1]);
            $worksheet4->write_number($row, 4, $unmodifProps{$anaID}->[2]);
            $row++;
        }
    }
    ### End of Peptide Proportions sheet

    ### Writing Info sheet
    $worksheet5->set_column("A:A", 40);
    $worksheet5->set_column("B:B", 100);
    foreach my $i (1..7) {
        $worksheet5->set_row($i, 35);
    }

    $worksheet5->write_string("A1", "Feature", $headerFormat);
    $worksheet5->write_string("B1", "Description", $headerFormat);
    $worksheet5->write_string("A2", "RT offset", $boldFormat);
    $worksheet5->write_string("B2", "The time (in minutes) at the beginning of the chromatographic run when the peptides' elution is not linear with the acetonitrile gradient. It corresponds roughly to the dead volume + the injection peak.", $multilineFormat);
    $worksheet5->write_string("A3", "Peptide Ions", $boldFormat);
    $worksheet5->write_string("B3", "Unique peptide ions found in the different analyses used in the computation", $multilineFormat);
    $worksheet5->write_string("A4", "Peptide Sequence", $boldFormat);
    $worksheet5->write_string("B4", "The Amino Acids sequence. Cysteines are in fact Carbamidomethyl Cysteines. Other modifications are not considered for the computation.", $multilineFormat);
    $worksheet5->write_string("A5", "Hydrophobicity Index (ACN%)", $boldFormat);
    $worksheet5->write_string("B5", "The peptide Hydrophobicity measured as the percentage of Acetonitrile at which the peptide is predicted to come out from the LC column. It is used to predict the retention time, assuming a linear ACN gradient.", $multilineFormat);
    $worksheet5->write_string("A6", "Shift", $boldFormat);
    $worksheet5->write_string("B6", "The time shift (in minutes) between the observed and predicted retention times.", $multilineFormat);
    $worksheet5->write_string("A7", "Obs. RT", $boldFormat);
    $worksheet5->write_string("B7", "The experimentally observed retention time (in minutes). It can be the corrected retention time if it comes from a peptidic quantification (with XIC extraction) or the time of the MS/MS spectrum used to identify this peptide when that is the only information available (check the 'Peptide Quantification' column in the Parameters worksheet).", $multilineFormat);
    $worksheet5->write_string("A8", "Pred. RT", $boldFormat);
    $worksheet5->write_string("B8", "The predicted retention time of the peptide (in minutes), according to the linear regression model of SSRCalc algorithm and the computed Hydrophobicity Index.", $multilineFormat);
    ### End of Info sheet

    $workbook->close();
}



####>Revision history<####
# 1.0.0 Created (VL 16/02/21)
