#!/usr/local/bin/perl -w

################################################################################
# displayGOAnalysis.cgi    1.2.3                                               #
# Authors: P. Poullet, G. Arras & F. Yvon (Institut Curie)                     #
# Contact: myproms@curie.fr                                                    #
# Display GO term lists and graphs of an existing GO analysis                  #
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
use Spreadsheet::WriteExcel;
use strict;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

####################
####>Parameters<####
####################
my $goAnaID=param('ID');
my $aspect=param('ASPECT');
my $projectID=&promsMod::getProjectID($dbh,$goAnaID,'GO_ANALYSIS');
my $resultDirUnix="$promsPath{go_unix}/project_$projectID/$goAnaID";
my $resultDirHtml="$promsPath{go}/project_$projectID/$goAnaID";

if(param('AJAX')){
	&getProtList(param('GOID'));
	$dbh->disconnect;
	exit;
}

if (param('XLS') && param('XLS') eq 'proteins') {
	&exportListXLS(param('goID'));
	$dbh->disconnect;
	exit;
}


my $view = (param('VIEW'))? lc(param('VIEW')): 'cloud';

my $sthParam = $dbh->prepare("SELECT NAME,PARAM_STRG FROM GO_ANALYSIS WHERE ID_GOANALYSIS=?");
$sthParam->execute($goAnaID);
my ($goAnaName,$paramStrg) = $sthParam->fetchrow_array;
$sthParam->finish;
$dbh->disconnect;
my $showAllNodes;
if($paramStrg =~ /showAllNodes\=([^;]+);/){ #deprecated parameter
	$showAllNodes = $1;
} else {
	my ($drawGraph) = ($paramStrg =~ /drawGraph\=([^;]+);/);
	$showAllNodes = ($drawGraph == 2)? 1:0;
}
my ($method) = ($paramStrg =~ /method\=([^;]+);/);

my %aspectName = ( 'P' => 'Biological Process',
                  'C' => 'Cellular Component',
                  'F' => 'Molecular Function'
                  );

#########################
####>Parsing results<####
#########################
my %pvalList;
my $countProt;
my $countPopProt;
my @significantTerms; # list containing all terms IDs that must be exclusively displayed in term list
open(RES,"$resultDirUnix/results_$aspect.txt") || die "Could not open result file: $!\n";
while(<RES>){

	# storing current line
	if(/^### (\d+)\t(\d+)/){ # 1st line in file
		$countProt = $1;
		$countPopProt = $2;
		next;
	}
	my ($goID,$goName,$pval,$numProt,$totNumProt,$protListStrg) = split /\t/;
	next if (!$numProt || $numProt !~ /\d/); # PP 07/01/16 to allow new generic plot view

	$goName = &formatGOName($goName);

	# storing all term informations in a hash table
	$pvalList{$goID}{'goName'}=$goName;
	$pvalList{$goID}{'numProt'} = $numProt;
	next if $goID eq 'unannotated' or $goID eq 'discarded';
	$pvalList{$goID}{'pval'}=$pval;
	next if $pval == 2; # means term is unsignificant but parent in graph (showAllNodes = 0) -> used in image map graph construction
	$pvalList{$goID}{'totNumProt'} = $totNumProt;
	$pvalList{$goID}{'ratio'} = (($numProt/$countProt)*100)/(($totNumProt/$countPopProt)*100);
	push @significantTerms, $goID;
}
close(RES);

#--------------------------------#
# Check for unannotated proteins #
#--------------------------------#
my $unannotatedProtString;
if(exists $pvalList{unannotated} && $pvalList{unannotated}{numProt}){
	$unannotatedProtString = "<FONT class=\"title3\"><A href=\"javascript:protTable(\'unannotated\');\">";
	$unannotatedProtString .= $pvalList{unannotated}{numProt}."/$countProt unannotated protein";
	$unannotatedProtString .= 's' if $pvalList{unannotated}{numProt} > 1;
	$unannotatedProtString .= ' (click for details)</A></FONT><BR><BR>';
}
else{
	$unannotatedProtString = '';
}

#------------------------------#
# Check for discarded proteins #
#------------------------------#
my $discardedProtString;
if(exists $pvalList{discarded} && $pvalList{discarded}{numProt}){
	$discardedProtString = "<FONT class=\"title3\"><A href=\"javascript:protTable(\'discarded\');\">";
	$discardedProtString .= $pvalList{discarded}{numProt}."/$countProt discarded protein";
	$discardedProtString .= 's' if $pvalList{discarded}{numProt} > 1;
	$discardedProtString .= ' (click for details)</A></FONT><BR><BR>';
}
else {
	$discardedProtString = '';
}

# Spread sheet case
if (param('XLS') && param('XLS') eq 'go'){
	&exportXLS(\%pvalList);
	exit;
}
#if(param('XLS')){
#	if (param('XLS') eq 'go'){
#		&exportXLS(\%pvalList);
#		exit;
#	}
#	else {
#		&exportListXLS(param('goID'));
#		exit;
#	}
#}

##############
####>HTML<####
##############
my ($light,$dark)=&promsConfig::getRowColors;
print header( -'content-encoding' => 'none');
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE>
	DIV.floatProtTableDiv
	{
		display:none;
		position:absolute;
		width:500px;
		height:350px;
		background-color:#FFFFFF;
		border: 1px solid #999999;
		padding: 10px;
	}
	DIV.closeButton
	{
		display:none;
		position:absolute;
		right:1px;
		top:1px;
		padding:0;
		border:solid 1px #999999;
		z-index:901;
	}
	.popup
	{
		width:1000px;
		height:350px;
		z-index:10;
		background-color:#FFFFFF;
		border:solid 3px #999999;
		padding:5px;
		position:absolute;
		display:none;
		box-shadow:10px 10px 20px #808080;
	}
	DIV.titleDiv
	{
		height:auto;
		padding-bottom:10px;
	}
	DIV.tableDiv
	{
		height:300px;
		overflow:auto;
	}
	A.cloudFont1
	{
		color:800080;
		text-decoration: none;
		/*white-space:nowrap;*/
	}
	A.cloudFont1:hover
	{
		color:#dd0000;
	}
	A.cloudFont2
	{
		color:191970;
		text-decoration: none;
		/*white-space:nowrap;*/
	}
	A.cloudFont2:hover
	{
		color:#dd0000;
	}
	A.cloudFont3
	{
		color:008000;
		text-decoration: none;
		/*white-space:nowrap;*/
	}
	A.cloudFont3:hover
	{
		color:#dd0000;
	}
	A.cloudFont4
	{
		color:CDCD00;
		text-decoration: none;
		/*white-space:nowrap;*/
	}
	A.cloudFont4:hover
	{
		color:#dd0000;
	}
	A.cloudFont5
	{
		color:FF8C00;
		text-decoration: none;
		/*white-space:nowrap;*/
	}
	A.cloudFont5:hover
	{
		color:#dd0000;
	}
	DIV.cloud{
		box-shadow:10px 10px 20px #808080;
		background-color:rgba(192,192,192,0.3);
		border:solid 3px #999999;
		padding:10px;
		width:80%;
		text-align:justify;
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
function floatProtTable(goID,coords){

	var coordTable = coords.split(',');

	var topRightX = parseInt(coordTable[2]);
	var topRightY = parseInt(coordTable[1]);


	var x = getDomOffset(document.getElementById('treeImg'), 'offsetLeft');
	var y = getDomOffset(document.getElementById('treeImg'), 'offsetTop');

	var protTable=document.getElementById("protTable");
	protTable.style.display='block';
	protTable.style.top=(y + topRightY)+'px';
	protTable.style.left=(x + topRightX)+'px';

	ajaxGetProteins(goID,protTable);
}
function protTable(goID){
	var protTable=document.getElementById("protTable");
	protTable.style.display='block';
	ajaxGetProteins(goID,protTable);
}
function getDomOffset( Obj, Prop ) {
    var iVal = 0;
    while (Obj && Obj.tagName != 'BODY') {
        eval('iVal += Obj.' + Prop + ';');
        Obj = Obj.offsetParent;
    }
    return iVal;
}
function closeFloatProtTable(){
	var protTable = document.getElementById('protTable');
	protTable.style.display = 'none';
}
// AJAX --->
var XHR=null;
function ajaxGetProteins(goID,protTable) {
	//wait image
	protTable.innerHTML='<IMG src="$promsPath{images}/scrollbarGreen.gif">';
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
	var URL = "$promsPath{cgi}/displayGOAnalysis.cgi?AJAX=1&ID=$goAnaID&ASPECT=$aspect&GOID="+goID;
	if(protTable.className == 'popup'){
		URL += '&FLOAT=1';
	}
	XHR.open("GET",URL,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText){
			protTable.innerHTML = XHR.responseText;
			if(protTable.className == 'popup'){
				putCloseButton(protTable);
			} else {
				protTable.scrollIntoView();
			}
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
function putCloseButton(protTable){
	var closeButton = document.getElementById('closeDiv');
	closeButton.style.display='block';
}
function sortTable(sortedElement){
	window.location="$promsPath{cgi}/displayGOAnalysis.cgi?ID=$goAnaID&ASPECT=$aspect&VIEW=table&SORT="+sortedElement;
}
function currentSort(sortedElement){
	document.getElementById(sortedElement).style.color="#DD0000";
}
|;
&promsMod::popupInfo();

if($view eq 'cloud' || $view eq 'table' || $view eq 'graphical'){
	print qq
|function sequenceView(id_protein,id_analysis){
	var winLocation="$promsPath{'cgi'}/sequenceView.cgi?id_ana="+id_analysis+"&id_prot="+id_protein+"&msdata="+top.showMSdata;
	top.openProtWindow(winLocation);
}
|;
	my $sort = (param('SORT'))?param('SORT'): 'pval';
	sub sortTable{
		if($sort eq 'pval'){
			$pvalList{$a}{$sort} <=> $pvalList{$b}{$sort}
		} else {
			$pvalList{$b}{$sort} <=> $pvalList{$a}{$sort}
		}
	}
	my $onLoadString = (scalar @significantTerms && $view eq 'table')? "onload=\"currentSort('$sort');\"": '';
	my $pvalueString = ($method eq 'bonferroni') ? 'Corrected<BR>P-value' : 'P-value';
	print qq
|</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" $onLoadString>
<CENTER>
<FONT class="title1">$aspectName{$aspect}</FONT><BR><BR>
|;
	if(scalar @significantTerms){
		if($view eq 'cloud'){
			print &tableButton, '&nbsp;', &treeButton, '&nbsp;', &graphicalButton;
			print qq
|&nbsp;<INPUT type="button" class="title3" value="Export" onclick="window.location='$promsPath{cgi}/displayGOAnalysis.cgi?XLS=go&ID=$goAnaID&ASPECT=$aspect';"><br>|;
			print $unannotatedProtString;
			print $discardedProtString;
			my $log10=log 10;
			my $MAXPVAL=1e-10;
			my $LOGMIN = 0; # pval = 1
			my $minFontSize = 9;
			my $maxFontSize = 30;

			my ($minGO,$maxGO) = (sort { $pvalList{$a}{pval} <=> $pvalList{$b}{pval}} @significantTerms)[0,-1];
			print $minGO."-".$maxGO;
			my ($maxLogPval, $x, $y);
			if ($pvalList{$minGO}{pval} < $MAXPVAL) {
				my $usedPval=($pvalList{$minGO}{pval} < 1e-307)? 1e-307 : $pvalList{$minGO}{pval};
				$maxLogPval = -log($usedPval)/$log10;
			} else {
				$maxLogPval = 10;
			}

			$x = ($maxFontSize - $minFontSize)/($maxLogPval - $LOGMIN);
			$y = $minFontSize - ($LOGMIN * $x);


			# Color code for GO colorization
			print qq
|<TABLE border=0 cellspacing=0 align=center>
<!--
<TR>
<TH align="right" class="title2">&nbsp;p-value :</TH>
|;
			#foreach my $font (30,22,18,14,9) {
			#	printf "<TD width=200 align=\"middle\" nowrap><FONT style=\"font-size:${font}pt\">&lt; %.0e ($font)</FONT></TH>\n",exp(-($font-$y)*$log10/$x);
			#}
			print qq
|</TR>
-->
<TR>
<TH align="right" class="title2">&nbsp;Enrichment factor :</TH>
<TH bgcolor="#800080" width=100><FONT color="#FFF">&lt; 2</FONT></TH>
<TH bgcolor="#191970" width=100><FONT color="#FFF">&lt; 4</FONT></TH>
<TH bgcolor="#008000" width=100><FONT color="#FFF">&lt; 6</FONT></TH>
<TH bgcolor="#CDCD00" width=100><FONT color="#FFF">&lt; 8</FONT></TH>
<TH bgcolor="#FF8C00" width=100><FONT color="#FFF">&ge; 8</FONT></TH>
</TR>
</TABLE><BR>

<DIV class="cloud">
|;
			my @cloudTerms;
			foreach my $goid (sort { $pvalList{$a}{goName} cmp $pvalList{$b}{goName}} @significantTerms){
				my $pval = sprintf ("%1.2e", $pvalList{$goid}{'pval'});
				my $ratio = sprintf ("%.2f", $pvalList{$goid}{'ratio'});
				my $fontSize = &getFontSize($pvalList{$goid}{'pval'}, $x, $y) . 'pt';
				my $popupString ="onmouseover=\"popup('<B>Enrichment factor:</B>&nbsp;$ratio<BR><B>p-value:</B>&nbsp;$pval')\" onmouseout=\"popout()\"";
				my $class=($ratio < 2)? 'cloudFont1' : ($ratio < 4) ? 'cloudFont2' : ($ratio < 6) ? 'cloudFont3' : ($ratio < 8)? 'cloudFont4' : 'cloudFont5';
				push @cloudTerms, "<A class=\"$class\" href=\"javascript:protTable('$goid');\" $popupString><FONT style=\"font-size:$fontSize;\">$pvalList{$goid}{'goName'}</FONT></A>";
			}
			my $cloudString = join " &bull;\n", @cloudTerms;

			print $cloudString;
			print "</DIV>\n<BR><BR>";

			sub getFontSize{
				my ($pval,$x,$y) = @_;
				$pval=1e-307 if $pval < 1e-307;
				my $logPval = -log ($pval)/$log10;
				#$pval = $LOGMAX if $pval > $LOGMAX;
				my $font = int ($logPval * $x + $y + 0.5);

				return $font;
			}
		}
		elsif ($view eq 'table') {
			print &cloudButton, '&nbsp;', &treeButton, '&nbsp;', &graphicalButton;
			print qq
|&nbsp;<INPUT type="button" class="title3" value="Export" onclick="window.location='$promsPath{cgi}/displayGOAnalysis.cgi?XLS=go&ID=$goAnaID&ASPECT=$aspect';">
<BR><BR>
$unannotatedProtString
$discardedProtString
<TABLE id="resultTable" border=0 cellspacing=0 align=center>
<TR bgcolor=$dark>
<TH class="rbBorder" rowspan=2>GO Term</TH>
<TH colspan=4 class="rbBorder">Proteins</TH>
<TH width=100 class="bBorder" rowspan=2><A href="javascript:sortTable('pval')"><FONT id="pval">$pvalueString</FONT></A></TH>
</TR>
<TR bgcolor=$dark>
<TH width=200 class="rbBorder">Set Frequency<BR>(total: $countProt)</TH>
<TH width=200 class="rbBorder">Background Frequency<BR>(total: $countPopProt)</TH>
<TH width=110 class="rbBorder"><A href="javascript:sortTable('ratio')"><FONT id="ratio">Enrichment factor</FONT></A></TH>
<TH width=110 class="rbBorder">List</TH>
</TR>
|;
			my $bgcolor = $light;
			foreach my $goID (sort sortTable @significantTerms){
				my $freq = sprintf("%3.1f", (($pvalList{$goID}{'numProt'}/$countProt)*100));
				my $popFreq = sprintf("%3.1f", (($pvalList{$goID}{'totNumProt'}/$countPopProt)*100));
				my $pvalF = sprintf ("%1.2e", $pvalList{$goID}{'pval'});
				my $ratio = sprintf ("%.2f", $pvalList{$goID}{'ratio'});
				print qq
|<TR bgcolor=$bgcolor valign=top>
<TD class="title3">&nbsp;<a href="http://amigo.geneontology.org/amigo/term/$goID" target=_blank>$pvalList{$goID}{'goName'}</a>&nbsp;</TD>
<TD align=center>$pvalList{$goID}{'numProt'}&nbsp;($freq% total)&nbsp;</TD>
<TD align=center>$pvalList{$goID}{'totNumProt'}&nbsp;($popFreq% total)&nbsp;</TD>
<TD align=center>$ratio</TD>
<TD align=center><INPUT type="button"  class="font11" value="Details" onclick="protTable('$goID');"></TD>
<TD align=center>$pvalF</TD>
</TR>
|;
				$bgcolor = ($bgcolor eq $light)?$dark : $light;
			}
			print '</TABLE>';
		}
		elsif ($view eq 'graphical') {
			print &cloudButton.'&nbsp;'.&treeButton.'&nbsp;'.&tableButton.'&nbsp;';
			print qq
|&nbsp;<INPUT type="button" class="title3" value="Export" onclick="window.location='$promsPath{cgi}/displayGOAnalysis.cgi?XLS=go&ID=$goAnaID&ASPECT=$aspect';"><br>|;
			my $log10=log(10);
			my $log2=log(2);
			my $idx=0;
			my $count=0;
			my $min=10000;
			my $max=0;
			my %pvalStrg;
			foreach my $goID ( keys %pvalList) {
				$min = $pvalList{$goID}{'numProt'} if ($min>$pvalList{$goID}{'numProt'});
				$max = $pvalList{$goID}{'numProt'} if ($max<$pvalList{$goID}{'numProt'});
			}
			print qq
|<SCRIPT LANGUAGE="JavaScript">
var GP;
function proteinLink(dSetIdx,identifier,goID) {
	//ajaxProteinData(null,protID,dataSetQuantif[dSetIdx],'ajaxPepRatio');
	protTable(goID)
}
window.onload=function() {
	GP=new genericPlot({div:'mainGraphDIV',name:'GP',width:700,height:700,
						axisX:{title:'Log2(enrichment factor)'}, // ,zoomable:true
						axisY:{title:'-Log10(p-value)'}, // ,zoomable:true
						axisClosure:true,
						zoomable:true,
						hideDataSets:true, // hide dataset label if only 1 displayed!
						pointLabel:myPointLabel,
						pointOnclick:proteinLink,
						convertValue:function(axis,val){return -Math.log(val)/$log10;} // only for threshold lines
					});
	GP.addThreshold({axis:'Y',label:'max. p-value',value:0.01,color:'#FF0000',keepAbove1:false,roundValue:2});
	GP.addDataSet(0,{name:'$aspectName{$aspect}',axisX:'X',axisY:'Y',sizeRule:{min:$min,max:$max,ratio:1}});
	GP.addDataAsString(0,'|;
				foreach my $goID (sort sortTable @significantTerms){
					my $freq = sprintf("%3.1f", (($pvalList{$goID}{'numProt'}/$countProt)*100));
					my $popFreq = sprintf("%3.1f", (($pvalList{$goID}{'totNumProt'}/$countPopProt)*100));
					#my $pvalF = ($pvalList{$goID}{'pval'})? -log(sprintf ("%1.2e", $pvalList{$goID}{'pval'}))/$log10 : 1e-300;
					my $pvalF = ($pvalList{$goID}{'pval'})? $pvalList{$goID}{'pval'} :  1e-300;
					$pvalStrg{$goID}=$pvalF;
					my $ratio = log($pvalList{$goID}{'ratio'})/$log2;
					my $numProt=$pvalList{$goID}{'numProt'};
					my $goName=$pvalList{$goID}{'goName'};
					$goName=~s/[';,]/ /g;
					if ($count != 0) {
						print ';';
					}
					my $pvalGP=-log(sprintf ("%1.2e", $pvalF))/$log10;
					print "$goName,$goID,$ratio,$pvalGP,$numProt";
					$count++;
				}
				print qq
|');
	GP.draw();
	//document.getElementById('mainGraphDIV').style.display='block';
}
var pvalObj=[];
|;
			foreach my $goID (keys %pvalStrg) {
				print "pvalObj['$goID']='$pvalStrg{$goID}';\n";
			}
			print qq
|function myPointLabel(dp,type) {
	var infoStrg=dp.label;
	if (type != 'min') {
		if (dp.dataSet.params.chart.dataSets.length > 1) {infoStrg+='\\nSet='+dp.dataSet.params.name;}
		var ef=(Math.pow(2,dp.getX())).toFixed(1);
		var pval=(pvalObj[dp.externalID]*1).toExponential(1);
		infoStrg+='\\nEnrich. factor='+ef+' \\np-value='+pval+'\\n# proteins='+dp.size;
	}
	return infoStrg;
}
</SCRIPT><br>
<DIV id="mainGraphDIV" style="background-color:#FFFFFF;border:solid 3px $dark;padding:2px;"></DIV>
|;

		}
	}
	else {
		print $unannotatedProtString;
		print $discardedProtString;
		print "<FONT class=\"title3\">No significant terms</FONT><BR>";
	}
	print qq
|<BR>
<DIV id="protTable" style="display:none"></DIV>
</CENTER>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup();</SCRIPT>
</BODY>
</HTML>
|;
}
elsif ($view eq 'tree'){
	print qq
|function sequenceView(id_protein,id_analysis){
	var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+id_analysis+"&id_prot="+id_protein+"&msdata="+opener.top.showMSdata;
	opener.top.openProtWindow(winLocation);
}
|;

	#---------------------------#
	# Generating image map code #
	#---------------------------#
	my $mapString = "<MAP name=\"tree\">\n";
	open(IMAP,"$resultDirUnix/imap_$aspect.tmp") or die $!;
	while(<IMAP>){
		next if /base referer/; # header
		chomp;
		my ($shape,$goID,$coords1,$coords2) = split /\s/;
		my $coords = ($shape eq 'rect')? "$coords1,$coords2" : $coords1;

		# significant node
		if(exists $pvalList{$goID}{'pval'} && $pvalList{$goID}{'pval'} != 2){
			$mapString.="<AREA shape=\"$shape\" onclick=\"floatProtTable('$goID','$coords');\" coords=\"$coords\"";
			my $pval = sprintf ("%1.2e", $pvalList{$goID}{'pval'});
			my $ratio = sprintf ("%.2f", $pvalList{$goID}{'ratio'});
			$mapString.=" onmouseover=\"popup('<B>Enrichment factor:</B>&nbsp;$ratio<BR><B>P-value:</B>&nbsp;$pval')\" onmouseout=\"popout()\">\n";
		}
		# unsignificant node (an ancestor of a significant node)
		else {
			$mapString.="<AREA shape=\"$shape\" coords=\"$coords\"";
			if($showAllNodes){
				$mapString.=" onmouseover=\"popup('Not significant term.')\" onmouseout=\"popout()\">\n";
			}else{
				$pvalList{$goID}{'goName'} =~ s/'/\\'/g; # protect quote char in javascript function
				$mapString.=" onmouseover=\"popup('$pvalList{$goID}{'goName'}')\" onmouseout=\"popout()\">\n";
			}
		}
	}
	$mapString.="</MAP>\n";
	close(IMAP);
	#--------------------#
	# Generating caption #
	#--------------------#
	### Each cell of the table has a background color given by the getPvalColor function, pval decreasing for each cell.
	my $captionString = "<TABLE border=0><TR><TD>p-value:&nbsp;</TD><TD><TABLE style=\"border:1px solid black\" cellspacing=0><TR>\n";
	my $pvalP = -2;
	my $bgcolorStrg = &convertRGBtoHTML(&convertHSVtoRGB(&getPvalColor(10**$pvalP)));
	my $colorString = "<TD bgcolor=$bgcolorStrg>10<SUP>$pvalP</SUP>&nbsp;</TD>";
	while($pvalP > -9){
		$colorString .= "<TD bgcolor=$bgcolorStrg>&nbsp;</TD>\n";
		$pvalP += -0.1;
		$bgcolorStrg = &convertRGBtoHTML(&convertHSVtoRGB(&getPvalColor(10**$pvalP)));
	}
	$colorString.="<TD bgcolor=$bgcolorStrg>10<SUP>-9</SUP>&nbsp;</TD>";
	$captionString .= $colorString."\n</TR></TABLE></TD></TR></TABLE>\n";

	print qq
|</SCRIPT>
<TITLE>$aspectName{$aspect} Tree View</TITLE>
</HEAD>
<BODY>
<CENTER>
<FONT class="title1">$aspectName{$aspect}</FONT><BR><BR>
$captionString<BR>
$mapString
<IMG name="treeImg" id="treeImg" usemap="#tree" src="$resultDirHtml/graph_$aspect.png"><BR><BR>
<DIV id="protTable" class="popup">
</DIV>
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


#########################################################
### Display protein list of a selected GO term (AJAX) ###
#########################################################
sub getProtList{

	my $goID = $_[0];
	my ($refProtList, $refMasterProtList)= @_[1,2];
	#print header;

	#my $refProtList=$_[1];
	#my $refMasterProtList=$_[2];
	my (%protInfo,%masterProteins);
	my $goName;
	my $sthProt = $dbh->prepare("SELECT ALIAS,ID_MASTER_PROTEIN,PROT_DES,ORGANISM FROM PROTEIN WHERE ID_PROTEIN=?");
	my ($uniAcID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='AC'");
	my $selUniprotAC=$dbh->prepare("SELECT VALUE from MASTERPROT_IDENTIFIER where ID_MASTER_PROTEIN=? and ID_IDENTIFIER=? ORDER BY RANK");
	# Preparing fetching analysis ID where each protein is detected by the maximum peptide number, to display sequenceView.
	# Changed on the 16/07/2014 because the experiment context was not considered -> The project best analysis was chosen and not the best analysis of the experiment (misleading in pepnum)!!!
	#$sthAna = $dbh->prepare("SELECT ANALYSIS_PROTEIN.ID_ANALYSIS,NUM_PEP
	#					FROM ANALYSIS_PROTEIN
	#					WHERE ANALYSIS_PROTEIN.ID_PROTEIN=?
	#					ORDER BY VISIBILITY DESC, NUM_PEP DESC LIMIT 0,1");
	my $sthAna = $dbh->prepare("SELECT AP.ID_ANALYSIS,NUM_PEP
					FROM GO_ANALYSIS GA
					INNER JOIN SAMPLE S ON GA.ID_EXPERIMENT=S.ID_EXPERIMENT
					INNER JOIN ANALYSIS A ON S.ID_SAMPLE=A.ID_SAMPLE
					INNER JOIN ANALYSIS_PROTEIN AP ON A.ID_ANALYSIS=AP.ID_ANALYSIS AND AP.ID_PROTEIN=? AND VISIBILITY > 0
					WHERE ID_GOANALYSIS=$goAnaID
					ORDER BY NUM_PEP DESC LIMIT 0,1");
	#$sthAna = $dbh->prepare("SELECT ANALYSIS_PROTEIN.ID_ANALYSIS,NUM_PEP
	#					FROM ANALYSIS_PROTEIN,GOANA_ANALYSIS
	#					WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=GOANA_ANALYSIS.ID_ANALYSIS
	#					AND ANALYSIS_PROTEIN.ID_PROTEIN=? AND VISIBILITY > 0
	#					AND ID_GOANALYSIS=?
	#					ORDER BY NUM_PEP DESC LIMIT 0,1;");

	# Fetching proteins IDs of corresponding GO ID in result file
	my ($ratio,$pvalue,$countProt,$countPopProt);
	open(RES,"$resultDirUnix/results_$aspect.txt") || die "Could not open result file";
	while(<RES>){
		if(/^### (\d+)\s(\d+)/){ # 1st line in file
			$countProt = $1;
			$countPopProt = $2;
			next;
		}
		my @info = split(/\t/, $_);
		if($info[0] eq $goID){
			if ($info[0] eq 'unannotated'){
				$goName = 'Unannotated proteins' ;
			} elsif ($info[0] eq 'discarded') {
				$goName = 'Discarded proteins';
			} else {
				$goName = &formatGOName($info[1]);
				$ratio = sprintf("%.2f",(($info[3]/$countProt)*100)/(($info[4]/$countPopProt)*100));
				$pvalue = sprintf("%1.2e",$info[2]);
			}
			foreach my $protID (split ';',$info[5]){
				$sthProt->execute($protID);
				($protInfo{$protID}{'Name'},my $masterProtID,$protInfo{$protID}{'Des'},$protInfo{$protID}{'Org'}) = $sthProt->fetchrow_array;

				$sthAna->execute($protID);
				($protInfo{$protID}{'AnaID'},$protInfo{$protID}{'NumPep'}) = $sthAna->fetchrow_array;
				if ($masterProtID) {
					@{$masterProteins{$masterProtID}{gene}}=();
					@{$masterProteins{$masterProtID}{uniAC}}=();
					$protInfo{$protID}{'masterProt'}=$masterProtID;
				}
			}
			last;
		}
	}
	close(RES);

	$sthProt->finish;
	$sthAna->finish;


	####>Fetching master proteins info<####
	my ($geneNameID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
	my $sthMP=$dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN=? AND ID_IDENTIFIER=$geneNameID ORDER BY RANK");
	foreach my $masterProtID (keys %masterProteins) {
		$sthMP->execute($masterProtID);
		while (my ($gene)=$sthMP->fetchrow_array) {
			push @{$masterProteins{$masterProtID}{gene}},$gene;
		}
		$selUniprotAC->execute($masterProtID,$uniAcID);
		while (my $uniprotAC=$selUniprotAC->fetchrow_array) {
			push @{$masterProteins{$masterProtID}{uniAC}},$uniprotAC;
		}
	}
	$sthMP->finish;

	if ($refProtList) {
		%{$refProtList}=%protInfo;
		%{$refMasterProtList}=%masterProteins;
		$dbh->disconnect;
		return;
	}



	####>HTML<####
	my ($light,$dark)=&promsConfig::getRowColors;
	print header(-'content-encoding'=>'no',-charset=>'UTF-8',-type=>'text/plain');
	warningsToBrowser(1);
	my ($titleClass,$tableClass);
	if(param('FLOAT')){
		($titleClass,$tableClass) = ('titleDiv','tableDiv');
	}
	else {
		($titleClass,$tableClass) = ('','');
	}
	print qq
|<DIV id="titleDiv" class="$titleClass">
<BR>
|;
	print &exportList($goID);
	print qq
|<BR><BR>
<CENTER>
	<FONT class="title1"><a href="http://amigo.geneontology.org/amigo/term/$goID" target=_blank>$goName</a></FONT><BR>
|;
	print qq|<FONT class="title3">Enrichment Factor: $ratio, p-value=$pvalue</FONT>| unless ($goID eq 'unannotated') or ($goID eq 'discarded');
	print "<BR><BR>" unless (param('FLOAT'));
	print qq
|</CENTER>
	<DIV id="closeDiv" class="closeButton">
	<A href="javascript:closeFloatProtTable();"><B>X</B></A>
	</DIV>
</DIV>
<DIV id="tableDiv" class="$tableClass"><CENTER><TABLE border=0 cellspacing=0>
|;
	print "<TR bgcolor=$dark><TH class=\"rbBorder\" align=\"left\">&nbsp;",scalar keys %protInfo," Proteins&nbsp;</TH><TH class=\"rbBorder\">&nbsp;Uniprot&nbsp;AC&nbsp;</TH><TH class=\"rbBorder\">&nbsp;Gene&nbsp;name&nbsp;</TH><TH class=\"rbBorder\">&nbsp;<FONT style=\"cursor:help;\" onmouseover=\"popup('Peptides in best Analysis');\" onmouseout=\"popout();\">Peptides</FONT>&nbsp;</TH><TH class=\"bBorder\" width=700>Description - Species</TH></TR>\n";
	my $bgcolor = $light;
	foreach my $protID (sort { $protInfo{$b}{'NumPep'} <=> $protInfo{$a}{'NumPep'} } keys %protInfo){
		my $geneStrg='-'; # default
		my $uniACStrg='-';
		if ($protInfo{$protID}{masterProt} && $masterProteins{$protInfo{$protID}{masterProt}}{gene}[0]) {
			if (scalar @{$masterProteins{$protInfo{$protID}{masterProt}}{gene}} > 1) {
				my $geneName1=shift(@{$masterProteins{$protInfo{$protID}{masterProt}}{gene}});
				$geneStrg="<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-".join('<BR>&nbsp;&nbsp;-',@{$masterProteins{$protInfo{$protID}{masterProt}}{gene}})."</FONT>')\" onmouseout=\"popout()\">$geneName1</A>";
			}
			else {$geneStrg=$masterProteins{$protInfo{$protID}{masterProt}}{gene}[0];}
		}
		if ($protInfo{$protID}{masterProt} && $masterProteins{$protInfo{$protID}{masterProt}}{uniAC}[0]) {
			my $uniAC=@{$masterProteins{$protInfo{$protID}{masterProt}}{uniAC}}[0];
			$uniACStrg=$uniAC;
		}
		print qq
|<TR bgcolor=$bgcolor><TH align="left" valign="top">&nbsp;<A href="javascript:sequenceView($protID,$protInfo{$protID}{AnaID})">$protInfo{$protID}{Name}</A>&nbsp;</TH>
<TH align="left" valign="top">&nbsp;$uniACStrg&nbsp;</TH>
<TH align="left" valign="top">&nbsp;$geneStrg&nbsp;</TH>
<TH valign="top">&nbsp;$protInfo{$protID}{NumPep}</TH>
<TD>&nbsp;$protInfo{$protID}{Des} <FONT class="org">$protInfo{$protID}{Org}</FONT></TD></TD>
</TR>
|;
		$bgcolor = ($bgcolor eq $light)? $dark : $light;
	}
	print "</TABLE><BR><BR></CENTER></DIV></BODY></HTML>\n";
}

sub treeButton{
	return qq
|<INPUT type="button" class="title3" value="Tree View" onclick="window.open('$promsPath{cgi}/displayGOAnalysis.cgi?VIEW=tree&ID=$goAnaID&ASPECT=$aspect','tree','menubar=no, scrollbars=yes, width=1000, height=750, resizable=yes');">
|;
}

sub cloudButton{
	return qq
|<INPUT type="button" class="title3" value="Cloud View" onclick="window.location='$promsPath{cgi}/displayGOAnalysis.cgi?VIEW=cloud&ID=$goAnaID&ASPECT=$aspect';">
|;
}

sub tableButton{
	return qq
|<INPUT type="button" class="title3" value="Table View" onclick="window.location='$promsPath{cgi}/displayGOAnalysis.cgi?VIEW=table&ID=$goAnaID&ASPECT=$aspect';">
|;
}

sub graphicalButton {
	return qq
|<INPUT type="button" class="title3" value="Graphical View" onclick="window.location='$promsPath{cgi}/displayGOAnalysis.cgi?VIEW=graphical&ID=$goAnaID&ASPECT=$aspect';">
|;
}
sub exportList {
	my $goID=$_[0];
	return qq
|<INPUT type="button" class="title3" value="Export List" onclick="window.location='$promsPath{cgi}/displayGOAnalysis.cgi?XLS=proteins&goID=$goID&VIEW=table&ID=$goAnaID&ASPECT=$aspect';">
|;
}

sub exportXLS{
	my %pvalList = %{$_[0]};
	my $workbook = Spreadsheet::WriteExcel->new("-");
	my $worksheet = $workbook->add_worksheet();
	print header(-type=>"application/vnd.ms-excel",-attachment=>"$goAnaName ($aspectName{$aspect}).xls");

	# Cell formatting
	my $headerMergeFormat = $workbook->add_format(align => 'center', valign => 'vcenter');
	$headerMergeFormat->set_bold();
	my $headerFormat = $workbook->add_format(align => 'center');
	$headerFormat->set_bold();
	$worksheet->set_column(0,0,30);
	$worksheet->set_column(1,2,20);

	# Writting column headers
	$worksheet->merge_range("A1:A2",'GO Term',$headerMergeFormat);
	$worksheet->merge_range("B1:D1",'Proteins',$headerMergeFormat);
	$worksheet->write_string(1,1,"Background Frequency",$headerFormat);
	$worksheet->write_string(1,2,"Set Frequency",$headerFormat);
	$worksheet->write_string(1,3,"Ratio",$headerFormat);
	$worksheet->merge_range("E1:E2",'P-value',$headerMergeFormat);
	my $formatExp = $workbook->add_format();
	$formatExp->set_num_format("0.00E+00");

	# Writting terms
	my $i = 1;
	foreach my $goID (sort {$pvalList{$a}{'pval'} <=> $pvalList{$b}{'pval'}} keys %pvalList){
		next if ($pvalList{$goID}{'pval'} == 2 or $goID eq 'unannotated' or $goID eq 'discarded');
		$i++;
		my $freq = sprintf("%3.1f", (($pvalList{$goID}{'numProt'}/$countProt)*100));
		my $popFreq = sprintf("%3.1f", (($pvalList{$goID}{'totNumProt'}/$countPopProt)*100));
		my $pvalF = sprintf ("%1.2e", $pvalList{$goID}{'pval'});
		my $ratio = sprintf ("%.3f", $pvalList{$goID}{'ratio'});
		$worksheet->write_string($i,0,$pvalList{$goID}{'goName'});
		$worksheet->write_number($i,1,$pvalList{$goID}{'totNumProt'});
		$worksheet->write_number($i,2,$pvalList{$goID}{'numProt'});
		$worksheet->write_number($i,3,$ratio);
		$worksheet->write_number($i,4,$pvalF,$formatExp);
	}
	$workbook->close();
	exit;
}

sub exportListXLS {

	my $goID=$_[0];
	my (%protList,%masterProtList);
	print header(-type=>"application/vnd.ms-excel",-attachment=>"$goID.xls");
	my $workbook = Spreadsheet::WriteExcel->new("-");
	my $worksheet = $workbook->add_worksheet();

	&getProtList($goID,\%protList,\%masterProtList);

	# Cell formatting
	my $headerMergeFormat = $workbook->add_format(align => 'center', valign => 'vcenter');
	$headerMergeFormat->set_bold();
	my $headerFormat = $workbook->add_format(align => 'center');
	$headerFormat->set_bold();

	# Writting column headers
	$worksheet->write_string("A1",'Proteins',$headerMergeFormat);
	$worksheet->write_string("B1",'UniprotAC',$headerMergeFormat);
	$worksheet->write_string("C1",'Gene Name',$headerMergeFormat);
	$worksheet->write_string("D1","Peptides",$headerFormat);
	$worksheet->write_string("E1","Description - Species",$headerFormat);
	my $i=1;
	foreach my $protID (sort { $protList{$b}{'NumPep'} <=> $protList{$a}{'NumPep'} } keys %protList){
		my $geneStrg='-';
		my $uniACStrg='-';
		if ($protList{$protID}{masterProt} && $masterProtList{$protList{$protID}{masterProt}}{gene}[0]) {
			if (scalar @{$masterProtList{$protList{$protID}{masterProt}}{gene}} > 1) {
				my $geneName1=shift(@{$masterProtList{$protList{$protID}{masterProt}}{gene}});
				$geneStrg=$geneName1;
			}
		}
		if ($protList{$protID}{masterProt} && $masterProtList{$protList{$protID}{masterProt}}{uniAC}[0]) {
			my $uniAC=@{$masterProtList{$protList{$protID}{masterProt}}{uniAC}}[0];
			$uniACStrg=$uniAC;
		}

		#$worksheet->write_number($i,0,$protList{$protID}{'Name'});
		$worksheet->write_string($i,0,$protList{$protID}{Name});
		$worksheet->write_string($i,1,$uniACStrg);
		$worksheet->write_string($i,2,$geneStrg);
		$worksheet->write_number($i,3,$protList{$protID}{NumPep});
		$worksheet->write_string($i,4,$protList{$protID}{Des}."-".$protList{$protID}{Org});
		$i++;
	}
	$workbook->close();
	exit;
}

sub formatGOName{
	my $goName = shift;

	# applying uppercase to term first letter, unless this letter must be lowercase (eg: mRNA, ncRNA...)
	unless($goName =~ /^\w{1,3}RNA/){
		$goName = ucfirst($goName);
	}

	return $goName;
}

#########################
# Graph color functions #
#########################

sub getPvalColor{
    # HSV gradient (Hue,Saturation,Lightness)
    my $pval = $_[0];
    my $logPval = -(log $pval)/(log 10);
    if($logPval >= 9){ # < 10^-9 -> red
        return "0, 1, 1";
    }
	elsif($logPval >= 5){ #10^-5 to 10^-9 -> yellow to red
        my $hue = ($logPval*(-0.0375))+0.3375;
        return "$hue, 1, 1";
    }
	elsif($logPval >= 2) {#10^-2 to 10^-5 -> white to yellow
        my $saturation = ($logPval*(1/3))-(2/3);
        return "0.15, $saturation, 1";
    }
	else {
        return "0, 0, 1";
    }
}

sub convertRGBtoHTML{
	my @RGB = @_;
	my @hexa = (0..9,'A','B','C','D','E','F');
	my $htmlColor = "#";
	foreach (@RGB){
		$_ *= 255;
		my $p1 = int($_/16);
		my $p0 = $_ % 16;
		my $hexaP1 = $hexa[$p1];
		my $hexaP0 = $hexa[$p0];
		$htmlColor.="$hexaP1$hexaP0";
	}
	return $htmlColor;
}

sub convertHSVtoRGB{
	my($h,$s,$v) = ($_[0] =~ /(.+)\, (.+)\, (\d)/);
	$h *= 360;
	my ($f,$p,$q,$t);
	if( $s == 0 ) {
		# achromatic (grey)
		return ($v,$v,$v);
	}
	else {
		$h /= 60;
		my $i = int($h);
		$f = $h - $i;
		$p = $v*(1-$s);
		$q = $v*(1-$s*$f);
		$t = $v * (1-$s*(1-$f));
		my ($r,$g,$b) = ($i == 0)?($v,$t,$p):($i == 1)?($v,$t,$p):($i == 2)?($v,$t,$p):($i == 3)?($v,$t,$p):($v,$t,$p);
		return ($r,$g,$b);
	}
}

####>Revision history<####
# 1.2.3 Updated link to AMIGO (PP 17/07/18)
# 1.2.2 Filter result file for missing numProt to prevent failure of new scatter plot view (PP 07/01/16)
# 1.2.1 replace graph by tree, add generic plot view and export List proteins (add column uniprotAC) (SL 19/10/15)
# 1.2.0 No longer uses category parameter for AJAX query (PP 18/07/14)
# 1.1.9 Minor bug correction in peptide number (GA 16/07/14)
# 1.1.8 Change in color code legende (PP 16/07/14)
# 1.1.7 Minor modification in Enrichment Factor table in cloud view (GA 16/06/14)
# 1.1.6 Update cloud view with colors and correct log error for too small p-values (GA 16/06/14)
# 1.1.5 Added Gene name to protein list (PP 29/05/13)
# 1.1.4 Bug fix when sorting table (FY 22/03/13)
# 1.1.3 Compatible with new data path .../results/project_ID/GoID (PP 13/09/12)
# 1.1.2 Adding notification for discarded proteins (absent from background)<BR>& Inverting background and set frequencies in table view (FY 17/04/12)
# 1.1.1 Changing font size management for cloud view (FY 16/04/12)
# 1.1.0 Adding a cloud style view for enriched terms +<BR>Managing view if no significant terms (FY 05/04/12)
# 1.0.1 Floating div for protein lists in graphical view (FY 27/01/12)
# 1.0.0 New script to display GO graphs and enrichment results (FY 08/07/11)