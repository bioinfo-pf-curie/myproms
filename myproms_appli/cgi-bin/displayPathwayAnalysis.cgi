#!/usr/local/bin/perl -w

##################################################################
# displayPathwayAnalysis.cgi    1.0.2              	 	 #
# Authors: P. Poullet, G. Arras & S. Liva (Institut Curie) 	 #
# Contact: myproms@curie.fr                               	 #
# Display in clound and table format Pathway Enrichment Analysis #
##################################################################
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
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use Spreadsheet::WriteExcel;
use Text::CSV;
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
my $pathwayID=param('ID');
my $view=param('VIEW')? param('VIEW') : 'cloud';
#my $button=($view eq 'cloud')? &tableButton : &cloudButton;
my $projectID=&promsMod::getProjectID($dbh,$pathwayID,'PATHWAY_ANALYSIS');
my $projectPath="$promsPath{pathAna}/project_$projectID";
my $pathwayPath=$projectPath."/".$pathwayID;
my $ajax = (param('AJAX'))? param('AJAX') : '';

my ($expID, $catID, $name, $desc, $param, $status, $recDate, $upDate, $user)=$dbh->selectrow_array("SELECT ID_EXPERIMENT, ID_CATEGORY, NAME, DES, PARAM_STRG, STATUS, RECORD_DATE, UPDATE_DATE, UPDATE_USER from PATHWAY_ANALYSIS where ID_PATHWAY_ANALYSIS=$pathwayID");
$desc='' unless $desc;
my ($strgFDR, $strgPVAL) = split(";", $param);
#my $FDRthreshold = (split("=",$strgFDR))[1];
my $PVALthreshold = (split("=",$strgPVAL))[1];

if ($ajax eq 'ajaxListProt') {
	&ajaxGetProteins;
	$dbh->disconnect;
	exit;
}
elsif ($ajax eq 'ajaxGetPathway') {
	&getPathwayInfo;
	$dbh->disconnect;
	exit;
}

if (param('XLS')) {
	if ( param('XLS') eq 'pathway'){
		&exportPathways;
		exit;
	}
	else {
		&exportProteins(param('pathIdent'),param('pathName'));
		exit;
	}
}

print header;warningsToBrowser(1);
print qq|
<HTML>
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
|;
&promsMod::popupInfo();
print qq|
function sortTable(sortItem) {
	window.location="./displayPathwayAnalysis.cgi?ID=$pathwayID&VIEW=table&SORT="+sortItem;
}
function sequenceView(listAnaID, protID, idType, ms) {
    var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+listAnaID+"&id_prot="+protID+"&msdata="+ms+"&id_type="+idType;
    top.openProtWindow(winLocation);
}

// AJAX --->
var XHR=null;
function ajaxGetProteins (reactNumber, type){
	var protTable=document.getElementById("protTable");
	protTable.style.display='block';
	protTable.innerHTML='<IMG src="$promsPath{images}/scrollbarGreen.gif">';

	//If XHR object already exists, the request is canceled & the object is deleted
	if(XHR && XHR.readyState != 0){
		XHR.abort();
		delete XHR;
	}
	var paramStrg="AJAX=ajaxListProt&ID=$pathwayID&VIEW=table&reactNumber="+reactNumber+"&TYPE="+type;
	//Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("POST","$promsPath{cgi}/displayPathwayAnalysis.cgi",true);
	//Send the proper header information along with the request
	XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
	XHR.setRequestHeader("Content-length", paramStrg.length);
	XHR.setRequestHeader("Connection", "close");
	XHR.onreadystatechange=function() {
		    if (XHR.readyState==4 && XHR.responseText) {
				protTable.innerHTML=XHR.responseText;
				protTable.scrollIntoView();
		    }
	}
	XHR.send(paramStrg);
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
</SCRIPT>
</HEAD>
<BODY>
<CENTER>
<FONT class="title">Display Pathway Enrichment Analysis for <FONT color="red">$name</FONT></FONT><BR><BR>|;
my (%pathwayInfo, @significantTerms);
@significantTerms=getPathwayInfo(\%pathwayInfo);#,$FDRthreshold,$PVALthreshold);
if (scalar @significantTerms) {
	my $strgButton = ($view eq "cloud")? &tableButton."&nbsp;".&graphicalButton."&nbsp;".&exportButtonPathways."<br>" : ($view eq "table")? &cloudButton."&nbsp;".&graphicalButton."&nbsp;".&exportButtonPathways."<br>" : &tableButton."&nbsp;".&cloudButton."&nbsp;".&exportButtonPathways."<br>";
	print $strgButton;
}


my %protIDs;
my $codeIdent="AC";
my ($identifierID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='$codeIdent'");
my ($strgFromSQL, $strgWhereSQL, $strgCondSQL)=($catID)? ("CATEGORY_PROTEIN C,","C.ID_CATEGORY=? and C.ID_PROTEIN=P.ID_PROTEIN and","") : ("PATHWAYANA_ANALYSIS PA,ANALYSIS_PROTEIN AP,","PA.ID_ANALYSIS = AP.ID_ANALYSIS and AP.ID_PROTEIN = P.ID_PROTEIN and", "and PA.ID_PATHWAY_ANALYSIS=? and AP.VISIBILITY>=1");
my $sthProt=$dbh->prepare("SELECT P.ID_PROTEIN, MI.VALUE  FROM $strgFromSQL PROTEIN P, MASTER_PROTEIN MP, MASTERPROT_IDENTIFIER MI
			where $strgWhereSQL
			P.ID_MASTER_PROTEIN = MP.ID_MASTER_PROTEIN and
			MP.ID_MASTER_PROTEIN = MI.ID_MASTER_PROTEIN and
			MI.ID_IDENTIFIER=? and MI.RANK=1 $strgCondSQL");
($catID)? $sthProt->execute($catID, $identifierID) : $sthProt->execute($identifierID, $pathwayID);
while(my ($protID, $uniprot) = $sthProt->fetchrow_array){
      $protIDs{$uniprot}{$protID}=1;
}
$sthProt->finish;
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my (%foundEntities, %notFoundEntities);
my $strgColor=$lightColor;

open(FOUND,"$pathwayPath/uniprotFound.csv" );
	while (my $lineFound=<FOUND>) {
		next if ($. == 1);
		chomp($lineFound);
		my($submit, $map)=split(",",$lineFound);
		$foundEntities{$submit}=1;
	}
close(FOUND);

open(NOT_FOUND,"$pathwayPath/uniprotNotFound.csv" );
	while (my $lineNotFound=<NOT_FOUND>) {
		next if ($. == 1);
		chomp($lineNotFound);
		#my($submit, $map)=split(",",$lineNotFound);
		$notFoundEntities{$lineNotFound}=1;
	}
close(NOT_FOUND);

my $totEntities = (scalar keys %foundEntities)+(scalar keys %notFoundEntities);
my $notFound = scalar(keys %notFoundEntities);
my $listNotFound=join("-", keys %notFoundEntities);
my $unannotatedProtString;
if(keys %notFoundEntities){
	$unannotatedProtString = "<FONT class=\"title3\"><A href=\"javascript:ajaxGetProteins('',\'unannotated\');\">";
	$unannotatedProtString .= $notFound."/$totEntities unannotated protein";
	$unannotatedProtString .= 's' if (scalar keys %notFoundEntities > 1);
	$unannotatedProtString .= ' (click for details)</A></FONT><BR><BR>';
}
else{
	$unannotatedProtString = '';
}
print qq|<BR><BR>
$unannotatedProtString<BR>|;
#my (%pathwayInfo, @significantTerms);
#@significantTerms=getPathwayInfo(\%pathwayInfo);#,$FDRthreshold,$PVALthreshold);
if (scalar @significantTerms) {
	if ($view eq 'cloud') {
		#print &tableButton, '&nbsp;', &graphicalButton;
print qq|
<TABLE border=0 cellspacing=0 align=center>
<TR>
    <TH align="right" class="title2">&nbsp;Enrichment factor :</TH>
    <TH bgcolor="#800080" width=100><FONT color="#FFF">&lt; 2</FONT></TH>
    <TH bgcolor="#191970" width=100><FONT color="#FFF">&lt; 4</FONT></TH>
    <TH bgcolor="#008000" width=100><FONT color="#FFF">&lt; 6</FONT></TH>
    <TH bgcolor="#CDCD00" width=100><FONT color="#FFF">&lt; 8</FONT></TH>
    <TH bgcolor="#FF8C00" width=100><FONT color="#FFF">&ge; 8</FONT></TH>
</TR>
</TABLE><br>|;

		my $log10=log 10;
		my $MAXPVAL=1e-10;
		my $LOGMIN = 0; # pval = 1
		my $minFontSize = 9;
		my $maxFontSize = 30;
		my ($minPath,$maxPath) = (sort { $pathwayInfo{$a}->[5] <=> $pathwayInfo{$b}->[5]} @significantTerms)[0,-1];
		my ($maxLogPval, $x, $y);
		if ($pathwayInfo{$minPath}->[5] < $MAXPVAL) {
			my $usedPval=($pathwayInfo{$minPath}->[5] < 1e-307)? 1e-307 : $pathwayInfo{$minPath}->[5];
			$maxLogPval = -log($usedPval)/$log10;
		} else {
			$maxLogPval = 10;
		}

		$x = ($maxFontSize - $minFontSize)/($maxLogPval - $LOGMIN);
		$y = $minFontSize - ($LOGMIN * $x);

		print qq|<DIV class="cloud">|;
		my @cloudTerms;
		foreach my $pathName (sort { $pathwayInfo{$a}->[5] <=> $pathwayInfo{$b}->[5]} keys %pathwayInfo){
			my $pval = sprintf ("%1.2e", $pathwayInfo{$pathName}->[5]);
			my $fdr = sprintf ("%.2f", $pathwayInfo{$pathName}->[6]*100);
			my $enrichment=$pathwayInfo{$pathName}->[4];
			my $reactNumber="#$pathwayInfo{$pathName}->[0]";
			my $strgSubmitEntities=$pathwayInfo{$pathName}->[12];
			my $strgMapEntities=$pathwayInfo{$pathName}->[13];
			my $fontSize = &getFontSize($pathwayInfo{$pathName}->[5], $x, $y) . 'pt';
			my $popupString ="onmouseover=\"popup('<B>Enrichment factor:</B>&nbsp;$enrichment<BR><B>p-value:</B>&nbsp;$pval<BR><B>FDR:</B>&nbsp$fdr<BR>')\" onmouseout=\"popout()\"";
			my $class=($enrichment < 2)? 'cloudFont1' : ($enrichment < 4) ? 'cloudFont2' : ($enrichment < 6) ? 'cloudFont3' : ($enrichment < 8)? 'cloudFont4' : 'cloudFont5';
			push @cloudTerms, "<A class=\"$class\" target=\"blank\" href=\"http://www.reactome.org/PathwayBrowser/$reactNumber\" $popupString><FONT style=\"font-size:$fontSize;\">$pathName</FONT></A>";
		}
		my $cloudString = join " &bull;\n", @cloudTerms;
		print $cloudString;
		print "</DIV>\n<BR><BR>";

		sub getFontSize {
			my ($pval,$x,$y) = @_;
			$pval=1e-307 if $pval < 1e-307;
			my $logPval = -log ($pval)/$log10;
			my $font = int ($logPval * $x + $y + 0.5);
			return $font;
		}
	}
	elsif ($view eq "table") {
		#print &cloudButton, '&nbsp;', &graphicalButton;
		my $strgColor=$lightColor;
		my $sort=param('SORT')? param('SORT') : 'pval';
		my ($pValColor,$enrichColor,$fdrColor)=($sort eq "pval")? ("red","black","black") : ($sort eq "ratio")? ("black","red","black") : ("black", "black","red");
		sub sortTable{
			if ($sort eq 'pval') {
				$pathwayInfo{$a}->[5] <=> $pathwayInfo{$b}->[5];
			}
			elsif ($sort eq 'ratio') {
				$pathwayInfo{$b}->[4] <=> $pathwayInfo{$a}->[4];
			}
			else {
				$pathwayInfo{$a}->[6] <=> $pathwayInfo{$b}->[6];
			}
		}
		my ($strgSort, $comp1, $comp2) = ($sort eq 'pval')? (5,$a, $b)  : ($sort eq 'ratio')? (4,$b,$a) : (6,$a,$b);
		print qq|
		<TABLE id="resultTable" border=0 cellspacing=0 align=center>
<TR bgcolor=$darkColor>
<TH class="rbBorder" rowspan=2>Pathway Term</TH>
<TH colspan=3 class="rbBorder">Proteins</TH>
<TH width=100 class="bBorder" rowspan=2><A href="javascript:sortTable('pval')"><FONT color="$pValColor">P-Value</FONT></A></TH>
</TR>
<TR bgcolor=$darkColor>
<TH width=110 class="rbBorder"><A href="javascript:sortTable('ratio')"><FONT id="ratio" color="$enrichColor">Enrichment factor</FONT></A></TH>
<TH width=110 class="rbBorder"><A href="javascript:sortTable('fdr')"><FONT id="fdr" color="$fdrColor">FDR (%)</FONT></A></TH>
<TH width=110 class="rbBorder">List</TH>
</TR>
|;
		foreach my $pathName (sort{ sortTable }  keys %pathwayInfo){
			my $pval = sprintf ("%1.2e", $pathwayInfo{$pathName}->[5]);
			my $fdr = sprintf("%.2f",$pathwayInfo{$pathName}->[6]*100);
			my $enrichment=$pathwayInfo{$pathName}->[4];
			my $reactFound=$pathwayInfo{$pathName}->[14];
			my $reactNumber="#$pathwayInfo{$pathName}->[0]";
			my $strgSubmitEntities=$pathwayInfo{$pathName}->[12];
			my $strgMapEntities=$pathwayInfo{$pathName}->[13];
			my $pathTerm="<A target=\"blank\" href=\"http://www.reactome.org/PathwayBrowser/$reactNumber\" >$pathName</A>";
			print qq|<TR bgcolor="$strgColor" valign="top">
					<TD nowrap class="title3">&nbsp;&nbsp;$pathTerm&nbsp;&nbsp;</TD>
					<TD align="center"><b>$enrichment</b></TD>
					<TD align="center">$fdr</TD>
					<TD align="center">
					<INPUT type="button"  class="font11" value="Details" onclick="ajaxGetProteins('$pathwayInfo{$pathName}->[0]','annotated');"></TD>
					<TD>&nbsp;&nbsp;$pval</TD>
				</TR>\n|;
			$strgColor=($strgColor eq $lightColor)? $darkColor : $lightColor;
		}
		print "</TABLE><BR>";
	}
	else {
		my $log10=log(10);
		my $idx=0;
		my $count=0;
		print qq
		|<SCRIPT LANGUAGE="JavaScript">
		var GP;
		function proteinLink(dSetIdx,identifier,pathIdent) {
			//ajaxProteinData(null,protID,dataSetQuantif[dSetIdx],'ajaxPepRatio');
			ajaxGetProteins(pathIdent,'annotated')
		}
		window.onload=function() {
			GP=new genericPlot({div:'mainGraphDIV',name:'GP',width:700,height:700,
						axisX:{title:'Log2(enrichment factor)'}, // ,zoomable:true
						axisY:{title:'-Log10(p-value)'}, // ,zoomable:true
						//sameScale:true,
						axisClosure:true,
						zoomable:true,
						hideDataSets:true, // hide dataset label if only 1 displayed!
						pointLabel:myPointLabel,
						pointOnclick:proteinLink,
						//pointOnList:ajaxListSelectedProteins,
						convertValue:function(axis,val){return -Math.log(val)/$log10;} // only for threshold lines
						});
			GP.addThreshold({axis:'Y',label:'max. p-value',value:$PVALthreshold,color:'#FF0000',keepAbove1:true,roundValue:1});
			GP.addDataSet(0,{name:'$name',axisX:'X',axisY:'Y',sizeRule:{min:1,max:400,ratio:1}});
			GP.addDataAsString(0,'|;
			foreach my $pathName (keys %pathwayInfo) {
				my $pathIdent= $pathwayInfo{$pathName}->[0];
				my $enrichment=log($pathwayInfo{$pathName}->[4])/log(2);
				my $pval = -log(sprintf ("%1.2e", $pathwayInfo{$pathName}->[5]))/$log10;
				#$pval=$pval/$log10;
				my $nbProt=$pathwayInfo{$pathName}->[15];
				if ($count != 0) {
					print qq|;|;
				}
				$pathName=~s/[';,]/ /g;
				print qq|$pathName,$pathIdent,$enrichment,$pval,$nbProt|;
				$count++;
				$idx++;
			}
			print qq|');|;
			#print qq|'CNOT1_HUMAN,58684,0.09,-0.21,113;FTO_HUMAN,58017,-0.05,0.26,47;FNBP4_HUMAN,62755,1.11,0.39,16;NUFP2_HUMAN,63406,1.13,0.54,27;UGPA_HUMAN,58406,-1.62,-0.20,61;CLP1L_HUMAN,58673,1.32,2.49,9;TPD52_HUMAN,65602,1.39,-2.08,21;DNJC2_HUMAN,63304,2.02,0.39,15;PHLB1_HUMAN,62793,-7.71,-2.01,5;BIG1_HUMAN,58608,-0.26,-1.45,23;PYRG2_HUMAN,60201,-0.90,-0.73,18;ARHGI_HUMAN,62850,-1.16,0.56,9;RFA1_HUMAN,59519,-0.67,-0.31,56;FND3A_HUMAN,62907,0.44,0.58,15;RT02_HUMAN,59623,1.33,-0.42,14;WFS1_HUMAN,59961,-2.24,-1.16,4;SPTA1_HUMAN,76442,-2.92,-1.68,16;HSP7C_HUMAN,58095,0.22,0.05,84;SC31A_HUMAN,62681,-0.36,-0.67,81;GT251_HUMAN,58979,0.60,-1.36,39;S100B_HUMAN,78799,1.82,0.72,24;AUXI_HUMAN,63061,-1.82,-1.32,7;AP1B1_HUMAN,58526,-0.76,0.16,48;GATB_HUMAN,58029,0.38,-1.28,6;E41L1_HUMAN,75048,-3.81,1.02,5;C43BP_HUMAN,63540,-1.14,-1.06,10;ASCC1_HUMAN,64607,1.38,-0.17,3;ATAT_HUMAN,66586,-1.21,0.91,15;SKP1_HUMAN,65205,-0.17,0.65,30;RBM19_HUMAN,62998,1.52,-0.24,9;T126A_HUMAN,66054,0.54,-0.09,14;PSMD2_HUMAN,59443,-0.08,0.28,137;SCYL2_HUMAN,63211,-0.72,0.97,22;NOTC2_HUMAN,67026,-0.34,-7.17,6;T2FA_HUMAN,63713,1.06,0.07,30;FOCAD_HUMAN,66124,0.33,-1.49,14;NSMA2_HUMAN,66364,0.63,1.29,15;GLYM_HUMAN,58035,-0.20,-1.93,68;SV2A_HUMAN,60246,-4.09,-1.63,7;EMC8_HUMAN,65773,0.34,0.28,3;GBLP_HUMAN,58935,0.70,-0.58,35;KCC2B_HUMAN,59093,3.77,1.64,8;1433F_HUMAN,58440,-0.73,-0.37,36;AT2C1_HUMAN,58581,-0.37,-0.82,21;LAMP2_HUMAN,62944,-0.67,-0.67,9;CYB5_HUMAN,65741,-0.85,-0.24,16;EST1A_HUMAN,66959,-1.15,0.20,3;CCD51_HUMAN,64626,1.65,-0.12,11;NDUS3_HUMAN,59250,0.36,-0.84,57;LASP1_HUMAN,65059,-0.09,1.67,22;MRRP1_HUMAN,64732,0.66,-0.83,40;ADHX_HUMAN,58488,-0.26,0.46,21;SDF2L_HUMAN,77656,0.04,-7.40,10;NUP85_HUMAN,59305,0.68,-0.16,40;EPHB4_HUMAN,58863,-7.08,-7.46,4;A2AP_HUMAN,66282,-2.54,-2.46,6;WRN_HUMAN,62838,2.29,-7.14,4;SDCB1_HUMAN,66670,-1.41,-1.00,13;NUD19_HUMAN,64749,1.98,-1.13,10;NDUBA_HUMAN,65463,0.23,-0.96,12;RL5_HUMAN,59564,0.72,-0.29,30;RHOA_HUMAN,59526,0.73,0.97,9;F213A_HUMAN,58890,-0.64,-0.12,16;DNPEP_HUMAN,64197,-0.42,-0.67,36;RPN1_HUMAN,59589,0.14,-0.26,82;SYTM_HUMAN,63479,0.88,0.04,28;EPS15_HUMAN,62896,-0.53,-1.66,18;TYB4_HUMAN,66886,-1.85,-0.58,16;NHEJ1_HUMAN,65112,0.28,-1.68,5;KAT1_HUMAN,64702,0.55,2.57,9;XPOT_HUMAN,59968,0.54,-0.47,41;CYC_HUMAN,58733,0.16,-0.65,38;PLXA2_HUMAN,66969,-3.88,1.18,9;LRC47_HUMAN,63378,0.04,0.48,55;1433B_HUMAN,58438,-0.52,0.37,32;NRP2_HUMAN,63169,1.77,0.90,5;TPPC3_HUMAN,66079,-0.91,-0.62,11;WDR4_HUMAN,59956,2.32,-0.97,8;P66A_HUMAN,63657,0.75,-0.89,42;2A5G_HUMAN,63523,-0.94,-0.38,16;PSMF1_HUMAN,60196,-1.08,-0.07,9;STA5B_HUMAN,59741,1.60,0.47,12;MIPEP_HUMAN,63637,0.94,-1.18,9;ANX11_HUMAN,57893,-0.89,-1.46,26;PRKRA_HUMAN,65152,-0.54,0.75,20;ES1_HUMAN,65375,0.01,-0.39,26;STAG1_HUMAN,62817,-0.42,0.58,35;ARFG3_HUMAN,63935,-1.41,-0.58,9;PPCE_HUMAN,59405,-0.27,0.15,60;RABX5_HUMAN,63857,-1.63,0.58,8;PP1R7_HUMAN,64493,-0.76,0.60,42;HTR5B_HUMAN,66962,-0.95,0.09,20;SEP11_HUMAN,58288,-0.95,-0.33,35;NUP53_HUMAN,65122,0.53,-0.09,21;ECI2_HUMAN,64653,-1.51,-0.68,8;RS6_HUMAN,59617,1.20,0.31,12;EIF2A_HUMAN,58836,0.81,-0.45,69;RHG18_HUMAN,63445,0.12,-7.66,8;PRS10_HUMAN,59422,-0.14,0.31,52;SDSL_HUMAN,65199,-2.40,2.27,8');

			print qq|GP.draw();
			//document.getElementById('mainGraphDIV').style.display='block';
		}
		function myPointLabel(dp,type) {
			var infoStrg=dp.label;
			if (type=='min') {
				return dp.label;
			}
			else {
				if (dp.dataSet.params.chart.dataSets.length > 1) {infoStrg+='\\nSet='+dp.dataSet.params.name;}
				var fc1=Math.pow(2,dp.getX());
				var fc2=Math.pow(10,dp.getY()*-1);
				infoStrg+='\\nEF='+fc1+' \\npVal= '+fc2;
			}
			infoStrg+='\\nnbProt='+dp.size;
			return infoStrg;
		}
		</SCRIPT>
		<DIV id="mainGraphDIV" style="float:left;background-color:#FFFFFF;border:solid 3px $darkColor;padding:2px;"></DIV>
		|;

	}
}
else {
	print "<FONT class=\"title3\">No significant terms</FONT><BR>";
}
print qq|<DIV id="protTable" style="float:left;display:none"></DIV>|;
print qq|

<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup();</SCRIPT>
</CENTER>
</BODY>
</HTML>
|;
sub getPathwayInfo {
	my ($refPathwayInfo)=@_;
	my $csv = Text::CSV->new({ sep_char => ',' });
	my @terms;
	open(PATHWAY,"$pathwayPath/pathway.csv" );
		while (my $line=<PATHWAY>) {
			chomp($line);
			next if ($. == 1);
			if ($csv->parse($line)) {
				my @fields=$csv->fields();
				foreach my $item(@fields) {
					my $pval=$fields[5];
					my $fdr=$fields[6];
					next if ($pval >= $PVALthreshold && $view ne 'graphical');
					#next if ($fdr >= $FDRthreshold);
					my $ratio = $fields[4];
					my $submitFoundEntities=$fields[12];
					my %submitFoundEntities;
					foreach my $uniprot (split(";",$submitFoundEntities)) {
						$submitFoundEntities{$uniprot}=1;
					}
					my $entitiesNumber=scalar(keys %submitFoundEntities);
					my $totalEntities=$fields[3];
					my $enrichment=sprintf("%.1f",($entitiesNumber/$totalEntities)/$ratio);
					next if ($enrichment <= 1);
					my $strgSubmitEntities= join("-",keys %submitFoundEntities);
					my $mapEntities=$fields[13];
					my %mapEntities;
					foreach my $mapUniprot (split(";", $mapEntities)) {
						$mapEntities{$mapUniprot}=1;
					}
					my $strgMapEntities= join("-",keys %mapEntities);
					$refPathwayInfo->{$fields[1]}=[$fields[0],$fields[1], $fields[2], $fields[3],$enrichment, $fields[5], $fields[6],$fields[7],$fields[8],$fields[9],$fields[10], $fields[11], $strgSubmitEntities, $strgMapEntities, $fields[14],$entitiesNumber];
					#0=Path Ident, 1=Path Name, 2=Entities found, 3=Entities total, 4=Ratio, 5=pVal, 6=FDR, 7=React Found, 8=React total, 9=React Ratio, 10=species ident, 11=species name, 12=Submit entities, 13=mapped entites, 14=Found react ident
					push @terms, $fields[1];
				}
			}
		}
	close(PATHWAY);
	if ($ajax eq 'ajaxGetPathway') {
		my $onchangeValue = (param('FROM') eq "cluster")? 'void' : 'ajaxGetPathwayProteinsList';
		print header(-type=>'text/plain'); warningsToBrowser(1);
		print "<SELECT id=\"pathTermsSELECT\" class=\"highlight\" onchange=\"$onchangeValue(this.value,$pathwayID,this.options[this.selectedIndex].text)\"><OPTION value=\"\">-= Select a pathway term =-</OPTION>\n";
		print "";
		my $numTerms=0;
		my $oldPathTerm="";
		foreach my $pathTerm (@terms) {
			next if ($pathTerm eq $oldPathTerm);
			$numTerms++;
			print "<OPTION value=\"$refPathwayInfo->{$pathTerm}[0]\">$pathTerm</OPTION>\n";
			$oldPathTerm=$pathTerm;
		}
		unless ($numTerms) {print "<OPTION value='' disabled>No pathway found</OPTION>\n";}
		print "</SELECT>\n";
	}
	else {
		return @terms;
	}
}
sub ajaxGetProteins {

#print header;warningsToBrowser(1);#debug
	my ($refProtInfo)=$_[0] if ($_[0]);
	my ($reactNumber,$type);
	if ($refProtInfo) {
		$reactNumber=$_[1];
		$type=$_[2];
	}
	else {
		$reactNumber=param('reactNumber');
		$type = (param('TYPE'))? param('TYPE') : "";
	}

	my $codeIdent="AC";
	my ($identifierID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='$codeIdent'");
	my ($strgFromSQL, $strgWhereSQL, $strgCondSQL)=($catID)? ("CATEGORY_PROTEIN C,","C.ID_CATEGORY=? and C.ID_PROTEIN=P.ID_PROTEIN and","") : ("PATHWAYANA_ANALYSIS PA,ANALYSIS_PROTEIN AP,","PA.ID_ANALYSIS = AP.ID_ANALYSIS and AP.ID_PROTEIN = P.ID_PROTEIN and", "and PA.ID_PATHWAY_ANALYSIS=? and AP.VISIBILITY>=1");
	my $sthProt=$dbh->prepare("SELECT P.ID_PROTEIN, MI.VALUE  FROM $strgFromSQL PROTEIN P, MASTER_PROTEIN MP, MASTERPROT_IDENTIFIER MI
				  where $strgWhereSQL
				  P.ID_MASTER_PROTEIN = MP.ID_MASTER_PROTEIN and
				  MP.ID_MASTER_PROTEIN = MI.ID_MASTER_PROTEIN and
				  MI.ID_IDENTIFIER=? and MI.RANK=1 $strgCondSQL");
	($catID)? $sthProt->execute($catID, $identifierID) : $sthProt->execute($identifierID, $pathwayID);
#my $count=0;
	while(my ($protID, $uniprot) = $sthProt->fetchrow_array) {
		next if ($protIDs{$uniprot});
#$count++;
#print "$count=$protID - $uniprot<br>\n";
		$protIDs{$uniprot}{$protID}=1;
	}
	$sthProt->finish;


	my %numEntities;
	open(FOUND,"$pathwayPath/uniprotFound.csv" );
		while (my $lineFound=<FOUND>) {
			next if ($. == 1);
			chomp($lineFound);
			my($submit, $map)=split(",",$lineFound);
			next if ($submit eq $map);
			push @{$numEntities{$submit}},$map;
		}
	close(FOUND);
	my ($pathName,$enrichment,$fdr,$pval,$strgSubmitEntities)="";

	my $csv = Text::CSV->new({ sep_char => ',' });
	open(PATHWAY,"$pathwayPath/pathway.csv" );
		while (my $line=<PATHWAY>) {
			chomp($line);
			next if ($. == 1);
			if ($csv->parse($line)) {
				my @fields=$csv->fields();
				foreach my $item(@fields) {
					if ($reactNumber eq $fields[0]) {
						$pathName=$fields[1];
						$pval=sprintf ("%1.2e", $fields[5]);
						$fdr=$fields[6];
						my $ratio = $fields[4];
						my $submitFoundEntities=$fields[12];
						my %submitFoundEntities;
						foreach my $uniprot (split(";",$submitFoundEntities)) {
							$submitFoundEntities{$uniprot}=1;
						}
						my $entitiesNumber=scalar(keys %submitFoundEntities);
						my $totalEntities=$fields[3];
						$enrichment=sprintf("%.1f",($entitiesNumber/$totalEntities)/$ratio);
						$strgSubmitEntities= join("-",keys %submitFoundEntities);
					}
				}
			}

		}
	close(PATHWAY);
	my %notFoundEntities;
	if ($type eq 'unannotated') {
		open(NOT_FOUND,"$pathwayPath/uniprotNotFound.csv" );
		while (my $lineNotFound=<NOT_FOUND>) {
			next if ($. == 1);
			chomp($lineNotFound);
			my($submit, $map)=split(",",$lineNotFound);
			$notFoundEntities{$submit}=1;
		}
		close(NOT_FOUND);
	}

	my $submitEntities = (keys %notFoundEntities)? join("-", keys %notFoundEntities) : $strgSubmitEntities;
	my $identGN="GN";
	my ($geneID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='$identGN'");
	my ($uniAcID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='AC'");
	my $selProtGene=$dbh->prepare("SELECT VALUE from MASTERPROT_IDENTIFIER where ID_MASTER_PROTEIN=? and ID_IDENTIFIER=? ORDER BY RANK");
	my $selUniprotAC=$dbh->prepare("SELECT VALUE from MASTERPROT_IDENTIFIER where ID_MASTER_PROTEIN=? and ID_IDENTIFIER=? ORDER BY RANK");
	my (%protInfo, %protGene, %protUniAC);
	foreach my $uniprotID (split('-', $submitEntities)) {
		if ($protIDs{$uniprotID}) {
			foreach my $protID (keys %{$protIDs{$uniprotID}}) {
				my ($alias, $masterID, $protDes, $species)=$dbh->selectrow_array("SELECT ALIAS, ID_MASTER_PROTEIN, PROT_DES, ORGANISM from PROTEIN where ID_PROTEIN=$protID");
				$selProtGene->execute($masterID, $geneID);
				while (my $gene = $selProtGene->fetchrow_array) {
					push @{$protGene{$protID}},$gene;
				}
				$selUniprotAC->execute($masterID,$uniAcID);
				while (my $uniprotAC=$selUniprotAC->fetchrow_array) {
					push @{$protUniAC{$protID}},$uniprotAC;
				}


				my ($anaID, $numPep)=$dbh->selectrow_array("SELECT ID_ANALYSIS, NUM_PEP FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=$protID AND VISIBILITY > 0 ORDER BY NUM_PEP LIMIT 0,1");
				my $sizeEntities = scalar(@{$numEntities{$uniprotID}}) if ($numEntities{$uniprotID});
				my $strgAlias = (!$numEntities{$uniprotID})? "<A href=\"javascript:sequenceView($anaID,$protID,'valid',0)\">$alias</A>" : "<A href=\"javascript:sequenceView($anaID,$protID,'valid',0)\" onmouseover=\"popup('<B>matches $sizeEntities entities in <br>reactome pathway database')\" onmouseout=\"popout()\">$alias</A>";
				$protInfo{$protID}=[$strgAlias, $protDes, $species, $anaID, $numPep];
				if (param('pathIdent')) {
					$refProtInfo->{$protID}=[$alias, $protDes, $species, $anaID, $numPep,@{$protGene{$protID}}[0],@{$protUniAC{$protID}}[0]];
				}
			}
		}
	}
	$selProtGene->finish;
	$dbh-> disconnect;

	if ($refProtInfo) {
		return;
	}


	my $numProt=scalar(keys %protInfo);
	my $strgNumProt = ($numProt == 1)? "Protein" : "Proteins";
	my ($lightColor,$darkColor)=&promsConfig::getRowColors;

	print header;warningsToBrowser(1);
	if (param('FROM')) {
		print join(";",keys %protInfo);
		return;
	}
	elsif ($type eq 'annotated') {
		print "<br>";
		#print &exportButton('list',$reactNumber);
		print &exportButtonProteins($reactNumber,$pathName);
		print qq|

		<br>
		<A class="title1" target=\"blank\" href=\"http://www.reactome.org/PathwayBrowser/$reactNumber\" >$pathName</A><br>
		<FONT class="title3">Enrichment Factor: $enrichment, p-value=$pval</FONT><BR><BR>
		|;
	}
	else {
		print qq|<FONT class="title1">Unannotated proteins</FONT><BR><BR>|;
	}
	print qq|
	<TABLE border=0 cellspacing=0>
		<TR bgcolor="$darkColor">
			<TH class="rbBorder" align="left">$numProt $strgNumProt</TH>
			<TH class="rbBorder">&nbsp;Uniprot&nbsp;AC&nbsp;</TH>
			<TH class="rbBorder">&nbsp;Gene&nbsp;name&nbsp;</TH>
			<TH class="rbBorder">&nbsp;Peptides&nbsp;</TH>
			<TH class="bBorder" width="700">&nbsp;Description&nbsp;-&nbsp;Species&nbsp;</TH>
		</TR>|;
		my $bgColor=$lightColor;
		foreach my $protID (sort{$protInfo{$b}->[4] <=> $protInfo{$a}->[4]} keys %protInfo) {
			my $anaID = $protInfo{$protID}->[3];
			my $geneStrg;
			my $uniACStrg;
			if (!$protGene{$protID}) {
				$geneStrg='-';
			}
			elsif (scalar @{$protGene{$protID}} > 1) {
				my $geneName1=shift(@{$protGene{$protID}});
				$geneStrg="<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-".join('<BR>&nbsp;&nbsp;-',@{$protGene{$protID}})."</FONT>')\" onmouseout=\"popout()\">$geneName1</A>";
			}
			else {
				$geneStrg=@{$protGene{$protID}}[0];
			}
			if (!$protUniAC{$protID}) {
				$uniACStrg='-';
			}
			else {
				$uniACStrg=@{$protUniAC{$protID}}[0];
			}

			print qq|
			<TR bgcolor="$bgColor">
				<TH align="left" valign="top">&nbsp;$protInfo{$protID}->[0]&nbsp;</TH>
				<TH align="left" valign="top">$uniACStrg</TH>
				<TH align="left" valign="top">$geneStrg</TH>
				<TH align="left" valign="top">$protInfo{$protID}->[4]</TH>
				<TD >$protInfo{$protID}->[1] &nbsp;&nbsp;<FONT class="org">$protInfo{$protID}->[2]</FONT></TD>
			</TR>
			|;
			$bgColor= ($bgColor eq $lightColor)? $darkColor : $lightColor;
		}
	print qq|
	</TABLE>
	</CENTER>
	|;
}
sub cloudButton{
	return qq|<INPUT type="button" class="title3" value="Cloud View" onclick="window.location='$promsPath{cgi}/displayPathwayAnalysis.cgi?VIEW=cloud&ID=$pathwayID';">|;
}
sub tableButton{
	return qq|<INPUT type="button" class="title3" value="Table View" onclick="window.location='$promsPath{cgi}/displayPathwayAnalysis.cgi?VIEW=table&ID=$pathwayID';">|;
}
sub exportButtonPathways {
	return qq|<INPUT type="button" class="title3" value="Export" onclick="window.location='$promsPath{cgi}/displayPathwayAnalysis.cgi?ID=$pathwayID&XLS=pathway';">|;
}
sub exportButtonProteins {
	my $reactNumber=$_[0];
	my $pathName=$_[1];
	return qq|<INPUT type="button" class="title3" value="Export" onclick="window.location='$promsPath{cgi}/displayPathwayAnalysis.cgi?ID=$pathwayID&XLS=proteins&pathIdent=$reactNumber&pathName=$pathName';">|;
}
sub graphicalButton {
	return qq|<INPUT type="button" class="title3" value="Graphical View" onclick="window.location='$promsPath{cgi}/displayPathwayAnalysis.cgi?VIEW=graphical&ID=$pathwayID';">|;
}
sub exportPathways {
	my $workbook = Spreadsheet::WriteExcel->new("-");
	my $worksheet = $workbook->add_worksheet();
	my %pathInfo;
	&getPathwayInfo(\%pathInfo);

	print header(-type=>"application/vnd.ms-excel",-attachment=>"$name.xls");
	#
	## Cell formatting
	my $headerMergeFormat = $workbook->add_format(align => 'center', valign => 'vcenter');
	$headerMergeFormat->set_bold();
	my $headerFormat = $workbook->add_format(align => 'center');
	$headerFormat->set_bold();
	$worksheet->set_column(0,0,30);
	$worksheet->set_column(1,2,20);

	## Writting column headers
	$worksheet->merge_range("A1:A2",'Pathway Term',$headerMergeFormat);
	$worksheet->merge_range("B1:C1",'Proteins',$headerMergeFormat);
	$worksheet->write_string(1,1,"Enrichment factor",$headerFormat);
	$worksheet->write_string(1,2,"FDR (%)",$headerFormat);
	#$worksheet->write_string(1,3,"Ratio",$headerFormat);
	$worksheet->merge_range("D1:D2",'P-value',$headerMergeFormat);
	my $formatExp = $workbook->add_format();
	$formatExp->set_num_format("0.00E+00");

	## Writting terms
	my $i = 1;
	foreach my $pathName (sort{ $pathInfo{$a}->[5] <=> $pathInfo{$b}->[5] }  keys %pathInfo){
		$i++;
		my $pval = sprintf ("%1.2e", $pathInfo{$pathName}->[5]);
		my $fdr = sprintf("%.2f",$pathInfo{$pathName}->[6]*100);
		my $enrichment=$pathInfo{$pathName}->[4];
		#my $reactFound=$pathInfo{$pathName}->[14];
		my $reactNumber="#$pathInfo{$pathName}->[0]";
		#my $strgSubmitEntities=$pathInfo{$pathName}->[12];
		#my $strgMapEntities=$pathInfo{$pathName}->[13];
		#my $pathTerm="$pathName";
		$worksheet->write_string($i,0,$pathName);
		$worksheet->write_number($i,1,$enrichment);
		$worksheet->write_number($i,2,$fdr);
		$worksheet->write_number($i,3,$pval);

	}
	$workbook->close();
	exit;
}

sub exportProteins {
	my $workbook = Spreadsheet::WriteExcel->new("-");
	my $worksheet = $workbook->add_worksheet();
	#print header;
	#print param('pathIdent');
	my %protInfo;
	my $pathIdent=param('pathIdent');
	my $pathName=param('pathName');
	$pathName=~s/ /_/g;
	&ajaxGetProteins(\%protInfo,$pathIdent,'annotated');
	print header(-type=>"application/vnd.ms-excel",-attachment=>"protein_List_For_$pathName.xls");
	my $headerMergeFormat = $workbook->add_format(align => 'center', valign => 'vcenter');
	$headerMergeFormat->set_bold();
	my $headerFormat = $workbook->add_format(align => 'center');
	$headerFormat->set_bold();

	$worksheet->write_string("A1",'Proteins',$headerMergeFormat);
	$worksheet->write_string("B1",'Uniprot AC',$headerMergeFormat);
	$worksheet->write_string("C1",'Gene Name',$headerMergeFormat);
	$worksheet->write_string("D1",'Peptide',$headerMergeFormat);
	$worksheet->write_string("E1",'Description - Species',$headerMergeFormat);
	my $i = 1;
	foreach my $protID (keys %protInfo){
		my $alias=$protInfo{$protID}->[0];
		my $protDes=$protInfo{$protID}->[1];
		my $species=$protInfo{$protID}->[2];
		my $anaID=$protInfo{$protID}->[3];
		my $numPep=$protInfo{$protID}->[4];
		my $gene=$protInfo{$protID}->[5];
		my $uniprot=$protInfo{$protID}->[6];
		$worksheet->write_string($i,0,$alias);
		$worksheet->write_string($i,1,$uniprot);
		$worksheet->write_string($i,2,$gene);
		$worksheet->write_string($i,3,$numPep);
		$worksheet->write_string($i,4,$protDes."-".$species);
		$i++;
	}
	$workbook->close();
	#print param('pathIdent');
	#foreach my $protID (keys %protInfo) {
	#	print $protInfo{$protID}->[0];
	#}
	#
	exit;
}

####>Revision history<####
# 1.0.2 add color to sort function (SL 14/03/18)
# 1.0.1 add export function and add graphical view (generic plot), comment FDR parameter and export protein list (add colimn uniprotAC) in excel format (SL 19/10/15)
# 1.0.0 New script to display Pathway enrichment analysis results in cloud and table format, replace runAndDisplayPathwayAnalysis.cgi (SL 09/03/15)