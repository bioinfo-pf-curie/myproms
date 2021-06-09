#!/usr/local/bin/perl -w

#############################################################################
# startGOQuantiAnalysis.cgi    1.3.2                                        #
# Authors: P. Poullet, G. Arras, F. Yvon & S. Liva (Institut Curie)         #
# Contact: myproms@curie.fr                                                 #
# Setting heatmap GO enrichment test parameters and performing it           #
#############################################################################
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

use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use POSIX qw(strftime);
use promsConfig;
use promsMod;
use promsQuantif;
use strict;
use goAnalysis;
use promsOboParser;
use GO::AnnotationProvider::AnnotationParser;
#use List::Util qw(sum);
use File::Copy;
use Encode qw(encode_utf8); # ... Encode needed if hard UTF8 chars in code!!!
use File::Path qw(rmtree);

my %promsPath = &promsConfig::getServerInfo;
my ($lightColor,$darkColor)=&promsConfig::getRowColors;

#---------------------------------#
# Fetching and checking form data #
#---------------------------------#
#my @binValues = (param('bin0'))? (param('bin0'),param('bin1'),param('bin2'),param('bin3'),param('bin4')): (5,15,60,15,5);
#die 'Incorrect bin values' unless sum(@binValues) == 100;

my @thrLogRatios = (defined param('logRatio0'))? (param('logRatio0'),param('logRatio1'),param('logRatio2'),param('logRatio3')):(-2,-1,1,2);
for (my $i=0;$i<$#thrLogRatios;$i++) {
    if ($thrLogRatios[$i] >= $thrLogRatios[$i+1]) {
        die 'Incorrect log ratio threshold values';
    }
}

#my $anaID = param('anaList') or die "No analysis";
die "No selected quantification" unless param("quanti");
my $quantiID = param("quanti");
my $ratioPos = param('ratio') // 1;

#---------------------------------#
# Fetching protein ids and ratios #
#---------------------------------#
my $dbh = &promsConfig::dbConnect;
my $minPep=param('minPep') || 1;
my $maxPvalue=param('ratioPvalue') || 1;
my $numPepCode='NUM_PEP_USED'; # could be 'DIST_PEP_USED'

##my $sthQuanti = $dbh->prepare("SELECT NAME,ID_MODIFICATION,Q.ID_QUANTIFICATION FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
my $sthQuanti = $dbh->prepare("SELECT Q.NAME,Q.ID_MODIFICATION,GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ','),Q.ID_QUANTIFICATION FROM QUANTIFICATION Q
								LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
								WHERE Q.ID_QUANTIFICATION=?
								GROUP BY Q.ID_QUANTIFICATION");
$sthQuanti->execute($quantiID);
my ($quantiName,$modifQuantifID,$multiModifStrg) = $sthQuanti->fetchrow_array;
my $isModifQuantif=($modifQuantifID || $multiModifStrg)? 1 : 0; 
$sthQuanti->finish;
my $proteinStrg=($isModifQuantif)? 'site' : 'protein';
my $upProteinStrg=ucfirst($proteinStrg);

my $tmpDir;
my $formView=0;
unless (param('tmpDir')) { # 1st call
	$formView=1;

	##>Clean older dirs if any
	&promsMod::cleanDirectory("$promsPath{tmp}/GO",'1h'); # 1 hour
	mkdir "$promsPath{tmp}/GO" unless -e "$promsPath{tmp}/GO";

	$tmpDir = strftime("%Y%m%d%H%M%S",localtime);
    mkdir "$promsPath{tmp}/GO/$tmpDir";

	#------#
	# HTML #
	#------#
	print header(-'content-encoding' => 'no',-charset=>'utf-8');
	warningsToBrowser(1);

	#my $image = &densityPlot(\@ratios,\@thrLogRatios);
	my $helpMessage = ucfirst($proteinStrg).'s are separated into 5 bins based on their quantitative ratios.<br>1 GO enrichment test is performed for each bin.';
	print qq
|<HTML>
<HEAD>
<TITLE>Gene Ontology Enrichment Analysis On Quantification Data</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT type="text/javascript">
|;
	&promsMod::popupInfo();
	print qq
|function checkForm(myForm){
    if (!myForm.anaName.value) {
        alert('Enter a name for the analysis.');
        return false;
    }
	else if (!myForm.ontology.value) {
        alert('Select an ontology file.');
        return false;
    }
	else if (!myForm.annotation.value) {
        alert('Select a species annotation.');
        return false;
    }
    else { // check bin size
		for (var i = 0;i<=4;i++) {
			var binValue = 1*(document.getElementById('bin'+i).innerHTML);
			if((binValue/100) * countRatios < 1){
				alert('Bin '+(1*i+1)+' is empty');
				return false;
			}
		}
	}
	document.enrichForm.start.disabled=true; // disable submit button to prevent multiple submit
    return true;
}
function selectStatistics(){
    var myForm = document.enrichForm;
    var value = getRadioValue(myForm.criteria);
    if(value == 'FDR'){
		myForm.pval.disabled = true;
		myForm.bonferroni.disabled = true;
		myForm.FDR.disabled = false;
		myForm.FDRmethod.disabled = false;
    }
	else if (value == 'pval'){
		myForm.FDR.disabled = true;
		myForm.FDRmethod.disabled = true;
		myForm.pval.disabled = false;
		myForm.bonferroni.disabled = false;
    }
}
function getRadioValue(radio){
    var value;
    var len = radio.length;
    for(var i = 0;i<len;i++){
		if(radio[i].checked){
			value = radio[i].value;
		}
    }
    return value;
}
// AJAX --->
var XHR=null;
function refreshPlot() {
	var plotDiv = document.getElementById('plotDiv');
	var img = document.getElementById('imgPlot');
	plotDiv.style.display='none';

	var scrollbarDiv=document.getElementById('scrollbarDiv');
	scrollbarDiv.style.display='';

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
	var GETstring = 'AJAX=1&quanti=$quantiID&tmpDir=$tmpDir&ratio='+document.getElementById('ratio').value;
	for(var i=0;i<4;i++){
		GETstring += '&logRatio'+i+'='+document.getElementById('logRatio'+i).value;
	}
	if(document.enrichForm.ratioPvalue.value){
		GETstring += '&ratioPvalue='+document.enrichForm.ratioPvalue.value;
	}
	if(document.enrichForm.minPep.disabled == false){
		GETstring += '&minPep='+document.enrichForm.minPep.value;
	}
	var URL = "$promsPath{cgi}/startGOQuantiAnalysis.cgi?"+GETstring;
	XHR.open("GET",URL,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText){
			img.src = XHR.responseText + '?' + new Date().getTime();
			scrollbarDiv.style.display='none';
			plotDiv.style.display='';
			document.getElementById('updateBUT').disabled=true; // No need to update chart again
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
var logRatios = [];
var countRatios = 0;
var totTrueProteins = 0;
var binSpans = [];
var inputLogRatios = [];
function initAll() {
    buildLogRatios(); // updates countRatios & totTrueProteins
    for (let i=0;i<=4;i++) {
        binSpans.push(document.getElementById('bin'+i));
    }
    for (let i=0;i<=3;i++) {
        inputLogRatios.push(document.getElementById('logRatio'+i));
    }
    updateBins();
}
function updateLogRatios(j) {
	var normRatio=document.getElementById('normRatio'+j).value;
	if (normRatio.match('^1/')) {normRatio=eval(normRatio);}
	if (!isNaN(normRatio*1)) {
		document.getElementById('logRatio'+j).value=1*((Math.log(normRatio)/Math.log(2)).toFixed(2));
	}
	updateBins(j);
}
function updateBins(j) { // j optional
    if (j !== undefined) {
		document.getElementById('updateBUT').disabled=false; // Allow update of chart
		correctRatios(j);
		/* Update all normal ratio values (to be safe) */
		for (let i=0;i<=3;i++) {
			let normRatioStrg=(inputLogRatios[i].value < 0)? '1/'+((1/(2**inputLogRatios[i].value)).toFixed(2))*1 : ((2**inputLogRatios[i].value).toFixed(2))*1;
			document.getElementById('normRatio'+i).value=normRatioStrg;
		}
	}
	else {
		document.getElementById('updateBUT').disabled=true; // No need to update chart
	}
    var binSize =[0,0,0,0,0];
	var trueProtBins=[{},{},{},{},{}]; // used for true protein count if PTM quantif
    for (let i=0;i<countRatios;i++) {
		let curLogRatio=quantificationData[logRatios[i]][0];
        var matchedBinIdx=4; // default
		for (let j=0; j<=3; j++) {
			if (curLogRatio <= inputLogRatios[j].value) {
				matchedBinIdx=j;
				break;
			}
		}
		binSize[matchedBinIdx]++;
		if ($isModifQuantif) {
			trueProtBins[matchedBinIdx][quantificationData[logRatios[i]][3]]=1; // array of hash of proteinIDs
		}
    }
    for (let i=0;i<=4;i++) {
        binSpans[i].innerHTML = 1*(((binSize[i] / countRatios) * 100).toFixed(1));
    }
	if ($isModifQuantif) {
		for (let i=0;i<=4;i++) {
			var binSizeProt=Object.keys(trueProtBins[i]).length;
			document.getElementById('binTrueProt'+i).innerHTML = 1*(((binSizeProt / totTrueProteins) * 100).toFixed(1)); // sum of bins > 100% due to same prot in multiple bins
		}

	}
}
function correctRatios(j) { // to make sure previous ratio threshold is always < newly set ratio
    for (var i=0;i<j;i++) {
        if (inputLogRatios[i].value*1 >= inputLogRatios[j].value*1) {
            inputLogRatios[i].value = inputLogRatios[j].value*1 + (i-j) * 0.1;
        }
    }
    for (var k=j+1;k<4;k++) {
        if (inputLogRatios[k].value*1 <= inputLogRatios[j].value*1) {
            inputLogRatios[k].value = inputLogRatios[j].value*1 + (k-j) * 0.1;
        }
    }
}
function buildLogRatios() {
    logRatios = [];
	var trueProteins={};
    var pValue = document.enrichForm.ratioPvalue.value;
    var minPepProt = document.enrichForm.minPep.value;
    if (pValue < 1 \|\| minPepProt > 1) {
        for (var i=0; i<quantificationData.length; i++) {
            if (quantificationData[i][1] !== undefined && quantificationData[i][1] < pValue && quantificationData[i][2] >= minPepProt) {
                logRatios.push(i);
				if ($isModifQuantif) {trueProteins[quantificationData[i][3]]=1;} // records list of proteinIDs
            }
        }
    }
	else {
		for (var i=0;i<quantificationData.length;i++) {
			logRatios.push(i);
			if ($isModifQuantif) {trueProteins[quantificationData[i][3]]=1;} // records list of proteinIDs
		}
    }
    countRatios = logRatios.length;
	document.getElementById('numProtSPAN').innerHTML=countRatios;
	if ($isModifQuantif) {
		totTrueProteins=Object.keys(trueProteins).length;
		document.getElementById('numTrueProtSPAN').innerHTML=totTrueProteins;
	}
}
function updateThreshold() {
    buildLogRatios();
    updateBins();
    refreshPlot();
}
function updateRatio(ratioPos){
    window.location="$promsPath{cgi}/startGOQuantiAnalysis.cgi?quanti=$quantiID&ratio="+ratioPos;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" onload="initAll();">
<CENTER>
<FONT class="title">Gene Ontology Enrichment Analysis On Quantification Data</FONT>&nbsp;<FONT class="help" onmouseover="popup('$helpMessage');" onmouseout="popout();">[?]</FONT><BR>
<DIV id="waitDIV">
<BR><BR><BR><FONT class="title3">Fetching data. Please wait...</FONT><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">
</DIV>
|;

}


################################>>PP 19/08/15
# 1.2.0 Attempt for multi-quantif GO analysis (PP ../08/15)
my (%quantifInfo,%quantifValues);
my $quantif=$quantiID.'_'.$ratioPos;
my %parameters=(QUANTIF_FAMILY=>'RATIO',VIEW=>'log2',NUM_PEP_CODE=>$numPepCode,QUANTIF_LIST=>[$quantif]); # STRICT_SELECT to match modProtID not loose protID
&promsQuantif::fetchQuantificationData($dbh,\%parameters,\%quantifInfo,\%quantifValues);
#my $countProt=0;
my $log2=log(2);
foreach my $protID0 (keys %{$quantifValues{$quantif}}) { # in case modif-quantif
	$quantifValues{$quantif}{$protID0}{RATIO}=log($quantifValues{$quantif}{$protID0}{RATIO})/$log2;
}

if ($formView) {
	print qq
|<SCRIPT type="text/javascript">
var quantificationData = [
|;
	my $jsCount=0;
	my $numProt=scalar keys %{$quantifValues{$quantif}};
	foreach my $protID0 (keys %{$quantifValues{$quantif}}) { # sort {$quantifValues{$quantif}{$a} <=> $quantifValues{$quantif}{$b}}
		my $pvalue = (defined $quantifValues{$quantif}{$protID0}{P_VALUE})? $quantifValues{$quantif}{$protID0}{P_VALUE} : ''; #'undefined';
		print "[$quantifValues{$quantif}{$protID0}{RATIO},$pvalue,$quantifValues{$quantif}{$protID0}{$numPepCode}";
		if ($isModifQuantif) {
			(my $protID=$protID0)=~s/-.+//;
			print ",$protID";
		}
		print "]";
		print ',' if ++$jsCount < $numProt;
		print "\n" unless $jsCount % 100;
	}
	print qq
|
];
</SCRIPT>
|;
}
# Excluding ratios with a p-value greater than threshold
my @ratios;
foreach my $protID0 (keys %{$quantifValues{$quantif}}) {
	if ($quantifValues{$quantif}{$protID0}{$numPepCode} < $minPep || ($maxPvalue < 1 && (!defined $quantifValues{$quantif}{$protID0}{P_VALUE} || $quantifValues{$quantif}{$protID0}{P_VALUE} > $maxPvalue))) {
	#if ($maxPvalue < 1 && (!defined $quantifValues{$quantif}{$protID}{P_VALUE} || $quantifValues{$quantif}{$protID}{P_VALUE} > $maxPvalue)) {
		delete $quantifValues{$quantif}{$protID0};
		next;
	}
	push @ratios,$quantifValues{$quantif}{$protID0}{RATIO};
}

if (param('AJAX')) {
    $dbh->disconnect;
    $tmpDir=param('tmpDir');
    &refreshPlot(\@ratios,\@thrLogRatios);
    exit;
}

if (param('anaName')) {
    $tmpDir=param('tmpDir');
    &processForm();
    exit;
}

#-------------------------------#
# Fetching Quantification Infos #
#-------------------------------#
my $rank = 0;
my $ratioString = "<SELECT name=\"ratio\" id=\"ratio\" onchange=\"updateRatio(this.value);\">\n";
foreach my $ratio (@{$quantifInfo{$quantiID}[1]{RATIOS}}) {
	$rank++;
	my $selected = ($rank == $ratioPos)? 'selected' : '';
	my ($testCond,$refCond)=split(/\//,$ratio);
	my $supRatioTag='';
	if ($testCond=~/%/) {
		$testCond=~s/%.+//;
		$refCond=~s/%.+//;
		$supRatioTag='°';
	}
	#my $supRatioTag=($quantifInfo{$quantiID}[1]{RATIO_TYPE}[0] eq 'SuperRatio' && $rank >= scalar keys %{$quantifInfo{$quantiID}[2]})? '°' : '';
    $ratioString .= "<OPTION value=\"$rank\" $selected>$quantifInfo{$quantiID}[2]{$testCond}{NAME}$supRatioTag/$quantifInfo{$quantiID}[2]{$refCond}{NAME}$supRatioTag</OPTION>\n";
}

$ratioString .= "</SELECT>\n";
##$sthCondition->finish;

#-------------------------#
# Fetching ontology files #
#-------------------------#
my $disableSubmitStrg = '';
my $ontologyString;
my $sthOBO = $dbh->prepare("SELECT ID_ONTOLOGY, OBO_FILE, NAME FROM ONTOLOGY WHERE STATUS>=1");
$sthOBO->execute;
my $refOboList = $sthOBO->fetchall_arrayref;
if (scalar @{$refOboList}) {
    $ontologyString =  "<SELECT name=\"ontology\">\n";
	$ontologyString .= "<OPTION value=\"\">-= Select =-</OPTION>\n" if scalar @{$refOboList} > 1;
    foreach my $oboEntry (@{$refOboList}) {
		my ($oboID,$oboFileName,$name) = @{$oboEntry};
        $oboFileName =~ s/^.+\/([^\/]+)/$1/; # removing FTP full URL
		$ontologyString .= "<OPTION value=$oboID>$name ($oboFileName)</OPTION>\n";
    }
    $ontologyString .= "</SELECT>";
}
else {
    $ontologyString = "&nbsp;<FONT color=\"red\">No ontology file found in database</FONT>\n";
    $disableSubmitStrg = 'disabled';
}

#---------------------------------------------------------#
# Fetching categories for background population selection #
#---------------------------------------------------------#
my $projectID = promsMod::getProjectID($dbh, $quantiID, 'QUANTIFICATION');
my $sthClass = $dbh->prepare("SELECT ID_CLASSIFICATION,NAME FROM CLASSIFICATION WHERE ID_PROJECT=$projectID ORDER BY NAME");
my $sthCategory = $dbh->prepare("SELECT ID_CATEGORY,NAME FROM CATEGORY WHERE ID_CLASSIFICATION=? ORDER BY DISPLAY_POS");
my $categoryStrg="<OPTION value=\"0\">-= Select =-</OPTION>\n";
$sthClass->execute;
while (my ($classID,$className) = $sthClass->fetchrow_array){
    $categoryStrg .= "<OPTGROUP label=\"Theme: $className\">\n";
    $sthCategory->execute($classID);
    while (my ($categoryID,$categoryName) = $sthCategory->fetchrow_array){
	    $categoryStrg .= "<OPTION value=$categoryID>$categoryName</OPTION>\n";
    }
}
$sthClass->finish;
$sthCategory->finish;

#---------------------------#
# Fetching annotation files #
#---------------------------#
my $annotationString;
my $sthGOA = $dbh->prepare("SELECT ID_GOANNOTATION, NAME, DES, ID_SPECIES, ANNOT_FILE FROM GOANNOTATION WHERE STATUS=1");
$sthGOA->execute;
my $refGOAList = $sthGOA->fetchall_arrayref;
if (scalar @{$refGOAList}) {
    $annotationString = "<SELECT name=\"annotation\">\n<OPTION value=\"\">-= Select =-</OPTION>\n";
	my $sthSp = $dbh->prepare("SELECT SCIENTIFIC_NAME, COMMON_NAME FROM SPECIES WHERE ID_SPECIES=?");
    foreach my $goaEntry (@{$refGOAList}){
		my ($id,$name,$description,$id_species,$annotFile) = @{$goaEntry};
		$sthSp->execute($id_species);
		my ($speciesSciName,$speciesComName) = $sthSp->fetchrow_array;
		$annotationString .= "<OPTION value=$id onmouseover=\"popup('$description')\" onmouseout=\"popout()\">$name - $speciesSciName";
		$annotationString .= " ($speciesComName)" if $speciesComName;
		$annotationString .= "</OPTION>\n";
    }
    $annotationString .= "</SELECT>";
	$sthSp->finish;
}
else {
    $annotationString = "&nbsp;<FONT color=\"red\">No annotation file found in database</FONT>\n";
    $disableSubmitStrg = 'disabled';
}


$dbh->disconnect;
print qq
|<SCRIPT type="text/javascript">
document.getElementById('waitDIV').style.display='none';
</SCRIPT>
<FORM name="enrichForm" id="enrichForm" method="post" enctype="multipart/form-data" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="quanti" value="$quantiID"/>
<BR>
<TABLE bgcolor=$darkColor border=0>
<TR><TH bgcolor=#FFF colspan=12>
<DIV id="scrollbarDiv"><BR><BR><FONT class="title3">Updating $proteinStrg ratios distribution</FONT><BR><IMG src="$promsPath{images}/scrollbarGreen.gif"><BR><BR></DIV>
<DIV id="plotDiv" style="display:none"><IMG id="imgPlot" src="$promsPath{images}/space.gif"></DIV>
</TH></TR>
<TR><TH align="right" colspan=2 nowrap>&nbsp;Log2 ratio thresholds:</TH>
|;
foreach my $idx (0..3) {
    print "<TD bgcolor=$lightColor colspan=2 align=\"center\"><INPUT type=\"text\" id=\"logRatio$idx\" name=\"logRatio$idx\" value=\"$thrLogRatios[$idx]\" style=\"text-align:center; width:50px\" onchange=\"updateBins($idx);\"/></TD>";
}
print qq
|<TH rowspan=2 colspan=2>&nbsp;<BUTTON type="button" id="updateBUT" onclick="refreshPlot();" disabled>Update<BR>chart</BUTTON>&nbsp;</TH></TR>
<TR><TH align="right" colspan=2>Normal scale ratios:</TH>
|;
foreach my $idx (0..3) {
	my $ratioStrg=($thrLogRatios[$idx] < 0)? '1/'.(1*(sprintf "%.2f",1/(2**$thrLogRatios[$idx]))) : 1*(sprintf "%.2f",2**$thrLogRatios[$idx]);
    print "<TD bgcolor=$lightColor colspan=2 align=\"center\"><INPUT type=\"text\" id=\"normRatio$idx\" name=\"normRatio$idx\" value=\"$ratioStrg\" style=\"text-align:center; width:50px\" onchange=\"updateLogRatios($idx);\"/></TD>";
}

print qq
|</TR>
<TR><TH align="right" nowrap>&nbsp;${upProteinStrg}s (%):</TH>|;
foreach my $num (0..4) {
    #print "<TD bgcolor=$lightColor colspan=2 align=\"center\"><INPUT type=\"text\" id=\"bin$_\" name=\"bin$_\" value=\"\" size=2/></TD>"; #  onchange=\"updateLogRatios($_);\"
    print "<TH bgcolor=$lightColor colspan=2><SPAN id=\"bin$num\">...</SPAN></TH>";
}
print qq
|<TH nowrap align=left>:<SPAN id="numProtSPAN">...</SPAN>&nbsp;tot.&nbsp;</TH>
</TR>
|;
if ($isModifQuantif) {
	print qq
|<TR><TH align="right" nowrap>&nbsp;Proteins (%):</TH>|;
	foreach my $num (0..4) {
		print "<TH bgcolor=$lightColor colspan=2><SPAN id=\"binTrueProt$num\">...</SPAN></TH>";
	}
	print qq
|<TH nowrap align=left>:<SPAN id="numTrueProtSPAN">...</SPAN>&nbsp;tot.&nbsp;</TH>
</TR>
|;

}
my $pValueType=(!$quantifInfo{$quantiID}[1]{FDR_CONTROL} || $quantifInfo{$quantiID}[1]{FDR_CONTROL}[0]=~/FALSE|NONE/i)? 'p-value' : 'adj. p-value';
print qq
|</TABLE>
<BR>
<TABLE bgcolor=$darkColor border=0>
        <TR>
            <TH valign='top' align='right'>&nbsp;&nbsp;Name : </TH><TD bgcolor=$lightColor><INPUT type='text' name='anaName' size="30"></TD>
        </TR>
        <TR>
            <TH valign='top' align='right'>&nbsp;&nbsp;Description : </TH><TD bgcolor=$lightColor><TEXTAREA rows="2" cols="60" name="description"></TEXTAREA></TD>
        </TR>
        <TR>
            <TH valign='top' align='right'>&nbsp;&nbsp;Quantification : </TH><TD bgcolor=$lightColor><B>$quantiName</B></TD>
        </TR>
        <TR>
            <TH valign='top' align='right'>&nbsp;&nbsp;Test : </TH><TD bgcolor=$lightColor>$ratioString</TD>
        </TR>
		<TR>
			<TH valign='top' align='right'>&nbsp;&nbsp;Protein filtering : </TH><TD bgcolor=$lightColor>
					&bull;At least <INPUT type='text' name="minPep" value="$minPep" size=2 onchange="updateThreshold();"> quantified peptide(s)/$proteinStrg<BR>
					&bull;Ratio $pValueType &le;<INPUT type="text" name="ratioPvalue" id="ratioPvalue" size=3 value="$maxPvalue" onchange="updateThreshold();">&nbsp;<I>(use <B>1</B> for no threshold)</I><br/>
                    &bull;Restrict to List: <SELECT name="restrictList">$categoryStrg</SELECT>
			</TD>
        </TR>
        <TR>
            <TH valign='top' align='right'>&nbsp;&nbsp;Ontology file : </TH>
            <TD bgcolor=$lightColor>
            $ontologyString
            &nbsp;<B>Depth:</B>&nbsp;<SELECT name="depth">
            <OPTION value=0 selected>All terms</OPTION>
|;
foreach (1..20) {
    print "<OPTION value=\"$_\">$_</OPTION>\n";
}
print qq
|            </SELECT>
            </TD>
        </TR>
        <TR>
            <TH valign='top' align='right'>&nbsp;&nbsp;Annotation : </TH>
            <TD  bgcolor=$lightColor>
                    $annotationString
            </TD>
        </TR>
        <TR>
            <TH valign='top' align='right'>&nbsp;&nbsp;Domain : </TH>
            <TD  bgcolor=$lightColor>
            <INPUT type='radio' name='aspect' value='P' checked> Biological Process &nbsp;&nbsp;
            <INPUT type='radio' name='aspect' value='C'> Cellular Component &nbsp;&nbsp;
            <INPUT type='radio' name='aspect' value='F'> Molecular Function &nbsp;&nbsp;
            </TD>
        </TR>
        <TR>
            <TH valign='top' align='right'>&nbsp;&nbsp;Advanced parameters : </TH><TD bgcolor=$lightColor>
		<TABLE border=0 cellspacing=0>
		<TR><TD>&nbsp;Background population:</TD><TD><INPUT type="radio" name="chPopSource" value="qproteome" checked> All <B>quantified</B> proteome (after <B>Protein filtering</B>)</TD></TR>
		<!--<TR><TD></TD><TD><INPUT type="radio" name="chPopSource" value="aaproteome"> All <B>annotated</B> proteome</TD></TR>-->
		<TR><TD></TD><TD><INPUT type="radio" name="chPopSource" value="category"> Select a List: <SELECT name="population">$categoryStrg</SELECT></TD></TR>
		<TR><TD></TD><TD onmouseover="popup('Background population identifiers must be the same<br>than the one used in the Annotation file.');" onmouseout="popout();"><INPUT type="radio" name="chPopSource" value="localFile"> Upload a local file: <INPUT type="file" name="localPopFile"></TD></TR>
		</TABLE>
                &nbsp;Enrichment test statistical settings:<BR>

                &nbsp;&nbsp;&nbsp;<INPUT type='radio' name='criteria' value='FDR' onClick="selectStatistics();" checked> Control FDR at <INPUT type='text' name='FDR' value='1' size=2>%
                with <SELECT name='FDRmethod'>
                            <OPTION value='BH' selected>Benjamini & Hochberg</OPTION>
                            <OPTION value='FDRsimulation'>Simulation</OPTION>
                </SELECT> method<BR>
                &nbsp;&nbsp;&nbsp;<INPUT type="radio" name="criteria" value="pval" onclick="selectStatistics();"> Use a p-value threshold: <INPUT type="text" name="pval" size=3 value="0.01" disabled>
                with <INPUT type="checkbox" name="bonferroni" value="1" checked disabled> Bonferroni correction<BR>
        </TD></TR>
        <TR>
            <TD align="center" colspan=2>
                <INPUT type="submit" name="start" value="Start" $disableSubmitStrg/> &nbsp;&nbsp;&nbsp;<INPUT type="reset" value="Clear" />
		<INPUT TYPE="hidden" name="tmpDir" value="$tmpDir" />
	    </TD>
        </TR>
    </TABLE>
</FORM>
</CENTER>
|;

##>Generate distribution image
my $image = &densityPlot(\@ratios,\@thrLogRatios);

print qq
|<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">
setPopup();
document.getElementById('scrollbarDiv').style.display='none';
document.getElementById('imgPlot').src='$image';
document.getElementById('plotDiv').style.display='';
</SCRIPT>
</BODY>
</HTML>
|;
exit;

sub refreshPlot{

    my $imageURL = &densityPlot(@_);

    print header(-'content-encoding' => 'none',-charset=>'utf-8');
    warningsToBrowser(1);
    print $imageURL;
}

sub densityPlot{
    my @data = sort {$a <=> $b} @{$_[0]};
    my @thrLogRatios = @{$_[1]};

    my $dataFile = "$promsPath{tmp}/GO/$tmpDir/rData.txt";

    my $size = scalar @data;
    my $percent = 0;
    my $abLineStrg;
    foreach my $thrValue (@thrLogRatios){
        $abLineStrg .= "abline(v='$thrValue', col='red')\n";
        #$abLineStrg .= 'text(x='.log2($data[$i]).", labels=".log2($data[$i]).")\n";
    }

    open Rdata, ">$dataFile" or die $!;
    print Rdata join("\n", @data);
    close Rdata;

    #my $userID = $ENV{'REMOTE_USER'};
    #my $fileName = $userID . '_densityPlot.png';
    #my $imagePath = "$promsPath{tmp}/$fileName";
    my $imagePath= "$promsPath{tmp}/GO/$tmpDir/densityPlot.png";
    unlink $imagePath if -e $imagePath;


    my $rScript = "$promsPath{tmp}/GO/$tmpDir/rScript.R";
    open Rscript, ">$rScript" or die $!;
    print Rscript qq
|data <- t(read.table("$dataFile"))
png(filename="$imagePath", width=500, height=300, units = "px")
hist(sapply(X=data, FUN=function(x) ifelse(x<5, ifelse(x>-5,x,-5.5), 5.5)), xlab="Log2(fold change)", ylab="$upProteinStrg count", main="$upProteinStrg distribution in bins", col='grey', xaxt='n', breaks=20, xlim=c(-6,6))
axis(side=1, at=c(-5,-4,-3,-2,-1,0,1,2,3,4,5), labels=c("<-5",-4,-3,-2,-1,0,1,2,3,4,"5<"))
$abLineStrg
dev.off()
|;
    close Rscript;

    system "cd $promsPath{tmp}/GO/$tmpDir; $promsPath{R}/R CMD BATCH --no-save --no-restore $rScript";

    unlink $dataFile;
    unlink $rScript;
    unlink $rScript . 'out';

    return "$promsPath{tmp_html}/GO/$tmpDir/densityPlot.png";
}

#sub log2{
#    my $x = shift;
#    return (log $x)/(log 2);
#}

sub processForm {
    #---------------------#
    # Fetching parameters #
    #---------------------#
    my $oboID = param('ontology');
    my $depth = (param('depth'))? param('depth') : 0;
    my $goaID = param('annotation');
    my $name = param('anaName');
    my $description = (param('description'))? param('description') : '';
    my $aspect = param('aspect');
    my $criteria = param('criteria');
    my $threshold = param($criteria) or die "unknown field $criteria";
    my $thres = ($criteria eq 'FDR')? param('FDR') : param('pval');
    my $method = ($criteria eq 'FDR')? param('FDRmethod'):(param('bonferroni'))? 'bonferroni':'none';
    my $bckgPop = (param('chPopSource') eq 'category') ? "popCategory=".param('population') : (param('chPopSource') eq 'localFile') ? "popFile=".param('localPopFile') : (param('chPopSource') eq 'aaproteome') ? "pop=allannotatedgenes" : "pop=allquantifiedproteome" ;

    #my $dbh = promsConfig::dbConnect;
    #my $projectID = promsMod::getProjectID($dbh, $anaID, 'ANALYSIS');
    my $projectID = promsMod::getProjectID($dbh, $quantiID, 'QUANTIFICATION');

    my $sthQMeth = $dbh->prepare("SELECT QM.NAME FROM QUANTIFICATION_METHOD QM, QUANTIFICATION Q
                                 WHERE QM.ID_QUANTIFICATION_METHOD=Q.ID_QUANTIFICATION_METHOD AND
                                 ID_QUANTIFICATION=?");
    $sthQMeth->execute($quantiID);
    my ($quantiName) = $sthQMeth->fetchrow_array;
    $sthQMeth->finish;

    print header( -'content-encoding' => 'none',-charset=>'utf-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<TITLE>On-going Quantification GO Analysis...</TITLE>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER><BR><FONT class="title">On-going Quantification GO Analysis...</FONT>
<BR><BR><FONT class="title3">Processing data. Please wait...</FONT><BR>
<IMG src="$promsPath{images}/scrollbarGreen.gif"><BR><BR><BR><BR>
</CENTER>
<FONT class=\"title3\">
|;

    #------------------#
    # Parsing OBO file #
    #------------------#
    print "Reading ontology...";
    my $oboFile = "$promsPath{'obo'}/$oboID.obo";
    my ($status) = $dbh->selectrow_array("SELECT STATUS FROM ONTOLOGY WHERE ID_ONTOLOGY=$oboID");
    my $slim = ($status == 1)? 1:0;
    my %ontologies;
    $ontologies{$aspect} = new promsOboParser( ontologyFile => $oboFile, aspect=> $aspect);
    if ($depth > 0) {
        $ontologies{$aspect} = $ontologies{$aspect}->sliceAtLevel($depth);
    }
    print " Done.<BR>\n";

    #------------------#
    # Parsing GOA file #
    #------------------#
    print "Reading annotation...";
    my $goaFile = "$promsPath{'goa'}/$goaID.goa";
    my $annotation = new GO::AnnotationProvider::AnnotationParser(annotationFile => $goaFile);
    goAnalysis::correctAnnotationForSlim($dbh, $annotation, \%ontologies, [$aspect]);
    print " Done.<BR>\n";

    my $sthGetExp;
    $sthGetExp = $dbh->prepare("SELECT ID_EXPERIMENT FROM DESIGN D,QUANTIFICATION Q WHERE D.ID_DESIGN=Q.ID_DESIGN AND ID_QUANTIFICATION=?");
    $sthGetExp->execute($quantiID);
    my ($expID) = $sthGetExp->fetchrow_array;
    $sthGetExp->finish;
    unless ($expID) {
		my $sthGetExp = $dbh->prepare("SELECT ID_EXPERIMENT FROM SAMPLE S, ANALYSIS A, ANA_QUANTIFICATION AQ
						WHERE S.ID_SAMPLE=A.ID_SAMPLE
						AND AQ.ID_ANALYSIS=A.ID_ANALYSIS
						AND AQ.ID_QUANTIFICATION=?
						LIMIT 0,1");
		$sthGetExp->execute($quantiID);
		($expID) = $sthGetExp->fetchrow_array;
		$sthGetExp->finish;
		die "Unsupported quantification method" unless $expID;
    }

    #-------------------#
    # Creating children #
    #-------------------#
    my $tmpName = strftime("%Y%m%d%H%M%S",localtime);

    my %protBinAlias;
    my %binWarnings;
    my ($identifierID, $identifierName) = $dbh->selectrow_array("SELECT GA.ID_IDENTIFIER,I.NAME FROM GOANNOTATION GA,IDENTIFIER I WHERE ID_GOANNOTATION=$goaID AND GA.ID_IDENTIFIER=I.ID_IDENTIFIER");
    $identifierID = 0 unless $identifierID;
	my $sthAlias = $dbh->prepare("SELECT ALIAS FROM PROTEIN WHERE ID_PROTEIN=?");

    #------------------------------#
    # Fetching foreground proteins #
    #------------------------------#
    my $sthListProt = $dbh->prepare("SELECT DISTINCT P.ID_PROTEIN,P.IDENTIFIER
                                    FROM CATEGORY_PROTEIN CP,PROTEIN P,ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S
                                    WHERE CP.ID_PROTEIN=P.ID_PROTEIN
                                    AND P.ID_PROTEIN=AP.ID_PROTEIN
                                    AND AP.ID_ANALYSIS=A.ID_ANALYSIS
                                    AND A.ID_SAMPLE = S.ID_SAMPLE
                                    AND CP.ID_CATEGORY=? AND S.ID_EXPERIMENT=$expID");
    
    
    print "# Fetching foregroung proteins...<br/>\n";
    my %protIDs; # ${uniprot} = DB_ID
    if(param('restrictList') && param('restrictList') != 0) {
		print "## Fetching protein set from Custom List...<br/>\n";
        my $catID = param('restrictList');
        $sthListProt->execute($catID);
		while(my ($protID, $protName) = $sthListProt->fetchrow_array){
            $protIDs{$protID} = $protName;
		}
    }
    $sthListProt->finish;

	############################################################
	####>Converting sites to proteins (in case PTM quantif)<####
	############################################################
	my %protInBins;
	foreach my $protID0 (keys %{$quantifValues{$quantif}}) {
		(my $protID=$protID0)=~s/-.+//;
        if(!param('restrictList') || (param('restrictList') && defined $protIDs{$protID})) {
            my $currentBin = ($quantifValues{$quantif}{$protID0}{RATIO} <= $thrLogRatios[0])? 0 :
                            ($quantifValues{$quantif}{$protID0}{RATIO} <= $thrLogRatios[1])? 1 :
                            ($quantifValues{$quantif}{$protID0}{RATIO} < $thrLogRatios[2])? 2 :
                            ($quantifValues{$quantif}{$protID0}{RATIO} < $thrLogRatios[3])? 3 : 4;
            $protInBins{$protID}{$currentBin}=1; # only 1 bin for prot quantif BUT can be >1 for site quantif
        }
	}
    my $protCount = scalar keys %protInBins; # %{$quantifValues{$quantif}}
    
	foreach my $protID (keys %protInBins) {
		$sthAlias->execute($protID);
		my ($protAlias)=$sthAlias->fetchrow_array;
		if ($identifierID) {
			my @ids = &goAnalysis::getProteinIds($dbh,$protID,$identifierID);
			my $selectedID = &goAnalysis::selectID(\@ids, $annotation);
			unless ($selectedID){
				$selectedID = $protAlias;
				foreach my $currentBin (sort{$a<=>$b} keys %{$protInBins{$protID}}) {
					push @{$binWarnings{$currentBin}}, "No $identifierName identifier found for $protAlias, this name was used instead";
				}
			}
			foreach my $currentBin (keys %{$protInBins{$protID}}) {
				$protBinAlias{$currentBin}{$selectedID} = $protID;
			}
		}
		else {
			foreach my $currentBin (keys %{$protInBins{$protID}}) {
				$protBinAlias{$currentBin}{$protAlias} = $protID;
			}
		}
	}
	$sthAlias->finish;

    #-----------------------------------#
    # Creating master GO analysis in DB #
    #-----------------------------------#
    my $masterParamStrg = "criteria=$criteria;threshold=$threshold;method=$method;aspect=$aspect;depth=$depth;$bckgPop;minPep=$minPep;logRatios=" . join(',', @thrLogRatios) . ';';
    $masterParamStrg .= "ratioPvalue=$maxPvalue;numProtSel=$protCount;";
	$masterParamStrg .= "numSiteSel=".(scalar keys %{$quantifValues{$quantif}}).';' if $isModifQuantif;
    
    my $userID = $ENV{'REMOTE_USER'};
    my $sthInsGOAna = $dbh->prepare("INSERT INTO GO_ANALYSIS (ID_EXPERIMENT,ID_ONTOLOGY,ID_GOANNOTATION,ID_PARENT_GOANA,NAME,DES,GOA_TYPE,ASPECT,PARAM_STRG,UPDATE_DATE,UPDATE_USER) VALUES (?,?,?,?,?,?,?,?,?,NOW(),?)");
    $sthInsGOAna->execute($expID,$oboID,$goaID,undef,$name,$description,'heatmap',$aspect,$masterParamStrg,$userID);
    my $masterID = $dbh->last_insert_id(undef, undef, 'GO_ANALYSIS', 'ID_GOANALYSIS');

    # and storing density plot image #
    mkdir "$promsPath{go_unix}/project_$projectID" unless -d "$promsPath{go_unix}/project_$projectID";
    mkdir "$promsPath{go_unix}/project_$projectID/$masterID";
    move("$promsPath{tmp}/GO/$tmpDir/densityPlot.png", "$promsPath{go_unix}/project_$projectID/$masterID/densityPlot.png");
    
	#foreach my $protID0 (keys %{$quantifValues{$quantif}}) {
	#	(my $protID=$protID0)=~s/-.+//;
	#	unless ($protAlias{$protID}) {
	#		$sthAlias->execute($protID);
	#		($protAlias{$protID})=$sthAlias->fetchrow_array;
	#	}
	#	my $currentBin = ($quantifValues{$quantif}{$protID0}{RATIO} <= $thrLogRatios[0])? 0 :
	#					($quantifValues{$quantif}{$protID0}{RATIO} <= $thrLogRatios[1])? 1 :
	#					($quantifValues{$quantif}{$protID0}{RATIO} < $thrLogRatios[2])? 2 :
	#					($quantifValues{$quantif}{$protID0}{RATIO} < $thrLogRatios[3])? 3 : 4;
	#
	#	if ($identifierID) {
	#		my @ids = &goAnalysis::getProteinIds($dbh, $protID, $identifierID);
	#		my $selectedID = &goAnalysis::selectID(\@ids, $annotation);
	#		unless ($selectedID){
	#			$selectedID = $protAlias{$protID};
	#			push @{$binWarnings{$currentBin}}, "No $identifierName identifier found for $protAlias{$protID}, this name was used instead";
	#		}
	#		$protBinAlias{$currentBin}{$selectedID} = $protID;
	#	} else {
	#		$protBinAlias{$currentBin}{$protAlias{$protID}} = $protID;
	#	}
	#}
	#$sthAlias->finish;

    #------------#
    # n,n tables #
    #------------#
    #my $sthGOAnaAnalysis = $dbh->prepare("INSERT INTO GOANA_ANALYSIS(ID_GOANALYSIS, ID_ANALYSIS) VALUES(?,?)");
    #$sthGOAnaAnalysis->execute($masterID,$anaID);
    #$sthGOAnaAnalysis->finish;

    my $sthGOAnaQuanti = $dbh->prepare("INSERT INTO GOANA_QUANTIFICATION (ID_GOANALYSIS,ID_QUANTIFICATION,RATIO) VALUES(?,?,?)");
    $sthGOAnaQuanti->execute($masterID,$quantiID,$ratioPos);
    $sthGOAnaQuanti->finish;

    #------------#
    # Population #
    #------------#
    my @population=(); # if left empty, then all annotated genes is defined as population by GO::TermFinder
    print "# Fetching background proteins...<br/>\n";
    if (param('chPopSource') eq 'category') {
		my $categoryID = param('population');
		my $sthProt = $dbh->prepare("SELECT CATEGORY_PROTEIN.ID_PROTEIN,IDENTIFIER FROM CATEGORY_PROTEIN,PROTEIN WHERE ID_CATEGORY=$categoryID
						AND CATEGORY_PROTEIN.ID_PROTEIN=PROTEIN.ID_PROTEIN");
		$sthProt->execute;
		while (my ($protID,$protName) = $sthProt->fetchrow_array){
			if ($identifierID) {
				my @identifiers = goAnalysis::getProteinIds($dbh,$protID,$identifierID);
				if (scalar @identifiers) {
					my $selectedID = goAnalysis::selectID(\@identifiers,$annotation);
					push @population,$selectedID;
				}
				else {
					push @population,$protName;
				}
			}
			else {
				push @population,$protName;
			}
		}
    }
    elsif (param('chPopSource') eq 'localFile') { # All annotated proteome
		my $fileName = tmpFileName(upload('localPopFile'));
        open (POPFILE,$fileName) or die $!;
        while (<POPFILE>) {
            chomp;
            push @population,$_;
        }
        close(POPFILE);
        unlink $fileName;
    }
	elsif (param('chPopSource') eq 'qproteome') {
		foreach my $bin (0..4) { # All quantified proteome (default option)
			push @population, keys %{$protBinAlias{$bin}};
		}
    }
    print "# Done.<BR>\n";


    #---------------------------#
    # Starting enrichment tests #
    #---------------------------#
    my %binData;
    my %significantTerms;
    my @unannotatedProteins;
    my @childToRemove;
    foreach my $bin (0..4){
        unless(scalar keys %{$protBinAlias{$bin}}){
            &processFailed('No protein in bin ' . ($bin+1));
        }

        my $fileName = "$tmpName"."_bin$bin";
		my $threshold = ($criteria eq 'FDR')? 100:1;
        my $size = sprintf "%.1f", ((scalar keys %{$protBinAlias{$bin}}) / $protCount) * 100;

        print 'Performing enrichment analysis for bin '. ($bin+1) . '...';
        &goAnalysis::termFinder(
            dbh => $dbh,
            ontologies => \%ontologies,
            annotation => $annotation,
            aspects => [$aspect],
            data => $protBinAlias{$bin},
            minPep => $minPep,
            criteria => $criteria,
            threshold => $threshold,
            method => $method,
            drawGraph => 0,
            fileName => $fileName,
            slim => 0,
            population => \@population,
            verbose => 0
        );
        my $childParamStrg = "bin=$bin;size=$size;";

        $sthInsGOAna->execute($expID, $oboID, $goaID, $masterID, undef, undef, undef, undef, $childParamStrg, $userID);
        my $goAnaID = $dbh->last_insert_id(undef, undef, 'GO_ANALYSIS', 'ID_GOANALYSIS');
		push @childToRemove, $goAnaID;
		$binData{$bin}{'ID'} = $goAnaID;
        if($binWarnings{$bin} && scalar @{$binWarnings{$bin}}){
            open ERR, ">>$promsPath{tmp}/GO/$fileName/errorLog.txt";
            foreach (@{$binWarnings{$bin}}){
                print ERR $_, "\n";
            }
            close ERR;
        }

		###> Add unannotated proteins percentage to global annotation
		my $resultDirUnix="";
		my %pvalList;
		open RESULTS, "$promsPath{tmp}/GO/$fileName/results_$aspect.txt" or push @unannotatedProteins, 0.00;
		my ($countProt,$countPopProt); # $countProt=total protein in bin <-> $countPopProt=total protein in population
		while (<RESULTS>) {
			if(/^### (\d+)\t(\d+)/){ # 1st line in file
				$countProt = $1;
				$countPopProt = $2;
				next;
			}
			my ($goID,$goName,$pval,$numProt,$totNumProt,$protListStrg,$FDR) = split /\t/;
			my @protList = split /;/, $protListStrg;
			$pvalList{$goID}{'numProt'} = $numProt;
			$pvalList{$goID}{'proteins'} = \@protList;
			if ($goID eq 'unannotated') {
				my $percent = sprintf '%.2f', ($numProt*100.0)/$countProt;
				push @unannotatedProteins,$percent;
				next;
				#last;
			}
			$pvalList{$goID}{'pval'}=$pval;
			next if $pval == 2; # means term is unsignificant but parent in graph (showAllNodes = 0) -> used in image map graph construction
			$pvalList{$goID}{'totNumProt'} = $totNumProt;
			#$significantTerms{$goID}++ if ($criteria eq 'pval' and $pval <= $threshold) or ($criteria eq 'FDR' and $FDR <= $thres/100);#($threshold/100));
			$significantTerms{$goID}++ if ($criteria eq 'pval' and $pval <= $thres) or ($criteria eq 'FDR' and $FDR <= $thres/100);#($threshold/100));
		}
		close RESULTS;
		$binData{$bin}{'pvalList'} = \%pvalList;
        &goAnalysis::storeData($projectID, $goAnaID, $fileName);
        print " Done.<BR>\n";
    }

    $sthInsGOAna->finish;
    my $unannotatedString=join(',', @unannotatedProteins);
    $dbh->do("UPDATE GO_ANALYSIS SET PARAM_STRG=CONCAT(PARAM_STRG,'unannotated=$unannotatedString;') WHERE ID_GOANALYSIS=$masterID");

    # To use hclust and compute a cluster, there must be n >= 2 objects to cluster (ie scalar significantTerms must be at least equal to 2)
    # otherwise, an R error will be launched and the script will fail to produce a heatmap
    unless(scalar keys %significantTerms > 1) {
		print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<CENTER>
<FONT class="title">Less than two significant terms in this Q.GO enrichment.<BR>No heatmap can be drawn with such a small set.</FONT>
</CENTER>
</BODY>
</HTML>
|;
		$dbh->rollback;
		$dbh->disconnect;
		rmtree("$promsPath{go_unix}/project_$projectID/$masterID");
		foreach my $chilID (@childToRemove) {
			rmtree("$promsPath{go_unix}/project_$projectID/$chilID");
		}
		exit;
    }

    $dbh->commit;
    $dbh->disconnect;

    #my $tmpDir = strftime("%Y%m%d%H%M%S",localtime)."_$masterID";
    #mkdir "$promsPath{tmp}/$tmpDir";

    # Merging all results #
    my $mergeFile = "$promsPath{tmp}/GO/$tmpDir/merge.txt";
    open MERGE, ">$mergeFile";
    my $i = 0;
	my $log10=log(10);
    foreach my $goid (keys %significantTerms){
		$i++;
		print MERGE $goid;
		foreach my $bin (sort {$a <=> $b} keys %binData){
			my $pval = (defined $binData{$bin}{'pvalList'}{$goid}{'pval'})? $binData{$bin}{'pvalList'}{$goid}{'pval'} : 1; # R does not accept NA values for clustering
			print MERGE "\t", $pval;
			$binData{$bin}{'pvalList'}{$goid}{'pval'} = ((defined $binData{$bin}{'pvalList'}{$goid}{'pval'}))? -1*log($binData{$bin}{'pvalList'}{$goid}{'pval'})/$log10: 0;
		}
		print MERGE "\n";
    }
    close MERGE;

    #&writeRscriptCluster($tmpDir,$mergeFile);
    system "cd $promsPath{tmp}/GO/$tmpDir; $promsPath{R}/R CMD BATCH --no-save --no-restore '--args $mergeFile' $promsPath{R_scripts}/goQuantiCluster.R";
    foreach(glob("$promsPath{tmp}/GO/$tmpDir/*")){
        #copy($_,"$promsPath{go_unix}/project_$projectID/$masterID/") || print STDERR "$_: $!";
        move($_,"$promsPath{go_unix}/project_$projectID/$masterID/") || print STDERR "$_: $!";
    }
    rmtree("$promsPath{tmp}/GO/$tmpDir");

    print qq
|Processing finished. Displaying result...</FONT>
|;
    sleep(3);
    print qq
|<SCRIPT type="text/javascript">
	parent.itemFrame.location='$promsPath{cgi}/openProject.cgi?ACT=experiment&VIEW=go&branchID=go_analysis:$masterID';
</SCRIPT>
</BODY>
</HTML>
|;
}

sub processFailed{
    my $errorText = shift;

    $dbh->rollback();
    $dbh->disconnect();

    print qq
|<FONT color="red">Processing failed : $errorText</FONT><BR>
<BR>
<CENTER>
<INPUT type="button" value="Go back" onclick="window.history.back();">
</CENTER>
</BODY>
</HTML>
|;
exit;
}

#sub writeRscriptCluster {
#   my ($tmpDir, $mergeFile)=@_;
#   open (R, ">$promsPath{tmp}/$tmpDir/goQuantiCluster.R");
#   print R "
#   library(\"gplots\")
#   data <- read.table(\"$mergeFile\",header=F,row.names=1,sep=\"\\t\")
#   data.log <- log10(as.matrix(data))
#   colnames(data.log) <- c(\"bin1\",\"bin2\",\"bin3\",\"bin4\",\"bin5\")
#
#   # Z-scoring
#   for(i in 1:nrow(data.log)){
#      mean.l <- mean(data.log[i,])
#      sd.l <- sd(data.log[i,])
#      data.log[i,] <- (data.log[i,]-mean.l)/sd.l
#   }
#
#   clust <- hclust(dist(data.log), method=\"average\")
#   write.table(data.frame(clust\$merge,sort(clust\$height)),file=\"extractDendro.txt\",quote=FALSE, sep=\"\\t\", col.names=NA )
#   hmCarpet <- heatmap.2(as.matrix(data.log), distfun = function(x) dist(x,method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'average'), key=T, symkey=FALSE, density.info=\"none\", trace=\"none\")\$carpet
#   writeLines(colnames(hmCarpet), con=\"extractOrder.txt\")
#   ";
#   close R;
#}

#sub log10{
#    my $x = shift;
#    return (log $x)/(log 10);
#}

####>Revision history
# 1.3.2 [ENHANCEMENT] Uses data from &promsConfig::fetchQuantificationData to display ratio info (PP 25/02/20)
# 1.3.1 Add foreground filtering option (VS 08/01/20)
# 1.3.0 Compatible with PTM quantification (PP 27/07/18)
# 1.2.0 Move tempDir for density plot inside $promsPath{tmp}/GO & added old tempDir clean up with promsMod::cleanDirectory (PP 06/02/18)
# 1.1.0 add cluster function from displayGOQuantiAnalysis (SL 02/12/14)
# 1.0.9 Update for Super-Ratio/Design selection (GA 17/07/14)
# 1.0.8 Add scrollbar waiting div while reloading of histogram (GA 19/05/14)
# 1.0.7 Add background population option and reloading of histogram<BR>TODO: check list efficiency in selection and add wait frame while reloading of histogram (GA 14/05/14)
# 1.0.6 Add checkMinPep() function that had been forgotten (GA 06/05/14)
# 1.0.5 Same modification than selectAnalyses.cgi 2.5.0 in annot pattern matching (GA 30/04/14)
# 1.0.4 Change syntax $x=$y // $z to old ()? one (PP 10/03/14)
# 1.0.3 Help popup near main title (FY 08/10/13)
# 1.0.2 Managing extreme ratios in R histogram (FY 11/09/13)
# 1.0.1 Reformatting param fetching to focus GO analysis on quantification and ratio items instead of analysis (FY 28/03/13)
# 1.0.0 New script to set GO Quanti analysis parameters (FY 07/12/12)
