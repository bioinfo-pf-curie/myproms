#!/usr/local/bin/perl -w

################################################################################
# displayPCA.cgi       1.2.21                                                  #
# Authors: P. Poullet, S.Liva (Institut Curie)      	                       #
# Contact: myproms@curie.fr                                                    #
# display and store the results of PCA analysis        	                       #
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
use promsQuantif;
use POSIX qw(strftime); # to get the time
use Text::CSV;
use Spreadsheet::WriteExcel;
use utf8; # Tells Perl that characters are UTF8. Necessary for Excel export to work with UTF-8 characters ...
#use Encode qw(encode_utf8); # ... Encode needed if hard UTF8 chars in code!!!

#######################
####>Configuration<####
#######################
#print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG (WARNING affects highlighting match)
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my ($explorID,$experimentID) = &promsMod::cleanNumericalParameters(param('explorID'),param('experimentID'));
my $ajax = param('AJAX') || '';
my $pcaType = param('selectPCA') || "quantif";
my $action = param('ACT') || '';
if ($action eq 'full3D') {
	&displayFull3D;
	exit;
}

####>Connexion to DB<####
my $dbh = &promsConfig::dbConnect;
my $projectID=param('PROJECT_ID') || &promsMod::getProjectID($dbh,$explorID,'EXPLORANA');
my $pathToFile = "$promsPath{explorAna}/project_$projectID/$explorID";
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
my $disabSave=($projectAccess eq 'guest')? ' disabled' : '';


#my $selQuantifModif=0;
my $quantifListStrg='';
my %quantificationIDs;
my $isPtmQuantif=0;
my %ptmUsed;

$dbh->do("SET SESSION group_concat_max_len = 1000000");
my $sthSelExplorQuantif=$dbh->prepare("SELECT Q.ID_QUANTIFICATION,GROUP_CONCAT(TARGET_POS ORDER BY TARGET_POS SEPARATOR '_'),Q.ID_MODIFICATION,GROUP_CONCAT(DISTINCT MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ',')
								FROM QUANTIFICATION Q
								LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
								INNER JOIN EXPLORANA_QUANTIF EQ ON EQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION
								WHERE ID_EXPLORANALYSIS=? GROUP BY Q.ID_QUANTIFICATION");
$sthSelExplorQuantif->execute($explorID);
while (my ($quantifID,$targetPosStrg,$quantifModID,$multiModifStrg) = $sthSelExplorQuantif->fetchrow_array) {
    $quantifListStrg.=':' if $quantifListStrg;
    $quantifListStrg.=$quantifID."_".$targetPosStrg;
    $quantificationIDs{$quantifID}=1;
	($isPtmQuantif,my $isMultiModifQuantif)=&promsQuantif::getQuantifModificationInfo($dbh,$quantifModID,$multiModifStrg,{},\%ptmUsed);
}
$sthSelExplorQuantif->finish;

my $featureItems='Proteins';
if (scalar keys %ptmUsed) {
    $isPtmQuantif=1;
	$featureItems='Sites';
}


####>EXPORT call<####
if ($action eq 'export') {
    &exportProteins();
    exit;
}

####>AJAX calls<####
if ($ajax eq 'changeDimensions') {
    &changeDimensions(param('dimX'),param('dimY'),param('dimZ'),'false');
    exit;
}
elsif ($ajax eq 'listSigProtDim') {
    &ajaxListSignificantProteins;
    exit;
}
elsif ($ajax=~/^property:/) {
    &ajaxGetPropertyValues;
    exit;
}
elsif ($ajax eq 'propDecorateGraph') {
    &ajaxPropDecorateGraph;
    exit;
}
elsif ($ajax eq 'customLists') {
    &ajaxGetCustomLists;
    exit;
}
elsif ($ajax eq 'listDecorateGraph') {
    &ajaxListDecorateGraph;
    exit;
}
elsif ($ajax eq 'saveAnnotations') {
    &ajaxSaveAnnotations;
    exit;
}

my %listThemes;
my $sthTh=$dbh->prepare("SELECT ID_CLASSIFICATION,NAME FROM CLASSIFICATION WHERE ID_PROJECT=$projectID");
$sthTh->execute;
while (my($themeID,$themeName)=$sthTh->fetchrow_array) {
    $listThemes{$themeID}=$themeName;
}
$sthTh->finish;
my (%goAnalyses, %parentGoAna, %pathwayAnalyses);

##GO ANALYSES
my $sthGO = $dbh -> prepare("SELECT ID_GOANALYSIS,NAME,ASPECT,ID_PARENT_GOANA FROM GO_ANALYSIS WHERE ID_EXPERIMENT=$experimentID");
$sthGO -> execute;
while (my ($goID,$name,$aspectStrg,$parentGoID) = $sthGO -> fetchrow_array) {
    if ($parentGoID) {
        push @{$parentGoAna{$parentGoID}},$goID;
    }
    else {
        @{$goAnalyses{$goID}} = ($name,$aspectStrg);
    }
}
$sthGO -> finish;

##PATHWAY ANALYSES
my $sthPA=$dbh->prepare("SELECT ID_PATHWAY_ANALYSIS, ID_CATEGORY, NAME, PARAM_STRG FROM PATHWAY_ANALYSIS WHERE ID_EXPERIMENT=$experimentID AND STATUS>0 AND ANALYSIS_TYPE='PATHWAY'");
$sthPA->execute;
while (my ($pathID, $catID, $name, $paramStrg)=$sthPA->fetchrow_array) {
    $catID = (!$catID)? "" : $catID;
    @{$pathwayAnalyses{$pathID}}=($catID, $name, $paramStrg);
}
$sthPA->finish;

my ($pcaName,$paramList,$filterList) = $dbh -> selectrow_array("SELECT NAME,PARAM_LIST,FILTER_LIST FROM EXPLORANALYSIS where ID_EXPLORANALYSIS = $explorID");
my ($quantifFam,$quantifMeasCode,$pepType)=('','',''); # prot PCA only
my %bioSampProperties;
if ($pcaType =~ /^quantif/) {
    my $sthProp=$dbh->prepare("SELECT DISTINCT P.ID_PROPERTY,P.NAME,P.PROPERTY_TYPE FROM PROPERTY P
                                INNER JOIN BIOSAMPLE_PROPERTY BP ON P.ID_PROPERTY=BP.ID_PROPERTY
                                INNER JOIN OBSERVATION O ON BP.ID_BIOSAMPLE=O.ID_BIOSAMPLE
                                INNER JOIN ANA_QUANTIFICATION AQ ON O.ID_ANALYSIS=AQ.ID_ANALYSIS
                                WHERE P.USE_IN_ANALYSIS=1 AND AQ.ID_QUANTIFICATION IN (".join(',',keys %quantificationIDs).")");
    $sthProp->execute;
    while (my($propID,$propName,$propType)=$sthProp->fetchrow_array) {
        $bioSampProperties{$propType}{$propID}=$propName;
    }
    $sthProp->finish;
}
else { # prot
	($quantifFam)=$paramList=~/quantFam=(\w+)/;
	$quantifMeasCode=$1 if $paramList=~/quantCode=(\w+)/;
	$pepType=$1 if ($filterList && $filterList=~/\/PEP=(\w+)/);
	$pepType=($pepType eq 'NUM_PEP_USED')? 'all' : ($pepType eq 'DIST_PEP_USED')? 'distinct' : $pepType;
}

#### START HTML
print header(-'content-encoding'=>'no',-charset=>'utf-8');
warningsToBrowser(1);
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
print qq
|<HTML>
<HEAD>
<TITLE>PCA</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.TD {font-weight:normal;}
.TH {font-weight:bold;}
.row_0{background-color:$darkColor;}
.row_1{background-color:$lightColor;}
.highlight{width:250px;}
.popup {z-index:999;background-color:#FFFFFF;border:solid 3px #999999;padding:5px;position:absolute;display:none;box-shadow:10px 10px 20px #808080;}
</STYLE>
<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/pcaPlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
<SCRIPT type="text/javascript">
|;
&promsMod::popupInfo();
print qq
|var full3dWindow;
function displayFull3D() {
	myForm=document.full3dForm;
	//Highlighting
	var objString = '';
    var concat = 0;
    for (let i in annotSetObject) {
        //var nameSubstr=i.replace(/\\[.+\|\\].+/,"");
        concat++;
        if (concat > 1) objString += ':%:'; // object separator
        objString += i+'&&'+annotSetObject[i].position+'&&'+annotSetObject[i].color+'&&'+annotSetObject[i].ajaxResp;
    }
	myForm.highlight.value=objString;
	//Display window
	full3dWindow=window.open('','full3dWindow','width=950,height=950,location=no,resizable=yes,scrollbars=yes');
	myForm.submit();
	full3dWindow.focus();
}
function editDeleteHighligth(name, action, newName) {
    if (action=='delete') {
        var pos = annotSetObject[name].position;
        for (let i in annotSetObject ) {
            if (annotSetObject[i].position > pos){
                annotSetObject[i].position --;
            }
        }
        rank--;
    }
    else {
        annotSetObject[newName]={};
        for (let j in annotSetObject[name]) {
            annotSetObject[newName][j]=annotSetObject[name][j];
        }
    }
    delete annotSetObject[name];
    document.getElementById('saveButton').style.display='';
}

var XHR = null;
function getXMLHTTP() {
    var xhr = null;
    if (window.XMLHttpRequest) {// Firefox & others
        xhr = new XMLHttpRequest();
    }
    else if (window.ActiveXObject) { // Internet Explorer
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

//AJAX TO CHANGE DIMENSIONS
function ajaxChangeDimensions() {
    var selDimX=document.getElementById('selectDimX');
    var selDimY=document.getElementById('selectDimY');
    var selDimZ=document.getElementById('selectDimZ');
    var dimX=selDimX.value;
    var dimY=selDimY.value;
	var dimZ=selDimZ.value;
    if (!dimX \|\| !dimY) {
        alert('Select valid values for X/Y axis dimension!');
        return;
    }
    var dimXlabel=selDimX.options[selDimX.selectedIndex].text;
    var dimYlabel=selDimY.options[selDimY.selectedIndex].text;

    //If XHR object already exists, the request is canceled & the object is deleted
    if (XHR && XHR.readyState != 0) {
        XHR.abort();
        delete XHR;
    }
    //Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if ( !XHR ) {
        return false;
    }
    XHR.open("GET","./displayPCA.cgi?AJAX=changeDimensions&dimX="+dimX+"&dimY="+dimY+"&dimZ="+dimZ+"&selectPCA="+selPCAtype+"&experimentID=$experimentID&explorID=$explorID",true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState == 4 && XHR.responseText) {
            PCA.resetData(XHR.responseText);
            PCA.redraw(dimXlabel,dimYlabel);
			//full3D form
			if (dimZ) {
				var myForm=document.full3dForm;
				myForm.axisTitles.value=dimXlabel+';'+dimYlabel+';'+selDimZ.options[selDimZ.selectedIndex].text;
				myForm.pointData.value=XHR.responseText;
				document.getElementById('full3dBUT').disabled=false;
			}
			else { // disable full3D button
				document.getElementById('full3dBUT').disabled=true;
			}
        }
    }
    XHR.send(null);
}

//DISPLAY PCA DATA
var selPCAtype='$pcaType';
function displayPCAdata(pcaType) {
	if (selPCAtype.replace('_sc','')==pcaType.replace('_sc','')) { // same PCA target -> ajax update
		selPCAtype=pcaType;
		ajaxChangeDimensions();
	}
	else {
	    window.location="./displayPCA.cgi?experimentID=$experimentID&explorID=$explorID&selectPCA="+pcaType;
	}
}
var annotSetObject = {};
|;
if ($pcaType =~/^quantif/) { # PCA ON QUANTIFICATIONS
    print qq
|
//AJAX FOR QUANTIFICATIONS
var sigProtSort='p-value';
function selectSigProtSort(newSort) {
    sigProtSort=newSort;
    ajaxListSignificantProteins();
}
function ajaxListSignificantProteins() {
    saveListGlobals.themesFetched=false; // resets ajaxManageSaveProteins mecanism
    var listDiv=document.getElementById('protListDIV');
    listDiv.innerHTML="<BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<FONT class=\\"title3\\">Fetching data. Please wait...</FONT><BR>&nbsp;&nbsp;&nbsp;&nbsp;<IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR>";
    listDiv.style.display='';
    listDiv.scrollIntoView({block:"start",inline:"nearest",behavior:"smooth"});

    var dim=document.getElementById('selectSigProtDim').value;
    //If XHR object already exists, the request is canceled & the object is deleted
    if (XHR && XHR.readyState != 0) {
        XHR.abort();
        delete XHR;
    }
    //Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if ( !XHR ) {
        return false;
    }
    XHR.open("GET","./displayPCA.cgi?AJAX=listSigProtDim&dim="+dim+"&sort="+sigProtSort+"&selectPCA="+selPCAtype+"&experimentID=$experimentID&explorID=$explorID",true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState == 4 && XHR.responseText) {
	    listDiv.innerHTML=XHR.responseText;
        }
    }
    XHR.send(null);
}
function exportProteins() {

    var selBox=document.getElementsByName('chkProt');
    var numProt=0;
    for (var i=0; i<selBox.length;i++) {
        if (selBox[i].checked) {
            numProt++;
            break;
        }
    }
    if (!numProt) {
        alert('ERROR: select at least 1 protein');
        return;
    }
    var defAction=document.protForm.ACT.value;
    document.protForm.ACT.value='export';
    document.protForm.submit();
    document.protForm.ACT.value=defAction;
}
function extendSelection(index,chk) {
    var auto=document.getElementById('autoExtend');
    if (!auto.checked){return;}
    var selBox=document.getElementsByName('chkProt');
    for (var i=index; i<selBox.length;i++) {
        selBox[i].checked=chk;
    }
}

function ajaxGetPropertyValues(selProp) {
    var selectDiv=document.getElementById('quantifHlDIV');
    if (!selProp) {
        selectDiv.innerHTML='';
        selectDiv.style.display = 'none';
        return;
    }

    //If XHR object already exists, the request is canceled & the object is deleted
    if (XHR && XHR.readyState != 0) {
        XHR.abort();
        delete XHR;
    }
    //Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if (!XHR) {
        return false;
    }

    XHR.open("GET","./displayPCA.cgi?AJAX="+selProp+"&experimentID=$experimentID&explorID=$explorID&selectPCA="+selPCAtype,true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState == 4 && XHR.responseText) {
            //var multiStrg=(selProp.match(':O:')? ' multiple' : '';
            selectDiv.innerHTML=XHR.responseText;
            selectDiv.style.display = '';
        }
    }
    XHR.send(null);
}

function prepareAjaxPropDecorateGraph(propValueData) {
    if (!propValueData) return;
    var valueSelect=document.getElementById('quantifHlSEL');
    var propSelect=document.getElementById('quantifHighlightType');
    ajaxPropDecorateGraph([ propValueData,valueSelect.options[valueSelect.selectedIndex].text,propSelect.value, propSelect.options[propSelect.selectedIndex].text ],false);
}
function ajaxPropDecorateGraph(params,saved,jobRank) {
    var [propValueData,propValueText,selPropInfo,propName]=params;
	if (!saved) saved=false;
    var paramStrg='AJAX=propDecorateGraph&explorID=$explorID&experimentID=$experimentID&qList=$quantifListStrg&propValueData='+encodeURIComponent(propValueData);

    //If XHR object already exists, the request is canceled & the object is deleted
    if (XHR && XHR.readyState != 0) {
        XHR.abort();
        delete XHR;
    }
	//Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if (!XHR) {
        return false;
    }
    XHR.open("POST","$promsPath{cgi}/displayPCA.cgi",true); //!saved Switches to synchronous for already saved highlights
    //Send the proper header information along with the request
    XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
    XHR.onreadystatechange=function() {
        if (XHR.readyState==4) {
            if (XHR.responseText) {
                colorIndex++;
                if (colorIndex >= colorList.length) colorIndex=0;
                var selData=selPropInfo.split(':'); // property:O:propID
                var propInfo=propValueData.split(':=:'); // propID:selected value
                var hlName=(selData[1]=='O')? propValueText+' ['+propName+']' : propName+': '+propValueText; // propInfo[1]+' ['+propName+']';
				if (addHighlighting(PCA,hlName,colorList[colorIndex],{'-1':XHR.responseText.split(';')},'^###(\$\|%)')) {
                    rank++;
                    annotSetObject[hlName] = {
						color:colorList[colorIndex],
						type:selData[0]+':'+selData[1],
						param:'#'+propValueData,
						position:rank,
						ajaxResp:XHR.responseText
					}; // # for ID tag
                    if (colorIndex==colorList.length) colorIndex=0;
                    if (!saved) document.getElementById('saveButton').style.display = '';
                }
                else {
                    colorIndex--;
                    if (colorIndex < 0) colorIndex=colorList.length-1;
                }
            }
            else {alert('No match found!');}
			
			/* Launch next job (if any) */
			launchNextJob(jobRank);
        }
    }
    XHR.send(paramStrg);
}

|;
}
else { # PCA ON PROTEINS
	print "var goAspects=new Object();\n";
	foreach my $goID (keys %goAnalyses) {
	    print "goAspects[$goID]='$goAnalyses{$goID}[1]';\n";
	}
	print qq
|function updateProtHighlightType(typeID) {
    if (!typeID) {
        document.getElementById('goAspect').style.display='none';
        document.getElementById('selTermsDIV').style.display='none';
        return;
    }
    var typeInfo=typeID.split(':');
    if (typeInfo[0]=='GOA') {
        updateGoAspects(typeInfo[1]);
    }
    else if (typeInfo[0]=='TH') {
        document.getElementById('goAspect').style.display='none';
        ajaxGetCustomLists(typeInfo[1]);
    }
    else if (typeInfo[0]=='PA') {
        document.getElementById('goAspect').style.display='none';
        document.getElementById('selTermsDIV').style.display='none';
        ajaxGetPathwayList(typeInfo[1]);
    }
}
//AJAX FOR GO
function updateGoAspects(goID) {
    var selAspect=document.getElementById('goAspect');
    selAspect.options.length=1;
    document.getElementById('selTermsDIV').innerHTML='';
    if (goAspects[goID].match('P')) {selAspect.options[1]=new Option('Biological process',goID+',P');}
    if (goAspects[goID].match('C')) {selAspect.options[2]=new Option('Cellular component',goID+',C');}
    if (goAspects[goID].match('F')) {selAspect.options[3]=new Option('Molecular function',goID+',F');}
    selAspect.style.display='';
}

function ajaxUpdateGoTermList(goIdStrg) {
    var termsDiv=document.getElementById('selTermsDIV');
    termsDiv.innerHTML='';
    if (!goIdStrg) {
        return;
    }

    //If XHR object already exists, the request is canceled & the object is deleted
    if (XHR && XHR.readyState != 0) {
        XHR.abort();
        delete XHR;
    }
    //Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if (!XHR) {
        return false;
    }
    XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxGoTerms&projectID=$projectID&goStrg="+goIdStrg,true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState==4 && XHR.responseText) {
            termsDiv.innerHTML=XHR.responseText;
            termsDiv.style.display='';
        }
    }
    XHR.send(null);
}

function ajaxGoDecorateGraph(termValue,termTxt) { // called ajax call to showProtQuantification.cgi
	ajaxMyGoDecorateGraph([termValue,termTxt]); // parameter conversion required for compatibility with highlightPCA()
}
function ajaxMyGoDecorateGraph(params,saved,jobRank) {
	var [termValue,termTxt]=params;
    var binArray=[];
    binArray=termValue.split(',');
    var nbElem=binArray.length;

    if (!saved) saved=false;
    if (!termValue) return;

    //If XHR object already exists, the request is canceled & the object is deleted
    if (XHR && XHR.readyState != 0) {
        XHR.abort();
        delete XHR;
    }
    //Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if (!XHR) {
        return false;
    }

    termTxt = termTxt.replace(/ \\[GO.+/,""); // not when saved is true

    if (nbElem==4) {
        var binInfo=binArray.pop();
        var binValues = binInfo.split(':');
        var binStrg = (binValues[0] == 0 \|\| binValues[0] == 1)?' ]'+binValues[1]+']' : (binValues[0] == 3 \|\| binValues[0] == 4)?' ['+binValues[1]+'[' :' ]'+binValues[1]+'[';
        termTxt=termTxt.concat(binStrg);
    }

    //var listParam = termValue.split(",");
    //var strgParam = listParam[0]+','+listParam[1];

    XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxTermProt&projectID=$projectID&goStrg="+termValue,true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState==4 && XHR.responseText) {
            colorIndex++;
			if (colorIndex >= colorList.length) colorIndex=0;
            var termData=termValue.split(',');
            if (addHighlighting(PCA,termTxt,colorList[colorIndex],{'-1':XHR.responseText.split(';')},'^###(\$\|-)')) {
                //colorList[colorIndex] = colorList[colorIndex].replace(/#/,"");
                rank++;
                annotSetObject[termTxt] = {color:colorList[colorIndex],type:'prot:GO',param:'#'+termValue,position:rank,ajaxResp:XHR.responseText}; // # for ID tag
                if (colorIndex==colorList.length) colorIndex=0;
                if (!saved) document.getElementById('saveButton').style.display = '';
            }
            else {colorIndex--;}
			
			/* Launch next job (if any) */
			launchNextJob(jobRank);
        }
    }
    XHR.send(null);
}
//AJAX FOR LIST HIGHLIGHTING
function ajaxGetCustomLists(themeID) {
    var termsDiv=document.getElementById('selTermsDIV');

    //If XHR object already exists, the request is canceled & the object is deleted
    if (XHR && XHR.readyState != 0) {
        XHR.abort();
        delete XHR;
    }
    //Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if (!XHR) {
        return false;
    }
    XHR.open("GET","$promsPath{cgi}/displayPCA.cgi?AJAX=customLists&experimentID=$experimentID&explorID=$explorID&themeID="+themeID,true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState==4 && XHR.responseText) {
            termsDiv.innerHTML=XHR.responseText;
            termsDiv.style.display='';
        }
    }
    XHR.send(null);
}
function ajaxListDecorateGraph(params,saved,jobRank) {
	var [listID]=params;
    if (!saved) saved=false;
    if (!listID) return;

    //If XHR object already exists, the request is canceled & the object is deleted
    if (XHR && XHR.readyState != 0) {
        XHR.abort();
        delete XHR;
    }
    //Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if (!XHR) {
        return false;
    }
    XHR.open("GET","$promsPath{cgi}/displayPCA.cgi?AJAX=listDecorateGraph&experimentID=$experimentID&explorID=$explorID&listID="+listID,true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState==4) {
            if (XHR.responseText) {
                colorIndex++;
                if (colorIndex >= colorList.length) colorIndex=0;
                var listData=XHR.responseText.split('::');
                if (addHighlighting(PCA,listData[0],colorList[colorIndex],{'-1':listData[1].split(';')},'^###(\$\|-)')) {
                    rank++;
                    annotSetObject[listData[0]] = {color:colorList[colorIndex],type:'prot:LIST',param:'#'+listID,position:rank,ajaxResp:XHR.responseText}; // # for ID tag
                    if (colorIndex==colorList.length) colorIndex=0;
                    if (!saved) document.getElementById('saveButton').style.display = '';
                }
                else {colorIndex--;}
            }
            else {alert('List is empty!');}
			
			/* Launch next job (if any) */
			launchNextJob(jobRank);
        }
    }
    XHR.send(null);
}

//AJAX FOR PATHWAY
function ajaxGetPathwayList(pathID) {
    var termsDiv=document.getElementById('selTermsDIV');
    //If XHR object already exists, the request is canceled & the object is deleted
    if (XHR && XHR.readyState != 0) {
        XHR.abort();
        delete XHR;
    }
    //Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if (!XHR) {
        return false;
    }
    XHR.open("GET","$promsPath{cgi}/displayPathwayAnalysis.cgi?AJAX=ajaxGetPathway&FROM=PCA&ID="+pathID,true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState==4 && XHR.responseText) {
            termsDiv.innerHTML=XHR.responseText;
            termsDiv.style.display='';
        }
    }
    XHR.send(null);
}
function ajaxGetPathwayProteinsList(params,saved,jobRank) {
	var [value,pathID,txt]=params;
    if (!saved) saved=false;
    //If XHR object already exists, the request is canceled & the object is deleted
    if (XHR && XHR.readyState != 0) {
        XHR.abort();
        delete XHR;
    }
    //Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if (!XHR) {
        return false;
    }
    XHR.open("GET","$promsPath{cgi}/displayPathwayAnalysis.cgi?AJAX=ajaxListProt&ID="+pathID+"&FROM=PCA&reactNumber="+value,true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState==4 && XHR.responseText) {
            colorIndex++;
            if (colorIndex >= colorList.length) colorIndex=0;
            if (addHighlighting(PCA,txt,colorList[colorIndex],{'-1':XHR.responseText.split(';')},'^###(\$\|-)')) {
                rank++;
                annotSetObject[txt] = {color:colorList[colorIndex],type:'prot:PA',param:'#'+pathID+','+value,position:rank,ajaxResp:XHR.responseText}; // # for ID tag
                if (colorIndex==colorList.length) colorIndex=0;
                if (!saved) document.getElementById('saveButton').style.display = '';
            }
            else {colorIndex--;}
			
			/* Launch next job (if any) */
			launchNextJob(jobRank);
        }
    }
    XHR.send(null);
}

//OTHER AJAX FOR PROTEINS
function ajaxListFromPCA(selectedPoints) {
    ajaxListSelectedProteins(selectedPoints[0].join(','),'protein');
}
function ajaxSearchConvertIdentifier(graphSearch,graphSearchArgs,searchTextIdx) { // (graph lib search function,array of function arguments,index of search text in array). To be called at end of convertion function
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxConvIdent&projectID=$projectID&quantifList=$quantifListStrg&TEXT="+encodeURIComponent(graphSearchArgs[searchTextIdx]),true); //search text
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			graphSearchArgs[searchTextIdx]=XHR.responseText; // replace old text with converted one
			graphSearch(...graphSearchArgs); // call graph lib search fuction & convert array to arguments
		}
	};
	XHR.send(null);
}
function ajaxListSelectedProteins(selectedPoints,sort='protein') {
    saveListGlobals.themesFetched=false; // resets ajaxManageSaveProteins mecanism
    var listDiv=document.getElementById('protListDIV');
    listDiv.innerHTML="<BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<FONT class=\\"title3\\">Fetching data. Please wait...</FONT><BR>&nbsp;&nbsp;&nbsp;&nbsp;<IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR>";
    listDiv.style.display='';
    listDiv.scrollIntoView({block:"start",inline:"nearest",behavior:"smooth"});
    var paramStrg="AJAX=ajaxListProt&CALL=PCA&id_project=$projectID&ACT=results&quantifFamily=$quantifFam&dispMeasure=$quantifMeasCode&pepType=$pepType&quantifList=$quantifListStrg&sort="+sort+"&selProt="+selectedPoints;

    //If XHR object already exists, the request is canceled & the object is deleted
    if (XHR && XHR.readyState != 0) {
        XHR.abort();
        delete XHR;
    }
    //Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if (!XHR) {
        return false;
    }
    XHR.open("POST","$promsPath{cgi}/compareQuantifications.cgi",true);
    //Send the proper header information along with the request
    XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
    XHR.onreadystatechange=function() {
        if (XHR.readyState==4 && XHR.responseText) {
            listDiv.innerHTML=XHR.responseText;
        }
    }
    XHR.send(paramStrg);
}
|;
}

print qq
|function selectSort(newSort,ajax) { // view=list (true or from Volcano ajax)
    var checkedProt=document.protForm.chkProt;
    var chkList=[];
    if (checkedProt.length) {
        for (var i=0; i<checkedProt.length; i++) {chkList.push(checkedProt[i].value);}
    }
    else {chkList.push(checkedProt.value);}
    ajaxListSelectedProteins(chkList.join(','),newSort);
}
function sequenceView(id_protein,id_analyses) {
    var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+id_analyses+"&id_prot="+id_protein+"&msdata="+top.showMSdata;
    top.openProtWindow(winLocation);
}

function checkAllProteins(chkStatus) {
    var checkBoxList=document.protForm.chkProt;
    if (checkBoxList.length) {
        for (var i=0; i < checkBoxList.length; i++) {checkBoxList[i].checked=chkStatus;}
    }
    else {checkBoxList.checked=chkStatus;}
}
|;
&promsMod::printAjaxManageSaveProteins($projectID,\%promsPath,'document.protForm.chkProt'); # no need for callback
print qq
|function ajaxSaveAnnotations() {
    //Object to string
    var objString = '';
    var concat = 0;
    for (let i in annotSetObject) {
        var nameSubstr=i.replace(/\\[.+\|\\].+/,"");
        concat++;
        if (concat > 1) objString += ':%:'; // object separator
        objString += 'name=='+nameSubstr; // == because annotSetObject[i][j] contains "=" if Treatment!
        for (let j in annotSetObject[i]) {
			if (j=='color' \|\| j=='ajaxResp') continue;
            objString += '//'+j+'=='+annotSetObject[i][j]; // == because annotSetObject[i][j] contains "=" if Treatment!
        }
    }
    var paramStrg="AJAX=saveAnnotations&experimentID=$experimentID&explorID=$explorID&selectPCA="+selPCAtype+"&string="+encodeURIComponent(objString);

    //If XHR object already exists, the request is canceled & the object is deleted
    if (XHR && XHR.readyState != 0) {
        XHR.abort();
        delete XHR;
    }
    //Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if (!XHR) {
        return false;
    }
    XHR.open("POST","$promsPath{cgi}/displayPCA.cgi",true);
    //Send the proper header information along with the request
    XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
    XHR.onreadystatechange=function() {
        if (XHR.readyState==4) {
            if (XHR.responseText.match('###OK###')) {
                document.getElementById('saveButton').style.display='none';
                var saveSpan=document.getElementById('saveSPAN');
                saveSpan.style.display='';
                setTimeout(function(){saveSpan.style.display='none';},5000);
            }
            else {
                alert('ERROR: Highlights could not be saved !');
            }
        }
    }
    XHR.send(paramStrg);
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">PCA : <FONT color=#DD0000>$pcaName</FONT></FONT><BR><BR>
<DIV id="waitDIV">
<BR><BR><FONT class="title3">Fetching data. Please wait...</FONT><BR>&nbsp;&nbsp;&nbsp;&nbsp;<IMG src="$promsPath{images}/scrollbarGreen.gif"><BR><BR>
</DIV>
</CENTER>
<DIV id="resultDIV" style="visibility:hidden"> <!-- display:none is not compatible with SVG text drawing with Raphaël in Chrome! -->
<INPUT type="button" id="full3dBUT" class="title3" value=" Full 3D view " onclick="displayFull3D()" disabled/>
|;

#### PCA ###
my %contribFiles=(prot=>'protContrib.txt',prot_sc=>'protContrib_sc.txt',quantif=>'quantifContrib.txt',quantif_sc=>'quantifContrib_sc.txt');
my %dimContribution;
open (CONTRIB, "$pathToFile/$contribFiles{$pcaType}" ) || &stopOnError("ERROR!<BR>$!: '$pathToFile/$contribFiles{$pcaType}'");
while (<CONTRIB>) {
    next if $. == 1;
    chomp;
    my ($dim, $contribValue) = split("\t",$_);
    #last if ($. > 4 && $contribValue < 5);
    $dimContribution{$dim} = sprintf"%.2f",$contribValue;
	last if $. == 11; # 10th dim max
}
close CONTRIB;
my $numDimensions=scalar keys %dimContribution;
my $stringPCA=&changeDimensions(1,2,3,'true'); # 1st dim, 2nd dim, 3rd dim, true:, transpo or not
my $pcaTarget=($pcaType=~/^quantif/)? 'Quantifications' : $featureItems; #'Proteins';
print qq
|<SCRIPT type="text/javascript">
var colorList = ['#0000FF','#4AA02C','#F660AB','#FBB917','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B'];
var colorIndex = -1;
var PCA;
var rank = 0;
function drawPCA() {
    PCA = new pcaPlot({div:'pcaDIV',width:500,height:500,
                       name:'PCA',
                       axisX:'Dim 1 ($dimContribution{1} %)',
                       axisY:'Dim 2 ($dimContribution{2} %)',
                       //pointOnclick:externalLink,
|;
print qq
|						pointOnList:ajaxListFromPCA,
						searchable:{text:'Extended search',externalSearch:ajaxSearchConvertIdentifier},
| if $pcaType=~/^prot/;
print qq
|                      allowHighlight:true,
                       updateHighlight:{callback:editDeleteHighligth},
					   exportAsImage:['Export as image','PCA_$pcaTarget\_$pcaName','./exportSVG.cgi']
					   });
    PCA.addDataAsString('$stringPCA');
    PCA.draw();
}
function highlightPCA(jobRank=1) {
	if (jobList[jobRank]) {
		let [hlFunction,params]=jobList[jobRank];
		hlFunction(params,true,jobRank);
	}
}
function launchNextJob(jobRank) {
	if (jobRank) {
		jobRank++; // next job
		if (jobList[jobRank]) {highlightPCA(jobRank);}
	}
}
var jobList={};
|;

###>Restoring saved highlights<###
my $sthSelectAnnot = $dbh -> prepare("SELECT NAME,ANNOT_RANK,ANNOT_TYPE,ANNOT_LIST FROM ANNOTATIONSET WHERE ID_EXPLORANALYSIS = $explorID ORDER BY ANNOT_RANK");
my $sthProperty=$dbh->prepare("SELECT NAME FROM PROPERTY WHERE ID_PROPERTY=?");
$sthSelectAnnot -> execute;
my $jobRank=0;
while (my ($name, $rank, $annotType, $annotList) = $sthSelectAnnot -> fetchrow_array) {
	$jobRank++;
    #my ($color, $selValue) = split("=", $annotList);
    if ($pcaType=~/^quantif/) {
        #print "displayParamLabel('$annotType','$annotList','$name',true);\n";
        if ($annotType =~ /^property:(.)/) {
            my $propType=$1;
            $annotList=~s/^#//; # remove propID tag
            my ($propID,$propValue)=split(':=:',$annotList);
            my $propValueText=($propType eq 'T')? &promsQuantif::convertTreatmentToText($propValue) : $propValue;
            $sthProperty -> execute($propID);
            my ($propName)=$sthProperty -> fetchrow_array;
            next unless $propName; # to be safe in case property has been deleted
			#print "\trestoreSavedHighlightings($jobRank,ajaxPropDecorateGraph,['$annotList','$propValueText','$annotType:$propID','$propName']);\n";
			print "jobList[$jobRank]=[ajaxPropDecorateGraph,['$annotList','$propValueText','$annotType:$propID','$propName']];\n";
        }
    }
    elsif ($pcaType =~ /^prot/ && $annotType=~ /^prot:/) {
        if ($annotType eq 'prot:GO') {
            my $termText=$annotList;
            $termText=~ s/\s+$//;
            $name=~s/^#//; # remove GOanaID tag
            my $termValue=$name;
			#print "\trestoreSavedHighlightings($jobRank,ajaxMyGoDecorateGraph,['$termValue','$termText']);\n";
			print "jobList[$jobRank]=[ajaxMyGoDecorateGraph,['$termValue','$termText']];\n";
         }
        elsif ($annotType eq 'prot:LIST') {
            my ($listID)=($name=~/^#(\d+)/); # extract cat ID
			#print "\trestoreSavedHighlightings($jobRank,ajaxListDecorateGraph,[$listID]);\n";
			print "jobList[$jobRank]=[ajaxListDecorateGraph,[$listID]];\n";
        }
        elsif ($annotType eq 'prot:PA') {
            $name=~s/^#//;
            my ($pathID, $reactNumber)=split(",", $name);
			#print "\trestoreSavedHighlightings($jobRank,ajaxGetPathwayProteinsList,['$reactNumber',$pathID,'$annotList']);\n";
			print "jobList[$jobRank]=[ajaxGetPathwayProteinsList,['$reactNumber',$pathID,'$annotList']];\n";
        }
    }
}
$sthSelectAnnot -> finish;
$sthProperty -> finish;

$dbh -> disconnect;
my ($selQuantif,$selQuantifSc,$selProt,$selProtSc)=($pcaType eq 'quantif')? ('selected','','','') : ($pcaType eq 'quantif_sc')? ('','selected','','') : ($pcaType eq 'prot')? ('','','selected','') : ('','','','selected');
print qq
|</SCRIPT>
<TABLE cellspacing=0>
    <TR class="row_0">
        <TH class="title3" nowrap>&nbsp;PCA type:
            <SELECT class="title3" name="selectPCA" onchange="displayPCAdata(this.value)">
			<!-- Protein/site PCAs are no longer computed
                <OPTION value="quantif" $selQuantif>Quantifications</OPTION>
                <OPTION value="quantif_sc" $selQuantifSc>Quantifications [scaled PCA]</OPTION>
                <OPTION value="prot" $selProt>$featureItems</OPTION>
                <OPTION value="prot_sc" $selProtSc>$featureItems [scaled PCA]</OPTION>
			-->
				<OPTION value="quantif" $selQuantif>Unscaled</OPTION>
                <OPTION value="quantif_sc" $selQuantifSc>Scaled</OPTION>
            </SELECT>
        </TH>
        <TH nowrap>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Show:</TH>
        <TH nowrap>
|;

foreach my $dimAxis ('X','Y','Z') {
    print qq
|<SELECT id="selectDim$dimAxis">
    <OPTION value="">-= $dimAxis =-</OPTION>
|;
    foreach my $dim (sort{$a <=> $b} keys %dimContribution) {
        print "<OPTION value=\"$dim\"";
        print ' selected' if (($dimAxis eq 'X' && $dim==1) || ($dimAxis eq 'Y' && $dim==2) || ($dimAxis eq 'Z' && $dim==3));
        print ">Dim $dim ($dimContribution{$dim} %)</OPTION>\n";
    }
    print "</SELECT>\n";
    print "vs." if $dimAxis eq 'X';
	print "size:" if $dimAxis eq 'Y';
}

print qq
|	</TH>
    <TD nowrap><INPUT type="button" onclick="ajaxChangeDimensions()" value="Change dimensions"></TD>
|;

if ($pcaType=~/^quantif/) {
    print "<TH nowrap>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$featureItems significant for:<SELECT id=\"selectSigProtDim\">";
    foreach my $dim (sort{$a <=> $b} keys %dimContribution) {print "<OPTION value=\"$dim\">Dim $dim</OPTION>\n";}
    print "</SELECT><INPUT type=\"button\" value=\"List\" onclick=\"ajaxListSignificantProteins()\"></TH>\n";
}

###>Data for full 3D view
my ($pointLabelStrg,$pointDataStrg)=('','');
if ($numDimensions >= 3) {
	foreach my $pointStrg (split(';',$stringPCA)) {
		my ($alias,@data)=split(',',$pointStrg);
		if ($pointLabelStrg) {
			$pointLabelStrg.=';';
			$pointDataStrg.=';';
		}
		$pointLabelStrg.=$alias;
		$pointDataStrg.=join(',',@data);
	}
}

print qq
|    </TR>
</TABLE>

<TABLE>
    <TR>
        <TD>
            <DIV id="pcaDIV" style="float:left;background-color:#FFFFFF;border:solid 3px $darkColor;padding:2px;"></DIV>
        </TD>
        <TD valign="top">
            <TABLE cellspacing=0>
                <TR bgcolor="$darkColor"><TD class="title3">&nbsp;Highlights PCA:&nbsp;</TD></TR>
|;
if ($pcaType =~/^quantif/) { # PCA on quantifications
		print qq
|               <TR>
                    <TD>
                        <SELECT name="quantifHighlightType" id="quantifHighlightType" class="highlight title3" onchange="ajaxGetPropertyValues(this.value)">
                            <OPTION value="">-=Select=-</OPTION>
							<!--
                                <OPTGROUP LABEL="Quantifications">
                                    <OPTION value="samp:quantif:name" disabled>Name</OPTION>
                                </OPTGROUP>
						    -->
                                <OPTGROUP LABEL="BioSample properties:">
|;
		if (scalar keys %{$bioSampProperties{'O'}}) {
                    foreach my $propID (sort{lc($bioSampProperties{'O'}{$a}) cmp lc($bioSampProperties{'O'}{$b})} keys %{$bioSampProperties{'O'}}) {
                        print "<OPTION value=\"property:O:$propID\">$bioSampProperties{O}{$propID}</OPTION>\n";
                    }
		}
		else {print "<OPTION value=\"\" disabled>None</OPTION>\n";}
		print qq|						</OPTGROUP>
								<OPTGROUP label="BioSample treatments:">
|;
		if (scalar keys %{$bioSampProperties{'T'}}) {
                    foreach my $propID (sort{lc($bioSampProperties{'T'}{$a}) cmp lc($bioSampProperties{'T'}{$b})} keys %{$bioSampProperties{'T'}}) {
                        print "<OPTION value=\"property:T:$propID\">$bioSampProperties{T}{$propID}</OPTION>\n";
                    }
		}
		else {print "<OPTION value=\"\" disabled>None</OPTION>\n";}
		print qq|						</OPTGROUP>
                        </SELECT>
                    </TD>
                </TR>
                <TR>
                    <TD>
                        <DIV id="quantifHlDIV">
						<!--
                            <B>Highlight name:</B><INPUT type="text" name="labelName" size="20" maxlength="50" value=""/><br>
                            <DIV id="addLabel" style="display:none"></DIV>
						-->
                        </DIV>
                    </TD>
                </TR>
|;
}
else { # PCA on proteins
	print qq
|                     <TR>
                        <TD>
                                <SELECT class="highlight title3" onchange="updateProtHighlightType(this.value)">
								<OPTION value="">-= Select =-</OPTION>
								<OPTGROUP label="GO Analyses:">
|;
	if (scalar keys %goAnalyses) {
		foreach my $goID (sort{lc($goAnalyses{$a}[0]) cmp lc($goAnalyses{$b}[0])} keys %goAnalyses) {
			print "<OPTION value=\"GOA:$goID\">$goAnalyses{$goID}[0]</OPTION>\n";
		}
	}
	else {print "<OPTION value=\"\" disabled>** No GO annotations **</OPTION>\n";}
	print qq
|								</OPTGROUP>
								<OPTGROUP label="Custom Lists:">
|;
		if (scalar keys %listThemes) {
			foreach my $themeID (sort{lc($listThemes{$a}) cmp lc($listThemes{$b})} keys %listThemes) {
				print "<OPTION value=\"TH:$themeID\">$listThemes{$themeID}</OPTION>\n";
			}
		}
		else {print "<OPTION value=\"\" disabled>** No List found **</OPTION>\n";}
		print qq
|								</OPTGROUP>
                                                                <OPTGROUP label="Pathway Analyses:">
|;
                if (scalar keys %pathwayAnalyses) {
                    foreach my $pathID (sort{lc($pathwayAnalyses{$a}[1]) cmp lc($pathwayAnalyses{$b}[1])} keys %pathwayAnalyses) {
                        print "<OPTION value=\"PA:$pathID\">$pathwayAnalyses{$pathID}[1]</OPTION>\n";
                    }
                }
                else {print "<OPTION value=\"\" disabled>** No Pathway annotations **</OPTION>\n";}
		print qq
|                                                                </OPTGROUP>
				    </SELECT>
                            </TD>
                        </TR>
                        <TR>
                            <TD>
                                <SELECT id="goAspect" class="highlight title3" style="display:none" onchange="ajaxUpdateGoTermList(this.value)">
                                    <OPTION value="">-= Select Aspect =-</OPTION>
                                </SELECT>
                            </TD>
                        </TR>
                        <TR>
                            <TD><DIV id="selTermsDIV"></DIV></TD>
                        </TR>
|;
}
print qq
|        				<TR><TH>
							<INPUT type="button" id="saveButton" style="font-weight:bold;display:none" onclick="ajaxSaveAnnotations()" value="Save highlights"$disabSave>
							<SPAN id="saveSPAN" class="title3" style="color:#FFF;background-color:#387D38;padding:1px 7px;display:none">Highlights saved !</SPAN>
						</TH></TR>
					</TABLE>
	</TD>
    </TR>
</TABLE>
</DIV>
<DIV id="protListDIV" style="position:absolute;height:500px;padding-right:25px;overflow:auto"></DIV>
<BR><BR>
<DIV id="divDescription" class="clDescriptionCont"></DIV>
<SCRIPT type="text/javascript">
setPopup();

/********* DRAWING PCA *********/
drawPCA();

document.getElementById('waitDIV').style.display='none';
document.getElementById('resultDIV').style.visibility='visible';

/********* HIGHLIGHTS PCA *********/
highlightPCA();

/********* FULL 3D PCA *********/
if ($numDimensions >= 3) {document.getElementById('full3dBUT').disabled=false;}

</SCRIPT>
|;
my @axisTitlesArray;
foreach my $i (1..$numDimensions){
	push @axisTitlesArray, "Dim $i ($dimContribution{$i} %)";
}
my $axisTitleStrg=join(";",@axisTitlesArray);

print qq|
<FORM name="full3dForm" method="POST" target="full3dWindow">
<INPUT type="hidden" name="ACT" value="full3D"/>
<INPUT type="hidden" name="explorID" value="$explorID"/>
<INPUT type="hidden" name="pcaName" value="$pcaName"/>
<INPUT type="hidden" name="axisTitles" value="$axisTitleStrg"/>
<INPUT type="hidden" name="pointLabels" value="$pointLabelStrg"/>
<INPUT type="hidden" name="pointData" value="$pointDataStrg"/>
<INPUT type="hidden" name="highlight" value=""/>

</FORM>

</BODY>
</HTML>
|;

#$dbh -> disconnect;
sub changeDimensions {
    my ($dimX, $dimY, $dimZ, $return) = @_;

    if ($return eq 'false') {
        print header(-'content-encoding' => 'no',-charset => 'utf-8');warningsToBrowser(1);
    }

    #my (%protList,$sthProt,%quantifNames);
    my %quantifNames;
    if ($pcaType=~/^quantif/) { # PCA on samples
        %quantifNames=&promsQuantif::getDistinctQuantifNames($dbh,$quantifListStrg);
    }
#    else { # protein view
#		$sthProt=$dbh->prepare("SELECT ALIAS FROM PROTEIN WHERE ID_PROTEIN=?");
#    }
	
	my $oldSingleModID=(keys %{$ptmUsed{NAME}})[0]; # in case of old single PTM-site format
	my %protIdList;
    my %coordFiles=(prot=>'protCoordinates.txt',prot_sc=>'protCoordinates_sc.txt',quantif=>'quantifCoordinates.txt',quantif_sc=>'quantifCoordinates_sc.txt');
    my @coord;
	my $coordIdx=0;
    open(COORD, "$pathToFile/$coordFiles{$pcaType}") || die "Error: $!";
    while (my $lineCoord = <COORD>) {
        chomp($lineCoord);
        my $line;
        next if ($. == 1);
        my (@infoLine) = split("\t", $lineCoord);
        my $alias;
        if ($pcaType=~/^prot/) {
			my ($protID,$modStrg)=($infoLine[0]=~/^(\d+)-?(.*)/);
			#unless ($protList{$protID}) {
			#	$sthProt->execute($protID);
			#	($protList{$protID})=$sthProt->fetchrow_array;
			#}
			if ($modStrg && $modStrg !~ /^-?\d+:/) { # no starting modifID => old single PTM format // potential Free res (-1:...)
				$modStrg=$oldSingleModID.':'.$modStrg;
				$infoLine[0]=$protID.'-'.$modStrg;
			}
            #$alias=($modStrg)? $protList{$protID}.'-'.$modStrg : $protList{$protID};
			if ($return eq 'true') {
				$alias=$protID; # temp. changed later
				push @{$protIdList{$protID}},[$coordIdx,$modStrg]; # undef for protein quantif
			}
        }
        else {
            if ($infoLine[0]=~/%/) {
                my ($test, $ref)=split(/%/, $infoLine[0]);
                #my $aliasTest=(split(/:/,$quantifNames{'OPTIMAL'}{$test}))[1];
                (my $aliasTest=$quantifNames{'OPTIMAL'}{$test})=~s/.*://;

                #my $aliasRef=(split(/:/,$quantifNames{'OPTIMAL'}{$ref}))[1];
                (my $aliasRef=$quantifNames{'OPTIMAL'}{$ref})=~s/.*://;
                $alias=$aliasTest." % ".$aliasRef;
                #$alias=$quantifNames{'OPTIMAL'}{$test}."_".$quantifNames{'OPTIMAL'}{$ref};
            }
            else {
                $alias=$quantifNames{'OPTIMAL'}{$infoLine[0]};
            }
        }
        if ($return eq 'true') {
            $line = join(",",$alias,$infoLine[0],$infoLine[$dimX],$infoLine[$dimY]) ;
        }
        else {
            $line = join(",",$infoLine[0],$infoLine[$dimX],$infoLine[$dimY]);
        }
		$line.=','.$infoLine[$dimZ] if ($dimZ && $infoLine[$dimZ]);
        push @coord,$line;
		$coordIdx++;
    }
    close(COORD);
	#$sthProt->finish if $pcaType=~/^prot/;

	if ($return eq 'true' && $pcaType=~/^prot/) {
		
		# my %quantifModifInfo;
		# if ($isPtmQuantif) {
		# 	my @quantifModifs=keys %ptmUsed;
		# 	&promsQuantif::getQuantifModificationInfo($dbh,\@quantifModifs,\%quantifModifInfo);
		# }
		
		my $protIdStrg=join(',',keys %protIdList);
		my $sthProt=$dbh->prepare("SELECT ID_PROTEIN,ALIAS FROM PROTEIN WHERE ID_PROTEIN IN ($protIdStrg)");
		$sthProt->execute;
		while (my ($protID,$alias)=$sthProt->fetchrow_array) {
			foreach my $refModCode (@{$protIdList{$protID}}) {
				my ($coordIdx,$modCode)=@{$refModCode};
				my $proteinName;
				if ($modCode) {
					#my ($formatCode,$displayCode)=&promsQuantif::formatProteinModificationSites($modCode,\%quantifModifInfo,'text');
					my $displayCode=&promsQuantif::displayModificationSites($modCode,$ptmUsed{DISPLAY},'text');
                    $proteinName=$alias.'-'.$displayCode;
				}
				else {$proteinName=$alias;}
				$coord[$coordIdx]=~s/^$protID/$proteinName/;
			}
		}
		$sthProt->finish;
	}

    my $stringScale = join(";",@coord);
	$stringScale=~s/'/./g; $stringScale=~s/"/./g;
    if ($return eq 'false') {print $stringScale;}
    else {return $stringScale;}
}

sub ajaxListSignificantProteins {
    print header(-type=>'text/plain', -'content-encoding' => 'no',-charset => 'utf-8');
    warningsToBrowser(1);

    my $dim=param('dim');
    my $sortOrder=(param('sort')) || 'p-value';
 
    my %protData;
	&fetchProteinsFromDimFile($dim,'html',\%protData);
	
	my $lowerFeatures=lc($featureItems);
	if (scalar keys %protData==0) {
		$dbh->disconnect;
		print qq|<FONT class="title2">No significant $lowerFeatures found for dimension $dim</FONT><BR>|;
		exit;
	}
		
    my $sthA=$dbh->prepare("SELECT DISTINCT ID_ANALYSIS FROM EXPLORANA_QUANTIF EQ,ANA_QUANTIFICATION AQ WHERE EQ.ID_QUANTIFICATION=AQ.ID_QUANTIFICATION AND ID_EXPLORANALYSIS=$explorID");
    my @anaList;
    $sthA->execute;
    while (my($anaID)=$sthA->fetchrow_array) {
        push @anaList,$anaID;
    }
    $sthA->finish;
    my $anaIDstrg=join(',',@anaList);

    $dbh->disconnect;

	my $numTotProteins=scalar keys %protData;
    my $proteinText=($sortOrder eq 'protein')? "<FONT color=\"#DD0000\">$lowerFeatures</FONT>" : $lowerFeatures;
    my $correlText=($sortOrder eq 'correlation')? '<FONT color="#DD0000">Correlation</FONT>' : 'Correlation';
    my $pvalueText=($sortOrder eq 'p-value')? '<FONT color="#DD0000">p-value</FONT>' : 'p-value';
	print qq
|<FORM name="protForm" method="POST">
<INPUT type="hidden" name="PROJECT_ID" value="$projectID"/>
<INPUT type="hidden" name="sort" value="$sortOrder"/>
<INPUT type="hidden" name="ACT" value="" />
<INPUT type="hidden" name="explorID" value="$explorID" />
<INPUT type="hidden" name="experimentID" value="$experimentID" />
<INPUT type="hidden" name="protList" value="" />
<INPUT type="hidden" name="dim" value="$dim" />

<FONT class="title2">List of significant $lowerFeatures for dimension $dim</FONT><BR>
<TABLE border=0 cellspacing=0 cellpadding=2 align=center>
<TR><TD colspan=8><DIV id="saveProtDIV" style="display:none;"></DIV></TD></TR>
<TR><TD colspan=8>
	<INPUT type="button" value="Clear list" style="font-weight:bold" onclick="document.getElementById('protListDIV').style.display='none'">
	<INPUT type="button" value="Check all" onclick="checkAllProteins(true)"/><INPUT type="button" value="Uncheck all" onclick="checkAllProteins(false)"/><INPUT type="button" id="saveFormBUTTON" value="Save proteins..." onclick="ajaxManageSaveProteins('getThemes','PROT')"$disabSave/>|;
	print qq|<INPUT type="button" id="saveSiteFormBUTTON" value="Save sites..." onclick="ajaxManageSaveProteins('getThemes','SITE')"$disabSave/>| if $isPtmQuantif;
	print qq
|<INPUT type="button" value="Export $lowerFeatures" onclick="exportProteins()"/>
</TD></TR>
<TR><INPUT type="checkbox" id="autoExtend" value="1" checked><FONT class="title3">Auto-extend selection</FONT></TR>
<TR class="row_0">
<TH class="rbBorder" align=left nowrap>&nbsp;$numTotProteins <A href="javascript:selectSigProtSort('protein')" onmouseover="popup('Click to sort proteins by <B>ascending name</B>.')" onmouseout="popout()">$proteinText</A>&nbsp;</TH>
<TH class="rbBorder" align="left" nowrap>&nbsp;Gene Name&nbsp;</TH>
<TH class="rbBorder">&nbsp;<A href="javascript:selectSigProtSort('correlation')" onmouseover="popup('Click to sort proteins by <B>decreasing correlation</B>.')" onmouseout="popout()">$correlText</A>&nbsp;</TH>
<TH class="rbBorder">&nbsp;<A href="javascript:selectSigProtSort('p-value')" onmouseover="popup('Click to sort proteins by <B>ascending p-value</B>.')" onmouseout="popout()">$pvalueText</A>&nbsp;</TH>
<TH class="rbBorder">&nbsp;MW<SUP>kd</SUP>&nbsp;</TH><TH class="bBorder">Description - Species</TH>
</TR>
|;
    my $numProt=0;
    my $index=-1;
    foreach my $modProtID (sort{&sortSigProteins($sortOrder,\%protData)} keys %protData) {
        $numProt++;
        $index++;
        my $trClass='row_'.($numProt % 2);
		my ($protID,$modStrg)=split('-',$modProtID);
		my ($protLabel,$correl,$pValue,$mw,$des,$species,$refGeneList)=@{$protData{$modProtID}};
        print qq
|<TR class="list $trClass" valign=top><TD class="TH" nowrap><INPUT type="checkbox" name="chkProt" value="$modProtID" onchange="extendSelection($index,this.checked)" /><A href="javascript:sequenceView($protID,'$anaIDstrg')">$protLabel</A>&nbsp;</TD>
<TH>|;
		my $numGenes=scalar @{$refGeneList};
		if ($numGenes==0) {print '-';}
        elsif ($numGenes==1) {print $refGeneList->[0];}
		else {
            print "<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-",join('<BR>&nbsp;&nbsp;-',@{$refGeneList}[1..$#{$refGeneList}]),"</B>')\" onmouseout=\"popout()\"><U>",$refGeneList->[0],"</U></A>&nbsp;";
        }
        print qq|</TH>
<TH>$correl</TH><TH>$pValue</TH>
<TD class="TH" align=right>$mw&nbsp;</TD><TD>$des <U><I>$species</I></U></TD>
</TR>
|;
    }
    print qq
|<TR><TD colspan=8><B>End of list.</B></TD></TR>
</TABLE>
</FORM>
|;
	exit;
}

sub ajaxGetPropertyValues {
    print header(-type=>'text/plain',-charset => 'utf-8');
    warningsToBrowser(1);

    #my $sthEQ = $dbh -> prepare("SELECT GROUP_CONCAT(DISTINCT ID_QUANTIFICATION SEPARATOR ',') FROM EXPLORANA_QUANTIF WHERE ID_EXPLORANALYSIS=$explorID GROUP BY ID_EXPLORANALYSIS");
    my ($propType,$propID)=($ajax=~/property:(.):(\d+)/);
    my $sthPV=$dbh->prepare("SELECT DISTINCT BP.PROPERTY_VALUE FROM BIOSAMPLE_PROPERTY BP
                                                            INNER JOIN OBSERVATION O ON BP.ID_BIOSAMPLE=O.ID_BIOSAMPLE
                                                            INNER JOIN ANA_QUANTIFICATION AQ ON O.ID_ANALYSIS=AQ.ID_ANALYSIS
                                                            INNER JOIN EXPLORANA_QUANTIF EQ ON AQ.ID_QUANTIFICATION=EQ.ID_QUANTIFICATION
                                                            WHERE BP.ID_PROPERTY=$propID AND EQ.ID_EXPLORANALYSIS=$explorID");
    #my ($propName)=($propType eq 'T')? $dbh->selectrow_array("SELECT NAME FROM PROPERTY WHERE ID_PROPERTY=$propID") : ('');
    $sthPV->execute;
    my %propValueCode;
    while (my ($propValue)=$sthPV->fetchrow_array) {
        my $propText=$propValue;
        if ($propType eq 'T') {
            $propValue=~s/^stepValue=\d:*//; # step independent: stepValue=1::quantity=xxx... -> quantity=xxx...
            $propValue='%' unless $propValue;
            $propText=&promsQuantif::convertTreatmentToText($propValue);
        }
        $propValueCode{$propValue}=$propText;
    }
    $sthPV->finish;
    $dbh->disconnect;

    print "<SELECT id=\"quantifHlSEL\" class=\"highlight title3\" onchange=\"prepareAjaxPropDecorateGraph(this.value)\">\n";
    if (scalar keys %propValueCode) {
        print "<OPTION value=\"\">-= Select =-</OPTION>\n";
        print "<OPTION value=\"$propID:=:%\">treated</OPTION>\n" if $propType eq 'T'; # will be skipped if already in list
        foreach my $propValue (sort{&promsMod::sortSmart(lc($propValueCode{$a}),lc($propValueCode{$b}))} keys %propValueCode) {
            next if $propValue eq '%';
            print "<OPTION value=\"$propID:=:$propValue\">$propValueCode{$propValue}</OPTION>\n";
        }
    }
    else {
        print "<OPTION value=\"\" disabled>No values</OPTION>\n";
    }
    print "</SELECT>\n";
    exit;
}

sub ajaxPropDecorateGraph {
    print header(-type=>'text/plain',-charset => 'utf-8');
    warningsToBrowser(1);

    my ($propID,$propertyValues)=split(':=:',param('propValueData'));
    my ($propType)=$dbh->selectrow_array("SELECT PROPERTY_TYPE FROM PROPERTY WHERE ID_PROPERTY=$propID");

    my %matchedQuantifs=&promsQuantif::getQuantificationsFromPropertyValue($dbh,$propID,$propType,$propertyValues,param('qList'));

    $dbh->disconnect;

    print join(';',keys %matchedQuantifs);
    exit;
}

sub ajaxGetCustomLists {
    print header(-type=>'text/plain',-charset => 'utf-8');
    warningsToBrowser(1);

    my $themeID=param('themeID');
    my $sthList=$dbh->prepare("SELECT ID_CATEGORY,NAME,LIST_TYPE FROM CATEGORY WHERE ID_CLASSIFICATION=$themeID ORDER BY DISPLAY_POS");
    $sthList->execute;
    my $htmlCode='';
	my %listTypeName=('PROT'=>'Proteins','SITE'=>'Sites');
    while (my ($listID,$listName,$listType)=$sthList->fetchrow_array) {
		$listType='PROT' unless $listType;
        $htmlCode.="<OPTION value=\"$listID\">$listName [$listTypeName{$listType}]</OPTION>\n";
    }
    $sthList->finish;
    $dbh->disconnect;

    print "<SELECT id=\"listHlSEL\" class=\"highlight title3\" onchange=\"ajaxListDecorateGraph([this.value])\">\n";
    if ($htmlCode) {
        print "<OPTION value=\"\">-= Select a List =-</OPTION>\n$htmlCode";
    }
    else {
        print "<OPTION value=\"\" disabled>No Lists</OPTION>\n";
    }
    print "</SELECT>\n";
    exit;
}
sub ajaxListDecorateGraph {
    print header(-type=>'text/plain',-charset => 'utf-8');
    warningsToBrowser(1);

    my $listID=param('listID');
    my ($listName)=$dbh->selectrow_array("SELECT NAME FROM CATEGORY WHERE ID_CATEGORY=$listID");
    my $sthProt=$dbh->prepare("SELECT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$listID");
    print $listName.'::';
    $sthProt->execute;
    my $okProt=0;
    while (my ($protID)=$sthProt->fetchrow_array) {
        print ';' if  $okProt;
        print $protID;
        $okProt=1;
    }
    $sthProt->finish;
    $dbh->disconnect;

    print 0 unless $okProt;
    exit;
}

sub ajaxSaveAnnotations {
    print header(-'content-encoding' => 'no',-charset => 'utf-8');
    warningsToBrowser(1);

    ##<Deleting all previous annot for PCA target
    if ($pcaType =~ /^quantif/) {
        $dbh -> do("DELETE FROM ANNOTATIONSET WHERE ID_EXPLORANALYSIS = $explorID AND ANNOT_TYPE LIKE 'property:%'");
    }
    else {
        $dbh -> do("DELETE FROM ANNOTATIONSET WHERE ID_EXPLORANALYSIS = $explorID AND ANNOT_TYPE LIKE 'prot:%'");
    }

    my ($annotSetID) = $dbh -> selectrow_array("SELECT MAX(ID_ANNOTATIONSET) FROM ANNOTATIONSET");
    my $sthInsertAnnot = $dbh -> prepare("INSERT INTO ANNOTATIONSET(ID_ANNOTATIONSET,ID_EXPLORANALYSIS,NAME,ANNOT_RANK,ANNOT_TYPE,ANNOT_LIST) VALUES (?,$explorID,?,?,?,?)");

    if (param('string')) {
        my @annotList;
        foreach my $annotStrg (split(":%:",param('string'))) {
            my %params;
            foreach my $dataStrg (split('//',$annotStrg)) {
                my @data=split('==',$dataStrg);
                $params{$data[0]}=$data[1];
            }
            push @annotList,\%params;
        }
        ##<Proteins
        if ($pcaType=~ /^prot/) {
            foreach my $refAnnot (@annotList) {
                $annotSetID++;
                #<GO annotation
                if ($refAnnot->{type} eq 'prot:GO' || $refAnnot->{type} eq 'prot:PA' ) {
                    $sthInsertAnnot -> execute($annotSetID,$refAnnot->{param},$refAnnot->{position},$refAnnot->{type},$refAnnot->{name}); # '#goAnaID,Aspect,Term' used as name
                }
                #<List annotation
                elsif ($refAnnot->{type} eq 'prot:LIST') {
                    $sthInsertAnnot -> execute($annotSetID,$refAnnot->{param},$refAnnot->{position},$refAnnot->{type},undef); # '#listID' used as name
                }
            }
        }
        ##<Quantifications
        else {
            foreach my $refAnnot (@annotList) {
                $annotSetID++;
                #<Property/treatment
                if ($refAnnot->{type} =~ /property:(.)/) {
                    my $propType=$1;
                    my $annotName;
                    if ($propType eq 'O') {$annotName=$refAnnot->{param};}
                    else { # Treatment
                        my ($propID,$propValue)=split(':=:',$refAnnot->{param});
                        $annotName=$propID.':=:'.&promsQuantif::convertTreatmentToText($propValue);
                    }
                    $sthInsertAnnot -> execute($annotSetID,$annotName, $refAnnot->{position}, $refAnnot->{type}, $refAnnot->{param}); # '#propID:=:value' used as name, Stored also in ANNOT_LIST in case too long for NAME field
                }
            }
        }
    }
    $sthInsertAnnot -> finish;
    $dbh -> commit;
    $dbh -> disconnect;

    print "###OK###";
    exit;
}

sub exportProteins {

    my ($pcaName) = $dbh -> selectrow_array("SELECT NAME from EXPLORANALYSIS where ID_EXPLORANALYSIS = $explorID");
    my $dim=param('dim');
    my $sortOrder=(param('sort')) || 'p-value';
    my %prot2Display;
    foreach my $modProtID (param('chkProt')) {
        $prot2Display{$modProtID}=1;
    }
    my %protData;
	&fetchProteinsFromDimFile($dim,'export',\%protData,\%prot2Display);
    my $numTotProteins=scalar keys %protData;
    my ($workbook, $worksheet, $formatLine, $title, $lineSize);
    my $iLine = my $yTitle = my $yLine = 0;
    $title = 'List proteins'.$pcaName;
    print header(-type => "application/vnd.ms-excel",-attachment => "$title.xls");
    $workbook = Spreadsheet::WriteExcel -> new("-");
    $worksheet = $workbook -> add_worksheet();
    my $formatTitle = $workbook -> add_format(bold=>1, color=>'red', size=>'12', align=>'vcenter');
    $worksheet -> set_row(0,25);
    $worksheet -> merge_range('A1:R1', "$title : $numTotProteins proteins (dimension $dim)", $formatTitle);
    $iLine++;

    my @titleHeader = ('Proteins', 'Gene Name', 'Correlation', 'p-value', 'MW(kDa)', 'Description','Species');
    my $formatHeader = $workbook -> add_format(bold => 1, size => 10, align => 'vcenter', color => 'white', bg_color => 'black');
    foreach my $title (@titleHeader) {
        $worksheet -> write($iLine, $yTitle, $title, $formatHeader);
        $yTitle++;
    }
    $iLine++;

    foreach my $modProtID (sort{&sortSigProteins($sortOrder,\%protData)} keys %protData) {
		#next unless $prot2Display{$modProtID};
        my ($protID,$modStrg)=split('-',$modProtID);
		my ($protLabel,$correl,$pValue,$mw,$des,$species,$refGeneList)=@{$protData{$modProtID}};
        $worksheet->write($iLine,$yLine,$protLabel);
        $yLine++;
		my $geneStrg=(scalar @{$refGeneList}==0)? '-' : join(',',@{$refGeneList});
		$worksheet->write($iLine,$yLine,$geneStrg);
        $yLine++;
        $worksheet->write($iLine,$yLine,$correl);
        $yLine++;
        $worksheet->write($iLine,$yLine,$pValue);
        $yLine++;
        $worksheet->write($iLine,$yLine,$mw);
        $yLine++;
        $worksheet->write($iLine,$yLine,$des);
        $yLine++;
        $worksheet->write($iLine,$yLine,$species);
        $iLine++;
        $yLine=0;
    }
    exit;
}


sub fetchProteinsFromDimFile { # GLOBALS: $dbh,$pathToFile,$isPtmQuantif,%ptmUsed
	my ($dim,$siteDisplayFormat,$refProtData,$refSelList)=@_;
	
	###>Get PTM info if any<###
	# my %quantifModifInfo;
	# if ($isPtmQuantif) {
	# 	my @quantifModifs=keys %ptmUsed;
	# 	&promsQuantif::getQuantifModificationInfo($dbh,\@quantifModifs,\%quantifModifInfo);
	# }

    my %protIdList;
    my $oldSingleModID=(keys %{$ptmUsed{NAME}})[0]; # in case of old single PTM-site format

    my $protFile="$pathToFile/quantifProtDim$dim";
    $protFile.=($pcaType eq 'quantif_sc')? '_sc.txt' : '.txt';
    open(DATA,$protFile);
    while (<DATA>) {
        #next if $.==1;
        chomp;
        my ($modProtID,$correl,$pValue)=split(/\t/,$_);
        next unless $modProtID;
		my ($protID,$modStrg)=($modProtID=~/^(\d+)-?(.*)/); # modCode is null for whole protein quantif
        if ($modStrg && $modStrg !~ /^-?\d+:/) { # potential Free res (-1:...)
			$modStrg=$oldSingleModID.':'.$modStrg;
			$modProtID=$protID.'-'.$modStrg;
		}
		next if ($refSelList && !$refSelList->{$modProtID});
		push @{$protIdList{$protID}},$modStrg; # undef for protein quantif
		$correl=sprintf "%.3f",$correl;
        $pValue=sprintf "%.1e",$pValue;
		@{$refProtData->{$modProtID}}=(undef,$correl,$pValue); # [0] for future aias
    }
    close DATA;
	
	if (scalar keys %protIdList > 0) {
		my $protIdStrg=join(',',keys %protIdList);
		my ($GNidentID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
		my $sthProtInfo=$dbh->prepare("SELECT P.ID_PROTEIN,ALIAS,MW,PROT_DES,ORGANISM,GROUP_CONCAT(DISTINCT MI.VALUE ORDER BY IDENT_RANK SEPARATOR ',')
										FROM PROTEIN P
										LEFT JOIN MASTERPROT_IDENTIFIER MI ON P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND MI.ID_IDENTIFIER=$GNidentID
										WHERE P.ID_PROTEIN IN ($protIdStrg) GROUP BY P.ID_PROTEIN");
		$sthProtInfo->execute;
		while (my ($protID,$alias,$mw,$des,$species,$geneStrg)=$sthProtInfo->fetchrow_array) {
			$alias=~s/ .*//; # Clean protein alias from badly parsed characters in case identifier conversion
			$alias=~s/[,;']/\./g; # Clean MaxQuant crappy contaminant identifiers
			$mw=sprintf "%.1f",$mw/1000;
			my @geneList=($geneStrg)? split(',',$geneStrg) : ();
			foreach my $modCode (@{$protIdList{$protID}}) {
				my $modProtID=$protID;
				my $protLabel=$alias;
				if ($modCode) {
					#my ($formatCode,$displayCode)=&promsQuantif::formatProteinModificationSites($modCode,\%quantifModifInfo,$siteDisplayFormat);
					my $displayCode=&promsQuantif::displayModificationSites($modCode,$ptmUsed{DISPLAY},$siteDisplayFormat);
                    $protLabel.='-'.$displayCode;
					$modProtID.='-'.$modCode;
				}
				$refProtData->{$modProtID}[0]=$protLabel;
				push @{$refProtData->{$modProtID}},($mw,$des,$species,\@geneList);
			}
		}
		$sthProtInfo->finish;
	}
}

sub sortSigProteins {
    my ($sort,$refProtData)=@_;
    if ($sort eq 'protein') {lc($refProtData->{$a}[0]) cmp lc($refProtData->{$b}[0])}
    elsif ($sort eq 'correlation') {$refProtData->{$b}[1]<=>$refProtData->{$a}[1] || $refProtData->{$a}[2]<=>$refProtData->{$b}[2] || lc($refProtData->{$a}[0]) cmp lc($refProtData->{$b}[0])}
    else {$refProtData->{$a}[2]<=>$refProtData->{$b}[2] || $refProtData->{$b}[1]<=>$refProtData->{$a}[1] || lc($refProtData->{$a}[0]) cmp lc($refProtData->{$b}[0])}
}

sub stopOnError {
    my ($error)=@_;
    print qq
|</DIV>
<BR><BR>
<FONT class="title3">$error</FONT>
</BODY>
</HTML>
|;
	exit;
}

sub displayFull3D {
	my $pcaName=param('pcaName');
	my @axisTitles=split(';',param('axisTitles'));
	my $pointLabelStrg=param('pointLabels');
	my $pointDataStrg=param('pointData');
	my $highlightStrg=param('highlight');
	
	####>Processing data<####
	my @pointLabels=split(';',$pointLabelStrg);
	my (%id2Index,@idList,%coord,%traces,%usedMatched,%norm2id);
	foreach my $pointStrg (split(';',$pointDataStrg)) {
		my ($id,$x,$y,$z)=split(',',$pointStrg);
		#$id=~s/%.+//; # remove normalizing quantif if any (WARNING: highlight will not work if based on normalizing quantif)
		($id,my $norm)=split('%',$id);
		$norm2id{$norm}=$id if $norm;
		push @{$coord{x}},$x;
		push @{$coord{y}},$y;
		push @{$coord{z}},$z;
		push @idList,$id;
		$id2Index{$id}=$#{$coord{x}};
		$usedMatched{$id}{0}=1; # default (no highlight)
	}

	##>W/o highlight (default)
	%{$traces{0}}=(name=>'Other',color=>'#555555',matched=>\@idList);
	
	##>With highlight: Keeping only the last highlight if point matches multiple
	foreach my $highlight (split(':%:',$highlightStrg)) {
		my ($name,$pos,$color,$matchStrg)=split('&&',$highlight);
		my @matched;
		foreach my $match (split(';',$matchStrg)) {
			#next unless defined($id2Index{$match}); # eg for GO: corresponding prot IDs are not used in PCA
			my $trueMatch=(defined($id2Index{$match}))? $match : ($norm2id{$match} && defined($id2Index{$norm2id{$match}}))? $norm2id{$match} : undef;
			next unless $trueMatch;
			$usedMatched{$trueMatch}{$pos}=1;
			push @matched,$trueMatch;
		}
		%{$traces{$pos}}=(name=>$name,color=>$color,matched=>\@matched);
	}
	foreach my $match (keys %usedMatched) {
		next if scalar keys %{$usedMatched{$match}}==1;
		my $finalPos;
		foreach my $pos (sort{$b<=>$a} keys %{$usedMatched{$match}}) {
			if ($finalPos) {
				foreach my $i (0..$#{$traces{$pos}{matched}}) {
					if ($traces{$pos}{matched}[$i] eq $match) {
						splice @{$traces{$pos}{matched}},$i,1;
						delete $traces{$pos} unless scalar @{$traces{$pos}{matched}};
						last;
					}
				}
			}
			else {$finalPos=$pos;} # keep this one
		}
	}
	
	####>Starting HTML<####
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>3D PCA</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">

<SCRIPT src="$promsPath{html}/js/plotly/plotly-gl3d-1.47.4.js"></SCRIPT>
<!--
<SCRIPT src="https://cdn.plot.ly/plotly-gl3d-latest.min.js"></SCRIPT>
-->
<SCRIPT type="text/javascript">
    var figure = {
    frames: [],
    layout: {
        title: '3D view of $pcaName',
		titlefont: {
			color: 'black',
			size: 22
		},
        autosize: true,
        dragmode: 'zoom',
        scene: {
            xaxis: {title: '$axisTitles[0]'},
            yaxis: {title: '$axisTitles[1]'},
            zaxis: {title: '$axisTitles[2]'},
            dragmode: 'orbit'
        },
        xaxis: {autorange: true},
        yaxis: {autorange: true},
        zaxis: {autorange: true},
		legend: {
			x: 1,
			y: 1,
			traceorder: 'normal',
			font: {
			  size: 14,
			  color: '#000'
			},
			bgcolor: '#E2E2E2',
			bordercolor: '#FFFFFF',
			borderwidth: 2
		}
    },
	data: [
|;
	my $firstTrace=1;
	foreach my $pos (sort{$a<=>$b} keys %traces) {
		print ",\n" unless $firstTrace;
		
		print qq
|		{
			name: '$traces{$pos}{name}',
            marker: {
				opacity: 1,
				color: '$traces{$pos}{color}',
				colorscale: []
			},
			type: 'scatter3d',
			mode: 'markers+text',
			hoverinfo: 'text',
			hovertext: [|;
		my $count=0;
		#>hoverText
		foreach my $match (@{$traces{$pos}{matched}}) {
			my $idx=$id2Index{$match};
			if ($count) {
				print ',' ;
				print "\n" unless $count % 100;
			}
			print "'$pointLabels[$idx]'";
			$count++;
		}
		print "],\n";
		#>Point coordinates
		foreach my $axis ('x','y','z') {
			print "			$axis: [";
			$count=0;
			foreach my $match (@{$traces{$pos}{matched}}) {
				my $idx=$id2Index{$match};
				if ($count) {
					print ',' ;
					print "\n" unless $count % 100;
				}
				print $coord{$axis}[$idx];
				$count++;
			}
			print ']';
			print ',' if $axis ne 'z';
			print "\n";
		}
		print '		}';
		$firstTrace=0;
	}
	my $showLegendStatus=(scalar keys %traces > 1)? 'true' : 'false';
	print qq
|	]
}
</SCRIPT>
<CENTER>
<!--
<FONT class="title">PCA : <FONT color=#DD0000>$pcaName</FONT></FONT><BR><BR>
-->
<DIV id="PCA_3dDIV" style="width:100%;height:100%;" class="plotly-graph-div"></DIV>
</CENTER>
<SCRIPt type="text/javascript">
(function(){
	window.PLOTLYENV={'BASE_URL': 'https://plot.ly'};

	var gd = document.getElementById('PCA_3dDIV')
	var resizeDebounce = null;

	function resizePlot() {
		var bb = gd.getBoundingClientRect();
		Plotly.relayout(gd, {
			width: bb.width,
			height: bb.height
		});
	}
	
	window.addEventListener('resize', function() {
		if (resizeDebounce) {
			window.clearTimeout(resizeDebounce);
		}
		resizeDebounce = window.setTimeout(resizePlot, 100);
	});
	
	Plotly.plot(gd, {
		data: figure.data,
		layout: figure.layout,
		frames: figure.frames,
		config: {displayModeBar: true,showlegend: $showLegendStatus} // displaylogo: false, linkText: 'Export to plot.ly',showLink: true
	});

}());
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

####>Revision history<####
# 1.2.21 [BUGFIX] Minor fix in display format used in &promsQuantif::displayModificationSites (PP 28/04/21)
# 1.2.20 [UPDATE] Compatibility with Free residues quantifications (PP 04/12/20)
# 1.2.19 [MODIF] Adapt query on PATHWAY_ANALYSIS to avoid GSEA entries (VL 18/11/20)
# 1.2.18 [BUGFIX] Correct calls to script runAndDisplayPathwayAnalysis replaced by displayPathwayAnalysis (VL 18/11/20)
# 1.2.17 [BUGFIX] Fixed truncated GROUP_CONCAT (PP 11/09/20)
# 1.2.16 [ENHANCEMENT] Better AJAX calls chaining during saved highlightings restoration (PP 03/09/20)
# 1.2.15 [CHANGE] Modified path to plotly script since it was changed to allow PCA 3D displaying when using a VPN browser access (VS 21/08/20)
# 1.2.14 [CHANGE] Code change to use new &promsQuantif::getQuantifModificationInfo (PP 02/07/20)
# 1.2.13 [BUGFIX] Remove forgotten temporary dev command (PP 27/04/20)
# 1.2.12 [BUGFIX] in 3D view highlighting for site quantifs with protein-level normalization (PP 22/04/20)
# 1.2.11 [ENHANCEMENT] Clean ' and " from strings given to JS (PP 23/03/20)
# 1.2.10 [CHANGE] Max # dimensions set to 10 regardless of contribution (PP 22/03/20)
# 1.2.9 [UPDATE] Changed RANK field to (ANNOT/IDENT)_RANK for compatibility with MySQL 8 (PP 04/03/20) 
# 1.2.8 [BUGFIX] in SQL query leading to truncated list of TARGET_POS when too many in same quantif (PP 22/02/20)
# 1.2.7 [CHANGE] Removed protein/site PCAs options and values distribution chart (PP 21/02/20)
# 1.2.6 [ENHANCEMENT] Better management of quantification parameters for ajaxListProt (PP 07/01/20)
# 1.2.5 [BUGFIX] in SQL query requiring GROUP_CONCAT for TARGET_POS & [UX] Smooth scroll to list of proteins (PP 20/11/19)
# 1.2.4 [ENHANCEMENT] Multi-site support (PP 28/10/19)
# 1.2.3 [ENHANCEMENT] Integrates R generated density plot (count normalized) as quality representation of imputed quantification values (VS 14/10/19)
# 1.2.2 Used local version of plotly-gl3d-1.47.4.min.js (PP 24/05/19)
# 1.2.1 check number of dimension to skip warning (SL 11/01/19)
# 1.2.0 Added full 3D view with plotly.js & avoid asynchronous ajax (PP 17/12/18)
# 1.1.2 Minor improvement in highlighting annotations sorting (PP 21/06/18)
# 1.1.1 Compatible with sites list (PP 14/12/17)
# 1.1.0 Minor bug fix on selection of max number dimensions (PP 18/01/17)
# 1.0.9 display quantifName for PCA with normalization (SL 19/07/16)
# 1.0.8 3rd dimension projected on point size (PP 09/11/15)
# 1.0.7 auto extended selection + add gene name + export excel for list in dim (SL 16/05/15)
# 1.0.6 Compatible with PTM quantifications (PP 16/06/15)
# 1.0.5 display modification data and quantification reference (SL 20/03/15)
# 1.0.4 add pathway annotation and manage bin for go analysis (SL 24/11/14)
# 1.0.3 JS ajaxManageSaveProteins moved to promsMod.pm (PP 22/08/14)
# 1.0.2 Quantification + protein highlithing & synchronized with new data file naming (PP 08/08/14)
# 1.0.1 Code cleaning & simplication. (PP 10/07/14)
# 1.0.0 New script to display and store the results of PCA (SL 08/04/14)
