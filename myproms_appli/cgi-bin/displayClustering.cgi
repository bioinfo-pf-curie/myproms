#!/usr/local/bin/perl -w

################################################################################
# displayClustering.cgi       1.2.6                                            #
# Authors: P. Poullet, S.Liva (Institut Curie)                                 #
# Contact: myproms@curie.fr                                                    #
# display and store the results of Clustering analysis        		           #
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

#######################
####>Configuration<####
#######################
#print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my ($explorID,$experimentID) = &promsMod::cleanNumericalParameters(param('explorID'),param('experimentID'));
my $ajax = param('AJAX') || '';
my $view=param('VIEW') || "1D";
my ($sel1D, $sel2D)= ($view eq "1D")? (" selected", "") : (""," selected");

####>Connexion to DB<####
my $dbh = &promsConfig::dbConnect;
my $projectID=param('PROJECT_ID') || &promsMod::getProjectID($dbh,$explorID,'EXPLORANA');
my $pathToFile = "$promsPath{explorAna}/project_$projectID/$explorID";

my $protAnnotMatchStrg=''; # required for ajax protein annotation calls
my $quantifListStrg='';
my %quantificationIDs;
my $isPtmQuantif=0;
my %ptmUsed;
my $sthSelExplorQuantif=$dbh->prepare("SELECT Q.ID_QUANTIFICATION,GROUP_CONCAT(DISTINCT TARGET_POS ORDER BY TARGET_POS SEPARATOR '_'),CONCAT(COALESCE(Q.ID_MODIFICATION,'0'),',',COALESCE(GROUP_CONCAT(DISTINCT MQ.ID_MODIFICATION SEPARATOR ','),'0'))
								FROM QUANTIFICATION Q
								LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
								INNER JOIN EXPLORANA_QUANTIF EQ ON EQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION
								WHERE ID_EXPLORANALYSIS=? GROUP BY Q.ID_QUANTIFICATION");
$sthSelExplorQuantif -> execute($explorID);
while (my($quantifID,$targetPos,$ptmStrg) = $sthSelExplorQuantif -> fetchrow_array) {
    $quantifListStrg.=':' if $quantifListStrg;
    $quantifListStrg.=$quantifID."_".$targetPos;
    $quantificationIDs{$quantifID}=1;
    next if $ptmStrg eq '0,0'; # no PTM quantif
    foreach my $modID (split(',',$ptmStrg)) {
        $ptmUsed{$modID}=1 if $modID;
    }
}
$sthSelExplorQuantif -> finish;

my $oldSingleModID=0; # to handle old single-site format stored in files (no modifID)
if (scalar keys %ptmUsed) {
    $isPtmQuantif=1;
	$oldSingleModID=(keys %ptmUsed)[0]; # just in case PTM
    $protAnnotMatchStrg=",'^###(\$|-)'";
}

##AJAX call
if ($ajax eq 'propAnnotateClustering') {
    &ajaxPropAnnotateClustering;
    exit;
}
elsif ($ajax eq 'themeAnnotateClustering') {
    &ajaxThemeAnnotateClustering;
    exit;
}
elsif ($ajax eq 'goAnnotateClustering') {
    $dbh->disconnect;
    &ajaxGoAnnotateClustering;
    exit;
}
elsif ($ajax eq 'saveAnnotations') {
    &ajaxSaveAnnotations;
    exit;
}
elsif ($ajax eq 'pathwayAnnotateClustering') {
    &ajaxPathwayAnnotateClustering;
    exit;
}
elsif ($ajax eq 'getListTheme') {
    &ajaxGetListTheme;
    exit;
}
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
my $disabSave=($projectAccess eq 'guest')? ' disabled' : '';

my ($clusteringName,$paramList,$filterList) = $dbh->selectrow_array("SELECT NAME,PARAM_LIST,FILTER_LIST FROM EXPLORANALYSIS WHERE ID_EXPLORANALYSIS = $explorID");
my ($quantifFam,$quantifMeasCode,$quantifMeasName)=('','',''); # for 2D only
if ($view eq '2D') {
	my %proteinQuantifFamilies=&promsQuantif::getProteinQuantifFamilies;
	($quantifFam)=$paramList=~/quantFam=(\w+)/;
	if ($paramList=~/quantCode=(\w+)/) {
		$quantifMeasCode=$1;
		foreach my $refMeasInfo (@{$proteinQuantifFamilies{MEASURES}{$quantifFam}}) {
			if ($refMeasInfo->[0] eq $quantifMeasCode) {
				$quantifMeasName=$refMeasInfo->[1];
				last;
			}
		}
	}
	else {$quantifMeasName=$proteinQuantifFamilies{NAME}{$quantifFam};}
}

my $pepType=($filterList && $filterList=~/\/PEP=(\w+)/)? $1 : ''; # needed in ajaxListSelectedProteins JS function
$pepType=($pepType eq 'NUM_PEP_USED')? 'all' : ($pepType eq 'DIST_PEP_USED')? 'distinct' : $pepType;

my %bioSampProperties;
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


my (%listThemes, %goAnalyses, %parentGoAna, %pathwayAnalyses);
if ($view eq '2D') {

    ##THEMES
    my $sthTh=$dbh->prepare("SELECT ID_CLASSIFICATION,NAME FROM CLASSIFICATION WHERE ID_PROJECT=$projectID");
    $sthTh->execute;
    while (my($themeID,$themeName)=$sthTh->fetchrow_array) {
        $listThemes{$themeID}=$themeName;
    }
    $sthTh->finish;

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
    my $sthPA=$dbh->prepare("SELECT ID_PATHWAY_ANALYSIS, ID_CATEGORY, NAME, PARAM_STRG FROM PATHWAY_ANALYSIS WHERE ID_EXPERIMENT=$experimentID AND STATUS>0");
    $sthPA->execute;
    while (my ($pathID, $catID, $name, $paramStrg)=$sthPA->fetchrow_array) {
        $catID = (!$catID)? "" : $catID;
        @{$pathwayAnalyses{$pathID}}=($catID, $name, $paramStrg);
    }
    $sthPA->finish;
}

#######################
####>Starting HTML<####
#######################
my $hmRowItem=($isPtmQuantif)? 'Site' : 'Protein';
my $hmRowItems=$hmRowItem.'s';
$quantifMeasName=~s/Protein/$hmRowItem/ if $view eq '2D'; # undef otherwise
print header(-'content-encoding'=>'no',-charset=>'utf-8');
warningsToBrowser(1);
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
print qq
|<HTML>
<HEAD>
<TITLE>Clustering</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.TD {font-weight:normal;}
.TH {font-weight:bold;}
.row_0{background-color:$darkColor;}
.row_1{background-color:$lightColor;}
.annotation{width:300px;}
.popup {z-index:999;background-color:#FFFFFF;border:solid 3px #999999;padding:5px;position:absolute;display:none;box-shadow:10px 10px 20px #808080;}
</STYLE>
<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/heatMap.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/peptidePlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
<SCRIPT language="Javascript">
|;
&promsMod::popupInfo();
if ($view eq '2D') {
	if ($quantifFam eq 'RATIO') {
		print qq
|function displayCellInfo(cellValue,hCell,rangeMinValue,rangeMaxValue) {
    var dispValue;
    if (cellValue==null) {dispValue='* no value *';}
    else {
        var overStrg=(cellValue <= rangeMinValue)? ' ≤ ' : (cellValue >= rangeMaxValue)? ' ≥ ' : '=';
        dispValue=Math.exp(cellValue*Math.log(2));
        if (dispValue < 1) {dispValue='1/'+(1/dispValue).toPrecision(2);}
        else {dispValue=dispValue.toPrecision(2);}
    }
    return '$quantifMeasName'+overStrg+dispValue; // Fold change
}
|;
	}
	else { # non-ratio quantif
		print qq
|function displayCellInfo(cellValue,hCell,rangeMinValue,rangeMaxValue) {
    var dispValue=(cellValue==null)? '* no value *' : Math.round(Math.exp(cellValue*Math.log(10)));
    return '$quantifMeasName='+dispValue;
}
|;
	}
	print qq
|function listProteins(selectedProteins) {
    ajaxListSelectedProteins(selectedProteins.join(','),'protein');
}

function ajaxListSelectedProteins(selectedPoints,sort='protein') { // 2D view only
    saveListGlobals.themesFetched=false; // resets ajaxManageSaveProteins mecanism
    var listDiv=document.getElementById('protListDIV');
    listDiv.innerHTML="<BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<FONT class=\\"title3\\">Fetching data. Please wait...</FONT><BR>&nbsp;&nbsp;&nbsp;&nbsp;<IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR>";
    listDiv.style.display='';
    listDiv.scrollIntoView({block:"start",inline:"nearest",behavior:"smooth"});
    var paramStrg="AJAX=ajaxListProt&CALL=PCA&id_project=$projectID&ACT=results&quantifFamily=$quantifFam&dispMeasure=$quantifMeasCode&pepType=$pepType&quantifList=$quantifListStrg&sort="+sort+"&selProt="+selectedPoints; // quantifFocus set to any modID
    var XHR = getXMLHTTP();
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
// -----> Block copied from compareQuantifications.cgi to complete selected proteins list functionalities
var isNav = (navigator.appName.indexOf("Netscape") !=-1);
function ajaxProteinData(e,protID,refData,action,midX,botY,divNum) { // action=ajaxPepRatio,ajaxProtStat
	var quantifData=(typeof(refData)=='object')? refData[0] : refData;
	var displayDiv=document.getElementById('displayDIV');
	var divX,divY;
	if (e) { // called from protein list
		divX = (isNav)? e.pageX : event.clientX + document.body.scrollLeft; //divX-=5;
		divY = (isNav)? e.pageY : event.clientY + document.body.scrollTop; //divY+=10;
	}
	else {
		var graphDiv=document.getElementById('mainGraphDIV');
		divX = getDomOffset(graphDiv,'offsetLeft') + graphDiv.offsetWidth +2;
		if (midX) { // called from heatmap
			var hmDiv=(HM.editMenu)? document.getElementById('mainGraphDIV_canvas') : document.getElementById('mainGraphDIV');
			divY= Math.max(getDomOffset(graphDiv,'offsetTop'),getDomOffset(hmDiv,'offsetTop') + botY -200);
		}
		else {divY = getDomOffset(graphDiv,'offsetTop');}
	}
	displayDiv.style.left = divX + 'px';
	displayDiv.style.top = divY + 'px';
	displayDiv.style.display='block';

	var infoDiv;
	if (!divNum \|\| divNum==1) {
		infoDiv=document.getElementById('infoDIV');
		document.getElementById('infoDIV2').innerHTML='';
	}
	else {infoDiv=document.getElementById('infoDIV2');}
	infoDiv.innerHTML="<BR><BR><IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR>";

	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	var qData=quantifData.split('_'); // id_quantif,ratio
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT="+action+"&id_quantif="+qData[0]+"&ratio="+qData[1]+"&id_prot="+protID+"&pepType=$pepType",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			if (action=='ajaxPepRatio') {
				var resCode;
				if (divNum && divNum==2) {
					resCode=XHR.responseText.replace(/pepPlot/g,'pepPlot2');
					resCode=resCode.replace(/pepList/g,'pepList2');
				}
				else {resCode=XHR.responseText;}
				var codeParts=resCode.split('#==========#');
				infoDiv.innerHTML=codeParts[0]; // HTML part
				eval(codeParts[1]); // javascript part
				// 2nd pep ratio (called from log2 plot)
				if (typeof(refData)=='object') {
					ajaxProteinData(null,protID,refData[1],'ajaxPepRatio',null,null,2);
				}
			}
			else {infoDiv.innerHTML=XHR.responseText;}
		}
	}
	XHR.send(null);

}
// <----- End of block

function selectSort(newSort,ajax) { // view=list (true or from Volcano ajax)
    var checkedProt=document.protForm.chkProt;
    var chkList=[];
    if (checkedProt.length) {
        for (let i=0; i<checkedProt.length; i++) {chkList.push(checkedProt[i].value);}
    }
    else {chkList.push(checkedProt[i].value);}
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
}
print qq
|
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

/******************** CLUSTERING ANNOTATION ********************/
var annotSetObject={quantif:{},prot:{}};
var annotRank = {quantif:0,prot:0};
var annotHasChanged=0;
function deleteAnnotObject (valType, name) { //valType => column / row
    annotHasChanged=1;
    if (valType == 'row') {
        delete annotSetObject.prot[name];
    }
    else {
        delete annotSetObject.quantif[name];
    }
    var selButton = document.getElementById('saveButton');
    selButton.style.display = '';
    return(true);
}

function selectAnnotTarget(target) {
    if (target=='quantif') {
        document.getElementById('protAnnotationType').style.display='none';
        document.getElementById('protAnnotationType').selectedIndex=0;
        document.getElementById('quantifAnnotationType').style.display='';
        document.getElementById('listDisplayDIV').style.display='none';
    }
    else { // prot
        document.getElementById('quantifAnnotationType').style.display='none';
        document.getElementById('quantifAnnotationType').selectedIndex=0;
        document.getElementById('protAnnotationType').style.display='';
        document.getElementById('listDisplayDIV').style.display='none';
    }
    document.getElementById('goAspect').style.display='none';
    document.getElementById('termDisplayDIV').style.display='none';
}

/*** QUANTIF/COND ANNOTATION ***/
// BIOSAMPLE PROPERTY / TREATMENT-BASED ANNOTATION
function ajaxPropAnnotateClustering(selProp,propName,saved) {
    if (!selProp) {
       return;
    }
    if (annotSetObject.quantif[propName]) {
        alert('ERROR: "'+propName+'" already used!');
        return;
    }
    if (!saved) saved=false;

    var paramStrg='AJAX=propAnnotateClustering&explorID=$explorID&experimentID=$experimentID&qList=$quantifListStrg&property='+selProp;

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
    XHR.open("POST","$promsPath{cgi}/displayClustering.cgi",!saved); // Switches to synchronous for already saved nnotations
	//Send the proper header information along with the request
    XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
    XHR.onreadystatechange=function() {
        if (XHR.readyState == 4 && XHR.responseText) {
            if (eval(XHR.responseText)) {
                var selData=selProp.split(':'); // property:O:propID
                annotRank.quantif++;
                annotSetObject.quantif[propName] = {type:selData[0]+':'+selData[1],param:'#'+selData[2],position:annotRank.quantif}; // # for ID tag
                if (!saved) document.getElementById('saveButton').style.display = '';
            }
        }
    }
    XHR.send(paramStrg);
}

/*** PROTEIN ANNOTATION ***/
function updateProtAnnotationType(protAnnotSel) {
    var typeID=protAnnotSel.value;
    if (!typeID) {
        document.getElementById('goAspect').style.display='none';
        document.getElementById('termDisplayDIV').style.display='none';
        document.getElementById('listDisplayDIV').style.display='none';
        return;
    }
    var typeInfo=typeID.split(':');
    if (typeInfo[0]=='GOA') {
        document.getElementById('listDisplayDIV').style.display='none';
        updateGoAspects(typeInfo[1]);
    }
    else if (typeInfo[0]=='TH') {
        document.getElementById('goAspect').style.display='none';
        document.getElementById('termDisplayDIV').style.display='none';
        document.getElementById('listDisplayDIV').style.display='none';
        ajaxGetListTheme(typeInfo[1]);
        //ajaxThemeAnnotateClustering(typeInfo[1],protAnnotSel.options[protAnnotSel.selectedIndex].text,false);
    }
    else if (typeInfo[0]=='PA') {
        document.getElementById('goAspect').style.display='none';
        document.getElementById('termDisplayDIV').style.display='none';
        document.getElementById('listDisplayDIV').style.display='block';
        ajaxGetPathwayList(typeInfo[1]);
    }
}

function selectBoxes(type, val) {
    var listSel=document.getElementsByName('list');
    var themeName
    var numChkList=0;
    if (type == 'TH') {
        var themeName = document.getElementById('protAnnotationType').options[document.getElementById('protAnnotationType').selectedIndex].text;
        var thCheck= (listSel[0].checked == true)? true : false;
        numChkList= (thCheck)? numChkList++ : 0;
        document.getElementById('annotListName').value = (listSel[0].checked == true)? themeName : '';
        document.getElementById('annotListName').disabled = (listSel[0].checked == true)? true : false;
        for (var i=1; i<listSel.length; i++) {
            var strgBoxes=listSel[i].value.split(':');
            listSel[i].checked= false;
            listSel[i].disabled= (thCheck)? true : false;
        }
    }
    else {
        for (var i=1; i<listSel.length; i++) {
            if (listSel[i].checked) {
                numChkList++;
            }
        }
        listSel[0].checked = false;
        listSel[0].disabled = (numChkList)? true : false;
    }

    document.getElementById('saveButton').style.display= (!numChkList && annotHasChanged)? 'block' : 'none';

}


//AJAX LIST THEME
function ajaxGetListTheme(themeID) {
    document.getElementById('listDisplayDIV').style.display='block';
    var listDiv=document.getElementById('listDisplayDIV');
    listDiv.innerHTML='';

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
    XHR.open("GET","$promsPath{cgi}/displayClustering.cgi?AJAX=getListTheme&explorID=$explorID&experimentID=$experimentID&PROJECT_ID=$projectID&themeID="+themeID,true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState==4 && XHR.responseText) {
            listDiv.innerHTML=XHR.responseText;
            //document.getElementById('annotList').style.display='block';
            document.getElementById('saveList').style.display='block';
        }
    }
    XHR.send(null);
}

function userListAnnotateClustering() {

    if (!document.getElementById('annotListName').value) {
        alert('ERROR: Provide a name for annotation!');
        return;
    }

    var listSel=document.getElementsByName('list');
    var numList=0;
    var listArray=new Array();
    for (var i=0; i<listSel.length; i++) {
        if (listSel[i].checked) {
            numList++;
            listArray.push(listSel[i].value);
        }
    }
    if (numList == 0) {
        alert('ERROR: Select at least 1 list');
        return;
    }

    var listID=new Array();
    var type;
    for (var list in listArray) {
        var strgList=listArray[list].split(':');
        if (strgList[0] == 'TH') {
            type = 'TH';
            annotHasChanged=1;
            ajaxThemeAnnotateClustering(strgList[1],strgList[2],false,'THEME');
        }
        else if (strgList[0] == 'LI') {
            type='LI';
            listID.push(strgList[1]);
        }
    }
    if (type == 'LI'){
        var listAnnot=(listID.length>1)? listID.join(':') : listID;
        annotHasChanged=1;
        ajaxThemeAnnotateClustering(listAnnot,document.getElementById('annotListName').value,false,'LIST');
    }

}

//AJAX PATHWAY ANNOTATION
function ajaxGetPathwayList(pathID) {
    document.getElementById('termDisplayDIV').style.display='block';
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
    XHR.open("GET","$promsPath{cgi}/runAndDisplayPathwayAnalysis.cgi?AJAX=ajaxGetPathway&FROM=cluster&ID="+pathID,true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState==4 && XHR.responseText) {
            termsDiv.innerHTML=XHR.responseText;
            var pathSel=document.getElementById('pathTermsSELECT');
            pathSel.multiple=true;
            pathSel.className = "";
            pathSel.style.height='200px';
            pathSel.remove(0); // removes "-= Select a GO term =-" option
            termsDiv.style.display='';
            document.getElementById('annotTermsDIV').innerHTML='<INPUT type="button" value="Annotate" onclick="userPathwayAnnotateClustering('+pathID+')"/>';
        }
    }
    XHR.send(null);
}
function ajaxGetPathwayProteinsList(name, termsStrg, pathID, saved) {
    if (annotSetObject.prot[name]) {
        alert('ERROR: "'+name+'" is already used!');
        return;
    }
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
    XHR.open("GET","$promsPath{cgi}/displayClustering.cgi?AJAX=pathwayAnnotateClustering&explorID=$explorID&experimentID=$experimentID&PROJECT_ID=$projectID&annotName="+encodeURIComponent(name)+"&pathID="+pathID+"&termsStrg="+termsStrg,true);

    XHR.onreadystatechange=function() {
        if (XHR.readyState==4 && XHR.responseText) {
            if (eval(XHR.responseText)) {
               annotRank.prot++;
               annotSetObject.prot[name] = {type:'prot:PA:M',param:'#'+pathID+':=:'+termsStrg,position:annotRank.prot}; // # for ID tag
               if (!saved) document.getElementById('saveButton').style.display = '';
            }
        }
    }
    XHR.send(null);
}
function userPathwayAnnotateClustering(pathID){
    var pathSel=document.getElementById('pathTermsSELECT');
    var termsArray=new Array();
    var numPath=0;
    for (var i=0; i<pathSel.options.length; i++) {
        if (pathSel.options[i].selected) {
            numPath++;
            var termStrg=pathSel.options[i].value+'\@'+pathSel.options[i].text;
            termsArray.push(termStrg);
        }
    }
    if (numPath == 0) {
        alert('ERROR: Select at least 1 pathway term!');
        return;
    }
    if (!document.getElementById('annotName').value) {
        alert('ERROR: Provide a name for annotation!');
        return;
    }

    var termsPathAnnot=termsArray.join(':\@:');
    var annotName=document.getElementById('annotName').value;
    ajaxGetPathwayProteinsList(annotName,termsPathAnnot,pathID,false);
}

// GO-BASED ANNOTATION
var goAspects=new Object();
|;
foreach my $goID (keys %goAnalyses) {
    print "goAspects[$goID]='$goAnalyses{$goID}[1]';\n";
}
print qq
|
function updateGoAspects(goID) {
    var selAspect=document.getElementById('goAspect');
    selAspect.options.length=1;
    //document.getElementById('selTermsDIV').innerHTML='';
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
    XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxGoTerms&projectID=$projectID&onChange=void&goStrg="+goIdStrg,true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState==4 && XHR.responseText) {
            termsDiv.innerHTML=XHR.responseText;
            var goSel=document.getElementById('goTermsSELECT');
            goSel.multiple=true;
            goSel.className = "";
            goSel.style.height='200px';
            goSel.remove(0); // removes "-= Select a GO term =-" option
            document.getElementById('termDisplayDIV').style.display='';
            document.getElementById('annotTermsDIV').innerHTML='<INPUT type="button" value="Annotate" onclick="userGoAnnotateClustering()"/>';
        }
    }
    XHR.send(null);
}

function userGoAnnotateClustering() {
    var goSel=document.getElementById('goTermsSELECT');
    var goID,goAspect,termsStrg='';
    for (var i=0; i<goSel.options.length; i++) {
        if (goSel.options[i].selected) {
            var termData=goSel.options[i].value.split(',');
            goID=termData[0];
            goAspect=termData[1];
            var termText=goSel.options[i].text.replace(/ \\[GO.+/,"");
            if (termsStrg) termsStrg+=':\@:';
            termsStrg+=termData[2]+'\@'+termText;
        }
    }
    if (!termsStrg) {
        alert('ERROR: Select at least 1 GO term!');
        return;
    }
    if (!document.getElementById('annotName').value) {
        alert('ERROR: Provide a name for annotation!');
        return;
    }
    var annotName=document.getElementById('annotName').value;
    ajaxGoAnnotateClustering(annotName,goID,goAspect,termsStrg,false);
}

function ajaxGoAnnotateClustering(annotName,goID,goAspect,termsStrg,saved) {
    if (annotSetObject.prot[annotName]) {
        alert('ERROR: "'+annotName+'" is already used!');
        return;
    }
    if (!saved) saved=false;

    var paramStrg='AJAX=goAnnotateClustering&explorID=$explorID&experimentID=$experimentID&PROJECT_ID=$projectID&annotName='+encodeURIComponent(annotName)+'&goID='+goID+'&goAspect='+goAspect+'&termsStrg='+termsStrg;

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
    XHR.open("POST","$promsPath{cgi}/displayClustering.cgi",!saved); // Switches to synchronous for already saved annotations
    //Send the proper header information along with the request
    XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
    XHR.onreadystatechange=function() {
        if (XHR.readyState == 4 && XHR.responseText) {
            if (eval(XHR.responseText)) {
                annotRank.prot++;
                annotSetObject.prot[annotName] = {type:'prot:GO:M',param:'#'+goID+':=:'+goAspect+':=:'+termsStrg,position:annotRank.prot}; // # for ID tag
                if (!saved) document.getElementById('saveButton').style.display = '';
            }
        }
    }
    XHR.send(paramStrg);
}

// LIST-BASED ANNOTATION
function ajaxThemeAnnotateClustering(themeID,themeName,saved,type) {

    if (annotSetObject.prot[themeName]) {
        alert('ERROR: "'+themeName+'" is already used!');
        return;
    }
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
    XHR.open("GET","$promsPath{cgi}/displayClustering.cgi?AJAX=themeAnnotateClustering&experimentID=$experimentID&explorID=$explorID&name="+themeName+"&type="+type+"&themeID="+themeID,!saved);
    XHR.onreadystatechange=function() {
        if (XHR.readyState==4 && XHR.responseText) {
            if (eval(XHR.responseText)) {
                annotRank.prot++;
                annotSetObject.prot[themeName] = {type:'prot:'+type,param:'#'+themeID,position:annotRank.prot}; // # for ID tag
                if (!saved) document.getElementById('saveButton').style.display = '';
            }
        }
    }
    XHR.send(null);
}
function selectClustDim (dim) {
    window.location="./displayClustering.cgi?VIEW="+dim+"&explorID=$explorID&experimentID=$experimentID&PROJECT_ID=$projectID";
}

//TODO: Save annotations
function ajaxSaveAnnotations() {
    //Object to string
    var objString = '';
    var concat = 0;
    for (var dim in annotSetObject) {
        for (var i in annotSetObject[dim]) {
            concat++;
            if (concat > 1) objString += ':%:'; // object separator
            objString += 'name=='+i;
            for (var j in annotSetObject[dim][i]) {
                objString += '//'+j+'=='+annotSetObject[dim][i][j];
            }
        }
    }

    var paramStrg="AJAX=saveAnnotations&experimentID=$experimentID&explorID=$explorID&string="+objString;
    //If XHR object already exists, the request is canceled & the object is deletsaveButtoned
    if (XHR && XHR.readyState != 0) {
        XHR.abort();
        delete XHR;
    }
    //Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if (!XHR) {
        return false;
    }
    XHR.open("POST","$promsPath{cgi}/displayClustering.cgi",true);
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
                alert('ERROR: Annotations could not be saved !');
            }
        }
    }
    XHR.send(paramStrg);
}

function toggleValueDistribution() {
	var graph = document.getElementById('valueDistribution');
	var isDisplayed = (graph.style.display == '');
	graph.style.display = (isDisplayed) ? 'none' : '';
	
	
	var displayButton = document.getElementById('valueDistributionButton');
	displayButton.innerHTML = (isDisplayed) ? 'Display values distribution' : 'Hide values distribution';
}

</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Clustering&nbsp;<SELECT class="title2" name="viewDim" onchange="selectClustDim(this.value)"><OPTION value="1D"$sel1D>1D</FONT></OPTION><OPTION value="2D"$sel2D>2D</OPTION></SELECT> : <FONT color=#DD0000>$clusteringName</FONT></FONT><BR>
<BR>
|;

print qq |
<DIV id="waitDIV">
<BR><BR><FONT class="title3">Fetching data. Please wait...</FONT><BR>&nbsp;&nbsp;&nbsp;&nbsp;<IMG src="$promsPath{images}/scrollbarGreen.gif"><BR><BR>
</DIV>
</CENTER>
<DIV id="resultDIV" style="visibility:hidden"> <!-- "display:none" is not compatible with SVG text drawing with Raphaël in Chrome! -->
|;

my (@proteinList, @proteinTab); #@proteinOrder,  @protein, , @proteinClust, @proteinNoOrder
my (@quantifOrder, @quantifInfo, @quantifTab); # @quantifList, , @quantifClust, @quantifNoOrder
my (%quantifNames, %proteinDendro, %proteinNames, %quantifDendro, %quantifProteinValues, %missingValues); # %quantifList, %quantifPrimaryName,

##>Quantifs/Conditions
%quantifNames=&promsQuantif::getDistinctQuantifNames($dbh,$quantifListStrg);

open (QUANTIF_ORDER,"$pathToFile/quantifCluster.txt") || &stopOnError("ERROR!<BR>$!: '$pathToFile/quantifCluster.txt'");
while (my $lineQuantif = <QUANTIF_ORDER>) {
    chomp($lineQuantif);
    next if ($. == 1);
    my ($index, $id, $quantif) = split("\t", $lineQuantif);
    my $strgQuantif;
    if ($quantif=~/%/) {
        my($testQuantif, $refQuantif)=split(/%/,$quantif);
        #my $aliasTest=(split(/:/,$quantifNames{'OPTIMAL'}{$testQuantif}))[0];
        (my $aliasTest =$quantifNames{'OPTIMAL'}{$testQuantif})=~s/.*://;
        #my $aliasRef=(split(/:/,$quantifNames{'OPTIMAL'}{$refQuantif}))[0];
        (my $aliasRef=$quantifNames{'OPTIMAL'}{$refQuantif})=~s/.*://;
        $strgQuantif=$aliasTest." % ".$aliasRef;
    }
    else {
        $strgQuantif=$quantifNames{OPTIMAL}{$quantif};
        #$strgQuantif=$quantifNames{OPTIMAL}{$quantif}  || $quantif;
        #$strgQuantif.='-'.$quantif;
    }
    push @quantifOrder, $quantif;
    push @quantifInfo, "['$strgQuantif','$quantif']";
    #push @quantifInfo, "['$quantifNames{OPTIMAL}{$quantif}','$quantif']";
    $quantifDendro{$id} = $quantif;
}
close (QUANTIF_ORDER);

open (QUANTIF_DENDRO,"$pathToFile/quantifDendro.txt") || &stopOnError("ERROR!<BR>$!: '$pathToFile/quantifDendro.txt'");
while (my $lineQuantifDendro = <QUANTIF_DENDRO>) {
    next if ($. == 1);
    chomp($lineQuantifDendro);
    my ($index, $dend1, $dend2, $height) = split("\t", $lineQuantifDendro);
    if ($dend1 =~ /-/){
        $dend1 =~ s/-//;
        $dend1 = "-".$quantifDendro{$dend1};
    }
    if ($dend2 =~ /-/){
        $dend2 =~ s/-//;
        $dend2 = "-".$quantifDendro{$dend2};
    }
    push @quantifTab, join(",", $index, $dend1, $dend2, $height);
}
close (QUANTIF_DENDRO);

##>Proteins
if ($view eq '2D') {
	
	###>Get PTM info if any<###
	my %quantifModifInfo;
	if ($isPtmQuantif) {
		my @quantifModifs=keys %ptmUsed;
		&promsQuantif::getQuantifModificationInfo($dbh,\@quantifModifs,\%quantifModifInfo);
	}

    my %protIdList;
    #my $oldSingleModID=(keys %ptmUsed)[0]; # in case of old single PTM-site format
    open (PROTEIN_ORDER,"$pathToFile/protCluster.txt") || &stopOnError("ERROR!<BR>$!: '$pathToFile/protCluster.txt'");
    while (my $lineProtein = <PROTEIN_ORDER>) {
        chomp($lineProtein);
        next if ($. == 1);
        my ($index, $pos, $modProtID) = split("\t", $lineProtein);
		my ($protID,$modCode)=$modProtID=~/^(\d+)-*(.*)/;
		if ($modCode && $modCode !~ /^\d+:/) { # no starting modifID => old single PTM format
			$modCode=$oldSingleModID.':'.$modCode;
			$modProtID=$protID.'-'.$modCode;
		}
        push @proteinList,$modProtID;
		push @{$protIdList{$protID}},$modCode; # undef for protein quantif
        $proteinDendro{$pos} = $modProtID;
    }
    close PROTEIN_ORDER;
    
	my $protIdStrg=join(',',keys %protIdList);
	my $sthProt=$dbh->prepare("SELECT ID_PROTEIN,ALIAS FROM PROTEIN WHERE ID_PROTEIN IN ($protIdStrg)");
	$sthProt->execute;
	while (my ($protID,$alias)=$sthProt->fetchrow_array) {
		foreach my $modCode (@{$protIdList{$protID}}) {
			if ($modCode) {
				my ($formatCode,$displayCode)=&promsQuantif::formatProteinModificationSites($modCode,\%quantifModifInfo,'text');
				$proteinNames{"$protID-$modCode"}=$alias.'-'.$displayCode;
			}
			else {$proteinNames{$protID}=$alias;}
		}
	}
	$sthProt->finish;
	
    open (PROTEIN_DENDRO,"$pathToFile/protDendro.txt") || &stopOnError("ERROR!<BR>$!: '$pathToFile/protDendro.txt'");
    while (my $lineProteinDendro = <PROTEIN_DENDRO>) {
        next if $. == 1;
        chomp($lineProteinDendro);
        my ($index, $dend1, $dend2, $height) = split("\t", $lineProteinDendro);
        if ($dend1 =~ /-/) {
            $dend1 =~ s/-//;
            $dend1 ="-".$proteinDendro{$dend1};
        }
        if ($dend2 =~ /-/) {
            $dend2 =~ s/-//;
            $dend2 = "-".$proteinDendro{$dend2};
        }
        push @proteinTab, join(",", $index, $dend1, $dend2, $height);
    }
    close PROTEIN_DENDRO;


    open (MATRIX,"$pathToFile/matrix.txt") || &stopOnError("ERROR!<BR>$!: '$pathToFile/matrix.txt'");
    my @quantifs;
    while (<MATRIX>) {
        chomp;
        if ($. == 1) {
            (my $empty,@quantifs) = split("\t",$_);
            #for (my $i = 1; $i <= $#quantif; $i++) {
            #    #$quantif[$i] =~ s/\"//g;
            #    push @quantifNoOrder, $quantif[$i];
            #}
            next;
        }
        my ($modProtID,@values) = split("\t",$_);
		my ($protID,$modCode)=$modProtID=~/^(\d+)-*(.*)/;
		if ($modCode && $modCode !~ /^\d+:/) {
			$modProtID=$protID.'-'.$oldSingleModID.':'.$modCode; # $oldSingleModID should be defined above
		}
        for (my $j=0; $j <= $#values; $j++) {
            if (!$values[$j] || $values[$j] eq 'NA') {$values[$j]='';} # NA imputed to 0 (before missMDA)
            $quantifProteinValues{$modProtID}{$quantifs[$j]} = $values[$j];
        }
    }
    close MATRIX;

    if (-e "$pathToFile/missingValues.txt") {
        my %quantifIndex;
        foreach my $i (0..$#quantifOrder) {$quantifIndex{$quantifOrder[$i]}=$i;}
        open (MISSING,"$pathToFile/missingValues.txt") || &stopOnError("ERROR!<BR>$!: '$pathToFile/missingValues.txt'");
        while (<MISSING>) {
            chomp;
            my ($modProtID,@quantifsNA)=split(/\t/,$_);
			my ($protID,$modCode)=$modProtID=~/^(\d+)-*(.*)/;
			if ($modCode && $modCode !~ /^\d+:/) {
				$modProtID=$protID.'-'.$oldSingleModID.':'.$modCode; # $oldSingleModID should be defined above
			}
            foreach my $quantif (@quantifsNA) {
                push @{$missingValues{$modProtID}},$quantifIndex{$quantif};
            }
        }
        close MISSING;
    }
}

my $strgProteinDendro = join(";", @proteinTab);
my $strgQuantifDendro = join(";", @quantifTab);
my $strgQuantif = join(",", @quantifInfo);

print qq
|<SCRIPT language="Javascript">
var HM;
var insertObject={};
var checkObject={};

function drawClustering() {
    HM = new heatMap({
        div:'heatMapDIV',
        editable:true,
        //colLabelOrientation:'diagonal',
        entities:{row:'$hmRowItem',column:'Quantification'},
        annotationOnDelete:deleteAnnotObject,
|;
if ($view eq "1D") {
    print qq
|       editable:false,
        colLabelPos:'bottom',
|;
}
else {
	if ($quantifFam eq 'RATIO') {
		print qq
|		editable:true, // scope & color
		normalization:{scope:'all',reference:'zero'},
|;
	}
 	else {print qq
|		editableItems:{scope:true,type:true,color:true},
		normalization:{scope:'all',reference:'min'},
|;
	}
   print qq
|		cellOnMouseOver:displayCellInfo,
		noMinCellHeight:true,
        labelsOnSelect:{row:{text:'List $hmRowItems',action:listProteins,ranking:true}},
        maxRange:{type:'quartile',value:7},
        flagText:'Imputed',
|;
}
print qq
|        exportAsImage:['Export as image','Clustering_$clusteringName','./exportSVG.cgi']
    });
    HM.setColumns([$strgQuantif]);
|;
if ($view eq '2D') {
    foreach my $modProtID (@proteinList) {
        my @values;
        foreach my $quantif (@quantifOrder) {
            push @values,$quantifProteinValues{$modProtID}{$quantif};
        }
        print "\tHM.addRow(['$proteinNames{$modProtID}','$modProtID','$proteinNames{$modProtID}'],[",join(',',@values),"]";
        print ",null,[",join(',',@{$missingValues{$modProtID}}),"]" if $missingValues{$modProtID}; # no excluded cell + list of index of flagged cells
        print ");\n";
    }
    print qq
|   HM.addTree('column','$strgQuantifDendro');\n
    HM.addTree('row','$strgProteinDendro');
|;
}
else {
    print qq|HM.addTree('column','$strgQuantifDendro');\n|;
}
print qq|
    HM.draw();
}
function annotateClustering() {
|;

###>Restoring saved annotations<###
my $sthSelectAnnot = $dbh -> prepare("SELECT NAME, RANK, ANNOT_TYPE, ANNOT_LIST FROM ANNOTATIONSET WHERE ID_EXPLORANALYSIS = $explorID ORDER BY RANK");
my $sthProperty=$dbh->prepare("SELECT NAME FROM PROPERTY WHERE ID_PROPERTY=?");
my $sthTheme=$dbh->prepare("SELECT NAME FROM CLASSIFICATION WHERE ID_CLASSIFICATION=?");
$sthSelectAnnot -> execute;
while (my ($name, $rank, $annotType, $annotList) = $sthSelectAnnot -> fetchrow_array) {
    if ($annotType =~ /^property:/) { # BioSample property/treatment
        my ($propID)=($name=~/^#(\d+)/); # extract property ID
        $sthProperty -> execute($propID);
        my ($propName)=$sthProperty -> fetchrow_array;
        next unless $propName; # to be safe in case property has been deleted
        print "\tajaxPropAnnotateClustering('$annotType:$propID','$propName',true);\n" ;
    }
    elsif ($annotType eq 'prot:THEME') { # Custom lists
        my ($themeID)=($name=~/^#(\d+)/); # extract classification ID
        $sthTheme -> execute($themeID);
        my ($themeName)=$sthTheme -> fetchrow_array;
        if ($view eq "1D") {
            print qq|
            annotRank.prot++;
            annotSetObject.prot['$themeName'] = {type:'prot:THEME',param:'#$themeID',position:annotRank.prot}; // # for ID tag
            |;
        }
        else {
            print "\tajaxThemeAnnotateClustering($themeID,'$themeName',true,'THEME');\n";
        }
    }
    elsif ($annotType eq 'prot:LIST') {
        $annotList=~s/#//g;
        if ($view eq "1D") {
            print qq|
            annotRank.prot++;
            annotSetObject.prot['$name'] = {type:'prot:LIST',param:'#$annotList',position:annotRank.prot}; // # for ID tag
            |;
        }
        else {
            print "\tajaxThemeAnnotateClustering('$annotList','$name',true,'LIST');\n";
        }

    }
    elsif ($annotType eq 'prot:GO:M') { # Multi-GO
        my ($goID,$goAspect,$termList)=split(':=:',$annotList);
        $goID=~s/^#//; # remove GOanaID tag
        if ($view eq "1D") {
            print qq|
            annotRank.prot++;
            annotSetObject.prot['$name'] = {type:'prot:GO:M',param:'$annotList',position:annotRank.prot}; // # for ID tag
            |;
        }
        else {
           print "\tajaxGoAnnotateClustering('$name',$goID,'$goAspect','$termList',true);\n";
        }
    }
    elsif ($annotType eq 'prot:PA:M') {
        my ($pathID, $termList)=split(':=:', $annotList);
        $pathID=~s/^#//;
        if ($view eq "1D") {
            print qq|
            annotRank.prot++;
            annotSetObject.prot['$name'] = {type:'prot:PA:M',param:'$annotList',position:annotRank.prot}; // # for ID tag
            |;
        }
        else {
           print "\tajaxGetPathwayProteinsList('$name','$termList',$pathID,true);\n";
        }
    }
}
$sthSelectAnnot -> finish;
$sthProperty -> finish;
$sthTheme -> finish;
$dbh -> disconnect;
print qq
|
}
</SCRIPT>|;
my $strgSelect = ($view eq "1D")? "Quantifications" : "<SELECT name=\"annotTarget\" class=\"title3\" onchange=\"selectAnnotTarget(this.value)\"><OPTION value=\"quantif\">Quantifications</OPTION><OPTION value=\"prot\">$hmRowItems</OPTION></SELECT>" ;
print qq|
<TABLE>
<TR>
    <TD rowspan=2><DIV id="heatMapDIV" style="float:left;background-color:#FFFFFF;border:solid 3px $darkColor;padding:2px;"></DIV></TD>
    <TD bgcolor="$darkColor" class="title3" style="height:30px" align="center">&nbsp;Annotation for:&nbsp;$strgSelect</TD>
</TR>
<TR><TD valign="top" nowrap><TABLE cellspacing=0 cellpadding=0 border=0>
    <TR><TD valign="top"><SELECT name="quantifAnnotationType" id="quantifAnnotationType" class="annotation title3" onchange="ajaxPropAnnotateClustering(this.value,this.options[this.selectedIndex].text,false)">
    <OPTION value="">-= Select =-</OPTION>
    <OPTGROUP label="BioSample properties:">
|;
if (scalar keys %{$bioSampProperties{'O'}}) {
    foreach my $propID (sort{lc($bioSampProperties{'O'}{$a}) cmp lc($bioSampProperties{'O'}{$b})} keys %{$bioSampProperties{'O'}}) {
        print "\t<OPTION value=\"property:O:$propID\">$bioSampProperties{O}{$propID}</OPTION>\n";
    }
}
else {print "\t<OPTION value=\"\" disabled>None</OPTION>\n";}
print qq
|	</OPTGROUP>
	<OPTGROUP label="BioSample treatments:">
|;
if (scalar keys %{$bioSampProperties{'T'}}) {
    foreach my $propID (sort{lc($bioSampProperties{'T'}{$a}) cmp lc($bioSampProperties{'T'}{$b})} keys %{$bioSampProperties{'T'}}) {
        print "\t<OPTION value=\"property:T:$propID\">$bioSampProperties{T}{$propID}</OPTION>\n";
    }
}
else {print "\t<OPTION value=\"\" disabled>None</OPTION>\n";}
print qq
|	</OPTGROUP>
    </SELECT></TD></TR>
    <TR><TD><SELECT id="protAnnotationType" class="annotation title3" onchange="updateProtAnnotationType(this)" style="display:none">
    <OPTION value="">-= Select =-</OPTION>
    <OPTGROUP label="GO Analyses:">
|;
if (scalar keys %goAnalyses) {
    foreach my $goID (sort{lc($goAnalyses{$a}[0]) cmp lc($goAnalyses{$b}[0])} keys %goAnalyses) {
        print "\t<OPTION value=\"GOA:$goID\">$goAnalyses{$goID}[0]</OPTION>\n";
    }
}
else {print "\t<OPTION value=\"\" disabled>** No GO Analyses found **</OPTION>\n";}
print qq
|   </OPTGROUP>
    <OPTGROUP label="Custom Lists:">
|;
if (scalar keys %listThemes) {
    foreach my $themeID (sort{lc($listThemes{$a}) cmp lc($listThemes{$b})} keys %listThemes) {
        print "\t<OPTION value=\"TH:$themeID\">$listThemes{$themeID}</OPTION>\n";
    }
}
else {print "\t<OPTION value=\"\" disabled>** No List found **</OPTION>\n";}
print qq
|   </OPTGROUP>
    <OPTGROUP label="Pathway Analyses:">
|;
if (scalar keys %pathwayAnalyses) {
    foreach my $pathID (sort{lc($pathwayAnalyses{$a}[1]) cmp lc($pathwayAnalyses{$b}[1])} keys %pathwayAnalyses) {
        print "<OPTION value=\"PA:$pathID\">$pathwayAnalyses{$pathID}[1]</OPTION>\n";
    }
}
else {print "<OPTION value=\"\" disabled>** No Pathway annotations **</OPTION>\n";}
print qq
|</OPTGROUP>
    </SELECT></TD></TR>
    <TR><TD>
    <SELECT id="goAspect" class="annotation title3" style="display:none" onchange="ajaxUpdateGoTermList(this.value)"><OPTION value="">-= Select Aspect =-</OPTION></SELECT></TD></TR>
    <TR><TD>
        <DIV id="termDisplayDIV" style="display:none">
        <DIV id="selTermsDIV"></DIV>
        <B>Annotation name:</B><INPUT type="text" id="annotName" style="width:160px"/><DIV id="annotTermsDIV" style="display:inline-block;"></DIV>
        <!--<INPUT type="button" value="Annotate" onclick="userGoAnnotateClustering()"/>-->
        </DIV>
    </TD></TR>
    <TR><TD><DIV id="listDisplayDIV" style="max-height:150px;max-width:300px;overflow:auto;display:none"></DIV><DIV id="saveList" style="display:none"><B>Annotation name:</B><INPUT type="text" id="annotListName" style="width:160px"/><br><INPUT type="button" value="Annotate" id="annotList" onclick="userListAnnotateClustering()"/></DIV></TD></TR>
    <TR><TD><INPUT type="button" id="saveButton" style="font-weight:bold;display:none" onclick="ajaxSaveAnnotations()" value="Save annotations"$disabSave></TD></TR>
    <TR><TD><SPAN id="saveSPAN" class="title3" style="color:#FFF;background-color:#387D38;padding:1px 7px;display:none">Annotations saved !</SPAN></TD></TR>
|;

my $distribPlotPath = "$promsPath{explorAna_html}/project_$projectID/$explorID/";
print qq |
    <TR><TH><br/>
        <button id='valueDistributionButton' onclick='toggleValueDistribution();' />Display values distribution</button><br/>
		<img src='$distribPlotPath/valueDistribution.png' id='valueDistribution' style='display:none' />
    </TH></TR>
| if(-e "$pathToFile/valueDistribution.png");


print qq |
    </TABLE></TD></TR>
</TABLE>
<DIV id="protListDIV" style="clear:both;height:535;overflow:auto"></DIV>
<DIV id="displayDIV" class="popup"> <!--filter:alpha(opacity=80);opacity:0.8;-->
	<DIV id="infoDIV"></DIV>
	<DIV id="infoDIV2"></DIV>
</DIV>

<DIV id="divDescription" class="clDescriptionCont"></DIV>
<SCRIPT type="text/javascript">
setPopup();

/********* DRAWING CLUSTERING *********/
drawClustering();
document.getElementById('waitDIV').style.display='none';
document.getElementById('resultDIV').style.visibility='visible';

/********* ANNOTATING CLUSTERING *********/
annotateClustering();

</SCRIPT>
</BODY>
</HTML>
|;


sub ajaxPropAnnotateClustering {
    print header(-type=>'text/plain',-charset => 'utf-8');
    warningsToBrowser(1);

    my ($propType,$propID)=(param('property')=~/property:(.):(\d+)/);

    my $sthPV=$dbh->prepare("SELECT DISTINCT BP.PROPERTY_VALUE FROM BIOSAMPLE_PROPERTY BP
                                INNER JOIN OBSERVATION O ON BP.ID_BIOSAMPLE=O.ID_BIOSAMPLE
                                INNER JOIN ANA_QUANTIFICATION AQ ON O.ID_ANALYSIS=AQ.ID_ANALYSIS
                                INNER JOIN EXPLORANA_QUANTIF EQ ON AQ.ID_QUANTIFICATION=EQ.ID_QUANTIFICATION
                                WHERE BP.ID_PROPERTY=$propID AND EQ.ID_EXPLORANALYSIS=$explorID");
    my ($propName)=$dbh->selectrow_array("SELECT NAME FROM PROPERTY WHERE ID_PROPERTY=$propID");
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

    my $jsAnnotStrg='';
    foreach my $propValue (sort{&promsMod::sortSmart(lc($a),lc($b))} keys %propValueCode) {
        my %matchedQuantifs=&promsQuantif::getQuantificationsFromPropertyValue($dbh,$propID,$propType,$propValue,param('qList'));
        if (scalar keys %matchedQuantifs) {
            $jsAnnotStrg.="," if $jsAnnotStrg;
            $jsAnnotStrg.="'$propValueCode{$propValue}':'".join(',',keys %matchedQuantifs)."'";
        }
    }
    print "HM.addAnnotation('column','$propName',{$jsAnnotStrg},'^###(\$|%)');";
    $dbh->disconnect;
    exit;
}
sub ajaxThemeAnnotateClustering {
    print header(-type=>'text/plain',-charset => 'utf-8');
    warningsToBrowser(1);

    my $themeID=param('themeID');
    my $themeName=param('name');
    my $type=param('type');

    my %protInCluster;
    &getProteinsInCluster(\%protInCluster);
    #my ($themeName)=$dbh->selectrow_array("SELECT NAME FROM CLASSIFICATION WHERE ID_CLASSIFICATION=$themeID");
    my $sthList=$dbh->prepare("SELECT ID_CATEGORY,NAME,LIST_TYPE FROM CATEGORY WHERE ID_CLASSIFICATION=? ORDER BY DISPLAY_POS");
	#my $sthProt=$dbh->prepare("SELECT DISTINCT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=?");
    $sthList->execute;
    my $jsAnnotStrg='';
    if ($type eq 'THEME') {
        $sthList->execute($themeID);
        while (my ($listID,$listName,$type)=$sthList->fetchrow_array) {
			if ($isPtmQuantif && $type eq 'SITE') {
				my %siteList;
				&promsQuantif::fetchCustomList($dbh,$listID,\%siteList);
				my $first=1;
				foreach my $modProtID (keys %siteList) {
					next unless $protInCluster{$modProtID};
					if ($first) {
						$jsAnnotStrg.=',' if $jsAnnotStrg;
						$jsAnnotStrg.="'$listName':'$modProtID";
						$first=0;
					}
					else {$jsAnnotStrg.=",$modProtID";}
				}
				$jsAnnotStrg.="'" unless $first;
			}
			else {
				#$sthProt->execute($listID);
				my %protList;
				&promsQuantif::fetchCustomList($dbh,$listID,\%protList,1);
				my $first=1;
				#while (my ($protID)=$sthProt->fetchrow_array) { #}
				foreach my $protID (keys %protList) {
					next unless $protInCluster{$protID};
					if ($first) {
						$jsAnnotStrg.=',' if $jsAnnotStrg;
						$jsAnnotStrg.="'$listName':'$protID";
						$first=0;
					}
					else {$jsAnnotStrg.=",$protID";}
				}
				$jsAnnotStrg.="'" unless $first;
			}
        }
    }
    else {
        foreach my $listID (split(':',$themeID)) {
            #$listID=~s/#//g;
            my ($listName)=$dbh->selectrow_array("SELECT NAME FROM CATEGORY WHERE ID_CATEGORY=$listID");
            #$sthProt->execute($listID);
            my $first=1;
			my %protList;
			&promsQuantif::fetchCustomList($dbh,$listID,\%protList,1);
            #while (my ($protID)=$sthProt->fetchrow_array) { #}
			foreach my $protID (keys %protList) {
                next unless $protInCluster{$protID};
                if ($first) {
                    $jsAnnotStrg.=',' if $jsAnnotStrg;
                    $jsAnnotStrg.="'$listName':'$protID";
                    $first=0;
                }
                else {$jsAnnotStrg.=",$protID";}
            }
            $jsAnnotStrg.="'" unless $first;
        }
    }
    $sthList->finish;
    #$sthProt->finish;
    $dbh->disconnect;
    print "HM.addAnnotation('row','$themeName',{$jsAnnotStrg}$protAnnotMatchStrg);"; # compatible with modif quantif
    exit;
}

sub ajaxGoAnnotateClustering {
    print header(-type=>'text/plain',-charset => 'utf-8');
    warningsToBrowser(1);

    my $annotName=param('annotName');
    my $goID=param('goID');
    my $aspect=param('goAspect');
    my (@termList,%termName);
    foreach my $termInfo (split(':@:',param('termsStrg'))) {
        my ($termID,$termText)=split('@',$termInfo);
        push @termList,$termID;
        $termName{$termID}=$termText;
    }
    my $numTerms=scalar @termList;
    my %protInCluster;
    &getProteinsInCluster(\%protInCluster);
    my $resultFile="$promsPath{go_unix}/project_$projectID/$goID/results_$aspect.txt";

    my (%matchedTerms,%matchedProt);
    my $numMatches=0;
    open (TERMS,$resultFile);
    TERM:while (<TERMS>) { # terms are sorted by ascending p-value
        foreach my $termID (@termList) {
            next if $matchedTerms{$termID};
            if (/^$termID\t/) {
                my @protList=split(';',(split(/\t/),$_)[5]);
                foreach my $protID (@protList) {
                    next unless $protInCluster{$protID}; # restrict to proteins used for clustering
                    next if $matchedProt{$protID}; # do not reuse a prot already matched to a previous (better) term
                    push @{$matchedTerms{$termID}},$protID;
                    $matchedProt{$protID}=1;
                }
                $numMatches++;
                last TERM if $numMatches==$numTerms;
            }
        }
    }
    close TERMS;
    if (scalar keys %matchedTerms) {
        print "HM.addAnnotation('row','$annotName',{";
        my $first=1;
        foreach my $termID (sort{&promsMod::sortSmart(lc($a),lc($b))} keys %matchedTerms) {
            if ($first) {$first=0;}
            else {print ',';}
            print "'$termName{$termID}':'",join(',',@{$matchedTerms{$termID}}),"'";
        }
        print "}$protAnnotMatchStrg);";
    }
    else {
        print "(function(){alert('No proteins matched!'); return false;})()";
    }
    exit;
}

sub ajaxPathwayAnnotateClustering {

    print header(-type=>'text/plain',-charset => 'utf-8');
    warningsToBrowser(1);

    my $annotName=param('annotName');
    my $termsStrg=param('termsStrg');
    my $pathwayID=param('pathID');
    my $projectPath="$promsPath{pathAna}/project_$projectID";
    my $pathwayPath=$projectPath."/".$pathwayID;
    my ($catID)=$dbh->selectrow_array("SELECT ID_CATEGORY from PATHWAY_ANALYSIS where ID_PATHWAY_ANALYSIS=$pathwayID");
    my %protInCluster;
    &getProteinsInCluster(\%protInCluster);
    my %reactItem;
    foreach my $item (split(":\@:",$termsStrg)) {
        foreach my $reactNum (split("\@", $item)) {
            $reactItem{$reactNum}=$reactNum;
        }
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
    my %protIDs;
    while(my ($protID, $uniprot) = $sthProt->fetchrow_array){
        $protIDs{$uniprot}{$protID}=1;
    }
    $sthProt->finish;
    my @uniprotID;
    my (%submitFoundEntities, %pathInfo);
    my $csv = Text::CSV->new({ sep_char => ',' });
    open(PATHWAY,"$pathwayPath/pathway.csv" );
        while (my $line=<PATHWAY>) {
            chomp($line);
            next if ($. == 1);
            if ($csv->parse($line)) {
                my @fields=$csv->fields();
                my %existUniprot;
                foreach my $item (@fields) {
                    if ($reactItem{$fields[0]}) {
                        $pathInfo{$reactItem{$fields[0]}} = $fields[1];
                        foreach my $uniprot (split(";",$fields[12])) {
                            next if ($existUniprot{$uniprot});
                            push @{$submitFoundEntities{$reactItem{$fields[0]}}}, $uniprot;
                            $existUniprot{$uniprot}=1;
                        }
                        next;
                    }
                }
            }
        }
    close(PATHWAY);

    my (%protInfo, %protGene);
    foreach my $reactNumber (keys %submitFoundEntities) {
        foreach my $uniprotID (@{$submitFoundEntities{$reactNumber}}) {
            if ($protIDs{$uniprotID}) {
                foreach my $protID (keys %{$protIDs{$uniprotID}}) {
                    next unless $protInCluster{$protID}; # restrict to proteins used for clustering
                    push @{$protInfo{$pathInfo{$reactNumber}}},$protID;
                }
            }
        }
    }
    if(scalar(keys %protInfo)){
        print "HM.addAnnotation('row','$annotName',{";
        my $first=1;
        foreach my $pathName (sort{&promsMod::sortSmart(lc($a),lc($b))} keys %protInfo) {
            if ($first) {$first=0;}
            else {print ',';}
            print "'$pathName':'",join(",", @{$protInfo{$pathName}}),"'";
        }
        print "}$protAnnotMatchStrg);";
    }
    else {
        print "(function(){alert('No proteins matched!'); return false;})()";
    }
    $dbh-> disconnect;
    exit;
}

sub ajaxSaveAnnotations {
    print header(-'content-encoding' => 'no',-charset => 'utf-8');
    warningsToBrowser(1);

    ##<Deleting all previous annot for Clustering
    $dbh -> do("DELETE FROM ANNOTATIONSET WHERE ID_EXPLORANALYSIS = $explorID");

    my ($annotSetID) = $dbh -> selectrow_array("SELECT MAX(ID_ANNOTATIONSET) FROM ANNOTATIONSET");
    my $sthInsertAnnot = $dbh -> prepare("INSERT INTO ANNOTATIONSET(ID_ANNOTATIONSET,ID_EXPLORANALYSIS,NAME,RANK,ANNOT_TYPE,ANNOT_LIST) values (?,$explorID,?,?,?,?)");

    if (param('string')) {
        my @annotList;
        foreach my $annotStrg (split(":%:",param('string'))) {
            print "annotStrg=$annotStrg\n";
            my %params;
            foreach my $dataStrg (split('//',$annotStrg)) {
                print "dataStrg=$dataStrg\n";
                my @data=split('==',$dataStrg);
                $params{$data[0]}=$data[1];
            }
            push @annotList,\%params;
        }

        foreach my $refAnnot (@annotList) {
            $annotSetID++;

            ##<Protein Annot
            #<GO annotation #<PATHWAY annotation
            if ($refAnnot->{type} eq 'prot:GO:M' || $refAnnot->{type} eq 'prot:PA:M') { # multiple GO terms in same annot
                $sthInsertAnnot -> execute($annotSetID,$refAnnot->{name},$refAnnot->{position},$refAnnot->{type},$refAnnot->{param});
            }
            ##<Theme annotation
            elsif ($refAnnot->{type} eq 'prot:THEME') {
                $sthInsertAnnot -> execute($annotSetID,$refAnnot->{param},$refAnnot->{position},$refAnnot->{type},undef); # '#themeID' used as name
            }
            #<List annotation
            elsif ($refAnnot->{type} eq 'prot:LIST') {
                $refAnnot->{param}=~s/#//g;
                my @listStrg=();
                foreach my $listID (split(':',$refAnnot->{param})) {
                    push @listStrg, "#$listID";
                }
                my $joinList = join(':',@listStrg);
                $sthInsertAnnot -> execute($annotSetID,$refAnnot->{name},$refAnnot->{position},$refAnnot->{type},$joinList);
            }
            ###<Quantif Annot
            ##<Property/treatment
            elsif ($refAnnot->{type} =~ /property:/) {
                $sthInsertAnnot -> execute($annotSetID,$refAnnot->{param},$refAnnot->{position},$refAnnot->{type},undef); # '#propID' used as name
            }
        }
    }
    $sthInsertAnnot -> finish;
    $dbh -> commit;
    $dbh -> disconnect;

    print "###OK###";
    exit;
}

sub getProteinsInCluster { # records both protID & (modProtID if any); GLOBALS: $explorID,$dbh,$pathToFile,$oldSingleModID
    my ($refProtList)=@_;
	open (PROTEIN_ORDER,"$pathToFile/protCluster.txt") || die "Error: $!";
    while (<PROTEIN_ORDER>) {
        chomp;
        next if $. == 1;
        my ($index, $pos, $modProtID) = split("\t",$_);
		my ($protID,$modCode)=$modProtID=~/^(\d+)-*(.*)/;
		$refProtList->{$protID}=1;
		if ($modCode) {
			if ($modCode !~ /^\d+:/) { # no starting modifID => old format
				$modProtID=$protID.'-'.$oldSingleModID.':'.$modCode;
			}
		}
        $refProtList->{$modProtID}=1;
    }
    close PROTEIN_ORDER;
}

sub ajaxGetListTheme {
    print header(-type=>'text/plain',-charset => 'utf-8');
    warningsToBrowser(1);
    #my ($lightColor,$darkColor)=&promsConfig::getRowColors;

    my $themeID=param('themeID');
	my %listTypeName=('PROT'=>'Proteins','SITE'=>'Sites');
    my ($themeName)=$dbh->selectrow_array("SELECT NAME FROM CLASSIFICATION WHERE ID_CLASSIFICATION=$themeID");
    my $sthList=$dbh->prepare("SELECT ID_CATEGORY,NAME,LIST_TYPE FROM CATEGORY WHERE ID_CLASSIFICATION=$themeID ORDER BY DISPLAY_POS");
    $sthList->execute();
    #print qq|<SELECT MULTIPLE>|;
    #my $bgColor=$lightColor;
    print qq|<INPUT type="checkbox" value="TH:$themeID:$themeName" name="list" onclick="selectBoxes('TH')">All lists in theme<br>|;
    while (my ($listID,$listName,$listType)=$sthList->fetchrow_array) {
		$listType='PROT' unless $listType;
        print qq|<INPUT type="checkbox" value="LI:$listID:$listName" name="list" onclick="selectBoxes('LI')">$listName [$listTypeName{$listType}]<br>|;
        #print qq|<OPTION value="$listID">$listName</OPTION>\n|;
        #$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
    }

    #print qq|</SELECT>|;
    $sthList->finish;
    $dbh->disconnect;
    exit;
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

####>Revision history<####
# 1.2.6 [ENHANCEMENT] Better management of quantification parameters for ajaxListProt (PP 07/01/20)
# 1.2.5 [BUGFIX] in SQL query requiring GROUP_CONCAT for TARGET_POS & [UX] Smooth scroll to list of proteins (PP 20/11/19)
# 1.2.4 [ENHANCEMENT] Multi-site support (PP 24/10/19)
# 1.2.3 [ENHANCEMENT] Integrates R generated density plot (count normalized) as quality representation of imputed quantification values (VS 14/10/19)
# 1.2.2 Minor improvement in highlighting annotations sorting (PP 20/06/18)
# 1.2.1 Compatible with modification sites list (PP 12/11/17)
# 1.2.0 Compatible with non-ratio (eg. MaxQuant) quantif (PP 11/01/17)
# 1.1.0 display quantifName for clustering with normalization (SL 19/07/16)
# 1.0.9 Add missing JS functions and DIV for displaying detailed protein quantification data (PP 13/07/16)
# 1.0.8 Minor change to allow 1px cell height on 2D clustering (PP 28/05/16)
# 1.0.7 Missing values flag and data max. range for 2D view (PP 29/01/16)
# 1.0.6 select multi list from theme for annotation (SL 11/06/15)
# 1.0.5 Compatibility between phospho-sites and custom list annotation (PP 22/05/15)
# 1.0.4 add pathway annotation and correct bin from GO analysis and add view 1D/2D (SL 24/11/14)
# 1.0.3 encodeURIComponent on GO annotName for AJAX (PP 28/08/14)
# 1.0.2 JS ajaxManageSaveProteins moved to promsMod.pm (PP 22/08/14)
# 1.0.1 Quantification + protein annotation & synchronized with new data file naming (PP 08/08/14)
# 1.0.0 New script to display and store the results of PCA (SL 08/04/14)
