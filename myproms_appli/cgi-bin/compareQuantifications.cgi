#!/usr/local/bin/perl -w

################################################################################
# compareQuantifications.cgi          1.6.5                                    #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Simple display of multiple quantifications                                   #
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
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use promsQuantif;
use strict;
use Spreadsheet::WriteExcel;
use utf8; # Tells Perl that characters are UTF8. Necessary for Excel export to work with UTF-8 characters ...
use Encode qw(encode_utf8); # ... Encode needed if hard UTF8 chars in code!!!
use File::Path qw(rmtree); # remove_tree

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my $userID=$ENV{'REMOTE_USER'};
my ($MIN_INF_RATIO,$MAX_INF_RATIO,$MAX_INF_RATIO_DB,$MIN_INF_RATIO_DB)=&promsQuantif::getExtremeRatios; # Max ratios allowed before switching to infinite
my $infChar=encode_utf8('∞');
my %proteinQuantifFamilies=&promsQuantif::getProteinQuantifFamilies;
my %maxQuantMeasures; #=('MQ_INT'=>'Intensity','MQ_IBAQ'=>'iBAQ','MQ_LFQ'=>'LFQ','MQ_SC'=>'MS/MS count');
foreach my $refMeas (@{$proteinQuantifFamilies{'MEASURES'}{'MQ'}}) {$maxQuantMeasures{$refMeas->[0]}=$refMeas->[1];}
my %emPAIMeasures; #=('EMPAI'=>'emPAI','EMPAI_MOL'=>'emPAI Mol. %','EMPAI_MR'=>'emPAI Mr. %');
foreach my $refMeas (@{$proteinQuantifFamilies{'MEASURES'}{'EMPAI'}}) {$emPAIMeasures{$refMeas->[0]}=$refMeas->[1];}
my %sinMeasures; #=('SIN_SIN'=>'SI<sub>N</sub>');
foreach my $refMeas (@{$proteinQuantifFamilies{'MEASURES'}{'SIN'}}) {$sinMeasures{$refMeas->[0]}=$refMeas->[1];}

####################
####>Parameters<####
####################
my $action=(param('ACT'))? param('ACT') : (param('AJAX'))? 'results' : 'frames';
my $selQuantifFamily=param('quantifFamily') || '';
my $selQuantifFocus=param('quantifFocus') || 0; # can be  0 (prot quantif)|modifID (modif quantifs only)|-modifID (modif & prot quantifs)
my $proteoform=($selQuantifFocus > 0)? 'site' : ($selQuantifFocus == 0)? 'protein' : 'protein/site';
my $Proteoform=ucfirst($proteoform);
my $selModifID=abs($selQuantifFocus);
my $projectID=param('id_project'); # not always defined
my $parentItem=(param('parentItem'))? lc(param('parentItem')) : ''; # not defined for AJAX calls
my ($item,$itemID)=split(/:/,$parentItem);
my $view=param('view') || 'list'; $view='list' if ($view eq 'correlMatrix' && $selQuantifFocus < 0); # not for mixte proteo/modif-proteo
my $sortOrder=param('sort') || 'protein';
my $restrictListID=param('restrictList') || 0; # exclusion/selection List ID
my $listAction=param('listAction') || 'restrict';
my $quantifList=param('quantifList');
my @selectedQuantifications=($quantifList)? split(':',$quantifList) : (); #param('usedQuantifs'); usedQuantifs requires selection of all options
my $numSelectedQuantifs=scalar @selectedQuantifications;
my $dispFoldChange=(param('foldChange'))? param('foldChange') : ($view eq 'volcano')? 2 : 1; $dispFoldChange=1 if $dispFoldChange < 1;
my $invDispFoldChange=1/$dispFoldChange;
my $dispPvalue=param('pValue') || 1; $dispPvalue=1 if ($dispPvalue < 0 || $dispPvalue > 1);  # 1 => will take quantif p-value as default
my $dispPepType=(param('pepType'))? param('pepType') : ($selQuantifFamily eq 'MQ')? 'RAZ_UNI_PEP' : 'all';
my $numPepCode=($dispPepType eq 'all')? 'NUM_PEP_USED' : ($dispPepType eq 'distinct')? 'DIST_PEP_USED' : $dispPepType;
my $dispNumPep=param('numPep') || 1; $dispNumPep=1 if $dispNumPep < 1;
my $empaiValue=param('empaiType') || 'EMPAI';
my $dispMeasure=param('dispMeasure') || 'MQ_INT';
my $strictFilter=param('strictFilter') || 0;
my $ajax=0;

if ($action eq 'frames') {&generateFrames; exit;}

####>Connecting to the database<####
my $dbh=&promsConfig::dbConnect;
unless ($projectID) {
	if ($parentItem) {
		$projectID=&promsMod::getProjectID($dbh,$itemID,$item);
	}
	elsif (param('branchID')) {
		my ($parentItem,$parentID)=split(':',param('branchID'));
		$projectID=&promsMod::getProjectID($dbh,$parentID,$parentItem);
	}
	else {
		$projectID=&promsMod::getProjectID($dbh,(split('_',$selectedQuantifications[0]))[0],'QUANTIFICATION');
	}
}
if (!$selQuantifFamily && $selectedQuantifications[0]) {
	my ($quantifID)=(split('_',$selectedQuantifications[0]))[0];
	my ($qMethCode)=$dbh->selectrow_array("SELECT QM.CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD QM WHERE Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$quantifID");
	$selQuantifFamily=($qMethCode=~/PROT_RATIO_PEP|TNPQ/)? 'RATIO' : $qMethCode;
}

my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
my $exportType=param('exportType') || '';
if (param('AJAX')) {
	$ajax=1;
	if (param('AJAX') eq 'getquantifList') {&ajaxGetQuantificationList; exit;}
	elsif (param('AJAX') eq 'ajaxListProt') {&ajaxListSelectedProteins; exit;} # also called by displayClustering.cgi & displayPCA.cgi
	exit;
}
elsif ($action eq 'selQuantifs') {&selectQuantifications; exit;}


###############
#     Main    #
###############
#print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
my ($workbook,%itemFormat); # globals for export to Excel
my (%selectedProteins,%excludedProteins); # proteins selected/excluded
if ($exportType) {
	#################################
	####>Prepare export to Excel<####
	#################################
	if ($exportType eq 'selected') {
		foreach my $modProtID (param('chkProt')) {
			$selectedProteins{$modProtID}=1;
		}
	}

	my $timeStamp1=strftime("%Y%m%d %H-%M",localtime);
	#my $timeStamp2=strftime("%Y-%m-%d %H:%M",localtime);

	$workbook=Spreadsheet::WriteExcel->new("-");
	$workbook->set_properties(title=>"Quantifications comparison",
							  author=>'myProMS server',
							  comments=>'Automatically generated with Perl and Spreadsheet::WriteExcel'
							  );
	$workbook->set_custom_color(40,189,215,255); # light light color  #BDD7FF (new myProMS colors V3.5+)
	#$workbook->set_custom_color(41,128,179,255); # dark blue color  #80B3FF
	%itemFormat=(
			title =>			$workbook->add_format(align=>'center',size=>18,bold=>1,border=>1),
			header =>			$workbook->add_format(align=>'center',size=>10,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			headerR =>			$workbook->add_format(align=>'right',valign=>'top',size=>10,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			mergeRowHeader =>	$workbook->add_format(align=>'center',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			mergeColHeader =>	$workbook->add_format(align=>'center',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			mergeColHeaderL =>	$workbook->add_format(align=>'left',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			text =>				$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>0,border=>1),
			textWrap =>			$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>1,border=>1),
			#mergeText =>		$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>0,border=>1),
			mergeColText =>		$workbook->add_format(align=>'left',size=>10,text_wrap=>1,border=>1),
			#numberVis =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,bold=>1),
			#numberHid =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1),
			#numberVisBC =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,bold=>1,color=>'grey'),
			#numberHidBC =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,color=>'grey'),
			#numberVisVP =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,bold=>1,italic=>1),
			#numberHidVP =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,italic=>1),
			number =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1),
			#mergeNumber =>		$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1),
			number1d =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',num_format=>'0.0',border=>1),
#mergeNumber1d =>	$workbook->add_format(size=>10,align=>'center',valign=>'top',num_format=>'0.0',border=>1)
			);

	print header(-type=>"application/vnd.ms-excel",-attachment=>"compare_quantifications_$timeStamp1.xls");
}
else {
	#######################
	####>Starting HTML<####
	#######################
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1); # writes PERL warnings as HTML comments <!-- --> must be activated after header
	my $viewScript=($view eq 'log2')? 'genericPlot.js' : 'volcanoPlot2.js'; # DO NOT LOAD BOTH: function conflicts!!!
	my $fontSize=($numSelectedQuantifs>=40)? '8px' : ($numSelectedQuantifs>=25)? '10px' : ($numSelectedQuantifs>=15)? '11px' : '12px';
	print qq
|<HEAD>
<TITLE>Compare Quantifications</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.TD {font-weight:normal;}
.TH {font-weight:bold;}
.LINK {cursor:pointer;}
.popup {z-index:100;background-color:#FFFFFF;border:solid 3px #999999;padding:5px;position:absolute;display:none;box-shadow:10px 10px 20px #808080;}
.row_0{background-color:$darkColor;}
.row_1{background-color:$lightColor;}
.rowHeader{background-color:$darkColor; min-height:50px; text-align:right;}
.colHeader{background-color:$darkColor; min-width:75px;}
TABLE.correlation TH{font-size:$fontSize;}
</STYLE>
<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/$viewScript"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/vennDiagram.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/peptidePlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/heatMap.js"></SCRIPT> <!-- Needed by ajaxPepView for DIA -->
<SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
<SCRIPT type="text/javascript">
|;
	&promsMod::popupInfo();
	print qq
|function sequenceView(id_protein,id_analyses){
	var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+id_analyses+"&id_prot="+id_protein+"&msdata="+top.showMSdata;
	top.openProtWindow(winLocation);
}
var spectWin,selectedPep;
function drawSpectrum(pepId,pepInfo) {
	if (selectedPep && document.getElementById(selectedPep)) { // selectedPep could have been remove by another AJAX call
		document.getElementById(selectedPep).style.color='#000000';
	}
	if (pepId) { // not defined if called from SVG SWATH heatmap
		document.getElementById(pepId).style.color='#DD0000'; //style.background = '#DD0000';
	}
	selectedPep=pepId;
	var paramString="RID="+pepInfo+"&CALL=pep";
	spectWin=window.open("$promsPath{cgi}/drawSpectrum.cgi?"+paramString,'SpectrumWindow','width=950,height=950,location=no,resizable=yes,scrollbars=yes');
	spectWin.focus();
}
function exportProteins(expType) { // only for list view
	if (!expType) {expType='all';}
	myForm=document.protForm;
	myForm.exportType.value=expType;
	if (expType=='selected') {
		var okProt=false;
		if (myForm.chkProt.length) {
			for (var i=0; i<myForm.chkProt.length; i++) {
				if (myForm.chkProt[i].checked) {okProt=true; break;}
			}
		}
		else if (myForm.chkProt.checked) {okProt=true;}
		if (!okProt) {
			alert('ERROR: No proteins selected!');
			return false;
		}
	}
	myForm.submit();
}
function selectSort(newSort,ajax) { // view=list (true or from Volcano ajax)
	//parent.selQuantifsFrame.setComparisonAsModified(2); // 2! because stored in top.promsFrame.selectedView
	top.promsFrame.selectedView=newSort; //openProject.cgi
	var myForm=parent.selQuantifsFrame.document.compForm;
	myForm.sort.value=newSort;
	var checkedProt=document.protForm.chkProt;
	var chkList=[];
	if (ajax) { // Take all proteins & recall ajaxListSelectedProteins
		if (checkedProt.length) {
			for (var i=0; i<checkedProt.length; i++) {chkList.push(checkedProt[i].value);}
		}
		else {chkList.push(checkedProt[i].value);}
		ajaxListSelectedProteins(chkList.join(','),null,newSort);
	}
	else { // Take selected & resumit comparison
		if (checkedProt.length) {
			for (var i=0; i<checkedProt.length; i++) {
				if (checkedProt[i].checked) {chkList.push(checkedProt[i].value);}
			}
		}
		else if (checkedProt.checked) {chkList.push(checkedProt[i].value);}
		myForm.chkProtList.value=chkList.join(':');
		myForm.submit();
	}
}
// Color list for plot
var hColors=['#E18B6B','#95B9C7','#7E2217','#9A9A9A','#8AFB17','#FBB917','#F660AB'];
var hColorIdx=0;
function getElementPosition(e) {
	var left=0;
	var top=0;
	while (e.offsetParent != undefined && e.offsetParent != null) { //Adding parent item position
		left += e.offsetLeft + (e.clientLeft != null ? e.clientLeft : 0);
		top += e.offsetTop + (e.clientTop != null ? e.clientTop : 0);
		e = e.offsetParent;
	}
	return [top,left];
}
function checkAllProteins(chkStatus) {
	var checkBoxList=document.protForm.chkProt;
	if (checkBoxList.length) {
		for (let i=0; i < checkBoxList.length; i++) {
			if (document.getElementById('tr_'+checkBoxList[i].value).style.display=='') checkBoxList[i].checked=chkStatus; // compatibility with selection from Venn Diagram
		}
	}
	else {checkBoxList.checked=chkStatus;}
}
// AJAX --->
function ajaxListSelectedProteins(selectedPoints,thresholds,sort) { // thresholds is not used
	var sortOrder=(sort)? sort : '$sortOrder'; // called from Volcano,heatmap (no sort) or Ajax list for resorting (sort)
	saveListGlobals.themesFetched=false; // used for themes management
	document.getElementById('displayDIV').style.display='none'; // in case popup is displayed
	// Adjusting protListDIV position to graph (graph can move due to summary option display)
	var graphDiv=document.getElementById('mainGraphDIV');
	var divX = getDomOffset(graphDiv,'offsetLeft') + graphDiv.offsetWidth + 2;
	var divY = getDomOffset(graphDiv,'offsetTop');
	var listDiv=document.getElementById('protListDIV');
	listDiv.style.left = divX + 'px';
	listDiv.style.top = divY + 'px';
	// Updating contents
	listDiv.innerHTML="<BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<FONT class=\\"title3\\">Fetching data. Please wait...</FONT><BR>&nbsp;&nbsp;&nbsp;&nbsp;<IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR>";
	listDiv.style.display='block';
	//Parameters 1
	var qSetList=[];
	if ('$view'=='log2') { // quantif are paired (2 dim)
		var usedQuantif={};
		for (var i=0; i<dataSetQuantif.length; i++) {
			if (!usedQuantif[dataSetQuantif[i][0]]) {
				qSetList.push(dataSetQuantif[i][0]);
				usedQuantif[dataSetQuantif[i][0]]=1;
			}
			if (!usedQuantif[dataSetQuantif[i][1]]) {
				qSetList.push(dataSetQuantif[i][1]);
				usedQuantif[dataSetQuantif[i][1]]=1;
			}
		}
	}
	else { // 1 dim
		qSetList=dataSetQuantif;
	}
	var paramStrg="AJAX=ajaxListProt&quantifFamily=$selQuantifFamily&quantifFocus=$selQuantifFocus&quantifList="+qSetList.join(':')+"&view=list&sort="+sortOrder+"&pepType=$dispPepType&numPep=$dispNumPep";
	if ('$selQuantifFamily'=='MQ') {
		paramStrg+="&dispMeasure=$dispMeasure";
	}
	else if ('$selQuantifFamily'=='EMPAI') {
		paramStrg+="&empaiType=$empaiValue";
	}
	else if ('$selQuantifFamily'=='RATIO') {
		paramStrg+="&foldChange=$dispFoldChange&pValue=$dispPvalue";
	}
	paramStrg+="&selProt=";

	//Parameters 2
	if (sort) { // Called from Ajax List => add string
		paramStrg+=selectedPoints; // string
	}
	else { // Called from Volcano => extracting no-duplicates list of proteins)
		var noDupProteins=new Object();
		for (var gr in selectedPoints) { // Object
			for (var i=0; i<selectedPoints[gr].length; i++) {
				noDupProteins[selectedPoints[gr][i]]=1;
			}
		}
		var p1=true;
		for (var prot in noDupProteins) {
			if (!p1) {paramStrg+=',';}
			paramStrg+=prot;
			p1=false;
		}
	}
	//paramStrg+='&thresholds='+thresholds.join(',');
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("POST","$promsPath{cgi}/compareQuantifications.cgi",true);
	//Send the proper header information along with the request
	XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
	XHR.setRequestHeader("Content-length", paramStrg.length);
	XHR.setRequestHeader("Connection", "close");
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			listDiv.innerHTML=XHR.responseText;
		}
	}
	XHR.send(paramStrg);
}
/*
function ajaxGoDecorateGraph(termIdStrg) {
	if (!termIdStrg) return;
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/compareQuantifications.cgi?ACT=ajaxTermProt&goStrg="+termIdStrg,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			var termData=termIdStrg.split(',');
			addHighlighting(VP,termData[2],hColors[hColorIdx],{'-1':XHR.responseText.split(';')});
			hColorIdx++;
			if (hColorIdx==hColors.length) hColorIdx=0;
		}
	}
	XHR.send(null);
}
*/
var isNav = (navigator.appName.indexOf("Netscape") !=-1);
function ajaxProteinData(e,protID,refData,action,midX,botY,divNum) { // action=ajaxPepData,ajaxProtStat
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
	var trueAction=(action=='ajaxPepData')? quantifAjaxPepAction[qData[0]] : action;
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT="+trueAction+"&id_quantif="+qData[0]+"&ratio="+qData[1]+"&id_prot="+protID+"&pepType=$dispPepType",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			if (action=='ajaxPepData') {
				var resCode;
				if (divNum && divNum==2) {
					resCode=XHR.responseText.replace(/pepPlot/g,'pepPlot2');
					resCode=resCode.replace(/pepList/g,'pepList2');
				}
				else {resCode=XHR.responseText;}
				var codeParts=resCode.split('#==========#');
				infoDiv.innerHTML=codeParts[0]; // HTML part
				eval(codeParts[1]); // javascript part
				// 2nd pep ratio (called from correlation plot)
				if (typeof(refData)=='object') {
					ajaxProteinData(null,protID,refData[1],action,null,null,2);
				}
			}
			else {infoDiv.innerHTML=XHR.responseText;}
		}
	}
	XHR.send(null);

}
/*
function ajaxUpdateGoTermList(goIdStrg) {
	var termsDiv=document.getElementById('goTermsDIV');
	if (!goIdStrg) {
		termsDiv.innerHTML='';
		return;
	}
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxGoTerms&goStrg="+goIdStrg,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			termsDiv.innerHTML=XHR.responseText;
		}
	}
	XHR.send(null);
}
*/
|;
&promsMod::printAjaxManageSaveProteins($projectID,\%promsPath,'document.protForm.chkProt','parent.selQuantifsFrame.ajaxUpdateRestrict');

print qq
|function getXMLHTTP() {
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
function getDomOffset( Obj, Prop ) {
	var iVal = 0;
	while (Obj && Obj.tagName != 'BODY') {
		eval('iVal += Obj.' + Prop + ';');
		Obj = Obj.offsetParent;
	}
	return iVal;
}

</SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
|;

	print "<DIV id=\"waitDIV\"><CENTER><BR><BR><BR><BR><BR><FONT class=\"title3\"><SPAN id=\"waitSPAN\">Fetching data.</SPAN> Please wait...</FONT><BR><IMG src=\"$promsPath{images}/scrollbarGreen.gif\"></CENTER></DIV>\n";
} # end of Html start document

####>Category filter<####
my ($filterClassName,$filterCatName); #,%restrictProteins
if ($restrictListID) {
	my $refprotList=($listAction eq 'restrict')? \%selectedProteins : \%excludedProteins;
	($filterClassName,$filterCatName,my $type)=$dbh->selectrow_array("SELECT CL.NAME,C.NAME,C.LIST_TYPE FROM CLASSIFICATION CL,CATEGORY C WHERE CL.ID_CLASSIFICATION=C.ID_CLASSIFICATION AND C.ID_CATEGORY=$restrictListID");
	if ($exportType ne 'selected') { # Not compatible with export=selected: same %selectedProteins is used!
		if ($type eq 'SITE') {
			&promsQuantif::fetchSiteList($dbh,$restrictListID,$refprotList,$selModifID);
		}
		else {
			my $sthCP=$dbh->prepare("SELECT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$restrictListID");
			$sthCP->execute;
			while (my ($protID)=$sthCP->fetchrow_array) {
				$refprotList->{$protID}=1;
			}
			$sthCP->finish;
		}
	}
}

####>Fetching protein quantification data<####
my (%quantifValues,%quantifInfo,%proteinInfo,%quantifAllNames,%quantifPTM); #,%pepInfo
#my ($minRatio,$maxRatio)=(10000,0.0001); # Globals updated by &fetchQuantificationData for view: log2 & heatmap
my %parameters=(QUANTIF_FAMILY=>$selQuantifFamily,VIEW=>$view,NUM_PEP_CODE=>$numPepCode,QUANTIF_LIST=>\@selectedQuantifications,SEL_MODIF_ID=>$selModifID);
$parameters{MEASURE}=$dispMeasure if $selQuantifFamily eq 'MQ';
if ($strictFilter) { # Not used by volcano/log2 views
	%{$parameters{'STRICT_FILTER'}}=(P_VALUE=>$dispPvalue,RATIO=>$dispFoldChange,INV_RATIO=>$invDispFoldChange) if $selQuantifFamily eq 'RATIO';
	$parameters{'STRICT_FILTER'}{$numPepCode}=$dispNumPep;
}
if ($view eq 'correlMatrix') {
	$parameters{VERBOSE}=1;
	print "<DIV id=\"correlProgressDIV\">\n<BR><FONT class=\"title3\">Fetching data...";
}
&promsQuantif::fetchQuantificationData($dbh,\%parameters,\%quantifInfo,\%quantifValues,\%proteinInfo,\%selectedProteins,\%excludedProteins);
print " Done.</FONT><BR>\n" if $view eq 'correlMatrix';
print "<FONT style=\"font-weight:bold;color:#DD0000\">$parameters{RETURN}</FONT><BR>\n" if $parameters{RETURN};
my ($minRatio,$maxRatio,$minPvalue)=(1,1,1);
if ($selQuantifFamily eq 'RATIO') {
	$minRatio=$parameters{MIN_RATIO}; # updated by &fetchQuantificationData
	$maxRatio=$parameters{MAX_RATIO};
	$minPvalue=$parameters{MIN_PVALUE};
}
#print "MIN=$minRatio, MAX=$maxRatio, PVAL=$minPvalue<BR> \n";
$minRatio/=1.1; # $minRatio/=($view eq 'heatmap')? 2 : 1.1;
$maxRatio*=1.1; #$maxRatio*=($view eq 'heatmap')? 2 : 1.1;
$minPvalue/=10;

my $dispMinRatio=sprintf "1/%.1f",1/$minRatio;
my $dispMaxRatio=sprintf "%.1f",$maxRatio;

####>Optimal naming<####
%quantifAllNames=&promsQuantif::getDistinctQuantifNames($dbh,$quantifList);
#foreach my $quantif (@selectedQuantifications) {}

####>Find quantif of modification<####
&findModificationQuantifs(\@selectedQuantifications,\%quantifInfo,\%quantifPTM); #updates %quantifPTM

####>Get peptide number for emPAI<####
##if ($selQuantifFamily eq 'EMPAI') {
##	my $sthGetPepNum=$dbh->prepare("SELECT ID_PROTEIN,NUM_PEP FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=?");
##	foreach my $quantiID (keys %quantifInfo) {
##		$sthGetPepNum->execute($quantifInfo{$quantiID}[4]);
##		while (my ($protID,$numPep)=$sthGetPepNum->fetchrow_array ) {
##			$pepInfo{$quantiID}{$protID}=$numPep;
##		}
##	}
##}
$dbh->disconnect;
#exit;

        ############################
############>Displaying results<############
        ############################

####>Recording quantif types<####
print qq
|<SCRIPT LANGUAGE="JavaScript">
var quantifAjaxPepAction={
|if !$exportType;
my $numQuantifications=scalar keys %quantifInfo;
my $quantifCount=0;
foreach my $quantifID (keys %quantifInfo) {
	$quantifCount++;
	my $quantifSotfware=($quantifInfo{$quantifID}[1]->{'SOFTWARE'})? $quantifInfo{$quantifID}[1]->{'SOFTWARE'}[0] : 'myProMS';
	my $pepAction=($quantifSotfware eq 'MQ')? 'ajaxPepMaxQuant' : ($quantifSotfware=~/SWATH|DIA/)? 'ajaxPepSwath' : 'ajaxPepRatio';
	my $comma=($quantifCount < $numQuantifications)? ',' : '';
	print "\t$quantifID:'$pepAction'$comma\n" if !$exportType;
}
print qq
|};
</SCRIPT>
| if !$exportType;

##################
####>Heat map<####
##################
my $filterInfoString='';
#if ($view eq 'heatmap') {
#	print qq
#|<SCRIPT LANGUAGE="JavaScript">
#var HM, dataSetQuantif=[];
#function proteinLink(protID) {
#    //console.log(rowID);
#	ajaxListSelectedProteins(protID,null,'$sortOrder'); // sort is necessary for protID format recognition
#}
#function cellLink(protID,quantifData,status,midX,botY) {
#	if (status == false) {
#		document.getElementById('displayDIV').style.display = 'none';
#		return;
#	}
#	else {
#		ajaxProteinData(null,protID,quantifData,'ajaxPepData',midX,botY);
#	}
#}
#var log2=Math.log(2);
#function displayCellInfo(cellValue) {
#	var dispValue;
#	if (cellValue==null) {dispValue='* no value *';}
#	else {
#		dispValue=Math.exp(cellValue*log2);
#		if (dispValue < 1) {dispValue='1/'+(1/dispValue).toPrecision(2);}
#		else {dispValue=dispValue.toPrecision(2);}
#	}
#    return dispValue;
#}
#window.onload=function() {
#	HM=new heatMap({div:'mainGraphDIV',name:'HM',editableItems:{scope:true,color:true}, //editable:true,
#					colLabelOrientation:'diagonal',
#					entities:{row:'Protein',column:'Quantification'},
#					moveLabel:true,
#					rowOnClick:proteinLink,
#					cellOnMouseOver:displayCellInfo,
#					cellOnClick:cellLink,
#					singleCellSelection:true,
#					normalization:{scope:'row',reference:'zero'}
#					});
#|;
#	###>Quantif labels<###
#	my $quantifStrg='';
#	foreach my $quantif (@selectedQuantifications) {
#		my ($quantifID,$ratioPos)=split('_',$quantif);
#		$quantifStrg.=',' if $quantifStrg;
#		#my ($testStatePos,$refStatePos)=split(/\//,$quantifInfo{$quantifID}[1]->{'RATIOS'}[$ratioPos-1]);
#		#$quantifStrg.="['".&promsMod::shortenName($quantifInfo{$quantifID}[0],15).":".$quantifInfo{$quantifID}[2]->{$testStatePos}{'NAME'}."/".$quantifInfo{$quantifID}[2]->{$refStatePos}{'NAME'}."','$quantif']";
#		$quantifStrg.="['".&promsMod::shortenName($quantifAllNames{OPTIMAL}{$quantif},30)."','$quantif','$quantifAllNames{FULL}{$quantif}']";
#		print "\tdataSetQuantif.push('$quantif');\n";
#	}
#	print "\tHM.setColumns([$quantifStrg]);\n";
#
#	###>Protein labels & data<###
#	$filterInfoString='For each protein, all quantifications are displayed if at least 1 passes filtering. Ratios are centered on fold-change=1'; # (&infin; ratios are ignored)'; # ∞
#	my $log2=log(2);
#	my $numQuantifs=scalar @selectedQuantifications;
#	#foreach my $protID (sort{lc($proteinInfo{$a}[0]) cmp lc($proteinInfo{$a}[0])} keys %proteinInfo) { #}
#	foreach my $protID (sort{&sortProteins($sortOrder,\%proteinInfo,\%quantifValues)} keys %proteinInfo) {
#		my $protDataStrg='';
#		my $i=0;
#		my $okFilters=0;
#		foreach my $quantif (@selectedQuantifications) {
#			#next if (!$quantifValues{$protID}{$quantif} || $quantifValues{$protID}{$quantif}{$numPepCode} < $dispNumPep || ($dispPvalue < 1 && (!$quantifValues{$protID}{$quantif}{'P_VALUE'} || $quantifValues{$protID}{$quantif}{'P_VALUE'} > $dispPvalue)));
#			#if ($quantifValues{$protID}{$quantif} && $quantifValues{$protID}{$quantif}{'RATIO'} > $MIN_INF_RATIO && $quantifValues{$protID}{$quantif}{'RATIO'} < $MAX_INF_RATIO) { #}
#			if ($quantifValues{$protID}{$quantif}) {
#				if ($quantifValues{$protID}{$quantif}{'RATIO'} >= $MAX_INF_RATIO) {$quantifValues{$protID}{$quantif}{'RATIO'}=$maxRatio;}
#				elsif ($quantifValues{$protID}{$quantif}{'RATIO'} < $MIN_INF_RATIO) {$quantifValues{$protID}{$quantif}{'RATIO'}=$minRatio;}
#				$okFilters=1 if (!$okFilters && ($quantifValues{$protID}{$quantif}{'RATIO'} <= $invDispFoldChange || $quantifValues{$protID}{$quantif}{'RATIO'} >= $dispFoldChange) && $quantifValues{$protID}{$quantif}{$numPepCode} >= $dispNumPep && ($dispPvalue==1 || ($quantifValues{$protID}{$quantif}{'P_VALUE'} && $quantifValues{$protID}{$quantif}{'P_VALUE'} <= $dispPvalue)));
#				$protDataStrg.=sprintf "%.2f",log($quantifValues{$protID}{$quantif}{'RATIO'})/$log2;
#			}
#			$protDataStrg.=',' if ++$i < $numQuantifs;
#		}
#		print "\tHM.addRow(['$proteinInfo{$protID}[0]',$protID],[$protDataStrg]);\n" if $okFilters;
#	}
#
#	print qq
#|	HM.draw();
#	document.getElementById('waitDIV').style.display='none';
#	document.getElementById('resultDIV').style.visibility='visible';
#
#}
#</SCRIPT>
#|;
#} # end of heatmap


######################
####>Volcano Plot<####
######################
if ($view eq 'volcano') { #elsif
	print qq
|<SCRIPT LANGUAGE="JavaScript">
var VP,dataSetQuantif=[];
function proteinLink(dSetIdx,identifier,protID) {
	ajaxProteinData(null,protID,dataSetQuantif[dSetIdx],'ajaxPepData');
}
window.onload=function() {
	VP=new volcanoPlot({div:'mainGraphDIV',width:700,height:500,name:'VP',
						foldChange:$dispFoldChange,pValue:$dispPvalue,
						pointOnclick:proteinLink,
						pointOnList:ajaxListSelectedProteins,
						allowHighlight:false,
						exportAsImage:['Export as image','VolcanoPlot','./exportSVG.cgi']
						});
|;
	###>Quantif sets<###
	my $quantifIdx=0;
	foreach my $quantif (@selectedQuantifications) {
		my ($quantifID,$ratioPos)=split('_',$quantif);
		my ($testStatePos,$refStatePos)=split(/\//,$quantifInfo{$quantifID}[1]->{'RATIOS'}[$ratioPos-1]);
		$testStatePos=~s/%\d+//; $refStatePos=~s/%\d+//;
		print "\tdataSetQuantif[$quantifIdx]='$quantif';\n";
		print "\tVP.addDataSet($quantifIdx,{name:'",&promsMod::shortenName($quantifInfo{$quantifID}[0],15).":".$quantifInfo{$quantifID}[2]->{$testStatePos}{'NAME'}."/".$quantifInfo{$quantifID}[2]->{$refStatePos}{'NAME'},"',sizeName:'Peptides'});\n";
		$quantifIdx++;
	}

	###>Quantif data<###
	$filterInfoString=($dispNumPep==1)? "All $proteoform ratios are shown" : "Only $proteoform ratios computed with at least $dispNumPep peptide(s) are shown";
	$quantifIdx=0;
	foreach my $quantif (@selectedQuantifications) {
		my $count=0;
		foreach my $modProtID (keys %{$quantifValues{$quantif}}) {
			#next if ($quantifValues{$quantif}{$modProtID}{$numPepCode} < $dispNumPep || !$quantifValues{$quantif}{$modProtID}{'P_VALUE'}); #!defined $quantifValues{$quantif}{$modProtID}{'P_VALUE'} || $quantifValues{$quantif}{$modProtID}{'P_VALUE'} > $dispPvalue ||
			next if $quantifValues{$quantif}{$modProtID}{$numPepCode} < $dispNumPep;
			if (!defined $quantifValues{$quantif}{$modProtID}{'P_VALUE'}) {$quantifValues{$quantif}{$modProtID}{'P_VALUE'}=1;}
			elsif ($quantifValues{$quantif}{$modProtID}{'P_VALUE'}==0) {$quantifValues{$quantif}{$modProtID}{'P_VALUE'}=$minPvalue;}
			my ($protID,$modStrg)=($modProtID=~/^(\d+)(.*)/); $modStrg='' unless $modStrg;
			# no filter on fold change
			if ($count==0) {
				print "\tVP.addDataAsString($quantifIdx,'";
			}
			else {print ";";}
			$count++;
			print "$proteinInfo{$protID}[0]$modStrg,$modProtID,";
			#>Set +/- infinite ratios
			if ($quantifValues{$quantif}{$modProtID}{RATIO} == $MAX_INF_RATIO_DB) {
				$proteinInfo{$protID}[-1]=1000 unless $proteinInfo{$protID}[-1]; # set prot length
				printf "+,%.2f",100*$quantifValues{$quantif}{$modProtID}{$numPepCode}/$proteinInfo{$protID}[-1]; # num pep/100 res
			}
			elsif ($quantifValues{$quantif}{$modProtID}{RATIO} == $MIN_INF_RATIO_DB) {
				$proteinInfo{$protID}[-1]=1000 unless $proteinInfo{$protID}[-1]; # set prot length
				printf "-,%.2f",100*$quantifValues{$quantif}{$modProtID}{$numPepCode}/$proteinInfo{$protID}[-1]; # num pep/100 res
			}
			else { # normal ratios
				$quantifValues{$quantif}{$modProtID}{RATIO}=1 if !$quantifValues{$quantif}{$modProtID}{P_VALUE};
				printf "$quantifValues{$quantif}{$modProtID}{RATIO},%.2e",$quantifValues{$quantif}{$modProtID}{P_VALUE};
			}
			print ",$quantifValues{$quantif}{$modProtID}{$numPepCode}";
			if ($count==100) {
				print "');\n";
				$count=0;
			}
		}
		if ($count > 0) {
			print "');\n";
		}
		$quantifIdx++;
	}
	print qq
|	VP.draw();
	document.getElementById('waitDIV').style.display='none';
	document.getElementById('resultDIV').style.visibility='visible';
}
</SCRIPT>
|;
} # end of volcano

###############################
####>Log2 correlation plot<####
###############################
elsif ($view eq 'log2') {

	###>Reordering quantifs (whole protein before ptm)<###
	if (scalar keys %quantifPTM) { # at least 1 PTM quantif
		my (@std,@ptm);
		foreach my $quantif (@selectedQuantifications) {
			if ($quantifPTM{$quantif}) {push @ptm,$quantif;}
			else {push @std,$quantif;}
		}
		@selectedQuantifications=@std;
		push @selectedQuantifications,@ptm;
	}
	my ($set1Name,$set2Name)=(scalar @selectedQuantifications > 2 || !$quantifAllNames{OPTIMAL}{$selectedQuantifications[0]})? ('dataset 1','dataset 2') : (&promsMod::shortenName($quantifAllNames{OPTIMAL}{$selectedQuantifications[0]},50),&promsMod::shortenName($quantifAllNames{OPTIMAL}{$selectedQuantifications[1]},50));
	$set1Name=~s/'/ /g; $set2Name=~s/'/ /g;
	my ($logText,$logBase,$protLink)=($selQuantifFamily eq 'RATIO')? ('Log2',log(2),'pointOnclick:proteinLink,'): ('Log10',log(10),''); # Avoided for emPAI and MaxQuant quanti values
	###>javaScript<###
	print qq
|<SCRIPT type="text/javascript">
var GP,dataSetQuantif=[];
function proteinLink(dSetIdx,identifier,protID) {
	ajaxProteinData(null,protID,dataSetQuantif[dSetIdx],'ajaxPepData');
}
window.onload=function() {
	GP=new genericPlot({div:'mainGraphDIV',name:'GP',width:600,height:600,
						axisX:{title:'$logText($set1Name)'}, // ,zoomable:true
						axisY:{title:'$logText($set2Name)'}, // ,zoomable:true
						sameScale:true,
						axisClosure:true,
						zoomable:true,
						//hideDataSets:true, // hide dataset label even if more than 1
						pointOpacity:0.6,
						pointLabel:myPointLabel,
						$protLink
						pointOnList:ajaxListSelectedProteins,
						convertValue:function(axis,val){return Math.log(val)/$logBase;}, // only for threshold lines
						exportAsImage:['Export as image','CorrelationPlot','./exportSVG.cgi']
						});
|;
	if ($selQuantifFamily eq 'RATIO') {
		print qq
|	GP.addThreshold({axis:'X',label:'max. fold change',value:$dispFoldChange,color:'#00A000',keepAbove1:true,roundValue:1});
	GP.addThreshold({axis:'X',label:'min. fold change',value:$invDispFoldChange,color:'#00A000',keepAbove1:true,roundValue:1});
	GP.addThreshold({axis:'Y',label:'max. fold change',value:$dispFoldChange,color:'#00A000',keepAbove1:true,roundValue:1});
	GP.addThreshold({axis:'Y',label:'min. fold change',value:$invDispFoldChange,color:'#00A000',keepAbove1:true,roundValue:1});
|;
	}

	###>Quantif sets<###
	my $minusInfX=my $minusInfY=my $plusInfX=my $plusInfY=0;
	my $quantifIdx=0;
	foreach my $quantifIdx1 (0..$#selectedQuantifications-1) {
		if ($selQuantifFamily eq 'RATIO') {
			$minusInfX=1 if $parameters{MINUS_INF}{$selectedQuantifications[$quantifIdx1]};
			$plusInfX=1 if $parameters{PLUS_INF}{$selectedQuantifications[$quantifIdx1]};
		}
		foreach my $quantifIdx2 (($quantifIdx1+1)..$#selectedQuantifications) {
			if ($selQuantifFamily eq 'RATIO') {
				$minusInfY=1 if $parameters{MINUS_INF}{$selectedQuantifications[$quantifIdx2]};
				$plusInfY=1 if $parameters{PLUS_INF}{$selectedQuantifications[$quantifIdx2]};
			}
			#my ($quantifID1,$ratioPos1)=split('_',$selectedQuantifications[$quantifIdx1]);
			#my ($quantifID2,$ratioPos2)=split('_',$selectedQuantifications[$quantifIdx2]);
			#my ($testStatePos1,$refStatePos1)=split(/\//,$quantifInfo{$quantifID1}[1]->{'RATIOS'}[$ratioPos1-1]);
			#my ($testStatePos2,$refStatePos2)=split(/\//,$quantifInfo{$quantifID2}[1]->{'RATIOS'}[$ratioPos2-1]);
			print "\tdataSetQuantif[$quantifIdx]=['$selectedQuantifications[$quantifIdx1]','$selectedQuantifications[$quantifIdx2]'];\n";
			#print "\tGP.addDataSet($quantifIdx,{name:'$quantifInfo{$quantifID1}[0] vs $quantifInfo{$quantifID2}[0]',axisX:'X',axisY:'Y'});\n";
			my $optimalName1=($quantifAllNames{OPTIMAL}{$selectedQuantifications[$quantifIdx1]})? $quantifAllNames{OPTIMAL}{$selectedQuantifications[$quantifIdx1]} : $quantifInfo{$selectedQuantifications[$quantifIdx1]}[0];
			my $optimalName2=($quantifAllNames{OPTIMAL}{$selectedQuantifications[$quantifIdx2]})? $quantifAllNames{OPTIMAL}{$selectedQuantifications[$quantifIdx2]} : $quantifInfo{$selectedQuantifications[$quantifIdx2]}[0];
			print "\tGP.addDataSet($quantifIdx,{name:'",&promsMod::shortenName($optimalName1,50)," vs ",&promsMod::shortenName($optimalName2,50),"',axisX:'X',axisY:'Y'});\n";
			$quantifIdx++;
		}
	}
	#>Infinite ratios thresholds
	if ($selQuantifFamily eq 'RATIO') {
		print "\tGP.addThreshold({axis:'X',label:'1/$infChar ratios',value:$dispMinRatio,color:'#FF0000',keepAbove1:true,roundValue:1,pattern:'- '});\n" if $minusInfX;
		print "\tGP.addThreshold({axis:'X',label:'$infChar ratios',value:$dispMaxRatio,color:'#FF0000',keepAbove1:true,roundValue:1,pattern:'- '});\n" if $plusInfX;
		print "\tGP.addThreshold({axis:'Y',label:'1/$infChar ratios',value:$dispMinRatio,color:'#FF0000',keepAbove1:true,roundValue:1,pattern:'- '});\n" if $minusInfY;
		print "\tGP.addThreshold({axis:'Y',label:'$infChar ratios',value:$dispMaxRatio,color:'#FF0000',keepAbove1:true,roundValue:1,pattern:'- '});\n" if $plusInfY;
	}

	###>Quantif data<###$selectedQuantifications[$quantifIdx1]
	#$filterInfoString=($dispNumPep==1)? 'All comparable protein ratios are shown' : "Protein ratios computed with less than $dispNumPep peptides are not shown";
	##$filterInfoString.=' (&infin; ratios are ignored)';
	##$filterInfoString.=" (1/&infin; and &infin; ratios are set to $dispMinRatio and $dispMaxRatio respectively)";
	$filterInfoString=($strictFilter)? "Strict filtering: A $proteoform point is displayed if both quantifications match red filters" : "Loose filtering: A $proteoform point is displayed if 1 of the 2 quantifications matches red filters";
	my (%protVenn,%proteinOK,%siteVenn);
	$quantifIdx=0;
	if ($selQuantifFamily eq 'RATIO') {
		foreach my $quantifIdx1 (0..$#selectedQuantifications-1) {
			my $quantif1=$selectedQuantifications[$quantifIdx1];
			foreach my $quantifIdx2 (($quantifIdx1+1)..$#selectedQuantifications) {
				my $quantif2=$selectedQuantifications[$quantifIdx2];
				my $count=0;
				foreach my $modProtID2 (keys %{$quantifValues{$quantif2}}) {
					my ($protID,$modStrg)=($modProtID2=~/^(\d+)(.*)/); $modStrg='' unless $modStrg;
					my $modProtID1;
					if ($modStrg) {
						if ($quantifPTM{$quantif1}) {$modProtID1=$modProtID2;} # ptm(test) vs ptm(ref)
						else {$modProtID1=$protID;} # ptm(test) vs whole(ref)
					}
					else {$modProtID1=$protID;} # must be whole(test) vs whole(ref)
					#>Prot #1
					if (!$proteinOK{$modProtID1} || !$proteinOK{$modProtID1}{$quantif1}) {
						if (!$quantifValues{$quantif1}{$modProtID1} || !$quantifValues{$quantif1}{$modProtID1}{'RATIO'}) {$proteinOK{$modProtID1}{$quantif1}=-2;}
						else {
							#>Ratio
							if ($quantifValues{$quantif1}{$modProtID1}{RATIO} == $MAX_INF_RATIO_DB) {$quantifValues{$quantif1}{$modProtID1}{RATIO}=$maxRatio*(1+rand(1));}
							elsif ($quantifValues{$quantif1}{$modProtID1}{RATIO} == $MIN_INF_RATIO_DB) {$quantifValues{$quantif1}{$modProtID1}{RATIO}=$minRatio/(1+rand(1));}
							#>P-value
							$proteinOK{$modProtID1}{$quantif1}=($quantifValues{$quantif1}{$modProtID1}{'P_VALUE'} && $quantifValues{$quantif1}{$modProtID1}{'P_VALUE'} > $dispPvalue)? -1 : 1;
							#>Num Peptides
							$proteinOK{$modProtID1}{$quantif1}=-1 if $quantifValues{$quantif1}{$modProtID1}{$numPepCode} < $dispNumPep;
						}
					}
					#>Prot #2
					if (!$proteinOK{$modProtID2} || !$proteinOK{$modProtID2}{$quantif2}) {
						if (!$quantifValues{$quantif2}{$modProtID2}{'RATIO'}) {$proteinOK{$modProtID2}{$quantif2}=-2;} # Just to be safe
						else {
							#>Ratio
							if ($quantifValues{$quantif2}{$modProtID2}{RATIO} == $MAX_INF_RATIO_DB) {$quantifValues{$quantif2}{$modProtID2}{RATIO}=$maxRatio*(1+rand(1));}
							elsif ($quantifValues{$quantif2}{$modProtID2}{RATIO} == $MIN_INF_RATIO_DB) {$quantifValues{$quantif2}{$modProtID2}{RATIO}=$minRatio/(1+rand(1));}
							#>P-value
							$proteinOK{$modProtID2}{$quantif2}=($quantifValues{$quantif2}{$modProtID2}{'P_VALUE'} && $quantifValues{$quantif2}{$modProtID2}{'P_VALUE'} > $dispPvalue)? -1 : 1;
							#>Num Peptides
							$proteinOK{$modProtID2}{$quantif2}=-1 if $quantifValues{$quantif2}{$modProtID2}{$numPepCode} < $dispNumPep;
						}
					}
					$protVenn{$protID}{$quantifIdx1}=1 if ($proteinOK{$modProtID1}{$quantif1} > 0 || (!$strictFilter && $proteinOK{$modProtID1}{$quantif1}>=-1 && $proteinOK{$modProtID2}{$quantif2} > 0));
					$protVenn{$protID}{$quantifIdx2}=1 if ($proteinOK{$modProtID2}{$quantif2} > 0 || (!$strictFilter && $proteinOK{$modProtID2}{$quantif2}>=-1 && $proteinOK{$modProtID1}{$quantif1} > 0));
					next if ($proteinOK{$modProtID1}{$quantif1}==-2 || $proteinOK{$modProtID2}{$quantif2}==-2);
					next if (($strictFilter && ($proteinOK{$modProtID1}{$quantif1} < 0 || $proteinOK{$modProtID2}{$quantif2} < 0)) || (!$strictFilter && $proteinOK{$modProtID1}{$quantif1} < 0 && $proteinOK{$modProtID2}{$quantif2} < 0));
					my $numPepMin=($quantifValues{$quantif1}{$modProtID1}{$numPepCode} < $quantifValues{$quantif2}{$modProtID2}{$numPepCode})? $quantifValues{$quantif1}{$modProtID1}{$numPepCode} : $quantifValues{$quantif2}{$modProtID2}{$numPepCode};
					if ($count==0) {
						print "\tGP.addDataAsString($quantifIdx,'";
					}
					else {print ";";}
					$count++;
					print "$proteinInfo{$protID}[0]$modStrg,$modProtID2,";
					printf "%.2f",log($quantifValues{$quantif1}{$modProtID1}{RATIO})/$logBase;
					print ',';
					printf "%.2f",log($quantifValues{$quantif2}{$modProtID2}{RATIO})/$logBase;
					print ",$numPepMin";
					#printf ("%.0f",$quantifValues{$selectedQuantifications[$quantifIdx1]}{$modProtID1}{$numPepCode}+$quantifValues{$selectedQuantifications[$quantifIdx2]}{$modProtID2}{$numPepCode});
					if ($count==100) {
						print "');\n";
						$count=0;
					}
				}
				if ($count > 0) {
					print "');\n";
				}
				$quantifIdx++;
			}
			#>Complete venn diagram for unique proteins in first quantif in list
			if ($quantifIdx1==0 && $numSelectedQuantifs<=5) {
				foreach my $modProtID1 (keys %{$quantifValues{$quantif1}}) {
					next if ($proteinOK{$modProtID1} && $proteinOK{$modProtID1}{$quantif1}); # never found in other quantif
					next if ($quantifValues{$quantif1}{$modProtID1}{'P_VALUE'} && $quantifValues{$quantif1}{$modProtID1}{'P_VALUE'} > $dispPvalue); # P-value
					next if $quantifValues{$quantif1}{$modProtID1}{$numPepCode} < $dispNumPep; # Num Peptides
					my ($protID,$modStrg)=($modProtID1=~/^(\d+)(.*)/); # $modStrg='' unless $modStrg;
					$protVenn{$protID}{$quantifIdx1}=1;
				}
			}
		}
	}
	elsif ($selQuantifFamily eq 'MQ') {
		foreach my $quantifIdx1 (0..$#selectedQuantifications-1) {
			my $quantif1=$selectedQuantifications[$quantifIdx1];
			foreach my $quantifIdx2 (($quantifIdx1+1)..$#selectedQuantifications) {
				my $quantif2=$selectedQuantifications[$quantifIdx2];
				my $count=0;
				foreach my $modProtID2 (keys %{$quantifValues{$quantif2}}) {
					my ($protID,$modStrg)=($modProtID2=~/^(\d+)(.*)/); $modStrg='' unless $modStrg;
					my $modProtID1;
					if ($modStrg) {
						if ($quantifPTM{$quantif1}) {$modProtID1=$modProtID2;} # ptm(test) vs ptm(ref)
						else {$modProtID1=$protID;} # ptm(test) vs whole(ref)
					}
					else {$modProtID1=$protID;} # must be whole(test) vs whole(ref)
					#>Prot #1
					if (!$proteinOK{$modProtID1} || !$proteinOK{$modProtID1}{$quantif1}) {
						if (!$quantifValues{$quantif1}{$modProtID1} || !$quantifValues{$quantif1}{$modProtID1}{$dispMeasure}) {$proteinOK{$modProtID1}{$quantif1}=-2;}
						else {$proteinOK{$modProtID1}{$quantif1}=($quantifValues{$quantif1}{$modProtID1}{$numPepCode} < $dispNumPep)? -1 : 1;}#>Num Peptides
					}
					#>Prot #2
					if (!$proteinOK{$modProtID2} || !$proteinOK{$modProtID2}{$quantif2}) {
						if (!$quantifValues{$quantif2}{$modProtID2} || !$quantifValues{$quantif2}{$modProtID2}{$dispMeasure}) {$proteinOK{$modProtID2}{$quantif2}=-2;}
						else {$proteinOK{$modProtID2}{$quantif2}=($quantifValues{$quantif2}{$modProtID2}{$numPepCode} < $dispNumPep)? -1 : 1;}#>Num Peptides
					}
					$protVenn{$protID}{$quantifIdx1}=1 if ($proteinOK{$modProtID1}{$quantif1} > 0 || (!$strictFilter && $proteinOK{$modProtID1}{$quantif1}>=-1 && $proteinOK{$modProtID2}{$quantif2} > 0));
					$protVenn{$protID}{$quantifIdx2}=1 if ($proteinOK{$modProtID2}{$quantif2} > 0 || (!$strictFilter && $proteinOK{$modProtID2}{$quantif2}>=-1 && $proteinOK{$modProtID1}{$quantif1} > 0));
					next if ($proteinOK{$modProtID1}{$quantif1}==-2 || $proteinOK{$modProtID2}{$quantif2}==-2);
					next if (($strictFilter && ($proteinOK{$modProtID1}{$quantif1} < 0 || $proteinOK{$modProtID2}{$quantif2} < 0)) || (!$strictFilter && $proteinOK{$modProtID1}{$quantif1} < 0 && $proteinOK{$modProtID2}{$quantif2} < 0));
					my $numPepMin=($quantifValues{$quantif1}{$modProtID1}{$numPepCode} < $quantifValues{$quantif2}{$modProtID2}{$numPepCode})? $quantifValues{$quantif1}{$modProtID1}{$numPepCode} : $quantifValues{$quantif2}{$modProtID2}{$numPepCode};
					if ($count==0) {
						print "\tGP.addDataAsString($quantifIdx,'";
					}
					else {print ";";}
					$count++;
					print "$proteinInfo{$protID}[0]$modStrg,$modProtID2,";
					printf "%.2f",log($quantifValues{$quantif1}{$modProtID1}{$dispMeasure})/$logBase;
					print ',';
					printf "%.2f",log($quantifValues{$quantif2}{$modProtID2}{$dispMeasure})/$logBase;
					print ",$numPepMin";
					if ($count==100) {
						print "');\n";
						$count=0;
					}
				}
				if ($count > 0) {
					print "');\n";
				}
				$quantifIdx++;
			}
			#>Complete venn diagram for unique proteins in first quantif in list
			if ($quantifIdx1==0 && $numSelectedQuantifs<=5) {
				foreach my $modProtID1 (keys %{$quantifValues{$quantif1}}) {
					next if ($proteinOK{$modProtID1} && $proteinOK{$modProtID1}{$quantif1}); # never found in other quantif
					next if $quantifValues{$quantif1}{$modProtID1}{$numPepCode} < $dispNumPep; # Num Peptides
					my ($protID,$modStrg)=($modProtID1=~/^(\d+)(.*)/); # $modStrg='' unless $modStrg;
					$protVenn{$protID}{$quantifIdx1}=1;
				}
			}
		}
	}
	elsif ($selQuantifFamily =~ /EMPAI|SIN/) {
		foreach my $quantifIdx1 (0..$#selectedQuantifications-1) {
			my $quantif1=$selectedQuantifications[$quantifIdx1];
			foreach my $quantifIdx2 (($quantifIdx1+1)..$#selectedQuantifications) {
				my $quantif2=$selectedQuantifications[$quantifIdx2];
				my $count=0;
				#foreach my $protID (keys %quantifValues) {#}
				foreach my $protID (keys %{$quantifValues{$quantif2}}) {
					if (!$proteinOK{$protID} || !$proteinOK{$protID}{$quantif1}) {
						if (!$quantifValues{$quantif1}{$protID} || !$quantifValues{$quantif1}{$protID}{$empaiValue}) {$proteinOK{$protID}{$quantif1}=-2;}
						else {
							$proteinOK{$protID}{$quantif1}=1;
							$protVenn{$protID}{$quantifIdx1}=1; # no data filtering?
						}
					}
					if (!$proteinOK{$protID}{$quantif2}) {
						if (!$quantifValues{$quantif2}{$protID} || !$quantifValues{$quantif2}{$protID}{$empaiValue}) {$proteinOK{$protID}{$quantif2}=-2;}
						else {
							$proteinOK{$protID}{$quantif2}=1;
							$protVenn{$protID}{$quantifIdx2}=1; # no data filtering?
						}
					}
					next if ($proteinOK{$protID}{$quantif1}==-2 || $proteinOK{$protID}{$quantif2}==-2);

					my $numPepMin=($quantifValues{$quantif1}{$protID}{'PEPTIDES'} < $quantifValues{$quantif2}{$protID}{'PEPTIDES'})? $quantifValues{$quantif1}{$protID}{'PEPTIDES'} : $quantifValues{$quantif2}{$protID}{'PEPTIDES'};

					####>Data
					###next if (!$quantifValues{$protID}{$selectedQuantifications[$quantifIdx1]}{$empaiValue} || !$quantifValues{$protID}{$selectedQuantifications[$quantifIdx2]}{$empaiValue});
					###my ($numPepQ1,$numPepQ2)=($pepInfo{$selectedQuantifications[$quantifIdx1]}{$protID},$pepInfo{$selectedQuantifications[$quantifIdx2]}{$protID});
					###next if (!$numPepQ1 || !$numPepQ2);
					###my $numPepMin=($numPepQ1 < $numPepQ2) ? $numPepQ1 : $numPepQ2 ;

					if ($count==0) {
						print "\tGP.addDataAsString($quantifIdx,'";
					}
					else {print ";";}
					$count++;
					print "$proteinInfo{$protID}[0],$protID,";
					printf "%.2f",log($quantifValues{$quantif1}{$protID}{$empaiValue})/$logBase;
					print ',';
					printf "%.2f",log($quantifValues{$quantif2}{$protID}{$empaiValue})/$logBase;
					print ",$numPepMin";
					if ($count==100) {
						print "');\n";
						$count=0;
					}
				}
				if ($count > 0) {
					print "');\n";
				}
				$quantifIdx++;
			}
			#>Complete venn diagram for unique proteins in first quantif in list
			if ($quantifIdx1==0 && $numSelectedQuantifs<=5) {
				foreach my $protID (keys %{$quantifValues{$quantif1}}) {
					next if ($proteinOK{$protID} && $proteinOK{$protID}{$quantif1}); # never found in other quantif
					$protVenn{$protID}{$quantifIdx1}=1;
				}
			}
		}
	}
	print qq
|	GP.draw();
	document.getElementById('waitDIV').style.display='none';
	document.getElementById('resultDIV').style.visibility='visible';
|;

	if ($numSelectedQuantifs <= 5) {
		my %vennDiagram;
		foreach my $protID (keys %protVenn) {
			my $vennStrg='';
			foreach my $quantifIdx (sort{$a<=>$b} keys %{$protVenn{$protID}}) {
				$vennStrg.=chr(65+$quantifIdx); # A,B,C,D,E
			}
			$vennDiagram{$vennStrg}++;
		}

		print "groupLabelsObj={";
		my $count=0;
		my $legendMaxSize=0;
		foreach my $quantif (@selectedQuantifications) {
			#print "$quantifCode{$quantif}:'$quantifLabel{$quantif}'";
			my $optimalName=($quantifAllNames{OPTIMAL}{$quantif})? $quantifAllNames{OPTIMAL}{$quantif} : $quantifInfo{$quantif}[0];
			print chr(65+$count).":'$optimalName'";
			$count++;
			print ',' if $count < $numSelectedQuantifs;
			if ($legendMaxSize < length($optimalName)) {$legendMaxSize=length($optimalName);}
		}
		my $legendPos=($legendMaxSize < 50)? 'side' : 'bottom';
		print "};\n\tnew VennDiagram({div:'vennDIV',legendPosition:'$legendPos',size:300,groupLabels:groupLabelsObj,setValues:{";
		$count=0;
		my $numSets=scalar keys %vennDiagram;
		foreach my $set (sort keys %vennDiagram) {
			$count++;
			print "$set:$vennDiagram{$set}";
			print ',' if $count < $numSets;
		}
		print "}});\n";
	}

	print qq
|}
function myPointLabel(dp,type) {
	var infoStrg=dp.label;
	if (type=='min') {
		return dp.label;
	}
	else {
		if (dp.dataSet.params.chart.dataSets.length > 1) {infoStrg+='\\nSet='+dp.dataSet.params.name;}
		var fc1=Math.exp(dp.getX()*$logBase);
		var fc2=Math.exp(dp.getY()*$logBase);
		//fc1=(fc1 < 1)? fc1='1/'+Math.round(100/fc1)/100 : Math.round(fc1*100)/100; // round with 2 decimals
		//fc2=(fc2 < 1)? fc2='1/'+Math.round(100/fc2)/100 : Math.round(fc2*100)/100; // round with 2 decimals
|;
		if ($selQuantifFamily =~ /EMPAI|MQ/) {
			my $description = ($selQuantifFamily eq 'EMPAI') ?  $emPAIMeasures{$empaiValue} : $maxQuantMeasures{$dispMeasure};
			print qq
|		infoStrg+='\\n$description='+fc1.toExponential(2)+' vs '+fc2.toExponential(2);
|;
		}
		else{
			print qq
|		if (fc1 < 1) {fc1=(fc1 < $minRatio*1.05)? fc1='1/$infChar' : '1/'+Math.round(100/fc1)/100;} // round with 2 decimals
		else {fc1=(fc1 > $maxRatio/1.05)? fc1='$infChar' : Math.round(fc1*100)/100;}
		if (fc2 < 1) {fc2=(fc2 < $minRatio*1.05)? fc2='1/$infChar' : '1/'+Math.round(100/fc2)/100;} // round with 2 decimals
		else {fc2=(fc2 > $maxRatio/1.05)? fc2='$infChar' : Math.round(fc2*100)/100;}

		infoStrg+='\\nFC='+fc1+' vs '+fc2;
|;
		}

			print qq
|	}
	infoStrg+='\\nMin. peptides='+dp.size;
	return infoStrg;
}
</SCRIPT>
|;
}


############################
####>Correlation values<####
############################
elsif ($view eq 'correlMatrix') {
	$filterInfoString=($strictFilter)? "Strict filtering: For each $proteoform, only quantifications that match red filters are included" : "Loose filtering: For each $proteoform, all quantifications are included if one of them matches red filters";
	print qq
|<SCRIPT type="text/javascript">
document.getElementById('waitSPAN').innerHTML='Computing correlations.';
</SCRIPT>
|;
	print "<BR><FONT class=\"title3\">Computing correlations and distances...";
	mkdir "$promsPath{tmp}/scratch" unless -d "$promsPath{tmp}/scratch";
	&promsMod::cleanDirectory("$promsPath{tmp}/scratch/correlation",'1h');

	my $tmpDir = strftime("%Y%m%d%H%M%S",localtime);
	my $workDir="$promsPath{tmp}/scratch/correlation/$tmpDir";
	mkdir $workDir;

	####>Writing data to files<####
	my $dataFile1="table1.txt"; # No +/-inf
	my $dataFile2="table2.txt";
	my $dataFile3="table3.txt";
	my %correlFiles=(1=>$dataFile1,3=>$dataFile3);
	my %jobInfo=(
		2=>{TITLE=>'Fold change tendancy correlation',
			DETAILS=>qq
|&bull;Values are splited in 3 bins (<B>]-&infin; ~ 1/$dispFoldChange]</B>, <B>]1/$dispFoldChange ~ $dispFoldChange\[</B>, <B>[$dispFoldChange ~ +&infin;[</B>) for correlation.<BR>
&bull;Based only on values from ${proteoform}s shared by both quantifications.<BR>
&bull;&plusmn;&infin; ratios are allowed.|
		},
		3=>{TITLE=>'Quantification coverage distance',
			DETAILS=>"&bull;Distance is based of presence or absence of value for each $proteoform in compared quantifications."
		}
	);
	my %protInQuantifs; #,%protUnion
	open (DATA1,">$workDir/$dataFile1");
	open (DATA3,">$workDir/$dataFile3");
	#>Header
	print DATA1 join("\t",('PROTEIN_ID',@selectedQuantifications)),"\n";
	print DATA3 join("\t",('PROTEIN_ID',@selectedQuantifications)),"\n";
	#>Data
	my ($log2,$log10)=(log(2),log(10));
	if ($selQuantifFamily eq 'RATIO') {
		$correlFiles{2}=$dataFile2;
		$jobInfo{1}={TITLE=>'Fold change correlation',
					 DETAILS=>qq
|&bull;Based only on values from ${proteoform}s shared by both quantifications.<BR>
&bull;&plusmn;&infin; ratios are excluded.|
					};
		$jobInfo{3}{DETAILS}.="<BR>\n&bull;&plusmn;&infin; ratios are allowed.";

		open (DATA2,">$workDir/$dataFile2");
		print DATA2 join("\t",('PROTEIN_ID',@selectedQuantifications)),"\n";
		PROT: foreach my $modProtID (keys %quantifValues) {
			#>Filtering data
			my $usedQuantifs=0;
			my $usedQuantifs2=0;
			my $goodQuantif=0;
			my $goodQuantif2=0;
			my (@ratioList,@catList);
			foreach my $quantif (@selectedQuantifications) {
				my $ratio=my $cat='NA';
				if (defined $quantifValues{$modProtID}{$quantif}) {
					my $okFilters=0;
					$okFilters++ if $quantifValues{$modProtID}{$quantif}{$numPepCode} >= $dispNumPep;
					$okFilters++ if ($dispPvalue == 1 || (defined $quantifValues{$modProtID}{$quantif}{'P_VALUE'} && $quantifValues{$modProtID}{$quantif}{'P_VALUE'} <= $dispPvalue));
					if (!$strictFilter || $okFilters==2) {
						#>Ratio (excluding +/-Inf)
						if ($quantifValues{$modProtID}{$quantif}{'RATIO'} != 0.001 && $quantifValues{$modProtID}{$quantif}{'RATIO'} != 1000) {
							$ratio=$quantifValues{$modProtID}{$quantif}{'RATIO'};
							$usedQuantifs++;
							$goodQuantif=1 if $okFilters==2;
						}
						#>Tendancy (including +/-Inf)
						$cat=($quantifValues{$modProtID}{$quantif}{'RATIO'} <= $invDispFoldChange)? 1 : ($quantifValues{$modProtID}{$quantif}{'RATIO'} < $dispFoldChange)? 2 : 3; # ]-Inf~1/dispFC] (1) or  ]1/dispFC-dispFC[ (2) or [dispFC~+Inf[ (3)
						$usedQuantifs2++;
						$goodQuantif2=1 if $okFilters==2;
					}
				}
				push @ratioList,$ratio;
				push @catList,$cat;
			}
			#>Printing data
			print DATA1 join("\t",($modProtID,@ratioList)),"\n" if ($goodQuantif && $usedQuantifs >= 2);
			print DATA2 join("\t",($modProtID,@catList)),"\n" if ($goodQuantif2 && $usedQuantifs2 >= 2);
			#next if $usedQuantifs2==0; # 1 quantif is enough (including +/-Inf)
			next if $goodQuantif2==0; # 1 quantif is enough (including +/-Inf)
			print DATA3 $modProtID;
			foreach my $idx (0..$#catList) {
				if ($catList[$idx] eq 'NA') {print DATA3 "\t0";}
				else {
					print DATA3 "\t1";
					#$protInQuantifs{$modProtID}{$selectedQuantifications[$idx]}=1;
					$protInQuantifs{$selectedQuantifications[$idx]}++;
				}
			}
			print DATA3 "\n";
		}
		close DATA2;
	}
	elsif ($selQuantifFamily eq 'MQ') { # MaxQuant States
		$jobInfo{1}={TITLE=>"$maxQuantMeasures{$dispMeasure} correlation",
					 DETAILS=>"&bull;Based only on values from ${proteoform}s shared by both quantifications."
					};
		PROT: foreach my $modProtID (keys %quantifValues) {
			#>Filtering data
			my $usedQuantifs=0;
			my $goodQuantif=0;
			my @valueList;
			foreach my $quantif (@selectedQuantifications) {
				my $value='NA';
				if (defined $quantifValues{$modProtID}{$quantif}) {
					my $okFilters=($quantifValues{$modProtID}{$quantif}{$numPepCode} >= $dispNumPep)? 1 : 0;
					if (!$strictFilter || $okFilters) {
						$value=log($quantifValues{$modProtID}{$quantif}{$dispMeasure})/$log10;
						$usedQuantifs++;
						$goodQuantif=1 if $okFilters;
					}
				}
				push @valueList,$value;
			}
			print DATA1 join("\t",($modProtID,@valueList)),"\n" if ($goodQuantif && $usedQuantifs >= 2);
			next if $usedQuantifs==0;
			print DATA3 $modProtID;
			foreach my $idx (0..$#valueList) {
				if ($valueList[$idx] eq 'NA') {print DATA3 "\t0";}
				else {
					print DATA3 "\t1";
					$protInQuantifs{$selectedQuantifications[$idx]}++;
				}
			}
			print DATA3 "\n";
		}
	}
	elsif ($selQuantifFamily =~ /EMPAI|SIN/) {
		my $qValueCode=($selQuantifFamily eq 'EMPAI')? $empaiValue : 'SIN_SIN';
		$jobInfo{1}={TITLE=>"$emPAIMeasures{$empaiValue} correlation",
					 DETAILS=>"&bull;Based only on values from ${proteoform}s shared by both quantifications."
					 };
		foreach my $modProtID (keys %quantifValues) {
			#>Filtering data
			my $goodQuantif=0;
			my @valueList;
			foreach my $quantif (@selectedQuantifications) {
				my $value='NA';
				if (defined $quantifValues{$modProtID}{$quantif}) {
					$value=log($quantifValues{$modProtID}{$quantif}{$qValueCode})/$log10;
					$goodQuantif++;
				}
				push @valueList,$value;
			}
			print DATA1 join("\t",($modProtID,@valueList)),"\n" if $goodQuantif >= 2;
			next if $goodQuantif==0;
			print DATA3 $modProtID;
			foreach my $idx (0..$#valueList) {
				if ($valueList[$idx] eq 'NA') {print DATA3 "\t0";}
				else {
					print DATA3 "\t1";
					$protInQuantifs{$selectedQuantifications[$idx]}++;
				}
			}
			print DATA3 "\n";
		}
	}
	close DATA1;
	close DATA3;

	##>Data error management
	unless (scalar keys %protInQuantifs) {
		rmtree($workDir);
		print qq
| Failed.</FONT>
</DIV>
<SCRIPT type="text/javascript">
document.getElementById('waitDIV').style.display='none';
</SCRIPT>
<BR><BR>
<FONT class="title2"><FONT style="color:#DD0000">ERROR:</FONT> No $proteoform data matched selected red filters.</FONT>
<BR><BR>
</BODY>
</HTML>
|;
		exit;
	}

	##>Shared proteins
	#foreach my $idx1 (0..$#selectedQuantifications) {
	#	foreach my $idx2 ($idx1..$#selectedQuantifications) {
	#		if ($idx1==$idx2) {
	#			$protUnion{$selectedQuantifications[$idx2]}{$selectedQuantifications[$idx1]}=1; # real value doesn't matter
	#			next;
	#		}
	#		foreach my $modProtID (keys %protInQuantifs) {
	#			$protUnion{$selectedQuantifications[$idx1]}{$selectedQuantifications[$idx2]}++ if ($protInQuantifs{$modProtID}{$selectedQuantifications[$idx1]} || $protInQuantifs{$modProtID}{$selectedQuantifications[$idx2]});
	#		}
	#		$protUnion{$selectedQuantifications[$idx2]}{$selectedQuantifications[$idx1]}=$protUnion{$selectedQuantifications[$idx1]}{$selectedQuantifications[$idx2]};
	#	}
	#}

	print '.';
	####>Running R<####
	my $bashFile="$workDir/command.sh";
	open (BASH,">$bashFile");
	print BASH qq
|#!/bin/bash

cd $workDir
$promsPath{R}/R CMD BATCH --no-save --no-restore \"--args correlation $dataFile1\" $promsPath{R_scripts}/correlation.R
|;
	print BASH "$promsPath{R}/R CMD BATCH --no-save --no-restore \"--args correlation $dataFile2\" $promsPath{R_scripts}/correlation.R\n" if $selQuantifFamily eq 'RATIO';
	print BASH qq
|$promsPath{R}/R CMD BATCH --no-save --no-restore \"--args distance $dataFile3\" $promsPath{R_scripts}/correlation.R

echo END > $workDir/end.txt
|;
	close BASH;
	chmod 0775, $bashFile;

	system "$bashFile &";

	my $nbWhile=0;
	my $maxNbWhile= 3 * ($numSelectedQuantifs * ($numSelectedQuantifs-1)) / 2; # 3*5sec per comparison
	while (!-e "$workDir/end.txt") {
		if ($nbWhile > $maxNbWhile) {
			#rmtree($workDir);
			print qq
| Failed.</FONT>
</DIV>
<SCRIPT type="text/javascript">
document.getElementById('waitDIV').style.display='none';
</SCRIPT>
<BR><BR>
<FONT class="title2"><FONT style="color:#DD0000">ERROR:</FONT> R is taking too long or died before completion.</FONT>
<BR><BR>
</BODY>
</HTML>
|;
			exit;
		}
		sleep 5; # 12 loops/min
		$nbWhile++;
		print '.';
	}

	##>R error management
	my $RoutStrg=`tail -3 $workDir/correlation.Rout`;
	unless ($RoutStrg=~/proc\.time\(\)/) {
		#rmtree($workDir);
		print qq
| Failed.</FONT>
</DIV>
<SCRIPT type="text/javascript">
document.getElementById('waitDIV').style.display='none';
</SCRIPT>
<BR><BR>
<FONT class="title2"><FONT style="color:#DD0000">ERROR:</FONT> R has generated an error.</FONT>
<BR><BR>
</BODY>
</HTML>
|;
		exit;
	}
	print '.';

	####>Parsing result files<####
	my %correlResults;
	foreach my $job (keys %correlFiles) {
		my @methods=($job==3)? ('Binary') : ('Pearson','Spearman'); #'Manhattan',
		foreach my $method (@methods) {
			open(RES,"$workDir/Results_$method"."_$correlFiles{$job}");
			my @results=<RES>;
			close RES;
			chomp($results[0]);
			$results[0]=~s/X//g; # remove X tag added by R
			my @quantifs=split(/\t/,$results[0]);
			foreach my $lineIdx (1..$#results) {
				chomp($results[$lineIdx]);
				my @data=split(/\t/,$results[$lineIdx]);
				#my $lineQuantif=shift(@data); # remove row name
				#$lineQuantif=~s/X//;
				my $lineQuantif=$quantifs[$lineIdx-1];
				foreach my $quantifIdx (0..$#data) {
					# 1-dist -> similarity
					#if ($method eq 'Manhattan') {$data[$quantifIdx]/=$protUnion{$lineQuantif}{ $quantifs[$quantifIdx] };} # "relative" Manhattan -> Soergel
					$correlResults{$job}{$method}{$lineQuantif}{ $quantifs[$quantifIdx] }=$data[$quantifIdx];
				}
			}
		}
	}

	####>Displaying results<####
	my %methodLinks=(
		'Pearson'=>'https://en.wikipedia.org/wiki/Pearson_correlation_coefficient',
		'Spearman'=>'https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient',
		'Binary'=>'https://en.wikipedia.org/wiki/Hamming_distance'
	);
	print qq
| Done.</FONT>
</DIV>
<SCRIPT type="text/javascript">
document.getElementById('waitDIV').style.display='none';
document.getElementById('correlProgressDIV').style.display='none';
</SCRIPT>
<BR>
<B>$filterInfoString.</B>
<BR>
<DIV> <!-- container for floating DIVs -->
|;

	my $span=$numSelectedQuantifs;
	my $footerSpan=$numSelectedQuantifs + 2;
	foreach my $job (sort{$a<=>$b} keys %correlResults) {
		my @methods=($job==3)? ('Binary') : ('Pearson','Spearman'); # 'Manhattan',
		print qq
|<DIV style="text-align:center; padding:5px 30px; float:left">
<FONT class="title3">$jobInfo{$job}{TITLE}</FONT>
<TABLE class="correlation" border=1>
<TABLE border=0 cellpadding=5>
<TR><TD rowspan=2 colspan=2></TD><TH class="colHeader" colspan=$span><A href="$methodLinks{$methods[0]}" target="_blank">$methods[0]</A></TH></TR>
<TR>
|;
		foreach my $quantif (@selectedQuantifications) {
			print "<TH class=\"colHeader\" nowrap><A href=\"javascript:void(null)\" onmouseover=\"popup('<B>$quantifAllNames{FULL}{$quantif}</B>')\" onmouseout=\"popout()\">",&promsMod::shortenName($quantifAllNames{OPTIMAL}{$quantif},25),"</A></TH>";
		}
		my $colSpan=($methods[1])? 1 : 2;
		foreach my $idx1 (0..$#selectedQuantifications) {
			print "<TR>";
			if ($idx1==0 && $methods[1]) {
				my $methodStrg=join('<BR>',(split(//,$methods[1])));
				print "<TH class=\"rowHeader\" rowSpan=\"$span\" align=\"center\"><A href=\"$methodLinks{$methods[1]}\" target=\"_blank\">$methodStrg</A></TH>\n";
			}
			my $numProtStrg=($job==3)? " ($protInQuantifs{$selectedQuantifications[$idx1]} ${proteoform}s)" : '';
			print "<TH class=\"rowHeader\" colspan=\"$colSpan\" nowrap><A href=\"javascript:void(null)\" onmouseover=\"popup('<B>$quantifAllNames{FULL}{$selectedQuantifications[$idx1]}</B>')\" onmouseout=\"popout()\">",&promsMod::shortenName($quantifAllNames{OPTIMAL}{$selectedQuantifications[$idx1]},25),"</A>$numProtStrg</TH>";
			foreach my $idx2 (0..$#selectedQuantifications) {
				my $method=($idx2 >= $idx1)? $methods[0] : $methods[1];
				if ($method) {
					my $correlColorS=($job==3)? sprintf '%.0f',100+155*$correlResults{$job}{$method}{ $selectedQuantifications[$idx1] }{ $selectedQuantifications[$idx2] }
											  : sprintf '%.0f',100+155*(1-$correlResults{$job}{$method}{ $selectedQuantifications[$idx1] }{ $selectedQuantifications[$idx2] });
					my $dispValue=1 * sprintf "%.3f",$correlResults{$job}{$method}{ $selectedQuantifications[$idx1] }{ $selectedQuantifications[$idx2] };
					print "<TH style=\"background-color:rgb(255,$correlColorS,$correlColorS)\"><A style=\"padding:30px\" href=\"javascript:void(null)\" onmouseover=\"popup('<B>$jobInfo{$job}{TITLE} ($method):</B><BR>$quantifAllNames{FULL}{$selectedQuantifications[$idx1]}<BR>vs.<BR>$quantifAllNames{FULL}{$selectedQuantifications[$idx2]}')\" onmouseout=\"popout()\">$dispValue</A></TH>";
				}
				else {print "<TH></TH>";}
			}
			print "</TR>\n";
		}
		print qq
|<TR><TD colspan="$footerSpan">$jobInfo{$job}{DETAILS}</TD></TR>
</TABLE>
</DIV>
|;
	}
	print "</DIV>\n"; # end of container

	rmtree($workDir);

}

######################
####>Protein List<####
######################
elsif ($view eq 'list') {
	$filterInfoString=($strictFilter)? "Strict filtering: For each $proteoform, only quantifications that match red filters are displayed" : "Loose filtering: For each $proteoform, all quantifications are displayed if at least 1 matches red filters";

	if ($exportType) {
		&exportProteinList(\@selectedQuantifications,\%quantifInfo,\%quantifValues,\%proteinInfo,\%quantifPTM);
		$workbook->close();
		exit;
	}

	print qq
|<SCRIPT type="text/javascript">
document.getElementById('waitDIV').style.display='none';
</SCRIPT>
<DIV style="float:top"><B>$filterInfoString.</B></DIV>
|;
	&printProteinList('COMPARE',\@selectedQuantifications,\%quantifInfo,\%quantifValues,\%proteinInfo,\%quantifAllNames,\%quantifPTM);
}

#########################
####>End of document<####
#########################
if ($view ne 'list' && $view ne 'correlMatrix') {
	print qq
|<DIV id="resultDIV" style="visibility:hidden">
<DIV style="float:top"><B>$filterInfoString.</B>
|;
	if ($numSelectedQuantifs <= 5) {
		print qq
|<DIV style="text-align:left;">
<BR>
<INPUT type="button" value="Export as image" style="font-weight:bold;font-size:10px" onclick="exportSVGtoImg('vennDIV','VennDiagram','./exportSVG.cgi')"/>
<DIV id="vennDIV" style="text-align:left;padding:3px;"></DIV>
</DIV>
|;
	}
	print qq
|</DIV>
<DIV id="mainGraphDIV" style="float:left;background-color:#FFFFFF;border:solid 3px $darkColor;padding:2px;"></DIV>
<DIV id="protListDIV" style="position:absolute;height:535;overflow:auto"></DIV>
</DIV>
|;
}
print qq
|<DIV id="displayDIV" class="popup"> <!--filter:alpha(opacity=80);opacity:0.8;-->
	<DIV id="infoDIV"></DIV>
	<DIV id="infoDIV2"></DIV>
</DIV>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup();</SCRIPT>
<DIV style="clear:both"></DIV>
<BR>
</BODY>
</HTML>
|;
exit;


#############################################<<<SUBROUTINES>>>###########################################

##########################################################################
#######################<Protein list routine>#############################
##########################################################################
sub printProteinList {
	# Globals: $dispFoldChange,$dispNumPep,$dispPvalue
	my ($call,$refSelQuantifs,$refQuantifInfo,$refQuantifValues,$refProteinInfo,$refQuantifAllNames,$refQuantifPTM)=@_; #,$dispStdDev
	my $numModifQuantifs=scalar keys %{$refQuantifPTM};
	my $numProtQuantifs=$numSelectedQuantifs-$numModifQuantifs;
	my $numTotIsoforms=0;
	if ($numModifQuantifs) {
		$numTotIsoforms=($numProtQuantifs)? scalar map {$_=~/^\d+-/? ($_) : ()} keys %{$refQuantifValues} : scalar keys %{$refQuantifValues};
	}
#my $numTotQuantifItems=scalar keys %{$refQuantifValues};
	my $numTotProteins=scalar keys %{$refProteinInfo};
	my ($minFoldChange,$maxFoldChange,$maxPvalue)=(0.5,2,0.05); # default values for arrow flag
	my $disabSave=($projectAccess eq 'guest')? ' disabled' : '';
	my (%quantifCode,%vennDiagram); #,%quantifLabel
	my $quantIdx=0;
	foreach my $quantif (@{$refSelQuantifs}) {
		$quantifCode{$quantif}=chr(65+$quantIdx++); # A,B,C,D,E
	}
	if ($numSelectedQuantifs <= 5 && !$ajax) {
		print qq
|<DIV style="text-align:center;">
<INPUT type="button" value="Export as image" style="font-weight:bold;font-size:10px" onclick="exportSVGtoImg('vennDIV','VennDiagram','./exportSVG.cgi')"/>
<DIV id="vennDIV" style="text-align:center;padding:3px;"></DIV>
</DIV>
|;
	}
	print qq
|<FORM name="protForm" method="POST">
<INPUT type="hidden" name="quantifFamily" value="$selQuantifFamily"/>
<INPUT type="hidden" name="ACT" value="$action"/>
<INPUT type="hidden" name="pepType" value="$dispPepType"/>
<INPUT type="hidden" name="id_project" value="$projectID"/>
<INPUT type="hidden" name="parentItem" value="$parentItem"/>
<INPUT type="hidden" name="view" value="$view"/>
<INPUT type="hidden" name="sort" value="$sortOrder"/>
<INPUT type="hidden" name="restrictList" value="$restrictListID"/>
<INPUT type="hidden" name="listAction" value="$listAction"/>
<INPUT type="hidden" name="quantifList" value="$quantifList"/>
<INPUT type="hidden" name="foldChange" value="$dispFoldChange"/>
<INPUT type="hidden" name="pValue" value="$dispPvalue"/>
<INPUT type="hidden" name="numPep" value="$dispNumPep"/>
<INPUT type="hidden" name="strictFilter" value="$strictFilter"/>
<INPUT type="hidden" name="exportType" value="all"/>
<INPUT type="hidden" name="empaiType" value="$empaiValue"/>


<TABLE border=0 cellspacing=0 cellpadding=2 align=center>
|;
	my $proteinText=($sortOrder eq 'protein')? '<FONT color=#DD0000>proteins</FONT>' : 'proteins';
	if ($ajax) {
		$proteinText="<A href=\"javascript:selectSort(\'protein\',$ajax)\" onmouseover=\"popup(\'Click to sort proteins by <B>ascending name</B>.\')\" onmouseout=\"popout()\">$proteinText</A>";
	}
	else {
		$proteinText="<A href=\"javascript:selectSort(\\'protein\\',$ajax)\" onmouseover=\"popup(\\'Click to sort proteins by <B>ascending name</B>.\\')\" onmouseout=\"popout()\">$proteinText</A>";
	}
	my $protString=($numModifQuantifs)? "&nbsp;$numTotIsoforms sites&nbsp;<BR>&nbsp;$numTotProteins $proteinText&nbsp;" : "&nbsp;$numTotProteins $proteinText&nbsp;";
	my ($numDispItems,$numDispIsoforms)=(0,0);
	my %dispProteins;

	if ($selQuantifFamily =~/RATIO|MQ/) {
		#my $invDispFoldChange=1/$dispFoldChange;
		my $clearButtonStrg=($view eq 'list' && !$ajax)? '' : '<INPUT type="button" value="Clear list" style="font-weight:bold" onclick="document.getElementById(\'protListDIV\').style.display=\'none\'">';
		#my ($ratioColspan,$stdDevColStrg)=($quantifMethod eq 'TNPQ')? (3,'') : (4,'<TH class="rbBorder" nowrap>&nbsp;Std. dev.&nbsp;</TH>'); # ratio/arrow/p-value/(+/-stdDev)
		my ($qValueCode,$qValueType,$qValueText)=($selQuantifFamily eq 'MQ')? ($dispMeasure,'state_',$maxQuantMeasures{$dispMeasure}) : ('RATIO','ratio_','Ratio');
		my $ratioColspan=3;
		my $numColumns=$numSelectedQuantifs*$ratioColspan+4;
		print qq
|<TR><TD colspan=$numColumns><DIV id="saveProtDIV" style="display:none;"></DIV></TD></TR>
<TR><TD colspan=$numColumns>$clearButtonStrg
<INPUT type="button" value="Check all" onclick="checkAllProteins(true)"/>
<INPUT type="button" value="Uncheck all" onclick="checkAllProteins(false)"/>
<INPUT type="button" id="saveFormBUTTON" value="Save proteins..." onclick="ajaxManageSaveProteins('getThemes','PROT')"$disabSave/>
|;
		#<Check if at least 1 modif quantif
		my $okModifQuantif=($selQuantifFocus > 0)? 1 : 0;
		if ($selQuantifFocus < 0) { # modif and/or prot quantifs
			foreach my $quantif (@{$refSelQuantifs}) {
				my ($quantifID,$ratioPos)=split('_',$quantif);
				if ($refQuantifInfo->{$quantifID}[4]) { # modif quantif
					$okModifQuantif=1;
					last;
				}
			}
		}
		print "<INPUT type=\"button\" id=\"saveSiteFormBUTTON\" value=\"Save sites...\" onclick=\"ajaxManageSaveProteins('getThemes','SITE',$selModifID)\"$disabSave/>\n" if $okModifQuantif;
		print "<INPUT type=\"button\" value=\"Export selected proteins\" onclick=\"exportProteins('selected')\"/>\n<INPUT type=\"button\" value=\"Export all proteins\" onclick=\"exportProteins('all')\"/></TD>\n" unless $ajax;
		print qq
|</TR>
<TR class="row_0">
<TH class="rbBorder" align=left rowspan=2 nowrap><DIV id="protCountDIV">$protString</DIV></TH> <!--$numTotProteins <A href="javascript:selectSort('protein',$ajax)" onmouseover="popup('Click to sort proteins by <B>ascending name</B>.')" onmouseout="popout()">$protString</A>&nbsp;</TH> -->
<TH class="rbBorder" rowspan=2 nowrap>&nbsp;Gene&nbsp;</TH>
|;
		my (%quantifAnalyses,%processedQuantifID);
		foreach my $quantif (@{$refSelQuantifs}) {
			my ($quantifID,$ratioPos)=split('_',$quantif);
			#$quantifAnalyses{$quantif}=$refQuantifInfo->{$quantifID}[5]; # anaID Strg 'ID1,ID2,...'
			unless ($processedQuantifID{$quantifID}) {
				foreach my $anaID (split(',',$refQuantifInfo->{$quantifID}[5])) { # anaID Strg 'ID1,ID2,...'
					$quantifAnalyses{$anaID}=1;
				}
				$processedQuantifID{$quantifID}=1;
			}
			my $parentItemStrg;
			for (my $i=1; $i<$#{$refQuantifInfo->{$quantifID}[3]}; $i++) { # skip project & quantif
				$parentItemStrg.=' > ' if $parentItemStrg;
				$parentItemStrg.=$refQuantifInfo->{$quantifID}[3][$i]->{'NAME'};
			}
			#my ($testStatePos,$refStatePos)=split(/\//,$refQuantifInfo->{$quantifID}[1]->{'RATIOS'}[$ratioPos-1]);
			#my $ratioTag='';
			#if ($testStatePos=~/%/) { # 2ndary ratio in SuperRatio
			#	$testStatePos=~s/%\d+//;
			#	$refStatePos=~s/%\d+//;
			#	$ratioTag=encode_utf8('°');
			#}
			#my $quantifFullName=$refQuantifInfo->{$quantifID}[0].":".$refQuantifInfo->{$quantifID}[2]->{$testStatePos}{'NAME'}.$ratioTag."/".$refQuantifInfo->{$quantifID}[2]->{$refStatePos}{'NAME'}.$ratioTag;
			my $quantifName=&promsMod::shortenName($refQuantifAllNames->{OPTIMAL}{$quantif},22); # $quantifFullName
			#$quantifLabel{$quantif}=$quantifName;
			print "<TH class=\"rbBorder\" colspan=$ratioColspan nowrap>&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<B>$parentItemStrg<BR>$refQuantifAllNames->{FULL}{$quantif}</B>')\" onmouseout=\"popout()\">$quantifName</A>&nbsp;</TH>\n";
		}
		my $massText=($sortOrder eq 'mass')? '<FONT color=#DD0000>MW<SMALL> kDa</SMALL></FONT>' : 'MW<SMALL> kDa</SMALL>';
		print qq
|<TH class="rbBorder" rowspan=2 nowrap>&nbsp;<A href="javascript:selectSort('mass',$ajax)" onmouseover="popup('Click to sort proteins by <B>decreasing mass</B>.')" onmouseout="popout()">$massText</A>&nbsp;</TH><TH class="bBorder" rowspan=2 width=700 nowrap>&nbsp;Description - Species&nbsp;</TH>
</TR>
<TR class="row_0">
|;
		foreach my $quantif (@{$refSelQuantifs}) {
			my $qValueLabel=($sortOrder eq $qValueType.$quantif)? "<FONT color=#DD0000>$qValueText</FONT>" : $qValueText;
			my $peptideLabel=($sortOrder eq "peptide_$quantif")? '<FONT color=#DD0000>Pep. used</FONT>' : 'Pep. used';
			print qq
|<TH class="rbBorder" colspan=2 nowrap>&nbsp;<A href="javascript:selectSort('$qValueType$quantif',$ajax)" onmouseover="popup('Click to sort proteins by <B>descending $qValueText</B> in this quantification.')" onmouseout="popout()">$qValueLabel</A>&nbsp;</TH>
<TH class="rbBorder" nowrap>&nbsp;<A href="javascript:selectSort('peptide_$quantif',$ajax)" onmouseover="popup('Click to sort proteins by <B>descending peptides</B> in this quantification.')" onmouseout="popout()">$peptideLabel</A>&nbsp;</TH>
|; # <TH class=\"rbBorder\" nowrap>&nbsp;p-value&nbsp;</TH>$stdDevColStrg
		}
		print "</TR>\n";
		my %modifiedProteins;
		&getModifiedProteinForms($refQuantifValues,\%modifiedProteins);
		my $anaIDstrg=join(',',sort{$a<=>$b} keys %quantifAnalyses);
		foreach my $modProtID (sort{&sortProteins($sortOrder,$refProteinInfo,$refQuantifValues)} keys %{$refQuantifValues}) { #$refProteinInfo
			my ($protID,$modStrg)=($modProtID=~/^(\d+)(.*)/);
			unless ($modStrg) {
				next if $modifiedProteins{$protID}; # this unmodified form co-exists with modified forms => skip because displayed before/later with modified form
				                                    # for this to work ok with strict filtering, filtering must be applied by &promsQuantif::fetchQuantificationData ($parameters{STRICT_FILTER})
				$modStrg='';
			}
			my (@protInQuantif); ## %anaIDlist,
			my $quantifDataStrg='';
			my $okFilters=($ajax)? 1 : 0;
			my $vennStrg='';
			my $vennInputStrg='';
			foreach my $quantif (@{$refSelQuantifs}) {
				my $usedProtID;
				if ($modStrg && $refQuantifPTM->{$quantif}) { # modif vs modif
					if (!$refQuantifValues->{$modProtID}{$quantif} || !$refQuantifValues->{$modProtID}{$quantif}{$qValueCode}) {
						$quantifDataStrg.="<TH colspan=2>-</TH><TH>-</TH>"; #<TH>-</TH>
						#$quantifDataStrg.="<TH>-</TH>" if $quantifMethod ne 'TNPQ'; # Std dev
						next;
					}
					$usedProtID=$modProtID;
				}
				else { # fall back to whole prot
					if (!$refQuantifValues->{$protID}{$quantif} || !$refQuantifValues->{$protID}{$quantif}{$qValueCode}) {
						$quantifDataStrg.="<TH colspan=2>-</TH><TH>-</TH>"; #<TH>-</TH>
						#$quantifDataStrg.="<TH>-</TH>" if $quantifMethod ne 'TNPQ'; # Std dev
						next;
					}
					$usedProtID=$protID;
				}
				##$anaIDlist{$quantifAnalyses{$quantif}}=1 unless $anaIDstrg;

				if ($selQuantifFamily eq 'RATIO') {
					my $ratio=$refQuantifValues->{$usedProtID}{$quantif}{'RATIO'};
					my $quantifClass=($ratio <= $invDispFoldChange || $ratio >= $dispFoldChange)? 'TH' : 'TD';
					$quantifClass=($quantifClass eq 'TH' && ($dispPvalue==1 || defined $refQuantifValues->{$usedProtID}{$quantif}{'P_VALUE'} && $refQuantifValues->{$usedProtID}{$quantif}{'P_VALUE'} <= $dispPvalue))? 'TH' : 'TD';
					$quantifClass=($quantifClass eq 'TH' && defined $refQuantifValues->{$usedProtID}{$quantif}{$numPepCode} && $refQuantifValues->{$usedProtID}{$quantif}{$numPepCode} >= $dispNumPep)? 'TH' : 'TD';
					$okFilters=(!$okFilters && $quantifClass eq 'TH')? 1 : $okFilters;
					my ($absRatio,$arrowDir)=($ratio >= 1)? ($ratio,'up_') : (1/$ratio,'down_');
					my $ratioStrg=($absRatio >= $MAX_INF_RATIO_DB)? '&infin;' : sprintf "%.2f",$absRatio; #'∞'
					$ratioStrg='1/'.$ratioStrg if ($ratio < 1 && $ratioStrg ne '1.00');
					# Arrow flag
					my $arrowStrg='';
					my $okRatio=($ratio <= $minFoldChange || $ratio >= $maxFoldChange)? 1 : 0;
					my $okPvalue=(defined $refQuantifValues->{$usedProtID}{$quantif}{'P_VALUE'} && $refQuantifValues->{$usedProtID}{$quantif}{'P_VALUE'} <= $maxPvalue)? 1 : 0;
					if ($okRatio || $okPvalue) {
						my $arrowColor=(!$okRatio || !$okPvalue)? 'gray' : ($ratio >= 1)? 'red' : 'green';
						$arrowStrg=(defined $refQuantifValues->{$usedProtID}{$quantif}{'P_VALUE'})? "<IMG class=\"LINK\" src=\"$promsPath{images}/$arrowDir$arrowColor.png\" onclick=\"ajaxProteinData(event,'$usedProtID','$quantif','ajaxProtStat')\">" : "<IMG src=\"$promsPath{images}/$arrowDir$arrowColor.png\">";
					}
					if ($strictFilter && !$ajax && $quantifClass eq 'TD') { # this quantif did not pass filter
						$quantifDataStrg.="<TH colspan=2>-</TH><TH>-</TH>";
					}
					else {
						if ($absRatio < $MAX_INF_RATIO_DB) {
							$quantifDataStrg.="<TD class=\"$quantifClass\" nowrap align=right>&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$usedProtID','$quantif','ajaxProtStat')\"/>$ratioStrg</A></TD>"; #<TH nowrap>&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$quantif','ajaxProtStat')\">$pValueStrg</A>&nbsp;</TH>";
						}
						else {$quantifDataStrg.="<TD class=\"$quantifClass\" align=right nowrap>&nbsp;$ratioStrg&nbsp;</TD>";}
						$quantifDataStrg.="<TD>$arrowStrg&nbsp;</TD>";
						$quantifDataStrg.="<TD class=\"$quantifClass\" align=\"center\">&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$usedProtID','$quantif','ajaxPepData')\"/>$refQuantifValues->{$usedProtID}{$quantif}{$numPepCode}</A>&nbsp;</TD>";
						$vennStrg.=$quantifCode{$quantif} if (!$ajax && $numSelectedQuantifs <= 5);
					}
				}
				else { # MaxQuant States
					my $quantifClass=(defined $refQuantifValues->{$usedProtID}{$quantif}{$numPepCode} && $refQuantifValues->{$usedProtID}{$quantif}{$numPepCode} >= $dispNumPep)? 'TH' : 'TD';
					if ($strictFilter && !$ajax && $quantifClass eq 'TD') { # this quantif did not pass filter
						$quantifDataStrg.="<TH colspan=2>-</TH><TH>-</TH>";
					}
					else {
						$quantifDataStrg.="<TD colspan=2 class=\"$quantifClass\" align=right nowrap>&nbsp;$refQuantifValues->{$usedProtID}{$quantif}{$qValueCode}&nbsp;</TD><TD class=\"$quantifClass\" align=\"center\">&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$usedProtID','$quantif','ajaxPepData')\"/>$refQuantifValues->{$usedProtID}{$quantif}{$numPepCode}</A>&nbsp;</TD>";
						$vennStrg.=$quantifCode{$quantif} if (!$ajax && $numSelectedQuantifs <= 5);
					}
					$okFilters=(!$okFilters && $quantifClass eq 'TH')? 1 : $okFilters;
				}
			}
			next unless $okFilters;
			$numDispItems++;
			$numDispIsoforms++ if $modStrg;
			$dispProteins{$protID}=1;
			if (!$ajax && $numSelectedQuantifs <= 5) {
				$vennDiagram{$vennStrg}++;
				$vennInputStrg="<INPUT type=\"hidden\" name=\"vennInput\" value=\"$vennStrg\"/>";
			}
			##$anaIDstrg=join(',',keys %anaIDlist) unless $anaIDstrg;
			my $trClass='list row_'.($numDispItems % 2);
			my $mass=(!$refProteinInfo->{$protID}[1])? '-' : $refProteinInfo->{$protID}[1];
			print qq
|<TR class="$trClass" valign=top id="tr_$modProtID"><TD class="TH" nowrap><INPUT type="checkbox" name="chkProt" value="$modProtID"/>$vennInputStrg<A href="javascript:sequenceView($protID,'$anaIDstrg')">$refProteinInfo->{$protID}[0]$modStrg</A>&nbsp;</TD>
<TH>&nbsp;|;
			if (scalar @{$refProteinInfo->{$protID}[5]} > 1) { # gene
				print "<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-",join('<BR>&nbsp;&nbsp;-',@{$refProteinInfo->{$protID}[5]}[1..$#{$refProteinInfo->{$protID}[5]}]),"</B>')\" onmouseout=\"popout()\">",$refProteinInfo->{$protID}[5][0],"</A>";
			}
			elsif ($refProteinInfo->{$protID}[5][0]) {print $refProteinInfo->{$protID}[5][0];}
			else {print '-';} # no gene mapped
			print qq
|&nbsp;</TH>
$quantifDataStrg
<TD class="TH" align=right>$mass&nbsp;</TD><TD>$refProteinInfo->{$protID}[2] <U><I>$refProteinInfo->{$protID}[3]</I></U></TD>
</TR>
|;
		}
		print "<TR><TD colspan=$numColumns><B>End of list.</B>";
		print qq|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<IMG src="$promsPath{images}/up_red.png"/><IMG src="$promsPath{images}/down_green.png"/>: Fold change &ge; 2 <B>and</B> p-value &le; 0.05&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<IMG src="$promsPath{images}/up_gray.png"/><IMG src="$promsPath{images}/down_gray.png"/>: Fold change &ge; 2 <B>or</B> p-value &le; 0.05| if $selQuantifFamily eq 'RATIO';
		print "</TD></TR>\n";
	}
	elsif ($selQuantifFamily =~ /EMPAI|SIN/) {
		#my $invDispFoldChange=1/$dispFoldChange;
		my $clearButtonStrg=($view eq 'list' && !$ajax)? '' : '<INPUT type="button" value="Clear list" style="font-weight:bold" onclick="document.getElementById(\'protListDIV\').style.display=\'none\'">';
		#my ($ratioColspan,$stdDevColStrg)=($quantifMethod eq 'TNPQ')? (3,'') : (4,'<TH class="rbBorder" nowrap>&nbsp;Std. dev.&nbsp;</TH>'); # ratio/arrow/p-value/(+/-stdDev)
		my $specColspan=2; # For emPAI & label free MaxQuant methods
		my $numColumns=$numSelectedQuantifs*$specColspan+4;
		print qq
|<TR><TD colspan=$numColumns><DIV id="saveProtDIV" style="display:none;"></DIV></TD></TR>
<TR><TD colspan=$numColumns>$clearButtonStrg<INPUT type="button" value="Check all" onclick="checkAllProteins(true)"/><INPUT type="button" value="Uncheck all" onclick="checkAllProteins(false)"/><INPUT type="button" id="saveFormBUTTON" value="Save proteins..." onclick="ajaxManageSaveProteins('getThemes')"$disabSave/>
|;
		print "&nbsp;&nbsp;<INPUT type=\"button\" value=\"Export selected proteins\" onclick=\"exportProteins('selected')\"/><INPUT type=\"button\" value=\"Export all proteins\" onclick=\"exportProteins('all')\"/></TD>\n" unless $ajax;
		print qq
|</TR>
<TR class="row_0">
<TH class="rbBorder" align=left rowspan=2 nowrap><DIV id="protCountDIV">$protString</DIV></TH> <!--$numTotProteins <A href="javascript:selectSort('protein',$ajax)" onmouseover="popup('Click to sort proteins by <B>ascending name</B>.')" onmouseout="popout()">$protString</A>&nbsp;</TH> -->
<TH class="rbBorder" rowspan=2 nowrap>&nbsp;Gene&nbsp;</TH>
|;
		my %quantifAnalyses;
		foreach my $quantif (@{$refSelQuantifs}) {
			#$quantifAnalyses{$quantif}=$refQuantifInfo->{$quantif}[4]; # anaID
			$quantifAnalyses{$refQuantifInfo->{$quantif}[5]}=1; # anaID
			my $parentItemStrg;
			for (my $i=1; $i<$#{$refQuantifInfo->{$quantif}[3]}; $i++) { # skip project & quantif
				$parentItemStrg.=' > ' if $parentItemStrg;
				$parentItemStrg.=$refQuantifInfo->{$quantif}[3][$i]->{'NAME'};
			}
			my $quantifName=&promsMod::shortenName($refQuantifInfo->{$quantif}[0],22); # $quantifName
			#print "<TH class=\"rbBorder\" colspan=$ratioColspan nowrap>&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<B>$parentItemStrg<BR>$refQuantifAllNames->{FULL}{$quantif}</B>')\" onmouseout=\"popout()\">$quantifName</A>&nbsp;</TH>\n";
			print "<TH class=\"rbBorder\" colspan=$specColspan nowrap>&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<B>$parentItemStrg<BR>$quantifName</B>')\" onmouseout=\"popout()\">$quantifName</A>&nbsp;</TH>\n";
		}
		my $massText=($sortOrder eq 'mass')? '<FONT color=#DD0000>MW<SMALL> kDa</SMALL></FONT>' : 'MW<SMALL> kDa</SMALL>';
		print qq
|<TH class="rbBorder" rowspan=2 nowrap>&nbsp;<A href="javascript:selectSort('mass',$ajax)" onmouseover="popup('Click to sort proteins by <B>decreasing mass</B>.')" onmouseout="popout()">$massText</A>&nbsp;</TH><TH class="bBorder" rowspan=2 width=700 nowrap>&nbsp;Description - Species&nbsp;</TH>
</TR>
<TR class="row_0">
|;
		my ($labelQCol,$labelQVal)=($selQuantifFamily eq 'EMPAI')? ($emPAIMeasures{$empaiValue},$empaiValue) : ('SI<sub>N</sub>','SIN_SIN');
		my $qFamlc=lc($selQuantifFamily);
		foreach my $quantif (@{$refSelQuantifs}) {
			my $colLabel=($sortOrder eq "${qFamlc}_$quantif")? "<FONT color=#DD0000>$labelQCol</FONT>" : $labelQCol;
			my $peptideLabel=($sortOrder eq "peptide_$quantif")? '<FONT color=#DD0000>Pep. found</FONT>' : 'Pep. found';
			print qq
|<TH class="rbBorder" nowrap>&nbsp;<A href="javascript:selectSort('${qFamlc}_$quantif',$ajax)" onmouseover="popup('Click to sort proteins by <B>descending value</B> in this quantification.')" onmouseout="popout()">$colLabel</A>&nbsp;</TH>
<TH class="rbBorder" nowrap>&nbsp;$peptideLabel&nbsp;</TH>
|; # <TH class=\"rbBorder\" nowrap>&nbsp;p-value&nbsp;</TH>$stdDevColStrg
		}
		print "</TR>\n";
		my $anaIDstrg=join(',',sort{$a<=>$b} keys %quantifAnalyses);
		foreach my $protID (sort{&sortProteins($sortOrder,$refProteinInfo,$refQuantifValues)} keys %{$refQuantifValues}) { #$refProteinInfo
			next if $refProteinInfo->{$protID}[0] =~/^CON__/;
			my $okFilters=($ajax)? 1 : 0;
			my $vennStrg='';
			my $vennInputStrg='';
			my (@protInQuantif); ## %anaIDlist,
			my $quantifDataStrg='';
			foreach my $quantif (@{$refSelQuantifs}) {
				if (!$refQuantifValues->{$protID}{$quantif} || !defined($refQuantifValues->{$protID}{$quantif}{$labelQVal})) { # defined because can be 0!
					$quantifDataStrg.="<TH>-</TH><TH>-</TH>"; #<TH>-</TH>
					next;
				}
				##$anaIDlist{$quantifAnalyses{$quantif}}=1 unless $anaIDstrg;
				my $quantifClass='TD';
				my $trClass='list row_'.($numDispItems % 2);
				if ($strictFilter && !$ajax && $quantifClass eq 'TD') { # this quantif did not pass filter
					$quantifDataStrg.="<TH>-</TH><TH>-</TH>";
				}
				else {
					my $qVal=$refQuantifValues->{$protID}{$quantif}{$labelQVal};
					$qVal=($qVal >= 0.1 && $qVal < 1000)? sprintf '%.2f',$qVal : ($qVal==0)? $qVal : sprintf '%.2e',$qVal;
					#my $numPep=($selQuantifFamily eq 'EMPAI')? $pepInfo{$quantif}{$protID} : $refQuantifValues->{$protID}{$quantif}{'NUM_PEP_TOTAL'};
					my $numPep=$refQuantifValues->{$protID}{$quantif}{'PEPTIDES'} || "-";
					$quantifDataStrg.="<TD class=\"$quantifClass\" align=right nowrap>&nbsp;$qVal&nbsp;</TD>";
					$quantifDataStrg.="<TD class=\"$quantifClass\" align=\"center\">&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$protID','$quantif','ajaxPepData')\"/>$numPep</A>&nbsp;</TD>";
					$vennStrg.=$quantifCode{$quantif} if (!$ajax && $numSelectedQuantifs <= 5);
				}
			}
			#next unless $okFilters;
			$numDispItems++;
			$dispProteins{$protID}=1;
			if (!$ajax && $numSelectedQuantifs <= 5) {
				$vennDiagram{$vennStrg}++;
				$vennInputStrg="<INPUT type=\"hidden\" name=\"vennInput\" value=\"$vennStrg\"/>";
			}
			##$anaIDstrg=join(',',keys %anaIDlist) unless $anaIDstrg;
			my $trClass='list row_'.($numDispItems % 2);
			my $mass=(!$refProteinInfo->{$protID}[1])? '-' : $refProteinInfo->{$protID}[1];
			print qq
|<TR class="$trClass" valign=top id="tr_$protID"><TD class="TH" nowrap><INPUT type="checkbox" name="chkProt" value="$protID"/>$vennInputStrg<A href="javascript:sequenceView($protID,'$anaIDstrg')">$refProteinInfo->{$protID}[0]</A>&nbsp;</TD>
<TH>&nbsp;|;
			if (scalar @{$refProteinInfo->{$protID}[5]} > 1) { # gene
				print "<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-",join('<BR>&nbsp;&nbsp;-',@{$refProteinInfo->{$protID}[5]}[1..$#{$refProteinInfo->{$protID}[5]}]),"</B>')\" onmouseout=\"popout()\">",$refProteinInfo->{$protID}[5][0],"</A>";
			}
			elsif ($refProteinInfo->{$protID}[5][0]) {print $refProteinInfo->{$protID}[5][0];}
			else {print '-';} # no gene mapped
			print qq
|&nbsp;</TH>
$quantifDataStrg
<TD class="TH" align=right>$mass&nbsp;</TD><TD>$refProteinInfo->{$protID}[2] <U><I>$refProteinInfo->{$protID}[3]</I></U></TD>
</TR>
|;
		}
		print qq
|<TR><TD colspan=$numColumns><B>End of list.</B></TD></TR>
|;
	}
	print qq
|</TABLE>
</FORM>
|;
	unless ($ajax) { # not possible if ajax call
		#my $protTypeStrg=(scalar keys %{$refQuantifPTM})? ' sites' : '';
		$protString=($numModifQuantifs)? "&nbsp;$numDispIsoforms/$numTotIsoforms sites<BR>&nbsp;".(scalar keys %dispProteins)."/$numTotProteins $proteinText&nbsp;" : "&nbsp;$numDispItems/$numTotProteins $proteinText&nbsp;";
		print qq
|<SCRIPT LANGUAGE="javascript">
document.getElementById('protCountDIV').innerHTML='$protString';
|;
		if ($numSelectedQuantifs <= 5) {
			print "window.onload=function() {\n\tgroupLabelsObj={";
			my $count=0;
			my $legendMaxSize=0;
			foreach my $quantif (@{$refSelQuantifs}) {
				$count++;
				#print "$quantifCode{$quantif}:'$quantifLabel{$quantif}'";
				my ($optimalName)=($refQuantifAllNames->{OPTIMAL}{$quantif})?$refQuantifAllNames->{OPTIMAL}{$quantif} : $refQuantifInfo->{$quantif}[0];
				print "$quantifCode{$quantif}:'$optimalName'";
				print ',' if $count < $numSelectedQuantifs;
				if ($legendMaxSize < length($optimalName)) {$legendMaxSize=length($optimalName);}
			}
			my $legendPos=($legendMaxSize < 50)? 'side' : 'bottom';
			print "};\n\tnew VennDiagram({div:'vennDIV',legendPosition:'$legendPos',size:300,clickAction:updateProteinList,groupLabels:groupLabelsObj,setValues:{";
			$count=0;
			my $numSets=scalar keys %vennDiagram;
			foreach my $set (sort keys %vennDiagram) {
				$count++;
				print "$set:$vennDiagram{$set}";
				print ',' if $count < $numSets;
			}
			print "}});\n}\n";
		}
		print qq
|function updateProteinList(set) {
	var myForm=document.protForm;
	var vennData=myForm.vennInput;
	var chkData=myForm.chkProt;
	// Scanning all proteins
	var existIsoforms=$numModifQuantifs;
	var protList={};
	var numMatch=0;
	var numProt=0;
	if (vennData.length) {
		for (var i=0; i<vennData.length; i++) {
			var chkStatus,displayStatus;
			if (set=='all') {numMatch++; chkStatus=true; displayStatus='';}
			else if (set.match(':full')) {
				var setInfo=set.split(':');
				var expr=new RegExp(setInfo[0]);
				if (vennData[i].value.match(expr)) {numMatch++; chkStatus=true; displayStatus='';}
				else {chkStatus=false; displayStatus='none';}
			}
			else {
				if (vennData[i].value==set) {numMatch++; chkStatus=true; displayStatus='';}
				else {chkStatus=false; displayStatus='none';}
			}
			chkData[i].checked=chkStatus;
			var tr=document.getElementById('tr_'+chkData[i].value);
			if (chkStatus) {
				tr.className='list row_'+(numMatch % 2);
				if (existIsoforms>=1) {
					protIdData=chkData[i].value.split('-');
					protList[protIdData[0]]=1;
				}
			}
			tr.style.display=displayStatus;
		}
		if (existIsoforms>=1) {
			for (var p in protList) {numProt++;}
		}
	}
	else {
		numMatch=1;
		numProt=1;
	}
	document.getElementById('protCountDIV').innerHTML=(existIsoforms>=1)? '&nbsp;'+numMatch+'/$numTotIsoforms sites<BR>&nbsp;'+numProt+'/$numTotProteins $proteinText&nbsp;' : '&nbsp;'+numMatch+'/$numTotProteins $proteinText&nbsp;';
}
</SCRIPT>
|;
	}
} # end of &printProteinList

#########################################################################
#######################<Export list routine>#############################
#########################################################################
sub exportProteinList { # Only for design or no-design labeled quantifs
	my ($refSelQuantifs,$refQuantifInfo,$refQuantifValues,$refProteinInfo,$refQuantifPTM)=@_; #,$dispStdDev
	#my $numTotProteins=scalar keys %{$refProteinInfo};
	#my $peptideTitle=($labelType eq 'SILAC')? 'Pept. sets' : 'Peptides'; # ??????????????

	####<Start printing>####
	my $worksheet2=$workbook->add_worksheet('Results');
	my $xlsRow=0;
	my $xlsCol=0;
	$worksheet2->set_column(0,0,30); # identifier col length

	if ($selQuantifFamily eq 'RATIO') {
		my $paramStrg="Parameters:";
		$paramStrg.=" •selected proteins only," if $exportType eq 'selected';
		$paramStrg.=" •abs. fold change ≥ $dispFoldChange, •p-value ≤ $dispPvalue, •$dispPepType peptides ≥ $dispNumPep, •$filterInfoString";
		$worksheet2->merge_range(0,0,0,4 + scalar @{$refSelQuantifs} * 4,$paramStrg,$itemFormat{'mergeColHeader'});
		$xlsRow++;
		# identifier header written at the end
		#<Gene
		$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow+1,$xlsCol,'Gene & Synonyms',$itemFormat{'mergeRowHeader'});
		#my $invDispFoldChange=1/$dispFoldChange;
		my $log2=log(2);
		foreach my $quantif (@{$refSelQuantifs}) {
			my ($quantifID,$ratioPos)=split('_',$quantif);
			my $parentItemStrg;
			for (my $i=1; $i<$#{$refQuantifInfo->{$quantifID}[3]}; $i++) { # skip project & quantif
				$parentItemStrg.=' > ' if $parentItemStrg;
				$parentItemStrg.=$refQuantifInfo->{$quantifID}[3][$i]->{'NAME'};
			}
			my ($testStatePos,$refStatePos)=split(/\//,$refQuantifInfo->{$quantifID}[1]->{'RATIOS'}[$ratioPos-1]);
			my $ratioTag='';
			if ($testStatePos=~/%/) { # 2ndary ratio in SuperRatio
				$testStatePos=~s/%\d+//;
				$refStatePos=~s/%\d+//;
				$ratioTag='°';
			}
			my $quantifName=$refQuantifInfo->{$quantifID}[0].":".$refQuantifInfo->{$quantifID}[2]->{$testStatePos}{'NAME'}.$ratioTag."/".$refQuantifInfo->{$quantifID}[2]->{$refStatePos}{'NAME'}.$ratioTag;
			$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow,$xlsCol+3,$parentItemStrg.' > '.$quantifName,$itemFormat{'mergeColHeader'});
			$worksheet2->write_string($xlsRow+1,$xlsCol,'Ratio',$itemFormat{'header'});
			$worksheet2->write_comment($xlsRow+1,$xlsCol,'Ratio of 1000 or 0.001 means protein was not detected in one condition.');
			$worksheet2->write_string($xlsRow+1,++$xlsCol,'Log2',$itemFormat{'header'});
			$worksheet2->write_string($xlsRow+1,++$xlsCol,'p-value',$itemFormat{'header'});
			$worksheet2->write_string($xlsRow+1,++$xlsCol,'Peptides used',$itemFormat{'header'});
		}
		#<MW, description & species
		$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow+1,$xlsCol,'MW (kDa)',$itemFormat{'mergeRowHeader'});
		$worksheet2->set_column(++$xlsCol,$xlsCol,80); # col length
		$worksheet2->merge_range($xlsRow,$xlsCol,$xlsRow+1,$xlsCol,'Description',$itemFormat{'mergeRowHeader'});
		$worksheet2->set_column(++$xlsCol,$xlsCol,30); # col length
		$worksheet2->merge_range($xlsRow,$xlsCol,$xlsRow+1,$xlsCol,'Species',$itemFormat{'mergeRowHeader'});

		###<Looping through proteins>###
		my %modifiedProteins;
		&getModifiedProteinForms($refQuantifValues,\%modifiedProteins);
		$xlsRow++;
		#my $numProt=0;
		my $numDispItems=0;
		my %dispProteins;
		foreach my $modProtID (sort{&sortProteins($sortOrder,$refProteinInfo,$refQuantifValues)} keys %{$refQuantifValues}) { #,\%pepInfo $refProteinInfo
			my ($protID,$modStrg)=($modProtID=~/^(\d+)(.*)/);
			unless ($modStrg) {
				next if $modifiedProteins{$protID}; #  unmodified prot & modified prot(s) exist(s) => skip because displayed before/later with modified prot
				$modStrg='';
			}
			$xlsCol=0;
			my $okFilters=0;
			my %okQuantif;
			foreach my $quantif (@{$refSelQuantifs}) {
				my $usedProtID;
				if ($modStrg && $refQuantifPTM->{$quantif}) { # modif vs modif
					next if (!$refQuantifValues->{$modProtID}{$quantif} || !$refQuantifValues->{$modProtID}{$quantif}{'RATIO'});
					$usedProtID=$modProtID;
				}
				else { # fall back to whole prot
					next if (!$refQuantifValues->{$protID}{$quantif} || !$refQuantifValues->{$protID}{$quantif}{'RATIO'});
					$usedProtID=$protID;
				}
				#<Ratio/p-value/num pept. filters
				unless ($okFilters) { # no ratioPos has passed filters
					$okQuantif{$quantif}=(($refQuantifValues->{$usedProtID}{$quantif}{'RATIO'} < 1 && $refQuantifValues->{$usedProtID}{$quantif}{'RATIO'} <= $invDispFoldChange) || ($refQuantifValues->{$usedProtID}{$quantif}{'RATIO'} >= 1 && $refQuantifValues->{$usedProtID}{$quantif}{'RATIO'} >= $dispFoldChange))? 1 : 0;
					$okQuantif{$quantif}=($okQuantif{$quantif} && ($dispPvalue>=1 || (defined $refQuantifValues->{$usedProtID}{$quantif}{'P_VALUE'} && $refQuantifValues->{$usedProtID}{$quantif}{'P_VALUE'} <= $dispPvalue)))? 1 : 0;
					$okQuantif{$quantif}=($okQuantif{$quantif} && $refQuantifValues->{$usedProtID}{$quantif}{$numPepCode} >= $dispNumPep)? 1 : 0;
				}
				$okFilters=1 if $okQuantif{$quantif};
			}
			next unless $okFilters;
			#$numProt++;
			$numDispItems++;
			$dispProteins{$protID}=1;
			#<Identifier
			$worksheet2->write_string(++$xlsRow,$xlsCol,$refProteinInfo->{$protID}[0].$modStrg,$itemFormat{'text'});
			#<Gene
			$worksheet2->write_string($xlsRow,++$xlsCol,join(',',@{$refProteinInfo->{$protID}[5]}),$itemFormat{'text'});
			#<Quantif
			foreach my $quantif (@{$refSelQuantifs}) {
				if (!$refQuantifValues->{$modProtID}{$quantif} || !$refQuantifValues->{$modProtID}{$quantif}{'RATIO'} || ($strictFilter && !$okQuantif{$quantif})) {
					foreach my $i (1..4) {
						$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
					}
					next;
				}
				# Ratio
				$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{$quantif}{'RATIO'},$itemFormat{'number'});
				# log2 ratio
				$worksheet2->write_number($xlsRow,++$xlsCol,log($refQuantifValues->{$modProtID}{$quantif}{'RATIO'})/$log2,$itemFormat{'number'});
				# p-value
				if (defined $refQuantifValues->{$modProtID}{$quantif}{'P_VALUE'}) {
					$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{$quantif}{'P_VALUE'},$itemFormat{'number'});
				}
				else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
				# Peptides
				$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{$quantif}{$numPepCode},$itemFormat{'number'});
			}
			# MW, desc, species
			my $mass=(!$refProteinInfo->{$protID}[1])? '-' : $refProteinInfo->{$protID}[1];
			$worksheet2->write_number($xlsRow,++$xlsCol,$mass,$itemFormat{'number1d'});
			$worksheet2->write_string($xlsRow,++$xlsCol,$refProteinInfo->{$protID}[2],$itemFormat{'textWrap'});
			$worksheet2->write_string($xlsRow,++$xlsCol,$refProteinInfo->{$protID}[3],$itemFormat{'text'});
		}
		#<Identifier header (written last to get remaining number of proteins)
		#my $protTypeStrg=(scalar keys %{$refQuantifPTM})? ' forms' : '';
		#$worksheet2->merge_range(1,0,2,0,"$numProt$protTypeStrg/$numTotProteins proteins",$itemFormat{'mergeRowHeader'});
		my $protString=(scalar keys %{$refQuantifPTM})? "$numDispItems/".(scalar keys %{$refQuantifValues})." sites\n".(scalar keys %dispProteins)."/".(scalar keys %{$refProteinInfo})." proteins" : "$numDispItems/".(scalar keys %{$refProteinInfo})." proteins";
		$worksheet2->merge_range(1,0,2,0,$protString,$itemFormat{'mergeRowHeader'});
	}
	elsif ($selQuantifFamily =~ /EMPAI|SIN/) {
		my ($paramValue,$paramDesc)=($selQuantifFamily eq 'EMPAI')? ($empaiValue,$emPAIMeasures{$empaiValue}) :($selQuantifFamily eq 'SIN')? ('SIN_SIN','SIN') : ($dispMeasure,$maxQuantMeasures{$dispMeasure});
		#<Gene
		$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow+1,$xlsCol,'Gene & Synonyms',$itemFormat{'mergeRowHeader'});
		#my $invDispFoldChange=1/$dispFoldChange;
		foreach my $quantifID (@{$refSelQuantifs}) {
			my $parentItemStrg;
			for (my $i=1; $i<$#{$refQuantifInfo->{$quantifID}[3]}; $i++) { # skip project & quantif
				$parentItemStrg.=' > ' if $parentItemStrg;
				$parentItemStrg.=$refQuantifInfo->{$quantifID}[3][$i]->{'NAME'};
			}
			my $quantifName=$refQuantifInfo->{$quantifID}[0];
			$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow,$xlsCol+1,$parentItemStrg.' > '.$quantifName,$itemFormat{'mergeColHeader'});
			$worksheet2->write_string($xlsRow+1,$xlsCol,$paramDesc,$itemFormat{'header'});
			$worksheet2->write_string($xlsRow+1,++$xlsCol,'Peptides found',$itemFormat{'header'});
		}
		#<MW, description & species
		$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow+1,$xlsCol,'MW (kDa)',$itemFormat{'mergeRowHeader'});
		$worksheet2->set_column(++$xlsCol,$xlsCol,80); # col length
		$worksheet2->merge_range($xlsRow,$xlsCol,$xlsRow+1,$xlsCol,'Description',$itemFormat{'mergeRowHeader'});
		$worksheet2->set_column(++$xlsCol,$xlsCol,30); # col length
		$worksheet2->merge_range($xlsRow,$xlsCol,$xlsRow+1,$xlsCol,'Species',$itemFormat{'mergeRowHeader'});

		###<Looping through proteins>###
		$xlsRow++;
		#my $numProt=0;
		my $numDispItems=0;
		my %dispProteins;
		foreach my $protID (sort{&sortProteins($sortOrder,$refProteinInfo,$refQuantifValues)} keys %{$refQuantifValues}) { #,\%pepInfo $refProteinInfo
			next if ($selQuantifFamily eq 'MQ' && $refProteinInfo->{$protID}[0] =~ /^CON__/);
			$xlsCol=0;
			my $okFilters=0;
			my %okQuantif;
			#foreach my $quantif (@{$refSelQuantifs}) {
			#	my $usedProtID;
			#	if ($modStrg && $refQuantifPTM->{$quantif}) { # modif vs modif
			#		next if (!$refQuantifValues->{$modProtID}{$quantif} || !$refQuantifValues->{$modProtID}{$quantif}{'RATIO'});
			#		$usedProtID=$modProtID;
			#	}
			#	else { # fall back to whole prot
			#		next if (!$refQuantifValues->{$protID}{$quantif} || !$refQuantifValues->{$protID}{$quantif}{'RATIO'});
			#		$usedProtID=$protID;
			#	}
			#	#<Ratio/p-value/num pept. filters
			#	unless ($okFilters) { # no ratioPos has passed filters
			#		$okQuantif{$quantif}=(($refQuantifValues->{$usedProtID}{$quantif}{'RATIO'} < 1 && $refQuantifValues->{$usedProtID}{$quantif}{'RATIO'} <= $invDispFoldChange) || ($refQuantifValues->{$usedProtID}{$quantif}{'RATIO'} >= 1 && $refQuantifValues->{$usedProtID}{$quantif}{'RATIO'} >= $dispFoldChange))? 1 : 0;
			#		$okQuantif{$quantif}=($okQuantif{$quantif} && ($dispPvalue>=1 || (defined $refQuantifValues->{$usedProtID}{$quantif}{'P_VALUE'} && $refQuantifValues->{$usedProtID}{$quantif}{'P_VALUE'} <= $dispPvalue)))? 1 : 0;
			#		$okQuantif{$quantif}=($okQuantif{$quantif} && $refQuantifValues->{$usedProtID}{$quantif}{$numPepCode} >= $dispNumPep)? 1 : 0;
			#	}
			#	$okFilters=1 if $okQuantif{$quantif};
			#}
			#next unless $okFilters;
			#$numProt++;
			$numDispItems++;
			$dispProteins{$protID}=1;
			#<Identifier
			$worksheet2->write_string(++$xlsRow,$xlsCol,$refProteinInfo->{$protID}[0],$itemFormat{'text'});
			#<Gene
			$worksheet2->write_string($xlsRow,++$xlsCol,join(',',@{$refProteinInfo->{$protID}[5]}),$itemFormat{'text'});
			#<Quantif
			foreach my $quantif (@{$refSelQuantifs}) {
				if (!$refQuantifValues->{$protID}{$quantif} || !$refQuantifValues->{$protID}{$quantif}{$paramValue} || ($strictFilter && !$okQuantif{$quantif})) {
					foreach my $i (1..2) {
						$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
					}
					next;
				}
				# emPAI / MaxQuant (iBAQ, Intensity,etc)
				#my $numPep=($selQuantifFamily eq 'EMPAI')? $pepInfo{$quantif}{$protID} : $refQuantifValues->{$protID}{$quantif}{'NUM_PEP_TOTAL'};
				my $numPep=($selQuantifFamily =~ /EMPAI|SIN/)? $refQuantifValues->{$protID}{$quantif}{'PEPTIDES'} : $refQuantifValues->{$protID}{$quantif}{'NUM_PEP_TOTAL'};
				$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$protID}{$quantif}{$paramValue},$itemFormat{'number'});
				# Peptide number
				$worksheet2->write_number($xlsRow,++$xlsCol,$numPep,$itemFormat{'number'});
			}
			# MW, desc, species
			my $mass=(!$refProteinInfo->{$protID}[1])? '-' : $refProteinInfo->{$protID}[1];
			$worksheet2->write_number($xlsRow,++$xlsCol,$mass,$itemFormat{'number1d'});
			$worksheet2->write_string($xlsRow,++$xlsCol,$refProteinInfo->{$protID}[2],$itemFormat{'textWrap'});
			$worksheet2->write_string($xlsRow,++$xlsCol,$refProteinInfo->{$protID}[3],$itemFormat{'text'});
		}
		#<Identifier header (written last to get remaining number of proteins)
		#my $protTypeStrg=(scalar keys %{$refQuantifPTM})? ' forms' : '';
		#$worksheet2->merge_range(1,0,2,0,"$numProt$protTypeStrg/$numTotProteins proteins",$itemFormat{'mergeRowHeader'});
		my $protString=(scalar keys %{$refQuantifPTM})? "$numDispItems/".(scalar keys %{$refQuantifValues})." sites\n".(scalar keys %dispProteins)."/".(scalar keys %{$refProteinInfo})." proteins" : "$numDispItems/".(scalar keys %{$refProteinInfo})." proteins";
		$worksheet2->merge_range(0,0,1,0,$protString,$itemFormat{'mergeRowHeader'});
	}

} # end of exportProteinList

####>Check for quantif of modification<####
sub findModificationQuantifs {
	my ($refSelQuantifs,$refQuantifInfo,$refQuantifPTM)=@_;
	foreach my $quantif (@{$refSelQuantifs}) {
		my ($quantifID,$ratioPos)=split('_',$quantif);
		if ($refQuantifInfo->{$quantifID}[4]) { # $modID: 0 if whole-protein quantif;
			$refQuantifPTM->{$quantif}=$refQuantifInfo->{$quantifID}[4];
		}
	}
}

sub getModifiedProteinForms {
	my ($refQuantifValues,$refModList)=@_;
	foreach my $modProtID (keys %{$refQuantifValues}) {
		my ($protID,$modStrg)=($modProtID=~/^(\d+)(.*)/);
		$refModList->{$protID}{$modProtID}=1 if $modStrg;
	}
}

############################################################
#######################<FRAMES>#############################
############################################################
sub generateFrames {
	print header(-charset=>'utf-8');
	# warningsToBrowser(1);
	print qq
|<HTML>
<FRAMESET rows="220,*" border=0 style="overflow:scroll;">
	<FRAME name="selQuantifsFrame" src="./compareQuantifications.cgi?ACT=selQuantifs&id_project=$projectID&parentItem=$item:$itemID&quantifFamily=$selQuantifFamily&view=$view&quantifFocus=$selQuantifFocus&sort=$sortOrder" >
	<FRAME name="listFrame" src="$promsPath{html}/nothing.html">
</FRAMESET>
</HTML>
|;
	exit;
}


##########################################################################
#######################<Items selection Form>#############################
##########################################################################
sub selectQuantifications {
#print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
	####<Connecting to the database>####
	my $dbh=&promsConfig::dbConnect;

	####<Project PTMs>####
	my $sthGetPM=$dbh->prepare("SELECT PM.ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES FROM PROJECT_MODIFICATION PM,MODIFICATION M WHERE PM.ID_MODIFICATION=M.ID_MODIFICATION AND PM.ID_PROJECT=$projectID");
	$sthGetPM->execute;
	my %projectVarMods;
	while (my ($modID,$psiName,$interName,$synName)=$sthGetPM->fetchrow_array) {
		$projectVarMods{$modID}=$psiName || $interName;
		unless ($projectVarMods{$modID}) {
			$synName=~s/^##//; $synName=~s/##$//; $synName=~s/##/, /;
			$projectVarMods{$modID}=$synName;
		}
	}
	$sthGetPM->finish;

	####<List of Quantifications>####
	my $focusFilterStrg=($selModifID==0)? 'ID_MODIFICATION IS NULL' : ($selQuantifFocus > 0)? "ID_MODIFICATION=$selModifID" : "(ID_MODIFICATION=$selModifID OR ID_MODIFICATION IS NULL)";

	my @queryList;
	#<Design quantifs
	push @queryList,qq |SELECT Q.ID_QUANTIFICATION,Q.ID_QUANTIFICATION_METHOD,CONCAT('EXPERIMENT:',E.ID_EXPERIMENT),E.NAME,CONCAT('DESIGN:',D.ID_DESIGN),D.NAME
						FROM EXPERIMENT E,DESIGN D,QUANTIFICATION Q
						WHERE E.ID_PROJECT=$projectID AND E.ID_EXPERIMENT=D.ID_EXPERIMENT AND D.ID_DESIGN=Q.ID_DESIGN AND Q.FOCUS='protein' AND Q.STATUS=1
							AND $focusFilterStrg
						ORDER BY Q.ID_QUANTIFICATION_METHOD ASC,E.DISPLAY_POS ASC,D.NAME ASC|;
	#<Internal quantifs with Samples
	push @queryList,qq |SELECT DISTINCT Q.ID_QUANTIFICATION,Q.ID_QUANTIFICATION_METHOD,CONCAT('EXPERIMENT:',E.ID_EXPERIMENT),E.NAME,CONCAT('SAMPLE:',S.ID_SAMPLE),S.NAME,CONCAT('ANALYSIS:',A.ID_ANALYSIS),A.NAME
						FROM EXPERIMENT E,SAMPLE S,ANALYSIS A,ANA_QUANTIFICATION AQ,QUANTIFICATION Q
						WHERE E.ID_PROJECT=$projectID AND E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SAMPLE=A.ID_SAMPLE AND S.ID_SPOT IS NULL
							AND A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.FOCUS='protein' AND Q.STATUS=1 AND ID_DESIGN IS NULL
							AND $focusFilterStrg
						ORDER BY Q.ID_QUANTIFICATION_METHOD ASC,E.DISPLAY_POS ASC,S.DISPLAY_POS ASC,A.DISPLAY_POS ASC|;
	#<Internal quantifs with Gels
	#,CONCAT('SPOT:',SP.ID_SPOT),SP.NAME,CONCAT('ANALYSIS:',A.ID_ANALYSIS),A.NAME
	push @queryList,qq |SELECT DISTINCT Q.ID_QUANTIFICATION,Q.ID_QUANTIFICATION_METHOD,CONCAT('EXPERIMENT:',E.ID_EXPERIMENT),E.NAME,CONCAT('GEL2D:',G.ID_GEL2D),G.NAME
						FROM EXPERIMENT E,GEL2D G,SPOT SP,SAMPLE S,ANALYSIS A,ANA_QUANTIFICATION AQ,QUANTIFICATION Q
						WHERE E.ID_PROJECT=$projectID AND E.ID_EXPERIMENT=G.ID_EXPERIMENT AND G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT AND S.ID_SAMPLE=A.ID_SAMPLE
							AND A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.FOCUS='protein' AND Q.STATUS=1 AND ID_DESIGN IS NULL
							AND $focusFilterStrg
						ORDER BY Q.ID_QUANTIFICATION_METHOD ASC,E.DISPLAY_POS ASC,G.DISPLAY_POS ASC,SP.NAME ASC,A.DISPLAY_POS ASC|;
	my $sthQM=$dbh->prepare('SELECT CODE FROM QUANTIFICATION_METHOD WHERE ID_QUANTIFICATION_METHOD=?');
	my (%seenQuantifs,%usedQuantifFamilies,%existQuantifFamilies,@itemList,%usedItems);
	foreach my $query (@queryList) {
		my $sthItem=$dbh->prepare($query);
		$sthItem->execute;
		while (my ($quantifID,$qmethID,@itemInfo)=$sthItem->fetchrow_array) {
			next if $seenQuantifs{$quantifID};
			$seenQuantifs{$quantifID}=1;
			$sthQM->execute($qmethID);
			my ($quantifCode)=$sthQM->fetchrow_array;
next if $quantifCode eq 'SSPA'; # skip SSPA... for now
			my $qFamily;
			FAM:foreach my $qFam (keys %{$proteinQuantifFamilies{'MEMBERS'}}) {
				foreach my $qmCode (@{$proteinQuantifFamilies{'MEMBERS'}{$qFam}}) {
					if ($qmCode eq $quantifCode) {
						$qFamily=$qFam;
						last FAM;
					}
				}
			}
			$existQuantifFamilies{$qFamily}=1;
			$selQuantifFamily=$qFamily unless $selQuantifFamily;
			next if $qFamily ne $selQuantifFamily;
			#$quantifMethods{$qmethID}='';
			#$selQuantifMethodID=$qmethID unless $selQuantifMethodID; # if not selected: take 1st one
			#next if $qmethID != $selQuantifMethodID;
			for (my $i=0; $i<$#itemInfo; $i+=2) {
				next if $usedItems{$itemInfo[$i]};
				my $labelString='';
				for (my $j=0; $j<=$i-2; $j+=2) {
					$labelString.=' > ' if $labelString;
					$labelString.=&promsMod::shortenName($itemInfo[$j+1],15);
				}
				$labelString.=' > ' if $labelString;
				$labelString.=$itemInfo[$i+1];
				push @itemList,[$itemInfo[$i],$labelString];
				$usedItems{$itemInfo[$i]}=1;
			}
		}
		$sthItem->finish;
	}
	$sthQM->finish;
	my $usedParent='';
	if ($item eq 'project') {
		$usedParent=$itemList[0][0] if $itemList[0];
	}
	else {$usedParent=uc($parentItem);}

	$dbh->disconnect;

	##############
	####<HTML>####
	##############
	print header(-charset=>'utf-8');
	warningsToBrowser(1); # writes PERL warnings as HTML comments <!-- --> must be activated after header
	my $ITEM=uc($item);
	print qq
|<HEAD>
<TITLE>Compare Multiple Quantifications</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE>
.trueFilter{color:#DD0000;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
|;
	&promsMod::popupInfo();
	print qq
|function changeQuantifFamily(newFamily) {
	//top.promsFrame.selectedComparison='items'; //openProject.cgi
	parent.listFrame.location="$promsPath{html}/nothing.html";
	window.location="$promsPath{cgi}/compareQuantifications.cgi?ACT=selQuantifs&id_project=$projectID&parentItem=$item:$itemID&quantifFamily="+newFamily+"&view=$view&quantifFocus=$selQuantifFocus";
}
function changeQuantifFocus(newFocus) {
	parent.listFrame.location="$promsPath{html}/nothing.html";
	window.location="$promsPath{cgi}/compareQuantifications.cgi?ACT=selQuantifs&id_project=$projectID&parentItem=$item:$itemID&quantifFamily=$selQuantifFamily&view=$view&quantifFocus="+newFocus;
}
function updateQuantifSelection(action) {
	var allSelect=document.compForm.allQuantifs;
	var usedSelect=document.compForm.usedQuantifs;
	var selected=false;
	if (action=='add') {
		for (var i=0; i<usedSelect.length; i++) {usedSelect.options[i].selected=false;}
		//Processing & adding parent label
		var selQuantifParent=document.compForm.quantifParent;
		var parentLabel=selQuantifParent.options[selQuantifParent.options.selectedIndex].text;
		var labelList=parentLabel.split(' > ');
		var lastLabel=labelList[labelList.length-1];
		if (lastLabel.length>11) {
			var begStrg=lastLabel.substr(0,4);
			var endStrg=lastLabel.substr(lastLabel.length-4,4);
			labelList[labelList.length-1]=begStrg+'...'+endStrg;
		}
		parentLabel=labelList.join(' > ');
		ALL:for (var i=0; i<allSelect.length; i++) {
			if (allSelect.options[i].selected) {
				selected=true;
				//Check if not already used
				for (var j=0; j<usedSelect.length; j++) {
					if (usedSelect.options[j].value==allSelect.options[i].value) continue ALL;
				}
				usedSelect.options[usedSelect.length]=new Option(parentLabel+' > '+allSelect.options[i].text,allSelect.options[i].value);
				usedSelect.options[usedSelect.length-1].selected=true;
				allSelect.options[i].selected=false;
			}
		}
		usedSelect.focus();
	}
	else { //remove
		var keepOptions=new Array();
		for (var i=0; i<usedSelect.length; i++) {
			if (!usedSelect.options[i].selected) {
				keepOptions.push(usedSelect.options[i].text);
				keepOptions.push(usedSelect.options[i].value);
			}
			else {selected=true;}
		}
		usedSelect.length=0;
		for (var i=0; i<keepOptions.length; i+=2) {
			usedSelect.options[usedSelect.length]=new Option(keepOptions[i],keepOptions[i+1]);
		}
	}
	if (!selected) {alert('No Analysis selected');}
}
function moveQuantifications(direction) {
	//setComparisonAsModified(2);
	var usedSelect=document.compForm.usedQuantifs;
	var selectOK=false;
	if (direction=='up') {
		for (var i=0; i<usedSelect.length; i++) {
			selectOK=switchQuantifications(usedSelect,i,-1,selectOK);
		}
	}
	else {
		for (var i=usedSelect.length-1; i>=0; i--) {
			selectOK=switchQuantifications(usedSelect,i,1,selectOK);
		}
	}
	usedSelect.focus();
	if (!selectOK) {alert('No Quantification selected');}
}
function clearQuantifications() {
	//setComparisonAsModified(2);
	document.compForm.usedQuantifs.length=0;
}
function switchQuantifications(usedSelect,i,delta,selectOK) {
	if (usedSelect.options[i].selected) {
		selectOK=true;
		if (i+delta<0 \|\| i+delta==usedSelect.length \|\| usedSelect.options[i+delta].selected) return selectOK;
		var prevQuantifID=usedSelect.options[i+delta].value;
		var prevText=usedSelect.options[i+delta].text;
		usedSelect.options[i+delta]=new Option(usedSelect.options[i].text,usedSelect.options[i].value);
		usedSelect.options[i]=new Option(prevText,prevQuantifID);
		usedSelect.options[i+delta].selected=true;
	}
	return selectOK;
}
var viewFilters={
	list:['ratio','pValue'],
	//heatmap:['ratio','pValue'],
	volcano:[],
	log2:['pValue'],
	correlMatrix:['pValue']
};
function updateView(view) {
	// Updating filters color
	var filters=['ratio','pValue'];
	var tags=['TEXT','INPUT'];
	for (var i=0; i<filters.length; i++) {
		for (var j=0; j<tags.length; j++) {
			if (document.getElementById(filters[i]+tags[j])) {document.getElementById(filters[i]+tags[j]).className='';}
		}
	}
	for (var i=0; i < viewFilters[view].length; i++) {
		for (var j=0; j<tags.length; j++) {
			if (document.getElementById(viewFilters[view][i]+tags[j])) {document.getElementById(viewFilters[view][i]+tags[j]).className='trueFilter';}
		}
	}
	// Updating strict filtering
	var myForm=document.compForm;
	if (myForm.strictFilter) {myForm.strictFilter.disabled=(view=='volcano')? true : false;}
	// Checking & submitting form
	var okSubmit=checkSelectedQuantifications(myForm);
	if (okSubmit) myForm.submit();
}
function checkSelectedQuantifications(myForm) {
	if (myForm.usedQuantifs.length<=1) {
		alert('Select at least 2 Quantifications.');
		return false;
	}

	myForm.quantifList.value='';
	for (var i=0; i<myForm.usedQuantifs.length; i++) {
		myForm.quantifList.value+=myForm.usedQuantifs[i].value;
		if (i<myForm.usedQuantifs.length-1) myForm.quantifList.value+=':';
	}
/*
	if (modCompType==1) { // modified by listFrame
		myForm.comparisonID.options[selCompIndex].text=selCompName; // removes '*' from name
	}
*/
	if (myForm.numPep && myForm.numPep.value <= 0) {myForm.numPep.value=1;}
	//if (myForm.numPep.value == 1 && myForm.view.value=='volcano') {myForm.numPep.value=2;} // no p-value if numPep < 2
	return true;
}
function autoSelectQuantifications() {
	var myForm=document.compForm;
	var searchString=myForm.selectText.value;
	if (!searchString) {
		alert('No search string typed.');
		return;
	}
	var bs=String.fromCharCode(92);
	var unsafe=bs+".+*?[^]\$(){}=!<>\|";
	for (var i=0;i<unsafe.length;++i) {
		searchString=searchString.replace(new RegExp("\\\\"+unsafe.charAt(i),"g"),bs+unsafe.charAt(i));
	}
	//Searching option text
	var selStatus=(myForm.autoAction.value=='select')? true : false;
	for (var i=0; i<myForm.allQuantifs.length; i++) {
		if (myForm.allQuantifs.options[i].text.match(searchString)) {
			myForm.allQuantifs.options[i].selected=selStatus;
		}
	}
}
function clearSelection() {
	var myForm=document.compForm;
	for (var i=0; i<myForm.allQuantifs.length; i++) {
		myForm.allQuantifs.options[i].selected=false;
	}
}
// AJAX --->
function ajaxUpdateRestrict() {
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxRestrictList&projectID=$projectID&restrictList=$restrictListID",true); // !!! Uses code from showProtQuantification.cgi !!!
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			document.getElementById('restrictDIV').innerHTML=XHR.responseText;
		}
	}
	XHR.send(null);
}
function ajaxGetQuantificationList(branchID) {
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/compareQuantifications.cgi?AJAX=getquantifList&quantifFamily=$selQuantifFamily&branchID="+branchID,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			listAllQuantifications(XHR.responseText);
		}
	}
	XHR.send(null);
}
/*
function ajaxUpdateCompGroupAnalyses(newGrNumber,compID) {
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/compareQuantifications.cgi?AJAX=getAnaText&compID="+compID,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			storeGroupAnalyses(newGrNumber,XHR.responseText);
		}
	}
	XHR.send(null);
}
*/
function getXMLHTTP() {
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
function listAllQuantifications(quantifText) {
	var quantifList=quantifText.split('\\n');
	var selAllQuantifs=document.compForm.allQuantifs;
	selAllQuantifs.length=0;
	for (var i=0; i<quantifList.length-1; i++) { // -1 because last line is empty
		var anaInfo=quantifList[i].split('#:#');
		selAllQuantifs.options[i]=new Option(anaInfo[1],anaInfo[0]);
		selAllQuantifs.options[i].selected=true;
	}
}
// <--- AJAX
</SCRIPT>
</HEAD>
<BODY topmargin=2 background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FORM name="compForm" method="post" target="listFrame" onsubmit="return checkSelectedQuantifications(this)">
<DIV style="white-space:nowrap">
<FONT class="title">Quantification method:</FONT><SELECT name="quantifFamily" onchange="changeQuantifFamily(this.value)" style="font-weight:bold;font-size:20px;color:#DD0000">
|;
	unless ($selQuantifFamily) {print '<OPTION value="">None found</OPTION>';}
	foreach my $qFamily (sort{lc($proteinQuantifFamilies{'NAME'}{$a}) cmp lc($proteinQuantifFamilies{'NAME'}{$b})} keys %{$proteinQuantifFamilies{'NAME'}}) {
		print "<OPTION value=\"$qFamily\"";
		print " selected" if $qFamily eq $selQuantifFamily;
		print " disabled" unless $existQuantifFamilies{$qFamily};
		print ">$proteinQuantifFamilies{NAME}{$qFamily}</OPTION>\n";
	}
	print "</SELECT>\n";
#print qq
#|&nbsp;&nbsp;<FONT class="title2">Comparison:<SELECT name="comparisonID" style="font-weight:bold;font-size:16px;width:250px" onchange="selectComparison(this.value,true)"><OPTION value="0">New</OPTION>
#|;
#	foreach my $compID (sort{lc($comparisonList{$a}{'NAME'}) cmp lc($comparisonList{$b}{'NAME'})} keys %comparisonList) {
#		print "<OPTION value=\"$compID\"";
#		print ' selected' if $compID==$comparisonID;
#		print ">$comparisonList{$compID}{NAME}</OPTION>\n";
#	}
#	print qq
#|</SELECT>&nbsp;<INPUT type="button" class="title3" value="View All" onclick="displayAllComparisons()"/>
#</FONT>
#|;
	if (!$selQuantifFamily || $selQuantifFamily eq 'RATIO') { # empty string if no quantif in selected family
		print qq
|&nbsp;&nbsp;&nbsp;&nbsp;<FONT class="title">Focus:<SELECT name="quantifFocus" onchange="changeQuantifFocus(this.value)" style="font-weight:bold;font-size:20px;color:#DD0000"><OPTION value="0">Proteins</OPTION>
|;
		foreach my $modID (sort{lc($projectVarMods{$a}) cmp lc($projectVarMods{$b})} keys %projectVarMods) {
			my ($sel1Strg,$sel2Strg)=($selQuantifFocus==$modID)? (' selected','') : ($selQuantifFocus==-$modID)? ('',' selected') : ('','');
			print "<OPTION value=\"$modID\"$sel1Strg>$projectVarMods{$modID}-proteins</OPTION>\n";
			print "<OPTION value=\"-$modID\"$sel2Strg>$projectVarMods{$modID}-proteins & proteins</OPTION>\n";
		}
		print "</SELECT>\n";
	}
	print qq
|</DIV>
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="parentItem" value="$parentItem">
<INPUT type="hidden" name="ACT" value="results">
<INPUT type="hidden" name="quantifList" value="">
<INPUT type="hidden" name="sort" value="$sortOrder">
<INPUT type="hidden" name="chkProtList" value="">
<FONT style="font-size:7px"><BR></FONT>
<TABLE bgcolor=$darkColor cellpadding=0 border=0><TR>
<TR>
<TH valign=top rowspan=2><FONT class="title2">Select :</FONT><SELECT name="quantifParent" class="title3" style="width:350px" onchange="ajaxGetQuantificationList(this.value)">
|;
	####<Looping through list of analyses>####
	my $matchedItem=0;
	foreach my $refItem (@itemList) {
		print "\t<OPTION value=\"$refItem->[0]\"";
		if ($refItem->[0] eq $usedParent) {
			print ' selected';
			$matchedItem=1;
		}
		print ">$refItem->[1]</OPTION>\n";
	}
	$usedParent=$itemList[0]->[0] unless $matchedItem; # Changed parent if no matching quantif
	my $volcanoOptionStrg='';
	if ($selQuantifFamily eq 'RATIO') {
		$volcanoOptionStrg='<OPTION value="volcano"';
		$volcanoOptionStrg.=' selected' if $view eq 'volcano';
		$volcanoOptionStrg.='>Volcano plot</OPTION>';
	}
	elsif ($view eq 'volcano') {$view='list';}
	my $correlMatOptionStrg='';
	if ($selQuantifFocus >= 0) {
		$correlMatOptionStrg='<OPTION value="correlMatrix"';
		$correlMatOptionStrg.=' selected' if $view eq 'correlMatrix';
		$correlMatOptionStrg.='>Correlation matrix</OPTION>';
	}
	#my $selHeat=($view eq 'heatmap')? ' selected' : '';
	my $selLog2=($view eq 'log2')? ' selected' : '';
	my $disabStrictFilter=($view eq 'volcano')? ' disabled' : '';
	print qq
|</SELECT><BR>
<SELECT multiple name="allQuantifs" style="width:600px;height:105px;font-weight:bold;font-size:14px;"></SELECT><BR>
<SELECT name="autoAction" style="font-weight:bold;"><OPTION value="select">Select</OPTION><OPTION value="unselect">Unselect</OPTION></SELECT> items containing
<INPUT type="text" name="selectText" size=10 value=""/>&nbsp;<INPUT type="button" value="Go" onclick="autoSelectQuantifications()"/>&nbsp;&nbsp<INPUT type="button" value="Clear" onclick="clearSelection()"/></TH>

<TH valign=middle><BR><INPUT type="button" name="add" value=">" style="width:30px" onclick="updateQuantifSelection('add')"/><BR><INPUT type="button" name="remove" value="<" style="width:30px" onclick="updateQuantifSelection('remove')"/></TH>

<TH valign=top><FONT class="title2">Datasets to compare:</FONT><BR>
<SELECT multiple name="usedQuantifs" style="width:600px;height:80px;font-weight:bold;font-size:14px;"></SELECT></TH>
<TD valign=middle><BR><INPUT type="button" name="up" value="Up" style="width:50px" onclick="moveQuantifications('up')"/><BR><INPUT type="button" name="down" value="Down" style="width:50px" onclick="moveQuantifications('down')"/><BR><INPUT type="button" value="Clear" style="width:50px" onclick="clearQuantifications()"/></TD>
</TR>
<TR><TH colspan=3><TABLE cellpadding=0 cellspacing=0 width=100% border=0>
|;
	if ($selQuantifFamily eq 'RATIO') {
		print qq
|<TR><TH align=left nowrap colspan=4>
&nbsp;<FONT id="ratioTEXT" class="trueFilter" onmouseover="popup('<B>Up &ge;</B> or <B>Down &le; 1/</B>')" onmouseout="popout()">&bull;Abs.<SUP>*</SUP> fold change &ge;</FONT><INPUT type="text" id="ratioINPUT" name="foldChange" class="trueFilter" value="2" size=2/>&nbsp;&nbsp;&nbsp;<FONT id="pValueTEXT" class="trueFilter">&bull;p-value &le;</FONT><INPUT type="text" id="pValueINPUT" name="pValue" class="trueFilter" value="0.05" size=5/>
&nbsp;&nbsp;&nbsp;<FONT class="trueFilter">&bull;</FONT><SELECT name="pepType" class="trueFilter"><OPTION value="all">All</OPTION><OPTION value="distinct">Dist.</OPTION></SELECT> <FONT class="trueFilter">pept. used &ge;</FONT><INPUT type="text" name="numPep" class="trueFilter" value="3" size=2/>
&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="strictFilter" value="1"$disabStrictFilter/>Strict filtering
</TH></TR>
|;
	}
	elsif ($selQuantifFamily eq 'MQ') {
		print qq
|<TR><TH align=left nowrap colspan=4>
&nbsp;<FONT class="title3">&bull;Measure :</FONT><SELECT name="dispMeasure" class="title3">
|;
		foreach my $refMeas (@{$proteinQuantifFamilies{'MEASURES'}{'MQ'}}) {print "<OPTION value=\"$refMeas->[0]\">$refMeas->[1]</OPTION>";}
		print qq
|</SELECT>
&nbsp;&nbsp;&nbsp;<FONT class="trueFilter">&bull;</FONT><SELECT name="pepType" class="trueFilter"><OPTION value="RAZ_UNI_PEP">Razor+unique</OPTION><OPTION value="UNIQUE_PEP">Unique</OPTION><OPTION value="PEPTIDES">All</OPTION></SELECT> <FONT class="trueFilter">peptides &ge;</FONT><INPUT type="text" name="numPep" class="trueFilter" value="3" size=2/>
&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="strictFilter" value="1"$disabStrictFilter/>Strict filtering
</TH></TR>
|;
	}
	elsif ($selQuantifFamily eq 'EMPAI') {
		print qq
|<TR><TH align=left nowrap colspan=4>
&nbsp;<FONT id="ratioTEXT" class="trueFilter">&bull;emPAI type :</FONT><SELECT name="empaiType" class="trueFilter">
|;
#<OPTION value="EMPAI_MR">$emPAIMeasures{EMPAI_MR}</OPTION><OPTION value="EMPAI_MOL">$emPAIMeasures{EMPAI_MOL}</OPTION><OPTION value="EMPAI">$emPAIMeasures{EMPAI}</OPTION>
		foreach my $refMeas (@{$proteinQuantifFamilies{'MEASURES'}{'EMPAI'}}) {print "<OPTION value=\"$refMeas->[0]\">$refMeas->[1]</OPTION>";}
		print qq
|</SELECT>
</TH></TR>
|;
	}
	print qq
|<TR><TH align=left nowrap width=10%>&nbsp;&bull;<SELECT name="listAction" class="title3"><OPTION value="restrict">Restrict to</OPTION><OPTION value="exclude">Exclude</OPTION></SELECT>:</TH><TD><DIV id="restrictDIV"><!--Custom list selection comes here with window.onload() --></DIV></TD>
<TH align=right><FONT class="title3">View:</FONT><SELECT name="view" class="title3" onchange="updateView(this.value)">
	<OPTION value="list">List</OPTION>
	<!--OPTION value="heatmap"\$selHeat>Heat map</OPTION-->
	$volcanoOptionStrg
	<OPTION value="log2"$selLog2>Correlation plot</OPTION>
	$correlMatOptionStrg
</SELECT>&nbsp;</TH>
<TH><INPUT type="submit" name="compare" value="Compare" class="title2" style="height:30px;width:120px"></TH></TR></TABLE></TH></TR>
</TABLE>
</FORM>
|;
	if ($selQuantifFamily) {
		print qq
|<SCRIPT LANGUAGE="JavaScript">
window.onload=function() {
	ajaxUpdateRestrict();
	ajaxGetQuantificationList('$usedParent');
}
</SCRIPT>
|;
	}
	print qq
|<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup();</SCRIPT>
</BODY>
</HTML>
|;
	exit;
} # end of &selectQuantifications


sub sortProteins {
	my ($sort,$refProteinInfo,$refQuantifValues)=@_; #,$refPepInfo
	my ($pID_a,$mStrg_a)=($a=~/^(\d+)(.*)/); $mStrg_a='' unless $mStrg_a;
	my ($pID_b,$mStrg_b)=($b=~/^(\d+)(.*)/); $mStrg_b='' unless $mStrg_b;
	if ($sort =~ /ratio_(\d+_\d+)/) {
		my $quantif=$1;
		if ($refQuantifValues->{$a}{$quantif} && $refQuantifValues->{$b}{$quantif}) {
			$refQuantifValues->{$b}{$quantif}{'RATIO'}<=>$refQuantifValues->{$a}{$quantif}{'RATIO'} || $refQuantifValues->{$b}{$quantif}{$numPepCode}<=>$refQuantifValues->{$a}{$quantif}{$numPepCode} || lc($refProteinInfo->{$pID_a}[0]) cmp lc($refProteinInfo->{$pID_b}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($mStrg_a),&promsMod::preparePtmString($mStrg_b))
		}
		elsif ($refQuantifValues->{$a}{$quantif}) {-1}
		elsif ($refQuantifValues->{$b}{$quantif}) {1}
		else {lc($refProteinInfo->{$pID_a}[0]) cmp lc($refProteinInfo->{$pID_b}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($mStrg_a),&promsMod::preparePtmString($mStrg_b))}
	}
	elsif ($sort =~ /state_(\d+_\d+)/) {
		my $quantif=$1;
		if ($refQuantifValues->{$a}{$quantif} && $refQuantifValues->{$b}{$quantif}) {
			$refQuantifValues->{$b}{$quantif}{$dispMeasure}<=>$refQuantifValues->{$a}{$quantif}{$dispMeasure} || $refQuantifValues->{$b}{$quantif}{$numPepCode}<=>$refQuantifValues->{$a}{$quantif}{$numPepCode} || lc($refProteinInfo->{$pID_a}[0]) cmp lc($refProteinInfo->{$pID_b}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($mStrg_a),&promsMod::preparePtmString($mStrg_b))
		}
		elsif ($refQuantifValues->{$a}{$quantif}) {-1}
		elsif ($refQuantifValues->{$b}{$quantif}) {1}
		else {lc($refProteinInfo->{$pID_a}[0]) cmp lc($refProteinInfo->{$pID_b}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($mStrg_a),&promsMod::preparePtmString($mStrg_b))}
	}
	# defaults to peptide NUM_PEP_USED or DIST_PEP_USED
	elsif ($sort =~ /peptide_(\d+_\d+)/) {
		my $quantif=$1;
		if ($refQuantifValues->{$a}{$quantif} && $refQuantifValues->{$b}{$quantif}) {
			$refQuantifValues->{$b}{$quantif}{$numPepCode}<=>$refQuantifValues->{$a}{$quantif}{$numPepCode} || $refQuantifValues->{$b}{$quantif}{'RATIO'}<=>$refQuantifValues->{$a}{$quantif}{'RATIO'} || lc($refProteinInfo->{$pID_a}[0]) cmp lc($refProteinInfo->{$pID_b}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($mStrg_a),&promsMod::preparePtmString($mStrg_b))
		}
		elsif ($refQuantifValues->{$a}{$quantif}) {-1}
		elsif ($refQuantifValues->{$b}{$quantif}) {1}
		else {lc($refProteinInfo->{$pID_a}[0]) cmp lc($refProteinInfo->{$pID_b}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($mStrg_a),&promsMod::preparePtmString($mStrg_b))}
	}
	elsif ($sort eq 'mass') {lc($refProteinInfo->{$b}[1]) <=> lc($refProteinInfo->{$a}[1])}
	elsif ($sort =~ /empai_(\d+)/ || $sort =~ /sin_(\d+)/){
		my $quantif=$1;
		my $qVal=($sort =~ /empai/)?$empaiValue:'SIN_SIN';
		if ($refQuantifValues->{$a}{$quantif} && $refQuantifValues->{$b}{$quantif}) {
			$refQuantifValues->{$b}{$quantif}{$qVal}<=>$refQuantifValues->{$a}{$quantif}{$qVal} || $refQuantifValues->{$b}{$quantif}{'PEPTIDES'}<=>$refQuantifValues->{$a}{$quantif}{'PEPTIDES'} || lc($refProteinInfo->{$pID_a}[0]) cmp lc($refProteinInfo->{$pID_b}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($mStrg_a),&promsMod::preparePtmString($mStrg_b))
		}
		elsif ($refQuantifValues->{$a}{$quantif}) {-1}
		elsif ($refQuantifValues->{$b}{$quantif}) {1}
		else {lc($refProteinInfo->{$pID_a}[0]) cmp lc($refProteinInfo->{$pID_b}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($mStrg_a),&promsMod::preparePtmString($mStrg_b))}
	}
	else { # default to protein identifier
		lc($refProteinInfo->{$pID_a}[0]) cmp lc($refProteinInfo->{$pID_b}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($mStrg_a),&promsMod::preparePtmString($mStrg_b))
	}
}


##########################################################################
#######################<AJAX-called routines>#############################
##########################################################################
sub ajaxGetQuantificationList {
	my ($parentItem,$parentID)=split(':',param('branchID'));

	my @queryList;
	if ($parentItem eq 'EXPERIMENT') {
		push @queryList,qq |SELECT Q.ID_QUANTIFICATION,Q.QUANTIF_ANNOT,D.NAME,Q.NAME
							FROM DESIGN D,QUANTIFICATION Q
							WHERE Q.ID_QUANTIFICATION_METHOD=? AND D.ID_EXPERIMENT=$parentID AND D.ID_DESIGN=Q.ID_DESIGN AND Q.FOCUS='protein' AND Q.STATUS=1
							ORDER BY D.NAME ASC,Q.NAME ASC|;
		push @queryList,qq |SELECT Q.ID_QUANTIFICATION,Q.QUANTIF_ANNOT,S.NAME,A.NAME,Q.NAME
							FROM SAMPLE S,ANALYSIS A,ANA_QUANTIFICATION AQ,QUANTIFICATION Q
							WHERE Q.ID_QUANTIFICATION_METHOD=? AND S.ID_EXPERIMENT=$parentID AND S.ID_SAMPLE=A.ID_SAMPLE AND S.ID_SPOT IS NULL
								AND A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.FOCUS='protein' AND Q.STATUS=1
							ORDER BY S.DISPLAY_POS ASC,A.DISPLAY_POS ASC,Q.NAME ASC|;
		push @queryList,qq |SELECT Q.ID_QUANTIFICATION,Q.QUANTIF_ANNOT,A.ID_ANALYSIS,G.NAME,SP.NAME,A.NAME,Q.NAME
							FROM GEL2D G,SPOT SP,SAMPLE S,ANALYSIS A,ANA_QUANTIFICATION AQ,QUANTIFICATION Q
							WHERE Q.ID_QUANTIFICATION_METHOD=? AND G.ID_EXPERIMENT=$parentID AND G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT AND S.ID_SAMPLE=A.ID_SAMPLE
								AND A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.FOCUS='protein' AND Q.STATUS=1
							ORDER BY G.DISPLAY_POS ASC,SP.NAME ASC,A.DISPLAY_POS ASC,Q.NAME ASC|;
	}
	elsif ($parentItem eq 'DESIGN') {
		push @queryList,qq |SELECT ID_QUANTIFICATION,QUANTIF_ANNOT,NAME
								FROM QUANTIFICATION
								WHERE ID_QUANTIFICATION_METHOD=? AND ID_DESIGN=$parentID AND FOCUS='protein' AND STATUS=1
								ORDER BY NAME ASC|;
	}
	elsif ($parentItem eq 'GEL2D') {
		push @queryList,qq |SELECT Q.ID_QUANTIFICATION,Q.QUANTIF_ANNOT,SP.NAME,A.NAME,Q.NAME
								FROM SPOT SP,SAMPLE S,ANALYSIS A,ANA_QUANTIFICATION AQ,QUANTIFICATION Q
								WHERE Q.ID_QUANTIFICATION_METHOD=? AND ID_DESIGN IS NULL AND SP.ID_GEL2D=$parentID AND SP.ID_SPOT=S.ID_SPOT AND S.ID_SAMPLE=A.ID_SAMPLE
									AND A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.FOCUS='protein' AND Q.STATUS=1
								ORDER BY SP.NAME ASC,A.DISPLAY_POS ASC,Q.NAME ASC|;
	}
	elsif ($parentItem eq 'SPOT') {
		push @queryList,qq |SELECT Q.ID_QUANTIFICATION,Q.QUANTIF_ANNOT,A.NAME,Q.NAME
							FROM SAMPLE S,ANALYSIS A,ANA_QUANTIFICATION AQ,QUANTIFICATION Q
							WHERE Q.ID_QUANTIFICATION_METHOD=? AND ID_DESIGN IS NULL AND ID_DESIGN IS NULL AND S.ID_SPOT=$parentID AND S.ID_SAMPLE=A.ID_SAMPLE
								AND A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.FOCUS='protein' AND Q.STATUS=1
							ORDER BY A.DISPLAY_POS ASC,Q.NAME ASC|;
	}
	elsif ($parentItem eq 'SAMPLE') {
		push @queryList,qq |SELECT Q.ID_QUANTIFICATION,Q.QUANTIF_ANNOT,A.NAME,Q.NAME
							FROM ANALYSIS A,ANA_QUANTIFICATION AQ,QUANTIFICATION Q
							WHERE Q.ID_QUANTIFICATION_METHOD=? AND ID_DESIGN IS NULL AND A.ID_SAMPLE=$parentID AND A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.FOCUS='protein' AND Q.STATUS=1
							ORDER BY A.DISPLAY_POS ASC,Q.NAME ASC|;
	}
	elsif ($parentItem eq 'ANALYSIS') {
		push @queryList,"SELECT Q.ID_QUANTIFICATION,Q.QUANTIF_ANNOT,Q.NAME FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE Q.ID_QUANTIFICATION_METHOD=? AND ID_DESIGN IS NULL AND AQ.ID_ANALYSIS=$parentID AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.STATUS=1 ORDER BY Q.NAME ASC";
	}


	####<Starting HTML
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);

	my %quantifCodeIDs;
	#my $sthQM=$dbh->prepare("SELECT CODE,ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD");
my $sthQM=$dbh->prepare("SELECT CODE,ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE != 'SSPA'"); # skip SSPA... for now
	$sthQM->execute;
	while (my ($qmCode,$qmID)=$sthQM->fetchrow_array) {$quantifCodeIDs{$qmCode}=$qmID;}
	$sthQM->finish;

	#my ($quantifMethod)=$dbh->selectrow_array("SELECT CODE FROM QUANTIFICATION_METHOD WHERE ID_QUANTIFICATION_METHOD=$selQuantifMethodID");

	my %seenQuantifs;
	foreach my $query (@queryList) {
		my $sthQ=$dbh->prepare($query);
		foreach my $qmCode (@{$proteinQuantifFamilies{'MEMBERS'}{$selQuantifFamily}}) {
			$sthQ->execute($quantifCodeIDs{$qmCode});
			while (my ($quantifID,$quantifAnnot,@quantifPath)=$sthQ->fetchrow_array) {
				next if $seenQuantifs{$quantifID};
				$seenQuantifs{$quantifID}=1;
				for (my $i=0; $i<$#quantifPath-1; $i++) {
					$quantifPath[$i]=&promsMod::shortenName($quantifPath[$i],15);
				}
				my $pathStrg=join (' > ',@quantifPath);

				##<Extracting ratios used
				if ($selQuantifFamily =~ /ratio/i) {
					my (%labelingInfo,%stateInfo);
					&promsQuantif::extractQuantificationParameters($dbh,$quantifAnnot,\%labelingInfo,\%stateInfo);
					my $ratioPos=0;
					foreach my $ratioData (@{$labelingInfo{'RATIOS'}}) {
						$ratioPos++;
						my ($testStatePos,$refStatePos)=split(/\//,$ratioData);
						my $normTag='';
						if ($ratioData=~/%/) { # Super ratio
							$normTag=($exportType)? '°' : encode_utf8('°');
							$testStatePos=~s/%\d+//;
							$refStatePos=~s/%\d+//;
						}
						print $quantifID,"_$ratioPos#:#$pathStrg : $stateInfo{$testStatePos}{NAME}$normTag/$stateInfo{$refStatePos}{NAME}$normTag\n";
					}
				}
				elsif ($selQuantifFamily eq 'MQ') { # MaxQuant intensities
					my (%labelingInfo,%stateInfo);
					&promsQuantif::extractQuantificationParameters($dbh,$quantifAnnot,\%labelingInfo,\%stateInfo);
					my $statePos=0;
					foreach my $stateData (@{$labelingInfo{'STATES'}}) {
						$statePos++;
						print $quantifID,"_$statePos#:#$pathStrg : $stateInfo{$statePos}{NAME}\n";
					}
				}
				else { # SIN, TNPQ emPAI
					print "$quantifID#:#$pathStrg\n";
				}
			}
		}
		$sthQ->finish;
	}

	$dbh->disconnect;
	exit;
}

sub ajaxListSelectedProteins { # also called by displayClustering.cgi & displayPCA.cgi
	my $call=param('CALL') || 'COMPARE'; # can be PCA if called from PCA interface
	my (%selectedProteins,%quantifInfo,%quantifValues,%proteinInfo,%quantifAllNames);
	foreach my $protID (split(',',param('selProt'))) {
		$selectedProteins{$protID}=1;
	}
	my %parameters=(QUANTIF_FAMILY=>$selQuantifFamily,VIEW=>$view,NUM_PEP_CODE=>$numPepCode,QUANTIF_LIST=>\@selectedQuantifications,SEL_MODIF_ID=>$selModifID); # STRICT_SELECT to match modProtID not loose protID
	$parameters{MEASURE}=$dispMeasure if $selQuantifFamily eq 'MQ';
	&promsQuantif::fetchQuantificationData($dbh,\%parameters,\%quantifInfo,\%quantifValues,\%proteinInfo,\%selectedProteins);
	%quantifAllNames=&promsQuantif::getDistinctQuantifNames($dbh,$quantifList);
	####>Find quantif of modification<####
	my %quantifPTM;
	&findModificationQuantifs(\@selectedQuantifications,\%quantifInfo,\%quantifPTM); #updates %quantifPTM
	###if ($selQuantifFamily eq 'EMPAI') {
	###	my $sthGetPepNum=$dbh->prepare("SELECT ID_PROTEIN,NUM_PEP FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=?");
	###	foreach my $quantiID (keys %quantifInfo) {
	###		$sthGetPepNum->execute($quantifInfo{$quantiID}[4]);
	###		while (my ($protID,$numPep)=$sthGetPepNum->fetchrow_array ) {
	###			$pepInfo{$quantiID}{$protID}=$numPep;
	###		}
	###	}
	###}


	$dbh->disconnect;

	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);

	&printProteinList($call,\@selectedQuantifications,\%quantifInfo,\%quantifValues,\%proteinInfo,\%quantifAllNames,\%quantifPTM);

	exit;
}


####>Revision history<####
# 1.6.5 Option to exclude list of proteins & [Fix]Minor display bug (PP 17/04/19)
# 1.6.4 DB +/-inf ratio value now fetched from promsQuantif.pm (PP 24/01/19)
# 1.6.3 Various improvements for EMPAI (PP 30/10/18)
# 1.6.2 Desactivate no-scroll option in FRAMESET designeb by &generateFrames fuinction (GA 29/06/18)
# 1.6.1 heatMap.js always loaded because needed by ajaxPepSwath in case of DIA quantif (PP 08/02/18)
# 1.6.0 New correlation matrix view (PP 06/02/18)
# 1.5.1 [FIX] Bug to allow modification focus selection even when to quantification family is not known (PP 26/01/18)
# 1.5.0 Improved protein/site quantifs comparison (PP 16/11/17)
# 1.4.3 Check software for SWATH or DIA for peptide ajax call (PP 10/11/17)
# 1.4.2 Bug fix in ajaxListProt call & updated Excel colors (PP 22/02/17)
# 1.4.1 Update emPAI type for quantifications (GA 31/01/17)
# 1.4.0 Handles MaxQuant and SWATH quantifications (PP 20/01/17)
# 1.3.6 Exclude MaxQuant from list of quantif methods (PP 16/12/16)
# 1.3.5 Compatible with ajaxPepSwath call for displaying transitions (PP 14/12/16)
# 1.3.4 Excludes SSPA from list of quantif methods & replaces p-values=0 for RATIO quatifs (PP 14/09/16)
# 1.3.3 Minor modif (GA 01/08/16)
# 1.3.2 Popup info for Abs; fold change (PP 13/04/16)
# 1.3.1 Update compare to make emPAI and MaxQuant compare (GA 29/02/16)
# 1.3.0 PTM quantification & heatmap view removed (PP 27/01/16)
# 1.2.5 JS ajaxManageSaveProteins moved to promsMod.pm (PP 22/08/14)
# 1.2.4 Requires heatMap.js v1.5.2 & moved &fetchQuantificationData to promsQuantif.pm (PP ../08/14)
# 1.2.3 Requires promsQuantif.pm (PP 27/06/14)
# 1.2.2 Sub-list export (PP 27/06/14)
# 1.2.1 Sub-list selection from Venn Diagram (PP 26/06/14)
# 1.2.0 Full SuperRatio support (PP 19/06/14)
# 1.1.4 &infin; ratios to heatmap, log2 views, venn diagram & data export (PP 02/06/14)
# 1.1.3 Fix for miss-spelling of heatMap.js (PP 23/04/14)
# 1.1.2 Restrict to protein quantification comparison & Full Log2 plot<BR>TODO: Handling of SIN & EMPAI (PP 08/04/14)
# 1.1.1 Add one option in compare quantifications for PD vs MCQ comparisons at peptide level : only for dev version ? (GA 18/03/14)
# 1.1.0 Minor corrections (PP 27/02/14)
# 1.0.9 Export SVG option (PP 12/02/14)
# 1.0.8 ajax call in JS ajaxProteinData() deported to showProtQuantification.cgi (PP 08/10/13)
# 1.0.7 Compatible with design SILAC & +/- infinite ratios (PP 03/09/13)
# 1.0.6 Compatible with label-free quantification (PP 03/07/13)
# 1.0.5 Minor big correction (PP 31/05/13)
# 1.0.4 -charset=>'utf-8' (PP 10/05/13)
# 1.0.3 Minor modification, related to neg PEP_BEG for Ghost Peptides (GA 25/02/13)
# 1.0.2 Prevents num peptide set to 1 when view is volcano (PP 29/01/13)
# 1.0.1 Requires heatMap.js and no longer heatMap2.js (PP 11/09/12)
# 1.0.0 First production version. Only for labeled quantifications (PP 30/07/12)
