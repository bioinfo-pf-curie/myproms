#!/usr/local/bin/perl -w

################################################################################
# startDesignProtRulerQuantif.cgi      1.0.8                                   #
# Authors: V. Laigle (Institut Curie)                                          #
# Contact: myproms@curie.fr                                                    #
################################################################################
#----------------------------------CeCILL License-------------------------------
# This file is part of myProMS
#
# Copyright Institut Curie 2019
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
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use promsQuantif;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};

###############################
####>Recovering parameters<####
###############################
my $designID=&promsMod::cleanNumericalParameters(param('ID'));
my $action=param('ACT') || '';
my $intParamID=param('intMetric') || '';

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;
my $projectID=&promsMod::getProjectID($dbh,$designID,'design');
#my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
#my $bioProjectAccessStrg=''; # for template selection
#if ($userInfo[1]=~/^(bio|manag)$/) { # for template selection
#	@userInfo=&promsMod::getUserInfo($dbh,$userID); # scan all accessible projects
#	$bioProjectAccessStrg=join(',',keys %{$userInfo[2]});
#}
my ($experimentID,$designName)=$dbh->selectrow_array("SELECT ID_EXPERIMENT,NAME FROM DESIGN WHERE ID_DESIGN=$designID");
my $titleString="Start Proteomic Ruler Quantification from Design <FONT color=#DD0000>$designName</FONT>";

################################
####> Getting info from DB <####
################################

my $quantifMethod = 'PROT_RULER';

# For now, Proteomic Ruler based only on MaxQuant quantifs with Intensity and/or LFQ data available
my %parentMethodParam = (
	'MQ' => {
		'params' => {
			'MQ_INT' => {},
			'MQ_LFQ' => {}
		}
	},
	'PROT_ABUNDANCE' => {
		'params' => {
			'MY_LFQ' => {}
		}
	}#,
#	'PROT_RATIO_PEP' => {
#		'params' => {
#			'MEAN_STATE' => {}
#		}
#	}
);
=begin comment
Final structure of the hash : 
key 1 = parent quantification method CODE
	key 2a = method ID
	key 2b = method Name
	key 2c = params
		key 3 = param CODE
			key 4a = param ID
			key 4b = param Name

=end comment
=cut

# Retrieve ID and name for each Quantification Method Parameter (parent quantif + current quantif)
my $sthGetMethod = $dbh->prepare("SELECT ID_QUANTIFICATION_METHOD, NAME FROM QUANTIFICATION_METHOD WHERE CODE=?");
my $sthParam 	= $dbh->prepare("SELECT ID_QUANTIF_PARAMETER, NAME FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=? AND CODE=?");

$sthGetMethod->execute($quantifMethod);
my ($quantifMethodID, $quantifMethodName) = $sthGetMethod->fetchrow_array;

foreach my $parentMethodCode (keys %parentMethodParam) {
	$sthGetMethod->execute($parentMethodCode);
	my ($parentMethodID, $parentMethodName) = $sthGetMethod->fetchrow_array;
	$parentMethodParam{$parentMethodCode}->{'methodID'} = $parentMethodID;
	$parentMethodParam{$parentMethodCode}->{'methodName'} = $parentMethodName;

	foreach my $parentParamCode (keys %{$parentMethodParam{$parentMethodCode}->{'params'}}) {
		$sthParam->execute($parentMethodID, $parentParamCode);
		my ($parentParamID, $parentParamName) = $sthParam->fetchrow_array;
		$parentMethodParam{$parentMethodCode}->{'params'}{$parentParamCode}{'paramID'}   = $parentParamID;
		$parentMethodParam{$parentMethodCode}->{'params'}{$parentParamCode}{'paramName'} = $parentParamName;
	}
}

$sthGetMethod->finish;
$sthParam->finish;

# Get all available quantif params for PROT_RULER quantif
my @quantifParamCodes;
my %quantifParams;
my $sthQuantifMethodParam=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER, NAME, CODE FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=$quantifMethodID ORDER BY ID_QUANTIF_PARAMETER");
$sthQuantifMethodParam->execute;
while (my ($quantifParamID, $quantifParamName, $quantifParamCode)=$sthQuantifMethodParam->fetchrow_array) {
	push @quantifParamCodes, $quantifParamCode;
	$quantifParams{$quantifParamCode}=[$quantifParamID,$quantifParamName];
}

####>Custom lists (Categories)<####
my %categoryList;
my $sthCat=$dbh->prepare("SELECT CL.ID_CLASSIFICATION,CL.NAME,ID_CATEGORY,CA.NAME,CA.LIST_TYPE FROM CLASSIFICATION CL,CATEGORY CA WHERE CL.ID_CLASSIFICATION=CA.ID_CLASSIFICATION AND ID_PROJECT=$projectID ORDER BY CL.NAME ASC,CA.DISPLAY_POS ASC");
$sthCat->execute;
while (my ($classID,$className,$catID,$catName,$listType)=$sthCat->fetchrow_array) {
	$listType='PROT' unless $listType;
	$categoryList{$classID}{'CL_NAME'}=$className unless $categoryList{$classID};
	push @{$categoryList{$classID}{'CAT'}},[$catID,$catName,$listType];
}
$sthCat->finish;

$dbh->disconnect;

if ($action eq 'ajaxDisplayQuantifs'){&ajaxDisplayQuantifs; exit;}

############################
####>Form was submitted<####
############################
if (param('launch')) {
	launchQuantifications();
	exit;
}

#######################
####>Starting HTML<####
#######################
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
print header(-charset=>'utf-8'); # because of ° in source selection
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Configure Prot Ruler</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TD.center {text-align:center}
</STYLE>
<SCRIPT type="text/javascript">
|;

print qq
|
// AJAX functions ------> //
function getXMLHTTP() {
	var xhr=null;
	if(window.XMLHttpRequest) {  // Firefox & others //
		xhr = new XMLHttpRequest();
	}
	else if(window.ActiveXObject){  // Internet Explorer //
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
	else {  // XMLHttpRequest not supported by browser //
		alert("Your browser does not support XMLHTTPRequest objects...");
	}
	return xhr;
}

function ajaxDisplayQuantifs(srcType,srcValue) {
	var quantifsCell = document.getElementById('QUANTIFS_CELL');
	if (srcType=='intMetric') {
		// Creation of the XMLHTTPRequest object //
		var XHR = getXMLHTTP();
		if (!XHR) {
			return false;
		}
		XHR.open("GET","$promsPath{cgi}/startDesignProtRulerQuantif.cgi?ACT=ajaxDisplayQuantifs&ID=$designID&intMetric="+srcValue,true);
		XHR.onreadystatechange=function() {
			if (XHR.readyState==4 && XHR.responseText) {
				if (srcValue) {
					quantifsCell.innerHTML=(XHR.responseText);
					displayFields('averagingMode', document.getElementById('AVG_MODE').value);
					displayFields('protRulerType', document.getElementById('PROT_RULER').value);
					displayFields('quantifScale', document.getElementById('QUANTIF_SCALE').value);
				} else {
					displayFields('averagingMode', "0");
					displayFields('protRulerType', "");
					displayFields('quantifScale', "sample");
					quantifsCell.innerHTML=(XHR.responseText);
				}
			}
		}
		XHR.send(null);
	}
}
// <------ End of AJAX functions //

function displayFields(srcType,srcValue) {
	var quantifBoxes = document.getElementsByName('quantifBoxes');
	var count;
	if (srcType=='averagingMode') {
		var groupSelect = document.getElementsByName('groups');
		if (srcValue=="2") {
			hideElement('grpHidden', isID=true, disable=false);
			showElement('grpLabel', isID=true, enable=false);
			count = 0;
			for (let groupSel of groupSelect) {
				if (quantifBoxes[count].checked) {
					showElement(groupSel, isID=false, enable=true);
				} else {
					hideElement(groupSel, isID=false, disable=true);
				}
				count++;
			}
		} else {
			hideElement('grpLabel', isID=true, disable=false);
			showElement('grpHidden', isID=true, enable=false);
			for (let groupSel of groupSelect) {
				hideElement(groupSel, isID=false, disable=true);
			}
		}
	} else if (srcType=='protRulerType') {
		var customQtySelect = document.getElementsByName('customQty');
		var askCellNb = (document.getElementById('QUANTIF_SCALE').value == 'cell')? false : true;
		var cellNbSelect = document.getElementsByName('cellNb');
		
		// Hide and disable everything, a bit less efficient but it avoids redondancy of hiding everything in each case
		hideElement('protRulerParams', isID=true, disable=false);
		hideElement('totalProteinAmount', isID=true, disable=true);
		hideElement('tpaLabel', isID=true, disable=false);
		hideElement('tpaCellNb', isID=true, disable=true);
		hideElement('tpaCellNbLabel', isID=true, disable=false);
		hideElement('histoneProteomicRuler', isID=true, disable=true);
		hideElement('hprLabel', isID=true, disable=false);
		hideElement('customProteins', isID=true, disable=false);
		hideElement('customProtSelection', isID=true, disable=true);
		hideElement('customProtFile', isID=true, disable=true);
		hideElement('customProtLabel', isID=true, disable=false);
		hideElement('customQtyLabel', isID=true, disable=false);
		hideElement('customQtyHidden', isID=true, disable=false);
		hideElement('cellNbLabel', isID=true, disable=false);
		hideElement('cellNbHidden', isID=true, disable=false);
		for (let customQtySel of customQtySelect) {
			hideElement(customQtySel, isID=false, disable=true);
		}
		for (let cellNbSel of cellNbSelect) {
			hideElement(cellNbSel, isID=false, disable=true);
		}
		// Show and enable what is needed in each case
		if (srcValue=="0") {
			showElement('protRulerParams', isID=true, enable=false);
			showElement('totalProteinAmount', isID=true, enable=true);
			showElement('tpaLabel', isID=true, enable=false);
			if (askCellNb) {
				showElement('tpaCellNb', isID=true, enable=true);
				showElement('tpaCellNbLabel', isID=true, enable=false);
			}
			showElement('customQtyHidden', isID=true, enable=false);
			showElement('cellNbHidden', isID=true, enable=false);
		} else if (srcValue=="1") {
			showElement('protRulerParams', isID=true, enable=false);
			showElement('histoneProteomicRuler', isID=true, enable=true);
			showElement('hprLabel', isID=true, enable=false);
			showElement('customQtyHidden', isID=true, enable=false);
			showElement('cellNbHidden', isID=true, enable=false);
		} else if (srcValue=="2") {
			showElement('protRulerParams', isID=true, enable=false);
			showElement('customProteins', isID=true, enable=false);
			showElement('customProtSelection', isID=true, enable=true);
			showElement('customProtFile', isID=true, enable=true);
			showElement('customProtLabel', isID=true, enable=false);
			showElement('customQtyLabel', isID=true, enable=false);
			count = 0;
			for (let customQtySel of customQtySelect) {
				if (quantifBoxes[count].checked) {
					showElement(customQtySel, isID=false, enable=true);
				} else {
					hideElement(customQtySel, isID=false, disable=true);
				}
				count++;
			}
			if (askCellNb) {
				showElement('cellNbLabel', isID=true, enable=false);
				count = 0;
				for (let cellNbSel of cellNbSelect) {
					if (quantifBoxes[count].checked) {
						showElement(cellNbSel, isID=false, enable=true);
					} else {
						hideElement(cellNbSel, isID=false, disable=true);
					}
					count++;
				}
			} else {
				showElement('cellNbHidden', isID=true, enable=false);
			}
		} else {
			showElement('customQtyHidden', isID=true, enable=false);
			showElement('cellNbHidden', isID=true, enable=false);
		}
	} else if (srcType == 'quantifScale') {
		var protRulerType = document.getElementById('PROT_RULER').value;
		var cellNbSelect = document.getElementsByName('cellNb');
		var customQtyLabel = document.getElementById('customQtyLabel');
		var tpaLabel = document.getElementById('tpaLabel');
		var headerEnd = '&nbsp;<BR>&nbsp;for each sample (pg)&nbsp;<BR><INPUT type="button" name="protQtyButton" id="protQtyButton" value="Same for all" onclick="reportValues(\\\'protQty\\\');">';
		
		if (srcValue == "cell") {
			tpaLabel.innerHTML = '&nbsp;Total protein quantity PER CELL (pg):&nbsp;';
			customQtyLabel.innerHTML = '&nbsp;Custom proteins quantities PER CELL' + headerEnd;
			for (let cellNbSel of cellNbSelect) {
				hideElement(cellNbSel, isID=false, disable=true);
			}
			hideElement('cellNbLabel', isID=true, disable=false);
			showElement('cellNbHidden', isID=true, enable=false);
			hideElement('tpaCellNb', isID=true, disable=true);
			hideElement('tpaCellNbLabel', isID=true, disable=false);
		} else {
			tpaLabel.innerHTML = '&nbsp;Total protein quantity per sample (pg):&nbsp;';
			customQtyLabel.innerHTML = '&nbsp;Custom proteins quantities' + headerEnd;
			
			if (protRulerType == "0") {
				showElement('tpaCellNb', isID=true, enable=true);
				showElement('tpaCellNbLabel', isID=true, enable=false);
				for (let cellNbSel of cellNbSelect) {
					hideElement(cellNbSel, isID=false, disable=true);
				}
				hideElement('cellNbLabel', isID=true, disable=false);
				showElement('cellNbHidden', isID=true, enable=false);
			} else if (protRulerType == "2") {
				count = 0;
				for (let cellNbSel of cellNbSelect) {
					if (quantifBoxes[count].checked) {
						showElement(cellNbSel, isID=false, enable=true);
					} else {
						hideElement(cellNbSel, isID=false, disable=true);
					}
					count++;
				}
				hideElement('cellNbHidden', isID=true, disable=false);
				showElement('cellNbLabel', isID=true, enable=false);
				hideElement('tpaCellNb', isID=true, disable=true);
				hideElement('tpaCellNbLabel', isID=true, disable=false);
			} else {
				for (let cellNbSel of cellNbSelect) {
					hideElement(cellNbSel, isID=false, disable=true);
				}
				hideElement('cellNbLabel', isID=true, disable=false);
				showElement('cellNbHidden', isID=true, enable=false);
				hideElement('tpaCellNb', isID=true, disable=true);
				hideElement('tpaCellNbLabel', isID=true, disable=false);
			}
		}
	}
}

function checkUncheckBox (myBox, boxFeature) {
	var myForm=document.setParamForm;
	var boxValue=myBox.checked;
	if (boxFeature == 'correction') {
		myForm.correctionFactor.disabled = (!boxValue)? true : false;
	} else if (boxFeature == 'quantif') {
		var boxNb = myBox.dataset.boxNb;
		var displayCustProt = (document.getElementById('PROT_RULER').value == '2')? true : false;
		var displayCellNb = (displayCustProt && document.getElementById('QUANTIF_SCALE').value != 'cell')? true : false;
		var displayGroups = (document.getElementById('AVG_MODE').value == '2')? true : false;

		if (boxValue) {
			if (displayCellNb) {  // Implies displayCustProt is true
				showElement(document.getElementsByName('customQty')[boxNb-1], isID=false, enable=true);
				showElement(document.getElementsByName('cellNb')[boxNb-1], isID=false, enable=true);
			} else if (displayCustProt) {
				showElement(document.getElementsByName('customQty')[boxNb-1], isID=false, enable=true);
				hideElement(document.getElementsByName('cellNb')[boxNb-1], isID=false, disable=true);
			} else {
				hideElement(document.getElementsByName('customQty')[boxNb-1], isID=false, disable=true);
				hideElement(document.getElementsByName('cellNb')[boxNb-1], isID=false, disable=true);
			}
			if (displayGroups) {
				showElement(document.getElementsByName('groups')[boxNb-1], isID=false, enable=true);
			} else {
				hideElement(document.getElementsByName('groups')[boxNb-1], isID=false, disable=true);
			}
		} else {
			hideElement(document.getElementsByName('customQty')[boxNb-1], isID=false, disable=true);
			hideElement(document.getElementsByName('cellNb')[boxNb-1], isID=false, disable=true);
			hideElement(document.getElementsByName('groups')[boxNb-1], isID=false, disable=true);
		}
	}
}

function reportValues(srcType) {
	var selection;
	if (srcType == 'cellNb') {
		selection = document.getElementsByName('cellNb');
	} else if (srcType == 'protQty') {
		selection = document.getElementsByName('customQty');
	}
	var count = 0;
	var valueToReport = 0;
	for (let cell of selection) {
		if (!valueToReport && !cell.disabled && cell.value) {
			valueToReport = cell.value;
			break;
		}
	}
	if (valueToReport) {
		for (let cell of selection) {
			if (!cell.disabled) {
				cell.value = valueToReport;
			}
		}
	}
}

function hideElement(elem, isID=false, disable=false) {
	var element;
	if (isID) {
		element = document.getElementById(elem);
	} else {
		element = elem;
	}
	element.style.display = 'none';
	if (disable) {
		element.disabled = true;
	}
}

function showElement(elem, isID=false, enable=false) {
	var element;
	if (isID) {
		element = document.getElementById(elem);
	} else {
		element = elem;
	}
	element.style.display = '';
	if (enable) {
		element.disabled = false;
	}
}

function cancelAction() {
	top.promsFrame.optionFrame.selectOption();
}

function checkForm(myForm) {
	var allQuantifs=document.getElementsByName('quantifBoxes');
	
	// name //
	if (!myForm.quantifName.value) {
		alert('ERROR: Missing name for quantification !');
		return false;
	}
	
	// prot. ruler type //
	if (!myForm.protRuler.value) {
		alert('ERROR: Choose a method for the Proteomic Ruler !');
		return false;
	} else if (myForm.protRuler.value=="0") {
		if (!myForm.totalProteinAmount.value){
			alert('ERROR: You chose the method "Total Protein Amount", please provide the total protein mass per sample (in pg, it is assumed to be the same in every sample).');
			return false;
		} else if (myForm.quantifScale.value != "cell" && (!myForm.tpaCellNb.value \|\| myForm.tpaCellNb.value == 0)) {
			alert('ERROR: You chose to provide protein quantities per sample + the number of cells, please provide the number of cells (it is assumed to be the same in every sample).');
			return false;
		}
	} else if (myForm.protRuler.value=="1") {
		if (!myForm.histoneProteomicRuler.value){
			alert('ERROR: You chose the method "Histone Proteomic Ruler" which needs the ploidy of your cell type, please provide it !');
			return false;
		}
	} else if (myForm.protRuler.value=="2") {
		if (!myForm.customProtSelection.value && !myForm.customProtFile.value) {
			alert('ERROR: You chose the method "Custom Proteins Ruler" which requires a list of your custom proteins (those of known quantity). Either select the corresponding list or provide a text file containing the Uniprot ACC of these proteins !');
			return false;
		} else if (myForm.customProtSelection.value && myForm.customProtFile.value) {
			alert('ERROR: You provided both a list and a file for your custom proteins, please provide only one or the other !');
			return false;
		}
	}

	// Quantification scale //
	if (!myForm.quantifScale.value) {
		alert('ERROR: Select at which scale (sample or cell) your values are set !');
		return false;
	}

	// Intensity metric selection //
	if (!myForm.intMetric.value){
		alert('ERROR: Please select a type of Intensity metric ! (If none of them is available it may be because you did not get Intensity measures from your previous quantifications)');
		return false;
	}

	// selected quantifs //
	var selectedQuantifs = [];
	for (let quantifBox of allQuantifs){
		if (quantifBox.checked){
			selectedQuantifs.push(quantifBox.dataset.boxNb);
		}
	}

	if (selectedQuantifs.length == 0){
		alert('ERROR: Select at least one quantification/sample to perform Proteomic Ruler on !');
		return false;
	}
	
	// custom prot quantity provided ? //
	if (myForm.protRuler.value=="2") {
		for (let i of selectedQuantifs) {
			let customQtyId='CUSTOM_QTY_' + i;
			if (!document.getElementById(customQtyId).value \|\| document.getElementById(customQtyId).value == 0){
				alert('ERROR: You chose the method "Custom Proteins Ruler", please provide the quantity of your custom proteins for every selected quantification/sample !');
				return false;
			}
		}
	}

	// cell number provided ? //
	if (myForm.protRuler.value == "2" && myForm.quantifScale.value != "cell") {
		for (let i of selectedQuantifs) {
			let cellNbId='CELL_NB_' + i;
			if (!document.getElementById(cellNbId).value \|\| document.getElementById(cellNbId).value == 0){
				alert('ERROR: You chose to provide protein quantities per sample + the number of cells, please provide the number of cells for every selected quantification/sample !');
				return false;
			}
		}
	}
	
	// groups selection //
	if (myForm.averagingMode.value=="2") {
		for (let i of selectedQuantifs) {
			let groupId='GROUP_' + i;
			if (!document.getElementById(groupId).value){
				alert('ERROR: You chose to normalize by groups, please provide the group for every selected quantification/sample !');
				return false;
			}
		}
	}
	
	// total protein concentration //
	if (!myForm.proteinConcentration.value) {
		alert('ERROR: Provide the total cellular protein concentration of your samples !');
		return false;
	}

	return true;
}
|;

print qq
|
</SCRIPT>
</HEAD>

<BODY background="$promsPath{images}/bgProMS.gif" >
<CENTER>
<FONT class="title">$titleString</FONT>
<BR><BR><BR>
<FORM name="setParamForm" id="PARAM_FORM" method="post" enctype="multipart/form-data" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="ID" value="$designID">
<TABLE border="0" bgcolor="$darkColor">
	<!--Quantification name-->
	<TR>
		<TH class="title3" align=right valign=top>&nbsp;Name:&nbsp;</TH>
		<TD colspan=4 bgcolor="$lightColor" style="width:700px">
		<INPUT type="text" name="quantifName" value="" placeholder="Name of quantification" class="title3" style="width:350px"/>
		</TD>
	</TR>
	
	<!--Proteomic Ruler Type-->
	<TR>	
		<TH align=right nowrap>&nbsp;Proteomic Ruler:&nbsp;</TH>
		<TD colspan=4 bgcolor="$lightColor">
		<SELECT name="protRuler" id="PROT_RULER" onchange="if (document.getElementById('intMetric').value) {displayFields('protRulerType',this.value);}" value="">
			<OPTION value="" selected>-= Select =-</OPTION>
			<OPTION value="0">Total Protein Approach</OPTION>
			<OPTION value="1">Histone Proteomic Ruler</OPTION>
			<OPTION value="2">Custom Proteins Ruler</OPTION>
		</SELECT>
		</TD>
	</TR>
	
	<!--Absolute quantification per cell or per sample-->
	<TR>
		<TH align=right nowrap>&nbsp;Type of values you will enter below:&nbsp;</TH>
		<TD colspan=4 bgcolor="$lightColor">
		<SELECT name="quantifScale" id="QUANTIF_SCALE" onchange="if (document.getElementById('intMetric').value) {displayFields('quantifScale', this.value);}" value="sample">
			<OPTION value="sample" selected>Quantity per sample + cell number</OPTION>
			<OPTION value="cell">Quantity per cell</OPTION>
		</SELECT>
		</TD>
	</TR>

	<!--Averaging mode-->
	<TR>	
		<TH align=right nowrap>&nbsp;Normalization mode:&nbsp;</TH>
		<TD colspan=4 bgcolor="$lightColor">
		<SELECT name="averagingMode" id="AVG_MODE" onchange="if (document.getElementById('intMetric').value) {displayFields('averagingMode', this.value);}" value="0">
			<OPTION value="0" selected>All samples separately</OPTION>
			<OPTION value="1">Same normalization for all</OPTION>
			<OPTION value="2">Same normalization within groups</OPTION>
			<OPTION value="3">Average all samples</OPTION>
		</SELECT>
		</TD>
	</TR>

	<!--Intensity Metric-->
	<TR>
		<TH align=right nowrap>
		<LABEL for="intMetric" id="metricLabel">&nbsp;Intensity metric:&nbsp;</LABEL>
		</TH>
		<TD colspan=4 bgcolor="$lightColor">
		<SELECT name="intMetric" id="intMetric" onchange="ajaxDisplayQuantifs('intMetric', this.value);">
			<OPTION value="" selected>-= Select =-</OPTION>
|;

foreach my $parentMethodCode (sort {$parentMethodParam{$a}->{'methodID'} <=> 
									$parentMethodParam{$b}->{'methodID'}
								   } keys %parentMethodParam) {
	print qq
|			<optgroup label="$parentMethodParam{$parentMethodCode}->{'methodName'} :">\n
|;
	foreach my $parentParamCode (sort { $parentMethodParam{$parentMethodCode}->{'params'}{$a}{'paramID'} <=> 
										$parentMethodParam{$parentMethodCode}->{'params'}{$b}{'paramID'}
									  } keys %{$parentMethodParam{$parentMethodCode}->{'params'}}) {
		print qq
|				<option value="$parentMethodParam{$parentMethodCode}->{'params'}{$parentParamCode}{'paramID'}">$parentMethodParam{$parentMethodCode}->{'params'}{$parentParamCode}{'paramName'}</option>\n
|;
	}
	print qq
|			</optgroup>\n
|;
}

print qq
|		</SELECT>
		</TD>
	</TR>
	
	<BR>
	<TR>
		<TD id=QUANTIFS_CELL colspan=5 bgcolor="$lightColor" style="text-align:center">
		<BR>Select <b>Intensity metric</b> to display available quantifications<BR>&nbsp;
		</TD>
	</TR>
	<BR>
	<!--Proteomic Ruler Parameters-->
	<TR id=protRulerParams style="display:none">
		<TH align=right nowrap>
		<LABEL for="totalProteinAmount" id="tpaLabel" style="display:none">&nbsp;Total protein quantity (pg):&nbsp;</LABEL>
		<LABEL for="histoneProteomicRuler" id="hprLabel" style="display:none">&nbsp;Ploidy of the cells:&nbsp;</LABEL>
		<LABEL for="customProteins" id="customProtLabel" style="display:none">&nbsp;Custom proteins (Uniprot ACC):&nbsp;</LABEL>
		</TH>
		<TD colspan=4 bgcolor="$lightColor">
		<INPUT type="number" min=0 step="any" name="totalProteinAmount" id="totalProteinAmount" placeholder="200" style="display:none">
		<LABEL for='tpaCellNb', id='tpaCellNbLabel' style="display:none">&nbsp;&nbsp;&nbsp;&nbsp;Number of cells per sample</LABEL>
		<INPUT type="number" min=0 step="any" name="tpaCellNb" id="tpaCellNb" placeholder="10000" style="display:none">
		<INPUT type="number" min=0 step="any" name="histoneProteomicRuler" id="histoneProteomicRuler" placeholder="2" style="display:none">		
		<DIV id="customProteins" style="display:none">
			<B>&nbsp;Choose a List:&nbsp;</B>
			<SELECT name="customProtSelection" id="customProtSelection" class="template">
				<OPTION value="">-= Select =-</OPTION>
|;
foreach my $classID (sort{lc($categoryList{$a}{'CL_NAME'}) cmp lc($categoryList{$b}{'CL_NAME'})} keys %categoryList) {
	print "<OPTGROUP label=\"$categoryList{$classID}{CL_NAME}:\">\n";
	foreach my $refCat (@{$categoryList{$classID}{'CAT'}}) {
		my $typeStrg=($refCat->[2] eq 'SITE')? ' [Sites]' : '';
		print '<OPTION value="',$refCat->[0],'">',$refCat->[1],"$typeStrg</OPTION>\n";
	}
	print "</OPTGROUP>\n";
}
print qq
|			</SELECT>
			<B>&nbsp;or Upload a file:&nbsp;</B>
			<INPUT type="file" name="customProtFile" id="customProtFile" accept=".txt,.doc,.docx,.xml,application/msword,application/vnd.openxmlformats-officedocument.wordprocessingml.document">
		</DIV>
		</TD>
	</TR>
	
	<!--Detectability correction-->
	<!--<TR>-->
	<!--	<TH align=right nowrap>&nbsp;Detectability correction:&nbsp;</TH>-->
	<!--	<TD colspan=4 bgcolor="$lightColor">-->
	<!--	<INPUT type="checkbox" name="correctionBox" onclick="checkUncheckBox(this, 'correction')">-->
	<!--	<LABEL for="CORR_FACTOR">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Correction factor:&nbsp;</LABEL>-->
	<!--	<INPUT type="text" name="correctionFactor" id="CORR_FACTOR" minlength="0" maxlength="25" size="26" placeholder="Number of theoretical peptides, ..." disabled>-->
	<!--	</TD>-->
	<!--</TR>-->

	<!--Total protein concentration-->
	<TR>	
		<TH align=right nowrap>&nbsp;Total cellular protein concentration (g/L):&nbsp;</TH>
		<TD colspan=4 bgcolor="$lightColor">
		<INPUT type="number" min=0 step="any" name="proteinConcentration" id="PROT_CONC" placeholder="200">
		</TD>
	</TR>
	
	<!--Protein selection (exclusion or restriction to a given list)-->
	<TR>
		<TH align=right nowrap>&nbsp;Protein selection:&nbsp;</TH>
		<TD colspan=4 bgcolor=$lightColor nowrap>
		<SELECT name="protSelType" class="template">
			<OPTION value="exclude">Exclude</OPTION>
			<OPTION value="restrict">Restrict to</OPTION>
		</SELECT>&nbsp;
		<B>proteins from List:&nbsp;</B>
		<SELECT name="protSelection" class="template">
			<OPTION value="">-= Select =-</OPTION>
|;
foreach my $classID (sort{lc($categoryList{$a}{'CL_NAME'}) cmp lc($categoryList{$b}{'CL_NAME'})} keys %categoryList) {
	print "<OPTGROUP label=\"$categoryList{$classID}{CL_NAME}:\">\n";
	foreach my $refCat (@{$categoryList{$classID}{'CAT'}}) {
		my $typeStrg=($refCat->[2] eq 'SITE')? ' [Sites]' : '';
		print '<OPTION value="',$refCat->[0],'">',$refCat->[1],"$typeStrg</OPTION>\n";
	}
	print "</OPTGROUP>\n";
}
print qq
|		</SELECT>
		</TD>
	</TR>

	<!--Organism name-->
	<TR>	
		<TH align=right nowrap>&nbsp;Organism name:&nbsp;</TH>
		<TD colspan=4 bgcolor="$lightColor">
		<SELECT name="organismName" id="ORG_NAME">
			<OPTION value="" selected>-= Select =-</OPTION>
			<OPTION value="Homo_sapiens">Homo sapiens</OPTION>
			<OPTION value="Mus_musculus">Mus musculus</OPTION>
			<OPTION value="Drosophila_melanogaster">Drosophila melanogaster</OPTION>
			<OPTION value="Caenorhabditis_elegans">Caenorhabditis elegans</OPTION>
			<OPTION value="Saccharomyces_cerevisiae">Saccharomyces cerevisiae</OPTION>
			<OPTION value="Schizosaccharomyces_pombe">Schizosaccharomyces pombe</OPTION>
			<OPTION value="Gallus_gallus">Gallus gallus</OPTION>
			<OPTION value="">Other</OPTION>
		</SELECT>
		</TD>
	</TR>
	
	<TR>	
		<!--Type of output-->
		<TH align=right nowrap>&nbsp;Desired output(s):&nbsp;</TH>
		<TD colspan=4 bgcolor="$lightColor">
		<SELECT name="output" id="OUTPUT" multiple>
|;
foreach my $quantifParamCode (@quantifParamCodes){
	print qq
|			<OPTION name="quantifParam" value="$quantifParams{$quantifParamCode}->[0]">$quantifParams{$quantifParamCode}->[1]</OPTION>
|;
}
print qq
|		</SELECT>
		</TD>
	</TR>
</TABLE>
<BR>
<INPUT type="submit" name="launch" value="Launch Quantification" class="title3">&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="cancelAction();">

</FORM>
</CENTER>
</BODY>
</HTML>
|;
#######################
#####>End of HTML<#####
#######################


sub launchQuantifications {
	####<Starting HTML>####
	print header(-'content-encoding'=>'no');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Launching Quantifications</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER><FONT class="title1">Launching Quantification(s)</FONT></CENTER>
<BR><BR><BR>
|;

	my ($projectID,$experimentID); # defined later for Design. needed for nav frame refresh
	
	################################
	####> Getting info from DB <####
	################################
	my $dbh=&promsConfig::dbConnect;
	($experimentID)=$dbh->selectrow_array("SELECT ID_EXPERIMENT FROM DESIGN WHERE ID_DESIGN=$designID");
	$projectID=&promsMod::getProjectID($dbh,$experimentID,'experiment');
	
	###########################################
	####>Getting parameters from HTML form<####
	###########################################
	my $protRulerType=param('protRuler');
	my $quantifScale=param('quantifScale');
	my $averagingMode=param('averagingMode');
	my $intMetric=param('intMetric');  # ID_QUANTIF_PARAMETER corresponding to the type of intensity used
	my @quantifBoxes=param("quantifBoxes");
	my (%availableStates, %statesPos);
	getAvailableStates(\%availableStates, $intMetric);
	
	# Compute target_pos for all available states
	my $sthQuantifAnnot = $dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
	my $sthTargetPos=$dbh->prepare("SELECT TARGET_POS FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=? AND ID_QUANTIF_PARAMETER=$intMetric GROUP BY TARGET_POS ORDER BY TARGET_POS ASC");

	foreach my $parentQuantifID (sort {$a <=> $b} keys %availableStates) {
		$sthQuantifAnnot->execute($parentQuantifID);
		my ($quantifAnnot) = $sthQuantifAnnot->fetchrow_array;
		my @annotFields = split(/::/, $quantifAnnot);
		my $statesAnnot;
		foreach my $annotField (@annotFields) {
			if ($annotField =~ /STATES=(.*)/) {
				$statesAnnot = $1;
				last;
			}
		}
		if ($statesAnnot) {
			getAnnotStatesOrder($statesAnnot, \@{$statesPos{$parentQuantifID}});
		} else {
			die "Could not get the STATES information from quantif annot for quantification $parentQuantifID : $!";
		}

		$sthTargetPos->execute($parentQuantifID);

		foreach my $parentStateID (@{$statesPos{$parentQuantifID}}) {
			my ($targetPos) = $sthTargetPos->fetchrow_array;  # Consider multiple targetPos
			$availableStates{$parentQuantifID}->{'states'}{$parentStateID}{'targetPos'} = $targetPos;
		}
	}
	$sthQuantifAnnot->finish;
	$sthTargetPos->finish;
	
	my @selectedStates;
	foreach my $parentQuantifID (sort {$a <=> $b} keys %availableStates) {
		foreach my $parentStateID (sort {$a <=> $b} keys %{$availableStates{$parentQuantifID}->{'states'}}) {
			if (grep(/^$parentQuantifID\_$parentStateID$/, @quantifBoxes)) {
				my $parentQuantifName 	= $availableStates{$parentQuantifID}->{'name'};
				my $parentStateName		= $availableStates{$parentQuantifID}->{'states'}{$parentStateID}{'name'};
				my $parentTargetPos		= $availableStates{$parentQuantifID}->{'states'}{$parentStateID}{'targetPos'};
				
				push @selectedStates, ["$parentQuantifID\_$parentTargetPos", 
									   "$parentQuantifName\_$parentStateName",
									   "$parentQuantifID,$parentStateID,$parentTargetPos"
									  ];
			}
		}
	}

	# Get quantif Method to write it in quantif_info
	my $quantifFamily;
	foreach my $parentMethodCode (keys %parentMethodParam) {	
		foreach my $parentParamCode (keys %{$parentMethodParam{$parentMethodCode}->{'params'}}) {
			if ($intMetric == $parentMethodParam{$parentMethodCode}->{'params'}{$parentParamCode}{'paramID'}) {
				$quantifFamily = $parentMethodCode;
			}
		}
	}
	$quantifFamily = ($quantifFamily eq 'PROT_RATIO_PEP')? 'RATIO:MEAN' : $quantifFamily;

	# Protein list restriction/exclusion ?
	my ($protSelType, $protSelList);
	if (param('protSelection')) {
		$protSelType=param('protSelType');
		$protSelList=param('protSelection');
	}
	
	my @groups;
	if ($averagingMode == 2) {
		@groups=param('groups');
	} else {
		push @groups,'None';
	}
	my $groupsString=join(';', @groups);
	my $totalProteinAmount = param('totalProteinAmount') || 'None';
	if ($totalProteinAmount ne 'None' && $quantifScale ne "cell") {  # tot prot amount per sample + nb of cells provided
		my $tpaCellNb = param('tpaCellNb');
		$totalProteinAmount = $totalProteinAmount / $tpaCellNb;
	}
	my $histoneProteomicRuler=param('histoneProteomicRuler') || 'None';

	# Get custom proteins ACC from list or from uploaded file and create string from it to pass as parameter
	my $customProtCatID=param('customProtSelection');
	my $custProtFile=tmpFileName(upload('customProtFile')) if upload('customProtFile');
	my $customProteins;
	my @custProtArr;

	if ($customProtCatID) {  # Get custom proteins Uniprot ACC as a semi-colon separated string from existing list with ID_CATEGORY
		($customProteins)=$dbh->selectrow_array("SELECT GROUP_CONCAT(MI.VALUE SEPARATOR ';') FROM CATEGORY_PROTEIN CP INNER JOIN PROTEIN P ON CP.ID_PROTEIN=P.ID_PROTEIN INNER JOIN MASTERPROT_IDENTIFIER MI ON P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN WHERE ID_CATEGORY=$customProtCatID AND MI.ID_IDENTIFIER=1 AND MI.RANK=1");
	} elsif ($custProtFile) {  # Get the same output but from user's text file
		open(CUST_PROT_FILE, "$custProtFile") or die "Could not open the custom proteins file $!";
		for my $line (<CUST_PROT_FILE>){
			chomp $line;
			if ($line =~ /\s/) {
				my @newCustomProts=split(/\s/, $line);
				push @custProtArr, @newCustomProts;
			} elsif ($line =~ /;/) {
				my @newCustomProts=split(';', $line);
				push @custProtArr, @newCustomProts;
			} elsif ($line =~ /[\w\d]/) {
				push @custProtArr, $line;
			}
		}
		unless (@custProtArr) {
			die "The custom proteins file was empty or could not be read. Try again and make sure your file contains Uniprot ACCs.";
		}
		$customProteins=join(';', @custProtArr);
	} else {  # Not using custom proteins
		$customProteins='None';
	}
	# Get custom proteins quantities and create corresponding string to pass as parameter
	my @customProtQty;
	if ($protRulerType == 2) {
		if ($quantifScale eq "cell") {
			@customProtQty=param('customQty');  # e.g. (6.5452, 6.5452, 6.5452)
		} else {
			my @qtyPerSample = param('customQty');
			my @cellNbPerSample = param('cellNb');
			for my $i (0..$#qtyPerSample) {
				push @customProtQty, $qtyPerSample[$i] / $cellNbPerSample[$i];
			}
		}
	} else {
		push @customProtQty,'None';
	}

	my $customQtyString=join(';', @customProtQty);
	my $logarithmized=0;
	my $logBase='None';
	my $detectabilityCorrection=(param('correctionBox') && param('correctionBox') eq 'on')? 1 : 0;
	my $correctionFactor=param('correctionFactor') || 'None';
	my $proteinConcentration=param('proteinConcentration');
	my $organismName=param('organismName') || 'None';
	my @quantifParamIDs=param('output');
	
	# Force copy number in output
	if (!(grep(/^$quantifParams{'COPY_NB'}->[0]$/, @quantifParamIDs))) {
		unshift @quantifParamIDs, $quantifParams{'COPY_NB'}->[0];
	}
	my $outputString=join(';',@quantifParamIDs);
	
	# Reformat selectedStates for quantif_info.txt and next scripts
	my $quantifIdsStrg;
	my $quantifNamesStrg;
	my $quantifStatesStrg;

	foreach my $selectedState (@selectedStates) {
		if ($quantifIdsStrg) {
			$quantifIdsStrg		.= ";$selectedState->[0]";
			$quantifNamesStrg 	.= ";$selectedState->[1]";
			$quantifStatesStrg	.= ";$selectedState->[2]";
		} else {
			$quantifIdsStrg		= "$selectedState->[0]";
			$quantifNamesStrg	= "$selectedState->[1]";
			$quantifStatesStrg	= "$selectedState->[2]";
		}
	}

	# Disconnect from DB before writing to quantif_info
	$dbh->disconnect;

	# Start writing to files
	mkdir "$promsPath{tmp}/quantification" unless -e "$promsPath{tmp}/quantification";
	my $currentQuantifDir="$promsPath{tmp}/quantification/current";
	mkdir $currentQuantifDir unless -e $currentQuantifDir;
	
	my $quantifName=param('quantifName');
	my $quantifDate=strftime("%Y%m%d%H%M%S",localtime);
	while (-e "$promsPath{tmp}/quantification/$quantifDate") { # to prevent multi-user launch collision
		sleep 2;
		$quantifDate=strftime("%Y%m%d%H%M%S",localtime);
	}

	print '<BR><FONT class="title3">Launching quantification process';
	print "[job #$quantifDate]...";
	
	my $jobDir=$quantifDate;
	my $quantifDir="$promsPath{tmp}/quantification/$jobDir";
	mkdir $quantifDir || die "ERROR detected: $!";
	
	my $labeling="FREE";
	my $algoType="PROT_RULER";
	
	###> Jobs info & data <###
	open (INFO,">$quantifDir/quantif_info.txt");  # item valueR valueDB
	print INFO "USER=$userID\n";
	print INFO "TYPE=DESIGN\n";
	print INFO "PARAMETERS:\n";
	print INFO "LABEL\t\t$labeling\n";
	print INFO "QUANTIF_NAME\t\t$quantifName\n";
	print INFO "ALGO_TYPE\t\t$algoType\n";
	print INFO "ID_DESIGN\t\t$designID\n";
	print INFO "QUANTIF_IDS\tquantif_ids\t$quantifIdsStrg\n";
	print INFO "QUANTIF_STATES\t\t$quantifStatesStrg\n";
	print INFO "QUANTIF_FAMILY\t\t$quantifFamily\n";

	if ($protSelList) {
		print INFO "PROTEINS\t\t$protSelType\t#$protSelList\n";
	}
	
	###> Data needed mostly by python script <###
	print INFO "\tprot_ruler_type\t$protRulerType\n";
	print INFO "\taveraging_mode\t$averagingMode\n";
	print INFO "\tquantif_names\t$quantifNamesStrg\n";
	print INFO "\tgroups\t$groupsString\n";
	print INFO "\ttotal_protein_amount\t$totalProteinAmount\n";
	print INFO "\thistone_proteomic_ruler\t$histoneProteomicRuler\n";
	print INFO "\tcustom_proteins\t$customProteins\n";
	print INFO "\tcustom_prot_qty\t$customQtyString\n";
	print INFO "\tint_metric\t$intMetric\n";
	print INFO "\tlogarithmized\t$logarithmized\n";
	print INFO "\tlog_base\t$logBase\n";
	print INFO "\tdetectability_correction\t$detectabilityCorrection\n";
	print INFO "\tcorrection_factor\t$correctionFactor\n";
	print INFO "\tprotein_concentration\t$proteinConcentration\n";
	print INFO "\torganism_name\t$organismName\n";
	print INFO "\toutput_params\t$outputString\n";

	open(FLAG,">$currentQuantifDir/$jobDir\_request.flag"); # flag file used by watchQuantif
	print FLAG "#";
	close FLAG;

	close INFO;

	my $fullQuantifType='DESIGN:'.$algoType;  # Just to be consistent with startDesignQuantification.cgi and launchQuantifications.pl
	###<Forking to launch quantifications>###
	my $childPid = fork;
	unless ($childPid) { # child here
		#>Disconnecting from server
		open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
		#open STDOUT, ">>$currentQuantifDir/test_debug.out" or die "Can't open test file for debug: $!";  # For debug
		#open STDERR, ">>$currentQuantifDir/test_debug.out" or die "Can't open test file for debug: $!";  # For debug
		open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
		system "./launchProtRulerQuantif.pl single $ENV{REMOTE_USER} $jobDir $fullQuantifType";
		exit;
	}
	print " Done.</FONT><BR>\n";

	print "<BR><FONT class=\"title3\"><BR>This page will refresh itself in a few seconds.</FONT>\n";

	###>Calling watch popup window<###
	sleep 3;
	print qq
|<SCRIPT type="text/javascript">
var monitorJobsWin=window.open("$promsPath{cgi}/monitorJobs.cgi?filterType=Quantification [DESIGN:PROT_RULER]&filterDateNumber=1&filterDateType=DAY&filterStatus=Queued&filterStatus=Running",'monitorJobsWindow','width=1200,height=500,scrollbars=yes,resizable=yes');
monitorJobsWin.focus();
</SCRIPT>
|;

	sleep 5;
	print qq
|<SCRIPT LANGUAGE="JavaScript">
// top.promsFrame.selectedAction='summary'; //
parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=design:$designID&ACT=experiment&EXPERIMENT=EXPERIMENT:$experimentID&ISNAVFRAME=0&VIEW=quanti";
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

sub ajaxDisplayQuantifs {
	
	my $dbh=&promsConfig::dbConnect;
	my ($lightColor,$darkColor)=&promsConfig::getRowColors;

	my %availableStates;
	my $statesNb;
	if ($intParamID) {
		$statesNb = getAvailableStates(\%availableStates, $intParamID);
	}

	print header(-charset=>'utf-8'); # because of ° in source selection
	
	if (%availableStates) {
		print qq
|		<TABLE border="0" cellspacing="0" bgcolor="$lightColor">
		<TR>
			<!--Parent quantif used for proteomic ruler-->
			<TH align=center scope="col" style="width:300px">&nbsp;Parent quantification(s)&nbsp;</TH>
			
			<!--Parent states used for proteomic ruler-->
			<TH align=center scope="col" style="width:150px">&nbsp;Parent state(s)&nbsp;</TH>

			<!--Custom protein quantity for each sample-->
			<TH align=center scope="col" id="customQtyLabel" style="display:none;width:300px">&nbsp;Custom proteins quantities (pg)
			<BR>
			<INPUT type="button" name="protQtyButton" id="protQtyButton" value="Same for all" onclick="reportValues('protQty');">
			</TH>
			<TH align=center scope="col" id="customQtyHidden" style="width:300px">&nbsp;<BR>&nbsp;</TH>

			<!--Cell number for each sample-->
			<TH align=center scope="col" id="cellNbLabel" style="display:none;width:300px">&nbsp;Cell number&nbsp;<BR>&nbsp;for each sample&nbsp;
			<BR>
			<INPUT type="button" name="cellNbButton" id="cellNbButton" value="Same for all" onclick="reportValues('cellNb');">
			</TH>
			<TH align=center scope="col" id="cellNbHidden" style="width:300px">&nbsp;<BR>&nbsp;</TH>

			<!--Groups-->
			<TH align=center scope="col" id="grpLabel" style="display:none;width:250px">&nbsp;Groups&nbsp;</TH>
			<TH align=center scope="col" id="grpHidden" style="width:250px">&nbsp;</TH>
		</TR>
|;
		my @lineColors=($lightColor, $darkColor);
		my $currentColor;
		my $lineCount = 1;
		my $firstState;

		foreach my $parentQuantifID (sort {$a <=> $b} keys %availableStates) {
			$firstState = 1;
			foreach my $parentStateID (sort {$a <=> $b} keys %{$availableStates{$parentQuantifID}->{'states'}}) {
				$currentColor=$lineColors[$lineCount % 2];
				print qq
|		<TR>
			<TD bgcolor="$currentColor" style="text-align:left">&nbsp;&nbsp;&nbsp;&nbsp;
|;
				if ($firstState == 1) {
					print qq
|			$availableStates{$parentQuantifID}->{'name'}</LABEL>
			&nbsp;&nbsp;>
|;
				$firstState = 0;
				}
				print qq
|			</TD>
			<TD bgcolor="$currentColor" style="text-align:left">&nbsp;&nbsp;&nbsp;&nbsp;
			<INPUT type="checkbox" name="quantifBoxes" id="QUANTIF_BOX_$lineCount" value="$parentQuantifID\_$parentStateID" data-box-nb=$lineCount onclick="checkUncheckBox(this, 'quantif');">
			<LABEL for="QUANTIF_BOX_$lineCount">
			&nbsp;&nbsp;&nbsp;&nbsp;$availableStates{$parentQuantifID}->{'states'}{$parentStateID}{'name'}&nbsp;&nbsp;
			</LABEL>
			</TD>
			<TD bgcolor="$currentColor" style="text-align:center">
			<INPUT type="number" min=0 step="any" name="customQty" id="CUSTOM_QTY_$lineCount" disabled style="display:none;width:150px">
			</TD>
			<TD bgcolor="$currentColor" style="text-align:center">
			<INPUT type="number" min=0 step="any" name="cellNb" id="CELL_NB_$lineCount" disabled style="display:none;width:150px">
			</TD>
			<TD bgcolor="$currentColor" style="text-align:center">
			<SELECT name="groups" id="GROUP_$lineCount" disabled style="display:none;width:100px">
				<OPTION value="" selected>-= Select =-</OPTION>
|;
				for my $i (1..$statesNb) {
					print qq
|				<OPTION value="group_$i">Group $i</OPTION>
|;
				}
				print qq
|			</SELECT>
			</TD>
		</TR>
|;
				$lineCount++;
			}
		}
		print qq
|		</TABLE>
|;	
	} else {
		print qq
|		<BR>Select <b>Intensity metric</b> to display available quantifications<BR>&nbsp;
|;
	}
}

sub getAvailableStates {
	my $dbh = &promsConfig::dbConnect;
	my ($refAvailableStates, $parentParamID) = @_;
	my $statesNumber = 0;
	my $sthString;
	if ($parentParamID) {
		$sthString = "SELECT EC.NAME, EC.ID_EXPCONDITION, Q.NAME, Q.ID_QUANTIFICATION
		FROM EXPCONDITION_QUANTIF EQ
		INNER JOIN EXPCONDITION EC ON EQ.ID_EXPCONDITION = EC.ID_EXPCONDITION
		INNER JOIN QUANTIFICATION Q ON EQ.ID_QUANTIFICATION = Q.ID_QUANTIFICATION
		LEFT JOIN PROTEIN_QUANTIFICATION PQ ON Q.ID_QUANTIFICATION=PQ.ID_QUANTIFICATION
		WHERE EC.ID_DESIGN=$designID AND PQ.ID_QUANTIF_PARAMETER=$parentParamID
		GROUP BY Q.ID_QUANTIFICATION, EC.ID_EXPCONDITION";
	} else {
		die "Cannot fetch quantifications / states without quantification parameter ID";
	}
	my $sthGetStates = $dbh->prepare($sthString);
	$sthGetStates->execute;
	while (my ($stateName, $stateID, $parentQuantifName, $parentQuantifID) = $sthGetStates->fetchrow_array) {
		$refAvailableStates->{$parentQuantifID}{'name'} = $parentQuantifName;
		$refAvailableStates->{$parentQuantifID}{'states'}{$stateID} = {'name' => $stateName};
		$statesNumber++;
	}
	$sthGetStates->finish;
	$dbh->disconnect;

	return $statesNumber;
}
=begin comment
Structure of hash of available states :
key 1 = parent quantification ID
	key 2a = 'name'
		value 2a = parent quantification NAME
	key 2b = 'states'
		key 3 = state (expcondition) ID
			key 4a = 'name'
				value 4a = state NAME
	(later)	key 4b = 'targetPos'
				value 4b = state TARGET_POS (if relevant)
	
=end comment
=cut

sub getAnnotStatesOrder {
	my ($statesAnnot, $refStates) = @_;

	my $statePos = 0;
	my @states = split(/;/, $statesAnnot);
    foreach my $stateData (@states) {
		$statePos++;
        $stateData=~s/#//g;  # removes id tags
        my ($numBioRep, $quantiObsIDs, $stateID) = split(/,/, $stateData);

		push @{$refStates}, $stateID;
    }
	return $statePos;
}

#####>Revision history<#####
# 1.0.8 [ENHANCEMENT] Add MyProMS LFQ to intensity metrics (VL 11/02/20)
# 1.0.7 [ENHANCEMENT] Add possibility to separate nb of cells also for "total protein amount" (VL 10/01/20)
# 1.0.6 [FEATURE] Add possibility to enter values per sample + cell number instead of uniquely per cell + minor modifs(VL 10/01/20)
# 1.0.5 [CHANGES] Use new job monitoring window opening parameters (VS 18/11/19)
# 1.0.4 [MODIF] Remove (temporarily) MEAN_STATE from accepted quantifs for Proteomic Ruler because it's not suitable as such (VL 23/10/19)
# 1.0.3 [MODIF] Switch from watchQuantification to monitorJobs script (VS 21/10/19)
# 1.0.2 [BUGFIX] Fix targetPos assignation when multiple states and targetPos from parent quantif (VL 07/10/19) 
# 1.0.1 [ENHANCEMENT] Add possibility to use other intensity metrics than MaxQuant (VL 24/09/19)
# 1.0.0 Created (VL 26/07/2019)
