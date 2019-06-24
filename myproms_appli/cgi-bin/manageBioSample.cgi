#!/usr/local/bin/perl -w

################################################################################
# manageBioSample.cgi               1.0.5                                      #
# Authors: P. Poullet, G. Arras, S.Liva (Institut Curie)              	       #
# Contact: myproms@curie.fr                                                    #
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
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use strict;

#print header; warningsToBrowser(1);#DEBUG

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my $maxTreatment = 5;
#################################
####>Fetching Parameters... <####
#################################
my $action=param ('ACT') || 'summary';
my $projectID = param('projectID');
my $biosampleID = param('biosampleID') || 0;

##AJAX : check sample name and sample code
if ($action eq 'ajaxCheckBiosample') {
	&checkBiosample(param('name'),$biosampleID);
	exit;
}


##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

##User Info and access
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
my $projectFullAccess=($projectAccess =~ /bioinfo/)? 1 : 0;

my ($titleStrg, %refBiosamples, %allSpecies, %allProperties, %biosampleInfo, %sampleTreatments, %sampleProperties);
my $numAllProperties = param('numAllProp') || 0;

if ($action=~/add|edit/ && param('submit')) {

	my $sampleName = param('bioSampName');
	my $speciesID = param('species');
	my $isRef = param('isReference') || 0;
	my $refSampID= param('refBioSample') || undef;
	my $desc = param('desc')? param('desc') : "";
	my $duplicate = param('valueDuplic') || 0;

	##Get all Information for treatment
	my %biosampleTreatment;
	foreach my $treatPos (1..$maxTreatment) {
		next unless param('treatment_'.$treatPos);
		my $strgTreatment=(param('step_'.$treatPos))? "stepValue=".param('step_'.$treatPos) : "stepValue=1"; # in case Step is disabled
		if (param('quantity_'.$treatPos)) {
			$strgTreatment .= "::quantity=".param('quantity_'.$treatPos)."::quantUnit=".param('quantUnit_'.$treatPos);
		}
		if (param('duration_'.$treatPos)) {
			$strgTreatment .= "::duration=".param('duration_'.$treatPos)."::durUnit=".param('durUnit_'.$treatPos);
		}
		#if (param('temperature_'.$treatPos)) {
		#	$strgTreatment .= "::temperature=".param('temperature_'.$treatPos);
		#}
		$biosampleTreatment{param('treatment_'.$treatPos)}{$treatPos} = $strgTreatment;
	}

	##Get all information for properties
	foreach my $propPos (1..$numAllProperties) {
		next unless param('propName_'.$propPos);
		$biosampleTreatment{param('propName_'.$propPos)}{$propPos} = param('propValue_'.$propPos);
	}

	my $sthInsertBiosamp = $dbh -> prepare("INSERT INTO BIOSAMPLE (ID_BIOSAMPLE, ID_REFBIOSAMPLE, ID_SPECIES, IS_REFERENCE, NAME, DES,RECORD_DATE, UPDATE_DATE, UPDATE_USER) VALUES (?,?,?,?,?,?,NOW(),NOW(),?)");
	my $sthInsertProjBiosamp = $dbh -> prepare("INSERT INTO PROJECT_BIOSAMPLE(ID_PROJECT, ID_BIOSAMPLE) VALUES (?,?)");
	my $sthInsertBiosampProperty = $dbh -> prepare("INSERT INTO BIOSAMPLE_PROPERTY(ID_BIOSAMPLE, ID_PROPERTY, RANK,  PROPERTY_VALUE) VALUES (?,?,?,?)");

	###>UPDATE biosample, biosample_property
	if ($action eq 'edit') {
		my $sthUpdateBiosamp = $dbh -> prepare("UPDATE BIOSAMPLE SET ID_REFBIOSAMPLE = ?, ID_SPECIES = $speciesID, NAME = ?, DES = ?, IS_REFERENCE = ?, UPDATE_DATE = NOW(), UPDATE_USER = '$userID' WHERE ID_BIOSAMPLE = $biosampleID");
		$sthUpdateBiosamp->execute($refSampID,$sampleName,$desc,$isRef);
		$sthUpdateBiosamp->finish;
		#>Properties
		$dbh -> do("DELETE from BIOSAMPLE_PROPERTY where ID_BIOSAMPLE = $biosampleID");
		foreach my $propID (sort{$biosampleTreatment{$a} <=> $biosampleTreatment{$b} } keys %biosampleTreatment) {
			foreach my $rank (keys %{$biosampleTreatment{$propID}}) {
				$sthInsertBiosampProperty -> execute($biosampleID,$propID,$rank,$biosampleTreatment{$propID}{$rank});
			}
		}
	}

	###>CREATE biosample, biosample_property
	else {
		($biosampleID) = $dbh -> selectrow_array("SELECT MAX(ID_BIOSAMPLE) FROM BIOSAMPLE");
		$biosampleID++;
		$sthInsertBiosamp -> execute($biosampleID,$refSampID,$speciesID,$isRef,$sampleName,$desc,$userID);
		$sthInsertProjBiosamp -> execute($projectID, $biosampleID);
		#>Properties
		foreach my $propID (sort{$biosampleTreatment{$a} <=> $biosampleTreatment{$b} } keys %biosampleTreatment) {
			foreach my $rank (keys %{$biosampleTreatment{$propID}}) {
				$sthInsertBiosampProperty -> execute($biosampleID,$propID,$rank,$biosampleTreatment{$propID}{$rank});
			}
		}
	}

	###>Duplicate biosample & biosample_property if needed
	if ($duplicate) {
		my ($dupBiosampID)=$dbh -> selectrow_array("SELECT MAX(ID_BIOSAMPLE) from BIOSAMPLE");
		foreach my $dup (1..$duplicate) {
			$dupBiosampID++;
			my $dupSampleName = $sampleName." #".$dupBiosampID;
			$sthInsertBiosamp -> execute($dupBiosampID,$refSampID,$speciesID,undef,$dupSampleName,$desc,$userID);
			$sthInsertProjBiosamp -> execute($projectID, $dupBiosampID);
			#>Properties
			foreach my $propID (sort{$biosampleTreatment{$a} <=> $biosampleTreatment{$b} } keys %biosampleTreatment) {
				foreach my $rank (keys %{$biosampleTreatment{$propID}}) {
					$sthInsertBiosampProperty -> execute($dupBiosampID,$propID,$rank,$biosampleTreatment{$propID}{$rank});
				}
			}
		}
	}
	$sthInsertBiosamp -> finish;
	$sthInsertProjBiosamp -> finish;

	$sthInsertBiosampProperty -> finish;

	$dbh -> commit;
	$dbh -> disconnect;

	print header(-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Update All Frames</TITLE>
<SCRIPT LANGUAGE="JavaScript">
top.promsFrame.selectedAction = 'summary';
parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=project&ID=$projectID&branchID=biosample:$biosampleID&VIEW="+parent.bioSampleView;
</SCRIPT>
</HEAD>
</HTML>
|;
    exit;
}

elsif ($action eq 'delete') {

	$dbh -> do("DELETE FROM PROJECT_BIOSAMPLE where ID_BIOSAMPLE = $biosampleID and ID_PROJECT = $projectID");
	my ($otherProj)=$dbh->selectrow_array("SELECT 1 FROM PROJECT_BIOSAMPLE WHERE ID_BIOSAMPLE = $biosampleID LIMIT 0,1");
	unless ($otherProj) { # not used in other projects
		$dbh -> do("DELETE FROM BIOSAMPLE_PROPERTY where ID_BIOSAMPLE = $biosampleID");
		$dbh -> do("DELETE FROM BIOSAMPLE where ID_BIOSAMPLE = $biosampleID");
	}
	$dbh -> commit;
	$dbh -> disconnect;

	print header(-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Update All Frames</TITLE>
<SCRIPT LANGUAGE="JavaScript">
top.promsFrame.selectedAction = 'properties';
parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=project&ID=$projectID&branchID=all:$projectID&VIEW="+parent.bioSampleView;
</SCRIPT>
</HEAD>
</HTML>
|;
	exit;
}

my $disableSteps=1;
if ($action eq 'summary' || $action eq 'edit') {

	##Get All info for specific sample
	my $sthSelBiosampleInfo = $dbh -> prepare("SELECT ID_REFBIOSAMPLE,ID_SPECIES, NAME, DES, IS_REFERENCE, RECORD_DATE, UPDATE_DATE, UPDATE_USER from BIOSAMPLE where ID_BIOSAMPLE = $biosampleID");
	$sthSelBiosampleInfo-> execute;
	while (my ($refBioSampID,$speciesID, $name, $desc, $isRef, $recDate, $upDate, $user) = $sthSelBiosampleInfo -> fetchrow_array) {
		if ($refBioSampID && $action eq 'summary') {
			($biosampleInfo{'REF_BIOSAMPLE'})=$dbh->selectrow_array("SELECT NAME FROM BIOSAMPLE WHERE ID_BIOSAMPLE=$refBioSampID");
		}
		else {$biosampleInfo{'REF_BIOSAMPLE'}='None';}
		my ($commonName, $scientificName) = $dbh -> selectrow_array("SELECT COMMON_NAME, SCIENTIFIC_NAME from SPECIES where ID_SPECIES = $speciesID");
		$biosampleInfo{'REF_BIOSAMPLE_ID'} = $refBioSampID || 0;
		$biosampleInfo{'SPECIES_ID'} = $speciesID;
		$biosampleInfo{'SPECIES'} = "$commonName ($scientificName)";
		$biosampleInfo{'NAME'} = $name;
		$biosampleInfo{'DESC'} =  $desc || "";
		$biosampleInfo{'IS_REF'} = $isRef || 0;
		$biosampleInfo{'RECORD_DATE'} = $recDate;
		($biosampleInfo{'REC_DATE'} = $recDate)=~s/ .+//;
		$biosampleInfo{'UPDATE_DATE'} = $upDate;
		($biosampleInfo{'UP_DATE'} = $upDate)=~s/ .+//;
		$biosampleInfo{'UPDATE_USER'} = $user;
	}

	##select properties and treatment for specific sample
	my $sthSelSampProp = $dbh -> prepare("SELECT BP.ID_PROPERTY, BP.RANK, BP.PROPERTY_VALUE, P.NAME, P.PROPERTY_TYPE FROM BIOSAMPLE_PROPERTY BP
					     INNER JOIN PROPERTY P ON BP.ID_PROPERTY = P.ID_PROPERTY
					     WHERE BP.ID_BIOSAMPLE = $biosampleID ORDER BY BP.RANK ");

	$sthSelSampProp -> execute;
	my %distinctSteps;
	while (my($propID, $rank, $propValue, $propName, $propType) = $sthSelSampProp -> fetchrow_array) {
		if ($propType eq 'T') {
			my @strgValue = split("::", $propValue);
			my $stepvalue = (split("=",$strgValue[0]))[1];
			my $strgQuantity = ($strgValue[1])? (split("=",$strgValue[1]))[1] : "";
			my $strgQuantUnit = ($strgValue[2])? (split("=",$strgValue[2]))[1] : "" ;
			my $strgDuration = ($strgValue[3])? (split("=",$strgValue[3]))[1] : "";
			my $strgDurUnit = ($strgValue[4])? (split("=",$strgValue[4]))[1] : "" ;
			@{$sampleTreatments{$rank}} = ($propID, $propName, $stepvalue, $strgQuantity, $strgQuantUnit, $strgDuration, $strgDurUnit);
			$distinctSteps{$stepvalue}=1;
		}
		else {
			@{$sampleProperties{$rank}} = ($propID, $propName, $propValue);
		}
	}
	$disableSteps=(scalar keys %distinctSteps <= 1)? 1 : 0;
	$titleStrg = ($action eq 'summary')? 'Viewing <FONT color="#DD0000">'.$biosampleInfo{'NAME'}.'</FONT>' : 'Editing <FONT color="#DD0000">'.$biosampleInfo{'NAME'}.'</FONT>';
}
else { # add
	$titleStrg = "Create a new biological sample";
	$biosampleInfo{'NAME'} = $biosampleInfo{'SPECIES'} = $biosampleInfo{'DESC'} = $biosampleInfo{'SAMPLE_CODE'} = "";
	$biosampleInfo{'IS_REF'} = $biosampleInfo{'REF_BIOSAMPLE_ID'} = $biosampleInfo{'SPECIES_ID'} = 0;
}

if ($action ne 'summary') {
	##>All Reference Biosamples
	if ($projectFullAccess) {
		my $sthRefBioSamp=$dbh->prepare("SELECT ID_BIOSAMPLE,NAME FROM BIOSAMPLE WHERE IS_REFERENCE=1 AND ID_BIOSAMPLE != $biosampleID");
		$sthRefBioSamp->execute;
		while (my ($refSampID,$refName)=$sthRefBioSamp->fetchrow_array) {
			$refBiosamples{$refSampID}=$refName;
		}
		$sthRefBioSamp->finish;
	}

	##>All Properties & treatments linked to current project or verified
	#my $sthSelProperties = $dbh -> prepare("SELECT ID_PROPERTY, NAME, PROPERTY_TYPE, POSSIBLE_VALUES FROM PROPERTY");
	my $sthSelProperties = $dbh -> prepare("SELECT P.ID_PROPERTY,NAME,PROPERTY_TYPE,POSSIBLE_VALUES
												FROM PROPERTY P
												LEFT OUTER JOIN PROJECT_PROPERTY PP ON P.ID_PROPERTY=PP.ID_PROPERTY
												WHERE PP.ID_PROJECT=$projectID OR P.IS_VERIFIED=1");
	$sthSelProperties -> execute;
	%{$allProperties{treatment}}=%{$allProperties{properties}}=%{$allProperties{values}}=(); # defined prim dim
	while (my($propID, $propName, $propType, $possValues) = $sthSelProperties -> fetchrow_array) {
		if ($propType eq 'T') {
			$allProperties{treatment}{$propID} = $propName;
		}
		else {
			$allProperties{properties}{$propID} = $propName;
			$allProperties{values}{$propID} = $possValues if $possValues;
		}
	}
	$numAllProperties=scalar keys %{$allProperties{properties}};
	$sthSelProperties -> finish;

	##All Species
	my $selSpecies = $dbh -> prepare("SELECT ID_SPECIES, COMMON_NAME, SCIENTIFIC_NAME FROM SPECIES WHERE IS_REFERENCE = 1 ORDER BY COMMON_NAME");
	$selSpecies -> execute;
	while (my($speciesID, $commonName,$scientifName) = $selSpecies -> fetchrow_array) {
		@{$allSpecies{$speciesID}} = ($commonName,$scientifName);
	}
	$selSpecies -> finish;
}


#######################
####>Starting HTML<####
#######################
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my $bgColor=$lightColor;

print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Biological Sample Management</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
|;
if ($action ne 'summary') {
	print qq
|var XHR = null;
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

function updateTreatmentSteps(chkStatus) {
	for (var pos=1; pos<=$maxTreatment; pos++) {
		document.getElementById('step_'+pos).disabled=!chkStatus;
	}
}
function checkTreatmentSteps(treatmentPos) {
	var treatmentStep=document.getElementById('step_'+treatmentPos).value;
	for (var pos=1; pos<=$maxTreatment; pos++) {
		if (pos < treatmentPos) {
			if (document.getElementById('step_'+pos).value > treatmentStep) {
				alert('ERROR: Steps must be in ascending order.');
				document.getElementById('step_'+treatmentPos).selectedIndex=document.getElementById('step_'+pos).selectedIndex;
				return;
			}
		}
		else if (pos > treatmentPos) {
			document.getElementById('step_'+pos).selectedIndex=document.getElementById('step_'+treatmentPos).selectedIndex;
		}
	}
	checkTreatments(treatmentPos)
}
function checkTreatments(treatmentPos) {
	var treatmentID=document.getElementById('treatment_'+treatmentPos).value;
	if (!treatmentID) return;
	var treatmentStep=document.getElementById('step_'+treatmentPos).value;
	for (var pos=1; pos<=$maxTreatment; pos++) {
		if (pos==treatmentPos) continue;
		if (treatmentID==document.getElementById('treatment_'+pos).value && treatmentStep==document.getElementById('step_'+pos).value) {
			alert('ERROR: This treatment if already used in step #'+treatmentStep);
			document.getElementById('treatment_'+treatmentPos).selectedIndex=0;
			return;
		}
	}
}
function numericOnly(e) {
	// deal with unicode character sets
	var unicode = e.charCode ? e.charCode : e.keyCode;
	//alert(unicode);
	// if the key is backspace, tab, point , del
	if (unicode == 8 \|\| unicode == 9 \|\| unicode == 46 \|\| unicode == 127
		// end, begin, arrows
		\|\| (unicode >=35  && unicode <= 40)
		// numbers
		\|\| (unicode >= 48 && unicode <= 57) )
	{
		// we allow the key press
		return true;
	}
	else
	{
		// otherwise we don't
		return false;
	}
}
var possValues = new Object();
var usedPropertiesPos=new Object();
|;
	foreach my $propID (keys %{$allProperties{values}}) {print "possValues[$propID] = '$allProperties{values}{$propID}';\n"; }
	foreach my $rank (keys %sampleProperties) {print "usedPropertiesPos[$rank]=$sampleProperties{$rank}[0];\n";}
	print qq
|function updateProperties(propertyID,propertyPos) {
	var strgInnerHtml;
	if (propertyID) {
		//Checking Property is already used
		for (var pos in usedPropertiesPos) {
			if (usedPropertiesPos[pos]==propertyID) {
				alert('ERROR: Property selected is already used at position #'+pos);
				var selPropName=document.getElementById('propName_'+propertyPos);
				if (usedPropertiesPos[propertyPos]) { // reverse to previous selection
					for (var i=0;i<selPropName.length;i++) {
						if (selPropName[i].value==usedPropertiesPos[propertyPos]) {
							selPropName.selectedIndex=i;
							break;
						}
					}
				}
				else {selPropName.selectedIndex=0;}
				return;
			}
		}
		//document.getElementById('propValue_'+propertyPos).disabled=false;
		nextPropPos=(propertyPos*1)+1;
		if (nextPropPos <= $numAllProperties) {
			document.getElementById('table_'+nextPropPos).style.display='';
		}
		usedPropertiesPos[propertyPos]=propertyID;

		if (possValues[propertyID]) {
			strgInnerHtml = '<SELECT NAME="propValue_'+propertyPos+'" id="propValue_'+propertyPos+'" style="width:200px"><OPTION value="">-= Select =-</OPTION>';
			var listParam = possValues[propertyID].split(':#:');
			for (var value in listParam) {
				strgInnerHtml += '<OPTION value="'+listParam[value]+'">'+listParam[value]+'</OPTION>';
			}
			strgInnerHtml += '</SELECT>';
			//document.getElementById('divInner_'+propertyPos).innerHTML = strgInnerHtml;
		}
		else {
			strgInnerHtml = '<INPUT type="text" name="propValue_'+propertyPos+'" id="propValue_'+propertyPos+'" style="width:200px" maxlength="100" value="">';
			//document.getElementById('propValue_'+propertyPos).disabled=false;
		}
	}
	else {
		//document.getElementById('divInner_'+propertyPos).innerHTML = '';
		strgInnerHtml='';
		delete usedPropertiesPos[propertyPos];
	}
	document.getElementById('divInner_'+propertyPos).innerHTML = strgInnerHtml;
}

var existSampName;
function ajaxCheckBiosample (itemName){
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	paramStrg="ACT=ajaxCheckBiosample&projectID=$projectID&biosampleID=$biosampleID&name="+itemName;
	XHR.open("POST","$promsPath{cgi}/manageBioSample.cgi",false);
	//Send the proper header information along with the request
	XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
	XHR.setRequestHeader("Content-length", paramStrg.length);
	XHR.setRequestHeader("Connection", "close");
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			existSampName = (XHR.responseText.match('###OK###'))? false : true;
		}
	}
	XHR.send(paramStrg);
}

function checkForm(myForm) {

	if (!myForm.bioSampName.value) {
		alert('Give a name to biological sample');
		return false;
	}
	if (!myForm.species.value){
		alert('Choose a species');
		return false;
	}

	ajaxCheckBiosample(myForm.bioSampName.value); // synchronous AJAX

	if (existSampName) {
		alert('The name used for this sample already exists. Please choose another one.');
		return false;
	}

	//treatments
	var okTreat;
	var okProp;
	for (var i=1; i<=$maxTreatment; i++) {
		if (document.getElementById('treatment_'+i).value) {
			if (document.getElementById('quantity_'+i).value) {
				if (document.getElementById('quantUnit_'+i).value) {
					if (document.getElementById('duration_'+i).value) {
						if (!document.getElementById('durUnit_'+i).value) {
							alert('Duration needs a unit');
							return false;
						}
					}
				}
				else {
					alert('Quantity needs a unit');
					return false;
				}
			}
			else if (document.getElementById('quantUnit_'+i).value \|\| document.getElementById('duration_'+i).value \|\| document.getElementById('durUnit_'+i).value) {
				alert('Quantity needs a value');
				return false;
			}
		}
	}

	//properties
	for (var i=1; i<=$numAllProperties; i++ ) {
		if (document.getElementById('propName_'+i).value) {
			if (!document.getElementById('propValue_'+i).value) {
				alert('Missing value for Property "'+document.getElementById('propName_'+i).value+'"!');
				return false;
			}
		}
	}

	//return false;
	return true;//for on submit form
}
|;
}
print qq
|</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" >
<CENTER>
<FONT class="title">$titleStrg</FONT><BR><BR>
|;

####>Summary<####
if ($action eq 'summary') {
	my $referenceStrg = ($biosampleInfo{'IS_REF'})? 'Is reference' : $biosampleInfo{'REF_BIOSAMPLE'};
	print qq
|<TABLE border=0  bgcolor="$darkColor">
	<TR><TH align="right" nowrap>Name :</TH><TD bgcolor=$lightColor width=500>$biosampleInfo{NAME}</TD></TR>
	<TR><TH align="right" valign="top" nowrap>Description :</TH><TD bgcolor=$lightColor>$biosampleInfo{DESC}</TD></TR>
	<TR><TH align="right" nowrap> Species :</TH><TD bgcolor=$lightColor>$biosampleInfo{SPECIES}</TD></TR>
	<TR><TH align="right" nowrap>&nbsp;Reference BioSample :</TH><TD bgcolor=$lightColor>$referenceStrg</TD></TR>
<TR><TH align="right" valign="top" nowrap>&nbsp;Distinctive treatments :</TH>
|;
	if (scalar keys %sampleTreatments) {
		my $stepHeader=($disableSteps)? '' : '<TH class="rbBorder">Step</TH>';
		print qq
|	<TD bgcolor="$lightColor" style="padding:2px"><TABLE width="100%" cellspacing="0">
		<TR bgcolor="$darkColor">$stepHeader<TH class="rbBorder">Name</TH><TH class="rbBorder">Quantity</TH><TH class="bBorder">Duration</TH></TR>
|;
		foreach my $rank (sort{$a<=>$b} keys %sampleTreatments) {
			print "<TR>\n";
			unless ($disableSteps) {
				print "\t<TD align=\"center\">$sampleTreatments{$rank}[2]</TD>\n" ;
			}
			my ($quantityStrg,$durationStrg)=('','');
			if ($sampleTreatments{$rank}[3]) { # value
				$quantityStrg=$sampleTreatments{$rank}[3];
				$quantityStrg.=' '.$sampleTreatments{$rank}[4] if $sampleTreatments{$rank}[4]; # unit
			}
			if ($sampleTreatments{$rank}[5]) { # value
				$durationStrg=$sampleTreatments{$rank}[5];
				$durationStrg.=' '.$sampleTreatments{$rank}[6] if $sampleTreatments{$rank}[6]; # unit
			}
			print qq
|	<TD align="left">$sampleTreatments{$rank}[1]</TD>
	<TD align="center">$quantityStrg</TD>
	<TD align="center">$durationStrg</TD>
	</TR>
|;
		}
		print "</TABLE></TD></TR>\n";
	}
	else {print "<TD bgcolor=\"$lightColor\">None</TD></TR>\n";}
	print qq
|<TR><TH align="right" valign="top" nowrap>Properties :</TH>
|;
	if (scalar keys %sampleProperties) {
		print qq
|	<TD bgcolor="$lightColor" style="padding:2px"><TABLE width="100%" cellspacing="0">
		<TR bgcolor="$darkColor"><TH class="rbBorder">Name</TH><TH class="bBorder" style="width:200px">Value</TH></TR>|;
		foreach my $rank (sort{$a<=>$b} keys %sampleProperties) {
			print "<TR bgcolor=\"$lightColor\"><TD>$sampleProperties{$rank}[1]</TD><TD>$sampleProperties{$rank}[2]</TD></TR>\n";
		}
		print "</TABLE></TD></TR>\n";
	}
	else {print "<TD bgcolor=\"$lightColor\">None</TD></TR>\n";}
	my $updateStrg=($biosampleInfo{'RECORD_DATE'} ne $biosampleInfo{'UPDATE_DATE'})? ". Updated on $biosampleInfo{'UP_DATE'}" : "";
	print qq
|<TR><TH align="right">History :</TH><TD bgcolor="$lightColor">Recorded on $biosampleInfo{REC_DATE}$updateStrg by $biosampleInfo{UPDATE_USER}.</TD></TR>
</TABLE>
<BR><BR>
|;
	my $sthObs=$dbh->prepare("SELECT O.ID_OBSERVATION,O.ID_ANALYSIS,O.TARGET_POS,GROUP_CONCAT(ID_MODIFICATION ORDER BY ID_MODIFICATION SEPARATOR ',')
								FROM OBSERVATION O
								LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION
								WHERE O.ID_BIOSAMPLE=$biosampleID GROUP BY O.ID_OBSERVATION");
	my @sthGetAnaInfo=( # 2d-gel or sample
		$dbh->prepare("SELECT 'experiment',E.NAME,'gel2d',G.NAME,'spot',SP.NAME,'analysis',A.NAME FROM ANALYSIS A,SAMPLE S,SPOT SP,GEL2D G,EXPERIMENT E WHERE G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT AND S.ID_SAMPLE=A.ID_SAMPLE AND VALID_STATUS > 0 AND G.ID_EXPERIMENT=E.ID_EXPERIMENT AND A.ID_ANALYSIS=?"),
		$dbh->prepare("SELECT 'experiment',E.NAME,'sample',S.NAME,'analysis',A.NAME FROM ANALYSIS A,SAMPLE S,EXPERIMENT E WHERE S.ID_SAMPLE=A.ID_SAMPLE AND ID_SPOT IS NULL AND VALID_STATUS > 0 AND S.ID_EXPERIMENT=E.ID_EXPERIMENT AND A.ID_ANALYSIS=?")
	);
	my $sthLQ=$dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND FOCUS='peptide' AND ID_ANALYSIS=? LIMIT 0,1");
	my $sthMod=$dbh->prepare("SELECT PSI_MS_NAME FROM MODIFICATION WHERE ID_MODIFICATION=?");

	my (%hierarchy,%obsAna,%anaList,%modifName,%anaPosition);
	$sthObs->execute;
	while (my ($obsID,$anaID,$targetPos,@mods)=$sthObs->fetchrow_array) {
		next unless $obsID; # LEFT JOIN return NULL if no match!!!
		$obsAna{$obsID}=$anaID;
		$anaList{$anaID}=1;
		foreach my $sthGAI (@sthGetAnaInfo) {
			$sthGAI->execute($anaID);
			@{$hierarchy{$obsID}}=$sthGAI->fetchrow_array;
		}
		if ($targetPos==0) { # Label-free
			push @{$hierarchy{$obsID}},('FREE','Label-free');
		}
		else { # Label
			$sthLQ->execute($anaID);
			my ($quantifAnnot)=$sthLQ->fetchrow_array;
			if ($quantifAnnot) {
				$quantifAnnot=~s/::SOFTWARE=\w+//; # remove software info for back compatibility
				my ($labelTypeStrg,@labelInfo)=split('::',$quantifAnnot);
				if ($targetPos > 0) { # iTRAQ
					my ($labType)=($quantifAnnot =~/TMT/i)?'TMT':'iTRAQ';
					my ($tgPos,$chanName,$labelStrg)=split(';',$labelInfo[$targetPos-1]);
					#push @{$hierarchy{$obsID}},('iTRAQ',$chanName);
					push @{$hierarchy{$obsID}},($labType,$chanName);
				}
				else { # SILAC
					if ($mods[0] && !$modifName{$mods[0]}) { # labeled channel
						$sthMod->execute($mods[0]);
						($modifName{$mods[0]})=$sthMod->fetchrow_array;
					}
					my ($matchedChanName,$matchedLabelStrg);
					INFO:foreach my $infoStrg (@labelInfo) {
						my ($tgPos,$chanName,$labelStrg)=split(';',$infoStrg);
						foreach my $modStrg (split('@',$labelStrg)) { # multiple mods for same label channel
							my @labelData=split('#',$modStrg);
							if ($mods[0]) {
								if ($labelData[4]) { # modID if massChroQ XIC
									foreach my $modID (@mods) { # Compare mod IDs
										if ($modID==$labelData[4]) {
											$matchedChanName=$chanName;
											$matchedLabelStrg=$labelStrg;
											last INFO;
										}
									}
								}
								if ($labelData[1] eq $modifName{$mods[0]}) { # Compare mod names
									$matchedChanName=$chanName;
									$matchedLabelStrg=$labelStrg;
									last INFO;
								}
							}
							elsif ($labelData[1] =~ /No label/i) { # unlabeled channel
								push @{$hierarchy{$obsID}},('SILAC',"$chanName [No label]");
							}
						}
					}
					if ($matchedLabelStrg) { # Labeled channel
						my $targetName="$matchedChanName ";
						my $first=1;
						foreach my $modStrg (split('@',$matchedLabelStrg)) {
							my @labelData=split('#',$modStrg);
							if ($first) {$first=0;}
							else {$targetName.='+';}
							$targetName.="[$labelData[1]&raquo;$labelData[2]]";
						}
						push @{$hierarchy{$obsID}},('SILAC',$targetName);
					}
				}
			}
		}
	}
	$sthObs->finish;
	foreach my $sthGAI (@sthGetAnaInfo) {$sthGAI->finish;}
	$sthLQ->finish;
	$sthMod->finish;

	if (scalar keys %anaList > 1) {
		my $anaStrg=join(',',keys %anaList);
		my @sthGetAnaPos=( # 2d-gel or sample
			$dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS A,SAMPLE S,SPOT SP,GEL2D G,EXPERIMENT E WHERE G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT AND S.ID_SAMPLE=A.ID_SAMPLE AND VALID_STATUS > 0 AND G.ID_EXPERIMENT=E.ID_EXPERIMENT AND A.ID_ANALYSIS IN ($anaStrg) ORDER BY E.DISPLAY_POS,G.DISPLAY_POS,SP.NAME,A.DISPLAY_POS"),
			$dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS A,SAMPLE S,EXPERIMENT E WHERE S.ID_SAMPLE=A.ID_SAMPLE AND ID_SPOT IS NULL AND VALID_STATUS > 0 AND S.ID_EXPERIMENT=E.ID_EXPERIMENT AND A.ID_ANALYSIS IN ($anaStrg) ORDER BY E.DISPLAY_POS,S.DISPLAY_POS,A.DISPLAY_POS")
		);
		my $anaRank=0;
		foreach my $sthGAP (@sthGetAnaPos) {
			$sthGAP->execute;
			while (my ($anaID)=$sthGAP->fetchrow_array) {
				$anaPosition{$anaID}=++$anaRank;
			}
			$sthGAP->finish;
		}
	}
	elsif (scalar keys %anaList==1) { # only 1 Analysis
		my $anaID = (keys %anaList)[0];
		$anaPosition{$anaID}=1;
	}

	my $maxHierarchyIndex=0;
	foreach my $refH (values %hierarchy) {$maxHierarchyIndex=$#{$refH} if $maxHierarchyIndex < $#{$refH};}
	my $totColSpan=($maxHierarchyIndex+1)/2;
	$totColSpan=1 if $totColSpan < 1; # in case no linked Obs

	my @prevHierarchy=();
	print qq
|<TABLE cellspacing=0>
<TR bgcolor="$darkColor"><TD class="title3 bBorder" colspan="$totColSpan" style="padding:2px">&nbsp;Linked Observations:&nbsp;</TD></TR>
|;
	my $numRow=0;
	foreach my $obsID (sort{$anaPosition{$obsAna{$a}}<=>$anaPosition{$obsAna{$b}} || $hierarchy{$a}[-1] cmp $hierarchy{$b}[-1]} keys %hierarchy) {
		#print "<TR><TD bgcolor=\"$lightColor\">",join(':',@{$hierarchy{$obsID}})"</TD></TR>\n";
		$numRow++;
		print "<TR bgcolor=\"$lightColor\">";
		for (my $i=0;$i<$#{$hierarchy{$obsID}}; $i+=2) {
			print "<TD nowrap>";
			my $hideItem=1;
			for (my $j=0; $j<=$i; $j+=2) {
				if (!$prevHierarchy[$j] || $prevHierarchy[$j+1] ne $hierarchy{$obsID}[$j+1]) {
					$hideItem=0;
					last;
				}
			}
			unless ($hideItem) {
				print "&nbsp;".$hierarchy{$obsID}[$i+1];
				print "&nbsp;>" if $i < $#{$hierarchy{$obsID}}-1;
			}
			print "&nbsp;" if $i==$#{$hierarchy{$obsID}}-1;
			print "</TD>";
			if ($i==$#{$hierarchy{$obsID}}-1) {
				for (my $j=$i+2; $j<$maxHierarchyIndex; $j+=2) {print "<TD></TD>";}
			}
		}
		print "</TR>\n";
		print "<TR bgcolor=\"$darkColor\"><TD colspan=$totColSpan></TD></TR>\n";
		@prevHierarchy=@{$hierarchy{$obsID}};
	}
	if ($numRow==0) { # no obs
		print qq
|<TR bgcolor="$lightColor"><TD>&nbsp;<B>None</B></TD></TR>
<TR bgcolor="$darkColor"><TD></TD></TR>
|;
	}
	print "</TABLE>\n";
}

else { # add edit
	print qq
|<FORM name="bioSampForm" method="post" onsubmit="return checkForm(this);">
<INPUT type="hidden" name="ACT" value="$action">
<INPUT type="hidden" name="projectID" value="$projectID">
<INPUT type="hidden" name="numAllProp" value="$numAllProperties">
|;
print "<INPUT type=\"hidden\" name=\"biosampleID\" value=\"$biosampleID\">\n" if $action eq 'edit';
print qq
|<TABLE border=0  bgcolor="$darkColor">
	<TR><TH align="right" nowrap> Name :</TH><TD bgcolor=$lightColor><INPUT type="text" name="bioSampName" size="25" maxlength="50" value="$biosampleInfo{NAME}"/></TD></TR>
	<TR><TH align="right" valign="top" nowrap>Description :</TH><TD bgcolor=$lightColor><TEXTAREA name="desc" rows="2" cols="50">$biosampleInfo{DESC}</TEXTAREA></TD></TR>
	<TR><TH align="right" nowrap> Species :</TH><TD bgcolor=$lightColor>
	<SELECT NAME="species"><OPTION value="">-= Select =-</OPTION>"|;
	foreach my $speciesID (sort{$allSpecies{$a}[0] cmp $allSpecies{$b}[0]} keys %allSpecies) {
		my $selSpecies = (defined($biosampleInfo{'SPECIES_ID'}) && $speciesID == $biosampleInfo{'SPECIES_ID'})? ' selected' : '';
		print qq|<OPTION value="$speciesID"$selSpecies>$allSpecies{$speciesID}[0]&nbsp;&nbsp;($allSpecies{$speciesID}[1])</OPTION>|;
	}
	print "</SELECT></TD></TR>\n";
	if ($projectFullAccess) {
		my $selRefStrg = ($biosampleInfo{'IS_REF'} == 1)? ' checked' : '';
		my $disabIsRefStrg=($biosampleInfo{'REF_BIOSAMPLE_ID'})? ' disabled' : '';
		print qq
|<TR><TH align="right" nowrap>Reference BioSample:</TH><TH align="left" bgcolor="$lightColor">
		<INPUT type="checkbox" name="isReference" id="isReference" value="1"$selRefStrg$disabIsRefStrg>Is reference&nbsp;&nbsp;&nbsp;or&nbsp;&nbsp;&nbsp;Select:
		<SELECT name="refBioSample" id="refBioSample" onchange="document.getElementById('isReference').disabled = (this.value)? true : false"><OPTION value="">-= None =-</OPTION>
|;
		foreach my $refSampID (sort{&promsMod::sortSmart(lc($refBiosamples{$a}),lc($refBiosamples{$b}))} keys %refBiosamples) {
			print "<OPTION value=\"$refSampID\"";
			print ' selected' if $refSampID==$biosampleInfo{'REF_BIOSAMPLE_ID'};
			print ">$refBiosamples{$refSampID}</OPTION>\n";
		}
		print "</TH></TR>\n";
	}
	elsif ($action eq 'edit') {
		my $referenceStrg = ($biosampleInfo{'IS_REF'})? 'Is reference' : $biosampleInfo{'REF_BIOSAMPLE'};
		print qq
|<TR><TH align="right" nowrap>Reference BioSample:</TH><TD bgcolor="$lightColor">$referenceStrg
<INPUT type="hidden" name="isReference" value="$biosampleInfo{IS_REF}"><INPUT type="hidden" name="refBioSample" value="$biosampleInfo{REF_BIOSAMPLE_ID}"></TD></TR>
|;
	}
	my ($chkStepStrg,$disabStepStrg)=($disableSteps)? ('',' disabled') : (' checked','');
	print qq
|	<TR><TH align="right" nowrap>Replicates :</TH><TH align="left" bgcolor=$lightColor>&nbsp;Create <INPUT type="text" style="width:30px" name="valueDuplic" onkeypress="return numericOnly(event);"> additional replicates</TH></TR>
	<TR><TH align="right" valign="top" nowrap>&nbsp;Distinctive treatments :</TH>
		<TD bgcolor="$lightColor" style="padding:2px"><INPUT type="checkbox" name="steps" value="1" onclick="updateTreatmentSteps(this.checked)"$chkStepStrg/><B>Multi-step treatments</B>
		<TABLE width="100%" cellspacing="0">
			<TR bgcolor="$darkColor"><TH class="rbBorder">Step</TH><TH class="rbBorder">Name</TH><TH class="rbBorder">Quantity</TH><TH class="bBorder">Duration</TH></TR>
|;
	foreach my $treatPos (1..$maxTreatment) {
		print qq|<TR><TD valign=top>#<SELECT name="step_$treatPos" id="step_$treatPos" onchange="checkTreatmentSteps($treatPos);"$disabStepStrg>\n|;
		foreach my $s (1..$maxTreatment) {
			my $selStep = ($sampleTreatments{$treatPos} && $s==$sampleTreatments{$treatPos}[2])? ' selected' : '';
			print "<OPTION value=\"$s\"$selStep>$s</OPTION>\n";
		}
		print qq
|		</SELECT></TD>
		<TD valign="top">
			<SELECT name="treatment_$treatPos" id="treatment_$treatPos" onchange="checkTreatments($treatPos);">\n\t
				<OPTION value="">-= Add treatment =-</OPTION>
|;
		if ($allProperties{treatment}) {
			my $usedTretID=($sampleTreatments{$treatPos})? $sampleTreatments{$treatPos}[0] : 0;
			foreach my $treatID (sort{$allProperties{treatment}{$a} cmp $allProperties{treatment}{$b}} keys %{$allProperties{treatment}}) {
				my $selTreat = ($treatID==$usedTretID)? ' selected' : '';
				print "<OPTION value=\"$treatID\"$selTreat>$allProperties{treatment}{$treatID}</OPTION>\n";
			}
		}
		else {print "<OPTION value=\"0\">No treatment</OPTION>\n";}
		print qq
|</SELECT>
</TD>
<TD valign=top nowrap>&nbsp;|;
		my $quantityValue = $sampleTreatments{$treatPos}[3] || "";
		print qq
|			<INPUT class="right" type="text" name="quantity_$treatPos" id="quantity_$treatPos" size="4" value="$quantityValue" onkeypress="return numericOnly(event);"/>
			<SELECT name="quantUnit_$treatPos" id="quantUnit_$treatPos">
				<OPTION value="">-</OPTION>
|;
		my $quantityUnit = $sampleTreatments{$treatPos}[4] || "";
		foreach my $unit (&getConcentrationUnits) {
			my $selQuantUnit = ($unit eq $quantityUnit)? ' selected' : '';
			print "\t<OPTION value=\"$unit\"$selQuantUnit>$unit</OPTION>\n";
		}
		print qq
|</SELECT>
</TD>
<TD valign=top nowrap>&nbsp;|;
		my $durationValue = $sampleTreatments{$treatPos}[5] || "";
		print qq
|		<INPUT class="right" type="text" name="duration_$treatPos" id="duration_$treatPos" size="4" value="$durationValue" onkeypress="return numericOnly(event);"/>
			<SELECT name="durUnit_$treatPos\" id="durUnit_$treatPos">
				<OPTION value="">-</OPTION>
|;
		my $durationUnit = $sampleTreatments{$treatPos}[6] || "";
		foreach my $unit ('sec','min','h','day') {
			my $selDurUnit = ($unit eq $durationUnit)? ' selected' : '';
			print qq "\t<OPTION value=\"$unit\"$selDurUnit>$unit</OPTION>\n";
		}
		print qq
|</SELECT>
</TD>
</TR>
|;
	}
	print qq
|		</TABLE>
		</TD>
	</TR>
	<TR>
		<TH align="right" valign="top" nowrap>Properties :</TH>
		<TD bgcolor=$lightColor style="padding:2px">
		<TABLE border=0 width="100%" cellspacing="0">
			<TR bgcolor="$darkColor">
				<TH class="rbBorder">Name</TH><TH class="bBorder" style="width:202px">Value</TH>
			</TR>
|;
	my $numUsedProperties = scalar keys %sampleProperties;
	my $usedPropIndex = 0;
	my $okProp = 1;
	foreach my $propPos (1..$numAllProperties) {
		my $usedPropID=($sampleProperties{$propPos})? $sampleProperties{$propPos}[0] : 0;
		my $displayStrg = ($propPos == 1 || $propPos <= $numUsedProperties+1)? '' : ' style="display:none;"';
		print qq
|				<TR id="table_$propPos"$displayStrg>
					<TD>
						<SELECT name="propName_$propPos" id="propName_$propPos" onchange="updateProperties(this.value,$propPos)" style="width:250px;">
							<OPTION value="">-= Add property =-</OPTION>
|;
		foreach my $propID (sort{$allProperties{properties}{$a} cmp $allProperties{properties}{$b}} keys %{$allProperties{properties}}) {
			my $selStrg = ($propID == $usedPropID)? ' selected' : '';
			print "<OPTION value=\"$propID\"$selStrg>$allProperties{properties}{$propID}</OPTION>\n";
		}
		my ($propValue,$disabStrg)=($propPos <= $numUsedProperties)? ($sampleProperties{$propPos}[2],'') : ('',' disabled');
		print qq
|</SELECT>
		</TD>
		<TD><DIV id="divInner_$propPos">|;
		if ($action eq 'edit'){
			if ($okProp > $numUsedProperties) {
				print "</SPAN>";
				next;
			}
			if ($usedPropID && $allProperties{values}{$usedPropID}) {
				print qq|<SELECT name="propValue_$propPos" id="propValue_$propPos" style="width:200px"><OPTION value="">-= Select =-</OPTION>|;
				foreach my $param (split(":#:",$allProperties{values}{$usedPropID})) {
					my $selStrg=($param eq $sampleProperties{$propPos}[2])? ' selected' : '';
					print "<OPTION value=\"$param\"$selStrg>$param</OPTION>\n";
				}
				print "</SELECT>";
			}
			else {
				print "<INPUT type=\"text\" name=\"propValue_$propPos\" id=\"propValue_$propPos\"  style=\"width:200px\" maxlength=\"100\" value=\"$propValue\"$disabStrg/>";
			}
			$okProp++;
		}
		print "</DIV></TD></TR>\n"; # </TABLE>
		$usedPropIndex++;#<INPUT type=\"text\" name=\"propValue_$propPos\" id=\"propValue_$propPos\" size=\"40\" maxlength=\"500\" value=\"$propValue\"$disabStrg/>#<DIV id="divInner_$propPos"></DIV>
	}
	print "<TR><TD colspan=2>No property recorded.</TD></TR>\n" unless $numAllProperties;
	print "</TABLE></TD></TR>\n";
	if ($action eq 'edit') {
		my $updateStrg=($biosampleInfo{'RECORD_DATE'} ne $biosampleInfo{'UPDATE_DATE'})? ". Updated on $biosampleInfo{'UP_DATE'}" : "";
		print qq
|<TR><TH align="right">History :</TH><TD bgcolor="$lightColor">Recorded on $biosampleInfo{REC_DATE}$updateStrg by $biosampleInfo{UPDATE_USER}.</TD></TR>
|;
	}
	print qq
|	<TR><TH colspan=2><INPUT type="submit" name="submit" id="buttonItem" value=" Save "></TD></TR>
</TABLE>
</FORM>
|;
}
print qq
|</CENTER>
</BODY>
</HTML>
|;
$dbh -> disconnect;

sub checkBiosample {
	my $dbh=&promsConfig::dbConnect;

	my ($itemName, $bioSampID) = @_;

	my $sthSelBiosampName = $dbh -> prepare("SELECT COUNT(*) FROM BIOSAMPLE B,PROJECT_BIOSAMPLE P WHERE B.ID_BIOSAMPLE=P.ID_BIOSAMPLE AND B.ID_BIOSAMPLE != $bioSampID AND UPPER(NAME) = ?");

	$sthSelBiosampName -> execute(uc($itemName));
	my ($existItemName) = $sthSelBiosampName -> fetchrow_array;
	$sthSelBiosampName -> finish;
	$dbh -> disconnect;

	print header(-charset=>'utf-8'); warningsToBrowser(1);
	if ($existItemName) {print '###BAD###';} else {print '###OK###';}
}

sub getConcentrationUnits {
    my @concUnits=(
	'none',
	'mM','µM','nM','pM','fM',
	'mg/ml','µg/ml','ng/ml','pg/ml','fg/ml',
	'mg','µg','ng','pg','fg',
	'mg/kg',
	'%',
	'gray'
    );
    return (@concUnits);
}



####>Revision history<####
# 1.0.5 Minor bug fix in label-free observation display (PP 11/04/19)
# 1.0.4 Minor modification for TMT (GA 03/04/17)
# 1.0.3 Restrict sample unique naming to project & minor display improvement for iTRAQ Observations (PP 15/04/15)
# 1.0.2 Minor change to maintain bioSample tree view in itemFrame (PP 27/08/14)
# 1.0.1 Updates and code cleaning (PP 15/07/14)
# 1.0.0 New script to create and manage biological sample (SL 24/06/14)