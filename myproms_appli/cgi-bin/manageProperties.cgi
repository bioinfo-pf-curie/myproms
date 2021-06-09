#!/usr/local/bin/perl -w

################################################################################
# manageProperties.cgi              1.0.5                                      #
# Authors: P. Poullet, G. Arras, S. Liva (Institut Curie)                      #
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

#################################
####>Fetching Parameters... <####
#################################
my $action=(param ('ACT'))? param ('ACT') : 'summary';
my $projectID=param('projectID');
my $propertyID=param('propertyID')? param('propertyID') : '';
my $type=param('type'); # O: properties, T: treatments

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
my $projectFullAccess=($projectAccess =~ /bioinfo/)? 1 : 0;
my $propertiesAccess=($projectAccess =~ /guest/)? 0 : 1;
my $disabButton=($propertiesAccess)? '' : ' disabled';
my $MAX_POSSIBLE_VALUES=20;

print header(-charset=>'utf-8'); warningsToBrowser(1);

if (param('submit')) { # add or edit or usage

	my (%allSamples,$rankTreat);

	my $sthSelectRank=$dbh->prepare("SELECT BP.ID_BIOSAMPLE,BP.ID_PROPERTY,BP.PROPERTY_RANK FROM BIOSAMPLE_PROPERTY BP INNER JOIN PROPERTY P ON BP.ID_PROPERTY = P.ID_PROPERTY
					WHERE P.PROPERTY_TYPE='$type' and BP.ID_BIOSAMPLE=? and BP.PROPERTY_RANK>? ORDER BY BP.PROPERTY_RANK");
	my $strgSQL = ($type eq "T")? " AND PROPERTY_RANK=?" : " ";
	my $sthInsertSampProp=$dbh->prepare("INSERT INTO BIOSAMPLE_PROPERTY(ID_BIOSAMPLE,ID_PROPERTY,PROPERTY_RANK,PROPERTY_VALUE) VALUES (?,?,?,?)");
	my $sthUpdateRank=$dbh->prepare("UPDATE BIOSAMPLE_PROPERTY SET PROPERTY_RANK=? WHERE ID_BIOSAMPLE=? AND ID_PROPERTY=?$strgSQL");
	my $sthDeleteSample=$dbh->prepare("DELETE FROM BIOSAMPLE_PROPERTY WHERE ID_BIOSAMPLE=? and ID_PROPERTY=?$strgSQL");
	my $sthUpdatePropValue=$dbh->prepare("UPDATE BIOSAMPLE_PROPERTY SET PROPERTY_VALUE=? WHERE ID_BIOSAMPLE=? and ID_PROPERTY=?$strgSQL");

	if ($action eq "usage") {
		my $sthSelectSample=$dbh->prepare("SELECT B.ID_BIOSAMPLE,B.NAME FROM BIOSAMPLE B INNER JOIN PROJECT_BIOSAMPLE PB ON B.ID_BIOSAMPLE = PB.ID_BIOSAMPLE WHERE PB.ID_PROJECT=$projectID");
		$sthSelectSample->execute();
		while (my($biosampID, $biosampName)=$sthSelectSample->fetchrow_array) {
			$allSamples{$biosampID} = $biosampName;
		}
		$sthSelectSample->finish;
	}

	if ($action eq "usage" && $type eq "O") { ## FROM PROPERTIES
		foreach my $sampID (sort{$a <=> $b} keys %allSamples) {
			my $currentValue=param('selValue_'.$sampID);
			my $oldValue=param('oldValue_'.$sampID);
			next if (!$currentValue && !$oldValue) || ($currentValue eq $oldValue);
			#print "old:$oldValue # current:$currentValue<br>\n";
			if (!$currentValue && $oldValue) {
				##DELETE
				my ($rank)= $dbh->selectrow_array("SELECT BP.PROPERTY_RANK FROM BIOSAMPLE_PROPERTY BP INNER JOIN PROPERTY P ON BP.ID_PROPERTY = P.ID_PROPERTY WHERE P.PROPERTY_TYPE='$type'
													AND BP.ID_BIOSAMPLE=$sampID AND BP.ID_PROPERTY=$propertyID ORDER BY BP.PROPERTY_RANK");
				#print "Deleting property: $allSamples{$sampID} ($sampID) deleted in rank $rank (oldValue:$oldValue currentValue:$currentValue)<BR>";
				$sthDeleteSample->execute($sampID,$propertyID);
				$dbh->commit;
				$sthSelectRank->execute($sampID, $rank);
				while (my($biosampID, $propID, $rankSQL)=$sthSelectRank->fetchrow_array) {
					#print "$biosampID, $propID, $rankSQL<br>";
					$sthUpdateRank->execute($rank,$biosampID, $propID);
					$dbh->commit;
					#print "Updating Rank: $allSamples{$sampID} ($sampID) in rank $rankSQL whith property $propID is now in rank $rank<br>";
					$rank++;
				}
			}
			elsif ($currentValue && !$oldValue) {
				##INSERT
				my ($maxRank)=$dbh->selectrow_array("SELECT MAX(PROPERTY_RANK) FROM BIOSAMPLE_PROPERTY BP INNER JOIN PROPERTY P ON  BP.ID_PROPERTY = P.ID_PROPERTY
							     		WHERE BP.ID_BIOSAMPLE=$sampID AND P.PROPERTY_TYPE = '$type'");
				$maxRank=0 if (!$maxRank);#skip warning if maxRank is null
				my $newRank=$maxRank+1;
				$sthInsertSampProp->execute($sampID,$propertyID,$newRank,$currentValue);
				$dbh->commit;
				#print "Inserting sample: $allSamples{$sampID} ($sampID) inserted with currentValue:$currentValue  in rank $newRank<br>";
			}
			else {
				##UPDATE
				$sthUpdatePropValue->execute($currentValue,$sampID,$propertyID);
				#print "Updating sample: $allSamples{$sampID} ($sampID) updated (oldValue:$oldValue currentValue:$currentValue propValue:$currentValue)<br>";
				$dbh->commit;
			}
		}

	}
	elsif ($action eq "usage" && $type eq "T") { ## FROM TREATMENT
		my (%oldList, %curList);
		foreach my $chkSampID (split(":", param('oldList'))) {
			$oldList{$chkSampID}=1;
		}
		my @currentList = param('list');
		foreach my $sampID (@currentList) {
			$curList{$sampID}=1;
		}

		foreach my $oldID (keys %oldList) {
			if (exists $curList{$oldID}) {
				delete $oldList{$oldID};
				delete $curList{$oldID};
			}
		}
		foreach my $sampID (sort{$a <=> $b} keys %allSamples) {
			my ($strgQuantity, $strgQuantUnit, $strgDuration, $strgDurUnit, @curValue, @propValue);

			my $oldValue=param('oldValue_'.$sampID);
			$rankTreat = (param('rank_'.$sampID))? param('rank_'.$sampID) : "";
			$strgQuantity = (param('quantity_'.$sampID))?  param('quantity_'.$sampID) : "";
			push @curValue, $strgQuantity;
			push @propValue, "quantity=$strgQuantity" if ($strgQuantity);
			$strgQuantUnit = (param('quantUnit_'.$sampID))? param('quantUnit_'.$sampID) : "";
			push @curValue, $strgQuantUnit;
			push @propValue, "quantUnit=$strgQuantUnit" if ($strgQuantUnit);
			$strgDuration = (param('duration_'.$sampID))? param('duration_'.$sampID) : "";
			push @curValue, $strgDuration;
			push @propValue, "duration=$strgDuration" if ($strgDuration);
			$strgDurUnit = (param('durUnit_'.$sampID))? param('durUnit_'.$sampID) : "";
			push @curValue, $strgDurUnit;
			push @propValue, "durUnit=$strgDurUnit" if ($strgDurUnit);
			my $currentValue=join("", @curValue);
			my $step=param('stepValue');
			unshift @propValue, "stepValue=$step";
			my $propertyValue=join("::", @propValue);

			if ($oldList{$sampID}) {
				##DELETE
				my ($rank)=  $rankTreat;
				#print "Deleting treatment: $allSamples{$sampID} ($sampID) deleted in rank $rank (oldValue:$oldValue currentValue:$currentValue)<BR>";
				$sthDeleteSample->execute($sampID,$propertyID,$rank);
				$dbh->commit;
				$sthSelectRank->execute($sampID, $rank);
				while (my($biosampID, $propID, $rankSQL)=$sthSelectRank->fetchrow_array) {
					$sthUpdateRank->execute($rank,$biosampID, $propID, $rankSQL);
					$dbh->commit;
					#print "Updating Rank: $allSamples{$sampID} ($sampID) in rank $rankSQL whith property $propID is now in rank $rank<br>";
					$rank++;
				}
			}
			elsif ($curList{$sampID}) {
				##INSERT
				my ($maxRank)=$dbh->selectrow_array("SELECT MAX(PROPERTY_RANK) FROM BIOSAMPLE_PROPERTY BP INNER JOIN PROPERTY P ON  BP.ID_PROPERTY = P.ID_PROPERTY
									WHERE BP.ID_BIOSAMPLE=$sampID AND P.PROPERTY_TYPE = '$type'");
				$maxRank=0 if (!$maxRank);#skip warning if maxRank is null
				my $newRank=$maxRank+1;
				$sthInsertSampProp->execute($sampID,$propertyID,$newRank,$propertyValue);
				$dbh->commit;
				#print "Inserting sample: $allSamples{$sampID} ($sampID) inserted (oldValue:$oldValue currentValue:$currentValue propValue:$propertyValue) in rank $newRank<br>";
			}
			else {
				##UPDATE
				next if (!$currentValue && !$oldValue) || ($currentValue eq $oldValue);
				$sthUpdatePropValue->execute($propertyValue,$sampID,$propertyID, $rankTreat);
				#print "Updating sample: $allSamples{$sampID} ($sampID) updated (oldValue:$oldValue currentValue:$currentValue propValue:$propertyValue) in rank $rankTreat<br>";
				$dbh->commit;
			}
		}
	}
	else {
		my $propertyName = param('propName');
		my $description = (param('desc'))? param('desc') : "";
		my $useInAna = (param('useInAna') eq 'yes')? 1 : 0;
		my $isVerif = (param('isVerif'))? param('isVerif') : 0;
		my @possibleValues;
		if (param('possibleValues')) {
			foreach my $i (1..$MAX_POSSIBLE_VALUES) {
				push @possibleValues,param('textValue_'.$i) if param('textValue_'.$i);
			}
		}
		my $strgPossValues = (scalar(@possibleValues) > 0)? join(":#:", @possibleValues) : undef;
		if ($action eq 'add') {
			my $insertProperty = $dbh -> prepare("INSERT INTO PROPERTY(ID_PROPERTY, NAME, DES, PROPERTY_TYPE,USE_IN_ANALYSIS, IS_VERIFIED, POSSIBLE_VALUES) VALUES (?,?,?,?,?,?,?)");
			my ($maxPropID) = $dbh -> selectrow_array("SELECT MAX(ID_PROPERTY) FROM PROPERTY");
			my $propID = ++$maxPropID;
			$insertProperty -> execute($propID, $propertyName, $description, $type, $useInAna, $isVerif, $strgPossValues);
			$insertProperty -> finish;
			unless($isVerif) {
				$dbh->do("INSERT INTO PROJECT_PROPERTY (ID_PROJECT,ID_PROPERTY) VALUES ($projectID,$propID)");
			}
		}
		elsif ($action eq 'edit') {
			my $updateProperty = $dbh -> prepare("UPDATE PROPERTY SET NAME=?, DES=?, USE_IN_ANALYSIS=?, IS_VERIFIED=?, POSSIBLE_VALUES=? where ID_PROPERTY=?");
			$updateProperty -> execute($propertyName,$description,$useInAna,$isVerif,$strgPossValues,$propertyID);
			$updateProperty -> finish;
			if ($isVerif) { # Accessible to all
				$dbh->do("DELETE FROM PROJECT_PROPERTY WHERE ID_PROPERTY=$propertyID");
			}
		}
		$dbh -> commit;
	}
	$sthSelectRank->finish;
	$sthUpdateRank->finish;
	$sthDeleteSample->finish;
	$sthInsertSampProp->finish;
	$dbh -> disconnect;
	print qq
|<HTML>
<BODY onload="window.location='./manageProperties.cgi?ACT=summary&projectID=$projectID&type=$type';">
</BODY>
</HTML>
|;
	exit;
}
elsif ($action eq 'delete') { # deletable only if not used for bioSamples of current project & user rights OK

	my $propID = param('propertyID');
	$dbh -> do("DELETE FROM PROJECT_PROPERTY WHERE ID_PROPERTY = $propID");
	my ($otherProj)=$dbh->selectrow_array("SELECT 1 FROM PROJECT_PROPERTY WHERE ID_PROPERTY = $propID LIMIT 0,1");
	unless ($otherProj) { # not used in other projects
		$dbh -> do("DELETE FROM PROPERTY WHERE ID_PROPERTY = $propID");
	}
	$dbh -> commit;
	$dbh -> disconnect;
	#print header;warningsToBrowser(1);
	print qq
|<HTML>
<BODY onload="window.location='./manageProperties.cgi?ACT=summary&projectID=$projectID&type=$type';">
</BODY>
</HTML>
|;
	exit;
}


my (%propList, %propIsUsed, %propEdit, %samplesUsage);
my $sthSelProperty = $dbh -> prepare("SELECT P.ID_PROPERTY,NAME,DES,USE_IN_ANALYSIS,IS_VERIFIED
					FROM PROPERTY P
					LEFT OUTER JOIN PROJECT_PROPERTY PP ON P.ID_PROPERTY=PP.ID_PROPERTY
					WHERE P.PROPERTY_TYPE='$type' AND (PP.ID_PROJECT=$projectID OR P.IS_VERIFIED=1)"); # linked to current project or verified
my $sthUsedProp=$dbh->prepare("SELECT 1 FROM BIOSAMPLE_PROPERTY PY
				JOIN PROJECT_BIOSAMPLE PT ON PY.ID_BIOSAMPLE=PT.ID_BIOSAMPLE
				WHERE PY.ID_PROPERTY=? AND PT.ID_PROJECT=$projectID LIMIT 0,1"); # not used for bioSamples of current project
my $sthSelSample=$dbh->prepare("SELECT B.ID_BIOSAMPLE,B.NAME FROM BIOSAMPLE B INNER JOIN PROJECT_BIOSAMPLE PB ON B.ID_BIOSAMPLE = PB.ID_BIOSAMPLE WHERE PB.ID_PROJECT=$projectID");# for sample usage
my $sthSampProp=$dbh->prepare("SELECT BP.ID_PROPERTY, BP.ID_BIOSAMPLE, BP.PROPERTY_VALUE, BP.PROPERTY_RANK FROM PROJECT_BIOSAMPLE PB
				  INNER JOIN BIOSAMPLE B ON PB.ID_BIOSAMPLE=B.ID_BIOSAMPLE
				  INNER JOIN BIOSAMPLE_PROPERTY BP ON B.ID_BIOSAMPLE=BP.ID_BIOSAMPLE
				  WHERE PB.ID_PROJECT=? AND BP.ID_PROPERTY=?"); # for javascript object
$sthSelProperty->execute;
while (my($propertyID, $propertyName, $propertyDes, $useInAnalysis, $isVerified) = $sthSelProperty -> fetchrow_array) {
	$propList{$propertyID} = [$propertyName, &promsMod::HTMLcompatible($propertyDes), $useInAnalysis, $isVerified];
	$sthUsedProp->execute($propertyID);
	my ($isUSed)=$sthUsedProp->fetchrow_array;
	$propIsUsed{$propertyID} = 1 if $isUSed;
}
$sthSelProperty->finish;
$sthUsedProp->finish;
my (%sampTreat, %sampProp, %distinctStep, %sampChecked);
if ($action eq 'edit' || $action eq 'usage') {
	my $propertyID=param('propertyID');
	my ($desc,$possValues) = $dbh -> selectrow_array("SELECT DES,POSSIBLE_VALUES FROM PROPERTY WHERE ID_PROPERTY = $propertyID");
	#$propEdit{propID} = $propertyID;
	$propEdit{name} = $propList{$propertyID}[0];
	$propEdit{use} = $propList{$propertyID}[2] || 0;
	$propEdit{verif} = $propList{$propertyID}[3] || 0;
	$propEdit{desc} = $desc || '';
	$propEdit{possValues} = $possValues || '';
	$propEdit{predefValues} = $possValues || ''; #for usage, to build select menu
	$propEdit{multiValues}=($propEdit{possValues})? 1 : 0;

	#>Used values (in case "Possible values" activated)
	if ($type eq 'O' && !$propEdit{possValues}) {
		my $sthVal=$dbh->prepare("SELECT DISTINCT PROPERTY_VALUE FROM BIOSAMPLE_PROPERTY WHERE ID_PROPERTY=$propertyID ORDER BY PROPERTY_VALUE");
		$sthVal->execute;
		while (my ($value)=$sthVal->fetchrow_array) {
			$propEdit{possValues}.=':#:' if $propEdit{possValues};
			$propEdit{possValues}.=$value;
		}
		$sthVal->finish;
	}

	#get all samples
	if ($action eq 'usage') {
		$sthSelSample->execute;
		while (my($sampleID, $sampleName) = $sthSelSample -> fetchrow_array) {
			$samplesUsage{$sampleID} = $sampleName;
		}
		$sthSelSample->finish;
		$sthSampProp->execute($projectID, $propertyID);
		while (my($propID, $sampID, $propValue, $rank)=$sthSampProp->fetchrow_array) {
			if ($type eq "O") {
				$sampProp{$sampID}="$propValue";
				#print "$sampID - $propValue<br>";
			}
			else {
				my @strgValue = split("::", $propValue);
				my $stepvalue = (split("=",$strgValue[0]))[1];
				$distinctStep{$stepvalue}=1;
				$sampTreat{$sampID}{$stepvalue}{$propID}= $rank."#".$propValue;
				$sampChecked{$sampID}=1 if (scalar(@strgValue));
			}
		}

	}
}
else {
	$propEdit{name} = '';
	$propEdit{use} = 0;
	$propEdit{verif} = 0;
	$propEdit{desc} = '';
	$propEdit{multiValues}=0;
}

#$dbh -> disconnect;

my ($strgProp1,$strgProp2) = ($type eq 'T')? ('Treatment','Treatments') : ('Property','Properties');
my $titleStrg = (scalar keys %propList)? "Available $strgProp2" : "No $strgProp2 Recorded";

print qq
|<HTML>
<HEAD>
<TITLE></TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
function addProperty() {
	window.location="./manageProperties.cgi?ACT=add&projectID=$projectID&type=$type";
}
function editProperty(propID) {
	window.location="./manageProperties.cgi?ACT=edit&projectID=$projectID&type=$type&propertyID="+propID;
}
function deleteProperty(propID) {
	window.location="./manageProperties.cgi?ACT=delete&projectID=$projectID&type=$type&propertyID="+propID;
}
function usageProperty(propID) {
	window.location="./manageProperties.cgi?ACT=usage&projectID=$projectID&type=$type&propertyID="+propID;
}
function displayTreatmentByStep(myForm, stepVal) {
	window.location="./manageProperties.cgi?ACT=usage&projectID=$projectID&type=$type&propertyID=$propertyID&step="+stepVal;
}
function cancelForm() {
	document.getElementById('createProperty').style.display = 'none';
}
function displayPossibleValues(predefVal) {
	if (predefVal==1) {
		for (let i = 1; i <= $MAX_POSSIBLE_VALUES; i++) {
			if (document.getElementById('textValue_'+i).value \|\| i <= 2) {
				document.getElementById('textValue_'+i).disabled = false;
				document.getElementById('dispDiv_'+i).style.display='block';
			}
		}
		document.getElementById('newValButton').style.display='';
	}
	else {
		for (let i = 1; i <= $MAX_POSSIBLE_VALUES; i++) {
			document.getElementById('textValue_'+i).disabled = false;
			document.getElementById('dispDiv_'+i).style.display='none';
		}
		document.getElementById('newValButton').style.display='none';
	}
}
function addRemoveInput(val,pos) {
	if (val == '+') {
		var lastVisPos=0;
		for (let i = 1; i <= $MAX_POSSIBLE_VALUES; i++) {
			if (document.getElementById('dispDiv_'+i).style.display != 'none') lastVisPos=i;
		}
		if (lastVisPos==$MAX_POSSIBLE_VALUES) {
			alert("The maximum number of possible values has been reached!");
			return;
		}
		lastVisPos++;
		document.getElementById('dispDiv_'+lastVisPos).style.display='block';
		document.getElementById('textValue_'+lastVisPos).disabled = false;
	}
	else if (val == '-') {
		document.getElementById('dispDiv_'+pos).style.display = 'none';
		document.getElementById('textValue_'+pos).value = '';
		document.getElementById('textValue_'+pos).disabled = true;
	}

}
function checkUncheck(chkVal, sampID ) {
	document.getElementById('quantity_'+sampID).disabled = !chkVal;
	document.getElementById('quantUnit_'+sampID).disabled = !chkVal;
	document.getElementById('duration_'+sampID).disabled = !chkVal;
	document.getElementById('durUnit_'+sampID).disabled = !chkVal;
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
function clearUsage() { // property/treatment usage
	/* checkbox (treatment only) */
	const treatChk=document.propertyForm.list;
	if (treatChk) {
		if (treatChk.length) {
			for (let i=0; i<treatChk.length; i++) {
				treatChk[i].checked=false;
			}
		}
		else {treatChk.checked=false;}
	}
	/* select */
	const sourceSelect=document.getElementsByClassName('sourceSelect');
	if (sourceSelect) {
		for (let i=0; i<sourceSelect.length; i++) {
			sourceSelect[i].selectedIndex=0;
		}
	}
	/* input */
	const sourceInput=document.getElementsByClassName('sourceInput');
	if (sourceInput) {
		for (let i=0; i<sourceInput.length; i++) {
			sourceInput[i].value='';
		}
	}
}
function checkForm(myForm) {
	var okForm = true;
|;
	if ($action ne 'usage') {
		print qq
|	if (!myForm.propName.value) {
		alert('ERROR: No name provided for property!');
		//return false;
		okForm = false;
	}
	if (myForm.possibleValues && myForm.possibleValues.value == 1) { // Not for treatments!
		var numValue=0;
		for (var i = 1; i <= $MAX_POSSIBLE_VALUES; i++) {
			if (document.getElementById('textValue_'+i).value) numValue++;
			if (numValue >= 2) break;
		}
		if (numValue < 2) {
			alert('ERROR: At least 2 predifined values are required!');
			//return false;
			okForm = false;
		}
	}
|;
	}
	elsif ($type eq "T") {
		print qq
|	var sampChkList=document.getElementsByName('list');
	for (let i=0; i<sampChkList.length; i++) {
		if (sampChkList[i].checked) {
			var sampID=sampChkList[i].value;
			var quantVal=document.getElementById('quantity_'+sampID).value;
			var quantUnitVal=document.getElementById('quantUnit_'+sampID).value;
			var durVal=document.getElementById('duration_'+sampID).value;
			var durUnitVal=document.getElementById('durUnit_'+sampID).value;
			if (!quantVal && quantUnitVal) {
				alert('ERROR: quantity unit has no value!');
				document.getElementById('quantity_'+sampID).scrollIntoView(1);
				okForm=false;
				return okForm;
			}
			else if(!durVal && durUnitVal) {
				alert('ERROR: duration unit has no value!');
				document.getElementById('duration_'+sampID).scrollIntoView(1);
				okForm=false;
				return okForm;
			}
		}
	}
|;
	}
	print qq
|	return okForm;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" >
<CENTER>
<FONT class="title">$titleStrg</FONT><BR><BR>
<FORM name="propertyForm" method="post" onsubmit="return checkForm(this);">
<INPUT type="hidden" name="type" value="$type">
<INPUT type="hidden" name="projectID" value="$projectID">
<INPUT type="hidden" name="ACT" value="$action">
|;
print "<INPUT type=\"hidden\" name=\"propertyID\" value=\"$propertyID\">\n" if ($action eq 'edit' || $action eq 'usage');
#$dbh -> disconnect;

if (scalar keys %propList) { #display available list of properties
	my $strgColor=$lightColor;
	print qq
|<TABLE border="0" bgcolor="$darkColor" cellspacing="0">
	<TR><TH class="rbBorder">&nbsp;Name&nbsp;</TH><TH class="rbBorder" width=400>&nbsp;Description&nbsp;</TH><TH class="bBorder">&nbsp;Use for analyses&nbsp;</TH><TH class="bBorder">&nbsp;</TH></TR>|;
	foreach my $propID (sort{$propList{$a}[0] cmp $propList{$b}[0]}  keys %propList ) {
		my $strgDisabDel = ($propIsUsed{$propID} || $propList{$propID}[3] == 1 || !$propertiesAccess)? ' disabled' : '';
		my $imgCheck = ($propList{$propID}[3] == 1)? "<img src=$promsPath{images}/good.gif>" : "";
		my $valueUseInAna = ($propList{$propID}[2] == 1)? "Yes" : "No";
		print qq
|	<TR bgcolor=$strgColor>
		<TD valign="top">&nbsp;$propList{$propID}[0]&nbsp;$imgCheck</TD>
		<TD valign="top">$propList{$propID}[1]</TD>
		<TD align="center" valign="top">$valueUseInAna</TD>
		<TD><INPUT type="button" value="Edit" style="height:23px" onclick="editProperty($propID)"$disabButton>
		    <INPUT type="button" value="Delete" style="height:23px" onclick="deleteProperty($propID)"$strgDisabDel>
		    <INPUT type="button" value="Usage" style="height:23px" onclick="usageProperty($propID)">

		</TD>
	</TR>
|;
		$strgColor = ($strgColor eq $lightColor)? $darkColor : $lightColor;
	}
	print qq
|	<TR>
		<TD class="tBorder" colspan="4" align="center">
			<INPUT type="button" value="New $strgProp1" onclick="addProperty()"$disabButton>&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Import from File" onclick="window.location='./importBioSampleData.cgi?projectID=$projectID&context=property'"$disabButton>
		</TD>
	</TR>
</TABLE>
|;
}
else {
	print qq
|<INPUT type="button" class="title2" value=" New $strgProp1 " onclick="addProperty()"$disabButton>&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Import from File" onclick="window.location='./importBioSampleData.cgi?projectID=$projectID&context=property'"$disabButton>
|;
}
if ($action eq 'summary') {
	$dbh->disconnect;
	print qq
|</FORM>
</CENTER>
</BODY>
</HTML>
|;
	exit;
}

#my $strgDIV = ($action eq 'edit' || $action eq 'add' || $action eq 'usage')? "<DIV id=\"createProperty\">" : "<DIV id=\"createProperty\" style=\"display:none\">";
#print "<BR><BR>$strgDIV\n";
print "<BR><BR><DIV id=\"createProperty\">\n"; # Need to hide form if "Cancel"

my $disabCheckbox = ($projectFullAccess)? '' : ' disabled';
my ($selNo, $selYes) = ($propEdit{use})? ('',' selected') : (' selected','');
my $selVerif = ($propEdit{verif})? ' checked' : '' ;
my ($selAny, $selDef) = ($propEdit{multiValues})? ('',' selected') : (' selected', '');
my $buttonValue = ($action eq 'edit')? 'Update' : 'Save';
my %editValues;
if (defined($propEdit{possValues})) {
	my $count = 0;
	foreach my $values (split(":#:", $propEdit{possValues})) {
		$count++;
		$editValues{$count} = $values;
	}
}
if ($action ne 'usage') {
	print qq
|<TABLE border=0  bgcolor="$darkColor">
	<TR><TH nowrap align="right" valign="top">&nbsp;Name :</TH><TD bgcolor="$lightColor"><INPUT TYPE="text" name="propName" value="$propEdit{name}" size="50"></TD></TR>
	<TR><TH nowrap align="right" valign="top">&nbsp;Description :</TH><TD bgcolor="$lightColor"><TEXTAREA name="desc" rows="2" cols="50">$propEdit{desc}</TEXTAREA></TD></TR>
	<TR><TH nowrap align="right" valign="top">&nbsp;Use for analyses :</TH><TD bgcolor="$lightColor"><SELECT name="useInAna"><OPTION value="no"$selNo>No</OPTION><OPTION value="yes"$selYes>Yes</OPTION></SELECT></TD></TR>
	<TR><TH nowrap align="right" valign="top">&nbsp;Is verified :</TH><TD bgcolor="$lightColor"><INPUT TYPE="checkbox" name="isVerif" value=1 $selVerif $disabCheckbox>
|;
	if ($type eq "O") {
		print qq
|	<TR><TH nowrap align="right" valign="top">&nbsp;Values :</TH><TD bgcolor="$lightColor"><SELECT name="possibleValues" onchange="displayPossibleValues(this.value)"><OPTION value="0"$selAny>Any</OPTION><OPTION value="1"$selDef>Predefined</OPTION></SELECT><BR>
|;
		for (my $i = 1; $i <= $MAX_POSSIBLE_VALUES; $i++) {
			my $value=($editValues{$i})? $editValues{$i} : '';
			my ($vis,$disab)=($propEdit{multiValues} && $editValues{$i})? ('','') : ('none','disabled');
			print qq
|<DIV name="displayDiv" id="dispDiv_$i" style="display:$vis">
<INPUT type="text" name="textValue_$i" id="textValue_$i" value="$value" size="30"><INPUT type="button" value="Delete" onclick="addRemoveInput('-',$i)">
</DIV>
|;
		}
		my $addDisplay=(!$propEdit{multiValues} || $editValues{$MAX_POSSIBLE_VALUES})? 'none' : '';
		print qq
|<INPUT type="button" id="newValButton" value="New value" onclick="addRemoveInput('+')" style="display:$addDisplay">
		</TD></TR>
|;
	}
	print qq
|	<TR><TD colspan="2" align="center"><INPUT type="submit" name="submit" value="$buttonValue">&nbsp;&nbsp;&nbsp;<INPUT type="reset" value="Cancel" onclick="cancelForm()"></TD></TR>
</TABLE>
|;
}
else { # usage
	$sthSampProp->execute($projectID, $propertyID);

	$sthSampProp->finish;
	my $step=param('step') || 1;
	my $strgUsage=($type eq "O")? "Add property '$propEdit{name}' to Biological Samples" : "Add Treatment '$propEdit{name}' to Biological Samples";

	print qq
|<FONT class="title">$strgUsage</FONT><BR><BR>
<DIV style="max-height:300px;overflow-y:scroll;display:inline-block">
<TABLE border="0" cellspacing="0">
<TR bgcolor="$darkColor">
|;
	if ($type eq "T") {
		print "<TH class=\"rbBorder\">&nbsp;Select&nbsp;</TH>";
	}
	print "<TH class=\"rbBorder\">&nbsp;Biosamples&nbsp;</TH>";
	if ($type eq "T") {
		print qq|<TH class="bBorder">&nbsp;Value for step&nbsp;<SELECT name="stepValue"style="font-weight:bold" onchange="javascript:displayTreatmentByStep(document.propertyForm,this.value)">|;
		foreach my $stepVal (sort{$a <=> $b} keys %distinctStep) {
			my $selStep= ($stepVal == $step)? ' selected' : ' ';
			print "<OPTION value=\"$stepVal\"$selStep>$stepVal</OPTION>\n";
		}
		print "</SELECT></TH>\n";
	}
	else {
		print "<TH class=\"bBorder\">&nbsp;Value&nbsp;</TH>\n";
	}
	print "</TR>\n";
	my ($disabledInfo, @checkList);
	my $colspan=($type eq "O")? 2 : 3;
	foreach my $sampID (sort{&promsMod::sortSmart(lc($samplesUsage{$a}),lc($samplesUsage{$b}))} keys %samplesUsage) {
		next if ($type eq "T" && $step > 1 && (!$sampTreat{$sampID} || !$sampTreat{$sampID}{$step-1}));
		print "<TR bgcolor=\"$lightColor\">";
		if ($type eq "T") {
			(my $checkSamp, $disabledInfo)=($sampChecked{$sampID} && $sampTreat{$sampID}{$step})? (" checked", " ") : (" "," disabled");
			push @checkList, $sampID if ($sampChecked{$sampID} && $sampTreat{$sampID}{$step});
			print qq|<TH>&nbsp;<INPUT TYPE="checkbox" name="list" onclick="checkUncheck(this.checked, $sampID)" value="$sampID"$checkSamp></TH>|;
		}
		print qq|<TH nowrap align="right">&nbsp;$samplesUsage{$sampID} :</TH>|;
		if ($propEdit{predefValues}) {
			print qq|<TD><SELECT name="selValue_$sampID" class="sourceSelect"><OPTION value="">-=Select=-</OPTION>\n|;
			foreach my $values (split(":#:|#", $propEdit{possValues})) {
				my $sel = (!$sampProp{$sampID})? '' : ($sampProp{$sampID} eq $values )? ' selected' : '';
				print qq|<OPTION value="$values"$sel>$values</OPTION>\n|;
			}
			print "</SELECT></TD>\n";
			my $strgOldVal = ($sampProp{$sampID})? $sampProp{$sampID} : "";
			print qq|<INPUT type="hidden" value="$strgOldVal" name="oldValue_$sampID"></TD>\n|;
		}
		else {
			if ($type eq "O") {
				my $valInput = ($sampProp{$sampID})? $sampProp{$sampID} : "";
				print qq
|<TD><INPUT TYPE="text" name="selValue_$sampID" value="$valInput" class="sourceInput"></TD><INPUT type="hidden" value="$valInput" name="oldValue_$sampID">
|;
			}
			else {
				print "<TD nowrap>";
				if ($sampTreat{$sampID}) {
					my $countStep=0;
					foreach my $stepVal (keys %{$sampTreat{$sampID}}) {
						if ($stepVal == $step) {
							my ($rank,$propValue)=split("#", $sampTreat{$sampID}{$stepVal}{$propertyID});
							my @strgValue = split("::", $propValue);
							my $strgStepVal = ($strgValue[0])? (split("=",$strgValue[0]))[1] : "";;
							my $strgQuantity = ($strgValue[1])? (split("=",$strgValue[1]))[1] : "";
							my $strgQuantUnit = ($strgValue[2])? (split("=",$strgValue[2]))[1] : "" ;
							my $strgDuration = ($strgValue[3])? (split("=",$strgValue[3]))[1] : "";
							my $strgDurUnit = ($strgValue[4])? (split("=",$strgValue[4]))[1] : "" ;
							my @oldValue = ($strgQuantity, $strgQuantUnit, $strgDuration, $strgDurUnit);
							my $strgOldvalue=join("", @oldValue);
							$rank = ($rank)? $rank : "no";
							print qq
|<INPUT TYPE="hidden" VALUE="$rank" NAME="rank_$sampID">
<INPUT TYPE="hidden" VALUE="$strgOldvalue" NAME="oldValue_$sampID">
<INPUT class="right sourceInput" type="text" name="quantity_$sampID" id="quantity_$sampID" size="4" value="$strgQuantity" onkeypress="return numericOnly(event);"$disabledInfo/>
<SELECT name="quantUnit_$sampID" id="quantUnit_$sampID" class="sourceSelect" $disabledInfo>
	<OPTION value="">-</OPTION>
|;
							foreach my $unit (&getConcentrationUnits) {
								my $selQuantUnit = ($unit eq $strgQuantUnit)? ' selected' : '';
								print "<OPTION value=\"$unit\"$selQuantUnit>$unit</OPTION>\n";
							}
							print qq
|</SELECT>
<INPUT class="right sourceInput" type="text" name="duration_$sampID" id="duration_$sampID" size="4" value="$strgDuration" onkeypress="return numericOnly(event);"$disabledInfo/>
<SELECT name="durUnit_$sampID" id="durUnit_$sampID" class="sourceSelect" $disabledInfo>
	<OPTION value="">-</OPTION>
|;
							foreach my $unit ('sec','min','h','day') {
								my $selDurUnit = ($unit eq $strgDurUnit)? ' selected' : '';
								print qq "\t<OPTION value=\"$unit\"$selDurUnit>$unit</OPTION>\n";
							}
							print "</SELECT>\n";
						}
						elsif (!$sampTreat{$sampID}{$step} &&  ($sampTreat{$sampID}{$step+1} || $sampTreat{$sampID}{$step-1})) {
							next if ($countStep==1);
							print qq
|<INPUT TYPE="hidden" VALUE="" NAME="oldValue_$sampID" ID="oldValue_$sampID">
<INPUT class="right sourceInput" type="text" name="quantity_$sampID" id="quantity_$sampID" size="4" value="" onkeypress="return numericOnly(event);"$disabledInfo/>
<SELECT name="quantUnit_$sampID" id="quantUnit_$sampID" class="sourceSelect" $disabledInfo>
	<OPTION value="">-</OPTION>
|;
							foreach my $unit (&getConcentrationUnits) {
								print "<OPTION value=\"$unit\">$unit</OPTION>\n";
							}
							print qq
|</SELECT>
<INPUT class="right sourceInput" type="text" name="duration_$sampID" id="duration_$sampID" size="4" value="" onkeypress="return numericOnly(event);"$disabledInfo/>
<SELECT name="durUnit_$sampID" id="durUnit_$sampID" class="sourceSelect" $disabledInfo>
	<OPTION value="">-</OPTION>
|;
							foreach my $unit ('sec','min','h','day') {
								print qq "<OPTION value=\"$unit\">$unit</OPTION>\n";
							}
							print "</SELECT>\n";
							$countStep++;
						}
					}
				}
				elsif ($step > 1) {
					#print qq||;
				}
				else {
					print qq
|<INPUT TYPE="hidden" VALUE="" NAME="oldValue_$sampID" ID="oldValue_$sampID">
<INPUT class="right sourceInput" type="text" name="quantity_$sampID" id="quantity_$sampID" size="4" value="" onkeypress="return numericOnly(event);"$disabledInfo/>
<SELECT name="quantUnit_$sampID" id="quantUnit_$sampID" class="sourceSelect" $disabledInfo>
	<OPTION value="">-</OPTION>
|;
					foreach my $unit (&getConcentrationUnits) {
						print "\t<OPTION value=\"$unit\">$unit</OPTION>\n";
					}
					print qq
|</SELECT>
<INPUT class="right sourceInput" type="text" name="duration_$sampID" id="duration_$sampID" size="4" value="" onkeypress="return numericOnly(event);"$disabledInfo/>
<SELECT name="durUnit_$sampID" id="durUnit_$sampID" class="sourceSelect" $disabledInfo>
	<OPTION value="">-</OPTION>
|;
					foreach my $unit ('sec','min','h','day') {
						print qq "<OPTION value=\"$unit\">$unit</OPTION>\n";
					}
					print "</SELECT>\n";
				}
				print "</TD>\n";
			}
		}
		print qq
|</TR>
<TR bgcolor="$darkColor"><TD colspan="$colspan"></TD></TR>
|;
	}
	my $disabledSubmit=($type eq "O" && !$propEdit{possValues} )? " disabled" : "";
	my $strgOldList=join(":", @checkList);
	print qq
|<TR bgcolor="$darkColor"><TH colspan="$colspan">&nbsp;<INPUT type="submit" value=" Save "$disabledSubmit name="submit">&nbsp;&nbsp;<INPUT type="button" value="Clear all" onclick="clearUsage()">&nbsp;&nbsp;<INPUT type="reset" value="Cancel" onclick="cancelForm()"><INPUT type="hidden" name="oldList" value="$strgOldList">&nbsp;
</TH></TR>
</TABLE>
</DIV>
|;
}

$dbh->disconnect;

print qq
|</DIV>
</FORM>
</CENTER>
</BODY>
</HTML>
|;

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
# 1.0.5 [FEATURE] Added "Clear all" option for Property/Treatment usage (PP 29/09/20)
# 1.0.4 [UPDATE] Changed RANK field to PROPERTY_RANK for compatibility with MySQL 8 (PP 17/03/20)
# 1.0.3 Display improvement (PP 06/11/15)
# 1.0.2 add one property/treatment to biosamples list, add usage item (SL 08/09/14)
# 1.0.1 Update for PROJECT_PROPERTY table (PP 21/08/14)
# 1.0.0 New script to handle biological sample properties (SL 25/06/14)
