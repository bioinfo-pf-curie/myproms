#!/usr/local/bin/perl -w

################################################################################
# manageModifications.cgi    1.2.1                                             #
# Authors: P. Poullet, G. Arras & F. Yvon (Institut Curie)                     #
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
use strict;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;

####################
####>Parameters<####
####################
my $action = (param('ACT'))? param('ACT'): 'list';
my $startLetter=(param('LET'))? param('LET') : 'A';

if ($action eq 'unimod') {
	&getModificationFromUnimod(param('unimodID'),param('ID'));
	exit;
}

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

if (param('submit')) { # after edit/add
    &processForm;
    exit;
}

#######################
####>Starting HTML<####
#######################
my ($light,$dark)=&promsConfig::getRowColors;
my $bgColor=$light;
print header(-'content-encoding'=>'no',-charset=>'utf-8');
warningsToBrowser(1);
if ($action eq 'list') {
    my $sthGetModifications =($startLetter eq '*')?
    $dbh->prepare("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DES,COMPOSITION,MONO_MASS,AVGE_MASS,UNIMOD_ACC,SPECIFICITY,DISPLAY_CODE,DISPLAY_COLOR,IS_SUBST,IS_LABEL,VALID_STATUS FROM MODIFICATION WHERE PSI_MS_NAME IS NULL OR PSI_MS_NAME NOT REGEXP '^[[:alnum:]]' OR INTERIM_NAME IS NULL OR INTERIM_NAME NOT REGEXP '^[[:alnum:]]' ORDER BY PSI_MS_NAME ASC, INTERIM_NAME ASC, SYNONYMES ASC")
	: ($startLetter eq '?')?
	$dbh->prepare("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DES,COMPOSITION,MONO_MASS,AVGE_MASS,UNIMOD_ACC,SPECIFICITY,DISPLAY_CODE,DISPLAY_COLOR,IS_SUBST,IS_LABEL,VALID_STATUS FROM MODIFICATION WHERE VALID_STATUS != 1 ORDER BY PSI_MS_NAME ASC, INTERIM_NAME ASC, SYNONYMES ASC")
	: ($startLetter eq '#')?
	$dbh->prepare("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DES,COMPOSITION,MONO_MASS,AVGE_MASS,UNIMOD_ACC,SPECIFICITY,DISPLAY_CODE,DISPLAY_COLOR,IS_SUBST,IS_LABEL,VALID_STATUS FROM MODIFICATION WHERE PSI_MS_NAME REGEXP '^[[:digit:]]' OR INTERIM_NAME REGEXP '^[[:digit:]]' OR SYNONYMES REGEXP '##[[:digit:]]' ORDER BY PSI_MS_NAME ASC, INTERIM_NAME ASC, SYNONYMES ASC")
	:
	$dbh->prepare("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DES,COMPOSITION,MONO_MASS,AVGE_MASS,UNIMOD_ACC,SPECIFICITY,DISPLAY_CODE,DISPLAY_COLOR,IS_SUBST,IS_LABEL,VALID_STATUS FROM MODIFICATION WHERE PSI_MS_NAME LIKE '$startLetter%' OR INTERIM_NAME LIKE '$startLetter%' OR SYNONYMES LIKE '##$startLetter%##' ORDER BY PSI_MS_NAME ASC, INTERIM_NAME ASC, SYNONYMES ASC");
    my $sthUsedAna=$dbh->prepare("SELECT 1 FROM ANALYSIS_MODIFICATION WHERE ID_MODIFICATION=? LIMIT 1");
    my $sthUsedProj=$dbh->prepare("SELECT 1 FROM PROJECT_MODIFICATION WHERE ID_MODIFICATION=? LIMIT 1");
	my @modificationList;
	$sthGetModifications->execute;
    #my $modifRef = $sthGetModifications->fetchall_arrayref();
    while (my ($modID,@modData) = $sthGetModifications->fetchrow_array) {
		my $used=0;
		$sthUsedAna->execute($modID);
		if ($sthUsedAna->fetchrow_array) {$used=1;}
		else {
			$sthUsedProj->execute($modID);
			$used=$sthUsedProj->fetchrow_array;
		}
		push @modificationList,[$modID,@modData,$used];
	}
    $sthGetModifications->finish;
	$sthUsedAna->finish;
	$sthUsedProj->finish;
    $dbh->disconnect;

    print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<TITLE>Modifications List</TITLE>
<SCRIPT language="Javascript">
|;
    &promsMod::popupInfo();
    print qq
|function manageModifications(action,modID,letter) {
	document.modifForm.ACT.value=action;
	document.modifForm.ID.value=modID;
	document.modifForm.LET.value=letter;
	document.modifForm.submit();
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<DIV style="float:top">
<BR>
<TABLE><TR><TH bgcolor="$dark">
<FONT class="title2">&nbsp;Go to:</FONT><SELECT class="title action" onchange="if (this.value) {window.location='$promsPath{cgi}/'+this.value;}">
	<OPTION value="">-= Select =-</OPTION>
	<OPTION value="promsMain.cgi">Main Window</OPTION>
	<OPTION value="selectProject.cgi">Project Selection</OPTION>
</SELECT>
</TH></TR></TABLE>
</DIV>
<CENTER>
<FONT class="title">List of Modifications</FONT>&nbsp;&nbsp;<INPUT type="button" class="title2" value="Add modification" onclick="window.location='$promsPath{cgi}/manageModifications.cgi?ACT=add';">
<BR><BR>
<FORM name="modifForm" method="POST">
<INPUT type="hidden" name="ACT" value="">
<INPUT type="hidden" name="ID" value="">
<INPUT type="hidden" name="LET" value="">
</FORM>
<TABLE bgcolor=$dark>
|;
    print "<TR><TH align=right class=\"title3\">&nbsp;Names starting with:</TH><TH bgcolor=$light>";
    foreach my $letter ('#','A'..'Z','*','?') {
		my $letterClass=($letter eq $startLetter)? 'selected' : 'selectable';
		my $mouseoverPopupStrg='';
		if ($letter !~ /\w/) {
			$mouseoverPopupStrg=' onmouseover="popup(\'<B>';
			$mouseoverPopupStrg.=($letter eq '#')? 'Names start with a digit' : ($letter eq '*')? 'Missing PSI-MS or Interim name' : 'Non-validated modifications';
			$mouseoverPopupStrg.='</B>\')" onmouseout="popout()"';
		}
		print "<A class=\"$letterClass\" href=\"javascript:manageModifications('list',0,'$letter')\"$mouseoverPopupStrg>&nbsp;$letter&nbsp;</A>";
		print "|" if $letter ne '?';
    }
    print "</TH></TR>\n";
    print qq
|</TABLE>
<BR>
|;
	if ($startLetter eq '#') {print "<FONT class=\"title2\">Names starting with a digit</FONT><BR>";}
	elsif ($startLetter eq '*') {print "<FONT class=\"title2\">Missing PSI-MS or Interim name</FONT><BR>";}
	elsif ($startLetter eq '?') {print "<FONT class=\"title2\">Non-validated modifications</FONT><BR>";}
	else {print "<FONT class=\"title2\">Names starting with ",lc($startLetter)," or $startLetter</FONT><BR>";}

    if (scalar @modificationList) {
        print "<TABLE border=0 cellspacing=0>\n";
        foreach my $entry (@modificationList) {
            my ($modificationID,$psiName,$interName,$altNames,$des,$compo,$mono,$avg,$unimodID,$spec,$dispCode,$dispColor,$isSubst,$isLabel,$validStatus,$isUsed)=@{$entry};
            my $disabStrg=($isUsed)? 'disabled' : '';
			###> Instentiate values for incomplete modifications
			$psiName='' unless $psiName;
			$interName='' unless $interName;
			$des='' unless $des;
			$compo='' unless $compo;
			$mono='' unless $mono;
			$avg=(!$avg)? '' : ($avg > 0)? '+'.$avg : $avg;
			$unimodID='' unless $unimodID;
			$validStatus=0 unless $validStatus;
			my ($altNameString,$specString)=('','');
			# Format synonyms of the modification
			if ($altNames) {
				$altNameString=$altNames;
				$altNameString=~ s/^##//;
				$altNameString=~ s/##$//;
				$altNameString=~ s/##/, /g;
			}
			# Format specificity of the modification
			if ($spec) {
				foreach my $aaOneLetter (sort{ $a cmp $b } (split(/,/, $spec))) {
					($aaOneLetter,my $context)=split(/;/,$aaOneLetter);
					next unless $aaOneLetter;
					$context='' unless $context;
					$aaOneLetter=($aaOneLetter eq '*') ? 'Any C-term' : ($aaOneLetter eq '+')? 'Protein C-term' : ($aaOneLetter eq '=')? 'Any N-term' : ($aaOneLetter eq '-')? 'Protein N-term' : ($aaOneLetter eq '.')? 'Any residue' : $aaOneLetter;
					$context=($context eq '*') ? ' (Any C-term)' : ($context eq '+')? ' (Protein C-term)' : ($context eq '=')? ' (Any N-term)' : ($context eq '-')? ' (Protein N-term)' : $context;
					$specString.="$aaOneLetter$context,";
				}
				chop($specString);
				$specString=~ s/,/, /g;
				$specString="<BR><B>&nbsp;&nbsp;&nbsp;Specificity:</B> $specString";
			}
			else{
				$specString="<BR><B>&nbsp;&nbsp;&nbsp;Specificity:</B> No residue";
			}
			my $titleName=$psiName || $interName || $altNameString;
			$altNameString="<BR><B>&nbsp;&nbsp;&nbsp;Alternative Name(s):</B> $altNameString" if $altNameString;
			$titleName.=" (<FONT color=\"$dispColor\">$dispCode</FONT>)" if $dispCode;
			my $imageSrc=($validStatus==1)? "$promsPath{images}/lightGreen1.gif" : "$promsPath{images}/lightRed1.gif";
			print qq
|<TR bgcolor=$bgColor>
<TD>&nbsp&nbsp</TD>
<TD width=600><FONT style="font-size:18px;font-weight:bold;">&nbsp;<IMG src='$imageSrc' width=11 height=11>&nbsp;$titleName</FONT>
<BR><B>&nbsp;&nbsp;&nbsp;PSI-MS name:</B> $psiName
<BR><B>&nbsp;&nbsp;&nbsp;Interim name:</B> $interName
$altNameString
<BR><B>&nbsp;&nbsp;&nbsp;Description:</B> $des
<BR><B>&nbsp;&nbsp;&nbsp;Average mass:</B> $avg
<BR><B>&nbsp;&nbsp;&nbsp;Unimod accession:</B> $unimodID
|;
#<BR><B>&nbsp;&nbsp;&nbsp;Composition:</B> $compo
#<BR><B>&nbsp;&nbsp;&nbsp;Monoisotopic mass:</B> $mono
			print qq
|$specString
</TD>
<TH width=100><INPUT type="button" value="Edit" style="width:75px" onclick="manageModifications('edit',$modificationID,'$startLetter')"><BR>
<INPUT type="button" value="Delete" style="width:75px" onclick="manageModifications('delete',$modificationID,'$startLetter')" $disabStrg>
</TH>
</TR>
|;
            $bgColor = ($bgColor eq $light)?$dark:$light;
        }
        print "</TABLE>\n";
    }
    else {print "<BR><FONT class=\"title2\">No matching modifications.</FONT><BR>\n";}
    print qq
|<BR><BR>
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
elsif ($action eq 'add' || $action eq 'edit') {
    my $modificationID = &promsMod::cleanNumericalParameters(param('ID')) || 0;
    my ($psiName,$interName,$altNames,$des,$compo,$mono,$avg,$unimodID,$spec,$dispCode,$dispColor,$isSubst,$isLabel,$validStatus);
	my $isProjectMod=0;
	if ($action eq 'edit') {
		($psiName,$interName,$altNames,$des,$compo,$mono,$avg,$unimodID,$spec,$dispCode,$dispColor,$isSubst,$isLabel,$validStatus)=$dbh->selectrow_array("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DES,COMPOSITION,MONO_MASS,AVGE_MASS,UNIMOD_ACC,SPECIFICITY,DISPLAY_CODE,DISPLAY_COLOR,IS_SUBST,IS_LABEL,VALID_STATUS FROM MODIFICATION WHERE ID_MODIFICATION=$modificationID");
		($isProjectMod)=$dbh->selectrow_array("SELECT 1 FROM PROJECT_MODIFICATION WHERE ID_MODIFICATION=$modificationID LIMIT 1");
		$isProjectMod=0 unless $isProjectMod;
	}
    ###> Instentiate values for incomplete modifications
    # Format synonyms of the modification
    my $altNameString='';
    if ($altNames) {
		$altNameString=$altNames;
		$altNameString=~ s/^##//;
		$altNameString=~ s/##$//;
		$altNameString=~ s/##/, /g;
    }
	my $modifName=$psiName || $interName || $altNameString;
	$psiName='' unless $psiName;
	$interName='' unless $interName;
    $des='' unless $des;
    $compo='' unless $compo;
    $mono='' unless $mono;
    $avg='' unless $avg;
	$spec='' unless $spec;
    $unimodID='' unless $unimodID;
    $dispCode='' unless $dispCode;
    $dispColor=($dispColor)?$dispColor:'FFFFFF';
    my ($isSubstYes,$isSubstNo)=($isSubst)? ('checked','') : ('','checked');
    my ($isLabelYes,$isLabelNo)=($isLabel)? ('checked','') : ('','checked');
    my ($isValidYes,$isValidNo)=($validStatus || $action eq 'add')? ('checked','') : ('','checked');

	my (%specDef,%usedSpec);
	my @AA=('=','-','*','+','A','C'..'I','K'..'N','P'..'T','V'..'Y');
	foreach my $aa (@AA) {
		next if $aa !~ /\w/;
		$specDef{$aa}=$aa;
    }
	$specDef{'-'}='Protein N-term';
	$specDef{'='}='Any N-term';
	$specDef{'+'}='Protein C-term';
	$specDef{'*'}='Any C-term';
	if ($modificationID) { # action=edit
		###<Specificity in projects (no modifications allowed)
		my $sthspecInAna=$dbh->prepare("SELECT DISTINCT(SPECIFICITY) FROM ANALYSIS_MODIFICATION WHERE ID_MODIFICATION=?");
		$sthspecInAna->execute($modificationID);
		while (my ($specString) = $sthspecInAna->fetchrow_array) {
			if ($specString) {
				foreach my $aa (split(',',$specString)) {
					next if $aa eq '?'; # ambiguity
					$usedSpec{$aa}=1; # set as used by Analyses
				}
			}
		}
		$sthspecInAna->finish;
	}
	
	my @destMods;
	if ($action eq 'edit') {
		#my $sthDestNames=$dbh->prepare("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,MONO_MASS FROM MODIFICATION WHERE ID_MODIFICATION != $modificationID AND (PSI_MS_NAME REGEXP '[[:alnum:]]' OR INTERIM_NAME REGEXP '[[:alnum:]]')");
		
		my $sthDestNames=($mono)? $dbh->prepare("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,MONO_MASS FROM MODIFICATION WHERE ID_MODIFICATION != ? AND MONO_MASS >= ".($mono-10)." AND MONO_MASS <= ".($mono+10))
								: $dbh->prepare("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,MONO_MASS FROM MODIFICATION WHERE ID_MODIFICATION != ? AND (COALESCE(PSI_MS_NAME,'') != '' OR COALESCE(INTERIM_NAME,'') != '') ORDER BY PSI_MS_NAME,INTERIM_NAME");
		$sthDestNames->execute($modificationID);
		while (my ($modID,$psiMsName,$interimName,$altNames,$monoMass) = $sthDestNames->fetchrow_array) {
			if ($psiMsName) {$psiMsName=~s/^\s+//; $psiMsName=~s/\s+$//;}
			unless ($psiMsName) {
				if ($interimName) {$interimName=~s/^\s+//; $interimName=~s/\s+$//;}
				unless ($psiMsName) {
					$altNames=~s/^##//;
					$altNames=(split('##',$altNames))[0];
				}
			}
			my $name=$psiMsName || $interimName || $altNames || '*No name*';
			$monoMass=0 unless $monoMass;
			push @destMods,[$modID,$name,$monoMass];
		}
		$sthDestNames->finish;
	}
	$dbh->disconnect;

 	my $actionTitle=($action eq 'add')? 'Adding New Modification' : "Editing Modification <FONT color='#DD0000'>$modifName</FONT>";
   print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<TITLE>Add/Edit Modification</TITLE>
<SCRIPT src="$promsPath{html}/js/local/jscolor.js"></SCRIPT>
<SCRIPT LANGUAGE="JavaScript">
|;
    #&promsMod::popupInfo();
    print qq
|function importFromUnimod() {
	var unimodID=document.getElementById('unimodID').value;
	if (!unimodID \|\| unimodID.match('\\D')) {
		alert('ERROR: Missing or invalid Unimod accession number');
		return false;
	}
	var waitDiv=document.getElementById('waitDIV');
	waitDiv.style.display='';
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/manageModifications.cgi?ACT=unimod&ID=$modificationID&unimodID="+unimodID,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			waitDiv.style.display='none';
//console.log(XHR.responseText);
			processUnimodData(unimodID,XHR.responseText);
		}
	}
	XHR.send(null);	
}
function processUnimodData(unimodID,resultStrg) {
	if (resultStrg.match('^DUPLICATE')) {
		var [dupStrg,dupVarMod]=resultStrg.split('::');
		alert('ERROR: This modification is already recorded ('+dupVarMod+').');
		return false;
	}
	else if (resultStrg.match('^NO_MATCH')) {
		alert('ERROR: Accession #'+unimodID+' is not found is Unimod file.');
		return false;
	}
	var resObject={};
	var resultLines=resultStrg.split('\\n');
	for (let i=0; i<resultLines.length; i++) {
		var [tag,valueStrg]=resultLines[i].split('::');
		resObject[tag]=valueStrg;
	}
//console.log(resObject);
	var myForm=document.editModForm;
	myForm.psiName.value=resObject.PSI_MS_NAME \|\| '';
	myForm.interName.value=resObject.INTERIM_NAME \|\| '';
	myForm.altName.value=resObject.SYNONYMES \|\| '';
	myForm.des.value=resObject.DES \|\| '';
	myForm.composition.value=resObject.COMPOSITION \|\| '';
	myForm.mono.value=resObject.MONO_MASS \|\| '';
	myForm.avg.value=resObject.AVGE_MASS \|\| '';
	if (resObject.IS_LABEL==1) {myForm.isLabel[0].checked=true;} else {myForm.isLabel[1].checked=true;}
	//myForm.isLabel[Math.abs(resObject.IS_LABEL*1-1)].checked=true; // 0 -> 1, 1 -> 0
	/*Specificity*/
	recordSpecificity(resObject.SPECIFICITY);
}
function recordSpecificity(specifString) {
	var myForm=document.editModForm;
	
	/*Check unmodifiable specificity (used in analyses)*/
	var unmodifSpecif={};
	if (myForm.aaspec1) {
		if (myForm.aaspec1.length) {
			for (let i=0; i<myForm.aaspec1.length; i++) {
				var [res,context]=myForm.aaspec1[i].value.split(';');
				unmodifSpecif[res]=1;
			}
		}
		else { // only 1 entry
			var [res,context]=myForm.aaspec1.value.split(';');
			unmodifSpecif[res]=1;
		}
	}
	/*Clear preselected specificity*/
	for (let i=0; i<myForm.aaspec2.length; i++) {
		myForm.aaspec2[i].selectedIndex=0;
		myForm.aaspec2[i].dataset.prevselidx=0;
		myForm.context[i].selectedIndex=0;
		for (let j=1; j<myForm.aaspec2[i].options.length; j++) {
			myForm.aaspec2[i].options[j].disabled=(unmodifSpecif[myForm.aaspec2[i].options[j].value])? true : false;
		}
		if (i>=1) {
			myForm.aaspec2[i].disabled=true;
			myForm.context[i].disabled=true;
			document.getElementById('specDIV_'+(i+1)).style.display='none';
		}
	}
	/*Inject editable specificity*/
	var specRes=specifString.split(',');
	var i=-1;
	for (let s=0; s<specRes.length; s++) {
		var [res,context]=specRes[s].split(';');
		if (unmodifSpecif[res]) continue;
		//if (res.match('[-=+*]')) context=res; // <- done by updateSpecificityDisplay()
		i++;
		myForm.aaspec2[i].disabled=false;
		var pos=i+1;
		//document.getElementById('specDIV_'+pos).style.display=''; // <- done by updateSpecificityDisplay()
		for (let j=1; j<myForm.aaspec2[i].options.length; j++) { // [0] <=> '-= Select =-'
			if (myForm.aaspec2[i].options[j].value==res) {
				//myForm.aaspec2[i].options[j].disabled=false;
				myForm.aaspec2[i].selectedIndex=j;
				if (context) {
					for (let k=1; k<myForm.context[i].options.length; k++) {
						if (myForm.context[i].options[k].value==context) {
							myForm.context[i].selectedIndex=k;
							break;
						}
					}
				}
				updateSpecificityDisplay(myForm.aaspec2[i],pos);
				break;
			}
		}
	}
}
// AJAX ------>
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
// <--- AJAX
function updateSpecificityDisplay(selObject,pos) {
	//Re-enable previously selected AA in all SELECTs
	var allSelects=document.getElementsByName('aaspec2');
	for (let i=0; i<allSelects.length; i++) {
		allSelects[i].options[selObject.dataset.prevselidx].disabled=false;
	}
	
	var nextSpecDiv=document.getElementById('specDIV_'+(pos+1));
	if (selObject.value) {		
		//Disable selected AA in all other SELECTs
		var allSelects=document.getElementsByName('aaspec2');
		for (let i=0; i<allSelects.length; i++) {
			if (i==pos-1) continue;
			allSelects[i].options[selObject.selectedIndex].disabled=true;
		}
		//Update data-prevselidx
		selObject.dataset.prevselidx=selObject.selectedIndex;
		//Handle special case of N/C-terms
		if (selObject.value.match('[-=+*]')) {
			var contextSel=document.getElementById('context_'+pos);
			for (let i=1; i<contextSel.length; i++) {
				if (contextSel.options[i].value==selObject.value) {
					contextSel.selectedIndex=i;
					break;
				}
			}
		}
		//Show next SELECT
		if (nextSpecDiv) {
			nextSpecDiv.style.display='';
			document.getElementById('aaspec2_'+(pos+1)).disabled=false;
			document.getElementById('context_'+(pos+1)).disabled=false;
		}
	}
	else {
		//Update data-prevselidx
		selObject.dataset.prevselidx=0;
		//Hide next SELECT if nothing selected
		if (nextSpecDiv && !document.getElementById('aaspec2_'+(pos+1)).value) {
			nextSpecDiv.style.display='none';
			document.getElementById('aaspec2_'+(pos+1)).disabled=true;
			document.getElementById('context_'+(pos+1)).disabled=true;
		}
	}
}
function resetColor(defaultColor) {
    var colorDisplay = new jscolor.color( document.getElementById( 'dispColor' ));
    colorDisplay.fromString( defaultColor );
}
function checkForm(myForm) {
    if (myForm.mergeMod && myForm.mergeMod.value != '') {
		var selText=myForm.mergeMod.options[myForm.mergeMod.selectedIndex].text;
		return confirm('Are you sure you want to merge this modification with '+selText+' ? ');
    }
	var unimod=myForm.unimodID.value;
	if (!unimod) {
		return confirm('No Unimod accession number recorded. Proceed anyway ?');
	}
	else if (unimod.match('\\D')) {
		alert('ERROR: Invalid Unimod accession number');
		return false;
	}
	if (!myForm.psiName && !myForm.interName && !myForm.altName) {
		alert('ERROR: Please provide a PSI-MS, Interim or Alternative name for Modification');
		return false;
	}
	var okSpecif=false;
	if (myForm.aaspec1) {okSpecif=true;}
	else {
		var allSelects=document.getElementsByName('aaspec2');
		for (let i=0; i<allSelects.length; i++) {
			if (allSelects[i].value && !allSelects[i].disabled) {
				okSpecif=true;
				break;
			}
		}
	}
	if (!okSpecif) {
		alert('ERROR: No site specificity defined');
		return false;
	}
    if ($isProjectMod > 0) {
		if (myForm.code.value == '') {
			alert('Modification already assigned as relevant in a Project! Project display code should not be empty!');
			return false;
		}
		/* Why not white?!!!! (PP)
		else if (myForm.dispColor.value == 'FFFFFF') {
			alert('Modification already assigned as relevant in a Project! Project display color should not be white!');
			return false;
		}
		*/
    }
    return true;
}
window.onload=function() {
	recordSpecificity('$spec');
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR>
<FONT class="title">$actionTitle</FONT><BR>

<FORM name="modifForm" method="POST">
<INPUT type="hidden" name="ACT" value="list">
<INPUT type="hidden" name="LET" value="$startLetter">
</FORM>
<FORM name="editModForm" method=POST onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="idModification" value="$modificationID">
<INPUT type="hidden" name="LET" value="$startLetter">

<TABLE border=0 cellpadding=2 bgcolor=$dark>
<TR><TH align=right class="title3">&nbsp;Unimod accession :</TH><TD nowrap bgcolor=$light>
	<INPUT type="text" name="unimodID" id="unimodID" value="$unimodID" placeholder="*None*" size=10>&nbsp;<INPUT type="button" value="Import from Unimod" onclick="importFromUnimod()"/>&nbsp;[<A href="http://www.unimod.org" target="_blank">Link to Unimod</A>]
	<DIV id="waitDIV" style="display:none"><IMG src="$promsPath{images}/scrollbarGreen.gif"></DIV>
</TD></TR>
<TR><TH align=right nowrap>PSI-MS name :</TH><TD nowrap bgcolor=$light><INPUT type="text" name="psiName" value="$psiName" placeholder="*None*" size=60></TD></TR>
<TR><TH align=right nowrap>Interim name :</TH><TD nowrap bgcolor=$light><INPUT type="text" name="interName" value="$interName" placeholder="*None*" size=60></TD></TR>
<TR><TH align=right valign=top nowrap>&nbsp;Alternative name(s) :</TH><TD nowrap bgcolor=$light><TEXTAREA name="altName" rows="1" cols="60">$altNameString</TEXTAREA></TD></TR>
<TR><TH align=right valign=top>Description :</TH><TD nowrap bgcolor=$light><TEXTAREA name="des" rows="2" cols="60">$des</TEXTAREA></TD></TR>
<TR><TH align=right>Composition :</TH><TD nowrap bgcolor=$light><INPUT type="text" name="composition" value="$compo" size=60></TD></TR>
<TR><TH align=right>Monoisotopic mass :</TH><TD nowrap bgcolor=$light><INPUT type="text" name="mono" value="$mono" size=10></TD></TD></TR>
<TR><TH align=right>Average mass :</TH><TD nowrap bgcolor=$light><INPUT type="text" name="avg" value="$avg" size=10></TD></TR>
<TR><TH align=right valign=top nowrap>Specificity :</TH><TD valign=top nowrap bgcolor=$light>
|;
	my %displayedContext=('-'=>'Protein N-term','='=>'Any N-term','+'=>'Protein C-term','*'=>'Any C-term',''=>'Anywhere');
	my %prevSelected;
	if (scalar keys %usedSpec) {
		print "<B>Not editable (used in Analyses):</B><BR>\n";
		foreach my $aa (@AA) {
			if ($usedSpec{$aa}) { # used in Ana => cannot be edited
				#my $contextStrg=($spec=~/$aa,(.)/)? ";$specDef{$aa}[1]" : '';
				my $context=($spec=~/$aa;(.)/)? $1 : '';
				print "&bull;$specDef{$aa}";
				#print " ($displayedContext{$specDef{$aa}[1]})" if $aa=~/\w/;
				my $contextStrg='';
				if ($context) {
					print " ($displayedContext{$context})";
					$contextStrg=";$context";
				}
				elsif ($aa=~/\w/) {
					print " (Anywhere)";
				}
				print "<INPUT type=\"hidden\" name=\"aaspec1\" value=\"$aa$contextStrg\"><BR>\n";
				$prevSelected{$aa}=1;
			}
		}
		print "<B>Editable specificity:</B><BR>\n";
	}
	
	#my $lastSelPos=0;
	#foreach my $i (1..15) {
	#	my $selAA=''; my $selAAIdx=0;
	#	my $specifStrg='';
	#	my $idx=0;
	#	foreach my $aa (@AA) {
	#		$idx++; # 0 => '-= Select =-'
	#		$specifStrg.="<OPTION value=\"$aa\"";
	#		if ($specDef{$aa}[2]==1 && !$selAA && !$prevSelected{$aa}) {
	#			$specifStrg.=' selected';
	#			$prevSelected{$aa}=1;
	#			$selAA=$aa;
	#			$selAAIdx=$idx;
	#			$lastSelPos=$i;
	#		}
	#		elsif ($prevSelected{$aa}) {
	#			$specifStrg.=' disabled';
	#		}
	#		$specifStrg.=">$specDef{$aa}[0]</OPTION>\n";
	#	}
	#	my $disabSelect=($lastSelPos+1 < $i)? ' disabled' : '';
	#	$specifStrg.="</SELECT>&nbsp;Position:<SELECT name=\"context\" id=\"context_$i\"$disabSelect><OPTION value=\"\">Anywhere</OPTION>\n";
	#	foreach my $context ('=','-','*','+') {
	#		$specifStrg.="<OPTION value=\"$context\"";
	#		$specifStrg.=' selected' if ($selAA && $context eq $specDef{$selAA}[1]);
	#		$specifStrg.=">$displayedContext{$context}</OPTION>\n";
	#	}
	#	$specifStrg.="</SELECT>\n";
	#	print "<DIV id=\"specDIV_$i\"";
	#	print ' style="display:none"' if $lastSelPos+1 < $i;
	#	print ">&bull;Site:<SELECT name=\"aaspec2\" id=\"aaspec2_$i\" onchange=\"updateSpecificityDisplay(this,$i)\" data-prevselidx=\"$selAAIdx\"$disabSelect><OPTION value=\"\">-= Select =-</OPTION>\n$specifStrg</DIV>\n"
	#}
	
	foreach my $i (1..15) { # selection is done by JS at window.onload()
		print qq
|<DIV id="specDIV_$i">&bull;Site:<SELECT name="aaspec2" id="aaspec2_$i" onchange="updateSpecificityDisplay(this,$i)" data-prevselidx="">
<OPTION value="">-= Select =-</OPTION>
|;
		foreach my $aa (@AA) {
			print "<OPTION value=\"$aa\">$specDef{$aa}</OPTION>\n";
		}
		print "</SELECT>&nbsp;Position:<SELECT name=\"context\" id=\"context_$i\"><OPTION value=\"\">Anywhere</OPTION>\n";
		foreach my $context ('=','-','*','+') {
			print "<OPTION value=\"$context\">$displayedContext{$context}</OPTION>\n";
		}
		print "</SELECT>\n</DIV>\n";
	}
	
	print qq
|</TD></TR>
<TR><TH align=right>Is label :</TH><TD nowrap bgcolor=$light>&nbsp;<INPUT type="radio" name="isLabel" value="1" $isLabelYes>Yes<INPUT type="radio" name="isLabel" value="0" $isLabelNo>No</TD></TR>
<TR><TH align=right>Is substitution :</TH><TD nowrap bgcolor=$light>&nbsp;<INPUT type="radio" name="isSubst" value="1" $isSubstYes>Yes<INPUT type="radio" name="isSubst" value="0" $isSubstNo>No</TD></TR>
<TR><TH align=right>Is valid :</TH><TD nowrap bgcolor=$light>&nbsp;<INPUT type="radio" name="isValid" value="1" $isValidYes>Yes<INPUT type="radio" name="isValid" value="0" $isValidNo>No</TD></TR>
<TR><TH align=right>Project display :</TH><TD valign=top nowrap bgcolor=$light>&nbsp;<B>-Set code:</B><INPUT type="text" name="code" value="$dispCode" size=3>&nbsp;<B>-Choose color:</B><INPUT id="dispColor" name="dispColor" class="color" value="$dispColor" autocomplete="on" style="background-image: none; background-color: rgb(41, 255, 191); color: rgb(0, 0, 0);"><INPUT type="button" value="Reset color" onclick="resetColor('$dispColor')"></TD></TR>
|;
	if ($action eq 'edit') {
		print qq
|<TR><TH align=right>Merge with :</TH><TD nowrap bgcolor=$light>
<SELECT name='mergeMod'>
|;
		if (scalar @destMods) {
			print "<OPTION value=\"\">-= Select =-</OPTION>\n";
			foreach my $refDest (sort{if ($mono) {abs($a->[2]-$mono) <=> abs($b->[2]-$mono) || $a->[2] <=> $b->[2]} else {lc($a->[1]) cmp lc($b->[1])}} @destMods) {
				my ($modID,$name,$monoMass)=@{$refDest};
				$monoMass='?' unless $monoMass;
				print "<OPTION value=\"$modID\">$name [Mono. mass: $monoMass]</OPTION>\n";
			}
		}
		else {
			my $infoStrg=($mono)? " [$mono &plusmn; 10 daltons]" : '';
			print "<OPTION value=\"\">** No match$infoStrg **</OPTION>\n";
		}
		print "</SELECT></TD></TR>\n";
	}
    print qq
|<TR><TD colspan=2 align=center>
<INPUT type="submit" name="submit" value=" Save ">
&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="document.modifForm.submit();">
</TD></TR>
</TABLE>
</FORM>
</CENTER>
</BODY>
</HTML>
|;

}
elsif ($action eq 'delete') {
	my $modificationID = &promsMod::cleanNumericalParameters(param('ID'));
	print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<TITLE>Deleting Modification</TITLE>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<BR><BR>
<FONT class="title2">Deleting modification...</FONT>
|;
	$dbh->do("DELETE FROM MODIFICATION WHERE ID_MODIFICATION=$modificationID");
	$dbh->commit;
	$dbh->disconnect;
	print "<FONT class=\"title2\"> Done.</FONT>\n";
	sleep 2;
	print qq
|<SCRIPT type="text/javascript">
window.location="./manageModifications.cgi?ACT=list&LET=$startLetter";
</SCRIPT>
</BODY>
</HTML>
|;
}

sub processForm {
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	my $modificationID=&promsMod::cleanNumericalParameters(param('idModification'));
	print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<TITLE>Processing Modification</TITLE>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
|;

	####<Merging 2 modifications>####
	if (param('mergeMod')) {

		print "<BR><FONT class='title3'>Merging two modifications...";
		my $newModID=param('mergeMod');
		my ($altNameOld,$specOld)=$dbh->selectrow_array("SELECT SYNONYMES,SPECIFICITY FROM MODIFICATION WHERE ID_MODIFICATION=$modificationID");
		my ($psiMsName,$interimName,$altName,$specNew,$validStatus)=$dbh->selectrow_array("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES,SPECIFICITY,VALID_STATUS FROM MODIFICATION WHERE ID_MODIFICATION=$newModID");
		if ($altNameOld) {
			if ($altName) {
				foreach my $syn (split(/##/,$altNameOld)) {
					next unless $syn;
					$altName.="$syn##" unless $altName=~/##$syn##/;
				}
			}
			else {$altName=$altNameOld;}
		}

		###<Merge specificity
		my %spec=();
		my $specAllStrg=($specOld && $specNew)? "$specOld,$specNew" : ($specOld)? $specOld : $specNew;
		if ($specAllStrg) {
			foreach my $aa (split(/,/,$specAllStrg)) {
				next unless $aa;
				$spec{$aa}=1;
			}
		}
		my $newSpec=(scalar keys %spec)? join(',', sort{$a cmp $b} keys %spec) : '';

		###<Add the alt. name to the synonyms of the merge modification and update the new specificity
		my $sthUpMod=$dbh->prepare("UPDATE MODIFICATION SET SYNONYMES=?,SPECIFICITY=? WHERE ID_MODIFICATION=?");
		$sthUpMod->execute($altName,$newSpec,$newModID);
		$sthUpMod->finish;

		print " Done.</FONT><BR>\n";

		###<Update Project / Analysis_Modification / Query_Modification

		###<PROJECT_MODIFICATION, ANALYSIS_MODIFICATION and SPECTRUM_MODIFICATION -> same way
		print "<BR><FONT class='title3'>Updating modification data for Project items...";
		my $count=0;
		foreach my $item ('PROJECT','ANALYSIS','SPECTRUM') {
			my $getItemID=$dbh->prepare("SELECT COUNT(ID_$item),ID_$item,ID_MODIFICATION FROM $item\_MODIFICATION WHERE ID_MODIFICATION=? OR ID_MODIFICATION=? GROUP BY ID_$item");
			my $sthUpItemmod=$dbh->prepare("UPDATE $item\_MODIFICATION SET ID_MODIFICATION=? WHERE ID_MODIFICATION=? AND ID_$item=?");
			$getItemID->execute($modificationID,$newModID);
			while (my ($numMod,$itemID,$modID)=$getItemID->fetchrow_array) {
				if ($numMod==1 && $modID == $modificationID) {
					$sthUpItemmod->execute($newModID,$modificationID,$itemID);
				}
				$count++;
				if ($count==100) {
					print '.';
					$count=0;
				}
			}
			$getItemID->finish;
			$sthUpItemmod->finish;
			$dbh->do("DELETE FROM $item\_MODIFICATION WHERE ID_MODIFICATION=$modificationID");
		}
		print " Done.</FONT><BR>\n";

		###<QUERY_MODIFICATION and PEPTIDE_MODIFICATION -> same way
		print "<BR><FONT class='title3'>Updating modification data for Peptides...";
		$count=0;
		foreach my $item ('QUERY','PEPTIDE') {
			my $extString=($item eq 'QUERY')? ',PEP_RANK' : '';
			my $sthItemInfo=$dbh->prepare("SELECT ID_MODIFICATION,ID_$item,POS_STRING,REF_POS_STRING$extString FROM $item\_MODIFICATION WHERE ID_MODIFICATION=? OR ID_MODIFICATION=?"); # ??? WHY $newModID ??? (PP 08/02/17)
			$sthItemInfo->execute($modificationID,$newModID);
			my %itemMods=();
			while (my ($modID,$itemID,$posStr,$refPosStr,$extInfo)=$sthItemInfo->fetchrow_array) {
				$extInfo='' unless $extInfo;
				$refPosStr='' unless $refPosStr;
				$itemMods{"$itemID:$extInfo"}{$modID}{'POS_STRING'}=$posStr;
				($itemMods{"$itemID:$extInfo"}{$modID}{'REF_POS_STRING'},$itemMods{"$itemID:$extInfo"}{$modID}{'PROB_STRING'})=split(/##/,$refPosStr); # PROB_STRING in case REF_POS_STRING with MaxQuant probability data
			}
			$sthItemInfo->finish;

			foreach my $pos ('POS_STRING','REF_POS_STRING') {
				my $sthUpItemPosMod=($item eq 'QUERY')?
				$dbh->prepare("UPDATE QUERY_MODIFICATION SET $pos=?,ID_MODIFICATION=$newModID WHERE ID_MODIFICATION=? AND ID_QUERY=? AND PEP_RANK=?"):
				$dbh->prepare("UPDATE PEPTIDE_MODIFICATION SET $pos=?,ID_MODIFICATION=$newModID WHERE ID_MODIFICATION=? AND ID_PEPTIDE=?");

				foreach my $itemInfo (keys %itemMods) {
					my ($itemID,$extInfo)=split(/:/,$itemInfo);
					$extInfo='' unless $extInfo;
					if (scalar keys %{$itemMods{$itemInfo}} == 2) {
						my $newPosStr=($itemMods{$itemInfo}{$modificationID}{$pos} && $itemMods{$itemInfo}{$newModID}{$pos})? $itemMods{$itemInfo}{$modificationID}{$pos}.$itemMods{$itemInfo}{$newModID}{$pos} : ($itemMods{$itemInfo}{$newModID}{$pos})? $itemMods{$itemInfo}{$newModID}{$pos} : $itemMods{$itemInfo}{$modificationID}{$pos} ;
						if ($newPosStr) {
							my @positions=split(/\./,$newPosStr);
							$newPosStr=join('.', sort{$a<=>$b} @positions );
						}
						$newPosStr.='##'.$itemMods{$itemInfo}{$modificationID}{'PROB_STRING'} if ($pos eq 'REF_POS_STRING' && $itemMods{$itemInfo}{$modificationID}{'PROB_STRING'});
						$newPosStr=undef unless $newPosStr;
						if ($item eq 'QUERY') {
							$sthUpItemPosMod->execute($newPosStr,$newModID,$itemID,$extInfo);
						}
						else {
							$sthUpItemPosMod->execute($newPosStr,$newModID,$itemID);
						}
					}
					else {
						next if $itemMods{$itemInfo}{$newModID}{$pos} || !defined($itemMods{$itemInfo}{$modificationID}{$pos});
						if ($item eq 'QUERY') {
							$sthUpItemPosMod->execute($itemMods{$itemInfo}{$modificationID}{$pos},$modificationID,$itemID,$extInfo) || die $dbh->errstr;;
						}
						else {
							$sthUpItemPosMod->execute($itemMods{$itemInfo}{$modificationID}{$pos},$modificationID,$itemID) || die $dbh->errstr;;
						}
					}
					$count++;
					if ($count>=1000) {
						print '.';
						$count=0;
					}
				}
				$sthUpItemPosMod->finish;
			}
			$dbh->do("DELETE FROM $item\_MODIFICATION WHERE ID_MODIFICATION=$modificationID");
		}
		print " Done.</FONT><BR>\n";

		$dbh->do("DELETE FROM MODIFICATION WHERE ID_MODIFICATION=$modificationID");

		$startLetter=($psiMsName)? (split(//,$psiMsName))[0] : ($interimName)? (split(//,$interimName))[0] : ($altName && $altName=~/^##(.)/)? $1 : ($validStatus)? '*' : '?';
	}

	####<Edit or Add a modification>####
	else {
		my $unimodID=param('unimodID') || undef;
		my $psiMsName=param('psiName') || undef;
		my $interimName=param('interName') || undef;
		my $altName=param('altName') || undef;
		if ($altName) {
			$altName=~s/,\s*/##/g;
			$altName='##'.$altName.'##';
		}
		my $des=param('des') || undef;
		my $composition=param('composition') || undef;
		my $mono=param('mono') || undef;
		my $avg=param('avg') || undef;
		my ($code,$color)=(param('code') && param('dispColor') && param('dispColor') !~ /^F+$/)? (param('code'),param('dispColor')) : (undef,undef);
		my $isSubst=param('isSubst');
		my $isLabel=param('isLabel');
		my $validStatus=param('isValid');

		###> Check specificity;
		my %specif;
		if (param('aaspec1')) {
			foreach my $spec (param('aaspec1')) {
				$specif{$spec}=1;
			}
		}
		if (param('aaspec2')) {
			my $specIdx=-1;
			foreach my $spec (param('aaspec2')) {
				$specIdx++;
				next unless $spec; # '-= Select =-'
				if ($spec=~/\w/) { # check context for normal residues only
					my $context=(param('context'))[$specIdx];
					$spec.=";$context" if $context; # Anywhere <=> empty context
				}
				$specif{$spec}=1;
			}
		}
		my $newSpec=join(',',sort{$a cmp $b} keys %specif);
		$newSpec=undef unless $newSpec;
		#my $validStatus=($mono && $avg && $newSpec)? 1 : 0;
		
		if ($modificationID) {
			my $sthUpMod=$dbh->prepare("UPDATE MODIFICATION SET UNIMOD_ACC=?,PSI_MS_NAME=?,INTERIM_NAME=?,DES=?,SYNONYMES=?,COMPOSITION=?,MONO_MASS=?,AVGE_MASS=?,SPECIFICITY=?,DISPLAY_CODE=?,DISPLAY_COLOR=?,IS_SUBST=?,IS_LABEL=?,VALID_STATUS=? WHERE ID_MODIFICATION=?");
			$sthUpMod->execute($unimodID,$psiMsName,$interimName,$des,$altName,$composition,$mono,$avg,$newSpec,$code,$color,$isSubst,$isLabel,$validStatus,$modificationID);
			$sthUpMod->finish;
		}
		else { # auto-increment
			my $sthInsMod=$dbh->prepare("INSERT INTO MODIFICATION (UNIMOD_ACC,PSI_MS_NAME,INTERIM_NAME,DES,SYNONYMES,COMPOSITION,MONO_MASS,AVGE_MASS,SPECIFICITY,DISPLAY_CODE,DISPLAY_COLOR,IS_SUBST,IS_LABEL,VALID_STATUS) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
			$sthInsMod->execute($unimodID,$psiMsName,$interimName,$des,$altName,$composition,$mono,$avg,$newSpec,$code,$color,$isSubst,$isLabel,$validStatus);
		}
		$startLetter=($psiMsName)? (split(//,$psiMsName))[0] : ($interimName)? (split(//,$interimName))[0] : ($altName && $altName=~/^##(.)/)? $1 : ($validStatus)? '*' : '?';
	}

	$dbh->commit;
	$dbh->disconnect;

	$startLetter=($startLetter=~/^\d/)? '#' : uc($startLetter);
	print qq
|<SCRIPT type="text/javascript">
window.location="./manageModifications.cgi?ACT=list&LET=$startLetter";
</SCRIPT>
</BODY>
</HTML>
|;
}

#sub createEditSpecificityString {
#    my ($dbh,$modificationID,$specStr,$refAA,$refSpecDef,$refUsUsedSpec)=@_;
#    my %specDef;
#    foreach my $aa (@{$refAA}) {
#		next if $aa !~ /\w/;
#		$specDef{$aa}=$aa;
#    }
#	$specDef{'-'}='Protein N-term';
#	$specDef{'='}='Any N-term';
#	$specDef{'+'}='Protein C-term';
#	$specDef{'*'}='Any C-term';
#	if ($modificationID) { # ! add
#		###<Specificity in projects (no modifications allowed)
#		my $sthspecInAna=$dbh->prepare("SELECT DISTINCT(SPECIFICITY) FROM ANALYSIS_MODIFICATION WHERE ID_MODIFICATION=?");
#		$sthspecInAna->execute($modificationID);
#		while (my ($specString) = $sthspecInAna->fetchrow_array) {
#			if ($specString) {
#				foreach my $aa (split(',',$specString)) {
#					next if $aa eq '?'; # ambiguity
#					$usedSpec{$aa}=1; # set as used by Analyses
#				}
#			}
#		}
#	}
#
#}


sub getModificationFromUnimod {
	print header(-type=>'text/plain',-'content-encoding'=>'no',-charset=>'utf-8');
	my ($unimodID,$modificationID)=&promsMod::cleanNumericalParameters(@_);
	
	my $dbh=&promsConfig::dbConnect;
	my ($matchedModID,$mPsiName,$mInterimName)=$dbh->selectrow_array("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME FROM MODIFICATION WHERE UNIMOD_ACC=$unimodID");
	if ($matchedModID && $matchedModID !=$modificationID) {
		my $name=$mPsiName || $mInterimName;
		print "DUPLICATE::$name";
		$dbh->disconnect;
		exit;
	}
	$dbh->disconnect;
	
	my $xmlFile="$promsPath{data}/unimod_tables.xml";
	my %vmodsInUnimod;
	my $handler = UnimodXMLHandler->new($xmlFile,\%vmodsInUnimod);
require XML::SAX::ParserFactory;
	my $xmlparser = XML::SAX::ParserFactory->parser(Handler => $handler );
	$xmlparser->parse_uri($xmlFile);

	if ($vmodsInUnimod{$unimodID}) {
		my $psiName=$vmodsInUnimod{$unimodID}{'PSI_MS_NAME'} || '';
		my $interimName=$vmodsInUnimod{$unimodID}{'INTERIM_NAME'} || '';
		my $synonymes=$vmodsInUnimod{$unimodID}{'SYNONYMES'} || ''; $synonymes=~s/^##//; $synonymes=~s/##$//; $synonymes=~s/##/, /g;
		my $des=$vmodsInUnimod{$unimodID}{'DES'} || '';
		my $composition=$vmodsInUnimod{$unimodID}{'COMPOSITION'} || '';
		my $monoMass=$vmodsInUnimod{$unimodID}{'MONO_MASS'} || '';
		my $avgeMass=$vmodsInUnimod{$unimodID}{'AVGE_MASS'} || '';
		my $isLabel=$vmodsInUnimod{$unimodID}{'IS_LABEL'} // 0;
		my $specificityString=($vmodsInUnimod{$unimodID}{'SPECIFICITY'})? join(',',keys %{$vmodsInUnimod{$unimodID}{'SPECIFICITY'}}) : '';
		print qq
|PSI_MS_NAME::$psiName
INTERIM_NAME::$interimName
SYNONYMES::$synonymes
DES::$des
COMPOSITION::$composition
MONO_MASS::$monoMass
AVGE_MASS::$avgeMass
IS_LABEL::$isLabel
SPECIFICITY::$specificityString
|;
	}
	else {
		print "NO_MATCH";
	}
	
	exit;
}

####>Revision history<####
# 1.2.1 [Fix] minor JavaScript bugs in form submission (PP 21/03/19)
# 1.2.0 New specificity form to handle site context & 'Add modification' option (PP 11/02/19)
# 1.1.1 Minor modification to avoid error when PSI_MS_NAME is not defined for merge option (GA 07/06/18)
# 1.1.0 Compatible with MaxQuant probality data in PEPTIDE_MODIFICATION.REF_POS_STRING & minro bug fixes (PP 10/02/17)
# 1.0.9 Improved listing options & code update & added X residue (PP 17/08/16)
# 1.0.8 Color on main "Go to" options (PP 25/03/16)
# 1.0.7 Print text for merging procedure so as to avoid timeout (GA 14/04/14)
# 1.0.6 Bug correction to make Reset Color button work<BR>Add valid-status update if enough information is given during editing i.e. mass and specificity (GA 10/09/13)
# 1.0.5 Prevents useless FFFFFF value in DISPLAY_COLOR field after edition (PP 06/09/13)
# 1.0.4 Bug correction to let DISPLAY_CODE null if empty string (GA 17/06/13)
# 1.0.3 Bug fix in createEditSpecificityString when specificity='?' in projects (PP 31/05/13)
# 1.0.2 Modification to avoid to set code/color as empty fields for already assigned relevant modifications in PROJECT_MODIFICATION (GA 21/05/13)
# 1.0.1 Minor change in GPL license (PP 10/05/13)
# 1.0.0 New script to manage modifications (GA 03/05/13)
