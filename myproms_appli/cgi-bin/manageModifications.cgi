#!/usr/local/bin/perl -w

################################################################################
# manageModifications.cgi    1.1.1                                             #
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

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

if (param('submit')) { # after edit
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
    $dbh->prepare("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DES,COMPOSITION,MONO_MASS,AVGE_MASS,UNIMOD_ACC,SPECIFICITY,DISPLAY_CODE,DISPLAY_COLOR,IS_SUBST,IS_LABEL,VALID_STATUS FROM MODIFICATION WHERE PSI_MS_NAME NOT REGEXP '^[[:alnum:]]' OR INTERIM_NAME NOT REGEXP '^[[:alnum:]]' ORDER BY PSI_MS_NAME ASC, INTERIM_NAME ASC , SYNONYMES ASC")
	: ($startLetter eq '?')?
	$dbh->prepare("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DES,COMPOSITION,MONO_MASS,AVGE_MASS,UNIMOD_ACC,SPECIFICITY,DISPLAY_CODE,DISPLAY_COLOR,IS_SUBST,IS_LABEL,VALID_STATUS FROM MODIFICATION WHERE VALID_STATUS != 1 ORDER BY PSI_MS_NAME ASC, INTERIM_NAME ASC , SYNONYMES ASC")
	: ($startLetter eq '#')?
	$dbh->prepare("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DES,COMPOSITION,MONO_MASS,AVGE_MASS,UNIMOD_ACC,SPECIFICITY,DISPLAY_CODE,DISPLAY_COLOR,IS_SUBST,IS_LABEL,VALID_STATUS FROM MODIFICATION WHERE PSI_MS_NAME REGEXP '^[[:digit:]]' OR INTERIM_NAME REGEXP '^[[:digit:]]' OR SYNONYMES REGEXP '##[[:digit:]]' ORDER BY PSI_MS_NAME ASC, INTERIM_NAME ASC , SYNONYMES ASC")
	:
	$dbh->prepare("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DES,COMPOSITION,MONO_MASS,AVGE_MASS,UNIMOD_ACC,SPECIFICITY,DISPLAY_CODE,DISPLAY_COLOR,IS_SUBST,IS_LABEL,VALID_STATUS FROM MODIFICATION WHERE PSI_MS_NAME LIKE '$startLetter%' OR INTERIM_NAME LIKE '$startLetter%' OR SYNONYMES LIKE '##$startLetter%##' ORDER BY PSI_MS_NAME ASC, INTERIM_NAME ASC , SYNONYMES ASC");
    $sthGetModifications->execute;
    my $modifRef = $sthGetModifications->fetchall_arrayref();
    $sthGetModifications->finish;
    $dbh->disconnect;

    print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<TITLE>Modifications List</TITLE>
<SCRIPT language="Javascript">
function manageModifications(action,modID,letter) {
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
<FONT class="title">List of Modifications</FONT><!--&nbsp;&nbsp;<INPUT type="button" class="title2" value="Add modifications" onclick="window.location='$promsPath{cgi}/manageModifications.cgi?ACT=add';" disabled>-->
<BR><BR>
<FORM name="modifForm" method="POST">
<INPUT type="hidden" name="ACT" value="">
<INPUT type="hidden" name="ID" value="">
<INPUT type="hidden" name="LET" value="">
</FORM>
<TABLE bgcolor=$dark>
|;
    print "<TR><TH align=right class=\"title3\">&nbsp;Name starts with:</TH><TH bgcolor=$light>";
    foreach my $letter ('#','A'..'Z','*','?') {
		my $letterClass=($letter eq $startLetter)? 'selected' : 'selectable';
		print "<A class=\"$letterClass\" href=\"javascript:manageModifications('list',0,'$letter')\">&nbsp;$letter&nbsp;</A>";
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

    if (scalar(@{$modifRef})) {
        print "<TABLE border=0 cellspacing=0>\n";
        foreach my $entry (@{$modifRef}) {
            my ($modificationID,$psiName,$interName,$altNames,$des,$compo,$mono,$avg,$unimodID,$spec,$dispCode,$dispColor,$isSubst,$isLabel,$validStatus)=@{$entry};
            #my $disabStrg=()? 'disabled' : '';
			my $disabStrg='disabled';
			###> Intentiate values for incomplete modifications
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
				$altNameString=~ s/##/,/g;
				$altNameString=substr $altNameString,1;
				chop($altNameString);
				$altNameString=~ s/,/,&nbsp;/g;
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
				$specString=~ s/,/,&nbsp;/g;
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
<INPUT type="button" value="Delete" style="width:75px" onclick="deleteModification($modificationID)" $disabStrg>
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
</BODY>
</HTML>
|;
}
elsif ($action eq 'add' || $action eq 'edit') {
    my $modificationID = param('ID');
    my ($psiName,$interName,$altNames,$des,$compo,$mono,$avg,$unimodID,$spec,$dispCode,$dispColor,$isSubst,$isLabel,$validStatus,$sthDetNames);
	my $nbmodInPM=0;
	if ($action eq 'edit') {
		($psiName,$interName,$altNames,$des,$compo,$mono,$avg,$unimodID,$spec,$dispCode,$dispColor,$isSubst,$isLabel,$validStatus)=$dbh->selectrow_array("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DES,COMPOSITION,MONO_MASS,AVGE_MASS,UNIMOD_ACC,SPECIFICITY,DISPLAY_CODE,DISPLAY_COLOR,IS_SUBST,IS_LABEL,VALID_STATUS FROM MODIFICATION WHERE ID_MODIFICATION=$modificationID");
		$sthDetNames=$dbh->prepare("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES FROM MODIFICATION WHERE ID_MODIFICATION!=$modificationID AND PSI_MS_NAME IS NOT NULL ORDER BY PSI_MS_NAME ASC");
		($nbmodInPM)=$dbh->selectrow_array("SELECT COUNT(*) FROM PROJECT_MODIFICATION WHERE ID_MODIFICATION=$modificationID");
	}
    ###> Intentiate values for incomplete modifications
    # Format synonyms of the modification
    my $altNameString='';
    if ($altNames) {
		$altNameString=$altNames;
		$altNameString=~ s/##/,/g;
		$altNameString=~ s/^,//;
		$altNameString=~ s/,$//;
		$altNameString=~ s/,/,&nbsp;/g;
    }
	my $modifName=$psiName || $interName || $altNameString;
	$psiName='*None*' unless $psiName;
	$interName='*None*' unless $interName;
    $des='' unless $des;
    $compo='' unless $compo;
    $mono='' unless $mono;
    $avg='' unless $avg;
    $unimodID='*None*' unless $unimodID;
    $dispCode='' unless $dispCode;
    $dispColor=($dispColor)?$dispColor:'FFFFFF';
    my ($isSubstYes,$isSubstNo)=($isSubst)? ('checked','') : ('','checked') ;
    my $isSubstStr="<INPUT type=\"radio\" name=\"isSubst\" value=\"1\" $isSubstYes>Yes<INPUT type=\"radio\" name=\"isSubst\" value=\"0\" $isSubstNo>No";
    my ($isLabelYes,$isLabelNo)=($isLabel)? ('checked','') : ('','checked') ;
    my $isLabelStr="<INPUT type=\"radio\" name=\"isLabel\" value=\"1\" $isLabelYes>Yes<INPUT type=\"radio\" name=\"isLabel\" value=\"0\" $isLabelNo>No";
    # Format specificity of the modification
    my $specString='';
	if ($spec) {
		foreach my $aaOneLetter (sort{ $a cmp $b } (split(/,+/, $spec))) {
			($aaOneLetter,my $context)=split(/;/,$aaOneLetter);
			$context='' unless $context;
			$aaOneLetter=($aaOneLetter eq '*') ? 'Any C-term' : ($aaOneLetter eq '+')? 'Protein C-term' : ($aaOneLetter eq '=')? 'Any N-term' : ($aaOneLetter eq '-')? 'Protein N-term' : $aaOneLetter;
			$context=($context eq '*') ? ' (Any C-term)' : ($context eq '+')? ' (Protein C-term)' : ($context eq '=')? ' (Any N-term)' : ($context eq '-')? ' (Protein N-term)' : $context;
			$specString.="$aaOneLetter$context,";
		}
		chop($specString);
		$specString=~ s/,/,&nbsp;/g;
    }
	else {
		($specString,$spec)=("No residue",'');
    }
    my ($editSpecString,$hiddenSpecString,%specDef)=&createEditSpecificityString($dbh,$modificationID);
    print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<TITLE>Modification edition</TITLE>
<style type="text/css">
.pushedButton {font-size:11px;font-weight:bold;color:red;}
.releasedButton {font-size:11px;font-weight:normal;color:black;}
</style>
<SCRIPT src="$promsPath{html}/js/local/jscolor.js"></SCRIPT>
<SCRIPT LANGUAGE="JavaScript">
|;
    &promsMod::popupInfo();
    print qq
|var nbmodInPM=$nbmodInPM;
function checkForm(myForm){
    if(myForm.mergeMod.value != ''){
	var selText=myForm.mergeMod.options[myForm.mergeMod.selectedIndex].text;
        return confirm('Are you sure you want to merge this modification with '+selText+' ? ');
    }
    if (nbmodInPM > 0) {
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
function showHideSpecificity(act) {
    if (act=='show') {
		document.getElementById('showSpec').style.display='none';
		document.getElementById('hideSpec').style.display='block';
		document.getElementById('specDIV').style.display='block';
    }
	else{
		document.getElementById('showSpec').style.display='block';
		document.getElementById('hideSpec').style.display='none';
		document.getElementById('specDIV').style.display='none';
    }

}
function changeAASpec(myB) {
    var myBState = document.getElementById('isClicked'+myB.id);
    if (myB.classList.contains('pushedButton')) {
		myB.classList.add('releasedButton');
		myB.classList.remove('pushedButton');
		myBState.value=0;
    }
	else{
		myB.classList.add('pushedButton');
		myB.classList.remove('releasedButton');
		myBState.value=1;
    }
}
var specAA=new Array();
|;
    my $pos=0;
    foreach my $specAA (sort{$a cmp $b} keys %specDef) {
		print "specAA[$pos]=['$specAA','$specDef{$specAA}{DEFINITION}'];\n";
		$pos++;
    }
    print qq
|
function resetColor(defaultColor) {
    var colorDisplay = new jscolor.color( document.getElementById( 'dispColor' ));
    colorDisplay.fromString( defaultColor );
}
function resetSpec() {
    for (var i=0 ; i < specAA.length ; i++) {
		var myB=document.getElementById('aaspec'+specAA[i][0]);
		if (!(myB.classList.contains(specAA[i][1]))) {
			if (myB.classList.contains('pushedButton')) {
				myB.classList.add('releasedButton');
				myB.classList.remove('pushedButton');
			}
			else{
				myB.classList.add('pushedButton');
				myB.classList.remove('releasedButton');
			}
		}
    }
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR>
<FONT class="title">Editing modification <FONT color='#DD0000'>$modifName</FONT></FONT>
<FORM name="modifForm" method="POST">
<INPUT type="hidden" name="ACT" value="list">
<INPUT type="hidden" name="LET" value="$startLetter">
</FORM>
<FORM name="editModForm" method=POST onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="idModification" value="$modificationID">
<INPUT type="hidden" name="specificity" value="$spec">
<INPUT type="hidden" name="LET" value="$startLetter">

$hiddenSpecString<BR>
<TABLE border=0 cellpadding=2 bgcolor=$dark>
<TR><TH align=right nowrap>PSI-MS name :</TH><TD nowrap bgcolor=$light>&nbsp;$psiName</TD></TR>
<TR><TH align=right nowrap>Interim name :</TH><TD nowrap bgcolor=$light>&nbsp;$interName</TD></TR>
<TR><TH align=right valign=top nowrap>&nbsp;Alternative name(s) :</TH><TD nowrap bgcolor=$light>&nbsp;$altNameString</TD></TR>
<TR><TH align=right valign=top>Description :</TH><TD nowrap bgcolor=$light><TEXTAREA name='des' rows='2' cols='65'>$des</TEXTAREA></TD></TR>
<TR><TH align=right>Monoisotopic :</TH><TD nowrap bgcolor=$light><INPUT type="text" name="mono" value="$mono" size=10></TD></TD></TR>
<TR><TH align=right>Average :</TH><TD nowrap bgcolor=$light><INPUT type="text" name="avg" value="$avg" size=10></TD></TR>
<TR><TH align=right>&nbsp;Unimod accession :</TH><TD nowrap bgcolor=$light>&nbsp;$unimodID</TD></TR>
<TR><TH align=right valign=top nowrap>Specificity :<INPUT type="button" id="showSpec" class="font11" value="Edit specificity" onclick="showHideSpecificity('show')"/><INPUT type="button" id="hideSpec" class="font11" value="Hide Specificity Editing" style="display:none" onclick="showHideSpecificity('hide')"/></TH><TD valign=top nowrap bgcolor=$light>&nbsp;$specString <BR>
<DIV id="specDIV" style="display:none">$editSpecString</DIV></TD></TR>
<TR><TH align=right>Project display :</TH><TD valign=top nowrap bgcolor=$light>&nbsp;<B>-Set code:</B><INPUT type="text" name="code" value="$dispCode" size=3>&nbsp;<B>-Choose color:</B><INPUT id="dispColor" name="dispColor" class="color" value="$dispColor" autocomplete="on" style="background-image: none; background-color: rgb(41, 255, 191); color: rgb(0, 0, 0);"><INPUT type="button" value="Reset color" onclick="resetColor('$dispColor')"></TD></TR>
<TR><TH align=right>Is label :</TH><TD nowrap bgcolor=$light>&nbsp;$isLabelStr</TD></TR>
<TR><TH align=right>Is substitution :</TH><TD nowrap bgcolor=$light>&nbsp;$isSubstStr</TD></TR>
|;
	if ($action eq 'edit') {
		print qq
|<TR><TH align=right>Merge with :</TH><TD nowrap bgcolor=$light>
<SELECT name='mergeMod'>
<OPTION value=''> -= Select =- </OPTION>
|;
		$sthDetNames->execute;
		while (my ($modID,$psiMsName,$interimName,$altNames) = $sthDetNames->fetchrow_array) {
			if ($altNames) {
				$altNames=~ s/##/,/g;
				$altNames=substr $altNames,1;
			}
			my $stringName=($psiMsName)? $psiMsName :($interimName)? $interimName : $altNames;
			print "<OPTION value='$modID'>$stringName</OPTION>\n";
		}
		$sthDetNames->finish;
		print "</SELECT></TD></TR>\n";
	}
    print qq
|<TR><TD colspan=2 align=center>
<INPUT type="submit" name="submit" value=" Save ">
&nbsp;&nbsp;&nbsp;<INPUT type="reset" value="  Cancel changes  " onclick="resetSpec()">
&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="document.modifForm.submit();">
</TD></TR>
</TABLE>
</FORM>
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

    $dbh->disconnect;
}
elsif ($action eq 'delete') {

}

sub processForm {
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	my $modificationID=param('idModification');
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
		my ($alNamesOld,$specOld)=$dbh->selectrow_array("SELECT SYNONYMES,SPECIFICITY FROM MODIFICATION WHERE ID_MODIFICATION=$modificationID");
		my ($psiMsName,$interiName,$alNamesNew,$specNew)=$dbh->selectrow_array("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES,SPECIFICITY FROM MODIFICATION WHERE ID_MODIFICATION=$newModID");
		my @psiLetters=($psiMsName)? split(//,$psiMsName) : ($interiName)? split(//,$interiName) : 'A' ;
		$startLetter=$psiLetters[0];
		if ($alNamesNew) {
			foreach my $syn (split(/##/,$alNamesOld)) {
				next unless $syn;
				$alNamesNew.="$syn##" unless $alNamesNew=~/##$syn##/;
			}
		}
		else {$alNamesNew=$alNamesOld;}

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
		$dbh->do("UPDATE MODIFICATION SET SYNONYMES='$alNamesNew',SPECIFICITY='$newSpec' WHERE ID_MODIFICATION=$newModID");

		print " Done.</FONT><BR>\n";

		###<Update Project / Analysis_Modification / Query_Modification

		###<PROJECT_MODIFICATION, ANALYSIS_MODIFICATION and SPECTRUM_MODIFICATION -> same way
		print "<BR><FONT class='title3'>Updating modification data for Project items...";
		my $count=0;
		foreach my $item ('PROJECT','ANALYSIS','SPECTRUM') {
			my $getItemID=$dbh->prepare("SELECT COUNT(ID_$item),ID_$item,ID_MODIFICATION FROM ${item}_MODIFICATION WHERE ID_MODIFICATION=$modificationID OR ID_MODIFICATION=$newModID GROUP BY ID_$item");
			my $sthUpItemmod=$dbh->prepare("UPDATE ${item}_MODIFICATION SET ID_MODIFICATION=$newModID WHERE ID_$item=? AND ID_MODIFICATION=$modificationID");
			$getItemID->execute;
			while (my ($numMod,$itemID,$modID)=$getItemID->fetchrow_array) {
				if ($numMod==1 && $modID == $modificationID) {
					$sthUpItemmod->execute($itemID);
				}
				$count++;
				if ($count==100) {
					print '.';
					$count=0;
				}
			}
			$getItemID->finish;
			$sthUpItemmod->finish;
			$dbh->do("DELETE FROM ${item}_MODIFICATION WHERE ID_MODIFICATION=$modificationID");
		}
		print " Done.</FONT><BR>\n";

		###<QUERY_MODIFICATION and PEPTIDE_MODIFICACTION -> same way
		print "<BR><FONT class='title3'>Updating modification data for Peptides...";
		$count=0;
		foreach my $item ('QUERY','PEPTIDE') {
			my $extString=($item eq 'QUERY')? ',PEP_RANK' : '';
			my $sthItemInfo=$dbh->prepare("SELECT ID_MODIFICATION,ID_$item,POS_STRING,REF_POS_STRING$extString FROM ${item}_MODIFICATION WHERE ID_MODIFICATION=$modificationID OR ID_MODIFICATION=$newModID"); # ??? WHY $newModID ??? (PP 08/02/17)
			$sthItemInfo->execute;
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
			$dbh->do("DELETE FROM ${item}_MODIFICATION WHERE ID_MODIFICATION=$modificationID");
		}
		print " Done.</FONT><BR>\n";

		$dbh->do("DELETE FROM MODIFICATION WHERE ID_MODIFICATION=$modificationID");

		sleep 2;
	}

	####<Edit or Add a modification>####
	else {
		my $mono=param('mono');
		my $avg=param('avg');
		#my $otherParams=(param('des'))? ",DES='".param('des')."'" : '';
		#$otherParams.=(defined(param('code')) && param('code') ne '')? ",DISPLAY_CODE='".param('code')."'" : '';
		#$otherParams.=(defined(param('dispColor')))? ",DISPLAY_COLOR='".param('dispColor')."'" : '';
		my ($code,$color)=(param('code') && param('dispColor'))? (param('code'),param('dispColor')) : (undef,undef);
		#$otherParams.=(param('isSubst') eq 'y')? ",IS_SUBST=1": ",IS_SUBST=0";
		#$otherParams.=(param('isLabel') eq 'y')? ",IS_LABEL=1": ",IS_LABEL=0";
		#$otherParams.=',IS_SUBST='.param('isSubst').',IS_LABEL='.param('isLabel');
		###> Check specificity;
		my @oldSpec=split(/,+/,param('specificity')) ; # To keep the context (ex: Trimethyl -> A (Protein N-term) only!!! <-> A;-)
		my @AA=('A','C'..'I','K'..'N','P'..'T','V'..'Y','nterm','cterm','pnterm','pcterm'); # added X
		my $newSpec='';
		foreach my $aa (@AA) {
			if(param("isClickedaaspec$aa")) {
				my $tag=1;
				foreach my $aaOld (@oldSpec) {
					next unless $aaOld;
					my $aa2=($aa eq 'nterm')?'=':($aa eq 'pnterm')?'-':($aa eq 'cterm')?'\*':($aa eq 'pcterm')?'\+':$aa;
					if ($aaOld =~ /^$aa2/) {
						$tag=0;
						$newSpec.="$aaOld,";
						last;
					}
				}
				$aa=($aa eq 'nterm')?'=':($aa eq 'pnterm')?'-':($aa eq 'cterm')?'*':($aa eq 'pcterm')?'+':$aa;
				if ($tag) {
					$newSpec.=',' if $newSpec;
					$newSpec.=$aa;
				}
			}
		}
		$newSpec=undef unless $newSpec;
		#$dbh->do("UPDATE MODIFICATION SET MONO_MASS=$mono, AVGE_MASS=$avg $otherParams WHERE ID_MODIFICATION=$modificationID");
		my $sthUpMod=$dbh->prepare("UPDATE MODIFICATION SET DES=?,MONO_MASS=?,AVGE_MASS=?,DISPLAY_CODE=?,DISPLAY_COLOR=?,IS_SUBST=?,IS_LABEL=?,SPECIFICITY=?,VALID_STATUS=? WHERE ID_MODIFICATION=$modificationID");
		my $validStatus=($mono && $avg && $newSpec)? 1 : 0;
		$sthUpMod->execute(param('des'),$mono,$avg,$code,$color,param('isSubst'),param('isLabel'),$newSpec,$validStatus);
		$sthUpMod->finish;
	}

	$dbh->commit;
	$dbh->disconnect;

	print qq
|<SCRIPT type="text/javascript">
window.location="./manageModifications.cgi?ACT=list&LET=$startLetter";
</SCRIPT>
</BODY>
</HTML>
|;
}

sub createEditSpecificityString {
    my ($dbh,$modificationID)=@_;
    my @AA=('A','C'..'I','K'..'N','P'..'T','V'..'Y'); # added X
    my %specDef;
    foreach my $aa (@AA) {
		#next unless $aa =~ /[A-Z]/;
		$specDef{$aa}{'DEFINITION'}='releasedButton';
		$specDef{$aa}{'IS_CLICKED'}=0;
    }
    $specDef{'nterm'}{'DEFINITION'}=$specDef{'pnterm'}{'DEFINITION'}=$specDef{'cterm'}{'DEFINITION'}=$specDef{'pcterm'}{'DEFINITION'}=0;
    $specDef{'nterm'}{'IS_CLICKED'}=$specDef{'pnterm'}{'IS_CLICKED'}=$specDef{'cterm'}{'IS_CLICKED'}=$specDef{'pcterm'}{'IS_CLICKED'}=0;
    ###> Specificity in Unimod/myProMS definition (could be modified)
    my ($specStr)=$dbh->selectrow_array("SELECT SPECIFICITY FROM MODIFICATION WHERE ID_MODIFICATION=$modificationID");
    if ($specStr) {
		foreach my $aa (split(/,+/,$specStr)) {
			$aa=($aa eq '=')?'nterm':($aa eq '-')?'pnterm':($aa eq '*')?'cterm':($aa eq '+')?'pcterm':$aa;
			if ($aa =~/;/) {
				my ($aaOnly,my $context)=split(/;/,$aa);
				my $contextString=($context eq '=')?'N-term':($context eq '-')?'Protein N-term':($context eq '*')?'C-term':'Protein C-term';
				$specDef{$aaOnly}{'POPUPINFO'}=" onmouseover=\"popup('<B>$contextString</B>')\" onmouseout=\"popout()\" ";
				#$specDef{$aaOnly}{'ANALYSIS'}=-1;
				$specDef{$aaOnly}{'DEFINITION'}='pushedButton';
				$specDef{$aaOnly}{'IS_CLICKED'}=1;
			}
			else {
				$specDef{$aa}{'DEFINITION'}='pushedButton';
				$specDef{$aa}{'IS_CLICKED'}=1;
			}
		}
    }
    ###> Specificity in projects (no modifications allowed)
    my $sthspecInAna=$dbh->prepare("SELECT DISTINCT(SPECIFICITY) FROM ANALYSIS_MODIFICATION WHERE ID_MODIFICATION=$modificationID");
    $sthspecInAna->execute;
    while (my ($specString) = $sthspecInAna->fetchrow_array) {
		if ($specString) {
			foreach my $aa (split(/,/,$specString)) {
				next if $aa eq '?'; # ambiguity
				$aa=($aa eq '=')?'nterm':($aa eq '-')?'pnterm':($aa eq '*')?'cterm':($aa eq '+')?'pcterm':$aa;
				$specDef{$aa}{'ANALYSIS'}=-1; # set as used by Analyses
			}
		}
    }
    $sthspecInAna->finish;
    my $numRes=0;
    my $editSpecString="<TABLE>\n";
    my $hiddenSpecString="";

    ###> N-TER SIDE
    my $disNterm=($specDef{"nterm"}{'ANALYSIS'} && $specDef{"nterm"}{'ANALYSIS'}==-1)?'disabled':'';
    my $disPNterm=($specDef{"pnterm"}{'ANALYSIS'} && $specDef{"pnterm"}{'ANALYSIS'}==-1)?'disabled':'';
    $editSpecString.="<INPUT name=\"aaspec\" id=\"aaspecnterm\" style=\"width: 90px\" type=\"button\" class=\"$specDef{nterm}{DEFINITION}\" value=\"Any N-term\" $disNterm onclick=\"changeAASpec(this)\">";
    $editSpecString.="<INPUT name=\"aaspec\" id=\"aaspecpnterm\" style=\"width: 110px\" type=\"button\" class=\"$specDef{pnterm}{DEFINITION}\" value=\"Protein N-term\" $disPNterm onclick=\"changeAASpec(this)\">";

    foreach my $aa (sort{$a cmp $b} @AA) {
		$specDef{$aa}{'POPUPINFO'}="onClick=\"changeAASpec(this)\"" unless $specDef{$aa}{'POPUPINFO'};
		$numRes++;
		my $disabled=($specDef{$aa}{'ANALYSIS'} && $specDef{$aa}{'ANALYSIS'}==-1) ? 'disabled':'';
		$editSpecString.="<INPUT name=\"aaspec\" id=\"aaspec$aa\" style=\"width: 30px\" type=\"button\" class=\"$specDef{$aa}{DEFINITION}\" value=\"$aa\" $disabled $specDef{$aa}{POPUPINFO}>";
		if ($numRes==11){
			$editSpecString.="<BR>";
			$numRes=-1;
		}
		$hiddenSpecString.="<INPUT type=\"hidden\" name=\"isClickedaaspec$aa\" id=\"isClickedaaspec$aa\" value=\"$specDef{$aa}{IS_CLICKED}\">\n";
    }

    ###> C-TER SIDE
    my $disCterm=($specDef{"cterm"}{'ANALYSIS'} && $specDef{"cterm"}{'ANALYSIS'}==-1)?'disabled':'';
    my $disPCterm=($specDef{"pcterm"}{'ANALYSIS'} && $specDef{"pcterm"}{'ANALYSIS'}==-1)?'disabled':'';
    $editSpecString.=qq
|<INPUT name="aaspec" id="aaspecpcterm" style="width:110px" type="button" class="$specDef{pcterm}{DEFINITION}" value="Protein C-term" $disPCterm onclick="changeAASpec(this)">
<INPUT name="aaspec" id="aaspeccterm" style="width:90px" type="button" class="$specDef{cterm}{DEFINITION}" value="Any C-term" $disCterm onclick="changeAASpec(this)">
</TR>
</TABLE>
|;

    ###> Hidden values
    $hiddenSpecString.=qq
|<INPUT type="hidden" name="isClickedaaspecnterm" id="isClickedaaspecnterm" value="$specDef{nterm}{IS_CLICKED}">
<INPUT type="hidden" name="isClickedaaspecterm" id="isClickedaaspecterm" value="$specDef{cterm}{IS_CLICKED}">
<INPUT type="hidden" name="isClickedaaspecpnterm" id="isClickedaaspecpnterm" value="$specDef{pnterm}{IS_CLICKED}">
<INPUT type="hidden" name="isClickedaaspecpcterm" id="isClickedaaspecpcterm" value="$specDef{pcterm}{IS_CLICKED}">
|;

    return($editSpecString,$hiddenSpecString,%specDef);
}

####>Revision history<####
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
