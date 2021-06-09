#!/usr/local/bin/perl -w

################################################################################
# proteinMining.cgi       1.0.4	                                               #
# Authors: P. Poullet, S.Liva (Institut Curie)	                               #
# Contact: myproms@curie.fr                                                    #
# Search for identifier and provide informations and peptide coverage          #
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
use POSIX qw(strftime); # to get the time
my %promsPath=&promsConfig::getServerInfo;
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my $dbh=&promsConfig::dbConnect;

my $ajax=param('AJAX') || '';

if ($ajax eq 'ajaxSearchProteins') {
    &ajaxSearchProteins(param('partialMatch'));
    exit;
}
elsif ($ajax eq 'ajaxDisplayEvidences') {
    &ajaxDisplayEvidences;
    exit;
}
elsif ($ajax eq 'ajaxDisplayMatchingProteins') {
    &ajaxDisplayMatchingProteins;
    exit;
}

my %species;
my $sthSelSpecies=$dbh->prepare("SELECT ID_SPECIES, SCIENTIFIC_NAME FROM SPECIES WHERE IS_REFERENCE=1");
$sthSelSpecies->execute;
while (my ($speciesID, $scientifName)=$sthSelSpecies->fetchrow_array) {
    $species{$speciesID}=$scientifName;
}

print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Protein Data Mining</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE>
/* FONT.notCompatible {text-decoration:line-through;} */
.dark {
    background-color: $darkColor;
    opacity: .999;
}
.light {
    background-color: $lightColor;
    opacity: .999;
}
</STYLE>
<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/peptidePlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
<SCRIPT LANGUAGE="JavaScript">
|;
&promsMod::popupInfo();
print qq
|var protIdStrg; // needed by ajaxDisplayMatchingProteins
var includeHiddenProt=0; // needed by all ajax calls
var selectedSpeciesID; // needed by ajaxDisplayMatchingProteins
function ajaxSearchProteins() {
	var myForm=document.searchProtForm;
    var searchString=myForm.searchString.value;
    var partial=(myForm.partialMatch.checked)? 1 : 0;
	includeHiddenProt=(myForm.hidden.checked)? 1 : 0;
	var species=myForm.species.value;

    if (!searchString) {
		alert('Please provide a protein/gene identifier');
		return;
    }
    if (!species) {
		alert('Please provide a species name');
		return;
    }

	document.getElementById('selectedProtDIV').innerHTML='';
	document.getElementById('matchingProtDIV').innerHTML='';
	var displayDiv=document.getElementById('searchProtDIV');
	displayDiv.style.display='';
	displayDiv.innerHTML='<BR><BR><BR><BR><BR><FONT class=\"title3\">Fetching data. Please wait...</FONT><BR><IMG src=\"$promsPath{images}/scrollbarGreen.gif\">';

    var params="AJAX=ajaxSearchProteins&searchString="+searchString+"&species="+species+"&partialMatch="+partial+"&hidden="+includeHiddenProt;

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
    XHR.open("POST","./proteinMining.cgi",true);
    XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=utf-8");
    XHR.setRequestHeader("Content-length", params.length);
    XHR.setRequestHeader("Connection", "close");
    XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			displayDiv.innerHTML = XHR.responseText;
			modifProtSelection=false;
		}
    }
    XHR.send(params);
}
var modifProtSelection=false;
function checkProtGroup(firstIdx,lastIdx,chkStatus) {
	var protStrg=document.getElementsByName('protChk');
	if (protStrg.length) {
		for (var i=firstIdx; i<=lastIdx; i++) {
			protStrg[i].checked=chkStatus;
		}
	}
	else {// only 1 checkbox
		protStrg.checked=chkStatus;
	}
	modifProtSelection=true;
}

function masterCheck(mastID) {
	modifProtSelection=true;
return; /*!!! protein group exclusion disabled !!!*/
    var protStrg=document.getElementsByName('protChk');
	if (!protStrg.length) return; // only 1 checkbox
    var numchk=0;
    for (var i=0; i<protStrg.length; i++) {
		if (protStrg[i].checked) {
			numchk++;
		}
		if (protStrg[i].id.match(mastID+':')) {
			protStrg[i].disabled=false;
		}
		else {
			protStrg[i].checked=false;
			protStrg[i].disabled=true;
			document.getElementById('protein:'+i).className='notCompatible';
		}
    }
    if (!numchk){
		for (var i=0; i<protStrg.length; i++) {
			protStrg[i].disabled=false;
			document.getElementById('protein:'+i).className='';
		}
    }
}
function ajaxDisplayEvidences(speciesID) {
	if (!modifProtSelection) {
		document.getElementById('searchProtDIV').style.display='none';
		document.getElementById('selectedProtDIV').style.display='';
		document.getElementById('matchingProtDIV').style.display='';
		return;
	}
	selectedSpeciesID=speciesID;
	selectedPepSeq=null;
    var protStrg=document.getElementsByName('protChk');
    var protArray=new Array();
    var isChk=false;
	if (protStrg.length) {
		for (var i=0; i<protStrg.length; i++) {
			if (protStrg[i].checked) {
				isChk=true;
				protArray.push(protStrg[i].value);
			}
		}
	}
	else if (protStrg.checked) { // only 1 checkbox
		isChk=true;
		protArray.push(protStrg.value);
	}
    if (!isChk) {
		alert('Select at leat one protein');
		return;
    }

	modifProtSelection=false;

	document.getElementById('searchProtDIV').style.display='none';
	document.getElementById('matchingProtDIV').innerHTML='';
	var displayDiv=document.getElementById('selectedProtDIV');
	displayDiv.style.display='';
	displayDiv.scrollIntoView(1);
	displayDiv.innerHTML='<BR><BR><BR><BR><BR><FONT class=\"title3\">Fetching data. Please wait...</FONT><BR><IMG src=\"$promsPath{images}/scrollbarGreen.gif\">';

    protIdStrg=protArray.join(',');
    var params="AJAX=ajaxDisplayEvidences&speciesID="+speciesID+"&species="+document.searchProtForm.species.value+"&hidden="+includeHiddenProt+"&protInfo="+protIdStrg; // species param in case speciesID=0

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
    XHR.open("POST","./proteinMining.cgi",true);
    XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=utf-8");
    XHR.setRequestHeader("Content-length", params.length);
    XHR.setRequestHeader("Connection", "close");
    XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			var codeParts=XHR.responseText.split('#==========#');
			displayDiv.innerHTML=codeParts[0]; // HTML part
			eval(codeParts[1]); // javascript part
		}
    }
    XHR.send(params);
}
var selectedPepSeq;
function prepareAjaxDisplayMatchingProteins(pepIdStrg,pepObj) { // called from peptidePlot
	ajaxDisplayMatchingProteins(pepObj.sequence,pepIdStrg);
}
function ajaxDisplayMatchingProteins(pepSeq,pepIdStrg) {
	if (selectedPepSeq) document.getElementById(selectedPepSeq).className='selectable';
	document.getElementById(pepSeq).className='selected';
	selectedPepSeq=pepSeq;
	var displayDiv=document.getElementById('matchingProtDIV');
	displayDiv.style.display='';
	displayDiv.innerHTML='<BR><BR><BR><BR><BR><FONT class=\"title3\">Fetching data. Please wait...</FONT><BR><IMG src=\"$promsPath{images}/scrollbarGreen.gif\">';
	displayDiv.scrollIntoView(1);
	var params="AJAX=ajaxDisplayMatchingProteins&speciesID="+selectedSpeciesID+"&species="+document.searchProtForm.species.value+"&hidden="+includeHiddenProt+"&pepSeq="+pepSeq+"&pepInfo="+pepIdStrg+"&protInfo="+protIdStrg;

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
    XHR.open("POST","./proteinMining.cgi",true);
    XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=utf-8");
    XHR.setRequestHeader("Content-length", params.length);
    XHR.setRequestHeader("Connection", "close");
    XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			displayDiv.innerHTML=XHR.responseText;
		}
    }
    XHR.send(params);
}
/*
function displayGraphicalView() {
    var dispDiv = (document.getElementById('graphicalDiv').style.display=='none')? 'block' : 'none';
    document.getElementById('graphicalDiv').style.display=dispDiv;
}
*/
var XHR=null;
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
</SCRIPT>
<BODY background="$promsPath{images}/bgProMS.gif">
<FORM method="post" name="searchProtForm">
<DIV style="float:top">
<BR>
<TABLE><TR><TH bgcolor="$darkColor">
<FONT class="title2">&nbsp;Go to:</FONT><SELECT class="title action" onchange="if (this.value) {window.location='$promsPath{cgi}/'+this.value;}">
	<OPTION value="">-= Select =-</OPTION>
	<OPTION value="promsMain.cgi">Main Window</OPTION>
	<OPTION value="selectProject.cgi">Project Selection</OPTION>
</SELECT>
</TH></TR></TABLE>
</DIV>
<CENTER>
<FONT class="title">Protein Identification Evidences</FONT>
<BR>(Find a protein of interest and display identified peptides statistics)
<BR><BR>
<TABLE border=0 bgcolor="$darkColor">
<TR>
	<TH align="right" valign="top" width="120px" class="title2">Search for :</TH>
	<TD bgcolor="$lightColor"><INPUT type="text" name="searchString"  value="" size="35" placeholder="Protein/gene identifier or keyword"/><BR>
		<FONT style="font-size:11px;font-weight:normal;">Match is <B>not</B> case sensitive. Non-[a-z 0-9] are ignored.</FONT><BR>
	    <INPUT type="checkbox" name="partialMatch"/>&nbsp;Partial word match.<BR>
		<INPUT type="checkbox" name="hidden"/>&nbsp;Include hidden proteins.
	</TD>
</TR>
<TR>
	<TH align="right" class="title2">Species :</TH>
	<TD bgcolor="$lightColor"><INPUT type="text" name="species" value="" size="25" placeholder="Scientific name" list="speciesList"/><BR>
	<DATALIST id="speciesList">
|;
	foreach my $specID (sort{$species{$a} cmp $species{$b}} keys %species) {
	    print "<OPTION value=\"$species{$specID}\">\n";
	}
	print qq
|	</DATALIST>
</TR>
<TR>
	<TH colspan="2"><INPUT type="button" class="title3" value="Search" onclick="ajaxSearchProteins()"></TH>
</TR>
</TABLE>
<BR><BR>
<DIV id="searchProtDIV" style="display:none"></DIV>
<DIV id="selectedProtDIV" style="display:none"></DIV>
<DIV id="matchingProtDIV" style="display:none"></DIV>
</BODY>
</FORM>
</CENTER>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup()</SCRIPT>
</HTML>
|;

sub ajaxSearchProteins {

    print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);

	my $searchString=param('searchString') || '';
	$searchString=~s/['",\.\*;\(\)\[\]]//g; # ignoring these characters
	$searchString = lc $searchString;
	my $speciesName=param('species') || '';
	my $includeHiddenProt=param('hidden') || 0;
    my ($partial)=@_;

    my ($speciesID)=$dbh->selectrow_array("SELECT ID_SPECIES FROM SPECIES WHERE SCIENTIFIC_NAME='$speciesName'");
	$speciesID=0 unless $speciesID;
	my %geneEntryCode;
	my $sthMap=$dbh->prepare("SELECT CODE,ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN' OR CODE='ID'");
	$sthMap->execute;
	while (my ($code,$identifier)=$sthMap->fetchrow_array) {$geneEntryCode{$code}=$identifier;}
	$sthMap->finish;
	my ($searchQuery1,$searchQuery2)=($partial)? (" LIKE '%$searchString%'"," LIKE '%$searchString%'") : ("='$searchString'"," REGEXP '(^$searchString\[=\\s\\\\-_;,\]|\[=\\s\\\\-_;,\]$searchString\[=\\s\\\\-_;,\]|\[=\\s\\\\-_;,\]$searchString\$)'"); # [[:<:]]word[[:>:]] word boundaries -> Not compatible with mysql8
	my $protVisQuery=($includeHiddenProt)? '' : 'AND VISIBILITY > 0';

    my (%protein,%masterIdent,%gene,%des,%date,%masterProtOrder);
	if ($speciesID) {
		my $sthSelMasterID=$dbh->prepare("SELECT DISTINCT MI.ID_MASTER_PROTEIN,PROT_DES,UPDATE_DATE FROM MASTER_PROTEIN M
											INNER JOIN MASTERPROT_IDENTIFIER MI ON M.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN
											WHERE M.ID_SPECIES=$speciesID
											AND (MI.VALUE$searchQuery1 OR PROT_DES$searchQuery2)");
		my $sthSelProt=$dbh->prepare("SELECT P.ID_PROTEIN,ALIAS,IDENTIFIER,PROT_LENGTH,MW,PROT_DES,MAX(VISIBILITY) FROM PROTEIN P,ANALYSIS_PROTEIN AP WHERE P.ID_PROTEIN=AP.ID_PROTEIN AND ID_MASTER_PROTEIN=? $protVisQuery GROUP BY P.ID_PROTEIN");
		my $sthSelGene=$dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN=? AND ID_IDENTIFIER=? ORDER BY IDENT_RANK");
		$sthSelMasterID->execute;
		while (my ($masterID,$masterDes,$updateDate)=$sthSelMasterID->fetchrow_array) {
			$sthSelProt->execute($masterID);
			while (my ($protID,$alias,$ident,$protLength,$mw,$protDes,$bestVis)=$sthSelProt->fetchrow_array) {
				$protein{$masterID}{$ident}{$protID}=1;
				if ($masterIdent{$masterID} && $masterIdent{$masterID}{$ident}) {
					$masterIdent{$masterID}{$ident}[4]=$bestVis if $masterIdent{$masterID}{$ident}[4] < $bestVis; # keep best visibility
				}
				else {
					@{$masterIdent{$masterID}{$ident}}=($alias,$protLength,$mw,$protDes,$bestVis);
				}
			}
			next unless scalar keys %{$protein{$masterID}}; # no protein with selected visibility
			$des{$masterID}=$masterDes || 'No description';
			$updateDate=~s/ .+//; # trim hours
			$date{$masterID}=$updateDate;
			while (my ($code,$identID) = each %geneEntryCode) {
				$sthSelGene->execute($masterID,$identID);
				while (my ($identValue)=$sthSelGene->fetchrow_array) {
					push @{$gene{$masterID}{$code}},$identValue;
				}
			}
			$masterProtOrder{$masterID}=1;
		}
		$sthSelMasterID->finish;
		$sthSelProt->finish;
		$sthSelGene->finish;
	}

	my $sthNoMasterProt=$dbh->prepare("SELECT P.ID_PROTEIN,ALIAS,IDENTIFIER,PROT_LENGTH,MW,PROT_DES,MAX(VISIBILITY),ORGANISM FROM PROTEIN P,ANALYSIS_PROTEIN AP WHERE P.ID_PROTEIN=AP.ID_PROTEIN AND ID_MASTER_PROTEIN IS NULL AND (IDENTIFIER$searchQuery2 OR PROT_DES$searchQuery2) AND ORGANISM LIKE '%$speciesName%' $protVisQuery GROUP BY P.ID_PROTEIN");
	$sthNoMasterProt->execute;
	while (my ($protID,$alias,$ident,$protLength,$mw,$protDes,$bestVis,$species)=$sthNoMasterProt->fetchrow_array) {
		$protein{0}{$ident}{$protID}=1;
		if ($masterIdent{0} && $masterIdent{0}{$ident}) {
			$masterIdent{0}{$ident}[4]=$bestVis if $masterIdent{0}{$ident}[4] < $bestVis; # keep best visibility
		}
		else {
			@{$masterIdent{0}{$ident}}=($alias,$protLength,$mw,$protDes,$bestVis,$species);
		}
	}
	$sthNoMasterProt->finish;
	$masterProtOrder{0}=2;

    $dbh->disconnect;

	if (scalar keys %masterIdent) {
		print qq
|<TABLE border=0 cellpadding=2 cellspacing=0>
<TR><TD colspan="6"><INPUT type="button" class="title3" value="Display evidences" onclick="ajaxDisplayEvidences($speciesID)"></TD></TR>
<TR bgcolor="$darkColor">
<TH class="rBorder" colspan="2">&nbsp;Name in myProMS&nbsp;</TH>
<TH class="rBorder">&nbsp;Identifier&nbsp;</TH>
<TH class="rBorder">&nbsp;Length&nbsp;</TH>
<TH class="rBorder">&nbsp;MW <SMALL>(kDa)</SMALL>&nbsp;</TH>
<TH width="500" align="left">&nbsp;Description&nbsp;</TH></TR>
|;
		my $i=0;
		my $m=1;
		foreach my $masterID (sort{$masterProtOrder{$a}<=>$masterProtOrder{$b} || $a<=>$b} keys %masterIdent) { # no master after
			my $geneStrg='No matching genes found';
			if ($gene{$masterID}{GN}) {
				$gene{$masterID}{GN}[0]=$gene{$masterID}{GN}[0].'</B>';
				$geneStrg='<B>Gene: '.join(', ',@{$gene{$masterID}{GN}});
			}
			my $spaceStrg=($i)? '<FONT class="font11">&nbsp;</FONT>' : '';
			my $lastProtIdx=$i + (scalar keys %{$masterIdent{$masterID}}) - 1;

			my $proteinStrg=($masterID)? "<B>#$m. Protein: <A href=\"http://www.uniprot.org/uniprot/$gene{$masterID}{ID}[0]\" target=\"_blank\" onmouseover=\"popup('<B>Description:</B><BR>$des{$masterID}<BR><B>Click for details</B>')\" onmouseout=\"popout()\">$gene{$masterID}{ID}[0]</A></B>, $geneStrg (Updated $date{$masterID})" : '<B>No cross-reference data available</B>';
			print qq
|<TR><TH colspan="6">$spaceStrg</TH></TR>
<TR bgcolor="$darkColor"><TD colspan="6" align="left"><INPUT type="checkbox" onclick="checkProtGroup($i,$lastProtIdx,this.checked)">&nbsp;$proteinStrg</TD></TR>
|;
			my $bgColor=$lightColor;
			foreach my $ident (sort{$masterIdent{$masterID}{$b}[1] <=> $masterIdent{$masterID}{$a}[1] || $masterIdent{$masterID}{$a}[0] cmp $masterIdent{$masterID}{$b}[0]} keys %{$masterIdent{$masterID}}) { # decreasing length then alias
				my $alias=$masterIdent{$masterID}{$ident}->[0];
				my $protLength=$masterIdent{$masterID}{$ident}[1];
				my $mw=sprintf("%.1f",$masterIdent{$masterID}{$ident}[2]/1000);
				my $protDes=$masterIdent{$masterID}{$ident}[3];
				$protDes.=" <FONT class=\"org\">$masterIdent{$masterID}{$ident}[5]</FONT>" if ($masterID==0 && $masterIdent{$masterID}{$ident}[4] && $masterIdent{$masterID}{$ident}[4] ne $speciesName);
				my $joinProt=join(',', keys %{$protein{$masterID}{$ident}});
				my $tag=($masterIdent{$masterID}{$ident}[4] > 0)? 'TH' : 'TD';
				my $mStrg=($masterID)? '' : "#$m.&nbsp;";
				print qq
|<TR bgcolor="$bgColor" valign="top" class="list">
<TD><INPUT type="checkbox" id="$masterID:$i" name="protChk" value="$joinProt" onclick="masterCheck($masterID)"></TD>
<$tag align="left"><FONT id="protein:$i">$mStrg$alias</FONT>&nbsp;</$tag>
<TD>&nbsp;$ident&nbsp;</TD>
<TD align="center">&nbsp;$protLength&nbsp;</TD><TD align="right">&nbsp;$mw&nbsp;</TD><TD>$protDes</TD></TR>
|;
				$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
				$i++;
				$m++ if $masterID==0;
			}
			$m++;
		}
		print qq
|<TR><TD colspan="6"><INPUT type="button" class="title3" value="Display evidences" onclick="ajaxDisplayEvidences($speciesID)"></TD></TR>
</TABLE>
|;
	}
	else {
		print "<FONT class=\"title3\">No match found.</FONT>\n";
	}
	print "<BR><BR>\n";
    exit;
}

sub ajaxDisplayEvidences {
	print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);

	my $speciesID=param('speciesID');
	my $speciesName=param('species') || ''; # in case speciesID=0
	my $includeHiddenProt=param('hidden') || 0;
	my $strgProtIDList=param('protInfo');
#print "\$strgProtIDList=$strgProtIDList<BR>\n";
	my @protListID=split(',',$strgProtIDList);
	my (%protPepInfo,%pepByAnalysis,%displayProt,%peptideIDs,%pepSeq2ID,%analyses,%searchTypes);

	my $protVisQuery=($includeHiddenProt)? '' : 'AND VISIBILITY > 0';
	my $sthSelProt=$dbh->prepare("SELECT ALIAS,IDENTIFIER,PROT_LENGTH,MW,PROT_DES,MAX(VISIBILITY) FROM PROTEIN P,ANALYSIS_PROTEIN AP WHERE P.ID_PROTEIN=AP.ID_PROTEIN AND P.ID_PROTEIN=? $protVisQuery GROUP BY P.ID_PROTEIN");
	my $protLength=0;
	foreach my $protID (@protListID) {
		$sthSelProt->execute($protID);
		my ($alias, $identifier,$pLength,$mw,$protDes,$bestVis)=$sthSelProt->fetchrow_array;
		$pLength=0 unless $pLength;
		my $kd=sprintf("%.1f",$mw/1000);
		if ($displayProt{$identifier}) {
			$displayProt{$identifier}[4]=$bestVis if $displayProt{$identifier}[4] < $bestVis;
			$displayProt{$identifier}[5]++; # num protIDs
		}
		else {
			@{$displayProt{$identifier}}=($alias,$pLength,$kd,$protDes,$bestVis,1);
		}
		$protLength=$pLength if $protLength < $pLength;
	}
	$sthSelProt->finish;

	#>Fetching all peptides no matter what is the protein visibility in each Analysis
    my $sthSelPeptideInfo=$dbh->prepare("SELECT PPA.ID_PEPTIDE,PPA.ID_ANALYSIS,ABS(PPA.PEP_BEG),ABS(PPA.PEP_END),PPA.IS_SPECIFIC,P.PEP_SEQ,P.PEP_LENGTH,A.FILE_FORMAT,P.SCORE
											FROM PEPTIDE_PROTEIN_ATTRIB PPA
											INNER JOIN PEPTIDE P ON PPA.ID_PEPTIDE=P.ID_PEPTIDE
											INNER JOIN ANALYSIS A ON PPA.ID_ANALYSIS=A.ID_ANALYSIS
											WHERE PPA.ID_PROTEIN=? AND P.SCORE IS NOT NULL"); # AND P.MISS_CUT=0 (no filter on missed-cleavage because some selected proteins might be missed if identified only miss-cut peptides)

	my $protVisQuery2=($includeHiddenProt)? '' : 'INNER JOIN ANALYSIS_PROTEIN AP ON PPA.ID_PROTEIN=AP.ID_PROTEIN AND AP.VISIBILITY > 0';
	my $sthMasterProt=$dbh->prepare("SELECT COUNT(DISTINCT P.ID_MASTER_PROTEIN) FROM PEPTIDE_PROTEIN_ATTRIB PPA
										INNER JOIN PROTEIN P ON PPA.ID_PROTEIN=P.ID_PROTEIN
										$protVisQuery2
										INNER JOIN MASTER_PROTEIN MP ON P.ID_MASTER_PROTEIN=MP.ID_MASTER_PROTEIN AND MP.ID_SPECIES=$speciesID
										WHERE PPA.ID_PEPTIDE=? AND PPA.ID_PROTEIN NOT IN ($strgProtIDList)");

	my $sthProt=$dbh->prepare("SELECT COUNT(DISTINCT IDENTIFIER) FROM PEPTIDE_PROTEIN_ATTRIB PPA
										INNER JOIN PROTEIN P ON PPA.ID_PROTEIN=P.ID_PROTEIN AND ID_MASTER_PROTEIN IS NULL AND (ORGANISM LIKE '%$speciesName%' OR ORGANISM='unknown organism')
										$protVisQuery2
										WHERE ID_PEPTIDE=? AND PPA.ID_PROTEIN NOT IN ($strgProtIDList)");

	my $lastPepStart=0;
	foreach my $protID (@protListID) {
		$sthSelPeptideInfo->execute($protID);
		while (my ($peptideID,$anaID,$pepBeg,$pepEnd,$isSpecific,$pepSeq,$pepLength,$fileFormat,$score)=$sthSelPeptideInfo->fetchrow_array) {
			$analyses{$anaID}=1;
			$isSpecific = ($isSpecific)? 1 : 2;
			$score=0 unless $score;
			if (!$pepByAnalysis{$pepSeq} || !$pepByAnalysis{$pepSeq}{$anaID}) { # do only once / analysis & pepSeq
				$pepByAnalysis{$pepSeq}{$anaID}=1; # consider all previously selected proteins as 1
				my $count=0;
				if ($speciesID) {
					$sthMasterProt->execute($peptideID);
					($count)=$sthMasterProt->fetchrow_array;
					$pepByAnalysis{$pepSeq}{$anaID}+=$count;
				}
				$sthProt->execute($peptideID); # in case also no master prot mapped => count individual proteins
				($count)=$sthProt->fetchrow_array;
				$pepByAnalysis{$pepSeq}{$anaID}+=$count;
			}
			my $countProt=$pepByAnalysis{$pepSeq}{$anaID};
			$peptideIDs{$peptideID}=1;
			$pepSeq2ID{$pepSeq}{$anaID}=$peptideID; # only 1 pepID/analysis
			$fileFormat=~s/\..+//;
			$searchTypes{$fileFormat}=1;
			if ($protPepInfo{$pepSeq}) {
				$protPepInfo{$pepSeq}[2] = $countProt if ($countProt > $protPepInfo{$pepSeq}[2]);
				$protPepInfo{$pepSeq}[4]{$fileFormat} = $score if (!$protPepInfo{$pepSeq}[4]{$fileFormat} || $score < $protPepInfo{$pepSeq}[4]{$fileFormat});
				$protPepInfo{$pepSeq}[5]{$fileFormat} = $score if (!$protPepInfo{$pepSeq}[5]{$fileFormat} || $score > $protPepInfo{$pepSeq}[5]{$fileFormat});
			}
			else {
				my %scMin=($fileFormat=>$score);
				my %scMax=($fileFormat=>$score);
				@{$protPepInfo{$pepSeq}}=($pepBeg,$pepEnd,$countProt,$pepLength,\%scMin,\%scMax);
			}
			$lastPepStart=$pepBeg if $pepBeg > $lastPepStart;
		}
    }
	$sthSelPeptideInfo->finish;
	$sthMasterProt->finish;
	$sthProt->finish;

	my $countProtAnalysis=scalar keys %analyses;

    print qq
|<INPUT type="button" class="title3" value="Display search results" onclick="document.getElementById('searchProtDIV').style.display=''; document.getElementById('selectedProtDIV').style.display='none'; document.getElementById('matchingProtDIV').style.display='none';"><BR><BR>
<TABLE border=0 cellpadding=2 cellspacing=0>
<CAPTION align="top"><B>Selected Protein(s)</B></CAPTION>
<TR bgcolor="$darkColor">
<TH class="rbBorder">&nbsp;Name in myProMS&nbsp;</TH>
<TH class="rbBorder">&nbsp;Identifier&nbsp;</TH>
<TH class="rbBorder">&nbsp;Length&nbsp;</TH>
<TH class="rbBorder">&nbsp;MW <SMALL>(kDa)</SMALL>&nbsp;</TH>
<TH width=500 align="left" class="bBorder">&nbsp;Description</TH></TR>
|;
    my $bgColor=$lightColor;
    foreach my $identifier (sort{$a cmp $b} keys %displayProt) {
		my $tag=($displayProt{$identifier}[4])? 'TH' : 'TD';
		my $numIDstrg=($displayProt{$identifier}[5]==1)? '' : "&nbsp;x$displayProt{$identifier}[5]";
		print qq
|<TR valign="top" bgcolor="$bgColor">
	<$tag align="left">&nbsp;$displayProt{$identifier}[0]$numIDstrg&nbsp;</$tag>
	<TD>&nbsp;$identifier&nbsp;</TD>
	<TD align="center">&nbsp;$displayProt{$identifier}[1]&nbsp;</TD>
	<TD align="right">&nbsp;$displayProt{$identifier}[2]&nbsp;</TD>
	<TD>$displayProt{$identifier}[3]&nbsp;</TD>
</TR>
|;
		$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
    }
	my $protStrg=(scalar keys %displayProt ==1)? 'protein is' : 'proteins are';
	my $analysisStrg=($countProtAnalysis==1)? 'Analysis' : 'Analyses';
	my $numSeqFound=scalar keys %protPepInfo;
	my $numPepFound=scalar keys %peptideIDs;
	my $numSearchTypes=scalar keys %searchTypes;
    print qq
|</TABLE><BR>

<DIV id="pepPlotDiv"></DIV><BR>

<TABLE border=0  cellpadding=2 cellspacing=0>
<CAPTION align="top"><B>Selected $protStrg identified in $countProtAnalysis $analysisStrg with $numSeqFound peptide sequences ($numPepFound identifications)</B></CAPTION>
<TR bgcolor="$darkColor">
	<TH class="rbBorder" rowspan=2>&nbsp;#&nbsp;</TH>
	<TH class="rbBorder" rowspan=2>&nbsp;Start&nbsp;</TH>
	<TH class="rbBorder" rowspan=2>Sequence</TH>
	<TH class="rbBorder" rowspan=2>&nbsp;Length&nbsp;</TH>
	<TH class="rbBorder" rowspan=2>&nbsp;Specificity&nbsp;<BR>(%)</TH>
	<TH class="rbBorder" rowspan=2>&nbsp;Occurence&nbsp;<BR>(%)</TH>
	<TH class="bBorder" colspan=$numSearchTypes>&nbsp;Scores [min~max]&nbsp;</TH>
</TR>
<TR bgcolor="$darkColor">
|;
	my $typeCount=0;
	foreach my $fileFormat (sort keys %searchTypes) {
		my $class=(++$typeCount==$numSearchTypes)? 'bBorder' : 'rbBorder';
		print "<TH class=\"$class\">$fileFormat</TH>\n";
	}
	print "</TR>\n";

	my $trClass='light';
	my %jsData;
	my $pepCount=0;
	my $prevBeg=0;
    foreach my $pepSeq (sort{$protPepInfo{$a}[0] <=> $protPepInfo{$b}[0] || $protPepInfo{$a}[3] <=> $protPepInfo{$b}[3]} keys %protPepInfo) {
		$pepCount++;
		#<Occurence
		my $countPepAnalysis=scalar( keys %{$pepByAnalysis{$pepSeq}});
		my $occur=sprintf("%.1f",($countPepAnalysis/$countProtAnalysis)*100); $occur*=1;
		my $occurStrg=($occur==100)? '<B>'.$occur.'</B>' : $occur;
		#<Specificity
		my @pepIDs;
		foreach my $anaID (keys %{$pepSeq2ID{$pepSeq}}) {push @pepIDs,$pepSeq2ID{$pepSeq}{$anaID};}
		my $pepIDstrg=join(',',@pepIDs);
		my $spec=sprintf("%.1f", 100/$protPepInfo{$pepSeq}[2]); $spec*=1;
		my $specStrg=($spec==100)?'<B>'.$spec.'</B>' : $spec;
		#<Other
		my ($minScore,$maxScore,%minScores,%maxScores);
		foreach my $fileFormat (sort keys %{$protPepInfo{$pepSeq}[4]}) {
			$minScores{$fileFormat}=sprintf("%.2f", $protPepInfo{$pepSeq}[4]{$fileFormat});
			$maxScores{$fileFormat}=sprintf("%.2f", $protPepInfo{$pepSeq}[5]{$fileFormat});
			$minScore=$minScores{$fileFormat} if (!$minScore || $minScore > $minScores{$fileFormat});
			$maxScore=$maxScores{$fileFormat} if (!$maxScore || $maxScore < $maxScores{$fileFormat});
		}
		my $beg=$protPepInfo{$pepSeq}[0];
		my $pepLength=@{$protPepInfo{$pepSeq}}[3];
		@{$jsData{$pepSeq}}=($beg,$occur,$spec,"'$minScore~$maxScore'",$pepIDstrg);
		my $begStrg=($beg==$prevBeg)? '' : $beg;
		print qq
|<TR class="list $trClass">
	<TH align="right">&nbsp;$pepCount&nbsp;</TH>
	<TD align="right">$begStrg&nbsp;</TD>
	<TH align="left">&nbsp;<A id="$pepSeq" class="selectable" href="javascript:ajaxDisplayMatchingProteins('$pepSeq','$pepIDstrg')" onmouseover="popup('Click to display matching proteins')" onmouseout="popout()">$pepSeq</A>&nbsp;</TH>
	<TD align="center">&nbsp;$pepLength&nbsp;</TD>
	<TD><DIV class="barContainerRight"><DIV class="barElementRight barGreen" style="width:$spec%"></DIV><DIV class="barValue">$specStrg&nbsp;</DIV></DIV></TD>
	<TD><DIV class="barContainerRight"><DIV class="barElementRight barRed" style="width:$occur%"></DIV><DIV class="barValue">$occurStrg&nbsp;</DIV></DIV></TD>
|;
		foreach my $fileFormat (sort keys %searchTypes) {
			if (defined $minScores{$fileFormat}) {print "<TD align=\"center\">&nbsp;$minScores{$fileFormat}~$maxScores{$fileFormat}&nbsp;</TD>\n";}
			else {print "<TD align=\"center\">&nbsp;-&nbsp;</TD>\n";}
		}
		print "</TR>\n";
		$trClass=($trClass eq 'light')? 'dark' : 'light';
		$prevBeg=$beg;
    }
	print qq
|</TABLE>
<BR>
<INPUT type="button" class="title3" value="Display search results" onclick="document.getElementById('searchProtDIV').style.display=''; document.getElementById('selectedProtDIV').style.display='none'; document.getElementById('matchingProtDIV').style.display='none';"><BR><BR>
|;

	###<Javascript for peptide plot>###
	my $plotWidth=($protLength)? $protLength*1.5 : $lastPepStart+50;
	$plotWidth=400 if $plotWidth < 400;
	$plotWidth=1000 if $plotWidth > 1000;
    print qq
|#==========#
var PP=new peptidePlot({div:'pepPlotDiv',width:$plotWidth,height:300,valueAxisLabel:'Occurence (%)',valueLabel:'Occurence (%)',
						protein:{length:$protLength},
						peptideOnClick:prepareAjaxDisplayMatchingProteins,
						peptideProperties: ['Specificity (%)','Score','id'],
						//convertValueDisplayed: function(v) {return v+'%'},
						peptideColor: {map:'Specificity (%)',type:'continuous',range:[0,100]}
						});
|;
	foreach my $pepSeq (sort{$protPepInfo{$a}[0] <=> $protPepInfo{$b}[0] || $protPepInfo{$a}[3] <=> $protPepInfo{$b}[3]} keys %protPepInfo) {
		print "PP.addPeptide([$jsData{$pepSeq}[0],'$pepSeq',null,null,$jsData{$pepSeq}[1],null,$jsData{$pepSeq}[2],$jsData{$pepSeq}[3],'$jsData{$pepSeq}[4]']);\n";
	}
	print "PP.draw();\n";
    exit;
}


sub ajaxDisplayMatchingProteins {
	print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);

	my $speciesID=param('speciesID') || 0;
	my $speciesName=param('species') || ''; # in case species not found in DB (no speciesID)
	my $includeHiddenProt=param('hidden') || 0;
	my $pepSeq=param('pepSeq');
	my @pepIdList=split(',',param('pepInfo'));
	my %selProtIdList;
	foreach my $protID (split(',',param('protInfo'))) {$selProtIdList{$protID}=1;}

	my %geneEntryCode;
	my $sthMap=$dbh->prepare("SELECT CODE,ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN' OR CODE='ID'");
	$sthMap->execute;
	while (my ($code,$identifier)=$sthMap->fetchrow_array) {$geneEntryCode{$code}=$identifier;}
	$sthMap->finish;

	my (%matchingProt,%masterIdent,%gene,%des,%date,%masterProtOrder);

	my $sthMatchProt=$dbh->prepare("SELECT DISTINCT ID_PROTEIN FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_PEPTIDE=?"); # could be a multi-species Analysis: No filter on species at this step!!!
	foreach my $pepID (@pepIdList) {
		$sthMatchProt->execute($pepID);
		while (my ($protID)=$sthMatchProt->fetchrow_array) {
			$matchingProt{$protID}=1;
		}
	}
	$sthMatchProt->finish;
	my $protVisQuery=($includeHiddenProt)? '' : 'AND VISIBILITY > 0';
	my $sthProt=$dbh->prepare("SELECT P.ID_MASTER_PROTEIN,MP.ID_MASTER_PROTEIN,P.PROT_DES,P.UPDATE_DATE,IDENTIFIER,ALIAS,P.PROT_LENGTH,P.MW,P.PROT_DES,ORGANISM,MAX(DISTINCT AP.VISIBILITY)
								FROM PROTEIN P
								LEFT JOIN MASTER_PROTEIN MP ON P.ID_MASTER_PROTEIN=MP.ID_MASTER_PROTEIN AND MP.ID_SPECIES=$speciesID
								JOIN ANALYSIS_PROTEIN AP ON P.ID_PROTEIN=AP.ID_PROTEIN $protVisQuery
								WHERE P.ID_PROTEIN=? GROUP BY P.ID_PROTEIN"
							);
	my $sthSelGene=$dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN=? AND ID_IDENTIFIER=? ORDER BY IDENT_RANK");

	foreach my $protID (keys %matchingProt) {
		$sthProt->execute($protID);
		my ($masterID,$master2ID,$masterDes,$updateDate,$identifier,@data)=$sthProt->fetchrow_array;
		next unless $data[0]; # in case no match
		next if ($masterID && !$master2ID); # Bad species prot can pass the LEFT JOIN on MASTER_PROTEIN but MP.ID_MASTER_PROTEIN is then NULL unlike P.ID_MASTER_PROTEIN!
		$masterID=0 unless $masterID;
		next if ($masterID==0 && $data[4] !~ /unknown/ && $data[4] !~ /$speciesName/); # species text filter if no master prot
		#next if ($masterID==0 && $data[4] !~ /$speciesName/); # species
		if ($masterIdent{$masterID} && $masterIdent{$masterID}{$identifier}) {
			push @{$masterIdent{$masterID}{$identifier}[0]},$protID;
			$masterIdent{$masterID}{$identifier}[6]=$data[5] if $masterIdent{$masterID}{$identifier}[6] < $data[5]; # keep best visibility
		}
		else {
			@{$masterIdent{$masterID}{$identifier}}=([$protID],@data); # keep list of protIDs in case further mining to be implemented
		}
		if ($masterID) {
			unless ($date{$masterID}) { # 1st time master prot is seen
				$des{$masterID}=$masterDes || 'No description';
				$updateDate=~s/ .+//; # trim hours
				$date{$masterID}=$updateDate;
				while (my ($code,$identID) = each %geneEntryCode) {
					$sthSelGene->execute($masterID,$identID);
					while (my ($identValue)=$sthSelGene->fetchrow_array) {
						push @{$gene{$masterID}{$code}},$identValue;
					}
				}
				$masterProtOrder{$masterID}=1;
			}
		}
		else {
			$masterProtOrder{0}=2;
		}
	}
	$sthProt->finish;
	$sthSelGene->finish;
	$dbh->disconnect;

	print qq
|<TABLE border=0 cellpadding=2 cellspacing=0>
<TR bgcolor="$darkColor">
<TH class="rBorder" colspan="2">&nbsp;Name in myProMS&nbsp;</TH>
<TH class="rBorder">&nbsp;Identifier&nbsp;</TH>
<TH class="rBorder">&nbsp;Length&nbsp;</TH>
<TH class="rBorder">&nbsp;MW <SMALL>(kDa)</SMALL>&nbsp;</TH>
<TH width="500" align="left">&nbsp;Description&nbsp;</TH></TR>
|;
	my $i=0;
	my $m=1;
	foreach my $masterID (sort{$masterProtOrder{$a}<=>$masterProtOrder{$b} || $a<=>$b} keys %masterIdent) { # no master after
		my $geneStrg='No matching genes found';
		if ($gene{$masterID}{GN}) {
			$gene{$masterID}{GN}[0]=$gene{$masterID}{GN}[0].'</B>';
			$geneStrg='<B>Gene: '.join(', ',@{$gene{$masterID}{GN}});
		}
		my $spaceStrg=($i)? '<FONT class="font11">&nbsp;</FONT>' : '';

		my $proteinStrg=($masterID)? "<B>#$m. Protein: <A href=\"http://www.uniprot.org/uniprot/$gene{$masterID}{ID}[0]\" target=\"_blank\" onmouseover=\"popup('<B>Description:</B><BR>$des{$masterID}<BR><B>Click for details</B>')\" onmouseout=\"popout()\">$gene{$masterID}{ID}[0]</A></B>, $geneStrg (Updated $date{$masterID})" : '<B>No cross-reference data available</B>';
		print qq
|<TR><TH colspan="6">$spaceStrg</TH></TR>
<TR bgcolor="$darkColor"><TD colspan="6" align="left">$proteinStrg</TD></TR>
|;
		my $bgColor=$lightColor;
		foreach my $ident (sort{$masterIdent{$masterID}{$b}[2] <=> $masterIdent{$masterID}{$a}[2] || $masterIdent{$masterID}{$a}[1] cmp $masterIdent{$masterID}{$b}[1]} keys %{$masterIdent{$masterID}}) { # decreasing length then alias
			my $alias=$masterIdent{$masterID}{$ident}[1];
			my $protLength=$masterIdent{$masterID}{$ident}[2];
			my $mw=sprintf("%.1f",$masterIdent{$masterID}{$ident}[3]/1000);
			my $protDes=$masterIdent{$masterID}{$ident}[4];
			$protDes.=" <FONT class=\"org\">$masterIdent{$masterID}{$ident}[5]</FONT>" if ($masterID==0 && $masterIdent{$masterID}{$ident}[5] && $masterIdent{$masterID}{$ident}[5] ne $speciesName);
			my $tag=($masterIdent{$masterID}{$ident}[6])? 'TH' : 'TD';
			my $selProtImgStrg=''; # default
			foreach my $protID (@{$masterIdent{$masterID}{$ident}[0]}) {
				if ($selProtIdList{$protID}) { # 1 of preselected proteins
					$selProtImgStrg="<IMG src=\"$promsPath{images}/good.gif\">";
					last;
				}
			}
			my $mStrg=($masterID)? '' : "#$m.";
			print qq
|<TR bgcolor="$bgColor" valign="top" class="list">
<TD>$selProtImgStrg</TD>
<$tag align="left">$mStrg&nbsp;$alias&nbsp;</$tag>
<TD>&nbsp;$ident&nbsp;</TD>
<TD align="center">&nbsp;$protLength&nbsp;</TD><TD align="right">&nbsp;$mw&nbsp;</TD><TD>$protDes</TD></TR>
|;
			$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
			$i++;
			$m++ if $masterID==0;
		}
		$m++;
	}
	print "</TABLE><BR>\n";
	if ($m > 10) {
		print qq
|<INPUT type="button" class="title3" value="Display search results" onclick="document.getElementById('searchProtDIV').style.display=''; document.getElementById('selectedProtDIV').style.display='none'; document.getElementById('matchingProtDIV').style.display='none';">
|;
	}
	print "<BR><BR>\n";
	exit;
}

####>Revision history<####
# 1.0.4 [BUGFIX] Updated search query string (VS 02/04/2021)
# 1.0.3 [UPDATE] Changed RANK field to IDENT_RANK for compatibility with MySQL 8 (PP 04/03/20) 
# 1.0.2 Minor change (PP 12/01/17)
# 1.0.1 More data mining and display improvement (PP 31/05/16)
# 1.0.0 Created from SL proteinsMining (PP 04/01/16)
