#!/usr/local/bin/perl -w

################################################################################
# processAnalyses.cgi    1.4.4                                                 #
# Authors: P. Poullet, G. Arras, F. Yvon & M. Le Picard (Institut Curie)       #
# Contact: myproms@curie.fr                                                    #
# Generates list of options available to manage multiple analyses at once      #
# Displayed in resultFrame                                                     #
# Called by selectOption.cgi (not accessible to guest)                         #
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

# print header; warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};

####################
####>Parameters<####
####################
my $branchID=param('ID');
my ($item,$itemID)=split(':',$branchID);
$item=uc($item);

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

####>Fetching Project ID<####
my $projectID=&promsMod::getProjectID($dbh,$itemID,$item);

####>Fetching user information<####
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
my $projectFullAccess=($projectAccess =~ /bioinfo|mass|manag|super/)? 1 : 0;
my $deleteAccess=($projectAccess =~ /bioinfo|mass|manag/)? 1 : 0;

####>Fetching analyses status<####
#my ($anaQuery,$validQuery,$filterQuery,$numMaldi,$validMascotQuery); #$numUnValidDat,
my ($anaQuery,$filterQuery);
my $anaFieldStrg='ID_ANALYSIS,VALID_STATUS,INSTRUMENT,MS_TYPE,FILE_FORMAT,LABELING';
if ($item eq 'EXPERIMENT') {
	$anaQuery="SELECT $anaFieldStrg FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_EXPERIMENT=$itemID";
	#$validQuery="SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_EXPERIMENT=$itemID AND VALID_STATUS=?";
	$filterQuery="SELECT COUNT(*) FROM PROTEIN_VALIDATION,ANALYSIS,SAMPLE WHERE PROTEIN_VALIDATION.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SEL_STATUS=-3 AND ID_EXPERIMENT=$itemID";
	#($numMaldi)=$dbh->selectrow_array("SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_EXPERIMENT=$itemID AND INSTRUMENT='MALDI-TOF-TOF'");
	#$validMascotQuery="SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_EXPERIMENT=$itemID AND VALID_STATUS>=1 AND FILE_FORMAT='MASCOT.DAT'";
}
elsif ($item eq 'SAMPLE') {
	$anaQuery="SELECT $anaFieldStrg FROM ANALYSIS WHERE ID_SAMPLE=$itemID";
	#$validQuery="SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS WHERE ID_SAMPLE=$itemID AND VALID_STATUS=?";
	$filterQuery="SELECT COUNT(*) FROM PROTEIN_VALIDATION,ANALYSIS WHERE PROTEIN_VALIDATION.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND SEL_STATUS=-3 AND ID_SAMPLE=$itemID";
	#($numMaldi)=$dbh->selectrow_array("SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS WHERE ID_SAMPLE=$itemID AND INSTRUMENT='MALDI-TOF-TOF'");
	#$validMascotQuery="SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS WHERE ID_SAMPLE=$itemID AND VALID_STATUS>=1 AND FILE_FORMAT='MASCOT.DAT'";
}
elsif ($item eq 'GEL2D') {
	$anaQuery="SELECT $anaFieldStrg FROM SPOT,SAMPLE,ANALYSIS WHERE SPOT.ID_SPOT=SAMPLE.ID_SPOT AND SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND ID_GEL2D=$itemID";
	#$validQuery="SELECT COUNT(ID_ANALYSIS) FROM SPOT,SAMPLE,ANALYSIS WHERE SPOT.ID_SPOT=SAMPLE.ID_SPOT AND SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND ID_GEL2D=$itemID AND VALID_STATUS=?";
	$filterQuery="SELECT COUNT(*) FROM PROTEIN_VALIDATION,ANALYSIS,SAMPLE,SPOT WHERE PROTEIN_VALIDATION.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SPOT.ID_SPOT=SAMPLE.ID_SPOT AND SEL_STATUS=-3 AND ID_GEL2D=$itemID";
	#($numMaldi)=$dbh->selectrow_array("SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS WHERE ID_SAMPLE=$itemID AND INSTRUMENT='MALDI-TOF-TOF'");
	#$validMascotQuery="SELECT COUNT(ID_ANALYSIS) FROM ANALYSIS WHERE ID_SAMPLE=$itemID AND VALID_STATUS=>=1 AND FILE_FORMAT='MASCOT.DAT'";
}
elsif ($item eq 'SPOT') {
	$anaQuery="SELECT $anaFieldStrg FROM SAMPLE,ANALYSIS WHERE SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND ID_SPOT=$itemID";
	#$validQuery="SELECT COUNT(ID_ANALYSIS) FROM SAMPLE,ANALYSIS WHERE SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND ID_SPOT=$itemID AND VALID_STATUS=?";
	$filterQuery="SELECT COUNT(*) FROM PROTEIN_VALIDATION,ANALYSIS,SAMPLE WHERE PROTEIN_VALIDATION.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SEL_STATUS=-3 AND ID_SPOT=$itemID";
	#($numMaldi)=$dbh->selectrow_array("SELECT COUNT(ID_ANALYSIS) FROM SAMPLE,ANALYSIS WHERE SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND ID_SPOT=$itemID AND INSTRUMENT='MALDI-TOF-TOF'");
	#$validMascotQuery="SELECT COUNT(ID_ANALYSIS) FROM SAMPLE,ANALYSIS WHERE SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND ID_SPOT=$itemID AND VALID_STATUS>=1 AND FILE_FORMAT='MASCOT.DAT'";
}
my $sthAna=$dbh->prepare($anaQuery);
my $sthLQ=$dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND ID_ANALYSIS=? AND FOCUS='peptide'");
$sthAna->execute;
my ($numNotImport,$numUnValid,$numPartValid,$numValid,$numMaldi,$numValidMascot,$numValidMascotDAT,$numValidPhenyx,$numValidMIS,$numValidSILAC,$numValidITRAQ,$numValidTMT)=(0,0,0,0,0,0,0,0,0,0,0,0);
while (my ($anaID,$validStatus,$instrument,$msType,$fileFormat,$labeling)=$sthAna->fetchrow_array) {
	$instrument='Default' unless $instrument;
	if ($validStatus==-1) {
		$numNotImport++;
		$numUnValid++;
	}
	elsif ($validStatus==0) {
		$numUnValid++;
	}
	elsif ($validStatus==1) {
		$numUnValid++;
		$numPartValid++;
	}
	if ($validStatus>=1) {
		$numValid++;
		$numValidMIS++ if $msType eq 'MIS';
		#if ($labeling) { (PP 17/08/17)
		#	if ($labeling=~/SILAC/i) {$numValidSILAC++;}
		#	elsif ($labeling=~/iTRAQ/i) {$numValidITRAQ++;}
		#	elsif ($labeling=~/TMT/i) {$numValidTMT++;}
		#	else { # Labeling type is not always mentioned in labeling name
		#		$sthLQ->execute($anaID);
		#		my ($okSILAC,$okITRAQ,$okTMT);
		#		while (my ($qAnnot)=$sthLQ->fetchrow_array) {
		#			my ($labelType)=($qAnnot=~/^LABEL=(\w+)/);
		#			if ($labelType eq 'SILAC') {$okSILAC=1;}
		#			elsif ($labelType eq 'ITRAQ') {$okITRAQ=1;}
		#			elsif ($labelType eq 'TMT') {$okTMT=1;}
		#		}
		#		$numValidSILAC++ if $okSILAC;
		#		$numValidITRAQ++ if $okITRAQ;
		#		$numValidTMT++ if $okTMT;
		#	}
		#}
	}
	$numMaldi++ if $instrument eq 'MALDI-TOF-TOF';
	$numValidMascot++ if (($fileFormat eq 'MASCOT.DAT' || $fileFormat eq 'MASCOT.PDM') && $msType eq 'MIS');
	$numValidMascotDAT++ if ( $fileFormat eq 'MASCOT.DAT' );
	$numValidPhenyx++ if ($fileFormat eq 'PHENYX.XML' && $msType eq 'MIS');
}
$sthAna->finish;
$sthLQ->finish;

#my $sthValid=$dbh->prepare($validQuery);
###>Not-validated + not-imported Analyses
#$sthValid->execute(-1);
#my ($numNotImport)=$sthValid->fetchrow_array;
###>Not-validated Analyses
#$sthValid->execute(0);
#my ($numUnValid)=$sthValid->fetchrow_array;
#$numUnValid+=$numNotImport;
###>Partially validated Analyses
#$sthValid->execute(1);
#my ($numPartValid)=$sthValid->fetchrow_array;
###>Validated Analyses
#$sthValid->execute(2);
#my ($numValid)=$sthValid->fetchrow_array;
#$sthValid->finish;
##my $numAnalyses=$numUnValid+$numPartValid+$numValid;
#$numUnValid+=$numPartValid;
#$numValid+=$numPartValid;
#my $numValidMascot=0;
#if ($numValid && $projectFullAccess) {
#	my $sthValMas=$dbh->prepare($validMascotQuery);
#	$sthValMas->execute;
#	($numValidMascot)=$sthValMas->fetchrow_array;
#	$sthValMas->finish;
#}

##>Import
my $disabImport=($projectFullAccess)? '' : ' disabled';
my $disabDecoy=($numUnValid && $projectFullAccess)? '' : ' disabled';
##>Combine
#my $disabCombine=($numUnValidDat>=2 && $projectAccess !~ /power/)? '' : ' disabled';
##>Filter + auto-selection
my $disabFilterSel=($numUnValid)? '' : ' disabled';
my $disabRemFilter=' disabled';
if ($numUnValid) {
	my ($filter)=$dbh->selectrow_array($filterQuery);
	print $filter;
	$disabRemFilter='' if $filter;
}

$dbh->disconnect;

##>Auto-selection
my $disabAutoSel=($numUnValid && $projectFullAccess)? '' : ' disabled';
my $disabValidation=($numUnValid-$numNotImport && $projectFullAccess)? '' : ' disabled';
my $disabElution=($numMaldi && $projectFullAccess)? '' : ' disabled';
###>Report validation
#my $disabReport=($numUnValid-$numNotImport && ($projectAccess eq 'bioinfo' || $projectAccess eq 'mass'))? '' : ' disabled';
##>Delete
my $disabDelete=(($numUnValid || $numValid) && $deleteAccess)? '' : ' disabled';
my $disabDuplicate=($numUnValid)? '' : ' disabled';
my $disabMaxquant=($item eq 'EXPERIMENT')?'' : ' disabled';

##>Quantification
my $disabemPAI=($numValidMascot && $projectFullAccess)? '' : ' disabled';
my $disabXIC=($numValidMIS && $projectFullAccess)? '' : ' disabled';
my $disabSIN=($numValidMascotDAT && $projectFullAccess)? '' : ' disabled';
my $disabSILAC=($numValidSILAC)? '' : ' disabled';
my $disabITRAQ=($numValidITRAQ)? '' : ' disabled';
my $disabTMT=($numValidTMT)? '' : ' disabled';
my $disabWatchQuantif=' disabled';
my $disabWatchPRS=($projectFullAccess)? '' : ' disabled';
my $disabPkv=($item eq 'EXPERIMENT')?'' : ' disabled';
my $disabOpenSwath=($item eq 'EXPERIMENT')?'' : ' disabled';
my $disabOpenSwathImport=($item eq 'EXPERIMENT')?'' : ' disabled';
my $disabSpectronautImport=($item eq 'EXPERIMENT')?'' : ' disabled';
my $disabTDA=($item eq 'EXPERIMENT')?'' : ' disabled';
if (-e "$promsPath{tmp}/quantification/current") {
	opendir (DIR, "$promsPath{tmp}/quantification/current");
	while (my $file = readdir (DIR))	{
		if ($file=~/^\d+\.*\d*_\d+_/) {
			$disabWatchQuantif='';
			last;
		}
	}
	close DIR;
}


####>Checking if files in private batch directory<####
#my $batchString='';
#my $batchFilesDir="$promsPath{tmp}/batch/$userID";
#if (-e $batchFilesDir && ($projectAccess eq 'bioinfo' || $projectAccess eq 'mass')) { # && ($item eq 'EXPERIMENT' || $item eq 'SAMPLE')
#	opendir (DIR, $batchFilesDir);
#	while (defined (my $file = readdir (DIR)))	{
#		next if -d "$batchFilesDir/$file"; # directory
#		$batchString='*';
#		last;
#	}
#	close DIR;
#}

#######################
####>Starting HTML<####
#######################
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my ($popBgColor,$popBdColor)=&promsConfig::getPopupColors;
print header(-charset=>'UTF-8',-'content-encoding'=>'no');
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Manage Multiple Analyses</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
var selType='';
function displayProcesses(type) {
	if (selType) document.getElementById(selType).style.display='none';
	if (type) document.getElementById(type).style.display='block';
	selType=type;
}
function watchQuantifications() {
	var watchQuantifWin=window.open("$promsPath{cgi}/watchQuantifications.cgi",'WatchQuantifWindow','width=1200,height=500,scrollbars=yes,resizable=yes');
	watchQuantifWin.focus();
	parent.optionFrame.selectOption(parent.optionFrame.document.getElementById('summary')); // refresh optionFrame with summary option
}
function watchPhosphoRS() {
	var watchPhosphoWin=window.open("$promsPath{cgi}/watchPhosphoAnalyses.cgi",'WatchPhosphoWindow','width=1000,height=500,scrollbars=yes,resizable=yes');
	watchPhosphoWin.focus();
	parent.optionFrame.selectOption(parent.optionFrame.document.getElementById('summary')); // refresh optionFrame with summary option
}
function selectAction(action) {
	switch(action) {
		case 'import':
			window.location="./importBatchAnalyses.cgi?ID=$branchID&ACT=start";
			break;
		case 'decoy':
			window.location="./importBatchAnalyses.cgi?ID=$branchID&ACT=start&decoy=1";
			break;
		case 'importElution':
			window.location="./selectAnalyses.cgi?ID=$branchID&id_project=$projectID&callType=impElution";
			break;
		case 'duplicate':
			//if (confirm('Duplicate from all partially/non-validated Analyses found in selected item?')) {
			//	window.location="./filterAnalysis.cgi?ITEM=$item&ID=$itemID&id_project=$projectID&duplicate=1&apply=1";
			//}
			window.location="./selectAnalyses.cgi?ID=$branchID&id_project=$projectID&callType=duplicate";
			break;
		case 'applyFilter':
			//window.location="./filterAnalysis.cgi?ITEM=$item&ID=$itemID&id_project=$projectID";
			window.location="./selectAnalyses.cgi?ID=$branchID&id_project=$projectID&callType=appFilter";
			break;
		case 'removeFilter':
			//if (confirm('Remove filter from all unvalidated Analyses found in selected item?')) {
			//	window.location="./filterAnalysis.cgi?id_project=$projectID&ITEM=$item&ID=$itemID&remFilter=1";
			//}
			window.location="./selectAnalyses.cgi?ID=$branchID&id_project=$projectID&callType=remFilter";
			break;
		case 'lowScores':
			//window.location="./autoSelect.cgi?id_project=$projectID&ITEM=$item&ID=$itemID&ACT=lowScores";
			window.location="./selectAnalyses.cgi?ID=$branchID&id_project=$projectID&callType=lowScores";
			break;
		case 'qualiSelection':
			window.location="./autoSelect.cgi?id_project=$projectID&ITEM=$item&ID=$itemID&ACT=select";
			break;
		case 'compSelection':
			window.location="./filterPeptideAnalysis.cgi?id_project=$projectID&ITEM=$item&ID=$itemID";
			break;
		case 'removeFlags':
			window.location="./selectAnalyses.cgi?ID=$branchID&callType=remFlag";
			break;
		case 'clearSelection':
			//window.location="./autoSelect.cgi?id_project=$projectID&ITEM=$item&ID=$itemID&ACT=clearSelection";
			window.location="./selectAnalyses.cgi?ID=$branchID&callType=clearSel";
			break;
		case 'reportEnd':
			window.location="./selectAnalyses.cgi?ID=$branchID&callType=report";
			break;
		case 'delete':
			window.location="./selectAnalyses.cgi?ID=$branchID&callType=delete";
			break;
		case 'phosphoRS':
			window.location="./selectAnalyses.cgi?ID=$branchID&callType=phosphoRS";
			break;
		// Quantifications
		case 'xic':
			window.location="./selAna4Quantification.cgi?ID=$branchID&quantifType=xic";
			break;
		case 'xicmcq':
			window.location="./selAna4Quantification.cgi?ID=$branchID&quantifType=xicmcq";
			break;
		case 'empai':
			window.location="./selAna4Quantification.cgi?ID=$branchID&quantifType=empai";
			break;
		case 'sin':
			window.location="./selAna4Quantification.cgi?ID=$branchID&quantifType=sin";
			break;
/* Label-based internal quantification no longer supported (PP 17/08/17)
		case 'silac':
			window.location="./selAna4Quantification.cgi?ID=$branchID&quantifType=silac";
			break;
		case 'itraq':
			window.location="./selAna4Quantification.cgi?ID=$branchID&quantifType=itraq";
			break;
		case 'tmt':
			window.location="./selAna4Quantification.cgi?ID=$branchID&quantifType=tmt";
			break;
*/
		case 'mzxml' :
			window.location="./managemzXMLFiles.cgi?id_project=$projectID&ID=$branchID";
			break;
		case 'maxquant' :
			window.location="./importMaxquant.cgi?id_project=$projectID&ID=$itemID";
			break;
		case 'pkv' :
			window.location="./importSwathData.cgi?id_project=$projectID&ID=$branchID&ACT=import&FORMAT=pkv";
			break;
		case 'openswath' :
			window.location="./importSwathDataRefonte.cgi?id_project=$projectID&ID=$branchID&ACT=quantification&FORMAT=openswath";
			break;
		case 'openswathImport' :
			window.location="./importSwathData.cgi?id_project=$projectID&ID=$branchID&ACT=import&FORMAT=openswath&USERID=$userID";
			break;
		case 'spectronautImport' :
			window.location="./importSwathData.cgi?id_project=$projectID&ID=$branchID&ACT=import&FORMAT=spectronaut";
			break;
		case 'prm' :
			window.location="./importTDAData.cgi?id_project=$projectID&ID=$branchID&ACT=import&FORMAT=prm";
			break;
		default:
			alert('"'+action+'" is not recognized');
			break;
	}
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Process Multiple Analyses</FONT>
<BR><BR>
<FONT class="title">Process type :</FONT><SELECT name="type" class="title" onchange="displayProcesses(this.value)">
	<OPTION value="">-= Select =-</OPTION>
	<OPTION value="anaMag">Analysis management</OPTION>
	<OPTION value="pepSel">Analysis validation</OPTION>
	<OPTION value="anaQuant">Analysis quantification</OPTION>
</SELECT>
<BR><BR><BR>
<DIV id="anaMag" style="display:none">
<TABLE bgcolor=$darkColor cellpadding=4 width=500>
<TR><TD colspan=2><FONT class="title2">&nbsp;Analysis management:</TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('import')"$disabImport/></TH><TD bgcolor=$lightColor width=100%>&nbsp;<FONT class="title3">Import multiple analyses</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('decoy')"$disabDecoy/></TH><TD bgcolor=$lightColor width=100%>&nbsp;<FONT class="title3">Import decoy data into multiple analyses</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('importElution')"$disabElution/></TH><TD bgcolor=$lightColor width=100%>&nbsp;<FONT class="title3">Import elution time into multiple analyses</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('delete')"$disabDelete/></TH><TD bgcolor=$lightColor>&nbsp;<FONT class="title3">Delete multiple analyses</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('duplicate')"$disabDuplicate/></TH><TD bgcolor=$lightColor>&nbsp;<FONT class="title3">Duplicate multiple analyses</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('maxquant')"$disabMaxquant/></TH><TD bgcolor=$lightColor>&nbsp;<FONT class="title3">Import <A href="http://maxquant.org/" target="_blank">MaxQuant</A> quantitation</FONT></TD></TR>
</TABLE>
</DIV>

<DIV id="pepSel" style="display:none">
<TABLE bgcolor=$darkColor cellpadding=4 width=500>
<TR><TD colspan=2><FONT class="title2">&nbsp;Peptide/protein selection:</TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('qualiSelection')"$disabAutoSel/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Qualitative auto-selection&nbsp;</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('compSelection')"$disabAutoSel/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Comparative auto-selection&nbsp;</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('clearSelection')"$disabAutoSel/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Clear peptide/protein selections&nbsp;</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('lowScores')"$disabAutoSel/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Activate lower-scoring peptides&nbsp;</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('reportEnd')"$disabValidation/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Report and/or end validation of analyses&nbsp;</FONT></TD></TR>
</TABLE>
<BR><BR>
<TABLE bgcolor=$darkColor cellpadding=4 width=500>
<TR><TD colspan=2><FONT class="title2">&nbsp;Filters and flags:</TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('applyFilter')"$disabFilterSel/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Apply a protein filter&nbsp;</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('removeFilter')"$disabRemFilter/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Remove protein filters&nbsp;</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('removeFlags')"$disabAutoSel/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Remove peptide flags&nbsp;</FONT></TD></TR>
</TABLE>
<BR><BR>
<TABLE bgcolor=$darkColor cellpadding=4 width=500>
<TR><TD colspan=2><FONT class="title2">&nbsp;Phosphorylation sites:</TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('phosphoRS')"/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Run PhosphoRS analysis&nbsp;</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="watchPhosphoRS()"$disabWatchPRS/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Monitor PhosphoRS analyses&nbsp;</FONT></TD></TR>
</TABLE>
</DIV>

<DIV id="anaQuant" style="display:none">
<INPUT type="button" value="Monitor on-going quantifications" onclick="watchQuantifications()"$disabWatchQuantif/>
<BR><BR>
<TABLE bgcolor=$darkColor cellpadding=4 width=500>
<TR><TD colspan=2><FONT class="title2">&nbsp;<A href="http://maxquant.org/" target="_blank">MaxQuant</A>:</TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('maxquant')"$disabMaxquant/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Import MaxQuant quantitation</FONT></TD></TR>
</TABLE>
<BR><BR>
<TABLE bgcolor=$darkColor cellpadding=4 width=500>
<TR><TD colspan=2><FONT class="title2">&nbsp;Peptide Quantification:</TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('mzxml')"/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Manage mzXML files</FONT></TD></TR>
<!-- <TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('xic')"$disabXIC/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">XIC extraction&nbsp;</FONT></TD></TR> -->
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('xicmcq')"$disabXIC/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">XIC extraction with <A href="http://pappso.inra.fr/bioinfo/masschroq/index.php" target="_blank">MassChroQ</A>&nbsp;</FONT></TD></TR>
</TABLE>
<BR><BR>
<TABLE bgcolor=$darkColor cellpadding=4 width=500>
<TR><TD colspan=2><FONT class="title2">&nbsp;Protein Label-free Quantification:</TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('empai')"$disabemPAI/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Import emPAI data&nbsp;</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('sin')"$disabSIN/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">SI<SUB>N</SUB> quantification&nbsp;</FONT></TD></TR>
</TABLE>
<BR><BR>
<!-- label-based internal quantification no longer supported (PP 17/08/17)
<TABLE bgcolor=$darkColor cellpadding=4 width=500>
<TR><TD colspan=2><FONT class="title2">&nbsp;Protein Label Quantification:</TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('silac')"$disabSILAC/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">SILAC-based quantification&nbsp;</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('itraq')"$disabITRAQ/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">iTRAQ-based quantification&nbsp;</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('tmt')"$disabTMT/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">TMT-based quantification&nbsp;</FONT></TD></TR>
</TABLE>
<BR><BR>
-->
<TABLE bgcolor=$darkColor cellpadding=4 width=500>
<TR><TD colspan=2><FONT class="title2">&nbsp;DIA/TDA Quantification:</TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('prm')"$disabTDA/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Import TDA data&nbsp;</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('pkv')"$disabPkv/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Import PeakView data&nbsp;</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('openswath')"$disabOpenSwath/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">OpenSwath based quantification&nbsp;</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('openswathImport')"$disabOpenSwathImport/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Import OpenSwath data&nbsp;</FONT></TD></TR>
<TR><TH><INPUT type="button" value=" Proceed " onclick="selectAction('spectronautImport')"$disabSpectronautImport/></TH><TD bgcolor=$lightColor width=100% nowrap>&nbsp;<FONT class="title3">Import Spectronaut data&nbsp;</FONT></TD></TR>
</TABLE>
</DIV>

</CENTER>
</BODY>
</HTML>
|;

#<TR><TH>&nbsp;<INPUT type="button" value=" Proceed " onclick="selectAction('combine')"$disabCombine/>&nbsp;</TH><TD bgcolor=$lightColor>&nbsp;<FONT class="title3">Combine analyses</FONT></TD></TR>

####>Revision history<####
# 1.4.4 Add userID for importSwathData.cgi call (GA 23/11/18)
# 1.4.3 Add "Monitor PhosphoRS" button (PP 08/11/18)
# 1.4.2 Change PRM Script path (VS 08/11/18)
# 1.4.1 Replace "PeakView based quantification" by "Import PeakView data" (MLP 09/02/18)
# 1.4.0 Add spectronaut import (MLP 12/01/18)
# 1.3.8	Add OpenSwath quantification (MLP 06/12/17)
# 1.3.7 Add OpenSwath import (MLP 12/09/17)
# 1.3.6 Label-based internal quantification no longer supported (PP 17/08/17)
# 1.3.5 Add TDA import (MLP 19/07/17)
# 1.3.4 Minor modification for TMT<BR>TODO: Internal quantifications -SILAC,iTRAQ,TMT- do not work (GA 03/04/17)
# 1.3.3 Minor modification for SIN that only works for MASCOT.DAT files (GA 03/02/17)
# 1.3.2 Added MaxQuant in Analysis quantification section (PP 10/01/17)
# 1.3.1 Minor modification on Swath quanti (MLP 13/07/16)
# 1.3.0 Add Swath import option through Analysis quantification  (MLP 09/05/16)
# 1.2.9 Add Maxquant import option through Analysis management (GA 24/09/15)
# 1.2.8 Restrict item "super deletion" to bioinfo/mass/manag users only (PP 05/12/14)
# 1.2.7 Set undefined $instrument to Default (PP 08/10/13)
# 1.2.6 Options display reorganization (PP 11/09/13)
# 1.2.5 Rename action 'mcq' to 'xicmcq' (PP 09/09/13)
# 1.2.4 Turn off old XIC extraction in Label Free Quantitation (GA 03/07/13)
# 1.2.3 Better detection of labeling method (PP 19/06/13)
# 1.2.2 Add MassChroQ extraction selection (GA 04/09/12)
# 1.2.1 Add PhosphoRS action selection (FY 01/03/12)
# 1.2.0 Merge 1.1.9GA & 1.1.7PP (02/04/12)
# 1.1.9 Add a new option in Analysis Quantification: import of mzXML files (GA 03/01/12)
# 1.1.8 Initiate the value $labeling to prevent from warnings (GA 24/11/11)
# 1.1.7 Resizable & scrollable quantification watch window (PP 22/11/11)
# 1.1.6 Manager has full access validation data (PP 01/10/11)
# 1.1.5 Correction for detection of labeled quantifications (PP 27/09/11)
# 1.1.4 Adds Label/label-free internal quantification options (PP 26/08/11)
# 1.1.3 Quantification option disabled (PP 22/04/11)
# 1.1.2 Merge of 1.1.1PP & 1.1.1FY (PP 11/04/11)
# 1.1.1 Adds super user/admin access to analysis management (15/03/11)<BR>Redirect to "selectAnalyses.cgi" for clear selection and clear flags options (FY 22/03/11)
# 1.1.0 Protein filtering moved to Pep/Prot selection (PP 25/02/11)
# 1.0.9 New display, quantification options (PP 11/2010)
