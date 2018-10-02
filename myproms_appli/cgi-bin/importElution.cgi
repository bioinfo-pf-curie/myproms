#!/usr/local/bin/perl -w

################################################################################
# importElution.cgi    1.0.2                                                   #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Imports elution time for MALDI-TOF-TOF analyses.                             #
# Called from editProjectItem.cgi or selectAnalyses.cgi                        #
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
use promsConfig;
use promsMod;
use strict;

#print header; warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;

####################
####>Parameters<####
####################
my $projectID=param('id_project');
my @analysisList=param('anaList');

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

#########################
####>Processing Form<####
#########################
if (param('save')) {

	print header(-'content-encoding'=>'no',-charset=>'utf-8'); # start_html,"\n"; # start_html required to force update
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Import Elution Data</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<CENTER>
<BR>
<FONT class="title">Importing Elution Data</FONT>
</CENTER>
<BR>
<BR>
|;

	####>Extracting Elution time from XML file<####
	print "<FONT class=\"title2\">+Reading template file...";
	my %labelElutionTime;
	my $spotCounter=0;
	my $XML=upload('templateFile');
	while (<$XML>) {
		if (/retentionTime="([^"]+)/) {
			my $retTime=sprintf "%.2f",$1;
			$retTime*=1;
			my ($label)=($_=~/ label="([^"]+)/);
			$labelElutionTime{$label}=$retTime;
		}
		$spotCounter++;
		if ($spotCounter==200) {
			print '.';
			$spotCounter=0;
		}
	}
	close $XML;
	print " DONE.</FONT><BR>";

	####>Looping through list of analyses<####
	my $sthAna=$dbh->prepare("SELECT NAME,VALID_STATUS,DATA_FILE FROM ANALYSIS WHERE ID_ANALYSIS=?");
	my $sthQ=$dbh->prepare("SELECT ID_QUERY,QUERY_NUM FROM QUERY_VALIDATION WHERE ID_ANALYSIS=?");
	my $sthP=$dbh->prepare("SELECT ID_PEPTIDE,QUERY_NUM FROM PEPTIDE WHERE ID_ANALYSIS=?");
	my $sthUpQ=$dbh->prepare("UPDATE QUERY_VALIDATION SET ELUTION_TIME=? WHERE ID_QUERY=?");
	my $sthUpP=$dbh->prepare("UPDATE PEPTIDE SET ELUTION_TIME=? WHERE ID_PEPTIDE=?");

	foreach my $anaID (@analysisList) {
		$sthAna->execute($anaID);
		my ($anaName,$validStatus,$dataFile)=$sthAna->fetchrow_array;
		print "<FONT class=\"title2\"><BR>+Processing Analysis <FONT color=\"#DD0000\">$anaName</FONT>...";

		##>Scanning .dat files for spot Labels
		my $usedDatFile;
		if ($validStatus<=1) {
			$usedDatFile="$promsPath{valid}/ana_$anaID/$dataFile";
		}
		else {
			$usedDatFile="$promsPath{peptide}/proj_$projectID/ana_$anaID/$dataFile";
			$usedDatFile=~s/\.dat\Z/_min\.dat/ unless -e $usedDatFile;  # assumes Fxxxxxx_min.dat file
		}
		unless (-e $usedDatFile) {
			print "<FONT color=\"#DD0000\">File $dataFile found!</FONT></FONT>\n";
			next;
		}
		open (DAT,$usedDatFile);
		my %queryLabels;
		my $curQuery;
		while (<DAT>) {
			if (/ name="query(\d+)"/) {
				$curQuery=$1;
				next;
			}
			if ($curQuery && /title=Label%3a%20(\w+)%2c/) {
				$queryLabels{$curQuery}=$1;
			}
		}
		close DAT;
		unless (scalar keys %queryLabels) {
			print "<FONT color=\"#DD0000\">No matching labels found!</FONT></FONT>\n";
			next;
		}
		print '.';

		##>Updating Elution time
		if ($validStatus<=1) {
			$sthQ->execute($anaID);
			while (my ($queryID,$queryNum)=$sthQ->fetchrow_array) {
				if ($queryLabels{$queryNum} && $labelElutionTime{$queryLabels{$queryNum}}) {
					$sthUpQ->execute($labelElutionTime{$queryLabels{$queryNum}},$queryID);
				}
			}
		}
		if ($validStatus>=1) {
			$sthP->execute($anaID);
			while (my ($peptideID,$queryNum)=$sthP->fetchrow_array) {
				if ($queryLabels{$queryNum} && $labelElutionTime{$queryLabels{$queryNum}}) {
					$sthUpP->execute($labelElutionTime{$queryLabels{$queryNum}},$peptideID);
				}
			}
		}
		print " DONE.</FONT><BR>\n";

		$dbh->commit;
	}

	$sthAna->finish;
	$sthQ->finish;
	$sthP->finish;
	$sthUpQ->finish;
	$sthUpP->finish;

	$dbh->disconnect;

	sleep 2;
	print qq
|<SCRIPT LANGUAGE="Javascript">
top.promsFrame.selectedAction='summary';
parent.optionFrame.selectOption();
</SCRIPT>
</HEAD></HTML>
|;
	exit;
}


#####################
####>Import Form<####
#####################
my ($analysisName)=$dbh->selectrow_array("SELECT NAME FROM ANALYSIS WHERE ID_ANALYSIS=$analysisList[0]");
$dbh->disconnect;

####>Starting HTML<####
my ($color1,$color2)=&promsConfig::getRowColors;
print header(-'content-encoding'=>'no',-charset=>'utf-8');
warningsToBrowser(1);
print qq
|<HEAD>
<TITLE>Import Elution Data</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
function checkForm(myForm) {
	if (!myForm.templateFile.value) {
		alert('ERROR: No template file provided!');
		return false;
	}
	return true;
}
</SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<CENTER>
<BR>
<FONT class="title">Import Elution Data for MALDI Analysis <FONT color="#DD0000">$analysisName</FONT></FONT>
<BR>
<BR>
<FORM name="importForm" method="post" onsubmit="return(checkForm(this));" enctype="multipart/form-data">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="anaList" value="$analysisList[0]">
<TABLE border=0 width=500><TR><TD bgcolor=$color2>
<TABLE border=0 width=100% cellpadding=2>
<TR>
	<TH align=right valign=top bgcolor=$color2 nowrap>Plate template file :</TH>
	<TD bgcolor=$color1><INPUT type="file" name="templateFile" size=80/><BR><SMALL><B><I>(XML format)</I></B></SMALL></TD>
</TR>
<TR>
	<TD colspan=2 align=center>
	<INPUT type="submit" name="save" value=" Save ">
	&nbsp;&nbsp;&nbsp;<INPUT type="button" value=" Cancel " onclick="parent.optionFrame.selectOption()">
	</TD>
</TR>
</TABLE>
</TD></TR></TABLE>
</FORM>
</BODY>
</HTML>
|;

####>Revision history<####
# 1.0.2 GPL license (PP 23/09/13)
# 1.0.1 Uses new data file path for validated .dat analyses (PP 06/03/13)
