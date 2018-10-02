#!/usr/local/bin/perl -w
################################################################################
# manageMotifAnalysis.cgi       1.0.1                                          #
# Authors: P. Poullet, S.Liva (Institut Curie)                                 #
# Contact: myproms@curie.fr                                                    #
# summary, edit, delete a motif enrichment analysis                            #
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
use promsQuantif;
use promsMod;
use POSIX qw(strftime); # to get the time
use File::Path qw(rmtree);

#print header(-charset=>'utf-8');warningsToBrowser(1);#DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $dbh=&promsConfig::dbConnect;
my $action=param('ACT') || 'summary';
my $motifID=param('ID');##can be motifID or heatMapMotifID
my $anaType=param('TYPE') || '';##MOTIF or HM_MOTIF
my $projectID=&promsMod::getProjectID($dbh,$motifID,'EXPLORANA');

my  ($catID, $expID, $name, $desc, $paramList, $filterList, $catExlusion, $status, $upDate, $user)=$dbh->selectrow_array("SELECT  ID_CATEGORY,ID_EXPERIMENT, NAME, DES, PARAM_LIST, FILTER_LIST, CAT_EXCLUSION, STATUS, UPDATE_DATE, UPDATE_USER from EXPLORANALYSIS where ID_EXPLORANALYSIS=$motifID");
$desc='' unless $desc;
$desc=&promsMod::HTMLcompatible($desc);

my ($pathToFile, %residues);
if ($anaType eq "MOTIF") {
  $pathToFile="$promsPath{explorAna}/project_".$projectID."/$motifID";
  %residues=("S"=>"Serine", "T"=>"Threonine", "Y"=>"Tyrosine", "C"=>"Cysteine", "H"=>"Histidine", "K"=>"Lysine", "R"=>"Arginine", "D"=>"Aspartic acid");
}

if ($action eq "delete") {
  print header(-'content-encoding'=>'no',-charset=>'utf-8');
  warningsToBrowser(1); # DEBUG

  if ($anaType eq "MOTIF") {
	if (!$catID) {
	  my $sthDeleteExplorQuantif=$dbh->do("DELETE from EXPLORANA_QUANTIF where ID_EXPLORANALYSIS=$motifID");
	}
	my $sthDeleteMotif=$dbh->do("DELETE from EXPLORANALYSIS where ID_EXPLORANALYSIS=$motifID");
	rmtree($pathToFile);
  }
  else {
	my $sthDeleteHMmotif=$dbh->do("DELETE from PARENT_EXPLORANALYSIS where ID_EXPLORANALYSIS=$motifID");
	my $sthDeleteMotif=$dbh->do("DELETE from EXPLORANALYSIS where ID_EXPLORANALYSIS=$motifID");
  }
  $dbh->commit;
  $dbh->disconnect;
  print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
    parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=experiment:$expID&branchID=motifanalysis:$motifID&VIEW=functAna";
</SCRIPT>
</HEAD>
</HTML>
|;
  exit;
}

if (param('submit')) {##Edit Action

  my $description = param('description')? param('description') : "";
  my $motifName = param('motifName');
  my $sthUpdateMotifAna=$dbh->do("UPDATE EXPLORANALYSIS set NAME = '$motifName', DES = '$description' where ID_EXPLORANALYSIS = $motifID") or die "Cannot prepare: " . $dbh->errstr();
  $dbh -> commit;
  $dbh -> disconnect;

  ####>Updating all frames<###
  print header(-'content-encoding'=>'no',-charset=>'utf-8');
  warningsToBrowser(1); # DEBUG
  print qq
|<HTML>
<HEAD>
<TITLE>Update All Frames</TITLE>
<SCRIPT LANGUAGE="JavaScript">
	top.promsFrame.selectedAction = 'summary';
	parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=experiment:$expID&branchID=motifanalysis:$motifID&VIEW=functAna";
</SCRIPT>
</HEAD>
</HTML>
|;
  exit;
}

my $statusStrg;
if ($anaType eq 'MOTIF') {
  $statusStrg= ($status == -2)? "<IMG src=$promsPath{images}/bad.gif>&nbsp;&nbsp;<FONT color=\"red\">***ERROR***</FONT>" : ($status == -1)? "<FONT color=\"orange\"><B>0n-going analysis</B></FONT>" : "<IMG src=$promsPath{images}/good.gif>&nbsp;<FONT color=\"green\">Finished</FONT>";
}

#### START HTML
print header(-charset=>'utf-8');
warningsToBrowser(1);
my $strgMotif = ($anaType eq "MOTIF")? "Motif Enrichment Analysis : <FONT color=\"red\">$name</FONT>" : "Heatmap Motif Enrichment Analysis : <FONT color=\"red\">$name</FONT>";
my $title = ($action eq 'summary')? $strgMotif :  "Editing $strgMotif" ;
print qq
|<HTML>
<HEAD>
<TITLE>Motif Analysis</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT langage="Javascript">
|;
if ($action eq 'edit') {
    print qq
|function checkForm(myForm) {
	if (!myForm.motifName.value) {
		alert('A name is expected for the analysis');
		return false;
	}
	return true;
}
|;
}
else { # summary
  if ($anaType eq "MOTIF" && $status == -1) {
	print qq
|//setTimeout(function(){window.location.reload(true);},5000);// Reload page every 5 seconds
setTimeout(function(){parent.optionFrame.location.reload(1);},5000);// Reload page every 5 seconds
|;
  }
}
print qq
|</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">$title</FONT><BR><BR>
|;
if ($action eq 'edit') {
  print qq
|<FORM NAME="displayMenu" method="POST" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="ID" value="$motifID" />
|;
}
my ($light,$dark)=&promsConfig::getRowColors;
print qq
|<TABLE border="0" width="500" cellspacing="2" cellpadding="2" bgcolor="$dark">
|;
if ($action eq 'edit') {
  print qq
|<TR><TH align="right" nowrap>Analysis name :</TH><TD bgcolor="$light"><INPUT type="text" name="motifName" size="50" maxlength="100" value="$name"/></TD></TR>
<TR><TH align="right" valign="top" nowrap>&nbsp;Description :</TH><TD bgcolor="$light"><TEXTAREA name="description" rows="2" cols="50">$desc</TEXTAREA></TD></TR>
<TR><TH colspan="2"><INPUT type="submit" name="submit" value="Update"></TD></TR>
|;
}
else {#summary
  if ($anaType eq "MOTIF") {
	#$desc=&promsMod::HTMLcompatible($desc);
	my ($strgDataSource,$strgCentralRes, $strgOccurence, $strgWidth, $strgBackground, $strgSignificance) = split("//",$paramList);
	my ($itemCentRes, $centralRes)=split("=",$strgCentralRes);
	my ($itemOcc, $occurence)=split("=",$strgOccurence);
	my ($itemWidth, $width)=split("=",$strgWidth);
	my ($itemBack, $backgroundInfo)=split("=",$strgBackground);
	my ($backgroundStrgValue,$backgroundValue)=split(":",$backgroundInfo);
	my ($itemSign,$significance)=split("=",$strgSignificance);
	my $strgTable = ($catID)? "List" : "Quantification";
	my ($foregroundName,$backgroundName);
	if ($catID) {
	  my ($catName, $className)=$dbh->selectrow_array("SELECT CA.NAME, CL.NAME FROM CATEGORY CA, CLASSIFICATION CL WHERE CA.ID_CATEGORY=$catID and CA.ID_CLASSIFICATION=CL.ID_CLASSIFICATION");
	  $foregroundName="&nbsp;$className > $catName";
	  if ($backgroundStrgValue eq "quanti") {
		my %groupQuantifNames=&promsQuantif::getDistinctQuantifNames($dbh,$backgroundValue);
		$backgroundName="&nbsp;Quantification&nbsp;:&nbsp;<b>$groupQuantifNames{'FULL'}{$backgroundValue}</b><br>";
	  }
	  else {
		$backgroundName="&bull;&nbsp;Type&nbsp;:&nbsp;<b>$backgroundStrgValue</b><br>&bull;&nbsp;Nb. seqs&nbsp;:&nbsp;<b>$backgroundValue</b>";
	  }
	}
	else {
	  my $quantifID=$dbh->selectrow_array("SELECT CONCAT(ID_QUANTIFICATION,\"_\",TARGET_POS) FROM EXPLORANA_QUANTIF WHERE ID_EXPLORANALYSIS=$motifID");
	  my %groupQuantifNames=&promsQuantif::getDistinctQuantifNames($dbh,$quantifID);
	  $foregroundName=$groupQuantifNames{'FULL'}{$quantifID};
	  if ($backgroundStrgValue eq "quanti") {
		$backgroundName="&nbsp;Type&nbsp;:&nbsp;<b>selected quantification</b><br>";
	  }
	  else {
		$backgroundName="&bull;&nbsp;Type&nbsp;:&nbsp;<b>$backgroundStrgValue</b><br>&bull;&nbsp;Nb. seqs&nbsp;:&nbsp;<b>$backgroundValue</b>";
	  }
	}

	my @logoMotif;
	my $strgMotifFound="";
	if ( -d $pathToFile ) {
	  opendir(my $dh, $pathToFile) || die "Can't opendir $pathToFile: $!";
	  @logoMotif = grep { /^logoMotif/  } readdir($dh);
	  closedir $dh;
	  $strgMotifFound= (scalar(@logoMotif) > 0)? "&nbsp;<b>".scalar(@logoMotif)." motif(s) found" : "<b>No motif found</b>";
	}

print qq
|<TR><TH align="right" width="150" nowrap>Analysis name :</TH><TD bgcolor="$light" nowrap>&nbsp;<b>$name</b></TD></TR>
<TR><TH align="right" width="150" nowrap>$strgTable :</TH><TD bgcolor="$light" nowrap><b>$foregroundName</b></TD></TR>
<TR><TH align="right" valign="top" nowrap>&nbsp;Description :</TH><TD bgcolor="$light">&nbsp;$desc</TD></TR>
<TR><TH align="right" valign="top" nowrap>&nbsp;Foreground selection :</TH>
  <TD nowrap bgcolor=$light>
	&bull;&nbsp;Central Residue&nbsp;:&nbsp;<b>$residues{$centralRes}</b><br>
	&bull;&nbsp;Min.seqs&nbsp;:&nbsp;<b>$occurence</b><br>
	&bull;&nbsp;Width&nbsp;:&nbsp;<b>$width</b><br>
	&bull;&nbsp;Significance&nbsp;:&nbsp;<b>$significance</b><br>
  </TD></TR>
  <TR><TH align="right" valign="top" nowrap>&nbsp;Background selection :</TH>
  <TD nowrap bgcolor=$light>$backgroundName</TD></TR>
<TR><TH align="right" valign="top" nowrap>&nbsp;Feature selection :</TH><TD nowrap bgcolor=$light>
|;
	if (!$catID) {
	  my ($strgPVAL, $strgFC, $strgTypeFC, $strgInfRatio) = split("//", $filterList);
	  my $FCthreshold = (split("=",$strgFC))[1];
	  my $PVALthreshold = (split("=",$strgPVAL))[1];
	  my $FCtype = (split("=",$strgTypeFC))[1];
	  my $infRatio=(split("=",$strgInfRatio))[1];
	  my $strgFoldChangeType = ($FCtype eq "abs")? "Up &ge;&nbsp;$FCthreshold or Down &le; 1/$FCthreshold" : ($FCtype eq "up")? "Up &ge; $FCthreshold" : "Down &le; 1/$FCthreshold";
	  my $strgRatio = ($infRatio eq "Y")? "infinite ratios excluded" : "infinite ratios not excluded";
		print qq|&bull;&nbsp;<b>$strgRatio</b><br>
				&bull;&nbsp;FC&nbsp;:&nbsp;<b>$strgFoldChangeType</b><br>
				&bull;&nbsp;p-value&nbsp;:&nbsp;<b>$PVALthreshold</b>
	  |;
	}
	else {
	  print qq|&nbsp;No feature available, list of proteins is involved|;
	}
	print qq|</TD></TR>
	<TR><TH align="right" bgcolor="$dark" nowrap>Motif :</TH><TD nowrap bgcolor="$light">$strgMotifFound</TD></TR>
	<TR><TH align="right" bgcolor="$dark" nowrap>Status :</TH><TD nowrap bgcolor="$light">&nbsp;&nbsp;$statusStrg</TD></TR>
	|;
  }
  else {
	my $sthSelParentExplo=$dbh->prepare("SELECT ID_PARENT_EXPLORANALYSIS, DISPLAY_POS FROM PARENT_EXPLORANALYSIS where ID_EXPLORANALYSIS=? ORDER BY DISPLAY_POS");
	$sthSelParentExplo->execute($motifID);
	my @fullName;
	while (my ($parentExplo, $displayPos)=$sthSelParentExplo->fetchrow_array) {
		my ($idCat, $explorName)=$dbh->selectrow_array("SELECT ID_CATEGORY, NAME FROM EXPLORANALYSIS WHERE ID_EXPLORANALYSIS=$parentExplo");
		push @fullName, "&bull;&nbsp;$explorName";
	}
	print qq|
	<TR><TH align="right" width="150" nowrap>Analysis name :</TH><TD bgcolor="$light" nowrap>$name</TD></TR>
	<TR><TH align="right" width="150" nowrap valign="top">Selection :</TH><TD bgcolor="$light" nowrap>|;
	print join("<br>",@fullName);
	print qq|</TD></TR>
	<TR><TH align="right" valign="top" nowrap>&nbsp;Description :</TH><TD bgcolor="$light">&nbsp;$desc</TD></TR>
	|;
  }
}

print "</TABLE><BR>\n";
if ($action eq 'edit'){
    print "</FORM>\n";
}

print qq|<DIV id="listProt" style="float:left;display:none"></DIV>
</CENTER>
</BODY>
</HTML>
|;
$dbh -> disconnect;

####>Revision history<####
# 1.0.1 add list option
# 1.0.0 new script to manage, edit, delete motif enrichment analysis (SL 21/10/14)
