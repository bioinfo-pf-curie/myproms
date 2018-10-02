#!/usr/local/bin/perl -w
################################################################################
# displayMotifAnalysis.cgi       1.0.1                                         #
# Authors: P. Poullet, S.Liva (Institut Curie)      	                       #
# Contact: myproms@curie.fr                                                    #
# display the results of Motif analysis        	                       		   #
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
use promsQuantif;
use POSIX qw(strftime); # to get the time
my %promsPath=&promsConfig::getServerInfo;

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;
my ($lightColor,$darkColor)=&promsConfig::getRowColors;

my $userID=$ENV{'REMOTE_USER'};
my $projectID=param('PROJECT_ID');
my $motifID=param('ID');
my $ajax = (param('AJAX'))? param('AJAX') : '';
my $pathToFile = "$promsPath{explorAna}/project_$projectID/$motifID";
my $imgMotif = "$promsPath{html}/data/exploratory_data/project_$projectID/$motifID";
my ($name)=$dbh->selectrow_array("SELECT NAME FROM EXPLORANALYSIS WHERE ID_EXPLORANALYSIS=$motifID");

if ($ajax eq 'ajaxListProt') {
  &ajaxGetProteins;
  $dbh->disconnect;
  exit;
}

my (%resultMotif, %order);
open(MOTIF,"$pathToFile/motifResult.txt");
while(my $line=<MOTIF>) {
  next if ($.==1);
  $line=~s/"//g;
  my ($index, $motif, $score, $fgMatches, $fgSize, $bgMatches, $bgSize, $foldChange)=split(/\t/,$line);
  $resultMotif{$motif}="<b>".sprintf("%.3f",$foldChange)."</b>-".sprintf("%.3f",$score)."-".$fgMatches."-".$fgSize."-".$bgMatches."-".$bgSize;
  $order{$motif}=$bgSize;
}
close(MOTIF);

#### START HTML
print header(-'content-encoding'=>'no',-charset=>'utf-8');
warningsToBrowser(1);
print qq|
<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">|;
&promsMod::popupInfo();
print qq|
function sequenceView(listAnaID, protID, idType, ms) {
  var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+listAnaID+"&id_prot="+protID+"&msdata="+ms+"&id_type="+idType;
  top.openProtWindow(winLocation);
}

// AJAX --->
var XHR=null;
function ajaxGetProteins (motif){
  var protTable=document.getElementById("protTable");
  protTable.style.display='block';
  protTable.innerHTML='<IMG src="$promsPath{images}/scrollbarGreen.gif">';

  //If XHR object already exists, the request is canceled & the object is deleted
  if(XHR && XHR.readyState != 0){
	XHR.abort();
	delete XHR;
  }
  var paramStrg="AJAX=ajaxListProt&ID=$motifID&PROJECT_ID=$projectID&MOTIF="+motif;
  //Creation of the XMLHTTPRequest object
  XHR = getXMLHTTP();
  if (!XHR) {
	return false;
  }
  XHR.open("POST","$promsPath{cgi}/displayMotifAnalysis.cgi",true);
  //Send the proper header information along with the request
  XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
  XHR.setRequestHeader("Content-length", paramStrg.length);
  XHR.setRequestHeader("Connection", "close");
  XHR.onreadystatechange=function() {
	if (XHR.readyState==4 && XHR.responseText) {
	  protTable.innerHTML=XHR.responseText;
	  protTable.scrollIntoView();
	}
  }
  XHR.send(paramStrg);
}

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
// <--- AJAX
</SCRIPT>
</HEAD>
<BODY>
<CENTER>
|;
if (!keys %resultMotif) {
  print qq|<br><br><b><font color="red">NO MOTIF FOUND</font></b>|;
}
else {
  print qq|
<FONT class="title">Display Motif Enrichment Analysis for <FONT color="red">$name</FONT></FONT><BR><BR>
<TABLE cellspacing=0 cellpadding=0 border=0>
<TR class="header" bgcolor=$darkColor><TH>&nbsp;Motif&nbsp;</TH><TH>&nbsp;Fold Change&nbsp;</TH><TH>&nbsp;Score&nbsp;</TH><TH>&nbsp;Fg. Matches&nbsp;</TH><TH>&nbsp;Fg. Size&nbsp;</TH><TH>&nbsp;Bg. Matches&nbsp;</TH><TH>&nbsp;Bg.Size&nbsp;</TH></TR>
|;
  my $bgColor=$lightColor;
  foreach my $motif (sort {$order{$b} <=> $order{$a}} keys %resultMotif) {
	print "<TR bgcolor=$bgColor nowrap><TD align=\"left\" ><A href=\"javascript:ajaxGetProteins('$motif')\"><IMG SRC=\"$imgMotif/logoMotif_".$motif."_.png\" width=120 heigth=auto border=\"1\"></A></TD>";
	foreach my $motifElem (split(/-/,$resultMotif{$motif})) {
	  print qq|<TD align="center" nowrap >$motifElem</TD>|;
	}
	print "</TR>";
	$bgColor= ($bgColor eq $lightColor)? $darkColor : $lightColor;
  }
  print "</TABLE><br><br>";
}
print qq|<DIV id="protTable" style="float:left;display:none"></DIV>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup();</SCRIPT>
</CENTER>
</BODY>
</HTML>|;

sub ajaxGetProteins {

  print header;warningsToBrowser(1);
  my $motif=param('MOTIF');

  ##for colorize the AA involvd i motif
  my %colorMotif;
  push my @seqMotif, split(//,$motif);
  for (my $i=0; $i<=$#seqMotif; $i++) {
	next if $seqMotif[$i] eq ".";
	$colorMotif{$motif}{$seqMotif[$i]}{$i}=$i;
  }

  my $sthGetSequence=$dbh->prepare("SELECT ID_MASTER_PROTEIN, ALIAS, PROT_DES, MW, ORGANISM FROM PROTEIN WHERE ID_PROTEIN=? and ID_PROJECT=$projectID");
  my ($geneNameID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
  my $sthMP=$dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN=? AND ID_IDENTIFIER=$geneNameID  ORDER BY RANK");

  my (%protInfo, %orderProt, %anaProt, %gene, %tmpGene);
  %gene=();
  open(PROT,"$pathToFile/foreground_protein.txt");
	while(my $line=<PROT>) {
	  chomp($line);
	  my ($protID)=(split(/\t/,$line))[0];
	  $sthGetSequence->execute($protID);
	  my ($masterProtID, $alias,$protDes, $mw, $organism)=$sthGetSequence->fetchrow_array;
	  $protDes="-" if (!$protDes);
	  $organism="-" if (!$organism);
	  $mw="-" if (!$mw);

	  $sthMP->execute($masterProtID);
	  while (my ($gene)=$sthMP->fetchrow_array) {
		next if $tmpGene{$protID}{$gene};
		push @{$gene{$protID}},$gene;
		$tmpGene{$protID}{$gene}=1;
	  }

	  my $strgGene = ($gene{$protID} && scalar @{$gene{$protID}} >= 1)? join(',',@{$gene{$protID}}) : "-";
	  if ($line=~/~/) {
		my ($idProt, $ambSite, $site, $sequence)=split(/\t/,$line);
		$protInfo{$sequence}=$protID."::".$alias."-".$ambSite."::".$site."::".$strgGene."::".$mw."::".$protDes."::".$organism;
		$orderProt{$sequence}=$alias."-".$ambSite;
		$anaProt{$idProt}=1;
	  }
	  else {
		my ($idProt, $site, $sequence)=split(/\t/,$line);
		$protInfo{$sequence}=$idProt."::".$alias."-".$site."::NA::".$strgGene."::".$mw."::".$protDes."::".$organism;
		$orderProt{$sequence}=$alias."-".$site;
		$anaProt{$idProt}=1;
	  }
	}
  close(PROT);

  my %sequenceMotif;
  open(MOTIF,"$pathToFile/sequenceMotif_".$motif."_.txt");
  while (my $line=<MOTIF>) {
	chomp($line);
	$line=~s/"//g;
	$sequenceMotif{$line}=1;
  }
  close(MOTIF);

  my (%proteinAnalysis);
  open(ANA,"$pathToFile/analysisList.txt");
  my $strgLine;
  while(my $line=<ANA>) {
	chomp($line);
	$strgLine=$line;
  }
  close(ANA);

  my $sthSelAnaProt=$dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=? and ID_ANALYSIS in ($strgLine) LIMIT 0,1");
  foreach my $protID (keys %anaProt) {
	next if $proteinAnalysis{$protID};
	$sthSelAnaProt->execute($protID);
	my ($anaID)=$sthSelAnaProt->fetchrow_array;
	if ($anaID) {
	  $proteinAnalysis{$protID}=$anaID;
	  next;
	}
  }
  my $massValue='MW<FONT style="font-size:9px;"> kDa</FONT>';

  print qq|<TABLE cellspacing=0 cellpadding=0 border=0 width=50%>|;
  print "<TR ><TH colspan=7><IMG SRC=\"$imgMotif/logoMotif_".$motif."_.png\" border=\"1\" width=600 heigth=auto></TH></TR>";
  print qq|<!--<TR ><TH >MOTIF $motif</TH></TR>-->
	  <TR class="header" bgcolor="$darkColor" ><TH nowrap>&nbsp;Isoforms&nbsp;</TH><TH nowrap>&nbsp;Sequence&nbsp;</TH><TH nowrap>&nbsp;Gene&nbsp;</TH><TH nowrap>&nbsp;$massValue&nbsp;</TH><TH nowrap>&nbsp;Description - Species&nbsp;</TH></TR>|;
  my $bgColor=$lightColor;
  foreach my $sequence (sort {$orderProt{$a} cmp $orderProt{$b} } keys %protInfo) {
	if ($sequenceMotif{$sequence}) {
	  print qq|<TR class="list" bgcolor="$bgColor"> |;
	  my ($idProt, $protein, $site, $gene, $mw, $desc,$organism)=split(/::/,$protInfo{$sequence});
	  my $shortenDesc = $desc;
	  my $shortenGene = $gene;
	  $mw=sprintf "%.1f",$mw/1000;

	  my @tmpGene = split(/,/,$gene);
	  my $geneName1=$tmpGene[0];
	  my $strgGene="";
	  if (scalar(@tmpGene)==1){
		$strgGene="<b>$geneName1</b>";
	  }
	  else {
		shift(@tmpGene);
		my $geneList=join('<br>',@tmpGene);
		$strgGene="<A href=\"javascript:void(null)\" onmouseover=\"popup('$geneList')\" onmouseout=\"popout()\"><u><b>$geneName1</b></u></A>";
	  }

	  my @seq=split(//,$sequence);
	  my $length=length($sequence);
	  my $index=($length-1)/2;
	  my $strgSequence="";
	  for (my $i=0; $i<=$#seq;$i++) {
		if ($i == $index) {
		  $strgSequence.="<font color=red>$seq[$i]</font>";
		}
		else {
		  if ($colorMotif{$motif}{$seq[$i]}{$i}) {
			$strgSequence.="<font color=orange>$seq[$i]</font>";
		  }
		  else {
			$strgSequence.=$seq[$i];
		  }
		}
	  }

	  my ($prot, $siteProt)=split(/-/,$protein);

	  my ($popupIso,$displayIso) = ($site eq "NA")? ($prot, $protein) : ($prot, "$protein ($site)");
	  print qq|<TD nowrap><A href="javascript:sequenceView('$proteinAnalysis{$idProt}','$idProt','valid',0)" onmouseover="popup('$popupIso')" onmouseout="popout()"><b>$displayIso</b></A>&nbsp;&nbsp;</TD>
		  <TD nowrap class="seq" align="center">&nbsp;$strgSequence&nbsp;</TD>
		  <TD nowrap align="left">&nbsp;&nbsp;$strgGene&nbsp;&nbsp;</TD>
		  <TD nowrap>&nbsp;&nbsp;$mw&nbsp;&nbsp;</TD>
		  <TD nowrap>&nbsp;&nbsp;<!--<A href="javascript:void(null)" onmouseover="popup('$desc')" onmouseout="popout()">-->$desc <!--</A>-->&nbsp;&nbsp;<i>$organism</i></TD>
		  <TR>
		  |;
	  $bgColor= ($bgColor eq $lightColor)? $darkColor : $lightColor;
	}
  }
	print qq|<TABLE>|;
}

####>Revision history<####
# 1.0.1 available for list (SL 26/07/17)
# 1.0.0 New script for displaying motif enrichment analysis (SL 27/06/17)