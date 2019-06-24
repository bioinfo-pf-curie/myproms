#!/usr/local/bin/perl -w

################################################################################
# manageReferencePeptides.cgi       1.0.0                                      #
# Authors: P. Poullet (Institut Curie)                                         #
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
use promsConfig;
use promsMod;
use POSIX qw(strftime); # to get the time
use XML::Simple;
use strict;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $MAX_NUM_PEPTIDES=50;

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

#############################
my $action= (param('ACT'))? param('ACT') : 'list' ;

if (param('save')) { # edit/add was submitted
	
	###>Fetching all parameters<###
	my $setRefID = &promsMod::cleanNumericalParameters(param('ID'));
	my $name = param('name');
	my $des = param('description');
	my $colType = param('colType');
	my $protIdent = param('protIdent');
	my $protSeq = param('protSeq');
	my $numPep=param('numPep');
	my @pepID=param('pepID');
	my @validPep=param('validPep');
	my @pepSeq=param('pepSeq');
	my @charge=param('charge');
	my @monoMass=param('monoMass');
	my @pepRT=param('pepRT');
	my $pepData="<PEPTIDE_DATA version=\"1.0.1\">\n";
	foreach my $idx (0..$numPep-1) {
		my $exclStrg=($validPep[$idx])? '' : ' excluded="1"';
		$charge[$idx]=~s/\D+//g; $charge[$idx]=2 unless $charge[$idx];
		$pepData.="<PEPTIDE pepId=\"$pepID[$idx]\" sequence=\"$pepSeq[$idx]\" monoMass=\"$monoMass[$idx]\" charge=\"$charge[$idx]\" iRT=\"$pepRT[$idx]\"$exclStrg/>\n";
	}
	$pepData.="</PEPTIDE_DATA>";
	
	###>Updating DB<###
	if ($setRefID) { # edit
		my $sthUp=$dbh->prepare("UPDATE REFERENCE_RT SET NAME=?,DES=?,COLUMN_TYPE=?,PROT_IDENTIFIER=?,PROT_SEQ=?,NUM_PEP=?,DATA=? WHERE ID_REFERENCE_RT=?") || die $dbh->errstr;
		$sthUp->execute($name,$des,$colType,$protIdent,$protSeq,$numPep,$pepData,$setRefID);
		$sthUp->finish;
	}
	else { # add
		my $sthIns=$dbh->prepare("INSERT INTO REFERENCE_RT (ID_REFERENCE_RT,NAME,DES,COLUMN_TYPE,PROT_IDENTIFIER,PROT_SEQ,NUM_PEP,DATA) VALUES (?,?,?,?,?,?,?,?)");
		my ($maxRefID)=$dbh->selectrow_array("SELECT MAX(ID_REFERENCE_RT) FROM REFERENCE_RT");
		$sthIns->execute(++$maxRefID,$name,$des,$colType,$protIdent,$protSeq,$numPep,$pepData);
		$sthIns->finish;
	}
	
	$dbh->commit;
	
	$action='list';
}
elsif ($action eq 'delete') {
	my $setRefID = &promsMod::cleanNumericalParameters(param('ID'));
	$dbh->do("DELETE FROM REFERENCE_RT WHERE ID_REFERENCE_RT=$setRefID") || die $dbh->errstr;
	$dbh->commit;
	$action='list';
}
if ($action eq 'list') {
	
	###>Fetching data<###
	my %referencePeptides;
	my $sthGetIRT = $dbh->prepare("SELECT ID_REFERENCE_RT,NAME,DES,COLUMN_TYPE,PROT_IDENTIFIER,NUM_PEP,DATA FROM REFERENCE_RT");
	my @sthUsed=(
		$dbh->prepare("SELECT 1 FROM ANALYSIS_REFRT WHERE ID_REFERENCE_RT=? LIMIT 1"),
		$dbh->prepare("SELECT 1 FROM SWATH_LIB WHERE ID_REFERENCE_RT=? LIMIT 1")
	);
	$sthGetIRT->execute;
	while (my ($setRefID,@pepData) = $sthGetIRT->fetchrow_array) {
		my $isUsed=0;
		foreach my $sth (@sthUsed) {
			$sth->execute($setRefID);
			($isUsed)=$sth->fetchrow_array;
			last if $isUsed;
		}
		push @pepData,$isUsed;
		@{$referencePeptides{$setRefID}}=@pepData;
	}
	$sthGetIRT->finish;
	foreach my $sth (@sthUsed) {$sth->finish;}
	
	$dbh->disconnect;

	###>HTML<###
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	my ($lightColor,$darkColor)=&promsConfig::getRowColors;
	my $bgColor=$lightColor;
	print qq
|<HTML><HEAD>
<TITLE>List of Reference Peptide Sets</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE>
UL.ulList {
	padding:0px 5px 0px 5px;
	margin: 0px 5px 0px 5px;
	list-style:none;
}
UL.ulList LI {margin:0; padding:0px 10px 0px 10px;}
</STYLE>
<SCRIPT type="text/javascript">
function updatePeptideDisplay(pepBut,pepID) {
	[document.getElementById('pepDIV_'+pepID).style.display,pepBut.value]=(pepBut.value=='more')? ['inline-block','less'] : ['none','more'];
}
function deleteSet(id,name){
	if (confirm('Delete '+name+' ?')) {
		window.location = "./manageReferencePeptides.cgi?ACT=delete&ID="+id;
	}
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
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
<INPUT type="button" class="title2" value=" Add a new Reference Set " onclick="window.location='./manageReferencePeptides.cgi?ACT=add'">
<BR><BR>
<TABLE border=0 cellspacing=0 cellpadding=4>
|;
	foreach my $setRefID (sort{$a<=>$b} keys %referencePeptides) {
		my ($setName,$des,$colType,$protIdent,$numPep,$pepData,$isUsed) = @{$referencePeptides{$setRefID}};
		$setName=~s/'/ /g;
		if ($des){
			$des = &promsMod::HTMLcompatible($des);
			$des =~ s/<BR>/<BR>&nbsp;&nbsp;&nbsp; /g;
			$des =~ s/&apos;*/\'/g; ### for IE
		}
		else {$des = 'No description';}
		my $disabDel=($isUsed)? 'disabled' : '';
		print qq
|<TR bgcolor=$bgColor>
<TD></TD>
<TD><FONT class="title2">$setName</FONT><BR>
<UL class="ulList">
<LI><B>Description: </B>$des</LI>
<LI><B>Column type: </B>$colType</LI>
<LI><B>Protein: </B>$protIdent</LI>
<LI><B>Number of peptides: </B>$numPep &nbsp;<INPUT type="button" id="pepBUT_$setRefID" value="more" onclick="updatePeptideDisplay(this,$setRefID)"/></LI>
<LI><DIV id="pepDIV_$setRefID" style="display:none;background-color:#FFF">
<TABLE cellspacing=0 bgcolor=$darkColor style="margin:2px">
<TR bgcolor=$darkColor><TH class="rbBorder">#</TH><TH class="rbBorder">&nbsp;Sequence&nbsp;</TH><TH class="rbBorder">&nbsp;Charge&nbsp;</TH><TH class="rbBorder">&nbsp;Mass<SUB>(mono)</SUB>&nbsp;</TH><TH class="bBorder">&nbsp;Ret.&nbsp;time&nbsp;(%)&nbsp;</TH></TR>
|;
	my $xmlRT = new XML::Simple();
	my $xmlData = $xmlRT->XMLin($pepData);
	my @pepList;
	foreach my $pepInfo (@{$xmlData->{PEPTIDE}}) {
		my $pepImg=($pepInfo->{excluded})? 'bad.gif' : 'good.gif';
		push @pepList,[$pepImg,$pepInfo->{sequence},$pepInfo->{charge},$pepInfo->{monoMass},$pepInfo->{iRT}];
	}
	my $pepBgColor=$lightColor;
	my $count=0;
	foreach my $refPepInfo (sort{$a->[-1]<=>$b->[-1]} @pepList) {
		$count++;
		my ($pepImg,$pepSeq,$charge,$monoMass,$pepRT)=@{$refPepInfo};
		print "<TR bgcolor=\"$pepBgColor\"><TD align=\"right\">$count&nbsp;</TD><TD><IMG src=\"$promsPath{images}/$pepImg\">$pepSeq</TD><TD align=\"center\">$charge</TD><TD align=\"center\">$monoMass</TD><TD align=\"center\">$pepRT</TD></TR>\n";
		$pepBgColor=($pepBgColor eq $lightColor)? $darkColor : $lightColor;
	}
	print qq
|</TABLE>
</DIV></LI>
</UL>
</TD>
<TH width=100>
<INPUT type="button" value="Edit" style="width:75px" onclick="window.location='./manageReferencePeptides.cgi?ACT=edit&ID=$setRefID'"><BR>
<INPUT type="button" value="Delete" style="width:75px" onclick="deleteSet($setRefID,'$setName')" $disabDel>
</TH>
</TR>
|;
		$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
	}
	if (scalar keys %referencePeptides == 0) {
		print "<TR><TH colspan=3><BR><FONT class=\"title2\">(No reference peptide sets recorded)</FONT></TH></TR>\n";
	}
	print qq
|</TABLE>
</CENTER>
<BR><BR><BR>
</BODY>
</HTML>
|;
}
elsif ($action=~/add|edit/) {
	
	my $setRefID = ($action eq 'add')? 0 : &promsMod::cleanNumericalParameters(param('ID'));
	
	###>Fetching data<###
	my ($setName,$des,$colType,$protIdent,$protSeq,$numPep,$pepData)=('','','','','',2,'');
	my @pepList;
	my $isUsed=0;
	if ($action eq 'edit') {
		($setName,$des,$colType,$protIdent,$protSeq,$numPep,$pepData) = $dbh->selectrow_array("SELECT NAME,DES,COLUMN_TYPE,PROT_IDENTIFIER,PROT_SEQ,NUM_PEP,DATA FROM REFERENCE_RT WHERE ID_REFERENCE_RT=$setRefID");
		unless($des){ $des = '' };
		my $xmlRT = new XML::Simple();
		my $xmlData = $xmlRT->XMLin($pepData);
		foreach my $pepInfo (@{$xmlData->{PEPTIDE}}) {
			push @pepList,[$pepInfo->{pepId},$pepInfo->{excluded},$pepInfo->{sequence},$pepInfo->{charge},$pepInfo->{monoMass},$pepInfo->{iRT}];	
		}
		$numPep=scalar @pepList unless $numPep;
		my @sthUsed=(
			$dbh->prepare("SELECT 1 FROM ANALYSIS_REFRT WHERE ID_REFERENCE_RT=? LIMIT 1"),
			$dbh->prepare("SELECT 1 FROM SWATH_LIB WHERE ID_REFERENCE_RT=? LIMIT 1")
		);
		foreach my $sth (@sthUsed) {
			$sth->execute($setRefID);
			($isUsed)=$sth->fetchrow_array;
			last if $isUsed;
		}
		foreach my $sth (@sthUsed) {$sth->finish;}
	}
	$dbh->disconnect;
	
	my $usedPepStrg=($isUsed)? '<FONT style="font-weight:bold;color:#DD0000">This Reference Set has been used. Peptide data cannot be fully modified.</FONT><BR>'
					: "&nbsp;<B>Number of peptides:</B><INPUT type=\"number\" name=\"numPep\" value=\"$numPep\" min=\"2\" max=\"$MAX_NUM_PEPTIDES\" style=\"width:45px\" onclick=\"updatePeptideDisplay(this.value)\"/><BR>";

	### HTML ###
	my $titleStrg=($action eq 'add')? 'Recording a new Reference Peptide Set' : 'Editing Reference Peptides Set';
	my ($lightColor,$darkColor)=&promsConfig::getRowColors;
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>$titleStrg</TITLE>
<LINK rel=stylesheet href="$promsPath{html}/promsStyle.css" type=text/css>
<STYLE type="text/css">
</STYLE>
<SCRIPT type="text/javascript">
function updatePeptideDisplay(numPep) {
	for (let i=1; i<=$MAX_NUM_PEPTIDES; i++) {
		document.getElementById('TR_'+i).style.display=(i<=numPep)? '' : 'none';
	}
}
function checkForm(myForm) {
	if (!myForm.name.value) {
		alert('Enter a name for this Set.');
		return false;
	}
	if (!myForm.colType.value) {
		alert('Enter a column type for this Set.');
		return false;
	}
	if (!myForm.protIdent.value) {
		alert('Enter a protein identifier for this Set.');
		return false;
	}
	if (!myForm.protSeq.value) {
		alert('Enter a protein sequence for this Set.');
		return false;
	}
	if (myForm.numPep.value <= 1) {
		alert('At least 2 peptide sequences are required.');
		return false;
	}
	for (let i=1; i<=myForm.numPep.value; i++) {
		var Idx=i-1;
		if (!myForm.pepSeq[Idx].value \|\| !myForm.charge[Idx].value \|\| !myForm.monoMass[Idx].value \|\| !myForm.pepRT[Idx].value)  {
			alert('Data are missing for peptide #'+i+'.');
			return false;
		}
		myForm.validPep[Idx].value=(myForm.validPepChk[Idx].checked)? 1 : 0;
	}
	return true;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR><BR>
<FONT class="title">$titleStrg</FONT><BR><BR>

<FORM name="setForm" method="post" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="ID" value="$setRefID">
<TABLE bgcolor=$darkColor>
<TR><TH align=right valign=top width=150>&nbsp;Name :</TH><TD nowrap bgcolor=$lightColor><INPUT type="text" name="name" style="width:250px" value="$setName"/></TD></TR>
<TR><TH align=right valign=top>&nbsp;Description :</TH><TD nowrap bgcolor=$lightColor><TEXTAREA name="description" rows="2" cols="65">$des</TEXTAREA></TD></TR>
<TR><TH align=right>&nbsp;Column type :</TH><TD nowrap bgcolor=$lightColor><INPUT type="text" name="colType" style="width:100px" value="$colType"/></TD></TR>
<TR><TH align=right>&nbsp;Protein identifier :</TH><TD nowrap bgcolor=$lightColor><INPUT type="text" name="protIdent" style="width:250px" value="$protIdent"/></TD></TR>
<TR><TH align=right valign=top>&nbsp;Protein Sequence :</TH><TD nowrap bgcolor=$lightColor><TEXTAREA name="protSeq" rows="4" cols="65">$protSeq</TEXTAREA></TD></TR>
<TR><TH align=right valign=top>&nbsp;Peptides :</TH><TD nowrap bgcolor=$lightColor>$usedPepStrg<DIV style="display:inline-block;background-color:#FFF">
<TABLE cellspacing=0 bgcolor=$darkColor style="margin:2px">
<TR bgcolor=$darkColor><TH class="rbBorder">#</TH><TH class="rbBorder">&nbsp;Valid&nbsp;</TH><TH class="rbBorder">&nbsp;Sequence&nbsp;</TH><TH class="rbBorder">&nbsp;Charge&nbsp;</TH><TH class="rbBorder">&nbsp;Mass<SUB>(mono)</SUB>&nbsp;</TH><TH class="bBorder">&nbsp;Ret.&nbsp;time&nbsp;(%)&nbsp;</TH></TR>
|;
	my $pepBgColor=$lightColor;
	my $curPepID=0;
	foreach my $i (1..$MAX_NUM_PEPTIDES) {
		last if ($isUsed && $i > @pepList);
		$curPepID++;
		my ($pepID,$isBad,$pepSeq,$charge,$monoMass,$pepRT)=($pepList[$i-1])? @{$pepList[$i-1]} : ($curPepID,'','','','','');
		my $checkStrg=($isBad)? '' : 'checked';
		my $pepStrg=($isUsed)? "$pepSeq<INPUT type=\"hidden\" name=\"pepSeq\" value=\"$pepSeq\"/>" : "<INPUT type=\"text\" name=\"pepSeq\" value=\"$pepSeq\" style=\"width:200px\"/>";
		my $rtStg=($isUsed)? "$pepRT<INPUT type=\"hidden\" name=\"pepRT\" value=\"$pepRT\"/>" : "<INPUT type=\"text\" name=\"pepRT\" value=\"$pepRT\" style=\"width:70px\"/>";
		my $trDispStrg=($i<=2 || $pepList[$i-1])? '' : 'style="display:none"';
		print qq
|<TR id="TR_$i" bgcolor=$pepBgColor $trDispStrg>
<TD align="right">$i&nbsp;<INPUT type="hidden" name="pepID" value="$pepID"/></TD>
<TD align="center"><INPUT type="checkbox" name="validPepChk" value=1 $checkStrg><INPUT type="hidden" name="validPep" value=1 $checkStrg></TH>
<TD>$pepStrg</TD>
<TD align="center"><INPUT type="text" name="charge" value="$charge" style="width:30px"/></TD>
<TD align="center"><INPUT type="text" name="monoMass" value="$monoMass" style="width:80px"/></TD>
<TD>&nbsp;$rtStg</TD>
</TR>
|;
		$curPepID=$pepID if $curPepID < $pepID; # to handle jumps is pepID increment
		$pepBgColor=($pepBgColor eq $lightColor)? $darkColor : $lightColor;
	}
	print qq
|</TABLE></DIV>
</TD></TR>
<TR><TH colspan=2><INPUT type="submit" name="save" value=" Save "/>&nbsp;<INPUT type="button" value=" Cancel " onclick="window.location='./manageReferencePeptides.cgi?ACT=list'"/></TH></TR>
</TABLE>
</FORM>
</CENTER>
</BODY>
</HTML>
|;
}

####>Revision history<####
# 1.0.0 New script to manage reference peptide sets (PP 30/04/19)
