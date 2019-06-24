#!/usr/local/bin/perl -w

################################################################################
# manageLabelReagents.cgi       1.0.0                                          #
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
use File::Path qw(remove_tree);
use File::Copy qw(move);
use strict;

#######################
####>Configuration<####
#######################
my $MAX_NUM_TAGS=25;
my %promsPath=&promsConfig::getServerInfo;
my $reagentDir="$promsPath{data}/label_reagents";
my $reagentDirHtml="$promsPath{html}/data/label_reagents";

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

#############################
my $action= (param('ACT'))? param('ACT') : 'list';

if (param('save')) { # edit/add was submitted
	
	###>Fetching all parameters<###
	my $productID = &promsMod::cleanNumericalParameters(param('ID'));
	my $name = param('name');
	my $company = param('company');
	my $prodNum = param('prodNum');
	my $lotNum = param('lotNum');
	my $useSatus=param('useStatus');
	
	my ($isUsed)=($productID)? $dbh->selectrow_array("SELECT 1 FROM QUANTIFICATION WHERE ID_PRODUCT=$productID LIMIT 1") : 0;

	###>Updating DB<###
	if ($isUsed) { # edit
		my $sthUp1=$dbh->prepare("UPDATE ISOTOPIC_CORRECTION SET NAME=?,COMPANY=?,PRODUCT_NUMBER=?,LOT_NUMBER=?,USE_STATUS=? WHERE ID_PRODUCT=?") || die $dbh->errstr;
		$sthUp1->execute($name,$company,$prodNum,$lotNum,$useSatus,$productID);
		$sthUp1->finish;
	}
	else { # add/edit
		my $label = param('label');
		my $numTags=param('numTags');
		my @tagName=param('tagName');
		my @minus2=param('minus2');
		my @minus1=param('minus1');
		my @plus1=param('plus1');
		my @plus2=param('plus2');
		my $tagData="<MASS_TAGS>\n";
		foreach my $tagPos (1..$numTags) { # <MASS_TAG RANK="1" NAME="TMT-126" MINUS_2="0" MINUS_1="0" PLUS_1="8.2" PLUS_2="0.4"/>
			my $idx=$tagPos-1;
			foreach my $refField (\@minus2,\@minus1,\@plus1,\@plus2) {
				$refField->[$idx]=0 unless $refField->[$idx];
				$refField->[$idx]=~s/,/\./; # just to be safe
			}
			$tagData.="<MASS_TAG RANK=\"$tagPos\" NAME=\"$tagName[$idx]\" MINUS_2=\"$minus2[$idx]\" MINUS_1=\"$minus1[$idx]\" PLUS_1=\"$plus1[$idx]\" PLUS_2=\"$plus2[$idx]\"/>\n";
		}
		$tagData.="</MASS_TAGS>";
		if ($productID) {
			my $sthUp2=$dbh->prepare("UPDATE ISOTOPIC_CORRECTION SET NAME=?,LABEL_TYPE=?,COMPANY=?,PRODUCT_NUMBER=?,LOT_NUMBER=?,ISOTOPIC_DISTRIBUTION=?,USE_STATUS=? WHERE ID_PRODUCT=?");
			$sthUp2->execute($name,$label.':'.$numTags,$company,$prodNum,$lotNum,$tagData,$useSatus,$productID);
			$sthUp2->finish;
		}
		else {
			my $sthIns=$dbh->prepare("INSERT INTO ISOTOPIC_CORRECTION (ID_PRODUCT,NAME,LABEL_TYPE,COMPANY,PRODUCT_NUMBER,LOT_NUMBER,ISOTOPIC_DISTRIBUTION,USE_STATUS,RECORD_DATE,RECORD_USER) VALUES (?,?,?,?,?,?,?,?,NOW(),?)");
			($productID)=$dbh->selectrow_array("SELECT MAX(ID_PRODUCT) FROM ISOTOPIC_CORRECTION");
			$sthIns->execute(++$productID,$name,$label.':'.$numTags,$company,$prodNum,$lotNum,$tagData,$useSatus,$ENV{REMOTE_USER});
			$sthIns->finish;
		}
	}
	##>Product sheet
	remove_tree("$reagentDir/reag_$productID") if param('deleteFile');
	my $prodFile='';
	if (param('prodSheet')) {
		mkdir $reagentDir unless -e $reagentDir;
		remove_tree("$reagentDir/reag_$productID") if -e "$reagentDir/reag_$productID"; # simpler than finding the file name
		mkdir "$reagentDir/reag_$productID";
		my $prodFile=(split(/[\\\/]/,param('prodSheet')))[-1];
		my $tmpfile = tmpFileName(upload('prodSheet')); # name of temp file being uploaded
		move($tmpfile,"$reagentDir/reag_$productID/$prodFile");
	}
	
	$dbh->commit;
	
	$action='list';
}
elsif ($action eq 'delete') {
	my $productID = &promsMod::cleanNumericalParameters(param('ID'));
	$dbh->do("DELETE FROM ISOTOPIC_CORRECTION WHERE ID_PRODUCT=$productID") || die $dbh->errstr;
	$dbh->commit;
	remove_tree("$reagentDir/reag_$productID") if -e "$reagentDir/reag_$productID";
	$action='list';
}
if ($action eq 'list') {
	
	###>Fetching data<###
	my %productList;
	my $sthGetProd = $dbh->prepare("SELECT ID_PRODUCT,U.ID_USER,LABEL_TYPE,NAME,COMPANY,PRODUCT_NUMBER,LOT_NUMBER,ISOTOPIC_DISTRIBUTION,USE_STATUS,RECORD_DATE,USER_NAME
								   FROM ISOTOPIC_CORRECTION I LEFT JOIN USER_LIST U ON I.RECORD_USER=U.ID_USER");
	my $sthUsed=$dbh->prepare("SELECT 1 FROM QUANTIFICATION WHERE ID_PRODUCT=? LIMIT 1");
	$sthGetProd->execute;
	while (my ($prodID,$recUserID,@prodData) = $sthGetProd->fetchrow_array) {
		unless ($prodData[-1]) { # record user name
			$prodData[-1]=$recUserID || 'Unknown';
		}
		@{$productList{$prodID}}=@prodData;
		my $prodFile='';
		if (-e $reagentDir && -e "$reagentDir/reag_$prodID") {
			my @files=glob("$reagentDir/reag_$prodID/*");
			$prodFile=(split(/[\\\/]/,$files[0]))[-1] if $files[0];
		}
		$sthUsed->execute($prodID);
		my ($isUsed)=$sthUsed->fetchrow_array;
		$isUsed=0 unless $isUsed;
		push @{$productList{$prodID}},($prodFile,$isUsed);
	}
	$sthGetProd->finish;
	$sthUsed->finish;
	
	$dbh->disconnect;

	###>HTML<###
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	my ($lightColor,$darkColor)=&promsConfig::getRowColors;
	my $bgColor=$lightColor;
	print qq
|<HTML><HEAD>
<TITLE>List of Label Reagents</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE>
UL.ulList {
	padding:0px 5px 0px 5px;
	margin: 0px 5px 0px 5px;
	list-style:none;
}
UL.ulList LI {margin:0; padding:0px 10px 0px 10px;}
TH.tag {min-width:65px;}
</STYLE>
<SCRIPT type="text/javascript">
function updateTagsDisplay(tagBut,tagID) {
	[document.getElementById('tagDIV_'+tagID).style.display,tagBut.value]=(tagBut.value=='more')? ['inline-block','less'] : ['none','more'];
}
function displayProductSheet(prodID,prodFile) {
	window.open('$reagentDirHtml/reag_'+prodID+'/'+prodFile,'fileSheetWindow','width=1000,height=750');
}
function deleteProduct(id,name){
	if (confirm('Delete '+name+' ?')) {
		window.location = "./manageLabelReagents.cgi?ACT=delete&ID="+id;
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
<INPUT type="button" class="title2" value=" Add a new Label Reagent " onclick="window.location='./manageLabelReagents.cgi?ACT=add'">
<BR><BR>
<TABLE border=0 cellspacing=0 cellpadding=4 style="min-width:600px">
|;
	foreach my $prodID (sort{$a<=>$b} keys %productList) {
		my ($labelType,$prodName,$company,$prodNum,$lotNum,$distribXML,$useStatus,$recDate,$recUser,$prodFile,$isUsed) = @{$productList{$prodID}};
		$prodName=~s/'/ /g;
		my ($label,$plex)=split(':',$labelType);
		$company='Unknown' unless $company;
		$prodNum='Unknown' unless $prodNum;
		$lotNum='Unknown' unless $lotNum;
		my $useStatusStrg=($useStatus eq 'NO')? '&nbsp;<FONT color="#DD0000">[Not usable]</FONT>' : '';
		my $disabDel=($isUsed)? 'disabled' : '';
		my $xmlTag = new XML::Simple();
		my $xmlData = $xmlTag->XMLin($distribXML); # <MASS_TAG RANK="1" NAME="TMT-126" MINUS_2="0" MINUS_1="0" PLUS_1="8.2" PLUS_2="0.4"/>
		my @tagList;
		foreach my $tagInfo (@{$xmlData->{MASS_TAG}}) {
			push @tagList,[$tagInfo->{RANK},$tagInfo->{NAME},$tagInfo->{MINUS_2},$tagInfo->{MINUS_1},$tagInfo->{PLUS_1},$tagInfo->{PLUS_2}];
		}
		my $numTags=scalar @tagList;
		print qq
|<TR bgcolor=$bgColor>
<TD></TD>
<TD><FONT class="title2">$prodName$useStatusStrg</FONT><BR>
<UL class="ulList">
<LI><B>Label: </B>$label ${plex}plex</LI>
<LI><B>Company: </B>$company</LI>
<LI><B>Product number: </B>$prodNum&nbsp;&nbsp;<B>Lot number: </B>$lotNum</LI>
<LI><B>Record date: </B>$recDate by $recUser</LI>
<LI><B>Number of tags: </B>$numTags &nbsp;<INPUT type="button" id="tagBUT_$prodID" value="more" onclick="updateTagsDisplay(this,$prodID)"/></LI>
<LI><DIV id="tagDIV_$prodID" style="display:none;background-color:#FFF">
<TABLE cellspacing=0 bgcolor=$darkColor style="margin:2px">
<TR bgcolor=$darkColor><TH class="rbBorder">#</TH><TH class="rbBorder">&nbsp;Mass tag&nbsp;</TH><TH class="rbBorder tag">-2</TH><TH class="rbBorder tag">-1</TH><TH class="rbBorder tag">Mono</TH><TH class="rbBorder tag">+1</TH><TH class="bBorder tag">+2</TH></TR>
|;
	my $tagBgColor=$lightColor;
	my $count=0;
	foreach my $refTagInfo (sort{$a->[0]<=>$b->[0]} @tagList) {
		$count++;
		my ($rank,$tagName,$minus2,$minus1,$plus1,$plus2)=@{$refTagInfo};
		$minus2=sprintf '%.1f',$minus2;
		$minus1=sprintf '%.1f',$minus1;
		$plus1=sprintf '%.1f',$plus1;
		$plus2=sprintf '%.1f',$plus2;
		print "<TR bgcolor=\"$tagBgColor\"><TD align=\"right\">$count&nbsp;</TD><TH align=\"left\">&nbsp;$tagName&nbsp;</TH><TD align=\"center\">&nbsp;$minus2%&nbsp;</TD><TD align=\"center\">&nbsp;$minus1%&nbsp;</TD><TH>&nbsp;100%&nbsp;</TH><TD align=\"center\">&nbsp;$plus1%&nbsp;</TD><TD align=\"center\">&nbsp;$plus2%&nbsp;</TD></TR>\n";
		$tagBgColor=($tagBgColor eq $lightColor)? $darkColor : $lightColor;
	}
	my $prodSheetStrg=($prodFile)? "$prodFile <INPUT type=\"button\" value=\"Display\" onclick=\"displayProductSheet($prodID,'$prodFile')\"/>" : 'Not recorded';
	print qq
|</TABLE>
</DIV></LI>
<LI><B>Product sheet:</B> $prodSheetStrg</LI>
</UL>
</TD>
<TH width=100>
<INPUT type="button" value="Edit" style="width:75px" onclick="window.location='./manageLabelReagents.cgi?ACT=edit&ID=$prodID'"><BR>
<INPUT type="button" value="Delete" style="width:75px" onclick="deleteProduct($prodID,'$prodName')" $disabDel>
</TH>
</TR>
|;
		$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
	}
	if (scalar keys %productList == 0) {
		print "<TR><TH colspan=3><BR><FONT class=\"title2\">(No label reagent sets recorded)</FONT></TH></TR>\n";
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
	
	my $productID = ($action eq 'add')? 0 : &promsMod::cleanNumericalParameters(param('ID'));
	
	###>Fetching data<###
	my ($label,$prodName,$company,$prodNum,$lotNum,$distribXML,$useStatus,$numTags,$prodFile)=('','','','','','','YES',2,'');
	my @tagList;
	my $isUsed=0;
	if ($action eq 'edit') {
		($prodName,my $labelType,$company,$prodNum,$lotNum,$distribXML,$useStatus) = $dbh->selectrow_array("SELECT NAME,LABEL_TYPE,COMPANY,PRODUCT_NUMBER,LOT_NUMBER,ISOTOPIC_DISTRIBUTION,USE_STATUS FROM ISOTOPIC_CORRECTION WHERE ID_PRODUCT=$productID");
		($label,my $plex)=split(':',$labelType);
		$company='Unknown' unless $company;
		$prodNum='Unknown' unless $prodNum;
		$lotNum='Unknown' unless $lotNum;
		my $xmlTag = new XML::Simple();
		my $xmlData = $xmlTag->XMLin($distribXML);
		foreach my $tagInfo (@{$xmlData->{MASS_TAG}}) {
			push @tagList,[$tagInfo->{RANK},$tagInfo->{NAME},$tagInfo->{MINUS_2},$tagInfo->{MINUS_1},$tagInfo->{PLUS_1},$tagInfo->{PLUS_2}];	
		}
		$numTags=scalar @tagList;
		($isUsed)=$dbh->selectrow_array("SELECT 1 FROM QUANTIFICATION WHERE ID_PRODUCT=$productID LIMIT 1");
		if (-e $reagentDir && -e "$reagentDir/reag_$productID") {
			my @files=glob("$reagentDir/reag_$productID/*");
			$prodFile=(split(/[\\\/]/,$files[0]))[-1] if $files[0];
		}
	}
	$dbh->disconnect;

	###>HTML<###
	my $titleStrg=($action eq 'add')? 'Recording a new Label Reagent' : 'Editing Label Reagent';
	my ($lightColor,$darkColor)=&promsConfig::getRowColors;
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>$titleStrg</TITLE>
<LINK rel=stylesheet href="$promsPath{html}/promsStyle.css" type=text/css>
<STYLE type="text/css">
TH.tag {min-width:65px;}
</STYLE>
<SCRIPT type="text/javascript">
function updateTagsDisplay(numTags) {
	for (let i=1; i<=$MAX_NUM_TAGS; i++) {
		document.getElementById('TR_'+i).style.display=(i<=numTags)? '' : 'none';
	}
}
function checkForm(myForm) {
	if (!myForm.name.value) {
		alert('Enter a name for this reagent.');
		return false;
	}
	if (!myForm.label.value) {
		alert('Selet a label type for this reagent.');
		return false;
	}
|;
	unless ($isUsed) {
		print qq
|	if (myForm.numTags.value <= 1) {
		alert('At least 2 tags are required.');
		return false;
	}
	for (let i=1; i<=myForm.numTags.value; i++) {
		var Idx=i-1;
		if (!myForm.tagName[Idx].value \|\| (!myForm.minus2[Idx].value && !myForm.minus1[Idx].value && !myForm.plus1[Idx].value && !myForm.plus2[Idx].value)) {
			alert('Data are missing for tag #'+i+'.');
			return false;
		}
	}
|;
	}
	print qq
|	return true;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR><BR>
<FONT class="title">$titleStrg</FONT><BR><BR>

<FORM name="reagentForm" method="post" enctype="multipart/form-data" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="ID" value="$productID">
<TABLE bgcolor=$darkColor>
<TR><TH align=right valign=top width=150>&nbsp;Name<SUP>*</SUP> :</TH><TD nowrap bgcolor=$lightColor><INPUT type="text" name="name" style="width:350px" value="$prodName"/></TD></TR>
<TR><TH align=right valign=top>&nbsp;Label :</TH><TD nowrap bgcolor=$lightColor>
|;
	my $usedProdStrg='';
	if ($isUsed) {
		$label='iTRAQ' if $label eq 'ITRAQ';
		print $label.'&nbsp;'.$numTags.'plex';
		$usedProdStrg="<FONT style=\"font-weight:bold;color:#DD0000\">This label reagent has been used. Tag data cannot be modified.</FONT><BR>";
	}
	else {
		my ($selITRAQ,$selTMT)=($label eq 'ITRAQ')? ('select','') : ($label eq 'TMT')? ('','select') : ('','');
		print qq
|<SELECT name="label"><OPTION value="">-= ! =-</OPTION><OPTION value="ITRAQ" $selITRAQ>iTRAQ</OPTION><OPTION value="TMT" $selTMT>TMT</OPTION></SELECT>&nbsp;
<INPUT type="number" name="numTags" value="$numTags" min="2" max="$MAX_NUM_TAGS" style="width:45px" onclick="updateTagsDisplay(this.value)"/>plex
|;
	}
	print qq
|</TD></TR>
<TR><TH align=right>&nbsp;Company :</TH><TD nowrap bgcolor=$lightColor><INPUT type="text" name="company" style="width:350px" value="$company"/></TD></TR>
<TR><TH align=right>&nbsp;Product number :</TH><TD nowrap bgcolor=$lightColor><INPUT type="text" name="prodNum" style="width:100px" value="$prodNum"/></TD></TR>
<TR><TH align=right>&nbsp;Lot number :</TH><TD nowrap bgcolor=$lightColor><INPUT type="text" name="lotNum" style="width:100px" value="$lotNum"/></TD></TR>
<TR><TH align=right valign=top>&nbsp;Tags :</TH><TD nowrap bgcolor=$lightColor>$usedProdStrg<DIV style="display:inline-block;background-color:#FFF">
<TABLE cellspacing=0 bgcolor=$darkColor style="margin:2px">
<TR bgcolor=$darkColor><TH class="rbBorder">#</TH><TH class="rbBorder">&nbsp;Mass tag&nbsp;</TH><TH class="rbBorder tag">-2</TH><TH class="rbBorder tag">-1</TH><TH class="rbBorder tag">Mono</TH><TH class="rbBorder tag">+1</TH><TH class="bBorder tag">+2</TH></TR>
|;
	my $tagBgColor=$lightColor;
	foreach my $i (1..$MAX_NUM_TAGS) {
		last if ($isUsed && $i > @tagList);
		my ($rank,$tagName,$minus2,$minus1,$plus1,$plus2)=($tagList[$i-1])? @{$tagList[$i-1]} : ($i,'','','','','');
		my $trDispStrg=($i<=2 || $tagList[$i-1])? '' : 'style="display:none"';
		print qq
|<TR id="TR_$i" bgcolor=$tagBgColor $trDispStrg>
<TD align="right">$i&nbsp;</TD>
|;
		if ($isUsed) {
			$minus2=sprintf '%.1f',$minus2;
			$minus1=sprintf '%.1f',$minus1;
			$plus1=sprintf '%.1f',$plus1;
			$plus2=sprintf '%.1f',$plus2;
			print qq
|<TH align="left">&nbsp;$tagName&nbsp;</TH>
<TD align="center">$minus2%</TD>
<TD align="center">$minus1%</TD>
<TH>100%</TD>
<TD align="center">$plus1%</TD>
<TD align="center">$plus2%</TD>
|;
		}
		else {
			print qq
|<TH align="left">&nbsp;<INPUT type="text" name="tagName" value="$tagName" style="width:100px"/>&nbsp;</TH>
<TD align="center"><INPUT type="text" name="minus2" value="$minus2" style="width:30px"/>%</TD>
<TD align="center"><INPUT type="text" name="minus1" value="$minus1" style="width:30px"/>%</TD>
<TH>100%</TD>
<TD align="center"><INPUT type="text" name="plus1" value="$plus1" style="width:30px"/>%</TD>
<TD align="center"><INPUT type="text" name="plus2" value="$plus2" style="width:30px"/>%</TD>
|;
		}
		print "</TR>\n";
		$tagBgColor=($tagBgColor eq $lightColor)? $darkColor : $lightColor;
	}
	print "<TR bgcolor=\"$lightColor\"><TD colspan=7>&nbsp;<I>Empty fields will be set to 0 (Mass tag field cannot be empty)</I></TD></TR>\n" unless $isUsed;
	print qq
|</TABLE></DIV>
</TD></TR>
<TR><TH align=right valign=top>&nbsp;Product sheet :</TH><TD bgcolor=$lightColor>
|;
	if ($prodFile) {
		print qq
|&nbsp;$prodFile<INPUT type="hidden" name="prevFile" value="$prodFile"/><BR>
<INPUT type="checkbox" name="deleteFile" value="1"/><B>Delete<BR>or
Replace with:</B>
|;
	}
	my $selNoStrg=($useStatus eq 'YES')? '' : 'selected';
	print qq
|<INPUT type="file" name="prodSheet"></TD></TR>
<TR><TH align=right valign=top>&nbsp;Usable :</TH><TD nowrap bgcolor=$lightColor><SELECT name="useStatus"><OPTION value="YES">Yes</OPTION><OPTION value="NO" $selNoStrg>No</OPTION></SELECT></TD></TR>
<TR><TH colspan=2><INPUT type="submit" name="save" value=" Save "/>&nbsp;<INPUT type="button" value=" Cancel " onclick="window.location='./manageLabelReagents.cgi?ACT=list'"/></TH></TR>
</TABLE>
<SUP>*</SUP><I>Mandatory field</I>
</FORM>
</CENTER>
</BODY>
</HTML>
|;
}

####>Revision history<####
# 1.0.0 New script for managing Label reagents that are required to allow TMT/iTRAQ peptide XIC quantification correction (PP 02/05/19)
