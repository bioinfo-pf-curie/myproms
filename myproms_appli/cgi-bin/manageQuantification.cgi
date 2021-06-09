#!/usr/local/bin/perl -w

################################################################################
# manageQuantification.cgi               1.2.8                                 #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
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
use promsQuantif;
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
my $action = param('ACT') || 'edit'; # edit, update, delete, (summary of quantifs moved to show[Pep/Prot]Quantification.cgi)
my $itemID = param('ID') || param('itemID') ;
my $projectID = param('PROJECT_ID');
my $branchID = param('branchID');


##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);

################
####>Delete<####
################
if ($action =~ /delete/) {
	my @quantifIdList = (param('quantifIdList')) ? param('quantifIdList') : ($itemID);
	my $title = ($action eq 'delete') ? 'Deleting ' : 'Delete ';
	$title .= ($action eq 'deleteMulti' || scalar @quantifIdList > 1) ? 'Quantifications' : 'Quantification';
	
	print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>$title</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<style>
#quantiTable tr {
	height: 25px;
}

#quantiTable td, #quantiTable th {
	padding: 0 10px 0 10px;
}
</style>
<script>
function checkall(checkStatus){
	var quantiBox = document.selQuantiForm.quantifIdList;
	if (!quantiBox) return; // no selectable quantifications
	if (quantiBox.length) { // more than 1 checkboxes
		for (let i=0; i < quantiBox.length; i++){
			if (quantiBox[i].disabled == false) {
				quantiBox[i].checked=checkStatus;
			}
		}
	}
	else {quantiBox.checked=checkStatus;} // Only 1 checkbox
}
function checkForm(myForm) {
	var quantiBox = myForm.quantifIdList;
	var okChecked=false;
	if (quantiBox.length) { // more than 1 checkboxes
		for (let i=0; i < quantiBox.length; i++){
			if (quantiBox[i].disabled == false){
				okChecked=quantiBox[i].checked;
				if (okChecked) break;
			}
		}
	}
	else {okChecked=quantiBox.checked;} // Only 1 checkbox
	if (okChecked==false) {
		alert("ERROR: No quantification selected!");
	}
	return okChecked;
}
</script>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" >
<CENTER><FONT class="title">$title</FONT>
<BR><BR><BR>
|;
	
	if ($action eq 'delete') {
		print "</CENTER>\n";
		my ($experimentID)=$dbh->selectrow_array("SELECT DISTINCT(ID_EXPERIMENT) FROM SAMPLE S, ANALYSIS A, ANA_QUANTIFICATION AQ WHERE S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=$quantifIdList[0]");
		my $designID;
		my $sthQ = $dbh->prepare("SELECT ID_DESIGN,NAME,FOCUS FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
		
		foreach my $qID (@quantifIdList) {
			$sthQ->execute($qID);
			($designID,my $qName,my $focus)=$sthQ->fetchrow_array;
			print "<FONT class=\"title2\">-Deleting '$qName'</FONT>";
			my ($verbose,$isTargeted)=($focus eq 'peptide' && $branchID eq "quantification:$qID")? (2,1) : (1,0); # single XIC deeletion or anything else
			print "<BR>" if $verbose==2;
			&promsQuantif::deleteQuantification($dbh,$projectID,$qID,{VERBOSE=>$verbose,IS_TARGETED=>$isTargeted}); # TARGETED only applies to XIC extraction
			$dbh->commit;
			print "<FONT class=\"title2\">Done.<BR></FONT>\n";
		}
		$sthQ->finish;
		$dbh->disconnect;
		sleep 2;
		
		my $branchIDStg=($designID)? "&DESIGN=DESIGN:$designID&branchID=design:$designID" : "&branchID=experiment:$experimentID";
		print qq
|</FONT>
<SCRIPT type="text/javascript">
	top.promsFrame.selectedAction='summary';
	parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID$branchIDStg&ACT=experiment&EXPERIMENT=EXPERIMENT:$experimentID&ISNAVFRAME=0&VIEW=quanti";
</SCRIPT>
</BODY>
</HTML>
|;
		exit;	
	}
	elsif($action eq 'deleteMulti') {
		print qq
|<form name="selQuantiForm" action="./manageQuantification.cgi" method="post" onsubmit="return (checkForm(this));" style=''>
<input type='hidden' name='PROJECT_ID' value='$projectID' />
<input type='hidden' name='ACT' value='delete' />
<input type='hidden' name='branchID' value='$branchID' />

<table id='quantiTable' border=0 cellspacing=0 cellpadding=0>
|;

		my $visFilterStrg=($userInfo[1] eq 'bioinfo')? '' : ($userInfo[1] eq 'mass')? 'AND Q.STATUS <= 2' : 'AND Q.STATUS <= 1'; # 3: hidden to mass &  bio, 2: hidden to biologists
		my $sthQ = $dbh->prepare("SELECT Q.ID_QUANTIFICATION, Q.NAME, Q.STATUS, QM.NAME, QM.DES FROM QUANTIFICATION Q INNER JOIN QUANTIFICATION_METHOD QM ON Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD WHERE Q.ID_DESIGN=$itemID $visFilterStrg");
        my @sthChild=(
					  $dbh->prepare("SELECT 1 FROM GOANA_QUANTIFICATION WHERE ID_QUANTIFICATION=? LIMIT 1"),
					  $dbh->prepare("SELECT 1 FROM EXPLORANA_QUANTIF WHERE ID_QUANTIFICATION=? LIMIT 1"),
					  $dbh->prepare("SELECT 1 FROM PATHWAYANA_QUANTIFICATION WHERE ID_QUANTIFICATION=? LIMIT 1"),
					  $dbh->prepare("SELECT 1 FROM PARENT_QUANTIFICATION WHERE ID_PARENT_QUANTIFICATION=? LIMIT 1")
					 );
		$sthQ->execute;
		my ($lightColor, $darkColor)=&promsConfig::getRowColors;
		my $rowColor = $lightColor;
		print qq
|<tr bgcolor="$darkColor">
	<th class="rbBorder" nowrap align=left><input id="checkAllQuantis" type="checkbox" onclick="checkall(this.checked)" />&nbsp;Name</th>
	<th class="bBorder" nowrap>Method</th>
</tr>
|;
        while (my ($qID, $qName, $qStatus, $qmName, $qmDes) = $sthQ->fetchrow_array) {
			my $hasChildren=0;
			if ($qStatus < 1) {
				my $statusStrg=($qStatus == -2)? 'Failed' : ($qStatus == -1)? 'On-going' : 'Not completed';
				$qName = "$qName <span style='color:red;font-weight:bold'>($statusStrg !)</span>";
			}
			else {
				foreach my $sth (@sthChild) {
					$sth->execute($qID);
					($hasChildren)=$sth->fetchrow_array;
					last if $hasChildren;
				}
			}
			my ($visStrg,$disabStatus)=($hasChildren)? ('style="visibility:hidden"','disabled') : ('','');
			print qq
|<tr bgcolor="$rowColor">
	<th align="left"><input name="quantifIdList" value='$qID' type="checkbox" $visStrg $disabStatus/>&nbsp;$qName</th>
	<td style='text-align:left'>$qmName ($qmDes)</td>
</tr>
|;
			$rowColor = ($rowColor eq $darkColor) ? $lightColor : $darkColor;
		}
		print qq
|</table><br/>
<input type="submit" name="Submit" value="Delete" class="title3">
</form>
</center>
</BODY>
</HTML>
|;
		$sthQ->finish;
		foreach my $sth (@sthChild) {$sth->finish;}
		$dbh->disconnect;
		
		exit;
	}
}
#####################
####>Update Name<#### (following 'edit')
#####################
elsif ($action eq 'update') {
	my $name=param('quantiName');
	my $status=param('quantiStatus');
	my $sthUp=$dbh->prepare("UPDATE QUANTIFICATION SET NAME=?,STATUS=? WHERE ID_QUANTIFICATION=$itemID");
	$sthUp->execute($name,$status);
	$sthUp->finish;
	$dbh->commit;

	my ($experimentID)=$dbh->selectrow_array("SELECT DISTINCT(ID_EXPERIMENT) FROM SAMPLE S, ANALYSIS A, ANA_QUANTIFICATION AQ WHERE S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=$itemID");
	my ($item,$parentID,$childID)=split(':',$branchID);# For XIC extractions linked to a TnPQ or Pep-Ratio quantification
	$childID=$parentID unless $childID;
	my ($designID)=$dbh->selectrow_array("SELECT ID_DESIGN FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$childID");

	$dbh->disconnect;

	my $branchIDStg=($designID)? "&DESIGN=DESIGN:$designID&branchID=quantification:$childID" : "&branchID=experiment:$experimentID";
	print header(-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<SCRIPT type="text/javascript">
	top.promsFrame.selectedAction='summary';
	parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID$branchIDStg&ACT=experiment&EXPERIMENT=EXPERIMENT:$experimentID&ISNAVFRAME=0&VIEW=quanti";
</SCRIPT>
</HEAD>
</HTML>
|;
	exit;
	
}
else {
	###################
	####>Edit form<#### No longer handles summary for peptide quantif
	###################
	my ($quantifName,$focus,$quantifStatus,$quantifMethodName,$quantifMethodDes)=$dbh->selectrow_array("SELECT Q.NAME,Q.FOCUS,Q.STATUS,QM.NAME,QM.DES FROM QUANTIFICATION Q,QUANTIFICATION_METHOD QM WHERE Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=$itemID");
	my $quantifString = ($quantifStatus==-2)? '<FONT color="#DD0000">Failed</FONT> (Click on "Monitor Quantification(s)" for more information)' : ($quantifStatus==-1)? 'Not launched yet' : ($quantifStatus==0)? 'On-going' : 'Finished';
	$focus=ucfirst($focus).'s';
	$focus.=' and fragment ions' if $quantifMethodName=~/Swath/i;
	
	$dbh->disconnect;
	
	#######################
	####>Starting HTML<####
	#######################
	my ($lightColor,$darkColor)=&promsConfig::getRowColors;
	my $bgColor=$lightColor;
	my %itemIcones=&promsConfig::getItemIcones;
	
	print header(-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Quantification</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
function cancelAction() {
	top.promsFrame.selectedAction='summary';
	top.promsFrame.optionFrame.selectOption();
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" >
<CENTER>
<FONT class="title">Editing Quantification <FONT color="#DD0000">$quantifName</FONT><BR><BR>
<FORM name="exportForm" method="post">
<INPUT type="hidden" name="itemID" value="$itemID">
<INPUT type="hidden" name="branchID" value="$branchID">
<INPUT type="hidden" name="PROJECT_ID" value="$projectID">
<INPUT type="hidden" name="ACT" value="update">

<TABLE bgcolor=$darkColor cellpadding=2 width=800>
	<TR><TH align=right valign=top>&nbsp;Name :</TH><TD bgcolor="$lightColor"><INPUT type="text" name="quantiName" value="$quantifName" size=50/></TD></TR>
	<TR><TH align=right valign=top>&nbsp;Method :</TH><TD bgcolor="$lightColor">&nbsp;$quantifMethodName ($quantifMethodDes)</TD></TR>
	<TR><TH align=right valign=top>&nbsp;Focus :</TH><TD bgcolor="$lightColor">&nbsp;$focus</TD></TR>
	<TR><TH align=right valign=top width=170>&nbsp;Status :</TH><TD align=left bgcolor="$lightColor">&nbsp;$quantifString</TD></TR>
	
|;
	if ($quantifStatus >= 1 && $userInfo[1]=~/bioinfo|mass/) {
		my ($selOpt2,$selOpt3)=($quantifStatus==2)? (' selected','') : ($quantifStatus==3)? ('',' selected') : ('','');
		print "<TR><TH align=right valign=top>&nbsp;Visibility :</TH><TD align=left bgcolor=\"$lightColor\"><SELECT name=\"quantiStatus\"><OPTION value=\"1\">Public</OPTION><OPTION value=\"2\"$selOpt2>Hide from collaborators</OPTION>";
		print "<OPTION value=\"3\"$selOpt3>For bioinformaticians only</OPTION>" if $userInfo[1] eq 'bioinfo';
		print "</SELECT></TD></TR>\n";
	}
	
	print qq
|	<TR><TH colspan=2><INPUT type="submit" name="save" value=" Save ">&nbsp;&nbsp;<INPUT type="button" value=" Cancel " onclick="cancelAction()"></TD></TR>
</TABLE>
</FORM>
<BR><BR>
</BODY>
</HTML>
|;
}

####>Revision history<####
# 1.2.8 [ENHANCEMENT] Added extra parameter to &deleteQuantification to indicated a targetted deletion (PP 28/05/20)
# 1.2.7 [BUGFIX] Restrict deletability to quantifications without children (PP 12/11/19)
# 1.2.6 [ENHANCEMENT] Add form check for multiple quantifications deletion (PP 09/11/19)
# 1.2.5 [ENHANCEMENT] Allow deletion of on-going quantifications (VS 21/10/19)
# 1.2.4 [FEATURE] Allow multiple quantifications deletion (VS 04/10/19)
# 1.2.3 [FEATURE] Editabble quantification visibility (PP 26/09/19)
# 1.2.2 [FEATURE] Removed peptide quantification summary (PP 29/08/19)
# 1.2.1 Compatible with deletion of multiple quantifications at once (PP 05/09/18)
# 1.2.0 Major code update and cleaning. 'summary' of design-quantif moved to showProtQuantification.cgi (PP 25/07/16)
# 1.1.5 Delete action updated for new BD tables used by modification quantification (PP 06/11/14)
# 1.1.4 Minor bug corrections for quantifications not linked to a DESIGN like multi XIC-extractions (GA 08/10/14)
# 1.1.3 Correct bug for error log file and move this option to watchQuantifications.cgi (GA 27/06/14)
# 1.1.2 Uses rmtree instead of remove_tree (PP 10/03/14)
# 1.1.1 system command removal (PP 08/11/13)
# 1.1.0 PEPTIDE_SET no longer used (PP 16/09/13)
# 1.0.9 Add edition of quantiName & UTF-8 header (GA 13/08/13)
# 1.0.8 Fix moved Rout file deletion (PP 03/07/13)
# 1.0.7 Adding export error option (GA 11/03/13)
# 1.0.6 Code cleaning (PP 04/01/13)
# 1.0.5 Minor modification ACT=experiment to make it homogene with openProject.cgi (GA 18/10/12)
# 1.0.4 Removal of redundant code / modification of edit option (GA 20/04/12)
# 1.0.3 No need to update ID_DESIGN for PARENT_QUANTIFICATION according to launchQuantification.cgi modification (GA 12/04/12)
# 1.0.2 Add nochmod.so in order to remove TNPQ generated files in deleteQuantification function (GA 16/12/11)
# 1.0.1 Update the presentation of the frame - removal of labelled part (GA 24/11/2011)
# 1.0.0 New script to handle quantifications (save, edit and delete options)
