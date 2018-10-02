#!/usr/local/bin/perl -w

################################################################################
# manageQuantification.cgi               1.2.1                                 #
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
my $action=(param ('ACT'))? param ('ACT') : 'edit'; # edit, update, delete, (summary of peptide quantif only) (summary of design quantifs moved to showProtQuantification.cgi)
my $quantifID=param('ID') ? param('ID') : param('itemID') ;
my $projectID=param('PROJECT_ID');
my $branchID=param('branchID');


##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

################
####>Delete<####
################
if ($action eq 'delete') {
	my @quantifIdList=param('quantifIdList') || ($quantifID);
	my ($experimentID)=$dbh->selectrow_array("SELECT DISTINCT(ID_EXPERIMENT) FROM SAMPLE S, ANALYSIS A, ANA_QUANTIFICATION AQ WHERE S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=$quantifIdList[0]");

	print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<HEAD>
<TITLE>Deleting Quantification</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" >
<CENTER><FONT class="title">Deleting Quantifications</FONT></CENTER>
<BR><BR><BR>
<FONT class="title2">
|;
	my $designID;
	my $sthQ=$dbh->prepare("SELECT ID_DESIGN,NAME FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
	foreach my $qID (@quantifIdList) {
		$sthQ->execute($qID);
		($designID,my $qName)=$sthQ->fetchrow_array;
		print "&nbsp;-Deleting '$qName'...";
		&promsQuantif::deleteQuantification($dbh,$projectID,$qID);
		$dbh->commit;
		print " Done.<BR>\n";
	}
	$sthQ->finish;
	$dbh->disconnect;
	print "</FONT>\n";
	sleep 2;
	my $branchIDStg=($designID)? "&DESIGN=DESIGN:$designID&branchID=design:$designID" : "&branchID=experiment:$experimentID";
	print qq
|<SCRIPT LANGUAGE="JavaScript">
top.promsFrame.selectedAction='summary';
parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID$branchIDStg&ACT=experiment&EXPERIMENT=EXPERIMENT:$experimentID&ISNAVFRAME=0&VIEW=quanti";
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}
#####################
####>Update Name<#### (following 'edit')
#####################
elsif ($action eq 'update') {
	my $name=param('quantiName');
	my $sthUp=$dbh->prepare("UPDATE QUANTIFICATION SET NAME=? WHERE ID_QUANTIFICATION=$quantifID");
	$sthUp->execute($name);
	$sthUp->finish;
	$dbh->commit;

	my ($experimentID)=$dbh->selectrow_array("SELECT DISTINCT(ID_EXPERIMENT) FROM SAMPLE S, ANALYSIS A, ANA_QUANTIFICATION AQ WHERE S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=$quantifID");
	my ($item,$parentID,$childID)=split(':',$branchID);# For XIC extractions linked to a TnPQ or Pep-Ratio quantification
	$childID=$parentID unless $childID;
	my ($designID)=$dbh->selectrow_array("SELECT ID_DESIGN FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$childID");

	$dbh->disconnect;

	my $branchIDStg=($designID)? "&DESIGN=DESIGN:$designID&branchID=quantification:$childID" : "&branchID=experiment:$experimentID";
	print header(-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<HEAD>
<SCRIPT LANGUAGE="JavaScript">
top.promsFrame.selectedAction='summary';
parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID$branchIDStg&ACT=experiment&EXPERIMENT=EXPERIMENT:$experimentID&ISNAVFRAME=0&VIEW=quanti";
</SCRIPT>
</HEAD>
</HTML>
|;
	exit;
}


##############################
####>Summary or Edit form<####
##############################
my ($quantifName,$focus,$quantifStatus,$quantifMethodName,$quantifMethodDes)=$dbh->selectrow_array("SELECT Q.NAME,Q.FOCUS,Q.STATUS,QM.NAME,QM.DES FROM QUANTIFICATION Q,QUANTIFICATION_METHOD QM WHERE Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=$quantifID");

$dbh->disconnect;

#######################
####>Starting HTML<####
#######################
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my $bgColor=$lightColor;
my %itemIcones=&promsConfig::getItemIcones;

print header(-charset=>'utf-8'); warningsToBrowser(1);
print qq
|<HEAD>
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
<FONT class="title">Quantification <FONT color="#DD0000">$quantifName</FONT></FONT><BR><BR>
|;

my $quantifString = ($quantifStatus==-2)? '<FONT color="#DD0000">Failed</FONT> (Click on "Monitor Quantification(s)" for more information)' : ($quantifStatus==-1)? 'Not launched yet' : ($quantifStatus==0)? 'On-going' : 'Finished';
my $nameString;
if ($action eq 'edit') {
	$nameString="<INPUT type=\"text\" name=\"quantiName\" value=\"$quantifName\" size=50/>";
	print qq
|<FORM name="exportForm" method="post">
<INPUT type="hidden" name="itemID" value="$quantifID">
<INPUT type="hidden" name="branchID" value="$branchID">
<INPUT type="hidden" name="PROJECT_ID" value="$projectID">
<INPUT type="hidden" name="ACT" value="update">
|;
}
else {$nameString='&nbsp;'.$quantifName;}
$focus=ucfirst($focus).'s';
$focus.=' and fragment ions' if $quantifMethodName=~/Swath/i;
print qq
|<TABLE bgcolor=$darkColor cellpadding=2 width=800>
<TR><TH align=right valign=top>&nbsp;Name :</TH><TD bgcolor="$lightColor">$nameString</TD></TR>
<TR><TH align=right valign=top>&nbsp;Method :</TH><TD bgcolor="$lightColor">&nbsp;$quantifMethodName ($quantifMethodDes)</TD></TR>
<TR><TH align=right valign=top>&nbsp;Focus :</TH><TD bgcolor="$lightColor">&nbsp;$focus</TD></TR>
<TR><TH align=right valign=top width=170>&nbsp;Status :</TH><TD align=left bgcolor="$lightColor">&nbsp;$quantifString</TD></TR>
|;

#	print qq
#|<TR><TH align=right valign=top nowrap>&nbsp;Quantifications involved</TH><TD bgcolor="$lightColor"><TABLE cellpadding=0 cellspacing=0>
#|;
#	###> Print the number of expCondition
#	#if(uc($quantiStyle) =~ /(TNPQ|PROT)/) {
#	my $sthgetQuantiInfo=$dbh->prepare("SELECT NAME FROM QUANTIFICATION Q,PARENT_QUANTIFICATION PQ WHERE PQ.ID_PARENT_QUANTIFICATION=Q.ID_QUANTIFICATION AND PQ.ID_QUANTIFICATION=$quantifID ORDER BY PAR_FUNCTION ASC");
#	$sthgetQuantiInfo->execute;
#	while (my ($parentName)=$sthgetQuantiInfo->fetchrow_array) {
#		#print "<TR><TD nowrap><TH align=right><IMG src=\"$promsPath{images}/$itemIcones{quantification}\">&nbsp;$parentName</TH></TD>";
#		print "<TR><TD nowrap><TH align=left><IMG src=\"$promsPath{images}/$itemIcones{quantification}\">&nbsp;$parentName</TH></TD>";
#		print "</TR>\n";
#	}
#	$sthgetQuantiInfo->finish;
#	#}else{
#	#	my $sthgetExpCondInfo=$dbh->prepare("SELECT NAME,EXPCONDITION_QUANTIF.QUANTIF_ELEMENT FROM EXPCONDITION_QUANTIF,EXPCONDITION WHERE EXPCONDITION_QUANTIF.ID_EXPCONDITION=EXPCONDITION.ID_EXPCONDITION AND EXPCONDITION_QUANTIF.ID_QUANTIFICATION=$quantifID ORDER BY QUANTIF_ELEMENT ASC");
#	#	$sthgetExpCondInfo->execute;
#	#	while (my ($expCondName,$quantifElement) = $sthgetExpCondInfo->fetchrow_array) {
#	#		print "<TR><TD nowrap><TH align=right><IMG src=\"$promsPath{images}/$itemIcones{expcondition}\">&nbsp;$expCondName</TH></TD>";
#	#		print "<TH align=left>&nbsp;&nbsp;&nbsp;&nbsp;&rarr;&nbsp;&nbsp;&nbsp;&nbsp;$quantifElement</TH>" if $quantifElement;
#	#		print "</TR>\n";
#	#	}
#	#	$sthgetExpCondInfo->finish;
#	#}
#	#print "\n";
#
#	print qq
#|</TABLE></TD></TR>|;

if ($action eq 'edit') {
	print "<TR><TH colspan=2><INPUT type=\"submit\" name=\"save\" value=\" Save \">&nbsp;&nbsp;<INPUT type=\"button\" value=\" Cancel \" onclick=\"cancelAction()\"></TD></TR>\n";
}
print "</TABLE>\n";
print "</FORM>\n" if $action eq 'edit';
print qq
|<BR><BR>
</BODY>
</HTML>
|;
#$dbh->disconnect;

####>Revision history<####
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