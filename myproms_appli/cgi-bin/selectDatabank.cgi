#!/usr/local/bin/perl -w

################################################################################
# selectDatabank.cgi          2.2.3                                            #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Lists the databanks available in myProMS                                     #
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
use strict;
use promsConfig;
use promsMod;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my %identifierTypes=&promsMod::getIdentifierTypes;

# print header(-charset=>'utf-8'),"DEBUG<BR>\n"; warningsToBrowser(1); # DEBUG
##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

############################################
####>Fetching list databanks in use... <####
############################################
my $sthDB=$dbh->prepare("SELECT ID_DATABANK,D.NAME,VERSION_NAME,D.DES,D.COMMENTS,IDENTIFIER_TYPE,DECOY_TAG,ORGANISM,FASTA_FILE,NUM_ENTRY,DT.NAME FROM DATABANK D,DATABANK_TYPE DT WHERE USE_STATUS='yes' AND D.ID_DBTYPE=DT.ID_DBTYPE ORDER BY D.NAME ASC");
$sthDB->execute;
my $refData=$sthDB->fetchall_arrayref; # reference to an array
$sthDB->finish;
my %notDeletable;
my $sthAD=$dbh->prepare("SELECT DISTINCT ID_DATABANK FROM ANALYSIS_DATABANK AB,ANALYSIS A WHERE A.ID_ANALYSIS=AB.ID_ANALYSIS AND VALID_STATUS=-1");
$sthAD->execute;
while (my $dbID=$sthAD->fetchrow_array) {
	$notDeletable{$dbID}=1;
}
$sthAD->finish;

$dbh->disconnect;


#######################
####>Starting HTML<####
#######################
my ($light,$dark)=&promsConfig::getRowColors;
my $bgColor=$light;
print header(-charset=>'utf-8');
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>List of Databanks</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TD {font-size:13px;}
</STYLE>
<SCRIPT LANGUAGE=\"JavaScript\">
function deleteDatabank(databankId,name) {
	if (confirm('Delete databank: '+name+'?   ')) {
		window.location='./editDatabank.cgi?ACT=delete&ID='+databankId;
	}
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
<FONT class="title">List of Sequence Databanks in Use</FONT>&nbsp;&nbsp;<INPUT type="button" class="title2" value="Add new Databank" onclick="window.location='./editDatabank.cgi?ACT=add'">
<BR><BR>
<TABLE border=0 cellspacing=0>
|;
# print scalar @{$refData};
###>Looping through databanks<###
foreach my $refRow (@{$refData}) {
	my ($databankID,$name,$versionName,$des,$comments,$identType,$decoyTag,$species,$fileString,$numEntry,$dbTypeName)=@{$refRow};
	$versionName='' unless $versionName;
	$des='' unless $des;
	$comments='' unless $comments;
	my $identTypeName=($identType)? $identifierTypes{$identType} : 'Undefined';
	$identTypeName=~s/ \(.+//; # example not needed
	my $decoyString=($decoyTag)? "&nbsp;&nbsp;(<B>Decoy</B> entries tagged with  <B>$decoyTag</B>)" : '';
	$species='Multiple species' unless $species;
	my $disabStrg=($notDeletable{$databankID})? 'disabled' : '';

	#>Checking fasta file
	my $fileInfo;
	if ($fileString=~/:/) { # mascot
		my ($mascotServer,$dbankDir,$fileName)=split(':',$fileString);
		$fileInfo="$mascotServer > $dbankDir > $fileName";
	}
	else {
		$fileInfo=$fileString;
		my $fastaFile="$promsPath{banks}/db_$databankID/$fileString";
		unless (-e $fastaFile) {$fileInfo.=" <FONT color=#DD0000>(File is missing)</FONT>";}
	}

	print qq
|<TR bgcolor=$bgColor>
<TD>&nbsp&nbsp</TD>
<TD width=600><FONT style="font-size:18px;font-weight:bold;">$name</FONT>
<BR><B>&nbsp;&nbsp;&nbsp;Version:</B> $versionName
<BR><B>&nbsp;&nbsp;&nbsp;Type:</B> $dbTypeName
<BR><B>&nbsp;&nbsp;&nbsp;Identifier type:</B> $identTypeName$decoyString
<BR><B>&nbsp;&nbsp;&nbsp;Species:</B> $species
|;
	print "<BR><B>&nbsp;&nbsp;&nbsp;Description:</B> ",&promsMod::HTMLcompatible($des) if $des;
	print "<BR><B>&nbsp;&nbsp;&nbsp;Comments:</B> ",&promsMod::HTMLcompatible($comments) if $comments;
	print qq
|<BR><B>&nbsp;&nbsp;&nbsp;Sequence file:</B> $fileInfo
<BR><B>&nbsp;&nbsp;&nbsp;Number of entries:</B> $numEntry
</TD>
<TH width=100><INPUT type="button" value="Edit" style="width:75px" onclick="window.location='./editDatabank.cgi?ACT=edit&ID=$databankID'"><BR>
<INPUT type="button" value="Delete" style="width:75px" onclick="deleteDatabank($databankID,'$name')" $disabStrg>
</TH>
</TR>
|;
	$bgColor=($bgColor eq $light)? $dark : $light;
}
if (scalar @{$refData} == 0) {
	print "<TR><TH colspan=3>(No databanks in use)</TH></TR>\n";
}
print qq
|</TABLE>
<BR><BR>
</CENTER>
</BODY>
</HTML>
|;

####>Revision history<####
# 2.2.3 Description & Comments displayed only if defined (PP 04/08/17)
# 2.2.2 Color on main "Go to" options (PP 25/03/16)
# 2.2.1 Moved &getIdentifierTypes from promsConfig.pm to promsMod.pm (PP 21/03/14)
# 2.2.0 Databanks are sorted by name (PP 31/10/13)
# 2.1.9 Correction in databank deletability (PP 20/08/13)
# 2.1.8 GPL license & UTF-8 header (PP 02/07/13)
# 2.1.7 Multi-databank search (PP 18/12/12)
# 2.1.6 Called by promsMain.cgi (PP 13/12/12)
# 2.1.5 web page title added (PP 10/03/2011)
