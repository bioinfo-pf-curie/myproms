#!/usr/local/bin/perl -w

################################################################################
# listInstruments.cgi               1.0.3                                      #
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
use promsConfig;
use promsMod;
use strict ;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;

#print header; warningsToBrowser(1); # DEBUG
#############################
####>Connecting to ProMS<####
#############################
my $dbh=&promsConfig::dbConnect;

############################################
####>Fetching list instrument in use... <####
############################################
my %instrumentList ;

my $sth=$dbh->prepare("SELECT ID_INSTRUMENT, NAME, DES, USE_STATUS FROM INSTRUMENT");
$sth->execute();
while (my ($id,$name,$des, $useStatus)=$sth->fetchrow_array()) {
	$instrumentList{$id}{'name'} = $name ;
	$instrumentList{$id}{'des'} = $des ;
	$instrumentList{$id}{'useStatus'} = $useStatus ;
}
$sth->finish();
$dbh->disconnect();

#######################
####>Starting HTML<####
#######################
my ($light,$dark)=&promsConfig::getRowColors;
my $bgColor=$light;
print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1); #NO DEBUG
print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TD {font-size:13px;}
</STYLE><SCRIPT LANGUAGE=\"JavaScript\">
function deleteInstrument(instrumentId,name) {
	if (confirm('Delete instrument: '+name+'?   ')) {
		window.location='./editInstrument.cgi?ACT=delete&ID='+instrumentId;
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
<FONT class="title">List of Instrument Definitions</FONT>&nbsp;&nbsp;<INPUT type="button" class="title2"value="Add new Instrument" onclick="window.location='./editInstrument.cgi?ACT=add'">
<BR><BR>
<TABLE border=0 cellspacing=0>
|;

###>Looping through databanks<###
foreach my $id (sort{lc($instrumentList{$a}{'name'}) cmp lc($instrumentList{$b}{'name'})} keys %instrumentList) {
	print qq
|<TR bgcolor=$bgColor>
<TD>&nbsp&nbsp</TD>
<TD width=600><FONT style="font-size:18px;font-weight:bold;">$instrumentList{$id}{name}</FONT>
<BR><B>&nbsp&nbsp Description :</B> $instrumentList{$id}{des}
<BR><B>&nbsp&nbsp In Use :</B> $instrumentList{$id}{useStatus}
</TD>
<TH width=100><INPUT type="button" value="Edit" style="width:75px" onclick="window.location='./editInstrument.cgi?ACT=edit&ID=$id'"><BR>
<INPUT type="button" value="Delete" style="width:75px" onclick="deleteInstrument($id,'$instrumentList{$id}{name}')">
</TH>
</TR>
|;
	$bgColor=($bgColor eq $light)? $dark : $light;
}

print qq
|</TABLE>
</CENTER>
</BODY>
</HTML>
|;

####>Revision history<####
# 1.0.3 Color on main "Go to" options (PP 25/03/16)
# 1.0.2 GPL license (PP 23/09/13)
# 1.0.1 Called by promsMain.cgi (PP 13/12/12)
