#!/usr/local/bin/perl -w

################################################################################
# editProtein.cgi      2.0.9                                                   #
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
$|=1;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use POSIX qw(strftime); # to get the time

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my $date=strftime("%Y-%m-%d %H:%M:%S", localtime);

####################
#    Parameters    #
####################
my $proteinID=param('id_prot');
my $analysisID=param('id_ana');
# my $type=param('type');
## if Clustering mode
#my $simGr=param('simGr');
my $projectID=param('id_project');

###############
#     Main    #
###############
## Connect to the database
my $dbh=&promsConfig::dbConnect;

if (param('save')) {
	my $alias=&promsMod::resize(param('alias'),50);
	my $description=&promsMod::resize(param('des'),250);
	my $comments=&promsMod::resize(param('comments'),250);
	my $organism=&promsMod::resize(param('organism'),100);
	$alias=$dbh->quote($alias);
	$description=$dbh->quote($description);
	$comments=$dbh->quote($comments);
	$organism=$dbh->quote($organism);
	$userID=$dbh->quote($userID);

	my $requete="UPDATE PROTEIN SET ALIAS=$alias,PROT_DES=$description,COMMENTS=$comments,ORGANISM=$organism,UPDATE_DATE='$date',UPDATE_USER=$userID WHERE ID_PROTEIN=$proteinID";
	$dbh->do($requete);
	$dbh->commit;
	$dbh->disconnect;

	print header;
	print qq
|<HTML><HEAD>
<SCRIPT LANGUAGE="Javascript">
window.location="./sequenceView.cgi?frame=info&id_ana=$analysisID&id_prot=$proteinID&msdata="+top.opener.showMSdata;
</SCRIPT>
</HEAD></HTML>
|;
	exit;
}

# Protein info
my ($protIdentifier,$protAlias,$protDes,$protComments,$protOrganism)=$dbh->selectrow_array("SELECT IDENTIFIER,ALIAS,PROT_DES,COMMENTS,ORGANISM FROM PROTEIN WHERE ID_PROTEIN=$proteinID");

$dbh->disconnect;

###########################
#      Starting HTML      #
###########################
my ($color1,$color2)=&promsConfig::getRowColors;
print header;
warningsToBrowser(1); # writes PERL warnings as HTML comments <!-- --> must be activated after header
print qq
|<HEAD>
<TITLE>Edit Protein</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
function checkForm(){
	if (document.editProt.alias.value){
		return true;
	}
	else{
		alert ("Type a name for protein.");
		return false;
	}
}
</SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<CENTER>
<FONT class="title">Edit Protein <FONT color="#DD0000">$protAlias</FONT></FONT>
<BR>
<BR>
<FORM name="editProt" action="./editProtein.cgi" onsubmit="return(checkForm());">
<INPUT type="hidden" name="id_ana" value="$analysisID">
<INPUT type="hidden" name="id_prot" value="$proteinID">
<INPUT type=\"hidden\" name=\"id_project\" value=\"$projectID\">
<TABLE border=0 width=744><TR><TD bgcolor=$color2>
<TABLE border=0 width=740 cellpadding=2>
<TR><TH align=right valign=top bgcolor=$color2>Identifier :</TH>
<TD bgcolor=$color1>$protIdentifier</TD>
</TR>
<TR><TH align=right valign=top bgcolor=$color2>Name in myProMS :</TH>
<TD bgcolor=$color1><INPUT type="text" name="alias" size="50" value="$protAlias"></TD>
</TR>
<TR><TH align=right valign=top bgcolor=$color2>Description :</TH>
<TD bgcolor=$color1><TEXTAREA name="des" rows="2" cols="65">$protDes</TEXTAREA></TD>
</TR>
|;
if (!defined($protComments)){$protComments='';}
print qq
|<TR><TH align=right valign=top bgcolor=$color2>Comments :</TH>
<TD bgcolor=$color1><TEXTAREA name="comments" rows="5" cols="65">$protComments</TEXTAREA></TD>
</TR>
<TR><TH align=right valign=top bgcolor=$color2>Organism :</TH>
<TD bgcolor=$color1><INPUT type="text" name="organism" size="50" value="$protOrganism"></TD>
</TR>
<TR><TD colspan=2 align=center>
<INPUT type="submit" name="save" value=" Save ">
&nbsp &nbsp &nbsp<INPUT type="button" value=" Cancel " onclick="history.back()">
</TD></TR>
</TABLE>
</TD></TR></TABLE>
</FORM>
</CENTER>
</BODY>
</HTML>
|;

####>Revision history<####
# 2.0.9 GPL license (PP 19/09/13)
