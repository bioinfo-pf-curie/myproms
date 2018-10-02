#!/usr/local/bin/perl -w

################################################################################
# editSpectrumComments.cgi             1.1.4                                   #
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
use strict ;
use CGI ':standard';
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use promsConfig;
use promsMod;


print header; warningsToBrowser(1);
##########################################
####>Configuration & Connection to DB<####
##########################################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my $dbh=&promsConfig::dbConnect;


####>Store new comment<####
if (param('editComment') ) {
	my $appel = param('editComment');
	my $spectrumID = param('ID');
	my $rank = param('hit') if param('hit');
	my $twoSpectra = (defined(param('twoSpectra')))? param('twoSpectra') : 0;
	my $newComment=param('editTextForm') ; #if defined param('editTextForm');
	chomp ($newComment);
	$newComment =~ s/,/;/g; # to avoid wrong parsing of comments
	$newComment=&promsMod::HTMLcompatible($newComment);
	if ($appel eq 'rank' ) {
		my $sthQuery= $dbh->prepare("SELECT INFO_PEP$rank FROM QUERY_VALIDATION where  ID_QUERY=$spectrumID");
		$sthQuery->execute;
		my $infoPep = $sthQuery->fetchrow_array;
		$sthQuery->finish;
		$infoPep =~ s/COM=[^,]*,?//; #remove old comments
		$infoPep .="COM=$newComment,";
		$infoPep=$dbh->quote($infoPep);
		$dbh->do("UPDATE QUERY_VALIDATION SET INFO_PEP$rank = $infoPep WHERE ID_QUERY=$spectrumID") || die $dbh->errstr();
		if ($twoSpectra>0) {
			$newComment=$dbh->quote($newComment);
			$dbh->do("UPDATE PEPTIDE  SET COMMENTS= $newComment WHERE ID_ANALYSIS=(SELECT ID_ANALYSIS FROM QUERY_VALIDATION WHERE ID_QUERY=$spectrumID ) AND PEP_RANK = $rank AND QUERY_NUM = (SELECT QUERY_NUM FROM QUERY_VALIDATION WHERE ID_QUERY=$spectrumID)") || die $dbh->errstr();

		}
	}
	elsif ($appel eq 'ref') {
		$newComment=$dbh->quote($newComment);
		$dbh->do("UPDATE SPECTRUM  SET COMMENTS= $newComment WHERE ID_SPECTRUM =$spectrumID") || die $dbh->errstr();
	}
	elsif ($appel eq 'pep') {
		my $newCommentQuoted=$dbh->quote($newComment);
		$dbh->do("UPDATE PEPTIDE  SET COMMENTS= $newCommentQuoted WHERE ID_PEPTIDE =$spectrumID") || die $dbh->errstr();
		if ($twoSpectra>0) {
			my $sthQueryRef= $dbh->prepare("SELECT  ID_ANALYSIS, QUERY_NUM, PEP_RANK FROM PEPTIDE WHERE  ID_PEPTIDE=$spectrumID");
			$sthQueryRef->execute;
			my ($anaID, $queryNum, $rank) = $sthQueryRef->fetchrow_array;
			my $sthQuery= $dbh->prepare("SELECT INFO_PEP$rank FROM QUERY_VALIDATION where ID_ANALYSIS = $anaID AND QUERY_NUM = $queryNum");
			$sthQuery->execute;
			my $infoPep = $sthQuery->fetchrow_array;
			$sthQuery->finish;
			$infoPep =~ s/COM=[^,]*,?//; #remove old comments
			$infoPep .="COM=$newComment,";
			$infoPep=$dbh->quote($infoPep);
			$dbh->do("UPDATE QUERY_VALIDATION SET INFO_PEP$rank = $infoPep WHERE ID_ANALYSIS = $anaID AND QUERY_NUM = $queryNum ") || die $dbh->errstr();
		}
	}

	print "<HTML>\n<SCRIPT LANGUAGE=\"JavaScript\">opener.location.reload(1); window.close();</SCRIPT>\n</HTML>\n";  # //commentEditWinbows.document.write('<p align="center" ><a href="javascript:window.close();">Close Window</a></p>') ;
	$dbh->commit;
	$dbh-> disconnect;
}

####>Dialog box edit Comment<####
elsif (param('EditComm')) {
	my $appel = param('EditComm');
	my $spectrumID = param('ID');
	my $rank = param('hit') if param('hit');
	my ($comments, $analysisStatus , $twocomments);
	my ($sthQuery, $sthStatusAnalysis, $sthOtherPeptide);

	$twocomments= 0;

	if ($appel eq 'ref') {
		$sthQuery= $dbh->prepare("SELECT COMMENTS FROM SPECTRUM where  ID_SPECTRUM=$spectrumID");
	}
	elsif ($appel eq 'pep') {
		$sthQuery= $dbh->prepare("SELECT COMMENTS FROM PEPTIDE where  ID_PEPTIDE=$spectrumID");

		$sthStatusAnalysis = $dbh->prepare(" SELECT VALID_STATUS FROM ANALYSIS WHERE ID_ANALYSIS=(SELECT ID_ANALYSIS FROM PEPTIDE WHERE ID_PEPTIDE=$spectrumID)");
		$sthStatusAnalysis->execute;
		$analysisStatus = $sthStatusAnalysis->fetchrow_array;
		$sthStatusAnalysis->finish;
		if ($analysisStatus==1) {
			$twocomments = 1;
		}
	}
	elsif ($appel eq 'rank') {
		$sthQuery= $dbh->prepare("SELECT INFO_PEP$rank FROM QUERY_VALIDATION where  ID_QUERY=$spectrumID");
		$sthStatusAnalysis = $dbh->prepare("SELECT VALID_STATUS FROM ANALYSIS WHERE ID_ANALYSIS=(SELECT ID_ANALYSIS FROM QUERY_VALIDATION WHERE ID_QUERY=$spectrumID)");
		$sthStatusAnalysis->execute;
		$analysisStatus = $sthStatusAnalysis->fetchrow_array;
		$sthStatusAnalysis->finish;

		if ($analysisStatus==1) {
			$sthOtherPeptide= $dbh->prepare("SELECT COUNT(*) FROM PEPTIDE WHERE ID_ANALYSIS=(SELECT ID_ANALYSIS FROM QUERY_VALIDATION WHERE ID_QUERY=$spectrumID ) AND QUERY_NUM = (SELECT QUERY_NUM FROM QUERY_VALIDATION WHERE ID_QUERY=$spectrumID) AND PEP_RANK = $rank");
			$sthOtherPeptide->execute;
			$twocomments = $sthOtherPeptide->fetchrow_array;
			$sthOtherPeptide->finish;
		}
	}

	$sthQuery->execute;

	if ($appel eq 'rank') {
		my $infoPep = $sthQuery->fetchrow_array;
		if ($infoPep=~ /COM=([^,]*),?/) {
			$comments = $1;
			chomp ($comments);
		}
	}
	else {$comments = $sthQuery->fetchrow_array;}
	#$comments =~ s/.*COM=(.+)/$1/ if ($appel eq 'rank') ;
	$comments ="" unless $comments;
	$sthQuery->finish;

	my $twoCommentWarning = ($twocomments==0)?"":"The two spectra will be affected <BR> <INPUT type=hidden name='twoSpectra' value=$twocomments>";

	my $hitString ="";
	$hitString = "<INPUT type=hidden name='hit' value=$rank> " if $appel eq 'rank';
	print qq
|<HTML>
<HEAD>
<TITLE>Edit Comment</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
function closeEdit() {
	window.close();
}
function checkForm(Form){
	if (Form.editTextForm.value.length>240) {
		alert("Comment length is limited at 240 caracters") ;
		return false ;
	}
}
</SCRIPT>
</HEAD>
<BODY>
<FORM  method=post onsubmit="return checkForm(this) ;" >
<TABLE>
<TR><TD>
<TEXTAREA cols=54 name='editTextForm' rows=6>$comments</TEXTAREA>
</TABLE>
<INPUT type=hidden name=ID value="$spectrumID">
<INPUT type=hidden name='editComment' value=$appel>
$hitString
$twoCommentWarning
<TABLE>
<TR><TD><BUTTON name=submitButton type=submit >Save Comment </BUTTON>	</TD>
<TD><BUTTON name=cancelButton type=button onclick=closeEdit();>Cancel </BUTTON>
</TD></TR>
</TABLE>
</FORM>
</BODY>
</HTML>
|;
	$dbh-> disconnect;
}

####>Revision history<####
# 1.1.4 GPL license (PP 23/09/13)
# 1.1.3 Fixed missing comma and parsing rules for COM attribute in INFO_PEP, ',' substitution to ';' in comments (FY 20/07/12)