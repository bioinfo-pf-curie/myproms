#!/usr/local/bin/perl -w

################################################################################
# listQueries.cgi      2.0.7                                                   #
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

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

####################
####>Parameters<####
####################
my $analysisID=param('ID');
my $msType=param('MSTYPE');
my $selectedSort=param('SORT');
my $selectedQueryID=param('SEL');
my $updateRank=(param('UPRANK'))? param('UPRANK') : 0;
my ($page,$maxPage)=split(/-/,param('PAGE'));
my ($lowQueryID,$upQueryID)=split(/-/,param('QR'));


##################################
####>Fetching list of queries<####
##################################
my ($sortedItem,$order)=split(/:/,$selectedSort);
my $sthQ=$dbh->prepare("SELECT ID_QUERY,QUERY_NUM,VALID_STATUS FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID AND ID_QUERY>=$lowQueryID AND ID_QUERY<=$upQueryID ORDER BY $sortedItem $order,ID_QUERY ASC");
$sthQ->execute || die $dbh->errstr;
my $refData=$sthQ->fetchall_arrayref; # reference to an array
$sthQ->finish;
$dbh->disconnect;

####>Checking selectedQueryID<#### in case !provided (page changed)
unless ($selectedQueryID) {
	my $firstQueryID;
	foreach my $refQuery (@{$refData}) {
		$firstQueryID=$refQuery->[0] unless $firstQueryID;
		if ($refQuery->[2]==-1 && !$selectedQueryID) {
			$selectedQueryID=$refQuery->[0];
			last;
		}
	}
	$selectedQueryID=$firstQueryID unless $selectedQueryID;
}


#######################
####>Starting HTML<####
#######################
# A {color: black; font-weight: bold; text-decoration:none}
# A:link {color: black} A:active {color: red; font-weight: bold} A:visited {color: black}
print header(-'content-encoding'=>'no',-charset=>'utf-8');
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Query List</TITLE>
<LINK rel=stylesheet href="$promsPath{html}/promsStyle.css" type=text/css>
<STYLE type="text/css">
TD,A{font-size:12px;}
TD {font-weight: bold;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
function openQuery(queryId,updateRank) {
	if (queryId==selectedQueryId) { // 1st time a query is selected
		document.getElementById(queryId).focus(); // scrolls window until query is visible
	}
	var oldQuery=top.promsFrame.queryFrame.document.getElementById(selectedQueryId);
	oldQuery.style.color='#000000';
	var newQuery=top.promsFrame.queryFrame.document.getElementById(queryId);
	newQuery.style.color='#DD0000'; //style.background = '#DD0000';
	selectedQueryId=queryId;
	if (updateRank) {
		parent.selectedQueryId=queryId;
		parent.rankFrame.location="$promsPath{cgi}/listRanks.cgi?ANAID=$analysisID&MSTYPE=$msType&ITEM=query&SORT=$selectedSort&PAGE=$page-$maxPage&RANGE=$lowQueryID-$upQueryID&ID="+queryId+"&SHOWFILT="+parent.showFiltered+"&varValue="+parent.varValue+"&multiAna="+parent.multiAnaScan;
	}
}
function changePage(newPage) {
	parent.selQueryPage=newPage;
	window.location="./listQueries.cgi?ID=$analysisID&MSTYPE=$msType&SEL=0&SORT=$selectedSort&UPRANK=1&PAGE="+newPage+"-$maxPage&QR="+parent.queryPages[newPage-1];
}
var selectedQueryId=$selectedQueryID;
</SCRIPT>
</HEAD>
<BODY topmargin="0" leftmargin="0" rightmargin="0" onload="openQuery($selectedQueryID,$updateRank)">
<CENTER>
|;
&promsMod::printPageMenu($page,$maxPage);
my ($color1,$color2)=&promsConfig::getRowColors;
my $bgColor=$color1;
print qq
|<TABLE border=0 cellspacing=0>
<TR bgcolor="$color2"><TH class="rbBorder" width=80>Query</TH><TH class="bBorder" width=75>Validated<BR>peptides</TH></TR>
|;
foreach my $refQuery (@{$refData}) {
	my ($queryID,$queryNum,$validStatus)=@{$refQuery};
	my $validPep=$validStatus;
	my $image='lightGreen1.gif';
	my $queryString="<A id=\"$queryID\" href=\"javascript:openQuery($queryID,1)\">q$queryNum</A>";
	if ($validStatus==-4) { # no interpretation
		$validPep='-';
		$image='lightBlue1.gif';
# 		$queryString="q$queryNum";
	}
	elsif ($validStatus==-3) { # better scores exist for all interpretations
		$validPep='-';
		$image='lightRed1.gif';
	}
	elsif ($validStatus==-2) { # rejected interpretation(s) + not verified one(s)
		$validPep=0;
		$image='lightGrayRed1.gif';
	}
	elsif ($validStatus==-1) { # not verified at all
		$validPep=0;
		$image='lightGray1.gif';
	}
	elsif ($validStatus==0) {
		$image='lightRed1.gif';
	}
	print "<TR class=\"list\" bgcolor=\"$bgColor\"><TD>&nbsp<IMG src='$promsPath{images}/$image' HSPACE=0 BORDER=0 height=11 width=11>&nbsp;$queryString</TD>\n";
	print "<TD align=center>$validPep</TD></TR>\n";
	$bgColor=($bgColor eq $color1)? $color2 : $color1;
}
print "</TABLE>\n";
&promsMod::printPageMenu($page,$maxPage);
print qq
|<HR width="90%">
</CENTER>
<B>
&nbsp<U>Query Status :</U><BR>
&nbsp<IMG src='$promsPath{images}/lightGray1.gif' width=11 height=11 border=noborder>&nbsp;: not verified.<BR>
&nbsp<IMG src='$promsPath{images}/lightGreen1.gif' width=11 height=11 border=noborder>&nbsp;: validated peptide(s).<BR>
&nbsp<IMG src='$promsPath{images}/lightGrayRed1.gif' width=11 height=11 border=noborder>&nbsp<IMG src='$promsPath{images}/lightRed1.gif' width=11 height=11 border=noborder>&nbsp;: rejected peptide(s).<BR>
&nbsp<IMG src='$promsPath{images}/lightBlue1.gif' width=11 height=11 border=noborder>&nbsp;: no interpretations.
</B>
<BR><BR>
</BODY>
</HTML>
|;

####>Revision history<####
# 2.0.7 Added multiAna as listRanks parameter (PP 26/09/13)
# 2.0.6 GPL license (PP 23/09/13)
# 2.0.5 header(-'content-encoding'=>'no') (PP 01/07/11)
