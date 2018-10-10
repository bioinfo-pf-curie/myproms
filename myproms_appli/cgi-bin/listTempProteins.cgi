#!/usr/local/bin/perl -w

################################################################################
# listTempProteins.cgi      2.2.1                                              #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Lists the proteins being validated                                           #
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
#use strict;
use promsConfig;
use promsMod;

#print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;

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
<TITLE>Protein List</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TD {font-weight:bold;}
A.good:link {color: black} A.good:visited {color:black}
A.poor:link {color: #707070} A.poor:visited {color:#707070}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
|;
&promsMod::popupInfo(200,'abs',0,'rel',25);
print qq
|</SCRIPT>
</HEAD>
<BODY topmargin="0" leftmargin="0" rightmargin="0">
<DIV id="waitDiv" style="position:absolute;top:40px;width:100%">
<TABLE align=center>
<TR><TH>Fetching data<SPAN id="waitSpan">...</SPAN></TH></TR>
<TR><TH><IMG src="$promsPath{images}/scrollbarGreen.gif" style="width:80%;height:4px"/></TH></TR>
</TABLE>
</DIV>
|;

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

########################
####>>>Parameters<<<####
########################
my $analysisID=param('ID');
my $msType=param('MSTYPE');
my $selectedProtID=param('SEL');
my $selectedSort=param('SORT');
my $showFiltered=param('SHOWFILT');
my $updateRank=(param('UPRANK'))? param('UPRANK') : 0;
my ($page,$maxPage)=split(/-/,param('PAGE'));
my $listMGr=param('MG');
my $displayPTM=param('DISPLAYPTM')? param('DISPLAYPTM') : 0;
my $showPTM=param('SHOWPTM')? param('SHOWPTM') : 'Indifferent';
my $PTMFilter=param('PTMFILTER')? param('PTMFILTER') : 'Indifferent';
(my $fixRes)=($PTMFilter =~ s/^(.):://g) ? $1 : '';

###################################
####>Fetching list of proteins<####
###################################

####>Converting match group string to list<####
my %matchGroups;
my @tmpList=split(',',$listMGr);
foreach my $subStrg (@tmpList) {
	if ($subStrg=~/-/) {
		my ($start,$end)=split(/-/,$subStrg);
		foreach my $mg ($start..$end) {$matchGroups{$mg}=1};
	}
	else {$matchGroups{$subStrg}=1;}
}

####>Finding proteins in range & corresponding match groups<####
#my $sthP=$dbh->prepare("SELECT ID_PROT_VALID,SEL_STATUS,MATCH_GROUP FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$analysisID AND IDENTIFIER NOT LIKE 'DECOY_%' ORDER BY MAX_SCORE DESC,ID_PROT_VALID ASC");
#$sthP->execute || die $sthP->errstr;
#my $refProtData=$sthP->fetchall_arrayref; # reference to an array
#$sthP->finish;

####>Finding number of databanks searched<####
my ($numDatabanks)=$dbh->selectrow_array("SELECT COUNT(*) FROM ANALYSIS_DATABANK WHERE ID_ANALYSIS=$analysisID");

####>Finding proteins in range & corresponding match groups<####
my @sortRules=split(/,/,$selectedSort);
my $sortString='';
foreach my $rule (@sortRules) {
	my ($sortedItem,$order)=split(/:/,$rule);
	$sortString.="$sortedItem $order,";
}
$sortString.='ID_PROT_VALID ASC';
my $filterString=($showFiltered)? '' : 'SEL_STATUS > -3 AND';
my $sthD=$dbh->prepare("SELECT ID_PROT_VALID,IDENTIFIER,PROT_DES,MW,PROT_LENGTH,ORGANISM,SEL_STATUS,NUM_MATCH,MAX_MATCH,CONF_LEVEL,MATCH_GROUP,DB_RANK FROM PROTEIN_VALIDATION WHERE $filterString ID_ANALYSIS=$analysisID AND IDENTIFIER NOT LIKE 'DECOY_%' ORDER BY $sortString");
$sthD->execute || die $sthD->errstr;
my $refData=$sthD->fetchall_arrayref; # reference to an array
$sthD->finish;

###################################
####>PTM filtering and/or info<####
###################################
my $qPTMFilter=quotemeta($PTMFilter);
my $numProtMatched=0;
my $usedDataList=$refData; # default; Changed if PTM filtering
if ($PTMFilter ne 'Indifferent' || $displayPTM) {
	my (%identiferMatch,%queryData);
	my ($numModifs)=$dbh->selectrow_array("SELECT COUNT(*) FROM ANALYSIS_MODIFICATION WHERE ID_ANALYSIS=$analysisID");
	my $sthQueryRank=$dbh->prepare("SELECT IDENTIFIER,QUERY_NUM,PEP_RANK FROM RANK_PROTEIN_MATCH WHERE ID_ANALYSIS=$analysisID");
	$sthQueryRank->execute;
	while (my ($identifier,$queryNum,$rank)=$sthQueryRank->fetchrow_array) {push @{$identiferMatch{$identifier}},[$queryNum,$rank];}
	$sthQueryRank->finish;
	my ($maxRank)=$dbh->selectrow_array("SELECT MAX_RANK FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
	$maxRank=($msType ne 'MIS')? 10 : ($maxRank)? $maxRank : &promsConfig::getMaxRank;
	my $rankStrg='INFO_PEP1';
	foreach my $rk (2..$maxRank) {$rankStrg.=",INFO_PEP$rk";}
	my $sthRankData=$dbh->prepare("SELECT QUERY_NUM,$rankStrg FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID AND QUERY_NUM > 0 AND VALID_STATUS > -4");
	$sthRankData->execute;
	while (my ($queryNum,@pepRanks)=$sthRankData->fetchrow_array) {@{$queryData{$queryNum}}=@pepRanks;}
	$sthRankData->finish;

	$usedDataList=[]; # reset to empty array ref
	PROT:foreach my $refProtID (@{$refData}) {
		next unless $matchGroups{$refProtID->[10]};
		my $PTMString='';
		my $hasSelPTM=0;
		my %varModifications;
		RANK: foreach my $refMatch (@{$identiferMatch{$refProtID->[1]}}) { # $identifier
			my ($queryNum,$rank)=@{$refMatch};
			if ($fixRes) { # fix moficiation only
				my ($sequence)=($queryData{$queryNum}[$rank-1]=~/SEQ=(\w+)/);
				$hasSelPTM=1 if $sequence =~ /$fixRes/;
				last RANK if (!$displayPTM && $hasSelPTM);
				next if (!$displayPTM);
			}
			my ($varMod)=($queryData{$queryNum}[$rank-1]=~/VMOD=([^,]+)/);
			next unless $varMod;
			next PROT if $PTMFilter eq 'None';
			my @modifs=split(/\s*\+\s*/,$varMod);
			foreach my $i (1..$#modifs) { # VMOD string starts with ' + ' => [0] is empty
				$modifs[$i]=~s/:[\d+\.]+\)/\)/; # \) needed to prevent ":13C" in "Label:13C(6)15N(2) (K)" to be matched
				$varModifications{$modifs[$i]}=1;
				$hasSelPTM=1 if ($PTMFilter eq 'Any' || ($PTMFilter ne 'Indifferent' && $modifs[$i]=~/$qPTMFilter/));
				last RANK if (!$displayPTM && $hasSelPTM); # selected modif is found. No need to keep looking
			}
			last if scalar keys %varModifications == $numModifs; # all possible modifs have been found. No need to keep looking
			last if (!$displayPTM && $PTMFilter eq 'Any');
		}
		my $hasPTMs=scalar keys %varModifications;
		if ($displayPTM) {
			$PTMString='<BR><B>PTM(s):</B>';
			if ($hasPTMs) {$PTMString.='<BR>&nbsp;-'.(join('<BR>&nbsp;-',sort keys %varModifications));}
			else {$PTMString.=' None';}
		}

		####>Print or not proteins according to their peptide PTMs
		if ($PTMFilter eq 'Any') {
			next unless $hasPTMs; # $PTMString;
		}
		elsif ($PTMFilter eq 'None') {
			next if $hasPTMs;
		}
		elsif ($PTMFilter ne 'Indifferent') {
			next unless $hasSelPTM;
		}
		push @{$refProtID},$PTMString;
		push @{$usedDataList},$refProtID;
	}
}

$dbh->disconnect;

&updateWaitBox;

##############################################################
####>Checking selectedProtID and corresponding link class<####	in case !provided (page changed) protein has become hidden: eg. hiding filtered proteins
##############################################################
my $sortByMatchGr=($selectedSort=~/MATCH_GROUP/)? 1 : 0;
my ($matchSelected,$selectedGroup); # expected selection
my ($firstProtID,$firstGroup); # 1st backup selection (1st unverified protein)
my ($newProtID,$newGroup); # 2nd backup selection (1st protein of the list)
foreach my $refProtID (@{$usedDataList}) { # $refData
	next unless $matchGroups{$refProtID->[10]};
	if ($refProtID->[0]==$selectedProtID) {
		$selectedGroup=$refProtID->[10];
		$matchSelected=1;
		last;
	}
	#>data for 1st protein in list
	unless ($firstProtID) {
		$firstProtID=$refProtID->[0];
		$firstGroup=$refProtID->[10];
	}
	#>data for 1st unverified protein in list
	if ($refProtID->[6]==-1 && !$newProtID) {
		$newProtID=$refProtID->[0];
		$newGroup=$refProtID->[10];
	}
}
unless ($matchSelected) { # default selected protein is not in list
	###>Finding a new selectedProtID<###
	if ($newProtID) { # another not-verified protein exists
		$selectedProtID=$newProtID;
		$selectedGroup=$newGroup;
	}
	else { # all proteins are verified: pick 1st in list
		($selectedProtID,$selectedGroup)=($firstProtID)? ($firstProtID,$firstGroup) : (0,0);
	}
	$updateRank=1; # set flag to update listRanks.cgi
}
&updateWaitBox;

################################
####>Re-starting JavaScript<####
################################
print "<SCRIPT language=\"JavaScript\">\n";

####>Passing Match Groups and Confidence Level to javascript<####
##>Match Group
print "var matchGroup = new Array;\n";
my $curGroup=0;
my $grSize;
my %pageDataRef;
foreach my $refProtID (sort{$a->[10] <=> $b->[10]} @{$usedDataList}) { # $refData
	next unless $matchGroups{$refProtID->[10]};
	$pageDataRef{$refProtID}=1;
	if ($refProtID->[10] > $curGroup) {
		$curGroup=$refProtID->[10];
		print "matchGroup['$curGroup'] = new Array();\n";
		$grSize=0;
	}
	print "matchGroup['$curGroup'][$grSize]=$refProtID->[0];\n";
	$grSize++;
}
print "matchGroup['0'] = new Array();\n" unless $curGroup;

##>Bad Confidence
print "var badConf = new Array();\n";
foreach my $refProtID (keys %pageDataRef) {
	print "badConf['$refProtID->[0]']=1;\n" if ($refProtID->[9] && $refProtID->[9]==1);
}

print qq
|function openProtein(protId,protGroup,updateRank) {
	if (!protId) {
		parent.rankFrame.location="$promsPath{html}/nothing.html";
		parent.spectrumFrame.location="$promsPath{html}/nothing.html";
		return;
	}
	if (protId==selectedProtId) { // 1st time a protein is selected
		document.getElementById(protId).focus(); // scrolls window until protein is visible
	}
	// Set previous match group members to default color(s)
	for (var i=0;i<matchGroup[selectedGroup].length;i++) {
		var oldGrProt=document.getElementById(matchGroup[selectedGroup][i]);
		if (oldGrProt) {
			if (badConf[matchGroup[selectedGroup][i]]) {oldGrProt.style.color='#707070';}
			else {oldGrProt.style.color='#000000';}
		}
	}
	// Highlight selected protein and other group members
	var newProt=parent.queryFrame.document.getElementById(protId);
	newProt.style.color='#DD0000'; //style.background = '#DD0000'; //value="Ca marche"; //
	for (var i=0;i<matchGroup[protGroup].length;i++) {
		var newGrProt=document.getElementById(matchGroup[protGroup][i]);
		if (newGrProt && newGrProt != newProt) {
			newGrProt.style.color='#0000DD'; //'#DD7777';
		}
	}
	selectedProtId=protId;
	selectedGroup=protGroup;
	// Update rankFrame
	if (updateRank) {
		parent.selectedProtId=protId;
		var pepRangeString=(parent.rankFrame.document.rankForm && protId==parent.rankFrame.document.rankForm.ID.value)? parent.rankFrame.pepRangeString : '';
		parent.rankFrame.location="$promsPath{cgi}/listRanks.cgi?ANAID=$analysisID&MSTYPE=$msType&ITEM=protein&SORT=$selectedSort&PAGE=$page-$maxPage&RANGE=$listMGr&ID="+protId+"&SHOWFILT=$showFiltered&sortType="+parent.listRankSort+"&Dlt="+parent.listRankDlt+"&maxDisp="+parent.maxDisp+"&maxSrch="+parent.maxSrch+"&varValue="+parent.varValue+pepRangeString+"&multiAna="+parent.multiAnaScan;
	}
}
function changePage(newPage) {
	parent.selProteinPage=newPage;
	window.location="./listTempProteins.cgi?ID=$analysisID&MSTYPE=$msType&SEL=0&SORT=$selectedSort&UPRANK=1&SHOWFILT=$showFiltered&PAGE="+newPage+"-$maxPage&MG="+parent.proteinPages[newPage-1]+"&DISPLAYPTM="+parent.displayPTM+"&PTMFILTER="+parent.ptmFilter;
}
function changeDisplayModif(newDisplay){
	parent.displayPTM=newDisplay; // 0 or 1
	window.location="./listTempProteins.cgi?ID=$analysisID&MSTYPE=$msType&SEL=0&SORT=$selectedSort&UPRANK=1&SHOWFILT=$showFiltered&PAGE="+parent.selProteinPage+"-$maxPage&MG="+parent.proteinPages[parent.selProteinPage-1]+"&DISPLAYPTM="+parent.displayPTM+"&PTMFILTER="+parent.ptmFilter;
}
var selectedProtId=$selectedProtID;
parent.selectedProtId=$selectedProtID; // in case coming from manual page change
var selectedGroup=$selectedGroup;
document.getElementById('waitDiv').style.display='none';
</SCRIPT>
<CENTER>
|;
&updateWaitBox;
&promsMod::printPageMenu($page,$maxPage);
my ($color1,$color2)=&promsConfig::getRowColors;
my $bgColor=$color2;
my ($newDisplay,$tag,$infoText)=($displayPTM)? (0,'-','hide') : (1,'+','show');
print qq
|<TABLE border=0 cellspacing=0>
<TR bgcolor="$color2"><TH class="rbBorder" width=190><A href="javascript:changeDisplayModif($newDisplay)" onmouseover=\"popup('<B>Click to $infoText PTMs in popups</B>')\" onmouseout=\"popout()\">Protein<SUP>$tag</SUP></A></TH><TH class="bBorder" width=60>&nbsp;Matches&nbsp;</TH></TR>
|;
my $curMatchGr=0;
my $numProtFetched=0;
my $numProtInPage=scalar keys %pageDataRef;
foreach my $refProtID (@{$usedDataList}) { # $refData
	next unless $pageDataRef{$refProtID}; # must be in current page
	$numProtFetched++;
	my ($protID,$identifier,$protDes,$protMW,$protLength,$protOrg,$selStatus,$numMatch,$maxMatch,$confLevel,$matchGroup,$dbRank,$PTMString)=@{$refProtID};
	$PTMString='' unless $PTMString;
	my $image;
	if ($selStatus<=-2) {
		$image='lightBlue1.gif';
	}
	elsif ($selStatus==-1) {
		$image='lightGray1.gif';
	}
	elsif ($selStatus==0) {
		$image='lightRed1.gif';
	}
	elsif ($selStatus==1) {
		$image='lightYellow1.gif';
	}
	elsif ($selStatus==2) {
		$image='lightGreen1.gif';
	}
	$protDes=&promsMod::HTMLcompatible($protDes);
	my $lengthString=($protLength)? "$protLength aa." : 'unknown.';
	my $mwString=($protMW)? sprintf "%.1f KDa.",$protMW/1000 : 'unknown.';
	my $class=($confLevel && $confLevel==1)? 'poor' : 'good';
	$protOrg=($protOrg)?$protOrg:'unknown';
	if ($sortByMatchGr) {
		if ($curMatchGr!=$matchGroup) {
			$bgColor=($bgColor eq $color1)? $color2 : $color1;
			$curMatchGr=$matchGroup;
		}
	}
	else {$bgColor=($bgColor eq $color1)? $color2 : $color1;}
	print "<TR class=\"list\" bgcolor=\"$bgColor\">";
	print "<TD>&nbsp<IMG src=\"$promsPath{images}/$image\" hspace=0 border=0 height=11 width=11>&nbsp";
	print "<A class=\"$class\" id=\"$protID\" href=\"javascript:openProtein($protID,$matchGroup,1)\" onmouseover=\"popup('<B>Protein:</B> $identifier<BR><B>Description:</B> $protDes.<BR><B>Size:</B> $lengthString<BR><B>Mass:</B> $mwString<BR><B>Organism:</B> <I><U>$protOrg</U></I>.$PTMString')\" onmouseout=\"popout()\">";
	my $identifierStrg=($numDatabanks>1)? $dbRank.'::'.$identifier : $identifier;
	$identifierStrg=&promsMod::resize($identifierStrg,20);
	if ($selStatus==-3) {print "<I>$identifierStrg</I>";} else {print "$identifierStrg";}
	print "</A></TD>\n";
	print "<TD align=center>$numMatch/$maxMatch</TD></TR>\n";
	last if $numProtFetched==$numProtInPage;
}
if (!$numProtFetched) {
	print "<TR bgcolor=\"$color1\"><TH align=left colspan=2>&nbsp;No proteins found&nbsp;</TH></TR>\n";
}
print "</TABLE>\n";
&promsMod::printPageMenu($page,$maxPage);
print qq
|<HR width="90%">
</CENTER>
<B>
&nbsp;<U>Validation status :</U><BR>
&nbsp;<IMG src="$promsPath{images}/lightGray1.gif">&nbsp;: not verified.<BR>
&nbsp;<IMG src="$promsPath{images}/lightRed1.gif">&nbsp;: no matching peptides.<BR>
&nbsp;<IMG src="$promsPath{images}/lightYellow1.gif">&nbsp;: partially verified.<BR>
&nbsp;<IMG src="$promsPath{images}/lightGreen1.gif">&nbsp;: fully verified.<BR>
&nbsp;<IMG src="$promsPath{images}/lightBlue1.gif">&nbsp;: excluded from list.<BR>
|;
print "&nbsp<IMG src='$promsPath{images}/lightBlue1.gif'>&nbsp;: <I>filtered out.</I><BR>\n" if $showFiltered;
print qq
|&nbsp;<FONT color="#DD000">Selected protein.</FONT><BR>
&nbsp;<FONT color="#0000DD">Same match group.</FONT><BR>
&nbsp;<FONT color="#707070">Poor confidence level.</FONT>
</B>

<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">
setPopup();
|;
if ($selectedProtID) { # Display protein peptides
	print qq
|openProtein($selectedProtID,$selectedGroup,$updateRank);
|;
}
elsif (scalar keys @{$refData}) { # Automatically go to next page
	my $newPage=$page+1;
	$newPage=1 if $newPage > $maxPage;
	print qq
|if (parent.ptmInitPage == $newPage) { // All pages were already scanned. Stop at initial page
	parent.rankFrame.location="$promsPath{html}/nothing.html";
	parent.spectrumFrame.location="$promsPath{html}/nothing.html";
}
else { // Go to next page
	window.location="./listTempProteins.cgi?ID=$analysisID&MSTYPE=$msType&SEL=0&SORT=$selectedSort&UPRANK=1&SHOWFILT=$showFiltered&PAGE=$newPage-$maxPage&MG="+parent.proteinPages[$newPage-1]+"&DISPLAYPTM="+parent.displayPTM+"&PTMFILTER="+parent.ptmFilter;
}
|;
}
else { # No proteins at all in Analysis: Do nothing
	print qq
|parent.rankFrame.location="$promsPath{html}/nothing.html";
parent.spectrumFrame.location="$promsPath{html}/nothing.html";
|;
}

print qq
|</SCRIPT>
</BODY>
</HTML>
|;


############## subroutine ###############
sub updateWaitBox {
	print "<SCRIPT language=\"JavaScript\">document.getElementById('waitSpan').innerHTML+='.';</SCRIPT>\n";
}
sub updateWaitBoxString {
	return "<SCRIPT language=\"JavaScript\">document.getElementById('waitSpan').innerHTML+='.';</SCRIPT>";
}

####>Revision history<####
# 2.2.1 Add fix modification filter (GA 18/08/17)
# 2.2.0 Code reorganisation to speed up PTM filtering/popup info (PP 17/06/16)
# 2.1.7 Minor modif to $protOrg to avoid warnings if unknown organism (GA 30/10/14)
# 2.1.6 Added multiAna as listRanks parameter (PP 26/09/13)
# 2.1.5 GPL license (PP 23/09/13)
# 2.1.4 Multi-databank flag on proteins (PP 11/12/12)
# 2.1.3 header(-'content-encoding'=>'no') (PP 01/07/11)
# 2.1.2 Added wait message during data processing (PP 09/06/11)
# 2.1.1 Correction bug direct link to a selected protein from sequenceView (PP 01/04/11)
# 2.1.0 Display PTM in the validation process + filter of proteins on the PTM
