#!/usr/local/bin/perl -w

################################################################################
# listMyTempProteins.cgi      1.1.1                                            #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Create, edit & display a user-defined list of proteins                       #
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

#print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

########################
####>>>Parameters<<<####
########################
my $action=param('ACT');
my $analysisID=param('ID');
my @itemInfo=&promsMod::getItemInfo($dbh,'analysis',$analysisID);
my $projectID=$itemInfo[0]{'ID'};
my $myProtFile="$promsPath{valid}/ana_$analysisID/my_proteins.txt";
if ($action eq 'form') {&displayForm;}
elsif ($action eq 'submit') {&processForm;}
elsif ($action eq 'delete') {
	unlink $myProtFile;
}

my $msType=param('MSTYPE');
my $selectedProtID=param('SEL');
my $selectedSort=param('SORT');
my $showFiltered=param('SHOWFILT');
my $pageStrg=param('PAGE');
my $updateRank=(param('UPRANK'))? param('UPRANK') : 0;

###################################
####>Fetching list of proteins<####
###################################
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
my $disableString=($projectAccess =~ /bioinfo|mass|manag|super/)? '' : 'disabled'; # have full access

####>Fetching proteins in personal list file<####
if (!-e $myProtFile) {
	$dbh->disconnect;
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Private Protein List</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
parent.rankFrame.location="$promsPath{html}/nothing.html";
parent.spectrumFrame.location="$promsPath{html}/nothing.html";
var selectedProtId=0;
</SCRIPT>
</HEAD>
<BODY leftmargin="0" rightmargin="0">
<CENTER>
<B>No protein list found.</B><BR>
<INPUT type="button" value="Create Protein List" onclick="parent.rankFrame.location='$promsPath{cgi}/listMyTempProteins.cgi?ACT=form&ID=$analysisID'" $disableString>
</BODY>
</HTML>
|;
	exit;
}

####>Fetching private list of proteins<####
my ($numDatabanks)=$dbh->selectrow_array("SELECT COUNT(*) FROM ANALYSIS_DATABANK WHERE ID_ANALYSIS=$analysisID");
my (%proteinList,%proteinPage);
my $badMG=9999999999;
my $sthD=$dbh->prepare("SELECT PROT_DES,MW,PROT_LENGTH,ORGANISM,SEL_STATUS,NUM_MATCH,MAX_MATCH,CONF_LEVEL,MATCH_GROUP,SCORE,MAX_SCORE,DB_RANK FROM PROTEIN_VALIDATION WHERE ID_PROT_VALID=?");
open (MYPROT,$myProtFile);
while (<MYPROT>) {
	my ($identifier,$protID,$page)=($_=~/^(\S+)\t(\d+)\t(\d+)/);
	if ($protID) {
		$sthD->execute($protID);
		my @data=$sthD->fetchrow_array;
		next if (!$showFiltered && $data[4]<=-3);
		push @{$proteinList{$identifier}},($protID,$identifier,@data);
		$proteinPage{$protID}=$page;
	}
	else {
		push @{$proteinList{$identifier}},($protID,$identifier,'',0,0,'',-4,0,0,2,++$badMG,0,0,0);
	}
}
close MYPROT;
$sthD->finish;

$dbh->disconnect;

########>Finding proteins in range & corresponding match groups<####
my @sortRules=split(',',$selectedSort);
my ($sortKey,$sortOrder)=split(':',$sortRules[0]);
my %sortKeyPos=(
	'IDENTIFIER'=>1,
	'SEL_STATUS'=>6,
	'NUM_MATCH'=>7,
	'MAX_MATCH'=>8,
	'MATCH_GROUP'=>10,
	'SCORE'=>11,
	'MAX_SCORE'=>12
);

##############################################################
####>Checking selectedProtID and corresponding link class<####	in case !provided (page changed) protein has become hidden: eg. hiding filtered proteins
##############################################################
my $sortByMatchGr=($selectedSort=~/MATCH_GROUP/)? 1 : 0;
my ($matchSelected,$selectedGroup)=(0,undef); # expected selection
my ($firstProtID,$firstGroup); # 1st backup selection (1st unverified protein)
my ($newProtID,$newGroup); # 2nd backup selection (1st protein of the list)
foreach my $identifier (sort{&sortProt($sortKeyPos{$sortKey},$sortOrder)} keys %proteinList) {
	next unless $proteinList{$identifier}[0];
	if ($proteinList{$identifier}[0]==$selectedProtID) {
		$selectedGroup=$proteinList{$identifier}[10];
		$matchSelected=1;
		last;
	}
	#>data for 1st protein in list
	unless ($firstProtID) {
		$firstProtID=$proteinList{$identifier}[0];
		$firstGroup=$proteinList{$identifier}[10];
	}
	#>data for 1st unverified protein in list
	if ($proteinList{$identifier}[6]==-1 && !$newProtID) {
		$newProtID=$proteinList{$identifier}[0];
		$newGroup=$proteinList{$identifier}[10];
	}
}
if ($matchSelected) {
	$updateRank=($updateRank==2)? 1 : 0; # 2: prevView=query => update; 0: prevView=protein =>no need for update
}
else { # default selected protein is not in list
	###>Finding a new selectedProtID<###
	if ($newProtID) { # another not-verified protein exists
		$selectedProtID=$newProtID;
		$selectedGroup=$newGroup;
	}
	elsif ($firstProtID) { # all proteins are verified: pick 1st in list
		$selectedProtID=$firstProtID;
		$selectedGroup=$firstGroup;
	}
}
my $onloadStrg=($selectedProtID)? " onload=\"openProtein($selectedProtID,$selectedGroup,$proteinPage{$selectedProtID},$updateRank)\"" : "";

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
FONT.noMatch {text-decoration:line-through}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
|;
&promsMod::popupInfo(200,'abs',0,'rel',25);

####>Passing Match Groups and Confidence Level to javascript<####
print "var matchGroup = new Array;\n";
my $curGroup=0;
my $grSize;
foreach my $identifier (sort{$proteinList{$a}[10] <=> $proteinList{$b}[10]} keys %proteinList) {
	next unless $proteinList{$identifier}[0];
	if ($proteinList{$identifier}[10] > $curGroup) {
		$curGroup=$proteinList{$identifier}[10];
		print "matchGroup['$curGroup'] = new Array();\n";
		$grSize=0;
	}
	print "matchGroup['$curGroup'][$grSize]=$proteinList{$identifier}[0];\n";
	$grSize++;
}

####>Bad Confidence<####
print "var badConf = new Array();\n";
foreach my $identifier (keys %proteinList) {
	next unless $proteinList{$identifier}[0];
	print "badConf['$proteinList{$identifier}[0]']=1;\n" if ($proteinList{$identifier}[9] && $proteinList{$identifier}[9]==1);
}

print qq
|function editProteinList() {
	parent.rankFrame.location="$promsPath{cgi}/listMyTempProteins.cgi?ACT=form&ID=$analysisID";
	parent.spectrumFrame.location="$promsPath{html}/nothing.html";
}
function deleteProteinList() {
	if (confirm('Delete protein list?')) {
		parent.rankFrame.location="$promsPath{html}/nothing.html";
		parent.spectrumFrame.location="$promsPath{html}/nothing.html";
		window.location="./listMyTempProteins.cgi?ACT=delete&ID=$analysisID";
	}
}
function openProtein(protId,protGroup,protPage,updateRank) {
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
	var newProt=document.getElementById(protId);
	newProt.style.color='#DD0000'; //style.background = '#DD0000'; //value="Ca marche"; //
	for (var i=0;i<matchGroup[protGroup].length;i++) {
		var newGrProt=parent.queryFrame.document.getElementById(matchGroup[protGroup][i]);
		if (newGrProt && newGrProt != newProt) {
			newGrProt.style.color='#0000DD'; //'#DD7777';
		}
	}
	selectedProtId=protId;
	selectedGroup=protGroup;
	// Update rankFrame
	if (updateRank) {
		parent.selectedProtId=protId;
		parent.selProteinPage=protPage;
		var pepRangeString=(parent.rankFrame.document.rankForm && protId==parent.rankFrame.document.rankForm.ID.value)? parent.rankFrame.pepRangeString : '';
		parent.rankFrame.location="$promsPath{cgi}/listRanks.cgi?ANAID=$analysisID&MSTYPE=$msType&ITEM=myList&SORT=$selectedSort&PAGE=$pageStrg&ID="+protId+"&SHOWFILT=$showFiltered&sortType="+parent.listRankSort+"&Dlt="+parent.listRankDlt+"&maxDisp="+parent.maxDisp+"&maxSrch="+parent.maxSrch+"&varValue="+parent.varValue+pepRangeString+"&multiAna="+parent.multiAnaScan;
	}
}
var selectedProtId=$selectedProtID;
var selectedGroup=$selectedGroup;
</SCRIPT>
</HEAD>
<BODY topmargin="0" leftmargin="0" rightmargin="0" $onloadStrg>
<CENTER>
<INPUT type="button" value="Edit" onclick="editProteinList()" $disableString><INPUT type="button" value="Delete" onclick="deleteProteinList()" $disableString>
|;
my ($color1,$color2)=&promsConfig::getRowColors;
my $bgColor=$color2;
print qq
|<TABLE border=0 cellspacing=0>
<TR bgcolor="$color2"><TH class="rbBorder" width=190>Protein</TH><TH class="bBorder" width=60>Matches</TH></TR>
|;
my $curMatchGr=0;
foreach my $identifier (sort{&sortProt($sortKeyPos{$sortKey},$sortOrder)} keys %proteinList) {
	my ($protID,$ident,$protDes,$protMW,$protLength,$protOrg,$selStatus,$numMatch,$maxMatch,$confLevel,$matchGroup,$score,$maxScore,$dbRank)=@{$proteinList{$identifier}};
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
	if ($sortByMatchGr) {
		if ($curMatchGr!=$matchGroup) {
			$bgColor=($bgColor eq $color1)? $color2 : $color1;
			$curMatchGr=$matchGroup;
		}
	}
	else {$bgColor=($bgColor eq $color1)? $color2 : $color1;}
	print "<TR class=\"list\" bgcolor=\"$bgColor\">";
	print "<TD>&nbsp<IMG src=\"$promsPath{images}/$image\" hspace=0 border=0 height=11 width=11>&nbsp";
	if ($protID) {
		print "<A class=\"$class\" id=\"$protID\" href=\"javascript:openProtein($protID,$matchGroup,$proteinPage{$protID},1)\" onmouseover=\"popup('<B>Protein:</B> $identifier<BR><B>Description:</B> $protDes.<BR><B>Size:</B> $lengthString<BR><B>Mass:</B> $mwString<BR><B>Organism:</B> <I><U>$protOrg</U></I>.')\" onmouseout=\"popout()\">";
		my $identifierStrg=($numDatabanks>1 && $dbRank)? $dbRank.'::'.$identifier : $identifier;
		$identifierStrg=&promsMod::resize($identifierStrg,20);
		if ($selStatus==-3) {print "<I>$identifierStrg</I>";} else {print "$identifierStrg";}
		print "</A></TD>\n<TD align=center>$numMatch/$maxMatch</TD>";
	}
	else {
		print "<FONT class=\"noMatch\"  onmouseover=\"popup('<B>Not found in current Analysis</B>')\" onmouseout=\"popout()\">$identifier</FONT></TD>\n<TD align=center>-</TD>";
	}
	print "</TR>\n";
}
print qq
|</TABLE>
<HR width="90%">
</CENTER>
<B>
&nbsp<U>Validation status :</U><BR>
&nbsp<IMG src='$promsPath{images}/lightGray1.gif' width=11 height=11>&nbsp;: not verified.<BR>
&nbsp<IMG src='$promsPath{images}/lightRed1.gif' width=11 height=11>&nbsp;: no matching peptides.<BR>
&nbsp<IMG src='$promsPath{images}/lightYellow1.gif' width=11 height=11>&nbsp;: partially verified.<BR>
&nbsp<IMG src='$promsPath{images}/lightGreen1.gif' width=11 height=11>&nbsp;: fully verified.<BR>
&nbsp<IMG src='$promsPath{images}/lightBlue1.gif' width=11 height=11>&nbsp;: excluded from list.<BR>
|;
print "&nbsp<IMG src='$promsPath{images}/lightBlue1.gif' width=11 height=11>&nbsp;: <I>filtered out.</I><BR>\n" if $showFiltered;
print qq
|&nbsp;<FONT color='#DD000'>Selected protein.</FONT><BR>
&nbsp;<FONT color='#0000DD'>Same match group.</FONT><BR>
&nbsp;<FONT class="noMatch">Not in current Analysis.</FONT><BR>
&nbsp;<FONT color='#707070'>Poor confidence level.</FONT>
</B>

<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup()</SCRIPT>

</BODY>
</HTML>
|;

############## Protein sort subroutines ###############
sub sortProt {
	my ($sortPos,$sortOrder)=@_;
	# Match group. (MG asc> MaxScore desc > Identifier asc)
	if ($sortPos==10) {$proteinList{$a}[10]<=>$proteinList{$b}[10] || $proteinList{$b}[12]<=>$proteinList{$a}[12] || lc($proteinList{$a}[1]) cmp lc($proteinList{$b}[1])}
	# Distribution asc. (Item asc > Identifier asc)
	elsif ($sortOrder eq 'ASC') {$proteinList{$a}[$sortPos]<=>$proteinList{$b}[$sortPos] || lc($proteinList{$a}[1]) cmp lc($proteinList{$b}[1])}
	# Distribution desc. (Item desc > Identifier asc)
	else {$proteinList{$b}[$sortPos]<=>$proteinList{$a}[$sortPos] || lc($proteinList{$a}[1]) cmp lc($proteinList{$b}[1])}
}


########################################
####<Displaying Form (in rankFrame)>####
########################################
sub displayForm {

	my $subAction=(-e $myProtFile)? 'edit' : 'add';

	####<Fetching list of validated analyses>####
	my @analysisList;
	my $anaQuery=qq |SELECT ID_ANALYSIS,EXPERIMENT.NAME,SAMPLE.NAME,ANALYSIS.NAME
							FROM EXPERIMENT,SAMPLE,ANALYSIS
							WHERE ID_PROJECT=$projectID AND EXPERIMENT.ID_EXPERIMENT=SAMPLE.ID_EXPERIMENT
							AND SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND VALID_STATUS>=1
							ORDER BY EXPERIMENT.DISPLAY_POS ASC,SAMPLE.DISPLAY_POS ASC,ANALYSIS.DISPLAY_POS ASC|;
	my $sthAna=$dbh->prepare($anaQuery);
	my $sthNPA=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=?");
	$sthAna->execute;
	while (my ($anaID,@anaInfo)=$sthAna->fetchrow_array) {
		next if $anaID==$analysisID;
		$sthNPA->execute($anaID);
		my ($numProt)=$sthNPA->fetchrow_array;
		if ($numProt) {
			for (my $i=0; $i<$#anaInfo-1; $i++) {
				$anaInfo[$i]=&promsMod::shortenName($anaInfo[$i],15);
			}
			push @analysisList,[$anaID,join (' > ',@anaInfo)." (x$numProt)"];
		}
	}
	$sthAna->finish;
	$sthNPA->finish;

	####<Fetching list of categories>####
	my @categoryList;
	my $sthCat=$dbh->prepare("SELECT ID_CATEGORY,CATEGORY.NAME,CLASSIFICATION.NAME FROM CLASSIFICATION,CATEGORY WHERE CLASSIFICATION.ID_CLASSIFICATION=CATEGORY.ID_CLASSIFICATION AND ID_PROJECT=$projectID ORDER BY CLASSIFICATION.NAME ASC,CATEGORY.NAME ASC");
	my $sthNPC=$dbh->prepare("SELECT COUNT(*) FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=?");
	$sthCat->execute;
	while (my ($catID,$catName,$className)=$sthCat->fetchrow_array) {
		$sthNPC->execute($catID);
		my ($numProt)=$sthNPC->fetchrow_array;
		if ($numProt) {
			$catName=&promsMod::shortenName($catName,19);
			$className=&promsMod::shortenName($className,19);
			push @categoryList,[$catID,join (' > ',($className,$catName))." (x$numProt)"];
		}
	}
	$sthCat->finish;
	$sthNPC->finish;
	$dbh->disconnect;

	####<Fetching private list of proteins>####
	my @proteinList;
	my $noMatch=0;
	if ($subAction eq 'edit') {
		open (MYPROT,$myProtFile);
		while (<MYPROT>) {
			my ($protIdent,$protID)=($_=~/^(\S+)\t(\d+)/);
			unless ($protID) {
				$protIdent="*$protIdent";
				$noMatch=1;
			}
			push @proteinList,$protIdent;
		}
		close MYPROT;
	}

	#######################
	####<Starting HTML>####
	#######################
	my $title=($subAction eq 'edit')? 'Editing Protein List' : 'Creating Protein List';
	my ($lightColor,$darkColor)=&promsConfig::getRowColors;
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Private Protein List Form</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE>
.checklist {
	border: 1px solid #ccc;
	list-style:none;
	height:240px;
	width:400px;
	overflow:auto;
}
.checklist, .checklist LI {margin:0; padding:0;}
.checklist LABEL {
    padding-left:3px;
    text-indent:-25px;
}
.liRow1 {background:#FFFFFF;}
.liRow2 {background:#F0F0F0;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
function changeImportSource(source) {
	document.getElementById('identListDiv').style.display='none';
	document.getElementById('selAnaDiv').style.display='none';
	document.getElementById('selCatDiv').style.display='none';
	if (source=='identList') {document.getElementById('combDiv').style.display='none';}
	else if ('$subAction'=='edit') {document.getElementById('combDiv').style.display='block';}
	document.getElementById(source+'Div').style.display='block';
}
function removeNotMatching() {
	var identArray=document.myListForm.identList.value.split(/[\\n,]/);
	var cleanArray=new Array();
	for (var i=0; i<identArray.length; i++) {
		if (identArray[i].match(/\\*/)) continue;
		cleanArray.push(identArray[i]);
	}
	document.myListForm.identList.value=cleanArray.join('\\n');
}
function cancelAction() {
	if ('$subAction'=='edit') {parent.writeMenu(parent.selectedProtID,2);}
	else {
		parent.selectedView='protein';
		parent.writeMenu(parent.selectedProtID,1);
	}
}
function checkForm(myForm) {
	//Identifier list
	if (myForm.source.value=='identList' && !myForm.identList.value) {
		alert('Type or paste a list of protein identifiers');
		return false;
	}
	//Analysis list
	else if (myForm.source.value=='selAna') {
		var okCheck=false;
		if (myForm.selAna.length) {
			for (var i=0; i<myForm.selAna.length; i++) {
				if (myForm.selAna[i].checked) {
					okCheck=true;
					break;
				}
			}
		}
		else if (!myForm.selAna.checked) okCheck=true;
		if (!okCheck) {
			alert('Select at least 1 Analysis');
			return false;
		}
	}
	//Category list
	else if (myForm.source.value=='selCat') {
		var okCheck=false;
		if (myForm.selCat.length) {
			for (var i=0; i<myForm.selCat.length; i++) {
				if (myForm.selCat[i].checked) {
					okCheck=true;
					break;
				}
			}
		}
		else if (!myForm.selCat.checked) okCheck=true;
		if (!okCheck) {
			alert('Select at least 1 Category');
			return false;
		}
	}
	return true;
}
</SCRIPT>
</HEAD>
<BODY topmargin="0" background="$promsPath{images}/bgProMS.gif">
<FORM name="myListForm" method="post" onsubmit="return(checkForm(this));" target="spectrumFrame">
<INPUT type="hidden" name="ACT" value="submit">
<INPUT type="hidden" name="ID" value="$analysisID">
<CENTER>
<FONT class="title">$title</FONT><BR>
<TABLE  bgcolor=$darkColor width=550>
<TR><TH align=right valign=top nowrap>&nbsp;Select source :</TH><TD bgcolor=$lightColor><SELECT class="title3" name="source" onchange="changeImportSource(this.value)">
	<OPTION value="identList">Type or paste list of identifiers</OPTION>
|;
	print "\t<OPTION value=\"selAna\">Generate from other Analyses</OPTION>\n" if scalar @analysisList;
	print "\t<OPTION value=\"selCat\">Generate from Classifications</OPTION>\n" if scalar @categoryList;
	print qq
|</SELECT>
</TD></TR>
<TR><TH align=right valign=top bgcolor=$darkColor>Identifiers :</TH><TD bgcolor=$lightColor valign=top>
<TABLE id="identListDiv" style="display:block" cellpadding=0 cellspacing=0><TR><TD><TEXTAREA name="identList" rows="15" cols="25">
|;
	foreach my $identifier (@proteinList) {print "$identifier\n";}
	print "</TEXTAREA></TD>";
	if ($noMatch) {
		print '<TD valign=top>* Not found in current Analysis<BR>&nbsp<INPUT type="button" value="Remove *Proteins" valign=top onclick="removeNotMatching();"></TD>';
	}
	print qq
|</TR></TABLE>
<DIV id="selAnaDiv" style="display:none">
<UL class="checklist">
|;
	my $countAna=0;
	my $liClass='liRow1';
	foreach my $refAna (@analysisList) {
		$countAna++;
		my $boxID="chkAna$countAna";
		print "\t<LI class=\"$liClass\"><INPUT type=\"checkbox\" name=\"selAna\" id=\"$boxID\" value=\"$refAna->[0]\"><LABEL FOR=\"$boxID\">$refAna->[1]</LABEL></LI>\n";
		$liClass=($liClass eq 'liRow1')? 'liRow2' : 'liRow1';
	}
	print "<UL>\n</DIV>\n<DIV id=\"selCatDiv\" style=\"display:none\">\n<UL class=\"checklist\">\n";

	my $countCat=0;
	$liClass='liRow1';
	foreach my $refCat (@categoryList) {
		$countCat++;
		my $boxID="chkCat$countCat";
		print "\t<LI class=\"$liClass\"><INPUT type=\"checkbox\" name=\"selCat\" id=\"$boxID\" value=\"$refCat->[0]\"><LABEL FOR=\"$boxID\">$refCat->[1]</LABEL></LI>\n";
		$liClass=($liClass eq 'liRow1')? 'liRow2' : 'liRow1';
	}
	print qq
|<UL></DIV>
<DIV id="combDiv" style="display:none"><INPUT type="checkbox" name="combine" value=1><B>Add to existing list</B></DIV></TD></TR>
<TR><TH colspan=2><INPUT type="submit" name="submit" value="Proceed">&nbsp;&nbsp;&nbsp<INPUT type="button" value="Cancel" onclick="cancelAction();"></TH></TR>
</TABLE>
</CENTER>
</FORM>
</BODY>
</HTML>
|;
	exit;
}

########################################
####<Processing Form (in rankFrame)>####
########################################
sub processForm {
	my $protSource=param('source');
	my $combine=param('combine');

	####<Starting HTML>####
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Private Protein List Form</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<BR>
<FONT class="title3">
|;
	my %identifierList;
	if ($protSource eq 'identList' || $combine) {
		print "Reading identifier list...";
		foreach my $identifier (split(/[\s+,]/,param('identList'))) {
			$identifier=~s/\*//; # not in analysis
			$identifier=~s/^\d+:://; # multi-databank search
			$identifierList{$identifier}=1 if $identifier=~/\S+/;
		}
		print " Done.<BR><BR>\n";
	}
	if ($protSource eq 'selAna') {
		print "Generating protein list from selected analyses...";
		my $sthAP=$dbh->prepare("SELECT IDENTIFIER FROM PROTEIN,ANALYSIS_PROTEIN WHERE ID_ANALYSIS=? AND PROTEIN.ID_PROTEIN=ANALYSIS_PROTEIN.ID_PROTEIN AND IDENTIFIER NOT LIKE 'DECOY_%'");
		foreach my $anaID (param('selAna')) {
			$sthAP->execute($anaID);
			while (my ($identifier)=$sthAP->fetchrow_array) {
				$identifierList{$identifier}=1;
			}
			print '.';
		}
		$sthAP->finish;
		print " Done.<BR><BR>\n";
	}
	elsif ($protSource eq 'selCat') {
		print "Generating protein list from selected categories...";
		my $sthCP=$dbh->prepare("SELECT IDENTIFIER FROM PROTEIN,CATEGORY_PROTEIN WHERE ID_CATEGORY=? AND PROTEIN.ID_PROTEIN=CATEGORY_PROTEIN.ID_PROTEIN");
		foreach my $catID (param('selCat')) {
			$sthCP->execute($catID);
			while (my ($identifier)=$sthCP->fetchrow_array) {
				$identifierList{$identifier}=1;
			}
			print '.';
		}
		$sthCP->finish;
		print " Done.<BR><BR>\n";
	}

	####<Splitting in pages>####
	print "Processing data...";
	my $sthPG=$dbh->prepare ("SELECT ID_PROT_VALID,IDENTIFIER,MATCH_GROUP FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$analysisID AND IDENTIFIER NOT LIKE 'DECOY_%' ORDER BY MAX_SCORE DESC,ID_PROT_VALID ASC");
	$sthPG->execute || die $sthPG->errstr;
	my $refData=$sthPG->fetchall_arrayref; # reference to an array
	$sthPG->finish;

	$dbh->disconnect;

	my (%matchGroups,%proteinData);
	foreach my $refProt (@{$refData}) {
		$matchGroups{$refProt->[2]}++;
		@{$proteinData{$refProt->[1]}}=($refProt->[0],$refProt->[2]) if $identifierList{$refProt->[1]};
	}
	my $protCount=0;
	my $page=1;
	my $pageSize=&promsConfig::getProteinPageSize;
	my (%usedMGr,%pageMGr);
	foreach my $refProt (@{$refData}) {
		next if $usedMGr{$refProt->[2]};
		$usedMGr{$refProt->[2]}=1;
		$protCount+=$matchGroups{$refProt->[2]};
		$pageMGr{$refProt->[2]}=$page;
		if ($protCount >= $pageSize) {
			$protCount=0;
			$page++;
		}
	}
	print " Done.<BR><BR>\n";

	unless (scalar keys %proteinData) {
		print qq
|<FONT color=#DD0000>No proteins were matched in current Analysis.</FONT>
</BODY>
</HTML>
|;
		exit;
	}
	if (scalar keys %proteinData > 200) {
		print qq
|<FONT color=#DD0000>Too many proteins were matched in current Analysis (>200).</FONT>
</BODY>
</HTML>
|;
		exit;
	}

	####<Writing list to file>####
	print "Storing data...";
	open (MYPROT, ">$myProtFile") || die $!; # overwrites if exists
	foreach my $identifier (sort{$a cmp $b} keys %identifierList) {
		if ($proteinData{$identifier}) {print MYPROT "$identifier\t$proteinData{$identifier}[0]\t$pageMGr{$proteinData{$identifier}[1]}\n";}
		else {print MYPROT "$identifier\t0\t0\n";}
	}
	close MYPROT;
	print " Done.<BR><BR>\n";

	####<Reloading personal list frame>####
	sleep 2;
	print qq
|<SCRIPT LANGUAGE="JavaScript">
parent.setView(parent.selectView);
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

####>Revision history<####
# 1.1.1 Added multiAna as listRanks parameter (PP 26/09/13)
# 1.1.0 GPL license (PP 23/09/13)
# 1.0.9 Multi-databank flag on proteins (PP 11/12/12)
# 1.0.8 Manager & super bio has full access validation data (PP 01/10/11)
