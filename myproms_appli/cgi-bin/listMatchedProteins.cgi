#!/usr/local/bin/perl -w

################################################################################
# listMatchedProteins.cgi         2.0.6                                        #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Lists the proteins matched by a given (set of) peptide(s)                    #
#  during a validation process.                                                #
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

# print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
#############################
####>Connecting to ProMS<####
#############################
my $dbh=&promsConfig::dbConnect;

####################
####>Parameters<####
####################
my (%matchedProt,%protDbRank,%matchMulti,%protDes,%protOrg,%protLength,%mWeight,%selStatus,%matchInfo,%score,%class,%protMGr);
my ($color1,$color2)=&promsConfig::getRowColors;
my $analysisID=param('ANAID');
my $msType=param('MSTYPE');
my $selIdentifier=param('PROT')? param('PROT') : '';
my $item=param('ITEM');
&excludeProt if param('EXCLUDE');
my $showFiltered=param('SHOWFILT');
my $filterString=($showFiltered)? '' : 'SEL_STATUS > -3 AND';
my $disableString=(param('DIS_EXCLUDE'))? param('DIS_EXCLUDE') : '';
if ($item eq 'rank') {&matchRank;}
else {&matchList;} # ITEM='list'
exit;

#################################<<<SUBROUTINES>>>######################################

##################################
####<< Subroutine matchRank >>#### Displays all proteins matched by a single peptide sequence
##################################
sub matchRank {
	my $match=param('MATCH');
	my ($queryID,$queryNum,$rank)=($match=~/prot_(\d+)_(\d+)_(\d+)/);
	my $sequence=param('SEQ');


	#########################################
	####<Finding matched proteins + info>####
	#########################################

	###<DB queries>###
	my ($numDatabanks)=$dbh->selectrow_array("SELECT COUNT(*) FROM ANALYSIS_DATABANK WHERE ID_ANALYSIS=$analysisID");
	my $sth1=$dbh->prepare("SELECT IDENTIFIER,MATCH_MULTI FROM RANK_PROTEIN_MATCH WHERE QUERY_NUM=$queryNum AND PEP_RANK=$rank AND ID_ANALYSIS=$analysisID ORDER BY IDENTIFIER ASC");
	my $sth2=$dbh->prepare("SELECT ID_PROT_VALID,DB_RANK,PROT_DES,ORGANISM,MW,PROT_LENGTH,SEL_STATUS,NUM_MATCH,MAX_MATCH,SCORE,CONF_LEVEL,MATCH_GROUP FROM PROTEIN_VALIDATION WHERE $filterString IDENTIFIER=? AND ID_ANALYSIS=$analysisID");

	###<Finding identifier for matched proteins from RANK_PROTEIN_MATCH table
	my $notAllExcluded=0;
	$sth1->execute();
	while (my ($identifier,$multi) = $sth1->fetchrow_array) {
		##<Fetching relevant data for matched proteins from PROTEIN_VALIDATION table
		$sth2->execute($identifier);
		my ($protID,$dbRank,$des,$organism,$weight,$size,$status,$numMatch,$maxMatch,$sc,$cl,$matchGr)=$sth2->fetchrow_array;
		next unless $protID; # There are non filtered proteins
		$matchedProt{$identifier}=$protID;
		$protDbRank{$identifier}=$dbRank if $numDatabanks>1;
		$matchMulti{$identifier}=$multi;
		$protDes{$identifier}=$des;
		$protOrg{$identifier}=$organism;
		$protLength{$identifier}=$size; $protLength{$identifier}='-' unless $protLength{$identifier};
		$mWeight{$identifier}=$weight;
		$selStatus{$identifier}=$status;
		$matchInfo{$identifier}="$numMatch/$maxMatch";
		$score{$identifier}=$sc;
		$class{$identifier}=($cl && $cl==1)? 'poor' : 'good';
		$notAllExcluded++ if $status>=-1;
		$protMGr{$identifier}=$matchGr;
	}

	$sth1->finish;
	$sth2->finish;
	$dbh->disconnect;


	##############
	####<HTML>####
	##############
	&startHTML;
	print qq
|<INPUT type="hidden" name="MATCH" value="$match">
<INPUT type="hidden" name="SEQ" value="$sequence">
<FONT class="title">Proteins matched by <FONT color="#DD0000">$sequence</FONT></FONT>
<BR><BR>
|;
if ($notAllExcluded) {
	print qq
|<TABLE border=0 width=1000><TR><TD>
<INPUT type="button" value="Check / Uncheck All" style="width:175px;" onclick="check(document.proteinList.id_prot)" $disableString/>&nbsp
<INPUT type="button" value="Exclude from list" style="width:150px;" onclick="excludeProteins(document.proteinList.id_prot)" $disableString/>
</TD></TR></TABLE>
|;
}
	&printTable(\%matchedProt);
	print "</FORM>\n</CENTER>\n";
	print end_html;
}

#####################################
####<< Subroutine matchList >>#### Displays proteins matched by a list of peptide sequences
#####################################
sub matchList {

	#############################
	####<List of SQL queries>####
	#############################
	my $sthSEL1=$dbh->prepare("SELECT QUERY_NUM,PEP_RANK FROM RANK_PROTEIN_MATCH WHERE IDENTIFIER=? AND ID_ANALYSIS=$analysisID");
	my $sthSEL2=$dbh->prepare("SELECT IDENTIFIER FROM RANK_PROTEIN_MATCH WHERE QUERY_NUM=? AND PEP_RANK=? AND ID_ANALYSIS=$analysisID");
	my $sthSEL3=$dbh->prepare("SELECT ID_PROT_VALID,PROT_DES,ORGANISM,MW,PROT_LENGTH,SEL_STATUS,NUM_MATCH,MAX_MATCH,SCORE,CONF_LEVEL,MATCH_GROUP FROM PROTEIN_VALIDATION WHERE $filterString IDENTIFIER=? AND ID_ANALYSIS=$analysisID");

	####################################################################
	####<Finding list of matching peptides and all proteins matched>####
	####################################################################
	my %matchList;
	my %tempProtList;
	my %proteinList;
	my $maxSame=0;
	$sthSEL1->execute($selIdentifier) || die $dbh->errstr;
	while (my ($queryNum,$rank) = $sthSEL1->fetchrow_array) {
		my ($rankInfo)=$dbh->selectrow_array("SELECT INFO_PEP$rank FROM QUERY_VALIDATION WHERE QUERY_NUM=$queryNum AND ID_ANALYSIS=$analysisID");
		($matchList{"$queryNum:$rank"})=($rankInfo=~/SEL=(-*\d)/);
		next if $matchList{"$queryNum:$rank"}==-1; # skip non selectable
		$maxSame++;
		$sthSEL2->execute($queryNum,$rank) || die $dbh->errstr;
		while (my ($identifier) = $sthSEL2->fetchrow_array) {
			$tempProtList{$identifier}=0;
		}
	}
	$tempProtList{$selIdentifier}=$maxSame;
	$sthSEL2->finish();

	################################################################
	####<Finding proteins strickly matched by peptides (sub)set>####
	################################################################
	foreach my $identifier (keys %tempProtList) {
		next if $identifier eq $selIdentifier;
		$sthSEL1->execute($identifier) || die $dbh->errstr;
		my $numSame=0;
		while (my ($queryNum,$rank) = $sthSEL1->fetchrow_array) {
			unless (defined($matchList{"$queryNum:$rank"})) {
				$numSame=0; # protein is matched by a non listed peptide => exclude
				last;
			}
			$numSame++ unless $matchList{"$queryNum:$rank"}==-1; # skip non selectable
		}
		$tempProtList{$identifier}=$numSame;
	}
	$sthSEL1->finish();

	###########################################
	####<Fetching info on matched proteins>####
	###########################################
	my $notAllExcluded;
	foreach my $identifier (keys %tempProtList) {
		next unless $tempProtList{$identifier}; # skip excluded proteins
		$sthSEL3->execute($identifier) || die $dbh->errstr;
		my ($protID,$des,$organism,$weight,$size,$status,$numMatch,$maxMatch,$sc,$confLevel,$matchGr)=$sthSEL3->fetchrow_array;
		next unless $protID; # There are non filtered proteins
		$proteinList{$tempProtList{$identifier}}{$identifier}=1;
		$matchedProt{$identifier}=$protID;
		$protDes{$identifier}=$des;
		$protOrg{$identifier}=$organism;
		$mWeight{$identifier}=$weight;
		$protLength{$identifier}=$size; $protLength{$identifier}='-' unless $protLength{$identifier};
		$selStatus{$identifier}=$status;
		$matchInfo{$identifier}="$numMatch/$maxMatch";
		$score{$identifier}=$sc;
		$class{$identifier}=($confLevel && $confLevel==1)? 'poor' : 'good';
		$notAllExcluded++ if $status>=-1;
		$protMGr{$identifier}=$matchGr;
	}
	$sthSEL3->finish;
	$dbh->disconnect;

	##############
	####<HTML>####
	##############
	&startHTML;
	print qq
|<FONT class="title">Proteins with similar peptide match pattern than <FONT color=#DD0000>$selIdentifier</FONT></FONT><BR>
<FONT class="title2">($maxSame peptides in list)</FONT>
<BR><BR>
|;
if ($notAllExcluded) {
	print qq
|<TABLE border=0 width=935><TR><TD>
<INPUT type="button" value="Check / Uncheck All" style="width:175px;" onclick="check(document.proteinList.id_prot)" $disableString/>&nbsp
<INPUT type="button" value="Exclude from list" style="width:150px;" onclick="excludeProteins(document.proteinList.id_prot)" $disableString/>
</TD></TR></TABLE>
|;
}
	foreach my $numSame (sort {$b<=>$a} keys %proteinList) {
		if ($numSame==$maxSame) {
			print "<FONT class=\"title2\">Proteins with identical match pattern :</FONT><BR>\n";
		}
		else {
			$pepString=($numSame==1)? 'peptide' : 'peptides';
			print "<FONT class=\"title2\">Proteins strickly matched by $numSame $pepString from list :</FONT><BR>\n";
		}
		&printTable(\%{$proteinList{$numSame}});
	}
	print "</FORM>\n</CENTER>\n";
	print end_html;
}


#######################
####<Starting HTML>####
#######################
sub startHTML {
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
A.good:link {color: black} A.good:visited {color:black}
A.poor:link {color: #707070} A.poor:visited {color:#707070}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
|;
&promsMod::popupInfo();
print qq
|function check(field){
	if (!field) {
		alert('There are no selectable proteins in list.');
		return;
	}
	if (checkflag=="false"){
		if (field.length){
			for (i = 0; i < field.length; i++){
				field[i].checked = true;
			}
		}
		else {field.checked = true;}
		checkflag="true";
	}
	else{
		if (field.length){
			for (i = 0; i < field.length; i++){
				field[i].checked = false;
			}
		}
		else {field.checked = false;}
		checkflag="false";
	}
}
function excludeProteins(field){
	if (!field) {
		alert('There are no selectable proteins in list.');
		return;
	}
	var checked=0;
	if (field.length){
		for (i = 0; i < field.length; i++){
			if (field[i].checked == true) {checked++;}
		}
	}
	else{
		if (field.checked == true){checked++;}
	}
	if (checked == 0){
		alert('You must select at least 1 protein');
	}
	else{
		if (confirm ("Exclude selected protein(s) from list of matched proteins ?")) {
			document.proteinList.EXCLUDE.value=1;
			document.proteinList.submit();
		}
	}
}
function findPageAndUpdate(protId,matchGr) {
	var page;
	PAGE:for (var i=0; i<parent.proteinPages.length; i++) {
		var pageBlocks=parent.proteinPages[i].split(/,/);
		for (var j=0; j<pageBlocks.length; j++) {
			var pageRange=pageBlocks[j].split(/-/);
			if (pageRange.length==2) {
				if (matchGr>=pageRange[0] && matchGr<=pageRange[1]) {
					page=i+1;
					break PAGE;
				}
			}
			else if (matchGr==pageRange[0]) {
				page=i+1;
				break PAGE;
			}
		}
	}
	parent.updateValidation(protId,0,0,0,0,'protein',page,1);
}
var checkflag="false";
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">

<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup()</SCRIPT>

<CENTER>
<FORM name="proteinList">
<INPUT type="hidden" name="EXCLUDE">
<INPUT type="hidden" name="ANAID" value="$analysisID">
<INPUT type="hidden" name="MSTYPE" value="$msType">
<INPUT type="hidden" name="PROT" value="$selIdentifier">
<INPUT type="hidden" name="ITEM" value="$item">
|;
}

####################
####<HTML table>####
####################
sub printTable {
	my $refList=$_[0];
	my $bgColor=$color1;
	print qq
|<TABLE border=0 cellpadding=0 cellspacing=0>
	<TR bgcolor="$color2">
	<TH class="rbBorder">Protein</TH>
|;
	print "<TH class=\"rbBorder\" width=65>Match<BR>freq.</TH>\n" if $item eq 'rank';
	my $scoreLabel=($msType eq 'PMF')? 'Max.<BR>Score' : 'Score';
	print qq
|	<TH class="rbBorder" width=65>Valid<BR>matches</TH>
	<TH class="rbBorder" width=60>$scoreLabel</TH>
	<TH class="rbBorder" width=70>Size (aa)</TH>
	<TH class="rbBorder" width=65>Mass (kDa)</TH>
	<TH class="bBorder" align=left width=480>&nbsp&nbsp;Description - Organism</TH>
</TR>
|;
	print "<TR bgcolor=\"$bgColor\"><TH>No valid proteins</TH></TR>\n" unless scalar keys %{$refList};
	foreach my $identifier (sort protSort keys %{$refList}){
	my $image;
		if ($selStatus{$identifier}<=-2) {
			$image='lightBlue1.gif';
		}
		elsif ($selStatus{$identifier}==-1) {
			$image='lightGray1.gif';
		}
		elsif ($selStatus{$identifier}==0) {
			$image='lightRed1.gif';
		}
		elsif ($selStatus{$identifier}==1) {
			$image='lightYellow1.gif';
		}
		elsif ($selStatus{$identifier}==2) {
			$image='lightGreen1.gif';
		}
		my $protString=($protDbRank{$identifier})? $protDbRank{$identifier}.'::'.$identifier : $identifier;
		#my $protString=&promsMod::resize($identifier,16);
		$protString="<FONT color=#DD0000>$protString</FONT>" if ($selIdentifier && $identifier eq $selIdentifier);
		$protString="<I>$protString</I>" if $selStatus{$identifier}==-3;
		print "<TR class=\"list\" bgcolor=\"$bgColor\">\n";
		print "<TH align=left valign=top>&nbsp";
		if ($selStatus{$identifier}>=-1) {
			print "<INPUT type=\"checkbox\" name=\"id_prot\" value=\"$matchedProt{$identifier}\" $disableString>&nbsp\n";
		}
		else {print "&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp\n";}
		print "<IMG src=\"$promsPath{images}/$image\" HSPACE=0 height=11 width=11>&nbsp\n";
		print "<A class=\"$class{$identifier}\" href=\"javascript:findPageAndUpdate($matchedProt{$identifier},$protMGr{$identifier})\" onmouseover=\"popup('<B>$identifier</B><BR>Click to display all matching peptides.')\" onmouseout=\"popout()\">$protString</A>&nbsp;</TH>\n";
		print "<TH valign=top>x$matchMulti{$identifier}</TH>\n" if $item eq 'rank';
		print "<TH valign=top>$matchInfo{$identifier}</TH>\n";
		print "<TH valign=top>$score{$identifier}</TH>\n";
		print "<TH align=right valign=top>$protLength{$identifier}&nbsp</TH>\n";
		if ($mWeight{$identifier}==0) {print "\t<TH align=right valign=top>-&nbsp</TH>\n";}
		else {printf "\t<TH align=right valign=top>%.1f&nbsp</TH>\n",$mWeight{$identifier}/1000;}
		print "\t<TD valign=top>$protDes{$identifier}. <I><U>$protOrg{$identifier}.</U></I></TD>\n";
		print "</TR>\n";
		$bgColor=($bgColor eq $color1)? $color2 : $color1;
	}
	print "</TABLE><BR>\n";
}
sub protSort {
	if ($score{$b}>$score{$a}) {return 1 ;} elsif ($score{$b}<$score{$a}){return -1 ;}
	elsif (($protDes{$b} =~/no\sdescription/) and  ($protDes{$a} !~/no\sdescription/)) {return -1}
	elsif (($protDes{$b} !~/no\sdescription/) and  ($protDes{$a} =~/no\sdescription/)) {return 1}
	elsif  (($protDes{$b} =~/unnamed/) and  ($protDes{$a} !~/unnamed/)) {return -1}
	elsif  (($protDes{$b} !~/unnamed/) and  ($protDes{$a} =~/unnamed/)) {return 1}
}




######################################
####<Exclude protein(s) from list>####
######################################
sub excludeProt {
	my @proteinList=param('id_prot');
	my $sthSEL=$dbh->prepare("SELECT SEL_STATUS FROM PROTEIN_VALIDATION WHERE ID_PROT_VALID=?");
	my $sthUP=$dbh->prepare("UPDATE PROTEIN_VALIDATION SET SEL_STATUS=-2 WHERE ID_PROT_VALID=?");
	my $numSelProt=0;
	foreach my $protID (@proteinList) {
		$sthSEL->execute($protID);
		my ($selStatus)=$sthSEL->fetchrow_array;
		$numSelProt-- if $selStatus > 0;
		$sthUP->execute($protID) || die $dbh->errstr;
	}
	$sthSEL->finish();
	$sthUP->finish();
	$dbh->commit();
	$dbh->disconnect();

	###<Starting HTML>###
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<SCRIPT LANGUAGE="JavaScript">
if (parent.selectedView=='protein') {
	parent.updateValidation(parent.selectedProtId,0,0,0,$numSelProt,parent.selectedView,parent.selProteinPage,1);
}
else { // view=query
	parent.updateValidation(parent.selectedQueryId,0,0,0,$numSelProt,parent.selectedView,parent.selQueryPage,1);
}
</SCRIPT>
</HEAD>
</HTML>
|;
exit;
}

####>Revision history<####
# 2.0.6 GPL license (PP 23/09/13)
# 2.0.5 Multi-databank flag on proteins (PP 11/12/12)
