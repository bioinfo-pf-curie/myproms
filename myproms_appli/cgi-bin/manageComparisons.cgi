#!/usr/local/bin/perl -w

################################################################################
# manageComparisons.cgi    1.0.4                                               #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Displays & edits Comparisons                                                 #
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
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use strict;

# print header; warningsToBrowser(1); #DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my $date=strftime("%Y-%m-%d %H:%M:%S", localtime);

####################
#    Parameters    #
####################
my $action=(param('ACT'))? param('ACT') : 'list';
my $projectID=param('id_project');

####>Connect to the database<####
my $dbh=&promsConfig::dbConnect;


if ($action eq 'edit') {&editComparison; exit;}
elsif ($action eq 'update') {&updateComparison; exit;}
elsif ($action eq 'delete') {&deleteComparison; exit;}


###############
#     Main    #
###############

###>Project PTMs
#my %allPostTransModifs=&promsConfig::getVariableModifications;
#my ($projectPtmString)=$dbh->selectrow_array("SELECT RELEVANT_PTMS FROM PROJECT WHERE ID_PROJECT=$projectID");
#my %projectVarMods;
#if ($projectPtmString) {
#	foreach my $varMod (split(';',$projectPtmString)) {$projectVarMods{$varMod}=$allPostTransModifs{$varMod};}
#}
my %allPostTransModifs=&promsMod::getVariableModifications($dbh);
my %projectVarMods;
my $sthGetPM=$dbh->prepare("SELECT ID_MODIFICATION FROM PROJECT_MODIFICATION WHERE ID_PROJECT=$projectID");
$sthGetPM->execute;
while (my ($modID)=$sthGetPM->fetchrow_array) {
	$projectVarMods{$modID}=$allPostTransModifs{$modID};
	$projectVarMods{$allPostTransModifs{$modID}[0]}=$allPostTransModifs{$modID}; # PSI name for compatibility with old string-based PTM managment
}
$sthGetPM->finish;

######################################
####>Fetching list of Comparisons<####
######################################
my $sthComp=$dbh->prepare("SELECT ID_COMPARISON,ID_CATEGORY,NAME,COMMENTS,NUM_GROUPS,CAT_EXCLUSION,PEPTIDE_PARAMS,AUTOCHECK_PARAMS,UPDATE_DATE,UPDATE_USER FROM COMPARISON WHERE ID_PROJECT=$projectID");
my $sthMod=$dbh->prepare("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES FROM MODIFICATION WHERE ID_MODIFICATION=?");
my $sthCat=$dbh->prepare("SELECT CONCAT(CL.NAME,' > ',CA.NAME) FROM CATEGORY CA,CLASSIFICATION CL WHERE CA.ID_CLASSIFICATION=CL.ID_CLASSIFICATION AND ID_CATEGORY=?");
my $sthUser=$dbh->prepare("SELECT USER_NAME FROM USER_LIST WHERE ID_USER=?");
my $sthAC=$dbh->prepare("SELECT ANA_POS,COMP_GROUP,ID_ANALYSIS FROM ANA_COMPARISON WHERE ID_COMPARISON=? ORDER BY COMP_GROUP ASC,ANA_POS ASC");
my $queryAS=qq |SELECT CONCAT(E.NAME,' > ',S.NAME),A.NAME,A.VALID_STATUS
					FROM EXPERIMENT E,SAMPLE S,ANALYSIS A
					WHERE E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SAMPLE=A.ID_SAMPLE AND S.ID_SPOT IS NULL AND A.ID_ANALYSIS=?|;
my $queryAG=qq |SELECT CONCAT(E.NAME,' > ',G.NAME,' > ',SP.NAME),A.NAME,A.VALID_STATUS
					FROM EXPERIMENT E,GEL2D G,SPOT SP,SAMPLE S,ANALYSIS A
					WHERE E.ID_EXPERIMENT=G.ID_EXPERIMENT AND G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT
					AND S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS=?|;

my @sthList=($dbh->prepare($queryAS),$dbh->prepare($queryAG));
my $sthLC=$dbh->prepare("SELECT CAT_POS,COMP_GROUP,CC.ID_COMPARISON,CL.NAME,CA.NAME FROM CAT_COMPARISON CC,CATEGORY CA,CLASSIFICATION CL WHERE CC.ID_CATEGORY=CA.ID_CATEGORY AND CA.ID_CLASSIFICATION=CL.ID_CLASSIFICATION AND ID_COMPARISON=?");

my %imageIcones=&promsConfig::getItemIcones;
my (%comparisonList,%pepModifications,%categoryList,%itemPath);
$sthComp->execute;
while (my ($compID,$catID,@compInfo)=$sthComp->fetchrow_array) {
	##>Peptide modification filters
	if ($compInfo[4]=~/modifFilters=([^;]+)/) {
		my $modFilterStrg=$1;
		foreach my $filterStrg (split',',$modFilterStrg) {
			$filterStrg=~s/.+://;
			foreach my $modID (split(/\./,$filterStrg)) {
				if ($allPostTransModifs{$modID}) {$pepModifications{$modID}=$allPostTransModifs{$modID}[0];}
				else {
					$sthMod->execute($modID);
					my ($psi,$inter,$syn)=$sthMod->fetchrow_array;
					$pepModifications{$modID}=$psi || $inter;
					unless ($pepModifications{$modID}) {
						$syn=~s/^##//; $syn=~s/##$//; $syn=~s/##/, /;
						$pepModifications{$modID}=$syn;
					}
					$pepModifications{$modID}="&lt;$modID?&gt;" unless $pepModifications{$modID}; # in case modif has been deleted/merged
				}
			}
		}
	}
	##>Category Filter
	if ($catID) {
		unless ($categoryList{$catID}) {
			$sthCat->execute($catID);
			($categoryList{$catID})=$sthCat->fetchrow_array;
		}
	}
	$sthUser->execute($compInfo[-1]); # login -> name
	($compInfo[-1])=$sthUser->fetchrow_array;
	@{$comparisonList{$compID}{'INFO'}}=($catID,@compInfo);

	##>Groups & Analyses/Lists
	#>Analyses
	$sthAC->execute($compID);
	while (my ($anaPos,$group,$anaID)=$sthAC->fetchrow_array) {
		$comparisonList{$compID}{'GROUP'}{$group}[--$anaPos]="A:$anaID";
		unless ($itemPath{"A:$anaID"}) {
			foreach my $sthAI (@sthList) {
				$sthAI->execute($anaID);
				my ($anaPath,$anaName,$valStatus)=$sthAI->fetchrow_array;
				if ($anaPath) {
					my $anaCode=($valStatus==1)? 'analysis:part_val' : 'analysis:val';
					$itemPath{"A:$anaID"}="$anaPath > <IMG src=\"$promsPath{images}/$imageIcones{$anaCode}\">&nbsp;$anaName";
					last;
				}
			}
		}
	}
	##<Lists>##
	$sthLC->execute($compID);
	while (my ($listPos,$group,$listID,$themeName,$listName)=$sthLC->fetchrow_array) {
		$comparisonList{$compID}{'GROUP'}{$group}[--$listPos]="C:$listID";
		$itemPath{"C:$listID"}="$themeName > <IMG src=\"$promsPath{images}/$imageIcones{category}\">&nbsp;$listName";
	}

}
$sthComp->finish;
$sthMod->finish;
$sthCat->finish;
$sthUser->finish;
$sthAC->finish;
$sthList[0]->finish;
$sthList[1]->finish;
$sthLC->finish;

$dbh->disconnect;


#######################
####>Starting HTML<####
#######################
my ($light,$dark)=&promsConfig::getRowColors;
print header(-'content-encoding'=>'no',-charset=>'utf-8');
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>List of Comparisons</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TD {font-size:13px;}
</STYLE>
<SCRIPT LANGUAGE=\"JavaScript\">
function deleteComparison(compID,name) {
	if (confirm("Delete Comparison: '"+name+"'?")) {
		window.location='./manageComparisons.cgi?ACT=delete&id_project=$projectID&ID='+compID;
	}
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">List of Comparisons</FONT>
<BR><BR>
<TABLE border=0 cellspacing=0>
|;
my %pepRulesNoGroups=('ba'=>'in Analysis','banr'=>'distinct in Analysis');
my %pepRulesGroups=('all'=>'in Group','anr'=>'distinct in Group','ba'=>'in best Analysis','banr'=>'distinct in best Analysis');
my %protHidRules=(1=>'hidden everywhere',2=>'where hidden'); # 0 means "Ignore proteins ..." not checked (noHidden=0)
####>Looping through comparisons<####
my $bgColor=$light;
foreach my $compID (sort{lc($comparisonList{$a}{'INFO'}[1]) cmp lc($comparisonList{$b}{'INFO'}[1])} keys %comparisonList) {
	my ($catID,$name,$comments,$numGroups,$catExclusion,$pepParams,$autocheckParams,$updateDate,$updateUser)=@{$comparisonList{$compID}{'INFO'}};
	$comments=&promsMod::HTMLcompatible($comments);
	$comments=~s/<BR>/<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/g;
	$updateDate=&promsMod::formatDate($updateDate);
	my %pepParam;
	foreach my $paramData (split(';',$pepParams)) {
		my ($pName,$pValue)=split('=',$paramData);
		$pepParam{$pName}=$pValue;
	}
	$pepParam{'pepThreshold'}=1 unless $pepParam{'pepThreshold'};
	my $virtualPepString=($pepParam{'virtualData'})? 'Yes' : 'No';
	my $protPepRuleString="&ge;$pepParam{pepThreshold} ";
	$protPepRuleString.=($numGroups)? $pepRulesGroups{$pepParam{'pepRule'}} : $pepRulesNoGroups{$pepParam{'pepRule'}};
	#$protPepRuleString.=' (Including peptides from XIC extractions)' if $pepParam{'virtualData'};
	my $protHidingStrg=($pepParam{'noHidden'}==1)? "Ignore proteins $protHidRules{$pepParam{hiddenRule}}" : 'Include visible and hidden';
	my $catFilterStrg;
	if ($catID) {
		$catFilterStrg=($catExclusion)? "Proteins <B>not</B> in List '$categoryList{$catID}'" : "Restricted to List '$categoryList{$catID}'";
	}
	else {$catFilterStrg='None';}
	my $missCutString=($pepParam{'noMissCut'})? 'Yes' : 'No';
	my $delocPhosphoString=($pepParam{'delocPhospho'})? 'Yes' : 'No';
	my %pepModifFilter=('exclude'=>'None','restrict'=>'None','ignore'=>'None');
	if ($pepParam{'modifFilters'}) {
		foreach my $filterStrg (split(',',$pepParam{'modifFilters'})) {
			my ($filter,$filterValue,$logic)=split(':',$filterStrg); # logic defined only for restrict (and/or)
			my $joinStrg=($filter ne 'restrict')? ', ' : (!$logic || $logic eq 'and')? ' <B>and</B> ' : ' <B>or</B> ';
			$pepModifFilter{$filter}='';
			foreach my $modID (split(/\./,$filterValue)) {
				$pepModifFilter{$filter}.=$joinStrg if $pepModifFilter{$filter};
				$pepModifFilter{$filter}.=$pepModifications{$modID};
			}
		}
	}

	my ($filterLogic,$pepFilterStrg,$selectedRefPTM,$ptmFilterStrg)=split('#',$autocheckParams); # $selectedRefPTM is PSI name (old version) or ptmID
	my $ptmSelectionStrg=(!$selectedRefPTM)? 'None' : ($projectVarMods{$selectedRefPTM})? $projectVarMods{$selectedRefPTM}[0] : "<FONT style=\"font-weight:bold;color:#DD0000\">$allPostTransModifs{$selectedRefPTM}[0]</FONT> <FONT class=\"font11\">(Warning: Not a project relevant PTM)</FONT>";
	my $groupStrg=($numGroups)? "List of compared Groups" : 'List of compared Analyses';
	print qq
|<TR bgcolor=$bgColor>
<TD rowspan=3>&nbsp&nbsp</TD>
<TD width=600><FONT style="font-size:18px;font-weight:bold;">$name</FONT>
<BR><B>&nbsp;&nbsp;&nbsp;Comments:</B> $comments
<BR><B>&nbsp;&nbsp;&nbsp;Include peptides rescued from XIC extractions:</B> $virtualPepString.
<BR><B>&nbsp;&nbsp;&nbsp;&bull;<U>Protein comparison rules:</U></B>
<BR><B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Peptide count rule:</B> $protPepRuleString.
<BR><B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Protein count rule:</B> $protHidingStrg.
<BR><B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Protein restriction:</B> $catFilterStrg.
<BR><B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;PTM selection:</B> $ptmSelectionStrg.

<BR><B>&nbsp;&nbsp;&nbsp;&bull;<U>Peptide comparison rules:</U></B>
<BR><B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Exclude missed cleavages:</B> $missCutString.
<BR><B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Restrict to PTM(s):</B> $pepModifFilter{restrict}.
<BR><B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Excluded PTM(s):</B> $pepModifFilter{exclude}.
<BR><B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Ignored PTM(s):</B> $pepModifFilter{ignore}.
<BR><B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Phosphorylation sites are delocalozed:</B> $delocPhosphoString.

<BR><B>&nbsp;&nbsp;&nbsp;&bull;<U>$groupStrg:</U></B>
</TD>
<TH width=100 rowspan=3><INPUT type="button" value="Edit" style="width:75px" onclick="window.location='./manageComparisons.cgi?ACT=edit&id_project=$projectID&ID=$compID'"><BR>
<INPUT type="button" value="Delete" style="width:75px" onclick="deleteComparison($compID,'$name')">
</TH>
</TR>
<TR bgcolor=$bgColor><TD>
|;
	my @pepFilters=split(';',$pepFilterStrg);
	my @ptmFilters;
	if ($ptmFilterStrg) {@ptmFilters=split(';',$ptmFilterStrg);}
	if ($numGroups) {
		foreach my $group (1..$numGroups) {
			#print "<HR width=90%>\n" if $group > 1;
			print "<TABLE border=0 cellpadding=0 cellspacing=0>\n<TR valign=top><TH align=left valign=top>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;#$group:&nbsp;</TH><TD nowrap>";
			foreach my $anaID (@{$comparisonList{$compID}{'GROUP'}{$group}}) {
				print "-$itemPath{$anaID}&nbsp;<BR>\n";
			}
			print "</TD><TD nowrap>";
			my $existFilter=0;
			my $filterString='&nbsp;&nbsp;&nbsp;<B>Filter:</B> ';
			my ($aGr,$rel,$numPep)=split(':',$pepFilters[$group-1]);
			if (($rel eq '>=' && $numPep > 0) || $rel eq '<=') {
				$existFilter=1;
				$filterString.="Peptides $rel $numPep";
			}
			if ($ptmFilterStrg) {
				my ($aGr,$rel,$numPtm)=split(':',$pepFilters[$group-1]);
				if (($rel eq '>=' && $numPtm > 0) || $rel eq '<=') {
					$filterString.=" & " if $existFilter;
					$existFilter=1;
					$filterString.="PTM $rel $numPtm";
				}
			}
			unless ($existFilter) {
				$filterString.='none';
			}
			print "$filterString&nbsp;</TD>\n</TR></TABLE>\n";
			print "<HR width=90%>\n";
		}
	}
	else { # no groups
		my $anaGroup=1;
		foreach my $anaID (@{$comparisonList{$compID}{'GROUP'}{0}}) {
			#print "<HR width=90%>\n" if $anaGroup > 1;
			#print "<TABLE border=0 cellpadding=0 cellspacing=0>\n<TR><TH align=left>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;#$anaGroup:&nbsp;</TH><TD>$itemPath{$anaID}</TD></TR></TABLE>\n";
			print "<B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;#$anaGroup:&nbsp;</B>$itemPath{$anaID}\n";
			my $existFilter=0;
			my $filterString='&nbsp;&nbsp;&nbsp;<B>Filter:</B> ';
			my ($aGr,$rel,$numPep)=split(':',$pepFilters[$anaGroup-1]);
			if (($rel eq '>=' && $numPep > 0) || $rel eq '<=') {
				$existFilter=1;
				$filterString.="Peptides $rel $numPep";
			}
			if ($ptmFilterStrg) {
				my ($aGr,$rel,$numPtm)=split(':',$pepFilters[$anaGroup-1]);
				if (($rel eq '>=' && $numPtm > 0) || $rel eq '<=') {
					$filterString.=" & " if $existFilter;
					$existFilter=1;
					$filterString.="PTM $rel $numPtm";
				}
			}
			unless ($existFilter) {
				$filterString.='none';
			}
			print "$filterString&nbsp;";
			print "<HR width=90%>\n";
			$anaGroup++;
		}
	}
	print "</TD></TR>\n";
	print "<TR bgcolor=$bgColor><TD colspan><B>&nbsp;&nbsp;&nbsp;Last modified:</B> $updateDate by $updateUser.</TD></TR>\n";
	$bgColor=($bgColor eq $light)? $dark : $light;
}
if (scalar keys %comparisonList == 0) {
	print "<TR><TH colspan=3>(No Comparisons recorded)</TH></TR>\n";
}
print qq
|</TABLE>
</CENTER>
</BODY>
</HTML>
|;


##############################
####<Editing a Comparison>####
##############################
sub editComparison {
	my $comparisonID=param('ID');
	my ($name,$comments)=$dbh->selectrow_array("SELECT NAME,COMMENTS FROM COMPARISON WHERE ID_COMPARISON=$comparisonID");
	$name=~s/"/&quot;/g;
	$comments='' unless $comments;
	$dbh->disconnect;

	print header; warningsToBrowser(1); # DEBUG
	my ($light,$dark)=&promsConfig::getRowColors;
	print qq
|<HTML>
<TITLE>Edit Comparison</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
function checkForm(myForm) {
	if (!myForm.name.value) {
		alert('ERROR: Provide a name for this Comparison.');
		return false;
	}
}
</SCRIPT>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Editing Comparison</FONT>
<BR>
<FORM name="editForm" method="post" onsubmit="return checkForm(this);">
<INPUT type="hidden" name="ACT" value="update">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="ID" value="$comparisonID">
<TABLE border=0><TR><TD bgcolor=$dark><TABLE border=0 cellpadding=2 width=100%>
	<TR><TH width=150 align=right>Name :</TH><TD bgcolor=$light><INPUT type="text" name="name" value="$name" size=65/></TD></TR>
	<TR><TH width=150 align=right valign=top>Comments :</TH><TD bgcolor=$light><TEXTAREA name="comments" rows=3 cols=65/>$comments</TEXTAREA></TD></TR>
	<TR><TH colspan=2><INPUT type="submit" name="save" value=" Save " />&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="window.location='./manageComparisons.cgi?id_project=$projectID'"></TH></TR>
</TD></TR></TABLE>
</BODY>
</HTML>
|;
	exit;
}

###############################
####<Updating a Comparison>####
###############################
sub updateComparison {
	my $comparisonID=param('ID');
	my $name=param('name');
	my $quotedName=$dbh->quote($name);
	#my $comments=(param('comments'))? param('comments') : 'NULL';
	my $quotedComments=$dbh->quote(param('comments'));
	$dbh->do("UPDATE COMPARISON SET NAME=$quotedName,COMMENTS=$quotedComments WHERE ID_COMPARISON=$comparisonID");
	$dbh->commit;
	$dbh->disconnect;

	print header; warningsToBrowser(1); # DEBUG
	print qq
|<HTML>
<SCRIPT LANGUAGE="JavaScript">
if (opener.selCompName && opener.projectID==$projectID) { // on-going comparison for same project => update comp name in list of options
	var compSEL=opener.document.compForm.comparisonID;
	for (var i=0; i<compSEL.options.length; i++) {
		if (compSEL.options[i].value==$comparisonID) {
			if (i==opener.selCompIndex) { // comp is selected
				opener.selCompName=$quotedName;
				var myStoreForm=opener.parent.listFrame.document.storeForm;
				myStoreForm.targetComp.options[1].text=$quotedName; // update target comp option name
				myStoreForm.compComments.value=$quotedComments;
			}
			compSEL.options[i].text=$quotedName;
			break;
		}
	}
}
window.location="./manageComparisons.cgi?id_project=$projectID";
</SCRIPT>
</HTML>
|;
	exit;
}

###############################
####<Deleting a Comparison>####
###############################
sub deleteComparison {
	my $comparisonID=param('ID');

	$dbh->do("DELETE FROM ANA_COMPARISON WHERE ID_COMPARISON=$comparisonID");
	$dbh->do("DELETE FROM CAT_COMPARISON WHERE ID_COMPARISON=$comparisonID");
	$dbh->do("DELETE FROM COMPARISON WHERE ID_COMPARISON=$comparisonID");

	$dbh->commit;
	$dbh->disconnect;

	print header; warningsToBrowser(1); # DEBUG
	print qq
|<HTML>
<SCRIPT LANGUAGE="JavaScript">
if (opener.selCompName && opener.projectID==$projectID) { // on-going comparison for same project => remove comp from list of options
	var compSEL=opener.document.compForm.comparisonID;
	for (var i=0; i<compSEL.options.length; i++) {
		if (compSEL.options[i].value==$comparisonID) {
			if (i==opener.selCompIndex) { // comp is selected
				compSEL.selectedIndex=0;
				opener.selCompIndex=0;
				opener.selCompName='New';
				// remove comp as target in listFrame
				var myListFrame=opener.parent.listFrame;
				myListFrame.displayTargetComparison(0);
				var tgtCompSEL=myListFrame.document.storeForm.targetComp;
				tgtCompSEL.removeChild(tgtCompSEL.options[1]);
			}
			compSEL.removeChild(compSEL.options[i]);
			break;
		}
	}
}
window.location="./manageComparisons.cgi?id_project=$projectID";
</SCRIPT>
</HTML>
|;
	exit;
}

####>Revision history<####
# 1.0.4 Update for '1 vs 1 peptide' comparison filters (PP 18/08/16)
# 1.0.3 Peptide threshold & list exclusion managment (PP 25/04/14)
# 1.0.2 Update for List comparison (PP 21/01/14)
# 1.0.1 Update with getVariableModifications from promsMod (GA 07/05/2013)
# 1.0.0 New script for listing, partial edition & deletion of comparisons (PP 09/03/2011)
