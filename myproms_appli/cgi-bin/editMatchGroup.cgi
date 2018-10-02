#!/usr/local/bin/perl -w

#############################################################################
# editMatchGroup.cgi         2.1.1                                        # 
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                   #
# Contact: myproms@curie.fr                                                 #
# Edits a validated Match Group of proteins                                 #
#############################################################################
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

#print header; warningsToBrowser(1); #DEBUG
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
my $itemID=param('id_item');
my $item=param('item');
my $projectID=&promsMod::getProjectID($dbh,$itemID,$item);
my ($totMatchGr,$modMatchGr);
($totMatchGr,$modMatchGr)=&processForm if param('save');
my $analysisID=param('id_analysis');
my $group=param('group');

#############################
####>Fetching data in DB<####
#############################
my @anaInfo=&promsMod::getItemInfo($dbh,'ANALYSIS',$analysisID);
#foreach my $refItem (@anaInfo) {
#	foreach my $keyType (keys %{$refItem}) {print ">$keyType=>$refItem->{$keyType}<BR>\n";}
#	print "*<BR>\n";
#}
#exit;
my %listProteins; #,%listClusters;
my $bestProtein;
my $sthGr=$dbh->prepare("SELECT PROTEIN.ID_PROTEIN,ALIAS,PROT_DES,ORGANISM,NUM_PEP,NUM_MATCH,SCORE,VISIBILITY,CONF_LEVEL,PROT_LENGTH FROM PROTEIN,ANALYSIS_PROTEIN WHERE PROTEIN.ID_PROTEIN=ANALYSIS_PROTEIN.ID_PROTEIN AND ID_ANALYSIS=$analysisID AND MATCH_GROUP=$group");
$sthGr->execute;
while (my ($protID,@info)=$sthGr->fetchrow_array) {
	@{$listProteins{$protID}}=@info;
	$bestProtein=$listProteins{$protID}[0] if $listProteins{$protID}[6]==2;
}
$sthGr->finish;

$dbh->disconnect;

#######################
####>Starting HTML<####
#######################
my ($color1,$color2)=&promsConfig::getRowColors;
print header;
warningsToBrowser(1); # writes PERL warnings as HTML comments <!-- --> must be activated after header
print qq
|<HEAD>
<TITLE>Edit Match Group</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
ACRONYM {cursor:help;} /*{border-bottom-width:1px;border-bottom-style:dashed;}*/
TH.big {font-size:15px;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
|;
&promsMod::popupInfo();
print qq
|function sequenceView(idProt){
	var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana=$analysisID&id_prot="+idProt+"&msdata="+top.showMSdata;
	top.openProtWindow(winLocation);
}
function clusterView(id_cluster){
	top.promsFrame.optionFrame.autoSelectButton('clusters'); //change button selection in optionFrame
	window.location="./listClusProteins.cgi?id_cluster="+id_cluster+"&ID=$itemID&ITEM=$item&id_project=$projectID";
}
function graphicalView() {
	top.openProtWindow('');
	document.groupForm.largWindow.value=(document.body)? top.protWindow.document.body.clientWidth : top.protWindow.innerWidth;
	document.groupForm.graphView.value=top.graphView;
	document.groupForm.action="./graphicalView.cgi";
	document.groupForm.target='ProteinWindow';
	document.groupForm.submit();
}
function checkBest(box) {
	document.getElementById('extendDiv').style.display='block';
	if (document.getElementById('result')) {document.getElementById('result').style.display='none';}
	if (box.checked==false) { // box was just unchecked
		box.checked=true;
		//alert('You must select 1 best matched protein.');
		return;
	}

	// A new check_best box has just been checked: uncheck old best matched protein
	var bestField=document.groupForm.check_best;
	var visField=document.groupForm.check_vis;
	for (i = 0; i < bestField.length; i++){
		if (bestField[i].value == box.value) {
			visField[i].checked = true;
			continue;
		}
		if (bestField[i].checked == true) {
			visField[i].checked = false; // old best matched: uncheck corresponding check_vis box
		}
		bestField[i].checked = false; // uncheck all other boxes
	}
}
function checkVisible(box) {
	document.getElementById('extendDiv').style.display='block';
	if (document.getElementById('result')) {document.getElementById('result').style.display='none';}
	if (box.checked==true) {return;} // box was just checked

	// Box was just unchecked: make sure it is NOT best matched protein
	var bestField=document.groupForm.check_best;
	for (i = 0; i < bestField.length; i++){
		if (bestField[i].checked == true && bestField[i].value == box.value) {
			box.checked=true;
			//alert('Best matched protein must remain visible.');
			break;
		}
	}
}
function setRule(b1,b2) {
	box1=document.getElementsByName(b1);
	box2=document.getElementsByName(b2);
	if (box1[0].checked==true) {box2[0].checked=false;}
}
</SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
|;
print '<INPUT type="button" value="<< Back" onClick="history.back()">',"\n" unless param('save');
print qq
|<CENTER>
<FONT class="title">Edit Match Group <FONT color="#DD0000">$bestProtein</FONT></FONT>
<BR><BR>
<FONT class="title2">1. Analysis <FONT color="#DD0000">$anaInfo[-1]{NAME}</FONT></FONT>
<BR>
<INPUT type="button" value="Peptide Distribution" onclick="javascript:graphicalView()\">
<BR>

<FORM name="groupForm" method="post">
<INPUT type="hidden" name="id_item" value=$itemID>
<INPUT type="hidden" name="item" value=$item>
<INPUT type="hidden" name="id_analysis" value=$analysisID>
<INPUT type="hidden" name="group" value=$group>
<INPUT type="hidden" name="largWindow">
<INPUT type="hidden" name="what" value="group">
<INPUT type="hidden" name="groupInfo" value="$analysisID:$group">
<INPUT type="hidden" name="graphView">

<TABLE border=0 cellspacing=0>
<TR bgcolor=$color2 valign=bottom>
<TH class="rbBorder" width=50><ACRONYM onmouseover="popup('A protein must be set as the <B>unique representative</B> of the <B>Match Group</B>.')" onmouseout="popout()">Alias</ACRONYM></TH>
<TH class="rbBorder" width=55>Visible</TH>
<TH class="rbBorder" width=200>Protein</TH>
<!-- <TH class="rbBorder" width=100>Cluster</TH> -->
<TH class="rbBorder" width=65>Peptides</TH>
<TH class="rbBorder" width=55>Score</TH>
<TH class="bBorder" align=left width=500>&nbsp;Description - Organism</TH>
</TR>
|;
# exit;
$bgColor=$color1;
# sort : numPep > score > length > id
foreach my $protID (sort {$listProteins{$b}[3]<=>$listProteins{$a}[3] || $listProteins{$b}[5]<=>$listProteins{$a}[5] || $listProteins{$b}[8]<=>$listProteins{$a}[8] || $a<=>$b} keys %listProteins){ # numPep>score>length>id
	print "<TR class=\"list\" bgcolor=\"$bgColor\" valign=top>\n";
	# Best
	my $chkBestString=($listProteins{$protID}[6]==2)? 'checked' : '';
	print "<TH><INPUT type=\"checkbox\" name=\"check_best\" value=\"$protID\" $chkBestString onClick=\"checkBest(this)\"></TH>\n";
	# Visible
	my $chkVisString=($listProteins{$protID}[6]>=1)? 'checked' : '';
	print "<TH><INPUT type=\"checkbox\" name=\"check_vis\" value=\"$protID\" $chkVisString onClick=\"checkVisible(this)\"></TH>\n";
	# Protein
	my $fWeight=($listProteins{$protID}[6]>=1)? 'bold' : 'normal';
	print "<TH align=left><A style=\"font-weight:$fWeight\" href=\"javascript:sequenceView($protID)\"";
	if ($listProteins{$protID}[7]==1){ # bad confidence
		print " onmouseover=\"popup('Massists have assigned a <B>bad confidence level</B> to this Protein.')\" onmouseout=\"popout()\"><FONT color=\"gray\">$listProteins{$protID}[0]</FONT>";
	}
	else{print ">$listProteins{$protID}[0]";}
	print "</A></TH>\n";
	# Peptides (index 3 = num peptides, index 4 = num matches)
	print "<TH>$listProteins{$protID}[3]";
	print "&nbsp<ACRONYM onmouseover=\"popup('$listProteins{$protID}[4] matches due to sequence repeats')\" onmouseout=\"popout()\">($listProteins{$protID}[4])</ACRONYM>" if $listProteins{$protID}[4]>$listProteins{$protID}[3];
	print "</TH>\n";
	# Score
	print "<TH><ACRONYM onmouseover=\"popup('<B>Score for $listProteins{$protID}[0]: $listProteins{$protID}[5]</B>')\" onmouseout=\"popout()\">";
	printf "%.0f</ACRONYM></TH>\n",$listProteins{$protID}[5];
# 	printf "<TH width=59>%.0f</TH>\n",$listProteins{$protID}[5];
	# Description - Organism
	print "<TD>$listProteins{$protID}[1] <FONT style=\"font-style:italic;text-decoration:underline;\">$listProteins{$protID}[2]</FONT></TD>\n";
	print "</TR>\n";
	$bgColor=($bgColor eq $color1)? $color2 : $color1;
}
print "</TABLE>\n";

####>Display results if form was just processed<####
if (param('save')) {
	my $str1=($modMatchGr>1)? 'Groups' : 'Group';
	my $str2=($modMatchGr>1)? 'were' : 'was';
	print "<FONT class=\"title2\" id=\"result\"><BR>Settings were saved ($modMatchGr of $totMatchGr Match $str1 $str2 modified).</FONT>\n";
}

####>Form to extend setting to other Match Groups<####
my $extendVis=(param('save'))? 'none' : 'block';
print qq
|<DIV id="extendDiv" style="display:$extendVis">
<BR>
<FONT class="title2">2. Extend proteins properties to all Analyses in</FONT>&nbsp
<SELECT name="parentItem" style="font-size:16px;font-weight:bold;">
|;
for (my $i=$#anaInfo;$i>=1;$i--) {
	#my $itemType=($i==2)? 'Sample' : ($i==4)? 'Experiment' : 'Project';
	print "<OPTION value=\"$anaInfo[$i]{ITEM}:$anaInfo[$i]{ID}\">$anaInfo[$i]{TYPE} '$anaInfo[$i]{NAME}'</OPTION>\n";
}
print qq
|</SELECT><BR>
<FONT style="font-weight:bold;font-style:italic">(Leave all boxes unchecked to restrict settings to current Match Group)</FONT><BR>
<TABLE cellspacing=0 width=800>
<TR bgcolor="$color2">
<TH class="big" align=left colspan=2>&nbsp&nbsp&nbsp&nbsp&nbsp;Rule for <U>Alias</U> protein:
<FONT style="font-size:13px;font-style:italic">&nbsp(Protein will follow <U>Visible</U> rule if no <U>Alias</U> rule is specified)</FONT></TH>
</TR>
<TR bgcolor="$color2" valign=top>
<TH align=right width=100><INPUT type="checkbox" name="best_all" value="1" onClick="setRule('best_all','best_top')"></TH>
<TH align=left nowrap>Protein is <U>Alias</U> in all Match Groups.</TH>
</TR>
<TR bgcolor="$color2" valign=top>
<TH align=right width=100>or <INPUT type="checkbox" name="best_top" value="1" onClick="setRule('best_top','best_all')"></TH>
<TH align=left nowrap>Protein is <U>Alias</U> only when matched by all peptides in Match Group.<BR>
<FONT style="font-size:11px;font-style:italic;">(Protein will otherwise follow <U>Visible</U> rule if any.)</FONT></TH>
</TR>
<TR bgcolor="$color2" valign=top>
<TH align=right width=100><INPUT type="checkbox" name="old_best" value="1"></TH>
<TH align=left nowrap>Former <U>Alias</U> proteins should remain <U>Visible</U>.
<FONT style="font-size:11px;font-style:italic;">(These proteins will be <U>Hidden</U> if box is not checked)</FONT></TH>
</TR>

<TR bgcolor="$color2">
<TH class="big" align=left colspan=2><BR>&nbsp&nbsp&nbsp&nbsp&nbsp;Rule for other <U>Visible</U> proteins :</TH></TR>
<TR bgcolor="$color2" valign=top><TH align=right width=100><INPUT type="checkbox" name="vis_all" value="1" onClick="setRule('vis_all','vis_top')"></TH>
<TH align=left nowrap>Proteins are <U>Visible</U> in all Match Groups.</TH></TR>
<TR bgcolor="$color2" valign=top><TH align=right width=100>or <INPUT type="checkbox" name="vis_top" value="1" onClick="setRule('vis_top','vis_all')"></TH>
<TH align=left nowrap>Proteins are <U>Visible</U> only when matched by all peptides in Match Group.</TH></TR>

<TR bgcolor="$color2">
<TH class="big" align=left colspan=2><BR>&nbsp&nbsp&nbsp&nbsp&nbsp;Rule for <U>Hidden</U> proteins :
&nbsp<FONT style="font-size:11px;font-style:italic;">(Proteins will not be <U>Hidden</U> where alone in Match Group)</FONT></TH></TR>
<TR bgcolor="$color2" valign=top><TH align=right width=100><INPUT type="checkbox" name="hide_all" value="1" onClick="setRule('hide_all','hide_low')"></TH>
<TH align=left nowrap>Proteins are <U>Hidden</U> everywhere.</TH></TR>
<TR bgcolor="$color2" valign=top><TH align=right width=100>or <INPUT type="checkbox" name="hide_low" value="1" onClick="setRule('hide_low','hide_all')"></TH>
<TH align=left nowrap>Proteins are <U>Hidden</U> only where another
<ACRONYM onmouseover="popup('Another protein with <B>same</B> or <B>more</B> matching peptides exists in Match Group.')" onmouseout="popout()">good candidate<SUP>*</SUP></ACRONYM> can be found.</TH></TR>
</TABLE>
<BR>
<FONT class="title2" color="red">3. Project-wide protein visibility rule may contradict or propagate the changes made.</FONT><BR>
<INPUT type="checkbox" name="overridePR" value="1" checked>&nbsp<FONT class="title3">Do not apply Project-wide rule on these changes.</FONT>

<BR><BR><INPUT type="submit" name="save" value="   Proceed   ">
</DIV>

</FORM>
</CENTER>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup()</SCRIPT>

</BODY>
</HTML>
|;


#################################################
####<<<Processing form: Storing data in DB>>>####
#################################################
sub processForm {

	####<Parameters>####
	my ($analysisID,$group)=split(/:/,param('groupInfo'));
	my $item=param('item');
	my $itemID=param('id_item');
	my ($bestProteinID)=param('check_best');
	my @visProteins=param('check_vis');
	my ($parentItem,$parentID)=split(/:/,lc(param('parentItem')));
	my $bestAction=(param('best_all'))? 2 : (param('best_top'))? 1 : 0;
	my $visAction=(param('vis_all'))? 2 : (param('vis_top'))? 1 : 0;
	my $hideAction=(param('hide_all'))? 2 : (param('hide_low'))? 1 : 0;
	my $oldBest2Vis=param('old_best');
	my $overridePR=(param('overridePR'))? param('overridePR') : 0;
	my (%totMatchGr,%modMatchGr);

	####<Common query>####
	my $sthUp=$dbh->prepare("UPDATE ANALYSIS_PROTEIN SET VISIBILITY=? WHERE ID_ANALYSIS=? AND ID_PROTEIN=?");

	####<Dealing with current Analysis>####
	$totMatchGr{"$analysisID:$group"}=1; # for current
	my %protVisGr;
	foreach my $protID (@visProteins) {
		if ($protID==$bestProteinID) {
			$protVisGr{$protID}=2;
			next;
		}
		$protVisGr{$protID}=1;
	}

	my $sthGr=$dbh->prepare("SELECT ID_PROTEIN,VISIBILITY FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=$analysisID AND MATCH_GROUP=$group");
	$sthGr->execute();
	while (my($protID,$vis)=$sthGr->fetchrow_array) {
		$protVisGr{$protID}=0 unless $protVisGr{$protID};
		next if $vis==$protVisGr{$protID}; # no change in visibility
		$sthUp->execute($protVisGr{$protID},$analysisID,$protID);
		$modMatchGr{"$analysisID:$group"}=1; # group has been updated
	}
	$sthGr->finish;

	unless ($bestAction || $visAction) {
		$sthUp->finish;

		####<Applying Project-wide protein visibility rule?>####
		&promsMod::applyVisibilityRule($dbh,$projectID) unless $overridePR; # updates classifications as well

		$dbh->commit;
		return (scalar keys %totMatchGr, scalar keys %modMatchGr);
	}


	####<Finding other Analyses>####
	my $query;
	if ($parentItem eq 'sample') {
		$query="SELECT ANALYSIS.ID_ANALYSIS,MATCH_GROUP,VISIBILITY,NUM_PEP FROM ANALYSIS_PROTEIN,ANALYSIS WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ID_SAMPLE=$parentID AND ID_PROTEIN=?";
	}
	elsif ($parentItem eq 'experiment') {
		$query="SELECT ANALYSIS.ID_ANALYSIS,MATCH_GROUP,VISIBILITY,NUM_PEP FROM ANALYSIS_PROTEIN,ANALYSIS,SAMPLE WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_EXPERIMENT=$parentID AND ID_PROTEIN=?";
	}
	else { # project
		$query="SELECT ID_ANALYSIS,MATCH_GROUP,VISIBILITY,NUM_PEP FROM ANALYSIS_PROTEIN,PROTEIN WHERE ANALYSIS_PROTEIN.ID_PROTEIN=PROTEIN.ID_PROTEIN AND ID_PROJECT=$parentID AND PROTEIN.ID_PROTEIN=?";
	}
	my $sthAna=$dbh->prepare($query);
	my (%visibility,%groupList);
	foreach my $protID (@visProteins) {
		$sthAna->execute($protID);
		while (my ($anaID,$matchGr,$vis,$numPep)=$sthAna->fetchrow_array) {
			next if $anaID==$analysisID; # skip current Analysis
			$groupList{$anaID}{$matchGr}++; # number of proteins in group
			@{$visibility{$protID}{$anaID}}=($vis,$matchGr,$numPep);
			$totMatchGr{"$anaID:$matchGr"}=1;
		}
	}
	$sthAna->finish;

	####<Finding old Best proteins>####
	my $sthOldBest=$dbh->prepare("SELECT ID_PROTEIN FROM ANALYSIS_PROTEIN WHERE VISIBILITY=2 AND ID_ANALYSIS=? AND MATCH_GROUP=?");
	my $sthMaxPep=$dbh->prepare("SELECT MAX(NUM_PEP) FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=? AND MATCH_GROUP=?");
	my (%oldBest,%maxNumPep);
	foreach my $anaID (keys %groupList) {
		foreach my $matchGr (keys %{$groupList{$anaID}}) {
			$sthOldBest->execute($anaID,$matchGr);
			($oldBest{$anaID}{$matchGr})=$sthOldBest->fetchrow_array;
			$sthMaxPep->execute($anaID,$matchGr);
			($maxNumPep{$anaID}{$matchGr})=$sthMaxPep->fetchrow_array;
		}
	}
	$sthOldBest->finish;
	$sthMaxPep->finish;

	####<Updating new Best protein>####
	if ($bestAction) {
		foreach my $anaID (keys %{$visibility{$bestProteinID}}) {
			my ($vis,$matchGr,$numPep)=@{$visibility{$bestProteinID}{$anaID}};
			if ($oldBest{$anaID}{$matchGr} != $bestProteinID) {
				if ($bestAction==2 || $maxNumPep{$anaID}{$matchGr}==$numPep) { # force best or new best matched by all pep
					my $newVisOld=($oldBest2Vis)? 1 : 0;
					$sthUp->execute($newVisOld,$anaID,$oldBest{$anaID}{$matchGr});
					$sthUp->execute(2,$anaID,$bestProteinID);
					$modMatchGr{"$anaID:$matchGr"}=1;
				}
				elsif ($bestAction==1 && $vis==0) { # set to visible only
					$sthUp->execute(1,$anaID,$bestProteinID);
					$modMatchGr{"$anaID:$matchGr"}=1;
				}
			}
		}
	}

	####<Updating Visible proteins>####
	if ($visAction) {
		foreach my $protID (@visProteins) {
			next if ($bestAction && $protID==$bestProteinID);
			foreach my $anaID (keys %{$visibility{$protID}}) {
				my ($vis,$matchGr,$numPep)=@{$visibility{$protID}{$anaID}};
				next if $vis != 0; # no need to update (could even be 2!)
				if ($visAction==2 || $maxNumPep{$anaID}{$matchGr}==$numPep) {
					$sthUp->execute(1,$anaID,$protID);
					$modMatchGr{$matchGr}=1;
				}
			}
		}
	}

	####<Updating Hidden proteins>####
	if ($hideAction) {
		my $sthCand=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=? AND MATCH_GROUP=? AND NUM_PEP>=?");
		my $sthGr2=$dbh->prepare("SELECT ID_PROTEIN FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=? AND MATCH_GROUP=? ORDER BY VISIBILITY DESC,NUM_PEP DESC,SCORE DESC,ID_PROTEIN ASC");
		foreach my $protID (sort{$protVisGr{$a}<=>$protVisGr{$b}} keys %protVisGr) {
			last if $protVisGr{$protID}>0; # proceed with hidden protein only
			foreach my $anaID (keys %{$visibility{$protID}}) {
				my ($vis,$matchGr,$numPep)=@{$visibility{$protID}{$anaID}};
				next unless $vis; # already hidden
				$sthCand->execute($anaID,$matchGr,$numPep);
				my ($numCand)=$sthCand->fetchrow_array; # number of best candidates (including protID itself)
				if (($hideAction==2 && $groupList{$anaID}{$matchGr}>1) || $numCand>=2) {
					if ($vis==2) { # protID is best of group => replace
						##<Checking if another best is available>## (No change occurs if all proteins in Match Gr must be hidden)
						$sthGr2->execute($anaID,$matchGr);
						while (my($pID)=$sthGr2->fetchrow_array) {
							next if (defined($protVisGr{$pID}) && $protVisGr{$pID}==0); # cannot be another 'to-be_hidden' nor protID
							$sthUp->execute(0,$anaID,$protID);
							$sthUp->execute(2,$anaID,$pID);
							$modMatchGr{$matchGr}=1;
							last;
						}
					}
					else { # vis==1
						$sthUp->execute(0,$anaID,$protID);
						$modMatchGr{$matchGr}=1;
					}
				}
			}
		}
		$sthCand->finish;
		$sthGr2->finish;
	}

	$sthUp->finish;

	####<Applying Project-wide protein visibility rule?>####
	&promsMod::applyVisibilityRule($dbh,$projectID) unless $overridePR; # updates classifications as well

	$dbh->commit;

	return (scalar keys %totMatchGr, scalar keys %modMatchGr);
}

####>Revision history
# 2.1.1 Red font and checking by default for project-wide rule option (FY 03/09/13)