#!/usr/local/bin/perl -w

################################################################################
# listConflicts.cgi    1.0.7                                                   #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Display protein conflicts into a set of analyses                             #
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
use strict;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;

#############################
####>Fetching parameters<####
#############################
if (param('AJAX')) {
	&ajaxGetProteinsMatchGroup;
	exit;
}
my $itemID = param('ID') or die 'Missing ID';
my $itemType = uc(param('TYPE')) or die 'Missing item';
my $action=param('ACT') || 'ambiguity'; # 'switch'; #

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

my $projectID = &promsMod::getProjectID($dbh,$itemID,$itemType);
my @userInfos = &promsMod::getUserInfo($dbh,$ENV{'REMOTE_USER'});
my $projectAccess=${$userInfos[2]}{$projectID};
my $userCanEditMatchGroups = ($projectAccess eq 'guest')?0:1;

my ($itemName) = $dbh->selectrow_array("SELECT NAME FROM $itemType WHERE ID_$itemType=$itemID");
my $displayType=&promsMod::getItemType($itemType);
my %protData;
&fetchProteins(\%protData); # Popuplate %protData

$dbh->disconnect;

#
# HTML #
print header;
warningsToBrowser(1);
my ($light,$dark)=&promsConfig::getRowColors;
my $bgColor = $dark;
print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.popup {background-color:#FFFFFF;border:solid 3px #999999;padding:5px;box-shadow:10px 10px 20px #808080;position:absolute;display:none;} //z-index:999; <- conflits with mouseover popup
</STYLE>
<SCRIPT language="Javascript">
var anaId;
var matchGroupId;
function setRowVisibility(rowId) {
    var img=document.getElementById('img_'+rowId);
    var rowBlock=document.getElementById('row_'+rowId);
    if (rowBlock.style.display=="none") { // block is hidden? => show block
		img.src='$promsPath{images}/minus1.gif';
		rowBlock.style.display="";
    }
    else { //block is shown => hide block
		img.src='$promsPath{images}/plus.gif';
		rowBlock.style.display="none";
    }
}
function sequenceView(id_analysis,id_protein){
    var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+id_analysis+"&id_prot="+id_protein+"&msdata="+top.showMSdata;
    top.openProtWindow(winLocation);
}
function editMatchGroup() {
    if (!anaId \|\| !matchGroupId){
		alert('Select an analysis');
    }
	else {
		window.location="./editMatchGroup.cgi?id_item=$itemID&item=$itemType&id_analysis="+anaId+"&group="+matchGroupId;
    }
}
function enableEditMatchGroup(analysis,matchGroup){
    anaId = analysis;
    matchGroupId = matchGroup;
    document.getElementById('editMatchGroupButton').disabled = false;
}
function disableEditMatchGroup(){
    document.getElementById('editMatchGroupButton').disabled = true;
}
function checkOnLoad(){
    var radios = document.getElementsByName('matchGroup');
    for(var i=0;i<radios.length;i++){
		if(radios[i].checked){
			setRowVisibility(radios[i].value);
			radios[i].click(); // to update global variables
		}
    }
}

// AJAX --->
var selectedImg;
var XHR=null;
function ajaxGetProteinsMatchGroup(anaID,matchGr,protID) {
	if (selectedImg) {selectedImg.src='$promsPath{images}/plus.gif';}
	var grDIV=document.getElementById('groupDIV');
	var linkImg=document.getElementById('img_'+protID+':'+anaID);
	if (grDIV.style.display == 'none' \|\| selectedImg.id != linkImg.id) {
		linkImg.src='$promsPath{images}/minus1.gif';
		grDIV.style.display='block';
	}
	else {
		linkImg.src='$promsPath{images}/plus.gif';
		grDIV.style.display='none';
		return;
	}
	selectedImg=linkImg;
	//wait image
	grDIV.innerHTML='<BR><TABLE><TR><TD>&nbsp;&nbsp;&nbsp;</TD><TD><IMG src="$promsPath{images}/scrollbarGreen.gif"></TD><TD>&nbsp;&nbsp;&nbsp;</TD></TR></TABLE><BR>';
	var divX=getDomOffset(linkImg,'offsetLeft') - 300;
    if (divX < 0) divX=0;
    var divY=getDomOffset(linkImg,'offsetTop') + linkImg.offsetHeight;
	grDIV.style.left = divX + 'px';
	grDIV.style.top = divY + 'px';

	//If XHR object already exists, the request is canceled & the object is deleted
	if(XHR && XHR.readyState != 0){
		XHR.abort();
		delete XHR;
	}
	//Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/listConflicts.cgi?AJAX=MGr&anaID="+anaID+"&matchGr="+matchGr+"&protID="+protID,true);
    XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			grDIV.innerHTML=(XHR.responseText);
		}
	}
	XHR.send(null);
}
function getXMLHTTP(){
	var xhr=null;
	if(window.XMLHttpRequest) {// Firefox & others
		xhr = new XMLHttpRequest();
	}
	else if(window.ActiveXObject){ // Internet Explorer
		try {
		  xhr = new ActiveXObject("Msxml2.XMLHTTP");
		} catch (e) {
			try {
				xhr = new ActiveXObject("Microsoft.XMLHTTP");
			} catch (e1) {
				xhr = null;
			}
		}
	}
	else { // XMLHttpRequest not supported by browser
		alert("Your browser does not support XMLHTTPRequest objects...");
	}
	return xhr;
}
// <--- AJAX
function getDomOffset( Obj, Prop ) {
	var iVal = 0;
	while (Obj && Obj.tagName != 'BODY') {
		eval('iVal += Obj.' + Prop + ';');
		Obj = Obj.offsetParent;
	}
	return iVal;
}
|;
&promsMod::popupInfo();
my $titlePiece =($action eq 'switch')? 'Changes in' : 'Ambiguities for';
print qq
|</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" onload="checkOnLoad();">
<CENTER>
<FONT class="title">List of $titlePiece Visibility Assignment in $displayType <FONT color="red">$itemName</FONT></FONT><BR><BR>
|;

if (scalar keys %protData) {
	print "<FONT class=\"title\">",scalar keys %protData," cases found.</FONT><BR><B>(Proteins identified with only 1 peptide are not considered)</B><BR>\n";
	if ($userCanEditMatchGroups) {
		print qq|<INPUT id="editMatchGroupButton" type="button" value="Edit match group" onclick="editMatchGroup();" disabled>|;
	}
	if ($action eq 'switch') {
		print qq
|<TABLE border=0 cellspacing=0 cellpadding=2 bgcolor=$dark>
<TR><TH class="rbBorder">&nbsp;Protein&nbsp;</TH>
<TH class="rbBorder">&nbsp;<ACRONYM onmouseover="popup('Number of Analyses where protein is <B>visible</B>')" onmouseout="popout()"># Visible (# Peptides)</ACRONYM>&nbsp;</TH>
<TH class="rbBorder">&nbsp;<ACRONYM onmouseover="popup('Number of Analyses where protein is <B>hidden</B>')" onmouseout="popout()"># Hidden (# Peptides)</ACRONYM>&nbsp;</TH><TH class="bBorder" width="600px">&nbsp;Description - Species&nbsp;</TH></TR>
|;
	}
	else { # ambiguity
		print qq
|<TABLE border=0 cellspacing=0 cellpadding=2 bgcolor=$dark>
<TR><TH class="rbBorder">&nbsp;Protein&nbsp;</TH>
<TH class="rbBorder">&nbsp;|;
		if ($itemType eq 'ANALYSIS') {print "# Peptides";}
		else {print qq |<ACRONYM onmouseover="popup('Number of Analyses where protein is <B>visible</B>')" onmouseout="popout()"># Analyses (# Peptides)</ACRONYM>|;}
		print qq
|&nbsp;</TH>
<TH class="bBorder" width="600px">&nbsp;Description - Species&nbsp;</TH></TR>
|;

	}

    #my $rowID = 0;
    foreach my $protID (sort {$protData{$a}{alias} cmp $protData{$b}{alias}} keys %protData){
		$bgColor = ($bgColor eq $dark)? $light: $dark;
		print "<TR bgcolor=$bgColor><TD>";
		if ($itemType eq 'ANALYSIS') {
			my $analysis=$protData{$protID}{visible}[0]; # only 1 analysis
			print "&nbsp;<B>$protData{$protID}{alias}</B>&nbsp;</TD><TD>",&getAnalysisHierarchyString($protID,$analysis);
			print "<IMG border=0 align=top id=\"img_$protID:$analysis->{ID}\" src=\"$promsPath{images}/plus.gif\" onclick=\"ajaxGetProteinsMatchGroup($analysis->{ID},$analysis->{matchGroup},$protID)\" onmouseover=\"popup('Click to display Match group')\" onmouseout=\"popout()\">" if $analysis->{protInMG} > 1;
			print "</TD>";
		}
		else {
			print "&nbsp;<IMG border=0 align=top id=\"img_$protID\" src=\"$promsPath{images}/plus.gif\" onclick=\"javascript:setRowVisibility($protID)\" onmouseover=\"popup('Click to display Analyses')\" onmouseout=\"popout()\">";
			print "<B>$protData{$protID}{alias}</B></TD><TD align=\"center\">", scalar @{$protData{$protID}{visible}};
		}
		print "</TD><TD align=\"center\">",scalar @{$protData{$protID}{invisible}},"</TD>" if $action eq 'switch';
		print "</TD><TD align=\"left\">",$protData{$protID}{description}, ' <FONT class="org">',$protData{$protID}{organism},"</FONT></TD></TR>\n";

		# Hidden row (contains details)
		if ($itemType ne 'ANALYSIS') {
			print "<TR id=\"row_$protID\" bgcolor=$bgColor style=\"display:none\"><TD></TD><TD valign=\"top\" nowrap>";
			foreach my $analysis (sort sortAnalyses @{$protData{$protID}{visible}}){
				print &getAnalysisHierarchyString($protID,$analysis);
				print "<IMG border=0 align=top id=\"img_$protID:$analysis->{ID}\" src=\"$promsPath{images}/plus.gif\" onclick=\"ajaxGetProteinsMatchGroup($analysis->{ID},$analysis->{matchGroup},$protID)\" onmouseover=\"popup('Click to display Match group')\" onmouseout=\"popout()\">" if $analysis->{protInMG} > 1;
				print "&nbsp;<BR>";
			}
			print "</TD>";
			if ($action eq 'switch') {
				print "<TD valign=\"top\" nowrap>";
				foreach my $analysis (sort sortAnalyses @{$protData{$protID}{invisible}}){
					print &getAnalysisHierarchyString($protID,$analysis),"<IMG border=0 align=top id=\"img_$protID:$analysis->{ID}\" src=\"$promsPath{images}/plus.gif\" onclick=\"ajaxGetProteinsMatchGroup($analysis->{ID},$analysis->{matchGroup},$protID)\" onmouseover=\"popup('Click to display Match group')\" onmouseout=\"popout()\">&nbsp;<BR>";
				}
				print "</TD>";
			}
			print "<TD></TD>";
		}
		print "</TR>\n";
    }

    print "</TABLE>\n";

}
else {
    print "<FONT class=\"title2\">No ambiguities detected</FONT>\n";
}

print qq
|</CENTER>
<DIV id="groupDIV" class="popup" style="display:none"></DIV>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">
setPopup();
</SCRIPT>
</BODY>
</HTML>
|;


sub fetchProteins{
	my ($refProtData)=@_;

    my $protDataSQL = 'P.ID_PROTEIN, P.ALIAS, P.PROT_DES, P.ORGANISM, A.ID_ANALYSIS, A.NAME, AP.VISIBILITY, AP.MATCH_GROUP, A.DISPLAY_POS, AP.CONF_LEVEL, AP.NUM_PEP';

    my @sthProt = ($itemType eq 'PROJECT')? ($dbh->prepare("SELECT $protDataSQL, E.NAME, E.DISPLAY_POS, S.NAME, S.DISPLAY_POS
							FROM EXPERIMENT E, SAMPLE S, ANALYSIS A, ANALYSIS_PROTEIN AP, PROTEIN P
							WHERE E.ID_PROJECT=?
							AND E.ID_EXPERIMENT=S.ID_EXPERIMENT
							AND S.ID_SAMPLE=A.ID_SAMPLE
							AND S.ID_SPOT IS NULL
							AND A.ID_ANALYSIS=AP.ID_ANALYSIS
							AND AP.ID_PROTEIN=P.ID_PROTEIN AND AP.NUM_PEP > 1"),
											$dbh->prepare("SELECT $protDataSQL, E.NAME, E.DISPLAY_POS, G.NAME, G.DISPLAY_POS, SP.NAME, SP.ID_SPOT
							FROM EXPERIMENT E, GEL2D G, SPOT SP, SAMPLE S, ANALYSIS A, ANALYSIS_PROTEIN AP, PROTEIN P
							WHERE E.ID_PROJECT=?
							AND G.ID_EXPERIMENT=E.ID_EXPERIMENT
							AND SP.ID_GEL2D=G.ID_GEL2D
							AND SP.ID_SPOT=S.ID_SPOT
							AND S.ID_SAMPLE=A.ID_SAMPLE
							AND A.ID_ANALYSIS=AP.ID_ANALYSIS
							AND AP.ID_PROTEIN=P.ID_PROTEIN AND AP.NUM_PEP > 1")) :
                ($itemType eq 'EXPERIMENT')? ($dbh->prepare("SELECT $protDataSQL, S.NAME, S.DISPLAY_POS
							FROM SAMPLE S, ANALYSIS A, ANALYSIS_PROTEIN AP, PROTEIN P
							WHERE S.ID_EXPERIMENT=?
							AND S.ID_SAMPLE=A.ID_SAMPLE
							AND A.ID_ANALYSIS=AP.ID_ANALYSIS
							AND AP.ID_PROTEIN=P.ID_PROTEIN
							AND S.ID_SPOT IS NULL AND AP.NUM_PEP > 1"),
											$dbh->prepare("SELECT $protDataSQL, G.NAME, G.DISPLAY_POS, SP.NAME, SP.ID_SPOT
							FROM GEL2D G, SPOT SP, SAMPLE S, ANALYSIS A, ANALYSIS_PROTEIN AP, PROTEIN P
							WHERE G.ID_EXPERIMENT=?
							AND SP.ID_GEL2D=G.ID_GEL2D
							AND SP.ID_SPOT=S.ID_SPOT
							AND S.ID_SAMPLE=A.ID_SAMPLE
							AND A.ID_ANALYSIS=AP.ID_ANALYSIS
							AND AP.ID_PROTEIN=P.ID_PROTEIN AND AP.NUM_PEP > 1")) :
                ($itemType eq 'SAMPLE')? ($dbh->prepare("SELECT $protDataSQL
							FROM ANALYSIS A, ANALYSIS_PROTEIN AP, PROTEIN P
							WHERE A.ID_SAMPLE=?
							AND A.ID_ANALYSIS=AP.ID_ANALYSIS
							AND AP.ID_PROTEIN=P.ID_PROTEIN AND AP.NUM_PEP > 1")) :
                ($itemType eq 'GEL2D')? ($dbh->prepare("SELECT $protDataSQL, SP.NAME, SP.ID_SPOT
						    FROM SPOT SP, SAMPLE S, ANALYSIS A, ANALYSIS_PROTEIN AP, PROTEIN P
							WHERE SP.ID_GEL2D=?
							AND SP.ID_SPOT=S.ID_SPOT
							AND S.ID_SAMPLE=A.ID_SAMPLE
							AND A.ID_ANALYSIS=AP.ID_ANALYSIS
							AND AP.ID_PROTEIN=P.ID_PROTEIN AND AP.NUM_PEP > 1")) :
                ($itemType eq 'SPOT')? ($dbh->prepare("SELECT $protDataSQL, S.NAME, S.DISPLAY_POS
							FROM SAMPLE S, ANALYSIS A, ANALYSIS_PROTEIN AP, PROTEIN P
							WHERE S.ID_SPOT=?
							AND S.ID_SAMPLE=A.ID_SAMPLE
							AND A.ID_ANALYSIS=AP.ID_ANALYSIS
							AND AP.ID_PROTEIN=P.ID_PROTEIN AND AP.NUM_PEP > 1")) :
				($itemType eq 'ANALYSIS')? ($dbh->prepare("SELECT $protDataSQL
							FROM ANALYSIS A, ANALYSIS_PROTEIN AP, PROTEIN P
							WHERE A.ID_ANALYSIS=?
							AND A.ID_ANALYSIS=AP.ID_ANALYSIS
							AND AP.ID_PROTEIN=P.ID_PROTEIN AND AP.NUM_PEP > 1")) :
                ();
	#my $sthMGA=$dbh->prepare("SELECT MATCH_GROUP,COUNT(*),MAX(NUM_PEP),ID_PROTEIN FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=? GROUP BY MATCH_GROUP");
	#my $sthMGA=$dbh->prepare("SELECT MATCH_GROUP,COUNT(*),MAX(NUM_PEP),ID_PROTEIN FROM
	#							(SELECT * FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=511 ORDER BY MATCH_GROUP, VISIBILITY DESC) AS O_AP
	#							GROUP BY MATCH_GROUP"); # to make sure ID_PROTEIN has best visibility (=2)
	my $sthMGA=($action eq 'switch')?
				$dbh->prepare("SELECT MATCH_GROUP,COUNT(*),MAX(NUM_PEP),ID_PROTEIN FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=? GROUP BY MATCH_GROUP") :
				$dbh->prepare("SELECT MATCH_GROUP,NUM_PEP,COUNT(*),MAX(VISIBILITY) FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=? GROUP BY MATCH_GROUP,NUM_PEP");

    my (%matchGrAna,%matchGrPepAna);
    foreach my $sth (@sthProt){
		$sth->execute($itemID);
		while (my ($protID,$protAlias,$protDes,$protOrganism,$anaID,$anaName,$visibility,$matchGroup,$anaDisplayPos,$confLevel,$numPep,@parentInfos) = $sth->fetchrow_array) {
			unless ($matchGrAna{$anaID}) {
				$sthMGA->execute($anaID);
				if ($action eq 'switch') {
					while (my ($mGr,$numProt,$maxPep)=$sthMGA->fetchrow_array) {
						@{$matchGrAna{$anaID}{$mGr}}=($numProt,$maxPep);
					}
				}
				else { # ambiguity
					my (%numProtMG,%maxPepMG);
					while (my ($mGr,$numPep,$numProt,$maxVis)=$sthMGA->fetchrow_array) {
						$numProtMG{$mGr}+=$numProt;
						if (!$maxPepMG{$mGr} || $maxPepMG{$mGr} < $numPep) {$maxPepMG{$mGr}=$numPep;}
						@{$matchGrPepAna{$anaID}{$mGr}{$numPep}}=($numProt,$maxVis);
					}
					foreach my $mGr (keys %numProtMG) {@{$matchGrAna{$anaID}{$mGr}}=($numProtMG{$mGr},$maxPepMG{$mGr});}
				}
			}
			next if ($visibility==0 && $numPep < $matchGrAna{$anaID}{$matchGroup}[1]); # OK for prot to be hidden
			my $analysis = { ID => $anaID, name => $anaName, visibility => $visibility, matchGroup => $matchGroup, displayPos => $anaDisplayPos, confLevel => $confLevel, protInMG => $matchGrAna{$anaID}{$matchGroup}[0], numPep => $numPep};
			while (my ($parentName, $parentDisplayPos) = splice @parentInfos,0,2) {
				push @{$analysis->{parents}}, { name => $parentName , displayPos => $parentDisplayPos };
			}
			push @{$refProtData->{$protID}{analyses}}, $analysis;
			$refProtData->{$protID}{alias} = $protAlias;
			$refProtData->{$protID}{description} = $protDes;
			$refProtData->{$protID}{organism} = $protOrganism;
		}
		$sth->finish;
    }
	$sthMGA->finish;


	foreach my $protID (keys %{$refProtData}) {

		my @visibleAna;
		my @invisibleAna;

		if ($action eq 'switch') {
			foreach my $analysis (@{$refProtData->{$protID}{analyses}}){
				if ($analysis->{visibility} > 0) {push @visibleAna, $analysis;}
				else {push @invisibleAna, $analysis;}
			}
			if (scalar @visibleAna && scalar @invisibleAna) { # there is switch
				$refProtData->{$protID}{visible} = \@visibleAna;
				$refProtData->{$protID}{invisible} = \@invisibleAna;
			}
			else {delete $refProtData->{$protID};}
		}

		elsif ($action eq 'ambiguity') {
			foreach my $analysis (@{$refProtData->{$protID}{analyses}}) {
				my $refMGrPepAnaData=$matchGrPepAna{$analysis->{ID}}{$analysis->{matchGroup}}{$analysis->{numPep}};
				if ($analysis->{visibility} > 0 && $refMGrPepAnaData->[0] > 1) { # visible && exist other prot with same num pep
					push @visibleAna,$analysis;
				}
				###elsif ($analysis->{visibility}==0 && $refMGrPepAnaData->[1] >= 1) { # hidden && exist visible prot with same num pep
				###	push @invisibleAna,$analysis;
				###} <-- COMMENTED!!! Ambiguities are ignored for hidden to prevent too many results
			}
			if (scalar @visibleAna || scalar @invisibleAna) { # there is ambuiguity
				$refProtData->{$protID}{visible} = \@visibleAna; # if scalar @visibleAna;
				$refProtData->{$protID}{invisible} = \@invisibleAna; # if scalar @invisibleAna;
			}
			else {delete $refProtData->{$protID};}
		}
	}
}

sub getAnalysisHierarchyString{
    my ($protID,$analysis) = @_;

    my $string='';
    if ($userCanEditMatchGroups && $analysis->{protInMG} > 1) {
		$string = "<INPUT type=\"radio\" name=\"matchGroup\" value=\"$protID\" onclick=\"enableEditMatchGroup($analysis->{ID},$analysis->{matchGroup})\">";
    }
    elsif ($analysis->{protInMG} == 1) {
		$string = '<INPUT type="radio" name="fakeMG" value="0" style="visibility:hidden">'; #"<IMG border=0 align=top src=\"$promsPath{images}/good.gif\">";
    }
    my ($protClass,$popupString) = promsMod::getProteinClass($analysis->{confLevel}, $analysis->{visibility});
	if ($itemType eq 'ANALYSIS') {
		$string .= "<A class=\"$protClass\" href=\"javascript:sequenceView($analysis->{ID},$protID)\" onmouseover=\"popup('$popupString');\" onmouseout=\"popout();\"><FONT class=\"$protClass\">$analysis->{numPep}</FONT></A>&nbsp";
	}
	else {
		if ($analysis->{parents}) {
			$string .= join (' > ', map {$_->{name}} @{$analysis->{parents}}) . ' > ';
		}
		$string .= "<A class=\"$protClass\" href=\"javascript:sequenceView($analysis->{ID},$protID)\" onmouseover=\"popup('$popupString');\" onmouseout=\"popout();\">$analysis->{name}</A>&nbsp;<FONT class=\"$protClass\">($analysis->{numPep})</FONT>&nbsp;";
	}

    return $string;
}

sub sortAnalyses{
    my $sort = 0;

    # Firstly sorting by parents if they exist
    if ($a->{parents} && $b->{parents}) {
		my $i = 0;
		while ($sort == 0) {
			if ($a->{parents}[$i] && $b->{parents}[$i]) {
			$sort = $a->{parents}[$i]{displayPos} <=> $b->{parents}[$i]{displayPos};
			}
			else { # 1 or both have no more parents. If only 1 have no more, this 1 is considered as first
				if ($a->{parents}[$i]) {return 1;}
				elsif ($b->{parents}[$i]){return -1;}
				else {last;} # go to sorting by names
			}
			$i++;
		}
    }

    # Else simply sorting by names
    if($sort == 0){
	$sort =  $a->{displayPos} <=> $b->{displayPos} ;
    }

    return $sort;
}

sub ajaxGetProteinsMatchGroup {
	print header(-charset=>'UTF-8');
	warningsToBrowser(1);

	my $analysisID=param('anaID');
	my $matchGroup=param('matchGr'); # can be undef
	my $refProtID=param('protID');

	my $dbh=&promsConfig::dbConnect;

	my $rootQuery="SELECT P.ID_PROTEIN,ALIAS,PROT_DES,ORGANISM,NUM_PEP,NUM_MATCH,SCORE,VISIBILITY,CONF_LEVEL,PROT_LENGTH FROM ANALYSIS_PROTEIN A,PROTEIN P WHERE A.ID_PROTEIN=P.ID_PROTEIN AND ID_ANALYSIS=$analysisID AND MATCH_GROUP";
	my $sthMG=($matchGroup)? $dbh->prepare("$rootQuery=$matchGroup") : $dbh->prepare("$rootQuery=(SELECT MATCH_GROUP FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=$analysisID AND ID_PROTEIN=$refProtID)");

	my %listProteins;
	$sthMG->execute;
	while (my ($protID,@protData)=$sthMG->fetchrow_array) {
		@{$listProteins{$protID}}=@protData;
	}
	$sthMG->finish;

	my ($anaName)=$dbh->selectrow_array("SELECT NAME FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");

	$dbh->disconnect;

	my ($lightColor,$darkColor)=&promsConfig::getRowColors;
	print qq
|<FONT class="title3">Match group of $listProteins{$refProtID}[0] in Analysis '$anaName'</FONT>&nbsp;&nbsp;<INPUT type="button" class="font11" value="Close" onclick="ajaxGetProteinsMatchGroup($analysisID,null,$refProtID)"/>
<TABLE border=0 cellspacing=0>
<TR bgcolor=$darkColor valign=bottom>
<TH class="rbBorder">&nbsp;Protein&nbsp;</TH>
<TH class="rbBorder">&nbsp;Visibility&nbsp;</TH>
<TH class="rbBorder">&nbsp;Peptides&nbsp;</TH>
<TH class="rbBorder">&nbsp;Score&nbsp;</TH>
<TH class="rbBorder">&nbsp;Length&nbsp;</TH>
<TH class="bBorder" align=left width=500>&nbsp;Description - Species&nbsp;</TH>
</TR>
|;
	my $bgColor=$lightColor;
	foreach my $protID (sort {$listProteins{$b}[6]<=>$listProteins{$a}[6] || $listProteins{$b}[3]<=>$listProteins{$a}[3] || $listProteins{$b}[5]<=>$listProteins{$a}[5] || $listProteins{$a}[8]<=>$listProteins{$b}[8] || $a<=>$b} keys %listProteins){ # numPep>score>length>id
		my ($protClass,$protPopup)=&promsMod::getProteinClass($listProteins{$protID}[7],$listProteins{$protID}[6]); # conf,vis
		my ($colorStrg1,$colorStrg2)=($protID==$refProtID)? ('<FONT color="#DD000">','</FONT>') : ('','');
		my $visibilityCell=($listProteins{$protID}[6]==2)? '<TH>alias</TH>' : ($listProteins{$protID}[6]==1)? '<TH>visible</TH>' : '<TD align="center">hidden</TD>';
		print qq
|<TR class="list" bgcolor="$bgColor" valign=top>
<TD nowrap>&nbsp;<A class="$protClass" href="javascript:sequenceView($analysisID,$protID)">$colorStrg1$listProteins{$protID}[0]$colorStrg2</A>&nbsp;</TD>
$visibilityCell
<TH>$listProteins{$protID}[3]|;
		print "&nbsp;<ACRONYM onmouseover=\"popup('$listProteins{$protID}[4] matches due to sequence repeats')\" onmouseout=\"popout()\">($listProteins{$protID}[4])</ACRONYM>" if $listProteins{$protID}[4]>$listProteins{$protID}[3];
		printf "</TH>\n<TH>&nbsp;%.1f&nbsp;</TH>",$listProteins{$protID}[5];
		print "<TH>$listProteins{$protID}[8]</FONT></TD>\n";
		print "<TD>$listProteins{$protID}[1] <FONT class=\"org\">$listProteins{$protID}[2]</FONT></TD>\n</TR>\n";
		$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
	}
	print "</TABLE>\n";
	exit;
}


####>Revision history
# 1.0.7 Minor modification due to call of listConflicts.cgy by ajaxGetProteinsMatchGroup in sequenceView.cgi (GA 04/10/17)
# 1.0.6 Minor change in score display in &ajaxGetProteinsMatchGroup (PP 17/06/15)
# 1.0.5 Scan for ambiguities or vis/hidden transitions (PP 04/12/14)
# 1.0.4 Changed sort and added Length to AJAX display (PP 22/09/14)
# 1.0.3 Stricter detection of ambiguities & AJAX call to display selected match group in popup div (PP 22/04/14)
# 1.0.2 Guests can not edit match groups anymore (FY 08/10/13)
# 1.0.1 Applying protein font class to analyses names (FY 09/09/13)
# 1.0.0 New script to display protein conflicts in a specified item (FY 02/09/13)
