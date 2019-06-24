#!/usr/local/bin/perl -w

################################################################################
# editClassification.cgi    2.3.1                                              #
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
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use strict;

#print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1); #DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my $date=strftime("%Y-%m-%d %H:%M:%S", localtime);

####################
#    Parameters    #
####################
my $themeID=param('id_theme');
my $projectID=param('id_project');
my $selectedList=0;

###############
#     Main    #
###############
## Connect to the database
my $dbh=&promsConfig::dbConnect;

##########################
####>Updating Changes<####
##########################
if (param('saveTheme')) { # processing theme form
	my $name=param('name');
	$name=~s/['"]//g; # not allowed (potential interference with Javascript)
	my $description=param('des');
	$name=$dbh->quote($name);
	$description=$dbh->quote($description);
	$dbh->do("UPDATE CLASSIFICATION SET NAME=$name,DES=$description,UPDATE_DATE=NOW(),UPDATE_USER='$userID' WHERE ID_CLASSIFICATION=$themeID");
	$dbh->commit;
}
elsif (param('saveList')) { # processing list form
	my $listID=param('id_list');
	my ($name,$description)=(param("name_$listID"),param("des_$listID"));
	$name=~s/['"]//g; # not allowed (potential interference with Javascript)
	$name=$dbh->quote($name);
	$description=$dbh->quote($description);
	if ($listID) { # edit old list
		$dbh->do("UPDATE CATEGORY SET NAME=$name,DES=$description,UPDATE_DATE='$date',UPDATE_USER='$userID' WHERE ID_CATEGORY=$listID");
		$selectedList=$listID;
	}
	else { # add new list
		my ($newListID)=$dbh->selectrow_array("SELECT MAX(ID_CATEGORY) FROM CATEGORY");
		$newListID++;
		my ($displayPos)=$dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM CATEGORY WHERE ID_CLASSIFICATION=$themeID");
		$displayPos++;
		$dbh->do("INSERT INTO CATEGORY (ID_CATEGORY,ID_CLASSIFICATION,NAME,DES,DISPLAY_POS,UPDATE_DATE,UPDATE_USER) VALUES ($newListID,$themeID,$name,$description,$displayPos,'$date','$userID')");
		$selectedList=$newListID;
	}
	$dbh->commit;

}
elsif (param('moveList')) {
	my $listID=abs(param('moveList'));
	my ($sort)=(param('moveList')>0)? ('DESC'): ('ASC');
	my $sthSel=$dbh->prepare("SELECT ID_CATEGORY,DISPLAY_POS FROM CATEGORY WHERE ID_CLASSIFICATION=$themeID ORDER BY DISPLAY_POS $sort");
	$sthSel->execute();
	my ($switchedListID,$newDisplayPos,$oldDisplayPos);
	while ((my $lID,$oldDisplayPos)=$sthSel->fetchrow_array) {
		last if $lID==$listID;
		$switchedListID=$lID;
		$newDisplayPos=$oldDisplayPos;
	}
	$sthSel->finish;
	if ($switchedListID) { # min/max displayPos has not been reached
		$dbh->do("UPDATE CATEGORY SET DISPLAY_POS=$newDisplayPos,UPDATE_DATE='$date',UPDATE_USER='$userID' WHERE ID_CATEGORY=$listID");
		$dbh->do("UPDATE CATEGORY SET DISPLAY_POS=$oldDisplayPos,UPDATE_DATE='$date',UPDATE_USER='$userID' WHERE ID_CATEGORY=$switchedListID");
		$dbh->commit;
	}
	$selectedList=$listID;
}
elsif (param('deleteList')) { # deleting a list
	my $strgParam=param('deleteList');
	my $listID;
	my $themeDel=0;
	if ($strgParam=~/:/) {
		foreach my $strgID (split(/,/, $strgParam)) {
			my ($type, $id)=split(/:/, $strgID);
			next if ($type ne "LI" );
			$listID=$id;
		}
		foreach my $strgID (split(/,/, $strgParam)) {
			my ($type, $id)=split(/:/, $strgID);
			if ($type eq "GO") {
				my ($goID,$item)=split(/-/,$id);
				my ($goParamStrg)=$dbh->selectrow_array("SELECT PARAM_STRG FROM GO_ANALYSIS where ID_GOANALYSIS=$goID");
				my ($strgUpdate,$strgDisp)=($item eq "prot")? ("protSetCategory=$listID", "protSetCategory=-1") : ("popCategory=$listID", "popCategory=-1");
				$goParamStrg=~s/$strgUpdate/$strgDisp/;
				$dbh->do("UPDATE GO_ANALYSIS SET PARAM_STRG='$goParamStrg' where ID_GOANALYSIS=$goID");
			}
			elsif ($type eq "AN") {
				my ($type,$item,$annotID,$rank)=split(/-/,$id);
				my $sthSelAnnot=$dbh->prepare("SELECT ID_ANNOTATIONSET,RANK FROM ANNOTATIONSET WHERE ID_EXPLORANALYSIS=? and ANNOT_TYPE like 'PROT%' ORDER BY RANK");
				my ($explorID)=$dbh->selectrow_array("SELECT ID_EXPLORANALYSIS FROM ANNOTATIONSET where ID_ANNOTATIONSET=$annotID");
				$sthSelAnnot->execute($explorID);
				if ($item eq 'LIST') {
					$dbh->do("DELETE FROM ANNOTATIONSET WHERE ID_ANNOTATIONSET=$annotID");
					my $i=0;
					while (my($annotSetID,$rank)=$sthSelAnnot->fetchrow_array) {
						next if ($annotSetID == $annotID);
						$i++;
						$dbh->do("UPDATE ANNOTATIONSET SET RANK=$i WHERE ID_ANNOTATIONSET=$annotSetID");
					}
				}
				else {
					next if $themeDel;
					$dbh->do("DELETE FROM MODIFICATION_SITE WHERE ID_CATEGORY=$listID");
					$dbh->do("DELETE FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$listID");
					$dbh->do("DELETE FROM CATEGORY WHERE ID_CATEGORY=$listID");
					$themeDel++;
				}
			}
		}
	}
	else {
		$listID=$strgParam;
		$dbh->do("DELETE FROM MODIFICATION_SITE WHERE ID_CATEGORY=$listID");
		$dbh->do("DELETE FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$listID");
		$dbh->do("DELETE FROM CATEGORY WHERE ID_CATEGORY=$listID");
	}
	$dbh->commit;
}
elsif (param('usageList')) { # display list involved in analysis
	&ajaxUsageList();
	exit;
}
elsif (param('import')) {
	&importList(param('import'));
	exit;
}

###############################
####>Fetching data from DB<####
###############################
my ($nameTheme,$desTheme,$updateDate,$updateUser)=$dbh->selectrow_array("SELECT NAME,DES,UPDATE_DATE,UPDATE_USER FROM CLASSIFICATION WHERE ID_CLASSIFICATION=$themeID");
($nameTheme,$desTheme)=&promsMod::chkDef($nameTheme,$desTheme);
my @listList=&promsMod::getChildrenList($dbh,$themeID,'CLASSIFICATION');
my (%listSize,%listType);
my $sthLS=$dbh->prepare("SELECT C.ID_CATEGORY,LIST_TYPE,COUNT(DISTINCT ID_PROTEIN),COUNT(ID_CATEGORY_PROTEIN) FROM CATEGORY C LEFT JOIN CATEGORY_PROTEIN CP ON CP.ID_CATEGORY=C.ID_CATEGORY WHERE C.ID_CLASSIFICATION=$themeID GROUP BY C.ID_CATEGORY ORDER BY DISPLAY_POS;");
$sthLS->execute;
while (my ($listID,$type,$numProt,$numProteoform)=$sthLS->fetchrow_array) {
	@{$listSize{$listID}}=($numProt,$numProteoform);
	$listType{$listID}=$type;
	$selectedList=$listID unless $selectedList;
}
$sthLS->finish;
#unless ($selectedList) {
#	foreach my $listID (sort{$listList[0]{'LIST'}{$a}[2] <=> $listList[0]{'LIST'}{$b}[2]} keys %{$listList[0]{'LIST'}}){
#		$selectedList=$listID;
#		last;
#	}
#}
my %notDelLists;
#my $sthDC=$dbh->prepare("SELECT DISTINCT(ID_CATEGORY) FROM COMPARISON WHERE ID_PROJECT=$projectID AND ID_CATEGORY IS NOT NULL");
#$sthDC->execute;
#while (my ($listID)=$sthDC->fetchrow_array) {
#	$notDelLists{$listID}=1; # all not deletable lists in project!
#}
#$sthDC->finish;
#my $sthCC=$dbh->prepare("SELECT DISTINCT(CC.ID_CATEGORY) FROM CAT_COMPARISON CC,COMPARISON CO WHERE CC.ID_COMPARISON=CO.ID_COMPARISON AND ID_PROJECT=$projectID");
#$sthCC->execute;
#while (my ($listID)=$sthCC->fetchrow_array) {
#	$notDelLists{$listID}=1; # all not deletable lists in project!
#}
#$sthCC->finish;

##>Deletability<##
my @sthNoDel=(
	$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,COMPARISON CO WHERE CO.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID"), # COMPARISON
	$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,CAT_COMPARISON CC WHERE CC.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID"), # CAT_COMPARISON
	#$dbh->prepare("SELECT CA.PARAM_STRG FROM CATEGORY CA,GO_ANALYSIS GA WHERE GA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID"), # GO_ANALYSIS
	$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,PATHWAY_ANALYSIS PA WHERE PA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID"), # PATHWAY_ANALYSIS
	$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,EXPLORANALYSIS EA WHERE EA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID") # EXPLORANALYSIS
	#$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,GO_ANALYSIS GA WHERE GA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID") # GO_ANALYSIS
);

foreach my $sth (@sthNoDel) {
	$sth->execute;
	while (my ($listID)=$sth->fetchrow_array) {
		$notDelLists{$listID}=1; # all not deletable lists in Theme
	}
}
foreach my $sth (@sthNoDel) {$sth->finish;}

##>GO ANALYSIS/ANNOTATIONSET Deletability (ID_CATEGORY in strg PARAM_STRG/ANNOT_LIST);
my $sthGOParam=$dbh->prepare("SELECT GA.PARAM_STRG FROM EXPERIMENT E
			      INNER JOIN GO_ANALYSIS GA on E.ID_EXPERIMENT=GA.ID_EXPERIMENT
			      where ID_PROJECT=$projectID and GOA_TYPE='graph'");
my $sthExploParam=$dbh->prepare("SELECT EA.ANA_TYPE,A.NAME, A.ANNOT_TYPE, A.ANNOT_LIST from EXPERIMENT E
				 INNER JOIN EXPLORANALYSIS EA on E.ID_EXPERIMENT=EA.ID_EXPERIMENT
				 INNER JOIN ANNOTATIONSET A on EA.ID_EXPLORANALYSIS=A.ID_EXPLORANALYSIS
				 where E.ID_PROJECT=$projectID AND ANNOT_TYPE LIKE 'prot:%'");
my $sthSelCat=$dbh->prepare("SELECT ID_CATEGORY from CATEGORY where ID_CLASSIFICATION=?");
$sthGOParam->execute;
my %usedListInGO;
while (my($paramStrg)=$sthGOParam->fetchrow_array) {
	next if ($paramStrg !~ /protSetCategory/);
	push my @tmpList, split(/;/,$paramStrg);
	foreach my $param(@tmpList) {
		if ($param=~/protSetCategory/) {
			my($name,$value)=split(/=/,$param);
			$usedListInGO{$value}=1;
		}
		if ($param=~/popCategory/) {
			my($name,$value)=split(/=/,$param);
			$usedListInGO{$value}=1 if ($value != "0");
		}
	}
}

$sthExploParam->execute;
while (my ($anaType,$name,$annotType,$annotList)=$sthExploParam->fetchrow_array) {
	if ($anaType eq "PCA") {#PCA
		if ($annotType eq "prot:LIST") {
			$name=~s/#//g;
			$usedListInGO{$name}=1;
		}
	}
	else {#cluster
		if ($annotType eq "prot:LIST") {
			foreach my $listStrg (split(':',$annotType)) {
				$listStrg=~s/#//g;
				$usedListInGO{$listStrg}=1;
			}
		}
		elsif ($annotType eq "prot:THEME") {
			$name=~s/#//g;
			if ($name == $themeID) {
				$sthSelCat->execute($name);
				while (my ($listID)=$sthSelCat->fetchrow_array) {
					$usedListInGO{$listID}=1;
				}
			}
		}
	}
}

foreach my $listID (sort{$listList[0]{'LIST'}{$a}[2] <=> $listList[0]{'LIST'}{$b}[2]} keys %{$listList[0]{'LIST'}}){
	if ($usedListInGO{$listID}) {
		$notDelLists{$listID}=1;
	}
	@{$listSize{$listID}}=(0,0) unless $listSize{$listID}; # avoid to launch warnings for empty lists
}

my $selListIsDeletable=($notDelLists{$selectedList})? 0 : 1;

$dbh->disconnect;
my $updateString=($updateDate)? "Last modified : $updateDate by $updateUser" : '';

#######################
####>Starting HTML<####
#######################
my ($color1,$color2)=&promsConfig::getRowColors;

print header(-'content-encoding'=>'no',-charset=>'utf-8');
warningsToBrowser(1); # writes PERL warnings as HTML comments <!-- --> must be activated after header
print qq
|<HEAD>
<TITLE> Theme : $nameTheme</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
function editTheme(action) {
	var vis1,vis2;
	if (action==1) {vis1='none';vis2='';}
	else {vis1='';vis2='none';}
	document.getElementById('displayTheme').style.display=vis1;
	document.getElementById('editTheme').style.display=vis2;
	document.getElementById('editThemeButton').style.display=vis1;
	document.getElementById('saveThemeButton').style.display=vis2;
}
function checkThemeForm(){
	if (document.editThemeForm.name.value){
		return true;
	}
	else{
		alert ("Type a name for theme.");
		return false;
	}
}
function selectList(id_list,isDeletable) {
	document.getElementById('displayList').style.display='none';
	if (editListMode \|\| addListMode \|\| !id_list) {return;}
	
	if(selList) {
		var oldObject=document.getElementById(selList);
		oldObject.style.color='#000000';
	}
	
	var newObject=document.getElementById(id_list);
	newObject.style.color='#DD0000';
	selList=id_list;
	document.getElementById('delButton').disabled=(isDeletable==1)? false : true;
	document.getElementById('useButton').disabled=(isDeletable==1)? true : false;
}
function editList(action) {
	var vis1,vis2;
	if (action>=0) { // Edit selected List
		vis1='none';vis2='';
		editListMode=1;
		var rowID1='list1_'+selList;
		var rowID2='list2_'+selList;
		document.getElementById(rowID1).style.display='none';
		document.getElementById(rowID2).style.display='';
	}
	else { // cancel
		vis1='';vis2='none';
		if (document.getElementById(selList)){
			document.getElementById(selList).style.color='#DD0000';
		}
		if (addListMode) { // new List
			addListMode=0;
			document.getElementById('list_0').style.display='none';
		}
		else { // editListMode
			editListMode=0;
			var rowID1='list1_'+selList;
			var rowID2='list2_'+selList;
			document.getElementById(rowID1).style.display='';
			document.getElementById(rowID2).style.display='none';
		}
	}
	if(document.getElementById('editListButton')) {
		document.getElementById('editListButton').style.display=vis1;
	}
	document.getElementById('saveListButton').style.display=vis2;
}
function deleteList(item){
	if (item == 'force') {
		var idToDelete=new Array();
		var selGO=document.getElementsByName('GO');
		var selAN=document.getElementsByName('AN');
		idToDelete.push('LI:'+selList);

		if (selGO){
			for (var i=0; i<selGO.length; i++) {
				idToDelete.push('GO:'+selGO[i].value);
			}
		}

		if (selAN) {
			for (var i=0; i<selAN.length; i++) {
				idToDelete.push('AN:'+selAN[i].value);
			}
		}

		var strgParam=idToDelete.join(',');
		window.location="./editClassification.cgi?id_project=$projectID&id_theme=$themeID&deleteList="+strgParam;

	}
	else if (!item) {
		if (confirm ("Do you really want to delete this List?")) {
			window.location="./editClassification.cgi?id_project=$projectID&id_theme=$themeID&deleteList="+selList;
		}
	}
}
function moveList(move){
	if (editListMode \|\| addListMode) {
		alert('Exit Edit Mode first.');
		return;
	}
	var moveList=move*selList;
	window.location="./editClassification.cgi?id_project=$projectID&id_theme=$themeID&moveList="+moveList;
}
function importList() {
	window.location="./editClassification.cgi?id_project=$projectID&id_theme=$themeID&import=1";
}

var XHR = null;
function getXMLHTTP() {
    var xhr = null;
    if (window.XMLHttpRequest) {// Firefox & others
        xhr = new XMLHttpRequest();
    }
    else if (window.ActiveXObject) { // Internet Explorer
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
function usageList() {
	var usageDiv=document.getElementById('displayList');
	usageDiv.style.display='none';
	usageDiv.innerHTML='';
	document.editListForm.id_list.value=selList;
	var paramStrg='id_project=$projectID&id_theme=$themeID&usageList='+selList;
	//If XHR object already exists, the request is canceled & the object is deleted
	if (XHR && XHR.readyState != 0) {
	    XHR.abort();
	    delete XHR;
	}
	    //Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if (!XHR) {
	    return false;
	}

	XHR.open("POST","$promsPath{cgi}/editClassification.cgi",true);
	//Send the proper header information along with the request
	XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
	XHR.setRequestHeader("Content-length", paramStrg.length);
	XHR.setRequestHeader("Connection", "close");
	XHR.onreadystatechange=function() {
	    if (XHR.readyState==4 && XHR.responseText) {
			usageDiv.style.display='';
			usageDiv.innerHTML=XHR.responseText;
			usageDiv.scrollIntoView(true);
	    }
	}
	XHR.send(paramStrg);
}
function checkListForm(){
	var modList;
	if (addListMode) {modList=0;}
	else {modList=selList;}
	var listName='name_'+modList;
	if (!document.getElementById(listName).value) {
		alert("Provide a name for List.");
		return false;
	}
	document.editListForm.id_list.value=modList;
}
var selList=$selectedList;
var editListMode=0;
var addListMode=0;
</SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif' onload="selectList($selectedList,$selListIsDeletable);">
<CENTER>
<FONT class="title">Editing Theme <FONT color="#DD0000">$nameTheme</FONT></FONT>
</CENTER>
<BR>
<FORM method="post" name="editThemeForm" onSubmit="return(checkThemeForm());">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="id_theme" value="$themeID">

<TABLE border=0 width=900>
<TR id="displayTheme"><TD bgcolor=$color2>
<TABLE width=896 cellpadding=2>
<TR height=35><TD width=150 class="title2" align=right valign=top bgcolor=$color2>Name :</TD>
<TD class="title2" bgcolor=$color1>$nameTheme</TD>
</TR>
<TR height=35><TD class="title2" align=right valign=top bgcolor=$color2>Description :</TD>
<TD class="title3" bgcolor=$color1>$desTheme</TD>
</TABLE>
</TD></TR>
<TR id="editTheme" style="display:none;"><TD bgcolor=$color2>
<TABLE width=896 cellpadding=2>
<TR height=35><TD width=150 class="title2" align=right valign=top bgcolor=$color2>Name :</TD>
<TD bgcolor=$color1><INPUT type="text" class="title2" name="name" size="40" value="$nameTheme" maxlength=50></TD>
</TR>
<TR height=35><TD class="title2" align=right valign=top bgcolor=$color2>Description :</TD>
<TD bgcolor=$color1><INPUT class="title3" type="text" name="des" size="72" value="$desTheme" maxlength=100></TD>
</TABLE>
</TD></TR>
<TR id="editThemeButton"><TD valign=top>&nbsp&nbsp&nbsp&nbsp
<INPUT type="button" value="Edit" class="title3" style="width:60px" onclick="editTheme(1)">
</TD></TR>
<TR id="saveThemeButton" style="display:none;"><TD valign=top>&nbsp&nbsp&nbsp&nbsp
<INPUT type="submit" name="saveTheme" value="Save" style="font-size:15px;width:60px;">
<INPUT type="button" value="Cancel" style="font-size:15px;width:70px;" onclick="editTheme(-1)">
&nbsp;&nbsp;&nbsp;&nbsp;$updateString
</TD></TR>
</TABLE>
</FORM>

<FORM method="post" name="editListForm" onsubmit="return(checkListForm());">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="id_theme" value="$themeID">
<INPUT type="hidden" name="id_list">
<TABLE border=0 cellspacing=0>
<TR><TH colspan=2 class="title">Lists in Theme <FONT color="#DD0000">$nameTheme</FONT>&nbsp;&nbsp;<INPUT type="button" class="title3" value="Create List from file" onclick="importList();"></TH></TR>
|;
if (scalar keys %{$listList[0]{'LIST'}}==0) {
	print qq
|<TR bgcolor="$color1"><TH colspan=2 class="title2" height=40>No Lists</TH></TR>
</TABLE>
|;
}
else {
	my $disabMove=(scalar keys %{$listList[0]{'LIST'}}==1)? ' disabled' : '';
	print qq
|<TR bgcolor=$color2><TH class="rBorder" valign=middle>
&nbsp;<INPUT type="button" style="width:100px" value="Move up"onclick="moveList(-1)"$disabMove>&nbsp;<BR>
&nbsp;<INPUT type="button" style="width:100px" value="Move down"onclick="moveList(1)"$disabMove>&nbsp;<BR>
</TH><TD>
<TABLE cellspacing=0>
<TR><TH class="rbBorder"><FONT class="title3">Name</FONT></TH><TH class="rbBorder" style="width:100px"><FONT class="title3">Type</FONT><TH class="rbBorder" style="min-width:100px"><FONT class="title3">Size</FONT></TH><TH class="bBorder" width=550 align=left><FONT class="title3">&nbsp;Description</FONT></TH></TR>
|;
	my $bgColor=$color1;
	foreach my $listID (sort{$listList[0]{'LIST'}{$a}[2] <=> $listList[0]{'LIST'}{$b}[2]} keys %{$listList[0]{'LIST'}}) {
		my $isDeletable=($notDelLists{$listID})? 0 : 1;
		my ($typeStrg,$sizeStrg)=($listType{$listID} eq 'SITE')? ('Sites',$listSize{$listID}[1].' / '.$listSize{$listID}[0]) : ('Proteins',$listSize{$listID}[0]);
		print qq
|<TR id="list1_$listID" bgcolor="$bgColor" class="list" height=25>
<TH align="left" nowrap>&nbsp;<A id="$listID" href="javascript:selectList($listID,$isDeletable);">$listList[0]{LIST}{$listID}[0]</A>&nbsp;</TH>
<TH>$typeStrg</TH>
<TH>$sizeStrg</TH>
<TH align="left">&nbsp;$listList[0]{LIST}{$listID}[1]</TH>
</TR>
<TR id="list2_$listID" bgcolor="$bgColor" class="list" height=25 style="display:none;">
<TD>&nbsp;<INPUT type="text" size=35 name="name_$listID" id="name_$listID" value="$listList[0]{LIST}{$listID}[0]" maxlength=50></TD>
<TH>$typeStrg</TH>
<TH>$sizeStrg</TH>
<TD>&nbsp;<INPUT type="text" size=65 name="des_$listID" value="$listList[0]{LIST}{$listID}[1]" maxlength=100></TD>
</TR>
|;
		$bgColor=($bgColor eq $color1)? $color2 : $color1;
	}
	print qq
|</TABLE></TD></TR>
<TR id="editListButton" bgcolor=$color2 height=30><TH colspan=2>
<INPUT type="button" value="Edit" style="width:60px;" onclick="editList(1);">
&nbsp;<INPUT type="button" id="delButton" value="Delete List" onclick="deleteList();">
&nbsp;<INPUT type="button" id="useButton" value="Usage" onclick="usageList();"></TH></TR>
<TR id="saveListButton" bgcolor=$color2 style="display:none;" height=30><TH colspan=2>
<INPUT type="submit" name="saveList" value="Save" style="width:60px;">
&nbsp;<INPUT type="button" value="Cancel" style="width:60px;" onclick="editList(-1);">
</TD></TR>
</TABLE>
|;
}
print qq
|<DIV id="displayList" style="display:none"></DIV>
</FORM>
</CENTER>
<BR>
<BR>
</BODY>
</HTML>
|;

sub ajaxUsageList {

	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	my $listID=param('usageList');

	##DISPLAY / NOT DELETE
	my $sthListComp=$dbh->prepare("SELECT ID_COMPARISON,NAME FROM COMPARISON WHERE ID_CATEGORY=$listID"); # COMPARISON
	my $sthCatComp=$dbh->prepare("SELECT ID_COMPARISON FROM CAT_COMPARISON WHERE ID_CATEGORY=$listID"); # CAT_COMPARISON
	my $sthListPA=$dbh->prepare("SELECT PA.ID_PATHWAY_ANALYSIS,PA.NAME,E.NAME FROM PATHWAY_ANALYSIS PA INNER JOIN EXPERIMENT E ON E.ID_EXPERIMENT=PA.ID_EXPERIMENT WHERE ID_CATEGORY=$listID"); # PATHWAY_ANALYSIS
	my $sthListEA=$dbh->prepare("SELECT EA.ID_EXPLORANALYSIS,EA.CAT_EXCLUSION,EA.NAME,EA.ANA_TYPE,E.NAME FROM EXPLORANALYSIS EA INNER JOIN EXPERIMENT E ON E.ID_EXPERIMENT=EA.ID_EXPERIMENT WHERE ID_CATEGORY=$listID"); # EXPLORANALYSIS

	##DISPLAY / UPDATE /DELETE
	#my $sthListGO=$dbh->prepare("SELECT CA.PARAM_STRG FROM CATEGORY CA,GO_ANALYSIS GA WHERE GA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID"); # GO_ANALYSIS
	my $sthGOParam=$dbh->prepare("SELECT GA.ID_GOANALYSIS,GA.NAME,GA.PARAM_STRG,E.NAME FROM EXPERIMENT E
				      INNER JOIN GO_ANALYSIS GA on E.ID_EXPERIMENT=GA.ID_EXPERIMENT
			              where E.ID_PROJECT=$projectID and GOA_TYPE='graph'");
	my $sthExploParam=$dbh->prepare("SELECT EA.NAME,EA.ANA_TYPE,A.ID_ANNOTATIONSET,A.NAME,A.ANNOT_TYPE,A.ANNOT_LIST,A.RANK,E.NAME FROM EXPERIMENT E
				         INNER JOIN EXPLORANALYSIS EA ON E.ID_EXPERIMENT=EA.ID_EXPERIMENT
				         INNER JOIN ANNOTATIONSET A ON EA.ID_EXPLORANALYSIS=A.ID_EXPLORANALYSIS
				         WHERE E.ID_PROJECT=$projectID AND ANNOT_TYPE LIKE 'prot:%' ORDER BY A.ANNOT_LIST");
	my ($catName)=$dbh->selectrow_array("SELECT NAME FROM CATEGORY WHERE ID_CATEGORY=$listID");
	##Search in comparison/cat_comparison
	my %listComp;
	$sthListComp->execute;
	while (my($compID, $compName)=$sthListComp->fetchrow_array) {
		$listComp{$compID}=$compName;
	}
	$sthListComp->finish;
	my %listCatComp;
	$sthCatComp->execute;
	while (my ($compID)=$sthCatComp->fetchrow_array) {
		my ($compName)=$dbh->selectrow_array("SELECT NAME FROM COMPARISON WHERE ID_COMPARISON=$compID");
		$listCatComp{$compID}=$compName;
	}
	$sthCatComp->finish;
	##Search in Pathway Analysis
	my %listPA;
	$sthListPA->execute;
	while (my ($paID,$paName,$expName)=$sthListPA->fetchrow_array) {
		@{$listPA{$paID}}=($expName,$paName);
	}
	$sthListPA->finish;
	##Search in GO Analysis
	my %listGO;
	$sthGOParam->execute;
	while (my ($goID,$goName,$paramStrg,$expName)=$sthGOParam->fetchrow_array) {
		#next if ($paramStrg !~ /protSetCategory/);
		foreach my $param(split(/;/,$paramStrg)) {
			if ($param=~/protSetCategory/) {
				my ($name,$catID)=split(/=/,$param);
				if ($catID == $listID) {
					@{$listGO{$goID}}=($expName,$goName,'prot');
				}
			}
			if ($param=~/popCategory/) {
				my ($name,$catID)=split(/=/,$param);
				if ($catID == $listID) {
					@{$listGO{$goID}}=($expName,$goName,'pop'); # if ($catID != "0");
				}
			}
		}
	}
	$sthGOParam->finish;
	## Search in ExplorAnalysis
	my %listEA;
	$sthListEA->execute;
	while (my($explorID, $exclusion, $explorName,$explorType,$expName)=$sthListEA->fetchrow_array) {
		$exclusion=0 unless $exclusion;
		@{$listEA{$explorType}{$explorID}}=($expName,$explorName,$exclusion);
	}
	$sthListEA->finish;
	## Search in Annotation Set
	my %listAnnot;
	$sthExploParam->execute;
	while (my ($explorName,$anaType,$annotID,$annotName,$annotType,$annotList,$rank,$expName)=$sthExploParam->fetchrow_array) {
		if ($anaType eq 'PCA') {#PCA
			if ($annotType eq "prot:LIST") {
				$annotName=~s/#//g;
				if ($annotName == $listID) {
					@{$listAnnot{$anaType}{$annotID}}=($expName,$explorName,$annotName,$rank);
				}
			}
		}
		elsif ($anaType eq 'cluster') {
			if ($annotType eq "prot:LIST") {
				foreach my $listStrg (split(':',$annotList)) {
					$listStrg=~s/#//g;
					if ($listStrg == $listID) {
						@{$listAnnot{$anaType}{$annotID}}=($expName,$explorName,$annotName,$rank); # $annotName."~~".$explorName."~~".$rank;
					}
				}
			}
			elsif ($annotType eq "prot:THEME") {
				$annotName=~s/#//g;
				if ($annotName == $themeID) {
					@{$listAnnot{$anaType}{$annotID}}=($expName,$explorName,'All lists of Theme',$rank); # "All lists~~".$explorName."~~".$rank;
				}
			}
		}
	}
	$sthExploParam->finish;

	my ($color1,$color2)=&promsConfig::getRowColors;
	if (keys %listComp || keys %listCatComp || keys %listPA || keys %listEA || keys %listGO || keys %listAnnot) {
		print qq
|<BR>
<FONT class="title3">The List <FONT color="#DD0000">$catName</FONT> is involved in the following data analyses</FONT>
<TABLE bgcolor="$color2">
|;
	}
	my $isDel=0;
	##Analyses to display not to delete
	if (keys %listComp || keys %listCatComp || keys %listPA || keys %listEA) {
		$isDel=1;
		if (keys %listComp) {
			print qq|<TR><TH align="right" valign="top">&nbsp;Comparison Filters:</TH><TD bgcolor=$color1>|;
			foreach my $compID (keys %listComp) {
				print "&nbsp;&nbsp;-&nbsp;$listComp{$compID}<br>\n";
			}
			print "</TD></TR>\n";
		}
		if (keys %listCatComp) {
			print qq|<TR><TH align="right" valign="top">&nbsp;Comparison Elements:</TH><TD bgcolor=$color1>|;
			foreach my $compID (keys %listCatComp) {
				print "&nbsp;&nbsp;-&nbsp;$listCatComp{$compID}<br>\n";
			}
			print "</TD></TR>\n";
		}
		if (keys %listPA) {
			print qq|<TR><TH align="right" valign="top">&nbsp;Pathway Analyses:</TH><TD bgcolor=$color1>|;
			foreach my $paID (keys %listPA) {
				#my ($expName)=$dbh->selectrow_array("SELECT E.NAME FROM PATHWAY_ANALYSIS PA
				#				     INNER JOIN EXPERIMENT E on PA.ID_EXPERIMENT=E.ID_EXPERIMENT
				#				     where PA.ID_PATHWAY_ANALYSIS=$paID");
				my ($expName,$paName)=@{$listPA{$paID}};
				print "&nbsp;&nbsp;-&nbsp;$expName > $paName<BR>\n";
			}
			print "</TD></TR>\n";
		}
		if (keys %listEA) {
			print qq|<TR ><TH align="right" valign="top">&nbsp;Exploratory Analyses :</TH><TD bgcolor=$color1 nowrap>|;
			foreach my $type (keys %listEA) {
				my $typeStrg=($type eq 'cluster')? 'Clustering' : ($type eq 'MOTIF')? 'Motif search' : ($type eq 'HM_MOTIF')? 'Motif clustering' : $type;
				print "&bull;<B>$typeStrg:</B><BR>\n";
				foreach my $explorID (keys %{$listEA{$type}}) {
					my ($expName,$anaName,$exclu)=@{$listEA{$type}{$explorID}};
					my $excluStrg = ($exclu == 1)? " (exclusion list)" : '';
					print "&nbsp;&nbsp;-&nbsp;$expName > $anaName$excluStrg&nbsp;<BR>\n";
				}
			}
			print "</TD></TR>\n";
		}
	}

	if (keys %listGO || keys %listAnnot) {
		print qq|<INPUT type="button" value="Force Delete" onclick="deleteList('force')">| if (!$isDel);
		print qq|<TR><TH align="right" valign="top">&nbsp;Annotations :</TH><TD bgcolor=$color1>|;
		if (keys %listGO) {
			print "&bull;&nbsp;<B>GO:</B><BR>";
			foreach my $goID (keys %listGO) {
				#my ($expName)=$dbh->selectrow_array("SELECT E.NAME FROM GO_ANALYSIS G, EXPERIMENT E where G.ID_EXPERIMENT=E.ID_EXPERIMENT and G.ID_GOANALYSIS=$goID");
				my ($expName,$goName,$item)=@{$listGO{$goID}};
				print qq |<INPUT type="hidden" name="GO" value="$goID-$item">&nbsp;&nbsp;-&nbsp;$expName > $goName<BR>\n|;
			}
		}

		if (keys %listAnnot) {
			foreach my $type (sort{$a cmp $b} keys %listAnnot) { # cluster or PCA
				my $typeStrg=($type eq 'cluster')? 'Clustering' : $type; # ($type eq 'MOTIF')? 'Motif search' : ($type eq 'HM_MOTIF')? 'Motif clustering' :
				print "&bull;<B>$typeStrg:</B><BR>";
				foreach my $annotID (keys %{$listAnnot{$type}}) {
					#my ($expName)=$dbh->selectrow_array("SELECT E.NAME FROM ANNOTATIONSET A
					#				     INNER JOIN EXPLORANALYSIS EA on A.ID_EXPLORANALYSIS=EA.ID_EXPLORANALYSIS
					#				     INNER JOIN EXPERIMENT E on EA.ID_EXPERIMENT=E.ID_EXPERIMENT
					#				     where A.ID_ANNOTATIONSET=$annotID");
					#my ($anName,$eaName,$rank)=split(/~~/,$listAnnot{$type}{$annotID});
					my ($expName,$eaName,$anName,$rank)=@{$listAnnot{$type}{$annotID}};
					if ($anName=~/All/) {
						if ($isDel) {
							print "&nbsp;&nbsp;-&nbsp;$anName involved in $eaName (not deletable)<br>";
						}
						else {
							print qq|<INPUT type="hidden" name="AN" value="$type-THEME-$annotID-$rank">&nbsp;&nbsp;-&nbsp;$expName > $eaName > $anName<BR>\n|;
						}
					}
					else {
						print qq|<INPUT type="hidden" name="AN" value="$type-LIST-$annotID-$rank">&nbsp;&nbsp;-&nbsp;$expName > $eaName > $anName<BR>\n|;
					}
				}
			}

		}
		print "</TD></TR>\n";
	}
	print "</TABLE>\n" if (keys %listComp || keys %listCatComp || keys %listPA || keys %listGO || keys %listEA || keys %listAnnot);
	$dbh->disconnect;
	exit;
}


sub importList {
	my $process=$_[0];

	#####################
	####<Import form>####
	#####################
	if ($process==1) {
		my ($nameTheme)=$dbh->selectrow_array("SELECT NAME FROM CLASSIFICATION WHERE ID_CLASSIFICATION=$themeID");
		my %identifierList;
		my $sthID=$dbh->prepare("SELECT ID_IDENTIFIER,NAME FROM IDENTIFIER");
		$sthID->execute;
		while (my($identID,$identName)=$sthID->fetchrow_array) {
			$identifierList{$identID}=$identName;
		}
		$sthID->finish;
		my %speciesList;
		my $sthSp=$dbh->prepare("SELECT ID_SPECIES,SCIENTIFIC_NAME,COMMON_NAME FROM SPECIES WHERE IS_REFERENCE=1");
		$sthSp->execute;
		while (my($speciesID,$scientName,$comName)=$sthSp->fetchrow_array) {
			$speciesList{$speciesID}=$scientName." ($comName)";
		}
		$sthSp->finish;

		####<Starting HTML>####
		my ($color1,$color2)=&promsConfig::getRowColors;
		print header(-'content-encoding'=>'no',-charset=>'utf-8');
		warningsToBrowser(1); # writes PERL warnings as HTML comments <!-- --> must be activated after header
		print qq
|<HEAD>
<TITLE>Theme : $nameTheme</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
var prevFileType='';
function updateForm(fileType) {
	var myForm=document.importListForm;
	if (fileType=='gmt') {
		myForm.listName.style.display='none';
		document.getElementById('nameSpan').style.display='';
		myForm.listDes.style.display='none';
		document.getElementById('desSpan').style.display='';
		document.getElementById('mergeSpan').style.display='';
	}
	else if (prevFileType=='gmt') {
		myForm.listName.style.display='';
		document.getElementById('nameSpan').style.display='none';
		myForm.listDes.style.display='';
		document.getElementById('desSpan').style.display='none';
		document.getElementById('mergeSpan').style.display='none';
	}
	prevFileType=fileType;
}
function checkListForm(myForm){
	if (!myForm.species.value) {
		alert("Select the species of origin.");
		return false;
	}
	if (!myForm.importSource.value) {
		alert("Select an import source.");
		return false;
	}
	if (myForm.importSource.value=='file') {
		if (!myForm.fileFormat.value) {
			alert("Select a file format.");
			return false;
		}
		if (myForm.fileFormat.value != 'gmt') {
			if (!myForm.listName.value) {
				alert("Provide a name for the List.");
				return false;
			}
		}
		if (!myForm.importFile.value) {
			alert("Provide a file to be imported.");
			return false;
		}
	}
	else if (!myForm.importText.value){ // pasted text
		alert("Paste a list of identifiers in the text box.");
		return false;
	}
	return true;
}
</SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'">
<CENTER>
<FONT class="title">Creating New List(s) in Theme <FONT color="#DD0000">$nameTheme</FONT> from File</FONT>
<FORM method="post" name="importListForm" enctype="multipart/form-data" onsubmit="return(checkListForm(this));" >
<INPUT type="hidden" name="id_project" value="$projectID"/>
<INPUT type="hidden" name="id_theme" value="$themeID"/>
<INPUT type="hidden" name="import" value="2"/>
<BR><BR>
<TABLE bgcolor="$color2">
<TR><TH align="right">&nbsp;Name :</TH><TD bgcolor="$color1"><INPUT type="text" name="listName" size="50" maxlength="50"/><SPAN id="nameSpan" style="display:none">&nbsp;From file</SPAN></TD></TR>
<TR><TH align="right" valign="top">&nbsp;Description :</TH><TD bgcolor="$color1"><TEXTAREA name="listDes" maxlength="100" cols="65" rows="2"></TEXTAREA><SPAN id="desSpan" style="display:none">&nbsp;From file</SPAN></TD></TR>
<TR><TH align="right" valign="top" nowrap>&nbsp;Identifiers source:</TH><TD bgcolor="$color1" valign="top" nowrap>
&bull;<B>Identifier type:</B><SELECT name="identifier"><OPTION value="">* Unspecified *</OPTION>
|;
		foreach my $identID (sort{lc($identifierList{$a}) cmp lc($identifierList{$b})} keys %identifierList) {
			print "<OPTION value=\"$identID\">$identifierList{$identID}</OPTION>\n";
		}
		print qq
|</SELECT>
&nbsp;&nbsp;&bull;<B>Species:</B><SELECT name="species">
	<OPTION value="">-= Select =-</OPTION>
	<OPTION value="-1">No species/Multi species</OPTION>
|;
		foreach my $speciesID (sort{lc($speciesList{$a}) cmp lc($speciesList{$b})} keys %speciesList) {
			print "<OPTION value=\"$speciesID\">$speciesList{$speciesID}</OPTION>\n";
		}
		print qq
|</SELECT><BR>
<FIELDSET style="padding:2px;white-space:nowrap"><LEGEND><INPUT type="radio" name="importSource" value="file"><B>Text file:</B></LEGEND>
&bull;<B>Format:</B><SELECT name="fileFormat" onchange="updateForm(this.value)"><OPTION value="">-= Select =-</OPTION><OPTION value="gmt">gmt</OPTION><OPTGROUP label="Custom:"><OPTION value="column">single column</OPTION><OPTION value="line">single line</OPTION></OPTGROUP></SELECT>
<SPAN id="mergeSpan" style="display:none">&nbsp;&nbsp;&bull;<INPUT type="checkbox" name="mergeGmt" value="1"><B>Merge all signatures into a single List</B></SPAN><BR>
&bull;<B>Choose file:</B><INPUT type="file" name="importFile" size="50">
</FIELDSET>
<B>Or</B><BR>
<FIELDSET style="padding:2px;white-space:nowrap"><LEGEND><INPUT type="radio" name="importSource" value="text"><B>Paste text:</B></LEGEND>
<TABLE cellpadding=0><TR>
<TD valign="top"><TEXTAREA name="importText" cols="20" rows="10"></TEXTAREA></TD><TD valign="top" class="font11"><B>Separator accepted:</B><BR>&bull; colon ','<BR>&bull; semi-colon ';'<BR>&bull; space ' '<BR>&bull; tab '\\t'<BR>&bull; new line '\\n'</TD>
</TR></TABLE>
</FIELDSET>
</TD></TR>
<TR><TH colspan=2><INPUT type="submit" value="Proceed"></TH></TR>
</TABLE>
</FORM>
</BODY>
<HTML>
|;

	}
	############################
	####<Submit import form>####
	############################
	else {
		my $listName=param('listName'); # ignored if file format is gmt
		my $listDes=param('listDes'); # ignored if file format is gmt
		my $identifierID=param('identifier');
		my $specieID=param('species');
		my $importSource=param('importSource');
		#>File import
		my $fileFormat=param('fileFormat');
		my $mergeList=param('mergeGmt') || 0;
		my $fileName = tmpFileName(upload('importFile')) if upload('importFile');
		#>Text import
		my $importText=param('importText');

		####<Starting HTML>####
		my ($color1,$color2)=&promsConfig::getRowColors;
		print header(-'content-encoding'=>'no',-charset=>'utf-8');
		warningsToBrowser(1); # writes PERL warnings as HTML comments <!-- --> must be activated after header
		print qq
|<HEAD>
<TITLE>Importing List(s) from File</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
function showIdentifiers(what,show) {
	var showStatus,hideStatus='';
	if (show) {
		showStatus='';
		hideStatus='none';
	}
	else {
		showStatus='none';
		hideStatus='';
	}
	if (what=='badVisDiv') {
		document.getElementById('badVisDiv').style.display=showStatus;
		document.getElementById('badVisBut1').style.display=hideStatus;
		document.getElementById('badVisBut2').style.display=showStatus;
	}
	else if (what=='notMappedDiv') {
		document.getElementById('notMappedDiv').style.display=showStatus;
		document.getElementById('notMappedBut1').style.display=hideStatus;
		document.getElementById('notMappedBut2').style.display=showStatus;
	}
}
</SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'">
<CENTER>
<FONT class="title">Importing List File</FONT>
<BR>
<IMG src="$promsPath{images}/engrenage.gif">
</CENTER>
<BR>
<FONT class="title2">
<BR>Reading file...
|;
		####<Reading file>####
		my (%listInfo,%identifierList,%allIdentifiers);
		if ($importSource eq 'file') {
			open(FILE,$fileName) or die $!;
			if ($fileFormat eq 'gmt') {
				my $listCount=0;
				while (<FILE>) {
					chomp;
					@{$listInfo{0}}=($listName,$listDes) if $mergeList;
					my ($name,$des,@identifers)=split(/\t+/,$_);
					if ($name) {
						unless ($mergeList) {
							$listCount++ ;
							$des=~s/^> //;
							@{$listInfo{$listCount}}=($name,$des);
						}
						foreach my $ident (@identifers) {
							$ident=~s/\s+\Z//;
							next unless $ident=~/\w/;
							$identifierList{$listCount}{$ident}=1;
							$allIdentifiers{$ident}=1;
						}
					}
				}
			}
			elsif ($fileFormat eq 'column') {
				@{$listInfo{0}}=($listName,$listDes);
				while (<FILE>) {
					next if $_=~/^#/;
					chomp;
					$_=~s/"//g;
					$_=~s/^\s+//g;
					$_=~s/\s+\Z//g;
					my ($ident)=($_=~/^(\S+)/);
					next unless $ident=~/\w/;
					$identifierList{0}{$ident}=1;
					$allIdentifiers{$ident}=1;
				}
			}
			elsif ($fileFormat eq 'line') {
				@{$listInfo{0}}=($listName,$listDes);
				while (<FILE>) {
					next if $_=~/^#/;
					chomp;
					$_=~s/"//g;
					$_=~s/^\s+//g;
					$_=~s/\s+\Z//g;
					next unless $_=~/\S/;
					foreach my $ident (split(/[\s+,;]/)) {
						next unless $ident=~/\w/;
						$identifierList{0}{$ident}=1;
						$allIdentifiers{$ident}=1;
					}
				}
			}
			close FILE;
		}
		else { # text import
			@{$listInfo{0}}=($listName,$listDes);
			foreach my $line (split(/\n/,$importText)) {
				$line=~s/"//g;
				$line=~s/^\s+//g;
				$line=~s/\s+\Z//g;
				foreach my $ident (split(/[\s+,;]/,$line)) {
					next unless $ident=~/\w/;
					$identifierList{0}{$ident}=1;
					$allIdentifiers{$ident}=1;
				}
			}
		}
		my $numLists=scalar keys %identifierList;
		my $numIdentifiers=scalar keys %allIdentifiers;
		my $listStrg=($numLists > 1)? 'Lists' : 'List';
		print qq
| Done ($numIdentifiers identifiers found in $numLists $listStrg).<BR>
<BR>Mapping identifiers: </FONT><B>0%|;

		####<Mapping identifiers>####
		my (%mappedIdentifiers,%badVisIdentifiers,%proteinList,%allProteins,%notMappedIdentifiers,%usedIdentifiers);
		my $addProtQuery=($identifierID)?" AND ID_IDENTIFIER=$identifierID":"";
		my $addSpeciesQuery=($specieID > -1) ? " AND MP.ID_SPECIES=$specieID" : ""; 
		my $sthProt=$dbh->prepare("SELECT P.ID_PROTEIN, MAX(AP.VISIBILITY) FROM PROTEIN P LEFT JOIN MASTER_PROTEIN MP ON P.ID_MASTER_PROTEIN = MP.ID_MASTER_PROTEIN$addSpeciesQuery LEFT JOIN MASTERPROT_IDENTIFIER MI ON MP.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN LEFT JOIN ANALYSIS_PROTEIN AP ON P.ID_PROTEIN=AP.ID_PROTEIN$addProtQuery WHERE P.ID_PROJECT=$projectID AND (MI.VALUE=? OR P.ALIAS=? OR P.IDENTIFIER=?) GROUP BY P.ID_PROTEIN");
		my @percent=(10,20,30,40,50,60,70,80,90,100);
		my @limitValue;
		foreach my $pc (@percent) {push @limitValue,int(0.5+($numIdentifiers*$pc/100));}
		my $pcIdx=0;
		my $count1=0;
		my $count2=0;
		my $maxCount2=int(0.5+($numIdentifiers/100));
		foreach my $listCount (keys %identifierList) {
			foreach my $ident (keys %{$identifierList{$listCount}}) {
				if ($usedIdentifiers{$ident}) { # already tested
					if ($mappedIdentifiers{$ident}) { # already mapped => update list contents with premapped proteins
						foreach my $protID (@{$mappedIdentifiers{$ident}}) {$proteinList{$listCount}{$protID}=1;}
					}
					next;
				}
				$sthProt->execute($ident, $ident, $ident);
				my $mapped=0;
				while (my ($protID,$bestVis)=$sthProt->fetchrow_array) {
					$mapped=1;
					#if ($bestVis >= 1) { !!! Commented to allow hidden proteins in List !!!
						push @{$mappedIdentifiers{$ident}},$protID;
						$proteinList{$listCount}{$protID}=1; # proteins found in current list
						$allProteins{$protID}=1;
					#}
					#else {
						$badVisIdentifiers{$ident}=1 if $bestVis==0;
					#}
				}
				$notMappedIdentifiers{$ident}=1 unless $mapped;
				$usedIdentifiers{$ident}=1;
				$count1++;
				if ($count1>=$limitValue[$pcIdx]) {
					print "$percent[$pcIdx]%";
					$pcIdx++;
				}
				$count2++;
				if ($count2==$maxCount2) {print '.'; $count2=0;}
			}
		}
		print '100%' if $pcIdx<10; # Just in case
		$sthProt->finish;
		my $numGood=scalar keys %mappedIdentifiers;
		my $numBadVis=scalar keys %badVisIdentifiers;
		my $numProt=scalar keys %allProteins;
		my $notMapped=scalar keys %notMappedIdentifiers;
		print qq
|</B><FONT class="title2"> Done.
<BR>&nbsp;&nbsp;&nbsp;-$numGood identifiers were mapped to $numProt proteins.
|;
		if ($numBadVis) {
			print qq
|<BR>&nbsp;&nbsp;&nbsp;-$numBadVis identifiers mapped to hidden proteins. <INPUT type="button" id="badVisBut1" value=" Show " onclick="showIdentifiers('badVisDiv',true)"/><INPUT type="button" id="badVisBut2" value=" Hide " onclick="showIdentifiers('badVisDiv',false)" style="display:none"/>
<DIV id="badVisDiv" class="title3" style="display:none">
|;
			print join(', ',sort{lc($a) cmp lc($b)} keys %badVisIdentifiers);
			print "</DIV>\n";
		}
		if ($notMapped) {
			print qq
|<BR>&nbsp;&nbsp;&nbsp;-$notMapped identifiers could not be mapped. <INPUT type="button" id="notMappedBut1" value=" Show " onclick="showIdentifiers('notMappedDiv',true)"/><INPUT type="button" id="notMappedBut2" value=" Hide " onclick="showIdentifiers('notMappedDiv',false)" style="display:none"/>
<DIV id="notMappedDiv" class="title3" style="display:none">
|;
			print join(', ',sort{lc($a) cmp lc($b)} keys %notMappedIdentifiers);
			print "</DIV>\n";
		}

		####<Storing results in DBs>####
		print qq
|<BR><BR>Storing results...|;
		my ($catID)=$dbh->selectrow_array("SELECT MAX(ID_CATEGORY) FROM CATEGORY");
		my ($displayPos)=$dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM CATEGORY WHERE ID_CLASSIFICATION=$themeID");
		my $sthInsCat=$dbh->prepare("INSERT INTO CATEGORY (ID_CATEGORY,ID_CLASSIFICATION,NAME,DES,DISPLAY_POS,LIST_TYPE,UPDATE_DATE,UPDATE_USER) VALUES (?,?,?,?,?,'PROT',NOW(),?)");
		my $sthInsProt=$dbh->prepare("INSERT INTO CATEGORY_PROTEIN (ID_CATEGORY,ID_PROTEIN) VALUES (?,?)");
		foreach my $listCount (sort{$a<=>$b} keys %listInfo) {
			$sthInsCat->execute(++$catID,$themeID,@{$listInfo{$listCount}},++$displayPos,$userID);
			#print "@{$listInfo{$listCount}}<BR>\n";
			foreach my $protID (keys %{$proteinList{$listCount}}) {
				$sthInsProt->execute($catID,$protID);
				#print "$protID<BR>\n";
			}
		}
		$sthInsCat->finish;
		$sthInsProt->finish;
		$dbh->commit;
		print qq
| Done.</FONT>
<BR><BR><INPUT type="button" class="title2" value=" Display Theme " onclick="window.location='./editClassification.cgi?id_project=$projectID&id_theme=$themeID'"/>
<BR><BR>
</BODY>
</HTML>
|;
	}

	$dbh->disconnect;
	exit;
}

####>Revision history<####
# 2.3.1 Add No/Multi species + Fix interactions with theme having no match (VS 14/02/2019)
# 2.3.0 Compatible with MODIFICATION_SITE & code optimization (PP 01/08/17)
# 2.2.2 Add unspecified identifier (GA 24/08/16)
# 2.2.1 Accepts also a pasted list of identifiers during import (PP 04/03/16)
# 2.2.0 Added list size & creation from imported file (PP 25/08/15)
# 2.1.0 add force list deletion (SL 15/06/15)
# 2.0.9 Extends non-deletability to list used in GO_ANALYSIS (PP 15/05/14)
# 2.0.8 Update for List comparison (PP 21/01/14)
# 2.0.7 Minor undef value correction (PP 29/11/13)
# 2.0.6 GPL license (PP 23/09/13)
# 2.0.5 ' and " not allowed due to interference with Javascript (PP 22/10/12)
# 2.0.4 Convert Classification to Theme & Category to List (PP 19/07/12)
# 2.0.3 Control for deletable Categories (PP 09/03/2011)
