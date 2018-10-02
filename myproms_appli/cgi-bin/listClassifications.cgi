#!/usr/local/bin/perl -w

################################################################################
# listClassifications.cgi                2.1.3                                 #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Lists Themes (former Classifications)                                        #
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

####################
#    Parameters    #
####################
my $projectID=param('id_project');

###############
#     Main    #
###############

# Connect to the database
my $dbh=&promsConfig::dbConnect;

if (param('save')) {&addTheme;}
elsif (param('delete')) {&deleteTheme;}

my ($projectName)=$dbh->selectrow_array("SELECT NAME FROM PROJECT WHERE ID_PROJECT=$projectID");

my %listTheme=&promsMod::getListClass($dbh,$projectID);

##>Deletability<## WARNING: Does not check for Lists used PARAMS of GO and Exploratory Analyses!!! # TODO: Add check or disable Delete option unless Theme is empty
my @sthDisab=(
	$dbh->prepare("SELECT 1 FROM CATEGORY CA,COMPARISON CO WHERE CO.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=? LIMIT 1"), # COMPARISON
	$dbh->prepare("SELECT 1 FROM CATEGORY CA,CAT_COMPARISON CC WHERE CC.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=? LIMIT 1"), # CAT_COMPARISON
	$dbh->prepare("SELECT 1 FROM CATEGORY CA,GO_ANALYSIS GA WHERE GA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=? LIMIT 1"), # GO_ANALYSIS
	$dbh->prepare("SELECT 1 FROM CATEGORY CA,EXPLORANALYSIS EA WHERE EA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=? LIMIT 1") # EXPLORANALYSIS
);

foreach my $themeID (keys %listTheme) {
	foreach my $sth (@sthDisab) {
		$sth->execute($themeID);
		if ($sth->fetchrow_array) {
			$listTheme{$themeID}[2]=1; # adds a 3rd index to array
			last;
		}
	}
}
foreach my $sth (@sthDisab) {$sth->finish;}

$dbh->disconnect;

###########################
#      Starting HTML      #
###########################
my ($color1,$color2)=&promsConfig::getRowColors;
my $bgColor=$color1;
print header(-'content-encoding'=>'no',-charset=>'utf-8');
warningsToBrowser(1); # writes PERL warnings as HTML comments <!-- --> must be activated after header
print qq
|<HEAD>
<TITLE>Manage of Themes</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
|;
&promsMod::popupInfo();
print qq
|function listThemeProteins(id_theme) {
	parent.selectedMode='classification:'+id_theme; // update promsFrame variable
	parent.selectedAction='proteins';
	parent.optionFrame.selectOption();
}
function editTheme(id_theme){
	window.location="./editClassification.cgi?id_project=$projectID&id_theme="+id_theme;
}
function deleteTheme(id_theme,name_theme){
	if (confirm ("Confirm deletion of Theme '"+name_theme+"'.")) {
		window.location="./listClassifications.cgi?delete=1&id_project=$projectID&id_theme="+id_theme;
	}
}
function addNewTheme(action){
	var vis1; var vis2;
	if (action>=0) {vis1='none';vis2='block';}
	else {vis1='block';vis2='none';}
	document.getElementById('addButton').style.display=vis1;
	document.getElementById('newThemeRow').style.display=vis2;
//	document.getElementById('saveButton').style.display=vis2;
}
function checkForm(){
	if (document.addThemeForm.name.value){
		return true;
	}
	else{
		alert ("Type a name for new Theme.");
		return false;
	}
}
</SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<CENTER>
<FONT class="title">List of Themes in Project <FONT color="#DD0000">$projectName</FONT></FONT>
</CENTER>
<BR>
<FORM method="post" name="addThemeForm" onSubmit="return(checkForm());">
<INPUT type="hidden" name="id_project" value=$projectID>
<TABLE border=0 cellspacing=0>
<TR bgcolor="$color2">
<TH class="rbBorder" width=400><FONT class="title3">Name</FONT></TH>
<TH class="bBorder" colspan=3 width=630 align=left><FONT class="title3">&nbsp;Description</FONT></TH>
</TR>
|;
foreach my $themeID (sort{lc($listTheme{$a}[0]) cmp lc($listTheme{$b}[0])} keys %listTheme){
	print "<TR class=\"list\" bgcolor=\"$bgColor\">\n";
	my $themeName=&promsMod::resize($listTheme{$themeID}[0],40);
	my $popupString=($themeName ne $listTheme{$themeID}[0])? "<B>Name:</B> $listTheme{$themeID}[0]<BR>" : "";
	$popupString.='Click to display protein contents.';
	print "<TH class=\"title3\"><A href=\"javascript:listThemeProteins($themeID);\" onmouseover=\"popup('$popupString')\" onmouseout=\"popout()\">$themeName</A></TH>\n";
	($listTheme{$themeID}[1])=&promsMod::chkDef($listTheme{$themeID}[1]);
	my $disabStrg=($listTheme{$themeID}[2])? ' disabled' : '';
	print qq
|<TH width=480 align=left>$listTheme{$themeID}[1]</TH>
<TH width=150 nowrap><INPUT type=\"button\" value=\"Manage\" onclick=\"editTheme($themeID)\">
<INPUT type=\"button\" value=\"Delete\" onclick=\"deleteTheme($themeID,'$listTheme{$themeID}[0]')\"$disabStrg></TH>
</TR>
|;
	$bgColor=($bgColor eq $color1)? $color2 : $color1;
}
if (scalar keys %listTheme==0) {
	print "<TR bgcolor=\"$color1\"><TH class=\"title2\" colspan=4>No Themes</TH></TR>\n";
	$bgColor=$color2;
}
print qq
|</TABLE>

<TABLE id="newThemeRow" border=0 cellspacing=0 style="display:none;">
<TR bgcolor=$bgColor>
<TD width=402>&nbsp<INPUT type="text" style="width:390px" name="name" value="" maxlength=50></TD>
<TD width=630>&nbsp<INPUT type="text" style="width:590px" name="des" value="" maxlength=100></TD>
</TR>
<TR height=30 bgcolor=$color2>
<TD align=center colspan=2><INPUT type="submit" name="save" value="Save" style="width:70px">
&nbsp;&nbsp;<INPUT type="button" value="Cancel" style="width:70px" onclick="addNewTheme(-1)"></TD>
</TR>
</TABLE>

<TABLE id="addButton" border=0 cellspacing=0>
<TR height=30 bgcolor=$color2>
<TD align=center width=1034><INPUT type="button" value="Add a new Theme" style="width:150px" onclick="addNewTheme(1)"></TD>
</TR>
</TABLE>

</FORM>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup()</SCRIPT>
</BODY>
</HTML>
|;


##############<SUBROUTINES>#############

####################
####<<<CREATE>>>####
####################
sub addTheme {
	my ($themeID)=$dbh->selectrow_array("SELECT Max(ID_CLASSIFICATION) FROM CLASSIFICATION");
	$themeID++;
	my $name=param('name');
	my $description=param('des');
	my $date=strftime("%Y-%m-%d %H:%M:%S",localtime);
	$name=$dbh->quote($name);
	$description=$dbh->quote($description);

	##<Insert the new theme in DB
	$dbh->do("INSERT INTO CLASSIFICATION (ID_CLASSIFICATION,ID_PROJECT,NAME,DES,UPDATE_DATE,UPDATE_USER) VALUES ($themeID,$projectID,$name,$description,'$date','$userID')");
	$dbh->commit; # not necessary with 'do'
}


####################
####<<<DELETE>>>####
####################
sub deleteTheme {
	my $themeID=param('id_theme');
	my $sthSelCat=$dbh->prepare("SELECT ID_CATEGORY FROM CATEGORY WHERE ID_CLASSIFICATION=$themeID");
	my $sthDelMS=$dbh->prepare("DELETE FROM MODIFICATION_SITE WHERE ID_CATEGORY=?");
	my $sthDelCat=$dbh->prepare("DELETE FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=?");
	$sthSelCat->execute();
	while (my ($idCat)=$sthSelCat->fetchrow_array) {
		$sthDelMS->execute($idCat);
		$sthDelCat->execute($idCat);
	}
	$sthSelCat->finish;
	$sthDelMS->finish;
	$sthDelCat->finish;
	$dbh->do("DELETE FROM CATEGORY WHERE ID_CLASSIFICATION=$themeID");
	$dbh->do("DELETE FROM CLASSIFICATION WHERE ID_CLASSIFICATION=$themeID");
	$dbh->commit;
}

####>Revision history<####
# 2.1.3 Compatible with MODIFICATION_SITE (PP 01/08/17)
# 2.1.2 Extends Theme non-deletability to list used in EXPLORANALYSIS (PP 26/10/15)
# 2.1.1 Extends Theme non-deletability to list used in GO_ANALYSIS (PP 15/05/14)
# 2.1.0 Update for List comparison (PP 21/01/14)
# 2.0.9 GPL license (PP 23/09/13)
# 2.0.8 Convert Classification to Theme (PP 19/07/12)
# 2.0.7 Rename of class deleteClass into deleteClassification (GA 20/06/2011)
# 2.0.6 Control for deletable Classifications (PP 09/03/2011)
# 2.0.5 new protein list management options (PP 22/02/2011)
