#!/usr/local/bin/perl -w

################################################################################
# manageTemplates.cgi       1.0.7                                              #
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
use POSIX qw(strftime); # to get the time
use strict;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID = $ENV{'REMOTE_USER'};


##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

#############################
my $action= (param('ACT'))? param('ACT') : 'list' ;
my ($userStatus,$userName) = $dbh->selectrow_array("SELECT USER_STATUS,USER_NAME FROM USER_LIST WHERE ID_USER='$userID'");

if (param('submit')) {
	&processForm;
	exit;
}
if ($action eq 'delete'){
	&deleteTemplate;
	exit;
}
elsif ($action eq 'list') {
	### Fetching data ###
	my $dispUserID = (param('DISPLAY_ID') && $userStatus=~/(bioinfo|mass|manag)/)? param('DISPLAY_ID') : $userID;
	my $sthGetTemplates = $dbh->prepare("SELECT ID_VAL_TEMPLATE,NAME,DES,VAL_TYPE,SEARCH_ENGINE,MS_TYPE,IS_DEFAULT FROM VALIDATION_TEMPLATE WHERE ID_USER='$dispUserID' ORDER BY VAL_TYPE DESC,SEARCH_ENGINE,NAME");
	$sthGetTemplates->execute;
	my $refTemplates = $sthGetTemplates->fetchall_arrayref;
	$sthGetTemplates->finish;

	### HTML ###
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	my ($light,$dark)=&promsConfig::getRowColors;
        my %msTypeFullName = &promsConfig::getMsType;
	my $bgColor=$light;
	print qq
|<HTML><HEAD>
<TITLE>List of User Validation Templates</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
function deleteTemplate(id,name){
	if (confirm('Confirm deleting '+name+' ?')){
		window.location = './manageTemplates.cgi?ACT=delete&ID='+id;
	}
}
|;
	if ($userStatus=~/bioinfo|mass|manag/) {
		print qq
|function changeUser(id){
	window.location='./manageTemplates.cgi?ACT=list&DISPLAY_ID='+id;
}
|;
	}
	print qq
|</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<DIV style="float:top">
<BR>
<TABLE><TR><TH bgcolor="$dark">
<FONT class="title2">&nbsp;Go to:</FONT><SELECT class="title action" onchange="if (this.value) {window.location='$promsPath{cgi}/'+this.value;}">
	<OPTION value="">-= Select =-</OPTION>
	<OPTION value="promsMain.cgi">Main Window</OPTION>
	<OPTION value="selectProject.cgi">Project Selection</OPTION>
</SELECT>
</TH></TR></TABLE>
</DIV>
<CENTER>
<FONT class="title">List of Validation Templates</FONT>
|;
	if ($userStatus =~ /bioinfo|mass|manag/) { # bioinfo mass & manag can see all/workgroup users templates
		my $wkgrStrg='';
		if ($userStatus eq 'manag') {
			my ($workgroup)=$dbh->selectrow_array("SELECT WORK_GROUP FROM USER_LIST WHERE ID_USER='$userID'");
			$wkgrStrg=(defined $workgroup)? 'AND WORK_GROUP='.$dbh->quote($workgroup) : 'AND WORK_GROUP IS NULL';
		}
		print "<FONT class=\"title\"> for:</FONT>&nbsp;<SELECT class=\"title\" onchange=\"changeUser(this.value);\">";

		my $sthGetUsers = $dbh->prepare("SELECT DISTINCT(VALIDATION_TEMPLATE.ID_USER),USER_LIST.USER_NAME FROM VALIDATION_TEMPLATE,USER_LIST WHERE VALIDATION_TEMPLATE.ID_USER=USER_LIST.ID_USER $wkgrStrg ORDER BY USER_NAME");
		$sthGetUsers->execute;

		my $hasTemplates=0;
		while (my ($user_ID,$user_name) = $sthGetUsers->fetchrow_array) {
			print "<OPTION value=\"$user_ID\"";
			if ($user_ID eq $dispUserID) {
				print " selected";
			}
			if ($user_ID eq $userID) {
				$hasTemplates=1;
			}
			print ">$user_name</OPTION>";
		}
		unless ($hasTemplates) {
			print "<OPTION value=\"$userID\"";
			print ' selected' if $userID eq $dispUserID;
			print ">$userName</OPTION>";
		}
		print "</SELECT>";
	}
	$dbh->disconnect;
	print qq
|<BR><BR>
<TABLE border=0 cellspacing=0>
|;
	my $lastValType = '';
	foreach my $refRow (@{$refTemplates}) {
		my ($IDTemplate,$templateName,$des,$valType,$searchEngine,$MSType,$isDefault) = @{$refRow};
		if ($des){
			$des = &promsMod::HTMLcompatible($des);
			$des =~ s/<BR>/<BR>&nbsp;&nbsp;&nbsp; /g;
			$des =~ s/&apos;*/\'/g; ### for IE
		}
		else {$des = 'No description';}
		my $defaultStrg = ($isDefault)? '(default)' : '';
		my $valTypeName = ($valType eq 'quali')? 'Qualitative' : ($valType eq 'compa')? 'Comparative' : 'Unknown';
		if ($lastValType ne $valType){
			print "<TR><TD colspan=2><BR><FONT style='font-size:18px;font-weight:bold;'>$valTypeName validation :</FONT></TD></TR>";
		}
		print qq
|<TR bgcolor=$bgColor>
<TD>&nbsp&nbsp</TD>
<TD width=400><FONT style="font-size:18px;font-weight:bold;">$templateName</FONT> <B>$defaultStrg</B><BR>
<B>&nbsp&nbsp&nbsp Validation type: </B>$valTypeName<BR>|;
		if($valType eq 'quali'){
			print qq
|<B>&nbsp&nbsp&nbsp Search engine: </B>$searchEngine<BR>
<B>&nbsp&nbsp&nbsp MS Type : </B>$msTypeFullName{$MSType}<BR>|;
		};
		print qq
|<B>&nbsp&nbsp&nbsp Description : </B>$des
</TD>
<TH width=100>
<INPUT type="button" value="Edit" style="width:75px" onclick="window.location='./manageTemplates.cgi?ACT=edit&ID=$IDTemplate'"><BR>
<INPUT type="button" value="Delete" style="width:75px" onclick="deleteTemplate($IDTemplate,'$templateName')">
</TH>
</TR>
|;
		$bgColor=($bgColor eq $light)? $dark : $light;
		$lastValType = $valType;
	}
	if (scalar @{$refTemplates} == 0) {
		print "<TR><TH colspan=3><BR>(No template saved)</TH></TR>\n";
	}
	print qq
|</TABLE>
</CENTER>
</BODY>
</HTML>
|;
	exit;
}
elsif ($action eq 'edit'){
	### Fetching data ###
	my $id = param('ID');
	my ($name,$des,$valType,$searchEngine,$MSType,$paramStrg,$isDefault) = $dbh->selectrow_array("SELECT NAME,DES,VAL_TYPE,SEARCH_ENGINE,MS_TYPE,PARAM_STRG,IS_DEFAULT FROM VALIDATION_TEMPLATE WHERE ID_VAL_TEMPLATE=$id");
	unless($des){ $des = '' };
	my $sthGetAllNames = $dbh->prepare("SELECT NAME FROM VALIDATION_TEMPLATE WHERE ID_USER='$userID' AND NAME!='$name'");
	$sthGetAllNames->execute;
	my @nameList;
	while (my $refName = $sthGetAllNames->fetchrow_array){
			push @nameList, "\"".$refName."\"";
	}
	my $searchEngineStrg = ($searchEngine)? "AND SEARCH_ENGINE=\'$searchEngine\'" : '';
	unless($searchEngine){ $searchEngine = ''};
	my $defaultTemplate = $dbh->selectrow_array("SELECT NAME FROM VALIDATION_TEMPLATE WHERE ID_USER='$userID' AND VAL_TYPE='$valType' $searchEngineStrg AND IS_DEFAULT=1");
	$dbh->disconnect;

	my ($lightColor,$darkColor)=&promsConfig::getRowColors;
	my $maxRank=&promsConfig::getMaxRank;

	if($valType eq 'quali'){
		### Parsing parameters ###
		my ($selGoodInt) = ($paramStrg =~ /selGoodInt:check:(\d*);/);  $selGoodInt = ($selGoodInt)? 'checked' : '';
		my ($rejBadInt) = ($paramStrg =~ /rejBadInt:check:(\d*);/); $rejBadInt = ($rejBadInt)? 'checked' : '';
		my ($minInt2Sc1,$minInt2Sc2,$minInt2Sc3,$minInt3Sc1,$minInt3Sc2,$minInt3Sc3,$minInt4Sc1,$minInt4Sc2,$minInt4Sc3);
		my ($minIntSc1,$minIntSc2,$minIntSc3);
		if($searchEngine =~ /SEQUEST/ ){
			($minInt2Sc1) = ($paramStrg =~ /minInt2Sc1:text:([\w\.]*);/);
			($minInt2Sc2) = ($paramStrg =~ /minInt2Sc2:text:([\w\.]*);/);
			($minInt2Sc3) = ($paramStrg =~ /minInt2Sc3:text:([\w\.]*);/);
			($minInt3Sc1) = ($paramStrg =~ /minInt3Sc1:text:([\w\.]*);/);
			($minInt3Sc2) = ($paramStrg =~ /minInt3Sc2:text:([\w\.]*);/);
			($minInt3Sc3) = ($paramStrg =~ /minInt3Sc3:text:([\w\.]*);/);
			($minInt4Sc1) = ($paramStrg =~ /minInt4Sc1:text:([\w\.]*);/);
			($minInt4Sc2) = ($paramStrg =~ /minInt4Sc2:text:([\w\.]*);/);
			($minInt4Sc3) = ($paramStrg =~ /minInt4Sc3:text:([\w\.]*);/);
		} else {
			($minIntSc1) = ($paramStrg =~ /minIntSc1:text:([\w\.]*);/);
			($minIntSc2) = ($paramStrg =~ /minIntSc2:text:([\w\.]*);/);
			($minIntSc3) = ($paramStrg =~ /minIntSc3:text:([\w\.]*);/);
		}
		my ($noReject) = ($paramStrg =~ /noReject:check:(\d*);/); $noReject = ($noReject)? 'checked' : '';
		my ($minSize) = ($paramStrg =~ /minSize:text:(\w*);/);
		my ($maxSize) = ($paramStrg =~ /maxSize:text:(\w*);/);
		my ($minDelta) = ($paramStrg =~ /minDelta:text:(\w*);/);
		my ($newMaxRankMS1) = ($paramStrg =~ /newMaxRankMS1:select:(\d*);/); $newMaxRankMS1 = $maxRank unless($newMaxRankMS1);
                my ($newMaxRankMS2) = ($paramStrg =~ /newMaxRankMS2:select:(\d*);/); $newMaxRankMS2 = $maxRank unless($newMaxRankMS2);
		my ($flaggedUp) = ($paramStrg =~ /flaggedUp:check:(\d*);/); $flaggedUp = ($flaggedUp)? 'checked' : '';
		my ($flaggedDown) = ($paramStrg =~ /flaggedDown:check:(\d*);/); $flaggedDown = ($flaggedDown)? 'checked' : '';
		my ($oneInt) = ($paramStrg =~ /oneInt:check:(\d*);/); $oneInt = ($oneInt)? 'checked' : '';
		my ($overPep) = ($paramStrg =~ /overPep:check:(\d*);/); $overPep = ($overPep)? 'checked' : '';
		my ($selProt) = ($paramStrg =~ /selProt:check:(\d*);/); $selProt = ($selProt)? 'checked' : '';
		my ($minPep) = ($paramStrg =~ /minPep:text:(\w*);/);
		my ($minProtSc) = ($paramStrg =~ /minProtSc:text:([\w\.]*);/);
		my ($minCov) = ($paramStrg =~ /minCov:text:(\w*);/);
		my ($bestMatch) = ($paramStrg =~ /bestMatch:check:(\d*);/); $bestMatch = ($bestMatch)? 'checked' : '';
		my ($overProt) = ($paramStrg =~ /overProt:check:(\d*);/); $overProt = ($overProt)? 'checked' : '';

		my $disableMIS;
		if ($MSType eq 'MIS') {
			$disableMIS='';
		}
		else {
			$disableMIS='disabled';
		}
		my $valTypeName = 'qualitative';
		my $selectedDefault = ($isDefault)? 'checked' : '';
		my $defaultTemplateStrg = ($defaultTemplate)? "(current is $defaultTemplate)" : '';

		### HTML ###
		print header(-'content-encoding'=>'no',-charset=>'utf-8');
		warningsToBrowser(1);
		print qq
|<HTML>
<HEAD>
<TITLE>Editing $name Template</TITLE>
<LINK rel=stylesheet href="$promsPath{html}/promsStyle.css" type=text/css>
<STYLE type="text/css">
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
|;
		print "var allNames = [".join(',',@nameList)."]\n";
		print qq
|function checkForm(form){
		if(!checkNameList(form.name.value)){
				return false;
		}
		if (!form.name.value){
				alert('Enter a template name.');
				return false;
		}
		return true;
}
function checkNameList(name){
		for(var i=0;i<allNames.length;i++){
				if(name.toLowerCase() == allNames[i].toLowerCase()){
						alert(name+' is already used.');
						return false;
				}
		}
		return true;
}
|;
		if ($searchEngine=~/SEQUEST/){
			print qq
|function updateScores() {
	var myForm=document.templateForm;
	var scProt=parseFloat(myForm.minProtSc.value);
|;
			for(my $i=2 ; $i<=4; $i++){
				for(my $j=1 ; $j<=3 ; $j++){
					print "\tvar sc$i$j=parseFloat(myForm.minInt"."$i"."Sc$j.value);\n";
				}
			}
			#for(my $i=2 ; $i<=4; $i++){
			#	for(my $j=1 ; $j<=3 ; $j++){
			#		print "        if (!sc$i$j \|\| sc$i"."1 < $absMinScore) {sc$i$j=$absMinScore;}\n";
			#	}
			#}
			for(my $i=2 ; $i<=4; $i++){
				for(my $j=2 ; $j<=3 ; $j++){
					my $previous=$j-1;
					print "\tif (sc$i$j > sc$i$previous) {sc$i$j=sc$i$previous;}\n";
				}
			}
			for(my $i=2 ; $i<=4; $i++){
				print "\tif (scProt < sc$i"."3) {scProt=sc$i"."3;}\n";
			}
			for(my $i=2 ; $i<=4; $i++){
				for(my $j=1 ; $j<=3 ; $j++){
					print "\tmyForm.minInt"."$i"."Sc$j.value=sc$i$j;\n";
				}
			}
			print qq
|	myForm.minProtSc.value=scProt;
}
function activatePepSel() {
	var myForm=document.templateForm;
	var disab;
	if (myForm.selGoodInt.checked==false && myForm.rejBadInt.checked==false){disab=true;}
	else {disab=false;}
|;
			for(my $i=2 ; $i<=4; $i++){
				for(my $j=1 ; $j<=3 ; $j++){
						print "\tmyForm.minInt"."$i"."Sc$j.disabled=disab;\n";
				}
			}
			print qq
|	myForm.noReject.disabled=disab;
	myForm.minSize.disabled=disab;
	myForm.maxSize.disabled=disab;
	myForm.minDelta.disabled=disab;
	myForm.newMaxRank.disabled=disab;
	myForm.flaggedUp.disabled=disab;
	myForm.flaggedDown.disabled=disab;
	myForm.oneInt.disabled=disab;
	myForm.overPep.disabled=disab;
	if (myForm.selGoodInt.checked==false) {
		myForm.oneInt.disabled=true;
	}
}
|;
		} else {
			print qq
|function updateScores() {
	var myForm=document.templateForm;
	var sc1=parseFloat(myForm.minIntSc1.value);
	var sc2=parseFloat(myForm.minIntSc2.value);
	var sc3=parseFloat(myForm.minIntSc3.value);
	var scProt=parseFloat(myForm.minProtSc.value);

	if (sc2 > sc1) {sc2=sc1;}
	if (sc3 > sc2) {sc3=sc2;}
	if (scProt < sc3) {scProt=sc3;}
	myForm.minIntSc1.value=sc1;
	myForm.minIntSc2.value=sc2;
	myForm.minIntSc3.value=sc3;
	myForm.minProtSc.value=scProt;
}
function activatePepSel() {
	var myForm=document.templateForm;
	var disab;
	if (myForm.selGoodInt.checked==false && myForm.rejBadInt.checked==false){disab=true;}
	else {disab=false;}
	myForm.minIntSc1.disabled=disab;
	myForm.minIntSc2.disabled=disab;
	myForm.minIntSc3.disabled=disab;
	myForm.noReject.disabled=disab;
	myForm.minSize.disabled=disab;
	myForm.maxSize.disabled=disab;
	myForm.minDelta.disabled=disab;
	myForm.newMaxRank.disabled=disab;
	myForm.flaggedUp.disabled=disab;
	myForm.flaggedDown.disabled=disab;
	myForm.oneInt.disabled=disab;
	myForm.overPep.disabled=disab;
	if (myForm.selGoodInt.checked==false) {
		myForm.oneInt.disabled=true;
	}
}
|;
		}
		print qq
|function activateProtSel() {
	var myForm=document.templateForm;
	var disab;
	if (myForm.selProt.checked==false){disab=true;} else {disab=false;}
	myForm.minPep.disabled=disab;
	myForm.minProtSc.disabled=disab;
	myForm.minCov.disabled=disab;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Editing Template <FONT color=red>$name</FONT></FONT><BR><BR>

<TABLE align=center>
<FORM name='templateForm' method="post" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="ID" value="$id">
<INPUT type="hidden" name="searchEngine" value="$searchEngine">
<INPUT type="hidden" name="valType" value="$valType">
<TR><TD valign=top><TABLE cellspacing=0>
<TR bgcolor=$darkColor><TH colspan=3><FONT class="title2">&nbsp&nbsp;Name and Description</FONT></TH></TR>
<TR bgcolor=$lightColor><TD width=30></TD><TH align="right" width=100>Name :&nbsp</TH>
<TH align=left><INPUT type="text" name="name" value=$name></TH></TR>
<TR bgcolor=$lightColor><TD width=30></TD><TH align="right" valign="top" width=100>Description :&nbsp</TH>
<TH align=left><TEXTAREA name="description" cols=50 rows=2 value=$des>$des</TEXTAREA></TH></TR>
<TR bgcolor=$lightColor><TD width=30></TD><TH colspan=2 align="left"><INPUT type="checkbox" name="setDefault" value=1 $selectedDefault> Set as default template for $searchEngine $valTypeName validation $defaultTemplateStrg</TH>
<TR bgcolor=$darkColor><TH colspan=3><FONT class="title2">&nbsp&nbsp;Interpretation Selection Rules</FONT></TH></TR>
<TR bgcolor=$lightColor><TD width=30></TD><TH nowrap align=left colspan=2>
<INPUT type="checkbox" name="selGoodInt" value=1 onclick="activatePepSel()" $selGoodInt > Select interpretations meeting the following criteria.<BR>
<INPUT type="checkbox" name="rejBadInt" value=1 onclick="activatePepSel()" $rejBadInt > Reject interpretations that <U>do not</U> meet these criteria.<BR>
<FONT style="font-size:4px"><BR></FONT>+ Criteria:<BR>
|;
		if ($searchEngine=~/SEQUEST/) {
			print qq
|<TABLE cellspacing=0 cellpadding=0>
		<TR><TH align=left>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- minimum score for:&nbsp;&nbsp;</TH>
				<TH align=left nowrap>2<SUP>+</SUP> peptides
1:<INPUT type="text" name="minInt2Sc1" id="minInt2Sc1" value="$minInt2Sc1" size="3" onchange="updateScores()" $disableMIS>,
2:<INPUT type="text" name="minInt2Sc2" id="minInt2Sc2" value="$minInt2Sc2" size="3" onchange="updateScores()" $disableMIS>,
<FONT class="title3">&#8805</FONT>3:<INPUT type="text" name="minInt2Sc3" id="minInt2Sc3" value="$minInt2Sc3" size="3" onchange="updateScores()" $disableMIS> peptides/protein.<TH>
		</TR>
		<TR><TH></TH><TH align=left nowrap>3<SUP>+</SUP> peptides
1:<INPUT type="text" name="minInt3Sc1" id="minInt3Sc1" value="$minInt3Sc1" size="3" onchange="updateScores()" $disableMIS>,
2:<INPUT type="text" name="minInt3Sc2" id="minInt3Sc2" value="$minInt3Sc2" size="3" onchange="updateScores()" $disableMIS>,
<FONT class="title3">&#8805</FONT>3:<INPUT type="text" name="minInt3Sc3" id="minInt3Sc3" value="$minInt3Sc3" size="3" onchange="updateScores()" $disableMIS> peptides/protein.<TH>
		</TR>
		<TR><TH></TH><TH align=left nowrap>4<SUP>+</SUP> peptides
1:<INPUT type="text" name="minInt4Sc1" id="minInt4Sc1" value="$minInt4Sc1" size="3" onchange="updateScores()" $disableMIS>,
2:<INPUT type="text" name="minInt4Sc2" id="minInt4Sc2" value="$minInt4Sc2" size="3" onchange="updateScores()" $disableMIS>,
<FONT class="title3">&#8805</FONT>3:<INPUT type="text" name="minInt4Sc3" id="minInt4Sc3" value="$minInt4Sc3" size="3" onchange="updateScores()" $disableMIS> peptides/protein.<TH>
		<TR></TABLE>
|;
		}
		else{
			print qq
|&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- minimum score if
1:<INPUT type="text" name="minIntSc1" value="$minIntSc1" size="3" onchange="updateScores()" $disableMIS>,
2:<INPUT type="text" name="minIntSc2" value="$minIntSc2" size="3" onchange="updateScores()" $disableMIS>,
<FONT class="title3">&#8805</FONT>3:<INPUT type="text" name="minIntSc3" value="$minIntSc3" size="3" onchange="updateScores()" $disableMIS> peptides/protein.<BR>
|;
		}
		my $optionString = '';
                if($MSType eq 'PMF' || $MSType eq 'SQ'){
                        $optionString .= '&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- and Rank ';
                        $optionString .= 'for PMF ' if($MSType eq 'SQ');
                        $optionString .= '<FONT class="title3">&#8804</FONT><SELECT name="newMaxRankMS1">';
                        foreach my $rank (1..10) {
                                $optionString.="<OPTION value=\"$rank\"";
                                $optionString.=" selected" if $rank==$newMaxRankMS1;
                                $optionString.=">$rank</OPTION>";
                        }
                        $optionString .= "</SELECT><BR>\n";
                }
                if($MSType eq 'MIS' || $MSType eq 'SQ'){
                        $optionString .= '&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- and Rank ';
                        $optionString .= 'for MS/MS ' if($MSType eq 'SQ');
                        $optionString .= '<FONT class="title3">&#8804</FONT><SELECT name="newMaxRankMS2">';
                        foreach my $rank (1..10) {
                                $optionString.="<OPTION value=\"$rank\"";
                                $optionString.=" selected" if $rank==$newMaxRankMS2;
                                $optionString.=">$rank</OPTION>";
                        }
                        $optionString .= "</SELECT><BR>\n";
                }
		print qq
|&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;
<INPUT type="checkbox" name="noReject" value=1 $noReject /> Exclude already rejected interpretations (if any) from count.</FONT><BR>
&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- and size from <INPUT type="text" name="minSize" value="$minSize" size="1">aa to <INPUT type="text" name="maxSize" value="$maxSize" size="1">aa.<BR>
&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- and mass error <FONT class="title3">&#8804</FONT><INPUT type="text" name="minDelta" value="$minDelta" size="4"> Da.<BR>
$optionString
&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;-<INPUT type="checkbox" name="flaggedUp" value=1  $flaggedUp > Interpretation was <IMG src="$promsPath{images}/lightYellow1.gif" hspace=0 border=0 height=11 width=11> flagged.<BR>
&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;-<INPUT type="checkbox" name="flaggedDown" value=1  $flaggedDown > Interpretation was not <IMG src="$promsPath{images}/lightOrange1.gif" hspace=0 border=0 height=11 width=11> flagged.<BR>


<FONT style="font-size:4px"><BR></FONT>
<INPUT type="checkbox" name="oneInt" value=1 $oneInt /> Select only 1 (best) interpretation/query.<BR>
<INPUT type="checkbox" name="overPep" value=1 $overPep /> Overwrite previous selections/rejections.<BR>
<FONT style="font-size:4px"><BR></FONT>
</TH></TR>
<TR bgcolor=$lightColor><TH colspan=3>
<FONT style="font-style:italic;">All proteins matched by selected peptides will be selected unless<BR>Protein Selection Rules are applied.</FONT>
</TH></TR>

<TR bgcolor=$darkColor><TH colspan=3><FONT class="title2">&nbsp&nbsp;Protein Selection Rules</FONT></TH></TR>
<TR bgcolor=$lightColor><TD width=30></TD><TH width=670 align=left colspan=2>
<INPUT type="checkbox" name="selProt" value=1 onclick="activateProtSel()" $selProt /> Select only proteins meeting the following criteria:<BR>
&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- containing at least <INPUT type="text" name="minPep" value="$minPep" size="1"> peptide(s).<BR>
&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- and with score <FONT class="title3">&#8805</FONT><INPUT type="text" name="minProtSc" value="$minProtSc" size="3" $disableMIS><BR>
&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;- and with peptide coverage <FONT class="title3">&#8805</FONT><INPUT type="text" name="minCov" value="$minCov" size="2">%.<BR>
<INPUT type="checkbox" name="bestMatch" value=1 $bestMatch /> Keep only best matching protein(s) for each match group.<BR>
<INPUT type="checkbox" name="overProt" value=1 $overProt /> Overwrite previous exclusions.
<FONT style="font-size=7px;"><BR></FONT>
</TH></TR>
<TR bgcolor=$darkColor><TH colspan=3><INPUT type="submit" name="submit" value="Save"/>&nbsp &nbsp &nbsp<INPUT type="reset" value="Cancel Changes"">&nbsp &nbsp &nbsp<INPUT type="button" value="Cancel" onclick="window.location='./manageTemplates.cgi?ACT=list'">
</TH>
</TR>
</TABLE>
</FORM>
</CENTER>
</BODY>
</HTML>
|;
		exit;
	}
	elsif ($valType eq 'compa'){
		my $analysisString = ' selected analyses';
		### Parsing parameters ###
		my ($selMaxRank) = ($paramStrg =~ /selMaxRank:select:(\d*);/);
		my ($actionPep) = ($paramStrg =~ /actionPep:radio:(\w*);/); my ($actionPep_selRej,$actionPep_flag) = ($actionPep eq 'selRej')? ('checked',''):($actionPep eq 'flag')?('','checked'):('','');
		my ($anaAction) = ($paramStrg =~ /anaAction:radio:(\w*);/); my ($anaAction_inter,$anaAction_diff) = ($anaAction eq 'inter')? ('checked',''):($anaAction eq 'diff')?('','checked'):('','');
		my ($selectPep) = ($paramStrg =~ /selectPep:check:(\d*);/); $selectPep = ($selectPep)? 'checked' : '';
		my ($rejectPep) = ($paramStrg =~ /rejectPep:check:(\d*);/); $rejectPep = ($rejectPep)? 'checked' : '';
		my ($FltAction) = ($paramStrg =~ /FltAction:radio:(\w*);/); my ($FltAction_inter,$FltAction_diff) = ($FltAction eq 'inter')? ('checked',''):($FltAction eq 'diff')?('','checked'):('','');
		my ($FlagUp) = ($paramStrg =~ /FlagUp:check:(\d*);/); $FlagUp = ($FlagUp)? 'checked' : '';
		my ($FlagDown) = ($paramStrg =~ /FlagDown:check:(\d*);/); $FlagDown = ($FlagDown)? 'checked' : '';
		my ($validPept) = ($paramStrg =~ /validPept:check:(\d*);/); $validPept = ($validPept)? 'checked' : '';
		my ($rejMinValPept) = ($paramStrg =~ /rejMinValPept:check:(\d*);/); $rejMinValPept = ($rejMinValPept)? 'checked' : '';
		my ($smPepSc) = ($paramStrg =~ /smPepSc:text:([\d\.]*);/);
		my ($discardEltTime) = ($paramStrg =~ /discardEltTime:check:(\d*);/); $discardEltTime = ($discardEltTime)? 'checked' : '';
		my ($maxEltTime) = ($paramStrg =~ /maxEltTime:text:(\d*);/);
		my ($ignRefMinScore) = ($paramStrg =~ /ignRefMinScore:check:(\d*);/); $ignRefMinScore = ($ignRefMinScore)? 'checked' : '';
		my ($refMinScore) = ($paramStrg =~ /refMinScore:text:([\w\.]*);/);

		my $valTypeName = 'comparative';
		my $selectedDefault = ($isDefault)? 'checked' : '';
		my $defaultTemplateStrg = ($defaultTemplate)? "(current is $defaultTemplate)" : '';


		### HTML ###
		print header(-'content-encoding'=>'no',-charset=>'utf-8');
		warningsToBrowser(1);
		print qq
	|<HTML><HEAD>
	<TITLE>Editing $name Template</TITLE>
	<LINK rel=stylesheet href="$promsPath{html}/promsStyle.css" type=text/css>
	<STYLE type="text/css">
	</STYLE>
	<SCRIPT LANGUAGE="JavaScript">
	|;
		print "var allNames = [".join(',',@nameList)."]\n";
		print qq
	|function checkForm(form){
		if(!checkNameList(form.name.value)){
			return false;
		}
		if (!form.name.value){
			alert('Enter a template name.');
			return false;
		}
		return true;
	}
	function checkNameList(name){
		for(var i=0;i<allNames.length;i++){
			if(name.toLowerCase() == allNames[i].toLowerCase()){
				alert(name+' is already used.');
				return false;
			}
		}
		return true;
	}
	</SCRIPT>
	</HEAD>
	<BODY background="$promsPath{images}/bgProMS.gif">
	<CENTER>
	<FONT class="title">Editing Template <FONT color=red>$name</FONT></FONT><BR><BR>
	<FORM name="filterForm" method="post" onsubmit="return(checkForm(this));">
		<INPUT type="hidden" name="ID" value="$id">
		<INPUT type="hidden" name="valType" value="$valType">
	<TABLE border=0 cellpadding=2>

	<TR>
	<TD nowrap valign=top bgcolor=$darkColor>
	<TABLE bgcolor=$darkColor border=0 cellspacing=0>
		<TR>
			<TH nowrap colspan=2><FONT class="title2">Name and Description</FONT></TH>
		</TR>
		<TR bgcolor=$lightColor>
			<TD align='right'>Name : </TD>
			<TD align='left'>
				<INPUT type='text' name='name' value="$name">
			</TD>
		</TR>
		<TR bgcolor=$lightColor>
			<TD align='right' valign='top'>&nbsp Description : </TD>
			<TD align='left'>
				<TEXTAREA name="description" cols=40 rows=3 value="$des">$des</TEXTAREA>
			</TD>
		</TR>
		<TR bgcolor=$lightColor>
			<TD align='left' colspan=2><INPUT type="checkbox" name="setDefault" value=1 $selectedDefault>Set as default template for $searchEngine $valTypeName validation $defaultTemplateStrg</TD>
		</TR>
		<TR>
			<TH nowrap colspan=2><FONT class="title2"><INPUT type="radio" name="actionPep" value="selRej" $actionPep_selRej>Select/Reject Peptides</FONT></TH>
		</TR>
		<TR bgcolor=$lightColor>
			<TD colspan=2><FONT style="line-height:5px;"><BR></FONT>
				<INPUT type="checkbox" name="selectPep" value=1 $selectPep>Select peptides matching the option below:<BR>
				<INPUT type="checkbox" name="rejectPep" value=1 $rejectPep>Reject peptides that <u>do not </u> match this option:<BR>
				<FONT style="line-height:5px;"><BR></FONT>
				&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;<INPUT type="radio" name="anaAction" value="inter" $anaAction_inter>&nbsp;common to Reference and$analysisString.<BR>
				&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;<INPUT type="radio" name="anaAction" value="diff" $anaAction_diff>&nbsp;unique to$analysisString.<BR>
				&nbsp;(Already selected/rejected peptides will not be affected)<BR>
				<FONT style="line-height:5px;"><BR></FONT>
			</TD>
		</TR>

		<TR>
			<TH nowrap colspan=2>
				<FONT class="title2">or&nbsp<INPUT type="radio" name="actionPep" value="flag" $actionPep_flag>Flag Peptides</FONT>
			</TH>
		</TR>
		<TR bgcolor=$lightColor>
			<TD colspan=2><FONT style="line-height:5px;"><BR></FONT>
				<INPUT type="checkbox" name="FlagUp" value=1 $FlagUp>&nbsp;Flag (<IMG src="$promsPath{images}/lightYellow1.gif" hspace=0 border=0 height=11 width=11>) peptides matching the option below:<BR>
				<INPUT type="checkbox" name="FlagDown" value=1 $FlagDown>&nbsp;Flag (<IMG src="$promsPath{images}/lightOrange1.gif" hspace=0 border=0 height=11 width=11>) peptides that <u>do not </u> match this option:<BR>
				<FONT style="line-height:5px;"><BR></FONT>
				&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;<INPUT type="radio" name="FltAction" value="inter" $FltAction_inter>&nbsp;common to Reference and$analysisString.<BR>
				&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp;<INPUT type="radio" name="FltAction" value="diff" $FltAction_diff>&nbsp;unique to$analysisString.<BR>
				<FONT style="line-height:5px;"><BR></FONT>
			</TD>
		</TR>

		<TR>
			<TH nowrap colspan=2><FONT class="title2">Options</FONT></TH>
		</TR>
		<TR bgcolor=$lightColor>
			<TD nowrap colspan=2><FONT style="line-height:5px;"><BR></FONT>
				&nbsp;Compare only peptides with rank &le;<SELECT name="selMaxRank" style="font-weight:bold;">
	|;
		foreach my $rank (1..10) {
				print "<OPTION value=\"$rank\"";
				print " selected" if $rank == $selMaxRank;
				print ">$rank</OPTION>";
		}
		print qq
	|</SELECT><BR>
				<INPUT type="checkbox" name="validPept" value=1 $validPept>Use only peptides already selected in Reference.<BR>
				<INPUT type="checkbox" name="rejMinValPept" value=1 $rejMinValPept>Ignore peptides in$analysisString with score below <INPUT type="text" name="smPepSc" value=$smPepSc size="3" >. <BR>
				<INPUT type="checkbox" name="ignRefMinScore" value=1 $ignRefMinScore>Ignore peptides in Reference with score below <INPUT type="text" name="refMinScore" value=$refMinScore size="3" >. <BR>
				<INPUT type="checkbox" name="discardEltTime" value=1 $discardEltTime>Ignore peptides with elution time differing by more than
				<INPUT type="text" name="maxEltTime" value=$maxEltTime size="2" > mn.<BR>
				<FONT style="line-height:5px;"><BR></FONT>
			</TD>
		</TR>
		<TR>
			<TH nowrap colspan=2><INPUT type="submit" name="submit" value=" Proceed "/>
				&nbsp &nbsp &nbsp<INPUT type="reset" value="  Clear  " />
				&nbsp &nbsp &nbsp<INPUT type="button" value="Cancel" onClick="window.location='./manageTemplates.cgi?ACT=list';">
			</TH>
		</TR>
	</TABLE>
	</TABLE>
	</FORM>
	</BODY>
	</HTML>
	|;
	}
	else {
		print header(-'content-encoding'=>'no',-charset=>'utf-8');
		print "<HTML><BODY> Unknown validation type !</BODY></HTML>";
		exit;
	}
}

sub processForm {
	### Fetching all parameters ###
	my $name = param('name');
	my $id = param('ID');
	my $des = param('description');
	my $setDefault = (param('setDefault'))? 1 : 0;
	my $valType = param('valType');
	my $paramStrg = '';
	my $searchEngine = param('searchEngine');

	if($valType eq 'quali'){
			my $selGoodInt = param('selGoodInt');
		my $rejBadInt = param('rejBadInt');
		$paramStrg .= "selGoodInt:check:$selGoodInt;rejBadInt:check:$rejBadInt;";
		my ($minInt2Sc1,$minInt2Sc2,$minInt2Sc3,$minInt3Sc1,$minInt3Sc2,$minInt3Sc3,$minInt4Sc1,$minInt4Sc2,$minInt4Sc3);
		my ($minIntSc1,$minIntSc2,$minIntSc3);
		if($searchEngine =~ /SEQUEST/ ){
			$minInt2Sc1 = param('minInt2Sc1');
			$minInt2Sc2 = param('minInt2Sc2');
			$minInt2Sc3 = param('minInt2Sc3');
			$minInt3Sc1 = param('minInt3Sc1');
			$minInt3Sc2 = param('minInt3Sc2');
			$minInt3Sc3 = param('minInt3Sc3');
			$minInt4Sc1 = param('minInt4Sc1');
			$minInt4Sc2 = param('minInt4Sc2');
			$minInt4Sc3 = param('minInt4Sc3');
			$paramStrg .= "minInt2Sc1:text:$minInt2Sc1;minInt2Sc2:text:$minInt2Sc2;minInt2Sc3:text:$minInt2Sc3;minInt3Sc1:text:$minInt3Sc1;minInt3Sc2:text:$minInt3Sc2;;minInt3Sc3:text:$minInt3Sc3;minInt4Sc1:text:$minInt4Sc1;minInt4Sc2:text:$minInt4Sc2;minInt4Sc3:text:$minInt4Sc3;"
		} else {
			$minIntSc1 =param('minIntSc1');
			$minIntSc2 =param('minIntSc2');
			$minIntSc3 =param('minIntSc3');
			$paramStrg .= "minIntSc1:text:$minIntSc1;minIntSc2:text:$minIntSc2;minIntSc3:text:$minIntSc3;"
		}

		my $noReject = param('noReject');
		my $minSize = param('minSize');
		my $maxSize = param('maxSize');
		my $minDelta = param('minDelta');
		my $newMaxRankMS1 =param('newMaxRankMS1');
                my $newMaxRankMS2 =param('newMaxRankMS2');
		my $flaggedUp = param('flaggedUp');
		my $flaggedDown =param('flaggedDown');
		my $oneInt =param('oneInt');
		my $overPep =param('overPep');
		my $selProt = param('selProt');
		my $minPep =param('minPep');
		my $minProtSc =param('minProtSc');
		my $minCov = param('minCov');
		my $bestMatch = param('bestMatch');
		my $overProt =param('overProt');
		$paramStrg .= "noReject:check:$noReject;minSize:text:$minSize;maxSize:text:$maxSize;minDelta:text:$minDelta;flaggedUp:check:$flaggedUp;flaggedDown:check:$flaggedDown;oneInt:check:$oneInt;overPep:check:$overPep;selProt:check:$selProt;minPep:text:$minPep;minProtSc:text:$minProtSc;minCov:text:$minCov;bestMatch:check:$bestMatch;overProt:check:$overProt;";
                $paramStrg .= "newMaxRankMS1:select:$newMaxRankMS1;" if($newMaxRankMS1);
                $paramStrg .= "newMaxRankMS2:select:$newMaxRankMS2;" if($newMaxRankMS2);
	}
	elsif($valType eq 'compa'){
		my $selMaxRank = param('selMaxRank');
		my $actionPep = param('actionPep');
		my $anaAction = param('anaAction');
		my $selectPep = param('selectPep');
		my $rejectPep = param('rejectPep');
		my $FltAction = param('FltAction');
		my $FlagUp = param('FlagUp');
		my $FlagDown = param('FlagDown');
		my $validPept = param('validPept');
		my $rejMinValPept = param('rejMinValPept');
		my $smPepSc = param('smPepSc');
		my $discardEltTime = param('discardEltTime');
		my $maxEltTime = param('maxEltTime');
		my $ignRefMinScore = param('ignRefMinScore');
		my $refMinScore = param('refMinScore');
		$paramStrg .= "selMaxRank:select:$selMaxRank;actionPep:radio:$actionPep;anaAction:radio:$anaAction;selectPep:check:$selectPep;rejectPep:check:$rejectPep;FltAction:radio:$FltAction;FlagUp:check:$FlagUp;FlagDown:check:$FlagDown;validPept:check:$validPept;rejMinValPept:check:$rejMinValPept;smPepSc:text:$smPepSc;discardEltTime:check:$discardEltTime;maxEltTime:text:$maxEltTime;ignRefMinScore:check:$ignRefMinScore;refMinScore:text:$refMinScore;";
	}
	### Updating DB ###
	my $searchEngineStrg = ($searchEngine)? "AND SEARCH_ENGINE=\'$searchEngine\'" : '';

	if($setDefault){
		$dbh->do("UPDATE VALIDATION_TEMPLATE SET IS_DEFAULT=0 WHERE ID_USER='$userID' AND VAL_TYPE='$valType' $searchEngineStrg AND IS_DEFAULT=1") || die $dbh->errstr;
	}
	$des = $dbh->quote($des);
	$name = $dbh->quote($name);
	$dbh->do("UPDATE VALIDATION_TEMPLATE SET NAME=$name, DES=$des, PARAM_STRG='$paramStrg', IS_DEFAULT=$setDefault WHERE ID_VAL_TEMPLATE=$id") || die $dbh->errstr ;
	$dbh->commit;
	$dbh->disconnect;
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	print qq
|<HTML>
<BODY onload="window.location='./manageTemplates.cgi?ACT=list';"></BODY>
</HTML>
|;
	exit;
}

sub deleteTemplate{
	my $id = param('ID');
	$dbh->do("DELETE FROM VALIDATION_TEMPLATE WHERE ID_VAL_TEMPLATE=$id") || die $dbh->errstr;
	$dbh->commit;
	$dbh->disconnect;
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	print qq
|<HTML>
<BODY onload="window.location='./manageTemplates.cgi?ACT=list'"></BODY>
</HTML>
|;
	exit;
}

####>Revision history<####
# 1.0.7 Color on main "Go to" options (PP 25/03/16)
# 1.0.6 GPL license (PP 23/09/13)
# 1.0.5 Called by promsMain.cgi (PP 13/12/12)
# 1.0.4 Fixes bug when logged user has no template: 1st user in list was displayed (PP 20/02/12)
# 1.0.3 Manager access management (PP 01/10/11)
# 1.0.2 SQ ranks management (FY 1/06/11)
# 1.0.1 Code reformatting & minor changes in user listing (PP 11/04/11)
# 1.0.0 New script to manage user validation templates (FY 03/03/2011)
