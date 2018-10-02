#!/usr/local/bin/perl -w

#############################################################################
# view2dGel.cgi                  1.1.6                                      #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                   #
# Contact: myproms@curie.fr                                                 #
# Display 2D gel with a Java applet                                         #
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
$| = 1;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use strict;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};

#print header; warningsToBrowser(1); # DEBUG
##################
####>Argument<####
##################
my $action=param('ACT');
my $gelID=param('id_gel');


############################
####>Auto-reload Window<#### required to get window true URL (VPN & https!) & pass it to Applet
############################
#if ($action eq 'open') {
#	my $projectID=param('id_project');
#
#	print header; warningsToBrowser(1);
#	print qq
#|<HTML>
#<HEAD>
#<SCRIPT LANGUAGE="JavaScript">
#function getWindowURL() {
#	alert('URL='+window.location.href);
#	var modURLs=window.location.href.split('5');
#	document.gelForm.ref_URL.value=modURLs.join('#');
#	alert('modURL='+document.gelForm.ref_URL.value);
#	document.gelForm.submit();
#}
#</SCRIPT>
#</HEAD>
#<BODY onload="getWindowURL()">
#<FORM name="gelForm" action="view2dGel.cgi" target="2D_Gel" method="POST">
#<INPUT type="hidden" name="id_project" value="$projectID"/>
#<INPUT type="hidden" name="id_gel" value="$gelID"/>
#<INPUT type="hidden" name="ACT" value="main"/>
#<INPUT type="hidden" name="ref_URL" value=""/>
#</FORM></BODY></HTML>
#|;
#	exit;
#}
#if ($action eq 'open2') {
#	my $projectID=param('id_project');
#
#	print header; warningsToBrowser(1);
#	print qq
#|<HTML>
#<HEAD>
#<SCRIPT LANGUAGE="JavaScript">
#function processURL(myForm) {
#	alert('pastedURL='+myForm.pastedURL.value);
#	var modURLs=myForm.pastedURL.value.split('5');
#	document.gelForm.ref_URL.value=modURLs.join('#');
#	alert('modURL='+myForm.ref_URL.value);
#	document.gelForm.submit();
#}
#</SCRIPT>
#</HEAD>
#<BODY>
#<FORM name="gelForm" action="view2dGel.cgi" target="2D_Gel" method="POST" onsubmit="processURL(this)">
#<INPUT type="hidden" name="id_project" value="$projectID"/>
#<INPUT type="hidden" name="id_gel" value="$gelID"/>
#<INPUT type="hidden" name="ACT" value="main"/>
#<B>Window URL:</B><INPUT type="text" name="pastedURL" value="" size=150/>
#<INPUT type="hidden" name="ref_URL" value=""/>
#&nbsp;<INPUT type="submit" value="GO"/>
#</FORM>
#<BR>
#HTTP_X_FORWARDED_FOR: $ENV{HTTP_X_FORWARDED_FOR}<BR>
#HTTP_REFERER: $ENV{HTTP_REFERER}<BR>
#<SCRIPT LANGUAGE="JavaScript">
#alert('HTTP_X_FORWARDED_FOR: $ENV{HTTP_X_FORWARDED_FOR}');
#alert('HTTP_REFERER: $ENV{HTTP_REFERER}');
#</SCRIPT>
#</BODY></HTML>
#|;
#	exit;
#}

#######################
####>Applet Window<####
#######################
if ($action eq 'main') {
	my $projectID=param('id_project');
	#my $refURL=param('ref_URL');
	#my $modURL=$refURL;
	#$refURL=~s/#/5/g;

	####>Connecting to DB<####
	my $dbh=&promsConfig::dbConnect;

	my ($gelName)=$dbh->selectrow_array("SELECT NAME FROM GEL2D WHERE ID_GEL2D=$gelID");
	my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
	my $editable=($userInfo[2]->{$projectID} eq 'guest')? 0 : 1;
	my ($imageFile)=$dbh->selectrow_array("SELECT IMAGE_FILE FROM GEL2D WHERE ID_GEL2D=$gelID");
	my ($cgiDirName)=($promsPath{'cgi'}=~/^\/([^\/]+)/);
	my $gelImageURL;
	if ($ENV{'HTTP_FRONT_END_HTTPS'} && $ENV{'HTTP_FRONT_END_HTTPS'} eq 'on') { # reverse proxy: https protocol
		$gelImageURL.="https://$ENV{HTTP_HOST}";
	}
	elsif ($promsPath{'vpn'} && $ENV{'HTTP_REFERER'} && $ENV{'HTTP_REFERER'}=~/vpn/) {$gelImageURL=$promsPath{'vpn'};} # VPN connection detected
	else { # intranet connection
		if ($ENV{'SCRIPT_URI'}=~/^http/) {($gelImageURL=$ENV{'SCRIPT_URI'})=~s/\/$cgiDirName.+//;}
		else {
			$gelImageURL=($ENV{'HTTPS'} && $ENV{'HTTPS'}=~/on/)? 'https' : 'http';
			$gelImageURL.="://$ENV{HTTP_HOST}"; # $ENV{SERVER_NAME} also works
		}
	}
	$gelImageURL.="$promsPath{gel}/project_$projectID/gel2D_$gelID.jpg";

	####>Starting HTML<####
	print header;
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>2D Gel</TITLE>
<LINK rel="icon" type="image/png" href="$promsPath{images}/myProMS_icon.png">
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
var id_project=$projectID;
var id_gel=$gelID;
//
var navFrame=parent.opener.parent.navFrame;
var itemFrame=parent.opener.parent.itemFrame;
var gelBranchID='gel2d:'+$gelID;
function updateFrame(action,idspot) {
	if (idspot>0) {
		if(action=='select') {
			if (navFrame.getSelectedBranchID() != gelBranchID) { // itemFrame must be reloaded
				navFrame.selItemBranchID='spot:'+idspot;
				navFrame.selectItem(navFrame.getItemIndex(gelBranchID));
				navFrame.showSelectedItem(gelBranchID);
			}
			else { // just select spot
				itemFrame.selectItem(itemFrame.getItemIndex('spot:'+idspot));
				itemFrame.showSelectedItem('spot:'+idspot);
			}
		}
		else {
			if (action=='delete') {
				navFrame.location='$promsPath{cgi}/openProject.cgi?&ACT=nav&ID='+parent.id_project+'&branchID='+gelBranchID;
			}
			else {
				navFrame.selItemBranchID='spot:'+idspot;
				navFrame.location='$promsPath{cgi}/openProject.cgi?&ACT=nav&ID='+parent.id_project+'&branchID='+gelBranchID+'&itemBranchID=spot:'+idspot;
			}
		}
	}
}
//Fetching & sending list of spots to Applet
function getSpotList() {
	for (var i=0; i< spotList.length; i++) {
		// 11 parameters
		document.gel2d.addSpot(spotList[i].id,spotList[i].name,spotList[i].x_pos,spotList[i].y_pos,spotList[i].ie_point,spotList[i].mw,spotList[i].intensity,spotList[i].ext_id,spotList[i].samp_id,spotList[i].samp_name,spotList[i].top_protein);
	}
}
//Free samples
function getFreeSamples(caller) {
	doAjaxQuery('WHAT=freeSamples&id_gel=$gelID',sendFreeSamples,caller);
}
function sendFreeSamples(sampleString,caller) {
	document.gel2d.sendFreeSamples(sampleString,caller);
}
function manageSpot(action,action2){
	var paramString='WHAT='+action2+'&id_gel=$gelID'+action;
	doAjaxQuery(paramString,updateSpot,action2);
}
function updateSpot(spotToAdd,action){
	if(action=='createSpot'){
		document.gel2d.addSpot(spotToAdd,action);
	}
	if(action=='editSpot'){
		document.gel2d.updateSpot(spotToAdd,action);
	}
}
function deleteSpot(ID_SPOT,ID_SAMPLE){
	var paramString='WHAT=deleteSpot&id_gel=$gelID&idspot='+ID_SPOT+'&idsample='+ID_SAMPLE;
	doAjaxQuery(paramString,deleteJavaSpot,'d');
}
function deleteJavaSpot(ID_SPOT,action){
	document.gel2d.deleteSpot(ID_SPOT);
}
function searchProtein(searchString,alias,description){
	var paramString='WHAT=searchProtein&id_gel=$gelID&searchString='+searchString+'&alias='+alias+'&description='+description;
	doAjaxQuery(paramString,printSearch,'d');
}
function printSearch(spotList,action){
	document.gel2d.printSearch(spotList);
}

// AJAX --->
var XHR=null;
function doAjaxQuery(paramStrg,processFunction,callFromJavaButton) {
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
	XHR.open("GET","$promsPath{cgi}/view2dGel.cgi?ACT=ajax&"+paramStrg,true);
    XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			processFunction(XHR.responseText,callFromJavaButton);
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
</SCRIPT>
</HEAD>
<BODY>
<CENTER>
<FONT class="title">2D Gel <FONT color=#DD0000>$gelName</FONT></FONT><BR>
<APPLET name="gel2d" code="ImageView.class" codebase="$promsPath{java_gel2d}" archive="gel2d.jar" width="1004" height="440" mayscript=true>
<PARAM name="gel2D_ID" value="$gelID">
<PARAM name="editable" value="$editable">
<PARAM name="gel2D_file" id="gel2D_file" value="$gelImageURL">
<PARAM name="image_file" value="$imageFile">
</APPLET>
|;
	print "<FONT style=\"font-size:11px;font-weight:bold;font-style:italic;\">Double-click on gel image to add spots</FONT>\n" if $editable;
	print qq
|</CENTER>
<SCRIPT LANGUAGE="JavaScript">
var spotList=new Array();
|;
	my $sthSpot=$dbh->prepare("SELECT ID_SPOT,NAME,X_POS,Y_POS,ISOELECTRIC_POINT,MOLECULAR_WEIGHT,INTENSITY,EXTERNAL_ID FROM SPOT WHERE ID_GEL2D=$gelID ORDER BY LOWER(NAME) ASC");
	my $sthSamp=$dbh->prepare("SELECT ID_SAMPLE,NAME FROM SAMPLE WHERE ID_SPOT=?");

	$sthSpot->execute;
	my $i=0;
	while (my ($spotID,$name,$xPos,$yPos,$iePoint,$mw,$intensity,$extID)=$sthSpot->fetchrow_array) {
		$sthSamp->execute($spotID);
		my ($sampID,$sampName)=$sthSamp->fetchrow_array;
		$intensity=-1 unless defined $intensity;
		$extID=($extID)? "'$extID'" : 'null';
		$sampID=-1 unless $sampID;
		$iePoint=-1 unless $iePoint;
		$mw=-1 unless $mw;
		$sampName=($sampName)? "'$sampName'" : 'null';
		my $majProt=&getMajProt($sampID,$dbh);
		$majProt=($majProt)? "'$majProt'" : 'null';
		print qq
|spotList[$i]=new Object();
	spotList[$i].id=$spotID;
	spotList[$i].name='$name';
	spotList[$i].x_pos=$xPos;
	spotList[$i].y_pos=$yPos;
	spotList[$i].ie_point=$iePoint;
	spotList[$i].mw=$mw;
	spotList[$i].intensity=$intensity;
	spotList[$i].ext_id=$extID;
	spotList[$i].samp_id=$sampID;
	spotList[$i].samp_name=$sampName;
	spotList[$i].top_protein=$majProt;
|;
		$i++;
	}
	$sthSpot->finish;
	$sthSamp->finish;

	$dbh->disconnect;
	print qq
|</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}
#}elsif($action eq 'createSpot') {
#
#}

#######################
####>AJAX Requests<####
#######################
elsif ($action eq 'ajax') {
	my $what=param('WHAT');

	####>Connecting to DB<####
	my $dbh=&promsConfig::dbConnect;

	####>Free Samples<####
	if ($what eq 'freeSamples') {
		my $sthFS=$dbh->prepare("SELECT ID_SAMPLE,NAME FROM SAMPLE WHERE ID_SPOT IS NULL AND ID_EXPERIMENT=(SELECT ID_EXPERIMENT FROM GEL2D WHERE ID_GEL2D=$gelID) ORDER BY NAME ASC");
		$sthFS->execute;
		my @sampleList;
		while (my @thisSample=$sthFS->fetchrow_array) {
			push @sampleList,join('@@',@thisSample);
		}
		$sthFS->finish;
		print header(-type=>'text/plain');
		if( $#sampleList eq -1){
			print '@#NONE#@';
		}else{
			print join('##',@sampleList);
		}
	}
	####>Create Spot<####
	elsif ($what eq 'createSpot') {
		my ($maxIdSpot)=$dbh->selectrow_array("SELECT MAX(ID_SPOT) FROM SPOT");
		my ($spotName,$x_pos,$y_pos,$pi,$pw,$intensity,$externalid,$sampleID,$sampleName)=(param('spotName'),param('x_pos'),param('y_pos'),param('pi'),param('pw'),param('intensity'),param('externalid'),param('sampleid'),param('sampleName'));
		my $sthIns=$dbh->prepare("INSERT INTO SPOT (ID_SPOT,ID_GEL2D,NAME,X_POS,Y_POS,ISOELECTRIC_POINT,MOLECULAR_WEIGHT,INTENSITY,EXTERNAL_ID,START_DATE,UPDATE_USER) VALUES (?,$gelID,?,?,?,?,?,?,?,NOW(),'$userID')");
		$intensity=undef unless $intensity;
		$pi=undef unless $pi;
		$pw=undef unless $pw;
		$sthIns->execute(++$maxIdSpot,param('spotName'),param('x_pos'),param('y_pos'),$pi,$pw,$intensity,param('externalid'));

		if((param('sampleid')<0) && !(param('sampleName') eq 'NONE')){
			#create a new sample
			my ($maxIdSample)=$dbh->selectrow_array("SELECT MAX(ID_SAMPLE) FROM SAMPLE");
			$maxIdSample++;
			my ($idexperiment)=$dbh->selectrow_array("SELECT ID_EXPERIMENT FROM GEL2D WHERE ID_GEL2D=$gelID");
			$dbh->do("INSERT INTO SAMPLE (ID_SAMPLE,ID_EXPERIMENT,ID_SPOT,NAME,START_DATE,UPDATE_USER) VALUES($maxIdSample,$idexperiment,$maxIdSpot,'$sampleName',NOW(),'$userID')");
			$dbh->commit;
			$sampleID=$maxIdSample;
		}elsif($sampleID>0){
			#sample already exists -> it has a majProt
			my ($idexperiment)=$dbh->selectrow_array("SELECT ID_EXPERIMENT FROM GEL2D WHERE ID_GEL2D=$gelID");
			my ($displayPos)=$dbh->selectrow_array("SELECT DISPLAY_POS FROM SAMPLE WHERE ID_SAMPLE=$sampleID");
			$dbh->do("UPDATE SAMPLE SET ID_SPOT=$maxIdSpot , DISPLAY_POS=null, UPDATE_DATE=NOW() , UPDATE_USER='$userID' WHERE ID_SAMPLE=$sampleID");
			# display_pos has to be updated
			my ($maxD)=$dbh->selectrow_array("SELECT max(DISPLAY_POS) FROM SAMPLE WHERE ID_EXPERIMENT=$idexperiment");
			if($maxD){#Is there still no spot associated
				for(my $i = $displayPos+1 ; $i <= $maxD ; $i++){
					my $newd=$i-1;
					$dbh->do("UPDATE SAMPLE SET DISPLAY_POS=$newd , UPDATE_DATE=NOW() , UPDATE_USER='$userID' WHERE ID_EXPERIMENT=$idexperiment AND DISPLAY_POS=$i");
				}
			}
			$dbh->commit;
		}
		my $protMaj=&getMajProt($sampleID,$dbh);
		
		($pi,$pw,$intensity,$externalid,$protMaj) = map { $_ // ''} ($pi,$pw,$intensity,$externalid,$protMaj);
		print header(-type=>'text/plain');
		print "$maxIdSpot##$spotName##$x_pos##$y_pos##$pi##$pw##$intensity##$externalid##$protMaj##$sampleID##$sampleName";

		$dbh->commit;
	}
	####>Delete Spot<####
	elsif($what eq 'deleteSpot'){
		if(param('idsample') > -1){# If the spot is associated to a real sample
			my $sampleID=param('idsample');
			$dbh->do("UPDATE SAMPLE SET ID_SPOT = null , UPDATE_DATE=NOW() , UPDATE_USER='$userID' WHERE ID_SAMPLE=$sampleID");
			$dbh->commit;
			my ($analysis)=$dbh->selectrow_array("SELECT COUNT(*) FROM ANALYSIS WHERE ID_SAMPLE=$sampleID");
			if($analysis eq 0){#No analysis related to the sample -> delete the sample
				$dbh->do("DELETE FROM SAMPLE WHERE ID_SAMPLE=$sampleID");
				$dbh->commit;
			}else{
				# Has to update the display pos if there is still at least one analysis related to this sample
				my ($idexperiment)=$dbh->selectrow_array("SELECT ID_EXPERIMENT FROM GEL2D WHERE ID_GEL2D=$gelID");
				my ($maxDisplayPos)=$dbh->selectrow_array("SELECT max(DISPLAY_POS) FROM SAMPLE WHERE ID_EXPERIMENT=$idexperiment");
				$maxDisplayPos++;
				$dbh->do("UPDATE SAMPLE SET DISPLAY_POS=$maxDisplayPos , UPDATE_DATE=NOW() , UPDATE_USER='$userID' WHERE ID_SAMPLE=$sampleID");
				$dbh->commit;
			}
		}
		my $idspot=param('idspot');
		$dbh->do("DELETE FROM SPOT WHERE ID_SPOT=$idspot");
		$dbh->commit;
		print header(-type=>'text/plain');
		print "$idspot";

	}

	####>Search for proteins<####
	elsif($what eq 'searchProtein') {
		my $query = "SELECT DISTINCT(SPOT.ID_SPOT) FROM SPOT,SAMPLE,ANALYSIS,ANALYSIS_PROTEIN,PROTEIN WHERE ID_GEL2D=$gelID AND SPOT.ID_SPOT=SAMPLE.ID_SPOT AND SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND ANALYSIS.ID_ANALYSIS=ANALYSIS_PROTEIN.ID_ANALYSIS AND ANALYSIS_PROTEIN.ID_PROTEIN=PROTEIN.ID_PROTEIN AND (";
		my $search=param('searchString');
		if(param('alias') == 1){
			$query.="lower(ALIAS)=lower(\"$search\") ";
		}

		if(param('alias') == 1 && param('description') == 1) {
			$query.="OR ";
		}

		if( param('description') == 1 ) {
			$search=lc($search);
			$query.="lower(PROT_DES) like \"%$search%\"";
		}

		$query.=")";

		my $sthSearch=$dbh->prepare($query);
		$sthSearch->execute;

		my $spotList="";
		while (my ($spotname)=$sthSearch->fetchrow_array) {
			$spotList.="$spotname@@";
		}
		$spotList='@@NONE@@' unless $spotList;
		print header(-type=>'text/plain');
		print "$spotList";


	}

	####>Edit spot information<####
	elsif($what eq 'editSpot'){
		#Get all the ids
		my ($idspot,$spotName,$x_pos,$y_pos,$pi,$pw,$intensity,$externalid,$sampleID,$sampleName)=(param('idspot'),param('spotName'),param('x_pos'),param('y_pos'),param('pi'),param('pw'),param('intensity'),param('externalid'),param('sampleid'),param('sampleName'));

		# Check if the association
		#print header(-type=>'text/plain');
		my ($oldSampleID)=$dbh->selectrow_array("SELECT ID_SAMPLE FROM SAMPLE WHERE ID_SPOT=$idspot");
		my ($maxIdSample)=$dbh->selectrow_array("SELECT MAX(ID_SAMPLE) FROM SAMPLE");
		if($maxIdSample){
			$maxIdSample++;
		}
		my ($idexperiment)=$dbh->selectrow_array("SELECT ID_EXPERIMENT FROM GEL2D WHERE ID_GEL2D=$gelID");

		##########################################################################
		####>Before updating the spot info, the sample info has to be updated>####
		##########################################################################
		if(!$oldSampleID){#The spot was not associated before
			#print "\nThe spot wasn't associated before";
			if(param('isnew') == 1){
				#print "\nThe spot is associated to a new sample maxiddsample=$maxIdSample";
				$dbh->do("INSERT INTO SAMPLE (ID_SAMPLE,ID_EXPERIMENT,ID_SPOT,NAME,START_DATE,UPDATE_USER) VALUES($maxIdSample,$idexperiment,$idspot,'$sampleName',NOW(),'$userID')");
				$sampleID=$maxIdSample;
				$dbh->commit;
			}else{
				if(!($sampleName eq 'NONE')){
					#print "\nThe spot is associated to an old sample";
					my ($displaypos)=$dbh->selectrow_array("SELECT DISPLAY_POS FROM SAMPLE WHERE ID_SAMPLE=$sampleID");
					$dbh->do("UPDATE SAMPLE SET ID_SPOT=$idspot , DISPLAY_POS=null, UPDATE_DATE=NOW() , UPDATE_USER='$userID' WHERE ID_SAMPLE=$sampleID");
					$dbh->commit;
					my ($maxDisplayPos)=$dbh->selectrow_array("SELECT max(DISPLAY_POS) FROM SAMPLE WHERE ID_EXPERIMENT=$idexperiment");
					#print "\nThe spot is associated to an old sample";
					if($maxDisplayPos){# Take care of the display_pos of samples not already associated
						$maxDisplayPos++;
						for(my $i = $displaypos+1 ; $i < $maxDisplayPos ; $i++){
							my $newd=$i-1;
							$dbh->do("UPDATE SAMPLE SET DISPLAY_POS=$newd , UPDATE_DATE=NOW() , UPDATE_USER='$userID' WHERE ID_EXPERIMENT=$idexperiment AND DISPLAY_POS=$i");
						}
					}
					$dbh->commit;
				}
			}
		}else{#The spot was associated before
			#print "\nThe spot was associated before";
			if($oldSampleID != $sampleID){
				if(!($sampleName eq 'NONE')){
					#print "\nBut it is associated to another sample";
					if(param('isnew') == 1){#Create a new sample
						#print "\nA new sample";
						$dbh->do("INSERT INTO SAMPLE (ID_SAMPLE,ID_EXPERIMENT,ID_SPOT,NAME,START_DATE,UPDATE_USER) VALUES($maxIdSample,$idexperiment,$idspot,'$sampleName',NOW(),'$userID')");
						$sampleID=$maxIdSample;
						$dbh->commit;
					}else{
						#print "\nAn old sample";
						my ($displaypos)=$dbh->selectrow_array("SELECT DISPLAY_POS FROM SAMPLE WHERE ID_SAMPLE=$sampleID");
						$dbh->do("UPDATE SAMPLE SET ID_SPOT=$idspot , DISPLAY_POS=null, UPDATE_DATE=NOW() , UPDATE_USER='$userID' WHERE ID_SAMPLE=$sampleID");
						my ($maxDisplayPos)=$dbh->selectrow_array("SELECT max(DISPLAY_POS) FROM SAMPLE WHERE ID_EXPERIMENT=$idexperiment");
						$maxDisplayPos++;
						if($maxDisplayPos){# Take care of the display_pos of samples not already associated
							for(my $i = $displaypos+1 ; $i < $maxDisplayPos ; $i++){
								my $newd=$i-1;
								$dbh->do("UPDATE SAMPLE SET DISPLAY_POS=$newd , UPDATE_DATE=NOW() , UPDATE_USER='$userID' WHERE ID_EXPERIMENT=$idexperiment AND DISPLAY_POS=$i");
								$dbh->commit;
							}
						}
						$dbh->commit;
					}
					my ($analysis)=$dbh->selectrow_array("SELECT COUNT(*) FROM ANALYSIS WHERE ID_SAMPLE=$oldSampleID");
					if($analysis eq 0){#No analysis related to the sample -> delete the sample
						$dbh->do("DELETE FROM SAMPLE WHERE ID_SAMPLE=$oldSampleID");
						$dbh->commit;
					}else{
						my ($maxDisplayPos)=$dbh->selectrow_array("SELECT max(DISPLAY_POS) FROM SAMPLE WHERE ID_EXPERIMENT=$idexperiment");
						$maxDisplayPos++;
						$dbh->do("UPDATE SAMPLE SET ID_SPOT=null , DISPLAY_POS=$maxDisplayPos, UPDATE_DATE=NOW() , UPDATE_USER='$userID' WHERE ID_SAMPLE=$oldSampleID");
						$dbh->commit;
					}
				}else{#it is associated to NONE, so, the only thing to do is to handle the display pos
					my ($analysis)=$dbh->selectrow_array("SELECT COUNT(*) FROM ANALYSIS WHERE ID_SAMPLE=$oldSampleID");
					if($analysis eq 0){#No analysis related to the sample -> delete the sample
						$dbh->do("DELETE FROM SAMPLE WHERE ID_SAMPLE=$oldSampleID");
						$dbh->commit;
					}else{
						my ($maxDisplayPos)=$dbh->selectrow_array("SELECT max(DISPLAY_POS) FROM SAMPLE WHERE ID_EXPERIMENT=$idexperiment");
						$maxDisplayPos++;
						$dbh->do("UPDATE SAMPLE SET ID_SPOT=null , DISPLAY_POS=$maxDisplayPos, UPDATE_DATE=NOW() , UPDATE_USER='$userID' WHERE ID_SAMPLE=$oldSampleID");
						$dbh->commit;
					}
				}
			}
		}

		######################################################
		####>Now, sending spot information to update java>####
		######################################################
		my $query="UPDATE SPOT SET NAME='$spotName' , UPDATE_DATE=NOW() , UPDATE_USER='$userID'";
		if($externalid){
			$query.=" , EXTERNAL_ID='$externalid'";
		}else{
			$externalid=undef;
		}
		if($intensity){
			$query.=" , INTENSITY=$intensity";
		}else{
			$intensity=undef;
		}
		if($pi){
			$query.=" , ISOELECTRIC_POINT=$pi";
		}else{
			$pi=undef;
		}
		if($pw){
			$query.=" , MOLECULAR_WEIGHT=$pw";
		}else{
			$pw=undef;
		}
		$query.=" WHERE ID_SPOT=$idspot";
		$dbh->do($query);
		$dbh->commit;
		my $protMaj=&getMajProt($sampleID,$dbh);
		
		($pi,$pw,$intensity,$externalid,$protMaj) = map { $_ // ''} ($pi,$pw,$intensity,$externalid,$protMaj);
		
		print header(-type=>'text/plain');
		print "$idspot##$spotName##$x_pos##$y_pos##$pi##$pw##$intensity##$externalid##$protMaj##$sampleID##$sampleName";
	}

	$dbh->disconnect;
	exit;
}

sub getMajProt {
	my ($sampleID,$dbh)=@_;
	my ($majProtAlias)=$dbh->selectrow_array("SELECT ALIAS FROM ANALYSIS_PROTEIN AP,PROTEIN P,ANALYSIS A WHERE AP.ID_PROTEIN=P.ID_PROTEIN AND A.ID_ANALYSIS=AP.ID_ANALYSIS AND A.ID_SAMPLE=$sampleID AND VALID_STATUS>=1 ORDER BY VISIBILITY DESC,NUM_PEP DESC,SCORE DESC,MATCH_GROUP ASC,P.ID_PROTEIN ASC LIMIT 0,1");
	#$majProtAlias=undef unless $majProtAlias;
	return $majProtAlias;
}
#sub getMajProt{
#	my ($sampleID,$dbh)=@_;
#	my @analysis;
#	my $protMaj;
#	if($sampleID>0){
#		my $sthAnalysis=$dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS WHERE ID_SAMPLE=$sampleID");
#		$sthAnalysis->execute;
#		while(my $idanalysis=$sthAnalysis->fetchrow_array){
#			push @analysis, $idanalysis;
#		}
#	}
#	if ($#analysis>=0) {
#		my $query="SELECT SCORE, NUM_PEP, MATCH_GROUP, VISIBILITY, ID_PROTEIN FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS IN (".join(',',@analysis).")";
#		my $sthProteins=$dbh->prepare($query);
#		$sthProteins->execute;
#
#		#Select the best protein between all the candidates
#		my $score=0.0;
#		my $num_pep=0;
#		my $match_group=0;
#		my $visibility=0;
#		my $id_prot=-1;
#		my ($scorec,$num_pepc,$match_groupc,$visibilityc,$id_protc);
#		while(($scorec,$num_pepc,$match_groupc,$visibilityc,$id_protc)=$sthProteins->fetchrow_array){
#			if($scorec>$score){
#				$score=$scorec;
#				$num_pep=$num_pepc;
#				$match_group=$match_groupc;
#				$visibility=$visibilityc;
#				$id_prot=$id_protc;
#			}elsif($scorec==$score){
#				if($num_pepc>$num_pep){
#					$score=$scorec;
#					$num_pep=$num_pepc;
#					$match_group=$match_groupc;
#					$visibility=$visibilityc;
#					$id_prot=$id_protc;
#				}elsif($num_pepc==$num_pep){
#					if($match_groupc<$match_group){
#						$score=$scorec;
#						$num_pep=$num_pepc;
#						$match_group=$match_groupc;
#						$visibility=$visibilityc;
#						$id_prot=$id_protc;
#					}elsif($match_groupc==$match_group){
#						if($visibilityc>$visibility){
#							$score=$scorec;
#							$num_pep=$num_pepc;
#							$match_group=$match_groupc;
#							$visibility=$visibilityc;
#							$id_prot=$id_protc;
#						}
#					}
#				}
#			}
#		}
#		#Get the alias of the protein selected as the majoriaty one
#		($protMaj)=$dbh->selectrow_array("SELECT ALIAS FROM PROTEIN WHERE ID_PROTEIN=$id_prot");
#	}
#	$protMaj=undef unless $protMaj;
#	return $protMaj;
#}

####>Revision history<####
# 1.1.6 Checking undefinition of spot parameters while creating or editing spot (FY 02/09/13)
# 1.1.5 Adding project_x subdirectory in gel URL (FY 02/09/13)
# 1.1.4 Reverse proxy (https) management (PP 17/10/11)
