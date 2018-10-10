#!/usr/local/bin/perl -w

################################################################################
# manageGOFiles.cgi    1.1.3                                                   #
# Authors: P. Poullet, G. Arras & F. Yvon (Institut Curie)                     #
# Contact: myproms@curie.fr                                                    #
# Script for managing annotation and ontology files                            #
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
use File::Copy;
use IO::Uncompress::Gunzip qw(gunzip);
use LWP::UserAgent;
use HTTP::Request;
use File::Listing qw(parse_dir);
use POSIX qw(strftime);
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
if ($userStatus!~/bioinfo|mass|manag/){
    $dbh->disconnect;
    die "Unauthorized access.\n";
}

if (param('submit')) {
    &processForm;
    exit;
}
elsif ($action eq 'checkUpdate') {
	&checkFileUpdate(param('ID'),param('TYPE'));
	exit;
}
elsif ($action eq 'previewGAF'){
	$dbh->disconnect;
    unless(param('ID')){
        die "Missing identifier";
    }
    &previewGAF(param('ID'));
    exit;
}
elsif ($action eq 'upGOA' || $action eq 'upOBO') {
	print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<TITLE>Processing GO File</TITLE>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<BR><BR><BR>
|;
    unless(param('ID')){
        die "Missing identifier";
    }
	if ($action eq 'upOBO'){
		&fetchRemoteFile('OBO',param('ID'));
	}
	else {
		&fetchRemoteFile('GOA',param('ID'));
	}
	$dbh->commit;
    $dbh->disconnect;
    sleep(3);
    print qq
|<SCRIPT language="JavaScript">
window.location="$promsPath{cgi}/manageGOFiles.cgi?ACT=list";
</SCRIPT>
</BODY></HTML>
|;
    exit;
}

elsif($action eq 'list'){

    my $sthOBO = $dbh->prepare("SELECT ID_ONTOLOGY,NAME,STATUS,OBO_FILE,DATA_VERSION,VERSION_DATE,UPDATE_DATE,UPDATE_USER FROM ONTOLOGY WHERE STATUS>=1");
    my $sthAnnot = $dbh->prepare("SELECT ID_GOANNOTATION, GOANNOTATION.ID_SPECIES, GOANNOTATION.NAME, DES, ANNOT_FILE, VERSION_DATE, UPDATE_DATE, UPDATE_USER,
                                SPECIES.COMMON_NAME, SPECIES.SCIENTIFIC_NAME, ID_IDENTIFIER
                                FROM GOANNOTATION, SPECIES
                                WHERE STATUS=1
                                AND SPECIES.ID_SPECIES=GOANNOTATION.ID_SPECIES
                                ORDER BY SPECIES.SCIENTIFIC_NAME,GOANNOTATION.NAME");
    my $sthIdentifier = $dbh->prepare("SELECT NAME FROM IDENTIFIER WHERE ID_IDENTIFIER=?");

    #### HTML ####
    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    my ($light,$dark)=&promsConfig::getRowColors;
    my $bgColor=$light;
    print qq
|<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<TITLE>List of Gene Ontology Files</TITLE>
<SCRIPT LANGUAGE="JavaScript">
function deleteFile(fileType,id,name){
	if(confirm('Delete '+name+' ?')){
		window.location = './manageGOFiles.cgi?ACT=delete&ID='+fileType+':'+id;
	}
}
</SCRIPT>
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
<FONT class="title">List of Gene Ontology Files</FONT>
<BR><BR>
<FONT class="title2">Ontology Files</FONT>&nbsp;&nbsp;<INPUT type="button" class="title3" value="Add new Gene Ontology File" onclick="window.location='$promsPath{cgi}/manageGOFiles.cgi?ACT=add&type=OBO';">
<BR>
<BR>
<TABLE cellspacing=0>
|;
    $sthOBO->execute;
    while(my ($id_ontology,$name,$status,$oboFile,$version,$versionDate,$updateDate,$updateUser) = $sthOBO->fetchrow_array){
		my $scopeText=($status==2)? 'Complete' : 'Slim';
		($version,$versionDate,$updateDate,$updateUser)=&promsMod::chkDef($version,$versionDate,$updateDate,$updateUser);
        print qq
|<TR bgcolor=$bgColor>
<TD>&nbsp;&nbsp;</TD>
<TD width=700><FONT style="font-size:18px;font-weight:bold;">$name [$scopeText]</FONT><BR>
<B>&nbsp;&nbsp;&nbsp;File name: </B>$oboFile<BR>
<B>&nbsp;&nbsp;&nbsp;Version: </B>$version<BR>
<B>&nbsp;&nbsp;&nbsp;File date: </B>$versionDate<BR>
<B>&nbsp;&nbsp;&nbsp;Update date: </B>$updateDate<BR>
</TD>
<TH width=100>
<INPUT type="button" value="Edit" style="width:75px" onclick="window.location='./manageGOFiles.cgi?ACT=edit&ID=OBO:$id_ontology'"><BR>
<INPUT type="button" value="Delete" style="width:75px" onclick="deleteFile('OBO',$id_ontology,'$name')">
</TH>
</TR>
|;
        $bgColor = ($bgColor eq $light)? $dark : $light ;
    }
    $sthOBO->finish;
    print qq
|</TABLE>
<BR>
<BR>
<BR>
<FONT class="title2">Annotation Files</FONT>&nbsp;&nbsp;<INPUT type="button" class="title3" value="Add new Annotation File" onclick="window.location='$promsPath{cgi}/manageGOFiles.cgi?ACT=add&type=GOA';">
<BR>
<BR>
<TABLE cellspacing=0>
|;
    $sthAnnot->execute;
    while(my ($id_goannotation,$id_species,$name,$description,$annotFile,$versionDate,$updateDate,$updateUser,$speciesCommonName,$speciesScientificName,$identifierID) = $sthAnnot->fetchrow_array){
	my $identifierName;
	if($identifierID){
	    $sthIdentifier->execute($identifierID);
	    ($identifierName) = $sthIdentifier->fetchrow_array;
	}
	else {
	    $identifierName = 'Analysis identifier';
	}
	($description,$versionDate,$updateDate,$updateUser)=&promsMod::chkDef($description,$versionDate,$updateDate,$updateUser);
    print qq
|<TR bgcolor=$bgColor>
<TD>&nbsp&nbsp</TD>
<TD width=700><FONT style="font-size:18px;font-weight:bold;">$name</FONT><BR>
<B>&nbsp;&nbsp;&nbsp;Description: </B>$description<BR>
<B>&nbsp;&nbsp;&nbsp;Species: </B>$speciesScientificName ($speciesCommonName)<BR>
<B>&nbsp;&nbsp;&nbsp;File name: </B>$annotFile<BR>
<B>&nbsp;&nbsp;&nbsp;File date: </B>$versionDate<BR>
<B>&nbsp;&nbsp;&nbsp;Indentifier type: </B>$identifierName<BR>
<B>&nbsp;&nbsp;&nbsp;Update date: </B>$updateDate<BR>
</TD>
<TH width=100>
<INPUT type="button" value="Edit" style="width:75px" onclick="window.location='./manageGOFiles.cgi?ACT=edit&ID=GOA:$id_goannotation'"><BR>
<INPUT type="button" value="Delete" style="width:75px" onclick="deleteFile('GOA',$id_goannotation,'$name')">
</TH>
</TR>
|;
        $bgColor = ($bgColor eq $light)? $dark : $light ;
    }
    $sthAnnot->finish;
    $sthIdentifier->finish;
    $dbh->disconnect;
    print qq|
</TABLE>
<BR>
</CENTER>
|;
}
elsif($action eq 'edit'){
    my ($fileType,$id) = split(/:/,param('ID'));
    my @fields;
    my ($name,$isFTP);
    if ($fileType eq 'OBO') {
        ($name,my $oboFile, my $dataVersion,my $versionDate,my $updateDate,my $updateUser,my $status) = $dbh->selectrow_array("SELECT NAME,OBO_FILE,DATA_VERSION,VERSION_DATE,UPDATE_DATE,UPDATE_USER, STATUS FROM ONTOLOGY WHERE ID_ONTOLOGY=$id");
		($dataVersion,$versionDate,$updateDate,$updateUser)=&promsMod::chkDef($dataVersion,$versionDate,$updateDate,$updateUser);
		# Check if local file or URL #
        my $fileString;
        my $scope = ($status == 2)? 'Complete' : 'Slim';
        if ($oboFile =~ /^(ftp|http):\/\//) { # url
			$isFTP=1;
			$fileString="<DIV id=\"fileUpdate\"><B>Checking for update...</B> <IMG src=\"$promsPath{images}/scrollbarGreen.gif\"></DIV>";
        }
		else { # local file
            $fileString = "<INPUT type=\"file\" name=\"localFile\" size=\"60\">";
        }
		my ($chkComplete,$chkSlim)=($status == 2)? ('checked','') : ('','checked');
        push @fields, (
                       ["Name", "<INPUT type=\"text\" name=\"name\" value=\"$name\" size=\"60\">"],
                       ["File", "$oboFile<BR>$fileString"],
                       ["Data version", $dataVersion],
                       ["Version date", $versionDate],
					   ["Scope","<INPUT type=\"radio\" name=\"scope\" value=\"complete\" $chkComplete>&nbsp;Complete&nbsp;<INPUT type=\"radio\" name=\"scope\" value=\"slim\" $chkSlim>&nbsp;Slim"]
                       );
    }
    elsif($fileType eq 'GOA') {
        (my $idSpecies,my $idIdentifier,$name,my $description,my $goaFile,my $versionDate,my $updateDate,my $updateUser) =
		$dbh->selectrow_array("SELECT ID_SPECIES,ID_IDENTIFIER,NAME,DES,ANNOT_FILE,VERSION_DATE,UPDATE_DATE,UPDATE_USER
									FROM GOANNOTATION
									WHERE ID_GOANNOTATION=$id");
		$idIdentifier = 0 unless $idIdentifier;

		my $speciesString = getSelectSpeciesHTML($dbh, $idSpecies);
		my $identifierString = getSelectIdentifierHTML($dbh, $idIdentifier);

        my $fileString = '';
        if($goaFile =~ /^ftp:\/\//){
			$isFTP=1;
			$fileString="<DIV id=\"fileUpdate\"><B>Checking for update...</B> <IMG src=\"$promsPath{images}/scrollbarGreen.gif\"></DIV>";
        }
		else { # local file
            $fileString = "<INPUT type=\"file\" name=\"localFile\" size=\"60\">";
        }

		my $fileContentString = qq
|<INPUT type="button" value="Preview" onclick="previewGAF();">
<DIV id="fileContent"></DIV>
|;
        push @fields, (
            ["Name", "<INPUT type='text' name='name' value=\"$name\">"],
            ["Description", "<TEXTAREA name='description' cols=60>$description</TEXTAREA>"],
			["Species", $speciesString],
            ["File", "$goaFile&nbsp;$fileString"],
			["Identifier used",$identifierString],
            ["Version date", $versionDate],
			["File content", $fileContentString]
        )
    }
    $dbh->disconnect;

    #### HTML ####
    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    my ($light,$dark)=&promsConfig::getRowColors;
    print qq
|<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<TITLE>Edit $name</TITLE>
<SCRIPT LANGUAGE="JavaScript">
var XHR=null;
|;
	if ($isFTP) {
		print qq
|window.onload=function() {
	var updateDiv = document.getElementById('fileUpdate');
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
    var URL = "$promsPath{cgi}/manageGOFiles.cgi?ID=$id&TYPE=$fileType&ACT=checkUpdate";
    XHR.open("GET",URL,true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState==4 && XHR.responseText){
			updateDiv.innerHTML = XHR.responseText;
        }
    }
    XHR.send(null);
}
|;
	}
	print qq
|function previewGAF(){
    var contentDiv = document.getElementById('fileContent');

    //wait image
    contentDiv.innerHTML='<IMG src="$promsPath{images}/scrollbarGreen.gif">';

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
    var URL = "$promsPath{cgi}/manageGOFiles.cgi?ID=$id&ACT=previewGAF";
    XHR.open("GET",URL,true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState==4 && XHR.responseText){
	contentDiv.innerHTML = XHR.responseText;
	//contentDiv.scrollIntoView();
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
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<BR>
<CENTER>
<FONT class="title">Edit $name</FONT>
<BR>
<BR>
<FORM name="editForm" method="post" enctype="multipart/form-data">
<INPUT type="hidden" name="id" value=$id>
<INPUT type="hidden" name="type" value=$fileType>
<TABLE border=0 cellpadding=2 bgcolor=$dark>
|;
    foreach my $row (@fields){
        print "<TR><TH align=right valign=top>&nbsp;",$row->[0]," :</TH><TD bgcolor=$light  width=\"600\">",$row->[1],"</TD></TH>\n";
    }
    print qq
|<TR>
<TH nowrap colspan=2><INPUT type="submit" name="submit" value=" Save "/>
		&nbsp;&nbsp;&nbsp;<INPUT type="reset" value="  Clear  " />
		&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="window.location='./manageGOFiles.cgi?ACT=list';">
</TH>
</TR>
</TABLE>
</FORM>
</BODY>
</HTML>
|;
}
elsif($action eq 'add'){
    my @fields;
    my $title;
    my $jsString = '';
    my $fileType = param('type');
    if($fileType eq 'OBO'){
        push @fields, (
            ["Name", "<INPUT type=\"text\" name=\"name\" size=\"60\">"],
            ["File", "<INPUT type=\"radio\" name=\"fileSource\" value=\"local\" onclick=\"selectFileSource(this);\"><B>Use a local file:</B><BR><INPUT type=\"file\" name=\"localFile\" size=\"60\" disabled><BR>
            <INPUT type=\"radio\" name=\"fileSource\" value=\"remote\" onclick=\"selectFileSource(this);\"><B>Use a remote file</B> (FTP/HTTP URL - .gz accepted):&nbsp;&nbsp;<I><A href=\"ftp://ftp.geneontology.org/pub/go/ontology\" target=\"_blank\">Download link</A> - <A href=\"http://www.geneontology.org\" target=\"_blank\">Info link</A></I><BR><INPUT type=\"text\" name=\"url\" size=\"80\" disabled>"],
            ["Scope","<INPUT type=\"radio\" name=\"scope\" value=\"complete\"><B>Complete</B>&nbsp;<INPUT type=\"radio\" name=\"scope\" value=\"slim\"><B>Slim</B>"]
        );
        $title = 'Add a new Gene Ontology File';
        $jsString .= qq
|function checkForm(myForm){
	if(!myForm.name.value){
		alert('Enter a name for Ontology');
		return false;
	}
	if (!myForm.fileSource[0].checked && !myForm.fileSource[1].checked) {
		alert('Provide a file for the Ontology');
		return false;
	}
	if (myForm.fileSource[0].checked && !myForm.localFile.value){
		alert('Select a local file for Ontology');
		return false;
	}
	if (myForm.fileSource[1].checked && !myForm.url.value){
		alert('Enter a valid URL for Ontology file');
		return false;
	}
	if (!myForm.scope[0].checked && !myForm.scope[1].checked) {
		alert('Define the scope of the Ontology (Complete or Slim?)');
		return false;
	}
	return true;
}
function selectFileSource(radioButton){
	var myForm = document.addForm;
	if(radioButton.value == 'local'){
		myForm.localFile.disabled = false;
		myForm.url.disabled = true;
	}
	else if(radioButton.value == 'remote'){
		myForm.localFile.disabled = true;
		myForm.url.disabled = false;
	}
}
|;
    }
    elsif ($fileType eq 'GOA'){
        # Prepare species file selection #
        #my $goaSelectString = "<SELECT name='species'>\n";
        #my $sthGetSpecies = $dbh->prepare("SELECT ID_SPECIES, COMMON_NAME, SCIENTIFIC_NAME FROM SPECIES");
        #$sthGetSpecies->execute;
        #while(my ($idSpecies,$comName,$sciName) = $sthGetSpecies->fetchrow_array){
        #        $goaSelectString .= "<OPTION value=$idSpecies>$sciName ($comName)</OPTION>\n";
        #}
        #$goaSelectString .= "</SELECT>\n";
        my $speciesString = &getSelectSpeciesHTML($dbh);
		my $identifierString = &getSelectIdentifierHTML($dbh);

        push @fields, (
            ["Name", "<INPUT type='text' name='name'>"],
            ["Description", "<TEXTAREA name='description' cols=60></TEXTAREA>"],
			["Species", "$speciesString&nbsp;<I>(Only reference species can be selected)</I>"],
            ["File", "<INPUT type=\"radio\" name=\"fileSource\" value=\"local\" onclick=\"selectFileSource(this);\"><B>Use a local file:</B><BR><INPUT type=\"file\" name='localFile' size=\"60\" disabled><BR>
            <INPUT type=\"radio\" name=\"fileSource\" value=\"remote\" onclick=\"selectFileSource(this);\"><B>Use a remote file</B> (FTP/HTTP URL - .gz accepted):<BR><INPUT type=\"text\" name=\"url\" size=\"80\" disabled><BR>
			<I>(e.g. <A href=\"ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN\" target=\"_blank\">ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/gene_association.goa_human.gz</A>)</I>"],
			["Identifier used", $identifierString]
        );
        $title = 'Add a new Annotation File';
        $jsString = qq
|function getRadioValue(radio){
	var value;
	var len = radio.length;
	for(i = 0;i<len;i++){
		if(radio[i].checked){
			value = radio[i].value;
		}
	}
	return value;
}
function selectFileSource(radioButton){
	var myForm = document.addForm;
	if(radioButton.value == 'local'){
		myForm.localFile.disabled = false;
		myForm.url.disabled = true;
	} else if(radioButton.value == 'remote'){
		myForm.localFile.disabled = true;
		myForm.url.disabled = false;
	}
}
function checkForm(myForm){
	if(!myForm.name.value){
		alert('Enter a name for GO Annotation');
		return false;
	}
	if (!myForm.fileSource[0].checked && !myForm.fileSource[1].checked) {
		alert('Provide a file for the Annotation');
		return false;
	}
	if (myForm.fileSource[0].checked && !myForm.localFile.value){
		alert('Select a local file for the Annotation');
		return false;
	}
	if (myForm.fileSource[1].checked && !myForm.url.value){
		alert('Enter a valid URL for Annotation file');
		return false;
	}
	if (myForm.identifier.value == -1) {
		alert('Select an identifier type.');
		return false;
	}
	return true;
}
|;
    }
    my ($light,$dark)=&promsConfig::getRowColors;

    # Starting HTML #
    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<TITLE>$title</TITLE>
<SCRIPT LANGUAGE="JavaScript">
$jsString
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<BR>
<CENTER>
<FONT class="title">$title</FONT>
<BR>
<BR>
<FORM name='addForm' method="post" enctype="multipart/form-data" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="type" value=$fileType>
<TABLE border=0 cellpadding=2 bgcolor=$dark>
|;
    foreach my $row (@fields){
        print "<TR><TH align=right valign=top>&nbsp;",$row->[0]," :</TH><TD bgcolor=$light width=\"600\">",$row->[1],"</TD></TH>\n";
    }
    print qq
|<TR>
<TH nowrap colspan=2><INPUT type="submit" name="submit" value=" Add "/>
		&nbsp;&nbsp;&nbsp;<INPUT type="reset" value="  Clear  "/>
		&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="window.location='./manageGOFiles.cgi?ACT=list';">
</TH>
</TR>
</TABLE>
</FORM>
</BODY>
</HTML>
|;
}
elsif($action eq 'delete'){

	my ($fileType,$id) = split(':',param('ID'));
    if ($fileType eq 'OBO') {

		my ($isUsed)=$dbh->selectrow_array("SELECT COUNT(*) FROM GO_ANALYSIS WHERE ID_ONTOLOGY = $id");
        if ($isUsed) {
			$dbh->do("UPDATE ONTOLOGY SET STATUS=0 WHERE ID_ONTOLOGY =$id");
		}
		else {$dbh->do("DELETE FROM ONTOLOGY WHERE ID_ONTOLOGY=$id");}
    }
	elsif ($fileType eq 'GOA') {
        my ($isUsed)=$dbh->selectrow_array("SELECT COUNT(*) FROM GO_ANALYSIS WHERE ID_GOANNOTATION=$id");
        if ($isUsed) {$dbh->do("UPDATE GOANNOTATION SET STATUS=0 WHERE ID_GOANNOTATION=$id");}
		else {$dbh->do("DELETE FROM GOANNOTATION WHERE ID_GOANNOTATION=$id");}
    }

    my $dir=lc($fileType);
    unlink "$promsPath{$dir}/$id.$dir";

    $dbh->commit;
    $dbh->disconnect;

    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<BODY onload="window.location='./manageGOFiles.cgi?ACT=list';"></BODY>
</HTML>
|;
}

sub processForm { #  edit or add processing

	print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<TITLE>Processing GO File</TITLE>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<BR><BR><BR>
|;
    my $name = param('name');
    my $id = (param('id'))?param('id') : 0; # no id means new file
    if(param('type') eq 'OBO') {
        my $status = (param('scope') eq 'complete')? 2 : 1;
        if (param('localFile')) {
			print "<FONT class=\"title2\">Processing local file...";
            mkdir $promsPath{'obo'} unless -d $promsPath{'obo'};
            # Fetching local file #
            $CGITempFile::TMPDIRECTORY=$promsPath{'obo'}; # folder where CGI temporary stores file being uploaded
            my $fileName = (split(/[\\\/]/,param('localFile')))[-1];
            my $tmpFile = tmpFileName(upload('localFile'));

			if ($fileName =~ /\.gz$/) {
				print "[Inflating file...";
				$fileName =~ s/\.gz$//;
				gunzip $tmpFile => "$promsPath{obo}/$fileName" or die "gzip error : $!";
				unlink $tmpFile;
				$tmpFile="$promsPath{obo}/$fileName";
				print " Done.]";
			}

			print " Done.</FONT><BR><BR><FONT class=\"title2\">Updating database...\n";

            # Parsing data version and date from this file #
            my ($dataVersion, $versionDate) = &getOBOVersion($tmpFile);

            my $sthNew = $dbh->prepare("INSERT INTO ONTOLOGY(NAME, OBO_FILE, DATA_VERSION, VERSION_DATE, UPDATE_DATE, UPDATE_USER,STATUS) VALUES (?,?,?,?,NOW(),?,?)");
            #my $status;
            #if(param('scope')){
            #    $status = (param('scope') eq 'complete')? 2 :1;
            #} else {
            #    $status = 1;
            #}
            $sthNew->execute($name,$fileName,$dataVersion,$versionDate,$userID,$status) || die $!;
            $sthNew->finish;

            if($id){
                # Updating database #
                my $sthUpdate = $dbh->prepare("UPDATE ONTOLOGY SET STATUS=0 WHERE ID_ONTOLOGY=?");
                $sthUpdate->execute($id) || die $!;
                $sthUpdate->finish;
                # Removing old file #
                unlink "$promsPath{obo}/$id.obo";
            }
            my $newID = $dbh->last_insert_id(undef,undef,'ONTOLOGY','ID_ONTOLOGY');
            move($tmpFile,"$promsPath{obo}/$newID.obo");
			print " Done.</FONT><BR>\n";
        }
        elsif (!$id) { # no id = new entry, and no localfile means URL was provided
            &fetchRemoteFile('OBO',0,param('url'));
        }
        else { # no local file upload or no URL change (can only be in edit context)
            my $sthUpdate = $dbh->prepare("UPDATE ONTOLOGY SET NAME=?,STATUS=?,UPDATE_DATE=NOW(),UPDATE_USER=? WHERE ID_ONTOLOGY=?");
            $sthUpdate->execute($name,$status,$userID,$id);
            $sthUpdate->finish;
        }
    }
    elsif (param('type') eq 'GOA') {
        my $description = param('description');
		my $idIdentifier = param("identifier");
		my $idSpecies = param('species');
        if (param('localFile')) {
            print "<FONT class=\"title2\">Processing local file...";
            mkdir $promsPath{'goa'} unless -d $promsPath{'goa'};
            # Importing and managing files
            $CGITempFile::TMPDIRECTORY=$promsPath{'goa'};
            my $fileName = (split(/[\\\/]/,param('localFile')))[-1];
            my $tmpFile = tmpFileName(upload('localFile'));
            if($id){
                # Database update
                my $sthUpPrevEntry = $dbh->prepare("UPDATE GOANNOTATION SET STATUS=0 WHERE ID_GOANNOTATION=?");
                $sthUpPrevEntry->execute($id);
                $sthUpPrevEntry->finish;
                unlink "$promsPath{goa}/$id.goa";
            }
			if ($fileName =~ /\.gz$/) {
				print "[Inflating file...";
				$fileName =~ s/\.gz$//;
				gunzip $tmpFile => "$promsPath{goa}/$fileName" or die "gzip error : $!";
				unlink $tmpFile;
				print " Done]";
				$tmpFile="$promsPath{goa}/$fileName";
			}

            print " Done.</FONT><BR><BR><FONT class=\"title2\">Updating database...\n";

            my $versionDate = (param('versionDate'))? param('versionDate') : strftime('%F',localtime);
			$idIdentifier=undef unless $idIdentifier; # 0 -> undef because of foreign key constraint
            my $sthNewEntry = $dbh->prepare("INSERT INTO GOANNOTATION(ID_SPECIES,ID_IDENTIFIER,NAME,DES,ANNOT_FILE,UPDATE_DATE,UPDATE_USER,VERSION_DATE,STATUS) VALUES (?,?,?,?,?,NOW(),?,?,1)");
            $sthNewEntry->execute($idSpecies,$idIdentifier,$name,$description,$fileName,$userID,$versionDate) || die $dbh->errstr;
            $sthNewEntry->finish;
            my $newID = $dbh->last_insert_id(undef,undef,'GOANNOTATION','ID_GOANNOTATION');
			move($tmpFile,"$promsPath{goa}/$newID.goa");
			print " Done.</FONT><BR>\n";
        }
		elsif (!$id){  # no id = new entry, and no localfile means URL was provided
            &fetchRemoteFile('GOA',0,param('url'));
        }
        else { # no local file
			# update GOA entry without local file update
			my $sthUpEntry = $dbh->prepare("UPDATE GOANNOTATION SET ID_SPECIES=?,ID_IDENTIFIER=?,NAME=?,DES=?,UPDATE_USER=?,UPDATE_DATE=NOW() WHERE ID_GOANNOTATION=?");
			if ($idIdentifier) {
				$sthUpEntry->execute($idSpecies,$idIdentifier,$name,$description,$userID,$id);
			}
			else {
				$sthUpEntry->execute($idSpecies,undef,$name,$description,$userID,$id);
			}
			$sthUpEntry->finish;
        }
    }
    $dbh->commit;
    $dbh->disconnect;
    sleep(3);
    print qq
|<SCRIPT language="JavaScript">
window.location="$promsPath{cgi}/manageGOFiles.cgi?ACT=list";
</SCRIPT>
</BODY></HTML>
|;
	exit;
}

sub fetchRemoteFile{
    # Id OR URL must be defined: URL and no id => new entry, Id but no URL => entry update
    my ($fileType,$id,$URL) = @_;
    if(!$id and !$URL){ die "Missing ID or URL"};

    # Fetching URL of an existing entry (ID but no URL provided)
    my ($name,$description,$idSpecies,$status,$idIdentifier);
    if(!$URL){
        if($fileType eq 'OBO'){
            ($name,$URL,$status) = $dbh->selectrow_array("SELECT NAME,OBO_FILE,STATUS FROM ONTOLOGY WHERE ID_ONTOLOGY=$id");
        }
        elsif ($fileType eq 'GOA'){
            ($name,$description,$idSpecies,$idIdentifier,$URL) = $dbh->selectrow_array("SELECT NAME,DES,ID_SPECIES,ID_IDENTIFIER,ANNOT_FILE FROM GOANNOTATION WHERE ID_GOANNOTATION=$id");
        }
    }
	else {
        $name = param('name');
        if($fileType eq 'GOA'){
            $description = param('description')? param('description'): '';
			$idIdentifier = param('identifier');
            $idSpecies = param('species') or die "Missing species identifier";
        }
		else {
            $status = (param('scope') eq 'complete')? 2 : 1;
        }
    }
    print "<BR><BR><FONT class=\"title2\">Downloading $fileType file $URL:</FONT>\n<BR><B>0%";
    # Test if remote file is gunziped #
    my $ext=($URL =~ /\.gz$/)? '.gz' : '';

    #----------------------#
    # Fetching remote file #
    #----------------------#
    my $tmpFileName=strftime("%Y%m%d%H%M%S",localtime).$userID;
	my $goFile="$promsPath{tmp}/$tmpFileName$ext";
	open (GO,">$goFile");
	my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
	$agent->timeout(360);
	if (!$promsPath{'http_proxy'}) {$agent->env_proxy;}
	elsif (lc($promsPath{'http_proxy'}) eq 'no') {$agent->no_proxy;}
	else {$agent->proxy(['http','ftp'],$promsPath{'http_proxy'});}
	my $prevPc=0;
	my $thPc=10;
	my $expectedLength;
	my $bytesReceived=0;
	my $response=$agent->request(HTTP::Request->new(GET => $URL),
					sub {
						my($chunk, $res) = @_;
						$bytesReceived += length($chunk);
						unless (defined $expectedLength) {
						   $expectedLength = $res->content_length || 0;
						}
						if ($expectedLength) {
							my $currPc=int(100 * $bytesReceived / $expectedLength);
							if ($currPc > $prevPc) {
								foreach my $pc ($prevPc+1..$currPc) {
									if ($pc == $thPc) {
										print "$thPc%";
										$thPc+=10;
									}
									else {print '.';}
								}
								$prevPc=$currPc;
							}
						}
						print GO $chunk;
					});

	close GO;

	if ($thPc <= 100) { # download was incolplete (complete= 110%!)
		unlink $goFile;
		$dbh->rollback;
		$dbh->disconnect;
		print "FAILED !</B><BR>\n",$response->status_line;
		print qq
|</B><BR><BR><CENTER><INPUT type="button" value="Go back" onclick="history.back();"></CENTER>
</BODY>
</HTML>
|;
		exit;
    }

	print "</B>\n<FONT class=\"title2\"> Done.</FONT><BR>\n" ;
	if ($ext eq '.gz') {
		print "<BR><FONT class=\"title2\">Inflating file...";
		gunzip $goFile => "$promsPath{tmp}/$tmpFileName" or die "gzip error : $!";
		unlink $goFile;
		$goFile="$promsPath{tmp}/$tmpFileName";
		print " Done.</FONT><BR>\n";
	}

    #-------------------#
    # Updating database #
    #-------------------#
    print "<BR><FONT class=\"title2\">Updating database...";

    #my $fileName = (split('/',$URL))[-1];
    #$fileName =~ s/\.gz$//;
    my $fileName = $URL;

    if($fileType eq 'OBO'){
        # Fetching data #
        my ($dataVersion,$versionDate) = &getOBOVersion($goFile); #
        my ($FTPFileDate,$response) = &getFTPFileDate($URL);

		my $addEntry=1;
		my $newID;
		if ($id) {
			($addEntry)=$dbh->selectrow_array("SELECT COUNT(*) FROM GO_ANALYSIS WHERE ID_ONTOLOGY=$id");
		}
        # New entry #
		if ($addEntry) {
			my $sthAddOBO = $dbh->prepare("INSERT INTO ONTOLOGY (NAME,OBO_FILE,DATA_VERSION,VERSION_DATE,UPDATE_DATE,UPDATE_USER,STATUS) VALUES (?,?,?,?,NOW(),?,?)");
			$sthAddOBO->execute($name,$URL,$dataVersion,$FTPFileDate,$userID,$status);
			$sthAddOBO->finish;

			# Moving file #
			$newID = $dbh->last_insert_id(undef,undef,'ONTOLOGY','ID_ONTOLOGY');
		}
		else {$newID=$id;} # $id is defined in this case

        # Updating previous entry #
        if ($id) {
			if ($addEntry) { # disable old entry
				$dbh->do("UPDATE ONTOLOGY SET STATUS=0 WHERE ID_ONTOLOGY=$id");
			}
			else { # entry not yet used => update
				my $sthUpOBO = $dbh->prepare("UPDATE ONTOLOGY SET DATA_VERSION=?,VERSION_DATE=?,UPDATE_DATE=NOW(),UPDATE_USER=? WHERE ID_ONTOLOGY=$id");
				$sthUpOBO->execute($dataVersion,$FTPFileDate,$userID);
				$sthUpOBO->finish;
			}
            unlink "$promsPath{obo}/$id.obo";
        }

		move($goFile,"$promsPath{obo}/$newID.obo");
    }
    elsif ($fileType eq 'GOA') {
        my ($fileDate,$response) = &getFTPFileDate($URL);

		my $addEntry=1;
		my $newID;
		if ($id) {
			($addEntry)=$dbh->selectrow_array("SELECT COUNT(*) FROM GO_ANALYSIS WHERE ID_GOANNOTATION=$id");
		}
		# New entry #
		if ($addEntry) {
			my $sthAddGOA = $dbh->prepare("INSERT INTO GOANNOTATION(ID_SPECIES,ID_IDENTIFIER,NAME,DES,ANNOT_FILE,VERSION_DATE,UPDATE_USER,UPDATE_DATE,STATUS) VALUES (?,?,?,?,?,?,?,NOW(),1)");
			$idIdentifier=undef unless $idIdentifier; # 0 -> undef because of foreign key constraint
			$sthAddGOA->execute($idSpecies,$idIdentifier,$name,$description,$fileName,$fileDate,$userID);
			$sthAddGOA->finish;
			$newID = $dbh->last_insert_id(undef,undef,'GOANNOTATION','ID_GOANNOTATION');
		}
		else {$newID=$id;} # $id is defined in this case

		# Updating previous entry #
        if ($id) {
			if ($addEntry) { # disable old entry
				$dbh->do("UPDATE GOANNOTATION SET STATUS=0 WHERE ID_GOANNOTATION=$id");
			}
			else { # entry not yet used => update
				my $sthUpGOA = $dbh->prepare("UPDATE GOANNOTATION SET VERSION_DATE=?,UPDATE_DATE=NOW(),UPDATE_USER=? WHERE ID_GOANNOTATION=$id");
				$sthUpGOA->execute($fileDate,$userID);
				$sthUpGOA->finish;
			}
            unlink "$promsPath{goa}/$id.goa";
        }

        move($goFile,"$promsPath{goa}/$newID.goa");
    }

	print " Done.</FONT>\n";

}

sub getFTPFileDate{
    # This function fetch the last-modif date from a file located in a FTP directory
    # Return the date and the response object of the FTP request
    # Return an undefined date value if the file is not found or if the request fails
	my $URL=$_[0];
    my @splittedURL = split('/',$URL);
    my $fileName = pop @splittedURL;
    my $dirURL = join('/',@splittedURL);
    my $fileDate;

	if ($URL=~/^http:\/\/www\.geneontology\.org/) { # convert http path into ftp
		if ($URL=~/GO_slims/i) { # slim (http://www.geneontology.org/GO_slims/goslim_generic.obo)
			$dirURL='ftp://ftp.geneontology.org/pub/go/ontology/subsets';
		}
		else { # complete (http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo)
			my $versionDir=pop @splittedURL;
			$dirURL="ftp://ftp.geneontology.org/pub/go/ontology/$versionDir";
		}
	}
#print header(-'content-encoding'=>'no',-charset=>'UTF-8',-type=>'text/plain');
    my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
    $agent->timeout(90);
	if (!$promsPath{'http_proxy'}) {$agent->env_proxy;}
	elsif (lc($promsPath{'http_proxy'}) eq 'no') {$agent->no_proxy;}
	else {$agent->proxy(['http','ftp'],$promsPath{'http_proxy'});}
    my $req = HTTP::Request->new(GET => $dirURL);
	#$req->content_type("text/ftp-dir-listing");
    my $response = $agent->request($req);
	#while (my ($key,$value)=each %{$response->headers}) {print "$key='$value' =====<BR>\n";}
    #$response->header('Content-Type' => "text/ftp-dir-listing");
	#$response->header(Accept => "text/plain");
#print "content-type:",$response->headers->{'content-type'},"<BR>----------------<BR>\n";
    if ($response->is_success){
		if ($response->headers->{'content-type'} eq 'text/html') {
			(my $textResponse=$response->content)=~s/<[^>]+>//g; # remove HTML tags
			foreach my $fileLine (split(/\n/,$textResponse)){
				$fileLine=~s/^\s+//;
				next unless $fileLine=~/^\d/;
				my @fileInfo=split(/\s+/,$fileLine); # expected format: size month day hour:min/year file
				if ($fileInfo[-1] eq $fileName) {
					my %months=('Jan'=>1,'Fev'=>2,'Mar'=>3,'Apr'=>4,'May'=>5,'Jun'=>6,'Jul'=>7,'Aug'=>8,'Sep'=>9,'Oct'=>10,'Nov'=>11,'Dec'=>12);
					$fileDate=($fileInfo[3]=~/:/)? strftime('%Y',localtime) : $fileInfo[3];
					$fileDate.="-$months{$fileInfo[1]}-$fileInfo[2]";
#print "DATE=$fileDate\n";
				}
			}
		}
		else { # assumes 'text/ftp-dir-listing' unix format: -rw-r--r--    1 82       82          70396 Oct 22 07:10 file
			foreach my $fileInfo (parse_dir($response->content)){ #,undef,'apache' 'unix''apache'
				if ($fileInfo->[0] eq $fileName){
					$fileDate = strftime('%F',localtime($fileInfo->[3]));
					last;
			   }
			}
		}
    }
	else {
		print "\n#ERROR: ",$response->status_line,"#<BR>\n";
	}
    return ($fileDate,$response);
}

sub getOBOVersion{
    my $fileName = shift;
    my $dataVersion;
    my $versionDate;

    open(OBO, $fileName) || die "There was an error during file reading: $! $fileName\n";
    while(<OBO>){
        if(/^data-version:.*(\d.{8}\d)/){ # /data-version: ([^\n\r]+)/
            $dataVersion = $1;
        }
		elsif (/^date: ([^\s]+)/){
            # fitting with mySQL date format
            my ($d,$m,$y) = split(':',$1);
            $versionDate = "$y-$m-$d";
        }
        last if /^\[Term\]/;
    }
    close(OBO);

    $dataVersion = 'N/A' unless $dataVersion;
    $versionDate = '' unless $versionDate;

    return ($dataVersion,$versionDate);
}

sub getSelectSpeciesHTML{
    my ($dbh, $selectedId) = @_;

    $selectedId = 0 unless $selectedId;

    my $speciesString = "<SELECT name=\"species\">\n";
    my $sthSpecies = $dbh->prepare("SELECT ID_SPECIES, SCIENTIFIC_NAME, COMMON_NAME FROM SPECIES WHERE IS_REFERENCE=1 ORDER BY SCIENTIFIC_NAME");
    $sthSpecies->execute;
    while(my ($idSpecies,$sciName,$comName) = $sthSpecies->fetchrow_array){
	my $selected = ($selectedId == $idSpecies)? 'selected':'';
	$speciesString .= "<OPTION value=$idSpecies $selected>$sciName ($comName)</OPTION>\n";
    }
    $sthSpecies->finish;
    $speciesString .= "</SELECT>\n";

    return $speciesString;
}

sub getSelectIdentifierHTML{
    my($dbh, $idIdentifier) = @_;

    #$idIdentifier = 0 unless $idIdentifier;

    my $identifierString = "<SELECT name=\"identifier\">\n";
    $identifierString .= "<OPTION value=\"-1\">-= Select an identifier type =-</OPTION>\n" if $action eq 'add';
    my $uniprotString;
    my $otherString;
    my $sthIdentifier = $dbh->prepare("SELECT ID_IDENTIFIER, NAME FROM IDENTIFIER ORDER BY NAME");
    $sthIdentifier->execute;
    while(my ($id, $name) = $sthIdentifier->fetchrow_array){
		my $selected = (defined $idIdentifier && $id == $idIdentifier)? 'selected' : '';
		my $option = "<OPTION value=\"$id\" $selected>$name</OPTION>\n";
		if ($name =~ /uniprot/i) {
			$uniprotString .= $option;
		}
		else {
			$otherString .= $option;
		}
    }
    $sthIdentifier->finish;


	if ($uniprotString) {
		#$identifierString .= "<OPTION disabled>--- Uniprot ---</OPTION>\n$uniprotString";
		$identifierString .= "<OPTGROUP label=\"UniProt:\">\n$uniprotString</OPTGROUP>\n";
    }
    my $selected = (defined $idIdentifier && $idIdentifier == 0)? 'selected' : '';
	$identifierString .= "<OPTGROUP label=\"Unspecified:\">\n<OPTION value=\"0\" $selected>Analysis identifier</OPTION>\n</OPTGROUP>\n";
    if ($otherString) {
		#$identifierString .= "<OPTION disabled>--- Other ---</OPTION>\n$otherString";
		$identifierString .= "<OPTGROUP label=\"Other:\">\n$otherString</OPTGROUP>\n";
    }
    #my $selected = (defined $idIdentifier && $idIdentifier == 0)? 'selected' : '';
    #$identifierString .= "<OPTION disabled>--- Unspecified ---</OPTION>\n<OPTION value=\"0\" $selected>Analysis identifier</OPTION>\n";

    $identifierString .= '</SELECT>';

    return $identifierString;
}

sub checkFileUpdate {
	my ($id,$fileType)=@_;

	my ($URL,$versionDate)=($fileType eq 'OBO')? $dbh->selectrow_array("SELECT OBO_FILE,VERSION_DATE FROM ONTOLOGY WHERE ID_ONTOLOGY=$id") :
									$dbh->selectrow_array("SELECT ANNOT_FILE,VERSION_DATE FROM GOANNOTATION WHERE ID_GOANNOTATION=$id");
	$dbh->disconnect;

	my $updateString="<INPUT type=\"button\" value=\"Update file\" onclick=\"window.location='$promsPath{cgi}/manageGOFiles.cgi?ACT=up$fileType&ID=$id';\"";
	my ($remoteFileDate,$response) = &getFTPFileDate($URL);
	if ($response->is_success) { # FTP was succesfully connected
		if ($remoteFileDate) { # defined value = file found
			if ($versionDate && $remoteFileDate eq $versionDate) { # remote version date = local version date
				$updateString = "<FONT style=\"font-weight:bold;color:#00AA00\">File is up to date</FONT>";
			}
			else { # remote version is more recent
				$updateString = "<FONT style=\"font-weight:bold;color:#0000DD\">New version available</FONT> <INPUT type=\"button\" value=\"Update file\" onclick=\"window.location='$promsPath{cgi}/manageGOFiles.cgi?ACT=up$fileType&ID=$id';\"/>";
			}
		}
		else { # file was not found
			$updateString = "<FONT style=\"font-weight:bold;color:#DD0000\">File not found on remote server</FONT>"; # at URL: $oboFile";
		}
	}
	else { # connexion to FTP failed
		$updateString = "&nbsp;<FONT style=\"font-weight:bold;color:#DD0000\">Unable to connect to remote server</FONT>";
	}

	print header(-'content-encoding'=>'no',-charset=>'UTF-8',-type=>'text/plain');
    warningsToBrowser(1);
	print $updateString;
}


sub previewGAF {
    my $goaID = shift;

    my $fileName = "$promsPath{goa}/$goaID.goa";

    my @data;

    open GAF, $fileName or die $!;
    while(<GAF>){
		next if /^!/;
		chomp;
		my @lineData = split /\t/;
		my @identifiers = ($lineData[1]);
		my @aliases = split /\|/, $lineData[10];
		foreach my $alias (@aliases){
			push @identifiers, $alias;
		}
		push @data, \@identifiers;
		last;
    }
    close GAF;

    print header(-'content-encoding'=>'no',-charset=>'UTF-8',-type=>'text/plain');
    warningsToBrowser(1);
    print "<B>Examples of identifiers and aliases used:</B><BR>\n";
    foreach my $idsRef (@data) {
		print '- ',join ', ', @$idsRef;
		print "<BR>\n";
    }
}

####>Revision history<####
# 1.1.3 Color on main "Go to" options (PP 25/03/16)
# 1.1.2 Better parsing of OBO version in &getOBOVersion subroutine (PP 18/04/14)
# 1.1.1 Added http support for remote file and gunzip support for local file (PP 21/11/13)
# 1.1.0 Checks define status of potentially uninitialized variables (PP 28/10/13)
# 1.0.9 Minor bug correction for GOA preview (FY 24/10/13)
# 1.0.8 Fixed critical bug when trying to update a GOA file by FTP (FY 21/10/13)
# 1.0.7 Restrict species selection to reference species (FY 16/09/13)
# 1.0.6 Better proxy declaration for LWP (PP 02/07/13)
# 1.0.5 Improved file update/download management & other minor fonctionalites (PP 28/05/13)
# 1.0.4 Reorganizing identifier-type select options for GOA files +<BR>Fixing foreign key constraint error when unselecting an identifier type (FY 23/04/13)
# 1.0.3 GOA files URL are now stored in GOANNOTATION table & add annotation file content preview in Edit context (FY 22/03/13)
# 1.0.2 Called by promsMain.cgi (PP 13/12/12)
# 1.0.1 Add remote file option for obo files<br>- creating GO directories if not exist<br>- Manage complete or slim obo files (FY 25/01/12)
# 1.0.0 New script to manage ontology and GO annotation files (FY 17/11/11)
