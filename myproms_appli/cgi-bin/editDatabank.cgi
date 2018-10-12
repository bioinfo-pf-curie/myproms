#!/usr/local/bin/perl -w

################################################################################
# editDatabank.cgi         2.4.2                                               #
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
use POSIX qw(strftime); # to get the time
use strict;
use LWP::UserAgent;
#use URI::Escape;
use IO::Uncompress::Gunzip qw(gunzip);
use IO::Uncompress::Unzip qw(unzip $UnzipError);
use File::Path qw(rmtree); # remove_tree
use File::Copy;
use promsConfig;
use promsMod;

#print header; # warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my %identifierTypes=&promsMod::getIdentifierTypes;
#$CGITempFile::TMPDIRECTORY=$promsPath{'banks'}; # folder where GCI temporary stores file being uploaded
my %mascotServers=&promsConfig::getMascotServers;

#################
####>Globals<####
#################
my (%numEntry,$totNumEntry);
my %speciesList=(
	'Arabidopsis thaliana'=>1,
	'Bos taurus'=>1,
	'Caenorhabditis elegans'=>1,
	'Chlamydomonas reinhardtii'=>1,
	'Danio rerio'=>1,
	'Dictyostelium discoideum'=>1,
	'Drosophila melanogaster'=>1,
	'Escherichia coli'=>1,
	'Gallus gallus'=>1,
	'Homo sapiens'=>1,
	'Mus musculus'=>1,
	'Mycoplasma pneumoniae'=>1,
	'Oryza sativa'=>1,
	'Plasmodium falciparum'=>1,
	'Pneumocystis carinii'=>1,
	'Rattus norvegicus'=>1,
	'Saccharomyces cerevisiae'=>1,
	'Schizosaccharomyces pombe'=>1,
	'Takifugu rubripes'=>1,
	'Xenopus laevis'=>1,
	'Zea mays'=>1
);

#############################
####>Fetching parameters<####
#############################
my $action=param('ACT');
my $databankID=(param('ID'))? param('ID') : 0; # 0 if ACT=add


##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;


&processForm if param('save');
if ($action eq 'delete') {&deleteDB;}
elsif ($action eq 'rules') {&ajaxTestParseRules;}
elsif ($action eq 'update') {&updateMascotDbFile;}


#######################
####>Fetching data<#### from tables ORGANISM, DATABANK and DATABANK_ORGANISM
#######################
my ($pageTitle,@nameText,@valueText);

####>Fetching list of used species
my $sthSp=$dbh->prepare("SELECT DISTINCT(ORGANISM) FROM DATABANK");
$sthSp->execute;
while (my ($species)=$sthSp->fetchrow_array) {
	next unless $species;
	$speciesList{$species}=1;
}
$sthSp->finish;

####>Fetching databank info
@nameText=('Name','Version Name','Version Date','Description','Type','Identifier type','Internal decoy','Species','Is contaminant DB','Sequence file');
my (%parseString,%defaultIdentType);
#my $mascotDisabStrg='';
#my $serverDisabStrg='';

if ($action eq 'edit') {

	push @nameText,('Number of entries','Comments');

	###>Fetching databank info
	my (@values)=$dbh->selectrow_array("SELECT D.NAME,VERSION_NAME,VERSION_DATE,D.DES,CONCAT(DT.ID_DBTYPE,':#:',DT.NAME),IDENTIFIER_TYPE,DECOY_TAG,ORGANISM,IS_CRAP,FASTA_FILE,NUM_ENTRY,COMMENTS FROM DATABANK D,DATABANK_TYPE DT WHERE ID_DATABANK=$databankID AND D.ID_DBTYPE=DT.ID_DBTYPE");
	$pageTitle="Editing Databank <FONT color='#DD0000'>$values[0]</FONT>";

	####>Formatting HTML code
	push @valueText,"<INPUT type='text' name='name' value='$values[0]' size='50'>"; # name
	$values[1]='' unless $values[1]; # version name
	push @valueText,"<INPUT type='text' name='versionName' value='$values[1]' size='20'>";
	$values[2]=~s/ .*//;  # version date (removing time)
	push @valueText,"<INPUT type='text' name='versionDate' value='$values[2]' size='10'> (yyyy-mm-dd)";
	$values[3]='' unless $values[3]; # description
	push @valueText,"<TEXTAREA name='des' rows='2' cols='65'>$values[3]</TEXTAREA>";
	my ($dbTypeID,$dbTypeName)=split(/:#:/,$values[4]);
	push @valueText,"&nbsp;$dbTypeName <INPUT type='hidden' name='dbType' id='dbType' value='$dbTypeID'/><INPUT type='button' id='defButton' style='font-size:10px' value='Test rules' onclick='ajaxTestParseRules()'/>"; # db type
	# identifier type
	my $identTypeStrg="<SELECT name='identifierType'>\n<OPTION value=''>* Unspecified *</OPTION>\n"; # identifier type
	foreach my $identType (sort{lc($identifierTypes{$a}) cmp lc($identifierTypes{$b})} keys %identifierTypes) {
		$identTypeStrg.="<OPTION value='$identType'";
		$identTypeStrg.=' selected' if ($values[5] && $values[5] eq $identType);
		$identTypeStrg.=">$identifierTypes{$identType}</OPTION>\n";
	}
	$identTypeStrg.="</SELECT>";
	push @valueText,$identTypeStrg;
	# internal decoy
	my ($usedDB)=$dbh->selectrow_array("SELECT COUNT(*) FROM ANALYSIS_DATABANK WHERE ID_DATABANK=$databankID");
	$values[6]='' unless $values[6];
	my $decoyString;
	if ($usedDB) { # db is used cannot edit decoy status
		$decoyString=($values[6])? "Yes&nbsp;<B>Decoy tag:</B> $values[6] <INPUT type=\"hidden\" name=\"decoyTag\" value=\"$values[6]\"/>" : 'No';
	}
	else {
		my ($selDecoy,$spanVis,$disTag)=($values[6])? (' selected','visibile','') : ('','hidden',' disabled');
		$decoyString="<SELECT name=\"hasDecoy\" onchange=\"activateDecoyTag(this.value)\"><OPTION value=\"no\">No</OPTION><OPTION value=\"yes\"$selDecoy>Yes</OPTION></SELECT>";
		$decoyString.="&nbsp;<SPAN id=\"decoySpan\" style=\"visibility:$spanVis\"><B>Decoy tag:</B><INPUT type=\"text\" name=\"decoyTag\" id=\"decoyTag\" size=10 value=\"$values[6]\"$disTag/>&nbsp;<I><B>eg.</B> REV_</I></SPAN>";
	}
	push @valueText,$decoyString;
	#>species
	my $spString="<SELECT name='species' onchange=\"activateNewSpecies(this.value)\">\n<OPTION value=''>Multiple species</OPTION>\n"; # species
	foreach my $species (sort{$a cmp $b} keys %speciesList) {
		my $selStrg=($values[7] && $values[7] eq $species)? ' selected' : '';
		$spString.="<OPTION value='$species'$selStrg>$species</OPTION>\n";
	}
	$spString.="<OPTION value='new'>-= New species =-</OPTION>\n";
	$spString.="</SELECT>\n";
	$spString.="&nbsp;<B>New species:</B><INPUT type='text' name='new_species' value='' size='25' disabled>";
	push @valueText,$spString;

	#>Contaminant
	my $crapInfo=($values[8])? 'Yes' : 'No';
	push @valueText,$crapInfo;

	#>Checking fasta file
	my $fileInfo;
	if ($values[9]=~/:/) { # mascot file
		#my ($mascotServer,$dbankDir,$fileName)=split(':',$values[8]);
		($fileInfo=$values[9])=~s/:/ > /g;
		$fileInfo.=" <INPUT type='button' style='font-size:10px' value='Check for update' onclick='updateMascotDbFile()'/>";
	}
	else {
		$fileInfo=$values[9];
		my $fastaFile="$promsPath{banks}/db_$databankID/$values[9]";
		unless (-e $fastaFile) {$fileInfo.=" <FONT color=#DD0000>(File is missing)</FONT>";}
	}
	push @valueText,"&nbsp;$fileInfo";

	$totNumEntry=$values[10]; # number of entries
	push @valueText,"&nbsp;$values[10]";
	$values[11]='' unless $values[11]; # comments
	push @valueText,"<TEXTAREA name='comments' rows='5' cols='65'>$values[11]</TEXTAREA>";
}
else { # add
	push @nameText,'Comments';
	#$databankID=0;

	####>Formatting HTML code<####
	$pageTitle='Adding new Databank to Server';
	push @valueText,"<INPUT type='text' name='name' value='' size='50'/>"; # name
	push @valueText,"<INPUT type='text' name='versionName' value='' size='20'/>"; # version name
	my $date=strftime("%Y-%m-%d", localtime);
	push @valueText,"<INPUT type='text' name='versionDate' value='$date' size='10'/> (yyyy-mm-dd)"; # version date
	push @valueText,"<TEXTAREA name='des' rows='2' cols='65'></TEXTAREA>"; # description
	# db type
	my $dbTypeString='<TABLE border=0 cellspan=0 cellpadding=0><TR>';
	$dbTypeString.="<TH valign=top><SELECT name='dbType' id='dbType' onchange='changeParseRule(this.value)'><OPTION selected value='0'>-= Choose from list =-</OPTION>";
	$parseString{0}='<BR><BR><BR><BR><BR><BR>';
	$defaultIdentType{0}='-';
	my $sthT=$dbh->prepare("SELECT ID_DBTYPE,NAME,DES,PARSE_STRING,DEF_IDENT_TYPE FROM DATABANK_TYPE ORDER BY LOWER(NAME) ASC");
	$sthT->execute || die $dbh->errstr;
	while (my ($dbTypeID,$name,$des,$string,$defIdentType)=$sthT->fetchrow_array) {
		$des='No description' unless $des;
		$string='<B>No info</B><BR><BR><BR><BR><BR><BR>' unless $string;
		$dbTypeString.="<OPTION value='$dbTypeID'>$name</OPTION>";
		$parseString{$dbTypeID}="<B>Description:</B> $des<BR><B>Parse rules:</B><BR>$string";
		$defaultIdentType{$dbTypeID}=($defIdentType)? $defIdentType : '-';
	}
	$sthT->finish;
	$dbTypeString.="</SELECT><BR><BR><INPUT type=\"button\" id=\"defButton\" style=\"font-size:10px\" value=\"Test rules\" onclick= \"ajaxTestParseRules()\" disabled/></TH>";
	$dbTypeString.="<TD><SPAN id='dbTypeSpan'><BR><BR><BR><BR><BR><BR></SPAN></TD></TR></TABLE>";
	push @valueText,$dbTypeString;
	# identifier type
	my $identTypeStrg="<SELECT name='identifierType'>\n<OPTION value=''>* Unspecified *</OPTION>\n";
	foreach my $identType (sort{lc($identifierTypes{$a}) cmp lc($identifierTypes{$b})} keys %identifierTypes) {
		$identTypeStrg.="<OPTION value='$identType'>$identifierTypes{$identType}</OPTION>\n";
	}
	$identTypeStrg.="</SELECT>";
	push @valueText,$identTypeStrg;
	# internal decoy
	my $decoyString='<SELECT name="hasDecoy" onchange="activateDecoyTag(this.value)"><OPTION value="no">No</OPTION><OPTION value="yes">Yes</OPTION></SELECT>';
	$decoyString.='&nbsp;<SPAN id="decoySpan" style="visibility:hidden"><B>Decoy tag:</B><INPUT type="text" name="decoyTag" id="decoyTag" size=10 value="" disabled/>&nbsp;<I><B>eg.</B> REV_</I></SPAN>';
	push @valueText,$decoyString;
	# species
	my $spString="<SELECT name='species' onchange=\"activateNewSpecies(this.value)\">\n<OPTION value='' selected>* Multiple species *</OPTION>\n"; # species
	foreach my $species (sort{$a cmp $b} keys %speciesList) {
		$spString.="<OPTION value='$species'>$species</OPTION>\n";
	}
	$spString.="<OPTION value='new'>-= New species =-</OPTION>\n";
	$spString.="</SELECT>\n";
	$spString.="&nbsp;<B>New species:</B><INPUT type='text' name='new_species' value='' size='25' disabled>";
	push @valueText,$spString;
	# contaminants
	push @valueText, "<INPUT type=\"checkbox\" name='contaminant' id='contaminant' value=\"1\"/>(check only if this databank is used in MaxQuant as contaminants database)";

	my $mascotDisabStatus=' disabled';
	my $mascotOptionsStrg="<OPTION value=\"\">-= Select =- </OPTION>\n";

	foreach my $mascotServer (sort{lc($a) cmp lc($b)} keys %mascotServers) {
		$mascotOptionsStrg.="<OPTGROUP label='$mascotServer'>\n";
		my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
		$agent->timeout(360);
		if ($mascotServers{$mascotServer}[1]) { # proxy settings
			if ($mascotServers{$mascotServer}[1] eq 'no') {$agent->no_proxy($mascotServers{$mascotServer}[0]);}
			else {$agent->proxy('http', $mascotServers{$mascotServer}[1]);}
		}
		else {$agent->env_proxy;}
		my $response = $agent->get("$mascotServers{$mascotServer}[0]/cgi/myproms4databanks.pl?ACT=list");
		while (my $wait = $response->header('Retry-After')) {
			sleep $wait;
			$response = $agent->get($response->base);
        }
		my @resultLines;
        if ($response->is_success) {
			@resultLines = split("\n",$response->content);
		}
		foreach (@resultLines) {
			next if /^#/;
			chomp;
			if (/^(\S+)\s+(\S+)/) {
				my $dbankName=$1;
				my $fullDbankfile=$2;
				my $dbFileName=(split(/\//,$fullDbankfile))[-1];
				$mascotOptionsStrg.="<OPTION value='$mascotServer##$dbankName##$fullDbankfile'>$dbankName ($dbFileName)</OPTION>\n";
				$mascotDisabStatus='';
			}
		}
		$mascotOptionsStrg.="</OPTGOUP>\n";
	}
	#$mascotDisabStrg=($mascotDisabStatus)? 'var mascotDisabStatus=true' : 'var mascotDisabStatus=false';

	#>Server data dir
	my $sharedDisabStatus=' disabled';
	my $selectShareStrg='Not available'; # default
	if ($promsPath{'shared'} && -e $promsPath{'shared'}) {
		$selectShareStrg=&promsMod::browseDirectory_getFiles($promsPath{'shared'},'sharedFile',{inputType=>'select',fileMatch=>qr/\.(fasta|gz|zip)$/i,returnString=>1,selectAttribStrg=>' id="sharedFile" onchange="updateRulesButton()" disabled'});
		$sharedDisabStatus='' if $selectShareStrg=~/^<SELECT/;
	}
	#$serverDisabStrg=($sharedDisabStatus)? 'var serverDisabStatus=true' : 'var serverDisabStatus=false';
	my $fileString=qq
|<INPUT type="radio" name="fileOption" value="mascotFile" onclick="updateFileOptions(this.value)"$mascotDisabStatus/><B>Use Mascot databank:</B>&nbsp;<SELECT name="mascotFile" id="mascotFile" onchange="updateRulesButton()" disabled>$mascotOptionsStrg</SELECT>
<BR><INPUT type="radio" name="fileOption" value="sharedFile" onclick="updateFileOptions(this.value)"$sharedDisabStatus/><B>Use shared directory:</B>&nbsp;$selectShareStrg
<BR><INPUT type="radio" name="fileOption" value="localFile" onclick="updateFileOptions(this.value)"/><B>Upload a local file:</B><BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="file" name="localFile" size="60" disabled/>
<BR><INPUT type="radio" name="fileOption" value="internetFile" onclick="updateFileOptions(this.value)"/><B>Download a file from Internet:<BR>URL:</B><INPUT type="text" name="internetFile" value="" size="60" disabled/>
<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<FONT style="font-size:11px;font-style:italic">http or ftp path - Normal or compressed (.gz,.zip) FASTA files are accepted.</FONT>
|;
	push @valueText,$fileString; # fasta file
	push @valueText,"<TEXTAREA name='comments' rows='5' cols='65'></TEXTAREA>"; # comments
}

$dbh->disconnect;


#######################
####>Starting HTML<####
#######################
print header(-'content-encoding'=>'no',-charset=>'utf-8'),"\n";
warningsToBrowser(1); # DEBUG
print qq
|<HTML>
<HEAD>
<TITLE>Managing Databanks</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
var action='$action';
var popupWin=null;
|;
if ($action eq 'add') {
	print "var parseString=new Object();\nvar defaultIdentifierType=new Object();\n";
	foreach my $dbTypeID (sort keys %parseString) {
		print "parseString['$dbTypeID']='$parseString{$dbTypeID}';\n";
		print "defaultIdentifierType['$dbTypeID']='$defaultIdentType{$dbTypeID}';\n";
	}
}
print qq
|function activateDecoyTag(decoyStatus) {
	var visSpan,disTag;
	if (decoyStatus=='yes') {
		visSpan='visible';
		disTag=false;
	}
	else {
		visSpan='hidden';
		disTag=true;
	}
	document.getElementById('decoySpan').style.visibility=visSpan;
	document.getElementById('decoyTag').disabled=disTag;
}
function activateNewSpecies(species) {
	var disabStatus;
	if (species=='new') {disabStatus=false;}
	else {disabStatus=true;}
	document.dbForm.new_species.disabled=disabStatus;
}
function changeParseRule(dbTypeID) {
	document.getElementById('dbTypeSpan').innerHTML=parseString[dbTypeID];
	if (!parseString[dbTypeID].match('<BR><BR>')) {
		document.getElementById('dbTypeSpan').innerHTML+='<BR>Strings in <I>italics</I> are optional.';
	}
	var selIdentType=document.dbForm.identifierType;
	if (defaultIdentifierType[dbTypeID]=='-') {
		selIdentType.selectedIndex=0;
	}
	else {
		for (var i=1;i<selIdentType.options.length;i++) {
			if (selIdentType.options[i].value==defaultIdentifierType[dbTypeID]) {
				selIdentType.selectedIndex=i;
				break;
			}
		}
	}
	updateRulesButton();
}
function updateRulesButton() {
	var disStatus=false,
		mascot=document.getElementById('mascotFile'),
		shared=document.getElementById('sharedFile'),
		mascotOK=(!mascot.disabled && mascot.value)? true : false,
		serverOK=(shared && !shared.disabled && shared.value && shared.value.match(/\\.fasta\$/))? true : false;
	if (document.getElementById('dbType').value==0 \|\| (!mascotOK && !serverOK)) {
		disStatus=true;
	}
	document.getElementById('defButton').disabled=disStatus;
	document.getElementById('rulesDIV').innerHTML=''; // clear previous results if any
}
function updateMascotDbFile() {
	window.location="$promsPath{cgi}/editDatabank.cgi?ID=$databankID&ACT=update";
}
var fileOpts=['mascotFile','sharedFile','localFile','internetFile'];
function updateFileOptions(selOption) {
	var myForm=document.dbForm;
	for (var i=0; i<fileOpts.length; i++) {
		if (!myForm[fileOpts[i]]) {continue;} // sharedFile may be missing
		myForm[fileOpts[i]].disabled=(fileOpts[i]==selOption)? false : true;
	}
	updateRulesButton();
}
function checkForm(myForm) {
	if (myForm.name.value=='') {
		alert('Type a name for this databank.');
		return false;
	}
	if (myForm.hasDecoy && myForm.hasDecoy.value=='yes' && !myForm.decoyTag.value) {
		alert('Provide a decoy tag for internal decoy entries.');
		return false;
	}
	if (myForm.species.value=='new' && !myForm.new_species.value) {
		alert('Provide a name for the new species used.');
		return false;
	}
	if (action=='add') {
		if (myForm.dbType.value==0) {
			alert('Select a databank type.');
			return false;
		}
		var okFile=true;
		if (!myForm.fileOption[0].checked && !myForm.fileOption[1].checked && !myForm.fileOption[2].checked && !myForm.fileOption[3].checked) {okFile=false;}
		else {
			if (myForm.fileOption[0].checked && !myForm.mascotFile.value) {okFile=false;}
			else if (myForm.fileOption[1].checked && !myForm.sharedFile.value) {okFile=false;}
			else if (myForm.fileOption[2].checked && !myForm.localFile.value) {okFile=false;}
			else if (myForm.fileOption[3].checked && !myForm.internetFile.value) {okFile=false;}
		}
		if (!okFile) {
			alert('Select a FASTA-formatted file for this databank.');
			return false;
		}
		if (myForm.fileOption[2].checked) { // local file selected
			popupWindow();
		}
	}
	myForm.save.style.visibility='hidden'; // prevent user from sending form twice!
	return true;
}
function popupWindow() {
	popupWin=window.open('','uploadWindow','width=500,height=200');
	var doc=popupWin.document;
	doc.open();
	doc.writeln('<HTML><HEAD>')
	doc.writeln('<TITLE>File Upload</TITLE>')
	doc.writeln('<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">');
	doc.writeln('</HEAD><BODY><CENTER>');
	doc.writeln('<BR><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">');
	doc.writeln('<FONT class="title2">Uploading sequence file...<BR>This could take a few minutes...</FONT>');
	doc.writeln('</CENTER></BODY></HTML>');
	doc.close;
}
function closeWindow() {
	if (popupWin != null) {
		if (!popupWin.closed) {
			popupWin.close();
		}
	}
}
// AJAX --->
function ajaxTestParseRules() {
	var defDiv=document.getElementById('rulesDIV');
	//defDiv.style.display='';
	defDiv.innerHTML="<BR><BR><FONT class=\\"title3\\">Fetching data. Please wait...</FONT><BR><IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR>";

	var paramStrg="ACT=rules";
	if ($databankID) {
		paramStrg+="&ID=$databankID";
	}
	else { // add : mascot or server file only!
		var dbTypeID=document.getElementById('dbType').value;
		paramStrg+="&TYPE_ID="+dbTypeID;
		if (document.dbForm.fileOption[0].checked) { //mascot file
			paramStrg+="&SRC=mascot&FILE="+document.getElementById('mascotFile').value;
		}
		else { // server
			paramStrg+="&SRC=server&FILE="+document.getElementById('sharedFile').value;
		}
	}

	//Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if (!XHR) {return false;}
	XHR.open("POST","$promsPath{cgi}/editDatabank.cgi",true);
	//Send the proper header information along with the request
	XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
	XHR.setRequestHeader("Content-length", paramStrg.length);
	XHR.setRequestHeader("Connection", "close");
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			defDiv.innerHTML=XHR.responseText;
			defDiv.scrollIntoView();
		}
	}
	XHR.send(paramStrg);
}
var XHR=null;
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
<BODY background="$promsPath{images}/bgProMS.gif" onunload="closeWindow()">
<CENTER>
<BR>
<FONT class="title">$pageTitle</FONT>
<BR><BR>
|;
####>FORM<####
print '<FORM name="dbForm" method="post" onsubmit="return checkForm(this);" action="./editDatabank.cgi" enctype="multipart/form-data">';
print hidden(-name=>'ACT',-default=>"$action"),"\n";
print hidden(-name=>'ID',-default=>"$databankID"),"\n" if $action eq 'edit';

####>TABLE
my ($light,$dark)=&promsConfig::getRowColors;
print qq
|<TABLE border=0>
	<TR><TD bgcolor=$dark><TABLE border=0 cellpadding=2>
|;
for my $i (0..$#nameText) {
	print "<TR>\n\t<TH align=right valign=top bgcolor=$dark width=150>$nameText[$i] :</TH>\n";
	print "\t<TD nowrap bgcolor=$light>$valueText[$i]</TD>\n</TR>\n";
}
my $resetString=($action eq 'add')? '  Clear  ' : ' Cancel Changes ';
print qq
|<TR><TD colspan=2 align=center>
<!-- SUBMIT button -->
<INPUT type="submit" name="save" value=" Save " />
<!-- CLEAR button -->
&nbsp &nbsp &nbsp<INPUT type="reset" value="$resetString" />
<!-- CANCEL button -->
&nbsp &nbsp &nbsp<INPUT type="button" value="Cancel" onclick="window.location='./selectDatabank.cgi'">
</TD></TR>
</TABLE>
</TD></TR>
</TABLE>
</FORM>
<DIV id="rulesDIV"></DIV><BR><BR>
</CENTER>
</BODY>
</HTML>
|;

############################<< SUBROUTINES >>###############################

###################################
####<Processing submitted form>####
###################################
sub processForm {
	print header(-'content-encoding'=>'no',-charset=>'utf-8'),"\n";
	warningsToBrowser(1); # DEBUG
	print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<CENTER>
<BR><BR><IMG src='$promsPath{images}/engrenage.gif'>
</CENTER>
|;

	####<Processing form parameters>####
	my $name=param('name'); $name=&promsMod::resize($name,50);
	my $versionName=param('versionName'); $versionName=&promsMod::resize($versionName,20);
	my $versionDate=param('versionDate');
	my $des=param('des'); $des=&promsMod::resize($des,100);
	my $identType=param('identifierType');
	my $decoyTag=param('decoyTag');
	my $species=(param('species') eq 'new')? param('new_species') : param('species');
	my $isContaminant=(param('contaminant'))? param('contaminant') : 0;
	my $fileOption=param('fileOption');
	my $comments=param('comments'); $comments=&promsMod::resize($comments,250);
	my @colName=('NAME','VERSION_NAME','VERSION_DATE','DES','IDENTIFIER_TYPE','DECOY_TAG','ORGANISM','COMMENTS');
	my @colValue=($dbh->quote($name),$dbh->quote($versionName),$dbh->quote($versionDate),$dbh->quote($des),$dbh->quote($identType),$dbh->quote($decoyTag),$dbh->quote($species),$dbh->quote($comments));

	####<Add>####
	if ($action eq 'add') {

		###<Computing new DB ID>###
		($databankID)=$dbh->selectrow_array("SELECT MAX(ID_DATABANK) FROM DATABANK");
		$databankID++;

		###<Copying file to dbPath/db_ID>###
		my ($mascotServer,$dbankName,$fileName,$fullDbankfile,$fileString);

		#<Mascot file
		if ($fileOption eq 'mascotFile') {
			($mascotServer,$dbankName,$fullDbankfile)=split('##',param('mascotFile'));
			$fileName=(split(/\//,$fullDbankfile))[-1];
			$fileString="$mascotServer:$dbankName:$fileName";

			if ($mascotServers{$mascotServer}[2]) { # Mascot is accessible from myProMS server
				$fullDbankfile=~s/.+\/sequence/$mascotServers{$mascotServer}[2]\/sequence/; # in case NFS mount of mascot root dir
			}
		}
		else { #  server,local,internet
			my $dbDir="$promsPath{banks}/db_$databankID";
			mkdir $promsPath{'banks'} unless -e $promsPath{'banks'};
			mkdir $dbDir unless -e $dbDir;

			#<server file
			if ($fileOption eq 'sharedFile') {
				my $sharedFile=param('sharedFile');
				$fileName=(split(/[\\\/]/,$sharedFile))[-1];
				$fullDbankfile="$dbDir/$fileName";
				move("$promsPath{shared}/$sharedFile",$fullDbankfile);
			}

			#<local file
			elsif ($fileOption eq 'localFile') {
				print "<FONT class=\"title2\">";
				$fileName=(split(/[\\\/]/,param('localFile')))[-1];
				$fullDbankfile="$dbDir/$fileName";
				my $tmpfile = tmpFileName(upload('localFile')); # name of temp file being uploaded
				move($tmpfile,$fullDbankfile);
# 			print "<FONT class=\"title2\">Writing $fileName to ProMS Server... ";
# 			my $uploadFile=param('localFile');	# also works  =upload('localFile');
# 			open (DB, ">$fullDbankfile");
# 			while (<$uploadFile>) {
# 				print DB $_;
# 			}
# 			close DB;
# 			print "Done.<BR>\n";
			}

			#<internet file
			elsif ($fileOption eq 'internetFile') {
				my $internetFile=param('internetFile');
				$fileName=(split(/\//,$internetFile))[-1];
				$fullDbankfile="$dbDir/$fileName";
				my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
				$agent->timeout(360);
				if (!$promsPath{'http_proxy'}) {$agent->env_proxy;}
				elsif (lc($promsPath{'http_proxy'}) eq 'no') {$agent->no_proxy;}
				else {$agent->proxy(['http','ftp'], $promsPath{'http_proxy'});}
				my $prevPc=0;
				my $thPc=10;
				print "<BR><BR><FONT class=\"title2\">Downloading file:</FONT>\n<BR><B>0%";
				open (DB,">$fullDbankfile");
				my $expectedLength;
				my $bytesReceived=0;
				my $response=$agent->request(HTTP::Request->new(GET => $internetFile),
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
									print DB $chunk;
								});

				close DB;
				if ($thPc <= 100) { # download was incomplete (complete= 110%!)
					unlink $fullDbankfile;
					$dbh->rollback;
					$dbh->disconnect;
					print "FAILED !<BR>\nError return: '",$response->status_line,"'";
					if ($internetFile=~/^ftp/) {
						$internetFile=&promsMod::resize($internetFile,200);
						$internetFile=~s/^ftp/<FONT color=#DD0000>http<\/FONT>/;
						print "<BR>Try again using http protocole instead ($internetFile)";
					}
					print qq
|</B><BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="button" value=" Try again " onclick="history.back()"/>&nbsp;&nbsp;<INPUT type="button" value=" Cancel " onclick="window.location='./selectDatabank.cgi'" />
</BODY>
</HTML>
|;
					exit;
				}
				print " Done</B>\n<BR><FONT class=\"title2\">Download is complete.<BR>\n";
			}

			###<Inflating file>###
			if ($fileName=~/\.(gz|zip)\Z/) {
				my $fullCompDbankfile=$fullDbankfile;
				my $compFileName=$fileName;
				$fullDbankfile=~s/\.(gz|zip)\Z//;
				$fileName=~s/\.(gz|zip)\Z//;
				print "<BR><FONT class=\"title2\">Inflating $compFileName...";
				if ($compFileName=~/\.gz\Z/) {
					#my @fileInfo=`gunzip -l $fullDbankfile`;
					#chop($fileInfo[1]);
					#($fileName)=($fileInfo[1]=~/([^\/\s]+)\Z/);
					#system "gunzip -q $fullDbankfile&";
					gunzip $fullCompDbankfile => $fullDbankfile or die "gzip error : $!";
				}
				elsif ($compFileName=~/\.zip\Z/) {
					#my @fileInfo=`unzip -l $fullDbankfile`;
					#chop($fileInfo[3]);
					#($fileName)=($fileInfo[3]=~/([^\/\s]+)\Z/);
					#system "unzip -qo -d $dbDir $fullDbankfile && rm $fullDbankfile&";
					unzip $fullCompDbankfile => $fullDbankfile or die "unzip failed: $UnzipError\n";
				}
				#while (-e $fullDbankfile) {
				#	print '.';
				#	sleep 3;
				#}
				#$fullDbankfile="$dbDir/$fileName"; # inflated file!!
				unlink $fullCompDbankfile;
				print "Done.</FONT><BR>\n";
			}
			$fileString=$fileName;
		}


		###<Checking file>###
		print "<BR><FONT class=\"title2\">Processing $fileName... ";
		if ($fileOption eq 'mascotFile' && !$mascotServers{$mascotServer}[2]) { # Distant Mascot server
			my $params={ACT=>'dbFile',DB=>$dbankName};
			my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
			$agent->timeout(360);
			if ($mascotServers{$mascotServer}[1]) { # proxy settings
				if ($mascotServers{$mascotServer}[1] eq 'no') {$agent->no_proxy($mascotServers{$mascotServer}[0]);}
				else {$agent->proxy('http', $mascotServers{$mascotServer}[1]);}
			}
			else {$agent->env_proxy;}
			my $response = $agent->post("$mascotServers{$mascotServer}[0]/cgi/myproms4databanks.pl",$params);
			while (my $wait = $response->header('Retry-After')) {
				sleep $wait;
				$response = $agent->get($response->base);
			}
			my @resultLines;
			if ($response->is_success) {
				@resultLines = split("\n",$response->content);
			}
			foreach (@resultLines) {
				next if /^#/;
				$totNumEntry=$1 if /^\S+\s+(\d+)/;
			}
		}
		else {$totNumEntry=`grep -c '>' $fullDbankfile`; chomp $totNumEntry;}

		if ($totNumEntry > 0) {
			print "Done ($totNumEntry entries).</FONT><BR>\n";

			###<Inserting data into tables DATABANK>###
			my $dbTypeID=param('dbType');
			push @colName,('FASTA_FILE','ID_DBTYPE');
			push @colValue,($dbh->quote($fileString),$dbTypeID);
			my $colNameString=join(",",@colName);
			my $colValueString=join(",",@colValue);
			$dbh->do("INSERT INTO DATABANK (ID_DATABANK,IS_CRAP,$colNameString,NUM_ENTRY,USE_STATUS) VALUES ($databankID,$isContaminant,$colValueString,$totNumEntry,'yes')") || die $dbh->errstr;
		}
		else {
			print "<BR><FONT color=#DD0000>ERROR: No entries found.</FONT></FONT>\n";
			#system "rm -r $promsPath{banks}/db_$databankID"; # not for mascot files
			#remove_tree("$promsPath{banks}/db_$databankID");
#rmtree("$promsPath{banks}/db_$databankID");
			print "&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" value=\" Try again \" onclick=\"history.back()\" />&nbsp;&nbsp;<INPUT type=\"button\" value=\" Cancel \" onclick=\"window.location='./selectDatabank.cgi'\" />\n";
			print "</BODY>\n</HTML>\n";
			exit;
		}
		sleep 2;
	}
	else { # edit
		my $updateString;
		for my $c (0..$#colName) { # looping through all conlumns
			$updateString.="$colName[$c]=$colValue[$c]";
			$updateString.=',' if $c<$#colName;
		}
		$dbh->do("UPDATE DATABANK SET $updateString WHERE ID_DATABANK=$databankID") || die $dbh->errstr;
	}

	$dbh->commit;
	$dbh->disconnect;

	print qq
|<SCRIPT LANGUAGE="JavaScript">
window.location ="./selectDatabank.cgi";
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}


#############################
####<Deleting a databank>####
#############################
sub deleteDB {

	my ($numAna)=$dbh->selectrow_array("SELECT count(*) FROM ANALYSIS_DATABANK WHERE ID_DATABANK=$databankID");
	if ($numAna) {
		####<DB was used => Setting use_status to 'no'>####
		$dbh->do("UPDATE DATABANK SET USE_STATUS='no' WHERE ID_DATABANK=$databankID") || die $dbh->errstr;
	}
	else {
		####<DB never used => Deleting entry>####
		$dbh->do("DELETE FROM DATABANK WHERE ID_DATABANK=$databankID");
	}

	####<Deleting all corresponding files>####
	#system "rm -r $promsPath{banks}/db_$databankID" if -e "$promsPath{banks}/db_$databankID"; # does not exist for Mascot DB
	#remove_tree("$promsPath{banks}/db_$databankID") if -e "$promsPath{banks}/db_$databankID"; # does not exist for Mascot DB
	rmtree("$promsPath{banks}/db_$databankID") if -e "$promsPath{banks}/db_$databankID"; # does not exist for Mascot DB

	$dbh->commit;
	$dbh->disconnect;

	print redirect("./selectDatabank.cgi");

	exit;
}

sub ajaxTestParseRules {

	my ($dbName,$dbTypeName,$dbFile,$numEntries,$parseRules,$fileSource);
	if ($databankID) { # edit
		($dbName,$dbTypeName,$dbFile,$numEntries,$parseRules)=$dbh->selectrow_array("SELECT D.NAME,DT.NAME,FASTA_FILE,NUM_ENTRY,PARSE_RULES FROM DATABANK D,DATABANK_TYPE DT WHERE D.ID_DBTYPE=DT.ID_DBTYPE AND ID_DATABANK=$databankID");
		$fileSource=($dbFile=~/:.+:/)? 'mascot' : 'local'; # mascot edit
	}
	else { # add
		my $dbTypeID=param('TYPE_ID');
		$fileSource=param('SRC'); # mascot or server
		$dbFile=param('FILE');
		($dbTypeName,$parseRules)=$dbh->selectrow_array("SELECT NAME,PARSE_RULES FROM DATABANK_TYPE WHERE ID_DBTYPE=$dbTypeID");
		if ($fileSource eq 'server') {$numEntries=`grep -c '^>' $promsPath{shared}/$dbFile`; chomp $numEntries;} # uncompressed fasta
		else {$numEntries=0;}
	}

	$dbh->disconnect;

	print header(-type=>'text/plain',-'content-encoding'=>'no',-charset=>'utf-8'),"\n";
	warningsToBrowser(1); # DEBUG
	print qq
|<FONT class="title">Testing databank rules <FONT color='#DD0000'>$dbTypeName</FONT></FONT>
<BR>
|;
	my @entryList;
	##<Mascot file
	if ($fileSource eq 'mascot') { # edit or add
		my ($mascotServer,$dbankDir,$dbFileName)=($databankID)? split(/:/,$dbFile) : split(/##/,$dbFile);
		my %mascotServers=&promsConfig::getMascotServers;
		my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
		$agent->timeout(360);
		#<Proxy
		if ($mascotServers{$mascotServer}[1]) {
			if ($mascotServers{$mascotServer}[1] eq 'no') {$agent->no_proxy($mascotServers{$mascotServer}[0]);}
			else {$agent->proxy('http', $mascotServers{$mascotServer}[1]);}
		}
		else {$agent->env_proxy;}
		my $response=$agent->post("$mascotServers{$mascotServer}[0]/cgi/myproms4databanks.pl",
									['ACT'=>'parse',
									'DB'=>$dbankDir,
									#'parseRules'=>uri_escape($parseRules), # URL encoding unsafe parameters
									'parseRules'=>$parseRules,
									'numEntries'=>$numEntries
									]
								   );
		if ($response->is_success) {
			if ($response->content=~/^#Error/) {
				(my $errorText=$response->content)=~s/^#Error: //;
				print "<FONT class=\"title2\"><FONT color='#DD0000'>**ERROR**</FONT>: Unexpected response from $mascotServer: \"$errorText\"</FONT><BR>\n";
			}
			else {
				@entryList=split(/\n/,$response->content);
				my $trueDbFileName=(split(/\/|\\/,shift(@entryList)))[-1]; # 1 line is db file, keep file name
				print "<FONT class=\"title2\">(Databank: $dbankDir, File: $trueDbFileName on $mascotServer)</FONT><BR><BR>\n";
			}
		}
		else { # Error
			print "<FONT class=\"title2\"><FONT color='#DD0000'>**ERROR**</FONT>: Bad anwser from $mascotServer!</FONT><BR>$!\n";
		}
	}

	##<Local/server file
	else {
		my ($idRule,$desRule,$orgRule)=&getParseRules($parseRules);
		my $trueDbFile=($fileSource eq 'local')? "$promsPath{banks}/db_$databankID/$dbFile" : "$promsPath{shared}/$dbFile";

		my $maxRange=(!$numEntries)? 10 : ($numEntries <= 1000)? $numEntries : 1000;
		my %entryPos;
		if ($maxRange==10) {
			foreach my $i (1..10) {$entryPos{$i}=1;} # pick first 10 positions
		}
		else {
			foreach my $i (1..10) {$entryPos{int(rand($maxRange)+0.5)}=1;} # pick 10 random positions
		}
		open (DB,$trueDbFile) || die "can't open $trueDbFile\n";
		my $curPos=0;
		while (<DB>) {
			last if scalar keys %entryPos==0;
			if (/^>/) {
				$curPos++;
				next unless $entryPos{$curPos};
				(my $entry=$_)=~s/^>\s*//;
				$entry =~s /.+//; # take 1 entry if multi-entry line (NCBI)
				$entry =~s /\s+\Z//; # chomp not always enough
				my ($identifier)=($entry=~/$idRule/);
				my $des; ($des)=($entry=~/$desRule/) if $desRule; # optional
				my $org; ($org)=($entry=~/$orgRule/) if $orgRule; # optional
				if ($org) {
					($des)=($entry=~/^\S+\s+(.+)$orgRule/) unless $desRule;
					$org=~s/\s+\Z//; # removing trailing spaces
				}
				else {
					($des)=($entry=~/^\S+\s+(.+)/) unless $desRule;
				}
				$des='' unless $des;
				$des=~s/\s+\Z//; # removing trailing spaces
				$org='' unless $org;
				push @entryList,"$identifier ##DES=$des ##ORG=$org ##HEAD=$entry\n";
				delete $entryPos{$curPos};
			}
		}
		close DB;
		if ($fileSource eq 'local') {
			print "<FONT class=\"title2\">(Databank: $dbName, Local file: $dbFile)</FONT><BR><BR>\n";
		}
		else { # shared folder
			print "<FONT class=\"title2\">(Shared directory file: $dbFile)</FONT><BR><BR>\n";
		}
	}

	if (scalar @entryList) {
		my ($light,$dark)=&promsConfig::getRowColors;
		print qq
|<TABLE cellspacing=0>
<TR bgcolor="$dark">
	<TH class="rbBorder" nowrap>Identifier</TH>
	<TH class="rbBorder" width=450>Description</TH>
	<TH class="rbBorder" nowrap>Species</TH>
	<TH class="bBorder" width=650>Full header</TH>
</TR>
|;
		my $color=$light;
		foreach my $entryStrg (@entryList) {
			my ($identifier,$des,$org,$fullEntry)=($entryStrg=~/^(\S+) ##DES=(.*) ##ORG=(.*) ##HEAD=(.+)/);
			$des='*Not found*' unless $des;
			my $orgStrg=($org)? "<FONT class=\"org\">$org</FONT>" : '*Not found*';
			print "<TR bgcolor=$color><TH align=left valign=top>&nbsp;$identifier&nbsp;</TH><TD valign=top>$des</TD><TD valign=top>&nbsp;$orgStrg&nbsp;</TD><TD valign=top>$fullEntry&nbsp;</TD></TR>\n";
			$color=($color eq $light)? $dark : $light;
		}
		print "</TABLE>\n";
	}

	exit;
}

sub updateMascotDbFile {

	my ($dbFile)=$dbh->selectrow_array("SELECT FASTA_FILE FROM DATABANK WHERE ID_DATABANK=$databankID");
	$dbh->disconnect;

	my ($mascotServer,$dbankDir,$dbFileName)=split(/:/,$dbFile);

	print header(-'content-encoding'=>'no',-charset=>'utf-8'),"\n";
	warningsToBrowser(1); # DEBUG
	print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<CENTER>
<BR><BR><IMG id="waitImg" src='$promsPath{images}/engrenage.gif'>
<BR><BR><FONT class="title2">Checking for update of databank <FONT color='#DD0000'>$dbankDir</FONT> on $mascotServer:</FONT>
<BR><BR>
<FONT class="title3">Scanning databank entries...|;

	my ($updated,$newDbFileName,$numEntry,$errorText)=&promsMod::updateMascotDB($databankID,1);
	print " Done</FONT><BR><BR>\n";

	if ($errorText) {
		print "<FONT class=\"title2\"><FONT color='#DD0000'>**ERROR**</FONT>: $errorText</FONT>\n";
	}
	elsif ($updated) {
		print "<FONT class=\"title2\">Update found: New file is '$newDbFileName' ($numEntry entries).</FONT><BR>\n";
	}
	else {
		print "<FONT class=\"title2\">No update found: File '$dbFileName' is still in use.</FONT><BR>\n";
	}

	print qq
|<BR><BR><INPUT class="title3" type="button" value="  OK  " onclick="window.location='./editDatabank.cgi?ACT=edit&ID=$databankID'"/>
</CENTER>
<SCRIPT LANGUAGE="JavaScript">
document.getElementById('waitImg').style.visibility='hidden';
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

sub getParseRules {
	my $parseRules=$_[0];
	my @rules=split(',:,',$parseRules);
	my ($idRule)=($rules[0]=~/ID=(.+)/);
	my ($desRule,$orgRule);
	if ($rules[1]) {
		if ($rules[1]=~/DES=/) {($desRule=$rules[1])=~s/DES=//;}
		else {($orgRule=$rules[1])=~s/ORG=//;}
	}
	if ($rules[2]) {($orgRule=$rules[2])=~s/ORG=//;}
	return ($idRule,$desRule,$orgRule);
}

####>Revision history<####
# 2.4.2 [Fix] bug in shared file path resolution (PP 11/10/18)
# 2.4.1 [Fix] Javascript bugs when no shared directory declared (PP 13/04/18)
# 2.4.0 Uses &promsMod::browseDirectory_getFiles for shared folder files (PP 04/12/17)
# 2.3.7 Bug fix in edit mode to handle contaminant field (PP 14/11/16)
# 2.3.6 Add contaminant checkbox for MaxQuant (GA 21/10/16)
# 2.3.5 Uses Perl modules to uncompress files (PP 20/10/16)
# 2.3.4 Moved &getIdentifierTypes from promsConfig.pm to promsMod.pm (PP 21/03/14)
# 2.3.3 Uses rmtree instead of remove_tree (PP 10/03/14)
# 2.3.2 Change syntax $x=$y // $z to old ()? one (PP 10/03/14)
# 2.3.1 Minor change in text (PP 20/11/13)
# 2.3.0 Moved Mascot databank update to &promsMod::updateMascotDB (PP 07/11/13)
# 2.2.9 Added parse rules test, databank file update for mascot & new file source (PP 04/11/13)
# 2.2.8 Added FTP to proxy declaration (PP 02/07/13)
# 2.2.7 Better proxy declaration for LWP (PP 02/07/13)
# 2.2.6 Uses LWP instead of wget (PP 14/05/13)
# 2.2.5 Compatible with multi-databank analyses (PP 18/12/12)
# 2.2.4 Added check on DB dir before delete (PP 05/07/12)
# 2.2.3 Minor bug correction (PP 26/03/12)
