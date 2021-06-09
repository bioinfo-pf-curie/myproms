#!/usr/local/bin/perl -w

################################################################################
# managemzXMLFiles.cgi      1.2.0                                              #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Script that allows to import (1 by 1) or delete mzXML files in a project     #
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

use strict;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use POSIX qw(strftime :sys_wait_h);  # Core module
use promsConfig;
use promsMod;
use File::Copy qw(copy move);
use File::Spec::Functions qw(splitpath); # Core module

# print header; warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
mkdir "$promsPath{tmp}/upload" unless -e "$promsPath{tmp}/upload";
$CGITempFile::TMPDIRECTORY="$promsPath{tmp}/upload";

###############################
####>Recovering parameters<####
###############################
my $projectID=param('id_project');
my $branchID=param('ID');
my $action=(param ('ACT'))? param ('ACT') : 'list';
my $mzXMLPath="$promsPath{tmp}/upload/project_$projectID";
my $maxUpFiles=10;

if ($action eq 'delete') {
	my @files = param('filebox');
	print header(-'content-encoding'=>'no'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<CENTER><FONT class="title">Deleting Files...</FONT></CENTER>
<BR><BR>
<FONT class="title3">
|;
	foreach my $file (@files) {
		unlink "$mzXMLPath/$file";
		print "&nbsp;$file... deleted.<BR>\n";
	}
	print "<BR><BR>All files have been deleted.</FONT>\n";

	sleep 5;
	print qq
|<SCRIPT type="text/javascript">
	window.location="./managemzXMLFiles.cgi?ID=$branchID&id_project=$projectID";
</SCRIPT>
</BODY>
</HTML>
|;
}

elsif ($action eq 'import') {
	####>Starting HTML<####
	print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>MZXML File Management</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY>
<BR>
<FONT class="title3">Retrieving files from your computer...</FONT><BR><BR>
|;
	my $form = new CGI;
	mkdir $mzXMLPath unless -e $mzXMLPath;
	my @upload_fh = $form->upload("mzXMLfiles");
	foreach my $fileName (@upload_fh) {
		my $tmpfile=tmpFileName($fileName);
		#print "tmpfile=$tmpfile<BR>\n$mzXMLPath/$fileName<BR>\n";
		print "<B>&nbsp;-$fileName...";
		move($tmpfile,"$mzXMLPath/$fileName");
		print " Done.</B><BR>\n";
	}
	sleep 3;
	print qq
|<SCRIPT type="text/javascript">
window.location="$promsPath{cgi}/managemzXMLFiles.cgi?ID=$branchID&id_project=$projectID";
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}
elsif ($action eq 'UseSharedDirectory') {

	####>Starting HTML<####
	print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>MZXML File Management</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY>
<BR>
<FONT class="title3">Retrieving files from shared data directory...</FONT><BR><BR>
|;
	mkdir $mzXMLPath unless -e $mzXMLPath;
	foreach my $sharedFile (param('sharedDirFiles')) {
		my (undef,$path,$fileName)=splitpath($sharedFile);
		print "<B>&nbsp;-$fileName...";
		my $newFile="$mzXMLPath/$fileName";
		&copyAndPrint("$promsPath{shared}/$sharedFile",$newFile);
		if ($fileName=~/\.zip$/) {
			print " Extracting files...<BR>\n";
			&promsMod::unzipArchive($newFile,$mzXMLPath,{move2Top=>1,mode=>'verbose',txtBefore=>'&nbsp;&nbsp;-',txtAfter=>"<BR>\n"});
			print "&nbsp;Done.</B><BR>\n";
		}
		elsif ($fileName =~ /\.gz$/) {
			print " Extracting files...";
			if ($newFile =~ /\.tar\.gz\Z/) {
				system "tar -zxf $newFile -C $mzXMLPath";
			}
			elsif ($newFile =~ /\.gz\Z/) {
				system "gunzip -c $newFile > $mzXMLPath";
			}
			print " Done.</B><BR>\n";
		}
		else {print " Done.</B><BR>\n";}
	}
	sleep 3;
	print qq
|<SCRIPT type="text/javascript">
window.location="$promsPath{cgi}/managemzXMLFiles.cgi?ID=$branchID&id_project=$projectID";
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

if ($action ne 'delete') {

	################
	####> HTML <####
	################
	my ($lightColor,$darkColor)=&promsConfig::getRowColors;

	print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<HEAD>
<TITLE>MZXML File Management</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TD.center {text-align:center} //{font-weight:bold;}
</STYLE>
<SCRIPT type="text/javascript">
var fileSource;
function cancelAction() {
	//top.promsFrame.optionFrame.selectOption();
	window.location="./processAnalyses.cgi?ID=$branchID";
}
function testCheckbox(myForm) {
	var selected=0;
	var anaBox=myForm.filebox;
	if (anaBox.length) { // more than 1 checkboxes
		for (i=0; i < anaBox.length; i++){
			if (anaBox[i].checked==true) {selected=1; break;}
		}
	}
	else if (anaBox.checked==true){selected=1;}
	return selected;
}
function checkFile(myFile) {
	if (myFile.files.length == 0) {
		alert('ERROR: "Please browse one or more files.');
	}else{
		for (var i = 0; i < myFile.files.length; i++) {
			if((!myFile.files[i].name.match('mzXML'))&&(!myFile.files[i].name.match('mzML'))){
				alert('ERROR: File(s) without mzXML or mzML extension detected: '+myFile.files[i].name+'.');
				myFile.value='';
				return;
			}
		}
		document.getElementById('submitButton').disabled=false;
	}
}
function checkAutoExtend(boxIdx,chkStatus) {
	var fileCheckBox=document.fileCleanForm.filebox;
	if (!document.fileCleanForm.autoExtend.checked \|\| !fileCheckBox.length) {return;} // no extension or only 1 file
	boxIdx++;
	for (let i=boxIdx; i<fileCheckBox.length; i++) {
		fileCheckBox[i].checked=chkStatus;
	}
}
function checkCleanForm(myForm) {
	var myForm=document.fileImportForm;
	myForm.ACT.value=fileSource;
	if (testCheckbox(myForm)==0) {
		alert('Error: No files selected !');
		return false;
	}
	if (confirm('Delete selected file(s) from server ?')) {return true;}
	return false;
}
function selectSource(source) {
	fileSource=source;
	var myForm=document.fileImportForm;
	document.getElementById('submitButton').disabled=false;

	/* Disable all file sources */
	if (document.getElementById('sharedDirDIV')) {document.getElementById('sharedDirDIV').style.display='none';} // sharedDIr unabled
	myForm.mzXMLfiles.disabled=true;

	/* Activate selected source */
	if (fileSource=='UseSharedDirectory') { // Shared directory
		document.getElementById('sharedDirDIV').style.display='';
		myForm.ACT.value=fileSource;
	}
	else if (fileSource=='UseUploadedFiles') { // Upload multiple files
		myForm.mzXMLfiles.disabled=false;
	}
}
|;
	&promsMod::browseDirectory_JavaScript() if $promsPath{'shared'};
	print qq
|</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" >
<CENTER>
<FONT class="title">Manage mzXML files</FONT><BR><BR>
|;

	##########################################
	####>Clean old files for ALL projects<####
	##########################################
	if ($action eq 'list') {
		print qq
|<DIV id="waitDIV">
<BR><BR><BR><FONT class="title3"><SPAN id="waitSPAN">Updating mzXML file list. Please wait...</SPAN></FONT><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">
<BR><BR><BR><BR>
</DIV>
|;
		opendir(UPLOAD,"$promsPath{tmp}/upload");
		while (defined (my $child = readdir(UPLOAD))) {
			next if ($child eq '.' || $child eq '..' || -l "$promsPath{tmp}/upload/$child"); # must be a valid file or directory
			my $remaining=&promsMod::cleanDirectory("$promsPath{tmp}/upload/$child",'3M'); # 3 months
			rmdir "$promsPath{tmp}/upload/$child" if (-d "$promsPath{tmp}/upload/$child" && !$remaining); # empty dir
			print qq
|<SCRIPT type="text/javascript">
document.getElementById('waitSPAN').innerHTML+='.';
</SCRIPT>
|;
		}
		close UPLOAD;
		print qq
|<SCRIPT type="text/javascript">
document.getElementById('waitSPAN').innerHTML+=' Done.';
setTimeout(function(){document.getElementById('waitDIV').style.display='none';},1000);
</SCRIPT>
|;
	}

	###################################################
	####>Get already mzXML files imported (if any)<####
	###################################################
	my %filesList;
	if (-e $mzXMLPath){
		opendir (DIR, $mzXMLPath) || print "ERROR: Unable to read '$mzXMLPath' !<BR>\n";
		while (defined (my $currentmzXMLFile = readdir (DIR))) {
			next unless $currentmzXMLFile =~ /.+\.mzX*ML\Z/;
			$filesList{$currentmzXMLFile}=(stat "$mzXMLPath/$currentmzXMLFile")[7]; # size
		}
		closedir DIR;
	}

	####> Option to add a new one
	print qq
|<FORM name="fileImportForm" method="post" action="managemzXMLFiles.cgi" enctype="multipart/form-data">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="ACT" value="import">
<INPUT type="hidden" name="ID" value="$branchID">

<TABLE bgcolor=$darkColor cellpadding=4 style="min-width:500px">
	<TR>
		<TH colspan="3"><FONT class="title2">&nbsp;Upload file(s) : select source</TH>
	</TR>
	<TR bgcolor=$darkColor></TD>
		<TH valign="top"><INPUT type="radio" name="selSource" value="UseUploadedFiles" onclick="selectSource(this.value);"></TH>
		<TH align="right" valign="top" nowrap>Upload multiple files :</TH>
		<TD bgcolor=$lightColor><INPUT type="file" id="mzXMLfiles" name="mzXMLfiles" multiple size=120 onchange="checkFile(this);"></TD>
	</TR>
|;

	if ($promsPath{'shared'}) {
			print qq
|	<TR bgcolor=$darkColor>
		<TH valign="top"><INPUT type="radio" name="selSource" value="UseSharedDirectory" onclick="selectSource(this.value);"></TH>
		<TH align="right" valign="top" nowrap>Shared data directory :</B></TH>
		<TD bgcolor=$lightColor><DIV id="sharedDirDIV" style="max-height:300px;overflow:auto;display:none">
|;
			&promsMod::browseDirectory_getFiles($promsPath{'shared'},'sharedDirFiles',{fileMatch=>qr/\.mzX*ML$/i});
			print qq
|			</DIV><BR>
			</TD>
|;
	}

	print qq
|<TR bgcolor=$darkColor>
<TH colspan="3">
    <INPUT type="submit" id="submitButton" value=" Proceed " disabled>&nbsp;&nbsp;&nbsp;
    <INPUT type="button" value=" Cancel " onclick="cancelAction();">
</TH></TR>

</TABLE><BR>
</FORM>
|;


	####> Already imported mzXML files
	print qq
|<FORM name="fileCleanForm" method="post" onsubmit="return checkCleanForm(this);" enctype="multipart/form-data">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="ACT" value="delete">
<INPUT type="hidden" name="ID" value="$branchID">

<TABLE bgcolor=$darkColor cellpadding=4 style="min-width:500px">
<TR><TH><FONT class="title2">&nbsp;Files already imported:&nbsp;</FONT></TH></TR>
<TR><TD><INPUT type="checkbox" name="autoExtend" value=1>Auto-extend selection</TD></TR>
|;

	my $bgColor=$lightColor;
	my $boxIdx=0;
	foreach my $dataFile (sort {lc($a) cmp lc($b)} keys %filesList) {
		my $fileSizeStrg=1*(sprintf "%.2f",$filesList{$dataFile} / (1024 ** 2))." Mo";
		print "<TR bgcolor=$bgColor><TH align=left nowrap><LABEL><INPUT type=\"checkbox\" name=\"filebox\" id=\"filebox\" value=\"$dataFile\" onchange=\"checkAutoExtend($boxIdx,this.checked);\">$dataFile</LABEL> - <I>$fileSizeStrg</I></TH></TR>\n";
		$boxIdx++;
		$bgColor=($bgColor eq $darkColor)? $lightColor : $darkColor;
	}
	my $deleteButtonStrg='';
	if (scalar keys %filesList) {
		$deleteButtonStrg='<INPUT type="submit" name="delete" value=" Delete ">&nbsp;&nbsp;&nbsp;';
	}
	else {
		print "<TR bgcolor=$lightColor><TH align=left nowrap>No files found.</TH></TR>\n";
	}
	print qq
|<TR><TH bgcolor=$darkColor>$deleteButtonStrg<INPUT type="button" value="Cancel" onclick="cancelAction();"></TH></TR>
</TABLE>
<BR><BR>
</CENTER>
</BODY>
</HTML>
|;
}


sub copyAndPrint {
	my ($sourceFile, $targetFile) = @_;
	my (undef, $path, $fileName) = splitpath($sourceFile);

	my $sourceFileSize = (stat $sourceFile)[7];
	my $sourceSizeGo = $sourceFileSize / (1024 ** 3);
	
	if ($sourceSizeGo < 0.1) {
		copy($sourceFile, $targetFile);
	}
	else { # fork and wait
		my $maxCopyTime = ($sourceSizeGo > 1)? $sourceSizeGo * 20 : 20;  # max 20min or 20min per Go
		my $targetFileSize = 0;
		my $loopNb = 0;  # To avoid infinite loop
	
		# Fork to copy files (makes it possible to keep the page active and avoid timeout for big files)
		my $childPid = fork;
		unless ($childPid) { # child here
			#>Disconnecting from server
			open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
			open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
			open STDERR, ">>/dev/null";
			copy($sourceFile, $targetFile);
			exit;
		}
		# print "<BR>&nbsp;&nbsp;&nbsp;&nbsp;- Storing $fileName.";
		sleep 1;  # To avoid sleeping 10 sec after in most cases
	
		# Start tracking target file and check if the copy is actually being done
		while ($loopNb < 6 * $maxCopyTime) {  # 10sec loops -> loopNb = 6 * nb of minutes passed
			my $res = waitpid($childPid, WNOHANG); # WNOHANG (from POSIX ":sys_wait_h"): parent doesn't hang during waitpid
			last if $res; # child has ended
	
			sleep 10;
			print ".";
			# print "<BR>" if ($loopNb > 30);  # End line every 5 minutes
			$loopNb++;
		}
		if ($loopNb >= 6 * $maxCopyTime) {
			print "Killing child process $childPid" if ($childPid);
			kill "SIGKILL", $childPid if ($childPid);
			die "Retrieving file is taking too long or process died silently before completion";
		} else {
			$targetFileSize = (-e "$targetFile") ? (stat $targetFile)[7] : 0;
			if ($targetFileSize == $sourceFileSize) {
				print("<BR>&nbsp;&nbsp;&nbsp;&nbsp;Done.");
			} else {
				die "Copy of file was not done properly. Exiting...";
			}
		}
	}
}

####> Revision history
# 1.2.0 [ENHANCEMENT] Copy no move from shared directory with fork to wait for large files & file size added (PP 08/06/21)
# 1.1.0 Uses &promsMod::cleanDirectory to delete all XML files older than 3 months & fixed JavaScript errors (PP 06/02/18)
# 1.0.11 Update to share directory (GA 20/11/17)
# 1.0.10 Set default CGI upload directory to $promsPath{tmp}/upload (PP 18/09/17)
# 1.0.9 Add mzML option import (GA 20/10/16)
# 1.0.8 Allow multiple file selection for mzXML (GA 15/10/14)
# 1.0.7 Uses mkdir instead of make_path (PP 10/03/14)
# 1.0.6 minor modif add File::Path qw(make_path) in header (GA 29/01/14)
# 1.0.5 system command removal (PP 08/11/13)
# 1.0.4 Minor changes (PP 06/09/13)
# 1.0.3 Minor add -> create upload directory (GA 07/12/12)
# 1.0.2 Minor bug correction -> remove printError call (GA 18/04/12)
# 1.0.1 Minor bug correction -> create mzXMLPath with system command (GA 11/04/12)
# 1.0.0 New script for mzXML files management (add or remove) for quantification purposes (GA 03/01/12)