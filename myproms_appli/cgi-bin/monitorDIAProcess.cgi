#!/usr/local/bin/perl -w

################################################################################
# monitorDIAProcess.cgi         1.0.4                                          #
# Authors: M. Le Picard (Institut Curie)                                       #
# Contact: myproms@curie.fr                                                    #
# monitor the libraries creation in myProMS                                    #
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

$|=1;       # buffer flush (allows instant printing)
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use strict;
use IO::Handle;
use promsConfig;
use File::Copy;
use File::Basename;
use XML::Simple;

my %promsPath=&promsConfig::getServerInfo;
my $dbh=&promsConfig::dbConnect;
my $userID=$ENV{'REMOTE_USER'};
my $workDir="$promsPath{data}/tmp/Swath";
my $action=param('ACT'); # library || import

print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);

if (param('AJAX')){
    if (param('AJAX') eq 'update'){&ajaxUpdateStatus;}
    if (param('AJAX') eq 'delete'){&ajaxDelete;}
    exit;
}

my $title=($action eq 'library')? "Monitoring libraries creation" : "Monitoring import of DIA data";

print qq |
<HTML>
<HEAD>
<TITLE>$title</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
function displayErrorDiv(jobDir,libID){
    var displayDiv=document.getElementById('error:'+jobDir+':'+libID);
    if (displayDiv.style.display==''){
        displayDiv.style.display='none';
    }
    else {
        displayDiv.style.display='';
    }
}

// -->AJAX

function deleteJob(jobDir,libID){
    autoUpdate=0;
    XHR=getXMLHTTP();
    XHR.open("GET","./monitorDIAProcess.cgi?AJAX=delete&ACT=$action&jobDir="+jobDir+"&libID="+libID,true);
    XHR.onreadystatechange=function (){
        if(XHR.readyState==4 && XHR.response){
            window.location.reload();
        }
    }
    XHR.send(null);
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

// <--AJAX
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">$title</FONT><BR><BR>
<INPUT type="button" value="Refresh window" onclick="window.location.reload()">&nbsp;&nbsp;<INPUT type="button" value="Close window" onclick="window.close()"><BR><BR><BR>
|;



###> recovering job names
my %jobDirList=&getJobList;
unless (scalar keys %{$jobDirList{$action}}){
    $dbh->disconnect;
	print qq
|<FONT class=\"title2\">No DIA process found.</FONT>
<SCRIPT language="JavaScript">
setTimeout('window.location.reload()',30000);
</SCRIPT>
</CENTER>
</BODY>
</HTML>
|;
	exit;
}

my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my @watchJob;
print "<TABLE>";
foreach my $jobDir (sort keys %{$jobDirList{$action}}){
    my $libID=$jobDirList{$action}{$jobDir};
    my ($libName,$libDesc,$useStatus,$user,$processType) = "";
    my $isDIA = 0;
    
    # GET INFO FILE
    my $fileInfo="$workDir/$jobDir/info_$libID.out";
    my %jobInfo;
    open(IN,"<$fileInfo");
    while(<IN>){
        my ($info,$value)=split(/=/,$_);
        $jobInfo{$info}=$value;
    }
    close IN;
    
    # GET LIBID (DIA) / DBID (TDA)
    if(length $libID && $jobInfo{'Software'} !~ /Skyline/) {
        $isDIA = 1;
    }
    
    if($libID) {
        if($isDIA) {
            ($libName,$useStatus) = $dbh->selectrow_array("SELECT NAME, USE_STATUS FROM SWATH_LIB WHERE ID_SWATH_LIB=$libID");
        } else {
            $libName = $dbh->selectrow_array("SELECT NAME FROM DATABANK WHERE ID_DATABANK=$libID");
        }
    }
    $libDesc = ($isDIA) ? "Library" : "Databank";
    
    if($action eq 'import') {
        $processType = ($isDIA) ? "Import DIA data" : "Import TDA data";
    } elsif($isDIA) {
        $processType = ($useStatus && $useStatus eq 'err')? 'Creating a library' : 'Updating a library';
    }
    
    # GET USER
    my $userName;
    if ($jobInfo{'User'}){
        chomp $jobInfo{'User'};
        $userName=$dbh->selectrow_array("SELECT USER_NAME FROM USER_LIST WHERE ID_USER=\"$jobInfo{'User'}\"");
    }
    
    print "<TR><TD><FONT class=\"title3\">&bull;&nbsp;&nbsp;Job: $jobDir&nbsp;&nbsp;-&nbsp;&nbsp;$processType";
    print "&nbsp;&nbsp;-&nbsp;&nbsp;Launched by user : $userName" if $userName;
    print "</FONT></TD></TR>";
    
    my ($status,$error)=&getJobStatus($jobDir,$libID);
    my $start=`head -1 $workDir/$jobDir/status_$libID.out`;
    
    print qq 
|<TR><TD nowrap>
    <TABLE cellspacing=0>
        <TR bgcolor="$darkColor">
            <TH class="rbBorder">#</TH>
            <TH class="rbBorder">&nbsp;&nbsp;$libDesc Name&nbsp;&nbsp;</TH>
            <TH class="rbBorder">&nbsp;&nbsp;Software&nbsp;&nbsp;</TH>
            <TH class="bBorder">&nbsp;&nbsp;Status&nbsp;&nbsp;</TH>
        </TR>
        <TR bgcolor="$lightColor">
            <TD>&nbsp;&nbsp;1&nbsp;&nbsp;</TD>
            <TD nowrap>&nbsp;&nbsp;$libName&nbsp;&nbsp;</TD>
            <TD nowrap>&nbsp;&nbsp;$jobInfo{'Software'}&nbsp;&nbsp;</TD>
            <TD nowrap style="padding:5px">
                &nbsp;&nbsp;<B>$start</B>&nbsp
                <SPAN id="status:$jobDir:$libID">$status</SPAN>&nbsp;&nbsp;
                <DIV id="error:$jobDir:$libID" style="display:none"><FIELDSET><LEGEND><B>Error message :</B></LEGEND>$error</FIELDSET></DIV>
            </TD>    
        </TR>
    </TABLE>
</TD></TR>    
|;
    
    push @watchJob,"$jobDir:$libID";
}
print "</TABLE></CENTER></BODY></HTML>";


if (@watchJob){
    my $watchJobStrg;
    foreach my $job (@watchJob){
        $watchJobStrg.=',' if $watchJobStrg;
        $watchJobStrg.=$job;
    }
    print qq |
<SCRIPT LANGUAGE="JavaScript">
var watchJobStrg="$watchJobStrg";
var autoUpdate=1;
function updateStatus(watchStatusStrg) {
    if (watchStatusStrg.match('##New job')){
        window.location.reload();
        return;
    }
    var statusStrg=watchStatusStrg.split('\\n');
    watchStatusStrg='';
    for (var i=0; i<statusStrg.length;i++){
        if (statusStrg[i].match('#')){
            var arrayStatus=statusStrg[i].split('##');
            if (arrayStatus[3].match('Finished')){
                window.location.reload();
                return;
            }
            document.getElementById('status:'+arrayStatus[1]+':'+arrayStatus[2]).innerHTML=arrayStatus[3];
            if (arrayStatus[4]){
                document.getElementById('error:'+arrayStatus[1]+':'+arrayStatus[2]).innerHTML='<FIELDSET><LEGEND><B>Error message :</B></LEGEND>'+arrayStatus[4]+'</FIELDSET>';
            }
            if (arrayStatus[0]==1){
                if(!watchStatusStrg.match(arrayStatus[1])){
                    if(watchStatusStrg){
                        watchStatusStrg+=',';
                    }
                    watchStatusStrg+=arrayStatus[1];       //jobDir
                }
                watchStatusStrg+=':'+arrayStatus[1];    //libID
            }
        }
    }
    if (watchStatusStrg && autoUpdate) {
		setTimeout('ajaxWatchStatus()',5000);
	}
}
// AJAX -->
function ajaxWatchStatus() {
    //If XHR object already exists, the request is canceled & the object is deleted
	if(XHR && XHR.readyState != 0){
		XHR.abort();
		delete XHR;
	}
    XHR=getXMLHTTP();
    if (!XHR){return false;}
    XHR.onreadystatechange=function() {
        if(XHR.readyState==4 && XHR.responseText){
            updateStatus(XHR.responseText);
        }
    }
    XHR.open("GET","./monitorDIAProcess.cgi?AJAX=update&ACT=$action&jobStrg="+watchJobStrg,true);
    XHR.send(null);
}
window.onload=setTimeout('ajaxWatchStatus()',5000);
// <-- AJAX
    
</SCRIPT>
    |; 
}


$dbh->disconnect;


sub ajaxUpdateStatus {
    my $jobStrg=param('jobStrg');
    my %watchJobList;
    foreach my $job (split(/,/,$jobStrg)){
        my ($jobDir,$libID)=split(/:/,$job);
        $watchJobList{$jobDir}=$libID;
    }
    
    ###> checking new job
    my %jobList=&getJobList;
    foreach my $jobDir (keys %{$jobList{$action}}){
        unless($watchJobList{$jobDir}){
            print "##New job ($jobDir)##\n";
            exit;
        }
    }
    
    foreach my $jobDir (keys %watchJobList){
        my $libID=$watchJobList{$jobDir};
        my ($status,$error)=&getJobStatus($jobDir,$libID);
        $error='' unless $error;
        if ($status eq 'Finished'){
            print "0##$jobDir##$libID##Finished\n";
        }
        else{
            if ($error && $error!~/WARNING/){
                print "0##$jobDir##$libID##$status##$error\n";
            }
            else{
                print "1##$jobDir##$libID##$status##$error\n";
            }
        }
    }
    exit;
}


sub getJobList{
    ###> recovering job names
    my %jobDirList;
    opendir (JOBDIR,$workDir);
    while (defined(my $jobDir = readdir(JOBDIR))){
        if ($jobDir=~/^\d{4}(_\d{2})+/){
            opendir(DIR,"$workDir/$jobDir");
            while (defined (my $jobFile= readdir (DIR))){
                if ($jobFile=~/status_(\d*)\.out/){
                    my $libID=$1;
                    $jobDirList{'library'}{$jobDir}=$libID;
                    last;
                }
            }
            close DIR;
        }
        elsif ($jobDir=~/^\d+/){
            opendir(DIR,"$workDir/$jobDir");
            while (defined (my $jobFile= readdir (DIR))){
                if ($jobFile=~/status_(\d*)\.out/){
                    my $libID=$1;
                    $jobDirList{'import'}{$jobDir}=$libID;
                    last;
                }
            }
            close DIR;
        }
    }
    close JOBDIR;
    return %jobDirList;
}


sub getJobStatus {
    my ($jobDir,$libID)=@_;
    my ($status,$warning);
    if (-s "$workDir/$jobDir/status_$libID.out"){
        if(`tail -1 $workDir/$jobDir/status_$libID.out`=~/Done/){
            $status='&nbsp;[&nbsp;<B>Finished</B>&nbsp;]&nbsp;';
        }
        elsif(`tail -1 $workDir/$jobDir/status_$libID.out`=~/Started/){
            $status='&nbsp;[&nbsp;Waiting for job to start...&nbsp;]&nbsp;';
        }
        elsif(`tail $workDir/$jobDir/status_$libID.out`=~/You had (\d+) "Peptide has unknown modification" errors/){
            $warning="Warning: You had $1 \"Peptide has unknown modification\" errors during SpectraST process.";
            $status='&nbsp;[&nbsp;Creating consensus library with unsplit mode in progress.&nbsp;]&nbsp;';
        }
        elsif(`tail -1 $workDir/$jobDir/status_$libID.out`=~/WARNING/){
            $warning=`tail -1 $workDir/$jobDir/status_$libID.out`;
            chomp $warning;
            $status='&nbsp;[&nbsp;';
            my $stat=`tail -2 $workDir/$jobDir/status_$libID.out`;
            $stat=~s/\n//g;
            my $qWarning=quotemeta($warning);
            ($status.=$stat)=~s/$qWarning//;
            $status.='&nbsp;]&nbsp;';
        }
        else{
            my $line=`tail -1 $workDir/$jobDir/status_$libID.out`;
            chomp $line;
            $status='&nbsp;[&nbsp;'.$line.'&nbsp;]&nbsp;';
        }
    }
    else{
        $status='Finished';
    }
    
    my $error='';
    if (-s "$workDir/$jobDir/status_$libID\_error.out"){
        open(ERROR, "$workDir/$jobDir/status_$libID\_error.out");
        while(my $line=<ERROR>){
            chomp $line;
            $error.=$line;
            $error.="<BR>";
        }
        close ERROR;
    }
    if ($error || $warning){
        my $tag=($error)? 'ERROR' : 'WARNING';
        $status.="<BR>&nbsp;&nbsp;<INPUT type=\"button\" value=\"Delete job\" onclick=\"deleteJob(\'$jobDir\','$libID')\">";
        $status.="&nbsp;<FONT color=\"red\"><B>**$tag**</B></FONT>&nbsp;";
        $status.="<INPUT type=\"button\" value=\"Show/Hide $tag\" onclick=\"displayErrorDiv(\'$jobDir\','$libID')\">";
        $error=$warning if $warning;
    }
    return ($status,$error);
}

sub ajaxDelete{
    my $jobDir=param('jobDir');
    my $libID=param('libID');
    if ($action eq 'library'){
        my $dbh=&promsConfig::dbConnect;
        $dbh->do("DELETE FROM SWATH_LIB WHERE ID_SWATH_LIB=$libID");
        $dbh->commit;
        $dbh->disconnect;
    }
    if (-e "$workDir/$jobDir"){
        system "rm -r $workDir/$jobDir";
    }
}

####>Revision history<####
# 1.0.4 Improved TDA data monitoring (VS 20/11/18)
# 1.0.3 Handling of TDA data monitoring  (VS 09/11/18)
# 1.0.2 Add import DIA data monitoring (MLP 08/02/18)
# 1.0.1 Allow peptide modification errors (MLP 30/01/18)
# 1.0.0 Creation (MLP 19/01/18)