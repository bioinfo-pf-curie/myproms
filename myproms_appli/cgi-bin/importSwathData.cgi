#!/usr/local/bin/perl -w

################################################################################
# importSwathData.cgi         1.3.8i                                           #
# Authors: M. Le Picard, P. Poullet & G. Arras (Institut Curie)                #
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

$|=1;       # buffer flush (allows instant printing)
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use strict;
use IO::Handle;
use promsConfig;
use promsMod;
use promsDIA;
use File::Copy;
use File::Basename;
use XML::Simple;
use POSIX qw(strftime); # to get the time
use File::Path qw(rmtree); # remove_tree
use String::Util qw(trim);
#use List::Util "first";
use File::Spec::Functions qw(splitpath); # Core module

my %promsPath=&promsConfig::getServerInfo('no_user');
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my $userID=(param('USERID'))?param('USERID'):'';

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect('no_user');
my $projectID=param('id_project');
my $branchID=param('ID');
my ($item,$experimentID)=split(':',$branchID);
print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);

my $action=(param('ACT'))? param('ACT') : "" ; ## import or quantification
my $format=(param('FORMAT'))? param('FORMAT') : ''; ##pkv, openswath, prm , spectronaut
my $submit=(param('submit'))? param('submit') : "";
my $fct=(param('FCT'))? param('FCT') : "" ;
my ($useFilteredFiles)=(param('FILTERED'))? param('FILTERED'):0;
my ($isCheckedFF,$dsExtension)=($useFilteredFiles==1)? ('checked','_filtered'):('','');
if ($action eq 'excludeLibraryModifications'){&excludeLibraryModifications; exit;}
if ($action eq "selectSamples"){&selectSamples; exit;}

my %clusterInfo=&promsConfig::getClusterInfo; #('debian'); # default is 'centos'
my $openMSPath=($clusterInfo{'on'}) ? $clusterInfo{'path'}{'openms'} : $promsPath{'openms'};
my $tricPath=($clusterInfo{'on'}) ? $clusterInfo{'path'}{'msproteomicstools'} : $promsPath{'msproteomicstools'};
my $pyprophetPath=($clusterInfo{'on'}) ? $clusterInfo{'path'}{'pyprophet'} : $promsPath{'pyprophet'};

if (param('CLUSTER') && param('CLUSTER')==1) {
    ###> Create new directory
    my $timeStamp=strftime("%Y%m%d%H%M%S",localtime);
    my $workdir=(param('workdir'))? param('workdir'):"$promsPath{data}/tmp/Swath/$timeStamp";
    mkdir $workdir unless -e $workdir;
    
    ####> Wait time
    print qq
    |<HTML><HEAD>
    <LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
    <STYLE type="text/css">
    </STYLE><BODY background="$promsPath{images}/bgProMS.gif"><CENTER>
        <BR><BR><IMG src='$promsPath{images}/engrenage.gif'><BR><BR>
    |;
    print "<FONT color=\"red\"><DIV id=\"waitDIV\"></DIV></FONT><BR>";
    print "<B>Processing</B><BR><BR>";
    
    my $scriptfile = "$workdir/openswathStoring.sh";
    my $currentDirectory=dirname $0;
    my $exportFileOut=(param('exportedlibrary'))? param('exportedlibrary') : '';
    my $tricMethod=(param('TRICMethod'))? param('TRICMethod') : 'BestOverall';
    my $paramFile=(param('exportparam'))?param('exportparam') : '';
    my $tsvfile=(param('uploadfile'))? param('uploadfile') : "feature_alignment_fileout.tsv";
    my $libID=(param('libID'))? param('libID') : '';
    my $softV=param('softversion');
    
    my ($fileName,$file);
    if ($tsvfile) {
        ###> Moving upload file
        my $newFile=tmpFileName($tsvfile);
        my $newFileDir="$workdir/$tsvfile";
        move($newFile,$newFileDir);
        $fileName = fileparse($newFileDir, qr/\.[^.]*/);
    }
    
    if ($exportFileOut){
        ###> upload export library file
        my $newLibExportFile=tmpFileName($exportFileOut);
        my $libExportFileDir="$workdir/$exportFileOut";
        move($newLibExportFile,$libExportFileDir);
    }
    
    if ($paramFile){
        my $newFile=tmpFileName($paramFile);
        my $newParamFile="$workdir/$paramFile";
        move ($newFile,$newParamFile);
    }
    
    open (BASH,"+>",$scriptfile);
	print BASH qq
|#!/bin/bash
export LC_ALL="C"
cd $currentDirectory
$clusterInfo{path}{perl}/perl $0 id_project=$projectID ID=$branchID ACT=$action FORMAT=$format submit=$submit workdir=$workdir CLUSTER=0 exportedlibrary=$exportFileOut TRICMethod=$tricMethod exportparam=$paramFile uploadfile=$tsvfile libID=$libID softversion=$softV USERID=$userID > $workdir/error.txt 2>&1
|;
	close BASH;
    my $modBash=0775;
	chmod $modBash, $scriptfile;
    my $size=(-e "$workdir/$tsvfile")?`stat -c "%s" $workdir/$tsvfile`: 1;
    $size=1073741824*15 unless $size;
    my $maxMem=($size/1073741824)*8;
    $maxMem=sprintf("%.0f",$maxMem);
    $maxMem=1 if ($maxMem < 1);
    $maxMem.='Gb';
    my $maxHours=48;
	my %jobParams=(
maxMem=>$maxMem,
numCPUs=>1,
maxHours=>$maxHours,
jobName=>"import_tric_$timeStamp",
outFile=>'PBSit.txt',
errorFile=>'PBSiterror.txt',
jobEndFlag=>"_END_$timeStamp",
	);

    my $childProcess=fork;
    unless ($childProcess){
        $clusterInfo{'runJob'}->($workdir,$scriptfile,\%jobParams);
        sleep 3;
        print qq |
<SCRIPT LANGUAGE="JavaScript">
var monitorWindow=window.open("$promsPath{cgi}/monitorDIAProcess.cgi?ACT=import",'Monitoring import of DIA data','width=1000,height=500,scrollbars=yes,resizable=yes');
monitorWindow.focus();
top.promsFrame.selectedAction='summary';
parent.navFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=experiment:$experimentID&ACT=nav";
</SCRIPT>
</BODY>
</HTML>
|;
        my $nbWhile=0;
        my $pbsError;
        while (`tail -1 $workdir/status_$libID.out`=~/Done/) {
            if ($nbWhile>$maxHours*60*2) {
                die "importswath is taking too long or died before completion";
            }
            sleep 30;
            $pbsError=$clusterInfo{'checkError'}->("$workdir/PBSiterror.txt");
            $nbWhile++;
        }
        exit;
    }    
    exit;
}

#######################
####>Starting HTML<####
#######################
###> Form to select files
if ($submit eq "") {
    print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.checklist {
    border: 1px solid #ccc;
    list-style:none;
    height:100px;
    min-width:250px;
    overflow:auto;

}
.checklist, .checklist LI {margin:0; padding:0;}
.checklist LABEL {
    padding-left:3px;
    text-indent:-25px;
}
.liRow1 {background:#FFFFFF;}
.liRow2 {background:#D8D8D8;}
.popup {background-color:#FFFFFF;border:solid 3px #999999;padding:5px;box-shadow:10px 10px 20px #808080;position:absolute;display:none;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">|;
    &promsMod::popupInfo();
    &promsMod::browseDirectory_JavaScript if $promsPath{'shared'};
    print qq |
function cancelAction(){
	window.location="./processAnalyses.cgi?ID=$branchID";
}

function selectExportOption (value){
    if (value == 'export'){
        document.getElementById('exportparameters').style.display='';
        document.getElementById('exportparamfile').style.display='none';
        document.getElementById('importlibexport').style.display='none';
    }
    else if (value == 'import'){
        document.getElementById('exportparameters').style.display='none';
        document.getElementById('exportparamfile').style.display='';
        document.getElementById('importlibexport').style.display='';
    }
}

function checkform (form){
    if(document.getElementById('exportparamfile').style.display=='' && form.exportparam.files.length==0){
        alert('ERROR: Select the export parameters file.');
        return false;
    }
    if(document.getElementById('importlibexport').style.display=='' && form.exportedlibrary.files.length==0){
        alert('ERROR: Select the library export file.');
        return false;
    }
    if(document.getElementById('exportparameters').style.display=='' && form.swathfile.files.length==0){
        alert('ERROR: Select the swath windows file.');
        return false;
    }

}
function selectSource (source){
    var myForm=document.parameters;
    myForm.mzXMLFiles.disabled=true;
    if (document.getElementById('sharedDirDIV')) {document.getElementById('sharedDirDIV').style.display='none';}

    if(source=='UseLocalDirectory'){
        myForm.mzXMLFiles.disabled=false;
    }
    else if(source=='UseSharedDirectory'){
        document.getElementById('sharedDirDIV').style.display='';
    }
}
function updateExpList(status){
    if(status==false){ // use “_with_dscore.csv” for TRIC
        window.location="./importSwathData.cgi?id_project=$projectID&ID=$branchID&ACT=quantification&FORMAT=openswath&FILTERED=0";
    }
    else{ // use “_with_dscore_filtered.csv” for TRIC
        window.location="./importSwathData.cgi?id_project=$projectID&ID=$branchID&ACT=quantification&FORMAT=openswath&FILTERED=1";
    }
}



// AJAX --->
    function getXMLHTTP() {
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
function selectLibrary (libID){
    document.getElementById('libraryModifications').style.display='';
    var XHR=getXMLHTTP();
    XHR.onreadystatechange=function () {
        if(XHR.readyState==4 && XHR.responseText ){
            document.getElementById('libraryModifications').innerHTML=XHR.responseText;
        }
    }
    XHR.open("GET","$promsPath{cgi}/importSwathData.cgi?ACT=excludeLibraryModifications&ID=$branchID&libID="+libID,true);
    XHR.send(null);
}

function ajaxSelectSamples (expID){
    var listSample=document.getElementById('listSamples'+expID);
    if (listSample.style.display==''){
        var i=1;
        while(document.getElementById('samples'+expID+i)){
            document.getElementById('samples'+expID+i).checked=false;
            i++;
        }
    }
    else{
        listSample.style.display='';
        var XHR=getXMLHTTP();
        XHR.onreadystatechange=function () {
            if(XHR.readyState==4 && XHR.responseText ){
                document.getElementById('listSamples'+expID).innerHTML=XHR.responseText;
            }
        }
        XHR.open("GET","$promsPath{cgi}/importSwathData.cgi?ACT=selectSamples&ID=$branchID&expID="+expID,true);
        XHR.send(null);
    }
}
// <--- AJAX

</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
    <CENTER>
|;
    if ($action eq 'import'){
        my $formTitle=($format eq 'pkv')? "Peakview" : ($format eq 'prm') ? "Skyline" : ($format eq 'openswath') ? "OpenSwath" : "Spectronaut" ;
        my $fileFormat=($format eq 'pkv')? ".xlsx" : ($format eq 'prm') ? ".csv,.tsv" : ($format eq 'openswath') ? ".tsv" : ".tsv" ;
        my $clusterValue=($format eq 'openswath')? 1:0;
        print qq |
            <FONT class="title1">Select options to import $formTitle data</FONT><BR><BR><BR><BR>
            <FORM method="POST" action="./importSwathData.cgi" name="parameters" enctype="multipart/form-data">
                <INPUT name="ID" type="hidden" value="$branchID">
                <INPUT name="id_project" type="hidden" value="$projectID">
                <INPUT name="ACT" type="hidden" value="$action">
                <INPUT name="FORMAT" type="hidden" value="$format">
                <INPUT name="CLUSTER" type="hidden" value="$clusterValue">
                <INPUT name="USERID" type="hidden" value="$userID">
                <TABLE bgcolor="$darkColor">
                <TR><TH align=right valign=top>Result file : </TH><TD bgcolor="$lightColor"><INPUT type="file" name="uploadfile" accept="$fileFormat" required>
        |;

        if($format eq 'pkv' || $format eq 'openswath' || $format eq 'spectronaut' ){
            print qq |
                </TD></TR><TR><TH align=right valign=top>Library name : </TH><TD bgcolor="$lightColor">|;
                my $sthLibList=$dbh->prepare("SELECT NAME,ID_SWATH_LIB,IDENTIFIER_TYPE FROM SWATH_LIB WHERE USE_STATUS='yes'") or die "Couldn't prepare statement: " . $dbh->errstr;
                $sthLibList->execute;
                my %libList;
                if ($sthLibList->rows==0) {print "No library available.";}
                else{
                    while (my($libName,$libID,$identType)=$sthLibList->fetchrow_array){
                        $identType=($identType)? $identType : 'Not defined';
                        $libList{$libName}="$libID&$identType";
                    }
                    $sthLibList->finish;
                    $dbh->commit;

                    print qq
                    |<SELECT name="libID" id="libID" required>
                         <option value="">-= Select Library =-</option>
                    |;
                    foreach my $libName (sort {lc($a) cmp lc($b)} keys %libList) {
                        my ($libID,$identType)=split('&',$libList{$libName});
                        print "<option value=\"$libID\">$libName [$identType]</option>";
                    }
                    print " </SELECT> ";
                }
                print qq |
                </TD></TR>
                <TR><TH align=right valign=top>Library export parameter file : </TH><TD bgcolor="$lightColor"><INPUT type="file" name="exportparam" required></TD></TR>
                |;
                print "<TR><TH align=right valign=top>Library export file : </TH><TD bgcolor=\"$lightColor\"><INPUT type=\"file\" name=\"exportedlibrary\" accept=\".tsv\" required ></TD></TR>" unless $format eq 'pkv';
                print qq |
                <TR><TH  align=right valign=top>FDR : </TH><TD bgcolor="$lightColor"><INPUT type="number" name="fdr" min="0" max="20" size="2" required> %</TD></TR>
                <TR><TH  align=right valign=top>PeakView parameters : </TH>
                    <TD bgcolor="$lightColor"><B>Peptide filter :</B><BR>
                    Number of peptides per protein : <INPUT type="text" name="nbpeptide" size="3"><BR>
                    Number of transitions per peptide : <INPUT type="text" name="nbtransition" size="3"><BR>
                    <INPUT type="checkbox" name="modpep">Exclude modified peptides<BR><BR>
                    <B>XIC options : </B><BR>
                    XIC extraction window (min) : <INPUT type="text" name="XICmin" size="3"><BR>
                    XIC width : <INPUT type="text" name="XIC" size="3">
                    <SELECT name="XICunit">
                        <OPTION value="DA">Da</OPTION>
                        <OPTION value="PPM">ppm</OPTION>
                    </SELECT>

                    </TD>
                </TR>| if $format eq 'pkv';
                print qq |
                <TR ><TH align=right valign=top>TRIC method used: </TH>
                    <TD bgcolor="$lightColor">
                        <SELECT name="TRICMethod" required >
                            <OPTION value="LocalMST" selected>LOCAL MST*</OPTION>
                            <OPTION value="BestOverall">BEST OVERALL</OPTION>
                        </SELECT>
                        <SMALL>*Recommanded option</SMALL>
                    </TD>
                </TR>| if $format eq 'openswath';
        }
        elsif($format eq 'prm'){
            print qq |
                <SMALL><BR><BR><U>The columns names have to be like</U> : <BR>
                Protein Name, Replicate Name, File Name, Peptide Modified Sequence, Precursor Mz<BR>
                Precursor Charge, Product Mz, Product Neutral Mass, Product Charge, Fragment Ion<BR>
                Retention Time, Missed Cleavages, Peptide Retention Time, Area, Begin Pos<BR>
                End Pos, Previous Aa, Next Aa<BR><BR></SMALL>
                </TD></TR>
            |;
            print "<TR><TH align=right valign=top>Fasta file : </TH><TD bgcolor=\"$lightColor\">";
            my $sthDBList=$dbh->prepare("SELECT D.NAME,ID_DATABANK,FASTA_FILE,DT.NAME FROM DATABANK D, DATABANK_TYPE DT WHERE DT.ID_DBTYPE=D.ID_DBTYPE AND USE_STATUS='yes'") or die "Couldn't prepare statement: " . $dbh->errstr;
            $sthDBList->execute;
            my %dbList;
            if ($sthDBList->rows==0) {print "No library available.";}
            else{
                while (my($dbName,$dbID,$fastaFile,$dbType)=$sthDBList->fetchrow_array){
                    my ($dbSource,$version);
                    if($fastaFile=~/:/){
                        my ($mascotServer,$dbankDir,$fileName)=split(':',$fastaFile);
                        $dbSource=$mascotServer;
                        $version=$dbankDir;
                    }
                    else{
                        if(! -e "$promsPath{banks}/db_$dbID/$fastaFile"){next;}
                        $dbSource="Local";
                    }
                    $dbList{$dbSource}{$dbID}=$dbName;
                    $dbList{$dbSource}{$dbID}.=" ($version)" if $version;
                    $dbList{$dbSource}{$dbID}.=" [$dbType]";
                }
                $sthDBList->finish;
                $dbh->commit;

                print qq
                    |<SELECT name="dbID" id="dbID" required>
                         <option value="">-= Select databank =-</option>
                |;
                foreach my $dbSource (sort {$a cmp $b} keys %dbList){
                    print "<OPTGROUP label=\"$dbSource :\">";
                    foreach my $dbID (sort {lc($dbList{$dbSource}{$a}) cmp lc($dbList{$dbSource}{$b})} keys %{$dbList{$dbSource}}) {
                        print "<OPTION value=\"$dbID\">$dbList{$dbSource}{$dbID}</OPTION>";
                    }
                    print "</OPTGROUP>";
                }
                print " </SELECT> ";
            }
            print "</TD></TR><TR><TH align=right valign=top>Search file :</TH><TD bgcolor=\"$lightColor\">";
            my %modificationList;
            my $sthModificationList=$dbh->prepare("SELECT PSI_MS_NAME,ID_MODIFICATION FROM MODIFICATION WHERE MONO_MASS IS NOT NULL AND PSI_MS_NAME IS NOT NULL");
            $sthModificationList->execute;
            while(my ($modName,$modID)=$sthModificationList->fetchrow_array){
                next unless $modName && $modID;
                $modificationList{$modID}=$modName;
            }
            $sthModificationList->finish;
            print "<SELECT name=\"modif\" id=\"modif\"  required multiple>";
            foreach my $modID (sort {$modificationList{$a} cmp $modificationList{$b}} keys %modificationList){
                print "<option value=\"$modID\">$modificationList{$modID}</option>";
            }
            print "</SELECT></TD></TR>";

        }
        print "<TR><TH align=right valign=top>Software version : </TH><TD bgcolor=\"$lightColor\"><INPUT type=\"text\" name=\"softversion\" required><SMALL>&nbsp;&nbsp;ex : 1.2, 2.1.3 ...</SMALL></TD></TR>";
        print qq |
                </SELECT></TD></TR>
                <TR ><TH colspan=2><input type="submit" name="submit" value="Submit">
                <!-- CLEAR button -->
                &nbsp &nbsp &nbsp<INPUT type="reset" value="Clear" />
                <!-- CANCEL button -->
                &nbsp &nbsp &nbsp<INPUT value="Cancel" onclick="cancelAction();" type="button" ></TH></TR>
            </TABLE>
            </FORM>
        |;
    }
    elsif ($action eq 'quantification'){
        print qq |
        <BR><FONT class="title1">OpenSwath quantification</FONT><BR><BR><BR>
        <FORM method="POST" name="parameters" action="./importSwathData.cgi" enctype="multipart/form-data" onsubmit="return(checkform(this));">
            <INPUT name="ID" type="hidden" value="$branchID">
            <INPUT name="id_project" type="hidden" value="$projectID">
            <INPUT name="ACT" type="hidden" value="$action">
            <INPUT name="FORMAT" type="hidden" value="$format">
            <TABLE bgcolor="$darkColor">
                <TR ><TH align=right valign=top>TRIC method: </TH>
                    <TD bgcolor="$lightColor">
                        <SELECT name="TRICMethod">
                            <OPTION value="LocalMST" selected>LOCAL MST*</OPTION>
                            <OPTION value="BestOverall">BEST OVERALL</OPTION>
                        </SELECT>
                        <SMALL>*Recommanded option</SMALL><BR>
                        &nbsp;<INPUT type="checkbox" value="$useFilteredFiles" onmouseover="popup('recommended for any large sample size')" onmouseout="popout()" name="FILTERED" id="useFilteredDScores" onchange="updateExpList(this.checked)" $isCheckedFF>Run TRIC on <I>“_with_dscore_filtered.csv”</I> files. 
                    </TD>
                </TR>
                <TR><TH align=right valign=top>Library name : </TH><TD bgcolor="$lightColor">&nbsp;|;
                my $sthLibList=$dbh->prepare("SELECT NAME,ID_SWATH_LIB,IDENTIFIER_TYPE FROM SWATH_LIB WHERE USE_STATUS='yes'") or die "Couldn't prepare statement: " . $dbh->errstr;
                $sthLibList->execute;
                my %libList;
                if ($sthLibList->rows==0) {print "No library available.";}
                else{
                    while (my($libName,$libID,$identType)=$sthLibList->fetchrow_array){
                        $identType=($identType)? $identType : 'Not defined';
                        $libList{$libName}="$libID&$identType";
                    }
                    $sthLibList->finish;
                    $dbh->commit;

                    print qq
                    |<SELECT name="libID" id="libID" onchange="selectLibrary(this.value);" required>
                         <option value="">-= Select Library =-</option>
                    |;
                    foreach my $libName (sort {lc($a) cmp lc($b)} keys %libList) {
                        my ($libID,$identType)=split('&',$libList{$libName});
                        print "<option value=\"$libID\">$libName [$identType]</option>";
                    }
                    print " </SELECT> ";
                }
                print qq |
                </TD></TR>
                <TR><TH align=right valign=top>Library export management :</TH>
                    <TD bgcolor="$lightColor">
                        <INPUT type="radio" name="libraryfileoption" value="export" onclick="selectExportOption(this.value)" required>Use the selected library<INPUT type="radio" name="libraryfileoption" value="import" onclick="selectExportOption(this.value)">Import your own library formated for Openswath<BR>
                    </TD>
                </TR>
                <TR id="exportparamfile" style="display:none"><TH align=right valign=top>Library export parameter file : </TH>
                    <TD bgcolor="$lightColor"><INPUT type="file" name="exportparam" ></TD>
                </TR>
                <TR id="importlibexport" style="display:none"><TH align=right valign=top>Library export file : </TH>
                    <TD bgcolor="$lightColor"><INPUT type="file" name="exportedlibrary" accept=".tsv"></TD>
                </TR>
                <TR id="exportparameters" style="display:none"><TH align=right valign=top>Library export options : </TH>
                    <TD bgcolor="$lightColor">
                    <B>&nbsp;- Mass range of fragment ions : </B><BR>
                        &nbsp;Min : <INPUT type="text" name="-lmin" size="5" value="350">&nbsp;&nbsp;Max : <INPUT type="text" name="-lmax" size="5" value="2000"><BR>

                    <BR><B>&nbsp;- Ion series and charge : </B><BR>
                        &nbsp;Ions* : <INPUT type="text" name="-s" size="5" value="b,y">
                        &nbsp;Charge : <INPUT type="text" name="-x" size="5" value="1,2"><BR><SMALL>&nbsp;*(separated by ',') for example : 'b,y' </SMALL><BR>

                    <BR><B>&nbsp;- Number of ions per peptide : </B><BR>
                        &nbsp;Min : <INPUT type="text" name="-o" size="4" value="3"> Max :<INPUT type="text" name="-n" value="6" size="4"><BR>

                    <BR><B>&nbsp;- Files : </B><BR>
                        &nbsp;Windows SWATH file : <INPUT onmouseover="popup('File containing the swath ranges.')" onmouseout="popout()" type="file" name="swathfile" ><BR>
                        &nbsp;File with modifications delta mass* : <INPUT onmouseover="popup('File with the modifications not specified by default.')" onmouseout="popout()" type="file" name="-m" ><BR>
                        &nbsp;<SMALL>*File headers : modified-AA TPP-nomenclature Unimod-Accession ProteinPilot-nomenclature is_a_labeling composition-dictionary<BR>
                        &nbsp;An example : S S[167] 21 [Pho] False {'H':1,'O':3,'P':1}</SMALL><BR>
                        &nbsp;Labelling file : <INPUT onmouseover="popup('File containing the amino acid isotopic labelling mass shifts.')" onmouseout="popout()" type="file" name="-i"><BR>
                        &nbsp;Fasta file : <INPUT onmouseover="popup('Fasta file to relate peptides to their proteins.')" onmouseout="popout()" type="file" name="-f" accept=".fasta"><BR>

                    <BR><B>&nbsp;- Other options : </B><BR>
                        &nbsp;<INPUT type="checkbox" name="other" value="-d">Remove duplicate masses from labelling<BR>
                        &nbsp;<INPUT type="checkbox" name="other" value="-e">Use theoretical mass<BR>
                        &nbsp;Time scale : <INPUT type="radio" name="-t" value="seconds" checked>seconds<INPUT type="radio" name="timescale" value="minutes">minutes<BR>
                        &nbsp;UIS order :  <INPUT onmouseover="popup('When using a switching modification, this determines the UIS order to be calculated.')" onmouseout="popout()" type="text" name="-y" value="2" size="3"><BR>
                        &nbsp;Maximum permissible error : <INPUT onmouseover="popup('Maximum error allowed at the annotation of a fragment ion.')" onmouseout="popout()" type="text" name="-p" value="0.05" size="4"><BR>
                        &nbsp;Allowed fragment mass modifications : <INPUT type="text" name="-g" size="10" onmouseover="popup('List of allowed fragment mass modifications. Useful for phosphorylations.')" onmouseout="popout()"><BR>

                    <DIV id="libraryModifications" style="display:none"></DIV>
                </TD></TR>
                <TR><TH align=right valign=top>OpenSwath parameters files : </TH>
                    <TD bgcolor="$lightColor">
                    &nbsp;iRT file : <INPUT type="file" name="iRTFile" accept=".traML" required><BR>
                    &nbsp;Windows file : <INPUT type="file" name="windowsFile" accept=".txt" required><BR>
                    </TD>
                </TR>
                <TR><TH align=right valign=top>mzXML files : </TH>
                    <TD bgcolor="$lightColor">
                        &nbsp;<INPUT type="radio" name="selSource" value="UseLocalDirectory" onclick="selectSource(this.value);"><B>Upload multiple files</B><BR>&nbsp;&nbsp;&nbsp;<INPUT type="file" name="mzXMLFiles" accept=".mzXML" required multiple="multiple" disabled><BR><BR>
                        &nbsp;<INPUT type="radio" name="selSource" value="UseSharedDirectory" onclick="selectSource(this.value);"><B>Import from shared data directory </B><BR>
                        &nbsp;&nbsp;&nbsp;&nbsp;<DIV id="sharedDirDIV" style="width:600px;max-height:300px;overflow:auto;display:none">
                        |;
                        &promsMod::browseDirectory_getFiles($promsPath{'shared'},'sharedDirFiles',{fileMatch=>qr/\.(dat|xml|mzXML)$/i});
                        print qq |
                        </DIV>
                    </TD>
                </TR>
                <TR><TH align=right valign=top>OpenSwath workflow options: </TH>
                    <TD bgcolor="$lightColor">
                    &nbsp;mz_threshold : <INPUT type="text" name="mzThreshold" size="3" value="0.05">
                    </TD>
                </TR>
                <TR><TH align=right valign=top>Pyprophet options: </TH>
                    <TD bgcolor="$lightColor">
                    &nbsp;d_score.cutoff : <INPUT type="text" name="dscore" size="3" value="1">
                    </TD>
                </TR>
                <TR id="TRICoption" ><TH align=right valign=top>Merge with other experiment: </TH>
                    <TD  bgcolor="$lightColor">|;
                    my $sthExperiment=$dbh->prepare("SELECT NAME,ID_EXPERIMENT FROM EXPERIMENT WHERE ID_PROJECT=$projectID");
                    $sthExperiment->execute();
                    if($sthExperiment->rows==0){
                        print "<SCRIPT langage=\"javascript\">document.getElementById('TRICoption').style.display='none';</SCRIPT>";
                    }
                    else{
                        my $sthAnaID=$dbh->prepare("SELECT ID_ANALYSIS,A.NAME FROM SAMPLE S,ANALYSIS A WHERE S.ID_SAMPLE=A.ID_SAMPLE AND S.ID_EXPERIMENT=?");
                        my $validExperiment=0;
                        while (my ($expName,$expID)=$sthExperiment->fetchrow_array){
                            $sthAnaID->execute($expID);
                            my $checkPyprophetFile=0;
                            while (my ($anaID,$anaName)=$sthAnaID->fetchrow_array){
                                my $anaDir="$promsPath{peptide}/proj_$projectID/ana_$anaID";
                                if (-s "$anaDir/$anaName\_with_dscore$dsExtension.csv"){
                                    $checkPyprophetFile=1;
                                    last;
                                }
                                elsif (-s "$anaDir/$anaName\_with_dscore$dsExtension.csv.tar.gz"){
                                    $checkPyprophetFile=1;
                                    last;
                                }
                            }
                            print "<INPUT type=\"checkbox\" name=\"tricmerge\" value=\"$expID\" onclick=\"ajaxSelectSamples(this.value)\">$expName&nbsp;<BR><DIV id=\"listSamples$expID\" style=\"display:none\" ></DIV>" if $checkPyprophetFile;
                            $validExperiment=1 if $checkPyprophetFile;
                        }
                        if (! $validExperiment){
                            print "<SCRIPT langage=\"javascript\">document.getElementById('TRICoption').style.display='none';</SCRIPT>";
                        }
                    }
                    print qq |
                    </TD>
                </TR>
                <TR ><TH colspan=2><input type="submit" name="submit" value="Submit">
                <!-- CLEAR button -->
                &nbsp &nbsp &nbsp<INPUT type="reset" value="Clear" />
                <!-- CANCEL button -->
                &nbsp &nbsp &nbsp<INPUT value="Cancel" onclick="cancelAction();" type="button" ></TH></TR>
            </TABLE>
        </FORM>
        <DIV id="divDescription" class="clDescriptionCont">
        <!--Empty div-->
        </DIV>
        <SCRIPT type="text/javascript">setPopup()</SCRIPT>
        |;
    }
    print qq |
    </CENTER>
</BODY>
</HTML>
|;
}
else{
    #############################
    ####>Fetching parameters<####
    #############################
    my %massATave=&promsConfig::getMassATave; # Average mass, needed for protein mass calculation
    my %massAAave=&promsConfig::getMassAAave; # Average mass, needed for protein mass calculation

    ###> Create new directory
    my $startTime=strftime("%s",localtime);
    my $time=strftime("%Y%m%d%H%M%S",localtime);
    my $workDir=(param('workdir'))? param('workdir'):"$promsPath{data}/tmp/Swath/$time";
    mkdir $workDir unless -e $workDir;

    ###> Wait time
    print qq
    |<HTML><HEAD>
    <LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
    <STYLE type="text/css">
    </STYLE><BODY background="$promsPath{images}/bgProMS.gif"><CENTER>
        <BR><BR><IMG src='$promsPath{images}/engrenage.gif'><BR><BR>
    |;
    print "<FONT color=\"red\"><DIV id=\"waitDIV\"></DIV></FONT><BR>";
    print "<B>Processing</B><BR><BR>";

    ############################
    ###>Fetching parameters <###
    ############################
    my $paramVersion=param('softversion');
    my ($softwareVersion)=&promsMod::cleanParameters($paramVersion);

    ###> Fetching project ID
    my ($protVisibility,$identProject)=$dbh->selectrow_array("SELECT PROT_VISIBILITY, ID_IDENTIFIER FROM PROJECT WHERE ID_PROJECT=$projectID");
    $protVisibility=0 unless $protVisibility;

    my ($fileName,$file);
    unless ($action eq 'quantification'){
        if (param('uploadfile')) {
            $file=param('uploadfile');
            ###> Moving upload file
            my $newFile=tmpFileName($file);
            my $newFileDir="$workDir/$file";
            move($newFile,$newFileDir);
            $fileName = fileparse($newFileDir, qr/\.[^.]*/);
        }
    }

    my $libID=(param('libID'))? param('libID') : '';

    my $exportFileOut=(param('exportedlibrary'))? param('exportedlibrary') : '';
    if ($exportFileOut){
        ###> upload export library file
        my $newLibExportFile=tmpFileName($exportFileOut);
        my $libExportFileDir="$workDir/$exportFileOut";
        move($newLibExportFile,$libExportFileDir);
    }


    ###> Fetching DB file
    my @dbID;
    if($libID){
        my $sthDBID=$dbh->prepare("SELECT ID_DATABANK,DB_RANK FROM DATABANK_SWATHLIB WHERE ID_SWATH_LIB=? ");
        $sthDBID->execute($libID);
        while (my ($dbID,$dbrank)=$sthDBID->fetchrow_array) {
            push @dbID,"$dbID&$dbrank";
        }
        $sthDBID->finish;
    }else{
        my @db=split(/&/,param('dbID'));
        push @dbID,"$db[0]&1";
    }

   my (@mzXmlFileList,$irtFile,$swathWindows,$w,$m);
   my $exportDir="$promsPath{data}/swath_lib/SwLib_$libID";
   if ($action eq 'quantification'){
        print "<BR><BR>Uploading files ...";
        ###> uploading files
        if (param('mzXMLFiles')){
            foreach my $mzXML (param('mzXMLFiles')){
                push @mzXmlFileList, $mzXML;
                my $newFile=tmpFileName($mzXML);
                my $newFileDir="$workDir/$mzXML";
                move($newFile,$newFileDir);
                
                my $size=`stat -c "%s" $newFileDir`;
                if ($size == 0 ) {
                    print "<HTML><BODY><SCRIPT LANGUAGE=\"JavaScript\">alert('Error during importation of upload files. Please retry.'); window.location=\"$promsPath{cgi}/importSwathData.cgi?id_project=$projectID&ID=$branchID&ACT=quantification&FORMAT=openswath\";</SCRIPT></BODY></HTML>";
                    exit;
                }
                print '.';
            }
        }
        elsif (param('sharedDirFiles')){
            foreach my $sharedFile (param('sharedDirFiles')){
                my (undef,$path,$fileName)=splitpath($sharedFile);
                 push @mzXmlFileList, $fileName;
                
                ##>move shared file in work directory
                my $NewFileDir="$workDir/$fileName";
                move("$promsPath{shared}/$sharedFile",$NewFileDir);
                
                my $size=`stat -c "%s" $NewFileDir`;
                if ($size == 0 ) {
                    print "<HTML><BODY><SCRIPT LANGUAGE=\"JavaScript\">alert('Error during importation of upload files. Please retry.'); window.location=\"$promsPath{cgi}/importSwathData.cgi?id_project=$projectID&ID=$branchID&ACT=quantification&FORMAT=openswath\";</SCRIPT></BODY></HTML>";
                    exit;
                }
                print '.';
            }
        }
            
        $irtFile=param('iRTFile');
        my $newIrtFile=tmpFileName($irtFile);
        my $newiRTFileDir="$workDir/$irtFile";
        move($newIrtFile,$newiRTFileDir);

        $swathWindows=param('windowsFile');
        my $swathWindowsFile=tmpFileName($swathWindows);
        my $swathWindowsNewDir="$workDir/$swathWindows";
        move($swathWindowsFile,$swathWindowsNewDir);
        
        
        ###> Move Swath windows file into workDir
        if (param('swathfile')){
            $w=upload('swathfile');
            my $swathFile=tmpFileName($w);
            my $swathFileDir="$exportDir/$w";
            move($swathFile,$swathFileDir);
        }
        ###> Move modifications file into workDir
        if (param('-m')){
            $m=upload('-m');
            my $newFile=tmpFileName($m);
            my $newDir="$exportDir/$m";
            move($newFile,$newDir);
        }
        
    }
    
    if (param('exportparam')){
        my $paramFile=param('exportparam');
        my $newFile=tmpFileName($paramFile);
        my $newParamFile="$workDir/$paramFile";
        move ($newFile,$newParamFile);
    }

    my $fileInfo="$workDir/info_$libID.out";
    my $software=($format eq 'pkv')? "Peakview" : ($format eq 'prm') ? "Skyline" : ($format eq 'openswath') ? "OpenSwath" : "Spectronaut" ;
	open(FILEINFO,">$fileInfo") || die "Error while opening $fileInfo";
	print FILEINFO "User=$userID\nSoftware=$software";
	close FILEINFO;



    my $fileStat="$workDir/status_$libID.out";
	open(FILESTAT,">$fileStat") || die "Error while opening $fileStat";
	print FILESTAT "Started ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
	close FILESTAT;

	my $fileError="$workDir/status_$libID\_error.out";

    my $anaReloadID;
#    my $childProcess=fork;
#    unless ($childProcess){
#        #>Disconnecting from server
#        open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
#		open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
#		open STDERR, '>/dev/null' or die "Can't open /dev/null: $!";

my $libraryVersion;
my $errorTRICMethod=0;
if ($action eq 'quantification'){
    $dbh=&promsConfig::dbConnect('no_user');
    $libraryVersion=$dbh->selectrow_array("SELECT VERSION_NAME FROM SWATH_LIB WHERE ID_SWATH_LIB=$libID") unless $libraryVersion;
    $softwareVersion="2.2";

    ###> recovering other pyprophet files
    my @pyprophetFilesOut;
    if(param('tricmerge') || param('samples')){
        open(FILESTAT,">>$fileStat");
        print FILESTAT "Recovering pyprophet files.\n";
        close FILESTAT;

        my @anaIDList;
        my $sthAnaID=$dbh->prepare("SELECT ID_ANALYSIS,NAME FROM ANALYSIS WHERE ID_SAMPLE=?");
        foreach my $sampleID (param('samples')){
            $sthAnaID->execute($sampleID);
            while (my ($anaID,$anaName)=$sthAnaID->fetchrow_array){
                push @anaIDList,$anaID;
                my $anaDir="$promsPath{peptide}/proj_$projectID/ana_$anaID";
                my $pyprophetFile="$anaName\_with_dscore$dsExtension.csv";
                if (-s "$anaDir/$pyprophetFile"){
                    copy("$anaDir/$pyprophetFile","$workDir/$pyprophetFile");
                    #push @pyprophetFilesOut,"$workDir/$pyprophetFile";
                    push @pyprophetFilesOut,$pyprophetFile;
                }
                elsif (-s "$anaDir/$pyprophetFile.tar.gz"){
                    copy("$anaDir/$pyprophetFile.tar.gz","$workDir/$pyprophetFile.tar.gz");
                    system "cd $workDir; tar -zxf  $pyprophetFile.tar.gz";
                    #push @pyprophetFilesOut,"$workDir/$pyprophetFile";
                    push @pyprophetFilesOut,$pyprophetFile;
                }
                else{
                    open(FILEERROR,">>$fileError");
                    print FILEERROR "No pyprophet file found for analysis $anaName ($anaID).";
                    close(FILEERROR);
                    exit;
                }
            }
        }
        $sthAnaID->finish;
    }

        ###> TRIC options
        open(FILESTAT,">>$fileStat");
        print FILESTAT "Fetching parameters.\n";
        close FILESTAT;
        my $tricOptions;
        my $tricMethod=param('TRICMethod');
        if ($tricMethod eq 'BestOverall'){
            $tricOptions='--method best_overall  --realign_method diRT';
        }
        else{
            $tricOptions='--method LocalMST --mst:useRTCorrection True --mst:Stdev_multiplier 3.0 --alignment_score 0.0001 --realign_method lowess --dscore_cutoff 1.0 --frac_selected 0 --disable_isotopic_grouping ';
        }

        ########################
        ###> Library Export <###
        ########################
        my $numValidProt;
        my ($libraryName)=$dbh->selectrow_array("SELECT NAME FROM SWATH_LIB WHERE ID_SWATH_LIB=$libID");
        if ($exportFileOut){
            $numValidProt=0;
        }
        else{
            my ($p,$g,$f,$t,$y,$labelling)=(param('-p'),param('-g'),param('-f'),param('-t'),param('-y'),param('-i'));
            my $other;
            foreach my $val (param('other')){
                $other.="$val ";
            }
            my ($lMax,$lMin,$s,$x,$n,$o);
            if ($other){
                ($lMax,$lMin,$s,$x,$n,$o,$other)=&promsMod::cleanParameters(param('-lmax'),param('-lmin'),param('-s'),param('-x'),param('-n'),param('-o'),$other) ;
            }
            else{($lMax,$lMin,$s,$x,$n,$o)=&promsMod::cleanParameters(param('-lmax'),param('-lmin'),param('-s'),param('-x'),param('-n'),param('-o'));}
            my $k='openswath';

            my $pepMod=param('pepMod');
            $pepMod=($pepMod)? $pepMod : '';

            ###> loading div
            open(FILESTAT,">>$fileStat");
            print FILESTAT "Library conversion for OpenSwath.\n";
            close FILESTAT;
            my $loadingDivID='none';
            my $loadingSPAN='none';
            ($exportFileOut,my $paramFile,$numValidProt,my $error)=&promsDIA::exportLibrarySCRIPT($dbh,$libID,$libraryName,$exportDir,\%promsPath,$startTime,$loadingDivID,$loadingSPAN,$w,$k,$p,$g,$f,$t,$y,$m,$labelling,$other,$lMax,$lMin,$s,$x,$n,$o,$pepMod);
            if ($error){
                open(FILEERROR,">>$fileError");
                print FILEERROR "ERROR: Library conversion failed : $error";
                close(FILEERROR);
                exit;
            }
            move("$exportDir/$exportFileOut","$workDir/$exportFileOut");
            move("$exportDir/$libraryName\_param","$workDir/$libraryName\_param");
        }
        $dbh->disconnect; # disconnect before long process without DB requirement
        
        ###########################
        ###> Running OpenSwath <###
        ###########################
        open(FILESTAT,">>$fileStat");
        print FILESTAT "Running OpenSwath quantification.\n";
        close FILESTAT;


        ###> fetching parameters
        my $mzThreshold=param('mzThreshold');
        my $dscore=param('dscore');
        my $outFile=$workDir.'/openSwathOUT.txt';



        ##########################################
        ###> library conversion for openswath <###
        ##########################################
        open(FILESTAT,">>$fileStat");
        print FILESTAT "Converting DIA library.\n";
        close FILESTAT;

        open(SCRIPT, "+>","$workDir/libConv.sh");
        print SCRIPT qq
|#!/bin/bash \n
$openMSPath/TargetedFileConverter -in $workDir/$exportFileOut -out $workDir/$libraryName.TraML  >>$outFile 2>&1 ;\n
$openMSPath/OpenSwathDecoyGenerator -method shuffle -decoy_tag DECOY_ -mz_threshold $mzThreshold -exclude_similar -similarity_threshold 0.05 -append -in $workDir/$libraryName.TraML -out $workDir/$libraryName\_DECOY.TraML -remove_unannotated >>$outFile 2>&1 ;\n
$openMSPath/TargetedFileConverter -in $workDir/$libraryName\_DECOY.TraML -out $workDir/$libraryName\_DECOY.tsv >>$outFile 2>&1;\n\n
|;
        close SCRIPT;
        if ($clusterInfo{'on'}) {
            my $clusterCommandString=$clusterInfo{'buildCommand'}->($workDir,"$workDir/libConv.sh");
            my $bashFile="$workDir/libraryConversion.sh";
            my $libSize=`stat -c "%s" $workDir/$exportFileOut`;
            my $maxHours=int(24);
            my $maxMem=($libSize/1073741824)*10;
            $maxMem=($maxMem < 20) ? 20 : ($maxMem > 50) ? 50 : sprintf("%.0f",$maxMem) ;
            $maxMem.='Gb';
            open(BASH1,">$bashFile");
            print BASH1 qq
|#!/bin/bash
##resources
#PBS -l mem=$maxMem
#PBS -l nodes=1:ppn=1
#PBS -l walltime=$maxHours:00:00
#PBS -q batch

##Information
#PBS -N libraryConversion_$time
#PBS -M marine.le-picard\@curie.fr
#PBS -m abe
#PBS -o $workDir/PBS.txt
#PBS -e $workDir/PBSerror.txt

## Command
$clusterCommandString
echo End_Library_Conversion
|;
            close BASH1;

            system "chmod 775 $workDir/libConv.sh";
            $clusterInfo{'sendToCluster'}->($bashFile);

            ### wait loop
            my $count=0;
            my ($error,$endProcess,$PBSerror);
            my $nbWhile=0;
            my $maxNbWhile=$maxHours*60*2;
            while (!$endProcess && !$error && !$PBSerror){
                if(-s $outFile){
                    $error=`grep -i 'error\\\|exception' $outFile`;
                }
                if(-s "$workDir/PBS.txt"){
                    $endProcess=`grep 'End_Library_Conversion' $workDir/PBS.txt`;
                }
                sleep 30;
                $nbWhile++;
                if ($nbWhile > $maxNbWhile) {
                    open(FILEERROR,">>$fileError");
                    print FILEERROR "Aborting: Library conversion is taking too long.";
                    close(FILEERROR);
                    exit;
                }
                $PBSerror=$clusterInfo{'checkError'}->("$workDir/PBSerror.txt") if -s "$workDir/PBSerror.txt";
            }
            $error=`grep -i 'error\\\|exception' $outFile` unless $error;
            if ($error){
                open(FILEERROR,">>$fileError");
                print FILEERROR $error;
                close(FILEERROR);
                exit;
            }
            if ($PBSerror){
                open(FILEERROR,">>$fileError");
                print FILEERROR $PBSerror;
                close(FILEERROR);
                exit;
            }
            system "rm $workDir/PBSerror.txt; rm $workDir/PBS.txt; rm $workDir/torqueID.txt; rm $workDir/ssh_connect.sh";
        }
        else{
            system "bash $workDir/libConv.sh&";
            my ($nbProcess,$count)=(0,0);
            my $pbsError;
            while ($nbProcess!=2 && !$pbsError){
                if(-s $outFile){
                    $pbsError=`grep -i 'error\\\|exception' $outFile`;
                    $nbProcess=`grep -c 'TargetedFileConverter took' $outFile`;
                }
                sleep 30;
            }
            if ($pbsError){
                open(FILEERROR,">>$fileError");
                print FILEERROR $pbsError;
                close(FILEERROR);
                exit;
            }
        }
        my $decoyLibrary=$libraryName.'_DECOY.tsv';


        ##########################################################
        ###> openswath workflow and pyprophet on search files <###
        ##########################################################
        open(FILESTAT,">>$fileStat");
        print FILESTAT "Processing OpenSwath and Pyprophet.\n";
        close FILESTAT;

        my $i=1;
        if ($clusterInfo{'on'}) {
            my %maxLaunch;
            my $maxHours=int(24);
            foreach my $mzXMLFile (sort {$a cmp $b} @mzXmlFileList) {
                (my $mzXMLout=$mzXMLFile)=~s/.mzXML$/.tsv/;
                my $mzxmlSize=`stat -c "%s" $workDir/$mzXMLFile`;
                open(FILESTAT,">>$fileStat");
                print FILESTAT "Analysing $mzXMLFile.\n";
                close FILESTAT;
                my $dir=$workDir.'/OpenSwath_'.$i;
                mkdir $dir unless -e $dir;
                $maxLaunch{$i}=0;
                my $outputFile=$dir.'/openSwathOUT.txt';
                open(SCRIPT, "+>","$dir/OpenSwath.sh");
                print SCRIPT qq
|#!/bin/bash

$openMSPath/OpenSwathWorkflow -in $workDir/$mzXMLFile  -tr $workDir/$decoyLibrary -sort_swath_maps -batchSize 1000 -tr_irt $workDir/$irtFile -swath_windows_file $workDir/$swathWindows -out_tsv $workDir/$mzXMLout >>$outputFile 2>&1 ;
$pyprophetPath/pyprophet --d_score.cutoff=\"$dscore\" --ignore.invalid_score_columns $workDir/$mzXMLout >>$outputFile 2>&1 ;\n
|;
                close SCRIPT;

                my $clusterCommandString=$clusterInfo{'buildCommand'}->($dir,"$dir/OpenSwath.sh");
                my $bashFile="$dir/OpenSwathRun.sh";
                #my $maxHours=int(24);
                my $maxMem=($mzxmlSize/1073741824)*5;
                $maxMem=($maxMem < 30) ? 30 : ($maxMem >100) ? 100 : sprintf("%.0f",$maxMem);


                $maxMem.='Gb';
                open(BASH1,">$bashFile");
                print BASH1 qq
|#!/bin/bash

##resources
#PBS -l mem=$maxMem
#PBS -l nodes=1:ppn=1
#PBS -l walltime=$maxHours:00:00
#PBS -q batch

##Information
#PBS -N OpenSwath_Workflow_$time\_$i
#PBS -M marine.le-picard\@curie.fr
#PBS -m abe
#PBS -o $dir/PBS.txt
#PBS -e $dir/PBSerror.txt

## Command
$clusterCommandString
echo End_OpenSwath_Workflow_$i
|;
                close BASH1;

                system "chmod 775 $dir/OpenSwath.sh";
                $clusterInfo{'sendToCluster'}->($bashFile);

                sleep 5;
                $i++;
            }

            ### wait loop
            my ($outError,$PBSerror);
            my %endProcess;
            my $count=0;
            my $nbWhile=0;
            my $maxNbWhile=$maxHours*60*2;
            while (scalar keys %endProcess < scalar @mzXmlFileList && !$PBSerror && !$outError){
                for (my $j=1; $j<=scalar @mzXmlFileList;$j++){
                    next if $endProcess{$j};
                    if (-s "$workDir/OpenSwath_$j/PBS.txt"){
                        my $process=`tail $workDir/OpenSwath_$j/PBS.txt`;
                        next unless $process;
                        $endProcess{$j}=1 if $process=~/End_OpenSwath_Workflow_$j/;
                    }
                    if (-s "$workDir/OpenSwath_$j/PBSerror.txt"){
                        $PBSerror=$clusterInfo{'checkError'}->("$workDir/OpenSwath_$j/PBSerror.txt") if -s "$workDir/OpenSwath_$j/PBSerror.txt";
                        if ($PBSerror=~/ERROR  : There was an error binding the path \/mchevail\/ACULM: No such file or directory/ && $maxLaunch{$j} < 5){
                            system "rm $workDir/OpenSwath_$j/PBSerror.txt; rm $workDir/OpenSwath_$j/PBS.txt; rm $workDir/OpenSwath_$j/ssh_connect.sh ; rm $workDir/OpenSwath_$j/torqueID.txt";
                            my $bashFile=$workDir.'/OpenSwath_'.$j.'/OpenSwathRun.sh';
                            $clusterInfo{'sendToCluster'}->($bashFile);
                            $PBSerror='';
                        }
                        $maxLaunch{$j}++;
                    }

                    $outError=`grep -i 'error\\\|exception' $workDir/OpenSwath_$j/openSwathOUT.txt` if (-s "$workDir/OpenSwath_$j/openSwathOUT.txt" && `grep -ci 'error\\\|exception' $workDir/OpenSwath_$j/openSwathOUT.txt`!=0);
                    last if ($PBSerror || $outError);
                }
                sleep 30;
                $nbWhile++;
                if ($nbWhile > $maxNbWhile) {
                    open(FILEERROR,">>$fileError");
                    print FILEERROR "Aborting: mzXML files processing is taking too long.\n";
                    close(FILEERROR);
                    exit;
                }
            }

            if ($outError){
                open(FILEERROR,">>$fileError");
                print FILEERROR $outError;
                close(FILEERROR);
                exit;
            }
            if ($PBSerror){
                open(FILEERROR,">>$fileError");
                print FILEERROR $PBSerror;
                close(FILEERROR);
                exit;
            }
        }
        else{
            ###> openswath workflow and pyprophet on search files
            my $i=1;
            open(SCRIPT, "+>","$workDir/OpenSwathRun.sh");
            foreach my $mzXMLFile (sort {$a cmp $b} @mzXmlFileList) {
                (my $mzXMLout=$mzXMLFile)=~s/.mzXML$/.tsv/;
                $i++;
                print SCRIPT qq
|$openMSPath/OpenSwathWorkflow -in $workDir/$mzXMLFile  -tr $workDir/$decoyLibrary -sort_swath_maps -batchSize 1000 -tr_irt $workDir/$irtFile -swath_windows_file $workDir/$swathWindows -out_tsv $workDir/$mzXMLout >>$outFile 2>&1 ;
$pyprophetPath/pyprophet --d_score.cutoff=\"$dscore\" --ignore.invalid_score_columns $workDir/$mzXMLout >>$outFile 2>&1 ;\n
|;
            }
            close SCRIPT;
            system "bash $workDir/OpenSwathRun.sh&";
            my $outError;
            my $endProcess=0;
            while (!$outError && !$endProcess){
                if (-s "$workDir/$outFile"){
                   $outError=`grep -i 'error\\\|exception' $workDir/$outFile`;
                }
                foreach my $mzXMLFile (@mzXmlFileList){
                    (my $mzXMLName=$mzXMLFile)=~s/.mzXML$//;
                    if (`ls $workDir/$mzXMLName* | wc -l` == 15 ){
                        $endProcess=1;
                    }
                    else{$endProcess=0;}
                    last if $endProcess==0;
                }
                sleep 30;
            }
            if ($outError){
                open(FILEERROR,">>$fileError");
                print FILEERROR $outError;
                close(FILEERROR);
                exit;
            }
        }
        sleep 30;

        foreach my $file (@mzXmlFileList) {
            $file=~s/.mzXML/_with_dscore$dsExtension.csv/;
            push @pyprophetFilesOut,"$file";
        }
        
        # Compute maxSize on all feature_alignment.py infiles
        my $maxSize=0;
        foreach my $file (@pyprophetFilesOut) {
            my $size=`stat -c "%s" $workDir/$file`;
            $maxSize+=$size;
        }
        sleep 30;

        ################################
        ###> TRIC on iprophet files <###
        ################################
        open(FILESTAT,">>$fileStat");
        print FILESTAT "Running TRIC on data.\n";
        close FILESTAT;
        my $command="cd $workDir; $tricPath/feature_alignment.py --in @pyprophetFilesOut --out feature_alignment_fileout.tsv --file_format openswath $tricOptions --max_rt_diff 30 --target_fdr 0.01 --max_fdr_quality 0.05 >>$outFile 2>&1";
        if ($clusterInfo{'on'}){
            my $clusterCommandString=$clusterInfo{'buildCommand'}->($workDir,$command);
            my $maxHours=int(24);
            my $maxMem=($maxSize/1073741824);
            #$maxMem=($maxMem < 20) ? 20 : ($maxMem > 90) ? 90 : sprintf("%.0f",$maxMem) ;
            $maxMem=($maxMem < 50) ? 50 : ($maxMem > 100) ? 100 : sprintf("%.0f",$maxMem);

            $maxMem.='Gb';
            my $bashFile="$workDir/TRICProcessus.sh";
            open(BASH1,">$bashFile");
            print BASH1 qq
|#!/bin/bash

##resources
#PBS -l mem=$maxMem
#PBS -l nodes=1:ppn=1
#PBS -l walltime=$maxHours:00:00
#PBS -q batch

##Information
#PBS -N OpenMS_TRIC_$time
#PBS -M marine.le-picard\@curie.fr
#PBS -m abe
#PBS -o $workDir/PBS.txt
#PBS -e $workDir/PBSerror.txt

## Command
$clusterCommandString
echo End_TRIC
|;
            close BASH1;

            ### lauch bash script
            $clusterInfo{'sendToCluster'}->($bashFile);

            my $count=0;
            my ($error,$PBSerror,$endProcess);
            my $nbWhile=0;
            my $maxNbWhile=$maxHours*60*2;
            while(!$endProcess && !$PBSerror && !$error){
                $PBSerror=$clusterInfo{'checkError'}->("$workDir/PBSerror.txt") if -s "$workDir/PBSerror.txt";
                if(-s $outFile){
                    $error=`grep -i 'error\\\|exception' $outFile`;
                }
                if(-s "$workDir/PBS.txt"){
                    $endProcess=`grep 'End_TRIC' $workDir/PBS.txt`;
                }
                sleep 30;
                $nbWhile++;
                if ($nbWhile > $maxNbWhile) {
                    open(FILEERROR,">>$fileError");
                    print FILEERROR "Aborting: TRIC is taking too long.\n";
                    close(FILEERROR);
                    exit;
                }
            }
            $error=`grep -i 'error\\\|exception' $outFile` unless $error;
            if ($error){
                if ($error=~/Exception: No data available for alignment \d_\d vs \d_\d/){
                    $errorTRICMethod=1;
                }else{
                    open(FILEERROR,">>$fileError");
                    print FILEERROR $error;
                    close(FILEERROR);
                    exit;
                }
            }
            if ($PBSerror){
                open(FILEERROR,">>$fileError");
                print FILEERROR $PBSerror;
                close(FILEERROR);
                exit;
            }
        }
        else{
            system "bash $command&";
            my $pbsError;
            my $count=0;
            while(!-s "$workDir/feature_alignment_fileout.tsv" && !$pbsError){
                if(-s $outFile){
                    $pbsError=`grep -i 'error\\\|exception' $outFile`;
                }
                sleep 30;
            }
            $pbsError=`grep -i 'error\\\|exception' $outFile` unless $pbsError;
            if ($pbsError){
                if ($pbsError=~/Exception: No data available for alignment 0_0 vs 0_1/){
                    $errorTRICMethod=1;
                }else{
                    open(FILEERROR,">>$fileError");
                    print FILEERROR $pbsError;
                    close(FILEERROR);
                    exit;
                }
            }
        }



        ####> If TRIC method selected by the user don't work -> lauch the other method
        if ($errorTRICMethod){
            if($clusterInfo{'on'}){
                system "rm $workDir/PBS.txt";
                system "rm $workDir/PBSerror.txt";
                system "rm $workDir/torqueID.txt";

                open(FILESTAT,">>$fileStat");
                print FILESTAT "WARNING : TRIC method ($tricMethod) did not work on data. ";
                print FILESTAT "Running the other method : BestOverall.\n" if $tricMethod eq 'LocalMST';
                print FILESTAT "Running the other method : LocalMST.\n" if $tricMethod eq 'BestOverall';
                close FILESTAT;

                my $outFile2=$workDir.'/Tric2OUT.txt';
                my $newTricOptions;
                if ($tricMethod eq 'BestOverall'){
                    $newTricOptions='--method LocalMST --mst:useRTCorrection True --mst:Stdev_multiplier 3.0 --alignment_score 0.0001 --realign_method lowess --dscore_cutoff 1.0 --frac_selected 0 --disable_isotopic_grouping ';
                }
                else{
                    $newTricOptions='--method best_overall  --realign_method diRT';
                }
                my $command="cd $workDir; $tricPath/feature_alignment.py --in @pyprophetFilesOut --out feature_alignment_fileout.tsv --file_format openswath $newTricOptions --max_rt_diff 30 --target_fdr 0.01 --max_fdr_quality 0.05 >>$outFile2 2>&1";

                my $clusterCommandString=$clusterInfo{'buildCommand'}->($workDir,$command);
                my $maxHours=24;
                my $maxMem=($maxSize/1073741824);
                $maxMem=($maxMem < 50) ? 50 : ($maxMem > 100) ? 100 : sprintf("%.0f",$maxMem);

                $maxMem.='Gb';
                my $bashFile="$workDir/scriptBashTRIC.sh";
                open(BASH2,">$bashFile");
                 print BASH2 qq
|#!/bin/bash

##resources
#PBS -l mem=$maxMem
#PBS -l nodes=1:ppn=1
#PBS -l walltime=$maxHours:00:00
#PBS -q batch

##Information
#PBS -N OpenMS_TRIC2_$time
#PBS -M marine.le-picard\@curie.fr
#PBS -m abe
#PBS -o $workDir/PBS.txt
#PBS -e $workDir/PBSerror.txt

## Command
$clusterCommandString
echo End_OpenMS
|;
                close BASH2;

                ### lauch bash script
                $clusterInfo{'sendToCluster'}->($bashFile);
                ### wait loop
                my ($error,$PBSerror,$endProcess);
                my $nbWhile=0;
                my $maxNbWhile=$maxHours*60*2;
                my $j=1;
                while(!$endProcess && !$PBSerror && !$error){
                    $PBSerror=$clusterInfo{'checkError'}->("$workDir/PBSerror.txt") if -s "$workDir/PBSerror.txt";
                    if(-s $outFile2){
                        $error=`grep -i 'error\\\|exception' $outFile2`;
                    }
                    if(-s "$workDir/PBS.txt"){
                        $endProcess=`grep 'End_OpenMS' $workDir/PBS.txt`;
                    }
                    sleep 30;
                    $nbWhile++;
                    if ($nbWhile > $maxNbWhile) {
                        open(FILEERROR,">>$fileError");
                        print FILEERROR "Aborting: TRIC is taking too long.";
                        close(FILEERROR);
                        exit;
                    }
                }
                if ($error || $PBSerror){
                    open(FILEERROR,">>$fileError");
                    print FILEERROR "TRIC method dit not work on the data. <BR>Select other TRIC options or method.";
                    close(FILEERROR);
                    exit;
                }
            }
        }
        
        

        ###> check if the final file exist
        if(-s "$workDir/feature_alignment_fileout.tsv"){$file="feature_alignment_fileout.tsv";}
        else{
            open(FILEERROR,">>$fileError");
            print FILEERROR "OpenSwath workflow did not work.";
            close(FILEERROR);
            exit;
        }

        open(FILESTAT,">>$fileStat");
        print FILESTAT "Import quantification results.\n";
        close FILESTAT;

    }


    ######################################################################
    ####>Fetching list of existing proteins for corresponding project<####
    ######################################################################
    $dbh=&promsConfig::dbConnect('no_user');
    ####>Fetching all project's proteins (they will be substracted from the new list => no duplicates)
    my (%projectProt,%peptideModif,%proteinLength,%incompleteProtein,%noDescriptionProtein,%proteinDescription,%noSequenceProtein,%noMWProtein);
    open(FILESTAT,">>$fileStat");
    print FILESTAT "Scanning project's proteins.\n";
    close FILESTAT;
    my $sthP=$dbh->prepare("SELECT ID_PROTEIN,IDENTIFIER,PROT_LENGTH,PROT_DES,PROT_SEQ,MW FROM PROTEIN WHERE ID_PROJECT=$projectID");
    $sthP->execute;
    my $refProtData=$sthP->fetchall_arrayref; # reference to an array
    $sthP->finish;
    $dbh->disconnect;
    my $nbLine=0;
    foreach my $refProt (@{$refProtData}) {
        next unless $refProt->[1]; # identifier is null (previous ID protection not removed)
        $projectProt{$refProt->[1]}=$refProt->[0]; # {identifier}=id
        $proteinLength{$refProt->[0]}=$refProt->[2]; #{id}=length
        $proteinDescription{$refProt->[0]}=$refProt->[3] unless !$refProt->[3]; #{id}=description

        if (!$refProt->[2]) { # no length
            $incompleteProtein{$refProt->[1]}=$refProt->[0]; # {identifier}=id
            if (!$refProt->[3] || $refProt->[3]=~/no\sdescription/i) { # no description
                $noDescriptionProtein{$refProt->[1]}=$refProt->[0]; # {identifier}=id
            }
        }
        if (!$refProt->[4]) { # no sequence
            $noSequenceProtein{$refProt->[1]}=$refProt->[0];    # {identifier}=id
        }
        if (!$refProt->[5]) { # no MW
            $noMWProtein{$refProt->[1]}=$refProt->[0];  # {identifier}=id
        }
    }


    my (%targetGhostNb,%decoyGhostNb);


    #########################################################
    ###> Fetching library and modifications informations <###
    #########################################################
    $dbh=&promsConfig::dbConnect('no_user');
    my (%massMods,%unimod2ID,%peakview2modID,%modifIDs,%peptideInfo,%fragMass,%modificationList);
    my ($fdr,$versionName,$taxonomy,$instrument,$quantifAnnot);
    if ($format eq 'prm'){
        ###> Fetching modifications mass and unimod ID
        my $sthMassMods=$dbh->prepare("SELECT ID_MODIFICATION, MONO_MASS FROM MODIFICATION WHERE MONO_MASS IS NOT NULL");
        $sthMassMods->execute;
        while(my ($modID,$massMod)=$sthMassMods->fetchrow_array){
            $massMods{$modID}=$massMod;
        }
        $sthMassMods->finish;
        ## recover modifications selected by the user
        foreach my $modID (param('modif')){
            $modificationList{$modID}=$massMods{$modID};
        }
    }
    else{
        ###> Fetching SWATH LIB modifications
        (my $libName,$instrument,$taxonomy,$versionName,my $split,my $stat)=$dbh->selectrow_array("SELECT NAME,INSTRUMENT,ORGANISM,VERSION_NAME,SPLIT,STATISTICS FROM SWATH_LIB WHERE ID_SWATH_LIB='$libID' ") or die "Couldn't prepare statement: " . $dbh->errstr;
        my $sthModifIDs=$dbh->prepare("SELECT M.ID_MODIFICATION, SLM.SPECIFICITY, MONO_MASS, UNIMOD_ACC, PEAKVIEW_CODE FROM MODIFICATION M, SWATH_LIB_MODIFICATION SLM WHERE M.ID_MODIFICATION=SLM.ID_MODIFICATION AND SLM.ID_SWATH_LIB=$libID");
        $sthModifIDs->execute();
        while (my ($modID,$specif,$massMod,$unimodID,$peakviewCode)=$sthModifIDs->fetchrow_array){
            my $modType="V";
            my @specifs=split(/,/,$specif);
            foreach my $aa (@specifs){
                $modifIDs{$modID}{$modType}{$aa}=1;
            }
            $massMods{$modID}=$massMod if $massMod;
            $peakview2modID{$peakviewCode}=$modID if $peakviewCode;
            next unless $unimodID && $modID;
            $unimod2ID{$unimodID}=$modID if $format eq 'openswath';
        }
        $sthModifIDs->finish;
        ###> Fetching SWATH LIB statistics
        my $xml=XML::Simple->new(KeepRoot=>1);
        my $xmlStat=$xml->XMLin($stat);
        my $pepNum=$xmlStat->{'NUM_ENTRY'}->{'NUM_PEP'};
        my $protNum=$xmlStat->{'NUM_ENTRY'}->{'NUM_PROT'};
        my $protSpeNum=$xmlStat->{'NUM_ENTRY'}->{'NUM_PROT_SPECIFIC'};
        if($format eq 'pkv'){
            my $XICunit=(param('XICunit'))? param('XICunit') : "N/A";
            ($fdr,my $nbPeptide,my $nbTransition,my $modifiedPeptide,my $XICwindow,my $XICwidth)=&promsMod::cleanParameters(param('fdr'),param('nbpeptide'),param('nbtransition'),param('modpep'),param('XICmin'),param('XIC'));
            $nbPeptide=$nbPeptide? $nbPeptide : 'N/A';
            $nbTransition=$nbTransition? $nbTransition : "N/A";
            $XICwindow=$XICwindow? $XICwindow : "N/A";
            $XICwidth=$XICwidth? $XICwidth : "N/A";
            my $exclusePep=($modifiedPeptide) ? 1 : 0;
            $quantifAnnot="LABEL=FREE::SOFTWARE=PKV;".$softwareVersion."::NB_PEPTIDES_PER_PROTEIN=".$nbPeptide."::NB_TRANSITION_PER_PEPTIDE=".$nbTransition."::EXCLUDE_MODIFIED_PEPTIDES=".$exclusePep."::XIC_EXTRACTION_WINDOW=".$XICwindow."::XIC_WIDTH=$XICwidth\_$XICunit";
            $quantifAnnot.="::LIBRARY_NAME=".$libName."::LIBRARY_PEPTIDES=".$pepNum."::LIBRARY_PROTEINS=".$protNum."::LIBRARY_SPECIFICS_PROTEINS=".$protSpeNum;
        }
        elsif ($format eq 'openswath'){
            $quantifAnnot="LABEL=FREE::SOFTWARE=OS;".$softwareVersion."::LIBRARY_NAME=".$libName."::LIBRARY_PEPTIDES=".$pepNum."::LIBRARY_PROTEINS=".$protNum."::LIBRARY_SPECIFICS_PROTEINS=".$protSpeNum;
            #  -> OS => OPENSWATH
        }
        else{
            $quantifAnnot="LABEL=FREE::SOFTWARE=SPC;".$softwareVersion."::LIBRARY_NAME=".$libName."::LIBRARY_PEPTIDES=".$pepNum."::LIBRARY_PROTEINS=".$protNum."::LIBRARY_SPECIFICS_PROTEINS=".$protSpeNum;
            #  -> SPC => SPECTRONAUT
        }
        #######################################
        ### Manage export library parameter ###
        #######################################
        my $paramFile="$workDir/$libName\_param";
        open(PARAM,"<$paramFile");
        while (<PARAM>){
            next if $_=~/(Desired proteins not found into the library|Export for)/;
            my ($paramName,$param)=split(/:/,$_);
            $param=~s/\s//g;
            $quantifAnnot.="::$paramName=$param";
        }
        close PARAM;
        if ($format eq 'openswath'){
            $quantifAnnot.=(param('TRICMethod') eq 'BestOverall' && $errorTRICMethod) ? "::TRIC_METHOD=LocalMST" : (param('TRICMethod') eq 'LocalMST' && $errorTRICMethod) ? "::TRIC_METHOD=BestOverall" : "::TRIC_METHOD=".param('TRICMethod');
        }
        ##########################################################
        ####>Recovering peptide miss cleavage and specificity<####
        ##########################################################
        my $count=0;
        open(FILESTAT,">>$fileStat");
        print FILESTAT "Scanning library.\n";
        close FILESTAT;
        print ".";
        if (-e "$promsPath{swath_lib}/SwLib_$libID/$libName.sptxt"){
            open(LIBFILE,"<","$promsPath{swath_lib}/SwLib_$libID/$libName.sptxt") or die ("open: $!");
            my $nameLine;
            my (%mod2Unimod,%pepModif);
            while (my $line=<LIBFILE>){
                if ($line=~/^Name: (.+)/) { #new peptide line
                    ($nameLine=$1)=~s/\s|n//g;
                    if($nameLine=~/\[\d*\]/ && $format eq 'pkv'){   #peptide mofification
                        while($nameLine=~/(\w)\[(\d*)\]/g){   #conversion of library modification code into peakview modification code (for peptide matches)
                            my $modMass=$2;
                            my $modifID;
                            my $aa=($1) ? $1 : '';
                            my $pos=$-[0];
                            foreach my $ID (keys %massMods){
                                if(($aa && sprintf("%.0f",$massMods{$ID}+$massAAave{$1}) == $2) || sprintf("%.0f",$massMods{$ID}) == $2){
                                    $modifID=$ID;
                                }
                                last if $modifID;
                            }

                            my $peakviewCode;
                            foreach my $pCode (keys %peakview2modID){
                                $peakviewCode=$pCode if $peakview2modID{$pCode}==$modifID;
                            }
                            $nameLine=~s/\[$2\]/\[$peakviewCode\]/;
                        }
                        $nameLine=~s/\//_/;
                    }
                    #elsif($nameLine=~/\[\d*\]/ && $format eq 'openswath' && $nameLine!~/^\[43\]/){    ### n-ter modifications are not take into account in peakview, or openswath
                    elsif($nameLine=~/\[\d*\]/ && ($format eq 'openswath' || $format eq 'spectronaut')){    ### n-ter modifications are not take into account in peakview, or openswath
                        while ($nameLine=~/(\w)*\[(\d*)\]/g){
                            my $modMass=$2;
                            my $modifID;
                            my $aa=($1) ? $1 : '';
                            foreach my $ID (keys %massMods){
                                if(($aa && sprintf("%.0f",$massMods{$ID}+$massAAave{$1}) == $2) || sprintf("%.0f",$massMods{$ID}) == $2){
                                    $modifID=$ID;
                                }
                                elsif ($modMass==43 && $massMods{$ID}==42.0106){
                                    $modifID=$ID;
                                }
                                last if $modifID;
                            }

                            my $modID=($mod2Unimod{$2}) ? $mod2Unimod{$2} : $dbh->selectrow_array("SELECT UNIMOD_ACC FROM MODIFICATION WHERE ID_MODIFICATION=$modifID");
                            $mod2Unimod{$2}=$modID unless $mod2Unimod{$2};
                            if ($format eq 'openswath'){
                                $nameLine=~s/\[$modMass\]/(UniMod:$modID)/;
                            }
                            else{
                                my $mass=sprintf("%.0f",$massMods{$modifID});
                                $nameLine=~s/\[$modMass\]/\[\+$mass\]/;
                            }

                        }
                        $nameLine=~s/\//_/;
                    }
                    else{$nameLine=~s/\//_/;}
                }
                elsif($line=~/^PrecursorMZ:(.*)/){  # M/Z
                    (my $mrObs=$1)=~s/\s//g;
                    push @{$peptideInfo{$nameLine}},$mrObs;
                }
                elsif($line=~/^Comment: (.+)/){
                    my @commentInfo=split(/ /,$1);
                    foreach my $info (@commentInfo){
                        if ($info=~/^NMC=(.*)/) {       #miss cleavage
                            my $missCut=($1 eq '') ? 'N/A' : $1;
                            push @{$peptideInfo{$nameLine}},$missCut;
                        }
                        elsif($info=~/^Protein=([0-9])\/(.+)/) {    #list of protein where the peptide was found
                            my $nbProt=($split==0) ? $1 : ($2=~/Subgroup_[0-9]_([0-9]).*/) ? $1 : '';
                            my $specificity=100/$nbProt;
                            $specificity=($specificity) ? (sprintf("%0.1f", $specificity))*1 : '-';
                            push @{$peptideInfo{$nameLine}},$specificity;
                        }
                        elsif($info=~/^Mods=(\d*)\/(.*)/){
                            my @result=&promsMod::convertVarModStringSwath($2);    #fetching peptide modification(s)
                            foreach my $res (@result){
                                my @resTab=split(/!/,$res);         # varMod!residues!positions
                                next if ($resTab[0]=~m/Label/);
                                my %hashMod;
                                my $modID=&promsMod::getModificationIDfromString($dbh,$resTab[0],$resTab[1],\%hashMod);       #fetching modification ID
                                $pepModif{$nameLine}{$modID}{$resTab[2]}=$massMods{$modID};
                            }
                        }
                    }
                }
                #elsif($format eq 'openswath' && $line=~/^\d/ && $nameLine!~/\[43\]/){
                #elsif($format eq 'openswath' && $line=~/^\d/){
                #    my ($fragMass,$fragIntensity,$fragInfo,@infos)=split(/\t/,$line);
                #    my ($fragmentType,$fragmentCharge);
                #    if($fragInfo=~/^\[?(b\d*|y\d*)\^?(\d)?\/.*/){
                #        $fragmentType=$1;
                #        $fragmentCharge=($2) ? $2 : 1;
                #    }
                #    $fragMass{$nameLine}{$fragmentType.'_'.$fragmentCharge}=$fragMass if ($fragmentCharge && $fragmentType);
                #}
            }
            close LIBFILE;
        }
        else{
             open(FILEERROR,">>$fileError");
            print FILEERROR "Library file not found.";
            close(FILEERROR);
            exit;
        }


        if (($format eq 'openswath') && ( $action eq 'import' || ($action eq 'quantification' && $exportFileOut))){
            if (-s "$workDir/$exportFileOut"){
                open (LIBEXPORT,"<$workDir/$exportFileOut");
                while (<LIBEXPORT>){
                    next if $.==1;
                    my ($mzPep,$fragMZ,$recalib,$transition,@other)=split(/\t/,$_);
                    my ($transitionID,$frag,$fragCharge,$peptide,$pepCharge)=split(/_/,$transition);
                    $fragMass{$peptide.'_'.$pepCharge}{$frag.'_'.$fragCharge}=$fragMZ;
                }
                close LIBEXPORT;
            }
            else{
                open(FILEERROR,">>$fileError");
                print FILEERROR "No library export file found.";
                close(FILEERROR);
                exit;
            }
        }
    }
    $dbh->disconnect;


    #############################
    ###> Parsing upload file <###
    #############################
    my (%peptideProt,%peptidePosition,%fragmentsInfos,%decoyNb,%targetNb,%sheetNames,%scoreList,%analysisFileNames,%sampleFileNames,%peptideList,%dbRankProt,%assosProtPeptide);
    if($format eq 'pkv'){
        ###> Split the upload file by sheet
        open(FILESTAT,">>$fileStat");
        print FILESTAT "Split the upload file by sheet.\n";
        close FILESTAT;
        #my $childConvert=fork;
        #unless($childConvert){
        #    #>Disconnecting from server
        #    open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
        #    open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
        #    open STDERR, '>/dev/null' or die "Can't open /dev/null: $!";
            system "$promsPath{python}/python $promsPath{python_scripts}/convert_xlsx2txt.py -f $workDir/$file >>$workDir/sortie.txt 2>&1";
        #    exit;
        #}
        my ($lastFile,$wait);
        while (!$wait) {
            sleep 20;
            if ( -s "$workDir/$fileName\_SheetNames.info" && !%sheetNames) {
                open(SHEETNAMES,"<$workDir/$fileName\_SheetNames.info");
                my $i=1;
                foreach (reverse(<SHEETNAMES>)) {
                    my ($sheet,$file)=split(/\t/,$_);
                    chomp $file;
                    $lastFile=$file if $i==1;
                    $sheetNames{'Area-ions'}=$file if $sheet=~/Area - ions/;
                    $sheetNames{'Area-peptides'}=$file if $sheet=~/Area - peptides/;
                    $sheetNames{'Area-proteins'}=$file if $sheet=~/Area - proteins/;
                    $sheetNames{'Score'}=$file if $sheet=~/Score/;
                    $sheetNames{'FDR'}=$file if $sheet=~/FDR/;
                    $sheetNames{'Observed_RT'}=$file if $sheet=~/Observed RT/;
                    $i++;
                }
                close SHEETNAMES;
            }
            $wait=(`tail -1 $workDir/sortie.txt`=~/Done/ && -s $lastFile)? 1:0;
        }
        sleep 3;

        #######################################################################
        ####>Selection of peptide with a fdr < $fdr (selected by the user)<####
        #######################################################################
        ###> In the FDR file
        $nbLine=0;
        open(FILESTAT,">>$fileStat");
        print FILESTAT "Scanning peakview file.\n";
        close FILESTAT;
        my %validFDR;
        my $fdrDec=$fdr/100;
        open(FDRFILE,"<","$sheetNames{'FDR'}") or die ("open: $!");
        while (my $line=<FDRFILE>){
            my @FDRFileLine=split(/\t/,$line);
            next if ($FDRFileLine[0]=~/RT/);
            if ($FDRFileLine[0] eq "Protein") {
                ###>Recovering of samples and analysis's names
                my %sampleNames;
                for (my $i=7;$i<scalar(@FDRFileLine);$i++){
                    $FDRFileLine[$i]=~s/FDR//;
                    $FDRFileLine[$i]=~s/\"//g;
                    my @columnName=split(/\(/,$FDRFileLine[$i]);
                    $analysisFileNames{$i}=$columnName[1];
                    $sampleFileNames{$i}=$columnName[0];
                }
            }
            else{
                for (my $i=7;$i<scalar @FDRFileLine;$i++){  # for each analysis
                    if ($FDRFileLine[$i]<$fdrDec) {
                        $validFDR{"$FDRFileLine[0]_$FDRFileLine[1]_$FDRFileLine[4]"}{$i}=$FDRFileLine[$i];      ## $FDRFileLine[0]= protein &&  $FDRFileLine[1]=peptide seq &&
                    }
                }
            }
        }
        close FDRFILE;


        ###> In the Observed RT file
        $nbLine=0;
        my %validPep;
        open(RTFILE,"<","$sheetNames{'Observed_RT'}") or die ("open: $!");
        while (my $line=<RTFILE>){
            my @RTFileLine=split(/\t/,$line);
            next if ($RTFileLine[0]=~/Protein/);
            my $decoy=$RTFileLine[6];
            for (my $i=7;$i<scalar @RTFileLine;$i++){  # for each analysis
                my $peptide=$RTFileLine[1].'_'.$RTFileLine[4];
                if (exists $validFDR{$RTFileLine[0].'_'.$peptide}{$i}) {
                    my $RT=$RTFileLine[$i];
                    chomp $RT;
                    if ((lc($decoy) eq 'false' || lc($decoy) eq 'faux')){   # $FDRFileLine[6]=DECOY (TRUE or FALSE)
                        ###>List of valid peptides
                        $validPep{'TARGET'}{$peptide}{$i}=$RT;
                        $targetNb{$i}++; #top target peptide number
                    }
                    elsif(lc($decoy) eq 'true' || lc($decoy) eq 'vrai'){
                        $validPep{'DECOY'}{$peptide}{$i}=$RT;
                        $decoyNb{$i}++;  #top decoy peptide number
                    }
                }
            }
        }
        close RTFILE;


        $dbh=&promsConfig::dbConnect('no_user');
        my $sthUpdateMod=$dbh->prepare("UPDATE MODIFICATION SET PEAKVIEW_CODE=? WHERE ID_MODIFICATION=?");
        ###> In the SCORE file
        open(SCOREFILE,"<","$sheetNames{'Score'}") or die ("open: $!");
        $nbLine=0;
        while (my $line=<SCOREFILE>){
            my @ScoreFileLine=split(/\t/,$line);
            next if ($ScoreFileLine[0] eq "Protein");
            ###>Peptide sequence
            my $peptide=$ScoreFileLine[1].'_'.$ScoreFileLine[4];

            ###>Protein identifier
            my @protIdentifiers=split("/",$ScoreFileLine[0]);
            shift @protIdentifiers;
            for (my $i=7;$i<scalar @ScoreFileLine;$i++){
                chomp $ScoreFileLine[$i];
                my $pepScore=$ScoreFileLine[$i];
                if (exists $validPep{'TARGET'}{$peptide}{$i}) {      #$ScoreFileLine[1]=peptide's sequence  and  $ScoreFileLine[4]=charge

                    ###>Recovering peptide RT and score
                    my $pepRT=$validPep{'TARGET'}{$peptide}{$i};
                    my $pepMZ=$ScoreFileLine[3];

                    ###> fetching peptide informations
                    @{$peptideList{$i}{$peptide}{$pepRT}}=($pepMZ,$pepScore);


                    ###> fetching modifications
                    if($ScoreFileLine[1]=~/\[\w+\]/){       ##pepseq
                        my $modifSeq=$ScoreFileLine[1];
                        while ($modifSeq=~/\[(\w+)\]/g){
                            my $peakviewCode=$1;
                            my $modificationID;
                            if($peakview2modID{$peakviewCode}){
                                $modificationID=$peakview2modID{$peakviewCode};
                            }
                            else{
                                open(MODFILE,"<","$promsPath{data}/Unified_Modification_PeakView.csv") or die ("open: $!"); #find peakview code in Unified_Modification_PeakView file
                                my $unimodID;
                                while (my $line=<MODFILE>) {
                                    my @lineStrg=split(/;/,$line);
                                    if ($lineStrg[6] eq $peakviewCode) {
                                        $unimodID=$lineStrg[7];
                                        last;
                                    }
                                }
                                close MODFILE;

                                $modificationID=$unimod2ID{$unimodID};

                                $sthUpdateMod->execute($peakviewCode,$modificationID);
                                $dbh->commit;
                            }
                            my $pos=$-[0];
                            $peptideModif{$i}{$peptide}{$modificationID}{$pos}=$massMods{$modificationID};
                            $modifSeq=~s/\[$1\]//;
                        }
                    }


                    ###>Recover analysis, proteins and peptides attribution
                    foreach my $prot (@protIdentifiers){
                        next if $prot=~/DECOY/;
                        @{$assosProtPeptide{$i}{$prot}{$peptide}}=($pepRT,$pepScore);   #RT,score
                    }
                    push @{$scoreList{$i}{'TARGET'}},$pepScore;        # list of score for TARGET peptide
                }
                elsif(exists $validPep{'DECOY'}{$peptide}{$i}){    # list of score for DECOY peptide
                    push @{$scoreList{$i}{'DECOY'}},$pepScore;
                }
            }
        }
        close SCOREFILE;
        $sthUpdateMod->finish;
        $dbh->disconnect;
    }
    elsif($format eq 'prm'){
        open(FILESTAT,">>$fileStat");
        print FILESTAT "Scanning prm file.\n";
        close FILESTAT;

        my (%colName2Index);
        open(IN,"<$workDir/$file");
        while(<IN>){
            if($.==1){
                $_=~s/\s*\Z//;
                my @colNames=split(/[;,\t]/,$_);
                my @refColNames=('Replicate Name','File Name','Precursor Mz','Precursor Charge','Product Mz','Product Neutral Mass','Product Charge','Fragment Ion','Retention Time','Area','Peptide Modified Sequence','Missed Cleavages','Begin Pos','End Pos','Previous Aa','Next Aa','Peptide Retention Time');
                #check columns names
                foreach my $i (0 .. $#colNames){
                    $colName2Index{$colNames[$i]}=$i;
                }
                foreach my $colName (@refColNames){
                    unless ($colName2Index{$colName}){
                        open(FILEERROR,">>$fileError");
                        print FILEERROR "Wrong columns names in file ! The columns names have to be like : <BR>Protein Name<BR>Replicate Name<BR>File Name<BR>Precursor Mz<BR>Precursor Charge<BR>Product Mz<BR>Product Neutral Mass<BR>Product Charge<BR>Fragment Ion<BR>Retention Time<BR>Area<BR>Peptide Modified Sequence<BR>Missed Cleavages<BR>Begin Pos<BR>End Pos<BR>Previous Aa<BR>Next Aa<BR>Peptide Retention Time";
                        close(FILEERROR);
                        exit;
                    }
                }
            }else{
                $_=~s/\s*\Z//;
                my @infos=split(/[;,\t]/,$_);
                my $prot=$infos[$colName2Index{'Protein Name'}];
                next unless $prot;
                my $anaName=$infos[$colName2Index{'Replicate Name'}];
                my $anaFile=$infos[$colName2Index{'File Name'}];
                my $pepMZ=$infos[$colName2Index{'Precursor Mz'}];
                my $pepCharge=$infos[$colName2Index{'Precursor Charge'}];
                my $fragmentMZ=$infos[$colName2Index{'Product Mz'}];
                my $fragmentNeutralMass=$infos[$colName2Index{'Product Neutral Mass'}];
                my $fragmentCharge=$infos[$colName2Index{'Product Charge'}];
                my $fragmentType=$infos[$colName2Index{'Fragment Ion'}];
                my $fragRT=$infos[$colName2Index{'Retention Time'}];
                my $area=$infos[$colName2Index{'Area'}];
                my $modifSeq=$infos[$colName2Index{'Peptide Modified Sequence'}];
                my $missCleav=$infos[$colName2Index{'Missed Cleavages'}];
                my $beg=$infos[$colName2Index{'Begin Pos'}];
                $beg=($beg eq '#N/A') ? '' : $beg;
                my $end=$infos[$colName2Index{'End Pos'}];
                $end=($end eq '#N/A') ? '' : $end;
                my $begFlank=$infos[$colName2Index{'Previous Aa'}];
                $begFlank=($begFlank eq 'X') ? '' : $begFlank;
                my $endFlank=$infos[$colName2Index{'Next Aa'}];
                $endFlank=($endFlank eq 'X') ? '' : $endFlank;
                my $pepRT=$infos[$colName2Index{'Peptide Retention Time'}];
                $analysisFileNames{$anaName}=$anaFile;
                $sampleFileNames{$anaName}=$anaName;

                my $peptide=$modifSeq.'_'.$pepCharge;
                my $fragment=$fragmentType.'_'.$fragmentCharge;

                ##peptide position
                $peptidePosition{$anaName}{$prot}{$peptide}{"$beg\_$begFlank"}="$end\_$endFlank" if ($beg && $begFlank && $end && $endFlank);


                $dbRankProt{$prot}=1 unless defined $dbRankProt{$prot};

                my $modSeq=$modifSeq;
                while($modSeq=~/(\w*)\[\+(\d*)\]/g){
                    my $modifID;
                    foreach my $modID (keys %modificationList){
                        $modifID=$modID if sprintf("%.0f",$modificationList{$modID})==$2;
                        last if $modifID;
                    }
                    my @aa=split("",$1);
                    my $aaPrev=($aa[-1])? $aa[-1] : '=';
                    $modifIDs{$modifID}{'V'}{$aaPrev}=1 if $modifID;
                    my $pos=length($1);
                    $modSeq=~s/\[\+(\d*)\]//;
                    $peptideModif{$anaName}{$peptide}{$modifID}{$pos}=$modificationList{$modifID};
                }

                ###>fetching peptide informations
                $peptideInfo{$peptide}[1]=$missCleav;
                @{$peptideList{$anaName}{$peptide}{"$pepRT"}}=($pepMZ,'');


                ###>fetching fragment informations
                @{$fragmentsInfos{$anaName}{$peptide}{$pepRT}{$fragment}}=($fragRT,$area,$fragmentMZ,$fragmentNeutralMass);


                ###>List of protein
                $peptideProt{$anaName}{$peptide}{$prot}=1;
                @{$assosProtPeptide{$anaName}{$prot}{$peptide}}=($pepRT);
            }
        }
        close IN;

        $quantifAnnot="LABEL=FREE::SOFTWARE=SKY;".$softwareVersion;
    }
    elsif ($format eq 'openswath'){
        open(FILESTAT,">>$fileStat");
        print FILESTAT "Scanning openswath file.\n";
        close FILESTAT;
        my %colName2Index;
        open (IN, "<$workDir/$file");
        $nbLine=0;
        my $lastPeptide='';
        my $lastAna;
        while (<IN>){
            if ($.==1){
                $_=~s/\s*\Z//;
                my @columns=split(/[,;\t]/,$_);
                my @refColNames=('ProteinName','filename','m/z','Charge','aggr_Fragment_Annotation','RT','aggr_Peak_Area','FullPeptideName');
                foreach my $i (0 .. $#columns){
                    $colName2Index{$columns[$i]}=$i;
                }

                ##check columns names
                foreach my $colName (@refColNames){
                    unless ($colName2Index{$colName}){
                        open(FILEERROR,">>$fileError");
                        print FILEERROR "Wrong columns names in file ! <BR>The columns names have to be like : <BR>ProteinName<BR>filename<BR>m/z<BR>Charge<BR>aggr_Fragment_Annotation<BR>RT<BR>aggr_Peak_Area<BR>FullPeptideName";
                        close(FILEERROR);
                        exit;
                    }
                }
            }
            else{
                $_=~s/\s*\Z//;
                my @infos=split(/[\t]/,$_);
                next if $infos[$colName2Index{'decoy'}]==1;
                next if $infos[$colName2Index{'ProteinName'}]=~/iRT/;
                next if $infos[$colName2Index{'ProteinName'}]=~/^1\/reverse_sp/;
                my $proteinsList=$infos[$colName2Index{'ProteinName'}];
                my $anaFile=$infos[$colName2Index{'filename'}];
                my $pepMZ=$infos[$colName2Index{'m/z'}];
                my $pepCharge=$infos[$colName2Index{'Charge'}];
                my $RT=$infos[$colName2Index{'RT'}];
                $RT=($RT/60)*1;
                my $areasList=$infos[$colName2Index{'aggr_Peak_Area'}];
                $areasList=~s/"//g;
                my $fragmentsList=$infos[$colName2Index{'aggr_Fragment_Annotation'}];
                $fragmentsList=~s/"//g;
                my $modifSeq=$infos[$colName2Index{'FullPeptideName'}];
                my $mScore=$infos[$colName2Index{'m_score'}];
                undef @infos;

                (my $anaName=$anaFile)=~s/.mzXML//;
                $anaName=fileparse($anaName, qr/\.[^.]*/);
                $anaFile=fileparse($anaFile, qr/\.[^.]*/);
                $analysisFileNames{$anaName}=$anaFile;
                $sampleFileNames{$anaName}=$anaName;


                my $peptide=$modifSeq.'_'.$pepCharge;

                ###>recovering peptide infos
                @{$peptideList{$anaName}{$peptide}{"$RT"}}=($pepMZ,'',$fragmentsList,$areasList);


                ##>recovering peptide/protein associations
                if($peptide eq $lastPeptide){
                    $peptideProt{$anaName}{$peptide}=\%{$peptideProt{$lastAna}{$peptide}};
                }
                else{
                    my @proteins=split(/\//,$proteinsList);
                    for my $i (1 .. $#proteins){
                        next if $proteins[$i]=~/reverse/;
                        @{$assosProtPeptide{$anaName}{$proteins[$i]}{$peptide}}=($RT);
                        $peptideProt{$anaName}{$peptide}{$proteins[$i]}=1;
                    }
                    undef @proteins;
                }


                ###>recovering peptide modifications
                if ($modifSeq=~/\(UniMod:\d+\)/){
                    my $peptideSeq=$modifSeq;
                    while ($peptideSeq=~/\(UniMod:(\d+)\)/g){
                        my $modificationID=$unimod2ID{$1};
                        my $pos=$-[0];
                        $peptideModif{$anaName}{$peptide}{$modificationID}{$pos}=$massMods{$modificationID};
                        $peptideSeq=~s/\(UniMod:$1\)//;
                    }
                    undef $peptideSeq;
                }

                $lastPeptide=$peptide;
                $lastAna=$anaName;
                undef $peptide;
                undef $anaName;
            }
        }
        close IN;
    }
    elsif ($format eq 'spectronaut'){
        open(FILESTAT,">>$fileStat");
        print FILESTAT "Scanning spectronaut file.\n";
        close FILESTAT;

        my %colName2Index;
        open (IN, "<$workDir/$file");
        $nbLine=0;
        my $lastPeptide='';
        my $lastAna;

        while (<IN>){
            if ($.==1){
                $_=~s/\s*\Z//;
                my @columns=split(/[,;\t\s]/,$_);
                my @refColNames=('R.Label','PG.UniProtIds','EG.ModifiedPeptide','EG.RTPredicted','FG.Charge','FG.PrecMz','F.Charge','F.FrgIon','F.FrgMz','F.PeakArea');
                foreach my $i (0 .. $#columns){
                    $colName2Index{$columns[$i]}=$i;
                }

                ##check columns names
                foreach my $colName (@refColNames){
                    unless (defined $colName2Index{$colName}){
                        open(FILEERROR,">>$fileError");
                        print FILEERROR "Wrong columns names in file ! <BR>The columns names have to be like : <BR>R.Label<BR>PG.UniProtIds<BR>EG.ModifiedPeptide<BR>EG.RTPredicted<BR>FG.Charge<BR>FG.PrecMz<BR>F.Charge<BR>F.FrgIon<BR>F.FrgMz<BR>F.PeakArea";
                        close(FILEERROR);
                        exit;
                    }
                }
            }else{
                $_=~s/\s*\Z//;
                my @infos=split(/[\t]/,$_);
                next if $infos[$colName2Index{'PG.UniProtIds'}]=~/iRT/;
                next if $infos[$colName2Index{'PG.UniProtIds'}]=~/^1\/reverse_sp/;
                #if ($infos[$colName2Index{'EG.ModifiedPeptide'}]!~/^_/){
                #    print $infos[$colName2Index{'EG.ModifiedPeptide'}],"<BR>";
                #}
                next if $infos[$colName2Index{'EG.ModifiedPeptide'}]!~/^_/;
                my $proteinsList=$infos[$colName2Index{'PG.UniProtIds'}];
                my $anaFile=$infos[$colName2Index{'R.Label'}];
                my $pepMZ=$infos[$colName2Index{'FG.PrecMz'}];
                $pepMZ=~s/,/\./;
                my $pepCharge=$infos[$colName2Index{'FG.Charge'}];
                my $RT=$infos[$colName2Index{'EG.RTPredicted'}];
                $RT=~s/,/\./;
                $RT=sprintf("%0.2f",$RT);
                my $area=$infos[$colName2Index{'F.PeakArea'}];
                if ($area!~/\d+\,?\d*/){
                    $area='NA';
                }
                else{
                    $area=~s/,/\./;
                }

                my $fragmentType=$infos[$colName2Index{'F.FrgIon'}];
                my $fragmentCharge=$infos[$colName2Index{'F.Charge'}];
                my $fragment=$fragmentType.'_'.$fragmentCharge;
                my $fragmentMZ=$infos[$colName2Index{'F.FrgMz'}];
                $fragmentMZ=~s/,/\./;
                my $modifSeq=$infos[$colName2Index{'EG.ModifiedPeptide'}];


                (my $anaName=$anaFile)=~s/\.raw//;
                $analysisFileNames{$anaName}=$anaFile;
                $sampleFileNames{$anaName}=$anaName;


                $modifSeq=~s/_//g;
                my $peptide=$modifSeq.'_'.$pepCharge;

                ###>recovering peptide infos
                @{$peptideList{$anaName}{$peptide}{"$RT"}}=($pepMZ,'');

                ###> recovering fragments infos
                @{$fragmentsInfos{$anaName}{$peptide}{"$RT"}{$fragment}}=('',$area,$fragmentMZ,'');


                ##>recovering peptide/protein associations
                if($peptide eq $lastPeptide){
                    $peptideProt{$anaName}{$peptide}=\%{$peptideProt{$lastAna}{$peptide}};
                }
                else{
                    my @proteins=split(/\//,$proteinsList);
                    for my $i (1 .. $#proteins){
                        next if $proteins[$i]=~/reverse/;
                        @{$assosProtPeptide{$anaName}{$proteins[$i]}{$peptide}}=($RT);
                        $peptideProt{$anaName}{$peptide}{$proteins[$i]}=1;
                    }
                }


                ###>recovering peptide modifications
                if ($modifSeq=~/\[\+(\d+)\]/){
                    my $peptideSeq=$modifSeq;
                    while ($peptideSeq=~/\[\+(\d+)\]/g){
                        my $modificationID;
                        foreach my $massID (keys %massMods){
                            if (sprintf("%0.f",$massMods{$massID}) == $1){
                                $modificationID=$massID;
                                last;
                            }
                        }
                        unless ($modificationID){
                            open(FILEERROR,">>$fileError");
                            print FILEERROR "The associated modification for the peptide $modifSeq is not found.";
                            close(FILEERROR);
                            exit;
                        }
                        my $pos=$-[0];
                        $peptideModif{$anaName}{$peptide}{$modificationID}{$pos}=$massMods{$modificationID};
                        $peptideSeq=~s/\[\+($1)\]//;
                    }
                }

                $lastPeptide=$peptide;
                $lastAna=$anaName;

            }
        }
        close IN;
    }

    ##########################
    ### RT-ordered peptide ###
    ##########################
    my (%peptideRTOrder,%peptideRT);
    foreach my $analysis (keys %peptideList){
        foreach my $peptide (keys %{$peptideList{$analysis}}){
            foreach my $rt (sort {$peptideList{$analysis}{$peptide}{$a} <=> $peptideList{$analysis}{$peptide}{$b}} keys %{$peptideList{$analysis}{$peptide}} ){
                $peptideRT{$peptide}=$rt;
            }
        }
        my $RTorder=1;
        foreach my $peptide (sort {$peptideRT{$a} <=> $peptideRT{$b}} keys %peptideRT){
            @{$peptideRTOrder{$analysis}{$peptide}}=($peptideRT{$peptide},$RTorder);
            $RTorder++;
        }
    }

    my (%numMatchProt,%protCoverage,%sequenceProtList,%protLength,%protDes,%protSeq,%protMW,%protOrg);
    unless ($format eq 'prm' && %peptidePosition){
        ##########################################################
        ####>Fetching protein sequence and mass from dataBank<####
        ##########################################################
        open(FILESTAT,">>$fileStat");
        print FILESTAT "Extracting protein sequences from sequence databank.\n";
        close FILESTAT;

        $dbh=&promsConfig::dbConnect('no_user');# if $format eq 'pkv'; # reconnect after long process
        foreach my $dbInfo (@dbID){
            my ($dbID,$dbRank)=split(/&/,$dbInfo);
            my $dbFasta=$dbh->selectrow_array("SELECT FASTA_FILE FROM DATABANK WHERE ID_DATABANK='$dbID' ");
            my $dbDir="$promsPath{data}/banks/db_$dbID";
            my $j=1;
            ###>Fetching parsing rules
            my ($parseRules,$identType)=$dbh->selectrow_array("SELECT PARSE_RULES,DEF_IDENT_TYPE FROM DATABANK_TYPE DT,DATABANK D WHERE D.ID_DBTYPE=DT.ID_DBTYPE AND ID_DATABANK=$dbID");
            $identType=($identType) ? $identType : '';
            my @rules=split(',:,',$parseRules);
            my ($idRule)=($rules[0]=~/ID=(.+)/);
            my ($desRule,$orgRule);
            if ($rules[1]) {
                if ($rules[1]=~/DES=/) {($desRule=$rules[1])=~s/DES=//;}
                else {($orgRule=$rules[1])=~s/ORG=//;}
            }
            if ($rules[2]) {($orgRule=$rules[2])=~s/ORG=//;}

            my %massATave=&promsConfig::getMassATave; # Average mass, needed for protein mass calculation
            my ($des,$org);
            my $count=0;
            if (-e "$dbDir/$dbFasta"){
                open(FASTAFILE,"<","$dbDir/$dbFasta") or die ("open: $!");
                my ($checkReverse,$identifier,$matched,$sequence,$info);
                my $nbmatch=0;
                while (my $line=<FASTAFILE>){
                    print ".";
                    if ($line=~/>reverse_/) {$checkReverse=1;}  #next if reverse = true
                    elsif ($line=~m/^>(.+)\n/) {       # prot description line
                        (my $newLine=$line)=~s/^>\s*//;
                        if ($identifier) {      # search prot
                            $sequence=~s/\s//;
                            $dbRankProt{$identifier}=$dbRank unless defined $dbRankProt{$identifier};
                            $sequenceProtList{$identifier}=$sequence; # {identifier}=sequence
                            $protLength{$identifier}=length($sequence);
                            $protDes{$identifier}=($des)? $des : $info;
                            $protOrg{$identifier}=($org)? $org : '';

                            my %countAA;
                            foreach my $aa (split(//,uc($sequence))) {
                                $countAA{$aa}++;
                            }
                            my $mass=0;
                            foreach my $aa (keys %countAA) {
                                $mass+=($massAAave{$aa}*$countAA{$aa}) if $massAAave{$aa}; # some characters are not amino-acids
                            }
                            $mass+=($massATave{H}+$massATave{H}+$massATave{O}) if $mass; # H(Nter) + OH(Cter)
                            $mass=sprintf "%.2f",$mass;
                            $mass=$mass*1;
                            $protMW{$identifier}=$mass;

                            $des='';
                            $org='';
                            $sequence='';
                        }
                        $checkReverse=0;
                        my @line=split(//,$newLine);
                        foreach my $entry (@line){
                            ($identifier)=($entry=~/$idRule/);
                            chomp $identifier;
                            if ($identType eq 'UNIPROT_ID') { # check for isoforms in 1st keyword before | & add to UNIPROT_ID
                                $entry=~s/^sp\|//;
                                if ($entry=~/^[^|]+-(\d+)\|/) {
                                    $identifier.="-$1";
                                }
                            }
                            #print $identifier,"<BR>";
                            ($des)=($entry=~/$desRule/) if $desRule; # optional
                            ($org)=($entry=~/$orgRule/) if $orgRule; # optional
                            if (not $des && not $org) {
                                $info=$newLine;
                            }
                        }
                    }
                    else{
                        if ($checkReverse==0) {
                            chomp $line;
                            $line=~s/\s//g;
                            $sequence.=$line;       # protein's sequence
                        }
                    }
                }
                 ## last prot
                $sequenceProtList{$identifier}=$sequence;
                $protLength{$identifier}=length($sequence);
                my %countAA;
                foreach my $aa (split(//,uc($sequence))) {
                    $countAA{$aa}++;
                }
                my $mass=0;
                foreach my $aa (keys %countAA) {
                    $mass+=($massAAave{$aa}*$countAA{$aa}) if $massAAave{$aa}; # some characters are not amino-acids
                }
                $mass+=($massATave{H}+$massATave{H}+$massATave{O}) if $mass; # H(Nter) + OH(Cter)
                $mass=sprintf "%.2f",$mass;
                $mass=$mass*1;
                $protMW{$identifier}=$mass;
                close FASTAFILE;
            }
            else {
                open(FILEERROR,">>$fileError");
                print FILEERROR "Fasta file not found! Skipping.";
                close(FILEERROR);
                exit;
            }
        }


        #####################################################################
        ####>Recovering peptide positions on protein and coverage calcul<####
        #####################################################################
        open(FILESTAT,">>$fileStat");
        print FILESTAT "Recovering peptide's positions on protein.\n";
        close FILESTAT;
        my $count=0;
        my %test;
        foreach my $analysisNumber (keys %assosProtPeptide){        # for each analysis
            foreach my $protein (keys $assosProtPeptide{$analysisNumber}){
                print ".";
                my %positionPeptide;
                my $sequence=$sequenceProtList{$protein};
                my $numMatch=0;
                my (@begPep,@endPep);
                foreach my $peptide (keys $assosProtPeptide{$analysisNumber}{$protein}){
                    my @pepSeq=split(/_/,$peptide);
                    my $pepSequence=$pepSeq[0];
                    $pepSequence=~s/\[\+?\w+\]//g if ($format eq 'pkv' || $format eq 'spectronaut');
                    $pepSequence=~s/\(UniMod:\d+\)//g if ($format eq 'openswath');

                    if ($sequence=~m/$pepSequence/ ) {  #peptide matches on protein
                        while ($sequence=~m/(.?)$pepSequence(.?)/g){
                            my $firstPosition=($1)? length($`)+2 :length($`)+1; #+2 : 1 aa before the peptide is matched
                            $positionPeptide{$firstPosition}++;
                            my $lastPosition=$firstPosition+length($pepSequence)-1;
                            $positionPeptide{$lastPosition}--;
                            my $begFlank=$1;
                            my $endFlank=$2;
                            if ($begFlank eq '') {   #N-term
                                $begFlank="-";
                            }
                            if ($endFlank eq '') {  #C-term
                                $endFlank="-";
                            }
                            $peptidePosition{$analysisNumber}{$protein}{$peptide}{"$firstPosition\_$begFlank"}="$lastPosition\_$endFlank";
                            $numMatch++;
                        }
                    }
                    else{   #for Leucine and Isoleucine inversion in peptide sequence
                        delete $assosProtPeptide{$analysisNumber}{$protein}{$peptide};
                    }
                }
                $numMatchProt{$analysisNumber}{$protein}=$numMatch; #number of peptide match
                delete $assosProtPeptide{$analysisNumber}{$protein} if scalar keys %{$assosProtPeptide{$analysisNumber}{$protein}} == 0;

                if ($assosProtPeptide{$analysisNumber}{$protein}) {
                    ###> Coverage calcul
                    my $pepCoverage;
                    my $hasPeptide=0;
                    my $matchBeg=0;
                    foreach my $position (sort {$a<=>$b} keys %positionPeptide) {
                        if ($hasPeptide==0) { ## start of peptide region
                            $matchBeg=$position;
                        }
                        $hasPeptide+=$positionPeptide{$position};
                        if ($hasPeptide==0 ) { ## end of peptide region
                            $pepCoverage+=($position-$matchBeg+1);
                        }
                    }
                    my $peptideCoverage=sprintf ("%.1f",$pepCoverage/length($sequence)*100);
                    $protCoverage{$analysisNumber}{$protein}=$peptideCoverage;
                }
            }
        }
    }


    ############################################################################################################
    ###Inserting data into tables SAMPLE,ANALYSIS,ANALYSIS_SWATH_LIB,ANALYSIS_DATABANK,ANALYSIS_MODIFICATION ###
    ############################################################################################################
    my $count=0;
    my %proteinList;
    open(FILESTAT,">>$fileStat");
    print FILESTAT "Storing analysis and samples data into database.\n";
    close FILESTAT;
    my $sthSampleID=$dbh->prepare("SELECT MAX(ID_SAMPLE) FROM SAMPLE");
    my $sthSample=$dbh->prepare("INSERT INTO SAMPLE (ID_SAMPLE,ID_EXPERIMENT,NAME,START_DATE,DISPLAY_POS,UPDATE_DATE,UPDATE_USER) VALUES (?,?,?,NOW(),?,NOW(),?)");
    my $sthAnalysis=$dbh->prepare("INSERT INTO ANALYSIS (ID_SAMPLE,NAME,START_DATE,DATA_FILE,VALID_STATUS,VALID_USER,LOWER_SCORES,FIRST_VALID_DATE,VALID_DATE,UPDATE_DATE,UPDATE_USER,WIFF_FILE,DECOY,FALSE_DISCOVERY,TAXONOMY,DISPLAY_POS,FILE_FORMAT,MS_TYPE,INSTRUMENT,VERIFIED_MG,MIN_SCORE,MAX_RANK) VALUES (?,?,NOW(),?,?,?,?,NOW(),NOW(),NOW(),?,?,?,?,?,?,?,?,?,?,?,?)");
    my $sthAnaSwath=$dbh->prepare("INSERT INTO ANALYSIS_SWATH_LIB (ID_ANALYSIS,ID_SWATH_LIB,VERSION_NAME) VALUES (?,?,?)");
    my $sthAnaDB=$dbh->prepare("INSERT INTO ANALYSIS_DATABANK (ID_DATABANK,ID_ANALYSIS,DB_RANK) VALUES (?,?,?)");
    my $sthAnaMod=$dbh->prepare("INSERT INTO ANALYSIS_MODIFICATION (ID_ANALYSIS,ID_MODIFICATION,MODIF_TYPE,SPECIFICITY) VALUES (?,?,?,?)");
    my (%analysisID);
    my $displayPos=$dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM SAMPLE WHERE ID_EXPERIMENT=$experimentID");
    foreach my $analysis (sort {$a cmp $b} keys %analysisFileNames){
        my $analysisFile=$analysisFileNames{$analysis};
        my $sampleName=$sampleFileNames{$analysis};
        my $analysisName=($format eq 'pkv') ? fileparse($analysisFile, qr/\.[^.]*/) : $analysis;

        $sthSampleID->execute;
        my $sampleID=$sthSampleID->fetchrow_array;

        $sampleID+=1;
        my $trueNb=($decoyGhostNb{$analysis}) ? $decoyGhostNb{$analysis} : 0;
        my $falseNb=($targetGhostNb{$analysis}) ? $targetGhostNb{$analysis} : 0;
        my $topPepNbDecoy=($decoyNb{$analysis}) ? $decoyNb{$analysis} : 0;
        my $topPepNbTarget=($targetNb{$analysis}) ? $targetNb{$analysis} : 0;

        my $fileFormat=($format eq 'pkv')? 'SWATH.PKV': ($format eq 'prm')?  'SKYLINE.CSV' : ($format eq 'openswath')? 'OPENSWATH.TSV' : 'SPECTRONAUT.XLS';
        my $msType=($format eq 'pkv')? 'SWATH': ($format eq 'prm')? 'TDA' : 'DIA';
        my $decoy=($format eq 'pkv') ? "INT:SEARCH,FDR=$fdr:precomputed" : '';

        $displayPos++;
        ###Insertion into SAMPLE
        $sthSample->execute($sampleID,$experimentID,$sampleName,$displayPos,$userID);
        $sthSample->finish;

        ###Insertion into ANALYSIS
        $sthAnalysis->execute($sampleID,$analysisName,$file,-1,$userID,0,$userID,$analysisFile,$decoy,"$topPepNbDecoy:$trueNb:$topPepNbTarget:$falseNb",$taxonomy,$displayPos,$fileFormat,$msType,$instrument,0,0,1);
        $sthAnalysis->finish;
        my $analysisID=$dbh->last_insert_id(undef, undef, 'ANALYSIS', 'ID_ANALYSIS');

        ###Insertion into ANALYSIS_SWATH_LIB
        unless ($format eq 'prm'){
            $sthAnaSwath->execute($analysisID,$libID,$versionName);
            $sthAnaSwath->finish;
        }

        ###Insertion into ANALYSIS_DATABANK
        foreach my $dbInfo (@dbID){
            my ($dbID,$dbRank)=split(/&/,$dbInfo);
            $sthAnaDB->execute($dbID,$analysisID,$dbRank);
        }
        $sthAnaDB->finish;


        foreach my $modifID (keys %modifIDs){
            foreach my $modType (keys %{$modifIDs{$modifID}}){
                my @aaList;
                foreach my $aa (keys %{$modifIDs{$modifID}{$modType}}){
                    push @aaList,$aa;
                }
                my $specif=join(",",@aaList);
                $sthAnaMod->execute($analysisID,$modifID,$modType,$specif);
                $sthAnaMod->finish;
            }
        }
        $analysisID{$analysis}=$analysisID;

        foreach my $protein (keys %{$assosProtPeptide{$analysis}}){
            $proteinList{$protein}=1;
        }
    }
    $sthSampleID->finish;


    ### Recovering protein informations
    if($format eq 'prm' || $format eq 'openswath'){
        my ($dbID,$dbRank)=split(/&/,$dbID[0]);
        my @analysisList=values (%analysisID);
        my %protList;
        foreach my $analysis (keys %assosProtPeptide){
            foreach my $protein (keys %{$assosProtPeptide{$analysis}}){
                $protList{$protein}{$analysisID{$analysis}}=1;
            }
        }
        &promsMod::getProtInfo('silent',$dbh,$dbID,\@analysisList,\%protDes,\%protMW,\%protOrg,\%protLength,\%protSeq,\%protList);


        ### compute peptide specificity and protein coverage
        foreach my $anaName (keys %peptideProt){
            foreach my $pep (keys %{$peptideProt{$anaName}}){
                if(scalar keys %{$peptideProt{$anaName}{$pep}} == 1){
                    $peptideInfo{$pep}[2]=100;
                }else{
                    $peptideInfo{$pep}[2]=(1/scalar keys %{$peptideProt{$anaName}{$pep}})*100;
                }
            }

            ## protein coverage
            foreach my $protein (keys %{$peptidePosition{$anaName}}){
                next unless $protLength{$protein};
                my %positionPeptide;
                foreach my $peptide (keys %{$peptidePosition{$anaName}{$protein}}){
                    foreach my $infoBeg (keys %{$peptidePosition{$anaName}{$protein}{$peptide}}){
                        my ($beg,$begFlank)=split('_',$infoBeg);
                        my ($end,$endFlank)=split('_',$peptidePosition{$anaName}{$protein}{$peptide}{$infoBeg});
                        $positionPeptide{$beg}++;
                        $positionPeptide{$end}--;
                        $numMatchProt{$anaName}{$protein}++;    #number of peptide match
                    }
                }
                my $hasPep=0;
                my ($pepCoverage,$matchBeg);
                foreach my $position (sort {$a <=> $b} keys %positionPeptide){
                    if($hasPep == 0){
                        $matchBeg=$position;
                    }
                    $hasPep+=$positionPeptide{$position};
                    if($hasPep==0){
                        $pepCoverage+=($position-$matchBeg+1);
                    }
                }
                my $coverage=sprintf ("%.1f",$pepCoverage/$protLength{$protein}*100);
                $protCoverage{$anaName}{$protein}=$coverage;
            }
        }
    }
    $dbh->commit;

    ###########################################
    #####Inserting data into table PROTEIN#####
    ###########################################
    open(FILESTAT,">>$fileStat");
    print FILESTAT "Storing proteins data into database.\n";
    close FILESTAT;
    my %proteinID;
    my $sthProtein=$dbh->prepare("INSERT INTO PROTEIN (ID_PROJECT,ID_MASTER_PROTEIN,ALIAS,IDENTIFIER,PROT_DES,PROT_LENGTH,PROT_SEQ,MW,COMMENTS,UPDATE_DATE,UPDATE_USER,ORGANISM) VALUES (?,NULL,?,?,?,?,?,?,NULL,NULL,NULL,?)");
    $count=0;
    foreach my $protein (keys %proteinList){
        print ".";
        ##recovering the alias of the protein in terms of the project identifier
        my $alias;
        if (! $identProject) {$alias=$protein;}
        else{
            if ($protein=~/[sp|tr]\|(\w+)\|(\w+_\w+)/) { #uniprot ALL
                if ($identProject==1) {$alias=$1;}
                elsif($identProject==2){$alias=$2;}
            }
            else{$alias=$protein;}
        }


        ### if this protein exists in PROTEIN table
        if(exists $projectProt{$protein} ){
            my $sthProteinUpdate;
            $proteinID{$protein}=$projectProt{$protein};
            # protein length = null (PROTEIN table)
            if (exists $incompleteProtein{$protein}){
                $sthProteinUpdate=$dbh->prepare("UPDATE PROTEIN SET PROT_LENGTH=? WHERE ID_PROTEIN=?");
                $sthProteinUpdate->execute($protLength{$protein},$projectProt{$protein});
            }
            # no protein description (PROTEIN table)
            if (exists $noDescriptionProtein{$protein}){
                $sthProteinUpdate=$dbh->prepare("UPDATE PROTEIN SET PROT_DES=? WHERE ID_PROTEIN=?");
                $sthProteinUpdate->execute($protDes{$protein},$projectProt{$protein});
            }
            # no protein sequence (PROTEIN table)
            if (exists $noSequenceProtein{$protein}){
                $sthProteinUpdate=$dbh->prepare("UPDATE PROTEIN SET PROT_SEQ=? WHERE ID_PROTEIN=?");
                $sthProteinUpdate->execute($sequenceProtList{$protein},$projectProt{$protein});
            }
            # no protein MW (PROTEIN table)
            if (exists $noMWProtein{$protein}){
                $sthProteinUpdate=$dbh->prepare("UPDATE PROTEIN SET MW=? WHERE ID_PROTEIN=?");
                $sthProteinUpdate->execute($protMW{$protein},$projectProt{$protein});
            }
        }
        else {   #no duplicates !!
            my $specie=($protOrg{$protein})? $protOrg{$protein} : $taxonomy;
            $sthProtein->execute($projectID,$alias,$protein,$protDes{$protein},$protLength{$protein},$sequenceProtList{$protein},$protMW{$protein},$specie) or die $dbh->errstr;
            my $protID=$dbh->last_insert_id(undef, undef, 'PROTEIN', 'ID_PROTEIN');
            $proteinID{$protein}=$protID;
        }
    }
    $sthProtein->finish;
    $dbh->commit;



    ##################################################
    ####Inserting data into table ANALYSIS_PROTEIN ###
    ##################################################
    my $sthTopProt=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=? AND VISIBILITY=2");
    my $sthBestVis=$dbh->prepare("SELECT MAX(VISIBILITY) FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=?");
    my $sthAnaProt=$dbh->prepare("INSERT INTO ANALYSIS_PROTEIN (ID_ANALYSIS,ID_PROTEIN,DB_RANK,CONF_LEVEL,SCORE,NUM_PEP,NUM_MATCH,PEP_COVERAGE,MATCH_GROUP,PEP_SPECIFICITY,VISIBILITY) VALUES (?,?,?,?,?,?,?,?,?,?,?)");
    ###>Fetching protein score, pep number and protein specificity
    foreach my $analysis (sort {$a cmp $b } keys %analysisID){
        my $analysisName=($format eq 'pkv') ? fileparse($analysisFileNames{$analysis}, qr/\.[^.]*/) : $analysis;
        my $analysisID=$analysisID{$analysis};
        my (%matchList,%proteinScore,%bestPepSpecificity,@pepInfo);
        $count=0;
        my $protScore;
        foreach my $protein (keys $assosProtPeptide{$analysis}){
            $protScore=0;
            my $pepSpecificity=0;
            foreach my $peptide (keys $assosProtPeptide{$analysis}{$protein}){
                $matchList{$protein}{$peptide}=1;
                ###protein score
                $protScore+=$assosProtPeptide{$analysis}{$protein}{$peptide}[1] if $assosProtPeptide{$analysis}{$protein}{$peptide}[1];

                ###max peptide specificity
                $pepSpecificity=($peptideInfo{$peptide}[2] && $pepSpecificity<$peptideInfo{$peptide}[2]) ? $peptideInfo{$peptide}[2] : $pepSpecificity;
            }
            $bestPepSpecificity{$protein}=$pepSpecificity;
            $proteinScore{$protein}=$protScore if $protScore;
        }

        ###>Finding number of times proteins are at top of match group hierarchy
        open(FILESTAT,">>$fileStat");
        print FILESTAT "Analysis $analysisName : Computing proteins visibility.\n";
        close FILESTAT;
        $count=0;
        my (%numProtTop,%bestProtVis,%numProtPeptides,%proteinPepDens);
        foreach my $protein (keys %matchList) {
            my $protID=$proteinID{$protein};
            if ($projectProt{$protein}) { # proteins are already in project
                $sthTopProt->execute($protID);
                ($numProtTop{$protein})=$sthTopProt->fetchrow_array;
                if ($protVisibility) {
                    $sthBestVis->execute($protID);
                    $bestProtVis{$protein}=$sthBestVis->fetchrow_array;
                }
            }
            else { # proteins are new in project
                $numProtTop{$protein}=0;
                $projectProt{$protein}=$protID; # add protein to project list
            }
            #>Peptide density
            $proteinPepDens{$protein}=$protLength{$protein}/(scalar (keys %{$assosProtPeptide{$analysis}{$protein}}));

        }
        $sthTopProt->finish;
        $sthBestVis->finish;


        #####>Finding match groups and protein visibility
        open(FILESTAT,">>$fileStat");
        print FILESTAT "Analysis $analysisName : Building match groups.\n";
        close FILESTAT;
        my (%visibility,%matchGroup);
        $count=0;
        #my @sortedProt=sort{$numProtPeptides{$b}<=>$numProtPeptides{$a} || $numProtTop{$b}<=>$numProtTop{$a} || $proteinScore{$b}<=>$proteinScore{$a} || &deltaLength($protLength{$a},$protLength{$b},$proteinPepDens{$a},$proteinPepDens{$b})  || $protLength{$a}<=>$protLength{$b} } keys %matchList;
        my @sortedProt;
        if(%proteinScore){
            @sortedProt=sort{scalar (keys %{$assosProtPeptide{$analysis}{$b}})<=>scalar (keys %{$assosProtPeptide{$analysis}{$a}}) || $numProtTop{$b}<=>$numProtTop{$a} || $proteinScore{$b}<=>$proteinScore{$a} || &deltaLength($protLength{$a},$protLength{$b},$proteinPepDens{$a},$proteinPepDens{$b})  || $protLength{$a}<=>$protLength{$b} } keys %matchList;
        }else{
            @sortedProt=sort{scalar (keys %{$assosProtPeptide{$analysis}{$b}})<=>scalar (keys %{$assosProtPeptide{$analysis}{$a}}) || $numProtTop{$b}<=>$numProtTop{$a} || &deltaLength($protLength{$a},$protLength{$b},$proteinPepDens{$a},$proteinPepDens{$b})  || $protLength{$a}<=>$protLength{$b} } keys %matchList;
        }
        my $numGroup=0;
        foreach my $i (0..$#sortedProt){
            my $prot1=$sortedProt[$i];
            next if $matchGroup{$prot1}; # already assigned to a match group
            $matchGroup{$prot1}=++$numGroup;
            $visibility{$prot1}=2;
            foreach my $j ($i+1..$#sortedProt){
                my $prot2=$sortedProt[$j];
                next if $matchGroup{$prot2};    # already assigned to a match group
                my $matchList=1;
                ##>Comparing peptide contents of prot1 and prot2<##
                foreach my $pep (keys %{$matchList{$prot2}}){
                    if (! $matchList{$prot1}{$pep}){
                        delete $matchList{$prot1}{$pep}; # to be safe
                        $matchList=0;
                        last;
                    }
                }
                if ($matchList) {
                    $matchGroup{$prot2}=$matchGroup{$prot1};
                    $visibility{$prot2}=($protVisibility && defined($bestProtVis{$prot2}) && ($bestProtVis{$prot2}==2 || ($protVisibility==2 && $bestProtVis{$prot2})))? 1 : 0; # visible or hidden
                }
            }
        }

        ###>Inserting data into table ANALYSIS_PROTEIN
        $count=0;
        foreach my $protein (keys %{$protCoverage{$analysis}}){
            my $dbRank=$dbRankProt{$protein};
            my $coverage=$protCoverage{$analysis}{$protein};
            my $proteinID=$proteinID{$protein};
            my $protScore=($proteinScore{$protein})? $proteinScore{$protein} : '';
            my $numPep=scalar keys %{$assosProtPeptide{$analysis}{$protein}};
            my $numMatch=$numMatchProt{$analysis}{$protein};
            my $matchGroup=$matchGroup{$protein};
            my $pepSpecificity=$bestPepSpecificity{$protein};
            my $visibility=$visibility{$protein};
            $sthAnaProt->execute($analysisID,$proteinID,$dbRank,2,$protScore,$numPep,$numMatch,$coverage,$matchGroup,$pepSpecificity,$visibility);

        }
        $sthAnaProt->finish;
    }
    $dbh->commit;



    my %peptideID;
    $count=0;
    my $sthPeptide=$dbh->prepare("INSERT INTO PEPTIDE (ID_ANALYSIS,PEP_SEQ,PEP_LENGTH,QUERY_NUM,PEP_RANK,SEARCH_RANK,SCORE,MISS_CUT,MR_EXP,MR_CALC,MR_OBS,MR_DELTA,COMMENTS,SUBST,CHARGE,ELUTION_TIME,VALID_STATUS,DATA,SPEC_COUNT) VALUES (?,?,?,?,?,?,?,?,?,NULL,?,NULL,NULL,NULL,?,?,?,NULL,?)");
    my $sthModification=$dbh->prepare("INSERT INTO PEPTIDE_MODIFICATION (ID_PEPTIDE,ID_MODIFICATION,MODIF_TYPE,POS_STRING,REF_POS_STRING) VALUES (?,?,?,?,NULL)");
    my $sthAttrib=$dbh->prepare("INSERT INTO PEPTIDE_PROTEIN_ATTRIB (ID_PROTEIN,ID_PEPTIDE,ID_ANALYSIS,PEP_BEG,PEP_END,FLANKING_AA,IS_SPECIFIC) VALUES (?,?,?,?,?,?,?)");
    my $sthUpdatePeptide=$dbh->prepare("UPDATE PEPTIDE SET MR_CALC=?, MR_DELTA=? WHERE ID_PEPTIDE=?");
    open(FILESTAT,">>$fileStat");
    print FILESTAT "Storing peptides data into database.\n";
    close FILESTAT;
    foreach my $analysis (keys %analysisID){
        my $analysisID=$analysisID{$analysis};
        ##############################################################
        ###Inserting data into tables PEPTIDE,PEPTIDE_MODIFICATION ###
        ##############################################################
        foreach my $peptide (keys %{$peptideList{$analysis}}){
            my $ghostPep=0;
            foreach my $rt (sort {$peptideList{$analysis}{$peptide}{$a} <=> $peptideList{$analysis}{$peptide}{$b}} keys %{$peptideList{$analysis}{$peptide}}){
                my @pepSeq=split(/_/,$peptide);
                my $pepSequence=$pepSeq[0];
                if($format eq 'openswath'){$pepSequence=~s/\(UniMod:\d+\)//g;}
                else{$pepSequence=~s/\[\+?\w+\]//g;}
                my $pepLength=length $pepSequence;
                my $pepCharge=$pepSeq[1];
                my $mrObs=($peptideInfo{$peptide}[0]) ? $peptideInfo{$peptide}[0] : $peptideList{$analysis}{$peptide}{$rt}[0] ;
                my $mrExp=($mrObs-1.007825032)*$pepCharge;
                my $missCut=($peptideInfo{$peptide}[1]) ? $peptideInfo{$peptide}[1] : '' ;
                my $pepScore=($ghostPep) ? 'NULL' : ($peptideList{$analysis}{$peptide}{$rt}[1] && $peptideList{$analysis}{$peptide}{$rt}[1]=~/\d+/)? $peptideList{$analysis}{$peptide}{$rt}[1] : 'NULL';
                $rt=sprintf("%0.2f",$rt);
                my $queryNum=($ghostPep) ? '' : $peptideRTOrder{$analysis}{$peptide}[1];
                my $rank=($ghostPep)? '' : 1;
                my $validStatus=($ghostPep) ? 0 :1;
                my $elutionTime="sc".$pepScore.";et".$rt;
                ###>Insertion into PEPTIDE
                $sthPeptide->execute($analysisID,$pepSequence,$pepLength,$queryNum,$rank,$rank,$pepScore,$missCut,$mrExp,$mrObs,$pepCharge,$elutionTime,$validStatus,1);
                my $peptideID=$dbh->last_insert_id(undef,undef,'PEPTIDE','ID_PEPTIDE');
                push @{$peptideID{$analysis}{$peptide}},$peptideID;
                ###>Insertion into PEPTIDE_MODIFICATION
                my %varMods;
                foreach my $modificationID (keys %{$peptideModif{$analysis}{$peptide}}){
                    my $position="";
                    foreach my $pos (sort {$a cmp $b} keys %{$peptideModif{$analysis}{$peptide}{$modificationID}}){
                        $pos='=' if $pos==0;
                        $position=($position)? "$position.$pos" : $pos;
                        $varMods{$pos}=($pos eq '=') ? $peptideModif{$analysis}{$peptide}{$modificationID}{0} : $peptideModif{$analysis}{$peptide}{$modificationID}{$pos};
                    }
                    $sthModification->execute($peptideID,$modificationID,'V',$position);
                }
                my $mrCalc=&mrCalc($pepSequence,\%varMods);     ## $pepSequence = peptide sequence without modification ; %varMod => {modification position}=modification mass
                my $delta=$mrCalc-$mrExp;
                $sthUpdatePeptide->execute($mrCalc,$delta,$peptideID);
            }
        }
        $sthPeptide->finish;
        $sthUpdatePeptide->finish;
        $sthModification->finish;
        $dbh->commit;
        ########################################################
        ####Inserting data into tables PEPTIDE_PROTEIN_ATTRIB###
        ########################################################
        $count=0;
        foreach my $protein (keys $assosProtPeptide{$analysis}){
            my $proteinID=$proteinID{$protein};
            foreach my $peptide (keys $assosProtPeptide{$analysis}{$protein}){
                my $ghostPep=0;
                foreach my $peptideID (@{$peptideID{$analysis}{$peptide}}){
                    ###>Max peptide specificity
                    my $specificity;
                    $specificity=($peptideInfo{$peptide}[2] && $peptideInfo{$peptide}[2]==100) ? 1 : 0;
                    ###>Peptide position on protein
                    foreach my $begPepInfo (keys %{$peptidePosition{$analysis}{$protein}{$peptide}}){
                        my @begPeptideInfo=split(/_/,$begPepInfo);
                        my $pepBeg=($ghostPep)? "-$begPeptideInfo[0]" : $begPeptideInfo[0];
                        my $endPepInfo=$peptidePosition{$analysis}{$protein}{$peptide}{$begPepInfo};
                        my @endPeptideInfo=split(/_/,$peptidePosition{$analysis}{$protein}{$peptide}{$begPepInfo});
                        my $pepEnd=($ghostPep)? "-$endPeptideInfo[0]" : $endPeptideInfo[0];
                        my $flanking="$begPeptideInfo[1]$endPeptideInfo[1]";
                        $sthAttrib->execute($proteinID,$peptideID,$analysisID,$pepBeg,$pepEnd,$flanking,$specificity);
                    }
                    $ghostPep=1;
                }
            }
        }
        $sthAttrib->finish;
        $dbh->commit;
    }




    ################################################
    ####> Update VALID_STATUS in ANALYSIS table <###
    ################################################
    my $sthAnalysisUpdate=$dbh->prepare("UPDATE ANALYSIS SET VALID_STATUS=2 WHERE ID_ANALYSIS=?");
    foreach my $analysis (sort {$a cmp $b } keys %analysisID){
        my $analysisID=$analysisID{$analysis};
        $sthAnalysisUpdate->execute($analysisID);
    }
    $sthAnalysisUpdate->finish;
    $dbh->commit;


    #######################################################################
    ####Inserting data into tables QUANTIFICATION and ANA_QUANTIFICATION###
    #######################################################################
    my $code=($format eq 'prm')? 'TDA' : 'DIA';
    my $quantiName=($format eq 'pkv')? 'PeakView_Quantification' : ($format eq 'prm')? 'Skyline_Quantification' : ($format eq 'openswath')? 'OpenSwath_Quantification' : 'Spectronaut_Quantification';

    my $quantiMethod=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE=\"$code\"");
    my $sthQuanti=$dbh->prepare("INSERT INTO QUANTIFICATION (ID_MODIFICATION,ID_QUANTIFICATION_METHOD,ID_DESIGN,NAME,FOCUS,QUANTIF_ANNOT,STATUS,UPDATE_DATE,UPDATE_USER) VALUES (NULL,?,NULL,?,?,?,?,NOW(),?)");
    $sthQuanti->execute($quantiMethod,$quantiName,'peptide',$quantifAnnot,1,$userID) or die $dbh->errstr;
    $sthQuanti->finish;
    my $quantiID=$dbh->last_insert_id(undef, undef, 'QUANTIFICATION', 'ID_QUANTIFICATION');
    my $sthAnaQuanti=$dbh->prepare("INSERT INTO ANA_QUANTIFICATION (ID_QUANTIFICATION,ID_ANALYSIS,QUANTIF_FILE,IS_REFERENCE) VALUES (?,?,NULL,NULL)");
    ##Inserting data into ANA_QUANTIFICATION
    foreach my $analysis (keys %analysisID) {
        $sthAnaQuanti->execute($quantiID,$analysisID{$analysis});
    }
    $sthAnaQuanti->finish;
    $dbh->commit;

    #####################################
    ####<Launches identifier mapping>####
    #####################################
    my $anaStrg=join(",",values %analysisID);
    system "./mapProteinIdentifiers.pl $userID $anaStrg";
    
    ###>  recovering memory used to launch TRIC processus on the cluster -> temporary!!!  <###
    if (-s "$workDir/torqueID.txt"){
        my $torqueID=`grep 'torque6\.curie\.fr' $workDir/torqueID.txt`;
        $torqueID=~s/\.torque6\.curie\.fr//;
        $torqueID=~s/\s//g;
        my $serverInstance=&promsConfig::detectServer;
        my $server=($serverInstance=~/prod/) ? 'http' : 'w-myproms-ud';
        my $day=strftime("%Y%m%d",localtime);
        my $logDir="/data/tmp/torque6_logs/$server/$day.E";
        if(-s $logDir){
            my ($askedMem,$usedMem);
            open (DIR,"<$logDir");
            while(<DIR>){
                if ($_=~/$torqueID/){
                    my @infos=split(/\s/,$_);
                    foreach my $info(@infos){
                        if ($info=~/Resource_List.mem\=/){
                            ($askedMem=$info)=~s/Resource_List\.mem=//;
                        }
                        elsif ($info=~/resources_used.mem=/){
                            ($usedMem=$info)=~s/resources_used\.mem=//;
                        }
                    }
                }
            }
            close DIR;
            my $tmpDir="/data/tmp/myproms/cluster_info_OpenSwath";
            mkdir $tmpDir unless -e $tmpDir;
            open(OUT,">$tmpDir/memInfo_$time.txt");
            print OUT "asked mem : ",$askedMem,"\nused mem : ",$usedMem;
            close OUT;
        }
    }
    
    #########################
    ###Separting data file###
    #########################
    $count=0;
    my $analysisID;
    my $warning;
    open(FILESTAT,">>$fileStat");
    print FILESTAT "Split each data files by analysis.\n";
    close FILESTAT;
    foreach my $analysis (keys %analysisFileNames){
        $analysisID=$analysisID{$analysis};
        $anaReloadID=($anaReloadID)? $anaReloadID: $analysisID;
        my $anaDir="$promsPath{peptide}/proj_$projectID/ana_$analysisID";
        mkdir $promsPath{'peptide'} unless -e $promsPath{'peptide'};
        mkdir "$promsPath{peptide}/proj_$projectID" unless -e "$promsPath{peptide}/proj_$projectID";
        mkdir $anaDir unless -e $anaDir;

        ###> storing pyprophet results files to launch TRIC on all analysis in a same time
        if ($action eq 'quantification'){
            my $pyprophetFile="$analysis\_with_dscore$dsExtension.csv";
            if (`stat -c "%s" $workDir/$pyprophetFile`>10000000){          ### compress pyprophet file if its size > 1Go
                system "cd $workDir; tar -zcf $pyprophetFile.tar.gz $pyprophetFile";
                move ("$workDir/$pyprophetFile.tar.gz","$anaDir/$pyprophetFile.tar.gz");
            }
            else{
                move ("$workDir/$pyprophetFile","$anaDir/$pyprophetFile");
            }
        }

        my $quantiPath = $promsPath{"quantification"}."/project_$projectID/quanti_$quantiID";
        system "mkdir -p $quantiPath";
        my $swathFile = "$quantiPath/swath_ana_$analysisID.txt";
        open(OUTFILE, ">", $swathFile) or die ("open: $!");
        if($format eq 'pkv'){
            open(INFILE,"<","$sheetNames{'Area-ions'}");
            while (my $line=<INFILE>){
                my @partLine=split(/\t/,$line);
                my (@columnNames,@newLine);
                if ($partLine[0] eq "Protein") {
                    push @columnNames, "ID_PEPTIDE";
                    push @columnNames, $partLine[5];        #Fragment MZ
                    push @columnNames, $partLine[6];        #Fragment Charge
                    push @columnNames, $partLine[7];        #Ion Type
                    push @columnNames, $partLine[8];        #Residue
                    push @columnNames, "Area";
                    push @columnNames, $partLine[4];        #RT
                    my $colNames=join('!',@columnNames);
                    print OUTFILE $colNames;
                    print OUTFILE "\n";
                }
                else{
                    if (exists $peptideList{$analysis}{$partLine[1].'_'.$partLine[3]}) {
                        chomp $partLine[$analysis+2];
                        my $peptideID=$peptideID{$analysis}{$partLine[1].'_'.$partLine[3]}[0];
                        push @newLine,$peptideID;
                        push @newLine,$partLine[5];     #Fragment MZ
                        push @newLine,$partLine[6];     #Fragment Charge
                        push @newLine,$partLine[7];     #Ion Type
                        push @newLine,$partLine[8];     #Residue
                        push @newLine,$partLine[$analysis+2]; #Area
                        #my $rt=($partLine[4]=~/None/) ? $partLine[4] : sprintf("%0.2f",$partLine[4]);
                        #push @newLine,$rt;        #RT
                        my $newLine=join('!',@newLine);
                        print OUTFILE $newLine;
                        print OUTFILE "\n";
                    }
                }
            }
            close INFILE;
            close OUTFILE;


            ###>Generating score files
            my $analysisName=fileparse($analysisFileNames{$analysis}, qr/\.[^.]*/);
            my $targetFile="$anaDir/$analysisName.target";
            my $decoyFile="$anaDir/$analysisName.decoy";
            open (TARGET,">$targetFile");
            open (DECOY,">$decoyFile");
            foreach my $type (keys %{$scoreList{$analysis}}){
                if ($type eq 'TARGET') {
                    foreach my $score (@{$scoreList{$analysis}{$type}}){
                        print TARGET "$score\n";
                    }
                }
                elsif($type eq 'DECOY'){
                    foreach my $score (@{$scoreList{$analysis}{$type}}){
                        print DECOY "$score\n";
                    }
                }
            }
            close TARGET;
            close DECOY;


            ###>Drawing density plot with R
            my $numCycles=0;
            my $numTarget=$targetNb{$analysis};
            my $numDecoy=$decoyNb{$analysis};
            my $Rscript="$anaDir/density.R";
            open(R,">$Rscript");
            print R qq
|target <- t(read.table("$targetFile"))
decoy <- t(read.table("$decoyFile"))
densT <- density(target)
densD <- density(decoy)

minX <- min(c(min(densT\$x,densD\$x)))
maxX <- max(c(max(densT\$x,densD\$x)))
maxY <- max(c(max(densT\$y*$numTarget,densD\$y*$numDecoy)))

png(filename="$anaDir/scores.png", width=500, height=500, units = "px")

plot(densT\$x,densT\$y*$numTarget, xlim=c(minX,maxX), ylim=c(0,maxY*1.1), yaxs='i', type='l', col='green', xlab="Scores", ylab="Abundance", main="$analysisName ($fdr% FDR precomputed)")
lines(densD\$x,densD\$y*$numDecoy, col='red')
legend("topright", c("Peptides","Decoy peptides"), cex=0.8, col=c('green','red'), lwd=2, bty="n")

dev.off()
|;
            close R;
            system "cd $anaDir; $promsPath{R}/R CMD BATCH --no-save --no-restore $Rscript";
            sleep 1;
            unlink $Rscript;
            unlink $Rscript.'out' if -e $Rscript.'out';

        }
        elsif($format eq 'openswath'){
            print OUTFILE "ID_PEPTIDE!Fragment MZ!Fragment Charge!Ion Type!Residue!Area!RT\n";
            foreach my $peptide (keys %{$peptideList{$analysis}}){
                my $ghostPep=0;
                foreach my $RT (sort {$peptideList{$analysis}{$peptide}{$a} <=> $peptideList{$analysis}{$peptide}{$b}} keys %{$peptideList{$analysis}{$peptide}}){
                    my ($pepMZ,$score,$fragmentsList,$areasList)=@{$peptideList{$analysis}{$peptide}{$RT}};
                    my $peptideID=($ghostPep) ? $peptideID{$analysis}{$peptide}[1] : $peptideID{$analysis}{$peptide}[0];
                    my @areas=split(/;/,$areasList);
                    my @fragments=split(/;/,$fragmentsList);
                    foreach my $i (0 .. $#fragments){
                        my ($fragmentID,$fragmentType,$fragmentCharge,@infos)=split(/_/,$fragments[$i]);
                        my $fragmentMZ=$fragMass{$peptide}{$fragmentType.'_'.$fragmentCharge};
                        unless ($fragmentMZ || $warning) {
                            $warning="Few fragments have no M/Z information. <BR>It might be due to a wrong selection of library or non exclusion of N-term modifications during export.";
                        }

                        $fragmentMZ=($fragmentMZ)?$fragmentMZ : '';
                        my ($ionType,$residue);
                        $fragmentType=~s/\s//;
                        if ($fragmentType=~/([a-z]{1})(.*)/){
                            $ionType=$1;
                            $residue=$2;
                        }
                        my $area=$areas[$i];
                        $area=~s/,/\./ if $area=~/,/;
                        print OUTFILE "$peptideID!$fragmentMZ!$fragmentCharge!$ionType!$residue!$area!\n";
                    }
                    $ghostPep=1;
                }
            }
            close OUTFILE;
        }
        else{
            print OUTFILE "ID_PEPTIDE!Fragment MZ!Fragment Charge!Ion Type!Residue!Area!RT\n";
            foreach my $peptide (keys %{$fragmentsInfos{$analysis}}){
                my $ghostpep=0;
                foreach my $rt (sort {$fragmentsInfos{$analysis}{$peptide}{$a} <=> $fragmentsInfos{$analysis}{$peptide}{$b}} keys %{$fragmentsInfos{$analysis}{$peptide}}){
                    foreach my $fragment (keys %{$fragmentsInfos{$analysis}{$peptide}{$rt}}){
                        my ($fragmentType,$fragmentCharge)=split(/_/,$fragment);
                        my ($ionType,$residue);
                        $fragmentType=~s/\s//;
                        if ($fragmentType=~/([a-z]{1})(.*)/){
                            $ionType=$1;
                            $residue=$2;
                        }
                        my $fragmentRT=($fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[0]) ? $fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[0] : '';
                        my $fragmentMZ=($fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[2]) ? $fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[2] : ($fragMass{$peptide}{$fragment})? $fragMass{$peptide}{$fragment} : '';
                        my $fragmentNeutralMass=$fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[3];
                        my $area=$fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[1];
                        my $peptideID=($ghostpep) ? $peptideID{$analysis}{$peptide}[1] :  $peptideID{$analysis}{$peptide}[0] ;
                        print OUTFILE "$peptideID!$fragmentMZ!$fragmentCharge!$ionType!$residue!$area!$fragmentRT\n";
                    }
                    $ghostpep=1;
                }
            }
            close OUTFILE;
        }
        symlink($swathFile,"$anaDir/$file");
    }
    open(FILESTAT,">>$fileStat");
    print FILESTAT "Done.\n";
    print FILESTAT "WARNING : ",$warning,"\n" if $warning;
    close(FILESTAT);

    system "rm $fileStat" if (-e $fileStat && !$warning) ;
#    } # END FORK

    print qq |
<SCRIPT LANGUAGE="JavaScript">
var monitorWindow=window.open("$promsPath{cgi}/monitorDIAProcess.cgi?ACT=import",'Monitoring import of DIA data','width=1000,height=500,scrollbars=yes,resizable=yes');
monitorWindow.focus();
top.promsFrame.selectedAction='summary';
parent.navFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=experiment:$experimentID&ACT=nav";
</SCRIPT>
	|;
    print "</CENTER></BODY></HTML>";
}



##########################
####Check delta length#### Compares delta lengthes between 2 proteins with Peptide density of the smaller one
########################## (delta is considered not significant if smaller than pep density)
sub deltaLength {
	my ($l_a,$l_b,$d_a,$d_b)=@_;
	my $pepDensVal=($l_b > $l_a)? $d_a : $d_b;
	if (abs($l_b-$l_a) > $pepDensVal) {return $l_a<=>$l_b} else {return 0}
}

sub mrCalc {
    my ($pepSeq,$refVarMods)=@_;    ## pepSeq => peptide sequence without modifications ; %varMods => {modification position}=modification mass
    my $mrCalc=0;
    my %massValueAA=&promsConfig::getMassAAmono;
    my @aaArray=split(//,$pepSeq);
    push @aaArray, 'C_term';
    unshift @aaArray, 'N_term';
    for (my $i=1; $i<$#aaArray; $i++){
        if($i==1){      ##N-term
            $mrCalc+=$massValueAA{$aaArray[0]};
            $mrCalc+=$refVarMods->{'='} if defined $refVarMods->{'='}; #N-term modification
        }
        if($i==$#aaArray-1){        ##C_term
            $mrCalc+=$massValueAA{$aaArray[$i+1]};
            $mrCalc+=$refVarMods->{'*'} if defined $refVarMods->{'*'}; # C-term modification
        }
        $mrCalc+=$massValueAA{$aaArray[$i]} if ($massValueAA{$aaArray[$i]});
        print "$aaArray[$i],$pepSeq,@aaArray" unless $massValueAA{$aaArray[$i]};
        exit unless $massValueAA{$aaArray[$i]};
        $mrCalc+=$refVarMods->{$i} if defined $refVarMods->{$i};
    }
    return $mrCalc;
}

sub excludeLibraryModifications {
    my $libID=param('libID');
    my $sthModification=$dbh->prepare("SELECT SLM.ID_MODIFICATION, PSI_MS_NAME FROM MODIFICATION M, SWATH_LIB_MODIFICATION SLM WHERE ID_SWATH_LIB=$libID AND SLM.ID_MODIFICATION=M.ID_MODIFICATION");
    $sthModification->execute;
    if ($sthModification->rows != 0){
        print "<BR><B>&nbsp;- Peptide modifications to be excluded :</B><BR>";
        while(my ($modifID,$modifName)=$sthModification->fetchrow_array){
            print "<INPUT type=\"checkbox\" name=\"pepMod\" value=\"$modifID\">$modifName&nbsp;";
        }
    }
	$dbh->disconnect;
}

sub selectSamples {
    my $expID=param('expID');
    my $sthListSamples=$dbh->prepare("SELECT ID_SAMPLE,NAME FROM SAMPLE WHERE ID_EXPERIMENT=$expID");
    $sthListSamples->execute;
    if ($sthListSamples->rows != 0){
        my $i=1;
        while (my ($sampleID,$name)=$sthListSamples->fetchrow_array){
            print "<INPUT type=\"checkbox\" name=\"samples\" id=\"samples$expID$i\" value=\"$sampleID\" checked>$name&nbsp;";
            $i++;
        }
    }
}


####>Revision history<####
# 1.3.8i Move to cluster end of TRIC process which is reading and storing feature_alignment_fileout.tsv into database (GA 23/11/18)
# 1.3.8h Minor modification on the exportLibrarySCRIPT call (VS 22/11/18)
# 1.3.8g Minor modification on the getProtInfo call (VS 16/11/2018)
# 1.3.8f Changed licence to CeCILL (VS 09/11/18)
# 1.3.8e Update swath files path (VS 08/11/18)
# 1.3.8d Test of TRIC version with filtered files (GA 20/09/18)
# 1.3.8c Set memory request to 200Go for TRIC (PP 01/09/18)
# 1.3.8b Modification in feature_alignment.py option equivalent to 1.4.1 dev script version (GA 01/08/18)
# 1.3.8 Minor modif to handle library export errors (MLP 17/04/18)
# 1.3.7 Recovering used memory for jobs launch on the cluster (MLP 26/03/18)
# 1.3.6 Update for auto-increment in table QUANTIFICATION (MLP 20/03/18) 
# 1.3.5 Min 40Gb for TRIC (PP 08/03/18)
# 1.3.4 minor modif (MLP 27/02/18)
# 1.3.3 Run import process in background -> monitorDIAProcess.cgi (MLP 08/02/18)
# 1.3.2 Minor modif on peakview file import (MLP 18/01/18)
# 1.3.1 Modification on samples selection for OpenSwath quantification (MLP 17/01/18)
# 1.3.0 Add spectronaut data import (MLP 12/01/18)
# 1.2.10 Minor modif (MLP 11/01/18)
# 1.2.9 Minor modif (MLP 11/01/18)
# 1.2.8 Added sleep in job-wait loops to prevent disk access overload & other minor changes (PP 28/12/17)
# 1.2.7 Minoir modif (MLP 20/12/17)
# 1.2.6 Modif to overcome server error during OpenSwath workflow process (MLP 18/12/17)
# 1.2.5 Modifications on fragment MZ recovering (MLP 15/12/17)
# 1.2.4 Minoir modif (MLP 14/12/17)
# 1.2.3 Minoir modif (MLP 14/12/17)
# 1.2.2 Launch multiple openswath workflow in the same time (MLP 14/12/17)
# 1.2.1 Add option : shared directory for files import (MLP 13/12/2017)
# 1.2.0 Add openswath quantification workflow (MLP 06/12/17)
# 1.1.4 Modif on the selection of modification for TDA data (MLP 08/11/17)
# 1.1.3 Minor Modif (24/10/17)
# 1.1.2 Minor Modif (03/10/17)
# 1.1.1 Minor modif (27/09/17)
# 1.1.0 Adapt script for OpenSwath (MLP 27/09/17)
# 1.0.9 Adapt script for TDA (MLP 27/07/17)
# 1.0.8 Minor modification on protein alias (MLP 22/03/17)
# 1.0.7 Modification for peptide's inversion (Leucine and Isoleucine) (MLP 03/02/17)
# 1.0.6 Minor beug correction (MLP 31/01/17)
# 1.0.5 Minor modif (MLP 30/01/17)
# 1.0.4 Minor update to store statistics library into QUANTIF_ANNOT (MLP 26/01/17)
# 1.0.3 Add modifications for multiple databank library (MLP 01/12/16)
# 1.0.2 Add &promsMod::cleanParameters verification and display IDENTIFIER_TYPE in SWATH library selection (MLP 20/10/16)
# 1.0.1 Minor beug correction (MLP 12/08/16)
# 1.0.0 Created (MLP 04/07/16)