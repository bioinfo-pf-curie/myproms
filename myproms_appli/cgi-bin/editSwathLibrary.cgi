#!/usr/local/bin/perl -w

################################################################################
# editSwathLibrary.cgi         1.8.9	                                       #
# Authors: M. Le Picard, V. Sabatet (Institut Curie)                           #
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
use File::Copy;
use File::Basename;
use XML::Simple;
use POSIX qw(strftime); # to get the time
use XML::SAX::ParserFactory;
use List::Util "first";
use File::Spec::Functions qw(splitpath); # Core module

#######################
####>Starting HTML<####
#######################
print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);

my %promsPath=&promsConfig::getServerInfo;
my %clusterInfo=&promsConfig::getClusterInfo;#('debian'); # default is 'centos'
my $tppPath=($clusterInfo{'on'}) ? $clusterInfo{'path'}{'tpp'} : $promsPath{'tpp'};
my $pythonPath=($clusterInfo{'on'})? $clusterInfo{'path'}{'python'} : $promsPath{'python'};
my $MAX_NB_THREAD = 4; # Number of threads to use when parallelization is available
		
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my $userID=$ENV{'REMOTE_USER'};

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

#############################
####>Fetching parameters<####
#############################
my $action=(param('ACT'))? param('ACT') : "" ;
my $libraryID=(param('ID'))? param('ID') : 0; # 0 if ACT=add
if ($action eq 'delete') {&deleteSwathLib;exit;}
if ($action eq 'archive') {&archiveSwathLib;exit;}
my $projectID=(param('projectID'))? param('projectID') : 'NULL';
my $experimentID=(param('experimentID'))? param('experimentID') : 0;
if ($action eq 'ajaxSelectExperiment') {&ajaxSelectExperiment; exit;}
if ($action eq 'ajaxSelectSample') {&ajaxSelectSample; exit;}
my $option=(param('option'))? param('option') : 0;
if ($action eq 'selectLibraryName') {&selectLibraryName; exit;}
if ($action eq 'selectSecondMergeLib') {&ajaxSelectSecondMergeLib; exit;}
my $submit=(param('submit'))? param('submit') : "";
if ($action eq 'rtdata') {&ajaxRTData; exit;}
if ($action eq 'selectDBMergeLib') {&selectDBMergeLib; exit;}
if ($action eq 'ajaxSelectSpecieDB') {&ajaxSelectSpecieDB; exit;}
if ($action eq 'restore') {&restorePreviousVersion;}
	

mkdir $promsPath{'data'} unless -e $promsPath{'data'};
mkdir "$promsPath{data}/tmp" unless -e "$promsPath{data}/tmp";
mkdir "$promsPath{data}/tmp/Swath" unless -e "$promsPath{data}/tmp/Swath";

#############################
####<Deleting a library>####
#############################
sub deleteSwathLib {
	my $sthSwathAna=$dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS_SWATH_LIB WHERE ID_SWATH_LIB=?");		#delete library if not used
	$sthSwathAna->execute($libraryID);
	if ($sthSwathAna->rows == 0) {
        $dbh->do("DELETE FROM DATABANK_SWATHLIB WHERE ID_SWATH_LIB=$libraryID")|| die $dbh->errstr;
		$dbh->do("DELETE FROM SWATH_LIB_MODIFICATION WHERE ID_SWATH_LIB=$libraryID")|| die $dbh->errstr;
		my $sthSwathParent1=$dbh->prepare("SELECT ID_PARENT_SWATH_LIB FROM PARENT_SWATH_LIB WHERE ID_SWATH_LIB=?");
		my $sthSwathParent2=$dbh->prepare("SELECT ID_SWATH_LIB FROM PARENT_SWATH_LIB WHERE ID_PARENT_SWATH_LIB=?");
		my $sthSwathNameParent=$dbh->prepare("SELECT NAME FROM SWATH_LIB WHERE ID_SWATH_LIB=?");
		$sthSwathParent1->execute($libraryID);
		$sthSwathParent2->execute($libraryID);
		if ($sthSwathParent1->rows != 0) {			#check if this library was created with an other library (merge, update)
			$dbh->do("DELETE FROM PARENT_SWATH_LIB WHERE ID_SWATH_LIB=$libraryID")|| die $dbh->errstr;
		}
		if ($sthSwathParent2->rows != 0) {
            while (my $swathParent=$sthSwathParent2->fetchrow_array){
                $sthSwathNameParent->execute($swathParent);
				if ($sthSwathNameParent->rows == 0) {
					$dbh->do("DELETE FROM PARENT_SWATH_LIB WHERE ID_SWATH_LIB=$libraryID")|| die $dbh->errstr;
				}
				else{
					my $libraryName=$dbh->selectrow_array("SELECT NAME FROM SWATH_LIB WHERE ID_SWATH_LIB=$libraryID");
					my $alert="You have to delete : ";
					while (my $name=$sthSwathNameParent->fetchrow_array){
                        $alert.="\"".$name."\" ";
                    }
                    $alert.="before \"$libraryName\".";
					print "<HTML><BODY><SCRIPT LANGUAGE=\"JavaScript\">alert('$alert'); window.location=\"". $promsPath{cgi} . "/listSwathLibraries.cgi\";</SCRIPT></BODY></HTML>";
					exit;
				}
            }
        }
		$sthSwathNameParent->finish;
		$sthSwathParent2->finish;
		$sthSwathParent1->finish;

		$dbh->do("DELETE FROM SWATH_LIB WHERE ID_SWATH_LIB=$libraryID")|| die $dbh->errstr;
    }
    else{
		$dbh->do("UPDATE SWATH_LIB SET USE_STATUS='no' WHERE ID_SWATH_LIB=$libraryID") || die $dbh->errstr;			# use_status='no' if library is used
	}

	
	$sthSwathAna->finish;
	$dbh->commit;
	$dbh->disconnect;
	
	
	####<Deleting all corresponding files>####
	system "rm -r $promsPath{swath_lib}/SwLib_$libraryID" if -e "$promsPath{swath_lib}/SwLib_$libraryID";
	print "<HTML><BODY><SCRIPT LANGUAGE=\"JavaScript\">window.location=\"$promsPath{cgi}/listSwathLibraries.cgi\";</SCRIPT></BODY></HTML>";
	exit;
}



###########################
####<Archive a library>####
###########################
sub archiveSwathLib {
	my $workDir="$promsPath{swath_lib}/SwLib_$libraryID";
	
	###> deleting archives for each library version (keep only the files of the consensus library)
	opendir (DIR, $workDir);
	while(my $file=readdir(DIR)){
		if ($file=~/\.tar\.gz/){
			system "rm $workDir/$file";
		}
	}
	close DIR;
	
	$dbh->do("UPDATE SWATH_LIB SET USE_STATUS='arc' WHERE ID_SWATH_LIB=$libraryID");
	$dbh->commit;
	$dbh->disconnect;
	
	print "<HTML><BODY><SCRIPT LANGUAGE=\"JavaScript\">window.location=\"$promsPath{cgi}/listSwathLibraries.cgi\";</SCRIPT></BODY></HTML>";
	
	exit;
}

#####################################################
###<Return to the previous version of the library>###
#####################################################
sub restorePreviousVersion{
	my ($des,$date,$user,$db,$param,$numPep,$numProt,$numProtSpe,$firstLine);
	
	my $libPath="$promsPath{swath_lib}/SwLib_$libraryID";
	my $versionPath="$promsPath{swath_lib}/SwLib_$libraryID/oldVersion";
	my ($libraryName,$versionLib)=$dbh->selectrow_array("SELECT NAME,VERSION_NAME FROM SWATH_LIB WHERE ID_SWATH_LIB=$libraryID");
	my ($dbInfoFile,$fileListFile);
	(my $prevVersion=$versionLib)=~s/v//;
	$prevVersion-=1;
	if (-e $versionPath){
		$dbInfoFile="dbinfo.txt";
		$fileListFile="filelist.txt";
	}
	else{
		
		print qq |
<HTML>
<HEAD>
<TITLE>Library</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR><BR><BR><BR><IMG src='$promsPath{images}/engrenage.gif'><BR>
<BR><BR><FONT class="title1">Restoring <FONT color="red"> $libraryName </FONT> previous version</FONT></BR></BR>
|;
		$dbInfoFile='dbinfo_Lib'.$libraryID.'_v'.$prevVersion.'.txt';
		$fileListFile='filelist_Lib'.$libraryID.'_v'.$prevVersion.'.txt';
		system "cd $libPath; rm Lib$libraryID\_$versionLib.tar.gz";
		system "cd $libPath; gunzip Lib$libraryID\_v$prevVersion.tar.gz ; tar -xf Lib$libraryID\_v$prevVersion.tar $dbInfoFile; tar -xf Lib$libraryID\_v$prevVersion.tar $fileListFile; gzip Lib$libraryID\_v$prevVersion.tar";
	}
	
	my $version;
	open(DBINFO,"<","$libPath/$dbInfoFile") or die ("open: $!");
	while (<DBINFO>) {
		if ($_=~/des=/) {
			$firstLine=$_;
			my @tabInfo=split(/\t/,$_);
			foreach my $info (@tabInfo){
				if ($info=~/des=/) {
					($des=$info)=~s/des=//;
				}
				elsif($info=~/version=/){
					($version=$info)=~s/version=//;
				}
				elsif($info=~/date=/){
					($date=$info)=~s/date=//;
				}
				elsif($info=~/user=/){
					($user=$info)=~s/user=//;
				}
				elsif($info=~/pep=/){
					($numPep=$info)=~s/pep=//;
				}
				elsif($info=~/prot=/){
					($numProt=$info)=~s/prot=//;
				}
				elsif($info=~/protspe=/){
					($numProtSpe=$info)=~s/protspe=//;
				}
				elsif($info=~/databank=/){
					($db=$info)=~s/databank=//;
				}
			}	
		}
		else{
			$param.=$_;
		}
		
	}
	close DBINFO;
	
	
	my $fragmentation;
	if ($param){
		($fragmentation)=($param=~/<option>-cI<\/option>\s*<opvalue>(\w+)<\/opvalue>/);
	}
	$fragmentation=($fragmentation)? $fragmentation : 'HCD';
	
	my $entry->{'NUM_PEP'}=$numPep;
	$entry->{'NUM_PROT'}=$numProt;
	$entry->{'NUM_PROT_SPECIFIC'}=$numProtSpe;
	my $xmlParser = XML::Simple->new( NoAttr=>1, RootName=>'NUM_ENTRY');
	my $stat=$xmlParser->XMLout($entry);
	
	
	my %fileList;
	
	if(-e $versionPath){		### for old libraries
		system "mv $versionPath/* $libPath";
		system "rm -r $versionPath";
		
		if (-e "$libPath/$fileListFile") {
			##### supprimer noms des fichiers actuels !!!
			open(FILELIST, "<","$libPath/$fileListFile") or die $!;
			while (my $line=<FILELIST>) {
				chomp $line;
				my @file=split(/&/,$line);
				foreach my $file (@file){
					my ($fileName,$versionFile)=split(/,/,$file);
					$versionFile=~s/v//;
					if ($versionFile<$prevVersion) {
						push @{$fileList{$versionFile}},$fileName;
					}
				}
			}
			close FILELIST;
			my %paramList;
			my $xmlParser = XML::Simple->new( NoAttr=>1, RootName=>'SWATH_LIB_DATA');
			my $xmlParams=$xmlParser->XMLin($param);
			foreach my $paramFiles (@{$xmlParams->{'FILES'}->{'FILE'}}){
				(my $versionFile=$paramFiles->{'version'})=~s/v//;
				if ($versionFile==$prevVersion) {
					push @{$paramList{$paramFiles->{'version'}}},$paramFiles->{'FileName'};
				}
			}
			@{$xmlParams->{'FILES'}->{'FILE'}}=();
			foreach my $versionFile (keys %paramList){
				foreach my $file (@{$paramList{$versionFile}}){
					push @{$xmlParams-> {'FILES'} -> {'FILE'}},{'FileName'=>"$file",'version'=>"$versionFile"};
				}
			}
			$param=$xmlParser->XMLout($xmlParams);
		}
		else{ 	### for older libraries !! (before the PARAM_STRG modification)
			my %paramList;
			my $xmlParser = XML::Simple->new( NoAttr=>1, RootName=>'SWATH_LIB_DATA');
			my $xmlParams=$xmlParser->XMLin($param);
			foreach my $paramFiles (@{$xmlParams->{'FILES'}->{'FILE'}}){
				(my $versionFile=$paramFiles->{'version'})=~s/v//;
				if ($versionFile<$prevVersion) {
					push @{$fileList{$versionFile}},$paramFiles->{'FileName'};
				}
				else{
					push @{$paramList{$paramFiles->{'version'}}},$paramFiles->{'FileName'};
				}
			}
			@{$xmlParams->{'FILES'}->{'FILE'}}=();
			foreach my $versionFile (keys %paramList){
				foreach my $file (@{$paramList{$versionFile}}){
					push @{$xmlParams-> {'FILES'} -> {'FILE'}},{'FileName'=>"$file",'version'=>"$versionFile"};
				}
			}
			$param=$xmlParser->XMLout($xmlParams);
			
			open(DBINFO,">","$libPath/$dbInfoFile") or die ("open: $!");
			print DBINFO "$firstLine\n";
			print DBINFO $param;
			close DBINFO;
		}
		if (%fileList) {
			open(FILELIST, ">","$libPath/$fileListFile") or die $!;
			foreach my $version (%fileList){
				foreach my $fileName (@{$fileList{$version}}){
					print FILELIST "$fileName,v$version&";
				}
			}
			close FILELIST;
		}
		else{
			system "rm $libPath/$fileListFile";			##first version of the library
		}
		
		system "rm $libPath/$dbInfoFile" if -e "$libPath/$dbInfoFile";
	}
	else{		
		system "rm -f $libPath/$libraryName.sptxt $libPath/$libraryName.spidx $libPath/$libraryName.splib $libPath/$libraryName.pepidx $libPath/$libraryName.mrm $libPath/sortie.txt $libPath/script.sh"; 
		print "<BR><BR> Fetching associated files ...";
		opendir (DIR,$libPath);
		my $count=0;
		while ( my $file=readdir(DIR) ){
			if ($file =~/\.tar\.gz/){
				system "cd $libPath; tar -zxf $file";
				(my $folderName=$file)=~s/\.tar\.gz//;
				system "rm $libPath/dbinfo_$folderName.txt;" if -e "$libPath/dbinfo_$folderName.txt";
				system "rm $libPath/filelist_$folderName.txt;" if -e "$libPath/filelist_$folderName.txt";
				
				if ($folderName eq "Lib$libraryID\_v$prevVersion"){
					move ("$libPath/script_$folderName.sh","$libPath/script.sh");
					move ("$libPath/sortie_$folderName.txt","$libPath/sortie.txt");
					system "cd $libPath; gunzip $file; tar --delete --file=$folderName.tar script_$folderName.sh; tar --delete --file=$folderName.tar sortie_$folderName.txt; tar --delete --file=$folderName.tar dbinfo_$folderName.txt; tar --delete --file=$folderName.tar filelist_$folderName.txt; gzip $folderName.tar";
				}
				else{
					system "rm $libPath/sortie_$folderName.txt;" if -e "$libPath/sortie_$folderName.txt";
					system "rm $libPath/script_$folderName.sh;" if -e "$libPath/script_$folderName.sh";
				}
			}
			print ".";
			$count++;
		}
		close DIR;
		
		
		my $outputFile=$libPath."/restaureOut.txt";
		print "<BR><BR> Processing concensus library generation ...";
		system "cd $libPath; $promsPath{'tpp'}/spectrast -cN$libPath/SpecLib_cons -cI$fragmentation -cAC $libPath/SpecLib*.splib >>$outputFile 2>&1 &";
		my ($wait,$error)=(0,0);
		$count=0;
		while (!$wait && !$error){
			if (-s "$libPath/SpecLib_cons.sptxt" && -s "$libPath/SpecLib_cons.splib" && -s "$libPath/SpecLib_cons.pepidx" && -s "$libPath/SpecLib_cons.spidx"){
				$wait=1;
			}
			if ($count==2){
				print ".";
				$count=0;
			}
			if (-s $outputFile){
				if (`tail -1 $outputFile`=~/SpectraST finished at/){
					$wait=1;
				}
				elsif (`tail -1 $outputFile`=~/error/){
					$error=1;
				}
			}
			$count++;
			sleep 30;
		}
		move ("$libPath/SpecLib_cons.sptxt","$libPath/$libraryName.sptxt");
		move ("$libPath/SpecLib_cons.splib","$libPath/$libraryName.splib");
		move ("$libPath/SpecLib_cons.pepidx","$libPath/$libraryName.pepidx");
		move ("$libPath/SpecLib_cons.spidx","$libPath/$libraryName.spidx");
		move ("$libPath/SpecLib_cons.mrm","$libPath/$libraryName.mrm") if(-s "$libPath/SpecLib_cons.spidx");
		system "cd $libPath;rm SpecLib_*; rm $outputFile; rm spectrast.log;";
		print "</CENTER></BODY></HTML>";
		
		print "<BR><BR> Storing data ...";
	}
	
	
	if ($date) {
		my $sthLibRestore=$dbh->prepare("UPDATE SWATH_LIB SET DES=?, PARAM_STRG=?, VERSION_NAME=?,UPDATE_DATE=?,UPDATE_USER=?,STATISTICS=?  WHERE ID_SWATH_LIB=?") or die "Couldn't prepare statement: " . $dbh->errstr;
		$sthLibRestore->execute($des,$param,"v$prevVersion",$date,$user,$stat,$libraryID);
		$sthLibRestore->finish;
	}
	else{
		my $sthLibRestore=$dbh->prepare("UPDATE SWATH_LIB SET DES=?, PARAM_STRG=?, VERSION_NAME=?,UPDATE_DATE=NULL ,UPDATE_USER=NULL, STATISTICS=?  WHERE ID_SWATH_LIB=?") or die "Couldn't prepare statement: " . $dbh->errstr;
		$sthLibRestore->execute($des,$param,"v$prevVersion",$stat,$libraryID);
		$sthLibRestore->finish;
	}
	
	my $sthLibDBRestore=$dbh->prepare("UPDATE DATABANK_SWATHLIB SET ID_DATABANK=? WHERE ID_SWATH_LIB=?") or die "Couldn't prepare statement: " . $dbh->errstr;
	$sthLibDBRestore->execute($db,$libraryID);
	$sthLibDBRestore->finish;
	$dbh->commit;
	$dbh->disconnect;
	
	
	print "<HTML><BODY><SCRIPT LANGUAGE=\"JavaScript\">window.location=\"$promsPath{cgi}/listSwathLibraries.cgi\";</SCRIPT></BODY></HTML>";
	exit;
}


#################################
####<Add or update a library>####
#################################

###> Form to select files
print qq
|<HTML>
<HEAD>
<TITLE>Library</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.row_0{background-color:$darkColor;}
.row_1{background-color:$lightColor;}
.highlight{width:250px;}
.checklist {
    border: 1px solid #ccc;
    list-style:none;
    height:120px;
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
<SCRIPT LANGUAGE="JavaScript">
    |;
    &promsMod::popupInfo();
	&promsMod::browseDirectory_JavaScript if $promsPath{'shared'};
    print qq |
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


    function ajaxSelectExperiment(projID) {
        document.getElementById('experimentSPAN').innerHTML="";
        document.getElementById('sampleSPAN').innerHTML="";
        if (projID != "") {
            //Creation of the XMLHTTPRequest object
            var XHR = getXMLHTTP();
            XHR.onreadystatechange=function() {
                if (XHR.readyState==4 && XHR.responseText) {
                   document.getElementById('experimentSPAN').innerHTML=XHR.responseText;
                }
            }
            XHR.open("GET","$promsPath{cgi}/editSwathLibrary.cgi?ACT=ajaxSelectExperiment&projectID="+projID,true);
            XHR.send(null);
        }

    }

    function ajaxSelectSample(expID) {
        document.getElementById('sampleSPAN').innerHTML="";
        if(expID != ""){
            //Creation of the XMLHTTPRequest object
            var XHR = getXMLHTTP();
            XHR.onreadystatechange=function() {
                if (XHR.readyState==4 && XHR.responseText) {
                   document.getElementById('sampleSPAN').innerHTML=XHR.responseText;
                }
            }
            XHR.open("GET","$promsPath{cgi}/editSwathLibrary.cgi?ACT=ajaxSelectSample&experimentID="+expID,true);
            XHR.send(null);
        }

    }

    function ajaxSelectLibrary(option){
        document.getElementById('libNameDIV').style.display='none';
        document.getElementById('libraryName').innerHTML = "";
        document.getElementById('chooseLibrary').innerHTML="";
        if (option == 'mergelib'){
            //Creation of the XMLHTTPRequest object
            var XHR = getXMLHTTP();
            XHR.onreadystatechange=function() {
                if (XHR.readyState==4 && XHR.responseText) {
                   document.getElementById('chooseLibrary').innerHTML=XHR.responseText;
                }
            }
            XHR.open("GET","$promsPath{cgi}/editSwathLibrary.cgi?ACT=selectLibraryName&option="+option,true);
            XHR.send(null);
        }
    }

    function selectSecondMergeLib(libID1){
        //Creation of the XMLHTTPRequest object
        var XHR = getXMLHTTP();
        XHR.onreadystatechange=function() {
            if (XHR.readyState==4 && XHR.responseText) {
               document.getElementById('secondSELECT').innerHTML=XHR.responseText;
            }
        }
        XHR.open("GET","$promsPath{cgi}/editSwathLibrary.cgi?ACT=selectSecondMergeLib&libID1="+libID1,true);
        XHR.send(null);
    }

	function ajaxSelectRT(rtID){
		document.getElementById('rtdata').innerHTML="";
		var XHR = getXMLHTTP();
		if(rtID!=0){
			XHR.onreadystatechange=function() {
				if (XHR.readyState==4 && XHR.responseText) {
				   document.getElementById('rtdata').innerHTML=XHR.responseText;
				}
			}
			XHR.open("GET","$promsPath{cgi}/editSwathLibrary.cgi?ACT=rtdata&rtID="+rtID,true);
			XHR.send(null);
		}
	}

	function ajaxSelectSpecie(species){
		if(species){
			document.getElementById('db').style.display='';
			document.getElementById('db').innerHTML='';
			var XHR = getXMLHTTP();
			XHR.onreadystatechange=function() {
				if (XHR.readyState==4 && XHR.responseText) {
				   document.getElementById('db').innerHTML=XHR.responseText;
				}
			}
			XHR.open("GET","$promsPath{cgi}/editSwathLibrary.cgi?ACT=ajaxSelectSpecieDB&species="+species,true);
			XHR.send(null);
		}
	}


    function selectLibraryName(option,split){
        document.getElementById('libNameDIV').style.display='';
        if ((option == 'new') \|\| (option == 'mergelib')){
            var libname="&nbsp;<INPUT size=45% type='text' name='libname' id='libname' pattern='[0-9a-zA-Z_-]+' required>(Accepted 0-9,a-z,A-Z,_,-)";
			document.getElementById('db').style.display='';
        }
		document.getElementById('libraryName').innerHTML=libname;
		if(option == 'mergelib'){
			if (document.getElementById('libfile').selectedOptions[0].text != "-= Select Library =-"){
				var libInfo=document.getElementById('libfile').selectedOptions[0].value;
				var splits=libInfo.split('&');
				var split=splits[1]
				var libID=splits[0];
				if (split == 1){
					document.getElementById("unsplit").checked = false;
					document.getElementById("split").checked = true;
				}
				else if (split == 0){
					document.getElementById("unsplit").checked = true;
					document.getElementById("split").checked = false;
				}
				document.getElementById('db').style.display='';
				var XHR = getXMLHTTP();
				XHR.onreadystatechange=function() {
					if (XHR.readyState==4 && XHR.responseText) {
					   document.getElementById('db').innerHTML=XHR.responseText;
					}
				}
				XHR.open("GET","$promsPath{cgi}/editSwathLibrary.cgi?ACT=selectDBMergeLib&ID="+libID,true);
				XHR.send(null);
			}
			else
                document.getElementById('libNameDIV').style.display='none';
		}
    }
	// <--- AJAX
	
	
	function checkFileForm (form){
		
		if (form.serverDir.value){
			return true;
		}
		else if (form.inputfiles1.files.length==0 && form.inputfiles2.files.length==0 && form.inputfiles3.files.length==0 && form.archfiles.files.length==0 && !form.sharedDirFiles){
			alert('ERROR: Select data files.');
			return false;
		}
		var nameMzXML=[],
				nameDat=[],
				nameTandem=[],
				nameSequest=[];
		if(form.inputfiles1.files.length \|\| form.inputfiles2.files.length \|\| form.inputfiles3.files.length){
			
			var nbMzXML=0,
				nbDat=0,
				nbTandem=0,
				nbSequest=0;
			for (var j=1; j<4 ; j++){
				var files = document.getElementById('inputfiles'+ j).files;
				for (var i=0; i<files.length; i++){
					var fileInfo=files[i].name.split(".");
					var fileType=fileInfo[fileInfo.length-1];
					if (fileType == 'mzXML'){
						nameMzXML.push(fileInfo[0]);
						nbMzXML++;
					}
					else if (fileType == 'dat'){
						nameDat.push(fileInfo[0]);
						nbDat++;
					}
					else if ((fileType == 'xml' && fileInfo[fileInfo.length-3] == 'tandem') \|\| fileType == 'tandem' \|\| (fileType == 'xml' && fileInfo[fileInfo.length-2] != 'pep')){			//tandem.pep.xml ; .tandem ; .xml -> xtandem
						nameTandem.push(fileInfo[0]);
						nbTandem++;
					}
					else if (fileType == 'xml' && fileInfo[fileInfo.length-2] == 'pep'){		//pep.xml -> sequest
						nameSequest.push(fileInfo[0]);
						nbSequest++;
					}
				}
			}
		}
		else if (form.sharedDirFiles.length){		// multiple files found
			for(let i=0 ; i<form.sharedDirFiles.length ; i++){
				var file=form.sharedDirFiles[i].value;
				if (form.sharedDirFiles[i].checked==true) {
					var fileInfo=file.split(".");
					var fileType=fileInfo[fileInfo.length-1];
					if(fileType == 'dat'){
						nameDat.push(fileInfo[0]);
					}
					else if ((fileType == 'xml' && fileInfo[fileInfo.length-3] == 'tandem') \|\| fileType == 'tandem' \|\| (fileType == 'xml' && fileInfo[fileInfo.length-2] != 'pep')){
						nameTandem.push(fileInfo[0]);
					}
					else if (fileType == 'xml' && fileInfo[fileInfo.length-2] == 'pep'){
						nameSequest.push(fileInfo[0]);
					}
					else if (fileType == 'mzXML'){
						nameMzXML.push(fileInfo[0]);
					}
				}
			}
			if(nameMzXML.length==0){
				alert("Wrong files number");
				return false;
			}
		}
		else if (form.sharedDirFiles.checked \|\| !form.sharedDirFiles.checked){	// single file found
			alert("Wrong files number");
			return false;
		}
		if ((nameDat.length && nameMzXML.length != nameDat.length)  \|\| (nameTandem.length && nameMzXML.length != nameTandem.length) \|\| (nameSequest.length && nameMzXML.length != nameSequest.length)){
				alert("Wrong files number");
				return false;
		}
		nameMzXML.sort();
		nameTandem.sort();
		nameDat.sort();
		nameSequest.sort();
		if ((nameDat.length && nameDat.join() !== nameMzXML.join()) \|\| (nameTandem.length && nameTandem.join() !== nameMzXML.join()) \|\| (nameSequest.length && nameSequest.join() !== nameMzXML.join())){
			alert("Wrong files names");
			return false;
		}
	}
	
	function selectSource (source){
		var myForm=document.parameters;
		
		myForm.projectFiles.disabled=true;
		myForm.archfiles.disabled=true;
		myForm.serverDir.disabled=true;
		for(var i=1; i<=3; i++){
			document.getElementById('inputfiles'+i).disabled=true;
		}
		if (document.getElementById('sharedDirDIV')) {document.getElementById('sharedDirDIV').style.display='none';}
		
		if(source=='UseLocalDirectory'){
			for(var i=1; i<=3; i++){
				document.getElementById('inputfiles'+i).disabled=false;
			}
		}
		else if(source=='UseProjectDirectory'){
			myForm.projectFiles.disabled=false;
		}
		else if(source=='UseArchiveFiles'){
			myForm.archfiles.disabled=false;
		}
		else if(source=='UseServerDirectory'){
			myForm.serverDir.disabled=false;
		}
		else if(source=='UseSharedDirectory'){
			document.getElementById('sharedDirDIV').style.display='';
		}
		
	}

</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
    <DIV style="float:top">
        <BR>
        <TABLE><TR><TH bgcolor="$darkColor">
            <FONT class="title2">&nbsp;Go to:</FONT>
            <SELECT class="title action" onchange="if (this.value) {window.location='$promsPath{cgi}/'+this.value;}">
                <OPTION value="">-= Select =-</OPTION>
                <OPTION value="promsMain.cgi">Main Window</OPTION>
                <OPTION value="selectProject.cgi">Project Selection</OPTION>
            </SELECT>
        </TH></TR></TABLE>
    </DIV>
    <CENTER>|;
if ($submit eq "") {
    if ($action eq "add" || $action eq "update" || $action eq "addDB2") {
		my ($splitLibUpdate,$selectLibraryName,$fdrType);
        if ($action eq 'update'){
            ($selectLibraryName,$splitLibUpdate,my $paramStrg)=$dbh->selectrow_array("SELECT NAME,SPLIT,PARAM_STRG FROM SWATH_LIB WHERE ID_SWATH_LIB='$libraryID' ") or die "Couldn't prepare statement: " . $dbh->errstr;
            print qq
            |<FONT class="title1">Update <FONT class="title1" color="#DD0000">$selectLibraryName</FONT></FONT><BR>
            <BR><BR><BR>
            <FORM method="POST" action="./editSwathLibrary.cgi" name="parameters" enctype="multipart/form-data" onsubmit="return(checkFileForm(this));">
                <TABLE bgcolor="$darkColor">

            |;
			
			my $xml=XML::Simple->new(KeepRoot=>1);
			my $param=$xml->XMLin($paramStrg);
			foreach my $command (@{$param->{'SWATH_LIB_DATA'}->{'COMMANDS_LINES'}->{'COMMAND'}}){
				if($command->{'CommandName'} eq 'Mayu.pl'){
					foreach my $params (@{$command->{'PARAMS'}}){
						if($params->{'option'} eq '-fdr'){
							$fdrType=$params->{'opvalue'};
						}
					}
				}
			}
			$fdrType=($fdrType)? $fdrType : '';
        }
        else{
            print qq
            |<FONT class="title1">Adding new spectral library to Server</FONT><BR><BR><BR><BR>
            <FORM method="POST" action="./editSwathLibrary.cgi" name="parameters" enctype="multipart/form-data" onsubmit="return(checkFileForm(this));">
            <TABLE bgcolor="$darkColor">
            <TR><TH align="right" valign="top">Task : </TH>
            <TD bgcolor="$lightColor">
            <INPUT type="radio" name="libcons" id="libcons" value="mergelib" onClick="ajaxSelectLibrary(this.value)" >Merge with an other library<BR>
            <INPUT type="radio" name="libcons" id="libcons" value="new" onClick="ajaxSelectLibrary(this.value)" onChange="selectLibraryName(this.value)">Create new library<BR>
            <SPAN id="chooseLibrary"></SPAN></TD></TR>
            <TR id="libNameDIV" style="display:none"><TH align=right>Library name : </TH><TD bgcolor='$lightColor'><SPAN id="libraryName"></SPAN></TD></TR>
			|;
			my %speciesHash;
            my $sthSpeciesList=$dbh->prepare("SELECT ID_SPECIES,COMMON_NAME,SCIENTIFIC_NAME FROM SPECIES WHERE IS_REFERENCE=1") or die "Couldn't prepare statement: " . $dbh->errstr;
            $sthSpeciesList->execute;
            if ($sthSpeciesList->rows==0) {print "No species.";}
            else{
                while (my($speciesID,$commonSpeciesName,$scientSpeciesName)=$sthSpeciesList->fetchrow_array){
                    @{$speciesHash{$speciesID}}=($commonSpeciesName,$scientSpeciesName);
                }
            }
            $sthSpeciesList->finish;
			print qq
			|<TR><TH align="right" valign="top">Species :</TH>
			<TD bgcolor='$lightColor'>
			&nbsp;<SELECT name="species" onChange="ajaxSelectSpecie(this.value)">
                <option value="">-= Select species =-</option>
            |;
            foreach my $speciesID (sort {$speciesHash{$a}[1] cmp $speciesHash{$b}[1]} keys %speciesHash){
				my $speciesName=($speciesHash{$speciesID}[1])? $speciesHash{$speciesID}[1]: $speciesHash{$speciesID}[0];
                print "<option value=\"$speciesHash{$speciesID}[0]_$speciesHash{$speciesID}[1]\">$speciesName</option>";
            }
            print qq
            |</SELECT><BR>&nbsp;<SMALL>(Update list of reference species if your species is not listed)</SMALL>
			</TD>
			</TR>
            |;
        }
		print "<INPUT type=\"hidden\" name=\"libfilename\" value=\"$selectLibraryName\">" if $action eq 'update';
        print qq
        |<INPUT type="hidden" name="ACT" id="ACT" value="$action"><INPUT type="hidden" name="libfile" value="$libraryID">
		 <TR><TH align="right" valign="top">&nbspConsensus library options : </TH>
           <TD bgcolor="$lightColor"><INPUT onmouseover="popup('This option additionnaly considers retention times when merging spectra.')" onmouseout="popout()" type="radio" name="liboption" id="split" value="split" |; if ($action eq 'update' && $splitLibUpdate==1) {print "checked";} else {print "required";} print qq |>Split
              <INPUT onmouseover="popup('This option assumes that all fragment ion spectra are correctly assigned.')" onmouseout="popout()" type="radio" name="liboption" value="unsplit" id="unsplit" |; if ($action eq 'update' && $splitLibUpdate==0) {print "checked";} print qq |>Unsplit<BR>
           </TD>
        </TR>
        <TR><TH align="right" valign="top">Files : </TH>
		<TD bgcolor="$lightColor"><B><INPUT type="radio" name="selSource" value="UseLocalDirectory" onclick="selectSource(this.value);"> Import files : </B><BR>
            &nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="file" accept=".dat,.mzXML,.xml" name="inputfiles1" id="inputfiles1" multiple="multiple" disabled ><BR>
            &nbsp;&nbsp;&nbsp;&nbsp;<INPUT accept=".xml,.dat,.mzXML" type="file" name="inputfiles2" id="inputfiles2" multiple="multiple" disabled ><BR>
            &nbsp;&nbsp;&nbsp;&nbsp;<INPUT accept=".xml,.dat,.mzXML" type="file" name="inputfiles3" id="inputfiles3" multiple="multiple" disabled><BR>
            &nbsp;&nbsp;&nbsp;&nbsp;(.dat, .mzXML and .tandem.pep.xml)<BR><BR>
            <B><INPUT type="radio" name="selSource" value="UseProjectDirectory" onclick="selectSource(this.value);"> Import from project :</B>
            |;
            my %projectHash;
            my $sthProjectList=$dbh->prepare("SELECT NAME, ID_PROJECT FROM PROJECT") or die "Couldn't prepare statement: " . $dbh->errstr;
            $sthProjectList->execute;
            if ($sthProjectList->rows==0) {print "No project. Create a new project.";}
            else{
                while (my($projectName,$projID)=$sthProjectList->fetchrow_array){
                    $projectHash{$projectName}=$projID;
                }
            }
            $sthProjectList->finish;

            print qq
            |<TABLE>
            <TR><TD valign="top">&nbsp;&nbsp;&nbsp;&nbsp;<SELECT name="projectFiles" onChange="ajaxSelectExperiment(this.value)" disabled>
                <option value="">-= Select Project =-</option>
            |;
            foreach my $projectName (sort {lc $a cmp lc $b} keys %projectHash){
                print "<option value=\"$projectHash{$projectName}\">$projectName</option>";
            }

            print qq
            |</SELECT>
			</TD>
            <TD valign="top"><SPAN id="experimentSPAN"></SPAN></TD>
            <TD valign="top"><SPAN id="sampleSPAN" ></SPAN></TD>
            </TR></TABLE>
			
			<BR><B><INPUT type="radio" name="selSource" value="UseArchiveFiles" onclick="selectSource(this.value);"> Import archive file :</B><BR>
			&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="file" accept=".zip,.gz" name="archfiles" id="archfiles" disabled><BR><SMALL>&nbsp;&nbsp;&nbsp;&nbsp;(zip or gz archive)</SMALL><BR>
			
			<BR><B><INPUT type="radio" name="selSource" value="UseServerDirectory" onclick="selectSource(this.value);"> Select directory from server :</B><BR>
			&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="text" name="serverDir" id="serverDir" disabled><BR><SMALL>&nbsp;&nbsp;&nbsp;&nbsp;(only for bioinformatician)</SMALL><BR>
			
			<BR><B><INPUT type="radio" name="selSource" value="UseSharedDirectory" onclick="selectSource(this.value);"> Shared data directory </B><BR>
			&nbsp;&nbsp;&nbsp;&nbsp;<DIV id="sharedDirDIV" style="width:600px;max-height:300px;overflow:auto;display:none">
			|;
			&promsMod::browseDirectory_getFiles($promsPath{'shared'},'sharedDirFiles',{fileMatch=>qr/\.(dat|xml|mzXML)$/i});
			print qq |
			</DIV><BR>
			</TD>
		</TR>	
        <TR><TH align="right">PTMProphet : </TH>
			<TD bgcolor="$lightColor">
				<label> <input type="checkbox" onclick="document.getElementById('PTMProphetOptions').style.display = (this.checked) ? '' : 'none';" name="PTMProphet" value="1" /> Use PTMs localization improvement (PTMProphet)</label><br/>
				<FIELDSET id="PTMProphetOptions" style='display:none; margin-top: 10px;'>
					<LEGEND><B>Options:</B></LEGEND>
					<label> Model to apply :&nbsp;
						<SELECT name ="PTMProphetEMModel" required>
							<option value="">-= Select fragmentation type =-</option>
							<option value="0">No EM Model</option>
							<option value="1">Intensity EM Model</option>
							<option value="2" selected>Intensity and matched peaks EM Model</option>
						</SELECT>
					</label><br/>
					<label> Min probability to evaluate peptide : <input type='text' name='PTMProphetMinProb' value='0.9' size='3' /></label></br>
					<label> MS1 PPM tolerance : <input stype='text' name='PTMProphetPPMTol' value='1' size='3' /></label><br/>
					<label> Fragments PPM tolerance : <input type='text' name='PTMProphetFragPPMTol' value='10' size='3' /></label><br/>
					Modifications to process :
							<label><input type='checkbox' name='PTMProphetModifs' value='STY:79.9663' />Phospho (S,T,Y)&nbsp;&nbsp;</label>
							<label><input type='checkbox' name='PTMProphetModifs' value='M:15.9949' />Oxydation (M)&nbsp;&nbsp;</label>
							<label><input type='checkbox' name='PTMProphetModifs' value='C:57.0214' />Carbamidomethyl (C)&nbsp;&nbsp;</label>
				</FIELDSET>
            </TD>
		</TR>
        <TR><TH align="right">Fragmentation type : </TH>
            <TD bgcolor="$lightColor">
            &nbsp;<SELECT name ="fragmentation" required>
			<option value="">-= Select fragmentation type =-</option>
				<option value="CID-QTOF">CID-QTOF</option>
				<option value="HCD">HCD</option>
				<option value="ETD">ETD</option>
				<option value="CID">CID</option>
            </SELECT>
            </TD>
		</TR>
		<TR><TH align="right">Instrument : </TH>
            <TD bgcolor="$lightColor">
            |;
            my $sthInstrumList=$dbh->prepare("SELECT NAME FROM INSTRUMENT") or die "Couldn't prepare statement: " . $dbh->errstr;
            $sthInstrumList->execute;
            print qq
            |&nbsp;<SELECT name ="instrument" >
                <option value="">-= Select Instrument =-</option>
            |;
            if ($sthInstrumList->rows==0) {print "No experiment in that project. Choose an other project.";}
            else{
                while (my ($instrumName)=$sthInstrumList->fetchrow_array){
                    print "<option value=\"$instrumName\"";
                    if ($instrumName eq 'ESI-TRAP' ){
                        print " selected";
                    }
                    print ">$instrumName</option>";
                }
            }
            $sthInstrumList->finish;
            print qq
            |</SELECT>
            </TD>
        </TR>|;
		if ($action eq 'update'){
			print "<TR><TH align=\"right\">Databank : </TH><TD bgcolor=\"$lightColor\">";
			my ($nameDBLib,$organismLib,$DBTypeLibUpdate)=$dbh->selectrow_array("SELECT D.NAME, ORGANISM, DT.ID_DBTYPE FROM DATABANK D, DATABANK_SWATHLIB DS, DATABANK_TYPE DT WHERE ID_SWATH_LIB=$libraryID AND DS.ID_DATABANK=D.ID_DATABANK AND D.ID_DBTYPE=DT.ID_DBTYPE");
			my $sthDBUpdateList=$dbh->prepare("SELECT D.NAME, D.ID_DATABANK,D.FASTA_FILE,D.DECOY_TAG,DT.NAME FROM DATABANK D,DATABANK_TYPE DT WHERE D.DECOY_TAG='yes' or D.DECOY_TAG LIKE '%rev%' AND D.USE_STATUS='yes' AND D.ID_DBTYPE=DT.ID_DBTYPE AND ORGANISM=\"$organismLib\" AND D.ID_DBTYPE=$DBTypeLibUpdate") or die "Couldn't prepare statement: " . $dbh->errstr;
			$sthDBUpdateList->execute;
			if ($sthDBUpdateList->rows==0) {
                print "No databank available for SWATH.";
            }
            else{
                my %dbHash;
                while (my($dbName,$dbID,$dbFileName,$decoyTag,$dbTypeName)=$sthDBUpdateList->fetchrow_array){
                    next if ($dbFileName=~m/:/) ;
                    next if ($decoyTag eq "No" || $decoyTag eq "");
                    $dbHash{"$dbName&nbsp;\[$dbTypeName\]"}=$dbID;
                }
                
                print qq
                |&nbsp;<SELECT name ="dbfile" required>
                    <option value="">-= Select Databank =-</option>
                |;
                foreach my $dbName (sort {lc $a cmp lc $b} keys %dbHash){
                    print "<option value=\"$dbHash{$dbName}\">$dbName</option>";
                }
				print "</SELECT></TD></TR>";
            }
			$sthDBUpdateList->finish;
		}
		else{
        print qq |<TR>

			<TH align="right">Databank : </TH>

            <TD bgcolor="$lightColor"><SPAN id="db">
            |;
            my $sthDBList=$dbh->prepare("SELECT D.NAME, D.ID_DATABANK,D.FASTA_FILE,D.DECOY_TAG,DT.NAME FROM DATABANK D,DATABANK_TYPE DT WHERE D.DECOY_TAG='yes' or D.DECOY_TAG LIKE '%rev%' AND D.USE_STATUS='yes' AND D.ID_DBTYPE=DT.ID_DBTYPE ") or die "Couldn't prepare statement: " . $dbh->errstr;
            $sthDBList->execute;

            if ($sthDBList->rows==0) {
                print "No databank available for SWATH.";
            }
            else{
                my %dbHash;
                while (my($dbName,$dbID,$dbFileName,$decoyTag,$dbTypeName)=$sthDBList->fetchrow_array){
                    next if ($dbFileName=~m/:/) ;
                    next if ($decoyTag eq "No" || $decoyTag eq "");
                    $dbHash{"$dbName&nbsp;\[$dbTypeName\]"}=$dbID;
                }
                
                print qq
                |&nbsp;<SELECT name ="dbfile" required>
                    <option value="">-= Select Databank =-</option>
                |;
                foreach my $dbName (sort {lc $a cmp lc $b} keys %dbHash){
                    print "<option value=\"$dbHash{$dbName}\">$dbName</option>";
                }
				print "</SELECT>";
            }
			$sthDBList->finish;
            print "</SPAN></TD></TR>";
		}

        print qq |<TR><TH align="right" valign="top">Mayu options: </TH>
            <TD bgcolor="$lightColor">&nbsp;FDR estimation with Mayu software. <BR>
                &nbsp;Missed cleavage :  <INPUT type="text" value="2" name="missedcleav" size="3" required><BR>
                &nbsp;FDR :  <INPUT type="text" value="0.01" name="fdr" size="4" required>
		|;
		if($action eq 'update' && $fdrType){
			print qq |
				Type : <SELECT name="fdrtype" required>
						<option value="$fdrType" selected>$fdrType</option>";
					   </SELECT>
			|;
		}else{
			print qq |
				Type : <SELECT name="fdrtype" required>
					<option value="mFDR">mFDR</option>
					<option value="pepFDR">pepFDR</option>
					<option value="protFDR" selected>protFDR</option>
					</SELECT>
			|;
		}
		print qq |		
            </TD>
        </TR>

        <TR><TH align="right" valign="top">RT file : </TH>
            <TD bgcolor="$lightColor">
            |;
            my $sthRTList=$dbh->prepare("SELECT NAME,ID_REFERENCE_RT FROM REFERENCE_RT") or die "Couldn't prepare statement: " . $dbh->errstr;
            $sthRTList->execute;
			my %rtHash;
            if ($sthRTList->rows==0) {print "No experiment in that project. Choose an other project.";}
            else{
                while (my($rtName,$rtID)=$sthRTList->fetchrow_array){
                    $rtHash{$rtName}=$rtID;
                }
                print qq
                |&nbsp;<SELECT name="refrtfile"  onChange="ajaxSelectRT(this.value)" required>
                    <option value="">-= Select RT File =-</option>
                |;
                foreach my $rtName (sort {lc $a cmp lc $b} keys %rtHash){
                    print "<option value=\"$rtHash{$rtName}\" >$rtName</option>";
				}
            }
			$sthRTList->finish;
            print "</SELECT>";
			print "<DIV id=\"rtdata\"></DIV>";
			print qq |
            </TD>
        </TR>


        <!--<TR><TH align="right">Convertion to Windows file : </TH>
           <TD bgcolor="$lightColor"><INPUT type="radio" name="convdos" value="oui"  >Oui<INPUT type="radio" name="convdos" value="non">Non</TD>
        </TR>--!>


        <TR><TH align="right" valign="top">Description : </TH>
           <TD bgcolor="$lightColor">&nbsp;<TEXTAREA rows="4" cols="43" name="descript" id="descript"></TEXTAREA></TD>
        </TR>
        <TR ><TH colspan=2><BR><input type="submit" name="submit" value="Submit" > 
        <!-- CLEAR button -->
        &nbsp &nbsp &nbsp<INPUT type="reset" value="Clear" />
        <!-- CANCEL button -->
        &nbsp &nbsp &nbsp<INPUT type="button" value="Cancel" onclick="window.location='./listSwathLibraries.cgi'"></TH></TR>
        </TABLE>
        </FORM>
        <DIV id="divDescription" class="clDescriptionCont">
        <!--Empty div-->
        </DIV>
        <SCRIPT type="text/javascript">setPopup()</SCRIPT>
        |;
    }

    elsif ($action eq "merge"){
        print qq |
        <FONT class="title1">Merge two SWATH libraries</FONT><BR><BR><BR><BR>
        <FORM method="POST" action="./editSwathLibrary.cgi" name="parameters" enctype="multipart/form-data">
        <TABLE bgcolor="$darkColor">
        <TR><TH align="right" valign="top">Select libraries : </TH>
        <TD bgcolor="$lightColor">
        |;
        my $sthLibList=$dbh->prepare("SELECT SL.NAME,SL.ID_SWATH_LIB,DT.NAME FROM SWATH_LIB SL, DATABANK_TYPE DT,DATABANK D, DATABANK_SWATHLIB DS WHERE SL.USE_STATUS='yes' AND DS.ID_SWATH_LIB=SL.ID_SWATH_LIB AND DS.ID_DATABANK=D.ID_DATABANK AND DS.DB_RANK=1 AND D.ID_DBTYPE=DT.ID_DBTYPE") or die "Couldn't prepare statement: " . $dbh->errstr;
        $sthLibList->execute;
        if ($sthLibList->rows==0) {print "No library available.";}
        else{
            my %libHash;
            while (my($libName,$libID,$dbType)=$sthLibList->fetchrow_array){
                $libHash{$libName}="$libID&$dbType";
            }
            
            print qq
            |<SELECT name="libfile1" id="libfile1" onChange="selectSecondMergeLib(this.value)" required>
                <option value="">-= Select Library =-</option>
            |;
            foreach my $libName (sort {lc $a cmp lc $b} keys %libHash){
				my ($libId,$dbType)=split(/&/,$libHash{$libName});
                print "<option value=\"$libId\">$libName [$dbType]</option>";
            }
            print "</SELECT><SPAN id=\"secondSELECT\"></SPAN></TD></TR>";
        }
		$sthLibList->finish;
		
        print qq |
        <INPUT type="hidden" name="ACT" id="ACT" value="$action">
        <TR><TH align="right" valign="top">Library name : </TH><TD bgcolor=\"$lightColor\">
        <INPUT onmouseover="popup('Accepted 0-9,a-z,A-Z,_,-')" onmouseout="popout()" size=45% type="text" name="libname" id="libname" pattern="[0-9a-zA-Z_-]+" required></TD></TR>
        <TR><TH align="right" valign="top">Description : </TH>
            <TD bgcolor="$lightColor"><TEXTAREA rows="4" cols="43" name="descript" id="descript"></TEXTAREA></TD>
        </TR>
        <TR ><TH colspan=2><BR><input type="submit" name="submit" value="Submit">
        <!-- CLEAR button -->
        &nbsp &nbsp &nbsp<INPUT type="reset" value="Clear" />
        <!-- CANCEL button -->
        &nbsp &nbsp &nbsp<INPUT type="button" value="Cancel" onclick="window.location='./listSwathLibraries.cgi'"></TH></TR>
        </TABLE></FORM>
        |;
    }
    elsif($action eq 'edit'){
        my ($selectLibraryName)=$dbh->selectrow_array("SELECT NAME FROM SWATH_LIB WHERE ID_SWATH_LIB='$libraryID' ") or die "Couldn't prepare statement: " . $dbh->errstr;
        print qq |
        <FONT class="title1">Edit  <FONT class="title1" color="#DD0000">$selectLibraryName</FONT></FONT><BR><BR><BR><BR>
        <FORM method="POST" action="./editSwathLibrary.cgi" name="parameters" enctype="multipart/form-data">
        <TABLE bgcolor="$darkColor">
        <INPUT type="hidden" name="ACT" id="ACT" value="$action">
        <INPUT type="hidden" name="ID" id="ID" value="$libraryID">
        <TR><TH align="right" valign="top">New name : </TH>
        <TD bgcolor="$lightColor"><INPUT onmouseover="popup('Accepted 0-9,a-z,A-Z,_,-')" onmouseout="popout()" size=45% type='text' name='libname' id='libname' pattern='[0-9a-zA-Z_-]+'></TD></TR>
        <TR><TH align="right" valign="top">Description : </TH>
            <TD bgcolor="$lightColor"><TEXTAREA rows="4" cols="43" name="descript" id="descript"></TEXTAREA></TD>
        </TR>
        <TR ><TH colspan=2><BR><input type="submit" name="submit" value="Submit">
        <!-- CLEAR button -->
        &nbsp &nbsp &nbsp<INPUT type="reset" value="Clear" />
        <!-- CANCEL button -->
        &nbsp &nbsp &nbsp<INPUT type="button" value="Cancel" onclick="window.location='./listSwathLibraries.cgi'"></TH></TR>
        </TABLE></FORM>
        |;
    }
}
else{
	if ($action eq 'edit') {
		print "<BR><BR><IMG src='$promsPath{images}/engrenage.gif'><BR><BR>";
		my ($libNewName,$des)=&promsMod::cleanParameters(param('libname'),param('descript'));
		
		my ($libName,$split)=$dbh->selectrow_array("SELECT NAME,SPLIT FROM SWATH_LIB WHERE ID_SWATH_LIB=$libraryID");
		my $sthLibUpdate=$dbh->prepare ("UPDATE SWATH_LIB SET NAME=? ,DES=? WHERE ID_SWATH_LIB=$libraryID");
		$sthLibUpdate->execute($libNewName,$des);
		$sthLibUpdate->finish;
		if ($split==0) {
			system "mv $promsPath{swath_lib}/SwLib_$libraryID/$libName.sptxt $promsPath{swath_lib}/SwLib_$libraryID/$libNewName.sptxt";
			system "mv $promsPath{swath_lib}/SwLib_$libraryID/$libName.splib $promsPath{swath_lib}/SwLib_$libraryID/$libNewName.splib";
			system "mv $promsPath{swath_lib}/SwLib_$libraryID/$libName.spidx $promsPath{swath_lib}/SwLib_$libraryID/$libNewName.spidx";
			system "mv $promsPath{swath_lib}/SwLib_$libraryID/$libName.pepidx $promsPath{swath_lib}/SwLib_$libraryID/$libNewName.pepidx";
			system "mv $promsPath{swath_lib}/SwLib_$libraryID/$libName.mrm $promsPath{swath_lib}/SwLib_$libraryID/$libNewName.mrm" if(-s "$promsPath{swath_lib}/SwLib_$libraryID/$libName.mrm");
		}
		elsif($split==1){
			system "mv $promsPath{swath_lib}/SwLib_$libraryID/$libName.sptxt $promsPath{swath_lib}/SwLib_$libraryID/$libNewName.sptxt";
		}
		print "<B>Done.</B>";
		$dbh->commit;
		$dbh->disconnect;
		print "<HTML><BODY><SCRIPT LANGUAGE=\"JavaScript\">window.location=\"$promsPath{cgi}/listSwathLibraries.cgi\";</SCRIPT></BODY></HTML>";
		exit;
	}
	
	my ($dbType,@dbID);
	my $startTime=strftime("%s",localtime);
	print "<BR><BR><IMG src='$promsPath{images}/engrenage.gif'><BR><BR>";
	print "<FONT color=\"red\"><DIV id=\"waitDIV\"></DIV></FONT><BR>";		
	## wait time
	my $divID="document.getElementById('waitDIV')";
	my $now=strftime("%s",localtime); # in sec
	my $waitTime=strftime("%Hh %Mm %Ss",localtime($now-$startTime-3600));
	my $status="Run time : $waitTime ";
	
		
	my (@inputFiles,$oldVersion,$instrument,$fdrType,$des,$workDir,$decoyTag,$missedCleavage,$identifierDB,$organismDB,$fastaDBFile,$fdr,$libName,$libPath,$libIDFile,$libFileName,$rtFileID,$split,$libName1,$versionLib1,$libName2,$versionLib2,$libID1,$libID2);	
	my $libOption=(param('liboption'))? param('liboption') :'' ;
	my $libCons=(param('libcons'))? param('libcons') :'' ;
	my $fragmentation=(param('fragmentation'))? param('fragmentation') : 'HCD';
	my ($time,$outputFile,$rtFileName,$dbFile,$dbDir,@inputFilesDat,@inputFilesMZXML,@inputFilesTandem,@inputFilesXML,@inputFilesSequest,@inputFilesPEPXML,@oldFiles,$maxMem,@taille,@tailleSort,$nbFile);
	if ($action eq "add" || $action eq "update" || $action eq "addDB2") {
		
		###> Create new directory
		$time=strftime("%Y_%m_%d_%H_%M_%S",localtime);
		
		##>delete older tmp folders
		foreach my $folder (glob "$promsPath{data}/tmp/Swath/*") {		
			next unless $folder;
			if (-M $folder && -M $folder>6) {
				system "rm -r $folder";
			}
			print "<B> .</B>";
		}
	
		if ($clusterInfo{'on'}) {
			print "<SCRIPT LANGUAGE=\"JavaScript\">$divID.innerHTML=\"\";$divID.innerHTML='$status';</SCRIPT>";
		}
		$workDir=(param('workDir'))? param('workDir') : "$promsPath{data}/tmp/Swath/$time";
		$libPath="$promsPath{swath_lib}";
		mkdir $libPath unless -e $libPath;
		mkdir $workDir unless -e $workDir;
		$outputFile=$workDir.'/sortie.txt';
		
		$libOption=param('liboption');
		$libCons=(param('libcons'))? param('libcons') :'' ;
		my $convDos=param('convdos');
		
		push @dbID,param('dbfile');
		
		$libFileName=param('libfilename') if $action eq 'update';
		if ($libCons eq 'mergelib') {
			($libIDFile,my $splitOption,$versionLib1,my $dbMergeID,$libFileName)=split(/&/,param('libfile'));
			push @dbID,$dbMergeID if "@dbID"!~/$dbMergeID/;
		}
		else{$libIDFile=param('libfile');}
		
		
		
		$rtFileID=param('refrtfile');
		$instrument=param('instrument');
		
		
		my $serverDir=(param('serverDir')) ? &promsMod::cleanParameters(param('serverDir')) : '';
		($des,$missedCleavage,$fdr,$libName)=&promsMod::cleanParameters(param('descript'),param('missedcleav'),param('fdr'),param('libname'));
		$fdrType=param('fdrtype');
		
		print "<BR><BR><BR><BR><B><FONT size=\"4pt\">Updating <FONT color=\"red\">$libFileName</FONT></FONT></B><BR><BR>" if $action eq "update";
		print "<BR><BR><BR><BR><B><FONT size=\"4pt\">Creation of <FONT color=\"red\">$libName</FONT></FONT></B><BR><BR>" unless $action eq "update";
	
		
		###################################
		####<Processing submitted form>####
		###################################
		my $script;
		$script="editSwathLibrary.cgi?ACT=add" if $action eq 'add';
		$script="editSwathLibrary.cgi?ACT=update&ID=$libIDFile" if $action eq 'update';
		
		
		##########################
		####>Recovering files<####
		##########################
		###> Recovering fasta file's name chosen by the user
		($fastaDBFile,$identifierDB,$organismDB,$decoyTag)=$dbh->selectrow_array("SELECT D.FASTA_FILE,DT.DEF_IDENT_TYPE,D.ORGANISM,D.DECOY_TAG FROM DATABANK D, DATABANK_TYPE DT WHERE ID_DATABANK='$dbID[0]' AND D.ID_DBTYPE=DT.ID_DBTYPE ") or die "Couldn't prepare statement: " . $dbh->errstr;
		$dbDir="$promsPath{data}/banks/db_$dbID[0]";
		$dbFile="$dbDir/$fastaDBFile";
		
		#$dbFileID2=$dbFileID;
		open(FASTA,"<",$dbFile) or die ("open :$!");
		my ($irtFastaSeq,$match);
		while (<FASTA>) {
			if ($_=~/>/) {
				if ($_=~/iRT/ && $_!~/reverse/) {$match=1;}
				else{$match=0;}
			}
			else{
				if ($match) {
					$_=~s/\s//;
					chomp ($_);
					$irtFastaSeq.=$_;
				}
			}
		}
		
		###> Recovering RT file in database to put it in workDir
		(my $rtData,$rtFileName)=$dbh->selectrow_array("SELECT DATA,NAME FROM REFERENCE_RT WHERE ID_REFERENCE_RT='$rtFileID' ") or die "Couldn't prepare statement: " . $dbh->errstr;
		open(RTFILE,">>","$workDir/$rtFileName.txt") or die ("open: $!");
		my $xmlRT =XML::Simple-> new (KeepRoot=>1);
		my $xmlData = $xmlRT->XMLin($rtData);
		foreach my $pepInfo (@{$xmlData->{'PEPTIDE_DATA'}->{'PEPTIDE'}}) {
			if (! $pepInfo->{'excluded'}) {
				my $pepiRTSeq=$pepInfo->{'sequence'};
				if ($irtFastaSeq!~/$pepiRTSeq/) {
					print "<SCRIPT LANGUAGE=\"JavaScript\">alert('Wrong iRT file selected. Select an other file.'); window.location=\"$promsPath{cgi}/$script\";</SCRIPT>";
					exit;
				}
				print RTFILE ($pepInfo->{'sequence'},"\t",$pepInfo->{'iRT'},"\n");
			}
			$now=strftime("%s",localtime); # in sec
			$waitTime=strftime("%Hh %Mm %Ss",localtime($now-$startTime-3600));
			$status="Run time : $waitTime ";
			print "<SCRIPT LANGUAGE=\"JavaScript\">$divID.innerHTML=\"\";$divID.innerHTML='$status';</SCRIPT>";
		}
		
		
		
		###>Fetching parsing rules
		my ($parseRules,$identType)=$dbh->selectrow_array("SELECT PARSE_RULES,DEF_IDENT_TYPE FROM DATABANK_TYPE DT,DATABANK D WHERE D.ID_DBTYPE=DT.ID_DBTYPE AND ID_DATABANK=$dbID[0]");
		my @rules=split(',:,',$parseRules);
		my ($idRule)=($rules[0]=~/ID=(.+)/);
		if ($action ne 'addDB2') {
			print "<BR>Uploading data<BR>";
			########################################################################################################
			###> Recovering files names selected by the user, upload input files and move them in work directory<###
			########################################################################################################
			###> From input 
			for (my $i=1; $i<=3; $i++){
				if (param("inputfiles$i")){
					foreach my $name (param("inputfiles$i")){
						if($name=~m/\.dat$/){push @inputFilesDat,$name;}
						if($name=~m/\.tandem\./){push @inputFilesTandem,$name;}
						if($name=~m/\.mzXML$/){push @inputFilesMZXML,$name;}
						if($name=~m/\.pep.xml$/ && $name!~/\.tandem\.pep\.xml$/){push @inputFilesPEPXML,$name;}
						if($name=~m/\.xml$/ && $name!~/\.tandem\.pep\.xml$/ && $name!~/\.pep\.xml$/) {push @inputFilesXML,$name;}
						
						##> allocates a temporary address at input files
						my $newFile=tmpFileName($name);
						my $NewFileDir="$workDir/$name";
						
						##>move them in work directory
						copy($newFile,$NewFileDir);
						
						push @inputFiles,$name unless $name=~/\.xml$/ && $name!~/\.tandem\.pep\.xml$/ && $name!~/\.pep\.xml$/; 
						my $taille=`stat -c "%s" $workDir/$name`;
						if ($taille == 0 ) {
							print "<SCRIPT LANGUAGE=\"JavaScript\">alert('Error during importation of upload files. Please retry.'); window.location=\"$promsPath{cgi}/$script\";</SCRIPT>";
							exit;
						}
						push(@taille,int($taille));
						
						$now=strftime("%s",localtime); # in sec
						$waitTime=strftime("%Hh %Mm %Ss",localtime($now-$startTime-3600));
						$status="Run time : $waitTime ";
						print "<SCRIPT LANGUAGE=\"JavaScript\">$divID.innerHTML=\"\";$divID.innerHTML='$status';</SCRIPT>";
					}
				}
			}
			###> From server directory
			if($serverDir){
				foreach my $file (glob ("$serverDir/*")){
					(my $name=$file)=~s/.*\///;
					next unless $name=~/[\.dat|\.xml|\.mzXML]/;
					if($name=~m/\.dat$/){push @inputFilesDat,$name;}
					if($name=~m/\.tandem\./){push @inputFilesTandem,$name;}
					if($name=~m/\.mzXML$/){push @inputFilesMZXML,$name;}
					if($name=~m/\.pep.xml$/ && $name!~/\.tandem\.pep\.xml$/){push @inputFilesPEPXML,$name;}
					if($name=~m/\.xml$/ && $name!~/\.tandem\.pep\.xml$/ && $name!~/\.pep\.xml$/) {push @inputFilesXML,$name;}
					
					##> allocates a temporary address at input files
					my $NewFileDir="$workDir/$name";
					
					##>move them in work directory
					copy("$serverDir/$name",$NewFileDir);
					
					push @inputFiles,$name unless $name=~/\.xml$/ && $name!~/\.tandem\.pep\.xml$/ && $name!~/\.pep\.xml$/; 
					my $taille=`stat -c "%s" $workDir/$name`;
					if ($taille == 0 ) {
						print "<SCRIPT LANGUAGE=\"JavaScript\">alert('Error during importation of upload files. Please retry.'); window.location=\"$promsPath{cgi}/$script\";</SCRIPT>";
						exit;
					}
					push(@taille,int($taille));
					
					$now=strftime("%s",localtime); # in sec
					$waitTime=strftime("%Hh %Mm %Ss",localtime($now-$startTime-3600));
					$status="Run time : $waitTime ";
					print "<SCRIPT LANGUAGE=\"JavaScript\">$divID.innerHTML=\"\";$divID.innerHTML='$status';</SCRIPT>";
				}
			}
			###>From projects
			my @analysisFiles;
			if (param('analysisFiles')) {
				foreach my $sample (param('analysisFiles')){
					system "ln -s $promsPath{valid}/ana_$sample $workDir";
					my @sample=split ("/",$sample);
					my $sampleName=$sample[1];
					push @inputFilesDat, "$sampleName";
					push @inputFiles, "$sampleName";
					push @analysisFiles,$sample;
				}
			}
			###> From archive
			if (param('archfiles')){
				###> extracting files from archive
				my $uplFile=(split(/[\\\/]/,param('archfiles')))[-1]; 
				my $newFile="$workDir/$uplFile";
				copy(tmpFileName(upload("archfiles")),$newFile);
				###>Inflating file
				if ($newFile =~/\.(gz|zip)\Z/) {
					print "<BR>Extracting files...<BR>";
					my $type=$1;
					if ($type eq 'gz') {
						my $totFile=`tar -tf $newFile | wc -l`;
						$totFile=~s/\s//;
						system "tar -zvxf $newFile -C $workDir --suffix=gz --index-file=$outputFile &" ;
						my ($wait,$nbFile,$lastNbLineLastFile)=(1,0,0);
						my $lastFile='';
						print "<DIV id=\"archDIV\"></DIV>";
						my $archDIV="document.getElementById('archDIV')";
						while($wait==1){
							sleep 1;
							my $file= `tail -1 $outputFile`;
							$file=~s/\s//;
							my $nbFile=`wc -l <$outputFile`;
							$nbFile=~s/\s//;
							if ($file ne $lastFile){
								print "<SCRIPT LANGUAGE=\"JavaScript\">$archDIV.innerHTML=\"\";$archDIV.innerHTML='$nbFile/$totFile : $file';</SCRIPT>";
							}
							elsif($nbFile==$totFile){
								my $nbLineLastFile=`wc -l <$workDir/$file`;	
								if($nbLineLastFile==$lastNbLineLastFile){
									$wait=0;
								}
								$lastNbLineLastFile=$nbLineLastFile;
							}
							$lastFile=$file;
							
							###>print wait time
							$now=strftime("%s",localtime); # in sec
							$waitTime=strftime("%Hh %Mm %Ss",localtime($now-$startTime-3600));
							$status="Run time : $waitTime ";
							print "<SCRIPT LANGUAGE=\"JavaScript\">$divID.innerHTML=\"\";$divID.innerHTML='$status';</SCRIPT>";
						}
						print "<SCRIPT LANGUAGE=\"JavaScript\">$archDIV.innerHTML=\"\";</SCRIPT>";
						print "<BR/>";
					}
					elsif ($type eq 'zip') {
						print "<BR><BR>\n";
						&promsMod::unzipArchive($newFile,$workDir,{mode=>'verbose',txtBefore=>'<B>&nbsp;&nbsp;-',txtAfter=>"</B><BR/>\n"});
					}
				}
				###>Deleting zip file
				unlink $newFile;
				
				###> recovering files names
				opendir (DIR,$workDir);
				while (my $file=readdir(DIR)){
					next if $file=~/^\.+$/; #  skip '.' & '..' directories
					next if $file=~/\.txt$/;
					if($file=~/\.dat$/){push @inputFilesDat,$file; }
					if($file=~/\.tandem\./){push @inputFilesTandem,$file;}
					if($file=~/\.mzXML$/){push @inputFilesMZXML,$file;}
					if($file=~/\.pep.xml$/ && $file!~/\.tandem\.pep\.xml$/){push @inputFilesPEPXML,$file;}
					if($file=~/\.xml$/ && $file!~/\.tandem\.pep\.xml$/ && $file!~/\.pep\.xml$/) {push @inputFilesXML,$file;}
					
					push @inputFiles,$file unless $file=~/\.xml$/ && $file!~/\.tandem\.pep\.xml$/ && $file!~/\.pep\.xml$/; 
					my $taille=`stat -c "%s" $workDir/$file`;
					if ($taille == 0 ) {
						print "<SCRIPT LANGUAGE=\"JavaScript\">alert('Error during importation of upload files. Please retry.'); window.location=\"$promsPath{cgi}/$script\";</SCRIPT>";
						exit;
					}
					push(@taille,int($taille));
				}
				if (scalar @inputFilesMZXML == 0){
					print "<SCRIPT LANGUAGE=\"JavaScript\">alert('No .mzXML files selected.'); window.location=\"$promsPath{cgi}/$script\";</SCRIPT>";
					exit;
				}
				
				
				####> check that for each search file there is a mzXML associated file
				if(@inputFilesDat){
					foreach my $datFile (@inputFilesDat){
						(my $datFileName=$datFile)=~s/.dat$//;
						my $msxmlRequired=$datFileName.'.mzXML';
						unless (first {$_ eq $msxmlRequired} @inputFilesMZXML) {
							print "<SCRIPT LANGUAGE=\"JavaScript\">alert('The names of the .dat and .mzXML selected files are not the same. Select other files or change the names.'); window.location=\"$promsPath{cgi}/$script\";</SCRIPT>";
							exit;
						}
					}
					if (scalar @inputFilesDat != scalar @inputFilesMZXML){
						print "<SCRIPT LANGUAGE=\"JavaScript\">alert('The number of the .dat and .mzXML files are not the same.'); window.location=\"$promsPath{cgi}/$script\";</SCRIPT>";
						exit;
					}	
				}
				if(@inputFilesTandem){
					foreach my $tandemFile (@inputFilesTandem){
						(my $tandemFileName=$tandemFile)=~s/.tandem.pep.xml$//;
						my $msxmlRequired=$tandemFileName.'.mzXML';
						unless (first {$_ eq $msxmlRequired} @inputFilesMZXML) {
							print "<SCRIPT LANGUAGE=\"JavaScript\">alert('The names of the .tandem and .mzXML selected files are not the same. Select other files or change the names.'); window.location=\"$promsPath{cgi}/$script\";</SCRIPT>";
							exit;
						}
					}
					if (scalar @inputFilesTandem != scalar @inputFilesMZXML){
						print "<SCRIPT LANGUAGE=\"JavaScript\">alert('The number of the .tandem and .mzXML files are not the same.'); window.location=\"$promsPath{cgi}/$script\";</SCRIPT>";
						exit;
					}
				}
				if(@inputFilesPEPXML){
					foreach my $sequestFile (@inputFilesPEPXML){
						(my $sequestFileName=$sequestFile)=~s/.pep.xml$//;
						my $msxmlRequired=$sequestFileName.'.mzXML';
						unless (first {$_ eq $msxmlRequired} @inputFilesMZXML) {
							print "<HTML><BODY><SCRIPT LANGUAGE=\"JavaScript\">alert('The names of the .pep.xml and .mzXML selected files are not the same. Select other files or change the names.'); window.location=\"$promsPath{cgi}/$script\";</SCRIPT></BODY></HTML>";
							exit;
						}
					}
					if (scalar @inputFilesPEPXML != scalar @inputFilesMZXML){
						print "<HTML><BODY><SCRIPT LANGUAGE=\"JavaScript\">alert('The number of the .pep.xml and .mzXML files are not the same.'); window.location=\"$promsPath{cgi}/$script\";</SCRIPT></BODY></HTML>";
						exit;
					}
				}
			}
			###> from shared directory
			if (param('sharedDirFiles')){
				foreach my $sharedFile (param('sharedDirFiles')){
					my (undef,$path,$fileName)=splitpath($sharedFile);
					if($fileName=~/\.dat$/){push @inputFilesDat,$fileName; }
					if($fileName=~/\.tandem\./){push @inputFilesTandem,$fileName;}
					if($fileName=~/\.mzXML$/){push @inputFilesMZXML,$fileName;}
					if($fileName=~/\.pep.xml$/ && $fileName!~/\.tandem\.pep\.xml$/){push @inputFilesPEPXML,$fileName;}
					if($fileName=~/\.xml$/ && $fileName!~/\.tandem\.pep\.xml$/ && $fileName!~/\.pep\.xml$/) {push @inputFilesXML,$fileName;}
					push @inputFiles,$fileName unless $fileName=~/\.xml$/ && $fileName!~/\.tandem\.pep\.xml$/ && $fileName!~/\.pep\.xml$/; 
					
					##>move shared file in work directory
					my $NewFileDir="$workDir/$fileName";
					copy("$promsPath{shared}/$sharedFile",$NewFileDir);
					
					my $taille=`stat -c "%s" $NewFileDir`;
					if ($taille == 0 ) {
						print "<SCRIPT LANGUAGE=\"JavaScript\">alert('Error during importation of upload files. Please retry.'); window.location=\"$promsPath{cgi}/$script\";</SCRIPT>";
						exit;
					}
					push(@taille,int($taille));
					
					###>wait time
					$now=strftime("%s",localtime); # in sec
					$waitTime=strftime("%Hh %Mm %Ss",localtime($now-$startTime-3600));
					$status="Run time : $waitTime ";
					print "<SCRIPT LANGUAGE=\"JavaScript\">$divID.innerHTML=\"\";$divID.innerHTML='$status';</SCRIPT>";
					
				}
			}
			
			####################################################
			###> Converting .xml files into .tandem.pep.xml <###
			####################################################
			if (scalar @inputFilesXML != 0) {
				foreach my $xmlFile (@inputFilesXML){
					(my $xmlFileName=$xmlFile)=~s/\.xml$//;
					open(XMLFILE,"<","$workDir/$xmlFile") or die ("open : $!");
					open(OUT,">","$workDir/Mod$xmlFile") or die ("open : $!");
					while (<XMLFILE>) {
						if ($_=~/protein, cleavage N-terminal mass change/) {
							(my $line=$_)=~s/protein, cleavage N-terminal mass change">.*<\/note>/protein, cleavage N-terminal mass change"><\/note>/;
							print OUT $line;
						} elsif($_=~/protein, cleavage C-terminal mass change/){
							(my $line=$_)=~s/protein, cleavage C-terminal mass change">.*<\/note>/protein, cleavage C-terminal mass change"><\/note>/;
							print OUT $line;
						} elsif($_=~/spectrum, path/){
							(my $line=$_)=~s/spectrum, path\">.*$xmlFileName\.mzXML/spectrum, path\">$workDir\/$xmlFileName\.mzXML/;
							print OUT $line;
						} else{ print OUT $_; }
					}
					close XMLFILE;
					close OUT;
					
					move("$workDir/Mod$xmlFile","$workDir/$xmlFileName.tandem");
					system "$tppPath/Tandem2XML $workDir/$xmlFileName.tandem $workDir/$xmlFileName.tandem.pep.xml";
					push @inputFilesTandem, "$xmlFileName.tandem.pep.xml";
					push @inputFiles,"$xmlFileName.tandem.pep.xml";
				}	
			}
		   


			#check compatibility of identifier types of databank and input files
			my $fileIdentType;
			my $matchIdent;
			foreach my $file (@inputFiles){
				$matchIdent=0;
				if (not defined ($identType)) {
					my $fileIdent;
					if($file=~m/\.dat$/){
						$fileIdent=`grep -m 1 -P "^q\\d+_p\\d+=\\d+," $workDir/$file | cut -d\\" -f2`;
						chomp $fileIdent;
						next unless $fileIdent;
						my $testDB=`grep -m 1 $fileIdent $dbFile`;
						chomp $testDB;
						if ($testDB) {
							$matchIdent=1;
							last;
						}
					}
					elsif($file=~m/\.tandem\./ || $file=~m/\.pep\.xml/ ){
						my $res=`grep -m 1 "<search_hit " $workDir/$file`; # ' | cut -d\\" -f10' is not safe (XML attributes not always in same order)
						($fileIdent)=($res=~/protein="([^"]+)/);
						chomp $fileIdent;
						next unless $fileIdent;
						my $testDB=`grep -m 1 $fileIdent $dbFile`;
						chomp $testDB;
						if ($testDB) {
							$matchIdent=1;
							last;
						}
					}
					else{next;}
				}
				else{
					my $fileIdent;
					if($file=~m/\.dat$/){$fileIdent=`grep -m 1 -P "^q\\d+_p\\d+=\\d+," $workDir/$file | cut -d\\" -f2`;}
					elsif($file=~m/\.tandem\./ || $file=~m/\.pep\.xml/){
						my $res=`grep -m 1 "<search_hit " $workDir/$file`; # ' | cut -d\\" -f10' is not safe (XML attributes not always in same order)
						($fileIdent)=($res=~/protein="([^"]+)/);
					}
					else{next;}
					chomp $fileIdent;
					next unless $fileIdent;
					$fileIdent=~s/$decoyTag//;
					if ($fileIdent=~/([A-Z][\w\.-]+\|[^\s\|]+_[^\s\|]+)/){
						$fileIdentType="UNIPROT_ALL";
						if($identType eq 'UNIPROT_ALL' || $identType eq 'NCBI_ALL'){
							$matchIdent=1;
							last;
						}
					}
					elsif($fileIdent=~/([^\s\|]+_[^\s\|]+)/){
						$fileIdentType="UNIPROT_ID";
						if ($identType eq 'UNIPROT_ID') {
							$matchIdent=1;
							last;
						}
					}
					elsif($fileIdent=~/([A-Z][\w\.-]+)/){
						$fileIdentType="UNIPROT_ACCESSION";
						if ($identType eq 'UNIPROT_ACCESSION') {
							$matchIdent=1;
							last;
						}
					}
					elsif($fileIdent=~/IPI:([^| .]*)/){
						$fileIdentType="IPI_ACCESSION";
						if ($identType eq 'IPI_ACCESSION') {
							$matchIdent=1;
							last;
						}
					}
					elsif($fileIdent=~/(gi\|\d+\.?\d*)/){
						$fileIdentType="GI_ACCESSION";
						if ($identType eq 'GI_ACCESSION') {
							$matchIdent=1;
							last;
						}
					}
				}
				$now=strftime("%s",localtime); # in sec
				$waitTime=strftime("%Hh %Mm %Ss",localtime($now-$startTime-3600));
				$status="Run time : $waitTime ";
				print "<SCRIPT LANGUAGE=\"JavaScript\">$divID.innerHTML=\"\";$divID.innerHTML='$status';</SCRIPT>";
			}
			$fileIdentType=($fileIdentType) ? $fileIdentType : '';
			if (!$matchIdent) {
				if ($libCons eq 'new') {
					print "<BR><BR><BR><B>***ERROR*** Identifiers from input files ($fileIdentType) and Databank ($identType) are different! <BR>Select an other Databank.</B><BR><BR><BR>";
					my $sthNewDBList=$dbh->prepare("SELECT D.NAME, D.ID_DATABANK,D.FASTA_FILE,D.DECOY_TAG,DT.NAME FROM DATABANK D,DATABANK_TYPE DT WHERE D.DECOY_TAG IS NOT NULL AND D.USE_STATUS='yes' AND DT.DEF_IDENT_TYPE=\"$fileIdentType\" AND D.ID_DBTYPE=DT.ID_DBTYPE  ") or die "Couldn't prepare statement: " . $dbh->errstr;
					$sthNewDBList->execute;
					if ($sthNewDBList->rows==0) {
						print "<INPUT type=\"button\" class=\"buttonadd\" value=\"Return to form.\" onclick=\"window.location='./editSwathLibrary.cgi?ACT=add'\">";
						exit;
					}
					else{
						my %dbHash;
						while (my($dbName,$dbID,$dbFileName,$decoyTag,$dbTypeName)=$sthNewDBList->fetchrow_array){
							next if ($dbFileName=~m/:/) ;
							next if ($decoyTag eq "No" || $decoyTag eq "");
							$dbHash{"$dbName&nbsp;\[$dbTypeName\]"}=$dbID;
						}
						$sthNewDBList->finish;
						print qq
						|<FORM method="POST" action="./editSwathLibrary.cgi" name="parameters" enctype="multipart/form-data"><CENTER><TABLE bgcolor="$darkColor">
							<INPUT type="hidden" name="ACT" value="addDB2">
							<INPUT type="hidden" name="ID" id="ID" value="$libraryID">
							<INPUT type="hidden" name="workDir" value="$workDir">
							<INPUT type="hidden" name="refrtfile" value="$rtFileID">
							<INPUT type="hidden" name="libcons" value="$libCons">
							<INPUT type="hidden" name="liboption" value="$libOption">
							<INPUT type="hidden" name="libname" value="$libName">
							<INPUT type="hidden" name="descript" value="$des">
							<INPUT type="hidden" name="missedcleav" value="$missedCleavage">
							<INPUT type="hidden" name="fdr" value="$fdr">
							<INPUT type="hidden" name="instrument" value="$instrument">
							<INPUT type="hidden" name="fragmentation" value="$fragmentation">

							<TR><TH align=right>Databank : </TH>
							<TD bgcolor="$lightColor">
							<SELECT name ="dbfile" id="dbfile" required>
							<option value="">-= Select Databank =-</option>
						|;
						foreach my $dbName (sort {lc $a cmp lc $b} keys %dbHash){
							print "<option value=\"$dbHash{$dbName}\">$dbName</option>";
						}

						print qq
						|</SELECT></TD></TR><TR><TH colspan=2><BR><input type="submit" name="submit" value="Submit">
						<!-- CLEAR button -->
						&nbsp &nbsp &nbsp<INPUT type="reset" value="Clear" />
						<!-- CANCEL button -->
						&nbsp &nbsp &nbsp<INPUT type="button" value="Cancel" onclick="window.location='./editSwathLibrary.cgi?ACT=add'">
						</TH></TR></TABLE></CENTER></FORM>
						|;
						exit;
					}
					$sthNewDBList->finish;
				}
				elsif($libCons eq 'mergelib'){
					print "<BR><BR><BR><B>***ERROR*** Identifiers from input files ($fileIdentType) and SWATH library ($identType) are different! <BR>Select other Files.</B><BR><BR><BR>";
					print "<INPUT type=\"button\" class=\"buttonadd\" value=\"Return to form.\" onclick=\"window.location='./editSwathLibrary.cgi?ACT=add'\">";
					exit;
				}
				elsif($action eq 'update'){
					print "<BR><BR><BR><B>***ERROR*** Identifiers from input files ($fileIdentType) and SWATH library ($identType) are different! <BR>Select other Files.</B><BR><BR><BR>";
					print "<INPUT type=\"button\" class=\"buttonadd\" value=\"Return to list.\" onclick=\"window.location='./listSwathLibraries.cgi'\">";
					exit;
				}
			}
		}
		else{
			###> Recovering files names selected by the user
			my $i;
			foreach my $file (glob("$workDir/*")){
				my @fileDir=split(/\//,$file);
				my $fileName=$fileDir[-1];
				#if($file=~ m/\.dat$/){push @inputFilesDat,$fileName; }
				#if($file=~m/\.pep.xml$/ && $file!~/\.tandem\.pep\.xml$/){push @inputFilesPEPXML,$fileName;}
				#if($file=~m/\.tandem(\.|$)/){push @inputFilesTandem,$fileName;}
				#
				if($file=~/\.dat$/){push @inputFilesDat,$file; }
				if($file=~/\.tandem\./){push @inputFilesTandem,$file;}
				if($file=~/\.mzXML$/){push @inputFilesMZXML,$file;}
				if($file=~/\.pep.xml$/ && $file!~/\.tandem\.pep\.xml$/){push @inputFilesPEPXML,$file;}
				if($file=~/\.xml$/ && $file!~/\.tandem\.pep\.xml$/ && $file!~/\.pep\.xml$/) {push @inputFilesXML,$file;}
				
				
				my $taille=`stat -c "%s" $workDir/$fileName`;
				push(@taille,int($taille));
				push(@inputFiles,$fileName) unless $fileName=~/iRT/;
			}
		}
		@tailleSort=sort @taille;
	}
	elsif ($action eq "merge"){
			$libID1=param('libfile1');
			$libID2=param('libfile2');
			($libName,$des)=&promsMod::cleanParameters(param('libname'),param('descript'));
			$startTime=strftime("%s",localtime);
			$time=strftime("%Y_%m_%d_%H_%M_%S",localtime);
			
			$workDir="$promsPath{data}/tmp/Swath/$time";
			$libPath="$promsPath{swath_lib}";
			mkdir $libPath unless -e $libPath;
			mkdir $workDir unless -e $workDir;
			$outputFile=$workDir.'/sortie.txt';
			
			
			###> recovering Databank IDs
			my $sthDBID=$dbh->prepare("SELECT DT.DEF_IDENT_TYPE,D.ORGANISM,DS.ID_DATABANK FROM DATABANK_SWATHLIB DS, DATABANK D, DATABANK_TYPE DT WHERE ID_SWATH_LIB IN (?,?) AND D.ID_DATABANK=DS.ID_DATABANK AND D.ID_DBTYPE=DT.ID_DBTYPE");
			$sthDBID->execute($libID1,$libID2);
			while (my($identDB,$orgDB,$dbID)=$sthDBID->fetchrow_array) {
				push @dbID, "$dbID" if "@dbID"!~/$dbID/;
				$identifierDB=$identDB;
				$organismDB=$orgDB;
			}
			
			###>Creating library
			($libName1,$versionLib1,$rtFileID,$split,$instrument)=$dbh->selectrow_array("SELECT NAME,VERSION_NAME,ID_REFERENCE_RT,SPLIT,INSTRUMENT FROM SWATH_LIB WHERE ID_SWATH_LIB=$libID1") or die "Couldn't prepare statement: " . $dbh->errstr;
			($libName2,$versionLib2)=$dbh->selectrow_array("SELECT NAME,VERSION_NAME FROM SWATH_LIB WHERE ID_SWATH_LIB=$libID2") or die "Couldn't prepare statement: " . $dbh->errstr;
	}
	
	my $libID;
	if ($action ne 'update'){
		my $sthLibUpdate2=$dbh->prepare("INSERT INTO SWATH_LIB (ID_REFERENCE_RT,NAME,DES,VERSION_NAME,IDENTIFIER_TYPE,ORGANISM,USE_STATUS,INSTRUMENT) values (?,?,?,?,?,?,?,?)") or die "Couldn't prepare statement: " . $dbh->errstr;
		$sthLibUpdate2->execute($rtFileID,$libName,$des,'v1',$identifierDB,$organismDB,"err",$instrument);
		$sthLibUpdate2->finish;
		$libID=$dbh->selectrow_array("SELECT MAX(ID_SWATH_LIB) FROM SWATH_LIB");
		$dbh->commit;
	}
	else{
		$libID=$libIDFile;
	}
	
	
	my $fileInfo="$workDir/info_$libID.out";
	open(FILEINFO,">$fileInfo") || die "Error while opening $fileInfo";
	print FILEINFO "User=$userID\nSoftware=TPP";
	close FILEINFO;
   
		
	my $fileStat="$workDir/status_$libID.out";
	open(FILESTAT,">$fileStat") || die "Error while opening $fileStat";
	print FILESTAT "Started ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
	close FILESTAT;
	
	my $fileError="$workDir/status_$libID\_error.out";

	# Create new job to monitor
	$dbh->do("INSERT INTO JOB_HISTORY (ID_JOB, ID_USER, ID_PROJECT, TYPE, STATUS, FEATURES, SRC_PATH, LOG_PATH, ERROR_PATH, STARTED) VALUES('$time', '$userID', $projectID, 'Import [TPP]', 'Queued', 'SOFTWARE=TPP;ID_LIBRARY=$libID;ACTION=$action', '$workDir', '$fileStat', '$fileError', NOW())");
	$dbh->commit;
	
	print qq |
		<SCRIPT type="text/javascript">
			var monitorJobsWin=window.open("$promsPath{cgi}/monitorJobs.cgi?filterType=Import [TPP]&filterDateNumber=1&filterDateType=DAY&filterStatus=Queued&filterStatus=Running",'monitorJobsWindow','width=1200,height=500,scrollbars=yes,resizable=yes');
			monitorJobsWin.focus();
		</SCRIPT>
	|;
	
	my $childConvert=fork;
	unless($childConvert){
		#>Disconnecting from server
		open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
		open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
		open STDERR, '>/dev/null' or die "Can't open /dev/null: $!";
		
		# Add process PID to current job in DB
		my $dbh = &promsConfig::dbConnect('no_user');
		$dbh->do("UPDATE JOB_HISTORY SET ID_JOB_CLUSTER='L$$' WHERE ID_JOB='$time'");
		$dbh->commit;
		
		if ($action eq "add" || $action eq "update" || $action eq "addDB2") {
			open(FILESTAT,">>$fileStat");
			print FILESTAT "Fetching files.\n";
			close(FILESTAT);
			###>Modification of sequest files
			if(@inputFilesPEPXML){
				foreach my $sequestFile (@inputFilesPEPXML){
					(my $fileName=$sequestFile)=~s/\.pep\.xml//;
					open(INFILE,"<$workDir/$sequestFile");
					open(OUTFILE,">$workDir/$fileName.sequest.pep.xml");
					my $match;
					while (<INFILE>){
						if($_=~/<search_summary>/){
							$match=1;
						}elsif($_=~/<\/WorkflowMessages>/){
							$match=0;
							next;
						}
						if($_=~/<msms_run_summary base_name/){
							$_=~s/base_name=".*(\w|\))\.msf"/base_name="$workDir\/$fileName"/;
							$_=~s/raw_data=".msf"/raw_data=".mzXML"/;
							print OUTFILE $_;
						}
						elsif($_=~/<search_summary base_name/){
							$_=~s/base_name=".*\.msf"/base_name="$workDir\/$fileName.msf"/;
							$_=~s/search_engine="Sequest HT"/search_engine="SEQUEST"/;
							print OUTFILE $_;
						}
						elsif($_=~/<sample_enzyme name="/){
							$_=~s/<sample_enzyme name=".*"/<sample_enzyme name="Trypsin"/;
							print OUTFILE $_;
						}
						elsif($_=~/<search_database/){
							$_=~s/local_path=""/local_path="$dbFile"/;
							print OUTFILE $_;
						}
						elsif($_=~/<parameter name="FastaDatabase"/){
							$_=~s/value="$fastaDBFile"/value="$dbFile"/;
							print OUTFILE $_;
						}
						elsif(!$match){print OUTFILE $_;}
					}
					close INFILE;
					close OUTFILE;
					push @inputFilesSequest, "$fileName.sequest.pep.xml";
				}
			}
		
			####> recovering non consensus library 
			if ($action eq 'update' || $libCons eq 'mergelib'){
				opendir (DIR, "$libPath/SwLib_$libIDFile/");
				while (my $file = readdir (DIR)){
					if ($file =~/\.tar\.gz/){
						(my $folderName=$file)=~s/\.tar\.gz//;
						system "cd $libPath/SwLib_$libIDFile; tar -zxf $file;";
						system "rm $libPath/SwLib_$libIDFile/dbinfo_$folderName.txt;" if -e "$libPath/SwLib_$libIDFile/dbinfo_$folderName.txt";
						system "rm $libPath/SwLib_$libIDFile/filelist_$folderName.txt;" if -e "$libPath/SwLib_$libIDFile/filelist_$folderName.txt";
						system "rm $libPath/SwLib_$libIDFile/sortie_$folderName.txt;" if -e "$libPath/SwLib_$libIDFile/sortie_$folderName.txt";
						system "rm $libPath/SwLib_$libIDFile/script_$folderName.sh;" if -e "$libPath/SwLib_$libIDFile/script_$folderName.sh";
					}
				}
				close DIR;
				opendir (CONS,"$libPath/SwLib_$libIDFile/");
				while (my $file2 = readdir (CONS)){
					if ($file2=~/SpecLib_Lib\d+_v\d+/){
						move ("$libPath/SwLib_$libIDFile/$file2",$workDir);
					}
				}
				close CONS;
			}
			
			
			
			##############################################
			###> Convert .dat files in .pep.xml files <###
			##############################################
			my $nbDat=1;
			my @inputMascotFiles;
			my $convDIV;
			my $last='Upload';
			if (@inputFilesDat) {
				open(FILESTAT,">>$fileStat");
				print FILESTAT "Conversion of dat files to pep.xml.\n";
				close(FILESTAT);
				$last="Conversion";
			}
			
			my $maxHours=int(48);
			foreach my $datFile (@inputFilesDat){
				my $mascotFile=$workDir.'/'.$datFile;
				(my $datName=$datFile)=~s/.dat//;
				if ($clusterInfo{'on'}){
					my $dir=$workDir."/datConversion_".$nbDat;
					mkdir $dir unless -e $dir;
					my $outFile=$dir.'/sortie.txt';
					my $datBashName='datConv.sh';
					open (BASH, ">$dir/$datBashName");
					print BASH "#!/bin/bash\n";
					print BASH "$tppPath/Mascot2XML $mascotFile -D$dbFile -Etrypsin -notgz >>$outFile 2>&1;";
					close BASH;
					my $bashName='datConvertion.sh';
					my $clusterCommandString=$clusterInfo{'buildCommand'}->($dir,"$dir/$datBashName");
					#my $maxHours=int(48);
					my $datSize=`stat -c "%s" $mascotFile`;
					$maxMem=($datSize/1073741824)*3;
					$maxMem=($maxMem < 5) ? 5 : ($maxMem > 30) ? 30 : sprintf("%.0f",$maxMem);
					$maxMem.='Gb';
					
					my $bashFile="$dir/$bashName";
					open (BASHDAT,">$bashFile");
					print BASHDAT qq
|#!/bin/bash

##resources
#PBS -l mem=$maxMem
#PBS -l nodes=1:ppn=1
#PBS -l walltime=$maxHours:00:00
#PBS -q batch

##Information
#PBS -N Dat_conversion_$nbDat
#PBS -M marine.le-picard\@curie.fr
#PBS -m abe
#PBS -o $dir/PBS.txt 
#PBS -e $dir/PBSerror.txt

## Command
$clusterCommandString
|;
					close BASHDAT;
					
					###> Execute bash file
					system "chmod 775 $dir/$datBashName";
					
					my ($jobClusterID) = $clusterInfo{'sendToCluster'}->($bashFile);
					
					$dbh->do("UPDATE JOB_HISTORY SET ID_JOB_CLUSTER = CONCAT(ID_JOB_CLUSTER, ';C$jobClusterID') WHERE ID_JOB='$time'");
					$dbh->commit;
				}
				else{
					my $outFile=$workDir.'/sortie'.$nbDat.'.txt';
					system "$tppPath/Mascot2XML $mascotFile -D$dbFile -Etrypsin -notgz >>$outFile 2>&1";
					if (-s $outFile && `tail $outFile`=~/error/){
						my $error=`tail $outFile`.'\n';
						open(FILEERROR,">>$fileError");
						print FILEERROR $error;
						close(FILEERROR);
						exit;
					}
				}
				$nbDat++;
				sleep 30;
			}
			
			my ($count,$nbFiles,$PBSerror)=(0,0,'');
			if ($clusterInfo{'on'}){
				my $nbWhile=0;
				my $maxNbWhile=$maxHours*60*2;
				while ($nbFiles != scalar @inputFilesDat && !$PBSerror){
					$nbFiles=`ls $workDir/*.pep.xml | grep -v .tandem.pep.xml |wc -l `;
					for (my $i=1;$i<=scalar @inputFilesDat ;$i++){
						$PBSerror=$clusterInfo{'checkError'}->("$workDir/datConversion_$i/PBSerror.txt") if -s "$workDir/datConversion_$i/PBSerror.txt";
						#$PBSerror=`grep -v 'InfluxDB\\|getcwd' $workDir/datConversion_$i/PBSerror.txt` if -s "$workDir/datConversion_$i/PBSerror.txt";
						if (-s "$workDir/datConversion_$i/sortie.txt" && `tail $workDir/datConversion_$i/sortie.txt`=~/Error/){
							$PBSerror=`tail $workDir/datConversion_$i/sortie.txt`;
						}
						last if $PBSerror;
					}
					sleep 30;
					$nbWhile++;
					if ($nbWhile > $maxNbWhile) {
						open(FILEERROR,">>$fileError");
						print FILEERROR "Aborting: Dat files conversion is taking too long.\n";
						close(FILEERROR);
						exit;
					}
				}
				if ($PBSerror){
					open(FILEERROR,">>$fileError");
					print FILEERROR "Dat files conversion failed.<BR>$PBSerror";
					close(FILEERROR);
					exit;
				}
			}
			sleep 30;
			
			
			foreach my $datFile (@inputFilesDat){
				(my $datName=$datFile)=~s/\.dat+//;
				my $pepXMLFile=$workDir.'/'.$datName.'.pep.xml';
				my $mascotPepXMLFile=$workDir.'/'.$datName.'.mascot.pep.xml';
				if (-s $pepXMLFile){
					push(@inputMascotFiles,$mascotPepXMLFile);                        #retrieve mascot files's names for xinteract
					move($pepXMLFile,$mascotPepXMLFile);
				}
				else{
					open(FILEERROR,">>$fileError");
					print FILEERROR "Dat file ($datFile) conversion failed.\n";
					close(FILEERROR);
					exit;
				}
			}
			
			
			###############################################################
			###> Change destination addresses in .tandem.pep.xml files <###
			###############################################################
			my $qServerPath=quotemeta($tppPath);
			my $qFastaDbFile=quotemeta($fastaDBFile);
			my $qDbDir=quotemeta($dbDir);
			my $qWorkDir=quotemeta($workDir);
			my $modifDIV;
			if (@inputFilesTandem){
				open(FILESTAT,">>$fileStat");
				print FILESTAT "Modification of X! TANDEM files.\n";
				close(FILESTAT);
				$last="Modification";
			}
			my $nbTandem=1;
			#my $maxHours=int(48);
			foreach my $tandemPepFile (@inputFilesTandem){
				(my $fileBaseName=$tandemPepFile)=~s/\.tandem.+//;
				my $qFileBaseName=quotemeta($fileBaseName);
				my $tandemFile=$workDir.'/'.$tandemPepFile;
				my $dir=$workDir."/tandemModification_".$nbTandem;
				mkdir $dir unless -e $dir;
				my $outFile=$dir.'/sortie.txt';
				my $tandemBashName='tandemModif.sh';
				open (BASH, ">$dir/$tandemBashName");
				print BASH "#!/bin/bash \n";
				print BASH "sed -r -i -e 's/[^\"]+\\\/$qFileBaseName\.?\"\/$qWorkDir\\\/$qFileBaseName\"/g' $workDir/$tandemPepFile >>$outFile 2>&1; \n";
					 #Taxonomy file
				print BASH "sed -r -i -e 's/[^\"]+\\\/taxonomy\.xml\/$qServerPath\\\/taxonomy\.xml/g' $workDir/$tandemPepFile >>$outFile 2>&1; \n";
					 #Database file
				print BASH "sed -r -i -e 's/local_path=\".*\.fasta\"/local_path=\"$qDbDir\\\/$qFastaDbFile\"/g' $workDir/$tandemPepFile >>$outFile 2>&1; \n";
				print BASH "sed -r -i -e 's/<parameter name=\"list path, sequence source #1\" value=\".*\.fasta\"/<parameter name=\"list path, sequence source #1\" value=\"$qDbDir\\\/$qFastaDbFile\"/g' $workDir/$tandemPepFile >>$outFile 2>&1; \n";
				print BASH "echo END >>$workDir/END.txt;";
				close BASH;
				
				if ($clusterInfo{'on'}){
					my $bashName='tandemModification.sh';
					my $clusterCommandString=$clusterInfo{'buildCommand'}->($dir,"$dir/$tandemBashName");
					#my $maxHours=int(48);
					my $datSize=`stat -c "%s" $tandemFile`;
					$maxMem=($datSize/1073741824)*3;
					$maxMem=($maxMem < 5) ? 5 : ($maxMem > 20) ? 20 : sprintf("%.0f",$maxMem);
					$maxMem.='Gb';
					
					my $bashFile="$dir/$bashName";
					open (BASHTANDEM,">$bashFile");
					print BASHTANDEM qq
|#!/bin/bash

##resources
#PBS -l mem=$maxMem
#PBS -l nodes=1:ppn=1
#PBS -l walltime=$maxHours:00:00
#PBS -q batch

##Information
#PBS -N Tandem_modification_$nbTandem
#PBS -M marine.le-picard\@curie.fr
#PBS -m abe
#PBS -o $dir/PBS.txt 
#PBS -e $dir/PBSerror.txt

## Command
$clusterCommandString
echo END_Tandem_modification_$nbTandem
|;
					close BASHTANDEM;
					
					###> Execute bash file
					system "chmod 775 $dir/$tandemBashName";
					my ($jobClusterID) = $clusterInfo{'sendToCluster'}->($bashFile);
					
					$dbh->do("UPDATE JOB_HISTORY SET ID_JOB_CLUSTER = CONCAT(ID_JOB_CLUSTER, ';C$jobClusterID') WHERE ID_JOB='$time'");
					$dbh->commit;
				}
				else{
					system "bash $dir/$tandemBashName";
				}
				$nbTandem++;
				sleep 30;
			}
			
			####> wait loop for cluster job
			($PBSerror,$count)=('',0);
			my $checkFiles=0;
			if ($clusterInfo{'on'}){
				my $nbWhile=0;
				my $maxNbWhile=$maxHours*60*2;
				while ($checkFiles != scalar @inputFilesTandem && !$PBSerror){
					for (my $i=1;$i<=scalar @inputFilesTandem ;$i++){
						$PBSerror=$clusterInfo{'checkError'}->("$workDir/tandemModification_$i/PBSerror.txt") if -s "$workDir/tandemModification_$i/PBSerror.txt";
						#$PBSerror=`grep -v 'InfluxDB\\|getcwd' $workDir/tandemModification_$i/PBSerror.txt` if -s "$workDir/tandemModification_$i/PBSerror.txt";
						last if $PBSerror;
						my $file=`head $workDir/tandemModification_$i/PBS.txt` if -s "$workDir/tandemModification_$i/PBS.txt";
						next unless $file;
						$checkFiles++ if $file=~/END_Tandem_modification_$i/;
					}
					$checkFiles=0 unless $checkFiles==scalar @inputFilesTandem;
					
					sleep 30;
					$nbWhile++;
					if ($nbWhile > $maxNbWhile) {
						open(FILEERROR,">>$fileError");
						print FILEERROR "Aborting: Tandem files conversion is taking too long.\n";
						close(FILEERROR);
						exit;
					}
				}
				if ($PBSerror){
					open(FILEERROR,">>$fileError");
					print FILEERROR "Tandem files conversion failed.<BR>$PBSerror.";
					close(FILEERROR);
					exit;
				}
				system "rm $workDir/END.txt;" if -e "$workDir/END.txt";
			}
			
			
			my @tandemFiles=map({$workDir.'/'.$_} @inputFilesTandem);
			my @sequestFiles=map({$workDir.'/'.$_} @inputFilesSequest);
			
			##############################
			####>Bash script commands<####
			##############################
			###> Create new bash file to store system commands
			open(BASH,"+>","$workDir/script.sh");
			print BASH "#!/bin/bash \n";
		
		
		
			############################
			###> iProphet xinteract <###
			############################
			
				#Xtandem files
			if (@tandemFiles) {	#if there are xtandem input files
				print BASH "echo \"Xinteract on X! TANDEM files.\" >>$fileStat; \n";
				print BASH "$tppPath/xinteract -OARPd -d$decoyTag -N$workDir/interact.tandem.pep.xml @tandemFiles >>$outputFile 2>&1;  \n ";
			}
			
				#Mascot files
			if (@inputMascotFiles) {   #if there are mascot input files
				print BASH "echo \"Xinteract on MASCOT files.\" >>$fileStat; \n";
				print BASH "$tppPath/xinteract -OARPd -d$decoyTag -N$workDir/interact.mascot.pep.xml @inputMascotFiles >>$outputFile 2>&1;  \n";
			}
			
				#Sequest files
			if(@sequestFiles){
				print BASH "echo \"Xinteract on SEQUEST files.\" >>$fileStat; \n";
				print BASH "$tppPath/xinteract -OARPd -d$decoyTag -N$workDir/interact.sequest.pep.xml @sequestFiles >>$outputFile 2>&1;  \n";
			}
			
			############################
			###> InterProphetParser <###		
			############################
			
			if (@tandemFiles){
				if(@inputMascotFiles){
					if(@sequestFiles){		## Tandem, Mascot and Sequest Files
						print BASH qq
|nbmascot=\`find $workDir/interact.mascot.pep* -type f \| wc -l\`;
nbtandem=\`find $workDir/interact.tandem.pep* -type f \| wc -l\`;
nbsequest=\`find $workDir/interact.sequest.pep* -type f \| wc -l\`;
if [ -e $workDir/interact.mascot.pep.xml ] && [ -e $workDir/interact.tandem.pep.xml ] && [ -e $workDir/interact.sequest.pep.xml ] && [ "\$nbtandem" -gt "4" ] && [ "\$nbmascot" -gt "4" ] && [ "\$nbsequest" -gt "4" ]
then
	echo "iProphet on X! TANDEM, MASCOT and SEQUEST files." >>$fileStat;
	$tppPath/InterProphetParser DECOY=$decoyTag $workDir/interact.tandem.pep.xml $workDir/interact.mascot.pep.xml $workDir/interact.sequest.pep.xml $workDir/iProphet.pep.xml >>$outputFile  2>&1;
else
	if [ -e $workDir/interact.mascot.pep.xml ]
	then
		if [ -e $workDir/interact.tandem.pep.xml ]
		then 
			echo "Xinteract on SEQUEST files did not work." >>$workDir/ERROR.txt;
			exit -1;
		else
			if [ -e $workDir/interact.sequest.pep.xml ]
			then
				echo "Xinteract on X! TANDEM files did not work." >>$workDir/ERROR.txt;
				exit -1;
			else
				echo "Xinteract on X! TANDEM and SEQUEST files did not work." >>$workDir/ERROR.txt;
				exit -1;
			fi
		fi
	else
		if [ -e $workDir/interact.tandem.pep.xml ]
		then
			if [ -e $workDir/interact.sequest.pep.xml ]
			then
				echo "Xinteract on MASCOT files did not work." >>$workDir/ERROR.txt;
				exit -1;
			else
				echo "Xinteract on MASCOT and SEQUEST files did not work." >>$workDir/ERROR.txt;
				exit -1;
			fi
		else
			if [ -e $workDir/interact.sequest.pep.xml ]
			then
				echo "Xinteract on MASCOT and X! TANDEM files did not work." >>$workDir/ERROR.txt;
				exit -1;
			else
				echo "Xinteract on MASCOT, X! TANDEM and SEQUEST files did not work." >>$workDir/ERROR.txt;
				exit -1;
			fi
		fi
	fi
fi
|;
					}
					else{						## Tandem and Mascot Files
						print BASH qq 
|nbmascot=\`find $workDir/interact.mascot.pep* -type f \| wc -l\`;
nbtandem=\`find $workDir/interact.tandem.pep* -type f \| wc -l\`;
if [ -e $workDir/interact.mascot.pep.xml ] && [ -e $workDir/interact.tandem.pep.xml ] && [ "\$nbtandem" -gt "4" ] && [ "\$nbmascot" -gt "4" ]
then
	echo "iProphet on X! TANDEM and MASCOT files." >>$fileStat;
	$tppPath/InterProphetParser DECOY=$decoyTag $workDir/interact.tandem.pep.xml $workDir/interact.mascot.pep.xml $workDir/iProphet.pep.xml >>$outputFile  2>&1;
else
	if [ -e $workDir/interact.mascot.pep.xml ]
	then
		echo "Xinteract on X! TANDEM files did not work." >>$workDir/ERROR.txt;
		exit -1;
	else
		if [ -e $workDir/interact.tandem.pep.xml ]
		then
			echo "Xinteract on MASCOT files did not work." >>$workDir/ERROR.txt;
			exit -1;
		else
			echo "Xinteract on MASCOT and X! TANDEM files did not work." >>$workDir/ERROR.txt;
			exit -1;
		fi
	fi
fi
|;
					}
				}
				elsif(@sequestFiles){		## Tandem and Sequest Files
					print BASH qq
|nbsequest=\`find $workDir/interact.sequest.pep* -type f \| wc -l\`;
nbtandem=\`find $workDir/interact.tandem.pep* -type f \| wc -l\`;
if [ -e $workDir/interact.sequest.pep.xml ] && [ -e $workDir/interact.tandem.pep.xml ] && [ "\$nbtandem" -gt "4" ] && [ "\$nbsequest" -gt "4" ]
then
	echo "iProphet on X! TANDEM and SEQUEST files." >>$fileStat;
	$tppPath/InterProphetParser DECOY=$decoyTag $workDir/interact.tandem.pep.xml $workDir/interact.sequest.pep.xml $workDir/iProphet.pep.xml >>$outputFile  2>&1;
else
	if [ -e $workDir/interact.sequest.pep.xml ]
	then
		echo "Xinteract on X! TANDEM files did not work." >>$workDir/ERROR.txt;
		exit -1;
	else
		if [ -e $workDir/interact.tandem.pep.xml ]
		then
			echo "Xinteract on SEQUEST files did not work." >>$workDir/ERROR.txt;
			exit -1;
		else
			echo "Xinteract on SEQUEST and X! TANDEM files did not work." >>$workDir/ERROR.txt;
			exit -1;
		fi
	fi
fi
|;
				}
				else{							## Tandem Files
					print BASH qq
|nbtandem=\`find $workDir/interact.tandem.pep* -type f \| wc -l\`;
if [ -e $workDir/interact.tandem.pep.xml ] && [ "\$nbtandem" -gt "4" ]
then
	echo "iProphet on X! TANDEM files." >>$fileStat;
	$tppPath/InterProphetParser DECOY=$decoyTag $workDir/interact.tandem.pep.xml $workDir/iProphet.pep.xml >>$outputFile  2>&1;
else
	echo "Xinteract on X! TANDEM files did not work." >>$workDir/ERROR.txt;
	exit -1;
fi
|;
				}
			}
			elsif(@inputMascotFiles){ 		
				if(@sequestFiles){			## Mascot and Sequest Files
					print BASH qq
|nbsequest=\`find $workDir/interact.sequest.pep* -type f \| wc -l\`;
nbmascot=\`find $workDir/interact.mascot.pep* -type f \| wc -l\`;
if [ -e $workDir/interact.sequest.pep.xml ] && [ -e $workDir/interact.mascot.pep.xml ] && [ "\$nbmascot" -gt "4" ] && [ "\$nbsequest" -gt "4" ]
then
	echo "iProphet on MASCOT and SEQUEST files." >>$fileStat;
	$tppPath/InterProphetParser DECOY=$decoyTag $workDir/interact.mascot.pep.xml $workDir/interact.sequest.pep.xml iProphet.pep.xml >>$outputFile  2>&1;
else
	if [ -e $workDir/interact.sequest.pep.xml ]
	then
		echo "Xinteract on MASCOT files did not work." >>$workDir/ERROR.txt;
		exit -1;
	else
		if [ -e $workDir/interact.mascot.pep.xml ]
		then
			echo "Xinteract on SEQUEST files did not work." >>$workDir/ERROR.txt;
			exit -1;
		else
			echo "Xinteract on SEQUEST and MASCOT files did not work." >>$workDir/ERROR.txt;
			exit -1;
		fi
	fi
fi
|;
				}
				else{							## Mascot Files
					print BASH qq
|nbmascot=\`find $workDir/interact.mascot.pep* -type f \| wc -l\`;
if [ -e $workDir/interact.mascot.pep.xml ] && [ "\$nbmascot" -gt "4" ]
then
	echo "iProphet on MASCOT files." >>$fileStat;
	$tppPath/InterProphetParser DECOY=$decoyTag $workDir/interact.mascot.pep.xml $workDir/iProphet.pep.xml >>$outputFile  2>&1;
else
	echo "Xinteract on MASCOT files did not work." >>$workDir/ERROR.txt;
	exit -1;
fi
|;
				}
			}
			else{								## Sequest Files
				print BASH qq
|nbsequest=\`find $workDir/interact.sequest.pep* -type f \| wc -l\`;
if [ -e $workDir/interact.sequest.pep.xml ] && [ "\$nbsequest" -gt "4" ]
then
	echo "iProphet on SEQUEST files." >>$fileStat;
	$tppPath/InterProphetParser DECOY=$decoyTag $workDir/interact.sequest.pep.xml $workDir/iProphet.pep.xml >>$outputFile  2>&1;
else
	echo "Xinteract on SEQUEST files did not work." >>$workDir/ERROR.txt;
	exit -1;
fi
|;
			}

			
if(param("PTMProphet")) {
	my $modifications = (param("PTMProphetModifs")) ? join(',', promsMod::cleanParameters(param("PTMProphetModifs"))) : "STY:79.9663,M:15.9949,C:57.0214";
	my $nIons = ($fragmentation eq 'ETD') ? 'c' : 'b';
	my $cIons = ($fragmentation eq 'ETD') ? 'z' : 'y';
	my $ms1PPMTol = (param('PTMProphetPPMTol')) ? (param('PTMProphetPPMTol')) : 1;
	my $fragPPMTol = (param('PTMProphetFragPPMTol')) ? (param('PTMProphetFragPPMTol')) : 10;
	my $EMmodel = (param('PTMProphetEMModel')) ? param('PTMProphetEMModel') : 2;
	my $minProb = (param('PTMProphetMinProb')) ? param('PTMProphetMinProb') : 0.9;
	
 			##########################
			###> PTMProphetParser <###
			##########################
			print BASH qq
|if [ -e $workDir/iProphet.pep.xml ]
then
	echo "PTMProphetParser: on iProphet file." >>$fileStat;
	cd $workDir;
	$tppPath/PTMProphetParser VERBOSE EM=$EMmodel MINPROB=$minProb PPMTOL=$ms1PPMTol FRAGPPMTOL=$fragPPMTol MAXTHREADS=$MAX_NB_THREAD NIONS=$nIons CIONS=$cIons $modifications $workDir/iProphet.pep.xml >> $outputFile  2>&1;
else
	echo "PTMProphetParser did not work." >>$workDir/ERROR.txt;
	exit -1;
fi
|;
}
	
			
			##############
			###> Mayu <###
			##############
			print BASH qq
|if [ -e $workDir/iProphet.pep.xml ]
then
	echo "Mayu on iProphet file." >>$fileStat;
	cd $workDir;
	$tppPath/Mayu.pl -A $workDir/iProphet.pep.xml -C $dbFile -E $decoyTag -G 0.01 -H 51 -I $missedCleavage -P $fdrType=$fdr:t >>$outputFile  2>&1;
else
	echo "InterProphetParser did not work." >>$workDir/ERROR.txt;
	exit -1;
fi
|;
			
			print BASH qq
|wait;
sleep 30;
if [ -e $workDir/*FDR0.01_t_1.07.csv ]
then
	cd $workDir;
	sort -t, -k5n $workDir/*FDR0.01_t_1.07.csv >>$workDir/tabFDR1.csv;
	wait;
	fdrMin=`sed '2q;d' $workDir/tabFDR1.csv \| cut -f5 -d,`;
	if [ \$fdrMin ]
	then
		echo "\$fdrMin" >>$workDir/FDRscore.txt;
		echo "The minimum score is : \$fdrMin." >>$fileStat;
	else
		echo "%Mayu : Not enough input files.\n" >>$workDir/ERROR.txt;
		exit -1;
	fi
|;
	
	
			###################
			###> Spectrast <###
			###################
			
			print BASH qq
|	sleep 30;
	echo "Creating spectras library." >>$fileStat;
	$tppPath/spectrast -cN$workDir/SpecLib -cI$fragmentation -cf \"Protein!~ $decoyTag\" -cP\$fdrMin -c_IRT$rtFileName.txt -c_IRR $workDir/iProphet.pep.xml >>$outputFile 2>&1;
|;
	
			###> Create consensus library
			if ($libOption eq 'unsplit'){ #unsplit mode
				print BASH "echo \"Creating consensus library with unsplit mode in progress.\" >>$fileStat;\n";
				if ($libCons eq 'mergelib' || $action eq 'update') {
					print BASH "cd $workDir; $tppPath/spectrast -cN$workDir/SpecLib_cons -cI$fragmentation -cAC -cM $workDir/SpecLib*.splib >>$outputFile 2>&1;\n";
				}
				elsif ($libCons eq 'new' ){
					print BASH "cd $workDir; $tppPath/spectrast -cN$workDir/SpecLib_cons -cI$fragmentation -cAC -cM $workDir/SpecLib.splib >>$outputFile 2>&1;\n";
				}
				$split=0;
			}
			else{	#split mode
				print BASH qq
|	echo "Creating consensus library with split mode in progress." >>$fileStat;
	$pythonPath/spectrast_cluster.py -d 2 $workDir/SpecLib.sptxt >>$outputFile 2>&1;
	for file in $workDir/SpecLib_*.sptxt; do $tppPath/spectrast -cNsplit-\${file%%.*} -cI$fragmentation \$file >>$outputFile 2>&1;  done 
	for file in $workDir/split-SpecLib_*.splib; do $tppPath/spectrast -cNcons-\${file%%.*} -cI$fragmentation -cAC \$file >>$outputFile 2>&1;  done 
|;
				if ($libCons eq 'mergelib'|| $action eq 'update') {
					print BASH "grep -hUv '###' $workDir/cons-split-SpecLib_*.sptxt $libPath/SwLib_$libIDFile/$libFileName.sptxt >>$workDir/SpecLib_cons_concat.sptxt;\n";
				}
				else{
					print BASH "grep -hUv '###' $workDir/cons-split-SpecLib_*.sptxt >>$workDir/SpecLib_cons_concat.sptxt;\n";
				}
				$split=1;
			}
			print BASH qq
|else
	echo "Mayu did not work, missing file : FDR0.01_t_1.07.csv." >>$workDir/ERROR.txt;
	exit -1;
fi
echo "END" >>$workDir/END.txt;
|;
			close BASH;
			
			###########################
			###> Execute bash file <###
			###########################
			
			if ($clusterInfo{'on'}) {
				my $clusterCommandString=$clusterInfo{'buildCommand'}->($workDir,"$workDir/script.sh");
				my $maxHours=int(48);
				$maxMem=($tailleSort[-1]/1073741824)*10;
				$maxMem*=2 if $clusterInfo{'name'} eq 'CentOS';
				$maxMem=($maxMem < 30) ? 30 : ($maxMem > 100) ? 100: sprintf("%.0f",$maxMem) ;
				$maxMem.='Gb';
				
				my $bashFile="$workDir/createLib.sh";
				open(BASH2,">$bashFile");
				print BASH2 qq
	|#!/bin/bash
	##resources
	#PBS -l mem=$maxMem
	#PBS -l nodes=1:ppn=1
	#PBS -l walltime=$maxHours:00:00
	#PBS -q batch
	
	##Information
	#PBS -N Swath_Lib_$time
	#PBS -M marine.le-picard\@curie.fr
	#PBS -m abe
	#PBS -o $workDir/PBS.txt
	#PBS -e $workDir/PBSerror.txt
	
	## Command
	$clusterCommandString
	echo End_Swath_Lib_$time >>$workDir/END.txt
	|;
				close BASH2;
				
				system "chmod 775 $workDir/script.sh";
				###> Execute bash file
				my ($jobClusterID) = $clusterInfo{'sendToCluster'}->($bashFile);
				
				# Add to DB
				$dbh->do("UPDATE JOB_HISTORY SET ID_JOB_CLUSTER = CONCAT(ID_JOB_CLUSTER, ';C$jobClusterID') WHERE ID_JOB='$time'");
				$dbh->commit;
				
				###>Waiting for job to run
				my $pbsError;
				my $j=1;
				my $nbWhile=0;
				my $maxNbWhile=$maxHours*60*2;
				while (!$pbsError  && !-s "$workDir/ERROR.txt" && !-s "$workDir/END.txt") {
					sleep 30;
					$nbWhile++;
					if ($nbWhile > $maxNbWhile) {
						open(FILEERROR,">>$fileError");
						print FILEERROR "Aborting: File processing is taking too long.";
						close(FILEERROR);
						exit;
					}
					$pbsError=$clusterInfo{'checkError'}->("$workDir/PBSerror.txt") if -s "$workDir/PBSerror.txt";
					#$pbsError=`grep -v 'InfluxDB\\|getcwd' $workDir/PBSerror.txt` if -s "$workDir/PBSerror.txt";
				}
				
				if ($pbsError) {
					open(FILEERROR,">>$fileError");
					print FILEERROR $pbsError;
					close(FILEERROR);
					exit;
				}
				
				if (-s "$workDir/ERROR.txt" ) {
					open(FILEERROR,">>$fileError");
					print FILEERROR `head -1 $workDir/ERROR.txt `;
					close(FILEERROR);
					if(-s $outputFile){
						if (`grep "Found 0 Decoys" $outputFile`){
							my $errorDecoy=`grep "Decoys" $outputFile`;
							my $nbErrorDecoy=`grep -c "Decoys" $outputFile`;
							if ($nbErrorDecoy>=2) {
								my ($errorTandem,$errorDat)=split(/\n/,$errorDecoy);
								if ($errorTandem=~/Found 0 Decoys/) {
									open(FILEERROR,">>$fileError");
									print FILEERROR "No DECOY found in X! Tandem data.";
									close(FILEERROR);
									exit;
								}
								elsif($errorDat=~/Found 0 Decoys/){
									open(FILEERROR,">>$fileError");
									print FILEERROR "No DECOY found in Mascot data (.dat).";
									close(FILEERROR);
									exit;
								}
							}
							else{
								if ($errorDecoy=~/Found 0 Decoys/) {
									open(FILEERROR,">>$fileError");
									print FILEERROR "No DECOY found in selected files.";
									close(FILEERROR);
									exit;
								}
							}
						}
						else{
							open(FILEERROR,">>$fileError");
							print FILEERROR `tail -3 $outputFile `;
							close(FILEERROR);
							exit;
						}
					}
				}
				elsif(-e $outputFile && `tail -1 $outputFile` !~ /SpectraST\sfinished\sat\s\D+\s\D+\s+\d+\s\d{2}:\d{2}:\d{2}\s\d{4}\swithout\serror.\s/){		## Consensus library error
					open(FILEERROR,">>$fileError");
					print FILEERROR "***WARNING***",`tail -1 $outputFile`;
					close(FILEERROR);
					last;
				}
				
			
			}
			else { ###>Run job on Web server
				system "bash $workDir/script.sh";
				if (-s "$workDir/ERROR.txt") {
					open(FILEERROR,">>$fileError");
					print FILEERROR `tail $workDir/ERROR.txt `;
					close(FILEERROR);
					if(-s $outputFile){
						if (`grep "Found 0 Decoys" $outputFile`){
							my $errorDecoy=`grep "Decoys" $outputFile`;
							my $nbErrorDecoy=`grep -c "Decoys" $outputFile`;
							if ($nbErrorDecoy>=2) {
								my ($errorTandem,$errorDat)=split(/\n/,$errorDecoy);
								if ($errorTandem=~/Found 0 Decoys/) {
									open(FILEERROR,">>$fileError");
									print FILEERROR "No DECOY found in X! Tandem data.";
									close(FILEERROR);
									exit;
								}
								elsif($errorDat=~/Found 0 Decoys/){
									open(FILEERROR,">>$fileError");
									print FILEERROR "No DECOY found in Mascot data (.dat).";
									close(FILEERROR);
									exit;
								}
							}
							else{
								if ($errorDecoy=~/Found 0 Decoys/) {
									open(FILEERROR,">>$fileError");
									print FILEERROR "No DECOY found in selected files.";
									close(FILEERROR);
									exit;
								}
							}
						}
					}
					else{
						open(FILEERROR,">>$fileError");
						print FILEERROR `tail -3 $outputFile`;
						close(FILEERROR);
						exit;
					}
				}
				elsif(-e $outputFile && `tail -1 $outputFile` !~ /SpectraST\sfinished\sat\s\D+\s\D+\s+\d+\s\d{2}:\d{2}:\d{2}\s\d{4}\swithout\serror.\s/){		## Consensus library error
					open(FILEERROR,">>$fileError");
					print FILEERROR "***WARNING***",`tail -1 $outputFile`;
					close(FILEERROR);
					last;
				}
			}
			sleep 15;			
		}
		elsif ($action eq "merge"){
			open(FILESTAT,">>$fileStat");
			print FILESTAT "Recovering library's files.\n";
			close(FILESTAT);
			foreach my $libIDFile ($libID1,$libID2){
				opendir (DIR, "$libPath/SwLib_$libIDFile/");
				while (my $file = readdir (DIR)){
					if ($file =~/\.tar\.gz/){
						(my $folderName=$file)=~s/\.tar\.gz//;
						system "cd $libPath/SwLib_$libIDFile; tar -zxf $file;";
						system "rm $libPath/SwLib_$libIDFile/dbinfo_$folderName.txt;" if -e "$libPath/SwLib_$libIDFile/dbinfo_$folderName.txt";
						system "rm $libPath/SwLib_$libIDFile/filelist_$folderName.txt;" if -e "$libPath/SwLib_$libIDFile/filelist_$folderName.txt";
						system "rm $libPath/SwLib_$libIDFile/sortie_$folderName.txt;" if -e "$libPath/SwLib_$libIDFile/sortie_$folderName.txt";
						system "rm $libPath/SwLib_$libIDFile/script_$folderName.sh;" if -e "$libPath/SwLib_$libIDFile/script_$folderName.sh";
					}
				}
				close DIR;
				opendir (CONS,"$libPath/SwLib_$libIDFile/");
				while (my $file2 = readdir (CONS)){
					if ($file2=~/SpecLib_Lib\d+_v\d+/){
						my $taille=`stat -c "%s" $libPath/SwLib_$libIDFile/$file2`;
						push(@taille,int($taille));
						move ("$libPath/SwLib_$libIDFile/$file2",$workDir);
					}
				}
				close CONS;
			}
			@tailleSort=sort @taille;
			
			open(FILESTAT,">>$fileStat");
			print FILESTAT "Creation of spectra librairies.\n";
			close(FILESTAT);
			my $command;
			if ($split == 1) {
				$command="grep -hUv '###' $libPath/SwLib_$libID1/$libName1.sptxt $libPath/SwLib_$libID2/$libName2.sptxt >>$workDir/SpecLib_cons_concat.sptxt";
			}
			elsif ($split == 0) {
				$command="$tppPath/spectrast -cN$workDir/SpecLib_cons -cI$fragmentation -cAC $workDir/*.splib >>$outputFile 2>&1";
			}
			
			my $clusterCommandString;
			if ($clusterInfo{'on'}){
				$clusterCommandString=$clusterInfo{'buildCommand'}->($workDir,$command);
			}
			else{
				$clusterCommandString=$command;
			}
			my $maxHours=20;
			my $bashFile="$workDir/script.sh";
			$maxMem=($tailleSort[-1]/1073741824)*10;
			$maxMem*=2 if $clusterInfo{'name'} eq 'CentOS';
			$maxMem=($maxMem < 20) ? 20 : ($maxMem > 100) ? 100: sprintf("%.0f",$maxMem) ;
			$maxMem.='Gb';
			open(BASH,"+>",$bashFile);
			print BASH qq
	|#!/bin/bash
	
	##resources
	#PBS -l mem=$maxMem
	#PBS -l nodes=1:ppn=1
	#PBS -l walltime=$maxHours:00:00
	#PBS -q batch
	
	##Information
	#PBS -N Swath_Lib_$time
	#PBS -M marine.le-picard\@curie.fr
	#PBS -m abe
	#PBS -o $workDir/PBS.txt
	#PBS -e $workDir/PBSerror.txt
	
	## Command
	$clusterCommandString
	echo "<BR><BR>End_Swath_Lib_$time" >>$workDir/END.txt
	|;
			close BASH;
			
			
			###> Execute bash 
			if ($clusterInfo{'on'}) {
				my ($jobClusterID) = $clusterInfo{'sendToCluster'}->($bashFile);
				
				# Add job to DB
				$dbh->do("UPDATE JOB_HISTORY SET ID_JOB_CLUSTER = CONCAT(ID_JOB_CLUSTER, ';C$jobClusterID') WHERE ID_JOB='$time'");
				$dbh->commit;
				
				###>Waiting for job to run
				my $pbsError;
				my $nbWhile=0;
				my $maxNbWhile=$maxHours*60*2;
				while (!$pbsError  && !-s "$workDir/END.txt") {
					sleep 30;
					$nbWhile++;
					if ($nbWhile > $maxNbWhile) {
						open(FILEERROR,">>$fileError");
						print FILEERROR "Aborting: File merge is taking too long.";
						close(FILEERROR);
						exit;
					}
					$pbsError=$clusterInfo{'checkError'}->("$workDir/PBSerror.txt") if -s "$workDir/PBSerror.txt";
					#$pbsError=`grep -v getcwd $workDir/PBSerror.txt` if -s "$workDir/PBSerror.txt";
				}
				if ($pbsError) {
					open(FILEERROR,">>$fileError");
					print FILEERROR $pbsError;
					close(FILEERROR);
					exit;
				}
			}
			else { ###>Run job on Web server
				system "bash $bashFile";
			}
			sleep 5;
		}
		
		## check errors in SpectraST
		if (-s "$workDir/spectrast.log") {
			my $numSpectraError=`grep -c ERROR $workDir/spectrast.log`;
			my $numNterModError=`grep -c 'Peptide ID has unknown modification: ' $workDir/spectrast.log`;		 ### skip N-term modifications : in sequest file N-term modification are tag with +42, in mascot or Xtandem : +43 ; spectrast do not recognize +42 modification
			chomp $numNterModError;
			if ($numSpectraError!=0) {
				if($numNterModError!=$numSpectraError){
					open(FILEERROR,">>$fileError");
					print FILEERROR "You had $numSpectraError errors during SpectraST process.";
					close(FILEERROR);
					exit;
				}
				else{
					open(FILESTAT,">>$fileStat");
					print FILESTAT "You had $numNterModError \"Peptide has unknown modification\" errors";
					close(FILESTAT);
				}
			}
		}
		
		
		if (-e "$workDir/SpecLib_cons_concat.sptxt" || -e "$workDir/SpecLib_cons.sptxt"){
			my ($numPep,$finalFileLib);
			####> Recovering the number of spectra in the library
			if (-e "$workDir/SpecLib_cons_concat.sptxt") {
				$numPep=`grep -c LibID  $workDir/SpecLib_cons_concat.sptxt`;
				open (INFILE,"<","$workDir/SpecLib_cons_concat.sptxt") or die("open: $!");
				$finalFileLib="$workDir/SpecLib_cons_concat.sptxt";
			}
			elsif(-e "$workDir/SpecLib_cons.sptxt"){
				$numPep=`grep -c LibID  $workDir/SpecLib_cons.sptxt`;
				open (INFILE,"<","$workDir/SpecLib_cons.sptxt") or die("open: $!");
				$finalFileLib="$workDir/SpecLib_cons.sptxt";
			}
			chomp $numPep;
			
	
			#########################################################
			####>Recover peptides modifications and protein list<####
			#########################################################
			open(FILESTAT,">>$fileStat");
			print FILESTAT "Recovering peptides modifications.\n";
			close(FILESTAT);
			my $line;
			my (%swathModifications,%protLibList,%protSpecificList,%numPepMod);
			while (($line=<INFILE>)){
				if ($line=~/Mods=(\d)\/(\S*)/){
					my($varModCode,$residues,$positions);
					my $modList=$2;
					my @result=&promsMod::convertVarModStringSwath($modList);
					foreach my $res (@result){
						my @resTab=split(/!/,$res);
						$varModCode=$resTab[0];
						$residues=$resTab[1];
						$positions=$resTab[2];
						$residues='=' if $positions eq '=';
						if ($varModCode eq 'GlyGly' || $varModCode eq 'Phospho' || $varModCode eq 'Acetyl'){
							$numPepMod{$varModCode}{'Multiple'}++;
						}
					}
					$numPepMod{'GlyGly'}{'Single'}++ if $modList=~/GlyGly/;
					$numPepMod{'Phospho'}{'Single'}++ if $modList=~/Phospho/;
					$numPepMod{'Acetyl'}{'Single'}++ if $modList=~/Acetyl/;
					$swathModifications{$varModCode}{$residues}=1;
				}
				if($line=~/Protein=\d\/(\S+)/){
					if ($split==0) {
						my @protIDs=split(/\//,$1);
						for (my $i=0; $i<@protIDs; $i++){
							next if $protIDs[$i]=~/reverse/;
							my $quoteProt=quotemeta($protIDs[$i]);
							$protLibList{$protIDs[$i]}=1;
							$protSpecificList{$protIDs[$i]}=1 unless scalar @protIDs>2;
						}
					}
					elsif ($split==1){
						my @protIDs=split(/\//,substr($1,10));
						for (my $i=0; $i<@protIDs; $i++){
							next if $protIDs[$i]=~/reverse/;
							my $quoteProt=quotemeta($protIDs[$i]);
							$protLibList{$protIDs[$i]}=1;
							$protSpecificList{$protIDs[$i]}=1 unless scalar @protIDs>2;
						}
					}
				}
			}
			close INFILE;
			
			###> Convertion of modifications that are not on a mass format (K[GlyGly] -> K[242])	#error during spectrast2tsv.py if K[GlyGly] !!!
			my %TPPModifCode=&promsConfig::getTPPModificationCode;
			foreach my $code (keys %TPPModifCode){
				if ($swathModifications{$code}){
					system "sed -i 's/\\\[$code\\\]/\\\[$TPPModifCode{$code}\\\]/g' $finalFileLib";			
				}
			}
			
			my $Stat;
			my $numProt=scalar keys %protLibList;
			my $numProtSpe=scalar keys %protSpecificList;
			my $numPepMod;
			if(%numPepMod){
				foreach my $mod (keys %numPepMod){
					$numPepMod.='&' if $numPepMod;
					$numPepMod.="$mod".':'.$numPepMod{$mod}{'Single'}.'/'.$numPepMod{$mod}{'Multiple'};
				}
			}
			my $entry->{'NUM_PEP'}=$numPep;
			$entry->{'NUM_PROT'}=$numProt;
			$entry->{'NUM_PROT_SPECIFIC'}=$numProtSpe;
			$entry->{'NUM_PEP_MOD'}=$numPepMod if $numPepMod;
	
			my $xmlParser = XML::Simple->new( NoAttr=>1, RootName=>'NUM_ENTRY');
			$Stat=$xmlParser->XMLout($entry);
	
			#################################
			####>Insertion into Database<####
			#################################
			$dbh=&promsConfig::dbConnect;
			open(FILESTAT,">>$fileStat");
			print FILESTAT "Storing data into database.\n";
			close(FILESTAT);
			my $versionFile=1;
			my ($oldEngine,$params,$nameFinalFiles,$oldFileNumberTot);
			my $version=1;
			## save the previous files version
			if ($action eq 'update') {
				my $finalDir="$libPath/SwLib_$libIDFile";
				
				my ($nameOldLib,$desS,$paramStrgS,$oldVersion,$updateDateS,$updateUserS,$statS,$dbS)=$dbh->selectrow_array("SELECT NAME, DES, PARAM_STRG, VERSION_NAME, UPDATE_DATE, UPDATE_USER, STATISTICS, ID_DATABANK FROM SWATH_LIB SL, DATABANK_SWATHLIB DS WHERE SL.ID_SWATH_LIB=$libIDFile AND DS.ID_SWATH_LIB=SL.ID_SWATH_LIB");
				my $xmlS =XML::Simple-> new (KeepRoot=>1);
				my $xmlOldStat = $xmlS->XMLin($statS);
				my $pepS=$xmlOldStat->{'NUM_ENTRY'}->{'NUM_PEP'};
				my $protS=$xmlOldStat->{'NUM_ENTRY'}->{'NUM_PROT'};
				my $protSpeS=$xmlOldStat->{'NUM_ENTRY'}->{'NUM_PROT_SPECIFIC'};
				$updateDateS=($updateDateS)? $updateDateS : '';
				$updateUserS=($updateUserS)? $updateUserS : '';
				open(SAVEFILE,">","$finalDir/dbinfo_Lib$libIDFile\_$oldVersion.txt");
				print SAVEFILE "des=$desS\tversion=$oldVersion\tdate=$updateDateS\tuser=$updateUserS\tdatabank=$dbS\tpep=$pepS\tprot=$protS\tprotspe=$protSpeS\n";
				print SAVEFILE "$paramStrgS";
				close SAVEFILE;
				
				$oldVersion=~s/v//;
				$version=$oldVersion+1;
				
				my $xmlOldParam = $xmlS->XMLin($paramStrgS);
				open(FILELIST,">>","$finalDir/filelist_Lib$libIDFile\_v$oldVersion.txt");
				foreach my $paramFiles (@{$xmlOldParam->{'SWATH_LIB_DATA'}->{'FILES'}->{'FILE'}}) {
					my $oldFile=$paramFiles->{'FileName'};
					my $olVersionFile=$paramFiles->{'version'};
					print FILELIST "$oldFile,$olVersionFile&";
					$olVersionFile=~s/v//;
					if ($olVersionFile>$versionFile) {
						$versionFile=$olVersionFile;
					}
					$oldEngine=$xmlOldParam->{'SWATH_LIB_DATA'}->{'SEARCH_ENGINES'}->{'NAME'};
				}
				close FILELIST;
				
				$oldFileNumberTot=$xmlOldParam->{'SWATH_LIB_DATA'}->{'FILE_NUMBER'}->{'NUMBER'};
				
				
				my $archiveGZ=$finalDir.'/Lib'.$libIDFile.'_v'.$oldVersion.'.tar.gz';
				my $archive=$finalDir.'/Lib'.$libIDFile.'_v'.$oldVersion.'.tar';
				move("$finalDir/sortie.txt","$finalDir/sortie_Lib$libIDFile\_v$oldVersion.txt");
				move("$finalDir/script.sh","$finalDir/script_Lib$libIDFile\_v$oldVersion.sh");
				system "cd $finalDir; gunzip $archiveGZ; tar -rf $archive filelist_Lib$libIDFile\_v$oldVersion.txt;tar -rf $archive dbinfo_Lib$libIDFile\_v$oldVersion.txt; tar -rf $archive sortie_Lib$libIDFile\_v$oldVersion.txt;tar -rf $archive script_Lib$libIDFile\_v$oldVersion.sh; gzip $archive; rm dbinfo_Lib$libIDFile\_v$oldVersion.txt; rm filelist_Lib$libIDFile\_v$oldVersion.txt; rm sortie_Lib$libIDFile\_v$oldVersion.txt; rm script_Lib$libIDFile\_v$oldVersion.sh";
				
				
				if (defined $xmlOldParam->{'SWATH_LIB_DATA'}->{'PARENT_LIB'}){			## for merge libraries && mergelib
					if (defined $xmlOldParam->{'SWATH_LIB_DATA'}->{'PARENT_LIB'}->{'PARENT'}){
						push @{$params->{'PARENT_LIB'}->{'PARENT'}},{'LibName'=>"$xmlOldParam->{'SWATH_LIB_DATA'}->{'PARENT_LIB'}->{'PARENT'}->{'LibName'}",'ID'=>"$xmlOldParam->{'SWATH_LIB_DATA'}->{'PARENT_LIB'}->{'PARENT'}->{'ID'}"};
					}
					else{
						push @{$params->{'PARENT_LIB'}->{'PARENT1'}},{'LibName'=>"$xmlOldParam->{'SWATH_LIB_DATA'}->{'PARENT_LIB'}->{'PARENT1'}->{'LibName'}",'ID'=>"$xmlOldParam->{'SWATH_LIB_DATA'}->{'PARENT_LIB'}->{'PARENT1'}->{'ID'}"};
						push @{$params->{'PARENT_LIB'}->{'PARENT2'}},{'LibName'=>"$xmlOldParam->{'SWATH_LIB_DATA'}->{'PARENT_LIB'}->{'PARENT2'}->{'LibName'}",'ID'=>"$xmlOldParam->{'SWATH_LIB_DATA'}->{'PARENT_LIB'}->{'PARENT2'}->{'ID'}"};
					}
				}
				
				$versionFile+=1;
				###> Adding library name used to create the new library (just for an update or a merge)
				push @{$params-> {'SOURCE_LIB'}},{'LibName'=>"$nameOldLib",'LibId'=>$libIDFile};
				$nameFinalFiles=$nameOldLib;                ##final files names= old library name
			}	
			
			###> Creating xml file that will be integrated to the database (SWATH_LIB : PARAM_STRG)
			my $xmlParams;
			if ($action eq "merge"){
				my $xmlS =XML::Simple-> new (KeepRoot=>1);
				my $paramStrg1=$dbh->selectrow_array("SELECT PARAM_STRG FROM SWATH_LIB WHERE ID_SWATH_LIB=$libID1");
				my $xmlOldParam1=$xmlS->XMLin($paramStrg1);
				my $nbFiles1=$xmlOldParam1->{'SWATH_LIB_DATA'}->{'FILE_NUMBER'}->{'NUMBER'};
				my $paramStrg2=$dbh->selectrow_array("SELECT PARAM_STRG FROM SWATH_LIB WHERE ID_SWATH_LIB=$libID2");
				my $xmlOldParam2=$xmlS->XMLin($paramStrg2);
				my $nbFiles2=$xmlOldParam2->{'SWATH_LIB_DATA'}->{'FILE_NUMBER'}->{'NUMBER'};
				my $fileNumberTot=$nbFiles2+$nbFiles1;
				
				push @{$params->{'PARENT_LIB'}->{'PARENT1'}},{'LibName'=>"$libName1",'ID'=>"$libID1"};
				push @{$params->{'PARENT_LIB'}->{'PARENT2'}},{'LibName'=>"$libName2",'ID'=>"$libID2"};
				if ($split==0) {
					push @{$params->{'COMMANDS_LINES'}->{'COMMAND'}},{'CommandName'=>'merge'};
					push @{$params->{'COMMANDS_LINES'}->{'COMMAND'}},{'CommandName'=>'spectrast','CommandRank'=>1,'PARAMS'=>[{'opvalue'=>'outputfile','selectby'=>'auto','option'=>'-cN'},{'opvalue'=>"$fragmentation",'selectby'=>'auto','option'=>'-cI'},{'opvalue'=>'no value','option'=>'-cAC'}]};				
				}
				$params-> {'FILE_NUMBER'}->{'NUMBER'}=$fileNumberTot;
				my $xmlParser = XML::Simple->new( NoAttr=>1, RootName=>'SWATH_LIB_DATA');
				$xmlParams=$xmlParser->XMLout($params);
			}
			else{
				###> Recovering files names
				foreach my $file (@inputFiles) {
					push @{$params-> {'FILES'} -> {'FILE'}},{'FileName'=>"$file",'version'=>"v$versionFile"};
				}
	
				open(FDR,"<","$workDir/FDRscore.txt");
				my $fdrMin;
				while (my $line=<FDR>) {
					$fdrMin=$line;
				}
				chomp $fdrMin;
				
				###> Recovering of all parameters used in bash script
				push @{$params-> {'COMMANDS_LINES'}->{'COMMAND'}},{'CommandName'=>'Mascot2XML','CommandRank'=>1,'PARAMS'=>[{'opvalue' => 'database','selectby' => 'user', 'option'=>'-D'},{'opvalue' => 'trypsin','selectby' => 'auto','option'=>'-E'},{'opvalue' => 'otgz','selectby' => 'auto','option'=>'-n'}]};
				push @{$params-> {'COMMANDS_LINES'}->{'COMMAND'}},{'CommandName'=>'xinteract','CommandRank'=>2,'PARAMS'=>[{'opvalue' => 'ARPd','selectby' => 'auto', 'option'=>'-O'},{'opvalue' => "$decoyTag",'selectby' => 'auto','option'=>'-d'},{'opvalue' => 'outputfile','selectby' => 'auto','option'=>'-N'}]};
				push @{$params-> {'COMMANDS_LINES'}->{'COMMAND'}},{'CommandName'=>'InterProphetParser','CommandRank'=>3,'PARAMS'=>{'opvalue' => "$decoyTag",'selectby' => 'auto', 'option'=>'DECOY'}};
				push @{$params-> {'COMMANDS_LINES'}->{'COMMAND'}},{'CommandName'=>'Mayu.pl','CommandRank'=>4,'PARAMS'=>[{'opvalue' => 'inputfile','selectby' => 'auto', 'option'=>'-A'},{'opvalue' => 'database','selectby' => 'user','option'=>'-C'},{'opvalue' => "$decoyTag",'selectby' => 'auto','option'=>'-E'},{'opvalue' => '0.01','selectby' => 'auto','option'=>'-G'},{'pvalue' => '51','selectby' => 'auto','option'=>'-H'},{'opvalue' => "$missedCleavage",'selectby' => 'user','option'=>'-I'},{'opvalue' => "$fdr:t",'selectby' => 'user','option'=>'-P'},{'opvalue' => "$fdrType",'selectby' => 'user','option'=>'-fdr'}]};
				push @{$params-> {'COMMANDS_LINES'}->{'COMMAND'}},{'CommandName'=>'spectrast','CommandRank'=>5,'PARAMS'=>[{'opvalue' => 'outputfile','selectby' => 'auto', 'option'=>'-cN'},{'opvalue' => "$fragmentation",'selectby' => 'auto','option'=>'-cI'},{'opvalue' => "Protein! ~ $decoyTag",'selectby' => 'auto','option'=>'-cf'},{'opvalue' => "$fdrMin",'selectby' => 'auto','option'=>'-cP'},{'opvalue' => 'irtfile','selectby' => 'user','option'=>'-c_IRT'},{'opvalue' => 'inputfile','selectby' => 'auto','option'=>'-c_IRR'}]};
	
				if ($libOption eq 'unsplit'){
					push @{$params-> {'COMMANDS_LINES'}->{'COMMAND'}},{'CommandName'=>'spectrast','CommandRank'=>6,'PARAMS'=>[{'opvalue' => 'outputfile','selectby' => 'auto', 'option'=>'-cN'},{'opvalue' => "$fragmentation",'selectby' => 'auto','option'=>'-cI'},{'opvalue' => 'no value','option'=>'-cAC'}]};
				}
				elsif ($libOption eq 'split'){
					push @{$params-> {'COMMANDS_LINES'}->{'COMMAND'}},{'CommandName'=>'spectrast_cluster.py','CommandRank'=>6,'PARAMS'=>{'opvalue' => '2','selectby' => 'auto', 'option'=>'-d'}};
					push @{$params-> {'COMMANDS_LINES'}->{'COMMAND'}},{'CommandName'=>'spectrast','CommandRank'=>7,'PARAMS'=>[{'opvalue' => 'outputfile','selectby' => 'auto', 'option'=>'-cN'},{'opvalue' => "$fragmentation",'selectby' => 'auto','option'=>'-cI'}]};
					push @{$params-> {'COMMANDS_LINES'}->{'COMMAND'}},{'CommandName'=>'spectrast','CommandRank'=>8,'PARAMS'=>[{'opvalue' => 'outputfile','selectby' => 'auto', 'option'=>'-cN'},{'opvalue' => "$fragmentation",'selectby' => 'auto','option'=>'-cI'},{'opvalue' => 'no value','option'=>'-cAC'}]};
				}
				if (("@inputFiles"=~ m/.dat/) && ("@inputFiles"=~ m/.tandem/)) {$params-> {'SEARCH_ENGINES'}->{'NAME'}='MASCOT,XTANDEM';}
				elsif ("@inputFiles"=~ m/.dat/) {
					if ($oldEngine) {
						$params-> {'SEARCH_ENGINES'}->{'NAME'}=($oldEngine eq 'MASCOT,XTANDEM' || $oldEngine eq 'XTANDEM') ? 'MASCOT,XTANDEM' : 'MASCOT';
					}else{$params-> {'SEARCH_ENGINES'}->{'NAME'}='MASCOT';}
				}
				elsif ("@inputFiles"=~ m/.tandem/) {
					if ($oldEngine) {
						$params-> {'SEARCH_ENGINES'}->{'NAME'}=($oldEngine eq 'MASCOT,XTANDEM' || $oldEngine eq 'MASCOT') ? 'MASCOT,XTANDEM' : 'XTANDEM';
					}else{$params-> {'SEARCH_ENGINES'}->{'NAME'}='XTANDEM';}
				}
				my $fileNumberTot;
				if ($libCons eq 'mergelib') {
					push @{$params->{'PARENT_LIB'}->{'PARENT'}},{'LibName'=>"$libFileName",'ID'=>"$libIDFile"};
					my $xmlS =XML::Simple-> new (KeepRoot=>1);
					my $paramStrg1=$dbh->selectrow_array("SELECT PARAM_STRG FROM SWATH_LIB WHERE ID_SWATH_LIB=$libIDFile");
					my $xmlOldParam1=$xmlS->XMLin($paramStrg1);
					my $nbFiles1=$xmlOldParam1->{'SWATH_LIB_DATA'}->{'FILE_NUMBER'}->{'NUMBER'};
					$fileNumberTot=$nbFiles1+scalar @inputFiles;
				}
				$fileNumberTot=($fileNumberTot) ? $fileNumberTot : ($oldFileNumberTot)? $oldFileNumberTot + scalar @inputFiles : scalar @inputFiles;
				$params-> {'FILE_NUMBER'}->{'NUMBER'}=$fileNumberTot;
				
				
				###>Transforming $params in xml
				my $xmlParser = XML::Simple->new( NoAttr=>1, RootName=>'SWATH_LIB_DATA');
				$xmlParams=$xmlParser->XMLout($params);
			}
			
			###> Insertion into table SWATH_LIB
			if ($action eq 'update') {
				my $sthLibUpdate2=$dbh->prepare("UPDATE SWATH_LIB SET DES=?, PARAM_STRG=?,USE_STATUS=?, VERSION_NAME=?,UPDATE_DATE=NOW(),UPDATE_USER=?,STATISTICS=?  WHERE ID_SWATH_LIB=?") or die "Couldn't prepare statement: " . $dbh->errstr;
				$sthLibUpdate2->execute($des,$xmlParams,'yes',"v$version",$userID,$Stat,$libIDFile);
				$sthLibUpdate2->finish;
				my $dbRank=0;
				my $match=0;
				my $sthOldDB=$dbh->prepare("SELECT ID_DATABANK,DB_RANK FROM DATABANK_SWATHLIB WHERE ID_SWATH_LIB=?");
				$sthOldDB->execute($libIDFile);
				while (my ($oldDBID,$oldDBRank)=$sthOldDB->fetchrow_array) {
					$dbRank=$oldDBRank if $dbRank<$oldDBRank;
					if ($oldDBID=!$dbID[0]) {
						$match=1;
						last;
					}
				}
				$sthOldDB->finish;
				$dbRank+=1;
				if ($match) {
					my $sthLibDBUpdate=$dbh->prepare("INSERT INTO DATABANK_SWATHLIB (ID_DATABANK,ID_SWATH_LIB,DB_RANK) VALUES (?,?,?)") or die "Couldn't prepare statement: " . $dbh->errstr;
					$sthLibDBUpdate->execute($dbID[0],$libIDFile,$dbRank);
					$sthLibDBUpdate->finish;
				}
			}
			else{
				my $sthLibUpdate2=$dbh->prepare("UPDATE SWATH_LIB SET PARAM_STRG=?, USE_STATUS=?, START_DATE=NOW(), STATISTICS=?, SPLIT=? WHERE ID_SWATH_LIB=?");
				$sthLibUpdate2->execute($xmlParams,"yes",$Stat,$split,$libID);
				$sthLibUpdate2->finish;
				$nameFinalFiles=$libName;    ##final files names= new library name
				##Insertion into DATABANK_SWATHLIB
				if($libCons eq 'new'){
					my $sthDatabankLib=$dbh->prepare("INSERT INTO DATABANK_SWATHLIB (ID_SWATH_LIB,ID_DATABANK,DB_RANK) values (?,?,?)");
					$sthDatabankLib->execute($libID,$dbID[0],1);
					$sthDatabankLib->finish;
				}
				elsif( $action eq 'merge'){
					my $dbRank=0;
					foreach my $dbID (@dbID){
						$dbRank++;
						my $sthDatabankLib=$dbh->prepare("INSERT INTO DATABANK_SWATHLIB (ID_SWATH_LIB,ID_DATABANK,DB_RANK) values (?,?,?)");
						$sthDatabankLib->execute($libID,$dbID,$dbRank);
						$sthDatabankLib->finish;
					}
					my $sthParentLib1=$dbh->prepare("INSERT INTO PARENT_SWATH_LIB (ID_SWATH_LIB,ID_PARENT_SWATH_LIB,VERSION_NAME) values (?,?,?)");
					$sthParentLib1->execute($libID,$libID1,$versionLib1);
					$sthParentLib1->finish;
					my $sthParentLib2=$dbh->prepare("INSERT INTO PARENT_SWATH_LIB (ID_SWATH_LIB,ID_PARENT_SWATH_LIB,VERSION_NAME) values (?,?,?)");
					$sthParentLib2->execute($libID,$libID2,$versionLib2);
					$sthParentLib2->finish;
				}
				elsif($libCons eq 'mergelib'){
					my $dbRank=0;
					foreach my $dbID (@dbID){
						next if $dbID eq '';
						$dbRank++;
						my $sthDatabankLib=$dbh->prepare("INSERT INTO DATABANK_SWATHLIB (ID_SWATH_LIB,ID_DATABANK,DB_RANK) values (?,?,?)");
						$sthDatabankLib->execute($libID,$dbID,$dbRank);
						$sthDatabankLib->finish;
					}
					my $sthParentLib=$dbh->prepare("INSERT INTO PARENT_SWATH_LIB (ID_SWATH_LIB,ID_PARENT_SWATH_LIB,VERSION_NAME) values (?,?,?)");
					$sthParentLib->execute($libID,$libIDFile,$versionLib1);
					$sthParentLib->finish;
				}
			}
			$dbh->commit;
	
			###>Insertion into SWATH_LIB_MODIFICATION
			foreach my $keys (keys(%swathModifications)){
				my $residuesList=join('',keys(%{$swathModifications{$keys}}));
				my $modID=&promsMod::getModificationIDfromString($dbh,$keys,$residuesList);
				my $sthLibMod;
				if ($action eq 'update') {
					my $modTest=$dbh->selectrow_array("SELECT COUNT(*) FROM SWATH_LIB_MODIFICATION WHERE ID_SWATH_LIB=$libIDFile AND ID_MODIFICATION=$modID");
					if ($modTest==0) {
						$sthLibMod=$dbh->prepare("INSERT INTO SWATH_LIB_MODIFICATION (ID_SWATH_LIB,ID_MODIFICATION,SPECIFICITY) VALUES (?,?,?)");
						$sthLibMod->execute($libIDFile,$modID,$residuesList);
						$sthLibMod->finish;
					}
					else{
						my $oldResiduesList=$dbh->selectrow_array("SELECT SPECIFICITY FROM SWATH_LIB_MODIFICATION WHERE ID_SWATH_LIB=$libIDFile AND ID_MODIFICATION=$modID");
						if ($residuesList ne $oldResiduesList){
							$sthLibMod=$dbh->prepare("UPDATE SWATH_LIB_MODIFICATION SET SPECIFICITY=? WHERE ID_SWATH_LIB=? AND ID_MODIFICATION=?");
							$sthLibMod->execute($residuesList,$libIDFile,$modID);
							$sthLibMod->finish;
						}
					}
				}
				elsif ($libCons eq 'new' || $libCons eq 'mergelib' || $action eq 'merge'){
					$sthLibMod=$dbh->prepare("INSERT INTO SWATH_LIB_MODIFICATION (ID_SWATH_LIB,ID_MODIFICATION,SPECIFICITY) values (?,?,?)");
					$sthLibMod->execute($libID,$modID,$residuesList);
					$sthLibMod->finish;
				}
			}
			$dbh->commit;
			
			############################
			####>Moving final files<####
			############################
			open(FILESTAT,">>$fileStat");
			print FILESTAT "Moving final files.\n";
			close(FILESTAT);
			
			my $finalDir="$libPath/SwLib_$libID";
			mkdir $finalDir unless -e $finalDir;
			my $versionDir='Lib'.$libID.'_v'.$version;
			
			if ($action eq 'merge'){
				system "cd $workDir; tar -czf $versionDir.tar.gz SpecLib_Lib*";
			}
			else{
				###> Moving final files to Swath_Lib and removing temp folder
				rename ("$workDir/SpecLib.sptxt","$workDir/SpecLib_Lib$libID\_v$version.sptxt");
				rename ("$workDir/SpecLib.splib","$workDir/SpecLib_Lib$libID\_v$version.splib");
				rename ("$workDir/SpecLib.pepidx","$workDir/SpecLib_Lib$libID\_v$version.pepidx");
				rename ("$workDir/SpecLib.spidx","$workDir/SpecLib_Lib$libID\_v$version.spidx");
				rename ("$workDir/SpecLib.mrm","$workDir/SpecLib_Lib$libID\_v$version.mrm") if(-s "$workDir/SpecLib.mrm");
				if ($libCons eq 'mergelib'){
					system "cd $workDir; tar -czf $versionDir.tar.gz SpecLib_Lib*";
				}
				else{
					system "cd $workDir; tar -cf $versionDir.tar SpecLib_Lib$libID\_v$version.sptxt";
					system "cd $workDir; tar -rf $versionDir.tar SpecLib_Lib$libID\_v$version.splib";
					system "cd $workDir; tar -rf $versionDir.tar SpecLib_Lib$libID\_v$version.pepidx";
					system "cd $workDir; tar -rf $versionDir.tar SpecLib_Lib$libID\_v$version.spidx";
					system "cd $workDir; tar -rf $versionDir.tar SpecLib_Lib$libID\_v$version.mrm" if(-s "SpecLib_Lib$libID\_v$version.mrm");
					system "cd $workDir; gzip $versionDir.tar";
				}
			}
			system "mv $workDir/$versionDir.tar.gz $finalDir";
			
			open(FILESTAT,">>$fileStat");
			print FILESTAT "Moving consensus library.\n";
			close(FILESTAT);
			
			#move ("$workDir/$versionDir.tar.gz","$finalDir/$versionDir.tar.gz");
			
			###> moving final files (consensus library)
			if ($libOption eq 'split' || ($action eq 'merge' && $split==1)){
				system "mv $workDir/SpecLib_cons_concat.sptxt $finalDir/$nameFinalFiles.sptxt";
			}
			elsif($libOption eq 'unsplit' || ($action eq 'merge' && $split==0)){
				system "mv $workDir/SpecLib_cons.sptxt $finalDir/$nameFinalFiles.sptxt";
				system "mv $workDir/SpecLib_cons.splib $finalDir/$nameFinalFiles.splib";
				system "mv $workDir/SpecLib_cons.pepidx $finalDir/$nameFinalFiles.pepidx";
				system "mv $workDir/SpecLib_cons.spidx $finalDir/$nameFinalFiles.spidx";
				system "mv $workDir/SpecLib_cons.mrm $finalDir/$nameFinalFiles.mrm" if(-s "$workDir/SpecLib_cons.mrm");
			}
			system "cp $workDir/sortie.txt $finalDir";
			system "cp $workDir/script.sh $finalDir" if (-e "$workDir/script.sh");
	
			open(FILESTAT,">>$fileStat");
			print FILESTAT "Ended";
			close(FILESTAT);
		}
		else{
			open(FILEERROR,">>$fileError");
			print FILEERROR "Spectrast did not work.";
			close(FILEERROR);
			exit;
		}
	}
}
print "</CENTER></BODY></HTML>";


sub ajaxSelectSecondMergeLib{
    my $libID1=param('libID1');
    my ($rtID1,$split1,$dbType)=$dbh->selectrow_array("SELECT SL.ID_REFERENCE_RT,SL.SPLIT,D.ID_DBTYPE FROM SWATH_LIB SL, DATABANK D, DATABANK_SWATHLIB DS WHERE SL.ID_SWATH_LIB=$libID1 AND DS.ID_SWATH_LIB=SL.ID_SWATH_LIB AND DS.ID_DATABANK=D.ID_DATABANK AND DS.DB_RANK=1");
    my $sthLibName=$dbh->prepare("SELECT S.NAME,S.ID_SWATH_LIB, DT.NAME FROM SWATH_LIB S, DATABANK_SWATHLIB DS, DATABANK D, DATABANK_TYPE DT WHERE S.ID_SWATH_LIB=DS.ID_SWATH_LIB AND DS.ID_DATABANK=D.ID_DATABANK AND D.ID_DBTYPE=$dbType AND DT.ID_DBTYPE=D.ID_DBTYPE AND S.SPLIT=$split1 AND S.ID_REFERENCE_RT=$rtID1 AND S.ID_SWATH_LIB!=$libID1 AND S.USE_STATUS='yes'") or die "Couldn't prepare statement: " . $dbh->errstr;
    $sthLibName->execute;
    if ($sthLibName->rows==0) {print "No library available.";}
    else{
        my %lib2Hash;
        while (my($libName2,$libID2,$dbType)=$sthLibName->fetchrow_array){
            $lib2Hash{$libName2}="$libID2&$dbType";
        }
        
        print qq
        |<SELECT name="libfile2" id="libfile2" required>
             <option value="">-= Select Library =-</option>
        |;
        foreach my $libName2 (sort {lc $a cmp lc $b} keys %lib2Hash){
			my ($libID,$dbType)=split(/&/,$lib2Hash{$libName2});
            print "<option value=\"$libID\">$libName2 [$dbType]</option>";
        }
        print "</SELECT>";
    }
	$sthLibName->finish;
}

sub selectLibraryName{
    my $sthLibName=$dbh->prepare("SELECT SL.NAME,SL.ID_SWATH_LIB,SPLIT,SL.VERSION_NAME,DT.NAME,D.ID_DATABANK FROM SWATH_LIB SL, DATABANK_SWATHLIB DS, DATABANK_TYPE DT, DATABANK D WHERE SL.USE_STATUS='yes' AND DS.ID_SWATH_LIB=SL.ID_SWATH_LIB AND DS.ID_DATABANK=D.ID_DATABANK AND D.ID_DBTYPE=DT.ID_DBTYPE AND DS.DB_RANK=1") or die "Couldn't prepare statement: " . $dbh->errstr;
    $sthLibName->execute;
    if ($sthLibName->rows==0) {print "No library available.";}
    else{
        my %libHash;
        while (my($libName,$libID,$splitLib,$version,$dbType,$dbID)=$sthLibName->fetchrow_array){
			#print "$libName,$libID,$splitLib,$version,$dbType,$dbID<BR>";
            $libHash{$libName}="$libID&$splitLib&$version&$dbType&$dbID";
        }
        
        print qq
        |<BR>&nbsp;<SELECT name="libfile" id="libfile" onChange="selectLibraryName('$option',this.value)" required>
             <option value="">-= Select Library =-</option>
        |;
        foreach my $libName (sort {lc $a cmp lc $b} keys %libHash){
			my ($libID,$splitLib,$version,$dbType,$dbID)=split(/&/,$libHash{$libName});
            print "<option value=\"$libID&$splitLib&$version&$dbID&$libName\">$libName [$dbType]</option>";
        }
        print "</SELECT>";
    }
	$sthLibName->finish;
}

sub ajaxSelectExperiment {
    my $sthExperimentList=$dbh->prepare("SELECT NAME, ID_EXPERIMENT FROM EXPERIMENT WHERE ID_PROJECT=$projectID ") or die "Couldn't prepare statement: " . $dbh->errstr;
    $sthExperimentList->execute;
    if ($sthExperimentList->rows==0) {print "No experiment in that project. Choose an other project.";}
    else{
        my %experimentHash;
        while (my($experimentName,$experimentID)=$sthExperimentList->fetchrow_array){
            $experimentHash{$experimentName}=$experimentID;
        }
        
        print qq
        |<SELECT name="experimentFiles" onChange="ajaxSelectSample(this.value)">
            <option value="">=-Select Experiment-=</option>
        |;
        foreach my $experimentName (sort {lc $a cmp lc $b} keys %experimentHash){
            print "<option value=\"$experimentHash{$experimentName}\">$experimentName</option>";
        }
        print "</SELECT>";
    }
    $sthExperimentList->finish;
}

sub ajaxSelectSample{
    my $sthSampleList=$dbh->prepare("SELECT S.NAME, S.ID_SAMPLE, A.DATA_FILE, A.ID_ANALYSIS ,D.DECOY_TAG FROM SAMPLE S LEFT JOIN ANALYSIS A ON S.ID_SAMPLE=A.ID_SAMPLE LEFT JOIN ANALYSIS_DATABANK AD ON AD.ID_ANALYSIS=A.ID_ANALYSIS LEFT JOIN DATABANK D ON AD.ID_DATABANK=D.ID_DATABANK WHERE ID_EXPERIMENT=$experimentID AND (D.DECOY_TAG='yes' or D.DECOY_TAG LIKE '%rev%')") or die "Couldn't prepare statement: " . $dbh->errstr;
    $sthSampleList->execute;
    if ($sthSampleList->rows==0){ print "No sample generated with a decoy databank found.";}
    else{
        my $liClass='liRow1';
        print "<UL valign=top class=\"checklist\">";
        my $lastID;
        my $i=0;
        while (my($sampleName,$sampleID,$analysisName,$analysisID,$decoyTag)=$sthSampleList->fetchrow_array){
            if ($i == 0){
                print "<LI class=\"$liClass\"><LABEL FOR=\"$sampleID\"><B>$sampleName\n</B><BR>";
                print "\t<INPUT type=\"checkbox\" name=\"analysisFiles\" value=\"$analysisID\/$analysisName\"><LABEL FOR=\"$analysisID\">$analysisName<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;DECOY=$decoyTag</LABEL>\n<BR>";
            }
            elsif ($sampleID==$lastID) {
                print "\t<INPUT type=\"checkbox\" name=\"analysisFiles\" value=\"$analysisID\/$analysisName\"><LABEL FOR=\"$analysisID\">$analysisName<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;DECOY=$decoyTag</LABEL>\n<BR>";
            }
            else{
                print "</LABEL><BR></LI>";
                $liClass=($liClass eq 'liRow1')? 'liRow2' : 'liRow1';
                print "<LI class=\"$liClass\"><LABEL FOR=\"$sampleID\"><B>$sampleName\n</B><BR>";
                print "\t<INPUT type=\"checkbox\" name=\"analysisFiles\" value=\"$analysisID\/$analysisName\"><LABEL FOR=\"$analysisID\">$analysisName<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;DECOY=$decoyTag</LABEL>\n<BR>";
            }
            $lastID=$sampleID;
            $i++;
        }
        print "</LABEL><BR></LI>";
        print "</UL>";
    }
    $sthSampleList->finish;
}

sub ajaxRTData{
	my $rtID=param('rtID');
	my $i=1;
	my $rtData=$dbh->selectrow_array("SELECT DATA FROM REFERENCE_RT WHERE ID_REFERENCE_RT=$rtID") or die "Couldn't prepare statement: " . $dbh->errstr;
	print "<TABLE><TR><TH align=\"left\">#</TH><TH align=\"left\">Sequence</TH><TH align=\"left\">Mass</TH><TH align=\"left\">iRT</TH></TR>";
	my $xml =XML::Simple-> new (KeepRoot=>1);
	my $xmlRTData = $xml->XMLin($rtData);
	foreach my $data (@{$xmlRTData->{PEPTIDE_DATA}->{PEPTIDE}}) {
		my $pepSeq=$data->{sequence};
		my $pepMass=$data->{monoMass};
		my $pepIrt=$data->{iRT};
		if (! $data->{excluded}) {
            print "<TR><TD>$i&nbsp&nbsp</TD><TD><B>$pepSeq</B></TD><TD>&nbsp$pepMass</TD><TD>&nbsp$pepIrt</TD></TR>";
        }
		else{
			print "<TR><TD><FONT color=\"#A4A4A4\">$i&nbsp;&nbsp;</FONT></TD><TD><FONT color=\"#A4A4A4\"><B>$pepSeq</B></FONT></TD><TD><FONT color=\"#A4A4A4\">&nbsp;$pepMass</FONT></TD><TD><FONT color=\"#A4A4A4\">&nbsp$pepIrt&nbsp;&nbsp;(excluded)</FONT></TD></TR>";
		}
		$i++;
	}
	print "</TABLE>";
}

sub ajaxSelectSpecieDB{
	my ($commonName,$scientName)=split(/_/,param('species'));
	($scientName)=quotemeta($scientName);
	($commonName)=quotemeta($commonName);
	my $sthDBSpeciesList=$dbh->prepare("SELECT D.NAME, D.ID_DATABANK,D.FASTA_FILE,D.DECOY_TAG,DT.NAME,D.ORGANISM FROM DATABANK D,DATABANK_TYPE DT WHERE D.DECOY_TAG='yes' or D.DECOY_TAG LIKE '%rev%' AND D.USE_STATUS='yes' AND D.ID_DBTYPE=DT.ID_DBTYPE AND (D.ORGANISM='$commonName' OR D.ORGANISM='$scientName' OR D.ORGANISM='')") or die "Couldn't prepare statement: " . $dbh->errstr;
	$sthDBSpeciesList->execute;

	if ($sthDBSpeciesList->rows==0) {
		print "No databank available for this organism.";
	}
	else{
		my %dbHash;
		while (my($dbName,$dbID,$dbFileName,$decoyTag,$dbTypeName,$organismDB)=$sthDBSpeciesList->fetchrow_array){
			next if ($dbFileName=~m/:/) ;
			next if ($decoyTag eq "No" || $decoyTag eq "");
			my $matchSpecies;
			if ($organismDB eq '' || !$organismDB) {$matchSpecies=0;}
            else{$matchSpecies=1;}
			$dbHash{$matchSpecies}{"$dbName&nbsp;\[$dbTypeName\]"}=$dbID;
		}
		
		print qq
		|&nbsp;<SELECT name ="dbfile" >
			<option value="">-= Select Databank =-</option>
		|;
		foreach my $optgroup (sort {$a cmp $b} keys %dbHash){
			if ($optgroup == 0) {
                print "<OPTGROUP label=\"Databank without species name.\">";
            }
            else{print "<OPTGROUP label=\"$scientName\">";}

			foreach my $dbName (sort {lc $a cmp lc $b} keys %{$dbHash{$optgroup}}){
				print "<option value=\"$dbHash{$optgroup}{$dbName}\">$dbName</option>";
			}
			print "</OPTGROUP>";
		}
		print "</SELECT>";
	}
	$sthDBSpeciesList->finish;
}

sub selectDBMergeLib{
	my $dbTypeID=$dbh->selectrow_array("SELECT ID_DBTYPE FROM DATABANK D,DATABANK_SWATHLIB DS WHERE ID_SWATH_LIB=$libraryID AND D.ID_DATABANK=DS.ID_DATABANK AND DS.DB_RANK=1") or die "Couldn't prepare statement: " . $dbh->errstr;
	my $sthDBUpdateList=$dbh->prepare("SELECT D.NAME, D.ID_DATABANK,D.FASTA_FILE,D.DECOY_TAG,DT.NAME FROM DATABANK D,DATABANK_TYPE DT WHERE D.DECOY_TAG='yes' or D.DECOY_TAG LIKE '%rev%' AND D.USE_STATUS='yes' AND DT.ID_DBTYPE=$dbTypeID AND D.ID_DBTYPE=DT.ID_DBTYPE ") or die "Couldn't prepare statement: " . $dbh->errstr;
    $sthDBUpdateList->execute;
    if ($sthDBUpdateList->rows==0){ print "No sample generated with a decoy databank found.";}
    else{
		my %dbHash;
		while (my($dbName,$dbID,$dbFileName,$decoyTag,$dbTypeName)=$sthDBUpdateList->fetchrow_array){
			next if ($dbFileName=~m/:/) ;
			next if ($decoyTag eq "No" || $decoyTag eq "");
			$dbHash{"$dbName&nbsp;\[$dbTypeName\]"}=$dbID;
		}
		
		print qq
		|&nbsp;<SELECT name ="dbfile">
			<option value="">-= Select Databank =-</option>
		|;
		foreach my $dbName (sort {lc $a cmp lc $b} keys %dbHash){
			print "<option value=\"$dbHash{$dbName}\">$dbName</option>";
		}
		print "</SELECT>";
	}
	$sthDBUpdateList->finish;
}



$dbh->disconnect;


####>Revision history<#####
# 1.8.9 [FIX] Fix output and error path of child processes (VS 08/01/20)
# 1.8.8 [ENHANCEMENT] Add PTMProphet to spectral library building (VS 18/11/19)
# 1.8.7 Handles new job monitoring system (VS 10/05/19)
# 1.8.6 Generalize SQL Queries to avoid specific iRT and Decoy filtering (VS 10/05/19)
# 1.8.5 Add acetylation recovering on importation (VS 22/10/2018) 
# 1.8.4 Add option to create .mrm library final file (MLP 19/17/18)
# 1.8.3 Minor modif : editting a library's name and PBSerror is now handled by &clusterInfo{checkError} (MLP 16/04/18)
# 1.8.2 Minor modif to manage cluster error (PBSerror.txt) (MLP 11/04/18)
# 1.8.1 Minor modif to manage cluster error (PBSerror.txt) (MLP 11/04/18)
# 1.8.0 Replace monitorLibrariesCreation by monitorDIAProcess (MLP 08/02/2018)
# 1.7.9 Minor modif to increase the time before temporary folder deletion (MLP 05/02/18)
# 1.7.8 Modif on job memory for merged swath librarie (MLP 02/02/18)
# 1.7.7 Minor modif to allow monitoring for creation of merged swath libraries (MLP 31/01/18)
# 1.7.6 Allow peptide modification errors (MLP 30/01/18)
# 1.7.5 Add monitoring Swath libraries creation (MLP 25/01/18)
# 1.7.4 Modif on database insertion to print more caracters to avoid timeout (MLP 11/01/18)
# 1.7.3 Add an option to archive a library (MLP 10/01/18)
# 1.7.2 New scripts for library update (MLP 08/01/18)
# 1.7.1 Added sleep in job-wait loops to prevent disk access overload & other minor changes (PP 26/12/17)
# 1.7.0 Modif on library update (save non consensus files) (MLP 20/12/17)
# 1.6.8 Minor modif (MLP 19/12/17)
# 1.6.7 Parallelisation of dat conversion and tandem modification (MLP 13/12/17)
# 1.6.6 Minor modif (MLP 08/12/17)
# 1.6.5 Add &promsConfig::getTPPModificationCode to convert modifications that are not on the available form for export (GlyGly => 242) (MLP 24/11/17)
# 1.6.4 Minor modif (MLP 17/11/17)
# 1.6.3 Minor modif to store the number of modified peptides (for GG or phospho) (MLP 15/11/17)
# 1.6.2 Modif on PBSerror to filter getcwd error (MLP 09/11/17)
# 1.6.1 Minor modif (MLP 08/11/17)
# 1.6.0 Add shared directory option to upload data files (MLP 08/11/17)
# 1.5.0 Modif to upload a data files archive and to launch job on CentOS (MLP 07/11/17)
# 1.4.3 Modif to allow the selection of FDR type (mFDR,pepFDR,protFDR) (MLP 24/10/17)
# 1.4.2 Minor modif to merge libraries (MLP 19/10/17)
# 1.4.1 Minor modif to select the acquisition data instrument (MLP 19/07/17) 
# 1.4.0 Minor modif (MLP 11/05/17)
# 1.3.9 Minor modification for edit form (MLP 28/02/2017)
# 1.3.8 Added XML format (PAPPSO X! Tandem) for input files (MLP 16/02/17) 
# 1.3.7 Minor modification (MLP 12/01/17)
# 1.3.6 Update to allow multiple databank library (MLP 05/12/2016)
# 1.3.5 Add explicit error if there are no Decoys in .dat or in .tandem.pep.xml files and add modifications for PARAM_STRG (table SWATH_LIB)  (MLP 29/11/2016)
# 1.3.4 Minor modifications (MLP 25/10/2016)
# 1.3.3 Add &promsMod::cleanParameters verification (MLP 28/09/2016)
# 1.3.3 Minor bug in file number test (MLP 27/09/2016)
# 1.3.2 Check the number and names of selected files and add restore option (MLP 21/09/2016)
# 1.3.1 Bug fix in identifier match test & restriction on species list (PP 17/08/16)
# 1.3.0 Minor corrections to run bash script on web server (MLP 12/08/2016)
# 1.2.9 Update 'deleteSwathLib' to totaly delete libraries that are not used (MLP 22/07/2016)
# 1.2.8 Verification of the compatibility between input files identifier type and databank type (MLP 28/06/2016)
# 1.2.7 Add num prot in table SWATH_LIB (MLP 18/05/2016)
# 1.2.6 Minor modifications (MLP 11/05/2016)
# 1.2.5 Minor modifications to increase asked memory (MLP 10/05/2016)
# 1.2.4 Minor modifications to print wait time (MLP 09/05/2016)
# 1.2.3 Minor modifications to allow to merge two librairies (MLP 02/05/2016)
# 1.2.2 Minor modifications to delete temporary folder (MLP 29/04/2016)
# 1.2.1 Minor update to insert peptides modifications in database in SWATH_LIB_MODIFICATION (MLP 28/04/2016)
# 1.2.0 Update to run job on cluster (MLP 25/04/2016)
# 1.1.0 Add split mode (MLP 19/04/2016)
# 1.0.2 Minor corrections (MLP 15/04/2016)
# 1.0.1 Change Swath_Lib path ($promsPath{data}/Swath_Lib to $promsPath{swath_lib}) (MLP 13/04/2016)
# 1.0.0 Allow to create and update swath library (MLP 11/04/2016)
