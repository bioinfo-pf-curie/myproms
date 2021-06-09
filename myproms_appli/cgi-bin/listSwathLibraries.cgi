#!/usr/local/bin/perl -w

################################################################################
# listSwathLibraries.cgi         1.2.2                                         #
# Authors: M. Le Picard, V. Sabatet (Institut Curie)                           #
# Contact: myproms@curie.fr                                                    #
# Lists the libraries available in myProMS                                     #
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
use promsMod;

my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my $bgColor=$lightColor;
my %promsPath=&promsConfig::getServerInfo;
my $userID = $ENV{'REMOTE_USER'};

##########################
####>Connecting to DB<####
##########################
print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);

my $dbh = &promsConfig::dbConnect;
my $act = param("ACT");

if($act) {
    if($act eq "showLibProteinContent") {
        my ($libID, $libName, $stats) = &promsMod::cleanParameters(param("libID"), param("libName"), param("stats"));
        if($libID and $libName and $stats) {
            showLibProteinContent($libID, $libName, $stats);
        }
        exit;
    } elsif($act eq "showAssociatedProject") {
        my $libID = &promsMod::cleanParameters(param("libID"));
        showAssociatedProject($libID) if($libID);
        exit;
    }
}

my $sthDataLibrary = $dbh->prepare("SELECT SL.ID_SWATH_LIB,RT.NAME,SL.NAME,SL.DES,SL.IDENTIFIER_TYPE,SL.ORGANISM,SL.VERSION_NAME,IF(SL.PARAM_STRG LIKE '%Spectronaut%', 'SPC', 'TPP') AS SOFTWARE,SL.START_DATE,SL.SPLIT,SL.STATISTICS,SL.USE_STATUS FROM SWATH_LIB SL LEFT JOIN REFERENCE_RT RT ON SL.ID_REFERENCE_RT=RT.ID_REFERENCE_RT WHERE SL.USE_STATUS!='no' AND SL.USE_STATUS!='err' ORDER BY SL.NAME ASC");
$sthDataLibrary->execute;
my $refDataLib = $sthDataLibrary->fetchall_arrayref;
$sthDataLibrary->finish;

#######################
#### Starting HTML ####
#######################

print qq |
<HTML>
    <HEAD>
        <TITLE>List of spectral libraries</TITLE>
        <LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
        <STYLE type="text/css">
            TD {font-size:13px;}
            .buttonlib{font-size:18px;font-weight:bold;}
            .buttonmerge{font-size:18px;font-weight:bold;margin-left:73%;margin-top:2%;position:fixed;}
            .tablelib{margin-top:4%;}
        </STYLE>
        <SCRIPT LANGUAGE="JavaScript">
            function getXMLHTTP() {
                var xhr=null;
                if(window.XMLHttpRequest) {// Firefox & others
                    xhr = new XMLHttpRequest();
                } else if(window.ActiveXObject){ // Internet Explorer
                    try {
                      xhr = new ActiveXObject("Msxml2.XMLHTTP");
                    } catch (e) {
                        try {
                            xhr = new ActiveXObject("Microsoft.XMLHTTP");
                        } catch (e1) {
                            xhr = null;
                        }
                    }
                } else { // XMLHttpRequest not supported by browser
                    alert("Your browser does not support XMLHTTPRequest objects...");
                }
                return xhr;
            }
            
            function deleteSwathLibrary(libraryID, name) {
                if (confirm('Delete library: '+name+'?')) {
                    window.location = './editSwathLibrary.cgi?ACT=delete&ID=' + libraryID;
                }
            }
            
            function archiveSwathLibrary(libraryID, name) {
                if (confirm('Archive library: ' + name + '?')) {
                    window.location = './editSwathLibrary.cgi?ACT=archive&ID=' + libraryID;
                }
            }
            
            function monitorLibCreation() {
                var monitorJobsWin=window.open("$promsPath{cgi}/monitorJobs.cgi?filterType=Spectral Library&filterDateNumber=1&filterDateType=DAY",'monitorJobsWindow','width=1200,height=500,scrollbars=yes,resizable=yes');
                monitorJobsWin.focus();
            }
            
            function showLibProteinContent(libraryID, libraryName, stats) {
                var libraryProteinContentDiv = document.getElementById('libraryProteinContent'+libraryID);
                var XHR = getXMLHTTP();
                libraryProteinContentDiv.innerHTML = '<IMG src="$promsPath{images}/scrollbarGreen.gif" \><br/>';
                XHR.onreadystatechange = function () {
                    if(XHR.readyState==4 && XHR.responseText ) {
                        libraryProteinContentDiv.innerHTML = XHR.responseText;
                    }
                }
                XHR.open("GET","$promsPath{cgi}/listSwathLibraries.cgi?ACT=showLibProteinContent&libID=" + libraryID + "&libName=" + libraryName + "&stats=" + stats, true);
                XHR.send(null);
            }
            
            function showAssociatedProjects(libraryID) {
                var libraryAssocProjectsContentDiv = document.getElementById('libraryAssocProjects'+libraryID);
                var XHR = getXMLHTTP();
                libraryAssocProjectsContentDiv.innerHTML = '<IMG src="$promsPath{images}/scrollbarGreen.gif" \><br/>';
                XHR.onreadystatechange = function () {
                    if(XHR.readyState==4 && XHR.responseText ) {
                        libraryAssocProjectsContentDiv.innerHTML = XHR.responseText;
                    }
                }
                XHR.open("GET","$promsPath{cgi}/listSwathLibraries.cgi?ACT=showAssociatedProject&libID=" + libraryID, true);
                XHR.send(null);
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
        <CENTER>
        <FONT class="title1">List of spectral libraries<BR></FONT><BR><BR><BR>
        
        <INPUT type="button" class="buttonlib" value="Add new spectral library" onclick="window.location='./editSwathLibrary.cgi?ACT=add'">
        <INPUT type="button" class="buttonlib" value="Merge two spectral libraries" onclick="window.location='./editSwathLibrary.cgi?ACT=merge'">
        <INPUT type="button" class="buttonlib" value="Monitor spectral libraries" onclick="monitorLibCreation()"><br/><br/><br/>
        
        <table>
    |;
    
    foreach my $dataRow (@{$refDataLib}) {
        my ($libraryID, $rtFile, $libraryName, $des, $identifierType, $organism, $versionName, $software, $startDate, $modeSplit, $stat, $useStatus)=@{$dataRow};
        my ($version, $dbFile);
        if ($versionName) {
            ($version=$versionName)=~s/v//;
        }
        
        $des='' unless $des;
        $stat =~ s/[\n\s]//g;
        my $split = (!$modeSplit) ? $modeSplit : ($modeSplit == 1) ? "Split" : "Unsplit";
        my $workDir = "$promsPath{swath_lib}/SwLib_$libraryID";
        if ((-e "$workDir/$libraryName\_peakview.tsv") && ((-M "$workDir/$libraryName\_peakview.tsv") >1)) {
            system "rm $workDir/$libraryName\_peakview.tsv";
        }
        if ((-e "$workDir/$libraryName\_peakview.txt") && ((-M "$workDir/$libraryName\_peakview.txt") >1)) {
            system "rm $workDir/$libraryName\_peakview.txt";
        }
        if ((-e "$workDir/$libraryName\_openswath.csv") && ((-M "$workDir/$libraryName\_openswath.csv") >1)) {
            system "rm $workDir/$libraryName\_openswath.csv";
        }
        if ((-e "$workDir/$libraryName\_param") && ((-M "$workDir/$libraryName\_param") >1)) {
            system "rm $workDir/$libraryName\_param";
        }
        if (-e "$workDir/sortieExport.txt") {
            system "rm $workDir/sortieExport.txt";
        }
        
        my $sthDatabankSwathLib = $dbh->prepare("SELECT D.NAME FROM DATABANK_SWATHLIB DS,DATABANK D WHERE DS.ID_DATABANK=D.ID_DATABANK AND DS.ID_SWATH_LIB=?");
        $sthDatabankSwathLib->execute($libraryID);
        while (my $db = $sthDatabankSwathLib->fetchrow_array) {
            $dbFile .= "; " if($dbFile);
            $dbFile .= "$db";
        }
        
        print qq |
            <tr id='lib$libraryID' bgcolor='$bgColor' style='margin-bottom:10px'><td>
                <table>
                    <td id='libInfos$libraryID' style='width:600px'>
                        <p style='font-size:18px;font-weight:bold;margin:3px 0 0 12px;'>$libraryName
                        
        |;
        print "<FONT color=\"red\" style=\"font-size:18px;font-weight:bold;\">&nbsp;(Archived)</FONT>" if ($useStatus eq 'arc');
        print "</p><p style='margin-left:12px;'>";
        print "<B>Version:</B> $versionName<br/>" if $versionName;
        if($split) { print "<B>Mode:</B> $split<br/>"; }
        print "<B>Identifier type:</B> $identifierType<br/>" if $identifierType;
        print qq |
               <B>Database(s):</B> $dbFile<br/>
        |;
        if($rtFile) { print "<B>RT:</B> $rtFile<br/>"; }
        if ($des) { print "<B>Description:</B> $des<br/>";}
        if ($organism){print "<B>Organism:</B> $organism<br/>";}
        if ($startDate){print "<B>Creation date:</B> $startDate<br/>";}
        print qq |
                        </p>
                    </td>
            
                    <td id="libActions$libraryID" style='width:225px;'>
                        <INPUT type="button" value="Delete" style="width:100px" onclick="deleteSwathLibrary($libraryID,'$libraryName')">
                        <INPUT type="button" value="Edit" style="width:100px" onclick="window.location='./editSwathLibrary.cgi?ACT=edit&ID=$libraryID'">
                        <INPUT type="button" value="Export library" style="width:100px" onclick="window.location='./exportSwathLibrary.cgi?ID=$libraryID'">
                        <INPUT type="button" value="Export proteins" style="width:100px" onclick="window.location='./exportLibraryProteins.cgi?ID=$libraryID'">
                        <INPUT type="button" value="Search" style="width:100px" onclick="window.location='./searchSwathLibrary.cgi?ID=$libraryID'">
                
        |;
        
        if (-e "$workDir/Lib$libraryID\_v$version.tar.gz" || $modeSplit==1) {
            print "<INPUT type=\"button\" value=\"Update\" style=\"width:100px\" onclick=\"window.location='./editSwathLibrary.cgi?ACT=update&ID=$libraryID'\">&nbsp;";
        }
        my $prevVersion=$version-1;
        if (-e "$workDir/Lib$libraryID\_v$version.tar.gz"){
            print "<INPUT type=\"button\" value=\"Archive\" style=\"width:100px\" onclick=\"archiveSwathLibrary($libraryID,'$libraryName')\">";
        }
        
        if ($version >1 and -e "$workDir/Lib$libraryID\_v$prevVersion.tar.gz" and $software eq 'TPP') {
            print "<INPUT type=\"button\" value=\"Restore previous version\" style=\"width:175px\" onclick=\"window.location='./editSwathLibrary.cgi?ACT=restore&ID=$libraryID'\"><BR>";
        }
        
        print qq |
                    </td>
                </tr>
                
                <tr>
                    <td style="font-weight:bold;" colspan=2>
                        <div id="libraryProteinContent$libraryID" style="margin:5px 0 0 18px;">
                            <INPUT type="button" value="Library statistics" style="width:161px" onclick="showLibProteinContent('$libraryID', '$libraryName', '$stat')">
                        </div>
            
                        <div id="libraryAssocProjects$libraryID" style="margin:5px 0 10px 18px;">
                            <INPUT type="button" value="Associated project(s)" style="width:161px" onclick="showAssociatedProjects('$libraryID')">
                        </div>
                    </td>
                </tr>
            </table>
        </td></tr>
        |;
          
        $bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
    }
    if (scalar @{$refDataLib} == 0) {
        print "<TR><TH colspan=3>(No SWATH libraries in use)</TH></TR>\n";
    }
    print qq
    |</TABLE>
    </CENTER>
</BODY>
</HTML>
|;
$dbh->disconnect;


sub showAssociatedProject {
    my ($libraryID, $libraryName) = @_;
    my $sthLibProj =$dbh->prepare("SELECT P.OWNER, P.ID_PROJECT, P.NAME, E.NAME, COUNT(A.ID_ANALYSIS) AS COUNT_ANA FROM ANALYSIS_SWATH_LIB ASL INNER JOIN ANALYSIS A ON A.ID_ANALYSIS=ASL.ID_ANALYSIS INNER JOIN SAMPLE S ON S.ID_SAMPLE=A.ID_SAMPLE INNER JOIN EXPERIMENT E ON E.ID_EXPERIMENT=S.ID_EXPERIMENT INNER JOIN PROJECT P ON P.ID_PROJECT=E.ID_PROJECT WHERE ASL.ID_SWATH_LIB=? GROUP BY P.OWNER, P.ID_PROJECT, E.ID_EXPERIMENT ORDER BY P.OWNER, P.NAME, E.NAME");
    my ($displayData, $currentUserName, $currentProjectName) = ('', '', '');
    
    $sthLibProj->execute($libraryID);
    while(my ($userName, $projectID, $projectName, $expName, $countAna) = $sthLibProj->fetchrow_array) {
        $displayData .= qq |
            <tr>
                <td align='center'>$userName</td>
                <td align='center'>$projectName</td>
                <td align='left'>$expName (used in $countAna analysis)</td>
            </tr>|;
    }
    
    if($displayData) {
        $displayData = qq |
            <table>
                <thead>
                    <th style='max-width: 200px;min-width: 170px;' align='center'>User</th>
                    <th style='max-width: 250px;min-width: 200px;' align='center'>Project</th>
                    <th style='max-width: 350px' align='center'>Experiment</th>
                </thead>
                <tbody>
                    $displayData
                </tbody>
            </table>
        |;
    } else {
        $displayData = "No related projects";
    }
    
    print qq |
        <br/><fieldset>
            <legend><B>Related Experiments:</B></legend>
            <div id="projectsLib:$libraryID" style="width: 100%; max-height: 200px !important; overflow-y: auto; overflow-x: hidden;">
                $displayData
            </div>
        </fieldset>
    |;
}


sub showLibProteinContent {
    my ($libraryID, $libraryName, $stat) = @_;
    my ($xmlStat, $numProt, $numPep, $numProtSpe, $numPepModList);
    my $xml =XML::Simple-> new (KeepRoot=>1);
    
    if ($stat) {
        $xmlStat = $xml->XMLin($stat);
        $numProt = $xmlStat->{'NUM_ENTRY'}->{'NUM_PROT'};
        $numPep = $xmlStat->{'NUM_ENTRY'}->{'NUM_PEP'};
        $numProtSpe = $xmlStat->{'NUM_ENTRY'}->{'NUM_PROT_SPECIFIC'};
        $numPepModList = $xmlStat->{'NUM_ENTRY'}->{'NUM_PEP_MOD'} if($xmlStat->{'NUM_ENTRY'}->{'NUM_PEP_MOD'});
    }
    
    unless($numPepModList){
        my (%modifLib, %numPepMod);
        my $sthModSwathLib=$dbh->prepare("SELECT SLM.ID_MODIFICATION, INTERIM_NAME FROM SWATH_LIB_MODIFICATION SLM INNER JOIN MODIFICATION M ON SLM.ID_MODIFICATION=M.ID_MODIFICATION WHERE ID_SWATH_LIB=?");

        $sthModSwathLib->execute($libraryID);
        while (my ($modID, $modName)=$sthModSwathLib->fetchrow_array){
            next if ($modName ne 'GG' && $modName ne 'Phospho' && $modName ne 'Acetyl');
            $modName='GlyGly' if $modName eq 'GG';
            $modifLib{$modName}{'Single'}=0;
            $modifLib{$modName}{'Multiple'}=0;
        }
       
        if (%modifLib and -s "$promsPath{swath_lib}/SwLib_$libraryID/$libraryName.sptxt") {
            open(LIB,"<$promsPath{swath_lib}/SwLib_$libraryID/$libraryName.sptxt");
            while (<LIB>){
                if ($_=~/Mods=(\d)\/(\S*)/){
                    my $modifList=$2;
                    my @modList=split(/\//,$modifList);
                    foreach my $mod (@modList){
                        my ($pos,$aa,$varModName)=split(/,/,$mod);
                        $numPepMod{$varModName}{'Multiple'}++ if ($varModName eq 'GlyGly' || $varModName eq 'Phospho' || $varModName eq 'Acetyl');
                    }
                    $numPepMod{'GlyGly'}{'Single'}++  if $modifList=~/GlyGly/;
                    $numPepMod{'Phospho'}{'Single'}++  if $modifList=~/Phospho/;      
                    $numPepMod{'Acetyl'}{'Single'}++  if $modifList=~/Acetyl/;
                }
            }
            close LIB;
        }
        
        if(%numPepMod){
            foreach my $mod (keys %numPepMod){
                $numPepModList.='&' if $numPepModList;
                $numPepModList.="$mod".':'.$numPepMod{$mod}{'Single'}.'/'.$numPepMod{$mod}{'Multiple'};
            }
        }
        
    }
    
    my $displayData = '';
    if ($numProt) { $displayData .= "<B>Number of proteins:</B> $numProt";}
    if ($numProtSpe) { $displayData .= "<BR><B>Number of unambigious proteins:</B> $numProtSpe";}
    if ($numPep) { $displayData .= "<BR><B>Number of peptides:</B> $numPep";}
    if ($numPepModList){
        my @modList=split(/\.\./,$numPepModList);
        foreach my $mod (@modList) {
            my ($modName,$numPepMod,$numModSite)=($mod=~/(.+):(\d+)\/(\d+)/);
            $displayData .= "<BR><B>Number of $modName sites:</B> $numModSite";
            $displayData .= "<BR><B>Number of peptides with $modName:</B> $numPepMod";
        }
    }
    
    $displayData = "No proteins information found" if(!$displayData);
    if($displayData) {
        print qq |
            <fieldset>
                <legend><B>Protein content:</B></legend>
                <div id="protInfoLib:$libraryID" style="width: 100%">
                    $displayData
                </div>
            </fieldset>
        |;
    }
}

####>Revision history<####
# 1.2.2 [CHANGE] Only archive TPP libraries (VS 22/10/20)
# 1.2.1 [CHANGE] Use new job monitoring window opening parameters (VS 18/11/19)
# 1.2.0 [ENHANCEMENT] Possibility to display projects and experiments associated to a library (VS 13/11/19)
# 1.1.11 [MODIF] Switch from monitorDIA to monitorJobs script (VS 21/10/19)
# 1.1.10 Add "Export proteins" option (VL 20/09/19)
# 1.1.9 Speed up page loading using asynchronous protein content computing (VS 25/03/19)
# 1.1.8 Add acetyl count displaying (VS 22/10/2018)
# 1.1.7 Minor modif to increase the time before temporary folder deletion (MLP 05/02/18)
# 1.1.6 Add monitoring libraries creation (MLP 25/01/18)
# 1.1.5 Minor modif (MLP 10/01/18)
# 1.1.4 Add an option to archive a library (MLP 10/01/18)
# 1.1.3 Minor modif (MLP 20/12/17)
# 1.1.2 Modif to display the number of modified peptides (for GG or phospho) (MLP 15/11/17)
# 1.1.1 Minor modif (MLP 22/03/17)
# 1.1.0 Added modifications to delete older tmp folders (MLP 15/12/2016)
# 1.0.9 Changing "SWATH libraries" into "spectral libraries" (MLP 24/11/2016)
# 1.0.8 Added "Restore" option to return to the previous library version (MLP 23/09/2016)
# 1.0.7 Added $sthDatabankSwathLib for Swath lib with different databanks (MLP 09/08/2016)
# 1.0.6 Minor modification of query names ($dataLibrary to $sthDataLibrary) (MLP 22/06/2016)
# 1.0.5 Remove entry NUM_ENTRY in table SWATH_LIB (MLP 03/06/2016)
# 1.0.4 Display number of protein and peptide (MLP 19/05/2016)
# 1.0.3 Minor modification to sort libraries list by library name (MLP 10/05/2016)
# 1.0.2 Added "Merge two Swath libraries" and "Edit" buttons (MLP 02/05/2016)
# 1.0.1 Added the suppression of export files after one day (MLP 20/04/2016)
# 1.0.0 Create listSwathLibraries.cgi to allow to list swath libraries (MLP 11/04/2016)


