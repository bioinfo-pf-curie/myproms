#!/usr/local/bin/perl -w

################################################################################
# listSwathLibraries.cgi         1.1.7                                         #
# Authors: M. Le Picard (Institut Curie)                                       #
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
my $userID=$ENV{'REMOTE_USER'};

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

my $sthDataLibrary=$dbh->prepare("SELECT SL.ID_SWATH_LIB,RT.NAME,SL.NAME,SL.DES,SL.IDENTIFIER_TYPE,SL.ORGANISM,SL.VERSION_NAME,SL.START_DATE,SL.SPLIT,SL.STATISTICS,SL.USE_STATUS FROM SWATH_LIB SL,REFERENCE_RT RT WHERE SL.USE_STATUS!='no' AND SL.USE_STATUS!='err' AND SL.ID_REFERENCE_RT=RT.ID_REFERENCE_RT ORDER BY SL.NAME ASC");
$sthDataLibrary->execute;
my $refDataLib=$sthDataLibrary->fetchall_arrayref;
$sthDataLibrary->finish;
my $sthDatabankSwathLib=$dbh->prepare("SELECT D.NAME FROM DATABANK_SWATHLIB DS,DATABANK D WHERE DS.ID_DATABANK=D.ID_DATABANK AND DS.ID_SWATH_LIB=?");
my $sthModSwathLib=$dbh->prepare("SELECT SLM.ID_MODIFICATION,PSI_MS_NAME FROM SWATH_LIB_MODIFICATION SLM, MODIFICATION M WHERE SLM.ID_MODIFICATION=M.ID_MODIFICATION AND ID_SWATH_LIB=?");

##>delete older tmp folders
foreach my $folder (glob "$promsPath{data}/tmp/Swath/*") {     
    if (-M $folder && -M $folder >6) {
        system "rm -r $folder";
    }
}


#######################
####>Starting HTML<####
#######################
print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);

print qq
|<HTML>
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
function deleteSwathLibrary(libraryID,name) {
	if (confirm('Delete library: '+name+'?   ')) {
		window.location='./editSwathLibrary.cgi?ACT=delete&ID='+libraryID;
	}
}
function archiveSwathLibrary(libraryID,name) {
	if (confirm('Archive library: '+name+'?   ')) {
		window.location='./editSwathLibrary.cgi?ACT=archive&ID='+libraryID;
	}
}
function monitorLibCreation() {
    var monitorLibWindows=window.open("$promsPath{cgi}/monitorDIAProcess.cgi?ACT=library",'Monitoring libraries creation','width=1000,height=500,scrollbars=yes,resizable=yes');
    monitorLibWindows.focus();
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
    <INPUT type="button" class="buttonlib" value="Monitor spectral libraries" onclick="monitorLibCreation()">
    <TABLE border=0 cellspacing=0 class="tablelib" >
    |;
    
    foreach my $dataRow (@{$refDataLib}) {
        my ($libraryID,$rtFile,$libraryName,$des,$identifierType,$organism,$versionName,$startDate,$modeSplit,$stat,$useStatus)=@{$dataRow};
        my $version;
        if ($versionName) {
            ($version=$versionName)=~s/v//;
        }
        
        my $split;
        $des='' unless $des;
        if ($modeSplit == 1) {
            $split="Split";
        }
        else{
            $split="Unsplit";
        }
        my $xml =XML::Simple-> new (KeepRoot=>1);
        my ($xmlStat,$numProt,$numPep,$numProtSpe,$numPepModList);
        if ($stat) {
            $xmlStat = $xml->XMLin($stat);
            $numProt=$xmlStat->{'NUM_ENTRY'}->{'NUM_PROT'};
            $numPep=$xmlStat->{'NUM_ENTRY'}->{'NUM_PEP'};
            $numProtSpe=$xmlStat->{'NUM_ENTRY'}->{'NUM_PROT_SPECIFIC'};
            $numPepModList=$xmlStat->{'NUM_ENTRY'}->{'NUM_PEP_MOD'} if $xmlStat->{'NUM_ENTRY'}->{'NUM_PEP_MOD'};
        }
        my $workDir="$promsPath{swath_lib}/SwLib_$libraryID";
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
        $sthDatabankSwathLib->execute($libraryID);
        my $dbFile;
        while (my $db=$sthDatabankSwathLib->fetchrow_array) {
            $dbFile.="$db\t";
        }
        
        
        unless($numPepModList){
            my %modifLib;
            $sthModSwathLib->execute($libraryID);
            while (my ($modID,$modName)=$sthModSwathLib->fetchrow_array){
                next if ($modName ne 'GG' && $modName ne 'Phospho');
                $modName='GlyGly' if $modName eq 'GG';
                $modifLib{$modName}{'Single'}=0;
                $modifLib{$modName}{'Multiple'}=0;
            }
           
            my %numPepMod;
            if (%modifLib){
                open(LIB,"<$promsPath{swath_lib}/SwLib_$libraryID/$libraryName.sptxt");
                while (<LIB>){
                    if ($_=~/Mods=(\d)\/(\S*)/){
                        my $modifList=$2;
                        my @modList=split(/\//,$modifList);
                        foreach my $mod (@modList){
                            my ($pos,$aa,$varModName)=split(/,/,$mod);
                            $numPepMod{$varModName}{'Multiple'}++ if ($varModName eq 'GlyGly' || $varModName eq 'Phospho');
                        }
                        $numPepMod{'GlyGly'}{'Single'}++  if $modifList=~/GlyGly/;
                        $numPepMod{'Phospho'}{'Single'}++  if $modifList=~/Phospho/;      
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
        
        print qq
        |<TR bgcolor=$bgColor>
            <TD>&nbsp&nbsp</TD>
            <TD width=600><FONT style="font-size:18px;font-weight:bold;">$libraryName</FONT>
        |;
        if ($useStatus eq 'arc') {
            print "<FONT color=\"red\" style=\"font-size:18px;font-weight:bold;\">&nbsp;(Archived)</FONT>";
        }
        print "<BR>";
        print "<BR><B>&nbsp;&nbsp;&nbsp;Version:</B> $versionName" if $versionName;
        print "<BR><B>&nbsp;&nbsp;&nbsp;Mode:</B> $split";
        print "<BR><B>&nbsp;&nbsp;&nbsp;Identifier type:</B> $identifierType" if $identifierType;
        print "<BR><B>&nbsp;&nbsp;&nbsp;Database(s):</B> $dbFile";
        print "<BR><B>&nbsp;&nbsp;&nbsp;RT:</B> $rtFile";    
               
        if ($des) { print "<BR><B>&nbsp;&nbsp;&nbsp;Description:</B> $des";}
        if ($numProt) { print "<BR><B>&nbsp;&nbsp;&nbsp;Number of proteins:</B> $numProt";}
        if ($numProtSpe) { print "<BR><B>&nbsp;&nbsp;&nbsp;Number of unambigious proteins:</B> $numProtSpe";}
        if ($numPep) { print "<BR><B>&nbsp;&nbsp;&nbsp;Number of peptides:</B> $numPep";}
        if ($numPepModList){
            my @modList=split(/&/,$numPepModList);
            foreach my $mod (@modList) {
                my ($modName,$numPepMod,$numModSite)=($mod=~/(.+):(\d+)\/(\d+)/);
                print "<BR><B>&nbsp;&nbsp;&nbsp;Number of $modName sites:</B> $numModSite";
                print "<BR><B>&nbsp;&nbsp;&nbsp;Number of peptides with $modName:</B> $numPepMod";
            }
        }
        if ($organism){print "<BR><B>&nbsp;&nbsp;&nbsp;Organism:</B> $organism";}
        if ($startDate){print "<BR><B>&nbsp;&nbsp;&nbsp;Creation date:</B> $startDate<BR><BR>";}
        print qq |
            </TD>
            <TH width=100>
                <INPUT type="button" value="Delete" style="width:75px" onclick="deleteSwathLibrary($libraryID,'$libraryName')">
                <INPUT type="button" value="Edit" style="width:75px" onclick="window.location='./editSwathLibrary.cgi?ACT=edit&ID=$libraryID'">
                <INPUT type="button" value="Export" style="width:75px" onclick="window.location='./exportSwathLibrary.cgi?ID=$libraryID'">
                <INPUT type="button" value="Search" style="width:75px" onclick="window.location='./searchSwathLibrary.cgi?ID=$libraryID'">
                
        |;
        if (-e "$workDir/Lib$libraryID\_v$version.tar.gz" || $modeSplit==1) {
            print "<INPUT type=\"button\" value=\"Update\" style=\"width:75px\" onclick=\"window.location='./editSwathLibrary.cgi?ACT=update&ID=$libraryID'\">&nbsp;";
        }
        my $prevVersion=$version-1;
        
        if (-e "$workDir/Lib$libraryID\_v$version.tar.gz"){
            print "<INPUT type=\"button\" value=\"Archive\" style=\"width:75px\" onclick=\"archiveSwathLibrary($libraryID,'$libraryName')\">";
        }
        
        if ($version >1 and -e "$workDir/Lib$libraryID\_v$prevVersion.tar.gz") {
            print "<INPUT type=\"button\" value=\"Restore previous version\" style=\"width:175px\" onclick=\"window.location='./editSwathLibrary.cgi?ACT=restore&ID=$libraryID'\"><BR>";
        }
        print "</TH></TR><TR><TH><BR></TH></TR>";
        
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

####>Revision history<####
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

