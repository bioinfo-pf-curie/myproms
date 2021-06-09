#!/usr/local/bin/perl -w

################################################################################
# exportSwathLibrary.cgi         1.2.2                                         #
# Authors: M. Le Picard (Institut Curie)                                       #
# Contact: myproms@curie.fr                                                    #
#Export a librarie available in myProMS  for peakview or openswath             #
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
use promsMod;
use promsDIARefonte;
use String::Util qw(trim);
use File::Compare;
use POSIX qw(strftime); # to get the time

my %promsPath=&promsConfig::getServerInfo;
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my $userID=$ENV{'REMOTE_USER'};

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

#######################
####>Starting HTML<####
#######################

my $libID = param('ID');
my ($libraryName, $libVersion, $software, $identType) = $dbh->selectrow_array("SELECT NAME, VERSION_NAME, IF(PARAM_STRG LIKE '%Spectronaut%', 'SPC', 'TPP'), IDENTIFIER_TYPE FROM SWATH_LIB WHERE ID_SWATH_LIB='$libID' ") or die "Couldn't prepare statement: " . $dbh->errstr;
my $action=(param('submit'))? param('submit') : param('ACT') ? param('ACT') : '';

if($action eq 'download') {
    my $format = &promsMod::cleanParameters(param("FORMAT"));
    downloadFile("$promsPath{swath_lib}/SwLib_$libID", "$libraryName.$format");
    exit;
}


print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);



###> Form to select files
print qq
|<HTML>
<HEAD>
<TITLE>Library</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.button{
    cursor:pointer;
    border-width:1px;
    border-style:outset;
    -moz-border-radius: 4px;
    -webkit-border-radius: 4px;
    border-radius: 4px;
    background : #EDDAA8;
    background : -webkit-gradient(linear, left top, left bottom, color-stop(0%,#EDDAA8), color-stop(100%,#EDBB52));
    background : -moz-linear-gradient(top, #EDDAA8 0%, #EDBB52 100%);
    background : -webkit-linear-gradient(top, #EDDAA8 0%, #EDBB52 100%);
    background : -o-linear-gradient(top, #EDDAA8 0%, #EDBB52 100%);
    background : -ms-linear-gradient(top, #EDDAA8 0%, #EDBB52 100%);
    background : linear-gradient(top, #EDDAA8 0%, #EDBB52 100%);
    filter : progid:DXImageTransform.Microsoft.gradient( startColorstr='#EDDAA8', endColorstr='#EDBB52',GradientType=0);
}
.popup {background-color:#FFFFFF;border:solid 3px #999999;padding:5px;box-shadow:10px 10px 20px #808080;position:absolute;display:none;} 
</STYLE>
<SCRIPT LANGUAGE="JavaScript">|;
    &promsMod::popupInfo();
    print qq |
    function selectFormat(format){
        document.getElementById('sptxt').style.display="none";
        document.getElementById('tsv').style.display="none";
        
        if (format == 'sptxt' \|\| (format == 'tsv' && '$software' == 'SPC')) {
            document.getElementById(format).style.display="";
        }
        else if (format == 'zip'){
            document.getElementById('test').style.display="none";
        }
    }
</SCRIPT>
</HEAD>|;


unless ($action){
    my $help;
    if ($identType) {
        $help="(List of desired proteins's identifiants separated by ',;' or enter/space)" if $identType=~/ID/; 
        $help="(List of desired proteins's accession numbers separated by ',;' or enter/space)" if $identType=~/ACCESSION/;
        $help="(List of desired proteins separated by ',;' or enter/space)" if $identType=~/ALL/;
    }else{$help="(List of desired proteins separated by ',;' or enter/space)";}
    
    print qq |
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
        <DIV id="form" style="display:" >
            <FONT class="title1">Export <FONT class="title1" color="#DD0000">$libraryName </FONT></FONT><BR><BR><BR><BR>
            <FORM method="POST" action="./exportSwathLibrary.cgi" name="export" enctype="multipart/form-data" >
            <TABLE bgcolor="$darkColor">
                <INPUT type="hidden" name="ID" id="ID" value="$libID">
                <TR><TH align=right valign="top">Export format : </TH>
                    <TD bgcolor="$lightColor">
                        <SELECT name="-k" onChange="selectFormat(this.value)" required>
                            <option value="">-= Select format =-</option>
        |;
        
        if($software eq 'TPP') {
            print qq |
                <option value="peakview">PeakView</option>
                <option value="openswath">OpenSwath*</option>
                <option value="spectronaut">Spectronaut</option>
                <option value="sptxt">Sptxt file</option>
                <option value="zip">Zip archive</option>
            |;
        } else {
            print qq |
                <option value="tsv">tsv</option>
            |;
        }
        
        print qq |
                        </SELECT><BR>
                    </TD>
                </TR>
        |;
        
        if($software eq 'TPP') {
            print qq |
                <TR><TH align=right valign="top">Mass range of fragment ions : </TH>
                    <TD bgcolor="$lightColor">
                        Min : <INPUT type="text" name="-lmin" size="5" value="350">&nbsp;&nbsp;Max : <INPUT type="text" name="-lmax" size="5" value="2000"><BR>
                    </TD>
                </TR>
                <TR><TH align=right valign="top">Ion series and charge : </TH>
                    <TD bgcolor="$lightColor">
                        Ions : <INPUT type="text" name="-s" size="5" value="b,y" required><SMALL> (separated by ',') for example : 'b,y' </SMALL><BR>
                        Charge : <INPUT type="text" name="-x" size="5" value="1,2"><BR>
                    </TD>
                </TR>
                <TR><TH align=right valign="top">Number of ions per peptide : </TH>
                    <TD bgcolor="$lightColor">
                        Min : <INPUT type="text" name="-o" size="4" value="3"> Max :<INPUT type="text" name="-n" value="20" size="4"><BR>
                    </TD>
                </TR>
                <TR><TH align=right valign="top">Files : </TH>
                    <TD bgcolor="$lightColor">
                        Windows SWATH file : <INPUT onmouseover="popup('File containing the swath ranges.')" onmouseout="popout()" type="file" name="swathfile" required><BR>
                        File with modifications delta mass* : <INPUT onmouseover="popup('File with the modifications not specified by default.')" onmouseout="popout()" type="file" name="-m" ><BR>
                        <SMALL>*File headers : modified-AA TPP-nomenclature Unimod-Accession ProteinPilot-nomenclature is_a_labeling composition-dictionary<BR>
                        An example : S S[167] 21 [Pho] False {'H':1,'O':3,'P':1}</SMALL><BR>
                        Labelling file : <INPUT onmouseover="popup('File containing the amino acid isotopic labelling mass shifts.')" onmouseout="popout()" type="file" name="-i"><BR>
                        Fasta file : <INPUT onmouseover="popup('Fasta file to relate peptides to their proteins.')" onmouseout="popout()" type="file" name="-f" accept=".fasta"><BR>
                    </TD>
                </TR>
                <TR><TH align=right valign="top">Other options : </TH>
                    <TD bgcolor="$lightColor">
                        <INPUT type="checkbox" name="other" value="-d">Remove duplicate masses from labelling<BR>
                        <INPUT type="checkbox" name="other" value="-e">Use theoretical mass<BR>
                        Time scale : <INPUT type="radio" name="-t" value="seconds" checked>seconds<INPUT type="radio" name="timescale" value="minutes">minutes<BR>
                        UIS order :  <INPUT onmouseover="popup('When using a switching modification, this determines the UIS order to be calculated.')" onmouseout="popout()" type="text" name="-y" value="2" size="3"><BR>
                        Maximum permissible error : <INPUT onmouseover="popup('Maximum error allowed at the annotation of a fragment ion.')" onmouseout="popout()" type="text" name="-p" value="0.05" size="4"><BR>
                        Allowed fragment mass modifications : <INPUT type="text" name="-g" size="10" onmouseover="popup('List of allowed fragment mass modifications. Useful for phosphorylations.')" onmouseout="popout()"><BR>
                    </TD>
                </TR>
                <TR><TH align="right" valign="top">Protein list :</TH>
                   <TD bgcolor="$lightColor"><INPUT type="file" name="protList"></INPUT><BR><SMALL>$help</SMALL></TD>
                </TR>|;
                my $sthModification=$dbh->prepare("SELECT SLM.ID_MODIFICATION, PSI_MS_NAME FROM MODIFICATION M, SWATH_LIB_MODIFICATION SLM WHERE ID_SWATH_LIB=$libID AND SLM.ID_MODIFICATION=M.ID_MODIFICATION");
                $sthModification->execute;
                if($sthModification->rows != 0){
                    print qq |
                <TR><TH align="right" valign="top">Peptide modifications to be excluded :</TH>
                    <TD bgcolor="$lightColor">
                    |;
                    while(my ($modifID,$modifName)=$sthModification->fetchrow_array){
                        print "<INPUT type=\"checkbox\" name=\"pepMod\" value=\"$modifID\">$modifName&nbsp;";
                    }                   
                    print qq |
                    </TD>
                </TR>|;
                }
        }
        
        print qq | <TR><TH colspan=2><br> |;
        print('<input type="submit" name="submit" value="Submit">') if($software eq 'TPP');
        print qq |
                   <!-- CLEAR button -->
                   &nbsp &nbsp &nbsp<INPUT type="reset" value="Clear" onclick="selectFormat('');" />
                   <!-- CANCEL button -->
                   &nbsp &nbsp &nbsp<INPUT type="button" value="Cancel" onclick="window.location='./listSwathLibraries.cgi'"></TH>
                </TR>
            </TABLE>
            </FORM>
        </DIV>
        <DIV id="divDescription" class="clDescriptionCont">
        <!--Empty div-->
        </DIV>
        <SCRIPT type="text/javascript">setPopup()</SCRIPT>
    </BODY>
    </CENTER>
    </HTML>
    |;
}

if($action){
    print qq
    |<HTML>
    <HEAD>
    <LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
    <STYLE type="text/css">
    </STYLE>
    </HEAD>
    
    <BODY  background="$promsPath{images}/bgProMS.gif">
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
        <DIV id="load">
            <BR><BR><IMG src='$promsPath{images}/engrenage.gif'>
    |;
    
    my $time = strftime("%Y%m%d%H%M%S", localtime);
    my $workDir = "$promsPath{tmp}/Swath/export_swath_lib_$libID\_$time";
    system("mkdir -p $workDir");
    
    #############################
    ###> fetching parameters <###
    #############################
    my ($swathWindows, $unimod, $w, $m, $fasta, $labelling);
    my $libraryID = param('ID');
    my $format = (param('-k')) ? param('-k') : 'openswath'; # peakview, openswath, spectronaut
    
    my @pepMod = (param('pepMod')) ? (param('pepMod')) : ();
    
    my $protList = param('protList');
    $protList = ($protList) ? $protList : '';
    if ($protList) {
        my $file = tmpFileName($protList);
        my $newDir = "$workDir/proteinList.txt";
        move($file, $newDir);
    }
    
    ### Upload TPP Window file
    $w = param('swathfile');
    if($w){
        uploadFile($workDir, $w);
    }
    
    ### Upload modifications file
    $m = param('-m');
    if($m) {
        uploadFile($workDir, $m);
    }
    
    ### Labelling file
    $labelling = param("-i");
    if($labelling) {
        uploadFile($workDir, $labelling);
    }
    
    ### Fasta file
    $fasta = param("-f");
    if($fasta) {
        uploadFile($workDir, $fasta);
    }
    
    ###> loading div
    print "<BR><BR><FONT color=\"red\"><DIV id=\"waitDIV\"></DIV></FONT><BR>";
    print "<BR><DIV id=\"loadingDIV\"></DIV><SPAN id=\"loadingSPAN\"></SPAN>";
    my $loadingDivID="document.getElementById('loadingDIV')";
    my $loadingSPAN="document.getElementById('loadingSPAN')";
    my $startTime=strftime("%s",localtime);
    my $processText="<B>Conversion to $format format ...</B>";
    
    my %libraryParams = (
        "ID"      => $libID,
        "name"    => $libraryName,
        "dir"     => "$promsPath{swath_lib}/SwLib_$libID",
        "version" => $libVersion,
    );
    
    DIAWorkflow::new($workDir, undef, \%promsPath, \%libraryParams);
    
    # Mass range (lmin/lmax), Ion series (s) and charge (x), Number of ions per peptide (min: o / max: n)
    my ($lMax, $lMin, $s, $x, $o, $n) = &promsMod::cleanParameters(param('-lmax'), param('-lmin'), param('-s'), param('-x'), param('-o'), param('-n'));
    
    # Maximum error (p), Allowed fragment mass modif (g), Fasta file (f), Time scale (t), UIS order (y), Labelling file (i == labelling)
    my ($p, $g, $t, $y) = &promsMod::cleanParameters(param('-p'), param('-g'), param('-t'), param('-y'));
    
    $labelling = "$workDir/$labelling" if($labelling);
    $w = "$workDir/$w" if($w);
    $fasta = "$workDir/$fasta" if($fasta);
    $m = "$workDir/$m" if($m);
    
    # Process other parameters (-d: remove duplicate masses / -e: use theorical mass)
    my $other = '';
    if (param('other')) {
        foreach my $value (param('other')) {
            $other .= " $value";
        }
        $other = &promsMod::cleanParameters($other);
    }
    
    my %exportParams = (
        "lmin" => $lMin, "lmax" => $lMax, "s" => $s, "x" => $x, "o" => $o, "n" => $n,
        "p" => $p, "g" => $g, "f" => $fasta, "t" => $t, "y" => $y, "i" => $labelling, "m" => $m, "other" => $other,
        "w" => $w, "format" => $format,
    );
    
    my ($finalFile, $paramFile) = DIAWorkflow::exportLibrary(\%exportParams, \@pepMod, 0, $protList, $processText, $loadingDivID, $loadingSPAN);
    
    ###> deleting temporary files
    system "rm $workDir/$w ;";
    system "rm $workDir/peakviewfile.txt" if -e "$workDir/peakviewfile.txt";
    system "rm $workDir/convertpeakview.sptxt" if  -e "$workDir/convertpeakview.sptxt";
    system "rm $workDir/convertopenswath.sptxt" if  -e "$workDir/convertopenswath.sptxt";
    system "rm $workDir/$libraryName\_peakview.tsv" if  -e "$workDir/$libraryName\_peakview.tsv";
    system "rm $workDir/proteinList.txt" if -e "$workDir/proteinList.txt";
    system "rm $m" if($m);
    system "rm $w" if($w);
    system "rm $fasta" if($fasta);
    system "rm $labelling" if($labelling);
    
    ###> Link to upload the final file : 
    print "<B>Done.</B><BR></DIV><BR>";
    print qq
    |   
        <SCRIPT LANGAGE="JavaScript">$loadingDivID.innerHTML="";$loadingSPAN.innerHTML="";</SCRIPT>
        <BR><BR><A class="button" href="$promsPath{tmp_html}/Swath/export_swath_lib_$libID/$finalFile" download><B>Download $format file</B></A>
        <BR><BR><A class="button" href="$promsPath{tmp_html}/Swath/export_swath_lib_$libID/$paramFile" download><B>Download param file</B></A><BR><BR><BR><BR><BR><BR>
        <INPUT type="button" class="buttonadd" value="Return to spectral libraries list." onclick="window.location='./listSwathLibraries.cgi'"></CENTER>
    </CENTER>
    </BODY>
    </HTML>
    |; 
}
print qq
|<DIV style="display:none" id="sptxt">
   <CENTER><BR><BR><A class="button" target="_blank" href="./exportSwathLibrary.cgi?ACT=download&ID=$libID&FORMAT=sptxt"><B>Download sptxt file</B></A></CENTER>
 </DIV>
 
 <DIV style="display:none" id="tsv">
   <CENTER><BR><BR><A class="button" target="_blank" href="./exportSwathLibrary.cgi?ACT=download&ID=$libID&FORMAT=tsv"><B>Download tsv file</B></A></CENTER>
 </DIV>
 |;

$dbh->disconnect; 


sub uploadFile {
    my ($destFolder, $fileName) = @_;
    
    if($fileName) {
        my $newFileDir = "$destFolder/".basename($fileName);
        
        if(not -e $newFileDir) {
            my $newFile = tmpFileName($fileName);
            move($newFile, $newFileDir);
        }
        
        return $newFileDir; 
    }
    
    return "";
}

sub downloadFile {
	my ($filePath, $fileName)=@_;
	print header(-type=>"application/octet-stream", -attachment=>$fileName);
	
	if(-e "$filePath/$fileName") {
        open(FILE,"$filePath/$fileName");
        while (<FILE>) {print $_;}
        close FILE;
    }
}

####>Revision history<####
# 1.2.2 [ENCHANCEMENT] Compatibility with Spectronaut spectral libraries (VS 22/10/2020)
# 1.2.1 [ENHANCEMENT] Handles spectronaut spectral library export format (VS 06/04/20)
# 1.2.0 [BUGFIX] Changed path to exported file: tmp_html directory (VS 16/12/19)
# 1.1.9 [ENHANCEMENT] Properly handles new DIA workflow for protein exportation (VS 18/11/19)
# 1.1.8 Update exportLibrary function to match the DIAWorkflow module (VS 07/01/19)
# 1.1.7 Minor modification on the exportLibrarySCRIPT call (VS 22/11/18)
# 1.1.6 Minor modif (MLP 13/12/17)
# 1.1.5 Create a function in promsDIA.pm to export library (MLP 06/12/17)
# 1.1.4 Minor modif on wait loop (MLP 31/10/17)
# 1.1.3 Update to allow the selection of peptide modifications to exlude (MLP 20/09/17)
# 1.1.2 Deleting N-ter acetylations in library file for openswath export (MLP 12/09/17)
# 1.1.1 Miror modif (MLP 27/07/17)
# 1.1.0 Changed the format of the exported file for OpenSwath (csv -> tsv) (MLP 11/05/17)
# 1.0.9 Minor corrections (MLP 17/03/17)
# 1.0.8 Minor modif on the form (MLP 28/02/17)
# 1.0.7 Add the selection of a protein list to export a part of the library (MLP 05/01/17)
# 1.0.6 Minor modif (MLP 29/11/16)
# 1.0.5 Bug correction for openswath conversion (MLP 28/09/16)
# 1.0.4 Minors modifications (MLP 08/07/16)
# 1.0.3 Change Swath_Lib path ($promsPath{data_html}/Swath_Lib to $promsPath{data_html}/swath_lib/) (MLP 07/07/16)
# 1.0.2 Deleting N-ter acetylations in sptxt file for peakview export (MLP 23/06/16)
# 1.0.1 Minor modifications to add pop up with definitions of export parameters (form) (MLP 02/05/16)
# 1.0.0 Create exportSwathLibrary.cgi to allow to export swath libraries for peakview or openswath (MLP 15/04/16)





