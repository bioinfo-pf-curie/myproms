#!/usr/local/bin/perl -w

################################################################################
# importSkyline.cgi         1.0.5                                              #
# Authors: V. Sabatet, M. Le Picard, V. Laigle (Institut Curie)                #
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

use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use strict;
use IO::Handle; 
use promsConfig;
use promsMod;
use promsImportDataManager;
use File::Copy;
use File::Basename;
use XML::Twig; 
use POSIX qw(strftime); # to get the time
use File::Path qw(rmtree); # remove_tree
use String::Util qw(trim);
use File::Spec::Functions qw(splitpath); # Core module

$|=1;       # buffer flush (allows instant printing)

#############################
####> Global parameters <####
#############################
my %promsPath = &promsConfig::getServerInfo;
my ($lightColor, $darkColor) = &promsConfig::getRowColors;
my $userID = $ENV{'REMOTE_USER'};

my $dbh = &promsConfig::dbConnect;
my $projectID = param('id_project');
my $branchID = param('ID');
my $submit = (param('submit')) ? param('submit') : "";

my ($item, $experimentID) = split(':', $branchID);

my ($fileError, $fileStat); # Log files

my (%peptideInfo, %peptideProt, %peptideModif, %peptidePosition, %fragmentsInfos, %analysisFileNames, %sampleFileNames, %peptideList, %dbRankProt, %assosProtPeptide);
my (%modifIDs, %searchModifs); # Store peptide modifications
my ($quantiAnnot, $isMS2, $hasMS1);

#######################
####>Starting HTML<####
#######################
print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);

## SHOULD DISPLAY IMPORT FORM
if ($submit eq "") {
    printImportForm();
    $dbh->disconnect;
    
## SHOULD PROCESS SUBMITED FORMS
} else {
    ### Print waiting screen
    print qq
    |<HTML>
        <HEAD>
            <LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
        </HEAD>
        <BODY style="background-image: url($promsPath{images}/bgProMS.gif); margin: auto; max-width: 900px;">
            <BR/><BR/><IMG src='$promsPath{images}/engrenage.gif'><BR/><BR/>
            <DIV id="waitDIV" style="color: red"></DIV><BR/>
            <B>Processing</B><BR/><BR/>
    |;

    ###############################
    ####> Fetching parameters <####
    ###############################
    my $file = &promsMod::cleanParameters(param('uploadfile'));
    my $dbID = &promsMod::cleanNumericalParameters(param('dbID'));
    
    my @dbID = ($dbID."&1");
    
    ### Creating new working directory
    my $startTime = strftime("%s",localtime);
    my $time = strftime("%Y%m%d%H%M%S",localtime);
    my $workDir = "$promsPath{data}/tmp/Swath/$time";
    mkdir $workDir unless -e $workDir;
    
    ### Declaring log files
    $fileStat="$workDir/status_$dbID.out";
    $fileError="$workDir/status_$dbID\_error.out";

    ### Moving uploaded file
    my $newFile = tmpFileName($file);
    my $newFileDir = "$workDir/$file";
    move($newFile, $newFileDir);
    my $fileName = fileparse($newFileDir, qr/\.[^.]*/);
    
    $dbh->disconnect;
    
    my $childProcess = fork;
    unless ($childProcess) {
        # Disconnecting from server
        #open STDOUT, ">$workDir/childprocess.out" or die "Can't open $workDir/childprocess.out: $!";
		open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
		open STDERR, ">$fileError" or die "Can't open $workDir/childprocess.err: $!";

        my ($fdr, $taxonomy, $instrument);
        my $dbh = &promsConfig::dbConnect;
        
        # Create new job to monitor
        $dbh->do("INSERT INTO JOB_HISTORY (ID_JOB, ID_USER, ID_PROJECT, ID_JOB_CLUSTER, TYPE, STATUS, FEATURES, SRC_PATH, LOG_PATH, ERROR_PATH, STARTED) VALUES('$time', '$userID', $projectID, 'L$$', 'Import [PRM]', 'Queued', 'SOFTWARE=Skyline', '$workDir', '$fileStat', '$fileError', NOW())");
        $dbh->commit;

        ###############################
        ###> Parsing uploaded file <###
        ###############################
        printToFile($fileStat, "Extracting information from Skyline file", 1);
        parseSkylineFile($dbh, $newFileDir, \%peptideInfo, \%peptideProt, \%peptideModif, \%peptidePosition, \%fragmentsInfos, \%analysisFileNames, \%sampleFileNames, \%peptideList, \%dbRankProt, \%assosProtPeptide, \%modifIDs);
        
        ################################
        ###> Initialize data import <###
        ################################
        my $importManager = promsImportDataManager->new($userID, $projectID, $experimentID, "prm", "Skyline_quantification", $quantiAnnot, $taxonomy, $instrument, $file, $fileStat, \@dbID);

        ###############################
        ### Import data in database ###
        ###############################
        $importManager->setFilesInfos(\%analysisFileNames, \%sampleFileNames);
        $importManager->importData(\%peptideInfo, \%peptideProt, \%peptideModif, \%peptidePosition, \%peptideList, \%fragmentsInfos, \%modifIDs, \%assosProtPeptide, \%dbRankProt);
        
        ############################
        ### Separating data file ###
        ############################
        my %analysisID = %{$importManager->getAnalysisID};
        my %peptidesID = %{$importManager->getPeptidesID};
        my $quantiID = $importManager->getQuantiID;
        my $hasGlobalPeptideXIC = 0;
        my ($analysisID, $quantiName, $quantiMethod, $fragmentQuantiID);  # For fragment quantification

        if ($isMS2) {  # Inserting fragments quantification data
            $quantiName = 'Skyline_Quantification (fragments)';
            $quantiMethod = $dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='TDA'");
                
            my $sthQuanti = $dbh->prepare("INSERT INTO QUANTIFICATION (ID_MODIFICATION,ID_QUANTIFICATION_METHOD,ID_DESIGN,NAME,FOCUS,QUANTIF_ANNOT,STATUS,UPDATE_DATE,UPDATE_USER) VALUES (NULL,?,NULL,?,?,?,?,NOW(),?)");
            $sthQuanti->execute($quantiMethod, $quantiName, 'peptide', $quantiAnnot, 1, $userID) or die $dbh->errstr;
            $sthQuanti->finish;
            $fragmentQuantiID=$dbh->last_insert_id(undef, undef, 'QUANTIFICATION', 'ID_QUANTIFICATION');
        }

        printToFile($fileStat, "Split each data files by analysis.", 1);
        foreach my $analysis (keys %analysisFileNames){
            $analysisID=$analysisID{$analysis};
            
            my $anaDir="$promsPath{peptide}/proj_$projectID/ana_$analysisID";
            mkdir $promsPath{'peptide'} unless -e $promsPath{'peptide'};
            mkdir "$promsPath{peptide}/proj_$projectID" unless -e "$promsPath{peptide}/proj_$projectID";
            mkdir $anaDir unless -e $anaDir;

            # Link analysis to fragments quantification
            if($isMS2 && $fragmentsInfos{$analysis}) {
                
                my $sthAnaQuanti = $dbh->prepare("INSERT INTO ANA_QUANTIFICATION (ID_QUANTIFICATION,ID_ANALYSIS,QUANTIF_FILE,IS_REFERENCE) VALUES (?,?,NULL,NULL)");
                $sthAnaQuanti->execute($fragmentQuantiID, $analysisID{$analysis});
                $sthAnaQuanti->finish;
                $dbh->commit;
                
                # Export fragment quantification data
                my $quantiPath = $promsPath{"quantification"}."/project_$projectID/quanti_$fragmentQuantiID";
                my $swathFile = "$quantiPath/swath_ana_".$analysisID{$analysis}.".txt";
                system "mkdir -p $quantiPath";
                open(FRAG_OUTFILE, ">", $swathFile) or die ("open: $!");
                print FRAG_OUTFILE "ID_PEPTIDE!Fragment MZ!Fragment Charge!Ion Type!Residue!Area!RT\n";
                
                foreach my $peptide (keys %{$fragmentsInfos{$analysis}}){
                    foreach my $rt (sort {$fragmentsInfos{$analysis}{$peptide}{$a} <=> $fragmentsInfos{$analysis}{$peptide}{$b}} keys %{$fragmentsInfos{$analysis}{$peptide}}){
                        foreach my $fragment (keys %{$fragmentsInfos{$analysis}{$peptide}{$rt}}){
                            my ($fragmentType,$fragmentCharge)=split(/_/,$fragment);
                            my ($ionType,$residue);
                            
                            $fragmentType =~ s/\s//;
                            if ($fragmentType =~ /([a-z]{1})(.*)/) {
                                $ionType=$1;
                                $residue=$2;
                            }
                            
                            my $fragmentRT = ($fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[0]) ? $fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[0] : '0.00';
                            my $fragmentMZ = ($fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[2]) ? $fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[2] : '0';
                            my $fragmentNeutralMass = $fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[3];
                            my $peptideID = $peptidesID{$analysis}{$peptide}{$fragmentRT}[0];
                            my $area = $fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[1];
                            
                            print FRAG_OUTFILE "$peptideID!$fragmentMZ!$fragmentCharge!$ionType!$residue!$area!$fragmentRT\n";
                        }
                    }
                }
                close FRAG_OUTFILE;
            }
            
            if($hasMS1) {
                # Export peptide quantification data
                my $quantiPath = $promsPath{"quantification"}."/project_$projectID/quanti_$quantiID";
                system "mkdir -p $quantiPath";
                
                open(PEP_OUTFILE, ">>", "$quantiPath/peptide_quantification.txt") or die ("open: $!");
                print PEP_OUTFILE "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n" if (-z "$quantiPath/peptide_quantification.txt");
                foreach my $peptide (keys %{$peptideList{$analysis}}) {
                    foreach my $rt (sort {$peptideList{$analysis}{$peptide}{$a} <=> $peptideList{$analysis}{$peptide}{$b}} keys %{$peptideList{$analysis}{$peptide}}){
                        my $pepArea = @{$peptideList{$analysis}{$peptide}{$rt}}[2];
                        my $peptideID = $peptidesID{$analysis}{$peptide}{$rt}[0];
                        
                        print PEP_OUTFILE "2\t$peptideID\t$pepArea\n";
                    }
                }
                close PEP_OUTFILE;
            }
        }
        
        if(!$hasMS1) {
            $dbh->do("DELETE FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$quantiID") || die $dbh->errstr;
            $dbh->commit;
            $dbh->do("DELETE FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantiID") || die $dbh->errstr;
            $dbh->commit;
        }
        
        printToFile($fileStat, "Import Ended", 1);
        
        $dbh->disconnect;
        exit;
    }
    
    print qq 
    |       <SCRIPT>
                var monitorJobsWin=window.open("$promsPath{cgi}/monitorJobs.cgi?filterType=Import [PRM]&filterDateNumber=1&filterDateType=DAY&filterStatus=Queued&filterStatus=Running",'monitorJobsWindow','width=1200,height=500,scrollbars=yes,resizable=yes');
                monitorJobsWin.focus();
                setTimeout(function() { parent.navFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=experiment:$experimentID&ACT=nav"; },7000);
            </SCRIPT>
        </BODY>
    </HTML>|;
}

=begin comment
sub parseSequenceModifications {
    my ($dbh, $modSeq ,$type, $peptideModif) = @_;
    my $errThreshold = 0.001; # Used for mass modifs
    my (%vmodsInUnimod, %modificationList, $varModName, $specificity); # Used for verbose modifs
    my %modifPatterns = (
        "mass"    => qr/(\w*)\[\+(\d+\.?\d*)\]/,
        "unimod"  => qr/(\w*)\(unimod:(\d+)\)/,
        "verbose" => qr/(\w*)\[([A-Za-z0-9 -]+)\s\((\w+)\)\]/,
    );
    
    # Fetching modifs mass and unimod ID informations <###
    my $sthMassMods = $dbh->prepare("SELECT ID_MODIFICATION, UNIMOD_ACC, MONO_MASS FROM MODIFICATION WHERE MONO_MASS IS NOT NULL ORDER BY UNIMOD_ACC");
    $sthMassMods->execute;
    while(my ($modID,$unimodID,$massMod) = $sthMassMods->fetchrow_array) {
        $modificationList{$modID}= {
            "mass" => $massMod,
            "unimod_id" => $unimodID,
        };
    }
    $sthMassMods->finish;
    
    # Parse current peptide for modifications
    while($modSeq =~ /$modifPatterns{$type}/g) {
        my $modifID;
        
        foreach my $modID (keys %modificationList) {
            if($type eq "mass") {
                next unless($modificationList{$modID}{"mass"});
                $modifID = $modID if ($2 >= $modificationList{$modID}{"mass"} - $errThreshold && $2 <= $modificationList{$modID}{"mass"} + $errThreshold);
            } elsif($type eq "unimod") {
                next unless($modificationList{$modID}{"unimod_id"});
                $modifID = $modID if ($modificationList{$modID}{"unimod_id"} == $2);
            } elsif($type eq "verbose") {
                $varModName = $2;
                $specificity = $3;
                $modifID = &promsMod::getModificationIDfromString($dbh, $varModName, $specificity);
            }
            last if $modifID;
        }
        
        unless($modifID) {
            printToFile($fileError, "ERROR: Couldn't find the $type peptide modification: $2", 0);
            exit;
        }
        
        my @aa=split("",$1);
        my $aaPrev=($aa[-1])? $aa[-1] : '=';
        my $pos=length($1);
        my $modMass = ($modificationList{$modifID}{"mass"}) ? $modificationList{$modifID}{"mass"} : "0";
        
        $modifIDs{$modifID}{'V'}{$aaPrev}=() if $modifID;
        $peptideModif->{$modifID}{$pos} = $modMass;
        
        if($type eq "mass") {
            my $patternContent = $2;
            $modSeq=~s/\[\+$patternContent\d*\]//;
        } elsif($type eq "unimod") {
            my $patternContent = "unimod:".$modificationList{$modifID}{"unimod_id"};
            $modSeq=~s/\($patternContent\)//;
        } elsif($type eq "verbose") {
            $varModName =~ s/-/\-/;
            $modSeq=~s/\[$varModName\s\($specificity\)\]//;
        }
    }
}
=end comment
=cut

sub printToFile {
    my ($filePath, $content, $append) = @_;
    my $openingMode = ($append) ? ">>" : ">";
    open(FILE, $openingMode, $filePath) || die "Error while opening $filePath\n";
	print FILE $content."\n";
	close FILE;
}

sub buildQuantiAnnot {
    my ($refAnnotSettings) = @_;

    my ($libCutOff, $libSearchFiles, $libIncludeAmbigous) = &promsMod::cleanParameters(param('libCutoff'), param('libSearchFiles'), param('libIncludeAmbigous'));
    $libIncludeAmbigous = ($libIncludeAmbigous && $libIncludeAmbigous eq "libIncludeAmbigous") ? 1 : 0;

    my ($fastaComment, $fastaProteinFilter, $fastaProteinFiltered) = &promsMod::cleanParameters(param('fastaComment'), param('fastaProteinFilter'), param('fastaProteinFiltered'));

    my $quantifAnnot = "LABEL=FREE::SOFTWARE=SKY;".$refAnnotSettings->{'softversion'};
    $quantifAnnot .= "::LIB_CUTOFF=".$libCutOff;
    $quantifAnnot .= ($libSearchFiles)? "::LIB_SEARCH_FILES=".$libSearchFiles : "::LIB_SEARCH_FILES=None";
    $quantifAnnot .= "::LIB_AMBIGOUS=".$libIncludeAmbigous;
    $quantifAnnot .= "::MS1_PRECURSOR_CHARGE=".$refAnnotSettings->{'ms1PrecursorCharge'};
    $quantifAnnot .= "::MS1_ISOTOPE_PEAKS=".$refAnnotSettings->{'ms1IsotopePeaks'};
    $quantifAnnot .= "::MS1_MASS_ANALYZER=".$refAnnotSettings->{'ms1MassAnalyzer'};
    $quantifAnnot .= "::MS1_PEAKS=".$refAnnotSettings->{'ms1Peaks'};
    $quantifAnnot .= "::MS1_POWER=".$refAnnotSettings->{'ms1Power'};
    $quantifAnnot .= "::MS1_POWER_CHARGE=".$refAnnotSettings->{'ms1PowerCharge'};
    
    if($refAnnotSettings->{'rtScans'} eq "UseFilteredScans") {
        $quantifAnnot .= "::RT_SCANS=FILTERED::RT_SCANS_FILTER=".$refAnnotSettings->{'rtScansFilter'};
    } else {
        $quantifAnnot .= "::RT_SCANS=ALL";
    }
    
    $quantifAnnot .= ($fastaComment)? "::FASTA_COMMENT=".$fastaComment : "::FASTA_COMMENT=None";
    $quantifAnnot .= "::FASTA_ENZYME=".$refAnnotSettings->{'fastaEnzyme'};
    $quantifAnnot .= "::FASTA_MAX_MC=".$refAnnotSettings->{'fastaMaxMC'};
    if ($fastaProteinFiltered && $fastaProteinFiltered eq "fastaProteinFiltered") {
        $quantifAnnot .= "::FASTA_FILTER=".$fastaProteinFilter;
    }

    $quantifAnnot .= ($refAnnotSettings->{'fixModifs'})? "::FIX_MODIFS=".$refAnnotSettings->{'fixModifs'} : "::FIX_MODIFS=None";
    $quantifAnnot .= ($refAnnotSettings->{'varModifs'})? "::VAR_MODIFS=".$refAnnotSettings->{'varModifs'} : "::VAR_MODIFS=None";
    
    return $quantifAnnot;
}

sub printImportForm {
    print qq
    |<!doctype html>
    <html lang="en-US">
        <head>
            <title>Import Skyline Data</title>
            <link rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
            <style>
                .popup { background-color: #FFFFFF; border: solid 3px #999999; padding: 5px; box-shadow: 10px 10px 20px #808080; position: absolute; display: none; }
                
                body { background-image: url($promsPath{images}/bgProMS.gif); margin: auto; max-width: 900px; }
                
                th { width: 145px; padding: 6px 5px 0 0; }
                
                input, select { display: inline-block; margin: 2px 0 1px 0; }
                
                td { line-height: 1.4em; padding: 2px 0 3px 5px; }
                
                td table td { padding: 0px; }
                
                input[type=checkbox], input[type=radio] { position: relative; top: 2px; margin-right: 2px; }
                
                input.number { text-align: right; }
            </style>
            <script type="text/javascript">
function cancelAction() {
    top.promsFrame.optionFrame.selectOption();
}
    |;
    &promsMod::popupInfo();
    &promsMod::browseDirectory_JavaScript if $promsPath{'shared'};
    print qq
    |        </script>
        </head>
        <body>
            <span class="title1" style="text-align:center; display:block; margin:25px 0 30px 0;">Select options to import Skyline data</span>
            <form method="POST" action="./importSkyline.cgi" name="parameters" enctype="multipart/form-data" style="margin: 0 auto;">
                <input name="ID" type="hidden" value="$branchID">
                <input name="id_project" type="hidden" value="$projectID">
                <table style="background-color:$darkColor">
                    <tr>
                        <th style="text-align:left; vertical-align:top; white-space:nowrap;">&nbsp;&nbsp;Type of quantification:&nbsp;&nbsp;</th>
                        <td style="background-color:$lightColor">
                            <select name="skylineQuantif" id="expType">
                                <option value="">-= Select =-</option>
                                <option value="ms1">MS1</option>
                                <option value="prm">PRM</option>
                            </select>
                        </td>
                    </tr>
                    <tr>
                        <th style="text-align:left; vertical-align:top; white-space:nowrap;">&nbsp;&nbsp;Skyline file (.sky):&nbsp;&nbsp;</th>
                        <td style="background-color:$lightColor">
                            <input type="file" name="uploadfile" accept=".sky" required />
                        </td>
                    </tr>
                    <tr>
                        <th style="text-align:left; vertical-align:top">&nbsp;&nbsp;Library:&nbsp;&nbsp;</th>
                        <td style="background-color:$lightColor">
                            cut-off score : 
                            <input type="text" class="number" value="95" size="3" value="0" name="libCutoff" />
                            <br/>
                            search files : 
                            <input type="text" size="30" name="libSearchFiles" /><br/>
                            <label>
                            <input type="checkbox" name="libIncludeAmbigous" value="libIncludeAmbigous" />
                            include ambigous matches
                            </label>
                        </td>
                    </tr>
                    <tr>
                        <th style="text-align:left; vertical-align:top; white-space:nowrap;">&nbsp;&nbsp;Fasta:&nbsp;&nbsp;</th>
                        <td style="background-color:$lightColor">
                            File : 
    |;

    my $sthDBList=$dbh->prepare("SELECT D.NAME,ID_DATABANK,FASTA_FILE,DT.NAME FROM DATABANK D, DATABANK_TYPE DT WHERE DT.ID_DBTYPE=D.ID_DBTYPE AND USE_STATUS='yes'") or die "Couldn't prepare statement: " . $dbh->errstr;
    $sthDBList->execute;
    if ($sthDBList->rows == 0) {
        print "No databank available.";
    } else {
        my %dbList;
        while (my($dbName,$dbID,$fastaFile,$dbType)=$sthDBList->fetchrow_array) {
            my ($dbSource,$version);
            if($fastaFile=~/:/) {
                my ($mascotServer,$dbankDir,$fileName)=split(':',$fastaFile);
                $dbSource=$mascotServer;
                $version=$dbankDir;
            } else {
                next if(! -e "$promsPath{banks}/db_$dbID/$fastaFile");
                $dbSource="Local";
            }
            $dbList{$dbSource}{$dbID}=$dbName;
            $dbList{$dbSource}{$dbID}.=" ($version)" if $version;
            $dbList{$dbSource}{$dbID}.=" [$dbType]";
        }
        $sthDBList->finish;
        $dbh->commit;

        print qq
        |                        <select name="dbID" id="dbID" required>
                                <option value="">-= Select databank =-</option>\n|;
        foreach my $dbSource (sort {$a cmp $b} keys %dbList){
            print "                                <optgroup label=\"$dbSource :\">\n";
            foreach my $dbID (sort {lc($dbList{$dbSource}{$a}) cmp lc($dbList{$dbSource}{$b})} keys %{$dbList{$dbSource}}) {
                print "                                    <option value=\"$dbID\">$dbList{$dbSource}{$dbID}</option>\n";
            }
            print "                                </optgroup>\n";
        }
        print "                            </select>\n";
    }
    
    print qq
    |                           <br/>
                            Comment on Database entries selection: 
                            <input type="text" value="" size="50" name="fastaComment" />
                            <br/>
                            <label><input type="checkbox" name="fastaProteinFiltered" value="fastaProteinFiltered">
                            Filter proteins having less than 
                            <input type="text" size="3" value="0" name="fastaProteinFilter" /> peptide(s)</label>
                        </td>
                    </tr>
                    <tr>
                        <th colspan=2>
                            <!-- SUBMIT button -->
                            <input type="submit" name="submit" value="Submit">&nbsp;&nbsp;&nbsp;&nbsp;
                            
                            <!-- CLEAR button -->
                            <input type="reset" value="Clear" />&nbsp;&nbsp;&nbsp;&nbsp;
                            
                            <!-- CANCEL button -->
                            <input type="button" value="Cancel" onclick="cancelAction();" >
                        </th>
                    </tr>
                </table>
            </form>

            <!-- Required for description popups -->
            <div id="divDescription" class="clDescriptionCont"></div>
            <script>setPopup();</script>
        </body>
    </html>
    |;
}


sub parseSkylineFile {
    my ($dbh, $skylineFile, $refPeptideInfo, $refPeptideProt, $refPeptideModif, $refPeptidePosition, $refFragmentsInfos, $refAnalysisFileNames, $refSampleFileNames, $refPeptideList, $refDbRankProt, $refAssosProtPeptide, $refModifIDs) = @_;

    # Parameters from form
    my $skylineQuantif = param('skylineQuantif') || 'ms1';
    if (!$skylineFile) {
        die "Unable to find or to open the Skyline file (.sky)"
    }
    if ($skylineQuantif eq 'prm') {
        $isMS2 = 1;
    }

    # Settings variables
    my ($softwareVersion);
    my ($libCutOff, $libSearchFiles, $libIncludeAmbigous);
    my ($ms1PrecursorCharge, $ms1IsotopePeaks, $ms1MassAnalyzer, $ms1Peaks, $ms1Power, $ms1PowerCharge);
    my ($rtScans, $rtScansFilter);
    my ($fastaEnzyme, $fastaMaxMC, $fastaComment, $fastaProteinFilter, $fastaProteinFiltered);
    my ($fixModifs, $varModifs);
    
    # Sub reference for XML parsing with XML::Twig
    my $settingsSub = sub {
        my ($innerTwig, $settings) = @_;  # twig handlers params
        my $peptideSettings     = $settings->first_child('peptide_settings');
        my $transitionSettings  = $settings->first_child('transition_settings');
        my $measuredResults     = $settings->first_child('measured_results');

        $fastaEnzyme = $peptideSettings->first_child('enzyme')->{'att'}->{'name'};
        $fastaMaxMC = $peptideSettings->first_child('digest_settings')->{'att'}->{'max_missed_cleavages'};

        # Possible modifications
        my (@allModifs, @fixedModifs, @variableModifs);
        push @allModifs, $peptideSettings->get_xpath('peptide_modifications/static_modifications/static_modification');
        push @allModifs, $peptideSettings->get_xpath('peptide_modifications/heavy_modifications/heavy_modification');
        
        foreach my $modif (@allModifs) {
            my $modifAttrs = $modif->atts; 
            
            my $modifFullName = $modifAttrs->{'name'};
            $modifFullName =~ s/\s/_/g;
            my $modifName = (split('_', $modifFullName))[0];
            my $modifType;
            if ($modifAttrs->{'variable'} && $modifAttrs->{'variable'} eq "true") {
                $modifType = 'V';
                push @variableModifs, $modifFullName;
            } else {
                $modifType = 'F';
                push @fixedModifs, $modifFullName;
            }
            my $specificity;
            if ($modifAttrs->{'aminoacid'} && $modifAttrs->{'terminus'}) {
                $specificity = $modifAttrs->{'terminus'};
                $specificity =~ s/N/=/;
                $specificity =~ s/C/\*/;
                $specificity .= ', '.$modifAttrs->{'aminoacid'};
            } elsif ($modifAttrs->{'aminoacid'}) {
                $specificity = $modifAttrs->{'aminoacid'};
            } elsif ($modifAttrs->{'terminus'}) {
                $specificity = $modifAttrs->{'terminus'};
                $specificity =~ s/N/=/;
                $specificity =~ s/C/\*/;
            }
            my $modifID = &promsMod::getModificationIDfromString($dbh, $modifName, $specificity);
            unless ($modifID) {
                printToFile($fileError, "ERROR: Couldn't retrieve any ID for the modification: $modifFullName", 0);
                exit;
            }

            my $modifMass = $dbh->selectrow_array("SELECT MONO_MASS FROM MODIFICATION WHERE ID_MODIFICATION=$modifID");
            
            $searchModifs{$modifFullName} = {
                'name'          => $modifName,
                'ID'            => $modifID,
                'monoMass'      => $modifMass,
                'unimod_id'     => $modifAttrs->{'unimod_id'},
                'aminoacid'     => $modifAttrs->{'aminoacid'},
                'terminus'      => $modifAttrs->{'terminus'},
                'type'          => $modifType,
                'specificity'   => $specificity
            }
        }

        $fixModifs = join(',', sort {$a cmp $b} @fixedModifs);
        $fixModifs =~ s/_/ /g;
        $varModifs = join(',', sort {$a cmp $b} @variableModifs);
        $varModifs =~ s/_/ /g;

        # Instrument/transition settings
        $ms1PrecursorCharge = $transitionSettings->first_child('transition_filter')->{'att'}->{'precursor_charges'};
        
        my $tranScanAttrs   = $transitionSettings->first_child('transition_full_scan')->atts;
        $ms1IsotopePeaks    = $tranScanAttrs->{'precursor_isotopes'};
        $ms1MassAnalyzer    = $tranScanAttrs->{'precursor_mass_analyzer'};
        $ms1MassAnalyzer    = ucfirst($ms1MassAnalyzer);
        $ms1Peaks           = $tranScanAttrs->{'precursor_isotope_filter'};
        $ms1Power           = $tranScanAttrs->{'precursor_res'};
        $ms1PowerCharge     = $tranScanAttrs->{'precursor_res_mz'};
        $rtScans            = ($tranScanAttrs->{'retention_time_filter_type'} && $tranScanAttrs->{'retention_time_filter_type'} eq "ms2_ids")? "UseFilteredScans" : "UseAllScans";
        $rtScansFilter      = $tranScanAttrs->{'retention_time_filter_length'};

        # Info from replicates
        my @replicateElmts = $measuredResults->children('replicate');
        foreach my $replicateElmt (@replicateElmts) {
            my $anaName = $replicateElmt->{'att'}->{'name'};
            
            # Get basename of raw file
            my $filePath = $replicateElmt->first_child('sample_file')->{'att'}->{'file_path'};
            my $ostype = ($filePath =~ /\\/) ? "MSDOS" : "UNIX";
            fileparse_set_fstype($ostype);
            my ($fileName, $path, $suffix) = fileparse($filePath, qr{\.\w+$});            
            my $anaFile = $fileName.$suffix;

            $refAnalysisFileNames->{$anaName} = $anaFile;
            $refSampleFileNames->{$anaName} = $anaName;
        }
        $innerTwig->purge;
    };
    
    # Sub reference for XML parsing with XML::Twig
    my $proteinSub = sub {
        my ($innerTwig, $protein) = @_;  # twig handlers params

        my $prot = $protein->{'att'}->{'name'};
        if ($prot) {
            $refDbRankProt->{$prot}=1;
        } else {
            die "Your file is missing protein name(s)";
        }

        my @peptides = $protein->children('peptide');  # Get all the peptides elements from the protein
        foreach my $pep (@peptides) {
            my $pepAttrs = $pep->atts;

            ### Peptide/Fragment content ###
            my $beg         = $pepAttrs->{'start'};
            my $end         = $pepAttrs->{'end'};
            my $begFlank    = $pepAttrs->{'prev_aa'};
            my $endFlank    = $pepAttrs->{'next_aa'};
            my $missCleav   = $pepAttrs->{'num_missed_cleavages'};
            my $pepSeq      = $pepAttrs->{'sequence'};
            my $modifSeq    = $pepAttrs->{'modified_sequence'};  # mass of modif in seq

            $beg            = (!$beg || $beg eq '#N/A') ? '' : $beg+1;
            $end            = (!$end || $end eq '#N/A') ? '' : $end;
            $begFlank       = (!$begFlank || $begFlank eq 'X') ? '' : $begFlank;
            $endFlank       = (!$endFlank || $endFlank eq 'X') ? '' : $endFlank;

            my %peptideModifs;
            my @allModifs;
            push @allModifs, $pep->get_xpath('variable_modifications/variable_modification');
            push @allModifs, $pep->get_xpath('implicit_modifications/implicit_static_modifications/implicit_modification');
            push @allModifs, $pep->get_xpath('implicit_modifications/implicit_heavy_modifications/implicit_modification');
            
            foreach my $modif (@allModifs) {
                my $modifPos        = $modif->{'att'}->{'index_aa'};
                my $modifFullName   = $modif->{'att'}->{'modification_name'};
                $modifFullName      =~ s/\s/_/g;
                
                my $modifID = $searchModifs{$modifFullName}{'ID'};
                unless ($modifID) {
                    printToFile($fileError, "ERROR: Couldn't retrieve any ID for the modification: $modifFullName", 0);
                    exit;
                }
                my $modifMass = $searchModifs{$modifFullName}{'monoMass'};

                my $aa;
                if ($modifPos == 0 && $searchModifs{$modifFullName}{'specificity'} =~ /[\-\=]/) {
                    $aa = '=';
                } elsif ($modifPos == (length($pepSeq)-1) && $searchModifs{$modifFullName}{'specificity'} =~ /[\*\+]/) {
                    $aa = '*';
                } else {
                    $aa = (split("", $pepSeq))[$modifPos];
                }
                $modifPos++;  # Switch position to 1-based for DB and display
                
                $modifIDs{$modifID}{$searchModifs{$modifFullName}{'type'}}{$aa}=();
                if ($searchModifs{$modifFullName}{'type'} eq 'V') {  # Only variable modifs stored in DB
                    $peptideModifs{$modifID}{$modifPos} = $modifMass;
                }
            }

            my %peptideRTs;
            my @peptideResults = $pep->get_xpath('peptide_results/peptide_result');
            foreach my $pepResult (@peptideResults) {
                my $rep = $pepResult->{'att'}->{'replicate'};
                my $rt  = $pepResult->{'att'}->{'retention_time'};
                $rt  = (!$rt || $rt eq '#N/A') ? '0.00' : sprintf("%.2f", $rt);
                $peptideRTs{$rep} = $rt;
            }

            my @precursors = $pep->children('precursor');
            foreach my $precursor (@precursors) {   
                my $precursorAttrs = $precursor->atts;  
            
                my $pepMZ     = $precursorAttrs->{'precursor_mz'};
                my $pepCharge = $precursorAttrs->{'charge'};
            
                my $peptide = $modifSeq.'_'.$pepCharge;
                $refPeptideInfo->{$peptide}[1] = $missCleav;

                my @precursorPeaks = $precursor->get_xpath('precursor_results/precursor_peak');
                foreach my $precursorPeak (@precursorPeaks) {
                    my $peakAttrs = $precursorPeak->atts;
                
                    my $anaName      = $peakAttrs->{'replicate'};
                    my $pepRT        = $peakAttrs->{'retention_time'};
                    my $totalArea    = $peakAttrs->{'area'};
                    my $library_dotp = $peakAttrs->{'library_dotp'};
                    my $isotope_dotp = $peakAttrs->{'isotope_dotp'};

                    $pepRT = (!$pepRT || $pepRT eq '#N/A') ? '0.00' : sprintf("%.2f", $pepRT);
                    next if (!$totalArea || $totalArea eq '#N/A');
                    if ($totalArea && !$hasMS1) {
                        $hasMS1 = 1;
                    }
                    $totalArea = sprintf("%.0f", $totalArea);
                    my $pepScore;
                    if ($isMS2) {
                        $pepScore = ($library_dotp)? sprintf("%.6f", $library_dotp) : ''; 
                    } else {
                        $pepScore = ($isotope_dotp)? sprintf("%.6f", $isotope_dotp) : '';
                    }

                    $refPeptideProt->{$anaName}{$peptide}{$prot}=1;
                    #parseSequenceModifications($dbh, $modifSeq, "mass", \%{$refPeptideModif->{$anaName}{$peptide}});
                    $refPeptideModif->{$anaName}{$peptide} = \%peptideModifs;
                    $refPeptidePosition->{$anaName}{$prot}{$peptide}{"$beg\_$begFlank"}="$end\_$endFlank" if ($beg && $begFlank && $end && $endFlank);
                    @{$refAssosProtPeptide->{$anaName}{$prot}{$peptide}}=($pepRT);
                    if (!$refPeptideList->{$anaName}{$peptide}{$pepRT}) { # Current line correspond to a peptide
                        @{$refPeptideList->{$anaName}{$peptide}{$pepRT}}=($pepMZ, $pepScore, $totalArea);
                    }
                }

                if ($isMS2) {
                    my @transitions = $precursor->children('transition');

                    foreach my $transition (@transitions) {
                        next if ($transition->{'att'}->{'fragment_type'} =~ /precursor/);
                    
                        my $transAttrs = $transition->atts;
                        my $fragmentType        = $transAttrs->{'fragment_type'};  # y, b, etc.
                        my $fragmentOrdinal     = $transAttrs->{'fragment_ordinal'};  # 1,2,...,8,9, etc
                        my $fragmentNeutralMass = $transAttrs->{'calc_neutral_mass'};
                        my $fragmentCharge      = $transAttrs->{'product_charge'};
                    
                        my $fragmentMZ          = $transition->first_child('product_mz')->text;
                    
                        my $fragmentName        = $fragmentType.$fragmentOrdinal;  # y1/y8/b1/b13/etc. 
                        my $fragment            = $fragmentName.'_'.$fragmentCharge;

                        my @transitionPeaks = $transition->get_xpath('transition_results/transition_peak');
                        foreach my $transPeak (@transitionPeaks) {
                            my $anaName      = $transPeak->{'att'}->{'replicate'};
                            #my $fragmentRT   = $transPeak->{'att'}->{'retention_time'};
                            my $fragmentArea = $transPeak->{'att'}->{'area'};
                            next if ($peptideRTs{$anaName} eq '0.00' || !$fragmentArea || $fragmentArea eq '#N/A');

                            $fragmentArea = sprintf("%.0f", $fragmentArea);
                            @{$refFragmentsInfos->{$anaName}{$peptide}{$peptideRTs{$anaName}}{$fragment}} = ($peptideRTs{$anaName},$fragmentArea,$fragmentMZ,$fragmentNeutralMass,$fragmentName);
                        }
                    }
                }
            }
        }
        $innerTwig->purge;
    };

    my $twig = XML::Twig->new(
        twig_handlers => {
            settings_summary    => $settingsSub,
            protein             => $proteinSub
        }
    );
    $twig->parsefile($skylineFile);
    
    my $root = $twig->root;
    $softwareVersion = $root->{'att'}{'software_version'};
    $softwareVersion =~ s/^.*\s(.+)$/$1/;
    
    my %annotSettings = (
        'softversion'           => $softwareVersion,
        #'libCutoff'             => undef,  # parameters still in form, wasn't able to get it from .sky file
        #'libSearchFiles'        => undef,
        #'libIncludeAmbigous'    => undef,
        'ms1PrecursorCharge'    => $ms1PrecursorCharge,
        'ms1IsotopePeaks'       => $ms1IsotopePeaks,
        'ms1MassAnalyzer'       => $ms1MassAnalyzer,
        'ms1Peaks'              => $ms1Peaks,
        'ms1Power'              => $ms1Power,
        'ms1PowerCharge'        => $ms1PowerCharge,
        'rtScans'               => $rtScans,
        'rtScansFilter'         => $rtScansFilter,
        'fastaEnzyme'           => $fastaEnzyme,
        'fastaMaxMC'            => $fastaMaxMC,
        #'fastaComment'          => undef,
        #'fastaProteinFilter'    => undef,
        #'fastaProteinFiltered'  => undef,
        'fixModifs'             => $fixModifs,
        'varModifs'             => $varModifs
    );

    $quantiAnnot = buildQuantiAnnot(\%annotSettings);
}

####>Revision history<####
# 1.0.5 [BUGFIX] Minor correction on quantifAnnot (VL 29/11/19)
# 1.0.4 [CHANGES] Use new job monitoring window opening parameters (VS 18/11/19)
# 1.0.3 Fits the new job monitoring system behavior (VS 09/10/19)
# 1.0.2 Create only one quantification for all analyses for MS2/fragments quantifications (VL 12/09/19)
# 1.0.1 Add handling of skyline scores for peptides (VL 11/09/19)
# 1.0.0 Created from parts of importTDA.cgi to replace it (VL 28/08/2019)
