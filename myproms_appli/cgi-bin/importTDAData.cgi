#!/usr/local/bin/perl -w

################################################################################
# importTDAData.cgi         1.0.3                                              #
# Authors: V. Sabatet, M. Le Picard (Institut Curie)                           #
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
use XML::Simple;
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
my $fileFormat = ".csv,.tsv";

my ($fileError, $fileInfo, $fileStat); # Log files
my (%modifIDs); # Store peptide modifications


#######################
####>Starting HTML<####
#######################
print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);

## SHOULD DISPLAY IMPORT FORM
if ($submit eq "") {
    printImportForm();
    
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
    my $quantiAnnot = buildQuantiAnnot();
    my @dbID = ($dbID."&1");
    
    ### Creating new working directory
    my $startTime = strftime("%s",localtime);
    my $time = strftime("%Y%m%d%H%M%S",localtime);
    my $workDir = "$promsPath{data}/tmp/Swath/$time";
    mkdir $workDir unless -e $workDir;
    
    ### Declaring log files
    $fileInfo = "$workDir/info_$dbID.out";
    $fileStat="$workDir/status_$dbID.out";
    $fileError="$workDir/status_$dbID\_error.out";

    printToFile($fileInfo, "User=$userID\nSoftware=Skyline", 0);
    printToFile($fileStat, "Started ".strftime("%H:%M:%S %d/%m/%Y",localtime), 0);

    ### Moving uploaded file
    my $newFile = tmpFileName($file);
    my $newFileDir = "$workDir/$file";
    move($newFile, $newFileDir);
    my $fileName = fileparse($newFileDir, qr/\.[^.]*/);

    my $childProcess=fork;
    unless ($childProcess) {
        # Disconnecting from server
        open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
		open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
		open STDERR, '>/dev/null' or die "Can't open /dev/null: $!";

        my ($fdr, $taxonomy, $instrument, $isMS2);
        $dbh = &promsConfig::dbConnect;
        
        
        ################################
        ###> Initialize data import <###
        ################################
        my $importManager = promsImportDataManager->new($userID, $projectID, $experimentID, "prm", "Skyline_quantification", $quantiAnnot, $taxonomy, $instrument, $file, $fileStat, \@dbID);

        
        ###############################
        ###> Parsing uploaded file <###
        ###############################
        my (%peptideInfo, %peptideProt, %peptideModif, %peptidePosition, %fragmentsInfos, %analysisFileNames, %sampleFileNames, %peptideList, %dbRankProt, %assosProtPeptide);
        my (%colName2Index, $colSeqName, $colRTName, $totFragmentRT, $nbFragmentPeptide);
        my %refColNames = (
            default => ['Protein Name', 'Replicate Name','File Name','Precursor Mz','Precursor Charge','Missed Cleavages','Begin Pos','End Pos','Previous Aa','Next Aa'],
            seq     => ['Peptide Modified Sequence Monoisotopic Masses', 'Peptide Modified Sequence Unimod Ids', 'Peptide Modified Sequence Average Masses'],
            RT      => ['Retention Time','Peptide Retention Time'],
            MS1     => ['Total Area MS1'],
            MS2     => ['Product Mz','Product Neutral Mass','Product Charge','Fragment Ion','Area','Total Area'],
        );

        printToFile($fileStat, "Scanning prm file.", 1);
        
        open(IN,"<$workDir/$file");
        while(<IN>){
            $_=~s/\s*\Z//;
            my @lineData = split(/[;,\t]/, $_);
            
            if($.==1) { # Check columns header
                my $errorMsg;
                $errorMsg = $colSeqName = $colRTName = "";
                
                foreach my $i (0 .. $#lineData) {
                    $colName2Index{$lineData[$i]} = $i;
                }

                # Check default columns names
                foreach my $colName (@{$refColNames{"default"}}) {
                    unless (exists $colName2Index{$colName}){
                        $errorMsg = "Required column '$colName' is missing from input file !";
                        last;
                    }
                }
                
                # Check for MS2 colums
                foreach my $colName (@{$refColNames{"MS2"}}) {
                    if(exists $colName2Index{$colName}) {
                        $isMS2 = 1;
                    } elsif($isMS2) {
                        $errorMsg = "Wrong columns names were provided for MS2 ! Following columns are required: ".join('<BR/>', @{$refColNames{"MS2"}});
                        last;
                    }
                }
                
                if(!$isMS2) {
                    # Check for MS1 colums
                    foreach my $colName (@{$refColNames{"MS1"}}) {
                        unless(exists $colName2Index{$colName}) {
                            $errorMsg = "Wrong or missing columns names were provided for MS1 ! Following columns are required: ".join('<BR/>', @{$refColNames{"MS1"}});
                            last;
                        }
                    }
                }
                
                # Check for peptide sequence and rt column
                foreach my $colName (@lineData) {
                    last if $colSeqName ne '' && $colRTName ne '';
                    if(grep(/$colName/, @{$refColNames{"seq"}}) && $colSeqName eq '') {
                        $colSeqName = $colName ;
                    } elsif(grep(/$colName/, @{$refColNames{"RT"}}) && $colRTName eq '') {
                        $colRTName = $colName;
                    }
                }
                
                if($colSeqName eq '') {
                    $errorMsg = "No sequence column in file ! The sequence columns names allowed are : ".join('<BR/>', @{$refColNames{"seq"}});
                } elsif($colRTName eq '') {
                    $errorMsg = "No retention time column in file ! The sequence columns names allowed are : ".join('<BR/>', @{$refColNames{"RT"}});
                }
                
                if($errorMsg) {
                    printToFile($fileError, $errorMsg, 0);
                    exit;
                }
            } else { # Peptide/Fragment content
                my $prot=$lineData[$colName2Index{'Protein Name'}];
                my $anaName     = $lineData[$colName2Index{'Replicate Name'}];
                my $anaFile     = $lineData[$colName2Index{'File Name'}];
                my $pepMZ       = $lineData[$colName2Index{'Precursor Mz'}];
                my $pepCharge   = $lineData[$colName2Index{'Precursor Charge'}];
                my $modifSeq    = $lineData[$colName2Index{$colSeqName}];
                my $missCleav   = $lineData[$colName2Index{'Missed Cleavages'}];
                my $pepRT       = $lineData[$colName2Index{$colRTName}];
                my $area        = ($isMS2) ? $lineData[$colName2Index{'Total Area'}] : $lineData[$colName2Index{'Total Area MS1'}];
                my $beg         = $lineData[$colName2Index{'Begin Pos'}];
                my $end         = $lineData[$colName2Index{'End Pos'}];
                my $begFlank    = $lineData[$colName2Index{'Previous Aa'}];
                my $endFlank    = $lineData[$colName2Index{'Next Aa'}];
                
                next if(!$prot || $area eq '#N/A');
                
                $pepRT          = '0.00' if ($pepRT eq '#N/A');
                $beg            = ($beg eq '#N/A') ? '' : $beg+1;
                $end            = ($end eq '#N/A') ? '' : $end+1;
                $begFlank       = ($begFlank eq 'X') ? '' : $begFlank;
                $endFlank       = ($endFlank eq 'X') ? '' : $endFlank;

                $analysisFileNames{$anaName} = $anaFile;
                $sampleFileNames{$anaName} = $anaName;

                my $peptide = $modifSeq.'_'.$pepCharge;
                
                ## peptide position
                $peptidePosition{$anaName}{$prot}{$peptide}{"$beg\_$begFlank"}="$end\_$endFlank" if ($beg && $begFlank && $end && $endFlank);
                $dbRankProt{$prot} = 1;

                # Fetching peptide modifications
                if($modifSeq=~/\w*\[\+\d+\.?\d*\]/) {
                    ## Checking for mass modifications : [+47] or [+47.00109]
                    parseSequenceModifications($modifSeq, "mass", \%{$peptideModif{$anaName}{$peptide}});
                } elsif ($modifSeq=~/\w*\(unimod:\d+\)/ || $modifSeq =~ /(\w*)\[([A-Za-z0-9 -]+)\s\((\w+)\)\]/) {
                    ## Checking unimod modifications : (unimod:ID)
                    parseSequenceModifications($modifSeq, "unimod", \%{$peptideModif{$anaName}{$peptide}});
                    ### With unimod export option, some PTMs match the following syntax: [Methyl-Propionyl (K)]
                    parseSequenceModifications($modifSeq, "verbose", \%{$peptideModif{$anaName}{$peptide}});
                } elsif ($modifSeq !~ /\w+/g) {
                    printToFile($fileError, "ERROR: Couldn't find the modification type in peptide: $modifSeq", 0);
                    exit;
                }

                ### fetching peptide informations
                $peptideInfo{$peptide}[1]=$missCleav;

                ### fetching fragment informations
                if($isMS2) {
                    my $fragmentMZ          = $lineData[$colName2Index{'Product Mz'}];
                    my $fragmentNeutralMass = $lineData[$colName2Index{'Product Neutral Mass'}];
                    my $fragmentCharge      = $lineData[$colName2Index{'Product Charge'}];
                    my $fragmentType        = $lineData[$colName2Index{'Fragment Ion'}];
                    my $fragment            = $fragmentType.'_'.$fragmentCharge;
                    
                    # current line correspond to a fragment
                    if($fragmentType !~ /precursor/) {
                        $area = $lineData[$colName2Index{"Area"}];
                        next if ($pepRT eq '0.00' || $area eq '#N/A');
                        @{$fragmentsInfos{$anaName}{$peptide}{$pepRT}{$fragment}} = ($pepRT,$area,$fragmentMZ,$fragmentNeutralMass,$fragmentType);
                    } elsif(!$peptideList{$anaName}{$peptide}{"$pepRT"}) { # Current line correspond to a peptide
                        @{$peptideList{$anaName}{$peptide}{"0.00"}}=($pepMZ, '', $area);
                    }
                } else {
                    @{$peptideList{$anaName}{$peptide}{$pepRT}} = ($pepMZ,'',$area);
                }
                
                ### List of protein
                $peptideProt{$anaName}{$peptide}{$prot}=1;
                @{$assosProtPeptide{$anaName}{$prot}{$peptide}}=($pepRT);
            }
        }
        close IN;
        
        
        ###############################
        ### Import data in database ###
        ###############################
        $importManager->setFilesInfos(\%analysisFileNames, \%sampleFileNames);
        $importManager->importData(\%peptideInfo, \%peptideProt, \%peptideModif, \%peptidePosition, \%peptideList, \%fragmentsInfos, \%modifIDs, \%assosProtPeptide, \%dbRankProt);
        
        
        ############################
        ### Separating data file ###
        ############################
        my ($analysisID, $warning);
        my %analysisID = %{$importManager->getAnalysisID};
        my %peptidesID = %{$importManager->getPeptidesID};
        my $quantiID = $importManager->getQuantiID;
        
        printToFile($fileStat, "Split each data files by analysis.", 1);
        foreach my $analysis (keys %analysisFileNames){
            $analysisID=$analysisID{$analysis};
            
            my $anaDir="$promsPath{peptide}/proj_$projectID/ana_$analysisID";
            mkdir $promsPath{'peptide'} unless -e $promsPath{'peptide'};
            mkdir "$promsPath{peptide}/proj_$projectID" unless -e "$promsPath{peptide}/proj_$projectID";
            mkdir $anaDir unless -e $anaDir;

            # Create a quantification for fragments quantification
            if($isMS2 && $fragmentsInfos{$analysis}) {
                my $quantiName = 'Skyline_Quantification (fragments)';
                my $quantiMethod = $dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='TDA'");
                
                # Inserting fragments quantification data
                my $sthQuanti = $dbh->prepare("INSERT INTO QUANTIFICATION (ID_MODIFICATION,ID_QUANTIFICATION_METHOD,ID_DESIGN,NAME,FOCUS,QUANTIF_ANNOT,STATUS,UPDATE_DATE,UPDATE_USER) VALUES (NULL,?,NULL,?,?,?,?,NOW(),?)");
                $sthQuanti->execute($quantiMethod, $quantiName, 'peptide', $quantiAnnot, 1, $userID) or die $dbh->errstr;
                $sthQuanti->finish;
                my $fragmentQuantiID=$dbh->last_insert_id(undef, undef, 'QUANTIFICATION', 'ID_QUANTIFICATION');
                
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
                            my $peptideID = $peptidesID{$analysis}{$peptide}[0];
                            my $area = $fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[1];
                            print FRAG_OUTFILE "$peptideID!$fragmentMZ!$fragmentCharge!$ionType!$residue!$area!$fragmentRT\n";
                        }
                    }
                }
                close FRAG_OUTFILE;
            }
            
            # Export peptide quantification data
            my $quantiPath = $promsPath{"quantification"}."/project_$projectID/quanti_$quantiID";
            system "mkdir -p $quantiPath";
            open(PEP_OUTFILE, ">>", "$quantiPath/peptide_quantification.txt") or die ("open: $!");
            print PEP_OUTFILE "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n" if (-z "$quantiPath/peptide_quantification.txt");
            foreach my $peptide (keys %{$peptideList{$analysis}}) {
                my $peptideID = $peptidesID{$analysis}{$peptide}[0];
                foreach my $RT (sort {$peptideList{$analysis}{$peptide}{$a} <=> $peptideList{$analysis}{$peptide}{$b}} keys %{$peptideList{$analysis}{$peptide}}){
                    my $pepArea = @{$peptideList{$analysis}{$peptide}{"$RT"}}[2];
                    print PEP_OUTFILE "2\t$peptideID\t$pepArea\n";
                    last;
                }
            }
            close PEP_OUTFILE;
        }
        
        printToFile($fileStat, "Done.", 1);
        printToFile($fileStat, "WARNING : $warning", 1) if $warning;
        system "rm $fileStat" if (-e $fileStat && !$warning);
    }
    
    print qq 
    |       <SCRIPT>
                monitorWindow=window.open("$promsPath{cgi}/monitorDIAProcess.cgi?ACT=import", 'Monitoring import of TDA data', 'width=1000,height=500,scrollbars=yes,resizable=yes');
                top.promsFrame.selectedAction='summary';
                setTimeout(function() { parent.navFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=experiment:$experimentID&ACT=nav"; },7000);
            </SCRIPT>
        </BODY>
    </HTML>|;
}

$dbh->disconnect;


sub parseSequenceModifications {
    my ($modSeq ,$type, $peptideModif) = @_;
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
                $modifID = &promsMod::getModificationIDfromString($dbh, $varModName, $specificity,\%vmodsInUnimod);
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

sub printToFile {
    my ($filePath, $content, $append) = @_;
    my $openingMode = ($append) ? ">>" : ">";
    open(FILE, $openingMode, $filePath) || die "Error while opening $filePath\n";
	print FILE $content."\n";
	close FILE;
}

sub buildQuantiAnnot {
    my $softwareVersion = &promsMod::cleanParameters(param('softversion'));
    my ($libCutOff, $libSearchFiles, $libIncludeAmbigous) = &promsMod::cleanParameters(param('libCutoff'), param('libSearchFiles'), param('libIncludeAmbigous'));
    my ($ms1PrecursorCharge, $ms1IsotopePeaks, $ms1MassAnalyzer, $ms1Peaks, $ms1Power, $ms1PowerCharge) = &promsMod::cleanParameters(param('ms1PrecursorCharge'), param('ms1IsotopePeaks'), param('ms1MassAnalyzer'),  param('ms1Peaks'),  param('ms1Power'),  param('ms1PowerCharge'));
    my ($rtScans, $rtScansFilter) = &promsMod::cleanParameters(param('rtFiltering'), param('rtFilterTime'));
    my ($fastaEnzyme, $fastaMaxMC, $fastaComment, $fastaProteinFilter, $fastaProteinFiltered) = &promsMod::cleanParameters(param('fastaEnzyme'), param('fastaMC'), param('fastaComment'), param('fastaProteinFilter'), param('fastaProteinFiltered'));
    my ($otherModifs, @modifs) = (param('otherModifs'), param('modifs'));
    
    $libIncludeAmbigous = ($libIncludeAmbigous && $libIncludeAmbigous eq "libIncludeAmbigous") ? 1 : 0;
    
    my $quantiAnnot = "LABEL=FREE::SOFTWARE=SKY;".$softwareVersion."::LIB_CUTOFF=".$libCutOff."::LIB_SEARCH_FILES=".$libSearchFiles."::LIB_AMBIGOUS=".$libIncludeAmbigous."::";
    $quantiAnnot .= "MS1_PRECURSOR_CHARGE=".$ms1PrecursorCharge."::MS1_ISOTOPE_PEAKS=".$ms1IsotopePeaks."::MS1_MASS_ANALYZER=".$ms1MassAnalyzer."::MS1_PEAKS=".$ms1Peaks."::MS1_POWER=".$ms1Power."::MS1_POWER_CHARGE=".$ms1PowerCharge."::";
    
    if($rtScans eq "UseFilteredScans") {
        $quantiAnnot .= "RT_SCANS=FILTERED::RT_SCANS_FILTER=".$rtScansFilter."::";
    } else {
        $quantiAnnot .= "RT_SCANS=ALL::";
    }
    
    $quantiAnnot .= "FASTA_COMMENT=".$fastaComment."::FASTA_ENZYME=".$fastaEnzyme."::FASTA_MAX_MC=".$fastaMaxMC."::";
    if($fastaProteinFiltered && $fastaProteinFiltered eq "fastaProteinFiltered") {
        $quantiAnnot .= "FASTA_FILTER=".$fastaProteinFilter."::";
    }
    
    my $modifs = join(',', @modifs);
    $modifs =~ s/,others//;
    $modifs .= ",".$otherModifs if($otherModifs);
    $quantiAnnot .= "MODIFS=$modifs";
    
    return $quantiAnnot;
}

sub printImportForm {
    print qq
    |<!DOCTYPE html>
    <HTML lang="fr">
        <HEAD>
            <TITLE>Import TDA Data</TITLE>
            <LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
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
            <SCRIPT>|;
    &promsMod::popupInfo();
    &promsMod::browseDirectory_JavaScript if $promsPath{'shared'};
    print qq 
    |
            </SCRIPT>
        </HEAD>
        <BODY>
            <SPAN class="title1" style="text-align: center; display: block; margin: 25px 0 30px 0;">Select options to import Skyline data</SPAN>
            <FORM method="POST" action="./importTDAData.cgi" name="parameters" enctype="multipart/form-data">
                <INPUT name="ID" type="hidden" value="$branchID">
                <INPUT name="id_project" type="hidden" value="$projectID">
                <TABLE style="background-color:$darkColor">
                    <TR>
                        <TH style="text-align:right; vertical-align:top">Result file : </TH>
                        <TD style="background-color:$lightColor">
                            <INPUT type="file" name="uploadfile" accept="$fileFormat" required /><BR/><BR/>
                            <U>Default required columns are</U> : <BR/>
                            "Protein Name", "Replicate Name", "File Name", "Peptide Modified Sequence Unimod Ids" OR "Peptide Modified Sequence Monoisotopic Masses"
                            "Precursor Mz", "Precursor Charge", "Missed Cleavages", "Peptide Retention Time", "Total Area", <BR/>
                            "Begin Pos", "End Pos", "Previous Aa", "Next Aa"<BR/><BR/>
                            <U>For MS/MS data, following columns are also required :</U><BR/>
                            "Product Mz", "Product Neutral Mass", "Product Charge", "Fragment Ion", "Retention Time", "Area"<BR/><BR/>
                        </TD>
                    </TR>
                    <TR>
                        <TH style="text-align:right; vertical-align:top">Software version : </TH>
                        <TD style="background-color:$lightColor">
                            <INPUT type="text" class="number" size="6" name="softversion" />
                            <SMALL>&nbsp;ex : 1.2, 2.1.3 ...</SMALL>
                        </TD>
                    </TR>
                    <TR>
                        <TH style="text-align:right; vertical-align:top">Library: </TH>
                        <TD style="background-color:$lightColor">
                            cut-off score : <INPUT type="text" class="number" value="95" size="3" value="0" name="libCutoff" /><br/>
                            search files : <INPUT type="text" size="30" name="libSearchFiles" /><br/>
                            <label><INPUT type="checkbox" name="libIncludeAmbigous" value="libIncludeAmbigous" />include ambigous matches</label><br/>
                        </TD>
                    </TR>
                    <TR>
                        <TH style="text-align:right; vertical-align:top">MS1 Filtering: </TH>
                        <TD style="background-color:$lightColor">
                            precursor charges: <INPUT type="text" class="number" size="7" value="2,3,4,5" name="ms1PrecursorCharge" /><br/>
                            isotope peaks included :
                            <SELECT name="ms1IsotopePeaks">
                                <OPTION value="None">None</OPTION>
                                <OPTION value="Count" selected>Count</OPTION>
                                <OPTION value="Percent">Percent</OPTION>
                            </SELECT><br/>
                            precursor mass analyzer :
                            <SELECT name="ms1MassAnalyzer">
                                <OPTION value="Centroided">Centroided</OPTION>
                                <OPTION value="QIT">QIT</OPTION>
                                <OPTION value="TOF">TOF</OPTION>
                                <OPTION value="Orbitrap" selected>Orbitrap</OPTION>
                                <OPTION value="FT-ICR">FT-ICR</OPTION>
                            </SELECT><br/>
                            peaks : <INPUT type="text" class="number" size="2" value="3" name="ms1Peaks" /><br/>
                            Resolving power <INPUT type="text" class="number" size="6" value="120000" name="ms1Power" /> at <INPUT type="text" class="number" size="3" value="200" name="ms1PowerCharge" /> m/z<br/>
                        </TD>
                    </TR>
                    <TR>
                        <TH style="text-align:right; vertical-align:top">Retention time: </TH>
                        <TD style="background-color:$lightColor">
                            <label><INPUT type="radio" name="rtFiltering" value="UseFilteredScans" onclick="" checked>Use only scans within <INPUT type="text" class="number" size="3" value="0" name="rtFilterTime" /> minutes of MS/MS IDs</label><br/>
                            <label><INPUT type="radio" name="rtFiltering" value="UseAllScans" onclick="">Include all matching scans</label><br/>
                        </TD>
                    </TR>
                    <TR>
                        <TH style="text-align:right; vertical-align:top">Fasta: </TH>
                        <TD style="background-color:$lightColor">
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
        |                   <SELECT name="dbID" id="dbID" required>
                                <OPTION value="">-= Select databank =-</OPTION>|;
        foreach my $dbSource (sort {$a cmp $b} keys %dbList){
            print "             <OPTGROUP label=\"$dbSource :\">";
            foreach my $dbID (sort {lc($dbList{$dbSource}{$a}) cmp lc($dbList{$dbSource}{$b})} keys %{$dbList{$dbSource}}) {
                print "             <OPTION value=\"$dbID\">$dbList{$dbSource}{$dbID}</OPTION>";
            }
            print "             </OPTGROUP>";
        }
        print "             </SELECT><br/>";
    }
    
    print qq
            |
                            Comment on Database entries selection: <INPUT type="text" value="" size="50" name="fastaComment" /><br/>
                            Enzyme : <INPUT type="text" value="trypsin" size="20" name="fastaEnzyme"  /><br/>
                            Max missed cleavage : <INPUT type="text" class="number" size="3" value="2" name="fastaMC" /><br/>
                            <label><INPUT type="checkbox" name="fastaProteinFiltered" value="fastaProteinFiltered" />Filter proteins having less than <INPUT type="text" size="3" value="0" name="fastaProteinFilter" /> peptide(s)</label><br/>
                        </TD>
                    </TR>
                    <TR>
                        <TH style="text-align:right; vertical-align:top">Modifications: </TH>
                        <TD style="background-color:$lightColor">
                            <table>
                                <tr>
                                    <td><label><INPUT type="checkbox" name="modifs" value="Methyl" />Methyl</label></td>
                                    <td><label><INPUT type="checkbox" name="modifs" value="Acetyl" />Acetyl</label></td>
                                </tr>
                                <tr>
                                    <td><label><INPUT type="checkbox" name="modifs" value="Dimethyl" />Dimethyl</label></td>
                                    <td><label><INPUT type="checkbox" name="modifs" value="Propionyl" />Propionyl</label></td>
                                </tr>
                                <tr>
                                    <td><label><INPUT type="checkbox" name="modifs" value="Trimethyl" />Trimethyl</td>
                                    <td><label><INPUT type="checkbox" name="modifs" value="Phospho" />Phospho</label><br/></td>
                                </tr>
                                <tr>
                                    <td colspan="2"><label><INPUT type="checkbox" name="modifs" value="others" />Others : <INPUT type="text" size="30" name="otherModifs" /></label></td>
                                </tr>
                            </table>
                        </TD>
                    </TR>
                    <TR>
                        <TH colspan=2>
                            <!-- SUBMIT button -->
                            <input type="submit" name="submit" value="Submit">&nbsp;&nbsp;&nbsp;
                            
                            <!-- CLEAR button -->
                            <INPUT type="reset" value="Clear" />&nbsp;&nbsp;&nbsp;
                            
                            <!-- CANCEL button -->
                            <INPUT value="Cancel" onclick="cancelAction();" type="button" >
                        </TH>
                    </TR>
                </TABLE>
            </FORM>
            
            <!-- Required for description popups -->
            <DIV id="divDescription" class="clDescriptionCont"></DIV>
            <SCRIPT>setPopup();</SCRIPT>
        </BODY>
    </HTML>
    |;
}


####>Revision history<####
# 1.0.3 Added RT for MS1 TDA (VS 27/11/18)
# 1.0.2 Fix output redirection (VS 23/11/18)
# 1.0.1 Added skyline parameters (VS 19/11/18)
# 1.0.0 Created (VS 12/11/18)