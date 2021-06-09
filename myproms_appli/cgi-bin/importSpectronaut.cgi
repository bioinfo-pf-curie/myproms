#!/usr/local/bin/perl -w

################################################################################
# importSpectronaut.cgi         1.1.6                                          #
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

use lib qw(.);
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use strict;
use IO::Handle;
use promsConfig;
use promsMod;
use promsDIARefonte; 
use promsImportDataManager;
use promsQuantif;
use XML::Simple;
use POSIX qw(strftime); # to get the time
use File::Path qw(rmtree); # remove_tree
use File::Spec::Functions qw(splitpath); # Core module
use File::Copy;
use File::Basename;
use String::Util qw(trim);
use List::Util qw(min max);
use Data::Dumper;

no warnings qw( once );
$CGI::LIST_CONTEXT_WARN = 0; # Disable warning of param values used as array

$|=1; # buffer flush (allows instant printing)

##########################
### LOAD CONFIGURATION ###
##########################
my %promsPath=&promsConfig::getServerInfo('no_user');
my %clusterInfo=&promsConfig::getClusterInfo; #('debian'); # default is 'centos'
my ($lightColor,$darkColor)=&promsConfig::getRowColors;


##########################
### PARSING PARAMETERS ###
##########################

# Global parameters
my $userID = (param('USERID')) ? param('USERID') : ($ENV{REMOTE_USER}) ? $ENV{REMOTE_USER} : '';
my $dbh = &promsConfig::dbConnect($userID);
my $action = (param('ACT')) ? param('ACT') : "" ; ## import or quantification

if($action eq 'ajaxUpdateDBList') {
    &ajaxUpdateDBList();
    exit;
}

my $projectID = param('id_project');
my $branchID = param('ID');
my @libIDs = (param('libID')) ? (param('libID')) : ();
@libIDs = split(/,/, $libIDs[0]) if(@libIDs && $libIDs[0] =~ /,/);

my ($item, $experimentID) = split(':', $branchID);
$experimentID = (param('expID')) ? param('expID') : $experimentID;

my %massATave = &promsConfig::getMassATave; # Average mass, needed for protein mass calculation
my %massAAave = &promsConfig::getMassAAave; # Average mass, needed for protein mass calculation

# Run options
my $runImport = (param('IMPORT')) ? param('IMPORT') : 0;
my $submit = (param('submit')) ? param('submit') : "";

# All import form fields
my $importRaw = (param('importRaw')) ? param('importRaw') : ""; # Should import MS1 quantification (raw + Normalized)
my $importCandidates = (param('importCandidates')) ? param('importCandidates') : ""; # Should import candidates table (from file)
my $createProtAbundQuanti = (param('createProtAbundQuanti')) ? param('createProtAbundQuanti') : ""; # Should create protein abundance quantification from MS1 normalized data
my $createProtRatio = (param('createProtRatio')) ? param('createProtRatio') : ""; # Should create protein ratio based on spectronaut exported candidates file
my $expParamFile = (param('expParamFile')) ? param('expParamFile') : '';
my $expParamFileShared = (param('sharedDirFilesParams')) ? param('sharedDirFilesParams') : ''; # Same as $expParamFile but coming from shared directory
my $fragmentsFile = (param('fragmentsFile')) ? param('fragmentsFile') : ''; # All fragments information
my $fragmentsFileShared = (param('sharedDirFilesFragments')) ? param('sharedDirFilesFragments') : ''; # Same as $fragmentsFile but coming from shared directory
my $candidatesFile = (param('candidatesFile')) ? param('candidatesFile') : ''; # All candidates protein quantity and ratios
my $candidatesFileShared = (param('sharedDirFilesCandidates')) ? param('sharedDirFilesCandidates') : ''; # Same as $candidatesFile but coming from shared directory 

# Data variables
my (%libsInfo); # Library statistics
my ($instrument, $taxonomy); # Library metadatas
my (%massMods, %modsInfosID, %modifIDs); # Libraries modifications infos and myProMS modification ID

my $quantifAnnot = ""; # All experiment information (extracted from libraries and experiment parameters file) [Used for MS1 and MS2 quantification data]
my (%spectronautQuantiParams); # Parameters to specify in the Label Free quantification "quantifAnnot" DB field
$spectronautQuantiParams{"intensity_metric"} = "MEAN_INT"; # Abundance metric used for minor grouping (values: mean/geom.mean/median/sum; default: mean intensity)
$spectronautQuantiParams{"quantity_level"} = "MS1"; # Level used to perform quantification (values: MS1/MS2; default: MS1)
$spectronautQuantiParams{"grouping_level"} = ""; # Level used to group post-analysis data (values: protein/peptide; default: protein)
$spectronautQuantiParams{"minor_grouping"} = 'modifSeq'; # Minor grouping strategy (values: modifSeq/strippedSeq/precursor; default: modifSeq)
$spectronautQuantiParams{"proteotypicity"} = "All"; # Wether to use proteotypic peptides or all of them (values: Only Proteotypic/All; default: All)
$spectronautQuantiParams{"probability_cutoff"} = 0; # PTM probability threshold used to do quantification (values: 0-100; default: 0; Optional in Spectronaut)

my (%peptideInfo); # All peptides information from fragment file content; All peptides PTM IDs

# Databanks information (IDs, rank, proteins content and statistics)
my @dbIDs = (param('dbID')) ? (param('dbID')) : ();
@dbIDs = split(/,/, $dbIDs[0]) if(@dbIDs && $dbIDs[0] =~ /,/);
my (%allProts, %proteinsDB, %protAnaCoverage, %dbRankProt); 

my (%quantiAbund, %stateQuantiRatio); # Quantification data at the protein/site level

my (%numMatchProt);

my (%designExperimentRatio, %designExperimentAbund, %peptideModif, %peptideProt, %peptidePosition, %fragmentsInfos, %decoyNb, %targetNb, %sheetNames, %scoreList, %analysisFileNames, %sampleFileNames, %peptideList, %assosProtPeptide, %pepExperimentDesign);


# Working environnement
my $time = strftime("%Y%m%d%H%M%S", localtime); # "20200514154020";
my $workDir = (param('workdir')) ? param('workdir') : "$promsPath{data}/tmp/Swath/$time";

# Log files
my $fileStat = "$workDir/status.out";
my $fileError = "$workDir/status\_error.out";
my $fileLog = "$workDir/log.out";
my $warning = "";

print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);

# Start import job on cluster given all form parameters
if($submit && $action eq "import") {
    ### Create import directory on server
    $CGITempFile::TMPDIRECTORY = $workDir; # folder where GCI temporary stores files
    mkdir $workDir unless -e $workDir;
    
    ### Print waiting screen
    print qq |
        <HTML>
            <HEAD>
                <LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
            </HEAD>
            <BODY background="$promsPath{images}/bgProMS.gif">
            <CENTER><BR/><BR/>
                <IMG src='$promsPath{images}/engrenage.gif'><BR/><BR/>
                <FONT color=\"red\"><DIV id=\"waitDIV\"></DIV></FONT><BR/>
                <B>Processing</B><BR/><BR/>
    |;
    
    ### Upload files on server
    ##### Fragment data file
    if ($fragmentsFile && !-e "$workDir/$fragmentsFile") {
        $fragmentsFile = uploadFile($workDir, $fragmentsFile);
    } elsif($fragmentsFileShared) {
        foreach my $file ($fragmentsFileShared) {
            $fragmentsFile = uploadFile($workDir, $file, 'shared');
            last;
        }
    }
    
    ##### Experiment parameters file
    if ($expParamFile  && !-e "$workDir/$expParamFile") {
        $expParamFile = uploadFile($workDir, $expParamFile);
    } elsif($expParamFileShared) {
        foreach my $file ($expParamFileShared) {
            $expParamFile = uploadFile($workDir, $file, 'shared');
            last;
        }
    }
    
    ### Candidates table file
    if ($candidatesFile && !-e "$workDir/$candidatesFile") {
        $candidatesFile = uploadFile($workDir, $candidatesFile);
    } elsif($candidatesFileShared) {
        foreach my $file ($candidatesFileShared) {
            $candidatesFile = uploadFile($workDir, $file, 'shared');
            last;
        }
    }
    
    ### Prepare bash script to run the import
    my $currentDirectory = dirname $0;
    my $scriptFile = "$workDir/spectronautImport.sh";
    
    ##### Build query options based on selected action
    my $queryOptions = "workdir=$workDir";
    my @paramNames = param;
    foreach my $paramName (@paramNames) {
        next if($paramName eq 'submit');
        my $paramValue = param($paramName);
        
        if($paramName eq 'fragmentsFile' || $paramName eq 'sharedDirFilesFragments') {
            $paramName = 'fragmentsFile';
            $paramValue = $fragmentsFile;
        } elsif($paramName eq 'candidatesFile' || $paramName eq 'sharedDirFilesCandidates') {
            $paramName = 'candidatesFile';
            $paramValue = $candidatesFile;
        } elsif($paramName eq 'expParamFile' || $paramName eq 'sharedDirFilesParams') {
            $paramName = 'expParamFile';
            $paramValue = $expParamFile;
        } elsif($paramName eq 'ACT') {
            $paramValue = "quantification";
        } elsif($paramName eq 'libID') {
            $paramValue = join(',', @libIDs);
        } elsif($paramName eq 'dbID') {
            $paramValue = join(',', @dbIDs);
        }
        $queryOptions .= ' '.$paramName.'="'.$paramValue.'"';
    }
    
    open (BASH,"+>", $scriptFile);
    print BASH "#!/bin/bash\n";
    print BASH "export LC_ALL='C'\n";
    print BASH "cd $currentDirectory\n";
    print BASH "echo \"$clusterInfo{path}{perl}/perl $0 IMPORT=1 $queryOptions >> $fileLog 2>&1\" > $workDir/query.txt\n";
    print BASH "$clusterInfo{path}{perl}/perl $0 IMPORT=1 $queryOptions > $workDir/log.txt 2>&1\n";
    close BASH;
    chmod 0775, $scriptFile;
    
    my $expTypeStr = (@libIDs) ? "ID_LIBRARY=".join(",", @libIDs) : "ID_DATABANK=".join(",", @dbIDs);
    
    ### Run either on cluster or locally
    if ($clusterInfo{'on'}) {
        my $sizeFragmentFile = (-e "$workDir/$fragmentsFile") ? `stat -c "%s" $workDir/$fragmentsFile`: 1073741824;
        my $sizeCandidateFile = (-e "$workDir/$candidatesFile") ? `stat -c "%s" $workDir/$candidatesFile`: 1073741824;
        my $maxMem = max(10, (($sizeFragmentFile+$sizeCandidateFile) / 1073741824) * 7);
        $maxMem = 120 if($maxMem > 120);
        
        my %jobParams = (
            maxMem => sprintf("%.0f", $maxMem).'Gb',
            numCPUs => 1,
            maxHours => 24,
            jobName => "importSpectronaut_exp$experimentID\_$time",
            outFile => 'PBSimport.txt',
            errorFile => 'PBSimportError.txt',
            jobEndFlag => "_END_".$action."Swath_$time",
            noWatch => '1',
        );
        my ($pbsError, $pbsErrorFile, $jobClusterID) = $clusterInfo{'runJob'}->($workDir, $scriptFile, \%jobParams);
        
        # Create new job to monitor
        $dbh->do("INSERT INTO JOB_HISTORY (ID_JOB, ID_USER, ID_PROJECT, ID_JOB_CLUSTER, TYPE, JOB_STATUS, FEATURES, SRC_PATH, LOG_PATH, ERROR_PATH, STARTED) VALUES('$time', '$userID', $projectID, 'C$jobClusterID', 'Import [Spectronaut]', 'Queued', 'SOFTWARE=SPECTRONAUT;$expTypeStr;FRAGMENT_FILE=$fragmentsFile;CANDIDATES_FILE=$candidatesFile', '$workDir', '$fileStat', '$fileError', NOW())");
        $dbh->commit;
    } elsif(!param('workdir')) { # If script was not started by the cluster job, create it in DB
        $dbh->do("INSERT INTO JOB_HISTORY (ID_JOB, ID_USER, ID_PROJECT, ID_JOB_CLUSTER, TYPE, JOB_STATUS, FEATURES, SRC_PATH, LOG_PATH, ERROR_PATH, STARTED) VALUES('$time', '$userID', $projectID, 'L$$', 'Import [Spectronaut]', 'Queued', 'SOFTWARE=SPECTRONAUT;$expTypeStr;FRAGMENT_FILE=$fragmentsFile;CANDIDATES_FILE=$candidatesFile', '$workDir', '$fileStat', '$fileError', NOW())");
        $dbh->commit;
        system("bash $scriptFile&");
    }
    $dbh->disconnect;
    
    print qq |
                    <SCRIPT LANGUAGE="JavaScript">
                        top.promsFrame.selectedAction='summary';
                        parent.navFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=experiment:$experimentID&ACT=nav";
                    </SCRIPT>
                </CENTER>
            </BODY>
        </HTML>
    |;
    exit;
}

# Import quantifications results
elsif($runImport && $action eq "quantification") {
    
    # RETRIEVE EXPERIMENT/LIBRARIES/QUANTIFICATION METADATA
    
    ## Parse all libraries for information
    if(@libIDs) {
        printToFile($fileStat, "Scanning libraries : Fetching content statistics", 1);
        foreach my $libID (@libIDs) {
            ### Fetching SWATH LIB modifications
            (my $libName,$instrument,my $libTaxonomy,my $paramStrg,my $version,my $split,my $stat) = $dbh->selectrow_array("SELECT NAME,INSTRUMENT,ORGANISM,PARAM_STRG,VERSION_NAME,SPLIT,STATISTICS FROM SWATH_LIB WHERE ID_SWATH_LIB='$libID' ") or die "Couldn't prepare statement: " . $dbh->errstr;
            if($libTaxonomy && $taxonomy && $taxonomy !~ /$libTaxonomy/) {
                $taxonomy .= ", $libTaxonomy";
            }
            my $sthModifIDs=$dbh->prepare("SELECT M.ID_MODIFICATION, SLM.SPECIFICITY, PSI_MS_NAME, INTERIM_NAME, SYNONYMES, MONO_MASS, UNIMOD_ACC, PEAKVIEW_CODE FROM MODIFICATION M, SWATH_LIB_MODIFICATION SLM WHERE M.ID_MODIFICATION=SLM.ID_MODIFICATION AND SLM.ID_SWATH_LIB=$libID");
            $sthModifIDs->execute();
            while (my ($modID,$specif,$psiName,$interimName,$synonyms,$massMod,$unimodID,$peakviewCode)=$sthModifIDs->fetchrow_array) {
                my $modType="V";
                my @specifs=split(/,/,$specif);
                foreach my $aa (@specifs){
                    $modifIDs{$modID}{$modType}{$aa}=1;
                }
                
                # Link mono mass to its modification ID
                $massMods{$modID}=$massMod if $massMod;
                
                # Link all potential names and unimodID to its modification ID
                $modsInfosID{$psiName} = $modID if($psiName);
                $modsInfosID{$interimName} = $modID if($interimName);
                $modsInfosID{$unimodID} = $modID if($interimName);
                
                if($synonyms) {
                    foreach my $synonym (split(/##/, $synonyms)) {
                        $modsInfosID{$synonym} = $modID if($synonym);
                    }
                }
                
                next unless $unimodID && $modID;
            }
            $sthModifIDs->finish;
            
            ### Fetching SWATH LIB statistics
            my $xml=XML::Simple->new(KeepRoot=>1);
            my $xmlStat=$xml->XMLin($stat);
            my $pepNum=$xmlStat->{'NUM_ENTRY'}->{'NUM_PEP'};
            my $protNum=$xmlStat->{'NUM_ENTRY'}->{'NUM_PROT'};
            my $protSpeNum=$xmlStat->{'NUM_ENTRY'}->{'NUM_PROT_SPECIFIC'};
            
            push(@{$libsInfo{"ID"}}, $libID);
            push(@{$libsInfo{"split"}}, ($split) ? 1 : 0);
            push(@{$libsInfo{"name"}}, $libName);
            push(@{$libsInfo{"pepNum"}}, $pepNum);
            push(@{$libsInfo{"protNum"}}, $protNum);
            push(@{$libsInfo{"protSpeNum"}}, $protSpeNum);
            push(@{$libsInfo{"version"}}, $version);
            push(@{$libsInfo{"software"}}, ($paramStrg =~ /Spectronaut/) ? "Spectronaut" : "TPP");
        } 
    } else {
        ### Fetching all modifications infos
         my $sthModifIDs=$dbh->prepare("SELECT ID_MODIFICATION, SPECIFICITY, PSI_MS_NAME, INTERIM_NAME, SYNONYMES, MONO_MASS, UNIMOD_ACC, PEAKVIEW_CODE FROM MODIFICATION");
        $sthModifIDs->execute();
        while (my ($modID,$specif,$psiName,$interimName,$synonyms,$massMod,$unimodID,$peakviewCode)=$sthModifIDs->fetchrow_array) {
            my $modType="V";
            my @specifs= ($specif) ? split(/,/,$specif) : ();
            
            # Link mono mass to its modification ID
            $massMods{$modID}=$massMod if $massMod;
            
            # Link all potential names to its modification ID
            $modsInfosID{$psiName} = $modID if($psiName);
            $modsInfosID{$interimName} = $modID if($interimName);
            $modsInfosID{$unimodID} = $modID if($interimName);
            
            if($synonyms) {
                foreach my $synonym (split(/##/, $synonyms)) {
                    $modsInfosID{$synonym} = $modID if($synonym);
                }
            }
            
            next unless $unimodID && $modID;
        }
        $sthModifIDs->finish;
    }
    $quantifAnnot = "";

    ## Parse experiment parameters file
    my $beginSettings = 0;
    if($expParamFile) {
        printToFile($fileStat, "Fetching experiment settings", 1);
        my $currentCategory = '';
        open(PARAM, "<$workDir/$expParamFile");
        while (<PARAM>) {
            my ($paramName, $paramValue);
            next if(!$beginSettings && $_ !~ /BEGIN-SETTINGS/ && $_ !~ /^Spectronaut/);
            last if($_ =~ /END-SETTINGS/);
            $beginSettings = 1 if(!$beginSettings && $_ =~ /BEGIN-SETTINGS/);
            if($_ =~ /Spectronaut\s+(.+)$/) {
                $paramValue = $1;
                $paramName = "SOFTWARE_VERSION";
            } elsif($_ =~ /(.+)─\s+(.+)/) {
                my ($tree, $treeContent) = ($1, $2);
                
                my $treeDepth = length($tree);
                next if($treeDepth < 11);
                
                ($paramName, $paramValue) = split(/:/, $treeContent);
                $paramName =~ s/(?!\s*$)|(?!^\s*)//g;
                if($treeDepth == 11) {
                    $currentCategory = $paramName;
                } elsif($treeDepth > 11) {
                    $paramName = "$currentCategory $paramName";
                }
            }
            
            if($paramValue) {
                $paramValue =~ s/[\t\n\r]//g;
                $paramValue =~ s/(?!\s*$)|(?!^\s*)//g;
            }
            
            if($paramName) {
                $quantifAnnot .= "::$paramName=$paramValue";
                
                # Store specific paramters for protein abundance quantif annot building
                if($paramName eq 'Proteotypicity Filter') {
                    $spectronautQuantiParams{"proteotypicity"} = ($paramValue eq 'Only Proteotypic') ? "unique" : ($paramValue eq 'Only Protein Group Specific') ? "shared" : "all";
                } elsif($paramName eq 'Cross Run Normalization Normalization Strategy') {
                    $spectronautQuantiParams{"normalization_strategy"} = $paramValue;
                } elsif($paramName eq 'Cross Run Normalization Normalize on') {
                    $spectronautQuantiParams{"normalization_on"} = $paramValue;
                } elsif($paramName eq 'SOFTWARE_VERSION') {
                    $spectronautQuantiParams{"version"} = $paramValue;
                } elsif($paramName eq 'Minor Group Quantity' && $createProtAbundQuanti) {
                    $spectronautQuantiParams{"intensity_metric"} = ($paramValue =~ /geom/i) ? "GEO_MEAN_INT" : ($paramValue =~ /mean/i) ? "MEAN_INT" : ($paramValue =~ /sum/i) ? "SUM_INT" : ($paramValue =~ /median/i) ? "MEDIAN_INT" : "MEAN_INT";
                } elsif($paramName eq 'Quantity MS-Level') {
                    $spectronautQuantiParams{"quantity_level"} = $paramValue;
                } elsif($paramName eq 'Differential Abundance Grouping') {
                    $spectronautQuantiParams{"grouping_level"} = ($paramValue =~ /Minor/) ? "peptide" : "protein";
                    $spectronautQuantiParams{"top_n"} = ($spectronautQuantiParams{"grouping_level"} eq 'protein') ? $spectronautQuantiParams{"major_top_n"} : $spectronautQuantiParams{"minor_top_n"};
                } elsif($paramName eq 'Minor (Peptide) Grouping') {
                    $spectronautQuantiParams{"minor_grouping"} = ($paramValue =~ /Precursor/) ? 'precursor' : ($paramValue =~ /Stripped/) ? 'strippedSeq' : 'modifSeq';
                } elsif($paramName eq 'PTM Localization Probability Cutoff') {
                    $paramValue =~ s/,/\./;
                    $paramValue *= 100;
                    $spectronautQuantiParams{"probability_cutoff"} = $paramValue;
                } elsif($paramName eq 'Minor Group Top N Max') {
                    $spectronautQuantiParams{"minor_top_n"} = $paramValue;
                } elsif($paramName eq 'Major Group Top N Max') {
                    $spectronautQuantiParams{"major_top_n"} = $paramValue;
                }
            }
        }
        close PARAM;
    }
    
    $quantifAnnot = "LABEL=FREE::SOFTWARE=SPC;$spectronautQuantiParams{version}::NB_LIB=".((@libIDs) ? scalar @libIDs : 0)."::$quantifAnnot";
    
    if(@libIDs) {
        $quantifAnnot = $quantifAnnot."::LIBRARY_NAME=".join(',', @{$libsInfo{"name"}})."::LIBRARY_PEPTIDES=".join(',', @{$libsInfo{"pepNum"}})."::LIBRARY_PROTEINS=".join(',', @{$libsInfo{"protNum"}})."::LIBRARY_SPECIFICS_PROTEINS=".join(',', @{$libsInfo{"protSpeNum"}})."::LIBRARY_VERSIONS=".join(',', @{$libsInfo{"version"}});
    }
    
    if(!$spectronautQuantiParams{"grouping_level"} || !$spectronautQuantiParams{"quantity_level"}) {
        printToFile($fileError, "Could not retrieve all parameters. Check the provided file for missing data.", 1);
        exit;
    }
    
    ##########################################################
    ###  Recovering peptide miss cleavage and specificity  ###
    ##########################################################
    #printToFile($fileStat, "Scanning library : Recovering peptide miss cleavage(s) and specificity", 1);
    #if(%libsInfo) {
    #    parseLibraryInfo(\%libsInfo);
    #}
    
    ###############################
    ###  Parsing fragment file  ###
    ###############################
    $dbh->disconnect; # Disconnect to avoid DB timeout
    printToFile($fileStat, "Parsing input files : Extracting fragments information", 1);
    if($fragmentsFile && -e "$workDir/$fragmentsFile") {
        extractFragmentFileContent("$workDir/$fragmentsFile");
    } else {
        printToFile($fileError, "Fragment file is not accessible", 1);
        exit;
    }

    # Update quantifAnnot with variable modifications found (either in DB for classic DIA or in peptides sequence for DirectDIA)
    $dbh = &promsConfig::dbConnect($userID);
    my $modIDsStr = "";
    foreach my $modID (sort keys %modifIDs) {
        my ($modifPsiName, $modifInterName, $modifSynonymes) = $dbh->selectrow_array("SELECT PSI_MS_NAME, INTERIM_NAME, SYNONYMES FROM MODIFICATION WHERE ID_MODIFICATION=$modID");
        my $residues = join(', ', keys %{$modifIDs{$modID}{'V'}});
        $residues =~ s/=/N-term/;
        $modIDsStr .= ', ' if($modIDsStr);
        $modIDsStr .= ($modifPsiName) ? $modifPsiName : ($modifInterName) ? $modifInterName : $modifSynonymes;
        $modIDsStr .= " ($residues)";
    }
    $quantifAnnot .= "::Variable Modifications=$modIDsStr";
    

    #################################
    ###  Parsing candidates file  ###
    #################################
    $dbh->disconnect; # Disconnect to avoid DB timeout
    printToFile($fileStat, "Parsing input files : Extracting candidates information", 1);
    if($candidatesFile) {
        if(-e "$workDir/$candidatesFile") {
            extractCandidatesFileContent("$workDir/$candidatesFile");
        } else {
            printToFile($fileError, "Candidates file is not accessible", 1);
            exit;
        }
    }
    
    ####################################################################
    ###  Fetching protein sequence and mass from libraries Databanks ###
    ####################################################################
    printToFile($fileStat, "Databank : Extracting proteins and statistics", 1);
    
    ### Fetch DBs file
    $dbh = &promsConfig::dbConnect($userID);
    my $dbRank = 0;
    foreach my $libID (@libIDs) {
        if($libID) {
            my $sthDBID = $dbh->prepare("SELECT DISTINCT ID_DATABANK FROM DATABANK_SWATHLIB WHERE ID_SWATH_LIB=?");
            $sthDBID->execute($libID);
            while (my ($dbID) = $sthDBID->fetchrow_array) {
                if(!grep(/^$dbID&/, @dbIDs)) {
                    $dbRank++;
                    push(@dbIDs, "$dbID&$dbRank");
                }
            }
            $sthDBID->finish;
        }
    }
    
    ### Extract proteins and statistics from DB
    for(my $i=0; $i<scalar @dbIDs; $i++) {
        my $dbID;
        if(@libIDs) {
            ($dbID, $dbRank) = split(/&/, $dbIDs[$i]);
        } else {
            ($dbID, $dbRank) = ($dbIDs[$i], $i+1);
            $dbIDs[$i] = "$dbIDs[$i]&$i";
        }
        
        extractDBContent($dbID, $dbRank, \%proteinsDB, \%dbRankProt); # needs $dbh
    }
    $dbh->disconnect;
    
    #####################################################################
    ####>Recovering peptide positions on protein and coverage calcul<####
    #####################################################################
    printToFile($fileStat, "Databank : Recovering peptides position on protein", 1);
    my %protNotFound = matchPeptidesToProteins(\%proteinsDB, \%assosProtPeptide);
    my $nbProtNotFound = scalar keys %protNotFound;
    if($nbProtNotFound > 0) {
        $warning = "$nbProtNotFound proteins in result files were not found in databank: ".join(', ', keys %protNotFound);
    }

    ###############################
    ### Import data in database ###
    ###############################
    my $quantiName = ($importRaw) ? 'Spectronaut Raw Quantification (MS2)' : '';
    my $importManager = promsImportDataManager->new($userID, $projectID, $experimentID, "spectronaut", $quantiName, $quantifAnnot, $taxonomy, $instrument, $fragmentsFile, $fileStat, \@dbIDs);
    $importManager->setFilesInfos(\%analysisFileNames, \%sampleFileNames);
    $importManager->setDIAInfos(join(',', @libIDs));
    $importManager->setProtCoverage(\%protAnaCoverage);
    $importManager->setDesignExperiment((%designExperimentRatio) ? \%designExperimentRatio : \%designExperimentAbund);
    
    printToFile($fileStat, "Storing : Insert results data in database", 1);
    $importManager->importData(\%peptideInfo, \%peptideProt, \%peptideModif, \%peptidePosition, \%peptideList, \%fragmentsInfos, \%modifIDs, \%assosProtPeptide, \%dbRankProt);
    
    ## Create experiment design    
    my $designRatioIDs = $importManager->createExperimentDesign(\%designExperimentRatio) if($importCandidates);
    my $designAbundIDs = $importManager->createExperimentDesign(\%designExperimentAbund) if($createProtAbundQuanti);
    
    ## Get all analysisID
    my %analysisNameID = %{$importManager->getAnalysisID};
    my @analysisIDs = ();
    foreach my $analysis (sort { lc($a) cmp lc($b)} keys %analysisNameID) {
        push(@analysisIDs, $analysisNameID{$analysis});
    }
    
    my %peptidesID = %{$importManager->getPeptidesID};
    my %proteinsID = %{$importManager->getProteinsID};
    my $quantiMS2RawID = $importManager->getQuantiID;
    
    
    ##################################
    ### Create all quantifications ###
    ##################################
    printToFile($fileStat, "Storing : Create quantifications", 1);
    my ($quantiProtAbundanceID, $quantiProtRatioID, %anaStateAbundInfo, %anaStateRatioInfo);
    $dbh = &promsConfig::dbConnect($userID);

    ## MS2 Normalized quantification (MS2 Raw quantification was created before with the $importManager->importData function)
    my $quantiMethodDIA = $dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='DIA'");
    $quantiName = 'Spectronaut Normalized Quantification (MS2)';
    my $quantiMS2NormID = $importManager->createQuantification($quantiName, 'peptide', $quantifAnnot, $quantiMethodDIA, \@analysisIDs);

    ## MS1 Quantification
    my $quantiMethodXIC = $dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='XIC'");
    my ($quantiRawMS1ID, $quantiNormMS1ID);
    
    ## Add MS1 normalized quantification data
    $quantiName = 'Spectronaut Normalized Quantification (MS1)';
    $quantiNormMS1ID = $importManager->createQuantification($quantiName, 'peptide', $quantifAnnot, $quantiMethodXIC, \@analysisIDs);
    
    if($importRaw) {
        ## Add MS1 raw quantification data
        $quantiName = 'Spectronaut Raw Quantification (MS1)';
        $quantiRawMS1ID = $importManager->createQuantification($quantiName, 'peptide', $quantifAnnot, $quantiMethodXIC, \@analysisIDs);
    }
    
    if($createProtAbundQuanti || $importCandidates) {
        # Queries global to abundance and ratio quantifications
        my $sthProtQ = $dbh->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_PROTEIN, ID_QUANTIFICATION, SEQ_BEG, SEQ_LENGTH, ID_QUANTIF_PARAMETER, QUANTIF_VALUE, TARGET_POS) VALUES (?,?,?,?,?,?,?)");
        my $intMetricID = $dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE='".$spectronautQuantiParams{"intensity_metric"}."'");
        my $sthQP = $dbh->prepare("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE=? AND ID_QUANTIFICATION_METHOD=?");
        my $sthMMod = $dbh->prepare("SELECT DISTINCT PM.ID_MODIFICATION FROM PEPTIDE_PROTEIN_ATTRIB PPA INNER JOIN PEPTIDE_MODIFICATION PM ON PM.ID_PEPTIDE=PPA.ID_PEPTIDE WHERE PPA.ID_ANALYSIS IN (".(join(', ', @analysisIDs)).") GROUP BY PM.ID_MODIFICATION ORDER BY PM.ID_MODIFICATION");
        my $sthParentQ = $dbh->prepare("INSERT INTO PARENT_QUANTIFICATION (ID_PARENT_QUANTIFICATION, ID_QUANTIFICATION) VALUES (?, ?)");
        
        # Queries only used if quantification grouping is on minor group
        my $sthQMod = $dbh->prepare("INSERT INTO MULTIMODIF_QUANTIFICATION (ID_MODIFICATION,ID_QUANTIFICATION, MODIF_RANK) VALUES (?,?,?)");

        # Query to retrieve proteins/peptides to check for related quantification
        my $proteotypityFilter = ($spectronautQuantiParams{"proteotypicity"} eq 'unique') ? 'AND AP.PEP_SPECIFICITY=100' : '';
        my $groupingFilter = ($spectronautQuantiParams{"grouping_level"} eq 'protein') ? 'GROUP BY P.ID_PROTEIN' : ($spectronautQuantiParams{"minor_grouping"} eq 'strippedSeq') ? 'GROUP BY P.ID_PROTEIN, PEP.PEP_SEQ, PPA.PEP_BEG, PPA.PEP_END' : 'GROUP BY P.ID_PROTEIN, PEP.ID_PEPTIDE, PPA.PEP_BEG, PPA.PEP_END';
        my $joinType = ($spectronautQuantiParams{"grouping_level"} eq 'protein') ? 'LEFT' : 'INNER';
        my $order = ($spectronautQuantiParams{"grouping_level"} eq 'protein') ? 'AP.ID_PROTEIN' : ($spectronautQuantiParams{"minor_grouping"} eq 'strippedSeq') ? 'P.ID_PROTEIN, PEP.PEP_SEQ, PPA.PEP_BEG, PPA.PEP_END' : 'P.ID_PROTEIN, PEP_SEQ, VMOD, PPA.PEP_BEG, PPA.PEP_END';
        my $queryGrpSelect = "SELECT P.ID_PROTEIN, P.ALIAS, PEP.PEP_SEQ, GROUP_CONCAT(DISTINCT PEPM.ID_MODIFICATION, ':', PEPM.POS_STRING SEPARATOR '&') AS VMOD, PEP_BEG AS PEP_PROT_POS
                    FROM ANALYSIS A 
                    INNER JOIN ANALYSIS_PROTEIN AP ON AP.ID_ANALYSIS=A.ID_ANALYSIS 
                    INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON PPA.ID_ANALYSIS=A.ID_ANALYSIS AND PPA.ID_PROTEIN=AP.ID_PROTEIN 
                    INNER JOIN PROTEIN P ON P.ID_PROTEIN=PPA.ID_PROTEIN 
                    INNER JOIN PEPTIDE PEP ON PEP.ID_PEPTIDE=PPA.ID_PEPTIDE 
                    $joinType JOIN PEPTIDE_MODIFICATION PEPM ON PEPM.ID_PEPTIDE=PEP.ID_PEPTIDE 
                    $joinType JOIN MODIFICATION M ON M.ID_MODIFICATION=PEPM.ID_MODIFICATION
                    WHERE AP.VISIBILITY=2 AND (CONF_LEVEL=0 || ((CONF_LEVEL=1 || CONF_LEVEL=2) && PEP.VALID_STATUS > 0)) $proteotypityFilter";
        
        ## Protein Abundance Quantification
        if($designAbundIDs && $createProtAbundQuanti) {
            printToFile($fileStat, "Storing : Create quantifications (".$spectronautQuantiParams{"grouping_level"}." abundance)", 1);
            
            my $quantiMethodProtAbundance = $dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='PROT_ABUNDANCE'");
            my $specQuantifAnnot = buildSpectronautQuantifAnnot(\%spectronautQuantiParams, $designAbundIDs, ($spectronautQuantiParams{"quantity_level"} eq 'MS2') ? $quantiMS2NormID : $quantiNormMS1ID, 'abundance');
            my $quantiLevel = ($spectronautQuantiParams{"grouping_level"} eq 'peptide') ? 'Peptide' : 'Protein';
            $quantiName = "Spectronaut $quantiLevel Abundance";
            $quantiProtAbundanceID = $importManager->createQuantification($quantiName, 'protein', $specQuantifAnnot, $quantiMethodProtAbundance, \@analysisIDs, $designAbundIDs->{"ID"});
            
            my $runDir="$promsPath{quantification}/project_$projectID/quanti_$quantiProtAbundanceID";
            system("mkdir -p $runDir") if(!-e $runDir);
            my $dbhLite=&promsQuantif::dbCreateProteinQuantification($quantiProtAbundanceID,$runDir); # SQLite
            my $sthInsProt=$dbhLite->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_PROTEIN,SITE_CODE,ID_QUANTIF_PARAMETER,QUANTIF_VALUE,TARGET_POS) VALUES (?,?,?,?,?)"); # PK autoincrement
            
            my $distPepMetricID = $dbh->selectrow_array($sthQP, {}, ('DIST_PEP_USED', $quantiMethodProtAbundance));
            my $numPepMetricID = $dbh->selectrow_array($sthQP, {}, ('NUM_PEP_USED', $quantiMethodProtAbundance));
            my $numPepTotMetricID = $dbh->selectrow_array($sthQP, {}, ('NUM_PEP_TOTAL', $quantiMethodProtAbundance));

            # Insert PTM rank into multi modif quantification table if required
            my (%modResIDs, %quantiModifsRank);
            if($spectronautQuantiParams{"grouping_level"} eq 'peptide') {
                $sthMMod->execute();
                my $modifRank = 1;
                while(my ($modID) = $sthMMod->fetchrow_array) {
                    $sthQMod->execute($modID, $quantiProtAbundanceID, $modifRank);
                    $quantiModifsRank{$modID} = $modifRank;
                    $modifRank++;
                }
            }
            
            # Get state ID by analysis
            foreach my $stateID (sort { $a <=> $b } keys %{$designAbundIDs->{"STATE"}}) {
                my $stateIndex = $designAbundIDs->{"STATE"}{$stateID}{"INDEX"};
                
                foreach my $obsID (sort { $a <=> $b } keys %{$designAbundIDs->{"STATE"}{$stateID}{"OBS"}}) {
                    my $anaID = $designAbundIDs->{"STATE"}{$stateID}{"OBS"}{$obsID};
                    @{$anaStateAbundInfo{$anaID}} = ($stateID, $stateIndex);
                }
            }
            
            foreach my $analysis (keys %analysisFileNames) {
                my $analysisID = $analysisNameID{$analysis};
                my $query = "$queryGrpSelect AND A.ID_ANALYSIS=$analysisID $groupingFilter ORDER BY $order";
                my $sthPepQM = $dbh->prepare($query);
                
                my $targetPos = @{$anaStateAbundInfo{$analysisID}}[1];
                $sthPepQM->execute();
                
                my %processedGroupItem;
                while(my ($protID, $protName, $pepSeq, $pepVMod, $pepPosProt) = $sthPepQM->fetchrow_array) {
                    if($spectronautQuantiParams{"grouping_level"} eq 'peptide') {
                        next if(!$pepVMod || $processedGroupItem{"$protID\_$pepSeq\_$pepVMod\_$pepPosProt"}); # Just in case, but should not happen
                        
                        my $pepSiteCodeStr = "";
                        my $modIndexStr = "";
                        my $pepLength = length($pepSeq);
                        
                        # Change pep vmod to SITE_CODE expected format
                        foreach my $modInfo (split(/\&/, $pepVMod)) {
                            my ($modID, $modPosAll) = split(/:/, $modInfo);
                            my $modRank = $quantiModifsRank{$modID};
                            
                            $pepSiteCodeStr .= "&" if($pepSiteCodeStr);
                            $pepSiteCodeStr .= "$modRank#";
                            
                            $modIndexStr .= "&" if($modIndexStr);
                            $modIndexStr .= "$modID#";
                            
                            my $modSiteStr = "";
                            my $modSiteIndexStr = "";
                            foreach my $modPos (split(/\./, $modPosAll)) {
                                my $residue = ($modPos =~ /^[-=+*]$/) ? $modPos : substr($pepSeq, $modPos-1, 1);
                                my $modProtPos = ($modPos =~ /^[-=+*]$/) ? $pepPosProt : $pepPosProt+$modPos-1;
                                
                                $modSiteStr .= '.' if($modSiteStr);
                                $modSiteStr .= ($residue =~ /^[-=+*]$/) ? $residue : "$residue$modProtPos";
                                
                                $modSiteIndexStr .= '.' if($modSiteIndexStr);
                                $modSiteIndexStr .= ($residue =~ /^[-=+*]$/) ? $residue : "$residue$modPos";
                            }
                            $pepSiteCodeStr .= $modSiteStr;
                            $modIndexStr .= $modSiteIndexStr;
                        }
                        
                        $pepSiteCodeStr .= "[$pepPosProt.$pepLength]";
                        $modIndexStr = "$pepSeq\_$modIndexStr";
                        
                        next if(!$quantiAbund{$analysis}{$protName}{$modIndexStr}); # No/Null intensity found for this peptide
                        
                        my %protQuantifItemsID = (
                            $intMetricID => $quantiAbund{"quantity"}{$analysis}{$protName}{$modIndexStr},
                            $numPepMetricID => scalar keys %{$quantiAbund{"allPep"}{$analysis}{$protName}{$modIndexStr}},
                            $numPepTotMetricID => scalar keys %{$quantiAbund{"allPepTot"}{$protName}{$modIndexStr}}
                            #$distPepMetricID => scalar keys %{$quantiAbund{"distPep"}{$analysis}{$protName}{$modIndexStr}}
                        );
                        
                        foreach my $metricsID (keys %protQuantifItemsID) {
                            $sthInsProt->execute($protID,$pepSiteCodeStr,$metricsID,$protQuantifItemsID{$metricsID},$targetPos);
                        }
                        $processedGroupItem{"$protID\_$pepSeq\_$pepVMod\_$pepPosProt"} = 1;
                        
                    } elsif($spectronautQuantiParams{"grouping_level"} eq 'protein') {
                        next if($processedGroupItem{$protID} || !$quantiAbund{"quantity"}{$analysis}{$protName}); # Already added this protein or it has no/null abundance quantification value
                        
                        my %protQuantifItemsID = (
                            $intMetricID => $quantiAbund{"quantity"}{$analysis}{$protName},
                            $numPepMetricID => scalar keys %{$quantiAbund{"allPep"}{$analysis}{$protName}},
                            $numPepTotMetricID => scalar keys %{$quantiAbund{"allPepTot"}{$protName}}
                            #$distPepMetricID => scalar keys %{$quantiAbund{"distPep"}{$analysis}{$protName}}
                        );
                        
                        foreach my $metricsID (keys %protQuantifItemsID) {
                            $sthInsProt->execute($protID,undef,$metricsID,$protQuantifItemsID{$metricsID},$targetPos);
                        }
                        
                        $processedGroupItem{$protID} = 1;
                    }
                }
                
                $sthPepQM->finish;
            }
            
            # Add peptide quantification to the abundance one based on quantity_level parameter defined in Spectronaut
            $sthParentQ->execute(($spectronautQuantiParams{"quantity_level"} eq 'MS2') ? $quantiMS2NormID : $quantiNormMS1ID, $quantiProtAbundanceID);
            
            # Link all states to the abundance quantification
            my $sthState = $dbh->prepare("INSERT INTO EXPCONDITION_QUANTIF (ID_EXPCONDITION,ID_QUANTIFICATION,COND_FUNCTION,QUANTIF_ELEMENT) VALUES (?,?,NULL,NULL)");
            foreach my $stateID (sort { $a <=> $b } keys %{$designAbundIDs->{"STATE"}}) {
                $sthState->execute($stateID, $quantiProtAbundanceID) || die $dbh->errstr();
            }
            $sthState->finish;
            
            $dbh->commit;
            $sthInsProt->finish;
            $dbhLite->commit;
            $dbhLite->disconnect;
        }
            
        ## Protein Ratio Quantification
        if($designRatioIDs && $importCandidates) {
            printToFile($fileStat, "Storing : Create quantifications (".$spectronautQuantiParams{"grouping_level"}." ratio)", 1);
            
            my $quantiMethodProtRatio = $dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='PROT_RATIO_PEP'");
            my $specQuantifAnnot = buildSpectronautQuantifAnnot(\%spectronautQuantiParams, $designRatioIDs, ($spectronautQuantiParams{"quantity_level"} eq 'MS2') ? $quantiMS2NormID : $quantiNormMS1ID, 'ratio');
            my $quantiLevel = ($spectronautQuantiParams{"grouping_level"} eq 'peptide') ? 'Peptide' : 'Protein';
            $quantiName = "Spectronaut $quantiLevel Ratio";
            $quantiProtRatioID = $importManager->createQuantification($quantiName, 'protein', $specQuantifAnnot, $quantiMethodProtRatio, \@analysisIDs, $designRatioIDs->{"ID"});
            
            # Insert PTM rank into multi modif quantification table if required
            my (%modResIDs, %quantiModifsRank);
            if($spectronautQuantiParams{"grouping_level"} eq 'peptide') {
                $sthMMod->execute();
                my $modifRank = 1;
                while(my ($modID) = $sthMMod->fetchrow_array) {
                    $sthQMod->execute($modID, $quantiProtRatioID, $modifRank);
                    $quantiModifsRank{$modID} = $modifRank;
                    $modifRank++;
                }
            }
            
            # Link peptide quanti to protein quanti based on Spectronaut quantity_level parameter (MS1/MS2)
            $sthParentQ->execute(($spectronautQuantiParams{"quantity_level"} eq 'MS2') ? $quantiMS2NormID : $quantiNormMS1ID, $quantiProtRatioID);
            
            my $runDir="$promsPath{quantification}/project_$projectID/quanti_$quantiProtRatioID";
            system("mkdir -p $runDir") if(!-e $runDir);
            my $dbhLite=&promsQuantif::dbCreateProteinQuantification($quantiProtRatioID,$runDir); # SQLite
            my $sthInsProt=$dbhLite->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_PROTEIN,SITE_CODE,ID_QUANTIF_PARAMETER,QUANTIF_VALUE,TARGET_POS) VALUES (?,?,?,?,?)"); # PK autoincrement
            
            # Compute Target pos by state index
            my %targetPos;
            foreach my $stateID1 (sort {$designRatioIDs->{"STATE"}{$a}{"INDEX"} <=> $designRatioIDs->{"STATE"}{$b}{"INDEX"}} keys %{$designRatioIDs->{"STATE"}}) {
                my $stateIndex1 = $designRatioIDs->{"STATE"}{$stateID1}{"INDEX"};
                my $stateName1 = $designRatioIDs->{"STATE"}{$stateID1}{"NAME"};
                
                foreach my $obsID (sort { $a <=> $b } keys %{$designRatioIDs->{"STATE"}{$stateID1}{"OBS"}}) {
                    my $anaID = $designRatioIDs->{"STATE"}{$stateID1}{"OBS"}{$obsID};
                    push(@{$targetPos{"ANA"}{$stateName1}}, $anaID);
                }
                
                foreach my $stateID2 (sort {$designRatioIDs->{"STATE"}{$a}{"INDEX"} <=> $designRatioIDs->{"STATE"}{$b}{"INDEX"}} keys %{$designRatioIDs->{"STATE"}}) {
                    my $stateIndex2 = $designRatioIDs->{"STATE"}{$stateID2}{"INDEX"};
                    my $stateName2 = $designRatioIDs->{"STATE"}{$stateID2}{"NAME"};
                    next if($stateID1 == $stateID2 || $targetPos{"POS"}{"$stateName2/$stateName1"} || $targetPos{"POS"}{"$stateName1/$stateName2"});
                    $targetPos{"POS"}{"$stateName2/$stateName1"} = (scalar keys %{$targetPos{"POS"}})+1;
                }
            }
            
            # Get state index by analysis
            my %anaStateIndex;
            foreach my $stateID (sort { $a <=> $b } keys %{$designRatioIDs->{"STATE"}}) {
                my $stateIndex = $designRatioIDs->{"STATE"}{$stateID}{"INDEX"};
                my $obsIndex = 1;
                
                foreach my $obsID (sort { $a <=> $b } keys %{$designRatioIDs->{"STATE"}{$stateID}{"OBS"}}) {
                    my $anaID = $designRatioIDs->{"STATE"}{$stateID}{"OBS"}{$obsID};
                    @{$anaStateIndex{$anaID}} = ($stateIndex, $obsIndex);
                    $obsIndex++;
                }
            }
            
            # TODO Make a proper resultPep.txt file to be able to display peptides quantity detail in showProtQuantification
            #my $quantiPath = $promsPath{"quantification"}."/project_$projectID/quanti_$quantiProtRatioID";    
            #system("mkdir -p $quantiPath/data");
            #open(DATA_TABLE, ">>", "$quantiPath/data/resultsPep.txt") or die ("open: $!");
            #open(DATA_EXCLUDED, ">>", "$quantiPath/data/excluded.txt") or die ("open: $!");
            #
            #print DATA_TABLE "ProteinName\tPeptideSequence\tPrecursorCharge\tFragmentIon\tProductCharge\tIsotopeLabelType\tCondition\tBioReplicate\tRun\tIntensity";
            #print DATA_EXCLUDED "ProteinName\tPeptideSequence\tPrecursorCharge\tFragmentIon\tProductCharge\tIsotopeLabelType\tCondition\tBioReplicate\tRun\tIntensity";
            #
            #foreach my $analysis(keys %analysisFileNames) {
            #    my $analysisID = $analysisNameID{$analysis};
            #    my ($stateIndex, $repIndex) = @{$anaStateIndex{$analysisID}};
            #    
            #    foreach my $protein (keys %{$assosProtPeptide{$analysis}}) {
            #        my $proteinID = $proteinsID{$protein};
            #        foreach my $peptide (keys %{$assosProtPeptide{$analysis}{$protein}}) {
            #            foreach my $rt (sort {$peptidesID{$analysis}{$peptide}{$a} <=> $peptidesID{$analysis}{$peptide}{$b}} keys %{$peptidesID{$analysis}{$peptide}}) {
            #                foreach my $peptideID (@{$peptidesID{$analysis}{$peptide}{$rt}}) {
            #    
            #                    # Get fragment and peptide charge
            #                    my $fragment=$fragmentType.'_'.$fragmentCharge;
            #                    $fragment .= "_$fragmentResidue" if($fragmentResidue);
            #                    
            #                    ### Store peptide infos
            #                    @{$peptideList{$anaName}{$peptide}{"$RT"}} = ($pepMZ, '', $pepRawArea, $pepArea);
            #                    
            #                    ### Store fragments infos
            #                    @{$fragmentsInfos{$anaName}{$peptide}{"$RT"}{$fragment}}=('',$area,$fragmentMZ,'',$areaNorm) if($areaNorm ne '1'); # Do not store fragments if it has no/null intensity (-> F.NormalizedPeakArea = 1 in Spectronaut)
            #                    
            #                    my $fragNormArea = @{$peptideList{$analysis}{$peptide}{$rt}}[3];
            #                    my $experiment = $pepExperimentDesign{$analysis}{$peptide}{$rt}{'experiment'};
            #                    my $condition = $pepExperimentDesign{$analysis}{$peptide}{$rt}{'condition'};
            #                    my $fraction = $pepExperimentDesign{$analysis}{$peptide}{$rt}{'fraction'};
            #                    my $replicate = $pepExperimentDesign{$analysis}{$peptide}{$rt}{'replicate'};
            #                    my $usedForQuanti = ($pepExperimentDesign{$analysis}{$peptide}{$rt}{'usedForPepQuantification'}) ? 'NA' : 'outMixBioPeptide';
            #                    (my $pepSeqStripped = $peptide) =~ s/[^A-Z]//;
            #                    
            #                    print MS1_RESULT_PEP_NORM "$protein\t$peptide\t$replicate\t1\t$proteinID\t1\t$pepSeqStripped\t$peptideID\t$fragNormArea\t$usedForQuanti\n";
            #                }
            #            }
            #        }
            #    }
            #}
            #close DATA_TABLE;
            #close DATA_EXCLUDED;
            
            foreach my $targetPosStates (sort { $targetPos{"POS"}{$a} <=> $targetPos{"POS"}{$b} } keys %{$targetPos{"POS"}}) {
                my ($stateName2, $stateName1) = split(/\//, $targetPosStates);
                my $targetPos = $targetPos{"POS"}{$targetPosStates};
                my @statesAna = @{$targetPos{"ANA"}{$stateName1}};
                push(@statesAna, @{$targetPos{"ANA"}{$stateName2}});
                
                my $query = "$queryGrpSelect AND A.ID_ANALYSIS IN (".join(',', @statesAna).") $groupingFilter ORDER BY $order";
                my $sthPepQM = $dbh->prepare($query);
                $sthPepQM->execute();
                
                my %processedGroupItem;
                while(my ($protID, $protName, $pepSeq, $pepVMod, $pepPosProt) = $sthPepQM->fetchrow_array) {
                    if($spectronautQuantiParams{"grouping_level"} eq 'peptide') {
                        next if(!$pepVMod || $processedGroupItem{"$protID\_$pepSeq\_$pepVMod\_$pepPosProt"});
                        
                        my $pepSiteCodeStr = "";
                        my $modIndexStr = "";
                        my $pepLength = length($pepSeq);
                        
                        # Change pep vmod to SITE_CODE expected format
                        foreach my $modInfo (split(/\&/, $pepVMod)) {
                            my ($modID, $modPosAll) = split(/:/, $modInfo);
                            my $modRank = $quantiModifsRank{$modID};
                            
                            $pepSiteCodeStr .= "&" if($pepSiteCodeStr);
                            $pepSiteCodeStr .= "$modRank#";
                            
                            $modIndexStr .= "&" if($modIndexStr);
                            $modIndexStr .= "$modID#";
                            
                            my $modSiteStr = "";
                            my $modSiteIndexStr = "";
                            foreach my $modPos (split(/\./, $modPosAll)) {
                                my $residue = ($modPos =~ /^[-=+*]$/) ? $modPos : substr($pepSeq, $modPos-1, 1);
                                my $modProtPos = ($modPos =~ /^[-=+*]$/) ? $pepPosProt : $pepPosProt+$modPos-1;
                                
                                $modSiteStr .= '.' if($modSiteStr);
                                $modSiteStr .= ($residue =~ /^[-=+*]$/) ? $residue : "$residue$modProtPos";
                                
                                $modSiteIndexStr .= '.' if($modSiteIndexStr);
                                $modSiteIndexStr .= ($residue =~ /^[-=+*]$/) ? $residue : "$residue$modPos";
                            }
                            $pepSiteCodeStr .= $modSiteStr;
                            $modIndexStr .= $modSiteIndexStr;
                        }
                        
                        $pepSiteCodeStr .= "[$pepPosProt.$pepLength]";
                        $modIndexStr = "$pepSeq\_$modIndexStr";
                        
                        if($stateQuantiRatio{"$stateName2/$stateName1"}{$protName}{$modIndexStr} || $stateQuantiRatio{"$stateName1/$stateName2"}{$protName}{$modIndexStr}) {
                            my %stateQuanti = ($stateQuantiRatio{"$stateName2/$stateName1"}{$protName}{$modIndexStr}) ? %{$stateQuantiRatio{"$stateName2/$stateName1"}{$protName}{$modIndexStr}} : %{$stateQuantiRatio{"$stateName1/$stateName2"}{$protName}{$modIndexStr}};
                            foreach my $metrics (keys %stateQuanti) {
                                my $metricsID = $dbh->selectrow_array($sthQP, {}, ($metrics, $quantiMethodProtRatio));
                                $sthInsProt->execute($protID,$pepSiteCodeStr,$metricsID,$stateQuanti{$metrics},$targetPos);
                            }
                            
                            $processedGroupItem{"$protID\_$pepSeq\_$pepVMod\_$pepPosProt"} = 1;
                        }
                        
                        if($stateQuantiRatio{"NUM_PEP_TOTAL"}{$protName}{$modIndexStr}) {
                            my $metricsID = $dbh->selectrow_array($sthQP, {}, ('NUM_PEP_TOTAL', $quantiMethodProtRatio));
                            $sthInsProt->execute($protID,$pepSiteCodeStr,$metricsID,$stateQuantiRatio{"NUM_PEP_TOTAL"}{$protName}{$modIndexStr},undef);
                            delete $stateQuantiRatio{"NUM_PEP_TOTAL"}{$protName}{$modIndexStr};
                        }
                        
                    } elsif($spectronautQuantiParams{"grouping_level"} eq 'protein') {
                        next if($processedGroupItem{$protID});
                        
                        if($stateQuantiRatio{"$stateName2/$stateName1"}{$protName} || $stateQuantiRatio{"$stateName1/$stateName2"}{$protName}) {
                            my %stateQuanti = ($stateQuantiRatio{"$stateName2/$stateName1"}{$protName}) ? %{$stateQuantiRatio{"$stateName2/$stateName1"}{$protName}} : %{$stateQuantiRatio{"$stateName1/$stateName2"}{$protName}};
                            foreach my $metrics (keys %stateQuanti) {
                                my $metricsID = $dbh->selectrow_array($sthQP, {}, ($metrics, $quantiMethodProtRatio));
                                $sthInsProt->execute($protID,undef,$metricsID,$stateQuanti{$metrics},$targetPos);
                            }
                            
                            $processedGroupItem{$protID} = 1;
                        }
                        
                        if($stateQuantiRatio{"NUM_PEP_TOTAL"}{$protName}) {
                            my $metricsID = $dbh->selectrow_array($sthQP, {}, ('NUM_PEP_TOTAL', $quantiMethodProtRatio));
                            $sthInsProt->execute($protID,undef,$metricsID,$stateQuantiRatio{"NUM_PEP_TOTAL"}{$protName},undef);
                            delete $stateQuantiRatio{"NUM_PEP_TOTAL"}{$protName};
                        }
                    }
                }
                
                $sthPepQM->finish;
            }
            
            # Link all states to the abundance quantification
            my $sthState = $dbh->prepare("INSERT INTO EXPCONDITION_QUANTIF (ID_EXPCONDITION,ID_QUANTIFICATION,COND_FUNCTION,QUANTIF_ELEMENT) VALUES (?,?,NULL,NULL)");
            foreach my $stateID (sort { $a <=> $b } keys %{$designRatioIDs->{"STATE"}}) {
                $sthState->execute($stateID, $quantiProtRatioID) || die $dbh->errstr();
            }
            $sthState->finish;
            
            $dbh->commit;
            $sthInsProt->finish;
            $dbhLite->commit;
            $dbhLite->disconnect;
        }
        
        $sthQMod->finish;
        $sthProtQ->finish;
        $sthParentQ->finish;
        $sthMMod->finish;
    }

    ##################################################################################
    ### Storing data files and associate them with their related quantification(s) ###
    ##################################################################################
    my $nbAnalysisFiles = scalar keys %analysisFileNames;
    my $iAnalysis = 0;
    printToFile($fileStat, "Storing : Split data files by analysis", 1);

    foreach my $analysis(keys %analysisFileNames) {
        my $analysisID = $analysisNameID{$analysis};
        my $anaDir = "$promsPath{peptide}/proj_$projectID/ana_$analysisID";
        my $expDir = "$promsPath{peptide}/proj_$projectID/exp_$experimentID";
        
        mkdir $promsPath{'peptide'} unless -e $promsPath{'peptide'};
        mkdir "$promsPath{peptide}/proj_$projectID" unless -e "$promsPath{peptide}/proj_$projectID";
        mkdir $anaDir unless -e $anaDir;
        mkdir $expDir unless -e $expDir;
        
        ## Store fragment file for all analysis
        if(-s "$workDir/$fragmentsFile" || -s "$expDir/$fragmentsFile") {
            printToFile($fileStat, "Storing : Fragment quantification file for $analysis ($iAnalysis/$nbAnalysisFiles processed)", 1);
            if(!-e "$expDir/$fragmentsFile") {
                move("$workDir/$fragmentsFile", "$expDir/$fragmentsFile");
            }
            
            symlink("$expDir/$fragmentsFile", "$anaDir/$fragmentsFile") if(-e "$expDir/$fragmentsFile");
        }
        
        ## Store analysis fragments XIC infos (Raw and Normalized data)
        printToFile($fileStat, "Storing : MS2 XIC data file for $analysis ($iAnalysis/$nbAnalysisFiles processed)", 1);
        
        ### Create MS2 quantification files (Raw + Normalized)
        my $quantiMS2RawPath = $promsPath{"quantification"}."/project_$projectID/quanti_$quantiMS2RawID";
        my $quantiMS2NormPath = $promsPath{"quantification"}."/project_$projectID/quanti_$quantiMS2NormID";
        my $swathFileRaw = "$quantiMS2RawPath/swath_ana_$analysisID.txt";
        my $swathFileNorm = "$quantiMS2NormPath/swath_ana_$analysisID.txt";
        system "mkdir -p $quantiMS2RawPath" if($importRaw);
        system "mkdir -p $quantiMS2NormPath";
        
        open(MS2_XIC_RAW, ">>", $swathFileRaw) or die ("open: $!") if($importRaw);
        open(MS2_XIC_NORM, ">>", $swathFileNorm) or die ("open: $!");
        
        print MS2_XIC_RAW "ID_PEPTIDE!Fragment MZ!Fragment Charge!Ion Type!Residue!Area!RT\n" if(-z $swathFileRaw && $importRaw);
        print MS2_XIC_NORM "ID_PEPTIDE!Fragment MZ!Fragment Charge!Ion Type!Residue!Area!RT\n" if(-z $swathFileNorm);
        foreach my $peptide (keys %{$fragmentsInfos{$analysis}}) {
            my $ghostpep=0;
            foreach my $rt (sort {$fragmentsInfos{$analysis}{$peptide}{$a} <=> $fragmentsInfos{$analysis}{$peptide}{$b}} keys %{$fragmentsInfos{$analysis}{$peptide}}){
                foreach my $fragment (keys %{$fragmentsInfos{$analysis}{$peptide}{$rt}}){
                    my ($fragmentType,$fragmentCharge,$fragmentLossType)=split(/_/,$fragment);
                    my ($ionType,$residue);
                    $fragmentType=~s/\s//;
                    if ($fragmentType=~/([a-z]{1})(.*)/) {
                        $ionType=$1;
                        $residue=$2;
                    }
                    $residue .= "-$fragmentLossType" if($fragmentLossType);
                    
                    my $fragmentRT=($fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[0]) ? $fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[0] : '';
                    my $fragmentMZ=($fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[2]) ? $fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[2] : '';
                    my $fragmentNeutralMass=$fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[3];
                    my $area=$fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[1];
                    my $areaNorm=$fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[4];
                    my $peptideID=($ghostpep) ? $peptidesID{$analysis}{$peptide}{$rt}[1] :  $peptidesID{$analysis}{$peptide}{$rt}[0] ;
                    my $fragUsedForPepQuanti = ($pepExperimentDesign{$analysis}{$peptide}{$rt}{'usedForPepQuantification'}) ? 1 : 0;
                    
                    print MS2_XIC_RAW "$peptideID!$fragmentMZ!$fragmentCharge!$ionType!$residue!$area!$fragmentRT\n" if($importRaw);
                    print MS2_XIC_NORM "$peptideID!$fragmentMZ!$fragmentCharge!$ionType!$residue!$areaNorm!$fragmentRT!$fragUsedForPepQuanti\n";
                }
                $ghostpep=1;
            }
        }
        close MS2_XIC_RAW if($importRaw);
        close MS2_XIC_NORM;
        
        symlink($swathFileRaw, "$anaDir/swath_ana_$analysisID.txt");
        symlink($swathFileNorm, "$anaDir/swath_ana_$analysisID\_norm.txt");

        
        #### Create MS1 normalized quantification file
        printToFile($fileStat, "Storing : MS1 Normalized XIC data file for $analysis ($iAnalysis/$nbAnalysisFiles processed)", 1);
        my $quantiPath = $promsPath{"quantification"}."/project_$projectID/quanti_$quantiNormMS1ID";
        system "mkdir -p $quantiPath";
        
        open(MS1_XIC_NORM, ">>", "$quantiPath/peptide_quantification.txt") or die ("open: $!");
        print MS1_XIC_NORM "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n" if (-z "$quantiPath/peptide_quantification.txt");
        foreach my $peptide (keys %{$peptideList{$analysis}}) {
            foreach my $rt (sort {$peptideList{$analysis}{$peptide}{$a} <=> $peptideList{$analysis}{$peptide}{$b}} keys %{$peptideList{$analysis}{$peptide}}){
                my $pepNormArea = @{$peptideList{$analysis}{$peptide}{$rt}}[3];
                my $peptideID = $peptidesID{$analysis}{$peptide}{$rt}[0];
                
                print MS1_XIC_NORM "2\t$peptideID\t$pepNormArea\n";
            }
        }
        close MS1_XIC_NORM;
        symlink("$quantiPath/peptide_quantification.txt", "$anaDir/MS1_norm_ana_$analysisID.txt");
        
        #### Create MS1 raw quantification file
        if($importRaw) {
            printToFile($fileStat, "Storing : MS1 Raw XIC data file for $analysis ($iAnalysis/$nbAnalysisFiles processed)", 1);
            $quantiPath = $promsPath{"quantification"}."/project_$projectID/quanti_$quantiRawMS1ID";
            system "mkdir -p $quantiPath";
            
            open(MS1_XIC_RAW, ">>", "$quantiPath/peptide_quantification.txt") or die ("open: $!");
            print MS1_XIC_RAW "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n" if (-z "$quantiPath/peptide_quantification.txt");
            foreach my $peptide (keys %{$peptideList{$analysis}}) {
                foreach my $rt (sort {$peptideList{$analysis}{$peptide}{$a} <=> $peptideList{$analysis}{$peptide}{$b}} keys %{$peptideList{$analysis}{$peptide}}){
                    my $pepRawArea = @{$peptideList{$analysis}{$peptide}{$rt}}[2];
                    my $peptideID = $peptidesID{$analysis}{$peptide}{$rt}[0];
                    
                    print MS1_XIC_RAW "2\t$peptideID\t$pepRawArea\n";
                }
            }
            close MS1_XIC_RAW;
            
            symlink("$quantiPath/peptide_quantification.txt", "$anaDir/MS1_raw_ana_$analysisID.txt");
        }
        
        $iAnalysis++;   
    }
    
    printToFile($fileStat, "Spectronaut import Ended.", 1);
    if($warning) {
        printToFile($fileStat, "WARNING : $warning", 1) if($warning);
    }
    $dbh->disconnect;
}

# No action: Display import form
else {
    printFormContent();
    exit;
}


### Functions
sub uploadFile {
    my ($destFolder, $fileName, $source) = @_;
    
    if($fileName) {
        my $newFileName = $fileName;
        $newFileName =~ s/\s/_/g;
        $newFileName =~ s/[()]//g;
        my $newFileDir = "$destFolder/".basename($newFileName);
        my $sourceFile = "";
        
        if(not -e $newFileDir) {
            $newFileName =~ s/^.+\///;
            if($source && $source eq 'shared') {
                $sourceFile = "$promsPath{shared}/$fileName";
                print("Retrieving file $newFileName from shared directory ...");
                system("cp '$sourceFile' '$newFileDir'&"); # async copy file to destination dir
            } else {
                $sourceFile = tmpFileName($fileName);
                print("Uploading file $fileName ...");
                system("mv '$sourceFile' '$newFileDir'&"); # async move file to destination dir
            }
            
            my $sourceFileSize = `stat -c %s '$sourceFile'`;
            my $targetFileSize = 0;
            while(!-e "$newFileDir" || $targetFileSize != $sourceFileSize) {
                sleep 10;
                $targetFileSize = (-e "$newFileDir") ? `stat -c %s '$newFileDir'`: 0;
                print(".");
            }
            print(" Done.</b><br/>");
        }
        
        return $newFileName; 
    }
    
    return "";
}

sub printToFile {
    my ($filePath, $content, $append) = @_;
    my $openingMode = ($append) ? ">>" : ">";
    open(FILE, $openingMode, $filePath) || die "Error while opening $filePath\n";
	print FILE $content."\n";
	close FILE;
}

sub printFormContent {
    print qq |
        <HTML>
            <HEAD>
                <LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
                <STYLE type="text/css">
                    .checklist {
                        border: 1px solid #ccc;
                        list-style: none;
                        height: 100px;
                        min-width: 250px;
                        overflow: auto;
                    }
                    .checklist, .checklist LI {margin:0; padding:0;}
                    .checklist LABEL {
                        padding-left: 3px;
                        text-indent: -25px;
                    }
                    
                    .liRow1 {
                        background: #FFFFFF;
                    }
                    .liRow2 {
                        background: #D8D8D8;
                    }
                    
                    .popup {
                        background-color: #FFFFFF;
                        border: solid 3px #999999;
                        padding: 5px;
                        box-shadow: 10px 10px 20px #808080;
                        position: absolute;
                        display: none;
                    }
                    
                    .listSamples {
                        margin-left: 20px;
                    }
                    
                    input[type="checkbox"], input[type="radio"] {
                        position: relative;
                        top: 2px;
                        margin-right: 6px;
                    }
                    
                    input[type="radio"] {
                        right: 2px;
                        margin-right: 8px;
                    }
                    
                    th {
                        min-width: 300px;
                    }
                    
                    td {
                        padding: 6px 0 5px 3px;
                    }
                </STYLE>
    |;
    
    printJS();
    
    print qq |
            </HEAD>
            <BODY background="$promsPath{images}/bgProMS.gif">
                <CENTER>
    |;
    
    printImportForm();
        
    print qq |
                </CENTER>
            </BODY>
        </HTML>
    |;
}

sub ajaxUpdateDBList {
    my $organism = (param('org')) ? param('org') : '';
    my %dbList;
    
    print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);
    
    my $sthDbList=$dbh->prepare("SELECT D.ID_DATABANK,D.NAME,VERSION_NAME,FASTA_FILE,DT.NAME,ORGANISM FROM DATABANK D,DATABANK_TYPE DT WHERE D.ID_DBTYPE=DT.ID_DBTYPE AND USE_STATUS='yes' AND (ORGANISM='$organism' OR ORGANISM='') AND (DT.NAME='UniProt - ACC' OR FASTA_FILE NOT LIKE '%mascot%')") or die "Couldn't prepare statement: " . $dbh->errstr;;
    $sthDbList->execute;
    if ($sthDbList->rows > 0) {
        while (my ($dbID,$dbName,$version,$fastaFile,$dbankType,$organism)=$sthDbList->fetchrow_array) {
            my $dbSource;
            if ($fastaFile=~/:/) {
                my ($mascotServer,$dbankDir,$fileName)=split(':',$fastaFile);
                $dbSource=$mascotServer;
                $version=$dbankDir;
            }
            else {
                if (!-e "$promsPath{banks}/db_$dbID/$fastaFile") {next;}
                $dbSource='Local';
            }
            
            $organism = 'Multi Species' if(!$organism);
            $dbList{$dbSource}{$organism}{$dbID} = $dbName;
            if($version && $dbSource eq 'Local') {
                $version = "v$version" if($version =~ /^[\d.,]+$/);
                $dbList{$dbSource}{$organism}{$dbID} .=" - $version";
            }
            $dbList{$dbSource}{$organism}{$dbID} .=" ($organism)";
            $dbList{$dbSource}{$organism}{$dbID} .=" [$dbankType]";
        }
    }
    $sthDbList->finish;
    
    my $databaseString="";
    foreach my $dbSource (sort{lc($a) cmp lc($b)} keys %dbList) {
        $databaseString.="<OPTGROUP label=\"$dbSource:\">\n";
        foreach my $organism (sort{lc($a) cmp lc($b)} keys %{$dbList{$dbSource}}) {
            foreach my $dbID (sort{lc($dbList{$dbSource}{$organism}{$a}) cmp lc($dbList{$dbSource}{$organism}{$b})} keys %{$dbList{$dbSource}{$organism}}) {
                $databaseString.="<OPTION value=\"$dbID\">$dbList{$dbSource}{$organism}{$dbID}</OPTION>\n";
            }
        }
        $databaseString.="</OPTGROUP>\n";
    }
    
    print($databaseString);
}

sub printImportForm {
    my $resultsFileFormat = ".tsv,.xls,.xlsx" ;
    my $paramsFileFormat = ".txt";
    
    print qq |
        <FONT class="title1">Select options to import Spectronaut data</FONT><BR/><BR/><BR/><BR/>
        <FORM method="POST" action="./importSpectronaut.cgi" name="parameters" enctype="multipart/form-data">
            <INPUT name="ID" type="hidden" value="$branchID">
            <INPUT name="id_project" type="hidden" value="$projectID">
            <INPUT name="ACT" type="hidden" value="import">
            <INPUT name="USERID" type="hidden" value="$userID">
            <TABLE bgcolor="$darkColor">
                <TR>
                    <TH align=right valign=top>Use a Spectral Library: </TH>
                    <TD bgcolor="$lightColor" style="padding: 0px 0px 5px 3px;">
                        <label><INPUT type="radio" name="selExpType" value="SWATH" onclick="selectExpType(this.value);" checked />Yes</label><br/>
                        <label><INPUT type="radio" name="selExpType" value="directDIA" onclick="selectExpType(this.value);" />No (DirectDIA)</label>
                    </TD>
                <TR>
                    <TH align=right valign=top><span id="expTypeLabel">Library</span>: </TH>
                    <TD bgcolor="$lightColor">
    |;
    
    my %libList;
    my $sthLibList=$dbh->prepare("SELECT NAME,ID_SWATH_LIB,PARAM_STRG,IDENTIFIER_TYPE FROM SWATH_LIB WHERE USE_STATUS='yes'") or die "Couldn't prepare statement: " . $dbh->errstr;
    $sthLibList->execute;
    if ($sthLibList->rows==0) {
        print "No library available.";
    } else {
        while (my($libName,$libID,$paramStr,$identType)=$sthLibList->fetchrow_array) {
            $identType=($identType)? $identType : 'Not defined';
            my $software = (!$paramStr) ? 'Unknown Software' : ($paramStr =~ /Spectronaut/) ? "Spectronaut" : "TPP";
            $libList{$software}{$libName}="$libID&$identType";
        }
        $sthLibList->finish;
        $dbh->commit;

        print qq |
                        <SELECT name="libID" id="libID" required multiple style="resize:vertical">
        |;
        foreach my $software (sort {lc($a) cmp lc($b)} keys %libList) {
            print "<OPTGROUP label='$software'>";
            foreach my $libName (sort {lc($a) cmp lc($b)} keys %{$libList{$software}}) {
                my ($libID,$identType)=split('&',$libList{$software}{$libName});
                print "         <option value=\"$libID\">$libName [$identType]</option>";
            }
            print "</OPTGROUP>";
        }
        
        print("</SELECT>");
    }
    
    print qq | <span id="orgSelect" style="display:none;"><SELECT name="orgID" id="orgID" required style="display:none" onchange="ajaxUpdateDBList(this.value)" disabled> |;
    my $sthOrg=$dbh->prepare("SELECT DISTINCT ORGANISM FROM DATABANK WHERE USE_STATUS='yes' AND (ORGANISM IS NOT NULL AND ORGANISM != '') ORDER BY ORGANISM") or die "Couldn't prepare statement: " . $dbh->errstr;;
    $sthOrg->execute;
    print("<OPTION value=''>--- Select a species ---</OPTION>");
    while (my ($organism)=$sthOrg->fetchrow_array) {
        print("<OPTION value='$organism'>$organism</OPTION>");
    }
    print "</SELECT></span>";
    $sthOrg->finish;
    
    print qq |
                        
                        </SELECT><br/><br/>
    
                        <SELECT name="dbID" id="dbID" required multiple style="resize:vertical;display:none" disabled>
                        </SELECT>
                    </TD>
                </TR>
                <TR>
                    <TH align=right valign=top>Fragments data file : <br/>($resultsFileFormat)</TH>
                    <TD bgcolor="$lightColor">
                        <label><INPUT type="radio" name="selSourceFragments" value="UseLocalDirectoryFragments" onclick="selectSource(this.value);" checked /><INPUT type="file" name="fragmentsFile" id="fragmentsFile" accept="$resultsFileFormat" required></label><br/>
                        <label><INPUT type="radio" name="selSourceFragments" value="UseSharedDirectoryFragments" onclick="selectSource(this.value);" /><B>Import from shared data directory </B></label><br/>
                        <DIV id="sharedDirDIVFragments" style="width:600px;max-height:300px;overflow:auto;display:none">
    |;

    (my $regexResultsFileFormat = $resultsFileFormat) =~ s/,/\|/g;
    &promsMod::browseDirectory_getFiles($promsPath{'shared'},'sharedDirFilesFragments',{fileMatch=>qr/($regexResultsFileFormat)$/i});
    
    print qq |
                        </DIV>
                    </TD>
                </TR>
                <TR>
                    <TH align=right valign=top>Experiment parameters file : </TH>
                    <TD bgcolor="$lightColor">
                    <label><INPUT type="radio" name="selSourceParams" value="UseLocalDirectoryParams" onclick="selectSource(this.value);" checked /><INPUT type="file" name="paramsFile" id="paramsFile" accept="$paramsFileFormat" required></label><br/>
                        <label><INPUT type="radio" name="selSourceParams" value="UseSharedDirectoryParams" onclick="selectSource(this.value);" /><B>Import from shared data directory </B></label><br/>
                        <DIV id="sharedDirDIVParams" style="width:600px;max-height:300px;overflow:auto;display:none">
    |;

    (my $regexParamsFileFormat = $paramsFileFormat) =~ s/,/\|/g;
    &promsMod::browseDirectory_getFiles($promsPath{'shared'},'sharedDirFilesParams',{fileMatch=>qr/($regexParamsFileFormat)$/i});
    
    print qq |
                        </DIV>
                </TR>
                <TR>
                    <TH align=right valign=top>Quantification import parameters: </TH>
                    <TD bgcolor="$lightColor">
                        <label id="createProtAbundQuanti"><INPUT type="checkbox" name="createProtAbundQuanti" checked> Import Proteins Abundance Quantification</label><br/>
                        <label><INPUT type="checkbox" id="importCandidates" name="importCandidates" onclick="toggleCandidatesSelection();"> Import Proteins Ratio</label><br/>
                        <FIELDSET id="candidatesSelectionDiv" style='display:none; margin-top: 10px;'>
                            <LEGEND><B>Candidates file ($resultsFileFormat):</B></LEGEND>
                            <span id="candidatesSelectionFile">
                                <label><INPUT type="radio" name="selSourceCandidates" value="UseLocalDirectoryCandidates" onclick="selectSource(this.value);" checked /><INPUT type="file" id="candidatesFile" name="candidatesFile" accept="$resultsFileFormat"></label><br/>
                                <label><INPUT type="radio" name="selSourceCandidates" value="UseSharedDirectoryCandidates" onclick="selectSource(this.value);" /><B>Import from shared data directory </B></label><br/>
                                <DIV id="sharedDirDIVCandidates" style="width:600px;max-height:300px;overflow:auto;display:none">
    |;

    &promsMod::browseDirectory_getFiles($promsPath{'shared'},'sharedDirFilesCandidates',{fileMatch=>qr/($regexResultsFileFormat)$/i});
    
    print qq |
                                </DIV>
                            </span>
                            <span id="candidatesSelectionTable" style="display:none"></span>
                        </FIELDSET>
                        
                        <label><INPUT type="checkbox" id="importRaw" name="importRaw"> Import Raw Quantifications</label><br/>
                    <TD>
                </TR>
                <TR>
                    <TH colspan=2>
                        <input type="submit" name="submit" value="Submit" />
                    
                        <!-- CLEAR button -->
                        &nbsp;&nbsp;&nbsp;<INPUT type="reset" value="Clear" />
                    
                        <!-- CANCEL button -->
                        &nbsp;&nbsp;&nbsp;<INPUT value="Cancel" onclick="window.location='./processAnalyses.cgi?ID=$branchID'" type="button" />
                    </TH>
                </TR>
            </TABLE>
        </FORM>
    |;
}

sub printJS {
    print qq |
        <SCRIPT LANGUAGE="JavaScript">
    |;
    
    &promsMod::popupInfo();
    &promsMod::browseDirectory_JavaScript if $promsPath{'shared'};
    
    print qq |
            function toggleCandidatesSelection() {
                var shouldDisplaySelectionTable = document.getElementById("importCandidates").checked;
                document.getElementById("candidatesSelectionDiv").style.display = (shouldDisplaySelectionTable) ? "" : "none";
            }

            function selectSource(source) {
                if(source == "UseLocalDirectoryCandidates") {
                    document.getElementById('candidatesFile').disabled = false;
                    if (document.getElementById('sharedDirDIVCandidates')) {
                        document.getElementById('sharedDirDIVCandidates').style.display='none';
                    }
                } else if(source == "UseSharedDirectoryCandidates") {
                    document.getElementById('candidatesFile').disabled = true;
                    if (document.getElementById('sharedDirDIVCandidates')) {
                        document.getElementById('sharedDirDIVCandidates').style.display='block';
                    }
                } else if(source == "UseLocalDirectoryFragments") {
                    document.getElementById('fragmentsFile').disabled = false;
                    if (document.getElementById('sharedDirDIVFragments')) {
                        document.getElementById('sharedDirDIVFragments').style.display='none';
                    }
                } else if(source == "UseSharedDirectoryFragments") {
                    document.getElementById('fragmentsFile').disabled = true;
                    if (document.getElementById('sharedDirDIVFragments')) {
                        document.getElementById('sharedDirDIVFragments').style.display='block';
                    }
                } else if(source == "UseLocalDirectoryParams") {
                    document.getElementById('paramsFile').disabled = false;
                    if (document.getElementById('sharedDirDIVParams')) {
                        document.getElementById('sharedDirDIVParams').style.display='none';
                    }
                } else if(source == "UseSharedDirectoryParams") {
                    document.getElementById('paramsFile').disabled = true;
                    if (document.getElementById('sharedDirDIVParams')) {
                        document.getElementById('sharedDirDIVParams').style.display='block';
                    }
                } 
            }
            
            function selectExpType(value) {
                if(value == "directDIA") {
                    document.getElementById('expTypeLabel').innerHTML = 'Databank';
                    document.getElementById('orgSelect').style.display = '';
                    document.getElementById('dbID').disabled = true;
                    document.getElementById('dbID').style.display = 'none';
                    document.getElementById('orgID').disabled = false;
                    document.getElementById('orgID').style.display = '';
                    document.getElementById('libID').disabled = true;
                    document.getElementById('libID').style.display = 'none';
                } else if(value == 'SWATH') {
                    document.getElementById('expTypeLabel').innerHTML = 'Library';
                    document.getElementById('libID').disabled = false;
                    document.getElementById('libID').style.display = '';
                    
                    document.getElementById('orgSelect').style.display = 'none';
                    document.getElementById('dbID').disabled = true;
                    document.getElementById('dbID').style.display = 'none';
                    document.getElementById('orgID').disabled = true;
                    document.getElementById('orgID').style.display = 'none';
                }
            }
            
            // AJAX --->
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
            
            function ajaxUpdateDBList(org) {
                var dbSelect = document.getElementById('dbID');

                if(org) {
                    var XHR = getXMLHTTP();
                    XHR.onreadystatechange=function () {
                        if(XHR.readyState==4 && XHR.responseText ){
                            dbSelect.innerHTML = XHR.responseText;
                            dbSelect.disabled = false;
                            dbSelect.style.display = "";
                        }
                    }
                    XHR.open("GET","$promsPath{cgi}/importSpectronaut.cgi?ACT=ajaxUpdateDBList&USERID=$userID&org="+ org, true);
                    XHR.send(null);
                } else {
                    dbSelect.disabled = true;
                    dbSelect.style.display = "none";
                }
            }
        </SCRIPT>
    |;
}

sub extractFragmentFileContent {
    my ($fragmentsFilePath) = @_;
    
    if(-e $fragmentsFilePath) {
        my (%colName2Index, %foundPep);
        my $protNotQuantified = 0;
        my $nFragments = `wc -l $fragmentsFilePath`;
        $nFragments = ($nFragments=~ /(\d+)/) ? $1 : 0;
        my $iFragment = 0;
        
        open (IN, "<$fragmentsFilePath");
        my $lastPeptide = '';
        my $lastAna;
        
        while (<IN>) {
            if ($.==1) {
                $_ =~ s/\s*\Z//;
                my @columns = split(/[\t]/,$_);
                my @refColNames = ('E.Name', 'R.Condition', 'R.FileName', 'R.Fraction', 'R.Label', 'R.Replicate', 'PG.Quantity', 'PEP.AllOccurringProteinAccessions', 'PEP.Quantity', 'PEP.IsProteotypic', 'PEP.UsedForProteinGroupQuantity', 'PEP.NrOfMissedCleavages', 'PEP.PeptidePosition', 'EG.IntPIMID', 'EG.IsImputed', 'EG.ModifiedSequence', 'EG.ApexRT', 'EG.PTMLocalizationProbabilities', 'EG.UsedForPeptideQuantity', 'FG.Charge', 'FG.PrecMz', 'FG.MS1Quantity', 'FG.MS1RawQuantity', 'FG.MS2Quantity', 'FG.MS2RawQuantity', 'F.Charge', 'F.FrgType', 'F.FrgNum', 'F.FrgMz', 'F.PeakArea', 'F.NormalizedPeakArea', 'F.FrgLossType', 'F.ExcludedFromQuantification');
                foreach my $i (0 .. $#columns) {
                    $colName2Index{$columns[$i]} = $i;
                }
    
                ##check columns names
                foreach my $colName (@refColNames) {
                    unless (defined $colName2Index{$colName}) {
                        printToFile($fileError, "Column '$colName' not found in fragments file !", 1);
                        exit;
                    }
                }
            } else {
                $_ =~ s/\s*\Z//;
                my @infos = split(/[\t]/,$_);
                next if(!$infos[$colName2Index{'PEP.AllOccurringProteinAccessions'}] || $infos[$colName2Index{'EG.IntPIMID'}] !~ /^_/);
    
                # Design info
                my $anaFile = $infos[$colName2Index{'R.Label'}];
                (my $anaName=$anaFile)=~s/\.raw//;
                $analysisFileNames{$anaName}=$anaFile;
                $sampleFileNames{$anaName}=$anaName;
                
                my $experiment = ($infos[$colName2Index{'E.Name'}]) ? $infos[$colName2Index{'E.Name'}] : "DIA Spectronaut";
                my $condition = ($infos[$colName2Index{'R.Condition'}] && $infos[$colName2Index{'R.Condition'}] ne 'NA') ? $infos[$colName2Index{'R.Condition'}] : "Condition 1";
                my $replicate = ($infos[$colName2Index{'R.Replicate'}]) ? $infos[$colName2Index{'R.Replicate'}] : "1";
                my $fraction = ($infos[$colName2Index{'R.Fraction'}] && $infos[$colName2Index{'R.Fraction'}] ne 'NA') ? $infos[$colName2Index{'R.Fraction'}] : undef;

                # Proteins Infos
                my $proteinsList = $infos[$colName2Index{'PEP.AllOccurringProteinAccessions'}];
                $proteinsList =~ s/"//g;
                my $proteinQuantity = $infos[$colName2Index{'PG.Quantity'}];
                if ($proteinQuantity !~ /\d+\,?\d*/) {
                    $proteinQuantity = 'NA';
                } else {
                    $proteinQuantity =~ s/,/\./;
                }
                
                # Peptides Infos
                my $pepUsedForProtQuanti = ($infos[$colName2Index{'PEP.UsedForProteinGroupQuantity'}] eq 'True') ? 1 : 0;
                my $pepUsedForPepGroupQuanti = ($infos[$colName2Index{'EG.UsedForPeptideQuantity'}] eq 'True') ? 1 : 0;
                my $pepMZ = $infos[$colName2Index{'FG.PrecMz'}];
                $pepMZ =~ s/,/\./;
                
                my $pepCharge = $infos[$colName2Index{'FG.Charge'}];
                my $RT = $infos[$colName2Index{'EG.ApexRT'}];
                $RT =~ s/,/\./;
                $RT = sprintf("%.2f", $RT);
                
                my $pepRawArea = ($spectronautQuantiParams{"quantity_level"} eq 'MS2') ? $infos[$colName2Index{'FG.MS2RawQuantity'}] : $infos[$colName2Index{'FG.MS1RawQuantity'}];
                if ($pepRawArea !~ /\d+\,?\d*/) {
                    $pepRawArea = 'NA';
                } else {
                    $pepRawArea =~ s/,/\./;
                }
                
                my $pepArea = ($spectronautQuantiParams{"quantity_level"} eq 'MS2') ? $infos[$colName2Index{'FG.MS2Quantity'}] : $infos[$colName2Index{'FG.MS1Quantity'}];
                if ($pepArea !~ /\d+\,?\d*/) {
                    $pepArea = 'NA';
                } else {
                    $pepArea =~ s/,/\./;
                }
                
                my $modifSeq=$infos[$colName2Index{'EG.IntPIMID'}];
                $modifSeq=~s/_//g;
                my $peptide=$modifSeq.'_'.$pepCharge;
                
                my $modifSeqStr=($infos[$colName2Index{'EG.ModifiedSequence'}]) ? $infos[$colName2Index{'EG.ModifiedSequence'}] : "";
                $modifSeqStr=~s/_//g;
                
                my $modifSeqStrScoreAll=($infos[$colName2Index{'EG.PTMLocalizationProbabilities'}]) ? $infos[$colName2Index{'EG.PTMLocalizationProbabilities'}] : "";
                $modifSeqStrScoreAll=~s/_//g;
                
                my $formattedSeq = ($colName2Index{'PEP.GroupingKey'} && $infos[$colName2Index{'PEP.GroupingKey'}]) ? $infos[$colName2Index{'PEP.GroupingKey'}] : ($spectronautQuantiParams{"minor_grouping"} eq 'modifSeq') ? $modifSeq : ($spectronautQuantiParams{"minor_grouping"} eq 'precursor') ? $peptide : ($modifSeq =~ s/\[[^\]]+\]|\(UniMod:\d+\)//g);
                
                my $pepMissCleavages = ($infos[$colName2Index{'PEP.NrOfMissedCleavages'}]) ? $infos[$colName2Index{'PEP.NrOfMissedCleavages'}] : 0;
                my $pepSpecificity = ($colName2Index{'PEP.IsProteotypic'} && $infos[$colName2Index{'PEP.IsProteotypic'}] eq 'True') ? 100 : 0;
                my $pepIsImputed = ($infos[$colName2Index{'EG.IsImputed'}] eq 'True') ? 1 : 0;
                my $pepScore = ($colName2Index{'EG.CScore'} && $infos[$colName2Index{'EG.CScore'}]) ? $infos[$colName2Index{'EG.CScore'}] : '';
                if ($pepScore !~ /\d+\,?\d*/) {
                    $pepScore = '';
                } else {
                    $pepScore =~ s/,/\./;
                }
                my $pepBegPositions = ($infos[$colName2Index{'PEP.PeptidePosition'}]) ? $infos[$colName2Index{'PEP.PeptidePosition'}] : "";
                $pepBegPositions =~ s/"//g;
                my @proteinsPepBegPos = split(/;/,$pepBegPositions);
                my $pepQuantity = ($infos[$colName2Index{'PEP.Quantity'}]) ? $infos[$colName2Index{'PEP.Quantity'}] : 0;
                if ($pepQuantity !~ /\d+\,?\d*/) {
                    $pepQuantity = 'NA';
                } else {
                    $pepQuantity =~ s/,/\./;
                }
                
                #next if($proteinQuantity eq '1' || $pepQuantity eq '1'); # Do not take into account peptides/proteins with no intensities (Spectronaut uses 1 as value)
                
                
                ### Set experiment design
                %{$pepExperimentDesign{$anaName}{$peptide}{"$RT"}} = (
                    'experiment' => $experiment,
                    'condition' => $condition,
                    'fraction' => $fraction,
                    'replicate' => $replicate,
                    'usedForProtQuantification' => $pepUsedForProtQuanti,
                    'usedForPepQuantification' => $pepUsedForPepGroupQuanti,
                );
                @{$peptideInfo{$peptide}} = ($pepMZ, $pepMissCleavages, $pepSpecificity) if(!$peptideInfo{$peptide});
                
                if($importCandidates) {
                    $designExperimentRatio{"name"} = "Spectronaut quantification (".$spectronautQuantiParams{"grouping_level"}." ratio)" if(!$designExperimentRatio{"name"});
                    
                    if($designExperimentRatio{"conditions"}{$condition}{$replicate} && $designExperimentRatio{"conditions"}{$condition}{$replicate}{'anaName'} ne $anaName) {
                        printToFile($fileError, "An inconsistency was found in the experimental design: $anaName is set with the same replicate index as ".$designExperimentRatio{"conditions"}{$condition}{$replicate}{'anaName'}. " (replicate $replicate)", 1);
                        exit;
                    }
                    
                    %{$designExperimentRatio{"conditions"}{$condition}{$replicate}} = (
                        'anaName'  => $anaName,
                        'techRep'  => undef,
                        'fraction' => $fraction
                    ) if($replicate && !$designExperimentRatio{"conditions"}{$condition}{$replicate});
                }
                
                if($createProtAbundQuanti) {
                    $designExperimentAbund{"name"} = "Spectronaut quantification (abundance)" if(!$designExperimentAbund{"name"});
                    
                    %{$designExperimentAbund{"conditions"}{$anaName}{1}} = (
                        'anaName'  => $anaName,
                        'techRep'  => undef,
                        'fraction' => $fraction
                    ) if($replicate && !$designExperimentAbund{"conditions"}{$anaName}{1});
                }
    
                # Fragments Infos
                my $area = $infos[$colName2Index{'F.PeakArea'}];
                if ($area !~ /\d+\,?\d*/) {
                    $area = 'NA';
                } else {
                    $area =~ s/,/\./;
                }
                
                my $areaNorm = $infos[$colName2Index{'F.NormalizedPeakArea'}];
                if ($areaNorm !~ /\d+\,?\d*/) {
                    $areaNorm = 'NA';
                } else {
                    $areaNorm =~ s/,/\./;
                }
                
                my $fragmentType=$infos[$colName2Index{'F.FrgType'}].$infos[$colName2Index{'F.FrgNum'}];
                my $fragmentCharge=$infos[$colName2Index{'F.Charge'}];
                my $fragmentResidue=($infos[$colName2Index{'F.FrgLossType'}] && $infos[$colName2Index{'F.FrgLossType'}] ne 'noloss') ? $infos[$colName2Index{'F.FrgLossType'}] : '';
                my $fragUsedForPepQuanti = ($infos[$colName2Index{'F.ExcludedFromQuantification'}] eq 'False') ? 1 : 0;
                my $fragmentMZ=$infos[$colName2Index{'F.FrgMz'}];
                $fragmentMZ=~s/,/\./;
    
                my $fragment=$fragmentType.'_'.$fragmentCharge;
                $fragment .= "_$fragmentResidue" if($fragmentResidue);
                
                ### Store peptide infos
                @{$peptideList{$anaName}{$peptide}{"$RT"}} = ($pepMZ, $pepScore, $pepRawArea, $pepArea, $pepIsImputed);
                
                ### Store fragments infos
                @{$fragmentsInfos{$anaName}{$peptide}{"$RT"}{$fragment}}=('',$area,$fragmentMZ,'',$areaNorm) if($areaNorm ne '1'); # Do not store fragments if it has no/null intensity (-> F.NormalizedPeakArea = 1 in Spectronaut)
                
                ### Store peptide/protein associations
                if($peptide eq $lastPeptide && $peptideProt{$lastAna}{$peptide}) {
                    $peptideProt{$anaName}{$peptide}=\%{$peptideProt{$lastAna}{$peptide}};
                } else {
                    my @proteins = split(/;/,$proteinsList);
                    for my $i (0 .. $#proteins) {
                        my $protein = $proteins[$i];
                        $allProts{$protein} = 1;
                        next if($protein=~/reverse/);
                        @{$assosProtPeptide{$anaName}{$protein}{$peptide}}=($RT, undef, $pepIsImputed);
                        $peptideProt{$anaName}{$peptide}{$protein}=1;
                        
                        if($spectronautQuantiParams{"grouping_level"} eq 'protein') {
                            if($pepUsedForProtQuanti && ($spectronautQuantiParams{"proteotypicity"} ne 'unique' || ($spectronautQuantiParams{"proteotypicity"} eq 'unique' && $pepSpecificity < 100)) && ($proteinQuantity && $proteinQuantity ne '1')) { # Do not store protein quantity if it has no/null intensity (-> PG.Quantity = 1 in Spectronaut)
                                $quantiAbund{"quantity"}{$anaName}{$protein} = $proteinQuantity;
                                $quantiAbund{"allPep"}{$anaName}{$protein}{"$formattedSeq"} = 1;
                                $quantiAbund{"allPepTot"}{$protein}{"$anaName.$formattedSeq"} = 1;
                                delete $spectronautQuantiParams{"nbLostGroupAbund"}{$protein} if($spectronautQuantiParams{"nbLostGroupAbund"}{$protein});
                            } elsif(!$spectronautQuantiParams{"nbLostGroupAbund"}{$protein} && !$quantiAbund{"allPepTot"}{$protein}) {
                                $spectronautQuantiParams{"nbLostGroupAbund"}{$protein} = 1;
                            }
                        }
                    }
                }
                
                ###>recovering peptide modifications
                if (index($modifSeqStr, '[') != -1 || index($modifSeqStr, '(') != -1) {
                    my $fileFormat = (index($modifSeqStr, '[') != -1) ? "SPC" : "OS";
                    my $peptideSeq = $modifSeqStr;
                    my $allPep = 0;
                    my %pepQuantiInfo;
                    my %pepMods;
                    
                    while ($peptideSeq=~/(\w?)\[([^(\]]+)(?:\)?\s\([^)]+\))?\]|(\w?)\(UniMod:(\d+)\)/g) {
                        my $modIdentifier = ($2) ? $2 : ($4) ? $4 : '';
                        $modIdentifier = 1 if($fileFormat eq 'OS' && $modIdentifier == 5); # OpenSwath mismatch Acetylation with Carbamyl
                        my $residue = ($1) ? $1 : ($3) ? $3 : '=';
                        my $pos = ($1 || $3) ? $-[0]+1 : 0;
                        
                        if(!$modsInfosID{$modIdentifier}) {
                            my $errMessage = "Modification $modIdentifier in peptide $modifSeqStr is not valid";
                            $errMessage .= (@libIDs) ? " (not referenced in libraries)" : " (not found in db)";
                            printToFile($fileError, $errMessage, 1);
                            exit;
                        }
                            
                        my $modificationID = $modsInfosID{$modIdentifier};
                        if(!@libIDs && !$modifIDs{$modificationID}) {
                            $modifIDs{$modificationID}{'V'}{$residue}=1;
                        }
                        
                        $peptideModif{$anaName}{$peptide}{$modificationID}{$pos} = $massMods{$modificationID};
                        push(@{$pepMods{$modificationID}}, ($residue eq '=') ? $residue : "$residue$pos");
                        
                        $peptideSeq=~s/\[[^\]]+\]// if($fileFormat eq 'SPC');
                        $peptideSeq=~s/\(UniMod:\d+\)// if($fileFormat eq 'OS');
                    }
                    
                    if($spectronautQuantiParams{"grouping_level"} eq 'peptide') { # Do not store peptide quantity if it has no/null intensity (-> PEP.Quantity = 1 in Spectronaut)
                        my $vModStr = "";
                        foreach my $modID (sort { $a <=> $b } keys %pepMods) {
                            $vModStr .= "&" if($vModStr);
                            $vModStr .= "$modID#".join('.', @{$pepMods{$modID}});
                        }
                        
                        my $vModIndex = $peptideSeq."_".$vModStr;
                        if($pepUsedForPepGroupQuanti && ($pepQuantity && $pepQuantity ne '1')) {
                            foreach my $protein(keys %{$peptideProt{$anaName}{$peptide}}) {
                                next if($quantiAbund{"quantity"}{$anaName}{$protein}{$vModIndex});
                                $quantiAbund{"quantity"}{$anaName}{$protein}{$vModIndex} = $pepQuantity;
                                $quantiAbund{"allPep"}{$anaName}{$protein}{$vModIndex}++;
                                $quantiAbund{"allPepTot"}{$protein}{"$anaName.$vModIndex"} = 1;
                            }
                            $foundPep{$vModIndex}=1;
                            delete $spectronautQuantiParams{"nbLostGroupAbund"}{$formattedSeq} if($spectronautQuantiParams{"nbLostGroupAbund"}{$formattedSeq});
                        } elsif(!$spectronautQuantiParams{"nbLostGroupAbund"}{$formattedSeq} && !$foundPep{$vModIndex}) {
                            $spectronautQuantiParams{"nbLostGroupAbund"}{$formattedSeq} = 1;
                        }
                    }
                }
                
                ####>recovering peptide modifications score if provided
                if ($modifSeqStrScoreAll && index($modifSeqStrScoreAll, '[') != -1) {
                    my $peptideSeq = $modifSeqStrScoreAll;
                    while ($peptideSeq=~/\[([^(\]]+)(?:\)?\s\([^)]+\)):\s([0-9,.]+)%\]/g) {
                        my $pos=$-[0];
                        my ($vModName, $vModScore) = ($1, $2);
                        
                        if(!$modsInfosID{$vModName}) {
                            my $errMessage = "Modification $1 in peptide $modifSeqStrScoreAll is not valid";
                            $errMessage .= (@libIDs) ? " (not referenced in libraries)" : " (not found in db)";
                            printToFile($fileError, $errMessage, 1);
                            exit;
                        }
                        
                        my $modificationID = $modsInfosID{$vModName};
                        if($vModScore) {
                            $vModScore =~ s/,/\./;
                            $vModScore /= 100;
                            if(!$peptideModif{$anaName}{$peptide}{$modificationID}{"-1"} || index($peptideModif{$anaName}{$peptide}{$modificationID}{"-1"}, "$pos:$vModScore") == -1) {
                                $peptideModif{$anaName}{$peptide}{$modificationID}{"-1"} .= "," if($peptideModif{$anaName}{$peptide}{$modificationID}{"-1"});
                                $peptideModif{$anaName}{$peptide}{$modificationID}{"-1"} .= "$pos:$vModScore";
                            }
                        }
                        
                        $peptideSeq=~s/\[[^\]]+\]//;
                    }
                }
                
                if($iFragment % 100000 == 0) {
                    printToFile($fileStat, "Parsing input files : Extracting fragments information ($iFragment/$nFragments)", 1);
                }
                
                $lastPeptide = $peptide;
                $lastAna = $anaName;
                $iFragment++;
                
                undef $peptide;
                undef $anaName;
            }
        }
    }
}

sub extractCandidatesFileContent {
    my ($candidatesFile) = @_;
    
    if(-e $candidatesFile) {
        my (%colName2Index);
        my $protNotQuantified = 0;
        my $nCandidates = `wc -l $candidatesFile`;
        $nCandidates = ($nCandidates=~ /(\d+)/) ? $1 : 0;
        my $iCandidate = 0;
        
        open (IN, "<$candidatesFile");
        my $lastPeptide = '';
        my $lastAna;
        
        while (<IN>) {
            if ($iCandidate == 0) {
                my @columns = split(/[,;\t\r\n]/,$_);
                my @refColNames = ('Valid', 'Condition Numerator', 'Condition Denominator', 'AVG Group Quantity Denominator', 'AVG Group Quantity Numerator', 'Ratio', '% Change', 'Standard Deviation', 'Qvalue', 'ProteinGroups', '# Unique Total Peptides', '# Unique Total EG.Id');
                foreach my $i (0 .. $#columns) {
                    $colName2Index{$columns[$i]} = $i;
                }
    
                ##check columns names
                foreach my $colName (@refColNames) {
                    unless (defined $colName2Index{$colName}) {
                        printToFile($fileError, "Column '$colName' not found in candidates file !", 1);
                        exit;
                    }
                }
            } else {
                my @infos = split(/\t/,$_);
                $spectronautQuantiParams{"nbLostProtRatio"}++ if(!$infos[$colName2Index{'Valid'}]);
                next if(!$infos[$colName2Index{'Valid'}] || $infos[$colName2Index{'Valid'}] eq 'False');
                
                #stateQuantiRatio
                my $numerator = $infos[$colName2Index{'Condition Numerator'}];
                my $denominator = $infos[$colName2Index{'Condition Denominator'}];
                my $group = $infos[$colName2Index{'Group'}];
                my @proteins = split(/;/, $infos[$colName2Index{'ProteinGroups'}]);
                    
                if($spectronautQuantiParams{"grouping_level"} eq 'peptide') {
                    my %pepMods;
                    my $vModStr = "";
                    (my $groupStr = $group) =~ s/_//g;
                    
                    while($groupStr =~ /(\w?)\[([^(\]]+)(?:\)?\s\([^)]+\))?\]/g) {
                        my $vModName = $2;
                        my $residue = ($1) ? $1 : '=';
                        my $pos = ($1) ? $-[0]+1 : 0;
                        
                        unless($modsInfosID{$vModName}) {
                            printToFile($fileError, "Modification $2 in peptide $group is not valid (not referenced in libraries).", 1);
                            exit;
                        }
                    
                        my $modificationID = $modsInfosID{$vModName};
                        push(@{$pepMods{$modificationID}}, ($residue eq '=') ? $residue : "$residue$pos");
                        
                        $groupStr=~s/\[[^\]]+\]//;
                    }
                    
                    foreach my $modID (sort { $a <=> $b } keys %pepMods) {
                        $vModStr .= "&" if($vModStr);
                        $vModStr .= "$modID#".join('.', @{$pepMods{$modID}});
                    }
                    
                    $group = $groupStr."_".$vModStr;
                }
                
                # Replace , by . in numeric values and transform scientific ones to decimal
                foreach my $fieldName (('Standard Deviation', 'Qvalue', 'Ratio')) {
                    my $numValue = $infos[$colName2Index{$fieldName}];
                    $numValue =~ s/,/\./;
                    
                    if ($numValue =~ /[\d\.]+E(-?)(\d+)/) {
                        my $suffixNb = ($1) ? $2 : 2;
                        $numValue = sprintf("%.".$suffixNb."f", $numValue);
                    }
                    
                    $infos[$colName2Index{$fieldName}] = $numValue;
                }
                
                my %groupQuantiInfo = (
                    "RATIO" => $infos[$colName2Index{'Ratio'}],
                    "PVAL_ADJ" => $infos[$colName2Index{'Qvalue'}],
                    "SD_GEO" => exp($infos[$colName2Index{'Standard Deviation'}]),
                    "DIST_PEP_USED" => $infos[$colName2Index{'# Unique Total Peptides'}],
                    "NUM_PEP_USED" => $infos[$colName2Index{'# Unique Total EG.Id'}]
                );
                    
                foreach my $protein (@proteins) {
                    if($spectronautQuantiParams{"grouping_level"} eq 'peptide') {
                        %{$stateQuantiRatio{"$numerator/$denominator"}{$protein}{$group}} = %groupQuantiInfo;
                        $stateQuantiRatio{"NUM_PEP_TOTAL"}{$protein}{$group} = $stateQuantiRatio{"$numerator/$denominator"}{$protein}{$group}{"NUM_PEP_USED"}; # Global amount of peptide used by protein
                    } else {
                        %{$stateQuantiRatio{"$numerator/$denominator"}{$protein}} = %groupQuantiInfo;
                        $stateQuantiRatio{"NUM_PEP_TOTAL"}{$protein} = $stateQuantiRatio{"$numerator/$denominator"}{$protein}{"NUM_PEP_USED"}; # Global amount of peptide used by protein
                    }
                }
                
                if($iCandidate % 100000 == 0) {
                    printToFile($fileStat, "Parsing input files : Extracting candidates information ($iCandidate/$nCandidates)", 1);
                }
            }
            
            $iCandidate++;
        }
    }
}

sub parseLibraryInfo {
    my ($libsInfoRef) = @_;
    my %libsInfo = %{$libsInfoRef};
    for(my $i=0; $i<scalar @{$libsInfo{"ID"}}; $i++) {
        my ($libID, $libName, $libType, $libSplit) = (@{$libsInfo{"ID"}}[$i], @{$libsInfo{"name"}}[$i], @{$libsInfo{"software"}}[$i], @{$libsInfo{"split"}}[$i]);
        
        if (-e "$promsPath{swath_lib}/SwLib_$libID/$libName.sptxt") {
            print("Looking for $promsPath{swath_lib}/SwLib_$libID/$libName.sptxt<br\>\n");
            open(LIBFILE,"<","$promsPath{swath_lib}/SwLib_$libID/$libName.sptxt") or die ("open: $!");
            my $nameLine;
            my (%mod2Unimod,%pepModif);
            while (my $line=<LIBFILE>) {
                if($libType eq 'TPP') {
                    if ($line=~/^Name: (.+)/) { #new peptide line
                        ($nameLine=$1)=~s/\s|n//g;
                        if($nameLine=~/\[\d*\]/) {
                            while ($nameLine=~/(\w)*\[(\d*)\]/g) {
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
                                my $mass=sprintf("%.0f",$massMods{$modifID});
                                my $sign = ($mass > 0) ? '+' : '-';
                                $nameLine=~s/\[$modMass\]/\[$sign$mass\]/;
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
                                my $nbProt=($libSplit==0) ? $1 : ($2=~/Subgroup_[0-9]_([0-9]).*/) ? $1 : '';
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
                }
            }
            close LIBFILE;
        } elsif (-e "$promsPath{swath_lib}/SwLib_$libID/$libName.tsv") {
            print("Looking for $promsPath{swath_lib}/SwLib_$libID/$libName.tsv<br\>\n");
			open (INFILE,"<","$promsPath{swath_lib}/SwLib_$libID/$libName.tsv") or die("open: $!");

            # Parse header to retrieve specific columns
			my $line=<INFILE>;
            my @headers = split(/[\t]/, $line);
            my %colName2Index;
            foreach my $i (0 .. $#headers) {
                $colName2Index{$headers[$i]} = $i;
            }

            while (($line=<INFILE>)) {
                my @lineContent = split(/[\t]/, $line);
                my $pepModifSeq = $lineContent[$colName2Index{'IntLabeledPeptide'}];
                $pepModifSeq = substr($pepModifSeq, 1, -1);
                my $pepMZ = $lineContent[$colName2Index{'PrecursorMz'}];
                $pepMZ =~ s/,/\./;
                my $pepCharge = $lineContent[$colName2Index{'PrecursorCharge'}];
                my @pepProtGroups = split(/;/, $lineContent[$colName2Index{'ProteinGroups'}]);
                my $pepSpecificity = (1/scalar @pepProtGroups)*100;
                @{$peptideInfo{"$pepModifSeq\_$pepCharge"}} = ($pepMZ, 0, $pepSpecificity);
            }
            exit;
        } else {
            printToFile($fileError, "Library file not found.", 1);
            exit;
        }
    }
}

sub extractDBContent {
    my ($dbID, $dbRank, $proteinsRef, $dbRankProtRef) = @_;
    my $dbFasta = $dbh->selectrow_array("SELECT FASTA_FILE FROM DATABANK WHERE ID_DATABANK='$dbID'");

    if($dbFasta =~ /prd-mascot/) { # Fetch from mascot server using &promsMod::getProtInfo
        my (%protDes, %protMW, %protOrg, %protLength, %protSeq, %protList);
        &promsMod::getProtInfo('silent', $dbh, $dbID, [-1, "$workDir/db$dbID.fasta", $clusterInfo{'on'}], \%protDes, \%protMW, \%protOrg, \%protLength, \%protSeq, \%allProts);
        
        if(%protSeq) {
            foreach my $identifier (keys %protSeq) {
                my $finalIdentifier = $identifier;
                
                # Force conversion of IDENTIFIER to UniProt - ACC so it matches spectronaut protein IDs
                if($identifier=~/^(?:tr\|sp)|([^\|]+)\|/) {
                    $finalIdentifier = $1;
                }
                
                $proteinsRef->{$finalIdentifier}{"sequence"} = ($protSeq{$identifier}) ? $protSeq{$identifier} : '';
                $dbRankProtRef->{$finalIdentifier}=$dbRank unless defined $dbRankProt{$identifier};
            }
        } else {
            printToFile($fileError, "Fasta file from databank $dbFasta ($dbID) was not found", 1);
            exit;
        }
    } else { # Fetch local DB with custome search (all to UniProt - ACC to match spectronaut proteins ID format)
        my $dbDir = "$promsPath{data}/banks/db_$dbID";
        my ($identType)=$dbh->selectrow_array("SELECT DEF_IDENT_TYPE FROM DATABANK_TYPE DT,DATABANK D WHERE D.ID_DBTYPE=DT.ID_DBTYPE AND ID_DATABANK=$dbID");
        $identType=($identType) ? $identType : '';
        
        if (-e "$dbDir/$dbFasta") {
            my ($noCheck, $identifier, $matched, $sequence);
            open(FASTAFILE,"<","$dbDir/$dbFasta") or die ("open: $!");
            while (my $line=<FASTAFILE>) {
                if ($line=~/>reverse_/) { $noCheck=1;}  #next if reverse = true
                elsif ($line=~m/^>(.+)\n/) {       # prot description line
                    (my $newLine=$line)=~s/^>\s*//;
                    if ($identifier) {      # search prot
                        $sequence=~s/\s//;
                        
                        ## Store protein information
                        $proteinsRef->{$identifier}{"sequence"} = $sequence;
                        $dbRankProtRef->{$identifier}=$dbRank unless defined $dbRankProt{$identifier};
                        
                        $sequence = '';
                        $identifier = '';
                    }
                    
                    $noCheck=0;
                    
                    ## Fetch protein information
                    my @line=split(/\s/,$newLine);
                    $identifier = $line[0];
                    if($identType eq 'UNIPROT_ALL' && $identifier =~ /\|(.+)\|/) {
                        $identifier = $1;
                    }
                    chomp $identifier;
                    
                    unless($allProts{$identifier}) {
                        $identifier = '';
                        $sequence = '';
                        $noCheck = 1;
                    }
                } else {
                    if($noCheck == 0) {
                        chomp $line;
                        $line=~s/\s//g;
                        $sequence .= $line;
                    }
                }
            }
            
            ## Store last protein information
            if($identifier) {
                $proteinsRef->{$identifier}{"sequence"} = $sequence;
                $dbRankProtRef->{$identifier}=$dbRank unless defined $dbRankProt{$identifier};
            }
            
            close FASTAFILE;
        } else {
            printToFile($fileError, "Fasta file not found", 1);
            exit;
        }
    }
    
    my $nbProt = scalar keys %{$proteinsRef};
    printToFile($fileStat, "Databank : Found $nbProt proteins in the DB fasta file $dbFasta", 1);
}

sub matchPeptidesToProteins {
    my ($proteinsDBRef, $assosProtPeptideRef) = @_;
    my %protNotFound;
    foreach my $analysisNumber (keys %{$assosProtPeptideRef}){        # for each analysis
        foreach my $protein (keys %{$assosProtPeptideRef->{$analysisNumber}}){
            my %positionPeptide;
            my $sequence = ($proteinsDBRef->{$protein}{"sequence"}) ? $proteinsDBRef->{$protein}{"sequence"} : '';
            my $numMatch=0;
            my (@begPep,@endPep);
            
            if(!$sequence || $sequence eq '') {
                $protNotFound{$protein} = 1;
            }
            
            foreach my $peptide (keys %{$assosProtPeptideRef->{$analysisNumber}{$protein}}){
                my @pepSeq=split(/_/,$peptide);
                my $pepSequence=$pepSeq[0];
                $pepSequence=~s/\[[+-]?\w+\]//g;
                
                if(!$pepSequence) {
                    printToFile($fileStat, "Databank : Peptide not found for $protein => $peptide", 1);
                }
                
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
                else{ # for Leucine and Isoleucine inversion in peptide sequence and for peptides of protein not found in DB fasta file
                    delete $assosProtPeptideRef->{$analysisNumber}{$protein}{$peptide};
                }
            }
            $numMatchProt{$analysisNumber}{$protein}=$numMatch; #number of peptide match
            delete $assosProtPeptideRef->{$analysisNumber}{$protein} if(scalar keys %{$assosProtPeptideRef->{$analysisNumber}{$protein}} == 0);

            if ($assosProtPeptideRef->{$analysisNumber}{$protein}) {
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
                $protAnaCoverage{$analysisNumber}{$protein}=$peptideCoverage;
            }
        }
    }
    
    # Return proteins found in results but not in DBs
    return %protNotFound;
}

sub buildSpectronautQuantifAnnot {  
    my ($spectronautQuantiParamsRef, $designExperimentIDsRef, $quantiPepID, $quantiType) = @_;
    my $quantifAnnot = "LABEL=FREE::SOFTWARE=Spectronaut;".$spectronautQuantiParamsRef->{"version"};
	$quantifAnnot .= "::ALGO_TYPE=DIA";
    $quantifAnnot .= "::TOP_N=".$spectronautQuantiParamsRef->{"top_n"} if($quantiType eq 'abundance' && $spectronautQuantiParamsRef->{"top_n"});
    $quantifAnnot .= "::ABUND_MEASURES=".$spectronautQuantiParamsRef->{"intensity_metric"} if($quantiType eq 'abundance');
    
    # Normalization
    my $normalization = "none.none";
    my $biasCorrection = "FALSE";
    if($spectronautQuantiParamsRef->{"normalization_strategy"}) {
        $biasCorrection = "TRUE";
        if($spectronautQuantiParamsRef->{"normalization_strategy"} eq 'Global Normalization') {
            if($spectronautQuantiParamsRef->{"normalization_on"} eq 'Median') {
                $normalization = "median.none";
            } elsif($spectronautQuantiParamsRef->{"normalization_on"} eq 'Average') {
                $normalization = "mean.none";
            } elsif($spectronautQuantiParamsRef->{"normalization_on"} eq 'Geometric Mean') {
                $normalization = "mean.geom";
            }
        } elsif($spectronautQuantiParamsRef->{"normalization_strategy"} eq 'Local Normalization') {
            $normalization = "rt.local.reg";
        }
    }
    $quantifAnnot .= "::BIAS_CORRECTION=$biasCorrection";
    $quantifAnnot .= "::NORMALIZATION_METHOD=$normalization";
    
    # Peptides selection
	$quantifAnnot .= "::PEPTIDES=".$spectronautQuantiParams{"proteotypicity"}.";1;1;all;all";
	$quantifAnnot .= "::PTM_POS=SPC:".$spectronautQuantiParamsRef->{"probability_cutoff"}.";exclude" if($spectronautQuantiParams{"grouping_level"} eq 'peptide'); # TODO: Give the possibility to the user to select the PTM score threshold to use to import data
    
	$quantifAnnot .= "::SINGLE_REF=0" if($quantiType eq 'ratio');
    $quantifAnnot .= "::RATIO_TYPE=";
    $quantifAnnot .= ($quantiType eq 'ratio') ? 'SimpleRatio' : 'None';
	$quantifAnnot .= "::VISIBILITY=1";
    
	$quantifAnnot .= "::STATES=";
    if($designExperimentIDsRef) {
        my $statesStr = "";
        foreach my $stateID (sort {$designExperimentIDsRef->{"STATE"}{$a}{"INDEX"} <=> $designExperimentIDsRef->{"STATE"}{$b}{"INDEX"}} keys %{$designExperimentIDsRef->{"STATE"}}) {
            my $stateIndex = $designExperimentIDsRef->{"STATE"}{$stateID}{"INDEX"};
            
            $statesStr .= ";" if($statesStr);
            $statesStr .= scalar keys %{$designExperimentIDsRef->{"STATE"}{$stateID}{"OBS"}}; # Num observations for this state
            $statesStr .= ",";
            my $obsStr = "";
            foreach my $obsID (sort { $a <=> $b } keys %{$designExperimentIDsRef->{"STATE"}{$stateID}{"OBS"}}) {
                my $anaID = $designExperimentIDsRef->{"STATE"}{$stateID}{"OBS"}{$obsID};
                $obsStr .= "." if($obsStr);
                $obsStr .= "#$obsID:#$quantiPepID:#$anaID:0";
            }
            $statesStr .= $obsStr;
            $statesStr .= ",#$stateID";
        }
        $quantifAnnot .= $statesStr;
    }
    
    if($quantiType eq 'ratio') {
        my %ratioCombination;
        my $ratioStr = "";
        foreach my $stateID1 (sort {$designExperimentIDsRef->{"STATE"}{$a}{"INDEX"} <=> $designExperimentIDsRef->{"STATE"}{$b}{"INDEX"}} keys %{$designExperimentIDsRef->{"STATE"}}) {
            my $stateIndex1 = $designExperimentIDsRef->{"STATE"}{$stateID1}{"INDEX"};
            
            foreach my $stateID2 (sort {$designExperimentIDsRef->{"STATE"}{$a}{"INDEX"} <=> $designExperimentIDsRef->{"STATE"}{$b}{"INDEX"}} keys %{$designExperimentIDsRef->{"STATE"}}) {
                my $stateIndex2 = $designExperimentIDsRef->{"STATE"}{$stateID1}{"INDEX"};
                
                if($stateID1 != $stateID2 && !($ratioCombination{"#$stateID1/#$stateID2"} || $ratioCombination{"#$stateID2/#$stateID1"})) {
                    $ratioCombination{"#$stateID2/#$stateID1"} = 1;
                    $ratioStr .= ";" if($ratioStr);
                    $ratioStr .= "#$stateID2/#$stateID1";
                }
            }
        }
        
        $quantifAnnot .= "::RATIOS=$ratioStr";
    }
    
    
    my $lostProteins = 0;
    if($quantiType eq 'abundance') {
        if($spectronautQuantiParams{"grouping_level"} eq 'peptide') {
            foreach my $grpSeq (keys %{$spectronautQuantiParamsRef->{"nbLostGroupAbund"}}) {
                $lostProteins += $spectronautQuantiParamsRef->{"nbLostGroupAbund"}{$grpSeq};
            }
        } elsif($spectronautQuantiParams{"grouping_level"} eq 'protein' && $spectronautQuantiParamsRef->{"nbLostGroupAbund"}) {
            $lostProteins = scalar keys %{$spectronautQuantiParamsRef->{"nbLostGroupAbund"}};
        }
    } elsif($quantiType eq 'ratio' && $spectronautQuantiParamsRef->{"nbLostProtRatio"}) {
        $lostProteins = $spectronautQuantiParamsRef->{"nbLostProtRatio"};
    }
	$quantifAnnot .= "::LOST_PROTEINS=$lostProteins";
}

####>Revision history<####
# 1.1.6 [CHANGE] Changed RT format from hours to minutes (VS 02/04/2021)
# 1.1.5 [CHANGE] Adapted import job's memory and time limit thresholds (VS 04/03/2021) 
# 1.1.4 [BUGFIX] Added forgotten dbConnect after preventive disconnection to avoid DB connection timeout (PP 18/02/21)
# 1.1.3 [BUGFIX] Disconnect from db during long processing steps to avoid DB timeout (VS 17/02/21)
# 1.1.2 [ENHANCEMENT] Async copy/move of uploaded/shared directory files (VS 16/02/21)
# 1.1.1 [BUGFIX] Fixed: DirectDIA associated all modifications to analysis (VS 02/02/21)
# 1.1.0 [CHANGE] Display total peptides instead of distinct peptides in abundance quantifications (VS 29/01/21) 
# 1.0.13 [BUGFIX] Fixed: handles large datasets (VS 29/01/21) 
# 1.0.12 [MINOR] Added traceability of Top N parameter used in Spectronaut identification process (VS 29/01/21)
# 1.0.11 [BUGFIX] Fixed: properly count amount of peptides by protein in quantification results (VS 29/01/21)
# 1.0.10 [BUGFIX] Fixed import using library: LABEL in quantif annot should be the first field (VS 11/12/2020)
# 1.0.9 [ENHANCEMENT] Handles imputed data (VS 26/11/20)
# 1.0.8 [BUGFIX] Added links between all states and its related ratio/abundance quantification (VS 03/11/20)
# 1.0.7 [ENHANCEMENT] Possibility to select a databank from mascot server when using Direct DIA (VS 27/10/20)
# 1.0.6 [ENHANCEMENT] Compatibility with DirectDIA (VS 22/10/2020)
# 1.0.5 [ENHANCEMENT] Store lost proteins/peptides amount in quantification (VS 10/09/20)
# 1.0.4 [CHANGE] Do not take into account peptides/proteins with no intensity in abundance quantification (VS 09/09/20)
# 1.0.3 [UPDATE] Changed JOB_HISTORY.STATUS to JOB_HISTORY.JOB_STATUS (PP 28/08/20)
# 1.0.2 [FEATURE] Create all-vs-all protein ratio quantification based on candidates file (VS 01/09/20)
# 1.0.1 [ENHANCEMENT] Handles input files from shared directory (VS 31/08/20)
# 1.0.0 [FEATURE] First version of Spectronaut data import (VS 24/08/20)
