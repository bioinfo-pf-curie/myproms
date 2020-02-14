#!/usr/local/bin/perl -w

################################################################################
# importSwathDataRefonte.cgi         1.4.4                                     #
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
use XML::Simple;
use POSIX qw(strftime); # to get the time
use File::Path qw(rmtree); # remove_tree
use File::Spec::Functions qw(splitpath); # Core module
use File::Copy;
use File::Basename;
use String::Util qw(trim);
use List::Util qw(min max);

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
my $userID = (param('USERID')) ? param('USERID') : ($ENV{REMOTE_USER}) ? $ENV{REMOTE_USER} : '';
my $projectID = param('id_project');
my $branchID = param('ID');
my $libID = (param('libID')) ? param('libID') : '';
my ($item, $experimentID) = split(':', $branchID);
$experimentID = (param('expID')) ? param('expID') : $experimentID;
my $dbh = &promsConfig::dbConnect($userID);

my $exportedLibFile = (param('exportedlibrary')) ? param('exportedlibrary') : '';
my $paramFile = (param('exportparam')) ? param('exportparam') : '';
my $decoyLibFile = (param('decoyLibrary')) ? param('decoyLibrary') : '';
my $alignmentFile = (param('uploadfile')) ? param('uploadfile') : ''; # TRIC alignment results as TSV file
my $alignmentFileName;

my $tricMethod = (param('TRICMethod')) ? param('TRICMethod') : 'best_overall';
my $enabledIPF = (param('enableIPF')) ? param('enableIPF') : 0;
my $softwareVersion = (param('softversion')) ? param('softversion') : '';
my $format = (param('FORMAT')) ? param('FORMAT') : 'openswath'; # peakview, openswath, spectronaut

my $isOnCluster = (param('CLUSTER')) ? param('CLUSTER') : 0;
my $shouldRunOnCluster = (param('RUNCLUSTER')) ? param('RUNCLUSTER') : 0;
my $sameParams = (param('sameParams')) ? 1 : 0;

my $submit = (param('submit')) ? param('submit') : "";
my $fct = (param('FCT')) ? param('FCT') : "" ;

print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);

my $action = (param('ACT')) ? param('ACT') : "" ; ## import or quantification
if ($action eq 'excludeLibraryModifications') {
    &excludeLibraryModifications; exit;
} elsif ($action eq "selectSamples") {
    &selectSamples; exit;
}


#####################
### Starting HTML ###
#####################

## SHOULD DISPLAY FORMS
if ($submit eq "") {
    printFormContent();

## SHOULD PROCESS DATA
} else {

    #############################
    ###  Fetching parameters  ###
    #############################
    my %massATave = &promsConfig::getMassATave; # Average mass, needed for protein mass calculation
    my %massAAave = &promsConfig::getMassAAave; # Average mass, needed for protein mass calculation

    my $software = ($format eq 'peakview') ? "Peakview" : ($format eq 'openswath') ? "OpenSwath" : "Spectronaut" ;
    my (@mzXmlFileList, @pyprophetFiles, $irtFile, $swathWindows, $unimod, $w, $m, $fasta, $labelling); # Useful variables -> files: mzXML, iRT, Windows, Unimod file, -W, -M, -f, -i
    my $unimodRef = "/bioinfo/projects_prod/myproms/tmp/Swath/unimod_ref.xml"; # TODO Update file with all unimod modifications
    my ($lMax, $lMin, $s, $x, $o, $n); # Useful variables for both library export and conversion
    my $errorTRICMethod = 0;
    
    ### Create new directory
    my $time = strftime("%Y%m%d%H%M%S", localtime); # "20190808181244";
    my $libDir = "$promsPath{data}/swath_lib/SwLib_$libID";
    my $workDir = (param('workdir')) ? param('workdir') : "$promsPath{data}/tmp/Swath/$time";
    $CGITempFile::TMPDIRECTORY = $workDir; # folder where GCI temporary stores files
    mkdir $workDir unless -e $workDir;

    ### Fetching parameters
    my $mzThreshold = param('mzThreshold');
    my $dscore = param('dscore');
    my $outFile = $workDir.'/openSwathOUT.txt';
    
    my $fileStat = "$workDir/status_$libID.out";
	my $fileError = "$workDir/status_$libID\_error.out";
    my $fileLog = "$workDir/log.out";
    
    ### Print waiting screen
    print qq |
        <HTML>
            <HEAD>
                <LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
                <BODY background="$promsPath{images}/bgProMS.gif">
                <CENTER><BR/><BR/>
                    <IMG src='$promsPath{images}/engrenage.gif'><BR/><BR/>
                    <FONT color=\"red\"><DIV id=\"waitDIV\"></DIV></FONT><BR/>
                    <B>Processing</B><BR/><BR/>
    |;
    
    
    #############################
    ###  Upload input files   ###
    #############################
    
    # Import quantification mode
    if($action eq 'quantification') {
        ### Upload mzXML data file(s)
        if(param('mzXMLFiles')) {
            push @mzXmlFileList, uploadFiles($workDir, 'mzXML');
        }
        
        ### Upload mzXML data file(s) from shared directory
        if(param('sharedDirFiles')) {
            push @mzXmlFileList, uploadFiles($workDir, 'sharedDir');
        }
        
        # Recover pyprophet files if any
        if(param('pyprophetFiles')) {
            push @pyprophetFiles, uploadFiles($workDir, 'pyprophet');
        }
        
        ### Check if it should add other analysis results (use other mzXML files in the workflow)
        if(param('samples')) {
            printToFile($fileStat, "Recovering results files to merge from previous experiment(s)", 1);
    
            ### Recovering mzXML / PyProphet files from other experiments
            my $sthAnaID=$dbh->prepare("SELECT ID_ANALYSIS,NAME FROM ANALYSIS WHERE ID_SAMPLE=?");
            foreach my $sampleID (param('samples')) {
                $sthAnaID->execute($sampleID);
                while (my ($anaID, $anaName) = $sthAnaID->fetchrow_array) {
                    my $anaDir = "$promsPath{peptide}/proj_$projectID/ana_$anaID";
                    my $mzXMLFile = "$anaName.mzXML";
                    my $pyprophetFile = "$mzXMLFile.tsv";
                    
                    if($sameParams) {
                        unless(-s "$workDir/$pyprophetFile") {
                            if (-s "$anaDir/$pyprophetFile") {
                                print("Recovering $pyprophetFile from $anaDir <br/>\n");
                                copy("$anaDir/$pyprophetFile", "$workDir/$pyprophetFile");
                            } elsif (-s "$anaDir/$pyprophetFile.tar.gz") {
                                print("Recovering $pyprophetFile.tar.gz from $anaDir<br/>\n");
                                system "cd $workDir; tar -zxf $anaDir/$pyprophetFile.tar.gz -C $workDir/;";
                            } else {
                                printToFile($fileError, "File $pyprophetFile not found for analysis $anaName ($anaID).", 1);
                                exit;
                            }
                        }
                        
                        push @pyprophetFiles, $pyprophetFile;
                    } else {
                        unless(-s "$workDir/$mzXMLFile") {
                            if (-s "$anaDir/$mzXMLFile") {
                                print("Recovering $mzXMLFile from $anaDir <br/>\n");
                                copy("$anaDir/$mzXMLFile", "$workDir/$mzXMLFile");
                            } elsif (-s "$anaDir/$mzXMLFile.tar.gz") {
                                print("Recovering $mzXMLFile.tar.gz from $mzXMLFile<br/>\n");
                                system "cd $workDir; tar -zxf $anaDir/$mzXMLFile.tar.gz -C $workDir/;";
                            } else {
                                printToFile($fileError, "File $mzXMLFile not found for analysis $anaName ($anaID).", 1);
                                exit;
                            }
                        }
                        
                        push @mzXmlFileList, $mzXMLFile;
                    }
                }
            }
            $sthAnaID->finish;
        }
        
        ### Upload iRT file
        $irtFile = param('iRTFile');
        if($irtFile) {
            uploadFile($workDir, $irtFile);
        }
        
        ### Upload export library file
        if ($exportedLibFile) {
            uploadFile($workDir, $exportedLibFile);
        }
        
        ### Upload decoy library
        if($decoyLibFile) {
            uploadFile($workDir, $decoyLibFile);
        }
    
        ### Upload OpenSwath Windows file
        $swathWindows = param('windowsFile');
        if($swathWindows) {
            uploadFile($workDir, $swathWindows);
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
        
        ### Upload unimod file
        $unimod = param("unimodFile");
        if($unimod) {
            uploadFile($workDir, $unimod);
            $enabledIPF = 1;
        }
        
        ### Labelling file
        $labelling = param("i");
        if($labelling) {
            uploadFile($workDir, $labelling);
        }
        
        ### Fasta file
        $fasta = param("f");
        if($fasta) {
            uploadFile($workDir, $fasta);
        }

    
    # Import OS Data mode
    } else {
        ### Upload alignment result file
        $alignmentFile = (param('uploadfile')) ? param('uploadfile') : '';
        if ($alignmentFile) {
            uploadFile($workDir, $alignmentFile);
            $alignmentFileName = fileparse("$workDir/$alignmentFile", qr/\.[^.]*/);
        }
        
        ### Upload export library file
        if ($exportedLibFile) {
            uploadFile($workDir, $exportedLibFile);
        }
        
        ### Upload export library param file
        if ($paramFile) {
            uploadFile($workDir, $paramFile);
        }
        
        ### Upload unimod file
        $unimod = param("unimodFile");
        if($unimod) {
            uploadFile($workDir, $unimod);
            $enabledIPF = 1;
        }
    }
    

    ###########################################
    ###  OPTION : Run full script on server ###
    ###########################################
    if ($shouldRunOnCluster && $clusterInfo{'on'}) {
        my $currentDirectory = dirname $0;
        my $scriptFile = "$workDir/openswathStoring.sh";
        
        # Build query options based on selected action
        my $queryOptions = "workdir=$workDir";
        my @paramNames = param;
        foreach my $paramName (@paramNames) {
            next if ($paramName eq "RUNCLUSTER" || $paramName eq "samples");
            next if($paramName eq 'mzXMLFiles' || $paramName eq 'sharedDirFiles');
            
            $queryOptions .= " $paramName=".param($paramName);
        }
        
        if($enabledIPF && !param("enableIPF")) {
            $queryOptions .= " enableIPF=1 ";
        }
        
        # Add mzXML Files
        if(@mzXmlFileList) {
            my $mzXMLOptionsStr = "";
            foreach my $mzXMLFile (@mzXmlFileList) {
                # Use custom delimiter OTHERFILE to separate each mzXML files provided (no space allowed)
                $mzXMLOptionsStr .= "OTHERFILE" if($mzXMLOptionsStr); 
                $mzXMLFile =~ s/\s/_/g;
                $mzXMLOptionsStr .= "$mzXMLFile";
            }
            $queryOptions .= " mzXMLFiles=$mzXMLOptionsStr";
        }
        
        # Add pyprophet files
        if(@pyprophetFiles) {
            my $pyprophetOptionsStr = "";
            foreach my $pyprophetFile (@pyprophetFiles) {
                # Use custom delimiter OTHERFILE to separate each mzXML files provided (no space allowed)
                $pyprophetOptionsStr .= "OTHERFILE" if($pyprophetOptionsStr); 
                $pyprophetFile =~ s/\s/_/g;
                $pyprophetOptionsStr .= "$pyprophetFile";
            }
            $queryOptions .= " pyprophetFiles=$pyprophetOptionsStr";
        }
        
        open (BASH,"+>", $scriptFile);
        print BASH "#!/bin/bash\n";
        print BASH "export LC_ALL='C'\n";
        print BASH "cd $currentDirectory\n";
        print BASH "echo \"$clusterInfo{path}{perl}/perl $0 CLUSTER=1 $queryOptions >> $fileLog 2>&1\" > $workDir/query.txt\n";
        print BASH "$clusterInfo{path}{perl}/perl $0 $queryOptions > $workDir/log.txt 2>&1\n";
        close BASH;
        chmod 0775, $scriptFile;
        
        my $maxMem;
        if($action eq 'quantification') {
            my $nbMzXMLFiles = (@mzXmlFileList) ? scalar @mzXmlFileList : 1;
            my $nbPyprophetFiles = (@pyprophetFiles) ? scalar @pyprophetFiles : 0;
            $maxMem = (($nbMzXMLFiles+$nbPyprophetFiles)/2 > 8) ? ($nbMzXMLFiles+$nbPyprophetFiles)/2 : ($enabledIPF) ? 30 : 8;
            $maxMem = 120 if($maxMem > 120);
        } else {
            my $size = (-e "$workDir/$alignmentFile") ? `stat -c "%s" $workDir/$alignmentFile`: 1073741824;
            $maxMem = max(20, ($size / 1073741824) * 8);
            $maxMem = 80 if($maxMem > 80);
        }
        
        my %jobParams = (
            maxMem => sprintf("%.0f", $maxMem).'Gb',
            numCPUs => 1,
            maxHours => ($action eq 'quantification') ? 120 : 8, # TODO Change to 240 hours when it is working properly
            jobName => $action."Swath_2_$time",
            outFile => 'PBSimport.txt',
            errorFile => 'PBSimportError.txt',
            jobEndFlag => "_END_".$action."Swath_$time",
            noWatch => '1',
        );
        my ($pbsError, $pbsErrorFile, $jobClusterID) = $clusterInfo{'runJob'}->($workDir, $scriptFile, \%jobParams);
        
        # Create new job to monitor
        $dbh->do("INSERT INTO JOB_HISTORY (ID_JOB, ID_USER, ID_PROJECT, ID_JOB_CLUSTER, TYPE, STATUS, FEATURES, SRC_PATH, LOG_PATH, ERROR_PATH, STARTED) VALUES('$time', '$userID', $projectID, 'C$jobClusterID', 'Import [$software]', 'Queued', 'SOFTWARE=$software;ID_LIBRARY=$libID', '$workDir', '$fileStat', '$fileError', NOW())");
        $dbh->commit;
        $dbh->disconnect;
        
        # Redirect to results
        print qq |
                        <SCRIPT LANGUAGE="JavaScript">
                            var monitorJobsWin=window.open("$promsPath{cgi}/monitorJobs.cgi?filterType=Import [$software]&filterDateNumber=3&filterDateType=DAY&filterStatus=Queued&filterStatus=Running",'monitorJobsWindow','width=1200,height=500,scrollbars=yes,resizable=yes');
                            monitorJobsWin.focus();
                            top.promsFrame.selectedAction='summary';
                            parent.navFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=experiment:$experimentID&ACT=nav";
                        </SCRIPT>
                    </CENTER>
                </BODY>
            </HTML>
        |;
        exit;
    } elsif(!param('workdir')) { # If script was not started by the cluster job, create it in DB
        $dbh->do("INSERT INTO JOB_HISTORY (ID_JOB, ID_USER, ID_PROJECT, ID_JOB_CLUSTER, TYPE, STATUS, FEATURES, SRC_PATH, LOG_PATH, ERROR_PATH, STARTED) VALUES('$time', '$userID', $projectID, 'L$$', 'Import [$software]', 'Queued', 'SOFTWARE=$software;ID_LIBRARY=$libID', '$workDir', '$fileStat', '$fileError', NOW())");
        $dbh->commit;
    }
    

    ########################
    ###  SWATH Workflow  ###
    ########################
    if ($action eq 'quantification') {

        ##########################
        ###  Initiate Workflow ###
        ##########################
        my ($libName, $libVersion) = $dbh->selectrow_array("SELECT NAME, VERSION_NAME FROM SWATH_LIB WHERE ID_SWATH_LIB=$libID");
        my %libraryParams = (
            "ID"      => $libID,
            "name"    => $libName,
            "dir"     => $libDir,
            "version" => $libVersion,
        );
    
        DIAWorkflow::new($workDir, $fileStat, \%promsPath, \%libraryParams, \@mzXmlFileList, $isOnCluster);
    
        ########################
        ###  Library Export  ###
        ########################
        if (!$exportedLibFile) { # Export selected library
            printToFile($fileStat, "Library conversion for OpenSwath", 1);
            
            ### Recover all exportation params
            # Excluded library modifications
            my $excludedPepMod = (param('pepMod')) ? param('pepMod') : '';
            
            # Mass range (lmin/lmax), Ion series (s) and charge (x), Number of ions per peptide (min: o / max: n)
            ($lMax, $lMin, $s, $x, $o, $n) = &promsMod::cleanParameters(param('-lmax'), param('-lmin'), param('-s'), param('-x'), param('-o'), param('-n'));
            
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
                "w" => $w, "format" => $format, "others" => $enabledIPF,
            );
            
            # TODO Uncomment
            if($enabledIPF) { # Only export param files for further database insertion, since the existing library mrm file will be used as input for the library conversion step
                $exportedLibFile = "$libName.mrm";
                $paramFile = $libName."_param";
                DIAWorkflow::exportParams(\%exportParams, "$workDir/$paramFile", "libExport");
                
                if (-e "$libDir/$exportedLibFile") {
                    copy("$libDir/$exportedLibFile", "$workDir/$exportedLibFile");
                } else {
                    printToFile($fileError, "ERROR: Required MRM library file does not exist. Run swath lib generation in Unsplit mode to create it.", 1);
                    exit;
                }
            } else { # Export library based on given params
                ($exportedLibFile, $paramFile, my $error) = DIAWorkflow::exportLibrary(\%exportParams, $excludedPepMod, $isOnCluster);
                if ($error) {
                    printToFile($fileError, "ERROR: Library conversion failed : $error", 1);
                    exit;
                }
                #$exportedLibFile = "$libName\_openswath.tsv"; # TODO Remove
                #$paramFile = "$libName\_param"; # TODO Remove
            }
            
            print("Exported library : $exportedLibFile<br/>"); # TODO Uncomment
        } 
    
        if(!$decoyLibFile) {
            ##########################################
            ###  library processing for openswath  ###
            ##########################################
            my ($productMZThresh, $precursorMZThresh, $precursorLowerMZ, $precursorUpperMZ) = &promsMod::cleanParameters(param('productMZThresh'), param('precursorMZThresh'), param('precursorLowerMZ'), param('precursorUpperMZ'));
            my %convertParams = (
                "lmin" => $lMin, "lmax" => $lMax, "s" => $s, "x" => $x, "o" => $o, "n" => $n,
                "productMZThresh" => $productMZThresh, "precursorMZThresh" => $precursorMZThresh, "precursorLowerMZ" => $precursorLowerMZ, "precursorUpperMZ" => $precursorUpperMZ,
                "w" => $swathWindows, "IPF" => $enabledIPF, "libFile" => $exportedLibFile, 
            );
            
            if($enabledIPF) {
                $convertParams{'maxNumLocalization'} = param("maxNumLocalization") if(param("maxNumLocalization"));
                $convertParams{'unimodFile'} = $unimod if($unimod);
            } else {
                $convertParams{'enableSpLosses'} = param('enableSpLosses') if(param('enableSpLosses'));
            }
            $convertParams{'enableUnspLosses'} = param('enableUnspLosses') if(param('enableUnspLosses'));
            
            
            printToFile($fileStat, "Converting DIA library", 1);
            
            # TODO Uncomment
            my $libConversionStatus = DIAWorkflow::convertLibrary(\%convertParams, $isOnCluster);
            if($libConversionStatus ne 'Done') {
                printToFile($fileError, $libConversionStatus, 1);
                exit;
            }
            
            $decoyLibFile = $libName.'_assay_DECOY.pqp';
        }
        
        if($enabledIPF) { # Phospho : Switch to a parsable tsv file for exported library, instead of a mrm one
            $exportedLibFile =~ s/\.mrm/\.tsv/;
        } else {
            $exportedLibFile = "$libName.tsv";
        }
        

        ##########################################################
        ###  openswath workflow and pyprophet on search files  ###
        ##########################################################
        printToFile($fileStat, "Running OpenSwath protein identification Workflow on mzXML files", 1);
        
        # OpensSwath options
        my ($mzExtractionWin, $ppm, $mzCorrectionFunction, $rtExtractionWin, $minPeakWidth) = &promsMod::cleanParameters(param('mzExtractionWin'), param('ppm'), param('mzCorrectionFunction'), param('rtExtractionWin'), param('minPeakWidth'));
        my %openSwathParams = (
            "decoyLibrary" => "$workDir/$decoyLibFile", "irtFile" => "$workDir/$irtFile", "swathWindows" => "$workDir/$swathWindows", "x" => $x, "o" => $o, "n" => $n,
            "mzExtractionWin" => $mzExtractionWin, "ppm" => $ppm, "mzCorrectionFunction" => $mzCorrectionFunction, "rtExtractionWin" => $rtExtractionWin, "minPeakWidth" => $minPeakWidth,
            "IPF" => $enabledIPF,
        );
         
        $openSwathParams{"enableUisScoring"} = param("enableUisScoring") if($enabledIPF || param("enableUisScoring"));
        $openSwathParams{"useMs1Traces"} = param("useMs1Traces") if($enabledIPF || param("useMs1Traces"));
        
        # TODO Uncomment
        my $openSwathStatus = DIAWorkflow::openswath(\%openSwathParams);
        if($openSwathStatus ne 'Done') {
            printToFile($fileError, $openSwathStatus, 1);
            exit;
        }
        

        ###################################
        ###  PyProphet on on OSW files  ###
        ###################################
        ## Add files name to merge in pyprophet analysis if any
        printToFile($fileStat, "Running error estimation on OpenSwath results", 1);
        
        # PyProphet options
        my $pyProphetMethod = (param('pyprophetMethod')) ? &promsMod::cleanParameters(param('pyprophetMethod')) : "experiment-wide";
        my %pyprophetParams = (
            "IPF" => $enabledIPF, "method" => $pyProphetMethod,
        );
        
        # TODO Uncomment
        my $pyprophetStatus = DIAWorkflow::pyprophet(\%pyprophetParams, 1); # TODO Uncomment
        my $nTries = 1;
        my $maxTries = 5;
        print("PyProphet status : $pyprophetStatus<br/>");
        while($pyprophetStatus ne 'Done' && $nTries < $maxTries) {
            $openSwathStatus = DIAWorkflow::openswath(\%openSwathParams);
            if($openSwathStatus ne 'Done') {
                printToFile($fileError, $openSwathStatus, 1);
                exit;
            }
            
            DIAWorkflow::pyprophet(\%pyprophetParams, 1);
            $nTries++;
        }
        
        # Compute pyprophet outfiles' name
        foreach my $mzXMLFile (@mzXmlFileList) {
            push @pyprophetFiles, "$mzXMLFile.tsv";
        }
        
        
        #################################
        ###  TRIC on pyprophet files  ###
        #################################
        printToFile($fileStat, "Running multi-alignment of peaks", 1);
        
        ### Recover TRIC options
        my ($fdrCutoff, $maxRTDiff, $maxFDRQuality, $dscoreCutoff) = &promsMod::cleanParameters(param('fdrCutoff'), param('maxRTDiff'), param('maxFDRQuality'), param('dscoreCutoff'));
        my %tricParams = (
            "fdrCutoff" => $fdrCutoff, "maxRTDiff" => $maxRTDiff, "maxFDRQuality" => $maxFDRQuality, "dscoreCutoff" => $dscoreCutoff,
            "format" => $format, "method" => $tricMethod, "inputFiles" => \@pyprophetFiles,
        );
        
        $alignmentFile = "feature_alignment_fileout.tsv";
        unless(-s "$workDir/$alignmentFile") {
            my $tricStatus = DIAWorkflow::tric(\%tricParams, $isOnCluster);
            print("TRIC status : $tricStatus<br/>");
            
            if($tricStatus ne 'Done') {
                $tricParams{"method"} = ($tricMethod ne 'LocalMST') ? 'LocalMST' : 'best_overall';
                
                printToFile($fileStat, "WARNING : TRIC method ($tricMethod) did not work on data", 1);
                printToFile($fileStat, "Running another method : $tricParams{method}", 1);
                $tricStatus = DIAWorkflow::tric(\%tricParams, $isOnCluster);
                
                print("TRIC status rerun : $tricStatus<br/>");
                if($tricStatus ne 'Done') {
                    printToFile($fileError, $tricStatus, 1);
                    exit;
                }
            }
        
            ### Check if the final file exists
            unless(-s "$workDir/$alignmentFile") {
                printToFile($fileError, "TRIC: Alignment of chromatograms failed.", 1);
                exit;
            }
        }
        
        $dbh->disconnect;
    }
    
    ### Transform PTMs text to unimod format : (Phospho) -> (UniMod:21) / (Acetyl) -> (UniMod:1)
    if ($enabledIPF) {
        castToUniMod("$workDir/$alignmentFile", ($unimod) ? "$workDir/$unimod" : $unimodRef);
    }

    #########################################################
    ###  Fetching library and modifications information   ###
    #########################################################
    printToFile($fileStat, "Scanning library : Fetching modifications information", 1);
    my (%peptideModif, %massMods, %unimod2ID, %peakview2modID, %modifIDs, %peptideInfo, %fragInfos, %modificationList);
    my ($fdr,$versionName,$taxonomy,$instrument,$quantifAnnot);
    my $dbh = &promsConfig::dbConnect($userID); # Reconnect after long workflow process
    
    ### Fetching SWATH LIB modifications
    (my $libName,$instrument,$taxonomy,$versionName,my $split,my $stat) = $dbh->selectrow_array("SELECT NAME,INSTRUMENT,ORGANISM,VERSION_NAME,SPLIT,STATISTICS FROM SWATH_LIB WHERE ID_SWATH_LIB='$libID' ") or die "Couldn't prepare statement: " . $dbh->errstr;
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
    
    ### Fetching SWATH LIB statistics
    my $xml=XML::Simple->new(KeepRoot=>1);
    my $xmlStat=$xml->XMLin($stat);
    my $pepNum=$xmlStat->{'NUM_ENTRY'}->{'NUM_PEP'};
    my $protNum=$xmlStat->{'NUM_ENTRY'}->{'NUM_PROT'};
    my $protSpeNum=$xmlStat->{'NUM_ENTRY'}->{'NUM_PROT_SPECIFIC'};
    if($format eq 'peakview'){
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
    open(PARAM, "<$workDir/$paramFile");
    while (<PARAM>) {
        next if($_ =~ /(Desired proteins not found into the library|Export for)/);
        my ($paramName, $param) = split(/:/ ,$_);
        #$param =~ s/\s//g;
        $quantifAnnot .= "::$paramName=$param";
    }
    close PARAM;
    
    ##########################################################
    ###  Recovering peptide miss cleavage and specificity  ###
    ##########################################################
    my $count=0;
    printToFile($fileStat, "Scanning library : Recovering peptide miss cleavage(s) and specificity", 1);
    if (-e "$promsPath{swath_lib}/SwLib_$libID/$libName.sptxt") {
        open(LIBFILE,"<","$promsPath{swath_lib}/SwLib_$libID/$libName.sptxt") or die ("open: $!");
        my $nameLine;
        my (%mod2Unimod,%pepModif);
        while (my $line=<LIBFILE>){
            if ($line=~/^Name: (.+)/) { #new peptide line
                ($nameLine=$1)=~s/\s|n//g;
                if($nameLine=~/\[\d*\]/ && $format eq 'peakview'){   #peptide mofification
                    while($nameLine=~/(\w)\[(\d*)\]/g){   #conversion of library modification code into peakview modification code (for peptide matches)
                        my $modMass=$2;
                        my $modifID;
                        my $aa=($1) ? $1 : '';
                        my $pos=$-[0];
                        foreach my $ID (keys %massMods) {
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
        }
        close LIBFILE;
    } else{
        printToFile($fileError, "Library file not found.", 1);
        exit;
    }
    $dbh->disconnect;
    
    ## Parse export library for fragments information
    if ($format eq 'openswath' && $exportedLibFile) {
        printToFile($fileStat, "Scanning library : Recovering fragments information in the library", 1);
        if (-s "$workDir/$exportedLibFile") {
            my (%colName2Index, %modifMass);
            my $nFrag = `wc -l $workDir/$exportedLibFile | cut -d' ' -f1`;
            $nFrag = ($nFrag) ? $nFrag-1 : 0;
            my $fragCount = 0;
            
            open (LIBEXPORT, "<$workDir/$exportedLibFile");
            while (<LIBEXPORT>) {
                if ($. == 1) {
                    $_ =~ s/\s*\Z//;
                    my @columns = split(/[,;\t]/, $_);
                    foreach my $i (0 .. $#columns) {
                        $colName2Index{$columns[$i]} = $i;
                    }
                } else {
                    if($fragCount % 100000 == 0) {
                        printToFile($fileStat, "Scanning library : Recovering fragments informations in the library ($fragCount/$nFrag)", 1);
                    }
                    $fragCount++;
                    
                    $_ =~ s/\s*\Z//;
                    next if($_ =~ /DECOY/);
                    my @infos = split(/[\t]/, $_);

                    my $fragIndex = $fragCount-1;
                    $fragInfos{$fragIndex} = {
                        "MZ"         => $infos[$colName2Index{'ProductMz'}],
                        "ionType"    => $infos[$colName2Index{'FragmentType'}],
                        "ionResidue" => $infos[$colName2Index{'FragmentSeriesNumber'}],
                        "charge"     => $infos[$colName2Index{'ProductCharge'}],
                        "peptide"    => $infos[$colName2Index{'ModifiedPeptideSequence'}]
                    };
                    
                    my $fragAnnot  = $infos[$colName2Index{'Annotation'}];
                    if($fragAnnot =~ /[bymap]\d{0,3}([-+])?([0-9A-Za-z]+)?/) {
                        if($2) {
                            if(!$modifMass{$2}) {
                                $modifMass{$2} = getModifMass($2, ($unimod) ? "$workDir/$unimod" : $unimodRef);
                            }
                            $fragInfos{$fragIndex}{"modif"} = sprintf("%.0f", $1.$modifMass{$2}) if($1 && $modifMass{2});
                        }
                    }
                }
            }
            
            close LIBEXPORT;
        } else {
            printToFile($fileError, "No library export file found.", 1);
            exit;
        }
    }
    
    ################################
    ###  Parsing alignment file  ###
    ################################
    my (%peptideProt, %peptidePosition, %fragmentsInfos, %decoyNb, %targetNb, %sheetNames, %scoreList, %analysisFileNames, %sampleFileNames, %peptideList, %dbRankProt, %assosProtPeptide);
    if($format eq 'peakview') {
        ### Split the upload file by sheet
        printToFile($fileStat, "Alignment file : Split by analysis", 1);
        system "$promsPath{python}/python $promsPath{python_scripts}/convert_xlsx2txt.py -f $workDir/$alignmentFile >>$workDir/sortie.txt 2>&1";
        
        my ($lastFile,$wait);
        while (!$wait) {
            sleep 20;
            if ( -s "$workDir/$alignmentFileName\_SheetNames.info" && !%sheetNames) {
                open(SHEETNAMES,"<$workDir/$alignmentFileName\_SheetNames.info");
                my $i=1;
                foreach (reverse(<SHEETNAMES>)) {
                    my ($sheet, $file) = split(/\t/,$_);
                    chomp $file;
                    $lastFile = $file if $i==1;
                    $sheetNames{'Area-ions'} = $file if $sheet=~/Area - ions/;
                    $sheetNames{'Area-peptides'} = $file if $sheet=~/Area - peptides/;
                    $sheetNames{'Area-proteins'} = $file if $sheet=~/Area - proteins/;
                    $sheetNames{'Score'} = $file if $sheet=~/Score/;
                    $sheetNames{'FDR'} = $file if $sheet=~/FDR/;
                    $sheetNames{'Observed_RT'} = $file if $sheet=~/Observed RT/;
                    $i++;
                }
                close SHEETNAMES;
            }
            $wait = (`tail -1 $workDir/sortie.txt`=~/Done/ && -s $lastFile) ? 1 : 0;
        }
        sleep 3;

        #######################################################################
        #### Selection of peptide with a fdr < $fdr (selected by the user) ####
        #######################################################################
        ### In the FDR file
        printToFile($fileStat, "Scanning peakview file", 1);
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
                    my $pepRT=sprintf("%.2f", $validPep{'TARGET'}{$peptide}{$i});
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
    } elsif ($format eq 'openswath') {
        printToFile($fileStat, "Alignment file : Openswath format scanning", 1);
        my %colName2Index;
        my $lastPeptide='';
        my $lastAna;
        my %prot;
        
        open (IN, "<$workDir/$alignmentFile");
        while (<IN>){
            if ($.==1){
                $_=~s/\s*\Z//;
                my @columns=split(/[,;\t]/,$_);
                my @refColNames=('ProteinName','filename','mz','Charge','aggr_Fragment_Annotation','RT','aggr_Peak_Area','FullPeptideName');
                foreach my $i (0 .. $#columns){
                    $colName2Index{$columns[$i]}=$i;
                }

                ##check columns names
                foreach my $colName (@refColNames){
                    unless ($colName2Index{$colName}){
                        printToFile($fileError, "Wrong columns names in file: $colName missing ! <BR/>The following columns are required : <BR/>ProteinName<BR/>filename<BR/>mz<BR/>Charge<BR/>aggr_Fragment_Annotation<BR/>RT<BR/>aggr_Peak_Area<BR/>FullPeptideName", 1);
                        exit;
                    }
                }
            } else {
                $_ =~ s/\s*\Z//;
                my @infos = split(/[\t]/,$_);
                
                next if ($infos[$colName2Index{'decoy'}] == 1);
                next if ($infos[$colName2Index{'ProteinName'}] =~ /iRT/);
                next if ($infos[$colName2Index{'ProteinName'}] =~ /^1\/reverse_sp/);
                
                my $proteinsList  = $infos[$colName2Index{'ProteinName'}];
                my $anaFile       = $infos[$colName2Index{'filename'}];
                my $pepMZ         = $infos[$colName2Index{'mz'}];
                my $pepCharge     = $infos[$colName2Index{'Charge'}];
                my $areasList     = $infos[$colName2Index{'aggr_Peak_Area'}];
                my $RT            = $infos[$colName2Index{'RT'}];
                my $modifSeq      = $infos[$colName2Index{'FullPeptideName'}];
                my $mScore        = $infos[$colName2Index{'m_score'}];
                my $fragmentsList = $infos[$colName2Index{'aggr_Fragment_Annotation'}];
                my $peptide       = $modifSeq.'_'.$pepCharge;
                
                #print("Extracting sequence : $infos[$colName2Index{FullPeptideName}] => $modifSeq\n");
                
                $RT = sprintf("%.2f", ($RT/60));
                $areasList =~ s/"//g;
                $fragmentsList =~ s/"//g;
                undef @infos;

                $anaFile =~ s/\\/\//g; # Change backslash by slash for windows format (not recognized by fileparse)
                (my $anaName = $anaFile) =~ s/.mzXML//;
                $anaName = fileparse($anaName, qr/\.[^.]*/);
                $anaFile = fileparse($anaFile, qr/\.[^.]*/);
                $analysisFileNames{$anaName} = $anaFile;
                $sampleFileNames{$anaName} = $anaName;
                
                # For DIA phospho, escape special characters (257097_sp|C0HKE3|H2A1D_MOUSE,257098_sp|C0HKE4|H2A1D_MOUSE; if multiple proteins)
                #my $quoteProteinsList=quotemeta($proteinsList);
                #$fragmentsList=~s/$quoteProteinsList//g;        

                # Recovering peptide infos
                @{$peptideList{$anaName}{$peptide}{"$RT"}} = ($pepMZ, '', $fragmentsList, $areasList);

                # Recovering peptide/protein associations
                if($peptide eq $lastPeptide && $peptideProt{$lastAna}{$peptide}) {
                    $peptideProt{$anaName}{$peptide} = \%{$peptideProt{$lastAna}{$peptide}};
                } else {
                    ## openswath normal => proteinsList like : 2/sp|Q6PIX5|RHDF1_MOUSE/sp|Q6PIX5|RHDF2_MOUSE ||| openswath phospho => proteinsList like :sp|Q6PIX5|RHDF1_MOUSE;sp|Q6PIX5|RHDF2_MOUSE
                    my @proteins = ($proteinsList=~/^\d/) ? split(/[\/]/, $proteinsList) : split(/[;]/,$proteinsList);
                    my $start = ($proteinsList=~/^\d/) ? 1 : 0;
                    for my $i ($start .. $#proteins) {
                        next if($proteins[$i] =~ /reverse/);
                        @{$assosProtPeptide{$anaName}{$proteins[$i]}{$peptide}} = ($RT);
                        $peptideProt{$anaName}{$peptide}{$proteins[$i]} = 1;
                        $prot{$proteins[$i]}=1;
                    }
                    undef @proteins;
                }

                # Recovering peptide modifications
                if ($modifSeq =~ /\(UniMod:\d+\)/) {
                    my $peptideSeq = $modifSeq;
                    while ($peptideSeq =~ /\(UniMod:(\d+)\)/g) {
                        my $modificationID = $unimod2ID{$1};
                        my $pos = (substr($peptideSeq, 0, 1) eq ".") ? 0 : $-[0]; # DIA phospho data can start and end with a dot: .(UniMod:X)CCCGTATAG.(UniMod:X) ('=' is N-ter modif | '*' is C-ter modif)
                        
                        $peptideModif{$anaName}{$peptide}{$modificationID}{$pos} = $massMods{$modificationID};
                        $peptideSeq =~ s/\(UniMod:$1\)//;
                        $peptideSeq =~ s/(^\.|\.$)//; #For DIA Phospho, remove dot of N-term and C-term modifs
                    }
                    undef $peptideSeq;
                }
                
                $lastPeptide = $peptide;
                $lastAna = $anaName;
                
                undef $peptide;
                undef $anaName;
            }
        }
        close IN;
    }
    elsif ($format eq 'spectronaut'){
        printToFile($fileStat, "Alignment file : Spectronaut format scanning", 1);

        my %colName2Index;
        open (IN, "<$workDir/$alignmentFile");
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
                        printToFile($fileError, "Wrong columns names in file ! <BR/>The columns names have to be like : <BR/>R.Label<BR/>PG.UniProtIds<BR/>EG.ModifiedPeptide<BR/>EG.RTPredicted<BR/>FG.Charge<BR/>FG.PrecMz<BR/>F.Charge<BR/>F.FrgIon<BR/>F.FrgMz<BR/>F.PeakArea", 1);
                        exit;
                    }
                }
            } else {
                $_=~s/\s*\Z//;
                my @infos=split(/[\t]/,$_);
                next if $infos[$colName2Index{'PG.UniProtIds'}]=~/iRT/;
                next if $infos[$colName2Index{'PG.UniProtIds'}]=~/^1\/reverse_sp/;
                #if ($infos[$colName2Index{'EG.ModifiedPeptide'}]!~/^_/){
                #    print $infos[$colName2Index{'EG.ModifiedPeptide'}],"<BR/>";
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
                            printToFile($fileError, "The associated modification for the peptide $modifSeq is not found.", 1);
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

    ##########################################################
    ###  Fetching protein sequence and mass from dataBank  ###
    ##########################################################
    $dbh = &promsConfig::dbConnect('no_user');# if $format eq 'peakview'; # reconnect after long process
    my (@dbID, %numMatchProt, %protCoverage, %sequenceProtList, %protLength, %protDes, %protMW, %protOrg);

    printToFile($fileStat, "Databank : Extracting protein sequences", 1);

    ### Fetch DB file
    if($libID) {
        my $sthDBID = $dbh->prepare("SELECT ID_DATABANK,DB_RANK FROM DATABANK_SWATHLIB WHERE ID_SWATH_LIB=?");
        $sthDBID->execute($libID);
        while (my ($dbID, $dbrank) = $sthDBID->fetchrow_array) {
            push @dbID, "$dbID&$dbrank";
        }
        $sthDBID->finish;
    }

    ### Parse all given databases
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

        my ($des,$org);
        my $count=0;
        if (-e "$dbDir/$dbFasta"){
            open(FASTAFILE,"<","$dbDir/$dbFasta") or die ("open: $!");
            my ($checkReverse,$identifier,$matched,$sequence,$info);
            my $nbmatch=0;
            while (my $line=<FASTAFILE>){
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
                        
                        #print("sequenceProtList : $identifier\n"); # TODO Remove

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
                        #print $identifier,"<BR/>";
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
        } else {
            printToFile($fileError, "Fasta file not found! Skipping.", 1);
            exit;
        }
        
        my $nbProt = scalar %protMW; # TODO Remove
        printToFile($fileStat, "Databank : Found $nbProt in the DB fasta file $dbDir/$dbFasta", 1);
    }


    #####################################################################
    ####>Recovering peptide positions on protein and coverage calcul<####
    #####################################################################
    printToFile($fileStat, "Databank : Recovering peptide's positions on protein", 1);

    foreach my $analysisNumber (keys %assosProtPeptide){        # for each analysis
        foreach my $protein (keys %{$assosProtPeptide{$analysisNumber}}){
            my %positionPeptide;
            my $sequence=($sequenceProtList{$protein}) ? $sequenceProtList{$protein} : '';
            my $numMatch=0;
            my (@begPep,@endPep);
            
            foreach my $peptide (keys %{$assosProtPeptide{$analysisNumber}{$protein}}){
                my @pepSeq=split(/_/,$peptide);
                my $pepSequence=$pepSeq[0];
                $pepSequence=~s/\[\+?\w+\]//g if ($format eq 'peakview' || $format eq 'spectronaut');
                $pepSequence=~s/\(UniMod:\d+\)//g if ($format eq 'openswath');
                $pepSequence=~s/\.//g if ($format eq 'openswath');
                
                if(!$pepSequence) {
                    printToFile($fileStat, "Databank : Peptide not found for $protein => $peptide", 1); # TODO Remove
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
    
    
    ###############################
    ### Import data in database ###
    ###############################
    my $quantiName=($format eq 'peakview')? 'PeakView_Quantification' : ($format eq 'openswath')? 'OpenSwath_Quantification' : 'Spectronaut_Quantification';
    my $importManager = promsImportDataManager->new($userID, $projectID, $experimentID, $format, $quantiName, $quantifAnnot, $taxonomy, $instrument, $alignmentFile, $fileStat, \@dbID);
    $importManager->setFilesInfos(\%analysisFileNames, \%sampleFileNames);
    $importManager->setDIAInfos($libID, \%decoyNb, \%targetNb, $fdr);
    
    printToFile($fileStat, "Storing : Insert results in Database", 1);
    $importManager->importData(\%peptideInfo, \%peptideProt, \%peptideModif, \%peptidePosition, \%peptideList, \%fragmentsInfos, \%modifIDs, \%assosProtPeptide, \%dbRankProt);

    my %analysisID = %{$importManager->getAnalysisID};
    my %peptidesID = %{$importManager->getPeptidesID};
    my $quantiID = $importManager->getQuantiID;
    

    ############################
    ### Separating data file ###
    ############################
    my $warning = "";
    my $mzXmlFileIndex = 0;
    my $nbAnalysisFiles = scalar keys %analysisFileNames;
    printToFile($fileStat, "Storing : Split data files by analysis", 1);

    foreach my $analysis(keys %analysisFileNames) {
        my $analysisID = $analysisID{$analysis};
        my $anaDir = "$promsPath{peptide}/proj_$projectID/ana_$analysisID";
        my $quantiPath = $promsPath{"quantification"}."/project_$projectID/quanti_$quantiID";
        my $swathFile = "$quantiPath/swath_ana_$analysisID.txt";

        mkdir $promsPath{'peptide'} unless -e $promsPath{'peptide'};
        mkdir "$promsPath{peptide}/proj_$projectID" unless -e "$promsPath{peptide}/proj_$projectID";
        mkdir $anaDir unless -e $anaDir;
        system "mkdir -p $quantiPath";

        # Move the mzXML file to the corresponding analysis folder
        my $mzXMLFile = "$analysis.mzXML";
        if(-s "$workDir/$mzXMLFile") {
            printToFile($fileStat, "Storing : Analysis file for $analysis ($mzXmlFileIndex/$nbAnalysisFiles processed)", 1);
            if(`stat -c "%s" $workDir/$mzXMLFile` > 10000000) { ### compress mzXML file if its size > 1Go
                system "cd $workDir; tar -zcf $mzXMLFile.tar.gz $mzXMLFile";
                move("$workDir/$mzXMLFile.tar.gz", "$anaDir/$mzXMLFile.tar.gz");
            } else {
                move("$workDir/$mzXMLFile", "$anaDir/$mzXMLFile");
            }
        }
        
        # Move the PyProphet file to the corresponding analysis folder
        my $pyprophetOutFile = "$mzXMLFile.tsv";
        if(-s "$workDir/$pyprophetOutFile") {
            printToFile($fileStat, "Storing : PyProphet file for $analysis ($mzXmlFileIndex/$nbAnalysisFiles processed)", 1);
            if(`stat -c "%s" $workDir/$pyprophetOutFile` > 10000000) { ### compress mzXML file if its size > 1Go
                system "cd $workDir; tar -zcf $pyprophetOutFile.tar.gz $pyprophetOutFile";
                move("$workDir/$pyprophetOutFile.tar.gz", "$anaDir/$pyprophetOutFile.tar.gz");
            } else {
                move("$workDir/$pyprophetOutFile", "$anaDir/$pyprophetOutFile");
            }
        }
        
        printToFile($fileStat, "Storing : XIC results file for $analysis ($mzXmlFileIndex/$nbAnalysisFiles processed)", 1);
        open(OUTFILE, ">", $swathFile) or die ("open: $!");
        
        if($format eq 'peakview') {
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
                        my $peptide;
                        foreach my $RT (keys %{$peptidesID{$analysis}{$partLine[1].'_'.$partLine[3]}}) {
                            my @peptides = $peptidesID{$analysis}{$partLine[1].'_'.$partLine[3]}{$RT};
                            if(@peptides) {
                                $peptide = $peptidesID{$analysis}{$partLine[1].'_'.$partLine[3]}{$RT}[0];
                                last;
                            }
                        }
                        push @newLine, ($peptide) ? $peptide : '';
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

        } elsif($format eq 'openswath') {
            my $fragNoMZ = 0;
            my $fragTot = 0;
            open(MISSMZ,">>$workDir/remaining_mz.csv"); # TODO Remove
            
            print OUTFILE "ID_PEPTIDE!Fragment MZ!Fragment Charge!Ion Type!Residue!Area!RT\n";
            
            foreach my $peptide (keys %{$peptideList{$analysis}}) {
                #my $ghostPep = 0;
                foreach my $rt (sort {$peptideList{$analysis}{$peptide}{$a} <=> $peptideList{$analysis}{$peptide}{$b}} keys %{$peptideList{$analysis}{$peptide}}) {
                    my ($pepMZ, $score, $fragmentsList, $areasList) = @{$peptideList{$analysis}{$peptide}{$rt}};
                    #my $peptideID = ($ghostPep) ? $peptidesID{$analysis}{$peptide}{$RT}[1] : $peptidesID{$analysis}{$peptide}{$RT}[0];
                    my $peptideID = $peptidesID{$analysis}{$peptide}{$rt}[0];
                    my @areas = split(/;/, $areasList);
                    my @fragments = split(/;/, $fragmentsList);
                    my ($fragmentIonResidue, $fragmentIonType, $fragmentCharge, $fragmentMZ, $fragmentID, $fragmentArea, $fragmentPeptide);
                    
                    foreach my $i (0 .. $#fragments) {
                        my @lineFragInfo = split(/_/, $fragments[$i]);
                        if(!@lineFragInfo) {
                            print("Can't find value for $peptide ($fragments[$i])\n<br/>");
                        }
                        
                        $fragmentID = ($lineFragInfo[0]) ? $lineFragInfo[0] : "";
                        next if(!$fragmentID);
                        
                        if(!$fragInfos{$fragmentID}) {
                            $fragNoMZ++;
                            print MISSMZ "$peptide & $fragments[$i] : $fragmentID not found\n"; # TODO Remove
                            $fragTot++;
                            next;
                        }
                        
                        my $lineFragIonResidue = $lineFragInfo[1];
                        my $lineFragIonCharge = $lineFragInfo[2];
                        
                        $fragmentIonType = $fragInfos{$fragmentID}{"ionType"};
                        $fragmentIonResidue = $fragInfos{$fragmentID}{"ionResidue"};
                        $fragmentIonResidue .= $fragInfos{$fragmentID}{"modif"} if($fragInfos{$fragmentID}{"modif"});
                        $fragmentCharge = $fragInfos{$fragmentID}{"charge"};
                        $fragmentPeptide = $fragInfos{$fragmentID}{"peptide"};
                        $fragmentMZ = $fragInfos{$fragmentID}{"MZ"};
                        
                        next if($fragmentIonResidue eq '-1');
                        
                        $fragmentArea = $areas[$i];
                        $fragmentArea =~ s/,/\./ if ($fragmentArea =~ /,/);
                        print OUTFILE "$peptideID!$fragmentMZ!$fragmentCharge!$fragmentIonType!$fragmentIonResidue!$fragmentArea!$rt\n";
                        
                        
                        if($fragTot % 10000 == 0) {
                            printToFile($fileStat, "Storing : XIC results file for $analysis ($mzXmlFileIndex/$nbAnalysisFiles file processed - $fragTot fragments stored)", 1);
                        }
                        $fragTot++;
                    }
                    #$ghostPep = 1;
                }
            }
            
            if($fragNoMZ) {
                $warning = "$analysis : $fragNoMZ/$fragTot fragments have no M/Z information.<br/>"; 
            }
            
            close OUTFILE;
            close MISSMZ; # TODO Remove
        } else {
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
                        my $fragmentMZ=($fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[2]) ? $fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[2] : '';
                        my $fragmentNeutralMass=$fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[3];
                        my $area=$fragmentsInfos{$analysis}{$peptide}{$rt}{$fragment}[1];
                        my $peptideID=($ghostpep) ? $peptidesID{$analysis}{$peptide}{$rt}[1] :  $peptidesID{$analysis}{$peptide}{$rt}[0] ;
                        print OUTFILE "$peptideID!$fragmentMZ!$fragmentCharge!$ionType!$residue!$area!$fragmentRT\n";
                    }
                    $ghostpep=1;
                }
            }
            close OUTFILE;
        }
        
        symlink($swathFile, "$anaDir/$alignmentFile");
        $mzXmlFileIndex++;
    }
    
    printToFile($fileStat, "OpenSwath Workflow Ended.", 1);
    if($warning) {
        printToFile($fileStat, "WARNING : $warning", 1) if($warning);
    }
    $dbh->disconnect;
}



#########################
###     METHODS       ###
#########################
sub excludeLibraryModifications {
    my $libID = param('libID');
    if($libID) {
        my $sthModification = $dbh->prepare("SELECT SLM.ID_MODIFICATION, PSI_MS_NAME FROM MODIFICATION M, SWATH_LIB_MODIFICATION SLM WHERE ID_SWATH_LIB=$libID AND SLM.ID_MODIFICATION=M.ID_MODIFICATION");
        $sthModification->execute;
        if ($sthModification->rows != 0) {
            print "<BR/><B>&nbsp;- Peptide modifications to be excluded :</B><BR/>";
            while(my ($modifID,$modifName)=$sthModification->fetchrow_array) {
                print "<INPUT type=\"checkbox\" name=\"pepMod\" value=\"$modifID\">$modifName&nbsp;";
            }
        }
        $dbh->disconnect;
    }
}

sub getModifMass {
    my ($modifFormula, $unimodFile) = @_;
    my $xml  = new XML::Simple;
    my $data = $xml->XMLin($unimodFile);

    for my $modif (@{$data->{"umod:modifications"}{"umod:mod"}}) {
        next unless($modif->{"umod:specificity"});
        my @specificities;
        if(ref($modif->{"umod:specificity"}) eq 'ARRAY') {
            @specificities = @{$modif->{"umod:specificity"}};
        } else {
            push @specificities, $modif->{"umod:specificity"};
        }

        for my $specificity (@specificities) {
            next unless($specificity->{"umod:NeutralLoss"});
            my @neutralLoss;
            if(ref($specificity->{"umod:NeutralLoss"}) eq 'ARRAY') {
                @neutralLoss = @{$specificity->{"umod:NeutralLoss"}};
            } else {
                push @neutralLoss, $specificity->{"umod:NeutralLoss"};
            }

            for my $neutralLoss (@neutralLoss) {
                next unless($neutralLoss->{"umod:element"});
                my @neutralLossElements;
                my $neutralLossFormula = "";
                my $remainingFormula = $modifFormula;
                if(ref($neutralLoss->{"umod:element"}) eq 'ARRAY') {
                    @neutralLossElements = @{$neutralLoss->{"umod:element"}};
                } else {
                    push @neutralLossElements, $neutralLoss->{"umod:element"};
                }
                
                for my $neutralLossElement (@neutralLossElements) {
                    my $element = $neutralLossElement->{"symbol"}.$neutralLossElement->{"number"};
                    $neutralLossFormula .= $element;
                    $remainingFormula =~ s/$element//;
                }
                
                if($remainingFormula eq "") {
                    return $neutralLoss->{"mono_mass"};
                }
            }
        }
        
    }
    return -1;
}

sub castToUniMod {
    my ($alignmentFile, $unimodFile) = @_;
    my $xml  = new XML::Simple;
    my $data = $xml->XMLin($unimodFile);

    for my $modif (@{ $data->{"umod:modifications"}{"umod:mod"} }) {
        my $type      = $modif->{"username_of_poster"};
        my $modifName = $modif->{"title"};
        my $uniModID  = $modif->{"record_id"};

        next unless $type eq 'unimod';
        
        # Change current PTM in alignment results to UniMod text format : (Phospho) -> (UniMod:21)
        system "sed -i 's/\($modifName\)/\(UniMod:$uniModID\)/g' $alignmentFile";
    }
}

sub selectSamples {
    if($experimentID) {
        my $sthListSamples = $dbh->prepare("SELECT ID_SAMPLE, NAME FROM SAMPLE WHERE ID_EXPERIMENT=?");
        $sthListSamples->execute($experimentID);
        if ($sthListSamples->rows > 0) {
            my $sthListAna = $dbh->prepare("SELECT ID_ANALYSIS, WIFF_FILE FROM ANALYSIS A WHERE ID_SAMPLE=?");
            
            while (my ($sampleID, $sampleName) = $sthListSamples->fetchrow_array) {
                $sthListAna->execute($sampleID);
                
                if ($sthListAna->rows > 0) {
                    while(my ($anaID, $wiffFile) = $sthListAna->fetchrow_array) {
                        my $anaDir = "$promsPath{peptide}/proj_$projectID/ana_$anaID";
                        my $mzXMLFile = "$wiffFile.mzXML";
                    
                        if((!$sameParams && (-s "$anaDir/$mzXMLFile" || -s "$anaDir/$mzXMLFile.tar.gz"))
                           || ($sameParams && (-s "$anaDir/$mzXMLFile.tsv" || -s "$anaDir/$mzXMLFile.tsv.tar.gz"))) {
                            print "<INPUT type=\"checkbox\" class='samples' name=\"samples\" onclick=\"checkSelectedSamples($experimentID)\" value=\"$sampleID\" checked>$sampleName&nbsp;";
                            last;
                        }
                    }
                }
            }
        }
    }
}

sub uploadFiles {
    my ($workDir, $typeFile) = @_;
    my (@srcFilesPaths, @files);
    
    ### Uploading files of defined type
    my $paramValue = $typeFile."Files";
    if(param($paramValue)) {
        if (ref(param($paramValue)) eq "Fh" || (!$isOnCluster && $shouldRunOnCluster)) {
            foreach my $file (param($paramValue)) {
                push(@srcFilesPaths, $file);
            }
        } else {
            @srcFilesPaths = split(/OTHERFILE/, param($paramValue));
        }
        
        foreach my $file (@srcFilesPaths) {
            my $newFileDir = "$workDir/".basename($file);
            
            if($paramValue ne 'sharedDirFiles') { # FROM LOCAL PATH
                uploadFile($workDir, $file);
            } elsif(-e "$promsPath{shared}/$file" && not -e $newFileDir) { # FROM SHARED PATH
                copy("$promsPath{shared}/$file", $newFileDir); # TODO Change with move
            }
            
            push(@files, basename($file));
            
            my $size=`stat -c "%s" $newFileDir`;
            if ($size == 0) {
                print "<HTML><BODY><SCRIPT LANGUAGE=\"JavaScript\">alert('Error during importation of upload files ($newFileDir) Please retry.'); window.location=\"$promsPath{cgi}/importSwathDataRefonte.cgi?id_project=$projectID&ID=$branchID&ACT=quantification&FORMAT=openswath\";</SCRIPT></BODY></HTML>";
                exit;
            }
            print("Imported $typeFile file $file.<br/>");
        }
    }

    return @files;
}

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
                    }
                    
                    th {
                        min-width: 300px;
                    }
                </STYLE>
    |;
    
    printJS();
    
    print qq |
            </HEAD>
            <BODY background="$promsPath{images}/bgProMS.gif">
                <CENTER>
    |;
    
    if ($action eq 'import') {
        printImportForm();
    }
    elsif ($action eq 'quantification') {
        printQuantiForm();
    }
    
    print qq |
                </CENTER>
            </BODY>
        </HTML>
    |;
}

sub printJS {
    print qq |
        <SCRIPT LANGUAGE="JavaScript">
    |;
    
    &promsMod::popupInfo();
    &promsMod::browseDirectory_JavaScript if $promsPath{'shared'};
    
    print qq |
            function cancelAction() {
                window.location="./processAnalyses.cgi?ID=$branchID";
            }
            
            function selectExportOption(value) {
                if (value == 'export') {
                    document.getElementById('exportparametersSpectrast2tsv').style.display = '';
                    document.getElementById('exportparamfile').style.display = 'none';
                    document.getElementById('importlibPQP').style.display = 'none';
                    document.getElementById('importlibTSV').style.display = 'none';
                    
                    if(!document.getElementById("enableIPF") \|\| !document.getElementById("enableIPF").checked) {
                        document.getElementById('exportLibNonIPFOptions').style.display = '';
                    }
                } else if (value == 'import') {
                    document.getElementById('exportparametersSpectrast2tsv').style.display = 'none';
                    document.getElementById('exportLibNonIPFOptions').style.display = 'none';
                    document.getElementById('exportparamfile').style.display = '';
                    document.getElementById('importlibPQP').style.display = '';
                    document.getElementById('importlibTSV').style.display = '';
                }
            }
            
            function checkform(form) {
                var enabledIPF = (document.getElementById('enableIPF') && document.getElementById('enableIPF').checked) ? 1 : 0;
                if(document.getElementById('exportparamfile').style.display=='' && form.exportparam.files.length==0){
                    alert('ERROR: Select the export parameters file.');
                    return false;
                }
                
                if(document.getElementById('importlibTSV').style.display=='' && form.exportedlibrary.files.length==0){
                    alert('ERROR: Select the library TSV file.');
                    return false;
                }
                
                if(document.getElementById('importlibPQP').style.display=='' && form.decoylibrary.files.length==0){
                    alert('ERROR: Select the library PQP file.');
                    return false;
                }
                
                if(!enabledIPF && document.getElementById('exportparametersSpectrast2tsv').style.display=='' && form.swathfile.files.length==0){
                    alert('ERROR: Select the swath windows file.');
                    return false;
                }
                
                if(document.getElementById("useMs1Traces")) {
                    document.getElementById("useMs1Traces").disabled = false;
                }
                
                if(document.getElementById("enableUisScoring")) {
                    document.getElementById("enableUisScoring").disabled = false;
                }
            }
            
            function selectSource(source) {
                var myForm = document.parameters;
                myForm.mzXMLFiles.disabled = true;
                if (document.getElementById('sharedDirDIV')) {
                    document.getElementById('sharedDirDIV').style.display='none';
                }
            
                if(source=='UseLocalDirectory') {
                    myForm.mzXMLFiles.disabled=false;
                } else if(source=='UseSharedDirectory') {
                    document.getElementById('sharedDirDIV').style.display='';
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
            
            function selectLibrary (libID) {
                var XHR = getXMLHTTP();
                XHR.onreadystatechange=function () {
                    if(XHR.readyState==4 && XHR.responseText) {
                        if(document.getElementById("enableIPF") && document.getElementById("enableIPF").checked) {
                            document.getElementById('libraryModifications').innerHTML="";
                        } else {
                            document.getElementById('libraryModifications').innerHTML=XHR.responseText;
                        }
                    } else {
                        document.getElementById('libraryModifications').innerHTML="";
                    }
                }
                XHR.open("GET","$promsPath{cgi}/importSwathDataRefonte.cgi?ACT=excludeLibraryModifications&ID=$branchID&libID="+libID,true);
                XHR.send(null);
            }
            
            function ajaxSelectSamples(expID) {
                var listSample=document.getElementById('listSamples'+expID);
                var listAllSamples = document.getElementById('allSamples'+expID);
                var sameParams = (document.getElementById('sameParams').checked) ? 'true' : '';
                
                if(listSample.innerHTML == '') {
                    listSample.style.display = '';
                    var XHR = getXMLHTTP();
                    XHR.onreadystatechange=function () {
                        if(XHR.readyState==4 && XHR.responseText ){
                            listSample.innerHTML = XHR.responseText;
                        }
                    }
                    XHR.open("GET","$promsPath{cgi}/importSwathDataRefonte.cgi?ACT=selectSamples&ID=$branchID&id_project=$projectID&sameParams=" + sameParams + "&expID=" + expID, true);
                    XHR.send(null);
                } else {
                    listSamplesCheckbox = listSample.getElementsByTagName("input");
                    for(var i=0; i<listSamplesCheckbox.length; i++) {
                        listSamplesCheckbox[i].checked = listAllSamples.checked;
                    }
                }
            }
            
            function resetExperiments() {
                var allListSamples = document.getElementsByClassName('listSamples');
                for (var i = 0; i < allListSamples.length; ++i) {
                    var listSample = allListSamples[i];  
                    listSample.innerHTML = '';
                }
                
                experiments = document.getElementsByClassName('experiment');
                for (var i=0; i < experiments.length; i++)  {
                    experiments[i].checked = false;
                }
            }
            
            function checkSelectedSamples(expID) {
                var listSample = document.getElementById('listSamples'+expID);
                var listAllSamples = document.getElementById('allSamples'+expID);
                var nbChecked = listSample.querySelectorAll('input[type="checkbox"]:checked').length;
                var nbTot = listSample.querySelectorAll('input[type="checkbox"]').length;
                listAllSamples.checked = (nbChecked == nbTot);
            }
            
            function checkIPF() {
                var checkedIPF = document.getElementById("enableIPF").checked;
                var checkedExportOption = document.getElementById("libraryExportOption").checked;
                 
                if(checkedIPF) {
                    if(checkedExportOption) {
                        document.getElementById('exportLibNonIPFOptions').style.display = 'none';
                    }
                    
                    document.getElementById("exportLibIPFOptions").innerHTML = '<B>- IPF options:</B><br/> Unimod modifications file : <INPUT onmouseover="popup(\\'File with the modifications not specified by default.\\')" onmouseout="popout()" type="file" name="unimodFile" ><br/> Maximum number of alternative localizations : <INPUT onmouseover="popup(\\'Maximum number of site-localization permutations.\\')" onmouseout="popout()" type="text" name="maxNumLocalization" value="10000" size="5">';
                    document.getElementById('libraryModifications').innerHTML = '';
                } else {
                    if(checkedExportOption) {
                        document.getElementById('exportLibNonIPFOptions').style.display = '';
                    }
                    
                    document.getElementById("exportLibIPFOptions").innerHTML = "";
                    selectLibrary(document.getElementById('libID').value);
                }
                
                document.getElementById("useMs1Traces").checked = checkedIPF;
                document.getElementById("useMs1Traces").disabled = checkedIPF;
                document.getElementById("enableSpLosses").checked = checkedIPF;
                document.getElementById("enableSpLosses").disabled = checkedIPF;
                document.getElementById("enableUisScoring").checked = checkedIPF;
                document.getElementById("enableUisScoring").disabled = checkedIPF;
            }
        </SCRIPT>
    |;
}

sub printImportForm {
    my $formTitle = ($format eq 'peakview') ? "Peakview" : ($format eq 'openswath') ? "OpenSwath" : "Spectronaut" ;
    my $fileFormat = ($format eq 'peakview') ? ".xlsx" : ($format eq 'openswath') ? ".tsv" : ".tsv" ;
    my $clusterValue = ($format eq 'openswath') ? 1 : 0;
    
    print qq |
        <FONT class="title1">Select options to import $formTitle data</FONT><BR/><BR/><BR/><BR/>
        <FORM method="POST" action="./importSwathDataRefonte.cgi" name="parameters" enctype="multipart/form-data">
            <INPUT name="ID" type="hidden" value="$branchID">
            <INPUT name="id_project" type="hidden" value="$projectID">
            <INPUT name="ACT" type="hidden" value="$action">
            <INPUT name="FORMAT" type="hidden" value="$format">
            <INPUT name="RUNCLUSTER" type="hidden" value="$clusterValue">
            <INPUT name="USERID" type="hidden" value="$userID">
            <TABLE bgcolor="$darkColor">
                <TR>
                    <TH align=right valign=top>Result file : </TH>
                    <TD bgcolor="$lightColor"> <INPUT type="file" name="uploadfile" accept="$fileFormat" required> </TD>
                </TR>
                <TR>
                    <TH align=right valign=top>Library name : </TH>
                    <TD bgcolor="$lightColor">
    |;
    
    my %libList;
    my $sthLibList=$dbh->prepare("SELECT NAME,ID_SWATH_LIB,IDENTIFIER_TYPE FROM SWATH_LIB WHERE USE_STATUS='yes'") or die "Couldn't prepare statement: " . $dbh->errstr;
    $sthLibList->execute;
    if ($sthLibList->rows==0) {
        print "No library available.";
    } else {
        while (my($libName,$libID,$identType)=$sthLibList->fetchrow_array) {
            $identType=($identType)? $identType : 'Not defined';
            $libList{$libName}="$libID&$identType";
        }
        $sthLibList->finish;
        $dbh->commit;

        print qq |
                        <SELECT name="libID" id="libID" required>
                            <option value="">-= Select Library =-</option>
        |;
        foreach my $libName (sort {lc($a) cmp lc($b)} keys %libList) {
            my ($libID,$identType)=split('&',$libList{$libName});
            print "         <option value=\"$libID\">$libName [$identType]</option>";
        }
        print "         </SELECT> ";
    }
    print qq |
                    </TD>
                </TR>
                <TR>
                    <TH align=right valign=top>Library export parameter file : </TH>
                    <TD bgcolor="$lightColor"><INPUT type="file" name="exportparam" required></TD>
                </TR>
    |;
    
    if($format eq 'peakview') {
        print qq |
                <TR>
                    <TH  align=right valign=top>FDR : </TH>
                    <TD bgcolor="$lightColor"><INPUT type="number" name="fdr" min="0" max="20" size="2" required> %</TD>
                </TR>
                <TR>
                    <TH align=right valign=top>PeakView parameters : </TH>
                    <TD bgcolor="$lightColor"><B>Peptide filter :</B><BR/>
                        Number of peptides per protein : <INPUT type="text" name="nbpeptide" size="3"><BR/>
                        Number of transitions per peptide : <INPUT type="text" name="nbtransition" size="3"><BR/>
                        <INPUT type="checkbox" name="modpep">Exclude modified peptides<BR/><BR/>
                        <B>XIC options : </B><BR/>
                        XIC extraction window (min) : <INPUT type="text" name="XICmin" size="3"><BR/>
                        XIC width : <INPUT type="text" name="XIC" size="3">
                        <SELECT name="XICunit">
                            <OPTION value="DA">Da</OPTION>
                            <OPTION value="PPM">ppm</OPTION>
                        </SELECT>
                    </TD>
                </TR>
        |;
    } else {
        print qq |
                <TR>
                    <TH align=right valign=top>Library export file : </TH>
                    <TD bgcolor=\"$lightColor\"><INPUT type=\"file\" name=\"exportedlibrary\" accept=\".tsv,.mrm\" required ></TD>
                </TR>
        |;
    }
    
    if ($format eq 'openswath') {
        print qq |
                <TR>
                    <TH align=right valign=top>TRIC method used: </TH>
                    <TD bgcolor="$lightColor">
                        <SELECT name="TRICMethod" required>
                            <OPTION value="LocalMST" selected>Local MST*</OPTION>
                            <OPTION value="best_overall">Best Overall</OPTION>
                            <OPTION value="best_cluster_score">Best Cluster Score</OPTION>
                            <OPTION value="global_best_overall">Global Best Overall</OPTION>
                            <OPTION value="global_best_cluster_score">Global Best Cluster Score</OPTION>
                        </SELECT>
                        <SMALL>*Recommanded option</SMALL>
                    </TD>
                </TR>
                <TR>
                    <TH align=right valign=top>Unimod File (IPF): </TH>
                    <TD bgcolor="$lightColor">
                        <INPUT onmouseover="popup(\\'File with the modifications not specified by default.\\')" onmouseout="popout()" type="file" name="unimodFile" >
                    </TD>
                </TR>
        |;
    }
    
    print qq |
                <TR>
                    <TH align=right valign=top>Software version : </TH>
                    <TD bgcolor=\"$lightColor\"><INPUT type=\"text\" name=\"softversion\" required><SMALL>&nbsp;&nbsp;ex : 1.2, 2.1.3 ...</SMALL></TD>
                </TR>
                <TR>
                    <TH colspan=2>
                        <input type="submit" name="submit" value="Submit">
                    
                        <!-- CLEAR button -->
                        &nbsp;&nbsp;&nbsp;<INPUT type="reset" value="Clear" />
                    
                        <!-- CANCEL button -->
                        &nbsp;&nbsp;&nbsp;<INPUT value="Cancel" onclick="cancelAction();" type="button" >
                    </TH>
                </TR>
            </TABLE>
        </FORM>
    |;
}

sub printQuantiForm {
    my $clusterValue = ($format eq 'openswath') ? 1 : 0;

    print qq |
        <BR/><FONT class="title1">OpenSwath quantification</FONT><BR/><BR/><BR/>
        <FORM method="POST" name="parameters" action="./importSwathDataRefonte.cgi" enctype="multipart/form-data" onsubmit="return(checkform(this));">
            <INPUT name="ID" type="hidden" value="$branchID">
            <INPUT name="id_project" type="hidden" value="$projectID">
            <INPUT name="ACT" type="hidden" value="$action">
            <INPUT name="FORMAT" type="hidden" value="$format">
            <INPUT name="RUNCLUSTER" type="hidden" value="$clusterValue">
            <INPUT name="USERID" type="hidden" value="$userID">
            <TABLE bgcolor="$darkColor">
            
                <!-- IPF SELECTION -->
                <TR>
                    <TH align=right valign=top>OpenSwath modification option : </TH>
                    <TD bgcolor="$lightColor">
                        <INPUT type="checkbox" name="enableIPF" id="enableIPF" onclick="checkIPF();">Use IPF (Inference peptido-form)
                    </TD>
                </TR>
                
                <!-- LIBRARY SELECTION -->
                <TR>
                    <TH align=right valign=top>Library name : </TH>
                    <TD bgcolor="$lightColor">&nbsp;
    |;
            
    my $sthLibList=$dbh->prepare("SELECT NAME,ID_SWATH_LIB,IDENTIFIER_TYPE FROM SWATH_LIB WHERE USE_STATUS='yes'") or die "Couldn't prepare statement: " . $dbh->errstr;
    $sthLibList->execute;
    my %libList;
    if ($sthLibList->rows==0) {
        print "No library available.";
    } else {
        while (my($libName,$libID,$identType)=$sthLibList->fetchrow_array) {
            $identType=($identType)? $identType : 'Not defined';
            $libList{$libName}="$libID&$identType";
        }
        $sthLibList->finish;
        $dbh->commit;

        print qq |
                        <SELECT name="libID" id="libID" onchange="selectLibrary(this.value);" required>
                            <option value="">-= Select Library =-</option>
        |;
        foreach my $libName (sort {lc($a) cmp lc($b)} keys %libList) {
            my ($libID,$identType)=split('&',$libList{$libName});
            print "         <option value=\"$libID\">$libName [$identType]</option>";
        }
        print "         </SELECT> ";
    }
    
    print qq |
                    </TD>
                </TR>
                
                <!-- EXPORT SELECTION -->
                <TR>
                    <TH align=right valign=top>Library export management :</TH>
                    <TD bgcolor="$lightColor">
                        <INPUT type="radio" id="libraryExportOption" name="libraryfileoption" value="export" onclick="selectExportOption(this.value)" required>Use the selected library
                        <INPUT type="radio" name="libraryfileoption" value="import" onclick="selectExportOption(this.value)">Import your own library formated for Openswath <BR/>
                    </TD>
                </TR>
                
                <!-- EXPORT OPTIONS (CUSTOM LIB) -->
                <TR id="exportparamfile" style="display:none">
                    <TH align=right valign=top>Library export parameter file : </TH>
                    <TD bgcolor="$lightColor"><INPUT type="file" name="exportparam" ></TD>
                </TR>
                
                <TR id="importlibTSV" style="display:none">
                    <TH align=right valign=top>Library export file (TSV): </TH>
                    <TD bgcolor="$lightColor"><INPUT type="file" name="exportedlibrary" accept=".tsv"></TD>
                </TR>
                
                <TR id="importlibPQP" style="display:none">
                    <TH align=right valign=top>Library decoy file (PQP): </TH>
                    <TD bgcolor="$lightColor"><INPUT type="file" name="decoyLibrary" accept=".pqp"></TD>
                </TR>
                
                
                <!-- EXPORT OPTIONS (SELECTED LIB) -->
                <TR id="exportparametersSpectrast2tsv" style="display:none">
                    <TH align=right valign=top>Library export global options : </TH>
                    <TD bgcolor="$lightColor">
                        <B>&nbsp;- Mass range of fragment ions : </B><BR/>
                            Min : <INPUT type="text" name="-lmin" size="5" value="350">&nbsp;&nbsp;Max : <INPUT type="text" name="-lmax" size="5" value="2000"><BR/>
    
                        <BR/><B>&nbsp;- Ion series and charge : </B><BR/>
                            Ions* : <INPUT type="text" name="-s" size="5" value="b,y">
                            Charge : <INPUT type="text" name="-x" size="5" value="1,2"><BR/><SMALL>&nbsp;*(separated by ',') for example : 'b,y' </SMALL><BR/>
    
                        <BR/><B>&nbsp;- Number of ions per peptide : </B><BR/>
                            Min : <INPUT type="text" name="-o" size="4" value="3"> Max :<INPUT type="text" name="-n" value="6" size="4"><BR/>
                    </TD>
                </TR>
                
                <TR id="exportLibNonIPFOptions" style="display:none">
                    <TH align=right valign=top>Library export specific options : </TH>
                    <TD bgcolor="$lightColor">
                        <B>&nbsp;- Files : </B><BR/>
                            Windows SWATH file : <INPUT onmouseover="popup('File containing the swath ranges.')" onmouseout="popout()" type="file" name="swathfile" ><BR/>
                            File with modifications delta mass* : <INPUT onmouseover="popup('File with the modifications not specified by default.')" onmouseout="popout()" type="file" name="-m" ><br/>
                            <SMALL>*File headers : modified-AA TPP-nomenclature Unimod-Accession ProteinPilot-nomenclature is_a_labeling composition-dictionary<br/>
                            An example : S S[167] 21 [Pho] False {'H':1,'O':3,'P':1}</SMALL><br/>
                            Labelling file : <INPUT onmouseover="popup('File containing the amino acid isotopic labelling mass shifts.')" onmouseout="popout()" type="file" name="-i"><br/>
                            Fasta file : <INPUT onmouseover="popup('Fasta file to relate peptides to their proteins.')" onmouseout="popout()" type="file" name="-f" accept=".fasta"><br/>
                        
                        <BR/><B>&nbsp;- Other options : </B><BR/>    
                            <INPUT type="checkbox" name="other" value="-d">Remove duplicate masses from labelling<br/>
                            <INPUT type="checkbox" name="other" value="-e">Use theoretical mass<br/>
                            Time scale : <INPUT type="radio" name="-t" value="seconds" checked>seconds<INPUT type="radio" name="-t" value="minutes">minutes<br/>
                            UIS order :  <INPUT onmouseover="popup('When using a switching modification, this determines the UIS order to be calculated.')" onmouseout="popout()" type="text" name="-y" value="2" size="3"><br/>
                            Maximum permissible error : <INPUT onmouseover="popup('Maximum error allowed at the annotation of a fragment ion.')" onmouseout="popout()" type="text" name="-p" value="0.05" size="4"><br/>
                            Allowed fragment mass modifications : <INPUT type="text" name="-g" size="10" onmouseover="popup('List of allowed fragment mass modifications. Useful for phosphorylations.')" onmouseout="popout()"><br/>
                            
                        <DIV id="libraryModifications"></DIV>
                    </TD>
                </TR>
                
                <TR id="exportparametersAssayDecoy">
                    <TH align=right valign=top>Assays and Decoys generation options : </TH>    
                    <TD bgcolor="$lightColor">
                        Product mz threshold : <INPUT type="text" name="productMZThresh" value="0.025" size="4"><br/>
                        Precursor mz threshold : <INPUT type="text" name="precursorMZThresh" value="0.025" size="4"><br/>
                        Precursor lower mz limit : <INPUT type="text" name="precursorLowerMZ" value="400" size="4"><br/>
                        Precursor upper mz limit : <INPUT type="text" name="precursorUpperMZ" value="1000" size="4"><br/>
                        
                        <INPUT type="checkbox" name="enableSpLosses" id='enableSpLosses' value="enableSpLosses" onmouseover="popup('Use specific losses (Required for Phosphorylation analyses and IPF)')">Enable specific losses<br/>
                        <INPUT type="checkbox" name="enableUnspLosses" value="enableUnspLosses">Enable unspecific losses<br/>
                        <br/>
                        <DIV id="exportLibIPFOptions"></DIV>
                    </TD>
                </TR>
                
                <TR>
                    <TH align=right valign=top>mzXML files : </TH>
                    <TD bgcolor="$lightColor">
                        <INPUT type="radio" name="selSource" value="UseLocalDirectory" onclick="selectSource(this.value);"><B>Upload multiple files</B><BR/>&nbsp;&nbsp;&nbsp;<INPUT type="file" name="mzXMLFiles" accept=".mzXML" required multiple="multiple" disabled><BR/><BR/>
                        <INPUT type="radio" name="selSource" value="UseSharedDirectory" onclick="selectSource(this.value);"><B>Import from shared data directory </B><BR/>
                        <DIV id="sharedDirDIV" style="width:600px;max-height:300px;overflow:auto;display:none">
    |;

    &promsMod::browseDirectory_getFiles($promsPath{'shared'},'sharedDirFiles',{fileMatch=>qr/\.(mzXML)$/i});
    
    print qq |
                        </DIV>
                    </TD>
                </TR>
                
                
                <!-- OPENSWATH OPTIONS -->
                <TR>
                    <TH align=right valign=top>OpenSwath workflow options: </TH>
                    <TD bgcolor="$lightColor">
                        <p>
                            iRT file : <INPUT type="file" name="iRTFile" accept=".traML" required><BR/>
                            Windows file for OpenSwath* : <INPUT type="file" name="windowsFile" accept=".txt" required><br/><SMALL>*With a blank header line.</SMALL>
                        </p>
                        
                        <p>
                            Minimal peak width (s) : <INPUT onmouseover="popup('Discard all peaks below this value (-1 means no action).')" onmouseout="popout()" type="text" name="minPeakWidth" size="3" value="20"><br/>
                            RT extraction window (s) : <INPUT onmouseover="popup('Only extract RT around this value (-1 means extract over the whole range, a value of 600 means to extract around +/- 300 s of the expected elution).')" onmouseout="popout()" type="text" name="rtExtractionWin" size="3" value="600"><br/>
                            MZ extraction window : <INPUT onmouseover="popup('Extraction window used')" onmouseout="popout()" type="text" name="mzExtractionWin" size="3" value="0.05">
                            <SELECT name="ppm" id="ppm">
                                <OPTION value="thomson" selected>Thomson</OPTION>
                                <OPTION value="ppm">ppm</OPTION>
                            </SELECT><br/>
                            MZ correction function :
                            <SELECT name="mzCorrectionFunction" id="mzCorrectionFunction" onmouseover="popup('Use the retention time normalization peptide MS2 masses to perform a mass correction (linear, weighted by intensity linear or quadratic) of all spectra.')" onmouseout="popout()">
                                <OPTION value="">-= Select value =-</OPTION>
                                <OPTION value="unweighted_regression">Unweighted regression</OPTION>
                                <OPTION value="weighted_regression">Weighted regression</OPTION>
                                <OPTION value="quadratic_regression">Quadratic regression</OPTION>
                                <OPTION value="weighted_quadratic_regression">Weighted quadratic regression</OPTION>
                                <OPTION value="weighted_quadratic_regression_delta_ppm">Weighted quadratic regression delta ppm</OPTION>
                                <OPTION value="quadratic_regression_delta_ppm">Quadratic regression delta ppm</OPTION>
                            </SELECT>
                        </p>
                        
                        <p>
                            <INPUT type="checkbox" onmouseover="popup('Use MS1 information.')" onmouseout="popout()" name="useMs1Traces" id="useMs1Traces">Use ms1 traces
                            <INPUT type="checkbox" onmouseover="popup('Use for IPF.')" onmouseout="popout()" name="enableUisScoring" id="enableUisScoring">Enable uis scoring<br/>
                        </p>
                    </TD>
                </TR>
                
                <!-- PYPROHET OPTIONS -->
                <TR>
                    <TH align=right valign=top>Pyprophet options: </TH>
                    <TD bgcolor="$lightColor">
                        <span onmouseover="popup('Compute pyprophet score on each run independently.')" onmouseout="popout()"><label><INPUT type="radio" name="pyprophetMethod" value="run-specific" checked required />Run specific </label></span>
                        <span onmouseover="popup('Take all experiment runs into account for the pyprophet score computation.')" onmouseout="popout()"><label><INPUT type="radio" name="pyprophetMethod" value="experiment-wide" required />Whole experiment</label></span>
                    </TD>
                </TR>
                
                <!-- TRIC OPTIONS -->
                <TR id="TRICoption">
                    <TH align=right valign=top>TRIC options: </TH>
                    <TD bgcolor="$lightColor">
                        <SELECT name="TRICMethod">
                            <OPTION value="LocalMST" selected>Local MST*</OPTION>
                            <OPTION value="best_overall">Best Overall</OPTION>
                            <OPTION value="best_cluster_score">Best Cluster Score</OPTION>
                            <OPTION value="global_best_overall">Global Best Overall</OPTION>
                            <OPTION value="global_best_cluster_score">Global Best Cluster Score</OPTION>
                        </SELECT>
                        <SMALL>*Recommanded option</SMALL><BR/>
                        
                        Maximal FDR quality : <INPUT onmouseover="popup('Extension m-score score cutoff, peakgroups of this quality will still be considered for alignment during extension.')" onmouseout="popout()" type="text" name="maxFDRQuality" size="3" value="-1"><br/>
                        FDR cutoff : <INPUT type="text" name="fdrCutoff" size="3" value="0.01"><br/>
                        d-score cutoff : <INPUT type="text" name="dscoreCutoff" size="3" value="1.96"><br/>
                        Maximal RT diff (s) : <INPUT onmouseover="popup('Maximal difference in RT for two aligned features.')" onmouseout="popout()" type="text" name="maxRTDiff" size="3" value="60"><br/>
                    </TD>
                </TR>
                
                <TR id="mergeExperiments" >
                    <TH align=right valign=top>
                        Merge with other experiment:<br/>
                    </TH>
                    <TD id="mergeExperimentsContent" bgcolor="$lightColor">
    |;

    my $sthExperiment = $dbh->prepare("SELECT NAME, ID_EXPERIMENT FROM EXPERIMENT WHERE ID_PROJECT=$projectID");
    $sthExperiment->execute();
    if($sthExperiment->rows > 0) {
        my $mergeExpHTML = "";
        my $sthAnaID = $dbh->prepare("SELECT ID_ANALYSIS,A.NAME FROM SAMPLE S,ANALYSIS A WHERE S.ID_SAMPLE=A.ID_SAMPLE AND S.ID_EXPERIMENT=?");
        
        while (my ($expName, $expID) = $sthExperiment->fetchrow_array) {
            $sthAnaID->execute($expID);
            
            if($sthAnaID->rows > 0) {
                while (my ($anaID, $anaName) = $sthAnaID->fetchrow_array) {
                    my $anaDir = "$promsPath{peptide}/proj_$projectID/ana_$anaID";
                    my $mzXMLFile = "$anaName.mzXML";
                    
                    if (-s "$anaDir/$mzXMLFile" || -s "$anaDir/$mzXMLFile.tar.gz") {
                        $mergeExpHTML .= "<INPUT type=\"checkbox\" name=\"merge\" class='experiment' id=\"allSamples$expID\" value=\"$expID\" onclick=\"ajaxSelectSamples(this.value)\">$expName<BR/><DIV class=\"listSamples\" id=\"listSamples$expID\" style=\"display:none\" ></DIV>\n";
                        last;
                    }
                }
            }
        }

        $mergeExpHTML = "<small><INPUT type='checkbox' id='sameParams' name='sameParams' value='sameParams' onclick='resetExperiments()' checked>All samples were processed with the exact same parameters<BR/></small><br/>\n".$mergeExpHTML if($mergeExpHTML);

        print $mergeExpHTML;
    }
    
    print qq |
                    </TD>
                </TR>
                
                <TR>
                    <TH colspan=2>
                        <input type="submit" name="submit" value="Submit">
                
                        <!-- CLEAR button -->
                        &nbsp;&nbsp;&nbsp;<INPUT type="reset" value="Clear" />
                
                        <!-- CANCEL button -->
                        &nbsp;&nbsp;&nbsp;<INPUT value="Cancel" onclick="cancelAction();" type="button" >
                    </TH>
                </TR>
            </TABLE>
        </FORM>
        
        <DIV id="divDescription" class="clDescriptionCont"></DIV>
        
        <SCRIPT type="text/javascript">setPopup()</SCRIPT>
    |;
}

####>Revision history<####
# 1.4.4 [BUGFIX] Fix import of mzXML files from previous experiments (VS 19/11/19)
# 1.4.3 [CHANGES] Use new job monitoring window opening parameters (VS 18/11/19)
# 1.4.2 [CHANGES] Apply changes according to promsDIA package modifications (VS 18/11/19)
# 1.4.1 [FEATURE] Handles new monitoring system (VS 29/10/19)
# 1.4.0 Overall script redesign to be more modular + separate promsDIA and promsImportManager (VS 21/10/19)
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
