#!/usr/local/bin/perl -w

################################################################################
# promsDIARefonte.pm     1.1.7                                                 #
# Authors: M. Le Picard (Institut Curie)                                       #
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

##################################################################
################### SWATH WORKFLOW BASED ON OPENMS ###############
##################################################################
package DIAWorkflow { 
    use POSIX qw(strftime);
    use strict;
    use DBD::mysql;	# MySQL
    use promsConfig;
    use List::Util qw(min max);
    use Exporter qw(import);
    use File::Basename;
    use File::stat;
    use Data::Dumper;
    
    our @EXPORT = qw(new exportLibrary convertLibrary openswath tric);
    
    # Default parameters ($foceLocal : Force workflow to run locally instead to create jobs on the cluster)
    my ($startTime, $timeStamp, $workDir, $logFile, $errorFile, $fileOut, $openMSPath, $pyprophetPath, $msproteomicstoolsPath, $forceLocalRun); 
    my (@mzXmlFileList, %library, %clusterInfo);
    my $NB_THREAD = 4; # Number of threads to use when parallelization is available
    my $MAX_PARALLEL_JOBS = 10; # Number of jobs allowed to run at the same time
    my $userID = ($ENV{REMOTE_USER}) ? $ENV{REMOTE_USER} : '';
    
    
    sub new {
        ($workDir, $logFile, my $promsPathRef, my $libraryRef, my $mzXmlFileListRef, my $isOnServer) = @_;
        @mzXmlFileList = ($mzXmlFileListRef) ? sort {$a cmp $b} @{$mzXmlFileListRef} : ();
        %library = %{$libraryRef}; # ID, Name, Version
        $library{"paramFile"} = $library{"name"}."_param";

        # Build paths
        %clusterInfo = &promsConfig::getClusterInfo; #('debian'); # default is 'centos'
        #$MAX_PARALLEL_JOBS = $clusterInfo{'maxJobs'} if($clusterInfo{'maxJobs'}); # TODO Uncomment
        $openMSPath = ($clusterInfo{'on'} || $isOnServer) ? $clusterInfo{'path'}{'openms'} : $promsPathRef->{'openms'};
        $pyprophetPath = ($clusterInfo{'on'}  || $isOnServer) ? $clusterInfo{'path'}{'pyprophet'} : $promsPathRef->{'pyprophet'};
        $msproteomicstoolsPath = ($clusterInfo{'on'} || $isOnServer) ? $clusterInfo{'path'}{'msproteomicstools'} : $promsPathRef->{'msproteomicstools'};
        
        $timeStamp = strftime("%Y%m%d%H%M%S", localtime);
        $startTime = strftime("%s", localtime);
        $fileOut = "$workDir/swathWorkflow.out";
        $errorFile = "$workDir/status_$library{ID}_error.out";
    }
    
    sub exportParams {
        my ($exportParamRef, $outFilePath, $step) = @_;
        $outFilePath = "$workDir/".$library{"paramFile"} if(!$outFilePath);
        my %exportParam = %{$exportParamRef};
        
        open(PARAMFILE,">>", $outFilePath) or die ("open: $!");
        
        ## Library export
        if(!$step || $step eq 'libExport') {
            print PARAMFILE "Export for $exportParam{format}\nION_MASS_LIMITS:$exportParam{lmin},$exportParam{lmax}\nION_TYPE:$exportParam{s}\nION_CHARGE:$exportParam{x}\nION_NUMBER_PER_PEPTIDE:$exportParam{o},$exportParam{n}\n";
            print PARAMFILE "MAXIMUM_ERROR_ALLOWED:$exportParam{p}\nTIME_SCALE:$exportParam{t}\nUIS_ORDER:$exportParam{y}\nSWATH_FILE:".basename($exportParam{'w'})."\n";
            print PARAMFILE "FASTA_FILE:$exportParam{f}\n" if($exportParam{"f"}); # Fasta file
            print PARAMFILE "FRAGMENT_MASS_MODIFICATIONS_ALLOWED:$exportParam{g}\n" if($exportParam{"g"}); # Allowed fragment modifications
            print PARAMFILE "MODIFICATION_FILE:".basename($exportParam{"m"})."\n" if($exportParam{"m"}); # Delta mass file
            print PARAMFILE "LABELLING_FILE:".basename($exportParam{"i"})."\n" if($exportParam{"i"}); # Labelling file
            print PARAMFILE "LIBRARY_VERSION:$library{version}\n";
            print PARAMFILE "OTHER:Use theorical mass\n" if($exportParam{"other"} =~ /-e/);
        }
        
        ## Library processing
        if(!$step || $step eq 'libProcessing') {
            print PARAMFILE "PRECURSOR_MZ_THRESHOLD:$exportParam{precursorMZThresh}\n" if($exportParam{"precursorMZThresh"});
            print PARAMFILE "PRODUCT_MZ_THRESHOLD:$exportParam{productMZThresh}\n" if($exportParam{"productMZThresh"});
            print PARAMFILE "PRECURSOR_MIN_MZ:$exportParam{precursorLowerMZ}\n" if($exportParam{"precursorLowerMZ"});
            print PARAMFILE "PRECURSOR_MAX_MZ:$exportParam{precursorUpperMZ}\n" if($exportParam{"precursorUpperMZ"});
            print PARAMFILE "IPF:Use IPF\n" if($exportParam{"IPF"});
            print PARAMFILE "UNIMOD_MOD_FILE:".basename($exportParam{"unimodFile"})."\n" if($exportParam{"unimodFile"});
            print PARAMFILE "MAX_NB_ALT_LOCALIZATIONS:$exportParam{maxNumLocalization}\n" if($exportParam{"maxNumLocalization"});
            print PARAMFILE "\nDesired proteins not found into the library : $exportParam{missingProt}\n" if($exportParam{"missingProt"});
        }
        
        ## OpenSwath workflow
        if(!$step || $step eq 'OpenSwath') {
            print PARAMFILE "IRT_FILE:".basename($exportParam{"irtFile"})."\n" if($exportParam{"irtFile"});
            print PARAMFILE "SWATH_WINDOWS:".basename($exportParam{"swathWindows"})."\n" if($exportParam{"swathWindows"});
            print PARAMFILE "FDR_CUTOFF:$exportParam{fdrCutoff}\n" if($exportParam{"fdrCutoff"});
            print PARAMFILE "MZ_EXTRACTION_WINDOWS:$exportParam{mzExtractionWin}\n" if($exportParam{"mzExtractionWin"});
            print PARAMFILE "MZ_EXTRACTION_WINDOWS_UNIT:$exportParam{ppm}\n" if($exportParam{"ppm"});
            print PARAMFILE "MZ_CORRECTION_FUNCTION:$exportParam{mzCorrectionFunction}\n" if($exportParam{"mzCorrectionFunction"});
            print PARAMFILE "RT_EXTRACTION_WINDOWS:$exportParam{rtExtractionWin}\n" if($exportParam{"rtExtractionWin"});
            print PARAMFILE "MIN_PEAK_WIDTH:$exportParam{minPeakWidth}\n" if($exportParam{"minPeakWidth"});
        }
        
        ## PyProphet
        if(!$step || $step eq 'PyProphet') {
            print PARAMFILE "PYPROPHET_METHOD:$exportParam{method}\n" if($exportParam{"method"});
        }
        
        ## TRIC alignment
        if(!$step || $step eq 'TRIC') {
            system "sed '/FDR_CUTOFF|MAX_RT_DIFF|MAX_FDR_QUALITY|D_SCORE_CUTOFF|TRIC_METHOD/d' $outFilePath";
            print PARAMFILE "FDR_CUTOFF:$exportParam{fdrCutoff}\n" if($exportParam{"fdrCutoff"});
            print PARAMFILE "MAX_RT_DIFF:$exportParam{maxRTDiff}\n" if($exportParam{"maxRTDiff"});
            print PARAMFILE "MAX_FDR_QUALITY:$exportParam{maxFDRQuality}\n" if($exportParam{"maxFDRQuality"});
            print PARAMFILE "D_SCORE_CUTOFF:$exportParam{dscoreCutoff}\n" if($exportParam{"dscoreCutoff"});
            print PARAMFILE "TRIC_METHOD:$exportParam{method}\n" if($exportParam{"method"});
        }
        
        close PARAMFILE;
    }
    
    sub exportLibrary {
        my ($exportParamRef, $pepMod, $forceLocalRun, $protList, $processText, $loadingDivID, $loadingSPAN) = @_;
        my %exportParam = %{$exportParamRef};
        my $libDir = $library{"dir"};
        my $libraryName = $library{"name"};
        my $dbh = &promsConfig::dbConnect('no_user');
        
        $processText = ($processText) ? $processText : '';
        
        ###################################
        #### Processing submitted form ####
        ###################################
        my ($divID, $now, $waitTime, $status, $loading);
        if ($processText) {
            $divID = "document.getElementById('waitDIV')";
            $now = strftime("%s", localtime); # in sec
            $waitTime = strftime("%Hh %Mm %Ss", localtime($now - $startTime - 3600));
            $status = "Updated $waitTime ago";
            print "<BR/><SCRIPT LANGUAGE=\"JavaScript\">$divID.innerHTML=\"\";$divID.innerHTML='$status';</SCRIPT>";
        }
        
        #################################################
        ###  Compute excluded peptides modifications  ###
        #################################################
        my %massAAave=&promsConfig::getMassAAave; 
        my @excludeMod;
        my $sthModInfo=$dbh->prepare("SELECT SLM.SPECIFICITY, MONO_MASS FROM MODIFICATION M, SWATH_LIB_MODIFICATION SLM WHERE SLM.ID_MODIFICATION=M.ID_MODIFICATION AND M.ID_MODIFICATION=? AND ID_SWATH_LIB=$library{ID}");
        if($pepMod) {
            foreach my $selectMod ($pepMod) {
                $sthModInfo->execute($selectMod);
                my ($specificity, $monoMass) = $sthModInfo->fetchrow_array;
                if($specificity) {
                    foreach my $aa (split(//, $specificity)) {
                        if ($aa eq '=') {
                            my $massExcluded = $monoMass + 1.00794;
                            push @excludeMod, sprintf("%0.f", $massExcluded);
                        } else {
                            my $massExcluded = $monoMass + $massAAave{$aa};
                            push @excludeMod, sprintf("%0.f", $massExcluded);
                        }
                    }
                }
            }
        }
        $sthModInfo->finish;
        $dbh->disconnect;
        
        ###############################################
        ### Recovering the protein filtering list   ###
        ###############################################
        my @protList;
        if ($protList) { # Only in exportLibrary.cgi
            open(PROTLIST, "<$protList") or die ("open: $!");
            while (<PROTLIST>) {
                foreach my $protein (split('[;,\s]', $_)){
                    push @protList, $protein;
                }
            }
            close PROTLIST;
        }
        
        ### Print progress
        if ($processText) {
            $processText = "<B>Step 1/2 :</B> Export library for PeakView ... " if ($exportParam{"format"} eq 'peakview');
            $loading = "<progress value=\"0\" max=\"100\" style=\"width:400px\"></progress>";
            print "<SCRIPT LANGAGE=\"JavaScript\">$loadingDivID.innerHTML=\"\";$loadingDivID.innerHTML='$processText<BR/><BR/> $loading';$loadingSPAN.innerHTML='0%';</SCRIPT>";
        }
        
        #######################################
        ###  Delete all N-ter acetylations  ###
        #######################################
        my (%mzHash, @validProt);
        open(LIBFILE, "<", "$libDir/$libraryName.sptxt") or die ("Could not open sptxt lib file at : $libDir/$libraryName.sptxt $!");
        open(OUTLIBFILE, ">", "$workDir/convert$exportParam{format}.sptxt") or die ("Could not write converted sptxt lib file: $workDir/convert$exportParam{format}.sptxt $!");
        my $saveLine;
        my ($match, $matchList) = (0, 0);
        while (my $line = <LIBFILE>) {
            if ($line =~ /^Name:(.*)/) {
                $match = 0;
                if ($matchList) {
                    print OUTLIBFILE $saveLine;
                }
                $matchList = 0;
                $saveLine = '';
                my $modOK = 1;
                if(@excludeMod) {
                    foreach my $modif (@excludeMod) {
                        if($1 =~ /\[$modif\]/) {
                            $modOK = 0;
                            last;
                        }
                    }
                }
                if ($modOK) {
                    $match = 1;
                    print OUTLIBFILE $line unless (@protList);        ## export all the library
                    $saveLine .= $line;
                }
            } else {
                if ($match == 1) {
                    $saveLine .= $line;
                    print OUTLIBFILE $line unless @protList;        ## export all the library
                    if ($line =~ /^PrecursorMZ: (.+)/) {
                        $mzHash{"$1"} = 1;
                    } elsif ($line =~ /^Comment: (.+)/ && @protList) {      ## export just the selected proteins
                        foreach my $protein (@protList) {
                            if ($1 =~ /$protein/) {
                                push @validProt, $protein unless ("@validProt" =~ /$protein/);
                                $matchList = 1;
                                last;
                            }
                        }
                    }
                }
            }     
        }
        close LIBFILE;
        close OUTLIBFILE;
        
        my $missingProt;
        if (@validProt && @protList) {
            foreach my $protein (@protList) {
                if ("@validProt" !~ /$protein/) {
                    $missingProt .= ',' unless (!$missingProt);
                    $missingProt .= $protein;
                }
            }
        }
        

        ########################################
        ### Generate library parameters file ###
        ########################################
        my $paramFile = $library{"paramFile"};
        system("rm $workDir/$paramFile") if(-s "$workDir/$paramFile"); 
        $exportParam{"missingProt"} = $missingProt if($missingProt); # Add missing proteins to parameters if any
        exportParams($exportParamRef, "$workDir/$paramFile", "libExport") if($paramFile);
        
        ########################
        ###  Export library  ###
        ########################
        my $output = "$workDir/sortieExport.txt";
        my $processFile = "$libraryName\_$exportParam{format}.tsv";
        my $finalFile = ($exportParam{format} eq "peakview") ? $libraryName."_peakview.txt" : $processFile;
        my $scriptFile = "$workDir/libExport.sh";
        
        # Build task parameters
        my $exportOptions = "";
        $exportOptions .= " -a $workDir/$processFile";
        $exportOptions .= " -f $exportParam{f} " if($exportParam{"f"}); # Fasta file
        $exportOptions .= " -g $exportParam{g} " if($exportParam{"g"}); # Allowed fragment modifications
        $exportOptions .= " -i $exportParam{i} " if($exportParam{"i"}); # Labelling file
        $exportOptions .= " -k ".$exportParam{"format"} if($exportParam{"format"});
        $exportOptions .= " -l ".$exportParam{"lmin"}.",".$exportParam{"lmax"} if($exportParam{"lmin"} && $exportParam{"lmax"});
        $exportOptions .= " -m ".$exportParam{"m"} if($exportParam{"m"}); # Delta mass file
        $exportOptions .= " -n ".$exportParam{"n"} if($exportParam{"n"});
        $exportOptions .= " -o ".$exportParam{"o"} if($exportParam{"o"});
        $exportOptions .= " -p ".$exportParam{"p"} if($exportParam{"p"});
        $exportOptions .= " -s ".$exportParam{"s"} if($exportParam{"s"});
        $exportOptions .= " -t ".$exportParam{"t"} if($exportParam{"t"});
        $exportOptions .= " -w ".$exportParam{"w"} if($exportParam{"w"});
        $exportOptions .= " -x ".$exportParam{"x"} if($exportParam{"x"});
        $exportOptions .= " -y ".$exportParam{"y"} if($exportParam{"y"});
        $exportOptions .= " ".$exportParam{"other"} if($exportParam{"other"});
        
        my $command = $msproteomicstoolsPath."/spectrast2tsv.py $exportOptions $workDir/convert$exportParam{format}.sptxt > $output 2>&1";
        if ($exportParam{"format"} eq "peakview") {
            $command .= "sed s/TRUE/FALSE/g $workDir/$processFile > $workDir/peakviewfile.txt;";
        }
        
        open (BASH,"+>", $scriptFile);
        print BASH "#!/bin/bash\n";
        print BASH "$command\n";
        close BASH;
        chmod 0775, $scriptFile;

        my $size = (-e "$libDir/$libraryName.sptxt") ? `stat -c "%s" $libDir/$libraryName.sptxt`: 1073741824;
        my $maxMem = max(1, $size/1073741824);
        
        # Run task
        runTask($workDir, $scriptFile, $maxMem, "libExport_$timeStamp", $forceLocalRun);

        #############################################
        ### Waiting for spectrast2tsv to finish   ###
        #############################################
        my $massNumTot = scalar keys %mzHash;
        $massNumTot *= $exportParam{"n"}; # ion number per peptide
        my $wait = 1;
        my $errorTxt = '';
        my $massNumber;
        my $prevNbLine = 0;
        my $waitNb = 0;
        
        while ($wait == 1) {
            sleep 30;
            $now = strftime("%s", localtime); # in sec
            $waitTime = strftime("%Hh %Mm %Ss", localtime($now-$startTime-3600));
            $status = "Updated $waitTime ago";
            
            ## loading process
            $massNumber = `tail -n +2 $workDir/$processFile | cut -f1 | uniq | wc -l` if (-s "$workDir/$processFile"); # Count amount of uniq precursor mass in the processed file
            $massNumber = ($massNumber) ? $massNumber : 1;
            chomp $massNumber;
            
            my $percent = $massNumber/$massNumTot*100;
            $percent = sprintf("%.0f", $percent);
            $percent = '100' if ($percent > 100);
            if ($processText) {
                print "<SCRIPT LANGUAGE=\"JavaScript\">$divID.innerHTML=\"\";$divID.innerHTML='$status';</SCRIPT>";
                $processText = "<B>Step 1/2 :</B> Export library for PeakView ... " if ($exportParam{"format"} eq 'peakview');
                $loading = "<progress value=\"$percent\" max=\"100\" style=\"width:400px\"></progress>";
                print "<SCRIPT LANGAGE=\"JavaScript\">$loadingDivID.innerHTML=\"\";$loadingDivID.innerHTML='$processText<BR/><BR/> $loading';$loadingSPAN.innerHTML='$percent%';</SCRIPT>";
            }
            
            open(EXPORTFILE, ">>$output");
            print EXPORTFILE $percent;
            close EXPORTFILE;
            
            if(-s $output) {
                my $processStatus = `grep -c 'Done.' $output`;
                chomp $processStatus;
                if ($processStatus) { # finished normally
                    $wait = 0;
                } else {
                    # Look for error
                    $processStatus = `grep -c 'Error' $output`;
                    chomp $processStatus;
                    if ($processStatus) {
                        $errorTxt = `tail $output`;
                        $wait = 0;
                    } else {
                        $processStatus = `grep -c 'This modification has not been recognized' $output`;
                        chomp $processStatus;
                        if ($processStatus) {
                            $errorTxt = '<BR/>A modification has not been recognized, select the modification file in the export form (\'File with modifications delta mass\').';
                            $wait = 0;
                        }
                    }
                }
            }
            
            # Parse processed file to check for added lines
            if(-s "$workDir/$processFile" && $waitNb == 2) {
                my @infoLine = split(/\s/, `wc -l $workDir/$processFile`);
                my $nbLine = $infoLine[0];
                $wait = ($nbLine > $prevNbLine) ? 1 : 0;
                $prevNbLine = $nbLine;
                $waitNb = 0;
            }
            
            $waitNb++;
        }

        # Remove unecessary files
        system "rm -f $workDir/PBSerror.txt $workDir/PBS.txt $workDir/torqueID.txt $workDir/ssh_connect.sh $workDir/PBScommand.sh $workDir/sortieExport.txt $scriptFile $workDir/convert$exportParam{format}.sptxt";
        
        if ($errorTxt) {
            if ($processText) {
                print("<BR/><BR/><BR/><BR/><BR/><FONT class='title2' color='#DD0000'>***ERROR: Data preprocessing failed : $errorTxt</FONT><BR/>");
            }
            
            return ('error', 'error', $errorTxt);
        }
        
        #############################
        ###> Add +50 to each iRT <###
        #############################
        if ($exportParam{"format"} eq "peakview") {
            my $numFileLine = `cat $workDir/peakviewfile.txt | wc -l`;
            my $numLine = 1;
            my $i = 1;
            
            open(INFILE, "<", "$workDir/peakviewfile.txt");
            open(OUTFILE, ">" , $finalFile);
            while (my $line=<INFILE>) {
                if ($line=~/^Q1/) {
                    print OUTFILE $line;
                } else {
                    my @lineInfo = split(/\t/, $line);
                    my $RT = $lineInfo[2] + 50;
                    my $iRT = $lineInfo[12] + 50;
                    $lineInfo[2] = $RT;
                    $lineInfo[12] = $iRT;
                    print OUTFILE join("\t", @lineInfo);
                }
                
                if ($i == 300) {
                    ## loading bar
                    my $percent = $numLine/$numFileLine*100;
                    $percent = sprintf("%.0f", $percent);
                    $percent ='100' if ($percent > 100);
                    
                    if($processText) {
                        my $loading = "<progress value=\"$percent\" max=\"100\" style=\"width:400px\"></progress>";
                        print "<SCRIPT LANGAGE=\"JavaScript\">$loadingDivID.innerHTML=\"\";$loadingDivID.innerHTML='<B>Step 2/2 :</B> Conversion for PeakView ... <BR/><BR/> $loading';$loadingSPAN.innerHTML='$percent%';</SCRIPT>";
                    }
                    $i = 0;
                }
                
                $i++;
                $numLine++; 
            }
            
            close OUTFILE;
            close INFILE;
        }

        system "rm -f $output $workDir/PBSerror.txt $workDir/PBS.txt $workDir/torqueID.txt $workDir/*.sh"; # $scriptFile"; # TODO Uncomment $scriptFile
        
        return ($finalFile, $paramFile);
    }
    
    sub convertLibrary {
        my ($convertParamRef, $forceLocalRun) = @_;
        my $targetLibName = $library{"name"};
        my $sourceLibName = $convertParamRef->{"libFile"};

        my $assayGeneratorOptions = " -swath_windows_file ".basename($convertParamRef->{"w"}) if($convertParamRef->{"w"});
        $assayGeneratorOptions .= " -precursor_mz_threshold ".$convertParamRef->{"precursorMZThresh"} if ($convertParamRef->{"precursorMZThresh"});
        $assayGeneratorOptions .= " -product_mz_threshold ".$convertParamRef->{"productMZThresh"} if($convertParamRef->{"productMZThresh"});
        $assayGeneratorOptions .= " -precursor_lower_mz_limit ".$convertParamRef->{"precursorLowerMZ"} if($convertParamRef->{"precursorLowerMZ"});
        $assayGeneratorOptions .= " -precursor_upper_mz_limit ".$convertParamRef->{"precursorUpperMZ"} if($convertParamRef->{"precursorUpperMZ"});
        $assayGeneratorOptions .= " -min_transitions ".$convertParamRef->{"o"} if($convertParamRef->{"o"});
        $assayGeneratorOptions .= " -max_transitions ".$convertParamRef->{"n"} if($convertParamRef->{"n"});
        $assayGeneratorOptions .= " -allowed_fragment_types ".$convertParamRef->{"s"} if($convertParamRef->{"s"});
        $assayGeneratorOptions .= " -allowed_fragment_charges ".$convertParamRef->{"x"} if($convertParamRef->{"x"});
        $assayGeneratorOptions .= " -product_lower_mz_limit ".$convertParamRef->{"lmin"} if($convertParamRef->{"lmin"});
        $assayGeneratorOptions .= " -product_upper_mz_limit ".$convertParamRef->{"lmax"} if($convertParamRef->{"lmax"});
        $assayGeneratorOptions .= " -threads $NB_THREAD";
        $assayGeneratorOptions .= ' -enable_detection_specific_losses' if($convertParamRef->{"IPF"} || $convertParamRef->{"enableSpLosses"});
        $assayGeneratorOptions .= ' -enable_detection_unspecific_losses' if($convertParamRef->{"enableUnspLosses"});
        
        if($convertParamRef->{"IPF"}) {
            $assayGeneratorOptions .= " -enable_ipf";
            $assayGeneratorOptions .= " -unimod_file ".$convertParamRef->{"unimodFile"} if($convertParamRef->{"unimodFile"});
            $assayGeneratorOptions .= " -max_num_alternative_localizations ".$convertParamRef->{"maxNumLocalization"} if($convertParamRef->{"maxNumLocalization"});
            $assayGeneratorOptions .= " ".$convertParamRef->{"others"} if($convertParamRef->{"others"});
            
        }
        
        my $decoyGeneratorOptions = " -product_mz_threshold ".$convertParamRef->{"productMZThresh"} if($convertParamRef->{"productMZThresh"});
        $decoyGeneratorOptions .= " -allowed_fragment_types ".$convertParamRef->{"s"} if($convertParamRef->{"s"});
        $decoyGeneratorOptions .= " -allowed_fragment_charges ".$convertParamRef->{"x"} if($convertParamRef->{"x"});
        $decoyGeneratorOptions .= " -method shuffle";
        $decoyGeneratorOptions .= " -min_decoy_fraction 0.6";
        $decoyGeneratorOptions .= " -shuffle_sequence_identity_threshold 0.8";
        $decoyGeneratorOptions .= " -shuffle_max_attempts 30";
        $decoyGeneratorOptions .= " -decoy_tag DECOY_";
        $decoyGeneratorOptions .= " -threads $NB_THREAD";
        $decoyGeneratorOptions .= ' -enable_detection_specific_losses' if($convertParamRef->{"IPF"} || $convertParamRef->{"enableSpLosses"});
        $decoyGeneratorOptions .= ' -enable_detection_unspecific_losses' if($convertParamRef->{"enableUnspLosses"});
        
        DIAWorkflow::exportParams($convertParamRef, "$workDir/".$library{"paramFile"}, "libProcessing") if($library{"paramFile"});
        
        # Build script file
        my $scriptFile = "$workDir/libConversion.sh";
        open (BASH,"+>", $scriptFile);
        print BASH "#!/bin/bash\n";
        print BASH "cd $workDir;\n";
        print BASH "OPENMS_DATA_PATH=$workDir;\n";
        print BASH "$openMSPath/OpenSwathAssayGenerator -in $sourceLibName -out $targetLibName\_assay.TraML $assayGeneratorOptions >> $fileOut 2>&1 ;\n";
        print BASH "rm $fileOut 2>&1 ;\n";
        print BASH "echo 'Library processing : Generating Decoys' >> $logFile 2>&1 ;\n" if($logFile);
        print BASH "$openMSPath/OpenSwathDecoyGenerator -in $targetLibName\_assay.TraML -out $targetLibName\_assay_DECOY.TraML $decoyGeneratorOptions >> $fileOut 2>&1 ;\n";
        print BASH "echo 'Library processing : Converting results to readable format' >> $logFile 2>&1 ;\n" if($logFile);
        print BASH "$openMSPath/TargetedFileConverter -in $targetLibName\_assay_DECOY.TraML -out $targetLibName\_assay_DECOY.pqp >> $fileOut 2>&1;\n" if(!$convertParamRef->{"targetFormat"} || $convertParamRef->{"targetFormat"} eq 'pqp');
        print BASH "$openMSPath/TargetedFileConverter -in $targetLibName\_assay_DECOY.TraML  -out $targetLibName.tsv >> $fileOut 2>&1;\n" if(!$convertParamRef->{"targetFormat"} || $convertParamRef->{"targetFormat"} eq 'tsv'); # Used to parse lib data
        close BASH;
        chmod 0775, $scriptFile;
        
        my $size = (-e "$workDir/$sourceLibName") ? `stat -c "%s" $workDir/$sourceLibName`: 1073741824;
        my $maxMem = max(10, ($size/1073741824)*10);
        $maxMem = min($maxMem, 60);
        
        # Run task
        runTask($workDir, $scriptFile, $maxMem, "libConversion_$timeStamp", $forceLocalRun, $NB_THREAD);

        # Watch for task completion        
        my $error = "";
        my $nbProcess = 0;
        my $processPattern = 'TargetedFileConverter took';
        my $nbProcessExpected = 1;
        if($convertParamRef->{"targetFormat"}) {
            if($convertParamRef->{"targetFormat"} eq 'TraML') {
                $processPattern = 'OpenSwathDecoyGenerator';
            }
        } elsif($convertParamRef->{"IPF"}) {
            $nbProcessExpected = 2;
        }
        
        my $nbWhile = 0;
        my $maxNbWhile = 48*60*2;
        my ($currentStatus, $currentPercent);
        while ($nbProcess != $nbProcessExpected && !$error) {
            sleep 30 if($nbWhile);
            $currentPercent = 0;
            
            if (-s $fileOut) {
                $error = `grep -i 'error\\\|exception' $fileOut`;
                
                $currentStatus = `grep -aoP "'.+'" $fileOut | tail -n 1`;
                $currentStatus = ($currentStatus) ? substr($currentStatus, 1, -2) : '';
                
                if(`tail -1 $fileOut | tr -d "\n"` !~ /done|took/) {
                    $currentPercent = `grep -aoP "[0-9]{2}[.][0-9]{2}+\\s." $fileOut | tail -n 1`;
                    $currentPercent = ($currentPercent && substr($currentPercent, 6, 1) eq '%') ? substr($currentPercent, 0, 5) : 0;
                }
                
                $nbProcess = `grep -c '$processPattern' $fileOut`;
            }
            
            if (-s "$workDir/PBSerror.txt") {
                $error = $clusterInfo{'checkError'}->("$workDir/PBSerror.txt");
            }
            
            my $tot = $nbWhile*30;
            my $main = sprintf("%0f", $tot/60);
            my $sec = $tot % 60;
            $currentStatus = "" if(!$currentStatus);
            if($logFile && ($currentPercent || $currentStatus)) {
                my $statusTxt = ($currentPercent) ? "($currentPercent%)" : "";
                $statusTxt  = ($currentStatus) ? ": $currentStatus $statusTxt" : $statusTxt;
                printToFile($logFile, "Library processing $statusTxt", 1);
            }
            print("Library processing $currentStatus => $error | $nbProcess ($tot sec total)<br/>");
            
            $nbWhile++;
            if ($nbWhile > $maxNbWhile) {
                $error = "Aborting: Library conversion is taking too long.";
            }
        }
        
        system "rm -f $workDir/PBSerror.txt $workDir/PBS.txt $workDir/torqueID.txt $workDir/*.sh" if(!$error); # $scriptFile"; # TODO Uncomment $scriptFile
        
        return ($error) ? $error : "Done";
    }
    
    sub openswath {
        my ($openSwathParamRef, $forceLocalRun) = @_;
        my ($MIN_MEMORY, $MAX_MEMORY, $MAX_TRIES) = (10, 30, 3);
        
        # OpensSwath options
        my $openSwathOptions = " -sort_swath_maps";
        $openSwathOptions .= " -batchSize 1000";
        $openSwathOptions .= " -tr ".$openSwathParamRef->{"decoyLibrary"};
        $openSwathOptions .= " -tr_irt ".$openSwathParamRef->{"irtFile"};
        $openSwathOptions .= " -swath_windows_file ".$openSwathParamRef->{"swathWindows"};
        $openSwathOptions .= " -mz_extraction_window ".$openSwathParamRef->{"mzExtractionWin"} if($openSwathParamRef->{"mzExtractionWin"});
        $openSwathOptions .= " -ppm " if($openSwathParamRef->{"ppm"} && $openSwathParamRef->{"ppm"} eq "ppm");
        $openSwathOptions .= " -mz_correction_function ".$openSwathParamRef->{"mzCorrectionFunction"} if($openSwathParamRef->{"mzCorrectionFunction"});
        $openSwathOptions .= " -rt_extraction_window ".$openSwathParamRef->{"rtExtractionWin"} if($openSwathParamRef->{"rtExtractionWin"});
        $openSwathOptions .= " -Scoring:TransitionGroupPicker:min_peak_width ".$openSwathParamRef->{"minPeakWidth"} if($openSwathParamRef->{"minPeakWidth"});
        #$openSwathOptions .= " -use_ms1_traces" if($openSwathParamRef->{"useMs1Traces"});
        #$openSwathOptions .= " -Scoring:TransitionGroupPicker:background_subtraction exact"; # TODO Check its impact on quantification values
        
        $openSwathOptions .= " -enable_uis_scoring" if($openSwathParamRef->{"IPF"} || $openSwathParamRef->{"enableUisScoring"});
        $openSwathOptions .= " -Scoring:Scores:use_uis_scores" if($openSwathParamRef->{"IPF"} || $openSwathParamRef->{"enableUisScoring"});
        #$openSwathOptions .= " -Scoring:Scores:use_ms1_mi" if($openSwathParamRef->{"IPF"}); # Check this parameter
        #$openSwathOptions .= " -Scoring:Scores:use_mi_score" if($openSwathParamRef->{"IPF"}); # Check this parameter
        #$openSwathOptions .= " -Scoring:Scores:use_total_mi_score" if($openSwathParamRef->{"IPF"}); # Check this parameter
        $openSwathOptions .= " -use_ms1_traces" if ($openSwathParamRef->{"IPF"} || $openSwathParamRef->{"useMs1Traces"});
        $openSwathOptions .= " -threads $NB_THREAD" if($NB_THREAD > 1);
        #$openSwathOptions .= " -readOptions cacheWorkingInMemory -tempDirectory $runDir/"; # TODO Check if it's really needed
        
        DIAWorkflow::exportParams($openSwathParamRef, "$workDir/".$library{"paramFile"}, "OpenSwath");
        
        ## Apply OpenSwath workflow on search files
        my (%processedFiles, %maxLaunch);
        
        # Watch for processing completion
        my $nbWhile = 0;
        my $maxNbWhile = 48*60*2;
        my $nProcessRunning = 0;
        my ($error, $PBSerror);
        
        while (scalar keys %processedFiles != scalar @mzXmlFileList && !$error && !$PBSerror) {
            sleep 30 if($nbWhile);
            my $pastTime = $nbWhile*30; # TODO Remove
            my $maxProcessPct = 0;
            
            for (my $i=1; $i <= scalar @mzXmlFileList; $i++) {
                my $mzXMLFile = $mzXmlFileList[$i-1];
                (my $OSWout = $mzXMLFile) =~ s/.mzXML$/.osw/;
                (my $mzXMLout = $mzXMLFile) =~ s/.mzXML$/.tsv/;
                (my $chromOut = $mzXMLFile) =~ s/.mzXML$/.mzML/;
                (my $anaName = $mzXMLFile) =~ s/.mzXML$//;
                
                my $runDir = $workDir.'/OpenSwath_'.$i;
                my $scriptFile = "$runDir/OpenSwath.sh";
                my $outputFile = "$runDir/openSwathOUT.txt";
    
                next if($processedFiles{$i});
                
                if(!$maxLaunch{$i}) { # Job not started yet
                    if(-s "$workDir/$OSWout") {
                        $processedFiles{$i} = 1;
                        next;
                    }
                    
                    if($nProcessRunning < $MAX_PARALLEL_JOBS) { # Start job
                        mkdir $runDir unless -e $runDir;
                        
                        open(SCRIPT, "+>", $scriptFile);
                        print SCRIPT "#!/bin/bash\n";
                        print SCRIPT "OPENMS_DATA_PATH=/usr/share/OpenMS;\n";
                        # ADD -out_chrom $workDir/compressed_$chromOut FOR EXPORTING CHROMATOGRAM
                        print SCRIPT "$openMSPath/OpenSwathWorkflow -in $workDir/$mzXMLFile -out_osw $workDir/$OSWout $openSwathOptions >> $outputFile 2>&1;\n";
                        #print SCRIPT "$openMSPath/FileConverter -in $workDir/compressed_$chromOut -in_type 'mzML' -out $workDir/$chromOut >> $outputFile 2>&1;\n"; # TODO uncomment;
                        #print SCRIPT "rm -f $workDir/compressed_$chromOut;\n"; # TODO uncomment;
                        print SCRIPT "echo 'OpenSwathWorkflow Done.' >> $outputFile 2>&1;\n";
                        close SCRIPT;
                        chmod 0775, $scriptFile;
                        
                        my $mzxmlSize = `stat -c "%s" $workDir/$mzXMLFile`;
                        my $maxMem = max($MIN_MEMORY, ($mzxmlSize/1073741824)*5);
                        $maxMem = min($maxMem, $MAX_MEMORY);
                        
                        # Run Task
                        $maxLaunch{$i} = 1;
                        runTask($runDir, $scriptFile, $maxMem, "openSwath_Workflow_$anaName\_$timeStamp\_".$maxLaunch{$i}, $forceLocalRun, $NB_THREAD);
                        $nProcessRunning++;
                        sleep 1;
                    }
                } else { # Check for process status
                    if (-s $outputFile) {
                        $error = `grep -i 'error\\\|exception' $outputFile`;
                        if($error) {
                            printToFile($logFile, "ERROR for $mzXMLFile in OpenSwath step : $error", 1) if($logFile);
                            
                            if($maxLaunch{$i} && $maxLaunch{$i} > $MAX_TRIES) {
                                printToFile($logFile, "ERROR for $mzXMLFile in OpenSwath step : limit of retries reached !", 1) if($logFile);
                                return $error;
                            } elsif($error !~ /does not have a corresponding chromatogram/) { # Error : could also be "OpenMS FATAL ERROR" and work with retry 
                                printToFile($logFile, "ERROR for $mzXMLFile : retrying !", 1) if($logFile);
                                system "rm -f $runDir/PBSerror.txt $runDir/PBS.txt $runDir/torqueID.txt $outputFile $workDir/$OSWout";
                                
                                $error = '';
                                $maxLaunch{$i}++;
                                
                                my $mzxmlSize = `stat -c "%s" $workDir/$mzXMLFile`;
                                my $maxMem = max($MIN_MEMORY, ($mzxmlSize/1073741824)*5);
                                $maxMem = min($maxMem, $MAX_MEMORY);
                                runTask($runDir, $scriptFile, $maxMem, "openSwath_Workflow_$anaName\_$timeStamp\_".$maxLaunch{$i}, $forceLocalRun, $NB_THREAD);
                                sleep 1;
                            } else {
                                return $error;
                            }
                        }
                        
                        # Looking for end of task in log file
                        if(-s "$runDir/PBS.txt") {
                            my $pastTimeM = $pastTime/60;
                            $processedFiles{$i} = 1;
                            my $tries = ($maxLaunch{$i}) ? $maxLaunch{$i}+1 : 1;
                            system "cat $runDir/openSwathOUT.txt >> $fileOut";
                            system "rm -rf  $runDir $workDir/*.pdf";
                            $nProcessRunning--;
                        } else {
                            my $currentProcessPct = `grep -aoP '[0-9.]+\\s\\%(?= \\s+Thread)' $runDir/openSwathOUT.txt | tail -1 | cut -d' ' -f1`;
                            $currentProcessPct = ($currentProcessPct) ? sprintf '%.2f', substr($currentProcessPct, 0, 5) : '0.00';
                            $maxProcessPct = $currentProcessPct if($currentProcessPct > $maxProcessPct);
                        }
                    }
                
                    # Looking for job error
                    if (-s "$runDir/PBSerror.txt") {
                        $PBSerror = $clusterInfo{'checkError'}->("$runDir/PBSerror.txt");
                        if($PBSerror && $maxLaunch{$i} <= $MAX_TRIES) {
                            $PBSerror = substr($PBSerror, 0, -1);
                            printToFile($logFile, "ERROR for $mzXMLFile ($PBSerror) : retrying !", 1) if($logFile);
                            system "rm -f $runDir/PBSerror.txt $runDir/PBS.txt $runDir/torqueID.txt $outputFile $workDir/$OSWout";
                            
                            $PBSerror = '';
                            $maxLaunch{$i}++;
                            
                            my $mzxmlSize = `stat -c "%s" $workDir/$mzXMLFile`;
                            my $maxMem = max($MIN_MEMORY, ($mzxmlSize/1073741824)*5);
                            $maxMem = min($maxMem, $MAX_MEMORY);
                            runTask($runDir, $scriptFile, $maxMem, "openSwath_Workflow_$anaName\_$timeStamp\_".$maxLaunch{$i}, $forceLocalRun, $NB_THREAD);
                            sleep 1;
                        }
                    }
                }
            }

            my $nFilesToProcces = scalar @mzXmlFileList;
            my $nFilesProcessed = scalar keys %processedFiles;
            
            if($logFile) {
                printToFile($logFile, "OpenSwath : $nFilesProcessed/$nFilesToProcces analysis processed ($maxProcessPct%)", 1);
            }
            
            $nbWhile++;
            if ($nbWhile > $maxNbWhile) {
                $error = "Aborting: OpenSwath processing is taking too long.";
            }
        }
        
        system "rm -fr OpenSwath_*/ *.pdf" if(!$error && !$PBSerror);
        
        return ($error) ? $error : ($PBSerror) ? $PBSerror : "Done";
    }
    

    sub pyprophet {
        my ($pyprophetParamRef, $allowErrorUntilEnd, $forceLocalRun) = @_;
            
         # PyProphet options
        my $pyProphetScoringOptions = ($pyprophetParamRef->{"IPF"}) ? " --ipf_max_peakgroup_rank 3" : "";
        $pyProphetScoringOptions .= " --threads $NB_THREAD" if($NB_THREAD > 1);
        $pyProphetScoringOptions .= " --tric_chromprob";
        my $pyProphetExportOptions = ($pyprophetParamRef->{"IPF"}) ? " --ipf peptidoform" : " --ipf disable";
        
        DIAWorkflow::exportParams($pyprophetParamRef, "$workDir/".$library{"paramFile"}, "PyProphet");
        
        # Watch for processing completion
        my ($error, $PBSerror, %processedFiles, %launched, $processTxt);
        my $nbWhile = 0;
        my $maxNbWhile = 48*60*2;
        my $nFilesToProcess = ($pyprophetParamRef->{"method"} eq "run-specific") ? scalar @mzXmlFileList : 1;
        my $nProcessRunning = 0;
        my $hasError = 0;
        
        while (scalar keys %processedFiles != $nFilesToProcess && !$error && !$PBSerror) {
            sleep 30 if($nbWhile);
            my $maxProcessPct = 0;
            $processTxt = '';
            
            if($pyprophetParamRef->{"method"} eq "run-specific") {
                for (my $i=1; $i <= scalar @mzXmlFileList; $i++) {
                    my $mzXMLFile = $mzXmlFileList[$i-1];
                    (my $anaName = $mzXMLFile) =~ s/.mzXML$//;
                    (my $pyprophetFile = $mzXMLFile) = "$mzXMLFile.tsv";
                    my $runDir = $workDir.'/PyProphet_'.$anaName;
                    my $pyprophetOut = "$runDir/pyprophet.out";
                    
                    next if($processedFiles{$i});
                    
                    if(!$launched{$i}) { # Start job
                        if(-s "$workDir/$pyprophetFile") {
                            $processedFiles{$i} = 1;
                            next;
                        }
                        
                        if($nProcessRunning < $MAX_PARALLEL_JOBS) {
                        
                            (my $OSWout = $mzXMLFile) = "$anaName.osw";
                            my $scriptFile = "$runDir/pyprophet_$anaName.sh";
                            
                            mkdir $runDir unless -e $runDir;
                            
                            open(SCRIPT, "+>", $scriptFile);
                            print SCRIPT "#!/bin/bash\n";
                            print SCRIPT "cd $workDir;\n";
                            
                            print SCRIPT "echo 'Computing score of $OSWout ...' >> $pyprophetOut 2>&1 ;\n";
                            if($pyprophetParamRef->{"IPF"}) {
                                print SCRIPT "$pyprophetPath/pyprophet score --in=$OSWout --level=ms1 $pyProphetScoringOptions >> $pyprophetOut 2>&1 ;\n";
                                print SCRIPT "$pyprophetPath/pyprophet score --in=$OSWout --level=ms2 $pyProphetScoringOptions >> $pyprophetOut 2>&1 ;\n";
                                print SCRIPT "$pyprophetPath/pyprophet score --in=$OSWout --level=transition $pyProphetScoringOptions >> $pyprophetOut 2>&1 ;\n";
                            } else {
                                print SCRIPT "$pyprophetPath/pyprophet score --in=$OSWout --level=ms1ms2 $pyProphetScoringOptions >> $pyprophetOut 2>&1 ;\n";
                            }
                            
                            print SCRIPT "echo 'Computing IPF score of $OSWout ...' >> $pyprophetOut 2>&1 ;\n" if($pyprophetParamRef->{"IPF"});
                            print SCRIPT "$pyprophetPath/pyprophet ipf --in=$OSWout >> $pyprophetOut 2>&1 ;\n" if($pyprophetParamRef->{"IPF"});
                            
                            print SCRIPT "echo 'Estimating peptides reliability for $OSWout ...' >> $pyprophetOut 2>&1 ;\n";
                            print SCRIPT "$pyprophetPath/pyprophet peptide --in=$OSWout --context=run-specific >> $pyprophetOut 2>&1 ;\n";
                            
                            print SCRIPT "echo 'Estimating proteins reliability for $OSWout ...' >> $pyprophetOut 2>&1 ;\n";
                            print SCRIPT "$pyprophetPath/pyprophet protein --in=$OSWout --context=run-specific >> $pyprophetOut 2>&1 ;\n";
                            
                            print SCRIPT "echo 'Export PyProphet results as $pyprophetFile ...' >> $pyprophetOut 2>&1 ;\n";
                            print SCRIPT "$pyprophetPath/pyprophet export --in=$OSWout $pyProphetExportOptions >> $pyprophetOut 2>&1 ;\n";
                            
                            print SCRIPT "echo 'PyProphet Done.' >> $pyprophetOut 2>&1 ;\n";
                            print SCRIPT "cat $pyprophetOut >> $fileOut;\n";
                            close SCRIPT;
                            chmod 0775, $scriptFile;
                            
                            # Run task
                            my $mzxmlSize = `stat -c "%s" $workDir/$mzXMLFile`;
                            my $maxMem = max(5, ($mzxmlSize/1073741824)*5);
                            $maxMem = min($maxMem, 20);
            
                            runTask($runDir, $scriptFile, $maxMem, "pyprophet\_$timeStamp\_$anaName", $forceLocalRun);
                            $launched{$i} = 1;
                            $nProcessRunning++;
                            sleep 1;
                        }
                    } else { # Monitor job status
                        if (-s "$runDir/PBSerror.txt") { # Check for process errors
                            $PBSerror = $clusterInfo{'checkError'}->("$runDir/PBSerror.txt");
                            last if($PBSerror);
                        }
                        
                        print("$pyprophetOut");
                        if(-s $pyprophetOut) {
                            $error = `grep -i 'error\\\|exception' $pyprophetOut`;
                            if($error) {
                                printToFile($logFile, "An error was thrown for $anaName : $error", 1) if($logFile);        
                                if($error !~ /IOError|does not have a corresponding chromatogram/) { # IOError => Could not write pdf file
                                    system "rm -rf $runDir $workDir/$anaName.osw $workDir/*.pdf";
                                    last if(!$allowErrorUntilEnd);
                                } else {
                                    $error = '';
                                }
                            }
                        }
                        
                        # Looking for end of task in log file
                        if(-s "$runDir/PBS.txt") {
                            system "rm -rf $runDir $workDir/*.pdf";
                            $processedFiles{$i} = 1;
                            $nProcessRunning--;
                        }
                    }
                }
                
                # Compute nb processed files
                my $nFilesProcessed = scalar keys %processedFiles;
                $processTxt = " : $nFilesProcessed/$nFilesToProcess analysis processed";
                
            } elsif($pyprophetParamRef->{"method"} eq "experiment-wide") {
                if(!$launched{1}) { # Start job
                    my $scriptFile = "$workDir/pyprophet.sh";
                    
                    open(SCRIPT, "+>", $scriptFile);
                    print SCRIPT "#!/bin/bash\n";
                    print SCRIPT "cd $workDir;\n";
                    print SCRIPT "echo 'Error estimation : Merging all osw files ...' >> $logFile 2>&1 ;\n" if($logFile);
                    print SCRIPT "$pyprophetPath/pyprophet merge --out=model.oswm --subsample_ratio=0.1 *.osw >> $fileOut 2>&1 ;\n";
                    
                    print SCRIPT "echo 'Error estimation : Computing scores of merged files entries ...' >> $logFile 2>&1 ;\n" if($logFile);
                    if($pyprophetParamRef->{"IPF"}) {
                        print SCRIPT "$pyprophetPath/pyprophet score --in=model.oswm --level=ms1 $pyProphetScoringOptions >> $fileOut 2>&1 ;\n";
                        print SCRIPT "$pyprophetPath/pyprophet score --in=model.oswm --level=ms2 $pyProphetScoringOptions >> $fileOut 2>&1 ;\n";
                        print SCRIPT "$pyprophetPath/pyprophet score --in=model.oswm --level=transition $pyProphetScoringOptions >> $fileOut 2>&1 ;\n";
                    } else {
                        print SCRIPT "$pyprophetPath/pyprophet score --in=model.oswm --level=ms1ms2 $pyProphetScoringOptions >> $fileOut 2>&1 ;\n";
                    }
                    
                    print SCRIPT "echo 'Error estimation : Apply score to individual runs (ms1) ...' >> $logFile 2>&1 ;\n" if($logFile);
                    if($pyprophetParamRef->{"IPF"}) {
                        print SCRIPT "for run in *.osw \ndo \n$pyprophetPath/pyprophet score --in=\$run --apply_weights=model.oswm --level=ms1 >> $fileOut \ndone 2>&1 ;\n";
                        print SCRIPT "for run in *.osw \ndo \n$pyprophetPath/pyprophet score --in=\$run --apply_weights=model.oswm --level=ms2 >> $fileOut \ndone 2>&1 ;\n";
                        print SCRIPT "for run in *.osw \ndo \n$pyprophetPath/pyprophet score --in=\$run --apply_weights=model.oswm --level=transition >> $fileOut \ndone 2>&1 ;\n";
                        
                        print SCRIPT "echo 'Error estimation : Apply IPF to individual runs ...' >> $logFile 2>&1 ;\n" if($logFile);
                        print SCRIPT "for run in *.osw \ndo \n$pyprophetPath/pyprophet ipf --in \$run \ndone >> $fileOut 2>&1 ;\n";
                    } else {
                        print SCRIPT "for run in *.osw \ndo \n$pyprophetPath/pyprophet score --in=\$run --apply_weights=model.oswm --level=ms1ms2 >> $fileOut \ndone 2>&1 ;\n";
                    }
                    
                    print SCRIPT "echo 'Error estimation : Generate tidy files from score results ...' >> $logFile 2>&1 ;\n" if($logFile);
                    print SCRIPT "for run in *.osw \ndo \n$pyprophetPath/pyprophet reduce --in=\$run --out=\${run}r >> $fileOut \ndone 2>&1 ;\n";
                    
                    print SCRIPT "echo 'Error estimation : Merge tidy files based on previous model ...' >> $logFile 2>&1 ;\n" if($logFile);
                    print SCRIPT "$pyprophetPath/pyprophet merge --template=model.oswm --out=model_global.oswm *.oswr >> $fileOut 2>&1 ;\n";
                    
                    print SCRIPT "echo 'Error estimation : Infer peptides in experiment-wide context ...' >> $logFile 2>&1 ;\n" if($logFile);
                    print SCRIPT "$pyprophetPath/pyprophet peptide --context=experiment-wide --in=model_global.oswm >> $fileOut 2>&1 ;\n";
                    
                    print SCRIPT "echo 'Error estimation : Infer proteins in experiment-wide context ...' >> $logFile 2>&1 ;\n" if($logFile);
                    print SCRIPT "$pyprophetPath/pyprophet protein --context=experiment-wide --in=model_global.oswm >> $fileOut 2>&1 ;\n";
                    
                    print SCRIPT "echo 'Error estimation : Backpropagate results on each run ...' >> $logFile 2>&1 ;\n" if($logFile);
                    print SCRIPT "for run in *.osw \ndo \n$pyprophetPath/pyprophet backpropagate --in=\$run --apply_scores=model_global.oswm >> $fileOut \ndone 2>&1 ;\n";
                    
                    print SCRIPT "echo 'Error estimation : Exporting results ...' >> $logFile 2>&1 ;\n" if($logFile);
                    print SCRIPT "for run in *.osw \ndo \n$pyprophetPath/pyprophet export --in \$run $pyProphetExportOptions >> $fileOut \ndone 2>&1 ;\n";
                    
                    print SCRIPT "rm -f model.oswm *.oswr >> $fileOut 2>&1 ;\n";
                    
                    print SCRIPT "echo 'PyProphet Done.' >> $fileOut 2>&1 ;\n";
                    close SCRIPT;
                    chmod 0775, $scriptFile;
                    
                    # Run task
                    my $totalSize = 0;
                    for (my $i=1; $i <= @mzXmlFileList; $i++) {
                        my $mzXMLFile = $mzXmlFileList[$i-1];
                        my $mzxmlSize = `stat -c "%s" $workDir/$mzXMLFile`;
                        $totalSize = $totalSize + $mzxmlSize;
                    }
                    my $maxMem = max(5, ($totalSize/1073741824)*5);
                    $maxMem = min($maxMem, 100);
                    runTask($workDir, $scriptFile, $maxMem, "pyprophet\_$timeStamp", $forceLocalRun);
                    $launched{1} = 1;
                    sleep 1;
                } else { # Monitor job completion status
                    # Check for process errors
                    if(-s "$workDir/PBSerror.txt") {
                        $PBSerror = $clusterInfo{'checkError'}->("$workDir/PBSerror.txt");
                    }
                    
                    $error = `grep -i 'error\\\|exception' $fileOut`;
                    if($error) {
                        printToFile($logFile, "An error was thrown for : $error", 1) if($logFile);
                        
                        if($error !~ /IOError|does not have a corresponding chromatogram/) {
                            system "rm -rf $workDir/*.osw $workDir/*.sh $workDir/*.sh $workDir/PBS*.txt $workDir/torqueID.txt $workDir/*.pdf $workDir/*.mzXML.tsv";
                        } else {# IOError => Could not write pdf file || KeyError => Error happenning during protein error estimation (Problem ?)
                            $error = '';
                        }
                    }
                    
                    # Looking for end of task in log file
                    if(-s "$workDir/PBS.txt") {
                        system "rm -rf $workDir/*.sh $workDir/PBS*.txt $workDir/torqueID.txt $workDir/*.pdf";
                        $processedFiles{1} = 1;
                    }
                }
            }
            
            if($logFile) {
                printToFile($logFile, "Error estimation computing $processTxt", 1);
            }
            
            my $pastTime = $nbWhile*30;
            print("Pyprophet : $pastTime sec elapsed (Jobs running : $nProcessRunning/$MAX_PARALLEL_JOBS");
            if ($nbWhile > $maxNbWhile) {
                $error = "Aborting: PyProphet is taking too long.";
            }
            $nbWhile++;
        }

        return ($error) ? $error : ($PBSerror) ? $PBSerror : "Done";
    }
    
    sub tric {
        my ($tricParamRef, $forceLocalRun) = @_;
        my @pyprophetFilesOut = @{$tricParamRef->{"inputFiles"}};
        
        # Compute task parameters
        my $scriptFile = "$workDir/TRICProcess.sh"; 
        my $format = (!$tricParamRef->{"format"} && $tricParamRef->{"format"} eq "pkv") ? "peakview" : "openswath";
        my $tricOptions;
        
        $tricOptions = " --method ".$tricParamRef->{"method"};
        if ($tricParamRef->{"method"} eq 'LocalMST') {
            $tricOptions .= " --realign_method lowess";
            $tricOptions .= " --mst:useRTCorrection True";
            $tricOptions .= " --mst:Stdev_multiplier 3.0";
        } else {
            $tricOptions .= " --realign_method diRT";
        }
        
        $tricOptions .= " --file_format ".$format;
        $tricOptions .= " --alignment_score 0.0001";
        $tricOptions .= " --fdr_cutoff ".$tricParamRef->{"fdrCutoff"} if($tricParamRef->{"fdrCutoff"});
        $tricOptions .= " --max_rt_diff ".$tricParamRef->{"maxRTDiff"} if($tricParamRef->{"maxRTDiff"});
        $tricOptions .= " --max_fdr_quality ".$tricParamRef->{"maxFDRQuality"} if($tricParamRef->{"maxFDRQuality"});
        $tricOptions .= " --use_dscore_filter --dscore_cutoff ".$tricParamRef->{"dscoreCutoff"} if($tricParamRef->{"dscoreCutoff"});
        $tricOptions .= " --disable_isotopic_grouping";
        
        DIAWorkflow::exportParams($tricParamRef, "$workDir/".$library{"paramFile"}, "TRIC");
        
        open(SCRIPT, "+>", $scriptFile);
        print SCRIPT "#!/bin/bash\n";
        print SCRIPT "cd $workDir;\n";
        print SCRIPT "$msproteomicstoolsPath/feature_alignment.py --in @pyprophetFilesOut --out feature_alignment_fileout.tsv $tricOptions >> $fileOut 2>&1\n";
        close SCRIPT;
        chmod 0775, $scriptFile;
        
        # Compute memory size to ask for TRIC, based on all feature_alignment.py infiles
        my $maxMem = max(30, scalar @pyprophetFilesOut / 2);
        $maxMem = 200 if($maxMem > 200);
        
        # Run task
        runTask($workDir, $scriptFile, $maxMem, "TRIC_$timeStamp", $forceLocalRun);
        
        # Watch for task completion
        my ($error, $PBSerror, $endProcess) = ('', '', '');
        my $nbWhile = 0;
        my $maxNbWhile = 48*60*2;
        while (!-s "$workDir/feature_alignment_fileout.tsv" && !$PBSerror && !$error) {
            sleep 30 if ($nbWhile);
            
            if (-s "$workDir/$fileOut") {
                $error = `grep -ia 'error\\\|exception' $workDir/$fileOut`;
            }
            
            if (-s "$workDir/PBSerror.txt") {
                $PBSerror = $clusterInfo{'checkError'}->("$workDir/PBSerror.txt");
            }
            
            my $tot = $nbWhile*30;
            my $main = sprintf("%0f", $tot/60);
            my $sec = $tot % 60;
            my $exist = (-s "$workDir/feature_alignment_fileout.tsv") ? "Exist" : "Not created";
            print("Waiting 30 sec more ($tot sec total) => $error | $PBSerror | $exist<br/>");
            
            $nbWhile++;
            if ($nbWhile > $maxNbWhile) {
                $error = "Aborting: TRIC is taking too long.";
            }
        }
        
        # Delete unnecessary files
        system "rm -f $workDir/*.tr $workDir/*.sh $workDir/PBS*.txt $workDir/torqueID.txt";
        
        if ($PBSerror) {
            if ($PBSerror =~ /Exception: No data available for alignment 0_0 vs 0_1/) {
                return "Error TRIC";
            } else {
                return $PBSerror;
            }
        } elsif($error) {
            return $error;
        } else {
            return "Done";
        }
    }
    
    sub runTask {
        my ($runDir, $scriptFile, $maxMem, $jobName, $forceLocalRun, $nCores) = @_;
        my ($jobID);
        ($jobID = $workDir) =~ s/.+\///g;
        my $dbh = &promsConfig::dbConnect('no_user');
        
        if ($clusterInfo{'on'} && !$forceLocalRun) { # Run on Cluster
            my %jobParams = (
                maxMem     => sprintf("%.0f", $maxMem).'Gb',
                numCPUs    => ($nCores) ? $nCores : 1,
                maxHours   => 48,
                jobName    => $jobName,
                outFile    => 'PBS.txt',
                errorFile  => 'PBSerror.txt',
                jobEndFlag => "End_$jobName",
                noWatch    => '1',
            );
            my ($pbsError, $pbsErrorFile, $jobClusterID) = $clusterInfo{'runJob'}->($runDir, $scriptFile, \%jobParams);
            
            $dbh->do("UPDATE JOB_HISTORY SET ID_JOB_CLUSTER = CONCAT(ID_JOB_CLUSTER, ';C".$jobClusterID."') WHERE ID_JOB='$jobID'");
            $dbh->commit;
            
        } else { # Run locally
            my $childConvert = fork;
            unless ($childConvert) {
                open STDOUT, ">$runDir/std.out" or die "Can't open $runDir/std.out: $!";
                open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
                open STDERR, ">$runDir/std.err" or die "Can't open $runDir/std.err: $!";
                
                $dbh->do("UPDATE JOB_HISTORY SET ID_JOB_CLUSTER = CONCAT(ID_JOB_CLUSTER, ';L".$$."') WHERE ID_JOB='$jobID'");
                $dbh->commit;
                
                system "bash $scriptFile";
                system "echo 'End_$jobName' > $runDir/PBS.txt";
            }
        }
        $dbh->disconnect;
    }
    
    sub printToFile {
        my ($filePath, $content, $append) = @_;
        my $openingMode = ($append) ? ">>" : ">";
        open(FILE, $openingMode, $filePath) || die "Error while opening $filePath\n";
        print FILE $content."\n";
        close FILE;
    }
}
1;

##########################################################
################### TANDEM XML PARSING ###################
##########################################################
package XTandemXMLHandler; {    ## for .xml and .tandem
    my (%infoRank,%matches,%protListNoSeq,%query2matches,%matchSpectrum);
    my (%modificationList,%bestScoreTarget,%bestScoreDecoy);
    my ($minScore,$maxRank,$dataPath,$dbID,$modifNTer,$type,$dbh,$msFilename,$analysisID,%info,%rankSeq,%prot,$refProtSeq,$refQueryInfo,$refProtDes,$refProtOrg,$refProtLength,$refCharge,$refMassObs,$refMassExp,$refElutionTime,$refNumValid,$refMatchList,$refMaxProtMatch,$refMaxProtScore,$refMaxQueryScore,$refAnalysisMods,$refRankProtMatch,$refProtMW,$refProtDbRank);
    my $i=1;
    
    sub new {
        ($type,$msFilename,$dataPath,$modifNTer,$dbh,$dbID,$minScore,$analysisID,$maxRank,$refProtSeq,$refQueryInfo,$refProtDes,$refProtOrg,$refProtLength,$refCharge,$refMassObs,$refMassExp,$refElutionTime,$refNumValid,$refMatchList,$refMaxProtMatch,$refMaxProtScore,$refMaxQueryScore,$refAnalysisMods,$refRankProtMatch,$refProtMW,$refProtDbRank) = @_;
        my $self = bless ({}, $type);
        
        ###>Fetching parsing rules
        my ($parseRules)=$dbh->selectrow_array("SELECT PARSE_RULES FROM DATABANK,DATABANK_TYPE WHERE DATABANK.ID_DBTYPE=DATABANK_TYPE.ID_DBTYPE AND ID_DATABANK=$dbID");
        my @rules=split(',:,',$parseRules);
        my ($idRule)=($rules[0]=~/ID=(.+)/);
        my ($orgRule,$desRule);
        if ($rules[1]) {
            if ($rules[1]=~/DES=/) {($desRule=$rules[1])=~s/DES=//;}
            else {($orgRule=$rules[1])=~s/ORG=//;}
        }
        if ($rules[2]) {($orgRule=$rules[2])=~s/ORG=//;}
    
        $self->{'idRule'}=$idRule;
        $self->{'desRule'}=$desRule if $desRule;
        $self->{'orgRule'}=$orgRule if $orgRule;
        
        ###>Fetching modifications IDs and names
        my $modificationSth=$dbh->prepare("SELECT PSI_MS_NAME,ID_MODIFICATION,MONO_MASS FROM MODIFICATION ");
        $modificationSth->execute;
        while (my ($modName,$modID,$monoMass)=$modificationSth->fetchrow_array) {
            $modificationList{$modID}="$monoMass&$modName" if ($monoMass && $modName);
        }
        
        return $self;    
    }
    sub start_document {
        my ($self,$doc)=@_;
        $self->{'queryNum'}=0;
    }
    sub end_document {
        my $self=$_[0];
    }
    
    sub start_element {
        my ($self, $element) = @_;
        $self->{'Element'} = $element->{'Name'};
        
        if ($element->{'Name'} eq 'group' && defined $element->{'Attributes'}{'{}rt'}) {            ## one group correspond to one spectra
            $self->{'queryNum'}++;
            $self->{'onGroup'}=1;
            $self->{'charge'}=$element->{'Attributes'}{'{}z'}{'Value'};
            $self->{'massexp'}=$element->{'Attributes'}{'{}mh'}{'Value'};
            $self->{'rttime'}=$element->{'Attributes'}{'{}rt'}{'Value'};
            $self->{'scan'}=$element->{'Attributes'}{'{}id'}{'Value'};
            
        }
        elsif($element->{'Name'} eq 'protein'){	
            $self->{'isOnProtein'}=1;
        }
        elsif($element->{'Name'} eq 'note'){	
            $self->{'isOnNote'}=1;
        }
        elsif($element->{'Name'} eq 'peptide'){
            $self->{'isOnPeptide'}=1;
            $self->{'protLength'}=$element->{'Attributes'}{'{}end'}{'Value'};
        }
        elsif ($element->{'Name'} eq 'domain') {    ## peptide informations
            $self->{'isOnDomain'}=1;
            $self->{'massCalc'}=$element->{'Attributes'}{'{}mh'}{'Value'};
            $self->{'delta'}=$element->{'Attributes'}{'{}delta'}{'Value'};
            $self->{'MIS'}=$element->{'Attributes'}{'{}missed_cleavages'}{'Value'};
            $self->{'sequence'}=$element->{'Attributes'}{'{}seq'}{'Value'};
            my @idDomain=split(/\./,$element->{'Attributes'}{'{}id'}{'Value'});
            $self->{'rankdomain'}="$idDomain[0]"."$idDomain[1]"."$idDomain[-1]";
            my $prev=$element->{'Attributes'}{'{}pre'}{'Value'};
            my @aaPrev=split(//,$prev);
            my $post=$element->{'Attributes'}{'{}post'}{'Value'};
            my @aaPost=split(//,$post);
            my $pos=$element->{'Attributes'}{'{}start'}{'Value'};
            my ($aaPost,$aaPrev);
            if ($aaPrev[-1]=~/\[/) {$aaPrev="-";}
            else{$aaPrev=$aaPrev[-1];}
            if($aaPost[0]=~/\]/){$aaPost="+";}
            else{$aaPost=$aaPost[0];}
            $self->{'protContext'}="$pos,$aaPrev,$aaPost";
            $self->{'start'}=$pos;
            #$self->{'score'}=$element->{'Attributes'}{'{}hyperscore'}{'Value'};
            $self->{'score'}=$element->{'Attributes'}{'{}expect'}{'Value'};
        }
        elsif ($element->{'Name'} eq 'aa') {        ## peptide modifications
            $self->{'mod'}{$element->{'Attributes'}{'{}modified'}{'Value'}}{$element->{'Attributes'}{'{}at'}{'Value'}}=$element->{'Attributes'}{'{}type'}{'Value'};
        }	
    }
    
    sub end_element {
        my ($self, $element) = @_;
        if($element->{'Name'} eq 'note'){	
            $self->{'isOnNote'}=0;
            if ($self->{'isOnProtein'}) {
                my $decoyTag="";
                my $isDecoyProtein;
                if ($self->{'protDes'}=~/reverse_/) {
                    $isDecoyProtein=1;
                    $decoyTag="DECOY_";
                }
                else{$isDecoyProtein=0;}
                my ($identifier)=($self->{'protDes'}=~/$self->{'idRule'}/);
                $identifier="$decoyTag$identifier";
                $self->{'matchDecoyProtein'}=($isDecoyProtein)? 1 : 0;	
                $self->{'matchTargetProtein'}=($isDecoyProtein)? 0 : 1;	
                $self->{'identifier'}=$identifier;
                ($refProtDes->{$self->{'identifier'}})=($self->{'protDes'}=~/$self->{'desRule'}/) if $self->{'desRule'} && !$refProtDes->{$self->{'identifier'}};
                ($refProtOrg->{$self->{'identifier'}})=($self->{'protDes'}=~/$self->{'orgRule'}/) if $self->{'orgRule'} && !$refProtOrg->{$self->{'identifier'}};
                $refProtDbRank->{$self->{'identifier'}}=1;
            }
        }
        elsif($element->{'Name'} eq 'domain'){
            $self->{'isOnDomain'}=0;
            if ($self->{'mod'}) {
                ###> Recovering peptide modifications 
                my %modif;
                foreach my $modifMass (sort {$a <=> $b} keys %{$self->{'mod'}}){
                    my $varMod;
                    my $match=0;
                    my ($numModDB,$unimod);
                    foreach my $modID (keys %modificationList){
                        my ($mass,$modName)=split(/&/,$modificationList{$modID});
                        if ($mass > $modifMass-0.005 && $mass < $modifMass+0.005) {
                            $varMod=$modName;
                            $numModDB=$modID;
                            $match=1;
                            last;
                        }
                    }
                    if (! $match) {
                        my $line;
                        open(MODFILE,"<","$dataPath/Unified_Modification_PeakView.csv") or die ("open: $!"); #find modification in Unified_Modification_PeakView file
                        while ($line=<MODFILE>) {
                            my @lineStrg=split(/;/,$line);
                            if ($lineStrg[11] eq $modifMass || $lineStrg[11] eq ($modifMass=~s/,/./) || (($lineStrg[11] > $modifMass-0.005) && ($lineStrg[11] < $modifMass+0.005)) ) {
                                $varMod=$lineStrg[14];
                                $unimod=$lineStrg[7];
                                last;
                            }
                        }
                        close MODFILE;
                        $numModDB=$dbh->fetchrow_array("SELECT ID_MODIFICATION FROM MODIFICATION WHERE UNIMOD_ACC=$unimod") if $unimod;
                    }
                    $refAnalysisMods->{'V'}{$varMod}{'numModDB'}=$numModDB;
                    $modif{$modifMass}=$varMod;
                }
                
                foreach my $modifMass (sort {$modif{$a} cmp $modif{$b}} keys %{$self->{'mod'}}){
                    my $varMod=$modif{$modifMass};
                    my ($posPeptide,$aaMod);
                    foreach my $posProt (sort {$a <=> $b} keys %{$self->{'mod'}{$modifMass}}){
                        $aaMod=$self->{'mod'}{$modifMass}{$posProt};
                        my $pos;
                        $pos=$posProt-$self->{'start'}+1;   ## peptide's modification position
                        $posPeptide=($posPeptide)? $posPeptide.".".$pos : $pos;
                        
                        if ($modifNTer && $varMod eq 'Acetyl' && $pos==1 && $aaMod ne 'K' ) {
                            $aaMod='N-term';
                        }
                    }
                     if ($aaMod eq 'N-term') {
                        $self->{'varModString'}.=" + $varMod ($aaMod)";
                    }
                    else{$self->{'varModString'}.=" + $varMod ($aaMod:$posPeptide)";}
                    $refAnalysisMods->{'V'}{$varMod}{'specificity'}=$aaMod if not defined $refAnalysisMods->{'V'}{$varMod}{'specificity'};
                    $refAnalysisMods->{'V'}{$varMod}{'specificity'}.=$aaMod if (defined $refAnalysisMods->{'V'}{$varMod}{'specificity'} && $refAnalysisMods->{'V'}{$varMod}{'specificity'}!~/$aaMod/);
                    
                }
                $self->{'actualSequence'}="$self->{'sequence'}$self->{'varModString'}";
                $self->{'mod'}=();
            }
            else{
                $self->{'actualSequence'}=$self->{'sequence'};
                $self->{'varModString'}='';
            }
            $info{"$self->{'sequence'}_$self->{'varModString'}"}{"$self->{'score'},$self->{'massCalc'},$self->{'delta'},$self->{'MIS'},$self->{'matchDecoyProtein'},$self->{'matchTargetProtein'},$self->{'varModString'}"}=1;
            $rankSeq{"$self->{'sequence'}_$self->{'varModString'}"}=$self->{'rankdomain'} unless $rankSeq{"$self->{'sequence'}_$self->{'varModString'}"};
            push @{$prot{"$self->{'sequence'}_$self->{'score'}"}{'DECOY'}},"$self->{'identifier'}&$self->{'protContext'}" if $self->{'matchDecoyProtein'};
            push @{$prot{"$self->{'sequence'}_$self->{'score'}"}{'TARGET'}},"$self->{'identifier'}&$self->{'protContext'}" if $self->{'matchTargetProtein'};
            $self->{'varModString'}='';
        }
        elsif($element->{'Name'} eq 'protein'){
            $i++;
            if ($i==2000) {
                print ".";
                $i=0;
            }
            $self->{'protContext'}='';
            $self->{'protDes'}='';
            $self->{'isOnProtein'}=0;
        }
        elsif ($element->{'Name'} eq 'peptide') {
            if (!$self->{'protSeq'} && not exists $refProtLength->{$self->{'identifier'}}) {
                my $ident;
                if ($self->{'identifier'}=~/DECOY_/){($ident=$self->{'identifier'})=~s/DECOY_//;}
                else{$ident=$self->{'identifier'}}
                $protListNoSeq{$ident}{$analysisID}=1;
            }
            else{
                my $mass=&computeAAMass($self->{'protSeq'},'protein');
                $mass=sprintf "%.2f",$mass; 
                $refProtMW->{$self->{'identifier'}}=$mass if !$refProtMW->{$self->{'identifier'}};
                $refProtLength->{$self->{'identifier'}}=$self->{'protLength'} if !$refProtLength->{$self->{'identifier'}};
                $refProtSeq->{$self->{'identifier'}}=$self->{'protSeq'} unless $refProtSeq->{$self->{'identifier'}};
            }
            $self->{'isOnPeptide'}=0;
            $self->{'matchProtSeq'}=0;
            $self->{'protSeq'}='';
        }
        elsif($element->{'Name'} eq 'group' && $self->{'onGroup'}){
           foreach my $seqMod (sort {$rankSeq{$a} <=>$rankSeq{$b}} keys %info){
                foreach my $info (keys %{$info{$seqMod}}){
                    my ($score,$massCalc,$delta,$mis,$matchDecoyProtein,$matchTargetProtein,$varModString)=split(/,/,$info);
                    my ($sequence,$mod)=split(/_/,$seqMod);
                    $varModString='' unless $varModString;
                    if ($score >= $minScore ){
                        $self->{'hitNum'}++;
                        ###>Computing max query score 
                        if ($matchTargetProtein && ((not defined $refMaxQueryScore->{$self->{'queryNum'}}) || ($refMaxQueryScore->{$self->{'queryNum'}}<$score))) {
                            $refMaxQueryScore->{$self->{'queryNum'}}=$score;
                        }
                        elsif($matchDecoyProtein && ((not defined $refMaxQueryScore->{"-$self->{'queryNum'}"}) || ($refMaxQueryScore->{"-$self->{'queryNum'}"}<$score))){
                            $refMaxQueryScore->{"-$self->{'queryNum'}"}=$score;
                        }
                        $matches{$self->{'hitNum'}}{'massCalc'}=$massCalc-1.00794;              # mh = calculated peptide mass + a proton
                        $matches{$self->{'hitNum'}}{'delta'}=$delta;
                        $matches{$self->{'hitNum'}}{'MIS'}=$mis;
                        $matches{$self->{'hitNum'}}{'sequence'}=$sequence;
                        $matches{$self->{'hitNum'}}{'vmod'}=$varModString;
                        $matches{$self->{'hitNum'}}{'score'}=$score;
                        $matches{$self->{'hitNum'}}{'matchDecoyProtein'}=$matchDecoyProtein;	
                        $matches{$self->{'hitNum'}}{'matchTargetProtein'}=$matchTargetProtein;
                        my $actualSequence="$sequence$varModString";
                        $matches{$self->{'hitNum'}}{'actualSequence'}=$actualSequence;
                        if ($matchDecoyProtein) {
                            $self->{'decoy'}++;
                            $matchSpectrum{"-$self->{'queryNum'}"}{$self->{'hitNum'}}=$score;
                            foreach my $protContext (@{$prot{"$sequence\_$score"}{'DECOY'}}){
                                my ($identifier,$beg)=split(/&/,$protContext);
                                $matches{$self->{'hitNum'}}{'protein'}=$identifier;
                                $refMatchList->{$identifier}{$actualSequence}=1;
                                $infoRank{$self->{'hitNum'}}{$identifier}{$beg}=1;
                                $matches{$self->{'hitNum'}}{'numMatchProt'}++;
                            }
                            $refCharge->{"-$self->{'queryNum'}"}=$self->{'charge'};
                            $refMassExp->{"-$self->{'queryNum'}"}=sprintf('%.4f',$self->{'massexp'}-1.00794);       # mh = parent ion mass (plus a proton) from the spectrum
                            $refMassObs->{"-$self->{'queryNum'}"}=sprintf('%.4f',($matches{$self->{'hitNum'}}{'massCalc'}+($refCharge->{"-$self->{'queryNum'}"}*1.00794))/$refCharge->{"-$self->{'queryNum'}"});
                            (my $rt=$self->{'rttime'})=~s/PT|S//g;
                            $rt=sprintf('%.2f',$rt/60);
                            $refElutionTime->{"-$self->{'queryNum'}"}="sc$self->{'scan'};et$rt";
                        }
                        elsif($matchTargetProtein){
                            $self->{'target'}++;
                            $matchSpectrum{$self->{'queryNum'}}{$self->{'hitNum'}}=$score;
                            foreach my $protContext (@{$prot{"$sequence\_$score"}{'TARGET'}}){
                                my ($identifier,$beg)=split(/&/,$protContext);
                                $matches{$self->{'hitNum'}}{'protein'}=$identifier;
                                $refMatchList->{$identifier}{$actualSequence}=1;
                                $infoRank{$self->{'hitNum'}}{$identifier}{$beg}=1;
                                $matches{$self->{'hitNum'}}{'numMatchProt'}++;
                            }
                            $refCharge->{$self->{'queryNum'}}=$self->{'charge'};
                            $refMassExp->{$self->{'queryNum'}}=sprintf('%.4f',$self->{'massexp'}-1.00794);
                            $refMassObs->{$self->{'queryNum'}}=sprintf('%.4f',($matches{$self->{'hitNum'}}{'massCalc'}+($refCharge->{$self->{'queryNum'}}*1.00794))/$refCharge->{$self->{'queryNum'}});
                            (my $rt=$self->{'rttime'})=~s/PT|S//g;
                            $rt=sprintf('%.2f',$rt/60);
                            $refElutionTime->{$self->{'queryNum'}}="sc$self->{'scan'};et$rt";
                        }
                    }
                    else{$matchSpectrum{$self->{'queryNum'}}=undef;}
                    
                }
            }
            $self->{'onGroup'}=0;
            %info=();
            %prot=();
            %rankSeq=();
        }
        elsif ($element->{'Name'} eq 'bioml') {     ## end of the xml file
            ##Compute query rank
            foreach my $queryNum (keys %matchSpectrum){
                my $rank=0;
                foreach my $hitNum (sort {$matchSpectrum{$queryNum}{$b} cmp $matchSpectrum{$queryNum}{$a} || $a <=> $b} keys %{$matchSpectrum{$queryNum}}){
                    $rank++;
                    $query2matches{$queryNum}{$rank}=$hitNum;
                }
            }
            ##Compute best peptide score
            foreach my $queryNum (keys %query2matches){
                foreach my $rank (keys %{$query2matches{$queryNum}}){
                    next if $rank>$maxRank;
                    my $hitNum=$query2matches{$queryNum}{$rank};
                    if ($matches{$hitNum}{'matchTargetProtein'}) {
                        @{$bestScoreTarget{$matches{$hitNum}{'actualSequence'}}}=($matches{$hitNum}{'score'},$hitNum) if (!$bestScoreTarget{$matches{$hitNum}{'actualSequence'}} || $bestScoreTarget{$matches{$hitNum}{'actualSequence'}}[0] < $matches{$hitNum}{'score'});
                    }
                    if ($matches{$hitNum}{'matchDecoyProtein'}) {
                        @{$bestScoreDecoy{$matches{$hitNum}{'actualSequence'}}}=($matches{$hitNum}{'score'},$hitNum) if (!$bestScoreDecoy{$matches{$hitNum}{'actualSequence'}} || $bestScoreDecoy{$matches{$hitNum}{'actualSequence'}}[0] < $matches{$hitNum}{'score'});
                    }
                }
            }
            ## Storing into %maxProtMatch, %maxProtScore, %rankProtMatch and %queryInfo
            foreach my $queryNum (keys %query2matches){
                my ($rankTarget,$rankDecoy)=(0,0);
                foreach my $rank (sort keys %{$query2matches{$queryNum}}){
                    next if $rank>$maxRank;
                    my $select;
                    my $usedRankID="";
                    my $hitNum=$query2matches{$queryNum}{$rank};
                    if ($matches{$hitNum}{'matchTargetProtein'}) {
                        $select=($bestScoreTarget{$matches{$hitNum}{'actualSequence'}}[1]==$hitNum)? 0:-1;
                        $refNumValid->{$queryNum}++ if ($select==0);
                        $rankTarget++;
                        $usedRankID="$queryNum:$rankTarget";
                        my $dataStrg="MATCH=$matches{$hitNum}{'numMatchProt'},SRK=$rank,SEL=$select,MIS=$matches{$hitNum}{'MIS'},CALC=$matches{$hitNum}{'massCalc'},SC=$matchSpectrum{$queryNum}{$hitNum},DELT=$matches{$hitNum}{'delta'},SEQ=$matches{$hitNum}{'sequence'},";
                        $dataStrg.="VMOD=$matches{$hitNum}{'vmod'}," if $matches{$hitNum}{'vmod'};
                        push @{$refQueryInfo->{$queryNum}},$dataStrg;
                    }
                    elsif($matches{$hitNum}{'matchDecoyProtein'}){
                        $select=($bestScoreDecoy{$matches{$hitNum}{'actualSequence'}}[1]==$hitNum)? 0:-1;
                        $refNumValid->{$queryNum}++ if ($select==0);
                        $rankDecoy++;
                        $usedRankID="$queryNum:$rankDecoy";
                        my $dataStrg="MATCH=$matches{$hitNum}{'numMatchProt'},SRK=$rank,SEL=$select,MIS=$matches{$hitNum}{'MIS'},CALC=$matches{$hitNum}{'massCalc'},SC=$matchSpectrum{$queryNum}{$hitNum},DELT=$matches{$hitNum}{'delta'},SEQ=$matches{$hitNum}{'sequence'},";
                        $dataStrg.="VMOD=$matches{$hitNum}{'vmod'}," if $matches{$hitNum}{'vmod'};
                        push @{$refQueryInfo->{$queryNum}},$dataStrg;
                    }
                    foreach my $identifier (keys %{$infoRank{$hitNum}}) {
                        next unless $refMatchList->{$identifier};
                        $refRankProtMatch->{$usedRankID}{$identifier}=();
                        $refMaxProtScore->{$identifier}+=$matches{$hitNum}{'score'} if $select==0;
                        foreach my $beg (keys %{$infoRank{$hitNum}{$identifier}}){
                            push @{$refRankProtMatch->{$usedRankID}{$identifier}} , $beg;
                            $refMaxProtMatch->{$identifier}{$matches{$hitNum}{'actualSequence'}}++;
                        }
                    }
                }
            }
            foreach my $identifier (keys %{$refMatchList}) {
                delete($refMatchList->{$identifier}) unless $refMaxProtScore->{$identifier};
                delete($refProtSeq->{$identifier}) unless $refMaxProtScore->{$identifier};
            }
            if (scalar keys %protListNoSeq) {
                my @analysisID;
                push @analysisID,$analysisID;
                &promsMod::getProtInfo('silent',$dbh,$dbID,\@analysisID,$refProtDes,$refProtMW,$refProtOrg,$refProtLength,undef,\%protListNoSeq);
            }
        }
        delete ($self->{'Elements'}{$element->{'Name'}});
        $self->{'Element'} ='';			
    }
    
    sub characters {
        my ($self, $element) = @_;
        if ($self->{'isOnNote'} && $self->{'isOnProtein'}) {
            $self->{'protDes'}.=$element->{'Data'};
        }
        elsif($self->{'isOnPeptide'} && !$self->{'isOnDomain'}){
            ($self->{'protSeq'}.=$element->{'Data'})=~s/\s//g;
        }
    }
    
    ###> Function that computes the mass of a peptide or a protein given its sequence
    sub computeAAMass {
        my ($sequence,$type)=@_;
        my (%massAAave,%massATave);
        if($type eq 'protein'){# Average mass for protein mass computation
            %massAAave=&promsConfig::getMassAAave; # Average mass, needed for protein mass calculation
            %massATave=&promsConfig::getMassATave; # Average mass, needed for protein mass calculation
        }elsif($type eq 'peptide'){# Isotopic mass for Peptide mass computation
            %massAAave=&promsConfig::getMassAAmono; # Average mass, needed for peptide mass calculation
            %massATave=&promsConfig::getMassATMono; # Average mass, needed for peptide mass calculation
        }
    
        my %countAA;
        my $mass=0.000000000;
        foreach my $aa (split(//,uc($sequence)) ) {
            $countAA{$aa}++;
        }
        foreach my $aa (keys %countAA) {
            $mass+=($massAAave{$aa}*$countAA{$aa}) if $massAAave{$aa}; # some characters are not amino-acids
        }
        $mass+=($massATave{H}+$massATave{H}+$massATave{O}) if $mass; # H(Nter) + OH(Cter)
        return $mass;
    }
}
1;

##################################################################
################### TANDEM.PEP.XML XML PARSING ###################
##################################################################
package XTandemPEPXMLHandler; { ## only for .tandem.pep.xml
    my %infoRank;
    my (%matches,%scorePep,$protDBSeq,%query2matches,%matchSpectrum);
    my (%modificationList,%bestScoreTarget,%bestScoreDecoy);
    my ($refProtSeq,$type,$msFilename,$dataPath,$maxRank,$dbh,$dbID,$minScore,$analysisID,$refQueryInfo,$refProtDes,$refProtOrg,$refProtLength,$refCharge,$refMassObs,$refMassExp,$refElutionTime,$refNumValid,$refMatchList,$refMaxProtMatch,$refMaxProtScore,$refMaxQueryScore,$refAnalysisMods,$refRankProtMatch,$refProtMW,$refProtDbRank);
    my $i=1;

    sub new {
        ($type,$msFilename,$dataPath,$protDBSeq,$dbh,$dbID,$minScore,$maxRank,$analysisID,$refProtSeq,$refQueryInfo,$refProtDes,$refProtOrg,$refProtLength,$refCharge,$refMassObs,$refMassExp,$refElutionTime,$refNumValid,$refMatchList,$refMaxProtMatch,$refMaxProtScore,$refMaxQueryScore,$refAnalysisMods,$refRankProtMatch,$refProtMW,$refProtDbRank) = @_;
        my $self = bless ({}, $type);
        $self->{'protDBSeq'}=$protDBSeq;
        
        ###>Fetching parsing rules
        my ($parseRules)=$dbh->selectrow_array("SELECT PARSE_RULES FROM DATABANK,DATABANK_TYPE WHERE DATABANK.ID_DBTYPE=DATABANK_TYPE.ID_DBTYPE AND ID_DATABANK=$dbID");
        my @rules=split(',:,',$parseRules);
        my ($idRule)=($rules[0]=~/ID=(.+)/);
        my ($orgRule,$desRule);
        if ($rules[1]) {
            if ($rules[1]=~/DES=/) {($desRule=$rules[1])=~s/DES=//;}
            else {($orgRule=$rules[1])=~s/ORG=//;}
        }
        if ($rules[2]) {($orgRule=$rules[2])=~s/ORG=//;}
    
        $self->{'idRule'}=$idRule;
        $self->{'desRule'}=$desRule if $desRule;
        $self->{'orgRule'}=$orgRule if $orgRule;
        
        ###>Fetching modifications IDs and names
        my $modificationSth=$dbh->prepare("SELECT PSI_MS_NAME,ID_MODIFICATION,MONO_MASS FROM MODIFICATION ");
        $modificationSth->execute;
        while (my ($modName,$modID,$monoMass)=$modificationSth->fetchrow_array) {
            $modificationList{$modID}="$monoMass&$modName" if ($monoMass && $modName);
        }
        return $self;    
    }
    sub start_document {
        my ($self, $doc) = @_;
        $self->{'queryNum'}=0;
    }
    sub end_document {
        my $self = $_[0];
    }
    
    sub start_element {
        my ($self, $element) = @_;
        $self->{'Element'} = $element->{'Name'};
        
        if ($element->{'Name'} eq 'spectrum_query') {
            my @queryInf=split(/\./,$element->{'Attributes'}{'{}spectrum'}{'Value'});
            $self->{'queryNum'}++;
            $self->{'charge'}=$element->{'Attributes'}{'{}assumed_charge'}{'Value'};
            $self->{'massexp'}=$element->{'Attributes'}{'{}precursor_neutral_mass'}{'Value'};
            $self->{'rttime'}=sprintf('%.2f',$element->{'Attributes'}{'{}retention_time_sec'}{'Value'}/60);
            $self->{'scan'}=$element->{'Attributes'}{'{}start_scan'}{'Value'};
        }
        elsif ($element->{'Name'} eq 'search_hit') {
            $self->{'hitNum'}++;
            my $entry=$element->{'Attributes'}{'{}protein'}{'Value'};
            my $decoyTag="";
            my $isDecoyProtein;
            if ($entry=~/reverse_/) {
                $isDecoyProtein=1;
                $decoyTag="DECOY_";
            }
            else{$isDecoyProtein=0;}
            my $identifier="$decoyTag$entry";
            $self->{'numMatchProt'}=$element->{'Attributes'}{'{}num_tot_proteins'}{'Value'};
            
            $self->{'masscalc'}=$element->{'Attributes'}{'{}calc_neutral_pep_mass'}{'Value'};
            $self->{'delta'}=$element->{'Attributes'}{'{}massdiff'}{'Value'};
            $self->{'MIS'}=$element->{'Attributes'}{'{}num_missed_cleavages'}{'Value'};
            $self->{'sequence'}=$element->{'Attributes'}{'{}peptide'}{'Value'};
            $self->{'matchDecoyProtein'}=($isDecoyProtein)? 1 : 0;	
            $self->{'matchTargetProtein'}=($isDecoyProtein)? 0 : 1;	
            $self->{'identifier'}=$identifier;
            
            ($refProtDes->{$identifier})=($element->{'Attributes'}{'{}protein_descr'}{'Value'}=~/(.*)\sOS=/);
            ($refProtOrg->{$identifier})=($element->{'Attributes'}{'{}protein_descr'}{'Value'}=~/$self->{'orgRule'}/) if $self->{'orgRule'} && !$refProtOrg->{$identifier};
            $refProtDbRank->{$self->{'identifier'}}=1;
            $refProtLength->{$identifier}=length($self->{'protDBSeq'}{$identifier});
            my $mass=&computeAAMass($self->{'protDBSeq'}{$identifier},'protein');
            $mass=sprintf "%.2f",$mass; 
            $refProtMW->{$identifier}=$mass;
            
            $self->{'before'}=($element->{'Attributes'}{'{}peptide_prev_aa'}{'Value'})?  $element->{'Attributes'}{'{}peptide_prev_aa'}{'Value'} : '-';
            $self->{'after'}=($element->{'Attributes'}{'{}peptide_next_aa'}{'Value'})? $element->{'Attributes'}{'{}peptide_next_aa'}{'Value'} : '-';
            
            if ($self->{'protDBSeq'}{$identifier}=~/(\w?)$self->{'sequence'}(\w?)/){
                while ($self->{'protDBSeq'}{$identifier}=~/(\w?)$self->{'sequence'}(\w?)/g) {
                    my ($pepBeg,$pepEnd)=($1,$2);
                    my $startPos=$-[0]+1;
                    my $endPos=$+[0];
                    if ($pepBeg) {$startPos++;}
                    else{$pepBeg='-';}
                    if ($pepEnd) {$endPos--;}
                    else{$pepEnd='-';}
                    push @{$self->{'protContext'}}, "$startPos,$pepBeg,$pepEnd:$identifier";
                }
            }
            else{push @{$self->{'protContext'}}, "-,$self->{'before'},$self->{'after'}:$identifier";}
        }
        elsif($element->{'Name'} eq 'alternative_protein'){
            my $entry=$element->{'Attributes'}{'{}protein'}{'Value'};
            my $decoyTag="";
            my $isDecoyProtein;
            if ($entry=~/reverse_/) {
                $isDecoyProtein=1;
                $decoyTag="DECOY_";
            }
            else{$isDecoyProtein=0;}
            my $identifier="$decoyTag$entry";
            push @{$self->{'altidentifier'}},$identifier;
            ($refProtDes->{$identifier})=($element->{'Attributes'}{'{}protein_descr'}{'Value'}=~/(.*)\sOS=/);
            ($refProtOrg->{$identifier})=($element->{'Attributes'}{'{}protein_descr'}{'Value'}=~/$self->{'orgRule'}/) if $self->{'orgRule'} && !$refProtOrg->{$identifier};
            $refProtDbRank->{$identifier}=1;
            $refProtLength->{$identifier}=length($self->{'protDBSeq'}{$identifier});
            my $mass=&computeAAMass($self->{'protDBSeq'}{$identifier},'protein');
            $mass=sprintf "%.2f",$mass; 
            $refProtMW->{$identifier}=$mass;
            
            if ($self->{'protDBSeq'}{$identifier}=~/(\w?)$self->{'sequence'}(\w?)/){
                while ($self->{'protDBSeq'}{$identifier}=~/(\w?)$self->{'sequence'}(\w?)/g) {
                    my ($pepBeg,$pepEnd)=($1,$2);
                    my $startPos=$-[0]+1;
                    my $endPos=$+[0];
                    if ($pepBeg) {$startPos++;}
                    else{$pepBeg='-';}
                    if ($pepEnd) {$endPos--;}
                    else{$pepEnd='-';}
                    push @{$self->{'protContext'}}, "$startPos,$pepBeg,$pepEnd:$identifier";
                }
            }
            else{push @{$self->{'protContext'}}, "-,$self->{'before'},$self->{'after'}:$identifier";}
        }
        elsif($element->{'Name'} eq 'search_score'){
            if ($element->{'Attributes'}{'{}name'}{'Value'} eq 'hyperscore') {
                $self->{'score'}=$element->{'Attributes'}{'{}value'}{'Value'};
                if ($self->{'score'} >= $minScore){
                    $matchSpectrum{$self->{'queryNum'}}{$self->{'hitNum'}}=$self->{'score'};
                    $scorePep{$self->{'hitNum'}}=$self->{'score'};
                    if ($self->{'matchTargetProtein'} && ((not defined($refMaxQueryScore->{$self->{'queryNum'}})) || (($refMaxQueryScore->{$self->{'queryNum'}})<($element->{'Attributes'}{'{}value'}{'Value'})))) {
                        $refMaxQueryScore->{$self->{'queryNum'}}=$element->{'Attributes'}{'{}value'}{'Value'} ;
                    }
                    elsif($self->{'matchDecoyProtein'} && ((not defined($refMaxQueryScore->{"-$self->{'queryNum'}"})) || (($refMaxQueryScore->{"-$self->{'queryNum'}"})<($element->{'Attributes'}{'{}value'}{'Value'})))){
                        $refMaxQueryScore->{"-$self->{'queryNum'}"}=$element->{'Attributes'}{'{}value'}{'Value'} ;
                    }
                }
                else{
                    $matchSpectrum{$self->{'queryNum'}}=undef;
                }
            }
        }
        elsif($element->{'Name'} eq 'mod_aminoacid_mass'){
            $self->{'mod'}{$element->{'Attributes'}{'{}mass'}{'Value'}}{$element->{'Attributes'}{'{}position'}{'Value'}}=1;
        }
    }
    
    sub end_element {
        my ($self, $element) = @_;
        if ($element->{'Name'} eq 'modification_info' && exists $self->{'mod'}) {
            ###>Recovering peptides modifications
            my %modifList;
            my @sequence=split(//,$self->{'sequence'});
            my %massAAave=&promsConfig::getMassAAmono;
            foreach my $mod ( keys %{$self->{'mod'}}){
                foreach my $pos (keys %{$self->{'mod'}{$mod}}){
                    my $aaMod;
                    if ($pos eq '[') {  # N-term modif
                        $pos='=';
                        $aaMod=$sequence[0];
                    }		
                    elsif ($pos eq ']'){# C-term modif
                        $pos='-';
                        $aaMod=$sequence[-1];
                    }		
                    else{$aaMod=$sequence[$pos-1];}
                    
                    my $varMod;
                    my $aaMass=$massAAave{$aaMod};
                    my $modifMass=$mod-$aaMass;
                    my $match=0;
                    my ($numModDB,$unimod);
                    foreach my $modID (keys %modificationList){
                        my ($mass,$modName)=split(/&/,$modificationList{$modID});
                        if ($mass > $modifMass-0.005 && $mass < $modifMass+0.005) {
                            $varMod=$modName;
                            $numModDB=$modID;
                            $match=1;
                            last;
                        }
                    }
                    if (! $match) {
                        my $line;
                        open(MODFILE,"<","$dataPath/Unified_Modification_PeakView.csv") or die ("open: $!"); #find modification in Unified_Modification_PeakView file
                        while ($line=<MODFILE>) {
                            my @lineStrg=split(/;/,$line);
                            if ($lineStrg[11] eq $modifMass || $lineStrg[11] eq ($modifMass=~s/,/./) || (($lineStrg[11] > $modifMass-0.005) && ($lineStrg[11] < $modifMass+0.005)) ) {
                                $varMod=$lineStrg[14];
                                $unimod=$lineStrg[7];
                                last;
                            }
                        }
                        close MODFILE;
                        $numModDB=$dbh->fetchrow_array("SELECT ID_MODIFICATION FROM MODIFICATION WHERE UNIMOD_ACC=$unimod") if $unimod;
                    }
                    $modifList{$varMod}{$aaMod}{$pos}=1;
                    $refAnalysisMods->{'V'}{$varMod}{'numModDB'}=$numModDB;
                }
            }
            foreach my $varMod (sort {$a cmp $b} keys %modifList){
                my $posPeptide;
                foreach my $aaMod (sort {$a cmp $b} keys %{$modifList{$varMod}}){
                    foreach my $pos (sort {$a <=> $b} keys %{$modifList{$varMod}{$aaMod}}){
                        $posPeptide=($posPeptide)? $posPeptide.".".$pos : $pos;
                    }
                    $self->{'varModString'}.=" + $varMod ($aaMod:$posPeptide)";
                    $refAnalysisMods->{'V'}{$varMod}{'specificity'}=$aaMod if not defined $refAnalysisMods->{'V'}{$varMod}{'specificity'};
                    $refAnalysisMods->{'V'}{$varMod}{'specificity'}.=$aaMod if (defined $refAnalysisMods->{'V'}{$varMod}{'specificity'} && $refAnalysisMods->{'V'}{$varMod}{'specificity'}!~/$aaMod/);
                }
            }
            $self->{'actualSequence'}="$self->{'sequence'}$self->{'varModString'}";
            $self->{'mod'}=();
        }
        elsif($element->{'Name'} eq 'search_hit'){
            if (defined $matchSpectrum{$self->{'queryNum'}}) {
                my $actualSequence=($self->{'actualSequence'})? $self->{'actualSequence'} : $self->{'sequence'};
                $matches{$self->{'hitNum'}}{'masscalc'}=$self->{'masscalc'};
                $matches{$self->{'hitNum'}}{'delta'}=$self->{'delta'};
                $matches{$self->{'hitNum'}}{'MIS'}=$self->{'MIS'};
                $matches{$self->{'hitNum'}}{'sequence'}=$self->{'sequence'};
                $matches{$self->{'hitNum'}}{'matchDecoyProtein'}=$self->{'matchDecoyProtein'};	
                $matches{$self->{'hitNum'}}{'matchTargetProtein'}=$self->{'matchTargetProtein'};
                $matches{$self->{'hitNum'}}{'protein'}=$self->{'identifier'};
                $matches{$self->{'hitNum'}}{'actualSequence'}=$actualSequence;
                $matches{$self->{'hitNum'}}{'vmod'}=$self->{'varModString'} if $self->{'varModString'};
                $refMatchList->{$self->{'identifier'}}{$actualSequence}=1;
                
                if (@{$self->{'altidentifier'}}) {
                    foreach my $identifier (@{$self->{'altidentifier'}}){
                        $refMatchList->{$identifier}{$actualSequence}=1;
                    }
                }
                foreach my $protContext (@{$self->{'protContext'}}){
                    my ($beg,$identifier)=split(/:/,$protContext);
                    $infoRank{$self->{'hitNum'}}{$identifier}{$beg}=1;
                }
                if ($self->{'matchDecoyProtein'}) {
                    $refCharge->{"-$self->{'queryNum'}"}=$self->{'charge'};
                    $refMassExp->{"-$self->{'queryNum'}"}=$self->{'massexp'};
                    $refMassObs->{"-$self->{'queryNum'}"}=sprintf('%.4f',($matches{$self->{'hitNum'}}{'masscalc'}+($refCharge->{"-$self->{'queryNum'}"}*1.00794))/$refCharge->{"-$self->{'queryNum'}"});
                    $refElutionTime->{"-$self->{'queryNum'}"}="sc$self->{'scan'};et$self->{'rttime'}";
                }
                elsif($self->{'matchTargetProtein'}){
                    $refCharge->{$self->{'queryNum'}}=$self->{'charge'};
                    $refMassExp->{$self->{'queryNum'}}=$self->{'massexp'};
                    $refMassObs->{$self->{'queryNum'}}=sprintf('%.4f',($matches{$self->{'hitNum'}}{'masscalc'}+($refCharge->{$self->{'queryNum'}}*1.00794))/$refCharge->{$self->{'queryNum'}});
                    $refElutionTime->{$self->{'queryNum'}}="sc$self->{'scan'};et$self->{'rttime'}";
                }
            }	
            $self->{'protContext'}=();
            $self->{'actualSequence'}='';
            $self->{'varModString'}='';
            $self->{'altidentifier'}=();
            $i++;
            if ($i==2000) {
                print ".";
                $i=0;
            }
        }
        elsif ($element->{'Name'} eq 'msms_pipeline_analysis') {
            ##Compute query rank
            foreach my $queryNum (keys %matchSpectrum){
                my $rank=0;
                foreach my $hitNum (sort {$matchSpectrum{$queryNum}{$b} cmp $matchSpectrum{$queryNum}{$a}} keys %{$matchSpectrum{$queryNum}}){
                    $rank++;
                    $query2matches{$queryNum}{$rank}=$hitNum;
                }
            }
            ##>Compute max peptide score
            foreach my $queryNum (keys %query2matches){
                foreach my $rank (keys %{$query2matches{$queryNum}}){
                    next if $rank>$maxRank;
                    my $hitNum=$query2matches{$queryNum}{$rank};
                    if ($matches{$hitNum}{'matchTargetProtein'}) {
                        @{$bestScoreTarget{$matches{$hitNum}{'actualSequence'}}}=($scorePep{$hitNum},$hitNum) if (!$bestScoreTarget{$matches{$hitNum}{'actualSequence'}} || $bestScoreTarget{$matches{$hitNum}{'actualSequence'}}[0] < $scorePep{$hitNum});
                    }
                    if ($matches{$hitNum}{'matchDecoyProtein'}) {
                        @{$bestScoreDecoy{$matches{$hitNum}{'actualSequence'}}}=($scorePep{$hitNum},$hitNum) if (!$bestScoreDecoy{$matches{$hitNum}{'actualSequence'}} || $bestScoreDecoy{$matches{$hitNum}{'actualSequence'}}[0] < $scorePep{$hitNum});
                    }
                }
            }
            ###>Storing into %maxProtScore, %maxProtMatch, %rankProtMatch, and %queryInfo
            foreach my $queryNum (keys %query2matches){
                my ($rankTarget,$rankDecoy)=(0,0);
                foreach my $rank (sort keys %{$query2matches{$queryNum}}){
                    next if $rank>$maxRank;
                    my $usedRankID="";
                    my $hitNum=$query2matches{$queryNum}{$rank};
                    if ($matches{$hitNum}{'matchTargetProtein'}) {
                        my $select=($bestScoreTarget{$matches{$hitNum}{'actualSequence'}}[1]==$hitNum)? 0:-1;
                        $refNumValid->{$queryNum}=1 if (!$refNumValid->{$queryNum} && $select==0);
                        $rankTarget++;
                        $usedRankID="$queryNum:$rankTarget";
                        my $dataStrg="MATCH=$self->{'numMatchProt'},SRK=$rank,SEL=$select,MIS=$matches{$hitNum}{'MIS'},CALC=$matches{$hitNum}{'masscalc'},SC=$matchSpectrum{$queryNum}{$hitNum},DELT=$matches{$hitNum}{'delta'},SEQ=$matches{$hitNum}{'sequence'},";
                        $dataStrg.="VMOD=$matches{$hitNum}{'vmod'}," if $matches{$hitNum}{'vmod'};
                        push @{$refQueryInfo->{$queryNum}},$dataStrg;
                    }
                    elsif($matches{$hitNum}{'matchDecoyProtein'}){
                        my $select=($bestScoreDecoy{$matches{$hitNum}{'actualSequence'}}[1]==$hitNum)? 0:-1;
                        $refNumValid->{"-$queryNum"}=1 if (!$refNumValid->{"-$queryNum"} && $select==0);
                        $rankDecoy++;
                        $usedRankID="-$queryNum:$rankDecoy";
                        my $dataStrg="MATCH=$self->{'numMatchProt'},SRK=$rank,SEL=$select,MIS=$matches{$hitNum}{'MIS'},CALC=$matches{$hitNum}{'masscalc'},SC=$matchSpectrum{$queryNum}{$hitNum},DELT=$matches{$hitNum}{'delta'},SEQ=$matches{$hitNum}{'sequence'},";
                        $dataStrg.="VMOD=$matches{$hitNum}{'vmod'}," if $matches{$hitNum}{'vmod'};
                        push @{$refQueryInfo->{"-$queryNum"}},$dataStrg;
                    }
                    foreach my $identifier (keys %{$infoRank{$hitNum}}) {
                        next unless $refMatchList->{$identifier};
                        $refRankProtMatch->{$usedRankID}{$identifier}=();
                        $refMaxProtScore->{$identifier}+=$scorePep{$hitNum};
                        foreach my $beg (keys %{$infoRank{$hitNum}{$identifier}}){
                            push @{$refRankProtMatch->{$usedRankID}{$identifier}} , $beg;
                            $refMaxProtMatch->{$identifier}{$matches{$hitNum}{'actualSequence'}}++;
                        }
                        $refProtSeq->{$identifier}=$self->{'protDBSeq'}{$identifier} unless $refProtSeq->{$identifier};
                    }
                }
            }
        }
        delete ($self->{'Elements'}{$element->{'Name'}});
        $self->{'Element'} ='';			
    }
    
    ###> Function that computes the mass of a peptide or a protein given its sequence
    sub computeAAMass {
        my ($sequence,$type)=@_;
        my (%massAAave,%massATave);
        if($type eq 'protein'){# Average mass for protein mass computation
            %massAAave=&promsConfig::getMassAAave; # Average mass, needed for protein mass calculation
            %massATave=&promsConfig::getMassATave; # Average mass, needed for protein mass calculation
        }elsif($type eq 'peptide'){# Isotopic mass for Peptide mass computation
            %massAAave=&promsConfig::getMassAAmono; # Average mass, needed for peptide mass calculation
            %massATave=&promsConfig::getMassATMono; # Average mass, needed for peptide mass calculation
        }

        my %countAA;
        my $mass=0.000000000;
        foreach my $aa (split(//,uc($sequence)) ) {
            $countAA{$aa}++;
        }
        foreach my $aa (keys %countAA) {
            $mass+=($massAAave{$aa}*$countAA{$aa}) if $massAAave{$aa}; # some characters are not amino-acids
        }
        $mass+=($massATave{H}+$massATave{H}+$massATave{O}) if $mass; # H(Nter) + OH(Cter)

        return $mass;
    }
}
1;

####>Revision history<####
# 1.1.7 [BUGFIX] Fix protein export issue (VS 18/11/19)
# 1.1.6 [FEATURES] Handles new monitoring system (VS 29/10/19)
# 1.1.5 Inner package : promsDIA, all steps of DIA Workflow (VS 11/06/19)
# 1.1.4 Encapsulated the promsPath within methods to avoid path retrieving issues when calling them by command line (VS 22/11/18)
# 1.1.3 Minor modification on the getProtInfo call (VS 16/11/2018)
# 1.1.2 Minor modif to handle library export errors (MLP 17/04/18)
# 1.1.1 Minor modif (MLP 10/01/18)
# 1.1.0 Add export function (MLP 06/12/17)
# 1.0.2 Minor correction (MLP 02/08/17)
# 1.0.1 Add modification to recover protein sequence (used to write analysis.fasta) (MLP 04/05/17)
# 1.0.0 Created (MLP 05/04/2017)