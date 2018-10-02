#!/usr/local/bin/perl -w

################################################################################
# storeAnalyses.cgi     3.7.2                                                  #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Stores analysis data into myProMS server and DB                              #
# Called by editProjectItem.cgi & importBatchAnalyses.cgi                      #
# $action eq 'add','anaBatch' or 'reval'                                       #
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
$| = 1;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use strict;
use XML::SAX;
use promsDIA;
#use XML::Simple;
#use File::Path qw(make_path remove_tree);
use File::Copy;
use Data::Dumper;
#######################
####>Configuration<####
#######################
my $MAX_RANK=10;
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
#my %massAAave=&promsConfig::getMassAAave; # Average mass, needed for protein mass calculation
#my %massATave=&promsConfig::getMassATave; # Average mass, needed for protein mass calculation
my @percent=(10,20,30,40,50,60,70,80,90,100); # needed to monitor process progression
my $maxRank; # Maximum number of interpretations allowed per query
my $minScore; # Minimum score allowed for interpretations
my %massValueAA = &promsConfig::getMassAAmono;
my %instrumentConver = ('ESI - QTOF'=> 'ESI-QUAD-TOF', 'ESI - ion trap'=>'ESI-TRAP', 'ESI-Linear trap'=>'ESI-TRAP', 'CID_FT-ICR'=>'ESI-FTICR', 'MALDI TOF-TOF'=>'MALDI-TOF-TOF');
my (%protSeq,%massExp,%massObs,%queryInfo,%matchedQueryMIS,%varMods,%analysisMods,%vmodsInUnimod,%matchedSeqMIS,%maxQueryScore,%numValid,%rankProtMatch,%maxProtMatch,%maxProtScore,%protDes,%protMW,%protOrg,%protLength,%protDbRank,%matchList,%matchGroup,%elutionTime,%charge,%extSpectrum,%spectralCount); # Globals (used for analysis by at least 2 subs)  @queryToSpectrum,@queryRT,@queryCharge,%protSelStatus
my (%proteinsInfo,%peptidesInfo,%spectrumPeptides); # for MSF to PDM file conversion
# Create a directory to save the MSF files
mkdir "$promsPath{valid}/multi_ana" unless -e "$promsPath{valid}/multi_ana/";

#############################
####>Fetching parameters<####
#############################
my ($batchItem,$batchItemID); # item used for batch import
my ($analysisID,$firstAnalysisID);
my $action=param('ACT');
if ($action=~/qvality\d|scanDB/) {&systemCommand; exit;}
my $decoyFile=param('decoy');
my $deleteFiles=(param('delUpFiles'))? param('delUpFiles') : 0;
my $splitMsfFiles=(param('splitMode'))? param('splitMode') : 0;
my $filePath=(param('path'))? param('path') : '';
$filePath=~s/\/\Z//; # removes ending '/' if exists
$filePath="/$filePath" if $filePath !~ /^\//; # adds '/' to path
my (@colName,@colValue);
my (%fileList,@temporaryFiles);
my $projectID=param('PROJECT_ID');

my @databankIDs; # add or reval
foreach my $dbParam ('databankID','databankID1','databankID2','databankID3') { # databankID for temp retro compatibility
	push @databankIDs,param($dbParam) if param($dbParam);
}
#my $scanDB=(param('scanDB') && param('scanDB') eq 'now')? 1 : 0; !!!OBSOLETE!!! Always performed in background (PP 09/08/12)
my $scanDB=0;
#my $selMinScore=(defined(param('minScore')) && param('minScore')=~/[^\d\.]/)? param('minScore') : 'default'; # Modified because '^' as the first character in the list negates the expressions and means - any character not in the list.
#my $selMinScore=(defined(param('minScore')) && param('minScore')=~/^\d+\.*\d*$/)? param('minScore') : 'default';
#my $maxFDR=(param('maxFDR') && param('maxFDR') > 0)? param('maxFDR') : 0; # maximum FDR allowed (%)
my $FDRalgo=param('FDRalgo') || 'qvality';
#my $des=param('des');
#my $comments=param('comments');
my ($paramScore,$paramFDR,$des,$comments)=&promsMod::cleanParameters(param('minScore'),param('maxFDR'),param('des'),param('comments'));
my $maxFDR=($paramFDR && $paramFDR > 0)? $paramFDR : 0; # maximum FDR allowed (%)
my $selMinScore=(defined($paramScore) && $paramScore=~/^\d+\.*\d*$/)? $paramScore : 'default';
my ($phenyxTaxoNumber,$phenyxTaxo,$phenyxSourceFile,$phenyxInstrument);
my $dbMixedDecoyFlag; # defined later if @databankIDs
my @databankScans; # list of analyses which proteins need to be annotated
my $isSpectra; # for X! Tandem -> used to know if the spectra option is selected


#######################
####>Starting HTML<####
#######################
my $titleAction=($action eq 'reval')? 'Reimporting' : 'Importing';
print header(-'content-encoding'=>'no',-charset=>'utf-8'); # start_html,"\n"; # start_html required to force update
warningsToBrowser(1);

print qq
|<HTML>
<HEAD>
<TITLE>Import Analyses</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT language="javascript">
var divID, streamLength, paramString;
function getWindowStream() {
	var streamText='';
	try {
		streamText=systemFrame.document.body.innerHTML;
		var newLength=streamText.length;
		if (newLength > streamLength) {
			document.getElementById(divID).innerHTML=streamText;
			streamLength=newLength;
		}
	}
	catch(e) {}
	if (systemFrame.document.readyState != 'complete') {setTimeout('getWindowStream()',500);}
	else if (streamText && (streamText.match('error') \|\| !streamText.match('__OK__'))) {
		systemFrame.location="./convertMsf2Pdm.cgi?error=1"+paramString;
	}
}
</SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<IFRAME name="systemFrame" style="display:none"></IFRAME>
<CENTER>
<FONT class="title">$titleAction Analysis Data
|;
#print " and Protein Annotations" if $scanDB; # OBSOLETE
print " in Batch Mode" if  $action eq 'anaBatch';
print "</FONT>\n</CENTER><BR><BR>\n";

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

############################
####>Checking databanks<####
############################
if (scalar @databankIDs) {
	my %decoyFlags;
	foreach my $dbID (@databankIDs) {
		my ($flag)=$dbh->selectrow_array("SELECT DECOY_TAG FROM DATABANK WHERE ID_DATABANK=$dbID");
		$decoyFlags{$flag}=1 if $flag;
	}
	$dbMixedDecoyFlag=join('|',keys %decoyFlags) if scalar keys %decoyFlags;
}
$dbh->disconnect;

##################################
####>Starting data processing<####
##################################
#if ($action eq 'add') { # OBSOLETE
#	my $dataFile=(split(/[\\\/]/,param('data_file')))[-1];
#	$fileList{$dataFile}{'anaName'} = param('name');
#	$fileList{$dataFile}{lc(param('PARENT'))} = param('PARENT_ID'); #  spot or sample
#	$fileList{$dataFile}{'file_path'} = param('data_file');
#	$fileList{$dataFile}{'pos'}=1;
#}
#els
if ($action eq 'reval') { #reval
	my $dataFile=(split(/[\\\/]/,param('data_file')))[-1];
	#$fileList{$dataFile}{'anaName'} = param('name');
	$fileList{$dataFile}{'ITEM_ID'} = param('ITEM_ID');
	$fileList{$dataFile}{'file_path'} = param('data_file');
	$fileList{$dataFile}{'pos'}=1;
}
elsif ($action eq 'anaBatch') {
	($batchItem,$batchItemID)=split(':',param('batchStart')); # item used for batch import
	my @files=param('file'); # selected files (including MSF internal requests)
	my @fileNames=param('fileName'); # all possible files (including MSF internal requests) in Batch directory
	my @anaNames=&promsMod::cleanParameters(param('anaName')); # normal import
	my @anaParents=param('selAnaParent'); # normal import
	my @targetAnalyses=param('selTargAnalysis'); # decoyFile import
	my %selectedFiles;


	print "<FONT class=\"title2\">+Pre-processing files :</FONT><BR>\n";
	foreach my $dataFile (@files) { # reformating selection list
		####--- tag INSERTION (start)---####
		#if ($dataFile =~ /:/ && $dataFile !~ /SAMPLE/) { #it is a MSF type file <-> Parsing is mandatory! (!!!SAMPLE????????)
		#	my ($dataFileName,$processingNodeAnalysis) = split(/:/,$dataFile);
		#	$selectedFiles{"$dataFileName"}=1;
		#	push @{$fileList{"$dataFileName"}{'analyses'}},$processingNodeAnalysis;
		#}
		#else {
			$selectedFiles{$dataFile}=1;
		#}
	}
	#if ($decoyFile) {
	#	my $fileRow=-1;
	#	for (my $i=0; $i<=$#fileNames; $i++) {
	#		my $fileName = $fileNames[$i];
	#		next unless defined ($selectedFiles{$fileName}); #keep only selected analyses
	#		die "Target Analysis not defined for $fileName" unless defined $targetAnalyses[$i];
	#		if (defined $fileList{$fileName}{'analyses'}) {
	#			my $j=0;# count msf file
	#			foreach my $analysis (@{$fileList{$fileName}{'analyses'}}){
	#				my $file=&parseProteomeDiscovererMSF($filePath,$fileName,$analysis,&getAllSearches("$filePath/$fileName"));
	#				push @temporaryFiles,$file; # file will be deleted from Batch dir after import
	#				my ($parItem,$parInfo)=split(':',$anaParents[$i+$j]);
	#				my @splitFile = split(/\// , $file);
	#				$fileList{$splitFile[-1]}{'anaName'} = $anaNames[$i+$j];
	#				$fileList{$splitFile[-1]}{lc($parItem)} = $parInfo;
	#				$fileList{$splitFile[-1]}{'file_path'} = $file;
	#				$fileList{$splitFile[-1]}{'msf_file'} = $fileName;
	#				$j++;
	#			}
	#			delete $fileList{$fileName}{'analyses'};
	#			delete $fileList{$fileName};
	#		}
	#		else {
	#			$fileList{$fileName}{'file_path'} = "$filePath/$fileName";
	#			$fileList{$fileName}{'target_analysis'} = $targetAnalyses[$i];
	#		}
	#	}
	#}
	#else {
	my $filePos=0;
	SEL_FILE:for (my $i=0; $i<=$#fileNames; $i++) { # looping to all files in batch dir
		next unless defined ($selectedFiles{$fileNames[$i]}); # keep only selected analyses
		if ($decoyFile) {
			die "Target Analysis not defined for $fileNames[$i]" unless defined $targetAnalyses[$i];
		}
		else {
			die "Analysis name not defined for $fileNames[$i]" unless defined $anaNames[$i]; # cheeking for parameters
			die "Parent not defined for $fileNames[$i]" unless defined $anaParents[$i]; # checking for parameters
		}
		my ($usedFileName,$usedFilePath);
		if ($fileNames[$i]=~/(.+):(\d+):(.+)\Z/) { # MSF file
			my ($msfFile,$searchNodeNumber,$fileIDstring)=($1,$2,$3);
			my $percolParamStrg=($FDRalgo eq 'precomputed' && $maxFDR)? "&percolThr=$maxFDR" : '';
			###>Launching msf to pdm conversion<###
			#my $pdmFile=&parseProteomeDiscovererMSF($dbh,$filePath,$msfFile,$searchNodeNumber,&getAllSearches("$filePath/$msfFile"));
			if ($splitMsfFiles) { # Split-mode
				my $nbFile=0;
				my @fileIDs=split(/;/,$fileIDstring);
				my $totalMergedFiles=scalar @fileIDs;
				foreach my $fileID (@fileIDs) {
					$nbFile++;
					print qq
|<FONT class=\"title3\">&nbsp;-Extracting search #$searchNodeNumber from $msfFile split-mode file #$nbFile/$totalMergedFiles:</FONT><BR>
<DIV id="streamDIV_${i}_$fileID">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>Waiting for process...</B></DIV>
<SCRIPT LANGUAGE="JavaScript">
systemFrame.location="./convertMsf2Pdm.cgi?minScore=$selMinScore$percolParamStrg&databankID=$databankIDs[0]&node=$searchNodeNumber&file=$msfFile&fileID=$fileID&path=$filePath";
divID='streamDIV_${i}_$fileID'; // reset
streamLength=0; // reset
paramString='&node=$searchNodeNumber&file=$msfFile&fileID=$fileID'; // reset
setTimeout('getWindowStream()',1500);
</SCRIPT>
|;
					#(my $pdmFile="$filePath/$msfFile")=~s/\.msf\Z/_$searchNodeNumber.pdm/;
					#$usedFileName=(split(/\//,$pdmFile))[-1];
					($usedFileName=$msfFile)=~s/\.msf\Z/_$searchNodeNumber\.$fileID.pdm/;
					$usedFilePath="$promsPath{tmp}/batch/$userID";
					my $pdmFile="$usedFilePath/$usedFileName"; # pdm generated in user's directory no matter where msf is

					#>Starting wait loop<#
					(my $databankFile="$promsPath{tmp}/$msfFile")=~s/\.msf\Z/_$searchNodeNumber\.$fileID\.$userID\.fasta/;
					(my $processEnd="$promsPath{tmp}/$msfFile")=~s/\.msf\Z/_$searchNodeNumber\.$fileID\.$userID\.end/;
					(my $processError="$promsPath{tmp}/$msfFile")=~s/\.msf\Z/_$searchNodeNumber\.$fileID\.$userID\.error/;
					my $count=my $totCount=0;
					my $pdmSize=0;
					while (!-e $processEnd && !-e $processError) {
						sleep 1;
						if (-e $pdmFile) {
							my $newSize=-s $pdmFile;
							if ($newSize > $pdmSize) {
								$count=0;
								$pdmSize=$newSize;
							}
						}
						if (++$count > 500) { # xxx cycles w/o change in pdm file
							print "<B>ERROR: The maximum duration allowed for process was exceeded. Skipping file import...</B>";
							unlink $pdmFile if -e $pdmFile;
							unlink $databankFile if -e $databankFile;
							next SEL_FILE;
						}
						print "<!--*-->\n" unless ++$totCount % 30; # keeps connection alive ~ every 30 sec
					}
					if (-e $processError) { # ERROR
						print "<B>ERROR: Premature ending of file convertion. Skipping file import...</B>";
						unlink $processError;
						unlink $pdmFile if -e $pdmFile;
						unlink $databankFile if -e $databankFile;
		#exit;
						next SEL_FILE;
					}
					unlink $processEnd if -e $processEnd;

					push @temporaryFiles,$pdmFile; # file will be deleted from Batch dir after import
					unless ($decoyFile) {
						## Copy the MSF file into the appropriate directory - Project Number, in multi_analysis directory

						# This MSF file has already been imported in this project before : need to check if it is the same MSF file or not
						if (-e "$promsPath{valid}/multi_ana/proj_$projectID/$msfFile"){
							###> Test if it is the same MSF file
							my $query="SELECT Date FROM SchemaInfo ORDER BY rowid DESC LIMIT 1";# Query to get the Date of creation of the MSF file: if it is the same name and the same date, then, it is for sure the same file!
							my $dbsqliteNew=DBI->connect("dbi:SQLite:$filePath/$msfFile", "", "", {PrintError => 1,RaiseError => 1});
							my ($dateNew)=($dbsqliteNew->selectrow_array($query));
							$dbsqliteNew->disconnect;

							my $dbsqliteOld=DBI->connect("dbi:SQLite:$promsPath{valid}/multi_ana/proj_$projectID/$msfFile", "", "", {PrintError => 1,RaiseError => 1});
							my ($dateOld)=($dbsqliteOld->selectrow_array($query));
							$dbsqliteOld->disconnect;

							my $newIndex=1;
							if ($dateNew ne $dateOld){# The already imported file is not the same... need to index the new one and check if it has not already been indexed before
								opendir (DIR, "$promsPath{valid}/multi_ana/proj_$projectID/");
								my $currentFile;
								while ( defined ($currentFile = readdir (DIR)) ) {
									next unless $currentFile =~/##\d+\.msf/;# interested in checking only MSF files that have been indexed...
									my $filewoExt=$currentFile;
									$filewoExt =~ s/\.msf//g; # File without extension
									my ($filewoExtName,$index)=split('##',$filewoExt);
									if($msfFile =~ /$filewoExtName/){ # File already indexed, check if it is the same or not...
										my $dbsqlitetemp=DBI->connect("dbi:SQLite:$promsPath{valid}/multi_ana/proj_$projectID/$currentFile", "", "", {PrintError => 1,RaiseError => 1});
										my ($dateTemp)=($dbsqlitetemp->selectrow_array($query));
										$dbsqlitetemp->disconnect;
										if ($dateTemp eq $dateNew){#it is the same file -> stop the loop
											$newIndex=$index;
											last;
										}
										else{# it is not the same -> keep the index that will rename the file after the loop ends
											$newIndex=$index+1 unless $newIndex>$index;
										}
									}
								}
								closedir DIR;

								###> Rename the MSF file with the corresponding index
								(my $msfFilewoExt=$msfFile)=~s/\.msf\Z//;
								$fileNames[$i]="$msfFilewoExt##$newIndex.msf";
								###> Copy the new indexed MSF file (only if it was not copied before).
								copy("$filePath/$msfFile","$promsPath{valid}/multi_ana/proj_$projectID/$fileNames[$i]") unless -e "$promsPath{valid}/multi_ana/proj_$projectID/$fileNames[$i]";
								###> Update the previously created PDM filename.
								my $newPDMFile="$usedFilePath/$msfFilewoExt##$newIndex"."_$searchNodeNumber.$fileID.pdm";
								move($pdmFile,$newPDMFile);
								###> Update the Fasta File that was created for fetching protein annotations and sequences
								###  (cf. &printPeptidesMSF)
								(my $rootFile=$usedFileName)=~s/\.pdm\Z//;
								my $tempOldDbFile="$promsPath{tmp}/$rootFile.$fileID.$userID.fasta";
								my $tempNewDbFile="$promsPath{tmp}/$msfFilewoExt##$newIndex"."_$searchNodeNumber.$fileID.$userID.fasta";
								move($tempOldDbFile,$tempNewDbFile);
								###> Update the variables containing the PDM filename.
								$temporaryFiles[-1]=$newPDMFile;
								$usedFileName=(split(/\//,$newPDMFile))[-1];
							}
						}
						else{
							mkdir "$promsPath{valid}/multi_ana/proj_$projectID" unless -e "$promsPath{valid}/multi_ana/proj_$projectID";
							copy("$filePath/$msfFile","$promsPath{valid}/multi_ana/proj_$projectID/$msfFile");
						}
					}
					print "<FONT class=\"title3\">&nbsp;&nbsp;Done.</FONT><BR>\n";
					$fileList{$usedFileName}{'anaName'} = "$anaNames[$i]_$fileID";
					$fileList{$usedFileName}{'file_path'} = "$usedFilePath/$usedFileName";
					$fileList{$usedFileName}{'file_id'} = $fileID;
					if ($decoyFile) {
						$fileList{$usedFileName}{'target_analysis'} = $targetAnalyses[$i];
					}
					else {
						my ($parItem,$parInfo)=split(':',$anaParents[$i]);
						$fileList{$usedFileName}{lc($parItem)} = $parInfo; # spot,sample=> id or new=sampName
					}
					$fileList{$usedFileName}{'pos'}=$filePos;
					$filePos++;
				}
				next;
			}else{ # Regulare mode: merge files stay together
				print qq
|<FONT class=\"title3\">&nbsp;-Extracting search #$searchNodeNumber from $msfFile:</FONT><BR>
<DIV id="streamDIV_$i">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>Waiting for process...</B></DIV>
<SCRIPT LANGUAGE="JavaScript">
systemFrame.location="./convertMsf2Pdm.cgi?minScore=$selMinScore$percolParamStrg&databankID=$databankIDs[0]&node=$searchNodeNumber&file=$msfFile&path=$filePath";
divID='streamDIV_$i'; // reset
streamLength=0; // reset
paramString='&node=$searchNodeNumber&file=$msfFile'; // reset
setTimeout('getWindowStream()',1500);
</SCRIPT>
|;
				#(my $pdmFile="$filePath/$msfFile")=~s/\.msf\Z/_$searchNodeNumber.pdm/;
				#$usedFileName=(split(/\//,$pdmFile))[-1];
				($usedFileName=$msfFile)=~s/\.msf\Z/_$searchNodeNumber.pdm/;
				$usedFilePath="$promsPath{tmp}/batch/$userID";
				my $pdmFile="$usedFilePath/$usedFileName"; # pdm generated in user's directory no matter where msf is

				#>Starting wait loop<#
				(my $databankFile="$promsPath{tmp}/$msfFile")=~s/\.msf\Z/_$searchNodeNumber\.$userID\.fasta/;
				(my $processEnd="$promsPath{tmp}/$msfFile")=~s/\.msf\Z/_$searchNodeNumber\.$userID\.end/;
				(my $processError="$promsPath{tmp}/$msfFile")=~s/\.msf\Z/_$searchNodeNumber\.$userID\.error/;
				my $count=my $totCount=0;
				my $pdmSize=0;
				while (!-e $processEnd && !-e $processError) {
					sleep 1;
					if (-e $pdmFile) {
						my $newSize=-s $pdmFile;
						if ($newSize > $pdmSize) {
							$count=0;
							$pdmSize=$newSize;
						}
					}
					if (++$count > 500) { # xxx cycles w/o change in pdm file
						print "<B>ERROR: The maximum duration allowed for process was exceeded. Skipping file import...</B>";
						unlink $pdmFile if -e $pdmFile;
						unlink $databankFile if -e $databankFile;
						next SEL_FILE;
					}
					print "<!--*-->\n" unless ++$totCount % 30; # keeps connection alive ~ every 30 sec
				}
				if (-e $processError) { # ERROR
					print "<B>ERROR: Premature ending of file convertion. Skipping file import...</B>";
					unlink $processError;
					unlink $pdmFile if -e $pdmFile;
					unlink $databankFile if -e $databankFile;
	#exit;
					next SEL_FILE;
				}
				unlink $processEnd if -e $processEnd;

				push @temporaryFiles,$pdmFile; # file will be deleted from Batch dir after import
				unless ($decoyFile) {
					## Copy the MSF file into the appropriate directory - Project Number, in multi_analysis directory

					# This MSF file has already been imported in this project before : need to check if it is the same MSF file or not
					if (-e "$promsPath{valid}/multi_ana/proj_$projectID/$msfFile"){
						###> Test if it is the same MSF file
						my $query="SELECT Date FROM SchemaInfo ORDER BY rowid DESC LIMIT 1";# Query to get the Date of creation of the MSF file: if it is the same name and the same date, then, it is for sure the same file!
						my $dbsqliteNew=DBI->connect("dbi:SQLite:$filePath/$msfFile", "", "", {PrintError => 1,RaiseError => 1});
						my ($dateNew)=($dbsqliteNew->selectrow_array($query));
						$dbsqliteNew->disconnect;

						my $dbsqliteOld=DBI->connect("dbi:SQLite:$promsPath{valid}/multi_ana/proj_$projectID/$msfFile", "", "", {PrintError => 1,RaiseError => 1});
						my ($dateOld)=($dbsqliteOld->selectrow_array($query));
						$dbsqliteOld->disconnect;

						my $newIndex=1;
						if ($dateNew ne $dateOld){# The already imported file is not the same... need to index the new one and check if it has not already been indexed before
							opendir (DIR, "$promsPath{valid}/multi_ana/proj_$projectID/");
							my $currentFile;
							while ( defined ($currentFile = readdir (DIR)) ) {
								next unless $currentFile =~/##\d+\.msf/;# interested in checking only MSF files that have been indexed...
								my $filewoExt=$currentFile;
								$filewoExt =~ s/\.msf//g; # File without extension
								my ($filewoExtName,$index)=split('##',$filewoExt);
								if($msfFile =~ /$filewoExtName/){ # File already indexed, check if it is the same or not...
									my $dbsqlitetemp=DBI->connect("dbi:SQLite:$promsPath{valid}/multi_ana/proj_$projectID/$currentFile", "", "", {PrintError => 1,RaiseError => 1});
									my ($dateTemp)=($dbsqlitetemp->selectrow_array($query));
									$dbsqlitetemp->disconnect;
									if ($dateTemp eq $dateNew){#it is the same file -> stop the loop
										$newIndex=$index;
										last;
									}
									else{# it is not the same -> keep the index that will rename the file after the loop ends
										$newIndex=$index+1 unless $newIndex>$index;
									}
								}
							}
							closedir DIR;

							###> Rename the MSF file with the corresponding index
							(my $msfFilewoExt=$msfFile)=~s/\.msf\Z//;
							$fileNames[$i]="$msfFilewoExt##$newIndex.msf";
							###> Copy the new indexed MSF file (only if it was not copied before).
							copy("$filePath/$msfFile","$promsPath{valid}/multi_ana/proj_$projectID/$fileNames[$i]") unless -e "$promsPath{valid}/multi_ana/proj_$projectID/$fileNames[$i]";
							###> Update the previously created PDM filename.
							my $newPDMFile="$usedFilePath/$msfFilewoExt##$newIndex"."_$searchNodeNumber.pdm";
							move($pdmFile,$newPDMFile);
							###> Update the Fasta File that was created for fetching protein annotations and sequences
							###  (cf. &printPeptidesMSF)
							(my $rootFile=$usedFileName)=~s/\.pdm\Z//;
							my $tempOldDbFile="$promsPath{tmp}/$rootFile.$userID.fasta";
							my $tempNewDbFile="$promsPath{tmp}/$msfFilewoExt##$newIndex"."_$searchNodeNumber.$userID.fasta";
							move($tempOldDbFile,$tempNewDbFile);
							###> Update the variables containing the PDM filename.
							$temporaryFiles[-1]=$newPDMFile;
							$usedFileName=(split(/\//,$newPDMFile))[-1];
						}
					}
					else{
						mkdir "$promsPath{valid}/multi_ana/proj_$projectID" unless -e "$promsPath{valid}/multi_ana/proj_$projectID";
						copy("$filePath/$msfFile","$promsPath{valid}/multi_ana/proj_$projectID/$msfFile");
					}
				}
				print "<FONT class=\"title3\">&nbsp;&nbsp;Done.</FONT><BR>\n";
			}
		}
		else{
			$usedFileName=$fileNames[$i];
			$usedFilePath=$filePath;
		}

		$fileList{$usedFileName}{'anaName'} = $anaNames[$i];
		$fileList{$usedFileName}{'file_path'} = "$usedFilePath/$usedFileName";
		(my $mzxmlFile=$usedFileName)=~s/\..*// ;
		$mzxmlFile.='.mzXML';
		$fileList{$usedFileName}{'mzxml'} = "$i";

		if ($decoyFile) {
			$fileList{$usedFileName}{'target_analysis'} = $targetAnalyses[$i];
		}
		else {
			my ($parItem,$parInfo)=split(':',$anaParents[$i]);
			$fileList{$usedFileName}{lc($parItem)} = $parInfo; # spot,sample=> id or new=sampName
		}
		$fileList{$usedFileName}{'pos'}=$filePos;
		$filePos++;
	}
}

####>Checking if remaining file(s)<####
unless (scalar keys %fileList) {
	print qq
|<BR><BR>
<FONT class="title2">&nbsp;&nbsp;File(s) not imported due to unexpected error(s).</FONT>
<BR><BR>
|;
	exit;
}

####################################
####>Looping through data files<####
####################################
$dbh=&promsConfig::dbConnect; # reconnect
my %newSamplesID;
my ($newSampleID,$newSamplePos);
my $numFiles=scalar keys %fileList;
my $filePos=0;
foreach my $dataFile (sort{$fileList{$a}{'pos'}<=>$fileList{$b}{'pos'}} keys %fileList) {
	$filePos++;

	my $computeFDR;
	if ($decoyFile) {print "<FONT class=\"title2\">+Processing decoy data for Analysis $filePos of $numFiles :</FONT><BR>\n";}
	else {print "<FONT class=\"title2\">+Processing Analysis $filePos of $numFiles :</FONT><BR>\n";}
	undef @colName;
	undef @colValue;
	push @colName,('UPDATE_DATE','UPDATE_USER');
	my $date=strftime("%Y-%m-%d %H:%M:%S", localtime);
	push @colValue,($dbh->quote($date),"'$userID'");

	undef %varMods; # moved out of &parseMascotDAT_PMF and &parseMascotDAT_MIS for common modifs detection (PP 16/07/13)
	undef %analysisMods; # new for MODIFICATION tables (PROJECT_MODIFICATION & ANALYSIS_MODIFICATION)

	if (!$decoyFile && ($action eq 'add' || $action eq 'anaBatch')) {
		###>Fetching a new ID for item
		($analysisID)=$dbh->selectrow_array("SELECT MAX(ID_ANALYSIS) FROM ANALYSIS");
		$analysisID++; #à décommenter
		push @colName,'START_DATE';
		push @colValue,$dbh->quote($date);

		##>On-the-fly creation of a new sample (startItem=experiment,gel2D,spot or parentItem=spot)
		push @colName,'ID_SAMPLE';

		#>Parent is Spot
		if ($fileList{$dataFile}{'spot'}) {
			my ($sampID)=$dbh->selectrow_array("SELECT ID_SAMPLE FROM SAMPLE WHERE ID_SPOT=$fileList{$dataFile}{spot}");
			unless ($sampID) { # create hidden sample
				my ($spotName)=$dbh->selectrow_array("SELECT NAME FROM SPOT WHERE ID_SPOT=$fileList{$dataFile}{spot}");
				unless ($newSampleID) {
					($newSampleID)=$dbh->selectrow_array("SELECT MAX(ID_SAMPLE) FROM SAMPLE");
				}
				$sampID=++$newSampleID;
				my $expID;
				if (lc($batchItem) eq 'experiment') {$expID=$batchItemID;}
				else {
					($expID)=$dbh->selectrow_array("SELECT ID_EXPERIMENT FROM SPOT,GEL2D WHERE SPOT.ID_GEL2D=GEL2D.ID_GEL2D AND ID_SPOT=$fileList{$dataFile}{spot}");
				}
				my $qSampName=$dbh->quote("Samp4:$spotName");
				my $qDate=$dbh->quote($date);
				$dbh->do("INSERT INTO SAMPLE (ID_SAMPLE,ID_SPOT,ID_EXPERIMENT,NAME,START_DATE,UPDATE_DATE,UPDATE_USER) VALUES ($sampID,$fileList{$dataFile}{spot},$expID,$qSampName,$qDate,$qDate,'$userID')");
			}
			$fileList{$dataFile}{'sample'}=$sampID;
		}
		#>Parent is free Sample (batch item is experiment)
		elsif ($fileList{$dataFile}{'sample'}=~/^new=(.+)/) {
			my $newSampName=$1;
			if ($newSamplesID{$newSampName}) { # new sample already created
				$fileList{$dataFile}{'sample'}=$newSamplesID{$newSampName};
			}
			else { # create new sample
				unless ($newSamplePos) {
					($newSamplePos)=$dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM SAMPLE WHERE ID_EXPERIMENT=$batchItemID");
				}
				$newSamplePos++;
				unless ($newSampleID) {
					($newSampleID)=$dbh->selectrow_array("SELECT MAX(ID_SAMPLE) FROM SAMPLE");
				}
				$newSampleID++;
				my $qSampName=&promsMod::resize($newSampName,48);
				$qSampName=$dbh->quote($qSampName);
				my $qDate=$dbh->quote($date);
				$dbh->do("INSERT INTO SAMPLE (ID_SAMPLE,ID_EXPERIMENT,NAME,DISPLAY_POS,START_DATE,UPDATE_DATE,UPDATE_USER) VALUES ($newSampleID,$batchItemID,$qSampName,$newSamplePos,$qDate,$qDate,'$userID')");

				$fileList{$dataFile}{'sample'}=$newSamplesID{$newSampName}=$newSampleID;
			}
		}
		push @colValue,$fileList{$dataFile}{'sample'};

		if ($dbMixedDecoyFlag) {
			push @colName,'DECOY';
			push @colValue,"'INT:BANK'";
			$computeFDR=1 if $maxFDR;
		}
	}
	elsif ($decoyFile) {
		push @colName,'DECOY';
		push @colValue,"'EXT:$dataFile'";
		$analysisID=$fileList{$dataFile}{'target_analysis'};
	}
	elsif ($action eq 'reval') { #reval
		$analysisID=$fileList{$dataFile}{'ITEM_ID'};
	}

	my $name;
	if ($action eq 'anaBatch') { # $action eq 'add' ||
		$name=$fileList{$dataFile}{'anaName'}; $name=&promsMod::resize($name,48);
		# my $des=param('des'); $des=&promsMod::resize($des,98);
		# my $comments=param('comments'); $comments=&promsMod::resize($comments,240); # < 250

		# push @colName,('DES','COMMENTS');
		# push @colValue,($dbh->quote($des),$dbh->quote($comments));
	}


	my ($fileFormat,$phenyxVersion,$wiffFile,$msType,$taxonomy,$DBFile,$instrument,$labeling,$decoy,$percolatorInPD,$percolatorInMascot,$modifNTer); #@fileData,$quantiMeth
	$phenyxTaxoNumber='';
	$phenyxTaxo='';
	$phenyxSourceFile='';
	$phenyxInstrument='';
	my ($data_fh,$uploadFile);
	my $linkedFile; # in case dataFile is symbolic link (to Mascot server)
	push @colName,'VALID_STATUS';
	my $validStatus=($action eq 'reval')? 1 : ($scanDB)? 0 : -1; # 1:Partial validation, -1:Protein data not imported, 0:Protein data imported
	push @colValue,$validStatus;

	if ($action eq 'anaBatch' || param('reval_file_ok')) { # reval_file_ok => full dat file in peptide_data
#print $fileList{$dataFile}{'file_path'}; #debug
		#open (DATAFILE, $fileList{$dataFile}{'file_path'});
		#@fileData=(<DATAFILE>);
		#close DATAFILE;
		open ($data_fh, $fileList{$dataFile}{'file_path'});

		if (-l $fileList{$dataFile}{'file_path'}) { # file is symbolic link
			$linkedFile=readlink($fileList{$dataFile}{'file_path'});
			unlink $fileList{$dataFile}{'file_path'} if param('reval_file_ok');
		}
	}
	else { # reval with old pepXXXXXX.dat file !!! NO LONGER SUPPORTED !!!
		$uploadFile=upload('data_file');	# also works : =param('data_file');
		####<Loading file in memory
		#@fileData=(<$uploadFile>);
		$data_fh=$uploadFile;
	}
	my @headData; # former global @fileData
	while(<$data_fh>) {
		push @headData,$_;
		last if $#headData==2;
	}
	close $data_fh;


	####>Checking if file is from MASCOT or Phenyx or Paragon (.group file converted to .xml via ProteinPilot sofware)
	if ($headData[0] =~ /MIME-Version: .+ Mascot/) {$fileFormat='MASCOT.DAT';} # default. Will be changed to [MASCOT/SEQUEST].PDM if from Protein Discover msf
	elsif ($headData[1] =~ /<mascot_search_results/) {$fileFormat='MASCOT.XML';}
# 	elsif ($headData[1] =~ /phenyx-ms/) {$fileFormat='PHENYX.XML';}
	elsif ($headData[2] =~ /IdentificationList/ && $headData[1] =~ /IdentificationResult.+version="([\d\.]+)"/) {$fileFormat='PHENYX.XML'; $phenyxVersion=$1;}
	elsif ($headData[0] =~ /<RESULTS/) {$fileFormat='PARAGON.XML';}
	elsif ($headData[0] =~ /<?xml/ && ($headData[2] =~ /<bioml/ || $headData[1] =~ /<bioml/)){$fileFormat='TDM.XML';}
	elsif ($headData[0] =~ /<?xml/ && $headData[2] =~ /<msms_pipeline_analysis /){$fileFormat='TDM.PEP.XML';}
	else {
		#&badFileName($uploadFile);
		print "<FONT class=\"title3\">ERROR: Unknown file format! Skipping file $dataFile.</FONT><BR>\n";
		delete $fileList{$dataFile};
		next; # next analysis
	}
	$firstAnalysisID=$analysisID unless $firstAnalysisID;

# 	if ($dataFile !~ /^F[0-9]{6}\.dat/) {&badFileName($uploadFile);}

	###>Copying file to validationPath/ana_$analysisID
	open (DATAFILE, $fileList{$dataFile}{'file_path'});
	mkdir "$promsPath{valid}/ana_$analysisID" unless -e "$promsPath{valid}/ana_$analysisID"; # <<<<<<<<<< ANALYSIS DIR

	####> upload mzXML files in anaDir (only for XTandem or dat files)
	if ($FDRalgo eq 'mayu') {
		(my $searchName=$dataFile)=~s/\.\D+//;
		if (param("mzxml$fileList{$dataFile}{'mzxml'}")) {
			my $mzxmlFile=upload("mzxml$fileList{$dataFile}{'mzxml'}");
			my $newFile="$promsPath{data}/tmp/upload/project_$projectID/$mzxmlFile";
			my $tmpfile = tmpFileName($mzxmlFile);
			move($tmpfile,$newFile);
		}
	}

	#next;
	if ($linkedFile) {symlink($linkedFile,"$promsPath{valid}/ana_$analysisID/$dataFile");} # anabatch || reval with symlink
	else {open (FILE,">$promsPath{valid}/ana_$analysisID/$dataFile");}
	if ($fileFormat eq 'MASCOT.DAT') { # default, can be changed to XXX.PDM below
		$percolatorInMascot=1 if ($FDRalgo eq 'precomputed' && $maxFDR && -e "$promsPath{tmp}/batch/$userID/.percolator/$dataFile.target.pop" );

		my %tempVarMods; #used for convertion of custom label mods to real unimod mods
		my $section='';
		#foreach my $line (@fileData) { #}
		while (my $line=<DATAFILE>) {
			print FILE $line unless $linkedFile;
			if ($line=~/^Content-Type: .+ name="(\w+)"/) {
				$section=$1;
				###if ($section eq 'peptides' && $maxFDR && $decoy eq 'INT:SEARCH') { # target
				###	$computeFDR=1;
				###	(my $rootName=$dataFile)=~s/\.[^\.]+\Z//;
				###	open (TARGET,">$promsPath{valid}/ana_$analysisID/$rootName.target");
				###	open (DECOY,">$promsPath{valid}/ana_$analysisID/$rootName.decoy");
				###}
				next;
			}
			if ($section eq 'parameters') {
				if ($line=~/^FILE=/) {
					if ($line=~/([^ =\/\\]+)\Z/) {
						$wiffFile=$1;
					}
					elsif ($line=~/^FILE=(.+)/) {
						$wiffFile=$1;
						$wiffFile=~s/\s+\Z//;
					}
					chomp($wiffFile) if $wiffFile; # otherwise \n at end of string
				}
				#elsif ($line=~/^IT_MODS=(.*Label:.+)/) { # SILAC: 'Label:xxx' mod not always used for search
				#	foreach my $mod (split(/,/,$1)) {
				#		$mod=~s/^_+//;
				#		$mod=~s/\s+\Z//;
				#		if ($mod=~/^(Label:.+)/) {
				#			$labeling.=',' if $labeling;
				#			$labeling.=$1;
				#		}
				#	}
				#}
				elsif ($line=~/^TAXONOMY=[\. ]*(.+)/) {$taxonomy=$1; $taxonomy=~s/\s*\Z//;} # remove trailing space,\n, etc..
				elsif ($line=~/^SEARCH=(\S+)/) {$msType=$1;}
				elsif ($line=~/^INSTRUMENT=(.+)/) {$instrument = $1; $instrument =~ s/\s+\Z//;}
				#elsif ($line=~/^QUANTITATION=(.+)/) {$quantiMeth = $1; $quantiMeth =~ s/\s+\Z//;}
				elsif ($line=~/^QUANTITATION=(.+)/ && $1 ne 'None' && $1 !~ /LABEL[-: ]FREE/i) {$labeling = $1; $labeling =~ s/\s+\Z//;} # overwrites data from IT_MODS / QUANTITATION=None in MALDI experiments...
				elsif ($line=~/^DECOY=1/) {
					unless ($dbMixedDecoyFlag) { # if mixed databank => OK for Qvality
						$decoy='INT:SEARCH';
						$computeFDR=1 if ($maxFDR && !$decoyFile); # && $FDRalgo ne 'precomputed'
					}
				}
				# Update the $minScore valor (Xcorr for Sequest goes from 0 to 5)
				elsif ($line=~ /^SEARCH_ALGO/) {
					if ($line=~/SEQUEST/) {
						$fileFormat='SEQUEST.PDM';
						#$minScore=&promsConfig::getMinScore('SEQUEST');
					}
					else {$fileFormat='MASCOT.PDM';}
				}
				# Percolator in Proteome Discoverer
				elsif ($line=~ /^FDR_ALGO=Percolator/) {$percolatorInPD=1;}

				#>Sequence modifications (moved here! (PP 16/07/13))
				# Possible fix modification (parameters section)
				elsif ($line=~/^MODS=(.+)\n/) {
					(my $modString=$1)=~s/\s*\Z//; # trailing \r !!!?
					foreach my $mod (split(/,/,$modString)) {
						my ($fixModName,$specificity)=&promsMod::convertVarModString($mod);
						$analysisMods{'F'}{$mod}{'specificity'}=$specificity;
						$analysisMods{'F'}{$mod}{'numModDB'}=&promsMod::getModificationIDfromString($dbh,$fixModName,$specificity,\%vmodsInUnimod);
					}
				}
				# Possible variable modification (parameters section)
				elsif ($line=~/^IT_MODS=(.+)\n/) {
					(my $modString=$1)=~s/\s*\Z//; # trailing \r !!!?
					my $numMod=1;
					foreach my $mod (split(/,/,$modString)) {
						$mod=~s/^_+//; # __Oxidation -> Oxidation !!! not a normal Oxydation (PP 07/05/12)!!
						$varMods{$numMod}=$mod;
						#$analysisMods{'V'}{$mod}{'numMod'}=$numMod; # Acetyl (K)
						my ($varModName,$specificity)=&promsMod::convertVarModString($mod);
						#$analysisMods{'V'}{$mod}{'specificity'}=$specificity;
						#$analysisMods{'V'}{$mod}{'numModDB'}=&promsMod::getModificationIDfromString($dbh,$varModName,$specificity,\%vmodsInUnimod);

						@{$tempVarMods{$mod}}=($numMod,$varModName,$specificity);
						$numMod++;
					}
				}
			}
			elsif ($section eq 'masses') {
				if ($line=~/^delta\d*=([^,]+),(.+)/) {
					my ($delta,$mod)=($1,$2);
					push @{$tempVarMods{$mod}},$delta if $tempVarMods{$mod}; # not for fix mods
				}
			}
			elsif ($section eq 'quantitation') {
				if ($fileFormat=~/\.PDM/ && $line=~/AminoAcids="([^\"]+)" Name="([^\"]+)" .+ DeltaMass="([^\"]+)/) {
					my ($newSpecif,$varModName,$deltaMass)=($1,$2,$3); # Same modif can be listed multiple times with diff specificity (eg Label:13C(6) -> K then Label:13C(6) -> R) => combine!!!
					foreach my $mod (keys %tempVarMods) {
						if (abs($tempVarMods{$mod}[3]-$deltaMass) < 0.2) { # match! => replace with new values
							my $numMod=$tempVarMods{$mod}[0];
							my $specificity=$tempVarMods{$mod}[2] || '';
							foreach my $res (split(//,$newSpecif)) { # just to safe. there should be only 1 res here
								next unless $res=~/[A-Z]/;
								next if $specificity=~/$res/;
								$specificity.=',' if $specificity;
								$specificity.=$res;
							}
							@{$tempVarMods{$mod}}=($numMod,$varModName,$specificity,$deltaMass);
							last;
						}
					}
				}
			}
			elsif ($section eq 'header') { # next section after quantitation => convert varmods into IDs
				foreach my $mod (keys %tempVarMods) {
					$analysisMods{'V'}{$mod}{'numMod'}=$tempVarMods{$mod}[0];
					$analysisMods{'V'}{$mod}{'specificity'}=$tempVarMods{$mod}[2];
					$analysisMods{'V'}{$mod}{'numModDB'}=&promsMod::getModificationIDfromString($dbh,$tempVarMods{$mod}[1],$tempVarMods{$mod}[2],\%vmodsInUnimod);
				}
			}

			###elsif ($computeFDR && $section=~/peptides/) { # target or decoy
			###	if ($line=~/^q\d+_p1=\d/) { # Only top match per query (rank1)
			###		my $sc=(split(/,/,$line))[7];
			###		if ($section eq 'peptides') {print TARGET "$sc\n";}
			###		else {print DECOY "$sc\n";}
			###	}
			###}
		}

		###if ($computeFDR) {
		###	close TARGET;
		###	close DECOY;
		###}
		$wiffFile='Not found' unless $wiffFile;
	}
	elsif ($fileFormat eq 'MASCOT.XML') {
		#foreach my $line (@fileData) { #}
		while (my $line=<DATAFILE>) {
			print FILE $line;
			next if $decoyFile;
			# if ($line=~/^<COM>(.+)<\/COM>/) {
				# my $tempFile=(split(/[\/\\]/,$1))[-1];
				# ($wiffFile)=($tempFile=~/^(\S+)/);
			# }
			if ($line=~/^<FILENAME>(.+)<\/FILENAME>/) {
				#my $tempFile=(split(/[\/\\]/,$1))[-1];
				($wiffFile)=($1=~/.*[\/|\\](\w+\.\w+)/);
			}
# 			if ($line=~/^<FILENAME>(.+)<\/FILENAME>/) {$wiffFile=(split(/[\/\\]/,$1))[-1]; $wiffFile=&promsMod::resize($wiffFile,20);}
			if ($line=~/^<TAXONOMY>\W*(.+)<\/TAXONOMY>/) {$taxonomy=$1;} # remove trailing space,\n, etc..
				#elsif ($line=~/<URI>.+\/data\/(.+)<\/URI>/) {$remoteFile=$1;}
			elsif ($line=~/<SEARCH>(\S+)<\/SEARCH>/) {$msType=$1;}
			elsif ($line=~/<INSTRUMENT>(\S+)<\/INSTRUMENT>/) {$instrument=$1;}
		}
		$wiffFile='Not found' unless $wiffFile;
	}
	elsif ($fileFormat eq 'PHENYX.XML') {
		$wiffFile='';
		$taxonomy='';
		$instrument = '';
		#$remoteFile='';
		$msType='MIS';
# 		foreach my $line (@fileData) {
# 			print FILE $line;
# 		}
		#foreach my $i (0..$#fileData) { #}
		while (my $line=<DATAFILE>) {
			$line=~s/<id[rl]:/</g; $line=~s/<\/id[rl]:/<\//g;	# Converts to v2.4+
			print FILE $line;
		}
	}
	elsif ($fileFormat eq 'PARAGON.XML') {
		my ($onSettings,$onSearch)=(0,0);
		$wiffFile='';
		$taxonomy='';
		$msType='MIS';
		#foreach my $line (@fileData) { #}
		while (my $line=<DATAFILE>) {
			print FILE $line;
			#print "$line<BR>";
			if ($line =~ /<PARAGON_METHOD_SETTINGS>/){
				$onSettings=1;
			}
			elsif ($line =~ /<\/PARAGON_METHOD_SETTINGS>/) {
				$onSettings=0;
			}

			if ($line =~ /<SEARCH/) {
				$onSearch=1;
			}
			elsif($line =~ /<\/SEARCH>/) {
				$onSearch=0;
			}

			if ($onSearch && $line =~ /<PARAMETERS .* do_reverse_search=\"(\w+)\"/) {
				my $decoyState=$1;
				if ($decoyState eq 'true'){
					$computeFDR=1 if $maxFDR;
					$decoy='INT:SEARCH';
				}
				$onSearch=0;
			}

			if ($onSettings) {
				if($line =~ /<UI_SPECIES>(.*)<\/UI_SPECIES>/){
					$taxonomy=$1;
				}
				elsif ($line =~ /<UI_INSTRUMENT>(.*)<\/UI_INSTRUMENT>/) {
					my $aliasInstrument=$1;
					($instrument)=$dbh->selectrow_array("SELECT NAME FROM INSTRUMENT WHERE ALIAS='$aliasInstrument'");
				}
				elsif ($line =~ /<UI_CYS_ALKYLATION>(.*)<\/UI_CYS_ALKYLATION>/) {
					my $mod=$1;
					$analysisMods{'F'}{$mod}{'specificity'}='C';
					$analysisMods{'F'}{$mod}{'numModDB'}=&promsMod::getModificationIDfromString($dbh,$mod,$analysisMods{'F'}{$mod}{'specificity'},\%vmodsInUnimod);
				}
				elsif ($line =~ /<UI_QUANT_TYPE>(.+)<\/UI_QUANT_TYPE>/) { # only if a quantification was set
					$labeling=$1;
				}
			}

			if (!$wiffFile && $line =~ /<PEAKLIST .* originalfilename=\"(.*?)\"/) {
				my $fileInfo=$1;
				$wiffFile=(split(/\\/,$fileInfo))[-1];# ex: C:\Documents and Settings\Administrator\Desktop\QStar\E5682BL.wiff
			}
		}# If not, the file created is empty...

		#if ($computeFDR) {
		#	push @colName,'DECOY';
		#	push @colValue,"'INT:SEARCH'";
		#}

	}
	elsif ($fileFormat eq 'TDM.PEP.XML' || $fileFormat eq 'TDM.XML'){
		while (my $line=<DATAFILE>) {
			print FILE $line;
			if ($line=~/spectrum, path/) {
				($wiffFile)=($line=~/>.*\/(.*)<\/note>/) if $fileFormat eq 'TDM.XML';
				($wiffFile)=($line=~/value=".*\/(.*)"\/>/) if $fileFormat eq 'TDM.PEP.XML';
			}
			elsif($line=~/output, spectra/){
				($isSpectra)=($line=~/>(\w+)<\/note>/) if $fileFormat eq 'TDM.XML';
				($isSpectra)=($line=~/value="(\w+)"\/>/) if $fileFormat eq 'TDM.PEP.XML';
			}
			elsif($line=~/potential N-terminus modifications">(.*)<\/note>/){
				$modifNTer=($1) ? 1 : 0;
			}
		}
		$wiffFile='' unless $wiffFile;
		($taxonomy,$DBFile)=$dbh->selectrow_array("SELECT ORGANISM,FASTA_FILE FROM DATABANK WHERE ID_DATABANK=$databankIDs[0]");
		$taxonomy='' unless $taxonomy;
		$instrument='ESI-QUAD-TOF';
		$msType='MIS';
	}
	close FILE unless $linkedFile;
	close DATAFILE;

	##>Min score & max rank
	if ($decoyFile || $action eq 'reval') {
		($minScore,$maxRank)=$dbh->selectrow_array("SELECT MIN_SCORE,MAX_RANK FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID"); # will be overwritten if decoyFile
	}
	else {
		$maxRank=($msType eq 'PMF')? $MAX_RANK : param('maxRank'); # No limit if PMF

		my $okPrecomputedFDR=0;
		if ($FDRalgo eq 'precomputed' && $maxFDR) {
			if ($percolatorInPD || $percolatorInMascot) {
				$minScore=0;
				$computeFDR=0;
				$okPrecomputedFDR=1;
			}
			else {
				$FDRalgo='qvality';
				print "<FONT class=\"title3\">WARNING: Precomputed FDR data not found! Computing FDR with Qvality...</FONT><BR>\n";
			}
		}
		if ($computeFDR) {$minScore=0;} # => Extract all data! -> will be updated after Qvality/DTcount execution
		elsif (!$okPrecomputedFDR) { # not if PD Percolator
			if ($selMinScore=~/def/i) {$minScore=&promsConfig::getMinScore($fileFormat);} # default
			else {$minScore=$selMinScore;}
		}###--> $minScore will be updated later if successful FDR computation
	}

	$maxRank=&promsConfig::getMaxRank unless $maxRank; # just to be safe
	$minScore=0 unless $minScore; # just to be safe
	# print "File Loaded (MS type=$msType)!<BR>\n";
# exit;

	if (!$decoyFile) {
		#push @colName,'ID_DATABANK'; OBSOLETE!!!!!!!!! # if reval + new DB
		#push @colValue,$databankID;

		if ($action eq 'anaBatch') { # $action eq 'add' ||
			push @colName,('DES','COMMENTS');
			push @colValue,($dbh->quote($des),$dbh->quote($comments));
			push @colName,'DATA_FILE';
			push @colValue,$dbh->quote($dataFile);
			push @colName,'FILE_FORMAT';
			push @colValue,$dbh->quote($fileFormat);
			push @colName,'MS_TYPE';
			push @colValue,$dbh->quote($msType);
			push @colName,'TAXONOMY';
			#$taxonomy=&promsMod::resize($taxonomy,50);
			push @colValue,$dbh->quote($taxonomy);
			push @colName,'WIFF_FILE';
			push @colValue,$dbh->quote($wiffFile);
			push @colName,'INSTRUMENT';
			#$instrument=&promsMod::resize($instrument,28);
			push @colValue,$dbh->quote($instrument);
			if ($fileFormat eq 'PHENYX.XML') {
				push @colName,'MAX_RANK';
				push @colValue,10;
			}
			elsif ($msType ne 'PMF') {
				push @colName,'MAX_RANK';
				push @colValue,$maxRank;
			}
			if ($msType ne 'PMF') {
				push @colName,'MIN_SCORE';
				push @colValue,$minScore;
			}
			#if ($quantiID) {
			#	push @colName,'ID_QUANTIF_METHOD';
			#	push @colValue,$quantiID;
			#}
			if ($labeling) {
				push @colName,'LABELING';
				push @colValue,$dbh->quote($labeling);
			}
			if ($decoy) {
				push @colName,'DECOY';
				my $decoyStrg=$decoy;
				###$decoyStrg.=",FDR:$maxFDR" if $computeFDR;
				push @colValue,$dbh->quote($decoyStrg);
			}
		}
	}

	################################
	####>Generating SQL Queries<####
	################################
	if (!$decoyFile && ($action eq 'add' || $action eq 'anaBatch')) {
		my $colNameString=join(",",@colName);
		my $colValueString=join(",",@colValue);

		###>Generating a new DISPLAY_POS default value<###
		my ($displayPos) = $dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM ANALYSIS WHERE ID_SAMPLE=$fileList{$dataFile}{sample}");
		$displayPos++;

			# no valid labels => only 1 entry
		$name=$dbh->quote($name);

		$dbh->do("INSERT INTO ANALYSIS (ID_ANALYSIS,NAME,DISPLAY_POS,LOWER_SCORES,$colNameString) VALUES ($analysisID,$name,$displayPos,0,$colValueString)") || die $dbh->errstr; #à décommenter

		my $dbRank=1;
		foreach my $dbID (@databankIDs) {
			$dbh->do("INSERT INTO ANALYSIS_DATABANK (ID_ANALYSIS,ID_DATABANK,DB_RANK) VALUES ($analysisID,$dbID,$dbRank)"); #à décommenter
			$dbRank++;
		}
		$dbh->commit; # id_analysis is protected as soon as possible!

	}
	elsif ($decoyFile || $action eq 'reval') {
		my $updateString;
		for (my $c=0;$c<=$#colName;$c++) { # looping through all conditions
			$updateString.="$colName[$c]=$colValue[$c]";
			$updateString.=',' if $c<$#colName;
		}
		$dbh->do("UPDATE ANALYSIS SET $updateString WHERE ID_ANALYSIS=$analysisID") || die $dbh->errstr;
	}
# $dbh->disconnect; exit; # DEBUG!!!

	######################################################
	####>Importing MS analysis data +/- databank data<####
	######################################################
	undef %massExp;
	undef %massObs;
	undef %queryInfo;
	undef %maxQueryScore;
	undef %matchedQueryMIS;
	undef %matchedSeqMIS;
	undef %numValid;
	undef %rankProtMatch;
	undef %maxProtMatch;
	undef %maxProtScore;
	undef %protDes;
	undef %protMW;
	undef %protOrg;
	#undef %protSelStatus;
	undef %protLength;
	undef %protDbRank;
	undef %matchList; # needed to compute match groups
	undef %matchGroup;
	undef %elutionTime;
	undef %extSpectrum;
	undef %protSeq;
	###undef %qvalityData;
	#undef @queryToSpectrum;# Globals (used for analysis by at least 2 subs)
	#undef @queryRT;
	#undef @queryCharge;
	#my ($queryID,$protValID); # Globals (used for analysis by at least 2 subs)

	####>Extracting data from Search Results file<####
#if ($msType ne 'PMF' && $msType ne 'MIS') { # (22/07/10)
#	print "<B>WARNING: Unknown search ($msType). Assuming search is MS/MS...</B><BR>\n";
#	$msType='MIS';
#}
	if ($fileFormat eq 'MASCOT.DAT' && $msType eq 'PMF') {&parseMascotDAT_PMF($dataFile);} # ,\@fileData
	elsif (($fileFormat=~/\.(DAT|PDM)/) && $msType eq 'MIS') {&parseMascotDAT_MIS($dataFile,0,$percolatorInMascot);} # ,\@fileData # Mascot or Protein Discover MSF
	elsif ($fileFormat=~/\.DAT/ && $msType eq 'SQ') { # Mascot mixed PMF/MIS search
		&parseMascotDAT_MIS($dataFile,1,$percolatorInMascot); # ,\@fileData
		&parseMascotDAT_PMF($dataFile,1); # ,\@fileData
	}
	elsif ($fileFormat eq 'MASCOT.XML' && $msType eq 'MIS') {&parseMascotXML_MIS($dataFile);}
	elsif ($fileFormat eq 'PHENYX.XML') {&parsePhenyxXML($dataFile,$phenyxVersion);}
	elsif ($fileFormat eq 'PARAGON.XML') {
		&parseParagonXML($dataFile,\%analysisMods);
		foreach my $anaVmod (keys %{$analysisMods{'V'}} ) {
			$analysisMods{'V'}{$anaVmod}{'numModDB'}=&promsMod::getModificationIDfromString($dbh,$anaVmod,$analysisMods{'V'}{$anaVmod}{'specificity'},\%vmodsInUnimod);
		}
	}
	elsif ($fileFormat eq 'TDM.PEP.XML') {
		my ($parseRules)=$dbh->selectrow_array("SELECT PARSE_RULES FROM DATABANK,DATABANK_TYPE WHERE DATABANK.ID_DBTYPE=DATABANK_TYPE.ID_DBTYPE AND ID_DATABANK=$databankIDs[0]");
		my @rules=split(',:,',$parseRules);
		my ($idRule)=($rules[0]=~/ID=(.+)/);
		my ($orgRule,$desRule);
		if ($rules[1]) {
			if ($rules[1]=~/DES=/) {($desRule=$rules[1])=~s/DES=//;}
			else {($orgRule=$rules[1])=~s/ORG=//;}
		}
		if ($rules[2]) {($orgRule=$rules[2])=~s/ORG=//;}

		##recovering of protein's infos
		open(DBFILE, "<$promsPath{data}/banks/db_$databankIDs[0]/$DBFile") or die ("open: $!");
		my ($identifier,$org,%protDBSeq);
		while (my $line=<DBFILE>) {
			if ($line=~/^>/) {
				(my $entry=$line)=~s/^>\s*//;
				my $decoyTag="";
				my $reverse="";
				if ($entry=~/reverse_/) {
					$decoyTag="DECOY_";
				}
				($identifier)=($entry=~/iRT/)? ($entry=~/^(\S+)/) : ($entry=~/$idRule/);
				$reverse="reverse_" if $identifier!~/reverse/ && $decoyTag eq "DECOY_";		## -> suppression of the "reverse_" tag with idRule in uniprot ACC
				$identifier=$decoyTag.$reverse.$identifier ;
			}
			else{
				chomp $line;
				$protDBSeq{$identifier}.=$line;
			}
		}
		close DBFILE;
		&parseTandemPEPXML($dataFile,\%protDBSeq);
	}
	elsif($fileFormat eq 'TDM.XML'){&parseTandemXML($dataFile,$modifNTer);}


# !!!!!!!!!!! PP (26/04/13): label Mod pas forcément: -unique (heavy,medium), -sur tous les AAs (SILAC,iTRAQ) !!!!!!!!!!!!!!!
	#if ($labeling && $labeling ne 'None' && $labeling !~ /LABEL[-: ]FREE/i) { # Add label modification in ANALYSIS_MODIFICATION
	#	$analysisMods{'L'}{$labeling}{'specificity'}='.';# All the residues
	#	$analysisMods{'L'}{$labeling}{'numModDB'}=&promsMod::getModificationIDfromString($dbh,$labeling,'.',\%vmodsInUnimod);;# All the residues
	#	$dbh->do("UPDATE MODIFICATION SET IS_LABEL=1 WHERE ID_MODIFICATION=$analysisMods{L}{$labeling}{numModDB}");
	#	$dbh->commit;
	#} Commented by PP 11/07/13

	####>Update Modification tables for this analysis and the project<####
	my %anaModDbData; # key = composite primary key (anaID:modID:modifType), value = specificity

	foreach my $modType (keys %analysisMods) {
		foreach my $mod (keys %{$analysisMods{$modType}}) {
			if ($mod =~ /\+/) {# Do not store into %anaModDbData double modification type : i.e. Label:2H(4)+GG (SILAC labeling and GlyGly on same residue !)
				print "<FONT class=\"title3\">WARNING: $mod is treated as a double modification on the same amino acid...</FONT><BR>\n";
				next;
			}
			my $modID = $analysisMods{$modType}{$mod}{'numModDB'};
			if (defined $anaModDbData{"$analysisID:$modID:$modType"}) {
				# concat existing specificity with new one
				my $specificity=quotemeta($analysisMods{$modType}{$mod}{'specificity'});
				next if ($anaModDbData{"$analysisID:$modID:$modType"} && $anaModDbData{"$analysisID:$modID:$modType"} =~ /$specificity/);
				$anaModDbData{"$analysisID:$modID:$modType"} .= ','.$analysisMods{$modType}{$mod}{'specificity'};
			}
			else {
				# new entry
				$anaModDbData{"$analysisID:$modID:$modType"} = $analysisMods{$modType}{$mod}{'specificity'};
			}
		}
	}

	my $sthInsAM=$dbh->prepare("INSERT INTO ANALYSIS_MODIFICATION (ID_ANALYSIS,ID_MODIFICATION,MODIF_TYPE,SPECIFICITY) VALUES (?,?,?,?)");
	while (my ($primaryKey, $specificity) = each %anaModDbData) {
		my ($anaID, $modID, $modType) = split(/:/,$primaryKey);
		$sthInsAM->execute($anaID,$modID,$modType,$specificity) or die $dbh->errstr;
	}
	$sthInsAM->finish;



	####>FDR Algorithm<####
	if ($computeFDR) {
		($computeFDR,$minScore)=&applyFDR($FDRalgo,$maxFDR,$analysisID,$dataFile,&promsConfig::getMinScore($fileFormat)); # $computeFDR & $minScore updated by function
		if ($computeFDR) { # ok or 0 decoy match
			$dbh->do("UPDATE ANALYSIS SET MIN_SCORE=$minScore,DECOY=CONCAT(DECOY,',FDR=$maxFDR:$FDRalgo') WHERE ID_ANALYSIS=$analysisID") or die $dbh->errstr; #à décommenter
		}
		my $fdrQueryStrg=($computeFDR)? ",DECOY=CONCAT(DECOY,',FDR=$maxFDR:$FDRalgo')" : '';
		$dbh->do("UPDATE ANALYSIS SET MIN_SCORE=$minScore$fdrQueryStrg WHERE ID_ANALYSIS=$analysisID");
	}
	elsif ($FDRalgo eq 'precomputed' && $maxFDR && ($percolatorInPD || $percolatorInMascot)) {
		&drawScoreDistribution($maxFDR,$analysisID,$dataFile);
		my $fdrQueryStrg=",DECOY=CONCAT(DECOY,',FDR=$maxFDR:$FDRalgo')";
		$dbh->do("UPDATE ANALYSIS SET MIN_SCORE=$minScore$fdrQueryStrg WHERE ID_ANALYSIS=$analysisID");
	}


	####>Creating Match groups<####
	&createMatchGroups;


	####>Protecting ID space in tables QUERY_VALIDATION and PROTEIN_VALIDATION<####
	#my ($maxQueryID,$maxProtValID);
	###>QUERY_VALIDATION
	#($queryID) = $dbh->selectrow_array("SELECT MAX(ID_QUERY) FROM QUERY_VALIDATION"); # fetching biggest ID in table
	#$maxQueryID=$queryID + scalar (keys %queryInfo) + 1; # protective ID: all IDs to be used +1 (prevents update on last real ID)
	#$dbh->do("INSERT INTO QUERY_VALIDATION (ID_QUERY,ID_ANALYSIS) VALUES ($maxQueryID,$analysisID)") || die $dbh->errstr(); # protecting ID space in table
	###>PROTEIN_VALIDATION
	#($protValID) = $dbh->selectrow_array("SELECT MAX(ID_PROT_VALID) FROM PROTEIN_VALIDATION");
	#$maxProtValID=$protValID + scalar (keys %maxProtMatch) + 1; # protective ID: all IDs to be used +1 (prevents update on last real ID)
	#$dbh->do("INSERT INTO PROTEIN_VALIDATION (ID_PROT_VALID,ID_ANALYSIS) VALUES ($maxProtValID,$analysisID)") || die $dbh->errstr(); # protecting ID space in table

	#$dbh->commit;

	####>Extracting data from databank file<####
	if (scalar keys %maxProtMatch) {
		if ($scanDB || (!$decoyFile && $fileFormat=~/\.PDM/)) {
			my $refProtMatch;
			if ($decoy) {
				my %trueProtMatch;
				foreach my $identifier (keys %maxProtMatch) {
					next if $identifier=~/^DECOY_/;
					$trueProtMatch{$identifier}=$maxProtMatch{$identifier}; # copying pointer!
				}
				$refProtMatch=\%trueProtMatch;
			}
			else {$refProtMatch=\%maxProtMatch;}
			#my $tempDbFile;

			if ($fileFormat=~/\.PDM/) {
				(my $tempDbFile="$promsPath{tmp}/$dataFile")=~s/\.pdm\Z/\.$userID.fasta/; # same dir as pdm file
				#if ($splitMsfFiles) {
				#	($tempDbFile="$promsPath{tmp}/$dataFile")=~s/\.pdm\Z/\.$fileList{$dataFile}{'file_id'}\.$userID.fasta/; # same dir as pdm file
				#}else{
				#
				#}
				&promsMod::getProtInfo('verbose',$dbh,$databankIDs[0],[$analysisID],\%protDes,\%protMW,\%protOrg,\%protLength,$refProtMatch,$tempDbFile); # only 1-db search allowed. $tempDbFile only for MSF
			}
			else {
				foreach my $dbID (@databankIDs) {
					#&promsMod::getProtInfo_old($dbh,$databankID,$analysisID,'verbose',\%protDes,\%protMW,\%protOrg,\%protLength,$refProtMatch,$tempDbFile); # $tempDbFile only for MSF
					&promsMod::getProtInfo('verbose',$dbh,$dbID,[$analysisID],\%protDes,\%protMW,\%protOrg,\%protLength,$refProtMatch);
					#unlink $tempDbFile if $tempDbFile;
				}
			}
		}
		elsif($fileFormat ne 'PARAGON.XML' && $fileFormat ne 'TDM.PEP.XML' && $fileFormat ne 'TDM.XML') {
			foreach my $identifier (keys %maxProtMatch) {
				$protDes{$identifier}='Data not imported';
				$protOrg{$identifier}='Data not imported';
				$protMW{$identifier}=0;
				$protLength{$identifier}=0;
			}
		}
	}


	############################################
	####>Loading data from file into tables<#### QUERY_VALIDATION, PROTEIN_VALIDATION and RANK_PROTEIN_MATCH
	############################################
	print "<FONT class=\"title3\">&nbsp;-Storing data into database:<BR>\n";
	my (@limitValue,$numEntry,$index,$counter1,$counter2,$maxCounter2);

	############Updating ANALYSIS TABLE IF NECESSARY#####################
	#unless ($decoyFile) {
	#	my @updateFields;
	#	push @updateFields,"WIFF_FILE='$phenyxSourceFile'" if $phenyxSourceFile;
	#	push @updateFields,"TAXONOMY='$phenyxTaxo'" if $phenyxTaxo;
	#	push @updateFields,"INSTRUMENT='$phenyxInstrument'" if $phenyxInstrument;
	#	push @updateFields,"VALID_STATUS=0" if $fileFormat=~/\.PDM/;
	#	if (scalar @updateFields) {
	#		my $queryStrg=join(',',@updateFields);
	#		$dbh->do("UPDATE ANALYSIS SET $queryStrg WHERE ID_ANALYSIS=$analysisID");
	#	}
	#}

	####<Loading parsed data into table QUERY_VALIDATION
	my $usedMaxRank=($msType eq 'MIS')? $maxRank : $MAX_RANK; # MAX_RANK for PMF & SQ
	print "&nbsp&nbsp;-Step 1 of 3:</FONT><BR>\n";
	my $pepString="INFO_PEP1";
	my $valueString='?';

	for (my $r=2;$r<=$usedMaxRank;$r++) {
		$pepString.=",INFO_PEP$r";
		$valueString.=',?';
	}

	my $sth1=$dbh->prepare("INSERT INTO QUERY_VALIDATION (ID_ANALYSIS,QUERY_NUM,EXT_SPECTRUMID,VALID_STATUS,MASS_DATA,MAX_SCORE,CHARGE,ELUTION_TIME,$pepString) VALUES ($analysisID,?,?,?,?,?,?,?,$valueString)");
	my $sth1bis=$dbh->prepare("INSERT INTO QUERY_MODIFICATION (ID_MODIFICATION,ID_QUERY,PEP_RANK,POS_STRING) VALUES (?,?,?,?)");

	$numEntry=scalar keys %queryInfo;
	undef @limitValue;
	foreach my $pc (@percent) {push @limitValue,int(0.5+($numEntry*$pc/100));}
	$index=0;
	$counter1=0;
	$counter2=0;

	$maxCounter2=int(0.5+($numEntry/100));
	print "<B>0%";
	foreach my $queryNum (sort{$a<=>$b} keys %queryInfo) {
		$counter1++;
		if ($counter1>=$limitValue[$index]) {
			print "$percent[$index]%"; # keeping connection alive
			$index++;
		}
		$counter2++;
		if ($counter2==$maxCounter2) {print '.'; $counter2=0;}
		next unless defined($maxQueryScore{$queryNum}); # skip valid status=-4! for performance and storage optimization

		my $numInsert=$#{$queryInfo{$queryNum}};
		my $queryValidStatus;


		if ($numInsert>=0 && $numValid{$queryNum}) {$queryValidStatus=-1;} # validatable interpretation(s)
		elsif ($numInsert>=0 && !$numValid{$queryNum}) {$queryValidStatus=-3;} # better scores exist for all interpretations
		else {$queryValidStatus=-4;} # no interpretations at all !!! OBSOLETE (skipped)

		my $massData="EXP=$massExp{$queryNum},OBS=$massObs{$queryNum}";
		for (my $i=$numInsert+1;$i<$usedMaxRank;$i++) {
			push @{$queryInfo{$queryNum}},undef;
		}
		#$queryID++;
# print ">id=$queryID: Q$queryNum SEL=$queryValidStatus mSC=$maxQueryScore{$queryNum} MASS=$massData PEP=@{$queryInfo{$queryNum}}<br>\n"; #,$mrExp{$queryNum},@{$queryInfo{$queryNum}}<br><br>\n";
		my $extSpecID=($extSpectrum{$queryNum})? $extSpectrum{$queryNum} : undef;
		my $elutTime=($elutionTime{$queryNum})? $elutionTime{$queryNum} : undef; # no elution time for PMF
		#my $extSpecID=($queryToSpectrum[$queryNum-1]) ? $queryToSpectrum[$queryNum-1] : undef;
		#$elutionTime{$queryNum}=($queryRT[$queryNum-1]) ? sprintf("%.2f",$queryRT[$queryNum-1]) : $elutionTime{$queryNum};
		#$charge{$queryNum}=($queryCharge[$queryNum-1]) ? $queryCharge[$queryNum-1] : $charge{$queryNum};

		$sth1->execute($queryNum,$extSpecID,$queryValidStatus,$massData,$maxQueryScore{$queryNum},$charge{$queryNum},$elutTime,@{$queryInfo{$queryNum}}) || die $sth1->errstr(); #à décommenter
		my $queryID = $dbh->last_insert_id(undef, undef, 'QUERY_VALIDATION', 'ID_QUERY');
		###> For modification table
		my $pepR=1;
		foreach my $pepInfo (@{$queryInfo{$queryNum}}) {
			next unless $pepInfo;
			if ($pepInfo=~ /VMOD=\s\+\s([^,]+)/) {
				my ($varModsString) = $1;
				my ($sequence) = ($pepInfo=~ /SEQ=(\w+)/);
				my @arrayVarMods = split (/\s\+\s/, $varModsString);
				my %modifs=();
				###> Read the all string
				foreach my $modif (@arrayVarMods) {
					my ($varModCode,$residues)=&promsMod::convertVarModString($modif);
					my @modifsPos;
					if ($modif=~ /:([\d\.]+)\)/) {
						@modifsPos=split(/\./,$1);
						#$modif=~ s/:[\d\.]+//;  # Update 10/06/2013 because Label:13C(6) (K:11) becomes LabelC(6) (K:11)
						$modif=~ s/:[\d\.]+\)/\)/; # Now Label:13C(6) (K:11) becomes Label:13C(6) (K)
					}
					#elsif ($residues =~ /\*/ || $residues =~ /\+/) {
					#	@modifsPos=(length($sequence));
					#}elsif ($residues =~ /=/ || $residues =~ /-/) {
					#	@modifsPos=(1);
					#}
					elsif ($residues =~ /(-|=|\+|\*)/) {
						@modifsPos =($1);
					}
					$modif=$varModCode if ($fileFormat eq 'PARAGON.XML' || $fileFormat eq 'TDM.PEP.XML' || $fileFormat eq 'TDM.XML');
					my $modificationID=($analysisMods{'V'}{$modif}{'numModDB'})? $analysisMods{'V'}{$modif}{'numModDB'} : ($analysisMods{'L'}{$modif}{'numModDB'}) ? $analysisMods{'L'}{$modif}{'numModDB'} : -1;
					foreach my $positionVMOD (@modifsPos) {
						$modifs{$modificationID}{$positionVMOD}=1;
					}
				}
				###> Put the Query_Modification into myProMS
				foreach my $modificationID (keys (%modifs)) {
					my $modifsPos=join('.',sort{if ($a=~/(-|=|\+|\*)\Z/ || $b=~/(-|=|\+|\*)\Z/) {return $a cmp $b} else {return $a<=>$b}} keys %{$modifs{$modificationID}});
					$modifsPos=-1 unless $modifsPos;
					$sth1bis->execute($modificationID,$queryID,$pepR,$modifsPos); #à décommenter
				}
			}
			$pepR++;
		}
	}
	$sth1->finish;
	$sth1bis->finish;

	print "</B></BR>\n";
# $dbh->disconnect(); print "END"; exit;

	####<Loading parsed data into table PROTEIN_VALIDATION
	print "<FONT class=\"title3\">&nbsp&nbsp;-Step 2 of 3:</FONT><BR>\n";
	my $sth2=$dbh->prepare("INSERT INTO PROTEIN_VALIDATION (ID_ANALYSIS,IDENTIFIER,DB_RANK,MW,PROT_LENGTH,PROT_DES,ORGANISM,NUM_MATCH,MAX_MATCH,SCORE,MAX_SCORE,SEL_STATUS,MATCH_GROUP) VALUES ($analysisID,?,?,?,?,?,?,0,?,?,?,-1,?)");
	$numEntry=scalar keys %maxProtMatch;
	undef @limitValue;
	foreach my $pc (@percent) {push @limitValue,int(0.5+($numEntry*$pc/100));}
	$index=0;
	$counter1=0;
	$counter2=0;
	$maxCounter2=int(0.5+($numEntry/100));
	print "<B>0%";
	foreach my $identifier (sort keys %maxProtMatch) {
		#next if $identifier=~/^DECOY_/; # skip decoy proteins
		#$protValID++;
		my $des=&promsMod::resize($protDes{$identifier},250); # max length allowed in table
		my $organism=&promsMod::resize($protOrg{$identifier},100); # max length allowed in table
# 		$protLength{$identifier}=0 unless $protLength{$identifier};
		my $score=($msType eq 'PMF')? $maxProtScore{$identifier} : 0;
# print ">($protValID) $identifier:<br>DES=$des<br>ORG=$organism<br>LEN=$protLength{$identifier} aa<br>MW=$protMW{$identifier}<br>MGr=$matchGroup{$identifier}<br>SC=$score/$maxProtScore{$identifier}<br>mMATCH=",scalar keys %{$maxProtMatch{$identifier}},"<BR>\n";

		$sth2->execute($identifier,$protDbRank{$identifier},$protMW{$identifier},$protLength{$identifier},$des,$organism,scalar keys %{$maxProtMatch{$identifier}},$score,$maxProtScore{$identifier},$matchGroup{$identifier}) || die $sth2->errstr(); #,$protSelStatus{$identifier}
		$counter1++;
		if ($counter1>=$limitValue[$index]) {
			print "$percent[$index]%"; # keeping connection alive
			$index++;
		}
		$counter2++;
		if ($counter2==$maxCounter2) {print '.'; $counter2=0;}
	}
	$sth2->finish;
	print "</B><BR>\n";
# $dbh->disconnect(); print "END"; exit;




	####<Loading parsed data into table RANK_PROTEIN_MATCH
	print "<FONT class=\"title3\">&nbsp&nbsp;-Step 3 of 3:</FONT><BR>\n";
	my $sth3=$dbh->prepare("INSERT INTO RANK_PROTEIN_MATCH (ID_ANALYSIS,QUERY_NUM,PEP_RANK,IDENTIFIER,MATCH_MULTI,MATCH_INFO) VALUES ($analysisID,?,?,?,?,?)");
	$numEntry=scalar keys %rankProtMatch;
	undef @limitValue;
	foreach my $pc (@percent) {push @limitValue,int(0.5+($pc*$numEntry/100));}
	$index=0;
	$counter1=0;
	$counter2=0;
	$maxCounter2=int(0.5+($numEntry/100));
	print "<B>0%";
	foreach my $rankID (keys %rankProtMatch) {
		my ($queryNum,$rank)=split(/:/,$rankID);
		foreach my $identifier (keys %{$rankProtMatch{$rankID}}) {
			my $numMatch=scalar @{$rankProtMatch{$rankID}{$identifier}};
			my $begString=join ':',@{$rankProtMatch{$rankID}{$identifier}};
			if (length($begString) > 255) { # Max size of DB field
				$numMatch=0;
				foreach my $beg (@{$rankProtMatch{$rankID}{$identifier}}) {
					if ($numMatch==0) {
						$begString="$beg";
						$numMatch=1;
						next;
					}
					last if length("$begString:$beg") > 255; # Max size of DB field
					$numMatch++;
					$begString.=":$beg";
				}
			}
			$sth3->execute($queryNum,$rank,$identifier,$numMatch,$begString) || die $sth3->errstr(); #à décommenter
		}
		$counter1++;
		if ($counter1>=$limitValue[$index]) {
			print "$percent[$index]%"; # keeping connection alive
			$index++;
		}
		$counter2++;
		if ($counter2==$maxCounter2) {print '.'; $counter2=0;}
	}
	$sth3->finish;
	print "</B><BR><FONT class=\"title3\">Done.</FONT><BR>\n";




	########################
	####>Ending process
	########################

	###########Updating ANALYSIS TABLE IF NECESSARY#####################
	unless ($decoyFile) {
		my @updateFields;
		push @updateFields,"WIFF_FILE='$phenyxSourceFile'" if $phenyxSourceFile;
		push @updateFields,"TAXONOMY='$phenyxTaxo'" if $phenyxTaxo;
		push @updateFields,"INSTRUMENT='$phenyxInstrument'" if $phenyxInstrument;
		if ($validStatus==-1 && ($numEntry==0 || $fileFormat=~/\.PDM/)) {
			push @updateFields,"VALID_STATUS=0";
			$validStatus=0;
		}
		if (scalar @updateFields) {
			my $queryStrg=join(',',@updateFields);
			$dbh->do("UPDATE ANALYSIS SET $queryStrg WHERE ID_ANALYSIS=$analysisID");
		}
		push @databankScans,$analysisID if ($validStatus==-1 && $fileFormat!~/TDM/);
		if($fileFormat=~/TDM/){
			$fileFormat.=($isSpectra eq 'yes' && $fileFormat eq 'TDM.XML')? ':Spec' : ($fileFormat eq 'TDM.XML') ? ':NoSpec' : '';
			$dbh->do("UPDATE ANALYSIS SET VALID_STATUS=0, FILE_FORMAT=\"$fileFormat\" WHERE ID_ANALYSIS=$analysisID");
		}
		#push @databankScans,$analysisID if $validStatus==-1;
	}

	####>Updating valid status if no proteins<####
	#if (!$decoyFile && $numEntry==0 && $validStatus==-1) { # no proteins
	#	$dbh->do("UPDATE ANALYSIS SET VALID_STATUS=0 WHERE ID_ANALYSIS=$analysisID");
	#}

	####>Deleting protective IDs<####
	#$dbh->do("DELETE FROM QUERY_VALIDATION WHERE ID_QUERY=$maxQueryID") || die $dbh->errstr();
	#$dbh->do("DELETE FROM PROTEIN_VALIDATION WHERE ID_PROT_VALID=$maxProtValID") || die $dbh->errstr();

	# print "OK load valid<BR>\n"; $dbh->disconnect(); exit;
	$dbh->commit;


	####>Deleting uploaded file<####
	unlink $fileList{$dataFile}{'file_path'} if $deleteFiles; # true only if $action eq 'anaBatch';
	if ($deleteFiles && $percolatorInMascot){
		unlink "$promsPath{tmp}/batch/$userID/.percolator/$dataFile.decoy.pop";
		unlink "$promsPath{tmp}/batch/$userID/.percolator/$dataFile.target.pop";
	}

	####>Moving qVality file<####
	if ($action eq 'reval') {
		if (-e "$promsPath{peptide}/proj_$projectID/ana_$analysisID/scores.png") {
			move("$promsPath{peptide}/proj_$projectID/ana_$analysisID/scores.png","$promsPath{valid}/ana_$analysisID/scores.png")
		}
	}

	print "<BR>\n";

}


##################################
####>Deleting temporary files<####
##################################
foreach my $tempFile (@temporaryFiles) {
	unlink $tempFile;
}
unlink glob "$promsPath{tmp}/*.$userID.fasta" if $action eq 'anaBatch';

#################################
####>Deleting original files<####
#################################
if ($deleteFiles && $action eq 'anaBatch') {
	foreach my $dataFile (keys %fileList) {
		if ($dataFile=~/(.+)_\d+\.pdm/) { # msf file
			unlink "$filePath/$1.msf";
		}
		else {
			unlink "$filePath/$dataFile";
		}
	}
}

########################################
####>Fetching info for frame reload<####
########################################
my @anaInfo=&promsMod::getItemInfo($dbh,'ANALYSIS',$firstAnalysisID) if $firstAnalysisID;

$dbh->disconnect; # safer to disconnect before fork

################################################
####>Launching databank scans in background<####
################################################
if (scalar @databankScans) {
	#my $analysisStrg=join(',',@databankScans);
	#system "./scanDatabank.pl $userID $analysisStrg";

	my $analysisStrg=join(',',@databankScans);

	print qq|<SCRIPT LANGUAGE="JavaScript">systemFrame.location="./storeAnalyses.cgi?ACT=scanDB&ANA_ID=$analysisStrg";</SCRIPT>|;
	my $endFile="$promsPath{logs}/scan_$userID.end";
	my $numCycles=0;
	while ($numCycles < 15 && !-e $endFile) { # 30 sec wait max
		sleep 2;
		$numCycles++;
	}
	sleep 1;
	unlink $endFile if -e $endFile;
}


##########################
####>Reloading frames<####
##########################
#exit; # DEBUG!
if ($firstAnalysisID) {
	print "<FONT class=\"title2\">Import is finished !</FONT>\n";
	print "<BR><FONT class=\"title3\">Protein annotations are being imported in background.</FONT>\n" if scalar @databankScans;
	sleep 3;

	my $anaBranchID="analysis:$firstAnalysisID";
	my $updateURL="parent.navFrame.location=\"$promsPath{cgi}/openProject.cgi?ID=$projectID&ACT=nav";

	my ($treeFrame,$targetTree);
	if ($anaInfo[2]{'ITEM'} eq 'SAMPLE') {
		$updateURL.="&branchID=$anaBranchID\"";
		$treeFrame='navFrame';
		$targetTree='projectTreeStatus';
	}
	else {
		$updateURL.="&branchID=".lc($anaInfo[2]{ITEM}).":$anaInfo[2]{ID}&itemBranchID=$anaBranchID\"";
		$treeFrame='itemFrame';
		$targetTree='itemTreeStatus';
	}

	print qq
|<SCRIPT LANGUAGE="JavaScript">
	top.promsFrame.selectedAction='summary';
	if ('$action' != 'reval' && top.$targetTree) {top.$targetTree=parent.$treeFrame.addItemToTree('$anaBranchID');}
	$updateURL;
</SCRIPT>
|;
}
else { # file formet error
	print "<FONT class=\"title2\">&nbsp;&nbsp;File(s) not imported due to unexpected error(s).</FONT>\n";
}

print "</BODY>\n</HTML>\n";

exit;


####################<<< SUBROUTINES >>>#########################



#################################################################
####<Parsing Mascot .dat file if creating a new PMF Analysis>####
#################################################################
sub parseMascotDAT_PMF {
	my ($msFileName,$mixedSearch)=@_; #,$refFiledata
	$mixedSearch=0 unless $mixedSearch;
	####<Parsing data from $msFileName into globals %massExp, %massObs, %queryInfo, %maxQueryScore, %numValid, %rankProtMatch, %maxProtMatch, %maxProtScore, %matchList, %matchGroup %spectralCount
	my %deltaMasses; # required to sort ranks properly
	my %tempRankProtMatch; # stores macthes with default temp ranks
	my %querySeq; # list of sequences+coded arMods=>(ranks,numProt) for same query
	#my %varMods; # list of possible variable modifications
	my %codedModSeq2varMods; # stores correspondance between pepSeq:codedVarMod and VarModStrg
	my %queryRankSeq; # list of query:Rank=>sequences+varMods
	my %bestDeltaMass; # required to record which query:rank best matches and interpretation (defined in PMF by min delta mass, not best score!)
	print "<FONT class=\"title3\">&nbsp;-Extracting PMF data from $msFileName...";

	my ($identifier,$queryNum,$tempRank);
	open (DATAFILE, $fileList{$msFileName}{'file_path'});
	#foreach my $line (@{$refFiledata}) { #}
	while (my $line=<DATAFILE>) {
		last if $line=~/name="(mixture|peptides)"/; # mixture for PMF, peptides for MIS

		###<Possible fix modification (parameters section)
		#if ($line=~/^MODS=(.+)\n/) {
		#	(my $modString=$1)=~s/\s*\Z//; # trailing \r !!!?
		#	foreach my $mod (split(/,/,$modString)) {
		#		my ($fixModName,$specificity)=&promsMod::convertVarModString($mod);
		#		$analysisMods{'F'}{$mod}{'specificity'}=$specificity;
		#		$analysisMods{'F'}{$mod}{'numModDB'}=&promsMod::getModificationIDfromString($dbh,$fixModName,$specificity,\%vmodsInUnimod);
		#	}
		#}
		####<Possible variable modification (parameters section)
		#elsif ($line=~/^IT_MODS=(.+)\n/) { # done twice if SQ search (it's OK or make %varMods global + reset between analyses)
		#	my $numMod=1;
		#	foreach my $mod (split(/,/,$1)) {
		#		#$mod=~s/^_+//; # __Oxidation -> Oxidation !!! not a normal Oxydation (PP 07/05/12)!!
		#		$varMods{$numMod}=$mod;
		#		my ($varModName,$specificity)=&promsMod::convertVarModString($mod);
		#		$analysisMods{'V'}{$mod}{'specificity'}=$specificity;
		#		$analysisMods{'V'}{$mod}{'numModDB'}=&promsMod::getModificationIDfromString($dbh,$varModName,$specificity,\%vmodsInUnimod);
		#
		#		$numMod++;
		#	}
		#	next;
		#}

		###<Mass info (summary section)
		unless ($mixedSearch) { # already done in MIS section
			if ($line=~/^qmass(\d+)=(\S+)/) {$massExp{$1}=$2;} # Mr(exp)
			elsif ($line=~/^qexp(\d+)=(\S+),(\d+)\+/) {
				$massObs{$1}=$2; # Observed
				$charge{$1}=$3; # charge
			}
			#next;
		}

		###<Match info for MS (PMF, summary section)
 		if ($line=~/^h\d+[=_]/) {
	 		if ($line=~/^h\d+=(\S+),(\S+),\S+,\S+/) {
		 		$identifier=$1;
		 		$maxProtScore{$identifier}=$2 unless $maxProtScore{$identifier}; # if already computed in MIS section (test MIS)
		 		#$protSelStatus{$identifier}=-1;
	 		}
	 		elsif ($line=~/^h\d+_q(\d+)=(\S+)/) {
		 		$queryNum=$1;
				next if $matchedQueryMIS{$queryNum}; # already matched in MIS section (test MIS)
				$maxQueryScore{$queryNum}=0;
		 		if ($2 eq '-1') { # hP_qQ=-1 no hit for query Q
					@{$queryInfo{$queryNum}}=() unless $queryInfo{$queryNum}; # set to empty if never encountered
				}
				else {
					chomp($line);
				 	my @data=split(/[=,]/,$line);

#					if ($data[10]>0) { # score > 0!!! for SQ search if interpretation is also found in MIS section (but not skipped in test MIS because < min score)
#						@{$queryInfo{$queryNum}}=() unless $queryInfo{$queryNum}; # set to empty if never encountered
#						$matchedQueryMIS{$queryNum}=1;
#print "*$queryNum ($data[10])*";
#						next;
#					}

					my $codedModSeq="$data[7]:$data[9]"; # $data[9] = coded varMod
				 	$tempRank=(!$querySeq{$queryNum})? 1 : ($querySeq{$queryNum}{$codedModSeq})? $querySeq{$queryNum}{$codedModSeq} : 1+scalar (keys %{$querySeq{$queryNum}});
				 	@{$tempRankProtMatch{"$queryNum:$tempRank"}{$identifier}}=($data[4]); # beg !!!!! ONLY 1 beg !!!!!!!!!
			 		#$maxProtMatch{$identifier}{$data[7]}++; # seq
			 		#$matchList{$identifier}{$data[7]}=1; # seq
			 		if (!$querySeq{$queryNum} || !$querySeq{$queryNum}{$codedModSeq}) { # new seq for query
						$querySeq{$queryNum}{$codedModSeq}=$tempRank; # 1st time seq is encountered
			 			#$numValid{$queryNum}++;
						$deltaMasses{$queryNum}{$tempRank}=$data[3];

				 		##<Variable modifications ($data[9]): list of 0..numMod: NterpepSeqCter eq 10002000
						my $varModString='';
						if ($codedModSeq2varMods{$codedModSeq}) {
							$varModString=($codedModSeq2varMods{$codedModSeq} eq 'none')? '' : $codedModSeq2varMods{$codedModSeq};
						}
						elsif ($data[9]=~/[^0]/) { # new $codedModSeq with a varMod
							foreach my $numMod (sort{$a<=>$b} keys %varMods) {
								my $position='';
								my $currentPosition=0;
								if ($data[9]=~/$numMod/) {
									my $currentVarMod=$varMods{$numMod};
									my @pos;
									while ($currentPosition != -1) {
										$currentPosition=index ($data[9],$numMod,$currentPosition+1);
										push @pos,$currentPosition if $currentPosition >= 0;
									}
									if (scalar @pos) {
										my $posString=':'.join('.', @pos);
										$currentVarMod=~s/\(([^\(]+)\)\Z/\($1$posString\)/;
									}
									$varModString.=" + $currentVarMod";
								}
							}
							$codedModSeq2varMods{$codedModSeq}=$varModString;
						}
						else {  # new $codedModSeq with no varMod
							$codedModSeq2varMods{$codedModSeq}='none';
							$varModString=''; # just to be safe
						}
						my $seqVmod=$data[7];
						$seqVmod.=$varModString if $varModString;
						#my $select=($matchedSeqMIS{$seqVmod})? -1 : 0; # priority to MIS
						my $absDeltaMass=abs($data[3]);
						my $select;
						if ($matchedSeqMIS{$seqVmod}) {$select=-1} # priority to MIS
						elsif (!$bestDeltaMass{$seqVmod}) {
							$select=0;
							@{$bestDeltaMass{$seqVmod}}=($absDeltaMass,$queryNum,$tempRank);
						}
						elsif ($absDeltaMass >= $bestDeltaMass{$seqVmod}[0]) {$select=-1}
						else { # smaller delta found
							$select=0;
							$queryInfo{$bestDeltaMass{$seqVmod}[1]}[$bestDeltaMass{$seqVmod}[2]-1]=~s/SEL=0,LSP=0/SEL=-1,LSP=1/; # set old best at -1
							@{$bestDeltaMass{$seqVmod}}=($absDeltaMass,$queryNum,$tempRank);
						}
						my $LSPstrg=($select==-1)? ',LSP=1' : ',LSP=0'; # Flag for lower-scoring peptide
						my $dataString="SEL=$select$LSPstrg,MIS=$data[1],CALC=$data[2],DELT=$data[3],SEQ=$data[7],SC=0,VMOD=$varModString,";
						#$varModString='' unless $varModString;
						push @{$queryInfo{$queryNum}},$dataString;
						#if ($spectralCount{"$data[7]:$varModString"}) {
						#	$spectralCount{"$data[7]:$varModString"}++;
						#}else{
						#	$spectralCount{"$data[7]:$varModString"}=1;
						#}
						$spectralCount{$seqVmod}++; # self initialization if not defined
					}
		 		}
	 		}
 			elsif ($line=~/^h\d+_q\d+_terms=(\S+)/) { # update match info to include boundary AA
				next if $matchedQueryMIS{$queryNum}; # already matched in MIS section
				my @boundAA=split(/:/,$1); # !!!! there should be only 1 set !!!!! unlike MIS
				$tempRankProtMatch{"$queryNum:$tempRank"}{$identifier}[0].=",$boundAA[0]";
			}
 		}
	}
	close DATAFILE;

	####<Computing MATCH for all peptides
	foreach my $queryNum (keys %queryInfo) {
		next if $matchedQueryMIS{$queryNum}; # already matched in MIS section
		next unless scalar @{$queryInfo{$queryNum}};
		foreach my $rank (values %{$querySeq{$queryNum}}) { # seq=>rank
			$queryInfo{$queryNum}[$rank-1].="MATCH=".scalar (keys %{$tempRankProtMatch{"$queryNum:$rank"}}).",";
			my ($seqVmod)=($queryInfo{$queryNum}[$rank-1]=~ /SEQ=(\w+)/);
			if ($queryInfo{$queryNum}[$rank-1]=~ /VMOD=([^,]+)/) {
				$seqVmod.=$1;
			}
			$queryRankSeq{"$queryNum:$rank"}=$seqVmod;
			$queryInfo{$queryNum}[$rank-1].="SPC=$spectralCount{$seqVmod},";
		}
	}

	####<Re-sorting ranks according to increasing mass delta (best has smallest delta)
	####<Converting tempRankProtMatch to rankProtMatch
	foreach my $queryNum (keys %queryInfo) {
		next if $matchedQueryMIS{$queryNum}; # already matched in MIS section
		next unless scalar @{$queryInfo{$queryNum}};
		if (scalar @{$queryInfo{$queryNum}}==1) {
			$numValid{$queryNum}++ if $queryInfo{$queryNum}[0]=~/SEL=0/;
			%{$rankProtMatch{"$queryNum:1"}}=%{$tempRankProtMatch{"$queryNum:1"}};
			next; # no need to resort queryInfo if only 1 rank
		}
		else {
			my @tempInfo;
			my $newRank=0;
			foreach my $oldRank (sort{abs($deltaMasses{$queryNum}{$a})<=>abs($deltaMasses{$queryNum}{$b})} keys %{$deltaMasses{$queryNum}}) {
				$newRank++;
				last if $newRank > $MAX_RANK; # discard ranks over absolute max Rank
				push @tempInfo,$queryInfo{$queryNum}[$oldRank-1];
				$numValid{$queryNum}++ if $queryInfo{$queryNum}[$oldRank-1]=~/SEL=0/;
				%{$rankProtMatch{"$queryNum:$newRank"}}=%{$tempRankProtMatch{"$queryNum:$oldRank"}};
			}
			@{$queryInfo{$queryNum}}=@tempInfo;
		}
		##<Computing peptide matches on proteins
		foreach my $queryRank (keys %rankProtMatch) {
			if ($mixedSearch) {
				my ($qNum,$rank)=split(/:/,$queryRank);
				next if $matchedQueryMIS{$qNum};
			}
			foreach my $identifier (keys %{$rankProtMatch{$queryRank}}) {
				$maxProtMatch{$identifier}{$queryRankSeq{$queryRank}}++; # seq+vMod
				$matchList{$identifier}{$queryRankSeq{$queryRank}}=1; # seq+vMod
			}
		}
	}

	####<Creating Match Groups
	#&createMatchGroups;

	print " Done.</FONT><BR>\n";
}

###################################################################
####<Parsing Mascot .dat file if creating a new MS/MS Analysis>####
###################################################################
sub parseMascotDAT_MIS {
	my ($msFileName,$mixedSearch,$percolatorInMascot)=@_; #,$refFiledata
	$mixedSearch=0 unless $mixedSearch;
	####<Parsing data from $msFileName into globals %massExp, %massObs, %queryInfo, %maxQueryScore, %numValid, %rankProtMatch, %maxProtMatch, %maxProtScore, %matchList, %matchGroup,  %elutionTime
	my (%bestScore,%matchedSeq); # locals required for calculation of protein max score
	my %seqQueryRank; # local required for setting interpretation at SEL=-1 when a better score is found
	my (%seqQueryRankBad,%bestScoreBad,%tempQueryInfo); # idem as above for mixed decoy DB
	#my $maxQueryNum=0;
	#my %matchQuality; # records rank and number of matching ranks (to exclude protein if matched by 1 maxRank)
	#my %varMods; # list of possible variable modifications
	my %begPosition; # stores list of match position foreach protein matched by given peptide
# 	my $onNoQueryNumber; # n° de query, needed for extrating elution Time
	my $currentQuery; # currently scanned query, needed for extracting elution Time
	my @databankRanks; # for multi-db search
	my ($prevQueryScore,$scoreRank); # needed for scoreRank count (22/07/10)
	my $section="";
	my ($qSign,$sFlag,$iFlag); # auto/file DECOY managment (not for mixed databank)
	my ($rankGood,$rankBad); # mixed decoy databank (22/07/10)
	my ($matchGood,$matchBad); # mixed decoy databank (made multi-line global 06/05/14)
	my %percolatorMascotScores=();
	&getPercolatorScores($msFileName,\%percolatorMascotScores) if $percolatorInMascot;

	print "<FONT class=\"title3\">&nbsp;-Extracting MIS data from $msFileName...";

	open (DATAFILE, $fileList{$msFileName}{'file_path'});
	#foreach my $line (@{$refFiledata}) { #}
	my $nbLine=0;
	while (my $line=<DATAFILE>) {
		#  last if $line=~/name="proteins"/;
		$nbLine++;
		if ($line =~ /gc0p4Jq0M2Yt08jU534c0p/) {
			$section = " ";
			($qSign,$sFlag,$iFlag)=($decoyFile)? ('-','_','DECOY_') : ('','','');
			next;
		}
		if ($section eq " " && $line =~ /name="(.+)"/) {
			$section=$1;
			if ($section =~ /query(\d+)/) {
				$currentQuery=$1;
				$section='query';
			}
			elsif ($section eq 'decoy_summary') {$qSign='-';}
			elsif ($section eq 'decoy_peptides') {($qSign,$sFlag,$iFlag)=('-','_','DECOY_');}
			next;
		}

		if ($section eq 'parameters') {
			####<Possible fix modification (parameters section)
			#if ($line=~/^MODS=(.+)\n/) {
			#	(my $modString=$1)=~s/\s*\Z//; # trailing \r !!!?
			#	foreach my $mod (split(/,/,$modString)) {
			#		my ($fixModName,$specificity)=&promsMod::convertVarModString($mod);
			#		$analysisMods{'F'}{$mod}{'specificity'}=$specificity;
			#		$analysisMods{'F'}{$mod}{'numModDB'}=&promsMod::getModificationIDfromString($dbh,$fixModName,$specificity,\%vmodsInUnimod);
			#	}
			#}
			####<Possible variable modification (parameters section)
			#elsif ($line=~/^IT_MODS=(.+)\n/) {
			#	(my $modString=$1)=~s/\s*\Z//; # trailing \r !!!?
			#	my $numMod=1;
			#	foreach my $mod (split(/,/,$modString)) {
			#		#$mod=~s/^_+//; # __Oxidation -> Oxidation !!! not a normal Oxydation (PP 07/05/12)!!
			#		$varMods{$numMod}=$mod;
			#		$analysisMods{'V'}{$mod}{'numMod'}=$numMod; # Acetyl (K)
			#		my ($varModName,$specificity)=&promsMod::convertVarModString($mod);
			#		$analysisMods{'V'}{$mod}{'specificity'}=$specificity;
			#		$analysisMods{'V'}{$mod}{'numModDB'}=&promsMod::getModificationIDfromString($dbh,$varModName,$specificity,\%vmodsInUnimod);
			#		$numMod++;
			#	}
			#}
			#els
			if ($line=~ /SEARCH_ALGO/ ) {
				# Only fake .DAT generated from MSF files have a SEARCH_ALGO line. It's a tag.
				# Update the file with XXX.ana to add the analysis
				my $anaFile = $msFileName;
				$anaFile =~ s/\_[0-9]*\.*[0-9]*\.pdm/\.ana/g;# reconstruct the name of the MSF file and change the extension
				open (ANAFILE, ">>$promsPath{valid}/multi_ana/proj_$projectID/$anaFile");
				print ANAFILE "$analysisID\n";
				close ANAFILE;
				(my $msfFile = $anaFile) =~ s/\.ana/\.msf/;
				my $dbsqlite = DBI->connect( "dbi:SQLite:$promsPath{valid}/multi_ana/proj_$projectID/$msfFile", "", "", {PrintError => 1,RaiseError => 1});
# COMMENTED by PP ----------->
				#my $getSpectrumIDs = $dbsqlite->prepare("SELECT SpectrumID, RetentionTime, MassPeaks.Charge FROM SpectrumHeaders, MassPeaks, Spectra WHERE Spectra.UniqueSpectrumID=SpectrumHeaders.UniqueSpectrumID AND SpectrumHeaders.MassPeakID=MassPeaks.MassPeakID ORDER BY SpectrumHeaders.Mass ASC");
				#$getSpectrumIDs->execute;
				#while (my ($spectrumID,$retentionTime,$charge) = $getSpectrumIDs->fetchrow_array) {
				#	push @queryToSpectrum, $spectrumID;
				#	push @queryRT, $retentionTime;#Already in minutes
				#	push @queryCharge, $charge;
				#}
				#my $sthSpec= $dbsqlite->prepare("SELECT SpectrumID, RetentionTime, ScanNumbers FROM SpectrumHeaders ORDER BY Charge ASC, Mass ASC");
				my ($proteomeDiscovererVersion)=$dbsqlite->selectrow_array("SELECT SoftwareVersion FROM SchemaInfo ORDER BY rowid ASC LIMIT 1"); # last record
				$proteomeDiscovererVersion=~s/^(\d\.\d+).*/$1/; # x.x.x.xx -> x.x (numerical value)
				my ($specTableN)=($proteomeDiscovererVersion >=2.2)?"MSnSpectrumInfo":"SpectrumHeaders";
				my $sthSpec= $dbsqlite->prepare("SELECT SpectrumID, RetentionTime, ScanNumbers FROM $specTableN ORDER BY Mass ASC, Charge ASC");
				$sthSpec->execute;
				my $queryNum=0;
				while (my ($spectrumID,$retentionTime,$scanNumber) = $sthSpec->fetchrow_array) {
					$queryNum++;
					$extSpectrum{$queryNum}=$extSpectrum{"-$queryNum"}=$spectrumID; # in case of decoy data
					$elutionTime{$queryNum}=$elutionTime{"-$queryNum"}=sprintf "et%.2f;sp$spectrumID;sc$scanNumber;",$retentionTime; #Already in minutes
				}
				$sthSpec->finish;
				$dbsqlite->disconnect;
			}

		}

# 		if ($line=~/^IT_MODS=(.+)\n/) {
# 			my $numMod=1;
# 			my $lineMod = $1 ;
# 			chomp ($lineMod);
# 			foreach my $mod (split(/,/,$lineMod)) {
# 				if ($mod=~/(.+\))/) {$varMods{$numMod}=$1;  	}
# 				else {$varMods{$numMod}=$mod;  	} # warning, may include new line
# 				$numMod++;
# 			}
# 		}


		elsif ($section eq 'summary' || (!$decoyFile && !$dbMixedDecoyFlag && $section eq 'decoy_summary')) { # skip decoy section if whole file is decoy
			###<Mass info (summary section)
			if ($line=~/^qmass(\d+)=(\S+)/) { # Mr(exp)
				$massExp{"$qSign$1"}=$2;
				$massExp{"-$1"}=$2 if $dbMixedDecoyFlag;
			}
			#elsif ($line=~/^qexp(\d+)=(\S+),/) {$massObs{"$qSign$1"}=$2;} # Observed
			elsif ($line=~/^qexp(\d+)=(\S+),(\d+)\+/) {
				$massObs{"$qSign$1"}=$2; # Observed
				$charge{"$qSign$1"}=$3; # charge
				if ($dbMixedDecoyFlag) {
					$massObs{"-$1"}=$2; # Observed
					$charge{"-$1"}=$3; # charge
				}
			}
		}

		elsif ($section eq 'peptides' || (!$decoyFile && !$dbMixedDecoyFlag && $section eq 'decoy_peptides')) { # skip decoy section if whole file is decoy or mixed decoy databank

			###<Match info for MS/MS (MIS, peptides section)
	 		if ($line=~/^q(\d+)_p(\d+)/) {

				my ($queryNum,$rank,$rankID)=("$qSign$1",$2,"$qSign$1:$2"); # qSign: decoy sign
		#$maxQueryNum=$queryNum if $maxQueryNum<$queryNum;  # only for mixed decoy DB
				next if $rank > $maxRank;
				next if ($percolatorInMascot && !$percolatorMascotScores{$section}{abs($queryNum)});
		#next if ($scoreRank && $scoreRank > $maxRank); # ??????????
				chomp $line;

				if ($line=~/^q\d+_p\d+_db=(\d+)/) { # multi-db search: eg. 010202010101
					my $dbStrg=$1;
					@databankRanks=($dbStrg=~ m/../g); # split every 2 characters -> (01,02,02,01,01,01)
					foreach my $i (0..$#databankRanks) {$databankRanks[$i]=~s/^0//;} # remove starting 0
				}
				elsif ($line=~/^q\d+_p\d+=/) {
					if ($rank==1) {
						#$maxQueryScore{$queryNum}=0; # commented 12/12/12 PP
						if ($dbMixedDecoyFlag) {
							%{$tempQueryInfo{$queryNum}}=();
							%{$tempQueryInfo{"-$queryNum"}}=();
							$maxQueryScore{"-$queryNum"}=0;
						}
						else {
							@{$queryInfo{$queryNum}}=();
						}
					}
	($matchGood,$matchBad)=(0,0);
					undef %begPosition; # reset for new peptide
					my ($pepData,$protData)=split(/;/,$line);

					if ($protData) { # qQ_pR != -1 => @data will be complete
						my @data=split(/[=,]/,$pepData);
						$matchedQueryMIS{$queryNum}=1; # query is used for MS/MS. Recorded even if score < min score (needed in PMF section if SQ search)

						next if ($data[8]==0 || $data[8] < $minScore); # $data[8] = peptide score (skip score==0!!!)

						##<Variable modifications ($data[7]): list of 0..numMod: NterpepSeqCter eq 10002000
						my $varModString='';
						# my $varPosModString ;
						foreach my $numMod (sort{$a<=>$b} keys %varMods) {
							my $position='';
							my $currentPosition=0;
							if ($data[7]=~/$numMod/) {
								if ($varMods{$numMod} =~/\+/) { # double modification on the same residue
										my @pos;
										while ($currentPosition != -1) {
											$currentPosition=index ($data[7],$numMod,$currentPosition+1);
											push @pos,$currentPosition if $currentPosition >= 0 ;
										}
										my $currentVarMod=$varMods{$numMod};
										my @singleVarMod=split(/\+| /,$currentVarMod);
										for(my $i =0; $i < $#singleVarMod ; $i++){
											$currentVarMod=$singleVarMod[$i]." ".$singleVarMod[$#singleVarMod];
											if (scalar @pos) {
												my $posString=':'.join('.', @pos);
										    	$currentVarMod=~s/\(([^\(]+)\)\Z/\($1$posString\)/;
										 	}
											$varModString.=" + $currentVarMod";
										}
								}
								else{
										my $currentVarMod=$varMods{$numMod};
										my @pos;
										while ($currentPosition != -1) {
											$currentPosition=index ($data[7],$numMod,$currentPosition+1);
											push @pos,$currentPosition if $currentPosition >= 0 ;
										}
										if (scalar @pos) {
											my $posString=':'.join('.', @pos);
											$currentVarMod=~s/\(([^\(]+)\)\Z/\($1$posString\)/;
										}
										$varModString.=" + $currentVarMod";
								}

							}
						}
						my $seqVmod = "$sFlag$data[5]";  # sFlag: decoy flag
						$seqVmod.=$varModString if $varModString;
						my $varModDataStrg=($varModString)? ",VMOD=$varModString" : ""; # ',' before!
						$matchedSeqMIS{$seqVmod}=1;

						###<Peptide selection
						my @protList=split(/,/,$protData);

						if ($rank==1) {$rankGood=$rankBad=0;} # reset for new query
						#my ($matchGood,$matchBad)=(0,0);
						if ($dbMixedDecoyFlag) {
							foreach my $identData (@protList) {
								my ($identifier)=($identData=~/"([^"]+)"/);
								if ($identifier=~/^$dbMixedDecoyFlag/) {$matchBad=1;} else {$matchGood=1;}
								last if ($matchGood && $matchBad);
							}
							$rankGood++ if $matchGood;
							$rankBad++ if $matchBad;
						}
						else { # all other cases
							$matchGood=1;
							$rankGood=$rank;
						}
#$matchedQueryMIS{$queryNum}=1 if $matchGood==1; # needed in PMF section if SQ search

						my (%numProt,%numFlaggedProt);
						for (my $m=0;$m<=$#protList;$m++) { # mixed & other DB
							my $usedRankID;
							my ($identifier,$beg)=($protList[$m]=~/"([^"]+)":\d:(\d+)/); # meaning of Value "multi" in data file is unclear !!!!!
							if ($dbMixedDecoyFlag && $identifier=~/^$dbMixedDecoyFlag/) { # eg. REV_Q9Z0P5
								$identifier="DECOY_$identifier"; # tagged as decoy (dbMixedDecoyFlag remains)
								$usedRankID="-$queryNum:$rankBad";
								$numFlaggedProt{$identifier}=1;
							}
							else {
								$identifier="$iFlag$identifier"; # iFlag: decoy flag
								$usedRankID="$queryNum:$rankGood";
								$numProt{$identifier}=1;
							}

							$matchList{$identifier}{$seqVmod}=1; # proteins can be repeated
							if (scalar keys %{$rankProtMatch{$usedRankID}}) { # $rankID
								my $badBeg=0;
								foreach my $prevBeg (@{$rankProtMatch{$usedRankID}{$identifier}}) {
									if ($beg==$prevBeg) {$badBeg=1; last;}
								}
								next if $badBeg; # beg already recorded (duplicated in data file!)
							}
			#if ($thisProtIsDecoy) {
			#	push @{$rankProtMatch{"-$rankID"}{$identifier}},$beg;
			#}
			#else {
							push @{$rankProtMatch{$usedRankID}{$identifier}},$beg; # identifier is listed as often as matched by peptide
			#}
							push @{$begPosition{$identifier}},"$m:$#{$rankProtMatch{$usedRankID}{$identifier}}";

							#$matchQuality{$identifier}{$rank}++; # records the type and number of matching ranks (to exclude protein if matched by 1 maxRank)
							$maxProtMatch{$identifier}{$seqVmod}++; # match from different queries with same seq are counted only once! (10/05/04)
							$matchedSeq{$identifier}{$seqVmod}=1; #$multi; # records list of peptide sequences matching identifier
							$protDbRank{$identifier}=($databankRanks[$m])? $databankRanks[$m] : 1;
						}

						##<Good or auto/file decoy matches
						if ($matchGood) {

							##<Peptide selection
							my $select;
							if (!defined($bestScore{$seqVmod})) { # $data[5]=peptide sequence
								$select=0;
								$bestScore{$seqVmod}=$data[8];
							}
							elsif ($bestScore{$seqVmod}>=$data[8]) { # There is already an interpretation with better score
								$select=-1;
							}
							else { # new best score
								$select=0;
								###<Setting previous interpretations to SEL=-1
								foreach my $rID (@{$seqQueryRank{$seqVmod}}) {
									my ($qNum,$rk)=split(/:/,$rID);
									if ($dbMixedDecoyFlag) {
										next unless $tempQueryInfo{$qNum}{$rk}=~/SEL=0/; # rk-1=index in array
										$tempQueryInfo{$qNum}{$rk}=~s/SEL=0,LSP=0/SEL=-1,LSP=1/;
									}
									else {
										next unless $queryInfo{$qNum}[$rk-1]=~/SEL=0/; # rk-1=index in array
										$queryInfo{$qNum}[$rk-1]=~s/SEL=0,LSP=0/SEL=-1,LSP=1/;
									}
									$numValid{$qNum}--;
								}
								$bestScore{$seqVmod}=$data[8];
							}

							#<Spectral count
							$spectralCount{$seqVmod}++;

							my $numGoodMatches=scalar (keys %numProt);
							my $lowerScFlag=($select==-1)? 1 : 0; # Flag for lower-scoring peptide
							if ($dbMixedDecoyFlag) { # use tempQueryInfo
								$tempQueryInfo{$queryNum}{$rankGood}="SEL=$select,LSP=$lowerScFlag,SRK=$rank,MIS=$data[1],CALC=$data[2],DELT=$data[3],SEQ=$data[5],SC=$data[8]$varModDataStrg,MATCH=$numGoodMatches,";
								$tempQueryInfo{$queryNum}{$rankGood}.=$percolatorMascotScores{$section}{$queryNum}{$rankGood} if $percolatorInMascot && $percolatorMascotScores{$section}{$queryNum}{$rankGood};
							}
							else { # use queryInfo directly
								push @{$queryInfo{$queryNum}},"SEL=$select,LSP=$lowerScFlag,MIS=$data[1],CALC=$data[2],DELT=$data[3],SEQ=$sFlag$data[5],SC=$data[8]$varModDataStrg,MATCH=$numGoodMatches,";
								$queryInfo{$queryNum}[-1].=$percolatorMascotScores{$section}{$queryNum}{$rankGood} if $percolatorInMascot && $percolatorMascotScores{$section}{$queryNum}{$rankGood};
								###$queryInfo{$queryNum}[0].="PEP=$qvalityData{$data[8]}{PEP},QVAL=$qvalityData{$data[8]}{QVALUE}," if ($rank==1 && $qvalityData{$data[8]}); # adding PEP & q-value
							}
							$maxQueryScore{$queryNum}=$data[8] if (!defined($maxQueryScore{$queryNum}) || $maxQueryScore{$queryNum}<$data[8]); # best score for query (when rank=1?)
							$numValid{$queryNum}++ unless $select==-1;
							push @{$seqQueryRank{$seqVmod}},"$queryNum:$rankGood"; # $rankID
						}
						#elsif ($rank==1) { # initialize
						#	@{$queryInfo{$queryNum}}=() unless $dbMixedDecoyFlag;
						#	$maxQueryScore{$queryNum}=0;
						#}

						##<Decoy matches from mixed databank
						if ($matchBad) { # only in mixed decoy DB context ($dbMixedDecoyFlag)

							##<Peptide selection
							my $select;
							if (!defined($bestScoreBad{$seqVmod})) { # $data[5]=peptide sequence
								$select=0;
								$bestScoreBad{$seqVmod}=$data[8];
							}
							elsif ($bestScoreBad{$seqVmod}>=$data[8]) { # There is already an interpretation with better score
								$select=-1;
							}
							else { # new best score
								$select=0;
								###<Setting previous interpretations to SEL=-1
								foreach my $rID (@{$seqQueryRankBad{$seqVmod}}) {
									my ($qNum,$rk)=split(/:/,$rID);
									next unless $tempQueryInfo{$qNum}{$rk}=~/SEL=0/; # rk-1=index in array
									$tempQueryInfo{$qNum}{$rk}=~s/SEL=0,LSP=0/SEL=-1,LSP=1/;
									$numValid{$qNum}--;
								}
								$bestScoreBad{$seqVmod}=$data[8];
							}
							my $numDecoyMatches=scalar (keys %numFlaggedProt);
							my $lowerScFlag=($select==-1)? 1 : 0; # Flag for lower-scoring peptide
							$tempQueryInfo{"-$queryNum"}{$rankBad}="SEL=$select,LSP=$lowerScFlag,SRK=$rank,MIS=$data[1],CALC=$data[2],DELT=$data[3],SEQ=_$data[5],SC=$data[8]$varModDataStrg,MATCH=$numDecoyMatches,";
							$maxQueryScore{"-$queryNum"}=$data[8] if (!defined($maxQueryScore{"-$queryNum"}) || $maxQueryScore{"-$queryNum"}<$data[8]); # best score for query
							$numValid{"-$queryNum"}++ unless $select==-1;
							push @{$seqQueryRankBad{$seqVmod}},"-$queryNum:$rankBad"; # $rankID

							$tempQueryInfo{"-$queryNum"}{$rankBad}.= $percolatorMascotScores{$section}{abs($queryNum)}{$rankBad} if $percolatorInMascot && $percolatorMascotScores{$section}{abs($queryNum)}{$rankBad};
#if ($queryNum==6) {
#	print "<BR>RK=$rank, BAD ($rankBad): $tempQueryInfo{-$queryNum}{$rankBad}<BR>\n";
#}
						}
						#elsif ($rank==1) {
						#	#@{$queryInfo{"-$queryNum"}}=(); # initialize
						#	%{$tempQueryInfo{"-$queryNum"}}=();
						#	$maxQueryScore{"-$queryNum"}=0;
						#}
					}
					#else { # qQ_pR=-1 only occurs with p1
					#	$maxQueryScore{$queryNum}=0;
					#	if ($dbMixedDecoyFlag) {
					#		#@{$queryInfo{"-$queryNum"}}=();
					#		%{$tempQueryInfo{$queryNum}}=();
					#		%{$tempQueryInfo{"-$queryNum"}}=(); # is it useful?
					#		$maxQueryScore{"-$queryNum"}=0;
					#	}
					#	else {
					#		@{$queryInfo{$queryNum}}=();
					#	}
					#}
				}
				elsif ($line=~/^q\d+_p\d+_terms=(\S+)/) { # update match info to include boundary AA
					my @boundAA=split(/:/,$1);
					foreach my $identifier (keys %begPosition) {
						my $usedRankID=($dbMixedDecoyFlag && $identifier=~/^DECOY_/)? "-$queryNum:$rankBad" : "$queryNum:$rankGood"; # eg. DECOY_REV_Q9Z0P5 ('DECOY_' already added!)
						foreach my $posSet (@{$begPosition{$identifier}}) {
							my ($allPos,$goodPos)=split(/:/,$posSet);
							$rankProtMatch{$usedRankID}{$identifier}[$goodPos].=",$boundAA[$allPos]";
						}
					}
				}
				elsif ($line=~/^q\d+_p\d+_percolator=(\S+)/) { # percolator from PD
					my ($qValue,$PEP)=split(/:/,$1);
					if ($dbMixedDecoyFlag) { # use tempQueryInfo
						$tempQueryInfo{$queryNum}{$rankGood}.="PEP=$PEP,QVAL=$qValue," if $matchGood;
						$tempQueryInfo{$queryNum}{$rankBad}.="PEP=$PEP,QVAL=$qValue," if $matchBad;
					}
					else { # use queryInfo directly
						$queryInfo{$queryNum}[-1].="PEP=$PEP,QVAL=$qValue," if $matchGood;
					}
				}
				elsif ($line=~/^q\d+_p\d+_subst=(\S+)/) { # keep substitution information for Mascot
					my @substitutions =split(/,/,$1);
					my $subString='';
					for (my $i = 0 ; $i <$#substitutions ; $i=$i+3) {
						my ($posRes,$res1,$res2)=($substitutions[$i],$substitutions[$i+1],$substitutions[$i+2]);
						$subString.=" + $res1->$res2 ($posRes)";
					}
					if ($dbMixedDecoyFlag) {
						$tempQueryInfo{$queryNum}{$rankGood}.="SUBST=$subString," if $matchGood;
						$tempQueryInfo{$queryNum}{$rankBad}.="SUBST=$subString," if $matchBad;
					}else{
						$queryInfo{$queryNum}[-1].="SUBST=$subString," if $matchGood;
					}
				}
			}
		}

		###<Query sections: Elution information
		elsif ($section eq 'query' && !$decoyFile) {
			if ($line =~ /^title=(.+)/) {
				my $titleStrg=$1;
				while ($titleStrg=~/%(\w{2})/) {
					my $code=$1;
					my $char=chr(hex($code));
					$titleStrg=~s/%$code/$char/g;
				}
				if ($titleStrg=~/Elution.*: (.+)[Pp]eriod/) {
					my $elutStrg=$1;
					if ($elutStrg=~/^(\S+) to (\S+)/) {$elutionTime{$currentQuery}="$1 to $2";}
					elsif ($elutStrg=~/^(\S+) min/) {$elutionTime{$currentQuery}="$1";}
				}
				elsif ($titleStrg=~/Spectrum(\d+) scans:\s*(\d+),/) { # Orbitrap in case no rtinseconds / $titleStrg='Spectrum887 scans:2142,' or $titleStrg='Spectrum9853 scans: 15264,'
					$elutionTime{$currentQuery}="sp$1;sc$2;";
				}
				elsif ($titleStrg=~/[Ss]can/) {
					#if ($titleStrg=~/Scan (\d+)/) {$elutionTime{$currentQuery}="sc$1;";} # Match only 'Scan 5684'
					if ($titleStrg=~/[Ss]can[= ](\d+)/) {$elutionTime{$currentQuery}="sc$1;";} # Match 'scan=5684' and 'Scan 5684'
					elsif ($titleStrg=~/ scans in range (\d+) \([^\)]+\) to (\d+)/) {$elutionTime{$currentQuery}="sc$1-$2;";}
				}
			}
			elsif ($line=~/^rtinseconds=(\d+)/) { # Orbitrap
				$elutionTime{$currentQuery}.=sprintf "et%.2f;",$1/60; # unless $elutionTime{$currentQuery};
			}#TODO: ADD SCAN/Parameter
			elsif ($line=~/^charge=(\d+)\+/) {
				$charge{$currentQuery}=$1;
			}
			elsif ($line=~/^scans=([\d\-]+)/) { # Orbitrap velos?
				my $sc=$1;
				$elutionTime{$currentQuery}=~s/sc[\d\-]+;// if $elutionTime{$currentQuery}; # in case seen in title
				$elutionTime{$currentQuery}.="sc$sc;";
			}

			if ($dbMixedDecoyFlag) { # for decoy
				$elutionTime{"-$currentQuery"}=$elutionTime{$currentQuery};
				$charge{"-$currentQuery"}=$charge{$currentQuery};
			}
		}
		if ($nbLine > 1000000) {
			$nbLine=0;
			print '.';
		}

	}
	close DATAFILE;

	####<Build %queryInfo (stored rank can be shifted compare to Mascot rank)
	if ($dbMixedDecoyFlag) {
		foreach my $qNum (keys %tempQueryInfo) {
			@{$queryInfo{$qNum}}=();
			foreach my $rk (sort{$a<=>$b} keys %{$tempQueryInfo{$qNum}}) {
				push @{$queryInfo{$qNum}},$tempQueryInfo{$qNum}{$rk};
				$matchedQueryMIS{$qNum}=1;  # needed in PMF section if SQ search
			}
		}
	}

	####<Calculating protein max score
	print '.';
	foreach my $identifier (keys %matchedSeq) {
		my $refBestScore=($dbMixedDecoyFlag && $identifier=~/DECOY_/)? \%bestScoreBad : \%bestScore;
		foreach my $pepSeq (keys %{$matchedSeq{$identifier}}) {
			$maxProtScore{$identifier}+=$refBestScore->{$pepSeq};
		}
	}
	####<Setting protein sel_status
#	print '.';
#	foreach my $identifier (keys %matchQuality) {
# 		if (scalar (keys %{$matchQuality{$identifier}})==1 && $matchQuality{$identifier}{$maxRank} && $matchQuality{$identifier}{$maxRank}==1) {
# 			$protSelStatus{$identifier}=-2; # Protein is matched by only 1 maxRank => Exclude
# 		}
# 		else {
#			$protSelStatus{$identifier}=-1;
# 		}
#	}

	####<Adding spectral count information
	foreach my $queryNum (keys %queryInfo) {
		#next if $queryNum < 0; # not for decoy !!! COMMENTED PP 16/05/14: also for decoy to better compute FDR when lower-scoring
		foreach my $rkIdx (0..$#{$queryInfo{$queryNum}}) { # seq=>rank
			my ($seqVmod)=($queryInfo{$queryNum}[$rkIdx]=~ /SEQ=(\w+)/);
			next unless $seqVmod;
			if ($queryInfo{$queryNum}[$rkIdx]=~ /VMOD=([^,]+)/) {
				$seqVmod.=$1;
			}
			#if ($spectralCount{$seqData}) { # not for decoy
				$queryInfo{$queryNum}[$rkIdx].="SPC=$spectralCount{$seqVmod}," if $spectralCount{$seqVmod};
			#}
		}
	}

	####<Creating Match Groups
	#&createMatchGroups unless $mixedSearch; # will be done in PMF section

	print " Done.</FONT><BR>\n";
}


###################################################################
####<Parsing Mascot .xml file if creating a new MS/MS Analysis>####
###################################################################
sub parseMascotXML_MIS { # scanDB=1 (action=add || action=edit)
	my $msFileName = $_[0] ;
	print "<FONT class=\"title3\">&nbsp;-Extracting data from Mascot XML $msFileName...";

	my $parser = XML::SAX::ParserFactory->parser(Handler => MascotXMLSAXHandler->new($msFileName)  );
	$parser->parse_uri("$promsPath{valid}/ana_$analysisID/$msFileName"); #all is done in this line !
	print '...';


	####<Creating Match Groups
	#&createMatchGroups;
	print " Done.</FONT><BR>\n";
}


################################################################
####<Parsing Phenyx .xml file if creating a new MS Analysis>####
################################################################
sub parsePhenyxXML {
	####<Parsing data from $msFileName into globals %elutionTime,  %massExp, %massObs,%maxQueryScore, %matchList,%maxProtMatch,%numValid,%maxProtScore,  %queryInfo,
	# initialise "%massExp"   massexp{numero du query}
	# initialise "%massObs"   massObs{numero du query}
	# initialise "%elutionTime"    elutionTime{numero du query}
	# initialise  et utilise  "%maxQueryScore"   maxQueryScore{ numero du query}
	# initialise "%matchList" :    matchList{  n° de gi }   {  sequence + modif  } =1;
	# altere "%maxProtMatch"   maxProtMatch{  n° de gi } {sequence + modif }   = ( nombre de match pour un( n° de gi) de (sequence +modif))
	# affecte "%protSelStatus"  protSelStatus{(n° de gi)} = -1
	# somme "%maxProtScore"   maxProtScore{(n° de gi)} = (somme des score des peptide retenus) ;
	# altere "%numValid" .  numValid[numero de query}  nombre d'interpretation valid du query
	# utilise et altere ajoute "%queryInfo"  #   $queryInfo{$numero du query} [rang]  = datastring =  SEL=(0 si best score, -1 si non best score) , MIS=(number of missed clivage),  CALC=(mass peptide calcule), DELT=(delta mass), SEQ=(peptide sequence),SC=(score),VMOD=(varmodstring),MATCH=(nombre de prot qui matchent)
	# utilise et altere et concatene "%rankProtMatch"     rankProtMatch{ (n° de query):(n° de rang) }  {  n° de gi   }   [     ]    = position de debut du peptide dans la prot ;
		#sub :     %matchGroup,
		# utilise et detruit "%matchList"
		# utilise "%maxProtScore"
		# initialise et utilise "%matchGroup"

	my ($msFileName,$version) = @_;
	print "<FONT class=\"title3\">&nbsp;-Extracting data from Phenyx XML (v. $version) $msFileName...";
	#my %modiftranslation = ('Oxidation_M' =>'Oxidation (M)', 'ACET_nterm'=>'Acetyl (N-term)') ;

	my $parser = XML::SAX::ParserFactory->parser(Handler => PhenyxSAXHandler->new($msFileName)  );
	$parser->parse_uri("$promsPath{valid}/ana_$analysisID/$msFileName"); #all is done in this line !
	print '...';

	####<Creating Match Groups
	#&createMatchGroups;


	#######parsing taxonomy_codes.txt to obtain taxonomy ########
	if ($phenyxTaxoNumber) {
# 		print " \n Parsing taxonomy file.<BR>\n";
		my $TaxoFile ="$promsPath{'banks'}/taxonomy_codes.txt";
		open (TAXOFILE, $TaxoFile);
		while (<TAXOFILE>) {
			my $line = $_ ;
			$line =~/^(\d+)\t/;
			if ($1==$phenyxTaxoNumber) {
				($phenyxTaxo)=($line =~/^\d+\t([^\t]*)/);
				last ;
			}
		}
		close TAXOFILE ;
		$phenyxTaxo = "All entries" if $phenyxTaxo eq "root" ;
	}

# 	print " \n parsing Phenyx XML Done.</FONT><BR><BR>\n";
	print " Done.</FONT><BR>\n";
}


#############################################################################
###<Parsing XML file got from PARAGON searches on ProteinPilot (AB Sciex)>###
### Group files converted with the executable group2xml                   ###
#############################################################################
sub parseParagonXML {
	my ($msFileName,$allVmods) = @_;
	#my (%bestScore,%matchedSeq); # locals required for calculation of protein max score

	print "<FONT class=\"title3\">&nbsp;-Extracting data from Paragon XML (from ProteinPilot) $msFileName...";
	my $handler = ParagonXMLHandler->new("$promsPath{valid}/ana_$analysisID/$msFileName",$allVmods);
	my $xmlparser = XML::SAX::ParserFactory->parser(Handler => $handler );
	$xmlparser->parse_uri("$promsPath{valid}/ana_$analysisID/$msFileName");

	print '...';

	print " Done.</FONT><BR>\n";
}


#######################################
####<Parsing .tandem.pep.xml file >####
#######################################
sub parseTandemPEPXML {
	my ($msFileName,$protDBSeq) = @_;
	print "<FONT class=\"title3\">&nbsp;-Extracting data from Tandem XML (from X! Tandem) $msFileName...";
	my $parser = XML::SAX::ParserFactory->parser(Handler => XTandemPEPXMLHandler->new($msFileName,$protDBSeq,$dbh,$databankIDs[0],$minScore,$maxRank,$analysisID,\%protSeq,\%queryInfo,\%protDes,\%protOrg,\%protLength,\%charge,\%massObs,\%massExp,\%elutionTime,\%numValid,\%matchList,\%maxProtMatch,\%maxProtScore,\%maxQueryScore,\%analysisMods,\%rankProtMatch,\%protMW,\%protDbRank)  );
	$parser->parse_uri("$promsPath{valid}/ana_$analysisID/$msFileName");

	print '...';
	print " Done.</FONT><BR>\n";
}


###############################
####<Parsing .tandem file >####
###############################
sub parseTandemXML {
	my ($msFileName,$modifNTer) = @_;
	print "<FONT class=\"title3\">&nbsp;-Extracting data from Tandem XML (from X! Tandem) $msFileName...";
	my $parser = XML::SAX::ParserFactory->parser(Handler => XTandemXMLHandler->new($msFileName,$modifNTer,$dbh,$databankIDs[0],$minScore,$analysisID,$maxRank,\%protSeq,\%queryInfo,\%protDes,\%protOrg,\%protLength,\%charge,\%massObs,\%massExp,\%elutionTime,\%numValid,\%matchList,\%maxProtMatch,\%maxProtScore,\%maxQueryScore,\%analysisMods,\%rankProtMatch,\%protMW,\%protDbRank)  );
	$parser->parse_uri("$promsPath{valid}/ana_$analysisID/$msFileName");

	print '...';
	print " Done.</FONT><BR>\n";
}

##########################################
####<Generates Temporary Match Groups>####
##########################################
sub createMatchGroups { # Globals %matchList, %maxProtScore, %matchGroup
	print "<FONT class=\"title3\">&nbsp;-Grouping proteins...";
	my $numGroup=0;
	my @sortedIdentifiers=(sort{scalar (keys %{$matchList{$b}})<=>scalar (keys %{$matchList{$a}}) || $maxProtScore{$b}<=>$maxProtScore{$a} || $a cmp $b} keys %matchList);
	my $count=0;
	foreach my $i (0..$#sortedIdentifiers) {
		$count++;
		if ($count==150) {print '.'; $count=0;}
		next if $matchGroup{$sortedIdentifiers[$i]}; # already assigned to a match group
		$matchGroup{$sortedIdentifiers[$i]}=++$numGroup;
#next; # SKIP grouping!!!
		foreach my $j ($i+1..$#sortedIdentifiers) {
			next if $matchGroup{$sortedIdentifiers[$j]}; # already assigned to a match group
			##<Comparing peptide contents of identifier#1 and identifier#2>## All peptides must match!
			my $matchOK=1;
			foreach my $seq (keys %{$matchList{$sortedIdentifiers[$j]}}) {
				if (!$matchList{$sortedIdentifiers[$i]}{$seq}) {
					delete $matchList{$sortedIdentifiers[$i]}{$seq}; # to be safe
					$matchOK=0;
					last;
				}
			}
			$matchGroup{$sortedIdentifiers[$j]}=$matchGroup{$sortedIdentifiers[$i]} if $matchOK;
		}
	}
	print " Done.</FONT><BR>\n";
}


################################
####<Mascot boundary String>####
################################
sub getBoundaryString {
    return '--gc0p4Jq0M2Yt08jU534c0p';
}

############################################################
####<Draw score distribution for precomputed percolator>#### Percolator in Proteome Discoverer
############################################################
sub drawScoreDistribution {
	# Globals: %queryInfo, %maxQueryScore, %rankProtMatch, %maxProtMatch, %matchList
	my ($fdr,$anaID,$dataFile)=@_;
	print "<FONT class=\"title3\">&nbsp;-Drawing score distribution for $fdr% FDR...";

	(my $rootName=$dataFile)=~s/\.[^\.]+\Z//;
	my $targetFile="$promsPath{valid}/ana_$anaID/$rootName.target";
	my $decoyFile="$promsPath{valid}/ana_$anaID/$rootName.decoy";
	my $qvalityFile="$promsPath{valid}/ana_$anaID/$rootName.qvality";

	###<Generating score files>###
	open (TARGET,">$targetFile");
	open (DECOY,">$decoyFile");
	my ($numTarget,$numDecoy)=(0,0);
	foreach my $queryNum (sort{$maxQueryScore{$b}<=>$maxQueryScore{$a}} keys %maxQueryScore) { # sort just for helping debug
		if ($queryNum > 0) {
			print TARGET "$maxQueryScore{$queryNum}\n";
			$numTarget++;
		}
		else {
			print DECOY "$maxQueryScore{$queryNum}\n";
			$numDecoy++;
		}
	}
	close TARGET;
	close DECOY;
	print '.';

	##<Drawing density plot with R
	my $Rscript="$promsPath{valid}/ana_$anaID/density.R";
	open(R,">$Rscript");
	print R qq
|target <- t(read.table("$targetFile"))
decoy <- t(read.table("$decoyFile"))
densT <- density(target)
densD <- density(decoy)

minX <- min(c(min(densT\$x,densD\$x)))
maxX <- max(c(max(densT\$x,densD\$x)))
maxY <- max(c(max(densT\$y*$numTarget,densD\$y*$numDecoy)))

png(filename="$promsPath{valid}/ana_$anaID/scores.png", width=500, height=500, units = "px")

plot(densT\$x,densT\$y*$numTarget, xlim=c(minX,maxX), ylim=c(0,maxY*1.1), yaxs='i', type='l', col='green', xlab="Scores", ylab="Abundance", main="$rootName ($fdr% FDR precomputed)")
lines(densD\$x,densD\$y*$numDecoy, col='red')
legend("topright", c("Peptides","Decoy peptides"), cex=0.8, col=c('green','red'), lwd=2, bty="n")

dev.off()
|;
	close R;
	system "cd $promsPath{valid}/ana_$anaID; $promsPath{R}/R CMD BATCH --no-save --no-restore $Rscript";
	sleep 1;
	unlink $Rscript;
	unlink $Rscript.'out' if -e $Rscript.'out';
	print " Done.</FONT><BR>\n";
}


#####################################################
####<Computing min score with qvality or DTcount>####
#####################################################
sub applyFDR {
	# Globals: %queryInfo, %maxQueryScore, %rankProtMatch, %maxProtMatch, %matchList
	my ($FDRalgo,$fdr,$anaID,$dataFile,$userScore)=@_;
	print "<FONT class=\"title3\">&nbsp;-Estimating threshold score for $fdr% FDR...";
	#my $rootName;
	(my $rootName=$dataFile)=~s/\.[^\.]+\Z//;

	my $targetFile="$promsPath{valid}/ana_$anaID/$rootName.target";
	my $decoyFile="$promsPath{valid}/ana_$anaID/$rootName.decoy";
	my $qvalityFile="$promsPath{valid}/ana_$anaID/$rootName.qvality";
	###<Generating score files>###
	open (TARGET,">$targetFile");
	open (DECOY,">$decoyFile");
	my $fdrScore=0;
	my $FDRsuccess=0;
	my ($numTarget,$numDecoy)=(0,0);
	if ($FDRalgo eq 'qvality') {
		###<Generating score files>###
		foreach my $queryNum (sort{abs($a)<=>abs($b)} keys %maxQueryScore) { # sort just for helping debug
			if ($queryNum > 0) {
				print TARGET "$maxQueryScore{$queryNum}\n";
				$numTarget++;
			}
			else {
				print DECOY "$maxQueryScore{$queryNum}\n";
				$numDecoy++;
			}
		}
	}
	elsif($FDRalgo eq 'mayu') {
		###<Generating score files>###
		foreach my $queryNum (sort{abs($a)<=>abs($b)} keys %maxQueryScore) { # sort just for helping debug
			if ($queryNum > 0) {
				print TARGET "$maxQueryScore{$queryNum}\n";
				$numTarget++;
			}
			else {
				print DECOY "$maxQueryScore{$queryNum}\n";
				$numDecoy++;
			}
		}

		###>processing Mayu
		(my $ext=$dataFile)=~s/\w+[^\.]//;
		($rootName=$dataFile)=~s/\..*//;

		my $anaDir="$promsPath{valid}/ana_$anaID";

		my $pepXMLFile;
		my ($DBFile,$decoyTag)=$dbh->selectrow_array("SELECT FASTA_FILE,DECOY_TAG FROM DATABANK WHERE ID_DATABANK=$databankIDs[0]");
		my $DBPath="$promsPath{data}/banks/db_$databankIDs[0]/$DBFile";
		my $mzxmlPath="$promsPath{data}/tmp/upload/project_$projectID";

		###> data files conversion into .pep.xml
		if ($ext eq '.xml' || $ext eq '.tandem') {
			if ($ext eq '.tandem') {
                system "cd $anaDir; mv $dataFile $rootName\_Mod.tandem;";
				open(XMLFILE,"<","$anaDir/$rootName\_Mod.tandem") or die ("open : $!");
            }
            else{open(XMLFILE,"<","$anaDir/$dataFile") or die ("open : $!");}

			open(OUT,">","$anaDir/$rootName.tandem") or die ("open : $!");
			while (<XMLFILE>) {
				if ($_=~/protein, cleavage N-terminal mass change/) {
					(my $line=$_)=~s/protein, cleavage N-terminal mass change">.*<\/note>/protein, cleavage N-terminal mass change"><\/note>/;
					print OUT $line;
				} elsif($_=~/protein, cleavage C-terminal mass change/){
					(my $line=$_)=~s/protein, cleavage C-terminal mass change">.*<\/note>/protein, cleavage C-terminal mass change"><\/note>/;
					print OUT $line;
				} elsif($_=~/spectrum, path/){
					#(my $line=$_)=~s/spectrum, path\">.*$rootName\.mzXML/spectrum, path\">$anaDir\/$rootName\.mzXML/;
					(my $line=$_)=~s/spectrum, path\">.*$rootName\.mzXML/spectrum, path\">$mzxmlPath\/$rootName\.mzXML/;
					print OUT $line;
				} else{ print OUT $_; }
			}
			close XMLFILE;
			close OUT;
			system "cd $anaDir;$promsPath{tpp}/Tandem2XML $rootName.tandem $rootName.tandem.pep.xml";
			$pepXMLFile="$rootName.tandem.pep.xml";
		}
		elsif($ext eq '.dat'){
			symlink("$mzxmlPath/$rootName.mzXML","$anaDir/$rootName.mzXML");   ## mzXML file is needed for .dat conversion
			system "cd $anaDir;$promsPath{tpp}/Mascot2XML $dataFile -D$DBPath -Etrypsin -notgz >>output.txt 2>&1 &";
			my $i=0;
			while (! -e "$anaDir/$rootName.pep.xml") {
				sleep 2;
				$i++;
				print '.' if $i==10;
				$i=0 if $i==10;
			}
			$pepXMLFile="$rootName.pep.xml";
		}
		elsif($ext eq '.tandem.pep.xml'){
			$pepXMLFile="$rootName.tandem.pep.xml";
		}

        ###> Modification of DB path into PepXML files
        unless ($ext eq '.dat'){
			my $qDBPath=quotemeta($DBPath);
			my $qMzxmlPath=quotemeta($mzxmlPath);
			system "sed -r -i -e 's/[^\"]+\\\/$rootName/$qMzxmlPath\\\/$rootName/g' $anaDir/$pepXMLFile" if $ext eq '.tandem.pep.xml';
			system "sed -r -i -e 's/local_path=\".*\.fasta\"/local_path=\"$qDBPath\"/g' $anaDir/$pepXMLFile";
            system "sed -r -i -e 's/<parameter name=\"list path, sequence source #1\" value=\".*\.fasta\"/<parameter name=\"list path, sequence source #1\" value=\"$qDBPath\"/g' $anaDir/$pepXMLFile";
        }

		###> xinteract
		system "cd $anaDir; $promsPath{tpp}/xinteract -OARPd -d$decoyTag -Ninteract.pep.xml $pepXMLFile >output.txt 2>&1 &";
		my $i=0;
        while (`tail -n 1 $anaDir/output.txt` !~ /job completed in \d* sec/) {
			if (`tail -n 1 $anaDir/output.txt` =~ /QUIT - the job is incomplete/) {
                print "<BR><BR><FONT color=\"red\">*ERROR* Xinteract don't run ! </FONT><BR>";
				exit;
            }
			sleep 2;
			$i++;
			print '.' if $i==10;
			$i=0 if $i==10;
		}


		###>InterProphet
		system "cd $anaDir; $promsPath{tpp}/InterProphetParser DECOY=$decoyTag interact.pep.xml iProphet.pep.xml >output.txt 2>&1 &";
		$i=0;
		while (!-e "$anaDir/iProphet.pep.xml") {
			sleep 2;
			$i++;
			print '.' if $i==10;
			$i=0 if $i==10;
		}

		###>Mayu
		system "cd $anaDir; $promsPath{tpp}/Mayu.pl -A iProphet.pep.xml -C $DBPath -E $decoyTag -G 0.01 -H 51 -I 2 -P protFDR=0.01:t >output.txt 2>&1";



		##>Get modifications ID and names from table MODIFICATION
		my %modificationList;
		my $modificationSth=$dbh->prepare("SELECT PSI_MS_NAME,ID_MODIFICATION,MONO_MASS FROM MODIFICATION");
		$modificationSth->execute;
		while (my ($modName,$modID,$monoMass)=$modificationSth->fetchrow_array) {
			$modificationList{$modID}="$monoMass&$modName" if ($monoMass && $modName);
		}
		#$dbh->disconnect;

		##>Recovering valid spectrum scans, peptides and proteins from Mayu result's file
		my (%validScan,%validProt,%validPepProt);
		my $mayuFile=`ls $anaDir/*psm_*`;
		chomp $mayuFile;
		if (-e $mayuFile) {
			open(MAYUOUT,"<","$mayuFile") or die ("open : $!");
			while (<MAYUOUT>) {
				next if $.==1;
				chomp $_;
				my ($spectrum,$pep,$prot,$mod,$score,$decoy,$mFDRsplit)=split(/,/,$_);
				$fdrScore=$score if (!$fdrScore || $fdrScore>=$score);

				## Recovering spectrum number
				my ($fileName,$scan,$id,$charge)=split(/\./,$spectrum);
				($scan=$scan)=~s/^0//;

				## Recovering modifications and peptides sequence
				my @sequence=split(//,$pep);
				my %modification;
				foreach my $modification (split(/:/,$mod)){
					my ($posMod,$mass)=split(/=/,$modification);
					my $varMod;
					my %massAAave=&promsConfig::getMassAAmono;
					my $aaMass=$massAAave{$sequence[$posMod-1]};
					my $modifMass=$mass-$aaMass;
					foreach my $modID (keys %modificationList){
						my ($mass,$modName)=split(/&/,$modificationList{$modID});
						if ($mass > $modifMass-0.005 && $mass < $modifMass+0.005) {
							$varMod=$modName;
							last;
						}
					}
					$modification{$varMod}{$sequence[$posMod-1]}{$posMod}=1;
				}
				my $varModString;
				foreach my $varMod (sort {$a cmp $b} keys %modification){
					foreach my $aaMod (keys %{$modification{$varMod}}){
						my $posModif;
						foreach my $posMod (sort {$a <=> $b} keys %{$modification{$varMod}{$aaMod}}){
							$posModif=($posModif)? $posModif.".".$posMod : $posMod;
						}
						$varModString.=" + $varMod ($aaMod:$posModif)";
					}
				}
				my $peptide="$pep";
				$peptide.="$varModString" if $varModString;
				$validProt{$prot}{$peptide}=1;
				$validPepProt{$peptide}{$prot}=1;
				$validScan{$scan}{$peptide}=1;
			}
			close MAYUOUT;
			print "] Done (Score=$fdrScore).</FONT><BR>\n";

			##>Updating data structures
			print "<FONT class=\"title3\">&nbsp;-Filtering data based on threshold score...";
			foreach my $queryNum (keys %maxQueryScore) {
				my ($scan)=($elutionTime{$queryNum}=~/sc(\d*);/);
				if($validScan{$scan}){				## valid spectrum
					for (my $i=0; $i<=$#{$queryInfo{$queryNum}}; $i++) {
						my ($sel)=($queryInfo{$queryNum}[$i]=~/SEL=([^,]+)/);
						my ($pepSeq)=($queryInfo{$queryNum}[$i]=~/SEQ=([^,]+)/);
						$pepSeq=~s/^_// if $pepSeq=~/^_/;
						my ($vMod)=($queryInfo{$queryNum}[$i]=~/VMOD=([^,]+)/);
						my ($sc)=($queryInfo{$queryNum}[$i]=~/SC=([^,]+)/);
						$pepSeq.=$vMod if $vMod ;
						(my $pepSeqMod=$pepSeq)=~s/ \+ Acetyl \(N-term\)// if $vMod && $vMod=~/N-term/;

						my $match=0;
						foreach my $peptide (keys %{$validScan{$scan}}){
							$match=($pepSeq eq $peptide) ? 1 : ($pepSeqMod && $pepSeqMod eq $peptide) ? 1 : 0;
							last if $match;
						}
						my $rankID="$queryNum:".($i+1);
						foreach my $identifier (keys %{$rankProtMatch{$rankID}}) {
							my $deletedPep=($pepSeqMod) ? $pepSeqMod : $pepSeq;
							if ($validPepProt{$deletedPep}{$identifier}){
								delete $rankProtMatch{$rankID}{$identifier} unless $match;
							}
							else{
								if ($sel==0) { # selectable interpretation
									delete $maxProtMatch{$identifier}{$pepSeq};
									delete $matchList{$identifier}{$pepSeq};
									$maxProtScore{$identifier}=sprintf('%.2f',$maxProtScore{$identifier}-$sc);
								}
								delete $rankProtMatch{$rankID}{$identifier};
								delete $maxProtMatch{$identifier} unless scalar keys %{$maxProtMatch{$identifier}};
							}
						}
						delete $rankProtMatch{$rankID} if $queryNum=~/-/;

						if (!$match || $queryNum=~/-/) {
							delete $queryInfo{$queryNum}[$i];
							$numValid{$queryNum}-- if $sel==0;
						}
						delete $rankProtMatch{$rankID} unless scalar keys %{$rankProtMatch{$rankID}};
					}
					@{$queryInfo{$queryNum}}=() unless $queryInfo{$queryNum};
				}
				else{
					for (my $i=0; $i<=$#{$queryInfo{$queryNum}}; $i++) {
						my ($sel)=($queryInfo{$queryNum}[$i]=~/SEL=([^,]+)/);
						my ($sc)=($queryInfo{$queryNum}[$i]=~/SC=([^,]+)/);
						my ($pepSeq)=($queryInfo{$queryNum}[$i]=~/SEQ=([^,]+)/);
						$pepSeq=~s/^_// if $pepSeq=~/^_/;
						my ($vMod)=($queryInfo{$queryNum}[$i]=~/VMOD=([^,]+)/);
						$pepSeq.=$vMod if $vMod;
						(my $pepSeqMod=$pepSeq)=~s/ \+ Acetyl \(N-term\)// if $vMod && $vMod=~/N-term/;

						my $rankID="$queryNum:".($i+1);
						foreach my $identifier (keys %{$rankProtMatch{$rankID}}) {
							my $deletedPep=($pepSeqMod) ? $pepSeqMod : $pepSeq;
							#print "noscan $identifier,$deletedPep,$queryNum<BR>" if $identifier eq 'DECOY_reverse_sp|A0A075B759|PAL4E_HUMAN';
							next if $validPepProt{$deletedPep}{$identifier};
							if ($sel==0) { # selectable interpretation
								delete $maxProtMatch{$identifier}{$pepSeq};
								delete $matchList{$identifier}{$pepSeq};
								$maxProtScore{$identifier}=sprintf('%.2f',$maxProtScore{$identifier}-$sc);
							}
							delete $maxProtMatch{$identifier} unless scalar keys %{$maxProtMatch{$identifier}};
						}
						delete $rankProtMatch{$rankID};
					}
					@{$queryInfo{$queryNum}}=();
					delete $maxQueryScore{$queryNum} if $maxQueryScore{$queryNum};
					delete $elutionTime{$queryNum};
					delete $numValid{$queryNum} if $numValid{$queryNum};
				}
			}
			$FDRsuccess=1;

			##analysis.fasta
			if (%protSeq) {
				open(ANAFASTA,">>$promsPath{valid}/ana_$analysisID/analysis.fasta");
				foreach my $identifier (keys %protSeq){
					next unless $maxProtMatch{$identifier};
					print ANAFASTA ">",$identifier,"\n$protSeq{$identifier}\n";
				}
				close ANAFASTA;
			}

			###>Deleting intermediate files
			push @temporaryFiles,"$anaDir/output.txt";
			push @temporaryFiles,$mayuFile;
			push @temporaryFiles,"$anaDir/$rootName.mzXML";
			push @temporaryFiles,"$anaDir/$rootName\_Mod.tandem" if -e "$anaDir/$rootName\_Mod.tandem";
			push @temporaryFiles,"$anaDir/$rootName.tandem" if -e "$anaDir/$rootName.tandem" && $ext ne '.tandem';
			push @temporaryFiles,"$anaDir/$rootName.tandem.pep.xml" unless $ext eq '.tandem.pep.xml';

			system "rm $anaDir/interact*";
			system "rm $anaDir/iProphet.pep.xml";
			system "rm $anaDir/*_main_*";

		}
		else{print "<BR><BR><FONT color=\"red\">*ERROR* Mayu don't run ! </FONT><BR>"; exit;}
	}
	else{# DTcount
		my $fracFDR=$fdr/100;
		my $upScore=0;
		my $match=0;
		my $nextScore=0;
#open (TEST,">$promsPath{valid}/ana_$anaID/test.txt");
#print TEST "PEP\tSC\tFDR\tUP\tNEXT\n";
		foreach my $queryNum (sort{$maxQueryScore{$b}<=>$maxQueryScore{$a}} keys %maxQueryScore) { # sort just for helping debug
			if ($queryNum > 0) {
				print TARGET "$maxQueryScore{$queryNum}\n";
				$numTarget++;
			}
			else {
				print DECOY "$maxQueryScore{$queryNum}\n";
				$numDecoy++;
			}
			#$fdrScore=$maxQueryScore{$queryNum} if $numDecoy/$numTarget <= $fracFDR;
			$upScore=$maxQueryScore{$queryNum} unless $upScore;
			if ($numDecoy/$numTarget <= $fracFDR) {
				$upScore=$fdrScore if ($fdrScore && $upScore > $fdrScore); # update upScore if new match with new score
				$fdrScore=$maxQueryScore{$queryNum};
				$match=1;
			}
			elsif ($match) { # next loop after match loop: decoy found (otherwise 'if' above is true)
				$nextScore=$maxQueryScore{$queryNum}; # record score
				$match=0;
			}
#print TEST "$pep\t$maxQueryScore{$queryNum}\t$fdrScore\t$upScore\t$nextScore\n";
		}
		$fdrScore=$upScore if ($upScore && $fdrScore==$nextScore); # Increase score if same score is found in following loops with additional decoys

	}
#close TEST;
	close TARGET;
	close DECOY;

	print '.';

	unless ($numDecoy) {
		print " <FONT color=\"#DD0000\">WARNING: No decoy data found!</FONT> (Score=$fdrScore).</FONT><BR>\n"; # 0!
		$FDRsuccess=1;
		return ($FDRsuccess,$fdrScore);
	}

	my $numCycles=0;
	my $endFile="$promsPath{valid}/ana_$anaID/fdr.end";
	my %qvalityData;
	###<Running Qvality>###
	if ($FDRalgo eq 'qvality') {

		system "$promsPath{qvality}/qvality -o \"$qvalityFile\" \"$targetFile\" \"$decoyFile\" 2> $promsPath{valid}/ana_$anaID/out.qvality&";
		#print qq|[<SCRIPT LANGUAGE="JavaScript">systemFrame.location="./storeAnalyses.cgi?ACT=qvality1&ANA_ID=$anaID&FILE=$rootName";</SCRIPT>|;
		#while ($numCycles < 60 && !-e $endFile) { # 2 min wait max
		while ($numCycles < 60 && !-e $qvalityFile) {
			sleep 2;
			$numCycles++;
			print '.';
		}
		sleep 2;
		#unlink $endFile if -e $endFile;
		print '.';

		###<Estimating threshold score>###
		my $maxQvalue=$fdr/100;

		if (-e $qvalityFile) {
			open (QVALITY,$qvalityFile);
			my ($prevSc,$prevQvalue);
			while (<QVALITY>) {
				next if $.==1; # skip header
				chomp;
				my ($sc,$pep,$qValue)=split(/\t/,$_);
				#last if $qValue > $maxQvalue;
				@{$qvalityData{$sc}}=($pep,$qValue); # inlcuding last qValue (>maxQvalue) in case rounded linear estimation = this value
				if ($qValue > $maxQvalue) { # linear estimation of fdr between flanking pairs of scores & qValues (PP 15/10/12)
					my $slope=($sc-$prevSc)/($qValue-$prevQvalue); # a in y=a*x+b (a=(y2-y1)/(x2-x1))
					my $zeroXvalue=$sc-$slope*$qValue; # b in y=a*x+b (b=y-a*x)
					$fdrScore=sprintf '%.5f',$slope*$maxQvalue+$zeroXvalue; # qVality rounds numbers to max 5 digits after .
					#print "$zeroXvalue,$slope<BR>";
					$fdrScore*=1; # 55.00000 -> 55
#print "<BR>**<BR>prevQvalue=$prevQvalue, qValue=$qValue<BR>prevSc=$prevSc, sc=$sc<BR>fdrScore=$fdrScore<BR>**<BR>\n";
					last;
				}
				#$fdrScore=$sc;
				$prevSc=$sc;
				$prevQvalue=$qValue;
			}
			close QVALITY;
			$fdrScore=$prevSc unless $fdrScore; # just to be safe
			print './';
		}
	}
	elsif($FDRalgo ne 'mayu') {print '[';}

	if ($FDRalgo ne 'qvality'  || -e $qvalityFile) {
		##<Drawing density plot with R
		my $Rscript="$promsPath{valid}/ana_$anaID/density.R";
		open(R,">$Rscript");
		print R qq
|target <- t(read.table("$targetFile"))
decoy <- t(read.table("$decoyFile"))
densT <- density(target)
densD <- density(decoy)

minX <- min(c(min(densT\$x,densD\$x)))
maxX <- max(c(max(densT\$x,densD\$x)))
maxY <- max(c(max(densT\$y*$numTarget,densD\$y*$numDecoy)))

png(filename="$promsPath{valid}/ana_$anaID/scores.png", width=500, height=500, units = "px")

plot(densT\$x,densT\$y*$numTarget, xlim=c(minX,maxX), ylim=c(0,maxY*1.1), yaxs='i', type='l', col='green', xlab="Scores", ylab="Abundance", main="$rootName ($fdr% FDR)")
lines(densD\$x,densD\$y*$numDecoy, col='red')
abline(v=$fdrScore, col='blue')
legend("topright", c("Peptides","Decoy peptides","Threshold score"), cex=0.8, col=c('green','red','blue'), lwd=2, bty="n")

dev.off()
|;
		close R;

		system "cd $promsPath{valid}/ana_$anaID; $promsPath{R}/R CMD BATCH --no-save --no-restore $Rscript&";
		#print qq|<SCRIPT LANGUAGE="JavaScript">systemFrame.location="./storeAnalyses.cgi?ACT=qvality2&ANA_ID=$anaID";</SCRIPT>|;

		$numCycles=0;
		#while ($numCycles < 30 && !-e $endFile) { # 1 min wait max
		while ($numCycles < 60 && !-e "$promsPath{valid}/ana_$anaID/scores.png") {
			sleep 2;
			$numCycles++;
			print '.';
		}
		sleep 2;

		unlink $Rscript;
		unlink $Rscript.'out' if -e $Rscript.'out';
		#unlink $endFile if -e $endFile;

		$FDRsuccess=1;
		print "] Done (Score=$fdrScore).</FONT><BR>\n" unless $FDRalgo eq 'mayu';
	}

	unless ($FDRsuccess) {
		$fdrScore=$userScore;
		print "] <FONT color=\"#DD0000\">ERROR!</FONT> Using user-defined score=$userScore</FONT><BR>\n";
	}
	if ($FDRalgo ne 'mayu') {

		###<Updating data structures>###
		print "<FONT class=\"title3\">&nbsp;-Filtering data based on threshold score...";
		##< %queryInfo & %rankProtMatch >##
		my %deletedProtMatches;
		foreach my $queryNum (keys %maxQueryScore) {
			RANK: for (my $i=0; $i<=$#{$queryInfo{$queryNum}}; $i++) {
				my ($rkScore0)=($queryInfo{$queryNum}[$i]=~/SC=([^,]+)/);
				my $rkScore=$rkScore0;
				if ($rkScore >= $fdrScore) {
					if ($FDRalgo eq 'qvality') { # Adding PEP & QVAL values
						$rkScore=sprintf "%.5f",$rkScore0; # qVality rounds numbers to max 5 digits after .
						$rkScore*=1;  # 55.00000 -> 55 (due to qvality output)
						if ($FDRsuccess && $queryNum>0 && $i==0) { # only for 1st rank
							unless ($qvalityData{$rkScore}) { # rounding to 5 digits failed: try 4
								$rkScore=sprintf "%.4f",$rkScore0;
								$rkScore*=1;  # 55.0000 -> 55 (due to qvality output)
							}
							$queryInfo{$queryNum}[$i].="PEP=$qvalityData{$rkScore}[0],QVAL=$qvalityData{$rkScore}[1],"; # if ($queryNum>0 && $i==0); # only for 1st rank
						}
					}
				}
				else { # reject interpretation
					for (my $j=$i; $j<=$#{$queryInfo{$queryNum}}; $j++) {
						my ($sel)=($queryInfo{$queryNum}[$j]=~/SEL=([^,]+)/);
						my $sc=0;
						my $seqVmod='';
						if ($sel==0) { # selectable interpretation
							$numValid{$queryNum}--;
							($sc)=($j==$i)? ($rkScore) : ($queryInfo{$queryNum}[$j]=~/SC=([^,]+)/);
							($seqVmod)=($queryInfo{$queryNum}[$j]=~/SEQ=([^,]+)/);
							my ($vMod)=($queryInfo{$queryNum}[$j]=~/VMOD=([^,]+)/);
							$seqVmod.=$vMod if $vMod;
						}
						$queryInfo{$queryNum}[$j]=undef; # clear interpretation

						my $rankID="$queryNum:".($j+1);
						foreach my $identifier (keys %{$rankProtMatch{$rankID}}) {
							push @{$deletedProtMatches{$identifier}},[$seqVmod,$sel,$sc];
							delete $rankProtMatch{$rankID}{$identifier};
						}
						delete $rankProtMatch{$rankID};
					}
					#$maxQueryScore{$queryNum}=0 if $i==0; # all ranks excluded
					delete $maxQueryScore{$queryNum} if $i==0; # all ranks excluded
					last RANK;
				}
			}
		}
		print '.';

		##< %maxProtMatch, %matchList >##
		foreach my $identifier (keys %deletedProtMatches) {
			foreach my $refRankData (@{$deletedProtMatches{$identifier}}) {
				my ($seqVmod,$sel,$sc)=@{$refRankData};
				if ($sel==0) { # selectable interpretation
					delete $maxProtMatch{$identifier}{$seqVmod};
					delete $matchList{$identifier}{$seqVmod}; # used by &createMatchGroups
					$maxProtScore{$identifier}=sprintf('%.1f',$maxProtScore{$identifier}-$sc);
				}
			}
			delete $maxProtMatch{$identifier} unless scalar keys %{$maxProtMatch{$identifier}};
		}
	}
	print " Done.</FONT><BR>\n";
	return ($FDRsuccess,$fdrScore);
}

sub systemCommand {

	print header(-'content-encoding'=>'no',-charset=>'utf-8'); # start_html,"\n"; # start_html required to force update
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>System Command</TITLE>
</HEAD>
<BODY>
|;
	###<Qvality-related commands>###
	if ($action=~/qvality(\d)/) {
		my $step=$1;
		my $anaID=param('ANA_ID');

		##<Running qvality>##
		if ($step==1) {
			my $rootName=param('FILE');
			my $targetFile="$promsPath{valid}/ana_$anaID/$rootName.target";
			my $decoyFile="$promsPath{valid}/ana_$anaID/$rootName.decoy";
			my $qvalityFile="$promsPath{valid}/ana_$anaID/$rootName.qvality";
			system "$promsPath{qvality}/qvality -o \"$qvalityFile\" \"$targetFile\" \"$decoyFile\" 2> $promsPath{valid}/ana_$anaID/out.qvality";
		}
		##<Generating score distribution image>##
		else {
			my $Rscript="$promsPath{valid}/ana_$anaID/density.R";
			system "cd $promsPath{valid}/ana_$anaID; $promsPath{R}/R CMD BATCH --no-save --no-restore $Rscript";
		}
		open(OUT,">$promsPath{valid}/ana_$anaID/fdr.end");
		print OUT "STEP $step";
		close OUT;
	}

	###<Launching background DB scan>###
	elsif ($action eq 'scanDB') {
		my $analysisStrg=param('ANA_ID');
		system "./scanDatabank.pl $userID $analysisStrg";
		open(OUT,">$promsPath{logs}/scan_$userID.end");
		print OUT "Scan launched by $userID";
		close OUT;
		sleep 1;
	}

	print qq
|</BODY>
</HTML>
|;
	exit;
}

sub getPercolatorScores {
		my ($anaFile,$refPercolatorScores)=@_;
		open(TARGET,"$promsPath{tmp}/batch/$userID/.percolator/$anaFile.target.pop") or die;
		while (my $line = <TARGET>) {
				next unless $line =~ /^query/;
				chomp($line);
				my ($queryString,$score,$qValue,$PEP,$peptide,$proteinIds)=split(/\t/,$line);
				my @queryInfos=split(/:|;/,$queryString); # $queryString='query:<$queryNum>;rank:1'
				$refPercolatorScores->{'peptides'}{$queryInfos[1]}{$queryInfos[-1]}="PEP=$PEP,QVAL=$qValue," if $qValue < $maxFDR/100.0;
		}
		close TARGET;
		open(DECOY,"$promsPath{tmp}/batch/$userID/.percolator/$anaFile.decoy.pop") or die;
		while (my $line = <DECOY>) {
				next unless $line =~ /^query/;
				chomp($line);
				my ($queryString,$score,$qValue,$PEP,$peptide,$proteinIds)=split(/\t/,$line);
				my @queryInfos=split(/:|;/,$queryString); # $queryString='query:<$queryNum>;rank:1'
				$refPercolatorScores->{'decoy_peptides'}{$queryInfos[1]}{$queryInfos[-1]}="PEP=$PEP,QVAL=$qValue," if $qValue < $maxFDR/100.0;
		}
		close DECOY;
}

###############################################################################
###########################     CLASS   #######################################
###############################################################################

###############################################################################
###########################    PHENYX PARSING   ###############################
###############################################################################



	#############    PeptideProteinedHandler   ################

	package PeptideProteinedHandler; {
		our %PeptideProteinedHandler;
		sub new {
			my $type = $_[0];
			return  bless {}, $type ;
	    }
		sub init_pepProt {
			my $self = $_[0] ;
			$self->{'key'} = $_[1];
		}
		sub getKey_PepProt {
			my $self = $_[0] ;
			return $self->{'key'};
		}
		sub setStart {
			my $self = $_[0] ;
			$self->{'start'} = $_[1]+1 ;
		}
		sub setLink {
			my $self = $_[0] ;
			$self->{'Link'} = $_[1] ; #Link to peptide MatchObjet
		}
		sub setAaBefore {
			my $self = $_[0] ;
			$self->{'aaBefore'} = $_[1] ;
		}
		sub setAaAfter{
			my $self = $_[0] ;
			$self->{'aaAfter'} = $_[1] ;
		}
		sub getTrueScore {
			my $self = $_[0] ;
			my $remoteScore =${$self->{'Link'}}->getScore();
			my $remoteSelStatus =${$self->{'Link'}}->getSelStatus();
			my $returnScore = ($remoteSelStatus == 0)?$remoteScore:0;
			return $returnScore;
		}
		sub getActualSequence {
			my $self = $_[0] ;
			return ${$self->{'Link'}}->getActualSequence();
		}
		sub getRankNumberInfo {
			my $self = $_[0] ;
			return ${$self->{'Link'}}->getRankNumber();
		}
			sub getQueryNumberInfo {
			my $self = $_[0] ;
			return  ${$self->{'Link'}}->getQueryNumber();
		}
		sub getPosInfo() {
			my $self = $_[0] ;
			my $returnString = $self->{'start'}.",";
			$returnString .=(defined ($self->{'aaBefore'}))?$self->{'aaBefore'}.",":"-," ;
			$returnString .=(defined ($self->{'aaAfter'}))?$self->{'aaAfter'}:"-" ;
			return $returnString ;
		}



	}

	#############    ProteinMatchedHandler   ################

	package ProteinMatchedHandler; {
		our %ProteinMatchedHandler;
		sub new {
			my $type = $_[0];
			return  bless {}, $type ;
	    }
		sub init_prot {
			my $self = $_[0] ;
			$self->{'key'} = $_[1];
		}
		sub getKey_prot {
			my $self = $_[0] ;
			return $self->{'key'};
		}
		sub addPeptide {
			my $self = $_[0] ;
			my $localPeptideProteined = $_[1] ;
			$self->{'peptideList'}{$localPeptideProteined->getKey_PepProt()}= \$localPeptideProteined ;
		}
		#sub updateGlobalProtSelStatus {
		#	my $self = $_[0] ;
		#	$protSelStatus{$self->{'identifier'}}=-1;
		#}
		sub updateGlobalMaxProtScore {
			my $self = $_[0] ;
			my $score ;
			foreach my $peptide (keys (%{$self->{'peptideList'}})) {
				$score += ${$self->{'peptideList'}{$peptide}}->getTrueScore() ;
			}
			$maxProtScore{$self->{'identifier'}}=$score;

		}# altere "%maxProtMatch"   maxProtMatch{  n° de gi } {sequence + modif }   = ( nombre de match pour un( n° de gi) de (sequence +modif))
		sub updateGlobalMaxProtMatch {
			my $self = $_[0] ;
			foreach my $peptide (keys (%{$self->{'peptideList'}})) {
				$maxProtMatch{$self->{'identifier'}}{${$self->{'peptideList'}{$peptide}}->getActualSequence()}++;
			}
		}# utilise et altere et concatene "%rankProtMatch"     rankProtMatch{ (n° de query):(n° de rang) }  {  n° de gi   }   [     ]    = position de debut du peptide dans la prot ;
		sub updateGlobalRankProtMatch {
			my $self = $_[0] ;
			foreach my $peptide (keys (%{$self->{'peptideList'}})) {
				my $rank = ${$self->{'peptideList'}{$peptide}}->getRankNumberInfo();
				next unless $rank<=$maxRank;
				my $query = ${$self->{'peptideList'}{$peptide}}->getQueryNumberInfo();
				my $posInfo = ${$self->{'peptideList'}{$peptide}}->getPosInfo();
				push @{$rankProtMatch{"$query:$rank"}{$self->{'identifier'}}}, $posInfo;
			}
		}
		sub setIdentifier{
			my $self = $_[0] ;
			$self->{'identifier'} = $_[1] ;
		}
	}

		#############    ProteinMatchedSupporter  ################

	package ProteinMatchedSupporter; {
		our  %ProteinMatchedSupporter;
		sub new {
			my $type = $_[0];
			return  bless {}, $type ;
	    }
		sub add_prot {
			my $self = $_[0] ;
			my $localProteinMatchHandler = $_[1] ;
			$self->{'ProteinList'}{$localProteinMatchHandler->getKey_prot()}= \$localProteinMatchHandler ;
		}
		sub loadGlobalsProteins {
			my $self = $_[0] ;
			foreach my $Protein (keys (%{$self->{'ProteinList'}})) {
				${$self->{'ProteinList'}{$Protein}}->updateGlobalRankProtMatch ();
				${$self->{'ProteinList'}{$Protein}}->updateGlobalMaxProtMatch ();
				${$self->{'ProteinList'}{$Protein}}->updateGlobalMaxProtScore ();
				#${$self->{'ProteinList'}{$Protein}}-> updateGlobalProtSelStatus();
			}
		}
	}

	#############    PeptideMatchHandler   ################

	package PeptideMatchHandler ; {
		our %PeptideMatchHandler;
		sub new {
			my $type = $_[0];
			return  bless {}, $type ;
	    }
		sub init_pep {
			my $self = $_[0] ;
			$self->{'key'} = $_[1];
			#($self->{'queryNumber'})= ($self->{'key'} =~ /%0%(\d+)/);
			$self->{'peptideMass'} = 18.01056 ; # water need to be added from residue to obtain peptide mass
	    }
		sub getQueryNumber {
			my $self = $_[0] ;
			return $self->{'queryNumber'};
		}
		sub getRankNumber{
			my $self = $_[0] ;
			return $self->{'rank'};
		}
		sub getKey_pep {
			my $self = $_[0] ;
			return $self->{'key'};
		}
		sub getScore {
			my $self = $_[0] ;
			return $self->{'score'}
		}
		sub getSelStatus {
			my $self = $_[0] ;
			return $self->{'SEL'}
		}
		sub getCharge {
			my $self = $_[0] ;
			return $self->{'charge'} ;
		}
		sub getActualSequence {
			my $self = $_[0] ;
			return $self->{'actualSequence'} ;
		}
		sub compute_pepMatchHandler {
			my $self = $_[0] ;
			my $sequence =$self->{'sequence'} ;
			foreach my $AA (split //, $sequence) {
				$self->{'peptideMass'} += $massValueAA{$AA};
			}
			my $i = -1;
			foreach my $modif (split /:/ ,$self->{'modif'}) {
				$i++ ;
				next unless ($modif) ;
				my ($modifName,$type,$delta)=${$self->{'generalInfoLink'}}->getModInfo($modif);
				$self->{'peptideMass'} += $delta;
				if ($modif ne $modifName) { #name changed for Mascot compatibility
					$modifName =~ s/\)/:$i\)/ if $modifName !~ /\(N-term\)/;
				}
				else {$modifName = "$modifName ($i)"; }
				$self->{'varModString'}.=" + $modifName" if $type eq "variable" ;


			}
			$self->{'actualSequence'} =   $self->{'sequence'} ;
			$self->{'actualSequence'}.=  $self->{'varModString'} if ($self->{'varModString'});
		}
		sub setDataString {
			my $self = $_[0] ;
			$self->{'nbrMatchedProt'} = ($self->{'protList'})? scalar (@{$self->{'protList'}}) : 0;
			$self->{'missCleavage'}=0 if !defined($self->{'missCleavage'}); # Missed cleavage is not provided for not-matching peptides!!!!
			 # datastring =  SEL=(0 si best score, -1 si non best score) , MIS=(number of missed clivage),  CALC=(mass peptide calcule), DELT=(delta mass), SEQ=(peptide sequence),SC=(score),VMOD=(varmodstring),MATCH=(nombre de prot qui matchent)
			my $LSPstrg=($self->{'SEL'}==-1)? ',LSP=1' : ',LSP=0'; # Flag for lower-scoring peptide
			$self->{'dataString'}= "SEL=".$self->{'SEL'}."$LSPstrg,MIS=".$self->{'missCleavage'}.",CALC=".$self->{'peptideMass'}.",DELT=".$self->{'deltaMass'} ;
			$self->{'dataString'}.= ",SEQ=".$self->{'sequence'}.",SC=".$self->{'score'}.",MATCH=".$self->{'nbrMatchedProt'}.",";
			$self->{'dataString'}.="VMOD=".$self->{'varModString'}."," if $self->{'varModString'};
		}
		sub updateGlobalQueryInfo {
			my $self = $_[0] ;
			$queryInfo{$self->{'queryNumber'}}[($self->{'rank'})-1]= $self->{'dataString'} if  $self->{'rank'}<=$maxRank;	#Global
		}
		sub updateGlobalMaxQueryScore {
			my $self = $_[0] ;
			$maxQueryScore{$self->{'queryNumber'}} = $self->{'score'} if (!defined $maxQueryScore{$self->{'queryNumber'}}) || ($maxQueryScore{$self->{'queryNumber'}} < $self->{'score'} ) ;  #Global
		}
		sub updateGlobalMatchList {
			my $self = $_[0] ;
			foreach my $identifier (@{$self->{'protList'}}) {
				$matchList{$identifier}{$self->{'actualSequence'}} = 1 ;  #Global
			}
		}
		sub setRankNumber{
			my $self = $_[0] ;
			$self->{'rank'}= $_[1] ;
		}
		sub setMisCleavage {
			my $self = $_[0] ;
			$self->{'missCleavage'} = $_[1] ;
		}
		sub setSequence {
			my $self = $_[0] ;
			$self->{'sequence'} = $_[1] ;
		}
		sub setScore {
			my $self = $_[0] ;
			$self->{'score'} = $_[1] ;
		}
		sub setDelta {
			my $self = $_[0] ;
			$self->{'deltaMass'} = $_[1] ;
		}
		sub setCharge {
			my $self = $_[0] ;
			$self->{'charge'} = $_[1] ;
		}
		sub setModif {
			my $self = $_[0] ;
			my $modifString = $_[1] ;
			$self->{'modif'}=$_[1] ;
		}
		sub setRefSpectra {
			my $self = $_[0] ;
			$self->{'RefSpectra'} = $_[1] ;
			($self->{'queryNumber'})=($self->{'RefSpectra'}=~/cmpd_(\d+)/);

		}
		sub setGeneralInfoLink {
			my $self = $_[0] ;
			$self->{'generalInfoLink'} = $_[1] ;
		}
		sub addProt {
			my $self = $_[0] ;
			my $newProt =  $_[1] ;
			push @{$self->{'protList'}}, $newProt ;
		}
	}

	#############    PeptideMatchSupporter   ################
	package PeptideMatchSupporter ; {
		our %PeptideMatchSupporter ;
		sub new {
			my $type = shift;
			return bless {}, $type;
	    }
		sub setGeneralInfoLink {
			my $self = $_[0] ;
			$self->{'generalInfoLink'} = $_[1] ;
		}
		sub add_peptide {
			my $self = $_[0] ;
			my $localPeptideMatchHandler = $_[1] ;
			$self->{'peptideList'}{$localPeptideMatchHandler->getKey_pep()}= \$localPeptideMatchHandler ;
		}
		sub setPeptideMatchHandlerMisCleavage {
			my $self = $_[0] ;
			my $MisCleavage = $_[1] ;
			my $PeptideMatchKey = $_[2];
			${$self->{'peptideList'}{$PeptideMatchKey}}->setMisCleavage($MisCleavage);
		}
		sub getCharge {
			my $self = $_[0] ;
			my $QueryNum = $_[1] ;
			my $charge ;
			foreach my $peptideKey (keys (%{$self->{'peptideList'}})) {
				return ${$self->{'peptideList'}{$peptideKey}} ->getCharge() if $QueryNum ==(${$self->{'peptideList'}{$peptideKey}} ->getQueryNumber()) ;
			}
		}
		sub getPeptideLink {
			my $self = $_[0] ;
			my $PeptideKey = $_[1];
			return $self->{'peptideList'}{$PeptideKey};
		}
		sub compute_pepMatchSupp {
			my $self = $_[0] ;
			my (%localPeptideSatus, %localRankStatus);
			foreach my $peptideKey (keys (%{$self->{'peptideList'}})) {
				my $queryNumber = ${$self->{'peptideList'}{$peptideKey}} ->getQueryNumber() ;
				$localRankStatus{$queryNumber}{$peptideKey}= ${$self->{'peptideList'}{$peptideKey}} ->getScore() ;
			}
			foreach my  $queryNumber (keys (%localRankStatus)) {
				my $rank=1 ;
				foreach my $peptideKey (sort {$localRankStatus{$queryNumber}{$b} <=>$localRankStatus{$queryNumber}{$a}|| $b cmp $a } keys %{$localRankStatus{$queryNumber}}) {
					${$self->{'peptideList'}{$peptideKey}} ->setRankNumber($rank) ;
					$rank++ ;
				}
			}

			foreach my $peptideKey (keys (%{$self->{'peptideList'}})) {
				${$self->{'peptideList'}{$peptideKey}} ->setGeneralInfoLink($self->{'generalInfoLink'}) ;
				${$self->{'peptideList'}{$peptideKey}} -> compute_pepMatchHandler () ;
			}
			foreach my $objet (keys (%{$self->{'peptideList'}})) {
				my $actualSequence = ${$self->{'peptideList'}{$objet}}->{'actualSequence'};
				my $localSel = 0;
				foreach my $keybis (keys %{$localPeptideSatus{$actualSequence}}) {
					if  ($localPeptideSatus{$actualSequence}{$keybis}{"score"}<  ${$self->{'peptideList'}{$objet}}->{'score'} ) {
						#$localPeptideSatus{$actualSequence}{$keybis}{"SEL"}=-1 ;
						${$self->{'peptideList'}{$keybis}}->{"SEL"} =-1;
				    }
					else {
						$localSel = -1;
						last ;
					}
				}
				$localPeptideSatus{$actualSequence}{$objet}{"score"} = ${$self->{'peptideList'}{$objet}}->{'score'} ;
				${$self->{'peptideList'}{$objet}}->{'SEL'} = $localSel ;
				$numValid{${$self->{'peptideList'}{$objet}}->{'queryNumber'}}++;	#GLOBAL
			}
			foreach my $pepTideKey (keys (%{$self->{'peptideList'}})) {
				${$self->{'peptideList'}{$pepTideKey}} ->setDataString() ;
			}
		}
		sub loadGlobalsPeptides {
			my $self = $_[0] ;
			foreach my $peptideKey (keys (%{$self->{'peptideList'}})) {
				${$self->{'peptideList'}{$peptideKey}} ->updateGlobalQueryInfo () ;
				${$self->{'peptideList'}{$peptideKey}} ->updateGlobalMaxQueryScore () ;
				${$self->{'peptideList'}{$peptideKey}} ->updateGlobalMatchList () ;
			}
		}
	}

	 #############    QueryInfoHandler   ################
	package QueryInfoHandler ; {
		our %QueryInfoHandler ;
		sub new {
			my $type = $_[0];
			return  bless {}, $type ;
	    }
		sub init_query {
			my $self = $_[0] ;
			#$self->{'queryNumber'} = $_[1];
			$self->{'key'} = $_[2];
			($self->{'queryNumber'}) = ($self->{'key'} =~ /cmpd_(\d+)/);
		}
		sub getKey_query {
			my $self = $_[0] ;
			return $self->{'key'};
		}
		sub fetchGlobal {
			my $self = $_[0] ;
			my $PeptideLink = $_[1] ;
			($elutionTime{$self->{'queryNumber'}}) = (($self->{'queryTittle'}) =~/Elution:\s(\d+\.?\d*)/ ) ; #global
			$massObs{$self->{'queryNumber'}} = sprintf("%.6f", $self->{'moz'});   #global
			if ($self->{'charge'} =~ /\D/) {
				$self->{'charge'}=$$PeptideLink->getCharge($self->{'queryNumber'}) }
			#$massExp{$self->{'queryNumber'}}=1;	}   #global
			#else
			{	$massExp{$self->{'queryNumber'}} =sprintf("%.6f", ($self->{'moz'} * $self->{'charge'} - $self->{'charge'}* 1.007825) ) ;} 	#global
			$charge{$self->{'queryNumber'}} = ($self->{'charge'}) ? $self->{'charge'} : $self->{'charges'};
		}
		sub setTittle {
			my $self = $_[0] ;
			$self->{'queryTittle'} = $_[1];
		}
		sub setMoz {
			my $self = $_[0] ;
			$self->{'moz'} = $_[1];
		}
		sub setCharge {
			my $self = $_[0] ;
			$self->{'charge'} = $_[1];
		}
		sub addInterpretation {
			my $self = $_[0] ;
			push @{$self->{'charge'}}, $_[1];
		}
	}

	 #############    QueryInfoSupporter   ################
	package QueryInfoSupporter ; {
		our %QueryInfoSupporter;
		sub new {
				my $type = $_[0] ;
				return bless {}, $type;
		    }
		sub add_query {
			my $self = $_[0] ;
			my $localQueryInfoHandler = $_[1] ;
			$self->{'QueryList'}{$localQueryInfoHandler->getKey_query()}= \$localQueryInfoHandler ;
		}
		sub compute_query {
			my $self = $_[0] ;
			foreach my $objet (keys (%{$self->{'QueryList'}})) {
				${$self->{'QueryList'}{$objet}} -> fetchGlobal($self->{'PeptideLink'}) ;
			}
		}
		sub setPeptideListLink {
			my $self = $_[0] ;
			$self->{'PeptideLink'} = $_[1];
		}
	}

	#############    SpectrumListHandler  ################
	package SpectrumListHandler; {
		our %SpectrumListHandler;
		sub new {
			my $type = shift;
			return bless {}, $type;
	    }
		sub init_spectrum {
			my $self = $_[0] ;
			$self->{'key'}  = $_[1] ;
			($self->{'queryNumber'}) = ($self->{'key'} =~ /sample_0%cmpd_(\d+)/) ;
		}
		sub getKey_spectrum {
			my $self = $_[0] ;
			return $self->{'key'};
		}
		sub getQueryNumber {
			my $self = $_[0] ;
			return $self->{'queryNumber'};
		}
		sub getTittle {
			my $self = $_[0] ;
			return $self->{'tittle'};
		}
		sub getPeakList {
			my $self = $_[0] ;
			return $self->{'peakList'};
		}
		sub setTittle {
			my $self = $_[0] ;
			 $self->{'tittle'} = $_[1] ;
		}
		sub setPeakList {
			my $self = $_[0] ;
			my $peakList = $self->{'peakString'} ;
			my $peakString = '';
			foreach my $line  (split /\n/, $peakList ) {
				$line =~ s/^\D+//;
				my ($moz, $intensity, $charge) = (split /\s+/, $line) ;
				next unless ($moz) ;
				$peakString .="$moz:$intensity,";
			}
			chop ($peakString) ;
			$self->{'peakList'} = $peakString;
		}

		sub addPeakList {
			my $self = $_[0] ;
			my $localPeakString = $_[1] ;
			$self->{'peakString'} .= $localPeakString;
		}
	}

	 #############    GeneralModResHandler  ################
	package SpectrumListSupporter; {
		our %SpectrumListSupporter;
		sub new {
			my ($type, $msfilename) = @_;
			my $self = bless {}, $type;
			$self->{'msFilename'} = $msfilename;
			$self->{'queriesNum'} = 0;
			return $self;
	    }
		sub add_spectrum {
			my $self = $_[0] ;
			my $localSpectrumListHandler = $_[1] ;
			$self->{'spectrumList'}{$localSpectrumListHandler->getKey_spectrum()}= \$localSpectrumListHandler ;
			($self->{'queriesNum'})++;
		}
		sub setGeneralInfoLink {
			my $self = $_[0] ;
			$self->{'generalInfoLink'} = $_[1] ;
		}
		sub getQueryNum {
			my $self = $_[0] ;
			return $self->{'queriesNum'};
		}
		sub createPgf {
			my $self = $_[0] ;
			my ($queryLikeNum , $tittle, $peakList);
			my $newFile = $self->{'msFilename'};
			$newFile =~ s/\.xml/\.pgf/;
			my %usedMod = ${$self->{'generalInfoLink'}}->getUsedModRes();
			my %parameters= ${$self->{'generalInfoLink'}}->getParametersInfo() ;
			open (DESCP, ">$promsPath{valid}/ana_$analysisID/$newFile") ; #et encore du Global ;

			########### Search parameter section #####################
			print DESCP "--gc0p4Jq0M2Yt08jU534c0p\n";
			print DESCP "Content-Type: peakListPhenyx; name=\"parameters\"\n";
			foreach my $parameter (keys (%parameters)) {
				my $parametersValue = $parameters{$parameter} ;
				print DESCP "$parameter=$parametersValue\n";
			}
			print DESCP "queries=".$self->getQueryNum()."\n";
			my $fixedModString ='' ;
			my $varModString  = '' ;
			foreach my $modRes (sort{lc($a) cmp lc($b)} keys %usedMod) { # sort added by PP 09/11/11
				$varModString  .=$usedMod{$modRes}{'name'}.", " if $usedMod{$modRes}{'type'} eq "variable" ;
				$fixedModString .=$usedMod{$modRes}{'name'}.", "  if $usedMod{$modRes}{'type'} eq "fixed" ;
			}
			chop ($varModString); chop ($varModString);
			chop ($fixedModString); chop ($fixedModString);
			print DESCP "MODS=$fixedModString\n";
			print DESCP "IT_MODS=$varModString\n";


			############### Mass section ##########################"
			print DESCP "--gc0p4Jq0M2Yt08jU534c0p\n";
			print DESCP "Content-Type: peakListPhenyx; name=\"masses\"\n";

			foreach my $modRes (keys (%usedMod)) {
				print DESCP "delta1=".$usedMod{$modRes}{'delta'}.",".$usedMod{$modRes}{'name'}."\n" if $usedMod{$modRes}{'type'} eq "variable" ;
				print DESCP "Fixed=".$usedMod{$modRes}{'delta'}.",".$usedMod{$modRes}{'name'}.",".$usedMod{$modRes}{'residu'}."\n"  if $usedMod{$modRes}{'type'} eq "fixed" ;
			}

			################# MS/MS peak section #################
			print DESCP "--gc0p4Jq0M2Yt08jU534c0p\n\n";

			foreach my $singleSpectrum (keys (%{$self->{'spectrumList'}})) {
				$queryLikeNum = ${$self->{'spectrumList'}{$singleSpectrum}} ->getQueryNumber ();
				$tittle = ${$self->{'spectrumList'}{$singleSpectrum}} ->getTittle ();
				$peakList = ${$self->{'spectrumList'}{$singleSpectrum}} ->getPeakList ();
				print DESCP "--gc0p4Jq0M2Yt08jU534c0p\n";
				print DESCP "Content-Type: peakListPhenyx; name=\"query$queryLikeNum\"\n";
				print DESCP "title=$tittle\n";
				print DESCP "Ions1=$peakList\n";
				print DESCP "\n";
			}
			close (DESCP) ;
		}
	}

	 #############    GeneralModResHandler  ################
	package GeneralModResHandler ;{
		our %GeneralModResHandler;

		sub new {
			my $type = shift;
			return bless {}, $type;
	    }
		sub init_mod {
			my $self = $_[0] ;
			$self->{'name'}  = $_[1] ;
		}
		sub setResidu {
			my $self = $_[0] ;
			$self->{'residu'}  = $_[1] ;
		}
		sub setType {
			my $self = $_[0] ;
			$self->{'type'}  = $_[1] ;
		}
		sub setDelta {
			my $self = $_[0] ;
			$self->{'delta'}  = $_[1] ;
		}
		sub getName {
			my $self = $_[0] ;
			return $self->{'name'}   ;
		}
		sub getDelta {
			my $self = $_[0] ;
			return $self->{'delta'}   ;
		}
		sub getMascotCompatibleName {
			my %modiftranslation = ('Oxidation_M' =>'Oxidation (M)', 'ACET_nterm'=>'Acetyl (N-term)', 'PHOS' => 'Phospho (STY)' ) ;
			my $self = $_[0] ;
			my $MascotCompatibleName = (exists ($modiftranslation{$self->{'name'}}))?$modiftranslation{$self->{'name'}}:$self->{'name'} ;
			return $MascotCompatibleName;
		}
		sub getResidu {
			my $self = $_[0] ;
			return $self->{'residu'}   ;
		}
	}

	 #############    GeneralInfoSupporter  ################
	package GeneralInfoSupporter; {
		our %GeneralInfoSupporter;
		sub new {
			my $type = shift;
			return bless {}, $type;
	    }
		sub addMosRes {
			my $self = $_[0] ;
			my $MosRes = $_[1] ;
			$self->{'ModResList'}{$MosRes->getName()}= \$MosRes ;
		}
		sub setModType {
			my $self = $_[0] ;
			my $ModType = $_[1] ;
			my $MosName = $_[2] ;
			${$self->{'ModResList'}{$MosName}}->setType($ModType) ;
			$self->{'ModResUsed'}{$MosName}{'name'}=${$self->{'ModResList'}{$MosName}}->getMascotCompatibleName() ;
			$self->{'ModResUsed'}{$MosName}{'type'}=$ModType;
			$self->{'ModResUsed'}{$MosName}{'delta'}=${$self->{'ModResList'}{$MosName}}->getDelta() ;
			$self->{'ModResUsed'}{$MosName}{'residu'}=${$self->{'ModResList'}{$MosName}}->getResidu() ;
		}
		sub addParameters {
			my $self = $_[0] ;
			my $parameters = $_[1] ;
			my $parametersValue = $_[2] ;
			$self->{'parameters'}{$parameters} = $parametersValue ;
		}


		sub getModInfo {
			my $self = $_[0] ;
			my $MosName = $_[1] ;
			return $self->{'ModResUsed'}{$MosName}{'name'},$self->{'ModResUsed'}{$MosName}{'type'},$self->{'ModResUsed'}{$MosName}{'delta'};
		}
		sub getUsedModRes {
			my $self = $_[0] ;
			return %{$self->{'ModResUsed'}} ;
		}
		sub getParametersInfo {
			my $self = $_[0] ;
			return %{$self->{'parameters'}} ;
		}
	}

	  #############    PhenyxSAXHandler   ################
	package PhenyxSAXHandler; {
		our %PhenyxSAXHandler;
		my  $Element;
		my %Elements ;
		my ($queryHandled, $queryListObj, $PeptideMatch, $PeptideListObj, $ProteinMatch, $proteinListObj, $peptideProt,$spectrumListObj,$spectrum, $generalInfo) ;  #Objets
		my ($modRes);
		sub new {
			my ($type, $msFilename) = @_;
			my $self = bless ({}, $type);
			$self->{'msFilename'} = $msFilename;
			return $self;
	    }
		sub start_document {
			my ($self, $doc) = @_;
			$queryListObj = QueryInfoSupporter -> new ();
			$PeptideListObj = PeptideMatchSupporter -> new ();
			$proteinListObj = ProteinMatchedSupporter->new ();
			$spectrumListObj = SpectrumListSupporter->new ($self->{'msFilename'});
			$generalInfo = GeneralInfoSupporter ->new () ;
			$spectrumListObj->setGeneralInfoLink(\$generalInfo);
			$PeptideListObj->setGeneralInfoLink(\$generalInfo);
			$queryListObj->setPeptideListLink(\$PeptideListObj);
		}
		sub end_document {
			print (".") ;
			$PeptideListObj->compute_pepMatchSupp () ; #data Processing
			$queryListObj->compute_query(); #data Processing
			$PeptideListObj->loadGlobalsPeptides () ; #updating global variables
			$proteinListObj->loadGlobalsProteins(); #updating global variables
			$spectrumListObj->createPgf () ;#creating PGF file, containing data for drawing spectrum
		}
		sub start_element {
			my ($self, $el) = @_;
			$Element = $el->{'Name'};
			$Elements{$el->{'Name'}}=1 ;
			$self->{'characters'}{'Data'} ='';
			if ($Element eq "compoundInfo") {
				$queryHandled  = QueryInfoHandler->new ();
				$queryHandled ->init_query($el->{'Attributes'}{'{}compoundNumber'}{'Value'},$el->{'Attributes'}{'{}key'}{'Value'}) ;

			}
			elsif ($Element eq "PeptideMatchDef") {
				$PeptideMatch  = PeptideMatchHandler-> new ();
				$PeptideMatch  -> init_pep($el->{'Attributes'}{'{}key'}{'Value'}) ;
			}
			elsif ($Element eq "DBMatch") {
				$ProteinMatch  = ProteinMatchedHandler-> new ();
				$ProteinMatch -> init_prot($el->{'Attributes'}{'{}key'}{'Value'}) ;
			}
			elsif ($Element eq "PeptideMatch") {
				$peptideProt  = PeptideProteinedHandler-> new ();
				$peptideProt -> init_pepProt($el->{'Attributes'}{'{}ref'}{'Value'}) ;
				$peptideProt -> setLink($PeptideListObj->getPeptideLink($el->{'Attributes'}{'{}ref'}{'Value'})) ;
			}
			elsif ($Element eq "ple:peptide") {
				$spectrum = SpectrumListHandler->new();
				$spectrum -> init_spectrum($el->{'Attributes'}{'{}key'}{'Value'}) ;
			}
			elsif (exists ($Elements{"OneDatabaseResult"})) {
				if ($Element eq "OneDatabaseResult") {
					$generalInfo->addParameters("release", ($el->{'Attributes'}{'{}name'}{'Value'})."_".($el->{'Attributes'}{'{}release'}{'Value'})) ;
				}
				elsif (exists ($Elements{"parameters"})) {
					if (($Element eq "OneParam")&& (($el->{'Attributes'}{'{}name'}{'Value'}) eq "fragTol")) {
						$generalInfo->addParameters("ITOL", ($el->{'Attributes'}{'{}value'}{'Value'})) ;
					}
				}
			}

			elsif (exists ($Elements{"inSilicoDefinitions"})) {
				if ($Element eq "oneModRes") {
					$modRes = GeneralModResHandler->new();
					$modRes -> init_mod($el->{'Attributes'}{'{}name'}{'Value'}) ;
				}
				elsif ((exists ($Elements{"oneModRes"})) && $Element eq "delta") { $modRes->setDelta($el->{'Attributes'}{'{}monoisotopic'}{'Value'}) }
				elsif ((exists ($Elements{"oneModRes"}))&& $Element eq "residues") {$modRes->setResidu($el->{'Attributes'}{'{}aa'}{'Value'})   }
			}
			elsif (exists ($Elements{"identificationAlgos"})) {
				if  ($Element eq "OneModif") {
					$generalInfo->setModType($el->{'Attributes'}{'{}mode'}{'Value'},$el->{'Attributes'}{'{}name'}{'Value'}) ;
				}
				elsif ($Element eq "scoring") {
					$generalInfo->addParameters("INSTRUMENT", $el->{'Attributes'}{'{}instrument'}{'Value'});
					$phenyxInstrument = $el->{'Attributes'}{'{}instrument'}{'Value'} ;
					if (exists ($instrumentConver{$phenyxInstrument})) {$phenyxInstrument =  $instrumentConver{$phenyxInstrument} ; };
				}
				elsif ($Element eq "cleavageEnzyme") {
					$generalInfo->addParameters("CLE", $el->{'Attributes'}{'{}name'}{'Value'});
					$generalInfo->addParameters("PFA", $el->{'Attributes'}{'{}missedCleavage'}{'Value'});
				}
				elsif ($Element eq "peptTolerance") {
					$generalInfo->addParameters("TOL", $el->{'Attributes'}{'{}value'}{'Value'});
					$generalInfo->addParameters("TOLU", $el->{'Attributes'}{'{}unit'}{'Value'});
				}

			}

		}

		sub end_element {
			my ($self, $el) = @_;
			$self->characters2 ($self->{'characters'});

			if ($el->{'Name'} eq "compoundInfo") {
				$queryListObj->add_query($queryHandled) ;
				$queryHandled ='' ;
			}
			elsif ($el->{'Name'} eq "oneModRes") {
				$generalInfo->addMosRes($modRes);
				$modRes ='' ;
			}

			elsif ($el->{'Name'} eq "PeptideMatchDef") {
				$PeptideListObj->add_peptide($PeptideMatch);
				$PeptideMatch ='' ;
			}

			elsif ($el->{'Name'} eq "DBMatch") {
				$proteinListObj->add_prot($ProteinMatch);
				$ProteinMatch ='' ;
			}
			elsif ($el->{'Name'} eq "PeptideMatch") {
				$ProteinMatch->addPeptide($peptideProt);
				$peptideProt ='' ;
			}
			elsif ($el->{'Name'} eq "ple:peaks") {
				$spectrum->setPeakList () ;
			}

			elsif ($el->{'Name'} eq "ple:peptide") {
				$spectrumListObj->add_spectrum($spectrum);
				$spectrum ='' ;
			}

			$self->{'characters'}{'Data'} ='';
			$Element = '' ;
			delete ($Elements{$el->{'Name'}});
		}

		sub characters {
			my ($self, $el) = @_;
			$self->{'characters'}{'Data'} .= $el->{'Data'} if defined ($el->{'Data'});
		}


		sub characters2 {
			my ($self, $el) = @_;

			if (exists ($Elements{"compoundInfo"})) {
				if ($Element eq "description" ) {$queryHandled->setTittle($el->{'Data'});}
				elsif ($Element eq "moz") {$queryHandled->setMoz($el->{'Data'});}
				elsif ($Element eq "charges") {$queryHandled->setCharge($el->{'Data'});}
				#elsif ($Element eq "onePeptideMatchesDefRef") {$queryHandled->addInterpretation($el->{'Data'});}

			}
			elsif (exists ($Elements{"header"})) {
				if ($Element eq "jobId" ) { $generalInfo->addParameters("ID", $el->{'Data'}); }
				elsif ($Element eq "title" ) { $generalInfo->addParameters("TITTLE", $el->{'Data'}); }
				elsif ($Element eq "date" ) { $generalInfo->addParameters("DATE", $el->{'Data'}); }

			}

			elsif (exists ($Elements{"olavJobSubmission"})) {
				if ($Element eq "clientPath" ) { $phenyxSourceFile = $el->{'Data'}; $generalInfo->addParameters("FILE", $el->{'Data'});}
				elsif ($Element eq "taxoCriterion" ) { $phenyxTaxoNumber = $el->{'Data'}; $generalInfo->addParameters("TAXO_ID", $el->{'Data'});}
			}
			elsif (exists ($Elements{"PeptideMatchDef"})) {
				if ($Element eq "sequence")	{$PeptideMatch->setSequence($el->{'Data'});}
				elsif ($Element eq "peptZScore")	{$PeptideMatch->setScore($el->{'Data'});}
				elsif ($Element eq "deltaMass")		{$PeptideMatch->setDelta($el->{'Data'}) unless exists $Elements{"ionicSeries"};}
				elsif ($Element eq "charge")		{$PeptideMatch->setCharge($el->{'Data'});}
				elsif ($Element eq "modifAt")		{$PeptideMatch->setModif($el->{'Data'});}
				elsif ($Element eq "spectrumRef")	{$PeptideMatch->setRefSpectra($el->{'Data'});}
				elsif ((exists ($Elements{"DBMatchesRefList"})) && ($Element eq "oneDBMatchesRef" )) {my ($fastaBase, $identifier)=split(/%/, $el->{'Data'}); $PeptideMatch->addProt($identifier);}
			}
			elsif (exists ($Elements{"DBMatch"})) {
				if (exists ($Elements{"PeptideMatch"})) {
					if ($Element eq "start" ) { $peptideProt->setStart($el->{'Data'})}
					elsif ($Element eq "missCleav" ) {$PeptideListObj->setPeptideMatchHandlerMisCleavage($el->{'Data'},$peptideProt->getKey_PepProt());}
					elsif ($Element eq "aaBefore" ) {$peptideProt->setAaBefore($el->{'Data'});}
					elsif ($Element eq "aaAfter" ) {$peptideProt->setAaAfter($el->{'Data'});}
				}
				elsif ($Element eq "AC" ) {$ProteinMatch->setIdentifier($el->{'Data'})}

			}
			elsif (exists ($Elements{"ple:peptide"})) {
				if ($Element eq "ple:PeptideDescr" ) {$spectrum->setTittle($el->{'Data'});}
				elsif ($Element eq "ple:peaks" ) {$spectrum->addPeakList($el->{'Data'});}
			}
		}
	}



#################################################################################
################################### MASCOT XML PARSING ##########################
#################################################################################


  #############    MascotXMLSAXHandler   ################
	package MascotXMLSAXHandler; {
		our %MascotXMLSAXHandler ;

		sub new {
			my ($type, $msFilename) = @_;
			my $self = bless ({}, $type);
			$self->{'msFilename'} = $msFilename;
			return $self;
		}
		sub start_document {
			my ($self, $doc) = @_;

		}
		sub end_document {
			my $self = $_[0];

		}
		sub start_element {
			my ($self, $el) = @_;
			$self->{'Element'} = $el->{'Name'};
			$self->{'Elements'}{$el->{'Name'}}=1 ;
			$self->{'characters'}{'Data'} = '';

			if (exists ( $self->{'Elements'}{'variable_mods'})) {
				if ($self->{'Element'} eq 'modification') {
					$self->{'localmodification'}{'identifier'} = $el->{'Attributes'}{'{}identifier'}{'Value'};
				}
			}
			elsif (exists ( $self->{'Elements'}{'masses'})) {
				if ($self->{'Element'} eq 'mass') {
					$self->{'localAA'}{'name'}= $el->{'Attributes'}{'{}name'}{'Value'} ;
				}
			}

			elsif (exists ( $self->{'Elements'}{'hits'})) {
				if ($self->{'Element'} eq 'protein') {
					$self->{'localprotein'}{'Access'} = $el->{'Attributes'}{'{}accession'}{'Value'} ;
					#$protSelStatus{$self->{'localprotein'}{'Access'}} = -1 ;
				}
				elsif ($self->{'Element'} eq 'peptide') {
					$self->{'localpeptide'}{'query'} = $el->{'Attributes'}{'{}query'}{'Value'} ;
					$self->{'localpeptide'}{'rank'} = $el->{'Attributes'}{'{}rank'}{'Value'} ;
				}
			}

			elsif (exists ( $self->{'Elements'}{'queries'})) {
				if 	($self->{'Element'} eq "queries") {
				################"Creating PGF file###############################
					my $newFile = $self->{'msFilename'};
					$newFile =~ s/\.xml/\.pgf/;
					open (DESCP, ">$promsPath{valid}/ana_$analysisID/$newFile") || die ("unable to create $promsPath{valid}/ana_$analysisID/$newFile \n") ; #et encore du Global ;

					###########"search parameter section #####################
					print DESCP "--gc0p4Jq0M2Yt08jU534c0p\n";
					print DESCP "Content-Type: peakListMascot; name=\"parameters\"\n";
					foreach my $parameter (keys (%{$self->{'Elements'}{'parameters'}})) {
						print DESCP "$parameter=".$self->{'Elements'}{'parameters'}{$parameter}."\n";
					}

					###############"masse section ##########################"
					print DESCP "--gc0p4Jq0M2Yt08jU534c0p\n";
					print DESCP "Content-Type: peakListMascot; name=\"masses\"\n";
					foreach my $modRes (keys (%{$self->{'VarMods'}})) {
						print DESCP "delta$modRes=".$self->{'VarMods'}{$modRes}{'delta'}.",".$self->{'VarMods'}{$modRes}{'name'}."\n" ;
					}

					###############"masseAA section ##########################"
					print DESCP "--gc0p4Jq0M2Yt08jU534c0p\n";
					print DESCP "Content-Type: peakListMascot; name=\"AAmass\"\n";
					foreach my $AA (keys ( %{$self->{'massesAA'}})) {
						print DESCP "$AA=".$self->{'massesAA'}{$AA}."\n" ;
					}

				}
				elsif ($self->{'Element'} eq "query") {
					print DESCP "\n--gc0p4Jq0M2Yt08jU534c0p\n";
					print DESCP "Content-Type: peakListMascot; name=\"query".$el->{'Attributes'}{'{}number'}{'Value'}."\"\n";

				}
			}
		}

		sub end_element {
			my ($self, $el) = @_;
			$self->characters2 ($self->{'characters'});

			$self->{'characters'}{'Data'} = '';
			if (exists ( $self->{'Elements'}{'variable_mods'})) {
				if ($el->{'Name'} eq 'modification') {
					$self->{'localmodification'}{'identifier'} = '';
				}
			}

			elsif (exists ( $self->{'Elements'}{'hits'})) {
				if ( $el->{'Name'} eq 'protein') {
					delete $self->{'localprotein'} ;
				}
				elsif ( $el->{'Name'} eq 'peptide') {
					if (($self->{'localpeptide'}{'rank'}<=$maxRank) && ($self->{'localpeptide'}{'score'} >= $minScore)) {
						if ($self->{'localpeptide'}{'var_mod_pos'}) {
							my ($modNter,$modCore, $modCter) = split (/\./, $self->{'localpeptide'}{'var_mod_pos'});
							my @arrayModCore = split (//, $modCore );
							my $varModString ;
							if ($modNter >0 ) {
								$varModString = " + ".$self->{'VarMods'}{$modNter}{'name'} ;
							}
							my $i = 0;
							foreach (@arrayModCore) {
								if ($arrayModCore[$i]>0) {
									my $localVarModsstring = $self->{'VarMods'}{$arrayModCore[$i]}{'name'} ;
									my $pos = $i+1;
									$localVarModsstring =~ s/\)/:$pos\)/;
									$varModString .=" + ". $localVarModsstring  ;
								}
							$i++ ;
							}
							if ($modCter >0 ) {
								$varModString .=" + ". $self->{'VarMods'}{$modCter}{'name'}  ;
							}

							$self->{'localpeptide'}{'varModString'} = $varModString ;
						}
						else {$self->{'localpeptide'}{'varModString'} = ''} ;

						$self->{'localpeptide'}{'actualSequence'} = $self->{'localpeptide'}{'sequence'};
						$self->{'localpeptide'}{'actualSequence'}.=" + ".$self->{'localpeptide'}{'varModString'} if $self->{'localpeptide'}{'varModString'} ;

						if (not defined ( $self->{'interpretedQueryList'}{$self->{'localpeptide'}{'query'}}{$self->{'localpeptide'}{'rank'}}{'info'})) {
							$self->{'interpretedQueryList'}{$self->{'localpeptide'}{'query'}}{$self->{'localpeptide'}{'rank'}}{'info'} = $self->{'localpeptide'} ;
							push (@{$self->{'peptideList'}{$self->{'localpeptide'}{'actualSequence'}}}, $self->{'localpeptide'}{'query'}.":".$self->{'localpeptide'}{'rank'}) ;

						}
						$self->{'proteinList'}{$self->{'localprotein'}{'Access'}}{'peptides'}{$self->{'localpeptide'}{'query'}.":".$self->{'localpeptide'}{'rank'}} =$self->{'localpeptide'}{'protInfo'};

						$self->{'interpretedQueryList'}{$self->{'localpeptide'}{'query'}}{$self->{'localpeptide'}{'rank'}}{'MATCH'}++ ;

						$matchList{$self->{'localprotein'}{'Access'}}{$self->{'localpeptide'}{'actualSequence'}}=1;
					}
					delete $self->{'localpeptide'} ;
				}

				if ( $el->{'Name'} eq 'hits') {

				##########SEL Status ###############"
					foreach my $actSequence (keys ( %{$self->{'peptideList'}})) {

						if (scalar (@{$self->{'peptideList'}{$actSequence}}) >1 ) {
							my $localBestScore ;
							foreach my $refs (@{$self->{'peptideList'}{$actSequence}}) {
								my ($query,$rank) = split (/:/, $refs);
								$localBestScore = $self->{'interpretedQueryList'}{$query}{$rank}{'info'}{'score'} if (not defined ($localBestScore)) || ($self->{'interpretedQueryList'}{$query}{$rank}{'info'}{'score'} > $localBestScore ) ;

							}

							foreach my $refs (@{$self->{'peptideList'}{$actSequence}}) {
								my ($query,$rank) = split (/:/, $refs);
								if ($self->{'interpretedQueryList'}{$query}{$rank}{'info'}{'score'} == $localBestScore) {
									$self->{'interpretedQueryList'}{$query}{$rank}{'SEL'}= 0 ;
								}
								else {
									$self->{'interpretedQueryList'}{$query}{$rank}{'SEL'}= -1 ;
								}
							}

						}
						else {
							my ($query,$rank) = split (/:/, $self->{'peptideList'}{$actSequence}[0]);
							$self->{'interpretedQueryList'}{$query}{$rank}{'SEL'}= 0 ;
						}
					}

				###########Setting queryInfo ############################
					foreach my $query (keys (%{$self->{'interpretedQueryList'}})) {
						foreach my $rank (keys (%{$self->{'interpretedQueryList'}{$query}})) {
							####setting datastirng ; ####
							my $datastring ;
							$datastring = "SEL=".$self->{'interpretedQueryList'}{$query}{$rank}{'SEL'};
							$datastring .= ($self->{'interpretedQueryList'}{$query}{$rank}{'SEL'}==-1)? ',LSP=1' : ',LSP=0'; # Flag for lower-scoring peptide
							$datastring .= ",MIS=".$self->{'interpretedQueryList'}{$query}{$rank}{'info'}{'MIS'};
							$datastring .= ",CALC=".$self->{'interpretedQueryList'}{$query}{$rank}{'info'}{'CALC'};
							$datastring .= ",DELT=".$self->{'interpretedQueryList'}{$query}{$rank}{'info'}{'DELTA'};
							$datastring .= ",SEQ=".$self->{'interpretedQueryList'}{$query}{$rank}{'info'}{'sequence'};
							$datastring .= ",SC=".$self->{'interpretedQueryList'}{$query}{$rank}{'info'}{'score'};
							$datastring .= ",VMOD=".$self->{'interpretedQueryList'}{$query}{$rank}{'info'}{'varModString'} if  $self->{'interpretedQueryList'}{$query}{$rank}{'info'}{'varModString'};
							$datastring .= ",MATCH=".$self->{'interpretedQueryList'}{$query}{$rank}{'MATCH'};
							####updating query Info  ####
							push (@{$queryInfo{$query}}, $datastring);
						}
					}
				###########Setting %numValid" ############################
					foreach my $query (keys (%{$self->{'interpretedQueryList'}})) {
						my $numValidInterpretation ;
							foreach my $rank (keys (%{$self->{'interpretedQueryList'}{$query}})) {
								$numValidInterpretation++ if $self->{'interpretedQueryList'}{$query}{$rank}{'SEL'} == 0;
							}
						$numValid{$query} = $numValidInterpretation ;
					}
				###########Setting %maxProtScore" ############################
					foreach my $identifier (keys( %{$self->{'proteinList'}})) {
						my $maxScore;
						foreach my $refs (keys %{$self->{'proteinList'}{$identifier}{'peptides'}}) {
							my ($query,$rank) = split (/:/, $refs);
							if ($self->{'interpretedQueryList'}{$query}{$rank}{'SEL'} == 0) {
								$maxScore += $self->{'interpretedQueryList'}{$query}{$rank}{'info'}{'score'} ;
							}
						}
						$maxProtScore{$identifier} = $maxScore ;

					}

				###########Setting %maxProtMatch" ############################
					foreach my $identifier (keys( %{$self->{'proteinList'}})) {
						foreach my $refs (keys %{$self->{'proteinList'}{$identifier}{'peptides'}}) {
							my ($query,$rank) = split (/:/, $refs);
							my $actualSequence =$self->{'interpretedQueryList'}{$query}{$rank}{'info'}{'actualSequence'} ;
							$maxProtMatch{$identifier}{$actualSequence} = 1; #sorry, I have no file with multiple !

						}
					}
				###########Setting%rankProtMatch" ############################
					foreach my $identifier (keys( %{$self->{'proteinList'}})) {
						foreach my $refs (keys %{$self->{'proteinList'}{$identifier}{'peptides'}}) {
							my ($query,$rank) = split (/:/, $refs);
							my $protRefString ;
							next unless $self->{'interpretedQueryList'}{$query}{$rank}{'SEL'} == 0 ;
							$protRefString = 		$self->{'proteinList'}{$identifier}{'peptides'}{$refs}{'position'};
							$protRefString .=",".	$self->{'proteinList'}{$identifier}{'peptides'}{$refs}{'before'};
							$protRefString .=",".	$self->{'proteinList'}{$identifier}{'peptides'}{$refs}{'after'};
							push (@{$rankProtMatch{$refs}{$identifier}}, $protRefString );

						}
					}
				}
			}

			elsif ( $el->{'Name'} eq 'queries') {
				close (DESCP) ;
			}
			delete ($self->{'Elements'}{$el->{'Name'}});
			$self->{'Element'} ='';
		}

		sub characters { #concatene Characters String, in case of multiple string return
			my ($self, $el) = @_;
			$self->{'characters'}{'Data'} .= $el->{'Data'};
		}

		sub characters2 { #processing of string between element, called at the beginning of end element
			my ($self, $el) = @_;

			if (exists ( $self->{'Elements'}{'variable_mods'})) {
				if (exists ( $self->{'Elements'}{'modification'})) {
					if ($self->{'Element'} eq 'name') {
						$self->{'VarMods'}{$self->{'localmodification'}{'identifier'}}{'name'}= $el->{'Data'};
					}
					if ($self->{'Element'} eq 'delta') {
						$self->{'VarMods'}{$self->{'localmodification'}{'identifier'}}{'delta'}= $el->{'Data'};
					}
				}
			}

			elsif (exists ( $self->{'Elements'}{'search_parameters'})) {
				$self->{'Elements'}{'parameters'}{$self->{'Element'}} = $el->{'Data'} ;
			}

			elsif (exists ( $self->{'Elements'}{'masses'})) {
				if ($self->{'Element'} eq 'mass') {
					$self->{'massesAA'}{$self->{'localAA'}{'name'}} = $el->{'Data'} ;
				}
			}

			elsif (exists ( $self->{'Elements'}{'hits'})) {

				if (exists( $self->{'Elements'}{'peptide'}) ) {
					if ($self->{'Element'} eq 'pep_exp_mz') {
						$massObs{$self->{'localpeptide'}{'query'}} = $el->{'Data'} if not defined ($massObs{$self->{'localpeptide'}{'query'}});
					}
					elsif ($self->{'Element'} eq 'pep_exp_mr') {
						$massExp{$self->{'localpeptide'}{'query'}} = $el->{'Data'}   if not defined ($massExp{$self->{'localpeptide'}{'query'}});
					}
					elsif ($self->{'Element'} eq 'pep_scan_title') {
						($elutionTime{$self->{'localpeptide'}{'query'}}) = ($el->{'Data'} =~ /Elution:\s(\d+\.?\d*)/) if not defined ($elutionTime{$self->{'localpeptide'}{'query'}});
					}
					elsif ($self->{'Element'} eq 'pep_score') {
						$self->{'localpeptide'}{'score'} = $el->{'Data'};
						$maxQueryScore{$self->{'localpeptide'}{'query'}} = $el->{'Data'}  if (not defined($maxQueryScore{$self->{'localpeptide'}{'query'}})) || (($maxQueryScore{$self->{'localpeptide'}{'query'}})<( $el->{'Data'}));
					}
					elsif ($self->{'Element'} eq 'pep_seq') {
						$self->{'localpeptide'}{'sequence'} = $el->{'Data'} ;
					}
					# elsif ($self->{'Element'} eq 'pep_var_mod') {
						# $self->{'localpeptide'}{'var_mod'} = $el->{'Data'} ;
					# }
					elsif ($self->{'Element'} eq 'pep_var_mod_pos') {
						$self->{'localpeptide'}{'var_mod_pos'} = $el->{'Data'} ;
					}
					elsif ($self->{'Element'} eq 'pep_miss') {
						$self->{'localpeptide'}{'MIS'} = $el->{'Data'} ;
					}
					elsif ($self->{'Element'} eq 'pep_calc_mr') {
						$self->{'localpeptide'}{'CALC'} = $el->{'Data'} ;
					}
					elsif ($self->{'Element'} eq 'pep_delta') {
						$self->{'localpeptide'}{'DELTA'} = $el->{'Data'} ;
					}
					elsif ($self->{'Element'} eq 'pep_res_before') {
						$self->{'localpeptide'}{'protInfo'}{'before'} = $el->{'Data'} ;
					}
					elsif ($self->{'Element'} eq 'pep_res_after') {
						$self->{'localpeptide'}{'protInfo'}{'after'} = $el->{'Data'} ;
					}
					elsif ($self->{'Element'} eq 'pep_start') {
						$self->{'localpeptide'}{'protInfo'}{'position'} = $el->{'Data'} ;
					}
				}

			}

			elsif (exists ( $self->{'Elements'}{'queries'})) {
				if 	($self->{'Element'} eq "StringTitle") {
					print DESCP "title=".$el->{'Data'}."\n";
				}
				elsif ($self->{'Element'} eq "StringIons1") {
					print DESCP "Ions1=".$el->{'Data'}."\n";
				}
			}
		}
	}

	###########################################################################################
	################################### PARAGON XML PARSING ###################################
	###########################################################################################
	package ParagonXMLHandler; {

		###> Be careful, those variables are global to that class... Need to be instantiate everytime a
		my (%matches,%anaInfo,%spectra,%repeatedSection);
		my (%bestScoreTarget,%bestScoreDecoy); # locals required for calculation of protein max score
		my %infoRank; # local to keep $infoRank{$matchID}{$identifier}{$beg} information for rank ProtMatch purposes
		my (%query2matches,%matches2proteins,%sequence2spectrum,%modFeatures);
		my (%allVmods,%elements,%specificity);# All the VARIABLE MODIFICATIONS that have been seen...
		my %search; # Keep all search information into that hash
		my %aaCode;

                sub new {#Constructor of the mzXMLHandler
                        my ($type,$msFilename,$anaVmods)= @_;
                        my $self = bless ({}, $type);
			$self->{'msFilename'}=$msFilename;
			$self->{'anaVmods'}=$anaVmods;
			my $dbh=&promsConfig::dbConnect;
			my ($parseRules)=$dbh->selectrow_array("SELECT PARSE_RULES FROM DATABANK,DATABANK_TYPE WHERE DATABANK.ID_DBTYPE=DATABANK_TYPE.ID_DBTYPE AND ID_DATABANK=$databankIDs[0]");
			$dbh->disconnect;
			my @rules=split(',:,',$parseRules);
			my ($idRule)=($rules[0]=~/ID=(.+)/);
			$self->{'idRule'}=$idRule;
                        return $self;
                }

                sub start_document {
                        my ($self) = @_;
			$self->{'queryNum'}=0;
			$self->{'tag'}=0;
			(%matches,%anaInfo,%spectra,%repeatedSection,%bestScoreTarget,%bestScoreDecoy,%infoRank,%query2matches,%matches2proteins,%sequence2spectrum,%allVmods,%elements,%search,%specificity)=((),(),(),(),(),(),(),(),(),(),(),(),());
                        #print "Starting reading XML file<BR>\n";
			%aaCode=('Alanine'=>'A',
			    'Arginine'=>'R',
			    'Asparagine'=>'N',
			    'AsnOrAsp'=>'B',
			    'Aspartic Acid'=>'D',
			    'Cysteine'=>'C',
			    'Glutamine'=>'Q',
			    'GlnOrGlu'=>'Z',
			    'Glutamic Acid'=>'E',
			    'Glycine'=>'G',
			    'Histidine'=>'H',
			    'Isoleucine'=>'I',
			    'Leucine'=>'L',
			    'Lysine'=>'K',
			    'Methionine'=>'M',
			    'Phenylalanine'=>'F',
			    'Proline'=>'P',
			    'Serine'=>'S',
			    'Threonine'=>'T',
			    'Tryptophan'=>'W',
			    'Tyrosine'=>'Y',
			    'Valine'=>'V',
			    'Selenocysteine'=>'U',
			    );
                }

                sub end_document {
                        my ($self) = @_;
                        #print "Finishing reading XML file<BR>\n";
                }

                sub start_element {# Read the elements of the XML
                        my ($self, $element) = @_;
                        if($element->{'Name'} eq "SPECTRUM"){ # Some spectrum are found with charge=0 but it is the precursor, then, match information do have a positive charge
			#if($element->{'Name'} eq "SPECTRUM" && $element->{'Attributes'}{'{}charge'}{'Value'} > 0){ # Some spectrum are found with charge=0
				$self->{'isOnSpectrum'}=1;
				$self->{'queryNum'}++;
				$self->{'rankNum'}=0;# For each query, the rank has to be initialised
                                $self->{'spectrumID'}=$element->{'Attributes'}{'{http://www.w3.org/XML/1998/namespace}id'}{'Value'};# ex: 2.1.1.2199.2.2
				#$extSpectrum{$self->{'queryNum'}}=$self->{'spectrumID'};
                                #$spectra{$self->{'spectrumID'}}{'spectrumCharge'}=$element->{'Attributes'}{'{}charge'}{'Value'};
                                #$spectra{$self->{'spectrumID'}}{'spectrumPrecursorElution'}=$element->{'Attributes'}{'{}precursorelution'}{'Value'};
                                #$spectra{$self->{'spectrumID'}}{'spectrumPrecursorMass'}=$element->{'Attributes'}{'{}precursormass'}{'Value'};
				###> Mixed-DB (decoy)... same values given to the queries and decoy ones
				$charge{$self->{'queryNum'}}=$element->{'Attributes'}{'{}charge'}{'Value'};
				$massObs{$self->{'queryNum'}}=$element->{'Attributes'}{'{}precursormass'}{'Value'};
				$massExp{$self->{'queryNum'}}=($element->{'Attributes'}{'{}precursormass'}{'Value'}-1.007825032)*$element->{'Attributes'}{'{}charge'}{'Value'};
				$elutionTime{$self->{'queryNum'}}=sprintf "et%.2f;sp$element->{'Attributes'}{'{http://www.w3.org/XML/1998/namespace}id'}{'Value'};",$element->{'Attributes'}{'{}elution'}{'Value'};# 2 type of elutions: precursorElution or Elution...
				if($self->{'isDecoySearch'}){
					$charge{"-$self->{'queryNum'}"}=$charge{$self->{'queryNum'}};
					$massObs{"-$self->{'queryNum'}"}=$massObs{$self->{'queryNum'}};
					$massExp{"-$self->{'queryNum'}"}=$massExp{$self->{'queryNum'}};
					$elutionTime{"-$self->{'queryNum'}"}=$elutionTime{$self->{'queryNum'}};
				}
				$self->{'tag'}++;
				if ($self->{'tag'} == 1000) {
					print ".";
					$self->{'tag'}=0;
				}
                        }elsif($element->{'Name'} eq "FEATURE_PROBABILITIES_SET" ){
				$self->{'isOnFeatures'}=1;
                        }elsif($self->{'isOnFeatures'} && $element->{'Name'} eq "MOD_FEATURE_SET" ){
				$self->{'isOnModFeature'}=1;
                        }elsif($self->{'isOnModFeature'} && $element->{'Name'} eq "MOD_FEATURE" ){
				$self->{'currentModDef'}=$element->{'Attributes'}{'{}mod'}{'Value'};
                        }elsif($self->{'isOnModFeature'} && $element->{'Name'} eq "OCCURRENCE" ){
				my $target=($element->{'Attributes'}{'{}target'}{'Value'})?$aaCode{$element->{'Attributes'}{'{}target'}{'Value'}}:'';
				if ($target) {
					$target.=($element->{'Attributes'}{'{}term_spec'}{'Value'} eq 'PepNTerm')? ";=":($element->{'Attributes'}{'{}term_spec'}{'Value'} eq 'ProtNTerm')? ";-": ($element->{'Attributes'}{'{}term_spec'}{'Value'} eq 'PepCTerm')? ";*" : ($element->{'Attributes'}{'{}term_spec'}{'Value'} eq 'ProtCTerm')? ";+" : "";
				}else{
					$target=($element->{'Attributes'}{'{}term_spec'}{'Value'} eq 'PepNTerm')? "=":($element->{'Attributes'}{'{}term_spec'}{'Value'} eq 'ProtNTerm')? "-": ($element->{'Attributes'}{'{}term_spec'}{'Value'} eq 'PepCTerm')? "*" : ($element->{'Attributes'}{'{}term_spec'}{'Value'} eq 'ProtCTerm')? "+" : "";
				}
				$specificity{$self->{'currentModDef'}}{$target}=1;
                        }elsif ($self->{'isOnSpectrum'} && $element->{'Name'} eq "MATCH"){#TODO -> use this information
				if ($element->{'Attributes'}{'{}confidence'}{'Value'}*100.0 >= $minScore) {
					$self->{'isOnMatch'}=1;
					$self->{'rankNum'}++; # New interpretation of this query (ie SPECTRUM) -> rank has to be updated.
					my $matchID=$element->{'Attributes'}{'{http://www.w3.org/XML/1998/namespace}id'}{'Value'};
					$query2matches{$self->{'queryNum'}}{$self->{'rankNum'}}=$matchID;
					$sequence2spectrum{"$element->{'Attributes'}{'{}seq'}{'Value'}:$element->{'Attributes'}{'{}charge'}{'Value'}"}{"$self->{'spectrumID'}"}=1;
					$self->{'matchID'}=$matchID;
					#$maxQueryScore{$self->{'queryNum'}}=$element->{'Attributes'}{'{}score'}{'Value'} unless $maxQueryScore{$self->{'queryNum'}};# macth are organized by decreasing scores
					$self->{'varMods'}=();
					$self->{'substMods'}=(); # for substitution

					$matches{$matchID}{'mz'}=$element->{'Attributes'}{'{}mz'}{'Value'};#m-obs
					$matches{$matchID}{'delta'}=$element->{'Attributes'}{'{}da_delta'}{'Value'};
					$matches{$matchID}{'sequence'}=$element->{'Attributes'}{'{}seq'}{'Value'};
					$matches{$matchID}{'score'}=$element->{'Attributes'}{'{}score'}{'Value'}; # Score of the peptide
					$matches{$matchID}{'confidence'}=$element->{'Attributes'}{'{}confidence'}{'Value'}*100.0; # Confidence: between 0 to 100 -> selection based on confidence
					$matches{$matchID}{'sequenceVmod'}=$element->{'Attributes'}{'{}ht'}{'Value'};
					$matches{$matchID}{'charge'}=$element->{'Attributes'}{'{}charge'}{'Value'};

					#$charge{$self->{'queryNum'}}=$charge{"-$self->{'queryNum'}"}=$element->{'Attributes'}{'{}charge'}{'Value'};
					#$matches{$matchID}{'mexp'}=($matches{$matchID}{'mz'}-1.007825032)*$matches{$matchID}{'charge'}; # M(Exp) -> Computed like this (Mobs - Masshydrogen )*Charge

					####> Update maxQueryScore with the first interpretation
					#if ($self->{'rankNum'}==1) {
					#	# Use of confidence instead of score which is not reliable...
					#	$maxQueryScore{$self->{'queryNum'}}=$matches{$matchID}{'confidence'};
					#	$maxQueryScore{"-$self->{'queryNum'}"}=$matches{$matchID}{'confidence'} if $self->{'isDecoySearch'};
					#	#$maxQueryScore{$self->{'queryNum'}}=$element->{'Attributes'}{'{}score'}{'Value'};
					#	#$maxQueryScore{"-$self->{'queryNum'}"}=$element->{'Attributes'}{'{}score'}{'Value'} if $self->{'isDecoySearch'};
					#}
				}
                        }elsif (($element->{'Name'} eq "MOD_FEATURE" || $element->{'Name'} eq "TERM_MOD_FEATURE") && $self->{'isOnMatch'}){# PTMs information
                                my $matchID=$self->{'matchID'};
                                my ($pos,$name)=($element->{'Attributes'}{'{}pos'}{'Value'},$element->{'Attributes'}{'{}mod'}{'Value'});
				#$self->{'varMods'}{$name}{$pos}=1;
				$self->{'varMods'}{$name}{$pos}=1 unless $name =~/^No /;
                        }elsif ($element->{'Name'} eq "SUBSTITUTION_FEATURE" && $self->{'isOnMatch'}) {# Substitution information
				my ($pos,$oneCodeLetter)=($element->{'Attributes'}{'{}pos'}{'Value'},$element->{'Attributes'}{'{}sub'}{'Value'});
				$self->{'substMods'}{$pos}{$oneCodeLetter}=1;
			}
#			elsif ($element->{'Name'} eq "TERM_MOD_FEATURE"){# Terminal PTMs information
#                                my $matchID=$self->{'matchID'};
#				my ($pos,$name)=($element->{'Attributes'}{'{}pos'}{'Value'},$element->{'Attributes'}{'{}mod'}{'Value'});
#				$self->{'varMods'}{$name}{$pos}=1;
#                        }elsif ($element->{'Name'} eq "SEARCH_STAT"){# Search parameters
#                                $anaInfo{'db_num_proteins'}=$element->{'Attributes'}{'{}num_proteins'}{'Value'};
#                                $anaInfo{'db_num_residues'}=$element->{'Attributes'}{'{}num_residues'}{'Value'};
#                        }elsif ($element->{'Name'} eq "FASTA"){# FASTA info
#                                $anaInfo{'db_filename'}=$element->{'Attributes'}{'{}filename'}{'Value'};
#                        }
			elsif ($element->{'Name'} eq "SEARCH" ){
				$self->{'isOnSearch'}=1;
				$self->{'searchID'}=$element->{'Attributes'}{'{http://www.w3.org/XML/1998/namespace}id'}{'Value'};# SEARCH:1:634817454795937500_0 : more than one search for same ID
                        }elsif ($self->{'isOnSearch'} && $element->{'Name'} eq "PARAMETERS" ){
				$search{$self->{'searchID'}}{'DECOY'}=($element->{'Attributes'}{'{}do_reverse_search'}{'Value'} eq 'true') ? 1 : 0;
				$self->{'isDecoySearch'}=$search{$self->{'searchID'}}{'DECOY'};
                        }
#			elsif ($element->{'Name'} eq "MSTOLERANCE" ){
#                                $anaInfo{'mstol_type'}=$element->{'Attributes'}{'{}type'}{'Value'};
#                                $anaInfo{'mstol_value'}=$element->{'Attributes'}{'{}value'}{'Value'};
#                        }elsif ($element->{'Name'} eq "MSMSTOLERANCE" ){
#                                $anaInfo{'msmstol_type'}=$element->{'Attributes'}{'{}type'}{'Value'};
#                                $anaInfo{'msmstol_value'}=$element->{'Attributes'}{'{}value'}{'Value'};
#                        }
			elsif ($element->{'Name'} eq "SEQUENCE" ){
                                # ex: <SEQUENCE xml:id="AADEAFAK" >
                                $self->{'sequenceID'}=$element->{'Attributes'}{'{http://www.w3.org/XML/1998/namespace}id'}{'Value'};
                                $self->{'isOnSequence'}=1;
				@{$self->{'matchIDs'}}=();
				@{$self->{'protContext'}}=();
				$self->{'tag'}++;
				if ($self->{'tag'} == 1000) {
					print ".";
					$self->{'tag'}=0;
				}
                        }
			elsif ($element->{'Name'} eq "PEPTIDE" && $self->{'isOnSequence'}==1){
                                # ex: <PEPTIDE matches="20571,20625,31229,47039,70067,87296" xml:id="PEPTIDE:9" />
				#$element->{'Attributes'}{'{http://www.w3.org/XML/1998/namespace}id'}{'Value'};
				push @{$self->{'matchIDs'}},$element->{'Attributes'}{'{}matches'}{'Value'}; # Important to group proteins!
#                                my @pepIDs=split(/,/,$element->{'Attributes'}{'{}matches'}{'Value'});
#				foreach my $pepID (@pepIDs) {
#					$peptides2matches{$self->{'sequenceID'}}=$pepID;# foreach peptide sequence in AA -> all the queries that match this sequence
#					$matches2peptides{$pepID}=$self->{'sequenceID'};# foreach pepID -> sequence in AA available in this hash
#				}
            }
			elsif ($element->{'Name'} eq "PROTEIN_CONTEXT" && $self->{'isOnSequence'}==1){
                                # ex: <PROTEIN_CONTEXT after="K" before="K" digest_prob="1" pos="530" protein="YDL229W" />
				my $entry=$element->{'Attributes'}{'{}protein'}{'Value'};
				my $decoyTag="";
				if ($entry=~/^RRRRR/) {
					$decoyTag="DECOY_";
					$entry=substr($entry, 5, length($entry));
				}
				my ($identifier)=($entry=~/$self->{'idRule'}/);
				$identifier="$decoyTag$identifier";
				my $before=($element->{'Attributes'}{'{}before'}{'Value'})? $element->{'Attributes'}{'{}before'}{'Value'} : '-';
				my $after=($element->{'Attributes'}{'{}after'}{'Value'})? $element->{'Attributes'}{'{}after'}{'Value'} : '-' ;
				my $pos=$element->{'Attributes'}{'{}pos'}{'Value'}+1;
				push @{$self->{'protContext'}}, "$pos,$before,$after:$identifier"; # pos,before,after:identifier -> ex 128,A,A:gi|229463044 for SYLDKVR peptide
            }
			elsif ($element->{'Name'} eq "PROTEIN") {#Protein information: score, sequence, name, description
                my $entry=$element->{'Attributes'}{'{http://www.w3.org/XML/1998/namespace}id'}{'Value'};
				my ($decoyTag,$isDecoyProtein)=('',0);
				if ($entry=~/^RRRRR/) {
					$decoyTag="DECOY_";
					$isDecoyProtein=1;
					$entry=substr($entry, 5, length($entry));
				}
				my ($identifier)=($entry=~/$self->{'idRule'}/);
				$self->{'identifier'}="$decoyTag$identifier";
				$self->{'description'}=$element->{'Attributes'}{'{}name'}{'Value'};
				$self->{'protScore'}=$element->{'Attributes'}{'{}protscore'}{'Value'};
				$self->{'protSequence'}=$element->{'Attributes'}{'{}sequence'}{'Value'};
				$self->{'protOrganism'}=$element->{'Attributes'}{'{}species'}{'Value'};
				$self->{'isOnProtein'}=1;
				#$proteins{$identifier};
				$protDes{$self->{'identifier'}}=$element->{'Attributes'}{'{}name'}{'Value'};
				$protOrg{$self->{'identifier'}}=$element->{'Attributes'}{'{}species'}{'Value'};
				$protLength{$self->{'identifier'}}=length($element->{'Attributes'}{'{}sequence'}{'Value'});
				$protDbRank{$self->{'identifier'}}=1; # One DB search
				my $mass=&computeAAMass($element->{'Attributes'}{'{}sequence'}{'Value'},'protein');
				$mass=sprintf "%.2f",$mass; # no need for more precision
				$protMW{$self->{'identifier'}}=$mass;
				###> To know if it is a reverse protein (ex: RRRRRgi|12585280) for DECOY computation
				if ($isDecoyProtein) {
					$self->{'isDecoyProtein'}=1;
					$self->{'isTargetProtein'}=0;
				}
				else{
					$self->{'isDecoyProtein'}=0;
					$self->{'isTargetProtein'}=1;
				}
            }
			elsif ($element->{'Name'} eq "PROTEIN_QUANT_RATIO" && $self->{'isOnProtein'}==1) {
				#<PROTEIN_QUANT_RATIO d="IT114" e="1.488459" lci="-1" n="IT115" p="1.763184e-014" r="0.02905569" uci="-1" />
				#$element->{'Attributes'}{'{http://www.w3.org/XML/1998/namespace}id'}{'Value'};
				#$element->{'Attributes'}{'{}d'}{'Value'};#Denominator
				#$element->{'Attributes'}{'{}e'}{'Value'};
				#$element->{'Attributes'}{'{}lci'}{'Value'};
				#$element->{'Attributes'}{'{}n'}{'Value'};#Numerator
				#$element->{'Attributes'}{'{}p'}{'Value'};#Probability?
				#$element->{'Attributes'}{'{}r'}{'Value'};#Ratio -> show in ProteinPilot interface, depending on which ratio the user want...
				#$element->{'Attributes'}{'{}uci'}{'Value'};
				#$element->{'Attributes'}{'{}'}{'Value'};
            }
			elsif ($element->{'Name'} eq "PROTEIN_QUANT_RAW_RATIO" && $self->{'isOnProtein'}) {
				#<PROTEIN_QUANT_RAW_RATIO d="IT114" e="1.488459" lci="-1" n="IT115" p="1.763184e-014" r="0.02905569" uci="-1" />
				#$element->{'Attributes'}{'{}'}{'Value'};
			}
			elsif ($element->{'Name'} eq "PROTEIN2MATCH" && $self->{'isOnProtein'}) {# Section for all the matching peptides against protein sequences...
					#<PROTEIN2MATCH background="0" digest_probability="1" match_id="68375" unused_protscore="0.99" use_quant="true" use_type="auto " />
					#use_type:  "auto " / "auto - low confidence" / "auto - discordant peptide type" / "no quant - weak signal"
					#$element->{'Attributes'}{'{}'}{'Value'};
				my $matchID=$element->{'Attributes'}{'{}match_id'}{'Value'};
				if ($matches{$matchID}) {# Maybe the score was too weak compared to minScore and this match was not recorded...
					my $actualSequence=$matches{$matchID}{'actualSequence'};

					$matchList{$self->{'identifier'}}{$actualSequence}=1;
					# Moved to another section
					#$maxProtScore{$self->{'identifier'}}+=$bestScore{$actualSequence};
					#$maxProtMatch{$self->{'identifier'}}{$actualSequence}++;
					$matches2proteins{$matchID}=1;
					if ($element->{'Attributes'}{'{}cleavage_features'}{'Value'}){# miss-cleavages information
						my @missSites=split(/;/,$element->{'Attributes'}{'{}cleavage_features'}{'Value'});
						$matches{$matchID}{'MIS'}=$#missSites+1; #ex: cleavage_features="missed K-G@4; missed R-S@13"
					}else{
						$matches{$matchID}{'MIS'}=0;
					}
					$matches{$matchID}{'matchDecoyProtein'}=1 if $self->{'isDecoyProtein'};
					$matches{$matchID}{'matchTargetProtein'}=1 if $self->{'isTargetProtein'};
				}
			}
			elsif ($element->{'Name'} eq "COVERAGE" && $self->{'isOnProtein'}==1) {# 3 type of coverage: with threshold value changing (0 / 0.5 / 0.95)
				#<COVERAGE coverage="0.1557" threshold="0" unique_peptides="7" />
				$self->{'coverage'}{$element->{'Attributes'}{'{}threshold'}{'Value'}}="$element->{'Attributes'}{'{}coverage'}{'Value'}:$element->{'Attributes'}{'{}unique_peptides'}{'Value'}";
			}
			elsif ($element->{'Name'} eq "MSMSPEAKS") {
				$self->{'isOnMSMS'}=1;
				$self->{'peaks'}="";
				# Name of the columns for the following peak-list: CENTROID,PEAK HEIGHT,PEAK AREA,CLUSTER AREA,CLUSTER AREA ERROR
				$self->{'peaks_attributes'}=$element->{'Attributes'}{'{}attributes'}{'Value'};
				# Number of lines for the following peak-list
				$self->{'peaks_size'}=$element->{'Attributes'}{'{}size'}{'Value'};
			}
			elsif ($element->{'Name'} eq "MSPEAKS") {
				$self->{'isOnMS'}=1;
				$self->{'peaks'}="";
				$self->{'peaks_attributes'}=$element->{'Attributes'}{'{}attributes'}{'Value'};
				$self->{'peaks_size'}=$element->{'Attributes'}{'{}size'}{'Value'};
			}
			elsif ($element->{'Name'} eq "ITRAQPEAKS") {
				$self->{'isOniTRAQ'}=1;
				$self->{'peaks'}="";
				$self->{'peaks_attributes'}=$element->{'Attributes'}{'{}attributes'}{'Value'};
				$self->{'peaks_size'}=$element->{'Attributes'}{'{}size'}{'Value'};
			#}elsif ($element->{'Name'} eq "PARAGON_METHOD_SETTINGS" && !$repeatedSection{$element->{'Name'}}) {
			###> TODO: for multiple search in PARAGON files.
			}
			elsif ($element->{'Name'} eq "PARAGON_METHOD_SETTINGS") {
                                $self->{'isOnPARAGONSettings'}=1;
                        }elsif ($self->{'isOnPARAGONSettings'}){
                                $self->{'paragagonSettingsName'}=$element->{'Name'};
                                $anaInfo{'PARAGONSETTINGS'}{$self->{'paragagonSettingsName'}}='';
                        }elsif ($self->{'isOnSearch'} && $element->{'Name'} eq "PEAKLIST" ) {
                                #<PEAKLIST count="4105" filenumber="8" original_count="4105" originalfilename="C:\Documents and Settings\Administrator\Desktop\QStar\E4208FD.wiff" partnumber="8" />
                                my @filename=split(/\\/,$element->{'Attributes'}{'{}originalfilename'}{'Value'});
				$search{$self->{'searchID'}}{'RAW_FILE'}=$filename[-1];
                        }elsif ($element->{'Name'} eq "DataDictionary") {
				$self->{'isOnDataDictionary'}=1;
			}
			elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "El"){
				$self->{'currentVmod'}=();
				$self->{'PryVal'}=0.00;
				$self->{'MassVal'}=0.00;
				$self->{'isOnEl'}=1;
			}
			elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "Nme"){# Name of an element, a modification, etc...
				$self->{'isOnNme'}=1;
				$self->{'Nme'}="";
			}
			elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "Sym"){# Symbol of the element
				$self->{'isOnSym'}=1;
				$self->{'Sym'}="";
			}
			elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "Iso"){# Define the isotope section
				$self->{'isOnIso'}=1;
			}
			elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "Mss"){# Mass of one of the isotope
				$self->{'isOnMss'}=1;
				$self->{'Mss'}="";
			}
			elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "Pry"){# Abundance of one of the isotope
				# In ESI Mass Spectrometry, we are focalised on the most abundant isotope
				# This is why in this XML-parsing, the most abundant Isotope is kept for mass computation
				$self->{'isOnPry'}=1;
				$self->{'Pry'}="";
			}
			elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "Mod"){
				$self->{'isOnMod'}=1;
                                $self->{'Fma'}='';
                                $self->{'RpF'}='';
                                $self->{'DisplayName'}='';
			}
			elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "Fma"){# Chemical formula of modification
				$self->{'isOnFma'}=1;
			}
			elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "RpF"){# Function that is replaced by the modification (could be H or more complex)
				$self->{'isOnRpF'}=1;
			}
			elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "DisplayName"){# Name of the modification that is displayed in ProteinPilot
				$self->{'isOnDisplayName'}=1;
			}

                        ###> Some redundant parts -> Do not read it again...
                        if (!$repeatedSection{$element->{'Name'}}){
                                $repeatedSection{$element->{'Name'}}=1;
                        }
                }

                sub end_element {
                        my ($self, $element) = @_;
                        if($element->{'Name'} eq "SPECTRUM"){#Keep order 1 spectrum to work only on MS information
                                $self->{'isOnSpectrum'}=0;
                        }elsif($element->{'Name'} eq "FEATURE_PROBABILITIES_SET" ){
				$self->{'isOnFeatures'}=0;
                        }elsif($self->{'isOnFeatures'} && $element->{'Name'} eq "MOD_FEATURE_SET" ){
				$self->{'isOnModFeature'}=0;
                        }elsif ($element->{'Name'} eq "SEQUENCE" ){
                                $self->{'isOnSequence'}=0;
				my $matchIDs=join(',',@{$self->{'matchIDs'}});
				foreach my $protContext (@{$self->{'protContext'}}) {
					my ($beg,$entry)=split(/:/,$protContext);
					foreach my $match (split(/,/,$matchIDs)) {
						$infoRank{$match}{$entry}{$beg}=1;
					}
				}
                        }elsif ($element->{'Name'} eq "PROTEIN") {
                                $self->{'isOnProtein'}=0;
                                #$proteins{$self->{'identifier'}}{'PROT_DESC'}=$self->{'description'};
                                #$proteins{$self->{'identifier'}}{'PROT_SCORE'}=$self->{'protScore'};
                                #$proteins{$self->{'identifier'}}{'PROT_SEQ'}=$self->{'protSequence'};
                        }elsif ($element->{'Name'} eq "SEARCH" ){
				$self->{'isOnSearch'}=0;
                        }elsif ($element->{'Name'} eq "MSMSPEAKS") {
                                #print $self->{'peaks'},"\n";
                                $self->{'isOnMSMS'}=0;
                        }elsif ($element->{'Name'} eq "MSPEAKS") {
                                $self->{'isOnMS'}=0;
                                #print $self->{'peaks'},"\n";
                        }elsif ($element->{'Name'} eq "ITRAQPEAKS") {
                                $self->{'isOniTRAQ'}=0;
                                #print $self->{'peaks'},"\n";
                        }elsif ($self->{'isOnSpectrum'} && $self->{'isOnMatch'} && $element->{'Name'} eq "MATCH") { # Do not consider the match for spectrum with charge=0
                                $self->{'isOnMatch'}=0;
				my ($varModString,$subString)=('','');
				my $matchID=$self->{'matchID'};
				if ($matches{$matchID}) {
					if ( $self->{'varMods'} ){
						foreach my $currentVarMod (sort{$a cmp $b} keys %{$self->{'varMods'}}) {
							if (!defined($self->{'anaVmods'}->{'V'}{$currentVarMod})) {
								$self->{'anaVmods'}->{'V'}{$currentVarMod}{'specificity'}=join(',',keys %{$specificity{$currentVarMod}});
								$self->{'anaVmods'}->{'V'}{$currentVarMod}{'numModDB'}=-1;
							}
							my $currentPosition=0;
							my @pos=(sort{$a<=>$b} keys %{$self->{'varMods'}{$currentVarMod}});
							my %aaVmod;# With PARAGON, same vmod could occur on multiple AA (ex: Oxidation -> M, P, R, etc...)
							my @sequence=split(//,$matches{$matchID}{'sequence'});
							foreach my $posRes (@pos){
								$aaVmod{$sequence[$posRes-1]}=1;
								#if (! $aaVmod{$sequence[$posRes-1]} =~ /$sequence[$posRes-1]/ ) {
								#	$self->{'anaVmods'}->{'V'}{$currentVarMod}{'specificity'}.=$sequence[$posRes-1];
								#}
							}
							my @aamod=(sort{$a cmp $b} keys %aaVmod);
							my $posString=join('', @aamod).':'.join('.', @pos);
							$currentVarMod.=" ($posString)";
							$varModString.=" + $currentVarMod";
						}
					}
					if ( $self->{'substMods'} ){ # AA substitution
						foreach my $currentPosSubst (sort{$self->{'substMods'}{$a}<=>$self->{'substMods'}{$b}} keys %{$self->{'substMods'}}) {
							my @residues=(sort{$a cmp $b} keys %{$self->{'substMods'}{$currentPosSubst}});
							my @sequence=split(//,$matches{$matchID}{'sequence'});
							foreach my $res (@residues){
								my $posRes=$currentPosSubst+1;
								$subString.=" + $sequence[$currentPosSubst]->$res ($posRes)";
							}
						}
					}
					#my $varModDataStrg=($varModString)? ",VMOD=$varModString":'';
					my $actualSequence="$matches{$matchID}{'sequence'}$varModString";
					$matches{$matchID}{'actualSequence'}=$actualSequence;

					#print "mz=$matches{$matchID}{'mz'} delta=$matches{$matchID}{'delta'} seq=$matches{$matchID}{'sequence'} score=$matches{$matchID}{'score'} seqVmod=$matches{$matchID}{'sequenceVmod'} charge=$matches{$matchID}{'charge'} confidence=$matches{$matchID}{'confidence'}\n";
					#print "$matchID '$actualSequence'<BR>\n";
					#my $select;
					# Change because score is not reliable for Paragon searches...
					#if(!$bestScore{$actualSequence}){
					#        $bestScore{$actualSequence}=$matches{$matchID}{'score'};
					#}elsif($bestScore{$actualSequence} < $matches{$matchID}{'score'}){
					#        $bestScore{$actualSequence}=$matches{$matchID}{'score'};
					#}
					####> Move on 09/03/2015 because bestScore should be filled considering maxrank : it is not known at this step if the Macth is for a decoy sequence or for a not decoy one
					#if(!$bestScore{$actualSequence} || ($bestScore{$actualSequence}[0] < $matches{$matchID}{'confidence'}) ){
					#	@{$bestScore{$actualSequence}}=($matches{$matchID}{'confidence'},$matchID,$self->{'rankNum'});
					#}
					#print "QUERY_NUM=$self->{'queryNum'} ACTUAL_SEQUENCE=$actualSequence MATCH_ID=$matchID mz=$matches{$matchID}{'mz'} delta=$matches{$matchID}{'delta'} seq=$matches{$matchID}{'sequence'} score=$matches{$matchID}{'score'} seqVmod=$matches{$matchID}{'sequenceVmod'} charge=$matches{$matchID}{'charge'} confidence=$matches{$matchID}{'confidence'} BS=$bestScore{$actualSequence}[0] BS=$bestScore{$actualSequence}[1] <BR>\n" if $matches{$matchID}{'sequence'} eq 'VLALAGFDVEK';
					#$numValid{$self->{'queryNum'}}++ unless $select==-1;
					#push @{$seqQueryRank{$matches{$matchID}{'sequence'}}},"$self->{'queryNum'}:$self->{'rankNum'}"; # $rankID
					#$query2matches{$self->{'queryNum'}}{$matches{$self->{'rankNum'}}}=$matchID;
					#$query2matches{$self->{'queryNum'}}{$self->{'rankNum'}}=$matchID; # Moved above
					$matches{$matchID}{'vmod'}=$varModString;
					$matches{$matchID}{'sub'}=$subString;
				}
                        }elsif ($element->{'Name'} eq "PARAGON_METHOD_SETTINGS") {
                                $self->{'isOnPARAGONSettings'}=0;
                        }elsif ($element->{'Name'} eq "RESULTS") {
				###> 1st loop to compute best score
				foreach my $queryNum (sort{$a<=>$b} keys %query2matches){
					my ($rankTarget,$rankDecoy,$rankMatch)=(0,0,0);
					foreach my $rank (sort{$a<=>$b} keys %{$query2matches{$queryNum}}) {
						my $matchID=$query2matches{$queryNum}{$rank};
						next unless $matches{$matchID};
						next unless $matches2proteins{$matchID};# Means that there is no match between this interpretation and the database
						if ($self->{'isDecoySearch'} && $matches{$matchID}{'matchDecoyProtein'}) {
							$rankDecoy++;
							if($rankDecoy <= $maxRank && (!$bestScoreDecoy{$matches{$matchID}{'actualSequence'}} || ($bestScoreDecoy{$matches{$matchID}{'actualSequence'}}[0] < $matches{$matchID}{'confidence'}) ) ){
								@{$bestScoreDecoy{$matches{$matchID}{'actualSequence'}}}=($matches{$matchID}{'confidence'},$matchID);
							}
						}
						if ( $matches{$matchID}{'matchTargetProtein'}) {
							$rankTarget++;
							if($rankTarget <= $maxRank && (!$bestScoreTarget{$matches{$matchID}{'actualSequence'}} || ($bestScoreTarget{$matches{$matchID}{'actualSequence'}}[0] < $matches{$matchID}{'confidence'}) ) ){
								@{$bestScoreTarget{$matches{$matchID}{'actualSequence'}}}=($matches{$matchID}{'confidence'},$matchID,$matches{$matchID}{'score'});
							}
						}
						last if ( ($rankTarget == $maxRank) && ($rankDecoy == $maxRank));# Avoid to get too many interpretations in INFO_PEP string
					}
				}
				###> 2nd loop to compute best score
				foreach my $queryNum (sort{$a<=>$b} keys %query2matches){
					my ($rankTarget,$rankDecoy,$rankMatch)=(0,0,0);
					#foreach my $rank (sort{$query2matches{$queryNum}{$a}<=>$query2matches{$queryNum}{$b}} keys %{$query2matches{$queryNum}}) {
					foreach my $rank (sort{$a<=>$b} keys %{$query2matches{$queryNum}}) {
						#if($rank==1){
						#	$queryInfo{$queryNum}=();
						#	$queryInfo{"-$queryNum"}=();
						#}
						my $matchID=$query2matches{$queryNum}{$rank};
						next unless $matches{$matchID};
						next unless $matches2proteins{$matchID};# Means that there is no match between this interpretation and the database
						next unless ($bestScoreTarget{$matches{$matchID}{'actualSequence'}} || $bestScoreDecoy{$matches{$matchID}{'actualSequence'}});
						my $vmod=$matches{$matchID}{'vmod'};
						$vmod ='' unless $vmod;
						my $actualSequence="$matches{$matchID}{'sequence'}$vmod";
						$vmod = ",VMOD=$vmod" if $vmod;
						my $sub=$matches{$matchID}{'sub'};
						$sub ='' unless $sub;
						$sub = ",SUBST=$sub" if $sub;

						my $mass=&computeAAMass($matches{$matchID}{'sequence'},'peptide');

                                                #print "Seq=$actualSequence BS=$bestScore{$actualSequence} MatchID=$matchID MatchScore=$matches{$matchID}{'score'}\n";
						#my $select=($matches{$matchID}{'score'}==$bestScore{$actualSequence})? 0:-1;
						#my $calc=scalar (keys %{$matches2proteins{$matchID}});
						my $usedRankID="";
						my $numProt=scalar (keys %{$infoRank{$matchID}});
						my $spc=scalar (keys %{$sequence2spectrum{"$matches{$matchID}{'sequence'}:$matches{$matchID}{'charge'}"}} );# Spectral-Count with no criterium of score -> all spectrum that match a sequence
						my $currentRank;
						if ($self->{'isDecoySearch'} && $matches{$matchID}{'matchDecoyProtein'}) {
							my $select=($bestScoreDecoy{$actualSequence}[1]==$matchID)? 0:-1;
							$queryInfo{"-$queryNum"}=() unless $queryInfo{"-$queryNum"};
							$rankDecoy++;
							$usedRankID="-$queryNum:$rankDecoy";
							$maxQueryScore{"-$queryNum"}=$matches{$matchID}{'confidence'} unless defined($maxQueryScore{"-$queryNum"});
							#push @{$queryInfo{"-$queryNum"}},"SEL=$select,SRK=$rank,MIS=$matches{$matchID}{'MIS'},CALC=$mass,DELT=$matches{$matchID}{'delta'},SEQ=$matches{$matchID}{'sequence'},SC=$matches{$matchID}{'score'}$vmod,MATCH=$numProt" if ($rankBad <= $maxRank);
							push @{$queryInfo{"-$queryNum"}},"SEL=$select,SRK=$rank,MIS=$matches{$matchID}{'MIS'},CALC=$mass,DELT=$matches{$matchID}{'delta'},SEQ=$matches{$matchID}{'sequence'},SC=$matches{$matchID}{'confidence'}$vmod$sub,MATCH=$numProt,SPC=$spc," if ($rankDecoy <= $maxRank);
							$numValid{"-$queryNum"}++ if (!$numValid{"-$queryNum"} && $select==0);
							$currentRank=$rankDecoy;
						}
						if ( $matches{$matchID}{'matchTargetProtein'}) {
							my $select=($bestScoreTarget{$actualSequence}[1]==$matchID)? 0:-1;
							$queryInfo{$queryNum}=() unless $queryInfo{$queryNum};
							$rankTarget++;
							$usedRankID="$queryNum:$rankTarget";
							$maxQueryScore{$queryNum}=$matches{$matchID}{'confidence'} unless defined($maxQueryScore{$queryNum});
							#push @{$queryInfo{$queryNum}},"SEL=$select,SRK=$rank,MIS=$matches{$matchID}{'MIS'},CALC=$mass,DELT=$matches{$matchID}{'delta'},SEQ=$matches{$matchID}{'sequence'},SC=$matches{$matchID}{'score'}$vmod,MATCH=$numProt" if ($rankGood <= $maxRank);
							push @{$queryInfo{$queryNum}},"SEL=$select,SRK=$rank,MIS=$matches{$matchID}{'MIS'},CALC=$mass,DELT=$matches{$matchID}{'delta'},SEQ=$matches{$matchID}{'sequence'},SC=$matches{$matchID}{'confidence'}$vmod$sub,MATCH=$numProt,SPC=$spc," if ($rankTarget <= $maxRank);
							$numValid{$queryNum}++ if (!$numValid{$queryNum} && $select==0);
							$currentRank=$rankTarget;
						}

						###> Last Step -> Fill rankProteinMatch information that is now available since the rank of the query is known
						if ($currentRank <= $maxRank) {
							foreach my $identifier (keys %{$infoRank{$matchID}}) {
								next unless $matchList{$identifier};
								$rankProtMatch{$usedRankID}{$identifier}=();
								$maxProtScore{$identifier}+=$bestScoreTarget{$actualSequence}[2] if $bestScoreTarget{$actualSequence}[2]; # Add individually score to compute a global protein score (do not consider confidence here)
								foreach my $beg (keys %{$infoRank{$matchID}{$identifier}}){
									push @{$rankProtMatch{$usedRankID}{$identifier}} , $beg;
									$maxProtMatch{$identifier}{$actualSequence}++;
								}
							}
						}
						last if ( ($rankTarget == $maxRank) && ($rankDecoy == $maxRank));# Avoid to get too many interpretations in INFO_PEP string
					}
					###> If no interpretation was added to that query, then, it does not have to be add in myProMS
					if (!$queryInfo{$queryNum}) {
						delete($charge{$queryNum});
						delete($massObs{$queryNum});
						delete($massExp{$queryNum});
						delete($elutionTime{$queryNum});
						delete($maxQueryScore{$queryNum});
					}
					if ($self->{'isDecoySearch'} && !$queryInfo{"-$queryNum"}) {
						delete($charge{"-$queryNum"});
						delete($massObs{"-$queryNum"});
						delete($massExp{"-$queryNum"});
						delete($elutionTime{"-$queryNum"});
						delete($maxQueryScore{"-$queryNum"});
					}
				}
				foreach my $identifier (keys %matchList) {
					delete($matchList{$identifier}) unless $maxProtScore{$identifier};
				}
				###> Store all the VMODS seen in this analysis so as to prevent from reading everytime the XML file to find the correct VMOD
				my $anaFile = $self->{'msFilename'};
				$anaFile =~ s/\.xml/\.pgf/g;# PARAGON generic file -> Write VMODS in a specific file.
				open(DESCP, ">$anaFile") || die ("unable to create $anaFile <BR>\n") ;
				print DESCP "--gc0p4Jq0M2Yt08jU534c0p\n";
				print DESCP "Content-Type: Paragon; name=\"search\"\n";
				my (%rawfileNames,$currentrawfileName);
				foreach my $searchID (keys %search) {
					$currentrawfileName=$search{$searchID}{'RAW_FILE'};
					$rawfileNames{$currentrawfileName}=1;
					print DESCP "$searchID,$currentrawfileName\n";
				}
				###> Update anaInfo information
				$anaInfo{'WIFF_FILE'}=(scalar (keys (%rawfileNames) ) == 1)? $currentrawfileName : "Merged Research File";
				print DESCP "--gc0p4Jq0M2Yt08jU534c0p\n";
				print DESCP "Content-Type: Paragon; name=\"masses\"\n";
				my $delta=1;
				foreach my $anaVmod (sort{ $a cmp $b} keys %{$self->{'anaVmods'}->{'V'}}){
					print DESCP "delta$delta=$allVmods{$anaVmod},$anaVmod\n";
					$delta++;
				}
				print DESCP "--gc0p4Jq0M2Yt08jU534c0p\n";
				close DESCP;
			}elsif ($element->{'Name'} eq "DataDictionary") {
				$self->{'isOnDataDictionary'}=0;
			}elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "El"){# End of Element part
				$self->{'isOnEl'}=0;
				# Ex: $elements{'He'}=4.00260325 , $elements{'H'}=1.007825032 , ...
				$elements{$self->{'Sym'}}=$self->{'MassVal'};
			}elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "Nme"){# Complete Name of the element
				$self->{'isOnNme'}=0;
			}elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "Sym"){# Symbol of the element
				$self->{'isOnSym'}=0;
			}elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "Iso"){# Define the isotope section
				$self->{'isOnIso'}=0;
			}elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "Mss"){# Mass of one of the isotope
				$self->{'isOnMss'}=0;
			}elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "Pry"){# Abundance of one of the isotope
				# In ESI Mass Spectrometry, we are focalised on the most abundant isotope
				# This is why in this XML-parsing, the most abundant Isotope is kept for mass computation
				$self->{'isOnPry'}=0;
				if ( $self->{'Pry'} > $self->{'PryVal'}) { # Keep the most abundant isotope for the current element
					$self->{'PryVal'}=$self->{'Pry'};
					$self->{'MassVal'}=$self->{'Mss'};
				}
			}elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "Mod"){
				###> Compute the delta-mass of the current VMOD...
				$self->{'isOnMod'}=0;
				$allVmods{$self->{'Nme'}}=&computeChemMass($self->{'Fma'},\%elements) - &computeChemMass($self->{'RpF'},\%elements);
				if($self->{'DisplayName'}){ $allVmods{$self->{'DisplayName'}}=$allVmods{$self->{'Nme'}};}# For some modifications, it has to be the DisplayName otherwise, an error would occur: ex: Amidated (for Terminal Amidated)
			}elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "RpF"){
				$self->{'isOnRpF'}=0;
			}elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "Fma"){
				$self->{'isOnFma'}=0;
			}elsif($self->{'isOnDataDictionary'} && $element->{'Name'} eq "DisplayName"){
				$self->{'isOnDisplayName'}=0;
			}
                }

                sub characters {
                        my ($self, $element) = @_;

                        if ($self->{'isOnMSMS'} && $self->{'isOnMSMS'} == 1 ) {
                                $self->{'peaks'}.=$element->{'Data'};
                        }elsif ($self->{'isOnPARAGONSettings'} && $self->{'paragagonSettingsName'}) {
                                $anaInfo{'PARAGONSETTINGS'}{$self->{'paragagonSettingsName'}}.=$element->{'Data'};
                        }elsif($self->{'isOnDataDictionary'} && $self->{'isOnSym'}){# Symbol of the element
				$self->{'Sym'}.=$element->{'Data'};
			}elsif($self->{'isOnDataDictionary'} && $self->{'isOnMss'}){# Mass of one of the isotope
				$self->{'Mss'}=$element->{'Data'};
			}elsif($self->{'isOnDataDictionary'} && $self->{'isOnPry'}){# Abundance of one of the isotope
				# In ESI Mass Spectrometry, we are focalised on the most abundant isotope
				# This is why in this XML-parsing, the most abundant Isotope is kept for mass computation
				$self->{'Pry'}.=$element->{'Data'};
			}elsif($self->{'isOnDataDictionary'} && $self->{'isOnRpF'}){
				$self->{'RpF'}.=$element->{'Data'};
			}elsif($self->{'isOnDataDictionary'} && $self->{'isOnFma'}){
				$self->{'Fma'}.=$element->{'Data'};
			}elsif($self->{'isOnDataDictionary'} && $self->{'isOnNme'}){
				$self->{'Nme'}.=$element->{'Data'};
			}elsif($self->{'isOnDataDictionary'} && $self->{'isOnDisplayName'}){
				$self->{'DisplayName'}.=$element->{'Data'};
			}
                }

                sub start_cdata {
                        my ($self, $element) = @_;
                        #$self->{'peaks'}="";
                }

                sub end_cdata {
                        my ($self, $element) = @_;
                        #if ($self->{'peaks'} && $self->{'isOnMSMS'}==1) {
                        #    #print "CDATA -> $self->{'peaks'}\n";
                        #}elsif ($self->{'peaks'} && $self->{'isOnMSMS'}==1) {
                        #
                        #}elsif ($self->{'peaks'} && $self->{'isOnMS'}==1) {
                        #
                        #}elsif ($self->{'peaks'} && $self->{'ITRAQPEAKS'}==1) {
                        #
                        #}
                }

                sub get_spectra {
                        my ($self)=@_;
                        return \%spectra;
                }

                sub get_matches {
                        my ($self)=@_;
                        return \%matches;
                }

		###> Function that computes the mass of a chemical formula based on the masses of the most abundant isotopes
		###> of all the elements that are contained in a PARAGON file
		sub computeChemMass {# Ex: C15H11O2, C1Cb6H13Nc2O, C2H3O, etc.
			my ($formula,$refElements)=@_;
			my ($element,$nb,$mass)=('','',0.000000000);

			chomp($formula);
			if ($formula){
				foreach my $letter (split(//,$formula)){
					if ($letter =~ /[A-Z]/) {
						if ($element){
							$nb=1 unless $nb;
							$mass+=$refElements->{$element}*$nb;
							$element='';
							$nb='';
						}
						$element=$letter;
					}elsif($letter =~ /[a-z]/) {
						$element.=$letter;
					}elsif($letter =~ /\d/) {
						$nb.=$letter;
					}
				}
				###> Last element that was missed in the loop has to be added in the mass computation
				$nb=1 unless $nb;
				$mass+=$refElements->{$element}*$nb;
			}
			return $mass;
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

		# Computation of Substitution : ex: A->V (A: origin / V: substitution -> returns 28.03130)
		sub computeSubstAA {
			my ($aaOrigin,$aaSubst)=@_;
			my %massAAave=&promsConfig::getMassAAave;
			return ($massAAave{$aaOrigin} - $massAAave{$aaSubst});
		}

	}




####>Revision history<####
# 3.7.2 Minor modif to rename promsTandemXML to promsDIA (MLP 06/12/17)
# 3.7.1 Minor bug correction (GA 30/11/17)
# 3.7.0 Update to PD 2_2 (GA 01/09/17)
# 3.6.5 Minor correction (MLP 03/08/17)
# 3.6.4 Minor modif (MLP 27/06/17)
# 3.6.3 Add symbolic link for mzxml file (in &applyFDR) (MLP 11/05/17)
# 3.6.2 Modification of Vmods parsing for double modification on same residue like Label:2H(4)+GG (GA 04/05/17)
# 3.6.1 Minor modif (MLP 04/05/17)
# 3.6.0 Added importation of X! Tandem files and Mayu option for data filtering (MLP 05/04/17)
# 3.5.9 Store substitutions from Mascot DAT files (GA 12/12/16)
# 3.5.8 Add Percolator for Mascot (GA 23/09/16)
# 3.5.7 Minor Change (PP 29/03/16)
# 3.5.6 Fix bug for modif specificity with msf file with quantif when labeling specificity is declared on multiple lines (PP 16/07/15)
# 3.5.5 Check SQlite queries for 2.0 PD MSF (GA 12/06/15)
# 3.5.4 Minor modification to get all decoy queries from Paragon XML (PP-GA 13/04/15)
# 3.5.3 Update BestScore computation for Paragon searches (GA 19/03/15)
# 3.5.2 Minor modification for Paragon searches for matches that count for decoy and target proteins at the same time (GA 19/03/15)
# 3.5.1 Labeling/Quantif detection for PARAGON.XML && bug fix in activation of some lower scoring peptides in PARAGON.XML<BR>& Modification for %bestScore of Paragon search to avoid to consider ranks not kept in the end (PP,GA 11/03/15)
# 3.5.0 Change for split-mode of msf files (GA 03/03/15)
# 3.4.9 Update ProteinPilot XML parsing:<BR>- to fill %protDbRank<BR>- to parse correctly db entries/identifiers with parsing rules<BR>- to correct a bug in %rankProteinMatch filling (GA 03/11/14)
# 3.4.8 Minor update for ProteinPilot XML (GA 21/10/14)
# 3.4.7 Comment line 278 which is an exit for dev purposes (GA 22/05/14)
# 3.4.6 Handles precomputed Percolator data in Proteome Discoverer (PP 12/05/14)
# 3.4.5 Minor improvement of scan (GA 15/04/14)
# 3.4.4 qvality command handles files with spaces & pdm always generated in user's directory (PP 15/04/14)
# 3.4.3 Improved parsing of scan/retention time (PP 26/03/14)
# 3.4.2 DT count FDR-based score estimation (PP 06/03/14)
# 3.4.1 Better file format error handling (PP 22/01/14)
# 3.4.0 MSF to PDM convertion moved to convertMsf2Pdm.cgi through IFRAME streaming (PP 17/01/14)
# 3.3.1 Minor pattern matching modification for $selMinScore (GA 27/11/13)
# 3.3.0 Replace most system commands with Perl commands (PP 05/11/13)
# 3.2.9 Removed 'next' command that aborted PMF-analyses import (FY 07/10/13)
# 3.2.8 Minor modification for decoy proteins printing if proteins are not kept in MSF (GA 04/09/13)
# 3.2.7 Decoy proteins ONLY for PD >= 1.4 (PP 28/08/13)
# 3.2.6 Adding decoy proteins in PDM files (GA 09/08/13)
# 3.2.5 Fitting with AUTO_INCREMENT for QUERY_VALIDATION and PROTEIN_VALIDATION primary keys (FY 25/07/13)
# 3.2.4 Merge of 3.2.3 & 3.2.2f (PP 24/07/13)
# 3.2.2f Converts custom label modifs name to true unimod name in pdm files (PP 17/07/13)
# 3.2.2d Improved label-free detection & labeling string no longer added to modification list (PP 11/07/13)
# 3.2.2c Fixed problem with qvality number rounding (PP 09/07/13)
# 3.2.2b Fixed typo when commenting &addSynonyms (PP 08/07/13)
# 3.2.2 More Proteome Discoverer version compatibility (PP 28/06/13)
# 3.2.1 Get misscleavage number from MSF PD 1.3 and 1.4 versions instead of computing it<BR>convertVarModString function is back again in promsMod (GA 26/06/13)
# 3.2.0 Minor modification to get scanNumber for all queries in DAT files (GA 14/06/13)
# 3.1.8 Minor modification in &getAllSearches to be consistent with new PD1.4 ProcessingNodes table<BR>+ modification in &convertVarModString<BR>+ remove too many calls of promsMod::getModificationIDfromString in peptide insertion (GA 27/05/13)
# 3.1.7 Bug correction, N-Term modifications were not kept and specificity for ANALYSIS_MODIFICATION were bad (GA 22/05/13)
# 3.1.6 unimod_tables.xml file handled by &promsMod::getModificationIDfromString (PP 21/05/13)
# 3.1.5 Getting specificity for modifications in Paragon searches to avoid empty statements in ANALYSIS_MODIFICATION +<BR>Bug correction in Paragon import +<BR>Minor modification in convertVarModString for residues formating (GA 21/05/13)
# 3.1.4 Fixing duplicate ANALYSIS_MODIFICATION entry in DB (FY 17/05/13)
# 3.1.3 Move locally &convertVarModString from promsMod and custom it for new PTMs handling (GA 15/05/13)
# 3.1.2 Lower-scoring peptide flag & commented PROJECT_MODIFICATION management (PP 26/04/13)
# 3.1.1 Get quantification information in Dat files and MSF files (GA 23/04/13)
# 3.1.0 Add code for MODIFICATION tables (GA 11/04/13)
# 3.0.9 Substitution in peptides: for paragon searches (GA 14/03/13)
# 3.0.8 Merge 3.0.7PP and 3.0.4GA (GA 14/03/13)
# 3.0.7 Handles PD 1.4 & minor bugs correction (PP 11/03/13)
# 3.0.6 Re-introduced [NC]-term check on terminal modifs assignment to peptides **Revalidation not tested** (PP 10/01/13)
# 3.0.5 Revalidation update. Allowed for MASCOT.DAT with no external decoy file (PP 18/12/12)
# 3.0.4b Minor modification on misscleavages count for PDM (GA 22/02/13) -> to discuss!
# 3.0.4 Linear estimation of FDR between flanking value pairs<BR>Multi-databank searches<BR>Skip peptides with score=0 (PP 11/12/12)
# 3.0.3 Change test on $maxQueryScore in applyFDR and elsewhere <BR>-> defined($maxQueryScore{$queryNum}): in Mascot, some queries could have a score=0 (GA 26/11/12)
# 3.0.2 Modification in varmods match to prevent misleading varmods such as Methyl, Dimethyl, Trimethyl etc. (GA 14/11/12)
# 3.0.1 Modify the query to get varmods for Sequest and Mascot - ProteomeDiscoverer (GA 02/10/12)
# 3.0.0 ProteinPilot group files (converted into XML) import for PARAGON data (GA 02/10/12)
# 2.9.1m2 Qvality for mixed target/decoy databanks & fix bug min score (PP 24/09/12)
# 2.9.1m Disconnects to DB before forking for DB scan (PP 11/09/12)
# 2.9.1 minor debug -> $FDRsuccess is defined 2 times in applyFDR function (GA 23/08/12)
# 2.9.0 Added Qvality for min score estimation based on preset FDR<BR>& immediate background databank scan (PP 13/08/12)
# 2.8.2 Modification for MSF redundancy names in same projects<BR>Minor-modification in PDM creation -> for merge, better file description (GA 24/07/12)
# 2.8.1 Minor modification of MSF parsing to get more information about the Fasta-File (GA 13/06/12)
# 2.8.0 Fix bug regular modifs used instead of terminal ones (PP 21/05/12)
# 2.7.9 Symbolic link is preserved if validation file is already a link (PP 03/05/12)
# 2.7.8 Modification of ORDER options in $sthSpInfo and $sthSpec to make it consistent one to another and identical to MASCOT.DAT file<BR>Modification of printParametersMSF for Neutral-Losses in Phosphorylations (GA 03/05/12)
# 2.7.7 Minor change in getAllSearches function to make it works for both PD 1.2 and 1.3 version<BR>Modification of QUERY for $sthParam2 to not get redundant Modifications (GA 27/04/12)
# 2.7.6 chmod commented (PP 20/03/12)
# 2.7.5 Bug on labeling corrected (QUANTITATION=None for MALDI experiments) (GA 13/01/12)
# 2.7.4 Optimization of SQLight queries during peptide data extraction from MSF file (PP 10/01/12)
# 2.7.3 Added sort on IT_MODS list for pdm files (PP 09/11/11)
# 2.7.2 correct bug in SEQUEST modification catching -> Add a new key {'PRINTED_NAME'} in modification hash table (GA 17/10/11)
# 2.7.1 better handling of Modif in MSF especially SILAC eg. Label:13C(6) (PP 05/07/11)
# 2.7.0 no longer uses DATA_FILE FROM QUERY_VALIDATION & PROTEIN_VALIDATION tables (PP 14/06/11)
# 2.6.9 Starting SQ search management (PP 12/05/11) & auto-correction of IT_MODS: __mod -> mod (PP 08/06/11)
# 2.6.8 multiple bugs in PMF mode & standardization of varMods position (PMF same as MIS now) (PP 13/05/11)
# 2.6.7 query $sthSpInfo modified (avoid to get a charge equals to 0 from MSF files) (GA 06/01/11)
# 2.6.6 add spectral count information in INFO_PEP at the end of the string (ex: ...,SPC=12,)
# 2.6.5 modification of ELUTION_TIME format for QUERY_VALIDATION ex: et32.28;sp3002;sc2993 <BR>et: elution time (min.) <BR>sp: spectrumNumber <BR>sc: scanNumber (information used to match a MS in T3PQ) <BR>->Impact on scripts that use this string (listRanks.cgi, sequenceView.cgi)
