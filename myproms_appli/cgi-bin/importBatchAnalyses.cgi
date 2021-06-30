#!/usr/local/bin/perl -w

################################################################################
# importBatchAnalyses.cgi         2.8.5                                        #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
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
$| = 1;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use POSIX qw(strftime); # to get the time
use File::Copy qw(copy move);
use File::Spec::Functions qw(splitpath); # Core module
use LWP::UserAgent;
use promsConfig;
use promsMod;
use strict;
use XML::Simple; # used for MSF 2.0 version

#print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
mkdir $promsPath{tmp} unless -e $promsPath{'tmp'};
my $batchDir="$promsPath{tmp}/batch";
mkdir $batchDir unless -e $batchDir;
$CGITempFile::TMPDIRECTORY=$batchDir;
my $userID=$ENV{'REMOTE_USER'};
mkdir "$batchDir/$userID" unless -e "$batchDir/$userID";
my ($color1,$color2)=&promsConfig::getRowColors;
my $maxUpFiles=10;
my $maxRank=&promsConfig::getMaxRank;
my %mascotServers=&promsConfig::getMascotServers;

#############################
####>Fetching parameters<####
#############################
my $ID=param('ID');
my ($item,$itemID)=split(':',$ID);
my $action=param('ACT'); # start,preview,delete,clean,UseUserDirectory,UseServerDirectory,UseZipFile,UseUploadedFiles,UseSharedDirectory
my $decoyFile=(param('decoy'))? 1 : 0;

if ($action eq 'UseMascot') {

	####>Parameters<####
	my $mascotAction=param('mascotAction');
	my $mascotServer=param('mascotServer');

	#################################
	####>Searching mascot server<####
	#################################
	if ($mascotAction eq 'list') {

		####>Parameters<####
		my $startDate=param('startDate');
		my $endDate=param('endDate');
		my $startJob=param('startJob');
		my $endJob=param('endJob');
		my $params={ACT=>'list',
					LOG=>join(':',param('searchLogs')),
					IDS=>param('mascotIDs')
					};
		$params->{DR}="$startDate:$endDate" if $startDate;
		$params->{JR}="$startJob:$endJob" if $startJob;
		$params->{STRG}=param('searchStrg') if param('searchStrg');

		####>Starting HTML<####
		print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);
		print qq
|<HTML>
<HEAD>
<TITLE>Searching Files on Mascot Server</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE>
A.file,A.badFile{font-size:11px;}
A.badFile{color:#707070;}
A.file:hover,A.badFile:hover{color:#FF6600; text-decoration:none;}
</STYLE>
</HEAD>
<BODY style="margin:0px;">
<DIV id="waitDiv">
<CENTER><BR><BR>
<IMG src="$promsPath{images}/engrenage.gif"><BR>
<B>Searching '$mascotServer'...</B>
</CENTER>
|;
		#foreach (1..200) {print "<!-- Forcing perl to print its buffer -->\n";}
		print "</DIV>\n";
		my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
		$agent->timeout(360);
		if ($mascotServers{$mascotServer}{proxy}) { # proxy settings
			if ($mascotServers{$mascotServer}{proxy} eq 'no') {$agent->no_proxy($mascotServers{$mascotServer}{url});}
			else {$agent->proxy('http', $mascotServers{$mascotServer}{proxy});}
		}
		else {$agent->env_proxy;}
        #push @{$agent->requests_redirectable},'POST';
        my $response = $agent->post("$mascotServers{$mascotServer}{url}/cgi/myproms4datFiles.pl", $params);
		while (my $wait = $response->header('Retry-After')) {
			sleep $wait;
			$response = $agent->get($response->base);
        }
		my @resultLines;
        if ($response->is_success) {
			@resultLines = split("\n",$response->content);
		}
		####>Building file tree<####
		print qq
|<SCRIPT type="text/javascript">
var fileInfo={};
|;
		my %treeOptions;
		$treeOptions{'BUTTON'}='font-size:9px;width:55px';
		%{$treeOptions{'CHECKBOX'}}=('folder'=>1,'file'=>1);
		@{$treeOptions{'SELECTION'}}=('mascot',0);
		my @fileTree=(0,'mascot',0,'','',1,$mascotServer,'');
		my $prevDirectory='';
		my @directoryTree;
		my $okServer=0;
		my $deniedRequest=-1;
		my %datNumbers;
		foreach (@resultLines) {
			$okServer=1 if /^#/;
			if (!$deniedRequest && /^# REQUEST DENIED/) {$deniedRequest=1; last;}
			$deniedRequest=0;
			next if /^#/;
			chomp;
			my ($filePath,$userID,$okFlag,$title)=split(/\t/,$_);
			my ($parent,$dataDir,$currentDirectory,$datFile)=($filePath=~/\//)? split(/\//,$filePath) : split(/\\/,$filePath); # Window or Unix/Linux
			if ($currentDirectory ne $prevDirectory) {
				if ($prevDirectory) {
					my @tempDirTree=@directoryTree;
					push @fileTree,\@tempDirTree; # cannot use \@directoryTree because reference will be updated
				}
				@directoryTree=(1,'folder',$currentDirectory,'','',1,$currentDirectory,''); # reset to current
				$prevDirectory=$currentDirectory;
			}
			my $fileType=($okFlag)? 'file' : 'badFile';
			my $fileClass=$fileType;
			$title='' unless $title;
			$title=&promsMod::HTMLcompatible($title);
			$title=quotemeta($title);
			$title=~s/(.{32,42}) /$1<BR>&nbsp;&nbsp;&nbsp;&nbsp;/g;
			#print "fileInfo['$fileType:$filePath']='$datFile:#:$okFlag:#:",quotemeta($title),"';\n";
			print "fileInfo['$fileType:$filePath']='$datFile:#:$userID:#:$okFlag:#:$title';\n";
			push @directoryTree,[2,$fileType,$filePath,$fileClass,'',1,$datFile,''];
			my ($datNum)=($datFile=~/(\d+)/);
			$datNumbers{$datFile}=$datNum;
		}
		push @fileTree,\@directoryTree if scalar @directoryTree; # last directory in list
		if ($deniedRequest) {push @fileTree,[1,'deniedRequest',0,'','',0,'Request denied.',''];} # connection was not allowed
		elsif (scalar @fileTree==8) {push @fileTree,[1,'noMatch',0,'','',0,'No matching files found.',''];} # no matching files

		###>Reordering files in tree: current order is finish time NOT start time (=dat number)
		foreach my $refDir (@fileTree[8..$#fileTree]) {
			my @fileList=@{$refDir}[8..$#{$refDir}];
			my $idx=7;
			foreach my $refFile (sort{$datNumbers{$a->[6]}<=>$datNumbers{$b->[6]}} @fileList) {
				$refDir->[++$idx]=$refFile;
			}
		}

		####>Tree JS functions<####
		&promsMod::writeJsTreeFunctions(\@fileTree,\%treeOptions);
		print qq
|//--->Local functions<---//
function actionOnSelect(tabIndex) {
	var infoTxt='';
	if (fileInfo[tableArray[tabIndex][4]]) {
		infoTxt+='<BR><FIELDSET><LEGEND><B>File info:</B></LEGEND>';
		var data=fileInfo[tableArray[tabIndex][4]].split(/:#:/);
		infoTxt+='&nbsp;&nbsp;<B>-Name:</B> '+data[0]+'<BR>&nbsp;&nbsp;<B>-Status:</B> ';
		infoTxt+=(data[2]==0)? '<FONT style="font-weight:bold;color:#DD0000">Removed</FONT>' : '<FONT style="font-weight:bold;color:#00DD00">Available</FONT>';
		infoTxt+='<BR>&nbsp;&nbsp;-<B>Search Title:</B> ';
		infoTxt+=(data[3])? data[3] : '';
		infoTxt+='<BR>&nbsp;&nbsp;-<B>User ID:</B> '+data[1];
		infoTxt+='</FIELDSET><BR><BR><BR><BR>';
	}
	parent.document.getElementById('infoSpan').innerHTML=infoTxt;
}
function actionOnCheck(tabIndex) {}
</SCRIPT>
<DIV id="treeDiv" style="display:none">
|;
		####>Drawing tree<####
		my $tableIndex=0;
		&promsMod::printTree(\@fileTree,\%treeOptions,\$tableIndex,1); # (refArray, is last child)

		print qq
|</DIV>

<SCRIPT type="text/javascript">
click2OpenClose(0,'open',1);
document.getElementById('waitDiv').style.display='none';
document.getElementById('treeDiv').style.display='block';
parent.document.fileAccessForm.mascotAction.value='get'; // switch search to get mode
</SCRIPT>
</BODY>
</HTML>
|;
	}

	#############################
	####>Fetching .dat files<####
	#############################
	elsif ($mascotAction eq 'get') {

		####>Parameters<####
		my $batchFilesDir="$batchDir/$userID";
		#mkdir $batchFilesDir unless -e $batchFilesDir;
		mkdir "$batchFilesDir/.percolator" unless -e "$batchFilesDir/.percolator"; # to store percolatot *.pop files created in Mascot DAT searches
		my ($itemType,@datFiles)=split(/[:,]/,param('datFiles'));

		####>Starting HTML<####
		print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);
		print qq
|<HTML>
<HEAD>
<TITLE>Fetching Files Mascot Server</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY>
<FONT class="title3">Retrieving files...</FONT><BR>
|;
		####>Looping through selected files<####
		my %alreadyRetrieved;
		foreach my $filePath (@datFiles) {
			my ($file)=($filePath=~/(F\d+\.dat)/);
			(my $spanStrg,$alreadyRetrieved{$file})=(-e "$batchFilesDir/$file")? ("<IMG src='$promsPath{images}/good.gif'>",1) : ('',0); # file already retrieved
			print "<B>-$file<SPAN id=\"$file\">$spanStrg</SPAN></B><BR>\n";
		}
		foreach my $filePath (@datFiles) {
			my ($file)=($filePath=~/(F\d+\.dat)/);
			next if $alreadyRetrieved{$file};
			print "<SCRIPT type=\"text/javascript\">document.getElementById('$file').innerHTML=\"...\";</SCRIPT>\n";
			#sleep 2;
			my $newDatFile="$batchFilesDir/$file";
			my ($popTargetFile,$popDecoyFile)=("$batchFilesDir/.percolator/$file.target.pop","$batchFilesDir/.percolator/$file.decoy.pop");
			unless (-e $newDatFile) {
				if ($mascotServers{$mascotServer}{data_local_path}) { # Mascot is on same computer or NFS
					$mascotServers{$mascotServer}{data_local_path}=~s/\/\Z//; # remove end '/' if any
					$filePath=~s/.+\/(\d{8}\/.+)/$1/; # remove everything before date folder (../data/<date>/Fxxxxxx.dat)
					if ($mascotServers{$mascotServer}{link_files}) { # symbolic link (no copy)
						symlink("$mascotServers{$mascotServer}{data_local_path}/$filePath",$newDatFile);
					}
					else { # copy
						copy("$mascotServers{$mascotServer}{data_local_path}/$filePath",$newDatFile);
					}
					# Copy *.pop files for percolator
					my ($year,$month)=$filePath=~/^(\d{4})(\d{2})/; # <date>/Fxxxxxx.dat
					my $mainCacheDir="$mascotServers{$mascotServer}{data_local_path}/cache/$year/$month"; # "$mascotServers{$mascotServer}[2]/data/cache/$year/$month"
					opendir(CACHEDIR,$mainCacheDir);
					my $foundFile=0;
					while ((my $cacheDir=readdir(CACHEDIR))) {
						next if ($cacheDir =~ /^\.+$/ || ! -d "$mainCacheDir/$cacheDir"); # skip . and ..
						opendir(DATDIR,"$mainCacheDir/$cacheDir");
						while ((my $cacheFile=readdir(DATDIR))) {
							next if $cacheFile =~ /^\.+$/; # skip . and ..
							if ($cacheFile =~ /$file.*\.target\.pop/) {
								if ($mascotServers{$mascotServer}{link_files}) {
									symlink("$mainCacheDir/$cacheDir/$cacheFile",$popTargetFile);
								}
								else{
									copy("$mainCacheDir/$cacheDir/$cacheFile",$popTargetFile);
								}
								$foundFile++;
							}
							elsif($cacheFile =~ /$file.*\.decoy\.pop/){
								if ($mascotServers{$mascotServer}{link_files}) {
									symlink("$mainCacheDir/$cacheDir/$cacheFile",$popDecoyFile);
								}
								else{
									copy("$mainCacheDir/$cacheDir/$cacheFile",$popDecoyFile);
								}
								$foundFile++;
							}
							last if $foundFile==2;
						}
						close DATDIR;
						last if $foundFile==2;
					}
					close CACHEDIR;
				}
				else { # Mascot is on another computer -> http
					my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
					if ($mascotServers{$mascotServer}{proxy}) { # proxy settings
						if ($mascotServers{$mascotServer}{proxy} eq 'no') {$agent->no_proxy($mascotServers{$mascotServer}{url});}
						else {$agent->proxy('http', $mascotServers{$mascotServer}{proxy});}
					}
					else {$agent->env_proxy;}
					#$agent->mirror("$mascotServers{$mascotServer}[0]/cgi/myproms4datFiles.pl?ACT=get&DAT=$filePath",$newDatFile);
					my $okFile=1;
					open(DAT,">$newDatFile") or $okFile=0;
					if ($okFile) {
						my $before=time;
						my $prevProgress=0;
						my $fileLength;
						my $bytesReceived = 0;
						$agent->request(HTTP::Request->new(GET => "$mascotServers{$mascotServer}{url}/cgi/myproms4datFiles.pl?ACT=get&DAT=$filePath"),
							sub {
								my ($chunk, $res) = @_;
								print DAT $chunk;
								unless (defined $fileLength) {
									$fileLength = $res->content_length || 0;
								}
								if ($fileLength) {
									$bytesReceived += length($chunk);
									my $progress=100 * $bytesReceived / $fileLength;
									if ($progress-$prevProgress > 5) { # in %
										my $progressStrg=sprintf "%d%%",$progress;
										print "<SCRIPT type=\"text/javascript\">document.getElementById('$file').innerHTML=\"...$progressStrg\";</SCRIPT>\n";
										$prevProgress=$progress;
									}
								}
								else {
									my $now=time;
									if ($now-$before > 30) { # in seconds
										print "<SCRIPT type=\"text/javascript\">document.getElementById('$file').innerHTML+=\".\";</SCRIPT>\n";
										$before=$now;
									}
								}
							}
						);
						close DAT;
					}
					###> Percolator part : retrieve the *.pop files from cache folder
					open(TARGET,">$popTargetFile") or $okFile=0;
					if ($okFile) {
						$agent->request(HTTP::Request->new(GET => "$mascotServers{$mascotServer}{url}/cgi/myproms4datFiles.pl?ACT=get&DAT=$filePath&POP=target"),
							sub {
								my ($chunk, $res) = @_;
								print TARGET $chunk;
								#unless (defined $fileLength) {
								#	$fileLength = $res->content_length || 0;
								#}
								#if ($fileLength) {
								#	$bytesReceived += length($chunk);
								#	my $progress=100 * $bytesReceived / $fileLength;
								#	if ($progress-$prevProgress > 5) { # in %
								#		my $progressStrg=sprintf "%d%%",$progress;
								#		print "<SCRIPT type=\"text/javascript\">document.getElementById('$file').innerHTML=\"...$progressStrg\";</SCRIPT>\n";
								#		$prevProgress=$progress;
								#	}
								#}
								#else {
								#	my $now=time;
								#	if ($now-$before > 30) { # in seconds
								#		print "<SCRIPT type=\"text/javascript\">document.getElementById('$file').innerHTML+=\".\";</SCRIPT>\n";
								#		$before=$now;
								#	}
								#}
							}
						);
						close TARGET;
					}
					open(DECOY,">$popDecoyFile") or $okFile=0;
					if ($okFile) {
						$agent->request(HTTP::Request->new(GET => "$mascotServers{$mascotServer}{url}/cgi/myproms4datFiles.pl?ACT=get&DAT=$filePath&POP=decoy"),
							sub {
								my ($chunk, $res) = @_;
								print DECOY $chunk;
								#unless (defined $fileLength) {
								#	$fileLength = $res->content_length || 0;
								#}
								#if ($fileLength) {
								#	$bytesReceived += length($chunk);
								#	my $progress=100 * $bytesReceived / $fileLength;
								#	if ($progress-$prevProgress > 5) { # in %
								#		my $progressStrg=sprintf "%d%%",$progress;
								#		print "<SCRIPT type=\"text/javascript\">document.getElementById('$file').innerHTML=\"...$progressStrg\";</SCRIPT>\n";
								#		$prevProgress=$progress;
								#	}
								#}
								#else {
								#	my $now=time;
								#	if ($now-$before > 30) { # in seconds
								#		print "<SCRIPT type=\"text/javascript\">document.getElementById('$file').innerHTML+=\".\";</SCRIPT>\n";
								#		$before=$now;
								#	}
								#}
							}
						);
						close DECOY;
					}
				}
			}
			#my $okFile=(-e $newDatFile)? `head -1 $newDatFile | grep -c Mascot` : 0;
			my $okFile=(-e $newDatFile)? `tail -1 $newDatFile | grep -c gc0p4Jq0M2Yt08jU534c0p` : 0;
			my $imageFlag='good.gif';
			unless ($okFile) {
				unlink $newDatFile;
				$imageFlag='bad.gif';
			}
			print "<SCRIPT type=\"text/javascript\">document.getElementById('$file').innerHTML=\"&nbsp;<IMG src='$promsPath{images}/$imageFlag'>\";</SCRIPT>\n";
		}
		sleep 3;
		print qq |
<SCRIPT type="text/javascript">
//alert('Done.');
parent.location="$promsPath{cgi}/importBatchAnalyses.cgi?ID=$ID&ACT=UseUserDirectory";
</SCRIPT>
</BODY>
</HTML>
|;
	}
	exit;
}

###################################################
####>Fetching files fromy shared data director<####
###################################################
elsif ($action eq 'UseSharedDirectory') {
	
	####>Starting HTML<####
	print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Fetching Shared Files</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY>
<BR>
<FONT class="title3">Retrieving files from shared data directory...</FONT><BR><BR>
|;
	my $batchFilesDir="$batchDir/$userID";
	foreach my $sharedFile (param('sharedDirFiles')) {
		#my $fileName=(split(/\//,$sharedFile))[-1];
		my (undef,$path,$fileName)=splitpath($sharedFile);
		print "<B>&nbsp;-$fileName...";
		my $newFile="$batchFilesDir/$fileName";
		move("$promsPath{shared}/$sharedFile",$newFile);
		if ($fileName=~/\.zip$/) {
			print " Extracting files...<BR>\n";
			#system "unzip -q -d $batchFilesDir $newFile"; # Inflating archive
			#unlink $newFile; # Deleting archive
			&promsMod::unzipArchive($newFile,$batchFilesDir,{move2Top=>1,mode=>'verbose',txtBefore=>'&nbsp;&nbsp;-',txtAfter=>"<BR>\n"});
			print "&nbsp;Done.</B><BR>\n";
		}
		elsif ($fileName =~ /\.gz$/) {
			print " Extracting files...";
			if ($newFile =~ /\.tar\.gz\Z/) {
				system "tar -zxf $newFile -C $batchFilesDir";
			}
			elsif ($newFile =~ /\.gz\Z/) {
				system "gunzip -c $newFile > $batchFilesDir";
			}
			print " Done.</B><BR>\n";
		}
		else {print " Done.</B><BR>\n";}
	}
	sleep 3;
	print qq |
<SCRIPT type="text/javascript">
window.location="$promsPath{cgi}/importBatchAnalyses.cgi?ID=$ID&ACT=UseUserDirectory";
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}


if ($action ne 'start' && $action ne 'delete') {

	#######################
	####>Starting HTML<####
	#######################
	my $bgColor = $color2;
	print header(-'content-encoding'=>'no',-charset=>'utf-8');warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Search file sources</TITLE>
<LINK rel=stylesheet href="$promsPath{html}/promsStyle.css" type=text/css>
<STYLE type="text/css">
.msf {
	border-style:solid;
	border-color:#999999;
}
.msf_tBorder {border-width: 2px 0px 0px 0px;}
.msf_trBorder {border-width: 2px 2px 0px 0px;}
.msf_rBorder {border-width: 0px 2px 0px 0px;}
.msf_rbBorder {border-width: 0px 2px 2px 0px;}
.msf_bBorder {border-width: 0px 0px 2px 0px;}
.msf_blBorder {border-width: 0px 0px 2px 2px;}
.msf_lBorder {border-width: 0px 0px 0px 2px;}
.msf_ltBorder {border-width: 2px 0px 0px 2px;}
.msf_ltrBorder {border-width: 2px 2px 0px 2px;}
.msf_lrBorder {border-width: 0px 2px 0px 2px;}
.msf_noBorder {}
</STYLE>
<SCRIPT type="text/javascript">
function cancelAction() {
	//top.promsFrame.selectedAction='summary'; // set default action to 'summary'
	top.promsFrame.optionFrame.selectOption();
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<DIV id="waitDiv" style="position:absolute;top:40%;width:100%">
<TABLE align=center>
<TR><TH>Scanning files. Please, wait<SPAN id="waitSpan">...</SPAN></TH></TR>
<TR><TH><IMG src="$promsPath{images}/scrollbarGreen.gif"/></TH></TR>
</TABLE>
</DIV>
|;

	####################################
	####>Fetching others parameters<####
	####################################
	my ($batchFilesDir,$batchDirLabel,$checkDeleteStrg);
	if ($action eq 'UseServerDirectory') {
		$batchFilesDir=param('newDir');
		$batchFilesDir=~s/^\s+//; # remove starting spaces
		$batchFilesDir=~s/\s+\Z//; # remove trailing spaces
		$batchFilesDir="/$batchFilesDir" unless $batchFilesDir=~/^\//; # path must start from root dir '/'
		if (!-e $batchFilesDir) {
			&printError("Directory '$batchFilesDir' not found on server!");
			exit;
		}
		$batchDirLabel="'$batchFilesDir'";
		$checkDeleteStrg='disabled';
	}
	
	else {
		my $selUserDir=($action eq 'UseUserDirectory' && param('userDir'))? param('userDir') : $userID;
		$batchFilesDir="$batchDir/$selUserDir";
		#mkdir $batchFilesDir unless -e $batchFilesDir;
		$batchDirLabel="'$selUserDir' directory";
		$checkDeleteStrg=($selUserDir eq $userID)? 'checked' : 'disabled';
	}

	##########################
	####>Connecting to DB<####
	##########################
	my $dbh=&promsConfig::dbConnect;

	########################################
	####>Recovering information from DB<####
	########################################
	my (%anaParentList,%DBlist,%DBmostUsed);  #%listSample,
	my @itemAnalyses; # decoyData only
	my $projectID=&promsMod::getProjectID($dbh,$itemID,$item);

	if ($action ne 'clean') {	# no need for DB info if selection from files to clean
		if ($decoyFile) { # 100% decoy analyses
			my @sthItem;
			if ($item eq 'sample') {
				$sthItem[0]=$dbh->prepare("SELECT ID_ANALYSIS,'',NAME FROM ANALYSIS WHERE ID_SAMPLE=$itemID AND VALID_STATUS<2 AND DECOY IS NULL ORDER BY DISPLAY_POS ASC");
			}
			elsif ($item eq 'spot') {
				$sthItem[0]=$dbh->prepare("SELECT ID_ANALYSIS,'',ANALYSIS.NAME FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_SPOT=$itemID AND VALID_STATUS<2 AND DECOY IS NULL ORDER BY ANALYSIS.DISPLAY_POS ASC");
			}
			elsif ($item eq 'gel2d') {
				$sthItem[0]=$dbh->prepare("SELECT ID_ANALYSIS,SPOT.NAME,ANALYSIS.NAME FROM ANALYSIS,SAMPLE,SPOT WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT=SPOT.ID_SPOT AND ID_GEL2D=$itemID AND VALID_STATUS<2 AND DECOY IS NULL ORDER BY SPOT.NAME ASC, ANALYSIS.DISPLAY_POS ASC");
			}
			elsif ($item eq 'experiment') {
				$sthItem[0]=$dbh->prepare("SELECT ID_ANALYSIS,CONCAT_WS('>',GEL2D.NAME,SPOT.NAME),ANALYSIS.NAME FROM ANALYSIS,SAMPLE,SPOT,GEL2D WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT=SPOT.ID_SPOT AND SPOT.ID_GEL2D=GEL2D.ID_GEL2D AND GEL2D.ID_EXPERIMENT=$itemID AND VALID_STATUS<2 AND DECOY IS NULL ORDER BY GEL2D.DISPLAY_POS ASC, SPOT.NAME ASC, ANALYSIS.DISPLAY_POS ASC");
				$sthItem[1]=$dbh->prepare("SELECT ID_ANALYSIS,SAMPLE.NAME,ANALYSIS.NAME FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT IS NULL AND ID_EXPERIMENT=$itemID AND VALID_STATUS<2 AND DECOY IS NULL ORDER BY SAMPLE.DISPLAY_POS ASC, ANALYSIS.DISPLAY_POS ASC");
			}
			foreach my $sth (@sthItem) {
				$sth->execute;
				while (my @anaData=$sth->fetchrow_array) {push @itemAnalyses,\@anaData;}
				$sth->finish;
			}
		}
		else { # normal analyses
			if ($item eq 'experiment') {
				##>Free Samples
				my $sthSamp=$dbh->prepare("SELECT ID_SAMPLE,NAME,DISPLAY_POS FROM SAMPLE WHERE ID_SPOT IS NULL AND ID_EXPERIMENT=$itemID");
				$sthSamp->execute;
				while (my ($sampID,$name,$pos)=$sthSamp->fetchrow_array) {
					#@{$listSample{$sampID}}=($name,$pos);
					@{$anaParentList{'SAMPLE'}{$sampID}}=($name,$pos);
				}
				$sthSamp->finish;

				##>Spots
				my $sthSpot=$dbh->prepare("SELECT ID_SPOT,GEL2D.NAME,SPOT.NAME FROM GEL2D,SPOT WHERE GEL2D.ID_GEL2D=SPOT.ID_GEL2D AND ID_EXPERIMENT=$itemID ORDER BY SPOT.NAME ASC");
				$sthSpot->execute;
				my $pos=0;
				while (my ($spotID,$gelName,$spotName)=$sthSpot->fetchrow_array) {
					@{$anaParentList{'SPOT'}{$spotID}}=("$gelName > $spotName",++$pos);
				}
				$sthSpot->finish;
			}
			elsif ($item eq 'gel2d') {
				##>Spots
				my $sthSpot=$dbh->prepare("SELECT ID_SPOT,NAME FROM SPOT WHERE ID_GEL2D=$itemID ORDER BY SPOT.NAME ASC");
				$sthSpot->execute;
				my $pos=0;
				while (my ($spotID,$name)=$sthSpot->fetchrow_array) {
					@{$anaParentList{'SPOT'}{$spotID}}=($name,++$pos);
				}
				$sthSpot->finish;
			}
			elsif ($item eq 'spot') {
				@{$anaParentList{'SPOT'}{$itemID}}=$dbh->selectrow_array("SELECT NAME,1 FROM SPOT WHERE ID_SPOT=$itemID");
			}
			elsif ($item eq 'sample') {
				#@{$listSample{$itemID}}=$dbh->selectrow_array("SELECT NAME,DISPLAY_POS FROM SAMPLE WHERE ID_SAMPLE=$itemID");
				@{$anaParentList{'SAMPLE'}{$itemID}}=$dbh->selectrow_array("SELECT NAME,DISPLAY_POS FROM SAMPLE WHERE ID_SAMPLE=$itemID");
			}

			####>Databanks<###
			my %DBsource;
			my $sthDB=$dbh->prepare("SELECT D.ID_DATABANK,D.NAME,VERSION_NAME,FASTA_FILE,DT.NAME FROM DATABANK D,DATABANK_TYPE DT WHERE D.ID_DBTYPE=DT.ID_DBTYPE AND USE_STATUS='yes'");
			$sthDB->execute;
			while (my ($dbID,$name,$version,$fastaFile,$dbankType)=$sthDB->fetchrow_array) {
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
				$DBlist{$dbSource}{$dbID}=$name;
				$DBlist{$dbSource}{$dbID}.=" ($version)" if $version;
				$DBlist{$dbSource}{$dbID}.=" [$dbankType]";
				$DBsource{$dbID}=$dbSource;
			}
			$sthDB->finish;
			###>Databanks in project<###
			if (scalar keys %DBsource > 5) {
				my $projectID=&promsMod::getProjectID($dbh,$itemID,$item);
				my $sthProjDB=$dbh->prepare("SELECT ID_DATABANK,COUNT(*) AS OCC,MAX(A.START_DATE) AS DT FROM ANALYSIS_DATABANK AD
												INNER JOIN ANALYSIS A ON A.ID_ANALYSIS=AD.ID_ANALYSIS
												INNER JOIN SAMPLE S ON S.ID_SAMPLE=A.ID_SAMPLE
												INNER JOIN EXPERIMENT E ON E.ID_EXPERIMENT=S.ID_EXPERIMENT
												WHERE ID_PROJECT=? GROUP BY ID_DATABANK ORDER BY DT DESC, OCC DESC");
				$sthProjDB->execute($projectID);
				my $rank=0;
				while (my ($dbID,$occ,$lastTimeUsed)=$sthProjDB->fetchrow_array) {
					next unless $DBsource{$dbID};
					@{$DBmostUsed{$dbID}}=($occ,++$rank,$DBsource{$dbID});
				}
				$sthProjDB->finish;
			}
		}
		&updateWaitBox;
	}

	$dbh->disconnect;

	if ($action eq 'UseZipFile') {
		my $uplFile=(split(/\\/,param('uplArch')))[-1]; # Windows
		$uplFile=(split(/\//,$uplFile))[-1]; # Others
		my $newFile="$batchFilesDir/$uplFile";
		my $tmpfile = tmpFileName(upload("uplArch")); # name of temp file being uploaded
		move($tmpfile,$newFile);
		system "unzip -q -d $batchFilesDir $newFile"; # Inflating archive
		unlink $newFile; # Deleting archive
		&updateWaitBox;
	}
	elsif ($action eq 'UseUploadedFiles') {
		foreach my $i (1..$maxUpFiles) {
			next unless param("uploaded_file_$i");
			my $uplFile=(split(/\\/,param("uploaded_file_$i")))[-1]; # Windows
			$uplFile=(split(/\//,$uplFile))[-1]; # Others
			my $newFile="$batchFilesDir/$uplFile";
			my $tmpfile = tmpFileName(upload("uploaded_file_$i")); # name of temp file being uploaded
			move($tmpfile,$newFile);
		}
		&updateWaitBox;
	}


	#######################
	####>Parsing Files<####
	#######################
	my %filesInfo;
	my %taxoLink;
	my %msTypes=&promsConfig::getMsType;

	opendir (DIR, $batchFilesDir) || &printError("Unable to read '$batchFilesDir'!");
	chdir ($batchFilesDir);
	my $currentDatFile;
	my $phenyxTaxo=0;
	while (defined ($currentDatFile = readdir (DIR))) {
		next if -d "$batchFilesDir/$currentDatFile"; # directory
		# next if ($action eq 'UseServerDirectory' && $currentDatFile !~ /\.(dat|xml|msf)$/i);
		#next if ($currentDatFile =~ /.+\.tmp\Z/); # avoid to read in a .tmp file that was previously created in the /tmp directory
		$filesInfo{$currentDatFile}{'OK'}=0; # default
		#my $localError = 0;

		###>MSF from Proteome Discoverer<###
		if ($currentDatFile =~ /.+\.msf\Z/i) {
			my $dbsqlite = DBI->connect( "dbi:SQLite:$currentDatFile", "", "", {PrintError=>0,RaiseError=>0}) || next; # if the file isn't opened, it's not a sqlite file!
			my ($proteomeDiscovererVersion)=$dbsqlite->selectrow_array("SELECT SoftwareVersion FROM SchemaInfo ORDER BY rowid ASC LIMIT 1"); # RaiseError=>1 causes SQLite error in next query before error handling!
			$proteomeDiscovererVersion=~s/^(\d\.\d+).*/$1/; # x.x.x.xx -> x.x (numerical value)
			#my ($myTable) = $dbsqlite->selectrow_array("SELECT sql FROM sqlite_master WHERE type='table' AND tbl_name='ProcessingNodes'") || next;
			my ($myTable) = $dbsqlite->selectrow_array("SELECT sql FROM sqlite_master WHERE type='table' AND tbl_name='ProcessingNodes'");
			next unless ((defined($myTable) && $myTable =~ /FriendlyName/) || $proteomeDiscovererVersion >= 2.0); # 2.x

			$filesInfo{$currentDatFile}{'OK'}=1; # version is OK!
			$filesInfo{$currentDatFile}{'FILE_FORMAT'}='PROTEOME_DISCOVER.MSF';
			$filesInfo{$currentDatFile}{'type'}="";
			%{$filesInfo{$currentDatFile}{'searches'}}=();
			my ($rawfilenumber,$currentSearch,$fileInfoTable);
			my %searches=();

			if ($proteomeDiscovererVersion >= 2.0) { # PD 2.x
				my ($wfXML) = $dbsqlite->selectrow_array("SELECT WorkflowXML FROM Workflows");
				my $xml = new XML::Simple();
				my $xmlData = $xml->XMLin($wfXML);
				my ($parenProcessingNumber,$processingNodeName,$processingNodeNumber);
				# Sort processingNodeParameters so as to get the rawfilenumber that comes first
				foreach my $processingNodeParameters (sort{$a->{ParentProcessingNodeNumber} cmp $b->{ParentProcessingNodeNumber}} @{$xmlData->{WorkflowTree}->{WorkflowNode}}) {
					$rawfilenumber=$processingNodeParameters->{ProcessingNodeNumber} if $processingNodeParameters->{ProcessingNodeName} eq 'SpectrumFilesNode';
					next unless $processingNodeParameters->{FriendlyName} eq 'Mascot' or uc($processingNodeParameters->{FriendlyName}) =~ 'SEQUEST';
					$processingNodeNumber=$processingNodeParameters->{ProcessingNodeNumber};
					$filesInfo{$currentDatFile}{'searches'}{$processingNodeNumber}{'DESCRIPTOR'}=$processingNodeParameters->{FriendlyName};
					$currentSearch=$processingNodeNumber;
					$filesInfo{$currentDatFile}{'type'}.="$currentSearch,$rawfilenumber\t";
					$filesInfo{$currentDatFile}{'SEARCHES'}{$processingNodeNumber}{'MOUSE_STRING'}.=&displaySearchParameters($dbsqlite,$currentSearch,$rawfilenumber,$xmlData); #,$descriptor
					&updateWaitBox;
				}
				$fileInfoTable='WorkflowInputFiles';
			}
			else {
				($rawfilenumber) = $dbsqlite->selectrow_array("SELECT ProcessingNodeNumber FROM ProcessingNodes WHERE FriendlyName='Spectrum Files'");
				###> For PD1.4 -> new name in the database in ProcessingNodes: NodeName='IseNode' & FriendlyName='Sequest HT' for SEQUEST searches
				my $sthGetParams;
				if ($proteomeDiscovererVersion >= 1.4) {
					# Need to work for Mascot, Sequest HT and SEQUEST : Sequest HT is multi-processing algorithm and SEQUEST is the old single processor one
					$sthGetParams = $dbsqlite->prepare("SELECT ProcessingNodeNumber, FriendlyName FROM ProcessingNodes WHERE (UPPER(FriendlyName) LIKE 'SEQUEST%' OR NodeName='Mascot')");
				}
				elsif ($proteomeDiscovererVersion >= 1.2) { # 1.2 ~ 1.3
					$sthGetParams = $dbsqlite->prepare("SELECT ProcessingNodeNumber, FriendlyName FROM ProcessingNodes WHERE (NodeName='SequestNode' OR NodeName='Mascot')");
				}
				$sthGetParams->execute;

				while(my ($processingNodeNumber,$descriptor) = $sthGetParams->fetchrow_array){
					#$descriptor=~s/Perform a //; $descriptor=~s/\.\Z//;
					$filesInfo{$currentDatFile}{'searches'}{$processingNodeNumber}{'DESCRIPTOR'}=$descriptor;
					$currentSearch=$processingNodeNumber;
					$filesInfo{$currentDatFile}{'type'}.="$currentSearch,$rawfilenumber\t";
					$filesInfo{$currentDatFile}{'SEARCHES'}{$processingNodeNumber}{'MOUSE_STRING'}.=&displaySearchParameters($dbsqlite,$currentSearch,$rawfilenumber); #,$descriptor
					&updateWaitBox;
				}
				$fileInfoTable='FileInfos';
				$sthGetParams->finish;
			}
			($filesInfo{$currentDatFile}{'MERGE_FILE_NUMBER'})=$dbsqlite->selectrow_array("SELECT COUNT(DISTINCT FileID) FROM $fileInfoTable");
			#($filesInfo{$currentDatFile}{'MERGE_FILE_IDS'})=$dbsqlite->selectrow_array("SELECT GROUP_CONCAT(FileID,';') FROM $fileInfoTable");
			###> Get Readable Name for MERGE Files 2018/12/20
			my $sthGetFileNames=$dbsqlite->prepare("SELECT FileID,Filename FROM $fileInfoTable");
			$sthGetFileNames->execute;
			my @tableFiles=();
			while(my ($fileID,$rawfileName) = $sthGetFileNames->fetchrow_array){
				my @filePath=split(/\\/,$rawfileName);
				my $searchFile=pop(@filePath);
				$searchFile=~s/.raw//g;
				push @tableFiles,"${fileID}££$searchFile";
			}
			$sthGetFileNames->finish;
			($filesInfo{$currentDatFile}{'MERGE_FILE_IDS'})=join(';',@tableFiles);
			$dbsqlite->disconnect;
			&updateWaitBox;
			next;
		}

		next if -B $currentDatFile; # only plain text files allowed from now on

		###>X! TANDEM<###
		if ($currentDatFile=~/.tandem.pep.xml/) {
			open(TDMFILE,"<",$currentDatFile) or die ("open: $!");
			my $DBpath;
			while (<TDMFILE>) {
                if ($_=~/spectrum, path" value="(.*)"\/>/) {
                    ($filesInfo{$currentDatFile}{'peakFile'}=$1)=~s/.+\///;
                }
                elsif($_=~/list path, sequence source #1" value="(.*)"\/>/){
					($DBpath=$1)=~s/.+\///;
					push @{$filesInfo{$currentDatFile}{'DB'}},$DBpath;
				}
				last if $DBpath && $filesInfo{$currentDatFile}{'peakFile'};
            }
            close TDMFILE;

			#my $xml=XML::Simple->new(KeepRoot=>1);
			#my $xmlData=$xml->XMLin($currentDatFile);
			#(my $peakFile=$xmlData->{'msms_pipeline_analysis'}->{'msms_run_summary'}->{'search_summary'}->{'parameter'}->{'spectrum, path'}->{'value'})=~s/.+\///;
			#$filesInfo{$currentDatFile}{'peakFile'}=$peakFile;
			#(my $DBpath=$xmlData->{'msms_pipeline_analysis'}->{'msms_run_summary'}->{'search_summary'}->{'search_database'}->{'local_path'})=~s/.+\///;
			#push @{$filesInfo{$currentDatFile}{'DB'}},$DBpath;

			$filesInfo{$currentDatFile}{'FILE_FORMAT'}='.TPX';
			$filesInfo{$currentDatFile}{'type'}='X! TANDEM MS/MS';
			$filesInfo{$currentDatFile}{'title'}='';
			$filesInfo{$currentDatFile}{'taxonomy'}='';
			$filesInfo{$currentDatFile}{'OK'}=1;
			next;
		}
		elsif ($currentDatFile=~/.(tandem|xml)/){
			open(TDMFILE,"<",$currentDatFile) or die ("open: $!");
			my $DBpath;
			while (<TDMFILE>) {
                if ($_=~/spectrum, path">(.*)<\/note>/) {
                    ($filesInfo{$currentDatFile}{'peakFile'}=$1)=~s/.+\///;
                }
                elsif($_=~/list path, sequence source #1">(.*)<\/note>/){
					($DBpath=$1)=~s/.+\///;
					push @{$filesInfo{$currentDatFile}{'DB'}},$DBpath;
				}
				last if $DBpath && $filesInfo{$currentDatFile}{'peakFile'};
            }
            close TDMFILE;


#			my $xml=XML::Simple->new(KeepRoot=>1);
#			my $xmlData=$xml->XMLin($currentDatFile);
#			foreach my $group (keys ($xmlData->{'bioml'}->{'group'})){
#				if ($xmlData->{'bioml'}->{'group'}[$group]->{'label'}=~/input parameters/ || $xmlData->{'bioml'}->{'group'}[$group]->{'label'}=~/performance parameters/){
#					foreach my $note (keys ($xmlData->{'bioml'}->{'group'}[$group]->{'note'})){
#						my $DBpath;
#						if ($xmlData->{'bioml'}->{'group'}[$group]->{'note'}[$note]->{'label'}=~/spectrum, path/) {
#							($filesInfo{$currentDatFile}{'peakFile'}=$xmlData->{'bioml'}->{'group'}[$group]->{'note'}[$note]->{'content'})=~s/.+\///;
#                        }
#                        elsif($xmlData->{'bioml'}->{'group'}[$group]->{'note'}[$note]->{'label'}=~/list path, sequence source #1/){
#							($DBpath=$xmlData->{'bioml'}->{'group'}[$group]->{'note'}[$note]->{'content'})=~s/.+\///;
#							push @{$filesInfo{$currentDatFile}{'DB'}},$DBpath;
#						}
#						last if $DBpath && $filesInfo{$currentDatFile}{'peakFile'};
#					}
#				}
#			}
			$filesInfo{$currentDatFile}{'FILE_FORMAT'}='.T';
			$filesInfo{$currentDatFile}{'type'}='X! TANDEM MS/MS';
			$filesInfo{$currentDatFile}{'title'}='';
			$filesInfo{$currentDatFile}{'taxonomy'}='';
			$filesInfo{$currentDatFile}{'OK'}=1;
			next;
		}

		###>MASCOT-LIKE<###
		open (DATFILE, $currentDatFile) || next;
		my $existType=0;
		foreach my $line (<DATFILE>) {
			$line=~s/\s*\Z//;# remove trailing space,\n, etc..
			if ($existType==0) {
				if ($line=~/MIME-Version: .+ Mascot/) {
					$filesInfo{$currentDatFile}{'FILE_FORMAT'}='MASCOT.DAT';
					$existType=1;
				}
				elsif ($line=~/<mascot_search_results/) {
					$filesInfo{$currentDatFile}{'FILE_FORMAT'}='MASCOT.XML';
					$existType=2;
				}
				elsif ($line=~/<IdentificationResult .*version=/) {
					$filesInfo{$currentDatFile}{'FILE_FORMAT'}='PHENYX.XML';
					$filesInfo{$currentDatFile}{'type'}='Phenyx MS/MS';
					$existType=3;
				}
				elsif ($line=~/<PARAGON_VERSION>/) { #<RESULTS
					$filesInfo{$currentDatFile}{'FILE_FORMAT'}='PARAGON.XML';
					$filesInfo{$currentDatFile}{'type'}='Paragon MS/MS';
					$filesInfo{$currentDatFile}{'title'}='';
					$existType=4;
				}
			}
			elsif ($existType==1) { # MASCOT.DAT
				if ($line=~ /^COM=(.*)/) {
					$filesInfo{$currentDatFile}{'title'} = $1;
				}
				elsif ($line=~ /^DB\d*=(.*)/) {
					push @{$filesInfo{$currentDatFile}{'DB'}},$1;
				}
				elsif ($line=~/^FILE=/) {
					if ($line=~/([^ =\/\\]+)\Z/) {
						$filesInfo{$currentDatFile}{'peakFile'} = $1;
					}
					elsif ($line=~/^FILE=(.+)/) {
						$filesInfo{$currentDatFile}{'peakFile'}=$1;
					}
					else {
						$filesInfo{$currentDatFile}{'peakFile'} = 'Unknown';
					}
				}
				#elsif ($line=~/^DECOY=(\S+)/) {
				#	$filesInfo{$currentDatFile}{'decoy'} = $1;
				#}
				elsif ($line=~/^TAXONOMY=[\. ]*(.+)/) {
					my $taxonomy=$1;
					$taxonomy=~s/ \(.*//;# remove eg. (human)
					$filesInfo{$currentDatFile}{'taxonomy'} = $taxonomy;
				}
				elsif ($line=~/^SEARCH=(\S+)/) {
					$filesInfo{$currentDatFile}{'type'} = $msTypes{$1};
				}
				elsif ($line=~/^MODS=(.+)/) {
					$filesInfo{$currentDatFile}{'otherInfo'}{'Fixed modifications'} = $msTypes{$1};
				}
				elsif ($line=~ /^IT_MODS=(.+)/) {
					$filesInfo{$currentDatFile}{'otherInfo'}{'Variable modifications'} = $1 ;
				}
				elsif ($line=~ /^CLE=(.+)/) {
					$filesInfo{$currentDatFile}{'otherInfo'}{'Enzyme'} = $1 ;
				}
				elsif ($line=~ /^PFA=(\d+)/) {
					$filesInfo{$currentDatFile}{'otherInfo'}{'Enzyme misscleavage'} = $1 ;
				}
				elsif ($line=~ /^INSTRUMENT=(.+)/ ){
					$filesInfo{$currentDatFile}{'otherInfo'}{'Instrument'} = $1 ;
				}
				elsif ($line=~ /^TOL=(.+)/ ){
					$filesInfo{$currentDatFile}{'otherInfo'}{'Peptide tolerance'} = $1 ;
				}
				elsif ($line=~ /^TOLU=(.+)/) {
					$filesInfo{$currentDatFile}{'otherInfo'}{'Peptide tolerance'} .= $1 ;
				}
				elsif ($line=~ /^ITOL=(.+)/) {
					$filesInfo{$currentDatFile}{'otherInfo'}{'Fragment tolerance'} = $1 ;
				}
				elsif ($line=~ /^ITOLU=(.+)/ ){
					$filesInfo{$currentDatFile}{'otherInfo'}{'Fragment tolerance'} .= $1 ;
				}
				elsif ($line=~ /^release=(.+)/) {
					$filesInfo{$currentDatFile}{'otherInfo'}{'Database'} = $1 ;
				}
				elsif ($line=~ /^queries=(\d+)/) {
					$filesInfo{$currentDatFile}{'otherInfo'}{'Queries'} = $1 ;
				}
				last if $line=~ /name="masses"/;
			}
			elsif ($existType==2) { # MASCOT.XML
				if ($line=~ /<COM>(.*)<\/COM>/) {
					$filesInfo{$currentDatFile}{'title'} = $1;
				}
				elsif ($line=~ /<DB>(.*)<\/DB>/) {
					push @{$filesInfo{$currentDatFile}{'DB'}},$1;
				}
				elsif ($line=~/<FILENAME>(.*)<\/FILENAME>/) {
					($filesInfo{$currentDatFile}{'peakFile'}) = ($1=~/.*[\/|\\](\w+\.\w+)/);
				}
				elsif ($line=~/<TAXONOMY>\W*(.*)<\/TAXONOMY>/) {
					my $taxonomy=$1;
					$taxonomy=~s/ \(.*//;# remove eg. (human)
					$taxonomy=~s/\s*\Z//;# remove trailing space,\n, etc..
					$filesInfo{$currentDatFile}{'taxonomy'} = $taxonomy;
					}
				elsif ($line=~/<SEARCH>(.*)<\/SEARCH>/) {
					$filesInfo{$currentDatFile}{'type'} = $msTypes{$1};
				}
				last if $line=~ /<\/search_parameters>/;
			}
			elsif ($existType==3) { # PHENYX.XML
				if ($line=~ /\Wtitle>.*CDATA\[(.*)\]\]/) {
					$filesInfo{$currentDatFile}{'title'} = $1;
				}
				elsif ($line=~ /<database db="([^"]+)"/) {
					push @{$filesInfo{$currentDatFile}{'DB'}},$1;
				}
				elsif ($line=~/<clientPath>.*CDATA\[(\w+\.\w+)\]/) {
					$filesInfo{$currentDatFile}{'peakFile'} = $1;
				}
				elsif ($line=~/<taxoCriterion>(\d+)/) {
					$filesInfo{$currentDatFile}{'taxonomy'} = $1;
					$taxoLink{$1}=1;
					$phenyxTaxo=1;
				}
				last if $line=~ /<\/[idl:]*SubmissionParam>/;
			}
			elsif ($existType==4) { # PARAGON.XML
				if ($line =~ /<PROGROUP .* fasta=\"(.*?)\"/) {
					my @filePath=split(/\\/,$1); # ex: C:\AB SCIEX\ProteinPilot Data\SearchDatabases\SwissProt_20110610.fasta
					my $toAdd=1;
					foreach my $dbName (@{$filesInfo{$currentDatFile}{'DB'}}) {
						$toAdd = 0 if $dbName eq $filePath[-1];
					}
					push @{$filesInfo{$currentDatFile}{'DB'}},$filePath[-1] if $toAdd; # Avoid to add several times the same DB because XML files are redundant
				}
				elsif ($line =~ /<PEAKLIST .* originalfilename=\"(.*?)\"/){
					my @filePath=split(/\\/,$1);# ex: C:\Documents and Settings\Administrator\Desktop\QStar\E5682BL.wiff
					$filesInfo{$currentDatFile}{'peakFile'}=$filePath[-1];
				}
				elsif($line =~ /<UI_SPECIES>(.*)<\/UI_SPECIES>/){
					$filesInfo{$currentDatFile}{'taxonomy'}=$1;
					last; # no need to read rest of file
				}
				elsif($line =~ /<UI_SPECIES\/>/){
					$filesInfo{$currentDatFile}{'taxonomy'}='All entries';
					last; # no need to read rest of file
				}
			}
		}
		close DATFILE;
		$filesInfo{$currentDatFile}{'OK'}=$existType; # version is OK!
		&updateWaitBox;
	}
	closedir DIR;

	if ($phenyxTaxo==1) {
		my $TaxoFile ="$promsPath{banks}/taxonomy_codes.txt";
		open (TAXOFILE, $TaxoFile) || &printError("Unable to read $TaxoFile!");
		while (<TAXOFILE>) {
			my $line = $_ ;
			if (($line =~/^(\d+)\t([^\t]*)/) && (exists ($taxoLink{$1}))) {
					$taxoLink{$1}=$2;
				}
		}
		close TAXOFILE ;
		foreach my $datFile (keys %filesInfo) {
			if ($filesInfo{$datFile}{'FILE_FORMAT'} eq 'PHENYX.XML') {
				$filesInfo{$datFile}{'taxonomy'}=$taxoLink{$filesInfo{$datFile}{'taxonomy'}} if $taxoLink{$filesInfo{$datFile}{'taxonomy'}} ne 'o' ;
			}
		}
		&updateWaitBox;
	}

	##############
	####>HTML<####
	print qq |
<SCRIPT type="text/javascript">
var checkStatus=true;
function cancelAction() {
	//top.promsFrame.optionFrame.selectOption();
	window.location="./importBatchAnalyses.cgi?ID=$ID&ACT=start";
}
function checkUncheckAll(checkStatus) {
	var anaBox=document.fileListForm.file;
	if (anaBox.length) { // more than 1 checkboxes
		for (i=0; i < anaBox.length; i++){anaBox[i].checked=checkStatus;}
	}
	else {anaBox.checked=checkStatus;} // Only 1 checkbox
}
function testCheckbox(myForm) {
	var selected=0;
	var anaBox=myForm.file;
	if (anaBox.length) { // more than 1 checkboxes
		for (i=0; i < anaBox.length; i++){
			if (anaBox[i].checked==true) {selected=1; break;}
		}
	}
	else if (anaBox.checked==true){selected=1;}
	return selected;
}
function showHideSearchParam(searchRequest,action) {
	if (action=='show') {
		document.getElementById(searchRequest+':show').style.display='block';
		document.getElementById(searchRequest+':hide').style.display='none';
	}
	else {
		document.getElementById(searchRequest+':show').style.display='none';
		document.getElementById(searchRequest+':hide').style.display='block';
	}
}
|;

	###################################
	####>File selection for import<####
	###################################
	if ($action ne 'clean') {
		print "var fileNames=[];\n";
		my $count=0;
		my $disabMZXML=1;
		foreach my $dataFile (sort {&promsMod::sortSmart(lc($a),lc($b))} keys %filesInfo) {

			if ($filesInfo{$dataFile}{'OK'}) {
				(my $searchName=$dataFile)=~s/\.[^\.]+//;
				($searchName=$dataFile)=~s/\.\D+// if $dataFile=~/tandem/;

				my $msName='';
				if ($filesInfo{$dataFile}{'peakFile'}) {
					($msName=$filesInfo{$dataFile}{'peakFile'})=~s/\.[^\.]+//;
					$msName=~ s/\s+\Z//;
				}
				my $numSearches=($filesInfo{$dataFile}{'searches'})? scalar keys %{$filesInfo{$dataFile}{'searches'}} : 1;
				print "fileNames[$count]=['$searchName','$msName',$numSearches];\n";
				$count++;

				if (-e "$promsPath{data}/tmp/upload/project_$projectID/$searchName.mzXML"){
                    $disabMZXML=1;
                }
                else{
					$disabMZXML=($disabMZXML==0) ? 0 : ($filesInfo{$dataFile}{'FILE_FORMAT'} eq '.T' || $filesInfo{$dataFile}{'FILE_FORMAT'} eq '.TPX' || $filesInfo{$dataFile}{'FILE_FORMAT'} eq 'MASCOT.DAT' ) ? 0 : 1 ; ##only for Xtandem and .dat files
				}
			}
			#else {print "fileNames[$count]=['','',0];\n";}
		}

		if ($decoyFile) {
			print qq
|function checkTargAnalysis(mySelTargAna) {
	if (!mySelTargAna.value) return; // no ana selected
	var allSelTargAna=document.fileListForm.selTargAnalysis;
	if (!allSelTargAna.length) return; // only 1 file
	for (var i=0;i<allSelTargAna.length;i++) {
		if (allSelTargAna[i]==mySelTargAna) continue;
		if (allSelTargAna[i].value==mySelTargAna.value) {
			alert('This analysis has already been selected for file: '+document.fileListForm.file[i].value);
			mySelTargAna.selectedIndex=0;
			break;
		}
	}
}
function checkFileForm(myForm) {
	if (testCheckbox(myForm)==0) {
		alert('ERROR: No files selected !');
		return false;
	}
	if (myForm.file.length) { // multiple files to import
		for (var i=0;i<myForm.file.length;i++) {
			if (myForm.file[i].checked && !myForm.selTargAnalysis[i].value) {
				alert('No target Analysis selected for file: '+myForm.file[i].value);
				return false;
			}
		}
	}
	else if (!myForm.selTargAnalysis.value) {
		alert('No target Analysis selected');
		return false;
	}
	return true;
}
|;
		}
		else { # not decoy files
			print qq |
function checkFileForm(myForm) {
	if (testCheckbox(myForm)==0) {
		alert('ERROR: No files selected !');
		return false;
	}
	var missAnaPos=testItemName(myForm.anaName);
	if (missAnaPos != 0) {
		alert('ERROR: Analysis name for file #'+missAnaPos+' is missing !');
		return false;
	}
	var missSampPos=testItemName(myForm.selAnaParent);
	if (missSampPos != 0) {
		alert('ERROR: Sample name for file #'+missSampPos+' is missing !');
		return false;
	}
	if (!myForm.databankID1.value) {
		alert('ERROR: Select at least 1 sequence databank to be used with files !');
		return false;
	}
	for (var i=1; i<=2; i++) {
		var prevDbValue=myForm['databankID'+i].value;
		if (myForm['databankID'+i].disabled==true \|\| !prevDbValue) {break;}
		for (var j=(i*1)+1; j<=3; j++) {
			var nextDbValue=myForm['databankID'+j].value;
			if (myForm['databankID'+j].disabled==true \|\| !nextDbValue) {break;}
			if (nextDbValue==prevDbValue) {
				alert('ERROR: Same databank selected at positions #'+i+' and #'+j);
				return false;
			}
		}
	}
	if (myForm.des.value && myForm.des.value.length>90) {
		alert('Description is too long.');
		return false;
	}
	if (myForm.comments.value && myForm.comments.value.length>240) {
		alert('Comment is too long.');
		return false;
	}
	return true;
}
function testItemName(itemNameArea) {
	var missingNamePos=0;
	var anaBox=document.fileListForm.file;
	if (anaBox.length) { // more than 1 checkboxes
		for (i=0; i < anaBox.length; i++){
			if (anaBox[i].disabled==false && anaBox[i].checked==true && !itemNameArea[i].value) {
				missingNamePos=i+1;
				break;
			}
		}
	}
	else if (anaBox.checked==true){
		if (!itemNameArea.value) {missingNamePos=1;}
	}
	return missingNamePos;
}
function selectAnaName(nameType) {
	var anaField=document.fileListForm.anaName;
	var rowIndex=-1;
	if (fileNames.length > 1) {
		for (var i=0; i<fileNames.length; i++) {
			if (fileNames[i][2]==0) {
				rowIndex++;
				continue;
			}
			for (var s=1; s<=fileNames[i][2]; s++) { // searches (Proteome Discover)
				rowIndex++;
				anaField[rowIndex].value=(nameType=='search')? fileNames[i][0] : (nameType=='ms')? fileNames[i][1] : '';
				if (fileNames[i][2]>1 && anaField[rowIndex].value) anaField[rowIndex].value+=' #'+s;
			}
		}
	}
	else if (fileNames[0][2]>0) { // only 1 file
		if (fileNames[0][2]>1) {
			for (var s=1; s<=fileNames[0][2]; s++) { // searches (Sequest)
				rowIndex++;
				anaField[rowIndex].value=(nameType=='search')? fileNames[0][0] : (nameType=='ms')? fileNames[0][1] : '';
				if (fileNames[0][2]>1 && anaField[rowIndex].value) anaField[rowIndex].value+=' #'+s;
			}
		}
		else { // only 1 search in data file
			anaField.value=(nameType=='search')? fileNames[0][0] : (nameType=='ms')? fileNames[0][1] : '';
		}
	}
}
function changeSampleName(sampleField) {
	if (sampleField.value != 'new') {return;}
	var sampleName='';
	while (sampleName=='') {
		sampleName=prompt('Enter name for new Sample :');
		if (sampleName=='') {alert('ERROR: Sample name is not valid !')};
		if (sampleName==null) {
			sampleField.selectedIndex=0;
			return;
		}
	}
	//Checking if name already exists
	for (var i=1;i<sampleField.length;i++) { // [0]=select
		if (sampleField.options[i].value.match('SPOT')) {continue;}
		if (sampleField.options[i].text==sampleName \|\| sampleField.options[i].text==sampleName+' (new Sample)') {
			if (confirm("WARNING: a Sample named '"+sampleName+"' exists already. Create anyway ?")) {
				break;
			}
			else {
				sampleField.selectedIndex=0;
				return;
			}
		}
	}
	//Adding name as new OPTION to all SELECT
	var sampSelect=document.fileListForm.selAnaParent;
	if (fileNames.length > 1) {
		for (var i=0;i<sampSelect.length;i++) {
			sampSelect[i][sampSelect[i].length]=new Option(sampleName+' (new Sample)','SAMPLE:new='+sampleName);
			if (sampleField==sampSelect[i]) {
				sampleField.selectedIndex=sampleField.length-1;
			}
		}
	}
	else { // only 1 file
		sampleField[sampleField.length]=new Option(sampleName,'SAMPLE:new='+sampleName);
		sampleField.selectedIndex=sampleField.length-1;
	}
}
function updateDatabankSelection(dbNum,dbValue) {
	var targetDisabled=(dbValue)? false : true;
	document.getElementById('db_'+(dbNum+1)).disabled=targetDisabled;
	if (dbNum==1) {
		if (!dbValue) {document.getElementById('db_3').disabled=true;}
		else if (document.getElementById('db_2').value) {document.getElementById('db_3').disabled=false;}
	}
}
|;
		}
		print qq
|document.getElementById('waitDiv').style.display='none';
</SCRIPT>
<CENTER>
<FONT class="title">Select Files to be Imported</FONT>
<BR><FONT class="title2">(From $batchDirLabel)</FONT>
<BR><BR>
|;
		unless (scalar keys %filesInfo) {
			&printError("No search files found!");
			exit;
		}
		print qq
|<FORM name="fileListForm" action="./storeAnalyses.cgi" method="post"  enctype="multipart/form-data" onsubmit="return checkFileForm(this);">
<INPUT type="hidden" name="batchStart" value="$item:$itemID">
<INPUT type="hidden" name="ACT" value="anaBatch">
<INPUT type="hidden" name="decoy" value="$decoyFile">
<INPUT type="hidden" name="fileAct" value="$action">
<INPUT type="hidden" name="path" value="$batchFilesDir">
<INPUT type="hidden" name="PROJECT_ID" value="$projectID">
<TABLE border=0 cellspacing=1 cellpadding=3 width=100%> <TR><TD colspan=3>
	<TABLE border=0 cellspacing=0 cellpadding=0 width=100%>
|;
		if ($action ne 'UseServerDirectory') {
			my $colspan=($decoyFile)? 6 : 7;
			print qq
|	<TR><TD colspan=$colspan valign=middle>
	<INPUT type="checkbox" name="delUpFiles" value=1 $checkDeleteStrg><FONT class="title3" style="color:#DD0000">Delete imported files afterwards.</FONT>
	</TD></TR>
	<TR><TD colspan=$colspan valign=middle>
	<INPUT type="checkbox" name="splitMode" value=1 unchecked><FONT class="title3" style="color:#DD0000">Split mode (needed to separate merge searches in msf files).</FONT>
	</TD></TR>
|;
		}
		print qq
|<TR class="header" bgcolor="$bgColor">
	<TH width=150 align=left><INPUT type="checkbox" onclick="checkUncheckAll(this.checked)">Search file</TH>
|;
		if ($decoyFile) {
			print "\t\t\t<TH nowrap width=200 nowrap>Target Analysis</TH>\n";
		}
		else {
			print qq
|	<TH nowrap width=100 nowrap>Analysis name<BR>
		<SELECT id="anaNameSel" onchange="selectAnaName(this.value)" style="font-size:9px">
		<OPTION value="manual">Manual</OPTION>
		<OPTION value="ms">MS file</OPTION>
		<OPTION value="search">Search file</OPTION>
		</SELECT>
	</TH>
	<TH nowrap width=100>Parents</TH>
|;
		}
		print qq
|	<TH nowrap width=100>MS file</TH>
	<TH nowrap width=100>Search type</TH>
	<TH nowrap width=140>Databank(s)<BR>Taxonomy</TH>

	|;
	print "<TH>mzXML file</TH>" unless $disabMZXML;
	print qq
|	<TH>Title</TH>
</TR>
|;
		my $selectParentString;
		if ($decoyFile) {
			$selectParentString="<SELECT name=\"selTargAnalysis\" style=\"width:190px##DISPLAY##\" onchange=\"checkTargAnalysis(this)\">\n<OPTION value=\"\">-= Select =-</OPTION>\n";
			my %prevItemName;
			my $prevParentStrg='';
			foreach my $refAnaData (@itemAnalyses) {
				my ($anaID,$parentStrg,$anaName)=@{$refAnaData};
				if ($parentStrg ne $prevParentStrg) {
					$selectParentString.="</OPTGROUP>\n" if $prevParentStrg;
					$selectParentString.="<OPTGROUP label=\"$parentStrg:\">\n";
				}
				$selectParentString.="<OPTION value=\"$anaID\">$anaName</OPTION>\n";
				$prevParentStrg=$parentStrg;
			}
			$selectParentString.="</OPTGROUP>\n" if $prevParentStrg;
			$selectParentString.="</SELECT>\n";
		}
		else {
			if ($item eq 'experiment' || $item eq 'gel2d') {
				$selectParentString="<SELECT name=\"selAnaParent\" style=\"width:100px##DISPLAY##\" onchange=\"changeSampleName(this)\">\n<OPTION value=\"\">-= Select =-</OPTION>\n";
				if ($anaParentList{'SPOT'}) {
					$selectParentString.="<OPTGROUP label=\"Spots:\">\n";
					foreach my $spotID (sort{$anaParentList{'SPOT'}{$a}[1]<=>$anaParentList{'SPOT'}{$b}[1]} keys %{$anaParentList{'SPOT'}}) {
						$selectParentString.="<OPTION value=\"SPOT:$spotID\">$anaParentList{SPOT}{$spotID}[0]</OPTION>\n";
					}
					$selectParentString.="</OPTGROUP>\n";
				}
				$selectParentString.="<OPTGROUP label=\"Samples:\">\n<OPTION value=\"new\">-= New =-</OPTION>\n";
				foreach my $sampID (sort{$anaParentList{'SAMPLE'}{$a}[1]<=>$anaParentList{'SAMPLE'}{$b}[1]} keys %{$anaParentList{'SAMPLE'}}) {
					$selectParentString.="<OPTION value=\"SAMPLE:$sampID\">$anaParentList{SAMPLE}{$sampID}[0]</OPTION>\n";
				}
				$selectParentString.="</OPTGROUP>\n";
				#$selectParentString="<OPTION value=\"\">-= Select =-</OPTION>\n<OPTION value=\"new\">-= New =-</OPTION>\n";
				#foreach my $sampleID (sort {$listSample{$a}[1]<=>$listSample{$b}[1]} keys %listSample) { # sort by display pos
				#	$selectParentString .= "<OPTION value=\"$sampleID\">$listSample{$sampleID}[0]</OPTION>\n";
				#}
				$selectParentString.="</SELECT>\n";
			}
			elsif ($item eq 'spot') {
				$selectParentString="<B>$anaParentList{SPOT}{$itemID}[0]</B><INPUT type=\"hidden\" name=\"selAnaParent\" value=\"SPOT:$itemID\">";
			}
			else { # sample
				$selectParentString="<B>$anaParentList{SAMPLE}{$itemID}[0]</B><INPUT type=\"hidden\" name=\"selAnaParent\" value=\"SAMPLE:$itemID\">";
			}
		}
		my $i=0;
		foreach my $dataFile (sort {&promsMod::sortSmart(lc($a),lc($b))} keys %filesInfo) {
			$bgColor=($bgColor eq $color1)? $color2 : $color1;
			#print "<TR class=\"list\" bgcolor=$bgColor valign=middle>\n";
			my ($disabString,$visString)=($filesInfo{$dataFile}{'OK'})? ('','') : ('disabled',';display:none');
			(my $localParentString=$selectParentString)=~s/##DISPLAY##/$visString/;

			#my $sampleSelString;
			#if ($item eq 'experiment' || $item eq 'gel2d') {
			#	$sampleSelString="<SELECT name=\"selSample\" onchange=\"changeSampleName(this)\" $disabString>$sampleOptString</SELECT>\n";
			#}
			#elsif ($item eq 'spot') {
			#	my $spotID=(keys %{$anaParentList{'SPOT'}})[0];
			#	$sampleSelString="<B>$anaParentList{'SPOT'}{$spotID}[0]</B><INPUT type=\"hidden\" name=\"selSample\" value=\"SPOT:$spotID\">";
			#}
			#else { # sample
			#	my $sampID=(keys %{$anaParentList{'SAMPLE'}})[0];
			#	$sampleSelString="<B>$anaParentList{'SAMPLE'}[0]</B><INPUT type=\"hidden\" name=\"selSample\" value=\"SAMPLE:$sampleID\">";
			#}
			###--- TAG INSERTION (start)---###
			if (!$filesInfo{$dataFile}{'OK'}) {
				#my $colspan=($decoyFile)? 5 : 6;
				my $colspan=($decoyFile && $disabMZXML) ? 5 : ($disabMZXML) ? 6 : 7;
				print qq |<TR class="list" bgcolor=$bgColor valign=middle>
	<TH nowrap align=left><INPUT type="checkbox" name="badFile" value="$dataFile" disabled>$dataFile</TH>
	<TH colspan=$colspan align=left><FONT style=\"color:#DD0000\">&nbsp;Unsupported file format or version.</FONT></TH>
</TR>
|;
			}
			elsif ($filesInfo{$dataFile}{'FILE_FORMAT'} eq 'PROTEOME_DISCOVER.MSF') { # SEQUEST
				my $colspan=($decoyFile)? 7 : 8;
				my $mergeFileString=($filesInfo{$dataFile}{'MERGE_FILE_NUMBER'}>1)? " (merge of $filesInfo{$dataFile}{'MERGE_FILE_NUMBER'} raw-files)" : "";
#-----------------------------------				
#				print qq
#|	<TR bgcolor=$bgColor valign=middle>
#		<TD nowrap colspan=$colspan><FIELDSET><LEGEND><B><U>$dataFile</U>$mergeFileString</B></LEGEND>
#		<TABLE>
#|;	#	<INPUT type="hidden" name="fileName" value="$dataFile">
#				#print "<TR class=\"list\" bgcolor=$bgColor valign=middle>\n";
#				foreach my $processingNodeNumber (sort{$a<=>$b} keys %{$filesInfo{$dataFile}{'SEARCHES'}}) {
#					print qq
#|		<TR bgcolor=$bgColor valign=top>
#		<TH nowrap align=left width=135><INPUT type="hidden" name="fileName" value="$dataFile:$processingNodeNumber:$filesInfo{$dataFile}{'MERGE_FILE_IDS'}"><INPUT type="checkbox" name="file" value="$dataFile:$processingNodeNumber:$filesInfo{$dataFile}{'MERGE_FILE_IDS'}">Request $processingNodeNumber</TH>
#|;
#					print "\t\t<TD align=center width=100><INPUT type='text' name='anaName' value='' size=14></TD>\n" unless $decoyFile;
#					print qq
#|		<TD align=center width=100>$localParentString</TD>
#		<TD>&nbsp;<B>$filesInfo{$dataFile}{searches}{$processingNodeNumber}{DESCRIPTOR}</B></TD>
#		<TD><DIV id="$dataFile:$processingNodeNumber:show" style="display:none"><FIELDSET><LEGEND><INPUT type="button" value="Hide Search Parameters" onclick="showHideSearchParam('$dataFile:$processingNodeNumber','hide')" style="font-size:11px"/></LEGEND>$filesInfo{$dataFile}{SEARCHES}{$processingNodeNumber}{MOUSE_STRING}</FIELDSET></DIV>
#			<DIV id="$dataFile:$processingNodeNumber:hide"><INPUT type="button" value="Show Search Parameters" onclick="showHideSearchParam('$dataFile:$processingNodeNumber','show')" style="font-size:11px"/></DIV>
#		</TD>
#		</TR>
#|;
#				}
#				print "\t\t</TABLE></FIELDSET></TD>\n\t</TR>\n";
#-----------------------------------
				
				print qq
|	<TR bgcolor=$bgColor valign=middle><TD class="msf msf_ltrBorder" nowrap colspan=$colspan>&nbsp;<B><U>$dataFile</U>$mergeFileString</B></TD></TR>
|;
				my $colspan2=($decoyFile)? 5 : 6;
				my $numSearches=scalar keys %{$filesInfo{$dataFile}{'SEARCHES'}};
				my $curSearch=0;
				foreach my $processingNodeNumber (sort{$a<=>$b} keys %{$filesInfo{$dataFile}{'SEARCHES'}}) {
					$curSearch++;
					my ($classLeft,$classMid,$classRight)=($curSearch==$numSearches)? ('msf_blBorder','msf_bBorder','msf_rbBorder') : ('msf_lBorder','msf_noBorder','msf_rBorder');
					print qq
|	<TR bgcolor=$bgColor valign=top>
		<TH class="msf $classLeft" nowrap align=left><INPUT type="hidden" name="fileName" value="$dataFile:$processingNodeNumber:$filesInfo{$dataFile}{'MERGE_FILE_IDS'}"><INPUT type="checkbox" name="file" value="$dataFile:$processingNodeNumber:$filesInfo{$dataFile}{'MERGE_FILE_IDS'}">Search #$processingNodeNumber ($filesInfo{$dataFile}{searches}{$processingNodeNumber}{DESCRIPTOR})</TH>
|;
					print "\t\t<TD class='msf $classMid' align=center><INPUT type='text' name='anaName' value='' size=14></TD>\n" unless $decoyFile;
					print qq
|		<TD class="msf $classMid" align=center>$localParentString</TD>
		<TD class="msf $classRight" colspan=$colspan2><DIV id="$dataFile:$processingNodeNumber:show" style="display:none"><FIELDSET><LEGEND><INPUT type="button" value="Hide Search Parameters" onclick="showHideSearchParam('$dataFile:$processingNodeNumber','hide')" style="font-size:11px"/></LEGEND>$filesInfo{$dataFile}{SEARCHES}{$processingNodeNumber}{MOUSE_STRING}</FIELDSET></DIV>
							  <DIV id="$dataFile:$processingNodeNumber:hide"><INPUT type="button" value="Show Search Parameters" onclick="showHideSearchParam('$dataFile:$processingNodeNumber','show')" style="font-size:11px"/></DIV>
		</TD>
	</TR>
|;
					print "\t<TR bgcolor=$bgColor><TH class='msf msf_lrBorder' colspan=$colspan><HR width=80%></TH></TR>\n" if $curSearch < $numSearches;
				}
#-----------------------------------

			###--- TAG INSERTION (end)---###
			}
			else { #(scalar keys %{$filesInfo{$dataFile}})
				print qq
|<TR class="list" bgcolor=$bgColor valign=middle>
	<TH nowrap align=left><INPUT type="hidden" name="fileName" value="$dataFile"><INPUT type="checkbox" name="file" value="$dataFile">$dataFile</TH>
|;
				if ($decoyFile) {print "\t\t<TD align=center>&nbsp;$localParentString&nbsp;</TD>\n";}
				else {
					print qq
|	<TD align=center>&nbsp;<INPUT type="text" name="anaName" value="" size="14">&nbsp;</TD>
	<TD align=center>&nbsp;$localParentString&nbsp;</TD>
|;
				}
				my $dbStrg=join('&bull;',@{$filesInfo{$dataFile}{'DB'}});
				print qq
|	<TD align=center>&nbsp;$filesInfo{$dataFile}{peakFile}&nbsp;</TD>
	<TD align=center>&nbsp;$filesInfo{$dataFile}{type}&nbsp;</TD>
	<TD align=center>&nbsp$dbStrg<BR>$filesInfo{$dataFile}{taxonomy}&nbsp</TD>

	|;

				if ($filesInfo{$dataFile}{'FILE_FORMAT'} eq '.T' || $filesInfo{$dataFile}{'FILE_FORMAT'} eq '.TPX' || $filesInfo{$dataFile}{'FILE_FORMAT'} eq 'MASCOT.DAT'){
					(my $mzxmlFile=$dataFile)=~s/\.\D+//;
					if (-e "$promsPath{data}/tmp/upload/project_$projectID/$mzxmlFile.mzXML"){
						print "<TD align=center><IMG src='$promsPath{images}/good.gif'></TD>"  unless $disabMZXML;
					}else{
						print "<TD align=center>&nbsp;<INPUT type=\"file\" name=\"mzxml$i\" id=\"mzxml$i\" accept=\".mzXML\">&nbsp;</TD>" ;
					}
				}
				else{print "<TD align=center>-</TD>" unless $disabMZXML;}
				$i++;
				print qq |
	<TD>&nbsp;$filesInfo{$dataFile}{title}</TD>
</TR>
|;
			}
			#else { # error on reading file
			#	print "<TH colspan=4 align=left>&nbsp;<FONT style=\"color:#DD0000\">&nbsp;Unknown file format</FONT></TH>\n</TR>\n";
			#}
		}
		print "</TABLE></TD></TR>\n";

		unless ($decoyFile) { # normal import
			my $databaseString="<OPTION selected value=\"\">-= Choose from list =-</OPTION>\n";
			if (scalar keys %DBmostUsed) {
				$databaseString.="<OPTGROUP label=\"Recently used:\">\n";
				my $count=0;
				foreach my $dbID (sort{$DBmostUsed{$b}[1]<=>$DBmostUsed{$a}[1]} keys %DBmostUsed) {
					$count++;
					my $dbSource=$DBmostUsed{$dbID}[2];
					$databaseString.="<OPTION value=\"$dbID\">$DBlist{$dbSource}{$dbID} ($dbSource)</OPTION>\n";
					last if $count==5;
				}
				$databaseString.="</OPTGROUP>\n";
				if (scalar keys %DBmostUsed > 1) {
					$databaseString.="<OPTGROUP label=\"Most used:\">\n";
					$count=0;
					foreach my $dbID (sort{$DBmostUsed{$b}[0]<=>$DBmostUsed{$a}[0]} keys %DBmostUsed) {
						$count++;
						my $dbSource=$DBmostUsed{$dbID}[2];
						$databaseString.="<OPTION value=\"$dbID\">$DBlist{$dbSource}{$dbID} ($dbSource)</OPTION>\n";
						last if $count==5;
					}
					$databaseString.="</OPTGROUP>\n";
				}
			}
			foreach my $dbSource (sort{lc($a) cmp lc($b)} keys %DBlist) {
				$databaseString.="<OPTGROUP label=\"All ".(($dbSource eq 'Local')? 'local' : $dbSource).":\">\n";
				foreach my $dbID (sort{lc($DBlist{$dbSource}{$a}) cmp lc($DBlist{$dbSource}{$b})} keys %{$DBlist{$dbSource}}) {
					$databaseString.="<OPTION value=\"$dbID\">$DBlist{$dbSource}{$dbID}</OPTION>\n";
				}
				$databaseString.="</OPTGROUP>\n";
			}

			my $maxRankOptionString ="";
			foreach my $rank (1..10) {
				#my $selStrg=($rank==$maxRank)? 'selected' : '';
				$maxRankOptionString .="<OPTION value='$rank'>$rank</OPTION>"; #$selStrg
			}
			my $defaultStrg="&nbsp;&nbsp;(<B>default</B> is set to: <B>".&promsConfig::getMinScore('MASCOT')."</B> for Mascot, <B>".&promsConfig::getMinScore('PARAGON')."</B> for Paragon, <B>".&promsConfig::getMinScore('PHENYX')."</B> for Phenyx, <B>".&promsConfig::getMinScore('TANDEM')."</B> for X! Tandem, <B>".&promsConfig::getMinScore('SEQUEST')."</B> for Sequest).";
			my $disabQval=(defined $promsPath{'qvality'})? '' : 'disabled';
			print qq
|</TABLE>
<TABLE width=100%>
	<TR><TD bgcolor=$color2>
	<TABLE width=100%>
	<TR>
		<TH align=right width=100>Databank(s) :</TH>
		<TD bgcolor=$color1 nowrap><B>#1:</B><SELECT id="db_1" name="databankID1" style="width:300px" onchange="updateDatabankSelection(1,this.value)">$databaseString</SELECT>&nbsp;&nbsp;<B>#2:</B><SELECT id="db_2" name="databankID2" style="width:300px" onchange="updateDatabankSelection(2,this.value)" disabled>$databaseString</SELECT>&nbsp;&nbsp;<B>#3:</B><SELECT id="db_3" name="databankID3" style="width:300px" disabled>$databaseString</SELECT></TD><!--&nbsp;&nbsp;<B>Scan:</B><SELECT name="scanDB"><OPTION value="now">now</OPTION><OPTION value="later" selected>later</OPTION></SELECT>-->
	</TR>
	<TR>
		<TH align=right valign=top>Description :</TH>
		<TD bgcolor=$color1><TEXTAREA name='des' rows='1' cols='70'></textarea></TD>
	</TR>
	<TR>
		<TH align=right valign=top nowrap>&nbsp;Data filtering<SUP>*</SUP> :</TH>
		<TD bgcolor=$color1>
			&bull;<B>If decoy search, filter data at</B> <INPUT type="text" name="maxFDR" value="1" size="3" style="text-align:right"><B>% FDR using</B> <SELECT name="FDRalgo"><OPTION value="qvality" $disabQval>Qvality</OPTION><OPTION value="DTcount">DT count</OPTION><OPTION value="precomputed" selected>*Precomputed*</OPTION><OPTION value="mayu">Mayu</OPTION></SELECT> (Use value <B>&le;0</B> to disable this option).<BR>
			&bull;<B>Otherwise use </B><INPUT type="text" name="minScore" value="default" size="8" style="text-align:right">$defaultStrg
		</TD>
	</TR>
	<TR>
		<TH align=right valign=top nowrap>Max. rank<SUP>*</SUP> :</TH>
		<TD bgcolor=$color1><SELECT name='maxRank'>$maxRankOptionString</SELECT></TD>
	</TR>
	<TR>
		<TH align=right valign=top>Comments :</TH>
		<TD bgcolor=$color1><TEXTAREA name='comments' rows='2' cols='70'></textarea></TD>
	</TR>
|;

		}
		my $colspan=($decoyFile)? 1 : 2;
		print qq
|	<TR><TH colspan=$colspan bgcolor=$color2>
		<INPUT type="submit" name="save" value=" Proceed ">
		&nbsp;&nbsp;&nbsp;<INPUT type='reset' value="  Clear  ">
		&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="cancelAction();"></TH>
	</TR>
	</TABLE>
|;
		print "</TD></TR></TABLE>\n*Applies only to interpretations based on MS/MS\n" unless $decoyFile;
		print qq |
<SCRIPT type="text/javascript">
const anaNameSEL=document.getElementById('anaNameSel');
if (anaNameSEL) {
	anaNameSEL.selectedIndex=2; // -> search file
	selectAnaName(anaNameSEL.value);
}
</SCRIPT>
|;
	}
	elsif ($action eq 'clean') {
		#######################################
		####>Cleaning Form, File Selection<####
		#######################################
		print qq
|function checkCleanForm(myForm) {
	if (testCheckbox(myForm)==0) {
		alert('Error: No files selected !');
		return false;
	}
	if (confirm('Delete selected file(s) from server ?')) {return true;}
	return false;
}
function checkUncheckAll(checkStatus) {
	var anaBox=document.fileCleanForm.file;
	if (anaBox.length) { // more than 1 checkboxes
		for (i=0; i < anaBox.length; i++){anaBox[i].checked=checkStatus;}
	}
	else {anaBox.checked=checkStatus;} // Only 1 checkbox
}
document.getElementById('waitDiv').style.display='none';
</SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<CENTER>
<FONT class="title">Select Files to be Deleted</FONT>
<BR><FONT class="title2">(From $batchDirLabel)</FONT>
<BR><BR>
|;
		unless (scalar keys %filesInfo) {
			&printError("No files found!");
			exit;
		}
		print qq
|<FORM name="fileCleanForm" method="post" onsubmit="return checkCleanForm(this);">
<INPUT type="hidden" name="ID" value="$ID" >
<INPUT type="hidden" name="ACT" value="delete">
<TABLE align=center border=0 cellspacing=0 cellpadding=0 >
<TR><TD>
<TABLE cellspacing=0 cellpadding=0 >
<TR class="header" bgcolor="$bgColor">
	<TH nowrap width=120 align=left><INPUT type="checkbox" onclick="checkUncheckAll(this.checked)">Search file</TH>
	<TH nowrap width=120 align=center>MS file</TH>
	<TH nowrap width=120>Analysis Type</TH>
	<TH width=140>Databank(s)<BR>Taxonomy</TH>
	<TH width=350>Title</TH>
</TR>
|;
		foreach my $dataFile (sort {lc($a) cmp lc($b)} keys %filesInfo) {
			$bgColor=($bgColor eq $color1)? $color2 : $color1;
			if (!$filesInfo{$dataFile}{'OK'}) {
				print qq
|<TR class="list" bgcolor=$bgColor valign=middle>
	<TH align=left nowrap><INPUT type="checkbox" name="file" value="$dataFile">$dataFile</TH>
	<TH colspan=4 align=left>&nbsp;<FONT style="color:#DD0000">Unsupported file format or version</FONT></TH>
</TR>
|;
				next;
			}
			elsif ($filesInfo{$dataFile}{'FILE_FORMAT'} eq 'PROTEOME_DISCOVER.MSF') {
				print qq
|<TR bgcolor=$bgColor valign=middle>
	<TD nowrap colspan=5><FIELDSET><LEGEND><INPUT type="checkbox" name="file" value="$dataFile"><B><U>$dataFile</U></B></LEGEND>
	<TABLE>
|;
				foreach my $processingNodeNumber (sort{$a<=>$b} keys %{$filesInfo{$dataFile}{'SEARCHES'}}) {
					print qq
|		<TR bgcolor=$bgColor valign=top>
		<TH nowrap align=left width=135>Request $processingNodeNumber</TH><TD>&nbsp;<B>$filesInfo{$dataFile}{searches}{$processingNodeNumber}{DESCRIPTOR}</B></TD>
		<TD><DIV id="$dataFile:$processingNodeNumber:show" style="display:none"><FIELDSET><LEGEND><INPUT type="button" value="Hide Search Parameters" onclick="showHideSearchParam('$dataFile:$processingNodeNumber','hide')" style="font-size:11px"/></LEGEND>$filesInfo{$dataFile}{SEARCHES}{$processingNodeNumber}{MOUSE_STRING}</FIELDSET></DIV>
			<DIV id="$dataFile:$processingNodeNumber:hide"><INPUT type="button" value="Show Search Parameters" onclick="showHideSearchParam('$dataFile:$processingNodeNumber','show')" style="font-size:11px"/></DIV>
		</TD>
		</TR>
|;
				}
				print "\t\t</TABLE></FIELDSET></TD>";
			}
			else {
				my $dbStrg=join('&bull;',@{$filesInfo{$dataFile}{'DB'}});
				print qq
|<TR class="list" bgcolor=$bgColor valign=middle>
	<TH align=left nowrap><INPUT type="checkbox" name="file" value="$dataFile">$dataFile</TH>
	<TD align=center>&nbsp;$filesInfo{$dataFile}{peakFile}&nbsp;</TD>
	<TD align=center>&nbsp;$filesInfo{$dataFile}{type}&nbsp;</TD>
	<TD align=center>$dbStrg<BR>$filesInfo{$dataFile}{taxonomy}</TD>
	<TD>&nbsp;$filesInfo{$dataFile}{title}&nbsp;</TD>
</TR>
|;
			}
		}
		print qq
|</TABLE></TD></TR>
<TR><TD><TABLE width=100%>
	<TR><TH bgcolor=$color2>
	<INPUT type="submit" name="delete" value=" Delete ">&nbsp &nbsp &nbsp
	<INPUT type="button" value="Cancel" onclick="cancelAction();">
	</TH></TR>
</TABLE></TD></TR>
</TABLE>
</FORM>
</CENTER>
|;
	}
	#####################
	####>Ending HTML<####
	#####################
	print qq
|</FORM>
</BODY>
</HTML>
|;
}


elsif ($action eq 'start') {

	##########################
	####>Connecting to DB<####
	##########################
	my $dbh=&promsConfig::dbConnect;
	my ($userStatus,$mascotIDs)=$dbh->selectrow_array("SELECT USER_STATUS,MASCOT_IDS FROM USER_LIST WHERE ID_USER='$userID'");
	$dbh->disconnect;

	my @userDirList=($userID);
	my $okMascot=0;
	if ($userStatus eq 'bio') {
		$mascotIDs='' unless defined $mascotIDs;
		$mascotIDs=~s/^,//; # in case no Mascot userID (PP 23/08/11)
		$okMascot=1 if $mascotIDs=~/^\d/;
	}
	else {
		$okMascot=1;
		$mascotIDs='0';
		opendir (DIR, $batchDir);
		while (my $thisDir = readdir (DIR)) {
			next if $thisDir eq $userID;
			next if $thisDir=~/^\./;
			push @userDirList, $thisDir if -d "$batchDir/$thisDir";
		}
		close DIR;
	}


	#########################################
	####>Starting Form, Option Selection<####
	#########################################
	my $today=strftime("%Y%m%d",localtime);
	####>Starting HTML<####
	print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);

	print qq |
<HTML>
<HEAD>
<TITLE>Search file sources</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT type="text/javascript">
|;
	&promsMod::browseDirectory_JavaScript if $promsPath{'shared'};
	print qq |
var fileSource;
function cancelAction() {
	//top.promsFrame.selectedAction='summary'; // set default action to 'summary'
	top.promsFrame.optionFrame.selectOption();
}
function selectSource(source) {
	fileSource=source;
	var myForm=document.fileAccessForm;
	document.getElementById('submitButton').disabled=false;

	/* Disable all file sources */
	myForm.userDir.disabled=true;
	if (myForm.newDir) {myForm.newDir.disabled=false;} // not defined if user is bio
	if (document.getElementById('sharedDirDIV')) {document.getElementById('sharedDirDIV').style.display='none';} // sharedDIr unabled
	myForm.uplArch.disabled=true;
	for (let i=1;i<=$maxUpFiles;i++) {
		document.getElementById('uploaded_file_'+i).disabled=true;
	}
	if ($okMascot) {document.getElementById('mascotTable').style.display='none';}
	
	/* Activate selected source */	
	if (fileSource=='UseUserDirectory') { // User directory
		myForm.userDir.disabled=false;
	}
	else if (fileSource=='UseServerDirectory') { // Custom directory
		myForm.newDir.disabled=false;
	}
	else if (fileSource=='UseSharedDirectory') { // Shared directory
		document.getElementById('sharedDirDIV').style.display='';
	}
	else if (fileSource=='UseZipFile') { // Upload archive
		myForm.uplArch.disabled=false;
	}
	else if (fileSource=='UseUploadedFiles') { // Upload multiple files
		for (let i=1;i<=$maxUpFiles;i++) {
			document.getElementById('uploaded_file_'+i).disabled=false;
		}
	}
	else if (fileSource=='UseMascot') { // Mascot server
		document.getElementById('mascotTable').style.display='';
		if (myForm.mascotServer.length==2) { // only 1 Mascot server => autoselect
			myForm.mascotServer.selectedIndex=1;
			updateLogList(myForm.mascotServer.value);
		}
	}
}
function displayNextFile(f) {
	// check for duplicates
	var fileValue=document.getElementById('uploaded_file_'+f).value;
	for (let i=1;i<=$maxUpFiles;i++) {
		if (i != f && document.getElementById('uploaded_file_'+i).value && document.getElementById('uploaded_file_'+i).value==fileValue) {
			alert('ERROR: This file (#'+f+') is already used at position #'+i+' !');
			document.getElementById('uploaded_file_'+f).value='';
			return;
		}
	}
	if (f==$maxUpFiles) {return;}
	f++;
	document.getElementById('tabFile_'+f).style.display='block';
}
var searchLogs={}; // updated after FORM
function updateLogList(mascotServer) {
	var selSearchLogs=document.getElementById('searchLogs');
	selSearchLogs.options.length=0;
	if (searchLogs[mascotServer]) {
		for (let i=0;i<searchLogs[mascotServer].length;i++) {
			selSearchLogs.options[i]=new Option(searchLogs[mascotServer][i],searchLogs[mascotServer][i]);
		}
		selSearchLogs.disabled=false;
		if (searchLogs[mascotServer].length==1) {selSearchLogs.selectedIndex=0;} // if only 1 search log => autoselect
	}
	else {
		selSearchLogs.disabled=true;
	}
}
function searchMascotServer() {
	document.getElementById('infoSpan').innerHTML='';
	var myForm=document.fileAccessForm;
	if (!myForm.mascotServer.value) {alert("ERROR: No Mascot server selected !"); return;}
	if (!myForm.searchLogs.value) {alert("ERROR: No Mascot log files selected !"); return;}
	//Date range check
	if (myForm.startDate.value) {
		if (!myForm.endDate.value) {
			alert('ERROR: No end date provided !');
			myForm.endDate.focus();
			return;
		}
		if (myForm.startDate.value.match(/\\D/) \|\| myForm.startDate.value*1 < 20000101 \|\| myForm.startDate.value*1 > $today) {
			alert("ERROR: Start date is not valid !");
			myForm.startDate.focus();
			return;
		}
	}
	if (myForm.endDate.value) {
		if (!myForm.startDate.value) {
			alert('ERROR: No start date provided !');
			myForm.startDate.focus();
			return;
		}
		if (myForm.endDate.value.match(/\\D/) \|\| myForm.endDate.value*1 < 20000101 \|\| myForm.endDate.value*1 > $today) {
			alert("ERROR: End date is not valid !");
			myForm.endDate.focus();
			return;
		}
	}
	if (myForm.startDate.value && myForm.endDate.value && myForm.startDate.value*1 > myForm.endDate.value*1) {
		alert('ERROR: Date range is not valid !');
		return;
	}
	//Job range check
	if (myForm.startJob.value) {
		if (!myForm.endJob.value) {
			alert('ERROR: No end job number provided !');
			myForm.endJob.focus();
			return;
		}
		if (myForm.startJob.value.match(/\\D/)) {
			alert("ERROR: Start job number is not valid !");
			myForm.startJob.focus();
			return;
		}
	}
	if (myForm.endJob.value) {
		if (!myForm.startJob.value) {
			alert('ERROR: No start job number provided !');
			myForm.startJob.focus();
			return;
		}
		if (myForm.endJob.value.match(/\\D/)) {
			alert("ERROR: End job number is not valid !");
			myForm.endJob.focus();
			return;
		}
	}
	if (myForm.startJob.value && myForm.endJob.value && myForm.startJob.value > myForm.endJob.value) {
		alert('ERROR: Job number range is not valid !');
		return;
	}
	//Search string

	//Submit
	if (!myForm.startDate.value && !myForm.startJob.value && !myForm.searchStrg.value && !confirm('No filters set. Proceed anyway?')) {return;}
	myForm.ACT.value='UseMascot';
	myForm.mascotAction.value='list';
	document.getElementById('mascotFrame').style.display='block';
	myForm.target='mascotFrame';
	myForm.submit();
}
function checkStartForm(myForm) {
	myForm.ACT.value=fileSource;
	myForm.target='resultFrame'; // in case changed by previous 'UseMascot'
	if (fileSource=='UseServerDirectory') {
		if (!myForm.newDir.value) {alert("ERROR: No server directory provided !"); return false;}
	}
	else if (fileSource=='UseSharedDirectory') {
		var failed=false;
		if (!myForm.sharedDirFiles) {failed=true;} // no file found
		else if (myForm.sharedDirFiles.length) { // multiple files found
			failed=true;
			for (let i=0; i<myForm.sharedDirFiles.length; i++) {
				if (myForm.sharedDirFiles[i].checked) {
					failed=false;
					break;
				}
			}
		}
		else if (!myForm.sharedDirFiles.checked) {failed=true;} // single file found
		if (failed) {alert("ERROR: No file selected from Shared data directory !"); return false;}
	}
	else if (fileSource=='UseZipFile') {
		if (!myForm.uplArch.value) {alert("ERROR: No archive file provided !"); return false;}
	}
	else if (fileSource=='UseUploadedFiles') {
		for (let i=1;i<=$maxUpFiles;i++) {
			if (document.getElementById('uploaded_file_'+i).value) {
				return true;
			}
		}
		alert("ERROR: No local file(s) provided !");
		return false;
	}
	else if (fileSource=='UseMascot') {
		if (myForm.mascotAction.value != 'get') {
			searchMascotServer();
			return false; // do not submit main form
		}
		//'get' mode: Fetch list of checked dat files
		var chkList=mascotFrame.getCheckedItems();
		if (chkList.match('folder')) { // transform 'folder' into 'file'
			var newList=[];
			var chkItems=chkList.split('+');
			for (let I=0;I<chkItems.length;I++) {
				// folders
				if (chkItems[I].match('folder')) {
					var itemBoxes=mascotFrame.document.treeForm.checkIndex;
					var chkFolders=chkItems[I].split(/[:,]/);
					for (let F=1;F<chkFolders.length;F++) {
						var tabIndex=mascotFrame.getItemIndex('folder:'+chkFolders[F]);
						for (let i=tabIndex+1;i<mascotFrame.tableArray.length;i++) {
							if (mascotFrame.tableArray[i][2]<=mascotFrame.tableArray[tabIndex][2]) {break;} // brother or parent is reached
							if (itemBoxes[i].disabled) {continue;}
							var branch=mascotFrame.tableArray[i][4].split(':');
							newList.push(branch[1]);
						}
					}
				}
				// files
				else {
					var chkFiles=chkItems[I].split(/[:,]/);
					for (let f=1;f<chkFiles.length;f++) {
					newList.push(chkFiles[f]);}
				}
			}
			myForm.datFiles.value='file:'+newList.join(',');
		}
		else {myForm.datFiles.value=chkList;}
		if (!myForm.datFiles.value) {alert("ERROR: No files selected !"); return false;}
		//else {alert(myForm.datFiles.value);}
		myForm.target='mascotFrame';
	}
//return false;
	return true;
}
function updateDisableStatus() { // needed after '<< BACK' button
	var radioSource=document.fileAccessForm.selSource;
	for (let i=0;i<radioSource.length;i++) {
		if (radioSource[i].checked) {
			selectSource(radioSource[i].value);
			break;
		}
	}
}
</SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif' onload="updateDisableStatus();">
<CENTER>
|;
	my $titleStrg='Select a File Source for Import of ';
	$titleStrg.='Decoy Data into ' if $decoyFile;
	$titleStrg.='Multiple Analyses';
	print qq
|<FONT class="title">$titleStrg</FONT>
<BR><BR>
<FORM name="fileAccessForm" method="post" onsubmit="return checkStartForm(this)" enctype="multipart/form-data">
<INPUT type="hidden" name="ID" value="$ID">
<INPUT type="hidden" name="ACT" value="preview">
<INPUT type="hidden" name="decoy" value="$decoyFile">
<INPUT type="hidden" name="mascotIDs" value="$mascotIDs">
<TABLE><TR><TD bgcolor=$color2>
	<TABLE>
	<TR>
		<TD colspan=2>&nbsp<FONT class="title2">Possible sources :</FONT></TD>
		<TD></TD>
	</TR>
	<TR>
		<TH><INPUT type="radio" name="selSource" id="selSource_uds" value="UseUserDirectory" onclick="selectSource(this.value);"></TH>
		<TH align=right><LABEL for="selSource_uds">A user directory on server</LABEL> :</TH>
		<TH bgcolor=$color1 align=left width=610px><SELECT name="userDir" style="font-weight:bold" disabled>
|;
	foreach my $userDir (sort{lc($a) cmp lc($b)} @userDirList) {
		print "<OPTION value=\"$userDir\"";
		print ' selected' if $userDir eq $userID;
		print ">$userDir</OPTION>\n";
	}
	print qq
|</SELECT>
		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="submit" name="clean" value="Clean My Directory" style="font-weight:bold;color:#DD0000;width:170px" onclick="fileSource='clean'">
		</TH>
	</TR>
|;
	##>shared directory<##
	if ($promsPath{'shared'}) {
		print qq
|	<TR valign=top>
		<TH><INPUT type="radio" name="selSource" id="selSource_sdd" value="UseSharedDirectory" onclick="selectSource(this.value);"></TH>
		<TH align=right><LABEL for="selSource_sdd">Shared data directory<LABEL> :</TH>
		<TD bgcolor=$color1><DIV id="sharedDirDIV" style="width:600px;max-height:300px;overflow:auto;display:none">
|;
		&promsMod::browseDirectory_getFiles($promsPath{'shared'},'sharedDirFiles',{fileMatch=>qr/\.(dat|msf|xml|zip|gz)$/i}); # List of handled file extensions
		print qq
|		</DIV></TD>
	</TR>	
|;			
	}
	if ($userStatus eq 'bioinfo') {
		print qq
|	<TR valign=top>
		<TH><INPUT type="radio" name="selSource" id="selSource_ads" value="UseServerDirectory" onclick="selectSource(this.value);"></TH>
		<TH align=right><LABEL for="selSource_ads">Any directory on server<SUP>*</SUP></LABEL> :</TH>
		<TD bgcolor=$color1><INPUT type="text" name="newDir" value="" size=70 disabled><BR><FONT style="font-style:italic"><SUP>*</SUP> for bioinformatician only</FONT>
		</TD>
	</TR>
|;
	}
	my $disabMascotFlag=($okMascot)? '' : ' disabled';
	print qq
|	<TR valign=top>
		<TH><INPUT type="radio" name="selSource" id="selSource_ms" value="UseMascot" onclick="selectSource(this.value);"$disabMascotFlag></TH>
		<TH align=right><LABEL for="selSource_ms">Mascot server</LABEL> :</TH>
		<TD bgcolor=$color1>
|;
		if ($okMascot) {
			print qq
|<TABLE cellpadding=0 id="mascotTable" width=100% style="display:none">
			<TR><TD valign=top><TABLE border=0 width=100%><TR>
					<TH align=left valign=top nowrap>Server :<SELECT name="mascotServer" onchange="updateLogList(this.value)"><OPTION value="">-=Select=-</OPTION>
|;
		foreach my $server (sort keys %mascotServers) {print "<OPTION value=\"$server\">$server</OPTION>";}
		print qq
|</SELECT><BR><BR><BR>Search filters :
					</TH>
					<TH align=right valign=top nowrap>&nbsp;&nbsp;Log files :</TH>
					<TD valign=top><SELECT multiple name="searchLogs" id="searchLogs" style="width:170px;height:60px" disabled></SELECT></TD>
				</TR>
				<TR>
					<TH align=left valign=top colspan=3 style="height:70px">
&nbsp;&nbsp;-Date range : from<INPUT type="text" name="startDate" style="width:80px" value=""/> to<INPUT type="text" name="endDate" style="width:80px" value=""/> <FONT style="font-size:9px">(yyyymmdd)</FONT><BR>
&nbsp;&nbsp;-Job range : from F<INPUT type="text" name="startJob" style="width:70px" value=""/>.dat to F<INPUT type="text" name="endJob" style="width:70px" value=""/>.dat<BR>
&nbsp;&nbsp;-Title contains :<INPUT type="text" name="searchStrg" style="width:200px" value=""/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="button" value=" Search " onclick="searchMascotServer()"/>
<INPUT type="hidden" name="mascotAction" value="">
<INPUT type="hidden" name="datFiles" value="">
					</TH>
				</TR></TABLE>
				</TD>
				<TD align=right rowspan=2 width=240><IFRAME name="mascotFrame" id="mascotFrame" width=230 height=300 frameborder=0 style="display:none"></IFRAME></TD>
			</TR>
			<TR>
				<TD valign=top><SPAN id="infoSpan"></SPAN></TD>
			</TR></TABLE>
|;
		}
		else { # bio with no Mascot IDs
			print qq|<FONT style="font-weight:bold;color:#DD0000;">No access</FONT> (No Mascot IDs recorded)|;
		}
		print qq
|		</TD>
	</TR>
	<TR>
		<TH><INPUT type="radio" name="selSource" id="selSource_uza" value="UseZipFile" onclick="selectSource(this.value);"></TH>
		<TH align=right><LABEL for="selSource_uza">Upload Zip archive</LABEL> :</TH>
		<TD bgcolor=$color1><INPUT type="file" name="uplArch" size=80 disabled></TD>
	</TR>
	<TR>
		<TH valign=top><INPUT type="radio" name="selSource" id="selSource_umf" value="UseUploadedFiles" onclick="selectSource(this.value);"></TH>
		<TH align=right valign=top><LABEL for="selSource_umf">Upload multiple files</LABEL> :</TH>
		<TD bgcolor=$color1>
|;
	print "\t\t\t<TABLE id=\"tabFile_1\"><TR><TH align=right width=70 nowrap>File #1 :</TH><TD><INPUT type=\"file\" id=\"uploaded_file_1\" name=\"uploaded_file_1\" size=80 onchange=\"displayNextFile(1)\" disabled></TD></TR></TABLE>\n";
	foreach my $i (2..$maxUpFiles) {
		print "\t\t\t<TABLE id=\"tabFile_$i\" style=\"display:none\"><TR><TH align=right width=70 nowrap>File #$i :</TH><TD><INPUT type=\"file\" id=\"uploaded_file_$i\" name=\"uploaded_file_$i\" size=80 onchange=\"displayNextFile($i)\" disabled></TD></TR></TABLE>\n";
	}
	print qq
|		</TD>
	</TR>
	<TR><TH colspan=3>
		<INPUT type="submit" id="submitButton" value=" Proceed " disabled>&nbsp&nbsp&nbsp
		<INPUT type="button" value=" Cancel " onclick="cancelAction();">
	</TH></TR>
	</TABLE></TR></TD>
</TABLE>
</FORM>
</CENTER>
<SCRIPT type="text/javascript">
|;
	####>Fetching list of search logs for all Mascot servers<####
	foreach my $mascotServer (keys %mascotServers) {
		my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
		$agent->timeout(360);
		if ($mascotServers{$mascotServer}{proxy}) { # proxy settings
			if ($mascotServers{$mascotServer}{proxy} eq 'no') {$agent->no_proxy($mascotServers{$mascotServer}{url});}
			else {$agent->proxy('http', $mascotServers{$mascotServer}{proxy});}
		}
		else {$agent->env_proxy;}
		my $response = $agent->get("$mascotServers{$mascotServer}{url}/cgi/myproms4datFiles.pl?ACT=log");
		while (my $wait = $response->header('Retry-After')) {
			sleep $wait;
			$response = $agent->get($response->base);
        }
		my @resultLines;
        if ($response->is_success) {
			@resultLines = split("\n",$response->content);
		}
		my $logString='';
		foreach (@resultLines) {
			next if /^#/;
			chomp;
			if (/\.log\Z/) {
				$logString.=',' if $logString;
				$logString.="'$_'";
			}
		}
		print "searchLogs['$mascotServer']=[$logString];\n" if $logString;
	}
	print qq
|</SCRIPT>
</BODY>
</HTML>
|;

}

elsif ($action eq 'delete') {
	my @files = param('file');
	my $path = "$batchDir/$userID";

	print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<CENTER><FONT class="title">Deleting Files...</FONT></CENTER>
<BR><BR>
<FONT class="title3">
|;
	foreach my $file (@files) {
		unlink "$path/$file";
		print "&nbsp;$file... deleted.<BR>\n";
	}
	print "<BR><BR>All files have been deleted.</FONT>\n";

	sleep 5;
	print qq
|<SCRIPT type="text/javascript">
	//top.promsFrame.selectedAction='summary'; // set default action to 'summary'
	top.promsFrame.optionFrame.selectOption();
</SCRIPT>
</BODY>
</HTML>
|;
}

sub printError {
	my ($text)=@_; #$header,
#	if ($header) {
#		print header(-charset=>'utf-8');
#		print qq
#|<HTML>
#<HEAD>
#<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
#<SCRIPT type="text/javascript">
#function cancelAction() {
#	//top.promsFrame.selectedAction='summary'; // set default action to 'summary'
#	top.promsFrame.optionFrame.selectOption();
#}
#</SCRIPT>
#</HEAD>
#<BODY background='$promsPath{images}/bgProMS.gif'>
#<CENTER>
#<BR>
#|;
#	}
	print qq
|<FONT class="title2" style="color:#DD0000">$text</FONT>
<BR><BR>
<INPUT type="button" value="<< Back" onclick="history.back()">&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="cancelAction();">
</CENTER>
</BODY>
</HTML>
|;
	exit;
}

#####---tag INSERTION---####
# Update the script so has to have a popup (modify ending html and add the &promsMod::popupInfo())
############################################################################################################
####>Display the search parameters input for MSF (ProteomeDiscover) for each search (Mascot or Sequest)<####
############################################################################################################
sub displaySearchParameters{
    my ($dbsqlite,$currentSearch,$rawfilenumber,$xmlData)=@_; #,$desc
    #my $req=$dbsqlite->prepare("SELECT ParameterName, ValueDisplayString, ParameterValue FROM ProcessingNodeParameters WHERE ProcessingNodeParameters.ProcessingNodeNumber=$currentSearch OR ProcessingNodeParameters.ProcessingNodeNumber=$rawfilenumber;");
    #my $req=$dbsqlite->prepare("SELECT ParameterName, ValueDisplayString, ParameterValue FROM ProcessingNodeParameters WHERE ProcessingNodeNumber=$currentSearch OR ProcessingNodeNumber=$rawfilenumber;");
    ###> Add ProcessingNodeParentNumber to get activation type -> CID, ETD, HCD... if available

    my ($ftol,$ptol,$miscl,$dbs,$enz,$activation);
    my ($static,$dynam,$activationString,$massFilterString)=("","","","");
    if ($xmlData) { # PD 2.0
		foreach my $processingNodeParameters (sort{$a->{ParentProcessingNodeNumber} cmp $b->{ParentProcessingNodeNumber}} @{$xmlData->{WorkflowTree}->{WorkflowNode}}) {
			next unless $processingNodeParameters->{ProcessingNodeNumber} == $currentSearch;
			foreach my $processingNodeParameter (@{$processingNodeParameters->{ProcessingNodeParameters}->{ProcessingNodeParameter}}){
				if($processingNodeParameter->{Name} eq 'FragmentTolerance'){$ftol=$processingNodeParameter->{content};}
				elsif($processingNodeParameter->{Name} eq 'PeptideTolerance'){$ptol=$processingNodeParameter->{content};}
				elsif($processingNodeParameter->{Name} eq 'MaxMissedCleavages' || $processingNodeParameter->{Name} eq 'MissedCleavages'){$miscl=$processingNodeParameter->{content};}
				elsif ($processingNodeParameter->{Name} eq 'FastaDatabase' || $processingNodeParameter->{Name} eq 'Database'){
					$dbs=$processingNodeParameter->{content};
				}
				elsif($processingNodeParameter->{IntendedPurpose} eq 'CleavageReagent'){
					my $enzymeXML=new XML::Simple();
					my $enzymeDesc=$enzymeXML->XMLin($processingNodeParameter->{content}); # <Enzyme Version="1" Name="Trypsin" CleavageSites="KR" CleavageInhibitors="P" Offset="1" CleavageSpecificity="SpecificAtBothEnds" />
					$enz=$enzymeDesc->{Name};
				}
				elsif($processingNodeParameter->{Name} eq 'ActivationTypeFilter'){
					$activationString="<TR><TH align=right nowrap valign=top>Activation Type Filter:</TH><TD width=100%>&nbsp;$processingNodeParameter->{content}</TD></TR>";
				}
				elsif($processingNodeParameter->{Name} eq 'MassAnalyzerFilter'){
					$massFilterString="\n<TR><TH align=right nowrap valign=top>Mass Analyser Filter:</TH><TD width=100%>&nbsp;$processingNodeParameter->{content}</TD></TR>";
				}
				elsif(($processingNodeParameter->{IntendedPurpose} =~ "Dynamic(Terminal)?Modification") && $processingNodeParameter->{IsValueSet} eq "True"){
					my $dynmodXML=new XML::Simple();
					my $dynmodDesc=$dynmodXML->XMLin($processingNodeParameter->{content}); #<Modification Version="2" AminoAcids="C" Name="Carbamidomethyl" Abbreviation="Carbamidomethyl" ID="8" UnimodAccession="4" DeltaMass="57.02146" DeltaAverageMass="57.05130" IsSubstitution="False" LeavingGroup="" Substitution="H(3) C(2) N O" PositionType="Any" />
					$dynam.=', ' if $dynam;
					$dynam.=$dynmodDesc->{Name};
				}
				elsif(($processingNodeParameter->{IntendedPurpose} =~ "Static(Terminal)?Modification") && $processingNodeParameter->{IsValueSet} eq "True"){
					my $statmodXML=new XML::Simple();
					my $statmodDesc=$statmodXML->XMLin($processingNodeParameter->{content}); #<Modification Version="2" AminoAcids="X" Name="MappingL" Abbreviation="L" ID="1167" UnimodAccession="-1" DeltaMass="113.08406" DeltaAverageMass="113.15890" IsSubstitution="True" LeavingGroup="" Substitution="C6H11NO" PositionType="Any" />
					$static.=', ' if $static;
					$static.=$statmodDesc->{Name} unless $statmodDesc->{Name} eq "MascotXValue";
				}
			}
		}
    }
    else {
		my $req=$dbsqlite->prepare("SELECT ParameterName, ValueDisplayString, ParameterValue FROM ProcessingNodeParameters WHERE ProcessingNodeNumber=$currentSearch OR ProcessingNodeNumber=$rawfilenumber OR ProcessingNodeNumber=(SELECT ProcessingNodeParentNumber FROM ProcessingNodes WHERE ProcessingNodeNumber=$currentSearch)");
		$req->execute;
		while(my ($paramN,$value,$paramV)=$req->fetchrow_array){
			$value=~s/\s*\Z//;
			if($paramN eq 'FragmentTolerance'){$ftol=$value;}
			elsif($paramN eq 'PeptideTolerance'){$ptol=$value;}
			elsif($paramN eq 'MaxMissedCleavages' || $paramN eq 'MissedCleavages'){$miscl=$value;}
			elsif($paramN eq 'FastaDatabase'){
				my $dbreq=$dbsqlite->prepare("SELECT VirtualFileName FROM FastaFiles WHERE FastaFileID=$paramV");
				$dbreq->execute;
				($dbs)=$dbreq->fetchrow_array;
				$dbs='' unless $dbs;
			}
			elsif($paramN eq 'Enzyme'){$enz=$value;}
			#elsif($paramN =~ /Stat/ && $paramN =~ /Mod/){
			elsif($paramN =~ /Stat/ && $paramN =~ /Mod_/){
				#$value =~ s/[0-9]*//g;
				#$value =~ s/\.//g;
				#$value =~ s/ \///g;
				#$value =~ s/\+//g;
				#$value =~ s/ Da //g;
				$value=~s/\/.+\(/\(/;
				$static.=', ' if $static;
				$static.=$value;
				#}elsif($paramN =~ /Dyn/ && $paramN =~ /Mod/){
			}
			elsif($paramN =~ /Dyn/ && $paramN =~ /Mod_/){
				#$value =~ s/[0-9]*//g;
				#$value =~ s/\.//g;
				#$value =~ s/ \///g;
				#$value =~ s/\+//g;
				#$value =~ s/ Da //g;
				$value=~s/\/.+\(/\(/;
				$dynam.=', ' if $dynam;
				$dynam.=$value;
			}
			elsif($paramN eq 'ActivationTypeFilter'){
				$activationString="<TR><TH align=right nowrap valign=top>Activation Type Filter:</TH><TD width=100%>&nbsp;$value</TD></TR>";
			}
			elsif($paramN eq 'MassAnalyzerFilter'){
				$massFilterString="\n<TR><TH align=right nowrap valign=top>Mass Analyser Filter:</TH><TD width=100%>&nbsp;$value</TD></TR>";
			}
		}
    }
    $static='None' unless $static;
    $dynam='None' unless $dynam;
    #my $mouseString="<TD colspan=6 align=right><FIELDSET><LEGEND><B>Search parameters:</B></LEGEND><TABLE border=0 cellpadding=0 width=100%>";
    my $mouseString=qq
|<TABLE border=0 cellpadding=0 width=100%>
<TR><TH align=right nowrap valign=top>Databank:</TH><TD width=100%>&nbsp;$dbs</TD></TR>
<TR><TH align=right nowrap valign=top>Enzyme:</TH><TD width=100%>&nbsp;$enz</TD></TR>
<TR><TH align=right nowrap valign=top>Enzyme miscleavages:</TH><TD width=100%>&nbsp;$miscl</TD></TR>
<TR><TH align=right nowrap valign=top>Fixed modifications:</TH><TD width=100%>&nbsp;$static</TD></TR>
<TR><TH align=right nowrap valign=top>Other modifications:</TH><TD width=100%>&nbsp;$dynam</TD></TR>
<TR><TH align=right nowrap valign=top>Fragment tolerance:</TH><TD width=100%>&nbsp;$ftol</TD></TR>
<TR><TH align=right nowrap valign=top>Peptide tolerance:</TH><TD width=100%>&nbsp;$ptol</TD></TR>
$activationString$massFilterString
</TABLE>
|; #<TR><TH align=right nowrap>Type of search:</TH><TD width=100%>&nbsp;$desc</TD></TR>

    return $mouseString;
}

sub updateWaitBox {
	print "<SCRIPT type=\"text/javascript\">document.getElementById('waitSpan').innerHTML+='.';</SCRIPT>\n";
}

####>Revision history<####
# 2.8.5 [BUGFIX] Fix wrong path to Mascot dat files when data_local_path is defined (PP 30/06/21)
# 2.8.4 [UX] Added Recently/Most used databanks in selection list & auto-select "search file" for Analysis name (PP 22/06/21)
# 2.8.3 Uses new &promsConfig::getMascotServers function (PP 25/06/19)
# 2.8.2 Modification of merge-file string (GA 08/01/19)
# 2.8.1 Creates batch/<userID> at start up if not exists (PP 04/12/18)
# 2.8.0 Added Shared directory as file source (requires promsMod 3.7.5+) & msf display update (PP 03/11/17)
# 2.7.1 Set default CGI file upload dir to $promsPath{tmp}/batch dir (PP 18/09/17)
# 2.7.0 Add "good.gif" in "mzXML file" column on the "Select Files to be Imported" form, if mzXML file is present in tmp dir (MLP 03/08/17)
# 2.6.9 Modification on the upload of mzXML file (used by Mayu algorithm) (MLP 04/05/17)
# 2.6.8 Add modifications in the "Select Files" form to select mzxml file used by Mayu (for X! Tandem and .dat files) (MLP 05/04/17)
# 2.6.7 Additional modifications for X! tandem files (MLP 16/12/16)
# 2.6.6 Added modifications to allow the import of X! tandem files (MLP 25/10/2016)
# 2.6.5 Retrieve Mascot *.pop files created in the cache folder when Percolator is run<BR> *.pop files are CSV files with score, q-value and PEP scores (GA 21/09/16)
# 2.6.4 More &updateWaitBox calls in case large merged PD files (PP 21/07/16)
# 2.6.3 Minor change to detect PD version >= 2.1 (PP 29/03/16)
# 2.6.2 Now uses $agent->request instead $agent->mirror to prevent server timeout on slow connection to Mascot (PP 23/03/16)
# 2.6.1 Check SQlite queries for 2.0 PD MSF (GA 12/06/15)
# 2.6.0 Minor changes to speed up PARAGON file scan (PP 18/03/15)
# 2.5.9 Modification for split-mode import (GA 03/03/15)
# 2.5.8 Minor modif to get Sequest HT or SEQUEST searches from PD1.4 (GA 28/11/14)
# 2.5.7 Modification for ProteinPilot XML for databank names (GA 29/10/14)
# 2.5.6 Minor update for ProteinPilot XML (GA 21/10/14)
# 2.5.5 Use of &promsMod::sortSmart for imported files listing (PP 20/10/14)
# 2.5.4 Checks for uninitialize $bds (PP 10/10/14)
# 2.5.3 Add Activation type and mass analyser filters in MSF files if the information is provided (GA 20/05/14)
# 2.5.2 Handles selection of percolator filtering performed by PD (PP 12/05/14)
# 2.5.1 Deletion allowed only if user directory is selected (PP 14/04/14)
# 2.5.0 Uses mkdir instead of make_path (PP 10/03/14)
# 2.4.9 Qvality vs DT count FDR estimations (PP 06/03/14)
# 2.4.8 Minor bug fix for double header print in &printError call (PP 06/09/13)
# 2.4.7 Forces Mascot dat files to be ordered by number (PP 07/06/13)
# 2.4.6 Update for new Proteome Discoverer 1.4 format (GA 27/05/13)
# 2.4.5 Uses LWP instead of wget (PP 13/05/13)
# 2.4.4 Make Paragon appear in the selection Threshold score line (GA 20/03/13)
# 2.4.3 Changed min date for Mascot dat file retrieval to 20000101 (PP 16/02/13)
# 2.4.2 Compatible with multi-databank searches: Mascot .dat only (PP 13/12/12)
# 2.4.1 Modification for PARAGON (GA 04/09/12)
# 2.4.0m Default rank selection is 1 (PP 11/09/12)
# 2.4.0 Added FDR-based score estimation option & removed DB Scan option (PP 03/08/12)
# 2.3.9 Minor improvement in symbolic link option to avoid multiple '/' in path (PP 03/08/12)
# 2.3.8 symbolic link option for .dat files from Mascot server (PP 03/05/12)
# 2.3.7 Minor change in sth for SQLite DB (ProcessingNodes table) to make it works for both PD 1.2 and 1.3 version (GA 27/04/12)
# 2.3.6 better handling of unrecognized file format (PP 19/03/12)
# 2.3.4 better handling of Modif in MSF especially SILAC eg. Label:13C(6) (PP 05/07/11)
# 2.3.3 Better handling of missformed MSF files for import & all MSF files for delete (PP 15/06/2011)
# 2.3.2 Compatibility with mixed (SQ) search (PP 16/04/2011)
# 2.3.1 Changes in bioinfo/mass/super user Server directory management (PP 15/03/2011)
