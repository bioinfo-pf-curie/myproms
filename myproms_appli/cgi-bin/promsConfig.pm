################################################################################
# promsConfig.pm           2.11.11D                                             #
# Authors: P. Poullet, G. Arras & V. Sabatet (Institut Curie)                  #
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

package promsConfig;

use CGI ':standard';
use DBI;
use Time::Piece;

use DBD::mysql;	# MySQL

require Exporter;

@ISA=qw(Exporter);
@EXPORT=qw();
@EXPORT_OK=qw();
$VERSION=1.00;
$CGI::LIST_CONTEXT_WARN=0; # Prevents security warning when CGI params are list of values

use strict;

####> DATABASE CONFIGURATION <####
my $host='DB_HOST'; # or IP address of MySQL server
my $port='DB_PORT'; # MySQL server port
my $dbName='DB_NAME'; # default: myproms
my $dbUser='DB_USER'; # default: myproms
my $password='DB_PASSWORD'; # default: myproms


####> WEB SERVER CONFIGURATION <####
sub getServerInfo {

	#>Checking user
	&checkUser unless $_[0];

	#>Server
	my %promsPath=('html'			=> '', # relative http path to html directory
					'cgi'			=> '/cgi-bin', # relative http path to cgi-bin directory
					'data'			=> 'DATA_DIR', # full unix path to data root directory
					'shared'			=> 'SHARED_DIR', # full unix path to shared directory (if any)
					'java'			=> 'JAVA_DIR', # full unix path to java (if any)
					'bin'			=> 'APPLI_DIR/bin', # full unix path to bin (e.g. phosphoRS)
					'R_scripts'		=> 'APPLI_DIR/scripts/R_scripts', # full unix path to R scripts
					'R'				=> 'R_DIR', # full path to R
					'qvality'		=> 'QVALITY_DIR', # full path to qvality (eg. /bioinfo/local/build/percolator_2_04/bin)
					'masschroq'		=> 'MASSCHROQ_DIR', # full path to masschroq (eg. /bioinfo/local/build/masschroq_1.2.1/bin)
					'tpp'           => 'TPP_DIR', # full path to TPP suite (TransProteomic Pipeline)
					'python'			=> 'PYTHON_DIR', # full path to default python
					'python2'		=> 'PYTHON2_DIR', # full path to python2
					'python3'		=> 'PYTHON3_DIR', # full path to python3
					'python_scripts'			=> 'APPLI_DIR/scripts/python',
					'msproteomicstools'		=> 'MSPROTEO_DIR', # full path to msproteomicstools (needed only if TPP is installed)
					'openms'			=> 'OPENMS_DIR', # full path to OpenMS tools
					'pyprophet'		=> 'PYPROPHET_DIR', # full path to pyprophet
					'http_proxy'	=> 'HTTP_PROXY' # URL:port (Leave empty if defined globally. Use 'no' if no proxy used)
					);

	my %autoPath=(
				'data_html'			=> "$promsPath{html}/data", # symlink to true data dir
				'images'			=> "$promsPath{html}/images",
				'banks'				=> "$promsPath{data}/banks",
				'swath_lib'			=> "$promsPath{data}/swath_lib",
				'swath_lib_html'	=> "$promsPath{html}/data/swath_lib",
				'valid'				=> "$promsPath{data}/validation",
				'logs'				=> "$promsPath{data}/logs",
                'metadata'          => "$promsPath{data}/metadata",
                'metadata_html'     => "$promsPath{html}/data/metadata",
				'peptide'			=> "$promsPath{data}/peptide_data",
				'spectrum'			=> "$promsPath{data}/spectrum_data",
				'gel'				=> "$promsPath{html}/data/gels",
				'gel_unix'			=> "$promsPath{data}/gels",
				'java_gel2d'		=> "$promsPath{html}/java/gel2d",
				'tmp'				=> "$promsPath{data}/tmp",
				'tmp_html'			=> "$promsPath{html}/data/tmp",
				'quantification'	=> "$promsPath{data}/quantification_data",
				'quanti_html'		=> "$promsPath{html}/data/quantification_data",
				'goa'				=> "$promsPath{data}/go/goa",
				'obo'				=> "$promsPath{data}/go/obo",
				'go_unix'			=> "$promsPath{data}/go/results",
				'go'				=> "$promsPath{html}/data/go/results",
				'genesets'			=> "$promsPath{data}/genesets",
				'gsea'				=> "$promsPath{data}/gsea",
				'gsea_html'			=> "$promsPath{html}/data/gsea",
				'explorAna'			=> "$promsPath{data}/exploratory_data",
				'explorAna_html'    => "$promsPath{html}/data/exploratory_data",
				'pathAna'			=> "$promsPath{data}/pathway_data"
				);

	@promsPath{keys %autoPath}=values %autoPath; # add %autoPath to %promsPath

	return (%promsPath);
}


####> CONNEXION TO DATABASE <####
sub dbConnect {

	#>Checking user
	&checkUser unless $_[0];
	if (!$host || !$dbName || !$dbUser) {
		die "ERROR: Database connection parameters not properly set. May be you forgot to call promsConfig::getServerInfo() first?";
	}

	#>DB connexion
	#my $dbh = DBI->connect("dbi:mysql:dbname=$dbName;host=$host;port=$port;mysql_init_command=SET NAMES utf8",$dbUser,$password,{RaiseError=>1,AutoCommit=>0}) || die "ERROR: ",$DBI::errstr;
	my $dbh = DBI->connect("dbi:mysql:dbname=$dbName;host=$host;port=$port",$dbUser,$password,
								{
									RaiseError			=> 1,
									AutoCommit			=> 0,
									mysql_init_command	=> "SET NAMES utf8, innodb_lock_wait_timeout = 600, sql_mode=''"
								}
							) || die "ERROR: ",$DBI::errstr;

	return $dbh;
}


####> CHECKING USER'S AUTHENTIFICATION <####
sub checkUser {
	my $badUser=0;
	if (cookie('myproms')) {
		my ($userID,$userKey)=(cookie('myproms')=~/USER=(.+),KEY=(.+)/);
		if ($userID && $userKey && $userKey eq crypt($userID,"myproms_$ENV{REMOTE_ADDR}")) {
			$ENV{'REMOTE_USER'}=$userID;
		}
		else {$badUser=1;}
	}
	else {$badUser=1;}

	##<Log again>##
	if ($badUser) {
		print header;
		print qq
|<HTML><HEAD>
<SCRIPT type="text/javascript">top.location="./login.cgi";</SCRIPT>
</HEAD></HTML>
|;
		exit;
	}
}

####> MAX PARALLEL JOBS (on web server) <####
sub getMaxParallelJobs {
	return 5;
}

####> CLUSTER SERVER <####
sub getClusterInfo {
	# my $clusterCommandPrefix='/path/to/myproms/singularity/image'; # TO BE EDITED!
	my %cluster; # Must be declared here even if only 1 cluster
	%cluster=(
### --Code below must be adapted to your own cluster. All keys of %cluster (pointing to values or subroutines) must be defined
		'on'=>0, # <------- set to 1 to activate
		# 'name'=>'CentOS',
		# 'maxCPUs'=>24,
		# 'maxJobs'=>50,
		# 'maxMemory'=>240, # Gb
		# 'singularityImage'=>$clusterCommandPrefix || '',
		'setSingularityImage' => sub {
			# $cluster{'singularityImage'}=$_[0] || '';
		},
		'buildCommand' => sub {
			# my ($runDir,$command)=@_; # $command: a command string or a bash file full name
			# my $userCommandFile;
			# if ($command !~ /[|>\n ]/ && -e $command) { # already a script (no 'bad' characters allowed in script path/name)
			# 	$userCommandFile=$command;
			# }
			# elsif ($command =~ /[;>]/) { # assume complex command => move to a bash file
			# 	$userCommandFile="$runDir/command.sh";
			# 	open (CMD,">$userCommandFile");
			# 	print CMD "#!/bin/bash\n$command\n";
			# 	close CMD;
			# 	my $modBash=0775;
			# 	chmod $modBash, $userCommandFile;
			# }
			# else { # assume simple command => run directly
			# 	$userCommandFile=$command;
			# }
			# return "$cluster{singularityImage} $userCommandFile";
		},
		'sendToCluster' => sub {
# 			my ($bashFile,$jobIdTag,$runDir)=@_;
# 			die "ERROR in sendToCluster: You must provide a full path for $bashFile" unless $bashFile=~/\//;
# 			$jobIdTag='' unless $jobIdTag;
# 			unless ($runDir) {($runDir=$bashFile)=~s/\/[^\/]+\Z//;}
# 			my $jobIdFile="$runDir/torqueID$jobIdTag.txt";
# 			my $connectFile="$runDir/ssh_connect.sh";
# 			open (CONN,">$connectFile");
# 			print CONN qq
# |#!/bin/bash
# /bioinfo/local/bin/sqsub_$serverInstance.sh $bashFile > $jobIdFile
# |;
# #ssh -o "UserKnownHostsFile=/dev/null" -o "StrictHostKeyChecking no" calcsub "$bashFile > $qsubFile" 2> $sshErrorFile
# 			close CONN;
# 			my $modBash=0775;
# 			chmod $modBash, $connectFile;
# 			system "$connectFile 2> $jobIdFile";

# 			my $jobID = `grep -Eo '[0-9]+' $jobIdFile | head -1 | tr -d '\n'`;
# 			return $jobID;
		},
		'checkError' => sub {
			# my ($errorFile,$refSkipErrors)=@_;
			# $refSkipErrors=[] unless $refSkipErrors;
			# my $userErrorStrg=(scalar @{$refSkipErrors})? " | grep -vi '".(join("' | grep -vi '",@{$refSkipErrors}))."'" : '';
			# my $error='';
			# if (-s $errorFile) {
			# 	# Remove:
			# 	# [0m[33mWARNING: Not mounting /bioinfo/projects_prod/myproms (already mounted in container)
			# 	# (Warning: )Permanently added...
			# 	# DIA InfluxDB error
			# 	$error=`grep -v "^\\[" $errorFile | grep -v "^perl: warning" | grep -vE "^    (LANG|LC_ALL|are supported)" | grep -v "Permanently " | grep -vi InfluxDB$userErrorStrg`;
			# 	$error=~s/\s//g unless $error=~/\S/; # remove empty and space-only lines if no text
			# }
			# return $error;
		},
		'killJob' => sub {
# 			my ($jobsIDRef, $waitForEnding)=@_;
			
# 			&promsMod::cleanDirectory("$promsPath{data}/tmp/delJobs",'15m'); # Remove previous deletion folders
			
# 			my $jobName = join('.', @{$jobsIDRef});
# 			my $runDir = "$promsPath{data}/tmp/delJobs/".join('.', @{$jobsIDRef});
# 			my $scriptPath = "$runDir/$jobName.sh";
# 			my %jobParams = (
# 				maxMem     => '20Mb',
# 				numCPUs    => 1,
# 				jobName    => $jobName."_script",
# 				outFile    => "PBS$jobName.txt",
# 				errorFile  => "PBSerror$jobName.txt",
# 			);
			
# 			system("mkdir -p $runDir");
# 			my $scriptContent = qq
# |#!/bin/bash

# ##resources
# #PBS -l mem=$jobParams{maxMem}
# #PBS -l nodes=$jobParams{numCPUs}:ppn=1
# #PBS -l walltime=00:30:00
# #PBS -q batch

# ##Information
# #PBS -N $jobParams{jobName}
# #PBS -m abe
# #PBS -o $runDir/$jobParams{outFile}
# #PBS -e $runDir/$jobParams{errorFile}

# ##Script Content
# |;
			
# 			foreach my $jobID (@{$jobsIDRef}) {
# 				$scriptContent .= "qdel $jobID\n";
# 			}
		
# 			open(FILE, '>', $scriptPath) || die "Error while opening $scriptPath\n";
# 			print FILE $scriptContent."\necho Done;";
# 			close FILE;
# 			chmod 0775, $scriptPath;
			
# 			system(($serverInstance =~ 'dev') ? "/bioinfo/local/bin/sqsub_dev.sh $scriptPath &" : "/bioinfo/local/bin/sqsub_prod.sh $scriptPath &"); # Start script
			
# 			# Remove generated files
# 			my $count = 0;
# 			my $maxCount = 60;
# 			while($waitForEnding && (!-e "$runDir/$jobParams{outFile}" || `cat $runDir/$jobParams{outFile}` !~ /Done/) && $count <= $maxCount) {
# 				$count++;
# 				sleep 2;
# 			}
		},
		'jobInfos' => sub {
			# my ($jobsIDRef, $realTime) = @_;
			# my @jobsID = @{$jobsIDRef};
			# my %jobInfos = ();
			# my $jobsLogDir = "/data/tmp/torque6_logs";
			# my $serverInstance = promsConfig::detectServer;
			# $jobsLogDir .= ($serverInstance =~ /prod/) ? '/http' : '/w-myproms-ud';
			# my $t = Time::Piece->new();
			
			# # Check if job is described in the exited job file summary
			# my $logFilesExistCmd = "find $jobsLogDir -name '".$t->year."*'";
			# my $getJobInfoExitedCmd = "grep -E '(".join('|', @jobsID).")\.torque6\.curie\.fr' $jobsLogDir/".$t->year."*.E";
			# my $jobsInfoStr = `$getJobInfoExitedCmd`;
			# if(`$logFilesExistCmd`) {
			# 	if($jobsInfoStr) {
			# 		foreach my $jobsInfoStr (split('\n', $jobsInfoStr)) {
			# 			my @allJobInfos = split(' ', $jobsInfoStr);
			# 			if($allJobInfos[1] =~ /([0-9]+).torque6.curie.fr/) {
			# 				my $jobID = $1;
			# 				foreach my $jobInfoField (@allJobInfos[2..$#allJobInfos]) {
			# 					my ($info, $value) = split('=', $jobInfoField);
			# 					$jobInfos{$jobID}{$info} = $value;
			# 				}
			# 				@jobsID = grep {!/$jobID/} @jobsID; # Remove job from queried list
			# 			}
			# 		}
			# 	}
			# }
			
			# # Try to find remaining jobs by using qstat on cluster
			# if(@jobsID && $realTime) {
			# 	my $jobName = "qstatJobs".join('.', @jobsID);
			# 	my $runDir = "$promsPath{data}/tmp/qstatJobs";
			# 	my $i = 0;
			# 	while(-s "$runDir".$i) {
			# 		$i++;
			# 	}
			# 	$runDir .= $i;
				
			# 	my $scriptPath = "$runDir/$jobName.sh";
			# 	my %jobParams = (
			# 		maxMem     => '20Mb',
			# 		numCPUs    => 1,
			# 		maxHours   => 1,
			# 		jobName    => $jobName."_script",
			# 		outFile    => "PBS$jobName.txt",
			# 		errorFile  => "PBSerror$jobName.txt",
			# 		noWatch    => '1',
			# 	);
				
			# 	mkdir $runDir;
			# 	my $scriptContent = qq |
			# 			#!/bin/bash
						
			# 			##resources
			# 			#PBS -l mem=$jobParams{maxMem}
			# 			#PBS -l nodes=$jobParams{numCPUs}:ppn=1
			# 			#PBS -l walltime=$jobParams{maxHours}:00:00
			# 			#PBS -q batch
						
			# 			##Information
			# 			#PBS -N $jobParams{jobName}
			# 			#PBS -m abe
			# 			#PBS -o $runDir/$jobParams{outFile}
			# 			#PBS -e $runDir/$jobParams{errorFile}
						
			# 			##Script Content
			# 	|;
				
			# 	$scriptContent .= "qstat -f ".join(' ', @jobsID)." > $runDir/qstat.out\n";
			
			# 	open(FILE, '>', $scriptPath) || die "Error while opening $scriptPath\n";
			# 	print FILE $scriptContent;
			# 	close FILE;
			# 	chmod 0775, $scriptPath;
				
			# 	system(($serverInstance =~ 'dev') ? "/bioinfo/local/bin/sqsub_dev.sh $scriptPath" : "/bioinfo/local/bin/sqsub_prod.sh $scriptPath"); # Start script
				
			# 	# Watch for task completion
			# 	my $nbWhile = 0;
			# 	while (!-s "$runDir/qstat.out" && $nbWhile < 120) {
			# 		sleep 1;
			# 		$nbWhile++;
			# 	}
				
			# 	if(-s "$runDir/qstat.out") {
			# 		open QSTAT, '<', "$runDir/qstat.out";
			# 		chomp(my @lines = <QSTAT>);
			# 		close QSTAT;
					
			# 		my $currentJobID = "";
			# 		foreach my $jobInfoField (@lines[0..$#lines]) {
			# 			chomp($jobInfoField);
			# 			if($jobInfoField =~ /([0-9]+).torque6.curie.fr/) {
			# 				$currentJobID = $1;
			# 				@jobsID = grep {!/$currentJobID/} @jobsID;
			# 			} else {
			# 				my ($info, $value) = split(' = ', $jobInfoField);
							
			# 				if($info && $value) {
			# 					$info =~ s/\s+//g;
			# 					if($info eq 'start_time') {
			# 						my $t = Time::Piece->strptime($value, "%a %h %d %H:%M:%S %Y");
			# 						$info = 'start';
			# 						$value = $t->epoch-3600;
			# 					}
			# 					$info = 'jobname' if($info eq 'Job_Name');
			# 					$jobInfos{$currentJobID}{$info} = $value
			# 				}
			# 			}
			# 		}
			# 	}
				
			# 	system("rm -rf $runDir");
			# }
			
			# # If jobs are remaining (only not started ones), check in the submitted jobs file summary
			# if(@jobsID) {
			# 	my $getJobInfoSubmitedCmd = "grep -E '(".join('|', @jobsID).")\.torque6\.curie\.fr' $jobsLogDir/".$t->year."*.S";
			# 	$jobsInfoStr = `$getJobInfoSubmitedCmd`;
			# 	if($jobsInfoStr) {
			# 		foreach my $jobsInfoStr (split('\n', $jobsInfoStr)) {
			# 			my @allJobInfos = split(' ', $jobsInfoStr);
			# 			if($allJobInfos[1] =~ /([0-9]+).torque6.curie.fr/) {
			# 				my $jobID = $1;
			# 				foreach my $jobInfoField (@allJobInfos[2..$#allJobInfos]) {
			# 					my ($info, $value) = split('=', $jobInfoField);
			# 					$jobInfos{$jobID}{$info} = $value;
			# 				}
							
			# 				@jobsID = grep {!/$jobID/} @jobsID; # Remove job from queried list
			# 			}
			# 		}
			# 	}
			# }
			
			# return %jobInfos;
		},
		'runJob' => sub {
# 			my ($runDir,$command,$refParam)=@_;
# 			# $refParam={
# 			#				-commandFile		(/path/to/)PBS command file. Optional, default: $runDir/PBScommand.sh
# 			#				-maxMem				Memory requested with unit. Optional, default: 1G
# 			#				-numCPUs			CPUs requested. Optional, default: 1
# 			#				-maxHours			Max job duration (in hours, no unit!). Optional, default: 24
# 			#				-jobName			Job name. Optional, default: myProMS_job
# 			#				-outFile			Cluster output file. Optional, default: $runDir/PBS.txt
# 			#				-errorFile			Cluster error output file. Optional, default: $runDir/PBSerror.txt
# 			#				-pbsRunDir			Cluster running directory. Optional, no default
# 			#				-commandBefore		Command to be run before main command. Optional, no default
# 			#				-commandAfter		Command to be run after main command. Optional, no default
# 			#				-jobIdTag		 	Extra tag for job id file. Optional, default ''
# 			#				-skipErrors			TODO: reference to array of error strings to be ignored in Cluster error output file (errorFile)
# 			#				-jobEndFlag			Flag to be echo(ed)	just before job ends, default '_JOB_END_', 'none' if nothing to echo (=>noWatch set to true!!!)
# 			#				-noWatch			Do not watch job: any non-(0/empty/undef) value, no default (undef: job is watched)
# 			# }

# 			##<Checking parameters>##
# 			$refParam={} unless $refParam;
# 			my $bashFile=$refParam->{'commandFile'} || 'PBScommand.sh';
# 			if ($bashFile !~ /\//) {$bashFile="$runDir/$bashFile";} # Full path required
# 			my $maxMem=$refParam->{'maxMem'} || '1Gb'; if ($maxMem=~/(\d+)Gb/ && $1 > $cluster{'maxMemory'}) {$maxMem=$cluster{'maxMemory'}.'Gb';}
# 			my $numCPUs=$refParam->{'numCPUs'} || 1; $numCPUs=$cluster{'maxCPUs'} if $numCPUs > $cluster{'maxCPUs'};
# 			my $maxHours=$refParam->{'maxHours'} || 24; $maxHours=~s/:.*//; # clean in case 'hh:mm:ss' format
# 			my $jobName=$refParam->{'jobName'} || 'myProMS_job';
# 			my $pbsFile=$refParam->{'outFile'} || 'PBS.txt';
# 			if ($pbsFile !~ /\//) {$pbsFile="$runDir/$pbsFile";} # Full path required
# 			my $pbsErrorFile=$refParam->{'errorFile'} || 'PBSerror.txt';
# 			if ($pbsErrorFile !~ /\//) {$pbsErrorFile="$runDir/$pbsErrorFile";} # Full path required
# 			my $clusterRunDirStrg='';
# 			if ($refParam->{'pbsRunDir'}) {
# 				if ($refParam->{'pbsRunDir'} !~ /\//) {$refParam->{'pbsRunDir'}="$runDir/".$refParam->{'pbsRunDir'};} # Full path required
# 				$clusterRunDirStrg='#PBS -d '.$refParam->{'pbsRunDir'};
# 			}
# 			my $commandBefore=$refParam->{'commandBefore'} || '';
# 			my $commandAfter=$refParam->{'commandAfter'} || '';
# 			my $jobIdTag=$refParam->{'jobIdTag'} || '';
# 			my $refSkipErrors=$refParam->{'skipErrors'} || [];
# 			my $jobEndFlag=$refParam->{'jobEndFlag'} || '_JOB_END_';
# 			my $noWatch=$refParam->{'noWatch'} || undef;
# 			$noWatch=1 if $jobEndFlag eq 'none';

# 			##<Building command (file)>##
# 			my $commandStrg=$cluster{'buildCommand'}->($runDir,$command);

# 			##<Writing PBS command file>##
# 			open (BASH,">$bashFile");
# 			print BASH qq
# |#!/bin/bash
# ##resources
# #PBS -l mem=$maxMem
# #PBS -l nodes=1:ppn=$numCPUs
# #PBS -l walltime=$maxHours:00:00
# #PBS -q batch
# ##Information
# #PBS -N $jobName
# #PBS -m abe
# #PBS -o $pbsFile
# #PBS -e $pbsErrorFile
# $clusterRunDirStrg

# ## Command(s)
# |;
# 			print BASH "export LC_ALL=\"C\"\n" if $cluster{'singularityImage'}=~/myproms_1\.1\.17\.img/;
# 			print BASH "$commandBefore\n" if $commandBefore;
# 			print BASH "$commandStrg\n";
# 			print BASH "$commandAfter\n" if $commandAfter;
# 			print BASH "\necho $jobEndFlag\n" if $jobEndFlag ne 'none';
# 			close BASH;
# 			chmod 0775, $bashFile;

# 			##<Sending job request to cluster>##
# 			my $jobID = $cluster{'sendToCluster'}->($bashFile,$jobIdTag,$runDir);

# 			##<Watch job>##
# 			my $pbsError='';
# 			unless ($noWatch) {
# 				sleep 10;
# 				my $nbWhile=0;
# 				my $maxNbWhile=$maxHours*60*2;
# 				while ((!-e $pbsFile || !`tail -3 $pbsFile | grep '$jobEndFlag'`) && !$pbsError) {
# 					if ($nbWhile > $maxNbWhile) {
# 						$pbsError="Aborting: Process is taking too long or died before completion";
# 						system "echo $pbsError >> $pbsErrorFile";
# 						$cluster{'killJob'}->([$jobID]);
# 						last;
# 					}
# 					sleep 30;

# 					##<Check for error
# 					$pbsError=$cluster{'checkError'}->($pbsErrorFile,$refSkipErrors);

# 					$nbWhile++;
# 				}
# 			}
			
# 			return ($pbsError, $pbsErrorFile, $jobID, $pbsFile);
		},

		#'getJobStatus' => sub {
		#	my ($jobID)=@_;
		#	my %jobStatus = (
		#		C => [1,'Job is completed after having run'],
		#		E => [1,'Job is exiting after having run'],
		#		H => [0,'Job is held'],
		#		Q => [0,'Job is queued, eligible to run or routed'],
		#		R => [0,'Job is running'],
		#		T => [0,'Job is being moved to new location'],
		#		W => [0,'Job is waiting for its execution time (-a option) to be reached'],
		#		S => [-1,'Job is suspended']
		#	);
		#	my $results=`qstat $jobID | grep ^$jobID`;
		#	my @returnInfo;
		#	if ($results) {
		#		my $status=(split(/\s+/,$results))[-2];
		#		@returnInfo=@{$jobStatus{$status}};
		#	}
		#	else {@returnInfo=(-1,'Job not found');}
		#	return @returnInfo;
		#},
		'path'		=>	{ # Singularity/docker
			# 'java'				=> '/usr/bin',
			# #'masschroq'			=> '/usr/local/bin',
			# 'masschroq'			=> '/usr/bin', # myproms_1.1.19-1.img!!!!!! && img 1.2.x
			# 'msproteomicstools'	=> '/usr/local/bin', #'/usr/local/lib/python2.7/dist-packages',
			# #'openms'			=> '/usr/local/bin', # img 1.1.x
			# 'openms'			=> '/usr/bin', # img 1.2.x
			# 'perl'				=> '/usr/local/bin',
			# 'pyprophet'			=> '/usr/local/bin',
			# 'python'			=> '/usr/bin', # same for python2 & python3
			# 'python2'			=> '/usr/bin', # python2
			# 'python3'			=> '/usr/bin', # python3
			# 'qvality'			=> '/usr/bin',
			# 'R'					=> '/usr/bin',
			# 'tpp'				=> '/usr/local/tpp/bin'
		}
	);

#	#<Cluster list># Multiple clusters? Defined specific %cluster for each of them & call &promsConfig::getClusterInfo('<clusterName>') in each script;
#	@{$cluster{'list'}}=('myCluster','my2ndCluster','...');

	return %cluster;
}

####> MASCOT SERVERS <####
sub getMascotServers {
	my %mascotServers=(
		##myServer=>{  Mascot server name
		##	url=>URL to Mascot server,
		##	proxy=>proxy string if any,
		##	data_local_path=>full path to Mascot dat files local directory if any,
		##	sequence_local_path=>full path to Mascot sequence files local directory if any,					
		##	link_files=>optional 0/1 flag for linking file
		##}
		## Config starts below --->
		#'MASCOT_SERVER_NAME'=>{url=>'MASCOT_SERVER_URL',proxy=>'MASCOT_SERVER_PROXY',data_local_path=>'MASCOT_DATA_DIR',sequence_local_path=>'MASCOT_SEQUENCE_DIR',link_files=>'MASCOT_LINK_FILES_FLAG'}
	);
	return %mascotServers;
}

####> FTP config <####
sub getFTPconfig {
	return (
		mode =>'passive' # or active
		#firewall => 'xxx'
	);
}

####> VALIDATION PARAMETERS <####
sub getMaxRank {	# Maximum number of query interpretations to be considered during MS analysis import
	return 10;		# Do not change this value if an imported analysis remains to be validated!!!
}
sub getMinScore {	# Default minimum score of query interpretations to be considered during MS analysis import
	my $fileFormat=$_[0];
	if (!$fileFormat || $fileFormat=~/MASCOT/) {return 20;} # Mascot search
	elsif ($fileFormat=~/PHENYX/) {return 5;} # Phenyx search
	elsif ($fileFormat=~/SEQUEST/) {return 1;} # Sequest search
	elsif ($fileFormat=~/PARAGON/) {return 90;} # Paragon search
    elsif ($fileFormat=~/TDM/) {return 1;} # X! Tandem search
	else {return 0;}
}
sub getMaxRefSpectra {	# maximum number of reference spectra allowed per peptide
	return 3;
}
sub getProteinPageSize {	# Average number of proteins displayed at once during validation
	return 150;
}
sub getQueryPageSize {	# Number of queries displayed at once during validation
	return 200;
}


####> DISPLAY CONFIGURATIONS <####
sub getRowColors { # HTML table row colors
	my ($lighter,$darker)=('#BDD7FF','#80B3FF');
	return ($lighter,$darker);
}
sub getPopupColors { # Popup layer colors
	my ($backgroundColor,$borderColor)=('#E5E5E5','#999999');
	return ($backgroundColor,$borderColor);
}
sub getItemIcones { # Navigation tree logos
	my %iconImg=(
		'project'=>'project.gif',
		'experiment'=>'experiment.gif',
		'experiment:locked' =>'experiment_locked.png',
		'sample'=>'sample.gif',
		'analysis:no_scan'=>'analysis_gray.gif',
		'analysis:no_val'=>'analysis_red.gif',
		'analysis:part_val'=>'analysis_yellow.gif',
		'analysis:val'=>'analysis_green.gif',
		'gel2d'=>'gel2D.gif',
		'spot'=>'spot.gif',
		'category'=>'category.gif',
		'mascot'=>'mascot.gif',
		'folder'=>'closed_folder.gif',
		'file'=>'file.gif',
		'badFile'=>'badFile.gif',
		'deniedRequest'=>'bad.gif',
		'design'=>'design.png',
		'expcondition'=>'category.gif',
		'quantification'=>'quantification.png', #'omega_blue.png',
		'quantification:peptide_no_label'=>'quanti_single.png',
		'quantification:peptide_no_label_swath'=>'swath.png',
		'quantification:peptide_label'=>'quantification.png',
		'quantification:peptide_label_swath'=>'swath.png',
        'quantification:protein_-2'=>'protein_none.gif',
        'quantification:protein_-1'=>'protein_error.gif',
        'quantification:protein_0'=>'protein_ongoing.gif',
        'quantification:protein_1'=>'protein_ended.gif',
        'quantification:protein_2'=>'protein_ended_massist.gif',
        'quantification:protein_3'=>'protein_ended_bioinfo.gif',
		'quantification:protein_label'=>'protein_labeled.gif',
		'no_label'=>'protein_unlabeled.gif', # not used
		'label'=>'protein_labeled.gif', # not used
		'go_analysis'=>'category.gif',
		'biosample' => 'cell.gif', # obsolete
		'biosample:not_used' => 'cell.gif',
		'biosample:used' => 'cell_color.gif',
		'exploranalysis:cluster' => 'clustering.png',
        'exploranalysis:clusterpep' => 'clustering.png',
		'exploranalysis:pca' => 'pca.png',
        'exploranalysis:pcapep' => 'pca.png',
		'pathwayanalysis' => 'pathway.gif',
		'motifanalysis:motif' => 'motifEnrichment.gif',
		'motifanalysis:hm_motif' => 'clustering.png'
		);
	return %iconImg;
}
sub getSpectrumColors { # Fragmentation spectrum color (added by FP)
	my ($higlightColor,$whiteColor,$blackColor,$colorA,$colorB)=('#99FFFF','#FFFFFF','#000000','#FF0000','#0000FF');
	return ($higlightColor,$whiteColor,$blackColor,$colorA,$colorB);
}

###>Maximum life-time of a partially validated Analysis before auto-end (in days)
sub getMaxPartialValidationDuration {
	return 183; # in days
}

################################################################################
##################> OTHER CONFIGURATIONS WITH DEFAULT VALUES <################## DO NOT CHANGE UNLESS YOU KNOW WHAT YOU ARE DOING !!!
################################################################################

sub getMsType { # Full names of MS types
	my %msType=('PMF'=>'Peptide Mass Fingerprint','MIS'=>'MS/MS Ions Search','SQ'=>'Sequence Query','sPMF'=>'PMF','sMIS'=>'MS/MS','sSQ'=>'SQ',
				'SWATH'=>'DIA','sSWATH'=>'SWATH-MS','TDA'=>'TDA','sTDA'=>'TDA','DIA'=>'DIA','sDIA'=>'DIA'
				);
	return (%msType);
}

sub getMassAAave { # Average mass, needed for protein mass calculation
	my %massAAave=(
		A=>71.0788,
		B=>114.59625,
		C=>103.1388,
		D=>115.0886,
		E=>129.1155,
		F=>147.1766,
		G=>57.0520,
		H=>137.1412,
		I=>113.1595,
		J=>0.0,
		K=>128.1742,
		L=>113.1595,
		M=>131.1925,
		N=>114.1039,
		O=>0.0,
		P=>97.1167,
		Q=>128.1308,
		R=>156.1876,
		S=>87.0782,
		T=>101.1051,
		U=>0.0,
		V=>99.1326,
		W=>186.2133,
		X=>111.0,
		Y=>163.1760,
		Z=>128.62315
	);
	return (%massAAave);
}

sub getMassATave { # Average mass, needed for protein mass calculation
	my %massATave=(
		H=>1.00794,
		C=>12.011,
		N=>14.0067,
		O=>15.9994
	);
	return (%massATave);
}

sub getTPPModificationCode {        # Modification code for Spectrast modifications conversion (editSwathLibrary.cgi)
    my %TPPModifCode=(
        'GlyGly'=>242
    );
    return %TPPModifCode;
}

#################################################################################
#################      Addition by FP            ################################
#################################################################################

sub getMassAAmono { # Monoisotopic Mass, needed for peptide Mass caculation #added by FP
	my %massValueAA = (
		A=>71.037110,
		B=>114.534930,
		C=>103.0091,
		D=>115.026940,
		E=>129.042590,
		F=>147.068410,
		G=>57.021460,
		H=>137.058910,
		I=>113.084060,
		J=>0.000000,
		K=>128.094960,
		L=>113.084060,
		M=>131.040490,
		N=>114.042930,
		O=>0.000000,
		P=>97.052760,
		Q=>128.058580,
		R=>156.101110,
		S=>87.032030,
		T=>101.047680,
		U=>150.953640,
		V=>99.068410,
		W=>186.079310,
		X=>111.000000,
		Y=>163.063330,
		Z=>128.550590,
		C_term=>17.002735,
		N_term=>1.007825
	);
	return (%massValueAA) ;
}

my %massAtMono = (
		H=>1.007825 ,
		C=>12.000000 ,
		N=>14.003070 ,
		O=>15.994910,
		CHO=>29.002735,
		CO=>27.99491,
		NH3=>17.026545,
		H2O => 18.01056
	);
sub getMassATMono { #added by FP
	return (%massAtMono);
}
sub getFragmentTolerance { #added by FP
	my $instrument = $_[0];
	my $tolerance ;
	if (($instrument=~/TRAP/) || ($instrument=~/QUAD/ && $instrument !~/TOF/)) {
		$tolerance=0.8;
	}
	else {
		$tolerance=0.2;
	}
	return $tolerance;
}
sub getNRLevel { # added by PP
	# my $instrument=$_[0];
	# add intrument management
	my $nrLevel=0.02;
	return $nrLevel;
}
sub getDeltaFragment {
	# my $instrument=$_[0];
	# add intrument management
	my $deltaFrag='Da';
	return $deltaFrag;
}
sub getDeltaParent {
	# my $instrument=$_[0];
	# add intrument management
	my $deltaParent='Da';
	return $deltaParent;
}
sub getFragmentationRules {
	my $instrument = $_[0];
	my %fragmentationRules;
	if ($instrument=~/ESI/) {
		%fragmentationRules = (
			'a'=>2,
			'a*'=>0,
			'a°'=>0,
			'b'=>2,
			'b*'=>1,
			'b°'=>1,
			'y'=>2,
			'y*'=>1,
			'y°'=>1,
			'ya'=>1,
			'yb'=>1,
			'+'=>1,
			'2+'=>1
		);
	}
	elsif ($instrument=~/MALDI/) {
		%fragmentationRules = (
			'a'=>2,
			'a*'=>1,
			'a°'=>1,
			'b'=>2,
			'b*'=>1,
			'b°'=>1,
			'y'=>2,
			'y*'=>1,
			'y°'=>1,
			'ya'=>2,
			'yb'=>2,
			'+'=>1,
			'2+'=>0
		);
	}
	else { #default
		%fragmentationRules = (
			'a'=>2,
			'a*'=>0,
			'a°'=>0,
			'b'=>2,
			'b*'=>1,
			'b°'=>1,
			'y'=>2,
			'y*'=>1,
			'y°'=>1,
			'ya'=>1,
			'yb'=>1,
			'+'=>1,
			'2+'=>1
		);
	}
	return %fragmentationRules;
}
sub getFragmentDef {
	my %fragmentDef=(
		'a' => -($massAtMono{"C"} + $massAtMono{"H"} + $massAtMono{"O"}),
		'b' => -$massAtMono{"H"},
		'c' => $massAtMono{"N"} + 2 * $massAtMono{"H"},
		'x' => $massAtMono{"C"} + $massAtMono{"O"} - $massAtMono{"H"},
		'y' => $massAtMono{"H"},
		'z' => -( $massAtMono{"N"} + 2 * $massAtMono{"H"} ),
		'z+1' => -( $massAtMono{"N"} + 2 * $massAtMono{"H"} ) + $massAtMono{"H"},
		'z+2' => -( $massAtMono{"N"} + 2 * $massAtMono{"H"} ) + 2 * $massAtMono{"H"},
		'°' => -( 2 * $massAtMono{"H"} + $massAtMono{"O"} ),
		'*' => -( 3 * $massAtMono{"H"} + $massAtMono{"N"} ),
		' ' => 0
	);
	$fragmentDef{"ya"} = $fragmentDef{"y"} + $fragmentDef{"a"};
	$fragmentDef{"yb"} = $fragmentDef{"y"} +$fragmentDef{"b"};
	return %fragmentDef;
}

sub getFragmentClassif {
	my %fragmentClassif;
	@{$fragmentClassif{"C_term"}} = ('x','y','z','z+1','z+2');
	@{$fragmentClassif{"N_term"}} = ('a','b','c');
	@{$fragmentClassif{"intern"}} = ('ya','yb');
	@{$fragmentClassif{"neutral_loss"}} = ('°','*') ;

	return %fragmentClassif;
}
1;

####>Revision history<####
# 2.11.11D Adapted for distribution (PP 04/06/21)
# 2.11.11 [UPDATE] sql_mode declaration in DB connection (PP 27/05/21)
# 2.11.10 [UPDATE] Minor change to allow spaces in user-defined cluster errors to be ignored (PP 22/04/21)
# 2.11.9 [UPDATE] Added "paper" instance for dev (PP 11/04/21)
# 2.11.8 [FEATURE] Add path of gsea folders (VL 29/10/20)
# 2.11.7 [CHANGE] Prod uses mypromsdb instead of bioinfo-db-prod (PP 07/03/21) 
# 2.11.6 [CHANGE] Uses Singularity image 1.3.5 on cluster (VS 29/01/21)
# 2.11.5 [CHANGE] Database lock time out set on 600s at connection (PP 29/01/21)
# 2.11.4 [UPDATE] $cluster{runJob} also returns cluster output file (PP 28/01/21)
# 2.11.3 [UPDATE] Modified cluster connection warning message "(Warning: )Permanently (added)" to be ignored has job error (PP 02/01/21)
# 2.11.2 [CHANGE] Uses default MySQL 8 sql_mode (no downgrade to MySQL 5) on connection to smorpym (PP 04/01/21)
# 2.11.1 [UPDATE] Added new message "Warning: Permanently added" to be ignored during cluster job error check (PP 11/12/20)
# 2.11.0 [CHANGE] Uses smorpym MySQL server for dev (PP 02/12/20)
# 2.10.8 [CHANGE] Added mascot IP url address to be able to use it on cluster (VS 27/10/20)
# 2.10.7 [CHANGES] Use Singularity image 1.3.4 on cluster (21/10/20)
# 2.10.6 [ENHANCEMENT] Declared paths for python, python2 and python3 & new path for singularity image 1.3.2 (PP 28/07/20)
# 2.10.5 [BUGFIX] Fix find command of job status folder (VS 07/04/20)
# 2.10.4 [CHANGE] Added status-based quantification icon (VS 21/02/20)
# 2.10.3 [CHANGE] Cluster max number of jobs sets to 50 instead of 20 (PP 07/02/20)
# 2.10.2 [ENHANCEMENT] Added alternative cluster status & singularity image for docker instance (PP 15/01/20)
# 2.10.1 [BUGFIX] In docker-based server detection (PP 14/01/20)
# 2.10.0 [FEATURE] Added server detection for docker container (PP 11/01/20)
# 2.9.22b [BUGFIX] Check torque6 log files path existence before grep (VS 08/01/20)
# 2.9.22 [CHANGES] Minor change in $ENV{REMOTE_USER} filtering in dev mode (PP 25/11/19)
# 2.9.21 [CHANGES] Change TPP path + Improve cluster killJob function + use Singularity image 1.3.2 on cluster (VS 18/11/19)
# 2.9.20 [CHANGES] Uses Singularity image 1.3.1 on cluster (PP 12/11/19)
# 2.9.19 [FEATURE] Added 'killJob' and 'jobInfos' to &getClusterInfo and (VS 16/10/19)
# 2.9.18 [FEATURE] Added cluster singularity image setter (PP 26/08/19)
# 2.9.17 [FEATURE] Added experiment lock system img (VS 08/08/19)
# 2.9.16 Uses new &promsConfig::getMascotServers function (PP 25/06/19)
# 2.9.15 &getFTPmode function now returns a hash, not a reference (PP 24/06/19)
# 2.9.14 Uses new mysql production server smorpym.curie.fr (PP 12/06/19)
# 2.9.13 &getFTPmode function (PP 11/06/19) 
# 2.9.12 Rename project_data folder into metadata (VS 11/06/19)
# 2.9.11 Add project data path for file metadata storage (VS 05/06/19)
# 2.9.10 Different icons for used/not_used BioSamples (PP 11/04/19) 
# 2.9.9 Added conditional BASH command 'export LC_ALL=\"C\"' in $cluster{runJob} function for image 1.1.17 (PP 07/01/19)
# 2.9.8 [Fix] bug in max. memory management for cluster & declare $CGI::LIST_CONTEXT_WARN=0 to prevent warning from using CGI param instead of multi_param in list context (PP 19/12/18)
# 2.9.7 Back to Singularity image myproms_1.1.17.img & filtering of Perl warning in &getClusterInfo{checkError} (PP 02/11/18)
# 2.9.6 Uses Singularity image myproms_1.1.19-1.img & updates in &getClusterInfo & better detection of git_demo server (PP 03/10/18)
# 2.9.5 New wrapping function $cluster{runJob}->() for job launch (PP 11/09/18)
# 2.9.4 Improved server config detection for scripts running on cluster & filters mounting warnings from singularity in $cluster->checkError (PP 30/08/18)
# 2.9.3 Minor modification in clusterInfo{path} (GA 25/06/18)
# 2.9.2 Handles webserver development sub-branch eg. ppoullet-refonte (PP 18/06/18)
# 2.9.1 Added &getMaxPartialValidationDuration for control of the auto-end validation process (PP 07/06/18)
# 2.9.0 Changes in &detectServer and &getClusterInfo to improve cluster support (PP 11/04/18)
# 2.8.9 Uses image myproms_1.1.17.img on CentOS cluster & added 'maxCPU' key in %clusterInfo (PP 24/01/18)
# 2.8.8 Another minor change in &getClusterInfo to improve file vs. command as parameter (PP 09/01/18)
# 2.8.7 Minor change in &getClusterInfo to improve file vs. command as parameter (PP 20/12/17)
# 2.8.6 Uses image myproms_1.1.16.img on CentOS cluster & added pyprophet path (PP 05/12/17)
# 2.8.5 Added &getTPPModificationCode (MLP 17/11/17)
# 2.8.4 Change database passwords for dev/valid/prod (PP 08/11/17)
# 2.8.3 Call of CentOS cluster is now instance-specific (PP 07/11/17)
# 2.8.2 Added img for peptide exploratory analysis (SL 06/11/17)
# 2.8.1 Added $promsPath{shared} in all instances (PP 31/10/17)
# 2.8.0 Added a dedicated function to retrieve computation cluster info: $promsPath{qsub} is depracated! & DIA type in &getMsType (PP 18/10/17)
# 2.7.4 Icons for Motif analysis (PP 16/08/17)
# 2.7.3 Add TDA type in &getMsType (PRM, MRM, SRM) (MLP 23/06/17)
# 2.7.2 Modification on &getMinScore : add TANDEM && add python_scripts path (MLP 13/01/17)
# 2.7.1 Bug fix in %fragmentDef for a and b masses & MySQL server port is now a variable (PP 07/11/16)
# 2.7.0 Added path for java & SWATH in &getMsType (PP 20/09/16)
# 2.6.9 Icon for SWATH peptide quantification (PP 02/08/16)
# 2.6.8 Minor code cleaning (PP 11/07/2016)
# 2.6.7 Minor modification on swath_lib path (MLP 07/07/2016)
# 2.6.6 Remove unused wget path in prod (PP 22/04/16)
# 2.6.5 Add paths for swath_lib, tpp, python & msproteomicstools (MLP 13/04/16)
# 2.6.4 Change myProMS form colors (PP 14/03/16)
# 2.6.3 New mascot server (GA 08/01/16)
# 2.6.1 add new Path tmp_html (SL 18/12/14)
# 2.6.0 add new image and new directories for pathway analysis (SL 21/10/14)
# 2.5.9 add new images cell, cluster and pca (SL 03/07/14)
# 2.5.8 add new path exploratory analyses (SL 13/06/14)
# 2.5.7 Add z+1 and z+2 for ETD fragmentation in &getFragmentClassif and &getFragmentDef (GA 15/04/14)
# 2.5.6 Check new version of MCQ 2.0.1 (GA 11/10/13)
# 2.5.5 %promsPath rearangement (PP 20/09/13)
# 2.5.4 Added qsub path in %promsPath (PP 06/09/13)
# 2.5.3 New icons for label(-free) quantification (PP 03/09/13)
# 2.5.2 Added no_label & label item in icon list (PP 04/06/13)
# 2.5.1 Move getVariableModifications to promsMod (GA 15/05/13)
# 2.5.0 Uses mysql_init_command to execute 'SET NAMES utf8' (PP 17/04/13)
# 2.4.9 Added 'SET NAMES utf8' at DB connexion (PP 17/04/13)
# 2.4.8 Bug fix for new valid/prod pipeline (PP 15/04/13)
# 2.4.7 New pipelines + %promsPath reorganization (PP 08/04/13)
# 2.4.6 &getMinScore add Paragon (GA 20/03/13)
# 2.4.5 Adding right fragment def for x ions (s/+CC/+CO/) (FY 19/03/13)
# 2.4.4 &getIdent2AliasModString no longer required (PP 12/03/13)
# 2.4.3 Updated Mac configuration (PP 21/03/13)
# 2.4.2 Added NCBI_ALL to &getIdentifierTypes and &getIdent2AliasModString (PP 28/01/13)
# 2.4.1 Added path to logs directory (PP 23/11/12)
# 2.4.0 Modifs for new pipeline (PP 03/10/12)
# 2.3.8 Added qvality & masschroq pathes in %promsPath (PP 03/08/12)
# 2.3.7 Added (Di/Tri)methylation in list of project VarMods (PP 11/05/12)
# 2.3.6 New icons for design & quantification (PP 19/04/12)
# 2.3.5 Change bin directory (FY 04/03/12)
# 2.3.4 Merge 2.2.3GA (Label-free) & 2.3.3PP (PP 02/04/12)
# 2.3.3 Added Flybase polypeptide id in &getIdentifierTypes (PP 26/03/12)
# 2.3.2 Added quanti & GO folders to production promsPath (PP 23/03/12)
# 2.3.1 Added GO data folders to promsPath (FY 19/01/12)
# 2.3.0 New prod server config params updated (PP 05/01/12)
# 2.2.9m Added 'quanti_html' key to %promsPath (PP 14/10/11)
# 2.2.8m databanks moved inside /bioinfo/projects_prod/myproms for both prod servers (PP 02/10/11)
# 2.2.7m DEV: tmp dir moved outside /bioinfo/projects_dev/tmp (PP 30/09/11)
# 2.2.7 Added new production server myproms.curie.fr (PP 30/09/11)
# 2.2.6 Handling SILAC 13C labeling (PP 13/07/11)
# 2.2.5 Updated for myproms-dev.curie.fr (PP 01/07/11)
# 2.2.4 Compatibility with mixed (SQ) search (PP 16/04/2011)
