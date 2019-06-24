#!/usr/local/bin/perl -w

################################################################################
# testMyProMS.cgi         2.5.2                                                #
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
use POSIX qw(strftime);
use strict;

###>Executables (binaries)<###
my @execPath=qw(
	java
	masschroq
	msproteomicstools
	openms
	pyprophet
	python
	qvality
	R
	tpp
);
#	qsub

###>List of R packages<###
my @Rpackages=qw(
	broom
	data.table
	dplyr
	FactoMineR
	ggdendro
	ggplot2
	ggseqlogo
	gplots
	gridExtra
	igraph
	limma
	missMDA
	MSstats
	multtest
	outliers
	plgem
	plyr
	preprocessCore
	purrr
	qvalue
	readr
	reshape2
	stringr
	tidyr
	tidyverse
);
#	affy
#	Biobase
#	BiocGenerics
#	ellipse
#	GGally
#	grid
#	imputeLCMD
#	parallel

@Rpackages=sort {lc($a) cmp lc($b)} @Rpackages;

#####>header<#####
print header(-'content-encoding'=>'no',-charset=>'utf-8');
warningsToBrowser(1);

my $action=param('ACT') || 'main';
my $call=param('CALL') || 'user';

#####>AJAX calls<#####
if (param('AJAX')) {
	require promsConfig;
	require File::Path; # qw(remove_tree); # qw(rmtree);
	my %promsPath=&promsConfig::getServerInfo;
	#if ($action eq 'R') {
	#	my $package=param('PACKAGE');
	#	my @RAnswer=`export LANG=en_US.UTF-8; $promsPath{R}/R -e "require($package); sessionInfo();" 2>&1`;
	#	#print "@RAnswer\n\n\n\n";
	#	my $responseOK=0;
	#	foreach (@RAnswer) {
	#		if (/^\s*\[\d+\].*\s$package\_(\S+)/) { # as 'other attached packages'
	#			print "&nbsp;Version $1 found\n";
	#			$responseOK=1;
	#			last;
	#		}
	#		elsif (/^\s*\[\d+\].*\s$package($|\s)/) { # as 'attached base packages'
	#			print "&nbsp;Found as base package\n";
	#			$responseOK=1;
	#			last;
	#		}
	#		elsif (/no package called \‘$package\’/) {
	#			print "&nbsp;<FONT color=red>Error: Package not found!</FONT>\n";
	#			$responseOK=1;
	#			last;
	#		}
	#	}
	#	unless ($responseOK) {
	#		print "&nbsp;<FONT color=blue>Unexpected response from server:</FONT><CODE style=\"display:block;white-space:pre-wrap;;background:#E0E0E0\">",join("\n",@RAnswer),"</CODE>\n";
	#	}
	#}
	if ($action eq 'job') {
		my $clusterIdx=param('CLUS');
		my $timeStamp=param('DIR');

		###>Run dir<###
		my $startJob=0;
		unless ($timeStamp) {
			mkdir "$promsPath{tmp}/test" unless -d "$promsPath{tmp}/test";
			$timeStamp=strftime("%Y%m%d%H%M%S",localtime);
			mkdir "$promsPath{tmp}/test/$timeStamp";
			$startJob=1;
		}
		my $runDir="$promsPath{tmp}/test/$timeStamp";

		###>local job<###
		if ($clusterIdx==-1) {
			my $Rscript="$runDir/Rscript.R";

			##>Launching job
			if ($startJob) {
				open (R,">$Rscript");
				foreach my $package (@Rpackages) {
					print R "require($package)\n";
				}
				print R qq
|#### Session info ####
sink(file.path("$runDir","RsessionInfo.txt"))
print(sessionInfo())
sink()
|;
				close R;

				my $childPid = fork;
				unless ($childPid) { # child here
					#>Disconnecting from server
					open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
					open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
					#>Command
					system "export LANG=en_US.UTF-8; cd $runDir; $promsPath{R}/R CMD BATCH --no-save --no-restore $Rscript; echo _END_ > end.txt";
					exit;
				}
				print "START#==========#$timeStamp#==========#";
			}
			else {
				if (-e "$runDir/end.txt") {
					my $sessionInfo="$runDir/RsessionInfo.txt";
					my $RoutScript=$Rscript.'out';
					if (-e $sessionInfo) {
						print "END#==========#";
						open(R_INFO,$sessionInfo);
						my @Rinfo=<R_INFO>;
						close R_INFO;
						my %sessionPackages;
						my ($version)=($Rinfo[0] && $Rinfo[0]=~/R version (\S+)/);
						foreach my $line (@Rinfo) {
							next unless $line =~ /\[.+_/;
							$line=~s/^.+\]\s*//;
							foreach my $packData (split/\s+/,$line) {
								my ($pack,$vers)=($packData=~/(.+)_([^_]+)$/);
								$sessionPackages{$pack}=$vers;
							}
						}
						&printRpackageVersion(\%sessionPackages);
						print "<BR>\n";
					}
					elsif (-e $RoutScript) {
						my $errorStrg=`tail -20 $RoutScript`;
						print "END#==========#<B>ERROR: Server responded:<BR>\n<PRE>$errorStrg</PRE>";
					}
					File::Path->remove_tree($runDir);
				}
				else {print "WAIT#==========#";}
			}
		}

		###>Cluster job<###
		else {
			my %clusterInfo=&promsConfig::getClusterInfo;
			my $refClusInfo={};
			if ($clusterInfo{'list'}) {%{$refClusInfo}=&promsConfig::getClusterInfo($clusterInfo{'list'}[$clusterIdx]);}
			else {$refClusInfo=\%clusterInfo;}

			if ($startJob) {

				#>Command file
				my $commandFile="$runDir/command.sh";
				open (COMM,">$commandFile");
				print COMM "#!/bin/bash\n";
				foreach my $pathName (sort{lc($a) cmp lc($b)} @execPath) {
					next unless defined $refClusInfo->{'path'}{$pathName};
					print COMM "\necho \">>>>$pathName:\"\n";
					if ($pathName eq 'java') {
						print COMM "$refClusInfo->{path}{$pathName}/java -version 2>&1\n";
					}
					elsif ($pathName eq 'masschroq') {
						print COMM "$refClusInfo->{path}{$pathName}/masschroq -v 2>&1\n";
					}
					elsif ($pathName eq 'msproteomicstools') {
						print COMM "pip show msproteomicstools 2>&1 | grep -P \"Version:\"\n";
						print COMM "grep msproteomicstools $refClusInfo->{path}{$pathName}/spectrast2tsv.py 2>&1\n";
					}
					elsif ($pathName eq 'openms') {
						print COMM "$refClusInfo->{path}{$pathName}/OpenMSInfo 2>&1\n";
					}

					elsif ($pathName eq 'pyprophet') {
						#print COMM "pip show pyprophet 2>&1 | grep -P \"Version:\"\n";
						print COMM "$refClusInfo->{path}{$pathName}/pyprophet --version 2>&1\n";
					}

					elsif ($pathName eq 'python') {
						print COMM "$refClusInfo->{path}{$pathName}/python -V 2>&1\n";
					}
					elsif ($pathName eq 'qvality') {
						print COMM "$refClusInfo->{path}{$pathName}/qvality -h 2>&1\n";
					}
					elsif ($pathName eq 'R') {
						my $Rscript="$runDir/Rscript.R";
						open (R,">$Rscript");
						foreach my $package (@Rpackages) {
							print R "require($package)\n";
						}
						print R qq
|#### Session info ####
sink(file.path("$runDir","RsessionInfo.txt"))
print(sessionInfo())
sink()
|;
						close R;
						print COMM "export LANG=en_US.UTF-8; cd $runDir; $refClusInfo->{path}{$pathName}/R CMD BATCH --no-save --no-restore $Rscript";
					}
					elsif ($pathName eq 'tpp') {
						print COMM "$refClusInfo->{path}{$pathName}/xinteract | grep TPP\n";
					}
				}
				close COMM;
				my $modBash=0775;
				chmod $modBash,$commandFile;

###				my $clusterCommand=$refClusInfo->{'buildCommand'}->($runDir,$commandFile);
###
###				#>qsub file
###				my $bashFile="$runDir/test_myproms.sh";
###				open (BASH,">$bashFile");
###				print BASH qq
###|#!/bin/bash
#####resources
####PBS -l mem=1Gb
####PBS -l nodes=1:ppn=1
####PBS -l walltime=1:00:00
####PBS -q batch
#####Information
####PBS -N myProMS_test_$timeStamp
#####PBS -M patrick.poullet\@curie.fr
####PBS -m abe
####PBS -o $runDir/PBS.txt
####PBS -e $runDir/PBSerror.txt
###
##### Command
###$clusterCommand
###
###echo _END_$timeStamp
###|;
###				close BASH;
###				#my $modBash=0775;
###				chmod $modBash, $bashFile;
###				$refClusInfo->{'sendToCluster'}->($bashFile);

my %jobParams=(
	maxMem=>'1Gb',
	numCPUs=>1,
	maxHours=>1,
	jobName=>"myProMS_test_$timeStamp",
	outFile=>'PBS.txt',
	errorFile=>'PBSerror.txt',
	jobEndFlag=>"_END_$timeStamp",
	noWatch=>1
);
$refClusInfo->{'runJob'}->($runDir,$commandFile,\%jobParams);

				print "START#==========#$timeStamp";
			}
			else { # check job
				#my $pbsError=(-e "$runDir/PBSerror.txt")? `less $runDir/PBSerror.txt` : '';
				my $pbsError=$refClusInfo->{'checkError'}->("$runDir/PBSerror.txt");
				if ($pbsError) {
					print "END#==========#<B>ERROR: Cluster responded:<BR>\n<PRE>$pbsError</PRE>";
					File::Path->remove_tree($runDir);
				}
				elsif (-e "$runDir/PBS.txt" && `tail -3 $runDir/PBS.txt | grep _END_$timeStamp`) {
					print "END#==========#";
					open (PBS,"$runDir/PBS.txt");
					my @results=<PBS>;
					close PBS;
					my $i=-1;
					while ($i < $#results) {
						$i++;
						next unless $results[$i]=~/^>>>>(.+):/;
						my $pathName=$1;
						my @response;
						foreach my $j ($i+1..$#results) {
							last if $results[$j]=~/^>>>>/;
							$i++;
							last if $results[$j]=~/^_END_/;
							push @response,$results[$j];
						}
						my $version;
						my %sessionPackages; # for R only
						if ($pathName eq 'java') {
							($version)=($response[0] && $response[0]=~/ version "([^"]+)"/);
						}
						elsif ($pathName eq 'R') {
							#($version)=($response[0] && $response[0]=~/R batch front end: (\S+)/);
							open(R_INFO,"$runDir/RsessionInfo.txt");
							my @Rinfo=<R_INFO>;
							close R_INFO;
							($version)=($Rinfo[0] && $Rinfo[0]=~/R version (\S+)/);
							foreach my $line (@Rinfo) {
								next unless $line =~ /\[.+_/;
								$line=~s/^.+\]\s*//;
								foreach my $packData (split/\s+/,$line) {
									my ($pack,$vers)=($packData=~/(.+)_([^_]+)$/);
									$sessionPackages{$pack}=$vers;
								}
							}
						}
						elsif ($pathName eq 'qvality') {
							($version)=($response[0] && $response[0]=~/qvality version ([\w\.-]+)/);
						}
						elsif ($pathName eq 'masschroq') {
							#($version)=($response[0] && $response[0]=~/MassChroQ version (\S+),/i);
							($version)=($response[0] && $response[0]=~/MassChroQ [version ]*([\d\.]+)/i);
						}
						elsif ($pathName eq 'tpp') {
							($version)=($response[0] && $response[0]=~/\(TPP ([^,]+)/);
						}
						elsif ($pathName eq 'python') {
							($version)=($response[0] && $response[0]=~/Python (.+)/);
						}
						elsif ($pathName eq 'msproteomicstools') {
							($version)=($response[0] && $response[0]=~/Version:\s*(\S+)/); # pip
							unless ($version) { # spectrast2tsv.py
								foreach my $line (@response) {
									($version)=($line && $line=~/msproteomicstools==([\d\.]+)/);
									last if $version;
								}
							}
						}
						elsif ($pathName eq 'openms') {
							($version)=($response[2] && $response[2]=~/Version\s*:\s*(\S+)/);
						}
						elsif ($pathName eq 'pyprophet') {
							#($version)=($response[0] && $response[0]=~/Version:\s*(\S+)/); # pip
							#unless ($version) {
							#	foreach my $line (@response) {
							#		($version)=($line && $line=~/^([\d\.]+)/);
							#		last if $version;
							#	}
							#}
							($version)=($response[0] && $response[0]=~/([\d\.]+)/);
						}
						&printVersionResponse($pathName,$refClusInfo->{path}{$pathName},$version,\@response,$clusterIdx,\%sessionPackages);
					}
					File::Path->remove_tree($runDir);
				}
				else {print "WAIT";}
			}
		}
	}
	exit;
}

###################> Main <###################
print qq
|<HTML>
<HEAD>
<TITLE>Testing myProMS</TITLE>
|;
my %promsPath;
if ($call eq 'server') {
	require promsConfig;
	%promsPath=&promsConfig::getServerInfo;
	print qq|<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">|;
}
if ($action=~/internet/) {
	&checkInternetConnection;
	exit;
}
elsif ($action eq 'lwp') {
	print "</HEAD><BODY>TEST IS OK</BODY></HTML>\n";
	exit;
}

print qq
|<SCRIPT LANGUAGE="JavaScript">
var startDate = new Date();
function updateVariablesDisplay() {
	if (document.getElementById('envDiv').style.display=='none') {
		document.getElementById('envButton').value='Hide server variables';
		document.getElementById('envDiv').style.display='block';
	}
	else {
		document.getElementById('envButton').value='Show server variables';
		document.getElementById('envDiv').style.display='none';
	}
}
function updateScriptDisplay() {
	if (document.getElementById('scriptDiv').style.display=='none') {
		document.getElementById('scriptButton').value='Hide versions';
		document.getElementById('scriptDiv').style.display='block';
	}
	else {
		document.getElementById('scriptButton').value='Show versions';
		document.getElementById('scriptDiv').style.display='none';
	}
}
function showHideRPackages(clusterIdx) {
	var packageDiv=document.getElementById('Rpackages_'+clusterIdx+'DIV');
	packageDiv.style.display=(packageDiv.style.display=='none')? '' : 'none';
}

//AJAX
/*** Job: local or cluster(s) ***/
function checkBinaries(clusterIdx,name) {
	var jobDiv,resSpanId,infoStrg;
	if (clusterIdx==-1) { // local binaries (R only)
		resSpanId='waitSPAN';
		document.getElementById('RpackagesBUT').style.display='none';
		jobDiv=document.getElementById('RpackagesDIV');
		infoStrg='Starting job';
	}
	else { // cluster
		resSpanId='wait'+clusterIdx+'SPAN';
		document.getElementById('cluster'+clusterIdx+'BUT').style.display='none';
		jobDiv=document.getElementById('cluster'+clusterIdx+'DIV');
		infoStrg='Launching job on cluster '+name;
	}
	jobDiv.innerHTML='<IMG src="$promsPath{images}/scrollbarGreen.gif"><BR>&nbsp;&nbsp;&nbsp;&nbsp;<B><SPAN id="'+resSpanId+'">'+infoStrg+'...</SPAN></B>';
	ajaxWatchJob(clusterIdx,null);
}
function ajaxWatchJob(clusterIdx,timeStamp) {
	//If XHR object already exists, the request is canceled & the object is deleted
	if (XHR && XHR.readyState != 0) {
		XHR.abort();
		delete XHR;
	}
	//Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}

	var resSpanId,buttonId,divId;
	if (clusterIdx==-1) {
		resSpanId='waitSPAN';
		buttonId='RpackagesBUT';
		divId='RpackagesDIV';
	}
	else {
		resSpanId='wait'+clusterIdx+'SPAN';
		buttonId='cluster'+clusterIdx+'BUT';
		divId='cluster'+clusterIdx+'DIV';
	}
	var tsParamStrg=(timeStamp)? "&DIR="+timeStamp : '';
	XHR.open("GET","$promsPath{cgi}/testMyProMS.cgi?AJAX=1&ACT=job&CLUS="+clusterIdx+tsParamStrg,true); // asynchronous
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			var info=XHR.responseText.split('#==========#');
			var resSpan=document.getElementById(resSpanId);
			if (info[0].match('START')) {
				timeStamp=info[1];
				resSpan.innerHTML='Job launched. Waiting for results...';
				setTimeout(function() {ajaxWatchJob(clusterIdx,timeStamp)},5000);
			}
			else if (info[0].match('WAIT')) {
				resSpan.innerHTML+='.';
				setTimeout(function() {ajaxWatchJob(clusterIdx,timeStamp)},5000); // timeStamp is defined
			}
			else if (info[0].match('END')) { // finished
				document.getElementById(buttonId).style.display='';
				document.getElementById(divId).innerHTML=info[1];
			}
		}
	}
	XHR.send(null);
}

var XHR = null;
function getXMLHTTP() {
    var xhr = null;
    if (window.XMLHttpRequest) {// Firefox & others
        xhr = new XMLHttpRequest();
    }
    else if (window.ActiveXObject) { // Internet Explorer
        try {
          xhr = new ActiveXObject("Msxml2.XMLHTTP");
        } catch (e) {
            try {
                xhr = new ActiveXObject("Microsoft.XMLHTTP");
            } catch (e1) {
                xhr = null;
            }
        }
    }
    else { // XMLHttpRequest not supported by browser
        alert("Your browser does not support XMLHTTPRequest objects...");
    }
    return xhr;
}
</SCRIPT>
</HEAD>
|;
if ($call eq 'server') {
	my ($light,$dark)=&promsConfig::getRowColors;
	print qq
|<BODY background="$promsPath{images}/bgProMS.gif">
<DIV style="float:top">
<BR>
<TABLE><TR><TH bgcolor="$dark">
<FONT class="title2">&nbsp;Go to:</FONT><SELECT class="title action" onchange="if (this.value) {window.location='$promsPath{cgi}/'+this.value;}">
	<OPTION value="">-= Select =-</OPTION>
	<OPTION value="promsMain.cgi">Main Window</OPTION>
	<OPTION value="selectProject.cgi">Project Selection</OPTION>
</SELECT>
</TH></TR></TABLE>
</DIV>
|;
}
else {
	print "<BODY>\n<BR>\n";
}

###>Detecting myProMS version
my $serverRootUnix;
if ($ENV{'WEBSERVER_ROOT'}) {
	$serverRootUnix=$ENV{'WEBSERVER_ROOT'};
}
else {
	$serverRootUnix=$ENV{'SCRIPT_FILENAME'} || $0;
	$serverRootUnix=~s/$promsPath{cgi}.+//;
}
my $indexUnix=(-d "$serverRootUnix$promsPath{html}" && -e "$serverRootUnix$promsPath{html}/index.html")? "$serverRootUnix$promsPath{html}/index.html" : (-d "$serverRootUnix/html$promsPath{html}")? "$serverRootUnix/html$promsPath{html}/index.html": undef;
my $mypromsVersion='Unknown';
if ($indexUnix) {
	my $versionLine=`grep '<!--version tag-->' $indexUnix`;
	($mypromsVersion)=$versionLine=~/ tag-->(\S+)/;
}
my $pwd=`pwd`; chomp($pwd);
print qq
|<H3>Checking myProMS installation: <INPUT type="button" value="Re-run test" onclick="window.location.reload()" style="font-weight:bold;font-size:14px"/></H3>
<BR>
<B>>myProMS version is $mypromsVersion.</B><BR><BR>
<B>>myProMS Server IP address is $ENV{SERVER_ADDR}.</B>&nbsp;<INPUT type="button" id="envButton" value="Show server variables" onclick="updateVariablesDisplay()"/>
<DIV id="envDiv" style="display:none">
<FIELDSET><LEGEND><FONT style="font-weight:bold;font-size:18px">Server variables:</FONT></LEGEND>
<B>>Script variables:</B><BR>
-Perl \$0: $0<BR>
-system `pwd`: $pwd<BR>
<BR><B>>Known pathes (\@INC):</B><BR>
|;
foreach my $v (@INC) {print "-$v<BR>\n";}
print "<BR><B>>Environment variables (\%ENV):</B><BR>\n";
foreach my $v (sort{lc($a) cmp lc($b)} keys %ENV) {print "$v: $ENV{$v}<BR>\n";}
print "</DIV><BR><BR>\n";

####>Checking Perl modules<####
my @requiredModules=qw(
	Archive::Tar
	CGI::Carp
	DBD::mysql
	DBD::SQLite
	DBI
	Encode
	File::Basename
	File::Compare
	File::Copy
	File::Copy::Recursive
	File::Find
	FileHandle
	File::Listing
	File::Path
	File::stat
	GD
	GO::AnnotationProvider::AnnotationParser
	GO::Node
	GO::OntologyProvider::OboParser
	GO::TermFinder
	GraphViz
	HTML::Entities
	HTTP::Request
	IO::Handle
	IO::Uncompress::Gunzip
	IO::Uncompress::Unzip
	JSON
	List::Util
	LWP::Simple
	LWP::UserAgent
	MIME::Base64
	Net::FTP
	Spreadsheet::WriteExcel
	Storable
	String::Util
	Text::CSV
	URI::Escape
	XML::SAX
	XML::SAX::ParserFactory
	XML::Simple
);	# CGI POSIX

print "<B>>Checking ",scalar @requiredModules," required Perl modules...</B><BR>\n";
my $lwpOK=0;
foreach my $mod0 (@requiredModules) {
	next if $mod0=~/^#/;
	(my $mod="$mod0.pm")=~s|::|/|g; # Foo::Bar::Baz => Foo/Bar/Baz.pm
	eval {require $mod};
	if ($@) {print "&nbsp;&nbsp;-<FONT color=red>Error: Module <B>$mod0</B> FAILED!</FONT><BR>\nERROR: $@<BR>\n";}
	else {
		print "&nbsp;&nbsp;-Module <B>$mod0</B> ",$mod0->VERSION," found<BR>\n";
		$lwpOK=1 if $mod0 eq 'LWP::UserAgent';
	}
}
print "<B>Done.</B><BR><BR>\n";

print "<B>>Checking myProMS modules...</B><BR>\n";
my @mypromsModules=qw(
	promsConfig
	promsMod
	promsQuantif
	promsOboParser
	goAnalysis
);
my $okConfig=0;
foreach my $mod0 (@mypromsModules) {
	next if $mod0=~/^#/;
	my $mod="$mod0.pm";
	eval {require $mod};
	if ($@) {print "&nbsp;&nbsp;-<FONT color=red>Error: Module <B>$mod0</B> FAILED!</FONT><BR>\nERROR: $@<BR>\n";}
	else {
		print "&nbsp;&nbsp;-Module <B>$mod0</B> found<BR>\n";
		$okConfig=1 if $mod0 eq 'promsConfig';
	}
}
print "<B>Done.</B><BR><BR>\n";

unless ($okConfig) {
	print "<BR><BR><B><FONT color=red>Aborting test until promsConfig.pm properly edited</FONT></B>\n";
	exit;
}

####>Checking myProMS data paths<####
print "<B>>Checking myProMS data paths...</B><BR>\n";
%promsPath=&promsConfig::getServerInfo('no_user') unless $call eq 'server'; # promsConfig loaded during eval
my @dataPath=qw(data banks swath_lib valid peptide spectrum gel_unix goa obo go_unix quantification explorAna pathAna tmp logs);
foreach my $pathName (@dataPath) {
	next unless $promsPath{$pathName};
	if ($pathName eq 'data') {
		if (-e $promsPath{$pathName}) {
			print "&nbsp;&nbsp;-Path for master <B>data</B> directory is found ($promsPath{data})<BR>\n";
		}
		else {
			print "&nbsp;&nbsp;<FONT color=red>Error: Path for master <B>data</B> directory was not FOUND!!! ($promsPath{data}). Secondary data directories cannot be checked.</FONT><BR>\n";
			last;
		}
		next;
	}
	if (-e $promsPath{$pathName}) {
		print "&nbsp;&nbsp;-Path for '<B>$pathName</B>' is found ($promsPath{$pathName})<BR>\n";
		&testWriteDirectory($promsPath{$pathName});
	}
	else { # create dir
		if ($promsPath{$pathName}=~/\/go\// && !-e "$promsPath{data}/go") { # goa,obo,go_unix
			mkdir "$promsPath{data}/go";
		}
		mkdir $promsPath{$pathName};
		if (-e $promsPath{$pathName}) {
			print "&nbsp;&nbsp;-Path for '<B>$pathName</B>' was succesfully created ($promsPath{$pathName})<BR>\n";
			&testWriteDirectory($promsPath{$pathName});
		}
		else {print "&nbsp;&nbsp;<FONT color=red>Error: Could not create path for '<B>$pathName</B>' ($promsPath{$pathName})</FONT><BR>\n";}
	}
}

print "<B>Done.</B><BR><BR>\n";

####>Checking myProMS server paths<####
print "<B>>Checking myProMS server paths...</B><BR>\n";
my %promsUnixPath=($call eq 'server')? &promsConfig::getServerInfo() : &promsConfig::getServerInfo('no_user'); # Reset to normal values in case used later
my @serverPath=qw(cgi html images java_gel2d data_html bin R_scripts shared);
#my $cgiPath=$promsPath{'cgi'}; # Needed unmodified for later LWP test
my $dir=`pwd`; chomp ($dir);
$dir =~ s/\/cgi.*//i;
foreach my $pathName (@serverPath) {
	if ($pathName=~/^(html|images|java_gel2d|\w+_html)\Z/) {
		my $unixRoot=$dir;
		if ($ENV{'DOCUMENT_ROOT'} eq '/Library/WebServer/Documents') { # Mac
			$unixRoot.='/Documents' unless $promsUnixPath{$pathName}=~/\/Documents/; # html is silent in URL
		}
		else { # Unix, Linux
			$unixRoot.='/html' unless $promsUnixPath{$pathName}=~/\/html/; # html is silent in URL
		}
		$promsUnixPath{$pathName}=~s/\/myproms// if $dir=~/myproms\Z/;
		$promsUnixPath{$pathName}=$unixRoot.$promsUnixPath{$pathName};
	}
	elsif ($pathName eq 'cgi') {
		$promsUnixPath{$pathName}=~s/\/myproms// if $dir=~/myproms\Z/;
		$promsUnixPath{$pathName}=~s/cgi-bin/CGI-Executables/ if $ENV{'DOCUMENT_ROOT'} eq '/Library/WebServer/Documents'; # Mac
		$promsUnixPath{$pathName}="$dir$promsPath{$pathName}";
	}
	elsif ($promsUnixPath{$pathName}) { # eg shared
		if (-e $promsUnixPath{$pathName}) {
			print "&nbsp;&nbsp;-Path for '<B>$pathName</B>' is found ($promsUnixPath{$pathName})<BR>\n";
		}
		else {print "&nbsp;&nbsp;<FONT color=red>Error: Path for '<B>$pathName</B>' was not found!!! ($promsUnixPath{$pathName})</FONT><BR>\n";}
	}
}
print "<B>Done.</B><BR><BR>\n";

print "<B>>Checking local binaries...</B><BR>\n";
foreach my $pathName (sort{lc($a) cmp lc($b)} @execPath) {
	next unless defined $promsPath{$pathName}; # in case exec no longer used
#print "*$pathName*<BR>";
	if ($pathName eq 'java') {
		my ($javaCommand,$javaPath);
		if ($promsPath{$pathName}) {
			$javaCommand="$promsPath{$pathName}/java";
			$javaPath=$promsPath{$pathName};
		}
		else {
			$javaCommand='java';
			$javaPath=`which java`;
			$javaPath=~s/\/java\s*\Z//;
		}
		my @response=`$javaCommand -showversion 2>&1`;
		my ($version)=($response[0] && $response[0]=~/ version "([^"]+)"/);
		&printVersionResponse($pathName,$javaPath,$version,\@response);
	}
	elsif ($pathName eq 'R') {
		my ($RCommand,$RPath);
		if ($promsPath{$pathName}) {
			$RCommand="$promsPath{$pathName}/R";
			$RPath=$promsPath{$pathName};
		}
		else {
			$RCommand='R';
			$RPath=`which R`;
			$RPath=~s/\/R\s*\Z//;
		}
		my @response=`$RCommand CMD BATCH --version 2>&1`;
		my ($version)=($response[0] && $response[0]=~/R batch front end: (\S+)/);
		&printVersionResponse($pathName,$RPath,$version,\@response);
	}
	elsif ($pathName eq 'qvality') {
		my ($qvalityCommand,$qvalityPath);
		if ($promsPath{$pathName}) {
			$qvalityCommand="$promsPath{$pathName}/qvality"; # no space before &1
			$qvalityPath=$promsPath{$pathName};
		}
		else {
			$qvalityCommand='qvality';
			$qvalityPath=`which qvality`;
			$qvalityPath=~s/\/qvality\s*\Z//;
		}
		my @response=`$qvalityCommand -h 2>&1`;
		my ($version)=($response[0] && $response[0]=~/qvality version ([\w\.-]+)/);
		&printVersionResponse($pathName,$qvalityPath,$version,\@response);
	}
	elsif ($pathName eq 'masschroq') {
		my ($masschroqCommand,$masschroqPath);
		if ($promsPath{$pathName}) {
			$masschroqCommand="$promsPath{$pathName}/masschroq";
			$masschroqPath=$promsPath{$pathName};
		}
		else {
			$masschroqCommand='masschroq';
			$masschroqPath=`which masschroq`;
			$masschroqPath=~s/\/masschroq\s*\Z//;
		}
		my @response=`$masschroqCommand -v 2>&1`;
		#my ($version)=($response[0] && $response[0]=~/MassChroQ version (\S+),/i);
		my ($version)=($response[0] && $response[0]=~/MassChroQ [version ]*([\d\.]+)/i);
		&printVersionResponse($pathName,$masschroqPath,$version,\@response);
	}
	elsif ($pathName eq 'tpp') {
		my @response=`$promsPath{tpp}/xinteract | grep TPP 2>&1`;
		my ($version)=($response[0] && $response[0]=~/\(TPP ([^,]+)/);
		&printVersionResponse($pathName,$promsPath{$pathName},$version,\@response);
	}
	elsif ($pathName eq 'python') {
		my ($pythonCommand,$pythonPath);
		if ($promsPath{'python'}) {
			$pythonPath=$promsPath{'python'};
			$pythonCommand="$pythonPath/python";
		}
		else {
			$pythonCommand='python';
			$pythonPath=`which python`;
			$pythonPath=~s/\/python\s*\Z//;
		}
		my @response=`$pythonCommand -V 2>&1`;
		my ($version)=($response[0] && $response[0]=~/Python (.+)/);
		&printVersionResponse($pathName,$pythonPath,$version,\@response);
	}
	elsif ($pathName eq 'msproteomicstools') {
		#1. Try pip
		my @response=`pip show msproteomicstools 2>&1 | grep -P "Version:"`;
		my ($version)=($response[0] && $response[0]=~/Version:\s*(\S+)/); # pip
		#2. Try spectrast2tsv.py
		unless ($version) {
			@response=`grep msproteomicstools $promsPath{$pathName}/spectrast2tsv.py 2>&1`;
			my $found;
			foreach my $line (@response) {
				if ($line=~/^\s*msproteomicstools/) {
					$found=1;
				}
				if ($line=~/msproteomicstools==([^']+)/) {
					$version=$1;
					last;
				}
			}
			$version='Unknown' if (!$version && $found);
		}
		&printVersionResponse($pathName,$promsPath{$pathName},$version,\@response);
	}
	elsif ($pathName eq 'pyprophet') {
		#my @response=`pip show pyprophet 2>&1 | grep -P "Version:"`;
		#my ($version)=($response[0] && $response[0]=~/Version:\s*(\S+)/); # pip
		#unless ($version) {
		#	@response=`$promsPath{$pathName}/pyprophet --version 2>&1`;
		#	foreach my $line (@response) {
		#		($version)=($line && $line=~/^([\d\.]+)/);
		#		last if $version;
		#	}
		#}
		my @response=`$promsPath{$pathName}/pyprophet --version 2>&1`;
		my ($version)=($response[0] && $response[0]=~/([\d\.]+)/);
		&printVersionResponse($pathName,$promsPath{$pathName},$version,\@response);
	}
	elsif ($pathName eq 'openms') {
		my @response=`$promsPath{$pathName}/OpenMSInfo 2>&1`;
		my ($version)=($response[2] && $response[2]=~/Version\s*:\s*(\S+)/);
		&printVersionResponse($pathName,$promsPath{$pathName},$version,\@response);
	}
	#elsif ($pathName eq 'qsub') {
	#	print "&nbsp;&nbsp;+<B>qsub</B> detected ($promsPath{$pathName}). Jobs are run on <B>cluster</B><BR>\n";
	#}
	# !!! TODO: qsub check does not work!!!
	#elsif ($pathName eq 'qsub') {
	#	next unless $promsPath{$pathName};
	#	my @qsubAnswer=`$promsPath{$pathName}/qsub --version`;
	#	my ($version)=($qsubAnswer[0]=~/Version: (\S+)/);
	#	if ($version) {print "&nbsp;&nbsp;-<B>qsub</B> $version found ($promsPath{$pathName})<BR>\n";}
	#	else {print "&nbsp;&nbsp;<FONT color=red>Error: <B>qsub</B> not found!!! (in '$promsPath{$pathName}')</FONT><BR>\n";}
	#}
}
print "<B>Done.</B><BR><BR>\n";

my %clusterInfo=&promsConfig::getClusterInfo; # exec on cluster
if ($clusterInfo{'on'} || ($clusterInfo{'list'} && scalar @{$clusterInfo{'list'}})) { # clusterInfo{list} is not mandatory in promsConfig.pm
	my $clusterStrg='cluster';
	$clusterStrg='clusters' if ($clusterInfo{'list'} && scalar @{$clusterInfo{'list'}} > 1);
	@{$clusterInfo{'list'}}=($clusterInfo{'name'}) if (!$clusterInfo{'list'} || scalar @{$clusterInfo{'list'}}==0);
	print qq
|<B>>Check binaries on $clusterStrg:</B><BR>
<TABLE><TR>
|;
	foreach my $i (0..$#{$clusterInfo{'list'}}) {
		my $refClusInfo={};
		if ($i==0) {$refClusInfo=\%clusterInfo;}
		else {%{$refClusInfo}=&promsConfig::getClusterInfo($clusterInfo{'list'}[$i]);}
		next unless $refClusInfo->{'on'};
		my $singularityImageStrg=($refClusInfo->{singularityImage})? " (using Singularity image '".$refClusInfo->{singularityImage}."')" : '';
		print "<TD valign='top'>&nbsp;&nbsp;<B>$refClusInfo->{name}$singularityImageStrg :</B>&nbsp;<INPUT type=\"button\" id=\"cluster${i}BUT\" value=\"Check\" onclick=\"checkBinaries($i,'$refClusInfo->{name}')\"/>&nbsp;&nbsp;<DIV id=\"cluster${i}DIV\"></DIV></TD>\n";
		print "<TD>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</TD>\n" if $i < $#{$clusterInfo{'list'}};
	}
	print "</TR></TABLE><BR>\n";
}
##>Cleaning test directory
if ($promsPath{'tmp'} && -d "$promsPath{tmp}/test") {
	my $okDir=1;
	opendir (DIR,"$promsPath{tmp}/test") or $okDir=0;
	if ($okDir) {
		my $now=strftime("%Y%m%d%H%M%S",localtime);
		while (defined (my $dir = readdir (DIR))) {
			#next if (!-d $dir || $dir =~/^\./);
			next unless $dir=~/\d{14}/;
			File::Path->remove_tree("$promsPath{tmp}/test/$dir") if $now-$dir >= 5000; # ~1/2h
		}
	}
	close DIR;
}

####>Data files<####
print "<B>>Checking other files...</B><BR>\n";
foreach my $refFileInfo ([$promsPath{data},'unimod_tables.xml'],[$promsPath{banks},'taxonomy_codes.txt'],[$promsPath{bin},'phosphoRS.jar']) {
	if (-e "$refFileInfo->[0]/$refFileInfo->[1]") {print "&nbsp;&nbsp;-<B>$refFileInfo->[1]</B> found ($refFileInfo->[0])<BR>\n";}
	else {print "&nbsp;&nbsp;<FONT color=red>Error: <B>$refFileInfo->[1]</B> not found in '$refFileInfo->[0]'</FONT><BR>\n";}
}
print "<B>Done.</B><BR><BR>\n";


####>Checking web page progressive display<####
print "<B>>Checking web page display mode...</B><BR><SPAN id=\"progDispSPAN\"></SPAN><BR><BR>\n";


####>Checking connection to database<####
print "<B>>Checking database...</B><BR>\n";
my $dbh=&promsConfig::dbConnect('no_user');
if ($dbh) {
	print "&nbsp;&nbsp;-Connection to database is OK<BR>\n";
	##>Read access
	my ($okSelect)=$dbh->selectrow_array("SELECT COUNT(*) FROM DATABANK_TYPE");
	if ($okSelect) {
		print "&nbsp;&nbsp;-Read access is OK<BR>\n";
	}
	else {
		print "&nbsp;&nbsp;<FONT color=red>Error: No read access to database or data may be missing from the database. Make sure you have run myProMS_start_values.sql</FONT><BR>\n";
	}
	##>Write access
	$dbh->do("INSERT INTO USER_PROFILE (ID_PROFILE,NAME) VALUES (-1,'ZZZZ')");
	my ($okInsert)=$dbh->selectrow_array("SELECT COUNT(*) FROM USER_PROFILE WHERE ID_PROFILE=-1");
	if ($okInsert) {
		print "&nbsp;&nbsp;-Write access is OK<BR>\n";
		$dbh->do("DELETE FROM USER_PROFILE WHERE ID_PROFILE=-1");
	}
	else {
		print "&nbsp;&nbsp;<FONT color=red>Error: No write access to database. Upgrade the database user's priviledges</FONT><BR>\n";
	}
	$dbh->disconnect;
}
print "<B>Done.</B><BR><BR>\n";

####>Checking connection to Mascot Server<####
my %mascotServers=&promsConfig::getMascotServers;
if (scalar keys %mascotServers && $lwpOK) {
	print "<B>>Checking communication with Mascot server(s)...</B><BR>\n";
	my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
	$agent->timeout(360);
	foreach my $mascotServer (sort{lc($a) cmp lc($b)} keys %mascotServers) {
		print "&nbsp;&nbsp;<B>&bull;$mascotServer</B> ($mascotServers{$mascotServer}[0]):<BR>\n";
		if ($mascotServers{$mascotServer}[1]) { # proxy settings
			if (lc($mascotServers{$mascotServer}[1]) eq 'no') {$agent->no_proxy($mascotServers{$mascotServer}[0]);}
			else {$agent->proxy('http', $mascotServers{$mascotServer}[1]);}
		}
		else {$agent->env_proxy;}
		#>Checking scripts
		foreach my $script ('myproms4datFiles.pl','myproms4databanks.pl','myproms4emPAI.pl') {
			my $response = $agent->get("$mascotServers{$mascotServer}[0]/cgi/$script?ACT=test");
			while (my $wait = $response->header('Retry-After')) {
				sleep $wait;
				$response = $agent->get($response->base);
			}
			my @mascotResponse;
			if ($response->is_success) {
				@mascotResponse = split("\n",$response->content);
			}
			my $responded=0;
			if ($mascotResponse[0]) {
				if ($mascotResponse[0]=~/# OK FROM (.+)/) {
					print "&nbsp;&nbsp;&nbsp;&nbsp;-$script on $1 is OK<BR>\n";
					$responded=1;
				}
				elsif ($mascotResponse[0]=~/.+DENIED.+HOST (.+)/) {
					print "<FONT color=red>&nbsp;&nbsp;&nbsp;&nbsp;-Server $1 is not allowed to connect to $mascotServer</FONT><BR>\n";
					$responded=1;
				}
			}
			print "<FONT color=red>&nbsp;&nbsp;&nbsp;&nbsp;-No response from $script</FONT><BR>\n" unless $responded;
		}
		#>Checking path to dat file if declared
		if ($mascotServers{$mascotServer}[2]) {
			if (-e "$mascotServers{$mascotServer}[2]/data") {
				print "&nbsp;&nbsp;&nbsp;&nbsp;-Local path to data is found: <B>$mascotServers{$mascotServer}[2]/data</B><BR>\n";
			}
			else {
				print "<FONT color=red>&nbsp;&nbsp;&nbsp;&nbsp;-Local path to data <B>not</B> found: <B>$mascotServers{$mascotServer}[2]/data</B></FONT><BR>\n";
			}
			if ($mascotServers{$mascotServer}[3]) {print "&nbsp;&nbsp;&nbsp;&nbsp;-Search files are <B>linked</B>, not copied<BR>\n";}
			else {print "&nbsp;&nbsp;&nbsp;&nbsp;-Search files are <B>copied</B>, not linked<BR>\n";}
		}
		else {print "&nbsp;&nbsp;&nbsp;&nbsp;-Search files are <B>copied</B><BR>\n";}
	}
	print "<B>Done.</B><BR><BR>\n";
}
elsif (scalar keys %mascotServers) {
	print "<B><FONT color=red>>Communication with Mascot server(s) cannot be checked (No LWP::UserAgent)</B></FONT><BR><BR>\n";
}
else {
	print "<B>>No Mascot servers declared.<B><BR><BR>\n";
}

####>Checking Internet connection to Uniprot<####
if ($lwpOK) {
	print "<B>>Check Internet connection: </B>";
	if ($promsPath{'http_proxy'} && $promsPath{'http_proxy'}) {
		my $configSettingStrg=(lc($promsPath{'http_proxy'}) eq 'no')? 'No proxy' : "Proxy: $promsPath{http_proxy}";
		print qq |<INPUT type="button" value="myProMS settings ($configSettingStrg)" onclick="window.location='./testMyProMS.cgi?ACT=internet_myproms&CALL=$call'"/>&nbsp;&nbsp;|;
	}
	print qq |<INPUT type="button" value="System's proxy" onclick="window.location='./testMyProMS.cgi?ACT=internet_system&CALL=$call'"/><BR><BR>|;
}
else {
	print "<B>>Cannot check Internet connection: LWP::UserAgent not found</B><BR><BR>\n";
}

####>Checking script versions<####
print qq
|<B>>Scripts:</B><INPUT type="button" id="scriptButton" value="Show versions" onclick="updateScriptDisplay()"/><BR>
<DIV id="scriptDiv" style="display:none">
<BR><FIELDSET><LEGEND><FONT style="font-weight:bold;font-size:24px">Scripts versions:</FONT></LEGEND>
<TABLE>
<TR><TD colspan=2><B>>HTML files ($promsUnixPath{html}):</B></TD></TR>
|;
opendir (HTML,$promsUnixPath{html}) || print "<TR><TD colspan=2>&nbsp;ERROR: Unable to read directory!</TD></TR>\n";
my @htmlFiles=readdir(HTML);
foreach my $currentFile (sort{lc($a) cmp lc($b)} @htmlFiles) {
	next if $currentFile !~ /\.(html|css)\Z/;
	&printScriptVersion($promsUnixPath{html},$currentFile);
}
close HTML;

print "<TR><TD colspan=2>&nbsp;</TD></TR><TR><TD colspan=2><B>>JavaScript files ($promsUnixPath{html}/js/local):</B></TD></TR>\n";
opendir (JS,"$promsUnixPath{html}/js/local") || print "<TR><TD colspan=2>&nbsp;ERROR: Unable to read directory!</TD></TR>\n";
my @jsFiles=readdir(JS);
foreach my $currentFile (sort{lc($a) cmp lc($b)} @jsFiles) {
	next if $currentFile !~ /\.js\Z/;
	&printScriptVersion("$promsUnixPath{html}/js/local",$currentFile);
}
close JS;

print "<TR><TD colspan=2>&nbsp;</TD></TR><TR><TD colspan=2><B>>Perl scripts ($promsUnixPath{cgi}):</B></TD></TR>\n";
opendir (CGI,$promsUnixPath{cgi}) || print "<TR><TD colspan=2>&nbsp;ERROR: Unable to read directory!</TD></TR>\n";
my @cgiFiles=readdir(CGI);
foreach my $currentFile (sort{lc($a) cmp lc($b)} @cgiFiles) {
	next if $currentFile !~ /\.(cgi|pl|pm)\Z/;
	&printScriptVersion($promsUnixPath{cgi},$currentFile);
}
close CGI;
print "<TR><TD colspan=2>&nbsp;</TD></TR><TR><TD colspan=2><B>>R scripts ($promsUnixPath{R_scripts}):</B></TD></TR>\n";
opendir (R,$promsUnixPath{R_scripts}) || print "<TR><TD colspan=2>&nbsp;ERROR: Unable to read directory!</TD></TR>\n";
my @rFiles=readdir(R);
foreach my $currentFile (sort{lc($a) cmp lc($b)} @rFiles) {
	next if $currentFile !~ /\.R\Z/;
	&printScriptVersion($promsUnixPath{R_scripts},$currentFile);
}
close R;
print "</TABLE></FIELDSET></DIV>\n";
sleep 1;
print qq
|<H3>Test is complete.</H3>
<BR><BR><BR>
<SCRIPT LANGUAGE="JavaScript">
var endDate = new Date();
var progSpan=document.getElementById('progDispSPAN');
if (endDate.getTime() - startDate.getTime() > 2000) {progSpan.innerHTML='&nbsp;&nbsp;<B>Progressive</B> display detected<BR><B>Done.</B>';}
else {progSpan.innerHTML='&nbsp;&nbsp;<FONT color=red>Error: Server is not configured for <B>progressive</B> display</FONT><BR><B>Done.</B>';}
</SCRIPT>
</BODY>
</HTML>
|;

sub testWriteDirectory {
	my ($dir)=@_;
	open (TEST,">$dir/testFile.txt") || print "$!<BR>\n";
	close TEST;
	my $ok=0;
	if (-e "$dir/testFile.txt") {
		print "&nbsp;&nbsp;&nbsp;&nbsp;-Write permission is OK<BR>\n";
		unlink "$dir/testFile.txt";
		$ok=1;
	}
	else {print "&nbsp;&nbsp;&nbsp;&nbsp;<FONT color=red>Error: Could not write into directory!</FONT><BR>\n";}
	return $ok;
}

sub printVersionResponse {
	my ($execName,$execPath,$version,$refResponse,$clusterIdx,$refSessionPackages)=@_;
	if ($version) {
		print "&nbsp;&nbsp;-<B>$execName</B> $version found ($execPath)"; # if (!defined($clusterIdx) || $clusterIdx >= 0); # Not for local R packages check
		if ($execName eq 'R') {
			if (defined($clusterIdx) && $clusterIdx >= 0) { # $refSessionPackages
				print ". <B>",scalar @Rpackages," packages</B> required <INPUT type=\"button\" value=\"Show/Hide packages\" onclick=\"showHideRPackages($clusterIdx)\"><DIV id=\"Rpackages_",$clusterIdx,"DIV\" style=\"display:none\">\n";
				&printRpackageVersion($refSessionPackages);
				print "</DIV><BR>\n";
			}
			else {
				print ". <B>",scalar @Rpackages," packages</B> required <INPUT type=\"button\" id=\"RpackagesBUT\" value=\"Check packages\" onclick=\"checkBinaries(-1)\"><DIV id=\"RpackagesDIV\"></DIV>\n";
				##foreach my $i (0..$#Rpackages) {
				##	print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>-$Rpackages[$i]:</B><SPAN id=\"SPAN_Rpackage$i\"></SPAN><BR>\n";
				##	#print "<BR>" if $i < $#Rpackages;
				##}
			}
		}
		else {print "<BR>\n";}
	}
	else {
		print "&nbsp;&nbsp;<FONT color=red>-Error detected for <B>$execName</B>!!! (in '$execPath')<BR>\n&nbsp;&nbsp;<B>Server responded:</B><BR>&nbsp;&nbsp;&nbsp;&nbsp;";
		if ($?==-1) {print $!;}
		else {print join('<BR>&nbsp;&nbsp;&nbsp;&nbsp;',@{$refResponse});}
		print "</FONT><BR>\n";
	}
}

sub printRpackageVersion {
	my ($refSessionPackages)=@_;
	foreach my $package (@Rpackages) {
		print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>-$package:</B>";
		if ($refSessionPackages->{$package}) {
			print "&nbsp;Version $refSessionPackages->{$package} found<BR>\n";
		}
		else {
			print "&nbsp;<FONT color=red>Error: Package not found!</FONT><BR>\n";
		}
	}
}

sub printScriptVersion {
	my ($dir,$file)=@_;
	my $fileName=(split(/[\/\\]/,$file))[-1];
	my $scriptVersionLine=`grep "# $fileName " $dir/$fileName`;
	my ($version)=($scriptVersionLine=~/ (\d\S+)/);
	$version='-' unless $version;
	print "<TR><TD>&nbsp;-$fileName:</TD><TD>$version</TD></TR>\n";
}


sub checkInternetConnection {
	require LWP::UserAgent;
	require Net::FTP;
	require promsConfig;
	%promsPath=&promsConfig::getServerInfo('no_user') unless $call eq 'server';

	my $proxySettingStrg=($action eq 'internet_system')? "System's proxy" : (lc($promsPath{'http_proxy'}) eq 'no')? 'No proxy' : "Proxy: $promsPath{http_proxy}";
	print "</HEAD>\n";
	if ($call eq 'server') {print qq|<BODY background="$promsPath{images}/bgProMS.gif">\n|;}
	else {print "<BODY>\n";}
	print qq
|<BR>
<INPUT type="button" value=" << Back" onclick="window.location='./testMyProMS.cgi?CALL=$call'"/>
<BR>
<H3>Checking myProMS' Internet connection with $proxySettingStrg:</H3>
|;
	my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
	#$agent->timeout(30);
	if ($action eq 'internet_myproms') {
		if (lc($promsPath{'http_proxy'}) eq 'no') {$agent->no_proxy;}
		else {$agent->proxy('http', $promsPath{'http_proxy'});}
	}
	else {$agent->env_proxy;} # system's proxy

	my $urlHTTP='http://www.uniprot.org/mapping?from=ACC&to=ID&query=P15311&format=tab';
	print "&nbsp;&nbsp;<B>+HTTP:</B> ($urlHTTP)<BR>\n";
	my $response = $agent->get($urlHTTP);
	while (my $wait = $response->header('Retry-After')) {
		sleep $wait;
		$response = $agent->get($response->base);
	}
	my @internetResponse;
	if ($response->is_success) {
		@internetResponse = split("\n",$response->content);
		my ($acc,$id);
		foreach (@internetResponse) {
			next if /^From/;
			chomp;
			($acc,$id)=split(/\s+/,$_);
			last if $acc;
		}
		if ($acc) {print "&nbsp;&nbsp;&nbsp;&nbsp;-HTTP connection is <B>OK</B> (<B>$acc</B> was successfully mapped to <B>$id</B> by the Uniprot resource).<BR>\n";}
		else {print "<FONT color=red>&nbsp;&nbsp;&nbsp;&nbsp;-No response from http://www.uniprot.org.</FONT><BR>\n";}
	}
	else {
		print "<FONT color=red>&nbsp;&nbsp;&nbsp;&nbsp;-ERROR: ",$response->status_line,"</FONT><BR>\n";
	}

	###>FTP connection
	print "<BR><H3>Checking FTP connection:</H3>\n";
	my %ftpConfig=&promsConfig::getFTPconfig;
	my $ftpHost='ftp.ebi.ac.uk';
	my $remoteDir='/pub/databases/fastafiles/uniprot';
	my $localDir=$promsPath{'tmp'};
	my $file='README.txt';
	my $localFile="$localDir/$file";
	unlink $localFile if -e $localFile;
	my %okModes;

	foreach my $refFTPmode (['Active',0],['Passive',1]) {
		unlink $localFile if -e $localFile;
		my ($modeName,$modValue)=@{$refFTPmode};
		my $error='';
		print "&nbsp;&nbsp;<B>+$modeName</B> mode:<BR>\n";
		my $ftp = Net::FTP->new($ftpHost, Passive => $modValue) or $error="Cannot connect to $ftpHost: $@";
		unless ($error) {
			$ftp->login() or $error=$ftp->message; #login("anonymous",'-anonymous@') <- default
			unless ($error) {
				$ftp->cwd($remoteDir) or $error=$ftp->message." ($remoteDir)";
				unless ($error) {
					$ftp->get($file,$localFile) or $error=$ftp->message." ($remoteDir/$file)";
					if (-e $localFile) {
						unlink $localFile;
						print "&nbsp;&nbsp;&nbsp;&nbsp;-FTP in ",lc($modeName)," mode is <B>OK</B> (file '<B>$remoteDir/$file</B>' was successfully retrieved from '<B>$ftpHost</B>' server).<BR>\n";
						$okModes{lc($modeName)}=1;
					}
				}
			}
		}
		$ftp->quit if $ftp;
		if ($error) {
			print "<FONT color=red>&nbsp;&nbsp;&nbsp;&nbsp;-ERROR: $error</FONT><BR>\n";
		}
	}
	my $mypromsFTPmode=($ftpConfig{mode} && $ftpConfig{mode} =~ /passive/i)? 'passive' : 'active';
	if (!$okModes{$mypromsFTPmode} && scalar keys %okModes) {
		my ($otherMode)=(keys %okModes)[0];
		print "&nbsp;&nbsp;<FONT color=red><B>+WARNING:</B></FONT> FTP mode is set to <B>$mypromsFTPmode</B> in promsConfig.pm. You must change it to <B>$otherMode</B> mode.<BR>\n";
	}
	else {
		print "&nbsp;&nbsp;+FTP mode is set to <B>$mypromsFTPmode</B> in promsConfig.pm<BR>\n";
	}
	
	print qq
|<H3>Done.</H3><BR><BR>
</BODY>
</HTML>
|;
	exit;
}


####>Revision history<####
# 2.5.2 [FIX] Bugs in version detection of local java & pyprophet (PP 24/06/19)
# 2.5.1 Added Singularity image info for cluster (PP 19/06/19)
# 2.5.0 Uses Net::FTP for testing default and passive FTP connection (PP 06/05/19)
# 2.4.14 Updated to docker image 1.2.x (PP 29/03/19)
# 2.4.13 [Fix] myProMS version detection for docker instance (PP 07/10/18)
# 2.4.12 Updated parsing rule for MassChroQ version also in local binaries section (PP 06/10/18)
# 2.4.11 Writes Perl $0 and `pwd` (PP 03/10/18)
# 2.4.10 $clusterInfo{'list'} is no longer mandatory if only 1 cluster is defined in promsConfig.pm (PP 27/09/18)
# 2.4.9 New parsing rule for MassChroQ version & uses $cluster{runJob} for launching cluster job (PP 17/09/18)
# 2.4.8 Minor improvement in version number detection (PP 23/04/18)
# 2.4.7 Extracts myProMS version number from index.html (PP 20/04/18)
# 2.4.6 Stricter check for pyprophet (PP 16/04/18)
# 2.4.5 Fix checks for msproteomicstools and pyprophet in local binaries context (PP 08/12/17)
# 2.4.4 Checks pyprophet & FTP URL with FTP and HTTP protocols (PP 05/12/17)
# 2.4.3 Minor change in cluster error management (PP 07/11/17)
# 2.4.2 Checks binaries & R packages on cluster(s) (PP 20/10/17)
# 2.4.1 Update for MassChroQ 2.2.x & OpenMS (PP 11/10/17)
# 2.4.0 New R packages with better AJAX check & FTP connection test (PP 28/08/17)
# 2.3.9 Minor fixes in display (PP 20/01/17)
# 2.3.8 Tests more "other" files (PP 26/10/16)
# 2.3.7 Bug fix in msproteomicstools path detection (PP 20/10/16)
# 2.3.6 Better test for msproteomicstools & checks web server progressive page display (PP 30/09/16)
# 2.3.5 Split tests by %promsPath directories types & MSstats package & java (PP 20/09/16)
# 2.3.4 Checks list of required Perl modules and R packages for myproMS v3.5 (PP 11/07/16)
# 2.3.3 Added test for tpp, python & msproteomicstools (PP 29/04/16)
# 2.3.2 Color on main "Go to" options (PP 25/03/16)
# 2.3.1 Added tmp_html path (PP 28/04/15)
# 2.3.0 Better handling of errors on binaries (PP 19/05/14)
# 2.2.9 Checks for List::Util module (PP 16/05/14)
# 2.2.8 Checks for Storable module (PP 10/03/14)
# 2.2.7 Display modules version & removal of local URL reconstruction test (PP 08/03/14)
# 2.2.6 Checks for File::Copy::Recursive module (PP 08/11/13)
# 2.2.5 Also checks for Mascot local data directory if declared (PP 28/10/13)
# 2.2.4 Added check of script absolute URL reconstruction (PP 23/10/13)
# 2.2.3 Added check for unimod_tables.xml file (PP 23/09/13)
# 2.2.2 Added check for myproms4emPAI.pl on Mascot server (PP 19/09/13)
# 2.2.1 Can be called from myProMS server (PP 26/09/13)
# 2.2.0 Uses LWP to test internet connection (PP 01/07/13)
# 2.1.9 Updated to test v3.0 (PP 24/06/13)
# 2.1.8 No write check on web server directories (PP 15/10/12)
# 2.1.7 User-triggered Internet connection check (PP 20/06/12)
# 2.1.6 Added check on internet connection (PP 29/05/12)
# 2.1.5 header(-'content-encoding'=>'no') (PP 17/10/11)
# 2.1.4 added better handling of mascot server error (PP 28/06/11)
# 2.1.3 added &lt;BR&gt;\n to Perl module error string (PP 16/06/11)
