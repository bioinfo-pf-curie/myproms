#!/usr/local/bin/perl -w
################################################################################
# myproms4datFiles.pl        1.0.9                                             #
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
$|=1;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use strict;


##################
####>Security<####
##################
my @remoteHostIPs=( # Uncomment (#) and add comma (,) after each server but last
	'10.2.0.30', # pangaea
	'10.2.200.39', # dev
	'10.2.0.193', # new prod
	'10.2.200.38', # dev/prod on bi-web02
	'10.200.10.172', # BIWS ppoullet
	'10.200.10.93' # BIWS fyvon
);
if (!param('SKIP') && scalar @remoteHostIPs && $ENV{'REMOTE_ADDR'}) {
	my $okRequest=0;
	foreach my $allowedIP (@remoteHostIPs) {
		if ($allowedIP eq $ENV{'REMOTE_ADDR'}) {
			$okRequest=1;
			last;
		}
	}
	unless ($okRequest) {
		print header(-type=>'text/plain'); warningsToBrowser(1);
		print "# REQUEST DENIED FOR HOST $ENV{REMOTE_ADDR}\n";
		exit;
	}
}

####################
####>Parameters<####
####################
my $action=param('ACT') || 'test';

if ($action ne 'get') { # different header for 'get'
	print header(-type=>'text/plain'); warningsToBrowser(1);
}

##############
####>Test<####
##############
if ($action eq 'test') {
	print "# OK FROM $ENV{SERVER_NAME}\n";
	exit;
}

#################################
####>Mascot server variables<####
#################################
elsif ($action eq 'mascotVar') {
	print qq |
<B>Script variable:</B><BR>
-Perl \$0: $0<BR><BR>
<B>Known paths (\@INC):</B><BR>
|;
	foreach my $v (@INC) {print "$v<BR>\n";}
	print "<BR><B>Environment variables (\%ENV):</B><BR>\n";
	foreach my $v (sort{lc($a) cmp lc($b)} keys %ENV) {print "$v: $ENV{$v}<BR>\n";}
	exit;
}

#############################
####>Listing directories<####
#############################
elsif ($action eq 'log') {
	print "# LOG\n";
	opendir(LOGDIR,'../logs');
	while ((my $file=readdir(LOGDIR))) {
		print "$file\n" if $file =~ /search.*\.log\Z/;
	}
	close LOGDIR;
	exit;
}

#############################
####>Listing directories<####
#############################
elsif ($action eq 'list') {
	my @logList=split(':',param('LOG'));
	my ($startDate,$endDate)=(param('DR'))? split(':',param('DR')) : (0,21001231);
	my ($startJob,$endJob)=(param('JR'))? split(':',param('JR')) : (0,999999);
	my $searchStrg=(param('STRG'))? quotemeta(param('STRG')) : '';
	my @userIDs=(param('IDS'))? split(',',param('IDS')) : (0);

	print "# LIST\n";
	my %jobNumbers;
	foreach my $searchLogFile (@logList) {
		open(LOG,"../logs/$searchLogFile"); # || die "Can't open Log file\n";
		while (<LOG>) {
			next if $.==1; # skip 1st line
			next unless /^\d+/; # skip empty line
			s/\s+\Z//; # chomp is not always enough!
			my @searchData=split(/\t/,$_);
			next if $jobNumbers{$searchData[0]}; # in case duplicates in different search logs
			if ($startJob) {
				next if $searchData[0]=~/\D/; # bad jobID
				next if $searchData[0]<$startJob;
				last if $searchData[0]>$endJob;
			}
			if ($startDate) {
				my ($dateDir)=($searchData[6]=~/(\d{8})/);
				next unless $dateDir; # non-numerical date
				next if $dateDir<$startDate;
				last if $dateDir>$endDate;
			}
			next if ($searchStrg && (!$searchData[5] || $searchData[5] !~ /$searchStrg/i));
			my $matchedID=0;
			#if ($searchData[-1]==0) { # no ID management
			#	$matchedID=1;
			#}
			#else {
				foreach my $ID (@userIDs) {
					if ($ID==0 || $ID==$searchData[-1]) { # 0: massist or bioinfo
						$matchedID=1;
						last;
					}
				}
			#}
			next unless $matchedID;
			$jobNumbers{$searchData[0]}=1;
			my $existFile=(-e $searchData[6])? 1 : 0;
			print "$searchData[6]\t$searchData[-1]\t$existFile\t$searchData[5]\n";
		}
		close LOG;
	}
	exit;
}

##########################
####>Sending dat file<####
##########################
elsif ($action eq 'get') {
	my $selFile=param('DAT');
	my $popFile=param('POP');
	if ($popFile) {
		my @infos=split(/\//,$selFile);
		my @date=split('',$infos[2]);
		my ($year,$month)=(join('',($date[0],$date[1],$date[2],$date[3])),join('',($date[4],$date[5])));
		#print "../data/$year/$month\n";
		opendir(CACHEDIR,"../data/cache/$year/$month");
		my $foundFile=0;
		while ((my $cacheDir=readdir(CACHEDIR))) {
			next if $cacheDir =~ /^\.\.?$/; # skip . and ..
			opendir(DATDIR,"../data/cache/$year/$month/$cacheDir") if -d "../data/cache/$year/$month/$cacheDir";
			while ((my $file=readdir(DATDIR))) {
				next if $cacheDir =~ /^\.\.?$/; # skip . and ..
				if ($file =~ /$infos[3].*\.$popFile\.pop/) {
					my $popFile="../data/cache/$year/$month/$cacheDir/$file";
					my $fileSize=-s $popFile;
					print header(-type=>'text/plain',-'content-length'=>$fileSize); warningsToBrowser(1);
					open(POP,$popFile);
					while (<POP>) {
						print $_;
					}
					close POP;
					$foundFile=1;
					last;
				}
			}
			close DATDIR;
			last if $foundFile;
		}
		close CACHEDIR;
	}
	else {
		# my $fileSize=-s $selFile;
		# print header(-type=>'text/plain',-'content-length'=>$fileSize); warningsToBrowser(1);
		print header(-type=>'text/plain'); warningsToBrowser(1);
		open(DAT,$selFile); # || die "Can't open file $selFile\n";
		while (<DAT>) {
			print $_;
		}
		close DAT;
	}
	exit;
}

####>Revision history<####
# 1.0.9 [BUGFIX] Removed content-length for file retrieval to prevent myProMS timeout (PP 02/07/21)
# 1.0.8 [FEATURE] Added ACT=mascotVar to retrieve Mascot server variables and check on $ENV{REMOTE_ADDR} & content-length for file retrieval (PP 28/06/21)
# 1.0.7 Send Mascot *.pop files created in the cache folder when Percolator is run (GA 22/09/16)
# 1.0.6 Cleaner line ending during search log scan (PP 28/10/13)
# 1.0.5 GPL license (PP 19/09/13)
# 1.0.4 management of Mascot User IDs (PP 30/03/11)
