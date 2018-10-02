#!/usr/local/bin/perl -w
################################################################################
# myproms4databanks.pl        1.0.8D                                            #
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
#use URI::Escape;
use strict;

print header(-type=>'text/plain'); warningsToBrowser(1);

##################
####>Security<####
##################
my @remoteHostIPs=(
	# '<IP address of myProMS web server>' (Add comma (,) after each entry but last)
);
if (scalar @remoteHostIPs) {
	my $okRequest=0;
	foreach my $allowedIP (@remoteHostIPs) {
		if ($allowedIP eq $ENV{'REMOTE_ADDR'}) {
			$okRequest=1;
			last;
		}
	}
	unless ($okRequest) {
		print "# REQUEST DENIED FOR HOST $ENV{REMOTE_ADDR}\n";
		exit;
	}
}

####################
####>Parameters<####
####################
my $action=(param('ACT'))? param('ACT') : 'test';

##############
####>Test<####
##############
if ($action eq 'test') {
	print "# OK FROM $ENV{SERVER_NAME}\n";
	exit;
}

###########################
####>Listing databanks<####
###########################
elsif ($action eq 'list') {
	open (MSDAT,"../config/mascot.dat");
	my %databanks;
	my $okDB=0;
	while (<MSDAT>) {
		next unless /\S/;
		next if /^#/;
		if (/^Databases/) {$okDB=1; next;}
		next unless $okDB;
		last if ($okDB && /^end/);
		my ($dbName,$filePathCode,$dbType,@dbData)=split(/\t/,$_);
		next if $dbType eq 'NA'; # Nucleic Acid

		##>Fetching real file name
		my $trueDbFile=&getTrueFile($filePathCode);
		$databanks{$dbName}=$trueDbFile if $trueDbFile;
	}
	close MSDAT;
	foreach my $dbName (sort{lc($a) cmp lc($b)} keys %databanks) {
		print "$dbName\t$databanks{$dbName}\n";
	}
	exit;
}

#######################################
####>Updating or adding a databank<####
#######################################
elsif ($action eq 'dbFile') {
	my $trueDbFile=&getDatabankFile(param('DB'));
	my $newFileName=(split(/\/|\\/,$trueDbFile))[-1];
	my $numEntries=0;
	if (!param('FILE') || $newFileName ne param('FILE')) { # no param('FILE') if adding a new DB to myProMS
		open (FAS,$trueDbFile);
		while (<FAS>) {
			if (/^>/) {
				$numEntries++;
				print "#$numEntries\n" unless $numEntries % 100000;
			}

		}
		close FAS;
	}
	print "$trueDbFile\t$numEntries";
	exit;
}

#############################
####>Testing parse rules<####
#############################
elsif ($action eq 'parse') {
	my $trueDbFile=&getDatabankFile(param('DB'));
	print "FILE=$trueDbFile\n";
	my ($idRule,$desRule,$orgRule)=&getParseRules(param('parseRules'));
	my $maxRange=(!param('numEntries'))? 10 : (param('numEntries') <= 1000)? param('numEntries') : 1000;
	my %entryPos;
	if ($maxRange==10) {
		foreach my $i (1..10) {$entryPos{$i}=1;} # pick first 10 positions
	}
	else {
		foreach my $i (1..10) {$entryPos{int(rand($maxRange)+0.5)}=1;} # pick 10 random positions
	}
	open (DB,$trueDbFile) || die "can't open $trueDbFile\n";
	my $curPos=0;
	while (<DB>) {
		last if scalar keys %entryPos==0;
		if (/^>/) {
			$curPos++;
			next unless $entryPos{$curPos};
			(my $entry=$_)=~s/^>\s*//;
			$entry =~s /.+//; # take 1 entry if multi-entry line (NCBI)
			$entry =~s /\s+\Z//; # chomp not always enough
			my ($identifier)=($entry=~/$idRule/);
			my $des; ($des)=($entry=~/$desRule/) if $desRule; # optional
			my $org; ($org)=($entry=~/$orgRule/) if $orgRule; # optional
			if ($org) {
				($des)=($entry=~/^\S+\s+(.+)$orgRule/) unless $desRule;
				$org=~s/\s+\Z//; # removing trailing spaces
			}
			else {
				($des)=($entry=~/^\S+\s+(.+)/) unless $desRule;
			}
			$des='' unless $des;
			$des=~s/\s+\Z//; # removing trailing spaces
			$org='' unless $org;
			print "$identifier ##DES=$des ##ORG=$org ##HEAD=$entry\n";
			delete $entryPos{$curPos};
		}
	}
	close DB;
	exit;
}

#############################
####>Scanning a databank<####
#############################
elsif ($action eq 'scan') {
	my %massAAave=( # Average mass, needed for protein mass calculation
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
	my %massATave=( # Average mass, needed for protein mass calculation
		H=>1.00794,
		C=>12.011,
		N=>14.0067,
		O=>15.9994
	);
	my $trueDbFile=&getDatabankFile(param('DB'));
	my (%protList,%protMW,%protLength,%protDes,%protOrg);
	foreach my $identifier (split(':',param('protList'))) {$protList{$identifier}=1;}
#my $dbOrganism=param('dbOrganism');
#	my $parseRules=uri_unescape(param('parseRules'));#deprotection of url information
##$parseRules=~s/�/\+/g;
#	my @rules=split(',:,',$parseRules);
##my @rules=split(',:,','ID=([^\|]+_[^\s\|]+),:,ORG=\s?- ([^\(]*),:,');
#	my ($idRule)=($rules[0]=~/ID=(.+)/);
#	my ($desRule,$orgRule);
#	if ($rules[1]) {
#		if ($rules[1]=~/DES=/) {($desRule=$rules[1])=~s/DES=//;}
#		else {($orgRule=$rules[1])=~s/ORG=//;}
#	}
#	if ($rules[2]) {($orgRule=$rules[2])=~s/ORG=//;}
	my ($idRule,$desRule,$orgRule)=&getParseRules(param('parseRules'));
	my $identType=param('identType') || '';

	my $maxProtID=scalar keys (%protList);

	###>Scanning DBank file
	my $counter=0;
	my $sequence='';
	my $matched;
	my @matchedID;
	open (DB,$trueDbFile) || die "can't open $trueDbFile\n";
	while (<DB>) {
		if (/^>/) {
			if ($matched) {
				my %countAA;
				foreach my $aa (split(//,uc($sequence))) {
					$countAA{$aa}++;
				}
				my $mass=0;
				foreach my $aa (keys %countAA) {
					$mass+=($massAAave{$aa}*$countAA{$aa}) if $massAAave{$aa}; # some characters are not amino-acids
				}
				$mass+=($massATave{H}+$massATave{H}+$massATave{O}) if $mass; # H(Nter) + OH(Cter)
				$mass=sprintf "%.2f",$mass; # no need for more precision

				###>Writing entry
				my @entryStrg;
				foreach my $identifier (@matchedID) {
					$protDes{$identifier}='no description' unless $protDes{$identifier};
					$protOrg{$identifier}='unknown organism' unless $protOrg{$identifier};
					$protMW{$identifier}=($mass)? $mass : 0;
					$protLength{$identifier}=($sequence)? length($sequence) : 0;
					push @entryStrg,"$identifier ##DES=$protDes{$identifier} ##ORG=$protOrg{$identifier} ##MW=$protMW{$identifier} ##LEN=$protLength{$identifier}";
				}
				print ">",join('',@entryStrg),"\n$sequence\n";
				@matchedID=(); %protDes=(); %protOrg=(); %protMW=(); %protLength=();
				$sequence='';
				$matched=0;
				last if $counter==$maxProtID; # all identifier were matched
			}
			chomp;
			(my $newLine=$_)=~s/^>\s*//;
			my @line=split(//,$newLine);
			foreach my $entry (@line) {
				my ($identifier)=($entry=~/$idRule/);
				next unless $identifier; # just to be safe

				if ($identType eq 'UNIPROT_ID') { # check for isoforms in 1st keyword before | & add to UNIPROT_ID
					$entry=~s/^sp\|//;
					if ($entry=~/^[^|]+-(\d+)\|/) {
						$identifier.="-$1";
					}
				}

				if (defined($protList{$identifier})) {
					my $des; ($des)=($entry=~/$desRule/) if $desRule; # optional
					my $org; ($org)=($entry=~/$orgRule/) if $orgRule; # optional
					if ($org) {
						($des)=($entry=~/^\S+\s+(.+)$orgRule/) unless $desRule;
						$org=~s/\s+\Z//; # removing trailing spaces
						$protOrg{$identifier}=$org;
					}
					else {
						($des)=($entry=~/^\S+\s+(.+)/) unless $desRule;
						#$protOrg{$identifier}=$dbOrganism if $dbOrganism; # Moved to calling script
					}
					$des='' unless $des;
					($protDes{$identifier}=$des)=~s/\s+\Z//; # removing trailing spaces

					$matched=1;
					push @matchedID,$identifier;

					$counter++;
				}
			}
		}
		elsif ($matched) {
			chomp;
			$sequence.=$_;
		}
	}
	if ($matched) { # last entry in DB file matches
		my %countAA;
		foreach my $aa (split(//,uc($sequence))) {
			$countAA{$aa}++;
		}
		my $mass=0;
		foreach my $aa (keys %countAA) {
			$mass+=($massAAave{$aa}*$countAA{$aa}) if $massAAave{$aa}; # some characters are not amino-acids
		}
		$mass+=($massATave{H}+$massATave{H}+$massATave{O}) if $mass; # H(Nter) + OH(Cter)
		$mass=sprintf "%.2f",$mass; # no need for more precision
		foreach my $identifier (@matchedID) {
			$protMW{$identifier}=$mass;
			$protLength{$identifier}=length($sequence);
		}
		###>Writing entry
		my @entryStrg;
		foreach my $identifier (@matchedID) {
			$protDes{$identifier}='no description' unless $protDes{$identifier};
			$protOrg{$identifier}='unknown organism' unless $protOrg{$identifier};
			$protMW{$identifier}=($mass)? $mass : 0;
			$protLength{$identifier}=($sequence)? length($sequence) : 0;
			push @entryStrg,"$identifier ##DES=$protDes{$identifier} ##ORG=$protOrg{$identifier} ##MW=$protMW{$identifier} ##LEN=$protLength{$identifier}";
		}
		print ">",join('',@entryStrg),"\n$sequence\n";
	}
	print "##END";
	close FAS;
	close DB;

	exit;
}

################################
####<Fetching databank file>####
################################
sub getDatabankFile {
	my ($dbName)=@_;
	my $trueDbFile;
	open (MSDAT,"../config/mascot.dat");
	while (<MSDAT>) {
		next unless /^$dbName\t([^\t]+)/;
		$trueDbFile=&getTrueFile($1);
		last;
	}
	close MSDAT;
	unless ($trueDbFile) {
		print "#Error: $dbName not found!";
		exit;
	}
	return $trueDbFile;
}

sub getTrueFile {
	my ($filePathCode)=@_;
	my ($filePath,$fileCode)=($filePathCode=~/^(.+)\/([^\/]+)\Z/);
	#$filePath=~s/.*sequence/$mascotSeqDir/; # only for moved mascot directory
	$fileCode=~s/\*/\.\*/;
	$fileCode=lc($fileCode); # Case mismatches allowed in Windows!!!
	opendir(DIR,$filePath);
	my $fileName;
	while ((my $file=readdir(DIR))) {
		if (lc($file)=~/^$fileCode\Z/) { # Case mismatches allowed in Windows!!!
			$fileName=$file;
			last;
		}
	}
	close DIR;
	if ($fileName) {return "$filePath/$fileName";}
	else {return undef;}
}

sub getParseRules {
	#my $parseRules=uri_unescape($_[0]); # deprotection of url information
	my $parseRules=$_[0];
	#$parseRules=uri_unescape($parseRules) if $parseRules=~/%/;
	$parseRules=~s/�/\+/g; # back comptatibility with myProMS 2.7.2
	my @rules=split(',:,',$parseRules);
#my @rules=split(',:,','ID=([^\|]+_[^\s\|]+),:,ORG=\s?- ([^\(]*),:,');
	my ($idRule)=($rules[0]=~/ID=(.+)/);
	my ($desRule,$orgRule);
	if ($rules[1]) {
		if ($rules[1]=~/DES=/) {($desRule=$rules[1])=~s/DES=//;}
		else {($orgRule=$rules[1])=~s/ORG=//;}
	}
	if ($rules[2]) {($orgRule=$rules[2])=~s/ORG=//;}
	return ($idRule,$desRule,$orgRule);
}

####>Revision history<####
# 1.0.8D Modified for distribution (PP 12/07/16)
# 1.0.8 Detects and moves UniProt isoform key from ACC to ID if ID is selected (PP 08/02/15)
# 1.0.7 Added parse rules test (PP 28/10/13)
# 1.0.6 GPL license (PP 19/09/13)
# 1.0.5 Error control on missing databank (PP 14/09/12)
# 1.0.4 URL encoding for parsing rules (GA 26/07/12)
# 1.0.3 Added ##END at end of fasta (PP 02/11/11)
