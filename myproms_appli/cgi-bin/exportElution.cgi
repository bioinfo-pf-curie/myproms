#!/usr/local/bin/perl -w

################################################################################
# exportElution.cgi    2.0.3                                                   #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Exports list of query with masses and elution times                          #
# for validated or non-validated MS/MS analysis.                               #
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
#############################################################

$| = 1;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use strict;

#print header; warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;

####################
####>Parameters<####
####################
my $analysisID=param('ana_ID');
my $projectID=param('proj_ID');
# my %convCode=('%20'=>' ','%28'=>'(','%29'=>')','%2c'=>',','%2d'=>'-','%2e'=>'.','%3a'=>':');

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

my ($analysisName,$validStatus,$dataFile,$wiffFile, $fileFormat)=$dbh->selectrow_array("SELECT NAME,VALID_STATUS,DATA_FILE,WIFF_FILE,FILE_FORMAT FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
$wiffFile='' unless $wiffFile;


############################################
####>Scanning database for elution data<####
############################################
my (%massObs,%charge,%elutionTime);
my $dataType;

if ($validStatus<=0) {
	my $noValidSth=$dbh->prepare("SELECT QUERY_NUM,ELUTION_TIME,MASS_DATA FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID AND QUERY_NUM>0");
	$noValidSth->execute;
	while (my ($queryNumber,$eluTime,$massData)=$noValidSth->fetchrow_array) {
		my ($ObsMass)=($massData=~/OBS=(\d+\.?\d*)/);
		my ($ExpMass)=($massData=~/EXP=(\d+\.?\d*)/);
		$massObs{$queryNumber}=sprintf "%.3f",$ObsMass;
		my $nonRoundedCharge=$ExpMass/$ObsMass +0.4;
		$charge{$queryNumber}=sprintf "%.0f",$nonRoundedCharge;
		($dataType,@{$elutionTime{$queryNumber}})=&getElutionTime($eluTime) if $eluTime;
	}
	$noValidSth->finish;
}
else {
	my $validSth=$dbh->prepare("SELECT QUERY_NUM,ELUTION_TIME,MR_OBS,MR_EXP FROM PEPTIDE WHERE ID_ANALYSIS=$analysisID");
	$validSth->execute;
	while (my ($queryNumber,$eluTime,$ObsMass,$ExpMass)=$validSth->fetchrow_array) {
		$massObs{$queryNumber}=sprintf "%.3f",$ObsMass;
		my $nonRoundedCharge=$ExpMass/$ObsMass +0.4 ;
		$charge{$queryNumber}=sprintf "%.0f",$nonRoundedCharge;
		($dataType,@{$elutionTime{$queryNumber}})=&getElutionTime($eluTime) if $eluTime;
	}
	$validSth->finish;
}
# my $currentQuery;
# open(ELU,$elutionFile) || die "Cannot find elution data file.";
# while (<ELU>) {
	#>Mass & charge
	# if (/^qexp(\d+)=(\S+),(\S+)/) { # Mass Obs + charge
		# $massObs{$1}=$2;
		# $charge{$1}=$3;
		# next;
	# }
	# if (/name="query(\d+)"/) {$currentQuery=$1; next;}
	# next unless $currentQuery;
	# last if /name="index"/;
# 	if (/^title=(\S+)/) {
# 		my $title=$1;
# 		foreach my $code (keys %convCode) {
# 			$title=~s/$code/$convCode{$code}/g;
# 		}
# 		if ($title=~/Elution: (\S+) to (\S+) min/) {@{$elutionTime{$currentQuery}}=($1,$2);}
# 		elsif ($title=~/Elution: (\S+) min/) {@{$elutionTime{$currentQuery}}=($1,'-');}
# 		else {@{$elutionTime{$currentQuery}}=(-1,-1);}
# 	}
	# if (/^title=.*Elution.*%3a%20(.+)[Pp]eriod/) {
		# my $elutStrg=$1;
		# while ($elutStrg=~/%(\w{2})/) {
			# my $code=$1;
			# my $char=chr(hex($code));
			# $elutStrg=~s/%$code/$char/g;
		# }
		# if ($elutStrg=~/^(\S+) to (\S+)/) {@{$elutionTime{$currentQuery}}=($1,$2);}
		# elsif ($elutStrg=~/^(\S+) min/) {@{$elutionTime{$currentQuery}}=($1,'-');}
		# else {@{$elutionTime{$currentQuery}}=(-1,-1);}
	# }
# }
# close ELU;
$dbh->disconnect;


########################
####>Starting EXCEL<####
########################
print header(-type=>"application/vnd.ms-excel",-attachment=>"$analysisName.xls");

################################
####>Printing data in table<####
################################
my ($col1Title,$col2Title)=($dataType=~/Spectrum/)? ('Spectrum','scan') : ('From','To');
print qq
|<TABLE border=1>
<TR><TH colspan=5><FONT style="font-size:18px">$analysisName</FONT><BR>($dataFile - $wiffFile)</TH></TR>
<TR><TH colspan=2 nowrap>$dataType</TH><TH rowspan=2>Mass (obs)</TH><TH rowspan=2>Charge</TH><TH rowspan=2>Query</TH></TR>
<TR><TH>$col1Title</TH><TH>$col2Title</TH></TR>
|;

foreach my $query (sort{$elutionTime{$a}[0]<=>$elutionTime{$b}[0] || $charge{$a} cmp $charge{$b}} keys %elutionTime) {
	print "<TR><TD>$elutionTime{$query}[0]</TD><TD>$elutionTime{$query}[1]</TD><TD>$massObs{$query}</TD><TD align=center>$charge{$query}</TD><TD>$query</TD></TR>\n";
}
print "</TABLE>\n";


####################
####<Subroutine>####
####################
sub getElutionTime {
	my ($eluTime)=@_;
	my $dataType;
	my @elutionTime;
	if ($eluTime=~/et(\d+\.\d+);/) { # new elution time format
		$elutionTime[0]=$1;
		if ($eluTime=~/\sto\s(\S+)/) {$elutionTime[1]=$1;}
		else {$elutionTime[1]='-';}
		$dataType='Elution time (min)';
	}
	elsif ($eluTime=~/^sp(\d+) sc(\d+)/) { # spectrum scan format
		@elutionTime=($1,$2);
		$dataType='Spectrum & Scan';
	}
	elsif ($eluTime=~/^rt(\S+)/) { # new retention time format
		@elutionTime=($1,'-');
		$dataType='Retention time (min)';
	}
	else { # old elution time format
		if ($eluTime=~/\sto\s/) {
			@elutionTime=split(/\sto\s/,$eluTime);
		}
		else {
			$elutionTime[0]=$eluTime;
			$elutionTime[1]='-';
		}
		$dataType='Elution time (min)';
	}
	return ($dataType,@elutionTime);
}

####>Revision history<####
# 2.0.3 GPL license (PP 23/09/13)
# 2.0.2 Minor update to export well the elution time (garras 06/01/11)<BR>See 2.6.5 modification in storeAnalyses.cgi
