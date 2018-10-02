#!/usr/local/bin/perl -w

################################################################################
# msms_gif.cgi                       2.1.8                                     #
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
use strict ;
use CGI ':standard';
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use POSIX qw(strftime); # to get the time
#use promsConfig;
#use promsMod;
use GD;

my $file = param('file');
my $spectrumType= param('TYPE');

#print header; warningsToBrowser(1);

#######################
####>recovering data<####
#######################
my %expIonTable ;
my %knowIon ;
my $IonList ;
open (DESC, "$file") ;
while (<DESC>) {
	my $line = $_ ;
	if ($line =~ /ionserie=(.+)/) {
		$IonList=$1 ;
		chomp ($IonList) ;
		my ($mz, $intens, $desc)  = split (/:/ , $IonList) ;
		$expIonTable{$mz}[0] = $intens ;
		$expIonTable{$mz}[1] = $desc if defined $desc ;
	}
}
close (DESC) ;
unlink ("$file") ; # remove tempory file ( .pgf format ! )


#######################
####>String Setup####
#######################
foreach my $mz (keys (%expIonTable)) {
	if (defined ($expIonTable{$mz}[1])) { #interpretation of fragment
			next unless $expIonTable{$mz}[1] ne " ";
			$expIonTable{$mz}[1] =~ s/-OxiM(\d+)/$1-64/g; #replacing neutral loss Ox M by -64
			$expIonTable{$mz}[1] =~ s/-PHOS(\d+)/$1-98/g; #replacing neutral loss Phospho  by -96
			$expIonTable{$mz}[1] =~ s/(,?)([^,]+?)\+\+/$1\($2\)\+\+/g; #setting (  ) if ++
			#$expIonTable{$mz}[1] =~ s/(.+),(y\d.*)/$2,$1/; #place y in first position
			#$expIonTable{$mz}[1] =~ s/\s,/ /;

			$knowIon{$mz}[0]=$expIonTable{$mz}[0]; #intensity


	}
}


#######################
####>creating spectrum viex<####
#######################

####>calculating<####
my $minMZ ;
my $maxMZ = 0 ;
my $higtIntensity =0;
foreach my $mz (keys (%expIonTable)) {
	$minMZ = $mz if !defined($minMZ) ;
	$minMZ = $mz if $mz<$minMZ ;
	$maxMZ = $mz if $mz>$maxMZ ;
	$higtIntensity = $expIonTable{$mz}[0] if $expIonTable{$mz}[0]>$higtIntensity ;
}


####>creating new Image<####
my $imgHight = param ('height');
my $imgWidht = param ('width');
my $leftMarging = 10 ;
my $rightMarging = 10;
my $topMarging = 40 ;
my $bottonMarging = 1 ;
my $bottomSpace = 20 ;

my $img = new GD::Image($imgWidht,$imgHight);

my $white = $img->colorAllocate(255,255,255);
my $black = $img->colorAllocate(0,0,0);
my $red = $img->colorAllocate(255,0,0);
my $green = $img->colorAllocate(0,255,0);
my $blue = $img->colorAllocate(0,0,255);
my $fushia = $img->colorAllocate(255,150,255);
my $lightblue = $img->colorAllocate(0,255,255);
$img->interlaced('true');
$img->transparent($white);

my $XRange =$maxMZ-$minMZ ;
my ($intervall)= ($XRange>800)?100:($XRange>400)?50:($XRange>80)?10:($XRange>16)?2:1 ;
my $Ystar = $imgHight -($bottonMarging+$bottomSpace) ;
my $Xstart = $minMZ-$intervall/5 ;
my $Xscale = 0.8*$imgWidht/$XRange ;
my $Yscale = 0.96*(($imgHight-($topMarging+$bottonMarging))/($higtIntensity));

my $lineColor = ($spectrumType eq 'ref')?$blue:$red;
my $fontColor = ($spectrumType eq 'ref')?$black:$blue;
my $dashedLineColor = ($spectrumType eq 'ref')?$fushia:$lightblue;

my $mzEcart = 10/$Xscale ;
my $YGap = 5/$Yscale ;

foreach my $mz ( keys (%knowIon)) { #sort pour debug
	my $hight =  $knowIon{$mz}[0];
	foreach my $localMz ( keys (%expIonTable)) {  #sort pour debug
		if (($localMz<$mz+$mzEcart) && ($localMz>$mz-$mzEcart) && ($expIonTable{$localMz}[0]>$hight)) {
			$knowIon{$mz}[2]= $expIonTable{$localMz}[0];
			if (defined ($expIonTable{$localMz}[1]) ){
				$knowIon{$mz}[2] += ($expIonTable{$localMz}[1] ne " ")?($YGap * (length ($expIonTable{$localMz}[1])+1)):$YGap;
				$knowIon{$mz}[2] += $YGap if $expIonTable{$localMz}[1] =~ /\+/ ;
			}
			$hight = $knowIon{$mz}[2];
		}
	}
 }
foreach my $mz (  keys (%knowIon)) { #sort pour debug
	next unless ($knowIon{$mz}[2]);
	my $hight = $knowIon{$mz}[2];
	my $exitValue ;
	do	{
		$exitValue = 0 ;
		foreach my $localMz ( keys (%knowIon)) {  #sort pour debug
			next unless ($localMz!=$mz) && ($localMz<$mz+$mzEcart) && ($localMz>$mz-$mzEcart) && defined($knowIon{$localMz}[2]);
			if (($knowIon{$mz}[2] < $knowIon{$localMz}[2] + $YGap * (length ($expIonTable{$localMz}[1])+1)) && ($knowIon{$localMz}[2] < $knowIon{$mz}[2] + $YGap * (length ($expIonTable{$mz}[1])+1))) {
				$knowIon{$mz}[2] = $knowIon{$localMz}[2] + $YGap * (length ($expIonTable{$localMz}[1])+1) ; #conflit resolution
				$knowIon{$mz}[2] += 2 * $YGap if $expIonTable{$localMz}[1] =~ /\+/ ;
				$exitValue = 1;
			}
		}
	} while ($exitValue == 1 ) ;
}

foreach my $mz (keys (%expIonTable)) {
	$higtIntensity = $knowIon{$mz}[2] if defined ($knowIon{$mz}[2]) && $knowIon{$mz}[2]>$higtIntensity ;
}
$Yscale = 0.96*(($imgHight-($topMarging+$bottonMarging))/($higtIntensity));

####>creating spectrum line<####
foreach my $mz ( keys (%expIonTable)) {
	my $Xpos = $leftMarging + ($mz-$Xstart)* $Xscale ;
	my $Yvalue = $Ystar -($expIonTable{$mz}[0]*$Yscale);
	$img->line($Xpos,$Ystar,$Xpos,$Yvalue,$lineColor);
	if (defined $expIonTable{$mz}[1]) {
		if ($expIonTable{$mz}[1] eq " ") {
			$img->stringUp(gdSmallFont,$Xpos-7,$Yvalue-2,"-",$fontColor);
		}
		elsif ($knowIon{$mz}[2]) {
			my $Yvaluebis = $Ystar -($knowIon{$mz}[2]*$Yscale);
			$img->stringUp(gdSmallFont,$Xpos-7,$Yvaluebis-4,$expIonTable{$mz}[1],$fontColor);
			$img->dashedLine($Xpos,$Yvalue-2,$Xpos,$Yvaluebis-2,$dashedLineColor);
		}
		else {
			$img->stringUp(gdSmallFont,$Xpos-7,$Yvalue-2,$expIonTable{$mz}[1],$fontColor);
		}
	}
}

####>creating axe<####
my $i= 0 ;
my $j=0 ;

while ($i<=($maxMZ+$intervall)) {
	$j++; 	last if $j>1000 ;
	if ($i>=($minMZ-$intervall)) {
		$img->line($leftMarging+($i-$Xstart)*$Xscale,$Ystar,$leftMarging+($i-$Xstart)*$Xscale,$Ystar+$bottomSpace/3,$black);
		$img->string(gdSmallFont,$leftMarging-8+($i-$Xstart)*$Xscale,$Ystar+$bottomSpace/3,"$i",$black);
	}
	$i += $intervall;
}
#$img->line($leftMarging,$Ystar,$imgWidht-$rightMarging,$Ystar,$black);
$img->line($leftMarging,$Ystar,$leftMarging+($i-($Xstart+$intervall))*$Xscale,$Ystar,$black);
####>printing Image<####
binmode STDOUT;
print "Content-type: image/png\n\n";
print $img->png;

####>Revision history<####
# 2.1.8  GPL license (PP 23/09/13)
