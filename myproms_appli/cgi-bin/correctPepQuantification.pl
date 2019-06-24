#!/usr/local/bin/perl -w

#############################################################################
# correctPepQuantification.pl      1.0.0                                    #
# Authors: P. Poullet (Institut Curie)                                      #
# Contact: myproms@curie.fr                                                 #
# Corrects isobaric XIC value based on isotope correction factors)          #
#############################################################################
#----------------------------------CeCILL License-------------------------------
# This file is part of myProMS
#
# Copyright Institut Curie 2019
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

use strict;
use POSIX qw(strftime); # to get the time
use XML::Simple;
use File::Path qw(rmtree);
use File::Copy qw(move);
use promsConfig;
use promsMod;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo('no_user');
my $dbh=&promsConfig::dbConnect('no_user');

####################
####>Parameters<####
####################
my ($pepQuantifID,$jobDir,$quantifItemID,$productID,$duplicateQuantif,$userID)=@ARGV;
($pepQuantifID,$productID)=&promsMod::cleanNumericalParameters($pepQuantifID,$productID);

my ($projectID)=&promsMod::getProjectID($dbh,$pepQuantifID,'quantification');

my $workDir="$promsPath{tmp}/quantification/$jobDir";
my $fileStat= "$workDir/status_$quantifItemID.out";
open(FILESTAT,">$fileStat") || die "Error while opening $fileStat";
print FILESTAT "Started ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
print FILESTAT "1/3 Generating data file\n";
close FILESTAT;

#####>Isotopic coefficients<######
my ($labelType,$distribXMLStrg)=$dbh->selectrow_array("SELECT LABEL_TYPE,ISOTOPIC_DISTRIBUTION FROM ISOTOPIC_CORRECTION WHERE ID_PRODUCT=$productID");
$labelType=~s/:.+//; # TMT:10 -> TMT
	
my $quantifDir="$promsPath{quantification}/project_$projectID/quanti_$pepQuantifID"; # original quanti

my ($quantifAnnot)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$pepQuantifID");
my (@usedTags,@usedTagPos);
foreach my $paramStrg (split('::',$quantifAnnot)) {
	next if $paramStrg !~ /^\d/;
	my ($pos,$tag,$mass)=split(';',$paramStrg);
	push @usedTags,$tag;
	push @usedTagPos,$pos;
}
my @coeffMatrix;
if ($labelType=~/TMT/i) {@coeffMatrix=&convertTMT2matrix($distribXMLStrg,\@usedTags);}
elsif ($labelType=~/ITRAQ/i) {
	#TODO: Implement iTRAQ coefficients
	$dbh->disconnect;
	die "ERROR: Unhandled label type: $labelType!!!";
}
else {
	$dbh->disconnect;
	die "ERROR: Unhandled label type: $labelType!!!";
}

##>Finding XIC parameter ID
my $refQuantifFile="$quantifDir/peptide_quantification_$usedTagPos[0].txt";
my $sthXP=$dbh->prepare("SELECT QP.CODE,QP.ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER QP,QUANTIFICATION_METHOD QM WHERE QP.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND QM.CODE=?");
$sthXP->execute($labelType);
my $xicParamID;
while (my ($code,$paramID)=$sthXP->fetchrow_array) {
	next if $code !~ /_(INTENSITY|AREA)$/;
	#>Head
	my $resp=`head -100 $refQuantifFile | grep -c -P "^$paramID\t"`;
	$resp=~s/\D//g;
	if ($resp && $resp > 0) {
		$xicParamID=$paramID;
		last;
	}
	#>Tail
	$resp=`tail -100 $refQuantifFile | grep -c -P "^$paramID\t"`;
	$resp=~s/\D//g;
	if ($resp && $resp > 0) {
		$xicParamID=$paramID;
		last;
	}
}
$sthXP->finish;

##>Extracting data from quantif files
my %xicData;
my $targetIdx=-1;
foreach my $targetPos (@usedTagPos) {
	$targetIdx++;
	my $pepQuantifFile="$quantifDir/peptide_quantification_$targetPos.txt";
	open(PEP,$pepQuantifFile) || die $!;
	while (<PEP>) {
		next if $.==1;
		next if $_ !~ /^$xicParamID\t/;
		chomp;
		my ($paramID,$pepID,$xicValue)=split(/\t/,$_);
		$xicData{$pepID}[$targetIdx]=$xicValue;
	}
	close PEP;
}
my @xicMatrix;
my $beforeMatrixFile="$workDir/beforeMatrix.txt";
open(MAT,">$beforeMatrixFile");
print MAT "pepID\tTag_",join("\tTag_",@usedTags),"\n"; # headers
foreach my $pepID (sort{$a<=>$b} keys %xicData) {
	print MAT 'p_'.$pepID;
	foreach my $i (0..$#usedTagPos) {
		my $val=$xicData{$pepID}[$i] || 0;
		print MAT "\t",$val;
	}
	print MAT "\n";
}
close MAT;
print " Done.</B>\n";

##>Writing R script
open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
print FILESTAT "2/3 Correcting data\n";
close FILESTAT;
my $afterMatrixFile="$workDir/afterMatrix.txt";
my $Rscript="$workDir/correctXIC.R";
open(R,">$Rscript");
print R qq
|#############################################################################
# Temporary R script for isobaric XIC correction                            #
# Created & called by correctPepQuantifications.cgi                         #
# Deleted at the end of job                                                 #
# Author: P. Poullet (Institut Curie)                                       #
#############################################################################

#>English
Sys.setenv( LANGUAGE = "en")

#>Coefficient matrix
|;
foreach my $i (0..$#coeffMatrix) {
	if ($i==0) {print R 'coeff <- rbind(';}
	else {print R '               ';}
	print R 'c(',join(',',@{$coeffMatrix[$i]}),')';
	print R ',' if $i < $#coeffMatrix;
	print R "\n";
}
print R "              )\n\n";

print R qq
|#>Loading XIC matrix
xicBefore <- read.table("$beforeMatrixFile",header=TRUE,row.names="pepID")

#>Applying correction coefficients
xicAfter <- t(solve(coeff,t(xicBefore)))

#>Writing corrected XIC matrix to file
write.table(xicAfter,file="$afterMatrixFile",quote=FALSE,sep="\\t",col.names=FALSE)
|;
close R;

##>Running R script
system "cd $workDir; $promsPath{R}/R CMD BATCH --no-save --no-restore $Rscript";

##>ERROR Management
my $RoutFile=$Rscript.'out';
my $RoutStrg=`tail -3 $RoutFile`;
unless ($RoutStrg=~/proc\.time\(\)/) {
	$dbh->disconnect;
	$RoutStrg=`tail -20 $RoutFile`;
	my $RerrorStrg="R script has generated an error!";
	my $inError=0;
	foreach my $line (split(/\n/,$RoutStrg)) {
		next if (!$inError && $line !~ /^Error in/); # skip everything before "Error in..."
		$inError=1;
		$RerrorStrg.="\n$line";
	}
	die $RerrorStrg;
}

##>Parsing results in afterMatrix.txt
open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
print FILESTAT "3/3 Importing corrected data\n";
close FILESTAT;

my %correctedData;
open(CORR,$afterMatrixFile) || die $!;
while(<CORR>) {
	chomp;
	my ($pepID,@correctValues)=split(/\t/,$_);
	$pepID=~s/\D+//; # remove "p_" prefix
	foreach my $i (0..$#usedTagPos) {
		$correctedData{$pepID}[$i]=($correctValues[$i] && $correctValues[$i] > 0)? $correctValues[$i] : 0;
	}
}
close CORR;

##>Writing final peptide quantif files
$targetIdx=-1;
foreach my $targetPos (@usedTagPos) {
	$targetIdx++;
	my $inQuantifFile="$quantifDir/peptide_quantification_$targetPos.txt";
	my $outQuantifFile="$workDir/peptide_quantification_$targetPos.txt";
	open(IN,$inQuantifFile) || die $!;
	open(OUT,">$outQuantifFile");
	while (<IN>) {
		if ($_ !~ /^$xicParamID\t/) {
			print OUT $_;
			next;
		}
		my ($paramID,$pepID,$xicValue)=split(/\t/,$_);
		print OUT "$paramID\t$pepID\t$correctedData{$pepID}[$targetIdx]\n";
	}
	close IN;
	close OUT;
}

####>Updating database<####
$quantifAnnot.="::CORRECTION=$duplicateQuantif;";
foreach my $i (0..$#coeffMatrix) {
	$quantifAnnot.='&' if $i > 0;
	$quantifAnnot.=join(',',@{$coeffMatrix[$i]});
}
my $usedQuantifDir;
if ($duplicateQuantif) {
	my %sthCorr=(SEL_Q=>$dbh->prepare("SELECT ID_QUANTIFICATION_METHOD,NAME FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?"),
			  INS_Q=>$dbh->prepare("INSERT INTO QUANTIFICATION (ID_QUANTIFICATION_METHOD,ID_PRODUCT,NAME,FOCUS,QUANTIF_ANNOT,STATUS,UPDATE_DATE,UPDATE_USER) VALUES (?,?,?,'peptide',?,1,NOW(),?)"),
			  INS_PAR=>$dbh->prepare("INSERT INTO PARENT_QUANTIFICATION (ID_PARENT_QUANTIFICATION,ID_QUANTIFICATION,FUNCTION) VALUES (?,?,'CORRECTION')"),
			  SEL_AQ=>$dbh->prepare("SELECT ID_ANALYSIS,QUANTIF_FILE,IS_REFERENCE FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=?"),
			  INS_AQ=>$dbh->prepare("INSERT INTO ANA_QUANTIFICATION (ID_QUANTIFICATION,ID_ANALYSIS,QUANTIF_FILE,IS_REFERENCE) VALUES (?,?,?,?)")
	);
	$sthCorr{SEL_Q}->execute($pepQuantifID);
	my ($quantifMethID,$name)=$sthCorr{SEL_Q}->fetchrow_array;
	my $newName=$name.' [Corrected '.(strftime("%Y/%m/%d %H:%M:%S",localtime)).']';
	$sthCorr{INS_Q}->execute($quantifMethID,$productID,$newName,$quantifAnnot,$userID);
	my $corrQuantifID=$dbh->last_insert_id(undef, undef,'QUANTIFICATION','ID_QUANTIFICATION');
	#$sthCorr{INS_PAR}->execute($pepQuantifID,$corrQuantifID);
	$sthCorr{SEL_AQ}->execute($pepQuantifID);
	while (my ($anaID,$qFile,$isRef)=$sthCorr{SEL_AQ}->fetchrow_array) {
		$sthCorr{INS_AQ}->execute($corrQuantifID,$anaID,$qFile,$isRef);
	}
	foreach my $sth (values %sthCorr) {
		$sth->finish;
	}
	
	$usedQuantifDir="$promsPath{quantification}/project_$projectID/quanti_$corrQuantifID";
	mkdir $usedQuantifDir || die $!;
}
else { # update source quantif (should be allowed only if quantif has never been used for prot quantif!)
	foreach my $targetPos (@usedTagPos) { # backuping original files
		my $srcQuantifFile="$quantifDir/peptide_quantification_$targetPos.txt";
		my $destQuantifFile="$quantifDir/backup_quantification_$targetPos.txt";
		move($srcQuantifFile,$destQuantifFile) || die $!;
	}
	my $sthUpCorr=$dbh->prepare("UPDATE QUANTIFICATION SET NAME=CONCAT(NAME,' [Corrected]'),ID_PRODUCT=?,QUANTIF_ANNOT=?,UPDATE_DATE=NOW(),UPDATE_USER=? WHERE ID_QUANTIFICATION=?");
	$sthUpCorr->execute($productID,$quantifAnnot,$userID,$pepQuantifID);
	$sthUpCorr->finish;
	
	$usedQuantifDir=$quantifDir;
}

##>Move/clean files & directory
foreach my $targetPos (@usedTagPos) {
	my $srcQuantifFile="$workDir/peptide_quantification_$targetPos.txt";
	my $destQuantifFile="$usedQuantifDir/peptide_quantification_$targetPos.txt";
	move($srcQuantifFile,$destQuantifFile) || die $!;
}
rmtree($workDir);
unlink $workDir if -e $workDir;

$dbh->commit;
$dbh->disconnect;

sub convertTMT2matrix {
	my ($distribXMLStrg,$refUsedTags)=@_;
	my %tagIndex=(C=>-1,N=>-1);
	my @header=('Tag',-2..2);
	my $xml = new XML::Simple();
	my $xmlData = $xml->XMLin($distribXMLStrg);
	
	my (@distribMatrix,%tags,@tagList,%tag2class);
	foreach my $tagData (@{$xmlData->{MASS_TAG}}) {
		push @distribMatrix,[$tagData->{NAME},$tagData->{MINUS_2},$tagData->{MINUS_1},100,$tagData->{PLUS_1},$tagData->{PLUS_2}];
		my $tag=$tagData->{NAME};
		my $tagClass=($tag eq 'TMT-126' || $tag=~/C$/)? 'C' : 'N'; # TMT-131 -> N
		push @{$tags{$tagClass}},$tag;
		$tag2class{$tag}=$tagClass;
		push @tagList,$tag;
	}
#TMT-Tag	-2	-1	0	+1	+2
#TMT-126	0	0	100	8.2	0.4
#TMT-127N	0	0.4	100	6.3	0
#TMT-127C	0	0.5	100	6.3	0
#TMT-128N	0	0.9	100	6.8	0
#TMT-128C	0.2	1.6	100	6.3	0
#TMT-129N	0	2.5	100	5	0
#TMT-129C	0	2.3	100	4.3	0
#TMT-130N	0	2.7	100	3.9	0
#TMT-130C	0	1.7	100	2.7	0
#TMT-131	0	3.9	100	3.7	0

	#my $l=-1;
	#foreach my $line (split(/\n/,$distribXMLStrg)) {
	#	$l++;
	#	my @lineData=split(/\t/,$line);
	#	if ($l==0) {
	#		@header=@lineData;
	#		next;
	#	}
	#	push @distribMatrix,\@lineData;
	#	my $tag=$lineData[0]; # remove tag column;
	#	my $tagClass=($tag eq '126' || $tag=~/C$/)? 'C' : 'N';
	#	push @{$tags{$tagClass}},$tag;
	#	$tag2class{$tag}=$tagClass;
	#	push @tagList,$tag;
	#}
	
	my %matrix;
	foreach my $i (0..$#distribMatrix) {
		my @lineData=@{$distribMatrix[$i]};
		my $tagClass=$tag2class{$lineData[0]};
		$tagIndex{$tagClass}++;	
		foreach my $j (1..$#lineData) {
			my $targetIndex=$tagIndex{$tagClass}+$header[$j];
			next if ($targetIndex < 0 || $targetIndex > $#{$tags{$tagClass}}); # outside tag range
			$matrix{$tagClass}[$targetIndex][$i]=$lineData[$j]/100;
		}	
	}
	##<Filling undef cells
	foreach my $tagClass (keys %matrix) {
		foreach my $i (0..$tagIndex{$tagClass}) {
			foreach my $j (0..$#distribMatrix) {
				$matrix{$tagClass}[$i][$j]=0 if (!$matrix{$tagClass}[$i] || !$matrix{$tagClass}[$i][$j]);
			}
		}
	}
	##<Merging sub-matrices C & N
	my @fullMatrix;
	foreach my $i (0..$tagIndex{C}) {
		foreach my $tagClass ('C','N') {
			push @fullMatrix,$matrix{$tagClass}[$i];
		}
	}
	##<Restricting matrix to used tags
	my @finalMatrix;
	if (!$refUsedTags || scalar @{$refUsedTags} == scalar @fullMatrix) {
		@finalMatrix=@fullMatrix;
	}
	else {
		my %usedTags;
		foreach my $tag (@{$refUsedTags}) {$usedTags{$tag}=1;}
		foreach my $i (0..$#fullMatrix) {
			next unless $usedTags{$tagList[$i]};
			my @coeffRow;
			foreach my $j (0..$#{$fullMatrix[$i]}) {
				push @coeffRow,$fullMatrix[$i][$j] if $usedTags{$tagList[$j]};
			}
			push @finalMatrix,\@coeffRow;
		}
	}

	return @finalMatrix;
}

# TODO: Implement iTRAQ coefficients
####>Revision history<####
# 1.0.0 New script for isobaric XIC correction (PP 19/02/19)