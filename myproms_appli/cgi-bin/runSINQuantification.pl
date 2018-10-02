#!/usr/local/bin/perl -w

################################################################################
# runSINQuantification.pl           1.2.0                                      #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Allows quantification of proteins by computing the Spectral Index            #
# Normalized as described by in the following reference :                      #
# Label-free, normalized quantification of complex mass spectrometry data      #
# for proteomic analysis. N.M. Griffin, J. Yu, F. Long, S. Shore, Y. Li,       #
# J. A Koziol, J. E Schnitzer. Nature Biotechnology 28 (2010) - 83-90.         #
# To get the intensity of validated spectra, it uses the peptide_view.cgi      #
# method adapted to the purpose of the analysis                                #
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
use promsConfig;
use promsMod;
use strict;
use POSIX qw(strftime);

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo('no_user');
my $minPepPerProt=3;# validate proteins
my $matchInterval=0.8;# user preference
my $intPercent=0.01;# user preference TODO --> update this criteria!!!y $useSameSet=1;#if 0, do not consider proteins in same group match

###############################
####>Recovering parameters<####
###############################
my ($quantiID,$quantifDate,$quantItemID,$userID)=@ARGV;
my ($analysisID,$qID)=split(/\./,$quantItemID);
my $quantifDir="$promsPath{tmp}/quantification/$quantifDate";
my $fileStat="$quantifDir/status_$quantItemID.out";
open(FILESTAT,">$fileStat") || die ("Error while opening $fileStat") ;
#print FILESTAT "#@ARGV#\n";
#print strftime("%H:%M:%S %d/%m/%Y",localtime);
print FILESTAT "Started ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";

###############################################################
###>Retrieve information from myProMS and compute SI values<###
###############################################################
my (%myQueries,%validatedPep,%FixMods,%massValueMods,%massValueFixedMods,%expIonTable,%filteredIon);
my $dbh=&promsConfig::dbConnect('no_user');
my ($dataFile,$minScore,$maxRank,$fileType)=$dbh->selectrow_array("SELECT DATA_FILE, MIN_SCORE, MAX_RANK, FILE_FORMAT FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
if($fileType !~ /MASCOT/) {
	open(FILESTAT,">>$fileStat");
	print FILESTAT "The analysis does not seem to habe been searched with Mascot software\n";
	print FILESTAT "Ended ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
	close(FILESTAT);
	exit;
}
my $DatFile="$promsPath{valid}/ana_$analysisID/$dataFile";
#print "$DatFile<BR>\n";
############################################
###>1st: Get all true validated peptides<###
############################################
#my $validPep=$dbh->prepare("SELECT ID_PEPTIDE,QUERY_NUM,PEP_SEQ,PEP_RANK FROM PEPTIDE WHERE ID_ANALYSIS=$analysisID AND PEP_RANK IS NOT NULL");
my $validPep=$dbh->prepare("SELECT ID_PEPTIDE,QUERY_NUM,PEP_SEQ,PEP_RANK FROM PEPTIDE WHERE ID_ANALYSIS=$analysisID AND PEP_RANK IS NOT NULL AND VALID_STATUS=1"); # Change on 2017/02/02 so as to not get ghost-peptides with no fragmentation spectrum
$validPep->execute;
print FILESTAT "1/5 Getting validated peptides from myProMS database\n";
close FILESTAT;
while( my ($peptideID,$queryNum,$pepSeq,$pepRank) = $validPep->fetchrow_array) {
	my $varMod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$peptideID,$analysisID,$pepSeq);
	$varMod='' unless $varMod;
	$validatedPep{"$pepSeq $varMod"}{'peptideID'}=$peptideID;
	$validatedPep{"$pepSeq $varMod"}{'query'}{$queryNum}=$pepRank;
	$myQueries{$analysisID}{$queryNum}=1;
}
$validPep->finish;

###--- Old version ---### -> It was not working for analysis that were shut
#
##########################################################################
#####>2nd: Get all the queries that match the peptide sequence peptide<###
##########################################################################
#my $pepString="INFO_PEP1";
#for(my $r=2;$r<=$maxRank;$r++){
#	$pepString.=",INFO_PEP$r";
#}
#my $sthallQueries=$dbh->prepare("SELECT QUERY_NUM,$pepString FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID AND QUERY_NUM>0");
#$sthallQueries->execute;
#
#
#open(FILESTAT,">>$fileStat");
#print FILESTAT "2/5 Getting all queries from myProMS database\n";
#close FILESTAT;
#while (my ($queryNum,@pepInfo)=$sthallQueries->fetchrow_array) {
#	for (my $ri=0;$ri<=$#pepInfo;$ri++) {
#		my $rankInfo=$pepInfo[$ri];
#		last unless $rankInfo;
#		my ($varMod)=($rankInfo=~/VMOD=([^,]+)/);
#		$varMod='' unless $varMod;
#		my ($pepSeq)=($rankInfo=~/SEQ=(\w+)/);
#		if($validatedPep{"$pepSeq $varMod"} && $rankInfo=~/SEL=[12]/){# 1: auto-selected / 2: manually selected
#			$validatedPep{"$pepSeq $varMod"}{'query'}{$queryNum}=$ri+1;
#			$myQueries{$analysisID}{$queryNum}=1;
#			#print "query=$queryNum"," rank=",$validatedPep{"$pepSeq $varMod"}{'query'}{$queryNum},"<BR>\n";
#			#print "$pepSeq $varMod","<BR>\n";
#		}
#	}
#}
#$sthallQueries->finish;
#$dbh->disconnect;
###--- Old version ---###

############################################################################################
###>2nd -> New (15/07/2011): Get all the queries that match the peptide sequence peptide<###
############################################################################################
###> Read the 1st part of the Mascot file (PDM or not) and get the queries
###> by reading the peptide part of
if (!-e $DatFile){
	my ($projectID)=&promsMod::getProjectID($dbh,$analysisID,'analysis');
	my $pepFileExt=($fileType=~/\.PDM/)? 'pdm' : 'dat';
	my $oldPepFileName=sprintf "P%06d.$pepFileExt",$analysisID;
	$DatFile="$promsPath{peptide}/proj_$projectID/$oldPepFileName";
}
$dbh->disconnect;

open (DATAFILE, $DatFile) || (print FILESTAT "Unable to open $DatFile file\n");
my ($section,$dateCode)=('','');
my (%varMods,%codeToquery,%queryToCode,%queryToSeq,%queryToVmod,%massQueryExp,%massQueryObs);
while (my $line=<DATAFILE>) {
	$section=$1 if $line=~/name="(\w+)"/;
	last if ($section eq 'proteins');
	next if ($section && $section ne 'peptides' && $section ne 'summary' && $section ne 'parameters');
	$line=~s/\s*\Z//; # delete \n and trailing spaces
	if ($section eq 'peptides'){
		if ($line=~/^q(\d+)_p(\d+)=/) {
			my ($queryNum,$rank)=($1,$2);
			my @queryInfo=split(/,/,$line);
			next unless $queryInfo[6];
			my $seqCode=$queryInfo[6];
			$codeToquery{$seqCode}{"$queryNum:$rank"}=1;
			$queryToCode{"$queryNum:$rank"}=$seqCode;
			$queryToSeq{"$queryNum:$rank"}=$queryInfo[4];
			my $varModString='';
			# my $varPosModString ;
			foreach my $numMod (sort{$a<=>$b} keys %varMods) {
				my $position='';
				my $currentPosition=0;
				if ($queryInfo[7]=~/$numMod/) {
					my $currentVarMod=$varMods{$numMod};
					my @pos;
					while ($currentPosition != -1) {
						$currentPosition=index ($queryInfo[7],$numMod,$currentPosition+1);
						push @pos,$currentPosition if $currentPosition >= 0 ;
					}
					if (scalar @pos) {
						my $posString=':'.join('.', @pos);
						$currentVarMod=~s/\(([^\(]+)\)\Z/\($1$posString\)/;
					}
					$varModString.=" + $currentVarMod";
				}
			}
			$queryToVmod{"$queryNum:$rank"}=$varModString;
		}
	}
	elsif ($section eq 'summary') {
		if ($line=~/^qmass(\d+)=(\S+)/) {$massQueryExp{$1}=$2;} # Mr(exp)
		elsif ($line=~/^qexp(\d+)=(\S+),(\d+)\+/) {
			$massQueryObs{$1}=$2; # Observed
		}
	}
	elsif ($section eq 'parameters') {
		if ($line=~/^IT_MODS=(.+)\n/) {
			(my $modString=$1)=~s/\s*\Z//; # trailing \r !!!?
			my $numMod=1;
			foreach my $mod (split(/,/,$modString)) {
				$mod=~s/^_+//; # __Oxidation -> Oxidation
				$varMods{$numMod}=$mod;
				$numMod++;
			}
		}
	}
}
close DATAFILE;

foreach my $valipPep (keys(%validatedPep)) {
	foreach my $queryNum (keys %{$validatedPep{$valipPep}{'query'}}){
		my $rank=$validatedPep{$valipPep}{'query'}{$queryNum};
		my $seqCode=$queryToCode{"$queryNum:$rank"};
		next unless $seqCode;# case query was like that: q1_p1=-1, no seqCode!
		foreach my $newQuery (keys %{$codeToquery{$seqCode}}) {
			my ($newQueryNum,$newRank)=split(/:/,$newQuery);
			if($queryToSeq{"$newQueryNum:$newRank"} eq $queryToSeq{"$queryNum:$rank"}){
				$validatedPep{$valipPep}{'query'}{$newQueryNum}=$newRank;
				$myQueries{$analysisID}{$newQueryNum}=1;
			}
		}
	}
}

################################################################################
###>3rd: Read the file and keep the ions of the queries that are interesting<###
################################################################################
open(FILESTAT,">>$fileStat");
print FILESTAT "3/5 Fetching ion information (m/z and intensity)\n";
close FILESTAT;
&copyIonsInformation($DatFile,$analysisID,\%FixMods,\%massValueMods,\%massValueFixedMods,\%expIonTable,\%filteredIon,\%myQueries);
######################################################
###>4th: Compute the SI for all top-match proteins<###
######################################################
open(FILESTAT,">>$fileStat");
print FILESTAT "4/5 Computing SI for peptides\n";
close(FILESTAT);
my $parentQuanti=&computeSIPeps($analysisID,\%massValueMods,\%massValueFixedMods,\%FixMods,\%expIonTable,\%filteredIon,\%validatedPep);
open(FILESTAT,">>$fileStat");
print FILESTAT "5/5 Computing SIN for proteins\n";
close(FILESTAT);
&computeSIProts($minPepPerProt,$analysisID,\%validatedPep,$parentQuanti);
open(FILESTAT,">>$fileStat");
print FILESTAT "Ended ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
close(FILESTAT);

$dbh=&promsConfig::dbConnect('no_user');
$dbh->do("UPDATE QUANTIFICATION SET STATUS=1 WHERE ID_QUANTIFICATION=$quantiID") || die $dbh->errstr();
$dbh->commit;
$dbh->disconnect;

sleep 2;
unlink $fileStat;

#----- Functions -----#

#################################################################
####>Function that copy into hash tables the ion information<####
#################################################################
sub copyIonsInformation{
	my ($DatFile,$analysisID,$refFixMods,$refmassValueMods,$refmassValueFixedMods,$refexpIonTable,$reffilteredIon,$refmyQueries)=@_;
	#print "datFile=$DatFile\n";
	open (FILE, $DatFile) || die "Unable to open $DatFile";
	my ($localIonSerie,$localQuery,$localIntMin,$localIntMax);
	my $onQuery=0;
	my $onMasses=0;

	while (my $line=<FILE>) {
		if ($line =~ /name=\"masses\"/) {
			$onMasses=1;
			####>creating %massValueFixedMods hash<####
			foreach my $numMod (keys (%{$refFixMods->{$analysisID}})) {
				$refmassValueFixedMods->{$analysisID}{$refFixMods->{$analysisID}{$numMod}[2]} = $refFixMods->{$analysisID}{$numMod}[1] ;
			}
		}
		elsif ($line =~ /name=\"query([0-9]+)\"/ && $refmyQueries->{$analysisID}{$1} ) {
			$localQuery=$1;
			$onQuery=1;
		}
		elsif ($line =~ /gc0p4Jq0M2Yt08jU534c0p/) {
			$onQuery=0;
			$onMasses=0;
		}
		elsif ($onMasses==1 && $line =~ /FixedMod(\d+)=(.+)\,(.+)/) {
			$refFixMods->{$analysisID}{$1}[0]= $3; #$3 description
			$refFixMods->{$analysisID}{$1}[1]= $2; #$2 delta mass
		}
		elsif ($onMasses==1 && $line =~ /FixedModResidues(\d+)=(\w+)/) {
			$refFixMods->{$analysisID}{$1}[2]= $2; #residue where fix mods occurs
		}
		elsif ($onMasses==1 && $line =~ /delta\d+=(.+)\,(.+\))/) {
			$refmassValueMods->{$analysisID}{$2} = $1;
		}
		elsif ($onQuery==1 && $line =~ /int_min=(.+)/) {
			$localIntMin = $1;
		}
		elsif ($onQuery==1 && $line =~ /int_max=(.+)/) {
			$localIntMax = $1;
		}
		elsif ($onQuery==1 && $line =~ /Ions1=(.+)/ && $fileType !~ /\.PDM/) {
			$localIonSerie = $1 ;
			chomp ($localIonSerie) ;
			my $higtIntensity=0;
			foreach my $Ion (split(/,/,$localIonSerie)) {
				my ($mz,$intens)=split(/:/,$Ion);
				$refexpIonTable->{$analysisID}{$localQuery}{$mz}[0]=sprintf "%.4f",$intens;
				#print "mz=$mz i=$refexpIonTable->{$analysisID}{$localQuery}{$mz}[0]<BR>\n";
			}
			my $minLevel2Match = $intPercent*$localIntMax;
			foreach my $mz (sort{$refexpIonTable->{$analysisID}{$localQuery}{$b}[0]<=>$refexpIonTable->{$analysisID}{$localQuery}{$a}[0]} keys %{$refexpIonTable->{$analysisID}{$localQuery}}){
				if ($refexpIonTable->{$analysisID}{$localQuery}{$mz}[0]>$minLevel2Match) {
					push @{$reffilteredIon->{$analysisID}{$localQuery}},$mz;
					#print "mz=$mz $localQuery $minLevel2Match<BR>\n";
				}else{
				    last;
				}
			}
		}
	}
	#print "<BR>\n";
	close FILE;
	#exit;
}

##########################################################################
####>Get SI information of all the peptides validated in the analysis<####
##########################################################################
sub computeSIPeps{
	my ($analysisID,$refmassValueMods,$refmassValueFixedMods,$refFixMods,$refexpIonTable,$reffilteredIon,$refvalidatedPep)=@_;
	$dbh=&promsConfig::dbConnect('no_user');

	my ($quantifParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE NAME='SI of a peptide'");#QUANTIFICATION_PARAMETER <-> ID_QUANTIF_PARAMETER of SI of a peptide
	my ($instrument) = $dbh->selectrow_array("SELECT INSTRUMENT FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID") ;
	$dbh->disconnect;
	foreach my $pep (keys(%{$refvalidatedPep})) {
		$refvalidatedPep->{$pep}{'SI'}=0;
		$refvalidatedPep->{$pep}{'MSMS'}=0;
		foreach my $queryNum (keys %{$refvalidatedPep->{$pep}{'query'}}){
			$refvalidatedPep->{$pep}{'SI'}+=&compareWithTheoreticalFragmentation($analysisID,$queryNum,$userID,$instrument,$refvalidatedPep->{$pep}{'query'}{$queryNum},$matchInterval,\%{$refexpIonTable->{$analysisID}{$queryNum}},\%{$refmassValueMods->{$analysisID}},\%{$refmassValueFixedMods->{$analysisID}},\%{$refFixMods->{$analysisID}},\@{$reffilteredIon->{$analysisID}{$queryNum}});
			$refvalidatedPep->{$pep}{'MSMS'}+=1;
		}
		#print "$refvalidatedPep->{$pep}{'peptideID'} $refvalidatedPep->{$pep}{'SI'}<BR>\n";
	}

	#$dbh=&promsConfig::dbConnect('no_user');
	#my $sthPepQuantif=$dbh->prepare("INSERT INTO PEPTIDE_QUANTIFICATION (ID_QUANTIFICATION,ID_PEPTIDE,ID_QUANTIF_PARAMETER,QUANTIF_VALUE) VALUES ($quantiID,?,$quantifParamID,?)");
	#foreach my $pep (keys %{$refvalidatedPep}) {
	#	$sthPepQuantif->execute($refvalidatedPep->{$pep}{'peptideID'},$refvalidatedPep->{$pep}{'SI'});
	#	#print "$refvalidatedPep->{$pep}{'peptideID'} $refvalidatedPep->{$pep}{'SI'}\n";
	#}
	#$sthPepQuantif->finish;
	#$dbh->commit;
	#$dbh->disconnect;
	my $projectID=&promsMod::getProjectID($dbh,$analysisID,'analysis');
	mkdir "$promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
	mkdir "$promsPath{quantification}/project_$projectID/quanti_$quantiID";
	open(QUANTI,">$promsPath{quantification}/project_$projectID/quanti_$quantiID/peptide_quantification.txt") || die $!;
	print QUANTI "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n";
	foreach my $pep (keys %{$refvalidatedPep}) {
		print QUANTI "$quantifParamID\t$refvalidatedPep->{$pep}{peptideID}\t$refvalidatedPep->{$pep}{SI}\n";
	}
	close QUANTI;

	return ($quantiID);
}

####################################################################
####>Compute SIN once the SI of all peptides have been computed<####
####################################################################
sub computeSIProts{
	my ($minPepPerProt,$analysisID,$refvalidatedPep,$parentQuantiID)=@_;
	$dbh=&promsConfig::dbConnect('no_user');
	my $sthProtValid=$dbh->prepare("SELECT PEP_SEQ,ID_PROTEIN,PEPTIDE.ID_PEPTIDE FROM PEPTIDE,PEPTIDE_PROTEIN_ATTRIB WHERE PEPTIDE.ID_PEPTIDE=PEPTIDE_PROTEIN_ATTRIB.ID_PEPTIDE AND PEPTIDE.ID_ANALYSIS=PEPTIDE_PROTEIN_ATTRIB.ID_ANALYSIS AND PEPTIDE.ID_ANALYSIS=$analysisID AND PEPTIDE.PEP_RANK IS NOT NULL");
	my $siGi=0;
	my %si;
	$sthProtValid->execute;
	my $visibility;
	while (my ($pepSeq,$proteinID,$peptideID)=$sthProtValid->fetchrow_array) {
		my $varMod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$peptideID,$analysisID,$pepSeq);
		$varMod='' unless $varMod;
		($visibility)=$dbh->selectrow_array("SELECT VISIBILITY FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=$analysisID AND ID_PROTEIN=$proteinID");
		if($si{$proteinID}){
			$si{$proteinID}{'SI'}+=$refvalidatedPep->{"$pepSeq $varMod"}{'SI'};
			$si{$proteinID}{'pepNumber'}+=1;
			if($visibility==2) {#Just top-match proteins contribute to the global SIGI
				$siGi+=$refvalidatedPep->{"$pepSeq $varMod"}{'SI'};
			}
		}
		else{
			$si{$proteinID}{'SI'}=$refvalidatedPep->{"$pepSeq $varMod"}{'SI'};
			$si{$proteinID}{'pepNumber'}=1;
			if($visibility==2) {#Just top-match proteins contribute to the global SIGI
				$siGi+=$refvalidatedPep->{"$pepSeq $varMod"}{'SI'};
			}
		}
	}
	$sthProtValid->finish;
	#print "<BR>SIGI=$siGi<BR>\n";
	if($siGi == 0) {
		$dbh->disconnect;
		open(FILESTAT,">>$fileStat");
		print FILESTAT "Problem during computation of SIGI (equals to 0...)\n";
		print FILESTAT "Ended ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
		close(FILESTAT);
		exit;
	}
	my $sthProtInfo=$dbh->prepare("SELECT IDENTIFIER,PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=?");
	my $sin;

	###> Create new instance of Quantification <-> The parent is the precedent quantification at the level peptide (SI values)
	my ($methodID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='SIN'") || die $dbh->errstr();
	my $quantiAnnot="minPepPerProt=$minPepPerProt;";
	$quantiAnnot=$dbh->quote($quantiAnnot);

	my ($quantiID)=$dbh->selectrow_array("SELECT MAX(ID_QUANTIFICATION)+1 FROM QUANTIFICATION");
	my ($anaName)=$dbh->selectrow_array("SELECT NAME FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
	#$quantiID++;
	$dbh->do("INSERT INTO QUANTIFICATION (ID_QUANTIFICATION,ID_QUANTIFICATION_METHOD,NAME,FOCUS,QUANTIF_ANNOT,UPDATE_DATE,UPDATE_USER,STATUS) VALUES ($quantiID,$methodID,'PROT_SIN_$anaName','protein',$quantiAnnot,NOW(),'$userID',-1)") || die $dbh->errstr();

	$dbh->do("INSERT INTO ANA_QUANTIFICATION (ID_QUANTIFICATION,ID_ANALYSIS) VALUES($quantiID,$analysisID)") || die $dbh->errstr();

	$dbh->do("INSERT INTO PARENT_QUANTIFICATION (ID_PARENT_QUANTIFICATION,ID_QUANTIFICATION) VALUES($parentQuantiID,$quantiID)") || die $dbh->errstr();

	my ($quantifParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE NAME='SIN'");#QUANTIFICATION_PARAMETER <-> ID_QUANTIF_PARAMETER of SIN of a protein

	###> Quantification is 'in-process' state
	$dbh->do("UPDATE QUANTIFICATION SET STATUS=0 WHERE ID_QUANTIFICATION=$quantiID") || die $dbh->errstr();
	$dbh->commit;

	my $sthProtQuantif=$dbh->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_PROTEIN,ID_QUANTIFICATION,ID_QUANTIF_PARAMETER,QUANTIF_VALUE) VALUES (?,$quantiID,$quantifParamID,?)");
	foreach my $proteinID (keys(%si)) {
		next unless $si{$proteinID}{'pepNumber'} >= $minPepPerProt;
		$sthProtInfo->execute($proteinID);
		my ($identifier,$protLength)=$sthProtInfo->fetchrow_array;
		next unless $protLength > 0;
		$sin=$si{$proteinID}{'SI'}/$siGi;
		$sin/=$protLength;
		$sthProtQuantif->execute($proteinID,$sin);
		#print "ID_PROTEIN=$proteinID SI=$si{$proteinID}{'SI'} PepNumber=$si{$proteinID}{'pepNumber'} ProtLength=$protLength SIN=$sin<BR>\n";
	}
	$sthProtInfo->finish;
	$sthProtQuantif->finish;

	###> Quantification is finished
	$dbh->do("UPDATE QUANTIFICATION SET STATUS=1 WHERE ID_QUANTIFICATION=$quantiID") || die $dbh->errstr();

	$dbh->commit;
	$dbh->disconnect;
}

##############################################################
####>Function that computes the theoretical fragmentation<####
##############################################################
sub compareWithTheoreticalFragmentation{
	my ($analysisID,$queryNum,$userID,$instrument,$rank,$matchInterval,$refIonTable,$mVM,$mVFM,$rFM,$refLimitedIonMZ)=@_;
        #print "$userID,$instrument,$rank,$matchInterval,$refIonTable,$mVM,$mVFM,$rFM,$refLimitedIonMZ<BR>\n";
	my %massValueMods=%{$mVM};
	my %massValueFixedMods=%{$mVFM};
	my %FixMods=%{$rFM};
	my %msTypeNames=&promsConfig::getMsType;
	my %massValueAA=&promsConfig::getMassAAmono;
	my %massValueAtom =&promsConfig::getMassATMono;
	my %fragmentDef=&promsConfig::getFragmentDef; #value of fragmentation
	my %fragmentClassif=&promsConfig::getFragmentClassif;
	my (%allowModif,%allowVARModif);
	my %massValueFragMod = (PHOS => 97.976896, OxiM => 64);
	#print "rank=$rank   idquery=$queryNum<BR>\n";
	$dbh=&promsConfig::dbConnect('no_user');
	my ($massData,$infoPep)=$dbh->selectrow_array("SELECT MASS_DATA,INFO_PEP$rank FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID AND QUERY_NUM=$queryNum") ;
	my ($sequence) = $queryToSeq{"$queryNum:$rank"};
	if(!$sequence){
		return 0;
	}
	if($sequence=~/X/){
		return (1,1);
	}
	#print "OK<BR>\n";
	#print "Sequence=$sequence<BR>\n";
	#my ($massObs) = ($massData=~/OBS=(\d+\.*\d*)/);
	#my ($massExp) = ($massData=~/EXP=(\d+\.*\d*)/);
	my $massObs=$massQueryObs{$queryNum};
	my $massExp=$massQueryExp{$queryNum};

	#my ($varModsString,$comments);
	#if ($infoPep=~ /VMOD=\s\+\s([^,]+)/) {
	#		$varModsString = $1;
	#	}
		#if ($infoPep=~ /COM=(.+)/) {
		#	$comments = $1;
		#	chomp ($comments);
		#}
	my $varModsString=$queryToVmod{"$queryNum:$rank"};
	#print "<BR>massObs=$massObs massExp=$massExp vmod=$varModsString\n";
	my ($usrParam,$refFragRules,$deltaParent,$deltaFrag,$tolFrag,$nrLevel) = &promsMod::getInstrParam ($dbh,$instrument,$userID); #$usrParam = -1 : instrument not defined; 0 : defined for all user, 1 user defined
	my %fragmentationRules=%{$refFragRules};
	####>user Preference
	my ($userPref) =$dbh->selectrow_array("SELECT USER_PREF FROM USER_LIST WHERE ID_USER='$userID'");
	$dbh->disconnect;
	$userPref='' unless $userPref;
	###################################
	####>recovering  variable Mods<####
	###################################
	my %varModsTable ; my $warningString;
	if (defined ($varModsString)) {
		$varModsString =~s/^\s\+\s// ;
		my @arrayVarMods = split (/\s\+\s/, $varModsString);
		foreach my $modif (@arrayVarMods) {
			if ($modif=~ /(.+)\s\(([^:]+):([\d\.]+)\)/) {
				my $modifName = $1;
				my $modiCompl = $2;
				my @modifsPos = split(/\./,$3);
				foreach my $modifPos (@modifsPos) {
					if (!defined($massValueMods{"$modifName ($modiCompl)"})) {$warningString .= "unhanded $modifName ($modiCompl) <BR> \n"; next;}
					$varModsTable{$modifPos} = "$modifName ($modiCompl)";
				}
			}
			elsif ($modif=~ /(\w+)\s\(N\-term\)/) {
				my $modifName = $1;
				$varModsTable{0} = "$modifName (N-term)";
			}
			elsif ($modif ne "") {
				$warningString .= "warning :unhanded $modif!<BR>\n"; #error with this modification, theorical calculation will not be good !!!
			}
		}
	}

	my %ionSerieN_term;
	foreach my $fragment (@{$fragmentClassif{"N_term"}}) {
		if ($fragmentationRules{$fragment} && $fragmentationRules{$fragment}>0) {
			$ionSerieN_term{$fragment} = $fragmentDef{$fragment};
			foreach my $neutral (@{$fragmentClassif{"neutral_loss"}}) {#° or *
				if ($fragmentationRules{$fragment.$neutral} && $fragmentationRules{$fragment.$neutral}>0) {
					$ionSerieN_term{$fragment.$neutral} = $fragmentDef{$fragment} + $fragmentDef{$neutral};
				}
			}
		}
	}
	my %ionSerieC_term;
	foreach my $fragment (@{$fragmentClassif{"C_term"}}) {
		if ($fragmentationRules{$fragment} && $fragmentationRules{$fragment}>0) {
			$ionSerieC_term{$fragment} = $fragmentDef{$fragment};
			foreach my $neutral (@{$fragmentClassif{"neutral_loss"}}) { #° or *
				if ($fragmentationRules{$fragment.$neutral} && $fragmentationRules{$fragment.$neutral}>0) {
					$ionSerieC_term{$fragment.$neutral} = $fragmentDef{$fragment} + $fragmentDef{$neutral};
				}
			}
		}
	}
	my %chargeSerie ;
	$chargeSerie{"+"} = 1 if ($fragmentationRules{'+'} && $fragmentationRules{'+'}>0);
	$chargeSerie{"2+"} = 2 if ($fragmentationRules{'2+'} && $fragmentationRules{'2+'}>0);

	###########################
	####>fetching Sequence<####
	###########################
	my @arraySequence ;
	my $allowFragOxiM = 0;
	my $allowFragPhospho = 0;
	$arraySequence[0] ="N_term";
	for (my $i=0; $i<length($sequence); $i++) {
		my $localAA =substr($sequence,$i,1);
		$arraySequence[$i+1]=$localAA;
	}
	push (@arraySequence, "C_term");
	$allowFragOxiM = 0;
	$allowFragPhospho = 0;

	for (my $i=length($sequence)-1; $i>=0 ; $i--) {
		$allowFragOxiM = 1 if (defined($varModsTable{$i+1}) && $varModsTable{$i+1} eq "Oxidation (M)");
		$allowFragPhospho = 1 if (defined($varModsTable{$i+1}) && $varModsTable{$i+1} =~ /Phospho\s\(ST/);
		$allowVARModif{"C_term"}{"PHOS"}[$i+1] = $allowFragPhospho;
		$allowVARModif{"C_term"}{"OxiM"}[$i+1] = $allowFragOxiM;
	}

	###########################################
	####>fetching N-terminal fragmentation<####
	###########################################
	my %massTable ; #  $massTable{ordre}{type de fragmentation}{charge}[]   ; [0] th mass, [2] infomatch [3] exp intensite [4] exp moz
	my $allowLossH2O = 0;
	my $allowLossNH3 = 0;
	my $peptideStartMass=0;
	for (my $i=1 ;$i<$#arraySequence ;$i++) {
		$peptideStartMass += $massValueAA{$arraySequence[$i]} ;
		$peptideStartMass  += $massValueMods{$varModsTable{$i}} if defined($varModsTable{$i}) ;
		$peptideStartMass  += $massValueFixedMods{$arraySequence[$i]} if defined $massValueFixedMods{$arraySequence[$i]};
		if ($i == 1){ # N-term
			$peptideStartMass += $massValueAA{$arraySequence[0]} ;
			$peptideStartMass  += $massValueMods{$varModsTable{0}} if defined($varModsTable{0}) ;
			$peptideStartMass  += $massValueFixedMods{$arraySequence[0]} if defined $massValueFixedMods{$arraySequence[0]};
		}
		if ($i == $#arraySequence-1) { # C-term
			$peptideStartMass  += $massValueMods{$varModsTable{$i+1}} if defined($varModsTable{$i+1}) ;
			$peptideStartMass  += $massValueFixedMods{$arraySequence[$i+1]} if defined $massValueFixedMods{$arraySequence[$i+1]};
		}
		#next if ($arraySequence[$i] eq "C_term" || $arraySequence[$i] eq "N_term");
		$allowLossH2O = 1 if ($allowLossH2O == 0 && $arraySequence[$i] =~ /[STED]/);
		$allowLossNH3 = 1 if ($allowLossNH3 == 0 && $arraySequence[$i] =~ /[NQRK]/);
		foreach my $fragment (keys (%ionSerieN_term)) {
			next if (($allowLossH2O == 0 && $fragment=~/°/) || ($allowLossNH3 == 0 && $fragment=~/\*/));
			foreach my $charge (keys %chargeSerie) {
				$massTable{$i}{$fragment}{$charge}[0] = ($peptideStartMass + $ionSerieN_term{$fragment} + $chargeSerie{$charge}*$massValueAtom{"H"})/$chargeSerie{$charge};
				foreach my $varMod (keys (%{$allowVARModif{"N_term"}})) {
					if ($allowVARModif{"N_term"}{$varMod}[$i] == 1 ) {
					$massTable{$i}{"$fragment-$varMod"}{$charge}[0] = $massTable{$i}{$fragment}{$charge}[0] - $massValueFragMod{$varMod}/$chargeSerie{$charge} }
				}
			}
		}
		if ($i == $#arraySequence-1) { # C-term
			$peptideStartMass += $massValueAA{$arraySequence[$i+1]} ;
		}
	}
	my $peptideMass = $peptideStartMass;
	###########################################
	####>fetching C-terminal fragmentation<####
	###########################################
	my $peptideEndmass =0;
	$allowLossH2O = 0;
	$allowLossNH3 = 0;
	for (my $i=$#arraySequence-1 ;$i>0 ;$i--)  {
		my $Endpos = (@arraySequence-1)-$i;
		$peptideEndmass += $massValueAA{$arraySequence[$i]} ;
		$peptideEndmass += $massValueMods{$varModsTable{$i}} if defined($varModsTable{$i}) ;
		$peptideEndmass += $massValueFixedMods{$arraySequence[$i]} if defined($massValueFixedMods{$arraySequence[$i]});
		if ($i==$#arraySequence-1) { # C-term
			$peptideEndmass += $massValueAA{$arraySequence[$i+1]};
			$peptideEndmass += $massValueMods{$varModsTable{$i+1}} if defined($varModsTable{$i+1});
			$peptideEndmass += $massValueFixedMods{$arraySequence[$i+1]} if defined($massValueFixedMods{$arraySequence[$i+1]});
		}
		if ($i==1) { # N-term
			#$peptideEndmass += $massValueAA{$arraySequence[0]};
			$peptideEndmass += $massValueMods{$varModsTable{0}} if defined($varModsTable{0}) ;
			$peptideEndmass += $massValueFixedMods{$arraySequence[0]} if defined($massValueFixedMods{$arraySequence[0]});
		}
		#next if ($arraySequence[$i] eq "C_term") || ($arraySequence[$i] eq "N_term") ;
		$allowLossH2O = 1 if ($allowLossH2O == 0 && $arraySequence[$i] =~ /[STED]/);
		$allowLossNH3 = 1 if ($allowLossNH3 == 0 && $arraySequence[$i] =~ /[NQRK]/);
		foreach my $fragment (keys %ionSerieC_term) {
			next if (($allowLossH2O == 0 && $fragment=~/°/) || ($allowLossNH3 == 0 && $fragment=~/\*/));
			foreach my $charge (keys %chargeSerie) {
				$massTable{$Endpos}{$fragment}{$charge}[0] = ($peptideEndmass+ $ionSerieC_term{$fragment} + $chargeSerie{$charge}*$massValueAtom{"H"})/$chargeSerie{$charge};
				foreach my $varMod (keys %{$allowVARModif{"C_term"}}) {
					if ($allowVARModif{"C_term"}{$varMod}[$i] == 1) {
						$massTable{$Endpos}{"$fragment-$varMod"}{$charge}[0] = $massTable{$Endpos}{$fragment}{$charge}[0] - $massValueFragMod{$varMod}/$chargeSerie{$charge};
					}
				}
			}
		}
		if ($i==1) { # N-term
			$peptideEndmass += $massValueAA{$arraySequence[0]};
		}
	}
	####>Looking for matched fragmentation in ionic serie<####
        #print "<BR>Autre Methode<BR>\n";
    my %used;
	foreach my $number (keys %massTable) {
		foreach my $modif (keys %{$massTable{$number}}) {
			foreach my $charge (keys %{$massTable{$number}{$modif}}) {
				foreach my $mz (@{$refLimitedIonMZ}) {
					if (abs($massTable{$number}{$modif}{$charge}[0]-$mz)<$matchInterval)  {
						$massTable{$number}{$modif}{$charge}[2] = 1;  # %massTable [2] infomatch
						if (!$massTable{$number}{$modif}{$charge}[3] || $massTable{$number}{$modif}{$charge}[3]<$refIonTable->{$mz}[0]) {
                            if(!$used{$mz}){#Avoid to use twice the intensity of an ion
								$massTable{$number}{$modif}{$charge}[3]=$refIonTable->{$mz}[0];
								$massTable{$number}{$modif}{$charge}[4]=$mz;
								$used{$mz}++;
								#print "Masse=$massTable{$number}{$modif}{$charge}[0]\tIntensite=$refIonTable->{$mz}[0]<BR>\n";
							}
						}
					}
				}
			}
		}
	}

	#print "Compute SIn<BR>\n";
	my $sumIntensities=0;#Sum of all the intensities of ions that match the sequence (blue one in peptide_view)
	my $matchIntensities=0;
	for (my $i=0 ;$i<@arraySequence-1 ;$i++) {
            #print "$i=$arraySequence[$i]<BR>\n";
		next if ($arraySequence[$i] eq "C_term" || $arraySequence[$i] eq "N_term");
		foreach my $modif (sort{if ($a =~ /b/ && $b !~ /b/) {return -1;} elsif ($b =~ /b/ && $a !~ /b/) {return 1 ; } else { $a cmp $b ;}} keys %ionSerieN_term) {
			#print "modif=$modif<BR>\n";
            foreach my $charge (sort{$a cmp $b} keys %chargeSerie) {
				if (defined($massTable{$i}{$modif}{$charge}[0]) && defined($massTable{$i}{$modif}{$charge}[2]) && defined($massTable{$i}{$modif}{$charge}[3])) {
					#printf$massTable{$i}{$modif}{$charge}[0], "%.2f\t",$massTable{$i}{$modif}{$charge}[0],"<BR>\n";
					$sumIntensities+=$massTable{$i}{$modif}{$charge}[3];
					$matchIntensities+=1;
					#print "N_TERM\t$massTable{$i}{$modif}{$charge}[0]\t$massTable{$i}{$modif}{$charge}[3]<BR>\n";
                    delete($massTable{$i}{$modif}{$charge}[3]);
				}
			}
		}
		my $endPos = (@arraySequence-1)- $i ;
		foreach my $modif (sort{$a cmp $b} keys %ionSerieC_term) {
                    #print "modif2=$modif<BR>\n";
			foreach my $charge (sort{$a cmp $b} keys %chargeSerie) {
				if (defined($massTable{$endPos}{$modif}{$charge}[0]) && defined($massTable{$endPos}{$modif}{$charge}[2]) && defined($massTable{$endPos}{$modif}{$charge}[3])) {
					$sumIntensities+=$massTable{$endPos}{$modif}{$charge}[3];
					$matchIntensities+=1;
					#print "C_TERM\t$massTable{$endPos}{$modif}{$charge}[0]\t$massTable{$endPos}{$modif}{$charge}[3]<BR>\n";
					delete($massTable{$i}{$modif}{$charge}[3]);
				}
			}
		}
	}
        #print "<BR>$sumIntensities\n";
	#exit();
	return ($sumIntensities);
}

####>Revision history<####
# 1.2.0 Peptide SI data now written to file $promsPath{quantification}/project_$projectID/quanti_$quantifID/peptide_quantification.txt (PP 11/05/18)
# 1.1.0 Change query to avoid retrieve ghost-peptides (GA 03/02/17)
# 1.0.9 Minor modification in args (GA 16/04/14)
# 1.0.8 system command removal (PP 13/11/13)
# 1.0.7 Remove ID_ANALYSIS value from insert in PROTEIN_QUANTIFICATION (PP 12/09/13)
# 1.0.6 Remove VAR_MOD from script (GA 27/05/13)
# 1.0.5 Minor changes to print in filestat some information (GA 17/10/12)
# 1.0.4 Minor changes in computeSIProts & checks for uninitialized $fragmentationRules{n+} (PP 30/04/12)
# 1.0.3 Minor update -> automatic Name for ProteinQuantification and STATUS handling (GA 30/03/12)
# 1.0.2 Update the script so as to avoid to use QUERY_VALIDATION table (GA 15/07/2011)<BR>Scan the mascot dat file instead to find the same peptides<BR>Modify the MYSQL in PEPTIDE_TABLE to avoid to select fake PEPTIDES (coming from XIC quantification)
# 1.0.1 Minor update: add ID_ANALYSIS in the insert of PROTEIN_QUANTIFICATION (GA 08/07/2011)
# 1.0.0 New script to compute SIN quantification
