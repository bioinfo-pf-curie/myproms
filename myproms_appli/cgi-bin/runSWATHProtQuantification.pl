#!/usr/local/bin/perl -w

################################################################################
# runSWATHProtQuantification.pl       1.3.0                                    #
# Component of site myProMS Web Server                                         #
# Authors: M.Le Picard (Institut Curie)                                        #
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
use promsConfig;
use promsMod;
use promsQuantif;
use strict;
#use XML::SAX;
use MIME::Base64;
use POSIX qw(strftime);
#use File::Copy;
use File::Copy::Recursive qw(dirmove dircopy);
#use File::Path qw(rmtree); # remove_tree
#use Data::Dumper;
#exit;
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo('no_user');
my %cluster=&promsConfig::getClusterInfo;#('debian'); # default is 'centos'


##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect('no_user');


###############################
####>Recovering parameters<####
###############################
my ($quantifID,$quantifDate)=@ARGV;
my $quantifDir="$promsPath{tmp}/quantification/$quantifDate";
my $projectID=&promsMod::getProjectID($dbh,$quantifID,'QUANTIFICATION');

$dbh->do("UPDATE QUANTIFICATION SET STATUS=0 WHERE ID_QUANTIFICATION=$quantifID");
$dbh->commit;

my $fileStat="$quantifDir/status_$quantifID.out";
open(FILESTAT,">$fileStat") || die "Error while opening $fileStat";
print FILESTAT "Started ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
print FILESTAT "1/3 Generating data files\n";
close FILESTAT;

my %quantifParameters=&promsQuantif::getQuantificationParameters("$quantifDir/quantif_info.txt");

my ($pepSpecificity,$pepMissedCleavage,$ptmFilter,$pepChargeState,$pepSource)=@{$quantifParameters{'DB'}{'PEPTIDES'}};


my $runDir="$quantifDir/quantif_$quantifID";
my $dataDir="$runDir/data";
my $resultDir="$runDir/results";
my $graphDir="$resultDir/graph";
mkdir $runDir unless -e $runDir;
mkdir $dataDir unless -e $dataDir;
mkdir $resultDir unless -e $resultDir;
mkdir $graphDir unless -e $graphDir;

&promsQuantif::writeQuantifParameterFiles($dataDir,$quantifParameters{'R'});    # $dataDir: full path to tmp quantif dir/data/


#####################################################
####>Protein lists for filtering & normalization<####
#####################################################
my ($protSelectionType,%selectExcludeProteins);
my $sthList=$dbh->prepare("SELECT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=? ORDER BY ID_PROTEIN");

###>Selection/exclusion list
if ($quantifParameters{'DB'}{'PROTEINS'}) {
    ($protSelectionType,my $listID)=@{$quantifParameters{'DB'}{'PROTEINS'}};
    $listID=~s/#//; # remove id flag
    $sthList->execute($listID);
    while (my ($protID)=$sthList->fetchrow_array) {
        $selectExcludeProteins{$protID}=1;
    }
}

###>Global standard
if ($quantifParameters{'DB'}{'NORMALIZATION_PROTEINS'}) { # globalStandards
    my $listID=$quantifParameters{'DB'}{'NORMALIZATION_PROTEINS'}[0];
    $listID=~s/#//; # remove id flag
    open(GLOBSTAND,">$dataDir/globalStandards.txt") or die $!;
    $sthList->execute($listID);
    while (my ($protID)=$sthList->fetchrow_array) {
        print GLOBSTAND $protID,"\n";
    }
    close GLOBSTAND;
}

$sthList->finish;

###################################################
####>Get the conditions associated with quanti<####
###################################################
my (%ana2Quant,%stateToCond,%conditionList,%condPeptideList); # %biorep2cond,
my ($pos,$run)=(0,0);
#foreach my $condID (@{$quantifParameters{'DB'}{'STATES'}}){ # DO NOT process Replicates here BUT from QUANTIF_ANNOT!!! PP
#    $pos++;
#    my $cond="State$pos";
#    $stateToCond{$condID}=$cond;
#    my $bioRepPos=0;
#    foreach my $bioReplicate (@{$quantifParameters{'DB'}{"QUANTITOCOND_$condID"}}){
#        $run++;
#        $bioRepPos++;
#        my ($observationID,$pepQuantifID,$analysisID)=split(/:/,$bioReplicate);
#        my $bioRepName=$cond.".BioRep".$bioRepPos;
#        $ana2Quant{$analysisID}{$pepQuantifID}=$bioRepName;
#        $conditionList{$cond}{$bioRepName}=$run;
#        $biorep2cond{$bioRepName}=$cond;
#    }
#}

my $numConditions=0;
my $ratioString;
my ($quantifAnnot,$quantifMethodID,$quantifiedModifID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_QUANTIFICATION_METHOD,ID_MODIFICATION FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID");
foreach my $annotStrg (split(/::/,$quantifAnnot)) {
    #### STATES ####
	if ($annotStrg =~ /^STATES=(.+)/) { # 'state1;state2;---;stateN',  stateX='nbBioRep,bioRep1.bioRep2.---.bioRepN,#condID', bioRepX='techRep1=techRep2=---=techRepN', techRepX='frac1+frac2+---+fracN'
        my $run=0;
        foreach my $state (split(/;/,$1)) {
            $numConditions++;
            if ($state=~/\+/) { # fractions!!! NOT HANDLED!
                $dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID"); # Failed
                $dbh->commit;
                $dbh->disconnect;
                die "ERROR in State #$pos: Fractions not compatible with DIA quantification. Correct your design.";
            }
            my $cond="State$numConditions";
            $state=~s/#//g; # remove all id tags
			my ($numBioRep,$quantiObsIDs,$condID)=split(/,/,$state);
            $stateToCond{$condID}=$cond;
            my $bioRepPos=0;
            foreach my $bioReplicate (split(/\./,$quantiObsIDs)) {
                $bioRepPos++;
                my $bioRepName=$cond.".BioRep".$bioRepPos;
                #$biorep2cond{$bioRepName}=$cond;
                my $techRepPos=0;
                foreach my $techReplicate (split(/&/,$bioReplicate)) {
                    $techRepPos++;
                    my $techRep='.TechRep'.$techRepPos;
                    $run++;
                    my ($observationID,$pepQuantifID,$analysisID)=split(/:/,$techReplicate);

                    $ana2Quant{$analysisID}{STATE}=$cond;
                    $ana2Quant{$analysisID}{PEP_QUANTIF}=$pepQuantifID;
                    $ana2Quant{$analysisID}{BIO_REP}=$bioRepName;
                    $ana2Quant{$analysisID}{TECH_REP}=$techRep;
                    $ana2Quant{$analysisID}{RUN}=$run;

                    $conditionList{$cond}{$bioRepName}{$techRep}=$run;
                }
            }
        }
    }
    #### RATIOS ####
	elsif ($annotStrg=~/^RATIOS=(.+)/) { # delay ratio processing because $stateToCond not yet defined
        $ratioString=$1;
        $ratioString=~s/#//g;
    }
}

###>Process ratio string
#my %ratioToPos;
#my $ratioPos=0;
#foreach my $ratio (split(/;/,$ratioString)){ # #test/#reference
#    $ratio=~s/#//g;
#    $ratioPos++;
#    my @cond=split(/\//,$ratio);
#    $ratioToPos{"$stateToCond{$cond[0]}/$stateToCond{$cond[1]}"}=$ratioPos;
#}


###>For modification based quantification
my ($ambiguousModifPos,$qQuantifModifStrg)=(0,'');
if ($quantifiedModifID){
    ###> Recovering targetted residues
    my $sthSpecif=$dbh->prepare("SELECT SPECIFICITY FROM ANALYSIS_MODIFICATION WHERE ID_ANALYSIS=? AND ID_MODIFICATION=$quantifiedModifID AND MODIF_TYPE='V'");
    my %resList;
    foreach my $anaID (keys %ana2Quant) {
        $sthSpecif->execute($anaID);
        my ($specifStrg)=$sthSpecif->fetchrow_array;
        foreach my $res (split(/,/,$specifStrg)){$resList{$res}=1;}
    }
    $sthSpecif->finish;
    my $quantifModifStrg=join('',sort keys %resList);
    $qQuantifModifStrg=quotemeta($quantifModifStrg); # \*\+
}


######################################
####>Fetching quantification data<####
######################################

###>Fetching protein and peptide data
my $selectedPepQuery=qq
|SELECT GROUP_CONCAT(DISTINCT PPA.ID_PROTEIN),P.ID_PEPTIDE,ABS(PEP_BEG),PEP_SEQ,GROUP_CONCAT(PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),CHARGE
FROM PEPTIDE P
LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
INNER JOIN ANALYSIS_PROTEIN AP ON PPA.ID_PROTEIN=AP.ID_PROTEIN
WHERE P.ID_ANALYSIS=? AND AP.ID_ANALYSIS=? AND AP.VISIBILITY=2|;
$selectedPepQuery.=($pepSpecificity eq 'unique')? ' AND PPA.IS_SPECIFIC=1' : ($pepSpecificity eq 'unique_shared')? ' AND AP.PEP_SPECIFICITY=100' : '';
$selectedPepQuery.=' AND MISS_CUT=0' unless $pepMissedCleavage;
$selectedPepQuery.=" GROUP BY PPA.ID_PROTEIN,P.ID_PEPTIDE ORDER BY AP.NUM_PEP DESC,ABS(PEP_BEG) ASC;";

my $sthSelectedPeptides=$dbh->prepare("$selectedPepQuery");

$ptmFilter=~s/#//g; # remove id flag
my ($ptmAllowed,@selectedPTMs)=split(/[:,]/,$ptmFilter);
my %allowedPtmID;
foreach my $modID (@selectedPTMs) {$allowedPtmID{$modID}=1;}
$allowedPtmID{$quantifiedModifID}=1 if $quantifiedModifID;

my (%validProteins,%pep2prot,%peptideList,%peptideBeg,%peptideUsed);
foreach my $analysisID (keys %ana2Quant){
    $sthSelectedPeptides->execute($analysisID,$analysisID);
    my %excludedSeq;
    while (my ($protID,$pepID,$pepBeg,$pepSeq,$varModStrg,$charge)=$sthSelectedPeptides->fetchrow_array) {
        #>Protein filtering
        if ($protSelectionType) {
            if ($protSelectionType eq 'exclude') {next if $selectExcludeProteins{$protID};}
            else {next unless $selectExcludeProteins{$protID};} # restrict
        }
        next if ($quantifiedModifID && (!$varModStrg || $varModStrg!~/(^|&)$quantifiedModifID:/));


        if ($pepSpecificity ne 'unique'){
            next if ($peptideUsed{$pepSeq} && $peptideUsed{$pepSeq}!=$protID);      #use peptide only once in non proteotypic case
            $peptideUsed{$pepSeq}=$protID;
        }

        #>Checking allowed mod
        if ($varModStrg && $ptmAllowed != 1) {
            my $hasPTM=0;
            foreach my $modif (split (/&/,$varModStrg)){
                my ($modifID)=($modif=~/^(\d+)/);
                if (!$allowedPtmID{$modifID}) {
                    $hasPTM=1;
                    last;
                }
            }
            if ($hasPTM) { # this modification is not allowed
                $excludedSeq{$pepSeq}=1 if $ptmAllowed <= -1; # -1 or -2 unmodified peptide not allowed if exists a modified version
                $excludedSeq{$pepID}=1 if $ptmAllowed >= 0; # 0 (exclude all) or 2 (Allow selected)
                next;
            }
        }
        if ($varModStrg) {
            $peptideList{$pepID}="$pepSeq&$varModStrg\_$charge" unless ($excludedSeq{$pepSeq} || $excludedSeq{$pepID});
        }
        else{$peptideList{$pepID}=$pepSeq.'_'.$charge unless ($excludedSeq{$pepSeq} || $excludedSeq{$pepID});}
        push @{$validProteins{$protID}{$analysisID}},$pepID unless ($excludedSeq{$pepSeq} || $excludedSeq{$pepID});
        push @{$pep2prot{$pepID}},$protID unless ($excludedSeq{$pepSeq} || $excludedSeq{$pepID});
        $peptideBeg{$protID}{$peptideList{$pepID}}=$pepBeg unless ($excludedSeq{$pepSeq} || $excludedSeq{$pepID});
    }
    $sthSelectedPeptides->finish;

    ###> Recover fragments info
    my $pepQuantifID=$ana2Quant{$analysisID}{PEP_QUANTIF};
    my $fragmentFile="$promsPath{peptide}/proj_$projectID/ana_$analysisID/swath_quanti_$pepQuantifID.txt";
    open(IN,$fragmentFile) or die "$fragmentFile : $!";
    while (<IN>) {
        next if $_=~m/PEP/;
        s/\s\Z//g; # chom is not enough. More hidden Windows character
        my ($peptideID,$fragMZ,$fragCharge,$fragType,$fragRes,$fragArea)=split(/!/,$_);
        if ($fragArea=~/^(N\/A|None|0)/) {
            $fragArea='NA';
        }
        next if (!$peptideList{$peptideID});

        my ($peptideSeq,$peptideCharge)=split(/_/,$peptideList{$peptideID});
        foreach my $protID (@{$pep2prot{$peptideID}}){
            #$condPeptideList{$protID}{"$peptideSeq\_$peptideCharge"}{"$fragType$fragRes\_$fragCharge"}{ $biorep2cond{$ana2Quant{$analysisID}{$pepQuantifID}} }{ $ana2Quant{$analysisID}{$pepQuantifID} }=$fragArea;
            # condPeptideList{protID}{pepIon}{fragIon}{State}{bioRep}{tehcRep}=fragArea
            $condPeptideList{$protID}{"$peptideSeq\_$peptideCharge"}{"$fragType$fragRes\_$fragCharge"}{ $ana2Quant{$analysisID}{STATE} }{ $ana2Quant{$analysisID}{BIO_REP} }{ $ana2Quant{$analysisID}{TECH_REP} }=$fragArea;
        }
    }
    close IN;
}

my $numDataLines=0;
########################################
####> Printing data into table.txt <####
########################################
my (%pepList, %pepTotal);
open(DATA,">$dataDir/table.txt") or die $!;     # Name of the table that is going to be processed by R scripts and statistically analyzed...
print DATA "ProteinName\tPeptideSequence\tPrecursorCharge\tFragmentIon\tProductCharge\tIsotopeLabelType\tCondition\tBioReplicate\tRun\tIntensity\n";
foreach my $proteinID (keys %condPeptideList){
    my $proteinModif=$proteinID;
    foreach my $peptide (keys %{$condPeptideList{$proteinID}}){
        my ($peptideSeq,$peptideCharge)=split(/_/,$peptide);
        my @seq=split(//,$peptideSeq);

        if ($quantifiedModifID){
            my @pos;
            my ($posStrg)=($peptideSeq=~/&$quantifiedModifID:([\d+\.\-\=\+\*]+)/);
            foreach my $pepPos (split(/\./,$posStrg)){
                if($pepPos=~/\d+/){push @pos, $seq[$pepPos-1].($peptideBeg{$proteinID}{$peptide}+$pepPos-1);}
                if($pepPos eq '-'){push @pos,'n0';}
                if($pepPos eq '='){push @pos,'n'.$peptideBeg{$proteinID}{$peptide};}
                if($pepPos eq '*'){push @pos,'c'.($peptideBeg{$proteinID}{$peptide}+$#seq);}
                if($pepPos eq '+'){push @pos,'c0';}
            }
            $proteinModif=$proteinID.'-'.join('.',@pos);
        }
        $pepTotal{$proteinID}++;
        foreach my $fragment (keys %{$condPeptideList{$proteinID}{$peptide}}){
            my ($fragType,$fragCharge)=split(/_/,$fragment);
            foreach my $condition (keys %conditionList){
                foreach my $bioReplicate (keys %{$conditionList{$condition}}){

                     foreach my $techRep (keys %{$conditionList{$condition}{$bioReplicate}}){
                        my $run=$conditionList{$condition}{$bioReplicate}{$techRep};
                        my $fragArea='NA'; # default
                        if ($condPeptideList{$proteinID}{$peptide}{$fragment}{$condition} && $condPeptideList{$proteinID}{$peptide}{$fragment}{$condition}{$bioReplicate} && $condPeptideList{$proteinID}{$peptide}{$fragment}{$condition}{$bioReplicate}{$techRep}) {
                            $fragArea=$condPeptideList{$proteinID}{$peptide}{$fragment}{$condition}{$bioReplicate}{$techRep} || 'NA';
                            #if ($fragArea eq '') {$fragArea='NA';}
                            #$numPepValues++;
                        }
                        print DATA "$proteinModif\t$peptideSeq\t$peptideCharge\t$fragType\t$fragCharge\tL\t$condition\t$bioReplicate\t$run\t$fragArea\n";
                        $pepList{$proteinModif}{$condition}{"$peptideSeq\_$peptideCharge"}{$fragment}=1 unless $fragArea eq 'NA';
                        $numDataLines++;
                      }

                }
            }
        }
    }
}
close DATA;


###########################
###> Running R Scripts <###
###########################
$dbh->disconnect;
open(FILESTAT,">>$fileStat");
print FILESTAT "2/3 Running quantification\n";
close(FILESTAT);
my $pathR=($cluster{'on'})? $cluster{'path'}{'R'} : $promsPath{'R'};
open(R_SCRIPT,">$runDir/quantiswath.R");
print R_SCRIPT qq
|
#########################################
# Launcher for quantification R scripts #
#########################################

filepath="$promsPath{R_scripts}/"
source(paste(filepath,"quantiSwath.R",sep=""))
|;
close R_SCRIPT;
my $RcommandString="export LANG=en_US.UTF-8; cd $runDir; $pathR/R CMD BATCH --no-save --no-restore quantiswath.R";
my $RoutFile='quantiswath.Rout';

### ALWAYS SKIP CLUSTER (PP 17/09/18)!!! ###
###if ($cluster{'on'}) {
###	my $numCPU=$cluster{'maxCPUs'} || 1;
###	my $clusterCommandString=$cluster{'buildCommand'}->($runDir,$RcommandString);
###	my $bashFile = "$runDir/runSWATHProtQuantification.sh";
###	my $maxHours=96;
###	#my $maxMem=int(1.5 + 1E-6 * $numPepValues);
###    #$maxMem*=5 if $cluster{'name'} eq 'CentOS';
###	my $maxMem=int(1E-5 * $numDataLines);
###    $maxMem=($maxMem<10)? 10 : ($maxMem > $cluster{'maxMem'})? $cluster{'maxMem'} : $maxMem;
###	$maxMem.='Gb';
###	open (BASH,">$bashFile");
###	print BASH qq
###|#!/bin/bash
#####resources
####PBS -l mem=$maxMem
####PBS -l nodes=1:ppn=$numCPU
####PBS -l walltime=$maxHours:00:00
####PBS -q batch
#####Information
####PBS -N protQuant_$quantifID
#####PBS -M marine.le-picard\@curie.fr
####PBS -m abe
####PBS -o $runDir/PBS.txt
####PBS -e $runDir/PBSerror.txt
###
##### Command
###$clusterCommandString
###echo _END_$quantifID
###|;
###
###	close BASH;
###	my $modBash=0775;
###	chmod $modBash, $bashFile;
###	$cluster{'sendToCluster'}->($bashFile);
###
###
###	###>Waiting for R job to run
###	my $pbsError;
###	my $nbWhile=0;
###	my $maxNbWhile=$maxHours*60*2;
###	while ((!-e "$runDir/PBS.txt" || !`tail -3 $runDir/PBS.txt | grep _END_$quantifID`) && !$pbsError) {
###		if ($nbWhile > $maxNbWhile) {
###            $dbh=&promsConfig::dbConnect('no_user'); # reconnect
###			$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID"); # Failed
###			$dbh->commit;
###			$dbh->disconnect;
###			die "Aborting quantification: R is taking too long or died before completion";
###		}
###		sleep 30;
###        $pbsError=$cluster{'checkError'}->("$runDir/PBSerror.txt");
###		#$pbsError=`head -5 $runDir/PBSerror.txt` if -e "$runDir/PBSerror.txt";
###		$nbWhile++;
###	}
###
###	# Wait 15 min max to make sure Rout file is visible from web server
###	foreach (1..30) {
###		if (-e "$runDir/$RoutFile") {last;}
###		sleep 30;
###	}
###}
###else { ###>Run job on Web server
	system $RcommandString;
###}
sleep 3;

if (-e "$resultDir/sessionInfo.txt") {
    my $versionLine=`grep MSstats_ $resultDir/sessionInfo.txt`;
    my ($version)=($versionLine=~/MSstats_(\S+)/);
    if ($version) {
        $dbh=&promsConfig::dbConnect('no_user'); # reconnect
        my $MSstatsStrg='MSstats;'.$version.'::';
        $dbh->do("UPDATE QUANTIFICATION SET QUANTIF_ANNOT=REPLACE(QUANTIF_ANNOT,'MSstats::','$MSstatsStrg') WHERE ID_QUANTIFICATION=$quantifID");
        $dbh->commit;
        $dbh->disconnect;
    }
}

####>ERROR Management<####
my $RoutStrg=`tail -3 $runDir/$RoutFile`;
unless ($RoutStrg=~/proc\.time\(\)/) {
	$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID");
	$dbh->commit;
	$dbh->disconnect;
	$RoutStrg=`tail -20 $runDir/$RoutFile`;
	my $RerrorStrg="R script has generated an error!";
	my $inError=0;
	foreach my $line (split(/\n/,$RoutStrg)) {
		next if (!$inError && $line !~ /^Error in/); # skip everything before "Error in..."
		$inError=1;
		$RerrorStrg.="\n$line";
	}
	die $RerrorStrg;
}



									############################################
			#############################> PARSING QUANTIFICATION RESULTS <##########################
									############################################

open(FILESTAT,">>$fileStat");
print FILESTAT "3/3 Parsing results\n";
close FILESTAT;

$dbh=&promsConfig::dbConnect('no_user'); # reconnect

my %quantifParamIDs;
my $sthQuantifParam=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,CODE FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=$quantifMethodID");
$sthQuantifParam->execute;
while (my ($paramID,$paramCode)=$sthQuantifParam->fetchrow_array) {
    $quantifParamIDs{$paramCode}=$paramID;
}
$sthQuantifParam->finish;

my ($sthInsModRes,$sthInsProtRes);
my $sthInsProt=$dbh->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_QUANTIFICATION,ID_PROTEIN,ID_QUANTIF_PARAMETER,QUANTIF_VALUE,TARGET_POS) VALUES ($quantifID,?,?,?,?)")  or die "Couldn't prepare statement: " . $dbh->errstr;
if ($quantifiedModifID) {
	$sthInsModRes=$dbh->prepare("INSERT INTO MODIFIED_RESIDUE (ID_QUANTIFICATION,RESIDUE,POSITION) VALUES ($quantifID,?,?)");
	$sthInsProtRes=$dbh->prepare("INSERT INTO PROTQUANTIF_MODRES (ID_MODIF_RES,ID_PROT_QUANTIF) VALUES (?,?)");
}

###>Generating ratio codes used by R<###
my (%measurePos2Code,%measureCode2RatioPos,%ratioPos2states); # a state can be involved in multiple ratios
my $rPos=0;
foreach my $s1 (1..$numConditions-1) {
	foreach my $s2 ($s1+1..$numConditions) {
		$rPos++;
		$measurePos2Code{$rPos}="State$s2-State$s1"; # for ResultsDAProt.txt
		@{$ratioPos2states{$rPos}}=("State$s1","State$s2");
	}
	last if $quantifParameters{'DB'}{'SINGLE_REF'};
}

##################################
####> Inserting data into DB <####
##################################
my %validProt;
#my %protModRes;
my %modResiduesID; # keeps track of modif residues already inseted in MODIFIED_RESIDUE to prevent duplicate inflation
if (-e "$resultDir/ResultsDAProt.txt") {
	my @RDbCodeList=(['pval','PVAL_ADJ'],['SD','SD_GEO']);
    my %protHeader2Idx;
	#my %refColumns;
    open(PROT,"$resultDir/ResultsDAProt.txt") or die "[$resultDir/ResultsDAProt.txt] $!";
    while (<PROT>) {
		$_=~s/\s*\Z//;
		my @lineData=split(/\t/,$_);
		####>Header------------------
		if ($.==1) { # 1st line of the file
			my $colIdx=0;
			foreach my $colName (@lineData) {
				$protHeader2Idx{$colName}=$colIdx;
				$colIdx++;
			}
			next;
		}
		####>Data------------------
		my $protID0=$lineData[0]; # $lineData[$protHeader2Idx{'Protein'}];
		next if !$protID0; # just to be safe
		next if ($quantifiedModifID && $protID0 !~/-/); # must be a mod prot
		my ($protID,@modResidues)=($quantifiedModifID)? split(/[-.]/,$protID0) : ($protID0,undef);
		&promsQuantif::insertModifiedResidues($protID,\@modResidues,\%modResiduesID,$dbh,$sthInsModRes) if $quantifiedModifID; # <-- %modResiduesID should be full at <PROT1>

		###>Conditions ratio data
		foreach my $targetPos (sort{$a<=>$b} keys %measurePos2Code) {
			#>Fold change
			my $log2FC=$lineData[$protHeader2Idx{'log2FC_'.$measurePos2Code{$targetPos}}];
			next if $log2FC=~/NA/i; # NA or NaN
			my $fcValue=($log2FC eq '-Inf')? 0.001 : ($log2FC eq 'Inf')? 1000 : 2**$log2FC; #exp($log2FC*$log2);
			$sthInsProt->execute($protID,$quantifParamIDs{'RATIO'},$fcValue,$targetPos);
			if ($quantifiedModifID) {
				my $protQuantifID=$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF');
				#&promsQuantif::insertModifiedResidues($protID,\@modResidues,\%modResiduesID,$dbh,$sthInsModRes);
				&promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$protQuantifID);
			}
			$validProt{$protID0}{$targetPos}=1;
			#>Associated parameters
			foreach my $refParam (@RDbCodeList) {
				my ($paramR,$paramDB)=@{$refParam};
				my $paramValue=$lineData[$protHeader2Idx{$paramR.'_'.$measurePos2Code{$targetPos}}];
				next if $paramValue=~/NA/i; # NA or NaN
				$sthInsProt->execute($protID,$quantifParamIDs{$paramDB},$paramValue,$targetPos);
				&promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $quantifiedModifID;
			}
		} # next ratioPos

	} # next line
	close PROT;
}
else { # File not read (does not exist)
    $dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID");
    $dbh->commit;
    $dbh->disconnect;
    die "Missing result file !\n";
}



################################
### Computing peptide counts ###
################################
####> read CVoutPep.txt to recover outliers
if (-e "$resultDir/CVoutPep.txt") {
    open(OUTLIERS, "<$resultDir/CVoutPep.txt") or die $!;
    my $line;
    while ($line=<OUTLIERS>) {
        next if $.==1;
        my ($protID,$pep,$state,$value)=split(/\s/,$line);
        my ($pepSeq,$pepCharge,$pepIon,$ionCharge)=split(/_/,$pep);
        delete $pepList{$protID}{$state}{"$pepSeq\_$pepCharge"}{"$pepIon\_$ionCharge"};
        delete $pepList{$protID}{$state}{"$pepSeq\_$pepCharge"} if scalar keys %{$pepList{$protID}{$state}{"$pepSeq\_$pepCharge"}} == 0;
    }
    close OUTLIERS;
}


####> read excluded.txt
if (-e "$dataDir/excluded.txt") {
    open(EXCLUDED, "<$dataDir/excluded.txt") or die $!;
    my $line;
    while ($line=<EXCLUDED>) {
        next if $.==1;
        my @lineInfo=split(/\s/,$line);
        delete $pepList{$lineInfo[0]}{$lineInfo[6]}{"$lineInfo[1]\_$lineInfo[2]"}{"$lineInfo[3]\_$lineInfo[4]"};    ##$lineInfo[0] -> protID ; $lineInfo[1] -> peptide seq ; $lineInfo[2] -> peptide charge ; $lineInfo[3] -> fragment ; $lineInfo[4] -> fragment charge ; $lineInfo[6] -> state
        delete $pepList{$lineInfo[0]}{$lineInfo[6]}{"$lineInfo[1]\_$lineInfo[2]"} if scalar keys %{$pepList{$lineInfo[0]}{$lineInfo[6]}{"$lineInfo[1]\_$lineInfo[2]"}} == 0;
    }
    close EXCLUDED;
}
else { # File not read (does not exist)
    $dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID");
    $dbh->commit;
    $dbh->disconnect;
    die "Missing file 'excluded.txt'!\n";
}


foreach my $protID0 (keys %pepList){
    next if !$validProt{$protID0};
    my ($protID,@modResidues)=($quantifiedModifID)? split(/[-.]/,$protID0) : ($protID0,undef);
    $sthInsProt->execute($protID,$quantifParamIDs{'NUM_PEP_TOTAL'},$pepTotal{$protID0},undef);
    &promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $quantifiedModifID;
	foreach my $ratioPos (sort {$a <=> $b} keys %ratioPos2states) {
		next if !$validProt{$protID0}{$ratioPos};
        my (%pepDistinct,%pepUsed);
        foreach my $state (@{$ratioPos2states{$ratioPos}}){
            foreach my $pep (keys %{$pepList{$protID0}{$state}}){
                $pepUsed{$pep}=1;
                my ($pepSeq,$pepCharge)=split(/_/,$pep);
                $pepDistinct{$pepSeq}=1;
            }
        }
        $sthInsProt->execute($protID,$quantifParamIDs{'NUM_PEP_USED'},scalar keys %pepUsed,$ratioPos);
        &promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $quantifiedModifID;
        $sthInsProt->execute($protID,$quantifParamIDs{'DIST_PEP_USED'},scalar keys %pepDistinct,$ratioPos);
        &promsQuantif::insertProteinModResidues(\@modResidues,$sthInsProtRes,$modResiduesID{$protID},$dbh->last_insert_id(undef,undef,'PROTEIN_QUANTIFICATION','ID_PROT_QUANTIF')) if $quantifiedModifID;
    }
}
$sthInsProt->finish;
$dbh->commit;


###> Quantification is finished.
$dbh->do("UPDATE QUANTIFICATION SET STATUS=1 WHERE ID_QUANTIFICATION=$quantifID");
$dbh->commit;

$dbh->disconnect;

open(FILESTAT,">>$fileStat");
print FILESTAT "Ended ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
close(FILESTAT);

###> Moving final files
mkdir "$promsPath{quantification}/project_$projectID/quanti_$quantifID" unless -e "$promsPath{quantification}/project_$projectID/quanti_$quantifID";
dirmove($runDir,"$promsPath{quantification}/project_$projectID/quanti_$quantifID");
sleep 2;
unlink $fileStat;


####>Revision history<####
# 1.3.0 Now runs on cluster itself: R is launcher by system command (PP 17/09/18)
# 1.2.1 Added 2 wait loops for file sync after job done on cluster & new maxMem estimation (PP 29/08/18)
# 1.2.0 Updated DB insertion to runXICProtQuantification level (PP 19/07/18)
# 1.1.8 Minor modif to allow non-proteotypic quantification with MSstats (MLP 19/04/18)
# 1.1.7 Modif on PBSerror management ($clusterInfo{'checkError'}) (MLP 16/04/18)
# 1.1.6 Modif to store results files columns's names in hash table (MLP 09/02/18)
# 1.1.5 Minor modif (MLP 24/01/18)
# 1.1.4 Updated to use multiple CPUs if available (PP 24/01/18)
# 1.1.3 Minor modif on missing values (MLP 16/01/18)
# 1.1.2 Minor modif (MLP 12/01/18)
# 1.1.1 Minor modif on dbh connection (MLP 12/01/18)
# 1.1.0 Change in design processing to handle technical replicates (PP 28/12/17)
# 1.0.11 Keep ratio from ResultsDAProt.txt even if no p-value (PP 26/12/17)
# 1.0.10 Minor modif (MLP 15/12/17)
# 1.0.9 Minor modif to launch job on CentOS cluster (MLP 13/12/17)
# 1.0.8 Add modifications based quantification (MLP 15/11/2017)
# 1.0.7 Update for OpenSwath data (MLP 19/09/17)
# 1.0.6 Minor modification (MLP 26/01/17)
# 1.0.5 Add corrections for peptides count (MLP 08/12/2016)
# 1.0.4 Add &promsMod::cleanParameters verification (MLP 25/10/2016)
# 1.0.3 Protein filtering & user-defined globalStandards protein list (PP 19/08/16)
# 1.0.2 Creating globalStandards.txt for Glabal Standards normalization (MLP 12/08/2016)
# 1.0.1 Modification of Rscripts path (MLP 11/08/2016)
# 1.0.0 Created (MLP 20/07/2016)