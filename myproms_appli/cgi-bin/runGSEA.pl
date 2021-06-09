#!/usr/local/bin/perl -w

################################################################################
# runGSEA.pl       1.0.0                                                       #
# Authors: V. Laigle (Institut Curie)                                          #
# Contact: myproms@curie.fr                                                    #
# Launch GSEA computation and store data                                       #
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
use strict;
use warnings;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use promsQuantif;

#####################
### Configuration ###
#####################
my %promsPath = &promsConfig::getServerInfo('no_user');
my %clusterInfo = &promsConfig::getClusterInfo;
my $dbh = &promsConfig::dbConnect('no_user');

############################
### Parameters and files ###
############################
my ($gseaID, $userID, $tmpDir) = @ARGV;
my $gseaDir     = "$promsPath{'tmp'}/gsea/$tmpDir/gsea_$gseaID";  # Not the same gseaDir as in startGSEA
mkdir $gseaDir unless (-d $gseaDir);
my $infoFile    = "$gseaDir/gsea_info.txt";
my $dataFile    = "$gseaDir/gsea_data.txt";
my $resultsFile = "$gseaDir/gsea_results.txt";
my $gUsedFile   = "$gseaDir/gsea_genes_used.txt";
my $graphFiles  = "$gseaDir/gsea_graph";
my $lostProts   = "$gseaDir/prot_not_used.txt";
my $logFile     = "$promsPath{'tmp'}/gsea/$tmpDir/status_$gseaID.out";
my $errorFile   = "$promsPath{'tmp'}/gsea/$tmpDir/error_$gseaID.txt";
my $projectID   = &promsMod::getProjectID($dbh, $gseaID,'PATHWAY_ANALYSIS');

######################################
### Fetch info about GSEA analysis ###
######################################
open(LOG, ">>$logFile");
print LOG "1/5 Fetching GSEA info and quantification data\n";
close LOG;

my ($gseaName, $gseaParams) = $dbh->selectrow_array("SELECT NAME, PARAM_STRG FROM PATHWAY_ANALYSIS WHERE ID_PATHWAY_ANALYSIS = $gseaID");
my ($quantifID, $targetPos, $isRatio);
my ($geneSetsID, $speciesID);
my ($protSelType, $protListID);
my ($pepWeights, $weightByPep, $pvalCutoff, $pvalCorrection, $maxRelevantGS);
foreach my $param (split(';', $gseaParams)) {
    $param =~ s/#//g;  # Clean ID tags
    if ($param =~ /^quantification=(\d+)$/) {
        $quantifID = $1;
    } elsif ($param =~ /^targetPos=(.*)$/) {
        $targetPos = $1;
    } elsif ($param =~ /^isRatio=(\d)$/) {
        $isRatio = $1;
    } elsif ($param =~ /^species=(.+)$/) {  # ID_SPECIES but may also equal "Custom"
        $speciesID = $1;
    } elsif ($param =~ /^geneSetsDB=(\d+)$/) {
        $geneSetsID = $1;
    } elsif ($param =~ /^protSel=(.+)$/) {
        ($protSelType, $protListID) = split(',', $1);
    } elsif ($param =~ /^weightByPep=(.+)$/) {
        $pepWeights = $1;
        $weightByPep = ($pepWeights eq "none")? 0 : 1;
    } elsif ($param =~ /^pvalCutoff=(.+)$/) {
        $pvalCutoff = $1;
    } elsif ($param =~ /^pvalCorr=(.+)$/) {
        $pvalCorrection = $1;
    } elsif ($param =~ /^maxGS=(\d+)$/) {
        $maxRelevantGS = $1;
    }
}

# Info about gene sets (.gmt file, identifier type)
my ($identID, $identType, $gmtFile) = $dbh->selectrow_array("SELECT I.ID_IDENTIFIER, I.NAME, GS.GENESETS_FILE FROM GENESETS GS LEFT JOIN IDENTIFIER I ON GS.ID_IDENTIFIER = I.ID_IDENTIFIER WHERE GS.ID_GENESETS = $geneSetsID");
$gmtFile = "$promsPath{'genesets'}/$gmtFile";
$identType = "User-defined" unless ($identType);  # $identID and $identType can be NULL in DB

# Species name
my $species;
if ($speciesID eq "Custom") {
    $species = $speciesID;
} else {
    ($species) = $dbh->selectrow_array("SELECT SCIENTIFIC_NAME FROM SPECIES WHERE ID_SPECIES = $speciesID");
}

# Identifier mapping to match identifiers in gene set file
my %identMap;
my %otherSpeciesProts;
my $sthIdentMap;
if ($identID) {  # Gene sets with usual identifiers -> get identifier with MASTER_PROTEIN (if it exists)
    $sthIdentMap = $dbh->prepare("SELECT P.ID_PROTEIN, MI.VALUE, P.ORGANISM FROM PROTEIN P LEFT JOIN MASTERPROT_IDENTIFIER MI ON P.ID_MASTER_PROTEIN = MI.ID_MASTER_PROTEIN WHERE P.ID_PROJECT = $projectID AND MI.ID_IDENTIFIER = $identID AND MI.IDENT_RANK = 1");
} else {  # Gene sets with user-defined identifiers -> just take the identifier from PROTEIN
    $sthIdentMap = $dbh->prepare("SELECT ID_PROTEIN, IDENTIFIER, ORGANISM FROM PROTEIN WHERE ID_PROJECT = $projectID");
}
$sthIdentMap->execute;
if ($species eq "Custom") {  # No species = custom gene sets -> use all proteins
    while (my ($protID, $identifier, $organism) = $sthIdentMap->fetchrow_array) {
        $identMap{$protID} = $identifier;
    }
} else {  # Otherwise take only proteins from the same organism as the gene sets
    while (my ($protID, $identifier, $organism) = $sthIdentMap->fetchrow_array) {
        if ($organism) {
            if ($organism eq $species) {
                $identMap{$protID} = ($identifier)? $identifier : '';  # In case no identifier found
            } else {
                $otherSpeciesProts{$protID} = $organism;
            }
        } else {
            $otherSpeciesProts{$protID} = "no species";
        }
    }
}
$sthIdentMap->finish;

#################################
### Fetch quantification data ###
#################################
my (%parameters, %quantifInfo, %quantifValues);
my (%dispModifSites, %proteinInfo, %selectedProteins, %excludedProteins);
my $quantif = $quantifID.'_'.$targetPos;

# Protein List exclusion / restriction
if ($protListID) {
    my $refprotList = ($protSelType eq 'restrict')? \%selectedProteins : \%excludedProteins;
    my ($type) = $dbh->selectrow_array("SELECT LIST_TYPE FROM CATEGORY WHERE ID_CATEGORY = $protListID");
    if ($type eq 'SITE') {  # Not used yet (in case GSEA is used on sites later)
        &promsQuantif::fetchCustomList($dbh, $protListID, $refprotList);
    } else {
        &promsQuantif::fetchCustomList($dbh, $protListID, $refprotList, 1);
    }
}

# Peptides numbers to get
my %pepTypeCodes = (
    none     => "NUM_PEP_USED",
    all      => "NUM_PEP_USED",
    distinct => "DIST_PEP_USED",
    msms     => "NUM_TRUE_USED"
);

# Quantification parameters
if ($isRatio) {
    %parameters = (
        QUANTIF_FAMILY  => 'RATIO',
        VIEW            => 'log2',
        NUM_PEP_CODE    => $pepTypeCodes{$pepWeights},
        QUANTIF_LIST    => [$quantif]
    );
} else {
    die "GSEA on non-ratio quantification is not supported yet";
}

&promsQuantif::fetchQuantificationData($dbh, \%parameters, \%quantifInfo, \%quantifValues, \%dispModifSites, \%proteinInfo, \%selectedProteins, \%excludedProteins);
my $log2 = log(2);
foreach my $protID (keys %{$quantifValues{$quantif}}) {
    $quantifValues{$quantif}{$protID}{RATIO} = log($quantifValues{$quantif}{$protID}{RATIO}) / $log2;
}

############################################
### Write data file expected by R script ###
############################################
open(LOG, ">>$logFile");
print LOG "2/5 Writing info for GSEA R script\n";
close LOG;

my %lostProteins;
open(DATA, ">$dataFile");
print DATA "Protein_ID\tGene_name\tValue\tPeptides_nb\tPvalue\n";
foreach my $protID (keys %{$quantifValues{$quantif}}) {
    if ($identMap{$protID}) {
        if ($isRatio) {
            $quantifValues{$quantif}{$protID}{P_VALUE} = "NA" unless $quantifValues{$quantif}{$protID}{P_VALUE};
            print DATA "$protID\t$identMap{$protID}\t$quantifValues{$quantif}{$protID}{RATIO}\t$quantifValues{$quantif}{$protID}{$pepTypeCodes{$pepWeights}}\t$quantifValues{$quantif}{$protID}{P_VALUE}\n";
        } else {
            die "GSEA on non-ratio quantification is not supported yet";
        }
    } elsif (!$otherSpeciesProts{$protID}) {
        $lostProteins{$protID} = 1;
    }
}
close DATA;

open(LOST, ">$lostProts");
print LOST "ProteinID\tWhy not used\n";
foreach my $protID (keys %otherSpeciesProts) {
    print LOST "$protID\t$otherSpeciesProts{$protID}\n";
}
foreach my $protID (keys %lostProteins) {
    print LOST "$protID\tNo $identType found\n";
}
close LOST;

####################################
### Write info file for R script ###
####################################
open(INFO, ">$infoFile");
print INFO "data_file\t$dataFile\n";
print INFO "is_ratio\t$isRatio\n";
print INFO "species\t$species\n";
print INFO "gmt_file\t$gmtFile\n";
print INFO "pep_weights\t$weightByPep\n";
print INFO "pval_cutoff\t$pvalCutoff\n";
print INFO "pval_correction\t$pvalCorrection\n";
print INFO "max_relevant_gs\t$maxRelevantGS\n";
print INFO "output_file\t$resultsFile\n";
print INFO "output_graph\t$graphFiles\n";
print INFO "genes_used_file\t$gUsedFile\n";
close INFO;

###############################
### Launch GSEA computation ###
###############################
open(LOG, ">>$logFile");
print LOG "3/5 Waiting for GSEA R script to finish\n";
close LOG;
my $pathR = ($clusterInfo{'on'})? $clusterInfo{'path'}{'R'} : $promsPath{'R'};
system "cd $gseaDir; $pathR/R CMD BATCH --no-save --no-restore '--args $gseaDir' $promsPath{R_scripts}/GSEA.R 2>> $errorFile";
&checkRerror("$gseaDir/GSEA.Rout");

# Second check for computation success
unless (-s "$gseaDir/gsea_results.txt") {
    $dbh->do("UPDATE PATHWAY_ANALYSIS SET STATUS = -1, UPDATE_DATE = NOW() WHERE ID_PATHWAY_ANALYSIS = $gseaID"); # Failed
    $dbh->commit;
    $dbh->disconnect;
    die "An error occured during GSEA computation, it failed to generate $gseaDir/gsea_results.txt or the file is empty.";
}


####################################
### Move results to final folder ###
####################################
open(LOG, ">>$logFile");
print LOG "4/5 Moving GSEA output and storing results\n";
close LOG;
unless (-d "$promsPath{'gsea'}/project_$projectID/gsea_$gseaID") {
    unless (-d "$promsPath{'gsea'}/project_$projectID") {
        mkdir "$promsPath{'gsea'}" unless (-d "$promsPath{'gsea'}");
        mkdir "$promsPath{'gsea'}/project_$projectID";
    }
    mkdir "$promsPath{'gsea'}/project_$projectID/gsea_$gseaID";
}
# Move result files to final directory
system "mv $gseaDir/*.txt $promsPath{'gsea'}/project_$projectID/gsea_$gseaID/.";
system "mv $gseaDir/*.png $promsPath{'gsea'}/project_$projectID/gsea_$gseaID/.";
system "mv $gseaDir/*.Rout $promsPath{'gsea'}/project_$projectID/gsea_$gseaID/.";
unlink $gseaDir;


# Update DB after success of GSEA computation
$dbh->do("UPDATE PATHWAY_ANALYSIS SET STATUS = 1, UPDATE_DATE = NOW() WHERE ID_PATHWAY_ANALYSIS = $gseaID");
$dbh->commit;
$dbh->disconnect;

open(LOG, ">>$logFile");
print LOG "5/5 Ended GSEA computation\n";
close LOG;
### End main computation, start subroutines ######


sub checkRerror {  # Globals: $dbh, $gseaID
    my ($RoutFile) = @_;
    my $gseaErrorStrg;
    if (-e $RoutFile) {
        my $RoutStrg = `tail -4 $RoutFile`;
        unless ($RoutStrg =~ /proc\.time\(\)/) {
            $RoutStrg = `tail -20 $RoutFile`;
            $gseaErrorStrg = "GSEA computation has generated the following error:\n";
            my $inError = 0;
            foreach my $line (split(/\n/, $RoutStrg)) {
                next if (!$inError && $line !~ /^Error in/ && $line !~ /^Error:/);  # skip everything before error
                $inError = 1;
                $gseaErrorStrg .= "\n$line";
            }
        }
    } else {
        $gseaErrorStrg = "The .Rout file $RoutFile was not even written ! Check that the corresponding R script has actually run.\n";
    }
    if ($gseaErrorStrg) {
        $dbh->do("UPDATE PATHWAY_ANALYSIS SET STATUS = -1, UPDATE_DATE = NOW() WHERE ID_PATHWAY_ANALYSIS = $gseaID");
        $dbh->commit;
        $dbh->disconnect;
        die $gseaErrorStrg;
    }
}

####>Revision history<####
# 1.0.0 Created (VL 21/10/20)
