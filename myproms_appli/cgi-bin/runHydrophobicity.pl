#!/usr/local/bin/perl -w

################################################################################
# runHydrophobicity.pl       1.0.0                                             #
# Authors: V. Laigle (Institut Curie)                                          #
# Contact: myproms@curie.fr                                                    #
# Launch Hydrophobicity computation and store data                             #
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

### myProMS software versions ($MYPROMS_HYDRO_VERSION) for Hydrophobicity Index computation (based on R script & data preprocessing in runHydrophobicity.pl)
# 1.0: based on SSRCalc algorithm 1st version (cf. script hydrophobicity.R for more info)

$| = 1;
use strict;
use warnings;
use POSIX qw(strftime);
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use promsQuantif;

#####################
### Configuration ###
#####################
my $MYPROMS_HYDRO_VERSION = '1.0';
my %promsPath = &promsConfig::getServerInfo('no_user');
my %clusterInfo = &promsConfig::getClusterInfo;
my $dbh = &promsConfig::dbConnect('no_user');

############################
### Parameters and files ###
############################
my ($quantifID, $quantifDir) = @ARGV;
my $projectID     = &promsMod::getProjectID($dbh, $quantifID, 'QUANTIFICATION');

my $workDir       = "$promsPath{tmp}/quantification/$quantifDir";
mkdir $workDir unless (-d $workDir);
my $infoFile      = "$workDir/hydro_info.txt";
my $dataFile      = "$workDir/hydro_data.txt";
my $resHydroFile  = "$workDir/hydro_results.txt";
my $resLinModFile = "$workDir/lin_mod_results.txt";
my $graphFile     = "$workDir/hydro_graph";
my $unmodifProps  = "$workDir/unmodified_proportions.txt";
my $logFile       = "$workDir/status_$quantifID.out";
my $errorFile     = "$workDir/Rerror.txt";

######################################################
### Fetch parameters and peptides data in analyses ###
######################################################
open(LOG, ">>$logFile");
print LOG "Started ", strftime("%H:%M:%S %d/%m/%Y", localtime), "\n";
print LOG "1/5 Fetching info and peptides data for hydrophobicity computation\n";
close LOG;

$dbh->do("UPDATE QUANTIFICATION SET STATUS = 0, QUANTIF_ANNOT = REPLACE(QUANTIF_ANNOT, 'SOFTWARE=myProMS','SOFTWARE=myProMS;$MYPROMS_HYDRO_VERSION') WHERE ID_QUANTIFICATION = $quantifID");
$dbh->commit;

my (%anaPepQuantifs, $rtOffset);
my ($quantifAnnot) = $dbh->selectrow_array("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION = $quantifID");
foreach my $quantifInfo (split('::', $quantifAnnot)) {
    next if ($quantifInfo =~ /LABEL=|SOFTWARE=/);
    if ($quantifInfo =~ /^ANA_PEPQUANTIF=/) {
        (my $anaPepQuantifStrg = $quantifInfo) =~ s/ANA_PEPQUANTIF=//;
        $anaPepQuantifStrg =~ s/#//g;  # Remove ID tags
        foreach my $anaPepQuantif (split(';', $anaPepQuantifStrg)) {
            my ($anaID, $pepQuantifID) = split(',', $anaPepQuantif);
            $anaPepQuantifs{$anaID} = $pepQuantifID;
        }
    } elsif ($quantifInfo =~ /^RT_OFFSET=/) {
        ($rtOffset = $quantifInfo) =~ s/RT_OFFSET=//;
    }
}

my (%anaInfo, %pepInfo, %modifiedPeps, %isGhost, %pepRTs);
my $sthGetAnaInfo = $dbh->prepare("SELECT A.NAME, S.ID_SAMPLE, S.NAME, A.VALID_STATUS 
    FROM ANALYSIS A 
    INNER JOIN SAMPLE S ON A.ID_SAMPLE = S.ID_SAMPLE 
    WHERE ID_ANALYSIS = ?"
);
foreach my $anaID (keys %anaPepQuantifs) {
    $sthGetAnaInfo->execute($anaID);
    my ($anaName, $sampleID, $sampleName, $valStatus) = $sthGetAnaInfo->fetchrow_array;
    $anaInfo{$anaID}{'anaName'}    = $anaName;
    $anaInfo{$anaID}{'sampleID'}   = $sampleID;
    $anaInfo{$anaID}{'sampleName'} = $sampleName;
    $anaInfo{$anaID}{'valStatus'}  = $valStatus;

    my $maxRT = 0;
    if ($valStatus >= 1) {
        $maxRT = &getAnaPepRT($anaID, \%pepInfo, \%modifiedPeps, \%isGhost);
    } else {
        $maxRT = &getQueryPepRT($anaID, \%pepInfo, \%modifiedPeps, \%isGhost);
    }

    if ($anaPepQuantifs{$anaID} eq "no_alignment") {
        if (exists $isGhost{$anaID}) {
            delete @{$pepInfo{$anaID}}{keys %{$isGhost{$anaID}}};
            delete @{$modifiedPeps{$anaID}}{keys %{$isGhost{$anaID}}};
        }
    } else {
        my ($labelType, $xicSoftCode) = &getLabelXicCode($anaPepQuantifs{$anaID});
        # Retrieve retention time according to quantif type when available in specific files
        if ($labelType eq 'FREE' && $xicSoftCode eq "MCQ") {
            &getMcqLfPepRT($anaPepQuantifs{$anaID}, \%pepRTs) unless ($pepRTs{$anaPepQuantifs{$anaID}});
        } elsif ($labelType eq 'SILAC' && $xicSoftCode eq "PD") {
            &getPdSilacPepRT($anaPepQuantifs{$anaID}, \%pepRTs) unless ($pepRTs{$anaPepQuantifs{$anaID}});
        }
        # Other cases -> RT is directly recorded in PEPTIDE table
        # $labelType = 'TMT', 'ITRAQ'
        # $xicSoftCode = "MQ", "OS", "SPC", "SKY"

        # Replace RTs from PEPTIDE table with reextraction RT
        if ($pepRTs{$anaPepQuantifs{$anaID}}) {
            foreach my $pepID (keys %{$pepInfo{$anaID}}) {
                $pepInfo{$anaID}{$pepID}{'RT'} = $pepRTs{$anaPepQuantifs{$anaID}}{$pepID} if ($pepRTs{$anaPepQuantifs{$anaID}}{$pepID});
            }
        }
    }

    if ($maxRT < 10) {  # Assume RT is in hours, need to convert it to minutes (especially for Spectronaut)
        foreach my $pepID (keys %{$pepInfo{$anaID}}) {
            $pepInfo{$anaID}{$pepID}{'RT'} = $pepInfo{$anaID}{$pepID}{'RT'} * 60;
        }        
    }
}
$sthGetAnaInfo->finish;

############################################
### Write data file expected by R script ###
############################################
open(LOG, ">>$logFile");
print LOG "2/5 Writing info for hydrophobicity R script\n";
close LOG;

open(DATA, ">$dataFile");
print DATA "Analysis_ID\tAnalysis\tSample_ID\tSample\tPeptide_ID\tSequence\tRT\tCharge\tValid_status\n";
foreach my $anaID (keys %anaInfo) {
    foreach my $pepID (keys %{$pepInfo{$anaID}}) {
        my $line = "$anaID\t$anaInfo{$anaID}{'anaName'}\t";
        $line .= "$anaInfo{$anaID}{'sampleID'}\t$anaInfo{$anaID}{'sampleName'}\t";
        $line .= "$pepID\t$pepInfo{$anaID}{$pepID}{'seq'}\t$pepInfo{$anaID}{$pepID}{'RT'}\t";
        $line .= "$pepInfo{$anaID}{$pepID}{'charge'}\t$pepInfo{$anaID}{$pepID}{'status'}\n";
        print DATA $line;
    }
}
close DATA;

open(PROPS, ">$unmodifProps");
print PROPS "Analysis_ID\tSample_ID\tUnmodified_pep_nb\tTotal_pep_nb\tUnmodified_proportion\n";
foreach my $anaID (keys %anaInfo) {    
    my $unmodifiedNb = ($pepInfo{$anaID}) ? scalar keys %{$pepInfo{$anaID}} : 0;
    my $modifNb      = ($modifiedPeps{$anaID})? scalar keys %{$modifiedPeps{$anaID}} : 0;
    my $totalNb      = $unmodifiedNb + $modifNb;
    my $unmodifProp  = $unmodifiedNb / $totalNb;

    print PROPS "$anaID\t$anaInfo{$anaID}{'sampleID'}\t$unmodifiedNb\t$totalNb\t$unmodifProp\n";
}
close PROPS;

####################################
### Write info file for R script ###
####################################
open(INFO, ">$infoFile");
print INFO "data_file\t$dataFile\n";
print INFO "output_hydro\t$resHydroFile\n";
print INFO "output_lin_mod\t$resLinModFile\n";
print INFO "graph_file\t$graphFile\n";
print INFO "rt_offset\t$rtOffset\n";
print INFO "plot_as\tpng\n";
print INFO "unmodified_props\t$unmodifProps\n";

close INFO;

#########################################
### Launch Hydrophobicity computation ###
#########################################
open(LOG, ">>$logFile");
print LOG "3/5 Waiting for Hydrophobicity R script to finish\n";
close LOG;
my $pathR = ($clusterInfo{'on'})? $clusterInfo{'path'}{'R'} : $promsPath{'R'};
system "cd $workDir; $pathR/R CMD BATCH --no-save --no-restore '--args $workDir' $promsPath{R_scripts}/hydrophobicity.R 2>> $errorFile";
&checkRerror("$workDir/hydrophobicity.Rout");

# Second check for computation success
unless (-s "$resHydroFile" && -s "$resLinModFile") {
    $dbh->do("UPDATE QUANTIFICATION SET STATUS = -1, UPDATE_DATE = NOW() WHERE ID_QUANTIFICATION = $quantifID"); # Failed
    $dbh->commit;
    $dbh->disconnect;
    die "An error occured during hydrophobicity computation, it failed to generate results files or the files are empty.";
}


####################################
### Move results to final folder ###
####################################
open(LOG, ">>$logFile");
print LOG "4/5 Moving Hydrophobicity output and storing results\n";
close LOG;
unless (-d "$promsPath{'quantification'}/project_$projectID/quanti_$quantifID") {
    unless (-d "$promsPath{'quantification'}/project_$projectID") {
        mkdir "$promsPath{'quantification'}" unless (-d "$promsPath{'quantification'}");
        mkdir "$promsPath{'quantification'}/project_$projectID";
    }
    mkdir "$promsPath{'quantification'}/project_$projectID/quanti_$quantifID";
}
# Move result files to final directory
system "mv $workDir/*.txt $promsPath{'quantification'}/project_$projectID/quanti_$quantifID/.";
system "mv $workDir/*.png $promsPath{'quantification'}/project_$projectID/quanti_$quantifID/.";
system "mv $workDir/*.Rout $promsPath{'quantification'}/project_$projectID/quanti_$quantifID/.";
unlink $workDir;


# Update DB after success of hydrophobicity computation
$dbh->do("UPDATE QUANTIFICATION SET STATUS = 1, UPDATE_DATE = NOW() WHERE ID_QUANTIFICATION = $quantifID");
$dbh->commit;
$dbh->disconnect;

open(LOG, ">>$logFile");
print LOG "5/5 Ended Hydrophobicity computation\n";
close LOG;
### End main computation, start subroutines ###


sub getAnaPepRT {  # Globals: %promsPath, $dbh, $projectID
    my ($anaID, $refPepInfo, $refModifiedPeps, $refIsGhost) = @_;

    my $maxRT = 0;
    my $sthPep = $dbh->prepare("SELECT P.ID_PEPTIDE, P.PEP_SEQ, P.CHARGE, P.ELUTION_TIME, P.VALID_STATUS,
        GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),
        GROUP_CONCAT(DISTINCT M.UNIMOD_ACC, ';', IS_LABEL ORDER BY PM.ID_MODIFICATION SEPARATOR '&')
        FROM PEPTIDE P 
        LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE = PM.ID_PEPTIDE
        LEFT JOIN MODIFICATION M ON PM.ID_MODIFICATION = M.ID_MODIFICATION
        WHERE P.ID_ANALYSIS = $anaID
        GROUP BY P.ID_PEPTIDE"
    );
    $sthPep->execute;
    while (my ($pepID, $pepSeq, $charge, $rtSc, $pepStatus, $modifs, $unimodsLabels) = $sthPep->fetchrow_array) {
        $refIsGhost->{$anaID}{$pepID} = 1 unless ($pepStatus);

        # Hydrophobicity not computed for modified peptides, except Carbamidomethyl Cysteine modifs (or labels)
        my $okModif = &checkLabelOrCarba($modifs, $unimodsLabels);
        if ($okModif) {
            my $pepRT = (!$rtSc) ? 0 : 
                            ($rtSc =~ /et([\d\.-]+)/) ? $1 : 
                            ($rtSc =~ /rt([\d\.-]+)/) ? $1 : 
                            ($rtSc =~ /^([\d\.-]+)/) ? $1 : 0;
            $refPepInfo->{$anaID}{$pepID}{'RT'}     = $pepRT;
            $refPepInfo->{$anaID}{$pepID}{'seq'}    = $pepSeq;
            $refPepInfo->{$anaID}{$pepID}{'status'} = $pepStatus;
            $refPepInfo->{$anaID}{$pepID}{'charge'} = $charge;
            # $refPepInfo->{$anaID}{$pepID}{'modifs'} = ($modifs)? $modifs : "";
            $maxRT = $pepRT if ($pepRT > $maxRT);
        } else {
            $refModifiedPeps->{$anaID}{$pepID} = 1;
        }
    }
    $sthPep->finish;

    return($maxRT);
}


sub getQueryPepRT {
    my ($anaID, $refPepInfo, $refModifiedPeps, $refIsGhost) = @_;
    my $maxRT = 0;
    my $sthPepQuery;  # = $dbh->prepare("SELECT * FROM ... ");

    return($maxRT);

}


sub getLabelXicCode {  # Globals: %promsPath, $dbh, $projectID
    my ($pepQuantifID) = @_;
    my ($labelType, $xicSoftCode);
    
    my ($pepQuantifAnnot) = $dbh->selectrow_array("SELECT QUANTIF_ANNOT FROM QUANTIFICATION Q WHERE Q.ID_QUANTIFICATION = $pepQuantifID");
	($labelType) = ($pepQuantifAnnot =~ /LABEL=([^:]+)/);
	$labelType = uc($labelType);
    if ($pepQuantifAnnot =~ /SOFTWARE=([^:|;]+);?([\d*\.?]*)/) {
        $xicSoftCode = $1;
    } elsif ($pepQuantifAnnot =~ /EXTRACTION_ALGO=/) {
    	$xicSoftCode = 'MCQ';
    } else {
    	my ($fileFormat) = $dbh->selectrow_array("SELECT FILE_FORMAT FROM ANALYSIS A, ANA_QUANTIFICATION Q WHERE A.ID_ANALYSIS = Q.ID_ANALYSIS AND ID_QUANTIFICATION = $pepQuantifID");
    	$xicSoftCode = ($fileFormat =~ /\.PDM\Z/)? 'PD' : 
                       ($fileFormat =~ /^MASCOT/)? 'MAS' : 
                       ($fileFormat =~ /PARAGON/)? 'PAR' : '?';
    }
    return ($labelType, $xicSoftCode);
}


sub getMcqLfPepRT {  # Globals: %promsPath, $dbh, $projectID
    my ($pepQuantifID, $refPepRTs) = @_;
    my $pepFile = "$promsPath{'quantification'}/project_$projectID/quanti_$pepQuantifID/peptide_quantification.txt";
    my ($rtParamID) = $dbh->selectrow_array("SELECT QP.ID_QUANTIF_PARAMETER
        FROM QUANTIFICATION Q
        INNER JOIN QUANTIFICATION_METHOD QM ON Q.ID_QUANTIFICATION_METHOD = QM.ID_QUANTIFICATION_METHOD
        INNER JOIN QUANTIFICATION_PARAMETER QP ON QM.ID_QUANTIFICATION_METHOD = QP.ID_QUANTIFICATION_METHOD
        WHERE Q.ID_QUANTIFICATION = $pepQuantifID AND QP.CODE = 'RT_APEX'"
    );
    open(PEP_QUANTIF, "$pepFile");
    while (my $line = <PEP_QUANTIF>) {
        next if ($line !~ /^$rtParamID/);  # Skip headers and lines with other quantif param
        chomp $line;
        my @values = split(/\t/, $line);
        my $pepID = $values[1];
        my $pepRT = $values[2] / 60;  # RT in seconds in MCQ output
        $refPepRTs->{$pepQuantifID}{$pepID} = $pepRT;
    }
    close PEP_QUANTIF;
}


sub getPdSilacPepRT {  # Globals: %promsPath, $dbh, $projectID
    my ($pepQuantifID, $refPepRTs) = @_;
    my @pepFiles = glob("$promsPath{'quantification'}/project_$projectID/quanti_$pepQuantifID/peptide_quantification*.txt");
    my ($rtParamID) = $dbh->selectrow_array("SELECT QP.ID_QUANTIF_PARAMETER
        FROM QUANTIFICATION Q
        INNER JOIN QUANTIFICATION_METHOD QM ON Q.ID_QUANTIFICATION_METHOD = QM.ID_QUANTIFICATION_METHOD
        INNER JOIN QUANTIFICATION_PARAMETER QP ON QM.ID_QUANTIFICATION_METHOD = QP.ID_QUANTIFICATION_METHOD
        WHERE Q.ID_QUANTIFICATION = $pepQuantifID AND QP.CODE = 'RT_APEX'"
    );
    foreach my $pepFile (@pepFiles) {
        open(PEP_QUANTIF, "$pepFile");
        while (my $line = <PEP_QUANTIF>) {
            next if ($line !~ /^$rtParamID/);  # Skip headers and lines with other quantif param
            chomp $line;
            my @values = split(/\t/, $line);
            my $pepID = $values[1];
            my $pepRT = $values[2];
            $refPepRTs->{$pepQuantifID}{$pepID} = $pepRT;
        }
        close PEP_QUANTIF;
    }
}


sub checkLabelOrCarba {
    my ($modifs, $unimodsLabels) = @_;
    if ($modifs) {
        if (!$unimodsLabels) {
            return 0;
        } else {
            my @modifsList = split(/&/, $modifs);
            my @unimodsLabelsList = split(/&/, $unimodsLabels);
            for my $i (0..$#modifsList) {
                my ($unimod, $label) = split(/;/, $unimodsLabelsList[$i]);
                if (!$unimod) {
                    return 0;
                } elsif ($unimod ne "4" && $label ne "1") {  # unimod 4 = carbamidomethyl or label
                    return 0;
                }
            }
        }
    }
    return 1;
}


sub checkRerror {  # Globals: $dbh, $quantifID
    my ($RoutFile) = @_;
    my $hydroErrorStrg;
    if (-e $RoutFile) {
        my $RoutStrg = `tail -4 $RoutFile`;
        unless ($RoutStrg =~ /proc\.time\(\)/) {
            $RoutStrg = `tail -20 $RoutFile`;
            $hydroErrorStrg = "Hydrophobicity computation has generated the following error:\n";
            my $inError = 0;
            foreach my $line (split(/\n/, $RoutStrg)) {
                next if (!$inError && $line !~ /^Error in/ && $line !~ /^Error:/);  # skip everything before error
                $inError = 1;
                $hydroErrorStrg .= "\n$line";
            }
        }
    } else {
        $hydroErrorStrg = "The .Rout file $RoutFile was not even written ! Check that the corresponding R script has actually run.\n";
    }
    if ($hydroErrorStrg) {
        $dbh->do("UPDATE QUANTIFICATION SET STATUS = -1, UPDATE_DATE = NOW() WHERE ID_QUANTIFICATION = $quantifID");
        $dbh->commit;
        $dbh->disconnect;
        die $hydroErrorStrg;
    }
}

####>Revision history<####
# 1.0.0 Created (VL 16/02/21)
