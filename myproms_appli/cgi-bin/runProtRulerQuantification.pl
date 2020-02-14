#!/usr/local/bin/perl -w

################################################################################
# runProtRulerQuantification.pl      1.0.2                                     #
# Authors: V. Laigle (Institut Curie)                                          #
# Contact: myproms@curie.fr                                                    #
################################################################################
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

$| = 1;
use promsConfig;
use promsMod;
use promsQuantif;
use strict;
use MIME::Base64;
use POSIX qw(strftime);
use File::Copy::Recursive qw(dirmove dircopy);


################################> TO DO <#######################################

# if correction factor, fetch data for it and add it to matrix input file

################################################################################

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo('no_user');
my %cluster=&promsConfig::getClusterInfo;  # ('debian'); # default is 'centos'

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

###> Get quantification parameters and other info from quantif_info.txt <###
my %quantifParameters=&promsQuantif::getQuantificationParameters("$quantifDir/quantif_info.txt");

my $runDir="$quantifDir/quantif_$quantifID";
my $dataDir="$runDir/data";
my $resultDir="$runDir/results";
my $graphDir="$resultDir/graph";
mkdir $runDir unless -e $runDir;
mkdir $dataDir unless -e $dataDir;
mkdir $resultDir unless -e $resultDir;
mkdir $graphDir unless -e $graphDir;

#########################################################################
### Insert input and output data for python script in parameters hash ###
### and keep consistent with the other value-types of the hash,       ###
### i.e. put every parameter in its own array (even if single value)  ###
#########################################################################
my @inputMatrix=("$dataDir/matrix_in.tsv",);
my @outFile=("$resultDir/matrix_out.tsv",);
my @groupsFile;
my @intensitiesIdx;
my %quantifAnnot;  # Will be filled gradually with parameters features

# Mapping column indexes to corresponding features to create input matrix
# and parse output from python script
my %finalColumnIndexes=(
	'Protein_ID'=> 0,
	'Protein_ACC'=> 1,  # Uniprot Accession Numbers (IDENTIFIERs)
	'Molecular_Weight'=> 2
);

#####> Mapping quantif parameter code to output nb in Proteomic Ruler <#####
my %quantifParamToOutput=('COPY_NB' => 0,           # Copy number
						  'CONCENTRATION' => 1,     # Concentration [nM]
						  'MASS_PER_CELL' => 2,     # Mass per cell [pg]
						  'MASS_ABUNDANCE' => 3,    # Abundance (mass/total mass) [*10^-6]
						  'MOL_ABUNDANCE' => 4,     # Abundance (molecules/total molecules) [*10^-6]
						  'COPY_RANK' => 5,         # Copy number rank
						  'REL_COPY_RANK' => 6      # Relative copy number rank
						 );

my %quantifParams;
my ($quantifMethodID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID");
my $sthQuantifParam=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,CODE FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=$quantifMethodID");
$sthQuantifParam->execute;
while (my ($paramID,$paramCode)=$sthQuantifParam->fetchrow_array) {
	$quantifParams{$paramID}=$paramCode;
}
$sthQuantifParam->finish;

my @quantifParamIDs=split(';', $quantifParameters{'R'}{'output_params'}[0]);
my @outputCodes;
foreach my $paramID (@quantifParamIDs){
	push @outputCodes, $quantifParamToOutput{$quantifParams{$paramID}};
}
my $output=join(';', @outputCodes);  # [0;3;5]
my @outputField=($output,);


# Resume mapping column indexes to corresponding features
my $columnIdx;
if ($quantifParameters{'R'}{'correction_factor'}[0] eq 'None') {
	$finalColumnIndexes{'Correction_factor'}='None';
	$quantifAnnot{'CORRECTION'}="None";
	$columnIdx=3;
} else {
	$finalColumnIndexes{'Correction_factor'}=3;
	$quantifAnnot{'CORRECTION'}=($quantifParameters{'R'}{'correction_factor'}[0]=~s/\s/_/g);
	$columnIdx=4;
}

my @selectedQuantifNames=split(';',$quantifParameters{'R'}{'quantif_names'}[0]);
my @sampleNames;
foreach my $parentQuantifName (@selectedQuantifNames){
	my $sampleName=($parentQuantifName=~/^(?:(?:LFQ\s)?[Ii]ntensity\s)?(.*?)(?:\s\[Intensities\])?$/)? $1 : $parentQuantifName;
	push @sampleNames, $sampleName;
	push @intensitiesIdx, $columnIdx;
	$finalColumnIndexes{$sampleName}=$columnIdx;
	$columnIdx++;
}
my $intensitiesStrg=join(';', @intensitiesIdx);
my @intensities=($intensitiesStrg,);

# Finished mapping indexes for input, start mapping output
# and fill information for targetPos
my %targetPos;
my @selectedQuantifications=split(';', $quantifParameters{'R'}{'quantif_ids'}[0]);
my @selectedStates;
foreach my $selectedState (split(';', $quantifParameters{'DB'}{'QUANTIF_STATES'}[0])) {
	my ($parQuantifID, $parStateID, $parPos) = split(',', $selectedState);
	push @selectedStates, [$parQuantifID, $parStateID, $parPos];
}
my %selectedQuantifIds;
for my $quantification (@selectedQuantifications){
	my ($parentQuantifId, $parentQuantifPos)=split('_', $quantification);
	$selectedQuantifIds{$parentQuantifId} = 1;
}
if ($quantifParameters{'R'}{'averaging_mode'}[0]==3){  # If 3, one quantification for all samples, no need for targetPos 
	foreach my $qParam (@quantifParamIDs){
		my $quantifIdParam="$quantifID\_$qParam";  # quantifID of the running quantif, not of the parent quantifs because only one quantif for all samples
		$finalColumnIndexes{$quantifIdParam}=$columnIdx;
		$targetPos{$quantifIdParam}=0;
		$columnIdx++;
	}
	foreach my $parent (@selectedStates){
		if ($quantifAnnot{'PARENT_Q'}){
			$quantifAnnot{'PARENT_Q'}.=";#$parent->[0],#$parent->[1],$parent->[2]:1";
		} else {
			$quantifAnnot{'PARENT_Q'}="#$parent->[0],#$parent->[1],$parent->[2]:1";
		}
	}
} else {
	my $posCount=0;
	foreach my $parent (@selectedStates){
		$posCount++;
		if ($quantifAnnot{'PARENT_Q'}){
			$quantifAnnot{'PARENT_Q'}.=";#$parent->[0],#$parent->[1],$parent->[2]:$posCount";
		} else {
			$quantifAnnot{'PARENT_Q'}="#$parent->[0],#$parent->[1],$parent->[2]:$posCount";
		}
		foreach my $qParam (@quantifParamIDs){
			my $quantifIdPosParam="$parent->[0]\_$parent->[2]\_$qParam";  # parent->[0] = quantifID of the parent quantif
			$finalColumnIndexes{$quantifIdPosParam}=$columnIdx;
			$targetPos{$quantifIdPosParam}=$posCount;
			$columnIdx++;
		}
	}
}

# Create groups' file if necessary and complete quantifAnnot
my @groups;
if ($quantifParameters{'R'}{'groups'}[0] eq 'None') {
	@groupsFile=("None",);
	$quantifAnnot{'GROUPING'}="None";
} else {
	@groupsFile=("$dataDir/groups.tsv",);
	open(GROUPS_FILE, ">$groupsFile[0]") or die "Couldn't open groups file: $!";
	@groups=split(';',$quantifParameters{'R'}{'groups'}[0]);
	for my $i (0..$#sampleNames){
		print GROUPS_FILE "$sampleNames[$i]\t$groups[$i]\n";
		my $groupNb=$groups[$i];
		$groupNb=~s/\D//g;
		if ($quantifAnnot{'GROUPING'}){
			$quantifAnnot{'GROUPING'}.=";#$selectedStates[$i]->[0],#$selectedStates[$i]->[1],$selectedStates[$i]->[2]:$groupNb";
		} else {
			$quantifAnnot{'GROUPING'}="#$selectedStates[$i]->[0],#$selectedStates[$i]->[1],$selectedStates[$i]->[2]:$groupNb";
		}
	}
	close GROUPS_FILE;
}

if ($quantifParameters{'R'}{'averaging_mode'}[0]==0){
	$quantifAnnot{'AVG_MODE'}="ALL_SEPARATE";
} elsif ($quantifParameters{'R'}{'averaging_mode'}[0]==1){
	$quantifAnnot{'AVG_MODE'}="ALL_SAME_NORM";
} elsif ($quantifParameters{'R'}{'averaging_mode'}[0]==2){
	$quantifAnnot{'AVG_MODE'}="GROUPS_SAME_NORM";
} elsif ($quantifParameters{'R'}{'averaging_mode'}[0]==3){
	$quantifAnnot{'AVG_MODE'}="AVG_ALL";
}

# Prot_Ruler type (int):
# 0 = Total protein amount,
# 1 = Histone proteoimic ruler,
# 2 = Custom proteins Ruler
$quantifAnnot{'RULER_TYPE'}=$quantifParameters{'R'}{'prot_ruler_type'}[0]; 

# Intensity metric (MaxQuant): Intensity, LFQ
$quantifAnnot{'METRIC'}=$quantifParameters{'R'}{'int_metric'}[0];

# Numerical parameters:
# Total_prot_concentration;Total_prot_amount;Ploidy;custom_prot_quantities_1,2,3...
# 2 out of the last 3 are "None"
my $customProtQtyStrg=$quantifParameters{'R'}{'custom_prot_qty'}[0];
$customProtQtyStrg=~s/;/,/g; 
$quantifAnnot{'NUM_PARAMS'}=join(";", ($quantifParameters{'R'}{'protein_concentration'}[0],
									   $quantifParameters{'R'}{'total_protein_amount'}[0],
									   $quantifParameters{'R'}{'histone_proteomic_ruler'}[0],
									   $customProtQtyStrg)
								);

my $quantifAnnotStrg=join("::", ("RULER_TYPE=$quantifAnnot{'RULER_TYPE'}",
								 "AVG_MODE=$quantifAnnot{'AVG_MODE'}",
								 "GROUPING=$quantifAnnot{'GROUPING'}",
								 "PARENT_Q=$quantifAnnot{'PARENT_Q'}",
								 "METRIC=$quantifAnnot{'METRIC'}",
								 "CORRECTION=$quantifAnnot{'CORRECTION'}",
								 "NUM_PARAMS=$quantifAnnot{'NUM_PARAMS'}")
						 );

$dbh->do("UPDATE QUANTIFICATION SET QUANTIF_ANNOT=CONCAT(QUANTIF_ANNOT, '::$quantifAnnotStrg') WHERE ID_QUANTIFICATION=$quantifID");
$dbh->commit;


# Add missing parameters to quantifParameters hash to write parameters' file for python script
$quantifParameters{'R'}{'input_matrix'}=\@inputMatrix;
$quantifParameters{'R'}{'out_file'}=\@outFile;
$quantifParameters{'R'}{'groups_file'}=\@groupsFile;
$quantifParameters{'R'}{'protein_acc'}=[$finalColumnIndexes{'Protein_ACC'},];
$quantifParameters{'R'}{'molecular_weights'}=[$finalColumnIndexes{'Molecular_Weight'},];
$quantifParameters{'R'}{'correction_factor_idx'}=[$finalColumnIndexes{'Correction_factor'},];
$quantifParameters{'R'}{'intensities'}=\@intensities;
$quantifParameters{'R'}{'output'}=\@outputField;

&promsQuantif::writeQuantifParameterFiles($dataDir,$quantifParameters{'R'});  # $dataDir: full path to tmp quantif dir/data/


######################################
####>Protein lists for filtering <####
######################################
my ($protSelectionType,%selectExcludeProteins);
my ($refSelectedProteins,$refExcludedProteins);
my $sthList=$dbh->prepare("SELECT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=? ORDER BY ID_PROTEIN");

###>Selection/exclusion list
if ($quantifParameters{'DB'}{'PROTEINS'}) {
	($protSelectionType,my $listID)=@{$quantifParameters{'DB'}{'PROTEINS'}};
	$listID=~s/#//; # remove id flag
	$sthList->execute($listID);
	while (my ($protID)=$sthList->fetchrow_array) {
		$selectExcludeProteins{$protID}=1;
	}
	($refSelectedProteins,$refExcludedProteins)=($protSelectionType eq 'restrict')? (\%selectExcludeProteins,undef) : (undef,\%selectExcludeProteins);
}

$sthList->finish;

################################################
####> Fetching protein quantification data <####
################################################

# Parameters to fetch data from previous MaxQuant quantification
my $selQuantifFamily=$quantifParameters{'DB'}{'QUANTIF_FAMILY'}[0];
my $view='list';
my $numPepCode='NUM_PEP_USED';
my $selModifID=0;

my ($intMetricCode)=$dbh->selectrow_array("SELECT CODE FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIF_PARAMETER=$quantifParameters{'R'}{'int_metric'}[0]");

my (%quantifValues,%quantifInfo,%proteinInfo);
my %parameters=(QUANTIF_FAMILY=>$selQuantifFamily,
				VIEW=>$view,
				NUM_PEP_CODE=>$numPepCode,
				QUANTIF_LIST=>\@selectedQuantifications,
				SEL_MODIF_ID=>$selModifID,
				MEASURE=>$intMetricCode
			   );

&promsQuantif::fetchQuantificationData($dbh,\%parameters,\%quantifInfo,\%quantifValues,undef,\%proteinInfo,$refSelectedProteins,$refExcludedProteins);

#####> Replace proteins ALIAS with Uniprot ACC in fetched data <#####
my $sthGetUniprotAcc=$dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER MI JOIN PROTEIN P ON MI.ID_MASTER_PROTEIN=P.ID_MASTER_PROTEIN WHERE P.ID_PROTEIN=? AND MI.ID_IDENTIFIER=1 AND MI.RANK=1");
foreach my $protID (keys %proteinInfo){
	$sthGetUniprotAcc->execute($protID);
	my ($protIdentifier)=$sthGetUniprotAcc->fetchrow_array;
	$proteinInfo{$protID}[0]=$protIdentifier;
}
$sthGetUniprotAcc->finish;


###################################################################
###> Writing values from quantifValues and proteinInfo to file <###
###> Building the input file for proteomic ruler python script <###
###################################################################
# Protein_ID Uniprot_ACC Molecular_Weight [Correction_factor] Intensity_sample_1 [Intensity_sample_2 ...]

my $header="Protein_IDs\tUniprot_ACC\tMolecular_Weights";
#if ($quantifParameters{'R'}{'detectability_correction'}[0]){
#	$header.="\tCorrection_factor";
#}
foreach my $name (@selectedQuantifNames){
	$header.="\t$name";
}
open(MATRIX_FILE, ">$quantifParameters{'R'}{'input_matrix'}[0]") or die "Couldn't open input matrix file $!";
print MATRIX_FILE "$header\n";

my $line;
foreach my $protID (keys %quantifValues){
	$line="$protID\t@{$proteinInfo{$protID}}[0]\t@{$proteinInfo{$protID}}[1]";   
	foreach my $quantif (@selectedQuantifications){
		if ($quantifValues{$protID} && $quantifValues{$protID}{$quantif} && defined($quantifValues{$protID}{$quantif}{$intMetricCode})){
			$line.="\t$quantifValues{$protID}{$quantif}{$intMetricCode}";
		} else {
			$line.="\tN/A";
		}
	}
	print MATRIX_FILE "$line\n";
}
close MATRIX_FILE;

$dbh->disconnect;


###############################
###> Running Python Script <###
###############################
open(FILESTAT,">>$fileStat");
print FILESTAT "2/3 Running quantification\n";
close(FILESTAT);
my $pathPython=($cluster{'on'})? $cluster{'path'}{'python'} : $promsPath{'python'};
open(PYTHON_SCRIPT,">$runDir/proteomic_ruler.txt");
print PYTHON_SCRIPT qq
|
##############################################
# Launcher for quantification python scripts #
##############################################

filepath="$promsPath{python_scripts}/"
source(paste(filepath,"proteomic_ruler.py",sep=""))
|;
close PYTHON_SCRIPT;
my $pythonCommandString="export LANG=en_US.UTF-8; cd $runDir; " .
	"$pathPython/python3 $promsPath{'python_scripts'}/proteomic_ruler.py -c $dataDir/param_char.txt"; # -n $dataDir/param_num.txt";

my $pythonOutFile=$quantifParameters{'R'}{'out_file'}[0];

system $pythonCommandString;

sleep 3;


########################################
###> PARSING QUANTIFICATION RESULTS <###
########################################
open(FILESTAT,">>$fileStat");
print FILESTAT "3/3 Parsing results\n";
close FILESTAT;

$dbh=&promsConfig::dbConnect('no_user'); # reconnect


# INSERT STATEMENTS
my $sthInsProtQ=$dbh->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_QUANTIFICATION,ID_PROTEIN,ID_QUANTIF_PARAMETER,QUANTIF_VALUE,TARGET_POS) VALUES ($quantifID,?,?,?,?)") or die "Couldn't prepare statement: " . $dbh->errstr;
my $sthInsParentQ=$dbh->prepare("INSERT INTO PARENT_QUANTIFICATION (ID_PARENT_QUANTIFICATION,ID_QUANTIFICATION,PAR_FUNCTION) VALUES (?,$quantifID,NULL)") or die "Couldn't prepare statement: " . $dbh->errstr;
my $sthInsAnaQ=$dbh->prepare("INSERT INTO ANA_QUANTIFICATION (ID_QUANTIFICATION,ID_ANALYSIS) VALUES ($quantifID,?)") or die "Couldn't prepare statement: " . $dbh->errstr;

# STATEMENTS TO LINK QUANTI TO PARENT ANALYSIS
my $allParentQStrg=join(',', (keys %selectedQuantifIds));
my ($allAnaStrg)=$dbh->selectrow_array("SELECT GROUP_CONCAT(DISTINCT ID_ANALYSIS) FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION IN ($allParentQStrg)");

# Parse results file and Insert data into DB
if (-e "$pythonOutFile") {
	open(RESULTS, $pythonOutFile) or die "Couldn't open results file: $!";
	foreach my $line (<RESULTS>){
		if ($line=~/Protein_IDs\tUniprot_ACC\tMolecular_Weights/) {
			next
		} else {
			chomp $line;
			my @fields=split("\t", $line);
			my $proteinID=$fields[$finalColumnIndexes{'Protein_ID'}];
			# my $uniprotACC=$fields[$finalColumnIndexes{'Protein_ACC'}];
			# my $molWeight=$fields[$finalColumnIndexes{'Molecular_Weight'}];
			# if ($quantifParameters{'R'}{'correction_factor'}[0] ne 'None'){
			#	my $correctFactor=$fields[$finalColumnIndexes{'Correction_factor'}];
			#}
			if ($quantifParameters{'R'}{'averaging_mode'}[0]==3){  # If 3, one quantification for all samples
				foreach my $qParam (@quantifParamIDs){
					my $quantifIdParam="$quantifID"."_"."$qParam";  # quantifID of the running quantif, not of the parent quantifs (only one quantif for all samples)
					my $quantifValue=$fields[$finalColumnIndexes{$quantifIdParam}];
					$sthInsProtQ->execute($proteinID,$qParam,$quantifValue,$targetPos{$quantifIdParam});
				}
			} else {
				foreach my $parentQuantif (@selectedQuantifications){
					foreach my $qParam (@quantifParamIDs){
						my $quantifIdPosParam="$parentQuantif\_$qParam";  # quantifID of the parent quantif, 
						my $quantifValue=$fields[$finalColumnIndexes{$quantifIdPosParam}];
						$sthInsProtQ->execute($proteinID,$qParam,$quantifValue,$targetPos{$quantifIdPosParam});
					}
				}
			}
		}
	}
	close RESULTS;
	
	foreach my $parentQuantifId (sort {$a <=> $b} keys %selectedQuantifIds){
		$sthInsParentQ->execute($parentQuantifId);
	}
	
	foreach my $anaID (split(',', $allAnaStrg)){
		$sthInsAnaQ->execute($anaID);
	}
	
} else { # File not read (does not exist)
	$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID");
	$dbh->commit;
	$dbh->disconnect;
	die "Missing results file !\n";
}

# Clean handles
$sthInsProtQ->finish;
$sthInsParentQ->finish;
$sthInsAnaQ->finish;
$dbh->commit;

# Quantification is finished.
$dbh->do("UPDATE QUANTIFICATION SET STATUS=1 WHERE ID_QUANTIFICATION=$quantifID");
$dbh->commit;
$dbh->disconnect;

open(FILESTAT,">>$fileStat");
print FILESTAT "Ended ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
close(FILESTAT);

# Moving final files
mkdir "$promsPath{quantification}/project_$projectID/quanti_$quantifID" unless -e "$promsPath{quantification}/project_$projectID/quanti_$quantifID";
dirmove($runDir,"$promsPath{quantification}/project_$projectID/quanti_$quantifID");
sleep 2;
unlink $fileStat;


#####>Revision history<#####
# 1.0.2 [ENHANCEMENT] Adapt code to accept new intensity metrics (VL 24/09/19)
# 1.0.1 [CODE] Adjusted parameter list to match change in &promsQuantif::fetchQuantificationData (PP 30/08/19)
# 1.0.0 Created (VL 26/07/2019)
