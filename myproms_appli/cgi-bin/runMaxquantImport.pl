#!/usr/local/bin/perl -w

################################################################################
# runMaxquantImport.pl       2.0.8                                             #
# Component of site myProMS Web Server                                         #
# Authors: P. Poullet, G. Arras, S. Liva (Institut Curie)                      #
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

use strict;
use POSIX qw(strftime);  # Core module
#use IO::Uncompress::Gunzip qw(gunzip);
#use Cwd;
use XML::Simple;
use File::Path qw(rmtree); # remove_tree
use File::Copy qw(copy move); # Core module
use File::Spec::Functions qw(splitpath); # Core module
use promsConfig;
use promsMod;
use promsQuantif;
# exit;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo('no_user');
my $MAX_DB=3;
my ($experimentID,$jobID,$onCluster)=@ARGV;
$experimentID=&promsMod::cleanNumericalParameters($experimentID);

my $currentQuantifDir="$promsPath{tmp}/quantification/current";
my $tmpFilesDir="$promsPath{tmp}/quantification/$jobID";
my $fileStat="$tmpFilesDir/status_$experimentID.out";
my $timeStamp=strftime("%H:%M:%S %d/%m/%Y",localtime);
&addToFileStat("Started $timeStamp\n");
my $detailFile="maxquant_import_$experimentID.log";
open(LOG,">$tmpFilesDir/$detailFile") || die "Error while opening $tmpFilesDir/$detailFile";
print LOG "Started $timeStamp";

my ($userID,$projectID,@databankIDs,$excludeCON,$matchGroupType,$importProtQuantif,$archiveFile); # $paramGrIdx,
my $contaminantDB=0; # default
open (INFO,"$tmpFilesDir/quantif_info.txt");
while(<INFO>) {
	chomp;
	my ($paramN,$paramV)=split('=',$_);
	if ($paramN eq 'USER') {$userID=$paramV;}
	elsif ($paramN eq 'ID_PROJECT') {$projectID=$paramV;}
	elsif ($paramN eq 'DATABANKS') {@databankIDs=split(' ',$paramV);}
	elsif ($paramN eq 'EXCLUDE_CON') {$excludeCON=$paramV;}
	elsif ($paramN eq 'ID_CONT_DB') {$contaminantDB=$paramV;}
	elsif ($paramN eq 'MG_TYPE') {$matchGroupType=$paramV;}
	#elsif ($paramN eq 'PARAM_GR') {$paramGrIdx=$paramV;}
	elsif ($paramN eq 'PROT_QUANTIF') {$importProtQuantif=$paramV;}
	elsif ($paramN eq 'ARCHIVE') {$archiveFile=$paramV;}
}
close INFO;
push @databankIDs,$contaminantDB if $contaminantDB;

my $importPepQuantif=1;
my $numSteps=($importPepQuantif)? 11 : 10;

####>Data Files<####
my $evidenceFile='evidence.txt'; # Contains all peptide identified and peptide quantitation
my $peptideFile='peptides.txt'; # Contains all peptide identified and peptide quantitation
my $proteinGroupsFile='proteinGroups.txt'; # Contains all quantitation information at protein level (LFQ, iBAQ, Intensity,...)
my $msmsFile='msms.txt'; # Contains m/z and intensities information
my $summaryFile='summary.txt'; # Contains some redundant information found in mqpar.xml (varmod, label type, experimental design)
my $parametersFile='parameters.txt'; # Contains some redundant information found in mqpar.xml (fixmod, peptide parameters, software version)
#my $multiGroupsFile='multi_groups.txt' # Provided my user if multi-group search without mqpar.xml file

###>Inflating file
if ($archiveFile) {
	&addToFileStat("0/$numSteps Extracting files\n");
	my $fullFile="$tmpFilesDir/$archiveFile";
	my ($type)=($archiveFile =~ /\.(gz|zip)\Z/);
	if ($type eq 'gz') {
		if ($archiveFile =~ /\.tar\.gz\Z/) {
			system "tar -zxf $fullFile -C $tmpFilesDir";
		}
		else { # gz only
			system "gunzip -c $fullFile > $tmpFilesDir";
		}
	}
	elsif ($type eq 'zip') {
		#system "unzip -q -d $tmpFilesDir $newFile";
		&promsMod::unzipArchive($fullFile,$tmpFilesDir);
	}

	###>Deleting archive file
	unlink $fullFile;
}

###>Moving all MQ files to top work directory (Only 1st-level sub-directories are checked!) <- in case extracted from archive with extra hierarchy
opendir (DIR,$tmpFilesDir);
while (my $item = readdir (DIR)) {
	next if $item=~/^\.+$/; #  skip '.' & '..' directories
	if (-d "$tmpFilesDir/$item") { # directory
		foreach my $mqFile (glob("$tmpFilesDir/$item/*.txt"),glob("$tmpFilesDir/*.xml")) {
			my (undef,$path,$fileName) = splitpath($mqFile);
			move($mqFile,"$tmpFilesDir/$fileName");
		}
		rmtree "$tmpFilesDir/$item"; rmdir "$tmpFilesDir/$item" if -e "$tmpFilesDir/$item";
	}
}
close DIR;

my $mqparFile;
my $fullMqparFile=(glob("$tmpFilesDir/*.xml"))[0];
(undef,my $tmpFilesPath,$mqparFile) = splitpath($fullMqparFile) if $fullMqparFile;
my @requiredFiles=($evidenceFile,$peptideFile);
push @requiredFiles,$mqparFile if $mqparFile;
push @requiredFiles,$proteinGroupsFile if $importProtQuantif;
push @requiredFiles,($summaryFile,$parametersFile) if !$mqparFile;
foreach my $file (@requiredFiles) {
	if (!-e "$tmpFilesDir/$file") {
		die "ERROR: File $file is missing!";
	}
}

####>Connectecting to DB<####
my $dbh=&promsConfig::dbConnect('no_user');

###>Link to a myProMS database in order to be consistent
my %idParseRules; # allows to used sub identifier instead of whole one (eg search with UNIPROT_ALL imported as UNIPROT_ID/ACC)
my $sthPR=$dbh->prepare("SELECT PARSE_RULES,IS_CRAP FROM DATABANK D,DATABANK_TYPE DT WHERE D.ID_DBTYPE=DT.ID_DBTYPE AND ID_DATABANK=?");
foreach my $dbID (@databankIDs) {
	$sthPR->execute($dbID);
	my ($parseRules,$isCrap)=$sthPR->fetchrow_array;
	my @rules=split(',:,',$parseRules);
	($idParseRules{$dbID})=($rules[0]=~/ID=(.+)/);
}
$sthPR->finish;

#push @databankIDs,$contaminantDB if $contaminantDB;


#my @summaryColumns=('Raw file','Experiment','Enzyme','Enzyme mode','Variable modifications','Max. missed cleavages','Labels0','Labels1','Labels2');
my $dataFile=$mqparFile || $parametersFile;
my $sthInsS=$dbh->prepare("INSERT INTO SAMPLE (ID_SAMPLE,ID_EXPERIMENT,NAME,DISPLAY_POS,START_DATE,UPDATE_DATE,UPDATE_USER) VALUES (?,$experimentID,?,?,NOW(),NOW(),?)") || die $dbh->errstr;
my $sthInsA=$dbh->prepare("INSERT INTO ANALYSIS (ID_SAMPLE,NAME,START_DATE,DATA_FILE,VALID_STATUS,VALID_USER,LOWER_SCORES,FIRST_VALID_DATE,VALID_DATE,LABELING,WIFF_FILE,DECOY,DISPLAY_POS,FILE_FORMAT,MS_TYPE,MAX_RANK,MIN_SCORE,UPDATE_DATE,UPDATE_USER) VALUES (?,?,NOW(),'$dataFile',2,?,0,NOW(),NOW(),?,?,?,?,'MAXQUANT.DIR','MIS',1,0,NOW(),?)");
my $sthInsAM=$dbh->prepare("INSERT INTO ANALYSIS_MODIFICATION (ID_ANALYSIS,ID_MODIFICATION,MODIF_TYPE,SPECIFICITY) VALUES (?,?,?,?)");
my $sthLabMod=$dbh->prepare("UPDATE MODIFICATION SET IS_LABEL=1 WHERE ID_MODIFICATION=?");
my $sthInsAD=$dbh->prepare("INSERT INTO ANALYSIS_DATABANK (ID_ANALYSIS,ID_DATABANK,DB_RANK) VALUES (?,?,?)");
my $sthInsO=$dbh->prepare("INSERT INTO OBSERVATION (ID_ANALYSIS,TARGET_POS) VALUES (?,?)");
my $sthInsOM=$dbh->prepare("INSERT INTO OBS_MODIFICATION (ID_OBSERVATION,ID_MODIFICATION) VALUES(?,?)");
my ($maxSampleID) = $dbh->selectrow_array("SELECT MAX(ID_SAMPLE) FROM SAMPLE");
my ($displayPosSamp) = $dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM SAMPLE WHERE ID_EXPERIMENT=$experimentID");


#my ($quantifMethodID,$areaParamID);
my $sthInsQ=$dbh->prepare("INSERT INTO QUANTIFICATION (ID_QUANTIFICATION_METHOD,ID_DESIGN,NAME,FOCUS,STATUS,QUANTIF_ANNOT,UPDATE_USER,UPDATE_DATE) VALUES (?,?,?,?,0,?,'$userID',NOW())");
my $sthInsAQ=$dbh->prepare("INSERT INTO ANA_QUANTIFICATION (ID_QUANTIFICATION,ID_ANALYSIS) VALUES (?,?)");
my (%mqRepIon) = &getModIDforReporterIon;

my (%anaInfo,%anaNames,%vmodsInUnimod,%modifName2ID,%designInfo);
my (%anaVmods,%labels);
my (%sampIDs,%ana2Observation,%samplePos);


####> Group/global parameters & design used for MaxQuant search <####
my (%fixedMods,%varMods,%designLabels);
#my (%labelInfos,%ana2pepQuantification);
my %ana2pepQuantification;
##my ($labeling,$labelingPlex);
my ($fdrPc,$pepUsed,$peptideFilterStrg,$versionStrg,$xmlParams,$useReporterIndex);
#my $isobarCorrectedStrg=''; # for isobaric labeling only
my $mqVersion=1; #default
my (@designAnalyses,@designRawFiles,@designExperiments,@designFractions,@designGroups);
my (%ana2ParamGroup,%group2AnaName);

&addToFileStat("1/$numSteps Importing experimental design\n");
if ($mqparFile) {
	
	###>mqpar.xml<###
	print LOG "1/$numSteps Importing experimental design from mqpar file\n";
	my $xml = new XML::Simple();
	$xmlParams = $xml->XMLin("$tmpFilesDir/$mqparFile",ForceArray=>['parameterGroup','string','short','int'],SuppressEmpty=>undef);
	
	##>Maxquant XIC software & version
	$versionStrg=($xmlParams->{maxQuantVersion})? ';'.$xmlParams->{maxQuantVersion} : '';
	
	##>Experimental design
	@designExperiments=@{$xmlParams->{experiments}{string}};
	@designFractions=@{$xmlParams->{fractions}{short}};
	#@designGroups=@{$xmlParams->{paramGroupIndices}{int}};
	my @paramGroupIndex=@{$xmlParams->{paramGroupIndices}{int}};
	foreach my $grIdx (@paramGroupIndex) {
		$designGroups[$grIdx]=$grIdx; # assumes gr name is gr index!!!!!!!
	}	
	#@designRawFiles=@{$xmlParams->{filePaths}{string}};
	#@designLabels=@{$xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{labelMods}{string}};
	foreach my $anaIdx (0..$#{$xmlParams->{filePaths}{string}}) {
		my ($rawFile)=$xmlParams->{filePaths}{string}[$anaIdx]=~/([^\\\/]+)$/;
		push @designRawFiles,$rawFile;
		(my $anaName=$rawFile)=~s/\.[^\.]+$//;
		push @designAnalyses,$anaName;
		push @{$group2AnaName{$paramGroupIndex[$anaIdx]}},$anaName;
		$ana2ParamGroup{$anaName}=$paramGroupIndex[$anaIdx];
	}
		
	$fdrPc=$xmlParams->{peptideFdr} * 100;
	
	##>Global fixed modifs (older MQ versions)<##
	if ($xmlParams->{fixedModifications}) {
		foreach my $modifStrg (@{$xmlParams->{fixedModifications}{string}}) {
			next unless $modifStrg;
			my ($fixModName,$specificity)=&promsMod::convertVarModString($modifStrg);
			my $modID=&promsMod::getModificationIDfromString($dbh,$fixModName,$specificity,\%vmodsInUnimod); # %vmodsInUnimod: prevents re-parsing of unimods_table file
			$modifName2ID{$fixModName}=$modID;
			foreach my $anaName (@designAnalyses) {
				push @{$fixedMods{$anaName}},[$modID,$specificity];
			}
		}
	}
	
	##>Group-specific parameters<##
	my (%allVarModStrgs,%grFixedMods,%grVarMods);
	foreach my $paramGrIdx (0..$#designGroups) {
		my $paramGr=$designGroups[$paramGrIdx];
		#>--Fixed modifs--<# Newer MQ versions (fixed mods are group-specific)
		if ($xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{fixedModifications}) { # 'string' attribute missing if no modif at all!
			foreach my $modifStrg (@{$xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{fixedModifications}{string}}) {
				next unless $modifStrg;
				my ($fixModName,$specificity)=&promsMod::convertVarModString($modifStrg);
				my $modID=&promsMod::getModificationIDfromString($dbh,$fixModName,$specificity,\%vmodsInUnimod); # %vmodsInUnimod: prevents re-parsing of unimods_table file
				$modifName2ID{$fixModName}=$modID;
				foreach my $anaName (@{$group2AnaName{$paramGr}}) {
					push @{$fixedMods{$anaName}},[$modID,$specificity];
				}
			}
		}
		#>--Variable modifs--<#
		if ($xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{variableModifications}) {
			foreach my $modifStrg (@{$xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{variableModifications}{string}}) {
				next unless $modifStrg;
				$allVarModStrgs{$modifStrg}=1;
				my ($varModName,$specificity)=&promsMod::convertVarModString($modifStrg);
				my $modID=($modifName2ID{$varModName})? $modifName2ID{$varModName} : &promsMod::getModificationIDfromString($dbh,$varModName,$specificity,\%vmodsInUnimod);
				$modifName2ID{$varModName}=$modID;
				foreach my $anaName (@{$group2AnaName{$paramGr}}) {
					push @{$varMods{$anaName}},[$modID,$specificity,$modifStrg];
				}
			}
		}
		#>--Labeling--<#
		if ($xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{labelMods}) {
			##foreach my $anaName (@{$group2AnaName{$paramGr}}) {
			##	@{$designLabels{$anaName}}=@{$xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{labelMods}{string}};
			##}
			#@{$designLabels{$paramGr}}=@{$xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{labelMods}{string}};
			my $okLabels=0;
			my @tempLabels;
			foreach my $label (@{$xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{labelMods}{string}}) { # matches multiplicity
				if ($label) {$okLabels++;}
				else {$label='';}
				push @tempLabels,$label;
			}
			if ($okLabels) {
				@{$designLabels{$paramGr}}=@tempLabels;
			}
		}
	}
	
	#>--Peptides used for quantif--<#
	$pepUsed=($xmlParams->{quantMode}==1)? 'razor' : ($xmlParams->{quantMode}==2)? 'unique' : 'all'; # 0: all, 1: unique+razor, 2: unique
	$peptideFilterStrg="$pepUsed;1;"; # peptide used;missedCut;
	if ($xmlParams->{restrictMods}{string}[0]) { # PTM used
		$peptideFilterStrg.='-' if lc($xmlParams->{useCounterparts}) eq 'false'; # exclude unmodified matching peptides (true/True according to MQ version!!!)
		#>Check if all varMods are listed
		my @matchedVmods;
		foreach my $modifStrg (@{$xmlParams->{restrictMods}{string}}) {
			##foreach my $refVarMod (@varMods) {
			##	if ($modifStrg eq $refVarMod->[2]) {
			##		push @matchedVmods,$modifStrg;
			##		last;
			##	}
			##}
			push @matchedVmods,$modifStrg if $allVarModStrgs{$modifStrg};
		}
		if (scalar @matchedVmods == scalar keys %allVarModStrgs) { # all allowed
			$peptideFilterStrg.=($peptideFilterStrg =~ /-$/)? '3' : '1'; # "-3" -> 'All modifications allowed (unmodified matching peptides not allowed)'
		}
		else {
			$peptideFilterStrg.='2:';
			my $firstPTM=1;
			foreach my $modifStrg (@matchedVmods) {
				if ($firstPTM) {$firstPTM=0;}
				else {$peptideFilterStrg.=',';}
				my ($varModName,$specificity)=&promsMod::convertVarModString($modifStrg);
				$peptideFilterStrg.='#'.$modifName2ID{$varModName};
			}
		}
	}
	else {$peptideFilterStrg.='0';} # no PTM allowed
	$peptideFilterStrg.=';all;all'; # charge;source used
}
else { # summary.txt + parameters.txt
	my %allVarModStrgs;
	print LOG "1/$numSteps Importing experimental design from summary.txt and parameters.txt files\n";
	
	open (SUMMARY,"$tmpFilesDir/$summaryFile") || die "Unable to open $tmpFilesDir/$summaryFile";
	my (%summaryColNum,@grKeyCols,%groupKeys);
	#my $paramGrIdx=0; # Force to 0...
	my $maxParamGr=-1; # same as idx here but can be different in mqpar.xml!
	while (my $line=<SUMMARY>) {
		$line=~s/\s+$//; # chomp is not enough <- Windows
		my @parameters=split(/\t/,$line);
		if ($parameters[0] eq 'Raw file') {
			my $ncol=0;
			foreach my $colName (@parameters) {
				$summaryColNum{$colName}=$ncol;
				$ncol++;
			}
			#>List of columns used for group key
			@grKeyCols=('Enzyme','Enzyme mode','Enzyme first search','Enzyme mode first search','Use enzyme first search','Variable modifications','Fixed modifications','Multi modifications','Variable modifications first search','Use variable modifications first search','Requantify','Multiplicity','Max. missed cleavages','Max. labeled AAs');
			my $idx=0;
			while ($summaryColNum{"Labels$idx"}) {
				push @grKeyCols,"Labels$idx";
				$idx++;
			}
		}
		elsif ($parameters[0] eq 'Total' || defined($designInfo{'Experiment'}{$parameters[0]})){
			# Experiment information is not kept so skip to the end !
			last;
		}
		else {
			if ($parameters[$summaryColNum{'LC-MS run type'}] =~ /Reporter ion MS/) {
				#code
				print LOG "WARNING: for Reporter ion MS2 or MS3, you need to provide the mqpar.xml file.\n";
				exit;
			}
			###> Add sample according to experiment rawname
			my $anaName=$parameters[$summaryColNum{'Raw file'}];
			push @designAnalyses,$anaName;
			push @designRawFiles,$anaName.'raw'; # assume .raw file format
			push @designExperiments,$parameters[$summaryColNum{'Experiment'}] if $summaryColNum{'Experiment'};
			push @designFractions,$parameters[$summaryColNum{'Fraction'}] if $summaryColNum{'Fraction'};
			#if ($#designLabels < 0) {
			my $okLabel=0;
			my @tempLabels;
			my $idx=0;
			while ($summaryColNum{"Labels$idx"}) {
				my $label='';
				if ($parameters[$summaryColNum{"Labels$idx"}]) {
					$label=$parameters[$summaryColNum{"Labels$idx"}];
					$okLabel++;
				}
				push @tempLabels,$label;
				$idx++;
				last if $idx == $parameters[$summaryColNum{'Multiplicity'}];
			}
			if ($okLabel) { # must have a value in at least 1 of label channels
				#if ($summaryColNum{'Labels0'}) {
				#	$parameters[$summaryColNum{'Labels0'}]='' unless $parameters[$summaryColNum{'Labels0'}];
				#}
				#push @{$designLabels{$anaName}},$parameters[$summaryColNum{'Labels0'}] if $summaryColNum{'Labels0'}; # @designLabels
				#push @{$designLabels{$anaName}},$parameters[$summaryColNum{'Labels1'}] if $summaryColNum{'Labels1'};
				#push @{$designLabels{$anaName}},$parameters[$summaryColNum{'Labels2'}] if $summaryColNum{'Labels2'};
				
				@{$designLabels{$anaName}}=@tempLabels;
			}
				
			if ($summaryColNum{'Fixed modifications'} && $parameters[$summaryColNum{'Fixed modifications'}]) {
				foreach my $modifStrg (split(/;/,$parameters[$summaryColNum{'Fixed modifications'}])) {
					next unless $modifStrg;
					my ($fixModName,$specificity)=&promsMod::convertVarModString($modifStrg);
					my $modID=($modifName2ID{$fixModName})? $modifName2ID{$fixModName} : &promsMod::getModificationIDfromString($dbh,$fixModName,$specificity,\%vmodsInUnimod);
					push @{$fixedMods{$anaName}},[$modID,$specificity];
					$modifName2ID{$fixModName}=$modID;
				}
			}
			
			if ($summaryColNum{'Variable modifications'} && $parameters[$summaryColNum{'Variable modifications'}]) {
				foreach my $modifStrg (split(/;/,$parameters[$summaryColNum{'Variable modifications'}])) {
					next unless $modifStrg;
					$allVarModStrgs{$modifStrg}=1;
					my ($varModName,$specificity)=&promsMod::convertVarModString($modifStrg);
					my $modID=($modifName2ID{$varModName})? $modifName2ID{$varModName} : &promsMod::getModificationIDfromString($dbh,$varModName,$specificity,\%vmodsInUnimod);
					push @{$varMods{$anaName}},[$modID,$specificity,$modifStrg];
					$modifName2ID{$varModName}=$modID;
				}
			}
			#}
			my ($sampleName)=($summaryColNum{'Experiment'}) ? $parameters[$summaryColNum{'Experiment'}] : "No_experiments_defined_in_maxquant";
			%{$designInfo{'Experiment'}{$sampleName}}=();
			
			#>Detecting parameter groups
			my @keyValues;
			foreach my $colName (@grKeyCols) {
				next unless $summaryColNum{$colName};
				my $value=defined($parameters[$summaryColNum{$colName}])? $parameters[$summaryColNum{$colName}] : '';
				push @keyValues,$value;
			}
			my $grKey=join('::',@keyValues);
			unless (defined $groupKeys{$grKey}) {
				$maxParamGr++;
				$groupKeys{$grKey}=$maxParamGr;
				push @designGroups,$maxParamGr;
			}
			my $paramGr=$groupKeys{$grKey};
			$ana2ParamGroup{$anaName}=$paramGr;
			push @{$group2AnaName{$paramGr}},$anaName;
			
			unless (defined $designLabels{$paramGr}) {
				my $label0=$parameters[$summaryColNum{'Labels0'}] if $summaryColNum{'Labels0'};
				my $label1=$parameters[$summaryColNum{'Labels1'}] if $summaryColNum{'Labels1'};
				my $label2=$parameters[$summaryColNum{'Labels2'}] if $summaryColNum{'Labels2'};
				if ($label0) {
					push @{$designLabels{$paramGr}},$label0;
				}
				if ($label1) {
					push @{$designLabels{$paramGr}},'' unless $label0;
					push @{$designLabels{$paramGr}},$label1;
				}
				push @{$designLabels{$paramGr}},$label2 if $label2;
				#if ($summaryColNum{'Labels0'}) {
				#	$parameters[$summaryColNum{'Labels0'}]='' unless $parameters[$summaryColNum{'Labels0'}];
				#}
				#push @{$designLabels{$paramGr}},$parameters[$summaryColNum{'Labels0'}] if $summaryColNum{'Labels0'}; # @designLabels
				#push @{$designLabels{$paramGr}},$parameters[$summaryColNum{'Labels1'}] if $summaryColNum{'Labels1'};
				#push @{$designLabels{$paramGr}},$parameters[$summaryColNum{'Labels2'}] if $summaryColNum{'Labels2'};
			}
		}
	}
	close SUMMARY;
	
	#if (scalar keys %groupKeys > 1) {
	#	###>Checking for multi_groups.txt file<###
	#	if (-e "$tmpFilesDir/$multiGroupsFile") {
	#		open (MULTIGR,"$tmpFilesDir/$multiGroupsFile") || die "Unable to open $tmpFilesDir/$multiGroupsFile";
	#		while (my $line=<MULTIGR>) {
	#			$line=~s/\s+$//; # chomp is not enough <- Windows
	#			#my ($rawFile,$grIdx)=split(/\t/,$line);
	#			my ($anaName,$grIdx)=$line=~/^\s*(.+\w)\W+(\d+)$/;
	#			$anaName=~s/\.\w{3}$//; # in case file extension was left
	#			unless (defined $anaParamGroup{$anaName}) {
	#				$dbh->disconnect;
	#				close LOG;
	#				close MULTIGR;
	#				die "ERROR: File '$anaName' is not part of search results";
	#			}
	#			$anaParamGroup{$anaName}=$grIdx;
	#		}
	#		close MULTIGR;
	#	}
	#	else {
	#		$dbh->disconnect;
	#		close LOG;
	#		die "ERROR: Import of multi-groups search without mqpar.xml file requires a user-made $multiGroupsFile file (see import form for details)";
	#	}
	#}
	
	my $addPepString='';
	open (PARAMFILE,"$tmpFilesDir/$parametersFile") || die "Unable to open $tmpFilesDir/$parametersFile";
	while (my $line=<PARAMFILE>) {
		$line=~s/\s+$//; # chomp is not enough <- Windows
		my @parameters=split(/\t/,$line);
		if ($parameters[0] eq 'Version') {
			##>Maxquant XIC software & version
			$parameters[1]=~s/\s+$//; # trailing spaces if any
			($mqVersion)=($parameters[1]=~/^(\d+\.\d+)/); # x.x.x.xx -> x.x (numerical value)
			$versionStrg=";$parameters[1]";
		}
		elsif ($parameters[0] eq 'PSM FDR') {
			$fdrPc=$parameters[1] * 100;
		}
		elsif ($parameters[0] eq 'Fixed modifications') { # old MQ version: Fixed mods are not group-specific
			next if (!$parameters[1] || $parameters[1] !~ /\w+/);
			my $fixStg=$parameters[1];
			$fixStg=~ s/\s+$//g; # Remove last whitespace: 'Carbamidomethyl (C) ' becomes 'Carbamidomethyl (C)'
			foreach my $modifStrg (split(/;/,$fixStg)){
				next unless $modifStrg;
				my ($fixModName,$specificity)=&promsMod::convertVarModString($modifStrg);
				my $modID=&promsMod::getModificationIDfromString($dbh,$fixModName,$specificity,\%vmodsInUnimod); # %vmodsInUnimod: prevents re-parsing of unimods_table file
				foreach my $anaName (@designAnalyses) {
					push @{$fixedMods{$anaName}},[$modID,$specificity];
				}
				$modifName2ID{$fixModName}=$modID;
			}
		}
		elsif ($parameters[0] eq 'Modifications included in protein quantification') {
			$peptideFilterStrg=';1;'; # ;missedCut;
			#>Check if all varMods are listed
			my @matchedVmods;
			foreach my $modifStrg (split(/;/,$parameters[1])) {
				##foreach my $refVarMod (@varMods) {
				##	if ($modifStrg eq $refVarMod->[2]) {
				##		push @matchedVmods,$modifStrg;
				##		last;
				##	}
				##}
				push @matchedVmods,$modifStrg if $allVarModStrgs{$modifStrg};
			}
			if (scalar @matchedVmods == scalar keys %allVarModStrgs) {$addPepString='1';} # all allowed
			else {
				$addPepString='2:';
				my $firstPTM=1;
				foreach my $modifStrg (@matchedVmods) {
					if ($firstPTM) {$firstPTM=0;}
					else {$addPepString.=',';}
					my ($varModName,$specificity)=&promsMod::convertVarModString($modifStrg);
					$addPepString.='#'.$modifName2ID{$varModName};
				}
			}
		}
		elsif ($parameters[0] eq 'Peptides used for protein quantification') {
			$pepUsed=($parameters[1]=~ /Razor/i)? 'razor' : ($parameters[1]=~ /Unique/i)? 'unique' : 'all';
		}
		elsif ($parameters[0] eq 'Discard unmodified counterpart peptides') {
			if ($parameters[1] =~ /True/i) {
				$peptideFilterStrg.='-';
				$addPepString='3' if $addPepString eq '1'; # overwritten to allow "-3" -> 'All modifications allowed (unmodified matching peptides not allowed)'
			}
			$peptideFilterStrg=$pepUsed.$peptideFilterStrg.$addPepString;
		}
	}
	close PARAMFILE;
	if (!$peptideFilterStrg) { $peptideFilterStrg="$pepUsed;1;0"; }
	$peptideFilterStrg.=';all;all'; # charge;source used
}

print LOG " -",scalar @designGroups," parameter group(s) detected\n";

##>Labeling<##
my %LABELCLASS_POS=('Label-free'=>1,'SILAC'=>2,'iTRAQ'=>3,'TMT'=>4,'Isobaric*'=>5);
my (%labeling,%labelingPlex,%distinctLabelings,%group2Labeling);
my (%labelModifSpecificity,%isotopeLabelInfo,%labelName,%labelList,%label2targetPos);
#my ($pepQuantifName,$labelingName,$pepQuantifAnnot,$quantifMethodID,$areaParamID,$intensityParamID); # intensity for isobaric
my (%labelingName,%peptideQuantifID,%areaParamID,%intensityParamID); # intensity for isobaric
my %isobarModifInfo; # for isobaric labeling only;
my $isobarCorrectedStrg=''; # if multiple isboric parameter groups: assume correction applies (or not) to all

foreach my $paramGrIdx (0..$#designGroups) {
	my $paramGr=$designGroups[$paramGrIdx];

	my ($pepQuantifName,$pepQuantifAnnot,$quantifMethodID,$labelClass);

	#>------SILAC-----<# !!! Actually any Isotopic labeling !!!
	#if (($xmlParams && $xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{multiplicity} > 1)||$#designLabels>1) { # $xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{labelMods}{string}[1]
	if ($designLabels{$paramGr}) { # not defined for Label-free
		$labelingPlex{$paramGr}=($xmlParams)? $xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{multiplicity} : 1+$#{$designLabels{$paramGr}}; #scalar @{$xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{labelMods}{string}};
		$labelClass='SILAC';
		$labeling{$paramGr}='SILAC';
		$labelingName{$paramGr}='SILAC '.$labelingPlex{$paramGr}.'plex';
		$pepQuantifName='SILAC made by MaxQuant';
		print LOG " -$labelingName{$paramGr} detected in group $paramGr\n";
	
		$pepQuantifAnnot='LABEL=SILAC::SOFTWARE=MQ'.$versionStrg;
		%{$labelName{$paramGr}}=('L'=>'Light','M'=>'Medium','H'=>'Heavy');
		my (@colNames,@labNames);
		if ($labelingPlex{$paramGr} == 1) { # unlikely but technically possible in MaxQuant
				@colNames=('H'); # or M?
				#@labNames=('Heavy');
				#%{$label2targetPos{$paramGr}}=('H'=>1);
		}
		elsif ($labelingPlex{$paramGr} == 2) { # possible ambiguity L/H or L/M or M/H? => read evidence.txt
			#my $headEvStrg=`head -1 $tmpFilesDir/$evidenceFile`;
			#$headEvStrg=~s/\s+$//; # chomp is not enough <- Windows
			#my @colHeaders=split(/\t/,$headEvStrg);
			#my $targetPos=0;
			#foreach my $colLabel ('L','M','H') {
			#	foreach my $header (@colHeaders) {
			#		if ($header eq "Intensity $colLabel") {
			#			push @colNames,$colLabel;
			#			push @labNames,$labelName{$paramGr}{$colLabel};
			#			$label2targetPos{$paramGr}{$colLabel}=++$targetPos;
			#			last;
			#		}
			#	}
			#}
			if ($designLabels{$paramGr}[0] eq '') { # L,H
				@colNames=('L','H');
				#@labNames=('Light','Heavy');
				#%{$label2targetPos{$paramGr}}=('L'=>1,'H'=>2);
			}
			else {
				@colNames=('M','H');
				#@labNames=('Medium','Heavy');
				#%{$label2targetPos{$paramGr}}=('M'=>1,'H'=>2);
			}
		}
		else { # triplex...
			@colNames=('L','M','H');
			#@labNames=('Light','Medium','Heavy');
			#%{$label2targetPos{$paramGr}}=('L'=>1,'M'=>2,'H'=>3);
		}
		@{$labelList{$paramGr}}=@colNames;
		my $targetPos=0;
		foreach my $colLabel (@colNames) {
			push @labNames,$labelName{$paramGr}{$colLabel};
			$label2targetPos{$paramGr}{$colLabel}=++$targetPos;
		}
	
		foreach my $labelIdx (0..$#{$designLabels{$paramGr}}) {
			#my $label=($labelIdx==0)? 'NA' :
			my @modifications=($designLabels{$paramGr}[$labelIdx] eq '')? ('None') : split(/;/,$designLabels{$paramGr}[$labelIdx]);
			my $labelCol=shift @colNames;
			my $labName=shift @labNames;
			my $targetPos=$labelIdx+1;
			$pepQuantifAnnot.="::$targetPos;$labName;";
			my @residueMods;
			foreach my $repIon (@modifications) {
				if ($repIon eq 'None') {
					push @residueMods,'None#No label##';
				}
				else {
					my $specificity=($repIon =~ /Lys/)? 'K' : ($repIon =~ /Arg/)? 'R' : ($repIon =~ /Nter/)? '=' : ($repIon =~ /Ile7/)? 'I': ($repIon =~ /Leu7/)? 'L' : ($repIon =~ /ICAT/)? 'C' : ($repIon =~ /18O/)? '*' : '';
					my ($unimodID)=$mqRepIon{$repIon}{'UNIMOD_ACC'};
					my ($modID,$psiName,$interimName)=$dbh->selectrow_array("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME FROM MODIFICATION WHERE UNIMOD_ACC=$mqRepIon{$repIon}{UNIMOD_ACC}");
					unless ($modID) { # not yet in DB
						$modID=&promsMod::getModificationIDfromString($dbh,$mqRepIon{$repIon}{'INTERIM_NAME'},$specificity,\%vmodsInUnimod);
						$sthLabMod->execute($modID); # set as label modif (just to be safe)
						($psiName,$interimName)=$dbh->selectrow_array("SELECT PSI_MS_NAME,INTERIM_NAME FROM MODIFICATION WHERE ID_MODIFICATION=$modID");
					}
					$psiName=$interimName unless $psiName;
					$psiName='' unless $psiName;
					push @residueMods,"$repIon#$psiName#$specificity#$psiName#$modID";
					$isotopeLabelInfo{$paramGr}{$labelCol}{$modID}{$specificity}=1;
					$labelModifSpecificity{$paramGr}{$modID}{$specificity}=1;
				}
			}
			$pepQuantifAnnot.=join('@',@residueMods);
		}
		($quantifMethodID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='SILAC'"); # SILAC
		($areaParamID{$paramGr})=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE='ISO_AREA' AND ID_QUANTIFICATION_METHOD=$quantifMethodID");
	}
	#>------iTRAQ or TMT-----<#
	elsif ($xmlParams && $xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{isobaricLabels}) { # iTRAQ, TMT
		my %isobarTag;
		my %aaThree2OneLetter=&promsMod::getAAThree2OneLetter;
		($mqVersion)=($xmlParams->{maxQuantVersion}=~/^(\d+\.\d+)/); # x.x.x.xx -> x.x (numerical value)
		my $tagPos=0;
		if ($mqVersion < 1.6) {
			foreach my $modifTgt (@{$xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{isobaricLabels}{string}}) {
				#my ($modif,$res,$tag)=($modifTgt=~/(.+)-([^-\d]+)(\d+)$/);
				###> Examples:
				###> TMT: TMT10plex-Nter126C,  TMT10plex-Lys127C, etc.
				###> iTRAQ: iTRAQ8plex-Nter113, iTRAQ8plex-Lys114, etc.
				my ($modif,$res,$tag)=($modifTgt=~/(.+)-([^-\d]+)(\d+\w?)$/);
				$isobarModifInfo{$paramGr}{NAME}=$modif;
				$res=($res eq 'Nter')? '=' : ($res eq 'Cter')? '*' : $aaThree2OneLetter{$res};
				$isobarModifInfo{$paramGr}{RES}{$res}=1;
				$isobarTag{$tag}=++$tagPos;
			}
		}
		else { # Seen for 1.6.3.4 maxQuantVersion
			foreach my $isobaricLabelInfo (@{$xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{isobaricLabels}{IsobaricLabelInfo}}) {
				#<IsobaricLabelInfo>
				#   <internalLabel>TMT10plex-Lys126C</internalLabel>
				#   <terminalLabel>TMT10plex-Nter126C</terminalLabel>
				#   <correctionFactorM2>0</correctionFactorM2>
				#   <correctionFactorM1>0</correctionFactorM2>
				#   <correctionFactorP1>8.2</correctionFactorP1>
				#   <correctionFactorP2>0.4</correctionFactorP2>
				#   <tmtLike>True</tmtLike>
				#</IsobaricLabelInfo>
				
				#>Internal / Terminal
				foreach my $location ('internalLabel','terminalLabel') {
					my $modifTgt=$isobaricLabelInfo->{$location};
					next unless $modifTgt;
					my ($modif,$res,$tag)=($modifTgt=~/(.+)-([^-\d]+)(\d+\w?)$/);
					$isobarModifInfo{$paramGr}{NAME}=$modif;
					$res=($res eq 'Nter')? '=' : ($res eq 'Cter')? '*' : $aaThree2OneLetter{$res};
					$isobarModifInfo{$paramGr}{RES}{$res}=1;
					$isobarTag{$tag}=++$tagPos unless $isobarTag{$tag}; # only once
				}
				
				#>Collect correction factors
				$isobarModifInfo{$paramGr}{CORR}{$tagPos}=join(',',($isobaricLabelInfo->{correctionFactorM2},$isobaricLabelInfo->{correctionFactorM1},$isobaricLabelInfo->{correctionFactorP1},$isobaricLabelInfo->{correctionFactorP2})) if defined $isobaricLabelInfo->{correctionFactorP1};
			}
		}
		if ($isobarModifInfo{$paramGr}{NAME}=~/iTRAQ/i) {
			$labeling{$paramGr}='ITRAQ';
			$pepQuantifName='iTRAQ made by MaxQuant';
			$labelClass='iTRAQ';
		}
		elsif ($isobarModifInfo{$paramGr}{NAME}=~/TMT/i) {
			$labeling{$paramGr}='TMT';
			$pepQuantifName='TMT made by MaxQuant';
			$labelClass='TMT';
		}
		else { # Warning: Should not happen!
			$labeling{$paramGr}='ISOBAR';
			$pepQuantifName='Isobaric made by MaxQuant';
			$labelClass='Isobaric*';
		}
		($labelingPlex{$paramGr})=($isobarModifInfo{$paramGr}{NAME}=~/(\d+)plex/i); # or scalar keys %isobarTag
		$labelingName{$paramGr}=$isobarModifInfo{$paramGr}{NAME};
		print LOG " -$labelingName{$paramGr} labeling detected in group $paramGr\n";
	
		$pepQuantifAnnot='LABEL='.$labeling{$paramGr}.'::SOFTWARE=MQ'.$versionStrg;
		@{$isobarModifInfo{$paramGr}{TAGS}}=();
		foreach my $tag (sort{$isobarTag{$a}<=>$isobarTag{$b}} keys %isobarTag) {
			push @{$isobarModifInfo{$paramGr}{TAGS}},$tag; # index!!!
			$pepQuantifAnnot.="::$isobarTag{$tag};$tag;"; # tagPos;tag
		}
		if ($isobarModifInfo{$paramGr}{CORR}) {
			$pepQuantifAnnot.='::CORRECTION=#0;'; # '#0' because not linked to a product sheet ID in myProMS & NOT the SAME format (Minus2,M1,Plus1,P2)!!!
			foreach my $tagPos (1..$labelingPlex{$paramGr}) {
				$pepQuantifAnnot.='&' if $tagPos > 1;
				$pepQuantifAnnot.=$isobarModifInfo{$paramGr}{CORR}{$tagPos} || '';
			}
		}
		my $specificity=join(',',sort keys %{$isobarModifInfo{$paramGr}{RES}});
		$isobarModifInfo{$paramGr}{ID}=&promsMod::getModificationIDfromString($dbh,$isobarModifInfo{$paramGr}{NAME},$specificity,\%vmodsInUnimod);
		$sthLabMod->execute($isobarModifInfo{$paramGr}{ID}); # set as label modif (just to be safe)
		$labelModifSpecificity{$paramGr}{$isobarModifInfo{$paramGr}{ID}}=$isobarModifInfo{$paramGr}{RES};
		($quantifMethodID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='$labeling{$paramGr}'");
		($intensityParamID{$paramGr})=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE='REP_INTENSITY' AND ID_QUANTIFICATION_METHOD=$quantifMethodID");
		($isobarModifInfo{$paramGr}{MONO_MASS})=$dbh->selectrow_array("SELECT MONO_MASS FROM MODIFICATION WHERE ID_MODIFICATION=$isobarModifInfo{$paramGr}{ID}");
		unless ($isobarModifInfo{$paramGr}{MONO_MASS}) {
			print LOG "[WARNING: Could not retrieve mass data for '$isobarModifInfo{$paramGr}{NAME}'! Complete this modification entry in myProMS for optimal data import.]";
		}
		foreach my $repIdx (0..($labelingPlex{$paramGr}-1)) {
			$label2targetPos{$paramGr}{$repIdx}=$repIdx+1;
			$labelName{$paramGr}{$repIdx}=$isobarModifInfo{$paramGr}{TAGS}[$repIdx];
			push @{$labelList{$paramGr}},$repIdx;
		}
		
		##>Checking index vs position (MQ version dependent: 'position' if v  >= 1.6?)
		$useReporterIndex=`head -1 $tmpFilesDir/$evidenceFile | grep -c 'Reporter intensity 0'`;
		$useReporterIndex=~s/\D//s; # 0 or 1
		$useReporterIndex*=1;
	
		##>Checking correction vs no correction
		unless ($isobarCorrectedStrg) {
			my $isCorrectedStrg=`head -1 $tmpFilesDir/$evidenceFile | grep -c 'Reporter intensity corrected'`;
			$isobarCorrectedStrg=($isCorrectedStrg=~/^0/)? 'not corrected' : 'corrected';
		}
	}
	#>------label-free-----<#
	else {
		$labeling{$paramGr}='FREE';
		$labelClass='Label-free';
		$labelingPlex{$paramGr}=0;
		$pepQuantifName='Label-free made by MaxQuant';
		$labelingName{$paramGr}=undef;
		print LOG " -No labeling detected in group $paramGr\n";
	
		$pepQuantifAnnot='LABEL=FREE::SOFTWARE=MQ'.$versionStrg;
		($quantifMethodID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='XIC'"); # label-free
		($areaParamID{$paramGr})=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE='XIC_AREA' AND ID_QUANTIFICATION_METHOD=$quantifMethodID");
	}
	
	##>Insert peptide QUANTIFICATION
	$sthInsQ->execute($quantifMethodID,undef,$pepQuantifName,'peptide',$pepQuantifAnnot);
	$peptideQuantifID{$paramGr}=$dbh->last_insert_id(undef,undef,'QUANTIFICATION','ID_QUANTIFICATION');
	
	push @{$distinctLabelings{$labelClass}},$paramGr;
	$group2Labeling{$paramGr}=$labelClass;
}

$sthLabMod->finish;

##>MaxQuant analysis design<##
my %displayPosAna;
#my $firstSampID;
my $firstAnaID;
my $expCondPos=0;
#my $cgiUnixDir=getcwd;
#chdir "$promsPath{peptide}/proj_$projectID"; # to make relative symlinks to parameter files
foreach my $anaIdx (0..$#designAnalyses) {
	##if ($xmlParams && $xmlParams->{paramGroupIndices}{int}[$anaIdx] != $paramGrIdx) {
	##	$dbh->disconnect;
	##	close LOG;
	##	die "ERROR: Import of more than 1 group is not supported (see paramGroupIndices in mqpar.xml)";
	##	#rmtree $tmpFilesDir;
	##	#exit;
	##}
	#my $paramGrIdx=($xmlParams)? $xmlParams->{paramGroupIndices}{int}[$anaIdx] : 0;
	
	my $anaName=$designAnalyses[$anaIdx];
	my $paramGr=$ana2ParamGroup{$anaName};
	my $labelClass=$group2Labeling{$paramGr};
	my ($sampleName,$realSampleName)=($#designExperiments>0 && $designExperiments[$anaIdx])? ($designExperiments[$anaIdx],' '.$designExperiments[$anaIdx]) : ("No experiment",'');
	$designInfo{'SamplesInGroup'}{$paramGr}{$sampleName}=1;
	
	my $fraction;
	$fraction=($#designFractions>0 && $designFractions[$anaIdx])? $designFractions[$anaIdx] : 0;
	$fraction=0 if $fraction==32767; # 32767 if no fraction

	unless ($sampIDs{$sampleName}) {
		$maxSampleID++;
		$displayPosSamp++;
		$sthInsS->execute($maxSampleID,$sampleName,$displayPosSamp,$userID) || die $dbh->errstr;
		$sampIDs{$sampleName}=$maxSampleID;
		#$firstSampID=$maxSampleID unless $firstSampID;
		$displayPosAna{$sampleName}=0; # sample-specific needed in case mixture of fractionnated & non-fractionnated samples (sample files may be mixted)
	}
	$displayPosAna{$sampleName}++;
	my $usedDisplayPos=$fraction || $displayPosAna{$sampleName}; # WARNING: There should not be mixture of fractions and no fraction for the same sample!!!

	#my ($rawFile)=($designRawFiles[$anaIdx]=~/([^\\\/]+)$/);
	my $rawFile=$designRawFiles[$anaIdx];
	#(my $anaName=$rawFile)=~s/\.[^\.]+$//;
	#my $anaName=$designAnalyses[$anaIdx];

	$sthInsA->execute($sampIDs{$sampleName},$anaName,$userID,$labelingName{$paramGr},$rawFile,"INT:SEARCH,FDR=$fdrPc:precomputed",$usedDisplayPos,$userID) || die $dbh->errstr;
	my $analysisID=$dbh->last_insert_id(undef,undef,'ANALYSIS','ID_ANALYSIS');

	$anaNames{$analysisID}=$anaName;
	$anaInfo{$anaName}{'ID_ANALYSIS'}=$analysisID;

	##>data directories
	mkdir "$promsPath{peptide}/proj_$projectID" unless -e "$promsPath{peptide}/proj_$projectID";
	mkdir "$promsPath{peptide}/proj_$projectID/ana_$analysisID";
	if ($mqparFile) {
		if ($anaIdx==0) {
			$firstAnaID=$analysisID;
			copy("$tmpFilesDir/$mqparFile","$promsPath{peptide}/proj_$projectID/ana_$analysisID/$mqparFile");
		}
		else {
			symlink("../ana_$firstAnaID/$mqparFile","$promsPath{peptide}/proj_$projectID/ana_$analysisID/$mqparFile"); # relative symbolic link with relative path to source
		}
	}
	else {
		if ($anaIdx==0) {
			$firstAnaID=$analysisID;
			copy("$tmpFilesDir/$parametersFile","$promsPath{peptide}/proj_$projectID/ana_$analysisID/$parametersFile");
			copy("$tmpFilesDir/$summaryFile","$promsPath{peptide}/proj_$projectID/ana_$analysisID/$summaryFile");
		}
		else {
			symlink("../ana_$firstAnaID/$parametersFile","$promsPath{peptide}/proj_$projectID/ana_$analysisID/$parametersFile"); # symbolic link with relative path to source
			symlink("../ana_$firstAnaID/$summaryFile","$promsPath{peptide}/proj_$projectID/ana_$analysisID/$summaryFile"); # symbolic link with relative path to source
		}
	}

	#>ANALYSIS_DATABANK
	my $dbRank=1;
	foreach my $dbID (@databankIDs) {
		$sthInsAD->execute($analysisID,$dbID,$dbRank);
		$dbRank++;
	}

	#>ANALYSIS_MODIFICATION
	#-Fixed
	if ($fixedMods{$anaName}) {
		foreach my $refMod (@{$fixedMods{$anaName}}) { # [modID,specif]
			$sthInsAM->execute($analysisID,$refMod->[0],'F',$refMod->[1]);
		}
	}
	#-Variable
	if ($varMods{$anaName}) {
		foreach my $refMod (@{$varMods{$anaName}}) { # [vModID,Specif,vModStrg]
			$sthInsAM->execute($analysisID,$refMod->[0],'V',$refMod->[1]);
			$anaVmods{$analysisID}{$refMod->[2]}=$refMod->[0];
		}
	}
	#-Labeling
	if ($labelModifSpecificity{$paramGr}) {
		foreach my $modID (keys %{$labelModifSpecificity{$paramGr}}) {
			$sthInsAM->execute($analysisID,$modID,'V',join(',',sort keys %{$labelModifSpecificity{$paramGr}{$modID}}));
		}
	}

	#>ANALYSIS_QUANTIFICATION
	$sthInsAQ->execute($peptideQuantifID{$paramGr},$analysisID);
	$ana2pepQuantification{$analysisID}=$peptideQuantifID{$paramGr};
#next unless $importProtQuantif;

	if ($labeling{$paramGr} eq 'SILAC') {
		#my @colNames=('L','H');
		##$labels{$analysisID}{'Labels0'}='NA'; # add label-free channel
		##$labels{$analysisID}{'Labels1'}=$parameters[$summaryColNum{'Labels1'}];
		#if ($labelingPlex{$paramGr}==3) {
		#	@colNames=('L','M','H');
		#	#$labels{$analysisID}{'Labels2'}=$parameters[$summaryColNum{'Labels2'}];
		#}
		foreach my $colName (@{$labelList{$paramGr}}) { # @colNames
			$sthInsO->execute($analysisID,-1); # $targetPosObs
			my $obsID=$dbh->last_insert_id(undef,undef,'OBSERVATION','ID_OBSERVATION');
			$ana2Observation{$analysisID}{$colName}=$obsID;
			if ($isotopeLabelInfo{$paramGr} && $isotopeLabelInfo{$paramGr}{$colName}) {
				foreach my $modID (keys %{$isotopeLabelInfo{$paramGr}{$colName}}) {
					$sthInsOM->execute($obsID,$modID);
				}
			}
			my $expCondName=($labelingPlex{$paramGr}==1)? $sampleName : $colName.$realSampleName;
			$designInfo{'Expcondition'}{$expCondName}{'ANALYSIS'}{$analysisID}=1;
			$designInfo{'ExpconditionOriginalName'}{$sampleName}{$expCondName}=1;
			$designInfo{'Fraction'}{$expCondName}{$fraction}{$analysisID}=1 if $fraction;
			$designInfo{'Expcondition2Label'}{$expCondName}=$colName; # Keep a memory of which label was used in order to link Observation and modification in the end
			$designInfo{'Expcondition'}{$expCondName}{'POS'}=++$expCondPos unless $designInfo{'Expcondition'}{$expCondName}{'POS'}; # only once/sample
		}
	}
	#>------Isobaric design------<#
	elsif ($isobarModifInfo{$paramGr}) { # $labeling=~/ITRAQ|TMT/
		foreach my $repIdx (0..($labelingPlex{$paramGr}-1)) {
			$sthInsO->execute($analysisID,$repIdx+1); # $targetPosObs
			my $obsID=$dbh->last_insert_id(undef,undef,'OBSERVATION','ID_OBSERVATION');
			$ana2Observation{$analysisID}{$repIdx}=$obsID;
			$sthInsOM->execute($obsID,$isobarModifInfo{$paramGr}{ID});
			#my $expCondName="Reporter intensity $repIdx$realSampleName";
			my $expCondName=($labelingPlex{$paramGr}==1)? $sampleName : $repIdx.$realSampleName; # not tested with 1-plex!!!
			$designInfo{'Expcondition'}{$expCondName}{'ANALYSIS'}{$analysisID}=1;
			$designInfo{'ExpconditionOriginalName'}{$sampleName}{$expCondName}=1;
			$designInfo{'Fraction'}{$expCondName}{$fraction}{$analysisID}=1 if $fraction;
			$designInfo{'Expcondition2Label'}{$expCondName}=$repIdx; #"Reporter intensity $repIdx"; # Keep a memory of which label was used in order to link Observation and modification in the end
			$designInfo{'Expcondition'}{$expCondName}{'POS'}=++$expCondPos unless $designInfo{'Expcondition'}{$expCondName}{'POS'}; # only once/sample
		}
	}
	elsif ($labeling{$paramGr} eq 'FREE') { # Label-free
		$sthInsO->execute($analysisID,0); # $targetPosObs
		$ana2Observation{$analysisID}{'NONE'}=$dbh->last_insert_id(undef, undef, 'OBSERVATION', 'ID_OBSERVATION');
		$designInfo{'Expcondition'}{$sampleName}{'ANALYSIS'}{$analysisID}=1;
		$designInfo{'ExpconditionOriginalName'}{$sampleName}{$sampleName}=1;
		$designInfo{'Fraction'}{$sampleName}{$fraction}{$analysisID}=1 if $fraction;
		$designInfo{'Expcondition'}{$sampleName}{'POS'}=++$expCondPos unless $designInfo{'Expcondition'}{$sampleName}{'POS'}; # only once/sample
#$designInfo{'Expcondition2Label'}{$ampleName}='NONE';
	}
	$designInfo{'Experiment'}{$sampleName}{$analysisID}=1;
}
#chdir $cgiUnixDir;
$sthInsAM->finish;
$sthInsOM->finish;
$sthInsS->finish;
$sthInsA->finish;
$sthInsAD->finish;
$sthInsO->finish;

my $numSamples=scalar keys %sampIDs;
my $sampStrg=($numSamples==1)? 'Sample' : 'Samples';
my $numAnalyses=scalar keys %anaNames;
my $anaStrg=($numAnalyses==1)? 'Analysis' : 'Analyses';
print LOG " -$numSamples $sampStrg and $numAnalyses $anaStrg created.\n";


###############################################
####>Check mqpar & txt files compatibility<####
###############################################
my $peptideHeadStrg=`head -1 $tmpFilesDir/$peptideFile`;
$peptideHeadStrg=~s/\s+$//; # chomp is not enough <- Windows
my ($numExpPeptide,$numExpMatched)=(0,0);
my $numExpParam=($numSamples==1 && (keys %sampIDs)[0] eq 'No experiment')? 0 : $numSamples;
foreach my $colHeader (split(/\t/,$peptideHeadStrg)) {
	$colHeader=~s/^\s+//; $colHeader=~s/\s+$//; # remove starting & trailing spaces if any
	if ($colHeader=~/^Experiment (.+)/) {
		my $expName=$1;
		$numExpMatched++ if $sampIDs{$expName};
		$numExpPeptide++;
	}
}
if ($numExpMatched < $numExpParam || $numExpPeptide > $numExpParam) {
	$dbh->rollback;
	$dbh->disconnect;
	#rmtree $tmpFilesDir;
	die "ERROR: Parameter file(s) do(es) not match data files!";
	#exit;
}

########################
####>myproMS Design<####
########################
#print "<BR><FONT class='title3'>Generating protein quantification design(s)...";
my %labelClass2Design;
my $importDate=strftime("%Y/%m/%d %H:%M:%S",localtime);
my ($maxDesignID)=$dbh->selectrow_array("SELECT MAX(ID_DESIGN) FROM DESIGN");
#$designID++;
#$dbh->do("INSERT INTO DESIGN (ID_DESIGN,ID_EXPERIMENT,NAME,DES,UPDATE_DATE,UPDATE_USER) VALUES ($designID,$experimentID,'MaxQuant [$importDate]','Automatically generated during MaxQuant data import',NOW(),'$userID')");
my $sthDesign=$dbh->prepare("INSERT INTO DESIGN (ID_DESIGN,ID_EXPERIMENT,NAME,DES,UPDATE_DATE,UPDATE_USER) VALUES (?,$experimentID,?,'Automatically generated during MaxQuant data import',NOW(),'$userID')");
my ($maxExpCondID)=$dbh->selectrow_array("SELECT MAX(ID_EXPCONDITION) FROM EXPCONDITION");
my $sthInsEC=$dbh->prepare("INSERT INTO EXPCONDITION (ID_EXPCONDITION,ID_DESIGN,NAME) VALUES (?,?,?)");
my $sthInsOE=$dbh->prepare("INSERT INTO OBS_EXPCONDITION (ID_EXPCONDITION,ID_OBSERVATION,FRACTION_GROUP) VALUES(?,?,?)");
my $numStates=0;
my $numLabelTypes=scalar keys %distinctLabelings;
foreach my $labelClass (sort{$LABELCLASS_POS{$a} <=> $LABELCLASS_POS{$b}} keys %distinctLabelings) {
	$maxDesignID++;
	my $designID=$maxDesignID;
	my $labelClassStrg=($numLabelTypes==1)? '' : ':'.$labelClass;
	$sthDesign->execute($designID,"MaxQuant$labelClassStrg [$importDate]");
	$labelClass2Design{$labelClass}=$designID;


#my (%label2targetPos,%labelName,@labelList);
#if ($labeling eq 'SILAC') {
#	if ($labelingPlex==2) {
#		%label2targetPos=('L'=>1,'H'=>2);
#		@labelList=('L','H');
#	}
#	else {
#		%label2targetPos=('L'=>1,'M'=>2,'H'=>3);
#		@labelList=('L','M','H');
#	}
#	%labelName=('L'=>'Light','M'=>'Medium','H'=>'Heavy');
#
#}
#elsif ($isobarModifInfo{ID}) { # $labeling=~/ITRAQ|TMT/
#	foreach my $repIdx (0..($labelingPlex-1)) {
#		$label2targetPos{$repIdx}=$repIdx+1;
#		$labelName{$repIdx}=$isobarModifInfo{TAGS}[$repIdx];
#		push @labelList,$repIdx;
#	}
#}
##else {%label2targetPos=('NONE'=>0);}

my %stateName2ID;
	foreach my $paramGr (@{$distinctLabelings{$labelClass}}) {
		if ($labeling{$paramGr} ne 'FREE') {
			#>EXPCONDITION
			foreach my $label (@{$labelList{$paramGr}}) {
#next if ($designInfo{'LABEL_EXPCOND'} && $designInfo{'LABEL_EXPCOND'}{$labelClass} && $designInfo{'LABEL_EXPCOND'}{$labelClass}{$label}); # to prevent state duplication in case multiple group with same labeling
if ($stateName2ID{ $labelName{$paramGr}{$label} }) { # to prevent state duplication in case multiple group with same labeling
	$designInfo{'LABEL_EXPCOND'}{$labelClass}{$label}{'ID'}=$stateName2ID{ $labelName{$paramGr}{$label} };
}
else {
				$maxExpCondID++;
				$sthInsEC->execute($maxExpCondID,$designID,$labelName{$paramGr}{$label});
				$designInfo{'LABEL_EXPCOND'}{$labelClass}{$label}{'ID'} = $stateName2ID{ $labelName{$paramGr}{$label} } = $maxExpCondID;
				$numStates++;
}
			}
		}
	
		my %labelReplicCount;
		foreach my $sampleName (sort{$sampIDs{$a}<=>$sampIDs{$b}} keys %{$designInfo{'SamplesInGroup'}{$paramGr}}) { # $designInfo{'ExpconditionOriginalName'}
			foreach my $expCondName (sort{$designInfo{'Expcondition'}{$a}{'POS'}<=>$designInfo{'Expcondition'}{$b}{'POS'}} keys %{$designInfo{'ExpconditionOriginalName'}{$sampleName}}) { # %{$designInfo{'Expcondition'}}
				my ($label,$targetPos);
				if ($labeling{$paramGr} eq 'FREE') { # Label-free -> each MQ experiment is a state (EXPCONDITION)
					$label='NONE';
					$targetPos=0;
					$maxExpCondID++;
					$designInfo{'Expcondition'}{$expCondName}{'ID'}=$maxExpCondID; # needed for EXPCONDITION_QUANTIF
					$sthInsEC->execute($maxExpCondID,$designID,$expCondName);
					$numStates++;
				}
				else {
					$label=$designInfo{'Expcondition2Label'}{$expCondName};
					$targetPos=$label2targetPos{$paramGr}{$label};
					$designInfo{'Expcondition'}{$expCondName}{'ID'}=$designInfo{'LABEL_EXPCOND'}{$labelClass}{$label}{'ID'}; # needed for EXPCONDITION_QUANTIF
				}
				#$targetPos=($expCondName=~/^([LMH]) /)? $label2targetPos{$1} : 0;
				#my ($label)=($designInfo{'Expcondition2Label'}{$expCondName})? $designInfo{'Expcondition2Label'}{$expCondName} : 'NONE';
				my %fractions;
				if ($designInfo{'Fraction'} && $designInfo{'Fraction'}{$expCondName}) {
					my %replicCount;
					foreach my $fraction (sort{$a<=>$b} keys %{$designInfo{'Fraction'}{$expCondName}}){
						#my $techRep=0;
						foreach my $analysisID (sort{$a<=>$b} keys %{$designInfo{'Fraction'}{$expCondName}{$fraction}}){ # only 1 possible!? (PP 21/11/16)
							#$techRep++;
							$replicCount{$fraction}++;
							$labelReplicCount{$label}{$fraction}++;
#next unless $ana2pepQuantification{$analysisID}; # No peptide quantif data for this analysis -> skip (PP 08/11/16)
							if ($labeling{$paramGr} eq 'FREE') {
								$sthInsOE->execute($maxExpCondID,$ana2Observation{$analysisID}{$label},$replicCount{$fraction});
							}
							else {
								$sthInsOE->execute($designInfo{'LABEL_EXPCOND'}{$labelClass}{$label}{'ID'},$ana2Observation{$analysisID}{$label},$labelReplicCount{$label}{$fraction}); # TODO: check $labelReplicCount{}{}
							}
							push @{$fractions{$replicCount{$fraction}}},"#$ana2Observation{$analysisID}{$label}:#$ana2pepQuantification{$analysisID}:#$analysisID:$targetPos";
						}
					}
				}
				else {
					my $replicNum=0;
					foreach my $analysisID (keys %{$designInfo{'Expcondition'}{$expCondName}{'ANALYSIS'}}) {
						$replicNum++;
						$labelReplicCount{$label}++;
#next unless $ana2pepQuantification{$analysisID}; # No peptide quantif data for this analysis -> skip (PP 08/11/16)
						if ($labeling{$paramGr} eq 'FREE') {
							$sthInsOE->execute($maxExpCondID,$ana2Observation{$analysisID}{$label},undef);
						}
						else {
							$sthInsOE->execute($designInfo{'LABEL_EXPCOND'}{$labelClass}{$label}{'ID'},$ana2Observation{$analysisID}{$label},$labelReplicCount{$label}); # TODO: check $labelReplicCount{}
						}
						push @{$fractions{$replicNum}},"#$ana2Observation{$analysisID}{$label}:#$ana2pepQuantification{$analysisID}:#$analysisID:$targetPos";
					}
				}
				my @replicates;
				foreach my $replic (sort{$a<=>$b} keys %fractions) {
					push @replicates,join('+',@{$fractions{$replic}});
				}
				#$designInfo{'ALL_STATES'}.=';' if $designInfo{'STATES'};
				#$designInfo{'ALL_STATES'}.=(scalar @replicates).','.join('.',@replicates).',#'.$maxExpCondID;
				#push @{$designInfo{'STATES'}{$sampleName}{$expCondName}},(scalar @replicates).','.join('.',@replicates).',#'.$maxExpCondID;
				push @{$designInfo{'STATES'}{$sampleName}},(scalar @replicates).','.join('.',@replicates).',#'.$designInfo{'Expcondition'}{$expCondName}{'ID'};
				push @{$designInfo{'ALL_STATES'}{$labelClass}{$targetPos}},@replicates;
		
			}
		}
	}
}
$sthDesign->finish;
$sthInsEC->finish;
$sthInsOE->finish;

my $designStrg=($numLabelTypes==1)? 'Design' : 'Designs';
my $stateStrg=($numStates==1)? 'State' : 'States';
print LOG " -$numLabelTypes $designStrg and $numStates $stateStrg generated (".strftime("%H:%M:%S %d/%m/%Y",localtime).").\n";

#$dbh->disconnect; close LOG; die "END"; ########## DEBUG

###################################
####>Protein Quantification(s)<####
###################################
# *** NOTE: Global(all samples together) intensities are not recorded: Not useful and require a single global design for all quantifs of all labeling types ***
my (%quantifParamIDs,%requiredParams,@protQuantifList,%protIntQuantifs,%recordedParams);
my ($MQMethID,$ratioQuantMethID);
my $numProtQuantif=0;
#my  $isobarCorrectedStrg=($labeling eq 'TMT')? 'not corrected' : 'corrected'; # for isobaric only. *******WARNING: Could be controlled by a parameter rather than by isobaric label type!!!********
if ($importProtQuantif) {

	####>Fetching list of quantification parameters<####
	($MQMethID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='MQ'");
	my $sthQP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,P.CODE FROM QUANTIFICATION_PARAMETER P,QUANTIFICATION_METHOD M WHERE P.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND M.ID_QUANTIFICATION_METHOD=?");
	$sthQP->execute($MQMethID);
	while (my ($paramID,$code)=$sthQP->fetchrow_array) {
		$quantifParamIDs{$code}=$paramID;
	}
	my $sthInsEQ=$dbh->prepare("INSERT INTO EXPCONDITION_QUANTIF (ID_EXPCONDITION,ID_QUANTIFICATION) VALUES (?,?)");
	my $sthInsParQ=$dbh->prepare("INSERT IGNORE INTO PARENT_QUANTIFICATION (ID_PARENT_QUANTIFICATION,ID_QUANTIFICATION) VALUES (?,?)");

	foreach my $labelClass (sort{$LABELCLASS_POS{$a} <=> $LABELCLASS_POS{$b}} keys %distinctLabelings) {
		my $designID=$labelClass2Design{$labelClass};
		
		my $globalQuantiID;
		foreach my $paramGr (@{$distinctLabelings{$labelClass}}) {
		
			if ($labeling{$paramGr} eq 'SILAC' && $labelingPlex{$paramGr} > 1 && !$ratioQuantMethID) {
				($ratioQuantMethID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='PROT_RATIO_PEP'");
				$sthQP->execute($ratioQuantMethID);
				while (my ($paramID,$code)=$sthQP->fetchrow_array) {
					$code="RATIO:$code" if $code=~/PEPTIDES|RAZ_UNI_PEP|UNIQUE_PEP/; # to distinguish from same codes used in MQ method
					$quantifParamIDs{$code}=$paramID;
				}
			}

			my %quantifParamCodes=(
				'MQ_INT'=>'Intensity','MQ_IBAQ'=>'iBAQ','MQ_LFQ'=>'LFQ intensity','MQ_SC'=>'MS/MS Count',
				'ISOBAR:MQ_INT'=>"Reporter intensity $isobarCorrectedStrg",
				'RATIO_RAW'=>'Ratio #%#','RATIO'=>'Ratio #%# normalized','RATIO_VAR'=>'Ratio #%# variability [%]',
				'PEPTIDES'=>'Peptides','RAZ_UNI_PEP'=>'Razor + unique peptides','UNIQUE_PEP'=>'Unique peptides',
				'NUM_PEP_USED'=>'Ratio #%# count'
			);
			#my @noRatioParamCodes=('MQ_INT','MQ_IBAQ','MQ_LFQ','MQ_SC','PEPTIDES','RAZ_UNI_PEP','UNIQUE_PEP');
			#my @ratioParamCodes=('RATIO_RAW','RATIO_NORM','RATIO_VAR','PEPTIDES','RAZ_UNI_PEP','UNIQUE_PEP','NUM_PEP_USED');
			my @stateIntensParamCodes=($isobarModifInfo{$paramGr})? ('ISOBAR:MQ_INT') : ('MQ_INT','MQ_IBAQ','MQ_LFQ','MQ_SC'); # MQ_INT must be 1st -> $measStrg
			my @ratioParamCodes=('RATIO_RAW','RATIO','RATIO_VAR','NUM_PEP_USED');
			my @allStatesParamCodes=('PEPTIDES','RAZ_UNI_PEP','UNIQUE_PEP');

			####
			####>For each sample: Create a non-ratio quantification (Intensity,iBAQ,LFQ,SpecCount) containing all (non-)label channels<####
			####
			foreach my $sampleName (sort{$sampIDs{$a}<=>$sampIDs{$b}} keys %{$designInfo{'SamplesInGroup'}{$paramGr}}) { # eg. label '123','456'  $designInfo{'ExpconditionOriginalName'}}
				my $sampleNameStrg=($sampleName eq 'No experiment')? '' : " $sampleName";
				$numProtQuantif++;
				my $quantifAnnot="LABEL=$labeling{$paramGr}";
				$quantifAnnot.='::SOFTWARE=MQ'.$versionStrg.'::PEPTIDES='.$peptideFilterStrg.'::RATIO_TYPE=None::STATES='.join(';',@{$designInfo{'STATES'}{$sampleName}});
				my $quantifName=($numSamples==1)? $sampleName : "$sampleName [Intensities]";
				$sthInsQ->execute($MQMethID,$designID,$quantifName,'protein',$quantifAnnot) || die $sthInsQ->errstr();
				my $quantifID=$dbh->last_insert_id(undef,undef,'QUANTIFICATION','ID_QUANTIFICATION');
				push @protQuantifList,$quantifID;
				$protIntQuantifs{$quantifID}=1;
				$sthInsParQ->execute($peptideQuantifID{$paramGr},$quantifID) || die $sthInsParQ->errstr();
				my %linkedAna;
				my $targetPos=0;
				foreach my $expCondName (sort{$designInfo{'Expcondition'}{$a}{'POS'}<=>$designInfo{'Expcondition'}{$b}{'POS'}} keys %{$designInfo{'ExpconditionOriginalName'}{$sampleName}} ) { # label: ' L 123','H 123',... label-free: '$sampleName'
					$targetPos++;
					my $expCondNameStrg=($expCondName eq 'No experiment')? '' : " $expCondName";
					#$quantifAnnot='LABEL=FREE::SOFTWARE=MQ'.$versionStrg.'::RATIO_TYPE=None::STATE='.$designInfo{'Expcondition'}{$expCondName}{'STATE'};
					#$sthinsQ->execute($MQMethID,$designID,$expCondName,'protein',$quantifAnnot) || die $sthinsQ->errstr();
					#my $quantifID=$dbh->last_insert_id(undef,undef,'QUANTIFICATION','ID_QUANTIFICATION');
					$sthInsEQ->execute($designInfo{'Expcondition'}{$expCondName}{'ID'},$quantifID) || die $sthInsEQ->errstr();
		
					if ($designInfo{'Fraction'} && $designInfo{'Fraction'}{$expCondName}) {
						foreach my $fraction (keys %{$designInfo{'Fraction'}{$expCondName}}){
							foreach my $analysisID (keys %{$designInfo{'Fraction'}{$expCondName}{$fraction}}) {
#next unless $ana2pepQuantification{$analysisID}; # No peptide quantif data for this analysis -> skip (PP 08/11/16)
								next if $linkedAna{$analysisID}; # do only once
								$sthInsAQ->execute($quantifID,$analysisID) || die $sthInsAQ->errstr();
								#$sthInsParQ->execute($ana2pepQuantification{$analysisID},$quantifID) || die $sthInsParQ->errstr();
								$linkedAna{$analysisID}=1;
							}
						}
					}
					else {
						#foreach my $analysisID (keys %{$designInfo{'Experiment'}{$expCondName}}){
						foreach my $analysisID (keys %{$designInfo{'Expcondition'}{$expCondName}{'ANALYSIS'}}) {
#next unless $ana2pepQuantification{$analysisID}; # No peptide quantif data for this analysis -> skip (PP 08/11/16)
							next if $linkedAna{$analysisID}; # do only once
							$sthInsAQ->execute($quantifID,$analysisID) || die $sthInsAQ->errstr();
							#$sthInsParQ->execute($ana2pepQuantification{$analysisID},$quantifID) || die $sthInsParQ->errstr();
							$linkedAna{$analysisID}=1;
						}
					}
		
					###>Make parameters state-specific: channel-spec for label ('Intensity L 123'), sample-spec for label-free ('Intensity')
					foreach my $paramCode (@stateIntensParamCodes) {
						#my $usedParamName=($isobarModifInfo{ID} && $paramCode eq 'MQ_INT')? "Reporter intensity $isobarCorrectedStrg" : $quantifParamCodes{$paramCode};
						(my $trueParamCode=$paramCode)=~s/ISOBAR://; # eg. ISOBAR:CODE -> CODE
						push @{$recordedParams{"$quantifParamCodes{$paramCode}$expCondNameStrg"}},[$quantifID,$trueParamCode,$targetPos]; # push because same param used for different quantifs (in case of ratios)
					}
					#$requiredParams{$quantifID}{$targetPos}=($isobarModifInfo{ID})? "Reporter intensity $isobarCorrectedStrg".$expCondNameStrg : 'Intensity'.$expCondNameStrg;
					$requiredParams{$quantifID}{$targetPos}=($isobarModifInfo{$paramGr})? $quantifParamCodes{'MQ_INT'}.$sampleNameStrg : $quantifParamCodes{'MQ_INT'}.$expCondNameStrg;
				}
		
				##>Make peptide params sample-specific: These peptide params not state-specific (global) for sample-level quantif
				foreach my $paramCode (@allStatesParamCodes) {
					push @{$recordedParams{"$quantifParamCodes{$paramCode}$sampleNameStrg"}},[$quantifID,$paramCode,0]; # push because same param used for different quantifs (in case of ratios)
				}
		
			}

			
			####
			####>Add a set of quantifications to keep labeling information (Ratio value)<####
			####
			if ($labeling{$paramGr} eq 'SILAC' && $labelingPlex{$paramGr} > 1) { # Add Ratio (No ratios for iTRAQ!!!)
		
				my $baseQuantifAnnot='LABEL='.$labeling{$paramGr}.'::SOFTWARE=MQ'.$versionStrg.'::BIAS_CORRECTION=TRUE::NORMALIZATION_METHOD=global.normalization.median::PEPTIDES='.$peptideFilterStrg.'::RATIO_TYPE=SimpleRatio';
		
				#my $globalQuantiID;
				if ($numSamples > 1 && !$globalQuantiID) { # -> Global ratio(s) != sample-specific ratio(s)
					my $quantifAnnot=$baseQuantifAnnot;
					#>RATIOS
					my $ratioStrg='';
					foreach my $refIdx (0..$#{$labelList{$paramGr}}-1) {
						foreach my $testIdx (1..$#{$labelList{$paramGr}}) {
							next if $testIdx==$refIdx;
							$ratioStrg.=';' if $ratioStrg;
							$ratioStrg.='#'.$designInfo{'LABEL_EXPCOND'}{$labelClass}{ $labelList{$paramGr}[$testIdx] }{'ID'}.'/#'.$designInfo{'LABEL_EXPCOND'}{$labelClass}{ $labelList{$paramGr}[$refIdx] }{'ID'};
						}
					}
					$quantifAnnot.='::RATIOS='.$ratioStrg;
					#>STATES
					my $labelStateStrg='';
					#$globalQuantiID=$quantifID+$numSamples+1; # use highest DB ID
					$numProtQuantif++;
					foreach my $targetPos (sort{$a<=>$b} keys %{$designInfo{'ALL_STATES'}{$labelClass}}) {
						my $label=$labelList{$paramGr}[$targetPos-1];
next unless $label;
						$labelStateStrg.=';' if $labelStateStrg;
						$labelStateStrg.=(scalar @{$designInfo{'ALL_STATES'}{$labelClass}{$targetPos}}).','.join('.',@{$designInfo{'ALL_STATES'}{$labelClass}{$targetPos}}).',#'.$designInfo{'LABEL_EXPCOND'}{$labelClass}{$label}{'ID'};
					}
					$quantifAnnot.='::STATES='.$labelStateStrg;
					$sthInsQ->execute($ratioQuantMethID,$designID,"Global ratios [$importDate]",'protein',$quantifAnnot) || die $sthInsQ->errstr();
					$globalQuantiID=$dbh->last_insert_id(undef,undef,'QUANTIFICATION','ID_QUANTIFICATION');
					#>EXPCONDITION_QUANTIF
					foreach my $label (keys %{$designInfo{'LABEL_EXPCOND'}{$labelClass}}) {
						$sthInsEQ->execute($designInfo{'LABEL_EXPCOND'}{$labelClass}{$label}{'ID'},$globalQuantiID) || die $sthInsEQ->errstr();
					}
					$sthInsParQ->execute($peptideQuantifID{$paramGr},$globalQuantiID) || die $sthInsParQ->errstr();
					push @protQuantifList,$globalQuantiID;
				}
		
				my @allStatesAnnot;
				foreach my $sampleName (sort{$sampIDs{$a}<=>$sampIDs{$b}} keys %{$designInfo{'SamplesInGroup'}{$paramGr}}) { # $designInfo{'ExpconditionOriginalName'}
					my $sampleNameStrg=($sampleName eq 'No experiment')? '' : " $sampleName";
		
					my $quantifAnnot=$baseQuantifAnnot;
					#>RATIOS
					my $ratioStrg='';
					foreach my $refIdx (0..$#{$labelList{$paramGr}}-1) {
						my $refExpCondName=$labelList{$paramGr}[$refIdx].$sampleNameStrg;
						foreach my $testIdx (1..$#{$labelList{$paramGr}}) {
							next if $testIdx==$refIdx;
							my $testExpCondName=$labelList{$paramGr}[$testIdx].$sampleNameStrg;
							$ratioStrg.=';' if $ratioStrg;
							$ratioStrg.='#'.$designInfo{'Expcondition'}{$testExpCondName}{'ID'}.'/#'.$designInfo{'Expcondition'}{$refExpCondName}{'ID'};
						}
					}
					$quantifAnnot.='::RATIOS='.$ratioStrg;
					#>STATES
					$quantifAnnot.='::STATES='.join(';',@{$designInfo{'STATES'}{$sampleName}});
		
					$numProtQuantif++;
					$sthInsQ->execute($ratioQuantMethID,$designID,"$sampleName [Ratios]",'protein',$quantifAnnot) || die $sthInsQ->errstr();
					my $quantifID=$dbh->last_insert_id(undef,undef,'QUANTIFICATION','ID_QUANTIFICATION');
					$sthInsParQ->execute($peptideQuantifID{$paramGr},$quantifID) || die $sthInsParQ->errstr();
		
					push @protQuantifList,$quantifID;
		
					my %linkedAna;
					foreach my $expCondName ( keys  %{$designInfo{'ExpconditionOriginalName'}{$sampleName}} ) {
						$sthInsEQ->execute($designInfo{'Expcondition'}{$expCondName}{'ID'},$quantifID) || die $sthInsEQ->errstr();
						if ($designInfo{'Fraction'} && $designInfo{'Fraction'}{$expCondName}) {
							foreach my $fraction (keys %{$designInfo{'Fraction'}{$expCondName}}){
								foreach my $analysisID (keys %{$designInfo{'Fraction'}{$expCondName}{$fraction}}){
#next unless $ana2pepQuantification{$analysisID}; # No peptide quantif data for this analysis -> skip (PP 08/11/16)
									next if $linkedAna{$analysisID}; # do only once
									$sthInsAQ->execute($quantifID,$analysisID) || die $sthInsAQ->errstr();
									#$sthInsParQ->execute($ana2pepQuantification{$analysisID},$quantifID) || die $sthInsParQ->errstr();
									if ($globalQuantiID) {
										$sthInsAQ->execute($globalQuantiID,$analysisID) || die $sthInsAQ->errstr();
										#$sthInsParQ->execute($ana2pepQuantification{$analysisID},$globalQuantiID) || die $sthInsParQ->errstr();
									}
									$linkedAna{$analysisID}=1;
								}
							}
						}
						else {
							#foreach my $analysisID (keys %{$designInfo{'Experiment'}{$expCondName}}){ #}
							foreach my $analysisID (keys %{$designInfo{'Expcondition'}{$expCondName}{'ANALYSIS'}}) {
#next unless $ana2pepQuantification{$analysisID}; # No peptide quantif data for this analysis -> skip (PP 08/11/16)
								next if $linkedAna{$analysisID}; # do only once
								$sthInsAQ->execute($quantifID,$analysisID) || die $sthInsAQ->errstr();
								$sthInsParQ->execute($ana2pepQuantification{$analysisID},$quantifID) || die $sthInsParQ->errstr();
								if ($globalQuantiID) {
									$sthInsAQ->execute($globalQuantiID,$analysisID) || die $sthInsAQ->errstr();
									#$sthInsParQ->execute($ana2pepQuantification{$analysisID},$globalQuantiID) || die $sthInsParQ->errstr();
								}
								$linkedAna{$analysisID}=1;
							}
						}
					}
		
					###> Make parameters ratio-specific (eg H/L)
					my $targetPos=0;
					foreach my $refIdx (0..$#{$labelList{$paramGr}}-1) {
						foreach my $testIdx (1..$#{$labelList{$paramGr}}) {
							next if $testIdx==$refIdx;
							$targetPos++;
							my $ratioName=$labelList{$paramGr}[$testIdx].'/'.$labelList{$paramGr}[$refIdx];
							foreach my $paramCode (@ratioParamCodes) {
								my $paramText="$quantifParamCodes{$paramCode}$sampleNameStrg";
								$paramText=~s/#%#/$ratioName/;
								push @{$recordedParams{$paramText}},[$quantifID,$paramCode,$targetPos]; # push because same param used for different quantifs (in case of ratios)
							}
							$requiredParams{$quantifID}{$targetPos}="Ratio $ratioName normalized";
						}
					}
					##>These peptide counts are not channel-specific
					foreach my $paramCode (@allStatesParamCodes) { # "RATIO:$paramCode" to distinguish PROT_RATIO_PEP parameters from MQ parameters
						push @{$recordedParams{"$quantifParamCodes{$paramCode}$sampleNameStrg"}},[$quantifID,"RATIO:$paramCode",0]; # push because same param used for different quantifs (in case of ratios)
					}
				}
		
				####>Global ratio quantif
				if ($globalQuantiID) {
					my $targetPos=0;
					foreach my $refIdx (0..$#{$labelList{$paramGr}}-1) {
						foreach my $testIdx (1..$#{$labelList{$paramGr}}) {
							next if $testIdx==$refIdx;
							$targetPos++;
							my $ratioName=$labelList{$paramGr}[$testIdx].'/'.$labelList{$paramGr}[$refIdx];
							foreach my $paramCode (@ratioParamCodes) {
								my $paramText=$quantifParamCodes{$paramCode};
								$paramText=~s/#%#/$ratioName/;
								push @{$recordedParams{$paramText}},[$globalQuantiID,$paramCode,$targetPos]; # push because same param used for different quantifs (in case of ratios)
							}
							$requiredParams{$globalQuantiID}{$targetPos}="Ratio $ratioName normalized";
						}
					}
					##>These peptide counts are not channel-specific
					foreach my $paramCode (@allStatesParamCodes) {
						push @{$recordedParams{$quantifParamCodes{$paramCode}}},[$globalQuantiID,"RATIO:$paramCode",0]; # push because same param used for different quantifs (in case of ratios)
					}
				}
			}
		}
	}
	$sthQP->finish;
	$sthInsEQ->finish;
	$sthInsParQ->finish;
	$sthInsEQ->finish;
	$sthInsParQ->finish;

}
$sthInsAQ->finish;
$sthInsQ->finish;
$dbh->commit;
#$dbh->rollback; # DEBUG
$dbh->disconnect; # Release connection

print LOG " -",(scalar @designGroups)," peptide and $numProtQuantif protein Quantifications to be imported.\n";
#print "Done.</FONT><BR>\n";


#exit; # DEBUG

########################
####>msmms.txt file<####
########################
my $counter=0;
my %msmsMissedCut; # missed cleavage data missing in evidence.txt before MaxQuant v~1.5 but exists in msms.txt
if (-e "$tmpFilesDir/$msmsFile") {
	&addToFileStat("2/$numSteps Reading msms file\n");
	print LOG "2/$numSteps Reading msms file...";
	open (MSMS,"$tmpFilesDir/$msmsFile") || die "Unable to open $tmpFilesDir/$msmsFile";
	my $firstmsmsLine=<MSMS>;
	my @parameters=split(/ *\t */,$firstmsmsLine);
	my %msmsColNum;
	my $ncol=0;
	foreach my $colName (@parameters) {
		$msmsColNum{$colName}=$ncol;
		$ncol++;
	}
	###> Copy the msmsFile in all peptide data directories
	foreach my $rawName (keys %anaInfo) {
		my $anaID=$anaInfo{$rawName}{'ID_ANALYSIS'};
		next unless $anaID;
		open($anaInfo{$rawName}{'INFILE'},'>',"$promsPath{peptide}/proj_$projectID/ana_$anaID/msms.txt") or die;
		my $fhavalue=$anaInfo{$rawName}{'INFILE'};
		print $fhavalue $firstmsmsLine;
	}
	while (<MSMS>) {
		$counter++;
		my @infos=split(/ *\t */,$_); # remove starting/trailing spaces
		#$infos[0]=~s/\s+$//; # remove trailing spaces
		my $fhavalue=$anaInfo{$infos[0]}{'INFILE'};
		$msmsMissedCut{$infos[$msmsColNum{'Sequence'}]}=$infos[$msmsColNum{'Missed cleavages'}];
		print $fhavalue $_;
		if ($counter > 10000) {
			print LOG '.';
			$counter=0;
		}
	}
	close MSMS;
	foreach my $rawName (keys %anaInfo) {
		my $fhavalue=$anaInfo{$rawName}{'INFILE'};
		next unless $fhavalue;
		close($fhavalue);
		delete($anaInfo{$rawName}{'INFILE'})
	}
	print LOG " Done (".strftime("%H:%M:%S %d/%m/%Y",localtime).").\n";
}


###########################
####>evidence.txt file<####
###########################
###> TODO : check DECOY + FDR + MAX_RANK for ANALYSIS
&addToFileStat("3/$numSteps Reading evidence file\n");
print LOG "3/$numSteps Reading evidence file...";

###>Filter "Only identified by site"
my %onlyBySite;
&filterOnlyIdentifiedBySite(\%onlyBySite); # uses proteinGroups.txt

open (EVIDENCE,"$tmpFilesDir/$evidenceFile")  || die "Unable to open $tmpFilesDir/$evidenceFile";
#my @evidenceColumns=('Sequence','Length','Modifications','Modified sequence','Missed cleavages','Proteins','Leading proteins','Leading razor protein',
#                     'Raw file','Charge','m/z','Mass','Mass Error [ppm]','Mass Error [Da]','Retention time','Calibrated retention time',
#                     'PEP','MS/MS Count','MS/MS Scan Number','Score','Intensity','Peptide ID');
my (%evidenceColNum,%pepInfo,%queryNum,%pepVmod);
my (%evidence2peptide,%maxProtMatch,%actualSeq2seq,%protDbRank,%matchList,%bestScore,%proteinRank); #,%bestScorebad ## ,%bestQuery
my (%seq2posString,%peptideXIC,%isGhost,%featureCount);
my %raw2UsedIdentifiers;
my %skipPeptide; # in case Contaminants are excluded

my %residue2varMod;
foreach my $anaName (keys %varMods) {
	foreach my $refVarMod (@{$varMods{$anaName}}) {
		foreach my $res (split(',',$refVarMod->[1])) {
			$residue2varMod{$res}=$refVarMod->[2]; # assumes only 1 varMod for a given residue
		}
	}
}
my $missedCutData=1;
$counter=0;
while (my $line=<EVIDENCE>) {
	$counter++;
	$line=~s/\s+$//; # chomp is not enough <- Windows
	my @parameters=split(/ *\t */,$line); # remove starting/trailing spaces
	if ($.==1) { # 1st line of the file
		my $ncol=0;
		foreach my $colName (@parameters) {
			$evidenceColNum{uc($colName)}=$ncol; # Upper case for all column headers!!!
			$ncol++;
		}
		$missedCutData=0 unless $evidenceColNum{uc('Missed cleavages')};
		next;
	}
	my $anaName=$parameters[$evidenceColNum{uc('Raw file')}];
	next unless $anaInfo{$anaName};
	my $analysisID=$anaInfo{$anaName}{'ID_ANALYSIS'};
	my $paramGr=$ana2ParamGroup{$anaName};
	$queryNum{$analysisID}++;
	###> $parameters[$evidenceColNum{'Proteins'}] can be empty if only one REV__ protein is matched !
	my $actualSequence=$parameters[$evidenceColNum{uc('Modified sequence')}];
	my $bearSequence=$parameters[$evidenceColNum{uc('Sequence')}];
	my ($matchGood,$matchBad)=(0,0); #,$matchCON,0
	my $proteinList=($parameters[$evidenceColNum{uc('Proteins')}])? $parameters[$evidenceColNum{uc('Proteins')}] : $parameters[$evidenceColNum{uc('Leading proteins')}];
#my $skipThisPeptide=0; # local var for speed
#if ($excludeCON && $proteinList=~/^CON__/) { # Leading protein!
#	$skipPeptide{$bearSequence}=1;
#	$skipThisPeptide=1;
#}
#my $protRank=1;
	my ($numTrueGood,$numCon)=(0,0);
	foreach my $rawIdentifier (split(/;/,$proteinList)) {
		last if $onlyBySite{$rawIdentifier};
		my $identifier=$rawIdentifier; # default
		$raw2UsedIdentifiers{$rawIdentifier}=$rawIdentifier; # default
		if ($identifier =~ /REV__/) {
			#$identifier =~ s/REV__/DECOY_/g;
			$matchBad=1;
		}
		#elsif ($excludeCON && $identifier=~/CON__/) {
		#	$matchCON=1;
		#}	
		#elsif ($skipThisPeptide) {
		#	$matchGood=1;
		#}
		elsif ($excludeCON && $identifier=~/CON__/) {
			$numCon++;
			$matchGood=1;
		}
		else {
			$numTrueGood++;
			$matchGood=1;
			my $tmpIdenfitier;
			foreach my $dbID (@databankIDs) {
				($tmpIdenfitier)=($rawIdentifier=~/$idParseRules{$dbID}/);
				if ($tmpIdenfitier) {
					$identifier=$raw2UsedIdentifiers{$rawIdentifier}=$tmpIdenfitier;
					last;
				}
			}

#$protDbRank{$identifier}=($identifier=~/CON__/ && $contaminantDB)? 2 : 1;
#$matchList{$protDbRank{$identifier}}{$identifier}{$analysisID}=1;
			$matchList{$identifier}{$analysisID}=1;
			$maxProtMatch{$analysisID}{$identifier}{$actualSequence}++;

#$proteinRank{$identifier}=$protRank if (!$proteinRank{$identifier} || $proteinRank{$identifier} > $protRank); # needed to select best MQ prot during MG re-creation (pseudo MG! Better to read proteinGroups.txt)
			unless ($proteinRank{$identifier}) {
				$proteinRank{$identifier}=($rawIdentifier eq $parameters[$evidenceColNum{uc('Leading Razor Protein')}])? 1 : 2; # needed to select best MQ prot during MG re-creation
			}
		}
#$protRank++;
#$protDbRank{$identifier}=($identifier=~/CON__/ && $contaminantDB)? 2 : 1;
#$matchList{$protDbRank{$identifier}}{$identifier}{$analysisID}=1;
#$maxProtMatch{$analysisID}{$identifier}{$actualSequence}++;
	}
	if ($matchBad) {
		$featureCount{'DECOY'}{$analysisID}{$actualSequence}++; # for FDR
#$bestScorebad{$analysisID}{$actualSequence}=$score if ($score && $bestScorebad{$analysisID} || !$bestScorebad{$analysisID}{$actualSequence});
#$bestScorebad{$analysisID}{$actualSequence}=$score if (!defined($bestScorebad{$analysisID}{$actualSequence}) || ($bestScorebad{$analysisID}{$actualSequence} && $score>$bestScorebad{$analysisID}{$actualSequence}));
# NOT NEEDED # @{$pepInfo{$analysisID}{"-$queryNum{$analysisID}"}}=($actualSequence,$parameters[$evidenceColNum{uc('Sequence')}],$parameters[$evidenceColNum{uc('Length')}],"-$queryNum{$analysisID}",$score,$missedCut,$massExp,$parameters[$evidenceColNum{uc('Mass')}],$parameters[$evidenceColNum{uc('m/z')}],$massErrorDa,$charge,$etString,$validStatus,"PEP=$pep",$parameters[$evidenceColNum{uc('MS/MS Count')}]);
# NOT NEEDED # @{$evidence2peptide{$parameters[$evidenceColNum{uc('id')}]}}=($analysisID,"-$queryNum{$analysisID},$actualSequence") unless $evidence2peptide{$actualSequence}{$parameters[$evidenceColNum{uc('id')}]};
#print "Bad '$actualSequence' $analysisID -$queryNum{$analysisID} score$bestScorebad{$analysisID}{$actualSequence}<BR>\n";
	}
	#elsif ($matchCON && !$matchGood) {
	#	$featureCount{'TARGET'}{$analysisID}{$actualSequence}++; # for FDR
	#}
	next unless $matchGood;
	#if ($matchGood) {
	$featureCount{'TARGET'}{$analysisID}{$actualSequence}++; # for FDR
#next if $skipThisPeptide;
if ($excludeCON && $numCon && !$numTrueGood) { # matches only CON__
	$skipPeptide{$bearSequence}=1;
	next;
}	
	#>Missed Cleavages
	my $missedCut=0;
	if ($missedCutData) {$missedCut=$parameters[$evidenceColNum{uc('Missed cleavages')}];}
	elsif ($msmsMissedCut{$bearSequence}) {$missedCut=$msmsMissedCut{$bearSequence};}
	#>Scores
	my ($score,$pepDataStrg);
	if ($parameters[$evidenceColNum{uc('Score')}] =~ /\d/) {
		$score=$parameters[$evidenceColNum{uc('Score')}]*1;
		$pepDataStrg='DELTA_SC='.($parameters[$evidenceColNum{uc('Delta score')}]*1);
		$pepDataStrg.='##PEP='.$parameters[$evidenceColNum{'PEP'}] if $parameters[$evidenceColNum{'PEP'}] =~ /\d/;
		$pepDataStrg.='##'; # for MQ_EVI below
	}
	else {$score=0;}
	$pepDataStrg.='MQ_EVI='.$parameters[$evidenceColNum{uc('id')}];
	my $charge=$parameters[$evidenceColNum{uc('Charge')}];
	#$queryNum{$analysisID}++;
	my $massExp=($parameters[$evidenceColNum{uc('m/z')}]-1.007825032)*$charge; # same way computed in storeAnalyses.cgi for Paragon
	my $massErrorDa=($evidenceColNum{uc('Mass Error [Da]')})? $parameters[$evidenceColNum{uc('Mass Error [Da]')}] : ($parameters[$evidenceColNum{uc('Mass Error [ppm]')}])? $parameters[$evidenceColNum{uc('Mass Error [ppm]')}]*$parameters[$evidenceColNum{uc('Mass')}]/1000000 : undef;
	$massErrorDa=undef if ($massErrorDa && $massErrorDa !~ /[^-\d\.]/); # eq 'NaN'
	#>Check for isobaric modif (Not listed as variable modif!!!!)
	if ($isobarModifInfo{$paramGr} && $isobarModifInfo{$paramGr}{MONO_MASS} && $massErrorDa && $massErrorDa > 10) {
		my $numModif=int(0.5+($massErrorDa*$charge/$isobarModifInfo{$paramGr}{MONO_MASS}));
		if ($numModif >= 1) {
			my ($Nterm,$protNterm,@positions);
			my $massShift=($isobarModifInfo{$paramGr}{MONO_MASS}/$charge); # +1.007825032;
			foreach my $res (sort keys %{$isobarModifInfo{$paramGr}{RES}}) {
				if ($res eq '=') {
					$Nterm=1;
					$numModif--;
					$massErrorDa-=$massShift;
				}
				elsif ($res eq '-') {
					$protNterm=1;
					$numModif--;
					$massErrorDa-=$massShift;
				}
				else {
					while ($bearSequence =~ /$res/g && $numModif > 0) { # 1st come, 1st served !
						push @positions,$-[0]+1;
						$massErrorDa-=$massShift;
						$numModif--;
					}
				}
				last if $numModif==0;
			}
			@positions=sort{$a<=>$b} @positions;
			if ($Nterm) {unshift @positions,'=';} # 1rst element
			elsif ($protNterm) {unshift @positions,'-';} # assumes N-term and Protein N-term are incompatible
			$isobarModifInfo{$paramGr}{ANA}{$analysisID}{ $queryNum{$analysisID} }=join('.', @positions);
		}
	}

	my $etString=($parameters[$evidenceColNum{uc('MS/MS Scan Number')}] =~ /\d/)? "et$parameters[$evidenceColNum{'CALIBRATED RETENTION TIME'}];sc$parameters[$evidenceColNum{'MS/MS SCAN NUMBER'}];":"et$parameters[$evidenceColNum{'CALIBRATED RETENTION TIME'}];";
	$actualSeq2seq{$actualSequence}=$bearSequence;
	#if ($matchGood) {#}
	##my $validStatus=0;
	(my $validStatus,$isGhost{$analysisID}{ $queryNum{$analysisID} })=($score)? (1,1) : (0,0);
	if (!$bestScore{$analysisID} || !defined($bestScore{$analysisID}{$actualSequence}) || $bestScore{$analysisID}{$actualSequence} < $score) {
		##$validStatus=1 if $score;
		##if ($bestScore{$analysisID}{$actualSequence}) { # > 0 (also $score > 0) Set previous best as ghost!
		##	my $prevBestQueryNum=$bestQuery{$analysisID}{$actualSequence};
		##	$isGhost{$analysisID}{$prevBestQueryNum}=1;
		##	$pepInfo{$analysisID}{$prevBestQueryNum}[12]=0; # validStatus
		##}
		$bestScore{$analysisID}{$actualSequence}=$score; # ranked by descending scores in file <= NOT TRUE (PP 25/10/19)
		##$bestQuery{$analysisID}{$actualSequence}=$queryNum{$analysisID};
	}
	##$isGhost{$analysisID}{ $queryNum{$analysisID} }=1 unless $validStatus;
	#$validStatus=1 if ($score && !$bestScore{$analysisID}{$actualSequence}); # if there is no score, then, it is a ghost peptide
	@{$pepInfo{$analysisID}{ $queryNum{$analysisID} }}=($actualSequence,$bearSequence,$parameters[$evidenceColNum{uc('Length')}],$queryNum{$analysisID},$score,$missedCut,$massExp,$parameters[$evidenceColNum{uc('Mass')}],$parameters[$evidenceColNum{uc('m/z')}],$massErrorDa,$charge,$etString,$validStatus,$pepDataStrg); # ,$parameters[$evidenceColNum{uc('MS/MS Count')}]
	@{$evidence2peptide{$parameters[$evidenceColNum{uc('id')}]}}=($analysisID,$queryNum{$analysisID},$actualSequence);
	if ($labeling{$paramGr} eq 'FREE') {
		my $eviColIdx=($distinctLabelings{'SILAC'})? $evidenceColNum{uc('Intensity L')} : $evidenceColNum{uc('Intensity')}; # data are in 'L' channel in evidence.txt if SILAC in another paramGr!!!!
		$peptideXIC{$analysisID}{ $queryNum{$analysisID} }{'XIC'}{'NONE'}=$parameters[$eviColIdx] if ($eviColIdx && $parameters[$eviColIdx] && $parameters[$eviColIdx] !~ /^(NaN|0)$/);
#print LOG "**$paramGr-$eviColIdx) $bearSequence:";
#print LOG " $parameters[$eviColIdx]\n" if $parameters[$eviColIdx];
	}
	elsif ($labeling{$paramGr} eq 'SILAC') {
		foreach my $colLabel (@{$labelList{$paramGr}}) { # L,M,H
			my $eviColIdx=($labelingPlex{$paramGr}==1)? $evidenceColNum{uc('Intensity L')} : $evidenceColNum{uc("Intensity $colLabel")}; # 'L' channel is used in evidence.txt if only 1-plex!!!!
			$peptideXIC{$analysisID}{ $queryNum{$analysisID} }{'XIC'}{$colLabel}=$parameters[$eviColIdx] if ($eviColIdx && $parameters[$eviColIdx] && $parameters[$eviColIdx] !~ /^(NaN|0)$/);
		}
	}
	elsif ($isobarModifInfo{$paramGr}) { # $labeling=~/ITRAQ|TMT/
		my $posShift=($useReporterIndex)? 0 : 1;
		foreach my $repIdx (@{$labelList{$paramGr}}) { # 0,1,2,...,n
			my $repTag=$repIdx+$posShift;
			my $eviColIdx=$evidenceColNum{uc("Reporter intensity $isobarCorrectedStrg $repTag")};
			$peptideXIC{$analysisID}{ $queryNum{$analysisID} }{'XIC'}{$repIdx}=$parameters[$eviColIdx] if ($eviColIdx && $parameters[$eviColIdx] && $parameters[$eviColIdx] !~ /^(NaN|0)$/);
		}
		#if ($useReporterIndex) {
		#	foreach my $repIdx (@{$labelList{$paramGr}}) { # 0,1,2,...,n
		#		$peptideXIC{$analysisID}{ $queryNum{$analysisID} }{'XIC'}{$repIdx}=$parameters[$evidenceColNum{uc("Reporter intensity $isobarCorrectedStrg $repIdx")}] if ($evidenceColNum{uc("Reporter intensity $isobarCorrectedStrg $repIdx")} && $parameters[$evidenceColNum{uc("Reporter intensity $isobarCorrectedStrg $repIdx")}] !~ /^(NaN|0)$/);
		#	}
		#}
		#else { # use Position
		#	foreach my $repPos (1..$labelingPlex{$paramGr}) {
		#		$peptideXIC{$analysisID}{ $queryNum{$analysisID} }{'XIC'}{$repPos-1}=$parameters[$evidenceColNum{uc("Reporter intensity $isobarCorrectedStrg $repPos")}] if ($evidenceColNum{uc("Reporter intensity $isobarCorrectedStrg $repPos")} && $parameters[$evidenceColNum{uc("Reporter intensity $isobarCorrectedStrg $repPos")}] !~ /^(NaN|0)$/);
		#	}
		#}
	}
	#$featureCount{'TARGET'}{$analysisID}{$actualSequence}++;
	#print "Good '$actualSequence' $analysisID $queryNum{$analysisID} $parameters[$evidenceColNum{'id'}] score=$bestScore{$analysisID}{$actualSequence}<BR>\n";
	#}
	foreach my $vmod (keys %{$anaVmods{$analysisID}}) { # non-label vmods
		if ($parameters[$evidenceColNum{uc($vmod)}] && $evidenceColNum{uc($vmod)} >= 1) {
			@{$pepVmod{$analysisID}{$queryNum{$analysisID}}{$vmod}}=(undef); # default
			my $nbMod=$parameters[$evidenceColNum{uc($vmod)}];
			if ($evidenceColNum{uc("$vmod Probabilities")} && $parameters[$evidenceColNum{uc("$vmod Probabilities")}]) {
				my $probSequence=$parameters[$evidenceColNum{uc("$vmod Probabilities")}]; # "C(0.837)DPRLGKYMAC(0.16)C(0.003)LLYR"
				my @probValues= $probSequence =~ /\(([^\)]+)\)/g; # @probValues = ((0.837) , (0.16) , (0.003))
				my @slices=split(/\([^\)]+\)/,$probSequence);
				my $pos=0;
				my (%probMod,@positions);
				my $probString='##PRB_MQ=';
				my $keepProb=0;
				for (my $i=0; $i <= $#slices; $i++) { # "<=" Change on 2017/12/06 because of GlyGly modification with probability on last residue !!!
					last if $i > $#probValues;# Be careful, 2 cases for GlyGly: 'AAAAAAALQAK(1)SDEK(1)AAVAGK' or 'AAAAAAALQAK(0.714)SDEK(0.286)'
					$pos+=length($slices[$i]);
					$probMod{$pos}=$probValues[$i];
					$probString.=',' if $i > 0;
					$probString.=$pos.':'.$probValues[$i];
					$keepProb=1 if ($probValues[$i] > 0 && $probValues[$i] < 1); # no need to keep prob with only 100% & 0%
				}
				$pepVmod{$analysisID}{$queryNum{$analysisID}}{$vmod}[0]=$probString if $keepProb;
				next if ($seq2posString{$actualSequence} && $seq2posString{$actualSequence}{$vmod}); # If already computed, no need to search for position...
				foreach my $pos (sort {$probMod{$b} <=> $probMod{$a} || $a<=>$b} keys %probMod) {
					last if $nbMod == 0;
					$nbMod--;
					push @positions,$pos;
				}
				$seq2posString{$actualSequence}{$vmod}=join('.', sort{$a<=>$b} @positions);
			}
			else {
				next if ($seq2posString{$actualSequence} && $seq2posString{$actualSequence}{$vmod}); # If already computed, no need to search for position...
				my $posString;
				if ($vmod =~ /C-TERM/i) {
					$posString=($vmod =~ /PROTEIN/i) ? '+' : '*';
				}
				elsif ($vmod =~ /N-TERM/i) {
					$posString=($vmod =~ /PROTEIN/i) ? '-' : '=';
				}
				else { # No prob data at all (No MS/MS) => Match between runs => extract pos info from actualSequence string: _SEQ(mod1)UENC(mod2)E_
					my @positions;
					my ($numExtraChars,$numFound)=(0,0);
					pos($actualSequence)=0; # reset the regexp match position (set by any previous vmod!!!) to begining of string => critical because of the "last if $numFound==$nbMod;" command!
					while ($actualSequence=~/(.)(\([^)]+\)+)/g) { # /('residue')('(modif (res))')/
						my ($res,$modLength,$pos)=($1,length($2),$-[0]-$numExtraChars);
						if ($residue2varMod{$res} && $residue2varMod{$res} eq $vmod) { # in case different vMods on same sequence
							push @positions,$pos;
							$numFound++;
							last if $numFound==$nbMod;
						}
						$numExtraChars+=($modLength+1); # [BUGFIX] "+1" PP 25/10/19
					}
					$posString=join('.',@positions);
					$pepVmod{$analysisID}{$queryNum{$analysisID}}{$vmod}[0]='##PRB_MQ=-1';
				}			
				$seq2posString{$actualSequence}{$vmod}=$posString;
			}
		}
	}
	if ($counter>100000) {
		print LOG '.';
		$counter=0;
	}
}
close EVIDENCE;

####>Reconnectecting to DB<####
$dbh=&promsConfig::dbConnect('no_user');

##>Update FDR info
my $sthUpAna=$dbh->prepare("UPDATE ANALYSIS SET FALSE_DISCOVERY=? WHERE ID_ANALYSIS=?");
foreach my $analysisID (keys %queryNum) { # %maxProtMatch
	my %pepCount=(TARGET=>0,DECOY=>0);
	my %specCount=(TARGET=>0,DECOY=>0);
	foreach my $type (keys %pepCount) {
		next if (!$featureCount{$type} || !$featureCount{$type}{$analysisID});
		$pepCount{$type}=scalar keys %{$featureCount{$type}{$analysisID}};
		foreach my $actualSequence (keys %{$featureCount{$type}{$analysisID}}) {$specCount{$type}+=$featureCount{$type}{$analysisID}{$actualSequence};}
		$specCount{$type}-=$pepCount{$type};
	}
	$sthUpAna->execute("$pepCount{DECOY}:$specCount{DECOY}:$pepCount{TARGET}:$specCount{TARGET}",$analysisID); # numDecoy:(specDecoy-numDecoy):numTarget:(spectarget-numTarget)
}
$sthUpAna->finish;
undef %queryNum;
print LOG " Done (".strftime("%H:%M:%S %d/%m/%Y",localtime).").\n";

&addToFileStat("4/$numSteps Preprocessing protein list\n");
print LOG "4/$numSteps Preprocessing protein list...";
my ($protVisibility,$projectIdentMapID)=$dbh->selectrow_array("SELECT PROT_VISIBILITY,ID_IDENTIFIER FROM PROJECT WHERE ID_PROJECT=$projectID");
$protVisibility=0 unless $protVisibility;
$projectIdentMapID=0 unless $projectIdentMapID;
my ($giIdentID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GI'");
my (%projectProt,%bestProtVis); #,%proteinLength,%incompleteProtein,%noDescriptionProtein);
#my $sthP=$dbh->prepare("SELECT ID_PROTEIN,IDENTIFIER,PROT_LENGTH,PROT_DES FROM PROTEIN WHERE ID_PROJECT=$projectID");
my $sthP=($protVisibility)?
	$dbh->prepare("SELECT P.ID_PROTEIN,IDENTIFIER,MAX(VISIBILITY) FROM PROTEIN P,ANALYSIS_PROTEIN AP WHERE P.ID_PROTEIN=AP.ID_PROTEIN AND ID_PROJECT=$projectID GROUP BY P.ID_PROTEIN")
	: $dbh->prepare("SELECT ID_PROTEIN,IDENTIFIER FROM PROTEIN WHERE ID_PROJECT=$projectID"); # no need for best vis of $protVisibility==0
$sthP->execute || die $sthP->errstr;
$counter=0;
while (my ($protID,$identifier,$maxVis)=$sthP->fetchrow_array) {
	next unless $identifier; # identifier is null (previous ID protection not removed)
	$counter++;
	$projectProt{$identifier}=$protID;
	$bestProtVis{$identifier}=$maxVis // 0;
	if ($counter>100000) {
		print LOG '.';
		$counter=0;
	}
}
$sthP->finish;
print LOG '/.';

$dbh->do("DELETE FROM PROTEIN WHERE ID_PROJECT=$projectID AND IDENTIFIER IS NULL") || $dbh->errstr; # in case of error in previous process (rare!)
$dbh->commit;
#my $insProt=$dbh->prepare("INSERT INTO PROTEIN (ID_PROTEIN,ID_PROJECT) VALUES (?,$projectID)") || $dbh->errstr;
##>Counting new proteins
my %newProteins;
$counter=0;
foreach my $identifier (keys %matchList) {
	foreach my $analysisID (keys %{$matchList{$identifier}}) {
		##if ($projectProt{$identifier}) { # new proteins added to %projectProt
		##	$matchList{$identifier}{$analysisID}=$projectProt{$identifier};
		##}
		##else {
		##	$newProteins{$identifier}=1;
		##	$bestProtVis{$identifier}=0;
		##	#$matchList{$identifier}{$analysisID}=1;
		##}
		unless ($projectProt{$identifier}) {
			$newProteins{$identifier}=1;
			$bestProtVis{$identifier}=0;
		}
		$counter++;
		if ($counter>100000) {
			print LOG '.';
			$counter=0;
		}
	}
}
##my ($maxProteinID)=$dbh->selectrow_array("SELECT MAX(ID_PROTEIN) FROM PROTEIN");
##my $protectProtID = $maxProteinID + 2 + scalar keys %newProteins; # +1 should be enough
##$dbh->do("INSERT INTO PROTEIN (ID_PROTEIN,ID_PROJECT) VALUES ($protectProtID,$projectID)") || $dbh->errstr; # ID protection
##$dbh->commit;

##foreach my $identifier (keys %matchList) {
##	foreach my $analysisID (keys %{$matchList{$identifier}}) {
##		if ($projectProt{$identifier}) { # new proteins added to %projectProt
##			$matchList{$identifier}{$analysisID}=$projectProt{$identifier};
##		}
##		else {
##			$matchList{$identifier}{$analysisID}=++$maxProteinID;
##			$projectProt{$identifier}=$maxProteinID;
##			$newProteins{$identifier}=$maxProteinID;
##			$bestProtVis{$identifier}=0;
##		}
##		$counter++;
##		if ($counter>100000) {
##			print '.';
##			$counter=0;
##		}
##	}
##}

print LOG " Done (".strftime("%H:%M:%S %d/%m/%Y",localtime).").\n";

&addToFileStat("5/$numSteps Scanning databank(s)\n");
print LOG "5/$numSteps Scanning databank(s)...";
my (%contMatchList,%protDes,%protMW,%protOrg,%protLength);
my (@anaIDs)=sort{$a<=>$b} keys %bestScore;
my @contextData = (-2, \@anaIDs, $onCluster);
my $prevRefMatchList;
my $dbRank=0;
my $annotFromDB=10; # flag for &promsMod::getProInfo
foreach my $dbID (@databankIDs) {
	$dbRank++;
	my $prefix;
	my $refMatchList={};
	if ($dbRank==1) {
		if ($contaminantDB) { # separate proteins from contaminants
			my $contDbRank=scalar @databankIDs; # last in list
			foreach my $identifier (keys %matchList) {
				if ($identifier=~/^CON__/) {
					$contMatchList{$identifier}=$matchList{$identifier};
					$protDbRank{$identifier}=$contDbRank;
				}
				else {$refMatchList->{$identifier}=$matchList{$identifier};} # add to list to be scanned
			}
		}
		else {$refMatchList=\%matchList;}
	}
	else {
		foreach my $identifier (keys %{$prevRefMatchList}) {
			if ($protLength{$identifier}) {$protDbRank{$identifier}=$dbRank-1;} # matched during previous db scan
			elsif ($dbID==$contaminantDB) {$protDbRank{$identifier}=1;} # set to 1 because no more db scan
			else {$refMatchList->{$identifier}=$matchList{$identifier};} # add to list to be scanned
		}
		if ($dbID==$contaminantDB) {
			$refMatchList=\%contMatchList;
			$prefix='CON__';
		}
	}
#print '.';
	#my $prefix=($dbID==$contaminantDB)? 'CON__' : undef;
	#&promsMod::getProtInfo('silent',$dbh,$dbID,\@anaIDs,\%protDes,\%protMW,\%protOrg,\%protLength,\%{$matchList{$dbRank}},undef,$prefix);
	&promsMod::getProtInfo('quiet',$dbh,$dbID,\@contextData,\%protDes,\%protMW,\%protOrg,\%protLength,undef,$refMatchList,undef,$annotFromDB,$prefix) if scalar keys %{$refMatchList};
	$prevRefMatchList=$refMatchList;
#print '.';
}
unless ($contaminantDB) {
	foreach my $identifier (keys %{$prevRefMatchList}) {
		$protDbRank{$identifier}=1 unless $protDbRank{$identifier};
	}
}
undef %contMatchList;
print LOG " Done (".strftime("%H:%M:%S %d/%m/%Y",localtime).").\n";

$dbh->disconnect; # Not needed for a while...

###> Compute Match Info based on analysis.fasta files
&addToFileStat("6/$numSteps Processing peptide/protein match data\n");
print LOG "6/$numSteps Processing peptide/protein match data...";
my (%numMatches,%sequenceList);
$counter=0;
foreach my $analysisID (@anaIDs) {
	if (-e "$promsPath{peptide}/proj_$projectID/ana_$analysisID/analysis.fasta") {
		my $identifier; #,$des,$org,$mw,$length;
		my $sequence='';
		open (FAS,"$promsPath{peptide}/proj_$projectID/ana_$analysisID/analysis.fasta") || die "Unable to open $promsPath{peptide}/proj_$projectID/ana_$analysisID/analysis.fasta";
		while(<FAS>) {
			$counter++;
			if (/^>(\S+)/) {
				my $newIdentifier=$1;
				if ($sequence && $maxProtMatch{$analysisID}{$identifier}) { # test
					$sequenceList{$identifier}=$sequence;
					if (!$numMatches{$analysisID} || !$numMatches{$analysisID}{$identifier}) {	
						foreach my $actualSequence (keys %{$maxProtMatch{$analysisID}{$identifier}}) {
							my $pepSeq=$actualSeq2seq{$actualSequence};
							#my ($seqBefore,$resEnd)=($sequence=~/(.*)$pepSeq(\w?)/);
							if (!$numMatches{$analysisID}{$identifier}{$pepSeq}) { # !$numMatches{$analysisID} || !$numMatches{$analysisID}{$identifier} || 
								while ($sequence=~/(\w?)$pepSeq(\w?)/g) { # potential multiple matches
									my ($resBeg,$resEnd)=($1,$2);
									my $startPos=$-[0] + 1;
									my $endPos=$+[0];
									if ($resBeg) {$startPos++;}
									else {$resBeg='-';}
									if ($resEnd) {$endPos--;}
									else {$resEnd='-';}
						#@{$matchInfos{$analysisID}{$identifier}{$actualSequence}{$startPos}}=($resBeg,$resEnd,$endPos);
									@{$numMatches{$analysisID}{$identifier}{$pepSeq}{$startPos}}=($resBeg,$resEnd,$endPos);
								}
							}
						}
					}
				}
				$identifier=$newIdentifier;
				$sequence='' if defined($identifier);
            }
			else {
				#chomp;
				$_=~s/\W+//g; # chomp not enough? (Windows!)
				$sequence.=$_;
			}
			if ($counter>10000) {
				print LOG '.';
				$counter=0;
			}
		}
		close FAS;

		if ($sequence && $maxProtMatch{$analysisID}{$identifier}) {
			$sequenceList{$identifier}=$sequence;
			if (!$numMatches{$analysisID} || !$numMatches{$analysisID}{$identifier}) {			
				foreach my $actualSequence (keys %{$maxProtMatch{$analysisID}{$identifier}}) {
					my $pepSeq=$actualSeq2seq{$actualSequence};
					if (!$numMatches{$analysisID}{$identifier}{$pepSeq}) { # !$numMatches{$analysisID} || !$numMatches{$analysisID}{$identifier} || 
						while ($sequence=~/(\w?)$pepSeq(\w?)/g) { # potential multiple matches
							my ($resBeg,$resEnd)=($1,$2);
							my $startPos=$-[0] + 1;
							my $endPos=$+[0];
							if ($resBeg) {$startPos++;}
							else {$resBeg='-';}
							if ($resEnd) {$endPos--;}
							else {$resEnd='-';}
							#@{$matchInfos{$analysisID}{$identifier}{$actualSequence}{$startPos}}=($resBeg,$resEnd,$endPos);
							@{$numMatches{$analysisID}{$identifier}{$pepSeq}{$startPos}}=($resBeg,$resEnd,$endPos);
						}
					}
				}
			}
		}
		unlink "$promsPath{peptide}/proj_$projectID/ana_$analysisID/analysis.fasta"; # no longer needed
    }
}
print LOG " Done (".strftime("%H:%M:%S %d/%m/%Y",localtime).").\n";
#die "DEBUG";

###> Read peptide file to get context and be able to write MATCH_MULTI and MATCH_INFO pos,AA-1,AA+1
###> 2nd step : create a Design, Condition, Observation, etc.
&addToFileStat("7/$numSteps Reading peptide file\n");
print LOG "7/$numSteps Reading peptide file...";
my (%peptideColNum,%pepProteins);
open (PEPTIDE,"$tmpFilesDir/$peptideFile") || die "Unable to open $tmpFilesDir/$peptideFile";
$counter=0;
while (my $line=<PEPTIDE>) {
	$counter++;
	$line=~s/\s+$//; # chomp is not enough <- Windows
	my @parameters=split(/ *\t */,$line); # remove starting/trailing spaces
	if ($.==1) { # 1st line of the file
		my $ncol=0;
		foreach my $colName (@parameters) {
			$peptideColNum{$colName}=$ncol;
			$ncol++;
		}
		next;
	}
	next if $skipPeptide{$parameters[$evidenceColNum{uc('Sequence')}]};
	#'Proteins'
	#'Evidence IDs' --> Get the evidence that gave this peptide
	my $rawRazorProtein=$parameters[$peptideColNum{'Leading razor protein'}];
	$rawRazorProtein =~ s/REV__/DECOY_/;
	my $razorProtein = $raw2UsedIdentifiers{$rawRazorProtein} || $rawRazorProtein;
	foreach my $rawIdentifier ( split(/;/,$parameters[$peptideColNum{'Proteins'}]) ) { # Be careful : CON__ prefix stands for contaminants and REV__ prefix stands for decoy/reverse db
		last if $onlyBySite{$rawIdentifier};
		#my $matchGood=($rawIdentifier !~ /REV__/ && ($rawIdentifier !~ /CON__/ || !$excludeCON))? 1 : 0;
		my $matchGood=($rawIdentifier =~ /REV__/)? 0 : 1;
		$rawIdentifier =~ s/REV__/DECOY_/;
		my $identifier=$raw2UsedIdentifiers{$rawIdentifier} || $rawIdentifier;
		foreach my $evID (split(/;/,$parameters[$peptideColNum{'Evidence IDs'}])) {
			next unless $evidence2peptide{$evID}; # matches a decoy or excluded contaminant (or skipped parameter group data????????????????????)
			my ($analysisID,$queryNum,$actualSequence)=@{$evidence2peptide{$evID}};
			my $pepSeq=$actualSeq2seq{$actualSequence};
			$pepProteins{$analysisID}{$pepSeq}{$identifier}=1 if $matchGood; # for peptide/protein specificity

			#if (!$numMatches{$analysisID} || !$numMatches{$analysisID}{$identifier} || !$numMatches{$analysisID}{$identifier}{$pepSeq}) {
			#	if ($identifier eq $razorProtein && $peptideColNum{'Start position'}) { # info available in file from v~1.5 only for $razorProtein
			#		my $startPos=$parameters[$peptideColNum{'Start position'}];
			#		my $endPos=$startPos + length($pepSeq)-1;
			#		@{$numMatches{$analysisID}{$identifier}{$pepSeq}{$startPos}}=($parameters[$peptideColNum{'Amino acid before'}],$parameters[$peptideColNum{'Amino acid after'}],$endPos);
			#	}
			#	else {@{$numMatches{$analysisID}{$identifier}{$pepSeq}{0}}=('X','X',0);}
			#}
			if ($identifier eq $razorProtein && $peptideColNum{'Start position'}) { # info available in file from v~1.5 only for $razorProtein
				%{$numMatches{$analysisID}{$identifier}{$pepSeq}}=(); # reset in case initialized with fasta file
				my $startPos=$parameters[$peptideColNum{'Start position'}];
				my $endPos=($peptideColNum{'End position'})? $parameters[$peptideColNum{'End position'}] : $startPos + length($pepSeq)-1;
				@{$numMatches{$analysisID}{$identifier}{$pepSeq}{$startPos}}=($parameters[$peptideColNum{'Amino acid before'}],$parameters[$peptideColNum{'Amino acid after'}],$endPos);
			}
			elsif (!$numMatches{$analysisID} || !$numMatches{$analysisID}{$identifier} || !$numMatches{$analysisID}{$identifier}{$pepSeq}) {
				#@{$numMatches{$analysisID}{$identifier}{$pepSeq}{0}}=('X','X',0); # non-razor not mapped in fasta file
				@{$numMatches{$analysisID}{$identifier}{$pepSeq}{1}}=('X','X',length($pepSeq)); # non-razor not mapped in fasta file
			}
		}
	}
	if ($counter>10000) {
		print LOG '.';
		$counter=0;
	}
}
close PEPTIDE;
undef %evidence2peptide;
print LOG " Done (".strftime("%H:%M:%S %d/%m/%Y",localtime).").\n";

####>Reconnectecting to DB<####
$dbh=&promsConfig::dbConnect('no_user');

print LOG "8/$numSteps Recording peptide data:\n"; # SLOW!!!!!!!!!!!!!
###> Store Peptide Information
my $sthInsPep=$dbh->prepare("INSERT INTO PEPTIDE (ID_ANALYSIS,PEP_SEQ,PEP_LENGTH,QUERY_NUM,PEP_RANK,SEARCH_RANK,SCORE,MISS_CUT,MR_EXP,MR_CALC,MR_OBS,MR_DELTA,CHARGE,ELUTION_TIME,VALID_STATUS,DATA,SPEC_COUNT) VALUES (?,?,?,?,1,1,?,?,?,?,?,?,?,?,?,?,?)");
my $sthInsPM=$dbh->prepare("INSERT INTO PEPTIDE_MODIFICATION (ID_MODIFICATION,ID_PEPTIDE,MODIF_TYPE,POS_STRING,REF_POS_STRING) VALUES (?,?,'V',?,?)");
my (%pepIDs,%specificSequences);
my $anaCounter=0;
foreach my $analysisID (sort{$a<=>$b} keys %pepInfo) {
	$anaCounter++;
	&addToFileStat("8/$numSteps Recording peptide data: Analysis $anaCounter/$numAnalyses\n");
	print LOG " -$anaNames{$analysisID} ($anaCounter/$numAnalyses)";
	my $paramGr=$ana2ParamGroup{$anaNames{$analysisID}};
	my $counter=0;
	my $numQueries=scalar keys %{$pepInfo{$analysisID}};
	foreach my $qNum (sort{$a<=>$b} keys %{$pepInfo{$analysisID}}) {
#next if $queryNum < 0 ; # decoy matches are not kept in myProMS DB
		$counter++;
		print LOG '.' unless $counter % 5000; # 10000
		my ($actualSequence,@data)=@{$pepInfo{$analysisID}{$qNum}};
###next if $skipPeptide{$data[0]}; # sequence
		##>Ghost?
		#$data[2]=0 if ($isGhost{$analysisID} && $isGhost{$analysisID}{$qNum}); # queryNum
		if ($isGhost{$analysisID} && $isGhost{$analysisID}{$qNum}) {
			next unless $peptideXIC{$analysisID}{$qNum}; # ghost without quantif data
			$data[2]=0;
		}
		##>Peptide is specific ?
		my ($isSpecific)=($pepProteins{$analysisID}{$actualSeq2seq{$actualSequence}} && scalar keys %{$pepProteins{$analysisID}{$actualSeq2seq{$actualSequence}}}==1)? 1:0;
		$specificSequences{$analysisID}{$data[0]}=1 if $isSpecific;
		#print "$peptideID,$analysisID,@data<BR>\n";
		#next unless ($data[3] && $data[3] == $bestScore{$analysisID}{$actualSequence}); # Do not consider lower scoring peptides
#print LOG ">$analysisID: @data\n";
		$sthInsPep->execute($analysisID,@data,$featureCount{'TARGET'}{$analysisID}{$actualSequence}) || die $sthInsPep->errstr();
		my $peptideID = $dbh->last_insert_id(undef,undef,'PEPTIDE','ID_PEPTIDE');
		#@{$pepIDs{$analysisID}{$actualSequence}{$peptideID}}=($data[0],$data[11]); # keep peqSeq and validStatus for further coverage computing
		$pepIDs{$analysisID}{$actualSequence}{$peptideID}=$data[11]; # validStatus
		$peptideXIC{$analysisID}{$qNum}{'ID_PEPTIDE'}=$peptideID;
		foreach my $vmod (keys %{$pepVmod{$analysisID}{$qNum}}) { # non-label vmods
			$sthInsPM->execute($anaVmods{$analysisID}{$vmod},$peptideID,$seq2posString{$actualSequence}{$vmod},$pepVmod{$analysisID}{$qNum}{$vmod}[0]);
		}
		##>Labeling modifs
		#>SILAC: only lightest version with XIC data is recorded as true peptide, heavier versions(s) is/are recorded later as ghost(s)
		#OBSOLETE comment -> SILAC, only light version is recorded as true peptide (even if no XIC data), heavier versions(s) is/are recorded later as ghost(s)
		if ($labeling{$paramGr} eq 'SILAC') {
			foreach my $colLabel (@{$labelList{$paramGr}}) {
				if ($peptideXIC{$analysisID}{$qNum}{'XIC'}{$colLabel}) { # use the 1srt channel with a XIC for true peptide
					#last if $peptideXIC{$analysisID}{$qNum}{'TRUE_PEP_LABEL'};
					$peptideXIC{$analysisID}{$qNum}{'TRUE_PEP_LABEL'}=$colLabel; # record the label channel used for true peptide (to be skipped for ghosts)
					last unless $isotopeLabelInfo{$paramGr}{$colLabel}; # Light channel: No modif to add
					foreach my $modID (keys %{$isotopeLabelInfo{$paramGr}{$colLabel}}) {
						my $posString=&getModifPosString($data[0],$isotopeLabelInfo{$paramGr}{$colLabel}{$modID}); # $data[0]=pepSeq
						next unless $posString; # eg "Arg10;Lys8" only 1 of the 2 modifs will match!
						$sthInsPM->execute($modID,$peptideID,$posString,undef);
					}
					last; # TODO: Correct peptide masses & mass error in case not a Light channel
				}
			}
			$peptideXIC{$analysisID}{$qNum}{'TRUE_PEP_LABEL'}='L' unless $peptideXIC{$analysisID}{$qNum}{'TRUE_PEP_LABEL'}; # just to be safe			
		}
		#>Isobaric
		if ($isobarModifInfo{$paramGr} && $isobarModifInfo{$paramGr}{ANA} && $isobarModifInfo{$paramGr}{ANA}{$analysisID} && $isobarModifInfo{$paramGr}{ANA}{$analysisID}{$qNum}) { # reconstructed variable isobaric modif positions
			$sthInsPM->execute($isobarModifInfo{$paramGr}{ID},$peptideID,$isobarModifInfo{$paramGr}{ANA}{$analysisID}{$qNum},undef);
		}
	}
	#print '.100%' if $currPC < 100;
	print LOG " Done.\n";
	$dbh->commit;
	sleep 3;
}
$sthInsPep->finish;
$sthInsPM->finish;
undef %isGhost;
undef %featureCount;
undef %skipPeptide;
print LOG " Done (".strftime("%H:%M:%S %d/%m/%Y",localtime).").\n";


print LOG "9/$numSteps Recording protein data:\n";
my (%maxProtScore,%ppa);
my $sthInsProt=$dbh->prepare("INSERT INTO PROTEIN (ID_PROJECT,IDENTIFIER,ALIAS,PROT_DES,ORGANISM,MW,PROT_SEQ,PROT_LENGTH) VALUES ($projectID,?,?,?,?,?,?,?)");
my $insAP=$dbh->prepare("INSERT INTO ANALYSIS_PROTEIN (ID_ANALYSIS,ID_PROTEIN,SCORE,CONF_LEVEL,DB_RANK,NUM_PEP,NUM_MATCH,PEP_COVERAGE,MATCH_GROUP,PEP_SPECIFICITY,VISIBILITY) VALUES (?,?,?,?,?,?,?,?,?,?,?)");
my $sthTop=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=? AND VISIBILITY=2");
my $sthAtt=$dbh->prepare("INSERT INTO PEPTIDE_PROTEIN_ATTRIB (ID_PEPTIDE,ID_PROTEIN,ID_ANALYSIS,PEP_BEG,PEP_END,FLANKING_AA,IS_SPECIFIC) VALUES (?,?,?,?,?,?,?)");

##>MQ match groups
my (%MQmatchGroup,%MQvisibility);
if ($matchGroupType eq 'MaxQuant') {
	&createMaxQuantMatchGroups($protVisibility,\%MQmatchGroup,\%MQvisibility,\%matchList,\%bestProtVis,\%raw2UsedIdentifiers,\%onlyBySite); # same for all analyses
}
my @newAnaMapping;
my %protTopMG; # prot is in top of a match group in project
$anaCounter=0;
foreach my $analysisID (sort{$a<=>$b} keys %maxProtMatch) {
	$anaCounter++;
	&addToFileStat("9/$numSteps Recording protein data: Analysis $anaCounter/$numAnalyses\n");
	print LOG " -$anaNames{$analysisID} ($anaCounter/$numAnalyses)...";

	####>Fetching starting ID and protecting ID space for table PROTEIN
	##my ($proteinID)=$dbh->selectrow_array("SELECT MAX(ID_PROTEIN) FROM PROTEIN");
	##$proteinID=0 unless $proteinID;
	##my $maxProteinID=$proteinID + scalar (keys %{$maxProtMatch{$analysisID}}) + 1;
	##$dbh->do("INSERT INTO PROTEIN (ID_PROTEIN,ID_PROJECT) VALUES ($maxProteinID,$projectID)") || $dbh->errstr;
	##$dbh->commit;

	foreach my $identifier (keys %{$maxProtMatch{$analysisID}}) {
		#my $refBestScore=($identifier=~/DECOY_/)? \%bestScorebad : \%bestScore; # no longer decoy in %maxProtMatch (PP 15/11/16)
		foreach my $actualSequence (keys %{$maxProtMatch{$analysisID}{$identifier}}) {
			#$maxProtScore{$analysisID}{$identifier}+=$refBestScore->{$analysisID}{$actualSequence};
			$maxProtScore{$analysisID}{$identifier}+=$bestScore{$analysisID}{$actualSequence};
		}
	}

	###################################
	####>Updating ANALYSIS_PROTEIN<####
	###################################
	####>Computing PEP_SPECIFICITY
	my %pepSpecificity;#$pepProteins{$actualSequence}{$identifier}
	foreach my $pepSeq (keys %{$pepProteins{$analysisID}}) {
		my $specificity=sprintf "%.2f",100/(scalar keys %{$pepProteins{$analysisID}{$pepSeq}});
		$specificity*=1; # 100.00 -> 100
		foreach my $identifier (keys %{$pepProteins{$analysisID}{$pepSeq}}) {
			$pepSpecificity{$identifier}=$specificity if (!$pepSpecificity{$identifier} || $pepSpecificity{$identifier}<$specificity);
		}
	}
	print LOG '.';

	####>Finding match groups
	print LOG '/.';
	my ($refmatchGroup,$refvisibility)=($matchGroupType eq 'MaxQuant')? (\%MQmatchGroup,\%MQvisibility) : &createMatchGroups($analysisID,$protVisibility,\%bestProtVis,\%maxProtScore,\%maxProtMatch,\%proteinRank);

	####>Storing data in DB
	print LOG '/';
	my $newProtein=0;
	$counter=0;
	my (%boundaryStatus,%maxCovLength,%seqEndMatched);
	foreach my $identifier (keys %{$maxProtMatch{$analysisID}}) {
#next if $identifier =~ /DECOY_/;
		$counter++;
		my $des=&promsMod::resize($protDes{$identifier},250); # max length allowed in table
		my $organism=&promsMod::resize($protOrg{$identifier},100); # max length allowed in table
		my $score=($maxProtScore{$analysisID}{$identifier})?$maxProtScore{$analysisID}{$identifier}:0;
		if ($newProteins{$identifier}) { # protein is new to project=> update values !!!
			my $alias=($projectIdentMapID==$giIdentID && $identifier=~/(gi\|\d+)/)? $1 : $identifier; # "GI" restriction on alias
			$sthInsProt->execute($identifier,$alias,$des,$organism,$protMW{$identifier},$sequenceList{$identifier},$protLength{$identifier}) || die $sthInsProt->errstr;
			my $proteinID=$dbh->last_insert_id(undef,undef,'PROTEIN','ID_PROTEIN');
			$projectProt{$identifier}=$proteinID; # update list of project prot IDs
#foreach my $anaID (keys %{$matchList{$identifier}}) { # update for all matching analyses
#	$matchList{$identifier}{$anaID}=$proteinID;
#}
			$newProtein=1;
			delete $newProteins{$identifier};
		}
		###
		###
		###
		#my $numMatch=0;
		foreach my $actualSequence (keys %{$maxProtMatch{$analysisID}{$identifier}}){
			#my %usedBeg; # moved locally to sequence to allow '0' value when no beg info (PP 10/11/16)
			###my ($beg,$flAA_Nter,$flAA_Cter)=split(/,/,$matchInfos{$analysisID}{$identifier}{$actualSequence});
			####next if $usedBeg{$beg};
			###my $end;
			my $isSpecific=($specificSequences{$analysisID}{$actualSeq2seq{$actualSequence}})? 1 : undef;
			foreach my $beg (sort{$a<=>$b} keys %{$numMatches{$analysisID}{$identifier}{$actualSeq2seq{$actualSequence}}}) {
				my ($flAA_Nter,$flAA_Cter,$end)=@{$numMatches{$analysisID}{$identifier}{$actualSeq2seq{$actualSequence}}{$beg}};
				foreach my $peptideID (keys %{$pepIDs{$analysisID}{$actualSequence}}){
					#my ($pepSeq,$validStatus)=@{$pepIDs{$analysisID}{$actualSequence}{$peptideID}};
					my $validStatus=$pepIDs{$analysisID}{$actualSequence}{$peptideID};
					#$end=$beg+length($pepSeq)-1; # Length recorded before $pepIDs{$analysisID}{$actualSequence}{$peptideID}
					#my $isSpecific=($specificSequences{$analysisID}{$pepSeq})? 1 : undef;
					if ($validStatus==0) { # for Ghost peptides,
						$sthAtt->execute($peptideID,$projectProt{$identifier},$analysisID,-abs($beg),-abs($end),$flAA_Nter.$flAA_Cter,$isSpecific) || die $sthAtt->errstr;
					}
					else {
						$sthAtt->execute($peptideID,$projectProt{$identifier},$analysisID,$beg,$end,$flAA_Nter.$flAA_Cter,$isSpecific) || die $sthAtt->errstr;
					}
					@{$ppa{$peptideID}{$projectProt{$identifier}}}=($beg,$end,"$flAA_Nter$flAA_Cter");
				}
				$boundaryStatus{$analysisID}{$identifier}{$beg}++;
				$boundaryStatus{$analysisID}{$identifier}{$end}--;
				$maxCovLength{$analysisID}{$identifier}=$end if (!$maxCovLength{$analysisID}{$identifier} || $maxCovLength{$analysisID}{$identifier} < $end);
				$flAA_Nter='' unless $flAA_Nter;
				$flAA_Cter='' unless $flAA_Cter;
				$seqEndMatched{$analysisID}{$identifier}=1 if $flAA_Cter eq '-'; # peptide matches end of sequence
				#$numMatch{$analysisID}{$identifier}++;
			}
			#$numMatch+=scalar keys %{$matchInfos{$analysisID}{$identifier}{$actualSequence}};
		}
		###>Computing number of peptides & matches
		my ($numPep,$numMatch)=(0,0);
		my %usedPepSeq;
		foreach my $actualSequence (keys %{$maxProtMatch{$analysisID}{$identifier}}) {
			next unless $bestScore{$analysisID}{$actualSequence}; # =0 if only ghost
			$numPep++;
			my $pepSeq=$actualSeq2seq{$actualSequence};
			next if $usedPepSeq{$pepSeq}; # Actual number of matches on protein. Can be < numPep because of PTM...
			$numMatch+=scalar keys %{$numMatches{$analysisID}{$identifier}{$pepSeq}};
			$usedPepSeq{$pepSeq}=1;
		}

		###>Computing PEP_COVERAGE
		my $coverage=0;
		my $hasPeptide=0;
		my $boundaryNter=0;
		my $pepCoverage;
		if ($protLength{$identifier}) {
			foreach my $boundary (sort{$a<=>$b} keys %{$boundaryStatus{$analysisID}{$identifier}}) {
				if ($hasPeptide==0) { # start of peptide region (cannot become <0!)
					$boundaryNter=$boundary;
				}
				$hasPeptide+=$boundaryStatus{$analysisID}{$identifier}{$boundary};
				if ($hasPeptide==0) { # end of peptide region (should be true for last boundary too)
					$coverage+=($boundary-$boundaryNter)+1;
				}
			}
			my $usedProtLength=($maxCovLength{$analysisID}{$identifier} <= $protLength{$identifier})? $protLength{$identifier} : ($seqEndMatched{$analysisID}{$identifier})? $maxCovLength{$analysisID}{$identifier} : -1*$maxCovLength{$analysisID}{$identifier}; # -: flag for protLength problem
			$pepCoverage=sprintf "%.1f",(100*$coverage)/$usedProtLength;
			$pepCoverage*=1; # 25.0 -> 25
		}

		###
		###
		###
		#my $numPep=scalar keys %{$maxProtMatch{$analysisID}{$identifier}};
		#print "'$identifier' NUMPEP=$numPep $protDbRank{$identifier},$protMW{$identifier},$protLength{$identifier},$des,$organism,$score,$maxProtScore{$analysisID}{$identifier},$refmatchGroup->{$identifier},$projectProt{$identifier}<BR>\n";
		my $confLevel = ($score)? 2 : 0;
		#print "$analysisID,$projectProt{$identifier},$score,$protDbRank{$identifier},$confLevel,",scalar keys %{$matchInfos{$analysisID}{$identifier}},",",scalar keys %{$maxProtMatch{$analysisID}{$identifier}},",$pepCoverage,$refmatchGroup->{$identifier},$pepSpecificity{$identifier},$refvisibility->{$identifier},<BR>\n";
		$insAP->execute($analysisID,$projectProt{$identifier},$score,$confLevel,$protDbRank{$identifier},$numPep,$numMatch,$pepCoverage,$refmatchGroup->{$identifier},$pepSpecificity{$identifier},$refvisibility->{$identifier}) || die $insAP->errstr();
		$protTopMG{$identifier}=1 if (defined($refvisibility->{$identifier}) && $refvisibility->{$identifier}==2);
		if ($counter > 300) {
			print LOG '.';
			$counter=0;
		}

	}
	push @newAnaMapping,$analysisID if $newProtein;
	
	$dbh->commit;
	print LOG " Done.\n";
	sleep 3;
}
$sthInsProt->finish;
$insAP->finish;
$sthTop->finish;
#$dbh->do("DELETE FROM PROTEIN WHERE ID_PROTEIN=$protectProtID") || $dbh->errstr; # in case of error in previous process (rare!)
#$dbh->commit;
undef %maxProtMatch;
undef %actualSeq2seq;
undef %protDbRank;
undef %matchList;
undef %bestScore;
##undef %bestQuery;
undef %proteinRank;
undef %bestProtVis;
undef %newProteins;
undef %protDes;
undef %protMW;
undef %protOrg;
undef %protLength;
undef %numMatches;
undef %sequenceList;
undef %pepProteins;
undef %pepIDs;
undef %maxProtScore;
undef %MQmatchGroup;
undef %MQvisibility;
print LOG " Done (".strftime("%H:%M:%S %d/%m/%Y",localtime).").\n";

#$dbh->disconnect; exit; # DEBUG


my $sthUpQuantif=$dbh->prepare("UPDATE QUANTIFICATION SET QUANTIF_ANNOT=CONCAT(QUANTIF_ANNOT,?),STATUS=1,UPDATE_DATE=NOW() WHERE ID_QUANTIFICATION=?");

#######################################
####> Peptide quantification data <####
#######################################
#$importPepQuantif=0; # DEBUG
if ($importPepQuantif) {

	&addToFileStat("10/$numSteps Recording XIC data: ...\n");
	print LOG "10/$numSteps Recording peptide XIC data:\n";
	mkdir "$promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
	$anaCounter=0;
	#foreach my $pepQuantifID (values %ana2pepQuantification) {#}
	foreach my $paramGr (sort{$a<=>$b} keys %peptideQuantifID) {
		my $pepQuantifID=$peptideQuantifID{$paramGr};
		
		###>Opening peptide quantif data file(s)<###
		my %pepQuantHandle;
		my $quantifDir="$promsPath{quantification}/project_$projectID/quanti_$pepQuantifID";
		mkdir $quantifDir;
		if ($labeling{$paramGr} eq 'FREE') {
			open($pepQuantHandle{0},">$quantifDir/peptide_quantification.txt") || die $!;
			print {$pepQuantHandle{0}} "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n";
		}
		else {
			foreach my $colLabel (@{$labelList{$paramGr}}) {
				my $targetPos=$label2targetPos{$paramGr}{$colLabel};
				open($pepQuantHandle{$targetPos},">$quantifDir/peptide_quantification_$targetPos.txt") || die $!;
				print {$pepQuantHandle{$targetPos}} "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n";
			}
		}
		
		##>List of analyses with peptide XIC data in current parameter group
		my @anaList;
		foreach my $anaName (@{$group2AnaName{$paramGr}}) {
			my $analysisID=$anaInfo{$anaName}{'ID_ANALYSIS'};
			push @anaList,$analysisID if $peptideXIC{$analysisID};
		}
		###>SILAC<###
		if ($labeling{$paramGr} eq 'SILAC') {
			my $sthAddVP=$dbh->prepare("INSERT INTO PEPTIDE (ID_ANALYSIS,PEP_SEQ,PEP_LENGTH,MISS_CUT,CHARGE,ELUTION_TIME,DATA,VALID_STATUS,QUERY_NUM,PEP_RANK) VALUES (?,?,?,?,?,?,?,0,0,0)");
			my $sthAddVPM=$dbh->prepare("INSERT INTO PEPTIDE_MODIFICATION (ID_PEPTIDE,ID_MODIFICATION,MODIF_TYPE,POS_STRING,REF_POS_STRING) VALUES (?,?,'V',?,?)");
			my $sthUpPep=$dbh->prepare("UPDATE PEPTIDE SET DATA=CONCAT(DATA,?) WHERE ID_PEPTIDE=?");
			foreach my $analysisID (sort{$a<=>$b} @anaList) { # keys %peptideXIC
				$anaCounter++;
				$counter=0;
				&addToFileStat("10/$numSteps Recording $labeling{$paramGr} XIC data: Analysis $anaCounter/$numAnalyses\n");
				print LOG " -$anaNames{$analysisID} ($anaCounter/$numAnalyses) $labeling{$paramGr}...";
				my $qsetNum=0;
				foreach my $queryNum (sort{$a<=>$b} keys %{$peptideXIC{$analysisID}}) {
					next unless $peptideXIC{$analysisID}{$queryNum}{'XIC'}; # no quantif data all
					$qsetNum++;
					foreach my $colLabel (@{$labelList{$paramGr}}) {
						$counter++;
						if ($counter > 500) {
							print LOG '.';
							$counter=0;
						}
						next unless $peptideXIC{$analysisID}{$queryNum}{'XIC'}{$colLabel};
						my $targetPos=$label2targetPos{$paramGr}{$colLabel};
						if ($isotopeLabelInfo{$paramGr} && $isotopeLabelInfo{$paramGr}{$colLabel} && (!$peptideXIC{$analysisID}{$queryNum}{'TRUE_PEP_LABEL'} || $peptideXIC{$analysisID}{$queryNum}{'TRUE_PEP_LABEL'} ne $colLabel) ) {
							#my $pepID=$peptideXIC{$analysisID}{$queryNum}{'ID_PEPTIDE'};
							my ($actualSequence,@data)=@{$pepInfo{$analysisID}{$queryNum}}; # @data=(pepSeq,Length,$anaID,$score,$missCut,$massExp,$Mass,$mz,$massError,$charge,$et,$vamidStatus,$data,$specCount)
							my $pepSeq=$data[0];
							###> Create a Ghost Peptide
							my $isSpecific=($specificSequences{$analysisID}{$pepSeq})? 1 : undef;
							$sthAddVP->execute($analysisID,$pepSeq,$data[1],$data[4],$data[9],$data[10],"QSET=$qsetNum");
							my $vpPepID=$dbh->last_insert_id(undef,undef,'PEPTIDE','ID_PEPTIDE');
							###> Add modifications
							my $posString='';
							###> Label ones
							foreach my $modID (keys %{$isotopeLabelInfo{$paramGr}{$colLabel}}) {
								#if ($isotopeLabelInfo{$paramGr}{$colLabel}{$modID}{'*'}) {
								#	$posString='*';
								#}
								#elsif ($isotopeLabelInfo{$paramGr}{$colLabel}{$modID}{'='}) {
								#	$posString='=';
								#}
								#else {
								#	my @aas=split(//,$pepSeq);
								#	my @pos;
								#	for (my $i = 0 ; $i <= $#aas ; $i++) {
								#		push @pos,$i+1 if $isotopeLabelInfo{$paramGr}{$colLabel}{$modID}{$aas[$i]}; # or $labelModifSpecificity{$paramGr}{$modID}{$aas[$i]}
								#	}
								#	next unless scalar @pos; # eg "Arg10;Lys8" only 1 of the 2 modifs will match!
								#	$posString=join('.',@pos);
								#}
								my $posString=&getModifPosString($pepSeq,$isotopeLabelInfo{$paramGr}{$colLabel}{$modID});
								next unless $posString; # eg "Arg10;Lys8" only 1 of the 2 modifs will match!
								$sthAddVPM->execute($vpPepID,$modID,$posString,undef);
							}
							###> Regular ones (vmod like oxidation, acetylation,...)
							foreach my $vmod (keys %{$pepVmod{$analysisID}{$queryNum}}) {
								$sthAddVPM->execute($vpPepID,$anaVmods{$analysisID}{$vmod},$seq2posString{$actualSequence}{$vmod},$pepVmod{$analysisID}{$queryNum}{$vmod}[0]); # use same site probability as validated peptide (good, bad?)
							}
							###> Create a link to PEPTIDE_PROTEIN_ATTRIB with this ghost peptide
							foreach my $protID (keys %{$ppa{$peptideXIC{$analysisID}{$queryNum}{'ID_PEPTIDE'}}}) {
								my ($beg,$end,$matchInfo)=@{$ppa{$peptideXIC{$analysisID}{$queryNum}{'ID_PEPTIDE'}}{$protID}};
								$sthAtt->execute($vpPepID,$protID,$analysisID,-abs($beg),-abs($end),$matchInfo,$isSpecific);
							}
							###> Add quantification for this specific Ghost peptide
							print {$pepQuantHandle{$targetPos}} "$areaParamID{$paramGr}\t$vpPepID\t$peptideXIC{$analysisID}{$queryNum}{XIC}{$colLabel}\n";
							#print "$quantiSILAC{$analysisID},$areaParamID,$vpPepID,$peptideXIC{$analysisID}{$queryNum}{\"XIC $colLabel\"},$labelInfos{$analysisID}{$colLabel}{'TARGET_POS'}<BR>\n";
						}
						else { # No modification ID (Light label) or Label used by true peptide
							#print "$quantiSILAC{$analysisID},$areaParamID,$peptideXIC{$analysisID}{$queryNum}{'ID_PEPTIDE'},$peptideXIC{$analysisID}{$queryNum}{\"XIC $colLabel\"},$labelInfos{$analysisID}{$colLabel}{'TARGET_POS'}<BR>";
							print {$pepQuantHandle{$targetPos}} "$areaParamID{$paramGr}\t$peptideXIC{$analysisID}{$queryNum}{ID_PEPTIDE}\t$peptideXIC{$analysisID}{$queryNum}{XIC}{$colLabel}\n";
							#my ($data)=$dbh->selectrow_array("SELECT DATA FROM PEPTIDE WHERE ID_PEPTIDE=$peptideXIC{$analysisID}{$queryNum}{'ID_PEPTIDE'}");
							#$data.="##QSET=$qsetNum";
							$sthUpPep->execute("##QSET=$qsetNum",$peptideXIC{$analysisID}{$queryNum}{'ID_PEPTIDE'});
						}
					}
				}
				$dbh->commit;
				print LOG " Done.\n";
			}
			$sthAddVP->finish;
			$sthAddVPM->finish;
			$sthUpPep->finish;
		}
		###>iTRAQ or TMT <###
		elsif ($isobarModifInfo{$paramGr}) { # $labeling=~/ITRAQ|TMT/
			foreach my $analysisID (sort{$a<=>$b} @anaList) { # keys %peptideXIC
				$anaCounter++;
				$counter=0;
				&addToFileStat("10/$numSteps Recording $labeling{$paramGr} XIC data: Analysis $anaCounter/$numAnalyses\n");
				print LOG " -$anaNames{$analysisID} ($anaCounter/$numAnalyses) $labeling{$paramGr}...";
				foreach my $queryNum (sort{$a<=>$b} keys %{$peptideXIC{$analysisID}}) {
					next unless $peptideXIC{$analysisID}{$queryNum}{'XIC'}; # no quantif data all
					foreach my $repIdx (@{$labelList{$paramGr}}) { #(sort{$a<=>$b} keys %{$peptideXIC{$analysisID}{$queryNum}{XIC}}) <<<<<<<<<<<<<<<
						next unless $peptideXIC{$analysisID}{$queryNum}{'XIC'}{$repIdx}; # skip 0/undef values
						$counter++;
						if ($counter > 500) {
							print LOG '.';
							$counter=0;
						}
						my $targetPos=$repIdx+1;
						print {$pepQuantHandle{$targetPos}} "$intensityParamID{$paramGr}\t$peptideXIC{$analysisID}{$queryNum}{ID_PEPTIDE}\t$peptideXIC{$analysisID}{$queryNum}{XIC}{$repIdx}\n";
					}
				}
				print LOG " Done.\n";
			}
		}
		###>Label-free<###
		elsif ($labeling{$paramGr} eq 'FREE') {
			foreach my $analysisID (sort{$a<=>$b} @anaList) { # keys %peptideXIC
				$anaCounter++;
				$counter=0;
				&addToFileStat("10/$numSteps Recording Label-free XIC data: Analysis $anaCounter/$numAnalyses\n");
				print LOG " -$anaNames{$analysisID} ($anaCounter/$numAnalyses) $labeling{$paramGr}...";
				foreach my $queryNum (sort{$a<=>$b} keys %{$peptideXIC{$analysisID}}) {
					$counter++;
					if ($counter > 500) {
						print LOG '.';
						$counter=0;
					}
					print {$pepQuantHandle{0}} "$areaParamID{$paramGr}\t$peptideXIC{$analysisID}{$queryNum}{ID_PEPTIDE}\t$peptideXIC{$analysisID}{$queryNum}{XIC}{NONE}\n" if $peptideXIC{$analysisID}{$queryNum}{'XIC'};
				}
				print LOG " Done.\n";
			}
		}

		foreach my $targetPos (keys %pepQuantHandle) {close $pepQuantHandle{$targetPos};}

		$sthUpQuantif->execute('',$pepQuantifID);
		
		$dbh->commit;

	}
	print LOG " Done (".strftime("%H:%M:%S %d/%m/%Y",localtime).").\n";

} # END of $importPepQuantif

$sthAtt->finish;
undef %pepInfo;
undef %pepVmod;
undef %seq2posString;
undef %peptideXIC;
undef %specificSequences;
undef %ppa;


#######################################
####> Protein quantification data <####
#######################################
if ($importProtQuantif) {
	
	&addToFileStat("11/$numSteps Importing protein quantification data\n");
	print LOG "11/$numSteps Importing protein quantification data from proteinGroups file...";
	my $pepTotStg=($pepUsed eq 'razor')? 'Peptide counts (razor+unique)' : ($pepUsed eq 'unique')? 'Peptide counts (unique)' : 'Peptide counts (all)';
	
	#my $sthInsProtQ=$dbh->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_PROTEIN,ID_QUANTIFICATION,ID_QUANTIF_PARAMETER,QUANTIF_VALUE,TARGET_POS) VALUES (?,?,?,?,?)");
	my (%dbhLite,%sthInsProtQ);
	foreach my $quantifID (@protQuantifList) {
		my $quantifDir="$promsPath{quantification}/project_$projectID/quanti_$quantifID";
		$dbhLite{$quantifID}=&promsQuantif::dbCreateProteinQuantification($quantifID,$quantifDir);
		$sthInsProtQ{$quantifID}=$dbhLite{$quantifID}->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_PROTEIN,ID_QUANTIF_PARAMETER,QUANTIF_VALUE,TARGET_POS) VALUES (?,?,?,?)");
	}
	
	open (PROTEIN,"$tmpFilesDir/$proteinGroupsFile") || die "Unable to open $tmpFilesDir/$proteinGroupsFile";

	my (%proteinColNum,%usedQuantifParams);
	$counter=0;
	while (my $line=<PROTEIN>) {
		$line=~s/\s+$//; # chomp is not enough <- Windows
		my @parameters=split(/ *\t */,$line); # remove starting/trailing spaces
		if ($.==1) { # 1st line of the file
			my $ncol=0;
			foreach my $colName (@parameters) {
				$proteinColNum{$colName}=$ncol;
				$ncol++;
			}
			next;
		}
		###> In proteinGroups.txt, you find in the column the distribution of peptides found for each entry
		###> The idea is to keep only quantification values made for the same set of peptides.
		###> ex : Protein IDs column => O43790;CON__O43790;CON__Q14533;Q14533
		###>      Peptide counts (all) column => 3;3;3;2
		###>      Then, only protein O43790, CON__O43790 and CON__Q14533 will be kept in MQ quantitation
		###> PP 24/11/16: Restrict rule above to proteins at top of a match group in at least 1 analysis
		my @pepTotal=split(/;/,$parameters[$proteinColNum{$pepTotStg}]);
		my @identifiers; #=split(/;/,$parameters[$proteinColNum{'Majority protein IDs'}]); # 'Protein IDs'
		foreach my $rawIdentifier (split(/;/,$parameters[$proteinColNum{'Majority protein IDs'}])) {
			last if $onlyBySite{$rawIdentifier};
			push @identifiers,$raw2UsedIdentifiers{$rawIdentifier} || $rawIdentifier;
		}
		next if (!$identifiers[0] || !$projectProt{$identifiers[0]}); # Skip "Only identified by site" & "dirty" lines... (PP 10/11/16)
		my %selProtFromMQgroup;
		$selProtFromMQgroup{$projectProt{$identifiers[0]}}=$pepTotal[0] if $protTopMG{$identifiers[0]}; # must be top of a MG in >= 1 MQ ana (bestVis=2) Should be always true for this identifier # if $projectProt{$identifiers[0]};
		for (my $i=1; $i<=$#identifiers; $i++){
			last if $pepTotal[$i] < $pepTotal[0]; # in case of very large groups, peptide counts are not recorded for all member
			$selProtFromMQgroup{$projectProt{$identifiers[$i]}}=$pepTotal[$i] if $protTopMG{$identifiers[$i]}; # must be top of a MG in >= 1 MQ ana (bestVis=2)
		}
		foreach my $paramText (keys %recordedParams) {  # %{$designInfo{'Quantification'}}
			next unless $proteinColNum{$paramText};
			my $ncol = $proteinColNum{$paramText};
			next unless $ncol>=0; # ?
			$parameters[$ncol]=~s/;.+//; # Global 'Peptide counts xxx' column ('X;Y;Z;..') is used when no design. Only 1st value is kept (PP 08/11/16)
			next if (!$parameters[$ncol] || $parameters[$ncol] eq 'NaN' || $parameters[$ncol] == 0); # avoid to store 0 values in myProMS db
			foreach my $refQuantif (@{$recordedParams{$paramText}}) {
				my ($qID,$qCode,$targetPos)=@{$refQuantif}; # @{$designInfo{'Quantification'}{$colName}};
				next unless $quantifParamIDs{$qCode};
				my $paramOK=0;
				if ($targetPos==0) { # required quantif value of at least 1 channel must be OK
					foreach my $tgPos (keys %{$requiredParams{$qID}}) {
						if ($parameters[$proteinColNum{$requiredParams{$qID}{$tgPos}}] && $parameters[$proteinColNum{$requiredParams{$qID}{$tgPos}}] ne 'NaN') {
							$paramOK=1;
							last;
						}
					}
					next unless $paramOK;
				}
				elsif ($parameters[$proteinColNum{$requiredParams{$qID}{$targetPos}}] && $parameters[$proteinColNum{$requiredParams{$qID}{$targetPos}}] ne 'NaN') { # do not record any values if required one fails
					$paramOK=1;
				}
				next unless $paramOK;
				my $usedTargetPos=$targetPos || undef; # 0 -> undef
				foreach my $protID (keys %selProtFromMQgroup) {
					#$sthInsProtQ->execute($protID,$qID,$quantifParamIDs{$qCode},$parameters[$ncol],$usedTargetPos);
					$sthInsProtQ{$qID}->execute($protID,$quantifParamIDs{$qCode},$parameters[$ncol],$usedTargetPos);
				}
				$usedQuantifParams{$qID}{$qCode}=1 if $protIntQuantifs{$qID}; # to record in quantifAnnot
			}
		}

		$counter++;
		if ($counter > 250) {
			print LOG '.';
			$counter=0;
		}
	}
	close PROTEIN;
	#$sthInsProtQ->finish;
	
	foreach my $quantifID (@protQuantifList) {
		$sthInsProtQ{$quantifID}->finish;
		$dbhLite{$quantifID}->commit;
		$dbhLite{$quantifID}->disconnect;

		#>Recording measures used<#
		my $measStrg='';
		if ($usedQuantifParams{$quantifID}) {
			foreach my $paramCode ('MQ_INT','MQ_IBAQ','MQ_LFQ','MQ_SC') {
				next unless $usedQuantifParams{$quantifID}{$paramCode};
				$measStrg.=';' if $measStrg;
				$measStrg.=$paramCode;
			}
			$measStrg='::ABUND_MEASURES='.$measStrg if $measStrg;
		}
		$sthUpQuantif->execute($measStrg,$quantifID); # measStrg='' if RATIO
	}
	$dbh->commit;
	print LOG " Done (".strftime("%H:%M:%S %d/%m/%Y",localtime).").\n";

}
$sthUpQuantif->finish;
undef %raw2UsedIdentifiers;
undef %projectProt;
undef %protTopMG;
$dbh->disconnect;

###>Move all files in peptide quantitation folder & link to protein quantification folder(s) !
&addToFileStat("12/$numSteps Moving data files\n");
print LOG "12/$numSteps Moving data files...";
mkdir "$promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
my $firstQuantifID=0;
foreach my $quantifID (sort{$a<=>$b} values %peptideQuantifID,@protQuantifList) {
	mkdir "$promsPath{quantification}/project_$projectID/quanti_$quantifID"; # unless -e "$promsPath{quantification}/project_$projectID/quanti_$quantifID";
	foreach my $file (@requiredFiles) {
		if ($firstQuantifID) {
			symlink("../quanti_$firstQuantifID/$file","$promsPath{quantification}/project_$projectID/quanti_$quantifID/$file"); # symbolic link with relative path to source
		}
		else {
			move("$tmpFilesDir/$file","$promsPath{quantification}/project_$projectID/quanti_$quantifID/$file");
#copy("$tmpFilesDir/$file","$promsPath{quantification}/project_$projectID/quanti_$quantifID/$file"); # DEBUG
		}
	}
	$firstQuantifID=$quantifID unless $firstQuantifID;
}
print LOG " Done (".strftime("%H:%M:%S %d/%m/%Y",localtime).").\n";
close LOG;
move("$tmpFilesDir/$detailFile","$promsPath{quantification}/project_$projectID/quanti_$firstQuantifID/$detailFile");
#sleep 2;
#rmtree $tmpFilesDir;

#print '**DONE**'; exit; #DEBUG !!!

################################
####>New identifiers to map<####
################################
if (scalar @newAnaMapping) { # true if new valid proteins in project
	open(MAPPING,">$tmpFilesDir/ana2map.txt") || die "Error while opening ana2map.txt: $!";
	my $anaStrg=join(',',@newAnaMapping);
	print MAPPING $anaStrg;
	close MAPPING;
}
&addToFileStat("Ended ".strftime("%H:%M:%S %d/%m/%Y",localtime)."\n");
sleep 30;

#########################
####<Update FILESTAT>####
#########################
sub addToFileStat {
	open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat: $!";
	print FILESTAT $_[0]; # SLOW!!!!!!!!!!!!!
	close FILESTAT;
}

################################
####<Get Modification position on sequence>####
################################
sub getModifPosString {
	my ($pepSeq,$refSpecificity)=@_;
	my (@posEnds,@posRes);
	foreach my $res (sort keys %{$refSpecificity}) {
		if ($res =~ /[\=\-\*\+]/) {push @posEnds,$res;}
		else {
			while ($pepSeq =~ /$res/g) {push @posRes,$-[0]+1;}
		}
	}
	return join('.',(@posEnds,sort{$a<=>$b} @posRes));
}

##########################################
####<Filter "Only identified by site">####
##########################################
sub filterOnlyIdentifiedBySite {
	my ($refOnlyBySite)=@_;
	open (PROTEIN,"$tmpFilesDir/$proteinGroupsFile") || die "Unable to open $tmpFilesDir/$proteinGroupsFile";
	my ($allProtColumnIdx,$onlyBySiteColumnIdx)=(1,0);
	while (my $line=<PROTEIN>) {
		my @parameters=split(/ *\t */,$line); # remove starting/trailing spaces
		if ($.==1) { # 1st line of the file
			my $ncol=0;
			foreach my $colName (@parameters) {
				if ($colName eq 'Proteins IDs') {$allProtColumnIdx=$ncol;} # all proteins in match group
				elsif ($colName eq 'Only identified by site') {$onlyBySiteColumnIdx=$ncol;}
				$ncol++;
			}
			next;
		}
		last if !$onlyBySiteColumnIdx;
		next unless $parameters[$onlyBySiteColumnIdx]; # Fully identified
		foreach my $rawIdentifier (split(/;/,$parameters[$allProtColumnIdx])) {
			$refOnlyBySite->{$rawIdentifier}=1;
		}
	}
	close PROTEIN;
}

################################
####<Generates Match Groups>####
################################
sub createMatchGroups {
	my ($analysisID,$protVisibility,$refBestProtVis,$refMaxScore,$refMaxProtMatch,$refProteinRank)=@_;
	my (%matchGroup,%visibility);
	my $numGroup=0;
	my @sortedIdentifiers=(sort{scalar (keys %{$refMaxProtMatch->{$analysisID}{$b}})<=>scalar (keys %{$refMaxProtMatch->{$analysisID}{$a}}) || $refMaxScore->{$analysisID}{$b}<=>$refMaxScore->{$analysisID}{$a} || $refProteinRank->{$a}<=>$refProteinRank->{$b} || $a cmp $b} keys %{$refMaxProtMatch->{$analysisID}});
	my $counter=0;
	foreach my $i (0..$#sortedIdentifiers) {
		next if $matchGroup{$sortedIdentifiers[$i]}; # already assigned to a match group
		$matchGroup{$sortedIdentifiers[$i]}=++$numGroup;
        $visibility{$sortedIdentifiers[$i]}=2; # Reference protein
		$refBestProtVis->{$sortedIdentifiers[$i]}=2; # update bestVis
	#next; # SKIP grouping!!!
		foreach my $j ($i+1..$#sortedIdentifiers) {
           next if $matchGroup{$sortedIdentifiers[$j]}; # already assigned to a match group
			##<Comparing peptide contents of identifier#1 and identifier#2>## All peptides must match!
			my $matchOK=1;
			foreach my $seq (keys %{$refMaxProtMatch->{$analysisID}{$sortedIdentifiers[$j]}}) {
				$counter++;
				if ($counter >= 100000) {
					print '.';
					$counter=0;
				}
 				if (!$refMaxProtMatch->{$analysisID}{$sortedIdentifiers[$i]}{$seq}) {
					delete $refMaxProtMatch->{$analysisID}{$sortedIdentifiers[$i]}{$seq}; # to be safe
					$matchOK=0;
					last;
				}
			}
			if ($matchOK) {
				$matchGroup{$sortedIdentifiers[$j]}=$matchGroup{$sortedIdentifiers[$i]};
				$visibility{$sortedIdentifiers[$j]}=($protVisibility && defined($refBestProtVis->{$sortedIdentifiers[$j]}) && ($refBestProtVis->{$sortedIdentifiers[$j]}==2 || ($protVisibility==2 && $refBestProtVis->{$sortedIdentifiers[$j]})))? 1 : 0;
				$refBestProtVis->{$sortedIdentifiers[$j]}=$visibility{$sortedIdentifiers[$j]} if (!defined($refBestProtVis->{$sortedIdentifiers[$j]}) || $visibility{$sortedIdentifiers[$j]} > $refBestProtVis->{$sortedIdentifiers[$j]}); # update bestVis
			}
		}
	}
	return(\%matchGroup,\%visibility);
}
sub createMaxQuantMatchGroups {
	my ($protVisibility,$refMQmatchGroup,$refMQvisibility,$refmatchList,$refBestProtVis,$refRaw2UsedIdentifiers,$refOnlyBySite)=@_;
	print LOG "&nbsp;-Applying MaxQuant match group rules...";
	open (PROTEIN,"$tmpFilesDir/$proteinGroupsFile") || die "Unable to open $tmpFilesDir/$proteinGroupsFile";
	my $counter=0;
	my $mgNum=0;
	my ($allProtColumnIdx,$topProtColumnIdx)=(0,1);
	while (my $line=<PROTEIN>) {
		my @parameters=split(/ *\t */,$line); # remove starting/trailing spaces
		if ($.==1) { # 1st line of the file
			my $ncol=0;
			foreach my $colName (@parameters) {
				if ($colName eq 'Proteins IDs') {$allProtColumnIdx=$ncol;} # all proteins in match group
				elsif ($colName eq 'Majority protein IDs') {$topProtColumnIdx=$ncol;} # most abundant proteins in match group (contain half of MG peptides)
				$ncol++;
			}
			next;
		}
		#my @topIdentifiers=split(/;/,$parameters[$topProtColumnIdx]);
		my $topIdentifier;
		foreach my $rawIdentifier (split(/;/,$parameters[$topProtColumnIdx])) {
			last if $refOnlyBySite->{$rawIdentifier}; # skip "Only identified by site"
			next if $rawIdentifier =~/REV__/;
			#$rawIdentifier =~ s/REV__/DECOY_/; # shouldn't be useful since DECOY are not in %{$refmatchList}
			my $identifier=$refRaw2UsedIdentifiers->{$rawIdentifier} || $rawIdentifier;
			if ($refmatchList->{$identifier}) { # good prot +/- CON__ depending on $excludeCON
				$topIdentifier=$identifier;
				last;
			}
		}
		next unless $topIdentifier; # Skip entire CON__-only MG if $excludeCON
		#next unless $refmatchList->{$topIdentifiers[0]}; # Skip "dirty" lines... (PP 10/11/16)
		$mgNum++;
		foreach my $rawIdentifier (split(/;/,$parameters[$allProtColumnIdx])) {
			my $identifier=$refRaw2UsedIdentifiers->{$rawIdentifier} || $rawIdentifier;
			next unless $refmatchList->{$identifier};
			$refMQmatchGroup->{$identifier}=$mgNum;
			$refMQvisibility->{$identifier}=($identifier eq $topIdentifier)? 2 : ($protVisibility && defined($refBestProtVis->{$identifier}) && ($refBestProtVis->{$identifier}==2 || ($protVisibility==2 && $refBestProtVis->{$identifier})))? 1 : 0;
		}
		$counter++;
		if ($counter >= 250) {
			print LOG '.';
			$counter=0;
		}
	}
	close PROTEIN;
	print LOG " Done.\n";
}

###########################################################
####<Table of Modifications for report ion in MaxQuant>####
###########################################################
sub getModIDforReporterIon { # Modification with corresponding Unimod ID as checked on Unimod website 09 of November 2015
	my %mqRepIon=();
	$mqRepIon{'18O'}{'UNIMOD_ACC'}=193;
	$mqRepIon{'Arg10'}{'UNIMOD_ACC'}=267;
	$mqRepIon{'Arg6'}{'UNIMOD_ACC'}=188;
	$mqRepIon{'DimethLys0'}{'UNIMOD_ACC'}=36;
	$mqRepIon{'DimethLys2'}{'UNIMOD_ACC'}=199;
	$mqRepIon{'DimethLys4'}{'UNIMOD_ACC'}=510;
	$mqRepIon{'DimethLys6'}{'UNIMOD_ACC'}=1291;
	$mqRepIon{'DimethLys8'}{'UNIMOD_ACC'}=330;
	$mqRepIon{'DimethNter0'}{'UNIMOD_ACC'}=36;
	$mqRepIon{'DimethNter2'}{'UNIMOD_ACC'}=199;
	$mqRepIon{'DimethNter4'}{'UNIMOD_ACC'}=510;
	$mqRepIon{'DimethNter6'}{'UNIMOD_ACC'}=1291;
	$mqRepIon{'DimethNter8'}{'UNIMOD_ACC'}=330;
	$mqRepIon{'ICAT-0'}{'UNIMOD_ACC'}=105;
	$mqRepIon{'ICAT-9'}{'UNIMOD_ACC'}=106;
	$mqRepIon{'ICPL-Lys0'}{'UNIMOD_ACC'}=365;
	$mqRepIon{'ICPL-Lys10'}{'UNIMOD_ACC'}=687;
	$mqRepIon{'ICPL-Lys4'}{'UNIMOD_ACC'}=364;
	$mqRepIon{'ICPL-Lys6'}{'UNIMOD_ACC'}=866;
	$mqRepIon{'ICPL-Nter0'}{'UNIMOD_ACC'}=365;
	$mqRepIon{'ICPL-Nter10'}{'UNIMOD_ACC'}=687;
	$mqRepIon{'ICPL-Nter4'}{'UNIMOD_ACC'}=364;
	$mqRepIon{'ICPL-Nter6'}{'UNIMOD_ACC'}=866;
	$mqRepIon{'Ile7'}{'UNIMOD_ACC'}=695;
	$mqRepIon{'Leu7'}{'UNIMOD_ACC'}=695;
	$mqRepIon{'Lys4'}{'UNIMOD_ACC'}=481;
	$mqRepIon{'Lys6'}{'UNIMOD_ACC'}=188;
	$mqRepIon{'Lys8'}{'UNIMOD_ACC'}=259;
	$mqRepIon{'mTRAQ-Lys0'}{'UNIMOD_ACC'}=888;
	$mqRepIon{'mTRAQ-Lys4'}{'UNIMOD_ACC'}=889;
	$mqRepIon{'mTRAQ-Lys8'}{'UNIMOD_ACC'}=1302;
	$mqRepIon{'mTRAQ-Nter0'}{'UNIMOD_ACC'}=888;
	$mqRepIon{'mTRAQ-Nter4'}{'UNIMOD_ACC'}=889;
	$mqRepIon{'mTRAQ-Nter8'}{'UNIMOD_ACC'}=1302;

	$mqRepIon{'18O'}{'INTERIM_NAME'}='double_O18';
	$mqRepIon{'Arg10'}{'INTERIM_NAME'}='13C6-15N4';
	$mqRepIon{'Arg6'}{'INTERIM_NAME'}='13C6';
	$mqRepIon{'DimethLys0'}{'INTERIM_NAME'}='di-Methylation ';
	$mqRepIon{'DimethLys2'}{'INTERIM_NAME'}='CHD2';
	$mqRepIon{'DimethLys4'}{'INTERIM_NAME'}='C13HD2';
	$mqRepIon{'DimethLys6'}{'INTERIM_NAME'}='Dimethyl:2H(6)';
	$mqRepIon{'DimethLys8'}{'INTERIM_NAME'}='Dimethyl:2H(6)13C(2)';
	$mqRepIon{'DimethNter0'}{'INTERIM_NAME'}='di-Methylation';
	$mqRepIon{'DimethNter2'}{'INTERIM_NAME'}='CHD2';
	$mqRepIon{'DimethNter4'}{'INTERIM_NAME'}='C13HD2';
	$mqRepIon{'DimethNter6'}{'INTERIM_NAME'}='Dimethyl:2H(6)';
	$mqRepIon{'DimethNter8'}{'INTERIM_NAME'}='Dimethyl:2H(6)13C(2)';
	$mqRepIon{'ICAT-0'}{'INTERIM_NAME'}='ICAT_light';
	$mqRepIon{'ICAT-9'}{'INTERIM_NAME'}='ICAT_heavy';
	$mqRepIon{'ICPL-Lys0'}{'INTERIM_NAME'}='ICPL_light';
	$mqRepIon{'ICPL-Lys10'}{'INTERIM_NAME'}='ICPL_medium';
	$mqRepIon{'ICPL-Lys4'}{'INTERIM_NAME'}='ICPL_heavy';
	$mqRepIon{'ICPL-Lys6'}{'INTERIM_NAME'}='ICPL:13C(6)2H(4)';
	$mqRepIon{'ICPL-Nter0'}{'INTERIM_NAME'}='ICPL_light';
	$mqRepIon{'ICPL-Nter10'}{'INTERIM_NAME'}='ICPL_medium';
	$mqRepIon{'ICPL-Nter4'}{'INTERIM_NAME'}='ICPL_heavy';
	$mqRepIon{'ICPL-Nter6'}{'INTERIM_NAME'}='ICPL:13C(6)2H(4)';
	$mqRepIon{'Ile7'}{'INTERIM_NAME'}='Label:13C(6)15N(1)';
	$mqRepIon{'Leu7'}{'INTERIM_NAME'}='Label:13C(6)15N(1)';
	$mqRepIon{'Lys4'}{'INTERIM_NAME'}='Lys4';
	$mqRepIon{'Lys6'}{'INTERIM_NAME'}='13C6';
	$mqRepIon{'Lys8'}{'INTERIM_NAME'}='13C6-15N2 ';
	$mqRepIon{'mTRAQ-Lys0'}{'INTERIM_NAME'}='mTRAQ';
	$mqRepIon{'mTRAQ-Lys4'}{'INTERIM_NAME'}='mTRAQ:13C(3)15N(1)';
	$mqRepIon{'mTRAQ-Lys8'}{'INTERIM_NAME'}='mTRAQ:13C(6)15N(2)';
	$mqRepIon{'mTRAQ-Nter0'}{'INTERIM_NAME'}='mTRAQ';
	$mqRepIon{'mTRAQ-Nter4'}{'INTERIM_NAME'}='mTRAQ:13C(3)15N(1)';
	$mqRepIon{'mTRAQ-Nter8'}{'INTERIM_NAME'}='mTRAQ:13C(6)15N(2)';
	return(%mqRepIon);
}


####>Revision history<####
# 2.0.8 [FEATURE] "Only identified by site" exclusion & more commits and disconnects to prevent database lock timeout (PP 26/05/21)
# 2.0.7 [BUGFIX] Change getProtInfo call to add cluster parameter (VL 11/05/21)
# 2.0.6 [UPDATE] Records measures used for protein intensity quantifications (PP 04/12/20)
# 2.0.5 [UPDATE] Uses SQLite database to store protein quantification data (PP 07/07/20)
# 2.0.4 [BUGFIX] Ghost status is now only based on score=0 and no longer on lower-scoring status (PP 27/04/20)
# 2.0.3 [BUGFIX] Minor fix to define contaminant DB parse rules (PP 12/02/20)
# 2.0.2 [CHANGE] Minor change in $contaminantDB initialization (PP 10/01/20)
# 2.0.1 [BUGFIX] Right number of parameterGroups detected in mqpar.xml & good 2nd+ PTM position in peptides (PP 25/10/19)
# 2.0.0 [FEATURE] Multi-parameter groups support, smarter SILAC ghost peptides creation & [BUGFIX] on data files symbolic links (PP 23/08/19)
# 1.0.2 Checks for non-numerical mass error in evidence.txt (PP 10/07/19)
# 1.0.1 Minor changes (PP 11/07/19)
# 1.0.0 Forked from importMaxquant.cgi to handle import in background/cluster (PP ../05/19)
