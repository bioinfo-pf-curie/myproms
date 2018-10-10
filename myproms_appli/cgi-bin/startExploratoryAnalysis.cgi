#!/usr/local/bin/perl -w

################################################################################
# startExploratoryAnalysis.cgi       2.1.8                                    #
# Authors: P. Poullet, S.Liva (Institut Curie)                                 #
# Contact: myproms@curie.fr                                                    #
# run PCA or Clustering exploratory analysis                                   #
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
use strict;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use promsQuantif;
use POSIX qw(strftime); # to get the time
use File::Copy;
use Archive::Tar;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my $dbh = &promsConfig::dbConnect;
my %proteinQuantifFamilies=&promsQuantif::getProteinQuantifFamilies;## For fetchQuantificationData function

my $userID=$ENV{'REMOTE_USER'};
my %pepTypeDesc=('NUM_PEP_USED'=>'All','DIST_PEP_USED'=>'Distinct','RAZ_UNI_PEP'=>'Razor + unique','UNIQUE_PEP'=>'Unique','IDENT_PEP'=>'Identified');
my $experimentID=param('ID');
my $action=param('ACT') || 'explorAna'; # or export
my $projectID=&promsMod::getProjectID($dbh,$experimentID,'EXPERIMENT');
my $ajax=param('AJAX') || "";

if ($ajax eq 'displaySelect') {
	my $fromMotif= (param('CALL'))? param('CALL') : "";
	&ajaxDisplayData(param('MODIF'),param('QUANTIF_FAM'),$fromMotif);
	exit;
}
elsif ($ajax eq 'copySelect') {
    &ajaxCopySelectionFromExpAna(param('EXPANA_ID'));
    exit;
}
elsif ($ajax eq 'propAnnotate') {
	&ajaxPropAnnotate;
	exit;
}
if (param('submitted')) {

    print header(-charset=>'utf-8', -'content-encoding' => 'no');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR><BR><BR><BR><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">
<BR><FONT class="title3">Fetching data. Please wait...|;

    ### Get Parameter for insert
    my $PCAName = param('PCAName') || undef;
    my $clusterName = param('ClusterName') || undef;
    my $valMod = param('modif') || 0;
    my ($quantifFam,$quantifMeasCode)=split(':',param('quantifFam'));
	my $quantifParamID=0; # only for non-ratio/multi-measure quantifs
    if ($quantifMeasCode) {
	    ($quantifParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER QP,QUANTIFICATION_METHOD QM WHERE QM.CODE='$quantifFam' AND QM.ID_QUANTIFICATION_METHOD=QP.ID_QUANTIFICATION_METHOD AND QP.CODE='$quantifMeasCode'");
    }
    else {$quantifMeasCode='';}

    my $listParam = param('list') || "";
    my $featureType = param('feature') || "";
    my $dataTransform=param('dataTransform') || 'NONE';
    my %transform=(NONE=>sub {return $_[0];},
		   LOG2=>sub {return log($_[0])/log(2);},
		   LOG10=>sub {return log($_[0])/log(10);}
			       );
    my @selQuantifs=param('test_quantif');
	my @groupQuantifs=param('group_quantif');

	my $isNormalized=0;
    if ($valMod < 0) {
		my @normQuantifs=param('ref_quantif');
		foreach my $i (0..$#selQuantifs) {
			$selQuantifs[$i].='%'.$normQuantifs[$i];
		}
		$isNormalized=1;
    }

    my $strgParamList="//feature=$featureType//quantFam=$quantifFam";
    $strgParamList.="//quantCode=$quantifMeasCode" if $quantifMeasCode;
    $strgParamList.="//dataTrans=$dataTransform";
    my ($metric,$itemMetric,$method);
    if ($clusterName) {
		($itemMetric,$metric) = split(/:/,param('metric'));
		$method = param('method');
    }
    my $missingValues=param('missingValues');
    my $topNvarProteins=param('topNvar') || 0;
    my $restrictListID = ($listParam)? $listParam : undef;
    my $condRestrict = ($listParam)? param('condRestrict') : "";
    my $catExclusion = ($condRestrict eq 'exclude')? 1 : 0;
    my $numPepType=param('pepType') || 'NUM_PEP_USED';
    my $minNumPeptides=param('numPep') || 1;
#>RATIO-specific parameters
    my $maxFoldChange = param('foldChange') || 1;
    my $minNumOkRatios = param('foldChangeOcc') || 1;
    my $maxPValue = param('pValue') || 1;
    my $minNumOkPValues = param('pValueOcc') || 1;
    my $infRatios=param('infRatios') || 0;
    my $minFoldChange = 1/$maxFoldChange;

####create parameter file for Analysis pipeline R scirpt
#>FILTERING
    my $delocalize = param('delocalize')? "TRUE" : "FALSE";
    my $geneNameR="FALSE";
	my $useGeneName=param('geneName') || 0;
	my $aggregate=param('aggregate')? "TRUE" : "FALSE";
	my $pValueAnovaChk=param('pvalChkAnova')? "TRUE" : "FALSE";
    my $pValueAnova=($pValueAnovaChk eq "TRUE")? param('pValueAnova') : "FALSE";
    my $nbProtChk=param('topNvarChk')? "TRUE" : "FALSE";
    my $nbProt=($nbProtChk eq "TRUE")? param('topNvar') : "FALSE";
#>STAT-specific parameters
	my $exportRprocessMatrices=($nbProtChk eq 'TRUE' || $aggregate eq 'TRUE')? 1 : 0; # default
    my $anova = (param('anovaSel'))? param('anovaSel') : (param('topNvarChk'))? "sample" : "none";
	
	if ($aggregate eq "TRUE" && $nbProtChk eq "FALSE") {
		$anova="sample";
		if ($action eq "export") {
			$nbProtChk="TRUE";
			$nbProt=-1;
		}
	}
	
	my ($groupList,$groupName,%groupQuantif);
    if ($anova eq "group") {
		$groupList=param('GROUP_LIST');
		$groupName=param('GROUP_ANNOT');
		my %groupQuantifNames=&promsQuantif::getDistinctQuantifNames($dbh,join(':',param('test_quantif')));
		foreach my $list (split('::',$groupList)) {
			my ($group, $quantifStrg)=split('=', $list);
			foreach my $quantif (split('\|', $quantifStrg)) {
				(my $quantifName=$groupQuantifNames{'EXTENDED'}{$quantif})=~s/\W+/\./g;
				$quantifName=~s/^\.+//; $quantifName=~s/\.+\Z//;
				push @{$groupQuantif{$group}},$quantifName;
			}
		}
    }

    my $strgFilterList="//PEP=$numPepType:$minNumPeptides//NA=$missingValues";
    $strgFilterList.="//TOPN_VAR=$topNvarProteins" if $topNvarProteins;
    if ($quantifFam eq 'RATIO') {
		$strgFilterList.="//PV=$maxPValue:$minNumOkPValues" if $valMod >= 0; # not for normalized data
		$strgFilterList.="//FC=$maxFoldChange:$minNumOkRatios//INF=$infRatios";
    }

    my (@selectedQuantifications, @testQuantifications, %quantifPos, %test2refQuantif, %isRefQuantif);
    my $strgNorm;
    foreach my $quantif (@selQuantifs) {
		my $testQuant=$quantif;
		if ($quantif=~/%/) { #test%ref
			$strgNorm.=':' if $strgNorm;
			$strgNorm.=$quantif;
			($testQuant, my $refQuant) = split(/%/,$quantif);
			$test2refQuantif{$testQuant}=$refQuant;
			my ($quantifRefID, $targetRefPos) = split("_",$refQuant);
			$targetRefPos=0 unless $targetRefPos;
			push @{$quantifPos{$quantifRefID}},$targetRefPos;
			$isRefQuantif{$refQuant}=1;
		}
		my ($quantifID,$targetPos) = split("_",$testQuant);
		$targetPos=0 unless $targetPos;
		push @{$quantifPos{$quantifID}},$targetPos;
		push @selectedQuantifications,$testQuant;
		push @testQuantifications,$testQuant;
    }
    $strgNorm='normalization='.$strgNorm if $strgNorm;
    foreach my $refQuant (keys %isRefQuantif) {push @selectedQuantifications,$refQuant;}

    my ($refSelectedProteins,$refExcludedProteins,$nbProtFilter);
    if ($restrictListID) {
		my %listProteins;
		my $sthCP = $dbh->prepare("SELECT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$restrictListID");
		$sthCP->execute;
		while (my ($protID) = $sthCP->fetchrow_array) {
			$listProteins{$protID} = 1;
		}
		$sthCP -> finish;
		($nbProtFilter)=scalar(keys %listProteins);
		$strgFilterList.="//nbProt=$nbProtFilter";
		($refSelectedProteins,$refExcludedProteins)=($condRestrict eq 'restrict')? (\%listProteins,undef) : (undef,\%listProteins);
    }

#####################################>TEST GO TERM 20/08/15
#if (1==2) {
#	my @goAnalyses=(89,90,91); # gr/SILAC-slim
#	#my @goAnalyses=(83,87,88); # gr/gr-slim
#	my $goID=undef;
#	my $aspect='P';
#	my %goProtList;
#	foreach my $goAnaID (@goAnalyses) {
#		my $resultDirUnix="$promsPath{go_unix}/project_$projectID/$goAnaID";
#		open(RES,"$resultDirUnix/results_$aspect.txt") || die "Could not open result file";
#		while(<RES>){
#			next if $_ =~ /^#/;
#			chomp;
#			my @info = split(/\t/, $_);
#			next if $info[0] =~ /unannotated|discarded/;
#			next if $info[5] eq '.';
#			if(!$goID || $info[0] eq $goID){
#				foreach my $protID (split ';',$info[5]){
#					$goProtList{$protID}=1;
#				}
#				last if $goID;
#			}
#		}
#		close RES;
#	}
#	$refSelectedProteins=\%goProtList;
#}
#####################################<<<<<<<<<<<<<<<<<<<<<

	my ($pcaID, $clusteringID,$sthUpdateFilterList);
	my ($exportDir,$exportPath,$quantifFile,$pepFile,$pepFileRef);
	my @file2Tar;
	if ($action eq 'explorAna') {

		####INSERT database
		my $sthInsertExplorAna = $dbh -> prepare("INSERT INTO EXPLORANALYSIS(ID_EXPLORANALYSIS,ID_CATEGORY,ID_EXPERIMENT,NAME,ANA_TYPE,PARAM_LIST,FILTER_LIST,CAT_EXCLUSION,STATUS,UPDATE_DATE,UPDATE_USER) VALUES (?,?,?,?,?,?,?,?,-1,NOW(),?)");
		#$pcaID = $clusteringID = 1; ## tocomment

		####Update Filter_list with nbProt
		#$sthUpdateFilterList=$dbh->prepare("UPDATE EXPLORANALYSIS SET FILTER_LIST=CONCAT(FILTER_LIST,?) WHERE ID_EXPLORANALYSIS=?");

		$strgParamList.="//$strgNorm" if $strgNorm;

		##PCA
		if (defined $PCAName) {
			($pcaID) = $dbh->selectrow_array("SELECT MAX(ID_EXPLORANALYSIS) FROM EXPLORANALYSIS");
			$pcaID++;
			$sthInsertExplorAna -> execute($pcaID,$restrictListID,$experimentID,$PCAName,'PCA',$strgParamList,$strgFilterList,$catExclusion,$userID);
			#$dbh -> commit;
		}

		##Clustering
		if (defined $clusterName) {
			$strgParamList .= "//method=$method//metric=$metric";
			($clusteringID) =  $dbh->selectrow_array("SELECT MAX(ID_EXPLORANALYSIS) FROM EXPLORANALYSIS");
			$clusteringID++;
			$sthInsertExplorAna -> execute($clusteringID,$restrictListID,$experimentID,$clusterName,'cluster',$strgParamList,$strgFilterList,$catExclusion,$userID);
			#$dbh -> commit;
		}
		$sthInsertExplorAna -> finish;

		##EXPLORANA_QUANTIF
		my $sthInsertExplorQuantif=$dbh->prepare("INSERT INTO EXPLORANA_QUANTIF(ID_QUANTIFICATION,ID_EXPLORANALYSIS,TARGET_POS) VALUES (?,?,?)");
		my %explorAnaQuatif;
		foreach my $explorAnaID ($pcaID,$clusteringID) {
			next if (!$explorAnaID);
			foreach my $quantifID (keys %quantifPos) {
				foreach my $targetPos (@{$quantifPos{$quantifID}}) {
					my $entryKey="$quantifID,$explorAnaID,$targetPos";
					next if $explorAnaQuatif{$entryKey}; # prevents duplicate entries
					$sthInsertExplorQuantif->execute($quantifID,$explorAnaID,$targetPos);
					$explorAnaQuatif{$entryKey}=1;
				}
			}
		}
		$dbh->commit;
		$sthInsertExplorQuantif->finish;

		###>Config<####
		mkdir "$promsPath{tmp}" unless -e "$promsPath{tmp}";
		mkdir "$promsPath{explorAna}/project_$projectID" unless -e "$promsPath{explorAna}/project_$projectID";
		mkdir "$promsPath{tmp}/exploratory_analysis" unless -e "$promsPath{tmp}/exploratory_analysis";
		my $numMissValue=scalar(@testQuantifications)*$missingValues/100;
		###>Directory for PCA
		if ($pcaID) {
			mkdir "$promsPath{tmp}/exploratory_analysis/$pcaID" unless -e "$promsPath{tmp}/exploratory_analysis/$pcaID";
			open(MAT_PCA,">$promsPath{tmp}/exploratory_analysis/$pcaID/matrix.txt");
			open(NA_PCA,">$promsPath{tmp}/exploratory_analysis/$pcaID/missingValues.txt");


			open(R_PARAMS,">$promsPath{tmp}/exploratory_analysis/$pcaID/R_parameters.txt");
			print R_PARAMS "EXCLUDE_AMB\tGENE_NAME\tPROTEIN_SELECTION\tAGGREGATE\tKEEP_PROT\tKEEP_PROT_NB\tP_VALUE_CHK\tP_VALUE\tMISSING_VALUE\tINF_RATIO\n";
			print R_PARAMS "$delocalize\t$geneNameR\t$anova\t$aggregate\t$nbProtChk\t$nbProt\t$pValueAnovaChk\t$pValueAnova\t$numMissValue\t$infRatios\n";
			close(R_PARAMS);
		}

		###>Directory for CLUSTER
		if ($clusteringID) {
			mkdir "$promsPath{tmp}/exploratory_analysis/$clusteringID" unless -e "$promsPath{tmp}/exploratory_analysis/$clusteringID";
			open(MAT_CLUSTER,">$promsPath{tmp}/exploratory_analysis/$clusteringID/matrix.txt");
			open(NA_CLUSTER,">$promsPath{tmp}/exploratory_analysis/$clusteringID/missingValues.txt");

			open(R_PARAMS,">$promsPath{tmp}/exploratory_analysis/$clusteringID/R_parameters.txt");
			print R_PARAMS "EXCLUDE_AMB\tGENE_NAME\tPROTEIN_SELECTION\tAGGREGATE\tKEEP_PROT\tKEEP_PROT_NB\tP_VALUE_CHK\tP_VALUE\tMISSING_VALUE\tINF_RATIO\n";
			print R_PARAMS "$delocalize\t$geneNameR\t$anova\t$aggregate\t$nbProtChk\t$nbProt\t$pValueAnovaChk\t$pValueAnova\t$numMissValue\t$infRatios\n";
			close(R_PARAMS);
		}
	}
	else { # export
		###>Config<####
		#mkdir "$promsPath{tmp}" unless -e "$promsPath{tmp}";
		mkdir "$promsPath{tmp}/scratch" unless -d "$promsPath{tmp}/scratch";
		&promsMod::cleanDirectory("$promsPath{tmp}/scratch/export",'15m');
		my $jobID=strftime("%Y%m%d%H%M%S",localtime);
		$exportDir=$jobID."_$userID";
		mkdir "$promsPath{tmp}/scratch/export/$exportDir";
		$exportPath="$promsPath{tmp}/scratch/export/$exportDir";
		$quantifFile=($quantifFam eq 'RATIO')? 'ratio.txt' : "$quantifMeasCode.txt"; # exported intensities are NOT logged
		push @file2Tar, ($exportRprocessMatrices)? "processed_$quantifFile" : $quantifFile;
		open(MAT_EXPORT_QUANTIF,">$exportPath/$quantifFile");
		if ($quantifFam eq 'RATIO' && !$isNormalized) {
			$pepFile=($numPepType eq 'NUM_PEP_USED')? 'peptide.txt' : 'dist_peptide.txt';
			open(MAT_EXPORT_PEP,">$exportPath/$pepFile");
			open(MAT_EXPORT_PVAL,">$exportPath/pvalue.txt");
			#my $pepFile2Tar = ($numPepType eq 'NUM_PEP_USED')? "matrix_processed.txt" : "matrix_processed_$pepFile";
			my $pepFile2Tar = "processed_$pepFile";
			push @file2Tar, ($exportRprocessMatrices)? ("processed_$pepFile","processed_pvalue.txt") : ($pepFile,'pvalue.txt');
		}
		push @file2Tar, ($exportRprocessMatrices)? "processed_annotation.txt" : "annotation.txt";
		push @file2Tar, "parameters.txt";
		open(ANNOT,">$exportPath/annotation.txt");
		open(PARAMS,">$exportPath/parameters.txt");
		my $focusStrg='Proteins';
		if ($valMod) {
		    my ($modifName)=$dbh->selectrow_array("SELECT PSi_MS_NAME FROM MODIFICATION WHERE ID_MODIFICATION=ABS($valMod)");
		    $focusStrg=$modifName.'-Proteins';
		}
		if ($exportRprocessMatrices) {
			push @file2Tar, "R_parameters.txt";
			my $numMissValue=scalar(@testQuantifications)*$missingValues/100;
			open(R_PARAMS,">$exportPath/R_parameters.txt");
			print R_PARAMS "MODIF\tQUANTIF_FAM\tPROTEIN_SELECTION\tAGGREGATE\tKEEP_PROT\tKEEP_PROT_NB\tP_VALUE_CHK\tP_VALUE\tMISSING_VALUE\n";
			print R_PARAMS "$valMod\t$quantifFam\t$anova\t$aggregate\t$nbProtChk\t$nbProt\t$pValueAnovaChk\t$pValueAnova\t$numMissValue\n";
			close(R_PARAMS);
		}

		my $listFilterStrg='None';
		if ($restrictListID) {
			my ($listName)=$dbh->selectrow_array("SELECT CONCAT(CL.NAME,' > ',CA.NAME) FROM CATEGORY CA,CLASSIFICATION CL WHERE CA.ID_CLASSIFICATION=CL.ID_CLASSIFICATION AND CA.ID_CATEGORY=$restrictListID");
			$listFilterStrg="$condRestrict: $listName ($nbProtFilter proteins)";
		}
		my $infRatioStrg=($infRatios < 0)? '=NA' : $infRatios;

		print PARAMS qq
|FOCUS\t$focusStrg
QUANTIF_METHOD\t$quantifFam
DATA_TRANSFORM=$dataTransform
%_MISSING_VALUES_ALLOWED\t$missingValues
LIST_FILTERING\t$listFilterStrg
MIN_NUM_PEPTIDES\t$minNumPeptides
PEPTIDE_TYPE\t$pepTypeDesc{$numPepType}
DATA_SELECTION\t$anova
|;
		if ($anova eq "group") {
			print PARAMS "GROUP\t$groupName\n";
		}
		if ($nbProtChk eq "TRUE") {
			print PARAMS "PROTEIN_SIGN\t$nbProt\n";
		}
		if ($pValueAnovaChk eq "TRUE") {
			print PARAMS "PVALUE_SIGN\t$pValueAnova\n";
		}
		if ($valMod) {
			print PARAMS qq
|AGGREGATE\t$aggregate
EXCLUDE_AMB\t$delocalize
|;
		}
		if ($quantifFam eq 'RATIO') {
			my $normStrg=($isNormalized)? 'Yes' : 'No';
			print PARAMS qq
|NORMALIZED_RATIOS\t$normStrg
MAX_ABS_RATIO\t$maxFoldChange
MIN_RATIO_OCCURENCE\t$minNumOkRatios
MINUS_INFINITE_RATIO\t-1000
PLUS_INFINITE_RATIO\t1000
MAX_P_VALUE\t$maxPValue
MIN_P_VALUE_OCCURENCE\t$minNumOkPValues
%_INFINITE_RATIOS_ALLOWED\t$infRatioStrg
|;
		}
		elsif ($quantifMeasCode) {
			my $quantifMeasName;
			foreach my $refMeas (@{$proteinQuantifFamilies{'MEASURES'}{$quantifFam}}) {
				if ($refMeas->[0] eq $quantifMeasCode) {
					$quantifMeasName=$refMeas->[1];
					last;
				}
			}
			print PARAMS "QUANTIF_CODE\t$quantifMeasName\n";
		}
		#elsif ($quantifFam eq 'EMPAI') {
		#	print PARAMS "QUANTIF_CODE=$empaiQuantifCode\n";
		#}
		close PARAMS;
	}

	my (%quantifValues,%quantifInfo,%proteinInfo, %groupQuantifNames);
	my ($view,$refProtInfo)=($action eq 'explorAna')? ('explorAna',undef) : ('export',\%proteinInfo);
	my $absValMod = ($isNormalized)? abs($valMod)  : $valMod;
    my %parameters=(QUANTIF_FAMILY=>$quantifFam,VIEW=>$view,NUM_PEP_CODE=>$numPepType,QUANTIF_LIST=>\@selectedQuantifications,SEL_MODIF_ID=>$absValMod,VERBOSE=>1);
	$parameters{MEASURE}=$quantifMeasCode if $quantifMeasCode;
	&promsQuantif::fetchQuantificationData($dbh,\%parameters,\%quantifInfo,\%quantifValues,$refProtInfo,$refSelectedProteins,$refExcludedProteins);
	my %quantifSoftwares;
	if ($quantifFam eq 'RATIO') {
		foreach my $quantif (@selectedQuantifications) {
			my ($quantifID,$targetPos)=split('_',$quantif);
			$quantifSoftwares{$quantif}=($quantifInfo{$quantifID}[1]->{SOFTWARE})? $quantifInfo{$quantifID}[1]->{SOFTWARE}[0] : 'myProMS';
			$quantifSoftwares{$quantif}='MaxQuant' if $quantifSoftwares{$quantif} eq 'MQ';
		}
	}

	my @header;
	if ($action eq 'explorAna') {
		foreach my $testQuant (@testQuantifications) {
			if ($test2refQuantif{$testQuant}) {push @header,$testQuant.'%'.$test2refQuantif{$testQuant};}
			else {push @header,$testQuant;}
		}
	}
	else { # export

		my %testQuantifNames=&promsQuantif::getDistinctQuantifNames($dbh,join(':',param('test_quantif')));
		my %refQuantifNames=($isNormalized)? &promsQuantif::getDistinctQuantifNames($dbh,join(':',param('ref_quantif'))) : ();

		foreach my $quantif (@testQuantifications) {
			(my $headerStrg=$testQuantifNames{'EXTENDED'}{$quantif})=~s/\W+/\./g;
			$headerStrg=~s/^\.+//; $headerStrg=~s/\.+\Z//;
			my $strgHeader=$headerStrg;
			$groupQuantifNames{$strgHeader}=();
			if ($isNormalized) {
				(my $refStrg=$refQuantifNames{'EXTENDED'}{$test2refQuantif{$quantif}})=~s/\W+/\./g;
				$refStrg=~s/^\.+//; $refStrg=~s/\.+\Z//;
				$headerStrg.='.VS.'.$refStrg;

			}
			$groupQuantifNames{$strgHeader}=$headerStrg;
			push @header,$headerStrg;
		}
	}

	if ($anova eq "group") {
		push @file2Tar, "group.txt";
		push @file2Tar, "anova_pvalue.txt";
		open(GROUP,">$exportPath/group.txt");
		print GROUP "sample\tgroup_snf\n";
		foreach my $group (sort{$a cmp $b} keys %groupQuantif) {
			foreach my $quantif (@{$groupQuantif{$group}}) {
				print GROUP "$groupQuantifNames{$quantif}\t$group\n";
			}
		}
		close(GROUP);
	}

	if ($anova eq "sample") {
		push @file2Tar, "sd.txt";
	}

	my $strgCond = join("\t",@header);
	my $numTestQuantifs=scalar @testQuantifications;
	my $numAllowedInf=($infRatios < 0)? $numTestQuantifs : $numTestQuantifs*$infRatios/100;
	my $numAllowedMissing=$numTestQuantifs*$missingValues/100;

	print MAT_PCA "\t$strgCond\n" if $pcaID;
	print MAT_CLUSTER "\t$strgCond\n" if $clusteringID;
	if ($action eq 'export') {
		print MAT_EXPORT_QUANTIF "\t$strgCond\n";
		if ($quantifFam eq 'RATIO' && !$isNormalized) {
			print MAT_EXPORT_PEP "\t$strgCond\n";
			print MAT_EXPORT_PVAL "\t$strgCond\n";
		}
		print ANNOT "PROTEIN\tIDENTIFIER\tGENE\tSYNONYMS\tDESCRIPTION\n";
	}

	print '/.';

	my $nbAllProt=0;
	my $count=0;
	my (%usedProteinInfo,%missingValuesProt);
	foreach my $modProtID (keys %quantifValues) { # skip sort because of modif quantif -> sort{lc($proteinInfo{$a}[0]) cmp lc($proteinInfo{$b}[0])} # sort by ALIAS
		my ($proteinID,$modStrg)=($modProtID=~/^(\d+)(.*)/);
		$modStrg='' unless $modStrg;
		#>Scan for +/- inf & missing values & peptide filter
		my (%infiniteRatios,%missingValueQuantif,%usedProtRatio); #,@trueValues
		my $numTrueValues=0;
		foreach my $quantif (@testQuantifications) {
			my $strgQuantif= ($isNormalized)? $quantif."%".$test2refQuantif{$quantif} : $quantif;
			$count++;
			if ($count==50000) {
				$count=0;
				print '.';
			}
			if ($quantifFam eq 'RATIO' && $test2refQuantif{$quantif}) { # no p-value filter if ratio normalization (final ratio is unrelated to recorded p-value!); +/-Inf ratio not allowed for reference
				if  (!$quantifValues{$proteinID} || !$quantifValues{$proteinID}{$test2refQuantif{$quantif}}{'RATIO'} || $quantifValues{$proteinID}{$test2refQuantif{$quantif}}{'RATIO'}==0.001 || $quantifValues{$proteinID}{$test2refQuantif{$quantif}}{'RATIO'}==1000 || !$quantifValues{$proteinID}{$test2refQuantif{$quantif}}{$numPepType} || $quantifValues{$proteinID}{$test2refQuantif{$quantif}}{$numPepType} < $minNumPeptides) {
					$missingValueQuantif{$strgQuantif}=1;
					next;
				}
			}
			if ($quantifValues{$modProtID}{$quantif}) {
				if (!$quantifValues{$modProtID}{$quantif}{$numPepType} || $quantifValues{$modProtID}{$quantif}{$numPepType} < $minNumPeptides) { # peptide filters
					$missingValueQuantif{$strgQuantif}=1;
					next;
				}
				if ($quantifFam eq 'RATIO') {
					if ($quantifValues{$modProtID}{$quantif}{'RATIO'} == 0.001) {
						if ($infRatios < 0) { # use as missing value
							$missingValueQuantif{$strgQuantif}=1;
						}
						else {
							$infiniteRatios{$quantif}=-1;
							$usedProtRatio{$quantif}=-1000;
						}
					}
					elsif ($quantifValues{$modProtID}{$quantif}{'RATIO'} == 1000) {
						if ($infRatios < 0) { # use as missing value
							$missingValueQuantif{$strgQuantif}=1;
						}
						else {
							$infiniteRatios{$quantif}=1;
							$usedProtRatio{$quantif}=1000;
						}
					}
					else { # usable value in test (and ref)
						$usedProtRatio{$quantif}=($test2refQuantif{$quantif})? $quantifValues{$modProtID}{$quantif}{'RATIO'}/$quantifValues{$proteinID}{$test2refQuantif{$quantif}}{'RATIO'} : $quantifValues{$modProtID}{$quantif}{'RATIO'};
						$numTrueValues++;
# TEMP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#$usedProtRatio{$quantif}=$quantifValues{$modProtID}{$quantif}{'NUM_PEP_USED'};
# <-------------
						#my $absRatio=($usedProtRatio{$quantif} < 1)? 1/$usedProtRatio{$quantif} : $usedProtRatio{$quantif};
						#$maxTrueRatio=$absRatio if $absRatio > $maxTrueRatio;
					}
				}
				else { # MQ/SIN/EMPAI
					$usedProtRatio{$quantif}=$quantifValues{$modProtID}{$quantif}{$quantifMeasCode};
					$numTrueValues++;
				}
			}
			else {$missingValueQuantif{$strgQuantif}=1;}
		}

		#my $numTrueValues=scalar @trueValues;
		my $numProtMissValues=scalar keys %missingValueQuantif;
		next if ($numProtMissValues > $numAllowedMissing || ($quantifFam eq 'RATIO' && scalar (keys %infiniteRatios) > $numAllowedInf) || ($action eq 'explorAna' && $numTrueValues < 2));
		%{$missingValuesProt{$modProtID}}=%missingValueQuantif if $numProtMissValues;

		#$maxTrueRatio*=2; # !!!!!!!!!!!!!!!!!!! +/-inf ratios converted into 2x max ratio for protein accross all quantifs !!!!!!!!!!!!!!!!!

		#>Re-scan for ratio & p-value filter
		if ($quantifFam eq 'RATIO') {
			my ($okRatio,$okPvalue)=(0,0);
			foreach my $quantif (@testQuantifications) {
				next if $missingValueQuantif{$quantif};
				$okRatio++ if ($maxFoldChange==1 || $usedProtRatio{$quantif} <= $minFoldChange || $usedProtRatio{$quantif} >= $maxFoldChange);
				$okPvalue++ if ($maxPValue==1 || $quantifSoftwares{$quantif} eq 'MaxQuant' || ($quantifValues{$modProtID}{$quantif}{'P_VALUE'} && $quantifValues{$modProtID}{$quantif}{'P_VALUE'} <= $maxPValue)); # !!!Does not rely on +/-inf to validate p-value filter!!! (+/-inf accepted only if another quantif has validated filter)
			}
			next if ($okRatio < $minNumOkRatios || $okPvalue < $minNumOkPValues);
		}

		$nbAllProt++;
		next if ($delocalize eq "TRUE" && $modStrg=~/:/);
		print MAT_PCA "$modProtID" if $pcaID;
		print MAT_CLUSTER "$modProtID" if $clusteringID;
		if ($action eq 'export') {
			my $entryName;
			if ($useGeneName) {$entryName=$proteinInfo{$proteinID}[5][0] || 'NO_GENE';} # gene
			else {$entryName=$proteinInfo{$proteinID}[0];} # ALIAS
			print MAT_EXPORT_QUANTIF "$entryName(uid$proteinID)$modStrg";
			if ($quantifFam eq 'RATIO' && !$isNormalized) {
				print MAT_EXPORT_PEP "$entryName(uid$proteinID)$modStrg";
				print MAT_EXPORT_PVAL "$entryName(uid$proteinID)$modStrg";
			}
		}
		foreach my $quantif (@testQuantifications) {
			if ($usedProtRatio{$quantif}) { # usable value
				if ($quantifFam eq 'RATIO') {
					if ($infiniteRatios{$quantif}) { # flag +/-inf ratios
						print MAT_PCA "\t$usedProtRatio{$quantif}" if $pcaID; # 1000 or -1000;
						print MAT_CLUSTER "\t$usedProtRatio{$quantif}" if $clusteringID;
						if ($action eq 'export') {
							print MAT_EXPORT_QUANTIF "\t$usedProtRatio{$quantif}";
							print MAT_EXPORT_PVAL "\tNA" unless $isNormalized;
						}
					}
					else {
						my $transValue=$transform{$dataTransform}->($usedProtRatio{$quantif});
						print MAT_PCA "\t$transValue" if $pcaID;
						print MAT_CLUSTER "\t$transValue" if $clusteringID;
						if ($action eq 'export') {
							print MAT_EXPORT_QUANTIF "\t$transValue";
							unless ($isNormalized) {
								if ($quantifValues{$modProtID}{$quantif}{'P_VALUE'}) {print MAT_EXPORT_PVAL "\t$quantifValues{$modProtID}{$quantif}{P_VALUE}";}
								else {print MAT_EXPORT_PVAL "\tNA";}
							}
						}
					}
					print MAT_EXPORT_PEP "\t$quantifValues{$modProtID}{$quantif}{$numPepType}" if ($action eq 'export' && !$isNormalized);
				}
				else { # Non-ratio quantifs MQ,SIN,EMPAI
					my $transValue=$transform{$dataTransform}->($usedProtRatio{$quantif});
					print MAT_PCA "\t$transValue" if $pcaID;
					print MAT_CLUSTER "\t$transValue" if $clusteringID;
					if ($action eq 'export') {print MAT_EXPORT_QUANTIF "\t$transValue";}
				}
			}
			else { # missing value
				if ($action eq 'explorAna') {
					print MAT_PCA "\tNA" if $pcaID;
					print MAT_CLUSTER "\tNA" if $clusteringID;
				}
				else {
					print MAT_EXPORT_QUANTIF "\tNA";
					if ($quantifFam eq 'RATIO' && !$isNormalized) {
						print MAT_EXPORT_PEP "\tNA";
						print MAT_EXPORT_PVAL "\tNA";
					}
				}
			}
		}
		if ($action eq 'explorAna') {
			if ($pcaID) {
				print MAT_PCA "\n";
			}
			if ($clusteringID) {
				print MAT_CLUSTER "\n";
			}
		}
		else { # export
			print MAT_EXPORT_QUANTIF "\n";
			if ($quantifFam eq 'RATIO' && !$isNormalized) {
				print MAT_EXPORT_PEP "\n";
				print MAT_EXPORT_PVAL "\n";
			}
			unless ($usedProteinInfo{$proteinID}) {
				my ($gene,$synStrg)=('','');
				if ($proteinInfo{$proteinID}[5][0]) { # gene
					$gene=shift @{$proteinInfo{$proteinID}[5]};
					$synStrg=join(',',@{$proteinInfo{$proteinID}[5]});
				}
				if ($useGeneName) {
					my $protName=$gene || 'NO_GENE';
					print ANNOT $protName;
				}
				else {print ANNOT $proteinInfo{$proteinID}[0];}
				print ANNOT "(uid$proteinID)\t$proteinInfo{$proteinID}[1]\t$gene\t$synStrg\t$proteinInfo{$proteinID}[2]\n"; # [1]=IDENTIFIER
				$usedProteinInfo{$proteinID}=1;
			}
		}
	}
	
	#exit;
	if ($action eq 'explorAna') {
		close MAT_CLUSTER if ($clusteringID);
		close MAT_PCA if ($pcaID);
	}
	else { # export
		close MAT_EXPORT_QUANTIF;
		close ANNOT;
		if ($quantifFam eq 'RATIO' && !$isNormalized) {
			close MAT_EXPORT_PEP;
			close MAT_EXPORT_PVAL;
		}
	}

	###>No proteins in matrix<###
#$nbAllProt=0; # TEST
	unless ($nbAllProt) {
		my ($buttonLabel,$onClickAction);
		if ($action eq 'explorAna') {
			my $sthUpdateStatus = $dbh -> prepare("UPDATE EXPLORANALYSIS SET STATUS=-2 WHERE ID_EXPLORANALYSIS = ?"); # set to error status
			$sthUpdateStatus->execute($clusteringID) if $clusteringID;
			$sthUpdateStatus->execute($pcaID) if $pcaID;
			$dbh->commit;

			my $explorID=($clusteringID)? $clusteringID : $pcaID;
			$buttonLabel='Continue';
			$onClickAction="top.promsFrame.selectedAction = 'summary'; parent.itemFrame.location='$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=experiment:$experimentID&branchID=exploranalysis:$explorID&VIEW=explorAna'";
		}
		else { # export
			unlink glob "$exportPath/*.txt";
			rmdir $exportPath;
			$buttonLabel='Try again';
			$onClickAction="parent.itemFrame.selectOption()";
		}
		$dbh->disconnect;

		print qq
| Done.</FONT>
<BR><BR><BR><BR><BR>
<FONT class="title2" color="#DD0000">***ERROR: No proteins matched selected filters!***</FONT>
<BR><INPUT type="button" class="title2" value="$buttonLabel" onclick="$onClickAction">
<BR><BR><BR>
</CENTER>
</BODY>
</HTML>
|;
		exit;
	}

	$dbh->disconnect;
	if ($action eq 'explorAna') {
		my $explorID=$pcaID || $clusteringID;

		###>Precrossing data (Imputing +/-inf & missing values, filtering top sd prot)<###
		print " Done.<BR>\nPreprocessing data...";
		my $explorDIR = "$promsPath{tmp}/exploratory_analysis/$explorID";
		my $childExplor = fork;
		unless ($childExplor) { # child here
			#>Disconnecting from server
			open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
			open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
			open STDERR, '>/dev/null' or die "Can't open /dev/null: $!";
			#$topNvarProteins=0 if $topNvarProteins > $nbAllProt;
			system "cd $explorDIR; $promsPath{R}/R CMD BATCH --no-save --no-restore $promsPath{R_scripts}/prepareExplorAna.R";
			exit;
		}

		##>Waiting for R to finish<##
		my $wait = 1;
		my $errorTxt = '';
		my $count1 = 0;
		my $count2 = 0;
		while ($wait == 1) {
			sleep 2;
			print '.';
			if (!-e "$explorDIR/prepareExplorAna.Rout") {
				if ($count1 > 300) { # <=> 10 min
					$errorTxt='No R output file found';
					$wait=0;
				}
				$count1++;
				next;
			}
			my $Rprocess = `grep -c '> proc.time()' $explorDIR/prepareExplorAna.Rout`;
			chomp $Rprocess;
			if ($Rprocess) { # finished normally
				$wait = 0;
			}
			else {
				$Rprocess = `grep -c '^Execution halted' $explorDIR/prepareExplorAna.Rout`;
				chomp $Rprocess;
				if ($Rprocess) {
					$errorTxt='Execution halted from R';
					$wait = 0;
				}
			}
			$count2++;
			if ($count2 > 900) { # <=> 30 min
				$errorTxt='Process duration has exceeded 30 min.';
				$wait=0;
			}
		}
		if ($errorTxt) {
			print qq
|</FONT>
<BR><BR><BR><BR><BR>
<FONT class="title2" color="#DD0000">***ERROR: Data preprocessing failed ($errorTxt)!***</FONT>
<BR><INPUT type="button" class="title2" value="Continue" onclick="top.promsFrame.selectedAction = 'summary'; parent.itemFrame.location='$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=experiment:$experimentID&branchID=exploranalysis:$explorID&VIEW=explorAna'">
<BR><BR><BR>
</CENTER>
</BODY>
</HTML>
|;
			exit;
		}
		##>Updating number of proteins used
		my %filteredProtList;
		my $filtering=0;
		if ($topNvarProteins && $topNvarProteins < $nbAllProt) {
			$filtering=1;
			#my $infoStrg=`wc -l $explorDIR/matrixProcessed.txt`;
			#my ($numLines)=$infoStrg=~/^(\d+)/;
			#$nbAllProt=$numLines-1; # remove header line
			open(MAT_PROC,"$explorDIR/matrixProcessed.txt") or die $!;
			while (<MAT_PROC>) {
				next if $.==1;
				my ($modProtID)=$_=~/^(\S+)/;
				$filteredProtList{$modProtID}=1;
			}
			close MAT_PROC;
			$nbAllProt=scalar keys %filteredProtList; # overwrites $nbAllProt
		}

		##>Counting number of missing values (after filtering if any)
		my $numMissingValues=0;
		foreach my $modProtID (keys %missingValuesProt) {
			next if ($filtering && !$filteredProtList{$modProtID});
			$numMissingValues+=scalar keys %{$missingValuesProt{$modProtID}};
			print NA_PCA "$modProtID\t",join("\t",keys %{$missingValuesProt{$modProtID}}),"\n" if $pcaID;
			print NA_CLUSTER "$modProtID\t",join("\t",keys %{$missingValuesProt{$modProtID}}),"\n" if $clusteringID;
		}
		close NA_PCA if $pcaID;
		close NA_CLUSTER if $clusteringID;

		####Update Filter_list with nbProt
		$dbh = &promsConfig::dbConnect; # reconnect after R
		$sthUpdateFilterList=$dbh->prepare("UPDATE EXPLORANALYSIS SET FILTER_LIST=CONCAT(FILTER_LIST,?) WHERE ID_EXPLORANALYSIS=?");
		$sthUpdateFilterList->execute("//nbAllProt=$nbAllProt//nbNA=$numMissingValues//imputeNA=PCA",$pcaID) if $pcaID;
		$sthUpdateFilterList->execute("//nbAllProt=$nbAllProt//nbNA=$numMissingValues//imputeNA=PCA",$clusteringID) if $clusteringID;
		$sthUpdateFilterList->finish;
		$dbh->commit;
		$dbh->disconnect;

		##>Replacing matrices<##
		if ($pcaID) {
			#move("$promsPath{tmp}/exploratory_analysis/$pcaID/matrix.txt","$promsPath{tmp}/exploratory_analysis/$pcaID/matrixBefore.txt");
			unlink "$promsPath{tmp}/exploratory_analysis/$pcaID/matrix.txt";
			move("$promsPath{tmp}/exploratory_analysis/$pcaID/matrixProcessed.txt","$promsPath{tmp}/exploratory_analysis/$pcaID/matrix.txt");
			if ($clusteringID) {
				#move("$promsPath{tmp}/exploratory_analysis/$clusteringID/matrix.txt","$promsPath{tmp}/exploratory_analysis/$clusteringID/matrixBefore.txt");
				unlink "$promsPath{tmp}/exploratory_analysis/$clusteringID/matrix.txt";
				copy("$promsPath{tmp}/exploratory_analysis/$pcaID/matrix.txt","$promsPath{tmp}/exploratory_analysis/$clusteringID/matrix.txt");
			}
		}
		elsif ($clusteringID) { # no PCA
			#move("$promsPath{tmp}/exploratory_analysis/$clusteringID/matrix.txt","$promsPath{tmp}/exploratory_analysis/$clusteringID/matrixBefore.txt");
			unlink "$promsPath{tmp}/exploratory_analysis/$clusteringID/matrix.txt";
			move("$promsPath{tmp}/exploratory_analysis/$clusteringID/matrixProcessed.txt","$promsPath{tmp}/exploratory_analysis/$clusteringID/matrix.txt");
		}
		print " Done.</FONT>\n<BR><BR><BR><BR><BR>\n";

		###>Forking to launch Analyses<###
		if ($pcaID) {
			my $childPidPCA = fork;
			unless ($childPidPCA) { # child here
				#>Disconnecting from server
				open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
				open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
				open STDERR, '>/dev/null' or die "Can't open /dev/null: $!";
				system "./launchExploratoryAnalyses.pl $pcaID PCA $projectID 2> $promsPath{tmp}/exploratory_analysis/$pcaID/error.txt";
				#system "./launchExploratoryAnalyses.pl $clusteringID Cluster 2> $promsPath{tmp}/exploratory_analysis/$clusteringID/error.txt";
				exit;
			}
		}
		if ($clusteringID) {
			my $childPidCluster = fork;
			unless ($childPidCluster) { # child here
				#>Disconnecting from server
				open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
				open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
				open STDERR, '>/dev/null' or die "Can't open /dev/null: $!";
				#system "./launchExploratoryAnalyses.pl $pcaID PCA 2> $promsPath{tmp}/exploratory_analysis/$pcaID/error.txt";
				system "./launchExploratoryAnalyses.pl $clusteringID cluster $projectID $metric $method $itemMetric 2> $promsPath{tmp}/exploratory_analysis/$clusteringID/error.txt";
				exit;
			}
		}

#print "<BR><BR>DONE!"; exit; # DEBUG
		#my $explorID=($pcaID)? $pcaID : $clusteringID;
		print qq
|<SCRIPT LANGUAGE="JavaScript">
top.promsFrame.selectedAction = 'summary';
parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=experiment:$experimentID&branchID=exploranalysis:$explorID&VIEW=explorAna";
</SCRIPT>
</BODY>
</HTML>
|;
	}
	else {# export
		my $errorTxt = '';
		my $exportDIR = "$promsPath{tmp}/export";

if ($exportRprocessMatrices) {

		my $childExport = fork;
		unless ($childExport) { # child here
		    #>Disconnecting from server
		    open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
		    open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
		    open STDERR, '>/dev/null' or die "Can't open /dev/null: $!";
		    system "cd $exportPath; $promsPath{R}/R CMD BATCH --no-save --no-restore '--args $quantifFile $pepFile' $promsPath{R_scripts}/analysisModifQuantif.R";
		    exit;
		}

		##>Waiting for R to finish<##
		my $wait = 1;
		my $count1 = 0;
		my $count2 = 0;
		while ($wait == 1) {
			sleep 2;
			print '.';
			if (!-e "$exportPath/analysisModifQuantif.Rout") {
				if ($count1 > 300) { # <=> 5 min
					$errorTxt='No R output file found';
					$wait=0;
				}
				$count1++;
				next;
			}
			my $Rprocess = `grep -c '> proc.time()' $exportPath/analysisModifQuantif.Rout`;
			chomp $Rprocess;
			if ($Rprocess) { # finished normally
				$wait = 0;
			}
			else {
				$Rprocess = `grep -c '^Execution halted' $exportPath/analysisModifQuantif.Rout`;
				chomp $Rprocess;
				if ($Rprocess) {
					$errorTxt='Execution halted from R';
					$wait = 0;
				}
			}
			$count2++;
			if ($count2 > 900) { # <=> 30 min
				$errorTxt='Process duration has exceeded 30 min.';
				$wait=0;
			}
		}

		if ($errorTxt) {
			print qq
|</FONT>
<BR><BR><BR><BR><BR>
<FONT class="title2" color="#DD0000">***ERROR: Data preprocessing failed ($errorTxt)!***</FONT>
</CENTER>
</BODY>
</HTML>
|;
			exit;
		}
		else {
			unlink "$exportPath/analysisModifQuantif.Rout";
		}

}

		# compress dir and link to donwload
		my $tar=Archive::Tar->new();
		chdir $exportPath;
		$tar->add_files(@file2Tar);
		my $archiveName=$exportDir.'.tgz';
		$tar->write("$promsPath{tmp}/scratch/export/$archiveName",COMPRESS_GZIP) or die "<B>**Error**: Archiving failed:</B> $!<BR>\n";
		sleep 2;
		unlink glob '*.txt';
		chdir $promsPath{cgi};
		rmdir $exportPath;
		
		#rmtree($exportPath); # delete original dir <==== Does not work 04/01/16

		##>Link to download
		print qq
| Done.<BR><BR><BR>
<INPUT type="button" class="title2" value="Download Dataset" onclick="window.location='$promsPath{tmp_html}/scratch/export/$archiveName'"/>
<SCRIPT LANGUAGE="JavaScript">
top.promsFrame.selectedAction = 'summary';
</SCRIPT>
</BODY>
</HTML>
|;
	}
    exit;
}

#### EXPLORATORY ANALYSIS FORM ####
my $sthSelectList = $dbh -> prepare("SELECT CL.ID_CLASSIFICATION,CL.NAME,CA.ID_CATEGORY,CA.NAME,CA.LIST_TYPE FROM CLASSIFICATION CL
					    INNER JOIN CATEGORY CA on CL.ID_CLASSIFICATION = CA.ID_CLASSIFICATION
					    INNER JOIN CATEGORY_PROTEIN CP on CA.ID_CATEGORY = CP.ID_CATEGORY
					    WHERE CL.ID_PROJECT = $projectID ORDER BY CA.DISPLAY_POS");
my $sthGetAllDesign = $dbh -> prepare("SELECT D.ID_DESIGN,D.NAME from EXPERIMENT E INNER JOIN DESIGN D on E.ID_EXPERIMENT = D.ID_EXPERIMENT where E.ID_EXPERIMENT = $experimentID");
my $sthSelModif=$dbh->prepare("SELECT PSI_MS_NAME FROM MODIFICATION where ID_MODIFICATION=?");
my $sthGetExpCondName = $dbh -> prepare("SELECT NAME from EXPCONDITION where ID_EXPCONDITION = ?");
my $sthGetAllInfo = $dbh -> prepare("SELECT Q.ID_QUANTIFICATION, Q.QUANTIF_ANNOT, Q.NAME, Q.ID_MODIFICATION from QUANTIFICATION Q where ID_DESIGN = ?");

#print header(-charset=>'utf-8', -'content-encoding' => 'no'); #debug
my (%listClassCat, %listClassification, %listCategory, %modification);

##Build the restrict/exclude menu
$sthSelectList -> execute;
while (my ($classificationID, $classificationName, $categoryID, $categoryName, $listType) = $sthSelectList -> fetchrow_array) {
	#$listProtein{$proteinID} = $protAlias;
	#print "list:$proteinID<br>";
	next if ($listCategory{$categoryID});
	push @{$listClassCat{$classificationID}}, $categoryID;
	$listClassification{$classificationID} = $classificationName;
	@{$listCategory{$categoryID}} = ($categoryName,$listType);
}
$sthSelectList -> finish;

##Build the Focus menu (protein/phospho and normalized)
$sthGetAllDesign -> execute;
while (my ($designID,$desName) = $sthGetAllDesign -> fetchrow_array) {
    #$designName{$designID}=$desName;
    $sthGetAllInfo -> execute($designID);
    while (my($quantifID, $quantifAnnot, $quantName, $modifID) = $sthGetAllInfo -> fetchrow_array) {
	#$quantifName{$quantifID} = $quantName;
		$modification{$modifID}=1 if $modifID;
    }
}

foreach my $modID (keys %modification) {
	$sthSelModif->execute($modID);
	($modification{$modID})=$sthSelModif->fetchrow_array;
}

$sthSelModif->finish;
$sthGetAllInfo -> finish;
$sthGetAllDesign -> finish;

$dbh -> disconnect;

#exit;
#my ($lightColor,$darkColor)=&promsConfig::getRowColors;
print header(-charset=>'utf-8', -'content-encoding' => 'no');
warningsToBrowser(1);

my ($title,$actionButton)=($action eq 'explorAna')? ("Start Principal Component and/or Clustering Analyses",'Run Analyses') : ('Export Multiple Quantifications From Current Experiment','Export dataset');
print qq
|<HTML>
<HEAD>
<TITLE>Exploratory Analyses</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT language="JavaScript">

function checkUncheckAE (myForm, aeName) {
    var pcaValueBox = myForm.PcaBox.checked;
    var clusterValueBox = myForm.ClusterBox.checked;
    if (aeName == 'Cluster') {
		myForm.ClusterName.disabled = (!clusterValueBox)? true : false;
		myForm.method.disabled = (!clusterValueBox)? true : false;
		myForm.metric.disabled = (!clusterValueBox)? true : false;
    }
    else {
		myForm.PCAName.disabled = (!pcaValueBox)? true : false;
    }
}

var numTestChecked=0;
function updateQuantifStatus(chkStatus,datasetIdx,refID,fromUser) {

	if (chkStatus) {numTestChecked++;} else {numTestChecked--;}
	document.getElementById('chkQuantifs').innerHTML=numTestChecked;
	//Normalization dataset
	if (document.expAnaForm.modif.value < 0) {
		document.getElementById('ref_quantif_'+refID).disabled=(chkStatus)? false : true;
		if (chkStatus == false) {
			document.getElementById('ref_quantif_'+refID).value="";
		}
	}|;
	if ($action eq "export"){
		print qq|
	if (document.getElementById('anovaSel').value == "group"){
//console.log('group_quantif_'+refID);
		document.getElementById('group_quantif_'+refID).disabled=(chkStatus)? false : true;
		if (chkStatus == false){  // && document.expAnaForm.modif.value > 0 )
			document.getElementById('span_quantif_'+datasetIdx).innerHTML= "<input type = text name=group_quantif id=group_quantif_"+refID+" placeholder = " + "'group for Anova'" + " disabled>";
		}
	}|;
	}
	print qq|
	//Auto-extend selection?
	if (fromUser && document.getElementById('autoExtend').checked) {
		var testQuantifs=document.expAnaForm.test_quantif;
		for (var i=datasetIdx+1; i < testQuantifs.length; i++) {
			if (testQuantifs[i].checked == chkStatus) {continue;}
			testQuantifs[i].checked=chkStatus;
			updateQuantifStatus(chkStatus,i,testQuantifs[i].value,false);
		}
	}
}
function copySelectionFromExpAna(expAnaID) {
	if (!expAnaID) return;
	var testQuantifs=document.expAnaForm.test_quantif;
	var refQuantifs=document.expAnaForm.ref_quantif; // undef if no normalization
	/* Deselection */
	if (expAnaID==-1) {
		for (var i=0; i < testQuantifs.length; i++) {
			if (testQuantifs[i].checked == false) {continue;}
			testQuantifs[i].checked=false;
			updateQuantifStatus(false,i,false);
		}
		if (refQuantifs) { // normalization
			for (var i=0; i < refQuantifs.length; i++) {
				refQuantifs[i].selectedIndex=0;
			}
		}
	}
	/* Copy selection from previous Explor Ana */
    //If XHR object already exists, the request is canceled & the object is deleted
    if (XHR && XHR.readyState != 0) {
        XHR.abort();
        delete XHR;
    }
    //Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if ( !XHR ) {
        return false;
    }

    XHR.open("GET","./startExploratoryAnalysis.cgi?AJAX=copySelect&ID=$experimentID&EXPANA_ID="+expAnaID,true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState == 4 && XHR.responseText) {
			var quantifList=XHR.responseText.split('\\n');
			var quantifObject={};
			for (var i=0; i<quantifList.length; i++) {
				var quantifs=quantifList[i].split('%'); // normalization ?
				quantifObject[quantifs[0]]=(quantifs[1])? quantifs[1] : '=';
			}
			for (var i=0; i < testQuantifs.length; i++) {
				var refQ=quantifObject[testQuantifs[i].value];
				if (refQ && !testQuantifs[i].checked) {
					testQuantifs[i].checked=true;
					updateQuantifStatus(true,i,false);
					if (refQ.match('_') && refQuantifs) { // normalization
						for (var j=1; j < refQuantifs[i].options.length; j++) {
							if (refQuantifs[i].options[j].value==refQ) {
								refQuantifs[i].selectedIndex=j;
								break;
							}
						}
					}
				}
			}
        }
    }
    XHR.send(null);
}
function checkForm (myForm) {
|;
if ($action eq 'explorAna') {
	print qq
|   var pcaValueBox = myForm.PcaBox.checked;
    var clusterValueBox = myForm.ClusterBox.checked;

    if (!pcaValueBox && !clusterValueBox) {
		alert('You must select PCA and/or Clustering to continue');
		return false;
    }

    if (!myForm.PCAName.value && pcaValueBox) {
		alert('Provide a name for PCA');
		return false;
    }
    if (!myForm.ClusterName.value && clusterValueBox) {
		alert('Provide a name for Clustering');
		return false;
    }
|;
}
print qq
|    if (!myForm.feature.value) {
		alert('Choose a feature type');
		return false;
    }
	var quantifFam=myForm.quantifFam.value;
	if (!quantifFam) {
		alert('Choose a quantification type');
		return false;
	}
	else if (quantifFam=='RATIO') {
		if (!myForm.foldChange.value \|\| isNaN(myForm.foldChange.value) \|\| myForm.foldChange.value < 1) {
			alert('Choose a valid Absolute fold change (>=1)');
			return false;
		}
		if (!myForm.foldChangeOcc.value \|\| isNaN(myForm.foldChangeOcc.value) \|\| myForm.foldChangeOcc.value < 1) {
			alert('Choose a valid fold change occurence (>=1)');
			return false;
		}
		if (!myForm.pValue.value \|\| isNaN(myForm.pValue.value) \|\| myForm.pValue.value <= 0 \|\| myForm.pValue.value > 1) {
			alert('Choose a valid p-value (]0-1])');
			return false;
		}
		if (!myForm.dataTransform.value.match('LOG') && !confirm('It is highly recommanded to log-transform protein ratios. Proceed anyway?')) {
			return false;
		}
	}
	if (!myForm.numPep.value \|\| isNaN(myForm.numPep.value) \|\| myForm.numPep.value < 1) {
		alert('Choose a valid number of peptides (>=1)');
		return false;
	}
	//Datasets
|;
if ($action eq 'explorAna') {
	print qq
|	if (numTestChecked < 3) {
		alert('At least 3 quantifications are needed to run the analyses');
		return false;
	}
|;
}
else {
	print qq
|	if (numTestChecked < 1) {
		alert('Select at least 1 quantification to be exported');
		return false;
	}
	if (myForm.anovaSel.value == "sample" && numTestChecked == 1) {
		alert('Select at least 2 quantification to be exported');
		return false;
	}
|;
}
print qq
|	var testQuantifs=myForm.test_quantif;
	var refQuantifs=myForm.ref_quantif; // may be null
	for (var i=0; i<testQuantifs.length; i++) {
		if (testQuantifs[i].checked && myForm.modif.value < 0 && !refQuantifs[i].value) { // normalized modif-proteome
			alert('Missing normalization quantification');
			return false;
		}
	}|;
	if ($action eq "export"){
		print qq|
	var objString='';
	var concat=0;
	if (myForm.anovaSel.value == "group") {
	    var groupByQuantifObject={};
	    var groupQuantifs=myForm.group_quantif;
	    for (var i=0; i<groupQuantifs.length;i++) {
			if (!groupQuantifs[i].disabled) {
//console.log(groupQuantifs[i]);
				if (!groupQuantifs[i].value) {alert('Missing value for group');return false;}
				if (!groupByQuantifObject[groupQuantifs[i].value]) {
					groupByQuantifObject[groupQuantifs[i].value]=[];
				}
				groupByQuantifObject[groupQuantifs[i].value].push(testQuantifs[i].value);
			}
	    }

	    for (var group in groupByQuantifObject) {
			if (groupByQuantifObject[group].length<2) {
				alert('2 quantifications are required by group');
				return false;
			}
	    }
	    for (var i in groupByQuantifObject) {
			concat++;
			if (concat == 1) {
				objString = i+'=';
			}
			else {
				objString += '::'+i+'=';
			}

			for (var j in groupByQuantifObject[i]) {
				objString += groupByQuantifObject[i][j]+'\|';
			}
	    }

	    myForm.GROUP_LIST.value=objString;
	}|;
	}
	print qq|
    document.getElementById('expAnaParams').style.display='none';
    return true;
}
//AJAX
var XHR = null;
function getXMLHTTP() {
    var xhr = null;
    if (window.XMLHttpRequest) {// Firefox & others
        xhr = new XMLHttpRequest();
    }
    else if (window.ActiveXObject) { // Internet Explorer
        try {
          xhr = new ActiveXObject("Msxml2.XMLHTTP");
        } catch (e) {
            try {
                xhr = new ActiveXObject("Microsoft.XMLHTTP");
            } catch (e1) {
                xhr = null;
            }
        }
    }
    else { // XMLHttpRequest not supported by browser
        alert("Your browser does not support XMLHTTPRequest objects...");
    }
    return xhr;
}
var selQuantifFam='', selModID=0;
var listDivData='';
function ajaxDisplayData(srcType,srcValue) {

	if (srcType=='focus') { // focus was changed
		var visibRatio;
		selModID=srcValue;
		var QFval=document.expAnaForm.QF.value;
		visibRatio=(QFval=="RATIO")?'':'none';

		//Display exclude amb site if modif
		var visibAmbStatus=(selModID == 0)? 'none' : '';
		var visibAggregStatus=(selModID == 0)? 'none' : '';
		document.getElementById('AMB').style.display=visibAmbStatus;
		document.getElementById('AGGREG').style.display=visibAggregStatus;
		if (selModID == 0) {
			document.expAnaForm.delocalize.checked=false;
			document.expAnaForm.aggregate.checked=false;
		}
|;
		if ($action eq "export") {
			print qq
|			//Back to protein selection with none by defaut (when changing focus)
			if (document.expAnaForm.anovaSel.value != 'none') {
				document.expAnaForm.geneName.checked=false;
				document.expAnaForm.pvalChkAnova.checked=false;
				document.expAnaForm.pValueAnova.disabled=true;
				document.expAnaForm.pValueAnova.value=0.05;
			}
|;
		}
		print qq
|		//document.expAnaForm.anovaSel.value='none';
		document.expAnaForm.aggregate.checked=false;
		document.expAnaForm.topNvarChk.checked=false;
		document.expAnaForm.topNvar.disabled=true;
		document.expAnaForm.topNvar.value=200;
		document.getElementById('ANOVA_PARAM').style.display=visibRatio;
		//document.getElementById('PROT_FILTER').style.display=visibRatio;

		if (selModID==0 ){
			//For phospho, enable just treat as missing value
			var infiniteRatio=document.getElementsByName('infRatios');
			for (var i=0;i<infiniteRatio.length;i++) {
				infiniteRatio[i].disabled = false;
			}
		}
	}
	else if (srcType=='quantifFam') { // quantif type was changed
		selQuantifFam=srcValue;
		var ratioVis,pepTypeOpts=[],dataTransIdx, disabAnova;
		if (!selQuantifFam) { // -= Select =-
			ratioVis='none';
			disabAnova = true;
			pepTypeOpts=[];
			dataTransIdx=0;
		}
		else if (selQuantifFam=='RATIO') {
			ratioVis='';
			disabAnova=false;
			pepTypeOpts=[['NUM_PEP_USED','$pepTypeDesc{NUM_PEP_USED}'],['DIST_PEP_USED','$pepTypeDesc{DIST_PEP_USED}']];
			dataTransIdx=1; // LOG2

			//For phospho, enable just treat as missing value
			var infiniteRatio=document.getElementsByName('infRatios');
			/*if (selModID != 0) {
			    for (var i=0;i<infiniteRatio.length;i++) {
					var disabOpt=(infiniteRatio[i].value == -1 \|\| infiniteRatio[i].value == 0 )? '' : 'disabled';
					infiniteRatio[i].disabled = disabOpt;
			    }
			}*/
		}
		else if (selQuantifFam.match('MQ:')) {
			ratioVis='none';
			pepTypeOpts=[['RAZ_UNI_PEP','$pepTypeDesc{RAZ_UNI_PEP}'],['UNIQUE_PEP','$pepTypeDesc{UNIQUE_PEP}']];
			dataTransIdx=2; // LOG10
			disabAnova = true;
		}
		else {
			ratioVis='none';
			pepTypeOpts=[['IDENT_PEP','$pepTypeDesc{IDENT_PEP}']];
			dataTransIdx=2;
			disabAnova = true;
		}
		document.getElementById('RATIO_FC').style.display=ratioVis;
		document.getElementById('RATIO_PV').style.display=ratioVis;
		document.getElementById('RATIO_INF').style.display=ratioVis;
		//document.getElementById('RATIO_PEP').style.display=ratioVis;
		//document.getElementById('anova_chk').disabled=disabAnova;
		document.getElementById('ANOVA_PARAM').style.display=ratioVis;
		//document.getElementById('PROT_FILTER').style.display=ratioVis;
//console.log('Ratio='+ratioVis);
		document.expAnaForm.dataTransform.options.selectedIndex=dataTransIdx;

		var pepTypeSEL=document.expAnaForm.pepType;
		pepTypeSEL.options.length=0;
		for (var i=0; i<pepTypeOpts.length; i++) {
			pepTypeSEL.options[i]=new Option(pepTypeOpts[i][1],pepTypeOpts[i][0]);
		}
	}
	var listDiv=document.getElementById('listQuantif');
	listDiv.innerHTML = '<B>Fetching data...</B>';
	if (!selQuantifFam) {
		listDiv.innerHTML = '<B>-= Select quantification type =-</B>';
		return;
	}

	//Ajax
    //If XHR object already exists, the request is canceled & the object is deleted
    if (XHR && XHR.readyState != 0) {
        XHR.abort();
        delete XHR;
    }
    //Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if ( !XHR ) {
        return false;
    }

    XHR.open("GET","./startExploratoryAnalysis.cgi?AJAX=displaySelect&MODIF="+selModID+"&QUANTIF_FAM="+selQuantifFam+"&ID=$experimentID",true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState == 4 && XHR.responseText) {
			listDiv.innerHTML=XHR.responseText;
			listDivData=XHR.responseText;
			numTestChecked=0;
        }
    }
    XHR.send(null);
}|;
if ($action eq "export") {
	print qq|
function updateLogTransform(isChk) {
	if (isChk == false) {
		if (document.getElementById('topNvarChk').checked==true){
			return;
		}
	}
	var dataT=document.getElementsByName('dataTransform');
	for (var i=0; i < dataT.length; i++) {
		if (!dataT[i].value) {
			dataT[i].disabled=isChk;
		}
	}
}	
function updateTopNParams(isChk) {
	var topItemIds=['topNvar','anovaSel','pvalChkAnova','pValueAnova'];
	for (let i=0; i<topItemIds.length; i++) {
		document.getElementById(topItemIds[i]).disabled=!isChk;
	}
	if (document.getElementById('aggregate').checked){return;}
	updateLogTransform(isChk);
}
function displayAnova(isChk) {
//console.log(listDivData);
	if (isChk == 'sample') {
		document.getElementById('GROUP_PARAM').style.display='none';
		var testQuantSelect = document.getElementsByName('test_quantif');
		var groupSelect = document.getElementsByName('group_quantif');
		for (var i=0; i < testQuantSelect.length; i++) {
			testQuantSelect[i].disabled=false;
		}
		for (var i=0; i < groupSelect.length; i++) {
			groupSelect[i].value="";
			groupSelect[i].style.visibility='hidden';
		}
		if (document.getElementById('quantifAnnotationType')) {
			document.getElementById('quantifAnnotationType').value="";
		}
	}
	if (isChk=="group") {
		var groupSelect = document.getElementsByName('group_quantif');
		var testQuantSelect = document.getElementsByName('test_quantif');
		for (var i=0; i < testQuantSelect.length; i++) {
			testQuantSelect[i].disabled=true;
		}

		for (var i = 0; i < groupSelect.length; i++) {
			groupSelect[i].style.display='';
		}
		if (document.getElementById('QF').value != "") {
			if (document.getElementById('groupHeaderSpan')){
				document.getElementById('groupHeaderSpan').style.display='';
			}
		}
		document.getElementById('GROUP_PARAM').style.display='';
    }
    else {
		if (document.getElementById('QF').value != "") {
			var refQuantifs=document.expAnaForm.ref_quantif;
			if (document.getElementById('groupHeaderSpan')){
				document.getElementById('groupHeaderSpan').style.display='none';
			}
			var groupSelect = document.getElementsByName('group_quantif');
			var testSelect = document.getElementsByName('test_quantif');
			var refSelect = document.getElementsByName('ref_quantif');
			for (var i = 0; i < groupSelect.length; i++) {
				if (testSelect[i].checked) {
					groupSelect[i].disabled=true;
					groupSelect[i].value='';
					testSelect[i].checked=false;
				}
				groupSelect[i].style.display='none';
				if (refQuantifs){
					if (refSelect[i].disabled == false)  {
						refSelect[i].value='';
						refSelect[i].disabled=true;
					}
				}
			}
			if (document.expAnaForm.allExpAnaSEL){
				document.expAnaForm.allExpAnaSEL.value="";
			}
		}
    }
}

function ajaxPropAnnotate(selProp,propName) {

	var groupQuantifList = document.getElementsByName('group_quantif');
	var testQuantifList = document.getElementsByName('test_quantif');

	if (!selProp) {
		for (var i=0;i<groupQuantifList.length;i++){
			groupQuantifList[i].style.visibility='hidden';
			groupQuantifList[i].value='';
			testQuantifList[i].checked=false;
			testQuantifList[i].disabled=true;
		}
		return;
	}

	if (selProp == 'custom') {
		for (var i=0;i<groupQuantifList.length;i++) {
			groupQuantifList[i].value='';
			groupQuantifList[i].style.visibility='';
			groupQuantifList[i].disabled=true;
			testQuantifList[i].checked=false;
			testQuantifList[i].disabled=false;
		}
		return;
	}
	else {
		for (var i=0;i<groupQuantifList.length;i++) {
			groupQuantifList[i].value='';
		}
	}

	var listQuantif = [];
	for (var i=0;i<testQuantifList.length;i++){
		listQuantif.push(testQuantifList[i].value);
	}

	var quantifList = listQuantif.join(":");
	var paramStrg='AJAX=propAnnotate&ID=$experimentID&qList='+quantifList+'&property='+selProp;

	document.expAnaForm.GROUP_ANNOT.value=propName;

	//If XHR object already exists, the request is canceled & the object is deleted
    if (XHR && XHR.readyState != 0) {
        XHR.abort();
        delete XHR;
    }
    //Creation of the XMLHTTPRequest object
    XHR = getXMLHTTP();
    if (!XHR) {
        return false;
    }
    XHR.open("POST","$promsPath{cgi}/startExploratoryAnalysis.cgi",true); // Switches to synchronous for already saved nnotations
	//Send the proper header information along with the request
    XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
    XHR.setRequestHeader("Content-length", paramStrg.length);
    XHR.setRequestHeader("Connection", "close");
    XHR.onreadystatechange=function() {
        if (XHR.readyState == 4 && XHR.responseText) {
			var resp=XHR.responseText.split('\\n');
			//var quantifProp=[];
			for (var i=0; i < resp.length-1; i++) {
				var list=resp[i].split(',');
				for (var j=1; j<list.length;j++) {
//console.log('j='+list[j]);
					//quantifProp.push('group_quantif_'+list[j]);
					document.getElementById('group_quantif_'+list[j]).value=list[0];
					document.getElementById('group_quantif_'+list[j]).style.visibility='';
					document.getElementById('group_quantif_'+list[j]).disabled=true;
					/*document.getElementById('test_quantif_'+list[j]).disabled=false;
					document.getElementById('test_quantif_'+list[j]).checked=false;*/
				}
			}
			for (var i=0;i<testQuantifList.length;i++){
				if (testQuantifList[i].disabled) {
					testQuantifList[i].disabled=false;
					groupQuantifList[i].style.visibility='';
					groupQuantifList[i].disabled=true;
				}
			}

        }
    }
    XHR.send(paramStrg);
}
	|;
}
else {
	print qq|
	function updateTopNParams(isChk) {
		
		document.getElementById('topNvar').disabled=!isChk;
	}
	|;
}
print qq|

function changeFeature (val) {
	if (val == "protQuant") {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/startExploratoryAnalysis.cgi?ID=$experimentID";
	}
	else {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/startExploratoryAnalysisPeptide.cgi?ID=$experimentID";
	}
}

</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">

<CENTER>
<FONT class="title">$title</FONT><BR><BR>
<DIV id="expAnaParams">
<FORM name="expAnaForm" method="POST" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="ID" value="$experimentID">
<INPUT type="hidden" name="ACT" value="$action">|;
if ($action eq "export") {
print qq|
<INPUT type="hidden" name="GROUP_LIST" value="">
<INPUT type="hidden" name="GROUP_ANNOT" value="">|;
}
print qq|

<TABLE border="0" bgcolor="$darkColor" width="1000px">
|;
if ($action eq 'explorAna') {
	print qq
|    <TR>
		<TD><INPUT type="checkbox" name="PcaBox" onclick="checkUncheckAE(document.expAnaForm, 'PCA')" checked></TD>
		<TH align=right nowrap>PCA name :</TH><TD bgcolor="$lightColor"><INPUT type="text" name="PCAName" size="35" maxlength="50" value=""/>&nbsp;</TD>
    </TR>
    <TR>
		<TD><INPUT type="checkbox" name="ClusterBox" onclick="javascript:checkUncheckAE(document.expAnaForm, 'Cluster')" checked></TD>
		<TH align=right nowrap>Clustering name :</TH>
		<TD bgcolor="$lightColor"><INPUT type="text" name="ClusterName" size="35" maxlength="50" value=""/>&nbsp;</TD>
    </TR>
|;
}
print qq
|	<TR><TH colspan="2" align="right" nowrap>Focus :</TH><TD bgcolor="$lightColor"><SELECT name="modif" id="modif" onchange="ajaxDisplayData('focus',this.value)">
	<OPTION value="0">Proteins</OPTION>
	<OPTGROUP label="Modifications:">
|;
my $disabFeature=($action eq "export")? " disabled" : "";
if (scalar keys %modification) {
	foreach my $modID (sort{lc($modification{$a}) cmp lc($modification{$b})} keys %modification) {
		print "<OPTION value=\"$modID\">$modification{$modID}-proteome</OPTION>\n";
	}
	print "\t</OPTGROUP>\n<OPTGROUP label=\"Normalized modifications:\">\n";
	foreach my $modID (sort{lc($modification{$a}) cmp lc($modification{$b})} keys %modification) {
		print "<OPTION value=\"-$modID\">Normalized $modification{$modID}-proteins</OPTION>\n";
	}
}
else {print "<OPTION value=\"\" disabled>*None found*</OPTION>\n";}
print qq
|	</OPTGROUP>
	</SELECT></TD>
	</TR>
    <TR>
		<TH colspan="2" align=right nowrap>Features :</TH>
		<TD nowrap bgcolor=$lightColor>
			<SELECT name="feature" onchange="javascript:changeFeature(this.value)">
			<OPTION value="">-= Select =-</OPTION>
			<OPTION value="protQuant" selected>Protein quantifications</OPTION>
			<OPTION value="pepCount"$disabFeature> Peptide counts</OPTION>
			</SELECT>
			<SPAN>&nbsp;&nbsp;<B>Quantification type:<SELECT name="quantifFam" id="QF" onchange="ajaxDisplayData('quantifFam',this.value)">
			<OPTION value="">-= Select =-</OPTION>
|;
foreach my $qFamily (sort{lc($proteinQuantifFamilies{'NAME'}{$a}) cmp lc($proteinQuantifFamilies{'NAME'}{$b})} keys %{$proteinQuantifFamilies{'NAME'}}) {
	my $isParent=0;
	foreach my $subFamily (@{$proteinQuantifFamilies{'MEMBERS'}{$qFamily}}) {
		if ($proteinQuantifFamilies{'MEASURES'}{$subFamily}) {
			$isParent=1;
			last;
		}
	}
	if ($isParent) {
		print "<OPTGROUP label=\"$proteinQuantifFamilies{NAME}{$qFamily}:\">";
		foreach my $subFamily (@{$proteinQuantifFamilies{'MEMBERS'}{$qFamily}}) {
			foreach my $refMeas (@{$proteinQuantifFamilies{'MEASURES'}{$subFamily}}) {print "<OPTION value=\"$subFamily:$refMeas->[0]\">$refMeas->[1]</OPTION>";}
		}
		print "</OPTGROUP>\n";
	}
	else {print "<OPTION value=\"$qFamily\">$proteinQuantifFamilies{NAME}{$qFamily}</OPTION>\n";}
	#if ($qFamily eq 'EMPAI') {
	#	print "<OPTION value=\"${qFamily}_MOL\">$proteinQuantifFamilies{NAME}{$qFamily} (Molar fraction)</OPTION>\n<OPTION value=\"${qFamily}_MR\">$proteinQuantifFamilies{NAME}{$qFamily} (Weight fraction)</OPTION>\n";
	#}

}
my $selAllInf='';
my $numSelPep=3;
my $sel34MissVal=' selected';
my $selAllMissVal='';
my $selGroupAnova='';
my $strgKeep='';
if ($action eq 'export') {
	$selAllInf=' selected';
	$numSelPep=1;
	$sel34MissVal='';
	$selAllMissVal=' selected';
	$strgKeep="onchange=\"updateLogTransform(this.checked)\"";
	#$strgAggregate="onchange=\"document.getElementById('geneName').disabled=!this.checked\"";
}
else {
	$selGroupAnova= ' disabled';
}
print qq
|			</SELECT></SPAN>
		</TD>
    </TR>
	<TR>
	<TH colspan="2" align=right valign=top nowrap>Data transform :</TH>
	<TD nowrap bgcolor=$lightColor>&bull;<SELECT name="dataTransform"><OPTION value="">No change to</OPTION><OPTION value="LOG2">Log2-transform</OPTION><OPTION value="LOG10">Log10-transform</OPTION></SELECT> <B>values</B>
|;
print qq|&nbsp;&nbsp;&nbsp;&nbsp;&bull;<INPUT type="checkbox" id="geneName" name="geneName" value="TRUE"><B>Use gene name as identifier</B>| if ($action eq "export");
print qq
|	</TD>
	</TR>
    <TR>
	<TH colspan="2" align=right valign=top nowrap>Data filtering :</TH>
	<TD nowrap bgcolor=$lightColor>
	    <TABLE cellspacing=0>
		<TR id="RATIO_FC" style="display:none">
		    <TH align="right" nowrap>&nbsp;Abs. fold change &ge;</TH>
		    <TD><INPUT type="text" name="foldChange" value="1" size="2">&nbsp;<B>in at least&nbsp;</B><INPUT type="number" name="foldChangeOcc" value="1"  style="width:50px">&nbsp;<B>quantification</B>&nbsp;</TD>
		</TR>
		<TR id="RATIO_INF" style="display:none">
		    <TH align="right">Infinite ratios:</TH>
			<TD><SELECT id="infRatios" name="infRatios"><OPTION value="-1">Treat as missing values</OPTION><OPTION value="0">None allowed</OPTION><OPTION value="5">Allow 5%</OPTION><OPTION value="10">Allow 10%</OPTION><OPTION value="25">Allow 25%</OPTION><OPTION value="34">Allow 34%</OPTION><OPTION value="50">Allow 50%</OPTION><OPTION value="100"$selAllInf>Allow all</OPTION></SELECT>&nbsp;per protein</TD>
			<!--<TD><SELECT id="infRatios" name="infRatios"><OPTION  name="infRatios" value="-1">Treat as missing values</OPTION><OPTION name="infRatios" value="0">None allowed</OPTION><OPTION name="infRatios" value="5">Allow 5%</OPTION><OPTION name="infRatios" value="10">Allow 10%</OPTION><OPTION name="infRatios" value="25">Allow 25%</OPTION><OPTION name="infRatios" value="34">Allow 34%</OPTION><OPTION name="infRatios" value="50">Allow 50%</OPTION><OPTION name="infRatios" value="100">Allow all</OPTION></SELECT>&nbsp;per protein</TD>-->

		</TR>
		<TR id="RATIO_PV" style="display:none">
		    <TH align="right">p-value &le;</TH>
		    <TD nowrap><INPUT type="text" name="pValue" value="1" size="5">&nbsp;<B>in at least&nbsp;</B><INPUT type="number" name="pValueOcc" value="1"  style="width:50px">&nbsp;<B>quantification</B>&nbsp;<SMALL>(Does not apply to normalized ratios)</SMALL></TD>
		</TR>
		<TR><!-- id="RATIO_PEP" style="display:none" -->
		    <TH align="right">Peptides:</TH>
		    <TD><SELECT name="pepType"><!-- *Updated by JavaScript* --></SELECT> &ge; <INPUT type="number" name="numPep" value="$numSelPep" style="width:50px"></TD>
		</TR>
		<TR>
		    <TH align="right">Missing values:</TH>
			<TD><SELECT name="missingValues"><OPTION value="0">None allowed</OPTION><OPTION value="5">Allow 5%</OPTION><OPTION value="10">Allow 10%</OPTION><OPTION value="25">Allow 25%</OPTION><OPTION value="34"$sel34MissVal>Allow 34%</OPTION><OPTION value="50">Allow 50%</OPTION><OPTION value="100"$selAllMissVal>Allow all</OPTION></SELECT>&nbsp;per protein</TD>
		</TR>
		<TR id="AMB" style="display:none">
		 <TH align="left" colspan=2><INPUT type="checkbox" id="delocalize" name="delocalize" >Exclude ambiguous sites</TH>
		<TR>
		<TR id="AGGREG" style="display:none">
		  <TH align="left" colspan=2><INPUT type="checkbox" id="aggregate" name="aggregate" $strgKeep>Keep only the most "informative" site per protein</TH>
		  <TD></TD>
		</TR>
		</TABLE>
	</TD>
	</TR>
	<TR id="PROT_FILTER">
		<TH colspan="2" align=right valign=top nowrap>Protein filtering :</TH>
		<TD nowrap bgcolor=$lightColor>
			<TABLE cellspacing=0>
				<TR>
					<TH align="right">Custom list:</TH>
					<TD nowrap>
						<SELECT name="condRestrict">
							<OPTION value="restrict">Restrict to </OPTION>
							<OPTION value="exclude">Exclude</OPTION>
						</SELECT>
						<SELECT name="list">
							<OPTION value="0">All proteins</OPTION>
							|;
							foreach my $classID (sort{lc($listClassification{$a}) cmp lc($listClassification{$b})} keys %listClassification) {
								print "<OPTGROUP LABEL=\"$listClassification{$classID}:\">\n";
								foreach my $catID (@{$listClassCat{$classID}}) { # ordered by displayPos from SQL query
									my $typeStrg=($listCategory{$catID}[1] eq 'SITE')? ' [Sites]' : '';
									print "<OPTION value=\"$catID\">$listCategory{$catID}[0]$typeStrg</OPTION>\n";
								}
								print "</OPTGROUP>\n";
							}
							print qq|
						</SELECT>
					</TD>
				</TR>
				<TR id="ANOVA_PARAM" style="display:none">
					<TH align="left">
						<INPUT type="checkbox" name="topNvarChk" id="topNvarChk" value="1" onchange="updateTopNParams(this.checked)">Keep only the:
					</TH>
					<TD nowrap>
						<INPUT type="number" name="topNvar" id="topNvar" value="200"  style="width:80px" disabled> <B>most changing proteins</B>|;
						if ($action eq "export") {
							print qq
							|&nbsp;<B>between:</B>
							<SELECT name="anovaSel" id="anovaSel" onchange="displayAnova(this.value)" disabled>
								<OPTION value="sample">individual quantifications (Std dev.)</OPTION>
								<OPTION value="group" $selGroupAnova>sets of quantifications (Anova)</OPTION>
							</SELECT>
							<!--<SPAN id="GROUP_PARAM" style="display:none">
								<INPUT type="checkbox" id="pvalChkAnova" name="pvalChkAnova" value="1" onchange="document.getElementById('pValueAnova').disabled=!this.checked"><B>p-value &le; </B><INPUT type="text" name="pValueAnova" id="pValueAnova" value="0.05" size="5" disabled>
							</SPAN>-->|;
						}
					print qq|
				   </TD>
				</TR>
				<TR id="GROUP_PARAM" style="display:none">
					<TH></TH>
					<TD><INPUT type="checkbox" id="pvalChkAnova" name="pvalChkAnova" value="1" onchange="document.getElementById('pValueAnova').disabled=!this.checked"><B>Further restrict to Anova p-value &le;</B>
						<INPUT type="text" name="pValueAnova" id="pValueAnova" value="0.05" size="5" disabled>
					</TD>
				</TR>
			</TABLE>
		</TD>
    </TR>
|;
if ($action eq 'explorAna') {
	print qq
|    <TR>
	<TH colspan="2" align=right nowrap>&nbsp;Clustering parameters :</TH>
	<TD nowrap bgcolor="$lightColor">
	    <B>Method:</B><SELECT name="method">
			<OPTION value="ward.D" selected>ward</OPTION>
			<OPTION value="average">average</OPTION>
			<OPTION value="single">single</OPTION>
			<OPTION value="complete">complete</OPTION>
	    </SELECT>&nbsp;
	    <B>Metric:</B><SELECT name="metric">
		    <OPTGROUP label="Correlation">
			<OPTION value="cor:pearson" selected>pearson</OPTION>
			<OPTION value="cor:pearsonabs">pearsonabs</OPTION>
			<OPTION value="cor:spearman">spearman</OPTION>
			<OPTION value="cor:spearmanabs">spearmanabs</OPTION>
		    </OPTGROUP>
		    <OPTGROUP label="Distance">
			<OPTION value="dist:euclidean">euclidean</OPTION>
			<OPTION value="dist:manhattan">manhattan</OPTION>
		    </OPTGROUP>
	    </SELECT>
	</TD>
    </TR>
|;
}
print qq
|	<TR><TH colspan=2 align="right" valign="top">Data selection :</TH><TD nowrap bgcolor="$lightColor"><DIV id="listQuantif" style="max-height:300px;overflow:auto"><B>-= Select quantification type =-</B></DIV></TD></TR>
    <TR>
	<TD colspan="3" align="center">
	    <INPUT type="submit" name="submitted" value="$actionButton">
	</TD>
    </TR>
</TABLE>
</FORM>
</DIV>
<BR><BR>
</BODY>
</HTML>
|;

sub ajaxDisplayData {

	my ($modificationID,$quantifFam,$fromMotif) = @_;

	my $normDataset=($modificationID < 0)? 1 : 0;
	my $phosphoMod = ($modificationID > 0)?1 : 0;
	$modificationID=abs($modificationID);

	my (%testQuantifications,%normQuantifications,%designName,%quantifName,%fetchedCond,%experimentInfo,%noDesignQuantifications,%sampleInfo,%analysisInfo);

	####<RATIO or MaxQuant Intensities>####
	if ($quantifFam =~ /RATIO|MQ/) {
		my $quantifFam0=$quantifFam;
		($quantifFam,my $quantifCode)=split(':',$quantifFam); # MQ:MQ_INT -> MQ
		my @quantifMethIDs;
		foreach my $subFamily (@{$proteinQuantifFamilies{'MEMBERS'}{$quantifFam}}) {
#print "Sub:$subFamily<br>";
			my ($quantifMethID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='$subFamily'");
			push @quantifMethIDs,$quantifMethID;
		}
		my $modStrg=($modificationID)? "AND ID_MODIFICATION=$modificationID" : 'AND ID_MODIFICATION IS NULL';

		###<Design quantifications>###
		my $sthGetAllDesign = $dbh -> prepare("SELECT D.ID_DESIGN,D.NAME FROM EXPERIMENT E INNER JOIN DESIGN D ON E.ID_EXPERIMENT = D.ID_EXPERIMENT WHERE E.ID_EXPERIMENT = ?");
		my $sthGetExpCondName = $dbh -> prepare("SELECT NAME FROM EXPCONDITION where ID_EXPCONDITION = ?");
		my $sthGetTestInfo = $dbh -> prepare("SELECT ID_QUANTIFICATION,QUANTIF_ANNOT,NAME,ID_MODIFICATION FROM QUANTIFICATION WHERE ID_DESIGN=? AND ID_QUANTIFICATION_METHOD IN (".join(',',@quantifMethIDs).") $modStrg");
		my $sthGetMeas;
		if ($quantifCode) {
			my ($quantifParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD IN (".join(',',@quantifMethIDs).") AND CODE='$quantifCode'");
			$sthGetMeas=$dbh->prepare("SELECT 1 FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIF_PARAMETER=$quantifParamID AND ID_QUANTIFICATION=? LIMIT 0,1");
		}

		&fetchDesignQuantifications($experimentID,$quantifFam,\%testQuantifications,\%designName,\%quantifName,$sthGetAllDesign,$sthGetExpCondName,$sthGetTestInfo,$sthGetMeas);

		$sthGetTestInfo -> finish;
		$sthGetMeas -> finish if $sthGetMeas;

		if ($normDataset) {
			my $sthExp=$dbh->prepare("SELECT ID_EXPERIMENT,NAME,DISPLAY_POS FROM EXPERIMENT WHERE ID_PROJECT=$projectID"); # $projectID is global
			my $sthGetNormInfo = $dbh->prepare("SELECT ID_QUANTIFICATION,QUANTIF_ANNOT,NAME,ID_MODIFICATION FROM QUANTIFICATION WHERE ID_DESIGN = ? AND ID_QUANTIFICATION_METHOD IN (".join(',',@quantifMethIDs).") AND ID_MODIFICATION IS NULL");
			$sthExp->execute;
			while (my ($expID,$expName,$dispPos)=$sthExp->fetchrow_array) {
				@{$experimentInfo{$expID}}=($expName,$dispPos);
				%{$normQuantifications{$expID}}=();
				&fetchDesignQuantifications($expID,$quantifFam,$normQuantifications{$expID},\%designName,\%quantifName,$sthGetAllDesign,$sthGetExpCondName,$sthGetNormInfo);
			}
			$sthExp->finish;
			$sthGetNormInfo->finish;
		}

		$sthGetAllDesign -> finish;
		$sthGetExpCondName->finish;

		###<Non-design quantifications>###
		#unless ($normDataset) { # skip non-design quantif if normalization
		if (!$normDataset && $quantifFam ne "MQ") {
			my ($quantifMethID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='PROT_RATIO_PEP'");
			my $sthND=$dbh->prepare("SELECT S.ID_SAMPLE,S.NAME,S.DISPLAY_POS,A.ID_ANALYSIS,A.NAME,A.DISPLAY_POS,Q.ID_QUANTIFICATION,Q.NAME,QUANTIF_ANNOT FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ,ANALYSIS A,SAMPLE S WHERE Q.ID_QUANTIFICATION=AQ.ID_QUANTIFICATION AND AQ.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND Q.ID_QUANTIFICATION_METHOD=$quantifMethID AND S.ID_EXPERIMENT=$experimentID AND FOCUS='protein' AND ID_DESIGN IS NULL $modStrg");
			$sthND->execute;
			while (my ($sampID,$sampName,$sampPos,$anaID,$anaName,$anaPos,$quantifID,$quantName,$quantifAnnot)=$sthND->fetchrow_array) {
				#push @{$noDesignQuantifications{$sampID}{$anaID}},$quantifID;
				@{$analysisInfo{$anaID}}=($anaName,$anaPos);
				@{$sampleInfo{$sampID}}=($sampName,$sampPos);
				$quantifName{$quantifID}=$quantName || 'No name';
				my (%labelingInfo,%stateInfo);
				&promsQuantif::extractQuantificationParameters($dbh,$quantifAnnot,\%labelingInfo,\%stateInfo);
				my $ratioPos = 0;
				foreach my $ratio (@{$labelingInfo{'RATIOS'}}) {
					$ratioPos++;
					my ($testStatePos,$refStatePos) = split(/\//,$ratio);
					push @{$noDesignQuantifications{$sampID}{$anaID}{$quantifID}},["$stateInfo{$testStatePos}{NAME}/$stateInfo{$refStatePos}{NAME}",$ratioPos];
				}
			}
			$sthND->finish;
		}
	}

	####<emPAI,SIN>####
	elsif ($quantifFam=~/^EMPAI|SIN/) {
		#my $quantifCode=$quantifFam;
		#$quantifFam=~s/_.+//;
		($quantifFam,my $quantifCode)=split(':',$quantifFam);
		my ($quantifMethID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='$quantifFam'");
		my $empaiQuantifParamID;
		if ($quantifCode=~/_/) { # EMPAI_MOL or EMPAI_MR
			($empaiQuantifParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=$quantifMethID AND CODE='$quantifCode'");
		}
		my $sthQ=$dbh->prepare("SELECT S.ID_SAMPLE,S.NAME,S.DISPLAY_POS,A.ID_ANALYSIS,A.NAME,A.DISPLAY_POS,Q.ID_QUANTIFICATION,Q.NAME FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ,ANALYSIS A,SAMPLE S WHERE Q.ID_QUANTIFICATION=AQ.ID_QUANTIFICATION AND AQ.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND Q.ID_QUANTIFICATION_METHOD=$quantifMethID AND S.ID_EXPERIMENT=$experimentID AND FOCUS='protein'");
		my $sthQ2=($empaiQuantifParamID)? $dbh->prepare("SELECT 1 FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=? AND ID_QUANTIF_PARAMETER=$empaiQuantifParamID LIMIT 0,1") : undef;
		$sthQ->execute;
		while (my ($sampID,$sampName,$sampPos,$anaID,$anaName,$anaPos,$quantifID,$quantName)=$sthQ->fetchrow_array) {
			if ($empaiQuantifParamID) {
				$sthQ2->execute($quantifID);
				my ($ok)=$sthQ2->fetchrow_array;
				next unless $ok;
			}
			push @{$noDesignQuantifications{$sampID}{$anaID}{$quantifID}},['',0];
			@{$analysisInfo{$anaID}}=($anaName,$anaPos);
			@{$sampleInfo{$sampID}}=($sampName,$sampPos);
			$quantifName{$quantifID}=$quantName || 'No name';
		}
		$sthQ->finish;
		$sthQ2->finish if $sthQ2;
	}

	####<Fetch ExplorAna with same quantifFam>####
	my %allExplorAna;
	my $sthAllExp=$dbh->prepare("SELECT ID_EXPLORANALYSIS,NAME,ANA_TYPE FROM EXPLORANALYSIS WHERE ID_EXPERIMENT=$experimentID AND PARAM_LIST LIKE '%quantFam=$quantifFam%'");
	$sthAllExp->execute;
	while (my ($eaID,$eaName,$anaType)=$sthAllExp->fetchrow_array) {
		$allExplorAna{$anaType}{$eaID}=$eaName;
	}
	$sthAllExp->finish;

	my (%quantificationIDs,%bioSampProperties);
	if (scalar %testQuantifications) {
		foreach my $designID (keys %testQuantifications) {
			foreach my $quantifID (keys %{$testQuantifications{$designID}}) {
				$quantificationIDs{$quantifID}=1;
			}
		}
		my $sthProp=$dbh->prepare("SELECT DISTINCT P.ID_PROPERTY,P.NAME,P.PROPERTY_TYPE FROM PROPERTY P
								INNER JOIN BIOSAMPLE_PROPERTY BP ON P.ID_PROPERTY=BP.ID_PROPERTY
								INNER JOIN OBSERVATION O ON BP.ID_BIOSAMPLE=O.ID_BIOSAMPLE
								INNER JOIN ANA_QUANTIFICATION AQ ON O.ID_ANALYSIS=AQ.ID_ANALYSIS
								WHERE P.USE_IN_ANALYSIS=1 AND AQ.ID_QUANTIFICATION IN (".join(',',keys %quantificationIDs).")");

		$sthProp->execute;
		while (my($propID,$propName,$propType)=$sthProp->fetchrow_array) {
			$bioSampProperties{$propType}{$propID}=$propName;
		}
		$sthProp->finish;
	}
	$dbh->disconnect;

    print header(-'content-encoding' => 'no',-charset => 'utf-8');
    warningsToBrowser(1);

	if (!scalar keys %testQuantifications && !scalar keys %noDesignQuantifications){
		print "<FONT class=\"title3\">No quantifications found.</FONT>\n";
		exit;
	}


	if (!$fromMotif) {
			print qq|
		<INPUT type="checkbox" id="autoExtend"><B>Auto-extend selection</B>
		<TABLE cellspacing=0 width=100% >
		<TR bgcolor="$darkColor"><TH>&nbsp;Quantifications (<SPAN id="chkQuantifs">0</SPAN> selected)&nbsp;<BR>
		&nbsp;Copy selection from:<SELECT name="allExpAnaSEL" id="allExpAnaSEL" onchange="copySelectionFromExpAna(this.value)"><OPTION value="">-= Select =-</OPTION><OPTION value="-1">** Deselect all **</OPTION>
		|;
			foreach my $anaType (sort{lc($a) cmp lc($b)} keys %allExplorAna) {
				my $typeLabel=($anaType=~/^clus/i)? 'Clustering' : 'PCA';
				print "<OPTGROUP label=\"$typeLabel:\">\n";
				foreach my $eaID (sort{lc($allExplorAna{$anaType}{$a}) cmp lc($allExplorAna{$anaType}{$b})} keys %{$allExplorAna{$anaType}}) {
					print "<OPTION value=\"$eaID\">$allExplorAna{$anaType}{$eaID}</OPTION>\n";
				}
				print "</OPTGROUP>\n";
			}
			print qq
		|</SELECT>&nbsp;
		</TH>
		<TH><SPAN id="groupHeaderSpan" style="display:none">
				<TABLE cellspacing=0 cellpadding=0 border=0>
				<TR><TD bgcolor="$darkColor" class="title3" style="height:30px" align="center">&nbsp;Sets based on:</TD>
				<TD valign="top">
				<SELECT name="quantifAnnotationType" id="quantifAnnotationType" class="annotation title3" onchange="ajaxPropAnnotate(this.value,this.options[this.selectedIndex].text,false)">
				<OPTION value="">-= Select =-</OPTION>
				<OPTGROUP label="Custom:">
					<OPTION value="custom">Manual selection</OPTION>
				</OPTGROUP>
				<OPTGROUP label="BioSample properties:">|;
				if (scalar keys %{$bioSampProperties{'O'}}) {
					foreach my $propID (sort{lc($bioSampProperties{'O'}{$a}) cmp lc($bioSampProperties{'O'}{$b})} keys %{$bioSampProperties{'O'}}) {
						print "\t<OPTION value=\"property:O:$propID\">$bioSampProperties{O}{$propID}</OPTION>\n";
					}
				}
				else {print "\t<OPTION value=\"\" disabled>None</OPTION>\n";}
			print qq
			|	</OPTGROUP>
				<OPTGROUP label="BioSample treatments:">
				|;
				if (scalar keys %{$bioSampProperties{'T'}}) {
					foreach my $propID (sort{lc($bioSampProperties{'T'}{$a}) cmp lc($bioSampProperties{'T'}{$b})} keys %{$bioSampProperties{'T'}}) {
						print "\t<OPTION value=\"property:T:$propID\">$bioSampProperties{T}{$propID}</OPTION>\n";
					}
				}
				else {print "\t<OPTION value=\"\" disabled>None</OPTION>\n";}
			print qq
			|	</OPTGROUP>
				</SELECT></TD></TR></TABLE>
								</SPAN></TH>|;
								print "<TH>&nbsp;Reference&nbsp;</TH>\n" if $normDataset;
								print "</TR>";
	}

	my $colSpan=3;#($normDataset && $anovaTag eq "group")? 3 : 2;
	my $bgColor=$lightColor;
	my $datasetIdx=0;
	#my $refQuantif=($modificationID)? \%testQuantifications : \%normQuantifications;
	my $ratioNameSep=($quantifFam=~/RATIO|MQ/)? ' : ' : '';

	###<Design quantifs>###
	if (scalar keys %testQuantifications) {
		if (!$fromMotif) {
			print qq
			|<TR><TH colspan="$colSpan" align="left"></TH></TR>
			<TR bgcolor="$darkColor"><TH colspan="$colSpan" align="left">Design-based quantifications:</TH></TR>
			|;
		}
		else {
			print qq|<SELECT name="quantif" id="quantif"><OPTION value="">--=Select=--</OPTION>|;
		}

		my $prevDesignID=0;
		foreach my $designID (sort{&promsMod::sortSmart(lc($designName{$a}),lc($designName{$b}))} keys %testQuantifications) {
			if ($fromMotif) {
				print qq|<OPTGROUP label="$designName{$designID} :">|;
			}

			my $prevQuantifID=0;
			foreach my $quantifID (sort{&promsMod::sortSmart(lc($quantifName{$a}),lc($quantifName{$b}))} keys %{$testQuantifications{$designID}}) {
				foreach my $refRatio (@{$testQuantifications{$designID}{$quantifID}}) {
					my $designNameStrg=($designID==$prevDesignID)? '<SPAN style="visibility:hidden">'.$designName{$designID}.' > </SPAN>' : $designName{$designID}.' > ';
					my $quantifNameStrg=($quantifID==$prevQuantifID)? '<SPAN style="visibility:hidden">'.$quantifName{$quantifID}.$ratioNameSep.'</SPAN>' : $quantifName{$quantifID}.$ratioNameSep;
					my ($ratioName,$rank)=@{$refRatio};
					my $valueStrg=$quantifID.'_'.$rank;

					if ($fromMotif) {
						print qq|<OPTION value="$valueStrg">&nbsp;&nbsp;$quantifNameStrg$ratioName</OPTION>	|;
					}
					else {
						print qq
						|<TR class="list" bgcolor="$bgColor">|;
						if (!$fromMotif) {
							print qq|<TD nowrap><INPUT type="checkbox" name="test_quantif" value="$valueStrg" onclick="updateQuantifStatus(this.checked,$datasetIdx,'$valueStrg',true)">$designNameStrg$quantifNameStrg$ratioName</TD>|;
						}
						else{
							print qq|<TD nowrap><INPUT type="checkbox" name="test_quantif" id="test_quantif_$valueStrg" value="$valueStrg" onclick="updateQuantifStatus(this.checked,'$valueStrg',true)">$designNameStrg$quantifNameStrg$ratioName</TD>|;
						}
					}

					#if ($phosphoMod && !$fromMotif) {
					#    print "<TD><INPUT TYPE=\"text\" name=\"group_quantif\" id=\"group_quantif_$datasetIdx\" placeholder=\"group for Anova\" disabled style=\"display:none\"></TD>";
					#}
					if (!$fromMotif ) {
					    print "<TD><INPUT TYPE=\"text\" name=\"group_quantif\" id=\"group_quantif_$valueStrg\" style=\"visibility:hidden\"></TD>";
					}
					if ($normDataset) {
						#print "<TD><SELECT name=\"ref_quantif\" id=\"ref_quantif_$datasetIdx\" disabled><OPTION value=\"\">-= Select =-</OPTION>\n";
						print "<TD><SELECT name=\"ref_quantif\" id=\"ref_quantif_$valueStrg\" disabled><OPTION value=\"\">-= Select =-</OPTION>\n";
						foreach my $expID (sort{$experimentInfo{$a}[1] <=> $experimentInfo{$b}[1]} keys %normQuantifications) {
							next unless scalar keys %{$normQuantifications{$expID}}; # no quantif
							print "<OPTGROUP label=\"Experiment $experimentInfo{$expID}[0]:\">\n";
							foreach my $desID (sort{&promsMod::sortSmart(lc($designName{$a}),lc($designName{$b}))} keys %{$normQuantifications{$expID}}) {
								print "<OPTGROUP label=\"&nbsp;&nbsp;Design $designName{$desID}:\">\n";
								foreach my $quanID (sort{&promsMod::sortSmart(lc($quantifName{$a}),lc($quantifName{$b}))} keys %{$normQuantifications{$expID}{$desID}}) {
									print qq"<OPTGROUP label=\"&nbsp;&nbsp;&nbsp;&nbsp;$quantifName{$quanID}\">\n";
									foreach my $refRatioRef (@{$normQuantifications{$expID}{$desID}{$quanID}}) {
										my ($ratioNameRef,$rankRef)=@{$refRatioRef};
										my $valueRefStrg=$quanID."_".$rankRef;
										print qq"<OPTION value=\"$valueRefStrg\">&nbsp;&nbsp;$ratioNameRef</OPTION>\n";
									}
									print "</OPTGROUP>\n";
								}
								print "</OPTGROUP>\n";
							}
							print "</OPTGROUP>\n";
						}
						print "</SELECT></TD>";
					}
					if (!$fromMotif) {print "</TR>\n";}
					$datasetIdx++;
					$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
					$prevDesignID=$designID;
					$prevQuantifID=$quantifID;
				}
			}
		}
		if ($fromMotif) {
				print qq|</SELECT>|;
			}
	}

	###<No Design quantifs>###
	if (scalar keys %noDesignQuantifications) {
		print qq
		|<TR><TH colspan="$colSpan" align="left"></TH></TR>
		<TR bgcolor="$darkColor"><TH align="left" colspan="$colSpan">Internal-analysis quantifications:</TH></TR>
		|;
		$bgColor=$lightColor;
		my $prevSampleID=0;
		foreach my $sampleID (sort{$sampleInfo{$a}[1]<=>$sampleInfo{$b}[1]} keys %noDesignQuantifications) {
			my $prevAnaID=0;
			foreach my $anaID (sort{$analysisInfo{$a}[1]<=>$analysisInfo{$b}[1]} keys %{$noDesignQuantifications{$sampleID}}) {
				my $prevQuantifID=0;
				foreach my $quantifID (sort{&promsMod::sortSmart(lc($quantifName{$a}),lc($quantifName{$b}))} keys %{$noDesignQuantifications{$sampleID}{$anaID}}) {
					foreach my $refRatio (@{$noDesignQuantifications{$sampleID}{$anaID}{$quantifID}}) {

						my $sampleNameStrg=($sampleID==$prevSampleID)? '<SPAN style="visibility:hidden">'.$sampleInfo{$sampleID}[0].' > </SPAN>' : $sampleInfo{$sampleID}[0].' > ';
						my $analysisNameStrg=($anaID==$prevAnaID)? '<SPAN style="visibility:hidden">'.$analysisInfo{$anaID}[0].' > </SPAN>' : $analysisInfo{$anaID}[0].' > ';
						my $quantifNameStrg=($quantifID==$prevQuantifID)? '<SPAN style="visibility:hidden">'.$quantifName{$quantifID}.$ratioNameSep.'</SPAN>' : $quantifName{$quantifID}.$ratioNameSep;
						my ($ratioName,$rank)=@{$refRatio};
						my $valueStrg=$quantifID.'_'.$rank;
						print qq
						|<TR class="list" bgcolor="$bgColor">\n<TD nowrap><INPUT type="checkbox" name="test_quantif" value="$valueStrg" onclick="updateQuantifStatus(this.checked,$datasetIdx,'$valueStrg',true)">$sampleNameStrg$analysisNameStrg$quantifNameStrg$ratioName</TD>\n
						|;

						if (!$fromMotif ) {
							print "<TD><SPAN id=\"span_quantif_$datasetIdx\"><INPUT TYPE=\"text\" name=\"group_quantif\" id=\"group_quantif_$valueStrg\" placeholder=\"group for Anova\" disabled style=\"display:none\"></SPAN></TD>\n";
							print "</TR>\n";
						}

						$datasetIdx++;
						$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
						$prevSampleID=$sampleID;
						$prevAnaID=$anaID;
						$prevQuantifID=$quantifID;
					}
				}
			}
		}
	}

	print "</TABLE>\n";
	exit;
}

sub ajaxCopySelectionFromExpAna {
	my $expAnaID=$_[0];
	print header(-'content-encoding' => 'no',-charset => 'utf-8');
    warningsToBrowser(1);

	my (%quantifList,%refQuantifList);
	my $sthEA=$dbh->prepare("SELECT PARAM_LIST FROM EXPLORANALYSIS WHERE ID_EXPLORANALYSIS=$expAnaID AND PARAM_LIST LIKE '%//normalization=%'");
	$sthEA->execute;
	my ($paramList)=$sthEA->fetchrow_array;
	$sthEA->finish;

	if ($paramList) { # normalization
		my ($normStrg)=($paramList=~/normalization=([^\/]+)/);
		foreach my $strg (split(/[;:]/,$normStrg)) {
			print "$strg\n"; # testQuant%refQuant
		}
	}
	else {
		my $sthSelEA=$dbh->prepare("SELECT ID_QUANTIFICATION,TARGET_POS FROM EXPLORANA_QUANTIF WHERE ID_EXPLORANALYSIS=$expAnaID");
		$sthSelEA->execute;
		while (my ($quantifID,$targetPos)=$sthSelEA->fetchrow_array) {
			$targetPos=0 unless $targetPos;
			print $quantifID,'_',$targetPos,"\n";
		}
		$sthSelEA->finish;
	}
	$dbh->disconnect;
	exit;
}


sub fetchDesignQuantifications {
	my ($expID,$quantifFam,$refQuantifications,$refDesignName,$refQuantifName,$sthGetAllDesign,$sthGetExpCondName,$sthGetQuantifInfo,$sthGetMeas)=@_;

	$sthGetAllDesign->execute($expID);
	while (my ($designID,$desName) = $sthGetAllDesign -> fetchrow_array) {
		$refDesignName->{$designID}=$desName;
		$sthGetQuantifInfo -> execute($designID);
		my %fetchedCond;
		while (my($quantifID, $quantifAnnot, $quantName, $modifID) = $sthGetQuantifInfo -> fetchrow_array) {
			if ($sthGetMeas) {
				$sthGetMeas->execute($quantifID);
				my ($okQuantif)=$sthGetMeas->fetchrow_array;
				next unless $okQuantif;
			}
			$refQuantifName->{$quantifID} = $quantName;
			my ($labelTypeStrg,@labelInfo) = split('::',$quantifAnnot);
			foreach my $infoStrg (@labelInfo) {
				my ($setting, $valueStrg)=split('=',$infoStrg);
				if ($quantifFam eq 'RATIO' && $setting eq 'RATIOS') {
					my $ratioPos = 0;
					foreach my $ratio (split(";",$valueStrg)) {
						$ratioPos++;
						$ratio=~s/\#//g; # for design experiments, RATIOS are displayed with '#' for condition IDs
						my ($testID,$refID) = split(/\//,$ratio);
						$testID=~s/(%)\d+//;
						my $ratioTag = ($1)? '°' : '';
						$refID=~s/%\d+//;
						unless ($fetchedCond{$testID}) {
							$sthGetExpCondName->execute($testID);
							($fetchedCond{$testID})=$sthGetExpCondName->fetchrow_array;
						}
						unless ($fetchedCond{$refID}) {
							$sthGetExpCondName->execute($refID);
							($fetchedCond{$refID})=$sthGetExpCondName->fetchrow_array;
						}
						my $testName=$fetchedCond{$testID};
						my $refName=$fetchedCond{$refID};
						my $ratioItem = "$testName$ratioTag/$refName$ratioTag";#"."_$ratioPos";
						push @{$refQuantifications->{$designID}{$quantifID}},[$ratioItem,$ratioPos];
					}
					last;
				}
				elsif ($quantifFam eq 'MQ' && $setting eq 'STATES') {
					my $statePos=0;
					foreach my $state (split(";",$valueStrg)) {
						$statePos++;
						$state=~s/\#//g;
						my ($numBioRep,$quantiObsIDs,$expCondID)= split(',',$state);
						$sthGetExpCondName->execute($expCondID);
						my ($stateName)=$sthGetExpCondName->fetchrow_array;
						push @{$refQuantifications->{$designID}{$quantifID}},[$stateName,$statePos];
					}
					last;
				}
			}
		}
    }
}

sub ajaxPropAnnotate {
	print header(-type=>'text/plain',-charset => 'utf-8');
    warningsToBrowser(1);

    my ($propType,$propID)=(param('property')=~/property:(.):(\d+)/);

	my $sthPV=$dbh->prepare("SELECT DISTINCT BP.PROPERTY_VALUE FROM BIOSAMPLE_PROPERTY BP
                                INNER JOIN OBSERVATION O ON BP.ID_BIOSAMPLE=O.ID_BIOSAMPLE
                                INNER JOIN ANA_QUANTIFICATION AQ ON O.ID_ANALYSIS=AQ.ID_ANALYSIS
                                WHERE BP.ID_PROPERTY=$propID");
    #my ($propName)=$dbh->selectrow_array("SELECT NAME FROM PROPERTY WHERE ID_PROPERTY=$propID");
    $sthPV->execute;
    my %propValueCode;
    while (my ($propValue)=$sthPV->fetchrow_array) {
        my $propText=$propValue;
        if ($propType eq 'T') {
            $propValue=~s/^stepValue=\d:*//; # step independent: stepValue=1::quantity=xxx... -> quantity=xxx...
            $propValue='%' unless $propValue;
            $propText=&promsQuantif::convertTreatmentToText($propValue);
        }
        $propValueCode{$propValue}=$propText;
    }
    $sthPV->finish;


    foreach my $propValue (keys %propValueCode) {
        my %matchedQuantifs=&promsQuantif::getQuantificationsFromPropertyValue($dbh,$propID,$propType,$propValue,param('qList'));
        if (scalar keys %matchedQuantifs) {
			print $propValueCode{$propValue}.",".join(',',keys %matchedQuantifs)."\n";
        }
    }
    $dbh->disconnect;
    exit;
}

####>Revision history<####
# 2.1.8 Adds proteinID in exported matrices to prevent row duplicates due to different proteins having same ALIAS (PP 26/09/18)
# 2.1.7 minor changes (SL ??/??/18)
# 2.1.6 minor changes, change protein selection (SL 07/03/18)
# 2.1.5 Uses &promsMod::cleanDirectory to delete old exported files (PP 05/02/18)
# 2.1.4 fix minor bugs (SL 31/01/18)
# 2.1.3 add automaticaly biosample treatment or properties for Anova by group, add modifID for fetchQuantif function (SL ../10/17)
# 2.1.2 can use Anova and other statistical parameter for export and explorAnalysis (SL 12/10/17)
# 2.1.1 add phospho analysis pipeline and add call from motif for ajaxDisplayData (SL 23/01/17)
# 2.1.0 Adapted for MaxQuant (PP 09/01/17)
# 2.0.6 Fix uninitialized reference name when exporting normalized quantifs (PP 22/07/16)
# 2.0.5 Bug fix in +/-inf ratio values used as NA & other minor changes (PP 23/05/16)
# 2.0.4 Bug fix in quantification naming in export mode (PP 24/03/16)
# 2.0.3 Added option to import selection from a previous ExplorAna (PP 14/03/16)
# 2.0.2 TopN most variable proteins & compatible with all non-design quantifs (PP 29/02/16)
# 2.0.1 Compatible with SIN/emPAI (beta) (PP 29/01/16)
# 2.0.0 Also used for quantification data export<BR>TODO: Handle non-design quantifs (PP ../01/16)
# 1.1.0 Added number of peptides filter & record quantif pairing if proteome normalisation (PP 18/06/15)
# 1.0.9 fix bugs, close </SELECT> (Sl 30/03/15)
# 1.0.8 add ajax function to display data with or without modification (SL 08/03/15)
# 1.0.7 +/-inf no longer used p-value filter & filter on number of missing values (PP 23/01/15)
# 1.0.6 Compatible with PTM quantification & minor bug fix (PP 09/01/15)
# 1.0.5 check/uncheck cluster or/and pca; add protein number in FILERT_LIST (SL 04/09/14)
# 1.0.4 Rename matrtice.txt to matrix.txt (PP 07/08/14)
# 1.0.3 Minor changes (PP 09/07/14)
# 1.0.2 Move R code in lauchExploratoryAnalyses.pl (GA 08/07/14)
# 1.0.1 Display & code update (PP 07/07/14)
# 1.0.0 New script for PCA or Clustering protein set (SL 27/03/14)
