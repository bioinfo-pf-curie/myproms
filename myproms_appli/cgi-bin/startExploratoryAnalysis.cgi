#!/usr/local/bin/perl -w

################################################################################
# startExploratoryAnalysis.cgi       2.5.11                                    #
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
#print header(-charset=>'utf-8', -'content-encoding' => 'no'); # DEBUG
my %promsPath=&promsConfig::getServerInfo;
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my $dbh = &promsConfig::dbConnect;
my %proteinQuantifFamilies=&promsQuantif::getProteinQuantifFamilies;## For fetchQuantificationData function
my $userID=$ENV{'REMOTE_USER'};

my %pepTypeDesc=('NUM_PEP_USED'=>'All','DIST_PEP_USED'=>'Distinct','RAZ_UNI_PEP'=>'Razor + unique','UNIQUE_PEP'=>'Unique','IDENT_PEP'=>'Identified');
my $experimentID=&promsMod::cleanNumericalParameters(param('ID'));
my $action=param('ACT') || 'explorAna'; # or export
my $projectID=&promsMod::getProjectID($dbh,$experimentID,'EXPERIMENT');
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $ajax=param('AJAX') || "";

if ($ajax eq 'displaySelect') {
	my $fromMotif=param('CALL') || "";
	my $focus=(param('FOCUS'))? param('FOCUS') : ($fromMotif && param('MODIF'))? 1 : 0;
	&ajaxDisplayData($focus,param('MODIF'),param('QUANTIF_FAM'),$fromMotif);
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
	my ($title,$actionStrg)=($action eq 'export')? ('Exporting Quantifications','Fetching data') : ('Starting Exploratory Analyses','Launching Analyses');
    print qq
|<HTML>
<HEAD>
<TITLE>Preparing Exploratory Analysis or Export</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR>
<FONT class="title">$title</FONT>
<BR><BR><BR><BR><BR><SPAN id="waitSPAN"><IMG src="$promsPath{images}/scrollbarGreen.gif"></SPAN>
<BR><FONT class="title3">$actionStrg. Please wait...|;

    ### Get Parameter for insert
    my $PCAName = param('PCAName') || undef;
    my $clusterName = param('ClusterName') || undef;
	my ($metric,$itemMetric,$method);
    if ($clusterName) {
		($itemMetric,$metric) = split(/:/,param('metric'));
		$method = param('method');
    }
	
	#---- GLOBALS --->
	my ($pcaID,$clusteringID,$explorID,$focus,$quantifFam,$quantifMeasCode,$dataTransform,$numQuantifs,$isNormalized,$missingValues,$topNvarProteins,$restrictListID,$condRestrict);
	my ($numPepType,$minNumPeptides,$maxFoldChange,$minNumOkRatios,$maxPValue,$minNumOkPValues,$infRatios,$minFoldChange);
	my ($delocalize,$geneNameR,$useGeneName,$aggregate,$pValueAnovaChk,$pValueAnova,$nbProtChk,$nbProt);
	my ($exportRprocessMatrices,$anova,$groupList,$groupName,$imputeData,$imputeDataR,$nameFormat);
	my ($refSelectedProteins,$refExcludedProteins,$nbProtFilter);
	my (@selModifications,@selectedQuantifications,@testQuantifications,%test2refQuantif);
	my (%quantifPos,%isRefQuantif);
    my %transform=(	NONE=>sub {return $_[0];},
					LOG2=>sub {return log($_[0])/log(2);},
					LOG10=>sub {return log($_[0])/log(10);}
			       );

	my $referenceExplorID=&promsMod::cleanNumericalParameters(param('datasetSelection')); # set to 0 if undef
	if ($referenceExplorID) {
		##>Fetch reference info from DB
		my ($restrictListID,$refType,$strgParamList,$strgFilterList,$catExclusion)=$dbh->selectrow_array("SELECT ID_CATEGORY,ANA_TYPE,PARAM_LIST,FILTER_LIST,CAT_EXCLUSION FROM EXPLORANALYSIS WHERE ID_EXPLORANALYSIS=$referenceExplorID");
		if ($refType eq 'cluster') { # reset method & metric
			$strgParamList=~s/\/\/(method|metric)=[^\/]+//;
		}
		
		##>INSERT new Ana Exp in DB
		my $sthInsertExplorAna=$dbh->prepare("INSERT INTO EXPLORANALYSIS(ID_EXPLORANALYSIS,ID_CATEGORY,ID_EXPERIMENT,NAME,ANA_TYPE,PARAM_LIST,FILTER_LIST,CAT_EXCLUSION,STATUS,UPDATE_DATE,UPDATE_USER) VALUES (?,?,?,?,?,?,?,?,-1,NOW(),?)");
		#>PCA
		if (defined $PCAName) {
			($pcaID) = $dbh->selectrow_array("SELECT MAX(ID_EXPLORANALYSIS) FROM EXPLORANALYSIS");
			$pcaID++;
			$sthInsertExplorAna->execute($pcaID,$restrictListID,$experimentID,$PCAName,'PCA',$strgParamList,$strgFilterList,$catExclusion,$userID);
			$dbh->commit;
		}
		#>Clustering
		if (defined $clusterName) {
			$strgParamList .= "//method=$method//metric=$metric";
			($clusteringID)=$dbh->selectrow_array("SELECT MAX(ID_EXPLORANALYSIS) FROM EXPLORANALYSIS");
			$clusteringID++;
			$sthInsertExplorAna->execute($clusteringID,$restrictListID,$experimentID,$clusterName,'cluster',$strgParamList,$strgFilterList,$catExclusion,$userID);
			$dbh->commit;
		}
		$sthInsertExplorAna->finish;
	
		$explorID=$pcaID || $clusteringID;
		
		##>EXPLORANA_QUANTIF
		my $sthRefExplorQuantif=$dbh->prepare("SELECT ID_QUANTIFICATION,TARGET_POS FROM EXPLORANA_QUANTIF WHERE ID_EXPLORANALYSIS=$referenceExplorID");
		my $sthInsertExplorQuantif=$dbh->prepare("INSERT INTO EXPLORANA_QUANTIF(ID_QUANTIFICATION,ID_EXPLORANALYSIS,TARGET_POS) VALUES (?,?,?)");
		$sthRefExplorQuantif->execute;
		while (my ($quantifID,$targetPos)=$sthRefExplorQuantif->fetchrow_array) {
			$sthInsertExplorQuantif->execute($quantifID,$pcaID,$targetPos) if $pcaID;
			$sthInsertExplorQuantif->execute($quantifID,$clusteringID,$targetPos) if $clusteringID;
		}
		$sthRefExplorQuantif->finish;
		$sthInsertExplorQuantif->finish;
		$dbh->commit;		
	}
	else { # usual cases
		$focus = param('focus') || 0;
		if ($focus) {
			@selModifications=($focus==1)? param('stdModif') : param('toNormModif');
		}
		($quantifFam,$quantifMeasCode)=split(':',param('quantifFam'));
		$quantifMeasCode='' unless $quantifMeasCode;
		my $listParam = param('list') || "";
		$dataTransform=param('dataTransform') || 'NONE';
		my @selQuantifs=param('test_quantif');
		$numQuantifs=scalar @selQuantifs;
		my @groupQuantifs=param('group_quantif');
	
		$isNormalized=0;
		if ($focus == 2) {##NORMALIZED PHOSPHO
			my @normQuantifs=param('ref_quantif');
			foreach my $i (0..$#selQuantifs) {
				$selQuantifs[$i].='%'.$normQuantifs[$i];
			}
			$isNormalized=1;
		}
		my $featureType=($focus)? 'siteQuant' : 'protQuant';
		my $strgParamList="//feature=$featureType//quantFam=$quantifFam";
		$strgParamList.="//quantCode=$quantifMeasCode" if $quantifMeasCode;
		$strgParamList.="//dataTrans=$dataTransform";
		$missingValues=param('missingValues');
		$topNvarProteins=param('topNvar') || 0;
		$restrictListID = ($listParam)? $listParam : undef;
		$condRestrict = ($listParam)? param('condRestrict') : "";
		my $catExclusion = ($condRestrict eq 'exclude')? 1 : 0;
		$numPepType=param('pepType') || 'NUM_PEP_USED';
		$minNumPeptides=param('numPep') || 1;
		#>RATIO-specific parameters
		$maxFoldChange = param('foldChange') || 1;
		$minNumOkRatios = param('foldChangeOcc') || 1;
		$maxPValue = param('pValue') || 1;
		$minNumOkPValues = param('pValueOcc') || 1;
		$infRatios=param('infRatios') || 0;
		$minFoldChange = 1/$maxFoldChange;
	
		####>Get parameters for data post-processing by R
		#>FILTERING
		$delocalize=param('delocalize')? "TRUE" : "FALSE";
		$geneNameR="FALSE";
		$useGeneName=param('geneName') || 0;
		$aggregate=param('aggregate')? "TRUE" : "FALSE";
		$pValueAnovaChk=param('pvalChkAnova')? "TRUE" : "FALSE";
		$pValueAnova=($pValueAnovaChk eq "TRUE")? param('pValueAnova') : "FALSE";
		$nbProtChk=param('topNvarChk')? "TRUE" : "FALSE";
		$nbProt=($nbProtChk eq "TRUE")? param('topNvar') : "FALSE";
		#>STAT-specific parameters
		$exportRprocessMatrices=($nbProtChk eq 'TRUE' || $aggregate eq 'TRUE')? 1 : 0; # default
		$anova = (param('anovaSel'))? param('anovaSel') : (param('topNvarChk'))? "sample" : "none";
	$groupList=param('GROUP_LIST'); # $anova = 'group'
	$groupName=param('GROUP_ANNOT'); # $anova = 'group'
		($imputeData,$imputeDataR)=($action eq 'explorAna' || param('imputeNA'))? (1,"TRUE") : (0,"FALSE"); # imputation is mandatory for explor analysis
		$nameFormat=param('nameFormat') || 'OPTIMAL'; # 'EXTENDED' # for export
		
		if ($aggregate eq "TRUE" && $nbProtChk eq "FALSE") {
			$anova="sample";
			if ($action eq "export") {
				$nbProtChk="TRUE";
				$nbProt=-1;
			}
		}
		
	#	my ($groupName,%groupQuantif);
	#    if ($anova eq "group") {
	#		my $groupList=param('GROUP_LIST');
	#		$groupName=param('GROUP_ANNOT');
	#		my %groupQuantifNames=&promsQuantif::getDistinctQuantifNames($dbh,join(':',param('test_quantif')));
	#		foreach my $list (split('::',$groupList)) {
	#			my ($group, $quantifStrg)=split('=', $list);
	#			foreach my $quantif (split('\|', $quantifStrg)) {
	#				(my $quantifName=$groupQuantifNames{$nameFormat}{$quantif})=~s/\W+/\./g;
	#				$quantifName=~s/^\.+//; $quantifName=~s/\.+\Z//;
	#				push @{$groupQuantif{$group}},$quantifName;
	#			}
	#		}
	#    }
	
		my $strgFilterList="//PEP=$numPepType:$minNumPeptides//NA=$missingValues";
		$strgFilterList.="//TOPN_VAR=$topNvarProteins" if $topNvarProteins;
		if ($quantifFam eq 'RATIO') {
			$strgFilterList.="//PV=$maxPValue:$minNumOkPValues" if $focus != 2; # not for normalized data
			$strgFilterList.="//FC=$maxFoldChange:$minNumOkRatios//INF=$infRatios";
		}
	
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
	
		if ($restrictListID) {
			my %listProteins;
			my $sthCP = $dbh->prepare("SELECT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$restrictListID");
			$sthCP->execute;
			while (my ($protID) = $sthCP->fetchrow_array) {
				$listProteins{$protID} = 1;
			}
			$sthCP -> finish;
			$nbProtFilter=scalar(keys %listProteins);
			$strgFilterList.="//nbProt=$nbProtFilter";
			($refSelectedProteins,$refExcludedProteins)=($condRestrict eq 'restrict')? (\%listProteins,undef) : (undef,\%listProteins);
		}
		
		if ($action eq 'explorAna') {
	
			##>INSERT into database
			my $sthInsertExplorAna = $dbh->prepare("INSERT INTO EXPLORANALYSIS(ID_EXPLORANALYSIS,ID_CATEGORY,ID_EXPERIMENT,NAME,ANA_TYPE,PARAM_LIST,FILTER_LIST,CAT_EXCLUSION,STATUS,UPDATE_DATE,UPDATE_USER) VALUES (?,?,?,?,?,?,?,?,-1,NOW(),?)");
			#$pcaID = $clusteringID = 1; ## tocomment
			$strgParamList.="//$strgNorm" if $strgNorm;
	
			##>PCA
			if (defined $PCAName) {
				($pcaID) = $dbh->selectrow_array("SELECT MAX(ID_EXPLORANALYSIS) FROM EXPLORANALYSIS");
				$pcaID++;
				$sthInsertExplorAna->execute($pcaID,$restrictListID,$experimentID,$PCAName,'PCA',$strgParamList,$strgFilterList,$catExclusion,$userID);
				$dbh->commit;
			}
			
			##>Clustering
			if (defined $clusterName) {
				$strgParamList .= "//method=$method//metric=$metric";
				($clusteringID) =  $dbh->selectrow_array("SELECT MAX(ID_EXPLORANALYSIS) FROM EXPLORANALYSIS");
				$clusteringID++;
				$sthInsertExplorAna -> execute($clusteringID,$restrictListID,$experimentID,$clusterName,'cluster',$strgParamList,$strgFilterList,$catExclusion,$userID);
				$dbh->commit;
			}
			$sthInsertExplorAna -> finish;
	
			$explorID=$pcaID || $clusteringID;
			
			##>EXPLORANA_QUANTIF
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
			$sthInsertExplorQuantif->finish;
			$dbh->commit;
		}
	
	} # end of ! $referenceExplorID

	$dbh->disconnect;
	
	if ($action eq 'explorAna') {

		####>Forking to continue as background process<####
		my $childExplor = fork;
		if ($childExplor) { # parent here
			print qq
| Done.</FONT><BR><BR><BR>
<FONT class="title2">Your exploratory analysis job is running. Analysis status will be updated when finished.</FONT>
<INPUT type="button" class="title2" value="Continue" onclick="displayExplorAna()"/>
<SCRIPT type="text/javascript">
document.getElementById('waitSPAN').style.display='none';
function displayExplorAna() {
	top.promsFrame.selectedAction = 'summary';
	parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=experiment:$experimentID&branchID=exploranalysis:$explorID&VIEW=explorAna";
}
</SCRIPT>
</BODY>
</HTML>
|;		
			exit;
		}
		
		#------> child here
		
		#>Disconnecting from server
		open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
		open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
		open STDERR, '>/dev/null' or die "Can't open /dev/null: $!";
		
		#>Reconnect to DB
		#$dbh=&promsConfig::dbConnect('no_user');
		
		###>Directories<####
		mkdir "$promsPath{tmp}" unless -e "$promsPath{tmp}";
		mkdir "$promsPath{explorAna}/project_$projectID" unless -e "$promsPath{explorAna}/project_$projectID";
		mkdir "$promsPath{tmp}/exploratory_analysis" unless -e "$promsPath{tmp}/exploratory_analysis";
		
		mkdir "$promsPath{tmp}/exploratory_analysis/$pcaID" if ($pcaID && !-e "$promsPath{tmp}/exploratory_analysis/$pcaID");
		mkdir "$promsPath{tmp}/exploratory_analysis/$clusteringID" if ($clusteringID && !-e "$promsPath{tmp}/exploratory_analysis/$clusteringID");
	}
	
#####>----------If exploratory analysis: background process from now on...

	my ($exportDir,$exportPath,$quantifFile,$pepFile,$pepFileRef,@file2Tar);
	
	if ($referenceExplorID) {
		##>Copy preprocessing files from reference Exploratory analysis
		foreach my $eaID ($pcaID,$clusteringID) {
			next unless $eaID;
			foreach my $dataFile ('matrix.txt','missingValues.txt','valueDistribution.png') {
				copy("$promsPath{data}/exploratory_data/project_$projectID/$referenceExplorID/$dataFile","$promsPath{tmp}/exploratory_analysis/$eaID/$dataFile") if -e "$promsPath{data}/exploratory_data/project_$projectID/$referenceExplorID/$dataFile";
			}
		}
	}
	else {
		#>Reconnect to DB
		$dbh=&promsConfig::dbConnect('no_user');
		
		####>Fetching quantification data (moved up to get software info. PP 27/03/19)<####
		my (%quantifValues,%quantifInfo,%proteinInfo,%groupQuantifNames,%dispModifSites);
		my ($view,$refProtInfo)=($action eq 'explorAna')? ('explorAna',undef) : ('export',\%proteinInfo);
		#my $absValMod = ($isNormalized)? abs($valMod)  : $valMod;
		my %parameters=(QUANTIF_FAMILY=>$quantifFam,VIEW=>$view,NUM_PEP_CODE=>$numPepType,QUANTIF_LIST=>\@selectedQuantifications,VERBOSE=>0); #,SEL_MODIF_ID=>$absValMod
		$parameters{MEASURE}=$quantifMeasCode if $quantifMeasCode; # only for non-RATIO quantifs
		&promsQuantif::fetchQuantificationData($dbh,\%parameters,\%quantifInfo,\%quantifValues,\%dispModifSites,$refProtInfo,$refSelectedProteins,$refExcludedProteins);
		my %quantifSoftwares;
		my ($okPvalue,$okSdGeo)=(0,0);
		if ($quantifFam eq 'RATIO') {
			foreach my $quantif (@selectedQuantifications) {
				my ($quantifID,$targetPos)=split('_',$quantif);
				@{$quantifSoftwares{$quantif}}=('myProMS',1);
				if ($quantifInfo{$quantifID}[1]->{SOFTWARE}) {
					$quantifSoftwares{$quantif}[0]=$quantifInfo{$quantifID}[1]->{SOFTWARE}[0];
					$quantifSoftwares{$quantif}[0]='MaxQuant' if $quantifSoftwares{$quantif} eq 'MQ';
					$quantifSoftwares{$quantif}[1]=$quantifInfo{$quantifID}[1]->{SOFTWARE}[1] if $quantifInfo{$quantifID}[1]->{SOFTWARE}[1];
				}
				if (!$isNormalized) {
					$okPvalue=1 if $quantifSoftwares{$quantif}[0] ne 'MaxQuant';
					$okSdGeo=1 if ($view eq 'export' && $quantifSoftwares{$quantif}[0] eq 'myProMS' && $quantifSoftwares{$quantif}[1] >= 2);
				}
			}
		}
		
		if ($action eq 'explorAna') {		
			my $numMissValue=scalar(@testQuantifications)*$missingValues/100;
			###>Directory for PCA
			if ($pcaID) {
				##mkdir "$promsPath{tmp}/exploratory_analysis/$pcaID" unless -e "$promsPath{tmp}/exploratory_analysis/$pcaID";
				open(MAT_PCA,">$promsPath{tmp}/exploratory_analysis/$pcaID/matrix.txt");
				#open(NA_PCA,">$promsPath{tmp}/exploratory_analysis/$pcaID/missingValues.txt");
				open(R_PARAMS,">$promsPath{tmp}/exploratory_analysis/$pcaID/R_parameters.txt");
				print R_PARAMS "EXCLUDE_AMB\tGENE_NAME\tPROTEIN_SELECTION\tAGGREGATE\tKEEP_PROT\tKEEP_PROT_NB\tANOVA_P_VALUE_CHK\tANOVA_P_VALUE\n";
				print R_PARAMS "$delocalize\t$geneNameR\t$anova\t$aggregate\t$nbProtChk\t$nbProt\t$pValueAnovaChk\t$pValueAnova\n";
				close R_PARAMS;
			}
			###>Directory for CLUSTER
			if ($clusteringID) {
				##mkdir "$promsPath{tmp}/exploratory_analysis/$clusteringID" unless -e "$promsPath{tmp}/exploratory_analysis/$clusteringID";
				open(MAT_CLUSTER,">$promsPath{tmp}/exploratory_analysis/$clusteringID/matrix.txt");
				#open(NA_CLUSTER,">$promsPath{tmp}/exploratory_analysis/$clusteringID/missingValues.txt");
				open(R_PARAMS,">$promsPath{tmp}/exploratory_analysis/$clusteringID/R_parameters.txt");
				print R_PARAMS "EXCLUDE_AMB\tGENE_NAME\tPROTEIN_SELECTION\tAGGREGATE\tKEEP_PROT\tKEEP_PROT_NB\tANOVA_P_VALUE_CHK\tANOVA_P_VALUE\n";
				print R_PARAMS "$delocalize\t$geneNameR\t$anova\t$aggregate\t$nbProtChk\t$nbProt\t$pValueAnovaChk\t$pValueAnova\n";
				close R_PARAMS;
			}
		}
		else { # export
			###>Config<####
			mkdir "$promsPath{tmp}/scratch" unless -d "$promsPath{tmp}/scratch";
			&promsMod::cleanDirectory("$promsPath{tmp}/scratch/export",'15m');
			my $jobID=strftime("%Y%m%d%H%M%S",localtime);
			$exportDir=$jobID."_$userID";
			mkdir "$promsPath{tmp}/scratch/export/$exportDir";
			$exportPath="$promsPath{tmp}/scratch/export/$exportDir";
			$quantifFile=($quantifFam eq 'RATIO')? 'ratio.txt' : ($quantifMeasCode eq 'MY_LFQ')? 'LFQ.txt' : "$quantifMeasCode.txt";
			my $processedPrefix=($exportRprocessMatrices)? 'processed_' : '';
			
			push @file2Tar, ($exportRprocessMatrices || $imputeData)? "processed_$quantifFile" : $quantifFile;
			open(MAT_EXPORT_QUANTIF,">$exportPath/$quantifFile");
			open(MAT_INFO_EXPORT,">$exportPath/value_info.txt") if ($exportRprocessMatrices || $imputeData);
			if (!$isNormalized) {
				$pepFile=($numPepType=~/NUM_PEP_USED|IDENT_PEP/)? 'peptide.txt' : 'dist_peptide.txt';
				open(MAT_EXPORT_PEP,">$exportPath/$pepFile");
				push @file2Tar,$processedPrefix.$pepFile;
				if ($quantifFam eq 'RATIO') {
					open(MAT_EXPORT_PVAL,">$exportPath/pvalue.txt") if $okPvalue;
					open(MAT_EXPORT_SDGEO,">$exportPath/sd_geo.txt") if $okSdGeo;
					push @file2Tar,$processedPrefix.'pvalue.txt' if $okPvalue;
					push @file2Tar,$processedPrefix.'sd_geo.txt' if $okSdGeo;
				}
			}
	
			push @file2Tar, "parameters.txt";
			push @file2Tar,$processedPrefix."annotation.txt";
			if ($exportRprocessMatrices || $imputeData) {
				push @file2Tar, "processed_value_info.txt";
				push @file2Tar, "processed_value_info.png";
			}
			open(ANNOT,">$exportPath/annotation.txt");
			open(PARAMS,">$exportPath/parameters.txt");
			
			my $focusStrg='Proteins';
			if ($focus) {
				$focusStrg='';
				#my ($modifName)=$dbh->selectrow_array("SELECT PSI_MS_NAME FROM MODIFICATION WHERE ID_MODIFICATION=$absValMod");
				my $sthSelModif=$dbh->prepare("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES FROM MODIFICATION WHERE ID_MODIFICATION=?");
				foreach my $modID (@selModifications) {
					$sthSelModif->execute($modID);
					my ($psiName,$interName,$synName)=$sthSelModif->fetchrow_array;
					my $modifName=$psiName || $interName || $synName;
					$modifName=~s/^##//; $modifName=~s/##.*$//;
					$focusStrg.='/' if $focusStrg;
					$focusStrg.=$modifName;
				}
				$focusStrg=$focusStrg.'-Sites';
			}
			if ($exportRprocessMatrices || $imputeData) {
				push @file2Tar, "R_parameters.txt";
				my $numMissValue=scalar(@testQuantifications)*$missingValues/100;
				open(R_PARAMS,">$exportPath/R_parameters.txt");
				print R_PARAMS "FOCUS\tQUANTIF_FAM\tPROTEIN_SELECTION\tAGGREGATE\tKEEP_PROT\tKEEP_PROT_NB\tANOVA_P_VALUE_CHK\tANOVA_P_VALUE\tMISSING_VALUE\tIMPUTE\n";
				print R_PARAMS "$focus\t$quantifFam\t$anova\t$aggregate\t$nbProtChk\t$nbProt\t$pValueAnovaChk\t$pValueAnova\t$numMissValue\t$imputeDataR\n";
				close(R_PARAMS);
			}
	
			my $listFilterStrg='None';
			if ($restrictListID) {
				my ($listName)=$dbh->selectrow_array("SELECT CONCAT(CL.NAME,' > ',CA.NAME) FROM CATEGORY CA,CLASSIFICATION CL WHERE CA.ID_CLASSIFICATION=CL.ID_CLASSIFICATION AND CA.ID_CATEGORY=$restrictListID");
				$listFilterStrg="$condRestrict: $listName ($nbProtFilter proteins)";
			}
			my $infRatioStrg=($infRatios < 0)? '100=MISSING_VALUE' : $infRatios;
	
			print PARAMS qq
|FOCUS\t$focusStrg
QUANTIF_METHOD\t$quantifFam
DATA_TRANSFORM\t$dataTransform
%_MISSING_VALUES_ALLOWED\t$missingValues
LIST_FILTERING\t$listFilterStrg
MIN_NUM_PEPTIDES\t$minNumPeptides
PEPTIDE_TYPE\t$pepTypeDesc{$numPepType}
DATA_SELECTION\t$anova
IMPUTE_MISSING_VALUES\t$imputeDataR
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
			if ($quantifFam eq 'RATIO') {
				print PARAMS qq
|MAX_ABS_RATIO\t$maxFoldChange
MIN_RATIO_OCCURENCE\t$minNumOkRatios
MINUS_INFINITE_RATIO\t-1000
PLUS_INFINITE_RATIO\t1000
MAX_P_VALUE\t$maxPValue
MIN_P_VALUE_OCCURENCE\t$minNumOkPValues
%_INFINITE_RATIOS_ALLOWED\t$infRatioStrg
|;
				if ($focus) {
					my $normStrg=($isNormalized)? 'Yes' : 'No';
					print PARAMS qq
|AGGREGATE_SITES\t$aggregate
EXCLUDE_AMBIGOUS_SITES\t$delocalize
RATIOS_NORMALIZED_DURING_EXPORT\t$normStrg
|;
				}
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

			##>Quantifications used with hierarchy (reference quantifs are obsolete)
			my $sthPH=$dbh->prepare("SELECT CONCAT (P.NAME,' > ',E.NAME,' > ',D.NAME) FROM QUANTIFICATION Q
												INNER JOIN DESIGN D ON Q.ID_DESIGN=D.ID_DESIGN
												INNER JOIN EXPERIMENT E ON D.ID_EXPERIMENT=E.ID_EXPERIMENT
												INNER JOIN PROJECT P ON E.ID_PROJECT=P.ID_PROJECT
												WHERE Q.ID_QUANTIFICATION=?");
			my $firstQuantifID=(split('_',$selectedQuantifications[0]))[0];
			$sthPH->execute($firstQuantifID);
			my ($hierarchyStrg)=$sthPH->fetchrow_array;
			$sthPH->finish;
			my $sthQN=$dbh->prepare("SELECT ID_QUANTIFICATION,NAME FROM QUANTIFICATION WHERE ID_QUANTIFICATION IN (".join(',',keys %quantifPos).")");
			$sthQN->execute;
			while (my ($quantifID,$qName)=$sthQN->fetchrow_array) {
				print PARAMS "QUANTIFICATION USED: $hierarchyStrg > $qName [UID: $quantifID]\n";
			}
			$sthQN->finish;

			close PARAMS;
		}
	
		my @header;
		if ($action eq 'explorAna') {
			foreach my $testQuant (@testQuantifications) {
				if ($test2refQuantif{$testQuant}) {
					push @header,$testQuant.'%'.$test2refQuantif{$testQuant};
				}
				else {
					push @header,$testQuant;
				}
			}
		}
		else { # export
			my %testQuantifNames=&promsQuantif::getDistinctQuantifNames($dbh,join(':',param('test_quantif')));
			my %refQuantifNames=($isNormalized)? &promsQuantif::getDistinctQuantifNames($dbh,join(':',param('ref_quantif'))) : ();
			my %usedHeaders;
			foreach my $quantif (@testQuantifications) {
				(my $headerStrg=$testQuantifNames{$nameFormat}{$quantif})=~s/\W+/\./g;
				$headerStrg=~s/^\.+//; $headerStrg=~s/\.+\Z//;
				my $refHeaderStrg=$headerStrg;
				if ($usedHeaders{$refHeaderStrg}) { # same label already used
					$headerStrg.='__uid'.$quantif ; # to prevent duplicates
					my ($firstIdx,$firstUid)=@{$usedHeaders{$refHeaderStrg}};
					if ($firstIdx > -1) { # Also update 1st instance of $refHeaderStrg (only once!)
						if ($isNormalized) {$header[$firstIdx]=~s/\.VS\./__uid$firstUid\.VS\./;}
						else {$header[$firstIdx].='__uid'.$firstUid;}
						$usedHeaders{$refHeaderStrg}[0]=-1; # set Idx to -1 once 1st instance has been updated
					}
				}
				else {@{$usedHeaders{$refHeaderStrg}}=(scalar @header,$quantif);} # record index in @header when seen for the 1rst time (in case duplicates follow)
				my $strgHeader=$headerStrg;
				if ($isNormalized) {
					(my $refStrg=$refQuantifNames{$nameFormat}{$test2refQuantif{$quantif}})=~s/\W+/\./g;
					$refStrg=~s/^\.+//; $refStrg=~s/\.+\Z//;
					$headerStrg.='.VS.'.$refStrg;
				}
				#$groupQuantifNames{$strgHeader}=$headerStrg;
				$groupQuantifNames{$quantif}=$headerStrg; # in case of anova group
				push @header,$headerStrg;
			}
	
			if ($anova eq "group") { # action = export
				my %groupQuantif;
				foreach my $list (split('::',$groupList)) {
					my ($group, $quantifStrg)=split('=', $list);
					foreach my $quantif (split('\|', $quantifStrg)) {
						push @{$groupQuantif{$group}},$quantif;
					}
				}
				push @file2Tar, "groups.txt";
				push @file2Tar, "anova_pvalue.txt";
				open(GROUP,">$exportPath/groups.txt");
				print GROUP "SAMPLE\tGROUP\n";
				foreach my $group (sort{$a cmp $b} keys %groupQuantif) {
					foreach my $quantif (@{$groupQuantif{$group}}) {
						print GROUP "$groupQuantifNames{$quantif}\t$group\n";
					}
				}
				close GROUP;
			}
		
			if ($anova eq "sample") {  # action = export
				push @file2Tar, "sd.txt";
			}
		}
		
		my $strgCond = join("\t",@header);
		my $numTestQuantifs=scalar @testQuantifications;
		my $numAllowedInf=($infRatios < 0)? $numTestQuantifs : $numTestQuantifs*$infRatios/100;
		my $numAllowedMissing=$numTestQuantifs*$missingValues/100;
		print MAT_PCA "\t$strgCond\n" if $pcaID;
		print MAT_CLUSTER "\t$strgCond\n" if $clusteringID;
		
		if ($action eq 'export') {
			print MAT_EXPORT_QUANTIF "\t$strgCond\n";
			print MAT_INFO_EXPORT "\t$strgCond\n" if ($exportRprocessMatrices || $imputeData);
			if (!$isNormalized) {
				print MAT_EXPORT_PEP "\t$strgCond\n";
				if ($quantifFam eq 'RATIO') {
					print MAT_EXPORT_PVAL "\t$strgCond\n";
					print MAT_EXPORT_SDGEO "\t$strgCond\n" if $okSdGeo;
				}
			}
			print ANNOT "\tIDENTIFIER\tENSEMBL_GENEID\tENTREZ_GENEID\tGENE_NAME\tSYNONYMS\tDESCRIPTION\n";
		}
	
		print '/.' if $action eq 'export';
	
		my (%protEntryLabel,%extraAnnotation);
		if ($action eq 'export') { # Pre-process row labels in case of duplicate entry due to identifier mapping
			my %seenProtEntry;
			foreach my $modProtID (keys %quantifValues) {
				my ($proteinID,$modStrg)=($modProtID=~/^(\d+)-*(.*)/);
				my $entryName;
				if ($useGeneName) {$entryName=$proteinInfo{$proteinID}[5][0] || 'NO_GENE';} # gene
				else {$entryName=$proteinInfo{$proteinID}[0];} # ALIAS
				$entryName.= '-'.$dispModifSites{$modProtID} if $modStrg; # site
				if ($seenProtEntry{$entryName}) { # already seen
					my $firstModProtID=$seenProtEntry{$entryName};
					$protEntryLabel{$firstModProtID}.='__uid'.$firstModProtID unless $protEntryLabel{$firstModProtID}=~/__uid$firstModProtID$/; # update original entry (only once!)
					$entryName.='__uid'.$modProtID;
				}
				else {$seenProtEntry{$entryName}=$modProtID;}
				$protEntryLabel{$modProtID}=$entryName;
			}
			#>Extra annotation
			my ($ensemblGeneID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='Ensembl'");
			my ($entrezGeneID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GeneID'");
			my $sthExAnnot=$dbh->prepare("SELECT P.ID_PROTEIN,MI.ID_IDENTIFIER,GROUP_CONCAT(MI.VALUE) FROM PROTEIN P
											INNER JOIN MASTERPROT_IDENTIFIER MI ON P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN
											WHERE P.ID_PROTEIN IN (".join(',',keys %proteinInfo).") AND MI.ID_IDENTIFIER IN ($ensemblGeneID,$entrezGeneID)
											GROUP BY P.ID_PROTEIN,MI.ID_IDENTIFIER");
			$sthExAnnot->execute;
			while (my($proteinID,$identID,$valueStrg)=$sthExAnnot->fetchrow_array) {
				@{$extraAnnotation{$proteinID}}=('','') unless $extraAnnotation{$proteinID};
				if ($identID==$ensemblGeneID) { # Ensembl GeneID
					foreach my $value (split(',',$valueStrg)) {
						if ($value=~/G\d+(\.|$)/) {
							$value=~s/\..*//; # Clean isoform info
							$extraAnnotation{$proteinID}[0]=$value;
							last;
						}
					}
				}
				else { # Entrez GeneID
					$extraAnnotation{$proteinID}[1]=(split(',',$valueStrg))[0]; # just to be safe (only 1 value there)
				}
			}
			$sthExAnnot->finish;
		}
	
		my $nbAllProt=0;
		my $count=0;
		my ($countMinusInf,$countPlusInf,$countMissingvalues,$countTrueValues)=0;
		my (%missingValuesProt); # %usedProteinInfo,
		foreach my $modProtID (keys %quantifValues) { # skip sort because of modif quantif -> sort{lc($proteinInfo{$a}[0]) cmp lc($proteinInfo{$b}[0])} # sort by ALIAS
			my ($proteinID,$modStrg)=($modProtID=~/^(\d+)-*(.*)/);
			$modStrg='' unless $modStrg;
			#>Scan for +/- inf & missing values & peptide filter
			my (%infiniteRatios,%missingValueQuantif,%usedProtRatio,%isInfiniteRatios); #,@trueValues
			my $numTrueValues=0;
			foreach my $quantif (@testQuantifications) {
				my $strgQuantif=($isNormalized)? $quantif."%".$test2refQuantif{$quantif} : $quantif;
				$count++;
				if ($count==150000) {
					$count=0;
					print '.' if $action eq 'export';
				}
				if ($test2refQuantif{$quantif}) { # Protein-level normalization!!! no p-value filter if ratio normalization (final ratio is unrelated to recorded p-value!); +/-Inf ratio not allowed for reference
					if (!$quantifValues{$proteinID} || !$quantifValues{$proteinID}{$test2refQuantif{$quantif}}{$quantifMeasCode}
						 || ($quantifFam eq 'RATIO' && ($quantifValues{$proteinID}{$test2refQuantif{$quantif}}{$quantifMeasCode}==0.001 || $quantifValues{$proteinID}{$test2refQuantif{$quantif}}{$quantifMeasCode}==1000))
						 || ($quantifFam !~ /SIN|EMPAI|PROT_RULER/ && (!$quantifValues{$proteinID}{$test2refQuantif{$quantif}}{$numPepType} || $quantifValues{$proteinID}{$test2refQuantif{$quantif}}{$numPepType} < $minNumPeptides))
					) {
						$missingValueQuantif{$strgQuantif}=1;
						next;
					}
				}
				if ($quantifValues{$modProtID}{$quantif}) {
					if (!$isNormalized) {
						if (!$quantifValues{$modProtID}{$quantif}{$numPepType} || $quantifValues{$modProtID}{$quantif}{$numPepType} < $minNumPeptides) { # peptide filters
							$missingValueQuantif{$strgQuantif}=1;
							next;
						}
					}
					if ($quantifFam eq 'RATIO') {
						if ($quantifValues{$modProtID}{$quantif}{'RATIO'} == 0.001) {
							if ($infRatios < 0) { # use as missing value
								$missingValueQuantif{$strgQuantif}=1;
								$isInfiniteRatios{$strgQuantif}="-Inf";
							}
							else {
								$infiniteRatios{$quantif}=-1;
								$usedProtRatio{$quantif}=-1000;
							}
						}
						elsif ($quantifValues{$modProtID}{$quantif}{'RATIO'} == 1000) {
							if ($infRatios < 0) { # use as missing value
								$missingValueQuantif{$strgQuantif}=1;
								$isInfiniteRatios{$strgQuantif}="+Inf";
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
					else { # MQ/PROT_RULER/PROT_ABUNDANCE/SIN/EMPAI (ok ref if any)
						#$usedProtRatio{$quantif}=$quantifValues{$modProtID}{$quantif}{$quantifMeasCode};
						$usedProtRatio{$quantif}=($test2refQuantif{$quantif})? $quantifValues{$modProtID}{$quantif}{$quantifMeasCode}/$quantifValues{$proteinID}{$test2refQuantif{$quantif}}{$quantifMeasCode} : $quantifValues{$modProtID}{$quantif}{$quantifMeasCode};
						$numTrueValues++;
					}
				}
				else {$missingValueQuantif{$strgQuantif}=1;}
			}
	
			#my $numTrueValues=scalar @trueValues;
			my $numProtMissValues=scalar keys %missingValueQuantif; # $numProtMissValues can equal $numQuantifs due to peptide filtering
			next if ($numProtMissValues==$numQuantifs || $numProtMissValues > $numAllowedMissing || ($quantifFam eq 'RATIO' && scalar (keys %infiniteRatios) > $numAllowedInf) || ($action eq 'explorAna' && $numTrueValues < 2));
			%{$missingValuesProt{$modProtID}}=%missingValueQuantif if $numProtMissValues;
	
			#$maxTrueRatio*=2; # !!!!!!!!!!!!!!!!!!! +/-inf ratios converted into 2x max ratio for protein accross all quantifs !!!!!!!!!!!!!!!!!
	
			#>Re-scan for ratio & p-value filter
			if ($quantifFam eq 'RATIO') {
				my ($okRatio,$okPvalue)=(0,0);
				foreach my $quantif (@testQuantifications) {
					next if $missingValueQuantif{$quantif};
					$okRatio++ if ($maxFoldChange==1 || $usedProtRatio{$quantif} <= $minFoldChange || $usedProtRatio{$quantif} >= $maxFoldChange);
					$okPvalue++ if ($maxPValue==1 || $quantifSoftwares{$quantif}[0] eq 'MaxQuant' || ($quantifValues{$modProtID}{$quantif}{'P_VALUE'} && $quantifValues{$modProtID}{$quantif}{'P_VALUE'} <= $maxPValue)); # !!!Does not rely on +/-inf to validate p-value filter!!! (+/-inf accepted only if another quantif has validated filter)
				}
				next if ($okRatio < $minNumOkRatios || $okPvalue < $minNumOkPValues);
			}
	
			$nbAllProt++;
			next if ($delocalize eq "TRUE" && $modStrg=~/~/);
			print MAT_PCA "$modProtID" if $pcaID;
			print MAT_CLUSTER "$modProtID" if $clusteringID;
			
			if ($action eq 'export') {
				#my $dispModStrg=($modStrg)? '-'.$dispModifSites{$modProtID} : '';
				#my $entryName;
				#if ($useGeneName) {$entryName=$proteinInfo{$proteinID}[5][0] || 'NO_GENE';} # gene
				#else {$entryName=$proteinInfo{$proteinID}[0];} # ALIAS
				#$entryName.=$dispModStrg.'__uid'.$modProtID;
				print MAT_EXPORT_QUANTIF $protEntryLabel{$modProtID};
				print MAT_INFO_EXPORT $protEntryLabel{$modProtID} if ($exportRprocessMatrices || $imputeData);
				if (!$isNormalized) {
					print MAT_EXPORT_PEP $protEntryLabel{$modProtID};
					if ($quantifFam eq 'RATIO') {
						print MAT_EXPORT_PVAL $protEntryLabel{$modProtID};
						print MAT_EXPORT_SDGEO $protEntryLabel{$modProtID} if $okSdGeo;
					}
				}
			}
			foreach my $quantif (@testQuantifications) {
				if ($usedProtRatio{$quantif}) { # usable value
					if ($quantifFam eq 'RATIO') {
						if ($infiniteRatios{$quantif}) { # flag +/-inf ratios
							$countMinusInf++ if ($usedProtRatio{$quantif} == -1000);
							$countPlusInf++ if ($usedProtRatio{$quantif} == 1000);
							my $infValue=($usedProtRatio{$quantif} == -1000)? "-Inf" : "+Inf";
							if ($pcaID){
								print MAT_PCA "\t$usedProtRatio{$quantif}";
							}
							if ($clusteringID){
								print MAT_CLUSTER "\t$usedProtRatio{$quantif}";
							}
							if ($action eq 'export') {
								print MAT_EXPORT_QUANTIF "\t$usedProtRatio{$quantif}";
								print MAT_INFO_EXPORT "\t$infValue" if ($exportRprocessMatrices || $imputeData);
								print MAT_EXPORT_PVAL "\tNA" unless $isNormalized;
								print MAT_EXPORT_SDGEO "\tNA" if $okSdGeo;
							}
						}
						else {
							my $transValue=$transform{$dataTransform}->($usedProtRatio{$quantif});
							print MAT_PCA "\t$transValue" if $pcaID;
							print MAT_CLUSTER "\t$transValue" if $clusteringID;
							if ($action eq 'export') {
								print MAT_EXPORT_QUANTIF "\t$transValue";
								print MAT_INFO_EXPORT "\tVal" if ($exportRprocessMatrices || $imputeData);
								$countTrueValues++;
								unless ($isNormalized) {
									if ($okPvalue) {
										if ($quantifValues{$modProtID}{$quantif}{'P_VALUE'}) {print MAT_EXPORT_PVAL "\t$quantifValues{$modProtID}{$quantif}{P_VALUE}";}
										else {print MAT_EXPORT_PVAL "\tNA";}
									}
									if ($okSdGeo) {
										if ($quantifValues{$modProtID}{$quantif}{'SD_GEO'}) {print MAT_EXPORT_SDGEO "\t".$quantifValues{$modProtID}{$quantif}{SD_GEO}*100/abs($transform{LOG2}->($usedProtRatio{$quantif}));}
										else {print MAT_EXPORT_SDGEO "\tNA";}
									}
								}
							}
						}
					}
					else { # Non-ratio quantifs
						my $transValue=$transform{$dataTransform}->($usedProtRatio{$quantif});
						print MAT_PCA "\t$transValue" if $pcaID;
						print MAT_CLUSTER "\t$transValue" if $clusteringID;
						if ($action eq 'export') {
							print MAT_EXPORT_QUANTIF "\t$transValue";
							print MAT_INFO_EXPORT "\tVal" if $imputeData;
						}
					}
					print MAT_EXPORT_PEP "\t$quantifValues{$modProtID}{$quantif}{$numPepType}" if ($action eq 'export' && !$isNormalized);
				}
				else { # missing value
					if ($action eq 'explorAna') {
						print MAT_PCA "\tNA" if $pcaID;
						print MAT_CLUSTER "\tNA" if $clusteringID;
					}
					else {
						$countMissingvalues++;
						print MAT_EXPORT_QUANTIF "\tNA";
						if ($exportRprocessMatrices || $imputeData) {
							if ($isInfiniteRatios{$quantif}) {
								print MAT_INFO_EXPORT "\t$isInfiniteRatios{$quantif}";
							}
							else {
								print MAT_INFO_EXPORT "\tMV"; #MV Missing Value
							}
						}
						if ($quantifFam eq 'RATIO' && !$isNormalized) {
							print MAT_EXPORT_PVAL "\tNA" if $okPvalue;
							print MAT_EXPORT_SDGEO "\tNA" if $okSdGeo;
						}
					}
					print MAT_EXPORT_PEP "\tNA" if ($action eq 'export' && !$isNormalized);
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
				print MAT_INFO_EXPORT "\n" if ($exportRprocessMatrices || $imputeData);
				if (!$isNormalized) {
					print MAT_EXPORT_PEP "\n";
					if ($quantifFam eq 'RATIO') {
						print MAT_EXPORT_PVAL "\n" if $okPvalue;
						print MAT_EXPORT_SDGEO "\n" if $okSdGeo;
					}
				}
				#unless ($usedProteinInfo{$proteinID}) {
					my ($gene,$synStrg)=("NO_GENE:$proteinInfo{$proteinID}[1]",''); # default
					if ($proteinInfo{$proteinID}[5][0]) { # gene
						$gene=shift @{$proteinInfo{$proteinID}[5]};
						$synStrg=join(',',@{$proteinInfo{$proteinID}[5]});
					}
					$proteinInfo{$proteinID}[2]=~s/\t/ /g; # just to be safe
					my ($ensemblGene,$entrezGene)=($extraAnnotation{$proteinID})? @{$extraAnnotation{$proteinID}} : ("NO_GENE:$proteinInfo{$proteinID}[1]","NO_GENE:$proteinInfo{$proteinID}[1]");
					print ANNOT "$protEntryLabel{$modProtID}\t$proteinInfo{$proteinID}[1]\t$ensemblGene\t$entrezGene\t$gene\t$synStrg\t$proteinInfo{$proteinID}[2]\n"; # [1]=IDENTIFIER // Now 1 line per entry in quantif table (PP 12/04/20)
					#$usedProteinInfo{$proteinID}=1;
				#}
			}
		}
		if ($action eq 'explorAna') {
			close MAT_CLUSTER if ($clusteringID);
			close MAT_PCA if ($pcaID);
		}
		else { # export
			close MAT_EXPORT_QUANTIF;
			close MAT_INFO_EXPORT if ($exportRprocessMatrices || $imputeData);
			close ANNOT;
			if (!$isNormalized) {
				close MAT_EXPORT_PEP;
				if ($quantifFam eq 'RATIO') {
					close MAT_EXPORT_PVAL if $okPvalue;
					close MAT_EXPORT_SDGEO if $okSdGeo;
				}
			}
		}

		###>No proteins in matrix<###
		unless ($nbAllProt) {
			my ($buttonLabel,$onClickAction);
			if ($action eq 'explorAna') {
				my $sthUpdateStatus = $dbh->prepare("UPDATE EXPLORANALYSIS SET STATUS=-2 WHERE ID_EXPLORANALYSIS = ?"); # set to error status
				$sthUpdateStatus->execute($clusteringID) if $clusteringID;
				$sthUpdateStatus->execute($pcaID) if $pcaID;
				$dbh->commit;
	
				#my $explorID=($clusteringID)? $clusteringID : $pcaID;
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
| if $action eq 'export';
			exit;
		}
	
		$dbh->disconnect;
		
		if ($action eq 'explorAna') {
			##my $explorID=$pcaID || $clusteringID;
	
			###>Precrossing data (Imputing +/-inf & missing values, filtering top sd prot)<###
			my $explorDIR = "$promsPath{tmp}/exploratory_analysis/$explorID";
			system "./launchExploratoryAnalyses.pl $explorID prepare &"; # 2> $errorPrepareFile &"; # Runs prepareExplorAna.R to preprocess data matrix
	
			##>Waiting for R to finish<##
			my $wait = 1;
			my $errorTxt = '';
			my $count1 = 0;
			my $count2 = 0;
			while ($wait == 1) {
				sleep 15;
				##print '.';
				#if (-s $errorPrepareFile) {
				#	$errorTxt=`cat $errorPrepareFile`;
				#	$wait=0;
				#	last;
				#}
				if (!-e "$explorDIR/prepareExplorAna.Rout") {
					if ($count1 > 40) { # <=> 10 min
						$errorTxt='No R output file found';
						$wait=0;
						last;
					}
					$count1++;
					next;
				}
				my $Rprocess = `grep -c '> proc.time()' $explorDIR/prepareExplorAna.Rout`;
				chomp $Rprocess;
				if ($Rprocess) { # finished normally
					$wait = 0;
					last;
				}
				else {
					$Rprocess = `grep -c '^Execution halted' $explorDIR/prepareExplorAna.Rout`;
					chomp $Rprocess;
					if ($Rprocess) {
						$errorTxt='Execution halted from R';
						$wait = 0;
						last;
					}
				}
				$count2++;
				if ($count2 > 480) { # * sleep 15sec <=> 2h
					$errorTxt='Process duration has exceeded 60 min.';
					$wait=0;
					last;
				}
			}
			if ($errorTxt) {
				
				#>Set status to error code -2
				$dbh = &promsConfig::dbConnect; # reconnect after R
				$dbh->do("UPDATE EXPLORANALYSIS SET STATUS=-2 WHERE ID_EXPLORANALYSIS=$pcaID") if $pcaID;
				$dbh->do("UPDATE EXPLORANALYSIS SET STATUS=-2 WHERE ID_EXPLORANALYSIS=$clusteringID") if $clusteringID;
				$dbh->commit;
				$dbh->disconnect;
			
##			print qq
##|</FONT>
##<BR><BR><BR><BR><BR>
##<FONT class="title2" color="#DD0000">***ERROR: Data preprocessing failed ($errorTxt)!***</FONT>
##<BR><INPUT type="button" class="title2" value="Continue" onclick="top.promsFrame.selectedAction='summary'; parent.itemFrame.location='$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=experiment:$experimentID&branchID=exploranalysis:$explorID&VIEW=explorAna'">
##<BR><BR><BR>
##</CENTER>
##</BODY>
##</HTML>
##|;
				exit;
			}
			#unlink $errorPrepareFile; # must be empty
			
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
			open(NA_PCA,">$promsPath{tmp}/exploratory_analysis/$pcaID/missingValues.txt") if $pcaID;
			open(NA_CLUSTER,">$promsPath{tmp}/exploratory_analysis/$clusteringID/missingValues.txt") if $clusteringID;
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
			$dbh=&promsConfig::dbConnect('no_user'); # reconnect after R
			my $sthUpdateFilterList=$dbh->prepare("UPDATE EXPLORANALYSIS SET FILTER_LIST=CONCAT(FILTER_LIST,?) WHERE ID_EXPLORANALYSIS=?");
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
					copy("$promsPath{tmp}/exploratory_analysis/$pcaID/valueDistribution.png","$promsPath{tmp}/exploratory_analysis/$clusteringID/valueDistribution.png") if -e "$promsPath{tmp}/exploratory_analysis/$pcaID/valueDistribution.png";
				}
			}
			elsif ($clusteringID) { # no PCA
				#move("$promsPath{tmp}/exploratory_analysis/$clusteringID/matrix.txt","$promsPath{tmp}/exploratory_analysis/$clusteringID/matrixBefore.txt");
				unlink "$promsPath{tmp}/exploratory_analysis/$clusteringID/matrix.txt";
				move("$promsPath{tmp}/exploratory_analysis/$clusteringID/matrixProcessed.txt","$promsPath{tmp}/exploratory_analysis/$clusteringID/matrix.txt");
			}
		}
		
	} # end of !$referenceExplorID

	if ($action eq 'explorAna') {
		###>Launch Analyses as unlinked process (&)<###
		if ($pcaID) {
			system "./launchExploratoryAnalyses.pl $pcaID PCA $projectID &";
		}
		if ($clusteringID) {
			system "./launchExploratoryAnalyses.pl $clusteringID cluster $projectID $metric $method $itemMetric &";
		}
	}
	else { # export
		my $errorTxt = '';
		my $exportDIR = "$promsPath{tmp}/export";

		if ($exportRprocessMatrices || $imputeData) {
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
				sleep 15;
				print '.';
				if (!-e "$exportPath/analysisModifQuantif.Rout") {
					if ($count1 > 40) { # <=> 10 min
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
					last;
				}
				else {
					$Rprocess = `grep -c '^Execution halted' $exportPath/analysisModifQuantif.Rout`;
					chomp $Rprocess;
					if ($Rprocess) {
						$errorTxt='Execution halted from R';
						$wait = 0;
						last;
					}
				}
				$count2++;
				if ($count2 > 240) { # *sleep 15sec <=> 60 min
					$errorTxt='Process duration has exceeded 60 min.';
					$wait=0;
					last;
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
		unlink glob '*.png';
		chdir $promsPath{cgi};
		rmdir $exportPath;
		
		#rmtree($exportPath); # delete original dir <==== Does not work 04/01/16

		##>Link to download
		print qq
| Done.<BR><BR><BR>
<INPUT type="button" class="title2" value="Download Dataset" onclick="window.location='$promsPath{tmp_html}/scratch/export/$archiveName'"/>
<SCRIPT type="text/javascript">
document.getElementById('waitSPAN').style.display='none';
top.promsFrame.selectedAction = 'summary';
</SCRIPT>
</BODY>
</HTML>
|;
	}
    exit;
}

#<<<--------------------------------END OF FORM SUBMISSION-----------------------------------#


#### EXPLORATORY ANALYSIS / EXPORT FORM ####
my $sthSelectList = $dbh->prepare("SELECT CL.ID_CLASSIFICATION,CL.NAME,CA.ID_CATEGORY,CA.NAME,CA.LIST_TYPE FROM CLASSIFICATION CL
					    INNER JOIN CATEGORY CA on CL.ID_CLASSIFICATION = CA.ID_CLASSIFICATION
					    INNER JOIN CATEGORY_PROTEIN CP on CA.ID_CATEGORY = CP.ID_CATEGORY
					    WHERE CL.ID_PROJECT = $projectID ORDER BY CA.DISPLAY_POS");
my $sthGetAllDesign = $dbh->prepare("SELECT ID_DESIGN FROM DESIGN WHERE ID_EXPERIMENT = $experimentID");
my $sthSelModif=$dbh->prepare("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DISPLAY_CODE,DISPLAY_COLOR FROM MODIFICATION WHERE ID_MODIFICATION=?");
my $sthGetExpCondName = $dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
my $sthGetAllModifs = $dbh->prepare("SELECT CONCAT(COALESCE(Q.ID_MODIFICATION,''),COALESCE(GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ','),'')),Q.QUANTIF_ANNOT,Q.ID_QUANTIFICATION
									FROM QUANTIFICATION Q
									LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
									WHERE ID_DESIGN=?
									GROUP BY Q.ID_QUANTIFICATION");

#print header(-charset=>'utf-8', -'content-encoding' => 'no'); #debug
my (%listClassCat, %listClassification, %listCategory, %modificationNotNorm, %modificationName, %freeResidues);

##Build the restrict/exclude menu
$sthSelectList->execute;
while (my ($classificationID, $classificationName, $categoryID, $categoryName, $listType) = $sthSelectList -> fetchrow_array) {
	#$listProtein{$proteinID} = $protAlias;
	#print "list:$proteinID<br>";
	next if ($listCategory{$categoryID});
	push @{$listClassCat{$classificationID}}, $categoryID;
	$listClassification{$classificationID} = $classificationName;
	@{$listCategory{$categoryID}} = ($categoryName,$listType);
}
$sthSelectList -> finish;

##Build the Focus menu (protein/sites and normalized)
my $hasFreeResQuantifs=0;
$sthGetAllDesign -> execute;
while (my ($designID) = $sthGetAllDesign -> fetchrow_array) {
    $sthGetAllModifs -> execute($designID);
    while (my ($multiModifStrg,$quantifAnnot) = $sthGetAllModifs -> fetchrow_array) {
		if ($multiModifStrg) {
			my $notNorm=1 unless $quantifAnnot=~/::INTRA_PROT_REF=/;
			foreach my $modifID (split(',',$multiModifStrg)) {
				if ($modifID==-1) {
					$hasFreeResQuantifs++;
					my ($context)=$quantifAnnot=~/::SITE_CONTEXT=#-1&([^:]+)/; # X-X-c/e/t-X{n,m}
					my ($targetRes)=$context=~/([a-z]|\/)/; # c/e/t cannot be x
					$targetRes=~s/\///g; # cet
					foreach my $res (split(//,uc($targetRes))) { # CET
						$freeResidues{$res}=$res;
					}
					next;
				}
				$modificationName{$modifID}=1;
				$modificationNotNorm{$modifID}=1 if $notNorm;
			}
		}
    }
}
foreach my $modID (keys %modificationName) {
	$sthSelModif->execute($modID);
	my ($psiName,$interName,$synName,$code,$color)=$sthSelModif->fetchrow_array;
	$modificationName{$modID}=$psiName || $interName || $synName;
	$modificationName{$modID}=~s/^##//; $modificationName{$modID}=~s/##.*$//;
	$modificationName{$modID}.=" [<FONT color='#$color'>$code</FONT>]";
}

$sthSelModif->finish;
$sthGetAllModifs->finish;
$sthGetAllDesign->finish;

###>Free residues<###
# if ($hasFreeResQuantifs) {
# 	my ($freeResSpecif)=$dbh->selectrow_array("SELECT SPECIFICITY FROM MODIFICATION WHERE ID_MODIFICATION=-1");
# 	foreach my $res (split(',',$freeResSpecif)) {
# 		$freeResidues{$res}=($res eq '-')? 'Protein N-term' : ($res eq '+')? 'Protein C-term' : $res;
# 	}
# $freeResidues{'C'}='C'; # TEMP!!!!!!
# }

####>Fetch ExplorAna for this Exp<####
my %allExplorAna;
if ($action eq 'explorAna') {
	my $sthAllExp=$dbh->prepare("SELECT E.ID_EXPLORANALYSIS,NAME,ANA_TYPE,COUNT(EQ.ID_QUANTIFICATION) FROM EXPLORANALYSIS E,EXPLORANA_QUANTIF EQ WHERE E.ID_EXPLORANALYSIS=EQ.ID_EXPLORANALYSIS AND  ID_EXPERIMENT=$experimentID GROUP BY E.ID_EXPLORANALYSIS");
	$sthAllExp->execute;
	while (my ($eaID,$eaName,$anaType,$numQuantif)=$sthAllExp->fetchrow_array) {
		@{$allExplorAna{$anaType}{$eaID}}=($eaName,$numQuantif);
	}
	$sthAllExp->finish;
}

$dbh->disconnect;

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
<STYLE type="text/css">
UL.ulList {
	padding:0px 5px 0px 5px;
	margin: 0px 5px 0px 5px;
	max-height:150px;
	list-style:none;
}
UL.ulList LI {margin:0; padding:0px 10px 0px 10px;}
</STYLE>
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

function updateDatasetSelection(eaID) {
	var dispStatus=(eaID=='0')? '' : 'none';
	var manualTRs=document.getElementsByClassName('manual');
	for (let i=0; i<manualTRs.length; i++) {
		manualTRs[i].style.display=dispStatus;
	}
}

function updateDesignTableDisplay(myImg,desTableId) {
	const desTable=document.getElementById(desTableId);
	if (desTable.style.display==='none') {
		desTable.style.display='';
		myImg.src='$promsPath{images}/minus1.gif';
	}
	else {
		desTable.style.display='none';
		myImg.src='$promsPath{images}/plus.gif';
	}
}
var numTestChecked=0;
var autoExtendSelection = 'project';
function updateQuantifStatus(chkStatus,datasetIdx,refID,fromUser) {

	if (chkStatus) {numTestChecked++;} else {numTestChecked--;}
	document.getElementById('chkQuantifs').innerHTML=numTestChecked;
	//Normalization dataset
	if (document.expAnaForm.focus.value=='2') {
		document.getElementById('ref_quantif_'+refID).disabled=(chkStatus)? false : true;
		if (chkStatus == false) {
			document.getElementById('ref_quantif_'+refID).value="";
		}
	}|;
	if ($action eq "export") {
		print qq|
	if (document.getElementById('anovaSel').value == "group"){
		//console.log('group_quantif_'+refID);
		document.getElementById('group_quantif_'+refID).disabled=(chkStatus)? false : true;
		if (chkStatus == false){  // && document.expAnaForm.focus.value.match('1\|2')
			document.getElementById('span_quantif_'+datasetIdx).innerHTML= "<input type = text name=group_quantif id=group_quantif_"+refID+" placeholder = " + "'group for Anova'" + " disabled>";
		}
	}|;
	}
	print qq|
	//Auto-extend selection?
	if (fromUser && document.getElementById('autoExtend').checked) {
		var testQuantifs=document.expAnaForm.test_quantif;
		//var currentChkBox = testQuantifs[datasetIdx];
		//var srcTable = currentChkBox.parentNode.parentNode.parentNode.parentNode;
		const designID=testQuantifs[datasetIdx].dataset.design;
		const [quantiID,tgtPos]=testQuantifs[datasetIdx].value.split('_');
		for (let i=datasetIdx+1; i < testQuantifs.length; i++) {
			if (autoExtendSelection == "design") {
				//var targetTable = testQuantifs[i].parentNode.parentNode.parentNode.parentNode;
				//if (targetTable != srcTable) { break; }
				if (testQuantifs[datasetIdx].dataset.design !== designID) {break;}
			}
			else if (autoExtendSelection == "quanti") {
				const [qID,tPos]=testQuantifs[i].value.split('_');
				if (qID !== quantiID) {break;}
			}
			if (testQuantifs[i].checked == chkStatus \|\| testQuantifs[i].disabled) {continue;}
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
		for (let i=0; i < testQuantifs.length; i++) {
			if (testQuantifs[i].checked == false) {continue;}
			testQuantifs[i].checked=false;
			updateQuantifStatus(false,i,testQuantifs[i].value,false);
		}
		if (refQuantifs) { // normalization
			for (let i=0; i < refQuantifs.length; i++) {
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
			for (let i=0; i<quantifList.length; i++) {
				var quantifs=quantifList[i].split('%'); // normalization ?
				quantifObject[quantifs[0]]=(quantifs[1])? quantifs[1] : '=';
			}
			for (let i=0; i < testQuantifs.length; i++) {
				var refQ=quantifObject[testQuantifs[i].value];
				if (refQ && !testQuantifs[i].checked) {
					testQuantifs[i].checked=true;
					updateQuantifStatus(true,i,testQuantifs[i].value,false);
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
	if (myForm.datasetSelection && myForm.datasetSelection.value != '0') { // nothing else to check
		document.getElementById('expAnaParams').style.display='none';
		return true;
	}
|;
}
print qq
|	var quantifFam=myForm.quantifFam.value;
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
			alert('Choose a valid (adj.) p-value (]0-1])');
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
	if (!myForm.anovaSel.disabled && myForm.anovaSel.value == "sample" && numTestChecked == 1) {
		alert('Select at least 2 quantification to be exported');
		return false;
	}
|;
}
print qq
|	var testQuantifs=myForm.test_quantif;
	var refQuantifs=myForm.ref_quantif; // may be null
	for (var i=0; i<testQuantifs.length; i++) {
		if (testQuantifs[i].checked && myForm.focus.value=='2' && !refQuantifs[i].value) { // normalized modif-proteome
			alert('Missing normalization quantification');
			return false;
		}
	}|;
	if ($action eq "export") {
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
function displayModifications(focusVal) {
	var [modifDisp,toNormModifDisp]=(focusVal==1)? ['','none'] : (focusVal==2)? ['none',''] : ['none','none'];
	document.getElementById('modifUL').style.display=modifDisp;
	document.getElementById('toNormModifUL').style.display=toNormModifDisp;
	if (focusVal==0 \|\| focusVal=='peptide') {
		ajaxDisplayData('focus',focusVal);
	}
	else {
		processModifSelection(focusVal,true);
	}
}
function processModifSelection(focusVal,doAjax) {
	var modifType=(focusVal==1)? document.expAnaForm.stdModif : document.expAnaForm.toNormModif;
	var selModifs=[];
	if (modifType.length) {
		for (let i=0; i<modifType.length; i++) {
			if (modifType[i].checked) selModifs.push(modifType[i].value);
		}
	}
	else if (modifType.checked) {selModifs=[modifType.value];}
	if (selModifs.length) {
		if (doAjax) {ajaxDisplayData('focus',selModifs.join(','));}
		else {return true;}
	}
	else {
		document.getElementById('listQuantif').innerHTML = '<B>-= Select at least 1 PTM =-</B>';
		if (!doAjax)  {return false;}
	}
}
function ajaxDisplayData(srcType,srcValue) {
	if (srcType=='focus') { // focus was changed
		if (srcValue=='peptide') {
			window.location="$promsPath{cgi}/startExploratoryAnalysisPeptide.cgi?ID=$experimentID";
		}
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
		else if (selQuantifFam.match('PROT_ABUNDANCE:')) {
			ratioVis='none';
			pepTypeOpts=[['NUM_PEP_USED','$pepTypeDesc{NUM_PEP_USED}'],['DIST_PEP_USED','$pepTypeDesc{DIST_PEP_USED}']];
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
		document.expAnaForm.dataTransform.options.selectedIndex=dataTransIdx;

		var pepTypeSEL=document.expAnaForm.pepType;
		pepTypeSEL.options.length=0;
		for (var i=0; i<pepTypeOpts.length; i++) {
			pepTypeSEL.options[i]=new Option(pepTypeOpts[i][1],pepTypeOpts[i][0]);
		}
	}
	var listDiv=document.getElementById('listQuantif');
	var focusVal=document.getElementById('focus').value;
	if (focusVal==1 \|\| focusVal==2) {
		if (!processModifSelection(focusVal,false)) return;
	}
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

    XHR.open("GET","./startExploratoryAnalysis.cgi?AJAX=displaySelect&FOCUS="+document.expAnaForm.focus.value+"&MODIF="+selModID+"&QUANTIF_FAM="+selQuantifFam+"&ID=$experimentID",true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState == 4 && XHR.responseText) {
			autoExtendSelection = 'project';
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
	print qq
|function updateTopNParams(isChk) {
	document.getElementById('topNvar').disabled=!isChk;
}
|;
}
print qq
|function checkImpute(type,chkStatus) {
	if (chkStatus) {
		 document.getElementById('infRatios').value=-1;
	}
	
}
function updateImpute(valStatus) {
	var disab=(valStatus==0)? true : false;
	if (document.getElementById('imputeNA')) {document.getElementById('imputeNA').disabled=disab;} // undef for explorAna

	var mvSel = document.getElementById("missingValues");
	var selectedId = mvSel.options[mvSel.selectedIndex].id;
	if(selectedId && selectedId == "optMissingX") {
		document.getElementById("spanMissingX").style.display =  "";
	} else {
		document.getElementById("spanMissingX").style.display =  "none";
	}
}
/*
function changeFeature(val) {
	if (val == "protQuant") {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/startExploratoryAnalysis.cgi?ID=$experimentID";
	}
	else {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/startExploratoryAnalysisPeptide.cgi?ID=$experimentID";
	}
}
*/
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">

<CENTER>
<FONT class="title">$title</FONT><BR><BR>
<DIV id="expAnaParams">
<FORM name="expAnaForm" method="POST" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="ID" value="$experimentID">
<INPUT type="hidden" name="ACT" value="$action">
|;
if ($action eq "export") {
	print qq
|<INPUT type="hidden" name="GROUP_LIST" value="">
<INPUT type="hidden" name="GROUP_ANNOT" value="">
|;
}
print qq
|<TABLE border="0" bgcolor="$darkColor">
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
	if (scalar keys %allExplorAna) {
		print qq
|	<TR>
		<TH colspan="2" align="right" valign="top" nowrap>Dataset selection :</TH><TD bgcolor="$lightColor" nowrap>
			<SELECT name="datasetSelection" onchange="updateDatasetSelection(this.value)"><OPTION value="0">Manual</OPTION><OPTGROUP label="Used Dataset from previous Analysis:">
|;
		foreach my $anaType (sort{lc($a) cmp lc($b)} keys %allExplorAna) {
			my $typeLabel=($anaType=~/^clus/i)? 'Clustering' : 'PCA';
			print "<OPTGROUP label=\"    -$typeLabel:\">\n";
			foreach my $eaID (sort{lc($allExplorAna{$anaType}{$a}[0]) cmp lc($allExplorAna{$anaType}{$b}[0])} keys %{$allExplorAna{$anaType}}) {
				print "<OPTION value=\"$eaID\">$allExplorAna{$anaType}{$eaID}[0] (x$allExplorAna{$anaType}{$eaID}[1] quantifications)</OPTION>\n";
			}
			print "</OPTGROUP>\n";
		}
		print qq
|</OPTGROUP></SELECT>&nbsp;</TD>
	</TR>
|;
	}
}
my $disabPepFocus=(1 || $action eq "export")? " disabled" : ""; # Disabled for now (PP 14/05/19)	
print qq
|	<TR class="manual"><TH colspan="2" align="right" valign="top" nowrap>Focus :</TH><TD bgcolor="$lightColor" style="min-width:600px" nowrap><SELECT name="focus" id="focus" onchange="displayModifications(this.value)">
	<OPTION value="peptide"$disabPepFocus>Peptides</OPTION>
	<OPTION value="0" selected>Proteins</OPTION>
	
|;
if (scalar keys %modificationName) {
	print "\t<OPTION value=\"1\">PTM-sites / Free residues:</OPTION>\n";
	# if (scalar keys %modificationNotNorm) { DISABLED TO SIMPLIFY (PP 02/12/20)
	# 	print "\t<OPTION value=\"2\">PTM-sites to normalize:</OPTION>\n";
	# }
}
print qq
|	</SELECT>
	<SPAN>&nbsp;&nbsp;<B>Quantification type:</B><SELECT name="quantifFam" id="QF" onchange="ajaxDisplayData('quantifFam',this.value)">
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
	print "<OPTGROUP label=\"$proteinQuantifFamilies{NAME}{$qFamily}:\">";
	if ($isParent) {
		foreach my $subFamily (@{$proteinQuantifFamilies{'MEMBERS'}{$qFamily}}) {
			foreach my $refMeas (@{$proteinQuantifFamilies{'MEASURES'}{$subFamily}}) {print "<OPTION value=\"$subFamily:$refMeas->[0]\">$refMeas->[1]</OPTION>";}
		}
	}
	else {print "<OPTION value=\"$qFamily\">$proteinQuantifFamilies{NAME}{$qFamily}</OPTION>\n";}
	print "</OPTGROUP>\n";
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
my $selImputeData='';
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
	#$selImputeData= ' disabled'; # commented because should be activated by default (PP)
}
print qq
|</SELECT></SPAN>
<DIV id="modifUL" style="display:none">
|;
if (scalar keys %modificationName) {
	print qq
|<B>&nbsp;PTM-sites:</B><BR>
<UL class="ulList">
|;
	foreach my $modID (sort{lc($modificationName{$a}) cmp lc($modificationName{$b})} keys %modificationName) {
		print "<LI><LABEL><INPUT type=\"checkbox\" name=\"stdModif\" value=\"$modID\" onclick=\"processModifSelection(1,true)\">$modificationName{$modID}-sites</LABEL></LI>\n";
	}
	print "</UL>\n";
}
if (scalar keys %modificationNotNorm) {
	print qq
|	<UL id="toNormModifUL" class="ulList" style="display:none">
|;
	foreach my $modID (sort{lc($modificationName{$a}) cmp lc($modificationName{$b})} keys %modificationNotNorm) {
		print "<LI><LABEL><INPUT type=\"checkbox\" name=\"toNormModif\" value=\"$modID\" onclick=\"processModifSelection(2,true)\">$modificationName{$modID}-sites</LABEL></LI>\n";
	}
	print "</UL>\n";
}
print qq
|<B>&nbsp;Free residues:</B><BR>
<UL class="ulList">
|;
	foreach my $res (sort keys %freeResidues) {
		print "<LI><LABEL><INPUT type=\"checkbox\" name=\"stdModif\" value=\"-1:$res\" onclick=\"processModifSelection(1,true)\">$freeResidues{$res}</LABEL></LI>\n";
	}
	print qq
|	</UL>
</DIV>	
		</TD>
    </TR>
	<TR class="manual">
	<TH colspan="2" align=right valign=top nowrap>Data transform :</TH>
	<TD nowrap bgcolor=$lightColor>&bull;<SELECT name="dataTransform"><OPTION value="">No change to</OPTION><OPTION value="LOG2">Log2-transform</OPTION><OPTION value="LOG10">Log10-transform</OPTION></SELECT> <B>values</B>
|;
print qq|&nbsp;&nbsp;&nbsp;&nbsp;&bull;<LABEL><INPUT type="checkbox" id="geneName" name="geneName" value="TRUE"><B>Use gene name as identifier</B></LABEL>| if ($action eq "export");
print qq
|	</TD>
	</TR>
    <TR class="manual">
	<TH colspan="2" align=right valign=top nowrap>Data filtering :</TH>
	<TD nowrap bgcolor=$lightColor>
	    <TABLE cellspacing=0>
		<TR id="RATIO_FC" style="display:none">
		    <TH align="right" nowrap>&nbsp;Abs. fold change &ge;</TH>
		    <TD><INPUT type="text" name="foldChange" value="1" size="2">&nbsp;<B>in at least&nbsp;</B><INPUT type="number" name="foldChangeOcc" value="1" style="width:50px">&nbsp;<B>quantification</B>&nbsp;</TD>
		</TR>
		<TR id="RATIO_INF" style="display:none">
		    <TH align="right">Infinite ratios:</TH>
			<TD><SELECT id="infRatios" name="infRatios"><OPTION value="-1">Treat as missing values</OPTION><OPTION value="0">None allowed</OPTION><OPTION value="5">Allow 5%</OPTION><OPTION value="10">Allow 10%</OPTION><OPTION value="25">Allow 25%</OPTION><OPTION value="34">Allow 34%</OPTION><OPTION value="50">Allow 50%</OPTION><OPTION value="100"$selAllInf>Allow all</OPTION></SELECT>&nbsp;per protein/site</TD>
			<!--<TD><SELECT id="infRatios" name="infRatios"><OPTION  name="infRatios" value="-1">Treat as missing values</OPTION><OPTION name="infRatios" value="0">None allowed</OPTION><OPTION name="infRatios" value="5">Allow 5%</OPTION><OPTION name="infRatios" value="10">Allow 10%</OPTION><OPTION name="infRatios" value="25">Allow 25%</OPTION><OPTION name="infRatios" value="34">Allow 34%</OPTION><OPTION name="infRatios" value="50">Allow 50%</OPTION><OPTION name="infRatios" value="100">Allow all</OPTION></SELECT>&nbsp;per protein</TD>-->

		</TR>
		<TR id="RATIO_PV" style="display:none">
		    <TH align="right">(Adj.) p-value &le;</TH>
		    <TD nowrap><INPUT type="text" name="pValue" value="1" size="5">&nbsp;<B>in at least&nbsp;</B><INPUT type="number" name="pValueOcc" value="1" style="width:50px">&nbsp;<B>quantification</B>&nbsp;<SMALL>(Does not apply to normalized ratios)</SMALL></TD>
		</TR>
		<TR><!-- id="RATIO_PEP" style="display:none" -->
		    <TH align="right">Peptides:</TH>
		    <TD><SELECT name="pepType"><!-- *Updated by JavaScript* --></SELECT> &ge; <INPUT type="number" name="numPep" value="$numSelPep" style="width:50px"></TD>
		</TR>
		<TR>
		    <TH align="right">Missing values:</TH>
			<TD><SELECT name="missingValues" id="missingValues" onchange="updateImpute(this.value)"><OPTION value="0">None allowed</OPTION><OPTION value="5">Allow 5%</OPTION><OPTION value="10">Allow 10%</OPTION><OPTION value="25">Allow 25%</OPTION><OPTION value="34"$sel34MissVal>Allow 34%</OPTION><OPTION value="50">Allow 50%</OPTION><OPTION id="optMissingX">Allow X%</OPTION><OPTION value="100"$selAllMissVal>Allow all</OPTION></SELECT><span id="spanMissingX" style="display:none;"> : <input name="inpMissingX" style="width:50px" id="inpMissingX" type="number" onfocusout="document.getElementById('optMissingX').value=this.value; document.getElementById('missingValues').value=this.value" length="3" max=100 placeholder="X" /> % </span>&nbsp;per protein
|;
print qq |&nbsp;&nbsp;<LABEL><INPUT type="checkbox" id="imputeNA" name="imputeNA" onchange="checkImpute('NA',this.checked)"$selImputeData><B>Impute missing values</B><LABEL>| if ($action eq "export");
print qq |</TD>
		</TR>
		<TR id="AMB" style="display:none">
		 <TH align="left" colspan=2><LABEL><INPUT type="checkbox" id="delocalize" name="delocalize" >Exclude ambiguous sites</LABEL></TH>
		<TR>
		<TR id="AGGREG" style="display:none">
		  <TH align="left" colspan=2><LABEL><INPUT type="checkbox" id="aggregate" name="aggregate" value="1" $strgKeep>Keep only the most "informative" site per protein</LABEL></TH>
		  <TD></TD>
		</TR>
		</TABLE>
	</TD>
	</TR>
	<TR id="PROT_FILTER" class="manual">
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
						<LABEL><INPUT type="checkbox" name="topNvarChk" id="topNvarChk" value="1" onchange="updateTopNParams(this.checked)">Keep only the<LABEL>:
					</TH>
					<TD nowrap>
						<INPUT type="number" name="topNvar" id="topNvar" value="200" style="width:80px" disabled> <B>most changing proteins</B>|;
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
|<TR class="manual"><TH colspan=2 align="right" valign="top">Data selection :</TH><TD nowrap bgcolor="$lightColor"><DIV id="listQuantif" style="max-height:300px;overflow:auto"><B>-= Select quantification type =-</B></DIV></TD></TR>
|;
if ($action eq 'export') {
	print qq
|<TR><TH colspan=2 align="right" valign="top">&nbsp;Quantification label format :</TH><TD nowrap bgcolor="$lightColor"><SELECT name="nameFormat">
<OPTION value="OPTIMAL">Auto-select distinct labels</OPTION>
<OPTION value="FULL">[Design if any.]Quantification name[.Ratio name if any]</OPTION>
<OPTION value="EXTENDED">Quantification name[.Ratio name if any]</OPTION>
<OPTION value="QUANTIF">Quantification name</OPTION>
<OPTION value="RATIO">Ratio name (Test.Reference states)</OPTION>
<OPTION value="TEST">Test state</OPTION>
</SELECT>&nbsp;<SPAN class="font11">(Non-alphanumerical characters will be changed to '.')</SPAN>&nbsp;</TD></TR>
|;
}
print qq
|	<TR>
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
#print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG

	my ($focus,$modifIDstrg0,$quantifFam,$fromMotif) = @_;
	my $normDataset=($focus==2)? 1 : 0;
	($modifIDstrg0,$quantifFam)=&promsMod::cleanParameters($modifIDstrg0,$quantifFam);
	my (%modifs,%freeRes,$modifIDstrg,$regExpModifIDstrg);
	if ($focus >= 1) {
		foreach my $modStrg (split(',',$modifIDstrg0)) {
			my ($modID,$res)=split(':',$modStrg); # clean Free Res meta data
			$modifs{$modID}=1;
			$freeRes{$res}=1 if ($modID==-1 && $res);
		}
		$modifIDstrg=join(',',keys %modifs);
		$regExpModifIDstrg=join('|',keys %modifs);
	}
	#if ($focus >= 1) {$regExpModifIDstrg=~s/,/\|/g;}
	
	###<Quantif visibility>###
	my $maxStatusStrg=($userInfo[1] eq 'bio' || $userInfo[1] eq 'manag')? ' = 1' : ($userInfo[1] eq 'mass')? ' IN (1,2)' : ' >= 1'; # filter for quantif visibility
	
	my (%testQuantifications,%normQuantifications,%designName,%quantifName,%fetchedCond,%experimentInfo,%noDesignQuantifications,%sampleInfo,%analysisInfo);
	####<RATIO or MaxQuant Intensities>####
	if ($quantifFam eq 'RATIO' || $quantifFam=~/^(MQ|PROT_RULER|PROT_ABUNDANCE):/) {
		my $quantifFam0=$quantifFam;
		($quantifFam,my $quantifCode)=split(':',$quantifFam); # MQ:MQ_INT -> MQ
		my @quantifMethIDs;
		foreach my $subFamily (@{$proteinQuantifFamilies{'MEMBERS'}{$quantifFam}}) {
			my ($quantifMethID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='$subFamily'");
			push @quantifMethIDs,$quantifMethID;
		}

		###<Modif info>###
		my %modifCodes;
		if ($focus >= 1) {
			# my $sthProjModif=$dbh->prepare("SELECT M.ID_MODIFICATION,DISPLAY_CODE,DISPLAY_COLOR FROM MODIFICATION M,PROJECT_MODIFICATION PM WHERE M.ID_MODIFICATION=PM.ID_MODIFICATION AND ID_PROJECT=$projectID");
			# $sthProjModif->execute;
			# while (my ($modID,$code,$color)=$sthProjModif->fetchrow_array) {
			# 	@{$modifCodes{$modID}}=("<FONT color='#$color'>$code</FONT>");
			# }
			# $sthProjModif->finish;
			
			my $sthModCode=$dbh->prepare("SELECT ID_MODIFICATION,DISPLAY_CODE,DISPLAY_COLOR FROM MODIFICATION WHERE DISPLAY_CODE IS NOT NULL");
			$sthModCode->execute;
			while (my ($modID,$code,$color)=$sthModCode->fetchrow_array) {
				@{$modifCodes{$modID}}=("<FONT color='#$color'>$code</FONT>");
				$modifCodes{$modID}[1]=\%freeRes if ($modID==-1 && scalar keys %freeRes);#<Free residues (with -1)
			}
			$sthModCode->finish;
		}
		
		###<Design quantifications>###
		my $sthGetAllDesign = $dbh -> prepare("SELECT ID_DESIGN,NAME FROM DESIGN WHERE ID_EXPERIMENT = ?");
		my $sthGetExpCondName = $dbh -> prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION = ?");
		my $sthGetTestInfo = ($focus==0)? 
										   $dbh->prepare("SELECT Q.ID_QUANTIFICATION,QUANTIF_ANNOT,NAME
																FROM QUANTIFICATION Q
																LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
																WHERE ID_DESIGN=? AND ID_QUANTIFICATION_METHOD IN (".join(',',@quantifMethIDs).")
																AND Q.ID_MODIFICATION IS NULL AND Q.STATUS $maxStatusStrg
																GROUP BY Q.ID_QUANTIFICATION
																HAVING GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ',') IS NULL")
										   
										: 
											# $dbh->prepare("SELECT Q.ID_QUANTIFICATION,QUANTIF_ANNOT,NAME,ID_MODIFICATION FROM QUANTIFICATION Q WHERE ID_DESIGN=? AND ID_QUANTIFICATION_METHOD IN (".join(',',@quantifMethIDs).") AND ID_MODIFICATION IN ($modifIDstrg) AND Q.STATUS $maxStatusStrg"),
											# $dbh->prepare("SELECT Q.ID_QUANTIFICATION,QUANTIF_ANNOT,NAME,GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ',') AS MULTI_MODIF
											# 					FROM QUANTIFICATION Q
											# 					LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
											# 					WHERE ID_DESIGN=? AND ID_QUANTIFICATION_METHOD IN (".join(',',@quantifMethIDs).")
											# 					AND Q.ID_MODIFICATION IS NULL AND Q.STATUS $maxStatusStrg
											# 					GROUP BY Q.ID_QUANTIFICATION
											# 					HAVING MULTI_MODIF REGEXP '(^|,)($regExpModifIDstrg)(,|\$)'")
											$dbh->prepare("SELECT Q.ID_QUANTIFICATION,QUANTIF_ANNOT,NAME,CONCAT(COALESCE(Q.ID_MODIFICATION,''),COALESCE(GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ','),'')) AS MULTI_MODIF
																FROM QUANTIFICATION Q
																LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
																WHERE ID_DESIGN=? AND ID_QUANTIFICATION_METHOD IN (".join(',',@quantifMethIDs).")
																AND Q.ID_MODIFICATION IS NULL AND Q.STATUS $maxStatusStrg
																GROUP BY Q.ID_QUANTIFICATION
																HAVING MULTI_MODIF REGEXP '(^|,)($regExpModifIDstrg)(,|\$)'");
										
		#my $sthChkTgPos=$dbh->prepare("SELECT 1 FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=? AND TARGET_POS=? LIMIT 1");
		
		my $quantifParamIdStrg; # Neede for fallback if mesaures are not listed in QUANTIF_ANNOT
		if ($quantifCode) {
			($quantifParamIdStrg)=$dbh->selectrow_array("SELECT GROUP_CONCAT(DISTINCT ID_QUANTIF_PARAMETER) FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD IN (".join(',',@quantifMethIDs).") AND CODE='$quantifCode'");
		}

		#foreach my $sthGetTestInfo (@sthGetAllTestInfo) {
			&fetchDesignQuantifications($experimentID,$quantifFam,\%testQuantifications,\%designName,\%quantifName,\%modifCodes,$sthGetAllDesign,$sthGetExpCondName,$sthGetTestInfo,$quantifCode,$quantifParamIdStrg); # $sthChkTgPos,$sthGetMeas
			$sthGetTestInfo->finish;
		#}
		#$sthGetMeas->finish if $sthGetMeas;

		if ($normDataset) {
			my $sthExp=$dbh->prepare("SELECT ID_EXPERIMENT,NAME,DISPLAY_POS FROM EXPERIMENT WHERE ID_PROJECT=$projectID"); # $projectID is global
			my $sthGetNormInfo = $dbh->prepare("SELECT Q.ID_QUANTIFICATION,QUANTIF_ANNOT,NAME
													FROM QUANTIFICATION Q
													LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
													WHERE ID_DESIGN=? AND ID_QUANTIFICATION_METHOD IN (".join(',',@quantifMethIDs).")
													AND Q.ID_MODIFICATION IS NULL AND Q.STATUS $maxStatusStrg
													GROUP BY Q.ID_QUANTIFICATION
													HAVING GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ',') IS NULL");
			$sthExp->execute;
			while (my ($expID,$expName,$dispPos)=$sthExp->fetchrow_array) {
				@{$experimentInfo{$expID}}=($expName,$dispPos);
				%{$normQuantifications{$expID}}=();
				&fetchDesignQuantifications($expID,$quantifFam,$normQuantifications{$expID},\%designName,\%quantifName,\%modifCodes,$sthGetAllDesign,$sthGetExpCondName,$sthGetNormInfo); # ,$sthChkTgPos
			}
			$sthExp->finish;
			$sthGetNormInfo->finish;
		}

		$sthGetAllDesign -> finish;
		$sthGetExpCondName->finish;
		# $sthChkTgPos->finish;

		###<Non-design quantifications>### (cannot be modif-quanitfs)
		#unless ($normDataset) { # skip non-design quantif if normalization
		if (!$normDataset && $quantifFam !~ /^(MQ|PROT_RULER|PROT_ABUNDANCE)$/) {
			my ($quantifMethID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='PROT_RATIO_PEP'");
			my $sthND=$dbh->prepare("SELECT S.ID_SAMPLE,S.NAME,S.DISPLAY_POS,A.ID_ANALYSIS,A.NAME,A.DISPLAY_POS,Q.ID_QUANTIFICATION,Q.NAME,QUANTIF_ANNOT
									FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ,ANALYSIS A,SAMPLE S
									WHERE Q.ID_QUANTIFICATION=AQ.ID_QUANTIFICATION AND AQ.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE
									AND Q.ID_QUANTIFICATION_METHOD=$quantifMethID AND S.ID_EXPERIMENT=$experimentID AND FOCUS='protein' AND ID_DESIGN IS NULL AND Q.STATUS $maxStatusStrg");
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
					my ($testCondID,$refCondID) = split(/\//,$ratio);
					push @{$noDesignQuantifications{$sampID}{$anaID}{$quantifID}},["$stateInfo{$testCondID}{NAME}/$stateInfo{$refCondID}{NAME}",$ratioPos];
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
		my $sthQ=$dbh->prepare("SELECT S.ID_SAMPLE,S.NAME,S.DISPLAY_POS,A.ID_ANALYSIS,A.NAME,A.DISPLAY_POS,Q.ID_QUANTIFICATION,Q.NAME FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ,ANALYSIS A,SAMPLE S
								WHERE Q.ID_QUANTIFICATION=AQ.ID_QUANTIFICATION AND AQ.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND Q.ID_QUANTIFICATION_METHOD=$quantifMethID AND S.ID_EXPERIMENT=$experimentID
								AND FOCUS='protein' AND Q.STATUS $maxStatusStrg");
		#my $sthQ2=($empaiQuantifParamID)? $dbh->prepare("SELECT 1 FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=? AND ID_QUANTIF_PARAMETER=$empaiQuantifParamID LIMIT 1") : 
		#								$dbh->prepare("SELECT 1 FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=? LIMIT 1");
		$sthQ->execute;
		while (my ($sampID,$sampName,$sampPos,$anaID,$anaName,$anaPos,$quantifID,$quantName)=$sthQ->fetchrow_array) {
			#$sthQ2->execute($quantifID);
			#my ($ok)=$sthQ2->fetchrow_array;
			my $dbhLite=&promsQuantif::dbConnectProteinQuantification($quantifID,$projectID);
			my $sthQ2=($empaiQuantifParamID)? $dbhLite->prepare("SELECT 1 FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIF_PARAMETER=$empaiQuantifParamID LIMIT 1") : 
											  $dbhLite->prepare("SELECT 1 FROM PROTEIN_QUANTIFICATION WHERE LIMIT 1");
			$sthQ2->execute;
			my ($ok)=$sthQ2->fetchrow_array;
			$sthQ2->finish;
			$dbhLite->disconnect;
			next if ($empaiQuantifParamID && !$ok);
			push @{$noDesignQuantifications{$sampID}{$anaID}{$quantifID}},['',0,$ok || 0];
			@{$analysisInfo{$anaID}}=($anaName,$anaPos);
			@{$sampleInfo{$sampID}}=($sampName,$sampPos);
			$quantifName{$quantifID}=$quantName || 'No name';
		}
		$sthQ->finish;
		#$sthQ2->finish if $sthQ2;
	}

	####<Fetch ExplorAna with same quantifFam>####
	my %allExplorAna;
	my $sthAllExp=$dbh->prepare("SELECT E.ID_EXPLORANALYSIS,NAME,ANA_TYPE,COUNT(EQ.ID_QUANTIFICATION) FROM EXPLORANALYSIS E,EXPLORANA_QUANTIF EQ WHERE E.ID_EXPLORANALYSIS=EQ.ID_EXPLORANALYSIS AND  ID_EXPERIMENT=$experimentID GROUP BY E.ID_EXPLORANALYSIS");
	$sthAllExp->execute;
	while (my ($eaID,$eaName,$anaType,$numQuantif)=$sthAllExp->fetchrow_array) {
		@{$allExplorAna{$anaType}{$eaID}}=($eaName,$numQuantif);
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

    print header(-type=>'text/plain',-charset => 'utf-8');
    warningsToBrowser(1);
	my %itemIcones=&promsConfig::getItemIcones;

	if (!scalar keys %testQuantifications && !scalar keys %noDesignQuantifications){
		print "<FONT class=\"title3\" color=\"#DD0000\">No quantifications found.</FONT>\n";
		exit;
	}


	if (!$fromMotif) {
		print qq|
<LABEL><INPUT type="checkbox" id="autoExtend"><B>Auto-extend selection within </B><SELECT onchange='autoExtendSelection = this.value;'><OPTION value='project' selected>Project</OPTION><OPTION value='design'>Design</OPTION><OPTION value='quanti'>Quantification</OPTION></SELECT></LABEL>
<TABLE cellspacing=0 width=100%>
<TR bgcolor="$darkColor"><TH>&nbsp;Quantifications (<SPAN id="chkQuantifs">0</SPAN> selected)&nbsp;<BR>
&nbsp;Copy selection from:<SELECT name="allExpAnaSEL" id="allExpAnaSEL" onchange="copySelectionFromExpAna(this.value)"><OPTION value="">-= Select =-</OPTION><OPTION value="-1">** Deselect all **</OPTION>
|;
		foreach my $anaType (sort{lc($a) cmp lc($b)} keys %allExplorAna) {
			my $typeLabel=($anaType=~/^clus/i)? 'Clustering' : 'PCA';
			print "<OPTGROUP label=\"$typeLabel:\">\n";
			foreach my $eaID (sort{lc($allExplorAna{$anaType}{$a}[0]) cmp lc($allExplorAna{$anaType}{$b}[0])} keys %{$allExplorAna{$anaType}}) {
				print "<OPTION value=\"$eaID\">$allExplorAna{$anaType}{$eaID}[0] (x$allExplorAna{$anaType}{$eaID}[1] quantifications)</OPTION>\n";
			}
			print "</OPTGROUP>\n";
		}
		print qq
|</SELECT>&nbsp;
</TH>
<TH><SPAN id="groupHeaderSpan" style="display:none">
		<TABLE cellspacing=0 cellpadding=0>
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
	my $datasetIdx=0;
	#my $ratioNameSep=($quantifFam=~/^(RATIO|MQ)$/)? ' : ' : '';
	my $ratioNameSep=($quantifFam=~/^(SIN|EMPAI)$/)? '' : ' : ';

	###<Design quantifs>###
	if (scalar keys %testQuantifications) {
		if (!$fromMotif) {
			print qq
|<TR><TH colspan="$colSpan" align="left"></TH></TR>
<TR bgcolor="$darkColor"><TH colspan="$colSpan" align="left">Design-based quantifications:</TH></TR>
|;
		}
		else {
			print qq|<SELECT name="quantif" id="quantif"><OPTION value="">-= Select =-</OPTION>|;
		}

		# my $prevDesignID=0;
		foreach my $designID (sort{&promsMod::sortSmart(lc($designName{$a}),lc($designName{$b}))} keys %testQuantifications) {
			if ($fromMotif) {
				print qq|<OPTGROUP label="$designName{$designID} :">|;
			}

			my $prevQuantifID=0;
			if (!$fromMotif) {
				#print qq |<TR style="background-color:white;height:.5em"><TD></TD><TD></TD></TR></TABLE>| if $prevDesignID;
				print qq 
|<TR><TH align="left" colspan=2 style="background-color:white"><IMG src="$promsPath{images}/minus1.gif" onclick="updateDesignTableDisplay(this,'table_des$designID')"><IMG src="$promsPath{images}/$itemIcones{design}">$designName{$designID}:</TH></TR>
<TR><TD colspan=2><TABLE id="table_des$designID" cellspacing=0 width=100%>
|;
			}

			my $bgColor=$lightColor;
			foreach my $quantifID (sort{&promsMod::sortSmart(lc($quantifName{$a}),lc($quantifName{$b}))} keys %{$testQuantifications{$designID}}) {
				foreach my $refRatio (@{$testQuantifications{$designID}{$quantifID}}) {
					# my $designNameStrg=($designID==$prevDesignID)? '<SPAN style="visibility:hidden">'.$designName{$designID}.' > </SPAN>' : $designName{$designID}.' > ';
					my $quantifNameStrg=($quantifID==$prevQuantifID)? '<SPAN style="visibility:hidden">'.$quantifName{$quantifID}.$ratioNameSep.'</SPAN>' : $quantifName{$quantifID}.$ratioNameSep;
					my ($ratioName,$rank,$hasData)=@{$refRatio};
					my $valueStrg=$quantifID.'_'.$rank;

					if ($fromMotif) {
						my ($disabStrg,$errorMsg)=($hasData)? ('','') : ('disabled',' [No data]');
						print qq|<OPTION value="$valueStrg" $disabStrg>&nbsp;&nbsp;$quantifNameStrg$ratioName$errorMsg</OPTION>|;
					}
					else {
						my ($chkBoxDisabStrg,$chkBoxVisStrg,$errorMsg)=($hasData)? ('','','') : ('disabled','style="visibility:hidden"',' [No data]');
						print qq|<TR class="list" bgcolor="$bgColor">|;
						#if (!$fromMotif) {
							print qq|<TD nowrap><INPUT type="checkbox" name="test_quantif" $chkBoxVisStrg value="$valueStrg" data-design="$designID" style="margin-left:30px" onclick="updateQuantifStatus(this.checked,$datasetIdx,'$valueStrg',true)" $chkBoxDisabStrg>$quantifNameStrg$ratioName$errorMsg</TD>|;
						#}
						#else{
						#	print qq|<TD nowrap><INPUT type="checkbox" name="$chkBoxName" $chkBoxVisStrg id="test_quantif_$valueStrg" value="$valueStrg" onclick="updateQuantifStatus(this.checked,'$valueStrg',true)">$designNameStrg$quantifNameStrg$ratioName</TD>|;
						#}
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
										my ($ratioNameRef,$rankRef,$hasData)=@{$refRatioRef};
										my $valueRefStrg=$quanID.'_'.$rankRef;
										my ($disabStrg,$errorMsg)=($hasData)? ('','') : ('disabled',' [No data]');
										print "<OPTION value=\"$valueRefStrg\" $disabStrg>&nbsp;&nbsp;$ratioNameRef$errorMsg</OPTION>\n";
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
					# $prevDesignID=$designID;
					$prevQuantifID=$quantifID;
				}
			}
			print "</TABLE></TD></TR>\n" if !$fromMotif;
		}
		if ($fromMotif) {
			print qq|</SELECT>|;
		}
	}

	###<Non-Design quantifs>###
	if (scalar keys %noDesignQuantifications) {
		print qq
		|<TR><TH colspan="$colSpan" align="left"></TH></TR>
		<TR bgcolor="$darkColor"><TH align="left" colspan="$colSpan">Internal-analysis quantifications:</TH></TR>
		|;
		my $bgColor=$lightColor;
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

						if (!$fromMotif) {
							print "<TD><SPAN id=\"span_quantif_$datasetIdx\"><INPUT type=\"text\" name=\"group_quantif\" id=\"group_quantif_$valueStrg\" placeholder=\"group for Anova\" disabled style=\"display:none\"></SPAN></TD>\n";
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
	my ($expID,$quantifFam,$refQuantifications,$refDesignName,$refQuantifName,$refModifCodes,$sthGetAllDesign,$sthGetExpCondName,$sthGetQuantifInfo,$quantifCode,$quantifParamIdStrg)=@_; # ,$sthChkTgPos,$sthGetMeas
	$quantifCode=$quantifFam unless $quantifCode; # eg; no $quantifCode for RATIO
	$sthGetAllDesign->execute($expID);
	while (my ($designID,$desName) = $sthGetAllDesign -> fetchrow_array) {
		$refDesignName->{$designID}=$desName;
		$sthGetQuantifInfo -> execute($designID);
		my %fetchedCond;
		while (my ($quantifID, $quantifAnnot, $quantName, $modifIdStrg) = $sthGetQuantifInfo -> fetchrow_array) {
			
			##<Quick check for available measures
			my (@measures,@usedFreeRes);
			my ($labelTypeStrg,@labelInfo) = split('::',$quantifAnnot);
			foreach my $infoStrg (@labelInfo) {
				my ($setting,$valueStrg)=split('=',$infoStrg);
				if ($quantifFam eq 'RATIO' && $setting eq 'RATIOS') {@measures=('RATIO');} #$numRatios=scalar (split(';',$valueStrg)); 
				#elsif ($setting eq 'STATES' && $valueStrg=~/\d,#/) {$numStates=scalar (split(';',$valueStrg));}
				elsif ($setting=~/^(ABUND_MEASURES|OUTPUT_MEAS)/) { # Abundance or Proteomic Ruler
					foreach my $meas (split(';',$valueStrg)) {
						next if $meas=~/^aLFQ_/;
						$meas=~s/,.+//; # clean extra parameters 
						push @measures,$meas;
					}
				}
	# Problem with MaxQuant TODO: update all MQ quantif with ABUND_MEASURES
				elsif ($setting eq 'SITE_CONTEXT') {
					# @usedFreeRes=split(';',$valueStrg);
					my ($targetRes)=$valueStrg=~/([a-z]|\/)/; # c/e/t cannot be x
					$targetRes=~s/\///g; # cet
					@usedFreeRes=split(//,uc($targetRes));
				}
			}
			my $okQuantif=0;
			if (scalar @measures) {
				foreach my $measCode (@measures) {
					if ($measCode eq $quantifCode) {
						$okQuantif=1;
						last;
					}
				}
			}
			elsif ($quantifParamIdStrg) { # Fallback
				my $dbhLite=&promsQuantif::dbConnectProteinQuantification($quantifID,$projectID);
				if ($quantifParamIdStrg) {
					($okQuantif)=$dbhLite->selectrow_array("SELECT 1 FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIF_PARAMETER IN ($quantifParamIdStrg) LIMIT 1");
				}
				$dbhLite->disconnect;
			}
			next unless $okQuantif;

			if ($modifIdStrg) {
				if ($refModifCodes->{-1}[1] && scalar @usedFreeRes) { # Check free residues compatibility
					my $okFreeRes=0;
					foreach my $res (@usedFreeRes) {
						if ($refModifCodes->{-1}[1]->{$res}) {
							$okFreeRes=1;
							last;
						}
					}
					next unless $okFreeRes;
				}
				my $modifCodeStrg='';
				foreach my $modID (split(',',$modifIdStrg)) {
					$modifCodeStrg.='/' if $modifCodeStrg;
					my @freeResNames=map{($_ eq '-')? 'Prot. Nter' : ($_ eq '+')? 'Prot. Cter' : $_} @usedFreeRes;
					$modifCodeStrg.=($modID==-1)? 'Free '.join(',',@freeResNames) : $refModifCodes->{$modID}[0];
				}
				$quantName.=' ['.$modifCodeStrg.']';
			}
			#my ($labelTypeStrg,@labelInfo) = split('::',$quantifAnnot);
			#my $okQuantif=0;
			
			# my %usedTgtPos;
			# my $sthChkTgPos=$dbhLite->prepare("SELECT DISTINCT TARGET_POS FROM PROTEIN_QUANTIFICATION WHERE TARGET_POS IS NOT NULL");
			# $sthChkTgPos->execute;
			# while (my ($tgtPos)=$sthChkTgPos->fetchrow_array) {$usedTgtPos{$tgtPos}=1;}
			# $sthChkTgPos->finish;
			
			foreach my $infoStrg (@labelInfo) {
				my ($setting, $valueStrg)=split('=',$infoStrg);
				if ($quantifFam eq 'RATIO' && $setting eq 'RATIOS') {
					my $ratioPos=0;
					foreach my $ratio (split(";",$valueStrg)) {
						$ratioPos++;
						#$sthChkTgPos->execute($quantifID,$ratioPos);
						#my ($okTgPosData)=$sthChkTgPos->fetchrow_array;
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
						push @{$refQuantifications->{$designID}{$quantifID}},[$ratioItem,$ratioPos,1]; # $usedTgtPos{$ratioPos} || 0]; # $okTgPosData
					}
					last;
				}
				elsif ($quantifFam =~ /^(MQ|PROT_ABUNDANCE)$/ && $setting eq 'STATES') {
					$valueStrg=~s/#//g;
					my $statePos=0;
					foreach my $state (split(";",$valueStrg)) {
						$statePos++;
						#$sthChkTgPos->execute($quantifID,$statePos);
						#my ($okTgPosData)=$sthChkTgPos->fetchrow_array;
						my ($numBioRep,$quantiObsIDs,$expCondID)= split(',',$state);
						$sthGetExpCondName->execute($expCondID);
						my ($stateName)=$sthGetExpCondName->fetchrow_array;
						push @{$refQuantifications->{$designID}{$quantifID}},[$stateName,$statePos,1]; # $usedTgtPos{$statePos} || 0]; # $okTgPosData
					}
					last;
				}
				elsif ($quantifFam eq 'PROT_RULER' && $setting eq 'PARENT_Q') {
					$valueStrg=~s/#//g;
					foreach my $parent (split(";",$valueStrg)) {
						my ($parQuantID,$expCondID,$condPos,$parentPos)=split(/[,:]/,$parent);
						#$sthChkTgPos->execute($quantifID,$parentPos);
						#my ($okTgPosData)=$sthChkTgPos->fetchrow_array;
						$sthGetExpCondName->execute($expCondID);
						my ($stateName)=$sthGetExpCondName->fetchrow_array;
						push @{$refQuantifications->{$designID}{$quantifID}},[$stateName,$parentPos,1]; # $usedTgtPos{$parentPos} || 0]; # $okTgPosData
					}
					last;
				}
			}
			$refQuantifName->{$quantifID} = $quantName; # if $okQuantif;
			#$dbhLite->disconnect;
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
# 2.5.11 [UPDATE] Export: parameters.txt report quantification(s) used & adapted for SITE_CONTEXT format instead of FREE_RESIDUES (PP 19/03/21)
# 2.5.10 [ENHANCEMENT] Show/hide design content in &ajaxDisplayData (PP 12/03/21)
# 2.5.9 [ENHANCEMENT] Allow selection of custom missing values threshold (VS 04/01/20)
# 2.5.8 [ENHANCEMENT] Removed SQLite scan to speed up matched quantifications listing (PP 03/12/20)
# 2.5.7 [ENHANCEMENT] Also allow auto-extend within a quantification (PP 01/12/20)
# 2.5.6 [ENHANCEMENT] Allow auto-extend based on the whole project or within the design of the selected quantification (VS 10/09/20)
# 2.5.5 [BUGFIX] Fix quantifications filtering based on selected QUANTIFICATION_PARAMETER (VS 10/09/20)
# 2.5.4 [ENHANCEMENT] Uses quantification data from SQLite file (PP 03/08/20)
# 2.5.3 [BUGFIX] Fix lack of protein-level normalization in non-ratio quantifications (PP 21/04/20)
# 2.5.2 [FEATURE] Add Ensembl GeneID and Entrez GeneID to annotation.txt (PP 13/04/20)
# 2.5.1 [BUGFIX] annotation.txt also exported for non-ratio quantifs with same number of entries than data files even for site (PP 12/04/20)
# 2.5.0 [FEATURE] Option to reuse a full dataset from an existing Exploratory analysis (PP 03/03/20)
# 2.4.2 [ENHANCEMENT] Early fork to run as background process before fetching quantification data (25/02/20)
# 2.4.1 [ENHANCEMENT] No longer uses fork to launch background commands to minimize memory usage (PP 21/02/20)
# 2.4.0 [FEATURE] Handles myProMS Protein Abundance (PP 10/01/20)
# 2.3.2 [FEATURE] Improvement in column/row header formats in export mode & [BUGFIX] in density plot image file management (PP 18/11/19)
# 2.3.1 [ENHANCEMENT] Explor analysis status is set to failed (-2) on R data preprocessing error (PP 06/11/19)
# 2.3.0 [FEATURE] Handles Proteomic Ruler & [ENHANCEMENT] Improved filter for unusable quantifications (PP 29/10/19)
# 2.2.2 [FEATURE] Handles quantification visibility & [BUGFIX] in response to AJAX call from startMotifAnalysis.cgi (PP 27/09/19)
# 2.2.1 [ENHANCEMENT] Handles R generated density plot (count normalized) as imputation quality representation (VS 14/10/19)
# 2.2.0 [FEATURE] Compatible with multi-modif quantifications (PP 15/07/19)
# 2.1.10 Multiple bug fixes & simplified protein/peptide switch [peptide option disabled for now] & wait loop improvement (PP 15/05/19)
# 2.1.9 Add impute data for export (15/01/19)
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
