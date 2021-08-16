#!/usr/local/bin/perl -w

################################################################################
# showProtQuantification.cgi     2.14.4                                        #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Displays protein quantification data                                         #
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
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use promsQuantif;
use phosphoRS;
use strict;
use Spreadsheet::WriteExcel;
#use utf8; # Tells Perl that characters are UTF8. Necessary for Excel export to work with UTF-8 characters ...
use Encode qw(encode_utf8 decode_utf8); # ... Encode needed if hard UTF8 chars in code!!!
#use File::Path qw(rmtree); # remove_tree

# print header(-charset=>'utf-8'); warningsToBrowser(1); #DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my $projectAccess; # must be global
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my %labelingName=('FREE'=>'Label-free','SILAC'=>'SILAC','ITRAQ'=>'iTRAQ','TMT'=>'TMT');
#my $MAX_INF_RATIO_DB=1000; # true infinite value in database
#my $MIN_INF_RATIO_DB=0.001;
# Limit ratios for +/-infinity switch in volcano plot
my ($MIN_INF_RATIO,$MAX_INF_RATIO,$MAX_INF_RATIO_DB,$MIN_INF_RATIO_DB)=&promsQuantif::getExtremeRatios; # Max ratios allowed before switching to infinite
# Default: $MAX_INF_RATIO & $MIN_INF_RATIO overwritten if $fullRatioRange
my $log2=log(2);
my %normalizationNames=&promsQuantif::getQuantifNormalizationName;
my %biasCorrLevelNames=('ion'=>'peptide','prot'=>'protein','ion_prot'=>'peptide & protein',
						'uf_ion'=>'Peptide','uf_prot'=>'Protein','uf_ion_prot'=>'Peptide & protein');
my %pValueCorrections=('none'=>'None','fdr'=>'Benjamini-Hochberg (FDR)','bonferroni'=>'Bonferroni'); # for algo v3 ([Simple/Super]Ratio)
my %fdrMethods=('FDR-BH'=>'Benjamini-Hochberg',
					'FDR-ABH'=>'Benjamini-Hochberg (Adaptative)',
					'FWER-BONF'=>'Bonferroni',
					'Qvalue'=>'Qvalue (Storey et al.)'
				   ); # Very old algo only (Ratio)
my %ptmProbSoft=('PRS'=>'PhosphoRS','MQ'=>'MaxQuant','SPC'=>'Spectronaut','PTMRS'=>'PtmRS');
my %absCIsCodes=(  # Additional uncertainty measures (confidence intervals) when displaying absolute quantification
	'MOL_PERCENT'	=> ['CI_INF_MOL_PCT', 'CI_SUP_MOL_PCT'],
	'MOL'			=> ['CI_INF_MOL', 'CI_SUP_MOL'],
	'MOL_CONC'		=> ['CI_INF_MOL', 'CI_SUP_MOL'],
	'MASS_PERCENT'	=> ['CI_INF_MASS_PCT', 'CI_SUP_MASS_PCT'],
	'MASS'			=> ['CI_INF_MASS', 'CI_SUP_MASS'],
	'MASS_CONC'		=> ['CI_INF_MASS', 'CI_SUP_MASS']
);
my %dispAbsUnits=(  # Units to display for absolute quantifications
	'MOL_PERCENT'	=> "%",
	'MOL'			=> "mmol",
	'MOL_CONC'		=> "mmol/mL",
	'MASS_PERCENT'	=> "%",
	'MASS'			=> "ng",
	'MASS_CONC'		=> "µg/mL"
);

####################
####>Parameters<####
####################
my $call=param('CALL') || 'quanti'; # ana or quanti
my $action=param('ACT') || 'select';
my $analysisID=param('id_ana') || 0; # 0=> call=quanti
my $selQuantifID=param('id_quantif') || 0; # 0=> call=ana
($analysisID,$selQuantifID)=&promsMod::cleanNumericalParameters($analysisID,$selQuantifID);
my $view=param('view') || 'graph'; # 'list';
my $fullRatioRange=param('fullRatioRange') || 0;
if ($fullRatioRange) {$MAX_INF_RATIO=1000; $MIN_INF_RATIO=1/$MAX_INF_RATIO;}
my $allowInfRatio=param('allowInfRatio') || 0; # allow infinite ratios in List view even with a p-value threshold (only for ratio)
my $foldChangeType=param('foldChgType') || 'abs'; $foldChangeType='abs' if $view eq 'graph';
my $dispFoldChange=param('foldChange') || 2; $dispFoldChange=1 if ($dispFoldChange=~/[^\d\.]/ || $dispFoldChange < 1);
my $filterProtIDs=(param('filterProtIDs'))? param('filterProtIDs') : '';
my $foldChangeRaw=param('foldChgRaw') || 0; # MaxQuant ratio only (0->NORM or 1->RAW)
my $ratioParamCode=($foldChangeRaw)? 'RATIO_RAW' : 'RATIO'; # raw ratio only for Software MaxQuant
#my $dispPvalue=(param('pValue'))? param('pValue') : 0.05; # undef at start
my $dispPvalue=param('pValue') || 0.05; $dispPvalue=1 if ($dispPvalue=~/[^\d\.e-]/ || $dispPvalue > 1);
#my $dispStdDev=param('stdDev') || 0; $dispStdDev=0 if $dispStdDev < 0; # fraction
#my $dispStdDevFrac=$dispStdDev/100;
my $dispCV=param('coefVar') // 0; $dispCV=0 if $dispCV=~/[^\d\.]/; # includes negative numbers
my $dispCVFrac=$dispCV/100;
my $dispPepType=param('pepType') || 'all'; # all/dist/msms(Num|Pc)
#my $numPepCode=($dispPepType eq 'all')? 'NUM_PEP_USED' : ($dispPepType eq 'distinct')? 'DIST_PEP_USED' : ($dispPepType eq 'allRep')? 'NUM_PEP_USED:REPLIC' : ($dispPepType eq 'distinctRep')? 'DIST_PEP_USED:REPLIC' : $dispPepType;
my $numPepCode=($dispPepType eq 'all')? 'NUM_PEP_USED' : ($dispPepType eq 'distinct')? 'DIST_PEP_USED' : ($dispPepType =~ /^msms/)? 'NUM_TRUE_USED' : $dispPepType;
my $numPepScope=param('numPepScope') || 'ratio'; # ratio/replicate(One|Both|Test|Ref)/msms(One|Both|Test|Ref)
my ($dispNumPep,$updateNumPep)=(param('numPep'))? (param('numPep'),0) : (3,1); # Can be changed to 1 later if modif-quantif & $updateNumPep
my $dispNumReplic=param('numReplic') || 1;
my (%dispRatios,%dispStates,%condToState,%ratioCode,%quantifParamInfo,%quantifParamID2code,$designID,$ratioType); # globals also for some ajax calls (%dispRatios & %dispStates mutually exclusive)
if (param('dispRatios')) {foreach my $pos (param('dispRatios')) {$dispRatios{$pos}=[];}} # Positions of ratios displayed in results
elsif (param('dispStates') || (defined(param('dispStates')) && param('dispStates') eq "0")) {foreach my $pos (param('dispStates')) {$dispStates{$pos}=1;}} # Positions of states displayed in results
my $dispSort=param('sort') || 'peptide'; # for list only
$dispSort='peptide' if (param('dispRatios') && $dispSort=~/ratio_(-?\d+)/ && !$dispRatios{$1});
my $dispMeasure=param('dispMeasure') || ''; # Defined later according to quantif method
my $restrictListID=param('restrictList') || 0; # exclusion/selection List ID
my $listAction=param('listAction') || 'restrict'; # exclude/restrict
my $quantifType=''; # never initialized if called from internal quantif
my $protRulerCase=0;  # Boolean, set to 1 if Proteomic Ruler is used (when quantifType is known) 
my $encodedDegStrg=($action eq 'export')? encode_utf8('°') : '°';
#if ($action eq 'delete') {&deleteQuantification($analysisID,$selQuantifID); exit;} # OBSOLETE only for label-based 2ndary quantifications
if ($action eq 'ajaxCompQuantifProt') {&ajaxComputeQuantifiedProteins; exit;}
elsif ($action eq 'ajaxCheckStatus') {&ajaxCheckQuantifStatus; exit;}
elsif ($action eq 'ajaxDispCorrel') {&ajaxDisplayCorrelationMatrix; exit;}
elsif ($action eq 'ajaxDispNormPep') {&ajaxDisplayNormPeptides; exit;}
elsif ($action eq 'ajaxDispMissValues') {&ajaxDisplayMissingValues; exit;}
elsif ($action eq 'ajaxProtStat') {&ajaxProteinStatistics; exit;}
elsif ($action eq 'ajaxDispBiasCoeffList') {&ajaxDisplayBiasCoefficientsList; exit;}
elsif ($action eq 'ajaxDispBiasCoeffPlot') {&ajaxDisplayBiasCoefficientsPlot; exit;}
elsif ($action eq 'ajaxPepRatio') {&ajaxPeptideRatios; exit;}
elsif ($action eq 'ajaxPepMaxQuant') {&ajaxPeptideMaxQuant; exit;}
elsif ($action eq 'ajaxPepSwath') {&ajaxPeptideSwath; exit;}
elsif ($action eq 'ajaxPepAbund') {&ajaxPeptideAbundance; exit;}
elsif ($action eq 'ajaxPepProtRuler') {&ajaxPeptideProtRuler; exit;}
elsif ($action eq 'ajaxListProt') {&ajaxListSelectedProteins; exit;}
elsif ($action eq 'ajaxListDecGraph') {&ajaxListDecorateGraph; exit;}
elsif ($action eq 'ajaxGoTerms') {&ajaxGetGoTerms; exit;}
elsif ($action eq 'ajaxTermProt') {&ajaxGetGoTermProteins; exit;}
elsif ($action eq 'ajaxSaveProt') {&ajaxManageSaveProteins; exit;}
elsif ($action eq 'ajaxRestrictList') {&ajaxRestrictProteinList; exit;}
elsif ($action eq 'ajaxManageTemplate') {&ajaxManageTemplate; exit;}
elsif ($action eq 'ajaxDispCustProt') {&ajaxDisplayCustomProteins; exit;}
elsif ($action eq 'ajaxConvIdent') {&ajaxSearchConvertIdentifier; exit;}
elsif ($action eq 'ajaxDispProtQty') {&ajaxDisplayProteinQuantities; exit;}

my $dispResults=param('displayRes') || 0; # display form was submitted

####>Connect to the database
my $dbh=&promsConfig::dbConnect;
my $projectID=($selQuantifID)? &promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION') : &promsMod::getProjectID($dbh,$analysisID,'ANALYSIS');
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
$projectAccess=${$userInfo[2]}{$projectID};

###>Edition submission (internal-analysis only)<###
if ($action eq 'edit') {
	my $qName=$dbh->quote(param('name'));
	$dbh->do("UPDATE QUANTIFICATION SET NAME=$qName WHERE ID_QUANTIFICATION=$selQuantifID");
	$dbh->commit;
	$action='display';
}


################
####>Main 1<####
################
my (%anaQuantifList,%quantificationTypes,@quantificationInfo,%methodList,%quantifModifInfo,%modificationContexts);
my ($quantifSoftware,$algoVersion,$isPepIntensity,$quantifMethDesc,$title,$isModifQuantif,$isMultiModifQuantif,$isFreeResQuantif,$isProteinByPTM); # ,$quantifModID
my ($focus,$quantifAnnot,$quantifStatus,$quantifMethodID,$updateDate,$updateUser);
#my $highlighMatchStrg='';
my $modifName='';
my $focusKeyWord='';
my $visibilityStrg='';
my ($workbook,%itemFormat); # globals for export to Excel

if ($call eq 'ana') {
	$title=($action eq 'export')? 'Protein Quantification' : 'List of Single-Analysis Quantifications';
	my $sthAQ=$dbh->prepare("SELECT Q.ID_QUANTIFICATION,Q.NAME,Q.FOCUS,STATUS,M.ID_QUANTIFICATION_METHOD,M.NAME,M.CODE,M.DES,Q.QUANTIF_ANNOT FROM ANA_QUANTIFICATION A,QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE ID_ANALYSIS=$analysisID AND A.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_DESIGN IS NULL AND Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD");
	$sthAQ->execute;
	while (my ($quantifID,$name,$focus,$status,$methodID,$methodName,$methodCode,$methodDesc,$qAnnot)=$sthAQ->fetchrow_array) {
		$name='No name' unless $name;
		@{$anaQuantifList{$quantifID}}=($name,$focus,$status,$methodID,$qAnnot);
		@{$methodList{$methodID}}=($methodName,$methodCode,$methodDesc);
		my $typeName=($focus eq 'peptide' && $methodCode=~/XIC|SIN|SILAC|ITRAQ|TMT|TDA|DIA/)? 'Peptide quantification'
			: ($focus eq 'peptide' && $methodCode eq 'DIA')? 'SWATH-MS/DIA'
			: ($focus eq 'peptide')? 'Peptide ratio'
			: ($focus eq 'protein' && $methodCode=~/PROT_RATIO/)? 'Protein ratio' # should never be TNPQ (call=ana)
			: ($focus eq 'protein')? 'Protein quantification' # EMPAI
			: 'Unknown';
		push @{$quantificationTypes{$typeName}},$quantifID;
		#push @{$quantificationTypes{'Peptide quantification'}},$quantifID if ($designID && $methodCode=~/PROT_RATIO|TNPQ/);
	}
	$sthAQ->finish;
	if ($action eq 'select') { #'display'; # $selQuantifID;
		$dbh->disconnect;
	}

	if ($selQuantifID) {
		($quantifType,$quantifMethDesc)=@{$methodList{$anaQuantifList{$selQuantifID}[3]}}[1,2];
		($focus,$quantifAnnot,$quantifStatus,$quantifMethodID,$updateDate,$updateUser)=$dbh->selectrow_array("SELECT FOCUS,QUANTIF_ANNOT,STATUS,ID_QUANTIFICATION_METHOD,UPDATE_DATE,UPDATE_USER FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	}
	$isPepIntensity=0;
}
else { # called from quanti frame
	
	#($quantifType,$designID,$quantifMethDesc,$quantifModID,$modifName,$interName,$synoName,$focus,$quantifAnnot,$quantifStatus,$quantifMethodID,$updateDate,$updateUser)=
	#$dbh->selectrow_array("SELECT QM.CODE,ID_DESIGN,QM.DES,M.ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,FOCUS,QUANTIF_ANNOT,STATUS,QM.ID_QUANTIFICATION_METHOD,UPDATE_DATE,UPDATE_USER FROM QUANTIFICATION Q
	#		LEFT JOIN MODIFICATION M ON Q.ID_MODIFICATION=M.ID_MODIFICATION
	#		INNER JOIN QUANTIFICATION_METHOD QM ON Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD
	#		WHERE ID_QUANTIFICATION=$selQuantifID");
	
	(my $qID,$quantifType,$designID,$quantifMethDesc,my $quantifModID,my $multiModifStrg,$focus,$quantifAnnot,$quantifStatus,$quantifMethodID,$updateDate,$updateUser)=
	$dbh->selectrow_array("SELECT Q.ID_QUANTIFICATION,QM.CODE,ID_DESIGN,QM.DES,Q.ID_MODIFICATION,GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ','),FOCUS,QUANTIF_ANNOT,STATUS,QM.ID_QUANTIFICATION_METHOD,UPDATE_DATE,UPDATE_USER FROM QUANTIFICATION Q
			LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
			INNER JOIN QUANTIFICATION_METHOD QM ON Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD
			WHERE Q.ID_QUANTIFICATION=$selQuantifID
			GROUP BY Q.ID_QUANTIFICATION");
	$protRulerCase=1 if $quantifType eq 'PROT_RULER';
	
	$action='proteinView' if ($action eq 'select' && $designID);
	$title=($designID)? 'Design-based Quantification of ' : 'Single-Analysis Quantification of ';
	
	my $ratioCorrStrg=($quantifAnnot=~/::INTRA_PROT_REF=/)? 'Fold Change-corrected ' : '';
	
	($isModifQuantif,$isMultiModifQuantif)=&promsQuantif::getQuantifModificationInfo($dbh,$quantifModID,$multiModifStrg,\%quantifModifInfo);
	$isFreeResQuantif=0; # default
	if ($isModifQuantif) {
		$isFreeResQuantif=1 if $quantifModifInfo{NAME}{-1};
		$biasCorrLevelNames{'prot'}='site'; $biasCorrLevelNames{'ion_prot'}='peptide & site';
		$modifName=($isMultiModifQuantif)? 'Multi' : $quantifModifInfo{NAME}{ (keys %{$quantifModifInfo{NAME}})[0] };
		$focusKeyWord=($isFreeResQuantif)? '' : '-sites';
		$title.=($action eq 'export')? "$ratioCorrStrg".(join('/',(sort{lc($a) cmp lc($b)} values %{$quantifModifInfo{EXPORT_NAME}})))."$focusKeyWord"
			  : "<FONT color=\"#DD0000\">$ratioCorrStrg".(join('/',(sort{lc($a) cmp lc($b)} values %{$quantifModifInfo{NAME}})))."</FONT>$focusKeyWord";
		#$highlighMatchStrg=",'^###-'"; # for matching modProtIDs with protIDs GO in volcanoPlot (for Lists:depends of the list type)
		$dispNumPep=1 if $updateNumPep;
		$title.=' (with Sequence context)' if $quantifAnnot=~/::KEEP_SEQ_CONTEXT=1/;
	}
	if ($quantifAnnot=~/::(PROTEIN_PTM|SITE_CONTEXT)=([^:]+)/) {
		my ($contextCode,$ptmContexts)=($1,$2);
		my @modifs;
		foreach my $ptmInfo (split(';',$ptmContexts)) {
			my ($modID,@pmtContexts)=split('&',$ptmInfo);
			$modID=~s/#//;
			push @modifs,$modID;
			$modificationContexts{$modID}=\@pmtContexts;
		}
		if ($contextCode eq 'PROTEIN_PTM') {
			$isProteinByPTM=1;
			my ($hasPtm,$isMultiPtm)=&promsQuantif::getQuantifModificationInfo($dbh,undef,join(',',@modifs),\%quantifModifInfo); # already done if modif quantif 
			$title.='Proteins enriched by ';
			$title.=($action eq 'export')? join('/',(sort{lc($a) cmp lc($b)} values %{$quantifModifInfo{EXPORT_NAME}}))
				: "<FONT color=\"#DD0000\">".(join('/',(sort{lc($a) cmp lc($b)} values %{$quantifModifInfo{NAME}})))."</FONT>";
		}
	}
	$title.='Proteins' if (!$isModifQuantif && !$isProteinByPTM);
	if ($action ne 'export' && $projectAccess=~/bioinfo|mass/) {
		$title.='&nbsp;<SPAN id="templateSpan">';
		if ($quantifAnnot=~/::TEMPLATE=1/) {
			$title.='<SPAN class="template">Template&nbsp;&nbsp;<INPUT type="button" class="font11" value="Remove" onclick="manageTemplate(\'remove\')"></SPAN>';
		}
		#elsif ($quantifAnnot=~/SOFTWARE=DIA/ || ($quantifAnnot=~/SOFTWARE=myProMS;(\d+)/ && $1 >= 3)) {
		#	$title.='<INPUT type="button" class="title3" value="Save as template" onclick="manageTemplate(\'add\')">';
		#}
		$title.='</SPAN>';
		if ($quantifStatus > 1) {
			$visibilityStrg='<BR><FONT class="title2" style="color:#DD0000">';
			$visibilityStrg.=($quantifStatus==2)? "Hidden from collaborators" : "For bioinformaticians only";
			$visibilityStrg.='</FONT>';
		}
	}
}

###> Initialize values according to quantif method <###
my (%abundanceMeasures,@aLFQparams);
if ($quantifType=~/^(PROT_ABUNDANCE|MQ|PROT_RULER)$/) {
	my ($measStrg)=($quantifType eq 'PROT_RULER')? $quantifAnnot=~/::OUTPUT_MEAS=([^:]+)/ : $quantifAnnot=~/::ABUND_MEASURES=([^:]+)/;
	if ($measStrg) {
		foreach my $measTypeStrg (split(';',$measStrg)) {
			my ($measType,@measParams)=split(',',$measTypeStrg);
			if ($measType=~/^aLFQ/) {
				@aLFQparams=($measType,@measParams);
				next;
			}
			$dispMeasure=$measType unless $dispMeasure;
			@{$abundanceMeasures{$measType}}=@measParams;
		}
	}
}
unless ($dispMeasure) {
	$dispMeasure=($quantifType eq 'PROT_ABUNDANCE')? 'MEAN_INT' : ($quantifType eq 'MQ')? 'MQ_INT' : ($quantifType eq 'PROT_RULER')? 'COPY_NB' : 'RATIO';
}

$quantifMethDesc=&promsMod::resize($quantifMethDesc,50);
#$ratioParamCode='RATIO' if $quantifType ne 'MQ';
#$quantifModID=0 unless $quantifModID;
$isModifQuantif=0 unless $isModifQuantif;
$isMultiModifQuantif=0 unless $isMultiModifQuantif;

if ($action eq 'export') {
	#################################
	####>Prepare export to Excel<####
	#################################
	my $timeStamp1=strftime("%Y-%m-%d-%H-%M",localtime);
	#my $timeStamp2=strftime("%Y-%m-%d %H:%M",localtime);

	$workbook=Spreadsheet::WriteExcel->new("-");
	$workbook->set_properties(title=>"Protein quantification data",
							  author=>'myProMS server',
							  comments=>'Automatically generated with Perl and Spreadsheet::WriteExcel'
							  );
	$workbook->set_custom_color(40,189,215,255); # light light color  #BDD7FF (new myProMS colors V3.5+)
	#$workbook->set_custom_color(41,128,179,255); # dark blue color  #80B3FF
	%itemFormat=(
			title =>			$workbook->add_format(align=>'center',size=>18,bold=>1,border=>1),
			header =>			$workbook->add_format(align=>'center',size=>10,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			headerR =>			$workbook->add_format(align=>'right',valign=>'top',size=>10,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			headerVCenter =>	$workbook->add_format(align=>'center',valign=>'vcenter',size=>10,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			mergeRowHeader =>	$workbook->add_format(align=>'center',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			mergeColHeader =>	$workbook->add_format(align=>'center',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			mergeColHeaderL =>	$workbook->add_format(align=>'left',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			text =>				$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>0,border=>1),
			textWrap =>			$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>1,border=>1),
			#mergeText =>		$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>0,border=>1),
			mergeColText =>		$workbook->add_format(align=>'left',size=>10,text_wrap=>1,border=>1),
			#numberVis =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,bold=>1),
			#numberHid =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1),
			#numberVisBC =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,bold=>1,color=>'grey'),
			#numberHidBC =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,color=>'grey'),
			#numberVisVP =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,bold=>1,italic=>1),
			#numberHidVP =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,italic=>1),
			number =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1),
			#mergeNumber =>		$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1),
			number1d =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',num_format=>'0.0',border=>1),
#mergeNumber1d =>	$workbook->add_format(size=>10,align=>'center',valign=>'top',num_format=>'0.0',border=>1)
			);

	print header(-type=>"application/vnd.ms-excel",-attachment=>"Protein_quantification_$timeStamp1.xls");

}
else {
	#######################
	####>Starting HTML<####
	#######################
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1); # writes PERL warnings as HTML comments <!-- --> must be activated after header
	print qq
|<HTML>
<HEAD>
<TITLE>Quantification results</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.template {font-size:15px;color:#F5F5F5;background-color:#00A000;border-radius:5px;padding:2px 2px 2px 10px;}
.TD {font-weight:normal;}
.TH {font-weight:bold;}
.LINK {cursor:pointer;}
.highlight{width:250px;}
.trueFilter {font-weight:bold;color:#DD0000;}
.noFilter {font-weight:bold;}
.ghostPeptide {font-weight:normal; font-style:italic;}
.popup {z-index:100;background-color:#FFFFFF;border:solid 3px #999999;padding:5px;position:absolute;display:none;box-shadow:10px 10px 20px #808080;}
</STYLE>
<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/heatMap.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/genericPlot.js"></SCRIPT>
|;
	if ($action eq 'summary') {
		print "<SCRIPT type=\"text/javascript\">\n";
		&promsMod::popupInfo;
	}
	else {
		print qq
|<SCRIPT src="$promsPath{html}/js/local/volcanoPlot2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/peptidePlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
<SCRIPT type="text/javascript">
|;
	&promsMod::popupInfo;
	print qq
|var view='$view';
function displayQuantification(quantifInfoStrg,action) {
	var quantifInfo=quantifInfoStrg.split(':'); // focus:quantifID
	var focus=quantifInfo[0], quantifID=quantifInfo[1];
	if (quantifID==0) {action='select';}
	if (action=='delete' && !confirm('Delete selected quantification?')) {return;} // Obsolete (PP 28/05/20)
	if (focus=='protein') {
		// window.location="$promsPath{cgi}/showProtQuantification.cgi?CALL=$call&id_ana=$analysisID&id_quantif="+quantifID+"&ACT="+action+"&view="+view+"&foldChange=$dispFoldChange&pValue=$dispPvalue&sort=$dispSort"; //&stdDev=\$dispStdDev
		window.location="$promsPath{cgi}/showProtQuantification.cgi?CALL=$call&id_ana=$analysisID&id_quantif="+quantifID+"&ACT="+action+"&view=$view"; //&foldChange=\$dispFoldChange&pValue=\$dispPvalue&sort=$dispSort"; //&stdDev=\$dispStdDev
	}
	else { // peptide
		window.location="$promsPath{cgi}/showPepQuantification.cgi?CALL=$call&id_ana=$analysisID&id_quantif="+quantifID+"&ACT="+action;
	}
}
function sequenceView(id_protein,anaIdStrg) {
	var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+anaIdStrg+"&id_prot="+id_protein+"&msdata="+top.showMSdata;
	top.openProtWindow(winLocation);
}
var spectWin,selectedPep;
function drawSpectrum(pepId,pepInfo) {
	//if (pepId != selectedPep \|\| (spectWin && spectWin.closed)) {
		if (selectedPep && document.getElementById(selectedPep)) { // selectedPep could have been remove by another AJAX call
			document.getElementById(selectedPep).style.color='#000000';
		}
		if (pepId) {document.getElementById(pepId).style.color='#DD0000';} // null if called from SVG;
		selectedPep=pepId;
		var paramString="RID="+pepInfo+"&CALL=pep";
		spectWin=window.open("$promsPath{cgi}/drawSpectrum.cgi?"+paramString,'SpectrumWindow','width=950,height=950,location=no,resizable=yes,scrollbars=yes');
	//}
	spectWin.focus();
}
|;
	}
	print qq
|</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif"> <!-- Do not add onload code here! Use window.onload instead. -->
<CENTER>
<FONT class="title">$title</FONT>$visibilityStrg
<BR>
|;

#print "*** ACTION=$action ***<BR>\n";

	######################################################
	####>Internal quantifications drop-down selection<####
	######################################################
	if ($call eq 'ana') {
		print qq
|<TABLE bgcolor="$darkColor">
<TR><TD class="title2" align=right>Select:</TD><TD><SELECT class="title2" onchange="displayQuantification(this.value,'display')"><OPTION value="protein:0">-= Select =-</OPTION>
|;
		my $prevFocus='';
		foreach my $typeName (sort keys %quantificationTypes) {
			my $nbValidXIC = 0;
			print "<OPTGROUP label=\"$typeName:\">\n";
			foreach my $quantifID (sort{lc($anaQuantifList{$a}[0]) cmp lc($anaQuantifList{$b}[0])} @{$quantificationTypes{$typeName}}) {
				print "<OPTION value=\"$anaQuantifList{$quantifID}[1]:$quantifID\"";
				my $errorFlag='';
				if ($anaQuantifList{$quantifID}[2] == -1 || $anaQuantifList{$quantifID}[2] == 0) {
					print ' disabled'; # just created or running
				}
				elsif ($anaQuantifList{$quantifID}[2] == -2) { # Error
					$errorFlag='***';
				}
				print ' selected' if $quantifID==$selQuantifID;
				my $methodName=$methodList{$anaQuantifList{$quantifID}[3]}[0];
				my $qAnnot=$anaQuantifList{$quantifID}[4];
				my $softwareExt='';
				if ($methodList{$anaQuantifList{$quantifID}[3]}[1] eq 'XIC'){ # Need to distinguish between MassChroQ, MaxQuant and Proteome-Discoverer
					$softwareExt=($qAnnot =~ /SOFTWARE=MCQ/) ? '-MCQ' : ($qAnnot =~ /SOFTWARE=MQ/) ? '-MQ' : ($qAnnot =~ /SOFTWARE=PD/) ? '-PD' : '';
				}
				print ">$errorFlag$anaQuantifList{$quantifID}[0] [$methodName]$softwareExt$errorFlag</OPTION>\n";
				$nbValidXIC++ if ($methodName ne 'TDA' && $methodName ne 'DIA');
			}
			#if ($typeName eq 'Peptide quantification') { #&& $nbValidXIC > 1) {
			#	print "<OPTION value=\"peptide:XIC\">Compare multiple XIC extractions</OPTION>";
			#}
			print "</OPTGROUP>\n";
		}
		#my $notDeletable=1;
		#if ($selQuantifID) {
		#	($notDeletable)=$dbh->selectrow_array("SELECT COUNT(*) FROM PARENT_QUANTIFICATION WHERE ID_PARENT_QUANTIFICATION=$selQuantifID");
		#}
		#my $disabDeleteStrg=($notDeletable || $designID)? ' disabled' : '';

		print qq
|</SELECT></TD>
</TR></TABLE><BR>
|;
		if ($action eq 'select') {
			print "</CENTER>\n</BODY>\n</HTML>\n";
			exit;
		}

		####>Hidden edition form<####
#		print qq
#|<FORM name="editForm" method="POST">
#<INPUT type="hidden" name="ACT" value="edit">
#<INPUT type="hidden" name="CALL" value="ana">
#<INPUT type="hidden" name="view" value="$view">
#<INPUT type="hidden" name="id_ana" value="$analysisID">
#<INPUT type="hidden" name="id_quantif" value="$selQuantifID">
#<INPUT type="hidden" name="foldChgType" value="$foldChangeType">
#<INPUT type="hidden" name="foldChange" value="$dispFoldChange">
#<!--INPUT type="hidden" name="stdDev" value="\$dispStdDev"-->
#<INPUT type="hidden" name="coefVar" value="$dispCV">
#<INPUT type="hidden" name="pValue" value="$dispPvalue">
#<INPUT type="hidden" name="pepType" value="$dispPepType">
#<INPUT type="hidden" name="numPep" value="$dispNumPep">
#<INPUT type="hidden" name="numPepScope" value="$numPepScope">
#<INPUT type="hidden" name="numReplic" value="$dispNumReplic">
#<INPUT type="hidden" name="sort" value="$dispSort">
#<DIV id="editDIV" style="display:none">
#<TABLE bgcolor=$darkColor border=0>
#<TR><TH align=right>&nbsp;Name :</TH><TD bgcolor=$lightColor><INPUT type="text" name="name" value="$anaQuantifList{$selQuantifID}[0]" size="50"></TD></TR>
#<TR><TH colspan=2><INPUT type="submit" value=" Save ">&nbsp;&nbsp;&nbsp;<INPUT type="button" value=" Cancel " onclick="document.getElementById('editDIV').style.display='none'"/></TR>
#</TABLE>
#</DIV>
#</FORM>
#|;
	}
	#print "<DIV id=\"waitDiv\"><BR><BR><BR><BR><BR><FONT class=\"title3\">Fetching data. Please wait...</FONT><BR><IMG src=\"$promsPath{images}/scrollbarGreen.gif\"><BR><BR><BR></DIV>\n";
}

###################################################
####>Fetching data for selected Quantification<####
###################################################
#($quantifType)=$dbh->selectrow_array("SELECT M.CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$selQuantifID") unless $quantifType; # internal quanti (SILAC, iTRAQ)
#if ($quantifType !~ /(PROT_RATIO|TNPQ|XIC|SIN|EMPAI)/) {
#my ($focus,$quantifAnnot,$quantifStatus,$quantifMethodID,$updateDate,$updateUser)=$dbh->selectrow_array("SELECT FOCUS,QUANTIF_ANNOT,STATUS,ID_QUANTIFICATION_METHOD,UPDATE_DATE,UPDATE_USER FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
if ($action=~/summary|export/) {
	if ($updateUser) {
		my $sthUser=$dbh->prepare("SELECT USER_NAME FROM USER_LIST WHERE ID_USER=?");
		$sthUser->execute($updateUser);
		my ($userName)=$sthUser->fetchrow_array;
		$updateUser=$userName || $updateUser;
	}
	else {$updateUser='?';}
}
if ($quantifType eq 'XIC' && $call eq 'quanti') { # To print the XIC in a good way
	$focus='protein';
	$action='xicView';
}
my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
my ($labelType)=($labelStrg)? ($labelStrg=~/LABEL=(.+)/) : ('FREE');
$labelType='FREE' unless $labelType;
$labelType=uc($labelType);

###>Number of proteins and sites<###
my $protQuantifiedStrg='';
if ($quantifStatus >= 1 && $quantifType !~ /EMPAI|SIN/) { # |MQ
	my $paramCode=($quantifType=~/^(MQ|PROT_RULER|PROT_ABUNDANCE)$/)? $dispMeasure : $ratioParamCode; # 'NUM_PEP_TOTAL';
	$protQuantifiedStrg="<SPAN id=\"protQuantifiedSPAN\"><INPUT type=\"button\" value=\"...\" onclick=\"ajaxComputeQuantifiedProteins('$paramCode')\"/> proteins quantified</SPAN>";
}

##################################
####>Protein restriction list<####
##################################
my (%restrictList,$isoformList);
if ($restrictListID) {
	$isoformList=&promsQuantif::fetchCustomList($dbh,$restrictListID,\%restrictList); # auto-detection site/prot if "no site" option is not specified
	#my $useSites;
	#if ($isModifQuantif) {
	#	($useSites)=$dbh->selectrow_array("SELECT 1 FROM CATEGORY WHERE ID_CATEGORY=$restrictListID AND LIST_TYPE='SITE'");
	#}
	#if ($useSites) {
	#	&promsQuantif::fetchSiteList($dbh,$restrictListID,\%restrictList);
	#}
	#else { # Normal proteins
	#	my $sthRP=$dbh->prepare("SELECT DISTINCT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$restrictListID");
	#	$sthRP->execute;
	#	while (my($protID)=$sthRP->fetchrow_array) {
	#		$restrictList{$protID}=1;
	#	}
	#	$sthRP->finish;
	#}
}

my $sthQP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,NAME,CODE FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=$quantifMethodID");
$sthQP->execute;
while (my ($paramID,$paramName,$paramCode)=$sthQP->fetchrow_array) {
	@{$quantifParamInfo{$paramCode}}=($paramID,$paramName);
	$quantifParamID2code{$paramID}=$paramCode;
}
$sthQP->finish;

########################
####>DESIGN QUANTIF<####
########################
my $numStates=1; # default
if ($designID) {

	####>View TNPQ/RATIO (Protein data)<####
	my ($usedAnaID)=($analysisID)? ($analysisID) : $dbh->selectrow_array("SELECT ANALYSIS.ID_ANALYSIS FROM ANALYSIS,ANA_QUANTIFICATION WHERE ANA_QUANTIFICATION.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ID_QUANTIFICATION IN (SELECT ID_PARENT_QUANTIFICATION FROM PARENT_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID) ORDER BY NAME ASC LIMIT 0,1");# Select one of the ANALYSIS related to that TNPQ quantification so as to make sequenceView method work...
#my ($quantifAnnot,$quantifMethodID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_QUANTIFICATION_METHOD FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	#my ($parentQuantifType)=$dbh->selectrow_array("SELECT CODE FROM QUANTIFICATION,PARENT_QUANTIFICATION,QUANTIFICATION_METHOD WHERE QUANTIFICATION.ID_QUANTIFICATION=PARENT_QUANTIFICATION.ID_PARENT_QUANTIFICATION AND QUANTIFICATION.ID_QUANTIFICATION_METHOD=QUANTIFICATION_METHOD.ID_QUANTIFICATION_METHOD AND PARENT_QUANTIFICATION.ID_QUANTIFICATION=$selQuantifID LIMIT 1");
#my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
	my %labelingInfo;
	#($labelingInfo{'LABELTYPE'})=($labelStrg)? $labelStrg=~/LABEL=(.+)/ : ('FREE');
	foreach my $infoStrg (@labelInfo) {
		my ($setting,$valueStrg)=split('=',$infoStrg);
		if ($setting eq 'RATIOS' && $designID) { # $parentQuantifType && $parentQuantifType eq 'XIC'
			$valueStrg=~s/\#//g; # for design experiments, RATIOS are displayed with '#' for condition IDs
		}
		if ($setting eq 'PEPTIDES') {$valueStrg=~s/\#//g;} # remove ID tag of selected modifications if any
		@{$labelingInfo{$setting}}=split(';',$valueStrg);
	}
	@{$labelingInfo{'IS_MODIF'}}=($isModifQuantif); # add modif status to be passed to &listProteinRatios

	##>Non-quantified proteins (myProMS v3)
	if ($quantifStatus >= 1 && $labelingInfo{'LOST_PROTEINS'}) {
		my $proteoStrg=($isModifQuantif)? "$modifName$focusKeyWord" : 'protein';
		my $numLostLFQ=''; # strg or number!
		if ($abundanceMeasures{'MY_LFQ'}) {
			$numLostLFQ=($labelingInfo{'LOST_PROTEINS_LFQ'})? $labelingInfo{'LOST_PROTEINS_LFQ'}[0] : '?'; # myProMS v3.5
		}
		if ($labelingInfo{'LOST_PROTEINS'}[0]) {
			$protQuantifiedStrg.=" / <B>$labelingInfo{LOST_PROTEINS}[0]</B> $proteoStrg";
			$protQuantifiedStrg.=($labelingInfo{'LOST_PROTEINS'}[0] == 1)? ' was ' : ($focusKeyWord)? 's were ' : ' were ';
			$protQuantifiedStrg.=' not quantified';
			$protQuantifiedStrg.=" and additional $numLostLFQ for LFQ)" if $numLostLFQ;
		}
		else { # 0
			$proteoStrg.='s' unless $proteoStrg=~/s$/;
			$protQuantifiedStrg.=" / <B>All</B> $proteoStrg were quantified";
			$protQuantifiedStrg.=" except $numLostLFQ excluded for LFQ" if $numLostLFQ;
		}
	}

	$ratioType=($labelingInfo{'RATIO_TYPE'})? $labelingInfo{'RATIO_TYPE'}[0] : (!$labelingInfo{'RATIOS'})? 'None' : 'Ratio';
	$numPepCode='PEPTIDES' if ($quantifType eq 'MQ' && $numPepCode=~/_PEP_USED/); # (NUM/DIST)_PEP_USED not available in no-ratio MaxQuant
	$numStates=($labelingInfo{'PARENT_Q'})? (scalar @{$labelingInfo{'PARENT_Q'}}) : (scalar @{$labelingInfo{'STATES'}});  # Take states info from parent quantif for Proteomic Ruler (no states in that case), otherwise it is directly available
	$algoVersion=($ratioType eq 'Ratio')? 1 : (!$labelingInfo{'SOFTWARE'})? 2 : ($labelingInfo{'SOFTWARE'}[0] eq 'myProMS')? $labelingInfo{'SOFTWARE'}[1] : 0; # myProMS version ONLY!
	$isPepIntensity=($quantifType eq 'PROT_ABUNDANCE' || ($labelingInfo{'MEAN_STATE'} && $labelingInfo{'MEAN_STATE'}[0]==1) || ($algoVersion==2 && $labelType eq 'FREE') || ($algoVersion>=3 && $labelingInfo{'ALGO_TYPE'}[0] eq 'PEP_INTENSITY'))? 1 : 0; # 0 for v1 Ratio, v2 Label & v3 PEP_RATIO
	#$quantifSoftware=($algoVersion > 0)? "myProMS v$algoVersion" : ($labelingInfo{'SOFTWARE'}[0] eq 'MQ')? 'MaxQuant' : $labelingInfo{'SOFTWARE'}[0];
	$quantifSoftware=($algoVersion > 0)? 'myProMS' : ($labelingInfo{'SOFTWARE'}[0] eq 'MQ')? 'MaxQuant' : $labelingInfo{'SOFTWARE'}[0];
	if ($algoVersion && $ratioType=~/Ratio/) { # myProMS version ONLY!
		$quantifMethDesc=($isPepIntensity)?	'Protein ratio based on peptide intensity' : 'Protein ratio based on peptide ratio'; # overwrites $quantifMethDesc
	}
	%ratioCode=&generateRatioCodes($labelType,$ratioType,$numStates,$isPepIntensity,$algoVersion) if $ratioType ne 'None';

	###>Gene ontology annotation
	my %goAnalyses;
	if ($view eq 'graph' && $dispResults) {
		my $sthGO=$dbh->prepare("SELECT ID_GOANALYSIS,NAME,ASPECT FROM GO_ANALYSIS WHERE ID_EXPERIMENT=(SELECT D.ID_EXPERIMENT FROM DESIGN D,QUANTIFICATION Q WHERE D.ID_DESIGN=Q.ID_DESIGN AND Q.ID_QUANTIFICATION=$selQuantifID) AND GOA_TYPE=\"graph\"");
		$sthGO->execute;
		while (my ($goID,$name,$aspectStrg)=$sthGO->fetchrow_array) {
			@{$goAnalyses{$goID}}=($name,$aspectStrg);
		}
	}

	###>Displaying summary
	#my $numRatios=scalar @{$labelingInfo{'RATIOS'}};
	my %stateInfo=&printQuantificationSummary($dbh,$usedAnaID,\%labelingInfo,\%goAnalyses);


	#### ---> From here action cannot be 'summary' --->
	$dbh->disconnect;
	print "<DIV id=\"waitDiv\"><BR><BR><BR><BR><BR><FONT class=\"title3\">Fetching data. Please wait...</FONT><BR><IMG src=\"$promsPath{images}/scrollbarGreen.gif\"><BR><BR><BR></DIV>\n" if $action ne 'export';


	###############################################
	#####>Fetching protein quantification data<####
	###############################################
	my (%quantifValues,%proteinInfo,%restrictProteins,%dispModifSites);
	if (scalar keys %restrictList) { # && $listAction eq 'restrict'
		foreach my $modProtID (keys %restrictList) {
			my ($protID0)=($modProtID=~/^(\d+)/);
			$restrictProteins{$protID0}=1;
			@{$proteinInfo{$protID0}}=() if $listAction eq 'restrict'; # in case prot in restrict are not quantified (=> total num proteins in list is displayed)
		}
	}

	my $refRestrictInfo=($restrictListID)? {ACTION=>$listAction,SITE=>$isoformList,PROT=>\%restrictProteins,LIST=>\%restrictList} : undef;
	my ($pvalueCode,$minPvalue)=&fetchDesignQuantificationData($projectID,$view,\%dispRatios,\%dispStates,\%ratioCode,\%labelingInfo,\%stateInfo,\%quantifParamID2code,\%quantifValues,\%proteinInfo,\%dispModifSites,\%quantifModifInfo,$refRestrictInfo);
		

	####################
	####>Export XLS<####
	####################
	if ($action eq 'export') {
		&exportProteinList($labelType,$ratioType,\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,\%dispModifSites,$pvalueCode);
		$workbook->close();
		exit;
	}
	#print qq |<DIV id="waitDiv"><BR><BR><BR><BR><BR><FONT class="title3">Fetching data. Please wait...</FONT><BR><IMG src="$promsPath{images}/scrollbarGreen.gif"></DIV><SCRIPT type="text/javascript">document.getElementById('waitDiv').style.display='none'</SCRIPT>|;
	print qq |<SCRIPT type="text/javascript">document.getElementById('waitDiv').style.display='none'</SCRIPT>|;
	#my $numRatios=scalar @{$labelingInfo{RATIOS}};

	########################
	####>Graphical view<####
	########################
	if ($view eq 'graph') {
		if ($protRulerCase) {
			&displayProtRulerPlot(\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo);
		}
		elsif ($ratioType eq 'None') {
			#&displayGenericPlot(\%quantifValues,\%proteinInfo,'MQ_INT','Intensity',$usedAnaID);
			&displayIntensityPlot(\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,\%dispModifSites);
		}
		elsif ($quantifSoftware eq 'MaxQuant') {# Ratios but no p-values!!!
			&displayRatioPlot(\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo);
		}
		else { # Ratios with p-values 
			&displayVolcanoPlot(\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,\%dispModifSites,$pvalueCode,$minPvalue);
		}
	}

	###################
	####>List view<####
	###################
	elsif ($view eq 'list') { # View type is list
		print qq
|<SCRIPT type="text/javascript">
window.onload=function() {
ajaxUpdateRestrict();
}
</SCRIPT>
|;
		if ($quantifSoftware eq 'MaxQuant' || $quantifType eq 'PROT_ABUNDANCE') {
			&listProteinIntensities($labelType,\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,\%dispModifSites);
		}
		elsif ($quantifType eq 'PROT_RULER') {
			&listProteomicRuler($labelType,\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo);
		}
		else {
			&listProteinRatios($labelType,$ratioType,\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,\%dispModifSites,$pvalueCode); #$labelingInfo{FDR_ALPHA}[0]/100,$dispStdDev
		}
	}
	&endHTML;
} # end of design quantif

######################
####>EMPAI or SIN<####
######################
elsif ($quantifType =~ /EMPAI|SIN/) {

	print "<DIV id=\"waitDiv\"><BR><BR><BR><BR><BR><FONT class=\"title3\">Fetching data. Please wait...</FONT><BR><IMG src=\"$promsPath{images}/scrollbarGreen.gif\"><BR><BR><BR></DIV>\n" if $action ne 'export';

	my @qpcodes=($quantifType eq 'EMPAI')? ('EMPAI','EMPAI_MOL','EMPAI_MR') : ('SIN_SIN');
	my ($usedAnaID)=($analysisID)? ($analysisID) : $dbh->selectrow_array("SELECT ID_ANALYSIS FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID"); # internal quantif
	##my $sthProtQ=$dbh->prepare("SELECT ID_PROTEIN,QUANTIF_VALUE,ID_QUANTIF_PARAMETER FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID AND QUANTIF_VALUE IS NOT NULL");
	my $dbhLite=&promsQuantif::dbConnectProteinQuantification($selQuantifID,$projectID);
	my $sthProtQ=$dbhLite->prepare("SELECT ID_PROTEIN,QUANTIF_VALUE,ID_QUANTIF_PARAMETER FROM PROTEIN_QUANTIFICATION WHERE QUANTIF_VALUE IS NOT NULL");		
	my ($GNidentID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
	my $sthProtI=$dbh->prepare("SELECT P.ID_PROTEIN,ALIAS,NUM_PEP,MW,PROT_DES,ORGANISM,GROUP_CONCAT(DISTINCT MI.VALUE ORDER BY IDENT_RANK SEPARATOR ',')
									FROM PROTEIN P
									INNER JOIN ANALYSIS_PROTEIN AP ON AP.ID_PROTEIN=P.ID_PROTEIN 
									LEFT JOIN MASTERPROT_IDENTIFIER MI ON P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND MI.ID_IDENTIFIER=$GNidentID
									WHERE VISIBILITY>=1 AND ID_ANALYSIS=$usedAnaID
									GROUP BY P.ID_PROTEIN");
	#my $sthGN=$dbh->prepare("SELECT MI.VALUE FROM PROTEIN P,MASTERPROT_IDENTIFIER MI WHERE P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND MI.ID_IDENTIFIER=$GNidentID AND P.ID_PROTEIN=? ORDER BY MI.IDENT_RANK");

	$sthProtI->execute;
	my (%proteinInfo,%quantifValues); #,@protOrder
	while (my ($protID,@info)=$sthProtI->fetchrow_array) {
		#next if ($restrictListID && !$restrictList{$protID});
		if ($restrictListID) {
			if ($listAction eq 'restrict') {next if !$restrictList{$protID};}
			else {next if $restrictList{$protID};} # exclude
		}
		$info[2]=($info[2])? 1*(sprintf "%.1f",$info[2]/1000) : 0; # MW
		my @geneList=split(',',$info[5]) if $info[5]; # genes: convert strg to ref of array
		$info[5]=\@geneList;
		@{$proteinInfo{$protID}}=@info;
		# my @geneList;
		# $sthGN->execute($protID);
		# while (my ($gene)=$sthGN->fetchrow_array) {push @geneList,$gene;}
		# push @{$proteinInfo{$protID}},\@geneList;
		foreach my $paramC (@qpcodes) { $quantifValues{$protID}{$paramC}=0; } # for later sort
	}
	$sthProtI->finish;
	#$sthGN->finish;
	$sthProtQ->execute;
	my %quantifProt;
	while (my ($protID,$qValue,$paramID)=$sthProtQ->fetchrow_array) {
		next unless $proteinInfo{$protID}; # skip hidden & restricted protein even if they were quantified!
		#push @protOrder,$protID;
		$quantifValues{$protID}{$quantifParamID2code{$paramID}}=$qValue;
		$quantifProt{$protID}=1;
	}
	$sthProtQ->finish;
	my $numQuantProteins=scalar keys %quantifProt;
	#&displayIntensityPlot(\%quantifValues,\%proteinInfo,'MQ_INT','Intensity',$usedAnaID) if ($quantifType eq 'MQ' && $action ne 'export');
	my $locationStrg='';
	if ($action eq 'export') {
		my @quantifInfo=&promsMod::getItemInfo($dbh,'QUANTIFICATION',$selQuantifID);
		foreach my $it (@quantifInfo) {
			$locationStrg.=' > ' if $locationStrg;
			$locationStrg.=$it->{'NAME'};
		}
	}

	$dbh->disconnect;
	$dbhLite->disconnect;

	my $numTotProteins=scalar keys %proteinInfo;
	my $proteinTitle="$numQuantProteins/$numTotProteins Proteins";
	my $colName=($quantifType eq 'EMPAI')? 'emPAI' : ($action eq 'export')? 'SIN' : 'SI<SUB>N</SUB>'; #:($quantifType eq 'MQ')? 'MaxQuant'  <BR>extracted  <BR>computed
	my $nbCol=($quantifType eq 'EMPAI')? 8 : 6;
	#my ($lightColor,$darkColor)=&promsConfig::getRowColors;

	if ($action eq 'export') {
		####<Start printing>####
		my $worksheet=$workbook->add_worksheet('Results');

		##>Table header<##
		my $xlsRow=0;
		my $xlsCol=0;
		$worksheet->set_column(0,0,30); # identifier col length
		$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow,$xlsCol+$nbCol,'myProMS '.$colName.' quantification : '.$locationStrg,$itemFormat{'mergeColHeaderL'});
		$worksheet->write_string(++$xlsRow,$xlsCol,$proteinTitle,$itemFormat{'header'});
		$worksheet->write_string($xlsRow,++$xlsCol,'Gene & Synonyms',$itemFormat{'header'});
		if ($quantifType eq 'EMPAI'){
			$worksheet->write_string($xlsRow,++$xlsCol,'emPAI',$itemFormat{'header'});
			$worksheet->write_string($xlsRow,++$xlsCol,'emPAI (Mol %)',$itemFormat{'header'});
			$worksheet->write_string($xlsRow,++$xlsCol,'emPAI (Mr %)',$itemFormat{'header'});
		}
		#elsif ($quantifType eq 'MQ'){
		#	$worksheet->write_string($xlsRow,++$xlsCol,'LFQ',$itemFormat{'header'});
		#	$worksheet->write_string($xlsRow,++$xlsCol,'Intensity',$itemFormat{'header'});
		#	$worksheet->write_string($xlsRow,++$xlsCol,'iBAQ',$itemFormat{'header'});
		#	$worksheet->write_string($xlsRow,++$xlsCol,'Spectral Count',$itemFormat{'header'});
		#}
		else{
			$worksheet->write_string($xlsRow,++$xlsCol,$colName,$itemFormat{'header'});
		}
		$worksheet->write_string($xlsRow,++$xlsCol,'Peptides',$itemFormat{'header'});
		$worksheet->write_string($xlsRow,++$xlsCol,'MW (kDa)',$itemFormat{'header'});
		$worksheet->set_column(++$xlsCol,$xlsCol,80); # col length
		$worksheet->write_string($xlsRow,$xlsCol,'Description',$itemFormat{'header'});
		$worksheet->set_column(++$xlsCol,$xlsCol,30); # col length
		$worksheet->write_string($xlsRow,$xlsCol,'Species',$itemFormat{'header'});

		##>Protein list
		foreach my $protID (sort{$quantifValues{$b}{$qpcodes[0]}<=>$quantifValues{$a}{$qpcodes[0]}} keys %proteinInfo) { # @protOrder
			$xlsCol=0;
			#>Identifier
			$worksheet->write_string(++$xlsRow,$xlsCol,$proteinInfo{$protID}[0],$itemFormat{'text'});
			#>Gene
			$worksheet->write_string($xlsRow,++$xlsCol,join(',',@{$proteinInfo{$protID}[5]}),$itemFormat{'text'});
			#>Quantif value
			foreach my $paramC (@qpcodes) {
				if ($quantifValues{$protID}{$paramC}) {$worksheet->write_number($xlsRow,++$xlsCol,$quantifValues{$protID}{$paramC},$itemFormat{'number'});}
				else {$worksheet->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
			}
			#>Num peptides
			$worksheet->write_number($xlsRow,++$xlsCol,$proteinInfo{$protID}[1],$itemFormat{'number'});
			#>MW
			$worksheet->write_number($xlsRow,++$xlsCol,$proteinInfo{$protID}[2],$itemFormat{'number1d'});
			#>Des
			$worksheet->write_string($xlsRow,++$xlsCol,$proteinInfo{$protID}[3],$itemFormat{'textWrap'});
			#>Species
			$worksheet->write_string($xlsRow,++$xlsCol,$proteinInfo{$protID}[4],$itemFormat{'text'});
		}
		$workbook->close();
		exit;
	}
	else {
		print qq
|<SCRIPT type="text/javascript">document.getElementById('waitDiv').style.display='none'</SCRIPT>
<TABLE border=0 cellspacing=0 cellpadding=2>
<TR><TD colspan=5><INPUT type="button" value="Export data" onclick="window.location='./showProtQuantification.cgi?CALL=$call&id_ana=$analysisID&id_quantif=$selQuantifID&ACT=export'"/></TD></TR>
<TR bgcolor="$darkColor">
<TH class="rbBorder" align=left nowrap>&nbsp;$proteinTitle&nbsp;</TH><TH class="rbBorder">&nbsp;Gene&nbsp;</TH>
|;
		if ($quantifType eq 'EMPAI') {
			print qq
|<TH class="rbBorder" nowrap>&nbsp;emPAI&nbsp;</TH>
<TH class="rbBorder" nowrap>&nbsp;emPAI (Mol %)&nbsp;</TH>
<TH class="rbBorder" nowrap>&nbsp;emPAI (Mr %)&nbsp;</TH>
|;
		}
		else {print "<TH class=\"rbBorder\" nowrap>&nbsp;$colName&nbsp;</TH>";}
		print qq
|<TH class="rbBorder" nowrap>&nbsp;Peptides&nbsp;</TH>
<TH class="rbBorder" nowrap>&nbsp;MW<SMALL> kDa</SMALL>&nbsp;</TH>
<TH class="bBorder" width=700>&nbsp;Description - Species&nbsp;</TH>
</TR>
|;
		my $bgColor=$lightColor;
		foreach my $protID (sort{$quantifValues{$b}{$qpcodes[0]}<=>$quantifValues{$a}{$qpcodes[0]}} keys %proteinInfo) { # @protOrder
			print "<TR valign=\"top\" bgcolor=\"$bgColor\">";
			print "<TD class=\"TH\" nowrap=\"nowrap\">&nbsp;<A href=\"javascript:sequenceView($protID,'$usedAnaID')\">$proteinInfo{$protID}[0]</A>&nbsp;</TH>"; # alias
			if (scalar @{$proteinInfo{$protID}[5]} > 1) { # gene
				print "<TH>&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-",join('<BR>&nbsp;&nbsp;-',@{$proteinInfo{$protID}[5]}[1..$#{$proteinInfo{$protID}[5]}]),"</B>')\" onmouseout=\"popout()\">",$proteinInfo{$protID}[5][0],"</A>&nbsp;</TH>";
			}
			elsif ($proteinInfo{$protID}[5][0]) {print "<TH>&nbsp;$proteinInfo{$protID}[5][0]&nbsp;</TH>";}
			else {print '<TH>&nbsp;-&nbsp;</TH>';} # no gene mapped
			foreach my $paramC (@qpcodes) {
				my $value;
				if (defined($quantifValues{$protID}{$paramC})) { # can be 0
					if ($quantifValues{$protID}{$paramC}>=0.1 && $quantifValues{$protID}{$paramC}<1000){
						$value=sprintf '%.2f',$quantifValues{$protID}{$paramC};
					}
					else {
						$value=sprintf '%.2e',$quantifValues{$protID}{$paramC};
					}
				}
				else {$value='-'};
				printf "<TH align=\"right\">&nbsp;$value&nbsp;</TH>";
			}
			print "<TH>&nbsp;$proteinInfo{$protID}[1]&nbsp;</TH><TD class=\"TH\" align=right>$proteinInfo{$protID}[2]&nbsp;</TD><TD>$proteinInfo{$protID}[3] <U><I>$proteinInfo{$protID}[4]</I></U></TD>";
			print "</TR>\n";
			$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
		}
		print "<TR><TD colspan=5><B>End of list.</B></TD></TR>\n</TABLE>\n";
	}
	&endHTML;
} # end of SIN/EMPAI HTML display


##################################
####>INTERNAL LABELED QUANTIF<#### Cannot be MaxQuant
##################################
else { #if ($focus eq 'protein') { #} protein

	if ($call eq 'quanti') { # possible since openProject v1.5.9
		($analysisID,@{$anaQuantifList{$selQuantifID}})=$dbh->selectrow_array("SELECT A.ID_ANALYSIS,Q.NAME,Q.FOCUS,STATUS,M.ID_QUANTIFICATION_METHOD FROM ANA_QUANTIFICATION A,QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE A.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=$selQuantifID");
	}
	my $siteDisplayFormat=($action eq 'export')? 'export' : ($view eq 'list')? 'html' : 'text';
	my %labelingInfo;
	##my $sthQP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,NAME,CODE FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=$anaQuantifList{$selQuantifID}[3]");
	##$sthQP->execute;
	##while (my ($paramID,$paramName,$paramCode)=$sthQP->fetchrow_array) {
	##	@{$quantifParamInfo{$paramCode}}=($paramID,$paramName);
	##}
	##$sthQP->finish;

	########################
	####>PROT_RATIO_PEP<####
	########################
	if ($quantifType eq 'PROT_RATIO_PEP') { # || $quantifType eq 'TNPQ'
		foreach my $infoStrg (@labelInfo) {
			my ($setting,$valueStrg)=split('=',$infoStrg);
			@{$labelingInfo{$setting}}=split(';',$valueStrg);
		}
		$ratioType=($labelingInfo{'RATIO_TYPE'})? $labelingInfo{'RATIO_TYPE'}[0] : 'Ratio';
		$quantifSoftware='myProMS';
		$algoVersion=($ratioType eq 'Ratio')? 1 : 2;
		$numStates=(scalar @{$labelingInfo{'STATES'}}); #-1;
		foreach my $ratioPos (1..$numStates*($numStates-1)/2) {$ratioCode{$ratioPos}=$ratioPos; $ratioCode{-$ratioPos}=-$ratioPos;}

		###>Gene ontology annotation
		my %goAnalyses;
		if ($view eq 'graph' && $dispResults) {
			my @anaParents=&promsMod::getItemInfo($dbh,'ANALYSIS',$analysisID);
			#my $sthGO=$dbh->prepare("SELECT GO.ID_GOANALYSIS,NAME,ASPECT FROM GO_ANALYSIS GO,GOANA_ANALYSIS GA WHERE GO.ID_GOANALYSIS=GA.ID_GOANALYSIS AND ID_ANALYSIS=$analysisID");
			my $sthGO=$dbh->prepare("SELECT ID_GOANALYSIS,NAME,ASPECT FROM GO_ANALYSIS WHERE ID_EXPERIMENT=$anaParents[1]{ID} AND GOA_TYPE=\"graph\"");
			$sthGO->execute;
			while (my ($goID,$name,$aspectStrg)=$sthGO->fetchrow_array) {
				@{$goAnalyses{$goID}}=($name,$aspectStrg);
			}
		}

		###>Displaying summary
		#my $numRatios=scalar @{$labelingInfo{'RATIOS'}};
		my %stateInfo=&printQuantificationSummary($dbh,$analysisID,\%labelingInfo,\%goAnalyses);

		###>Error in Quantification
		if ($action ne 'export') {
			if ($anaQuantifList{$selQuantifID}[2] == -1) {
				$dbh->disconnect;
				print qq
|<SCRIPT type="text/javascript">document.getElementById('waitDiv').style.display='none'</SCRIPT>
<BR><FONT class=\"title2\" color=\"#DD0000\">***An error has occured during quantification***</FONT>
</CENTER>
</BODY>
</HTML>
|;
				exit;
			}
			print "<DIV id=\"waitDiv\"><BR><BR><BR><BR><BR><FONT class=\"title3\">Fetching data. Please wait...</FONT><BR><IMG src=\"$promsPath{images}/scrollbarGreen.gif\"><BR><BR><BR></DIV>\n";
		}

		###>Fetching protein quantification data
		my (%quantifValues,%proteinInfo,%usedIsoforms,%dispModifSites);
		if (scalar keys %restrictList && $listAction eq 'restrict') { # cannot be sites
			foreach my $protID (keys %restrictList) {@{$proteinInfo{$protID}}=();} # in case prot in restrict are not quantified (=>correct num of proteins displayed)
		}
		#my $sthProtQ=$dbh->prepare("SELECT ID_PROTEIN,TARGET_POS,QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID AND ID_QUANTIF_PARAMETER=?");
		##my $sthProtQ=$dbh->prepare("SELECT PQ.ID_PROTEIN,GROUP_CONCAT(RESIDUE,POSITION ORDER BY POSITION SEPARATOR '.'),TARGET_POS,QUANTIF_VALUE
		##								FROM PROTEIN_QUANTIFICATION PQ
		##								LEFT JOIN PROTQUANTIF_MODRES PQMR ON PQ.ID_PROT_QUANTIF=PQMR.ID_PROT_QUANTIF
		##								LEFT JOIN MODIFIED_RESIDUE MR ON PQMR.ID_MODIF_RES=MR.ID_MODIF_RES AND MR.ID_QUANTIFICATION=PQ.ID_QUANTIFICATION
		##								WHERE PQ.ID_QUANTIFICATION=$selQuantifID AND ID_QUANTIF_PARAMETER=? GROUP BY PQ.ID_PROT_QUANTIF");
		my $dbhLite=&promsQuantif::dbConnectProteinQuantification($selQuantifID,$projectID);
		my $sthProtQ=$dbhLite->prepare("SELECT ID_PROTEIN,SITE_CODE,TARGET_POS,QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIF_PARAMETER=?");
		my $pvalueCode=($labelingInfo{'FDR_CONTROL'}[0] eq 'TRUE')? 'PVAL_ADJ' : 'PVAL';
		my @extraParams=($action eq 'export')? ('CONF_LOW','CONF_UP') : ();
		my %formattedModRes;
		my $minPvalue=1;
		foreach my $paramCode ('RATIO','SD_GEO',$pvalueCode,@extraParams) {
			$sthProtQ->execute($quantifParamInfo{$paramCode}[0]);
			while (my ($protID,$modResStrg,$ratioPos,$qValue)=$sthProtQ->fetchrow_array) {
				next unless $dispRatios{$ratioPos};
				#next if ($restrictListID && !$restrictList{$protID});
				if ($restrictListID) {
					if ($listAction eq 'restrict') {next if !$restrictList{$protID};}
					else {next if $restrictList{$protID};} # exclude
				}
				
				@{$proteinInfo{$protID}}=(); # actual number of quantified proteins (only list filtering)
				if ($modResStrg) { # quantif of modification
					#$modResStrg=~s/#//g if $modResStrg=~/^#/; # remove all # if not a multi-modif quantif
					unless ($formattedModRes{$modResStrg}) {
						#@{$formattedModRes{$modResStrg}}=&promsQuantif::formatProteinModificationSites($modResStrg,\%quantifModifInfo,$siteDisplayFormat);
						$formattedModRes{$modResStrg}[0]=&promsQuantif::standardizeSiteCode($modResStrg,$quantifModifInfo{RANK2ID});
						$formattedModRes{$modResStrg}[1]=&promsQuantif::displayModificationSites($formattedModRes{$modResStrg}[0],$quantifModifInfo{DISPLAY},$siteDisplayFormat);
					}
					$protID.='-'.$formattedModRes{$modResStrg}[0];
					$dispModifSites{$protID}=$formattedModRes{$modResStrg}[1];
				}
				#$quantifValues{$protID}{$paramCode}{$ratioPos}=($paramCode eq 'SD_GEO')? $qValue*100 : $qValue;
				$quantifValues{$protID}{$paramCode}{$ratioPos}=$qValue;
				$minPvalue=$qValue if ($paramCode eq $pvalueCode && $qValue > 0 && $qValue < $minPvalue);
				$usedIsoforms{$protID}=1;
			}
		}
		$minPvalue/=10;

		foreach my $paramCode ('DIST_PEP_USED','NUM_PEP_USED','NUM_PEP_TOTAL') { # NO SUPER RATIO EXPECTED FOR INTERNAL QUANTIF!!!
			$sthProtQ->execute($quantifParamInfo{$paramCode}[0]);
			while (my ($protID,$modResStrg,$ratioPos,$qValue)=$sthProtQ->fetchrow_array) {
				next if ($ratioPos && !$dispRatios{$ratioPos}); # no ratioPos for NUM_PEP_TOTAL
				next unless $proteinInfo{$protID}; # NUM_PEP_TOTAL has no ratioPos & could bring additional protID
				#next if ($restrictListID && !$restrictList{$protID});
				if ($modResStrg) { # quantif of modification
					#$modResStrg=~s/#//g if $modResStrg=~/^#/; # remove all # if not a multi-modif quantif
					unless ($formattedModRes{$modResStrg}) {
						@{$formattedModRes{$modResStrg}}=&promsQuantif::formatProteinModificationSites($modResStrg,\%quantifModifInfo,$siteDisplayFormat);
						$formattedModRes{$modResStrg}[0]=&promsQuantif::standardizeSiteCode($modResStrg,$quantifModifInfo{RANK2ID});
						$formattedModRes{$modResStrg}[1]=&promsQuantif::displayModificationSites($formattedModRes{$modResStrg}[0],$quantifModifInfo{DISPLAY},$siteDisplayFormat);
					}
					$protID.='-'.$formattedModRes{$modResStrg}[0];
					next unless $usedIsoforms{$protID};
					$dispModifSites{$protID}=$formattedModRes{$modResStrg}[1];
				}
				if ($ratioPos) { # defined for +/-inf ratios
					$quantifValues{$protID}{$paramCode}{$ratioPos}=$qValue if $dispRatios{$ratioPos};
				}
				else { # undef for quantified proteins
					foreach my $r (keys %dispRatios) { # extend to all displayed ratios
						$quantifValues{$protID}{$paramCode}{$r}=$qValue;
					}
				}
			}
		}
		$sthProtQ->finish;

		###>Fetching protein info
		#my %proteinInfo;
		#my $extraFieldStrg=($view eq 'graph' && $action ne 'export')? ',PROT_LENGTH' : ',MW,PROT_DES,ORGANISM';
		#my $sthProtI=$dbh->prepare("SELECT P.ID_PROTEIN,ALIAS$extraFieldStrg FROM ANALYSIS_PROTEIN A,PROTEIN P WHERE A.ID_PROTEIN=P.ID_PROTEIN AND VISIBILITY>=1 AND ID_ANALYSIS=$analysisID");
		##my $sthProtI=$dbh->prepare("SELECT DISTINCT P.ID_PROTEIN,ALIAS,MW,PROT_DES,ORGANISM,PROT_LENGTH FROM PROTEIN_QUANTIFICATION PQ,PROTEIN P WHERE PQ.ID_PROTEIN=P.ID_PROTEIN AND ID_QUANTIFICATION=$selQuantifID");
		my ($GNidentID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
		my $sthGN=$dbh->prepare("SELECT MI.VALUE FROM PROTEIN P,MASTERPROT_IDENTIFIER MI WHERE P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND MI.ID_IDENTIFIER=$GNidentID AND P.ID_PROTEIN=? ORDER BY MI.IDENT_RANK");
		my @protIdList=keys %proteinInfo;
		while (my @subProtList=splice(@protIdList,0,1000)) { # chuncks of 1000 proteins
			my $protIdStrg=join(',',@subProtList);
			my $sthProtI=$dbh->prepare("SELECT DISTINCT ID_PROTEIN,ALIAS,MW,PROT_DES,ORGANISM,PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN IN ($protIdStrg)");
			$sthProtI->execute;
			while (my ($protID,@info)=$sthProtI->fetchrow_array) {
				#next if ($restrictListID && !$restrictList{$protID});
				if ($restrictListID) {
					if ($listAction eq 'restrict') {next if !$restrictList{$protID};}
					else {next if $restrictList{$protID};} # exclude
				}
				next unless $quantifValues{$protID}; # only quantified proteines
				my @geneList;
				if ($view eq 'list' || $action eq 'export') {
					$info[1]=sprintf "%.1f",$info[1]/1000;
					$sthGN->execute($protID);
					while (my ($gene)=$sthGN->fetchrow_array) {push @geneList,$gene;}
				}
				$info[0]=~s/[,;']/\./g; # Clean identifier: MaxQuant contaminant can be crappy!!!
				@{$proteinInfo{$protID}}=(@info,\@geneList,$analysisID); # to matched multi-ana quantif
				#$proteinInfo{$protID}[1]=sprintf "%.1f",$proteinInfo{$protID}[1]/1000 if ($view eq 'list' || $action eq 'export'); # MW
			}
			$sthProtI->finish;
		}
		$sthGN->finish;

		$dbh->disconnect;
		$dbhLite->disconnect;

		if ($action eq 'export') {
			&exportProteinList($labelType,$ratioType,\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,\%dispModifSites,$pvalueCode);
			$workbook->close();
			exit;
		}

		print qq |<SCRIPT type="text/javascript">document.getElementById('waitDiv').style.display='none'</SCRIPT>|;
		#my $numRatios=scalar @{$labelingInfo{RATIOS}};
		#my $numPepCode=($dispPepType eq 'all')? 'NUM_PEP_USED' : 'DIST_PEP_USED';

		########################
		####>Graphical view<####
		########################
		if ($view eq 'graph') {
			&displayVolcanoPlot(\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,\%dispModifSites,$pvalueCode,$minPvalue);
		}

		###################
		####>List view<####
		###################
		else {
			print qq
|<SCRIPT type="text/javascript">
window.onload=function() {
	ajaxUpdateRestrict();
}
</SCRIPT>
|;
			&listProteinRatios($labelType,$ratioType,\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,\%dispModifSites,$pvalueCode); #$labelingInfo{FDR_ALPHA}[0]/100,$dispStdDev
		}
	}
	&endHTML;
}


#############################################<<<SUBROUTINES>>>###########################################

sub endHTML {
	print qq
|</CENTER>
<DIV id="displayDIV" class="popup"> <!--filter:alpha(opacity=80);opacity:0.8;-->
	<DIV id="infoDIV"></DIV>
</DIV>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">
setPopup();
</SCRIPT>
</BODY>
</HTML>
|;
}

sub fetchDesignQuantificationData { # GLOBAL: $selQuantifID, $ratioType, $dispMeasure, $action, $quantifType, $dispPepType, $numPepScope, $numPepCode, $protRulerCase, $restrictListID
	# returns: pvalueCode,$minPvalue,
	my ($projectID,$finalView,$refDispRatios,$refDispStates,$refRatioCode,$refLabelingInfo,$refStateInfo,$refQuantifParamID2code,$refQuantifValues,$refProteinInfo,$refDispModifSites,$refQuantifModifInfo,$refRestrictInfo)=@_;
	my $siteDisplayFormat=($action eq 'export')? 'export' : ($finalView eq 'list')? 'html' : 'text';
	my ($actionList,$refRestrictProteins,$refRestrictList,$isoformList)=($refRestrictInfo)? ($refRestrictInfo->{ACTION},$refRestrictInfo->{PROT},$refRestrictInfo->{LIST} || {},$refRestrictInfo->{SITE} || 0) : ('nothing',{},{},0);
	my %usedIsoforms;
	my $quantifSoftware=($refLabelingInfo->{'SOFTWARE'})? $refLabelingInfo->{'SOFTWARE'}[0] : 'myProMS'; $quantifSoftware='MaxQuant' if $quantifSoftware eq 'MQ';
	my $pvalueCode='';
	my $minPvalue=1;
	my @params;
	if ($ratioType eq 'None') { # $quantifType eq 'PROT_RULER' or 'PROT_ABUNDANCE'
		@params=($dispMeasure);
		push @params, @{$absCIsCodes{$dispMeasure}} if ($quantifType eq 'PROT_ABUNDANCE' && exists($absCIsCodes{$dispMeasure}));
	}
	elsif ($quantifSoftware eq 'MaxQuant') { # MaxQuant with ratio
		@params=($ratioParamCode,'RATIO_VAR'); # filtering can apply on RATIO_VAR in list view
	}
	else {
		$pvalueCode=($ratioType=~/S\w+Ratio/ || $refLabelingInfo->{'FDR_CONTROL'}[0] eq 'TRUE')? 'PVAL_ADJ' : 'PVAL';
		@params=('RATIO',$pvalueCode);
		push @params,'SD_GEO' if $quantifType eq 'PROT_RATIO_PEP'; # && $quantifSoftware ne 'MaxQuant'); # no SD_GEO for TNPQ!!!
		push @params,('CONF_LOW','CONF_UP') if ($action eq 'export' && $ratioType eq 'Ratio');
		#push @params,'DIST_PEP_USED' if ($ratioType=~/S\w+Ratio/ && ($numPepCode eq 'DIST_PEP_USED' || $action eq 'export')); # for all ratios
	}
	
	my $numStates=($refLabelingInfo->{'PARENT_Q'})? scalar @{$refLabelingInfo->{'PARENT_Q'}} : scalar @{$refLabelingInfo->{'STATES'}};
	my (%replicatesPos,%msmsStatePos);
	if ($numPepScope =~ /^replicate/) {
		&getReplicateTargetPos($numStates,$refDispRatios,$refRatioCode,$refLabelingInfo,$refStateInfo,\%replicatesPos); # stores replicate tgtPos in %dispRatios (dispRatios{rPos}[0,1]) & %replicatesPos
	}
	elsif ($dispPepType =~ /^msms/ && $quantifType eq 'PROT_RATIO_PEP') { # $numPepScope= ratio/msms(One|Both)
		&getMsmsStateTargetPos($refDispRatios,$refRatioCode,$refLabelingInfo,\%msmsStatePos); # stores correponding state tgtPos in %dispRatios (dispRatios{rPos}[0,1]) & %msmsStatePos, 
		&getReplicateTargetPos($numStates,$refDispRatios,$refRatioCode,$refLabelingInfo,$refStateInfo,\%replicatesPos) if $numPepScope ne 'ratio'; # updates %dispRatios with replicates tgtPos (dispRatios{ratioPos}[2,3]) <= Needed to compute #NUM_PEP_USED/state (not in BD!)
	}
	my (@pepParams,		# peptide measures for non-ratio or global to all targetPos
		@pepSupParams);	# peptide measures for ratio
	if ($ratioType=~/S\w+Ratio/) { # No peptide data for super ratios before 18/02/15: use those from linked primary ratios
		@pepSupParams=('NUM_PEP_USED'); # NUM_PEP_USED always needed (point size in graph, displayed in list)
		push @pepSupParams,'DIST_PEP_USED' if ($numPepCode eq 'DIST_PEP_USED' || ($quantifType ne 'MQ' && $action eq 'export')); # for all ratios
		push @pepSupParams,$numPepCode if $numPepCode eq 'NUM_TRUE_USED';
		#my @pepSupParams; # channel-specific params
		# @pepSupParams=($numPepCode) if $numPepCode=~/(NUM|DIST)_(PEP|TRUE)_USED/; # NUM_PEP_USED,DIST_PEP_USED,NUM_TRUE_USED
		# push @pepSupParams,'NUM_PEP_USED' if ($numPepCode ne 'NUM_PEP_USED' && $quantifSoftware ne 'MaxQuant'); # always displayed from now on (PP 15/04/21)
		# if ($finalView eq 'list' || $dispPepType eq 'msmsPc' || $action eq 'export') { # msmsPc requires NUM_PEP_USED for replicates
		# 	push @pepSupParams,'NUM_PEP_USED' if $numPepCode ne 'NUM_PEP_USED'; 
		# 	push @pepSupParams,'DIST_PEP_USED' if ($numPepCode ne 'DIST_PEP_USED' && $quantifSoftware ne 'MaxQuant'); # why needed (PP 13/12/20)?
		# }
		push @pepParams,$numPepCode if $numPepCode !~ /(NUM|DIST)_(PEP|TRUE)_USED/; # MaxQuant PEPTIDES, RAZ_UNI_PEP or UNIQUE_PEP (not targetPos-specific)
		if ($finalView eq 'list' || $action eq 'export') {
			if ($quantifSoftware eq 'MaxQuant') {
				foreach my $param ('PEPTIDES','RAZ_UNI_PEP','UNIQUE_PEP') {
					push @pepParams,$param if $param ne $numPepCode;
				}
			}
			else {@pepParams=('NUM_PEP_TOTAL');}
		}
#print "PEP='@pepParams', SUP='@pepSupParams' SEL=$numPepCode,$quantifParamInfo{$numPepCode}[0]<BR>\n";
	}
	elsif (!$protRulerCase) { # No ratio (No peptide params for Proteomic Ruler)
		if ($quantifSoftware eq 'MaxQuant') {
			if ($action eq 'export' || $finalView eq 'list') {
				@pepParams=('PEPTIDES','RAZ_UNI_PEP','UNIQUE_PEP');
			}
			elsif ($finalView eq 'graph') {
				@pepParams=('PEPTIDES');
				push @pepParams,$numPepCode if $numPepCode ne 'PEPTIDES';
			}
		}
		else {
			if ($action eq 'export') {
				@pepParams=('NUM_PEP_USED','DIST_PEP_USED','NUM_PEP_TOTAL');
				push @pepParams,$numPepCode if $numPepCode !~ /(DIST|NUM)_PEP_USED/;
			}
			else {
				if ($finalView eq 'list') {
					@pepParams=('NUM_PEP_USED','NUM_PEP_TOTAL');
				}
				elsif ($finalView eq 'graph') {
					@pepParams=('NUM_PEP_USED');
				}
				push @pepParams,$numPepCode if $numPepCode ne 'NUM_PEP_USED';
			}
		}
		# if ($finalView eq 'graph') {
		# 	if ($action eq 'export' && $quantifType ne 'MQ') {
		# 		@pepParams=($numPepCode =~ /NUM_PEP_USED/)? ($numPepCode,'DIST_PEP_USED','NUM_PEP_TOTAL') : ($numPepCode =~ /(DIST_PEP|NUM_TRUE)_USED/)? ($numPepCode,'NUM_PEP_USED','NUM_PEP_TOTAL') : ('NUM_PEP_USED','NUM_PEP_TOTAL');
		# 	} else {
		# 		@pepParams=($numPepCode);
		# 	}
		# } else {
		# 	if ($quantifType eq 'MQ') {
		# 		@pepParams=('PEPTIDES','RAZ_UNI_PEP','UNIQUE_PEP');
		# 	}
		# 	else {
		# 		@pepParams=($numPepCode =~ /NUM_PEP_USED/ && $action eq 'export')? ($numPepCode,'DIST_PEP_USED','NUM_PEP_TOTAL') : ($numPepCode =~ /(DIST_PEP|NUM_TRUE)_USED/ || ($action eq 'export'))? ($numPepCode,'NUM_PEP_USED','NUM_PEP_TOTAL') : ('NUM_PEP_USED','NUM_PEP_TOTAL');
		# 	}
		# }
	}
# print "**<BR>\@params=@params,<BR>\@pepParams=@pepParams,<BR>\@pepSupParams=@pepSupParams,<BR>**";
	my $paramIdStrg='';
	foreach my $paramCode (@params,@pepParams,@pepSupParams) {
		$paramIdStrg.=',' if $paramIdStrg;
		$paramIdStrg.=$quantifParamInfo{$paramCode}[0];
		#$quantifParamID2code{$quantifParamInfo{$paramCode}[0]}=$paramCode;
	}

	my %neededTargetPos; # combine all targetPos required & remove duplicates
	if ($ratioType=~/Ratio/) { # $ratioType=~/S\w+Ratio/
		foreach my $rPos (keys %{$refDispRatios}) { # %{$refDispRatios}
			$neededTargetPos{abs($rPos)}=1;
			foreach my $refTgtPos (@{$refDispRatios->{$rPos}}) {
				foreach my $tgtPos (@{$refTgtPos}) {$neededTargetPos{$tgtPos}=1;}
			}
		}
	}
	else {@neededTargetPos{keys %{$refDispStates}}=values %{$refDispStates};}

# print "*PEPSCOPE=$numPepScope PEPCODE=$numPepCode PARAM=$paramIdStrg (P=@params, PEP=@pepParams, SUP=@pepSupParams)<BR>TGT_POS=",join(',',keys %neededTargetPos),"<BR>\n";
	my $targetPosStrg='AND (TARGET_POS IS NULL OR ABS(TARGET_POS) IN ('.join(',',keys %neededTargetPos).'))'; # ABS because some replicate tgtPos can be < 0 to distinguish context!!!

# my $bugProtID=0; #227550; #'236807-0:C877'; #227077; #227223;  #226811; #1/INF        #226662; #INF
	my %allQuantifData;
	my $dbhLite=&promsQuantif::dbConnectProteinQuantification($selQuantifID,$projectID);
	if ($actionList eq 'restrictSQL' && $isoformList) { # Called from AJAX list with sites quantif
		my $sthProtQS=$dbhLite->prepare("SELECT ID_QUANTIF_PARAMETER,ID_PROTEIN,SITE_CODE,IFNULL(TARGET_POS,0),QUANTIF_VALUE
											FROM PROTEIN_QUANTIFICATION
											WHERE ID_PROTEIN=? AND SITE_CODE=? AND ID_QUANTIF_PARAMETER IN ($paramIdStrg) $targetPosStrg");
		foreach my $modProtID (keys %{$refRestrictList}) {
			my ($protID,$modResCode)=($modProtID=~/^(\d+)-(.+)/);
			$sthProtQS->execute($protID,$modResCode);
			while (my ($paramID,@data)=$sthProtQS->fetchrow_array) {
				push @{$allQuantifData{ $refQuantifParamID2code->{$paramID} }},\@data;
# print "*$paramID, @data[2,3]<BR>\n" if $data[0] eq $bugProtID;
			}
		}
		$sthProtQS->finish;
	}
	else { # Non-AJAX call or protein quantif
		my $restrictProtStrg=($actionList eq 'restrictSQL' && !$isoformList)? 'ID_PROTEIN IN ('.join(',',keys %{$refRestrictList}).') AND' : '';
		my $sthProtQ=$dbhLite->prepare("SELECT ID_QUANTIF_PARAMETER,ID_PROTEIN,SITE_CODE,IFNULL(TARGET_POS,0),QUANTIF_VALUE
											FROM PROTEIN_QUANTIFICATION
											WHERE $restrictProtStrg ID_QUANTIF_PARAMETER IN ($paramIdStrg) $targetPosStrg");
										
		$sthProtQ->execute;
		while (my ($paramID,@data)=$sthProtQ->fetchrow_array) {
			push @{$allQuantifData{ $refQuantifParamID2code->{$paramID} }},\@data;
# print "*$bugProtID: $paramID, @data[2,3]<BR>\n" if $data[0] eq $bugProtID;
		}
		$sthProtQ->finish;
	}
	$dbhLite->disconnect;

#print "PARAM='@params', ",scalar @{$allQuantifData{ 'MQ_INT' }},"<BR>\n";

	####<Reshape data (protein/site->param->targetPos=value)>####
	my %formattedModRes;
	foreach my $paramCode (@params) {
		foreach my $refData (@{$allQuantifData{$paramCode}}) {
			my ($protID,$modResStrg,$targetPos,$qValue)=@{$refData};
			#next unless $dispRatios{$targetPos};
			if ($ratioType eq 'None') { # no-ratio MaxQuant or Proteomic Ruler
				if ($protRulerCase && $refLabelingInfo->{'AVG_MODE'}[0] eq 'AVG_ALL') {
					next unless (scalar keys %{$refDispStates});
				} else {
					next unless $refDispStates->{$targetPos};
				}
			}
			else {
				if ($paramCode eq 'RATIO') {
					if ($qValue < $MIN_INF_RATIO_DB) {$qValue=$MIN_INF_RATIO_DB;}
					elsif ($qValue > $MAX_INF_RATIO_DB) {$qValue=$MAX_INF_RATIO_DB;}
				}
				if ($refDispRatios->{-$targetPos}) { # reverse ratio management
					$targetPos=-$targetPos;
					$qValue=1/$qValue if $paramCode=~/^RATIO/;
				}
				elsif (!$refDispRatios->{$targetPos}) {next;}
			}
			my $protID0=$protID;
			if ($modResStrg) { # quantif of modification
				#$modResStrg=~s/#//g if $modResStrg=~/^#/; # remove all # if not a multi-modif quantif
				unless ($formattedModRes{$modResStrg}) {
					#@{$formattedModRes{$modResStrg}}=&promsQuantif::formatProteinModificationSites($modResStrg,\%quantifModifInfo,$siteDisplayFormat);
					$formattedModRes{$modResStrg}[0]=&promsQuantif::standardizeSiteCode($modResStrg,$refQuantifModifInfo->{RANK2ID});
					$formattedModRes{$modResStrg}[1]=&promsQuantif::displayModificationSites($formattedModRes{$modResStrg}[0],$refQuantifModifInfo->{DISPLAY},$siteDisplayFormat);
				}
				$protID.='-'.$formattedModRes{$modResStrg}[0];
				$refDispModifSites->{$protID}=$formattedModRes{$modResStrg}[1];
			}
			if ($refRestrictInfo) {
				if ($actionList eq 'restrict') {next if (($isoformList && !$refRestrictList->{$protID}) || (!$isoformList && !$refRestrictProteins->{$protID0}));}
				else {next if (($isoformList && $refRestrictProteins->{$protID}) || (!$isoformList && $refRestrictProteins->{$protID0}));} # exclude
			}
			@{$refProteinInfo->{$protID0}}=(); # actual number of quantified proteins (only list filtering)
			#$quantifValues{$protID}{$paramCode}{$targetPos}=($paramCode eq 'SD_GEO')? $qValue*100 : $qValue;
			$refQuantifValues->{$protID}{$paramCode}{$targetPos}=$qValue;
			$minPvalue=$qValue if ($paramCode eq $pvalueCode && $qValue > 0 && $qValue < $minPvalue);
			$usedIsoforms{$protID}=1;
		}
	}
	$minPvalue/=10;
	if ($ratioType=~/S\w+Ratio/) { # No peptide data for super ratios before 18/02/15: use those from linked primary ratios
		# TODO: Simplify (remove linkedRatios) to make compatible with only quantif >18/02/15
		my %linkedRatios; # for SuperRatio only
		if ($ratioType eq 'SuperRatio' && $numPepScope !~ /^replicate/ && $dispPepType !~ /^msms/) { # if scope is replicateXXX/msmsXXX: all required extra targetPos are already in %replicatesPos/%msmsStatePos
			my $rPos=0;
			foreach my $ratioData (@{$refLabelingInfo->{'RATIOS'}}) { # Find linked primary ratios
				$rPos++;
				my ($testStateCode,$refStateCode)=split(/\//,$ratioData);
				if ($testStateCode=~/%(\d+)/) { # ratio of ratio
					#next unless $dispRatios{$rPos};
					next if (!$refDispRatios->{$rPos} && !$refDispRatios->{-$rPos});
					my $linkedTestPos=$1;
					my ($linkedRefPos)=($refStateCode=~/%(\d+)/);
					$testStateCode=~s/%\d+//;
					$refStateCode=~s/%\d+//;
					push @{$linkedRatios{$linkedTestPos}},$rPos;
					push @{$linkedRatios{$linkedRefPos}},$rPos;
				}
			}
		}
		foreach my $paramCode (@pepSupParams) {
			foreach my $refData (@{$allQuantifData{$paramCode}}) {
				my ($protID,$modResStrg,$targetPos,$qValue)=@{$refData};
				next unless $refProteinInfo->{$protID}; # $linkedRatios{$targetPos} could bring additional protID NOT in $proteinInfo{$protID}
				my $trueTargetPos=$targetPos; # reverse ratio management
				#next if (!$refDispRatios->{$targetPos} && !$linkedRatios{$targetPos});
				if ($ratioType eq 'None') { # no-ratio MaxQuant
					next unless $refDispStates->{$targetPos};
				}
				else {
					$targetPos=-$targetPos if $refDispRatios->{-$targetPos}; # -N
					if ($numPepScope =~ /^replicate/) {
						# next unless $replicatesPos{abs($trueTargetPos)}; # abs() for +/-Inf context if any
						next if (!$refDispRatios->{$targetPos} && !$replicatesPos{abs($trueTargetPos)}); # allows (NUM|DIST)_PEP_USED for ratio and/or replicates
					}
					elsif ($dispPepType =~ /^msms/) {
						next if (!$msmsStatePos{abs($trueTargetPos)} && !$refDispRatios->{$targetPos} && !$replicatesPos{abs($trueTargetPos)}); # NUM_TRUE_PEP & NUM_PEP_USED (list/export)
#print "*OK_SUP $paramCode TGT:$targetPos<BR>\n" if $protID==$bugProtID;
					}
					else { # usual case
						next if (!$refDispRatios->{$targetPos} && !$linkedRatios{$trueTargetPos}); # allows (NUM|DIST)_PEP_USED
					}
					#next if ($restrictListID && !$restrictList{$protID});
			
					if ($modResStrg) { # quantif of modification: $linkedRatios{$targetPos} could bring additional sites
						#$modResStrg=~s/#//g if $modResStrg=~/^#/; # remove all # if not a multi-modif quantif
						unless ($formattedModRes{$modResStrg}) {
							#@{$formattedModRes{$modResStrg}}=&promsQuantif::formatProteinModificationSites($modResStrg,\%quantifModifInfo,$siteDisplayFormat);
							$formattedModRes{$modResStrg}[0]=&promsQuantif::standardizeSiteCode($modResStrg,$quantifModifInfo{RANK2ID});
							$formattedModRes{$modResStrg}[1]=&promsQuantif::displayModificationSites($formattedModRes{$modResStrg}[0],$quantifModifInfo{DISPLAY},$siteDisplayFormat);
						}
						$protID.='-'.$formattedModRes{$modResStrg}[0];
						next unless $usedIsoforms{$protID}; # $linkedRatios{$targetPos} could bring additional sites
						$refDispModifSites->{$protID}=$formattedModRes{$modResStrg}[1];
					}
				}
				$refQuantifValues->{$protID}{$paramCode}{$targetPos}=$qValue; # if $refDispRatios->{$targetPos};
#print "*VAL #$paramCode: (PROT=$protID, TGT=$targetPos) $qValue*<BR>\n" if $targetPos==1; #$protID eq $bugProtID;
				if ($linkedRatios{$trueTargetPos}) { # SuperRatio
					foreach my $rPos (@{$linkedRatios{$trueTargetPos}}) {
						#$refQuantifValues->{$protID}{$paramCode}{$rPos}+=$qValue;
						##$refQuantifValues->{$protID}{$paramCode}{$rPos}=$qValue if ($refQuantifValues->{$protID}{$paramCode} && $refQuantifValues->{$protID}{$paramCode}{$rPos} && $refQuantifValues->{$protID}{$paramCode}{$rPos} > $qValue); # keep smallest value
#if ($numPepScope eq 'replicate') {
#	next unless $replicatesPos{abs($rPos)};
#}
#else {
						$refQuantifValues->{$protID}{$paramCode}{$rPos}=$qValue if ($refDispRatios->{$rPos} && (!$refQuantifValues->{$protID}{$paramCode} || !$refQuantifValues->{$protID}{$paramCode}{$rPos} || $refQuantifValues->{$protID}{$paramCode}{$rPos} > $qValue)); # keep smallest value
						$refQuantifValues->{$protID}{$paramCode}{-$rPos}=$qValue if ($refDispRatios->{-$rPos} && (!$refQuantifValues->{$protID}{$paramCode} || !$refQuantifValues->{$protID}{$paramCode}{-$rPos} || $refQuantifValues->{$protID}{$paramCode}{-$rPos} > $qValue)); # keep smallest value
#}
					}
				}
			}
		}
	}

	my $refDispSets=($ratioType eq 'None')? $refDispStates : $refDispRatios;
	foreach my $paramCode (@pepParams) { # empty for SuperRatio in graph view
		foreach my $refData (@{$allQuantifData{$paramCode}}) {
			my ($protID,$modResStrg,$targetPos,$qValue)=@{$refData};
			next unless $refProteinInfo->{$protID}; # NUM_PEP_TOTAL has no ratioPos & could bring additional protID
			#next if ($restrictListID && !$restrictList{$protID});
			if ($modResStrg) { # quantif of modification
				#$modResStrg=~s/#//g if $modResStrg=~/^#/; # remove all # if not a multi-modif quantif
				unless ($formattedModRes{$modResStrg}) {
					#@{$formattedModRes{$modResStrg}}=&promsQuantif::formatProteinModificationSites($modResStrg,\%quantifModifInfo,$siteDisplayFormat);
					$formattedModRes{$modResStrg}[0]=&promsQuantif::standardizeSiteCode($modResStrg,$quantifModifInfo{RANK2ID});
					$formattedModRes{$modResStrg}[1]=&promsQuantif::displayModificationSites($formattedModRes{$modResStrg}[0],$quantifModifInfo{DISPLAY},$siteDisplayFormat);
				}
				$protID.='-'.$formattedModRes{$modResStrg}[0];
				next unless $usedIsoforms{$protID};
				$refDispModifSites->{$protID}=$formattedModRes{$modResStrg}[1];
			}
			#next if ($targetPos && !$refDispRatios->{$targetPos}); # no ratioPos for NUM_PEP_TOTAL
			if ($targetPos) { # no ratioPos for NUM_PEP_TOTAL or MQ intensities
				if ($ratioType eq 'None') { # no-ratio MaxQuant
					next unless $refDispStates->{$targetPos};
				}
				else {
					if ($numPepScope =~ /^replicate/) {
						next unless $replicatesPos{abs($targetPos)}; # abs to allow neg tgt pos from Inf ratios
					}
					elsif ($dispPepType =~ /^msms/) {
#print "?_PEP $paramCode TGT:$targetPos<BR>\n" if $protID eq $bugProtID;
						next unless $msmsStatePos{abs($targetPos)};
#print "*OK_PEP $paramCode TGT:$targetPos<BR>\n" if $protID eq $bugProtID;
					}
					else {
						if ($refDispRatios->{-$targetPos}) {$targetPos=-$targetPos;} # reverse ratio management
						elsif (!$refDispRatios->{$targetPos}) {next;}
					}
				}
				$refQuantifValues->{$protID}{$paramCode}{$targetPos}=$qValue if $refDispSets->{$targetPos}; # NUM_PEP_TOTAL has no ratioPos & could bring additional sites
			}
			else { # undef for quantified proteins
				foreach my $tPos (keys %{$refDispSets}) { # extend to all displayed ratios ($targetPos not always defined for single-ratio quantif)
					$refQuantifValues->{$protID}{$paramCode}{$tPos}=$qValue;
				}
			}
		}
	}
		
	####<Compute 2ndary peptide counts for selected scope (ratio/state) based on replicate/msms (Super|Simple)Ratio only)>####
	if ($ratioType=~/S\w+Ratio/ && ($numPepScope =~ /^replicate/ || $dispPepType =~ /^msms/)) {
		my $usedNumPepCode=$numPepCode.'_SUBSET'; # to allow ratio-level computation without overwriting true ratio numPepCode
		foreach my $modProtID (keys %{$refQuantifValues}) {
			foreach my $ratioPos (keys %{$refDispRatios}) {
				next unless $refQuantifValues->{$modProtID}{RATIO}{$ratioPos}; # no ratio in this ratioPos
				my @stateIdx=($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} >= $MAX_INF_RATIO)? (1) : ($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} <= $MIN_INF_RATIO)? (0) : (0,1); # REF,TEST
				my $trueInfRatio=($refQuantifValues->{$modProtID}{RATIO}{$ratioPos}==0.001 || $refQuantifValues->{$modProtID}{RATIO}{$ratioPos}==1000)? 1 : 0;
				#my (%minNumPep,%okState,%numAnyPepState);
				my %minNumPep=(0=>0,1=>0); # initialized because summed when $numPepScope eq 'ratio' even for +/-INF ratios
				my (%okState,%numAnyPepState);
				#my $numPep4ratio=0; # for $dispPepType='msmsXXX' and $numPepScope=ratio
# print "ALL_ST='@stateIdx' Ratio=$refQuantifValues->{$modProtID}{RATIO}{$ratioPos}<BR>" if $modProtID eq $bugProtID;
				#my $tgtPosSwap=($trueInfRatio)? -1 : 1; # Inf ratio ALWAYS takes negative tgt pos since 21/09/20

				###<State-level computation>###
				foreach my $s (@stateIdx) {
					# $minNumPep{$s}=0;
					# $numAnyPepState{$s}=0;
					my $numOkReplic=0; # state for msmsXXX
					foreach my $i (0..$#{$refDispRatios->{$ratioPos}[$s]}) { # replicate: list replicates in state / msms: state itself ($i=0 only)
						my $targetPos=$refDispRatios->{$ratioPos}[$s][$i];
						my $tgtPosSwap=($trueInfRatio && defined $refQuantifValues->{$modProtID}{$numPepCode}{-$targetPos})? -1 : 1; # if Inf ratio: take negative tgt pos if exists / defined because msms can be 0!!!
						$targetPos*=$tgtPosSwap;
# print "ST_R:$s idx:$i SWAP=$tgtPosSwap => TGT=$targetPos ($numPepCode: ",($refQuantifValues->{$modProtID}{$numPepCode}{$targetPos} || 0),"  or ",($refQuantifValues->{$modProtID}{$numPepCode}{-$targetPos} || 0),"<BR>" if $modProtID eq $bugProtID;
						if ($numPepScope =~ /^replicate/) {
							if ($refQuantifValues->{$modProtID}{$numPepCode}{$targetPos} && $refQuantifValues->{$modProtID}{$numPepCode}{$targetPos} >= $dispNumPep) {
								$minNumPep{$s}=$refQuantifValues->{$modProtID}{$numPepCode}{$targetPos} if (!$minNumPep{$s} || $minNumPep{$s} > $refQuantifValues->{$modProtID}{$numPepCode}{$targetPos}); # keep smallest non-0
								$numOkReplic++;
							}
						}
						elsif (defined $refQuantifValues->{$modProtID}{$numPepCode}{$targetPos}) { # msmsXXX
							$minNumPep{$s}=$refQuantifValues->{$modProtID}{$numPepCode}{$targetPos};
							#$numPep4ratio+=$minNumPep{$s};
							foreach my $replicTgPos (@{$refDispRatios->{$ratioPos}[$s+2]}) { # s+2: [0/1]->[2/3] # @{$refDispRatios->MsmsReplic{$ratioPos}[$s]}
								#$numAnyPepState{$s}+=$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$replicTgPos*$tgtPosSwap} if $refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$replicTgPos*$tgtPosSwap};
								if ($trueInfRatio) {
									if ($refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{-$replicTgPos}) { #*$tgtPosSwap
										$numAnyPepState{$s}+=$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{-$replicTgPos}; #*$tgtPosSwap
# print "OK $replicTgPos*swap<BR>" if $modProtID eq $bugProtID;
									}
# 	else { # Strange case for +/-INF where NUM_TRUE_USED assigned to tgtPos>0 and replicate NUM_PEP to tgtPos<0
# 		my $newReplicTgPos=($tgtPosSwap==1)? -$replicTgPos : $replicTgPos;
# 		if ($refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$newReplicTgPos}) {
# 			$numAnyPepState{$s}+=$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$newReplicTgPos};
# print "OK $replicTgPos -> $newReplicTgPos<BR>" if $modProtID eq $bugProtID;
# 		}
# 		else {
# print "BAD $replicTgPos, MSMS_TGT=$targetPos (VAL=$refQuantifValues->{$modProtID}{$numPepCode}{$targetPos})<BR>" if $modProtID eq $bugProtID;
# 		}
# 	}
								}
								elsif ($refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$replicTgPos}) { # normal ratio
									$numAnyPepState{$s}+=$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$replicTgPos};
								}
# print "*$s: $replicTgPos: ",$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$replicTgPos*$tgtPosSwap} || 0,"*<BR>\n" if $modProtID eq $bugProtID;
							}
# print "BUG PROT=$modProtID (TGT_SWAP=$tgtPosSwap)<BR>\n" if (!$numAnyPepState{$s} && $numPepScope ne 'ratio' && $dispPepType eq 'msmsPc');
							$refQuantifValues->{$modProtID}{'FRAC_TRUE_USED'}{$targetPos}=$minNumPep{$s}/$numAnyPepState{$s} if ($numPepScope ne 'ratio' && $dispPepType eq 'msmsPc' && $finalView eq 'list'); # state-level fraction (needed only for display in list)
						}
					}
					if ($numPepScope =~ /^replicate/) {
						$okState{$s}=1 if ($numOkReplic >= $dispNumReplic || $numOkReplic == scalar @{$refDispRatios->{$ratioPos}[$s]}); # filter cannot be > number replicates
					}
				}

				###<Ratio-level computation>###
# $refQuantifValues->{$modProtID}{'NUM_TRUE_USED_RATIO'}{$ratioPos}=$minNumPep{0}+$minNumPep{1} if $dispPepType =~ /^msms/;
				if ($numPepScope eq 'ratio') { # msms with scope=ratio
					#$refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos}=$numPep4ratio;
					$refQuantifValues->{$modProtID}{$usedNumPepCode}{$ratioPos}=$minNumPep{0}+$minNumPep{1};
					$refQuantifValues->{$modProtID}{'FRAC_TRUE_USED'}{$ratioPos}=($minNumPep{0}+$minNumPep{1})/$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$ratioPos} if $dispPepType eq 'msmsPc'; # ($numAnyPepState{0}+$numAnyPepState{1})
#print "Frac: {$modProtID}{$numPepCode}{$ratioPos}=",$refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos},"<BR>\n";
				}
				elsif ($numPepScope =~ /^replicate(Both|One)$/) {
					if (scalar @stateIdx==2) { # normal ratio
						if ($numPepScope eq 'replicateBoth') {
							$refQuantifValues->{$modProtID}{$usedNumPepCode}{$ratioPos}=(!$okState{0} || !$okState{1})? 0 : ($minNumPep{0} <= $minNumPep{1})? $minNumPep{0} : $minNumPep{1}; # keep smallest
						}
						else { # replicateOne
							$refQuantifValues->{$modProtID}{$usedNumPepCode}{$ratioPos}=(!$okState{0} && !$okState{1})? 0 : (!$okState{0})? $minNumPep{1} : (!$okState{1})? $minNumPep{0} : ($minNumPep{0} <= $minNumPep{1})? $minNumPep{1} : $minNumPep{0}; # keep largest if both states pass threshold
						}
					}
					else { # +/-INF ratio: Use state where protein is found
						$refQuantifValues->{$modProtID}{$usedNumPepCode}{$ratioPos}=($okState{$stateIdx[0]})? $minNumPep{$stateIdx[0]} : 0;
					}
				}
				elsif ($numPepScope =~ /^replicate(Ref|Test)$/) { 
					my $usedStateIdx=($numPepScope eq 'replicateRef')? 0 : 1;
					$refQuantifValues->{$modProtID}{$usedNumPepCode}{$ratioPos}=($okState{$usedStateIdx})? $minNumPep{$usedStateIdx} : 0;
				}
				elsif ($numPepScope =~ /^msms(Both|One)$/) {
					if (scalar @stateIdx==2) { # normal ratio
						my $usedStateIdx;
						if ($dispPepType eq 'msmsNum') {
							if ($numPepScope eq 'msmsBoth') {
								$usedStateIdx=($minNumPep{0} <= $minNumPep{1})? 0 : 1; # keep smallest
							}
							else { # msmsOne
								$usedStateIdx=($minNumPep{0} <= $minNumPep{1})? 1 : 0; # keep largest
							}
						}
						else { # msmsPc
# print "ERROR NON-INF: prot=$modProtID ($numAnyPepState{0},$numAnyPepState{1})<BR>\n" if (!$numAnyPepState{0} || !$numAnyPepState{1});
							my @fracPep=($minNumPep{0}/$numAnyPepState{0},$minNumPep{1}/$numAnyPepState{1});
							#my @fracPep=(($numAnyPepState{0})? $minNumPep{0}/$numAnyPepState{0} : 0,($numAnyPepState{1})? $minNumPep{1}/$numAnyPepState{1} : 0);
							if ($numPepScope eq 'msmsBoth') {
								$usedStateIdx=($fracPep[0] <= $fracPep[1])? 0 : 1; # keep smallest
							}
							else { # msmsOne
								$usedStateIdx=($fracPep[0] <= $fracPep[1])? 1 : 0; # keep largest
							}
							$refQuantifValues->{$modProtID}{'FRAC_TRUE_USED'}{$ratioPos}=$fracPep[$usedStateIdx];
						}
						$refQuantifValues->{$modProtID}{$usedNumPepCode}{$ratioPos}=$minNumPep{$usedStateIdx};
						# $refQuantifValues->{$modProtID}{'NUM_PEP_USED_STATE'}{$ratioPos}=$numAnyPepState{$usedStateIdx}; # WARNING: overwrites default ratio-wise NUM_PEP_USED with state-wise NUM_PEP_USED
					}
					else { # +/-INF ratio: Use state where protein is found
						$refQuantifValues->{$modProtID}{$usedNumPepCode}{$ratioPos}=$minNumPep{$stateIdx[0]}; # $stateIdx[0]=0 or 1
# print "BAD! P=$modProtID, R=$ratioPos, V=$minNumPep{$stateIdx[0]} ($stateIdx[0])<BR>\n" if !$numAnyPepState{$stateIdx[0]};
						if ($dispPepType eq 'msmsPc') {$refQuantifValues->{$modProtID}{'FRAC_TRUE_USED'}{$ratioPos}=($numAnyPepState{$stateIdx[0]})? $minNumPep{$stateIdx[0]}/$numAnyPepState{$stateIdx[0]} : 0;} # 0 => major error in data => protein will be skipped by peptide filtering
# print "P=$modProtID, $minNumPep{$stateIdx[0]}/$numAnyPepState{$stateIdx[0]} = $refQuantifValues->{$modProtID}{FRAC_TRUE_USED}{$ratioPos}<BR>\n" if ($modProtID eq "$bugProtID" && $dispPepType eq 'msmsPc');
						# $refQuantifValues->{$modProtID}{'NUM_PEP_USED_STATE'}{$ratioPos}=$numAnyPepState{$stateIdx[0]}; # WARNING: overwrites default ratio-wise NUM_PEP_USED with state-wise NUM_PEP_USED (just to safe)
					}
				}
				elsif ($numPepScope =~ /^msms(Ref|Test)$/) { 
					my $usedStateIdx=($numPepScope eq 'msmsRef')? 0 : 1;
					$refQuantifValues->{$modProtID}{$usedNumPepCode}{$ratioPos}=$minNumPep{$usedStateIdx} || 0; # undef if +/-INF
					if ($dispPepType eq 'msmsPc') {$refQuantifValues->{$modProtID}{'FRAC_TRUE_USED'}{$ratioPos}=($numAnyPepState{$usedStateIdx})? $minNumPep{$usedStateIdx}/$numAnyPepState{$usedStateIdx} : 0;} # 0 => major error in data => protein will be skipped by peptide filtering
					# $refQuantifValues->{$modProtID}{'NUM_PEP_USED_STATE'}{$ratioPos}=$numAnyPepState{$usedStateIdx}; # WARNING: overwrites default ratio-wise NUM_PEP_USED with state-wise NUM_PEP_USED
				}
			}
		}
	}
	elsif ($ratioType eq 'None' && $dispPepType eq 'msmsPc') { # compute 'FRAC_TRUE_USED' for ABUNDANCE
		foreach my $modProtID (keys %{$refQuantifValues}) {
			foreach my $targetPos (keys %{$refQuantifValues->{$modProtID}{'NUM_TRUE_USED'}}) {	
				$refQuantifValues->{$modProtID}{'FRAC_TRUE_USED'}{$targetPos}=$refQuantifValues->{$modProtID}{'NUM_TRUE_USED'}{$targetPos}/$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$targetPos};
			}
		}
	}


	###>Fetching protein info
	if (scalar keys %{$refProteinInfo}) {
		my $dbh=&promsConfig::dbConnect;
		my ($GNidentID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
		my ($allAnaStrg)=$dbh->selectrow_array("SELECT GROUP_CONCAT(ID_ANALYSIS) FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
		my ($geneQuery1,$geneQuery2)=($finalView eq 'list' || $action eq 'export')? (",GROUP_CONCAT(DISTINCT MI.VALUE ORDER BY IDENT_RANK SEPARATOR ',')","LEFT JOIN MASTERPROT_IDENTIFIER MI ON P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND MI.ID_IDENTIFIER=$GNidentID") : ('','');
		my ($anaQuery1,$anaQuery2)=($finalView eq 'list')? (",GROUP_CONCAT(DISTINCT ID_ANALYSIS SEPARATOR ',')","LEFT JOIN ANALYSIS_PROTEIN AP ON AP.ID_PROTEIN=P.ID_PROTEIN AND ID_ANALYSIS IN ($allAnaStrg)") : ('',''); #  "LEFT JOIN" to handle rare MaxQuant bug where protein is quantified but no identified in selected Analysi(e)s
		#my $protIdStrg=join(',',keys %{$refProteinInfo});
		my @protIdList=keys %{$refProteinInfo};
		while (my @subProtList=splice(@protIdList,0,1000)) { # chuncks of 1000 proteins
			my $protIdStrg=join(',',@subProtList);
			my $sthProtInfo=$dbh->prepare("SELECT P.ID_PROTEIN,ALIAS,MW,PROT_DES,ORGANISM,PROT_LENGTH$geneQuery1$anaQuery1
											FROM PROTEIN P
											$geneQuery2
											$anaQuery2
											WHERE P.ID_PROTEIN IN ($protIdStrg) GROUP BY P.ID_PROTEIN"); # P.ID_PROJECT=$projectID
			$sthProtInfo->execute;
			while (my ($protID,$alias,$mw,$protDes,$organism,$protLength,$geneStrg,$anaStrg)=$sthProtInfo->fetchrow_array) {
				next unless $refProteinInfo->{$protID}; # just to be safe
				$alias=~s/ .*//; # clean alias from badly parsed characters (identifier conversion)
				$mw=($mw)? 1*(sprintf "%.1f",$mw/1000) : 0;
				my @geneList=split(',',$geneStrg) if $geneStrg;
				if ($finalView eq 'list') {$anaStrg=0 unless $anaStrg;}
				else {$anaStrg=$allAnaStrg;}
				@{$refProteinInfo->{$protID}}=($alias,$mw,$protDes,$organism,$protLength,\@geneList,$anaStrg); # [5] for genelist
			}	
			$sthProtInfo->finish;
		}
		$dbh->disconnect;	
	}

	return ($pvalueCode,$minPvalue);

}

sub generateRatioCodes { # WARNING: Computes all possible ratios even if "All vs State1"!!!
	my ($labelType,$ratioType,$numStates,$isPepIntensity,$algoVersion)=@_;
	my %ratioCode;
	if ($ratioType=~/S\w+Ratio/) {
		if ($isPepIntensity) {
			my $ratioPos=0;
			foreach my $refStatePos (1..$numStates-1) {
				foreach my $testStatePos ($refStatePos+1..$numStates) {
					next if $testStatePos==$refStatePos;
					$ratioPos++;
					$ratioCode{$ratioPos}=($algoVersion>=3)? "State$testStatePos.State$refStatePos" : "State$testStatePos-State$refStatePos";
					$ratioCode{-$ratioPos}='-'.$ratioCode{$ratioPos};
				}
			}
		}
		elsif ($algoVersion >= 3) {
			my $numPrimRatios=$numStates-1;
			foreach my $ratioPos (1..$numPrimRatios) {
				$ratioCode{$ratioPos}='State'.($ratioPos+1).'.State1';
				$ratioCode{-$ratioPos}='-'.$ratioCode{$ratioPos};
			}
			if ($ratioType eq 'SuperRatio') {
				my $ratioPos=$numPrimRatios;
				foreach my $x (2..$numPrimRatios) { # ratios of ratios
					foreach my $y ($x+1..$numPrimRatios+1){
						next if $y==$x;
						$ratioPos++;
						$ratioCode{$ratioPos}="State$y.State1.State$x.State1";
						$ratioCode{-$ratioPos}='-'.$ratioCode{$ratioPos};
					}
				}
			}
		}
		else {
			if ($ratioType eq 'SuperRatio') {
				my $numPrimRatios=$numStates-1;
				foreach my $ratioPos (1..$numPrimRatios) {
					$ratioCode{$ratioPos}=chr(64+$ratioPos);
					$ratioCode{-$ratioPos}='-'.$ratioCode{$ratioPos};
				}
				my $ratioPos=$numPrimRatios;
				foreach my $x (1..$numPrimRatios-1) { # ratios of ratios
					foreach my $y ($x+1..$numPrimRatios){
						next if $y==$x;
						$ratioPos++;
						$ratioCode{$ratioPos}=chr(64+$y).'-'.chr(64+$x);
						$ratioCode{-$ratioPos}='-'.$ratioCode{$ratioPos};
					}
				}
			}
			else {
				my $numPrimRatios=$numStates*($numStates-1)/2;
				foreach my $ratioPos (1..$numPrimRatios) {
					$ratioCode{$ratioPos}=chr(64+$ratioPos);
					$ratioCode{-$ratioPos}='-'.$ratioCode{$ratioPos};
				}
			}
		}
	}
	else {
		foreach my $ratioPos (1..$numStates*($numStates-1)/2) {$ratioCode{$ratioPos}=$ratioPos; $ratioCode{-$ratioPos}=-$ratioPos;}
	}
	return %ratioCode;
}

sub getReplicateTargetPos { # NO GLOBALS
	# Adds to @{$refDispRatios->{$ratioPos}} the list of replicPos from both test/ref states
	# Fills %{$refReplicatesPos}: replicPos => state
	# NO NEGATIVE REPLIC POS HERE
	my ($numStates,$refDispRatios,$refRatioCode,$refLabelingInfo,$refStateInfo,$refReplicatesPos)=@_;
	my $usedTargetPos=scalar @{$refLabelingInfo->{'RATIOS'}};
	$usedTargetPos+=$numStates if ($refLabelingInfo->{'MEAN_STATE'} && $refLabelingInfo->{'MEAN_STATE'}[0]==1);

	foreach my $ratioPos (keys %{$refDispRatios}) {
		my $trueRatioPos=abs($ratioPos); # in case of reversed ratio
		my %usedStates;
		while ($refRatioCode->{$trueRatioPos}=~/State(\d+)/g) {$usedStates{$1}=1;}
		delete $usedStates{'1'} if scalar keys %usedStates > 2; # SuperRatio: skip State1
		foreach my $statePos (sort{if ($ratioPos < 0) {$b<=>$a} else {$a<=>$b}} keys %usedStates) { # order REF then TEST
			my $currTargetPos=$usedTargetPos;
			foreach my $prevStatePos (1..$statePos-1) {$currTargetPos+=$refStateInfo->{$prevStatePos}{'NUM_REPLIC'};}
			foreach my $replicPos ($currTargetPos+1..$currTargetPos+$refStateInfo->{$statePos}{'NUM_REPLIC'}) {
				#push @{$refDispRatios->{$ratioPos}},$replicPos;
				$refReplicatesPos->{$replicPos}=$statePos;
			}
			push @{$refDispRatios->{$ratioPos}},[$currTargetPos+1..$currTargetPos+$refStateInfo->{$statePos}{'NUM_REPLIC'}]; # REF then TEST
		}
#print "R=$ratioPos: REF0=@{$refDispRatios->{$ratioPos}[0]} TEST1=@{$refDispRatios->{$ratioPos}[1]}<BR>\n"; # last 2 indexes for replicate pos
	}
}

sub getMsmsStateTargetPos { # NO GLOBALS
	# Adds to @{$refDispRatios->{$ratioPos}} the 2 test/ref states
	# Fills %{$refReplicatesPos}: replicPos => state
	# NO NEGATIVE REPLIC POS HERE
	my ($refDispRatios,$refRatioCode,$refLabelingInfo,$refMsmsStatePos)=@_;
	$refMsmsStatePos={} unless $refMsmsStatePos;
	my $numRatios=scalar @{$refLabelingInfo->{'RATIOS'}};

	foreach my $ratioPos (keys %{$refDispRatios}) {
		my $trueRatioPos=abs($ratioPos); # in case of reversed ratio
		my %usedStates;
		while ($refRatioCode->{$trueRatioPos}=~/State(\d+)/g) {$usedStates{$1}=1;}
		delete $usedStates{'1'} if scalar keys %usedStates > 2; # SuperRatio: skip State1
		my @stTgtPos=(sort{if ($ratioPos < 0) {$b<=>$a} else {$a<=>$b}} keys %usedStates); # order REF then TEST
		push @{$refDispRatios->{$ratioPos}},([$numRatios+$stTgtPos[0]],[$numRatios+$stTgtPos[1]]); # same structure than Replicate
#print "MSMS: RPOS=$ratioPos, STPOS=",$numRatios+$stTgtPos[0],",",$numRatios+$stTgtPos[1],"<BR>";
		$refMsmsStatePos->{$numRatios+$stTgtPos[0]}=$stTgtPos[0];
		$refMsmsStatePos->{$numRatios+$stTgtPos[1]}=$stTgtPos[1];
	}
}

sub printQuantificationSummary { # GLOBALS: $selQuantifID, $view, $projectID, $designID, $ratioType, $lightColor, $darkColor, $quantifType, $quantifSoftware, $isPepIntensity
	my ($dbh,$anaID,$refInfo,$refGoAnalyses)=@_;
	#my ($lightColor,$darkColor)=&promsConfig::getRowColors;
	my %proteinQuantifFamilies=&promsQuantif::getProteinQuantifFamilies;
	my $resultDir="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/results";
	my $graphDirHtml="$promsPath{quanti_html}/project_$projectID/quanti_$selQuantifID/results/graph";
	my $protNormDirHtml="$promsPath{quanti_html}/project_$projectID/quanti_$selQuantifID/results/resultsNorm"; # protein-level bias correction only
	my ($newLine,$startB,$endB)=($action eq 'export')? ("\n",'','') : ('<BR>','<B>','</B>');

	my $ratioTypeStrg=($ratioType eq 'Ratio')? 'Simple ratios (old version)' : ($ratioType eq 'SuperRatio')? 'Super ratios' : ($ratioType eq 'None')? 'No ratio' : 'Simple ratios';
	$ratioTypeStrg=$startB.$ratioTypeStrg.$endB;
	#my $quantifSoftware=($refInfo->{'SOFTWARE'})? $refInfo->{'SOFTWARE'}[0] : 'myProMS';
	#$quantifSoftware='MaxQuant' if $quantifSoftware eq 'MQ';
	my $softNameStrg=($quantifSoftware eq 'myProMS' && $quantifType eq 'PROT_RULER')? 'myProMS-Ruler' : ($quantifSoftware eq 'myProMS')? 'myProMS-Quant' : $quantifSoftware;
	my $softVersionStrg=($algoVersion)? ' v'.$algoVersion : ($refInfo->{'SOFTWARE'} && $refInfo->{'SOFTWARE'}[1])? ' v'.$refInfo->{'SOFTWARE'}[1] : ''; # any Software

	my $topNString='';
	if ($isPepIntensity && !$isModifQuantif && $quantifSoftware =~ /myProMS|Spectronaut/) {
		$topNString=($refInfo->{'MATCH_PEP'} && $refInfo->{'MATCH_PEP'}[0])? ' - '.$startB.'Matching'.$endB : ' - '.$startB.'Any'.$endB;
		$topNString.=($refInfo->{'TOP_N'})? ' top'.$startB.$refInfo->{'TOP_N'}[0].$endB.' peptide' : ' top'.$startB.'N'.$endB.' peptide';
		$topNString.='s' if ($topNString && (!$refInfo->{'TOP_N'} || $refInfo->{'TOP_N'}[0]>1));
	}
	my $swathStrg=($quantifSoftware=~/SWATH/)? ' (SWATH)' : ($quantifSoftware=~/DIA/) ? ' (DIA)' : '';

	#<Conditions & replicates
	my (%stateInfo,%replicateName,%numTechReplicates,%bioRepLabel);
	my ($protRulerTypeStrg, $intMetric, $protRulerNormStrg, $protRulerParamsStrg);  # Initialize variables for proteomic ruler
	#my $existReplicates=0;
	#my $supRatioReplicates=0;
	my $numStates=0;
	###> Get parent quantifType to distinguish from label/label-free quantification
	#my ($parentQuantifType)=$dbh->selectrow_array("SELECT CODE FROM QUANTIFICATION,PARENT_QUANTIFICATION,QUANTIFICATION_METHOD WHERE QUANTIFICATION.ID_QUANTIFICATION=PARENT_QUANTIFICATION.ID_PARENT_QUANTIFICATION AND QUANTIFICATION.ID_QUANTIFICATION_METHOD=QUANTIFICATION_METHOD.ID_QUANTIFICATION_METHOD AND PARENT_QUANTIFICATION.ID_QUANTIFICATION=$selQuantifID LIMIT 1");
	#if ($parentQuantifType && $parentQuantifType eq 'XIC'){#} # Be careful with Pep-Ratio from label vs. label-free experiments... STATES are not formated the same!!!
	if ($designID && $refInfo->{STATES}) {  # refInfo{STATES} undef in Proteomic Ruler case, cannot work (VL 30/07/19)
		my $sthAnaName=$dbh->prepare("SELECT NAME FROM ANALYSIS WHERE ID_ANALYSIS=?");
		my $sthExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
		my $sthPepQuantif=$dbh->prepare("SELECT NAME,QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
		my $sthObsBioSamp=$dbh->prepare("SELECT ID_BIOSAMPLE FROM OBSERVATION WHERE ID_OBSERVATION=?");
		my $sthObsMsSamp=$dbh->prepare("SELECT ID_SAMPLE FROM ANALYSIS A WHERE ID_ANALYSIS=?");
		my $sthBsName=$dbh->prepare("SELECT NAME FROM BIOSAMPLE WHERE ID_BIOSAMPLE=?");
		my $sthMsName=$dbh->prepare("SELECT NAME FROM SAMPLE WHERE ID_SAMPLE=?");
		my (%anaName,%parQuantifInfo);
		foreach my $stateData (@{$refInfo->{STATES}}) {
			$stateData=~s/#//g; # removes id tags
			$numStates++;
			my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$stateData);
			#if ($ratioType eq 'SuperRatio') { #=~/S\w+Ratio/
			#	if ($numStates==1) {$supRatioReplicates=$numBioRep;}
			#	elsif ($numBioRep==1) {$supRatioReplicates--;}
			#}
			#$expCondID=~ s/^.//; # To remove the first character (which is a '#' to show that it is an ID)
			$sthExpCondName->execute($expCondID);
			($stateInfo{$numStates}{'NAME'})=$sthExpCondName->fetchrow_array;
			$condToState{$expCondID}=$numStates;
			$stateInfo{$numStates}{'POPUP'}='';
			#my (%replicateInfo,%replicateStrg);
			my %bioRepHeader; #%channelName
			my ($currBioRep,$allTechRepCount,$numPools,$numObs)=(0,0,0,0);
			foreach my $bioReplicate (split(/\./,$quantiObsIDs)) { # bio rep separator
				$currBioRep++;
				my (%fracAnaNames,%fracAnaNamesPopup,%numPoolObs,%bioRepNames);
				my $numTechRep=0;
				my $numObsInBioRep=0;
				foreach my $techReplicate (split(/&/,$bioReplicate)) { # tech rep separator
					$numTechRep++;
					$numPoolObs{$numTechRep}=0;
					foreach my $pooledObs (split(/\+/,$techReplicate)) {
						$numPoolObs{$numTechRep}++;
						my ($obsID,$parQuantifID,$anaID,$targetPos);
						my @obsData=split(/:/,$pooledObs);
						if (scalar @obsData==4) {($obsID,$parQuantifID,$anaID,$targetPos)=@obsData;} # design post 07/2013 -> obsID,parQuantifID,anaID,targetPos
						else { # older design -> parQuantifID,anaID
							($parQuantifID,$anaID)=@obsData;
							$obsID=0;
							$targetPos=0;
						}
						unless ($anaName{$anaID}) {
							$sthAnaName->execute($anaID);
							($anaName{$anaID})=$sthAnaName->fetchrow_array;
						}
						if ($obsID) {
							$sthObsBioSamp->execute($obsID);
							my ($bsID)=$sthObsBioSamp->fetchrow_array;
							$bioRepNames{BIOSAMPLE}{$bsID}++ if $bsID;
						}
						$sthObsMsSamp->execute($anaID);
						my ($sampID)=$sthObsMsSamp->fetchrow_array;
						$bioRepNames{MSSAMPLE}{$sampID}++;
						
						unless($parQuantifInfo{$parQuantifID}) {
							$sthPepQuantif->execute($parQuantifID);
							my ($name,$qAnnot)=$sthPepQuantif->fetchrow_array;
							$parQuantifInfo{$parQuantifID}{NAME}=$name;
							my @quantifInfo=split(/::/,$qAnnot);
							foreach my $channelInfo (@quantifInfo[1..$#quantifInfo]) {
								my ($chanPos,$chanName)=split(/;/,$channelInfo); # + other unsued elements
								$parQuantifInfo{$parQuantifID}{CHANNELS}{$chanPos}=$chanName;
							}
						}
						my $targetName=($targetPos)? ' > '.$parQuantifInfo{$parQuantifID}{CHANNELS}{$targetPos} : '';
						push @{$fracAnaNames{$numTechRep}},$anaName{$anaID}.$targetName;
						push @{$fracAnaNamesPopup{$numTechRep}},$anaName{$anaID}.' : '.$parQuantifInfo{$parQuantifID}{NAME}.$targetName;
						# my $targetName='';
						# if ($targetPos) {
						# 	unless ($channelName{$parQuantifID}) {
						# 		$sthPepQuantif->execute($parQuantifID);
						# 		my ($qAnnot)=$sthPepQuantif->fetchrow_array;
						# 		my @quantifInfo=split(/::/,$qAnnot);
						# 		foreach my $channelInfo (@quantifInfo[1..$#quantifInfo]) {
						# 			my ($chanPos,$chanName)=split(/;/,$channelInfo); # + other unsued elements
						# 			$channelName{$parQuantifID}{$chanPos}=$chanName;
						# 		}
						# 	}
						# 	$targetName=' > '.$channelName{$parQuantifID}{$targetPos};
						# }
						# push @{$fracAnaNames{$numTechRep}},$anaName{$anaID}.$targetName;
					}
					$numObsInBioRep+=$numPoolObs{$numTechRep};
				}
				$allTechRepCount+=$numTechRep;
				$replicateName{$numStates}{$currBioRep}=($action eq 'export')? $stateInfo{$numStates}{'NAME'} : $startB.$stateInfo{$numStates}{'NAME'}.$endB;
				$numTechReplicates{"$numStates:$currBioRep"}=$numTechRep;

				##<Guess most suitable name for bioRep
				if ($bioRepNames{BIOSAMPLE} && scalar keys %{$bioRepNames{BIOSAMPLE}}==1) { # Bio-sample
					my $bsID=(keys %{$bioRepNames{BIOSAMPLE}})[0];
					#if ($bioRepNames{BIOSAMPLE}{$bsID}==$numObsInBioRep) {
						$sthBsName->execute($bsID);
						($bioRepLabel{"$numStates:$currBioRep"})=$sthBsName->fetchrow_array;
					#}
				}
				if (!$bioRepLabel{"$numStates:$currBioRep"} && scalar keys %{$bioRepNames{MSSAMPLE}}==1) { # MS Sample
					my $sampID=(keys %{$bioRepNames{MSSAMPLE}})[0];
					#if ($bioRepNames{MSSAMPLE}{$sampID}==$numObsInBioRep) {
						$sthMsName->execute($sampID);
						($bioRepLabel{"$numStates:$currBioRep"})=$sthMsName->fetchrow_array;
					#}
				}
				if (!$bioRepLabel{"$numStates:$currBioRep"}) {
					$bioRepLabel{"$numStates:$currBioRep"}='Replicate #'.$currBioRep;
				}

				if ($numObsInBioRep > 1) { # 2 or more fractions and/or tech rep
					if ($numBioRep > 1) {
						$replicateName{$numStates}{$currBioRep}.="$startB:$endB ".$bioRepLabel{"$numStates:$currBioRep"};
						if ($action eq 'export') {$stateInfo{$numStates}{'POPUP'}.=$startB.$bioRepLabel{"$numStates:$currBioRep"}.$endB.$newLine;}
						else {@{$bioRepHeader{$currBioRep}}=($bioRepLabel{"$numStates:$currBioRep"},$numTechRep);}
					}
					foreach my $techRep (1..$numTechRep) {
						$numPools++;
						if ($numTechRep > 1) {
							$stateInfo{$numStates}{'POPUP'}.=($action eq 'export')? "•Tech. rep. #$techRep" : "<TD valign=top nowrap>&bull;<U>Tech. rep. #$techRep";
							$stateInfo{$numStates}{'POPUP'}.=" (Pool)" if $numPoolObs{$techRep} > 1;
							$stateInfo{$numStates}{'POPUP'}.=($action eq 'export')? ":$newLine" : "</U>:$newLine";
						}
						elsif ($numPoolObs{$techRep} > 1) {
							$stateInfo{$numStates}{'POPUP'}.=($action eq 'export')? "•Pool:$newLine" : "<TD valign=top nowrap>&bull;<U>Pool:</U>$newLine";
						}
						if ($numTechRep > 1 || $numPoolObs{$techRep} > 1) {
							$stateInfo{$numStates}{'POPUP'}.=($action eq 'export')?
															"  +".join("$newLine  +",@{$fracAnaNamesPopup{$techRep}}).$newLine
															: "&nbsp;&nbsp;+".join("$newLine&nbsp;&nbsp;+",@{$fracAnaNames{$techRep}})."<BR>";
						}
						else {
							$replicateName{$numStates}{$currBioRep}.="$startB:$endB $fracAnaNames{$techRep}[0]";
							$stateInfo{$numStates}{'POPUP'}.=($action eq 'export')? "•$fracAnaNamesPopup{$techRep}[0]$newLine" : "<TD valign=top nowrap>&bull;$fracAnaNamesPopup{$techRep}[0]<BR>";
						}
					}
				}
				else { # only 1 obs in replicate
					$replicateName{$numStates}{$currBioRep}.="$startB:$endB $fracAnaNames{1}[0]";
					$stateInfo{$numStates}{'POPUP'}.=($action eq 'export')? "•$fracAnaNamesPopup{1}[0]$newLine" : "<TD valign=top nowrap>&bull;$fracAnaNamesPopup{1}[0]<BR>";
				}
				$numObs+=$numObsInBioRep;
				#$stateString.=" $replicateName{$numStates}{$currBioRep},";
			}
			#$stateInfo{$numStates}{'NUM_REPLIC'}=($currBioRep > 1)? $currBioRep : $allTechRepCount; # needed for targetPos selection for numPep filtering
			$stateInfo{$numStates}{'NUM_REPLIC'}=( ($refInfo->{'RESIDUAL_VAR'} && $refInfo->{'RESIDUAL_VAR'}[0] eq 'biological') || (!$refInfo->{'RESIDUAL_VAR'} && $currBioRep > 1) )? $currBioRep : $allTechRepCount; # needed for targetPos selection for numPep filtering
			#chop($stateString);
			#$stateInfo{$numStates}{'NAME_LONG'}=$stateString;
			$stateInfo{$numStates}{'NAME_LONG'}=($action eq 'export')? ': ' : ':&nbsp;';
			if ($numBioRep>1) {
				#$existReplicates=1;
				$stateInfo{$numStates}{'REPLIC'}=($action eq 'export' || $ratioType=~/S\w+Ratio/ || $quantifType=~/^(MQ|PROT_ABUNDANCE)$/)? '' : "&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Correlation\" onclick=\"displayReplicCorrelation(this,$numStates,'$stateInfo{$numStates}{NAME}')\">"; # not true replicates if SuperRatio
				$stateInfo{$numStates}{'NAME_LONG'}.=($ratioType eq 'SuperRatio' && $numStates==1 && $numBioRep==scalar @{$refInfo->{STATES}})? "$numBioRep sets / " : "$numBioRep bio. rep. / ";
			}
			else {$stateInfo{$numStates}{'REPLIC'}='';}
			if ($allTechRepCount > $numBioRep) {
				$stateInfo{$numStates}{'NAME_LONG'}.="$allTechRepCount tech. rep. / ";
			}
			elsif ($numPools) {
				$stateInfo{$numStates}{'NAME_LONG'}.=($numPools > 1)? "$numPools pools / " : "$numPools pool / ";
			}
			$stateInfo{$numStates}{'NAME_LONG'}.="$numObs observation";
			$stateInfo{$numStates}{'NAME_LONG'}.='s' if $numObs > 1;

			if (scalar keys %bioRepHeader) {
				my $popupStrg='<TABLE border=1 cellspacing=0><TR>';
				foreach my $currBioRep (sort{$a<=>$b} keys %bioRepHeader) {
					$popupStrg.="<TH colspan=$bioRepHeader{$currBioRep}[1]>$bioRepHeader{$currBioRep}[0]</TH>";
				}
				$popupStrg.='</TR><TR>';
				$stateInfo{$numStates}{'POPUP'}=$popupStrg.'</TR><TR>'.$stateInfo{$numStates}{'POPUP'}.'</TR></TABLE>';
			}
#print "#$numStates: $numBioRep, $supRatioReplicates<BR>\n";
		}
		$sthAnaName->finish;
		$sthExpCondName->finish;
		$sthPepQuantif->finish;
		$sthObsBioSamp->finish;
		$sthObsMsSamp->finish;
		$sthBsName->finish;
		$sthMsName->finish;
	}
	elsif ($designID && $protRulerCase) {  # Proteomic Ruler case (no states)
		my $sthParentQuantifName=$dbh->prepare("SELECT NAME FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
		my $sthParentStateName = $dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
		my $sthQuantifAnnot=$dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
		my $sthAnaInfo = $dbh->prepare("SELECT A.NAME, A.ID_ANALYSIS FROM ANALYSIS A 
										INNER JOIN OBSERVATION O ON A.ID_ANALYSIS=O.ID_ANALYSIS
										INNER JOIN OBS_EXPCONDITION OE ON O.ID_OBSERVATION=OE.ID_OBSERVATION
										INNER JOIN EXPCONDITION E ON OE.ID_EXPCONDITION=E.ID_EXPCONDITION
										WHERE E.ID_EXPCONDITION=? ORDER BY A.DISPLAY_POS ASC");
		
		my (%parentQuantifName,%anaName);
		my %grouping;
		if ($refInfo->{'AVG_MODE'}[0] eq "GROUPS_SAME_NORM") {
			foreach my $groupData (@{$refInfo->{GROUPING}}) {
				$groupData=~s/#//g; # removes id tags	
				my ($groupQuantif,$groupNb)=split(':',$groupData);
				$grouping{$groupQuantif}=$groupNb;
			}
		}
		if ($refInfo->{'AVG_MODE'}[0] eq "AVG_ALL") {
			$stateInfo{0}{'NAME'}="Average on all states";
			$stateInfo{0}{'POPUP'}='';
			$replicateName{0}{0}.="$startB:$endB $stateInfo{0}{'NAME'}";
				
			foreach my $stateData (@{$refInfo->{PARENT_Q}}) {
				$stateData=~s/#//g; # removes id tags
				my ($parentQuantif,$targetPos)=split(':',$stateData);
				my ($parentQuantifID, $parentStateID, $parentPos) = split(',', $parentQuantif);
				unless ($parentQuantifName{$parentQuantifID}) {
					$sthParentQuantifName->execute($parentQuantifID);
					($parentQuantifName{$parentQuantifID}) = $sthParentQuantifName->fetchrow_array;
				}
				$sthParentStateName->execute($parentStateID);
				my ($parentStateName) = $sthParentStateName->fetchrow_array;
				my $stateName = "$parentQuantifName{$parentQuantifID} - $parentStateName";
				$sthAnaInfo->execute($parentStateID);
				$stateInfo{0}{'POPUP'} .= ($action eq 'export')? $stateName.$newLine :
										   "<TD valign=top nowrap>$stateName"."$newLine";

				while (my ($analysisName, $analysisID) = $sthAnaInfo->fetchrow_array) {
					unless ($anaName{$analysisID}) {
						$anaName{$analysisID} = $analysisName;
					}
					$stateInfo{0}{'POPUP'} .= ($action eq 'export')? "•$anaName{$analysisID}$newLine" :
											   "<TD valign=top nowrap>&bull;$anaName{$analysisID}$newLine";
				}
			}
		}
		else {
			foreach my $stateData (@{$refInfo->{PARENT_Q}}) {
				$stateData=~s/#//g; # removes id tags
				$numStates++;
				my ($parentQuantif,$targetPos)=split(':',$stateData);
				my ($parentQuantifID, $parentStateID, $parentPos) = split(',', $parentQuantif);
				unless ($parentQuantifName{$parentQuantifID}) {
					$sthParentQuantifName->execute($parentQuantifID);
					($parentQuantifName{$parentQuantifID}) = $sthParentQuantifName->fetchrow_array;
				}
				$sthParentStateName->execute($parentStateID);
				my ($parentStateName) = $sthParentStateName->fetchrow_array;
				
				my $stateName = "$parentQuantifName{$parentQuantifID} - $parentStateName";
				$stateInfo{$numStates}{'NAME'} = $stateName;
				
				$stateInfo{$numStates}{'POPUP'}='';
				my $currBioRep=0;
				$replicateName{$numStates}{$currBioRep}.="$startB:$endB $stateInfo{$numStates}{'NAME'}";

				$sthAnaInfo->execute($parentStateID);
				while (my ($analysisName, $analysisID) = $sthAnaInfo->fetchrow_array) {
					unless ($anaName{$analysisID}) {
						$anaName{$analysisID} = $analysisName;
					}
					$stateInfo{$numStates}{'POPUP'} .= ($action eq 'export')? "•$anaName{$analysisID}$newLine" :
													   "<TD valign=top nowrap>&bull;$anaName{$analysisID}$newLine";
				}
				$stateInfo{$numStates}{'REPLIC'}='';
				if (%grouping) {
					$stateInfo{$numStates}{'NAME_LONG'}.=($action eq 'export')? ": Group $grouping{$parentQuantif}" : ":&nbsp;Group&nbsp;$grouping{$parentQuantif}";
				}
			}
		}
		$sthParentQuantifName->finish;
		$sthParentStateName->finish;
		$sthQuantifAnnot->finish;
		$sthAnaInfo->finish;

		# Set additional information strings specific to proteomic ruler
		$protRulerParamsStrg = "Total cellular protein concentration: $refInfo->{'NUM_PARAMS'}[0] g/L";
		if ($refInfo->{'RULER_TYPE'}[0]==0){
			$protRulerTypeStrg="Total protein amount";
			$protRulerParamsStrg.=" - Total protein quantity per sample: $refInfo->{'NUM_PARAMS'}[1] pg";
		} elsif ($refInfo->{'RULER_TYPE'}[0]==1){
			$protRulerTypeStrg="Histone Proteomic Ruler";
			$protRulerParamsStrg.=" - Ploidy of cell type: $refInfo->{'NUM_PARAMS'}[2]";
		} elsif ($refInfo->{'RULER_TYPE'}[0]==2){
			$protRulerTypeStrg="Custom Proteins Ruler";
			$protRulerParamsStrg.="&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Display Custom Proteins\" style=\"width:150px\" onclick=\"displayCustomProts(this)\">&nbsp;" unless ($action eq 'export');
		}
		if ($refInfo->{'METRIC'}) {
			# Need to fetch DB because METRIC is a quantif parameter from parent quantification method,
			# so it is not present in %quantifParamInfo (results in empty string)
			($intMetric)=$dbh->selectrow_array("SELECT NAME FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIF_PARAMETER=$refInfo->{'METRIC'}[0]");
		}
		if ($refInfo->{'AVG_MODE'}) {
			if ($refInfo->{'AVG_MODE'}[0] eq "ALL_SEPARATE") {
				$protRulerNormStrg="All samples separately";
			} elsif ($refInfo->{'AVG_MODE'}[0] eq "ALL_SAME_NORM") {
				$protRulerNormStrg="Same normalization for all";
			} elsif ($refInfo->{'AVG_MODE'}[0] eq "GROUPS_SAME_NORM") {
				$protRulerNormStrg="Same normalization within groups";
			} elsif ($refInfo->{'AVG_MODE'}[0] eq "AVG_ALL") {
				$protRulerNormStrg="Average all samples";
			}
		}
	}
	else { # no design
		foreach my $stateData (@{$refInfo->{STATES}}) {
			$numStates++;
			(my $numBioRep,$stateInfo{$numStates}{'NAME'},my $repPosStrg)=split(',',$stateData);
			my $currBioRep=0;
			foreach my $replicName (split(/\./,$stateInfo{$numStates}{'NAME'})) {
				$currBioRep++;
				$replicateName{$numStates}{$currBioRep}=$replicName;
				$numTechReplicates{"$numStates:$currBioRep"}=1;
			}
			$stateInfo{$numStates}{'NAME'}=~s/\./\+/g;
			if ($numBioRep>1) {
				$stateInfo{$numStates}{'REPLIC'}=($action eq 'export')? " ($numBioRep replicates)" : "&nbsp;($numBioRep replicates)&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Correlation\" onclick=\"displayReplicCorrelation(this,$numStates,'$stateInfo{$numStates}{NAME}')\">";
				#$existReplicates=1;
			}
			else {$stateInfo{$numStates}{'REPLIC'}='';}
		}
	}

	#<Setting default displayed ratios (0 for MaxQuant no-ratio)
	if ($ratioType ne 'None') {
		my $maxDefaultRatios=($numStates <= 4)? $numStates-1 : 3; # all ratios with 1st reference
		unless (scalar keys %dispRatios) {
			foreach my $pos (1..$maxDefaultRatios) {$dispRatios{$pos}=[];}
		}
	}

	#<Peptides (retro-compatibility in case no miss-cut param)
	my ($pepRangeStrg,$pepMissCutStrg,$pepPtmStrg,$pepChargeStrg,$pepSourceStrg);
	if ($refInfo->{'PEPTIDES'}) {
		if ($quantifSoftware eq 'MaxQuant') {
			$pepRangeStrg=($refInfo->{'PEPTIDES'}[0] eq 'unique')? 'Unique peptides only' : ($refInfo->{'PEPTIDES'}[0] eq 'razor')? 'Razor+unique peptides allowed' : 'All peptides allowed';
		}
		else {
			$pepRangeStrg=($refInfo->{'PEPTIDES'}[0] eq 'unique')? 'Proteotypic peptides only' : ($refInfo->{'PEPTIDES'}[0] eq 'unique_shared')? 'Shared peptides allowed if exist proteotypic' : 'Non-proteotypic peptides allowed';
		}
		$pepMissCutStrg='Missed cleavage allowed';
		my $pepIdx=0;
		if (scalar @{$refInfo->{'PEPTIDES'}}==5) { # Missed cleavage option exists
			$pepMissCutStrg=($refInfo->{'PEPTIDES'}[1]==0)? 'Missed cleavage not allowed' : ($refInfo->{'PEPTIDES'}[1]==-1)? 'Miss-cleaved & cleaved counterparts not allowed' : 'Missed cleavage allowed';
			$pepIdx++;
		}
		$pepIdx++;
		my ($ptmAllowed,@selectedPTMs)=split(/[:,]/,$refInfo->{'PEPTIDES'}[$pepIdx]);
		if (abs($ptmAllowed)<=1) { # old options
			$pepPtmStrg=($ptmAllowed==1)? 'All modifications allowed' : ($ptmAllowed==-1)? 'Modified and unmodified matching peptides not allowed' : 'No modifications allowed';
		}
		elsif ($ptmAllowed==-3) { # Only MaxQuant so far
			$pepPtmStrg='All modifications allowed (unmodified matching peptides not allowed)';
		}
		else { # new custom selection options (2 or -2)
			my @ptmNames;
			my $sthVM=$dbh->prepare('SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DES FROM MODIFICATION WHERE ID_MODIFICATION=?');
			foreach my $modID (@selectedPTMs) {
				$sthVM->execute($modID);
				my ($psiName,$interName,$synName,$des)=$sthVM->fetchrow_array;
				my $ptmName=$psiName || $interName || $des;
				unless ($ptmName) {
					$synName=~s/^##//; $synName=~s/##$//; $synName=~s/##/, /;
					$ptmName=$synName;
				}
				push @ptmNames,$ptmName;
			}
			$sthVM->finish;
			if ($action eq 'export') {$pepPtmStrg='Modifications allowed: '.join(', ',sort{lc($a) cmp lc($b)} @ptmNames);}
			else {
				$pepPtmStrg="<A href=\"javascript:void(null)\" onmouseover=\"popup('<B>Modifications allowed:</B><BR>&nbsp;&bull;".join('<BR>&nbsp;&bull;',sort{lc($a) cmp lc($b)} @ptmNames)."')\" onmouseout=\"popout()\">Some modifications allowed<SUP>*</SUP></A>";
			}
			$pepPtmStrg.=' (Other modified and unmodified forms not allowed)' if $ptmAllowed==-2;
		}
		#if ($quantifType !~ /(PROT_RATIO|TNPQ|XIC|SIN|EMPAI)/) { #} Modified in 29/10/12
		#if ($quantifType !~ /XIC|SIN|EMPAI/) {
			$pepIdx++;
			$pepChargeStrg=($refInfo->{PEPTIDES}[$pepIdx] eq 'best')? 'Best charge state' : 'All charge states';
			$pepIdx++;
			$pepSourceStrg=($refInfo->{PEPTIDES}[$pepIdx] && $refInfo->{PEPTIDES}[$pepIdx] eq 'best')? 'Best source' : 'All sources';
		#}
	}
	my %convPtmParam=('ambiguous'=>'delocalized','exclude'=>'excluded','valid'=>'not delocalized');
	
	my $ptmQuantifStrg='';
	my ($boldTag1,$boldTag2)=($action eq 'export')? ('','') : ('<B>','</B>');
	if (scalar keys %modificationContexts) {
		my $modifContextStrg='';
		my $hasContext=0;
		foreach my $modID (keys %modificationContexts) {
			$modifContextStrg.=', ' if $modifContextStrg;
			$modifContextStrg.="$boldTag1$quantifModifInfo{NAME}{$modID}$boldTag2: ";
			if (scalar @{$modificationContexts{$modID}}) {
				if ($action ne 'export') {
					foreach my $context (@{$modificationContexts{$modID}}) {$context=~s/([a-z])/<B>$1<\/B>/g;} # emphasize targetted residue(s)
				}
				$modifContextStrg.="'".join("' '",@{$modificationContexts{$modID}})."'";
				$hasContext=1;
			}
			else {$modifContextStrg.='No context';}
		}
		if ($hasContext) {$ptmQuantifStrg="Site context(s) for ".$modifContextStrg;}
		elsif ($isProteinByPTM) {$ptmQuantifStrg="No context for enrichment PTM(s)";} # context if mandatory for Free Residue
	}
	elsif ($refInfo->{'PTM_POS'}) {
		$ptmQuantifStrg=($refInfo->{'PTM_POS'}[0]=~/(\w+):(.+)/)? "$boldTag1$modifName$boldTag2-sites positions are $boldTag1"."confirmed$boldTag2 if $ptmProbSoft{$1} probability &ge; $boldTag1$2%$boldTag2, others are $boldTag1$convPtmParam{$refInfo->{PTM_POS}[1]}$boldTag2" : "$boldTag1$modifName$boldTag2-sites are $convPtmParam{$refInfo->{PTM_POS}[0]}";
	}

	#<Bias correction
	my ($biasCorrFlag,$biasCorrLevel)=($refInfo->{'BIAS_CORRECTION'})? split(',',lc($refInfo->{'BIAS_CORRECTION'}[0])) : ('false','');
	unless ($biasCorrLevel) { $biasCorrLevel=($biasCorrFlag eq 'true')? 'ion' : ''; }
	my ($biasCorrIons,$biasCorrProteins)=($biasCorrLevel eq 'ion')? (1,0) : ($biasCorrLevel eq 'ion_prot')? (1,1) : ($biasCorrLevel eq 'prot')? (0,1) : (0,0);
	my ($biasCorrectStrg,$biasCorrectStrg2)=('','');
	if ($action eq 'export') {
		if ($refInfo->{'BIAS_CORRECTION'}) {		
			if ($refInfo->{'BIAS_CORRECTION'}[0] eq 'FALSE') {$biasCorrectStrg='None';}
			else {
				$biasCorrectStrg='Unknown'; # default
				if ($ratioType=~/S\w+Ratio/ || $quantifType eq 'PROT_ABUNDANCE') {
					if ($refInfo->{'NORMALIZATION_METHOD'}[0]) {
						$biasCorrectStrg=($normalizationNames{$refInfo->{'NORMALIZATION_METHOD'}[0]})? $normalizationNames{$refInfo->{'NORMALIZATION_METHOD'}[0]} : $refInfo->{'NORMALIZATION_METHOD'}[0];
						$biasCorrectStrg.=" at $biasCorrLevelNames{$biasCorrLevel} level";
						$biasCorrectStrg.='s' if $biasCorrLevel eq 'ion_prot';
					}
					else {$biasCorrectStrg=$normalizationNames{'global.mad.normalization'};}
				}
				else {$biasCorrectStrg="Scale normalization";} # ?
				if ($refInfo->{'BIAS_CORRECTION'}[1]) { # Prot list used for bias correction
					#if ($refInfo->{'BIAS_CORRECTION'}[1] =~ /shared/) { # proteins/sites/peptides shared by all states
					#	my $item=($refInfo->{'BIAS_CORRECTION'}[1] eq 'sharedPep')? 'peptides' : ($isModifQuantif)? 'site' : 'protein';
					#	my $numItemStrg='';
					#	if ($item ne 'peptides') {
					#		my $numItems=$refInfo->{'BIAS_CORRECTION'}[2]; # not defined for 'sharedPep'
					#		$item.='s' if ($numItems > 1);
					#		$numItemStrg=" $numItems";
					#	}
					#	$biasCorrectStrg.=" based on$numItemStrg $item shared by all States";
					#}
					$refInfo->{'BIAS_CORRECTION'}[1]=~s/#//;
					if ($refInfo->{'BIAS_CORRECTION'}[1]=~/^\d+$/) { # Custom list
						my ($listName)=$dbh->selectrow_array("SELECT CONCAT(CL.NAME,' > ',C.NAME) FROM CATEGORY C,CLASSIFICATION CL WHERE C.ID_CLASSIFICATION=CL.ID_CLASSIFICATION AND ID_CATEGORY=$refInfo->{BIAS_CORRECTION}[1]");
						my $notStrg=($refInfo->{'BIAS_CORRECTION'}[2] && $refInfo->{'BIAS_CORRECTION'}[2] eq 'not')? ' not' : '';
						$biasCorrectStrg.=" based on proteins$notStrg in List: $listName";
					}
					$biasCorrectStrg.=" using peptides shared by all States" if $refInfo->{'BIAS_CORRECTION'}[-1] eq 'sharedPep';
				}
				if ($refInfo->{'PEPTIDES_NORM'}) { # only 'manual' possible for now
					$biasCorrectStrg.=" based on set of peptides manually selected";
				}

				#$biasCorrectStrg.="</TD></TR>\n<TR><TD></TD><TD>";
				if ($refInfo->{'BIAS_COEFF'}) { # algo v1
					my $coeffIdx=0;
					foreach my $statePos (sort{$a<=>$b} keys %replicateName) {
						foreach my $replicPos (sort{$a<=>$b} keys %{$replicateName{$statePos}}) {
							foreach my $numTechRep (1..$numTechReplicates{"$statePos:$replicPos"}) {
								$biasCorrectStrg2.=$newLine if $coeffIdx;
								$biasCorrectStrg2.=$replicateName{$statePos}{$replicPos};
								$biasCorrectStrg2.=' : Tech. Rep. #'.$numTechRep if $numTechReplicates{"$statePos:$replicPos"} > 1;
								$biasCorrectStrg2.=($refInfo->{'BIAS_COEFF'}[$coeffIdx])? ' / '.$refInfo->{'BIAS_COEFF'}[$coeffIdx] : ' / ?';
								$coeffIdx++;
							}
						}
					}
				}
			}
		}
	}
	else { # HTML display
		if ($refInfo->{'BIAS_CORRECTION'}) {
			my $protNormFlag=($refInfo->{'INTRA_PROT_REF'})? 'true' : 'false';
			my $biasCorrButtonStrg='';
			if ($quantifStatus > 0 && $quantifSoftware !~ /MaxQuant|Spectronaut/) {
				my $numBoxImg=1; # check if normalization boxplots are split into multiple image files
				while (-e "$resultDir/graph/Beforeallpeptide$numBoxImg.jpeg") {
					$numBoxImg++;
				}
				$numBoxImg-- unless $numBoxImg==1;
				$biasCorrButtonStrg.="<INPUT type=\"button\" class=\"font11\" value=\"Box plots\" style=\"width:80px\" onclick=\"displayNormalizationBoxPlots(this,$biasCorrIons,$biasCorrProteins,$protNormFlag,$numBoxImg)\"/>";
			}
			#my $biasCorrButtonStrg=($quantifStatus<=0 || $quantifSoftware eq 'MaxQuant')? '' : "<INPUT type=\"button\" class=\"font11\" value=\"Box plots\" onclick=\"displayNormalizationBoxPlots(this,$biasCorrFlag,$protNormFlag)\"/>";
			if ($refInfo->{'BIAS_CORRECTION'}[0] eq 'FALSE') {$biasCorrectStrg="<TD valign=top>None</TD><TD>&nbsp;&nbsp;$biasCorrButtonStrg</TD></TR>\n";}
			else {
$biasCorrButtonStrg.="<INPUT type=\"button\" class=\"font11\" value=\"Coeff. list\" style=\"width:80px\" onclick=\"displayBiasCoefficientsList(this,$biasCorrIons,$biasCorrProteins)\"/>";
				#my $normalizationMethod=($ratioType=~/S\w+Ratio/ && $refInfo->{'NORMALIZATION_METHOD'}[0])? $normalizationNames{$refInfo->{'NORMALIZATION_METHOD'}[0]} : ($ratioType=~/S\w+Ratio/)? $normalizationNames{'global.mad.normalization'} : "Scale normalization";
				my $normalizationMethod='Unknown'; # default
				if ($ratioType=~/S\w+Ratio/ || $quantifType eq 'PROT_ABUNDANCE') {
					if ($refInfo->{'NORMALIZATION_METHOD'}[0]) {
						$normalizationMethod=($normalizationNames{$refInfo->{'NORMALIZATION_METHOD'}[0]})? $normalizationNames{$refInfo->{'NORMALIZATION_METHOD'}[0]} : $refInfo->{'NORMALIZATION_METHOD'}[0];
					}
					else {$normalizationMethod=$normalizationNames{'global.mad.normalization'};} # ?
				}
				else {$normalizationMethod="Scale normalization";}
				$biasCorrectStrg="<TD colspan=2 nowrap><B>$normalizationMethod</B>";
				if ($biasCorrLevel) {
					$biasCorrectStrg.=" at <B>$biasCorrLevelNames{$biasCorrLevel}</B> level";
					$biasCorrectStrg.='s' if $biasCorrLevel eq 'ion_prot';
				}
				if ($refInfo->{'BIAS_CORRECTION'}[1]) { # Prot list used/not used for bias correction
					#if ($refInfo->{'BIAS_CORRECTION'}[1] =~ /^shared/) { # proteins/sites/peptides shared by all states
					#	my $item=($refInfo->{'BIAS_CORRECTION'}[1] eq 'sharedPep')? 'peptides' : ($isModifQuantif)? 'site' : 'protein';
					#	my $numItemStrg='';
					#	if ($item ne 'peptides') {
					#		my $numItems=$refInfo->{'BIAS_CORRECTION'}[2]; # not defined for 'sharedPep'
					#		$item.='s' if ($numItems > 1);
					#		$numItemStrg=" $numItems";
					#	}
					#	$biasCorrectStrg.=" based on$numItemStrg <B>$item shared</B> by all States";
					#}
					$refInfo->{'BIAS_CORRECTION'}[1]=~s/#//;
					if ($refInfo->{'BIAS_CORRECTION'}[1]=~/^\d+$/) { # Custom list
						my ($listName)=$dbh->selectrow_array("SELECT CONCAT(CL.NAME,' > ',C.NAME) FROM CATEGORY C,CLASSIFICATION CL WHERE C.ID_CLASSIFICATION=CL.ID_CLASSIFICATION AND ID_CATEGORY=$refInfo->{BIAS_CORRECTION}[1]");
						my $notStrg=($refInfo->{'BIAS_CORRECTION'}[2] && $refInfo->{'BIAS_CORRECTION'}[2] eq 'not')? '<B>not</B> in' : '<B>in</B>';
						$biasCorrectStrg.=" based on proteins $notStrg List: <B>$listName</B>";
					}
					$biasCorrectStrg.=" using <B>peptides shared</B> by all States" if $refInfo->{'BIAS_CORRECTION'}[-1] eq 'sharedPep';
				}
				if ($refInfo->{'PEPTIDES_NORM'}) { # only 'manual' possible for now
					$biasCorrectStrg.=" based on <B>peptides manually selected</B> <INPUT type=\"button\" class=\"font11\" value=\" List \" onclick=\"ajaxDisplayNormPeptides(this)\">";
				}

				$biasCorrectStrg.="&nbsp;&nbsp;$biasCorrButtonStrg" if ($quantifStatus > 0 && $algoVersion < 3 && !$refInfo->{'BIAS_COEFF'}); # display only boxplot before normalization has occured

				$biasCorrectStrg.="</TD></TR>\n<TR><TD></TD><TD nowrap>";
				if ($quantifStatus > 0) {
					if ($refInfo->{'BIAS_COEFF'}) {
						# my $valueStrg='';
						# my $coeffIdx=0;
						# foreach my $statePos (sort{$a<=>$b} keys %replicateName) {
						# 	foreach my $replicPos (sort{$a<=>$b} keys %{$replicateName{$statePos}}) {
						# 		foreach my $numTechRep (1..$numTechReplicates{"$statePos:$replicPos"}) {
						# 			$valueStrg.='<BR>' if $coeffIdx;
						# 			$valueStrg.='&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&bull;'.$replicateName{$statePos}{$replicPos};
						# 			$valueStrg.=' : Tech. Rep. #'.$numTechRep if $numTechReplicates{"$statePos:$replicPos"} > 1;
						# 			$valueStrg.=($refInfo->{'BIAS_COEFF'}[$coeffIdx])? ' / '.$refInfo->{'BIAS_COEFF'}[$coeffIdx] : ' / ?';
						# 			$valueStrg.='&nbsp;';
						# 			$coeffIdx++;
						# 		}
						# 	}
						# }
						# $valueStrg="<DIV style=\"max-height:100px;overflow-y:scroll;display:inline-block\">$valueStrg</DIV>" if $coeffIdx > 4;
						# $biasCorrectStrg.="$valueStrg</TD><TD valign=\"middle\">&nbsp;&nbsp;&nbsp;$biasCorrButtonStrg</TD></TR>\n";
						$biasCorrectStrg.="$biasCorrButtonStrg</TD></TR>\n";
					}
					elsif ($algoVersion >= 3) { # fetch data from file
						# my %biasData;
						# &fetchBiasData($resultDir,\%biasData);
						# my $valueStrg='';
						# my $numAllReplics=0;
						# my $okData=0;
						# foreach my $statePos (sort{$a<=>$b} keys %replicateName) {
						# 	foreach my $replicPos (sort{$a<=>$b} keys %{$replicateName{$statePos}}) {
						# 		$valueStrg.='<BR>' if ($valueStrg && $okData);
						# 		$okData=0;
						# 		$numAllReplics++;
						# 		if ($isPepIntensity) {
						# 			$valueStrg.=$replicateName{$statePos}{$replicPos}.': 1/'.((sprintf "%.3f",$biasData{"$statePos:$replicPos"}[0])*1).' ~ '.((sprintf "%.3f",$biasData{"$statePos:$replicPos"}[1])*1);
						# 			$okData=1;
						# 		}
						# 		else { # PEP_RATIO
						# 			$okData=0;
						# 			next if $statePos==1; # State1 not listed
						# 			next unless $biasData{"$statePos:$replicPos"}; # SuperRatio: data are global to State but tagged as xxx_rep1_techRep1
						# 			$valueStrg.=($ratioType eq 'SuperRatio')?
						# 				'<B>'.$stateInfo{$statePos}{'NAME'}.' / '.$stateInfo{1}{'NAME'}.'</B>'
						# 				: $replicateName{$statePos}{$replicPos}.' <B>/</B> '.$replicateName{1}{$replicPos};
						# 			$valueStrg.=': 1/'.$biasData{"$statePos:$replicPos"}[0].' ~ '.$biasData{"$statePos:$replicPos"}[1];
						# 			$okData=1;
						# 		}
						# 	}
						# }
						# $biasCorrButtonStrg.="<BR>&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Coeff. plot\" style=\"width:80px\" onclick=\"displayBiasCoefficientsPlot(this,$biasCorrIons,$biasCorrProteins)\"/>";
						# $valueStrg="<DIV style=\"max-height:100px;overflow-y:scroll;display:inline-block\">$valueStrg</DIV>" if $numAllReplics > 5;
						# $biasCorrectStrg.="$valueStrg</TD><TD width=100% valign=\"middle\">&nbsp;&nbsp;&nbsp;$biasCorrButtonStrg</TD></TR>\n";
						$biasCorrButtonStrg.="<INPUT type=\"button\" class=\"font11\" value=\"Coeff. plot\" style=\"width:80px\" onclick=\"displayBiasCoefficientsPlot(this,$biasCorrIons,$biasCorrProteins)\"/>";
						$biasCorrectStrg.="$biasCorrButtonStrg</TD></TR>\n";
					}
				}
			}
		}
	}
	
	#<P-value correction
	my $pValueType=(!$refInfo->{'FDR_CONTROL'} || $refInfo->{'FDR_CONTROL'}[0]=~/FALSE|none/i)? 'p-value' : 'Adj. p-value';
	my ($fdrStrg,$pvalStrg,$alterStrg)=('','',''); # (very old algo only!)
	if ($ratioType eq 'Ratio') {
		($fdrStrg,$pvalStrg)=($refInfo->{'FDR_CONTROL'}[0] eq 'FALSE')? ('None','Not given') : ($refInfo->{'FDR_ALPHA'}[0].'% ('.$fdrMethods{$refInfo->{'FDR_METHOD'}[0]}.')' , $refInfo->{P_VALUE_OUTLIER}[0] );
		$alterStrg=(!$refInfo->{ALTER})? '' : ($refInfo->{ALTER}[0] eq 'two.sided')? 'Two sided' : ($refInfo->{ALTER}[0] eq 'less')? 'Lesser' : 'Greater';
	}

	####################################
	####<Exporting summary to Excel>####
	####################################
	if ($action eq 'export') {
		my $worksheet1=$workbook->add_worksheet('Summary');
		my $xlsRow=0;
		$worksheet1->set_column(0,0,50);
		$worksheet1->set_column(1,1,40);
		$worksheet1->set_column(2,2,50);
		$worksheet1->set_row($xlsRow,26);
		$worksheet1->merge_range($xlsRow,0,$xlsRow,2,decode_utf8($title),$itemFormat{'title'});

		##<Source>##
		my @quantifInfo=&promsMod::getItemInfo($dbh,'QUANTIFICATION',$selQuantifID);
		$xlsRow++;
		$worksheet1->merge_range(++$xlsRow,0,$xlsRow,2,'Source',$itemFormat{'mergeColHeader'});
		#<Name
		$worksheet1->write_string(++$xlsRow,0,'Name :',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$quantifInfo[-1]{'NAME'},$itemFormat{'mergeColText'});
		#<Location
		$worksheet1->write_string(++$xlsRow,0,'Location in myProMS :',$itemFormat{'headerR'});
		my $locationStrg='';
		foreach my $it (@quantifInfo) {
			last if $it->{'ITEM'} eq 'QUANTIFICATION';
			$locationStrg.=' > ' if $locationStrg;
			$locationStrg.=$it->{'NAME'};
		}
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$locationStrg,$itemFormat{'mergeColText'});
		#<Design
		$worksheet1->write_string(++$xlsRow,0,'Design level :',$itemFormat{'headerR'});
		my $designStrg=($designID)? 'Experiment' : 'Analysis';
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$designStrg,$itemFormat{'mergeColText'});
		if ($quantifType ne "PROT_RULER") {
			#<Peptide quantif
			$worksheet1->write_string(++$xlsRow,0,'Peptide XIC extraction engine :',$itemFormat{'headerR'});
			my ($isMCQ)=$dbh->selectrow_array("SELECT 1 FROM PARENT_QUANTIFICATION PQ,QUANTIFICATION P WHERE PQ.ID_QUANTIFICATION=$selQuantifID && PQ.ID_PARENT_QUANTIFICATION=P.ID_QUANTIFICATION AND P.FOCUS='peptide' AND P.QUANTIF_ANNOT LIKE '%::ALLCHARGESTATES=%' LIMIT 0,1");  # Keep for back compatibility ?
			my ($parQuantifAnnot)=$dbh->selectrow_array("SELECT P.QUANTIF_ANNOT FROM PARENT_QUANTIFICATION PQ,QUANTIFICATION P WHERE PQ.ID_QUANTIFICATION=$selQuantifID AND PQ.ID_PARENT_QUANTIFICATION=P.ID_QUANTIFICATION AND P.FOCUS='peptide' LIMIT 0,1");
			my ($xicSoftware) = ($parQuantifAnnot =~ /SOFTWARE=([^:]+)/);
			my ($xicSoftStrg, $xicSoftVersion) = split(';', $xicSoftware);
			my $pepQuantStrg=($isMCQ || ($xicSoftStrg && $xicSoftStrg eq 'MCQ'))? 'MassChroQ' : ($xicSoftStrg && $xicSoftStrg eq 'MQ')? 'MaxQuant' : 'Proteome Discoverer';
			$pepQuantStrg.=" (v$xicSoftVersion)" if ($xicSoftStrg && $xicSoftStrg =~ /MCQ|MQ|PD/);
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$pepQuantStrg,$itemFormat{'mergeColText'});
			#<Labeling
			$worksheet1->write_string(++$xlsRow,0,'Labeling :',$itemFormat{'headerR'});
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$labelingName{$labelType}.$swathStrg,$itemFormat{'mergeColText'});
		}
		##<Parameters>##
		$xlsRow++;
		$worksheet1->merge_range(++$xlsRow,0,$xlsRow,2,'Quantification parameters',$itemFormat{'mergeColHeader'});
		$worksheet1->write_string(++$xlsRow,0,'Quantification method :',$itemFormat{'headerR'});
		if ($quantifType eq "PROT_RULER") {
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$quantifMethDesc.' ['.$protRulerTypeStrg.', based on '.$intMetric.' values - Software: '.$softNameStrg.$softVersionStrg.']',$itemFormat{'mergeColText'});
			$worksheet1->write_string(++$xlsRow,0,'Normalization mode :',$itemFormat{'headerR'});
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$protRulerNormStrg,$itemFormat{'mergeColText'});
			$worksheet1->write_string(++$xlsRow,0,'Numerical parameters :',$itemFormat{'headerR'});
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$protRulerParamsStrg,$itemFormat{'mergeColText'});
		} else {
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$quantifMethDesc.' ['.$ratioTypeStrg.$topNString.' - Software: '.$softNameStrg.$softVersionStrg.']',$itemFormat{'mergeColText'});
		}
		#<States
		foreach my $statePos (sort{$a<=>$b} keys %stateInfo) {
			$worksheet1->write_string(++$xlsRow,0,"State #$statePos :",$itemFormat{'headerR'});
			my $replicateNames=($stateInfo{$statePos}{'NAME_LONG'})? $stateInfo{$statePos}{'NAME_LONG'} : '';
			if ($stateInfo{$statePos}{'POPUP'}) {
				$worksheet1->write_string($xlsRow,1,decode_utf8("$stateInfo{$statePos}{NAME}$replicateNames"),$itemFormat{'text'});
				$worksheet1->write_string($xlsRow,2,decode_utf8($stateInfo{$statePos}{'POPUP'}),$itemFormat{'textWrap'});
			}
			else {$worksheet1->merge_range($xlsRow,1,$xlsRow,2,decode_utf8("$stateInfo{$statePos}{NAME}$replicateNames"),$itemFormat{'mergeColText'});}
		}
		#<Peptide selection
		if ($quantifType ne "PROT_RULER") {
			my ($pepString,$numPepLines)=($quantifType !~ /SIN|EMPAI/ && $quantifSoftware !~ /SWATH|DIA/ )? ("\n•$pepChargeStrg.\n•$pepSourceStrg.",5) : ('',3);
			if ($ptmQuantifStrg) {
				$numPepLines++;
				$ptmQuantifStrg="\n•$ptmQuantifStrg.";
			}
			if (scalar keys %modificationContexts) {
				my $modifContextStrg='';
				my $hasContext=0;
				foreach my $modID (keys %modificationContexts) {
					$modifContextStrg.=', ' if $modifContextStrg;
					$modifContextStrg.="$quantifModifInfo{'EXPORT_NAME'}{$modID}: ";
					if (scalar @{$modificationContexts{$modID}}) {
						$modifContextStrg.="'".join("' '",@{$modificationContexts{$modID}})."'";
						$hasContext=1;
					}
					else {$modifContextStrg.='No context';}
				}
				if ($hasContext) {$pepString.="\n•Site context(s) for ".$modifContextStrg;}
				elsif ($isProteinByPTM) {$pepString.="\n•No context for enrichment PTM(s)";} # context if mandatory for Free Residue
			}
			if ($quantifType eq 'PROT_ABUNDANCE' && $refInfo->{'PEPTIDES'}[0] ne 'unique' && !$isModifQuantif) {
				my $MGSharedRule=(!$refInfo->{'MG_SHARED_PEP'} || $refInfo->{'MG_SHARED_PEP'}[0] eq 'best')? 'assigned to best protein' : ($refInfo->{'MG_SHARED_PEP'}[0] eq 'share')? 'used for each protein' : 'excluded from dataset'; #'exclude';
				$numPepLines++;
				$pepString.="\n•Peptides shared by multiple Match Groups are $MGSharedRule";
			}
			$worksheet1->set_row(++$xlsRow,13*$numPepLines);
			$worksheet1->write_string($xlsRow,0,'Peptide selection :',$itemFormat{'headerR'});
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,decode_utf8("•$pepRangeStrg.\n•$pepMissCutStrg.\n•$pepPtmStrg.$pepString$ptmQuantifStrg"),$itemFormat{'mergeColText'});
		}
		#<Protein selection
		my $protSelectionStrg='All visible proteins'; # default
		my $contaminantLink='without';
		if ($refInfo->{'PROTEINS'}) {
			($protSelectionStrg,$contaminantLink)=($refInfo->{'PROTEINS'}[0] eq 'exclude')? ('Exclude proteins in List ','and') : ('Restrict to proteins in List ','without');
			(my $listID=$refInfo->{'PROTEINS'}[1])=~s/#//;
			my ($listStrg)=$dbh->selectrow_array("SELECT CONCAT(T.NAME,' > ',L.NAME) FROM CATEGORY L,CLASSIFICATION T WHERE L.ID_CLASSIFICATION=T.ID_CLASSIFICATION AND ID_CATEGORY=$listID");
			$protSelectionStrg.="$listStrg";
		}
		if ($refInfo->{'EXCL_CONTAMINANTS'} && $refInfo->{'EXCL_CONTAMINANTS'}[0]) {
			$protSelectionStrg.=" $contaminantLink contaminants";
		}
		$worksheet1->write_string(++$xlsRow,0,'Protein selection :',$itemFormat{'headerR'});
		#$protSelectionStrg.=' ('.$protQuantifiedStrg.')' if $protQuantifiedStrg; <- contains ajax call to num. quantified prot
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$protSelectionStrg,$itemFormat{'mergeColText'});
		if ($quantifType ne "PROT_RULER") {
			#<Bias
			$worksheet1->write_string(++$xlsRow,0,'Bias correction :',$itemFormat{'headerR'});
			if ($biasCorrectStrg2) {
				$worksheet1->write_string($xlsRow,1,$biasCorrectStrg,$itemFormat{'text'});
				$worksheet1->write_string($xlsRow,2,$biasCorrectStrg2,$itemFormat{'textWrap'});
			}
			else {$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$biasCorrectStrg,$itemFormat{'mergeColText'});}
			if ($ratioType ne "None") {
				#<Inf ratios
				$worksheet1->write_string(++$xlsRow,0,'Infinite ratios :',$itemFormat{'headerR'});
				my $infRatioString=($refInfo->{'MIN_INF_RATIOS'} && $refInfo->{'MIN_INF_RATIOS'}[0])? 'Avoided whenever possible' : 'Not avoided';
				$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$infRatioString,$itemFormat{'mergeColText'});
			}
			#<Advanced settings
			if ($ratioType eq 'Ratio') {
				if ($refInfo->{'THRESHOLD_CV'}) { # replicates used TODO: Make sure THRESHOLD_CV is reported for design quantif
					$worksheet1->write_string(++$xlsRow,0,'Variation coefficient threshold between replicates :',$itemFormat{'headerR'});
					my $biasCorrStrg.=($refInfo->{'THRESHOLD_CV'}[0] eq 'FALSE')? 'Auto' : $refInfo->{'THRESHOLD_CV'}[1];
					$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$biasCorrStrg,$itemFormat{'mergeColText'});
				}
				$worksheet1->write_string(++$xlsRow,0,'FDR control :',$itemFormat{'headerR'});
				$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$fdrStrg,$itemFormat{'mergeColText'});
				$worksheet1->write_string(++$xlsRow,0,"$pValueType threshold for outlier detection :",$itemFormat{'headerR'});
				$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$pvalStrg,$itemFormat{'mergeColText'});
				if ($refInfo->{'ALTER'}) {
					$worksheet1->write_string(++$xlsRow,0,'Alternative hypothesis for comparison :',$itemFormat{'headerR'});
					$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$alterStrg,$itemFormat{'mergeColText'});
				}
				if ($refInfo->{'CONFIDENCE_LEVEL'}) {
					$worksheet1->write_string(++$xlsRow,0,'Confidence interval on protein abundance :',$itemFormat{'headerR'});
					$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$refInfo->{'CONFIDENCE_LEVEL'}[0],$itemFormat{'mergeColText'});
				}
			}
		}

		##<Creation date & user
		$updateDate=&promsMod::formatDate($updateDate);
		$worksheet1->write_string(++$xlsRow,0,'Created :',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,decode_utf8("$updateDate by $updateUser"),$itemFormat{'mergeColText'});

		##<Export settings>##
		$xlsRow++;
		$worksheet1->merge_range(++$xlsRow,0,$xlsRow,2,'Export settings',$itemFormat{'mergeColHeader'});
		if ($ratioType eq 'None') {  # Abundance quantifications
			$worksheet1->set_row(++$xlsRow, 13 * scalar keys %dispStates);
			$worksheet1->write_string($xlsRow, 0, 'States exported :', $itemFormat{'headerR'});
			my $disStatesStrg='';
			foreach my $statePos (sort {abs($a) <=> abs($b)} keys %dispStates) {
				$disStatesStrg .= "\n" if $disStatesStrg;
				my ($repNb, $bioReps, $state);
				if ($quantifType eq 'PROT_RULER') {
					$state = $statePos;
				} else {
					($repNb, $bioReps, $state) = split(',', $refInfo->{STATES}[$statePos - 1]);
					if ($designID) {
						$state = $condToState{$state};
					}
				}
				$disStatesStrg .= "•$stateInfo{$state}{NAME}";
			}
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,decode_utf8($disStatesStrg),$itemFormat{'mergeColText'});

		}
		else {  # Ratio quantifications
			$worksheet1->set_row(++$xlsRow,13 * scalar keys %dispRatios);
			$worksheet1->write_string($xlsRow,0,'Ratios exported :',$itemFormat{'headerR'});
			my $disRatiosStrg='';
			foreach my $rPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
				my $rPos0=abs($rPos);
				$disRatiosStrg.="\n" if $disRatiosStrg;
				my ($testStatePos,$refStatePos)=split(/\//,$refInfo->{RATIOS}[$rPos0-1]);
				my $ratioTag='';
				if ($designID) {
					if ($testStatePos=~/%/) { # SuperRatio
						$ratioTag='°';
						$testStatePos=~s/%\d+//;
						$refStatePos=~s/%\d+//;
					}
					$testStatePos=$condToState{$testStatePos};
					$refStatePos=$condToState{$refStatePos};
				}
				$disRatiosStrg.=($rPos>0)? "•$stateInfo{$testStatePos}{NAME}$ratioTag/$stateInfo{$refStatePos}{NAME}$ratioTag" : "•$stateInfo{$refStatePos}{NAME}$ratioTag/$stateInfo{$testStatePos}{NAME}$ratioTag";
			}
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,decode_utf8($disRatiosStrg),$itemFormat{'mergeColText'});
			my $dispFoldChangeStrg=($filterProtIDs)? '*User selection*' : $dispFoldChange;
			my ($foldChgTypeStrg,$dispFcStrg)=($foldChangeType eq 'abs')? ('Absolute fold change ≥',"$dispFoldChangeStrg (+INF=1000, 1/INF=0.001)") : ($foldChangeType eq 'up')? ('Fold change ≥',"$dispFoldChange (INF=1000)") : ('Fold change ≤',"1/$dispFoldChange (1/INF=0.001)");
			$worksheet1->write_string(++$xlsRow,0,decode_utf8($foldChgTypeStrg),$itemFormat{'headerR'});
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$dispFcStrg,$itemFormat{'mergeColText'});
			$worksheet1->write_string(++$xlsRow,0,decode_utf8("$pValueType ≤"),$itemFormat{'headerR'});
			my $pValueStrg=($filterProtIDs)? '*User selection*' : ($allowInfRatio)? $dispPvalue.' (+/-INF allowed)' : $dispPvalue;
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$pValueStrg,$itemFormat{'mergeColText'});
			if ($quantifType eq 'PROT_RATIO_PEP') {
				#$worksheet1->write_string(++$xlsRow,0,'Standard deviation ≤',$itemFormat{'headerR'});
				#$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$dispStdDev,$itemFormat{'mergeColText'});
				$worksheet1->write_string(++$xlsRow,0,decode_utf8('Coefficient of variation ≤'),$itemFormat{'headerR'});
				my $cvStrg=($filterProtIDs)? '*User selection*' : ($dispCV)? $dispCV.'%' : '*No filter*';
				$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$cvStrg,$itemFormat{'mergeColText'});
			}
		}
		unless ($quantifType eq 'PROT_RULER') {
			#my $numPepStrg=($dispPepType eq 'distinct')? 'Distinct peptides used ≥' : ($dispPepType eq 'distinctRep')? 'Distinct peptides used per replicate ≥' : ($dispPepType eq 'allRep')? 'Peptides used per replicate ≥' : 'Peptides used ≥';
			# my $numPepStrg=($dispPepType eq 'distinct')? 'Distinct p' : ($dispPepType =~ /^msms/ && $numPepScope eq 'ratio')? 'Identified p' : ($dispPepType =~ /^msms/)? 'Min. ident. p' : 'P';
			my ($pcRatioStrg,$pcStateStrg)=($dispPepType eq 'msmsPc')? ('% identified','Min. %') : ('Identified','Min.');
			my $numPepStrg=($dispPepType eq 'distinct')? 'Distinct p' : ($dispPepType =~ /^msms/ && $numPepScope eq 'ratio')? "$pcRatioStrg p" : ($dispPepType =~ /^msms/)? "$pcStateStrg ident. p" : 'P';
			$numPepStrg.='eptides used';
			if ($ratioType eq 'None') {
				$numPepStrg.=' in at least one state';
			}
			elsif ($numPepScope =~ /^replicate/ || $dispPepType =~ /^msms/) {
				$numPepStrg.=' per replicate' if $numPepScope =~ /^replicate/;
				$numPepStrg.=($numPepScope =~ /One$/)? ' in one state' : ($numPepScope =~ /Both$/)? ' in both states' : ($numPepScope =~ /Test$/)? ' in Test state' : ' in Reference state';
			}
			else {$numPepStrg.=' for ratio';}
			$numPepStrg.=' ≥';
			$worksheet1->write_string(++$xlsRow,0,decode_utf8($numPepStrg),$itemFormat{'headerR'});
			my $numPepValueStrg=$dispNumPep;
			if ($numPepScope =~ /^replicate/) {
				$numPepValueStrg.=" in at least $dispNumReplic replicate";
				$numPepValueStrg.='s' if $dispNumReplic > 1;
			}
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$numPepValueStrg,$itemFormat{'mergeColText'});
		}
		$worksheet1->write_string(++$xlsRow,0,'Sort by :',$itemFormat{'headerR'});
		my $sortOption;
		if ($dispSort=~/ratio_(\d+)/) {
			my $rPos=$1;
			my ($testStatePos,$refStatePos)=split(/\//,$refInfo->{RATIOS}[$rPos-1]);
			my $ratioTag='';
			if ($designID) {
				if ($testStatePos=~/%/) { # SuperRatio
					$ratioTag=$encodedDegStrg;
					$testStatePos=~s/%\d+//;
					$refStatePos=~s/%\d+//;
				}
				$testStatePos=$condToState{$testStatePos};
				$refStatePos=$condToState{$refStatePos};
			}
			$sortOption="Ratio: '$stateInfo{$testStatePos}{NAME}$ratioTag/$stateInfo{$refStatePos}{NAME}$ratioTag'";
		}
		if ($dispSort=~/state_(\d+)/) {
			$sortOption="$quantifParamInfo{$dispMeasure}[1]: $stateInfo{$1}{NAME}";
		}
		else {$sortOption=($dispSort eq 'identifier')? 'Identifiers' : ($dispSort eq 'peptide')? 'Peptides used' : ($dispSort eq 'mw')? 'Molecular weight' : '';}
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$sortOption,$itemFormat{'mergeColText'});
		if ($restrictListID) {
			my ($listStrg)=$dbh->selectrow_array("SELECT CONCAT(CL.NAME,' > ',C.NAME) FROM CATEGORY C,CLASSIFICATION CL WHERE C.ID_CLASSIFICATION=CL.ID_CLASSIFICATION AND ID_CATEGORY=$restrictListID");
			my $actionStrg=($listAction eq 'restrict')? 'Restrict to' : 'Exclude';
			$worksheet1->write_string(++$xlsRow,0,"$actionStrg proteins in List :",$itemFormat{'headerR'});
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$listStrg,$itemFormat{'mergeColText'});
		}

		return %stateInfo;
	} # end of &printQuantificationSummary when action='export'


	####################################
	####<Displaying summary in HTML>####
	####################################
	my $correlImgFileName=($quantifSoftware eq 'myProMS')? 'Beforematrixcorrelation.jpeg' : 'AllProtein_Before_matrixcorrelation.jpeg'; # !!! Temp! Remove "AllProtein_" for SWATH!!!!
	#my $correlMatFileName=($quantifSoftware eq 'myProMS')? 'Beforematrixcorrelation.txt' : 'AllProtein_Before_matrixcorrelation.txt'; # !!! Temp! Remove "AllProtein_" for SWATH!!!!
	my $correlMatFileName=($algoVersion >= 2)? 'Beforematrixcorrelation.txt' : 'AllProtein_Before_matrixcorrelation.txt'; # !!! Temp! Remove "AllProtein_" for SWATH!!!!
	my $ajaxPepAction=($quantifSoftware=~/SWATH|DIA|Spectronaut/)? 'ajaxPepSwath' : ($quantifSoftware eq 'MaxQuant')? 'ajaxPepMaxQuant' : ($protRulerCase)? 'ajaxPepProtRuler' : ($ratioType eq 'None')? 'ajaxPepAbund' : 'ajaxPepRatio';
	#my $dispMeasParamStrg=(($quantifSoftware =~ /MaxQuant|Spectronaut/ && $ratioType eq 'None') || $protRulerCase)? "&dispMeasure=$dispMeasure" : '';
	my $dispMeasParamStrg=($ratioType eq 'None' || $protRulerCase)? "&dispMeasure=$dispMeasure" : '';
	print qq
|<SCRIPT type="text/javascript">
// AJAX ------>
function getXMLHTTP() {
	var xhr=null;
	if(window.XMLHttpRequest) {// Firefox & others
		xhr = new XMLHttpRequest();
	}
	else if(window.ActiveXObject){ // Internet Explorer
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
// <--- AJAX
function getDomOffset( Obj, Prop ) {
	var iVal = 0;
	while (Obj && Obj.tagName != 'BODY') {
		eval('iVal += Obj.' + Prop + ';');
		Obj = Obj.offsetParent;
	}
	return iVal;
}
function updateSettings(act) {
	if (act=='more') {
		document.getElementById('moreSettings').style.display='none';
		document.getElementById('lessSettings').style.display='block';
		document.getElementById('advancedSetDIV').style.display='block';
	}
	else {
		document.getElementById('moreSettings').style.display='block';
		document.getElementById('lessSettings').style.display='none';
		document.getElementById('advancedSetDIV').style.display='none';
	}
}
function getElementPosition(e) {
	var left=0;
	var top=0;
	while (e.offsetParent != undefined && e.offsetParent != null) { //Adding parent item position//
		left += e.offsetLeft + (e.clientLeft != null ? e.clientLeft : 0);
		top += e.offsetTop + (e.clientTop != null ? e.clientTop : 0);
		e = e.offsetParent;
	}
	return [top,left];
}
function checkQuantifStatus() { // for on-going quantif only
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxCheckStatus&id_quantif=$selQuantifID",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			var newStatus=parseInt(XHR.responseText);
			if (newStatus != $quantifStatus) { // quantif status has changed
				window.location.reload();
			}
			else {
				setTimeout(checkQuantifStatus,10000);
			}
		}
	}
	XHR.send(null);
}
function manageTemplate(value){
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxManageTemplate&id_quantif=$selQuantifID&projectID=$projectID&TEMPLATE="+value,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			/*document.getElementById('templateSpan').innerHTML='';*/
			document.getElementById('templateSpan').innerHTML=XHR.responseText;
			if (value == 'add') {
				alert("Succesfully saved as template.");
			}
			else {
				alert("Succesfully removed from templates.");
			}
		}
	}
	XHR.send(null);
}

/* Color list for plots */
var hColors=['#E18B6B','#95B9C7','#7E2217','#9A9A9A','#8AFB17','#FBB917','#F660AB','#000000','#4AA02C']; // reverse of volcanoPlot library
|;
	if ($quantifSoftware ne 'MaxQuant') {
		my $biasCorrProtOptStrg='';
		if ($biasCorrLevel =~ /prot/) { # ABUNDANCE
			foreach my $refMeas (@{$proteinQuantifFamilies{'MEASURES'}{$quantifType}}) {
				my ($measCode,$measName,$isOptional)=@{$refMeas};
				next unless $abundanceMeasures{$measCode};
				next unless $measCode=~/_INT|LFQ/;
				$biasCorrProtOptStrg="<OPTGROUP label='$biasCorrLevelNames{'uf_prot'}:'>" if (!$biasCorrProtOptStrg && $biasCorrLevel eq 'ion_prot'); # multi-level
				$biasCorrProtOptStrg.="<OPTION value='$measCode'>$measName</OPTION>";
			}
			$biasCorrProtOptStrg.='</OPTGROUP>' if ($biasCorrProtOptStrg && $biasCorrLevel eq 'ion_prot'); # multi-level
		}
		print qq
|function displayNormalizationBoxPlots(biasButton,biasCorrIons,biasCorrProteins,isProtNorm,numBoxImg) { // isProtNorm -> sites normalized by prot
	var hasCorrection=biasCorrIons \|\| biasCorrProteins;
	var imageStrg="<TABLE cellpadding=0 cellspacing=0><TR>";
	var hSpan=(hasCorrection && isProtNorm)? 3 : (hasCorrection \|\| isProtNorm)? 2 : 1;
	if ('$quantifSoftware'.match('MSstats')) {imageStrg+="<TH class=\\"title2\\">Fragment intensities";}
	else if (hasCorrection \|\| isProtNorm) {
		var corrLevelStrg=(!hasCorrection)? '' : (biasCorrIons && !biasCorrProteins)? '$biasCorrLevelNames{uf_ion}' : "<SELECT id='boxPlotCorrMeasSEL' onchange='updateCorrBoxPlotImages("+biasCorrIons+","+biasCorrProteins+","+isProtNorm+",this,document.getElementById(\\"boxPlotImgNumSEL\\"))' class='title3'>";
		if (biasCorrIons && biasCorrProteins) {corrLevelStrg+="<OPTION value='ion'>$biasCorrLevelNames{uf_ion}</OPTION>";}
		if (biasCorrProteins) {corrLevelStrg+="$biasCorrProtOptStrg</SELECT>";}
		imageStrg+="<TH colspan="+hSpan+" class=\\"title2\\">"+corrLevelStrg+" normalization";
		if (isProtNorm) {
			imageStrg+="&nbsp;for&nbsp;<SELECT class=\\"title3\\" onchange=\\"updateNormalizationDisplay(this.value)\\"><OPTION value=\\"sites\\">Sites</OPTION><OPTION value=\\"protRef\\">Proteins</OPTION></SELECT>";
		}
	}
	else {imageStrg+="<TH colspan="+hSpan+" class=\\"title2\\">Peptide intensities";}
	imageStrg+="&nbsp;&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></TH></TR>";
	
	if (numBoxImg > 1) {
		imageStrg+="<TR><TH colspan="+hSpan+" align='left'>Sample set:<SELECT id='boxPlotImgNumSEL' onchange='updateCorrBoxPlotImages("+biasCorrIons+","+biasCorrProteins+","+isProtNorm+",document.getElementById(\\"boxPlotCorrMeasSEL\\"),this)' class='title3'>";
		for (let i=1; i<=numBoxImg; i++) {imageStrg+="<OPTION value='"+i+"'>"+i+"</OPTION>";}
		imageStrg+="</SELECT> of "+numBoxImg+"</TH></TR>";
	}
	if (hasCorrection \|\| isProtNorm) {
		imageStrg+="<TR><TH>Raw data</TH>"; // Before
		if (hasCorrection) {imageStrg+="<TH>Bias correction</TH>";} // After
		if (isProtNorm) {imageStrg+="<TH id='pepSitesTH'>Protein-level correction</TH>";}
		imageStrg+="</TR>";
	}

	/*var pepBoxRepRaw,pepBoxRepBias;*/
	var [pepBoxRepRawStrg,pepBoxRepBiasStrg,pepBoxRepFinalStrg,pepBoxRefBeforeStrg,pepBoxRefAfterStrg]=['','','','',''];
	if ('$ratioType'==='Ratio') {
		/*pepBoxRepRaw='pep_box_repl_before.png';*/
		pepBoxRepRawStrg="<IMG src='$graphDirHtml/pep_box_repl_before.png' width=480 heigth=480>";
		/*pepBoxRepBias='pep_box_repl_after.png';*/
		pepBoxRepBiasStrg="<IMG src='$graphDirHtml/pep_box_repl_after.png' width=480 heigth=480>";
	}
	else { //(Super/Simple)Ratio,ABUNDANCE,MSstats [Files: Beforeallpeptide.jpeg (raw) > Afternormallpeptide.jpeg (+/-bias) > Afterallpeptide.jpeg (+/-protRef,outliers)]//
		var pepBoxRepRaw='Beforeallpeptide';
		var pepBoxRepBias=(isProtNorm)? 'Afternormallpeptide' : 'Afterallpeptide';
		pepBoxRepRawStrg+="<IMG id='boxplotImg_A' src='' width=480 heigth=480>";
		if (hasCorrection) {pepBoxRepBiasStrg+="<IMG id='boxplotImg_B' src='' width=480 heigth=480>";}
		if (isProtNorm) {
			pepBoxRepFinalStrg+="<IMG id='boxplotImg_C' src='' width=480 heigth=480>";
			pepBoxRefBeforeStrg="<IMG id='boxplotImg_D' src='' width=480 heigth=480>";
			if (hasCorrection) {pepBoxRefAfterStrg="<IMG id='boxplotImg_E' src='' width=480 heigth=480>";}
		}
	}
	/* Current dataset */
	imageStrg+="<TR id='pepSitesTR'><TD>"+pepBoxRepRawStrg+"</TD>";
	if (hasCorrection) {imageStrg+="<TD>"+pepBoxRepBiasStrg+"</TD>";}
	if (isProtNorm) {imageStrg+="<TD>"+pepBoxRepFinalStrg+"</TD>";}
	imageStrg+="</TR>";
	/* Reference dataset */
	if (isProtNorm) {
		imageStrg+="<TR id='pepProtRefTR' style='display:none'><TD>"+pepBoxRefBeforeStrg+"</TD>";
		if (hasCorrection) {imageStrg+="<TD>"+pepBoxRefAfterStrg+"</TD>";}
		else {imageStrg+="<TD></TD>";}
		imageStrg+="<TD></TD></TR>";
	}
	imageStrg+="</TABLE>";
	
	document.getElementById('infoDIV').innerHTML=imageStrg;
	if ('$ratioType'!=='Ratio') {updateCorrBoxPlotImages(biasCorrIons,biasCorrProteins,isProtNorm,document.getElementById('boxPlotCorrMeasSEL'),document.getElementById('boxPlotImgNumSEL'));}
	var [top,left]=getElementPosition(biasButton);
	top+=30;
	left-=(hasCorrection)? 460 : 230;
// alert('window.innerWidth='+window.innerWidth+', document.body.clientWidth='+document.body.clientWidth+' ,document.documentElement.clientWidth='+document.documentElement.clientWidth);
	// var left=(document.body.clientWidth)? (document.body.clientWidth/2)-500 : (window.innerWidth)? (window.innerWidth/2)-500 : (document.documentElement.clientWidth)? (document.documentElement.clientWidth/2)-500 : ;
	if (left < 0) {left=0;}
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = left+'px';
	displayDiv.style.top = top+'px';
	displayDiv.style.display='block';
}
function updateCorrBoxPlotImages(biasCorrIons,biasCorrProteins,isProtNorm,corrMeasSel,setNumSel) {
	var hasCorrection=biasCorrIons \|\| biasCorrProteins;
	var imgPath;
	if (corrMeasSel \|\| (!biasCorrIons && biasCorrProteins)) { // protein/site-level norm
		imgPath=(corrMeasSel.value==='ion')? '$graphDirHtml' : "$protNormDirHtml/"+corrMeasSel.value+'/results/graph';
	}
	else {imgPath='$graphDirHtml';} // peptide-level norm
	const numberTag=(setNumSel)? setNumSel.value : '';
	var pepBoxRepRaw='Beforeallpeptide';
	var pepBoxRepBias=(isProtNorm)? 'Afternormallpeptide' : 'Afterallpeptide';

	document.getElementById('boxplotImg_A').src=imgPath+'/'+pepBoxRepRaw+numberTag+'.jpeg';
	if (hasCorrection) {document.getElementById('boxplotImg_B').src=imgPath+'/'+pepBoxRepBias+numberTag+'.jpeg';}
	if (isProtNorm) {
		document.getElementById('boxplotImg_C').src=imgPath+'/Afterallpeptide'+numberTag+'.jpeg';
		document.getElementById('boxplotImg_D').src=imgPath+'/BeforeREFallpeptide'+numberTag+'.jpeg';
		if (hasCorrection) {document.getElementById('boxplotImg_E').src=imgPath+'/AfternormREFallpeptide'+numberTag+'.jpeg';}
	}
}

function displayBiasCoefficientsList(coeffButton,biasCorrIons,biasCorrProteins) {
	var [top,left]=getElementPosition(coeffButton);
	top+=20;
	left-=250;
	if (left < 0) {left=0;}
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = left+'px';
	displayDiv.style.top = top+'px';
	displayDiv.style.display='block';

	var [coeffStrg,startLevel]=(biasCorrIons && !biasCorrProteins)? ['<CENTER><FONT class="title3">$biasCorrLevelNames{uf_ion}','ion'] : ["<CENTER><SELECT id='biasCorrListSEL' onchange='updateBiasCoefficientsList(this.value)' class='title3'>",null];
	if (biasCorrIons && biasCorrProteins) {coeffStrg+="<OPTION value='ion'>$biasCorrLevelNames{uf_ion}</OPTION>";}
	if (biasCorrProteins) {coeffStrg+="$biasCorrProtOptStrg</SELECT><FONT class=\\"title3\\">";}
	coeffStrg+="-level bias correction coefficients</FONT>&nbsp;&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></CENTER>";
	coeffStrg+="<DIV id='coeffListDIV' style='max-height:400; overflow:auto; padding:10px 30px 10px 30px'></DIV>";
	document.getElementById('infoDIV').innerHTML=coeffStrg;
	// initialize default view:
	updateBiasCoefficientsList(startLevel \|\| document.getElementById('biasCorrListSEL').value);
}
function updateBiasCoefficientsList(biasLevel) {
	var coeffListDiv=document.getElementById('coeffListDIV');
	coeffListDiv.innerHTML='<BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;<IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\">&nbsp;&nbsp;&nbsp;&nbsp;<BR><BR>'; // Clean any previous view
	
	/*Creation of the XMLHTTPRequest object*/
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxDispBiasCoeffList&id_quantif=$selQuantifID&bias_level="+biasLevel+"&projectID=$projectID",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			coeffListDiv.innerHTML=XHR.responseText;
		}
	}
	XHR.send(null);	
}

function displayBiasCoefficientsPlot(coeffButton,biasCorrIons,biasCorrProteins) {
	var [top,left]=getElementPosition(coeffButton);
	top+=20;
	left-=250;
	if (left < 0) {left=0;}
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = left+'px';
	displayDiv.style.top = top+'px';
	displayDiv.style.display='block';
	
	var [coeffStrg,startLevel]=(biasCorrIons && !biasCorrProteins)? ['<CENTER><FONT class="title3">$biasCorrLevelNames{uf_ion}','ion'] : ["<CENTER><SELECT id='biasCorrPlotSEL' onchange='updateBiasCoefficientsPlot(this.value)' class='title3'>",null];
	if (biasCorrIons && biasCorrProteins) {coeffStrg+="<OPTION value='ion'>$biasCorrLevelNames{uf_ion}</OPTION>";}
	if (biasCorrProteins) {coeffStrg+="$biasCorrProtOptStrg</SELECT><FONT class=\\"title3\\">";}
	coeffStrg+="-level bias correction coefficients</FONT>&nbsp;&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></CENTER>";
	coeffStrg+="<DIV id='coeffPlotDIV'></DIV>";
	document.getElementById('infoDIV').innerHTML=coeffStrg;
	// initialize default view:
	updateBiasCoefficientsPlot(startLevel \|\| document.getElementById('biasCorrPlotSEL').value);
}
function updateBiasCoefficientsPlot(biasLevel) {
	document.getElementById('coeffPlotDIV').innerHTML='<BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;<IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\">&nbsp;&nbsp;&nbsp;&nbsp;<BR><BR>'; // Clean any previous view
	
	/*Creation of the XMLHTTPRequest object*/
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxDispBiasCoeffPlot&id_quantif=$selQuantifID&bias_level="+biasLevel+"&projectID=$projectID",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			eval(XHR.responseText); // pure JS
		}
	}
	XHR.send(null);	
}

function updateNormalizationDisplay(type) {
	var [dispSites,dispProtRef]=(type=='sites')? ['','none'] : ['none',''];
	document.getElementById('pepSitesTR').style.display=dispSites;
	document.getElementById('pepSitesTH').style.display=dispSites;
	document.getElementById('pepProtRefTR').style.display=dispProtRef;
}
function ajaxDisplayNormPeptides(pepButton) {
	var [top,left]=getElementPosition(pepButton);
	top+=20;
	left-=600;
	if (left < 0) {left=0;}
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = left+'px';
	displayDiv.style.top = top+'px';
	displayDiv.style.display='block';
	var infoDiv=document.getElementById('infoDIV');
	infoDiv.innerHTML='<BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;<IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\">&nbsp;&nbsp;&nbsp;&nbsp;<BR><BR>';
	// Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxDispNormPep&id_quantif=$selQuantifID&projectID=$projectID",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			infoDiv.innerHTML=(XHR.responseText);
		}
	}
	XHR.send(null);
}
function ajaxComputeQuantifiedProteins(paramCode) {
	var protSpan=document.getElementById('protQuantifiedSPAN');
	protSpan.innerHTML='<B>Computing...</B>';
	/*Creation of the XMLHTTPRequest object*/
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxCompQuantifProt&id_quantif=$selQuantifID&isModif=$isModifQuantif&qMethID=$quantifMethodID&paramCode="+paramCode,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			var [numProt,numSites]=XHR.responseText.split(':');
			var protQuantifiedStrg='<B>'+numProt+'</B> protein';
			if (numProt*1 > 1) protQuantifiedStrg+='s';
			if (numSites*1) {
				protQuantifiedStrg+=', <B>'+numSites+"</B> $modifName$focusKeyWord";
				if (numSites > 1 && !protQuantifiedStrg.match(/s\$/)) protQuantifiedStrg+='s'; //WARNING: $modifName contains single quote!
			}
			protQuantifiedStrg+=' quantified';
			protSpan.innerHTML=protQuantifiedStrg;
		}
	}
	XHR.send(null);
}
function displayMissingValues(missButton) {
	var [top,left]=getElementPosition(missButton);
	top+=20;
	left-=150;
	if (left < 0) {left=0;}
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = left+'px';
	displayDiv.style.top = top+'px';
	displayDiv.style.display='block';
	var infoDiv=document.getElementById('infoDIV');
	infoDiv.innerHTML='<BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;<IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\">&nbsp;&nbsp;&nbsp;&nbsp;<BR><BR>';
	// Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxDispMissValues&id_quantif=$selQuantifID&projectID=$projectID",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			infoDiv.innerHTML=(XHR.responseText);
		}
	}
	XHR.send(null);
}
function displayReplicCorrelation(corrButton,statePos,replicName) {
	var [top,left]=getElementPosition(corrButton);
	top+=20;
	left-=250;
	if (left < 0) {left=0;}
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = left+'px';
	displayDiv.style.top = top+'px';
	displayDiv.style.display='block';
	var infoDiv=document.getElementById('infoDIV');
	var htmlStrg;
	if (statePos) {
		htmlStrg="<TABLE cellpadding=0 cellspacing=0><TR><TH colspan=2 class=\\"title2\\">Correlation between replicates for '"+replicName+"'&nbsp;&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></TH></TR>";
		htmlStrg+='<TR><TD><IMG src="$graphDirHtml/pep_corr_repl_cond'+statePos+'.png"></TD></TR></TABLE>';
		infoDiv.innerHTML=htmlStrg;
	}
	else { // S(uper/imple)Ratio
		htmlStrg="<TABLE cellpadding=0 cellspacing=0><TR><TH class=\\"title2\\" colspan=2>Correlation between all replicates&nbsp;&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></TH></TR>";
		if ($biasCorrFlag) {
			htmlStrg+='<TR><TD></TD><TH class="darkBg">Before normalization</TH></TR>';
			htmlStrg+='<TR><TH class="darkBg">A<BR>f<BR>t<BR>e<BR>r<BR><BR>n<BR>o<BR>r<BR>m<BR>a<BR>l<BR>i<BR>z<BR>a<BR>t<BR>i<BR>o<BR>n</TH>';
		}
		else {
			htmlStrg+='<TR><TD></TD><TH class="darkBg">No normalization</TH></TR>';
			htmlStrg+='<TR><TD></TD>';
		}
		htmlStrg+='<TD><DIV id="globCorrDIV"><BR><BR><IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR></DIV></TD></TR></TABLE>';
		infoDiv.innerHTML=htmlStrg;
		//Creation of the XMLHTTPRequest object
		var XHR = getXMLHTTP();
		if (!XHR) {
			return false;
		}
		XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxDispCorrel&id_quantif=$selQuantifID&projectID=$projectID",true);
		XHR.onreadystatechange=function() {
			if (XHR.readyState==4 && XHR.responseText) {
				document.getElementById('globCorrDIV').innerHTML='';
				eval(XHR.responseText); // javascript code
			}
		}
		XHR.send(null);
	}
}
function displayProteinQuantities(qtyButton) {
	var [top,left]=getElementPosition(qtyButton);
	top+=20;
	left-=350;
	if (left < 0) {left=0;}
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = left+'px';
	displayDiv.style.top = top+'px';
	displayDiv.style.display='block';
	var infoDiv=document.getElementById('infoDIV');
	infoDiv.innerHTML='<BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;<IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\">&nbsp;&nbsp;&nbsp;&nbsp;<BR><BR>';
	// Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxDispProtQty&id_quantif=$selQuantifID&projectID=$projectID",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			infoDiv.innerHTML=(XHR.responseText);
		}
	}
	XHR.send(null);
}
function displayCustomProts(custProtButton) {
	var [top,left]=getElementPosition(custProtButton);
	top+=20;
	left-=250;
	if (left < 0) {left=0;}
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = left+'px';
	displayDiv.style.top = top+'px';
	displayDiv.style.display='block';
	var infoDiv=document.getElementById('infoDIV');
	var htmlStrg;
	htmlStrg="<TABLE cellpadding=0 cellspacing=0><TR><TH class=\\"title2\\" colspan=2>Custom proteins&nbsp;&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></TH></TR>";
	htmlStrg+='<TD><DIV id="customProtDIV"><BR><BR></DIV></TD></TR></TABLE>';
	infoDiv.innerHTML=htmlStrg;
	// Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxDispCustProt&id_quantif=$selQuantifID&projectID=$projectID",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			document.getElementById('customProtDIV').innerHTML=(XHR.responseText);
		}
	}
	XHR.send(null);
}
/*
function displayReplicCorrelation_old(corrButton,statePos,replicName) {
	var imageStrg;
	if (statePos) {
		imageStrg="<TABLE cellpadding=0 cellspacing=0><TR><TH colspan=2 class=\\"title2\\">Correlation between replicates for '"+replicName+"'&nbsp;&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></TH></TR>";
		imageStrg+='<TR><TD><IMG src="$graphDirHtml/pep_corr_repl_cond'+statePos+'.png"></TD></TR></TABLE>';
	}
	else { // S(uper/imple)Ratio
		imageStrg="<TABLE cellpadding=0 cellspacing=0><TR><TH class=\\"title2\\">Correlation between all replicates&nbsp;&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></TH></TR>";
		imageStrg+='<TR><TD><IMG src="$graphDirHtml/$correlImgFileName" width=500 heigth=500></TD></TR></TABLE>';
	}
	var infoDiv=document.getElementById('infoDIV').innerHTML=imageStrg;
	var [top,left]=getElementPosition(corrButton);
	top+=20;
	left-=250;
	if (left < 0) {left=0;}
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = left+'px';
	displayDiv.style.top = top+'px';
	displayDiv.style.display='block';
}
*/
function displayImageFromButton(buttonObj,title,image,imgWidth) {
	if (!imgWidth) {imgWidth=700;}
	var imageStrg="<TABLE cellpadding=0 cellspacing=0><TR><TH colspan=2 class=\\"title2\\">"+title+"&nbsp;&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></TH></TR>";
	imageStrg+="<TR><TD><IMG src='$graphDirHtml/"+image+"' width="+imgWidth+"></TD></TR></TABLE>";
	var infoDiv=document.getElementById('infoDIV').innerHTML=imageStrg;
	var [top,left]=getElementPosition(buttonObj);
	top+=20;
	left-=350;
	if (left < 0) {left=0;}
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = left+'px';
	displayDiv.style.top = top+'px';
	displayDiv.style.display='block';
}
function displayRatioStat(linkID,ratioCode,ratioName) {
	var ratioStrg='Ratio: '+ratioName;
	if (ratioCode.match(/^-/)) { // reverse ratio
		ratioCode=ratioCode.replace('-',''); // true ratioCode (replace w/o regExp -> 1st '-' only)
		var rName=ratioName.split('/');
		ratioName=rName[1]+'/'+rName[0]+' (Reversed)';
		ratioStrg='Reverse ratio: '+rName[1]+'/'+rName[0]
	}
	var ratioLink=document.getElementById(linkID);
	var imageStrg="<TABLE cellpadding=0 cellspacing=0><TR><TH colspan=2 class=\\"title2\\">"+ratioStrg+"'&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></TH></TR>";

	if ('$ratioType'=='Ratio') {
		imageStrg+="<TR><TH>Boxplot</TH><TH>Distribution</TH></TR>";
		imageStrg+="<TR><TD><IMG src=\\"$graphDirHtml/prot_ratio_box_"+ratioCode+".png\\"></TD><TD><IMG src=\\"$graphDirHtml/prot_ratio_hist_"+ratioCode+".png\\"></TD></TR>";
	}
	else if ($algoVersion>=3) {
		imageStrg+="<TR><TH colspan=2><BR>Fold Change distribution<BR><IMG src=\\"$graphDirHtml/distriblog2FC_"+ratioCode+".jpeg\\" width=600></TH></TR>";
		imageStrg+="<TR><TH colspan=2><BR>$pValueType distribution<BR><IMG src=\\"$graphDirHtml/distribPValue_"+ratioCode+".jpeg\\" width=600></TH></TR>";
	}
	else { // (Super/Simple)Ratio
		if ('$labelType'=='FREE') {
			imageStrg+="<TR><TH colspan=2>Distribution</TH></TR>";
			imageStrg+="<TR><TD colspan=2><IMG src=\\"$graphDirHtml/Graph_Box_Hist_ratioprot_"+ratioCode+".png\\"></TD></TR>";
		}
		else {
			if (ratioCode.match('-') \|\| '$biasCorrFlag'=='false') { // ratio of ratio or no normalization
				imageStrg+="<TR><TD colspan=2><IMG src=\\"$graphDirHtml/Graph_Box_Hist_ratioprot_"+ratioCode+".png\\"></TD></TR>";
			}
			else { // primary ratio with normalization
				// imageStrg+="<TR><TD><TABLE cellpadding=0 cellspacing=0><TR><TH>Before normalization</TH><TH>After normalization</TH></TR><TR><TD><IMG src=\\"$graphDirHtml/MAplotPep_"+ratioCode+"Before.png\\" width=240></TD><TD><IMG src=\\"$graphDirHtml/MAplotPep_"+ratioCode+"After.png\\" width=240></TD></TR></TABLE></TD><TD><IMG src=\\"$graphDirHtml/Graph_Box_Hist_ratioprot_"+ratioCode+".png\\"></TD></TR>";
				imageStrg+="<TR><TH>Before normalization<BR><IMG src=\\"$graphDirHtml/MAplotPep_"+ratioCode+"Before.png\\" width=320><BR>After normalization<BR><IMG src=\\"$graphDirHtml/MAplotPep_"+ratioCode+"After.png\\" width=320></TH><TD valign=\\"top\\"><IMG src=\\"$graphDirHtml/Graph_Box_Hist_ratioprot_"+ratioCode+".png\\"></TD></TR>";
			}
		}
	}
	imageStrg+="</TABLE>";
	var infoDiv=document.getElementById('infoDIV').innerHTML=imageStrg;
	var [top,left]=getElementPosition(ratioLink);
	top+=20;
	left-=500;
	if (left < 0) {left=0;}
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = left+'px';
	displayDiv.style.top = top+'px';
	displayDiv.style.display='block';
}
|;
	}
	if ($action eq 'summary') {
		print qq
|function resumeQuantification() {
	if (confirm("WARNING: Resuming a failed quantification will succeed only if failure was due to a non data-dependent cause (eg system failure).\\n                   Proceed?")) {
		window.location="$promsPath{cgi}/startDesignQuantification.cgi?RESUME=$selQuantifID&PROJ_ID=$projectID&ID=$designID";
	}
}
</SCRIPT>
<BR>
|;
	}
	else {
		print qq
|var filters={
	all:['foldChgRawFONT1','foldChgType','foldChgRawFONT2','foldChgType','foldChange','pValueFONT','pValue','coefVarFONT1','coefVar','coefVarFONT2'],
	graph:['coefVarFONT1','coefVar','coefVarFONT2'],
	list:['foldChgRawFONT1','foldChgType','foldChgRawFONT2','foldChgType','foldChange','pValueFONT','pValue','coefVarFONT1','coefVar','coefVarFONT2']
};
function updateFiltersClass(view) {
	for (let i=0; i<filters.all.length; i++) {
		let filter=document.getElementById(filters.all[i]);
		if (filter) filter.className='noFilter';
	}
	for (let i=0; i<filters[view].length; i++) {
		let filter=document.getElementById(filters[view][i]);
		if (filter) filter.className='trueFilter';
	}
	let fullRatioSpan=document.getElementById('fullRatioSPAN');
	if (fullRatioSpan) { // Missing if not ratio-quantif
		fullRatioSpan.style.display=(view==='graph')? '' : 'none';
		let allowInf=document.getElementById('allowInfSPAN');
		if (allowInf) {allowInf.style.display=(view==='graph')? 'none' : '';}
	}
	alterUpdateButton();
}
function updateReplicStatus(pepType) {
	let numPepSc=document.displayForm.numPepScope;
	if (!numPepSc) return; // no scope available
	let pepReplicSpan=document.getElementById('pepReplicSPAN'); // Replicate data?
	let replicDisab=false; //default
	if (pepType.match('msms')) { // truely identified peptides selected
		if (pepReplicSpan) {pepReplicSpan.style.display='none';}
		replicDisab=true; // disable replicate filtering because not compatible with true peptide selection)
	}
	let scOptLength=numPepSc.options.length;
	for (let i=1; i<scOptLength; i++) {
		numPepSc.options[i].disabled=(numPepSc.options[i].value.match('replicate'))? replicDisab : !replicDisab; // mutually exclusive with true pep.
	}
	if (numPepSc.options[numPepSc.selectedIndex].disabled) {
		numPepSc.selectedIndex=0; // set to 'ratio' if current has bee disabled
		if (pepReplicSpan) {pepReplicSpan.style.display='none';}
	}
}
function alterUpdateButton() {
	document.getElementById('submitBUT').className=($dispResults)? 'trueFilter' : 'noFilter';
}
function updateSortBy() {
	if ('$view'=='graph') {
		if (document.getElementById('listFromGraph') && document.getElementById('protListDIV').style.display != 'none') { // a protein list is displayed from graph
			ajaxListSelectedProteins(lastListFromGraph[0],lastListFromGraph[1]); // re-fetch list to update sort order
		}
	}
	else {document.displayForm.submit();} // view=list => resubmit form
}
function updatedListDisplay() {
	var foldChange=document.getElementById('foldChange').value;
	if (foldChange < 1) {
		alert('ERROR: Fold change must be >= 1!');
		return;
	}
/*
	var stdDev=document.getElementById('stdDev').value;
	if (stdDev < 1) {
		alert('ERROR: Standard deviation must be >= 1!');
		return;
	}
*/
	var pValue=document.getElementById('pValue').value;
	if (pValue < 0) {
		alert('ERROR: $pValueType must be >= 0!');
		return;
	}
	var sort=document.getElementById('sort').value;
	window.location="$promsPath{cgi}/showProtQuantification.cgi?CALL=$call&id_ana=$analysisID&id_quantif=$selQuantifID&ACT=$action&view=$view&foldChange="+foldChange+"&pValue="+pValue+"&sort="+sort; //+"&stdDev="+stdDev
}
var selGoID=0;
var goAspects={};
|;
		foreach my $goID (keys %{$refGoAnalyses}) {
			print "goAspects[$goID]='$refGoAnalyses->{$goID}[1]';\n";
		}
		print qq
|function updateGoAspects(goID) {
	var selAspect=document.getElementById('goAspect');
	selAspect.options.length=1;
	document.getElementById('goTermsDIV').innerHTML='';
	if (goID==0) {
		selAspect.style.visibility='hidden';
	}
	else {
		var i=1;
		if (goAspects[goID].match('P')) {selAspect.options[i]=new Option('Biological process',goID+',P'); i++;}
		if (goAspects[goID].match('C')) {selAspect.options[i]=new Option('Cellular component',goID+',C'); i++;}
		if (goAspects[goID].match('F')) {selAspect.options[i]=new Option('Molecular function',goID+',F');}
		selAspect.style.visibility='visible';
	}
}
function checkAllProteins(chkStatus) {
	var checkBoxList=document.protForm.chkProt;
	if (checkBoxList.length) {
		for (var i=0; i < checkBoxList.length; i++) {checkBoxList[i].checked=chkStatus;}
	}
	else {checkBoxList.checked=chkStatus;}
}

function checkProtToExport(required=false) {
	var protChkBox=document.getElementsByName('chkProt');
	if (!protChkBox.length) {protChkBox=[protChkBox];} // force to array if single element
	var filterProtIDs = [];
	protChkBox.forEach(function(chkBox) {
		if(chkBox.checked) {
			filterProtIDs.push(chkBox.value);
		}
	});
	
	if (filterProtIDs.length) {
		exportProteins(filterProtIDs);
	} else if(required) {
		alert('You must select the proteins/sites to export.');
	} else {
		exportProteins();
	}
}

function exportProteins(filterProtIDs=[]) {
	document.displayForm.filterProtIDs.value = '';
	if (filterProtIDs.length) {
		document.displayForm.filterProtIDs.value = filterProtIDs.join(',');
	}
	var defAction=document.displayForm.ACT.value;
	document.displayForm.ACT.value='export';
	document.displayForm.submit();
	document.displayForm.ACT.value=defAction;
}
|;
# function exportIsoforms() {
# 	//TODO: Check for selected proteins
# 	document.protForm.action='../sliva/createPhosphoSequence.cgi';
# 	document.protForm.submit();
# 	document.protForm.action='./showProtQuantification.cgi';
# }
# |;
		my $dispTargetPosStrg=($ratioType eq 'None')? join(',',sort{$a<=>$b} keys %dispStates) : join(',',sort{abs($a)<=>abs($b)} keys %dispRatios);
		if ($view eq 'graph') {
			#my $maxRatioAllowed=($ratioType eq 'Ratio')? scalar @{$refInfo->{RATIOS}} : (scalar @{$refInfo->{STATES}})-1;
			print qq
|//datasetIdx to targetPos
var dataset2TargetPos=[$dispTargetPosStrg];
var hColorIdx=0;
var usedHighLights={},usedColorIdx=[];
for (let idx=0; idx<hColors.length; idx++) {usedColorIdx[idx]=0;}
function updateUsedColors(hlName,action) { // callback for chart highlighting
	if (action != 'delete') return;
	var idx=usedHighLights[hlName];
	delete usedHighLights[hlName];
	usedColorIdx[idx]--;
	if (usedColorIdx[idx]==0) {hColorIdx=idx;}
}
// AJAX --->
function proteinLink(dSetIdx,identifier,protID) {
	ajaxProteinData(null,protID,dataset2TargetPos[dSetIdx],'$ajaxPepAction');
}
var lastListFromGraph=[]; // needed to recall ajaxListSelectedProteins if "Sort by" is changed in graph view
function ajaxListSelectedProteins(selectedPoints,thresholds) {
	lastListFromGraph=[selectedPoints,thresholds];
	saveListGlobals.themesFetched=false; // used for themes management
	document.getElementById('displayDIV').style.display='none'; // in case popup is displayed
	// Adjusting protListDIV position to graph (graph can move due to summary option display)
	var graphDiv=document.getElementById('mainGraphDIV');
	var divX = getDomOffset(graphDiv,'offsetLeft') + graphDiv.offsetWidth + 2;
	var divY = getDomOffset(graphDiv,'offsetTop');
	var listDiv=document.getElementById('protListDIV');
	listDiv.style.left = divX + 'px';
	listDiv.style.top = divY + 'px';
	// Updating contents
	listDiv.innerHTML="<BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<FONT class=\\"title3\\">Fetching data. Please wait...</FONT><BR>&nbsp;&nbsp;&nbsp;&nbsp;<IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR><BR>";
	listDiv.style.display='block';
	//Parameters (extracting no-duplicates list of proteins)
	var noDupProteins=new Object();
	for (var gr in selectedPoints) {
		for (var i=0; i<selectedPoints[gr].length; i++) {
			noDupProteins[selectedPoints[gr][i]]=1;
		}
	}
	var paramStrg='ACT=ajaxListProt&id_ana=$anaID&id_quantif=$selQuantifID&view=$view&dispTargetPosAjax=$dispTargetPosStrg&pepType=$dispPepType&numPep=$dispNumPep&numPepScope=$numPepScope&numReplic=$dispNumReplic&sort='+document.getElementById('sort').value;
	if (document.displayForm.dispMeasure) { // for MaxQuant intensities or Proteomic Ruler
		paramStrg+='&dispMeasure='+document.displayForm.dispMeasure.value+'&selProt=';
	}
	else { // Ratio
		paramStrg+='&foldChange=$dispFoldChange&pValue=$dispPvalue&coefVar=$dispCV&selProt='; //&stdDev=\$dispStdDev
	}
	var p1=true;
	for (var prot in noDupProteins) {
		if (!p1) {paramStrg+=',';}
		paramStrg+=encodeURIComponent(prot); // because of some modProtID
		p1=false;
	}
//paramStrg+='&thresholds='+thresholds.join(','); // no longer used
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("POST","$promsPath{cgi}/showProtQuantification.cgi",true);
	//Send the proper header information along with the request
	XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			listDiv.innerHTML='<INPUT type="hidden" id="listFromGraph" value="1"/>'+XHR.responseText; // listFromGraph: flag
		}
	}
	XHR.send(paramStrg);
}
function ajaxListDecorateGraph(listID) {
	if (listID==0) return;
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	var noSiteParam=($isModifQuantif==0)? '&noSite=1' : '&noSite=0'; // 1: asks for a list of protIDs without sites
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxListDecGraph"+noSiteParam+"&listID="+listID,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			var listData=XHR.responseText.split('::'); // list_name::hasPTMs(0/1)::prot1,prot2,...,protN
			var matchPattern=(listData[1]=='0' && $isModifQuantif==1)? '^###-' : null;
			addHighlighting(mainChart,listData[0],hColors[hColorIdx],{'-1':listData[2].split(',')},matchPattern);
			usedHighLights[listData[0]]=hColorIdx; usedColorIdx[hColorIdx]++;
			hColorIdx++;
			if (hColorIdx==hColors.length) hColorIdx=0;
		}
	}
	XHR.send(null);
}
function ajaxGoDecorateGraph(termIdStrg) {
	if (!termIdStrg) return;
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxTermProt&projectID=$projectID&goStrg="+termIdStrg,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			var termData=termIdStrg.split(',');
			var matchPattern=($isModifQuantif==1)? '^###-' : null;
			addHighlighting(mainChart,termData[2],hColors[hColorIdx],{'-1':XHR.responseText.split(';')},matchPattern); // \$highlighMatchStrg
			usedHighLights[listData[0]]=hColorIdx; usedColorIdx[hColorIdx]++;
			hColorIdx++;
			if (hColorIdx==hColors.length) hColorIdx=0;
		}
	}
	XHR.send(null);
}
function ajaxSearchConvertIdentifier(graphSearch,graphSearchArgs,searchTextIdx) { // (graph lib search function,array of function arguments,index of search text in array). To be called at end of convertion function
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxConvIdent&projectID=$projectID&quantifList=$selQuantifID&id_ana=$analysisID&TEXT="+encodeURIComponent(graphSearchArgs[searchTextIdx]),true); //search text
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			graphSearchArgs[searchTextIdx]=XHR.responseText; // replace old text with converted one
			graphSearch(...graphSearchArgs); // call graph lib search fuction & convert array to arguments
		}
	};
	XHR.send(null);
}
|;
		}
		else {print "// AJAX --->\n";}
		#my $maxRatioAllowed=($ratioType eq 'Ratio')? scalar @{$refInfo->{RATIOS}} : (scalar @{$refInfo->{STATES}})-1;
		print qq
|var isNav = (navigator.appName.indexOf("Netscape") !=-1);
function ajaxProteinData(e,protID,ratioPos,action) {
	var displayDiv=document.getElementById('displayDIV');
	var divX,divY;
	if (e) { // called from protein list
		divX = (isNav)? e.pageX : event.clientX + document.body.scrollLeft; //divX-=5;
		divY = (isNav)? e.pageY : event.clientY + document.body.scrollTop; //divY+=10;
	}
	else {
		var graphDiv=document.getElementById('mainGraphDIV');
		divX = getDomOffset(graphDiv,'offsetLeft') + graphDiv.offsetWidth +2;
		divY = getDomOffset(graphDiv,'offsetTop');
	}
	displayDiv.style.left = divX + 'px';
	displayDiv.style.top = divY + 'px';
	displayDiv.style.display='block';

	var infoDiv=document.getElementById('infoDIV');
	infoDiv.innerHTML="<BR><BR><IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR>";

	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	//var extraParams=(action=='$ajaxPepAction')? '&pepType=$dispPepType' : '';
	//XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT="+action+"&id_ana=$anaID&id_quantif=$selQuantifID&view=$view&sort=$dispSort&foldChange=$dispFoldChange&pValue=$dispPvalue&coefVar=$dispCV&pepType=$dispPepType&numPep=$dispNumPep&id_prot="+encodeURIComponent(protID)+"&ratio="+ratioPos,true); //+extraParams    &stdDev=\$dispStdDev
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT="+action+"&id_ana=$anaID&id_quantif=$selQuantifID&pepType=$dispPepType&id_prot="+encodeURIComponent(protID)+"$dispMeasParamStrg&dispTargetPosAjax=$dispTargetPosStrg&ratio="+ratioPos,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			if (action=='$ajaxPepAction') {
				var codeParts=XHR.responseText.split('#==========#');
				infoDiv.innerHTML=codeParts[0]; // HTML part
				eval(codeParts[1]); // javascript part
			}
			else {infoDiv.innerHTML=XHR.responseText;}
		}
	}
	XHR.send(null);
}
function ajaxUpdateGoTermList(goIdStrg) {
	var termsDiv=document.getElementById('goTermsDIV');
	if (!goIdStrg) {
		termsDiv.innerHTML='';
		return;
	}
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxGoTerms&projectID=$projectID&goStrg="+goIdStrg,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			termsDiv.innerHTML=XHR.responseText;
		}
	}
	XHR.send(null);
}
|;
		&promsMod::printAjaxManageSaveProteins($projectID,\%promsPath,'document.protForm.chkProt','ajaxUpdateRestrict');
		print qq
|function ajaxUpdateRestrict() {
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxRestrictList&projectID=$projectID&restrictList=$restrictListID&listAction=$listAction&noSelect=1",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			document.getElementById('restrictSPAN').innerHTML='<SELECT name="restrictList" class="trueFilter" onchange="alterUpdateButton()"><OPTION value="">-= Select =-</OPTION>'+XHR.responseText+'</SELECT>';
			if ('$view'=='graph' && $dispResults==1) {
				document.getElementById('listDecoSPAN').innerHTML='<SELECT id="listDeco" onchange="ajaxListDecorateGraph(this.value)"><OPTION value="0">-= None =-</OPTION>'+XHR.responseText;
			}
		}
	}
	XHR.send(null);
}
</SCRIPT>
<BR>
<FORM name="displayForm" method="POST">
<INPUT type="hidden" name="ACT" value="$action">
<INPUT type="hidden" name="CALL" value="$call">
<INPUT type="hidden" name="id_ana" value="$anaID">
<INPUT type="hidden" name="id_quantif" value="$selQuantifID">
<INPUT type="hidden" name="displayRes" value="1">
<INPUT type="hidden" name="filterProtIDs" value="">
|;
	}
	my $statusImage=($quantifStatus==-2)? 'lightRed1.gif' : ($quantifStatus==-1)? 'lightGray1.gif' : ($quantifStatus==0)? 'lightYellow1.gif' : 'lightGreen1.gif';
	#my $colspanStrg=($ratioType eq 'None')? '' : ($quantifSoftware eq 'MaxQuant')? ' colspan=2' : ' colspan=3'; #($existReplicates)? ' colspan=3' : ' colspan=2';
	my $colspanStrg=($ratioType=~/Ratio/)? ' colspan=3' : ($quantifType=~/^(MQ|PROT_ABUNDANCE)$/)? ' colspan=2' : '';
	print qq
|<TABLE bgcolor="$darkColor" border=0>
<TR><TH nowrap align="right">Labeling :</TH><TD bgcolor="$lightColor"$colspanStrg valign=top>&nbsp;<B>$labelingName{$labelType}$swathStrg</B></TD></TR>
|;
	#my $maxDefaultRatios=$numStates-1; # all ratios with 1st reference
	#unless (scalar keys %dispRatios) { # setting default displayed ratios
	#	foreach my $pos (1..$maxDefaultRatios) {$dispRatios{$pos}=1;}
	#}

	my $colSpanCorrel=($quantifSoftware eq 'MaxQuant' || $quantifStatus > 0)? 1 : 2; # +/- correlation buttons
	my $stateLabel=($numStates==1)? 'State' : 'States';
	print qq
|<TR><TH nowrap align="right" valign="top">$stateLabel :</TH><TD bgcolor="$lightColor" colspan=$colSpanCorrel valign="top" style="min-width:600px">
<DIV style="max-height:150px;overflow:auto"><TABLE bgcolor="$darkColor" width=100%>
|;
	foreach my $statePos (sort{$a<=>$b} keys %stateInfo) {
		#print "<TR><TH nowrap align=right>&nbsp;State #$statePos:</TH><TD bgcolor=$lightColor width=450><B>$stateInfo{$statePos}{NAME}</B>$stateInfo{$statePos}{REPLIC}</TD></TR>\n";
		my $replicateNames=($stateInfo{$statePos}{'NAME_LONG'})? $stateInfo{$statePos}{'NAME_LONG'} : '';
		print "<TR><TH nowrap align=right>&nbsp;#$statePos:</TH><TD bgcolor=$lightColor width=100%>";
		if ($ratioType eq 'None') {
			if ($action eq 'summary') {print "&nbsp;<IMG src=\"$promsPath{images}/$statusImage\"/>";}
			elsif ($numStates > 1) {
				unless (scalar keys %dispStates) {%dispStates=(1=>1);} # set default
				my $selStatus=($dispStates{$statePos})? ' checked' : '';
				print "<INPUT type=\"checkbox\" name=\"dispStates\" value=\"$statePos\"$selStatus>";
			}
			else {print "<INPUT type=\"hidden\" name=\"dispStates\" value=\"$statePos\">"}
		}
		print "&nbsp;";
		print "<A href=\"javascript:void(null)\" onmouseover=\"popup('$stateInfo{$statePos}{POPUP}')\" onmouseout=\"popout()\">" if $stateInfo{$statePos}{'POPUP'};
		print "<B>$stateInfo{$statePos}{NAME}$replicateNames</B>";
		print "</A>" if $stateInfo{$statePos}{'POPUP'};
		print $stateInfo{$statePos}{REPLIC};

		#if ($quantifStatus > 0 && ($ratioType=~/S\w+Ratio/ || $quantifType eq 'PROT_ABUNDANCE') && $isPepIntensity && $quantifSoftware eq 'myProMS') { #}
		if ($quantifStatus > 0 && $quantifSoftware eq 'myProMS' && $isPepIntensity) {
			my $imageFile=($algoVersion >= 3)? "distriblog2FC_State$statePos.jpeg" : "Graph_Box_Hist_ratioprot_mean_State$statePos.png";
			print "&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Mean distrib.\" onclick=\"displayImageFromButton(this,'Protein Mean Distribution for State $stateInfo{$statePos}{NAME}','$imageFile')\">" if -e "$resultDir/graph/$imageFile";
		}
		print "&nbsp;</TD></TR>\n";
	}
	print "</TABLE></DIV></TD>\n";

	if ($quantifStatus > 0 && $quantifType eq 'PROT_ABUNDANCE') {
		my $disabCorrStrg=(-e "$resultDir/$correlMatFileName")? '' : 'disabled';
		print qq
|<TD bgcolor="$lightColor" align="center" valign="middle"><INPUT type="button" class="font11" value="Missing values" style="width:120px" onclick="displayMissingValues(this)"><BR>
<INPUT type="button" class="font11" value="Global correlation" style="width:120px" onclick="displayReplicCorrelation(this)" $disabCorrStrg>
|;

		if ($refInfo->{'ABUND_MEASURES'} && grep(/aLFQ.*(?:amounts|calibration)/, @{$refInfo->{'ABUND_MEASURES'}})) {
			print qq
|<BR>
<INPUT type="button" class="font11" value="Protein quantities" style="width:120px" onclick="displayProteinQuantities(this)">
</TD>
|;
		}
		else {print "</TD>";}
	}
	elsif ($refInfo->{RATIOS}) {
		if ($quantifStatus > 0 && $quantifSoftware ne 'MaxQuant') {
			print qq |<TD bgcolor=$lightColor valign=middle>|;
			if ($ratioType=~/S\w+Ratio/) {
				#my $disabStrg=(-e "$resultDir/graph/$correlImgFileName")? '' : 'disabled';
				my $disabCorrStrg=(-e "$resultDir/$correlMatFileName")? '' : 'disabled';
				print qq
|&nbsp;<INPUT type="button" class="font11" value="Missing values" style="width:120px" onclick="displayMissingValues(this)">&nbsp;<BR>
&nbsp;<INPUT type="button" class="font11" value="Global correlation" style="width:120px" onclick="displayReplicCorrelation(this)" $disabCorrStrg>&nbsp;<BR>
|;
			}
			print qq |&nbsp;<INPUT type="button" class="font11" value="Global coeff. var." style="width:120px" onclick="displayImageFromButton(this,'Coefficient of variation between all replicates','pep_cv_hist.png')">&nbsp;| if -e "$resultDir/graph/pep_cv_hist.png"; # before myProMS v3! ($ratioType eq 'Ratio' || $supRatioReplicates); #displayPepVariation(this)
			print "</TD>\n";
		}
		if (scalar @{$refInfo->{RATIOS}} > 1) { # multiple ratios
			unless (scalar keys %dispRatios) {%dispRatios=(1=>[],2=>[],3=>[]);} # set default
			print "<TH bgcolor=$lightColor valign=top><TABLE bgcolor=$darkColor cellpadding=1><TR>";
			if ($action eq 'summary') {
				print "<TH bgcolor=$lightColor nowrap>Ratios computed:</TH><TH>St.</TH>";
				foreach my $stPosT (2..$numStates) {print "<TH>#$stPosT</TH>"}
				print "</TR>\n";
				#my $numReferences=(($ratioType eq 'SimpleRatio' && $labelType ne 'FREE' && $quantifSoftware ne 'MaxQuant') || ($refInfo->{'SINGLE_REF'} && $refInfo->{'SINGLE_REF'}[0]==1))? 1 : $numStates-1;
				my $numReferences=($refInfo->{'SINGLE_REF'} && $refInfo->{'SINGLE_REF'}[0]==1)? 1 : $numStates-1;
				my $ratioPos=0;
				foreach my $stPosR (1..$numReferences) {
					print "<TR>";
					if ($stPosR==1) {
						my $extraRefStrg.=($numReferences >= 4)? '<BR>e&nbsp;<BR>r&nbsp;<BR>e&nbsp;<BR>n&nbsp;<BR>c&nbsp;<BR>e&nbsp;<BR>s&nbsp;' : '';
						print "<TH bgcolor=\"$lightColor\" valign=\"middle\" rowspan=$numReferences align=right><FONT style=\"line-height:90%\">R&nbsp;<BR>e&nbsp;<BR>f&nbsp;$extraRefStrg</FONT></TH>";
					}
					print "<TH>#$stPosR</TH>";
					foreach my $stPosT (2..$numStates) {
						if ($stPosT<=$stPosR) {print "<TH bgcolor=$lightColor></TH>"; next;}
						$ratioPos++;
						my $ratioTag=($ratioType eq 'SuperRatio' && $ratioPos > $numReferences)? $encodedDegStrg : '';
						my $ratioName="$stateInfo{$stPosT}{NAME}$ratioTag/$stateInfo{$stPosR}{NAME}$ratioTag";
						print "<TH bgcolor=$lightColor onmouseover=\"popup('<B>$ratioName</B>')\" onmouseout=\"popout()\"><IMG id=\"menuRatio_$ratioPos\" src=\"$promsPath{images}/$statusImage\"";
						print " onclick=\"displayRatioStat('menuRatio_$ratioPos','$ratioCode{$ratioPos}','$ratioName')\"" if ($quantifSoftware eq 'myProMS' && $quantifStatus > 0);
						print "></TH>\n";
					}
					print "</TR>\n";
				}
			}
			elsif ($refInfo->{'SINGLE_REF'} && $refInfo->{'SINGLE_REF'}[0]==1) { # all vs ref
				print "<TH bgcolor=$lightColor nowrap colspan=$numStates>Ratios displayed</TH></TR><TR><TH bgcolor=$lightColor nowrap rowspan=2 align=right><FONT color=\"#00BB00\">&nbsp;Normal</FONT></TH>";
				my $colspan=$numStates-1;
				my $ratioPos=0;
				foreach my $stPosT (2..$numStates) {
					$ratioPos++;
					print "<TH>&nbsp;#$ratioPos&nbsp;";
					if ($quantifSoftware eq 'myProMS') {
						my $ratioName="$stateInfo{$stPosT}{NAME}/$stateInfo{1}{NAME}";
						print "<BR><IMG id=\"menuRatio_$ratioPos\" src=\"$promsPath{images}/plus.gif\" onclick=\"displayRatioStat('menuRatio_$ratioPos','$ratioCode{$ratioPos}','$ratioName')\">";
						print "</TH>\n";
					}
				}
				print "</TR><TR>";
				$ratioPos=0;
				foreach my $stPosT (2..$numStates) {
					$ratioPos++;
					my $ratioName="$stateInfo{$stPosT}{NAME}/$stateInfo{1}{NAME}";
					print "<TH bgcolor=\"#00DD00\" onmouseover=\"popup('<B>$ratioName</B>')\" onmouseout=\"popout()\"><INPUT type=\"checkbox\" name=\"dispRatios\" id=\"$ratioPos\" value=\"$ratioPos\" onclick=\"alterUpdateButton(); if (this.checked) {document.getElementById('-$ratioPos').checked=false;}\"";
					print ' checked' if $dispRatios{$ratioPos};
					print "/></TH>\n";
				}
				
				print "</TR>\n<TR><TH bgcolor=$lightColor nowrap align=right valign=top nowrap><FONT color=\"#DD0000\">&nbsp;Reversed</FONT></TH>";
				my $revRatioPos=0;
				foreach my $stPosT (2..$numStates) {
					$revRatioPos--;
					my $ratioName="$stateInfo{1}{NAME}/$stateInfo{$stPosT}{NAME}";
					print "<TH bgcolor=\"#DD0000\" onmouseover=\"popup('<B>$ratioName</B>')\" onmouseout=\"popout()\"><INPUT type=\"checkbox\" name=\"dispRatios\" id=\"$revRatioPos\" value=\"$revRatioPos\" onclick=\"alterUpdateButton(); if (this.checked) {document.getElementById((-1*$revRatioPos)+'').checked=false;}\"";
					print ' checked' if $dispRatios{$revRatioPos};
					print "/></TH>\n";
				}
				print "</TR>\n";
			}
			else { # all vs all
				print "<TH bgcolor=$lightColor nowrap>Ratios displayed:</TH><TH>St.</TH>";
				foreach my $stPosT (1..$numStates) {print "<TH>#$stPosT</TH>"}
				print "</TR>\n";
				my $ratioPos=0;
				my %revRatioPos;
				foreach my $stPosR (1..$numStates) {
					print "<TR>";
					if ($stPosR==1) {print "<TH bgcolor=$lightColor rowspan=$numStates align=right valign=top nowrap><FONT color=\"#00BB00\">&nbsp;Normal</FONT>&nbsp;<BR><FONT color=\"#DD0000\">&nbsp;Reversed</FONT>&nbsp;</TH>";}
					print "<TH>#$stPosR</TH>";
					foreach my $stPosT (1..$numStates) {
						if ($stPosT==$stPosR) {print "<TH bgcolor=$lightColor></TH>\n";}
						elsif ($stPosT > $stPosR) {
							$ratioPos++;
							$revRatioPos{"$stPosT:$stPosR"}=-$ratioPos;
							my $ratioTag=($ratioType eq 'SuperRatio' && $stPosR > 1)? $encodedDegStrg : '';
							my $ratioName="$stateInfo{$stPosT}{NAME}$ratioTag/$stateInfo{$stPosR}{NAME}$ratioTag";
							print "<TH bgcolor=\"#00DD00\" onmouseover=\"popup('<B>$ratioName</B>')\" onmouseout=\"popout()\"><INPUT type=\"checkbox\" name=\"dispRatios\" id=\"$ratioPos\" value=\"$ratioPos\" onclick=\"alterUpdateButton(); if (this.checked) {document.getElementById('-$ratioPos').checked=false;}\"";
							print ' checked' if $dispRatios{$ratioPos};
							print "/>";
							print "<IMG id=\"menuRatio_$ratioPos\" src=\"$promsPath{images}/plus.gif\" onclick=\"displayRatioStat('menuRatio_$ratioPos','$ratioCode{$ratioPos}','$ratioName')\">" if $quantifSoftware eq 'myProMS';
							print "</TH>\n";
						}
						else { # reverse ratio
							my $revRatPos=$revRatioPos{"$stPosR:$stPosT"};
							my $ratioTag=($ratioType eq 'SuperRatio' && $stPosT > 1)? $encodedDegStrg : '';
							my $ratioName="$stateInfo{$stPosT}{NAME}$ratioTag/$stateInfo{$stPosR}{NAME}$ratioTag";
							print "<TH bgcolor=\"#DD0000\" onmouseover=\"popup('<B>$ratioName</B>')\" onmouseout=\"popout()\"><INPUT type=\"checkbox\" name=\"dispRatios\" id=\"$revRatPos\" value=\"$revRatPos\" onclick=\"alterUpdateButton(); if (this.checked) {document.getElementById((-1*$revRatPos)+'').checked=false;}\"";
							print ' checked' if $dispRatios{$revRatPos};
							print "/></TH>\n";
						}
					}
					print "</TR>\n";
				}
			}
			print "</TABLE>\n";
		}
		else { # single ratio
			print "<TD bgcolor=$lightColor class=\"title3\">";
			if ($action eq 'summary') {
				print "&nbsp;<IMG src=\"$promsPath{images}/$statusImage\">";
				if ($quantifSoftware eq 'myProMS') {
					print "<A id=\"menuRatio_1\" href=\"javascript:displayRatioStat('menuRatio_1','$ratioCode{1}','$stateInfo{2}{NAME}/$stateInfo{1}{NAME}')\">$stateInfo{2}{NAME}/$stateInfo{1}{NAME}</A>&nbsp;";
				}
				else {print "$stateInfo{2}{NAME}/$stateInfo{1}{NAME}&nbsp;";}
			}
			else {
				print "<TABLE cellspacing=0 cellpadding=0 ><TR><TD align=right class=\"title3\">&nbsp;<FONT color=\"#00BB00\">Computed ratio:</FONT></TD><TD class=\"title3\"><INPUT type=\"checkbox\" name=\"dispRatios\" id=\"1\" value=\"1\" onclick=\"if (this.checked) {document.getElementById('-1').checked=false;}\"";
				print ' checked' if $dispRatios{1};
				print "/>";
				if ($quantifSoftware eq 'myProMS') {
					print "<A id=\"menuRatio_1\" href=\"javascript:displayRatioStat('menuRatio_1','$ratioCode{1}','$stateInfo{2}{NAME}/$stateInfo{1}{NAME}')\">$stateInfo{2}{NAME}/$stateInfo{1}{NAME}</A>&nbsp;";
				}
				else {print "$stateInfo{2}{NAME}/$stateInfo{1}{NAME}&nbsp;";}
				print "</TD></TR><TR><TD align=right class=\"title3\">&nbsp;<FONT color=\"#DD0000\">Reversed ratio:</FONT></TD><TD class=\"title3\"><INPUT type=\"checkbox\" name=\"dispRatios\" id=\"-1\" value=\"-1\" onclick=\"if (this.checked) {document.getElementById('1').checked=false;}\"";
				print ' checked' if $dispRatios{-1};
				print "/>$stateInfo{1}{NAME}/$stateInfo{2}{NAME}&nbsp;</TD></TR></TABLE>\n";
			}
			print "</TD>";
		}
	}
	print "</TR>\n";

	if ($isModifQuantif) {
		print "<TR><TH align=right nowrap valign=top>&nbsp;Fold change correction :</TH><TD nowrap bgcolor=$lightColor$colspanStrg>&nbsp;";
		if ($refInfo->{'INTRA_PROT_REF'}) {
			print "Using protein fold changes in ";
			if ($refInfo->{'INTRA_PROT_REF'}[0] eq '-1') {
				print "<B>current</B> dataset";
				if ($refInfo->{'PEPTIDES_REF'}[0] eq 'manual') {print " (<B>Manual</B> peptide selection <INPUT type=\"button\" class=\"font11\" value=\" List \" onclick=\"ajaxDisplayNormPeptides(this)\">)";}
			}
			else {
				my ($refQuantifID,$rPos,$testRefCondID,$refRefCondID)=split(/[_:]/,$refInfo->{'INTRA_PROT_REF'}[0]);
				my ($refQuantifPath)=$dbh->selectrow_array("SELECT CONCAT(E.NAME,' > ',D.NAME,' > ',Q.NAME) FROM QUANTIFICATION Q,DESIGN D,EXPERIMENT E WHERE Q.ID_DESIGN=D.ID_DESIGN AND D.ID_EXPERIMENT=E.ID_EXPERIMENT AND Q.ID_QUANTIFICATION=$refQuantifID");
				my $sthCondName = $dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
				$sthCondName->execute($testRefCondID);
				my ($testRefCondName)=$sthCondName->fetchrow_array;
				$sthCondName->execute($refRefCondID);
				my ($refRefCondName)=$sthCondName->fetchrow_array;
				$sthCondName->finish;
				print "dataset from <B>$refQuantifPath : $testRefCondName/$refRefCondName</B>";
			}
		}
		else {print 'None';}
		print "</TD></TR>\n";
	}

	#my $pepString=($quantifType !~ /XIC|SIN|EMPAI/)? "&nbsp;&nbsp;&bull;$pepChargeStrg&nbsp;&nbsp;&nbsp;&bull;$pepSourceStrg&nbsp;" : ''; # Change on 29/10/12
	my $pepString=($quantifSoftware eq 'myProMS' && !$protRulerCase)? "&nbsp;&nbsp;&bull;$pepChargeStrg&nbsp;&nbsp;&nbsp;&bull;$pepSourceStrg&nbsp;" : '';
	#my $pValueMeaningStrg=($algoVersion < 3)? '' : ($refInfo->{'RESIDUAL_VAR'}[0] eq 'biological')? ''
	my $infRatioString=($refInfo->{'MIN_INF_RATIOS'} && $refInfo->{'MIN_INF_RATIOS'}[0])? 'Avoided whenever possible' : 'Not avoided';
	$ptmQuantifStrg='<BR>&nbsp;&bull;'.$ptmQuantifStrg if $ptmQuantifStrg;
	if (!$protRulerCase) {	
		print "<TR><TH align=right nowrap valign=top>&nbsp;Peptide selection :</TH><TD nowrap bgcolor=$lightColor$colspanStrg>&nbsp;";
		if ($refInfo->{'PEPTIDES'}) {print "&bull;$pepRangeStrg&nbsp;&nbsp;&nbsp;&bull;$pepMissCutStrg&nbsp;&nbsp;&nbsp;&bull;$pepPtmStrg&nbsp;$pepString&nbsp;$ptmQuantifStrg";}
		else {print 'All';}
		if ($quantifType eq 'PROT_ABUNDANCE' && $refInfo->{'PEPTIDES'}[0] ne 'unique' && !$isModifQuantif) {
			my $MGSharedRule=(!$refInfo->{'MG_SHARED_PEP'} || $refInfo->{'MG_SHARED_PEP'}[0] eq 'best')? 'assigned to <B>best</B> protein' : ($refInfo->{'MG_SHARED_PEP'}[0] eq 'share')? 'used for <B>each</B> protein' : '<B>excluded</B> from dataset'; #'exclude';
			print "<BR>&nbsp;&bull;Peptides shared by multiple Match Groups are $MGSharedRule";
		}
		print "</TD></TR>\n";
	}
	my $protSelectionStrg='<B>All</B> visible proteins'; # default
	my $contaminantLink='without';
	if ($refInfo->{'PROTEINS'}) {
		($protSelectionStrg,$contaminantLink)=($refInfo->{'PROTEINS'}[0] eq 'exclude')? ('<B>Exclude</B> proteins in List ','and') : ('<B>Restrict</B> to proteins in List ','without');
		(my $listID=$refInfo->{'PROTEINS'}[1])=~s/#//;
		my ($listStrg)=$dbh->selectrow_array("SELECT CONCAT(T.NAME,' > ',L.NAME) FROM CATEGORY L,CLASSIFICATION T WHERE L.ID_CLASSIFICATION=T.ID_CLASSIFICATION AND ID_CATEGORY=$listID");
		$protSelectionStrg.="<B>$listStrg</B>";
	}
	if ($refInfo->{'EXCL_CONTAMINANTS'} && $refInfo->{'EXCL_CONTAMINANTS'}[0]) {
		$protSelectionStrg.=" <B>$contaminantLink</B> contaminants";
	}
	$protSelectionStrg.=" ($protQuantifiedStrg)" if $protQuantifiedStrg;
	print qq
|<TR><TH align=right nowrap valign=top>&nbsp;Protein selection :</TH><TD nowrap bgcolor=$lightColor$colspanStrg>&nbsp;$protSelectionStrg</TD></TR>
<TR><TH align=right nowrap valign=top>&nbsp;Quantification settings :|;
	if ($ratioType eq 'Ratio') {
		print qq
|<BR><INPUT type="button" id="moreSettings" class="font11" value="More settings" onclick="updateSettings('more')"/><INPUT type="button" id="lessSettings" class="font11" value="Less settings" style="display:none" onclick="updateSettings('less')"/>|;
	}
	if ($quantifType eq "PROT_RULER") {
		print qq
|</TH><TD bgcolor=$lightColor$colspanStrg valign=top>
	&nbsp;-<B>Quantification method:</B> $quantifMethDesc [$protRulerTypeStrg, based on "$intMetric" values - Software: $softNameStrg$softVersionStrg]<BR>
	&nbsp;-<B>Normalization mode:</B> $protRulerNormStrg<BR>
	&nbsp;-<B>Numerical parameters:</B> $protRulerParamsStrg<BR>
|;
	}
	else {
		#my $quantifSoftwareStrg=($quantifSoftware eq 'myProMS')? 'myProMS v'.$algoVersion : $quantifSoftware;
		print qq
|</TH><TD bgcolor=$lightColor$colspanStrg valign=top>
	&nbsp;-<B>Quantification method:</B> $quantifMethDesc [$ratioTypeStrg$topNString - Software: $softNameStrg$softVersionStrg]<BR>
|;
		if ($refInfo->{'RATIOS'}) {
			if ($algoVersion>=3) {
				my $pvalTypeStrg=($refInfo->{'RESIDUAL_VAR'}[0] eq 'technical')? '<FONT color="#DD0000">'.$refInfo->{'RESIDUAL_VAR'}[0].'</FONT>' : $refInfo->{'RESIDUAL_VAR'}[0];
				print "&nbsp;-<B>$pValueType</B> estimates <B>$pvalTypeStrg</B> variation [Correction: <B>",$pValueCorrections{$refInfo->{'FDR_CONTROL'}[0]},"</B>]<BR>\n";
			}
			if ($quantifSoftware ne 'MaxQuant') {print "&nbsp;-<B>Infinite ratios:</B> $infRatioString\n";}
		}
		elsif ($quantifType eq 'PROT_ABUNDANCE') { # myProMS abundance
			print "&nbsp;-<B>Available measures:</B>";
			my $previous=0;
			foreach my $refMeas (@{$proteinQuantifFamilies{'MEASURES'}{$quantifType}}) {
				my ($measCode,$measName,$isOptional)=@{$refMeas};
				next unless $abundanceMeasures{$measCode};
				print ',' if $previous;
				$measName=~s/ intensity//;
				print " $measName";
				if ($measCode eq 'MY_LFQ') {
					my $lfqParamStrg=" [&ge;$abundanceMeasures{$measCode}[0] pep. ratios, ";
					$lfqParamStrg.=($abundanceMeasures{$measCode}[1])? 'large ratios stabilized]' : 'large ratios <B>not</B> stabilized]';
					print $lfqParamStrg;
				}
				$previous=1;
			}
			print "<BR>\n";
			if (scalar @aLFQparams) {
				my $aLFQmeas=shift(@aLFQparams);
				my ($pqi) = $aLFQmeas=~/aLFQ_(iBAQ|TOP)/;
				my  ($absModel, $pqiNorm, $outType, $topX, $topStrictness, $qtyUnit);
				if ($pqi eq 'TOP') {
					($absModel, $pqiNorm, $outType, $topX, $topStrictness, $qtyUnit) = @aLFQparams;
					$pqi = ucfirst(lc($pqi));
				} else {
					($absModel, $pqiNorm, $outType, $qtyUnit) = @aLFQparams;
				}
				my %unitsMatching = (  # internal code to name displayed
					"mol" 			 => "mol",
					"mmol" 			 => "mmol",
					"copy_number" 	 => "Copy Number",
					"conc_mol_l" 	 => "mol/L or mmol/mL",
					"conc_mmol_l" 	 => "mmol/L",
					"mass_g" 		 => "g",
					"mass_mg" 		 => "mg",
					"mass_ug" 		 => "μg",
					"mass_conc_g_l"  => "g/L or mg/mL or μg/μL",
					"mass_conc_mg_l" => "mg/L or μg/mL",
					"mass_conc_ug_l" => "μg/L"
				);
				$absModel = ucfirst($absModel);
				print "&nbsp;-<B>Parameters for absolute measures:</B>";
				print "&nbsp;&nbsp;PQI: $pqi";
				print "$topX ($topStrictness)" if ($pqi eq "Top");
				print ",&nbsp;Model: $absModel";
				print ",&nbsp;Initial units : $unitsMatching{$qtyUnit}" if ($qtyUnit);
				print "<BR>\n";
			}
		}
		print "<TABLE cellspacing=0><TR valign=top><TD nowrap>&nbsp;-<B>Bias correction:</B></TD>$biasCorrectStrg</TABLE>\n" if $quantifSoftware =~ /myProMS|MSstats|Spectronaut/i;
	}
	if ($ratioType eq 'Ratio') {
		print qq
|	<DIV id="advancedSetDIV" style="display:none">
	&nbsp;&bull;<B><U>Advanced settings:</U></B>
|;
		if ($refInfo->{'THRESHOLD_CV'}) { # replicates used TODO: Make sure THRESHOLD_CV is reported for design quantif
			my $biasCorrStrg='<B>Variation coefficient threshold between replicates:</B> ';
			$biasCorrStrg.=($refInfo->{'THRESHOLD_CV'}[0] eq 'FALSE')? 'Auto' : $refInfo->{'THRESHOLD_CV'}[1];
			print "<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-$biasCorrStrg\n";
		}
		print qq
|	<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-<B>FDR control:</B> $fdrStrg
	<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-<B>$pValueType threshold for outlier detection:</B> $pvalStrg
|;
		print "<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-<B>Alternative hypothesis for comparison:</B> $alterStrg\n" if $refInfo->{'ALTER'};
		print "<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-<B>Confidence interval on protein abundance:</B> $refInfo->{CONFIDENCE_LEVEL}[0]\n" if $refInfo->{'CONFIDENCE_LEVEL'};
		print "</DIV>\n";
	}
	print qq
|</TD></TR>
|;
	if ($action eq 'summary') {
		my $resumeStrg=($quantifStatus==-2 && $quantifType eq 'PROT_ABUNDANCE')? ' <INPUT type="button" value="Resume quantification..." onclick="resumeQuantification()"/>' : '';
		my $statusStrg=($quantifStatus==-2)? '<FONT color="#DD0000">Failed</FONT> (Click on "Monitor Jobs" for more information)'.$resumeStrg : ($quantifStatus==-1)? 'Not launched yet' : ($quantifStatus==0)? 'On-going' : 'Finished';
		$updateDate=&promsMod::formatDate($updateDate);
		print qq
|<TR><TH align=right nowrap valign=top>&nbsp;Creation date :</TH><TD bgcolor=$lightColor$colspanStrg valign=top>&nbsp;$updateDate by <B>$updateUser</B></TD></TR>
<TR><TH align=right nowrap valign=top>&nbsp;Status :</TH><TH align="left" bgcolor=$lightColor$colspanStrg valign=top>&nbsp;$statusStrg</TH></TR>
</TABLE>
<BR>
|;
		$dbh->disconnect;

		if ($quantifStatus==-1 || $quantifStatus==0) { # on-going: reload every 10 sec
			print qq
|<SCRIPT type="text/javascript">
setTimeout(checkQuantifStatus,10000);
</SCRIPT>
|;
		}

		&endHTML;
		exit;
	}
	#my ($selGraphView,$selListView,$listDivVis)=($view eq 'graph')? (' selected','','none') : ('',' selected','block');
	my ($selListView,$foldChPvalClass)=($view eq 'list')? (' selected','trueFilter') : ('','noFilter'); # 'graph' is default
	my ($selFcUp,$selFcDown)=($foldChangeType eq 'up')? (' selected','') : ($foldChangeType eq 'down')? ('',' selected') : ('','');
	$dispPvalue=($dispPvalue)? $dispPvalue : $refInfo->{'FDR_ALPHA'}[0]/100;
	my $hasTruePepCount=($refInfo->{'NUM_TRUE_USED'})? $refInfo->{'NUM_TRUE_USED'}[0] : 0;
	#my ($selDistPepType,$selRazPepType,$selUniPepType)=($numPepCode eq 'DIST_PEP_USED')? (' selected','','') : ($dispPepType eq 'RAZ_UNI_PEP')? ('',' selected','') :  ($dispPepType eq 'UNIQUE_PEP')? ('','',' selected') : ('','','');
	#my ($selDistPepType,$selRazPepType,$selUniPepType,$selRatioCountType)=($numPepCode eq 'DIST_PEP_USED')? (' selected','','','') : ($numPepCode eq 'RAZ_UNI_PEP')? ('',' selected','','') :  ($numPepCode eq 'UNIQUE_PEP')? ('','',' selected','') :  ($numPepCode eq 'NUM_PEP_USED')? ('','','',' selected') : ('','','','');
	my ($selDistPepType,$selNumTruePepType,$selPcTruePepType,$selAllRepPepType,$selDistRepPepType,$selRazPepType,$selUniPepType,$selRatioCountType)=('','','','','','','','');
	if ($numPepCode eq 'DIST_PEP_USED') {$selDistPepType=' selected';}
	elsif ($numPepCode eq 'NUM_TRUE_USED') {
		if ($dispPepType eq 'msmsNum') {$selNumTruePepType=' selected';} else {$selPcTruePepType=' selected';}
	}
	#elsif ($numPepCode eq 'NUM_PEP_USED:REPLIC') {$selAllRepPepType=' selected';}
	#elsif ($numPepCode eq 'DIST_PEP_USED:REPLIC') {$selDistRepPepType=' selected';}
	elsif ($numPepCode eq 'RAZ_UNI_PEP') {$selRazPepType=' selected';}
	elsif ($numPepCode eq 'UNIQUE_PEP') {$selUniPepType=' selected';}
	elsif ($numPepCode eq 'NUM_PEP_USED') {$selRatioCountType=' selected';}
	#$selDistPepType=($dispPepType eq 'distinct')? ' selected' : ''; # 'all' is default

	print qq
|<TR><TH align=right nowrap valign=top>&nbsp;<FONT class="title2">Display :</FONT></TH><TD bgcolor="$lightColor" $colspanStrg><TABLE border=0 cellspacing=0 cellpadding=0>
	<TR><TD rowspan=2 nowrap><FONT class="title3">&bull;View:</FONT><SELECT name="view" class="title3" onchange="updateFiltersClass(this.value); /*document.displayForm.submit()*/"><OPTION value="graph">Graphical</OPTION><OPTION value="list"$selListView>List</OPTION></SELECT>&nbsp;&nbsp;&nbsp;
|;
	if ($ratioType eq 'None') {
		print "</TD><TD colspan=2 nowrap><FONT style=\"font-weight:bold\">&bull;Measure:</FONT><SELECT name=\"dispMeasure\" style=\"font-weight:bold\" onchange=\"alterUpdateButton(); /*document.displayForm.submit()*/\">";
		my $dbhLite=&promsQuantif::dbConnectProteinQuantification($selQuantifID,$projectID);
		my $sthMeas=$dbhLite->prepare("SELECT 1 FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIF_PARAMETER=? LIMIT 1");
		foreach my $refMeas (@{$proteinQuantifFamilies{'MEASURES'}{$quantifType}}) {
			my ($measCode,$measName,$isOptional)=@{$refMeas};
			my $disabStrg='';
			unless ($abundanceMeasures{$measCode}) {
				if ($isOptional) { # check if exist data
					$sthMeas->execute($quantifParamInfo{$measCode}[0]);
					my ($OK)=$sthMeas->fetchrow_array;
					# $disabStrg='disabled' unless $OK; # flag for disabled
					next unless $OK;
				}
				else {
					$disabStrg='disabled';
				}
			}
			my $selStrg=($measCode eq $dispMeasure && !$disabStrg)? 'selected' : '';
			print "<OPTION value=\"$measCode\" $disabStrg $selStrg>$measName</OPTION>";
		}
		$sthMeas->finish;
		$dbhLite->disconnect;
		print "</SELECT>\n";
	}
	else {
		my ($dispFullRatioStrg,$dispAllowInf)=($view eq 'graph')? ('','style="display:none"') : ('style="display:none"','');
	 	my $chkFullStrg=($fullRatioRange)? 'checked' : '';
		my $chkInfStrg=($allowInfRatio)? 'checked' : '';
		print qq
|<SPAN id="fullRatioSPAN" $dispFullRatioStrg><BR><B>&bull;<LABEL><INPUT type="checkbox" name="fullRatioRange" value=1 $chkFullStrg>Full ratio range</LABEL></B>&nbsp;&nbsp;&nbsp;</SPAN></TD>
<TD colspan=2 nowrap><FONT id="foldChgRawFONT1" class="$foldChPvalClass">&bull;|;
		if ($quantifSoftware eq 'MaxQuant') {
			my $selRawStrg=($foldChangeRaw)? ' selected' : '';
			print "</FONT><SELECT name=\"foldChgRaw\" id=\"foldChgRaw\" class=\"$foldChPvalClass\" onchange=\"alterUpdateButton()\"><OPTION value=\"0\">Norm.</OPTION><OPTION value=\"1\"$selRawStrg>Raw</OPTION></SELECT><FONT id=\"foldChgRawFONT2\" class=\"$foldChPvalClass\" onchange=\"alterUpdateButton()\">";
		}
		print qq
|Ratio:</FONT><SELECT name="foldChgType" id="foldChgType" class="$foldChPvalClass" onchange="alterUpdateButton()"><OPTION value="abs">Up &ge; or Down &le; 1/</OPTION><OPTION value="up"$selFcUp>Up &ge;</OPTION><OPTION value="down"$selFcDown>Down &le; 1/</OPTION></SELECT><INPUT type="text" name="foldChange" id="foldChange" value="$dispFoldChange" size=2 class="$foldChPvalClass" onchange="alterUpdateButton()"/>
|;
		if ($quantifSoftware eq 'MaxQuant') {print "<INPUT type=\"hidden\" name=\"pValue\" id=\"pValue\" value=\"1\" onchange=\"alterUpdateButton()\"/>";}
		else { # p-value
			print qq
|&nbsp;&nbsp;&nbsp;<FONT id="pValueFONT" class="$foldChPvalClass">&bull;$pValueType &le;</FONT><INPUT type="text" name="pValue" id="pValue" value="$dispPvalue" size=5 class="$foldChPvalClass" onchange="alterUpdateButton()"/>
<SPAN id="allowInfSPAN" $dispAllowInf><FONT class="trueFilter">(</FONT><LABEL class="trueFilter"><INPUT type="checkbox" name="allowInfRatio" value=1 onclick="alterUpdateButton()" $chkInfStrg>Allow infinite ratios<SUP onmouseover="popup('<B>p-value</B> threshold &lt;1 will not apply to <B>infinite ratios</B>')" onmouseout="popout()">?</SUP>)</LABEL></SPAN>
|;
		}
		if ($quantifType eq 'PROT_RATIO_PEP') {
			if ($quantifSoftware eq 'MaxQuant') {print "<INPUT type=\"hidden\" name=\"coefVar\" id=\"coefVar\" value=\"0\"/>";}
			else {print qq |<SPAN style="display:none">&nbsp;&nbsp;&nbsp;<FONT id="coefVarFONT1" class="trueFilter">&bull;Coeff. var.<SUP onmouseover="popup('Use <B>0%</B> for no filtering')" onmouseout="popout()">?</SUP> &le;</FONT><INPUT type="text" name="coefVar" id="coefVar" value="$dispCV" size=2 class="trueFilter" onchange="alterUpdateButton()"/><FONT id="coefVarFONT2" class="trueFilter">%</FONT></SPAN>|;} # hidden: not useful (PP 10/12/20)
		}
		#else {} # TNPQ
	}
	if (!$protRulerCase) {  # Peptide selection -> not relevant for Proteomic Ruler (VL 30/07/19)
		# Peptide selection
		print qq |&nbsp;&nbsp;&nbsp;<FONT class="trueFilter">&bull;</FONT><SELECT name="pepType" class="trueFilter" onchange="updateReplicStatus(this.value); alterUpdateButton()">|; # onchange="document.displayForm.submit()"
		if ($quantifSoftware eq 'MaxQuant') {
			my $ratioPepOptStrg=($ratioType eq 'None')? '' : "<OPTION value=\"NUM_PEP_USED\"$selRatioCountType>Ratio counts</OPTION>";
			print qq
|<OPTION value="PEPTIDES">All peptides</OPTION><OPTION value="RAZ_UNI_PEP"$selRazPepType>Razor+unique peptides</OPTION><OPTION value="UNIQUE_PEP"$selUniPepType>Unique peptides</OPTION>$ratioPepOptStrg
</SELECT> <FONT class="trueFilter">&ge;</FONT>
|;
		}
		else {
			print qq |<OPTION value="all">All</OPTION><OPTION value="distinct"$selDistPepType>Dist.</OPTION>|;
			print qq |<OPTION value="msmsNum"$selNumTruePepType># Valid.</OPTION><OPTION value="msmsPc"$selPcTruePepType>% Valid.</OPTION>| if $hasTruePepCount;
			#if ($quantifSoftware eq 'myProMS' && $algoVersion >= 2) {
			#	print qq |<OPTION value="allRep"$selAllRepPepType>All/replic.</OPTION><OPTION value="distinctRep"$selDistRepPepType>Dist./replic.</OPTION>|;
			#}
			print "</SELECT>";
			print "<SUP class=\"trueFilter\" onmouseover=\"popup('<B>Validated:</B> peptides formally identified by <B>MSMS</B><BR>(Unlike those identified by <B>MBWR</B>)')\" onmouseout=\"popout()\">?</SUP>" if $hasTruePepCount;
			print " <FONT class=\"trueFilter\">peptides &ge;</FONT>\n";
		}
		my $numPepWidth=($hasTruePepCount)? '50px' : '50px';
		print qq
|<INPUT type="number" name="numPep" min="1" max="100" value="$dispNumPep" style="width:$numPepWidth" class="trueFilter" onchange="alterUpdateButton()"/>|;
		if ($quantifSoftware eq 'myProMS' && $algoVersion >= 3 && $quantifType ne 'PROT_ABUNDANCE') {
			my $maxNumReplicInState=1;
			foreach my $statePos (keys %stateInfo) {
				$maxNumReplicInState=$stateInfo{$numStates}{'NUM_REPLIC'} if $maxNumReplicInState < $stateInfo{$numStates}{'NUM_REPLIC'};
			}
if ($hasTruePepCount|| $maxNumReplicInState > 1) {
			print qq
|<FONT class="trueFilter"> for </FONT>
|;
			if ($maxNumReplicInState > 1) {
				my $dispRepSpan=($numPepScope=~/^replicate/)? '' : 'none';
				my $replicHelpStrg=qq
|&bull;If <B>replic. in both states</B> only 1 state is considered for 1/&infin; and &infin; ratios.<BR>&bull;Max. number of replic./state is used if filter is set higher.|;
				print qq
|<SPAN id="pepReplicSPAN" style="display:$dispRepSpan"><FONT class="trueFilter"> &ge;</FONT>
<INPUT type="number" name="numReplic" min=1 max=$maxNumReplicInState value="$dispNumReplic" style="width:40px" class="trueFilter" onchange="alterUpdateButton()"/><SUP class="trueFilter" onmouseover="popup('$replicHelpStrg')" onmouseout="popout()">?</SUP></SPAN>
|;
			}
			print qq
|<SELECT name="numPepScope" class="trueFilter" onchange="alterUpdateButton(); document.getElementById('pepReplicSPAN').style.display=(this.value.match('replicate'))? '' : 'none'">
<OPTION value="ratio">ratio</OPTION>
|;		
				my @filterScope=(['One','one state'],['Both','both states'],['Test','test state'],['Ref','reference state']);
			if ($maxNumReplicInState > 1) {
					#print "<OPTGROUP value=\""
					my $disabReplic=($dispPepType =~ /^msms/)? ' disabled' : '';
					foreach my $refScope (@filterScope) {
						my ($scope,$scopeStrg)=@{$refScope};
						print "<OPTION value=\"replicate$scope\"";
						print ' selected' if $numPepScope eq 'replicate'.$scope;
						print "$disabReplic>replic. in $scopeStrg</OPTION>";
					}

			}
			if ($hasTruePepCount) {
					my $disabMsms=($dispPepType =~ /^msms/)? '' : ' disabled';
					foreach my $refScope (@filterScope) {
						my ($scope,$scopeStrg)=@{$refScope};
						print "<OPTION value=\"msms$scope\"";
						print ' selected' if $numPepScope eq 'msms'.$scope;
						print "$disabMsms>$scopeStrg</OPTION>";
					}
			}
			print "</SELECT>\n";
}
		}
	}
	my $selExclStrg=($listAction eq 'exclude')? 'selected' : '';
	print qq
|&nbsp;&nbsp;&nbsp;</TD>
<TD rowspan=2><INPUT type="submit" id="submitBUT" value="Update" class="noFilter"/></TD></TR>
	<TR><TD nowrap width=10%><FONT class="trueFilter">&bull;</FONT><SELECT name="listAction" class="trueFilter" onchange="alterUpdateButton()"><OPTION value="restrict">Restrict to</OPTION><OPTION value="exclude" $selExclStrg>Exclude</OPTION></SELECT><FONT class="trueFilter">&nbsp;proteins in List:&nbsp;</FONT></TD><TD><SPAN id="restrictSPAN"><!--Custom list selection comes here with window.onload --></SPAN>
	&nbsp;&nbsp;&nbsp;<B>&bull;Sort by:</B><SELECT name="sort" id="sort" onchange="updateSortBy()" style="font-weight:bold">
|;
	my @sortOptions=(['identifier','Identifiers']);
	if ($ratioType eq 'None') {
		if (scalar keys %stateInfo == 1) {
			if ($protRulerCase && $refInfo->{'AVG_MODE'}[0] eq "AVG_ALL"){
				push @sortOptions,['state_0','Measure'];
			} else {
				push @sortOptions,['state_1','Measure'];
			}
		}
		else {
			foreach my $statePos (sort{$a<=>$b} keys %stateInfo) {
				push @sortOptions,["state_$statePos","Measure: '$stateInfo{$statePos}{NAME}'"];
			}
		}
	}
	else {
		if (scalar @{$refInfo->{RATIOS}} == 1) {
			push @sortOptions,['ratio_1','Ratio'];
		}
		else {
			my $ratioPos=0;
			foreach my $ratioData (@{$refInfo->{RATIOS}}) {
				$ratioPos++;
				next unless $dispRatios{$ratioPos};
				my ($testStatePos,$refStatePos)=split(/\//,$ratioData);
				#if ($parentQuantifType && $parentQuantifType eq 'XIC'){#}
				my $ratioTag='';
				if ($designID) {
					if ($testStatePos=~/%/) {
						$ratioTag=$encodedDegStrg;
						$testStatePos=~s/%\d+//;
						$refStatePos=~s/%\d+//;
					}
					$testStatePos=$condToState{$testStatePos};
					$refStatePos=$condToState{$refStatePos};
				}
				push @sortOptions,["ratio_$ratioPos","Ratio: '$stateInfo{$testStatePos}{NAME}$ratioTag/$stateInfo{$refStatePos}{NAME}$ratioTag'"];
			}
		}
	}
	if ($protRulerCase) {
		push @sortOptions,(['mw','Molecular weight']);
		$dispSort=(!$dispSort || $dispSort eq 'peptide')? (($refInfo->{'AVG_MODE'}[0] eq "AVG_ALL")? "state_0" : "state_1") : $dispSort;
	} else {
		push @sortOptions,(['peptide','Peptide filter'],['mw','Molecular weight']);
	}	
	foreach my $refOption (@sortOptions) {
		print "<OPTION value=\"$refOption->[0]\"";
		print ' selected' if $refOption->[0] eq $dispSort;
		print ">$refOption->[1]</OPTION>";
	}
	print qq
|</SELECT>&nbsp;&nbsp;&nbsp;</TD>
</TR>
</TABLE></TD></TR>
|;
##	if (scalar keys %{$refGoAnalyses}) {
##		print qq
##|<TR><TH align=right nowrap valign=top>&nbsp;<FONT class="title3">GO decoration :</FONT></TH><TD bgcolor=$lightColor$colspanStrg><TABLE border=0 cellspacing=0 cellpadding=0>
##	<TR><TD nowrap><SELECT id="goAna" class="title3" onchange="updateGoAspects(this.value)"><OPTION value="0">-= None =-</OPTION>
##|;
##		foreach my $goID (sort{lc($refGoAnalyses->{$a}[0]) cmp lc($refGoAnalyses->{$b}[0])} keys %{$refGoAnalyses}) {
##			print "<OPTION value=\"$goID\">$refGoAnalyses->{$goID}[0]</OPTION>";
##		}
##		print qq
##|</SELECT><FONT class="title3">:&nbsp;</FONT><SELECT id="goAspect" class="title3" style="visibility:hidden" onchange="ajaxUpdateGoTermList(this.value)"><OPTION value="">-= Select Aspect =-</OPTION></SELECT>&nbsp;</TD>
##<TD><DIV id="goTermsDIV"></DIV></TD></TR></TD></TABLE></TD></TR>
##|;
##	}

	if ($view eq 'graph' && $dispResults) {
		print qq |<TR><TH align=right nowrap valign=top>&nbsp;<FONT class="title3">Chart highlighting :</FONT></TH><TD bgcolor=$lightColor$colspanStrg><TABLE border=0 cellspacing=0 cellpadding=0>
	<TR><TH align=right nowrap>&nbsp;Custom lists:</TH><TD nowrap colspan=3><SPAN id="listDecoSPAN"></SPAN></TD></TR>
|;
		if (scalar keys %{$refGoAnalyses}) {
			print qq|	<TR><TH align=right nowrap>&nbsp;Gene Ontology:</TH><TD nowrap><SELECT id="goAna" onchange="updateGoAspects(this.value)"><OPTION value="0">-= None =-</OPTION>
|;
			foreach my $goID (sort{lc($refGoAnalyses->{$a}[0]) cmp lc($refGoAnalyses->{$b}[0])} keys %{$refGoAnalyses}) {
				print "<OPTION value=\"$goID\">$refGoAnalyses->{$goID}[0]</OPTION>";
			}
			print qq |</SELECT><FONT class="title3">:&nbsp;</FONT><SELECT id="goAspect" style="visibility:hidden" onchange="ajaxUpdateGoTermList(this.value)"><OPTION value="">-= Select Aspect =-</OPTION></SELECT>&nbsp;</TD>
<TD><DIV id="goTermsDIV"></DIV></TD></TR>
|;
		}
		print "</TABLE></TD></TR>\n";
	}

	print qq
|</TABLE>
</FORM>
<BR>
|;
	unless ($dispResults) { # form was not submitted (1st time displayed)
		$dbh->disconnect;
		print qq
|<SCRIPT type="text/javascript">
window.onload=function() {
	ajaxUpdateRestrict();
}
</SCRIPT>
<FONT class="title1"><BR><BR><BR>Choose settings and click on "Update" to display data</FONT>
|;
		&endHTML;
		exit;
	}
	return %stateInfo;
} # end of &printQuantificationSummary



sub sortProteins {
	my ($sortItem,$refProteinInfo,$refQuantifValues,$firstRatioUsed)=@_;
	my ($pID1,$modStrg1)=($a=~/^(\d+)-?(.*)/); $modStrg1='' unless $modStrg1;
	my ($pID2,$modStrg2)=($b=~/^(\d+)-?(.*)/); $modStrg2='' unless $modStrg2;
	if ($sortItem eq 'identifier') {lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))}
	elsif ($sortItem eq 'mw') {lc($refProteinInfo->{$pID2}[1]) <=> lc($refProteinInfo->{$pID1}[1]) || lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))}
	elsif ($sortItem =~ /ratio_(-?\d+)/) {
		my $ratioPos=$1;
		if ($refQuantifValues->{$a} && $refQuantifValues->{$b}) {
			#$refQuantifValues->{$b}{'RATIO'}{$ratioPos}<=>$refQuantifValues->{$a}{'RATIO'}{$ratioPos} || $refQuantifValues->{$b}{'NUM_PEP_USED'}{$ratioPos}<=>$refQuantifValues->{$a}{'NUM_PEP_USED'}{$ratioPos} || lc($refProteinInfo->{$a}[0]) cmp lc($refProteinInfo->{$b}[0])
			&sortCheckDefined($ratioParamCode,$refQuantifValues,$ratioPos,-1) || &sortCheckDefined('NUM_PEP_USED',$refQuantifValues,$ratioPos,-1) || &sortCheckDefined('NUM_PEP_TOTAL',$refQuantifValues,$ratioPos,-1) || lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))
		}
		elsif ($refQuantifValues->{$a}) {-1}
		elsif ($refQuantifValues->{$b}) {1}
		else {lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))}
	}
	elsif ($sortItem =~ /state_(\d+)/) {
		my $statePos=$1;
		if ($refQuantifValues->{$a} && $refQuantifValues->{$b}) {
			#$refQuantifValues->{$b}{'RATIO'}{$ratioPos}<=>$refQuantifValues->{$a}{'RATIO'}{$ratioPos} || $refQuantifValues->{$b}{'NUM_PEP_USED'}{$ratioPos}<=>$refQuantifValues->{$a}{'NUM_PEP_USED'}{$ratioPos} || lc($refProteinInfo->{$a}[0]) cmp lc($refProteinInfo->{$b}[0])
			&sortCheckDefined($dispMeasure,$refQuantifValues,$statePos,-1) || &sortCheckDefined('PEPTIDES',$refQuantifValues,$statePos,-1) || &sortCheckDefined('RAZ_UNI_PEP',$refQuantifValues,$statePos,-1) || &sortCheckDefined('UNIQUE_PEP',$refQuantifValues,$statePos,-1) || lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))
		}
		elsif ($refQuantifValues->{$a}) {-1}
		elsif ($refQuantifValues->{$b}) {1}
		else {lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))}
	}
	# defaults to peptide count (NUM_PEP_USED / DIST_PEP_USED, ...)
	else {
		if ($refQuantifValues->{$a} && $refQuantifValues->{$b}) { # $sortItem || NUM_PEP_USED because they can be different
			#$refQuantifValues->{$b}{$numPepCode}{1}<=>$refQuantifValues->{$a}{$numPepCode}{1} || $refQuantifValues->{$b}{'RATIO'}{1}<=>$refQuantifValues->{$a}{'RATIO'}{1} || lc($refProteinInfo->{$a}[0]) cmp lc($refProteinInfo->{$b}[0])
			&sortCheckDefined($sortItem,$refQuantifValues,$firstRatioUsed,-1) || &sortCheckDefined('NUM_PEP_USED',$refQuantifValues,$firstRatioUsed,-1) || &sortCheckDefined($ratioParamCode,$refQuantifValues,$firstRatioUsed,-1) || lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))
		}
		elsif ($refQuantifValues->{$a}) {-1}
		elsif ($refQuantifValues->{$b}) {1}
		else {lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))}
	}
}

sub sortCheckDefined {
	my ($qCode,$refQuantifValues,$ratioPos,$order)=@_;
	$order=1 if !$order;
	if (defined $refQuantifValues->{$a}{$qCode}{$ratioPos} && defined $refQuantifValues->{$b}{$qCode}{$ratioPos}) {
		if ($order==1) {$refQuantifValues->{$a}{$qCode}{$ratioPos}<=>$refQuantifValues->{$b}{$qCode}{$ratioPos}}
		else {$refQuantifValues->{$b}{$qCode}{$ratioPos}<=>$refQuantifValues->{$a}{$qCode}{$ratioPos}}
	}
	elsif (defined $refQuantifValues->{$a}{$qCode}{$ratioPos}) {1*$order}
	elsif (defined $refQuantifValues->{$b}{$qCode}{$ratioPos}) {-1*$order}
	else {0}
}

sub displayVolcanoPlot {
	my ($refLabelingInfo,$refStateInfo,$refQuantifValues,$refProteinInfo,$refDispModifSites,$pvalueCode,$minPvalue)=@_;
	my $pValueType=(!$refLabelingInfo->{'FDR_CONTROL'} || $refLabelingInfo->{'FDR_CONTROL'}[0]=~/FALSE|none/i)? 'p-value' : 'Adj. p-value';
	print qq
|<SCRIPT type="text/javascript">
var mainChart,VP;
window.onload=function() {
	ajaxUpdateRestrict();
	VP=new volcanoPlot({div:'mainGraphDIV',
						width:500,height:450,
						foldChange:$dispFoldChange,
						pValue:$dispPvalue,
						pValueType:'$pValueType',
						allowHighlight:true,
						searchable:{text:'Extended search',externalSearch:ajaxSearchConvertIdentifier},
						updateHighlight:{callback:updateUsedColors,editable:false},
						pointOnclick:proteinLink,
						pointOnList:ajaxListSelectedProteins,
						exportAsImage:['Export as image','VolcanoPlot','./exportSVG.cgi']
						});
	mainChart=VP; // needed for highlighting
|;
	my %stateScopeLabel=(
		'Both'=>' used',
		'One'=>' (best st.)',
		'Test'=>' (test st.)',
		'Ref'=>' (ref. st.)'
	);
	# my $pepSizeStrg;
	# if ($numPepScope=~/^replicate(.+)/) {
	# 	my $scope=$1;
	# 	my $distStrg=($dispPepType eq 'all')? '' : ' dist.';
	# 	$pepSizeStrg="Min.$distStrg pep./rep.$stateScopeLabel{$scope}";
	# 	#$pepSizeStrg=($numPepScope eq 'replicateBoth')? "Min.$distStrg pep./rep." : "Min.$distStrg pep./rep. (best st.)";
	# }
	# elsif ($dispPepType =~ /^msms/) {
	# 	if ($numPepScope=~/^msms(.+)/) {
	# 		my $scope=$1;
	# 		$pepSizeStrg=($numPepScope eq 'msmsBoth')? 'Min. valid. ' : 'Valid. ';
	# 		$pepSizeStrg.="pep.$stateScopeLabel{$scope}";
	# 	}
	# 	else {$pepSizeStrg="Valid. pep. used"} # ratio
	# 	#$pepSizeStrg=($numPepScope eq 'ratio')? "Valid. pep. used" : ($numPepScope eq 'msmsBoth')? 'Min. valid. pep. both st.' : 'Min. valid. pep. best st.';
	# }
	# else {
	# 	$pepSizeStrg=($dispPepType eq 'all')? 'All pep. used' : 'Dist. pep. used';
	# }
	my $pepSizeStrg='All pep. used';
	my $ratioIdx=0;
	my (%ratioPos2index,%supRatios);
	foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
		my ($testStatePos,$refStatePos)=($ratioPos > 0)? split(/\//,$refLabelingInfo->{'RATIOS'}[$ratioPos-1]) : reverse(split(/\//,$refLabelingInfo->{'RATIOS'}[abs($ratioPos)-1]));
		my $ratioTag='';
		if ($designID) {
			if ($testStatePos=~/%/) {
				$ratioTag=$encodedDegStrg;
				$testStatePos=~s/%\d+//;
				$refStatePos=~s/%\d+//;
				$supRatios{$ratioPos}=1;
			}
			$testStatePos=$condToState{$testStatePos};
			$refStatePos=$condToState{$refStatePos};
		}
		print "\tVP.addDataSet($ratioIdx,{name:'$refStateInfo->{$testStatePos}{NAME}$ratioTag/$refStateInfo->{$refStatePos}{NAME}$ratioTag',sizeName:'$pepSizeStrg'});\n";
		$ratioPos2index{$ratioPos}=$ratioIdx;
		$ratioIdx++;
	}
	my (%dispProteins,%dispIsoforms);
	my ($pepCode4Filter,$pepFilterValue)=($dispPepType eq 'msmsPc')? ('FRAC_TRUE_USED',$dispNumPep/100) : ($numPepScope =~ /^replicate/ || $dispPepType =~ /^msms/)? ($numPepCode.'_SUBSET',$dispNumPep) : ($numPepCode,$dispNumPep); # numPepCode is GLOBAL
	my $pepCode4Size=($quantifSoftware eq 'myProMS' || $numPepCode eq 'DIST_PEP_USED')? 'NUM_PEP_USED' : $numPepCode; # Always show all pep.
	foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
		my $count=0;
		my $ratioIdx=$ratioPos2index{$ratioPos};
		foreach my $modProtID (keys %{$refQuantifValues}) {
			my ($protID,$modStrg)=($modProtID=~/^(\d+)(.*)/); $modStrg='' unless $modStrg;
			#next if (!$quantifValues{$modProtID}{$pvalueCode} || !$quantifValues{$modProtID}{$pvalueCode}{$ratioPos} || !$quantifValues{$modProtID}{'RATIO'} || !$quantifValues{$modProtID}{'RATIO'}{$ratioPos});
			next if (!$refQuantifValues->{$modProtID}{'RATIO'} || !$refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos});
			
			next if $refQuantifValues->{$modProtID}{$pepCode4Filter}{$ratioPos} < $pepFilterValue;

			#next if ($dispConfInt && defined $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos} && 1.96*$refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos} > $dispConfInt);
			if ($dispCV && defined $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}) {
				if ($ratioType eq 'Ratio') {
					next if abs(log($refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos})/log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos})) > $dispCVFrac;
				}
				else {
					next if $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}/abs(log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos})/$log2) > $dispCVFrac;
				}
			}
			if ($count==0) {
				print "\tVP.addDataAsString($ratioIdx,'";
			}
			else {print ";";}
			$count++;
#$modStrg=~s/,/./g; # for compatibility with volcano-plot
#(my $usedProtID=$modProtID)=~s/,/./g; # for compatibility with volcano-plot
			my $dispModIdent=($refDispModifSites->{$modProtID})? $refProteinInfo->{$protID}[0].'-'.$refDispModifSites->{$modProtID} : $refProteinInfo->{$protID}[0].$modStrg;
			print "$dispModIdent,$modProtID,";
			#>Set missing p-values=1!!!
			if (!$refQuantifValues->{$modProtID}{$pvalueCode} || !defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}) {
				$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}=1;
			}
			elsif ($refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}==0) {
				$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}=$minPvalue;
			}
			#>Set +/- infinite ratios
			if ($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos} == 1000 || $refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos} > $MAX_INF_RATIO) {
				$refProteinInfo->{$protID}[4]=1000 unless $refProteinInfo->{$protID}[4]; # set prot length
				printf "+,%.2f",100*$refQuantifValues->{$modProtID}{$pepCode4Size}{$ratioPos}/$refProteinInfo->{$protID}[4]; # num pep/100 res
			}
			elsif ($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos} == 0.001 || $refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos} < $MIN_INF_RATIO) {
				$refProteinInfo->{$protID}[4]=1000 unless $refProteinInfo->{$protID}[4]; # set prot length
				printf "-,%.2f",100*$refQuantifValues->{$modProtID}{$pepCode4Size}{$ratioPos}/$refProteinInfo->{$protID}[4]; # num pep/100 res
			}
			else { # normal ratios within [$MIN_INF_RATIO-$MAX_INF_RATIO] boundaries
				#$quantifValues{$modProtID}{$pvalueCode}{$ratioPos}=1 if (!$quantifValues{$modProtID}{$pvalueCode} || !$quantifValues{$modProtID}{$pvalueCode}{$ratioPos});
				#if ($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos} >= $MAX_INF_RATIO) {$refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos}='+';}
				#elsif ($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos} <= $MIN_INF_RATIO) {$refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos}='-';}
				printf "$refQuantifValues->{$modProtID}{RATIO}{$ratioPos},%.2e",$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos};
			}
			print ",$refQuantifValues->{$modProtID}{$pepCode4Size}{$ratioPos}";
			if ($count==100) {
				print "');\n";
				$count=0;
			}
			$dispProteins{$protID}=1;
			$dispIsoforms{$modProtID}=1 if $modStrg;
		}
		if ($count > 0) {
			print "');\n";
		}
	}
	my $numProtDispStrg=(scalar keys %dispProteins).' proteins';
	$numProtDispStrg.=' / '.(scalar keys %dispIsoforms).' sites' if scalar keys %dispIsoforms;
	$numProtDispStrg.=' displayed';
	print qq
|	VP.draw();
}
</SCRIPT>
</CENTER>
<TABLE cellspacing=0 cellpadding=0>
<TR><TD><FONT class="title3">$numProtDispStrg</FONT></TD></TR>
<TR><TD><DIV id="mainGraphDIV" style="float:left;background-color:#FFFFFF;border:solid 3px $darkColor;padding:2px;"></DIV></TD></TR>
</TABLE>
<DIV id="protListDIV" style="position:absolute;height:535;overflow:auto"></DIV>
|;
}

sub displayIntensityPlot {
	my ($refLabelingInfo,$refStateInfo,$refQuantifValues,$refProteinInfo,$refDispModifSites)=@_; # globals: $dispMeasure,$numPepCode,%quantifParamInfo
	my $measureName=$quantifParamInfo{$dispMeasure}[1];
	my $pepCode4Display=($quantifType eq 'PROT_ABUNDANCE')? 'NUM_PEP_USED' : ($quantifType eq 'MQ')? 'PEPTIDES' : $numPepCode;
	my $numPepName=$quantifParamInfo{$pepCode4Display}[1];
	$numPepName.=' pep.' if ($quantifType eq 'MQ' && $pepCode4Display ne 'PEPTIDES');
	my $isLogged=1;
	my $axisYtitle=($isLogged)? (1,"Log10($measureName)") : (0,$measureName);
	my $log10=log(10);
#print '**',join(':',sort{$a<=>$b} keys %dispStates),"**<BR>\n";
	print qq
|</CENTER>
<SCRIPT type="text/javascript">
function ipPointLabel(dp,type) {
	var infoStrg=dp.label;
	if (type != 'min') {
		if (dp.dataSet.params.chart.dataSets.length > 1) {infoStrg+='\\nState: '+dp.dataSet.params.name;}
		var intensity=($isLogged==1)? (1*(Math.pow(10,dp.getY())).toFixed(1)) : dp.getY();
		infoStrg+='\\n$measureName='+intensity+'\\nrank='+dp.getX()+'\\n$numPepName='+dp.size;
	}
	return infoStrg;
}
var mainChart,IP;
window.onload=function() {
	ajaxUpdateRestrict();
	IP=new genericPlot({div:'mainGraphDIV',width:600,height:450,
						axisX:{title:'Rank'}, // ,zoomable:true
						axisY:{title:'$axisYtitle'}, // ,zoomable:true
						axisClosure:true,
						//hideDataSets:true, // hide dataset label if only 1 displayed!
						zoomable:true,
						pointLabel:ipPointLabel,
						pointOnclick:proteinLink,
						pointOnList:ajaxListSelectedProteins,
						allowHighlight:true,
						searchable:{text:'Extended search',externalSearch:ajaxSearchConvertIdentifier},
						updateHighlight:{callback:updateUsedColors,editable:false},
						exportAsImage:['Export as image','IntensityPlot','./exportSVG.cgi']
					});
	mainChart=IP; // needed for highlighting
|;
	my $dsetIdx=0;
	foreach my $statePos (sort{$a<=>$b} keys %dispStates) {
		print "\tIP.addDataSet($dsetIdx,{name:'$refStateInfo->{$statePos}{NAME}'});\n";
		$dsetIdx++;
	}
	my $protCount;
	$dsetIdx=0;
	my (%dispProteins, %dispIsoforms);
	my ($pepCode4Filter,$pepFilterValue)=($dispPepType eq 'msmsPc')? ('FRAC_TRUE_USED',$dispNumPep/100) : ($numPepCode,$dispNumPep); # numPepCode is GLOBAL
	foreach my $statePos (sort{$a<=>$b} keys %dispStates) {
		my $count=0;
		$protCount=0;
		foreach my $modProtID (sort{&sortCheckDefined($dispMeasure,$refQuantifValues,$statePos,1) || &sortCheckDefined($pepCode4Display,$refQuantifValues,$statePos,1)} keys %{$refQuantifValues}) {
			next if (!$refQuantifValues->{$modProtID}{$dispMeasure} || !$refQuantifValues->{$modProtID}{$dispMeasure}{$statePos});
			next if (!$refQuantifValues->{$modProtID}{$pepCode4Filter} || !$refQuantifValues->{$modProtID}{$pepCode4Filter}{$statePos} || $refQuantifValues->{$modProtID}{$pepCode4Filter}{$statePos} < $pepFilterValue); # Not in DB if ==0 (eg. UNIQUE_PEP)
			#next if ($dispConfInt && defined $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos} && 1.96*$refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos} > $dispConfInt);
			if ($count==0) {
				print "\tIP.addDataAsString($dsetIdx,'";
			}
			else {print ";";}
			$count++;
			$protCount++;
			my ($protID,$modStrg)=($modProtID=~/^(\d+)(.*)/); $modStrg='' unless $modStrg;
			$refProteinInfo->{$protID}[0]=~s/[,;']/\./g; # Clean MaxQuant crappy contaminant identifiers
			my $dispModIdent=($refDispModifSites->{$modProtID})? $refProteinInfo->{$protID}[0].'-'.$refDispModifSites->{$modProtID} : $refProteinInfo->{$protID}[0].$modStrg;
			print "$dispModIdent,$modProtID,$protCount,";
			if ($isLogged) {print log($refQuantifValues->{$modProtID}{$dispMeasure}{$statePos})/$log10;} else {print $refQuantifValues->{$modProtID}{$dispMeasure}{$statePos};}
			print ',',$refQuantifValues->{$modProtID}{$pepCode4Display}{$statePos};
			if ($count==100) {
				print "');\n";
				$count=0;
			}
			$dispProteins{$protID}=1;
			$dispIsoforms{$modProtID}=1 if $modStrg;
		}
		if ($count > 0) {
			print "');\n";
		}
		$dsetIdx++;
	}
	my $numProtDispStrg=(scalar keys %dispProteins).' proteins';
	$numProtDispStrg.=' / '.(scalar keys %dispIsoforms).' sites' if scalar keys %dispIsoforms;
	$numProtDispStrg.=' displayed';
	print qq
|	IP.draw();
}
</SCRIPT>
<TABLE cellspacing=0 cellpadding=0>
<TR><TD><FONT class="title3">$numProtDispStrg</FONT></TR></TD>
<TR><TD><DIV id="mainGraphDIV" style="float:left;background-color:#FFFFFF;border:solid 3px $darkColor;padding:2px;"></DIV></TR></TD>
</TABLE>
<DIV id="protListDIV" style="position:absolute;height:535;overflow:auto"></DIV>
|;
}

sub displayRatioPlot {
	my ($refLabelingInfo,$refStateInfo,$refQuantifValues,$refProteinInfo)=@_; # globals: $dispMeasure,$numPepCode,%quantifParamInfo
	my $numPepName=$quantifParamInfo{$numPepCode}[1];
	$numPepName.=' peptides' if $numPepCode ne 'PEPTIDES';
	my $revFoldChange=1/$dispFoldChange;
	my $log2=log(2);
	print qq
|</CENTER>
<SCRIPT type="text/javascript">
function myPointLabel(dp,type) {
	var infoStrg=dp.label;
	if (type != 'min') {
		if (dp.dataSet.params.chart.dataSets.length > 1) {infoStrg+='\\nState: '+dp.dataSet.params.name;}
		var ratio=(dp.getX() >= 0)? 1*((Math.pow(2,dp.getX())).toFixed(2)) : '1/'+(1*((1/(Math.pow(2,dp.getX()))).toFixed(2)));
		infoStrg+='\\nFC='+ratio+'\\nCoef. var.='+dp.getY()+'%\\n$numPepName='+dp.size;
	}
	return infoStrg;
}
function convertThreshold(axis,value) {
	if (axis=='X') {value=Math.log2(value);}
	return value;
}
var mainChart,RP;
window.onload=function() {
	ajaxUpdateRestrict();
	RP=new genericPlot({div:'mainGraphDIV',width:600,height:450,
						axisX:{title:'Log2(fold change)'}, // ,zoomable:true
						axisY:{title:'Coefficient variation (%)'}, // ,zoomable:true
						axisClosure:true,
						//hideDataSets:true,
						zoomable:true,
						editThreshold:true,
						convertValue:convertThreshold, // threshold lines only
						pointLabel:myPointLabel,
						pointOnclick:proteinLink,
						pointOnList:ajaxListSelectedProteins,
						allowHighlight:true,
						searchable:{text:'Extended search',externalSearch:ajaxSearchConvertIdentifier},
						updateHighlight:{callback:updateUsedColors,editable:false},
						exportAsImage:['Export as image','RatioPlot','./exportSVG.cgi']
					});
	mainChart=RP; // needed for highlighting
	RP.addThreshold({axis:'Y',label:'CV thresold',value:$dispCV,color:'#FF0000',editable:true});
	RP.addThreshold({axis:'X',label:'min. fold change',value:$revFoldChange,color:'#00A000',keepAbove1:true,editable:true});
	RP.addThreshold({axis:'X',label:'max. fold change',value:$dispFoldChange,color:'#00A000',keepAbove1:true,editable:true});
|;
	my $ratioIdx=0;
	my %ratioPos2index;
	foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
		my ($testStateID,$refStateID)=($ratioPos > 0)? split(/\//,$refLabelingInfo->{'RATIOS'}[$ratioPos-1]) : reverse(split(/\//,$refLabelingInfo->{'RATIOS'}[abs($ratioPos)-1]));
		my $testStatePos=$condToState{$testStateID};
		my $refStatePos=$condToState{$refStateID};
		print "\tRP.addDataSet($ratioIdx,{name:'$refStateInfo->{$testStatePos}{NAME}/$refStateInfo->{$refStatePos}{NAME}'});\n";
		$ratioPos2index{$ratioPos}=$ratioIdx;
		$ratioIdx++;
	}
	my %dispProteins;
	foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
		my $count=0;
		foreach my $modProtID (keys %{$refQuantifValues}) {
			my ($protID,$modStrg)=($modProtID=~/^(\d+)(.*)/); $modStrg='' unless $modStrg;
next if (!$refQuantifValues->{$modProtID}{$ratioParamCode} || !$refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} || !$refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos});
			next if $refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos} < $dispNumPep;
#next if ($dispCV && defined $refQuantifValues->{$modProtID}{'RATIO_VAR'}{$ratioPos} && $refQuantifValues->{$modProtID}{'RATIO_VAR'}{$ratioPos} > $dispCV);
			if ($count==0) {
				print "\tRP.addDataAsString($ratioPos2index{$ratioPos},'";
			}
			else {print ";";}
			$count++;
			$refProteinInfo->{$protID}[0]=~s/[,;']/\./g; # Clean MaxQuant crappy contaminant identifiers
			my $coefVar=$refQuantifValues->{$modProtID}{RATIO_VAR}{$ratioPos} || 0; # not defined for < 3 peptides
			print "$refProteinInfo->{$protID}[0]$modStrg,$modProtID,",log($refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos})/$log2,",$coefVar,$refQuantifValues->{$protID}{$numPepCode}{$ratioPos}";
			if ($count==100) {
				print "');\n";
				$count=0;
			}
			$dispProteins{$protID}=1;
		}
		if ($count > 0) {
			print "');\n";
		}
	}
	my $numProtDisp=scalar keys %dispProteins;
	print qq
|	RP.draw();
}
</SCRIPT>
<TABLE cellspacing=0 cellpadding=0>
<TR><TD><FONT class="title3">$numProtDisp proteins displayed</FONT></TD></TR>
<TR><TD><DIV id="mainGraphDIV" style="float:left;background-color:#FFFFFF;border:solid 3px $darkColor;padding:2px;"></DIV></TD></TR>
</TABLE>
<DIV id="protListDIV" style="position:absolute;height:535;overflow:auto"></DIV>
|;
}

sub listProteinRatios {
	my ($labelType,$ratioType,$refLabelingInfo,$refStateInfo,$refQuantifValues,$refProteinInfo,$refDispModifSites,$pvalueCode)=@_; # $labelType not global for all calls
	my ($minFoldChange,$maxFoldChange,$maxPvalue)=(0.5,2,0.05); # default values for arrow flag
	my $invDispFoldChange=1/$dispFoldChange;
	my $quantifSoftware=($refLabelingInfo->{'SOFTWARE'})? $refLabelingInfo->{'SOFTWARE'}[0] : 'myProMS'; $quantifSoftware='MaxQuant' if $quantifSoftware eq 'MQ';
	my $pValueType=(!$refLabelingInfo->{'FDR_CONTROL'} || $refLabelingInfo->{'FDR_CONTROL'}[0]=~/FALSE|none/i)? 'p-value' : 'Adj. p-value';
	my $numTotQuantifItems=scalar keys %{$refQuantifValues};
	my $numTotProteins=scalar keys %{$refProteinInfo};
	my $isModifQuantif=$refLabelingInfo->{'IS_MODIF'}[0];
	my $protString=($isModifQuantif)? "&nbsp;$numTotQuantifItems sites&nbsp;<BR>&nbsp;$numTotProteins proteins&nbsp;" : "&nbsp;$numTotProteins proteins&nbsp;";
	#my $peptideTitle=($numPepCode eq 'DIST_PEP_USED')? 'Dist. P' : 'P';
	my $peptideTitle=($dispPepType eq 'all')? 'Pep.' : ($dispPepType =~ /^msms/)? 'Valid. pep.' : 'Dist. pep.';
	$peptideTitle.=(($refLabelingInfo->{'ALGO_TYPE'} && $refLabelingInfo->{'ALGO_TYPE'}[0] eq 'PEP_RATIO') || ($refLabelingInfo->{'QUANTIFICATION_METHOD'} && $refLabelingInfo->{'QUANTIFICATION_METHOD'}[0] eq 'Ratio'))? ' sets' : '';
	my $disabSave=($projectAccess eq 'guest')? ' disabled' : '';
	my $clearButtonStrg=($view eq 'list')? '' : '<INPUT type="button" value="Clear list" style="font-weight:bold" onclick="document.getElementById(\'protListDIV\').style.display=\'none\'">';
	my ($ratioColspan,$cvColStrg)=($quantifType eq 'TNPQ')? (3,'') : (4,'<TH class="rbBorder" nowrap>&nbsp;CV (%)&nbsp;</TH>'); # ratio/arrow/p-value/(+/-stdDev)
my $pepTypeFiltering=0;
if ($dispPepType ne 'all' || $numPepScope ne 'ratio') { # includes Distinct peptides!!!!
	$ratioColspan++;
	$pepTypeFiltering=1;
}
#my $pepFilterTitle=($dispPepType =~ /^msms/)? 'Valid. pep.' : ($dispPepType eq 'all' && $numPepScope ne 'ratio')? 'Dist. pep.';
	my $selProtRequired = ($action eq 'ajaxListProt') ? 'true' : 'false';
	#my $ratioType=($refLabelingInfo->{'RATIO_TYPE'})? $refLabelingInfo->{'RATIO_TYPE'}[0] : 'Ratio';
	my %supRatios;
	if ($ratioType eq 'SuperRatio') {
		foreach my $ratioPos (keys %dispRatios) {
			$supRatios{$ratioPos}=1 if $refLabelingInfo->{'RATIOS'}[abs($ratioPos)-1]=~/%/;
		}
	}
	my $numColumns=scalar (keys %dispRatios) * $ratioColspan + 5;
	print qq
|<FORM name="protForm" method="POST">
<TABLE border=0 cellspacing=0 cellpadding=2>
<TR><TD colspan=$numColumns><DIV id="saveProtDIV" style="display:none;"></DIV></TD></TR>
<TR><TD colspan=$numColumns>$clearButtonStrg
<INPUT type="button" value="Check all" onclick="checkAllProteins(true)"/>
<INPUT type="button" value="Uncheck all" onclick="checkAllProteins(false)"/>
<INPUT type="button" id="saveFormBUTTON" value="Save proteins..." onclick="ajaxManageSaveProteins('getThemes','PROT')"$disabSave/>
|;
	print "<INPUT type=\"button\" id=\"saveSiteFormBUTTON\" value=\"Save sites...\" onclick=\"ajaxManageSaveProteins('getThemes','SITE')\"$disabSave/>\n" if $isModifQuantif;
	#print "<INPUT type=\"button\" value=\"Export data\" onclick=\"exportProteins()\"/>\n" if $action ne 'ajaxListProt';
	print "<INPUT type=\"button\" value=\"Export data\" onclick=\"checkProtToExport($selProtRequired)\"/>\n"; #if ($action ne 'ajaxListProt'); # && $algoVersion != 3); # !!!TEMP!!!
	#print "<INPUT type=\"button\" value=\"Export site fragments\" onclick=\"exportIsoforms()\"/>\n" if ($quantifModID && $userID eq 'sliva'); # !!!!! TEMP !!!!
	print qq
|</TD></TR>
<TR bgcolor="$darkColor">
<TH class="rbBorder" align=left rowspan=2><DIV id="protCountDIV">$protString</DIV></TH><TH class="rbBorder" rowspan=2>Gene</TH>
|; #&nbsp;$proteinTitle&nbsp;
	foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
		my ($testStatePos,$refStatePos)=($ratioPos > 0)? split(/\//,$refLabelingInfo->{'RATIOS'}[$ratioPos-1]) : reverse(split(/\//,$refLabelingInfo->{'RATIOS'}[abs($ratioPos)-1]));
		my $ratioTag='';
		if ($designID) {
			if ($supRatios{$ratioPos}) {
				$ratioTag=$encodedDegStrg;
				$testStatePos=~s/%\d+//;
				$refStatePos=~s/%\d+//;
			}
			$testStatePos=$condToState{$testStatePos};
			$refStatePos=$condToState{$refStatePos};
		}
		print "<TH class=\"rbBorder\" colspan=$ratioColspan nowrap>&nbsp;<A id=\"listRatio_$ratioPos\" href=\"javascript:displayRatioStat('listRatio_$ratioPos','$ratioCode{$ratioPos}','$refStateInfo->{$testStatePos}{NAME}$ratioTag/$refStateInfo->{$refStatePos}{NAME}$ratioTag')\">$refStateInfo->{$testStatePos}{NAME}$ratioTag/$refStateInfo->{$refStatePos}{NAME}$ratioTag</A>&nbsp;</TH>\n";
	}
	if ($ratioType eq 'Ratio') { # old algo
		my $popupInfo=($dispPepType eq 'all')? 'used/total' : 'distinct/used/total';
		print "<TH class=\"rbBorder\" rowspan=2 nowrap onmouseover=\"popup('<B>$popupInfo</B>')\" onmouseout=\"popout()\">&nbsp;$peptideTitle&nbsp;<BR>used</TH>\n";
	}
	else { # Super/Simple Ratio
		print "<TH class=\"rbBorder\" rowspan=2 nowrap>&nbsp;Peptides&nbsp<BR>in set</TH>\n";
	}
	print qq
|<TH class="rbBorder" rowspan=2 nowrap>&nbsp;MW<SMALL> kDa</SMALL>&nbsp;</TH><TH class="bBorder" rowspan=2 width=700 nowrap>&nbsp;Description - Species&nbsp;</TH>
</TR>
<TR bgcolor="$darkColor">
|;
	my $popupInfo='&bull;';
	if ($ratioType=~/S\w+Ratio/) {
		my %stateScopeLabel=(
			'Both'=>'in both states',
			'One'=>'in best state',
			'Test'=>'in test state',
			'Ref'=>'in reference state',
			'ratio'=>'in ratio'
		);
		if ($numPepScope=~/^replicate(.+)/) {
			my $scope=$1;
			my $distStrg=($dispPepType eq 'all')? 'pep.' : 'dist. pep.';
			$popupInfo.="From replic. with min. #$distStrg &ge; $dispNumPep $stateScopeLabel{$scope}";
			#$popupInfo=($numPepScope eq 'replicateBoth')? "Min.$distStrg pep. &ge; thres./replicate" : "Min.$distStrg pep. &ge; thres./replicate in best state";
		}
		elsif ($dispPepType =~ /^msms/) {
			my $scope='ratio'; # default
			if ($numPepScope=~/^msms(.+)/) {
				$scope=$1;
				$popupInfo.=($numPepScope eq 'msmsBoth' && $dispPepType eq 'msmsPc')? 'Min. % valid' : ($numPepScope eq 'msmsBoth')? 'Min. valid.' : ($dispPepType eq 'msmsPc')? '% valid.' : 'Valid.';
			}
			else {$popupInfo.='Valid.';} # ratio
			# $popupInfo.="+<FONT class=\\\'ghostPeptide\\\'>MBWR</FONT> pep. used $stateScopeLabel{$scope}";
			$popupInfo.=" pep. used $stateScopeLabel{$scope}";
		}
		else { # ratio => must be distinct
			# $popupInfo=($dispPepType eq 'all')? 'used' : 'distinct/used';
			$popupInfo.='Distinct';
		}
	}
	my $firstRatioUsed=0;
	foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
		$firstRatioUsed=$ratioPos unless $firstRatioUsed;
		print "<TH class=\"rbBorder\" colspan=2 nowrap>&nbsp;Ratio&nbsp;</TH><TH class=\"rbBorder\" nowrap>&nbsp;$pValueType&nbsp;</TH>";
		if ($ratioType=~/S\w+Ratio/) {
			print "<TH class=\"rbBorder\" nowrap onmouseover=\"popup('<B>All peptides used for ratio</B>')\" onmouseout=\"popout()\">&nbsp;Pep.&nbsp;</TH>";
			print "<TH class=\"rbBorder\" nowrap onmouseover=\"popup('<B>Peptides used for filtering:<BR>$popupInfo</B>')\" onmouseout=\"popout()\">&nbsp;Pep. filt.&nbsp;</TH>" if $pepTypeFiltering;
		}
		else {print $cvColStrg;}
	}
	print "</TR>\n";

	my $bgColor=$lightColor;
	my $ajaxPepAction=($quantifSoftware=~/SWATH|DIA|Spectronaut/)? 'ajaxPepSwath' : ($quantifSoftware eq 'MaxQuant')? 'ajaxPepMaxQuant' : ($protRulerCase)? 'ajaxPepProtRuler' : ($ratioType eq 'None')? 'ajaxPepAbund' : 'ajaxPepRatio';
	my $numDispItems=0;
	my %dispProteins;
	my ($pepCode4Filter,$pepFilterValue)=($dispPepType eq 'msmsPc')? ('FRAC_TRUE_USED',$dispNumPep/100) : ($numPepScope =~ /^replicate/ || $dispPepType =~ /^msms/)? ($numPepCode.'_SUBSET',$dispNumPep) : ($numPepCode,$dispNumPep); # numPepCode is GLOBAL
	my $sortItem=($dispSort eq 'peptide')? $pepCode4Filter : $dispSort;
	foreach my $modProtID (sort{&sortProteins($sortItem,$refProteinInfo,$refQuantifValues,$firstRatioUsed)} keys %{$refQuantifValues}) { #%{$refProteinInfo}
#last if $numDispItems >= 5;
		my ($protID,$modStrg)=($modProtID=~/^(\d+)(.*)/); # eg. 12345-Y40,T152 -> 12345, -Y40,T152
		$modStrg='' unless $modStrg;
		#>Peptide filter
		my $refUsedRatioPos=0;
		if ($view eq 'list') { # || $numPepScope eq 'replicate' # can be 'graph' if ajax call from graphical view
			my $okPeptide=0;
			foreach my $ratioPos (keys %dispRatios) {
				next if (!$refQuantifValues->{$modProtID}{$pepCode4Filter}{$ratioPos} || $refQuantifValues->{$modProtID}{$pepCode4Filter}{$ratioPos} < $pepFilterValue); # may not be defined for +/-inf ratios
				$okPeptide=1; # at least 1 ratio must pass the filter
				$refUsedRatioPos=$ratioPos;
				last; # unless $numPepScope eq 'replicate';
			}
			next unless $okPeptide;
		}
		my $mass=(!$refProteinInfo->{$protID}[1])? '-' : $refProteinInfo->{$protID}[1];
		my $anaList = $refProteinInfo->{$protID}[-1];
		my $quantifDataStrg='';
		my $okFilters=($view eq 'list')? 0 : 1; # 1 in case ajax call from graphical view
		my %protCV;
		foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
			$refUsedRatioPos=$ratioPos unless $refUsedRatioPos;
			if (!$refQuantifValues->{$modProtID} || !$refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos}) {
				$quantifDataStrg.="<TH colspan=2>-</TH><TH>-</TH>";
				$quantifDataStrg.="<TH>-</TH>" if $quantifType ne 'TNPQ'; # Std dev or SuperRatio prim ratio
				next;
			}
			# Ratio
			my $ratio=$refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos};
			my ($absRatio,$arrowDir)=($ratio >= 1)? ($ratio,'up_') : (1/$ratio,'down_');
			my $ratioStrg=($absRatio == 1000)? '&infin;' : sprintf "%.2f",$absRatio; # '∞'
			$ratioStrg='1/'.$ratioStrg if ($ratio < 1 && $ratioStrg ne '1.00');
			# Arrow flag
			my $arrowStrg='';
			my $okRatio=($ratio <= $minFoldChange || $ratio >= $maxFoldChange)? 1 : 0;
			my $okPvalue=(defined $refQuantifValues->{$modProtID}{$pvalueCode} && defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} && $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} <= $maxPvalue)? 1 : 0;
			if ($okRatio || $okPvalue) {
				my $arrowColor=(!$okRatio || !$okPvalue)? 'gray' : ($ratio >= 1)? 'red' : 'green';
				$arrowStrg="<IMG class=\"LINK\" src=\"$promsPath{images}/$arrowDir$arrowColor.png\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'$ajaxPepAction')\">";
			}
			#if ($ratio <= $minFoldChange || $ratio >= $maxFoldChange) {
			#	my $arrowColor=(!$refQuantifValues->{$modProtID}{$pvalueCode} || !defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} || $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} > $maxPvalue)? 'gray' : ($ratio >= 1)? 'red' : 'green';
			#	$arrowStrg="<IMG class=\"LINK\" src=\"$promsPath{images}/$arrowDir$arrowColor.png\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'$ajaxPepAction')\">";
			#}
			$quantifDataStrg.="<TH nowrap align=right>&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'$ajaxPepAction')\"/>$ratioStrg</A></TH><TD>$arrowStrg&nbsp;</TD>";
			# p-value
			#my $pValueStrg=(!defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos})? '-' : ($refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}>=0.01)? sprintf "%.2f",$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} : sprintf "%.1e",$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos};

			if (defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}) {
				my $pValueStrg=($refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}>=0.01)? sprintf "%.2f",$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} : sprintf "%.1e",$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos};
				$quantifDataStrg.="<TH nowrap>&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'ajaxProtStat')\">$pValueStrg</A>&nbsp;</TH>";
			}
			else {$quantifDataStrg.="<TH nowrap>&nbsp;-&nbsp;</TH>";}
			# Std dev OR num Pep
			if ($ratioType eq 'Ratio') { # Algo v1
				if ($quantifType ne 'TNPQ') { # Std dev
					if ($refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}) {
						#$protCV{$ratioPos}=($ratioType eq 'Ratio')? abs(log($refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos})/log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos})) : $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}/abs(log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos})/$log2);
						$protCV{$ratioPos}=abs(log($refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos})/log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos}));
						my $cvStrg=sprintf "%.1f",($protCV{$ratioPos}*100);
						$quantifDataStrg.="<TH nowrap>&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'ajaxProtStat')\">$cvStrg</A>&nbsp;</TH>";
					}
					else {$quantifDataStrg.="<TH nowrap>&nbsp;-&nbsp;</TH>";}
				}
			}
			else {
				my $cellTag=($refQuantifValues->{$modProtID}{$pepCode4Filter} && $refQuantifValues->{$modProtID}{$pepCode4Filter}{$ratioPos} >= $pepFilterValue)? 'TH' : 'TD';
				$quantifDataStrg.="<$cellTag align=center nowrap>&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'$ajaxPepAction')\"/>".$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$ratioPos}."&nbsp;</$cellTag>\n";
				if ($pepTypeFiltering) { # 1 extra column
					$quantifDataStrg.="<$cellTag align=center nowrap>&nbsp;";
					if ($dispPepType =~ /^msms/ && !$refQuantifValues->{$modProtID}{$pepCode4Filter}{$ratioPos}) { # not stored if 0
						$quantifDataStrg.='-';
					}
					elsif ($numPepScope =~ /^replicate/ || $dispPepType =~ /^msms/) {
						$quantifDataStrg.='<SPAN '; #'<A href="javascript:void(null)"';
						#my $tgtPosSwap=($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos}==0.001 || $refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos}==1000)? -1 : 1; # Inf ratio ALWAYS takes negative tgt pos since 21/09/20
						if ($numPepScope eq 'ratio') { # => $dispPepType=msms(Num|Pc)
							$quantifDataStrg.=" onmouseover=\"popup('<B>Valid. pep. in ratio: ".$refQuantifValues->{$modProtID}{$numPepCode.'_SUBSET'}{$ratioPos}."</B>')\" onmouseout=\"popout()\"" if $dispPepType eq 'msmsPc'; # no popup for msmsNum
						}
						else { # state-level peptide filtering
							my %pepRepStrg=(0=>'',1=>'');
							foreach my $s (0,1) {
								foreach my $i (0..$#{$dispRatios{$ratioPos}[$s]}) { # only [0] for msms
									#my $targetPos=$dispRatios{$ratioPos}[$s][$i]*$tgtPosSwap;
									my $targetPos=$dispRatios{$ratioPos}[$s][$i];
									$targetPos*=-1 if (($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos}==0.001 || $refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos}==1000) && $refQuantifValues->{$modProtID}{$numPepCode} && defined $refQuantifValues->{$modProtID}{$numPepCode}{-$targetPos}); # if Inf ratio: take negative tgt pos if exists ()
									$pepRepStrg{$s}.=', ' if $pepRepStrg{$s} ne '';
									if ($dispPepType eq 'msmsPc') {
										my ($numPep,$fracPep)=($refQuantifValues->{$modProtID}{$numPepCode} && $refQuantifValues->{$modProtID}{$numPepCode}{$targetPos})? ($refQuantifValues->{$modProtID}{$numPepCode}{$targetPos},$refQuantifValues->{$modProtID}{'FRAC_TRUE_USED'}{$targetPos}) : (0,0);
										my $pcPep=1*(sprintf "%.1f",$fracPep*100);
										my $pcPepStrg=($pcPep)? ' ('.$pcPep.'%)' : '';
										$pepRepStrg{$s}.=($fracPep >= $pepFilterValue)? "<B>$numPep$pcPepStrg</B>" : "$numPep$pcPepStrg";
									}
									else {
										my $numPep=($refQuantifValues->{$modProtID}{$numPepCode} && $refQuantifValues->{$modProtID}{$numPepCode}{$targetPos})? $refQuantifValues->{$modProtID}{$numPepCode}{$targetPos} : 0;
										$pepRepStrg{$s}.=($numPep >= $dispNumPep)? '<B>'.$numPep.'</B>' : $numPep;
									}
								}
							}
							my $item=($numPepScope =~ /^replicate/)? 'replicates' : 'states';
							$quantifDataStrg.=" onmouseover=\"popup('<B>$peptideTitle in $item:</B><BR>&nbsp;$pepRepStrg{1} / $pepRepStrg{0}')\" onmouseout=\"popout()\"";
						}
						# $quantifDataStrg.='/>'.$refQuantifValues->{$modProtID}{$pepCodeDisplay}{$ratioPos};
						my $dispPepValue=($dispPepType eq 'msmsPc')? 1*(sprintf "%.1f",$refQuantifValues->{$modProtID}{$pepCode4Filter}{$ratioPos}*100) : $refQuantifValues->{$modProtID}{$pepCode4Filter}{$ratioPos};
						$quantifDataStrg.="/>$dispPepValue</SPAN>";
		# 				if ($dispPepType =~ /^msms/) {
		# $quantifDataStrg.='('.($refQuantifValues->{$modProtID}{$numPepCode.'_SUBSET'}{$ratioPos} || 0).')';
		# 					#my $numGhostPep=$refQuantifValues->{$modProtID}{'NUM_PEP_USED_STATE'}{$ratioPos}-($refQuantifValues->{$modProtID}{$numPepCode.'_SUBSET'}{$ratioPos} || 0);
		# 					my $numGhostPep=$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$ratioPos}-($refQuantifValues->{$modProtID}{$pepCodeDisplay}{$ratioPos} || 0);
		# 					$quantifDataStrg.='+<FONT class="ghostPeptide">'.$numGhostPep.'</FONT>' if $numGhostPep;
		# 				}
		# $quantifDataStrg.='['.$refQuantifValues->{$modProtID}{$pepCode4Filter}{$ratioPos}.']';
						#$quantifDataStrg.='/'.$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$ratioPos} if $numPepCode eq 'NUM_TRUE_USED';
					}
					else { # normal dist. pep
						$quantifDataStrg.=$refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos}
					}


					# 	$quantifDataStrg.="/>";
					# 	$quantifDataStrg.=($numPepCode eq 'DIST_PEP_USED')? $refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos}.'/'.$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$ratioPos} : $refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos};
					# }
					$quantifDataStrg.="&nbsp;</$cellTag>";
				}
			}
			unless ($okFilters) { # no ratioPos has yet passed filters
				if ($foldChangeType eq 'abs') {
					$okFilters=(($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} < 1 && $refQuantifValues->{$modProtID}{RATIO}{$ratioPos} <= $invDispFoldChange) || ($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos} >= 1 && $refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos} >= $dispFoldChange))? 1 : 0;
				}
				elsif ($foldChangeType eq 'up') {
					$okFilters=($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} >= 1 && $refQuantifValues->{$modProtID}{RATIO}{$ratioPos} >= $dispFoldChange)? 1 : 0;
				}
				else { # down
					$okFilters=($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} <= 1 && $refQuantifValues->{$modProtID}{RATIO}{$ratioPos} <= $invDispFoldChange)? 1 : 0;
				}
				$okFilters=($okFilters && ($dispPvalue>=1 || ($allowInfRatio && $absRatio==1000) || (defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} && $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} <= $dispPvalue)))? 1 : 0;
				#if ($okFilters && $dispConfInt && defined $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}) { # undefined always passes filter!
				#	$okFilters=(1.96*$refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos} <= $dispConfInt)? 1 : 0;
				#}
				if ($okFilters && $dispCV && defined $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}) { # undefined always passes filter!
					#$okFilters=0 if ($refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}/abs(log($refQuantifValues->{$modProtID}{RATIO}{$ratioPos}/$log2)) > $dispCVFrac);
					$okFilters=0 if $protCV{$ratioPos} > $dispCVFrac;
				}
			}
		}
#print "$modProtID<BR>\n" unless $okFilters;
		next unless $okFilters;
		$numDispItems++;
		$dispProteins{$protID}=1;
		my $dispModIdent=($refDispModifSites->{$modProtID})? $refProteinInfo->{$protID}[0].'-'.$refDispModifSites->{$modProtID} : $refProteinInfo->{$protID}[0].$modStrg;
		print "<TR bgcolor=\"$bgColor\" valign=top><TD class=\"TH\" nowrap><INPUT type=\"checkbox\" name=\"chkProt\" value=\"$modProtID\"/><A href=\"javascript:sequenceView($protID,'$anaList')\">$dispModIdent</A>&nbsp;</TD>\n";
		# Gene(s)
		if (scalar @{$refProteinInfo->{$protID}[5]} > 1) {
			print "<TH>&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-",join('<BR>&nbsp;&nbsp;-',@{$refProteinInfo->{$protID}[5]}[1..$#{$refProteinInfo->{$protID}[5]}]),"</B>')\" onmouseout=\"popout()\">",$refProteinInfo->{$protID}[5][0],"</A>&nbsp;</TH>";
		}
		elsif ($refProteinInfo->{$protID}[5][0]) {print '<TH>&nbsp;',$refProteinInfo->{$protID}[5][0],'&nbsp;</TH>';}
		else {print '<TH>-</TH>';} # no gene mapped
		# quantif data
		print $quantifDataStrg;
		if ($refQuantifValues->{$modProtID} && $quantifType ne 'MQ') {# For SILAC quantifications made with MaxQuant, only redundant peptides are provided in proteinGroups.txt (column called 'Ratio H/L count')
			#my $distPepStrg=($refQuantifValues->{$modProtID}{'DIST_PEP_USED'}==$refQuantifValues->{$modProtID}{'NUM_PEP_USED'})? '' : $refQuantifValues->{$modProtID}{'DIST_PEP_USED'}.'/';
			#my $distPepUsedStrg=($refQuantifValues->{$modProtID}{DIST_PEP_USED})? $refQuantifValues->{$modProtID}{DIST_PEP_USED} : '-'; #
			my $distPepUsedStrg=($ratioType=~/S\w+Ratio/)? $refQuantifValues->{$modProtID}{NUM_PEP_TOTAL}{$refUsedRatioPos} : ($numPepCode eq 'DIST_PEP_USED')? "$refQuantifValues->{$modProtID}{DIST_PEP_USED}{$refUsedRatioPos}/$refQuantifValues->{$modProtID}{NUM_PEP_USED}{$refUsedRatioPos}/$refQuantifValues->{$modProtID}{NUM_PEP_TOTAL}{$refUsedRatioPos}": "$refQuantifValues->{$modProtID}{NUM_PEP_USED}{$refUsedRatioPos}/$refQuantifValues->{$modProtID}{NUM_PEP_TOTAL}{$refUsedRatioPos}"; #
			print "<TH>&nbsp;$distPepUsedStrg&nbsp;</TH>";
		}
		else {print "<TH>-</TH>"}
		# mw, des & species
		print "<TD class=\"TH\" align=right>$mass&nbsp;</TD><TD>$refProteinInfo->{$protID}[2] <U><I>$refProteinInfo->{$protID}[3]</I></U></TD>\n";
		print "</TR>\n";
		$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
	}
	$protString=($isModifQuantif)? "&nbsp;$numDispItems/$numTotQuantifItems sites<BR>&nbsp;".(scalar keys %dispProteins)."/$numTotProteins proteins&nbsp;" : "&nbsp;$numDispItems/$numTotProteins proteins&nbsp;";
	print qq
|<TR><TD colspan=$numColumns><B>End of list.</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<IMG src="$promsPath{images}/up_red.png"/><IMG src="$promsPath{images}/down_green.png"/>: Fold change &ge; 2 <B>and</B> $pValueType &le; 0.05&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<IMG src="$promsPath{images}/up_gray.png"/><IMG src="$promsPath{images}/down_gray.png"/>: Fold change &ge; 2 <B>or</B> $pValueType &le; 0.05</TD></TR>
</TABLE>
</FORM>
<SCRIPT type="text/javascript"><!-- Not executed if called from AJAX -->
document.getElementById('protCountDIV').innerHTML='$protString';
</SCRIPT>
|;
}

sub listProteinIntensities { # Abundance or Maxquant (including MaxQuant SILAC ratio!!!)
	my ($labelType,$refLabelingInfo,$refStateInfo,$refQuantifValues,$refProteinInfo,$refDispModifSites)=@_; # $labelType  not global for all calls
	my $quantifSoftware=($refLabelingInfo->{'SOFTWARE'})? $refLabelingInfo->{'SOFTWARE'}[0] : 'myProMS'; $quantifSoftware='MaxQuant' if $quantifSoftware eq 'MQ';
	my $measureName=$quantifParamInfo{$dispMeasure}[1]; # only for no-ratio
	my $invDispFoldChange=1/$dispFoldChange;
	my $numTotQuantifItems=scalar keys %{$refQuantifValues};
	my $numTotProteins=scalar keys %{$refProteinInfo};
	my $protString=($isModifQuantif)? "&nbsp;$numTotQuantifItems sites&nbsp;<BR>&nbsp;$numTotProteins proteins&nbsp;" : "&nbsp;$numTotProteins proteins&nbsp;";
	my $disabSave=($projectAccess eq 'guest')? ' disabled' : '';
	my $clearButtonStrg=($view eq 'list')? '' : '<INPUT type="button" value="Clear list" style="font-weight:bold" onclick="document.getElementById(\'protListDIV\').style.display=\'none\'">';
	my $refDispSets=($ratioType eq 'None')? \%dispStates : \%dispRatios;
	my $numDispSets=scalar keys %{$refDispSets};
	my $numColumns=($quantifType eq 'PROT_ABUNDANCE' && exists($absCIsCodes{$dispMeasure}))? $numDispSets * 3 + 5 :
				   ($quantifType eq 'PROT_ABUNDANCE')? $numDispSets * 2 + 5 :
				   ($ratioType eq 'None')? $numDispSets + 7 : 
				   $numDispSets * 3 + 6;
	my $selProtRequired = ($action eq 'ajaxListProt') ? 'true' : 'false';
	print qq
|<FORM name="protForm" method="POST">
<TABLE border=0 cellspacing=0 cellpadding=2>
<TR><TD colspan=$numColumns><DIV id="saveProtDIV" style="display:none;"></DIV></TD></TR>
<TR><TD colspan=$numColumns>$clearButtonStrg
<INPUT type="button" value="Check all" onclick="checkAllProteins(true)"/>
<INPUT type="button" value="Uncheck all" onclick="checkAllProteins(false)"/>
<INPUT type="button" id="saveFormBUTTON" value="Save proteins..." onclick="ajaxManageSaveProteins('getThemes','PROT')"$disabSave/>
|;
	print "<INPUT type=\"button\" id=\"saveSiteFormBUTTON\" value=\"Save sites...\" onclick=\"ajaxManageSaveProteins('getThemes','SITE')\"$disabSave/>\n" if $isModifQuantif;
	print "<INPUT type=\"button\" value=\"Export data\" onclick=\"checkProtToExport($selProtRequired)\"/>\n";
	print qq
|</TD></TR>
<TR bgcolor="$darkColor">
<TH class="rbBorder" align=left rowspan=2><DIV id="protCountDIV">$protString</DIV></TH><TH class="rbBorder" rowspan=2>Gene</TH>
|;
	if ($ratioType eq 'None') {
		if ($quantifType eq 'PROT_ABUNDANCE') {
			my $colSpanNb = (exists($absCIsCodes{$dispMeasure}))? 3 : 2;
			$colSpanNb++ if $dispPepType ne 'all'; # always displayed
			foreach my $statePos (sort{$a<=>$b} keys %dispStates) {
				print "<TH class=\"rbBorder\" nowrap colspan=$colSpanNb>&nbsp;$refStateInfo->{$statePos}{NAME}&nbsp;</TH>";
			}
		}
		else {
			my $stateStrg=($numDispSets > 1)? 'States' : 'State';
			print "<TH class=\"rbBorder\" colspan=$numDispSets>&nbsp$stateStrg&nbsp;</TH>";
		}
	}
	else { # MaxQuant ratio
		foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
			my ($testStateID,$refStateID)=($ratioPos > 0)? split(/\//,$refLabelingInfo->{'RATIOS'}[$ratioPos-1]) : reverse(split(/\//,$refLabelingInfo->{'RATIOS'}[abs($ratioPos)-1]));
			print "<TH class=\"rbBorder\" colspan=2 nowrap>&nbsp;$refStateInfo->{$condToState{$testStateID}}{NAME}/$refStateInfo->{$condToState{$refStateID}}{NAME}&nbsp;</TH>\n";
		}
	}
	if ($quantifType eq 'PROT_ABUNDANCE') {
		print qq |<TH class="rbBorder" rowspan=2 nowrap>All<BR>&nbsp;peptides&nbsp;</TH>|;
	}
	else {
		print qq |<TH class="rbBorder" colspan=3 nowrap>&nbsp;Peptides&nbsp;</TH>|;
	}
	print qq
|<TH class="rbBorder" rowspan=2 nowrap>&nbsp;MW<SMALL> kDa</SMALL>&nbsp;</TH><TH class="bBorder" rowspan=2 width=700 nowrap>&nbsp;Description - Species&nbsp;</TH>
</TR>
<TR bgcolor="$darkColor">
|;
	my $firstSetUsed;
	if ($ratioType eq 'None') {
		my ($pepStrg,$popupInfo)=($dispPepType eq 'all')? ('Pep.','All peptides used for state') : ($dispPepType  eq 'msmsPc')? ('% valid. pep.','% validated peptides used') : ($dispPepType  =~ /^msms/)? ('Valid. pep.','Validated peptides used') : ('Dist. pep.','Distinct peptides used');
		foreach my $statePos (sort{$a<=>$b} keys %dispStates) {
			$firstSetUsed=$statePos unless $firstSetUsed;
			if ($quantifType eq 'PROT_ABUNDANCE') {
				my $measName=($dispMeasure eq 'MY_LFQ')? 'LFQ' : ucfirst(lc($dispMeasure)).'.'; $measName=~s/_/ /;
				print "<TH class=\"rbBorder\" nowrap>&nbsp;$measName&nbsp;</TH>\n";
				print "<TH class=\"rbBorder\" nowrap>&nbsp;95% CI&nbsp;</TH>\n" if (exists($absCIsCodes{$dispMeasure}));
				print "<TH class=\"rbBorder\" nowrap onmouseover=\"popup('<B>All peptides used for state</B>')\" onmouseout=\"popout()\">&nbsp;Pep.&nbsp;</TH>\n";
				print "<TH class=\"rbBorder\" nowrap onmouseover=\"popup('<B>$popupInfo</B>')\" onmouseout=\"popout()\">&nbsp;$pepStrg&nbsp;</TH>\n" if $dispPepType ne 'all';
			}
			else {print "<TH class=\"rbBorder\" nowrap>&nbsp;$refStateInfo->{$statePos}{NAME}&nbsp;</TH>";}
		}
	}
	else {
		my $ratioTitle=($foldChangeRaw)? 'Raw ratio' : 'Norm. ratio';
		foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
			$firstSetUsed=$ratioPos unless $firstSetUsed;
			print "<TH class=\"rbBorder\" nowrap>&nbsp;$ratioTitle&nbsp;</TH><TH class=\"rbBorder\" nowrap>&nbsp;Ratio count&nbsp;</TH>\n";
		}
	}
	# my $quantifSoftware=($refLabelingInfo->{'SOFTWARE'})? $refLabelingInfo->{'SOFTWARE'}[0] : 'myProMS';
	# $quantifSoftware='MaxQuant' if $quantifSoftware eq 'MQ';
	if ($quantifSoftware eq 'MaxQuant') {
		print qq
|<TH class="rbBorder"nowrap>&nbsp;Unique&nbsp;</TH><TH class="rbBorder" nowrap>&nbsp;Razor+uniq.&nbsp;</TH><TH class="rbBorder" nowrap>&nbsp;All&nbsp;</TH>
|;
	}
	print "</TR>\n";

	my $bgColor=$lightColor;
	my $numDispItems=0;
	my %dispProteins;
	my ($pepCode4Filter,$pepFilterValue)=($dispPepType eq 'msmsPc')? ('FRAC_TRUE_USED',$dispNumPep/100) : ($numPepCode,$dispNumPep); # numPepCode is GLOBAL
	my $sortItem=($dispSort eq 'peptide')? $pepCode4Filter : $dispSort;
	foreach my $modProtID (sort{&sortProteins($sortItem,$refProteinInfo,$refQuantifValues,$firstSetUsed)} keys %{$refQuantifValues}) { #%{$refProteinInfo}
		my ($protID,$modResStrg)=($modProtID=~/^(\d+)(.*)/); # eg. 12345-Y40,T152 -> 12345, -Y40,T152
		$modResStrg='' unless $modResStrg;
		#>Peptide filter
		if ($view eq 'list') { # can also be called from graph view by ajax on subset of proteins => no filtering
			my $okPeptide=0;
			foreach my $setPos (keys %{$refDispSets}) {
				next if (!$refQuantifValues->{$modProtID}{$pepCode4Filter} || !$refQuantifValues->{$modProtID}{$pepCode4Filter}{$setPos} || $refQuantifValues->{$modProtID}{$pepCode4Filter}{$setPos} < $pepFilterValue); # may not be defined for +/-inf ratios
				$okPeptide=1; # at least 1 state must pass the filter
				last;
			}
			next unless $okPeptide;
		}
		my $mass=(!$refProteinInfo->{$protID}[1])? '-' : $refProteinInfo->{$protID}[1];
		my $anaList = $refProteinInfo->{$protID}[-1];
		my $quantifDataStrg='';
		if ($ratioType eq 'None') {
			foreach my $statePos (sort{$a<=>$b} keys %dispStates) {
				if ($refQuantifValues->{$modProtID} && $refQuantifValues->{$modProtID}{$dispMeasure}{$statePos}) {
					$quantifDataStrg.=sprintf "<TH nowrap>&nbsp;%.2g&nbsp;</TH>",$refQuantifValues->{$modProtID}{$dispMeasure}{$statePos};
					my ($dispPepValue,$numPepValue)=($refQuantifValues->{$modProtID}{$numPepCode}{$statePos})? ($refQuantifValues->{$modProtID}{$pepCode4Filter}{$statePos},$refQuantifValues->{$modProtID}{$numPepCode}{$statePos}) : (0,0);
					my $cellTag=($dispPepValue >= $pepFilterValue)? 'TH' : 'TD';
					if ($quantifType eq 'PROT_ABUNDANCE') {
						if (exists($absCIsCodes{$dispMeasure})){  # Display 95% confidence interval
							if ($refQuantifValues->{$modProtID}{$absCIsCodes{$dispMeasure}->[0]}{$statePos}) {
								$quantifDataStrg.=sprintf "<TH nowrap>&nbsp;%.2g - %.2g&nbsp;</TH>",
								$refQuantifValues->{$modProtID}{$absCIsCodes{$dispMeasure}->[0]}{$statePos},
								$refQuantifValues->{$modProtID}{$absCIsCodes{$dispMeasure}->[1]}{$statePos};
							} else {
								$quantifDataStrg.="<TH>&nbsp;-&nbsp;</TH>";
							}
						}
						my $cellTag1=($dispPepType eq 'all')? $cellTag : 'TH';
						$quantifDataStrg.="<$cellTag1 align=\"center\" nowrap>&nbsp;".$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$statePos}."&nbsp;</$cellTag1>";
						if ($dispPepType ne 'all') {
							$quantifDataStrg.="<$cellTag align=\"center\" nowrap>&nbsp;";
							if ($dispPepType eq 'msmsPc' && $numPepValue) {
								my $pcPep=1*(sprintf "%.1f",$dispPepValue*100);
								my ($bold1,$bold2)=($dispPepValue >= $pepFilterValue)? ('<B>','</B>') : ('','');
								$quantifDataStrg.="<SPAN onmouseover=\"popup('".$bold1."Valid. pep. used: $numPepValue$bold2')\" onmouseout=\"popout()\">$pcPep</SPAN>";
							}
							else {$quantifDataStrg.=($dispPepValue)? $dispPepValue : '-';}
							$quantifDataStrg.="&nbsp;</$cellTag>";
						}
						
						# my $extraPepStrg=($numPepCode =~ /(DIST_PEP|NUM_TRUE)_USED/)? '/'.$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$statePos} : '';
						# $quantifDataStrg.="<$cellTag align=\"center\">&nbsp;$numPepValue$extraPepStrg&nbsp;</$cellTag>";
					}
				}
				else {
					$quantifDataStrg.='<TH>-</TH>'; # meas
					if ($quantifType eq 'PROT_ABUNDANCE') {
						$quantifDataStrg.='<TH>-</TH>' if exists($absCIsCodes{$dispMeasure});
						$quantifDataStrg.='<TD align="center">-</TD>'; # NUM_PEP_USED
						$quantifDataStrg.='<TD align="center">-</TD>' if $dispPepType ne 'all';
					}
				}
			}
		}
		else { # MQ ratio
			my $okFilters=($view eq 'list')? 0 : 1; # 1 in case ajax call from graphical view
			foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
				if (!$refQuantifValues->{$modProtID} || !$refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos}) {
					$quantifDataStrg.="<TH>-</TH><TH>-</TH>";
					next;
				}
				unless ($okFilters) { # no ratioPos has yet passed filters
					if ($foldChangeType eq 'abs') {
						$okFilters=(($refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} < 1 && $refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} <= $invDispFoldChange) || ($refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} >= 1 && $refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} >= $dispFoldChange))? 1 : 0;
					}
					elsif ($foldChangeType eq 'up') {
						$okFilters=($refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} >= 1 && $refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} >= $dispFoldChange)? 1 : 0;
					}
					else { # down
						$okFilters=($refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} <= 1 && $refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} <= $invDispFoldChange)? 1 : 0;
					}
					if ($okFilters && $dispCV && $refQuantifValues->{$modProtID}{'RATIO_VAR'} &&  defined $refQuantifValues->{$modProtID}{'RATIO_VAR'}{$ratioPos}) { # undefined always passes filter!
						$okFilters=0 if $refQuantifValues->{$modProtID}{'RATIO_VAR'}{$ratioPos} > $dispCV;
					}
					if ($okFilters && $refQuantifValues->{$modProtID}{$numPepCode} && $refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos}) {
						$okFilters=0 if $refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos} < $dispNumPep;
					}
				}
				my $ratioStrg=($refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} >= 1)? sprintf "%.2f",$refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} : '1/'.sprintf "%.2f",1/$refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos};
				#$quantifDataStrg.="<TH nowrap>&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'ajaxProtStat')\">$ratioStrg</A>&nbsp;</TH><TH>$refQuantifValues->{$modProtID}{NUM_PEP_USED}{$ratioPos}</TH>";
				$quantifDataStrg.="<TH nowrap>&nbsp;$ratioStrg&nbsp;</TH><TH>$refQuantifValues->{$modProtID}{NUM_PEP_USED}{$ratioPos}</TH>";
			}
			next unless $okFilters;
		}
		$numDispItems++;
		$dispProteins{$protID}=1;
		my $dispModIdent=($refDispModifSites->{$modProtID})? $refProteinInfo->{$protID}[0].'-'.$refDispModifSites->{$modProtID} : $refProteinInfo->{$protID}[0].$modResStrg;
		my $class=($anaList)? 'TH' : 'TD'; # quantified MaxQuant proteins can be hidden in myProMS
		print "<TR bgcolor=\"$bgColor\" valign=top><TD class=\"$class\" nowrap><INPUT type=\"checkbox\" name=\"chkProt\" value=\"$modProtID\"/><A href=\"javascript:sequenceView($protID,'$anaList')\">$dispModIdent</A>&nbsp;</TD>\n";
		# Gene(s)
		if (scalar @{$refProteinInfo->{$protID}[5]} > 1) {
			print "<TH>&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-",join('<BR>&nbsp;&nbsp;-',@{$refProteinInfo->{$protID}[5]}[1..$#{$refProteinInfo->{$protID}[5]}]),"</B>')\" onmouseout=\"popout()\">",$refProteinInfo->{$protID}[5][0],"</A>&nbsp;</TH>";
		}
		elsif ($refProteinInfo->{$protID}[5][0]) {print '<TH>&nbsp;',$refProteinInfo->{$protID}[5][0],'&nbsp;</TH>';}
		else {print '<TH>-</TH>';} # no gene mapped
		# quantif data
		print $quantifDataStrg;
		# peptides (all states)
		if ($quantifType eq 'PROT_ABUNDANCE') {
			my $allPeptideStrg=$refQuantifValues->{$modProtID}{'NUM_PEP_TOTAL'}{$firstSetUsed} || '-';
			print "<TH>$allPeptideStrg</TH>";
		}
		else {
			my $uniqStrg=$refQuantifValues->{$modProtID}{'UNIQUE_PEP'}{$firstSetUsed} || '-';
			my $razUniStrg=$refQuantifValues->{$modProtID}{'RAZ_UNI_PEP'}{$firstSetUsed} || '-';
			my $peptideStrg=$refQuantifValues->{$modProtID}{'PEPTIDES'}{$firstSetUsed} || '-';
			print "<TH>$uniqStrg</TH><TH>$razUniStrg</TH><TH>$peptideStrg</TH>";
		}
		# mw, des & species
		print "<TD class=\"TH\" align=right>$mass&nbsp;</TD><TD>$refProteinInfo->{$protID}[2] <U><I>$refProteinInfo->{$protID}[3]</I></U></TD>\n";
		print "</TR>\n";
		$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
	}
	$protString=($isModifQuantif)? "&nbsp;$numDispItems/$numTotQuantifItems sites<BR>&nbsp;".(scalar keys %dispProteins)."/$numTotProteins proteins&nbsp;" : "&nbsp;$numDispItems/$numTotProteins proteins&nbsp;";
	print qq
|</TABLE>
</FORM>
<SCRIPT type="text/javascript"><!-- Not executed if called from AJAX -->
document.getElementById('protCountDIV').innerHTML='$protString';
</SCRIPT>
|;
}


sub displayProtRulerPlot {
	my ($refLabelingInfo,$refStateInfo,$refQuantifValues,$refProteinInfo)=@_; # globals: $dispMeasure,$numPepCode,%quantifParamInfo
	my $measureName=$quantifParamInfo{$dispMeasure}[1];
	my $isLogged=1;
	my $axisYtitle=($isLogged)? (1,"Log10($measureName)") : (0,$measureName);
	my $log10=log(10);
	print qq
|</CENTER>
<SCRIPT type="text/javascript">
function myPointLabel(dp,type) {
	var infoStrg=dp.label;
	if (type != 'min') {
		if (dp.dataSet.params.chart.dataSets.length > 1) {infoStrg+='\\nState: '+dp.dataSet.params.name;}
		var intensity=($isLogged==1)? (1*(Math.pow(10,dp.getY())).toFixed(1)) : dp.getY();
		infoStrg+='\\n$measureName='+intensity+'\\nrank='+dp.getX();
	}
	return infoStrg;
}
var mainChart,IP;
window.onload=function() {
	ajaxUpdateRestrict();
	IP=new genericPlot({div:'mainGraphDIV',width:600,height:450,
						axisX:{title:'Rank'}, // ,zoomable:true
						axisY:{title:'$axisYtitle'}, // ,zoomable:true
						axisClosure:true,
						//hideDataSets:true, // hide dataset label if only 1 displayed!
						zoomable:true,
						pointLabel:myPointLabel,
						pointOnclick:proteinLink,
						pointOnList:ajaxListSelectedProteins,
						allowHighlight:true,
						searchable:{text:'Extended search',externalSearch:ajaxSearchConvertIdentifier},
						updateHighlight:{callback:updateUsedColors,editable:false},
						exportAsImage:['Export as image','IntensityPlot','./exportSVG.cgi']
					});
	mainChart=IP; // needed for highlighting
|;
	my $dsetIdx=0;
	if ($refLabelingInfo->{'AVG_MODE'}[0] eq 'AVG_ALL') {
		print "\tIP.addDataSet($dsetIdx,{name:'Average on all states'});\n";
	} else {
		foreach my $statePos (sort{$a<=>$b} keys %dispStates) {
			print "\tIP.addDataSet($dsetIdx,{name:'$refStateInfo->{$statePos}{NAME}'});\n";
			$dsetIdx++;
		}
	}
	my $protCount;
	$dsetIdx=0;
	my %dispProteins;

	if ($refLabelingInfo->{'AVG_MODE'}[0] eq 'AVG_ALL') {
		my $statePos=0;
		my $count=0;
		$protCount=0;
		foreach my $protID (sort{&sortCheckDefined($dispMeasure,$refQuantifValues,$statePos,1) || &sortCheckDefined($numPepCode,$refQuantifValues,$statePos,1)} keys %{$refQuantifValues}) {
			next if (!$refQuantifValues->{$protID}{$dispMeasure} || !$refQuantifValues->{$protID}{$dispMeasure}{$statePos});
			if ($count==0) {
				print "\tIP.addDataAsString($dsetIdx,'";
			}
			else {print ";";}
			$count++;
			$protCount++;
			$refProteinInfo->{$protID}[0]=~s/[,;']/\./g; # Clean MaxQuant crappy contaminant identifiers
			print "$refProteinInfo->{$protID}[0],$protID,$protCount,";
			if ($isLogged) {
				print log($refQuantifValues->{$protID}{$dispMeasure}{$statePos})/$log10;
			} else {
				print $refQuantifValues->{$protID}{$dispMeasure}{$statePos};
			}
			print ',10';  # Size of points
			if ($count==100) {
				print "');\n";
				$count=0;
			}
			$dispProteins{$protID}=1;
		}
		if ($count > 0) {
			print "');\n";
		}
	} else {
		foreach my $statePos (sort{$a<=>$b} keys %dispStates) {
			my $count=0;
			$protCount=0;
			foreach my $protID (sort{&sortCheckDefined($dispMeasure,$refQuantifValues,$statePos,1) || &sortCheckDefined($numPepCode,$refQuantifValues,$statePos,1)} keys %{$refQuantifValues}) {
				next if (!$refQuantifValues->{$protID}{$dispMeasure} || !$refQuantifValues->{$protID}{$dispMeasure}{$statePos});
				#next if (!$refQuantifValues->{$protID}{$numPepCode} || !$refQuantifValues->{$protID}{$numPepCode}{$statePos} || $refQuantifValues->{$protID}{$numPepCode}{$statePos} < $dispNumPep); # Not in DB if ==0 (eg. UNIQUE_PEP)
				#next if ($dispConfInt && defined $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos} && 1.96*$refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos} > $dispConfInt);
				if ($count==0) {
					print "\tIP.addDataAsString($dsetIdx,'";
				}
				else {print ";";}
				$count++;
				$protCount++;
				$refProteinInfo->{$protID}[0]=~s/[,;']/\./g; # Clean MaxQuant crappy contaminant identifiers
				print "$refProteinInfo->{$protID}[0],$protID,$protCount,";
				if ($isLogged) {print log($refQuantifValues->{$protID}{$dispMeasure}{$statePos})/$log10;} else {print $refQuantifValues->{$protID}{$dispMeasure}{$statePos};}
				#print ','.$refQuantifValues->{$protID}{$numPepCode}{$statePos};  # Number of peptides for size of points
				print ',10';  # Size of points
				if ($count==100) {
					print "');\n";
					$count=0;
				}
				$dispProteins{$protID}=1;
			}
			if ($count > 0) {
				print "');\n";
			}
			$dsetIdx++;
		}
	}
	my $numProtDisp=scalar keys %dispProteins;
	print qq
|	IP.draw();
}
</SCRIPT>
<TABLE cellspacing=0 cellpadding=0>
<TR><TD><FONT class="title3">$numProtDisp proteins displayed</FONT></TR></TD>
<TR><TD><DIV id="mainGraphDIV" style="float:left;background-color:#FFFFFF;border:solid 3px $darkColor;padding:2px;"></DIV></TR></TD>
</TABLE>
<DIV id="protListDIV" style="position:absolute;height:535;overflow:auto"></DIV>
|;
}

sub listProteomicRuler {
	my ($labelType,$refLabelingInfo,$refStateInfo,$refQuantifValues,$refProteinInfo)=@_;
	my $measureName=$quantifParamInfo{$dispMeasure}[1];
	my $numTotQuantifItems=scalar keys %{$refQuantifValues};
	my $numTotProteins=scalar keys %{$refProteinInfo};
	my $protString=($isModifQuantif)? "&nbsp;$numTotQuantifItems sites&nbsp;<BR>&nbsp;$numTotProteins proteins&nbsp;" : "&nbsp;$numTotProteins proteins&nbsp;";
	my $disabSave=($projectAccess eq 'guest')? ' disabled' : '';
	my $clearButtonStrg=($view eq 'list')? '' : '<INPUT type="button" value="Clear list" style="font-weight:bold" onclick="document.getElementById(\'protListDIV\').style.display=\'none\'">';
	my $numColumns=($refLabelingInfo->{'AVG_MODE'}[0] eq 'AVG_ALL')? 7 : scalar (keys %dispStates) + 6;
	my %protRulerMeasures=(	'COPY_NB'			=> 'Copy number per cell',
							'CONCENTRATION'		=> 'Concentration [nM]',
							'MASS_PER_CELL'		=> 'Mass per Cell [pg]',
							'MASS_ABUNDANCE'	=> 'Mass Abundance (mass/total mass) [*10^-6]',
							'MOL_ABUNDANCE'		=> 'Mol. Abundance (molecules/total mol.) [*10^-6]',
							'COPY_RANK'			=> 'Copy Number Rank',
							'REL_COPY_RANK'		=> 'Relative Copy Number Rank'
						  );
	my $selProtRequired = ($action eq 'ajaxListProt') ? 'true' : 'false';
	if (%{$refQuantifValues}) {
		print qq
|<FORM name="protForm" method="POST">
<TABLE border=0 cellspacing=0 cellpadding=2>
<TR><TD colspan=$numColumns><DIV id="saveProtDIV" style="display:none;"></DIV></TD></TR>
<TR><TD colspan=$numColumns>$clearButtonStrg
<INPUT type="button" value="Check all" onclick="checkAllProteins(true)"/>
<INPUT type="button" value="Uncheck all" onclick="checkAllProteins(false)"/>
<INPUT type="button" id="saveFormBUTTON" value="Save proteins..." onclick="ajaxManageSaveProteins('getThemes','PROT')"$disabSave/>
|;
		print "<INPUT type=\"button\" value=\"Export data\" onclick=\"checkProtToExport($selProtRequired)\"/>\n";
	}
	else {
		print qq
|<TABLE border=0 cellspacing=0 cellpadding=2>
|;
	}

	print qq
|</TD></TR>
<TR bgcolor="$darkColor">
<TH class="rbBorder" align=left rowspan=2><DIV id="protCountDIV">$protString</DIV></TH><TH class="rbBorder" rowspan=2>Gene</TH>
|;
	my $numStates=($refLabelingInfo->{'AVG_MODE'}[0] eq 'AVG_ALL')? 1 : scalar keys %dispStates;
	my $quantifParamStrg=$protRulerMeasures{$dispMeasure};
	if (%{$refQuantifValues}) {
		print "<TH class=\"rbBorder\" colspan=$numStates>&nbsp$quantifParamStrg&nbsp;</TH>";
	}

	print qq
|<TH class="rbBorder" rowspan=2 nowrap>&nbsp;MW<SMALL> kDa</SMALL>&nbsp;</TH><TH class="bBorder" rowspan=2 width=700 nowrap>&nbsp;Description - Species&nbsp;</TH>
</TR>
<TR bgcolor="$darkColor">
|;
	my $firstSetUsed;
	if (%{$refQuantifValues}) {
		if ($refLabelingInfo->{'AVG_MODE'}[0] eq 'AVG_ALL') {
			$firstSetUsed=1;
			print "<TH class=\"rbBorder\" nowrap>Average on all states&nbsp;</TH>";
		} else {
			foreach my $statePos (sort{$a<=>$b} keys %dispStates) {
				$firstSetUsed=$statePos unless $firstSetUsed;
				print "<TH class=\"rbBorder\" nowrap>&nbsp;$refStateInfo->{$statePos}{NAME}&nbsp;</TH>";
			}
		}
	}

	my $bgColor=$lightColor;
	my $numDispItems=0;
	my %dispProteins;
	foreach my $modProtID (sort{&sortProteins($dispSort,$refProteinInfo,$refQuantifValues,$firstSetUsed)} keys %{$refProteinInfo}) { #%{$refQuantifValues} Changed to keys of refProteinInfo to display protein information even when quantif values are not available (typically, custom proteins list) (VL 22/08/19)
		my ($protID,$modResStrg)=($modProtID=~/^(\d+)(.*)/); # eg. 12345-Y40,T152 -> 12345, -Y40,T152
		$modResStrg='' unless $modResStrg;
		if ($refProteinInfo->{$protID}[0]) {  # display line only if protein has an identifier (VL 22/08/19)
			my $okPeptide=1;
			my $mass=(!$refProteinInfo->{$protID}[1])? '-' : $refProteinInfo->{$protID}[1];
			my $anaList = $refProteinInfo->{$protID}[-1];
			my $quantifDataStrg='';
			if (%{$refQuantifValues}) {
				if ($refLabelingInfo->{'AVG_MODE'}[0] eq 'AVG_ALL') {
					my $statePos=0;
					$quantifDataStrg.=($refQuantifValues->{$modProtID} && $refQuantifValues->{$modProtID}{$dispMeasure}{$statePos})? "<TH align=right>&nbsp;$refQuantifValues->{$modProtID}{$dispMeasure}{$statePos}&nbsp;</TH>" : '<TH align=right>-&nbsp;</TH>';
				} else {
					foreach my $statePos (sort{$a<=>$b} keys %dispStates) {
						$quantifDataStrg.=($refQuantifValues->{$modProtID} && $refQuantifValues->{$modProtID}{$dispMeasure}{$statePos})? "<TH align=right>&nbsp;$refQuantifValues->{$modProtID}{$dispMeasure}{$statePos}&nbsp;</TH>" : '<TH align=right>-&nbsp;</TH>';
					}
				}
			}
				
			$numDispItems++;
			$dispProteins{$protID}=1;
			my $class=($anaList)? 'TH' : 'TD'; # quantified MaxQuant proteins can be hidden in myProMS
			print "<TR bgcolor=\"$bgColor\" valign=top><TD class=\"$class\" nowrap>";
			if (%{$refQuantifValues}) {
				print "<INPUT type=\"checkbox\" name=\"chkProt\" value=\"$modProtID\"/>";
			}
			print "<A href=\"javascript:sequenceView($protID,'$anaList')\">$refProteinInfo->{$protID}[0]$modResStrg</A>&nbsp;</TD>\n";
			# Gene(s)
			if ($refProteinInfo->{$protID}[5] && scalar @{$refProteinInfo->{$protID}[5]} > 1) {
				print "<TH>&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-",join('<BR>&nbsp;&nbsp;-',@{$refProteinInfo->{$protID}[5]}[1..$#{$refProteinInfo->{$protID}[5]}]),"</B>')\" onmouseout=\"popout()\">",$refProteinInfo->{$protID}[5][0],"</A>&nbsp;</TH>";
			}
			elsif ($refProteinInfo->{$protID}[5] && $refProteinInfo->{$protID}[5][0]) {print '<TH>&nbsp;',$refProteinInfo->{$protID}[5][0],'&nbsp;</TH>';}
			else {print '<TH>-</TH>';} # no gene mapped
			# quantif data
			print $quantifDataStrg;
			# mw, des & species
			print "<TD class=\"TH\" align=right>$mass&nbsp;</TD><TD>$refProteinInfo->{$protID}[2] <U><I>$refProteinInfo->{$protID}[3]</I></U></TD>\n";
			print "</TR>\n";
			$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
		}
	}
	$protString=($isModifQuantif)? "&nbsp;$numDispItems/$numTotQuantifItems sites<BR>&nbsp;".(scalar keys %dispProteins)."/$numTotProteins proteins&nbsp;" : "&nbsp;$numDispItems/$numTotProteins proteins&nbsp;";
	print qq
|</TABLE>
</FORM>
<SCRIPT type="text/javascript"><!-- Not executed if called from AJAX -->
document.getElementById('protCountDIV').innerHTML='$protString';
</SCRIPT>
|;
	
}


sub exportProteinList { # Only for design or no-design labeled quantifs
	my ($labelType,$ratioType,$refLabelingInfo,$refStateInfo,$refQuantifValues,$refProteinInfo,$refDispModifSites,$pvalueCode)=@_; #,$dispStdDev
	my $invDispFoldChange=1/$dispFoldChange;
	my $pValueType=(!$refLabelingInfo->{'FDR_CONTROL'} || $refLabelingInfo->{'FDR_CONTROL'}[0]=~/FALSE|none/i)? 'p-value' : 'Adj. p-value';
	my $quantifSoftware=($refLabelingInfo->{'SOFTWARE'})? $refLabelingInfo->{'SOFTWARE'}[0] : 'myProMS'; $quantifSoftware='MaxQuant' if $quantifSoftware eq 'MQ';
	my $algoVersion=0; # myProMS only
	if ($quantifSoftware eq 'myProMS') {
		$algoVersion=($refLabelingInfo->{'SOFTWARE'} && $refLabelingInfo->{'SOFTWARE'}[1])? $refLabelingInfo->{'SOFTWARE'}[1] : ($ratioType eq 'Ratio')? 1 : 2;
	}
	my $numTotQuantifItems=scalar keys %{$refQuantifValues};
	my $numTotProteins=scalar keys %{$refProteinInfo};
	my $peptideTitle=($numPepCode eq 'DIST_PEP_USED')? 'Dist. P' : 'P';
	$peptideTitle.=(($refLabelingInfo->{'ALGO_TYPE'} && $refLabelingInfo->{'ALGO_TYPE'}[0] eq 'PEP_RATIO') || ($refLabelingInfo->{'QUANTIFICATION_METHOD'} && $refLabelingInfo->{'QUANTIFICATION_METHOD'}[0] eq 'Ratio'))? 'ept. sets' : 'eptides';

	####<Start printing>####
	my $worksheet2=$workbook->add_worksheet('Results');
	my $xlsRow=0;
	my $xlsCol=0;

	###<Headers>###
	#<Identifiers
	$worksheet2->set_column(0,0,30); # identifier col width
	$worksheet2->set_column(1,1,30); # Gene names col width
	$worksheet2->set_row(1,39); # Increase 2nd headers row height
	$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow+1,$xlsCol,'Gene & Synonyms',$itemFormat{'mergeRowHeader'});
	#$worksheet2->merge_range($xlsRow,$xlsCol,$xlsRow+1,$xlsCol,"Proteins",$itemFormat{'mergeRowHeader'});
	
	my $pepFilterComment='';
	unless ($quantifType eq 'PROT_RULER') {
		my ($pcRatioStrg,$pcStateStrg)=($dispPepType eq 'msmsPc')? ('% identified','Min. %') : ('Identified','Min.');
		$pepFilterComment=($dispPepType eq 'distinct')? 'Distinct p' : ($dispPepType =~ /^msms/ && $numPepScope eq 'ratio')? "$pcRatioStrg p" : ($dispPepType =~ /^msms/)? "$pcStateStrg ident. p" : 'P';
		$pepFilterComment.='eptides used';
		if ($ratioType eq 'None') {
			$pepFilterComment.=' in at least one state';
		}
		elsif ($numPepScope =~ /^replicate/ || $dispPepType =~ /^msms/) {
			$pepFilterComment.=' per replicate' if $numPepScope =~ /^replicate/;
			$pepFilterComment.=($numPepScope =~ /One$/)? ' in one state' : ($numPepScope =~ /Both$/)? ' in both states' : ($numPepScope =~ /Test$/)? ' in Test state' : ' in Reference state';
		}
		else {$pepFilterComment.=' for ratio';}
		$pepFilterComment.=' ≥ '.$dispNumPep;
		if ($numPepScope =~ /^replicate/) {
			$pepFilterComment.=" in at least $dispNumReplic replicate";
			$pepFilterComment.='s' if $dispNumReplic > 1;
		}
	}
	my $fullPepTypeFiltering=($dispPepType =~ /^msms/ || $numPepScope ne 'ratio')? 1 : 0; # simple dist for ratio is NOT included here!
	my %deltaColSpan;
	my $firstRatioUsed=0;
	if ($ratioType eq "None") {  # Abundance quantifications
		my $measStrg = ($dispMeasure eq 'MY_LFQ')? 'LFQ' : ucfirst(lc($dispMeasure)).'.';
		$measStrg =~ s/_/ /;
		$measStrg .= " ($dispAbsUnits{$dispMeasure})" if ($dispAbsUnits{$dispMeasure});
		foreach my $statePos (sort {$a <=> $b} keys %dispStates) {
			$firstRatioUsed = $statePos unless $firstRatioUsed;  # For compatibility with sortProteins and ratios
			my ($repNb, $bioReps, $state);
			if ($quantifType eq 'PROT_RULER') {
				$state = $statePos;
			} else {
				($repNb, $bioReps, $state) = split(',', $refLabelingInfo->{STATES}[$statePos - 1]);
				$state =~ s/#//g;
				if ($designID) {
					$state = $condToState{$state};
				}
			}
			$deltaColSpan{$statePos} = ($quantifType eq 'PROT_RULER')? 1 :
									   ($dispMeasure =~ /^(MOL|MASS)(_(PERCENT|CONC))?$/)? 5 :
									   3;  # Value [CI_low CI_up] Log10value [Peptides * 2] (-1 col)
			$deltaColSpan{$statePos}++ if $dispPepType =~/^msms/;
			$worksheet2->merge_range($xlsRow, ++$xlsCol, $xlsRow, $xlsCol+$deltaColSpan{$statePos}, decode_utf8($refStateInfo->{$state}{NAME}), $itemFormat{'mergeColHeader'});
			$worksheet2->write_string($xlsRow+1, $xlsCol, $measStrg, $itemFormat{'headerVCenter'});
			if ($quantifType eq 'PROT_ABUNDANCE' && $dispMeasure =~ /^(MOL|MASS)(_(PERCENT|CONC))?$/) {  # Need to display CIs
				$worksheet2->write_string($xlsRow+1, ++$xlsCol, 'Lower bound.', $itemFormat{'headerVCenter'});
				$worksheet2->write_comment($xlsRow+1, $xlsCol, "Lower boundary of abundance confidence interval: There is a ".($refLabelingInfo->{'CONFIDENCE_LEVEL'}[0]*100)."% certainty that true mean is included between lower and upper boundaries.");
				$worksheet2->write_string($xlsRow+1, ++$xlsCol, 'Upper bound.', $itemFormat{'headerVCenter'});
				$worksheet2->write_comment($xlsRow+1, $xlsCol, "Upper boundary of abundance confidence interval: There is a ".($refLabelingInfo->{'CONFIDENCE_LEVEL'}[0]*100)."% certainty that true mean is included between lower and upper boundaries.");
			}
			$worksheet2->write_string($xlsRow+1, ++$xlsCol, "Log10($measStrg)", $itemFormat{'headerVCenter'});
			unless ($quantifType eq 'PROT_RULER') {  # Don't consider peptides number for Proteomic Ruler (not recorded)
				$worksheet2->write_string($xlsRow+1,++$xlsCol,'Dist. pept. used',$itemFormat{'headerVCenter'});
				$worksheet2->write_string($xlsRow+1,++$xlsCol,'Pept. used',$itemFormat{'headerVCenter'});
				if ($dispPepType  =~ /^msms/) {
					my $pepTypeStrg = 'Valid. pep. used';
					$worksheet2->write_string($xlsRow+1,++$xlsCol,$pepTypeStrg,$itemFormat{'headerVCenter'});
				}
			}
		}
	}
	else {  # Ratio quantifications
		# if ($dispPepType ne 'all' || $numPepScope ne 'ratio') { # includes Distinct peptides!!!!
		# 	$pepTypeFiltering=1;
		# }
		foreach my $rPos0 (sort{abs($a)<=>abs($b)} keys %dispRatios) {
			my $rPos=abs($rPos0);
			$firstRatioUsed=$rPos0 unless $firstRatioUsed;
			my ($testStatePos,$refStatePos)=split(/\//,$refLabelingInfo->{RATIOS}[$rPos-1]);
			my $ratioTag='';
			if ($designID) {
				$testStatePos=~s/#//;
				$refStatePos=~s/#//;
				if ($refLabelingInfo->{RATIOS}[$rPos-1]=~/%/) {
					#$supRatios{$rPos0}=1;
					$ratioTag='°';
					$testStatePos=~s/%\d+//;
					$refStatePos=~s/%\d+//;
				}
				$testStatePos=$condToState{$testStatePos};
				$refStatePos=$condToState{$refStatePos};
			}
			if ($rPos0 < 0) { # reverse ratio
				my $tempPos=$testStatePos;
				$testStatePos=$refStatePos;
				$refStatePos=$tempPos;
			}
			$deltaColSpan{$rPos0}=($quantifType eq 'TNPQ')? 4 : ($quantifSoftware eq 'myProMS' && $fullPepTypeFiltering)? 7 : 5; # 5 before
			$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow,$xlsCol+$deltaColSpan{$rPos0},decode_utf8($refStateInfo->{$testStatePos}{'NAME'}.$ratioTag.'/'.$refStateInfo->{$refStatePos}{'NAME'}.$ratioTag),$itemFormat{'mergeColHeader'});
			$worksheet2->write_string($xlsRow+1,$xlsCol,'Ratio',$itemFormat{'header'});
			$worksheet2->write_comment($xlsRow+1,$xlsCol,'Ratio of 1000 or 0.001 means protein was not detected in one condition.');
			$worksheet2->write_string($xlsRow+1,++$xlsCol,'Log2(Ratio)',$itemFormat{'header'});
			$worksheet2->write_comment($xlsRow+1,$xlsCol,'Log2(Ratio) of 1000 or -1000 means protein was not detected in one condition.');
			#if ($ratioType eq 'Ratio') {
			#	$worksheet2->write_string($xlsRow+1,++$xlsCol,'Lower bound.',$itemFormat{'header'});
			#	$worksheet2->write_comment($xlsRow+1,$xlsCol,"Lower boundary of ratio confidence interval: There a ".($refLabelingInfo->{'CONFIDENCE_LEVEL'}[0]*100)."% certainty that true ratio is included between lower and upper boundaries.");
			#	$worksheet2->write_string($xlsRow+1,++$xlsCol,'Upper bound.',$itemFormat{'header'});
			#	$worksheet2->write_comment($xlsRow+1,$xlsCol,"Upper boundary of ratio confidence interval: There a ".($refLabelingInfo->{'CONFIDENCE_LEVEL'}[0]*100)."% certainty that true ratio is included between lower and upper boundaries.");
			#}
			$worksheet2->write_string($xlsRow+1,++$xlsCol,$pValueType,$itemFormat{'header'});
			$worksheet2->write_string($xlsRow+1,++$xlsCol,'CV %',$itemFormat{'header'}) if $quantifSoftware eq 'myProMS';
			if ($ratioType=~/S\w+Ratio/) {
				$worksheet2->write_string($xlsRow+1,++$xlsCol,'Dist. pept. used',$itemFormat{'header'});
				$worksheet2->write_string($xlsRow+1,++$xlsCol,'Pept. used',$itemFormat{'header'});
				if ($fullPepTypeFiltering) {
					if ($numPepScope =~ /^replicate(.+)/) {
						my $scope=$1;
						my $distStrg=($dispPepType eq 'all')? 'pep.' : 'dist. pep.';
						$worksheet2->write_string($xlsRow+1,++$xlsCol,"Pep. filter",$itemFormat{'header'});
						$worksheet2->write_string($xlsRow+1,++$xlsCol,"Distrib. $distStrg in replic.",$itemFormat{'header'});
						# my ($pepTypeStrg1,$pepTypeStrg2)=($dispPepType eq 'distinct')? ('Min. dist. ','Dist. p') : ('Min.','P');
						# $pepTypeStrg1.=' peptides';
						# $pepTypeStrg2.='eptide';
						# $worksheet2->write_string($xlsRow+1,++$xlsCol,"$pepTypeStrg1/replicate",$itemFormat{'header'});
						# $worksheet2->write_string($xlsRow+1,++$xlsCol,"$pepTypeStrg2 distrib. in replicates",$itemFormat{'header'});
					}
					elsif ($dispPepType =~ /^msms/) {
						# my $pepTypeStrg=($numPepScope eq 'ratio')? 'Valid. pep. used' : ($numPepScope eq 'msmsBoth')? 'Min. valid. pep. both st./all pep. used' : ($numPepScope eq 'msmsOne')? 'Min. valid. pep. best st.' : ($numPepScope eq 'msmsTest')? 'Min. valid. pep. test st.' : 'Min. valid. pep. ref. st.';
						$worksheet2->write_string($xlsRow+1,++$xlsCol,"Pep. filter",$itemFormat{'header'});
						$worksheet2->write_string($xlsRow+1,++$xlsCol,"Distrib. valid. pep. in states",$itemFormat{'header'});
					}
					$worksheet2->write_comment($xlsRow+1,$xlsCol-1,decode_utf8($pepFilterComment));
				}
			}
		}
	}
	#<Dataset-wide peptides
	if ($ratioType eq 'Ratio') {
		$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow,$xlsCol+2,$peptideTitle,$itemFormat{'mergeColHeader'});
		$worksheet2->write_string($xlsRow+1,$xlsCol,'Distinct used',$itemFormat{'header'});
		$worksheet2->write_string($xlsRow+1,++$xlsCol,'Used',$itemFormat{'header'});
		$worksheet2->write_string($xlsRow+1,++$xlsCol,'Total',$itemFormat{'header'});
	}
	elsif ($ratioType ne 'None') {
		$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow+1,$xlsCol,"Total peptides in set",$itemFormat{'mergeRowHeader'});
	}
	elsif ($quantifType ne 'PROT_RULER') {
		$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow+1,$xlsCol,"All peptides",$itemFormat{'mergeRowHeader'});
		$worksheet2->set_column($xlsCol,$xlsCol,20); # increase col width
	}
	#<MW, description & species
	$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow+1,$xlsCol,'MW (kDa)',$itemFormat{'mergeRowHeader'});
	$worksheet2->set_column(++$xlsCol,$xlsCol,80); # col length
	$worksheet2->merge_range($xlsRow,$xlsCol,$xlsRow+1,$xlsCol,'Description',$itemFormat{'mergeRowHeader'});
	$worksheet2->set_column(++$xlsCol,$xlsCol,30); # col length
	$worksheet2->merge_range($xlsRow,$xlsCol,$xlsRow+1,$xlsCol,'Species',$itemFormat{'mergeRowHeader'});

	###<Looping through proteins>###
	$xlsRow++;
	my $numDispItems=0;
	my %dispProteins;
	my %filterProteins;
	if ($filterProtIDs) {
		foreach my $modProtID (split(',',$filterProtIDs)) {$filterProteins{$modProtID}=1;}
	}
	# my ($pepCode4Filter,$pepFilterValue)=($dispPepType eq 'msmsPc')? ('FRAC_TRUE_USED',$dispNumPep/100) : ($numPepCode,$dispNumPep); # numPepCode is GLOBAL
	my ($pepCode4Filter,$pepFilterValue)=($dispPepType eq 'msmsPc')? ('FRAC_TRUE_USED',$dispNumPep/100) : ($ratioType ne 'None' && ($numPepScope =~ /^replicate/ || $dispPepType =~ /^msms/))? ($numPepCode.'_SUBSET',$dispNumPep) : ($numPepCode,$dispNumPep); # numPepCode is GLOBAL
	my $pepCode4Display;
	if ($ratioType eq 'None') {
		$pepCode4Display=($quantifType eq 'PROT_ABUNDANCE')? 'NUM_PEP_USED' : ($quantifType eq 'MQ')? 'PEPTIDES' : $numPepCode;
	}
	else { # ratio
		$pepCode4Display=($dispPepType eq 'msmsPc')? $numPepCode.'_SUBSET' : $pepCode4Filter; # msmsPc: displays counts instead of fraction
	}
	my $sortItem=($dispSort eq 'peptide')? $pepCode4Filter : $dispSort;
	my $log10=log(10);
	foreach my $modProtID (sort{&sortProteins($sortItem,$refProteinInfo,$refQuantifValues,$firstRatioUsed)} keys %{$refQuantifValues}) { #%{$refProteinInfo}
		#my $qModProtID=quotemeta($modProtID);
		#next if ($filterProtIDs && $filterProtIDs !~ /(^|,)$qModProtID(,|$)/);
		next if ($filterProtIDs && !$filterProteins{$modProtID});
		my ($protID,$modResStrg)=($modProtID=~/^(\d+)(.*)/); # eg. 12345-Y40,T152 -> 12345, -Y40,T152
		$modResStrg='' unless $modResStrg;
		$xlsCol=0;
		my $refUsedRatioPos=0;
		
		if (!$filterProtIDs) {
			#<Peptide filter
			my $okPeptide=0;
			if ($ratioType eq 'None') {
				foreach my $statePos (keys %dispStates) {
					next if ($quantifType ne 'PROT_RULER' && (
							!$refQuantifValues->{$modProtID}{$pepCode4Filter}{$statePos} ||
							$refQuantifValues->{$modProtID}{$pepCode4Filter}{$statePos} < $pepFilterValue
					));
					$okPeptide = 1; # at least 1 state must pass the filter (except prot ruler, peptides not recorded)
					$refUsedRatioPos = $statePos;  # Consistency with ratios
					last;
				}
			}
			else {
				foreach my $ratioPos (keys %dispRatios) {
				next if (!$refQuantifValues->{$modProtID}{$pepCode4Filter}{$ratioPos} || $refQuantifValues->{$modProtID}{$pepCode4Filter}{$ratioPos} < $pepFilterValue); # may not be defined for +/-inf ratios
					$okPeptide=1; # at least 1 ratio must pass the filter
					$refUsedRatioPos=$ratioPos;
					last;
				}
			}
			next unless $okPeptide;

			#<Ratio/p-value/Std dev filters
			my $okFilters = 0;
			if ($ratioType eq 'None') {
				$okFilters = 1;
			}
			else {
				foreach my $ratioPos (keys %dispRatios) {
					unless ($okFilters) { # no ratioPos has passed filters
						my $absRatio=0;
						if (defined $refQuantifValues->{$modProtID}{RATIO}{$ratioPos}) {
							$absRatio=($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} >= 1)? $refQuantifValues->{$modProtID}{RATIO}{$ratioPos} : 1/$refQuantifValues->{$modProtID}{RATIO}{$ratioPos};
							if ($foldChangeType eq 'abs') {
								#$okFilters=(($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} < 1 && $refQuantifValues->{$modProtID}{RATIO}{$ratioPos} <= $invDispFoldChange) || ($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} >= 1 && $refQuantifValues->{$modProtID}{RATIO}{$ratioPos} >= $dispFoldChange))? 1 : 0;
								$okFilters=($absRatio >= $dispFoldChange)? 1 : 0;
							}
							elsif ($foldChangeType eq 'up') {
								$okFilters=($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} >= 1 && $refQuantifValues->{$modProtID}{RATIO}{$ratioPos} >= $dispFoldChange)? 1 : 0;
							}
							else { # down
								$okFilters=($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} <= 1 && $refQuantifValues->{$modProtID}{RATIO}{$ratioPos} <= $invDispFoldChange)? 1 : 0;
							}
						}
						$okFilters=($okFilters && ($dispPvalue>=1 || ($allowInfRatio && $absRatio==1000) || (defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} && $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} <= $dispPvalue)))? 1 : 0;
						if ($okFilters && $dispCV && defined $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}) { # undefined always passes filter!
							if ($ratioType eq 'Ratio') { # v1
								$okFilters=0 if abs(log($refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos})/log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos})) > $dispCVFrac;
							}
							else { # no CV filter for TNPQ or MaxQuant
								$okFilters=0 if $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}/abs(log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos})/$log2) > $dispCVFrac;
							}
						}
					}
				}
			}
			next unless $okFilters;
		}
		
		$numDispItems++;
		$dispProteins{$protID}=1;
		# Identifier
		my $dispModIdentifier=($refDispModifSites->{$modProtID})? $refProteinInfo->{$protID}[0].'-'.$refDispModifSites->{$modProtID} : $refProteinInfo->{$protID}[0];
		$worksheet2->write_string(++$xlsRow,$xlsCol,$dispModIdentifier,$itemFormat{'text'});
		# Gene(s)
		$worksheet2->write_string($xlsRow,++$xlsCol,join(',',@{$refProteinInfo->{$protID}[5]}),$itemFormat{'text'});

		if ($ratioType eq 'None') {  # Abundance quantifications
			foreach my $statePos (sort {$a <=> $b} keys %dispStates) {
				$refUsedRatioPos=$statePos unless $refUsedRatioPos;  # For consistency with ratios
				if (!$refQuantifValues->{$modProtID}{$dispMeasure}{$statePos}) {
					foreach my $i (0..$deltaColSpan{$statePos}) {
						$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
					}
					next;
				}
				# Abundance measure
				$worksheet2->write_number($xlsRow, ++$xlsCol, $refQuantifValues->{$modProtID}{$dispMeasure}{$statePos}, $itemFormat{'number'});
				# Log10 Abundance measure
				$worksheet2->write_number($xlsRow, ++$xlsCol, (log($refQuantifValues->{$modProtID}{$dispMeasure}{$statePos}) / $log10), $itemFormat{'number'});
				if ($dispMeasure =~ /^(MOL|MASS)(_(PERCENT|CONC))?$/) {  # Display Confidence intervals
					if ($refQuantifValues->{$modProtID}{$absCIsCodes{$dispMeasure}[0]}{$statePos} && $refQuantifValues->{$modProtID}{$absCIsCodes{$dispMeasure}[1]}{$statePos}) {
						$worksheet2->write_number($xlsRow, ++$xlsCol, $refQuantifValues->{$modProtID}{$absCIsCodes{$dispMeasure}[0]}{$statePos}, $itemFormat{'number'});
						$worksheet2->write_number($xlsRow, ++$xlsCol, $refQuantifValues->{$modProtID}{$absCIsCodes{$dispMeasure}[1]}{$statePos}, $itemFormat{'number'});
					} else {
						$worksheet2->write_blank($xlsRow, ++$xlsCol, $itemFormat{'number'});
						$worksheet2->write_blank($xlsRow, ++$xlsCol, $itemFormat{'number'});
					}
				}
				# Peptides in state
				unless ($quantifType eq 'PROT_RULER') {  # Don't consider peptides number for Proteomic Ruler (not recorded)
					$worksheet2->write_number($xlsRow, ++$xlsCol, $refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$statePos}, $itemFormat{'number'});
					$worksheet2->write_number($xlsRow, ++$xlsCol, $refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$statePos}, $itemFormat{'number'});		
					if ($dispPepType =~ /^msms/) {
						if ($refQuantifValues->{$modProtID}{$numPepCode}{$statePos}) {
							# my $numGhostPep = $refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$statePos} - $refQuantifValues->{$modProtID}{$numPepCode}{$statePos};
							# $numGhostPep = ($numGhostPep)? '+' . $numGhostPep : '';
							$worksheet2->write_number($xlsRow, ++$xlsCol, $refQuantifValues->{$modProtID}{$numPepCode}{$statePos}, $itemFormat{'number'}); # . $numGhostPep
						} else {
							$worksheet2->write_blank($xlsRow, ++$xlsCol, $itemFormat{'number'});
						}
					}
					# else {
					# 	# DIST_PEP_USED
					# 	if ($refQuantifValues->{$modProtID}{'DIST_PEP_USED'} && $refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$statePos}) {
					# 		$worksheet2->write_number($xlsRow, ++$xlsCol, $refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$statePos}, $itemFormat{'number'});
					# 	} else {
					# 		$worksheet2->write_blank($xlsRow, ++$xlsCol, $itemFormat{'number'});
					# 	}
					# 	# NUM_PEP_USED
					# 	if ($refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$statePos} && $refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$statePos}) {
					# 		$worksheet2->write_number($xlsRow, ++$xlsCol, $refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$statePos}, $itemFormat{'number'});
					# 	} else {
					# 		$worksheet2->write_blank($xlsRow, ++$xlsCol, $itemFormat{'number'});
					# 	}
					# }
				}
			}
		}
		else { # Ratio quantifications
			foreach my $ratioPos (sort{$a<=>$b} keys %dispRatios) {
				$refUsedRatioPos=$ratioPos unless $refUsedRatioPos;
				if (!$refQuantifValues->{$modProtID}{RATIO}{$ratioPos}) { # !$refQuantifValues->{$protID}{'RATIO'} ||
					foreach my $i (0..$deltaColSpan{$ratioPos}) {
						#$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
						$xlsCol++;
					}
					next;
				}
				# Ratio
				$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{RATIO}{$ratioPos},$itemFormat{'number'});
				# log2 ratio
				my $log2Ratio=($refQuantifValues->{$modProtID}{RATIO}{$ratioPos}==1000)? 1000 : ($refQuantifValues->{$modProtID}{RATIO}{$ratioPos}==0.001)? -1000 : log($refQuantifValues->{$modProtID}{RATIO}{$ratioPos})/$log2;
				$worksheet2->write_number($xlsRow,++$xlsCol,$log2Ratio,$itemFormat{'number'});
				#if ($ratioType eq 'Ratio') {
				#	# low conf limit
				#	if ($refQuantifValues->{$protID}{'CONF_LOW'} && defined $refQuantifValues->{$modProtID}{'CONF_LOW'}{$ratioPos}) {
				#		$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'CONF_LOW'}{$ratioPos},$itemFormat{'number'});
				#	}
				#	else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
				#	# upper conf limit
				#	if ($refQuantifValues->{$protID}{'CONF_UP'} && defined $refQuantifValues->{$modProtID}{'CONF_UP'}{$ratioPos}) {
				#		$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'CONF_UP'}{$ratioPos},$itemFormat{'number'});
				#	}
				#	else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
				#}
				# p-value
				if (defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}) {
					$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos},$itemFormat{'number'});
				}
				else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}

				if ($quantifSoftware eq 'myProMS') {
					# Std dev -> CV
					if ($refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}) {
						my $cv=($ratioType eq 'Ratio')? abs(log($refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos})/log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos})) : $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}/abs(log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos})/$log2);
						$worksheet2->write_number($xlsRow,++$xlsCol,100*$cv,$itemFormat{'number'});
					}
					else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
				}
				# Simple/SuperRatio peptides
				if ($ratioType=~/S\w+Ratio/) {
					# DIST_PEP_USED
					if ($refQuantifValues->{$modProtID}{'DIST_PEP_USED'} && $refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$ratioPos}) {
						$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$ratioPos},$itemFormat{'number'});
					}
					else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
					# NUM_PEP_USED
					if ($refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$ratioPos} && $refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$ratioPos}) {
						$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$ratioPos},$itemFormat{'number'});
					}
					else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}

					if ($fullPepTypeFiltering) { # 2 extra columns
						$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{$pepCode4Display}{$ratioPos},$itemFormat{'number'});
						if ($numPepScope eq 'ratio') { # => $dispPepType=msms(Num|Pc)
							$worksheet2->write_number($xlsRow,++$xlsCol,1*(sprintf "%.1f",$refQuantifValues->{$modProtID}{'FRAC_TRUE_USED'}{$ratioPos}*100),$itemFormat{'number'}) if $dispPepType eq 'msmsPc';
						}
						else { # state-level peptide filtering
							my %pepRepStrg=(0=>'',1=>'');
							foreach my $s (0,1) {
								foreach my $i (0..$#{$dispRatios{$ratioPos}[$s]}) { # only [0] for msms
									#my $targetPos=$dispRatios{$ratioPos}[$s][$i]*$tgtPosSwap;
									my $targetPos=$dispRatios{$ratioPos}[$s][$i];
									$targetPos*=-1 if (($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos}==0.001 || $refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos}==1000) && $refQuantifValues->{$modProtID}{$numPepCode} && defined $refQuantifValues->{$modProtID}{$numPepCode}{-$targetPos}); # if Inf ratio: take negative tgt pos if exists ()
									$pepRepStrg{$s}.=', ' if $pepRepStrg{$s} ne '';
									if ($dispPepType eq 'msmsPc') {
										my ($numPep,$fracPep)=($refQuantifValues->{$modProtID}{$numPepCode} && $refQuantifValues->{$modProtID}{$numPepCode}{$targetPos})? ($refQuantifValues->{$modProtID}{$numPepCode}{$targetPos},$refQuantifValues->{$modProtID}{'FRAC_TRUE_USED'}{$targetPos}) : (0,0);
										my $pcPep=1*(sprintf "%.1f",$fracPep*100);
										my $pcPepStrg=($pcPep)? ' ('.$pcPep.'%)' : '';
										$pepRepStrg{$s}.="$numPep$pcPepStrg";
									}
									else {
										my $numPep=($refQuantifValues->{$modProtID}{$numPepCode} && $refQuantifValues->{$modProtID}{$numPepCode}{$targetPos})? $refQuantifValues->{$modProtID}{$numPepCode}{$targetPos} : 0;
										$pepRepStrg{$s}.=$numPep;
									}
								}
							}
							$worksheet2->write_string($xlsRow,++$xlsCol,$pepRepStrg{1}.'/'.$pepRepStrg{0},$itemFormat{'text'});
						}
					}
					# if ($numPepScope =~ /^replicate/) {
					# 	if ($refQuantifValues->{$modProtID}{$numPepCode} && $refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos}) {
					# 		$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos},$itemFormat{'number'});
					# 		my %pepRepStrg=(0=>'',1=>'');
					# 		foreach my $s (0,1) {
					# 			foreach my $i (0..$#{$dispRatios{$ratioPos}[$s]}) {
					# 				my $replicPos=$dispRatios{$ratioPos}[$s][$i];
					# 				$replicPos*=-1 if (($refQuantifValues->{$modProtID}{RATIO}{$ratioPos}==0.001 || $refQuantifValues->{$modProtID}{RATIO}{$ratioPos}==1000) && $refQuantifValues->{$modProtID}{$numPepCode} && $refQuantifValues->{$modProtID}{$numPepCode}{-$replicPos}); # if Inf ratio: take negative tgt pos if exists
					# 				my $numPep=($refQuantifValues->{$modProtID}{$numPepCode} && $refQuantifValues->{$modProtID}{$numPepCode}{$replicPos})? $refQuantifValues->{$modProtID}{$numPepCode}{$replicPos} : 0;
					# 				$pepRepStrg{$s}.=',' if $pepRepStrg{$s} ne '';
					# 				$pepRepStrg{$s}.=$numPep;
					# 			}
					# 		}
					# 		$worksheet2->write_string($xlsRow,++$xlsCol,$pepRepStrg{1}.'/'.$pepRepStrg{0},$itemFormat{'text'});
					# 	}
					# 	else {
					# 		$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
					# 		$worksheet2->write_string($xlsRow,++$xlsCol,'',$itemFormat{'text'});
					# 	}
					# }
					# elsif ($dispPepType =~ /^msms/) {
					# 	if ($refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$ratioPos}) {
					# 		my %identSt;
					# 		foreach my $s (0,1) {
					# 			my $stPos=$dispRatios{$ratioPos}[$s][0];
					# 			$stPos*=-1 if (($refQuantifValues->{$modProtID}{RATIO}{$ratioPos}==0.001 || $refQuantifValues->{$modProtID}{RATIO}{$ratioPos}==1000) && $refQuantifValues->{$modProtID}{$numPepCode} && $refQuantifValues->{$modProtID}{$numPepCode}{-$stPos}); # if Inf ratio: take negative tgt pos if exists
					# 			$identSt{$s}=($refQuantifValues->{$modProtID}{$numPepCode} && $refQuantifValues->{$modProtID}{$numPepCode}{$stPos})? $refQuantifValues->{$modProtID}{$numPepCode}{$stPos} : 0;
					# 		}
					# 		my $numGhostPep=$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$ratioPos}-$refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos};
					# 		$numGhostPep=($numGhostPep)? '+'.$numGhostPep : '';
					# 		$worksheet2->write_string($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos}.$numGhostPep,$itemFormat{'number'});
					# 		$worksheet2->write_string($xlsRow,++$xlsCol,$identSt{1}.'/'.$identSt{0},$itemFormat{'text'});
					# 	}
					# 	else {
					# 		$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
					# 		$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
					# 	}
					# }
					# else {
					# }
				}
			}
		}
		# Peptides
		unless ($quantifType eq 'PROT_RULER') {
			if ($refQuantifValues->{$modProtID}) {
				if ($ratioType eq 'Ratio') {
					# DIST_PEP_USED
					if ($refQuantifValues->{$modProtID}{'DIST_PEP_USED'} && $refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$refUsedRatioPos}) {
						$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$refUsedRatioPos},$itemFormat{'number'});
					}
					else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
					# NUM_PEP_USED
					if ($refQuantifValues->{$modProtID}{'NUM_PEP_USED'} && $refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$refUsedRatioPos}) {
						$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$refUsedRatioPos},$itemFormat{'number'});
					}
					else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
				}
				# NUM_PEP_TOTAL
				if ($refQuantifValues->{$modProtID}{'NUM_PEP_TOTAL'} && $refQuantifValues->{$modProtID}{'NUM_PEP_TOTAL'}{$refUsedRatioPos}) {
					$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'NUM_PEP_TOTAL'}{$refUsedRatioPos},$itemFormat{'number'});
				}
				else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
			}
			else {
				if ($ratioType eq 'Ratio') {
					$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
					$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
				}
				$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
			}
		}
		# MW, desc, species
		my $mass=(!$refProteinInfo->{$protID}[1])? 0 : $refProteinInfo->{$protID}[1];
		$worksheet2->write_number($xlsRow,++$xlsCol,$mass,$itemFormat{'number1d'});
		$worksheet2->write_string($xlsRow,++$xlsCol,$refProteinInfo->{$protID}[2],$itemFormat{'textWrap'});
		$worksheet2->write_string($xlsRow,++$xlsCol,$refProteinInfo->{$protID}[3],$itemFormat{'text'});
	}
	#<Identifier header (written last to get remaining number of proteins)
	#my $numTotQuantifItems=scalar keys %{$refQuantifValues};
	my $protString=($isModifQuantif)? "$numDispItems/$numTotQuantifItems sites\n".(scalar keys %dispProteins)."/$numTotProteins proteins" : "$numDispItems/$numTotProteins proteins";
	$worksheet2->merge_range(0,0,1,0,$protString,$itemFormat{'mergeRowHeader'});
}

sub ajaxListSelectedProteins { # GLOBALS: $ratioType,%dispStates,%dispRatios (all still undef!)
# print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
	my ($refList,$refTgPosList)=@_; # optional argument (Proteomic ruler)
	my (@selProteins,@targetPosList);
	if ($refList) {  # Sub called from another sub (ajaxDisplayCustomProteins) and parameters don't come from get/post request
		@selProteins=@{$refList};
		@targetPosList=@{$refTgPosList};
	}
	else {
		@selProteins=split(',',param('selProt'));
		#my ($minFoldChange,$maxFoldChange,$maxPvalue)=split(',',param('thresholds'));
		@targetPosList=split(',',param('dispTargetPosAjax'));
		#foreach my $tgtPos (@targetPosList) {$dispRatios{$tgtPos}=1;} # %dispRatios is global but not defined yet
	}

	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;
	my $projectID=&promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION');

	my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
	$projectAccess=${$userInfo[2]}{$projectID}; # global


	###>Fetching protein quantification data
	my %quantifParamID2code;
	my $sthQP2=$dbh->prepare("SELECT M.CODE,ID_QUANTIF_PARAMETER,P.NAME,P.CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD M,QUANTIFICATION_PARAMETER P WHERE Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND M.ID_QUANTIFICATION_METHOD=P.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$selQuantifID");
	$sthQP2->execute;
	while (my ($qType,$paramID,$paramName,$paramCode)=$sthQP2->fetchrow_array) {
		$quantifType=$qType;
		@{$quantifParamInfo{$paramCode}}=($paramID,$paramName);
		$quantifParamID2code{$paramID}=$paramCode;
	}
	$sthQP2->finish;

	####>Protein ratio (from peptides)<####
	my %labelingInfo;
	# WARNING: $desingID is GLOBAL!!!!!!!!
	#(my $quantifAnnot,$designID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_DESIGN FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	(my $quantifAnnot,$designID,my $quantifModID,my $multiModifStrg)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_DESIGN,Q.ID_MODIFICATION,GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ',')
																	 FROM QUANTIFICATION Q
																	 LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
																	 WHERE Q.ID_QUANTIFICATION=$selQuantifID GROUP BY Q.ID_QUANTIFICATION");
	#my ($isModifQuantif,$isMultiModifQuantif)=(0,0);
	#my %quantifModifInfo;
	#if ($quantifModID || $multiModifStrg) {
	#	$isModifQuantif=1;
	#	my $modifStrg=$quantifModID || $multiModifStrg;
	#	$quantifModID=$modifStrg if $modifStrg !~ /,/; # <- new single PTM quantif (also uses MULTIMODIF_Q table!)
	#	my $sthMod=$dbh->prepare("SELECT ID_MODIFICATION,DISPLAY_CODE,DISPLAY_COLOR FROM MODIFICATION WHERE ID_MODIFICATION IN ($modifStrg) ORDER BY FIELD(ID_MODIFICATION,$modifStrg)");
	#	my $modifRank=0;
	#	$sthMod->execute;
	#	while (my ($modID,$displayCode,$displayColor)=$sthMod->fetchrow_array) {
	#		$modifRank++;
	#		$quantifModifInfo{ID2RANK}{$modID}=$modifRank;
	#		#$quantifModifInfo{RANK2ID}{$modifRank}=$modID;
	#		$quantifModifInfo{DISPLAY}{$modID}=[$displayCode,$displayColor];
	#	}
	#	$sthMod->finish;
	#	$isMultiModifQuantif=1 if $modifRank > 1;
	#}
	my %quantifModifInfo;
	my ($isModifQuantif,$isMultiModifQuantif)=&promsQuantif::getQuantifModificationInfo($dbh,$quantifModID,$multiModifStrg,\%quantifModifInfo);
	
	my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
	my ($labelType)=($labelStrg)? ($labelStrg=~/LABEL=(.+)/) : ('FREE');
	$labelType=uc($labelType);
	$labelingInfo{'LABELTYPE'}=$labelType;
	foreach my $infoStrg (@labelInfo) {
		my ($setting,$valueStrg)=split('=',$infoStrg);
		$valueStrg=~s/#//g;
		@{$labelingInfo{$setting}}=split(';',$valueStrg);
	}
	$ratioType=($labelingInfo{'RATIO_TYPE'})? $labelingInfo{'RATIO_TYPE'}[0] : (!$labelingInfo{'RATIOS'})? 'None' : 'Ratio'; # GLOBAL!
	my $algoVersion=($ratioType eq 'Ratio')? 1 : (!$labelingInfo{'SOFTWARE'})? 2 : ($labelingInfo{'SOFTWARE'}[0] eq 'myProMS')? $labelingInfo{'SOFTWARE'}[1] : 0; # myProMS version ONLY!
	my $isPepIntensity=($labelingInfo{'MEAN_STATE'} && $labelingInfo{'MEAN_STATE'}[0]==1)? 1 : ($algoVersion==2 && $labelType eq 'FREE')? 1 : ($algoVersion>=3 && $labelingInfo{'ALGO_TYPE'}[0] eq 'PEP_INTENSITY')? 1 : 0; # 0 for v1 Ratio, v2 Label & v3 PEP_RATIO

	my $numStates;
	if ($quantifType eq 'PROT_RULER') {
		$protRulerCase=1;
		$numStates=scalar @{$labelingInfo{'PARENT_Q'}};
	}
	else {
		@{$labelingInfo{'IS_MODIF'}}=($isModifQuantif); # add modif status to be passed to &listProteinRatios
		$ratioType=($labelingInfo{'RATIO_TYPE'})? $labelingInfo{'RATIO_TYPE'}[0] : 'Ratio';
		$numStates=scalar @{$labelingInfo{'STATES'}};
		%ratioCode=&generateRatioCodes($labelType,$ratioType,$numStates,$isPepIntensity,$algoVersion) if $ratioType ne 'None';
	}
	
	my %stateInfo;
	my $statePos=0;
	
	if ($labelingInfo{PARENT_Q}) {
		my $sthParentQuantifName=$dbh->prepare("SELECT NAME FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
		my $sthParentStateName = $dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
		my (%parentQuantifName,%stateName);
		foreach my $stateData (@{$labelingInfo{PARENT_Q}}) {
			$statePos++;
			my ($parentQuantif,$targetPos)=split(':',$stateData);
			my ($parentQuantifID, $parentStateID, $parentPos) = split(',', $parentQuantif);
			unless ($parentQuantifName{$parentQuantifID}) {
				$sthParentQuantifName->execute($parentQuantifID);
				($parentQuantifName{$parentQuantifID}) = $sthParentQuantifName->fetchrow_array;
			}
			unless ($stateName{$parentStateID}) {
				$sthParentStateName->execute($parentStateID);
				($stateName{$parentStateID}) = $sthParentStateName->fetchrow_array;
			}
			$stateInfo{$statePos}{'NAME'} = "$parentQuantifName{$parentQuantifID} - $stateName{$parentStateID}";
		}
		$sthParentQuantifName->finish;
		$sthParentStateName->finish;
	}
	else {
		my $sthgetExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
		foreach my $stateData (@{$labelingInfo{STATES}}) {
			$statePos++;
			my ($numBioRep,$stateContent,$repPosStrg)=split(',',$stateData);
			if ($designID){
				my $expCondID=$repPosStrg;
				#$expCondID=~ s/#//;
				$condToState{$expCondID}=$statePos;
				$sthgetExpCondName->execute($expCondID);
				($stateInfo{$statePos}{'NAME'})=$sthgetExpCondName->fetchrow_array;
				if ($numPepScope =~ /^(replicate|msms)/) {
					my ($numBioRep,$allTechRepCount)=(0,0);
					foreach my $bioReplicate (split(/\./,$stateContent)) {
						$numBioRep++;
						$allTechRepCount+=scalar (split(/&/,$bioReplicate));
					}
					#$stateInfo{$statePos}{'NUM_REPLIC'}=($numBioRep > 1)? $numBioRep : $allTechRepCount; # needed for targetPos selection for numPep filtering
					$stateInfo{$statePos}{'NUM_REPLIC'}=( ($labelingInfo{'RESIDUAL_VAR'} && $labelingInfo{'RESIDUAL_VAR'}[0] eq 'biological') || (!$labelingInfo{'RESIDUAL_VAR'} && $numBioRep > 1) )? $numBioRep : $allTechRepCount; # needed for targetPos selection for numPep filtering
				}
			}
			else {
				$stateInfo{$statePos}{'NAME'}=$stateContent;
				$stateInfo{$statePos}{'NAME'}=~s/\./\+/g;
			}
		}
		$sthgetExpCondName->finish;
	}

	$dbh->disconnect;

	my (%quantifValues,%proteinInfo);
	my (%dispModifSites);

	my (%fileFormattedList,%queryModifSites);
	foreach my $modProtID (@selProteins) {
		my ($protID,$modResStrg)=($modProtID=~/^(\d+)-?(.*)/); $modResStrg='' unless $modResStrg;
		if ($modResStrg) {
			$queryModifSites{$modProtID}=&promsQuantif::siteCode2QuantifFileCode($modResStrg,$isMultiModifQuantif,$quantifModifInfo{ID2RANK});
			#$dispModifSites{$modProtID}=&promsQuantif::displayModificationSites($modResStrg,$quantifModifInfo{DISPLAY},'html'); <- done by fetchDesignQuantificationData
			$fileFormattedList{$protID.'-'.$queryModifSites{$modProtID}}=1;
		}
		else {$fileFormattedList{$protID}=1;}
	}
	if ($ratioType eq 'None') {
		foreach my $statePos (@targetPosList) {$dispStates{$statePos}=1;} # %dispStates is global but not defined yet
	}
	else {
		foreach my $ratioPos (@targetPosList) {$dispRatios{$ratioPos}=[];} # %dispRatios is global but not defined yet
	}

	my $refRestrictInfo={ACTION=>'restrictSQL',SITE=>$isModifQuantif,LIST=>\%fileFormattedList};
	my ($pvalueCode)=&fetchDesignQuantificationData($projectID,'list',\%dispRatios,\%dispStates,\%ratioCode,\%labelingInfo,\%stateInfo,\%quantifParamID2code,\%quantifValues,\%proteinInfo,\%dispModifSites,\%quantifModifInfo,$refRestrictInfo);

	####<Starting HTML>###
	if (!@_) {  # if @_, it means that this sub is called from another sub and the header is already printed
		print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);
	}
	my $quantifSoftware=($labelingInfo{'SOFTWARE'})? $labelingInfo{'SOFTWARE'}[0] : 'myProMS'; $quantifSoftware='MaxQuant' if $quantifSoftware eq 'MQ';
	if ($quantifSoftware eq 'MaxQuant' || $quantifType eq 'PROT_ABUNDANCE') {
		&listProteinIntensities($labelType,\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,\%dispModifSites);
	}
	elsif ($quantifType eq 'PROT_RULER') {
		&listProteomicRuler($labelType,\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo);
	}
	else {
		&listProteinRatios($labelType,$ratioType,\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,\%dispModifSites,$pvalueCode); #,10
	}
	if (!@_) {  # if @_, it means that this sub is called from another sub and we need to perform other things after
		exit;
	}
}

####################<<<ajaxPeptideRatios>>>#####################
sub ajaxPeptideRatios {
# print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
	# TODO: +/- infinite ratio: Add missing filters for peptides (charge, source)
	my $modProtID=param('id_prot');
	my ($protID,$modStrg)=($modProtID=~/^(\d+)-?(.*)/); $modStrg='' unless $modStrg;
	my $selRatioPos=param('ratio');
	my $trueRatioPos=abs($selRatioPos);
	#my $selRatioIdx=$selRatioPos-1;
#print "PROT=$modProtID<BR>\n";
	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;
	my $projectID=&promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION');
	my $dbhLite=&promsQuantif::dbConnectProteinQuantification($selQuantifID,$projectID);

	#my ($quantifAnnot,$designID,$quantifType,$quantifModID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_DESIGN,M.CODE,Q.ID_MODIFICATION FROM QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$selQuantifID");
	my ($quantifAnnot,$designID,$quantifMethodID,$quantifType,$quantifModID,$multiModifStrg)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_DESIGN,M.ID_QUANTIFICATION_METHOD,M.CODE,Q.ID_MODIFICATION,GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ',')
																	 FROM QUANTIFICATION Q
																	 LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
																	 INNER JOIN QUANTIFICATION_METHOD M ON Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD
																	 WHERE Q.ID_QUANTIFICATION=$selQuantifID GROUP BY Q.ID_QUANTIFICATION");
	#my ($isModifQuantif,$isMultiModifQuantif)=(0,0);
	#my (%quantifModifInfo,%quantifModifNames); #%quantifModifRanks,
	#my ($textModRes,$htmlModRes)=('','');
	#if ($quantifModID || $multiModifStrg) {
	#	$isModifQuantif=1;
	#	my $modifStrg=$quantifModID || $multiModifStrg;
	#	$quantifModID=$modifStrg if $modifStrg !~ /,/; # <- New single PTM quantif (also uses MULTIMODIF_Q table!)
	#	my $sthMod=$dbh->prepare("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DES,DISPLAY_CODE,DISPLAY_COLOR FROM MODIFICATION WHERE ID_MODIFICATION IN ($modifStrg) ORDER BY FIELD(ID_MODIFICATION,$modifStrg)");
	#	my $modifRank=0;
	#	$sthMod->execute;
	#	while (my ($modID,$psiName,$interName,$synName,$des,$displayCode,$displayColor)=$sthMod->fetchrow_array) {
	#			$modifRank++;
	#			#$quantifModifInfo{RANK2ID}{$modifRank}=$modID;
	#			$quantifModifInfo{ID2RANK}{$modID}=$modifRank;
	#			$quantifModifInfo{DISPLAY}{$modID}=[$displayCode,$displayColor];
	#			my $ptmName=$psiName || $interName || $des;
	#			unless ($ptmName) {
	#				$synName=~s/^##//; $synName=~s/##$//; $synName=~s/##/, /;
	#				$ptmName=$synName;
	#			}
	#			$quantifModifNames{$modID}=$ptmName;
	#	}
	#	$sthMod->finish;
	#	$isMultiModifQuantif=1 if $modifRank > 1;
	#}
	my %quantifModifInfo;
	my ($isModifQuantif,$isMultiModifQuantif)=&promsQuantif::getQuantifModificationInfo($dbh,$quantifModID,$multiModifStrg,\%quantifModifInfo);
	
	my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
	my ($labelType)=($labelStrg)? $labelStrg=~/LABEL=(.+)/ : ('FREE');
	$labelType=uc($labelType); # iTRAQ to ITRAQ
	my (%labelingInfo,%condToState,%condToObs);
	$labelingInfo{'LABELTYPE'}=$labelType;
	#if ($designID || $labelType =~ /SILAC|ITRAQ/) {
		foreach my $infoStrg (@labelInfo) {
			my ($setting,$valueStrg)=split('=',$infoStrg);
			$valueStrg=~s/#//g; # remove all ID tags
			@{$labelingInfo{$setting}}=split(';',$valueStrg);
		}
	#}
	my $ratioType=($labelingInfo{RATIO_TYPE})? $labelingInfo{RATIO_TYPE}[0] : 'Ratio';
	my $algoVersion=($ratioType eq 'Ratio')? 1 : ($labelingInfo{SOFTWARE} && $labelingInfo{SOFTWARE}[1])? $labelingInfo{SOFTWARE}[1] : ($labelingInfo{SOFTWARE} && $labelingInfo{SOFTWARE}[0] eq 'myProMS')? 3 : 2;
	my $isPepIntensity=($algoVersion==2 && $labelType eq 'FREE')? 1 : ($algoVersion>=3 && $labelingInfo{ALGO_TYPE}[0]=~/^(PEP_INTENSITY|TDA|DIA)$/)? 1 : 0; # 0 for v1 Ratio, v2 Label & v3 PEP_RATIO
	my $numRatios=scalar @{$labelingInfo{RATIOS}};
	my $numStates=scalar @{$labelingInfo{STATES}};
	my $idxShift=($quantifType eq 'PROT_RATIO_PEP')? 0 : 2;

	if ($modStrg && !$isModifQuantif) { # fall back to whole protID if called with a modProtID for a non-modif quantif (eg. called from log2 plot)
		$modProtID=$protID;
		$modStrg='';
	}
	my ($protAlias,$protLength)=$dbh->selectrow_array("SELECT ALIAS,PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=$protID");
	$protLength=0 unless $protLength;

	my $sthQP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,NAME,CODE FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=$quantifMethodID");
	$sthQP->execute;
	while (my ($paramID,$paramName,$paramCode)=$sthQP->fetchrow_array) {
		@{$quantifParamInfo{$paramCode}}=($paramID,$paramName); # %quantifParamInfo undef until now
	}
	$sthQP->finish;
	
	#my $sthPR;
	my ($textModRes,$htmlModRes,$siteQueryStrg)=('','','');
	if ($modStrg) { # modification quantif
		my $fileModStrg=&promsQuantif::siteCode2QuantifFileCode($modStrg,$isMultiModifQuantif,$quantifModifInfo{ID2RANK});
		$modProtID=$protID.'-'.$fileModStrg; # redefine modProtID to match quantif file format !!!!!!!!!!!!!!!!!!!!!!!!
		$siteQueryStrg="AND SITE_CODE='$fileModStrg'";
		$textModRes=&promsQuantif::displayModificationSites($modStrg,$quantifModifInfo{DISPLAY},'text');
		$htmlModRes=&promsQuantif::displayModificationSites($modStrg,$quantifModifInfo{DISPLAY},'html');
		
		##$sthPR=$dbh->prepare("SELECT QUANTIF_VALUE
		##						FROM PROTEIN_QUANTIFICATION P
		##						LEFT JOIN PROTQUANTIF_MODRES PQMR ON P.ID_PROT_QUANTIF=PQMR.ID_PROT_QUANTIF
		##						LEFT JOIN MODIFIED_RESIDUE MR ON PQMR.ID_MODIF_RES=MR.ID_MODIF_RES
		##						INNER JOIN QUANTIFICATION_PARAMETER Q ON P.ID_QUANTIF_PARAMETER=Q.ID_QUANTIF_PARAMETER AND Q.CODE=?
		##						WHERE P.ID_QUANTIFICATION=$selQuantifID AND ID_PROTEIN=$protID AND TARGET_POS=? GROUP BY P.ID_PROT_QUANTIF
		##						HAVING GROUP_CONCAT(COALESCE(MODIF_RANK,''),':',RESIDUE,POSITION ORDER BY MODIF_RANK,POSITION SEPARATOR '.')='$dbModStrg'");
	}
	# else { # whole protein quantif
	# 	##$sthPR=$dbh->prepare("SELECT QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION P,QUANTIFICATION_PARAMETER Q WHERE P.ID_QUANTIF_PARAMETER=Q.ID_QUANTIF_PARAMETER AND Q.CODE=? AND ID_QUANTIFICATION=$selQuantifID AND ID_PROTEIN=$protID AND TARGET_POS=?");
	# }
	my $sthPR=$dbhLite->prepare("SELECT QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIF_PARAMETER=? AND ID_PROTEIN=$protID AND TARGET_POS=? $siteQueryStrg");

	##<Fetching labeling type
	#my ($quantifAnnot,$designID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_DESIGN FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	#my ($pepQuantifIDs)=$dbh->selectrow_array("SELECT GROUP_CONCAT(ID_PARENT_QUANTIFICATION SEPARATOR ',') FROM PARENT_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID AND PAR_FUNCTION IS NULL GROUP BY ID_QUANTIFICATION");
	#my ($apexRTParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER QP,QUANTIFICATION Q WHERE QP.ID_QUANTIFICATION_METHOD=Q.ID_QUANTIFICATION_METHOD AND CODE='RT_APEX' AND ID_QUANTIFICATION IN ($pepQuantifIDs) LIMIT 1");

	my ($trueProtRatio,$trueTestStateName,$trueRefStateName,$is2ndaryRatio);
	my (%protRatio,%testStatePos,%refStatePos,%testStateName,%refStateName,%usedExpConds,%anaList,%pepQuantifList,%protMean); #,%usedReplicates
	my (%anaTg2StatePos,%usedTargetPos); # algo v1 (ratio)
	my %replicDesign; # v3
	my %formattedModRes; my $dispModSites='';
	if ($designID) { # Design quanti ***statePos=CondID but converted later in true pos by %condToState***
		my ($testStateID,$refStateID)=($selRatioPos > 0)? split(/\//,$labelingInfo{RATIOS}[$trueRatioPos-1]) : reverse (split(/\//,$labelingInfo{RATIOS}[$trueRatioPos-1]));
		my ($testRatioPos,$refRatioPos); # for 2ndary ratio only
		if ($testStateID=~/%/) {  # && $algoVersion < 3) { #} 2NDARY RATIO!!!
			$is2ndaryRatio=1;
			#$testStatePos=~s/%\d+//;
			#$refStatePos=~s/%\d+//;
			($testRatioPos)=($testStateID=~/%(\d+)/);
			($refRatioPos)=($refStateID=~/%(\d+)/);
			foreach my $ratioPos ($testRatioPos,$refRatioPos) {
				($testStatePos{$ratioPos},$refStatePos{$ratioPos})=split(/\//,$labelingInfo{RATIOS}[$ratioPos-1]); # Converted to real pos below
				$usedExpConds{$refStatePos{$ratioPos}}=$usedExpConds{$testStatePos{$ratioPos}}=1;
			}
			$sthPR->execute($quantifParamInfo{'RATIO'}[0],$trueRatioPos);
			($trueProtRatio)=$sthPR->fetchrow_array;
			$trueProtRatio=1/$trueProtRatio if $selRatioPos < 0; # reverse
		}
		else { # PRIMARY RATIO or FREE
			$testStatePos{$selRatioPos}=$testStateID; # Converted to real pos below
			$refStatePos{$selRatioPos}=$refStateID; # Converted to real pos below
			$usedExpConds{$refStateID}=$usedExpConds{$testStateID}=1;
		}
		my $statePos=0;
		foreach my $stateInfo (@{$labelingInfo{STATES}}) {
			$statePos++;
			my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$stateInfo);
			next unless $usedExpConds{$expCondID}; # Restrict to states used for selected (2ndary) ratio (PP 16/01/19)
			my $curBioRep=0;
			foreach my $bioReplicate (split(/\./,$quantiObsIDs)) {
				$curBioRep++;
				$replicDesign{$statePos}{$curBioRep}=($quantiObsIDs=~/&/)? 1 : 0;
				next if (!$isPepIntensity && $statePos==1);	# Use only test state for label design to prevent use of State1 (possibly more obs than used by selRatioPos)			
				foreach my $techReplicate (split(/&/,$bioReplicate)) { # tech rep separator
					foreach my $fraction (split(/\+/,$techReplicate)) {
						my ($obsID,$pepQuantID,$anaID,$tgPos)=split(':',$fraction);
						$pepQuantifList{$pepQuantID}=1;	
						$anaList{$anaID}=1;
						
						#(Ratio) Algo v1
						$tgPos=0 unless $tgPos;
						$anaTg2StatePos{"$anaID.$tgPos"}=$statePos; # absolute pos in all states list
						$usedTargetPos{$tgPos}=1;

					}
				}
			}			
			#$expCondID=~ s/^.//;
			$condToState{$expCondID}=$statePos; # absolute pos in all states list
			$condToObs{$expCondID}=$quantiObsIDs;
		}
		my $sthExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");

		foreach my $ratioPos (keys %testStatePos) {
			$sthPR->execute($quantifParamInfo{'RATIO'}[0],abs($ratioPos));
			my ($pRatio)=$sthPR->fetchrow_array;
			if ($is2ndaryRatio || $selRatioPos > 0) {$protRatio{$ratioPos}=$pRatio;}
			else {$protRatio{$ratioPos}=1/$pRatio;}
			my $testStateExpCondID=$testStatePos{$ratioPos};
			$sthExpCondName->execute($testStateExpCondID);
			($testStateName{$ratioPos})=$sthExpCondName->fetchrow_array;
			my $refStateExpCondID=$refStatePos{$ratioPos};
			$sthExpCondName->execute($refStateExpCondID);
			($refStateName{$ratioPos})=$sthExpCondName->fetchrow_array;
			$testStatePos{$ratioPos}=$condToState{$testStateExpCondID}; # absolute pos in all states list
			$refStatePos{$ratioPos}=$condToState{$refStateExpCondID}; # absolute pos in all states list
			#if ($protRatio{$ratioPos} >= $MIN_INF_RATIO_DB) {$usedReplicates{$ratioPos}=$condToObs{$refStateExpCondID};}
			#elsif ($protRatio{$ratioPos} <= $MAX_INF_RATIO_DB) {$usedReplicates{$ratioPos}=$condToObs{$testStateExpCondID};}
			if ($isPepIntensity) { # $labelType eq 'FREE' && $ratioType ne 'Ratio'
				$sthPR->execute($quantifParamInfo{'MEAN_STATE'}[0],$numRatios+$testStatePos{$ratioPos});
				($protMean{$testStatePos{$ratioPos}})=$sthPR->fetchrow_array;
				$sthPR->execute($quantifParamInfo{'MEAN_STATE'}[0],$numRatios+$refStatePos{$ratioPos});
				($protMean{$refStatePos{$ratioPos}})=$sthPR->fetchrow_array;
#< Check for bug in state mean allocation => switch if necessary! (except if +/-INF)				
if ($protMean{$testStatePos{$ratioPos}} && $protMean{$refStatePos{$ratioPos}} && (($pRatio < 1 && $protMean{$testStatePos{$ratioPos}} > $protMean{$refStatePos{$ratioPos}}) || ($pRatio > 1 && $protMean{$testStatePos{$ratioPos}} < $protMean{$refStatePos{$ratioPos}})) ) {
	my $meanProt=$protMean{$testStatePos{$ratioPos}};
	$protMean{$testStatePos{$ratioPos}}=$protMean{$refStatePos{$ratioPos}};
	$protMean{$refStatePos{$ratioPos}}=$meanProt;
}
#print "TEST=$protMean{$testStatePos{$ratioPos}}, REF=$protMean{$refStatePos{$ratioPos}}<BR>\n";
			}
		}
		$sthExpCondName->finish;
		if ($is2ndaryRatio) { # SuperRatio
			$trueTestStateName=$testStateName{$testRatioPos}.$encodedDegStrg;
			$trueRefStateName=$testStateName{$refRatioPos}.$encodedDegStrg;
		}
		else {
			$trueTestStateName=$testStateName{$selRatioPos};
			$trueRefStateName=$refStateName{$selRatioPos};
		}
	}
	else { # internal quanti
		$sthPR->execute($quantifParamInfo{'RATIO'}[0],$trueRatioPos);
		($protRatio{$selRatioPos})=$sthPR->fetchrow_array;
		($testStatePos{$selRatioPos},$refStatePos{$selRatioPos})=($selRatioPos > 0)? split(/\//,$labelingInfo{RATIOS}[$trueRatioPos-1]) : reverse(split(/\//,$labelingInfo{RATIOS}[$trueRatioPos-1]));
		($testStateName{$selRatioPos},my $testChanNum)=$labelingInfo{STATES}[$testStatePos{$selRatioPos}-1]=~/,(.+),(\d+)/; $testStateName{$selRatioPos}=~s/\./\+/g;
		($refStateName{$selRatioPos},my $refChanNum)=$labelingInfo{STATES}[$refStatePos{$selRatioPos}-1]=~/,(.+),(\d+)/; $refStateName{$selRatioPos}=~s/\./\+/g;
my ($anaID)=$dbh->selectrow_array("SELECT ID_ANALYSIS FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
$anaTg2StatePos{"$anaID.$testChanNum"}=$testStatePos{$selRatioPos};
$anaTg2StatePos{"$anaID.$refChanNum"}=$refStatePos{$selRatioPos};
$usedTargetPos{$refChanNum}=$usedTargetPos{$testChanNum}=1;
$anaList{$anaID}=1;
#print "\$selRatioPos=$selRatioPos; \$testStatePos{\$selRatioPos}=$testStatePos{$selRatioPos}, \$testStateName{\$selRatioPos}=$testStateName{$selRatioPos}<BR>\n";
# Define $usedReplicates for infinite ratio !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		#if ($protRatio <= 0.001 || $protRatio >= 1000) { # infinite ratio
		#	my ($stateChannelStrg)=($protRatio <= 0.001)? $labelingInfo{STATES}[$refStatePos-1]=~/,([^,]+)\Z/ : $labelingInfo{STATES}[$testStatePos-1]=~/,([^,]+)\Z/;
		#	# get state data from parent peptide quanti
		#	# sthQuery....
		#	foreach my $channel (split(/\./,$stateChannelStrg)) {
		#		#...
		#		# execute($channel);
		#
		#	}
		#	# ...
		#}
		#($usedReplicates{$selRatioPos})=($protRatio{$selRatioPos} <= 0.001)? $labelingInfo{STATES}[$refStatePos{$selRatioPos}-1]=~/,([^,]+)\Z/ : $labelingInfo{STATES}[$testStatePos{$selRatioPos}-1]=~/,([^,]+)\Z/;
		my ($pepQuantifID)=$dbh->selectrow_array("SELECT ID_PARENT_QUANTIFICATION FROM PARENT_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
		$pepQuantifList{$pepQuantifID}=1;
	}
	$sthPR->finish;
	
	if ($modStrg) {
		$modStrg='-'.$modStrg;
		$htmlModRes='-'.$htmlModRes;
		$textModRes='-'.$textModRes;
	}

	#my $sthPB=$dbh->prepare("SELECT ABS(PEP_BEG) FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_PROTEIN=$protID AND ID_PEPTIDE=? ORDER BY ABS(PEP_BEG) ASC");
	#my $pepQuery=($apexRTParamID)? "SELECT SCORE,DATA,ID_ANALYSIS,QUANTIF_VALUE,QUERY_NUM,PEP_RANK FROM PEPTIDE P
	#									LEFT JOIN PEPTIDE_QUANTIFICATION PQ ON P.ID_PEPTIDE=PQ.ID_PEPTIDE AND ID_QUANTIF_PARAMETER=$apexRTParamID AND ID_QUANTIFICATION IN ($pepQuantifIDs)
	#									WHERE P.ID_PEPTIDE=?"
	#								: "SELECT SCORE,DATA,ID_ANALYSIS,NULL,QUERY_NUM,PEP_RANK FROM PEPTIDE P,PEPTIDE_QUANTIFICATION PQ WHERE P.ID_PEPTIDE=PQ.ID_PEPTIDE AND P.ID_PEPTIDE=? AND ID_QUANTIFICATION IN ($pepQuantifIDs)";
	#my $sthPD=$dbh->prepare($pepQuery);

	my $dataDir="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/data";
	my $resultDir="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/results";
	
	my $pepQuantifIDs=join(',',sort{$a<=>$b} keys %pepQuantifList);
#print "'$pepQuantifIDs'<BR>\n";
	my ($apexRTParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER QP,QUANTIFICATION Q WHERE QP.ID_QUANTIFICATION_METHOD=Q.ID_QUANTIFICATION_METHOD AND CODE='RT_APEX' AND ID_QUANTIFICATION IN ($pepQuantifIDs) LIMIT 1");
	my (%protPeptideData,%peptideApexRT);
	my $sthProtPep=$dbh->prepare(qq
|SELECT P.ID_PEPTIDE,GROUP_CONCAT(DISTINCT ABS(PEP_BEG) ORDER BY ABS(PEP_BEG) ASC SEPARATOR ','),SCORE,DATA,P.ID_ANALYSIS,QUERY_NUM,PEP_RANK,VALID_STATUS
	FROM PEPTIDE P,PEPTIDE_PROTEIN_ATTRIB PPA,ANA_QUANTIFICATION AQ
	WHERE P.ID_PEPTIDE=PPA.ID_PEPTIDE AND PPA.ID_PROTEIN=? AND P.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION IN ($pepQuantifIDs) GROUP BY P.ID_PEPTIDE
|);
	my ($usePtmProb,%ptmProb,$ptmSoft);
	my $existSpectrumFile=1; # default
	my $fileFormat = "";
	$sthProtPep->execute($protID);
	while (my ($pepID,$begStrg,$score,$pepData,$anaID,$qNum,$rank,$pepStatus)=$sthProtPep->fetchrow_array) {
		$score=0 unless $score;
		$pepData='' unless $pepData;
		@{$protPeptideData{$pepID}{BEG}}=split(',',$begStrg);
		@{$protPeptideData{$pepID}{INFO}}=($score,$pepData,$anaID,$qNum,$rank,$pepStatus);
		if (!defined $usePtmProb) { # do only once!
			($fileFormat)=$dbh->selectrow_array("SELECT FILE_FORMAT FROM ANALYSIS WHERE ID_ANALYSIS=$anaID") unless($fileFormat);
			if ($fileFormat =~ /MAXQUANT.DIR|SPECTRONAUT|SEQUEST/) {
				$usePtmProb=($isModifQuantif)? 1 : 0;
				$existSpectrumFile=0 if($fileFormat eq "MAXQUANT.DIR" && !-e "$promsPath{peptide}/proj_$projectID/ana_$anaID/msms.txt");
			}
			else {$usePtmProb=0;}
		}
		%{$ptmProb{$pepID}}=() if $usePtmProb;
	}
	$sthProtPep->finish;

	if ($usePtmProb) { # MaxQuant/Spectronaut/ptmRS PTM probability
		#($quantifModName)=$dbh->selectrow_array("SELECT PSI_MS_NAME FROM MODIFICATION WHERE ID_MODIFICATION=$quantifModID"); # !!! TODO: HANDLE MULTI-MODIF !!!
		my $sthMQP=$dbh->prepare("SELECT ID_MODIFICATION,POS_STRING,REF_POS_STRING FROM PEPTIDE_MODIFICATION WHERE ID_MODIFICATION IN (".join(',',keys %{$quantifModifInfo{NAME}}).") AND ID_PEPTIDE=?");
		foreach my $pepID (keys %ptmProb) {
			$sthMQP->execute($pepID);
			while (my ($modID,$modPosStrg,$refPosStrg)=$sthMQP->fetchrow_array) {
				if ($refPosStrg && $refPosStrg=~/##PRB_(MQ|PTMRS|SPC)=([^#]+)/) {
					($ptmSoft, my $modProbStrg)=($ptmProbSoft{$1}, $2);
					while ($modProbStrg =~ /([^,:]+):([^,]+)/g) {
						push @{$ptmProb{$pepID}{$modID}},[$1,$2];
					}
				}
				elsif($fileFormat eq "MAXQUANT.DIR") { # all 100%
					foreach my $pos (split(/\./,$modPosStrg)) {
						push @{$ptmProb{$pepID}{$modID}},[$pos,1]; # set to 100%
					}
				}
			}
			delete $ptmProb{$pepID} unless scalar keys %{$ptmProb{$pepID}};
		}
		$sthMQP->finish;
	}

	if ($apexRTParamID) {
		QUANTIF:foreach my $pepQuantifID (split(',',$pepQuantifIDs)) {
			foreach my $quantFile (glob ("$promsPath{quantification}/project_$projectID/quanti_$pepQuantifID/peptide_quantification*")) {
				open(QUANTI,$quantFile);
				while (<QUANTI>) {
					next if $.==1;
					my ($paramID,$pepID,$qValue)=split(/\t/,$_);
					$peptideApexRT{$pepID}=$qValue if ($paramID==$apexRTParamID && $protPeptideData{$pepID});
				}
				close QUANTI;
				last QUANTIF if scalar keys %peptideApexRT == scalar keys %protPeptideData;
			}
		}
	}

	my (%peptideSets,%sequenceBeg,%dataSources);
my %dataSrcAlias; # v3
	my (%rawPeptideSets,%pepOccurence,%varModName); #,%sequenceBeg,%dataSources

	######################
	####<Normal ratio>####
	######################
	my %normalRatio;
my %usedAnalyses; # anaID used for selected ratio (needed for filtering table.txt for SuperRatio. No longer use Experiment) v3
my %usedStates; # (Super/Simple)Ratio (v2 v3)
	my $numNormalRatios=0;
	###my ($pepFile,$pepSetFactor)=($ratioType eq 'Ratio')? ("$resultDir/pep_mean_cv_ratio.txt",2) : ($labelType eq 'FREE')? ("$resultDir/Log2MeanPep.txt",1) : ("$resultDir/RatioPep.txt",2);
	#my $pepSetFactor=($ratioType eq 'Ratio' || $labelType ne 'FREE')? 2 : 1;	
	my $pepFile=($algoVersion>=3)? 'resultsPep.txt' : ($ratioType eq 'Ratio')? 'pep_mean_cv_ratio.txt' : ($labelType eq 'FREE')? 'Log2MeanPep.txt' : 'RatioPep.txt';
	$pepFile="$resultDir/$pepFile";
	foreach my $ratioPos (keys %protRatio) {
		$normalRatio{$ratioPos}=0;
		if ($protRatio{$ratioPos} && $protRatio{$ratioPos} != 0.001 && $protRatio{$ratioPos} != 1000) { # test on defined protRatio because of rare possible exclusion of all peptides by global CVout
			$normalRatio{$ratioPos}=1;
			$numNormalRatios++;
			if ($ratioType eq 'Ratio') { # v1
				my %varModName;
				my $prevProtID=0;
				open(PEP,$pepFile) || die "ERROR: $!";
				while(<PEP>) {
					next if $.==1;
					chomp;
					my ($pID,$pepSeq,$varMod,$charge,$pepIdKey,$protName,$protValid,@values)=split(/\t/,$_); # ratio	mean.cond1	mean.cond2	PERM.Repl_Light	PERM.Repl_Heavy	PERM_Out	Pep_Out_0.05
					my ($trueProtID)=($pID=~/^(\d+)/);
					last if ($prevProtID==$protID && $trueProtID != $protID); # sorted by true protID but not by modProtID
					next if $pID ne $modProtID;
					$prevProtID=$trueProtID;
					my $seqVarMod=$pepSeq;
					if ($varMod && $varMod ne 'NA') {
						if ($varMod=~/^&/) { # var mod code
							$varMod=&promsMod::decodeVarMod($dbh,$pepSeq,$varMod,\%varModName);
						}
						else { # var mod text
							$varMod=~s/^ \+ //;
						}
						$seqVarMod.=' + '.$varMod;
					}
					else {$varMod='';}
					unless ($sequenceBeg{$pepSeq}) {
						my ($pep1ID)=$pepIdKey=~/^(\d+)/;
						#$sthPB->execute($pep1ID);
						#while (my ($beg)=$sthPB->fetchrow_array) {
						#	push @{$sequenceBeg{$pepSeq}},$beg;
						#}
@{$sequenceBeg{$pepSeq}}=@{$protPeptideData{$pep1ID}{BEG}};
					}
					my ($dataSrc,$scoreStrg,$scoreMin,$scoreMax); #,$prsStrg,$prsMin,$prsMax,$rtStrg,$rtMin,$rtMax;
					$pepIdKey=~s/\s//g; # very old quantif with spaces around ';'
					foreach my $pepID (split(/\D/,$pepIdKey)) {
						#$sthPD->execute($pepID);
						#my ($score,$pepData,$anaID)=$sthPD->fetchrow_array;
						#$score=0 unless $score;
						#$pepData='' unless $pepData;
my ($score,$pepData,$anaID)=@{$protPeptideData{$pepID}{INFO}}[0..2];
						#$scoreStrg.='/' if $scoreStrg;
						#$scoreStrg.=($score)? $score : '-';
						unless ($dataSrc) {
							$dataSrc='#'.$anaID;
							$dataSrc.='.'.$1 if $pepData=~/SOURCE_RANK=(\d+)/;
						}
						if ($score) {
							$scoreMin=$score unless ($scoreMin && $score>$scoreMin);
							$scoreMax=$score unless ($scoreMax && $score<$scoreMax);
						}
					#if ($quantifModID && $pepData=~/PRS=\d;([^;]+);/) {
					#	my $prs=$1;
					#	$prsMin=$prs unless ($prsMin && $prsMin>$prs);
					#	$prsMax=$prs unless ($prsMax && $prsMax<$prs);
					#}
						#$dataSrc=$dSrc || '-';
	#$anaList{$anaID}=1;
					}
					$scoreStrg='-';
					if ($scoreMin && $scoreMax) {
						if ($scoreMin == $scoreMax) {$scoreStrg=1*(sprintf '%.2f',$scoreMin);}
						else {$scoreStrg='['.(1*(sprintf '%.2f',$scoreMin)).'-'.(1*(sprintf '%.2f',$scoreMax)).']';}
					}
				#if ($quantifModID) {
				#	$prsStrg='-';
				#	if ($prsMin && $prsMax) {
				#		if ($prsMin == $prsMax) {$prsStrg=1*(sprintf '%.1f',$prsMin);}
				#		else {$prsStrg='['.(1*(sprintf '%.1f',$prsMin)).'-'.(1*(sprintf '%.1f',$prsMax)).']';}
				#	}
				#}
					$dataSrc='' unless $dataSrc;
					$dataSources{$dataSrc}=$dataSrc;
					@{$peptideSets{$ratioPos}{$pepIdKey}{INFO}}=($sequenceBeg{$pepSeq}[0],$pepSeq,$varMod,$charge,$scoreStrg,{NAME=>'',DATA=>[$dataSrc]});
					my $pepRatio=($selRatioPos > 0)? $values[$trueRatioPos-1] : 1/$values[$trueRatioPos-1];
					push @{$peptideSets{$ratioPos}{$pepIdKey}{DATA}},$pepRatio,
						$values[$numRatios+$testStatePos{$ratioPos}-1],$values[$numRatios+$refStatePos{$ratioPos}-1], # mean_up, mean_down
						$values[$numRatios+$numStates+$testStatePos{$ratioPos}-1],$values[$numRatios+$numStates+$refStatePos{$ratioPos}-1],$values[$numRatios+($numStates*2)], # cv_out_up, cv_out_down, cv_out
						$values[$numRatios+($numStates*2)+$trueRatioPos+$idxShift]; # pep_out_ratio
				}
				close PEP;
			}
			else { # (Super/Simple)Ratio
				my %selStateCode; # v3
				my $rPos=0;
				PR:foreach my $x (1..$numStates-1) {
					foreach my $y ($x+1..$numStates) {
						$rPos++;
						if ($rPos==$trueRatioPos) {
							if ($algoVersion>=3) {
								if ($isPepIntensity) {
									$selStateCode{"State$x"}=$x; #$ratioPos;
									$selStateCode{"State$y"}=$y; #$ratioPos;
								}
								else {$selStateCode{"State$y.State$x"}=$ratioPos;}
							}
							$usedStates{"State$x"}=$usedStates{"State$y"}=1;
							last PR;
						}
					}
					last if $is2ndaryRatio; # no need to check $labelInfo{SINGLE_REF}
				}
				if ($is2ndaryRatio) {
					SR:foreach my $x (2..$numStates-1) {
						foreach my $y ($x+1..$numStates) {
							$rPos++;
							if ($rPos==$trueRatioPos) {
								if ($algoVersion>=3) {
									$selStateCode{"State$x.State1"}=$x-1; # primary ratioPos
									$selStateCode{"State$y.State1"}=$y-1; # primary ratioPos
								}
								$usedStates{"State$x"}=$usedStates{"State$y"}=$usedStates{'State1'}=1;
								last SR;
							}
						}
					}
				}

				my (%peptideData, %seenExcluded); # simpleRatio with isPepIntensity
				open(PEP,$pepFile) || die "ERROR: $!";
				if ($algoVersion>=3) {
					my (%pepHeader2Idx,%replicates);
					while(<PEP>) {
						#next if $.==1;
						chomp;
						#my ($condCode,$experiment,$bioRep,$techRep,$pID,$pepCode,$log2Value,$pepIdKey,$outlier)=split(/\t/,$_);
						#ProteinID	Peptide	Experiment	replicate	repTech	Condition	proteinName	PeptideId	proteinValidity	log2Measure	normProtein	out
						#my ($pID,$pepCode,$experiment,$bioRep,$techRep,$condCode,$protName,$pepIdKey,$protValid,$log2Value,$protNorm,$outlier)=split(/\t/,$_);
						#ProteinID	Peptide	Experiment	Condition	proteinValidity	log2Measure	replicate	repTech	proteinName	PeptideId	normProtein	out # v3b
						my @lineData=split(/\t/,$_);
						####>Header------------------
						if ($.==1) { # 1st line of the file
							my $colIdx=0;
							foreach my $colName (@lineData) {
								$pepHeader2Idx{$colName}=$colIdx;
								$colIdx++;
							}
							next;
						}
						####>Data------------------
						next if $lineData[$pepHeader2Idx{'ProteinID'}] ne $modProtID;
						next unless $selStateCode{$lineData[$pepHeader2Idx{'Condition'}]};
						my $rPos=$selStateCode{$lineData[$pepHeader2Idx{'Condition'}]}; # always primary ratio or meanProt!
						my $isOutlier=($lineData[$pepHeader2Idx{'out'}] eq 'NA')? 0 : 1;
						@{$peptideData{$rPos}{$lineData[$pepHeader2Idx{'Peptide'}]}{$lineData[$pepHeader2Idx{'replicate'}].$lineData[$pepHeader2Idx{'repTech'}]}}=($lineData[$pepHeader2Idx{'log2Measure'}],$lineData[$pepHeader2Idx{'PeptideId'}],$isOutlier);
my ($statePos)=$lineData[$pepHeader2Idx{'Condition'}]=~/^State(\d+)/;
my ($bioRepPos)=$lineData[$pepHeader2Idx{'replicate'}]=~/(\d+)/;
my ($techRepPos)=$lineData[$pepHeader2Idx{'repTech'}]=~/(\d+)/;
my $numBioRep=scalar keys %{$replicDesign{$statePos}};
$dataSrcAlias{$lineData[$pepHeader2Idx{'replicate'}].$lineData[$pepHeader2Idx{'repTech'}]}=($numBioRep==1 && $replicDesign{$statePos}{$bioRepPos}==0)? 'No replic.': ($numBioRep==1)? "tech.rep$techRepPos" : ($replicDesign{$statePos}{$bioRepPos}==0)? "bio.rep$bioRepPos" : "bio.rep$bioRepPos/tech.rep$techRepPos";
					}
				}
				else { # v2+ No state recorded: only experiment!
					my $selExp=chr(64+abs($ratioPos)); # cannot be 2ndary ratio!!! OK for iTRAQ in results files
					while(<PEP>) {
						next if $.==1;
						chomp;
						my ($expState,$pID,$pepCode,$log2Value,$pepIdKey,$pepExcluded,$nbPep,$log2PepStrg)=split(/\t/,$_); # $log2PepStrg defined only for label-free since FunctionLimma.R v2.4.1
						##my ($trueProtID)=($pID=~/^(\d+)/);
						##last if ($prevProtID==$protID && $trueProtID != $protID); # sorted by true protID but maybe not by modProtID <- NEW: Not even true for LABEL-FREE!!!
						next if $pID ne $modProtID;
						##$prevProtID=$trueProtID;
						next if $log2Value eq 'NA'; # No $log2Value for excluded peptide!!!!!! TODO: record a value here
						my $statePos;
						my $usedPos;
						if ($labelType eq 'FREE') {
							($statePos)=($expState=~/(\d+)/);
							next unless $protMean{$statePos};
							$usedPos=$statePos;
						}
						else {
							if ($expState ne $selExp) {next;}
							$usedPos=$ratioPos;
						}
						@{$peptideData{$usedPos}{$pepCode}{ALL_REP}}=($log2Value,$pepIdKey,$pepExcluded);
					}
				}

				close PEP;

				foreach my $usedPos (keys %peptideData) {
					foreach my $pepCode (keys %{$peptideData{$usedPos}}) {

						my ($seqModCode,$charge)=split(/_/,$pepCode);
						my ($pepSeq,$varMod);
						if ($seqModCode=~/(\w+)&(.+)/) {
							$pepSeq=$1;
							$varMod=&promsMod::decodeVarMod($dbh,$pepSeq,$2,\%varModName);
						}
						else {
							$pepSeq=$seqModCode;
							$varMod='';
						}

						foreach my $repCode (keys %{$peptideData{$usedPos}{$pepCode}}) {
							my ($log2Value,$pepIdKey,$pepExcluded)=@{$peptideData{$usedPos}{$pepCode}{$repCode}};
				my @pepIdList=split(/[+_;\.=]/,$pepIdKey);
				next unless $protPeptideData{$pepIdList[0]}; # peptide does not belong to selected dataset
							unless ($sequenceBeg{$pepIdList[0]}) {
#print "$pepIdKey: $pepIdList[0]<BR>\n";
								#$sthPB->execute($pepIdList[0]);
								#while (my ($beg)=$sthPB->fetchrow_array) {
								#	push @{$sequenceBeg{$pepSeq}},$beg;
								#}
@{$sequenceBeg{$pepSeq}}=@{$protPeptideData{$pepIdList[0]}{BEG}};
							}
						#my @pepIdList=($pepIdKey);
						#my @pepLog2Mean=($log2Value);
						#if ($log2PepStrg) { # Add unagreggated peptides
						#	$log2PepStrg=~s/^_//;
						#	if ($pepIdKey=~/_/) { # or $nbPep > 1?
						#		push @pepIdList,split(/_/,$pepIdKey);
						#		push @pepLog2Mean,split(/_/,$log2PepStrg);
						#	}
						#}
						#my %dataSrcPepIndex; # only if $log2PepStrg is defined
							my (%dataSrc,$scoreStrg,$scoreMin,$scoreMax); #,$prsStrg,$prsMin,$prsMax,$rtStrg,$rtMin,$rtMax;

							foreach my $pepID (@pepIdList) {
								#$sthPD->execute($pepID);
								#my ($score,$pepData,$anaID)=$sthPD->fetchrow_array;
								#$score=0 unless $score;
								#$pepData='' unless $pepData;
my ($score,$pepData,$anaID)=@{$protPeptideData{$pepID}{INFO}}[0..2];
								#$scoreStrg.='/' if $scoreStrg;
								#$scoreStrg.=($score)? $score : '-';
								if ($score) {
									$scoreMin=$score unless ($scoreMin && $score>$scoreMin);
									$scoreMax=$score unless ($scoreMax && $score<$scoreMax);
								}
							#if ($quantifModID && $pepData=~/PRS=\d;([^;]+);/) {
							#	my $prs=$1;
							#	$prsMin=$prs unless ($prsMin && $prsMin>$prs);
							#	$prsMax=$prs unless ($prsMax && $prsMax<$prs);
							#}
								#$dSrc='-' unless $dSrc;
								my $dSrc='#'.$anaID;
								$dSrc.='.'.$1 if $pepData=~/SOURCE_RANK=(\d+)/;
								$dataSrc{$dSrc}=$dSrc;
								$dataSources{$dSrc}=$dSrc;
								#$anaList{$anaID}=1;
							}
							$scoreStrg='-';
							if ($scoreMin && $scoreMax) {
								if ($scoreMin == $scoreMax) {$scoreStrg=1*(sprintf '%.2f',$scoreMin);}
								else {$scoreStrg='['.(1*(sprintf '%.2f',$scoreMin)).'-'.(1*(sprintf '%.2f',$scoreMax)).']';}
							}
						#if ($quantifModID) {
						#	$prsStrg='-';
						#	if ($prsMin && $prsMax) {
						#		if ($prsMin == $prsMax) {$prsStrg=1*(sprintf '%.1f',$prsMin);}
						#		else {$prsStrg='['.(1*(sprintf '%.1f',$prsMin)).'-'.(1*(sprintf '%.1f',$prsMax)).']';}
						#	}
						#}
							my $pepUsed=($pepExcluded)? 0 : 1;
							my @srcList=keys %dataSrc;
							my %srcInfo=(NAME=>$dataSrcAlias{$repCode} || $repCode,DATA=>\@srcList);
							#if ($labelType eq 'FREE') { #}
							if ($isPepIntensity) {
								if ($pepExcluded && $algoVersion==2) { # v2: special case: same occurences are on different lines
									push @{$seenExcluded{$usedPos}{$pepCode}},[$pepIdKey,$sequenceBeg{$pepSeq}[0],$pepSeq,$varMod,$charge,$scoreStrg,\%srcInfo,2**$log2Value]; # (keys %dataSrc)[0]
								}
								else {
									@{$peptideSets{$usedPos}{$pepIdKey}{INFO}}=($sequenceBeg{$pepSeq}[0],$pepSeq,$varMod,$charge,$scoreStrg,\%srcInfo); # join('<BR>',sort keys %dataSrc)
									push @{$peptideSets{$usedPos}{$pepIdKey}{DATA}},(2**$log2Value,scalar (split(/[+_;\.=]/,$pepIdKey)),$pepUsed);
								}
							}
							else {
								@{$peptideSets{$usedPos}{$pepIdKey}{INFO}}=($sequenceBeg{$pepSeq}[0],$pepSeq,$varMod,$charge,$scoreStrg,\%srcInfo); # join('<BR>',sort keys %dataSrc)
								my $ratioVal=($ratioPos > 0)? 2**$log2Value : 2**(-$log2Value);
								push @{$peptideSets{$usedPos}{$pepIdKey}{DATA}},($ratioVal,(scalar (split(/[+_;\.=]/,$pepIdKey)))/2,$pepUsed);
							}
						}
					}
				}
				#$dataSources{'-'}='...'; # no datasource distinction for (Super/simple)Ratio
				if (scalar keys %seenExcluded) { # Label-free with SimpleRatio (algo v2 only)
					foreach my $statePos (keys %seenExcluded) {
						foreach my $pepCode (keys %{$seenExcluded{$statePos}}) {
							my (@pepID,$seqBeg,$pepSeq,$varMod,$charge,@scores,%dataSrc,$meanValue,$occ);
							foreach my $refPep (@{$seenExcluded{$statePos}{$pepCode}}) {
								push @pepID,$refPep->[0];
								($seqBeg,$pepSeq,$varMod,$charge)=@{$refPep}[1..4];
								push @scores,$refPep->[5] if $refPep->[5] ne '-';
								foreach my $src (keys %{$refPep->[6]}) {$dataSrc{$src}=1;}
								$meanValue+=$refPep->[7];
								$occ++;
							}
							my $pepIdKey=join('_',sort{$a<=>$b} @pepID);
							my $scoreStrg;
							my $numScores=scalar @scores;
							if ($numScores==0) {$scoreStrg='-';}
							elsif ($numScores==1) {$scoreStrg=$scores[0];}
							else {
								my ($scoreMin,$scoreMax)=(sort{$a<=>$b}@scores)[0,-1];
								$scoreStrg='['.(1*(sprintf '%.2f',$scoreMin)).'-'.(1*(sprintf '%.2f',$scoreMax)).']';
							}
							my @srcList=keys %dataSrc;
							@{$peptideSets{$statePos}{$pepIdKey}{INFO}}=($seqBeg,$pepSeq,$varMod,$charge,$scoreStrg,{NAME=>'',DATA=>\@srcList}); # join('<BR>',sort keys %dataSrc)
							push @{$peptideSets{$statePos}{$pepIdKey}{DATA}},($meanValue/$occ,$occ,0);
						}
					}
				}

			}
		}
	}

	##########################################
	####<Raw data (former Infinite ratio)>####
	##########################################
	if ($ratioType eq 'Ratio') {
		###<1: Assume data is in table.txt (new quantif procedure: runXICProtQuantification v2.2.3+)
		###open(DATA,"$dataDir/table.txt") || die "$dataDir<BR>ERROR: $!";
		###my $prevProtID=0;
		###while(<DATA>) {
		###	next if $.==1;
		###	chomp;
		###	my ($mpID,$pepSeq,$varMod,$charge,$pepIdKey,$identifier,$validity,@values)=split(/\t/,$_);
		###	my ($pID)=($mpID=~/^(\d+)/);
		###	last if ($pID != $protID && $prevProtID==$protID); # modProtID not always in continuous order (but protID is!)
		###	$prevProtID=$pID;
		###	next if $mpID ne $modProtID;
		###	if ($varMod) {
		###		if ($varMod=~/^&/) { # var mod code
		###			$varMod=&decodeVarMod($dbh,$pepSeq,$varMod,\%varModName);
		###		}
		###		else { # var mod text
		###			$varMod=~s/^ \+ //;
		###		}
		###		#$varMod=' + '.$varMod;
		###	}
		###	else {$varMod='';}
		###	unless ($sequenceBeg{$pepSeq}) {
		###		my ($pep1ID)=$pepIdKey=~/^(\d+)/;
		###		$sthPB->execute($pep1ID);
		###		while (my ($beg)=$sthPB->fetchrow_array) {
		###			push @{$sequenceBeg{$pepSeq}},$beg;
		###		}
		###	}
		###	my ($dataSrc,$scoreMin,$scoreMax,$prsMin,$prsMax,$rtMin,$rtMax);
		###	foreach my $pepID (split(/\s*;\s*/,$pepIdKey)) {
		###		$sthPD->execute($pepID);
		###		my ($score,$pepData,$anaID,$rtApex,$qNum,$rk,$tgPos)=$sthPD->fetchrow_array;
		###		$score=0 unless $score;
		###		$pepData='' unless $pepData;
		###		#$scoreStrg.='/' if $scoreStrg;
		###		#$scoreStrg.=($score)? $score : '-';
		###		if ($score) {
		###			$scoreMin=$score if (!$scoreMin || $scoreMin > $score);
		###			$scoreMax=$score if (!$scoreMax || $scoreMax < $score);
		###		}
		###		unless ($dataSrc) {
		###			$dataSrc='#'.$anaID;
		###			$dataSrc.='.'.$1 if $pepData=~/SOURCE_RANK=(\d+)/;
		###		}
		###		if ($quantifModID && $pepData=~/PRS=\d;([^;]+);/) {
		###			my $prs=$1;
		###			$prsMin=$prs if (!$prsMin || $prsMin > $prs);
		###			$prsMax=$prs if (!$prsMax || $prsMax < $prs);
		###		}
		###		if ($rtApex) {
		###			$rtMin=$rtApex if (!$rtMin || $rtMin > $rtApex);
		###			$rtMax=$rtApex if (!$rtMax || $rtMax < $rtApex);
		###		}
		###		#$dataSrc=$dSrc || '-'; # Not good for LABEL-FREE!!!! (!= for refState & testState)
		###		$anaList{$anaID}=1;
		###	}
		###	#$dataSrc='' unless $dataSrc;
		###	$dataSources{$dataSrc}=$dataSrc;
		###	my $uniquePep=$pepSeq.$varMod.'_'.$charge;
		###	foreach my $xic (@values) {
		###		$pepOccurence{$selRatioPos}{$uniquePep}++ if $xic ne 'NA';
		###	}
		###	my $pepCode=$pepSeq.':'.$varMod.'_'.$charge;
		###	$pepCode.="_$dataSrc" if $labelType ne 'FREE';
		###	@{$rawPeptideSets{$selRatioPos}{$pepCode}{INFO}}=($sequenceBeg{$pepSeq}[0],$pepSeq,$varMod,$charge,$uniquePep);
		###	push @{$rawPeptideSets{$selRatioPos}{$pepCode}{DATA}{$refStatePos{$selRatioPos}}},[$scoreMin,$values[$refStatePos{$selRatioPos}-1],$rtMax,$dataSrc,undef,undef,undef];
		###	push @{$rawPeptideSets{$selRatioPos}{$pepCode}{DATA}{$testStatePos{$selRatioPos}}},[$scoreMax,$values[$testStatePos{$selRatioPos}-1],$rtMin,$dataSrc,undef,undef,undef];
		###}


#		my %anaTg2StatePos;
#my %usedTargetPos;
#		if ($designID) {
#			foreach my $expCondID (keys %condToObs) {
#next if ($condToState{$expCondID} != $testStatePos{$selRatioPos} && $condToState{$expCondID} != $refStatePos{$selRatioPos}); # Algo v1: No 2ndary ratios
#				foreach my $replic (split(/\./,$condToObs{$expCondID})) {
#					foreach my $fraction (split(/\+/,$replic)) { # there should be only 1 fraction in label-free...
#						my ($obsID,$pepQuantID,$anaID,$tgPos)=split(':',$fraction);
#						$tgPos=0 unless $tgPos;
#						$anaTg2StatePos{"$anaID.$tgPos"}=$condToState{$expCondID}; # absolute pos in all states list
#$usedTargetPos{$tgPos}=1;
#					}
#				}
#			}
#		}
#		else { # internal
#			my ($anaID)=$dbh->selectrow_array("SELECT ID_ANALYSIS FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
#
#			my ($testChanNum)=$labelingInfo{STATES}[$testStatePos{$selRatioPos}-1]=~/,(\d+)$/;
#			my ($refChanNum)=$labelingInfo{STATES}[$refStatePos{$selRatioPos}-1]=~/,(\d+)$/;
#			$anaTg2StatePos{"$anaID.$testChanNum"}=$testStatePos{$selRatioPos};
#			$anaTg2StatePos{"$anaID.$refChanNum"}=$refStatePos{$selRatioPos};
#$usedTargetPos{$refChanNum}=$usedTargetPos{$testChanNum}=1;
##			my ($parQuantiAnnot,$parQuantifID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,PQ.ID_PARENT_QUANTIFICATION FROM QUANTIFICATION Q,PARENT_QUANTIFICATION PQ WHERE Q.ID_QUANTIFICATION=PQ.ID_PARENT_QUANTIFICATION AND PQ.ID_QUANTIFICATION=$selQuantifID LIMIT 0,1"); # there should be only one
##			my ($labelStrg,@labelInfo)=split('::',$parQuantiAnnot);
##			foreach my $infoStrg (@labelInfo) {
##				last if $infoStrg !~ /^\d/; # entering MassChroQ params if from MassChroQ XIC
##				my ($chanNum,$chanName,$labelDataStrg)=split(';',$infoStrg);
##print "$chanNum,$chanName,$labelDataStrg<BR>\n";
##				#$anaTg2StatePos{"$anaID.$chanNum"}=$testPos;
##
##			}
#			#my ($testPos,$refPos)=split(/\//,$labelingInfo{'RATIOS'}[$selRatioPos-1]);
#			#$anaTg2StatePos{"$anaID.$testPos"}=$testPos;
#			#$anaTg2StatePos{"$anaID.$refPos"}=$refPos;
#		}

		###my $sthXic=$dbh->prepare("SELECT QUANTIF_VALUE,TARGET_POS FROM PEPTIDE_QUANTIFICATION WHERE ID_PEPTIDE=? AND ID_QUANTIFICATION IN ($pepQuantifIDs) LIMIT 0,1");
		#my $sthTgPos=$dbh->prepare("SELECT TARGET_POS FROM PEPTIDE_QUANTIFICATION WHERE ID_PEPTIDE=? AND ID_QUANTIFICATION IN ($pepQuantifIDs)");

		##>Fetching targetPos from quantif files
		my %pepTargetPos;
		foreach my $targetPos (keys %usedTargetPos) {
			my $pepQuantifFile=($targetPos==0)? 'peptide_quantification.txt' : "peptide_quantification_$targetPos.txt";
			open (PEP_QUANTIF,"$promsPath{quantification}/project_$projectID/quanti_$pepQuantifIDs/$pepQuantifFile");
			while(<PEP_QUANTIF>) {
				next if $.==1;
				chomp;
				my ($paramID,$pepID,$quantifValue)=split(/\t/,$_);
				if ($protPeptideData{$pepID}) {push @{$pepTargetPos{$pepID}},$targetPos;}
			}
			close PEP_QUANTIF;
		}
		
		my %targetPos2StatePos;
		open(DATA,"$dataDir/table.txt") || die "$dataDir<BR>ERROR: $!";
		my $prevProtID=0;
		while(<DATA>) {
			next if $.==1;
			chomp;
			my ($mpID,$pepSeq,$varMod,$charge,$pepIdKey,$identifier,$validity,@values)=split(/\t/,$_);
			my ($pID)=($mpID=~/^(\d+)/);
			last if ($pID != $protID && $prevProtID==$protID); # modProtID not always in continuous order (but protID is!)
			$prevProtID=$pID;
			next if $mpID ne $modProtID;
			if ($varMod) {
				if ($varMod=~/^&/) { # var mod code
					$varMod=&promsMod::decodeVarMod($dbh,$pepSeq,$varMod,\%varModName);
				}
				else { # var mod text
					$varMod=~s/^ \+ //;
				}
				#$varMod=' + '.$varMod;
			}
			else {$varMod='';}

			my (%multiPepData,%xicValue,%qSetList,$refAnaID);
			foreach my $statePos (values %anaTg2StatePos) {
				$xicValue{$statePos}=$values[$statePos-1];
			}
			my @peptideList=split(/\s*;\s*/,$pepIdKey);
			unless ($sequenceBeg{$pepSeq}) {
				#$sthPB->execute($peptideList[0]);
				#while (my ($beg)=$sthPB->fetchrow_array) {
				#	push @{$sequenceBeg{$pepSeq}},$beg;
				#}
@{$sequenceBeg{$pepSeq}}=@{$protPeptideData{$peptideList[0]}{BEG}};
			}
			my $uniquePep=$pepSeq.$varMod.'_'.$charge;
			#foreach my $xic (@values) {
			#	$pepOccurence{$selRatioPos}{$uniquePep}++ if $xic ne 'NA';
			#}
			my $basePepCode=$pepSeq.':'.$varMod.'_'.$charge;
			#my @targetList;

			MAIN_PEP:foreach my $pepIdStrg (@peptideList) { # No longer checks
				my $statePos;
				###my $xicValue=0;
my $dataSrc;
				foreach my $pepID (split(/\D/,$pepIdStrg)) { # compatible with XIC summing
					#$sthPD->execute($pepID);
					#my ($score,$pepData,$anaID,$rtApex,$qNum,$rank)=$sthPD->fetchrow_array;
					#$score=0 unless $score;
					#$pepData='' unless $pepData;
my ($score,$pepData,$anaID,$qNum,$rank)=@{$protPeptideData{$pepID}{INFO}};
my $rtApex=$peptideApexRT{$pepID};
					###$sthXic->execute($pepID);
					###my ($xicVal,$tgPos)=$sthXic->fetchrow_array;
					###next unless $xicVal;
					###$tgPos=0 unless $tgPos;
					###next MAIN_PEP unless $anaTg2StatePos{"$anaID.$tgPos"}; # peptide does not belong to test nor ref conditions (multiple ratios in same quantif)
					###$statePos=$anaTg2StatePos{"$anaID.$tgPos"};
					my %matchedTgPos;
					#$sthTgPos->execute($pepID);
					#while (my ($tgPos)=$sthTgPos->fetchrow_array) { # multiple matches for same peptideID for iTRAQ. Only 1 for other (non-)labelling methods
					#	$tgPos=0 unless $tgPos;
					#	$matchedTgPos{$anaTg2StatePos{"$anaID.$tgPos"}}=1 if $anaTg2StatePos{"$anaID.$tgPos"};
					#}
foreach my $tgPos (@{$pepTargetPos{$pepID}}) {
	$matchedTgPos{$anaTg2StatePos{"$anaID.$tgPos"}}=1 if $anaTg2StatePos{"$anaID.$tgPos"};
}
					next MAIN_PEP unless scalar keys %matchedTgPos; # peptide does not belong to test nor ref conditions (multiple ratios in same quantif)
		#$anaList{$anaID}=1;
					$refAnaID=$anaID;
					###$xicValue{$statePos}+=$xicVal; # recompute sum
unless ($dataSrc) {
					$dataSrc='#'.$anaID; # my
					$dataSrc.='.'.$1 if $pepData=~/SOURCE_RANK=(\d+)/;
	#$dataSources{$dataSrc}=$dataSrc unless $dataSources{$dataSrc};	
}
					my $qSet;
					if ($labelType eq 'FREE') {$qSet='-';}
					else {
						$qSet=($pepData=~/QSET=([^#]+)/)? $1 : 0;
						$qSetList{$qSet}=1;
					}
					my $prsStrg=$1 if ($isModifQuantif && $pepData=~/PRS=([^#]+)/);
					#$pepCode.="_$dataSrc" if $labelType ne 'FREE';
					foreach my $statePos (keys %matchedTgPos) {
						@{$multiPepData{$statePos}{$qSet}}=($pepID,$score,$dataSrc,$rtApex,$prsStrg,$qNum,$rank) if (!$multiPepData{$statePos} || !$multiPepData{$statePos}{$qSet} || $multiPepData{$statePos}{$qSet}[1]<$score); # qset seen 1+ for same state!!! => prevents overwrite of best scoring pep in qset (may not be the best fix) PP 01/02/16
#print "+$statePos:$qSet: @{$multiPepData{$statePos}{$qSet}}<BR>\n";
					}
				}
				$pepOccurence{$selRatioPos}{$uniquePep}++;
			}
			my $quantifSet=($labelType eq 'FREE')? '-' : 'A'.$refAnaID.'_'.join('+',sort{$a<=>$b} keys %qSetList);
			my $pepCode=$basePepCode;
			@{$rawPeptideSets{$selRatioPos}{$pepCode}{INFO}}=($sequenceBeg{$pepSeq}[0],$pepSeq,$varMod,$charge,$uniquePep);
				#push @{$rawPeptideSets{$selRatioPos}{$pepCode}{DATA}{$qSet}{$ana2StatePos{$anaID}}},[$score,(sprintf "%.0f",$xicValue),$rtApex,$dataSrc,$prsStrg,$pepID,$qNum,$rank];
			foreach my $statePos (keys %multiPepData) {
				push @{$rawPeptideSets{$selRatioPos}{$pepCode}{DATA}{$quantifSet}{$statePos}},[$xicValue{$statePos},$multiPepData{$statePos}]; # push because of label-free (only 1 index for labeled because of dataSrc in pepCode)
			}
		}
		close DATA;
		###$sthXic->finish;
		#$sthTgPos->finish;
	}
	else { # (Super/Simple)Ratio
		my (%statePos2RatioPos); # %selExp,
		foreach my $ratioPos (keys %protRatio) { # in case different ratios use the same (eg ref) state
			if ($isPepIntensity) { # $labelType eq 'FREE'
				foreach my $statePos (keys %protMean) {@{$statePos2RatioPos{$statePos}}=($ratioPos);}
			}
			else {
##$selExp{chr(64+$ratioPos)}=$ratioPos; # not for iTRAQ (exp='A' always)
				push @{$statePos2RatioPos{$testStatePos{$ratioPos}}},$ratioPos;
				push @{$statePos2RatioPos{$refStatePos{$ratioPos}}},$ratioPos;
			}
			%{$rawPeptideSets{$ratioPos}}=();
		}

		open(DATA,"$dataDir/table.txt") || die "$dataDir<BR>ERROR: $!";
		#my $selExp=chr(64+$selRatioPos); # cannot be 2ndary ratio!!!
		my $prevProtID=0;
		while(<DATA>) {
			next if $.==1;
			chomp;
			my ($mpID,$pepCode,$state,@data)=split(/\t/,$_);
#next unless $usedStates{$state}; # this state is not used for selected ratio
#print "$mpID<BR>\n";
			my ($pID)=($mpID=~/^(\d+)/);
			#last if ($pID != $protID && $prevProtID==$protID); # modProtID not always in continuous order (but protID is!)
			$prevProtID=$pID;
			next if $mpID ne $modProtID;
my $quantifSet=($labelType eq 'FREE')? '-' : $data[-5]; # count from end because +/- techRep
# if ($ratioType eq 'SuperRatio' && $state eq 'State1') { # restrict State1 to selected ratio context (anaIDs)
if (!$isPepIntensity && $state eq 'State1') { # restrict State1 to selected ratio context (anaIDs)
	my ($anaID)=$quantifSet=~/^A(\d+)/;
	next unless $anaList{$anaID};
}
			my ($statePos)=($state=~/(\d+)/);
#if ($labelType eq 'FREE') {next unless $protMean{$statePos};} # state filter
#else {next unless $selExp{$data[-6]};} # Experiment: index from end because +/- techRep
	###next if ($labelType ne 'FREE' && $labelType ne 'ITRAQ' && $data[-6] ne 'NA' && !$selExp{$data[-6]}); # Experiment: index from end because +/- techRep
	###my $quantifSet=($labelType eq 'FREE')? '-' : $data[-5]; # count from end because +/- techRep
			my $pepIdKey=$data[-4]; # count from end because +/- techRep
			my ($pep1ID)=$pepIdKey=~/^(\d+)/;
#next unless $protPeptideData{$pep1ID}; # also filter on states used
			#next unless ($exp eq 'NA' || $exp eq $selExp);
			#my $xicValue=sprintf "%.0f",$data[-1]; # count from end because +/- techRep
			my $xicValue=$data[-1]; # count from end because +/- techRep
			my ($seqModCode,$charge,$dataSrc)=split(/_/,$pepCode);
			$pepCode=~s/_[^_]+\Z//;# if $labelType eq 'FREE'; # remove dataSource for label-free quantif
			my $uniquePep=$seqModCode.'_'.$charge;

			foreach my $ratioPos (@{$statePos2RatioPos{$statePos}}) { # selected ratio(s) that use(s) this state
				$pepOccurence{$ratioPos}{$uniquePep}++;
				unless ($rawPeptideSets{$ratioPos}{$pepCode}) { # 1st time seen
					$dataSrc='-' unless $dataSrc;
					my ($pepSeq,$varMod);
					if ($seqModCode=~/(\w+)&(.+)/) {
						$pepSeq=$1;
						$varMod=&promsMod::decodeVarMod($dbh,$pepSeq,$2,\%varModName);
					}
					else {
						$pepSeq=$seqModCode;
						$varMod='';
					}
					unless ($sequenceBeg{$pepSeq}) {
						#$sthPB->execute($pep1ID);
						#while (my ($beg)=$sthPB->fetchrow_array) {
						#	push @{$sequenceBeg{$pepSeq}},$beg;
						#}
@{$sequenceBeg{$pepSeq}}=@{$protPeptideData{$pep1ID}{BEG}};
					}
					@{$rawPeptideSets{$ratioPos}{$pepCode}{INFO}}=($sequenceBeg{$pepSeq}[0],$pepSeq,$varMod,$charge,$uniquePep);
				}
				$dataSources{$dataSrc}=$dataSrc; # old quantif with true data source file in table.xtx

				my %multiPepData;
				foreach my $pepIdStrg (split(/\+/,$pepIdKey)) {
					my ($pepID)=$pepIdStrg=~/^(\d+)/; # because of early quantifs: pepID1;peppID2;...
					#$sthPD->execute($pepID);
					#my ($score,$pepData,$anaID,$rtApex,$qNum,$rank)=$sthPD->fetchrow_array;
					#$score=0 unless $score;
					#$pepData='' unless $pepData;
unless ($protPeptideData{$pepID}) {print "*$pepID*";}
my ($score,$pepData,$anaID,$qNum,$rank,$validStatus)=@{$protPeptideData{$pepID}{INFO}};
my $rtApex=$peptideApexRT{$pepID};
					#$anaList{$anaID}=1;
					$dataSrc='#'.$anaID;
					$dataSrc.='.'.$1 if $pepData=~/SOURCE_RANK=(\d+)/;
					my $qSet=($labelType eq 'FREE')? '-' : ($pepData=~/QSET=([^#]+)/)? $anaID.':'.$1 : $anaID.':0'; # compatible with merged files?????
					my $prsStrg=$1 if ($isModifQuantif && $pepData=~/PRS=([^#]+)/);
					@{$multiPepData{$qSet}}=($pepID,$score,$dataSrc,$rtApex,$prsStrg,$qNum,$rank,$validStatus);
				}
				#push @{$rawPeptideSets{$ratioPos}{$pepCode}{DATA}{$quantifSet}{$statePos}},[$score,$xicValue,$rtApex,$dataSrc,$prsStrg,$pep1ID,$qNum,$rank]; # push because of label-free (only 1 index for labeled because of dataSrc in pepCode)
				push @{$rawPeptideSets{$ratioPos}{$pepCode}{DATA}{$quantifSet}{$statePos}},[$xicValue,\%multiPepData]; # push because of label-free (only 1 index for labeled because of dataSrc in pepCode)
			}
		}
		close DATA;
	}
	#$sthPB->finish;
	#$sthPD->finish;

	##<Converting dataSources>##
	my $numSources=scalar keys %dataSources;

	my %analysisSrc;
	my $sthDS=$dbh->prepare("SELECT VALID_STATUS,WIFF_FILE FROM ANALYSIS WHERE ID_ANALYSIS=?");
	foreach my $anaID (keys %anaList) {
		$sthDS->execute($anaID);
		my ($validStatus,$wiffFile)=$sthDS->fetchrow_array;
		my $anaDir=($validStatus==2)? "$promsPath{peptide}/proj_$projectID/ana_$anaID" : "$promsPath{valid}/ana_$anaID";
		if (-e "$anaDir/mergedFiles.txt") {
			open(MERGE,"$anaDir/mergedFiles.txt");
			while (<MERGE>) {
				my ($srcRk,$srcName,$srcID)=split(/\t/,$_);
				$analysisSrc{'#'.$anaID.'.'.$srcRk}=$wiffFile.'>'.$srcName;
			}
			close MERGE;
		}
		else {$analysisSrc{'#'.$anaID}=$wiffFile;}
	}
	$sthDS->finish;
	foreach my $srcKey (keys %analysisSrc) {
		$dataSources{$srcKey}=$analysisSrc{$srcKey} if $srcKey=~/^#/;
#print "-$srcKey: $dataSources{$srcKey}<BR>\n";
	}

	$dbh->disconnect;
	$dbhLite->disconnect;

	######################
	####<Starting HTML>###
	######################
	#my ($lightColor,$darkColor)=&promsConfig::getRowColors;
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);
	my ($biasStrg1,$biasStrg2)=($ratioType=~/S\w+Ratio/)? ('','') : ($labelingInfo{'BIAS_CORRECTION'}[0] eq 'TRUE')? ('<SUP>*</SUP>','<I>Intensities after data normalization</I>') : ('','<BR>');
	my $anaStrg=join(',',keys %anaList);

	if ($is2ndaryRatio) { # SuperRatio
		#my $numNormalRatios=0;
		#foreach my $ratioPos (keys %protRatio) {$numNormalRatios++ if $normalRatio{$ratioPos};}
		my $protRatioStrg;
		if ($trueProtRatio < 1) {$protRatioStrg=($numNormalRatios==2)? '1/'.sprintf "%.2f",1/$trueProtRatio : '1/&infin;';} #'∞'
		else {$protRatioStrg=($numNormalRatios==2)? sprintf "%.2f",$trueProtRatio : '&infin;';}
		print qq
|<TABLE><TR><TD nowrap>
<FONT class="title2">Peptide data for <A href="javascript:sequenceView($protID,'$anaStrg')">$protAlias$htmlModRes</A> in $trueTestStateName/$trueRefStateName (ratio=$protRatioStrg):</FONT>
</TD></TR></TABLE>
|;
	}
	my %targetList;
	my $lastPepStart=0;
	foreach my $ratioPos (sort{if ($selRatioPos > 0) {abs($b)<=>abs($a)} else {abs($a)<=>abs($b)}} keys %protRatio) { # test then ref
		my $protRatioStrg;
		if (!$protRatio{$ratioPos}) {$protRatioStrg='unknown <FONT class="title3">[All peptides excluded]</FONT>';} # In case global CVout excludes all peptides of this ratio
		elsif ($protRatio{$ratioPos} < 1) {$protRatioStrg=($normalRatio{$ratioPos})? '1/'.sprintf "%.2f",1/$protRatio{$ratioPos} : '1/&infin;';} #'∞'
		else {$protRatioStrg=($normalRatio{$ratioPos})? sprintf "%.2f",$protRatio{$ratioPos} : '&infin;';}
		if ($is2ndaryRatio) {
			print qq
|<HR width=100%>
<TABLE><TR><TD nowrap>
&nbsp;&nbsp;&nbsp;&bull;<FONT class="title2">Primary ratio: $testStateName{$ratioPos}/$refStateName{$ratioPos} (ratio=$protRatioStrg):</FONT>
|;
		}
		else {
			print qq
|<TABLE><TR><TD nowrap>
<FONT class="title2">Peptide data for <A href="javascript:sequenceView($protID,'$anaStrg')">$protAlias$htmlModRes</A> in $testStateName{$ratioPos}/$refStateName{$ratioPos} (ratio=$protRatioStrg):</FONT>
|;
		}
		print '&nbsp;&nbsp;<FONT class="title3">View:';
		my $rawDataVis='';
		if ($normalRatio{$ratioPos}) {
			$rawDataVis='none';
			print qq
|</FONT><SELECT class="title3" onchange="var divList=['pepPlot:$selQuantifID:$ratioPos','pepList:$selQuantifID:$ratioPos','rawData:$selQuantifID:$ratioPos']; for (var i=0; i<3; i++) {var vis=(divList[i]==this.value)? '' : 'none'; document.getElementById(divList[i]).style.display=vis;}">
<OPTION value="pepPlot:$selQuantifID:$ratioPos">Graphical</OPTION><OPTION value="pepList:$selQuantifID:$ratioPos">List</OPTION><OPTION value="rawData:$selQuantifID:$ratioPos">Raw data</OPTION>
</SELECT>
|;
		}
		else {print " Raw data</FONT>";}
		print qq
|&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'"/>
</TD></TR></TABLE>
<DIV id="rawData:$selQuantifID:$ratioPos" style="display:$rawDataVis">
|;
		my ($numColumns,$scColClass,$srcColStrg)=($isPepIntensity)? (11,'bBorder','') : ($numSources > 1)? (12,'rbBorder','<TH class="bBorder">&nbsp;Source&nbsp;</TH>') : (11,'bBorder',''); # $labelType eq 'FREE'
		#my $rtColStrg=($labelType eq 'FREE')? '' : "<TH class=\"$rtColClass\" nowrap>&nbsp;Ret. time&nbsp;</TH>";
		#if ($labelType eq 'ITRAQ') { # no need for RT column
		#	$rtColStrg='';
		#	$numColumns--;
		#}
		if ($ratioType eq 'Ratio' && !scalar keys %{$rawPeptideSets{$ratioPos}}) { # TODO: temp until update
			print "<BR>&nbsp;&nbsp;<FONT class=\"title3\">This quantification is too old for infinite ratios display.<BR>&nbsp;&nbsp;Delete it and run it again if necessary.<BR><BR><BR>\n";
		}
		else {
			my $scLabelStrg=($labelType=~/ITRAQ|TMT/)? 'Score' : 'Scores';
			my $hasGhosts=0;
			print qq
|<TABLE cellspacing=0>
<TR bgcolor=$darkColor>
<TH class="rbBorder">#</TH><TH class="rbBorder">&nbsp;Pos.&nbsp;<TH class="rbBorder">Peptides</TH><TH class="rbBorder">&nbsp;Charge&nbsp;</TH><TH class="rbBorder">&nbsp;Occur.&nbsp;</TH><TH class="rbBorder" colspan=3>&nbsp;XIC values&nbsp;</TH><TH class="$scColClass" colspan=3>&nbsp;$scLabelStrg&nbsp;</TH>$srcColStrg
</TR>
|;
			my $bgColor=$lightColor;
			my $numPep=0;
			my $peptide;
			my $prevPeptide='';
			my $prevUniquePep='';
			my $prevBeg=0;
			foreach my $pepCode (sort{$rawPeptideSets{$ratioPos}{$a}{INFO}[0]<=>$rawPeptideSets{$ratioPos}{$b}{INFO}[0] || $rawPeptideSets{$ratioPos}{$a}{INFO}[1] cmp $rawPeptideSets{$ratioPos}{$b}{INFO}[1] || $rawPeptideSets{$ratioPos}{$a}{INFO}[2] cmp $rawPeptideSets{$ratioPos}{$b}{INFO}[2] || $rawPeptideSets{$ratioPos}{$a}{INFO}[3]<=>$rawPeptideSets{$ratioPos}{$b}{INFO}[3]} keys %{$rawPeptideSets{$ratioPos}}) {
				my $varMod=($rawPeptideSets{$ratioPos}{$pepCode}{INFO}[2])? ' + '.$rawPeptideSets{$ratioPos}{$pepCode}{INFO}[2] : '';
				my $uniquePep=$rawPeptideSets{$ratioPos}{$pepCode}{INFO}[4];
				my $numPepStrg='';
				$peptide=$rawPeptideSets{$ratioPos}{$pepCode}{INFO}[1].$varMod;
				my @residues=($usePtmProb)? split('',$rawPeptideSets{$ratioPos}{$pepCode}{INFO}[1]) : (); # pep seq
				foreach my $quantifSet (sort{&promsMod::sortSmart($rawPeptideSets{$ratioPos}{$pepCode}{DATA}{$a}{$testStatePos{$ratioPos}}[0][3],$rawPeptideSets{$ratioPos}{$pepCode}{DATA}{$b}{$testStatePos{$ratioPos}}[0][3])} keys %{$rawPeptideSets{$ratioPos}{$pepCode}{DATA}}) {
					if ($dispPepType eq 'distinct') {
						if ($peptide ne $prevPeptide) {
							$numPep++;
							$numPepStrg=$numPep;
						}
					}
					else {
						$numPep++;
						$numPepStrg=$numPep;
					}
					my $peptideStrg=($peptide ne $prevPeptide)? $rawPeptideSets{$ratioPos}{$pepCode}{INFO}[1].$varMod : '';
					my $begStrg=($prevBeg != $rawPeptideSets{$ratioPos}{$pepCode}{INFO}[0])? "$rawPeptideSets{$ratioPos}{$pepCode}{INFO}[0]&nbsp;" : '';
					my ($chargeStrg,$occStrg)=($prevUniquePep eq $uniquePep)? ('','') : ($rawPeptideSets{$ratioPos}{$pepCode}{INFO}[3].'<SUP>+</SUP>','x'.$pepOccurence{$ratioPos}{$uniquePep});
					print qq
|<TR bgcolor="$bgColor" valign="top">
	<TD class="rBorder" align="right">&nbsp;$numPepStrg&nbsp;</TD><TD align="right" valign="top">$begStrg</TD><TH class="font11" align="left" valign="top" nowrap>$peptideStrg&nbsp;</TH>
	<TH valign="top">$chargeStrg</TH><TD align="center" valign="top">&nbsp;$occStrg&nbsp;</TD>
|;
					#my ($sumRetTime,$numRetTime)=(0,0); # labeled
					my %numPeps=($testStatePos{$ratioPos}=>0,$refStatePos{$ratioPos}=>0); # FREE
					my %numPeps2=($testStatePos{$ratioPos}=>0,$refStatePos{$ratioPos}=>0);
					my (%peptID,%sc,%dtsrc,%rt,%prs,%qNum,%rank);
					foreach my $statePos ($testStatePos{$ratioPos},$refStatePos{$ratioPos}) {
						#my ($xicValue,$retTime)=($rawPeptideSets{$pepCode}{DATA}{$statePos} && $rawPeptideSets{$pepCode}{DATA}{$statePos}[1])? ($rawPeptideSets{$pepCode}{DATA}{$statePos}[1],$rawPeptideSets{$pepCode}{DATA}{$statePos}[2]) : ('-','-');
						my $xicStrg='';
						#my $prsStrg; # labeled
#print "*$statePos: ",scalar @{$rawPeptideSets{$ratioPos}{$pepCode}{DATA}{$statePos}},"<BR>\n";
						foreach my $refData (@{$rawPeptideSets{$ratioPos}{$pepCode}{DATA}{$quantifSet}{$statePos}}) {
							
							#>Multi-peptide data
							my $isTruePep=0;
							foreach my $qSet (keys %{$refData->[1]}) {
#print "*$qSet: @{$refData->[1]{$qSet}}<BR>\n";
								#($peptID{$qSet}{$statePos},$sc{$qSet}{$statePos},$rt{$qSet}{$statePos},$prs{$qSet}{$statePos},$qNum{$qSet}{$statePos},$rank{$qSet}{$statePos})=@{$refData->[2]{$qSet}};
								push @{$peptID{$qSet}{$statePos}},$refData->[1]{$qSet}[0]; # push because of new Label-Free (but only 1 qSet per dataScr)
								push @{$sc{$qSet}{$statePos}},$refData->[1]{$qSet}[1];
								push @{$dtsrc{$qSet}{$statePos}},$refData->[1]{$qSet}[2];
								push @{$rt{$qSet}{$statePos}},$refData->[1]{$qSet}[3];
								push @{$prs{$qSet}{$statePos}},$refData->[1]{$qSet}[4];
								push @{$qNum{$qSet}{$statePos}},$refData->[1]{$qSet}[5];
								push @{$rank{$qSet}{$statePos}},$refData->[1]{$qSet}[6];
								$isTruePep=1 if ($refData->[1]{$qSet}[7] || $refData->[1]{$qSet}[1]);
								$hasGhosts=1 if (!$hasGhosts && !$refData->[1]{$qSet}[7] && !$refData->[1]{$qSet}[1]);
								$numPeps2{$statePos}++;
							}
							my ($ghostFlag1,$ghostFlag2)=($isTruePep)? ('','') : ('<SPAN class="virtualProt">','</SPAN>');
							if ($isPepIntensity || $labelType eq 'FREE') { # $labelType eq 'FREE'
								$xicStrg.='<BR>' if $xicStrg;
					#my $rt=($refData->[2])? sprintf "%.2f",$refData->[2] : '-';
					#my $dtSrc=($refData->[0])? $dataSources{$refData->[0]} : '-';
								#$xicStrg.="&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<TABLE cellspacing=0><TR><TD align=right nowrap><U>Ret. time:</U></TD><TD>$rt<TD></TR><TR><TD align=right><U>Source:</U></TD><TD>$dtSrc</TD></TR></TABLE>')\" onmouseout=\"popout()\">$refData->[1]</A>&nbsp;";
								$xicStrg.=($refData->[0] && $refData->[0] ne 'NA')? sprintf "&nbsp;$ghostFlag1%.2e$ghostFlag2&nbsp;",$refData->[0] : '-';
								$numPeps{$statePos}++;
							}
							else {
								#if ($refData->[2]) {
								#	$sumRetTime+=$refData->[2];
								#	$numRetTime++;
								#}
								#($dataSrc,$prsStrg)=@{$refData}[3..4];
								#$xicStrg=sprintf "&nbsp;%1.2e&nbsp;",$refData->[1] if $refData->[1];
								if ($refData->[0] && $refData->[0] ne 'NA') { # NA if algo Ratio
									$xicStrg=sprintf "&nbsp;$ghostFlag1%.2e$ghostFlag2&nbsp;",$refData->[0];
									#if (!$refData->[6]) {$xicStrg="&nbsp;<I>$xicValue</I>";} # !query number -> ghost
									#else {
									#	$xicStrg="&nbsp;<A id=\"pep_$refData->[5]\" href=\"javascript:drawSpectrum('pep_$refData->[5]','pep_$refData->[5]_$refData->[6]_$refData->[7]')\" onmouseover=\"popup('Click to display <B>Fragmentation Spectrum</B>.')\" onmouseout=\"popout()\">$xicValue</A>";
									#}
									#if ($prsStrg) {$xicStrg.=&phosphoRS::printIcon($prsStrg,{format=>'string'});}
									#else {$xicStrg.='&nbsp;'}
								}
							}
							
						}		
						$xicStrg='&nbsp;-&nbsp;' unless $xicStrg;
						my $sepStrg;
						if ($isPepIntensity || $labelType eq 'FREE') { # $labelType eq 'FREE' for 'Ratio'
							my $numSep=(sort{$b<=>$a} values %numPeps)[0];
							$sepStrg='|<BR>' x $numSep;
						}
						else {$sepStrg='/';}
						#$retTimeStrg.=($retTime)? sprintf "%.2f",$retTime : '-';
						#$retTimeStrg.=($retTime)? &phosphoRS::printIcon($prsStrg,{format=>'string',text=>'&nbsp;'.(sprintf "%.2f",$retTime).'&nbsp;'}) : '-';
						#$retTimeStrg.=&phosphoRS::printIcon($prsStrg,{format=>'string'}) if $prsStrg;
						my $align;
						if ($statePos==$testStatePos{$ratioPos}) {
							$align='right';
							#$retTimeStrg.=($ratioType eq 'Ratio')? '-' : ' / '; # labeled
							#$retTimeStrg.='/';
						}
						else {
							$align='left';
							print "<TD>$sepStrg</TD>";
						}
						print "<TD align=\"$align\">$xicStrg</TD>";
					}

					#>Other columns
					my $dataSrc='-';
					my $align='right';
					my $numSep2=(sort{$b<=>$a} values %numPeps2)[0];
					foreach my $statePos ($testStatePos{$ratioPos},$refStatePos{$ratioPos}) {
						my $scoreStrg;
						my $plusStrg='';
						foreach my $qSet (sort keys %peptID) {
							if ($peptID{$qSet}{$statePos}) {
								foreach my $i (0..$#{$peptID{$qSet}{$statePos}}) { # i > 0 for Label-Free only!
									$scoreStrg.="<BR>" if $scoreStrg;
									($dataSrc,my $dtSrcStrg)=($dtsrc{$qSet}{$statePos}[$i])? ($dtsrc{$qSet}{$statePos}[$i],$dataSources{$dtsrc{$qSet}{$statePos}[$i]}) : ('-','-'); # changes if label-free (unique for labeled)
									if ((!$dtSrcStrg || $dtSrcStrg eq $dataSrc) && $dataSrc=~/(#\d+)\.(\d+)/) { # for Algo v1
										$dtSrcStrg=$dataSources{$dataSrc}=$dataSources{$1}.'>file#'.$2;
									}
									my $dataPopStrg='';
									$dataPopStrg='<TABLE cellspacing=0>';
									if ($rt{$qSet}{$statePos}[$i]) {
										my $ret=1*(sprintf "%.2f",$rt{$qSet}{$statePos}[$i]);
										$dataPopStrg.="<TR><TD align=right nowrap><U>Ret. time:</U></TD><TD>$ret min.</TD></TR>";
									}
									$dataPopStrg.="<TR><TD align=right><U>Source:</U></TD><TD>$dtSrcStrg</TD></TR></TABLE>";
									my $modProbPopStrg="";
									if ($usePtmProb && $ptmProb{$peptID{$qSet}{$statePos}[$i]}) {
										foreach my $modID (sort{$quantifModifInfo{ID2RANK}{$a}<=>$quantifModifInfo{ID2RANK}{$b}} keys %{$ptmProb{$peptID{$qSet}{$statePos}[$i]}}) {
											$modProbPopStrg.="<br/><U>$ptmSoft&nbsp;$quantifModifInfo{NAME}{$modID}-site&nbsp;probabilities:</U><BR>";
											foreach my $refPos (@{$ptmProb{$peptID{$qSet}{$statePos}[$i]}{$modID}}) {
												my $siteStrg=($refPos->[0]=~/\d/)? $residues[$refPos->[0]-1].$refPos->[0] : ($refPos->[0] eq '-')? 'Protein N-term' : ($refPos->[0] eq '=')? 'N-term' : ($refPos->[0] eq '+')? 'Protein C-term' : ($refPos->[0] eq '*')? 'C-term' : $refPos->[0];
												$modProbPopStrg.="&nbsp;&nbsp;-".$siteStrg.":".&promsMod::PTMProbIcon($refPos->[1],{text=>'&nbsp;'.($refPos->[1]*100).'%&nbsp;',inPopup=>1})."<BR>";
											}
										}
										$modProbPopStrg .= "<br/>" if($modProbPopStrg);
									}
									my ($score,$linkStrg,$linkPopStrg)=('?','void(null)','');
									if ($sc{$qSet}{$statePos}[$i]) {
										$score=1 * (sprintf '%.1f',$sc{$qSet}{$statePos}[$i]);
										if ($existSpectrumFile) {
											$linkStrg="drawSpectrum('pep_$peptID{$qSet}{$statePos}[$i]','pep_$peptID{$qSet}{$statePos}[$i]_$qNum{$qSet}{$statePos}[$i]_$rank{$qSet}{$statePos}[$i]')";
											$linkPopStrg='Click&nbsp;to&nbsp;display&nbsp;<B>Fragmentation&nbsp;Spectrum</B>';
										}
									}
									
									my $popupContent = "$dataPopStrg$modProbPopStrg$linkPopStrg";
									$popupContent =~ s/([^\\])'/$1\\'/ig;
									$scoreStrg.="&nbsp;$plusStrg<A id=\"pep_$peptID{$qSet}{$statePos}[$i]\" href=\"javascript:$linkStrg\" onmouseover=\"popup('$popupContent')\" onmouseout=\"popout()\">$score</A>";
									$scoreStrg.=($prs{$qSet}{$statePos}[$i])? &phosphoRS::printIcon($prs{$qSet}{$statePos}[$i],{format=>'string'}) : '&nbsp;';
								}
							}
							else {
								$scoreStrg.="<BR>" if $scoreStrg;
								$scoreStrg.='&nbsp;-&nbsp;';
							}
							$plusStrg='+';
						}
						if ($labelType=~/ITRAQ|TMT/) {
							print "<TD colspan=3 align=\"center\">$scoreStrg</TD>"; # same score for all states
							last;
						}
						else {
							print "<TD align=\"$align\">$scoreStrg</TD>";
							if ($statePos==$testStatePos{$ratioPos}) {
								my $sepStrg=($isPepIntensity || $labelType eq 'FREE')? '|<BR>' : '/<BR>'; # $labelType eq 'FREE' for Ratio
								print '<TD>',$sepStrg x $numSep2,'</TD>';
							}
							$align='left';
						}
					}

					print "<TD>&nbsp;$dataSources{$dataSrc}&nbsp;</TD>" if (!$isPepIntensity && $numSources > 1); # $labelType ne 'FREE'
					print "</TR>\n";
					$prevBeg=$rawPeptideSets{$ratioPos}{$pepCode}{INFO}[0];
					$prevPeptide=$peptide;
					$prevUniquePep=$uniquePep;
					$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
				}
			}
			print "<TR bgcolor=\"$darkColor\"><TD colspan=$numColumns></TD></TR>\n" if $bgColor eq $darkColor;
			print "</TABLE>\n";
			print qq |<SPAN class="font11" style="font-weight:bold">Peptide ions identified by MSMS, <FONT class="virtualProt">Peptide ions rescued by Match Between or Within Runs (MBWR).</FONT></SPAN>| if $hasGhosts;
			print "<BR><BR>\n";
		}
		print "</DIV>\n";

		####<Normal ratio>####
		if ($normalRatio{$ratioPos}) {
			my ($refClass,$sourceColStrg,$colspan)=($numSources > 1)? ('rbBorder','<TH class="bBorder">&nbsp;Source&nbsp;</TH>',8) : ('bBorder','',7);
			my ($labelFreeDivStrg,%titleStrg,$valueHeader);
			if ($isPepIntensity) { # $labelType eq 'FREE' && $ratioType ne 'Ratio'
				$valueHeader='Mean';
				foreach my $targetPos (sort{if ($selRatioPos > 0) {$b<=>$a} else {$a<=>$b}} keys %protMean) {
					push @{$targetList{$ratioPos}},$targetPos;
					my $stateName=($labelFreeDivStrg)? $trueRefStateName : $trueTestStateName;
					my $disProtValue=sprintf "%.3e",$protMean{$targetPos};
					$titleStrg{$targetPos}="&nbsp;&nbsp;&nbsp;&bull;<FONT class=\"title2\">Mean values for $stateName: $disProtValue</FONT><BR>";
					$labelFreeDivStrg.="$titleStrg{$targetPos}<DIV id=\"pepPlot:$selQuantifID:$ratioPos:$targetPos\"></DIV>\n";
				}
			}
			else {
				$valueHeader='Ratio';
				@{$targetList{$ratioPos}}=($ratioPos);
				$labelFreeDivStrg=$titleStrg{$ratioPos}='';
			}
			print qq
|<DIV id="pepPlot:$selQuantifID:$ratioPos">$labelFreeDivStrg</DIV>
<DIV id="pepList:$selQuantifID:$ratioPos" style="display:none">
|;
			#foreach my $targetPos (@{$targetList{$ratioPos}}) { #}
			foreach my $tpIdx (0..$#{$targetList{$ratioPos}}) {
				my $targetPos=$targetList{$ratioPos}[$tpIdx];
				print qq
|$titleStrg{$targetPos}
<TABLE cellspacing=0>
<TR bgcolor=$darkColor>
<TH class="rbBorder">#</TH><TH class="rbBorder">&nbsp;Pos.&nbsp;<TH class="rbBorder">Peptides</TH><TH class="rbBorder">&nbsp;Charge&nbsp;</TH><TH class="rbBorder">&nbsp;$valueHeader&nbsp;</TH>
|;
				if ($ratioType eq 'Ratio') {$colspan++; print "<TH class=\"rbBorder\">&nbsp;$testStateName{$targetPos}$biasStrg1&nbsp;</TH><TH class=\"rbBorder\">&nbsp;$refStateName{$targetPos}$biasStrg1&nbsp;</TH>";}
				else {print "<TH class=\"rbBorder\">&nbsp;Occur.&nbsp;</TH>";}
				print "<TH class=\"$refClass\">&nbsp;Scores&nbsp;</TH>$sourceColStrg\n";
				print "</TR>\n";
				my $bgColor=$lightColor;
				my $numPep=0;
				my $prevPeptide='';
				my $prevBeg=0;
				my $maxBeg=0;
				foreach my $pepIdKey (sort{$peptideSets{$targetPos}{$a}{INFO}[0]<=>$peptideSets{$targetPos}{$b}{INFO}[0] || $peptideSets{$targetPos}{$a}{INFO}[1] cmp $peptideSets{$targetPos}{$b}{INFO}[1] || $peptideSets{$targetPos}{$a}{INFO}[2] cmp $peptideSets{$targetPos}{$b}{INFO}[2] || $peptideSets{$targetPos}{$a}{INFO}[3]<=>$peptideSets{$targetPos}{$b}{INFO}[3]} keys %{$peptideSets{$targetPos}}) { # beg,seq,vMod,charge
					$maxBeg=$sequenceBeg{$peptideSets{$targetPos}{$pepIdKey}{INFO}[1]}[-1];
					my $varMod=($peptideSets{$targetPos}{$pepIdKey}{INFO}[2])? " + $peptideSets{$targetPos}{$pepIdKey}{INFO}[2]" : '';
					my $numPepStrg='';
					if ($dispPepType eq 'distinct') {
						my $peptide="$peptideSets{$targetPos}{$pepIdKey}{INFO}[1]$varMod";
						if ($peptide ne $prevPeptide) {
							$numPep++;
							$numPepStrg=$numPep;
							$prevPeptide=$peptide;
						}
					}
					else {
						$numPep++;
						$numPepStrg=$numPep;
					}
					my $begStrg=($prevBeg != $peptideSets{$targetPos}{$pepIdKey}{INFO}[0])? "$peptideSets{$targetPos}{$pepIdKey}{INFO}[0]&nbsp;" : '';
					my $flagImgStrg=($peptideSets{$targetPos}{$pepIdKey}{DATA}[-1] eq 'NA')? '' : ($peptideSets{$targetPos}{$pepIdKey}{DATA}[-1]==1)? "<IMG src=\"$promsPath{images}/good.gif\">" : "<IMG src=\"$promsPath{images}/bad.gif\">";
					#my $ratio=($labelType eq 'FREE' && $ratioType eq 'SimpleRatio')? sprintf "%.0f",$peptideSets{$targetPos}{$pepIdKey}{DATA}[0] : ($peptideSets{$targetPos}{$pepIdKey}{DATA}[0]>=1)? sprintf "%.2f",$peptideSets{$targetPos}{$pepIdKey}{DATA}[0] : '1/'.sprintf "%.2f",1/$peptideSets{$targetPos}{$pepIdKey}{DATA}[0];
					my $ratio=($isPepIntensity)? sprintf "%.0f",$peptideSets{$targetPos}{$pepIdKey}{DATA}[0] : ($peptideSets{$targetPos}{$pepIdKey}{DATA}[0]>=1)? sprintf "%.2f",$peptideSets{$targetPos}{$pepIdKey}{DATA}[0] : '1/'.sprintf "%.2f",1/$peptideSets{$targetPos}{$pepIdKey}{DATA}[0];
					print qq
|<TR bgcolor="$bgColor">
<TD class="rBorder" align=right>&nbsp;$numPepStrg&nbsp;</TD><TD align=right>$begStrg</TD><TH class="font11" align=left nowrap>$flagImgStrg$peptideSets{$targetPos}{$pepIdKey}{INFO}[1]$varMod&nbsp;</TH>
<TH>$peptideSets{$targetPos}{$pepIdKey}{INFO}[3]<SUP>+</SUP></TH>
<TD align=center>&nbsp;$ratio&nbsp;</TD>
|;
					if ($ratioType eq 'Ratio') {
						my $testVal=sprintf "%.2e",$peptideSets{$targetPos}{$pepIdKey}{DATA}[1];
						my $refVal=sprintf "%.2e",$peptideSets{$targetPos}{$pepIdKey}{DATA}[2];
						print "<TD align=center>&nbsp;$testVal&nbsp;</TD><TD align=center>&nbsp;$refVal&nbsp;</TD>\n";
					}
					else {print "<TD align=center>&nbsp;x$peptideSets{$targetPos}{$pepIdKey}{DATA}[1]&nbsp;</TD>";} # (Super/Simple)Ratio: pep set occurence
					print "<TD align=center>&nbsp;$peptideSets{$targetPos}{$pepIdKey}{INFO}[4]&nbsp;</TD>\n";
					if ($numSources > 1) {
						my @srcFiles;
						foreach my $srcKey (@{$peptideSets{$targetPos}{$pepIdKey}{INFO}[5]{DATA}}) {push @srcFiles,$dataSources{$srcKey} || $srcKey;}
						if ($ratioType eq 'Ratio') {print "<TD>&nbsp;",join('<BR>',sort @srcFiles),"&nbsp;</TD>";} # $dataSources{$peptideSets{$targetPos}{$pepIdKey}{INFO}[5]}
						else {
							my $srcName=($algoVersion>=3)? $peptideSets{$targetPos}{$pepIdKey}{INFO}[5]{NAME} : '...';
							print "<TD>&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<U>Source(s):</U><BR>",join('<BR>',sort @srcFiles),"')\" onmouseout=\"popout()\">$srcName</A>&nbsp;</TD>";
						} # $peptideSets{$targetPos}{$pepIdKey}{INFO}[5]
					}
					print "</TR>\n";
					$prevBeg=$peptideSets{$targetPos}{$pepIdKey}{INFO}[0];
					$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
				}
				print "<TR bgcolor=\"$lightColor\"><TD colspan=\"$colspan\">&nbsp;<B>No peptide data available</B></TD></TR>\n" unless $numPep;
				print "<TR bgcolor=\"$darkColor\"><TD colspan=\"$colspan\"></TD></TR>\n" if $bgColor eq $darkColor;
				$lastPepStart=$maxBeg if $maxBeg > $lastPepStart;
				print "</TABLE>\n";
				print "<BR>\n" if $tpIdx < $#{$targetList{$ratioPos}};
			}
			print qq
|$biasStrg1$biasStrg2
</DIV><BR>
|;
		}
	}

	exit unless $numNormalRatios; # !!! CONDITIONAL EXIT !!!

	my $plotWidth=($protLength)? $protLength : $lastPepStart+50;
	$plotWidth=400 if $plotWidth < 400;
	$plotWidth=1000 if $plotWidth > 1000;
	#my $plotReferenceFlag=($labelType eq 'FREE' && $ratioType ne 'Ratio')? 'false' : 'true';

	###<javascript for peptidePlot>###
	print "#==========#\n";
	foreach my $ratioPos (sort{if ($selRatioPos > 0) {abs($b)<=>abs($a)} else {abs($a)<=>abs($b)}} keys %protRatio) { # test then ref
		next unless $protRatio{$ratioPos}; # Undefined if all peptides were excluded by global CVout
		if ($normalRatio{$ratioPos}) {
			foreach my $targetPos (@{$targetList{$ratioPos}}) {
				my $trueTargetPos=abs($targetPos);
				my ($extraIdStrg,$protValueType,$protValue)=($isPepIntensity)? (':'.$targetPos,'protMean',$protMean{$targetPos}) : ('','protRatio',$protRatio{$ratioPos}); # $labelType eq 'FREE' && $ratioType ne 'Ratio'
				if (($isPepIntensity || $labelType eq 'FREE') && !defined($protValue)) { # no peptide at all for this state
					print qq
|document.getElementById('pepPlot:$selQuantifID:$ratioPos$extraIdStrg').innerHTML='&nbsp;&nbsp;&nbsp;<FONT class="title3">No peptide data available</FONT>'
|;
					next;
				}
				if ($isPepIntensity) { # $labelType eq 'FREE' && $ratioType ne 'Ratio'
					my $ratioPropStrg=($algoVersion>=3)? ",'Source'" : '';
					my $pepColorStrg=($algoVersion>=3)? "peptideColor: {map:'Source',type:'discrete'}" : '';
					print qq
|var PP_$trueTargetPos=new peptidePlot({div:'pepPlot:$selQuantifID:$ratioPos$extraIdStrg',width:$plotWidth,height:300,valueAxisLabel:'Log2(mean XIC)',valueLabel:'Mean XIC',
									protein: {name:'$protAlias$textModRes',length:$protLength,value:$protValue,relativeThresholds:[0.5,2]},
									peptideProperties: ['Score'$ratioPropStrg],
									convertValue: function(v) {return Math.log(v)/$log2},
									convertValueDisplayed: function(v) {return (v*1).toPrecision(3)},
									//convertValueDisplayed: function(v) {return Math.round(v*1)},
									$pepColorStrg
									});
|;
				}
				else {
					my $ratioPropStrg=($ratioType eq 'Ratio')? ",'Test','Source'" : ($algoVersion>=3)? ",'Source'" : '';
					my $pepColorStrg=($ratioType eq 'Ratio' || $algoVersion>=3)? "peptideColor: {map:'Source',type:'discrete'}" : '';
					print qq
|var PP_$trueTargetPos=new peptidePlot({div:'pepPlot:$selQuantifID:$ratioPos$extraIdStrg',width:$plotWidth,height:300,valueAxisLabel:'Log2(ratio)',valueLabel:'Ratio',
									protein: {name:'$protAlias$textModRes',length:$protLength,value:$protValue,relativeThresholds:[0.5,2]},
									reference: {label:'Reference: No fold change',value:1},
									peptideProperties: ['Score'$ratioPropStrg],
									convertValue: function(v) {return Math.log(v)/$log2},
									convertValueDisplayed: function(v) {return (v<1)? '1/'+(1/v).toPrecision(3) : (v*1).toPrecision(3)},
									$pepColorStrg
									});
|;
				}
				foreach my $pepIdKey (sort{$peptideSets{$targetPos}{$a}{INFO}[0]<=>$peptideSets{$targetPos}{$b}{INFO}[0] || $peptideSets{$targetPos}{$a}{INFO}[3]<=>$peptideSets{$targetPos}{$b}{INFO}[3]} keys %{$peptideSets{$targetPos}}) { # beg,charge
					my $vModStrg=($peptideSets{$targetPos}{$pepIdKey}{INFO}[2])? "'$peptideSets{$targetPos}{$pepIdKey}{INFO}[2]'" : 'null';
					my $excluded=($peptideSets{$targetPos}{$pepIdKey}{DATA}[-1] eq '0')? 1 : 0; # 'eq' because 'NA' possible
					print "PP_$trueTargetPos.addPeptide([[",join(',',@{$sequenceBeg{$peptideSets{$targetPos}{$pepIdKey}{INFO}[1]}}),"],'$peptideSets{$targetPos}{$pepIdKey}{INFO}[1]',$vModStrg,$peptideSets{$targetPos}{$pepIdKey}{INFO}[3],$peptideSets{$targetPos}{$pepIdKey}{DATA}[0],$excluded,'$peptideSets{$targetPos}{$pepIdKey}{INFO}[4]'"; # [beg],seq,mod,charge,ratio,exclu,score
					# source
					if ($algoVersion>=3) {
						print ",'$peptideSets{$targetPos}{$pepIdKey}{INFO}[5]{NAME}'";
					}
					elsif ($ratioType eq 'Ratio') {
						my $testVal=sprintf "%.2e",$peptideSets{$targetPos}{$pepIdKey}{DATA}[1];
						my $refVal=sprintf "%.2e",$peptideSets{$targetPos}{$pepIdKey}{DATA}[2];
						print ",'$testVal/$refVal','",$dataSources{$peptideSets{$targetPos}{$pepIdKey}{INFO}[5]{DATA}->[0]},"'";
					}
					print "]);\n";
				}
				print "PP_$trueTargetPos.draw();\n";
			}
		}
	}
	exit;
}


####################<<<ajaxPeptideMaxQuant>>>#####################
sub ajaxPeptideMaxQuant { # Data display not yet implemented
	my $modProtID=param('id_prot');
	my ($protID,$modStrg)=($modProtID=~/^(\d+)-?(.*)/);
	my $dispMeasure=param('dispMeasure') || '';
	my $selTargetPos=param('ratio');
	my $trueTargetPos=abs($selTargetPos);

	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;
	my $projectID=&promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION');
	my $dbhLite=&promsQuantif::dbConnectProteinQuantification($selQuantifID,$projectID);

	#my ($quantifAnnot,$designID,$quantifType,$quantifModID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_DESIGN,M.CODE,Q.ID_MODIFICATION FROM QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$selQuantifID");
	my ($quantifAnnot,$designID,$quantifMethodID,$quantifType,$quantifModID,$multiModifStrg)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_DESIGN,M.ID_QUANTIFICATION_METHOD,M.CODE,Q.ID_MODIFICATION,GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ',')
																	 FROM QUANTIFICATION Q
																	 LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
																	 INNER JOIN QUANTIFICATION_METHOD M ON Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD
																	 WHERE Q.ID_QUANTIFICATION=$selQuantifID GROUP BY Q.ID_QUANTIFICATION");
	#my $isModifQuantif=($quantifModID || $multiModifStrg)? 1 : 0;
	#my $isMultiModifQuantif=($multiModifStrg)? 1 : 0;
	my %quantifModifInfo;
	my ($isModifQuantif,$isMultiModifQuantif)=&promsQuantif::getQuantifModificationInfo($dbh,$quantifModID,$multiModifStrg,\%quantifModifInfo);
	my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
	my ($labelType)=($labelStrg)? $labelStrg=~/LABEL=(.+)/ : ('FREE');
	$labelType=uc($labelType); # iTRAQ to ITRAQ
	my (%labelingInfo,%condToState,%condToObs);
	$labelingInfo{'LABELTYPE'}=$labelType;
	#if ($designID || $labelType =~ /SILAC|ITRAQ/) {
		foreach my $infoStrg (@labelInfo) {
			my ($setting,$valueStrg)=split('=',$infoStrg);
			$valueStrg=~s/#//g; # remove all ID tags
			@{$labelingInfo{$setting}}=split(';',$valueStrg);
		}
	#}
	my $ratioType=$labelingInfo{RATIO_TYPE}[0];

	if ($modStrg && !$isModifQuantif) { # fall back to whole protID if called with a modProtID for a non-modif quantif (eg. called from log2 plot)
		$modProtID=$protID;
		$modStrg='';
	}
	my ($protAlias,$protLength)=$dbh->selectrow_array("SELECT ALIAS,PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=$protID");
	$protLength=0 unless $protLength;

	my $paramCode=$dispMeasure || 'RATIO';
	my ($paramID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER,NAME,CODE FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=$quantifMethodID AND CODE='$paramCode'");

	my ($displayModStrg,$siteQueryStrg)=('','');
	if ($modStrg) { # modification quantif
		my $fileModStrg=&promsQuantif::siteCode2QuantifFileCode($modStrg,$isMultiModifQuantif,$quantifModifInfo{ID2RANK});
		$displayModStrg=&promsQuantif::displayModificationSites($modStrg,$quantifModifInfo{DISPLAY},'html');
		$modProtID=$protID.'-'.$fileModStrg; # redefine modProtID to match quantif file format !!!!!!!!!!!!!!!!!!!!!!!!
		$siteQueryStrg="AND SITE_CODE='$fileModStrg'";
	}
	my $sthPR=$dbhLite->prepare("SELECT QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIF_PARAMETER=$paramID AND ID_PROTEIN=$protID AND TARGET_POS=$trueTargetPos $siteQueryStrg");
	$sthPR->execute;
	my ($protValue)=$sthPR->fetchrow_array;
	$sthPR->finish;

	#########################
	###<Processing design>###
	#########################
	my ($quantifTargetName,$testStateID,$refStateID,$selStateID);
	my $sthExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
	if ($ratioType eq 'None') {
		(my $numBioRep,my $quantiObsIDs,$selStateID)=split(',',$labelingInfo{STATES}[$selTargetPos-1]);
		($quantifTargetName)=$sthExpCondName->execute($selStateID);
	}
	else {
		($testStateID,$refStateID)=($selTargetPos > 0)? split(/\//,$labelingInfo{RATIOS}[$trueTargetPos-1]) : reverse (split(/\//,$labelingInfo{RATIOS}[$trueTargetPos-1]));
		$sthExpCondName->execute($testStateID);
		my ($testStateName)=$sthExpCondName->fetchrow_array;
		$sthExpCondName->execute($refStateID);
		my ($refStateName)=$sthExpCondName->fetchrow_array;
		$quantifTargetName="$testStateName/$refStateName";
	}
	$sthExpCondName->finish;

	my (%anaList); #,%pepQuantifList
	foreach my $stateInfo (@{$labelingInfo{STATES}}) {
		my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$stateInfo);
		if ($ratioType eq 'None') {next if $expCondID != $selStateID;}
		else {next if $expCondID != $testStateID && $expCondID != $refStateID;}
		foreach my $bioRep (split(/\./,$quantiObsIDs)) {
			foreach my $techRep (split(/&/,$bioRep)) { # no fraction for SWATH (so far...)
				foreach my $frac (split(/\+/,$techRep)) {
					my ($obsID,$pepQuantID,$anaID,$tgPos)=split(':',$frac);
					$anaList{$anaID}=1;
				}
			}
		}
	}
	my $anaString=join(',',keys %anaList);
	#my $pepQuantifString=join(',',keys %pepQuantifList);

	#my %runInfo;
	#my $sthRun=$dbh->prepare("SELECT A.ID_ANALYSIS,A.NAME,Q.NAME FROM ANALYSIS A,ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION IN ($pepQuantifString)");
	#$sthRun->execute;
	#while (my ($anaID,$anaName,$quantName)=$sthRun->fetchrow_array) {
	#	@{$runInfo{$anaID}}=($anaName,$quantName);
	#}
	#$sthRun->finish;

	$dbh->disconnect;
	$dbhLite->disconnect;

	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);

	my $protValueStrg;
	if ($ratioType eq 'None') {
		$protValueStrg=($dispMeasure eq 'MQ_INT')? 'Intensity' : ($dispMeasure eq 'MQ_IBAQ')? 'iBAQ'  : ($dispMeasure eq 'MQ_LFQ')? 'LFQ' : 'MS/MS count';
		$protValueStrg.='='.sprintf "%.0f",$protValue;
	}
	else {
		$protValueStrg='Ratio=';
		if ($protValue < 1) {$protValueStrg.='1/'.sprintf "%.2f",1/$protValue;}
		else {$protValueStrg.=sprintf "%.2f",$protValue;}
	}
	$modStrg=($modStrg)? '-'.$modStrg : '';

	print qq
|<TABLE><TR><TD nowrap><FONT class="title2">Peptide raw data for <A href="javascript:sequenceView($protID,'$anaString')">$protAlias$modStrg</A> in $quantifTargetName ($protValueStrg):</FONT>
&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'"/>
</TD></TR></TABLE>
<BR>
<BR><FONT class="title3">Not yet implemented...</FONT><BR><BR>
|;
	exit;
}


####################<<<ajaxPeptideProtRuler>>>#####################
sub ajaxPeptideProtRuler { # Data display not yet implemented
	my $modProtID = param('id_prot');
	my ($protID, $modStrg) = ($modProtID =~ /^(\d+)-?(.*)/);
	my $dispMeasure = param('dispMeasure') || '';
	my $selTargetPos = param('ratio');
	my $trueTargetPos = abs($selTargetPos);

	###<Connecting to the database>###
	my $dbh = &promsConfig::dbConnect;
	my $projectID = &promsMod::getProjectID($dbh, $selQuantifID, 'QUANTIFICATION');
	my $dbhLite = &promsQuantif::dbConnectProteinQuantification($selQuantifID, $projectID);
	my %proteinQuantifFamilies = &promsQuantif::getProteinQuantifFamilies;
	my $dispMeasName;
	foreach my $refMeas (@{$proteinQuantifFamilies{'MEASURES'}{'PROT_RULER'}}) {
		my ($measCode, $measName, $isOptional) = @{$refMeas};
		next if ($measCode ne $dispMeasure);
		$dispMeasName = $measName;
	}

	my ($quantifAnnot, $designID, $quantifMethodID, $quantifType, $quantifModID, $multiModifStrg) = $dbh->selectrow_array("SELECT QUANTIF_ANNOT, ID_DESIGN, M.ID_QUANTIFICATION_METHOD, M.CODE, Q.ID_MODIFICATION,
		GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ',')
		FROM QUANTIFICATION Q
		LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION = MQ.ID_QUANTIFICATION
		INNER JOIN QUANTIFICATION_METHOD M ON Q.ID_QUANTIFICATION_METHOD = M.ID_QUANTIFICATION_METHOD
		WHERE Q.ID_QUANTIFICATION = $selQuantifID GROUP BY Q.ID_QUANTIFICATION"
	);

	my %labelingInfo;
	my ($labelStrg, @labelInfo) = split('::', $quantifAnnot);
	my ($labelType) = ($labelStrg)? $labelStrg =~ /LABEL=(.+)/ : ('FREE');
	$labelType = uc($labelType);
	$labelingInfo{'LABELTYPE'} = $labelType;
	foreach my $infoStrg (@labelInfo) {
		my ($setting, $valueStrg) = split('=', $infoStrg);
		$valueStrg =~ s/#//g; # remove all ID tags
		@{$labelingInfo{$setting}} = split(';', $valueStrg);
	}

	# Get protein data
	my ($protAlias) = $dbh->selectrow_array("SELECT ALIAS FROM PROTEIN WHERE ID_PROTEIN = $protID");
	my ($paramID) = $dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER, NAME, CODE FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD = $quantifMethodID AND CODE = '$dispMeasure'");
	my ($protValue) = $dbhLite->selectrow_array("SELECT QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIF_PARAMETER = $paramID AND ID_PROTEIN = $protID AND TARGET_POS = $trueTargetPos");

	# Process design info
	my ($selParentQuantif, $selParTargetPos) = split(':', $labelingInfo{PARENT_Q}[$selTargetPos-1]);
	my ($selParentQuantifID, $selParentStateID, $selParentPos) = split(',', $selParentQuantif);
	my ($quantifTargetName) = $dbh->selectrow_array("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION = $selParentStateID");

	my $sthAnaInfo = $dbh->prepare("SELECT A.ID_ANALYSIS FROM ANALYSIS A 
		INNER JOIN OBSERVATION O ON A.ID_ANALYSIS=O.ID_ANALYSIS
		INNER JOIN OBS_EXPCONDITION OE ON O.ID_OBSERVATION=OE.ID_OBSERVATION
		INNER JOIN EXPCONDITION E ON OE.ID_EXPCONDITION=E.ID_EXPCONDITION
		WHERE E.ID_EXPCONDITION=? ORDER BY A.DISPLAY_POS ASC"
	);
	my %anaList;
	$sthAnaInfo->execute($selParentStateID);
	while (my ($analysisID) = $sthAnaInfo->fetchrow_array) {
		unless ($anaList{$analysisID}) {
			$anaList{$analysisID} = 1;
		}
	}
	my $anaString = join(',', keys %anaList);

	$dbh->disconnect;
	$dbhLite->disconnect;

	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);
	my $protValueStrg .= $dispMeasName . ' = ' . sprintf("%.2E", $protValue);
	print qq
|<TABLE><TR><TD nowrap><FONT class="title2">Peptide raw data for <A href="javascript:sequenceView($protID, '$anaString')">$protAlias</A> in $quantifTargetName ($protValueStrg):</FONT>
&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'"/>
</TD></TR></TABLE>
<BR>
<BR><FONT class="title3">Not yet implemented...</FONT><BR><BR>
|;
	exit;
}


sub ajaxCheckQuantifStatus {
	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;
	my ($status)=$dbh->selectrow_array("SELECT STATUS FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	$dbh->disconnect;
	
	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);
	print $status;
	exit;
}

####################<<<ajaxPeptideSwath>>>#####################
sub ajaxPeptideSwath {
#print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
	my $modProtID=param('id_prot');
	my ($protID,$modStrg)=($modProtID=~/^(\d+)-?(.*)/); $modStrg='' unless $modStrg;
	my $selRatioPos=param('ratio');
	my $trueRatioPos=abs($selRatioPos);
	#my $selRatioIdx=$selRatioPos-1;

	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;
	my $projectID=&promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION');
	my $dbhLite=&promsQuantif::dbConnectProteinQuantification($selQuantifID,$projectID);

	#my ($quantifAnnot,$quantifModID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_MODIFICATION FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	my ($quantifAnnot,$quantifMethodID,$quantifModID,$multiModifStrg)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,Q.ID_QUANTIFICATION_METHOD,Q.ID_MODIFICATION,GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ',')
																	 FROM QUANTIFICATION Q
																	 LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
																	 WHERE Q.ID_QUANTIFICATION=$selQuantifID GROUP BY Q.ID_QUANTIFICATION");
	#my ($textModRes,$htmlModRes)=('','');
	#my ($isModifQuantif,$isMultiModifQuantif)=(0,0);
	#my %quantifModifInfo;
	#if ($quantifModID || $multiModifStrg) {
	#	$isModifQuantif=1;
	#	my $modifStrg=$quantifModID || $multiModifStrg;
	#	$quantifModID=$modifStrg if $modifStrg !~ /,/; # <- new single PTM quantif (also uses MULTIMODIF_Q table!)
	#	my $sthMod=$dbh->prepare("SELECT ID_MODIFICATION,DISPLAY_CODE,DISPLAY_COLOR FROM MODIFICATION WHERE ID_MODIFICATION IN ($modifStrg) ORDER BY FIELD(ID_MODIFICATION,$modifStrg)");
	#	my $modifRank=0;
	#	$sthMod->execute;
	#	while (my ($modID,$displayCode,$displayColor)=$sthMod->fetchrow_array) {
	#		$modifRank++;
	#		$quantifModifInfo{ID2RANK}{$modID}=$modifRank;
	#		$quantifModifInfo{$modifRank}=[$modID,$displayCode,$displayColor];
	#	}
	#	$sthMod->finish;
	#	$isMultiModifQuantif=1 if $modifRank > 1;
	#}
	my %quantifModifInfo;
	my ($isModifQuantif,$isMultiModifQuantif)=&promsQuantif::getQuantifModificationInfo($dbh,$quantifModID,$multiModifStrg,\%quantifModifInfo);
	my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
	my ($labelType)=($labelStrg)? $labelStrg=~/LABEL=(.+)/ : ('FREE');
	$labelType=uc($labelType); # iTRAQ to ITRAQ
	my (%labelingInfo,%condToState,%condToObs);
	$labelingInfo{'LABELTYPE'}=$labelType;
	foreach my $infoStrg (@labelInfo) {
		my ($setting,$valueStrg)=split('=',$infoStrg);
		$valueStrg=~s/#//g; # remove all ID tags
		@{$labelingInfo{$setting}}=split(';',$valueStrg);
	}
	#my $ratioType=$labelingInfo{RATIO_TYPE}[0]; # PROT_RATIO_PEP
	my $numRatios=scalar @{$labelingInfo{RATIOS}};
	my $numStates=scalar @{$labelingInfo{STATES}};

	if ($modStrg && !$isModifQuantif) { # fall back to whole protID if called with a modProtID for a non-modif quantif (eg. called from log2 plot)
		$modProtID=$protID;
		$modStrg='';
	}
	my ($protAlias,$protLength)=$dbh->selectrow_array("SELECT ALIAS,PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=$protID");
	$protLength=0 unless $protLength;

	my ($ratioParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=$quantifMethodID AND CODE='RATIO'");
	#my $trueProtRatio;
	my ($textModRes,$htmlModRes,$siteQueryStrg)=('','','');
	if ($modStrg) { # modification quantif
		my $fileModStrg=&promsQuantif::siteCode2QuantifFileCode($modStrg,$isMultiModifQuantif,$quantifModifInfo{ID2RANK});
		$siteQueryStrg="AND SITE_CODE='$fileModStrg'";
		$textModRes=&promsQuantif::displayModificationSites($modStrg,$quantifModifInfo{DISPLAY},'text');
		$htmlModRes=&promsQuantif::displayModificationSites($modStrg,$quantifModifInfo{DISPLAY},'html');
		##($trueProtRatio)=$dbh->selectrow_array("SELECT QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION P
		##									LEFT JOIN PROTQUANTIF_MODRES PQMR ON P.ID_PROT_QUANTIF=PQMR.ID_PROT_QUANTIF
		##									LEFT JOIN MODIFIED_RESIDUE MR ON PQMR.ID_MODIF_RES=MR.ID_MODIF_RES
		##									INNER JOIN QUANTIFICATION_PARAMETER Q ON P.ID_QUANTIF_PARAMETER=Q.ID_QUANTIF_PARAMETER AND Q.CODE='RATIO'
		##									WHERE P.ID_QUANTIFICATION=$selQuantifID AND ID_PROTEIN=$protID AND TARGET_POS=$trueRatioPos GROUP BY P.ID_PROT_QUANTIF
		##									HAVING GROUP_CONCAT(RESIDUE,POSITION ORDER BY POSITION SEPARATOR '.')='$dbModStrg'");
	}
	# else { # whole protein quantif
	# 	#($trueProtRatio)=$dbh->selectrow_array("SELECT QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION P,QUANTIFICATION_PARAMETER Q WHERE P.ID_QUANTIF_PARAMETER=Q.ID_QUANTIF_PARAMETER AND Q.CODE='RATIO' AND ID_QUANTIFICATION=$selQuantifID AND ID_PROTEIN=$protID AND TARGET_POS=$trueRatioPos");
	# }
	my ($trueProtRatio)=$dbhLite->selectrow_array("SELECT QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIF_PARAMETER=$ratioParamID AND ID_PROTEIN=$protID AND TARGET_POS=$trueRatioPos $siteQueryStrg");
	my $protRatio=($selRatioPos < 0)? 1/$trueProtRatio : $trueProtRatio;

	#########################
	###<Processing design>###
	#########################
	my (%anaList,%pepQuantifList); #,%usedReplicates
	my ($testStateID,$refStateID)=($selRatioPos > 0)? split(/\//,$labelingInfo{RATIOS}[$trueRatioPos-1]) : reverse (split(/\//,$labelingInfo{RATIOS}[$trueRatioPos-1]));
	my $sthExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
	$sthExpCondName->execute($testStateID);
	my ($testStateName)=$sthExpCondName->fetchrow_array;
	$sthExpCondName->execute($refStateID);
	my ($refStateName)=$sthExpCondName->fetchrow_array;
	$sthExpCondName->finish;

	my ($testStatePos,$refStatePos);
	my $statePos=0;
	my $runCount=0;
	foreach my $stateInfo (@{$labelingInfo{STATES}}) {
		$statePos++;
		my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$stateInfo);
		if ($expCondID == $testStateID) {$testStatePos=$statePos;}
		elsif ($expCondID == $refStateID) {$refStatePos=$statePos;}
		foreach my $bioRep (split(/\./,$quantiObsIDs)) {
			foreach my $techRep (split(/&/,$bioRep)) { # no fraction for SWATH (so far...)
				$runCount++;
				next if ($expCondID != $testStateID && $expCondID != $refStateID);
				my ($obsID,$pepQuantID,$anaID,$tgPos)=split(':',$techRep);
				$anaList{$anaID}=$runCount;
				$pepQuantifList{$pepQuantID}=1;
			}
		}
	}
	my $anaString=join(',',keys %anaList);
	my $pepQuantifString=join(',',keys %pepQuantifList);

	my %runInfo;
	my $sthRun=$dbh->prepare("SELECT A.ID_ANALYSIS,A.NAME,Q.NAME FROM ANALYSIS A,ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION IN ($pepQuantifString)");
	$sthRun->execute;
	while (my ($anaID,$anaName,$quantName)=$sthRun->fetchrow_array) {
		@{$runInfo{$anaID}}=($anaName,$quantName);
	}
	$sthRun->finish;
#print "T=$testStatePos, R=$refStatePos<BR>\n";
	##>Processing ratio design
	my (%columnStrg,%design);
	#my @trueStatesPos=($selRatioPos > 0)? ($testStatePos,$refStatePos) : ($refStatePos,$testStatePos);
	foreach my $statePos ($testStatePos,$refStatePos) {
		my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$labelingInfo{STATES}[$statePos-1]);
		my $bioRepCount=0;
		foreach my $bioRep (split(/\./,$quantiObsIDs)) {
			$bioRepCount++;
			foreach my $techRep (split(/&/,$bioRep)) { # no fraction for SWATH (so far...)
				my ($obsID,$pepQuantID,$anaID,$tgPos)=split(':',$techRep);
				push @{$columnStrg{$statePos}},"['$runInfo{$anaID}[0]',$anaID,'$runInfo{$anaID}[0] > $runInfo{$anaID}[1]']";
				push @{$design{$statePos}},["State$statePos:BioRep$bioRepCount:$anaList{$anaID}",$anaID];
#print "State$statePos:BioRep$bioRepCount:$anaList{$anaID}<BR>\n";
			}
		}
	}


	####################
	###<Peptide data>###
	####################

	###<Loading peptide quantif data from table.txt & excluded.txt>###
	my $dataDir="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/data";
	my $resultDir="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/results";
	my %peptideData;
	foreach my $file ('table.txt','excluded.txt') {
		next unless -e "$dataDir/$file";
		my $excluded=($file eq 'excluded.txt')? 1 : 0;
		open(PEP,"$dataDir/$file") || die "ERROR: $!";
		my $protFound=0;
		while(<PEP>) {
			next if $.==1;
			if (!/^$modProtID\t/) {
				next unless $protFound;
				last; # file is ordered by prot
			}
			$protFound=1;
			chomp;
			my ($prID,$pepSeqVarMod,$charge,$fragIon,$fragIonCharge,$labelType,$state,$bioRep,$run,$intensity)=split(/\t/,$_);
			next if $state !~ /State($testStatePos|$refStatePos)/; # restrict to selected states
			next if $intensity eq 'NA';
			$bioRep=~s/$state.//; # '.' -> any single-character separator
			@{$peptideData{$pepSeqVarMod.'_'.$charge}{$fragIon.'_'.$fragIonCharge}{"$state:$bioRep:$run"}}=($intensity,$excluded);
#print "$_<BR>\n" if ($pepSeqVarMod =~ /^MIKPFFHSLSEK.+/ && $fragIon eq 'b5');
		}
		close PEP;
	}

	if (-e "$resultDir/CVoutPep.txt") {
		open(PEP,"$resultDir/CVoutPep.txt") || die "ERROR: $!";
		my $protFound=0;
		while(<PEP>) { # 75817	DSDQTTIAEDLMK_2_y5_1	State1	0.145431865367481
			next if $.==1;
			if (!/^$modProtID\t/) {
				next unless $protFound;
				last; # file is ordered by prot
			}
			$protFound=1;
			chomp;
			my ($prID,$pepCode,$state,$intensity)=split(/\t/,$_);
			my ($pepSeqVarMod,$charge,$fragIon,$fragIonCharge)=split('_',$pepCode);
			$pepSeqVarMod=~s/\*/:/g; # : converted to * because incompatible with MSstats
			#foreach my $bioRep (keys %{$peptideData{$pepSeqVarMod.'_'.$charge}{$state}}) { # bioReps & runs are not specified => assume all are excluded
			#	foreach my $run (keys %{$peptideData{$pepSeqVarMod.'_'.$charge}{$state}{$bioRep}}) {
			#		$peptideData{$pepSeqVarMod.'_'.$charge}{$fragIon.'_'.$fragIonCharge}{"$state:$bioRep:$run"}[1]=2 if $peptideData{$pepSeqVarMod.'_'.$charge}{$fragIon.'_'.$fragIonCharge}{"$state:$bioRep:$run"}; # set exclusion to level 2
			#	}
			#}
			foreach my $runKey (keys %{$peptideData{$pepSeqVarMod.'_'.$charge}{$fragIon.'_'.$fragIonCharge}}) {
				$peptideData{$pepSeqVarMod.'_'.$charge}{$fragIon.'_'.$fragIonCharge}{$runKey}[1]=2 if $runKey=~/^$state:/; # set exclusion to level 2
			}

		}
		close PEP;
	}

	###<Fetching peptides info from DB>###
	my (%peptideInfo,%peptideBeg,%varModName,%varModText);
	my $sthPI=$dbh->prepare("SELECT PEP_SEQ,GROUP_CONCAT(PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),CHARGE,ABS(PEP_BEG),P.ID_ANALYSIS,P.ID_PEPTIDE,QUERY_NUM,PEP_RANK,SCORE
							FROM PEPTIDE P
							LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
							INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE AND P.ID_ANALYSIS IN ($anaString)
							WHERE ID_PROTEIN=? AND PEP_RANK=1 GROUP BY P.ID_PEPTIDE");
	$sthPI->execute($protID);
	while (my ($pepSeq,$varModCode,$charge,$beg,$anaID,@info)=$sthPI->fetchrow_array) {
		my $vModStrg;
		if ($varModCode) {$vModStrg='&'.$varModCode;} else {$varModCode=$vModStrg='';}
		if ($peptideData{$pepSeq.$vModStrg.'_'.$charge}) { # match on seq+vMod+charge because ids are not recorded in table.txt)
			$peptideInfo{$pepSeq}{$varModCode}{$charge}{$anaID}=\@info; # pepID,qNum,rank,score
			$peptideBeg{$pepSeq}=$beg;
			$varModText{$varModCode}=&promsMod::decodeVarMod($dbh,$pepSeq,$varModCode,\%varModName) if $varModCode;
		}
	}
	$sthPI->finish;

	$dbh->disconnect;
	$dbhLite->disconnect;

	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);
	my $normalRatio=($protRatio==1000 || $protRatio==0.001)? 0 : 1; # inf or no inf
	my $protRatioStrg;
	if ($protRatio < 1) {$protRatioStrg=($normalRatio)? '1/'.sprintf "%.2f",1/$protRatio : '1/&infin;';} #'∞'
	else {$protRatioStrg=($normalRatio)? sprintf "%.2f",$protRatio : '&infin;';}
	if ($modStrg) {
		$htmlModRes='-'.$htmlModRes;
		$textModRes='-'.$textModRes;
	}
	my $numTestRep=scalar @{$design{$testStatePos}};
	my $numRefRep=scalar @{$design{$refStatePos}};
	print qq
|<TABLE><TR><TD nowrap>
<FONT class="title2">Peptide fragments raw data for <A href="javascript:sequenceView($protID,'$anaString')">$protAlias$htmlModRes</A> in $testStateName/$refStateName (ratio=$protRatioStrg):</FONT>
&nbsp;&nbsp;<FONT class="title3">View:</FONT><SELECT class="title3" onchange="var divList=['pepHeat:$selQuantifID:$selRatioPos','pepList:$selQuantifID:$selRatioPos']; for (var i=0; i<2; i++) {var vis=(divList[i]==this.value)? '' : 'none'; document.getElementById(divList[i]).style.display=vis;}">
<OPTION value="pepHeat:$selQuantifID:$selRatioPos">Graphical</OPTION><OPTION value="pepList:$selQuantifID:$selRatioPos">List</OPTION>
</SELECT>
&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'"/>
</TD></TR></TABLE>
<DIV id="pepHeat:$selQuantifID:$selRatioPos"></DIV>
<DIV id="pepList:$selQuantifID:$selRatioPos" style="display:none">
<TABLE cellspacing=0>
<TR bgcolor="$darkColor">
<TH class="rbBorder">#</TH><TH class="rbBorder">&nbsp;Pos.&nbsp;<TH class="rbBorder">&nbsp;Peptides&nbsp;</TH><TH class="rbBorder">&nbsp;Charge&nbsp;</TH><TH class="rbBorder">&nbsp;Fragment&nbsp;</TH><TH class="rbBorder" colspan="$numTestRep">&nbsp;$testStateName&nbsp;</TH><TH class="bBorder" colspan="$numRefRep">&nbsp;$refStateName&nbsp;</TH></TR>
|;
	my $bgColor=$lightColor;
	my $pepCount=0;
	foreach my $pepSeq (sort{$peptideBeg{$a}<=>$peptideBeg{$b} || length($a)<=>length($b)} keys %peptideInfo) {
		foreach my $varMod (sort keys %{$peptideInfo{$pepSeq}}) {
			my $rootPepKey=($varMod)? $pepSeq.'&'.$varMod : $pepSeq;
			my $pepSeqVarMod=($varMod)? $pepSeq.'+'.$varModText{$varMod} : $pepSeq;
			foreach my $charge (sort{$a<=>$b} keys %{$peptideInfo{$pepSeq}{$varMod}}) {
				$pepCount++;
				my $pepKey=$rootPepKey.'_'.$charge;
				my $numFrag=scalar keys %{$peptideData{$pepKey}};
				print "<TR bgcolor=\"$bgColor\"><TD class=\"rBorder\" align=\"right\" valign=\"top\" rowspan=$numFrag>&nbsp;$pepCount&nbsp;</TD><TD align=\"right\" valign=\"top\" rowspan=$numFrag>&nbsp;$peptideBeg{$pepSeq}&nbsp;</TD><TH align=\"right\" valign=\"top\" nowrap rowspan=$numFrag>&nbsp;$pepSeqVarMod</TH><TH valign=\"top\"rowspan=$numFrag>&nbsp;$charge</SUP>+</SUP></TH>";
				my $firstFrag=1;
				foreach my $fragCharge (sort{&promsMod::sortSmart($a,$b)} keys %{$peptideData{$pepKey}}) {
					if ($firstFrag) {$firstFrag=0;}
					else {print "<TR bgcolor=\"$bgColor\">"}
					my ($frag,$ch)=split('_',$fragCharge);
					my $chStrg=($ch==1)? '+' : $ch.'+';
					print "<TH>&nbsp;$frag<SUP>$chStrg</SUP>&nbsp;</TH>";
					#>Looping through columns (runs)
					foreach my $statePos ($testStatePos,$refStatePos) {
						foreach my $refRun (@{$design{$statePos}}) {
							my ($runKey,$anaID)=@{$refRun};
							if ($peptideData{$pepKey}{$fragCharge}{$runKey}) {
								my $fragValue=($peptideData{$pepKey}{$fragCharge}{$runKey}[0]=~/None/)? 'None' : sprintf '%.0f',$peptideData{$pepKey}{$fragCharge}{$runKey}[0];
								print "<TD align=\"right\"><A href=\"javascript:drawSpectrum(null,'pep_",(join('_',@{$peptideInfo{$pepSeq}{$varMod}{$charge}{$anaID}}[0..2])),"')\">$fragValue&nbsp;</A>&nbsp;</TD>";
							}
							else {print "<TD align=\"right\">-&nbsp;</TD>";}
						}
					}
					print "</TR>\n";
				}
				$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
			}
		}
	}

	print qq
|</TABLE>
<BR>
</DIV>
<BR>
#==========#
var HMF=new heatMap({
	div:'pepHeat:$selQuantifID:$selRatioPos',
	// cellHeight:21,
	// editable:true,
	editableItems:{scope:true,type:true,color:true},
	entities:{row:'Fragments',column:'Run'},
	flagText:'excluded fragments',
	cellOnMouseOver: function(value,hCell) {
		var valueStrg;
		if (value==null) {valueStrg='* no value *';}
		else {valueStrg='Intensity: '+(Math.round(10*value)/10);}
		if (hCell) {
			var rowIdData=rowList[hCell.rowIdx].id.split(':'), colId=columnList[hCell.columnIdx].id;
			if (peptideId[rowIdData[0]][colId]) {valueStrg+='\\nPeptide score: '+(Math.round(100*peptideId[rowIdData[0]][colId][1])/100);}
			else {valueStrg+='\\n* no peptide identified *';}
		}
		return valueStrg;
	},
	singleCellSelection:true,
	cellOnClick:function(rowId,colId,status) {
		if (!status) return;
		var rowIdData=rowId.split(':');
		if (peptideId[rowIdData[0]][colId]) {drawSpectrum(null,peptideId[rowIdData[0]][colId][0]);}
		else {alert('Peptide was not identified in this run');}
	},
	normalization:{scope:'row',reference:'user',limitValues:{ref:0,min:0}},
	exportAsImage:['Export as image','Peptide_fragments_heat_map','$promsPath{cgi}/exportSVG.cgi']
});
|;

	##>Setting columns: runs and States
	my %stateGroups;
	my $anaCount=0; # index
	print 'HMF.setColumns([';
	foreach my $statePos ($testStatePos,$refStatePos) { # Test then Ref
		$stateGroups{$statePos}[0]=$anaCount; # index of 1st run in state
		print ',' if $statePos==$refStatePos;
		print join(',',@{$columnStrg{$statePos}});
		$anaCount+=scalar @{$columnStrg{$statePos}};
		$stateGroups{$statePos}[1]=$anaCount-1; # index of last run in state
	}
	print "])\n";

	##>Adding rows: peptide fragments
	my %peptideGroups;
	my $jsPeptideIdStrg="var peptideId={\n";
	$pepCount=0; # already used for List view
	my $fragCount=0;
	foreach my $pepSeq (sort{$peptideBeg{$a}<=>$peptideBeg{$b} || length($a)<=>length($b)} keys %peptideInfo) {
		foreach my $varMod (sort keys %{$peptideInfo{$pepSeq}}) {
			my $rootPepKey=($varMod)? $pepSeq.'&'.$varMod : $pepSeq;
			foreach my $charge (sort{$a<=>$b} keys %{$peptideInfo{$pepSeq}{$varMod}}) {
				$pepCount++;
				my $pepKey=$rootPepKey.'_'.$charge;
				@{$peptideGroups{$pepCount}{'PEP'}}=($pepSeq,$varMod,$charge);
				$peptideGroups{$pepCount}{'FRAG'}[0]=$fragCount; # index    unless $peptideGroups{$pepCount}{'FRAG'};
				my $firstFrag=1;
				foreach my $fragCharge (sort{&promsMod::sortSmart($a,$b)} keys %{$peptideData{$pepKey}}) {
					$peptideGroups{$pepCount}{'FRAG'}[1]=$fragCount;
					my ($frag,$ch)=split('_',$fragCharge);
					print "HMF.addRow(['$frag','$pepCount:$fragCount','Charge: $ch+'],[";
					if ($firstFrag) { # same peptideID for all fragments
						$jsPeptideIdStrg.=",\n" if $pepCount > 1;
						$jsPeptideIdStrg.="\t'$pepCount':{";
					}
					#>Looping through columns (runs)
					my @excluded;
					my $localFragCount=0; # index
					foreach my $statePos ($testStatePos,$refStatePos) {
						foreach my $refRun (@{$design{$statePos}}) {
							my ($runKey,$anaID)=@{$refRun};
							if ($firstFrag) { # generate colID 2 pep
								$jsPeptideIdStrg.=',' if $localFragCount;
								if ($peptideInfo{$pepSeq}{$varMod}{$charge}{$anaID}) {$jsPeptideIdStrg.="'$anaID':['pep_".(join('_',@{$peptideInfo{$pepSeq}{$varMod}{$charge}{$anaID}}[0..2]))."',$peptideInfo{$pepSeq}{$varMod}{$charge}{$anaID}[3]]";}
								else {$jsPeptideIdStrg.="'$anaID':null";}
							}
							print ',' if $localFragCount;

							if ($peptideData{$pepKey}{$fragCharge}{$runKey}) {
								if($peptideData{$pepKey}{$fragCharge}{$runKey}[0]=~/None/){print '';}
								else{print $peptideData{$pepKey}{$fragCharge}{$runKey}[0]}
								push @excluded,$localFragCount if $peptideData{$pepKey}{$fragCharge}{$runKey}[1];
							}
							$localFragCount++;
						}
					}
					print ']';
					if (scalar @excluded) {
						print ',[],[',join(',',@excluded),']';
					}
					print ");\n";
					$fragCount++;

					if ($firstFrag) {
						$jsPeptideIdStrg.='}';
						$firstFrag=0;
					}
				}
			}
		}
	}
	$jsPeptideIdStrg.="\n};\n";

	##>Defining groups
	print "HMF.defineGroups('column',[['$refStateName',1,'Reference state',$stateGroups{$refStatePos}[0],$stateGroups{$refStatePos}[1]],['$testStateName',2,'Test state',$stateGroups{$testStatePos}[0],$stateGroups{$testStatePos}[1]]]);\n";
	print "HMF.defineGroups('row',[";
	foreach my $pepCount (sort{$a<=>$b} keys %peptideGroups) {
		print ',' if $pepCount > 1;
		my ($pepSeq,$varMod,$charge)=@{$peptideGroups{$pepCount}{'PEP'}};
		my $pepSeqVarMod=($varMod)? $pepSeq.'+'.$varModText{$varMod} : $pepSeq;
		#my $rowID='pep_'.(join('_',@{$peptideInfo{$pepSeq}{$varMod}{$charge}}[0..2]));
		#print "['$pepSeqVarMod',$pepCount,'Charge: $charge+\\nScore: $peptideInfo{$pepSeq}{$varMod}{$charge}[3]\\nPosition: $peptideBeg{$pepSeq}-",($peptideBeg{$pepSeq}+length($pepSeq)-1),"',$peptideGroups{$pepCount}{FRAG}[0],$peptideGroups{$pepCount}{FRAG}[1]]";
		print "['$pepSeqVarMod',$pepCount,'Charge: $charge+\\nPosition: $peptideBeg{$pepSeq}-",($peptideBeg{$pepSeq}+length($pepSeq)-1),"',$peptideGroups{$pepCount}{FRAG}[0],$peptideGroups{$pepCount}{FRAG}[1]]";
	}
	print qq
|]);
HMF.draw();
var rowList=HMF.getRowList();
var columnList=HMF.getColumnList();
$jsPeptideIdStrg
|;

	exit;
}

####################<<<ajaxPeptideAbundance>>>#####################
sub ajaxPeptideAbundance {

	my $modProtID=param('id_prot');
	my ($protID,$modStrg)=($modProtID=~/^(\d+)-?(.*)/);
	my $dispMeasure=param('dispMeasure') || '';
	my $selTargetPos=param('ratio');
	#my @dispTargetPos=split(',',param('dispTargetPosAjax'));

	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;
	my $projectID=&promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION');
	my $dbhLite=&promsQuantif::dbConnectProteinQuantification($selQuantifID,$projectID);

	#my ($quantifAnnot,$designID,$quantifType,$quantifModID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_DESIGN,M.CODE,Q.ID_MODIFICATION FROM QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$selQuantifID");
	my ($quantifAnnot,$designID,$quantifMethodID,$quantifType,$quantifModID,$multiModifStrg)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_DESIGN,M.ID_QUANTIFICATION_METHOD,M.CODE,Q.ID_MODIFICATION,GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ',')
																	 FROM QUANTIFICATION Q
																	 LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
																	 INNER JOIN QUANTIFICATION_METHOD M ON Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD
																	 WHERE Q.ID_QUANTIFICATION=$selQuantifID GROUP BY Q.ID_QUANTIFICATION");
	#my $isModifQuantif=($quantifModID || $multiModifStrg)? 1 : 0;
	#my $isMultiModifQuantif=($multiModifStrg)? 1 : 0;
	my %quantifModifInfo;
	my ($isModifQuantif,$isMultiModifQuantif)=&promsQuantif::getQuantifModificationInfo($dbh,$quantifModID,$multiModifStrg,\%quantifModifInfo);
	my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
	my ($labelType)=($labelStrg)? $labelStrg=~/LABEL=(.+)/ : ('FREE');
	$labelType=uc($labelType); # iTRAQ to ITRAQ
	my (%labelingInfo,%condToState,%condToObs);
	$labelingInfo{'LABELTYPE'}=$labelType;
	#if ($designID || $labelType =~ /SILAC|ITRAQ/) {
		foreach my $infoStrg (@labelInfo) {
			my ($setting,$valueStrg)=split('=',$infoStrg);
			$valueStrg=~s/#//g; # remove all ID tags
			@{$labelingInfo{$setting}}=split(';',$valueStrg);
		}
	#}

	if ($modStrg && !$isModifQuantif) { # fall back to whole protID if called with a modProtID for a non-modif quantif (eg. called from log2 plot)
		$modProtID=$protID;
		$modStrg='';
	}
	my ($protAlias,$protLength)=$dbh->selectrow_array("SELECT ALIAS,PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=$protID");
	$protLength=0 unless $protLength;

	my ($paramID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER,NAME,CODE FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=$quantifMethodID AND CODE='$dispMeasure'");

	my ($displayModStrg,$siteQueryStrg)=('','');
	if ($modStrg) { # modification quantif
		my $fileModStrg=&promsQuantif::siteCode2QuantifFileCode($modStrg,$isMultiModifQuantif,$quantifModifInfo{ID2RANK});
		$displayModStrg='-'.&promsQuantif::displayModificationSites($modStrg,$quantifModifInfo{DISPLAY},'html');
		$modProtID=$protID.'-'.$fileModStrg; # redefine modProtID to match quantif file format !!!!!!!!!!!!!!!!!!!!!!!!
		$siteQueryStrg="AND SITE_CODE='$fileModStrg'";
	}
	my $sthPR=$dbhLite->prepare("SELECT QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIF_PARAMETER=$paramID AND ID_PROTEIN=$protID AND TARGET_POS=$selTargetPos $siteQueryStrg");
	$sthPR->execute;
	my ($protValue)=$sthPR->fetchrow_array;
	$sthPR->finish;
	$dbhLite->disconnect;

	#########################
	###<Processing design>###
	#########################
	my ($numBioRep,$quantiObsIDs,$selStateID)=split(',',$labelingInfo{STATES}[$selTargetPos-1]);
	my ($quantifTargetName)=$dbh->selectrow_array("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=$selStateID");
	my %anaList;
	foreach my $bioRep (split(/\./,$quantiObsIDs)) {
		foreach my $techRep (split(/&/,$bioRep)) {
			foreach my $frac (split(/\+/,$techRep)) {
				my ($obsID,$pepQuantID,$anaID,$tgPos)=split(':',$frac);
				$anaList{$anaID}=1;
			}
		}
	}
	my $anaString=join(',',keys %anaList);
	$dbh->disconnect;

	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);
	my $protValueStrg=($dispMeasure eq 'SUM_INT')? 'Sum' : ($dispMeasure eq 'MEAN_INT')? 'Mean' : ($dispMeasure eq 'GEO_MEAN_INT')? 'Geometric mean' : ($dispMeasure eq 'MEDIAN_INT')? 'Median' : 'LFQ';
	$protValueStrg.='='.sprintf "%.2E",$protValue;
	print qq
|<TABLE><TR><TD nowrap><FONT class="title2">Abundance for <A href="javascript:sequenceView($protID,'$anaString')">$protAlias$displayModStrg</A> in $quantifTargetName ($protValueStrg):</FONT>
&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'"/>
</TD></TR></TABLE>
<BR>
<BR><FONT class="title3">Not yet implemented...</FONT><BR><BR>
|;
	exit;
}

sub ajaxComputeQuantifiedProteins { # globals: $selQuantifID
	my $isModifQuantif=param('isModif');
	my $quantifMethodID=param('qMethID');
	my $paramCode=param('paramCode');
	
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);
	
	my $dbh=&promsConfig::dbConnect;
	my $projectID=&promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION');
	my $dbhLite=&promsQuantif::dbConnectProteinQuantification($selQuantifID,$projectID);

	my ($paramID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=$quantifMethodID AND CODE='$paramCode'");
	##my ($numProteinsQuantified)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID AND ID_QUANTIF_PARAMETER=$paramID");
	my ($numProteinsQuantified)=$dbhLite->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIF_PARAMETER=$paramID");
	#my ($numIsoFormsQuantified)=($isModifQuantif)? $dbh->selectrow_array("SELECT COUNT(ID_PROT_QUANTIF) FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID AND ID_QUANTIF_PARAMETER=$paramID") : 0;
	my ($numIsoFormsQuantified)=($isModifQuantif)? $dbhLite->selectrow_array("SELECT COUNT(*) FROM (SELECT DISTINCT ID_PROTEIN,SITE_CODE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIF_PARAMETER=$paramID)") : 0;
	
	$dbh->disconnect;
	$dbhLite->disconnect;
	print "$numProteinsQuantified:$numIsoFormsQuantified";
	exit;
}

sub ajaxProteinStatistics { # ONLY for ratios. If multiple ratios were computed, ratio reversion (if any) applies only to selected ratio
#print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
	my $modProtID=param('id_prot');
	my $selRatioPos=param('ratio');
	my @dispTgPosAjax=(param('dispTargetPosAjax'))? split(',',param('dispTargetPosAjax')) : ($selRatioPos);
	foreach my $ratioPos (@dispTgPosAjax) {$dispRatios{$ratioPos}=[];} # %dispRatios is global but not defined yet
	my ($protID,$modResStrg)=($modProtID=~/^(\d+)-?(.*)/);

	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;
	my $projectID=&promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION');
	my $dbhLite=&promsQuantif::dbConnectProteinQuantification($selQuantifID,$projectID);

	my ($protAlias)=$dbh->selectrow_array("SELECT ALIAS FROM PROTEIN WHERE ID_PROTEIN=$protID");
	#my ($quantifAnnot,$designID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_DESIGN FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	my ($quantifAnnot,$designID,$quantifMethodID,$quantifModID,$multiModifStrg)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_DESIGN,ID_QUANTIFICATION_METHOD,Q.ID_MODIFICATION,GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ',')
																	 FROM QUANTIFICATION Q
																	 LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
																	 WHERE Q.ID_QUANTIFICATION=$selQuantifID GROUP BY Q.ID_QUANTIFICATION");
	#my ($isModifQuantif,$isMultiModifQuantif)=(0,0);
	##my (%quantifModifInfo,%quantifModifRanks);
	#my %quantifModifInfo;
	#if ($quantifModID || $multiModifStrg) {
	#	$isModifQuantif=1;
	#	my $modifStrg=$quantifModID || $multiModifStrg;
	#	$quantifModID=$modifStrg if $modifStrg !~ /,/; # <- new single PTM quantif (also uses MULTIMODIF_Q table!)
	#	my $sthMod=$dbh->prepare("SELECT ID_MODIFICATION,DISPLAY_CODE,DISPLAY_COLOR FROM MODIFICATION WHERE ID_MODIFICATION IN ($modifStrg) ORDER BY FIELD(ID_MODIFICATION,$modifStrg)");
	#	my $modifRank=0;
	#	$sthMod->execute;
	#	while (my ($modID,$displayCode,$displayColor)=$sthMod->fetchrow_array) {
	#		$modifRank++;
	#		$quantifModifInfo{ID2RANK}{$modID}=$modifRank;
	#		$quantifModifInfo{DISPLAY}{$modID}=[$displayCode,$displayColor];
	#	}
	#	$sthMod->finish;
	#	$isMultiModifQuantif=1 if $modifRank > 1;
	#}
	my %quantifModifInfo;
	my ($isModifQuantif,$isMultiModifQuantif)=&promsQuantif::getQuantifModificationInfo($dbh,$quantifModID,$multiModifStrg,\%quantifModifInfo);
	my (%stateName,%condToState,@ratioInfo,$alterHyp,$adjStrg,$trueAlpha);
	my ($labelType,$numStates,%labelingInfo,%stateInfo,$isPepIntensity);
	my $ratioType='Ratio'; # default
	my ($quantifSoftware,$algoVersion)=('myProMS',0); # default
	my $residualVar='biological'; # default
	my $hasMsmsCount=0; # dafault
	foreach my $infoStrg (split('::',$quantifAnnot)) {
		my ($setting,$valueStrg)=split('=',$infoStrg);
		@{$labelingInfo{$setting}}=split(/;/,$valueStrg);
		if ($setting eq 'LABEL') {$labelType=$valueStrg; $labelType=uc($labelType);}
		elsif ($setting eq 'SOFTWARE') {
			($quantifSoftware,$algoVersion)=split(';',$valueStrg);
			$quantifSoftware='MaxQuant' if $quantifSoftware eq 'MQ';
		}
		elsif ($setting eq 'RESIDUAL_VAR') {$residualVar=$valueStrg;} # Must be parsed before 'STATES'!!!
		elsif ($setting eq 'STATES') {
			my $statePos=0;
			if ($designID) {
				$valueStrg=~s/\#//g; # clean ID tags
				my $sthExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
				foreach my $stateData (split(';',$valueStrg)) {
					$statePos++;
					my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$stateData);
					$sthExpCondName->execute($expCondID);
					$condToState{$expCondID}=$statePos;
					($stateName{$statePos})=$sthExpCondName->fetchrow_array;
					#$stateInfo{$statePos}{'NUM_REPLIC'}=($numBioRep > 1)? $numBioRep : scalar (split(/&/,$quantiObsIDs)); # needed for targetPos selection for numPep filtering
					$stateInfo{$statePos}{'NUM_REPLIC'}=($residualVar eq 'biological')? $numBioRep : scalar (split(/&/,$quantiObsIDs)); # needed for targetPos selection for numPep filtering
				}
			}
			else {
				foreach my $stateData (split(';',$valueStrg)) {
					$statePos++;
					(my $numBioRep,$stateName{$statePos},my $repPosStrg)=split(',',$stateData);
					$stateName{$statePos}=~s/\./\+/g;
				}
			}
			$numStates=$statePos;
		}
		elsif ($setting eq 'RATIOS') {
			$valueStrg=~s/#//g if $designID; # for label-free experiments, RATIOS are displayed with '#' to tell
			@ratioInfo=split(';',$valueStrg);
		}
		elsif ($setting eq 'FDR_CONTROL') {
			$adjStrg=($valueStrg eq 'TRUE')? '_ADJ' : '';
		}
		elsif ($setting eq 'FDR_ALPHA') { # not always defined
			$trueAlpha=$valueStrg/100;
		}
		elsif ($setting eq 'ALTER') {
			$alterHyp=$valueStrg;
		}
		elsif ($setting eq 'RATIO_TYPE') {
			$ratioType=$valueStrg;
			($quantifSoftware,$algoVersion)=('myProMS',1) if $ratioType eq 'Ratio';
		}
		elsif ($setting eq 'NUM_TRUE_USED') {$hasMsmsCount=1;}
	}
	$algoVersion=2 if ($quantifSoftware eq 'myProMS' && !$algoVersion);
	$trueAlpha=0.05 unless $trueAlpha;
	$alterHyp='two.sided' unless $alterHyp;
	$adjStrg='_ADJ' if $ratioType=~/S\w+Ratio/; # overwrites FDR_CONTROL (just to be safe)

	my %replicatesPos;
	if ($quantifSoftware eq 'myProMS' && $algoVersion>=2) { # 'replicate'
		my $isPepIntensity=($labelingInfo{'MEAN_STATE'} && $labelingInfo{'MEAN_STATE'}[0]==1)? 1 : ($algoVersion==2 && $labelType eq 'FREE')? 1 : ($algoVersion>=3 && $labelingInfo{'ALGO_TYPE'}[0]=~/^(PEP_INTENSITY|TDA|DIA)$/)? 1 : 0; # 0 for v1 Ratio, v2 Label & v3 PEP_RATIO
		#foreach my $i (0..$#ratioInfo) {$dispRatios{$i+1}=[];}
		my %ratioCode=generateRatioCodes($labelType,$ratioType,$numStates,$isPepIntensity,$algoVersion);
		&getReplicateTargetPos($numStates,\%dispRatios,\%ratioCode,\%labelingInfo,\%stateInfo,\%replicatesPos); # populates $dispRatios{ratioPos}[0,1]
		&getMsmsStateTargetPos(\%dispRatios,\%ratioCode,\%labelingInfo) if $hasMsmsCount; # populates $dispRatios{ratioPos}[2,3]
	}
	
	my %paramInfo;
	my $sthQP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,NAME,CODE FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=$quantifMethodID");
	$sthQP->execute;
	while (my ($paramID,$parName,$parCode)=$sthQP->fetchrow_array) {
		@{$paramInfo{$paramID}}=($parCode,$parName);
	}
	$sthQP->finish;
	my (%statData,%paramName);
	#my $sthPS;
	my ($dispModSites,$siteQueryStrg)=('','');
	if ($modResStrg) { # modification quantif
		$dispModSites='-'.&promsQuantif::displayModificationSites($modResStrg,$quantifModifInfo{DISPLAY},'html');
		my $dbModStrg=&promsQuantif::siteCode2QuantifFileCode($modResStrg,$isMultiModifQuantif,$quantifModifInfo{ID2RANK});
		$siteQueryStrg="AND SITE_CODE='$dbModStrg'";
		##$sthPS=$dbh->prepare("SELECT GROUP_CONCAT(COALESCE(MODIF_RANK,''),':',RESIDUE,POSITION ORDER BY MODIF_RANK,POSITION SEPARATOR '.') AS MODIF_RES,CODE,NAME,TARGET_POS,QUANTIF_VALUE
		##						FROM PROTEIN_QUANTIFICATION P
		##						LEFT JOIN PROTQUANTIF_MODRES PQMR ON P.ID_PROT_QUANTIF=PQMR.ID_PROT_QUANTIF
		##						LEFT JOIN MODIFIED_RESIDUE MR ON PQMR.ID_MODIF_RES=MR.ID_MODIF_RES
		##						INNER JOIN QUANTIFICATION_PARAMETER Q ON P.ID_QUANTIF_PARAMETER=Q.ID_QUANTIF_PARAMETER
		##						WHERE P.ID_QUANTIFICATION=$selQuantifID AND ID_PROTEIN=$protID GROUP BY P.ID_PROT_QUANTIF
		##						HAVING MODIF_RES='$dbModStrg'");
		# $sthPS=$dbhLite->prepare("SELECT GROUP_CONCAT(COALESCE(MODIF_RANK,''),':',RESIDUE,POSITION ORDER BY MODIF_RANK,POSITION SEPARATOR '.') AS MODIF_RES,ID_QUANTIF_PARAMETER,TARGET_POS,QUANTIF_VALUE
		# 						FROM PROTEIN_QUANTIFICATION P
		# 						LEFT JOIN PROTQUANTIF_MODRES PQMR ON P.ID_PROT_QUANTIF=PQMR.ID_PROT_QUANTIF
		# 						LEFT JOIN MODIFIED_RESIDUE MR ON PQMR.ID_MODIF_RES=MR.ID_MODIF_RES
		# 						WHERE ID_PROTEIN=$protID GROUP BY P.ID_PROT_QUANTIF
		# 						HAVING MODIF_RES='$dbModStrg'");
	}
	# else {
	# 	##$sthPS=$dbh->prepare("SELECT NULL,CODE,NAME,TARGET_POS,QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION PQ,QUANTIFICATION_PARAMETER QP WHERE PQ.ID_QUANTIF_PARAMETER=QP.ID_QUANTIF_PARAMETER AND ID_QUANTIFICATION=$selQuantifID AND ID_PROTEIN=$protID"); # AND TARGET_POS=$ratioPos
	# 	$sthPS=$dbhLite->prepare("SELECT NULL,ID_QUANTIF_PARAMETER,TARGET_POS,QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION WHERE ID_PROTEIN=$protID"); # AND TARGET_POS=$ratioPos
	# }
	my $sthPS=$dbhLite->prepare("SELECT ID_QUANTIF_PARAMETER,TARGET_POS,QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION WHERE ID_PROTEIN=$protID $siteQueryStrg");
	$sthPS->execute;
	while (my ($paramID,$targetPos,$quantifValue)=$sthPS->fetchrow_array) {
		my ($parCode,$parName)=@{$paramInfo{$paramID}};
		$targetPos=0 unless $targetPos; # multi-ratio data
		$statData{$parCode}{$targetPos}=$quantifValue;
		if ($parCode eq 'RATIO' && $dispRatios{-$targetPos}) { # reverse ratio
			$statData{$parCode}{$targetPos}=1/$statData{$parCode}{$targetPos};
		}
		if ($parCode eq 'SD_GEO') {$paramName{$parCode}='Coeff. variation';} # overwrites SD_GEO definition
		else {$paramName{$parCode}=$parName;}
	}
	$sthPS->finish;

	$dbh->disconnect;
	$dbhLite->disconnect;

	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);
	my $numRatios=scalar @ratioInfo;
	my ($pepLabel,$distPepUsedStrg);
	if ($ratioType eq 'Ratio') {
		$pepLabel='Peptides used';
		$distPepUsedStrg=($numPepCode eq 'DIST_PEP_USED')? "$statData{DIST_PEP_USED}{0}/$statData{NUM_PEP_USED}{0}/$statData{NUM_PEP_TOTAL}{0}" : "$statData{NUM_PEP_USED}{0}/$statData{NUM_PEP_TOTAL}{0}";
	}
	else { # Super/Simple Ratio
		$pepLabel='All peptides';
		$distPepUsedStrg=$statData{NUM_PEP_TOTAL}{0};
	}
	print qq
|<FONT class="title2">All statistics for $protAlias$dispModSites&nbsp;</FONT>&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'"/>
<TABLE bgcolor=$darkColor>
<TR><TH nowrap align=right>$pepLabel :</TH><TH bgcolor=$lightColor align=left colspan=$numRatios>$distPepUsedStrg</TH></TR>
<TR><TH></TH>
|;
	my $ratioPos=0;
	my %linkedRatios;
	foreach my $ratioData (@ratioInfo) {
		$ratioPos++;
		my ($testStatePos,$refStatePos)=split(/\//,$ratioData);
		my $ratioTag='';
		if ($designID) {
			if ($testStatePos=~/%/) {
				$ratioTag=$encodedDegStrg;
				my ($linkedTestPos)=$testStatePos=~/%(\d+)/;
				my ($linkedRefPos)=$refStatePos=~/%(\d+)/;
				@{$linkedRatios{$ratioPos}}=($linkedTestPos,$linkedRefPos); # needed for NUM_PEP_USED
				$testStatePos=~s/%\d+//;
				$refStatePos=~s/%\d+//;
			}
			$testStatePos=$condToState{$testStatePos};
			$refStatePos=$condToState{$refStatePos};
		}
		next unless ($dispRatios{$ratioPos} || $dispRatios{-$ratioPos});
		my ($selectTag1,$selectTag2)=($ratioPos==abs($selRatioPos))? ('<FONT color=#DD0000>','</FONT>') : ('','');
		if ($dispRatios{$ratioPos}) {print "<TH nowrap>&nbsp;$selectTag1$stateName{$testStatePos}$ratioTag/$stateName{$refStatePos}$ratioTag$selectTag2&nbsp;</TH>";} # normal ratio
		else {print "<TH nowrap>&nbsp;$selectTag1$stateName{$refStatePos}$ratioTag/$stateName{$testStatePos}$ratioTag$selectTag2&nbsp;</TH>";} # reverse ratio
	}
	print "</TR>\n";
	my @quantifParams=('RATIO',"PVAL$adjStrg","NORM_PVAL$adjStrg",'SD_GEO');
	push @quantifParams,('CONF_LOW','CONF_UP') if $algoVersion>=3;
	push @quantifParams,('DIST_PEP_USED','NUM_PEP_USED') if $ratioType=~/S\w+Ratio/;
	push @quantifParams,'NUM_TRUE_USED' if $hasMsmsCount;
	my %log2Ratio;
	foreach my $parCode (@quantifParams) {
		next if $parCode =~ /^CONF_(LOW|UP)/; # written with Coef Var
		next unless $paramName{$parCode}; # NORM_PVAL$adjStrg not always defined
		$paramName{$parCode} ='CV [Ratio 95% conf.]' if $parCode eq 'SD_GEO';
		$paramName{$parCode}.='<SUP onmouseover="popup(\'Peptides formally identified by <B>MSMS</B> in each state\')" onmouseout="popout()">?</SUP>' if $parCode eq 'NUM_TRUE_USED';
		print "<TR><TH nowrap align=right>&nbsp;$paramName{$parCode} :</TH>";
		foreach my $rPos0 (sort{abs($a)<=>abs($b)} keys %dispRatios) { # (1..$ratioPos)
			my $rPos=abs($rPos0);
			my $valueStrg='-';
			if (defined $statData{$parCode}{$rPos}) {
				if ($parCode eq 'RATIO') {
					my $absRatio=($statData{$parCode}{$rPos} >= 1)? $statData{$parCode}{$rPos} : 1/$statData{$parCode}{$rPos};
					$valueStrg=sprintf "%.2f",$absRatio;
					$valueStrg='1/'.$valueStrg if ($statData{$parCode}{$rPos} < 1 && $valueStrg > 1);
					$log2Ratio{$rPos}=log($statData{$parCode}{$rPos})/log(2);
					$valueStrg.=' (Log2: '.(sprintf "%.2f",$log2Ratio{$rPos}).')';
				}
				elsif ($parCode=~/PVAL/) {
					$valueStrg=($statData{$parCode}{$rPos})? sprintf "%.2e",$statData{$parCode}{$rPos} : '~0';
					if ($parCode eq "PVAL$adjStrg" && $statData{$parCode}{$rPos} <= $trueAlpha) {
						$valueStrg.="<IMG src=\"$promsPath{images}/good.gif\">";
					}
				}
				elsif ($parCode eq 'SD_GEO') {
					#$valueStrg=($ratioType eq 'Ratio')? sprintf "%.1f\%",100 * abs(log($statData{$parCode}{$rPos})/log($statData{RATIO}{$rPos})) : sprintf "%.1f",100*$statData{$parCode}{$rPos}/abs(log($statData{RATIO}{$rPos})/log(2));
					#$valueStrg.=($statData{'CONF_LOW'})? '% ['.(sprintf "%.2f",$statData{'CONF_LOW'}{$rPos}).' - '.(sprintf "%.2f",$statData{'CONF_UP'}{$rPos}).']' : '% (Conf. int.: &plusmn;'.(sprintf "%.2f",1.96*$statData{$parCode}{$rPos}).' Log2 unit)'; # +/- conf interval
					my $cvPc=($ratioType eq 'Ratio')? 100*abs(log($statData{$parCode}{$rPos})/log($statData{RATIO}{$rPos})) : 100*$statData{$parCode}{$rPos}/abs($log2Ratio{$rPos}); # v1 : v2,v3
					my ($ratioLow,$ratioUp)=($ratioType eq 'Ratio')? ($statData{'CONF_LOW'}{$rPos},$statData{'CONF_UP'}{$rPos}) : ($statData{'CONF_LOW'})? (2**$statData{'CONF_LOW'}{$rPos},2**$statData{'CONF_UP'}{$rPos}) : (2**($log2Ratio{$rPos} - 1.96*$statData{$parCode}{$rPos}),2**($log2Ratio{$rPos} + 1.96*$statData{$parCode}{$rPos})); # v1 : v3 : v2
					if ($ratioLow > $ratioUp) { # Bug in data output from R in S*Ratio ???
						my $tempVal=$ratioLow;
						$ratioLow=$ratioUp;
						$ratioUp=$tempVal;
					}
# BUG? NACB1_YEAST in UPS1 Skiping
					my $ratioLowStrg=($ratioLow >= 1)? sprintf "%.2f",$ratioLow : sprintf "1/%.2f",1/$ratioLow;
					my $ratioUpStrg=($ratioUp >= 1)? sprintf "%.2f",$ratioUp : sprintf "1/%.2f",1/$ratioUp;
					$valueStrg=(sprintf "%.1f",$cvPc)."% [$ratioLowStrg - $ratioUpStrg]"; #.' '.$statData{$parCode}{$rPos};
#$valueStrg.=' ('.(sprintf "%.2f",2**($log2Ratio{$rPos} - 1.96*$statData{$parCode}{$rPos})).'-'.(sprintf "%.2f",2**($log2Ratio{$rPos} + 1.96*$statData{$parCode}{$rPos})).')'; 
				}
				else {
					$valueStrg=$statData{$parCode}{$rPos};
					if ($parCode=~/(NUM|DIST)_PEP_USED/) {
						my @pepInRep;
						my $okReplic=0;
						my @stateIdx=(-$rPos==$rPos0)? (0,1) : (1,0); # St2/St1 +/-reverse
						foreach my $s (@stateIdx) {
							$pepInRep[$s]=[];
							foreach my $replicPos (@{$dispRatios{$rPos0}[$s]}) {
								if ($statData{$parCode}{$replicPos}) {
									push @{$pepInRep[$s]},$statData{$parCode}{$replicPos};
									$okReplic=1;
								}
								else {push @{$pepInRep[$s]},0;}
							}
						}
						$valueStrg.=' ['.join(', ',@{$pepInRep[$stateIdx[0]]}).' / '.join(', ',@{$pepInRep[$stateIdx[1]]}).']<SUP onmouseover="popup(\'Distribution in <B>replicates</B> of each state\')" onmouseout="popout()">?</SUP>' if $okReplic; # TEST/REF
					}
				}
			}
			elsif ($parCode eq 'NUM_TRUE_USED') {
				my @stateIdx=(-$rPos==$rPos0)? (2,3) : (3,2); # [0,1] used for list of replicate targetPos!!!
				my @numTruePep;
				foreach my $s (@stateIdx) {
					#push @numTruePep,($dispRatios{$rPos0} && $dispRatios{$rPos0}[$s])? $dispRatios{$rPos0}[$s][0] : 0;
					push @numTruePep,($dispRatios{$rPos0} && $dispRatios{$rPos0}[$s] && $statData{$parCode}{ $dispRatios{$rPos0}[$s][0] })? $statData{$parCode}{ $dispRatios{$rPos0}[$s][0] } : 0;
				}
				$valueStrg=join(' / ',@numTruePep);
			}
			elsif ($linkedRatios{$rPos} && $parCode=~/(NUM|DIST)_PEP_USED/) {
				$valueStrg=($statData{$parCode}{$linkedRatios{$rPos}[0]} <= $statData{$parCode}{$linkedRatios{$rPos}[1]})? $statData{$parCode}{$linkedRatios{$rPos}[0]} : $statData{$parCode}{$linkedRatios{$rPos}[1]};
			}
			print "<TH bgcolor=$lightColor align=left nowrap>&nbsp;$valueStrg&nbsp;</TH>";
		}
		print "</TR>\n";
	}
	#if ($ratioType eq 'Ratio') { # old algo
	#	print "<TR><TH nowrap align=right>&nbsp;Confidence interval :</TH>";
	#	foreach my $rPos (1..$ratioPos) {
	#		my ($lowConf,$upConf);
	#		if ($alterHyp ne 'less') {
	#			$lowConf=(!defined $statData{'CONF_LOW'}{$rPos})? '?' : ($statData{'CONF_LOW'}{$rPos} >= 1)? sprintf "%.2f",$statData{'CONF_LOW'}{$rPos} : sprintf "1/%.2f",1/$statData{'CONF_LOW'}{$rPos};
	#			$lowConf=~s/0+\Z//;
	#		}
	#		if ($alterHyp ne 'greater') {
	#			$upConf=(!defined $statData{'CONF_UP'}{$rPos})? '?' : ($statData{'CONF_UP'}{$rPos} >= 1)? sprintf "%.2f",$statData{'CONF_UP'}{$rPos} : sprintf "1/%.2f",1/$statData{'CONF_UP'}{$rPos};
	#			$upConf=~s/0+\Z//;
	#		}
	#		my $confIntStrg=($alterHyp eq 'less')? "0-$upConf" : ($alterHyp eq 'greater')? "&ge;$lowConf" : "$lowConf-$upConf";
	#		print "<TH bgcolor=$lightColor align=left nowrap>&nbsp;$confIntStrg&nbsp;</TH>";
	#	}
	#	print "</TR>\n";
	#}
	print "</TABLE>\n";

	exit;
}

sub ajaxSearchConvertIdentifier {
	my $searchText=lc((split(/\s+/,param('TEXT')))[0]);
	my $quantifList=param('quantifList'); # can be quantifID or string of qID1_tgPos1[:qID2_tgPos2:...]
	my $quantifIdStrg='';
	if ($quantifList) {
		foreach my $quantif (split(':',$quantifList)) {
			$quantifIdStrg.=',' if $quantifIdStrg;
			$quantifIdStrg.=(split('_',$quantif))[0];
		}
	}
	my $projectID=&promsMod::cleanNumericalParameters(param('projectID'));
	
	####<Starting HTML>###
	print header(-type=>'text/plain'); warningsToBrowser(1);
	
	my $dbh=&promsConfig::dbConnect;
	
	####<Finding identifier(s) & species used for proteins>#### (needed for search keyword conversion from graphical library)
	my %allIdentifierCodes;
	my $sthAll=$dbh->prepare("SELECT CODE,ID_IDENTIFIER FROM IDENTIFIER");
	$sthAll->execute;
	while (my ($code,$id)=$sthAll->fetchrow_array) {$allIdentifierCodes{$code}=$id;}
	$sthAll->finish;
	
	unless ($projectID) {
		my ($item,$itemID)=($quantifIdStrg)? ('quantification',(split(',',$quantifIdStrg))[0]) : ('analysis',$analysisID);
		$projectID=&promsMod::getProjectID($dbh,$itemID,$item);
	}
	my ($usedIdentierIdStrg)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM PROJECT WHERE ID_PROJECT=$projectID");
	unless ($usedIdentierIdStrg) {
		my $sthIdent=($quantifIdStrg)? $dbh->prepare("SELECT DISTINCT D.IDENTIFIER_TYPE FROM ANA_QUANTIFICATION AQ
														INNER JOIN ANALYSIS_DATABANK AD ON AQ.ID_ANALYSIS=AD.ID_ANALYSIS
														INNER JOIN DATABANK D ON AD.ID_DATABANK=D.ID_DATABANK
														WHERE ID_QUANTIFICATION IN ($quantifIdStrg)")
									: $dbh->prepare("SELECT DISTINCT D.IDENTIFIER_TYPE FROM ANALYSIS_DATABANK AD
														INNER JOIN DATABANK D ON AD.ID_DATABANK=D.ID_DATABANK
														WHERE ID_ANALYSIS=$analysisID");
		$sthIdent->execute;
		my %usedIdentifiers;
		while (my ($identType)=$sthIdent->fetchrow_array) {
			my $identCode;
			if ($identType=~/UNIPROT_A/) {$identCode='AC';}
			elsif ($identType eq 'UNIPROT_ID') {$identCode='ID';}
			elsif ($identType eq 'FLYBASE_ID') {$identCode='FlyBase';}
			elsif ($identType eq 'GI_ACCESSION') {$identCode='GI';}
			elsif ($identType eq 'IPI_ACCESSION') {$identCode='IP';}
			else {next;}
			$usedIdentifiers{ $allIdentifierCodes{$identCode} }=1;
		}
		$sthIdent->finish;
		$usedIdentierIdStrg=join(',',keys %usedIdentifiers);
	}
	#<Species used in project
	my ($usedSpecies,$gr)=$dbh->selectrow_array("SELECT GROUP_CONCAT(DISTINCT ID_SPECIES SEPARATOR ','),1 AS GR FROM MASTER_PROTEIN MP
													INNER JOIN PROTEIN P ON MP.ID_MASTER_PROTEIN=P.ID_MASTER_PROTEIN
													WHERE P.ID_PROJECT=$projectID GROUP BY GR");
	my ($speciesRestrict1,$speciesRestrict2)=($usedSpecies)? ("INNER JOIN MASTER_PROTEIN MP ON MI.ID_MASTER_PROTEIN=MP.ID_MASTER_PROTEIN AND MP.ID_SPECIES IN ($usedSpecies)",
															  "AND ID_SPECIES IN ($usedSpecies)") : ('','');

	####<Searching MASTERPROT_IDENTIFIER & MASTER_PROTEIN>####
	my @sthSearch=($dbh->prepare("SELECT DISTINCT MI.ID_MASTER_PROTEIN FROM MASTERPROT_IDENTIFIER MI $speciesRestrict1 WHERE MI.VALUE LIKE '$searchText%'"), # partial at end of word
				   $dbh->prepare("SELECT ID_MASTER_PROTEIN FROM MASTER_PROTEIN WHERE PROT_DES LIKE '%$searchText%' $speciesRestrict2")); # full partial
	my %matchedIDs;
	foreach my $sthS (@sthSearch) {
		$sthS->execute;
		while (my ($masterProtID)=$sthS->fetchrow_array) {$matchedIDs{$masterProtID}=1;}
		$sthS->finish;
	}

	my %matchedIdentifiers;
	$matchedIdentifiers{uc($searchText)}=1;
	if (scalar keys %matchedIDs) {
		my $masterIdStrg=join(',',keys %matchedIDs);
		
		#<Match mapped identifiers
		#my $sthMatch=$dbh->prepare("SELECT GROUP_CONCAT(VALUE SEPARATOR '|'),1 AS GR FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN IN ($masterIdStrg) AND ID_IDENTIFIER IN ($usedIdentierIdStrg) GROUP BY GR");			
		my $sthMatch1=$dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN IN ($masterIdStrg) AND ID_IDENTIFIER IN ($usedIdentierIdStrg)");			
		$sthMatch1->execute;
		while (my ($identifier)=$sthMatch1->fetchrow_array) {$matchedIdentifiers{$identifier}=1;}
		$sthMatch1->finish;
		
		#<Match proteins linked to matched master protein(s)
		my $sthMatch2=$dbh->prepare("SELECT CONCAT(IDENTIFIER,'|',ALIAS) FROM PROTEIN WHERE ID_MASTER_PROTEIN IN ($masterIdStrg) AND ID_PROJECT=$projectID");			
		$sthMatch2->execute;
		while (my ($identAlias)=$sthMatch2->fetchrow_array) {
			foreach my $identifier (split(/\|/,$identAlias)) {
				next if (length($identifier) < 3 || (length($identifier) < 4 && $identifier!~/\d/)); # try to avoid sp, trembl, ...
				$matchedIdentifiers{$identifier}=1;
			}
		}
		$sthMatch2->finish;
	}
	$dbh->disconnect;
	
	print join('|',keys %matchedIdentifiers);
	
	exit;
}

sub ajaxGetGoTerms {
	my ($goID,$aspect)=split(',',param('goStrg'));
	my $projectID=param('projectID');
	my $dbh=&promsConfig::dbConnect;
	my $selGOfromParentID=$dbh->prepare("SELECT ID_GOANALYSIS, PARAM_STRG from GO_ANALYSIS where ID_PARENT_GOANA=?");
	$selGOfromParentID->execute($goID);
	my %chilGoID;
	while (my($goAnaID, $paramStrg)=$selGOfromParentID->fetchrow_array) {
		my ($binStrg,$sizeStrg)=split(/;/,$paramStrg);
		my $binValue = (split(/=/,$binStrg))[1];
		$chilGoID{$goAnaID}=$binValue;
	}
	$selGOfromParentID->finish;
	my $onChange=param('onChange') || 'ajaxGoDecorateGraph';

	####<Starting HTML>###
	print header(-type=>'text/plain'); warningsToBrowser(1);
	print "<SELECT id=\"goTermsSELECT\" class=\"highlight\" onchange=\"$onChange(this.value,this.options[this.selectedIndex].text)\"><OPTION value=\"\">-= Select a GO term =-</OPTION>\n";

	if (!scalar(keys %chilGoID)) {
		my $resultFile="$promsPath{go_unix}/project_$projectID/$goID/results_$aspect.txt";
		open (TERMS,$resultFile);
		my $numTerms=0;
		while (<TERMS>) {
			if (/^(GO:\d+)\t([^\t]+)\t([^\t]+)/) {
				next if $3 > 1;
				my $pValue=sprintf "%.2e",$3;
				print "<OPTION value=\"$goID,$aspect,$1\">$2 [$1] ($pValue)</OPTION>\n";
				$numTerms++;
			}
		}
		close TERMS;
		unless ($numTerms) {print "<OPTION value='' disabled>No terms found</OPTION>\n";}
	}
	else {##there is a parent
		my @logRatioStrg=split(/;/,$dbh->selectrow_array("SELECT PARAM_STRG from GO_ANALYSIS where ID_GOANALYSIS=$goID"));
		my %ratioValue;
		foreach my $item (@logRatioStrg) {
			if ($item =~/logRatios/) {
				my ($value)=(split(/=/,$item))[1];
				my $count=0;
				my ($oldRatio, $strgRatioEnd);
				foreach my $ratio (split(/,/,$value)) {
					my $strgRatio="";
					if ($count == 0) {
						$strgRatio = "1/&infin;&nbsp;".&revertLog2($ratio);
					}
					elsif ($count==3) {
						$strgRatio = "$oldRatio&nbsp;".&revertLog2($ratio);
						$ratioValue{$count}=$strgRatio;
						$strgRatioEnd = &revertLog2($ratio)."&nbsp;&infin;";
						$ratioValue{$count+1}=$strgRatioEnd;
						last;
					}
					else {
						$strgRatio = "$oldRatio&nbsp;".&revertLog2($ratio);
					}
					$ratioValue{$count}=$strgRatio;
					$oldRatio=&revertLog2($ratio);
					$count++;
				}
			}

		}

		foreach my $goAnaID (sort{$chilGoID{$a} <=>$chilGoID{$b}} keys %chilGoID) {
			print "<OPTGROUP label=\"bin $ratioValue{$chilGoID{$goAnaID}} \">";
			my $resultFile="$promsPath{go_unix}/project_$projectID/$goAnaID/results_$aspect.txt";
			open (TERMS,$resultFile);
			my $numTerms=0;
			while (<TERMS>) {
				if (/^(GO:\d+)\t([^\t]+)\t([^\t]+)/) {
					next if $3 > 1;
					my $pValue=sprintf "%.2e",$3;
					print "<OPTION value=\"$goAnaID,$aspect,$1,$chilGoID{$goAnaID}:$ratioValue{$chilGoID{$goAnaID}}\">$2 [$1] ($pValue)</OPTION>\n";
					$numTerms++;
				}
			}
			close TERMS;
			print "</OPTGROUP>";
			unless ($numTerms) {print "<OPTION value='' disabled>No terms found</OPTION>\n";}
		}
	}
	print "</SELECT>\n";
	$dbh->disconnect;
}

sub ajaxGetGoTermProteins {
	my ($goID,$aspect,$termID)=split(',',param('goStrg'));
	my $projectID=param('projectID');

	####<Starting HTML>###
	print header(-type=>'text/plain',-charset => 'utf-8'); warningsToBrowser(1);
	my $resultFile="$promsPath{go_unix}/project_$projectID/$goID/results_$aspect.txt";
	my $protIdStrg;
	open (TERMS,$resultFile);
	while (<TERMS>) {
		if (/^$termID\t/) {
			$protIdStrg=(split(/\t/,$_))[5];
			last;
		}
	}
	close TERMS;
	print $protIdStrg;
}

sub ajaxListDecorateGraph {
    my $listID=param('listID');
	my $noSites=param('noSite') || 0;

    ###<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);

	###<Fetching list of proteins/sites>###
	my $dbh=&promsConfig::dbConnect;

    my $sthList=$dbh->prepare("SELECT NAME,LIST_TYPE FROM CATEGORY WHERE ID_CATEGORY=?");
	$sthList->execute($listID);
    my ($listName,$listType)=$sthList->fetchrow_array;
	$listName=&promsMod::shortenName($listName,30);
	$listName=~s/['"]/\./g;
	$sthList->finish;
	
	my %protList;
	my $numPTMs=&promsQuantif::fetchCustomList($dbh,$listID,\%protList,$noSites);
	my $hasPTMs=($numPTMs)? 1 : 0;
	print $listName.'::'.$hasPTMs.'::';
	if (scalar keys %protList) {print join(',',keys %protList);}
    else {print 0;}
	
	$dbh->disconnect;
    exit;
}

sub ajaxManageSaveProteins {
	my $job=param('job');
	my $listType=param('listType') || 'PROT';
	my $themeID=param('themeID') || 0;

	###<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);

	my $dbh=&promsConfig::dbConnect;

	###<Fetching list of available themes/classifications>###
	if ($job eq 'getThemes') {
		my $projectID=param('projID');
		my %themes;
		my $sthTh=$dbh->prepare("SELECT ID_CLASSIFICATION,NAME FROM CLASSIFICATION WHERE ID_PROJECT=$projectID");
		$sthTh->execute;
		while (my ($thID,$thName)=$sthTh->fetchrow_array) {
			$themes{$thID}=$thName;
		}
		$sthTh->finish;
		my $entities=($listType eq 'SITE')? 'sites' : 'proteins'; # used between SPAN because can be changed by JS ajaxManageSaveProteins() in promsMod
		print qq
|<TABLE border=0 bgcolor=$darkColor id="storeList"><TR>
<TH valign=top nowrap>&nbsp;<FONT style="font-size:15px;">Save checked <SPAN id="listTypeEntity">$entities</SPAN> in</FONT></TH>
<TD valign=top>
<SELECT id="themeSELECT" style="width:250px;font-weight:bold;" onchange="ajaxManageSaveProteins(this.value,saveListGlobals.listType)">
|; # saveListGlobals.listType: global JS variable because of possible swtich between PROT and SITE without new AJAX call
		my $displayStatus;
		if (scalar keys %themes) { # exists at least 1 classification
			print "<OPTION selected value=\"getLists:0\">-=Select a Theme=-</OPTION>\n";
			foreach my $thID (sort{lc($themes{$a}) cmp lc($themes{$b})} keys %themes) {
				print "<OPTION value=\"getLists:$thID\">$themes{$thID}</OPTION>\n";
			}
			$displayStatus='none';
		}
		else {$displayStatus='block';}
		print qq
|<OPTION value="getLists:-1">-=Create a new Theme=-</OPTION>
</SELECT><BR>
<INPUT type="text" id="themeNameIN" value="" placeholder="New Theme name" style="width:250px;display:$displayStatus"/>
</TD>
<TD valign=top><DIV id="listSelectionDIV">|;
		my $disabSave='disabled';
		unless (scalar keys %themes) { # no Theme
			$disabSave='';
			print qq
|<SELECT id="listSelectionSELECT" style="width:250px;font-weight:bold" onchange="ajaxManageSaveProteins(this.value,'$listType')">
<OPTION value="selectList:-1:-1">-=Create a new List=-</OPTION>
</SELECT><BR>
<INPUT type="text" id="listNameIN" value="" placeholder="New List name" style="width:250px;"/>
|;
		}
		print qq
|</DIV></TD>
<TD nowrap valign=top><INPUT type="button" id="saveProtBUTTON" value=" Save " onclick="ajaxManageSaveProteins('saveProt')" $disabSave>&nbsp;<INPUT type="button" value="Cancel" onclick="ajaxManageSaveProteins('cancel')"></TD>
<TR><TD></TD>
<TH align=left colspan=3><SELECT name="saveProtAct" id="saveProtAct"><OPTION value="add">Add to List</OPTION><OPTION value="remove">Remove from List</OPTION><OPTION value="replace">Replace List content</OPTION></SELECT>
</TH>
</TR></TABLE>
|;
	}
	####<Fetching list of available lists/categories>####
	elsif ($job eq 'getLists') {
		my %lists;
		if ($themeID > 0) {
			my $sthL=$dbh->prepare("SELECT L.ID_CATEGORY,NAME,DISPLAY_POS,COUNT(CP.ID_CATEGORY_PROTEIN) FROM CATEGORY L,CATEGORY_PROTEIN CP WHERE L.ID_CATEGORY=CP.ID_CATEGORY AND ID_CLASSIFICATION=$themeID AND LIST_TYPE='$listType' GROUP BY L.ID_CATEGORY");
			my @sthNoMod=(
				$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,COMPARISON CO WHERE CO.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID AND LIST_TYPE='$listType'"), # COMPARISON
				$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,CAT_COMPARISON CC WHERE CC.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID AND LIST_TYPE='$listType'"), # CAT_COMPARISON
				$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,GO_ANALYSIS GA WHERE GA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID AND LIST_TYPE='$listType'"), # GO_ANALYSIS
				$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,PATHWAY_ANALYSIS PA WHERE PA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID AND LIST_TYPE='$listType'"), # PATHWAY_ANALYSIS
				$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,EXPLORANALYSIS EA WHERE EA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID AND LIST_TYPE='$listType'") # EXPLORANALYSIS
			);
			$sthL->execute;
			while (my ($listID,$listName,$displayPos,$numItems)=$sthL->fetchrow_array) {
				@{$lists{$listID}}=($listName,$displayPos,$numItems,0);
			}
			foreach my $sth (@sthNoMod) {
				$sth->execute;
				while (my ($listID)=$sth->fetchrow_array) {
					$lists{$listID}[3]=1 if $lists{$listID}; # update modif flag
				}
			}
			$sthL->finish;
			foreach my $sth (@sthNoMod) {$sth->finish;}
		}
		my $displayStatus;
		print "<SELECT id=\"listSelectionSELECT\" style=\"width:250px;font-weight:bold\" onchange=\"ajaxManageSaveProteins(this.value,'$listType')\">\n";
		if (scalar keys %lists) {
			print "<OPTION value=\"selectList:$themeID:0\">-=Select a List=-</OPTION>\n";
			foreach my $listID (sort{$lists{$a}[1] cmp $lists{$b}[1]} keys %lists) {
				print "<OPTION value=\"selectList:$themeID:$listID\"";
				print ' disabled' if $lists{$listID}[3]; # list is used and cannot be modified
				print ">$lists{$listID}[0] [x$lists{$listID}[2]]</OPTION>\n";
			}
			$displayStatus='none';
		}
		else {$displayStatus='block';}
		print qq
|<OPTION value="selectList:$themeID:-1">-=Create a new List=-</OPTION>
</SELECT>
<INPUT type="text" id="listNameIN" value="" placeholder="New List name" style="width:250px;display:$displayStatus"/>
|;
	}
	####<Saving proteins in List>####
	else {
		my $listID=param('listID');
		#my $modifID=param('modifID');
		#my $replace=(param('replace'))? param('replace') : 0;
		my $saveProtAction=param('saveProtAct') || 'add' ; # add/remove/replace
		my ($themeName,$listName);
		if ($themeID==-1) { # Adding new theme
			my $projectID=param('projID');
			$themeName=param('themeName');
			($themeID)=$dbh->selectrow_array("SELECT MAX(ID_CLASSIFICATION) FROM CLASSIFICATION");
			$themeID++;
			my $sthInsCL=$dbh->prepare("INSERT INTO CLASSIFICATION (ID_CLASSIFICATION,NAME,ID_PROJECT) VALUES ($themeID,?,?)");
			$sthInsCL->execute($themeName,$projectID);
			$sthInsCL->finish;
			$listID=-1; # to be safe
		}
		else {
			$themeID=~s/\D+//g; # clean to prevent prepare
			($themeName)=$dbh->selectrow_array("SELECT NAME FROM CLASSIFICATION WHERE ID_CLASSIFICATION=$themeID");
		}
		my %oldList;
		#my $checkOld=0;
		if ($listID==-1) { # Creating new list
			#$replace=0; # no need to clean list
			$saveProtAction='add' if $saveProtAction eq 'replace'; # no need to clean list
			$listName=param('listName');
			($listID)=$dbh->selectrow_array("SELECT MAX(ID_CATEGORY) FROM CATEGORY");
			$listID++;
			my ($displayPos)=$dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM CATEGORY WHERE ID_CLASSIFICATION=$themeID");
			$displayPos++;
			my $sthInsCA=$dbh->prepare("INSERT INTO CATEGORY (ID_CATEGORY,NAME,LIST_TYPE,DISPLAY_POS,ID_CLASSIFICATION,UPDATE_USER,UPDATE_DATE) VALUES ($listID,?,?,$displayPos,?,'$userID',NOW())");
			$sthInsCA->execute($listName,$listType,$themeID);
			$sthInsCA->finish;
		}
		else { # use existing list
			$listID=&promsMod::cleanNumericalParameters($listID);
			($listName)=$dbh->selectrow_array("SELECT NAME FROM CATEGORY WHERE ID_CATEGORY=$listID");

			if ($saveProtAction eq 'replace') {
				###<Delete previous category contents
				$dbh->do("DELETE FROM MODIFICATION_SITE WHERE ID_CATEGORY=$listID");
				$dbh->do("DELETE FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$listID");
				$saveProtAction='add'; # replace -> add
			}
			else { # add/remove => fetch list content
				###<Fetching existing protein/site contents to prevent duplicates (TODO: use INSERT IGNORE instead)
				my $noSiteFlag=($listType eq 'SITE')? 0 : 1;
				&promsQuantif::fetchCustomList($dbh,$listID,\%oldList,$noSiteFlag); # 1 => no site
				#$checkOld=1 if scalar keys %oldList;
			}
		}

		if ($listType eq 'PROT') {
			###<Removing duplicates (eg. if from modification quantif) and skipping already recorded prot
			my (%newProtID,%deleteProtID);
			foreach my $modProtID (split(',',param('protID'))) {
				my ($protID,$modCode)=split('-',$modProtID);
				if ($saveProtAction eq 'add') {
					#next if ($checkOld && $oldList{$protID}); # protein is already in category
					$newProtID{$protID}=1 unless $oldList{$protID};
				}
				else { # remove
					$deleteProtID{$protID}=1 if $oldList{$protID};
				}
			}
			###<Updating protein list
			if ($saveProtAction eq 'add') {
				my $sthInsCP=$dbh->prepare("INSERT INTO CATEGORY_PROTEIN (ID_PROTEIN,ID_CATEGORY) VALUES (?,$listID)");
				foreach my $protID (keys %newProtID) {
					$sthInsCP->execute($protID);
					#$oldList{$protID}=1; # In case multiple occurences of same protID in list (eg called from listProteins.cgi)
				}
				$sthInsCP->finish;
			}
			else { # remove
				$dbh->do("DELETE FROM CATEGORY_PROTEIN WHERE ID_PROTEIN IN (".join(',',keys %deleteProtID).") AND ID_CATEGORY=$listID");
			}
		}
		elsif ($listType eq 'SITE') {
			###<Updating site list
			if ($saveProtAction eq 'add') {
				my $sthInsCP=$dbh->prepare("INSERT INTO CATEGORY_PROTEIN (ID_PROTEIN,SEQ_BEG,SEQ_LENGTH,ID_CATEGORY) VALUES (?,?,?,$listID)");
				my $sthInsMS=$dbh->prepare("INSERT INTO MODIFICATION_SITE (ID_CATEGORY,ID_CATEGORY_PROTEIN,ID_MODIFICATION,RESIDUE,POSITION) VALUES (?,?,?,?,?)");
				foreach my $modProtID (split(',',param('protID'))) {
					#next if ($checkOld && $oldList{$modProtID});
					next if $oldList{$modProtID};
					my ($protID,$allFullModCode)=split('-',$modProtID);
					my ($allModCode,$seqContext)=$allFullModCode=~/^([^\[]+)(.*)/; # sequence context-compatible
					next unless $allModCode; # in case of mixte of sites and prot (eg from compareQuantifications.cgi)
					my ($seqBeg,$seqLength)=($seqContext && $seqContext=~/(\d+)\.(\d+)/)? ($1,$2) : (undef,undef);
					$sthInsCP->execute($protID,$seqBeg,$seqLength);
					my ($catProtID)=$dbh->last_insert_id(undef,undef,'CATEGORY_PROTEIN','ID_CATEGORY_PROTEIN');
					foreach my $modifCode (split(/\+/,$allModCode)) { # SAME CODE IN editClassification.cgi for list import ====>
						my ($modID,$modCode)=$modifCode=~/^(\d+):(.+)/;
						$modID=-1 if $modID==0; # Free residues
						#<Encode modCode for DB storage
						if ($modCode=~/~/) {
							$modCode=~s/^ProtNt~/=0./; # position was absent
							$modCode=~s/^PepNt(\d+)~/=$1./;
							$modCode=~s/^(\d+)~/-$1./; # normal res
							$modCode=~s/ProtCt:/\*99999./; # position was absent
							$modCode=~s/PepCt(\d+):/\*$1./;
							$modCode=~s/(\d+):/+$1./; # normal res
							$modCode=~s/\///;
						}
						else {
							$modCode=~s/ProtNt/n0/; # position was absent
							$modCode=~s/ProtCt/c99999/; # position was absent
							$modCode=~s/PepNt/n/;
							$modCode=~s/PepCt/c/;
						}
						#<Separate individual site
						foreach my $site (split(/\./,$modCode)) {
							$site=~/(.)(\d+)/; # (res)(position)
							$sthInsMS->execute($listID,$catProtID,$modID,$1,$2);
						}
					} # <====== SAME CODE
					$sthInsMS->finish;
					$sthInsCP->finish;
				}
			}
			else { # remove
				my $sthDelMS=$dbh->prepare("DELETE FROM MODIFICATION_SITE WHERE ID_CATEGORY_PROTEIN=?");
				my $sthDelCP=$dbh->prepare("DELETE FROM CATEGORY_PROTEIN WHERE ID_CATEGORY_PROTEIN=?");
				foreach my $modProtID (split(',',param('protID'))) {
					if ($modProtID=~/-/) { # site
						next unless $oldList{$modProtID}; # not in List
						$sthDelMS->execute($oldList{$modProtID});
						$sthDelCP->execute($oldList{$modProtID});
					}
					else { # protein => remove all matching sites
						foreach my $oldModProtID (keys %oldList) {
							if ($oldModProtID=~/^$modProtID-/) {
								$sthDelMS->execute($oldList{$oldModProtID});
								$sthDelCP->execute($oldList{$oldModProtID});
							}
						}
					}
				}
				$sthDelMS->finish;
				$sthDelCP->finish;
			}
		}

		###<Updating Theme
		$dbh->do("UPDATE CLASSIFICATION SET UPDATE_USER='$userID',UPDATE_DATE=NOW() WHERE ID_CLASSIFICATION=$themeID");
		$dbh->commit;

		print qq
|<B>List <FONT color="#DD0000">$listName</FONT> of Theme <FONT color="#DD0000">$themeName</FONT> was successfully updated.</B>
<!-- THEME_ID=$themeID THEME_NAME=$themeName### LIST_ID=$listID LIST_NAME=${listName}::: (parsed by JS when needed) -->
|;

	}

	$dbh->disconnect;
	exit;
}

sub ajaxRestrictProteinList {
	my $projectID=param('projectID');
	my $submitForm=param('submitForm') || '';
	my $noSelect=param('noSelect') || 0;

	my $dbh=&promsConfig::dbConnect;
	# my $sthL=$dbh->prepare("SELECT T.ID_CLASSIFICATION,T.NAME,L.ID_CATEGORY,L.NAME,L.DISPLAY_POS,LIST_TYPE FROM CLASSIFICATION T,CATEGORY L WHERE T.ID_CLASSIFICATION=L.ID_CLASSIFICATION AND T.ID_PROJECT=$projectID");
	my $sthL=$dbh->prepare("SELECT T.ID_CLASSIFICATION,T.NAME,L.ID_CATEGORY,L.NAME,L.DISPLAY_POS,LIST_TYPE,COUNT(CP.ID_CATEGORY_PROTEIN) FROM CLASSIFICATION T,CATEGORY L,CATEGORY_PROTEIN CP WHERE T.ID_CLASSIFICATION=L.ID_CLASSIFICATION AND L.ID_CATEGORY=CP.ID_CATEGORY AND T.ID_PROJECT=$projectID GROUP BY L.ID_CATEGORY");
	$sthL->execute;
	my (%savedLists,%themeInfo);
	while (my ($themeID,$themeName,$listID,$listName,$listPos,$type,$numItems)=$sthL->fetchrow_array) {
		$type='PROT' unless $type;
		$themeInfo{$themeID}=$themeName;
		@{$savedLists{$themeID}{$listID}}=($listName,$listPos,$type,$numItems);
	}
	$sthL->finish;
	$dbh->disconnect;

	###<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
	unless ($noSelect) { # otherwise <SELECT> and 1st "-= Selec =-" option => written by JS function
		print "<SELECT name=\"restrictList\" class=\"title3\"";
		print " onchange=\"document.$submitForm.submit()\"" if $submitForm;
		if (scalar keys %savedLists) {print ">\n<OPTION value=\"\">-= Select =-</OPTION>\n";}
		else {print ">\n<OPTION value=\"\">***No List found***</OPTION>\n";}
	}
	foreach my $themeID (sort{lc($themeInfo{$a}) cmp lc($themeInfo{$b})} keys %themeInfo) {
		print "<OPTGROUP label=\"$themeInfo{$themeID}\">\n";
		foreach my $listID (sort{lc($savedLists{$themeID}{$a}[1]) cmp lc($savedLists{$themeID}{$b}[1])} keys %{$savedLists{$themeID}}) {
			print "<OPTION value=\"$listID\"";
			print ' selected' if $listID==$restrictListID; # in case a list was already selected
			my $typeStrg=($savedLists{$themeID}{$listID}[2] eq 'SITE')? 'Sites ' : '';
			print ">$savedLists{$themeID}{$listID}[0] [$typeStrg","x$savedLists{$themeID}{$listID}[3]]</OPTION>\n";
		}
		print "</OPTGROUP>\n";
	}
	print "</SELECT>\n" unless $noSelect;
	exit;
}


sub ajaxDisplayCorrelationMatrix {
#print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1); # DEBUG
	my $projectID=param('projectID');

	my (%states,%labelsName);
	my ($quantifSoftware,$algoVersion,$drawGroups,$ratioType,$hasNormalization)=&getSamplesInfo(\%states,\%labelsName);
	
	my (%labels,%labelData);
	my $minCorrelation=1;
	if ($quantifSoftware eq 'myProMS') {
		my @correlStages=(['BEFORE','Beforematrixcorrelation.txt']);
		if ($hasNormalization) {push @correlStages,['AFTER','Aftermatrixcorrelation.txt'];}
		foreach my $refStage (@correlStages) {
			my ($stage,$corFile)=@{$refStage};
			my $corrMatrixFile="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/results/$corFile";
			#my $corrMatrixFile="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/results/Afteroutmatrixcorrelation.txt";
			open(CORR,$corrMatrixFile) or die "ERROR: Cannot open file '$corrMatrixFile': $!";
			my (%convertLabels,@labelArray);
			my @rawLabelList; # Temp
			while (<CORR>) {
				chomp;
				if ($.==1) {
					my %bioRepCountShift=(A=>0);
					my %numBioRepInState;
					my $prevExpCode='';
					foreach my $rawLabel (split(/\t/,$_)) {
						next unless $rawLabel; # top left cell is empty
						push @rawLabelList,$rawLabel;
						my ($expCode,$stateNum,$bioRepNum,$techRepNum);
						if ($algoVersion==2) {
							($expCode,$stateNum,$bioRepNum,$techRepNum)=($rawLabel=~/^(\w*)-?State(\d+)-biorep(\d+)-techRep(\d+)/); # "\w-" only for SuperRatio
						}
						else { # v3+
							#($stateNum,$expCode,$bioRepNum,$techRepNum)=($rawLabel=~/^[\D]+(\d+)_([^_]+)_[\D]+(\d+)_[\D]+(\d+)/); # v3b
							($stateNum,$expCode,$bioRepNum,$techRepNum)=($rawLabel=~/^State(\d+)_([^_]+)_rep(\d+)_techRep(\d+)/); # v3b
						}

						#if ($ratioType eq 'SuperRatio' && $stateNum==1) { # Cummulate BioReps across all Exp for Ref SILAC (State1)
						if ($stateNum==1) { # Cummulate BioReps across all Exp for Ref SILAC (State1)
							my $nextCode=chr(ord($expCode)+1); # A -> B
							$bioRepCountShift{$nextCode}=$bioRepCountShift{$expCode} if $expCode ne $prevExpCode;
							$bioRepCountShift{$nextCode}++;
							$bioRepNum+=$bioRepCountShift{$expCode};
							$prevExpCode=$expCode;
						}
						my $labelCode="State$stateNum-bioRep$bioRepNum-techRep$techRepNum";
						my $label=$labelsName{$labelCode} || "*$labelCode*";
						@{$labels{$labelCode}}=($stateNum,$bioRepNum,$techRepNum,$label);
						$convertLabels{$rawLabel}=$labelCode;
						push @labelArray,$labelCode;
					}
					next;
				}
				my @data=split(/\t/,$_);
				my $rawLabel=(scalar @data > scalar keys %convertLabels)? shift @data : shift @rawLabelList; # 1st col is row label OR no raw label (assume same order than in header)
				foreach my $i (0..$#data) {
					#if ($stage eq 'BEFORE') {next if $i < $.-1;} else {next if $i >= $.-1;} # TOP-RIGHT=BEFORE, BOTTOM-LEFT=AFTER
					if ($stage eq 'BEFORE') {
						if ($i < $.-2) {
							$labelData{$convertLabels{$rawLabel}}{$labelArray[$i]}='NA'; # Initialize 'AFTER' with NA in case no Normalization
							next;
						}
					}
					else {next if $i > $.-1;} # TOP-RIGHT=BEFORE, BOTTOM-LEFT=AFTER
					$labelData{$convertLabels{$rawLabel}}{$labelArray[$i]}=$data[$i];
					$minCorrelation=$data[$i] if ($data[$i] ne 'NA' && $minCorrelation > $data[$i]);
				}
			}
			close CORR;
		}
	}
	elsif ($quantifSoftware =~ /MSstats/) {

	}

	###<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
	my $numLabels=scalar keys %labels;
	#my $matSize=($numLabels < 5)? 100 : ($numLabels < 50)?  $numLabels*20 : 1000;
	#JS: width:$matSize,height:$matSize,
	my $cellSize=($numLabels == 2)? 100 : ($numLabels <= 4)? 50 : ($numLabels <= 25)? 25 : ($numLabels <= 50)? 15 : 10;
	print qq
|var GCM=new heatMap({
div:'globCorrDIV',cellWidth:$cellSize,cellHeight:$cellSize,
// editable:true,
colLabelOrientation:'vertical',
moveLabel:true,
// normalization:{scope:'row',reference:'z-score',colors:'GYR'}
normalization:{scope:'row',reference:'user',limitValues:{ref:0,max:1}} //$minCorrelation
});
|;
	my @sortedLabels=(sort{$labels{$a}[0]<=>$labels{$b}[0] || $labels{$a}[1]<=>$labels{$b}[1] ||$labels{$a}[2]<=>$labels{$b}[2]} keys %labels);

	#>Columns
	#print "GCM.setColumns([['",join("'],['",@sortedLabels),"']]);\n";
	print "GCM.setColumns([";
	my $firstColLabel=1;
	foreach my $colLabelCode (@sortedLabels) {
		print ',' unless $firstColLabel;
		# print "['$labels{$colLabelCode}[3]',null,'$states{$labels{$colLabelCode}[0]}[0]: $labels{$colLabelCode}[3]']";
		my $repLabel=($states{$labels{$colLabelCode}[0]}[1] > 1)? $labels{$colLabelCode}[3] : $states{$labels{$colLabelCode}[0]}[0]; # Display State if only 1 replic
		print "['$repLabel',null,'$states{$labels{$colLabelCode}[0]}[0]: $labels{$colLabelCode}[3]']";
		$firstColLabel=0;
	}
	print "]);\n";

	#>Rows
	foreach my $rowLabelCode (@sortedLabels) {
		# print "GCM.addRow(['$labels{$rowLabelCode}[3]',null,'$states{$labels{$rowLabelCode}[0]}[0]: $labels{$rowLabelCode}[3]'],[";
		my $repLabel=($states{$labels{$rowLabelCode}[0]}[1] > 1)? $labels{$rowLabelCode}[3] : $states{$labels{$rowLabelCode}[0]}[0]; # Display State if only 1 replic
		print "GCM.addRow(['$repLabel',null,'$states{$labels{$rowLabelCode}[0]}[0]: $labels{$rowLabelCode}[3]'],[";
		my $firstLabel=1;
		foreach my $colLabelCode (@sortedLabels) {
			print ',' unless $firstLabel;
			print $labelData{$rowLabelCode}{$colLabelCode} if $labelData{$rowLabelCode}{$colLabelCode} ne 'NA';
			$firstLabel=0;
		}
		print "]);\n";
	}
	if ($drawGroups) {
		my $groupStrg='';
		my $prevStatePos=0;
		my $begGrIdx=0;
		foreach my $i (0..$#sortedLabels) {
			my $statePos=$labels{$sortedLabels[$i]}[0];
			if ($i > 0 && $statePos > $prevStatePos) {
				if ($states{$prevStatePos}[1] > 1) { # No group if State has only 1 replic
					my $endGrIdx=$i-1;
					$groupStrg.=',' if $groupStrg;
					$groupStrg.="['$states{$prevStatePos}[0]',$prevStatePos,,$begGrIdx,$endGrIdx]";
				}
				$begGrIdx=$i;
			}
			$prevStatePos=$statePos;
		}
		if ($states{$prevStatePos}[1] > 1) {
			$groupStrg.=',' if $groupStrg;
			$groupStrg.="['$states{$prevStatePos}[0]',$prevStatePos,,$begGrIdx,$#sortedLabels]";
		}
		print qq
|GCM.defineGroups('row',[$groupStrg]);
GCM.defineGroups('column',[$groupStrg]);
|;
	}
	print "GCM.draw();\n";
	exit;
}

sub ajaxDisplayNormPeptides {
#print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1); # DEBUG
	my $projectID=param('projectID');
	my (%states,%labelsName);
	my ($quantifSoftware,$algoVersion,$replicInStates,$ratioType,$biasCorrection,$quantifAnnot)=&getSamplesInfo(\%states,\%labelsName);
	my $pepNormType=($quantifAnnot=~/PEPTIDES_NORM=manual/)? 'biasCorr' : ($quantifAnnot=~/PEPTIDES_REF=manual/)? 'protCorr' : '';
	unless ($pepNormType) {
		print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
		print qq
|<CENTER>
<FONT class="title3">Normalizing peptides</FONT>&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'"><BR><BR>
<B>Sorry, no data recorded for this quantification.</B>
<BR><BR>
</CENTER>
|;
		exit;
	}

	my %labels;
	foreach my $labelCode (keys %labelsName) {
		my ($stateNum,$bioRepNum,$techRepNum)=$labelCode=~/^State(\d+)-bioRep(\d+)-techRep(\d+)/;
		@{$labels{$labelCode}}=($stateNum,$bioRepNum,$techRepNum);
	}

	my (%peptideList,%protAnalyses);
	if ($pepNormType eq 'protCorr') { # Protein-level correction
		open (DATA_REF,"$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/data/tableRef.txt");
		while(<DATA_REF>) {
			next if $.==1;
			chomp;
			my ($modProtID,$fullPepCode,$state,$treatment,$biorep,$techRep,$exp,$obsKey,$pepIdKey,$protIdent,$validity,$xicValue)=split(/\t/,$_);
			my ($protID)=($modProtID=~/^(\d+)/);
			my @pepIDs=split(/\D/,$pepIdKey);
			my ($pepCode,$anaID)=split('_#',$fullPepCode); # seqModcharge
			$protAnalyses{$protID}{$anaID}=1;
			$biorep='bio'.ucfirst($biorep); # repX => bioRepX
			$peptideList{$protID}{$pepCode}{XIC}{"$state-$biorep-$techRep"}=$xicValue;
			push @{$peptideList{$protID}{$pepCode}{INFO}[0]},@pepIDs;
		}
		close DATA_REF;
	}
	else { # Bias correction
		#open (PROT_NORM,"$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/data/normProtein.txt");
		#while(<PROT_NORM>) {
		#	next if $.==1;
		#	chomp;
		#	%{$peptideList{$_}}=(); # protID, NOT a modProtID!
		#}
		#close PROT_NORM;
		open (DATA,"$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/data/table.txt");
		while(<DATA>) {
			next if $.==1;
			chomp;
			my ($protID,$fullPepCode,$state,$treatment,$biorep,$techRep,$exp,$obsKey,$pepIdKey,$protIdent,$validity,$xicValue)=split(/\t/,$_);
			#next unless $peptideList{$protID}; # non-modified protID(s) is/are expected at the end of file
			next if $protID=~/^\d+-/; # sites
			my @pepIDs=split(/\D/,$pepIdKey);
			my ($pepCode,$anaID)=split('_#',$fullPepCode); # seqModcharge
			$protAnalyses{$protID}{$anaID}=1;
			$biorep='bio'.ucfirst($biorep); # repX => bioRepX
			$peptideList{$protID}{$pepCode}{XIC}{"$state-$biorep-$techRep"}=$xicValue;
			push @{$peptideList{$protID}{$pepCode}{INFO}[0]},@pepIDs;
		}
		close DATA;
	}

	my $dbh=&promsConfig::dbConnect;
	my $sthProt=$dbh->prepare("SELECT ID_PROTEIN,ALIAS,PROT_DES,PROT_LENGTH,ORGANISM FROM PROTEIN WHERE ID_PROTEIN IN (".join(',',keys %peptideList).")");
	my %protInfo;
	$sthProt->execute;
	while (my ($protID,$alias,$des,$length,$organism)=$sthProt->fetchrow_array) {@{$protInfo{$protID}}=($alias,$des,$length,$organism);}
	$sthProt->finish;

	my (%varModName);
	foreach my $protID (keys %peptideList) {
		my (@distPeptides,%pepID2pepCode);
		foreach my $pepCode (keys %{$peptideList{$protID}}) {
			my ($seqModCode,$charge)=split('_',$pepCode); #  KSAPATGGVKKPHR&5:10&8:1_2_#3625
			my $seqVarMod;
			if ($seqModCode=~/(\w+)&(.+)/) { # with PTM(s)
				$seqVarMod=$1.' + '.&promsMod::decodeVarMod($dbh,$1,$2,\%varModName);
			}
			else { # no PTM
				$seqVarMod=$seqModCode;
			}
			push @{$peptideList{$protID}{$pepCode}{INFO}},($seqVarMod,$charge);
			my $firstPepID=$peptideList{$protID}{$pepCode}{INFO}[0][0];
			push @distPeptides,$firstPepID;
			$pepID2pepCode{$firstPepID}=$pepCode;
		}
		my $sthPep=$dbh->prepare("SELECT P.ID_PEPTIDE,ABS(PEP_BEG),PPA.IS_SPECIFIC,MISS_CUT
									FROM PEPTIDE P
									INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
									WHERE PPA.ID_PROTEIN=$protID AND P.ID_PEPTIDE IN (".join(',',@distPeptides).")");
		$sthPep->execute;
		while (my ($pepID,$pepBeg,$isSpecif,$missCut)=$sthPep->fetchrow_array) { # only 1 pepID per pepCode
			push @{$peptideList{$protID}{ $pepID2pepCode{$pepID} }{INFO}},($pepBeg,$isSpecif,$missCut); # (\@pepIDs,$seqVarMod,$charge,$pepBeg,$isSpecif,$missCut)
		}
	}
	$dbh->disconnect;

	###<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
	print qq
|<CENTER><FONT class="title3">Normalizing peptides</FONT>&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'"></CENTER>
<BR>
<TABLE cellspacing=0>
|;
	my $numStates=scalar keys %states;
	my $numSamples=scalar keys %labels;
	my @sortedSamples=sort{$labels{$a}[0]<=>$labels{$b}[0] || $labels{$a}[1]<=>$labels{$b}[1] || $labels{$a}[2]<=>$labels{$b}[2]} keys %labels; # State > bioRep > techRep
	my $numCols=7+$numSamples;
	my $numProt=scalar keys %peptideList;
	my $protCount=0;
	foreach my $protID (keys %peptideList) {
		$protCount++;
		my $anaStrg=join(',',keys %{$protAnalyses{$protID}});
		my ($alias,$protDes,$protLength,$organism)=@{$protInfo{$protID}};
		print qq
|<TR bgcolor="$darkColor"><TD class="bBorder" colspan=$numCols><TABLE>
	<TR><TH valign=top>#$protCount.<A href="javascript:sequenceView($protID,'$anaStrg')">$alias</A>:</TH><TD bgcolor="$lightColor" width=100%>$protDes <FONT class="org">$organism</FONT> ($protLength aa)</TD></TR>
	</TABLE></TD></TR>
<TR><TH width=20 rowspan=2>&nbsp;&nbsp;</TH>
<TH bgcolor="$darkColor" class="rbBorder" width=25 rowspan=2>#</TH>
<TH bgcolor="$darkColor" class="rbBorder" width=50 rowspan=2>&nbsp;Start&nbsp;</TH>
<TH bgcolor="$darkColor" class="rbBorder" rowspan=2>Peptide&nbsp;&nbsp;</TH>
<TH bgcolor="$darkColor" class="rbBorder" width=50 rowspan=2>&nbsp;Charge&nbsp;</TH>
<TH bgcolor="$darkColor" class="rbBorder" width=100 rowspan=2>&nbsp;Proteotypic&nbsp;</TH>
<TH bgcolor="$darkColor" class="rbBorder" width=60 rowspan=2>&nbsp;Miss.&nbsp;cut&nbsp;</TH>
|;
		#<States
		foreach my $stateRank (sort{$a<=>$b} keys %states) {
			my $cellClass=($stateRank==$numStates)? 'bBorder' : 'rbBorder';
			print "<TH bgcolor=\"$darkColor\" class=\"$cellClass\" colspan=\"$states{$stateRank}[1]\">$states{$stateRank}[0]</TH>\n";
		}
		print "</TR><TR>\n";
		#<Replicates
		my $sampCount=0;
		foreach my $labelCode (@sortedSamples) {
			$sampCount++;
			my $cellClass=($sampCount==$numSamples)? 'bBorder' : 'rbBorder';
			print "<TH bgcolor=\"$darkColor\" class=\"$cellClass\">&nbsp;$labelsName{$labelCode}&nbsp;</TH>\n";
		}
		print "</TR>\n";

		my $bgColor=$lightColor;
		my $numPep=0;
		foreach my $pepCode (sort{$peptideList{$protID}{$a}{INFO}[3]<=>$peptideList{$protID}{$b}{INFO}[3] || $peptideList{$protID}{$a}{INFO}[2]<=>$peptideList{$protID}{$b}{INFO}[2] || $peptideList{$protID}{$a}{INFO}[5]<=>$peptideList{$protID}{$b}{INFO}[5]} keys %{$peptideList{$protID}}) { # pepBeg || charge || missCut
			$numPep++;
			my ($refPepIDs,$seqVarMod,$charge,$pepBeg,$isSpecif,$missCut)=@{$peptideList{$protID}{$pepCode}{INFO}};
			#my $pepIdStrg=$protID.':'.join(',',@{$refPepIDs});
			my $specifStrg=($isSpecif)? 'yes' : 'no';
			my $missCutStrg=$missCut || 'no';
			print qq
|<TR>
	<TD></TD>
	<TD bgcolor="$bgColor" class="rBorder" align=right>&nbsp;$numPep&nbsp;</TD>
	<TD bgcolor="$bgColor" align=right>$pepBeg&nbsp;</TD>
	<TH bgcolor="$bgColor" class="font11" align=left nowrap>$seqVarMod&nbsp;</TH>
	<TD bgcolor="$bgColor" align=center>$charge<SUP>+</SUP></TD>
	<TD bgcolor="$bgColor" align=center>$specifStrg</TD>
	<TD bgcolor="$bgColor" align=center>$missCutStrg</TD>
|;
			foreach my $labelCode (@sortedSamples) {
				my $valueStrg=(defined $peptideList{$protID}{$pepCode}{XIC}{$labelCode})? sprintf "%.2e",$peptideList{$protID}{$pepCode}{XIC}{$labelCode} : '-';
				print "<TD align=\"right\" bgcolor=\"$bgColor\">&nbsp;$valueStrg&nbsp;</TH>\n";
			}
			print "</TR>\n";
			$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
		}
		print "<TR><TD colspan=8>&nbsp;</TD></TR>\n" if $protCount < $numProt;
	}
	print "</TABLE><BR><BR>\n";
	exit;
}

sub ajaxDisplayMissingValues {

	my $projectID=param('projectID');
	
	my (%states,%labelsName);
	my ($quantifSoftware,$algoVersion,$replicInStates,$ratioType)=&getSamplesInfo(\%states,\%labelsName);
	
	my (%labels,%missingValues,%missingValuesProt);
	if ($quantifSoftware =~ /myProMS|MSstats/) {
		my $resultDir="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/results";
		my $fileName=($quantifSoftware eq 'myProMS')? 'Brut_percentageNAsample.txt' : 'brut_percentageNAsample.txt';
		unless (-e "$resultDir/$fileName") {
			print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
			print qq
|<CENTER>
<FONT class="title3">Missing values (%)</FONT>&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'"><BR><BR>
<B>Sorry, missing values were not recorded for this quantification.</B>
<BR><BR>
</CENTER>
|;
			exit;
		}
		open(MISS,"$resultDir/$fileName");
		
		if ($quantifSoftware eq 'myProMS') {
			if ($algoVersion>=3) {
				while (<MISS>) {
					chomp;
					my ($sampLabel,$pcMissValue)=split(/\t/,$_); # A/State1/rep1/techRep1
					my ($expCode,$stateNum,$bioRepNum,$techRepNum)=$sampLabel=~/^(\w+)\/State(\d+)\/rep(\d+)\/techRep(\d+)/;
					my $labelCode="State$stateNum-bioRep$bioRepNum-techRep$techRepNum";
					my $label=$labelsName{$labelCode} || "*$labelCode*";
					@{$labels{$labelCode}}=($stateNum,$bioRepNum,$techRepNum,$label);
					$missingValues{$labelCode}=$pcMissValue;
				}
				close MISS;
			}
			elsif ($algoVersion==2) {
				my %lines;
				while (<MISS>) {
					chomp;
					@{$lines{$.}}=split(/\t/,$_);
				}
				close MISS;
				foreach my $i (0..$#{$lines{1}}) { # A_State1_rep1_techRep1
					my ($expCode,$stateNum,$bioRepNum,$techRepNum)=($lines{1}[$i]=~/^(\w+)_State(\d+)_rep(\d+)_techRep(\d+)/);
					my $labelCode="State$stateNum-bioRep$bioRepNum-techRep$techRepNum";
					my $label=$labelsName{$labelCode} || "*$labelCode*";
					@{$labels{$labelCode}}=($stateNum,$bioRepNum,$techRepNum,$label);
					$missingValues{$labelCode}=$lines{2}[$i];
				}
			}

			##<Abundance protein-level missing values
			my %fileTypes=('LFQ'=>'LFQ','Other'=>'');
			foreach my $type (keys %fileTypes) {
				if (-e "$resultDir/percentageNAProt$fileTypes{$type}.txt") {
					open(MISSPROT,"$resultDir/percentageNAProt$fileTypes{$type}.txt");
					while (<MISSPROT>) {
						next if $.==1;
						chomp;
						my ($state,$pcValue)=split(/\t/,$_);
						$missingValuesProt{$type}{$state}=$pcValue;
					}
					close MISSPROT;
				}
			}
		}
		elsif ($quantifSoftware =~ /MSstats/) {
			my @missData=<MISS>;
			close MISS;
			my @replic=split(/\t/,$missData[0]);
			my @missVal=split(/\t/,$missData[1]);
			my %techRep;
			foreach my $i (0..$#replic) {
				my ($stateNum,$bioRepNum)=$replic[$i]=~/State(\d+)\.BioRep(\d+)/;
				#my $techRepNum=1; # for now...
				my $techRepNum=++$techRep{"$stateNum:$bioRepNum"};
				my $labelCode="State$stateNum-bioRep$bioRepNum-techRep$techRepNum";
				my $label=$labelsName{$labelCode} || "*$labelCode*";
				@{$labels{$labelCode}}=($stateNum,$bioRepNum,$techRepNum,$label);
				$missingValues{$labelCode}=$missVal[$i];
			}
		}
	}
	
	###<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
	print qq
|<CENTER><FONT class="title3">Missing values (%)</FONT>&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'"></CENTER>
<BR>
<TABLE cellspacing=0>
|;
	my $spaceColSpan=2;
	my $numProtTypes=scalar keys %missingValuesProt;
	if ($numProtTypes) {
		print "<TR><TH class='darkBg rBorder'>State</TH>";
		if ($replicInStates) {
			print "<TH class='darkBg rBorder'>Replicate</TH>";
			$spaceColSpan++;
		}
		print "<TH class='darkBg rBorder'>&nbsp;Ions&nbsp;</TH>";
		$spaceColSpan+=$numProtTypes;
		my ($otherClass,$lfqClass)=($numProtTypes==2)? ('rBorder',''): ('','');
		print "<TH class='darkBg $otherClass'>Abundance</TH>" if $missingValuesProt{'Other'};
		print "<TH class='darkBg $lfqClass'>LFQ</TH>" if $missingValuesProt{'LFQ'};
		print "</TR>\n";
	}
	my $prevStateNum=0;
	my $bgClass='lightBg';
	foreach my $labelCode (sort{$labels{$a}[0]<=>$labels{$b}[0] || $labels{$a}[1]<=>$labels{$b}[1] ||$labels{$a}[2]<=>$labels{$b}[2]} keys %labels) {
		my $statePos=$labels{$labelCode}[0];
		my $dispValue=1*(sprintf "%.1f",$missingValues{$labelCode});
		my ($stateColStrg,$spaceRowStrg)=($replicInStates && $prevStateNum < $statePos)? ("<TH class=\"darkBg rBorder\" align=\"right\" rowspan=$states{$statePos}[1] nowrap>&nbsp;$states{$statePos}[0]&nbsp;</TH>","<TR><TD colspan=$spaceColSpan></TD></TR>\n") : ('','');
		print qq
|$spaceRowStrg<TR>
	$stateColStrg<TH class="$bgClass" align="right" nowrap>&nbsp;$labels{$labelCode}[3] :</TH>
	<TD class="$bgClass" width="100"><DIV class="barContainerLeft"><DIV class="barElementLeft barRed" style="width:$dispValue%"></DIV><DIV class="barValue">&nbsp;$dispValue&nbsp;</DIV></DIV></TD>
|;
		my $protClass=($replicInStates)? 'darkBg' : $bgClass;
		if ($prevStateNum < $statePos) {
			if ($missingValuesProt{'Other'}) {
				my $dispValueProt=1*(sprintf "%.1f",$missingValuesProt{'Other'}{"State$statePos"});
				print qq |<TD class="$protClass" width="100" rowspan=$states{$statePos}[1] nowrap><DIV class="barContainerLeft"><DIV class="barElementLeft barRed" style="width:$dispValueProt%"></DIV><DIV class="barValue">&nbsp;$dispValueProt&nbsp;</DIV></DIV></TD>\n|;
			}
			if ($missingValuesProt{'LFQ'}) {
				my $dispValueProt=1*(sprintf "%.1f",$missingValuesProt{'LFQ'}{"State$statePos"});
				print qq |<TD class="$protClass" width="100" rowspan=$states{$statePos}[1] nowrap><DIV class="barContainerLeft"><DIV class="barElementLeft barRed" style="width:$dispValueProt%"></DIV><DIV class="barValue">&nbsp;$dispValueProt&nbsp;</DIV></DIV></TD>\n|;
			}
		}
		print "</TR>\n";
		$prevStateNum=$labels{$labelCode}[0];
		$bgClass=($bgClass eq 'lightBg')? 'darkBg' : 'lightBg';
	}
	print "</TABLE>\n";
	exit;
}

sub ajaxDisplayBiasCoefficientsList {
	
	my $biasLevel=param('bias_level') || 'ion';
	my $projectID=param('projectID');
	
	my (%states,%labelsName,@replicOrder);
	my ($quantifSoftware,$algoVersion,$replicInStates,$ratioType,$hasBiasCorr,$quantifAnnot)=&getSamplesInfo(\%states,\%labelsName,\@replicOrder);
	my ($algoType)=($quantifAnnot=~/ALGO_TYPE=([^:]+)/); $algoType='' unless $algoType;
	my $refStateStrg=($algoType eq 'PEP_RATIO')? '/'.$states{1}[0] : ''; # SILAC reference
	if ($replicInStates && $algoType eq 'PEP_RATIO') { # check if replicates in Test states
		$replicInStates=-1; # flag for relic in reference SIALC only
		my $numStates=scalar keys %states;
		foreach my $statePos (2..$numStates) {
			if ($states{$statePos}[1] > 1) {
				$replicInStates=1;
				last;
			}
		}
	}
	
	my %biasData;
	if ($algoVersion >= 3) {
		my $biasPath="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/results"; # default
		if ($biasLevel ne 'ion') {
			$biasPath.="/resultsNorm/$biasLevel/results";
		}
		&fetchBiasData($biasPath,\%biasData,'long');
	}
	elsif ($quantifAnnot=~/BIAS_COEFF=([^:]+)/) {
		my @coeffData=split(';',$1);
		my $coeffIdx=-1;
		foreach my $replicIdx (0..$#replicOrder) {
			$biasData{$replicOrder[$replicIdx]}[0]=$coeffData[$replicIdx]; # [0] because only 1 bias correction value (median shift?)
		}
	}
	
	###<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
	my ($itemType,$indexShift)=($algoType eq 'PEP_RATIO')? ('Ratio',1) : ('State',0);
	print "<TABLE cellspacing=0>\n<TR bgcolor=\"$darkColor\"><TH class=\"rbBorder\" colspan=2>&nbsp;$itemType&nbsp;</TH>";
	print "<TH class=\"rbBorder\">&nbsp;Replicate&nbsp;</TH>" if $replicInStates==1;
	my $coeffHelpStrg=($algoVersion >= 3)? 'Median shift ~ Distribution stretch' : 'Median shift';
	print "<TH class=\"bBorder\">&nbsp;Bias coeff.<SUP onmouseover=\"popup('$coeffHelpStrg')\" onmouseout=\"popout()\">?</SUP>&nbsp;</TH></TR>\n";
	my $bgColor=$darkColor;
	my $prevState=0;
	foreach my $replicCode (@replicOrder) {
		my ($statePos)=$replicCode=~/State(\d+)/;
		next if ($algoType eq 'PEP_RATIO' && $statePos==1); # SILAC reference
		my ($itemPosStrg,$stateLabelStrg)=('','');
		if ($statePos != $prevState) {
			$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
			$itemPosStrg='&nbsp;#'.($statePos-$indexShift).'&nbsp;';
			$stateLabelStrg="$states{$statePos}[0]$refStateStrg:";
		}
		print "<TR bgcolor=\"$bgColor\"><TH align=\"right\">$itemPosStrg</TH><TH align=\"right\">&nbsp;";
		if ($replicInStates) {
			print "$stateLabelStrg</TH>";
			print "<TD>$labelsName{$replicCode}</TD>" if $replicInStates==1;
		}
		else {print "$labelsName{$replicCode}$refStateStrg:</TH>";}
		# $replicLabel=~s/["',]/./g;
		print "<TD>";
		if ($algoVersion >= 3) {
			printf "1/%.3f",1/$biasData{$replicCode}[0];
			printf " ~ %.3f",$biasData{$replicCode}[1];
		}
		elsif ($biasData{$replicCode}[0]) {printf "1/%.3g",$biasData{$replicCode}[0];}
		else {print "1/?";}
		print "&nbsp;</TD></TR>\n";
		$prevState=$statePos;
	}
	print "</TABLE>\n";
	exit;
}

sub ajaxDisplayBiasCoefficientsPlot {
	
	my $biasLevel=param('bias_level') || 'ion';
	my $projectID=param('projectID');
	
	my (%states,%labelsName);
	my ($quantifSoftware,$algoVersion,$replicInStates,$ratioType,$hasBiasCorr,$quantifAnnot)=&getSamplesInfo(\%states,\%labelsName);
	my ($algoType)=($quantifAnnot=~/ALGO_TYPE=([^:]+)/); $algoType='' unless $algoType;
	my $refStateStrg=($algoType eq 'PEP_RATIO')? '/'.$states{1}[0] : ''; # SILAC reference
	
	my $biasPath="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/results"; # default
	if ($biasLevel ne 'ion') {
		$biasPath.="/resultsNorm/$biasLevel/results";
	}
	my %biasData;
	&fetchBiasData($biasPath,\%biasData,'long');
	
	###<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
	print qq
|var BCP=new genericPlot({div:'coeffPlotDIV',width:350,height:350,
	axisX:{title:'Log2(Median shift)'},
	axisY:{title:'Log2(Distribution stretch)'},
	axisClosure:true,
	zoomable:true,
	pointLabel:clpPointLabel,
	allowHighlight:true,
	editHighlight:false,
	searchable:true, //{text:'Extended search',externalSearch:ajaxSearchConvertIdentifier},
	exportAsImage:['Export as image','BiasCoefficientsPlot','./exportSVG.cgi']
});
BCP.addThreshold({axis:'X',value:0,color:'#AAA',pattern:'- ',editable:false});
BCP.addThreshold({axis:'Y',value:0,color:'#AAA',pattern:'- ',editable:false});
BCP.addDataSet(0,{name:'Bias correction'});
|;
	print "BCP.addDataAsString(0,'";
	my (%numBioRep,%numTechRep);
	my ($maxStateBioRep,$maxStateTechRep)=(0,0);
	my $others=0;
	my $log2=log(2);
	foreach my $replicCode (keys %biasData) {
		#$replicCode=~s/biorep/bioRep/;
		my ($statePos)=$replicCode=~/State(\d+)/;
		my $replicLabel=($replicInStates)? "$states{$statePos}[0]$refStateStrg: $labelsName{$replicCode}" : $labelsName{$replicCode}.$refStateStrg;
		$replicLabel=~s/["',]/./g;
		my ($coefAdd,$coefMul)=@{$biasData{$replicCode}};
		print ';' if $others;
		print "$replicLabel,$replicCode,",log($coefAdd)/$log2,",",log($coefMul)/$log2;
		$others=1;
		$numBioRep{$statePos}{$replicCode}=1;
		$numTechRep{$statePos}{$replicCode}=1;
		$maxStateBioRep=scalar keys %{$numBioRep{$statePos}} if $maxStateBioRep < scalar keys %{$numBioRep{$statePos}};
		$maxStateTechRep=scalar keys %{$numTechRep{$statePos}} if $maxStateTechRep < scalar keys %{$numTechRep{$statePos}};
	}
	print "');\n";
	print qq
|
BCP.draw();
function clpPointLabel(dp,type) {
	var infoStrg=dp.label;
	if (type != 'min') {
		if (dp.dataSet.params.chart.dataSets.length > 1) {infoStrg+='\\nState: '+dp.dataSet.params.name;}
		// infoStrg+='\\nMedian shift=1/'+1*(dp.getX()).toFixed(3)+'\\nDist. stretch='+1*(dp.getY()).toFixed(3);
		infoStrg+='\\nMedian shift=1/'+(2**(1*(dp.getX()))).toFixed(3)+'\\nDist. stretch='+(2**(1*(dp.getY()))).toFixed(3);
	}
	return infoStrg;
}
|;
	my $numStates=scalar keys %states;
	if ($replicInStates && $numStates <= 20) {
		my $colorIdx=-1;
		my $refRep=($maxStateBioRep > $numStates)? \%numBioRep : \%numTechRep; # bioRep or techRep
		foreach my $statePos (sort{$a<=>$b} keys %states) {
			next if ($algoType eq 'PEP_RATIO' && $statePos==1);
			$states{$statePos}[0]=~s/[',]/./g;
			$colorIdx++;
			my $repArrayListStrg="['".join("','",keys %{$refRep->{$statePos}})."']";
			print qq |addHighlighting(BCP,'$states{$statePos}[0]$refStateStrg',hColors[$colorIdx % hColors.length],{'-1':$repArrayListStrg});\n|;
		}
	}
	exit;
}

sub fetchBiasData {
	my ($resultDir,$refBiasData,$format)=@_;
	$format='short' unless $format;
	open(BIAS,"$resultDir/allbias_coef.txt");
	my (%header2Idx,%biasType);
	my $lastStateIdx=-1;
	while (<BIAS>) {
		chomp;
		my ($rowName,@lineData)=split(/\s+/,$_);
		if ($.==1) {
			my $colIdx=-1;
			foreach my $colName (@lineData) {
				$colIdx++;
				$header2Idx{$colName}=$colIdx;
				$lastStateIdx++ if $colName=~/State/; # not yet in Experiment section
			}
		}
		@{$biasType{$rowName}}=@lineData; # includes header
	}
	close BIAS;

	foreach my $i (0..$lastStateIdx) {
		my ($exp,$statePos,$bioRepPos,$techRepPos)=$biasType{'run'}[$i]=~/^([^_]+)_State(\d+).*_rep(\d+)_techRep(\d+)/;
		my $replicCode=($format eq 'short')? "$statePos:$bioRepPos" : "State$statePos-bioRep$bioRepPos-techRep$techRepPos";
		foreach my $biasName ('coefAdd','coefMul') {
			my $biasValue=(!$biasType{$biasName})? 1 : ($biasName eq 'coefAdd')? 2**($biasType{$biasName}[$i]+$biasType{$biasName}[$header2Idx{$exp}]) : ($biasName eq 'coefAdd')? 1/(2**($biasType{$biasName}[$i]+$biasType{$biasName}[$header2Idx{$exp}])) : $biasType{$biasName}[$i]*$biasType{$biasName}[$header2Idx{$exp}];
			push @{$refBiasData->{$replicCode}},$biasValue;
		}
	}
}

sub ajaxDisplayCustomProteins {

	my $projectID=param('projectID');
	my (@quantifNames, @customProteins, @customProtIDs, @customQty);

	unless (-e "$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/quantif_info.txt") {
		print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
		print qq
|<CENTER>
<B>Sorry, information about your proteins is not available for this quantification.</B>
<BR><BR>
</CENTER>
|;
		exit;
	}
	
	### Retrieve custom prots names/identifiers/etc. and quantities from quantif_info
	open(INFO,"$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/quantif_info.txt");
	while (<INFO>) {
		chomp;
		$_=~s/^\s+//;
		my ($label, $value);
		if (/quantif_names/) {
			($label, $value)=split(/\t/,$_);
			@quantifNames=split(';', $value);
		} elsif (/custom_proteins/) {
			($label, $value)=split(/\t/,$_);
			@customProteins=split(';', $value);
		} elsif (/custom_prot_qty/) {
			($label, $value)=split(/\t/,$_);
			@customQty=split(';', $value);
		}
	}
	my $colNb=(scalar @customQty > 9)? 10 : (scalar @customQty);
	my (@customProtUsed, @customProtNotUsed);
	open(PROT_ID_MAP, "$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/data/matrix_in.tsv");
	while (<PROT_ID_MAP>) {
		my @fields=split(/\t/, $_);
		if (grep(/^$fields[1]$/, @customProteins)) {
			push @customProtIDs, $fields[0];
			push @customProtUsed, $fields[1];
		}
	}
	#my $customProtIdStrg=join(',', @customProtIDs);
	for my $custProt (@customProteins) {
		if (!grep(/$custProt/, @customProtUsed)) {
			push @customProtNotUsed, $custProt;
		}
	}
	
	###<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
	print qq
|<BR>
<DIV align="center">
<TABLE bgcolor="$darkColor" width=100%>
	<TR>
		<TH></TH><TH align=center colspan=$colNb>Total amounts</TH>
	</TR>
	<TR>
		<TH align="left" nowrap>Parent quantification name</TH>
|;
	for my $quantifName (@quantifNames) {
		print qq
|		<TD align="center" bgcolor="$lightColor">$quantifName</TD>
|;
	}		
	print qq
|	</TR>
	<TR>
		<TH align="left" nowrap>Quantity per cell (pg)</TH>
|;
	for my $custQty (@customQty) {
		print qq
|		<TD align="center" bgcolor="$lightColor">$custQty</TD>
|;
	}		
	print qq
|	</TR>
</TABLE>
<BR>
|;
	&ajaxListSelectedProteins(\@customProtIDs,[1]);
	print qq
|<TABLE bgcolor="$lightColor" border=0 width=100%>
	<TR>
		<TH align="center" colspan=15 bgcolor="$darkColor">Custom proteins not found in analysis / quantifications, or not used (because of match group for example)</TH>
	</TR>
	<TR>
|;
	my $count=0;
	for my $custProt (@customProtNotUsed) {
		print qq
|		<TD bgcolor="$lightColor">&nbsp;$custProt&nbsp;</TD>
|;
		$count++;
		if ($count==15) {
			$count=0;
			print "</TR><TR>"
		}
	}
	print qq
|	</TR>
	<BR>
</TABLE>
|;

	print "</DIV>";
	exit;
}

sub getSamplesInfo {
	my ($refStates,$refLabelsName,$refReplicOrder)=@_; #$refReplicOrder optional
	
	my $dbh=&promsConfig::dbConnect;
	my $sthQuantif=$dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
	my $sthState=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
	my $sthAnaInfo=$dbh->prepare("SELECT A.NAME,S.NAME FROM ANALYSIS A,SAMPLE S WHERE A.ID_SAMPLE=S.ID_SAMPLE AND A.ID_ANALYSIS=?");

	$sthQuantif->execute($selQuantifID);
	my ($quantifAnnot0)=$sthQuantif->fetchrow_array;
	my $quantifAnnot=$quantifAnnot0;
	$quantifAnnot=~s/::/\$/g; # easier to match
	my ($labelType)=($quantifAnnot=~/LABEL=([^\$]+)/);
	$labelType=uc($labelType);
	#my ($quantifSoftware,$algoVersion)=($quantifAnnot=~/=([^=]+MSstats)/)? ($1,0) : ($quantifAnnot=~/=myProMS;(\d+)/)? ('myProMS',$1) : ('myProMS',2);
	my ($quantifSoftware,$algoVersion)=($quantifAnnot=~/SOFTWARE=([^;\$]+);*([^\$]*)/);
	$quantifSoftware='myProMS' unless $quantifSoftware;
	$quantifSoftware='MaxQuant' if $quantifSoftware eq 'MQ';
	unless ($algoVersion) {
		$algoVersion=($quantifSoftware eq 'myProMS')? 2 : 0;
	}
	my ($ratioType)=($quantifAnnot=~/RATIO_TYPE=([^\$]+)/); $ratioType='Ratio' unless $ratioType;
	my $biasCorrection=($quantifAnnot=~/BIAS_CORRECTION=FALSE/)? 0 : 1;
	my ($stateStrg)=($quantifAnnot=~/STATES=([^\$]+)/);
	$stateStrg=~s/#//g;

	my $hasBioReplicates=($stateStrg=~/\./)? 1 : 0;
	my $hasTechReplicates=($stateStrg=~/&/)? 1 : 0;
	my $hasReplicates=$hasBioReplicates || $hasTechReplicates;
	my $hasFractions=($stateStrg=~/\+/)? 1 : 0;
	my $replicInStates=0;
	my $statePos=0;
	my (%labelsInfo,%sampChannelOcc);
	foreach my $stateData (split(';',$stateStrg)) {
		$statePos++;
		my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$stateData);
		$sthState->execute($expCondID);
		($refStates->{$statePos}[0])=$sthState->fetchrow_array;
		my (%parentSampleName,%anaName,%channelName);
		my $currBioRep=0;
		my $numAllRep=0;
		foreach my $bioReplicate (split(/\./,$quantiObsIDs)) { # bio rep separator
			$currBioRep++;
			my $multiTechRep=($bioReplicate=~/&/)? 1 : 0;
			my $numTechRep=0;
			foreach my $techReplicate (split(/&/,$bioReplicate)) { # tech rep separator
				$numTechRep++;
				my $replicCode="State$statePos-bioRep$currBioRep-techRep$numTechRep";
				push @{$refReplicOrder},$replicCode if $refReplicOrder;
				if ($hasFractions) { # Uses BioRep > TechRep > Fractions hierarchy (TODO: Use BioSample if any)
					my $sampName='';
					if ($hasReplicates) {
						$sampName='BioRep'.$currBioRep if $hasBioReplicates;
						if ($hasTechReplicates) {
							$sampName.=' > ' if $sampName;
							$sampName.='TechRep'.$numTechRep;
						}
					}
					else {
						my $numFrac=scalar (split(/\+/,$techReplicate));
						$sampName='Fract.1-'.$numFrac;
					}
					$refLabelsName->{$replicCode}=$sampName;
				}
				else { # Uses Project hierarchy (Sample > Analysis) because more precise?
					my ($firstFrac)=($techReplicate=~/^([^\+]+)/);
					my @obsData=split(/:/,$firstFrac);
					my ($parQuantifID,$anaID,$targetPos);
					if (scalar @obsData==4) {($parQuantifID,$anaID,$targetPos)=@obsData[1..3];} # design post 07/2013 -> obsID,parQuantifID,anaID,targetPos
					else { # older design -> parQuantifID,anaID
						($parQuantifID,$anaID)=@obsData;
						$targetPos=0;
					}
					unless ($anaName{$anaID}) {
						$sthAnaInfo->execute($anaID);
						($anaName{$anaID},$parentSampleName{$anaID})=$sthAnaInfo->fetchrow_array;
					}
					my $targetName='';
					if ($targetPos) { # labeling
						unless ($channelName{$parQuantifID}) {
							$sthQuantif->execute($parQuantifID);
							my ($qAnnot)=$sthQuantif->fetchrow_array;
							my @quantifInfo=split(/::/,$qAnnot);
							foreach my $channelInfo (@quantifInfo[1..$#quantifInfo]) {
								my ($chanPos,$chanName)=split(/;/,$channelInfo); # + other unused elements
								$channelName{$parQuantifID}{$chanPos}=$chanName;
							}
						}
						$targetName=' > '.$channelName{$parQuantifID}{$targetPos};
					}
					@{$labelsInfo{$replicCode}}=($statePos,$parentSampleName{$anaID},$anaName{$anaID},$targetName);
					$sampChannelOcc{$parentSampleName{$anaID}.$targetName}++;
				}
			}
			$numAllRep+=$numTechRep;
		}
		$refStates->{$statePos}[1]=$numAllRep; # [0]: name, [1]: number of replicates
	}

	unless ($hasFractions) { # Analysis hierarchy
		foreach my $replicCode (keys %labelsInfo) {
			$refLabelsName->{$replicCode}=($hasReplicates || $refStates->{$labelsInfo{$replicCode}[0]}[0] eq $labelsInfo{$replicCode}[1])? $labelsInfo{$replicCode}[1] : "[$refStates->{$labelsInfo{$replicCode}[0]}[0]] $labelsInfo{$replicCode}[1]";
			my $sampChan=$labelsInfo{$replicCode}[1].$labelsInfo{$replicCode}[3];
			$refLabelsName->{$replicCode}.=($sampChannelOcc{$sampChan}==1)? $labelsInfo{$replicCode}[3] : ' > '.$labelsInfo{$replicCode}[2].$labelsInfo{$replicCode}[3];
		}
	#$labelsName{$replicCode}.=" ($replicCode)"; # DEBUG
	}

	$sthState->finish;
	$sthQuantif->finish;
	$sthAnaInfo->finish;
	$dbh->disconnect;
	
	return ($quantifSoftware,$algoVersion,$hasReplicates,$ratioType,$biasCorrection,$quantifAnnot0);
}

# sub getSamplesInfo_old { # PROBLEM: No fraction management
# 	my ($refStates,$refLabelsName)=@_;
	
# 	my $dbh=&promsConfig::dbConnect;
# 	my $sthQuantif=$dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
# 	my $sthState=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
# 	my $sthAnaInfo=$dbh->prepare("SELECT A.NAME,S.NAME FROM ANALYSIS A,SAMPLE S WHERE A.ID_SAMPLE=S.ID_SAMPLE AND A.ID_ANALYSIS=?");

# 	$sthQuantif->execute($selQuantifID);
# 	my ($quantifAnnot0)=$sthQuantif->fetchrow_array;
# 	my $quantifAnnot=$quantifAnnot0;
# 	$quantifAnnot=~s/::/\$/g; # easier to match
# 	my ($labelType)=($quantifAnnot=~/LABEL=([^\$]+)/);
# 	$labelType=uc($labelType);
# 	#my ($quantifSoftware,$algoVersion)=($quantifAnnot=~/=([^=]+MSstats)/)? ($1,0) : ($quantifAnnot=~/=myProMS;(\d+)/)? ('myProMS',$1) : ('myProMS',2);
# 	my ($quantifSoftware,$algoVersion)=($quantifAnnot=~/SOFTWARE=([^;\$]+);*([^\$]*)/);
# 	$quantifSoftware='myProMS' unless $quantifSoftware;
# 	$quantifSoftware='MaxQuant' if $quantifSoftware eq 'MQ';
# 	unless ($algoVersion) {
# 		$algoVersion=($quantifSoftware eq 'myProMS')? 2 : 0;
# 	}
# 	my ($ratioType)=($quantifAnnot=~/RATIO_TYPE=([^\$]+)/); $ratioType='Ratio' unless $ratioType;
# 	my $biasCorrection=($quantifAnnot=~/BIAS_CORRECTION=FALSE/)? 0 : 1;
# 	my ($stateStrg)=($quantifAnnot=~/STATES=([^\$]+)/);
# 	$stateStrg=~s/#//g;

# 	my $replicInStates=0;
# 	my $statePos=0;
# 	my (%labelsInfo,%sampChannelOcc);
# 	foreach my $stateData (split(';',$stateStrg)) {
# 		$statePos++;
# 		my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$stateData);
# 		$sthState->execute($expCondID);
# 		($refStates->{$statePos}[0])=$sthState->fetchrow_array;
# 		my (%parentSampleName,%anaName,%channelName);
# 		my $currBioRep=0;
# 		my $numAllRep=0;
# 		foreach my $bioReplicate (split(/\./,$quantiObsIDs)) { # bio rep separator
# 			$currBioRep++;
# 			my $multiTechRep=($bioReplicate=~/&/)? 1 : 0;
# 			my $numTechRep=0;
# 			foreach my $techReplicate (split(/&/,$bioReplicate)) { # tech rep separator
# 				$numTechRep++;
# 				#my $multiFrac=($techReplicate=~/\+/)? 1 : 0;
# 				my ($firstFrac)=($techReplicate=~/^([^\+]+)/);
# 				my @obsData=split(/:/,$firstFrac);
# 				my ($parQuantifID,$anaID,$targetPos);
# 				if (scalar @obsData==4) {($parQuantifID,$anaID,$targetPos)=@obsData[1..3];} # design post 07/2013 -> obsID,parQuantifID,anaID,targetPos
# 				else { # older design -> parQuantifID,anaID
# 					($parQuantifID,$anaID)=@obsData;
# 					$targetPos=0;
# 				}
# 				unless ($anaName{$anaID}) {
# 					$sthAnaInfo->execute($anaID);
# 					($anaName{$anaID},$parentSampleName{$anaID})=$sthAnaInfo->fetchrow_array;
# 				}
# 				my $targetName='';
# 				if ($targetPos) { # labeling
# 					unless ($channelName{$parQuantifID}) {
# 						$sthQuantif->execute($parQuantifID);
# 						my ($qAnnot)=$sthQuantif->fetchrow_array;
# 						my @quantifInfo=split(/::/,$qAnnot);
# 						foreach my $channelInfo (@quantifInfo[1..$#quantifInfo]) {
# 							my ($chanPos,$chanName)=split(/;/,$channelInfo); # + other unused elements
# 							$channelName{$parQuantifID}{$chanPos}=$chanName;
# 						}
# 					}
# 					$targetName=' > '.$channelName{$parQuantifID}{$targetPos};
# 				}
# 				@{$labelsInfo{"State$statePos-biorep$currBioRep-techRep$numTechRep"}}=($statePos,$parentSampleName{$anaID},$anaName{$anaID},$targetName);
# 				$sampChannelOcc{$parentSampleName{$anaID}.$targetName}++;
# 			}
# 			$numAllRep+=$numTechRep;
# 		}
# 		$refStates->{$statePos}[1]=$numAllRep; # [0]: name, [1]: number of replicates
# 	}

# 	$replicInStates=1 if $statePos < scalar keys %labelsInfo; # bio|tech replicates => draw label groups
# 	foreach my $label (keys %labelsInfo) {
# 		$refLabelsName->{$label}=($replicInStates)? $labelsInfo{$label}[1] : "[$refStates->{$labelsInfo{$label}[0]}[0]] $labelsInfo{$label}[1]";
# 		my $sampChan=$labelsInfo{$label}[1].$labelsInfo{$label}[3];
# 		$refLabelsName->{$label}.=($sampChannelOcc{$sampChan}==1)? $labelsInfo{$label}[3] : ' > '.$labelsInfo{$label}[2].$labelsInfo{$label}[3];
# 	#$labelsName{$label}.=" ($label)"; # DEBUG
# 	}

# 	$sthState->finish;
# 	$sthQuantif->finish;
# 	$sthAnaInfo->finish;
# 	$dbh->disconnect;
	
# 	return ($quantifSoftware,$algoVersion,$replicInStates,$ratioType,$biasCorrection,$quantifAnnot0);
# }

sub ajaxManageTemplate {
	my $templateAction=param('TEMPLATE');
	my $dbh=&promsConfig::dbConnect;
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
	if ($templateAction eq 'add') {
		$dbh->do("UPDATE QUANTIFICATION SET QUANTIF_ANNOT=CONCAT(QUANTIF_ANNOT,'::TEMPLATE=1') WHERE ID_QUANTIFICATION=$selQuantifID");
		#print "&nbsp;Yes&nbsp;&nbsp;<INPUT type=\"button\" id=\"remove\" value=\"Remove from templates\" onclick=\"manageTemplate(this.id)\">";
		print '<SPAN class="template">Template&nbsp;&nbsp;<INPUT type="button" class="font11" value="Remove" onclick="manageTemplate(\'remove\')"></SPAN>';
	}
	elsif ($templateAction eq 'remove') {
		$dbh->do("UPDATE QUANTIFICATION SET QUANTIF_ANNOT=REPLACE(QUANTIF_ANNOT,'::TEMPLATE=1','') WHERE ID_QUANTIFICATION=$selQuantifID");
		#print "&nbsp;No&nbsp;&nbsp;<INPUT type=\"button\" id=\"add\" value=\"Save as template\" onclick=\"manageTemplate(this.id)\">";
		print '<INPUT type="button" class="title3" value="Save as template" onclick="manageTemplate(\'add\')">';
	}
	$dbh->commit;
	$dbh->disconnect;
	exit;
}

#sub dbEncodeModifPosition {
#	my ($fullModStrg,$refModRanks)=@_;
#	return '' unless $fullModStrg;
#	my $isMultiModif=(!$refModRanks || scalar keys %{$refModRanks} == 1)? 0 : 1; 
#	my %rankModStrg;
#	foreach my $modStrg (split(/\+/,$fullModStrg)) {
#		my ($modID,$resStrg)=$modStrg=~/^(\d+):(.+)/;
#		my $modifRank=($isMultiModif)? $refModRanks->{$modID} : 1;
#		my $modRkStrg=($isMultiModif)? $modifRank : '';
#		if ($modStrg=~/~/) { # ambiguity format => convert to DB output format
#			my ($start,$end,$numMods,$numSites)=$resStrg=~/^(\d+)~(\d+):(\d+)\/(\d+)/;
#			#$rankModStrg{$modifRank}=$modRkStrg.':'.$numMods.$numSites.'.'.$modRkStrg.':-'.$start.'.'.$modRkStrg.':+'.$end;
#			$rankModStrg{$modifRank}="$modRkStrg:$numMods$numSites.$modRkStrg:-$start.$modRkStrg:+$end";
#		}
#		else {
#			$resStrg=~s/\./\.$modRkStrg:/g;
#			$rankModStrg{$modifRank}=$modRkStrg.':'.$resStrg;
#		}
#	}
#	my $fullDbModStrg='';
#	foreach my $modifRank (sort{$a<=>$b} keys %rankModStrg) {
#		$fullDbModStrg.='.' if $fullDbModStrg;
#		$fullDbModStrg.=$rankModStrg{$modifRank};
#	}
#	
#	return $fullDbModStrg;
#}

sub revertLog2 {
	my $l2Val=shift;
	#my $fc=exp($l2Val*log(2));
	my $fc=2**$l2Val;
	my $formatFc=($fc > 1)? sprintf "%.2f",$fc : '1/'.sprintf "%.2f",1/$fc;
	$formatFc=~s/\.0+\Z//;
	$formatFc=1 if $formatFc eq '1/1';
	return $formatFc;
}


sub ajaxDisplayProteinQuantities {  # Globals: %promsPath, $selQuantifID, $lightColor, $darkColor
	my $dbh=&promsConfig::dbConnect;
	my $projectID = param('projectID');

	my (%states, %labelsName);
	my ($quantifSoftware, $algoVersion, $replicInStates, $ratioType) = &getSamplesInfo(\%states, \%labelsName);

	my (%proteins, %initialQty, %totalQty, %protMW);
	my $protQtyFile = "$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/data/prot_quantities.csv";
	my $totalQtyFile = "$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/results/absolute/aLFQ_abs_total.txt";
	my $mwFile = "$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/data/prot_molecular_weights.tsv";

	my %unitsMatching = (  # internal code to name displayed # TODO: remove code redundancy
		"mol" 			 => "mol",
		"mmol" 			 => "mmol",
		"copy_number" 	 => "Copy Number",
		"conc_mol_l" 	 => "mmol/mL",
		"conc_mmol_l" 	 => "mmol/L",
		"mass_g" 		 => "g",
		"mass_mg" 		 => "mg",
		"mass_ug" 		 => "μg",
		"mass_conc_g_l"  => "μg/μL",
		"mass_conc_mg_l" => "μg/mL",
		"mass_conc_ug_l" => "μg/L"
	);

	# Find which case we are working with (with parameters used)
	my ($quantifAnnot) = $dbh->selectrow_array("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION = $selQuantifID");
	my ($aLFQAnnot) = $quantifAnnot =~ /ABUND_MEASURES=.*(aLFQ[^;]+)/;
	my @aLFQAnnots = split(',', $aLFQAnnot);

	my $isCalibration 	= ($aLFQAnnots[1] =~ /calibration/)? 1 : 0;
	my $hasTotalQty 	= ($aLFQAnnots[3] =~ /amounts/)? 1 : 0;
	my $hasInitQty 		= (exists($unitsMatching{$aLFQAnnots[-1]}))? 1 : 0;
	my $initialUnits 	= ($hasInitQty)? $aLFQAnnots[-1] : '';
	my $isConcentration = ($initialUnits && $initialUnits =~ /conc/)? 1 : 0;
	my $isMass			= ($initialUnits && $initialUnits =~ /mass/)? 1 : 0;
	my $qtyMultiplier;  # To convert mol, mol/L, g or g/L recorded back to initial units

	if ($hasTotalQty && -f $totalQtyFile) {  # Exists only if known quantities provided
		open(TOTAL_QTY, "$totalQtyFile");
		while (my $line = <TOTAL_QTY>) {
			next if ($. == 1);
			my ($statePos, $molAmount, $massAmount) = split("\t", $line);
			$statePos =~ s/^State//;
			$totalQty{$statePos}{'MOL'} = ($molAmount =~ /[Nn]\/?[Aa]/)? $molAmount : sprintf("%.3g", $molAmount);
			$totalQty{$statePos}{'MASS'} = ($massAmount =~ /[Nn]\/?[Aa]/)? $massAmount : sprintf("%.3g", $massAmount);
		}
		close TOTAL_QTY;
	}

	if ($initialUnits) {
		$qtyMultiplier = ($initialUnits =~ /^mmol|^conc_mmol_l|^mass_mg|^mass_conc_mg_l/)? 1e03 :
						 ($initialUnits =~ /^mass_ug|^mass_conc_ug_l/)? 1e06 :
						 ($initialUnits =~ /^copy_number/)? 6.02e+23 :
						 1;  # default is 1 for /^mol|^conc_mol_l|^mass_g|^mass_conc_g_l/
		if ($isMass && -f $mwFile) {
			open(MW_FILE, "$mwFile");
			while (my $line = <MW_FILE>) {
				next if ($. == 1);
				chomp $line;
				my @values = split("\t", $line);
				$protMW{$values[0]} = $values[1];
			}
			close MW_FILE;
		}
	}

	if ($hasInitQty && -f $protQtyFile) {  # Same, exists only if known quantities provided	
		my $sthProtAlias = $dbh->prepare("SELECT ALIAS FROM PROTEIN WHERE ID_PROTEIN = ?");  # TODO: Improve efficiency for this query by retrieving all prots at once and storing in hash
		open(INIT_QTY, "$protQtyFile");
		while (my $line = <INIT_QTY>) {
			next if ($. == 1);
			my @values = split(',', $line);
			my ($statePos, $bioRep, $techRep) = split('_', $values[0]);
			$statePos =~ s/^State//;
			$bioRep =~ s/^rep//;
			$techRep =~ s/^techRep//;
			if ($isCalibration) {
				(my $protID = $values[1]) =~ s/\"//g;  # remove quotes introduced around protID for aLFQ
				unless ($proteins{$protID}) {
					$sthProtAlias->execute($protID);
					($proteins{$protID}) = $sthProtAlias->fetchrow_array;
				}
				if ($isMass) {
					$initialQty{$statePos}{$bioRep}{$techRep}{$protID} = ($values[2] =~ /[Nn]\/?[Aa]/)? $values[2] : sprintf("%.3g", $values[2] * $protMW{$protID} * $qtyMultiplier);
				} else {
					$initialQty{$statePos}{$bioRep}{$techRep}{$protID} = ($values[2] =~ /[Nn]\/?[Aa]/)? $values[2] : sprintf("%.3g", $values[2] * $qtyMultiplier);
				}
			} else {
				$initialQty{$statePos}{$bioRep}{$techRep} = ($values[1] =~ /[Nn]\/?[Aa]/)? $values[1] : sprintf("%.3g", $values[1] * $qtyMultiplier);
			}
		}
		close INIT_QTY;
		$sthProtAlias->finish;
	}

	# Start HTML
	print header(-type=>'text/plain',-charset=>'UTF-8');
	warningsToBrowser(1);
	print qq
|<CENTER>
<FONT class="title3">Protein quantities per sample</FONT>
&nbsp;&nbsp;
<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'">
</CENTER>
<BR>
|;

	unless (%totalQty || %initialQty) {
		print qq
|<CENTER>
<B>Information about the protein quantities is not available for this quantification.</B>
<BR><BR>
</CENTER>
|;
		exit;
	}

	my $headerRowNb = ($isCalibration)? 2 : 1;
	my $headerColNb = ($isCalibration)? scalar keys %proteins : 1;
	
	print "<CENTER><TABLE border=0 cellspacing=0 width=100%><TR>\n";
	print "<TH class=\"darkBg rbBorder\" align=\"center\" rowspan=\"$headerRowNb\" nowrap>&nbsp;State&nbsp;</TH>\n";

	if ($hasTotalQty) {
		if ($isConcentration) {
			print "<TH class=\"darkBg rbBorder\" align=\"center\" rowspan=\"$headerRowNb\" nowrap>&nbsp;Total molar&nbsp;<BR>&nbsp;conc. (mmol/mL)&nbsp;</TH>\n";
			print "<TH class=\"darkBg rbBorder\" align=\"center\" rowspan=\"$headerRowNb\" nowrap>&nbsp;Total mass&nbsp;<BR>&nbsp;conc. (µg/mL)&nbsp;</TH>\n";
		} else {
			print "<TH class=\"darkBg rbBorder\" align=\"center\" rowspan=\"$headerRowNb\" nowrap>&nbsp;Total molar&nbsp;<BR>&nbsp;amount (mmol)&nbsp;</TH>\n";
			print "<TH class=\"darkBg rbBorder\" align=\"center\" rowspan=\"$headerRowNb\" nowrap>&nbsp;Total mass&nbsp;<BR>&nbsp;amount (ng)&nbsp;</TH>\n";
		}
	}
	print "<TH class=\"darkBg rbBorder\" align=\"center\" rowspan=\"$headerRowNb\" nowrap>&nbsp;Sample&nbsp;</TH>\n";
	print "<TH class=\"darkBg rbBorder\" align=\"center\" colspan=\"$headerColNb\" nowrap>&nbsp;Provided quantities ($unitsMatching{$initialUnits})&nbsp;</TH>\n";

	if ($isCalibration) {
		print "</TR><TR>\n";
		foreach my $protID (sort {$a <=> $b} keys %proteins) {
			print "<TH class=\"darkBg rbBorder\" align=\"center\" nowrap>&nbsp;$proteins{$protID}&nbsp;</TH>\n";
		}
		print "</TR>\n";
	}

	my $bgColor = "lightBg";
	foreach my $statePos (sort{$a <=> $b} keys %states) {
		my $stateName = $states{$statePos}[0];
		my $repNb = $states{$statePos}[1];
		my $repNum = 0; 
		print "<TR><TH class=\"darkBg rbBorder\" valign=\"center\" align=\"center\" rowspan=\"$repNb\" nowrap>&nbsp;$stateName&nbsp;</TH>\n";
		
		if ($hasTotalQty) {
			print "<TD class=\"lightBg rbBorder\" valign=\"center\" align=\"center\" rowspan=\"$repNb\" nowrap>&nbsp;$totalQty{$statePos}{'MOL'}&nbsp;</TD>\n";
			print "<TD class=\"lightBg rbBorder\" valign=\"center\" align=\"center\" rowspan=\"$repNb\" nowrap>&nbsp;$totalQty{$statePos}{'MASS'}&nbsp;</TD>\n";
		}

		foreach my $bioRep (sort{$a <=> $b} keys %{$initialQty{$statePos}}) {
			foreach my $techRep (sort{$a <=> $b} keys %{$initialQty{$statePos}{$bioRep}}) {
				$repNum ++;
				print "<TR>\n" unless ($bioRep == 1 && $techRep == 1);
				my $label = "State$statePos-bioRep$bioRep-techRep$techRep";
				my $tdClass = ($repNum == $repNb)? "class=\"$bgColor rbBorder\"" : "class=\"$bgColor rBorder\"";
				print "<TD $tdClass nowrap>&nbsp;$labelsName{$label}&nbsp;</TD>\n";
				if ($isCalibration) {
					foreach my $protID (sort{$a <=> $b} keys %proteins) {
						if ($initialQty{$statePos}{$bioRep}{$techRep}{$protID}) {
							print "<TD $tdClass align=\"center\" nowrap>&nbsp;$initialQty{$statePos}{$bioRep}{$techRep}{$protID}&nbsp;</TD>\n";
						} else {
							print "<TD $tdClass align=\"center\" nowrap> - </TD>\n";
						}
					}
				} else {
					print "<TD $tdClass align=\"center\" nowrap>&nbsp;$initialQty{$statePos}{$bioRep}{$techRep}&nbsp;</TD>\n";
				}
				print "</TR>\n";
				$bgColor = ($bgColor eq "lightBg") ? "darkBg" : "lightBg";
			}
		}
	}
	print "</TABLE></CENTER>\n";
	$dbh->disconnect;
}


####>Revision history<#####
# 2.14.4 [BUGFIX] Handles rare MaxQuant bug where a quantified protein is not identified in selected Analyses (PP 06/08/21)
# 2.14.3 [BUGFIX] Removed forgotten debug text & fixed global correlation heatMap (PP 26/06/21)
# 2.14.2 [BUGFIX] Fix minor bug in AJAX protein list table header for MaxQuant (PP 16/06/21)
# 2.14.1 [BUGFIX] Fix error when trying to display raw peptide data for Proteomic Ruler (VL 10/06/21)
# 2.14.0 [CHANGE] Always uses max num peptides used for point size & multiple bug fixes for PEP_RATIO (PP 26/05/21)
# 2.13.7 [BUGFIX] Fix listing of manually selected peptide for bias correction (PP 03/05/21)
# 2.13.6 [BUGFIX] Fix non-implemented peptide view popup for Abundance (PP 15/03/21)
# 2.13.5 [FEATURE] Displays missing values at protein-level for myProMS Abundance quantification (PP 12/03/21) 
# 2.13.4 [FEATURE] Compatible with context for Free residues and PTM-enriched proteins (PP 08/03/21) 
# 2.13.3 [FEATURE] Displays MY_LFQ-specific lost proteins/sites (PP 12/02/21)
# 2.13.2 [Merge] master into ppoullet (PP 02/02/21)
# 2.13.1b [FEATURE] Abundance multi-level bias correction & MG shared peptides rule && contaminants exclusion (PP 01/02/21)
# 2.13.1 [MINOR] Displays TopN parameters for Spectronaut quantifications (VS 29/01/21)
# 2.13.0 [FEATURE] Added "% valid peptides" & moved design-quantif data fetch to a &fetchDesignQuantificationData (PP 04/01/21)
# 2.12.4 [FEATURE] Add export of absolute quantif + Fix export from list view for other measures (VL 16/12/20)
# 2.12.3 [ENHANCEMENT] Add display of 95% confidence intervals for absolute quantification (VL 07/12/20)
# 2.12.2 [ENHANCEMENT] Add display of total and initial protein quantities for absolute quantification (VL 30/11/20)
# 2.12.1 [MINOR] Changes for sample names in bias correction plot to handle PEP_RATIO quantifications (PP 08/12/20)
# 2.12.0 [FEATURE] Displays reference peptides data for TDA & compatible with Free residues quantification (PP 23/11/20)
# 2.11.13 [ENHANCEMENT] Improved extended search (PP 13/10/20)
# 2.11.12b [UPDATE] Parent XIC quantification is displayed in observations popup (PP 13/10/20)
# 2.11.12 [ENHANCEMENT] Add parameters for absolute quantifications (VL 13/10/20)
# 2.11.11 [MINOR] Fix proteomic ruler parent quantification parameter display (VL 19/10/20)
# 2.11.10 [FEATURE] New scatter plot to display bias coefficients (PP 07/10/20)
# 2.11.9 [ENHANCEMENT] Improved handling of negative targetPos (PP 21/09/20)
# 2.11.8b [BUGFIX] Negative replicate targetPos are now retrieved by SQLite query (PP 17/09/20)
# 2.11.8 [ENHANCEMENT] Add displaying of spectronaut/ptmRS score in peptide view (VS 16/09/20) 
# 2.11.7 [BUGFIX] Fix filtering based on num. peptides in one replicate (PP 16/09/20)
# 2.11.6 [BUGFIX] Fix typo in & promsMod subroutine call (PP 11/09/20)
# 2.11.5 [BUGFIX] Removed forgotten reference to ID_QUANTIFICATION field in SQLite query in &ajaxPeptideSwath (PP 02/09/20)
# 2.11.4 [UPDATE] Change in "shared peptides" parameter handling & remove "shared proteins" (PP 17/08/20)
# 2.11.3 [BUGFIX] JS fix for display of multi-file normalization boxplots (PP 16/08/20)
# 2.11.2 [FEATURE] Multi-file boxplot for normalization & compatible with "Shared peptides" normalization (PP 14/08/20)
# 2.11.1 [BUGFIX] Fix display of missing values summary for MSstats (PP 13/08/20)
# 2.11.0 [ENHANCEMENT] Uses quantification data from SQLite file & default init for $quantifType & site sequence context compatible (PP 30/07/20)
# 2.10.4 [ENHANCEMENT] Added site information displaying for protein/peptide intensity graph/table (VS 16/07/2020)
# 2.10.3.1 [FEATURE] Option to display full ratio range & sites are correctly displayed in abundance graphical view (PP 29/06/20)
# 2.10.3 [BUGFIX] Added forgotten check for defined values for truely identified peptides (PP 18/06/20)
# 2.10.2 [FEATURE] Normalization with shared proteins and LFQ 2 parameters and ghost peptide flagging (PP 08/06/20)
# 2.10.1 [MINOR] Changed MaxQuantProbIcon to PTMProbIcon (VS 06/06/20)
# 2.10.0 [FEATURE] Non-ghost peptides filtering & [BUGFIX] Fix issues #97,#100 (PP 20/05/20)
# 2.9.4 [CHANGE] Uses 'sites' instead of 'isoforms' for PTM-quantifications (PP 06/04/20)
# 2.9.3 [BUGFIX] Clean List name before passing to JS in &ajaxListDecorateGraph (PP 02/04/20)
# 2.9.2 [UPDATE] Changed RANK field to IDENT_RANK for compatibility with MySQL 8 (PP 04/03/20) 
# 2.9.1 [BUGFIX] In intensity graph for PTMs & restores display of MSstats bias correction info (PP 03/03/20)
# 2.9.0 [ENHANCEMENT] Optimizations to speed up summary display and data extraction from DB & display p-value +/- adjusted (PP 25/02/20)
# 2.8.4 [BUGFIX] When removing proteins from list, only apply deletion on selected list (VS 17/02/20)
# 2.8.3 [BUGFIX] in &ajaxListSelectedProteins & [ENHANCEMENT] in auto-refresh of on-going quantif (PP 07/02/20)
# 2.8.2 [BUGFIX] Fix minor display bugs for DIA:Abundance (PP 29/01/20)
# 2.8.1 [BUGFIX] When highlighting of a site-quantif with a custom list of proteins (PP 23/01/20)
# 2.8.0 [FEATURE] Handles myProMS Protein Abundance (PP 03/01/20)
# 2.7.8 [BUGFIX] Replicate filtering, list filtering +/- sites, DIA missing value display, selected protein export & [ENHANCEMENT] Peptide-based normalization (PP 19/12/19)
# 2.7.7 [CHANGE] Change "best" state to "one" state for replicate peptide filtering (PP 25/11/19)
# 2.7.6 [FEATURE] Compatible with option to extend missed cleavage exclusion to overlapping peptides & [BUGFIX] in PTM-site management (PP 21/11/19)
# 2.7.5 [BUGFIX] in &ajaxListSelectedProteins in replicate peptide filtering (PP 07/11/19)
# 2.7.4 [ENHANCEMENT] Deletion of non-design quantif uses &promsQuantif::deleteQuantification (PP 06/11/19)
# 2.7.3 [FEATURE] Full support to multi-modifs quantif & TDA (PP 31/10/19)
# 2.7.2 [BUGFIX] Add \%dispModifSites which was missing in arguments of &listProteinRatios and pValueCode was taken instead (VL 28/10/19)
# 2.7.1 [ENHANCEMENT] Adapt code to new intensity metrics in Proteomic Ruler (VL 25/09/19)
# 2.7.0 [FEATURE] Merge with Proteomic Ruler code (PP 29/08/19)
# 2.6.6 Add display of custom proteins for Proteomic Ruler (VL 21/08/19)
# 2.6.5 Add possibility to display Proteomic Ruler quantifications (VL 20/08/19)
# 2.6.4.1 [FEATURE] Handles multi-modif & partial support for TDA by myProMS algo (PP 26/08/19)
# 2.6.4 [FEATURE] Export proteins selected from volcano plot (VS 09/08/19)
# 2.6.3 Switch protein means when necessary in &ajaxPeptideRatios to handle R AnaDiff bug for PEP_INTENSITY (PP 17/04/19)
# 2.6.2 Compatible with recreated sites (PP 04/02/19)
# 2.6.1 Minor change in setting default value for num peptide threshold (PP 25/01/19)
# 2.6.0 Displays MaxQuant probabilities as popup in peptide view for any PTM-quantif (PP 09/01/19)
# 2.5.0 Compatible with protein-level normalization of PTMs & missing values display (PP 19/12/18) 
# 2.4.7 Improvement in emPAI/SIN display (PP 27/11/18)
# 2.4.6 Added TDA displaying (VS 26/11/18)
# 2.4.5 [Fix] Bug in &ajaxRestrictProteinList (PP 11/09/18)
# 2.4.4 Chart highlighting now possible from custom lists (PP 07/09/18)
# 2.4.3 Clean protein alias from badly parsed characters in case identifier conversion (PP 09/07/18)
# 2.4.2 Updated quanti selection in ana call for compatibility with TMT and DIA.
# 2.4.1 [Fix] bug display of negative normalization proteins list & added = in pepID string split (PP 15/06/18)
# 2.4.0 Reads peptide quantif data from file(s) not DB (PP 11/05/18)
# 2.3.8 [Fix] Minor JS bug when displaying ratios distribution (PP 07/05/18)
# 2.3.7 [Fix] Proper image file selection for label-free v2 vs v3 algo (PP 16/04/18)
# 2.3.6 Modif to save a quanti as template (MLP 16/04/18)
# 2.3.5 Minor change is text display (PP 06/04/18)
# 2.3.4 Quantifications with myProMS algo v3 are now exportable in XLS<BR>TODO: Export MaxQuant (PP 16/02/18)
# 2.3.3 Minor update to handle NAs in algo v3 correlation matrix (PP 13/12/17)
# 2.3.2 Fix minor bug in &ajaxPeptideRatios for algo v2 %seenExcluded dataSrc (PP 11/12/17)
# 2.3.1 Uses color-codes for peptide source in peptide view for algo v3 (04/12/17)
# 2.3.0 Compatible with quantif algo myProMS v3 (PP 27/11/17)
# 2.2.3 Modifs to display DIA quantification from MLP (PP 16/11/17)
# 2.2.2 Now uses &promsQuantif::fetchSiteList, update in &ajaxManageSaveProteins function & auto-reload for on-going quantification (PP 11/08/17)
# 2.2.1 Filters sites if restriction list and quantif are sites (PP 07/08/17)
# 2.2.0 Compatible with MODIFICATION_SITE, improved BioRep naming & include GA 2.1.4b update (PP 03/08/17)
# 2.1.4 Updated Excel colors (PP 22/02/17)
# 2.1.3 Minor modification (MLP 26/01/17)
# 2.1.2 Improved handling of MaxQuant quantif (PP 20/01/17)
# 2.1.1 Fill missing format for 2 merge_range functions in XLS export (GA 20/12/16)
# 2.1.0 Major changes in MaxQuant data display<BR>TODO: Handle export (PP 13/12/16)
# 2.0.2 Minor modification for emPAI (GA 25/11/16)
# 2.0.1 Bug fix in SWATH fragment display for multi-ratio quantif (PP 26/10/16)
# 2.0.0 Reverse ratios handling & Coef. variation filtering & SWATH peptide display (PP 19/10/16)
# 1.8.4 Minor fix to keep +/-inf flag in log2 ratios during export (PP 10/08/16)
# 1.8.3 Uses SINGLE_REF flag for multi-ratio quantifs (PP 04/08/16)
# 1.8.2 Handles SWATH & now used by 'summary' call of protein ratio quantifs (PP 02/08/16)
# 1.8.1 Minor modification for MaxQuant (GA 01/08/16)
# 1.8.0 Bug fix for rare cases of missing ratio due to exclusion of all peptides by global CV outlier detection (PP ../05/16)
# 1.7.9 Display MQ quantifications in myProMS (intensity, iBAQ, LFQ, MS/MS count) (GA 09/05/16)
# 1.7.8 Compatibility with peptidePlot.js v1.1.0 & bug fix for raw peptide data display of iTRAQ & other minor updates (PP 28/04/16)
# 1.7.7 Added matching peptides option for label-free quantification (PP 25/02/16)
# 1.7.6 Minor bug fix in &ajaxListSelectedProteins for label-free (PP 19/02/16)
# 1.7.5 Minor change to prevent undefined quantif name (PP 11/02/16)
# 1.7.4 Bug fix in peptide display of internal quantif & +/-inf conversion limit for volcano only (PP 02/02/16)
# 1.7.3 Bug fix in $sthPR->execute for internal quantif (PP 30/01/16)
# 1.7.2 Compatible with quantifSet XIC summing (PP 19/01/16)
# 1.7.1 Compatible with old FC quantifs (PP 15/10/15)
# 1.7.0 Compatible with PTM quantification & called openProject v1.5.9 (PP 13/05/15)
# 1.6.5 display goTerm from child Go Analysis and pathway term from Pathway analysis, add revert2log function (SL 05/11/14)
# 1.6.4 Minor bug correction numStates for internal quantifs (PP 23/09/14)
# 1.6.3 Uses &promsMod::printAjaxManageSaveProteins function for List save management (PP 25/08/14)
# 1.6.2 Exploratory Analysis constraint on List deletion & fix in R image display (PP 25/08/14)
# 1.6.1 Update for EMPAI in export and printing (GA 26/06/14)
# 1.6.0 Support for SuperRatio, default view is now List & change "Category" to "List" (PP 19/06/14)
# 1.5.1 Default display is now 3 peptide/protein (PP 29/04/14)
# 1.5.0 More links to ratios statistics (PP 14/04/14)
# 1.4.9 Add log2 option to compare XIC extractions when call=analysis<BR>TODO: change color definition (GA 24/03/14)
# 1.4.8 Uses rmtree instead of remove_tree (PP 10/03/14)
# 1.4.7 ratio display management for internal quantif & graph view merge for design/internal quantif (PP 27/02/13)
# 1.4.6 Export to Excel (PP 20/02/14)
# 1.4.5 Export SVG & improvement of ajaxPepRatio action (PP 11/02/14)
# 1.4.4 system command removal (PP 12/11/13)
# 1.4.3 Displays minimum infinite ratios choice (PP 21/10/13)
# 1.4.2 Improved handling of infinite ratios (PP 09/10/13)
# 1.4.1 Minor modification to avoid to skip $labelingInfo information (GA 23/09/13)
# 1.4.0 Generated from listAnaQuantifications.cgi 1.3.3PP2 (PP 02/09/13)
# 1.3.3PP2 Volcano plot with missing p-values set to 1 & infinite ratios detection (PP 21/08/13)
# 1.3.3 Minor modification in XLS formating for XIC-export (GA 20/08/13)
# 1.3.2 Add export option for XIC-quantification results (GA 14/08/13)
# 1.3.1 Another improvement of new design display (PP ../07/13)
# 1.3.0 Improved display of new design options (PP 11/07/13)
# 1.2.9 Minor change (PP 08/07/13)
# 1.2.8 Fix missing parameters for ms2 alignment during MCQ extraction (GA 05/07/13)
# 1.2.7 extended design management (PP 04/07/13)
# 1.2.6 Better label detection & fix bug ghost peptide modification (PP 04/07/13)
# 1.2.5 Remove VAR_MOD from script (GA 28/05/13)
# 1.2.4 Modifs for unified naming of files between label and unlabel quantifs (PP 16/04/13)
# 1.2.3 Restrict GO decoration selection to qualitative GO analyses (FY 04/04/13)
# 1.2.2 ABS(PEP_BEG) on Label-quantif because of negative values for ghost-peptides (PP 03/04/13)
# 1.2.1 Minor modification in &deleteQuantification -> delete virtual proteins (GA 08/03/13)
# 1.2.0 Minor modification -> ABS(PEP_BEG) because of negative values for ghost-peptides (GA 08/02/13)
# 1.1.9 Ratio selection & compatibility with new peptide filters (PP 15/01/13)
# 1.1.8 Minor display bugs correction (PP 09/01/13)
# 1.1.7 Minor adding column Mr(obs) for XIC printing (GA 07/12/12)
# 1.1.6 Correction of bug in ajaxListSelectedProteins that was not working properly + print XIC (GA 23/11/12)
# 1.1.5 Aesthetic modifications for 'Peptides used' column (GA 16/11/12)
# 1.1.4 Add quantification parameters for XIC extractions + add deletion of ghost peptides associated to a quantification (GA 22/10/12)
# 1.1.3 Compatible with new data path .../results/project_ID/GoID (PP 13/09/12)
# 1.1.2 Minor bug fixes & improvements (PP 30/07/12)
# 1.1.1 Added custom protein list filtering on display (PP 18/07/12)
# 1.1.0 Find all GO analyses linked to Experiment & true number of filtered proteins (PP 08/06/12)
# 1.0.9 Added save proteins to List. Still not fully stable (PP 07/06/12)
# 1.0.8 More code factorisation & minor display changes: Not fully stable (PP 23/04/12)
# 1.0.7 Bug correction & code factorization (PP 19/04/12)
# 1.0.6 Merge 1.0.5GA & 1.0.4c (PP 06/04/12)
# 1.0.4c PEPTIDE.QUANTI_DATA change to DATA (PP 05/04/12)
# 1.0.4b minor change in peptide selection display (PP 27/03/12)
# 1.0.4 added graphical display of peptide for labeled quanti (PP 27/02/12)
# 1.0.3 Added selected sources (merged MSF) & graphical displays (PP 02/02/12)
# 1.0.2 Protein-level internal analysis quantification (PP 19/09/11)
# 1.0.1 No longer uses CHANNEL as quantif parameter (PP 29/07/11)
# 1.0.0 Created (PP 25/07/11)
