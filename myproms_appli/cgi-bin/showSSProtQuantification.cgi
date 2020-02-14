#!/usr/local/bin/perl -w

################################################################################
# showSSProtQuantification.cgi     1.0.5                                       #
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
use strict;
use Spreadsheet::WriteExcel;
use utf8; # Tells Perl that characters are UTF8. Necessary for Excel export to work with UTF-8 characters ...
use Encode qw(encode_utf8); # ... Encode needed if hard UTF8 chars in code!!!

#print header(-charset=>'utf-8'); warningsToBrowser(1); #DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my $projectAccess; # must be global
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my %labelingName=('FREE'=>'Label-free','SILAC'=>'SILAC','ITRAQ'=>'iTRAQ');
my %normalizationNames=&promsQuantif::getQuantifNormalizationName;
my %featureCodes=(sp_count=>['Peptide/spectrum matches','PSMs'],xic=>['Peptide abundance','XICs'],
				  all_ion=>['All peptide ions','ions'],all_pep=>['All peptides','peptides'],
				  dist_ion=>['Distinct peptide ions','dist. ions'],dist_pep=>['Distinct peptides','dist. pep.'],dist_seq=>['Distinct peptide sequences','dist. seq.']
				  );

####################
####>Parameters<####
####################
my $action=(param('ACT'))? param('ACT') : 'select';
if ($action eq 'exportProtData') {&exportProteinData; exit;}
my $selQuantifID=(param('id_quantif'))? param('id_quantif') : 0; # 0=> call=ana
my $view=param('view') || 'graph'; # list
my $dispPvalue=param('pValue') // 0.05; $dispPvalue=1 if ($dispPvalue <= 0 || $dispPvalue > 1);
my $dispDelta=param('delta') // 50; $dispDelta=0 if ($dispDelta < 0 || $dispDelta > 100);
my $dispMeanPep=param('meanPep') || 0;
my $dispSort=param('sort') || 'delta'; # for list only
my %sortOptions=('set'=>['State specificity',1],'delta'=>['Best delta',2],'p-value'=>['Adj. p-value',3],'mean'=>['Mean feature count',4],'identifier'=>['Identifiers',5],'mw'=>['Molecular weight',6]);
my (%dispSets,%condToState); # globals also for some ajax calls
if (param('dispSets')) {foreach my $setPos (param('dispSets')) {@{$dispSets{$setPos}}=();}} # Positions of subsets displayed in results
my $restrictListID=param('restrictList') || 0;

if ($action eq 'ajaxListProt') {&ajaxListSelectedProteins; exit;}
elsif ($action eq 'ajaxDrawProt') {&ajaxDrawSelectedProteins; exit;}
elsif ($action eq 'ajaxPepDistrib') {&ajaxPeptideDistribution; exit;}

my $dispResults=param('displayRes') || 0; # display form was submitted

####>Connect to the database
my $dbh=&promsConfig::dbConnect;
my $projectID=&promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION');
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
$projectAccess=${$userInfo[2]}{$projectID};

################
####>Main 1<####
################
my (%anaQuantifList,%quantificationTypes,@quantificationInfo,%methodList);
#my $highlighMatchStrg='';
my ($workbook,%itemFormat); # globals for export to Excel

#my ($designID,$quantifModID,$modifName)=$dbh->selectrow_array("SELECT ID_DESIGN,M.ID_MODIFICATION,PSI_MS_NAME FROM QUANTIFICATION Q
#																	LEFT JOIN MODIFICATION M ON Q.ID_MODIFICATION=M.ID_MODIFICATION
#																	WHERE ID_QUANTIFICATION=$selQuantifID");
my ($designID)=$dbh->selectrow_array("SELECT ID_DESIGN FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");

#my $title='State-Specific ';
#if ($quantifModID) {
#	$title.=$modifName.'-';
#	$highlighMatchStrg=",'^###-'"; # for matching modProtID with protID in graph
#}
#$title.='Protein Analysis';
my $title='State-Specific Protein Analysis';


if ($action eq 'export') {
	#################################
	####>Prepare export to Excel<####
	#################################
	my $timeStamp1=strftime("%Y%m%d %H-%M",localtime);
	#my $timeStamp2=strftime("%Y-%m-%d %H:%M",localtime);

	$workbook=Spreadsheet::WriteExcel->new("-");
	$workbook->set_properties(title=>"Protein quantification data",
							  author=>'myProMS server',
							  comments=>'Automatically generated with Perl and Spreadsheet::WriteExcel'
							  );
	#$workbook->set_custom_color(40,224,224,255); # light light color #E0E0FF
	#$workbook->set_custom_color(41,176,176,255); # dark blue color  #B0B0FF
	$workbook->set_custom_color(40,189,215,255); # light light color  #BDD7FF (new myProMS colors V3.5+)
	%itemFormat=(
			title =>			$workbook->add_format(align=>'center',size=>18,bold=>1,border=>1),
			header =>			$workbook->add_format(align=>'center',size=>10,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			headerR =>			$workbook->add_format(align=>'right',valign=>'top',size=>10,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
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
<TITLE>List of Quantifications</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.TD {font-weight:normal;}
.TH {font-weight:bold;}
.LINK {cursor:pointer;}
.highlight{width:250px;}
.popup {z-index:100;background-color:#FFFFFF;border:solid 3px #999999;padding:5px;position:absolute;display:none;box-shadow:10px 10px 20px #808080;}
</STYLE>
|;
	if ($action eq 'summary') {
		print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
		&promsMod::popupInfo;
	}
	else {
		print qq
|<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/genericPlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/boxPlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/peptidePlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
<SCRIPT LANGUAGE="JavaScript">
|;
	&promsMod::popupInfo;
	print qq
|var view='$view';
function sequenceView(id_protein,anaIdStrg){
	var winLocation="$promsPath{cgi}/sequenceView.cgi?id_prot="+id_protein+"&msdata="+top.showMSdata+"&id_ana="+anaIdStrg;
	top.openProtWindow(winLocation);
}
|;
	}
	print qq
|</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">$title</FONT>
<BR>
|;
}

#print "*** ACTION=$action ***<BR>\n";


###################################################
####>Fetching data for selected Quantification<####
###################################################
my ($quantifAnnot,$quantifStatus,$quantifMethodID,$updateDate,$updateUser)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,STATUS,ID_QUANTIFICATION_METHOD,UPDATE_DATE,UPDATE_USER FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
my ($labelType)=($labelStrg)? ($labelStrg=~/LABEL=(.+)/) : ('FREE');
$labelType='FREE' unless $labelType;
if ($action=~/summary|export/) {
	if ($updateUser) {
		my $sthUser=$dbh->prepare("SELECT USER_NAME FROM USER_LIST WHERE ID_USER=?");
		$sthUser->execute($updateUser);
		my ($userName)=$sthUser->fetchrow_array;
		$updateUser=$userName || $updateUser;
	}
	else {$updateUser='?';}
}

##################################
####>Protein restriction list<####
##################################
my %restrictProteins;
if ($restrictListID) {
	my $sthRP=$dbh->prepare("SELECT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$restrictListID");
	$sthRP->execute;
	while (my($protID)=$sthRP->fetchrow_array) {
		$restrictProteins{$protID}=1;
	}
	$sthRP->finish;
}

######################
####>SSPA QUANTIF<####
######################
my $minPvalue=1;

my ($usedAnaID)=$dbh->selectrow_array("SELECT A.ID_ANALYSIS FROM ANALYSIS A,ANA_QUANTIFICATION AQ WHERE AQ.ID_ANALYSIS=A.ID_ANALYSIS AND ID_QUANTIFICATION IN (SELECT ID_PARENT_QUANTIFICATION FROM PARENT_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID) ORDER BY NAME ASC LIMIT 0,1");# Select one of the ANALYSIS related to that TNPQ quantification so as to make sequenceView method work...

my (%labelingInfo,%quantifParamInfo,%quantifParamID2Code);
my $sthQP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,NAME,CODE FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=$quantifMethodID");
$sthQP->execute;
while (my ($paramID,$paramName,$paramCode)=$sthQP->fetchrow_array) {
	@{$quantifParamInfo{$paramCode}}=($paramID,$paramName);
	$quantifParamID2Code{$paramID}=$paramCode;
}
$sthQP->finish;

foreach my $infoStrg (@labelInfo) {
	my ($setting,$valueStrg)=split('=',$infoStrg);
	if ($setting eq 'PEPTIDES') {$valueStrg=~s/\#//g;} # remove ID tag of selected modifications if any
	@{$labelingInfo{$setting}}=split(';',$valueStrg);
}
my $numStates=scalar @{$labelingInfo{'STATES'}};

if ($action ne 'summary' && !scalar keys %dispSets) { # set default dispSets selection to single set(s)
	foreach my $setPos (1..$numStates) {
		#last if $labelingInfo{'SETS'}[$setPos-1]=~/\+/; # in case some single-sets are not used (used < numStates)
		@{$dispSets{$setPos}}=();
	}
}

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

#### ---> Form here on action cannot be 'summary' --->
print "<DIV id=\"waitDiv\"><BR><BR><BR><BR><BR><FONT class=\"title3\">Fetching data. Please wait...</FONT><BR><IMG src=\"$promsPath{images}/scrollbarGreen.gif\"></DIV>\n" if $action ne 'export';

###>Fetching protein quantification data
my (%quantifValues,%proteinInfo);
if (scalar keys %restrictProteins) {
	foreach my $protID (keys %restrictProteins) {@{$proteinInfo{$protID}}=();} # in case prot in restrict are not quantified (=>correct num of proteins displayed)
}

##>DELTA_PC & PVAL_ADJ
#my $sthProtQ1=$dbh->prepare("SELECT PQ.ID_PROTEIN,GROUP_CONCAT(RESIDUE,POSITION ORDER BY POSITION SEPARATOR '.'),ID_QUANTIF_PARAMETER,TARGET_POS,QUANTIF_VALUE
#								FROM PROTEIN_QUANTIFICATION PQ
#								LEFT JOIN PROTQUANTIF_MODRES PQMR ON PQ.ID_PROT_QUANTIF=PQMR.ID_PROT_QUANTIF
#								LEFT JOIN MODIFIED_RESIDUE MR ON PQMR.ID_MODIF_RES=MR.ID_MODIF_RES AND MR.ID_QUANTIFICATION=PQ.ID_QUANTIFICATION
#								WHERE PQ.ID_QUANTIFICATION=$selQuantifID AND ID_QUANTIF_PARAMETER IN ($quantifParamInfo{DELTA_PC}[0],$quantifParamInfo{PVAL_ADJ}[0]) GROUP BY PQ.ID_PROT_QUANTIF");
my $sthProtQ1=$dbh->prepare("SELECT PQ.ID_PROTEIN,ID_QUANTIF_PARAMETER,TARGET_POS,QUANTIF_VALUE
								FROM PROTEIN_QUANTIFICATION PQ
								WHERE PQ.ID_QUANTIFICATION=$selQuantifID AND ID_QUANTIF_PARAMETER IN ($quantifParamInfo{DELTA_PC}[0],$quantifParamInfo{PVAL_ADJ}[0]) GROUP BY PQ.ID_PROT_QUANTIF");
$sthProtQ1->execute;
my (%bestSet,%decodedMods);
while (my ($protID,$paramID,$setPos,$qValue)=$sthProtQ1->fetchrow_array) { #,$modResStrg
	next if ($restrictListID && !$restrictProteins{$protID}); # skip if list restriction
	next unless $dispSets{$setPos}; # skip if setPos not selected
	#@{$proteinInfo{$protID}}=(); # actual number of quantified proteins (only list filtering)
	#if ($modResStrg) {
	#	unless (defined $decodedMods{$modResStrg}) {
	#		$decodedMods{$modResStrg}=($modResStrg)? '-'.&promsQuantif::decodeModifPosAmbiguity($modResStrg) : ''; # quantif of modification
	#	}
	#	$protID.=$decodedMods{$modResStrg};
	#}
	#else {$modResStrg='';}
	if ($view eq 'graph') {$quantifValues{$setPos}{$protID}{$quantifParamID2Code{$paramID}}=$qValue;}
	else {
		$quantifValues{$protID}{$quantifParamID2Code{$paramID}}=$qValue;
	}
	$bestSet{$protID}=$setPos;
	$minPvalue=$qValue if ($quantifParamID2Code{$paramID} eq 'PVAL_ADJ' && $qValue > 0 && $qValue < $minPvalue);
}
$minPvalue/=10;

$sthProtQ1->finish;
#print "<BR>--DATA#1 OK!<BR>\n";

##>Mean
#my $sthProtQ2=$dbh->prepare("SELECT PQ.ID_PROTEIN,GROUP_CONCAT(RESIDUE,POSITION ORDER BY POSITION SEPARATOR '.'),TARGET_POS,QUANTIF_VALUE
#								FROM PROTEIN_QUANTIFICATION PQ
#								LEFT JOIN PROTQUANTIF_MODRES PQMR ON PQ.ID_PROT_QUANTIF=PQMR.ID_PROT_QUANTIF
#								LEFT JOIN MODIFIED_RESIDUE MR ON PQMR.ID_MODIF_RES=MR.ID_MODIF_RES AND MR.ID_QUANTIFICATION=PQ.ID_QUANTIFICATION
#								WHERE PQ.ID_QUANTIFICATION=$selQuantifID AND ID_QUANTIF_PARAMETER=$quantifParamInfo{MEAN}[0] GROUP BY PQ.ID_PROT_QUANTIF");
my $sthProtQ2=$dbh->prepare("SELECT PQ.ID_PROTEIN,TARGET_POS,QUANTIF_VALUE
								FROM PROTEIN_QUANTIFICATION PQ
								WHERE PQ.ID_QUANTIFICATION=$selQuantifID AND ID_QUANTIF_PARAMETER=$quantifParamInfo{MEAN}[0] GROUP BY PQ.ID_PROT_QUANTIF");
$sthProtQ2->execute;
while (my ($protID,$statePos,$qValue)=$sthProtQ2->fetchrow_array) { #,$modResStrg statePos 1 to numStates are retrieved
	#if ($modResStrg) {
	#	next unless defined $decodedMods{$modResStrg}; # mod never seen => prot never seen
	#	$protID.=$decodedMods{$modResStrg};
	#}
	#else {$modResStrg='';}
	next unless $bestSet{$protID}; # skip if protein not is selected sets
	if ($view eq 'graph') {$quantifValues{$statePos}{$protID}{MEAN}=$qValue;} # WARNING: Not the same structure for MEAN
	else {$quantifValues{$protID}{MEAN}{$statePos}=$qValue;}
}
$sthProtQ2->finish;

#print "<BR>--DATA#2 OK!<BR>\n";

###>Compute test State mean & filter on selected thresholds
if ($view eq 'graph') {
	foreach my $protID (keys %bestSet) {
		#my $setPos=$bestSet{$protID};
		#next unless defined $quantifValues{$setPos}{$protID}{DELTA_PC};
		my $minSetMean;
		#foreach my $statePos (1..$numStates) { #} find smallest mean in set
		foreach my $statePos (split(/\+/,$labelingInfo{SETS}[$bestSet{$protID}-1])) { # state(s) involved in set => find smallest mean in best set
			$minSetMean=$quantifValues{$statePos}{$protID}{MEAN} if (!defined $minSetMean || $minSetMean > $quantifValues{$statePos}{$protID}{MEAN});
		}
		if ($minSetMean < $dispMeanPep) {
			foreach my $statePos (1..$numStates) { # MEAN data
				delete $quantifValues{$statePos}{$protID};
			}
			delete $quantifValues{$bestSet{$protID}}{$protID}; # OTHER data
			delete $bestSet{$protID};
			next;
		}
		$quantifValues{$bestSet{$protID}}{$protID}{SET_MEAN}=$minSetMean;
		my ($protID,$modResStrg)=($protID=~/^(\d+)(.*)/);
		@{$proteinInfo{$protID}}=(); # actual number of quantified proteins (only list filtering)
	}
}
else { # list
	foreach my $protID (keys %bestSet) {
		#next if !defined $quantifValues{$protID}{DELTA_PC};
		if ($quantifValues{$protID}{DELTA_PC} < $dispDelta || $quantifValues{$protID}{PVAL_ADJ} > $dispPvalue) {
			delete $quantifValues{$protID};
			delete $bestSet{$protID};
			next;
		}
		my $minSetMean;
		foreach my $statePos (split(/\+/,$labelingInfo{SETS}[$bestSet{$protID}-1])) { # state(s) involved in set => find smallest mean in best set
			$minSetMean=$quantifValues{$protID}{MEAN}{$statePos} if (!defined $minSetMean || $minSetMean > $quantifValues{$protID}{MEAN}{$statePos});
		}
		if ($minSetMean < $dispMeanPep) {
			delete $quantifValues{$protID};
			delete $bestSet{$protID};
			next;
		}
		$quantifValues{$protID}{SET_MEAN}=$minSetMean;
		$quantifValues{$protID}{SET}=$bestSet{$protID};
		my ($protID,$modResStrg)=($protID=~/^(\d+)(.*)/);
		@{$proteinInfo{$protID}}=(); # actual number of quantified proteins (only list filtering)
	}
}

###>Fetching protein info
my $sthPI=$dbh->prepare("SELECT DISTINCT P.ID_PROTEIN,ALIAS,MW,PROT_DES,ORGANISM,PROT_LENGTH FROM PROTEIN_QUANTIFICATION Q,PROTEIN P WHERE Q.ID_PROTEIN=P.ID_PROTEIN AND ID_QUANTIFICATION=$selQuantifID");
my ($GNidentID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
my $sthGN=$dbh->prepare("SELECT MI.VALUE FROM PROTEIN P,MASTERPROT_IDENTIFIER MI WHERE P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND MI.ID_IDENTIFIER=$GNidentID AND P.ID_PROTEIN=? ORDER BY MI.RANK");
my $sthPA=$dbh->prepare("SELECT AP.ID_ANALYSIS FROM ANALYSIS_PROTEIN AP,ANA_QUANTIFICATION AQ WHERE AP.ID_ANALYSIS=AQ.ID_ANALYSIS AND ID_PROTEIN=? AND AQ.ID_QUANTIFICATION=$selQuantifID ORDER BY VISIBILITY DESC,NUM_PEP DESC LIMIT 0,1");# To be sure that the PROTEIN is in the ANALYSIS

$sthPI->execute;
while (my ($protID,@info)=$sthPI->fetchrow_array) {
	next unless $proteinInfo{$protID};
	#next if ($restrictListID && !$restrictProteins{$protID});
	my @geneList=();
	if ($view eq 'list' || $action eq 'export') {
		$info[1]=sprintf "%.1f",$info[1]/1000;
		$sthGN->execute($protID);
		while (my ($gene)=$sthGN->fetchrow_array) {push @geneList,$gene;}
	}
	$sthPA->execute($protID);
	my ($anaID)=$sthPA->fetchrow_array;
	@{$proteinInfo{$protID}}=(@info,\@geneList,$anaID); # [5] for genelist
}
$sthPI->finish;
$sthGN->finish;
$sthPA->finish;

#print "<BR>--PROT OK!<BR>\n";

$dbh->disconnect;

####################
####>Export XLS<####
####################
if ($action eq 'export') {
	&exportProteinList(\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,'PVAL_ADJ',scalar keys %bestSet);
	$workbook->close();
	exit;
}
#print qq |<DIV id="waitDiv"><BR><BR><BR><BR><BR><FONT class="title3">Fetching data. Please wait...</FONT><BR><IMG src="$promsPath{images}/scrollbarGreen.gif"></DIV><SCRIPT type="text/javascript">document.getElementById('waitDiv').style.display='none'</SCRIPT>|;
print qq |<SCRIPT type="text/javascript">document.getElementById('waitDiv').style.display='none';</SCRIPT>|;
#my $numRatios=scalar @{$labelingInfo{RATIOS}};

########################
####>Graphical view<####
########################
if ($view eq 'graph') {
	my $allowHlight=(scalar keys %goAnalyses)? 'true' : 'false';
	&displaySSPAPlot(\%quantifValues,\%proteinInfo,$allowHlight,'PVAL_ADJ',scalar keys %bestSet,"$featureCodes{$labelingInfo{PEPTIDES}[3]}[1]");
}

###################
####>List view<####
###################
else { # View type is list
	print qq
|<SCRIPT type="text/javascript">
window.onload=function() {
	ajaxUpdateRestrict();
}
</SCRIPT>
|;
	&printProteinList(\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,'PVAL_ADJ',scalar keys %bestSet); #$labelingInfo{FDR_ALPHA}[0]/100,$dispStdDev
}
&endHTML;


#############################################<<<SUBROUTINES>>>###########################################

sub endHTML {
	print qq
|</CENTER>
<BR><BR>
<DIV id="displayDIV" class="popup"> <!--filter:alpha(opacity=80);opacity:0.8;-->
	<DIV id="infoDIV"></DIV>
</DIV>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT LANGUAGE="javascript">
setPopup();
</SCRIPT>
</BODY>
</HTML>
|;
}

sub printQuantificationSummary { # GLOBALS: $selQuantifID, $view, $projectID, $designID, $ratioType, $lightColor, $darkColor
	my ($dbh,$anaID,$refInfo,$refGoAnalyses)=@_;
	#my ($lightColor,$darkColor)=&promsConfig::getRowColors;
	my $resultDir="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/results";

	my ($newLine,$startB,$endB)=($action eq 'export')? ("\n",'','') : ('<BR>','<B>','</B>');

	#<Conditions & replicates
	my (%stateInfo,%replicateName,%numTechReplicates,%bioRepLabel,%anaName,%anaInStates,%numCondInSample);
	my $numStates=0;
	my $replicPos=0;
	###< Get parent quantifType to distinguish from label/label-free quantification
	my $sthAnaName=$dbh->prepare("SELECT NAME FROM ANALYSIS WHERE ID_ANALYSIS=?");
	my $sthExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
	my $sthPepQuantif=$dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
	my $sthObsBioSamp=$dbh->prepare("SELECT ID_BIOSAMPLE FROM OBSERVATION WHERE ID_OBSERVATION=?");
	my $sthObsMsSamp=$dbh->prepare("SELECT ID_SAMPLE FROM ANALYSIS A WHERE ID_ANALYSIS=?");
	my $sthCondMsSamp=$dbh->prepare("SELECT COUNT(E.ID_EXPCONDITION) FROM ANALYSIS A,OBSERVATION O,OBS_EXPCONDITION OE,EXPCONDITION E WHERE A.ID_ANALYSIS=O.ID_ANALYSIS AND O.ID_OBSERVATION=OE.ID_OBSERVATION AND OE.ID_EXPCONDITION=E.ID_EXPCONDITION AND E.ID_DESIGN=$designID AND A.ID_SAMPLE=?");
	my $sthBsName=$dbh->prepare("SELECT NAME FROM BIOSAMPLE WHERE ID_BIOSAMPLE=?");
	my $sthMsName=$dbh->prepare("SELECT NAME FROM SAMPLE WHERE ID_SAMPLE=?");
	foreach my $stateData (@{$refInfo->{STATES}}) {
		$stateData=~s/#//g; # removes id tags
		$numStates++;
		my %anaList;
		my $anaCount=0;
		my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$stateData);
		$sthExpCondName->execute($expCondID);
		($stateInfo{$numStates}{'NAME'})=$sthExpCondName->fetchrow_array;
		$condToState{$expCondID}=$numStates;
		$stateInfo{$numStates}{'POPUP'}='';
		my (%replicateInfo,%replicateStrg,%channelName,%bioRepHeader);
		my ($currBioRep,$allTechRepCount,$numPools,$numObs)=(0,0,0,0);
		foreach my $bioReplicate (split(/\./,$quantiObsIDs)) { # bio rep separator
			$currBioRep++;
			my (%fracAnaNames,%numPoolObs,%bioRepNames); # %numObs
			my $numTechRep=0;
			my $numObsInBioRep=0;
			foreach my $techReplicate (split(/&/,$bioReplicate)) { # tech rep separator
				$numTechRep++;
				#my @fracAnaNames;
				$numPoolObs{$numTechRep}=0;
				foreach my $pooledObs (split(/\+/,$techReplicate)) {
					$numPoolObs{$numTechRep}++;
					#my ($parQuantifID,$anaID,$targetPos);
					#my @obsData=split(/:/,$pooledObs);
					#if (scalar @obsData==4) {($parQuantifID,$anaID,$targetPos)=@obsData[1..3];} # design post 07/2013 -> obsID,parQuantifID,anaID,targetPos
					#else { # older design -> parQuantifID,anaID
					#	($parQuantifID,$anaID)=@obsData;
					#	$targetPos=0;
					#}
					my ($obsID,$parQuantifID,$anaID,$targetPos)=split(/:/,$pooledObs);
					unless ($anaName{$anaID}) {
						$sthAnaName->execute($anaID);
						($anaName{$anaID})=$sthAnaName->fetchrow_array;
					}
					$anaList{$anaID}=++$anaCount;
					$sthObsBioSamp->execute($obsID);
					my ($bsID)=$sthObsBioSamp->fetchrow_array;
					$bioRepNames{BIOSAMPLE}{$bsID}++ if $bsID;
					$sthObsMsSamp->execute($anaID);
					my ($sampID)=$sthObsMsSamp->fetchrow_array;
					$bioRepNames{MSSAMPLE}{$sampID}++;
					my $targetName='';
					if ($targetPos) {
						unless ($channelName{$parQuantifID}) {
							$sthPepQuantif->execute($parQuantifID);
							my ($qAnnot)=$sthPepQuantif->fetchrow_array;
							my @quantifInfo=split(/::/,$qAnnot);
							foreach my $channelInfo (@quantifInfo[1..$#quantifInfo]) {
								my ($chanPos,$chanName)=split(/;/,$channelInfo); # + other unsued elements
								$channelName{$parQuantifID}{$chanPos}=$chanName;
							}
						}
						$targetName=' > '.$channelName{$parQuantifID}{$targetPos};
					}
					push @{$fracAnaNames{$numTechRep}},$anaName{$anaID}.$targetName;
				}
				$numObsInBioRep+=$numPoolObs{$numTechRep};
			}
			$allTechRepCount+=$numTechRep;

			$replicateName{++$replicPos}=($action eq 'export')? $stateInfo{$numStates}{'NAME'} : $startB.$stateInfo{$numStates}{'NAME'}.$endB;
			$numTechReplicates{$replicPos}=$numTechRep;

			##<Guess most suitable name for bioRep>##
			##<Bio-sample
			if ($bioRepNames{BIOSAMPLE} && scalar keys %{$bioRepNames{BIOSAMPLE}}==1) {
				my $bsID=(keys %{$bioRepNames{BIOSAMPLE}})[0];
				#if ($bioRepNames{BIOSAMPLE}{$bsID}==$numObsInBioRep) {
					$sthBsName->execute($bsID);
					($bioRepLabel{"$numStates:$currBioRep"})=$sthBsName->fetchrow_array;
				#}
			}
			##<MS Sample or MS Analysis
			if (!$bioRepLabel{"$numStates:$currBioRep"} && scalar keys %{$bioRepNames{MSSAMPLE}}==1) { # BioRep linked to only 1 MS Sample
				my $sampID=(keys %{$bioRepNames{MSSAMPLE}})[0];
				#>Check if MS Sample is linked to only 1 condition in current design
				unless ($numCondInSample{$sampID}) {
					$sthCondMsSamp->execute($sampID);
					($numCondInSample{$sampID})=$sthCondMsSamp->fetchrow_array;
				}
				if ($numCondInSample{$sampID}==1) { # MS Sample name can be used as proxy for bioRep
				#if ($bioRepNames{MSSAMPLE}{$sampID}==$numObsInBioRep) {
					$sthMsName->execute($sampID);
					($bioRepLabel{"$numStates:$currBioRep"})=$sthMsName->fetchrow_array;
				#}
				}
				##<MS Analysis +/- label channel
				elsif ($numObsInBioRep==1) {
					$bioRepLabel{"$numStates:$currBioRep"}=$fracAnaNames{1}[0];
				}
			}
			##<Replicate Rank
			if (!$bioRepLabel{"$numStates:$currBioRep"}) {
				$bioRepLabel{"$numStates:$currBioRep"}='Replicate #'.$currBioRep;
			}

			if ($numObsInBioRep > 1) { # 2 or more fractions and/or tech rep
				if ($numBioRep > 1) {
					$replicateName{$replicPos}.="$startB:$endB ".$bioRepLabel{"$numStates:$currBioRep"};
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
														"  +".join("$newLine  +",@{$fracAnaNames{$techRep}}).$newLine
														: "&nbsp;&nbsp;+".join("$newLine&nbsp;&nbsp;+",@{$fracAnaNames{$techRep}})."<BR>";
					}
					else {
						$replicateName{$replicPos}.="$startB:$endB $fracAnaNames{$techRep}[0]";
						$stateInfo{$numStates}{'POPUP'}.=($action eq 'export')? "•$fracAnaNames{$techRep}[0]$newLine" : "<TD valign=top nowrap>&bull;$fracAnaNames{$techRep}[0]<BR>";
					}
				}
			}
			else { # only 1 obs in replicate
				$replicateName{$replicPos}.="$startB:$endB $fracAnaNames{1}[0]";
				$stateInfo{$numStates}{'POPUP'}.=($action eq 'export')? "•$fracAnaNames{1}[0]$newLine" : "<TD valign=top nowrap>&bull;$fracAnaNames{1}[0]<BR>";
			}
			$numObs+=$numObsInBioRep;
			#$stateString.=" $replicateName{$replicPos},";
		} # end of bioReplicate
		#chop($stateString);
		#$stateInfo{$numStates}{'NAME_LONG'}=$stateString;
		$stateInfo{$numStates}{'NAME_LONG'}=($action eq 'export')? ': ' : ':&nbsp;';
		if ($numBioRep>1) {
			#$existReplicates=1;
			$stateInfo{$numStates}{'NAME_LONG'}.="$numBioRep bio. rep. / ";
		}
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
		$anaInStates{$numStates}=join(',',sort{$anaList{$a}<=>$anaList{$b}} keys %anaList);
	} # end of state
	$sthAnaName->finish;
	$sthExpCondName->finish;
	$sthPepQuantif->finish;
	$sthObsBioSamp->finish;
	$sthObsMsSamp->finish;
	$sthBsName->finish;
	$sthMsName->finish;
	$sthCondMsSamp->finish;

	#<Peptides (retro-compatibility in case no miss-cut param)
	my $pepRangeStrg=($refInfo->{PEPTIDES}[0] eq 'all')? 'Non-proteotypic peptides allowed' : ($refInfo->{PEPTIDES}[0] eq 'unique_shared')? 'Shared peptides allowed if exist proteotypic' : 'Proteotypic peptides only';
	my $pepMissCutStrg=($refInfo->{PEPTIDES}[1]==0)? 'Missed cleavage not allowed' : ($refInfo->{PEPTIDES}[1]==-1)? 'Miss-cleaved & cleaved counterparts not allowed' : 'Missed cleavage allowed';
	my $featureStrg=$featureCodes{$refInfo->{PEPTIDES}[3]}[0] || 'Unknown';
	my $rescuedStrg=($refInfo->{PEPTIDES}[4] || $refInfo->{PEPTIDES}[3] eq 'sp_count')? 'MBR-rescued not allowed' : 'MBR-rescued allowed'; # spectral count=0 for MBR peptides

	my $pepPtmStrg; #=($refInfo->{PEPTIDES}[$pepIdx]==1)? 'PTMs allowed' : ($refInfo->{PEPTIDES}[$pepIdx]==-1)? 'Exclude sequence if PTM found' : 'PTMs not allowed';
	my ($ptmAllowed,@selectedPTMs)=split(/[:,]/,$refInfo->{PEPTIDES}[2]);
	if (abs($ptmAllowed)<=1) { # old options
		$pepPtmStrg=($ptmAllowed==1)? 'All modifications allowed' : ($ptmAllowed==-1)? 'Modified and unmodified matching peptides not allowed' : 'No modifications allowed';
	}
	else { # new custom selection options (2 or -2)
		my @ptmNames;
		my $sthVM=$dbh->prepare('SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DES FROM MODIFICATION WHERE ID_MODIFICATION=?');
		foreach my $modID (@selectedPTMs) {
			$sthVM->execute($modID);
			my ($psiName,$interName,$synName,$des)=$sthVM->fetchrow_array;
			my $ptmName=$psiName || $interName || $synName || $des;
			push @ptmNames,$ptmName;
		}
		$sthVM->finish;
		if ($action eq 'export') {$pepPtmStrg='Modifications allowed: '.join(', ',sort{lc($a) cmp lc($b)} @ptmNames);}
		else {
			$pepPtmStrg="<A href=\"javascript:void(null)\" onmouseover=\"popup('<B>Modifications allowed:</B><BR>&nbsp;&bull;".join('<BR>&nbsp;&bull;',sort{lc($a) cmp lc($b)} @ptmNames)."')\" onmouseout=\"popout()\">Some modifications allowed<SUP>*</SUP></A>";
		}
		$pepPtmStrg.=' (Other modified and unmodified forms not allowed)' if $ptmAllowed==-2;
	}
	#}
	my %convPtmParam=('ambiguous'=>'delocalized','exclude'=>'excluded');
	#my $ptmQuantifStrg=(!$refInfo->{PTM_POS})? '' : ($refInfo->{PTM_POS}[0]=~/PRS:(.+)/)? "$modifName-site positions are <B>confirmed</B> if PRS probability &ge; <B>$1%</B>, others are <B>$convPtmParam{$refInfo->{PTM_POS}[1]}</B>" : ($refInfo->{PTM_POS}[0]==1)? "$modifName-sites are delocalized" : "$modifName-sites are not delocalized";

	#<Bias correction
	my %fdrMethods=('FDR-BH'=>'Benjamini-Hochberg',
					'FDR-ABH'=>'Benjamini-Hochberg (Adaptative)',
					'FWER-BONF'=>'Bonferroni',
					'Qvalue'=>'Qvalue (Storey et al.)'
				   );


	####################################
	####<Exporting summary to Excel>####
	####################################
	if ($action eq 'export') {
		my $worksheet1=$workbook->add_worksheet('Summary');
		my $xlsRow=0;
		$worksheet1->set_column(0,0,50);
		$worksheet1->set_column(1,1,40);
		$worksheet1->set_column(2,2,30);
		$worksheet1->merge_range($xlsRow,0,$xlsRow,2,$title,$itemFormat{'title'});

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
		#<Labeling
		$worksheet1->write_string(++$xlsRow,0,'Labeling :',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$labelingName{$labelType},$itemFormat{'mergeColText'});

		##<Parameters>##
		$xlsRow++;
		$worksheet1->merge_range(++$xlsRow,0,$xlsRow,2,'Quantification parameters',$itemFormat{'mergeColHeader'});
		$worksheet1->write_string(++$xlsRow,0,'Quantification method :',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,'SSP Analysis',$itemFormat{'mergeColText'});
		#<States
		foreach my $statePos (sort{$a<=>$b} keys %stateInfo) {
			$worksheet1->write_string(++$xlsRow,0,"State #$statePos :",$itemFormat{'headerR'});
			my $replicateNames=($stateInfo{$statePos}{'NAME_LONG'})? $stateInfo{$statePos}{'NAME_LONG'} : '';
			if ($stateInfo{$statePos}{'POPUP'}) {
				$worksheet1->write_string($xlsRow,1,"$stateInfo{$statePos}{NAME}$replicateNames",$itemFormat{'text'});
				$worksheet1->write_string($xlsRow,2,$stateInfo{$statePos}{'POPUP'},$itemFormat{'textWrap'});
			}
			else {$worksheet1->merge_range($xlsRow,1,$xlsRow,2,"$stateInfo{$statePos}{NAME}$replicateNames",$itemFormat{'mergeColText'});}
		}
		#<Features counted
		$worksheet1->write_string(++$xlsRow,0,'Feature counted :',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$featureStrg,$itemFormat{'mergeColText'});

		#<Peptide selection
		my $numPepLines=3;
		#if ($ptmQuantifStrg) {
		#	$numPepLines++;
		#	$ptmQuantifStrg="\n•".$ptmQuantifStrg;
		#}
		$worksheet1->set_row(++$xlsRow,13*$numPepLines);
		$worksheet1->write_string($xlsRow,0,'Peptide selection :',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,"•$pepRangeStrg.\n•$pepMissCutStrg.\n•$pepPtmStrg.\n•$rescuedStrg",$itemFormat{'mergeColText'});

		##<Export settings>##
		$xlsRow++;
		$worksheet1->merge_range(++$xlsRow,0,$xlsRow,2,'Export settings',$itemFormat{'mergeColHeader'});
		$worksheet1->set_row(++$xlsRow,13 * scalar keys %dispSets);
		$worksheet1->write_string($xlsRow,0,'Sets exported :',$itemFormat{'headerR'});
		my $disSetsStrg='';
		foreach my $sPos (sort{$a<=>$b} keys %dispSets) {
			$disSetsStrg.="\n" if $disSetsStrg;
			my $setLabel='';
			my $setSize=0;
			foreach my $pos (split('\+',$refInfo->{SETS}[$sPos-1])) {
				$setLabel.=' + ' if $setSize;
				$setLabel.=$stateInfo{$pos}{'NAME'};
				$setSize++;
			}
			$disSetsStrg.="•$setLabel";
		}
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$disSetsStrg,$itemFormat{'mergeColText'});
		$worksheet1->write_string(++$xlsRow,0,'Best delta (%) ≥',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,"$dispDelta",$itemFormat{'mergeColText'});
		$worksheet1->write_string(++$xlsRow,0,'p-value ≤',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$dispPvalue,$itemFormat{'mergeColText'});
		$worksheet1->write_string(++$xlsRow,0,'Feature count',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$dispMeanPep,$itemFormat{'mergeColText'});
		$worksheet1->write_string(++$xlsRow,0,'Sort by :',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$sortOptions{$dispSort}->[0],$itemFormat{'mergeColText'});
		if ($restrictListID) {
			my ($listStrg)=$dbh->selectrow_array("SELECT CONCAT(CL.NAME,' > ',C.NAME) FROM CATEGORY C,CLASSIFICATION CL WHERE C.ID_CLASSIFICATION=CL.ID_CLASSIFICATION AND ID_CATEGORY=$restrictListID");
			$worksheet1->write_string(++$xlsRow,0,'Restrict to proteins in List :',$itemFormat{'headerR'});
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$listStrg,$itemFormat{'mergeColText'});
		}

		return %stateInfo;
	} # end of &printQuantificationSummary when action='export'


	####################################
	####<Displaying summary in HTML>####
	####################################
	print qq
|<SCRIPT type="text/javascript">
var statesName={
|;
	foreach my $statePos (1..$numStates) {
		print "\t$statePos:'$stateInfo{$statePos}{NAME}'";
		print ',' if $statePos < $numStates;
		print "\n";
	}
	print qq
|};
// Names used for replicates (used by ajaxDrawSelectedProteins and showPeptideCount)
var replicatesLabel={
|;
	my $replicCount=0;
	foreach my $repKey (keys %bioRepLabel) {
		$replicCount++;
		$bioRepLabel{$repKey}=~s/'/\\'/g;		
		print "\t'$repKey':'$bioRepLabel{$repKey}'";
		print ',' if $replicCount < $replicPos;
		print "\n";
	}
	print qq
|};
// List of analyses corresponding to each state
var anaInStates={
|;
	foreach my $statePos (1..$numStates) {
		print "\t$statePos:'$anaInStates{$statePos}'";
		print ',' if $statePos < $numStates;
		print "\n";
	}
	print qq
|};

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
/* OBSOLETE -->
function getDomOffset( Obj, Prop ) {
	var iVal = 0;
	while (Obj && Obj.tagName != 'BODY') {
		eval('iVal += Obj.' + Prop + ';');
		Obj = Obj.offsetParent;
	}
	return iVal;
}
*/
function getElementPosition(e) {
	var left=0;
	var top=0;
	while (e.offsetParent != undefined && e.offsetParent != null) { //Adding parent item position
		left += e.offsetLeft + (e.clientLeft != null ? e.clientLeft : 0);
		top += e.offsetTop + (e.clientTop != null ? e.clientTop : 0);
		e = e.offsetParent;
	}
	return [top,left];
}
|;
	if ($action eq 'summary') {
		print qq
|</SCRIPT>
<BR>
|;
	}
	else {
		print qq
|function extendSelection(chkBox,selIdx) {
	for (let i=selIdx+1; i<chkBox.length; i++) {
		if (chkBox[i].checked != chkBox[selIdx].checked) {chkBox[i].checked = chkBox[selIdx].checked;}
		else {break;} // stop propagation if a box is already in selected check status
	}
}
function updateSortBy() {
	if ('$view'=='graph') {
		if (document.getElementById('listFromGraph') && document.getElementById('protListDIV').style.display != 'none') { // a protein list is displayed from graph
			ajaxListSelectedProteins(lastListFromGraph[0],lastListFromGraph[1]); // re-fetch list to update sort order
		}
	}
	else {document.displayForm.submit();} // view=list => resubmit form
}
var selGoID=0;
var goAspects=new Object();
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
		for (let i=0; i < checkBoxList.length; i++) {checkBoxList[i].checked=chkStatus;}
	}
	else {checkBoxList.checked=chkStatus;}
}
function exportProteins() {
	var defAction=document.displayForm.ACT.value;
	document.displayForm.ACT.value='export';
	document.displayForm.submit();
	document.displayForm.ACT.value=defAction;
}
|;
		if ($view eq 'graph') {
			my $dispSetPosStrg=join(',',sort{$a<=>$b} keys %dispSets);
			#my $maxRatioAllowed=($ratioType eq 'Ratio')? scalar @{$refInfo->{RATIOS}} : (scalar @{$refInfo->{STATES}})-1;
			print qq
|//GP datasetIdx to
var dataset2SetPos=[$dispSetPosStrg];
// Color list for plot
var hColors=['#E18B6B','#95B9C7','#7E2217','#9A9A9A','#8AFB17','#FBB917','#F660AB'];
var hColorIdx=0;
function proteinLink(dSetIdx,identifier,protID) {
	//var setPos=dataset2SetPos[dSetIdx];
	//sequenceView(protID,anaInStates[statePos]);
	ajaxDrawSelectedProteins({dSetIdx:[protID]});
}
// AJAX --->
var lastListFromGraph=[]; // needed to recall ajaxListSelectedProteins if "Sort by" is changed in graph view
function ajaxListSelectedProteins(selectedPoints,thresholds) {
	lastListFromGraph=[selectedPoints,thresholds];
	ajaxProcessSelectedProteins('ajaxListProt',document.getElementById('protListDIV'),selectedPoints);
}
function ajaxDrawSelectedProteins(selectedPoints,thresholds) {
	document.getElementById('displayDIV').style.display='none'; // in case popup is displayed
	// Adjusting displayDIV position to graph (graph can move due to summary option display)
	var mainGraphDiv=document.getElementById('mainGraphDIV');
	//var divX = getDomOffset(mainGraphDiv,'offsetLeft') + mainGraphDiv.offsetWidth + 2;
	//var divY = getDomOffset(mainGraphDiv,'offsetTop');
	var rect = mainGraphDiv.getBoundingClientRect();
	var divX=rect.right + document.body.scrollLeft + 2;
	var divY=rect.top + document.body.scrollTop;
	var graphDiv=document.getElementById('protGraphDIV');
	graphDiv.style.left = divX + 'px';
	graphDiv.style.top = divY + 'px';
	ajaxProcessSelectedProteins('ajaxDrawProt',graphDiv,selectedPoints);
}
function ajaxProcessSelectedProteins(action,targetDiv,selectedPoints) {
	// Updating contents
	prevImgID=null; // reset any [-] image from previous call
	targetDiv.innerHTML="<BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<FONT class=\\"title3\\">Fetching data. Please wait...</FONT><BR>&nbsp;&nbsp;&nbsp;&nbsp;<IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR>";
	targetDiv.style.display='block';
	//Parameters (extracting list of proteins)
	var proteinList=[];
	for (let gr in selectedPoints) {
		for (let i=0; i<selectedPoints[gr].length; i++) {
			proteinList.push(selectedPoints[gr][i]);
		}
	}
	var paramStrg='ACT='+action+'&id_quantif=$selQuantifID&view=$view&dispSetPosAjax=$dispSetPosStrg&sort='+document.getElementById('sort').value+'&pValue=$dispPvalue&delta=$dispDelta&meanPep=$dispMeanPep&selProt='+proteinList.join(',');

	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("POST","$promsPath{cgi}/showSSProtQuantification.cgi",true);
	//Send the proper header information along with the request
	XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
	//XHR.setRequestHeader("Content-length", paramStrg.length);
	//XHR.setRequestHeader("Connection", "close");
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			targetDiv.innerHTML='<INPUT type="hidden" id="listFromGraph" value="1"/>'; // listFromGraph: flag
			if (action=='ajaxDrawProt') {
				var codeParts=XHR.responseText.split('#==========#');
				targetDiv.innerHTML+=codeParts[0]; // HTML part
				eval(codeParts[1]); // javascript part
			}
			else {targetDiv.innerHTML+=XHR.responseText;}
			if (action=='ajaxListProt') {
				targetDiv.scrollIntoView({block:"start",inline:"nearest",behavior:"smooth"});
			}
		}
	}
	XHR.send(paramStrg);
}
function ajaxSearchConvertIdentifier(graphSearch,graphSearchArgs,searchTextIdx) { // (graph lib search function,array of function arguments,index of search text in array). To be called at end of convertion function
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxConvIdent&projectID=$projectID&quantifList=$selQuantifID&TEXT="+encodeURIComponent(graphSearchArgs[searchTextIdx]),true); //search text
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
var prevImgID=null;
function peptideDistribution(e,img,protID) {
	var displayDiv=document.getElementById('displayDIV');
	if (img.id==prevImgID) { // second click on img => close ajax popup window
		displayDiv.style.display='none';
		img.src='$promsPath{images}/plus.gif';
		prevImgID=null;
		return;
	}
	if (prevImgID && document.getElementById(prevImgID)) {document.getElementById(prevImgID).src='$promsPath{images}/plus.gif';}
	img.src="$promsPath{images}/minus1.gif";
	prevImgID=img.id;
	var divX,divY;
	if (e) { // called from protein list
		divX = (isNav)? e.pageX : event.clientX + document.body.scrollLeft; //divX-=5;
		divY = (isNav)? e.pageY : event.clientY + document.body.scrollTop; //divY+=10;
	}
	else {
		var graphDiv=document.getElementById('bpDIV_'+protID);
		//divX = getDomOffset(graphDiv,'offsetLeft') + graphDiv.offsetWidth + 2;
		//divY = getDomOffset(graphDiv,'offsetTop');
		var rect = graphDiv.getBoundingClientRect();
		divX=rect.right + document.body.scrollLeft + 2;
		divY=rect.top + document.body.scrollTop;
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
	XHR.open("GET","$promsPath{cgi}/showSSProtQuantification.cgi?ACT=ajaxPepDistrib&id_quantif=$selQuantifID&id_prot="+protID,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			var codeParts=XHR.responseText.split('#==========#');
			infoDiv.innerHTML=codeParts[0]; // HTML part
			eval(codeParts[1]); // javascript part
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
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxRestrictList&projectID=$projectID&restrictList=$restrictListID&submitForm=displayForm",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			document.getElementById('restrictDIV').innerHTML=XHR.responseText;
		}
	}
	XHR.send(null);
}
</SCRIPT>
<FORM name="displayForm" method="POST">
<INPUT type="hidden" name="ACT" value="result">
<INPUT type="hidden" name="id_ana" value="$anaID">
<INPUT type="hidden" name="id_quantif" value="$selQuantifID">
<INPUT type="hidden" name="displayRes" value="1">
|;
	}
	my $colspanStrg=($action eq 'summary')? '' : 'colspan=2'; #($existReplicates)? ' colspan=3' : ' colspan=2';
	print qq
|<TABLE bgcolor=$darkColor border=0>
<TR><TH nowrap align="right">Labeling :</TH><TH bgcolor="$lightColor"$colspanStrg align="left" valign="top">&nbsp;$labelingName{$labelType}</TH></TR>
|;
	my $colSpanCorrel=($quantifStatus > 0)? 1 : 2; # +/- correlation buttons
	foreach my $statePos (sort{$a<=>$b} keys %stateInfo) {
		#print "<TR><TH nowrap align=right>&nbsp;State #$statePos:</TH><TD bgcolor=$lightColor width=450><B>$stateInfo{$statePos}{NAME}</B>$stateInfo{$statePos}{REPLIC}</TD></TR>\n";
		my $replicateNames=($stateInfo{$statePos}{'NAME_LONG'})? $stateInfo{$statePos}{'NAME_LONG'} : '';
		print "<TR><TH nowrap align=right>&nbsp;State #$statePos :</TH><TD bgcolor=$lightColor width=450 nowrap>&nbsp;";
		print "<A href=\"javascript:void(null)\" onmouseover=\"popup('$stateInfo{$statePos}{POPUP}')\" onmouseout=\"popout()\">" if $stateInfo{$statePos}{'POPUP'};
		print "<B>$stateInfo{$statePos}{NAME}$replicateNames</B>";
		print "</A>" if $stateInfo{$statePos}{'POPUP'};
		print "&nbsp;</TD>";
		if ($statePos==1 && $action ne 'summary') {
			print "<TH bgcolor=$lightColor rowspan=$numStates valign=top align=right><DIV style=\"max-height:100px;overflow-y:scroll\"><TABLE bgColor=$darkColor><TR><TH>State specificity</TH></TR>\n";
			my $setPos=0;
			my $prevSize=0;
			foreach my $set (@{$refInfo->{SETS}}) {
				$setPos++;
				my $setLabel='';
				my $setSize=0;
				my %setPosList;
				foreach my $pos (split('\+',$set)) {
					$setLabel.=' + ' if $setLabel;
					$setLabel.=$stateInfo{$pos}{'NAME'};
					$setSize++;
					$setPosList{$pos}=1; # hash instead of array because of view=list
				}
				my $setBgColor=($setPos % 2)? $lightColor : $darkColor;
				if ($setSize != $prevSize) {
					print "</TD></TR>\n" if $setSize > 1;
					print "<TR><TD bgcolor=\"$lightColor\" class=\"TH\">";
				}
				else {print "<BR>\n";}
				my $chkStrg='';
				if ($dispSets{$setPos}) {
					$chkStrg='checked';
					@{$dispSets{$setPos}}=($setLabel,\%setPosList); # record set label
				}
				print "<INPUT type=\"checkbox\" name=\"dispSets\" value=\"$setPos\" onclick=\"extendSelection(document.displayForm.dispSets,",($setPos-1),")\" $chkStrg>$setLabel&nbsp;";
				$prevSize=$setSize;
			}
			print "</TD></TR>\n</TABLE></DIV></TH>\n";
		}
		print "</TR>\n";
	}

	#$ptmQuantifStrg='<BR>&nbsp;&bull;'.$ptmQuantifStrg if $ptmQuantifStrg;
	print qq
|<TR><TH align=right nowrap valign=top>&nbsp;Feature counted :</TH><TD nowrap bgcolor="$lightColor" $colspanStrg>&nbsp;$featureStrg&nbsp;</TD></TR>
<TR><TH align=right nowrap valign=top>&nbsp;Peptide selection :</TH><TD nowrap bgcolor="$lightColor" $colspanStrg>&nbsp;&bull;$pepRangeStrg&nbsp;&nbsp;&nbsp;&bull;$pepMissCutStrg&nbsp;&nbsp;&nbsp;&bull;$pepPtmStrg&nbsp;&nbsp;&nbsp;&bull;$rescuedStrg&nbsp;</TD></TR>
|; # $ptmQuantifStrg
	if ($action eq 'summary') {
		my $statusStrg=($quantifStatus==-2)? '<FONT color=#DD0000>Failed</FONT> (Click on "Monitor Quantifications" for more information)' : ($quantifStatus==-1)? 'Not launched yet' : ($quantifStatus==0)? 'On-going' : 'Finished';
		$updateDate=&promsMod::formatDate($updateDate);
		print qq
|<TR><TH align=right nowrap valign=top>&nbsp;Creation date :</TH><TD bgcolor="$lightColor" $colspanStrg valign=top>&nbsp;$updateDate by <B>$updateUser</B></TD></TR>
<TR><TH align=right nowrap valign=top>&nbsp;Status :</TH><TH align="left" bgcolor="$lightColor" $colspanStrg valign=top>&nbsp;$statusStrg</TD></TR>
</TABLE>
<BR>
|;
		$dbh->disconnect;

		if ($quantifStatus==-1 || $quantifStatus==0) { # on-going: reload every 5 sec
			print qq
|<SCRIPT type="text/javascript">
setTimeout(function(){parent.optionFrame.location.reload();},5000);
</SCRIPT>
|;
		}
		&endHTML;
		exit;
	}
	my $trueFilterStyle='style="font-weight:bold;color:#DD0000"';
	my ($selListView,$foldChPvalStyle)=($view eq 'list')? (' selected',$trueFilterStyle) : ('','style="font-weight:bold"'); # 'graph' is default
	$dispPvalue=($dispPvalue)? $dispPvalue : 0.05;
	print qq
|<TR><TH align=right nowrap valign=top>&nbsp;<FONT class="title2">Display :</FONT></TH><TD bgcolor="$lightColor" $colspanStrg><TABLE border=0 cellspacing=0 cellpadding=0>
	<TR><TD colspan=2 nowrap><FONT class="title3">&bull;View:</FONT><SELECT name="view" class="title3" onchange="document.displayForm.submit()"><OPTION value="graph">Graphical</OPTION><OPTION value="list"$selListView>List</OPTION></SELECT>
&nbsp;&nbsp;&nbsp;<FONT $foldChPvalStyle>&bull;Best delta:</FONT><INPUT type="text" name="delta" value="$dispDelta" size=2 $foldChPvalStyle/><FONT $foldChPvalStyle>%</FONT>&nbsp;&nbsp;&nbsp;<FONT $foldChPvalStyle>&bull;p-value &le;</FONT><INPUT type="text" name="pValue" value="$dispPvalue" size=5 $foldChPvalStyle/>
&nbsp;&nbsp;&nbsp;<FONT $trueFilterStyle>&bull;Feature count &ge;</FONT><INPUT type="text" name="meanPep" value="$dispMeanPep" size=2 $trueFilterStyle/>
&nbsp;&nbsp;&nbsp;<B>&bull;Sort by:</B><SELECT name="sort" id="sort" onchange="updateSortBy()" style="font-weight:bold">
|;
	#my @sortOptions=(['identifier','Identifiers'],['delta','Best delta'],['p-value','Adj. p-value'],['mean','Mean # peptides'],['mw','Molecular weight']);
	foreach my $option (sort{$sortOptions{$a}->[1]<=>$sortOptions{$b}->[1]} keys %sortOptions) {
		print "<OPTION value=\"$option\"";
		print ' selected' if $option eq $dispSort;
		print ">$sortOptions{$option}->[0]</OPTION>";
	}
	print qq
|</SELECT>&nbsp;&nbsp;&nbsp;</TD><TD rowspan=2><INPUT type="submit" value="Update" style="font-weight:bold"/></TD></TR>
	<TR><TD nowrap width=10%><FONT class="title3">&bull;Restrict to proteins in List:&nbsp;</FONT></TD><TD><DIV id="restrictDIV"><!--Custom list selection comes here with window.onload() --></DIV></TD></TR>
</TABLE></TD></TR>
|;
	if (scalar keys %{$refGoAnalyses}) {
		print qq
|<TR><TH align=right nowrap valign=top>&nbsp;<FONT class="title3">GO decoration :</FONT></TH><TD bgcolor="$lightColor" $colspanStrg><TABLE border=0 cellspacing=0 cellpadding=0>
	<TR><TD nowrap><SELECT id="goAna" class="title3" onchange="updateGoAspects(this.value)"><OPTION value="0">-= None =-</OPTION>
|;
		foreach my $goID (sort{lc($refGoAnalyses->{$a}[0]) cmp lc($refGoAnalyses->{$b}[0])} keys %{$refGoAnalyses}) {
			print "<OPTION value=\"$goID\">$refGoAnalyses->{$goID}[0]</OPTION>";
		}
		print qq
|</SELECT><FONT class="title3">:&nbsp;</FONT><SELECT id="goAspect" class="title3" style="visibility:hidden" onchange="ajaxUpdateGoTermList(this.value)"><OPTION value="">-= Select Aspect =-</OPTION></SELECT>&nbsp;</TD>
<TD><DIV id="goTermsDIV"></DIV></TD></TR></TD></TABLE></TD></TR>
|;
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



#sub getAnalysesInStates {
#	my ($refInfoStates,$format)=@_;
#	$format='array' unless $format;
#	my %anaInStates;
#	my $statePos=0;
#	foreach my $stateData (@{$refInfoStates}) {
#		$statePos++;
#		my %anaList;
#		my $anaCount=0;
#		my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$stateData);
#		foreach my $bioReplicate (split(/\./,$quantiObsIDs)) { # bio rep separator
#			foreach my $techReplicate (split(/&/,$bioReplicate)) { # tech rep separator
#				foreach my $pooledObs (split(/\+/,$techReplicate)) {
#					my @obsData=split(/:/,$pooledObs);
#					my $anaID;
#					if (scalar @obsData==4) {$anaID=$obsData[2];} # design post 07/2013 -> obsID,parQuantifID,anaID,targetPos
#					else {$anaID=$obsData[1];} # older design -> parQuantifID,anaID
#					$anaList{$anaID}=++$anaCount;
#				}
#			}
#		}
#		if ($format=~/str/) {$anaInStates{$statePos}=join(',',sort{$anaList{$a}<=>$anaList{$b}} keys %anaList);}
#		else {@{$anaInStates{$statePos}}=sort{$anaList{$a}<=>$anaList{$b}} keys %anaList;}
#	}
#	return %anaInStates;
#}
#
sub sortProteins {
	my ($sort,$refProteinInfo,$refQuantifValues)=@_;
	my ($pID1,$modStrg1)=($a=~/^(\d+)(.*)/); $modStrg1='' unless $modStrg1;
	my ($pID2,$modStrg2)=($b=~/^(\d+)(.*)/); $modStrg2='' unless $modStrg2;
	if ($sort eq 'identifier') {lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))}
	elsif ($sort eq 'mw') {$refProteinInfo->{$pID2}[1] <=> $refProteinInfo->{$pID1}[1] || lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))}
	elsif ($sort eq 'delta') {$refQuantifValues->{$b}{DELTA_PC} <=> $refQuantifValues->{$a}{DELTA_PC} || $refQuantifValues->{$a}{PVAL_ADJ} <=> $refQuantifValues->{$b}{PVAL_ADJ} || lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))}
	elsif ($sort eq 'p-value') {$refQuantifValues->{$a}{PVAL_ADJ} <=> $refQuantifValues->{$b}{PVAL_ADJ} ||$refQuantifValues->{$b}{DELTA_PC} <=> $refQuantifValues->{$a}{DELTA_PC} ||  lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))}
	elsif ($sort eq 'mean') {$refQuantifValues->{$b}{SET_MEAN} <=> $refQuantifValues->{$a}{SET_MEAN} || lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))}
	# default to 'set'
	else {length($refQuantifValues->{$a}{SET}) <=> length($refQuantifValues->{$a}{SET}) || &promsMod::sortSmart($refQuantifValues->{$a}{SET},$refQuantifValues->{$b}{SET}) || lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))}
}

sub displaySSPAPlot {
	my ($refQuantifValues,$refProteinInfo,$allowHlight,$pvalueCode,$numProtUsed,$featureName)=@_;
	print qq
|<SCRIPT type="text/javascript">
var GP;
window.onload=function() {
	ajaxUpdateRestrict();
	GP=new genericPlot({div:'mainGraphDIV',width:600,height:550,
						axisX:{title:'Best delta (%)',forceTo0:2}, // ,zoomable:true
						axisY:{title:'-Log10(p-value)',forceTo0:2}, // ,zoomable:true
						zoomable:true,
						axisClosure:true,
						pointOpacity:0.7,
						pointLabel:setPointPopupText,
						pointOnclick:proteinLink,
						pointOnList:[['Draw',ajaxDrawSelectedProteins],['List',ajaxListSelectedProteins]],
						allowHighlight:$allowHlight,
						searchable:{text:'Extended search',externalSearch:ajaxSearchConvertIdentifier},
						convertValue:function(axis,thVal) {if (axis=='X') {return thVal} else {return -1*Math.log(thVal)/2.302585093;}} // -log10
						});
	GP.addThreshold({axis:'X',label:'Best delta threshold',value:$dispDelta,color:'#0B0',editable:true});
	GP.addThreshold({axis:'Y',label:'p-value threshold',value:$dispPvalue,color:'#F00',editable:true});
|;
	my $setIdx=0;
	foreach my $setPos (sort{$a<=>$b} keys %dispSets) {
		print "\tGP.addDataSet($setIdx,{name:'$dispSets{$setPos}[0]',sizeName:'Features in set',sizeRule:{min:1,max:200,ratio:1}});\n";
		$setIdx++;
	}
	my $log10=log(10);
	$setIdx=0;
	foreach my $setPos (sort{$a<=>$b} keys %dispSets) {
		my $count=0;
		foreach my $protID (keys %{$refQuantifValues->{$setPos}}) {
			next unless $refQuantifValues->{$setPos}{$protID}{DELTA_PC}; # $refQuantifValues->{$setPos}{$protID} is always defined for all primary states due to {MEAN} value
			my ($protID,$modStrg)=($protID=~/^(\d+)(.*)/); $modStrg='' unless $modStrg;
			if ($count==0) {
				print "\tGP.addDataAsString($setIdx,'";
			}
			else {print ";";}
			$count++;
			#>Set missing p-values=1!!!
			if (!defined $refQuantifValues->{$setPos}{$protID}{$pvalueCode}) {
				$refQuantifValues->{$setPos}{$protID}{$pvalueCode}=1;
			}
			elsif ($refQuantifValues->{$setPos}{$protID}{$pvalueCode}==0) {
				$refQuantifValues->{$setPos}{$protID}{$pvalueCode}=$minPvalue;
			}
			my $mean=(sprintf "%.1f",$refQuantifValues->{$setPos}{$protID}{SET_MEAN})*1;
			print "$refProteinInfo->{$protID}[0]$modStrg,$protID,$refQuantifValues->{$setPos}{$protID}{DELTA_PC},",(-log($refQuantifValues->{$setPos}{$protID}{$pvalueCode})/$log10),",",$mean;
			if ($count==100) {
				print "');\n";
				$count=0;
			}
		}
		if ($count > 0) {
			print "');\n";
		}
		$setIdx++;
	}
	print qq
|	GP.draw();
}
function setPointPopupText(dp,type) {
	var infoStrg=dp.label;
	if (type=='max') {
		if (dp.dataSet.params.chart.dataSets.length > 1) infoStrg+='\\nSet='+dp.dataSet.params.name;
		infoStrg+='\\nDelta='+dp.x+'%\\np-value='+(Math.pow(10,-1*dp.y*1).toPrecision(4))+'\\nMean='+dp.size+' $featureName';
	}
	return infoStrg;
}
var peptideCount={};
function showPeptideCount(grID) { // identifier:protID:statePos
	if (prevImgID) { // in case previous use of displayDIV was for peptidePlot
		document.getElementById(prevImgID).src='$promsPath{images}/plus.gif';
		prevImgID=null;
	}
	var grInfo=grID.split(':');
	var pepStrg='<INPUT type="button" value="Close" class="font11" onclick="document.getElementById(\\'displayDIV\\').style.display=\\'none\\';"><TABLE bgcolor=$darkColor><TR><TH colspan=2>'+grInfo[0]+' in '+statesName[grInfo[2]]+'</TH></TR>';
	var dataList=peptideCount[grID].split(',');
	for (let i=0; i<dataList.length; i++) {
		var replicInfo=dataList[i].split('='); // replicCode=value
		pepStrg+='<TR><TH align=right nowrap>&nbsp;'+replicatesLabel[replicInfo[0]]+'</TH><TD bgcolor=$lightColor>&nbsp;'+replicInfo[1]+'</TD></TR>';
	}
	pepStrg+='</TABLE>';
	document.getElementById('infoDIV').innerHTML=pepStrg;
	// Adjusting displayDIV position to boxPlot DIV
	var barPlotDiv=document.getElementById('bpDIV_'+grInfo[1]);
	//var divX = getDomOffset(barPlotDiv,'offsetLeft') + barPlotDiv.offsetWidth + 2;
	//var divY = getDomOffset(barPlotDiv,'offsetTop');
	var rect = barPlotDiv.getBoundingClientRect();
	var divX=rect.right + document.body.scrollLeft + 2;
	var divY=rect.top + document.body.scrollTop;
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = divX + 'px';
	displayDiv.style.top = divY + 'px';
	displayDiv.style.display='block';
}
function exportProteinData(format,protDivID,protAlias) { // protDivID,protAlias defined only for image
	if (format=='xls') { // data export
		/* Serialize peptideCount */
		var protData='';
		for (let grID in peptideCount) {
			if (protData != '') {protData+=';';}
			protData+=grID+'>'+peptideCount[grID];
		}
		document.exportProtDataForm.DATA.value=protData;
		/* Serialize statesName & replicatesLabel (only once!) */
		if (!document.exportProtDataForm.STATES.value) {
			var states='';
			for (let statePos in statesName) {
				if (states != '') {states+=';';}
				states+=statePos+'='+statesName[statePos];
			}
			document.exportProtDataForm.STATES.value=states;

			var replicLabel='';
			for (let replic in replicatesLabel) {
				if (replicLabel != '') {replicLabel+=';';}
				replicLabel+=replic+'='+replicatesLabel[replic];
			}
			document.exportProtDataForm.REPLICATES.value=replicLabel;
		}
		document.exportProtDataForm.submit();
	}
	else { // image export
		exportSVGtoImg(protDivID,protAlias,'./exportSVG.cgi',format); // WARNING: defined in chartLibrary2.js
	}
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
			addHighlighting(GP,termData[2],hColors[hColorIdx],{'-1':XHR.responseText.split(';')});
			hColorIdx++;
			if (hColorIdx==hColors.length) hColorIdx=0;
		}
	}
	XHR.send(null);
}
</SCRIPT>
</CENTER>
<FORM name="exportProtDataForm" method="POST" target="_blank"><!-- Form used to export counts for a single protein to Excel -->
<INPUT type="hidden" name="ACT" value="exportProtData">
<INPUT type="hidden" name="STATES" value="">
<INPUT type="hidden" name="REPLICATES" value="">
<INPUT type="hidden" name="DATA" value="">
</FORM>
<FONT class="title3">$numProtUsed proteins selected.</FONT><BR>
<DIV id="mainGraphDIV" style="float:left;background-color:#FFFFFF;border:solid 3px $darkColor;padding:2px;"></DIV>
<DIV id="protListDIV" style="clear:both;max-height:600;overflow:auto"></DIV>
<DIV id="protGraphDIV" style="position:absolute;max-height:600;overflow:auto"></DIV>
|;
}

sub printProteinList { # Global: %dispSets
	my ($refLabelingInfo,$refStateInfo,$refQuantifValues,$refProteinInfo,$pvalueCode,$numProtUsed)=@_;
	my $numStates=scalar @{$refLabelingInfo->{'STATES'}};
	my $numTotProteins=scalar keys %{$refProteinInfo};
	#my $protString=($quantifModID)? "&nbsp;$numProtUsed isoforms&nbsp;<BR>&nbsp;$numTotProteins proteins&nbsp;" : "&nbsp;$numProtUsed proteins&nbsp;";
	my $protString="&nbsp;$numProtUsed proteins&nbsp;";
	my $disabSave=($projectAccess eq 'guest')? ' disabled' : '';
	my $clearButtonStrg=($view eq 'list')? '' : '<INPUT type="button" value="Clear list" style="font-weight:bold" onclick="document.getElementById(\'protListDIV\').style.display=\'none\'">';
	my $numColumns= 7 + $numStates;
	print qq
|<FORM name="protForm" method="POST">
<TABLE border=0 cellspacing=0 cellpadding=2>
<TR><TD colspan=$numColumns><DIV id="saveProtDIV" style="display:none;"></DIV></TD></TR>
<TR><TD colspan=$numColumns>$clearButtonStrg
<INPUT type="button" value="Check all" onclick="checkAllProteins(true)"/>
<INPUT type="button" value="Uncheck all" onclick="checkAllProteins(false)"/>
<INPUT type="button" id="saveFormBUTTON" value="Save proteins..." onclick="ajaxManageSaveProteins('getThemes')"$disabSave/>
|;
	print "<INPUT type=\"button\" value=\"Export data\" onclick=\"exportProteins()\"/>\n" unless $action eq 'ajaxListProt';
	print qq
|</TD></TR>
<TR bgcolor="$darkColor">
<TH class="rbBorder" align=left rowspan=2>&nbsp;$protString&nbsp;</TH><TH class="rbBorder" rowspan=2>&nbsp;Gene&nbsp;</TH>
<TH class="rbBorder" rowspan=2>&nbsp;Specificity&nbsp;</TH><TH class="rbBorder" rowspan=2>&nbsp;Best delta (%)&nbsp;</TH><TH class="rbBorder" rowspan=2>&nbsp;p-value&nbsp;</TH>
<TH class="rbBorder" colspan=$numStates>&nbsp;Feature count&nbsp;</TH>
<TH class="rbBorder" rowspan=2 nowrap>&nbsp;MW<SMALL> kDa</SMALL>&nbsp;</TH><TH class="bBorder" rowspan=2 width=700 nowrap>&nbsp;Description - Species&nbsp;</TH>
</TR>
<TR bgcolor="$darkColor">
|;
	foreach my $statePos (1..$numStates) {
		print "<TH class=\"rbBorder\" nowrap>&nbsp;$refStateInfo->{$statePos}{NAME}&nbsp;</TH>\n";
	}
	print qq
|</TR>
<TR bgcolor="$darkColor">
|;
	my $bgColor=$lightColor;
	my $numDispItems=0;
	foreach my $protID (sort{&sortProteins($dispSort,$refProteinInfo,$refQuantifValues)} keys %{$refQuantifValues}) { #%{$refProteinInfo}
		my $anaID = $refProteinInfo->{$protID}[-1];
		$numDispItems++;
		print "<TR bgcolor=\"$bgColor\" valign=top><TD class=\"TH\" nowrap><INPUT type=\"checkbox\" name=\"chkProt\" value=\"$protID\"/><A href=\"javascript:sequenceView($protID,$anaID)\">$refProteinInfo->{$protID}[0]</A><IMG id=\"pepDataList_$protID\" src=\"$promsPath{images}/plus.gif\" onclick=\"peptideDistribution(event,this,$protID)\">&nbsp;</TD>\n";
		# Gene(s)
		if (scalar @{$refProteinInfo->{$protID}[5]} > 1) {
			print "<TH>&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-",join('<BR>&nbsp;&nbsp;-',@{$refProteinInfo->{$protID}[5]}[1..$#{$refProteinInfo->{$protID}[5]}]),"</B>')\" onmouseout=\"popout()\">",$refProteinInfo->{$protID}[5][0],"</A>&nbsp;</TH>";
		}
		elsif ($refProteinInfo->{$protID}[5][0]) {print '<TH>&nbsp;',$refProteinInfo->{$protID}[5][0],'&nbsp;</TH>';}
		else {print '<TH>-</TH>';} # no gene mapped
		my $setPos=$refQuantifValues->{$protID}{SET};
		# Specificity
		print "<TH>&nbsp;$dispSets{$setPos}[0]&nbsp;</TH>";
		# Best delta
		my $delta=1*(sprintf "%.1f",$refQuantifValues->{$protID}{DELTA_PC});
		print "<TH>&nbsp;$delta&nbsp;</TH>";
		# p-value
		my $pValueStrg=($refQuantifValues->{$protID}{$pvalueCode}==0)? '~0' : ($refQuantifValues->{$protID}{$pvalueCode}>=0.01)? sprintf "%.2f",$refQuantifValues->{$protID}{$pvalueCode} : sprintf "%.1e",$refQuantifValues->{$protID}{$pvalueCode};
		if ($refQuantifValues->{$protID}{$pvalueCode}<0.05) {print "<TH>&nbsp;$pValueStrg&nbsp;</TH>";}
		else {print "<TD align=\"center\">&nbsp;$pValueStrg&nbsp;</TD>";}
		# States mean (feature count)
		foreach my $statePos (1..$numStates) {
			my $numFeat=1*(sprintf "%.1f",$refQuantifValues->{$protID}{MEAN}{$statePos});
			if ($numFeat) {
				if ($dispSets{$setPos}[1]{$statePos}) {print "<TH>&nbsp;<A href=\"javascript:sequenceView($protID,anaInStates[$statePos])\">$numFeat</A>&nbsp;</TH>";} # state is part of set
				else {print "<TD align=\"center\">&nbsp;<A href=\"javascript:sequenceView($protID,anaInStates[$statePos])\">$numFeat</A>&nbsp;</TD>";}
			}
			else {print "<TD align=\"center\">&nbsp;$numFeat&nbsp;</TD>";}
		}
		# mw, des & species
		my $mass=(!$refProteinInfo->{$protID}[1])? '-' : $refProteinInfo->{$protID}[1];
		print "<TD class=\"TH\" align=right>$mass&nbsp;</TD><TD>$refProteinInfo->{$protID}[2] <U><I>$refProteinInfo->{$protID}[3]</I></U></TD>\n";
		print "</TR>\n";
		$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
	}
	print qq
|<TR><TH colspan=$numColumns align="left">End of list.</TH></TR>
</TABLE>
</FORM>
|;
} # end of printProteinList

sub exportProteinList { # Only for design or no-design labeled quantifs
###	my ($refLabelingInfo,$refStateInfo,$refQuantifValues,$refProteinInfo,$pvalueCode)=@_; #,$dispStdDev
###	my $invDispFoldChange=1/$dispFoldChange;
###	my $log2=log(2);
###	my $numTotProteins=scalar keys %{$refProteinInfo};
###	my $peptideTitle=($labelType eq 'SILAC')? 'Pept. sets' : 'Peptides';
###
###	####<Start printing>####
###	my $worksheet2=$workbook->add_worksheet('Results');
###	my $xlsRow=0;
###	my $xlsCol=0;
###
###	###<Headers>###
###	#<Identifiers
###	$worksheet2->set_column(0,0,30); # identifier col length
###	$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow+1,$xlsCol,'Gene & Synonyms',$itemFormat{'mergeRowHeader'});
###	#$worksheet2->merge_range($xlsRow,$xlsCol,$xlsRow+1,$xlsCol,"Proteins",$itemFormat{'mergeRowHeader'});
###	#<Ratios
###	my $ratioType=($refLabelingInfo->{'RATIO_TYPE'})? $refLabelingInfo->{'RATIO_TYPE'}[0] : 'Ratio';
###	my (%ratioColDelta); # %supRatios,
###	my $curRatioPos=0;
###	my $firstRatioUsed=0;
###	foreach my $ratioData (@{$refLabelingInfo->{RATIOS}}) {
###		$curRatioPos++;
###		next unless $dispSets{$curRatioPos};
###		$firstRatioUsed=$curRatioPos unless $firstRatioUsed;
###		my ($testStatePos,$refStatePos)=split(/\//,$ratioData);
###		my $ratioTag='';
###		if ($designID) {
###			$testStatePos=~s/#//;
###			$refStatePos=~s/#//;
###			if ($ratioData=~/%/) {
###				#$supRatios{$curRatioPos}=1;
###				$ratioTag='°';
###				$testStatePos=~s/%\d+//;
###				$refStatePos=~s/%\d+//;
###			}
###			$testStatePos=$condToState{$testStatePos};
###			$refStatePos=$condToState{$refStatePos};
###		}
###		$ratioColDelta{$curRatioPos}=($quantifType eq 'TNPQ')? 4 : ($ratioType eq 'Ratio')? 5 : 5; # S(imple/uper)Ratio    ratio/Log2 ratio/(low. bound./up. bound.)*/p-value/(+/-stdDev)
###		$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow,$xlsCol+$ratioColDelta{$curRatioPos},$refStateInfo->{$testStatePos}{'NAME'}.$ratioTag.'/'.$refStateInfo->{$refStatePos}{'NAME'}.$ratioTag,$itemFormat{'mergeColHeader'});
###		$worksheet2->write_string($xlsRow+1,$xlsCol,'Ratio',$itemFormat{'header'});
###		$worksheet2->write_comment($xlsRow+1,$xlsCol,'Ratio of 1000 or 0.001 means protein was not detected in one condition.');
###		$worksheet2->write_string($xlsRow+1,++$xlsCol,'Log2',$itemFormat{'header'});
###		if ($ratioType eq 'Ratio') {
###			$worksheet2->write_string($xlsRow+1,++$xlsCol,'Lower bound.',$itemFormat{'header'});
###			$worksheet2->write_comment($xlsRow+1,$xlsCol,"Lower boundary of ratio confidence interval: There a ".($refLabelingInfo->{'CONFIDENCE_LEVEL'}[0]*100)."% certainty that true ratio is included between lower and upper boundaries.");
###			$worksheet2->write_string($xlsRow+1,++$xlsCol,'Upper bound.',$itemFormat{'header'});
###			$worksheet2->write_comment($xlsRow+1,$xlsCol,"Upper boundary of ratio confidence interval: There a ".($refLabelingInfo->{'CONFIDENCE_LEVEL'}[0]*100)."% certainty that true ratio is included between lower and upper boundaries.");
###		}
###		$worksheet2->write_string($xlsRow+1,++$xlsCol,'p-value',$itemFormat{'header'});
###		$worksheet2->write_string($xlsRow+1,++$xlsCol,'Std. dev.',$itemFormat{'header'}) if $quantifType ne 'TNPQ';
###		if ($ratioType=~/S\w+Ratio/) {
###			$worksheet2->write_string($xlsRow+1,++$xlsCol,'Dist. pept. used',$itemFormat{'header'});
###			$worksheet2->write_string($xlsRow+1,++$xlsCol,'Pept. used',$itemFormat{'header'});
###		}
###	}
###	#<Peptides
###	if ($ratioType eq 'Ratio') {
###		$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow,$xlsCol+2,$peptideTitle,$itemFormat{'mergeColHeader'});
###		$worksheet2->write_string($xlsRow+1,$xlsCol,'Distinct used',$itemFormat{'header'});
###		$worksheet2->write_string($xlsRow+1,++$xlsCol,'Used',$itemFormat{'header'});
###		$worksheet2->write_string($xlsRow+1,++$xlsCol,'Total',$itemFormat{'header'});
###	}
###	else {
###		$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow+1,$xlsCol,"Total $peptideTitle",$itemFormat{'mergeRowHeader'});
###	}
###	#<MW, description & species
###	$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow+1,$xlsCol,'MW (kDa)',$itemFormat{'mergeRowHeader'});
###	$worksheet2->set_column(++$xlsCol,$xlsCol,80); # col length
###	$worksheet2->merge_range($xlsRow,$xlsCol,$xlsRow+1,$xlsCol,'Description',$itemFormat{'mergeRowHeader'});
###	$worksheet2->set_column(++$xlsCol,$xlsCol,30); # col length
###	$worksheet2->merge_range($xlsRow,$xlsCol,$xlsRow+1,$xlsCol,'Species',$itemFormat{'mergeRowHeader'});
###
###	###<Looping through proteins>###
###	$xlsRow++;
###	my $numDispItems=0;
###	my %dispProteins;
###	foreach my $modProtID (sort{&sortProteins($dispSort,$refProteinInfo,$refQuantifValues,$firstRatioUsed)} keys %{$refQuantifValues}) { #%{$refProteinInfo}
####last if $numDispItems >= 5;
###		my ($protID,$modResStrg)=($modProtID=~/^(\d+)(.*)/); # eg. 12345-Y40,T152 -> 12345, -Y40,T152
###		$modResStrg='' unless $modResStrg;
###		$xlsCol=0;
###		my $refUsedRatioPos=0;
###		#<Peptide filter
###		my $okPeptide=0;
###		foreach my $ratioPos (keys %dispSets) {
###			next if (!$refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos} || $refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos} < $dispNumPep); # may not be defined for +/-inf ratios
###			$okPeptide=1; # at least 1 ratio must pass the filter
###			$refUsedRatioPos=$ratioPos;
###			last;
###		}
###		next unless $okPeptide;
###
###		#<Ratio/p-value/Std dev filters
###		my $okFilters=0;
###		foreach my $ratioPos (keys %dispSets) {
###			unless ($okFilters) { # no ratioPos has passed filters
###				if (defined $refQuantifValues->{$modProtID}{"RATIO$ratioMQTag"}{$ratioPos}) {
###					if ($foldChangeType eq 'abs') {
###						$okFilters=(($refQuantifValues->{$modProtID}{"RATIO$ratioMQTag"}{$ratioPos} < 1 && $refQuantifValues->{$modProtID}{"RATIO$ratioMQTag"}{$ratioPos} <= $invDispFoldChange) || ($refQuantifValues->{$modProtID}{"RATIO$ratioMQTag"}{$ratioPos} >= 1 && $refQuantifValues->{$modProtID}{"RATIO$ratioMQTag"}{$ratioPos} >= $dispFoldChange))? 1 : 0;
###					}
###					elsif ($foldChangeType eq 'up') {
###						$okFilters=($refQuantifValues->{$modProtID}{"RATIO$ratioMQTag"}{$ratioPos} >= 1 && $refQuantifValues->{$modProtID}{"RATIO$ratioMQTag"}{$ratioPos} >= $dispFoldChange)? 1 : 0;
###					}
###					else { # down
###						$okFilters=($refQuantifValues->{$modProtID}{"RATIO$ratioMQTag"}{$ratioPos} <= 1 && $refQuantifValues->{$modProtID}{"RATIO$ratioMQTag"}{$ratioPos} <= $invDispFoldChange)? 1 : 0;
###					}
###				}
###				$okFilters=($okFilters && ($dispPvalue>=1 || (defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} && $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} <= $dispPvalue)))? 1 : 0;
###				if ($okFilters && $dispConfInt && defined $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}) { # undefined always passes filter!
###					$okFilters=(1.96*$refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos} <= $dispConfInt)? 1 : 0;
###				}
###			}
###		}
###		next unless $okFilters;
###		$numDispItems++;
###		$dispProteins{$protID}=1;
###		# Identifier
###		$worksheet2->write_string(++$xlsRow,$xlsCol,$refProteinInfo->{$protID}[0].$modResStrg,$itemFormat{'text'});
###		# Gene(s)
###		$worksheet2->write_string($xlsRow,++$xlsCol,join(',',@{$refProteinInfo->{$protID}[5]}),$itemFormat{'text'});
###		#<Ratios
###		foreach my $ratioPos (sort{$a<=>$b} keys %dispSets) {
###			$refUsedRatioPos=$ratioPos unless $refUsedRatioPos;
###			if (!$refQuantifValues->{$modProtID}{"RATIO$ratioMQTag"}{$ratioPos}) { # !$refQuantifValues->{$protID}{'RATIO'} ||
###				foreach my $i (0..$ratioColDelta{$ratioPos}) {
###					$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
###				}
###				next;
###			}
###			# Ratio
###			$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{"RATIO$ratioMQTag"}{$ratioPos},$itemFormat{'number'});
###			# log2 ratio
###			$worksheet2->write_number($xlsRow,++$xlsCol,log($refQuantifValues->{$modProtID}{"RATIO$ratioMQTag"}{$ratioPos})/$log2,$itemFormat{'number'});
###			if ($ratioType eq 'Ratio') {
###				# low conf limit
###				if ($refQuantifValues->{$protID}{'CONF_LOW'} && defined $refQuantifValues->{$modProtID}{'CONF_LOW'}{$ratioPos}) {
###					$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'CONF_LOW'}{$ratioPos},$itemFormat{'number'});
###				}
###				else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
###				# upper conf limit
###				if ($refQuantifValues->{$protID}{'CONF_UP'} && defined $refQuantifValues->{$modProtID}{'CONF_UP'}{$ratioPos}) {
###					$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'CONF_UP'}{$ratioPos},$itemFormat{'number'});
###				}
###				else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
###			}
###			# p-value
###			if (defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}) {
###				$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos},$itemFormat{'number'});
###			}
###			else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
###
###			if ($quantifType ne 'TNPQ') {
###				# Std dev
###				if ($refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}) {
###					$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos},$itemFormat{'number'});
###				}
###				else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
###			}
###			# SuperRatio peptides
###			if ($ratioType=~/S\w+Ratio/) {
###				# DIST_PEP_USED
###				if ($refQuantifValues->{$modProtID}{'DIST_PEP_USED'} && $refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$ratioPos} && $refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$ratioPos}) {
###					$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$ratioPos},$itemFormat{'number'});
###				}
###				else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
###				# NUM_PEP_USED
###				if ($refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$ratioPos} && $refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$ratioPos}) {
###					$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$ratioPos},$itemFormat{'number'});
###				}
###				else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
###			}
###		}
###		# Peptides
###		if ($refQuantifValues->{$modProtID}) {
###			if ($ratioType eq 'Ratio') {
###				# DIST_PEP_USED
###				if ($refQuantifValues->{$modProtID}{'DIST_PEP_USED'} && $refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$refUsedRatioPos} && $refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$refUsedRatioPos}) {
###					$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$refUsedRatioPos},$itemFormat{'number'});
###				}
###				else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
###				# NUM_PEP_USED
###				if ($refQuantifValues->{$modProtID}{'NUM_PEP_USED'} && $refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$refUsedRatioPos}) {
###					$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$refUsedRatioPos},$itemFormat{'number'});
###				}
###				else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
###			}
###			# NUM_PEP_TOTAL
###			if ($refQuantifValues->{$modProtID}{'NUM_PEP_TOTAL'} && $refQuantifValues->{$modProtID}{'NUM_PEP_TOTAL'}{$refUsedRatioPos}) {
###				$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'NUM_PEP_TOTAL'}{$refUsedRatioPos},$itemFormat{'number'});
###			}
###			else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
###		}
###		else {
###			if ($ratioType eq 'Ratio') {
###				$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
###				$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
###			}
###			$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
###		}
###		# MW, desc, species
###		my $mass=(!$refProteinInfo->{$protID}[1])? 0 : $refProteinInfo->{$protID}[1];
###		$worksheet2->write_number($xlsRow,++$xlsCol,$mass,$itemFormat{'number1d'});
###		$worksheet2->write_string($xlsRow,++$xlsCol,$refProteinInfo->{$protID}[2],$itemFormat{'textWrap'});
###		$worksheet2->write_string($xlsRow,++$xlsCol,$refProteinInfo->{$protID}[3],$itemFormat{'text'});
###	}
###	#<Identifier header (written last to get remaining number of proteins)
###	my $numTotQuantifItems=scalar keys %{$refQuantifValues};
###	my $protString=($quantifModID)? "$numDispItems/$numTotQuantifItems isoforms\n".(scalar keys %dispProteins)."/$numTotProteins proteins" : "$numDispItems/$numTotProteins proteins";
###	$worksheet2->merge_range(0,0,1,0,$protString,$itemFormat{'mergeRowHeader'});
}

sub exportProteinData {
	my (%stateNames,%bioRepLabel,%peptideCount,%numReplicates);
	foreach my $state (split(';',param('STATES'))) {
		my ($statePos,$stateLabel)=split('=',$state);
		$stateNames{$statePos}=$stateLabel;
	}
	foreach my $replic (split(';',param('REPLICATES'))) {
		my ($repKey,$repLabel)=split('=',$replic);
		$bioRepLabel{$repKey}=$repLabel;
	}
	my %seenState;
	foreach my $grData (split(';',param('DATA'))) {
		my ($grID,$grData)=split('>',$grData);
		my ($identifier,$protID,$statePos)=split(':',$grID);
		foreach my $repData (split(',',$grData)) {
			my ($repKey,$repCount)=split('=',$repData);
			push @{$peptideCount{$identifier}{$statePos}},$repCount;
			push @{$numReplicates{$statePos}},$repKey unless $seenState{$statePos}; # same list for all proteins
		}
		$seenState{$statePos}=1;
	}

	#################################
	####<Prepare export to Excel>####
	#################################
	my $timeStamp=strftime("%Y%m%d %H-%M",localtime);
	my $workbook=Spreadsheet::WriteExcel->new("-");
	$workbook->set_properties(title=>"State-Specific Protein Analysis",
							  author=>'myProMS server',
							  comments=>'Automatically generated with Perl and Spreadsheet::WriteExcel'
							  );
	$workbook->set_custom_color(40,189,215,255); # light light color  #BDD7FF (new myProMS colors V3.5+)
	my %itemFormat=(
			#title =>			$workbook->add_format(align=>'center',size=>18,bold=>1,border=>1),
			header =>			$workbook->add_format(align=>'center',size=>10,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			headerR =>			$workbook->add_format(align=>'right',valign=>'top',size=>10,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			mergeRowHeader =>	$workbook->add_format(align=>'center',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			mergeColHeader =>	$workbook->add_format(align=>'center',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			text =>				$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>0,border=>1),
			number =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1)
	);
	print header(-type=>"application/vnd.ms-excel",-attachment=>"SSPA_raw_counts_$timeStamp.xls");

	my $worksheet=$workbook->add_worksheet('SSPA - Raw counts');
	my $xlsRow=0;
	my $xlsCol=0;

	##<Headers (2 lines)
	$worksheet->set_column(0,0,30); # identifier col length
	$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow+1,$xlsCol,'Identifier',$itemFormat{'mergeRowHeader'});
	foreach my $statePos (sort{$a<=>$b} keys %stateNames) {
		my $endCol=$xlsCol+scalar @{$numReplicates{$statePos}};
		$worksheet->merge_range($xlsRow,$xlsCol+1,$xlsRow,$endCol,$stateNames{$statePos},$itemFormat{'mergeColHeader'}); # states name
		foreach my $repKey (@{$numReplicates{$statePos}}) {
			$worksheet->write_string($xlsRow+1,++$xlsCol,$bioRepLabel{$repKey},$itemFormat{'header'}); # replicates name
		}
		$xlsCol=$endCol;
	}
	$xlsRow++;

	##<Data
	foreach my $identifier (sort keys %peptideCount) {
		$xlsRow++;
		$xlsCol=0;
		$worksheet->write_string($xlsRow,$xlsCol,$identifier,$itemFormat{'text'});
		foreach my $statePos (sort{$a<=>$b} keys %stateNames) {
			foreach my $repCount (@{$peptideCount{$identifier}{$statePos}}) {
				$worksheet->write_number($xlsRow,++$xlsCol,$repCount,$itemFormat{'number'}); # counts
			}
		}
	}

	$workbook->close();

	#print header(-charset=>'utf-8'); warningsToBrowser(1);
	#print "<TABLE border=1>\n<TR><TH rowspan=2>Identifier</TH>";
	#foreach my $statePos (sort{$a<=>$b} keys %stateNames) {
	#	print "<TH colspan=",scalar @{$numReplicates{$statePos}},">$stateNames{$statePos}</TH>\n";
	#}
	#print "</TR>\n<TR>";
	#foreach my $statePos (sort{$a<=>$b} keys %stateNames) {
	#	foreach my $repKey (@{$numReplicates{$statePos}}) {
	#		print "<TH>$bioRepLabel{$repKey}</TH>\n";
	#	}
	#}
	#print "</TR>\n";
	#foreach my $identifier (sort keys %peptideCount) {
	#	print "<TR><TH>$identifier</TH>";
	#	foreach my $statePos (sort{$a<=>$b} keys %stateNames) {
	#		foreach my $repCount (@{$peptideCount{$identifier}{$statePos}}) {
	#			print "<TD align='center'>$repCount</TD>";
	#		}
	#	}
	#	print "</TR>\n";
	#}
	#print "</TABLE>\n";

	exit;
}


sub ajaxListSelectedProteins {
	my @selProteins=split(',',param('selProt'));
	my @setPosList=split(',',param('dispSetPosAjax'));
	foreach my $setPos (@setPosList) {@{$dispSets{$setPos}}=();}
#print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG

	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;
	my $projectID=&promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION');
	my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
	$projectAccess=${$userInfo[2]}{$projectID};

	###<Fetching quantification parameters>###
	my $sthQP2=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,P.NAME,P.CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD M,QUANTIFICATION_PARAMETER P WHERE Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND M.ID_QUANTIFICATION_METHOD=P.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$selQuantifID");
	my (%quantifParamInfo,%quantifParamCode);
	$sthQP2->execute;
	while (my ($paramID,$paramName,$paramCode)=$sthQP2->fetchrow_array) {
		@{$quantifParamInfo{$paramCode}}=($paramID,$paramName);
		$quantifParamCode{$paramID}=$paramCode;
	}
	$sthQP2->finish;

	###<Fetching quantification design>###
	my %labelingInfo;
	(my $quantifAnnot)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
	foreach my $infoStrg (@labelInfo) {
		my ($setting,$valueStrg)=split('=',$infoStrg);
		$valueStrg=~s/#//g;
		@{$labelingInfo{$setting}}=split(';',$valueStrg);
	}
	my %stateInfo;
	my $statePos=0;
	my $sthgetExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
	foreach my $stateData (@{$labelingInfo{STATES}}) {
		$statePos++;
		my ($numBioRep,$sampData,$expCondID)=split(',',$stateData);
		$sthgetExpCondName->execute($expCondID);
		($stateInfo{$statePos}{'NAME'})=$sthgetExpCondName->fetchrow_array;
	}
	$sthgetExpCondName->finish;

	#<Update global %setPosList
	my $setPos=0;
	foreach my $set (@{$labelingInfo{SETS}}) {
		$setPos++;
		my $setLabel='';
		my %setPosList;
		foreach my $pos (split('\+',$set)) {
			$setLabel.=' + ' if $setLabel;
			$setLabel.=$stateInfo{$pos}{'NAME'};
			$setPosList{$pos}=1; # hash instead of array because of view=list
		}
		@{$dispSets{$setPos}}=($setLabel,\%setPosList); # record set label
	}

	###<Protein data>###
	my %quantifValues;
	my $sthPD=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,TARGET_POS,QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID AND ID_PROTEIN=?");
	foreach my $protID (@selProteins) {
		$sthPD->execute($protID);
		while (my ($paramID,$setPos,$qValue)=$sthPD->fetchrow_array) {
			if ($quantifParamCode{$paramID} eq 'MEAN') {$quantifValues{$protID}{'MEAN'}{$setPos}=$qValue;} # WARNING: Not the same structure for MEAN
			else {$quantifValues{$protID}{$quantifParamCode{$paramID}}=$qValue;}
			$quantifValues{$protID}{SET}=$setPos if $quantifParamCode{$paramID} eq 'DELTA_PC';
		}
	}
	$sthPD->finish;

	###<Protein info>###
	my %proteinInfo;
	my $sthPI=$dbh->prepare("SELECT ALIAS,MW,PROT_DES,ORGANISM,PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=?");
	my ($GNidentID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
	my $sthGN=$dbh->prepare("SELECT MI.VALUE FROM PROTEIN P,MASTERPROT_IDENTIFIER MI WHERE P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND MI.ID_IDENTIFIER=$GNidentID AND P.ID_PROTEIN=? ORDER BY MI.RANK");
	my $sthPA=$dbh->prepare("SELECT AP.ID_ANALYSIS FROM ANALYSIS_PROTEIN AP,ANA_QUANTIFICATION AQ WHERE AP.ID_ANALYSIS=AQ.ID_ANALYSIS AND ID_PROTEIN=? AND AQ.ID_QUANTIFICATION=$selQuantifID ORDER BY VISIBILITY DESC,NUM_PEP DESC LIMIT 0,1");# To be sure that the PROTEIN is in the ANALYSIS
	foreach my $protID (@selProteins) {
		$sthPI->execute($protID);
		my @info=$sthPI->fetchrow_array;
		my @geneList=();
		$info[1]=($info[1])? sprintf "%.1f",$info[1]/1000 : '-';
		$sthGN->execute($protID);
		while (my ($gene)=$sthGN->fetchrow_array) {push @geneList,$gene;}
		$sthPA->execute($protID);
		my ($anaID)=$sthPA->fetchrow_array;
		@{$proteinInfo{$protID}}=(@info,\@geneList,$anaID); # [5] for genelist
	}
	$sthPI->finish;
	$sthGN->finish;
	$sthPA->finish;

	$dbh->disconnect;

	####<Starting HTML>####
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);

	&printProteinList(\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,'PVAL_ADJ',scalar @selProteins);

	exit;
}

sub ajaxDrawSelectedProteins {

	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;
	my $projectID=&promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION');
	#my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
	#$projectAccess=${$userInfo[2]}{$projectID};

	###<Protein info>###
	my %proteinAlias;
	my $sthPI=$dbh->prepare("SELECT ALIAS FROM PROTEIN WHERE ID_PROTEIN=?");
	foreach my $protID (split(',',param('selProt'))) {
		$sthPI->execute($protID);
		($proteinAlias{$protID})=$sthPI->fetchrow_array;
	}
	$sthPI->finish;
	my $numSelProt=scalar keys %proteinAlias;


	###<Fetching quantification design>###
	my (%stateNames,%refAnalysis);
	(my $quantifAnnot)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	my $sthgetExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
	my $sthProt=$dbh->prepare("SELECT 1 FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=? AND ID_PROTEIN=? AND VISIBILITY=2");
	my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
	my $numStates=0;
	my $featureStrg='';
	foreach my $infoStrg (@labelInfo) {
		my ($setting,$valueStrg)=split('=',$infoStrg);
		if ($setting eq 'STATES') {
			$valueStrg=~s/#//g;
			STATE:foreach my $stateData (split(';',$valueStrg)) {
				$numStates++;
				my ($numBioRep,$sampData,$expCondID)=split(',',$stateData);
				$sthgetExpCondName->execute($expCondID);
				($stateNames{$numStates})=$sthgetExpCondName->fetchrow_array;
				#<Finding reference anaID where protein is found
				if (scalar keys %refAnalysis < $numSelProt) {
					foreach my $bioReplicate (split(/\./,$sampData)) {
						foreach my $techReplicate (split(/&/,$bioReplicate)) {
							foreach my $fraction (split(/\+/,$techReplicate)) {
								my ($obsID,$parQuantifID,$anaID,$targetPos)=split(/:/,$fraction);
								foreach my $protID (keys %proteinAlias) {
									next if $refAnalysis{$protID};
									$sthProt->execute($anaID,$protID);
									my ($okProt)=$sthProt->fetchrow_array;
									if ($okProt) {
										$refAnalysis{$protID}=$anaID;
									}
								}
								next STATE if scalar keys %refAnalysis == $numSelProt;
							}
						}
					}
				}
			}
		}
		elsif ($setting eq 'PEPTIDES') {
			my ($featCode)=(split(';',$valueStrg))[3];
			$featureStrg=$featureCodes{$featCode}[0];
		}
	}
	$sthgetExpCondName->finish;
	$sthProt->finish;

	$dbh->disconnect;

	###<Reading table.txt>###
	my %proteinData;
	open(DATA,"$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/data/table.txt");
	my @colIndex2Set;
	my $numMatched=0;
	while (<DATA>) {
		chomp;
		my ($protID,@data)=split(/\t/,$_);
		if ($.==1) {
			foreach my $colLabel (@data) {
				$colLabel=~/State(\d+)\.BioRep(\d+)/;
				push @colIndex2Set,[$1,$2];
			}
			next;
		}
		next unless $proteinAlias{$protID};
		foreach my $idx (0..$#data) {
			my ($statePos,$repPos)=@{$colIndex2Set[$idx]};
			push @{$proteinData{$protID}{$statePos}},["$statePos:$repPos",$data[$idx]]; # repCode,count value
		}
		$numMatched++;
		last if $numMatched==$numSelProt;
	}
	close DATA;

	####<Starting HTML>####
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);

	my $numCols=($numSelProt <= 4)? $numSelProt : 5;
	print qq
|<TABLE border=0 cellpadding=0 cellspacing=0>
<TR><TH align="left" colspan=$numCols><INPUT type="button" value="Export raw counts" class="font11x" style="font-weight:bold" onclick="exportProteinData('xls')"/></TH></TR>
<TR>
|;
	my $colRank=0;
	my $protCount=0;
	foreach my $protID (sort{uc($proteinAlias{$a}) cmp uc($proteinAlias{$b})} keys %proteinData) {
		$colRank++;
		$protCount++;
		print qq
|<TH>&nbsp;<A href="javascript:sequenceView($protID,$refAnalysis{$protID})">$proteinAlias{$protID}</A><IMG id="pepDataDraw_$protID" src="$promsPath{images}/plus.gif" onclick="peptideDistribution(null,this,$protID)">&nbsp;<BR>
<FONT class="font11">Export as:</FONT><SELECT id="exportProtDataSEL" class="font11"><OPTION value="png">PNG</OPTION><OPTION value="svg">SVG</OPTION></SELECT><INPUT type="button" value="Go" class="font11" onclick="exportProteinData(document.getElementById('exportProtDataSEL').value,'bpDIV_$protID','$proteinAlias{$protID}')"/>
<!--
<INPUT type="button" value="Export as PNG" class="font11" onclick="exportSVGtoImg('bpDIV_$protID','$proteinAlias{$protID}','./exportSVG.cgi','png')"/>
<INPUT type="button" value="Export as SVG" class="font11" onclick="exportSVGtoImg('bpDIV_$protID','$proteinAlias{$protID}','./exportSVG.cgi','svg')"/>
-->
<DIV id="bpDIV_$protID"></DIV>
</TH>
|;
		if ($colRank == $numCols) {
			print "\n</TR>\n";
			print "<TR>\n" if $protCount < $numSelProt;
			$colRank=0;
		}
	}
	if ($protCount < $numSelProt) {
		print "<TD colspan=",($numCols-$colRank),"></TD>\n";
	}
	print qq
|</TR></TABLE>
<BR><BR>
#==========#
peptideCount={}; // reset
|;
	my %groupColor=(1=>'#0000FF',2=>'#4AA02C',3=>'#000000',4=>'#F660AB',5=>'#FBB917',6=>'#8AFB17',7=>'#9A9A9A',8=>'#7E2217',9=>'#95B9C7',10=>'#E18B6B'); # palette from genericPlot.js
	foreach my $protID (sort{uc($proteinAlias{$a}) cmp uc($proteinAlias{$b})} keys %proteinData) {
		print qq
|var BP_$protID=new boxPlot({
		div:'bpDIV_$protID',width:150,height:400,
		forceTo0:2, // ~0
		boxesTitle:' ',
		valuesTitle:'$featureStrg for $proteinAlias{$protID}',
		boxOnClick:showPeptideCount
    });
|;
		foreach my $statePos (1..$numStates) {
			my $colorRk=$statePos % 10; $colorRk=10 if $colorRk==0; # 10 colors max
			print "BP_$protID.addBox({label:'$stateNames{$statePos}',id:'$proteinAlias{$protID}:$protID:$statePos',color:'$groupColor{$colorRk}',data:[";
			my $dataStrg='';
			my $first=1;
			foreach my $refRep (@{$proteinData{$protID}{$statePos}}) {
				print ',' unless $first;
				print "'$refRep->[1]:'+replicatesLabel['$refRep->[0]']"; # 'value:label' replicatesLabel is a global JS variable
				$dataStrg.=',' unless $first;
				$dataStrg.=join('=',@{$refRep});
				$first=0;
			}
			print "]});\n";
			print "peptideCount['$proteinAlias{$protID}:$protID:$statePos']='$dataStrg';\n";
		}
		print "BP_$protID.draw();\n";
	}

	exit;
}


sub ajaxPeptideDistribution {
	my $protID=param('id_prot');

	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;
	#my $projectID=&promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION');
	my ($quantifAnnot)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	#my $sthgetExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
	my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
	my ($labelType)=($labelStrg)? ($labelStrg=~/LABEL=(.+)/) : ('FREE');
	$labelType='FREE' unless $labelType;
	my %labelingInfo;
	foreach my $infoStrg (@labelInfo) {
		my ($setting,$valueStrg)=split('=',$infoStrg);
		$valueStrg=~s/\#//g; # remove ID tag of selected modifications if any
		@{$labelingInfo{$setting}}=split(';',$valueStrg);
	}

	###<States composition>###
	my (%quantifDesign,%obs2Ana); # ,@stateNames
	my (%ana2Obs,%labelModifs,%anaLabelModifs);

	my $sthALM=$dbh->prepare("SELECT M.ID_MODIFICATION FROM ANALYSIS_MODIFICATION AM,MODIFICATION M WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND AM.ID_ANALYSIS=? AND M.IS_LABEL=1");
	my $sthObsMod=$dbh->prepare("SELECT M.ID_MODIFICATION FROM OBSERVATION O LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION WHERE O.ID_OBSERVATION=?"); # returns NULL if no modif!!!

	my $numStates=0;
	foreach my $state (@{$labelingInfo{STATES}}) { # 'state1;state2;---;stateN',  stateX='nbBioRep,bioRep1.bioRep2.---.bioRepN,#condID', bioRepX='techRep1=techRep2=---=techRepN', techRepX='frac1+frac2+---+fracN'
		$numStates++;
		my ($numBioRep,$quantiObsIDs,$expCondID)=split(/,/,$state);
		#$sthgetExpCondName->execute($expCondID);
		#($stateNames[$numStates-1])=$sthgetExpCondName->fetchrow_array;
		@{$quantifDesign{$numStates}}=();
		foreach my $bioReplicate (split(/\./,$quantiObsIDs)) {
			my @bioRep;
			my $numBioRepObs=0;
			foreach my $techReplicate (split(/&/,$bioReplicate)) {
				my @fractionObs;
				foreach my $fraction (split(/\+/,$techReplicate)) {
					$numBioRepObs++;
					my ($obsID,$parQuantifID,$anaID,$targetPos)=split(/:/,$fraction); # $parQuantifID can be 0 for SSPA
					$obs2Ana{$obsID}=$anaID;
					push @fractionObs,$obsID;
					%{$labelModifs{$obsID}}=();
					push @{$ana2Obs{$anaID}},$obsID;
					if ($labelType eq 'SILAC') { # SILAC QUANTI
						unless ($anaLabelModifs{$anaID}) { # get list of labeling modif for Analysis
							$sthALM->execute($anaID);
							while (my ($modID)=$sthALM->fetchrow_array) {$anaLabelModifs{$anaID}{$modID}=1;}
						}
						$sthObsMod->execute($obsID);
						while (my ($modID)=$sthObsMod->fetchrow_array) {
							$modID=0 unless $modID; # undef for no-label channel
							$labelModifs{$obsID}{$modID}=1; # {0} for no-label obs
							#$anaLabelModifs{$anaID}{$modID}=1 if $modID;
						}
					}
					# iTRAQ case is handled when fetching peptide quanti data (<- ??? PP 20/03/15)
					elsif ($labelType =~ /ITRAQ|TMT/) {
						my $modID;
						if ($anaLabelModifs{$anaID}) {$modID=(keys %{$anaLabelModifs{$anaID}})[0];}
						else {
							$sthALM->execute($anaID);
							($modID)=$sthALM->fetchrow_array;
next unless $modID; # !!!TEMP!!! PARAGON quantif with no validated protein (PP 30/03/15)
							$anaLabelModifs{$anaID}{$modID}=1;
						}
						$labelModifs{$obsID}{$modID}=1; # if $modID;
					}
					elsif ($labelType eq 'FREE') {$anaLabelModifs{$anaID}{0}=1; $labelModifs{$obsID}{0}=1;}
					else {die "ERROR: Unrecognized labeling ('$labelType')!";}
				} # end of fraction
				push @bioRep,\@fractionObs;
			} # end of tech replicate (sample)
			push @{$quantifDesign{$numStates}},\@bioRep;
		} # end of bio-replicate
	} # end of cond/state

	$sthALM->finish;
	$sthObsMod->finish;
	#$sthgetExpCondName->finish;

	###<Fetching and filtering peptides associated with samples>###
	# Fetch all peptide even where protein is hidden(PP 16/05/17)
	my ($pepSpecificity,$pepMissedCleavage,$ptmFilter,$pepFocus)=@{$labelingInfo{PEPTIDES}};
	#my $chargeStrg=($pepFocus=~/sp_count|_ion/)? 'CHARGE' : '0'; # sp_count|all_ion|all_pep|dist_ion|dist_pep|dist_seq (before: peptide or pepSeq)
	#my $spCountField=($pepFocus eq 'sp_count')? 'SPEC_COUNT' : 1;
	#my $peptideQuery=qq|SELECT PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),ABS(PEP_BEG),CHARGE,SPEC_COUNT
	#	FROM PEPTIDE P
	#	LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
	#	INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
	#	INNER JOIN ANALYSIS_PROTEIN AP ON PPA.ID_PROTEIN=AP.ID_PROTEIN AND P.ID_ANALYSIS=AP.ID_ANALYSIS
	#	WHERE P.ID_ANALYSIS=? AND AP.VISIBILITY=2 AND PPA.ID_PROTEIN=?|; # No need for group_concat on PEP_BEG since only 1 position is usable
	my $peptideQuery=qq|SELECT PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),ABS(PEP_BEG),CHARGE,SPEC_COUNT
		FROM PEPTIDE P
		LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
		INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
		WHERE P.ID_ANALYSIS=? AND PPA.ID_PROTEIN=?|; # No need for group_concat on PEP_BEG since only 1 position is usable
	$peptideQuery.=($pepSpecificity eq 'unique')? ' AND PPA.IS_SPECIFIC=1' : ($pepSpecificity eq 'unique_shared')? ' AND AP.PEP_SPECIFICITY=100' : ''; # Filter at DB level
	$peptideQuery.=' AND MISS_CUT=0' unless $pepMissedCleavage;
	$peptideQuery.=' GROUP BY P.ID_PEPTIDE';

	my $sthGetPep=$dbh->prepare($peptideQuery);

	$ptmFilter=~s/#//g; # remove id flag
	my ($ptmAllowed,@selectedPTMs)=split(/[:,]/,$ptmFilter); # @selectedPTMs defined only if $ptmAllowed >= 2
	my %allowedPtmID;
	foreach my $modID (@selectedPTMs) {$allowedPtmID{$modID}=1;}
	#$allowedPtmID{$quantifiedModifID}=1 if $quantifiedModifID; # to allow modif quantification

	##<Running design
	my (%anaData,%peptideData,%maxPeptideData,%excludedSeq,%peptidePos,%varModName,%varModCode2Text);
	my $refAnaID; # for sequenceView
	my $lastPepStart=0;
	foreach my $statePos (sort{$a<=>$b} keys %quantifDesign) {
		my %bioRepData;
		my $numBioRep=0;
		foreach my $refBioRep (@{$quantifDesign{$statePos}}) {
			$numBioRep++;
			my %techRepData;
			my $numTechRep=0;
			foreach my $refTechRep (@{$refBioRep}) {
				$numTechRep++;
				my %fracSpCount;
				foreach my $obsID (@{$refTechRep}) { # 1 obs = 1 fraction
					%{$fracSpCount{$obsID}}=(); # defined so no needed to checked later
					my $anaID=$obs2Ana{$obsID};
					if (!$anaData{$anaID}) {
						$sthGetPep->execute($anaID,$protID);
						PEP:while (my ($pepSeq,$varModStrg,$pepBeg,$charge,$spCount)=$sthGetPep->fetchrow_array) {
							next if $excludedSeq{$pepSeq}; # skip any future seq
							$spCount=1 unless $spCount; # not defined if ghost from match-between-run
							#next unless $spCount;
							#$varModStrg='' unless $varModStrg;
#next if ($quantifiedModifID && (!$varModStrg || $varModStrg !~ /(^|&)$quantifiedModifID:/)); # skip peptides not containing quantified modif
							#>Processing var mods & matching label channel (checking allowed & removing label if labeled quantif)
							$refAnaID=$anaID unless $refAnaID;
							my @varMods=($varModStrg)? split(/&/,$varModStrg) : ();
							my $varModCode='';
							$lastPepStart=$pepBeg if $pepBeg > $lastPepStart;
							##if ($labeling eq 'SILAC') {
							my $matchedChannel=0;
							foreach my $anaObsID (@{$ana2Obs{$anaID}}) {
								my $matchedObs=0;
								my $isLabeled=0;
								foreach my $vMod (@varMods) {
									my ($modID)=($vMod=~/^(\d+)/);
									if ($labelModifs{$anaObsID}{$modID}) { # matched labeled channel
										$matchedObs=1;
										$isLabeled=1;
										##last;
									}
									elsif ($anaLabelModifs{$anaID}{$modID}) { # un-matched labeled channel
										$isLabeled=1;
									}
								}
								if (!$isLabeled && $labelModifs{$anaObsID}{0}) { # Non-labeled channel
									$matchedObs=1;
								}
								if ($matchedObs) {
									my $hasExcludedPTM=0;
									foreach my $vMod (@varMods) {
										my ($modID)=($vMod=~/^(\d+)/);
										next if $labelModifs{$anaObsID}{$modID}; # skip labeling modifs
										if ($ptmAllowed != 1 && !$allowedPtmID{$modID}) {
											$hasExcludedPTM=1;
											last;
										}
										$varModCode.='&'.$vMod;
									}
									if ($hasExcludedPTM) { # this modification is not allowed
										push @{$excludedSeq{$pepSeq}},$protID if $ptmAllowed <= -1; # -1 or -2 unmodified peptide not allowed if exists a modified version (unmodif pep may have been recorded before detection => clean below)
										next PEP; # skip all other obs for this ana
									}
									if ($anaObsID == $obsID) {$matchedChannel=1;}
									else { # stored to be recalled when proper obsID is wanted
										push @{$anaData{$anaID}{$anaObsID}},[$pepSeq,$varModCode,$pepBeg,$charge,$spCount];
									}
									last; # obs in Ana: current anaObs was matched by peptide: no need to scan others
								}
							}

							##<Recording peptides
							if ($matchedChannel) {
								#$varModCode='' if $pepFocus eq 'dist_seq';
								$varModCode2Text{"$pepSeq:$varModCode"}=&promsMod::decodeVarMod($dbh,$pepSeq,$varModCode,\%varModName) if ($varModCode && !$varModCode2Text{"$pepSeq:$varModCode"});
								my $ionKey="$pepSeq#$varModCode#$charge";
								#if ($pepFocus eq 'sp_count') { # record best spCount/ion in each fraction
									$fracSpCount{$obsID}{$ionKey}=$spCount if (!$fracSpCount{$obsID} || !$fracSpCount{$obsID}{$ionKey} || $fracSpCount{$obsID}{$ionKey} < $spCount);
								#}
								#else { # all|distinct ion|pep|seq
								#	$techRepData{$numTechRep}{$ionKey}++;
								#}
								$peptidePos{$ionKey}=$pepBeg;
							}
						}
					}
					else { # analysis data already fetched & assigned to an Obs =>  $anaData{$anaID}{$obsID} defined
						foreach my $refPep (@{$anaData{$anaID}{$obsID}}) {
							my ($pepSeq,$varModCode,$pepBeg,$charge,$spCount)=@{$refPep};
							#$varModCode='' if $pepFocus eq 'dist_seq';
							$varModCode2Text{"$pepSeq:$varModCode"}=&promsMod::decodeVarMod($dbh,$pepSeq,$varModCode,\%varModName) if ($varModCode && !$varModCode2Text{"$pepSeq:$varModCode"});
							my $ionKey="$pepSeq#$varModCode#$charge";
							#if ($pepFocus eq 'sp_count') { # record best spCount for each ion in each fraction (cannot be directly summed because of multiple occurences of same ion with same spCount)
								$fracSpCount{$obsID}{$ionKey}=$spCount if (!$fracSpCount{$obsID} || !$fracSpCount{$obsID}{$ionKey} || $fracSpCount{$obsID}{$ionKey} < $spCount);
							#}
							#else { # all|distinct ion|pep|seq
							#	$techRepData{$numTechRep}{$ionKey}++;
							#}
							$peptidePos{$ionKey}=$pepBeg;
						}
					}
				} # End of obs/fraction

				###>Sum each ion's spCount across fractions
				#if ($pepFocus eq 'sp_count') {
					foreach my $obsID (keys %fracSpCount) {
						foreach my $ionKey (keys %{$fracSpCount{$obsID}}) {
							$techRepData{$numTechRep}{$ionKey}+=$fracSpCount{$obsID}{$ionKey};
						}
					}
				#}

				###>Clean excluded pepSeq ($ptmAllowed <= -1)
				foreach my $pepSeq (keys %excludedSeq) {
						my @peps=grep{/^$pepSeq#/} keys %{$techRepData{$numTechRep}};
						foreach my $pep (@peps) {
							delete $techRepData{$numTechRep}{$pep};
							delete $techRepData{$numTechRep} unless scalar keys %{$techRepData{$numTechRep}};
						}
					@{$excludedSeq{$pepSeq}}=(); # reset list of matching prot to 0 before next techRep
				}

			} # End of techRep

			###>Counting peptides (aggregate techReps -> mean of pepCount)
			#if ($pepFocus=~/^(sp|all)_/) { # sp_count | all_(ion|pep)
				foreach my $techRepPos (keys %techRepData) {
					while (my ($ionKey,$count) = each %{$techRepData{$techRepPos}}) {
						$bioRepData{$numBioRep}{$ionKey}+=$count;
					}
				}
				foreach my $ionKey (keys %{$bioRepData{$numBioRep}}) {
					$bioRepData{$numBioRep}{$ionKey}/=$numTechRep;
				}
			#}
			#else { # ion or peptide or pepSeq (distinct +/- charge)
			#	foreach my $techRepPos (keys %techRepData) {
			#		foreach my $ionKey (keys %{$techRepData{$techRepPos}}) {
			#			$bioRepData{$numBioRep}{$ionKey}+=1 - (($numStates-1)*0.0.5); # 1, 0.95, 0.9, 0.85, ...
			#		}
			#		foreach my $ionKey (keys %{$bioRepData{$numBioRep}}) {
			#			$bioRepData{$numBioRep}/=$numTechRep;
			#		}
			#	}
			#}
		} # End of bioRep

		foreach my $bioRepPos (keys %bioRepData) {
			while (my ($ionKey,$count) = each %{$bioRepData{$bioRepPos}}) { # (key,value)
				$peptideData{$statePos}{$ionKey}+=$count;
			}
		}
		foreach my $ionKey (keys %{$peptideData{$statePos}}) {
			$peptideData{$statePos}{$ionKey}/=$numBioRep;
			$maxPeptideData{$ionKey}=$peptideData{$statePos}{$ionKey} if (!defined $maxPeptideData{$ionKey} || $maxPeptideData{$ionKey} < $peptideData{$statePos}{$ionKey});
		}

	} # End of State
	$sthGetPep->finish;

	my $sthPI=$dbh->prepare("SELECT ALIAS,PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=?");
	$sthPI->execute($protID);
	my ($proteinAlias,$protLength)=$sthPI->fetchrow_array;
	$sthPI->finish;

	$dbh->disconnect;

	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);
	my $plotWidth=($protLength)? $protLength*1.5 : $lastPepStart+50;
	$plotWidth=400 if $plotWidth < 400;
	$plotWidth=1000 if $plotWidth > 1000;
	my @sortedStatePos=sort{$a<=>$b} keys %peptideData;
	my $plotHeight=70*(1+$sortedStatePos[-1]-$sortedStatePos[0]);
	#my $JSstateStrg="'".join("','",@stateNames)."'";
	print qq
|<TABLE><TR>
<TD nowrap><FONT class="title2">Peptide distribution for <A href="javascript:sequenceView($protID,'$refAnaID')">$proteinAlias</A> in all States:</FONT></TD>
<TD nowrap>&nbsp;&nbsp;<FONT class="title3">(View:</FONT><SELECT class="title3" onchange="for (let i=1; i<=2; i++) {var pepDiv=document.getElementById('pepPlot'+i+'DIV'); pepDiv.style.display=(pepDiv.style.display=='none')? '' : 'none';}"><OPTION value="1">Normalized counts</OPTION><OPTION value="2">Absolute counts</OPTION></SELECT><FONT class="title2">)</FONT></TD>
<TD nowrap>&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'; document.getElementById(prevImgID).src='$promsPath{images}/plus.gif'; prevImgID=null;"></TD>
</TR></TABLE>
<DIV id="pepPlot1DIV"></DIV>
<DIV id="pepPlot2DIV"></DIV>
#==========#
var PP1=new peptidePlot({div:'pepPlot1DIV',width:$plotWidth,height:$plotHeight,valueAxisLabel:'State position',valueLabel:'State',
	protein:{length:$protLength},
	convertValueDisplayed: function(v) {return statesName[Math.floor(v)]}, //state pos -> state name
	peptideProperties:['Mean # features','Norm. mean # features (%)'],
	peptideColor:{map:'Norm. mean # features (%)',type:'continuous',range:[0,100]}
});
|;
	foreach my $statePos (@sortedStatePos) {
		my $stateIdx=$statePos-1;
		foreach my $ionKey (sort{$peptidePos{$a}<=>$peptidePos{$b}} keys %{$peptideData{$statePos}}) {
			my ($pepSeq,$varModCode,$charge)=split('#',$ionKey);
			my $modStrg=($varModCode)? "'".$varModCode2Text{"$pepSeq:$varModCode"}."'" : 'null';
			#print "PP.addPeptide([$peptidePos{$ionKey},'$pepSeq',$modStrg,$charge,$peptideData{$statePos}{$ionKey},0,'$stateNames[$stateIdx]']);\n";
			print "PP1.addPeptide([$peptidePos{$ionKey},'$pepSeq',$modStrg,$charge,",$statePos+(0.1*($charge-2)),",0,",1*(sprintf "%.3f",$peptideData{$statePos}{$ionKey}),",",1*(sprintf "%.3f",100*$peptideData{$statePos}{$ionKey}/$maxPeptideData{$ionKey}),"]);\n";
		}
	}
	print "PP1.draw();\n";
	#---2nd view (absolute counts)
	print qq
|var PP2=new peptidePlot({div:'pepPlot2DIV',width:$plotWidth,height:$plotHeight,valueAxisLabel:'State position',valueLabel:'State',
	protein:{length:$protLength},
	convertValueDisplayed: function(v) {return statesName[Math.floor(v)]}, //state pos -> state name
	peptideProperties:['Mean # features','Norm. mean # features (%)'],
	peptideColor:{map:'Mean # features',type:'continuous'}
});
|;
	foreach my $statePos (@sortedStatePos) {
		my $stateIdx=$statePos-1;
		foreach my $ionKey (sort{$peptidePos{$a}<=>$peptidePos{$b}} keys %{$peptideData{$statePos}}) {
			my ($pepSeq,$varModCode,$charge)=split('#',$ionKey);
			my $modStrg=($varModCode)? "'".$varModCode2Text{"$pepSeq:$varModCode"}."'" : 'null';
			print "PP2.addPeptide([$peptidePos{$ionKey},'$pepSeq',$modStrg,$charge,",$statePos+(0.1*($charge-2)),",0,",1*(sprintf "%.3f",$peptideData{$statePos}{$ionKey}),",",1*(sprintf "%.3f",100*$peptideData{$statePos}{$ionKey}/$maxPeptideData{$ionKey}),"]);\n";
		}
	}
	print "PP2.draw();\ndocument.getElementById('pepPlot2DIV').style.display='none';\n";
	exit;
}


# TODO: Export to Excel
####>Revision history<####
# 1.0.5 [FEATURE] Compatible with peptide xic feature (PP 28/02/20)
# 1.0.4 [FEATURE] Compatible with MBR-peptide filter (PP 22/01/20)
# 1.0.3 [FEATURE] Added extended search option for graphical view (PP 19/12/19)
# 1.0.2 [FEATURE] Compatible with option to extend missed cleavage exclusion to overlapping peptides (PP 19/11/19)
# 1.0.1 Multiple bug fixes (PP 19/01/18)
# 1.0.0 Tested with SILAC & label-free quantif (PP 01/12/17)
# 0.9.3 Working beta version. TODO: Export to Excel (PP 16/05/17)
# 0.0.1 Created (PP 04/08/16)