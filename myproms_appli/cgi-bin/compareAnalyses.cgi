#!/usr/local/bin/perl -w

################################################################################
# compareAnalyses.cgi               2.3.1                                      #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Compares the protein & peptide contents of multiple (groups of) analyses     #
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
use Encode 'decode_utf8';
#use utf8; # needed only if utf8 char from script itself must be printed to Excel
use POSIX qw(strftime); # to get the time
use Spreadsheet::WriteExcel;
use promsConfig;
use promsMod;
use phosphoRS;
use strict;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my $MAX_NUM_GROUPS=50;
my ($color1,$color2)=&promsConfig::getRowColors;
my $CONTEXT_SIZE=15; # size of sequence window used for modif site export

####################
####>Parameters<####
####################
if (param('AJAX')) {
	if (param('AJAX') eq 'getAnaList') {&ajaxGetAnalysesList;}
	elsif(param('AJAX') eq 'ajaxRestrictList') {&ajaxRestrictProteinList;}
	elsif(param('AJAX') eq 'getAnaText') {&ajaxUpdateCompGroupAnalyses;}
	elsif (param('AJAX') eq 'getProtDetails') {&ajaxGetProteinDetails;}
	elsif (param('AJAX') eq 'getSiteDetails') {&ajaxGetModificationSiteDetails;}
	exit;
}
#print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG;
my $projectID=param('id_project');
my $parentItem=lc(param('parentItem'));
my $sortOrder=(param('sort'))? param('sort') : 'protein';
my ($item,$itemID)=split(/:/,$parentItem);
my $comparisonID=(param('comparisonID'))? param('comparisonID') : (param('targetComp'))? param('targetComp') : 0;
my $action=(param('ACT'))? param('ACT') : 'frames';
my $compType=param('compType') || 'full_prot';

####>Master Frame<####
if ($action eq 'frames') {&generateFrames; exit;}

my $numGroups=param('numGroups') || 0;
my $groupPosStrg=param('groupPos') || ''; # only if full_prot called from 1v1_prot
my $anaList=param('anaList');
my $virtualData=(param('virtualData'))? 1 : 0; # XIC-based identifications (not used by 1v1_prot)
my $pepSpecificity=param('pepSpecificity') || 'all'; # unique,unique_shared,all
#>Protein-comparison parameters (defined even if peptide comp)----->
my $peptideRule=param('pepRule') || 'all'; # all,nr,ba,banr (not used by 1v1_pep)
my $pepThreshold=param('pepThreshold') || 1; # 1 to 5 (not used by 1v1_pep)
my $catFilter=param('catFilter') || 0; # (not used by 1v1_pep)
my $catFilterRule=param('catFilterRule') || 'restrict'; # (not used by 1v1_pep)
my $autochkStrg=param('autocheck') || ''; # (used only by full_prot)
# Visibility of unchecked proteins
my $hideShowStatus=(param('hideShowStatus'))? param('hideShowStatus') : ''; # (used only by full_prot)
my $hideShowStrg=($hideShowStatus eq '')? 'Hide Unchecked' : 'Show Unchecked'; # (used only by full_prot)
#>Peptide-comparison parameters (defined even if protein comparison)--->
my $noMissCut=param('noMissCut') || 0;
my $delocPhospho=param('delocPhospho') || 0;
my @modifList=(param('modifList'))? split(',',param('modifList')) : (); # not defined before submit
my $restrictModifLogic=param('restrictModifLogic') || 'all';
my %modifFilters=(
	exclude=>{},
	restrict=>{},
	ignore=>{},
	active=>0
);
foreach my $modID (@modifList) {
	my $filter=param("filterMod_$modID");
	next unless $filter; # param("filterMod_$modID") not always defined. eg action='savComp'
	if ($filter ne 'allow') { # no need to record 'allow'
		$modifFilters{$filter}{$modID}=1;
		$modifFilters{active}=1;
	}
}
#>Modif site-specific parameters (defined for all compare modes)--->
my $targetedRes=param('siteRes') || '';
my $siteMeasure=param('siteMeas') || 'occ'; # specAbs, specRel

####>Analysis Selection Frame<####
if ($action eq 'selAna') {&selectAnalyses; exit;}
my $call=param('call') || 'selAna'; # only if full_prot called from 1v1_prot

####>Comparison Result Frame<####
my (@analysisGroups,@groupPos);
if ($numGroups) { # real groups
	foreach my $grDataStrg (split(',',$anaList)) {
		push @analysisGroups,[split(/\./,$grDataStrg)]; # :
	}
	if ($groupPosStrg) { # only if full_prot called from 1v1_prot
		foreach my $pos (split(',',$groupPosStrg)) {
			push @groupPos,$pos;
		}
	}
}
else { # fake groups (1 ana / gr)
	foreach my $codeID (split(/\./,$anaList)) { # :
		push @analysisGroups,[$codeID];
	}
}
my $numAnaGroups=scalar @analysisGroups;
my $listComparison=1 if $anaList=~/C:/; # list comparison

#>Checking sortOrder compatibility
if ($sortOrder=~/group:(\d+)/) {$sortOrder="group:$numAnaGroups" if $1>$numAnaGroups;}
elsif ($sortOrder=~/Delta/ && $numAnaGroups != 2) {$sortOrder='protein';}
my $noHidden=(param('noHidden'))? 1 : 0; # (not used by 1v1_pep)
my $hiddenRule=($noHidden)? param('hiddenRule') : 0; # show all hidden, 1=>skip if hidden everywhere, 2=>skip any hidden (not used by 1v1_pep)

####>Saving a comparison<####
if ($action eq 'saveComp') {&storeComparison; exit;}

my $exportType=param('exportType');
#my $checkBoxName='check_'.uc($item).'_'.$itemID;

####>Connecting to the database<####
my $dbh=&promsConfig::dbConnect;
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};

###>Project PTMs & status
my %allPostTransModifs=&promsMod::getVariableModifications($dbh);
#my ($projectPtmString,$projectStatus)=$dbh->selectrow_array("SELECT RELEVANT_PTMS,STATUS FROM PROJECT WHERE ID_PROJECT=$projectID");
my ($projectStatus)=$dbh->selectrow_array("SELECT STATUS FROM PROJECT WHERE ID_PROJECT=$projectID");
$projectStatus=0 unless $projectStatus;

my $sthGetPM=$dbh->prepare("SELECT ID_MODIFICATION FROM PROJECT_MODIFICATION WHERE ID_PROJECT=$projectID");
$sthGetPM->execute;
my %projectVarMods;
while (my ($modID)=$sthGetPM->fetchrow_array) {
	$projectVarMods{$modID}=$allPostTransModifs{$modID};
}
$sthGetPM->finish;


my ($compName,$compComments);
if ($comparisonID) {
	($compName,$compComments)=$dbh->selectrow_array("SELECT NAME,COMMENTS FROM COMPARISON WHERE ID_COMPARISON=$comparisonID");
}
$compName='' unless $compName;
$compComments='' unless $compComments;
#$compComments=quotemeta($compComments);

############################
####>Peptide comparison<####
############################
if ($compType eq '1v1_pep') {
	&comparePeptides;
	exit;
}

####<Selected proteins/sites>####
my %selectedProteins; # proteins checked (Export selected or change in list display)
if (param('chkProt')) {
	foreach my $protInfo (param('chkProt')) {
		my $protID=(split(':',$protInfo))[-1];
		$protID=~s/-.+// unless $compType =~ /^modif_sites/; # in case coming for modif site
		$selectedProteins{$protID}=1;
#print "$protInfo '$protID'<BR>\n";
	}
}

#########################
####>Site comparison<####
#########################
if ($compType =~ /^modif_sites/) {
	$sortOrder='peptide' if $sortOrder=~/Delta/;
	&compareModificationSites;
	exit;
}

############################
####>Protein comparison<####
############################

###############
#     Main    #
###############

###################################
####>Protein comparison Header<####
###################################
if ($compType eq '1v1_prot') {
	#&compareProteins;
	#exit;
	&display1v1ComparisonHead;
}
else {
	if ($exportType) {
		my $timeStamp1=strftime("%Y%m%d %H-%M",localtime);
		print header(-type=>"application/vnd.ms-excel",-attachment=>"Compare items_$timeStamp1.xls");
	}
	else {
		print header(-'content-encoding'=>'no',-charset=>'utf-8');
		warningsToBrowser(1);
		print qq
|<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TABLE.summary TH {font-size:15px;}
/*.item {font-size:17px;color:#DD0000;}*/
.small {font-size:10px;}
.TD {font-weight:normal;}
.TH {font-weight:bold;}
ACRONYM {cursor:help;}
.popPTM {font-weight:bold;background-color:$color1;}
.badPTM {color:#505050;}
.row_0 {background-color:$color2;}
.row_1 {background-color:$color1;}
|;
		foreach my $varMod (sort{$a cmp $b} keys %projectVarMods) {
			print ".$projectVarMods{$varMod}[2] \{$projectVarMods{$varMod}[3]\}\n";
		}
		print qq
|/* PTMs */
.visibleLinkElement {visibility:visible;cursor:pointer;}
</STYLE>
|;
		if ($numAnaGroups<=5) {
			print qq
|<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/vennDiagram.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
|;
		}
		print qq
|<SCRIPT type="text/javascript">
function updateLists(themeSel,listSel,XHR) {
	parent.selItemsFrame.ajaxUpdateRestrict();
	if (themeSel.value=='getLists:-1') { // update anaParent in selItemsFrame
		var idData=XHR.responseText.match(/THEME_ID=(\\d+)/);
		var nameData=XHR.responseText.match(/THEME_NAME=(.+)###/);
		var selParent=parent.selItemsFrame.compForm.anaParent;
		var idx=(selParent.options[selParent.length-1].text=='None')? selParent.length-1 : selParent.length;
		selParent[idx]=new Option(nameData[1],'CLASSIFICATION:'+idData[1]);
	}
}
|;
		##>Save proteins to custom Lists
		&promsMod::printAjaxManageSaveProteins($projectID,\%promsPath,'document.protView.chkProt','updateLists');

		&promsMod::popupInfo();

		print qq
|</SCRIPT>
</HEAD>
<BODY style="margin:2px;" background="$promsPath{images}/bgProMS.gif">
<CENTER>
<DIV id="waitDIV">
<BR><BR><BR><BR><BR><FONT class="title3">Fetching data...</FONT><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">
<BR><BR><BR><FONT class="title3">Status: <SPAN id="waitSPAN"><SPAN>...</FONT>
</DIV>
</CENTER>
|;
	}
} # end of full comparison header

####>Looping through analyses<####
###<Queries
my @sthAnaList;
#if ($item eq 'project') {
my $projQuery1=qq |SELECT E.NAME,G.NAME,SP.NAME,A.NAME FROM EXPERIMENT E,GEL2D G,SPOT SP,SAMPLE S,ANALYSIS A
						WHERE A.ID_ANALYSIS=? AND S.ID_SAMPLE=A.ID_SAMPLE AND SP.ID_SPOT=S.ID_SPOT AND G.ID_GEL2D=SP.ID_GEL2D AND E.ID_EXPERIMENT=G.ID_EXPERIMENT|;
my $projQuery2=qq |SELECT E.NAME,S.NAME,A.NAME FROM EXPERIMENT E,SAMPLE S,ANALYSIS A
						WHERE A.ID_ANALYSIS=? AND S.ID_SAMPLE=A.ID_SAMPLE AND E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SPOT IS NULL|;
push @sthAnaList,($dbh->prepare($projQuery1),$dbh->prepare($projQuery2)); # free sample then gel

my $sthCatInfo=$dbh->prepare("SELECT CLASSIFICATION.NAME,CATEGORY.NAME FROM CLASSIFICATION,CATEGORY WHERE CLASSIFICATION.ID_CLASSIFICATION=CATEGORY.ID_CLASSIFICATION AND ID_CATEGORY=?");
my $sthCatProt=$dbh->prepare("SELECT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=?");
#}
#elsif ($item eq 'experiment') {
#	push @sthAnaList,$dbh->prepare("SELECT SAMPLE.NAME,ANALYSIS.NAME FROM SAMPLE,ANALYSIS WHERE ID_ANALYSIS=? AND SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND SAMPLE.ID_SPOT IS NULL"); # free sample
#	push @sthAnaList,$dbh->prepare("SELECT GEL2D.NAME,SPOT.NAME,ANALYSIS.NAME FROM GEL2D,SPOT,SAMPLE,ANALYSIS WHERE ID_ANALYSIS=? AND SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND SPOT.ID_SPOT=SAMPLE.ID_SPOT AND GEL2D.ID_GEL2D=SPOT.ID_GEL2D"); # gel
#}
#elsif ($item eq 'gel2d') {
#	push @sthAnaList,$dbh->prepare("SELECT SPOT.NAME,ANALYSIS.NAME FROM SPOT,SAMPLE,ANALYSIS WHERE ID_ANALYSIS=? AND SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND SPOT.ID_SPOT=SAMPLE.ID_SPOT");
#}
#elsif ($item eq 'sample' || $item eq 'spot') {
#	push @sthAnaList,$dbh->prepare("SELECT NAME FROM ANALYSIS WHERE ID_ANALYSIS=?");
#}
my $visibilityStrg=($hiddenRule==2)? 'AND VISIBILITY>=1' : '';
my ($confLevelStrg,$pepValidStrg)=($virtualData)? ('','') : ('AND CONF_LEVEL > 0','AND P.VALID_STATUS > 0'); # virtual proteins created by MassChroQ alignment
my $protPepSpecifStrg=($pepSpecificity=~/unique/)? 'AND PEP_SPECIFICITY=100' : '';
my $isSpecificStrg=($pepSpecificity eq 'unique')? 'AND IS_SPECIFIC=1' : '';
###my $sthAP=$dbh->prepare("SELECT ID_PROTEIN,VISIBILITY,CONF_LEVEL,SCORE FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=? $confLevelStrg $visibilityStrg"); # skip virtual proteins
my $sthAP=$dbh->prepare("SELECT AP.ID_PROTEIN,VISIBILITY,CONF_LEVEL,AP.SCORE,P.ID_PEPTIDE,PEP_SEQ FROM ANALYSIS_PROTEIN AP,PEPTIDE_PROTEIN_ATTRIB PPA,PEPTIDE P WHERE AP.ID_ANALYSIS=? AND PPA.ID_ANALYSIS=? AND AP.ID_PROTEIN=PPA.ID_PROTEIN AND PPA.ID_PEPTIDE=P.ID_PEPTIDE $confLevelStrg $visibilityStrg $pepValidStrg $protPepSpecifStrg $isSpecificStrg ORDER BY AP.ID_PROTEIN");
my $sthPep=$dbh->prepare("SELECT P.ID_PEPTIDE,PEP_SEQ FROM PEPTIDE_PROTEIN_ATTRIB PPA,PEPTIDE P WHERE ID_PROTEIN=? AND PPA.ID_ANALYSIS=? AND PPA.ID_PEPTIDE=P.ID_PEPTIDE $pepValidStrg $isSpecificStrg"); # skip ghosts
my $sthPepMod=$dbh->prepare("SELECT PM.ID_MODIFICATION,P.ID_PEPTIDE,POS_STRING FROM PEPTIDE_MODIFICATION PM,PEPTIDE P WHERE PM.ID_PEPTIDE=P.ID_PEPTIDE AND P.ID_ANALYSIS=? $pepValidStrg"); # <- no pepSpecificity filter applied (not needed)
my (@groupLabels,@groupPeptides,%proteinList,%bestVisibility,%bestAnalysis,%protAnalyses,%processedProtAna,%peptideAnaDistrib,%ptmAnalyses,%catProteins,%anaPepVmods);

my (@parentItemAna,%groupListCompProt,%listInGroup,$contextAna,@groupListAna);
my $listContextStrg='';
if ($listComparison) { # fetch list of analyses in selected parent item
	my ($contextName)=$dbh->selectrow_array("SELECT NAME FROM ".uc($item)." WHERE ID_".uc($item)."=$itemID");
	$listContextStrg='List context is '.&promsMod::getItemType($item)." '$contextName'";
	$listContextStrg="&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<FONT class=\"title3\" color=\"#DD0000\">$listContextStrg</FONT>" unless $exportType;

	###>Analyses in selected parent item<###
	my @queryList;
	my $ITEM=uc($item);
	if ($ITEM eq 'PROJECT') {
		#push @queryList,"SELECT DISTINCT ID_ANALYSIS FROM PROTEIN P,ANALYSIS_PROTEIN AP WHERE P.ID_PROTEIN=AP.ID_PROTEIN AND ID_PROJECT=$itemID";
		my $projQS=qq |SELECT A.ID_ANALYSIS
								FROM EXPERIMENT E,SAMPLE S,ANALYSIS A
								WHERE E.ID_PROJECT=$itemID AND E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SAMPLE=A.ID_SAMPLE AND S.ID_SPOT IS NULL AND A.VALID_STATUS>=1 ORDER BY E.DISPLAY_POS,S.DISPLAY_POS,A.DISPLAY_POS|;
		my $projQG=qq |SELECT A.ID_ANALYSIS
								FROM EXPERIMENT E,GEL2D G,SPOT SP,SAMPLE S,ANALYSIS A
								WHERE E.ID_PROJECT=$itemID AND E.ID_EXPERIMENT=G.ID_EXPERIMENT AND G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT
								AND S.ID_SAMPLE=A.ID_SAMPLE AND A.VALID_STATUS>=1 ORDER BY E.DISPLAY_POS,G.DISPLAY_POS,SP.NAME,A.DISPLAY_POS|;
		push @queryList,($projQS,$projQG);
	}
	elsif ($ITEM eq 'EXPERIMENT') {
		my $expQS=qq |SELECT A.ID_ANALYSIS
								FROM SAMPLE S,ANALYSIS A
								WHERE S.ID_EXPERIMENT=$itemID AND S.ID_SAMPLE=A.ID_SAMPLE AND S.ID_SPOT IS NULL AND A.VALID_STATUS>=1 ORDER BY S.DISPLAY_POS,A.DISPLAY_POS|;
		my $expQG=qq |SELECT A.ID_ANALYSIS
								FROM GEL2D G,SPOT SP,SAMPLE S,ANALYSIS A
								WHERE G.ID_EXPERIMENT=$itemID AND G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT
								AND S.ID_SAMPLE=A.ID_SAMPLE AND A.VALID_STATUS>=1 ORDER BY G.DISPLAY_POS,SP.NAME,A.DISPLAY_POS|;
		push @queryList,($expQS,$expQG);
	}
	elsif ($ITEM eq 'GEL2D') {
		my $gelQuery=qq |SELECT A.ID_ANALYSIS
								FROM SPOT SP,SAMPLE S,ANALYSIS A
								WHERE SP.ID_GEL2D=$itemID AND SP.ID_SPOT=S.ID_SPOT
								AND S.ID_SAMPLE=A.ID_SAMPLE AND A.VALID_STATUS>=1 ORDER BY SP.NAME,A.DISPLAY_POS|;
		push @queryList,$gelQuery;
	}
	#elsif ($parentItem eq 'SPOT') {
	#	push @queryList,"SELECT A.ID_ANALYSIS FROM SAMPLE S,ANALYSIS A WHERE S.ID_SPOT=$parentID AND S.ID_SAMPLE=A.ID_SAMPLE AND A.VALID_STATUS>=1";
	#}
	elsif ($ITEM eq 'SAMPLE') {
		push @queryList,"SELECT ID_ANALYSIS FROM ANALYSIS WHERE ID_SAMPLE=$itemID AND VALID_STATUS>=1 ORDER BY DISPLAY_POS";
	}
	elsif ($ITEM eq 'ANALYSIS') {
		push @queryList,"SELECT $itemID";
	}
	foreach my $query (@queryList) {
		my $sthAna=$dbh->prepare($query);
		$sthAna->execute;
		while (my ($anaID)=$sthAna->fetchrow_array) {
			push @parentItemAna,$anaID;
		}
		$sthAna->finish;
	}

	###>Proteins in List comparison<###
	foreach my $g (0..$#analysisGroups) {
		foreach my $codeID (@{$analysisGroups[$g]}) {
			my ($type,$typeID)=split(':',$codeID); # A or C
			next if $type eq 'A';
			foreach my $anaID (@parentItemAna) {$groupListAna[$g]{$anaID}=1;} # for process progress info
			$listInGroup{$g}++;
			%{$catProteins{$typeID}}=();
			$sthCatProt->execute($typeID);
			while (my ($protID)=$sthCatProt->fetchrow_array) {
				$catProteins{$typeID}{$protID}=1;
				$groupListCompProt{$g}{$protID}=1;
			}
		}
	}
}
$contextAna=join(',',@parentItemAna);

###>Category filter<###
my ($filterClassName,$filterCatName);
if ($catFilter) {
	#my $sthCP=$dbh->prepare("SELECT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$catFilter");
	%{$catProteins{$catFilter}}=();
	$sthCatProt->execute($catFilter);
	while (my ($protID)=$sthCatProt->fetchrow_array) {
		$catProteins{$catFilter}{$protID}=1;
	}
	#$sthCP->finish;
	#($filterClassName,$filterCatName)=$dbh->selectrow_array("SELECT CLASSIFICATION.NAME,CATEGORY.NAME FROM CLASSIFICATION,CATEGORY WHERE CLASSIFICATION.ID_CLASSIFICATION=CATEGORY.ID_CLASSIFICATION AND ID_CATEGORY=$catFilter");
	$sthCatInfo->execute($catFilter);
	($filterClassName,$filterCatName)=$sthCatInfo->fetchrow_array;
}
$sthCatProt->finish;


###>Primary scan<###
my $newLine=($exportType)? "\n" : '<BR>';
my $maxItemInGr=1;
foreach my $g (0..$#analysisGroups) {
#$trueAnaList.=',' if $trueAnaList;
	my $grStrg=($numGroups)? 'Group '.($g+1).'/'.$numGroups.', ' : '';
	my $numItemInGr=scalar @{$analysisGroups[$g]};
	$maxItemInGr=$numItemInGr if $numItemInGr > $maxItemInGr;
	my $numAnaInGr=($groupListAna[$g])? scalar keys %{$groupListAna[$g]} : $numItemInGr;
	my (%protPeptides,%groupProteins,%bestScore,%groupBestVis);
	my %itemAnaScanIdx; # if C only scan parent item analyses once/protein
	my $codeIndex=-1;
	my $numAna=0;
	foreach my $codeID (@{$analysisGroups[$g]}) {
		$codeIndex++;
		my ($type,$typeID)=split(':',$codeID); # A or C
#$trueAnaList.=':' if $codeIndex;
		my (@typeInfo,@localAnaList);
		if ($type eq 'A') { # Analysis
			foreach my $sthAna (@sthAnaList) {
				$sthAna->execute($typeID);
				@typeInfo=$sthAna->fetchrow_array;
				last if $typeInfo[0]; # 2 loops in case of project or experiment
			}
			@localAnaList=($typeID);
		}
		else { # C -> List
			$grStrg.='[List] ';
			$sthCatInfo->execute($typeID); # catID in fact!
			@typeInfo=$sthCatInfo->fetchrow_array;
			@localAnaList=@parentItemAna;
		}
		if ($numGroups) {
			unless ($groupLabels[$g][0]) {$groupLabels[$g][0]=($groupPos[$g])? 'Group #'.$groupPos[$g] : 'Group #'.($g+1);}
			$groupLabels[$g][1].=$newLine if $groupLabels[$g][1];
			$groupLabels[$g][1].='-'.join(' > ',@typeInfo);
		}
		else {
			@{$groupLabels[$g]}=(&promsMod::shortenName($typeInfo[-1],20),join(' > ',@typeInfo));
		}
		if ($type eq 'C') { #  && defined($itemAnaScanIdx)
			foreach my $protID (keys %itemAnaScanIdx) { # ana list already scanned for protein  // %peptideAnaDistrib
				next unless $catProteins{$typeID}{$protID}; # protein is not in List
				$peptideAnaDistrib{$protID}[$g][$codeIndex]=$peptideAnaDistrib{$protID}[$g][$itemAnaScanIdx{$protID}]; # replicate ana pattern: copy reference NOT data @{}!
			}
			#next;
		}
#$trueAnaList.=join(':',@localAnaList);
		my $anaIndex=-1; # index of Analyses in group (not A/C items)
		foreach my $anaID (@localAnaList) {
			$anaIndex++;
			unless ($exportType) {
				if ($numGroups) {
					$numAna++;
					&updateProgress('Processing '.$grStrg."Analysis $numAna/$numAnaInGr..."); # $numAna can be > $numAnaInGr if >1 lists in gr!
				}
				else {
					&updateProgress('Processing Analysis '.($g+1)."/$numAnaGroups..."); # $numAna can be > $numAnaInGr if >1 lists in gr!
				}
			}

			###> PTMs
			unless ($anaPepVmods{$anaID}) {
				%{$anaPepVmods{$anaID}}=();
				##$sthPepMod->execute($anaID,$anaID);
				$sthPepMod->execute($anaID);
				while (my ($modID,$pepID,$posString) = $sthPepMod->fetchrow_array) {
					$anaPepVmods{$anaID}{$pepID}{$modID}=$posString;
				}
			}
			$sthAP->execute($anaID,$anaID);
			my $prevProtID=0;
			my ($prevVis,$prevConfLevel,$prevScore);
			my (%distinctPepAna,%allPepAna);
			while (my ($protID,$vis,$confLevel,$score,$pepID,$pepSeq)=$sthAP->fetchrow_array) { # No longer uses ANALYSIS_PROTEIN.NUM_PEP (PP 11/08/15)
				next if ($exportType && $exportType eq 'sel' && !$selectedProteins{$protID});
				next if ($catFilter && ($catFilterRule eq 'restrict' && !$catProteins{$catFilter}{$protID}) || ($catFilterRule eq 'exclude' && $catProteins{$catFilter}{$protID}));
				#next if ($type eq 'C' && !$groupListCompProt{$g}{$protID}); # List comparison => skip if prot is not in (group of) list(s)
				next if ($type eq 'C' && !$catProteins{$typeID}{$protID}); # List comparison => skip if prot is not in List /// (group of) list(s)   $groupListCompProt{$g}{$protID}
				if ($protID != $prevProtID) {
					$processedProtAna{$protID}{$anaID}=0;
					@{$protAnalyses{$protID}{$anaID}}=(-1,$vis,$confLevel,$score); # $protAnalyses{$protID}{$anaID}[0] tagged in case skipped at this stage
					#$protAnalyses{$protID}{$anaID}[4]=($noHidden)? $vis : 1; # used for checkbox value for graphicalView & to record hidden proteins if hiddenRule=1
				}
				next if ($vis==0 && $hiddenRule); # hiddenRule=1 (cannot be 2 if vis=0)
				$groupProteins{$protID}=1;

				##>True numPep because of +/-distinct or +/-virtual peptides
				if ($prevProtID && $protID != $prevProtID) {
					my $numPep=($peptideRule=~/nr/)? scalar keys %distinctPepAna : scalar keys %allPepAna;
					$protAnalyses{$prevProtID}{$anaID}[0]=$numPep; # overwrites $protAnalyses{$protID}{$anaID}}[0]
					&updateProteinData($g,$codeIndex,$type,$anaID,$anaIndex,$prevProtID,$prevVis,$prevConfLevel,$prevScore,$numPep,\%bestScore,\%groupBestVis,\%itemAnaScanIdx);
					%distinctPepAna=%allPepAna=();
				}

				###>Compute true $numPep (PP 11/08/15)
				$allPepAna{$pepID}=1;
				my ($varModStrg,@modifications)=('',());
				if ($anaPepVmods{$anaID}{$pepID}) {
					foreach my $modID (sort{$a<=>$b} keys %{$anaPepVmods{$anaID}{$pepID}}) {
						$varModStrg.="+$modID:$anaPepVmods{$anaID}{$pepID}{$modID}";
						if ($projectVarMods{$modID}) {
							$ptmAnalyses{$protID}{$g}{$anaID}{$modID}{'ALL'}++;
							push @modifications,$modID;
						}
					}
					if ($peptideRule=~/nr/) {# prevents INTRA-analysis redundancy
						foreach my $modID (@modifications) {
							$ptmAnalyses{$protID}{$g}{$anaID}{$modID}{'NR'}{"$pepSeq$varModStrg"}=1;
						}
					}
				}
				if ($peptideRule=~/nr/) {
					$protPeptides{$protID}{"$pepSeq$varModStrg"}=1;
					$distinctPepAna{"$pepSeq$varModStrg"}=1;
				}

				($prevProtID,$prevVis,$prevConfLevel,$prevScore)=($protID,$vis,$confLevel,$score);

			}
			##>True numPep because of +/-distinct or +/-virtual peptides (Again for last protein in loop)
			if ($prevProtID) {
				my $numPep=($peptideRule=~/nr/)? scalar keys %distinctPepAna : scalar keys %allPepAna;
				$protAnalyses{$prevProtID}{$anaID}[0]=$numPep; # overwrites $protAnalyses{$protID}{$anaID}}[0]
				&updateProteinData($g,$codeIndex,$type,$anaID,$anaIndex,$prevProtID,$prevVis,$prevConfLevel,$prevScore,$numPep,\%bestScore,\%groupBestVis,\%itemAnaScanIdx);
				#%distinctPepAna=%allPepAna=();
			}
		} # end of localAnaList
	#$itemAnaScanIdx=$codeIndex if $type eq 'C';
	}
	##>Computing peptides distribution in groups & global best visibility
	foreach my $protID (keys %groupProteins) {
		#>Correcting numPep in %proteinList for nr cases
		#$proteinList{$protID}[$g][0]=scalar keys %{$protPeptides{$protID}} if $peptideRule=~/nr/; # Correcting numPep in %proteinList
		if ($peptideRule eq 'anr') {
			$proteinList{$protID}[$g][0]=scalar keys %{$protPeptides{$protID}};
		}
		elsif ($peptideRule eq 'banr') {
			$proteinList{$protID}[$g][0]=$protAnalyses{$protID}{$groupBestVis{$protID}[0]}[0];
		}
		$groupPeptides[$g]{$proteinList{$protID}[$g][0]}=1; # to select num peptides for auto-check

		@{$bestVisibility{$protID}}=@{$groupBestVis{$protID}} if (!$bestVisibility{$protID} || $bestVisibility{$protID}[1]<$groupBestVis{$protID}[1]);
	}

}

foreach my $sthAna (@sthAnaList) {$sthAna->finish;}
$sthAP->finish;
$sthCatInfo->finish;

###>Secondary scan: Filling blanks<###
&updateProgress('Checking protein distribution...') unless $exportType;
foreach my $protID (keys %proteinList) {
	foreach my $g (0..$#analysisGroups) {
		#next if defined($proteinList{$protID}[$g]); #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if ($hiddenRule) {
			##>Scanning group's analyses for hidden protein
			my $predefinedProtGr=(defined($proteinList{$protID}[$g]))? 1 : 0;
			my (%protPeptides,%bestScore);
			my $codeIndex=-1;
			foreach my $codeID (@{$analysisGroups[$g]}) {
				my ($type,$typeID)=split(':',$codeID); # A or C
				$codeIndex++;
				next if ($type eq 'C' && !$catProteins{$typeID}{$protID});

				my @localAnaList=($type eq 'A')? ($typeID) : @parentItemAna;
				my $anaIndex=-1; # index of Analyses in group (not A/C items)
				foreach my $anaID (@localAnaList) {
					$anaIndex++;

					next unless $protAnalyses{$protID}{$anaID}; # protein must be in analysis
					next if $processedProtAna{$protID}{$anaID}; # protein must be not processed
					if ($type eq 'C' && !$groupListCompProt{$g}{$protID}) { # List comparison => skip if prot is not in (group of) list(s)
						#$peptideAnaDistrib{$protID}[$g][$codeIndex][$anaIndex]=0;
						next;
					}
					my (%distinctPepAna,%allPepAna);

					if ($protAnalyses{$protID}{$anaID}[0]<0) { # not yet updated
						################ v1.3.2
						#if ($peptideRule=~/nr/) {
						#	$sthPep->execute($protID,$anaID);
						#	while (my ($pepSeq,$varMod)=$sthPep->fetchrow_array) {
						#		my $varModStrg=($varMod)? " $varMod" : '';
						#		$distinctPepAna{"$pepSeq$varModStrg"}=1;
						#	}
						#	$protAnalyses{$protID}{$anaID}[0]=scalar keys %distinctPepAna;
						#}
						#else {$protAnalyses{$protID}{$anaID}[0]=abs($protAnalyses{$protID}{$anaID}[0]);}

						$sthPep->execute($protID,$anaID);
						#while (my ($pepSeq,$varMod)=$sthPep->fetchrow_array) {
						#	my $varModStrg='';
						#	if ($varMod) {
						#		$varModStrg=" $varMod";
						#		&getVarModCodes(\%{$ptmAnalyses{$protID}{$g}{$anaID}},$varModStrg,"$pepSeq$varModStrg") if scalar keys %projectVarMods>=1;
						#	}
						#	$distinctPepAna{"$pepSeq$varModStrg"}=1 if $peptideRule=~/nr/;
						#}
						while (my ($pepID,$pepSeq)=$sthPep->fetchrow_array) {
							$allPepAna{$pepID}=1;
							my ($varModStrg,@modifications)=('',());
							if ($anaPepVmods{$anaID}{$pepID}) {
								foreach my $modID (sort{$a<=>$b} keys %{$anaPepVmods{$anaID}{$pepID}}) {
									$varModStrg.="+$modID:$anaPepVmods{$anaID}{$pepID}{$modID}";
									if ($projectVarMods{$modID}) {
										$ptmAnalyses{$protID}{$g}{$anaID}{$modID}{'ALL'}++;
										push @modifications,$modID;
									}
								}
							}
							if ($peptideRule=~/nr/) {# prevents INTRA-analysis redundancy
								foreach my $modID (@modifications) {
									$ptmAnalyses{$protID}{$g}{$anaID}{$modID}{'NR'}{"$pepSeq$varModStrg"}=1;
								}
							}
							$distinctPepAna{"$pepSeq$varModStrg"}=1 if $peptideRule=~/nr/;
						}
						if ($peptideRule=~/nr/) {
							$protAnalyses{$protID}{$anaID}[0]=scalar keys %distinctPepAna;
						}
						else {
							#$protAnalyses{$protID}{$anaID}[0]=abs($protAnalyses{$protID}{$anaID}[0]);
							$protAnalyses{$protID}{$anaID}[0]=scalar keys %allPepAna;
						}
						###############################
					}

					if ($type eq 'A') {
						$peptideAnaDistrib{$protID}[$g][$codeIndex]=$protAnalyses{$protID}{$anaID}[0] if ($numGroups || $listInGroup{$g});
#print "ANA2 $protID:$g:$codeIndex => $peptideAnaDistrib{$protID}[$g][$codeIndex]<BR>\n";
					}
					else { # C
						$peptideAnaDistrib{$protID}[$g][$codeIndex][$anaIndex]=$protAnalyses{$protID}{$anaID}[0]; # 1 more dimension!!!
#print "CAT2 $protID:$g:$codeIndex:$anaIndex => $peptideAnaDistrib{$protID}[$g][$codeIndex][$anaIndex]<BR>\n";
					}
#print "*$protID:$g<BR>\n";
					next if $predefinedProtGr;
#print "**Pass:&nbsp;";
					#next if ($proteinList{$protID}[$g] && $proteinList{$protID}[$g][1]); # prot is visible in group

					if (!defined($proteinList{$protID}[$g])) { # 1st ana in group
						@{$proteinList{$protID}[$g]}=($protAnalyses{$protID}{$anaID}[0],0,$protAnalyses{$protID}{$anaID}[2]);
						$bestScore{$protID}=$protAnalyses{$protID}{$anaID}[3];
						$bestAnalysis{$protID}{$g}=$anaID if $peptideRule=~/ba/;
					}
					elsif ($peptideRule=~/ba/) { # prot is hidden in >1 ana
#print "*ba: $protID, $g, $anaID<BR>\n";
						if ($protAnalyses{$protID}{$anaID}[0]>$proteinList{$protID}[$g][0] || ($protAnalyses{$protID}{$anaID}[0]==$proteinList{$protID}[$g][0] && $protAnalyses{$protID}{$anaID}[3]>$bestScore{$protID})) {
							$proteinList{$protID}[$g][0]=$protAnalyses{$protID}{$anaID}[0]; # numpep
							$proteinList{$protID}[$g][2]=$protAnalyses{$protID}{$anaID}[2]; # confLevel
							$bestScore{$protID}=$protAnalyses{$protID}{$anaID}[3];
							$bestAnalysis{$protID}{$g}=$anaID;
						}
					}
					else { # prot is hidden in >1 ana
#print "*else: $protID, $g, $anaID<BR>\n";
						$proteinList{$protID}[$g][0]+=$protAnalyses{$protID}{$anaID}[0]; # numpep
						$proteinList{$protID}[$g][2]=$protAnalyses{$protID}{$anaID}[2] if $proteinList{$protID}[$g][2]<$protAnalyses{$protID}{$anaID}[2]; # confLevel
					}
#print "-$protID:$g => $proteinList{$protID}[$g][0]<BR>\n";
					if ($peptideRule=~/nr/) {
						if (scalar keys %distinctPepAna) { # peptides already scanned for $protAnalyses{$protID}{$anaID}[0]
							foreach my $pepMod (keys %distinctPepAna) {
								$protPeptides{$protID}{$pepMod}=1;
							}
						}
						#else { #Commented in v1.3.2
						#	$sthPep->execute($protID,$anaID);
						#	while (my ($pepSeq,$varMod)=$sthPep->fetchrow_array) {
						#		my $varModStrg=($varMod)? " $varMod" : '';
						#		$protPeptides{$protID}{"$pepSeq$varModStrg"}=1;
						#	}
						#}
					}
				}
			} # end of anaID loop
			if (!$predefinedProtGr && $proteinList{$protID}[$g]) { # prot is hidden in at least 1 analysis
				$proteinList{$protID}[$g][0]=scalar keys %{$protPeptides{$protID}} if $peptideRule=~/nr/; # Correcting numPep in %proteinList
#print "+$protID:$g => $proteinList{$protID}[$g][0]<BR>\n";
				$groupPeptides[$g]{$proteinList{$protID}[$g][0]}=1;# if $proteinList{$protID}[$g][0]>0; # to select num peptides for auto-check
			}
			#else {@{$proteinList{$protID}[$g]}=(0,0,0);}
			elsif (!$proteinList{$protID}[$g]) {@{$proteinList{$protID}[$g]}=(0,0,0);}
		} # end of $hiddenRule
		elsif (!$proteinList{$protID}[$g]) {@{$proteinList{$protID}[$g]}=(0,0,0);}
	}
}
$sthPep->finish;
$sthPepMod->finish;


####>Peptide threshold<####
#my $pepThreshold=5;
if ($pepThreshold > 1) {
	&updateProgress("Filtering data (nb. peptides &ge; $pepThreshold)...") unless $exportType;
	foreach my $protID (keys %proteinList) {
		my $okProt=0;
		foreach my $g (0..$#{$proteinList{$protID}}) {
			if ($proteinList{$protID}[$g][0] >= $pepThreshold) {
				$okProt++;
			}
			else {@{$proteinList{$protID}[$g]}=(0,0,0);}
		}
		unless ($okProt) { ## Delete prot everywhere
			delete $proteinList{$protID};
			delete $peptideAnaDistrib{$protID};
			delete $selectedProteins{$protID};

			#delete $ptmAnalyses{$protID};
			#delete $processedProtAna{$protID};
			#delete $protAnalyses{$protID};
			#delete $groupProteins{$protID};
			#delete $bestScore{$protID};
			#delete $bestAnalysis{$protID};
			#delete $groupBestVis{$protID};
			#delete $protPeptides{$protID};
			#delete $itemAnaScanIdx{$protID};
		}
	}
}

###################################
####>1 vs 1 Protein Comparison<####
###################################
if ($compType eq '1v1_prot') {

	$dbh->disconnect;
	#&updateProgress('Post processing data...') unless $exportType;

	my (%proteinDistrib,%groupDistrib);
	foreach my $protID (keys %proteinList) {
		foreach my $g (0..$#analysisGroups) {
			if ($proteinList{$protID}[$g][0]) {
				$proteinDistrib{$g}{$protID}=1;
				$groupDistrib{$protID}{$g}=1;
			}
		}
	}
	&display1v1ComparisonBody('protein',\@groupLabels,\%proteinDistrib,\%groupDistrib);
	exit;
}


&updateProgress('Computing peptide distribution...') unless $exportType;
####>Peptide distrib in grouped Ana<####
foreach my $protID (keys %peptideAnaDistrib) {
	foreach my $g (0..$#{$peptideAnaDistrib{$protID}}) {
		foreach my $codeIndex (0..$#{$analysisGroups[$g]}) {
			#$peptideAnaDistrib{$protID}[$g][$anaIndex]=0 unless $peptideAnaDistrib{$protID}[$g][$anaIndex];
			if ($peptideAnaDistrib{$protID}[$g][$codeIndex]) {
				if (ref($peptideAnaDistrib{$protID}[$g][$codeIndex])) { # List comparison => 1 more dimension
					foreach my $anaIndex (0..$#parentItemAna) {
						$peptideAnaDistrib{$protID}[$g][$codeIndex][$anaIndex]=0 unless $peptideAnaDistrib{$protID}[$g][$codeIndex][$anaIndex];
					}
				}
			}
			else {$peptideAnaDistrib{$protID}[$g][$codeIndex]=0;}
		}
	}
}

####>PTMs distrib in grouped Ana<####
my (%protPtmGroups,%ptmOccurences,%ptmProteins,%proteinPtms);
foreach my $protID (keys %ptmAnalyses) { # $projectPtmString must be not null
	foreach my $g (keys %{$ptmAnalyses{$protID}}) {
		my $ptmInGr=0;
		my %nrPepPtms;
		foreach my $anaID (keys %{$ptmAnalyses{$protID}{$g}}) {
			next if ($peptideRule=~/^ba/ && $anaID != $bestAnalysis{$protID}{$g});
			foreach my $varMod (keys %{$ptmAnalyses{$protID}{$g}{$anaID}}) {
				if ($peptideRule=~/nr/) {
					foreach my $pep (keys %{$ptmAnalyses{$protID}{$g}{$anaID}{$varMod}{'NR'}}) {
						$nrPepPtms{$varMod}{$pep}=1; # prevents INTER-analysis redundancy
					}
				}
				else {
					$protPtmGroups{$protID}{$g}{$varMod}+=$ptmAnalyses{$protID}{$g}{$anaID}{$varMod}{'ALL'};
				}
			}
		}
		if ($peptideRule=~/nr/) {
			foreach my $varMod (keys %nrPepPtms) {
				$protPtmGroups{$protID}{$g}{$varMod}=scalar keys %{$nrPepPtms{$varMod}};
			}
		}
		foreach my $varMod (keys %{$protPtmGroups{$protID}{$g}}) {
			$ptmOccurences{$g}{$protPtmGroups{$protID}{$g}{$varMod}}=1;
			$ptmProteins{$varMod}{$protID}=1;
			$proteinPtms{$protID}{$varMod}=1;
		}
	}
}

#print "-->(1:6): $peptideAnaDistrib{202}[1][0]<BR>\n";
####>Computing Relative Delta % (only if 2 ana/groups)<####
my (%absoluteDelta,%relativeDelta);
if ($numAnaGroups==2) {
	foreach my $protID (sort keys %proteinList) {
		$absoluteDelta{$protID}=$proteinList{$protID}[1][0]-$proteinList{$protID}[0][0];
		$relativeDelta{$protID}=$absoluteDelta{$protID}*100/($proteinList{$protID}[0][0]+$proteinList{$protID}[1][0]);
		#$relativeDelta{$protID}=($proteinList{$protID}[1][0]-$proteinList{$protID}[0][0])*100/($proteinList{$protID}[0][0]+$proteinList{$protID}[1][0]);
	}
}

####>Fetching protein annotation & distribution<####
&updateProgress('Fetching protein annotation...') unless $exportType;
my (%protAnnotation,%protDistribution,%masterProteins);
my $sthProt=$dbh->prepare("SELECT ALIAS,PROT_DES,MW,ORGANISM,ID_MASTER_PROTEIN FROM PROTEIN WHERE ID_PROTEIN=?");
my ($geneNameID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
my $sthMPG=$dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN=? AND ID_IDENTIFIER=$geneNameID ORDER BY RANK");
foreach my $protID (keys %proteinList) {
	$sthProt->execute($protID);
	@{$protAnnotation{$protID}}=$sthProt->fetchrow_array;
	$protAnnotation{$protID}[2]=sprintf "%.1f",$protAnnotation{$protID}[2]/1000; # MW
	if ($sortOrder eq 'peptide') {
		@{$protDistribution{$protID}}=(0,0,0);
		foreach my $g (0..$#analysisGroups) {
			$protDistribution{$protID}[0]++ if $proteinList{$protID}[$g][0]; # numb ana where protein is found
			$protDistribution{$protID}[2]+=$proteinList{$protID}[$g][0]; # total num peptides
		}
		my $avePep=$protDistribution{$protID}[2]/$numAnaGroups;
		foreach my $g (0..$#analysisGroups) {
			$protDistribution{$protID}[1]+=abs($proteinList{$protID}[$g][0]-$avePep); # sum of deltas to average
		}
		$protDistribution{$protID}[1]=sprintf "%.10f",$protDistribution{$protID}[1]/$protDistribution{$protID}[2]; # normalized sum of deltas to average
	}
	if ($protAnnotation{$protID}[4]) {
		my $masterProtID=$protAnnotation{$protID}[4];
		unless ($masterProteins{$masterProtID}) {
			@{$masterProteins{$masterProtID}}=();
			$sthMPG->execute($masterProtID);
			while (my ($gene)=$sthMPG->fetchrow_array) {
				push @{$masterProteins{$masterProtID}},$gene;
			}
		}
	}
}
$sthProt->finish;
$sthMPG->finish;


#####################
####>Export mode<####
#####################
if ($exportType) {

	$dbh->disconnect;
	#my $timeStamp1=strftime("%Y%m%d %H-%M",localtime);
	my $timeStamp2=strftime("%Y-%m-%d %H:%M",localtime);

	my $workbook=Spreadsheet::WriteExcel->new("-");
	eval { # Perl 5.8 compatibility
		$workbook->set_properties(title=>"Compare Project Items",
								  author=>'myProMS server',
								  comments=>'Automatically generated with Perl and Spreadsheet::WriteExcel'
								  );
	};
	$workbook->set_custom_color(40,189,215,255); # light light color  #BDD7FF (new myProMS colors V3.5+)
	#$workbook->set_custom_color(41,128,179,255); # dark blue color  #80B3FF
	my %itemFormat=(
			title =>			$workbook->add_format(align=>'center',size=>18,bold=>1,border=>1),
			header =>			$workbook->add_format(align=>'center',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			mergeRowHeader =>	$workbook->add_format(align=>'center',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			mergeColHeader =>	$workbook->add_format(align=>'center',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			identVis =>			$workbook->add_format(align=>'left',valign=>'top',size=>10,border=>1,bold=>1),
			identHid =>			$workbook->add_format(align=>'left',valign=>'top',size=>10,border=>1),
			identVisBC =>			$workbook->add_format(align=>'left',valign=>'top',size=>10,border=>1,bold=>1,color=>'grey'),
			identHidBC =>			$workbook->add_format(align=>'left',valign=>'top',size=>10,border=>1,color=>'grey'),
			identVisVP =>			$workbook->add_format(align=>'left',valign=>'top',size=>10,border=>1,bold=>1,italic=>1),
			identHidVP =>			$workbook->add_format(align=>'left',valign=>'top',size=>10,border=>1,italic=>1),
#mergeIdent =>		$workbook->add_format(align=>'left',valign=>'top',size=>10,border=>1),
			identBC =>			$workbook->add_format(align=>'left',size=>10,color=>'grey',border=>1),
#mergeIdentBC =>		$workbook->add_format(align=>'left',valign=>'top',size=>10,color=>'grey',border=>1),
			identVP =>			$workbook->add_format(align=>'left',size=>10,italic=>1,border=>1),
#mergeIdentVP =>		$workbook->add_format(align=>'left',valign=>'top',size=>10,italic=>1,border=>1),
			text =>				$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>0,border=>1),
			textWrap =>				$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>1,border=>1),
			mergeText =>		$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>0,border=>1),
			mergeColText =>		$workbook->add_format(align=>'left',size=>10,text_wrap=>0,border=>1),
			numberVis =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,bold=>1),
			numberHid =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1),
			numberVisBC =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,bold=>1,color=>'grey'),
			numberHidBC =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,color=>'grey'),
			numberVisVP =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,bold=>1,italic=>1),
			numberHidVP =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,italic=>1),
			number =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1),
			mergeNumber =>		$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1),
			number1d =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',num_format=>'0.0',border=>1),
#mergeNumber1d =>	$workbook->add_format(size=>10,align=>'center',valign=>'top',num_format=>'0.0',border=>1)
			);

	my $worksheet=$workbook->add_worksheet('Results');

	#print header(-type=>"application/vnd.ms-excel",-attachment=>"Compare items_$timeStamp1.xls");

	my ($ptmFactor,$mainRowSpan,$itemString)=(scalar keys %projectVarMods>=1)? (1,3,'Peptides/PTMs') : (0,2,'Peptides');
	my $deltaRowSpan=$mainRowSpan-1;

	my (@pepColSpans,@ptmColSpans,@grColSpans,$peptideString,$itemColSpan);
	if ($numGroups || !$listComparison) { # gr comparison OR simple Ana comparison
		(my $pepColSpan,my $ptmColSpan,$peptideString)=($numGroups)? (3,4,"$itemString, Analyses, Distribution") : (1,2,$itemString);
		my $grColSpan=$pepColSpan+($ptmColSpan*$ptmFactor);
		foreach my $g (0..$#analysisGroups) {
			push @pepColSpans,$pepColSpan;
			push @ptmColSpans,$ptmColSpan;
			push @grColSpans,$grColSpan;
			$itemColSpan+=$grColSpan;
		}
		$itemColSpan+=2 if ($numAnaGroups==2 && !$listComparison); # for (abs/rel)Delta
	}
	elsif ($listComparison) { # list +/- Ana comparison
		$peptideString="$itemString, Analyses, Distribution";
		$itemColSpan=0;
		foreach my $g (0..$#analysisGroups) {
			my ($pepColSpan,$ptmColSpan)=($listInGroup{$g})? (3,4) : (1,2);
			my $grColSpan=$pepColSpan+($ptmColSpan*$ptmFactor);
			push @pepColSpans,$pepColSpan;
			push @ptmColSpans,$ptmColSpan;
			push @grColSpans,$grColSpan;
			$itemColSpan+=$grColSpan;
		}
		#$itemColSpan+=2 if $numAnaGroups==2;
	}
	my $lastColIndex=$itemColSpan+4;

	####>Start printing<####
	my $xlsRow=0;
	my $xlsCol=0;
	$worksheet->set_column(0,0,30); # identifier col length

	##>Title & list context
	$worksheet->merge_range($xlsRow,0,$xlsRow,$lastColIndex,"Compare Analyses and/or Lists ($timeStamp2)",$itemFormat{'title'});
	$worksheet->merge_range(++$xlsRow,0,$xlsRow,$lastColIndex,$listContextStrg,$itemFormat{'mergeColText'}) if $listComparison;
	##>List filter
	my $catFilterStrg;
	if ($catFilter) {
		if ($catFilterRule eq 'exclude') {$catFilterStrg="Proteins not in Theme '$filterClassName' > List '$filterCatName'";}
		else {$catFilterStrg="Proteins restricted to Theme '$filterClassName' > List '$filterCatName'";}
	}
	else {$catFilterStrg='No List restriction';}
	$worksheet->merge_range(++$xlsRow,0,$xlsRow,$lastColIndex,$catFilterStrg,$itemFormat{'mergeColText'});
	##>Peptide rules
	if ($peptideRule) {
		my $peptideRuleStrg='';
		$peptideRuleStrg="Peptide rule: At least $pepThreshold ";
		my $pepStrg=($pepThreshold==1)? 'peptide' : 'peptides';
		if ($numGroups || $listComparison) {
			$peptideRuleStrg.=($peptideRule eq 'all')? "$pepStrg in Group" : ($peptideRule eq 'anr')? "distinct $pepStrg in Group" : ($peptideRule eq 'ba')? "$pepStrg in best Analysis" : "distinct $pepStrg in best Analysis";
		}
		else {
			$peptideRuleStrg.=($peptideRule=~/nr/)? "$pepStrg distinct in Analysis" : "$pepStrg in Analysis";
		}
		$worksheet->merge_range(++$xlsRow,0,$xlsRow,$lastColIndex,$peptideRuleStrg,$itemFormat{'mergeColText'});
	}
	##>Table header
	my $proteinStrg=($exportType eq 'all')? 'All proteins' : 'Selected proteins';
	$worksheet->merge_range(++$xlsRow,0,$xlsRow+$mainRowSpan-1,0,$proteinStrg,$itemFormat{'mergeRowHeader'});
	$worksheet->merge_range($xlsRow,1,$xlsRow+$mainRowSpan-1,1,'Gene',$itemFormat{'mergeRowHeader'});
	$xlsCol=$itemColSpan+1;
	$worksheet->merge_range($xlsRow,2,$xlsRow,$xlsCol,$peptideString,$itemFormat{'mergeColHeader'});
	$worksheet->merge_range($xlsRow,++$xlsCol,$xlsRow+$mainRowSpan-1,$xlsCol,'MW (kDa)',$itemFormat{'mergeRowHeader'});
	$worksheet->set_column(++$xlsCol,$xlsCol,80); # col length
	$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow+$mainRowSpan-1,$xlsCol,'Description',$itemFormat{'mergeRowHeader'});
	$worksheet->set_column(++$xlsCol,$xlsCol,25); # col length
	$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow+$mainRowSpan-1,$xlsCol,'Species',$itemFormat{'mergeRowHeader'});
	$xlsRow++;
	$xlsCol=2;
	$worksheet->set_row($xlsRow,15*(1+$maxItemInGr)) if $maxItemInGr > 1; # row height
#my $smiley = "\x{263A}";
	foreach my $g (0..$#analysisGroups) {
		my $labelStrg=($numGroups)? "$groupLabels[$g][0]:$newLine$groupLabels[$g][1]" : $groupLabels[$g][1];
		my $endCol=$xlsCol+$grColSpans[$g]-1;
		if ($grColSpans[$g] > 1) {$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow,$endCol,decode_utf8($labelStrg),$itemFormat{'mergeColHeader'});}
		else {$worksheet->write_string($xlsRow,$xlsCol,decode_utf8($labelStrg),$itemFormat{'header'});}
		$xlsCol=$endCol+1;
	}
	if ($numAnaGroups==2 && !$listComparison) {
		if ($deltaRowSpan==1) { # in case no relevant PTMs in porject (scalar keys %projectVarMods=0)
			$worksheet->write_string($xlsRow,$xlsCol,'Abs. Delta (%)',$itemFormat{'header'});
			$worksheet->write_string($xlsRow,++$xlsCol,'Rel. Delta (%)',$itemFormat{'header'});
		}
		else {
			$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow+$deltaRowSpan-1,$xlsCol,'Abs. Delta (%)',$itemFormat{'mergeRowHeader'});
			$worksheet->merge_range($xlsRow,++$xlsCol,$xlsRow+$deltaRowSpan-1,$xlsCol,'Rel. Delta (%)',$itemFormat{'mergeRowHeader'});
		}
	}
	if (scalar keys %projectVarMods>=1) {
		$xlsRow++;
		$xlsCol=2;
		foreach my $g (0..$#analysisGroups) {
			my $endCol=$xlsCol+$pepColSpans[$g]-1;
			if ($pepColSpans[$g] > 1) {$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow,$endCol,'Peptides',$itemFormat{'mergeColHeader'});}
			else {$worksheet->write_string($xlsRow,$xlsCol,'Peptides',$itemFormat{'header'});}
			$xlsCol=$endCol+1;
			$endCol=$xlsCol+$ptmColSpans[$g]-1;
			if ($ptmColSpans[$g] > 1) {$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow,$endCol,'PTMs',$itemFormat{'mergeColHeader'});}
			else {$worksheet->write_string($xlsRow,$xlsCol,'Peptides',$itemFormat{'header'});}
			$xlsCol=$endCol+1;
		}
	}

	###>Looping through proteins<###
	foreach my $protID (sort{&sortProt($sortOrder)} keys %protAnnotation) {
		$xlsCol=0;
		##>identifier
		my $visTag=($bestVisibility{$protID}[1])? 'Vis' : 'Hid';
		#my $confTag=($confLevel{$protID}==0)? 'VP' : ($confLevel{$protID}==1)? 'BC' : '';
		$worksheet->write_string(++$xlsRow,$xlsCol,$protAnnotation{$protID}[0],$itemFormat{'ident'.$visTag}); #.$confTag
		##<Gene
		my $geneName=($protAnnotation{$protID}[4] && $masterProteins{$protAnnotation{$protID}[4]}[0])? $masterProteins{$protAnnotation{$protID}[4]}[0] : '-';
		$worksheet->write_string($xlsRow,++$xlsCol,$geneName,$itemFormat{'ident'.$visTag});
		##>Analyses/groups
		foreach my $g (0..$#analysisGroups) {
			if ($proteinList{$protID}[$g][0]) { # prot exist in gr
				#>peptide
				my $gVis=($proteinList{$protID}[$g][1])? 'Vis' : 'Hid';
				my $gConf=($proteinList{$protID}[$g][2]==0)? 'VP' : ($proteinList{$protID}[$g][2]==1)? 'BC' : '';
				$worksheet->write_number($xlsRow,++$xlsCol,$proteinList{$protID}[$g][0],$itemFormat{'number'.$gVis.$gConf});
				#>Distribution
				my $numAnaExistPep=0;
				my $pepAnaDistribStrg='';
				my $itemAnaDistribStrg='';
				my $hiddenTag='';
				if (($numGroups || $listInGroup{$g}) && $peptideAnaDistrib{$protID}[$g]) { # Groups and/or Lists
					my $itemAnaScan=0; # if C only scan parent item analyses once
					my $codeIndex=0; #my $anaOnLine=0;
					foreach my $codeID (@{$analysisGroups[$g]}) {
						$pepAnaDistribStrg.='-' if $codeIndex;
						my ($type,$typeID)=split(':',$codeID); # A or C
						if ($type eq 'A') {
							$numAnaExistPep++ if $peptideAnaDistrib{$protID}[$g][$codeIndex];
							$pepAnaDistribStrg.=$peptideAnaDistrib{$protID}[$g][$codeIndex];
							$pepAnaDistribStrg.='*' if ($peptideAnaDistrib{$protID}[$g][$codeIndex] && !$protAnalyses{$protID}{$typeID}[1]); # * = tag for hidden
							$hiddenTag='*' if ($hiddenRule==1 && $proteinList{$protID}[$g][1] && $peptideAnaDistrib{$protID}[$g][$codeIndex] && !$protAnalyses{$protID}{$typeID}[1]); # prot is vis in group but hidden in some analyses
						}
						else { # C
							if ($peptideAnaDistrib{$protID}[$g][$codeIndex]) { # protein exists in analysis/list   $peptideAnaDistrib{$protID}[$g][$codeIndex]
								$numAnaExistPep++;
								if ($itemAnaScan) { # ana list already scanned
									$pepAnaDistribStrg.='L'; # at least 2 Lists in same group
									$codeIndex++;
									next;
								}
								my $anaIndex=0;
								foreach my $anaID (@parentItemAna) {
									$itemAnaDistribStrg.='-' if $anaIndex;
									$itemAnaDistribStrg.=$peptideAnaDistrib{$protID}[$g][$codeIndex][$anaIndex];
									$itemAnaDistribStrg.='*' if ($peptideAnaDistrib{$protID}[$g][$codeIndex][$anaIndex] && !$protAnalyses{$protID}{$anaID}[1]);
									$hiddenTag='*' if ($hiddenRule==1 && $proteinList{$protID}[$g][1] && $peptideAnaDistrib{$protID}[$g][$codeIndex][$anaIndex] && !$protAnalyses{$protID}{$anaID}[1]); # prot is vis in group but hidden in some analyses
									$anaIndex++;
								}
								if ($listInGroup{$g} > 1) {$pepAnaDistribStrg.='L';}
								else {$pepAnaDistribStrg.="[$itemAnaDistribStrg]";}
								$itemAnaScan=1;
							}
							else {$pepAnaDistribStrg.='0';}
						}
						$codeIndex++;
					}
					$worksheet->write_number($xlsRow,++$xlsCol,$numAnaExistPep,$itemFormat{'number'.$gVis.$gConf});
					$worksheet->write_string($xlsRow,++$xlsCol,$pepAnaDistribStrg,$itemFormat{'text'});
					my $commentStrg='';
					$commentStrg.="L=$itemAnaDistribStrg." if ($listInGroup{$g} && $listInGroup{$g} > 1 && $itemAnaDistribStrg);
					if ($hiddenTag) {
						$commentStrg.=$newLine if $commentStrg;
						$commentStrg.='*: Protein is hidden in some Analyses.';
					}
					$worksheet->write_comment($xlsRow,$xlsCol,$commentStrg) if $commentStrg;
				}
				#>PTMs
				next unless scalar keys %projectVarMods>=1;
				if ($protPtmGroups{$protID} && $protPtmGroups{$protID}{$g}) {
					my $textFormat=(scalar keys %{$protPtmGroups{$protID}{$g}} > 1)? 'textWrap' : 'text';
					my $ptmNameStrg='';
					foreach my $varMod (sort{$a cmp $b} keys %{$protPtmGroups{$protID}{$g}}) {
						$ptmNameStrg.=$newLine if $ptmNameStrg;
						$ptmNameStrg.="$projectVarMods{$varMod}[0]:";
					}
					$worksheet->write_string($xlsRow,++$xlsCol,$ptmNameStrg,$itemFormat{$textFormat});
					my $ptmNumStrg='';
					foreach my $varMod (sort{$a cmp $b} keys %{$protPtmGroups{$protID}{$g}}) {
						$ptmNumStrg.=$newLine if $ptmNumStrg;
						$ptmNumStrg.=$protPtmGroups{$protID}{$g}{$varMod};
					}
					$worksheet->write_string($xlsRow,++$xlsCol,$ptmNumStrg,$itemFormat{$textFormat});

					if (($numGroups || $listInGroup{$g}) && $peptideAnaDistrib{$protID}[$g]) { # PTM distrib data
						my $allPtmOccString='';
						my $allPtmDistribString='';
						foreach my $varMod (sort{$a cmp $b} keys %{$protPtmGroups{$protID}{$g}}) {
							my %numAnaExistPtm;
							my $ptmAnaDistribStrg='';
							my $firstCode=1;
							foreach my $codeID (@{$analysisGroups[$g]}) {
								$ptmAnaDistribStrg.='-' unless $firstCode;
								my ($type,$typeID)=split(':',$codeID); # A or C
								my @localAnaList;
								if ($type eq 'A') {
									@localAnaList=($typeID);
								}
								else { # C
									@localAnaList=@parentItemAna;
									$ptmAnaDistribStrg.='[';
								}
								my $firstAna=1;
								foreach my $anaID (@localAnaList) {
									my $numProtPtmGr=0;
									if ($ptmAnalyses{$protID}{$g}{$anaID} && $ptmAnalyses{$protID}{$g}{$anaID}{$varMod}) {
										$numAnaExistPtm{$anaID}=1;
										$numProtPtmGr=($peptideRule=~/nr/)? scalar keys %{$ptmAnalyses{$protID}{$g}{$anaID}{$varMod}{'NR'}} : $ptmAnalyses{$protID}{$g}{$anaID}{$varMod}{'ALL'};
									}
									$ptmAnaDistribStrg.='-' unless $firstAna;
									$ptmAnaDistribStrg.=$numProtPtmGr;
									$firstAna=0;
								}
								$ptmAnaDistribStrg.=']' if $type eq 'C';
								$firstCode=0;
							}
							$allPtmOccString.=$newLine if $allPtmOccString;
							$allPtmOccString.=scalar keys %numAnaExistPtm;
							$allPtmDistribString.=$newLine if $allPtmDistribString;
							$allPtmDistribString.=$ptmAnaDistribStrg;
						}
						$worksheet->write_string($xlsRow,++$xlsCol,$allPtmOccString,$itemFormat{$textFormat});
						$worksheet->write_string($xlsRow,++$xlsCol,$allPtmDistribString,$itemFormat{$textFormat});
					}
				}
				else {
					my $endCol=++$xlsCol+$ptmColSpans[$g]-1;
					$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow,$endCol,'-',$itemFormat{'mergeNumber'});
					$xlsCol=$endCol;
				}
			}
			else {
				my $endCol=++$xlsCol+$grColSpans[$g]-1;
				if ($endCol-$xlsCol > 1) {$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow,$endCol,'-',$itemFormat{'mergeNumber'});}
				else {$worksheet->write_string($xlsRow,$xlsCol,'-',$itemFormat{'number'});}
				$xlsCol=$endCol;
			}
		}
		##>Absolute & Relative Delta
		if ($numAnaGroups==2 && !$listComparison) {
			$worksheet->write_number($xlsRow,++$xlsCol,$absoluteDelta{$protID},$itemFormat{'number'});
			$worksheet->write_number($xlsRow,++$xlsCol,$relativeDelta{$protID},$itemFormat{'number1d'});
		}
		##>Mass
		if ($protAnnotation{$protID}[2]) {
			$worksheet->write_number($xlsRow,++$xlsCol,$protAnnotation{$protID}[2],$itemFormat{'number1d'});
		}
		else {
			$worksheet->write_string($xlsRow,++$xlsCol,'-',$itemFormat{'text'});
		}
		##>Description
		$worksheet->write_string($xlsRow,++$xlsCol,$protAnnotation{$protID}[1],$itemFormat{'text'});
		##>Species
		$worksheet->write_string($xlsRow,++$xlsCol,$protAnnotation{$protID}[3],$itemFormat{'text'});
	}

	$workbook->close();
	exit;
}

####>Fetching list of available classification and categories<####
my (%classifications,%categories);
my $sthClass=$dbh->prepare("SELECT ID_CLASSIFICATION,NAME FROM CLASSIFICATION WHERE ID_PROJECT=$projectID");
my $sthCat=$dbh->prepare("SELECT ID_CATEGORY,NAME,DISPLAY_POS FROM CATEGORY WHERE ID_CLASSIFICATION=?");
$sthClass->execute;
while (my ($classID,$className)=$sthClass->fetchrow_array) {
	$classifications{$classID}=$className;
	$sthCat->execute($classID);
	while (my ($catID,$catName,$displayPos)=$sthCat->fetchrow_array) {
		@{$categories{$classID}{$catID}}=($catName,$displayPos);
	}
}
$sthClass->finish;
$sthCat->finish;

$dbh->disconnect;

####>Fetching Comparison autocheck parameters (if any)<####
my ($filterLogic,%autocheckParams,$selectedRefPTM,$okAutocheck);
if ($autochkStrg) {
	#my $autochkStrg='1:>=:3;2:>=:4';
	($filterLogic,my $pepFilterStrg,$selectedRefPTM,my $ptmFilterStrg)=split('#',$autochkStrg);
	$selectedRefPTM=undef if ($selectedRefPTM && !$projectVarMods{$selectedRefPTM}); # in case project PTMs list has changed after comparison was saved
	if ($filterLogic eq 'frequency') {
		($autocheckParams{'FREQ_ANAGR'},$autocheckParams{'FREQ_PEP'})=split(';',$pepFilterStrg);
		$autocheckParams{'FREQ_PTM'}=$ptmFilterStrg if $selectedRefPTM;
		$okAutocheck=1;
	}
	else { # and,or
		foreach my $grData (split(';',$pepFilterStrg)) {
			my ($aGr,$rel,$numPep)=split(':',$grData);
			@{$autocheckParams{'PEP'}{$aGr}}=($rel,$numPep);
			$okAutocheck=1 if ($rel ne '>=' || $numPep>0); # do not check all proteins if ">=0" for all anaGroups
		}
		if ($selectedRefPTM) {
			foreach my $grData (split(';',$ptmFilterStrg)) {
				my ($aGr,$rel,$numPtm)=split(':',$grData);
				@{$autocheckParams{'PTM'}{$aGr}}=($rel,$numPtm);
				$okAutocheck=1 if ($rel ne '>=' || $numPtm>0); # do not check all proteins if ">=0" for all anaGroups
			}
		}
	}
}
else {$filterLogic='and';}
my $numProtChecked=scalar keys %selectedProteins;
$okAutocheck=0 if $numProtChecked; # check by Perl

####>Computing Venn diagram<####
my %vennDiagram;
if ($numAnaGroups<=5) {
	foreach my $protID (keys %proteinList) {
		my  @protSet;
		foreach my $g (0..$#analysisGroups) {
			$protSet[$g]=($proteinList{$protID}[$g][0])? chr(65+$g) : '';
		}
		$vennDiagram{join('',@protSet)}++;
	}
}

#######################
####>Starting BODY<####
#######################
my $numAllProteins=scalar keys %proteinList;
if ($numAllProteins==0) {
	print qq
|<SCRIPT type="text/javascript">document.getElementById('waitDIV').style.display='none';</SCRIPT>
<BR><BR><FONT class="title2">No proteins matched!</FONT>
</CENTER>
</BODY>
</HTML>
|;
	exit;
}

#print "SETS:<BR>"; foreach my $set (keys %vennDiagram) {print ">$set => $vennDiagram{$set}<BR>\n";}; exit;
my $onloadStrg=''; #onload="refreshChkProteins()"';
#$onloadStrg.="autoCheck(0);" if ($okAutocheck && scalar keys %proteinList);
$onloadStrg.=($okAutocheck && $numAllProteins)? "updateComparisonStrategy('$filterLogic'); autoCheck(0); hideShowUnchecked(0);" : ($numProtChecked)? 'update2ndVennDiagram();' : ''; # $numProtChecked after Export?, save list or change in list display
#$onloadStrg.="displayPtmImages(document.getElementById('refPTM').value);" if $projectPtmString;

###>PTMs
print qq
|<SCRIPT type="text/javascript">
//PTMs data
var usedPTMs=new Object();
|;
if (scalar keys %ptmProteins) {
	foreach my $varMod (sort{lc($a) cmp lc($b)} keys %ptmProteins) {
		print "usedPTMs['$varMod']=[",join(',',keys %{$ptmProteins{$varMod}}),"];\n";
	}
}
print qq
|var selectedProtID;
var isNav = (navigator.appName.indexOf("Netscape") !=-1);
function unselectProtein() {
	if (selectedProtID) {
		document.getElementById('protDetailsDIV').style.display='none';
		document.getElementById('IMG_'+selectedProtID).src="$promsPath{images}/plus.gif";
	}
	selectedProtID=null;
}
function ajaxGetProteinDetails(e,protID,bestAnaStrg) {
	var protDIV=document.getElementById('protDetailsDIV');
	if (!protDIV) { // page is not fully loaded yet
		alert('Display process is still on-going. Try again when completed.');
		return;
	}
	if (protID==selectedProtID) { // image is minus1.gif
		unselectProtein();
		return;
	}
	unselectProtein();
	selectedProtID=protID;
	document.getElementById('IMG_'+protID).src="$promsPath{images}/minus1.gif";
	var divX = (isNav)? e.pageX : event.clientX + document.body.scrollLeft; divX-=5;
	var divY = (isNav)? e.pageY : event.clientY + document.body.scrollTop; divY+=10;
	protDIV.style.left = divX + 'px';
	protDIV.style.top = divY + 'px';
	protDIV.innerHTML="<CENTER><BR><FONT class=\\"title3\\">Fetching data...</FONT><BR><IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><INPUT type=\\"button\\" value=\\" Cancel \\" class=\\"font11\\" onclick=\\"document.getElementById('protDetailsDIV').style.display='none'\\"><BR><BR>";
	protDIV.style.display='';
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	/* POST */
	var paramStrg="AJAX=getProtDetails&id_project=$projectID&protID="+protID+"&hiddenRule=$hiddenRule&pepRule=$peptideRule&pepThreshold=$pepThreshold&virtualData=$virtualData&pepSpecificity=$pepSpecificity&numGroups=$numGroups&groupPos=$groupPosStrg&anaList=$anaList&anaContext=$contextAna&bestAna="+bestAnaStrg;
	XHR.open("POST","$promsPath{cgi}/compareAnalyses.cgi",true);
	//Send the proper header information along with the request
	XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
	XHR.setRequestHeader("Content-length", paramStrg.length);
	XHR.setRequestHeader("Connection", "close");
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			protDIV.innerHTML=XHR.responseText;
		}
	}
	XHR.send(paramStrg);
}
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

function generateAutoCheckString() {
	var protForm=document.protView;
	//logic;
	var autoChkStrg=protForm.logic.value+'#';
	if (protForm.logic.value=='frequency') {
		autoChkStrg+=protForm.freqAnaGr.value+';'+protForm.freqPep.value;
		if (protForm.refPTM) { // relevant PTMs in project
			var refPtm=protForm.refPTM.value;
			if (refPtm) {
				autoChkStrg+='#'+refPtm+'#'+protForm.freqPtm.value;
			}
		}
	}
	else { // and,or
		//peptides
		for (var g=0; g<$numAnaGroups; g++) {
			if (g>0) {autoChkStrg+=';'}
			var aGr=g+1;
			autoChkStrg+=aGr+':'+protForm['pepCompSign_'+g].value+':'+protForm['numPep_'+g].value;
		}
		//PTMs
		if (protForm.refPTM) { // relevant PTMs in project
			var refPtm=protForm.refPTM.value;
			if (refPtm) {
				autoChkStrg+='#'+refPtm+'#';
				for (var g=0; g<$numAnaGroups; g++) {
					if (g>0) {autoChkStrg+=';'}
					var aGr=g+1;
					autoChkStrg+=aGr+':'+protForm['ptmCompSign_'+g].value+':'+protForm['numPtm_'+g].value;
				}
			}
		}
	}
	return autoChkStrg;
}
function selectSort(newSort) {
	parent.selItemsFrame.setComparisonAsModified(2); // 2! because stored in top.promsFrame.selectedView
	top.promsFrame.selectedView=newSort; //openProject.cgi
	parent.selItemsFrame.document.compForm.sort.value=newSort;
	var myForm=document.protView;
	myForm.exportType.value=null;
	myForm.sort.value=newSort;
	myForm.autocheck.value=generateAutoCheckString();
	myForm.action="./compareAnalyses.cgi";
	myForm.target='listFrame';
	myForm.submit();
}
function exportResult(newExportType) {
	if (newExportType=='sel' && numProtChecked==0) {
		alert('ERROR: No proteins selected!');
		return;
	}
	var myForm=document.protView;
	myForm.exportType.value=newExportType;
	myForm.action="./compareAnalyses.cgi";
	myForm.target='listFrame';
	myForm.submit();
}
function sequenceView(proteinId,analysisId){
	var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+analysisId+"&id_prot="+proteinId+"&msdata="+top.showMSdata;
	top.openProtWindow(winLocation);
}
function graphicalView() {
	if (numProtChecked==0) {
		alert('ERROR: No proteins selected!');
		return;
	}
	top.openProtWindow('');
	document.protView.largWindow.value=(document.body)? top.protWindow.document.body.clientWidth : top.protWindow.innerWidth;
	document.protView.graphView.value=top.graphView;
	document.protView.action="./graphicalView.cgi";
	document.protView.target='ProteinWindow';
	document.protView.submit();
}
function uncheckAll(checkBoxList) {
	parent.selItemsFrame.setComparisonAsModified(1);
	if (checkBoxList.length) {
		for (var i=0; i < checkBoxList.length; i++) {checkBoxList[i].checked=false;}
	}
	else {checkBoxList.checked=false;}
	updateNumChecked(-numProtChecked);
	clear2ndVennDiagram();
	hideShowUnchecked(0);
}
function manualCheck(chkStatus) {
	parent.selItemsFrame.setComparisonAsModified(1);
	var newChecked=(chkStatus)? 1 : -1;
	updateNumChecked(newChecked);
	//document.getElementById('manChkBUT').disabled=(numProtChecked > 0)? false : true;
	update2ndVennDiagram();
}
function activatePtmFiltering(refPtm) {
	parent.selItemsFrame.setComparisonAsModified(1);
	var logicPtmDisabStatus=(refPtm && document.getElementById('logic').value != 'frequency')? false : true;
	for (let a=0; a < $numAnaGroups; a++) {
		document.getElementById('ptmCompSign_'+a).disabled=logicPtmDisabStatus;
		document.getElementById('numPtm_'+a).disabled=logicPtmDisabStatus;
	}
	document.getElementById('freqPtm').disabled=(refPtm)? !logicPtmDisabStatus : true; // Warning: not defined if no relevant PTMs in project}
}
function updateComparisonStrategy(strategy) {
	var logicPepDisabStatus=(strategy=='frequency')? true : false;
	var refPTM=document.getElementById('refPTM');
	var logicPtmDisabStatus=(logicPepDisabStatus \|\| (refPTM && !refPTM.value))? true : false;
	for (let a=0; a < $numAnaGroups; a++) {
		document.getElementById('pepCompSign_'+a).disabled=logicPepDisabStatus;
		document.getElementById('numPep_'+a).disabled=logicPepDisabStatus;
		if (refPTM) {
			document.getElementById('ptmCompSign_'+a).disabled=logicPtmDisabStatus;
			document.getElementById('numPtm_'+a).disabled=logicPtmDisabStatus;
		}
	}
	document.getElementById('freqAnaGr').disabled = !logicPepDisabStatus;
	document.getElementById('freqPep').disabled = !logicPepDisabStatus;
	if (refPTM) {
		document.getElementById('freqPtm').disabled = (refPTM.value)? !logicPtmDisabStatus : true;
	}
}
function autoCheck(warn) {
	//document.getElementById('manChkBUT').disabled=true;
	var protForm=document.protView;
	var checkBoxList=protForm.chkProt;
	//Peptides
	var pepDistribList=protForm.pepDistrib;
	var pepCompSign=[];
	var numPep=[];
	for (let a=0; a < $numAnaGroups; a++) {
		pepCompSign[a]=document.getElementById('pepCompSign_'+a).value;
		numPep[a]=document.getElementById('numPep_'+a).value*1; //convert to number
	}
	var logic=document.getElementById('logic').value; //inter-group logic
	var freqAnaGr=protForm.freqAnaGr;
	var freqPep=protForm.freqPep;
	var freqPtm=protForm.freqPtm;
	
	//PTMs
	var refPTM;
	var refPtmIndex;
	var ptmDistribList;
	var ptmCompSign=[];
	var numPtm=[];
	if (protForm.ptmDistrib) {
		refPTM=document.getElementById('refPTM').value;
		refPtmIndex=document.getElementById('refPTM').selectedIndex - 1; // -=Select=- -> 0
		ptmDistribList=protForm.ptmDistrib;
		for (let a=0; a < $numAnaGroups; a++) {
			ptmCompSign[a]=document.getElementById('ptmCompSign_'+a).value;
			numPtm[a]=document.getElementById('numPtm_'+a).value*1; //convert to number
		}
	}
	//Scanning protein list
	var numMatch=0;
	var allMatches=0;
	PROT: for (let i=0; i<pepDistribList.length; i++) {
		var anaDistrib=pepDistribList[i].value.split(':');
		var allPtmDistib;
		var refPtmDistrib;
		if (refPTM) {
			allPtmDistib=ptmDistribList[i].value.split(';');
			refPtmDistrib=allPtmDistib[refPtmIndex].split(':');
		}
		if (logic=='and') {
			//Peptides
			for (let a=0; a<pepCompSign.length; a++) {
				if ((pepCompSign[a]=='>=' && anaDistrib[a]*1<numPep[a]) \|\| (pepCompSign[a]=='<=' && anaDistrib[a]*1>numPep[a])) {
					continue PROT;
				}
			}
			//PTMs
			if (refPTM) {
				for (let a=0; a<ptmCompSign.length; a++) {
					if ((ptmCompSign[a]=='>=' && refPtmDistrib[a]*1<numPtm[a]) \|\| (ptmCompSign[a]=='<=' && refPtmDistrib[a]*1>numPtm[a])) {
						continue PROT;
					}
				}
			}
		}
		else if (logic=='or') {
			var okMatch;
			for (let a=0; a<pepCompSign.length; a++) { // same length than ptmCompSign
				okMatch=true; // reset for each group
				//Peptides
				if ((pepCompSign[a]=='>=' && anaDistrib[a]*1<numPep[a]) \|\| (pepCompSign[a]=='<=' && anaDistrib[a]*1>numPep[a])) {
					okMatch=false;
				}
				//PTMs
				if (okMatch==true && refPTM) { // logic is always 'and' between pep & ptm
					if ((ptmCompSign[a]=='>=' && refPtmDistrib[a]*1<numPtm[a]) \|\| (ptmCompSign[a]=='<=' && refPtmDistrib[a]*1>numPtm[a])) {
						okMatch=false;
					}
				}
				if (okMatch==true) {break;}
			}
			if (okMatch==false) {continue PROT;}
		}
		else { // frequency
			//Peptides
			var numMatchPep=0;
			for (let a=0; a<anaDistrib.length; a++) {
				if (anaDistrib[a]*1 >= freqPep.value) {numMatchPep++;}
			}
			if (numMatchPep < freqAnaGr.value) {continue PROT;}
			//PTMs
			if (refPTM) {
				var numMatchPtm=0;
				for (let a=0; a<refPtmDistrib.length; a++) {
					if (refPtmDistrib[a]*1 >= freqPtm.value) {numMatchPtm++;}
				}
				if (numMatchPtm < freqAnaGr.value) {continue PROT;}
			}
		}
		//Result
		if (!checkBoxList[i].checked) { //okMatch &&
			numMatch++;
			checkBoxList[i].checked=true;
			var chkData=checkBoxList[i].value.split(':');
			var tr=document.getElementById('prot_'+chkData[chkData.length-1]);
			tr.className='list row_'+(numMatch % 2);
			tr.style.display=''; // make sure checked prot is visible
		}
		allMatches++;
	}
	if (numMatch) {updateNumChecked(numMatch);}
	update2ndVennDiagram();
	if (warn==1) {
		if (numMatch) {parent.selItemsFrame.setComparisonAsModified(1);}
		alert(allMatches+' proteins matched selected criteria ('+numMatch+' new)!');
	}
}
function update2ndVennDiagram() {
	if ($numAnaGroups > 5) {return;}
	clear2ndVennDiagram();
	if (numProtChecked==0) {return;}
	var protForm=document.protView;
	var checkBoxList=protForm.chkProt;
	if (!checkBoxList.length) {return;} // 1 protein only
	var pepDistribList=protForm.pepDistrib;
	var newVennValues={};
	var labelCode=['A','B','C','D','E'];
	for (var i=0; i<checkBoxList.length; i++) {
		if (!checkBoxList[i].checked) continue;
		//Peptides
		var label=null; // important
		var anaDistrib=pepDistribList[i].value.split(':');
		for (var a=0; a < anaDistrib.length; a++) {
			if (anaDistrib[a]*1 > 0) {
				if (label) {label+=labelCode[a];}
				else {label=labelCode[a];}
			}
		}
		if (newVennValues[label]) {newVennValues[label]++;}
		else {newVennValues[label]=1;}
	}
	document.getElementById('filterTextDIV').style.display='';
	document.getElementById('filterVennDIV').style.display='';
	new VennDiagram({div:'venn2DIV',size:300,groupLabels:groupLabelsObj,setValues:newVennValues});
}
function clear2ndVennDiagram() {
	if ($numAnaGroups > 5) {return;}
	document.getElementById('filterTextDIV').style.display='none';
	document.getElementById('filterVennDIV').style.display='none';
	document.getElementById('venn2DIV').innerHTML='';
}
function hideShowUnchecked(warn) {
	var myForm=document.protView;
	if (numProtChecked==0 && myForm.hideShowButton.value=='Hide Unchecked') {
		if (warn) alert('There are no proteins checked!');
		return;
	}
	if (myForm.hideShowStatus.value=='none') {
		myForm.hideShowStatus.value=''; // show all
		myForm.hideShowButton.value='Hide Unchecked';
	}
	else { // '' ~ block
		myForm.hideShowStatus.value='none'; // show checked only
		myForm.hideShowButton.value='Show Unchecked';
	}
	var checkBoxList=myForm.chkProt;
	if (checkBoxList.length) {
		var numVisible=0;
		for (var i=0; i < checkBoxList.length; i++) {
			var chkData=checkBoxList[i].value.split(':');
			var tr=document.getElementById('prot_'+chkData[chkData.length-1]);
			if (myForm.hideShowStatus.value=='none') { // show checked only
				if (checkBoxList[i].checked) {
					numVisible++;
					tr.className='list row_'+(numVisible % 2);
					tr.style.display='';
				}
				else {
					tr.style.display='none';
				}
			}
			else { // show all
				numVisible++;
				tr.className='list row_'+(numVisible % 2);
				tr.style.display='';
			}
		}
	}
	else if (!checkBoxList.checked) { // only 1 prot in comparison
		var chkData=checkBoxList.value.split(':');
		document.getElementById('prot_'+chkData[chkData.length-1]).style.display=myForm.hideShowStatus.value;
	}
}
function updateNumChecked(newChecked) {
	numProtChecked+=newChecked;
	document.getElementById('numProtChecked').innerHTML=numProtChecked;
}
var numProtChecked=$numProtChecked;
var groupLabelsObj;
function whenWindowLoaded() {
|;
	my $vennDisplayStrg;
	if ($numAnaGroups <= 5) {
		$vennDisplayStrg='';
		print "\tdocument.getElementById('vennTABLE').style.display='';\n\tgroupLabelsObj={";
		foreach my $g (0..$#analysisGroups) {
			print chr(65+$g),":'$groupLabels[$g][0]'";
			print ',' if $g < $#analysisGroups;
		}
		print "};\n\tnew VennDiagram({div:'vennDIV',size:300,clickAction:updateProteinList,groupLabels:groupLabelsObj,setValues:{";
		my $numSets=scalar keys %vennDiagram;
		my $count=0;
		foreach my $set (sort keys %vennDiagram) {
			print "$set:$vennDiagram{$set}";
			$count++;
			print ',' if $count < $numSets;
		}
		print "}});\n";
	}
	else {$vennDisplayStrg='display:none';}
	print qq
|	$onloadStrg
}
function updateProteinList(set) { // called when a set is selected in Venn diagram
	/* Clear all selections */
	resetSelection();
	/* Auto selection */
	if (set=='all') {return;}
	else if (set.match(':full')) {
		var setInfo=set.split(':');
		var a=set.charCodeAt(0)-65; // ascii code
		//document.getElementById('numPep_'+a).selectedIndex=1; // numPep>=1
		var numPepLC=document.getElementById('numPep_'+a);
		if (numPepLC.options.length > 1) {numPepLC.selectedIndex=1;} // numPep>=1 unless no protein in gr because of pep threshold
	}
	else {
		var labelCode=['A','B','C','D','E'];
		for (var a=0; a < $numAnaGroups; a++) {
			if (set.match(labelCode[a])) {
				//document.getElementById('numPep_'+a).selectedIndex=1;} // numPep>=1
				var numPepLC=document.getElementById('numPep_'+a);
				if (numPepLC.options.length > 1) {numPepLC.selectedIndex=1;} // numPep>=1 unless no protein in gr because of pep threshold
			}
			else {
				document.getElementById('pepCompSign_'+a).selectedIndex=1; // numPep<=0
			}
		}
	}
	autoCheck(0);
	hideShowUnchecked(0);
}
function resetSelection() {
	var myForm=document.protView;
	if (myForm.refPTM) { // not defined if no PTMs
		myForm.refPTM.selectedIndex=0;
		activatePtmFiltering('');
	}
	document.getElementById('logic').selectedIndex=0; // logic = 'and'
	updateComparisonStrategy('and');
	uncheckAll(myForm.chkProt);
	//hideShowUnchecked(0);
	for (var a=0; a < $numAnaGroups; a++) {
		document.getElementById('pepCompSign_'+a).selectedIndex=0;
		document.getElementById('numPep_'+a).selectedIndex=0;
	}
	if (myForm.ptmDistrib) {
		for (var a=0; a < $numAnaGroups; a++) {
			document.getElementById('ptmCompSign_'+a).selectedIndex=0;
			document.getElementById('numPtm_'+a).selectedIndex=0;
		}
	}
}
document.getElementById('waitDIV').style.display='none';
</SCRIPT>

<TABLE id="vennTABLE" border=0 align="center" cellspacing=0 cellpadding=0 style="display:none"><TR>
<TH><DIV style="text-align:center;$vennDisplayStrg">
	<FONT class="title2">Complete set</FONT>&nbsp;&nbsp;<INPUT type="button" value="Export as image" style="font-weight:bold;font-size:10px" onclick="exportSVGtoImg('vennDIV','VennDiagram','./exportSVG.cgi')"/>
	<DIV id="vennDIV" style="text-align:center;padding:3px;"></DIV>
</DIV></TH>
<TH><DIV id="filterTextDIV" style="text-align:center;display:none"><FONT class="title2">&gt;&gt;&gt;Filtering&gt;&gt;&gt;</FONT></DIV></TH>
<TH><DIV id="filterVennDIV" style="text-align:center;display:none">
	<FONT class="title2">Filtered set</FONT>&nbsp;&nbsp;<INPUT type="button" value="Export as image" style="font-weight:bold;font-size:10px" onclick="exportSVGtoImg('venn2DIV','VennDiagram','./exportSVG.cgi')"/>
	<DIV id="venn2DIV" style="text-align:center;padding:3px"></DIV>
</DIV></TH>
</TR></TABLE>
|;

####>Printing Store form<####
&printStoreForm;

####>Printing protein list<####
print qq
|<FORM name="protView" method="post" style="margin:0">
<!--parameters for Export -->
<INPUT type="hidden" name="ACT" value="list">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="parentItem" value="$parentItem">
<INPUT type="hidden" name="numGroups" value="$numGroups">
<INPUT type="hidden" name="catFilter" value="$catFilter">
<INPUT type="hidden" name="catFilterRule" value="$catFilterRule">
<INPUT type="hidden" name="pepRule" value="$peptideRule">
<INPUT type="hidden" name="pepThreshold" value="$pepThreshold">
<INPUT type="hidden" name="noHidden" value="$noHidden">
<INPUT type="hidden" name="hiddenRule" value="$hiddenRule">
<INPUT type="hidden" name="comparisonID" value="$comparisonID">
<INPUT type="hidden" name="virtualData" value="$virtualData">
<INPUT type="hidden" name="autocheck" value="">
<INPUT type="hidden" name="hideShowStatus" value="$hideShowStatus">
<INPUT type="hidden" name="sort" value="$sortOrder">
<INPUT type="hidden" name="anaList" value="$anaList">
<INPUT type="hidden" name="exportType" value="">
<!-- No need for peptide modification filters -->
<!--parameters for Peptide Distribution -->
<INPUT type="hidden" name="what" value="checkbox">
<INPUT type="hidden" name="item" value="$item">
<INPUT type="hidden" name="id_item" value="$itemID">
<INPUT type="hidden" name="largWindow">
<INPUT type="hidden" name="graphView">

<TABLE border=0 cellspacing=0 cellpadding=2><TR><TD><DIV id="saveProtDIV" style="display:none;"></DIV></TD></TR></TABLE>
|;

###>Action table<###
my ($actGrWidth,$grColSpan)=(scalar keys %projectVarMods>=1)? (170,2) : (80,1);
my $disPtmStrg=' disabled';
print qq
|$listContextStrg
<TABLE border=0 cellpadding=0 cellspacing=0>
<TR><TH class="rbBorder" bgcolor=$color2 class="font11" rowspan=2 colspan=2 valign=top nowrap>&nbsp;<FONT class="title3">Auto-check :</FONT>|;
if (scalar keys %projectVarMods>=1) {
	print "<BR>&nbsp;PTM:<SELECT name=\"refPTM\" id=\"refPTM\" class=\"font11\" onchange=\"activatePtmFiltering(this.value)\"><OPTION value=\"\">-=Select=-</OPTION>";
	my $defSelStrg=((scalar keys %projectVarMods)==1)? ' selected' : '';
	foreach my $varMod (sort{$a cmp $b} keys %projectVarMods) {
		my $selStrg=($selectedRefPTM && $varMod eq $selectedRefPTM)? ' selected' : $defSelStrg;
		print "<OPTION value=\"$varMod\"$selStrg>$projectVarMods{$varMod}[0]</OPTION>";
		$disPtmStrg='' if $selStrg;
		$selStrg='';
	}
	print "</SELECT>&nbsp;";
}
print "<BR>(<SPAN id=\"numProtChecked\">$numProtChecked</SPAN> checked)</TH>\n";
foreach my $g (0..$#analysisGroups) {
	print "<TH class=\"rbBorder\" width=$actGrWidth bgcolor=$color2 colspan=$grColSpan nowrap><FONT onmouseover=\"popup('<B>$groupLabels[$g][1]</B>')\" onmouseout=\"popout()\">$groupLabels[$g][0]</FONT></TH>\n";
}
my $disSaveCompStrg=($projectStatus > 0 || $projectAccess eq 'guest' || $call eq '1v1_prot')? 'disabled' : ''; # project was ended/archived
my $disSaveProtStrg=($projectStatus > 0 || $projectAccess eq 'guest')? 'disabled' : ''; # project was ended/archived
my ($selOrStatus,$selFreqStatus)=($filterLogic eq 'or')? (' selected','') : ($filterLogic eq 'frequency')? ('',' selected') : ('','');
my $thTableColSpan=($numAnaGroups==2 && !$listComparison)? 4 : 2;
my $freqStrg='';
foreach my $g (0..$#analysisGroups) {
	my $gNum=$g+1;
	$freqStrg.="<OPTION value=\"$gNum\"";
	$freqStrg.=' selected' if ($autocheckParams{'FREQ_ANAGR'} && $autocheckParams{'FREQ_ANAGR'}==$gNum);
	$freqStrg.=">$gNum</OPTION>";
}
my (%allNumPeptides,%allNumPepStrg);
foreach my $g (0..$#analysisGroups) {
	foreach my $numPep (keys %{$groupPeptides[$g]}) {$allNumPeptides{PEP}{$numPep}=1;}
	if (scalar keys %projectVarMods>=1) {
		if ($ptmOccurences{$g}) {
			foreach my $occ (keys %{$ptmOccurences{$g}}) {$allNumPeptides{PTM}{$occ}=1;}
		}
	}
}
foreach my $type (keys %allNumPeptides) {
	$allNumPepStrg{$type}='';
	foreach my $num (sort{$a<=>$b} keys %{$allNumPeptides{$type}}) {
		$allNumPepStrg{$type}.="<OPTION value=\"$num\"";
		$allNumPepStrg{$type}.=' selected' if ($autocheckParams{'FREQ_'.$type} && $autocheckParams{'FREQ_'.$type}==$num);
		$allNumPepStrg{$type}.=">$num</OPTION>";
	}
}
$allNumPepStrg{PTM}='' unless $allNumPepStrg{PTM};
print qq
|<TH class="rbBorder" rowspan=2 colspan=$thTableColSpan height=100% valign=top><TABLE border=0 cellpadding=0 cellspacing=0 height=100%><TR>
	<TH class="rBorder" bgcolor=$color2 nowrap valign="top">&nbsp;Frequency:<SELECT name="freqAnaGr" id="freqAnaGr" class="font11" disabled>$freqStrg</SELECT>/$numAnaGroups&nbsp;
<BR><FONT class="font11">&nbsp;Num. Pept.&ge;<SELECT name="freqPep" id="freqPep" class="font11" disabled>$allNumPepStrg{PEP}</SELECT>&nbsp;
|;
print qq |<BR>&nbsp;Num. PTM&ge;<SELECT name="freqPtm" id="freqPtm" class="font11" disabled><OPTION value="0">0</OPTION>$allNumPepStrg{PTM}</SELECT>&nbsp;| if scalar keys %projectVarMods; #$allNumPepStrg{PTM};
print qq
|</FONT>
	</TH>
	<TH class="rBorder" bgcolor=$color2 nowrap>
&nbsp;Strategy<SUP onmouseover="popup('Select<BR>&nbsp;&nbsp;&nbsp;-A <B>logic</B> between Group filters: <B>And</B> (intersection) or <B>Or</B> (union)<BR>&nbsp;&nbsp;&nbsp;or<BR>&nbsp;&nbsp;&nbsp;-a global <B>frequency</B> of detection')\" onmouseout=\"popout()\">*</SUP>:
<SELECT name="logic" id="logic" class="font11" onchange="updateComparisonStrategy(this.value)"><OPTION value="and">And</OPTION><OPTION value="or"$selOrStatus>Or</OPTION><OPTION value="frequency"$selFreqStatus>Frequency</OPTION></SELECT>&nbsp;<BR>
&nbsp;<INPUT type="button" value="Check Matching" class="font11" style="width:110px" onclick="autoCheck(1)">&nbsp;<INPUT type="button" value="Uncheck All" class="font11" style="width:85px" onclick="uncheckAll(document.protView.chkProt)">&nbsp;
<BR>&nbsp;<INPUT type="button" name="hideShowButton" value="$hideShowStrg" class="font11" style="width:110px" onclick="hideShowUnchecked(1)">&nbsp;
	</TH>
	<TH nowrap valign=top><INPUT type="button" value="Save Comparison..." style="width:150px" onclick="showSaveCompForm('show')" $disSaveCompStrg><BR>
<INPUT type="button" id="saveFormBUTTON" value="Save Selected..." style="width:150px" onclick="ajaxManageSaveProteins('getThemes');" $disSaveProtStrg></TH>
	<TH nowrap valign=top><INPUT type="button" value="Export Selected" onclick="exportResult('sel')" style="width:120px"><BR><INPUT type="button" value="Export All" onclick="exportResult('all')" style="width:120px"></TD>
	<TH nowrap valign=top><INPUT type="button" value="Peptide Distribution" style="width:150px" onclick="graphicalView()"></TH>
</TR></TABLE></TH></TR>
<TR>
|;
foreach my $g (0..$#analysisGroups) {
	my $aGr=$g+1;
	# Peptides
	my ($selSupPep,$selInfPep)=($autocheckParams{'PEP'} && $autocheckParams{'PEP'}{$aGr}[0] eq '<=')? ('',' selected') : (' selected','');
	print "<TH class=\"rbBorder\" bgcolor=$color2 nowrap>";
	print "<FONT class=\"font11\">Peptides</FONT><BR>" if scalar keys %projectVarMods>=1;
	print "<SELECT name=\"pepCompSign_$g\" id=\"pepCompSign_$g\" class=\"font11\"><OPTION value=\">=\"$selSupPep>&ge;</OPTION><OPTION value=\"<=\"$selInfPep>&le;</OPTION></SELECT>";
	print "<SELECT name=\"numPep_$g\" id=\"numPep_$g\" class=\"font11\"><OPTION value=\"0\">0</OPTION>";
	foreach my $numPep (sort{$a<=>$b} keys %{$groupPeptides[$g]}) {
		print "<OPTION value=\"$numPep\"";
		print ' selected' if ($autocheckParams{'PEP'} && $autocheckParams{'PEP'}{$aGr}[1]==$numPep);
		print ">$numPep</OPTION>";
	}
	print "</SELECT></TH>\n";
	# PTMs
	if (scalar keys %projectVarMods>=1) {
		my ($selSupPtm,$selInfPtm)=($autocheckParams{'PTM'} && $autocheckParams{'PTM'}{$aGr}[0] eq '<=')? ('',' selected') : (' selected','');
		print "<TH class=\"rbBorder\" bgcolor=$color2 nowrap><FONT class=\"font11\">PTM<BR></FONT>";
		print "<SELECT name=\"ptmCompSign_$g\" id=\"ptmCompSign_$g\" class=\"font11\"$disPtmStrg><OPTION value=\">=\"$selSupPtm>&ge;</OPTION><OPTION value=\"<=\"$selInfPtm>&le;</OPTION></SELECT>";
		print "<SELECT name=\"numPtm_$g\" id=\"numPtm_$g\" class=\"font11\"$disPtmStrg><OPTION value=\"0\">0</OPTION>";
		if ($ptmOccurences{$g}) {
			foreach my $occ (sort{$a<=>$b} keys %{$ptmOccurences{$g}}) {
				print "<OPTION value=\"$occ\"";
				print ' selected' if ($autocheckParams{'PTM'} && $autocheckParams{'PTM'}{$aGr}[1]==$occ);
				print ">$occ</OPTION>";
			}
		}
		print "</SELECT></TH>\n";
	}
}
print qq
|</TR>
<TR><TH style="font-size:4px">&nbsp;</TH></TR>
|;

###>List head table<###
# Printing according to selected sort
my $proteinText=scalar(keys %protAnnotation).' Proteins';
my $massText='MW <FONT class="small">(kDa)</FONT>';
my $organismText='Species';
my $peptideText='Peptides';
my $absDeltaText='Abs. &#916';
my $relDeltaText='Rel. &#916 %';
foreach my $g (0..$#analysisGroups) {
	$groupLabels[$g][0]="<FONT color=#DD0000>$groupLabels[$g][0]</FONT>" if $sortOrder eq 'group:'.($g+1);
}
if ($sortOrder eq 'protein') {
	$proteinText="<FONT color=#DD0000>$proteinText</FONT>";
}
elsif ($sortOrder eq 'mass') {
	$massText="<FONT color=#DD0000>$massText</FONT>";
}
elsif ($sortOrder eq 'organism') {
	$organismText="<FONT color=#DD0000>$organismText</FONT>";
}
elsif ($sortOrder eq 'peptide'){
	$peptideText="<FONT color=#DD0000>$peptideText</FONT>";
}
elsif ($sortOrder eq 'absDelta'){
	$absDeltaText="<FONT color=#DD0000>$absDeltaText</FONT>";
}
elsif ($sortOrder eq 'relDelta'){
	$relDeltaText="<FONT color=#DD0000>$relDeltaText</FONT>";
}
my $pepStrg;
if ($numGroups || $listComparison) {
	$pepStrg=($peptideRule eq 'all')? "&ge;$pepThreshold in Group" : ($peptideRule eq 'anr')? "&ge;$pepThreshold distinct in Group" : ($peptideRule eq 'ba')? "&ge;$pepThreshold in best Analysis" : "&ge;$pepThreshold distinct in best Analysis";
}
else {
	$pepStrg=($peptideRule=~/nr/)? "&ge;$pepThreshold distinct in Analysis" : "&ge;$pepThreshold in Analysis";
}
$peptideText.=" <FONT class=\"small\">($pepStrg)</FONT>";
my ($pepColSize,$ptmColSize)=(scalar keys %projectVarMods>=1 && ($numGroups || $listComparison))? (60,60) : (scalar keys %projectVarMods>=1)? (50,50) : ($numGroups || $listComparison)? (120,0) : (100,0);
my $grColSizeH=$pepColSize+$ptmColSize-4;
my $grColSizeR=$pepColSize+$ptmColSize;
my $numPepCols=($numAnaGroups==2 && !$listComparison)? ($numAnaGroups*$grColSpan)+2 : $numAnaGroups*$grColSpan; # +2 for (abs/rel) Delta
#my $tableWidth=929+($grColSizeR*$numAnaGroups);
#$tableWidth+=144 if ($numAnaGroups==2 && !$listComparison);
my $ptmTitleStrg='';
if (scalar keys %projectVarMods>=1) {
	my $ptmInfoStrg='';
	foreach my $varMod (sort{$a cmp $b} keys %projectVarMods) {
		$ptmInfoStrg.="<BR><FONT class=\\'popPTM $varMod\\'>$projectVarMods{$varMod}[1]</FONT>: $projectVarMods{$varMod}[0]";
	}
	$ptmTitleStrg=" - <A href=\"javascript:void(null)\" onmouseover=\"popup('<B>Post-Trans. Modifications:</B>$ptmInfoStrg')\" onmouseout=\"popout()\">PTMs</A>";
}
#>Drawing list head table
#my $manualChkVis=($numAnaGroups > 5)? 'display:none' : '';
#<TABLE border=0 cellpadding=0 cellspacing=0 width=$tableWidth>
print qq
|<TR class="row_0">
<TH class="rbBorder" rowspan=2><A href="javascript:selectSort('protein')" onmouseover="popup('Click to sort proteins by <B>ascending name</B>.')" onmouseout="popout()">$proteinText</A>&nbsp;<!--<INPUT type="button" id="manChkBUT" value="Update Checked" class="font11" style="width:110px;\$manualChkVis" onclick="update2ndVennDiagram()" disabled>&nbsp;--></TH>
<TH class="rbBorder" rowspan=2>&nbsp;Gene&nbsp;</TH>
<TH class="rbBorder" colspan=$numPepCols><A href="javascript:selectSort('peptide')" onmouseover="popup('Click to sort proteins by <B>best distribution</B>.')" onmouseout="popout()">$peptideText</A>$ptmTitleStrg</TH>
<TH class="rbBorder" width=75 rowspan=2><A href="javascript:selectSort('mass')" onmouseover="popup('Click to sort proteins by <B>decreasing mass</B>.')" onmouseout="popout()">$massText</A></TH>
<TH class="bBorder" width=650 rowspan=2>Description - <A href="javascript:selectSort('organism')" onmouseover="popup('Click to sort proteins by <B>ascending species</B>.')" onmouseout="popout()">$organismText</A></TH>
</TR>
<TR class="row_0">
|;
foreach my $g (0..$#analysisGroups) {
	my $numAnaStrg=($numGroups)? scalar @{$analysisGroups[$g]}.' items:<BR>' : '';
	print "<TH class=\"rbBorder\" width=$grColSizeH colspan=$grColSpan><A href=\"javascript:selectSort('group:",$g+1,"')\" onmouseover=\"popup('<B>$numAnaStrg$groupLabels[$g][1]</B><BR>Click to sort proteins by <B>decreasing number of matching peptides</B>.')\" onmouseout=\"popout()\">$groupLabels[$g][0]</A></TH>\n";
}
if ($numAnaGroups==2 && !$listComparison) {
	print qq
|<TH class="rbBorder" width=70><A href=\"javascript:selectSort('absDelta')\" onmouseover=\"popup('<B>Absolute delta:</B> #2-#1<BR>Click to sort proteins by <B>decreasing absolute delta</B>.')\" onmouseout=\"popout()\">$absDeltaText</A></TH>
<TH class="rbBorder" width=70><A href=\"javascript:selectSort('relDelta')\" onmouseover=\"popup('<B>Relative delta:</B> (#2-#1)/(#1+#2)%<BR>Click to sort proteins by <B>decreasing relative delta</B>.')\" onmouseout=\"popout()\">$relDeltaText</A></TH>
|;
}
print "</TR>\n";

my $numProt=0;
foreach my $protID (sort{&sortProt($sortOrder)} keys %protAnnotation) {
	my $chkStatus=($selectedProteins{$protID})? ' checked' : '';
	#my $displayStatus=($chkStatus)? 'block' : $hideShowStatus;
	my $displayStatusStrg=($chkStatus)? '' : " style=\"display:$hideShowStatus\"";
	my $trClass='row_'.(++$numProt % 2);
	#print "<TABLE id=\"prot_$protID\" cellspacing=0 cellpadding=0 width=$tableWidth border=0 style=\"display:$displayStatus\"><TR class=\"list $trClass\" valign=\"top\">\n";
	print "<TR id=\"prot_$protID\" class=\"list $trClass\" valign=\"top\"$displayStatusStrg>\n";
	#>Protein
	my $bestAnaStrg='';
	my %distinctAna;
	foreach my $g (0..$#analysisGroups) {
		foreach my $codeID (@{$analysisGroups[$g]}) {
			my ($type,$typeID)=split(':',$codeID); # A or C
			my @localAnaList=($type eq 'A')? ($typeID) : @parentItemAna;
			foreach my $anaID (@localAnaList) {
				$distinctAna{$anaID}=1 if ($protAnalyses{$protID}{$anaID} && $protAnalyses{$protID}{$anaID}[3]);# protein peptides are listed
			}
		}
		if ($peptideRule=~/ba/) {
			$bestAnaStrg.=':' if $bestAnaStrg ne '';
			$bestAnaStrg.=($bestAnalysis{$protID}{$g})? $bestAnalysis{$protID}{$g} : 0;
		}
	}
	my $chkBoxValue=join(':',keys %distinctAna,$protID);
	my $distAnaStrg=join(',',keys %distinctAna);
	my ($protTag,$classStrg)=($bestVisibility{$protID}[1])? ('TH','') : ('TD',' class="hiddenProt"'); # $classStrg is important for proper recording to custom Lists!!!
	print qq
|<TD width=252 nowrap>
<INPUT type="checkbox" name="chkProt"$classStrg value="$chkBoxValue" onclick="manualCheck(this.checked)"$chkStatus><A class="$protTag" href="javascript:sequenceView($protID,'$distAnaStrg')">$protAnnotation{$protID}[0]</A><IMG id="IMG_$protID" src="$promsPath{images}/plus.gif" align="top" onclick="ajaxGetProteinDetails(event,$protID,'$bestAnaStrg')">
</TD>
|;
	#>Gene strg
	my $masterProtID=$protAnnotation{$protID}[4];
	if ($masterProtID) {
		if (scalar @{$masterProteins{$masterProtID}} > 1) {
			print "<TH>&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-".(join('<BR>&nbsp;&nbsp;-',@{$masterProteins{$masterProtID}}[1..$#{$masterProteins{$masterProtID}}]))."</B>')\" onmouseout=\"popout()\">$masterProteins{$masterProtID}[0]</A>&nbsp;</TH>";
		}
		elsif ($masterProteins{$masterProtID}[0]) {print "<TH>&nbsp;$masterProteins{$masterProtID}[0]&nbsp;</TH>";}
		else {print "<TH>&nbsp;-&nbsp;</TH>";} # no gene mapped
	}
	else {print "<TH>&nbsp;-&nbsp;</TH>";} # no mapping at all

	#>(Groups of) Analyses/Lists
	my (@pepDistrib,%ptmDistrib);
	foreach my $g (0..$#analysisGroups) {
		my $pepTag=($proteinList{$protID}[$g][1])? 'TH' : 'TD'; # protein global(best) visibility
		#>Converting subgroup data to string
		my $numAnaExistPep=0;
		my $pepAnaDistribStrg='';
		my $itemAnaDistribStrg='';
		my $visInAnaTag='TD'; # Prot is visible in ana distribution
		my $hiddenTag='';
		if (($numGroups || $listInGroup{$g}) && $peptideAnaDistrib{$protID}[$g]) { # prot exists in group
			$pepAnaDistribStrg='<B>Distr.:</B> ';
			my $itemAnaScan=0; # if C only scan parent item analyses once
			my $codeIndex=0;
			foreach my $codeID (@{$analysisGroups[$g]}) {
				$pepAnaDistribStrg.='-' if $codeIndex;
				my ($type,$typeID)=split(':',$codeID); # A or C
				if ($type eq 'A') {
					$numAnaExistPep++ if $peptideAnaDistrib{$protID}[$g][$codeIndex];
					#$pepAnaDistribStrg.=($protAnalyses{$protID}{$typeID}[1])? "<B>$peptideAnaDistrib{$protID}[$g][$codeIndex]</B>" : $peptideAnaDistrib{$protID}[$g][$codeIndex];
					if ($protAnalyses{$protID}{$typeID}[1]) {
						$pepAnaDistribStrg.="<B>$peptideAnaDistrib{$protID}[$g][$codeIndex]</B>";
						$visInAnaTag='TH';
					}
					else {$pepAnaDistribStrg.=$peptideAnaDistrib{$protID}[$g][$codeIndex];}
					$hiddenTag='-' if ($hiddenRule==1 && $proteinList{$protID}[$g][1] && $peptideAnaDistrib{$protID}[$g][$codeIndex] && !$protAnalyses{$protID}{$typeID}[1]); # prot is vis in group but hidden in some analyses
				}
				else { # C
					if ($peptideAnaDistrib{$protID}[$g][$codeIndex]) { # protein exists in analysis/list   $peptideAnaDistrib{$protID}[$g][$codeIndex]
						$numAnaExistPep++;
						if ($itemAnaScan) { # ana list already scanned
							$pepAnaDistribStrg.='<B>L</B>'; # at least 2 Lists in same group
							$codeIndex++;
							next;
						}
						my $anaIndex=0;
						foreach my $anaID (@parentItemAna) {
							$itemAnaDistribStrg.='-' if $anaIndex;
							#if ($anaOnLine==12) {$pepAnaDistribStrg.='<BR>'; $anaOnLine=0;}
							#$numAnaExistPeptide{$anaID}=1 if $peptideAnaDistrib{$protID}[$g][$codeIndex][$anaIndex];
		#$itemAnaDistribStrg.=($protAnalyses{$protID}{$anaID}[1])? "<B>$peptideAnaDistrib{$protID}[$g][$codeIndex][$anaIndex]</B>" : $peptideAnaDistrib{$protID}[$g][$codeIndex][$anaIndex];
		if ($protAnalyses{$protID}{$anaID}[1]) {
			$itemAnaDistribStrg.="<B>$peptideAnaDistrib{$protID}[$g][$codeIndex][$anaIndex]</B>";
			$visInAnaTag='TH';
		}
		else {$itemAnaDistribStrg.=$peptideAnaDistrib{$protID}[$g][$codeIndex][$anaIndex];}
							$hiddenTag='-' if ($hiddenRule==1 && $proteinList{$protID}[$g][1] && $peptideAnaDistrib{$protID}[$g][$codeIndex][$anaIndex] && !$protAnalyses{$protID}{$anaID}[1]); # prot is vis in group but hidden in some analyses
							$anaIndex++;
							#$anaOnLine++;
						}
						if ($listInGroup{$g} > 1) {$pepAnaDistribStrg.='<B>L</B>';}
						else {$pepAnaDistribStrg.="[$itemAnaDistribStrg]";}
						$itemAnaScan=1;
					}
					else {$pepAnaDistribStrg.='0';}
				}
				$codeIndex++;
			}
			$pepAnaDistribStrg.="<BR><B>L</B>=$itemAnaDistribStrg" if ($listInGroup{$g} && $listInGroup{$g} > 1 && $itemAnaDistribStrg);
			$pepAnaDistribStrg.='<BR><B>-:</B> Protein is hidden in some Analyses' if $hiddenTag;
		}
		my ($fontStrg1,$fontStrg2)=($proteinList{$protID}[$g][2]==1)? ('<FONT color="gray">','</FONT>') : ('','');
		my $numAnaExistStrg=(($numGroups || $listInGroup{$g}) && $numAnaExistPep)? "<SUP><A href=\"javascript:void(null)\" onmouseover=\"popup('$pepAnaDistribStrg')\" onmouseout=\"popout()\">$fontStrg1$numAnaExistPep$hiddenTag$fontStrg2</A></SUP>" : ''; # " <FONT class=\"small\">($numAnaExistPep)</FONT>"
		push @pepDistrib,$proteinList{$protID}[$g][0];
		#my %numProtPtmGr;
		if ($proteinList{$protID}[$g][0]) {
			#peptides
			#print "<TD class=\"$pepTag\" width=$pepColSize align=\"center\">$fontStrg1$proteinList{$protID}[$g][0]$fontStrg2$numAnaExistStrg</TD>\n";
			print "<TD width=$pepColSize align=\"center\">$fontStrg1<FONT class=\"$pepTag\">$proteinList{$protID}[$g][0]</FONT>$fontStrg2<FONT class=\"$visInAnaTag\">$numAnaExistStrg</FONT></TD>\n";

			#PTMs
			next unless scalar keys %projectVarMods>=1;
			my $varModText='-';
			if ($protPtmGroups{$protID} && $protPtmGroups{$protID}{$g}) {
				$varModText='';
				foreach my $varMod (sort{$a cmp $b} keys %{$protPtmGroups{$protID}{$g}}) {
					$varModText.="<FONT class=\"$projectVarMods{$varMod}[2]\">";
					$varModText.="<FONT class=\"font11\">$protPtmGroups{$protID}{$g}{$varMod}</FONT>" if $protPtmGroups{$protID}{$g}{$varMod} > 1;
					$varModText.="$projectVarMods{$varMod}[1]</FONT>";
					if ($numGroups || $listInGroup{$g}) {
						my %numAnaExistPtm;
						my $ptmAnaDistribStrg="<B>Distr. </B><FONT class=\\'popPTM $varMod\\'>$projectVarMods{$varMod}[1]</FONT>: ";
						my $itemPtmAnaDistribStrg='';
						my $itemAnaScan=0; # if C only scan parent item analyses once
						#my $anaOnLine=0;
						my $firstCode=1;
						foreach my $codeID (@{$analysisGroups[$g]}) {
							$ptmAnaDistribStrg.='-' unless $firstCode;
							my ($type,$typeID)=split(':',$codeID); # A or C
							my @localAnaList;
							if ($type eq 'A') {
								@localAnaList=($typeID);
							}
							else { # C
								if ($itemAnaScan) { # ana list already scanned
									$ptmAnaDistribStrg.='<B>L</B>'; # at least 2 Lists in same group
									next;
								}
								@localAnaList=@parentItemAna;
								#$ptmAnaDistribStrg.='[';
							}
							my $firstAna=1;
							foreach my $anaID (@localAnaList) {
								my $numProtPtmGr=0;
								if ($ptmAnalyses{$protID}{$g}{$anaID} && $ptmAnalyses{$protID}{$g}{$anaID}{$varMod}) {
									$numAnaExistPtm{$anaID}=1;
									$numProtPtmGr=($peptideRule=~/nr/)? scalar keys %{$ptmAnalyses{$protID}{$g}{$anaID}{$varMod}{'NR'}} : $ptmAnalyses{$protID}{$g}{$anaID}{$varMod}{'ALL'};
								}
								#if ($anaOnLine==12) {$ptmAnaDistribStrg.='<BR>'; $anaOnLine=0;}
								if ($type eq 'A') {
									$ptmAnaDistribStrg.=($protAnalyses{$protID}{$anaID}[1])? "<B>$numProtPtmGr</B>" : $numProtPtmGr;
								}
								else {
									$itemPtmAnaDistribStrg.='-' unless $firstAna;
									$itemPtmAnaDistribStrg.=($protAnalyses{$protID}{$anaID}[1])? "<B>$numProtPtmGr</B>" : $numProtPtmGr;
								}
								$firstAna=0;
								#$anaOnLine++;
							}
							if ($type eq 'C') {
								if ($listInGroup{$g} > 1) {$ptmAnaDistribStrg.='<B>L</B>';}
								else {$ptmAnaDistribStrg.="[$itemPtmAnaDistribStrg]";}
								$itemAnaScan=1;
							}
							$firstCode=0;
						}
						#foreach my $anaID (@{$analysisGroups[$g]}) {
						#	my $numProtPtmGr=0;
						#	if ($ptmAnalyses{$protID}{$g}{$anaID} && $ptmAnalyses{$protID}{$g}{$anaID}{$varMod}) {
						#		$numAnaExistPtm++;
						#		$numProtPtmGr=($peptideRule=~/nr/)? scalar keys %{$ptmAnalyses{$protID}{$g}{$anaID}{$varMod}{'NR'}} : $ptmAnalyses{$protID}{$g}{$anaID}{$varMod}{'ALL'};
						#	}
						#	if ($anaOnLine==12) {$ptmAnaDistribStrg.='<BR>'; $anaOnLine=0;}
						#	$ptmAnaDistribStrg.='-' unless $firstAna;
						#	$ptmAnaDistribStrg.=($protAnalyses{$protID}{$anaID}[1])? "<B>$numProtPtmGr</B>" : $numProtPtmGr;
						#	$firstAna=0;
						#	$anaOnLine++;
						#}
						$ptmAnaDistribStrg.="<BR><B>L</B>=$itemPtmAnaDistribStrg" if ($listInGroup{$g} && $listInGroup{$g} > 1 && $itemPtmAnaDistribStrg);
						$varModText.="<SUP><A class=\"$varMod\" href=\"javascript:void(null)\" onmouseover=\"popup('$ptmAnaDistribStrg')\" onmouseout=\"popout()\">".(scalar keys %numAnaExistPtm)."</A></SUP>";
					}
				}
				$varModText='-' unless $varModText;
			}
			print "<TD class=\"$pepTag\" width=$pepColSize align=\"center\">$varModText</TD>\n";
		}
		else {print "<TD colspan=$grColSpan class=\"$pepTag\" width=$grColSizeR align=\"center\">-</TD>\n";}
		if (scalar keys %projectVarMods>=1) {
			foreach my $varMod (keys %projectVarMods) {
				if ($protPtmGroups{$protID} && $protPtmGroups{$protID}{$g} && $protPtmGroups{$protID}{$g}{$varMod}) {push @{$ptmDistrib{$varMod}},$protPtmGroups{$protID}{$g}{$varMod};}
				else {push @{$ptmDistrib{$varMod}},0;}
			}
		}
	}
	print "<INPUT type=\"hidden\" name=\"pepDistrib\" value=\"",join(':',@pepDistrib),"\"/>\n";
	if (scalar keys %projectVarMods>=1) {
		my @protPtmDistrib;
		foreach my $varMod (sort{$a cmp $b} keys %projectVarMods) {
			#push @protPtmDistrib,"$projectVarMods{$varMod}[2]=".join(':',@{$ptmDistrib{$varMod}});
			push @protPtmDistrib,join(':',@{$ptmDistrib{$varMod}});
		}
		print "<INPUT type=\"hidden\" name=\"ptmDistrib\" value=\"",join(';',@protPtmDistrib),"\"/>\n";
	}
	#>Relative delta
	if ($numAnaGroups==2 && !$listComparison) {
		my $relDelta=sprintf "%.1f",$relativeDelta{$protID};
		$relDelta=~s/\.0\Z//;
		print "<TD class=\"$protTag\" width=72 align=\"center\">$absoluteDelta{$protID}</TD><TD class=\"$protTag\" width=72 align=\"center\">$relDelta</TD>\n";
	}
	#>Mass
	my $mass=($protAnnotation{$protID}[2]==0)? '-' : $protAnnotation{$protID}[2];
	print "<TD class=\"$protTag\" width=77 align=\"right\">$mass&nbsp;&nbsp;</TD>\n";
	#>Description - Species
	print "<TD width=600>$protAnnotation{$protID}[1] <FONT class=\"org\">$protAnnotation{$protID}[3]</FONT></TD>\n";
	#print "</TR></TABLE>\n";
	print "</TR>\n";
	#$bgColor=($bgColor eq $color1)? $color2 : $color1;
}
print qq
|</TABLE>
<B><I>End of list.</I></B>
</FORM>
<BR><BR>
<DIV id="protDetailsDIV" style="z-index:200;background-color:#FFFFFF;border:solid 3px #999999;padding:5px;position:absolute;display:none;"></DIV> <!--filter:alpha(opacity=80);opacity:0.8;-->
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">
setPopup();
whenWindowLoaded();
</SCRIPT>
</BODY>
</HTML>
|;


sub printStoreForm {
	print qq
|<SCRIPT type="text/javascript">
function showSaveCompForm(action) {
	document.getElementById('storeComp').style.display=(action=='show')? '' : 'none';
}
function displayTargetComparison(tgtID) {
	if (tgtID==0) {
		document.getElementById('newNameSPAN').style.display='';
		document.getElementById('newCompComments').style.display='';
		document.getElementById('compComments').style.display='none';
	}
	else {
		document.getElementById('newNameSPAN').style.display='none';
		document.getElementById('newCompComments').style.display='none';
		document.getElementById('compComments').style.display='';
	}
}
// Saving comparison
function checkSaveComparison() {
	var myForm=document.storeForm;
	if (myForm.targetComp.value==0 && !myForm.newCompName.value) {
		alert('ERROR: Provide a name for Comparison');
		return;
	}
	if (myForm.targetComp.value != 0 && !confirm('Overwrite current Comparison?')) {return;}
	if ('$compType'=='full_prot') {myForm.autocheck.value=generateAutoCheckString();}
	myForm.submit();
}
</SCRIPT>
<FORM name="storeForm" method="post" style="margin:0">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="parentItem" value="$parentItem">
<INPUT type="hidden" name="compType" value="$compType">
<INPUT type="hidden" name="ACT" value="saveComp">
<INPUT type="hidden" name="numGroups" value="$numGroups">
<INPUT type="hidden" name="catFilter" value="$catFilter">
<INPUT type="hidden" name="catFilterRule" value="$catFilterRule">
<INPUT type="hidden" name="anaList" value="$anaList">
<INPUT type="hidden" name="sort" value="$sortOrder">
<INPUT type="hidden" name="pepRule" value="$peptideRule">
<INPUT type="hidden" name="pepThreshold" value="$pepThreshold">
<INPUT type="hidden" name="virtualData" value="$virtualData">
<INPUT type="hidden" name="noHidden" value="$noHidden">
<INPUT type="hidden" name="hiddenRule" value="$hiddenRule">
<INPUT type="hidden" name="noMissCut" value="$noMissCut">
|;
	print "<INPUT type=\"hidden\" name=\"modifList\" value=\"",join(',',@modifList),"\">\n";
	if ($modifFilters{active}) {
		foreach my $filter ('exclude','restrict','ignore') {
			foreach my $modID (sort{$a<=>$b} keys %{$modifFilters{$filter}}) {print "<INPUT type=\"hidden\" name=\"filterMod_$modID\" value=\"$filter\">\n";}
		}
	}
	print qq
|<INPUT type="hidden" name="delocPhospho" value="$delocPhospho">
<INPUT type="hidden" name="autocheck" value="=">
<!--INPUT type="hidden" name="protID" value="-"-->

<TABLE border=0 id="storeComp" style="display:none">
<TR><TD bgcolor=$color2><TABLE cellpadding=2>
	<TR><TH align=left colspan=2><FONT style="font-size:15px;">&nbsp;&nbsp;Saving Comparison:</FONT></TH></TR>
	<TR><TH align=right width=100>Target :</TH><TD bgcolor=$color1><TABLE cellpadding=0 cellspacing=0><TR><TD><SELECT name="targetComp" onchange="displayTargetComparison(this.value)"><OPTION value="0">New</OPTION>
|;
	my ($visOldSPAN,$visNewSPAN)=('none','');
	if ($comparisonID) {
		print "<OPTION value=\"$comparisonID\" selected>$compName</OPTION>\n";
		$visOldSPAN='';
		$visNewSPAN='none';
	}
	print qq
|</SELECT></TD><TD><SPAN id="newNameSPAN" style="display:$visNewSPAN">&nbsp;<B>Name:</B>&nbsp;<INPUT type="text" name="newCompName" value="" style="width:450px"/></SPAN></TD></TR></TABLE>
</TD></TR>
	<TR><TH align=right valign=top>Comments :</TH><TD bgcolor=$color1><TEXTAREA name="newCompComments" id="newCompComments" rows="2" cols="85" style="display:$visNewSPAN"></TEXTAREA><TEXTAREA name="compComments" id="compComments" rows="2" cols="85" style="display:$visOldSPAN">$compComments</TEXTAREA></TD>
	<TR><TH colspan=2><INPUT type="button" name="saveComp" value=" Save " onclick="checkSaveComparison()">&nbsp;<INPUT type="button" value="Cancel" onclick="showSaveCompForm('hide')"></TR>
</TABLE></TD></TR>
</TABLE>
</FORM>
|;

}

sub updateProgress {
	my ($text,$action)=@_;
	$text='' unless $text;
	$action='=' unless $action; # '=' or '+='
	print qq
|<SCRIPT type="text/javascript">document.getElementById('waitSPAN').innerHTML$action'$text';</SCRIPT>
|;
}

############################################################
#######################<FRAMES>#############################
############################################################
sub generateFrames {
	print header(-charset=>'utf-8');
	# warningsToBrowser(1);
	print qq
|<HTML>
<FRAMESET rows="205,*" scrolling=no border=0>
	<FRAME name="selItemsFrame" src="./compareAnalyses.cgi?ACT=selAna&id_project=$projectID&parentItem=$parentItem&sort=$sortOrder&compType=$compType&comparisonID=$comparisonID" scrolling=no>
	<FRAME name="listFrame" src="$promsPath{html}/nothing.html">
</FRAMESET>
</HTML>
|;
	exit;
}


##########################################################################
#######################<Items selection Form>#############################
##########################################################################
sub selectAnalyses {

	####<Connecting to the database>####
	my $dbh=&promsConfig::dbConnect;

	#################################
	####<Fetching item hierarchy>####
	#################################

	my %allPostTransModifs=&promsMod::getVariableModifications($dbh); # usable PTMs
	my $sthGetPM=$dbh->prepare("SELECT ID_MODIFICATION FROM PROJECT_MODIFICATION WHERE ID_PROJECT=$projectID");
	$sthGetPM->execute;
	my %projectVarMods;
	while (my ($modID)=$sthGetPM->fetchrow_array) {
		$projectVarMods{$modID}=$allPostTransModifs{$modID};
	}
	$sthGetPM->finish;

	####<List of Comparisons>####
	my %comparisonList;
	my $sthC=$dbh->prepare("SELECT ID_COMPARISON,NAME,NUM_GROUPS,PEPTIDE_PARAMS,ID_CATEGORY,CAT_EXCLUSION,SORT_ORDER,AUTOCHECK_PARAMS FROM COMPARISON WHERE ID_PROJECT=$projectID");
	$sthC->execute;
	while (my($compID,$name,$numGroups,$pepParams,$catFilterID,$catExclusion,$sort,$autoChk)=$sthC->fetchrow_array) {
		$pepParams.=';virtualData=0' if $pepParams !~ /virtualData/; # compatibility with comparisons before update 1.4.0
		#if ($cpType ne '1v1_pep') {
		$pepParams.=';pepThreshold=1' if $pepParams !~ /pepThreshold/; # compatibility with comparisons before update 1.7.0
		#}
		$pepParams.=";delocPhospho=0" if $pepParams !~ /delocPhospho/; # compatibility with comparisons before update 2.0.0
		$pepParams.=";noMissCut=0" if $pepParams !~ /noMissCut/; # compatibility with comparisons before update 2.0.6
		$pepParams.=";pepSpecificity=all" if $pepParams !~ /pepSpecificity/; # compatibility with comparisons before update 2.2.0
		$comparisonList{$compID}{'NAME'}=$name;
#$comparisonList{$compID}{'COMP_TYPE'}="'$cpType'"; # TODO: add field TYPE to COMPARISON table
		$comparisonList{$compID}{'NUM_GROUPS'}=$numGroups;
		$comparisonList{$compID}{'PEP_PARAMS'}="'$pepParams'";
		$comparisonList{$compID}{'FILTER'}=$catFilterID || 0;
		$comparisonList{$compID}{'FILTER_RULE'}=($catFilterID && $catExclusion)? "'exclude'" : "'restrict'";
		$comparisonList{$compID}{'SORT'}=($sort)? "'$sort'" : '';
		#$comparisonList{$compID}{'AUTOCHECK'}=($autoChk)? "'$autoChk'" : "''";
		if ($autoChk) {
			my $ptmStrg=(split('#',$autoChk))[2];
			if ($ptmStrg && $ptmStrg=~/\D/) { # PSI name NOT ptmID (compatibility with old string-based version of PTM managment)
				foreach my $ptmID (keys %allPostTransModifs) {
					if ($ptmStrg eq $allPostTransModifs{$ptmID}[0]) {
						$autoChk=~s/#$ptmStrg#/#$ptmID#/; # convert to ID
						last;
					}
				}
			}
		}
		else {$autoChk='';}
		$comparisonList{$compID}{'AUTOCHECK'}="'$autoChk'";
	}
	$sthC->finish;

	####<List of Analyses>####
	my @queryList;
	push @queryList,qq |SELECT DISTINCT CONCAT('EXPERIMENT:',E.ID_EXPERIMENT),E.NAME,CONCAT('GEL2D:',G.ID_GEL2D),G.NAME
							FROM EXPERIMENT E,GEL2D G,SPOT SP,SAMPLE S,ANALYSIS A
							WHERE E.ID_PROJECT=$projectID AND E.ID_EXPERIMENT=G.ID_EXPERIMENT
							AND G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT
							AND S.ID_SAMPLE=A.ID_SAMPLE AND A.VALID_STATUS>=1
							ORDER BY E.DISPLAY_POS ASC,G.DISPLAY_POS ASC|; #,SP.NAME ASC
	push @queryList,qq |SELECT DISTINCT CONCAT('EXPERIMENT:',E.ID_EXPERIMENT),E.NAME,CONCAT('SAMPLE:',S.ID_SAMPLE),S.NAME
							FROM EXPERIMENT E,SAMPLE S,ANALYSIS A
							WHERE E.ID_PROJECT=$projectID AND E.ID_EXPERIMENT=S.ID_EXPERIMENT
							AND S.ID_SAMPLE=A.ID_SAMPLE AND S.ID_SPOT IS NULL AND A.VALID_STATUS>=1
							ORDER BY E.DISPLAY_POS ASC,S.DISPLAY_POS ASC|;
	my (@itemList,%usedItems);
	foreach my $query (@queryList) {
		my $sthItem=$dbh->prepare($query);
		$sthItem->execute;
		while (my @itemInfo=$sthItem->fetchrow_array) {
			for (my $i=0; $i<$#itemInfo; $i+=2) {
				next if $usedItems{$itemInfo[$i]};
				my $labelString='';
				for (my $j=0; $j<=$i-2; $j+=2) {
					$labelString.=' > ' if $labelString;
					$labelString.=&promsMod::shortenName($itemInfo[$j+1],21);
				}
				$labelString.=' > ' if $labelString;
				$labelString.=$itemInfo[$i+1];
				push @itemList,[$itemInfo[$i],$labelString];
				$usedItems{$itemInfo[$i]}=1;
			}
		}
		$sthItem->finish;
	}
	my $usedParent;
	if ($item eq 'project') {
		$usedParent=$itemList[0][0];
	}
	elsif ($item eq 'spot') {
		($usedParent)=$dbh->selectrow_array("SELECT CONCAT('GEL2D:',G.ID_GEL2D) FROM GEL2D G,SPOT SP WHERE G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=$itemID");
	}
	elsif ($item eq 'analysis') {
		($usedParent)=$dbh->selectrow_array("SELECT CONCAT('SAMPLE:',S.ID_SAMPLE) FROM SAMPLE S,ANALYSIS A WHERE S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS=$itemID AND S.ID_SPOT IS NULL");
		unless ($usedParent) { # GEL
			my $query=qq |SELECT CONCAT('GEL2D:',G.ID_GEL2D) FROM GEL2D G,SPOT SP,SAMPLE S,ANALYSIS A
								WHERE G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT
								AND S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS=$itemID|;
			($usedParent)=$dbh->selectrow_array($query);
		}
	}
	else {$usedParent=uc($parentItem);}

	####<List of Modifications used in Analyses in Project for filter>####
	my %modifications;
	my $sthProjMod=$dbh->prepare("SELECT DISTINCT ID_MODIFICATION FROM ANALYSIS_MODIFICATION AM,ANALYSIS A,SAMPLE S,EXPERIMENT E WHERE AM.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_EXPERIMENT=E.ID_EXPERIMENT AND MODIF_TYPE='V' AND ID_PROJECT=$projectID");
	my $sthMod=$dbh->prepare("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES FROM MODIFICATION WHERE ID_MODIFICATION=?");
	$sthProjMod->execute;
	while (my ($modID)=$sthProjMod->fetchrow_array) {
		$sthMod->execute($modID);
		my ($psiName,$interName,$syn)=$sthMod->fetchrow_array;
		$modifications{$modID}=$psiName || $interName;
		unless ($modifications{$modID}) {
			$syn=~s/^##//; $syn=~s/##$//; $syn=~s/##/, /;
			$modifications{$modID}=$syn;
		}
		$modifications{$modID}=&promsMod::resize($modifications{$modID},50);
	}
	$sthMod->finish;
	$sthProjMod->finish;

	my ($phosphoID)=$dbh->selectrow_array('SELECT ID_MODIFICATION FROM MODIFICATION WHERE UNIMOD_ACC=21');

	####<List of Classification/Category for filter>####
	my %classificationList=&promsMod::getListClass($dbh,$projectID) if $compType ne '1v1_pep';

	$dbh->disconnect;

	##############
	####<HTML>####
	##############
	print header(-charset=>'utf-8');
	warningsToBrowser(1); # writes PERL warnings as HTML comments <!-- --> must be activated after header
	my ($lightColor1,$darkColor2)=&promsConfig::getRowColors;
	my $ITEM=uc($item);
	print qq
|<HEAD>
<TITLE>Compare Multiple Analyses</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE>
.formButton {font-weight:bold;font-size:10px}
.formSelect {font-weight:bold;font-size:11px}
.formText {font-weight:bold;font-size:11px}
.popup {z-index:100;max-height:140px;min-width:200px;overflow-y:auto;background-color:#FFFFFF;border:solid 3px #999999;padding:5px;position:absolute;display:none;box-shadow:10px 10px 20px #808080;}
</STYLE>
<SCRIPT type="text/javascript">
|;
	&promsMod::popupInfo();
	print qq
|//----> Saved comparisons
var projectID=$projectID; // needed in case edition/deletion of comparison in popup window
var comparisonList=new Object();
comparisonList['0']=new Object(); // no comparisons saved
//comparisonList['0'].COMP_TYPE='$compType';
comparisonList['0'].NUM_GROUPS=$numGroups; // 0..n
comparisonList['0'].PEP_PARAMS='noHidden=1;hiddenRule=2;pepRule=$peptideRule;pepThreshold=$pepThreshold;virtualData=$virtualData;noMissCut=$noMissCut;delocPhospho=$delocPhospho;pepSpecificity=$pepSpecificity';
comparisonList['0'].FILTER_RULE='$catFilterRule';
comparisonList['0'].FILTER=0; // 0/ID List
comparisonList['0'].SORT='$sortOrder'; // only needed for full_prot
comparisonList['0'].AUTOCHECK='$autochkStrg';
|;
	my $selCompName=($comparisonID)? $comparisonList{$comparisonID}{NAME} : 'New';
	my $selCompIndex=0;
	my $count=0;
	foreach my $compID (sort{$a<=>$b} keys %comparisonList) {
		$count++;
		$selCompIndex=$count if $compID==$comparisonID;
		print "comparisonList['$compID']=new Object();\n";
		foreach my $attrib ('NUM_GROUPS','PEP_PARAMS','FILTER','FILTER_RULE','SORT','AUTOCHECK') { # ,'VIRTUAL_DATA'
			print "comparisonList['$compID'].$attrib=$comparisonList{$compID}{$attrib};\n";
			#my $attribValue=($comparisonList{$compID}{$attrib}=~/\D/ || $comparisonList{$compID}{$attrib} eq '')? "'$comparisonList{$compID}{$attrib}'" : $comparisonList{$compID}{$attrib};
			#print "comparisonList['$compID'].$attrib=$attribValue;\n";
		}
	}
	print qq
|//<-----
//------>Modifications specificity
var modifSpecificity={
|;
	my $numModifs=scalar keys %projectVarMods;
	my $countM=0;
	foreach my $modID (keys %projectVarMods) {
		print "\t$modID:[";
		my @resList=split(',',$projectVarMods{$modID}[4]);
		my $numRes=scalar @resList;
		my $countR=0;
		foreach my $specif (@resList) {
			next if $specif=~/[=\*]/; # irrelevant for "Any-term" modifs
			#my $res=($specif eq '-')? 'Protein N-term' : ($specif eq '=')? 'Any N-term' : ($specif eq '+')? 'Protein C-term' : ($specif eq '*')? 'Any C-term' : $specif;
			my $res=($specif eq '-')? 'Protein N-term' : ($specif eq '+')? 'Protein C-term' : $specif;
			print "['$specif','$res']";
			print ',' if ++$countR < $numRes;
		}
		print ']';
		print ',' if ++$countM < $numModifs;
		print "\n";
	}
	print qq
|};
//<-----
var selCompName='$selCompName';
var selCompIndex=$selCompIndex; // index of selected comparison in document.compForm.comparisonID!!!
var modCompType=0;
function selectComparisonType(newCompType) {
	var myForm=document.compForm;
	if (newCompType=='1v1_pep') {
		//Peptide vs Protein options
		document.getElementById('pepCompOptions').style.display='';
		document.getElementById('protCompOptions').style.display='none';
		disableListComparison();
/*
		//Disable all lists
		for (var i=0; i<myForm.anaParent.options.length; i++) {
			if (myForm.anaParent.options[i].value.match('CLASSIFICATION:')) {
				if (i==myForm.anaParent.selectedIndex) { // deselect Classification opt
					myForm.anaParent.selectedIndex=1;
					ajaxGetAnalysesList(myForm.anaParent.value);
				}
				myForm.anaParent.options[i].disabled=true;
			}
		}
		//Remove lists from selected items
		for (var i=myForm.usedAna.options.length-1; i>=0; i--) {
			if (myForm.usedAna.options[i].value.match(/^C:/)) {
				myForm.usedAna.removeChild(myForm.usedAna.options[i]);
			}
		}
		//Remove lists from existing groups
		var modGroup=false;
		for (var g=0; g<analysisGroups.length; g++) {
			for (var i=analysisGroups[g].length-2; i>=0; i-=2) {
				if (analysisGroups[g][i+1].match(/^C:/)) { // remove 2 indexes text,value
					analysisGroups[g].splice(i,2);
					modGroup=true;
				}
			}
		}
		if (modGroup) {alert('Warning: Custom lists were removed from all groups');}
*/
	}
	else {
		//Peptide vs Protein options
		document.getElementById('pepCompOptions').style.display='none';
		document.getElementById('protCompOptions').style.display='';
		//Activate all lists
		for (var i=0; i<myForm.anaParent.options.length; i++) {
			if (myForm.anaParent.options[i].value.match('CLASSIFICATION:')) {myForm.anaParent.options[i].disabled=false;}
		}
		if (newCompType.match('modif_sites')) {
			document.getElementById('pepNumSPAN').style.display='none';
			var compInfo=newCompType.split(':');
			document.getElementById('modifResSPAN').style.display='';
			disableListComparison();
			/* Update SELECT according to selected modif */
			myForm.siteRes.options.length=1;
			for (var i=0; i<modifSpecificity[compInfo[1]].length; i++) {
				myForm.siteRes.options[i+1]=new Option(modifSpecificity[compInfo[1]][i][1],modifSpecificity[compInfo[1]][i][0]);
			}
		}
		else {
			document.getElementById('pepNumSPAN').style.display='';
			document.getElementById('modifResSPAN').style.display='none';
		}
	}
	showModifFilterForm(false);
}
function disableListComparison() {
	var myForm=document.compForm;
	//Disable all lists
	for (var i=0; i<myForm.anaParent.options.length; i++) {
		if (myForm.anaParent.options[i].value.match('CLASSIFICATION:')) {
			if (i==myForm.anaParent.selectedIndex) { // deselect Classification opt
				myForm.anaParent.selectedIndex=1;
				ajaxGetAnalysesList(myForm.anaParent.value);
			}
			myForm.anaParent.options[i].disabled=true;
		}
	}
	var listIsUsed=false;
	//Remove lists from selected items
	for (var i=myForm.usedAna.options.length-1; i>=0; i--) {
		if (myForm.usedAna.options[i].value.match(/^C:/)) {
			myForm.usedAna.removeChild(myForm.usedAna.options[i]);
			listIsUsed=true;
		}
	}
	//Remove lists from existing groups
	for (var g=0; g<analysisGroups.length; g++) {
		for (var i=analysisGroups[g].length-2; i>=0; i-=2) {
			if (analysisGroups[g][i+1].match(/^C:/)) { // remove 2 indexes text,value
				analysisGroups[g].splice(i,2);
				listIsUsed=true;
			}
		}
	}
	if (listIsUsed) {alert('Warning: Custom lists were removed from selected items');}
}
function selectComparison(compID,onChange) {
	if (onChange && compID==0) { // comparison was changed to 'New' => Reset resultFrame
		modCompType=0;
		parent.listFrame.location='$promsPath{html}/nothing.html';
	}
	var myForm=document.compForm;
	var selComp=myForm.comparisonID;
	selComp.options[selCompIndex].text=selCompName; // !corrects name in case modified '*'
	selCompName=selComp.options[selComp.selectedIndex].text;
	selCompIndex=selComp.selectedIndex;

	var newGrNumber=comparisonList[compID].NUM_GROUPS;
	myForm.numGroups.selectedIndex=(newGrNumber==0)? 0 : (newGrNumber*1)-1;
	myForm.currentGroup.length=1; // restart from scratch
	myForm.usedAna.length=0;
	analysisGroups.length=0;
	updateNumGroups(newGrNumber);

	/*** Peptide options ***/
	//Reset all modif fiters to 'allow'
	var allMods=(myForm.modifList.value+'').split(',');
	for (var m=0; m<allMods.length; m++) {myForm['filterMod_'+allMods[m]].selectedIndex=0;}
	myForm.restrictModifLogic.selectedIndex=0; // restrictModifLogic
	//Loop through all peptide options
	var pepParams=comparisonList[compID].PEP_PARAMS.split(';');
	for (var p=0; p<pepParams.length; p++) {
		if (!pepParams[p]) continue; // in case empty (;;)
		var param=pepParams[p].split('=');
		if (param[0].match('noHidden\|virtualData\|noMissCut\|delocPhospho')) { // checkboxes
			if (!myForm[param[0]]) continue; // delocPhospho defined only if Phospho is used in project
			myForm[param[0]].checked=(param[1]=='1')? true : false;
		}
		else if (param[0]=='modifFilters') {
			var filters=param[1].split(',');
			for (var f=0; f<filters.length; f++) { // 1 select per modif
				var filterMods=filters[f].split(':');
				var selIndex=(filterMods[0]=='exclude')? 1 : (filterMods[0]=='restrict')? 2 : 3; // cannot be 0 (allow)
				var modList=filterMods[1].split('.');
				for (var m=0; m<modList.length; m++) {
					if (myForm['filterMod_'+modList[m]]) { // Make sure this mod is still used in project (just to be safe)
						myForm['filterMod_'+modList[m]].selectedIndex=selIndex;
					}
				}
				if (filterMods[2] && filterMods[2]=='any') { // restrictModifLogic
					myForm.restrictModifLogic.selectedIndex=1;
				}
			}
		}
		else { // hiddenRule, pepRule, pepThreshold, pepSpecificity (select)
			myForm[param[0]].selectedIndex=0; // default
			for (var i=0; i<myForm[param[0]].options.length; i++) {
				if (param[1] && param[1]==myForm[param[0]].options[i].value) {
					myForm[param[0]].selectedIndex=i;
					break;
				}
			}
		}
	}
	//Sort Order & Autocheck (used after submission)
	myForm.sort.value=comparisonList[compID].SORT;
	myForm.autocheck.value=comparisonList[compID].AUTOCHECK;

	//Category filter
	for (var i=0; i<myForm.catFilter.options.length; i++) {
		if (myForm.catFilter.options[i].value==comparisonList[compID].FILTER) {
			myForm.catFilter.selectedIndex=i;
			break;
		}
	}
	// Category filter rule
	myForm.catFilterRule.selectedIndex=(comparisonList[compID].FILTER_RULE=='exclude')? 1 : 0;

	//AJAX-dependent processes (last because of server delay)
	if (compID != 0) {
		ajaxUpdateCompGroupAnalyses(newGrNumber,compID);
	}
	else {
		selectedGroup=0;
	}
}
function showModifFilterForm(status) {
	var modifDiv=document.getElementById('modifDIV');
	if (!status) {
		modifDiv.style.display='none';
		return;
	}
	// Show
	modifDiv.style.display='block'; // must be 'block'!
	var coord=getElementPosition(document.getElementById('showModif'));
	var top=15;
	var left=Math.max(0,coord[1]);
	modifDiv.style.left = left+'px';
	modifDiv.style.top = top+'px'
}
function checkCatFilterRule(catID) {
	if (!catID) document.compForm.catFilterRule.selectedIndex=0; // no filter -> set rule to 'restrict'
}
function setComparisonAsModified(modType) {
	if (selCompIndex != 0) {
		document.compForm.comparisonID.options[selCompIndex].text=selCompName+'*';
		if (modCompType < modType) {modCompType=modType;}
	}
}
function displayAllComparisons() {
	var compWin=window.open("$promsPath{cgi}/manageComparisons.cgi?id_project=$projectID",'ComparisonsWindow','width=900,height=700,location=no,resizable=yes,scrollbars=yes');
	compWin.focus();
}
var selectedGroup=0;
var analysisGroups=new Array();
function updateNumGroups(newGrNumber) {
	var myForm=document.compForm;
	var selPepRule=(myForm.pepRule)? myForm.pepRule.value : 0; // defined only for full_prot
	if (newGrNumber > myForm.currentGroup.length) {
		if (myForm.currentGroup.length==1) { // never 0
			myForm.currentGroup[0]=new Option(1,1);
			analysisGroups[0]=new Array();
			selectedGroup=1;
			if (myForm.pepRule) { // defined only for full_prot
				myForm.pepRule.length=0;
				myForm.pepRule[0]=new Option('in Group','all');
				myForm.pepRule[1]=new Option('distinct in Group','anr');
				myForm.pepRule[2]=new Option('in best Analysis','ba');
				myForm.pepRule[3]=new Option('distinct in best Analysis','banr');
				myForm.pepRule.selectedIndex=(selPepRule=='ba')? 0 : 1;
			}
		}
		for (var i=myForm.currentGroup.length; i<newGrNumber; i++) {
			myForm.currentGroup[i]=new Option(i+1,i+1);
			analysisGroups[i]=new Array();
		}
	}
	else {
		myForm.currentGroup.length=newGrNumber;
		analysisGroups.length=newGrNumber;
		if (newGrNumber==0) {
			myForm.currentGroup[0]=new Option('-',0);
			analysisGroups=new Array();
			selectedGroup=0;
			if (myForm.pepRule) { // defined only for full_prot
				myForm.pepRule.length=0;
				myForm.pepRule[0]=new Option('in Analysis','ba');
				myForm.pepRule[1]=new Option('distinct in Analysis','banr');
				myForm.pepRule.selectedIndex=(selPepRule.match('nr'))? 1 : 0;
			}
		}
		else if (newGrNumber<selectedGroup) {displayGroup(1);}
	}
}
function setCurrentGroup(curGroup) {
	var myForm=document.compForm;
	//Storing analyses in selectedGroup
	analysisGroups[selectedGroup-1].length=0;
	for (var i=0; i<myForm.usedAna.options.length; i++) {
		analysisGroups[selectedGroup-1].push(myForm.usedAna.options[i].text);
		analysisGroups[selectedGroup-1].push(myForm.usedAna.options[i].value); // 2x indexes
	}
	displayGroup(curGroup);
}
function displayGroup(selGroup) { //Display selGroup analyses
	//selGroup*=1;
	var myForm=document.compForm;
	myForm.usedAna.options.length=0;
	var anaCount=0;
	for (var i=0; i<analysisGroups[selGroup-1].length; i+=2) { // 2x indexes
		myForm.usedAna[anaCount]=new Option(analysisGroups[selGroup-1][i],analysisGroups[selGroup-1][i+1]);
		anaCount++;
	}
	selectedGroup=selGroup;
}
function updateAnalyses(action) {
	setComparisonAsModified(2);
	var allSelect=document.compForm.allAna;
	var usedSelect=document.compForm.usedAna;
	var selected=false;
	if (action=='add') {
		for (var i=0; i<usedSelect.options.length; i++) {usedSelect.options[i].selected=false;}
		//Processing & adding parent label
		var selAnaParent=document.compForm.anaParent;
		var parentLabel=selAnaParent.options[selAnaParent.options.selectedIndex].text;
		var labelList=parentLabel.split(' > ');
		var lastLabel=labelList[labelList.length-1];
		if (lastLabel.length>21) {
			var begStrg=lastLabel.substr(0,9);
			var endStrg=lastLabel.substr(lastLabel.length-9,9);
			labelList[labelList.length-1]=begStrg+'...'+endStrg;
		}
		parentLabel=labelList.join(' > ');
		ALL:for (var i=0; i<allSelect.length; i++) {
			if (allSelect.options[i].selected) {
				selected=true;
				//Check if not already used
				for (var j=0; j<usedSelect.length; j++) {
					if (usedSelect.options[j].value==allSelect.options[i].value) continue ALL;
				}
				usedSelect.options[usedSelect.length]=new Option(parentLabel+' > '+allSelect.options[i].text,allSelect.options[i].value);
				usedSelect.options[usedSelect.length-1].selected=true;
				allSelect.options[i].selected=false;
			}
		}
		usedSelect.focus();
	}
	else { //remove
		var keepOptions=new Array();
		for (var i=0; i<usedSelect.length; i++) {
			if (!usedSelect.options[i].selected) {
				keepOptions.push(usedSelect.options[i].text);
				keepOptions.push(usedSelect.options[i].value);
			}
			else {selected=true;}
		}
		usedSelect.length=0;
		for (var i=0; i<keepOptions.length; i+=2) {
			usedSelect.options[usedSelect.length]=new Option(keepOptions[i],keepOptions[i+1]);
		}
	}
	if (!selected) {alert('No Analysis selected');}
}
function moveAnalyses(direction) {
	setComparisonAsModified(2);
	var usedSelect=document.compForm.usedAna;
	var selectOK=false;
	if (direction=='up') {
		for (var i=0; i<usedSelect.length; i++) {
			selectOK=switchAnalyses(usedSelect,i,-1,selectOK);
		}
	}
	else {
		for (var i=usedSelect.length-1; i>=0; i--) {
			selectOK=switchAnalyses(usedSelect,i,1,selectOK);
		}
	}
	usedSelect.focus();
	if (!selectOK) {alert('No Analysis selected');}
}
function clearAnalyses() {
	setComparisonAsModified(2);
	document.compForm.usedAna.length=0;
}
function switchAnalyses(usedSelect,i,delta,selectOK) {
	if (usedSelect.options[i].selected) {
		selectOK=true;
		if (i+delta<0 \|\| i+delta==usedSelect.length \|\| usedSelect.options[i+delta].selected) return selectOK;
		var prevAnaID=usedSelect.options[i+delta].value;
		var prevText=usedSelect.options[i+delta].text;
		usedSelect.options[i+delta]=new Option(usedSelect.options[i].text,usedSelect.options[i].value);
		usedSelect.options[i]=new Option(prevText,prevAnaID);
		usedSelect.options[i+delta].selected=true;
	}
	return selectOK;
}
var OneVsOneWin;
function checkSelectedAnalyses(myForm) {
	if (selectedGroup==0) { // no groups
		if (myForm.usedAna.length<=1) {
			alert('Select at least 2 Analyses or Lists.');
			return false;
		}
		myForm.anaList.value='';
		for (var i=0; i<myForm.usedAna.length; i++) {
			myForm.anaList.value+=myForm.usedAna[i].value;
			if (i<myForm.usedAna.length-1) myForm.anaList.value+='.'; // :
		}
	}
	else { //groups
		setCurrentGroup(selectedGroup); // recorded currently displayed group
		myForm.anaList.value='';
		for (var g=0; g<analysisGroups.length; g++) {
			if (analysisGroups[g].length==0) {
				alert('No analyses selected for Group #'+(g+1));
				return false;
			}
			if (g>0) myForm.anaList.value+=',';
			//myForm.anaList.value+=g+'=';
			for (var i=1; i<analysisGroups[g].length; i+=2) { // 2x indexes
				if (i>1) myForm.anaList.value+='.'; // :
				myForm.anaList.value+=analysisGroups[g][i];
			}
		}
	}
	if (modCompType==1) { // modified by listFrame
		myForm.comparisonID.options[selCompIndex].text=selCompName; // removes '*' from name
	}
	/* Popup (or not!) window */
	if (myForm.compType.value=='1v1_prot') {
		if (OneVsOneWin && !OneVsOneWin.closed) {OneVsOneWin.focus();}
		else {OneVsOneWin=window.open(null,'OneVsOneWindow','width=1000,height=500,location=no,resizable=yes,scrollbars=yes');}
		myForm.target='OneVsOneWindow';
		parent.listFrame.location='$promsPath{html}/nothing.html';
	}
	else {
		myForm.target='listFrame';
	}
	return true;
}
function autoSelectAnalyses() {
	var myForm=document.compForm;
	var searchString=myForm.selectText.value;
	if (!searchString) {
		alert('No search string typed.');
		return;
	}
	var bs=String.fromCharCode(92);
	var unsafe=bs+".+*?[^]\$(){}=!<>\|";
	for (var i=0;i<unsafe.length;++i) {
		searchString=searchString.replace(new RegExp("\\\\"+unsafe.charAt(i),"g"),bs+unsafe.charAt(i));
	}
	//Searching option text
	var selStatus=(myForm.autoAction.value=='select')? true : false;
	for (var i=0; i<myForm.allAna.length; i++) {
		if (myForm.allAna.options[i].text.match(searchString)) {
			myForm.allAna.options[i].selected=selStatus;
		}
	}
}
function clearSelection() {
	var myForm=document.compForm;
	for (var i=0; i<myForm.allAna.length; i++) {
		myForm.allAna.options[i].selected=false;
	}
}
// AJAX --->
function ajaxGetAnalysesList(branchID) {
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/compareAnalyses.cgi?AJAX=getAnaList&branchID="+branchID,true);
	//XHR.open("POST","$promsPath{cgi}/compareAnalyses.cgi",true);
	//XHR.send("AJAX=getAnaList&branchID="+branchID);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			listAllAnalyses(XHR.responseText);
		}
	}
	XHR.send(null);
}
function ajaxUpdateRestrict(compID) {
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	var catFilterID=(compID)? comparisonList[compID].FILTER : document.compForm.catFilter.value;
	var asyncFlag=(compID)? false : true;
	XHR.open("GET","$promsPath{cgi}/compareAnalyses.cgi?AJAX=ajaxRestrictList&projectID=$projectID&catFilter="+catFilterID,asyncFlag); //true
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			document.getElementById('restrictSPAN').innerHTML=XHR.responseText;
		}
	}
	XHR.send(null);
}
function ajaxUpdateCompGroupAnalyses(newGrNumber,compID) {
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/compareAnalyses.cgi?AJAX=getAnaText&compID="+compID,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			storeGroupAnalyses(newGrNumber,XHR.responseText);
		}
	}
	XHR.send(null);
}
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
function listAllAnalyses(anaText) {
	var anaList=anaText.split('\\n');
	var selAllAna=document.compForm.allAna;
	selAllAna.length=0;
	var o=0;
	for (var i=0; i<anaList.length; i++) {
		var anaInfo=anaList[i].split('#:#');
		if (anaInfo.length==2) { // skip warning & last empty line
			selAllAna.options[o]=new Option(anaInfo[1],anaInfo[0]);
			selAllAna.options[o].selected=true;
			o++;
		}
	}
}
function storeGroupAnalyses(newGrNumber,anaText) {
	var anaList=anaText.split('\\n');
	var myForm=document.compForm;
	if (newGrNumber==0) { // no groups
		selectedGroup=0;
		for (var i=0; i<anaList.length-1; i++) { // -1 because last line is empty
			var grData=anaList[i].split('#:#'); // (gr,codeID,codeText)
			myForm.usedAna[i]=new Option(grData[2],grData[1]);
		}
	}
	else {
		for (var i=0; i<anaList.length-1; i++) { // -1 because last line is empty
			var grData=anaList[i].split('#:#'); // (gr,codeID,codeText)
			grIdx=(grData[0]*1)-1;
			analysisGroups[grIdx].push(grData[2]); // code Text
			analysisGroups[grIdx].push(grData[1]); // codeID (2x indexes)
		}
		displayGroup(myForm.currentGroup.value); // defines selectedGroup
	}
	checkSelectedAnalyses(myForm);
	myForm.submit();
}
// <--- AJAX
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
window.onload=function() {
	ajaxGetAnalysesList('$usedParent');
	ajaxUpdateRestrict($comparisonID);
	selectComparisonType('$compType');
	selectComparison($comparisonID,false);
}
</SCRIPT>
</HEAD>
<BODY topmargin=2 background="$promsPath{images}/bgProMS.gif"> <!--refreshSelAnalyses(document.compForm)-->
<CENTER>
<IFRAME name="storeFrame" style="display:none"></IFRAME>
<FORM name="compForm" method="post" target="listFrame" onsubmit="return checkSelectedAnalyses(this)">
<DIV style="white-space:nowrap">
<SELECT name="compType" id="compType" style="font-weight:bold;font-size:18px;color:#DD0000" onchange="selectComparisonType(this.value)"><OPTION value="1v1_pep"
|;
	print ' selected' if $compType eq '1v1_pep';
	print ">1 vs 1 Peptide</OPTION><OPTION value=\"1v1_prot\"";
	print ' selected' if $compType eq '1v1_prot';
	print ">1 vs 1 Protein</OPTION><OPTION value=\"full_prot\"";
	print ' selected' if $compType eq 'full_prot';
	print ">Full Protein</OPTION><OPTGROUP label=\"Modification sites:\">";
	foreach my $modID (sort{lc($projectVarMods{$a}[0]) cmp lc($projectVarMods{$b}[0])} keys %projectVarMods) {
		my $selStrg=($compType eq "modif_sites:$modID")? ' selected' : '';
		print "<OPTION value=\"modif_sites:$modID\"$selStrg>$projectVarMods{$modID}[0]-sites</OPTION>";
	}
	print qq
|</SELECT>&nbsp;<FONT class="title2">Comparison with&nbsp;</FONT><FONT class="title2"><SELECT name="numGroups" style="font-weight:bold;font-size:16px" onchange="updateNumGroups(this.value); setComparisonAsModified(2)"><OPTION value="0">No</OPTION>
|;
	#Print number of groups that can be compared (from 2 to $MAX_NUM_GROUPS)
	for (my $nbGroups=2; $nbGroups<=$MAX_NUM_GROUPS; $nbGroups++){
		print '<OPTION value="',$nbGroups,'">',$nbGroups,'</OPTION>';
	}
	print qq
|</SELECT> Groups</FONT>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<FONT class="title2">(Source:</FONT><SELECT name="comparisonID" style="font-weight:bold;font-size:14px;width:250px" onchange="selectComparison(this.value,true)"><OPTION value="0">New</OPTION>
|;
	foreach my $compID (sort{lc($comparisonList{$a}{'NAME'}) cmp lc($comparisonList{$b}{'NAME'})} keys %comparisonList) {
		print "<OPTION value=\"$compID\"";
		print ' selected' if $compID==$comparisonID;
		print ">$comparisonList{$compID}{NAME}</OPTION>\n";
	}
	print qq
|</SELECT>&nbsp;<INPUT type="button" class="title4" value="View All" onclick="displayAllComparisons()"/><FONT class="title2">)</FONT>
</DIV>
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="parentItem" value="$parentItem">
<INPUT type="hidden" name="ACT" value="list">
<INPUT type="hidden" name="anaList" value="">
<INPUT type="hidden" name="sort" value="$sortOrder">
<INPUT type="hidden" name="autocheck" value="">
<FONT style="font-size:2px"><BR></FONT>
<TABLE bgcolor=$darkColor2 cellpadding=0 border=0><TR>
<TR>
<TH valign=top rowspan=2><FONT class="title3">Select :</FONT><SELECT name="anaParent" class="title4" style="width:350px;font-weight:bold;font-size:12px" onchange="ajaxGetAnalysesList(this.value)">
<OPTION value="" disabled style="color:black;">----- Project items -----</OPTION>
|;
	####<Looping through list of analyses>####
	foreach my $refItem (@itemList) {
		print "\t<OPTION value=\"$refItem->[0]\"";
		print ' selected' if $refItem->[0] eq $usedParent;
		print ">$refItem->[1]</OPTION>\n";
	}
	print "<OPTION value=\"\" disabled style=\"color:black;\">----- Custom lists -----</OPTION>\n";
	if (scalar %classificationList) {
		foreach my $classID (sort{lc($classificationList{$a}[0]) cmp lc($classificationList{$b}[0])} keys %classificationList) {
			print "<OPTION value=\"CLASSIFICATION:$classID\">$classificationList{$classID}[0]</OPTION>\n";
			print ' selected' if "CLASSIFICATION:$classID" eq $usedParent;
		}
	}
	else {print "<OPTION value=\"\" disabled>None</OPTION>\n";}
	my ($selSpecAbs,$selSpecRel)=($siteMeasure eq 'specAbs')? (' selected','') : ($siteMeasure eq 'specRel')? ('',' selected') : ('','');
	my ($selUniquePep,$selUniqueSharedPep)=($pepSpecificity eq 'unique')? (' selected','') : ($pepSpecificity eq 'unique_shared')? ('',' selected') : ('','');
	print qq
|</SELECT><BR>
<SELECT multiple name="allAna" class="formSelect" style="width:450px;height:110px"></SELECT><BR>
<SELECT name="autoAction" class="formSelect"><OPTION value="select">Select</OPTION><OPTION value="unselect">Unselect</OPTION></SELECT> <FONT class="formText">items containing</FONT>
<INPUT type="text" name="selectText" style="font-size:11px;width:120px" value=""/>&nbsp;<INPUT type="button" class="formButton" value="Go" onclick="autoSelectAnalyses()"/>&nbsp;&nbsp;<INPUT type="button" class="formButton" value="Clear" onclick="clearSelection()"/></TH>
<TH valign=middle><BR><INPUT type="button" class="formButton" name="add" value=">" style="width:30px" onclick="updateAnalyses('add')"/><BR><INPUT type="button" class="formButton" name="remove" value="<" style="width:30px" onclick="updateAnalyses('remove')"/></TH>
<TH valign=top><FONT class="title3">Compared analyses - Group:</FONT><SELECT name="currentGroup" style="font-weight:bold;font-size:12px;" onchange="setCurrentGroup(this.value)"><OPTION value="">-</OPTION></SELECT><BR>
<SELECT multiple name="usedAna" class="formSelect" style="width:450px;height:70px"></SELECT></TH>
<TD valign=middle><BR><INPUT type="button" class="formButton" name="up" value="Up" style="width:50px" onclick="moveAnalyses('up')"/><BR><INPUT type="button" class="formButton" name="down" value="Down" style="width:50px" onclick="moveAnalyses('down')"/><BR><INPUT type="button" class="formButton" value="Clear" style="width:50px" onclick="clearAnalyses()"/></TD>
</TR>
<TR><TH colspan=3 align=left>
<INPUT type="checkbox" name="virtualData" value="1" /><FONT class="formText" onmouseover="popup('Peptides added by quantification algorithms<BR>using the \\'<B>M</B>atch <B>B</B>etween <B>R</B>uns\\' feature')" onmouseout="popout()">Include MBR-rescued peptides<SUP>*</SUP></FONT>
&nbsp;&nbsp;&nbsp;&nbsp;<FONT class="formText" onmouseover="popup('-<B>Proteotypic:</B> Peptides found in only 1 protein.<BR>-<B>Proteotypic + shared:</B> Use all peptides found for a protein if at least 1 is proteotypic.')" onmouseout="popout()">Peptide specificity<SUP>*</SUP>:</FONT><SELECT name="pepSpecificity"><OPTION value="all">All</OPTION><OPTION value="unique"$selUniquePep>Proteotypic</OPTION><OPTION value="unique_shared"$selUniqueSharedPep>Proteotypic + shared</OPTION></SELECT>

<DIV id="pepCompOptions">
<INPUT type="checkbox" name="noMissCut" value="1" /><FONT class="formText">Exclude peptides with missed cleavages</FONT>
<BR>&nbsp;<INPUT type="checkbox" name="fake" value="0" style="visibility:hidden"/><FONT class="formText">Peptide modification filters</FONT>&nbsp;<INPUT type="button" id="showModif" class="formButton" value="Open form" onclick="showModifFilterForm(true)" />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="submit" name="compare" value="Compare Analyses" style="width:180px;font-size:14px;font-weight:bold;">
</DIV>

<DIV id="protCompOptions"><INPUT type="checkbox" name="noHidden" value="1" checked onclick="setComparisonAsModified(2);document.getElementById('hiddenRule').disabled=(this.checked)?false:true"/><FONT class="formText" style="color:#DD0000">Ignore proteins </FONT><SELECT name="hiddenRule" id="hiddenRule" class="formSelect" style="color:#DD0000" onchange="setComparisonAsModified(2)">
<OPTION value="1">hidden everywhere</OPTION><OPTION value="2">where hidden</OPTION></SELECT>
&nbsp;&nbsp;<SPAN id="pepNumSPAN"><FONT class="formText">Peptides: &ge;</FONT><SELECT name="pepThreshold" class="formSelect" style="color:#DD0000" onchange="setComparisonAsModified(2)">
	<OPTION value="1">1</OPTION><OPTION value="2">2</OPTION><OPTION value="3">3</OPTION><OPTION value="4">4</OPTION><OPTION value="5">5</OPTION>
</SELECT><SELECT name="pepRule" class="formSelect" style="width:160px" onchange="setComparisonAsModified(2)">
	<OPTION value="ba">in Analysis</OPTION><OPTION value="banr">distinct in Analysis</OPTION>
</SELECT></SPAN><SPAN id="modifResSPAN"><FONT class="formText">Compare:</FONT><SELECT name="siteMeas" class="formSelect" onchange="setComparisonAsModified(2)"><OPTION value="occ">Occurence</OPTION><OPTION value="specAbs"$selSpecAbs>Sp. count</OPTION><OPTION value="specRel"$selSpecAbs>Sp. count %</OPTION></SELECT>&nbsp;&nbsp;<FONT class="formText">Res.:</FONT><SELECT name="siteRes" class="formSelect" onchange="setComparisonAsModified(2)"><OPTION value="">all</OPTION></SELECT></SPAN>
<BR><SELECT name="catFilterRule" class="formSelect" onchange="setComparisonAsModified(2)"><OPTION value="restrict">Restrict to</OPTION><OPTION value="exclude">Exclude</OPTION>
</SELECT>:<SPAN id="restrictSPAN" style="width:270px"><SELECT name="catFilter" style="width:270px"><OPTION value="">All proteins</OPTION>
</SELECT></SPAN>
&nbsp;<INPUT type="submit" name="compare" value="Compare Analyses" style="width:180px;font-size:14px;font-weight:bold;">
</DIV>

</TH></TR>
</TABLE>
<DIV id="modifDIV" class="popup">
<TABLE cellspacing=0 style="width:400px">
<TR><TD><INPUT type="button" class="formButton" value="Close form" onclick="showModifFilterForm(false)"/></TD></TR>
<TR bgcolor="$darkColor2"><TD>&nbsp;<B>Modification filters</B>&nbsp;&nbsp;<FONT class="formText">(multi-restrict logic:</FONT><SELECT name="restrictModifLogic" class="formSelect"><OPTION value="all">all</OPTION><OPTION value="any">any</OPTION></SELECT><FONT class="formText">)</FONT><B>:</B>&nbsp;</TD></TR>
|;
	my $bgColor=$lightColor1;
	foreach my $modID (sort{lc($modifications{$a}) cmp lc($modifications{$b})} keys %modifications) {
		print qq
|<TR bgcolor="$bgColor">
<TD nowrap><SELECT name="filterMod_$modID" class="formSelect"><OPTION value="allow">Allow</OPTION>
|;
		foreach my $filter ('exclude','restrict','ignore') {print "<OPTION value=\"$filter\">",ucfirst($filter),"</OPTION>";} # filter order must be matched in JS selectComparison function
		print qq
|</SELECT>&nbsp;<FONT class="formText">$modifications{$modID}</FONT>|;
		print '&nbsp;<INPUT type="checkbox" name="delocPhospho" value="1" /><FONT class="formText">Delocalize Phosphorylation sites</FONT>' if $modID==$phosphoID;
		print qq
|</TD>
</TR>
|;
		$bgColor = ($bgColor eq $lightColor1)? $darkColor2 : $lightColor1;
	}
	print "</TABLE>\n";
	print "<INPUT type=\"hidden\" name=\"modifList\" value=\"",join(',',sort{$a<=>$b} keys %modifications),"\">\n";
	print qq
|</DIV>
</FORM>
</CENTER>
<DIV id="divDescription" class="clDescriptionCont"></DIV>
<SCRIPT type="text/javascript">setPopup();</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

##########################################################################
#####################<1 vs 1 Peptide comparison>##########################
##########################################################################
sub comparePeptides {
	#my $delocPhospho=param('delocPhospho') || 0;

	&display1v1ComparisonHead;

	my $phosphoID=0;
	my $phosphoSpecif='STY';
	if ($delocPhospho) {
		($phosphoID)=$dbh->selectrow_array('SELECT ID_MODIFICATION FROM MODIFICATION WHERE UNIMOD_ACC=21');
	}

	####>Looping through analyses<####
	###<Queries
	my @sthAnaList;
	my $projQuery1=qq |SELECT E.NAME,G.NAME,SP.NAME,A.NAME FROM EXPERIMENT E,GEL2D G,SPOT SP,SAMPLE S,ANALYSIS A
							WHERE A.ID_ANALYSIS=? AND S.ID_SAMPLE=A.ID_SAMPLE AND SP.ID_SPOT=S.ID_SPOT AND G.ID_GEL2D=SP.ID_GEL2D AND E.ID_EXPERIMENT=G.ID_EXPERIMENT|;
	my $projQuery2=qq |SELECT E.NAME,S.NAME,A.NAME FROM EXPERIMENT E,SAMPLE S,ANALYSIS A
							WHERE A.ID_ANALYSIS=? AND S.ID_SAMPLE=A.ID_SAMPLE AND E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SPOT IS NULL|;
	push @sthAnaList,($dbh->prepare($projQuery1),$dbh->prepare($projQuery2)); # free sample then gel
	my $pepValidStrg=($virtualData)? '' : 'AND P.VALID_STATUS > 0';
	my $missCutStrg=($noMissCut)? 'AND MISS_CUT != 1' : '';
	my $sthPep=($pepSpecificity eq 'unique')?
		$dbh->prepare("SELECT PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&')
								FROM PEPTIDE P
								LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
								INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE AND IS_SPECIFIC=1
								WHERE P.ID_ANALYSIS=? $pepValidStrg $missCutStrg
								GROUP BY P.ID_PEPTIDE")
		: ($pepSpecificity eq 'unique_shared')?
		$dbh->prepare("SELECT PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&')
								FROM PEPTIDE P
								LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
								INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
								INNER JOIN ANALYSIS_PROTEIN AP ON P.ID_ANALYSIS=AP.ID_ANALYSIS AND AP.ID_PROTEIN=PPA.ID_PROTEIN AND PEP_SPECIFICITY=100
								WHERE P.ID_ANALYSIS=? $pepValidStrg $missCutStrg
								GROUP BY P.ID_PEPTIDE")
		: $dbh->prepare("SELECT PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&')
								FROM PEPTIDE P
								LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
								WHERE ID_ANALYSIS=? $pepValidStrg $missCutStrg
								GROUP BY P.ID_PEPTIDE");
	my $restrictModifs=scalar keys %{$modifFilters{restrict}};

	####<Scanning analyses>####
	my (@groupLabels,%peptideDistrib,%groupDistrib); #,%peptidesInGroup
	my $newLine=($exportType)? "\n" : '<BR>';
	foreach my $g (0..$#analysisGroups) {
		my $numAna=0;
		foreach my $codeID (@{$analysisGroups[$g]}) {
			unless ($exportType) {
				if ($numGroups) {
					$numAna++;
					my $numAnaInGr=scalar @{$analysisGroups[$g]};
					&updateProgress('Processing Group #'.($g+1)." / Analysis $numAna/$numAnaInGr..."); # $numAna can be > $numAnaInGr if >1 lists in gr!
				}
				else {
					&updateProgress('Processing Analysis '.($g+1)."/$numAnaGroups..."); # $numAna can be > $numAnaInGr if >1 lists in gr!
				}
			}

			my ($type,$anaID)=split(':',$codeID); # Can only be 'A'
			my @typeInfo;
			foreach my $sthAna (@sthAnaList) {
				$sthAna->execute($anaID);
				@typeInfo=$sthAna->fetchrow_array;
				last if $typeInfo[0]; # 2 loops in case of project or experiment
			}
			if ($numGroups) {
				unless ($groupLabels[$g][0]) {$groupLabels[$g][0]='Group #'.($g+1);} # groupPos is irrelevent here
				$groupLabels[$g][1].=$newLine if $groupLabels[$g][1];
				$groupLabels[$g][1].='-'.join(' > ',@typeInfo);
			}
			else {
				@{$groupLabels[$g]}=(&promsMod::shortenName($typeInfo[-1],20),join(' > ',@typeInfo));
			}

			###<peptides
			$sthPep->execute($anaID);
			while (my ($pepSeq,$varModStrg)=$sthPep->fetchrow_array) {
				my $ambigCode='';
				my $newVarModStrg='';
				my $okPeptide=1; # default
				if ($varModStrg) {
					#<Restrict filter
					my $okRestrict=0;
					if ($restrictModifs) {
						foreach my $modID (keys %{$modifFilters{restrict}}) {
							$okRestrict++ if $varModStrg=~/(^|&)$modID:/;
						}
						$okRestrict=0 if ($restrictModifLogic eq 'all' && $okRestrict < $restrictModifs); # all: all restrict modifs must be in peptide, any: peptide OK if any is found
					}
					else {$okRestrict=1;}
					next unless $okRestrict; # next peptide

					if (($delocPhospho && $varModStrg=~/(^|&)$phosphoID:/) || $modifFilters{active}) {
						foreach my $modifStrg (split('&',$varModStrg)) {
							my ($modID,$posStrg)=split(':',$modifStrg);

							#<Exclude filter
							if ($modifFilters{exclude}{$modID}) {
								$okPeptide=0;
								last;
							}

							#<Ignore filter
							next if $modifFilters{ignore}{$modID}; # next modif

							#<Delocalize Phospho
							if ($delocPhospho && $modID==$phosphoID) {
								my $numModifRes = () = $posStrg =~ /(\d+)/g;
								my @targetResIdx;
								while ($pepSeq =~ /[$phosphoSpecif]/g) {push @targetResIdx,$-[0];} # position of match (1st=0)
								my $numTargetRes=scalar @targetResIdx;
								if ($numModifRes==$numTargetRes) {
									$ambigCode='-P'.$posStrg;
								}
								else { # ambiguity
									$ambigCode='-P'.($targetResIdx[0]+1).'~'.($targetResIdx[-1]+1).':'.$numModifRes.'/'.$numTargetRes;
								}
							}
							else {
								$newVarModStrg.='&'.$modifStrg;
							}
						}
					}
					else {
						$newVarModStrg=$varModStrg;
					}
				}
				elsif ($restrictModifs) { # no modif
					$okPeptide=0;
				}
				next unless $okPeptide;
				$peptideDistrib{$g}{$pepSeq.$newVarModStrg.$ambigCode}=1;
				$groupDistrib{$pepSeq.$newVarModStrg.$ambigCode}{$g}=1;
			}
		}
	}
	foreach my $sth (@sthAnaList) {$sth->finish;}
	$sthPep->finish;
	$dbh->disconnect;

	&display1v1ComparisonBody('peptide',\@groupLabels,\%peptideDistrib,\%groupDistrib);

	exit;
}




##########################################################################
###################<Modification sites comparison>########################
##########################################################################
sub compareModificationSites {
	my ($selModifID)=($compType =~ /^modif_sites:(\d+)/);
	$CONTEXT_SIZE++ unless $CONTEXT_SIZE % 2; # must be odd number
	my $MID_CONTEXT=int($CONTEXT_SIZE/2);

	####<Fetching Comparison autocheck parameters (if any)>####
	my ($filterLogic,%autocheckParams,$okAutocheck);
	if ($autochkStrg) {
		#my $autochkStrg='1:>=:3;2:>=:4';
		($filterLogic,my $pepFilterStrg)=split('#',$autochkStrg);
		if ($filterLogic eq 'frequency') {
			($autocheckParams{'FREQ_ANAGR'})=split(';',$pepFilterStrg); # just in case usi
		}
		else { # and,or
			foreach my $grData (split(';',$pepFilterStrg)) {
				my ($aGr,$rel,$numPep)=split(':',$grData);
				@{$autocheckParams{'PEP'}{$aGr}}=($rel,$numPep);
				$okAutocheck=1 if ($rel ne '>=' || $numPep>0); # do not check all proteins if ">=0" for all anaGroups
			}
		}
	}
	else {$filterLogic='and';}
	my $numProtChecked=scalar keys %selectedProteins;
	$okAutocheck=0 if $numProtChecked; # check by Perl


	my ($jobDir,$timeStamp2,$exportTitleStrg); # for export only
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
|;
	if ($exportType) {
		mkdir "$promsPath{tmp}/scratch" unless -d "$promsPath{tmp}/scratch";
		&promsMod::cleanDirectory("$promsPath{tmp}/scratch/export",'15m');
		$jobDir=strftime("%Y%m%d%H%M%S",localtime).'_'.$userID;
		$timeStamp2=strftime("%Y-%m-%d %H:%M",localtime);
		mkdir "$promsPath{tmp}/scratch/export/$jobDir";
		$exportTitleStrg="<FONT class=\"title\">Exporting $projectVarMods{$selModifID}[0]-Sites Comparison</FONT>";
	}
	else {
		$exportTitleStrg='';
		print qq
|<STYLE type="text/css">
.small {font-size:10px;}
.TD {font-weight:normal;}
.TH {font-weight:bold;}
ACRONYM {cursor:help;}
.row_0 {background-color:$color2;}
.row_1 {background-color:$color1;}
.modifSite {font-family:monospace;font-weight:bold;font-size:18px;}
.$projectVarMods{$selModifID}[2] {$projectVarMods{$selModifID}[3]}
</STYLE>
|;
		if ($numAnaGroups<=5) {
			print qq
|<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/vennDiagram.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
|;
		}
		print qq
|<SCRIPT type="text/javascript">
var selectedSite;
var isNav = (navigator.appName.indexOf("Netscape") !=-1);
function unselectSite() {
	if (selectedSite) {
		document.getElementById('siteModifDetailsDIV').style.display='none';
		document.getElementById('IMG_'+selectedSite).src="$promsPath{images}/plus.gif";
	}
	selectedSite=null;
}
function ajaxGetModifSiteDetails(e,site,pepIdStrg) {
	var siteDIV=document.getElementById('siteModifDetailsDIV');
	if (!siteDIV) { // page is not fully loaded yet
		alert('Display process is still on-going. Try again when completed.');
		return;
	}
	if (site==selectedSite) { // image is minus1.gif
		unselectSite();
		return;
	}
	unselectSite();
	selectedSite=site;
	document.getElementById('IMG_'+site).src="$promsPath{images}/minus1.gif";
	var divX = (isNav)? e.pageX : event.clientX + document.body.scrollLeft; divX-=5;
	var divY = (isNav)? e.pageY : event.clientY + document.body.scrollTop; divY+=10;
	siteDIV.style.left = divX + 'px';
	siteDIV.style.top = divY + 'px';
	siteDIV.innerHTML="<CENTER><BR><FONT class=\\"title3\\">Fetching data...</FONT><BR><IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><INPUT type=\\"button\\" value=\\" Cancel \\" class=\\"font11\\" onclick=\\"document.getElementById('siteModifDetailsDIV').style.display='none'\\"><BR><BR>";
	siteDIV.style.display='';
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	/* POST */
	var paramStrg="AJAX=getSiteDetails&modID=$selModifID&site="+site+"&numGroups=$numGroups&anaList=$anaList&pepIDs="+pepIdStrg;
	XHR.open("POST","$promsPath{cgi}/compareAnalyses.cgi",true);
	//Send the proper header information along with the request
	XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
	XHR.setRequestHeader("Content-length", paramStrg.length);
	XHR.setRequestHeader("Connection", "close");
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			siteDIV.innerHTML=XHR.responseText;
		}
	}
	XHR.send(paramStrg);
}
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
function generateAutoCheckString() {
	var protForm=document.protView;
	//logic;
	var autoChkStrg=protForm.logic.value+'#';
	if (protForm.logic.value=='frequency') {
		autoChkStrg+=protForm.freqAnaGr.value+';1';
	}
	else {
		//peptides
		for (var g=0; g<$numAnaGroups; g++) {
			if (g>0) {autoChkStrg+=';'}
			//var aGr=g+1;
			//autoChkStrg+=aGr+':'+protForm['pepCompSign_'+g].value+':'+protForm['numPep_'+g].value;
			autoChkStrg+=(g+1)+':';
			autoChkStrg+=(protForm['numPep_'+g].value*1 >= 0)? '>=:'+protForm['numPep_'+g].value : '<=:0'; // 0/1='Any,Found' : -1='Not found'
		}
	}
	return autoChkStrg;
}
function selectSort(newSort) {
	parent.selItemsFrame.setComparisonAsModified(2); // 2! because stored in top.promsFrame.selectedView
	top.promsFrame.selectedView=newSort; //openProject.cgi
	parent.selItemsFrame.document.compForm.sort.value=newSort;
	var myForm=document.protView;
	myForm.exportType.value=null;
	myForm.sort.value=newSort;
	myForm.autocheck.value=generateAutoCheckString();
	myForm.action="./compareAnalyses.cgi";
	myForm.target='listFrame';
	myForm.submit();
}
function exportSites(newExportType) {
	if (newExportType=='sel' && numProtChecked==0) {
		alert('ERROR: No site selected!');
		return;
	}
	window.open('','ExportWindow','width=600,height=200,location=no,resizable=yes,scrollbars=yes');
	var myForm=document.protView;
	myForm.exportType.value=newExportType;
	myForm.action="./compareAnalyses.cgi";
	myForm.target='ExportWindow';
	myForm.submit();
}
function sequenceView(proteinId,analysisId){
	var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+analysisId+"&id_prot="+proteinId+"&msdata="+top.showMSdata;
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
function uncheckAll(checkBoxList) {
	parent.selItemsFrame.setComparisonAsModified(1);
	if (checkBoxList.length) {
		for (var i=0; i < checkBoxList.length; i++) {checkBoxList[i].checked=false;}
	}
	else {checkBoxList.checked=false;}
	updateNumChecked(-numProtChecked);
	hideShowUnchecked(0);
}
function manualCheck(chkStatus) {
	parent.selItemsFrame.setComparisonAsModified(1);
	var newChecked=(chkStatus)? 1 : -1;
	updateNumChecked(newChecked);
}
function updateComparisonStrategy(strategy) {
	var logicPepDisabStatus=(strategy=='frequency')? true : false;
	for (let a=0; a < $numAnaGroups; a++) {
		document.getElementById('numPep_'+a).disabled=logicPepDisabStatus;
	}
	document.getElementById('freqAnaGr').disabled = !logicPepDisabStatus;
}
function autoCheck(warn) {
	var protForm=document.protView;
	var checkBoxList=protForm.chkProt;
	//Peptides
	var pepDistribList=protForm.pepDistrib;
	var numPep=[];
	for (let a=0; a < $numAnaGroups; a++) {
		numPep[a]=document.getElementById('numPep_'+a).value*1; //convert to number
	}
	var logic=document.getElementById('logic').value; //inter-group logic
	var freqAnaGr=protForm.freqAnaGr;

	//Scanning protein list
	var numMatch=0;
	var allMatches=0;
	PROT: for (let i=0; i<pepDistribList.length; i++) {
		var anaDistrib=pepDistribList[i].value.split(':');
		if (logic=='and') {
			//Peptides
			for (let a=0; a<anaDistrib.length; a++) {
				if ((anaDistrib[a]*1==1 && numPep[a]==-1) \|\| (anaDistrib[a]*1==0 && numPep[a]==1)) {
					continue PROT;
				}
			}
		}
		else if (logic=='or') {
			var okMatch;
			for (let a=0; a<anaDistrib.length; a++) {
				okMatch=true; // reset for each group
				if ((anaDistrib[a]*1==1 && numPep[a]==-1) \|\| (anaDistrib[a]*1==0 && numPep[a]==1)) {
					okMatch=false;
				}
				if (okMatch==true) {break;}
			}
			if (okMatch==false) {continue PROT;}
		}
		else { //logic = frequency
			var numMatchGr=0;
			for (let a=0; a<anaDistrib.length; a++) {
				if (anaDistrib[a]*1==1) {numMatchGr++;}
			}
			if (numMatchGr < freqAnaGr.value) {continue PROT;}
		}
		//Result
		if (!checkBoxList[i].checked) { //okMatch &&
			numMatch++;
			checkBoxList[i].checked=true;
			var chkData=checkBoxList[i].value.split(':');
			var tr=document.getElementById('prot_'+chkData[chkData.length-1]);
			tr.className='list row_'+(numMatch % 2);
			tr.style.display=''; // make sure checked prot is visible
		}
		allMatches++;
	}
	if (numMatch) {updateNumChecked(numMatch);}
	if (warn==1) {
		if (numMatch) {parent.selItemsFrame.setComparisonAsModified(1);}
		alert(allMatches+' sites matched selected criteria ('+numMatch+' new)!');
	}
}
function hideShowUnchecked(warn) {
	var myForm=document.protView;
	if (numProtChecked==0 && myForm.hideShowButton.value=='Hide Unchecked') {
		if (warn) alert('There are no sites checked!');
		return;
	}
	if (myForm.hideShowStatus.value=='none') {
		myForm.hideShowStatus.value=''; // show all
		myForm.hideShowButton.value='Hide Unchecked';
	}
	else { // '' ~ block
		myForm.hideShowStatus.value='none'; // show checked only
		myForm.hideShowButton.value='Show Unchecked';
	}
	var checkBoxList=myForm.chkProt;
	if (checkBoxList.length) {
		var numVisible=0;
		for (var i=0; i < checkBoxList.length; i++) {
			var chkData=checkBoxList[i].value.split(':');
			var tr=document.getElementById('prot_'+chkData[chkData.length-1]);
			if (myForm.hideShowStatus.value=='none') { // show checked only
				if (checkBoxList[i].checked) {
					numVisible++;
					tr.className='list row_'+(numVisible % 2);
					tr.style.display='';
				}
				else {
					tr.style.display='none';
				}
			}
			else { // show all
				numVisible++;
				tr.className='list row_'+(numVisible % 2);
				tr.style.display='';
			}
		}
	}
	else if (!checkBoxList.checked) { // only 1 prot in comparison
		var chkData=checkBoxList.value.split(':');
		document.getElementById('prot_'+chkData[chkData.length-1]).style.display=myForm.hideShowStatus.value;
	}
}
function updateNumChecked(newChecked) {
	numProtChecked+=newChecked;
	document.getElementById('numProtChecked').innerHTML=numProtChecked;
}
var numProtChecked=$numProtChecked;
|;
		##>Save proteins to custom Lists
		&promsMod::printAjaxManageSaveProteins($projectID,\%promsPath,'document.protView.chkProt','updateLists');

		&promsMod::popupInfo();
	} # end of no export
	print qq
|</SCRIPT>
</HEAD>
<BODY style="margin:2px;" background="$promsPath{images}/bgProMS.gif">
<CENTER>
$exportTitleStrg
<DIV id="waitDIV">
<BR><BR><BR><FONT class="title3">Fetching data...</FONT><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">
<BR><BR><BR><FONT class="title3">Status: <SPAN id="waitSPAN"><SPAN>...</FONT>
</DIV>
</CENTER>
|;

	####<Category filter>####
	my ($filterClassName,$filterCatName);
	if ($catFilter) {
		my $sthCP=$dbh->prepare("SELECT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$catFilter");
		%{$catProteins{$catFilter}}=();
		$sthCP->execute;
		while (my ($protID)=$sthCP->fetchrow_array) {
			$catProteins{$catFilter}{$protID}=1;
		}
		$sthCP->finish;
		($filterClassName,$filterCatName)=$dbh->selectrow_array("SELECT CLASSIFICATION.NAME,CATEGORY.NAME FROM CLASSIFICATION,CATEGORY WHERE CLASSIFICATION.ID_CLASSIFICATION=CATEGORY.ID_CLASSIFICATION AND ID_CATEGORY=$catFilter");
	}

	####<Looping through analyses>####
	###<Queries
	my @sthAnaList;
	my $projQuery1=qq |SELECT FILE_FORMAT,E.NAME,G.NAME,SP.NAME,A.NAME FROM EXPERIMENT E,GEL2D G,SPOT SP,SAMPLE S,ANALYSIS A
							WHERE A.ID_ANALYSIS=? AND S.ID_SAMPLE=A.ID_SAMPLE AND SP.ID_SPOT=S.ID_SPOT AND G.ID_GEL2D=SP.ID_GEL2D AND E.ID_EXPERIMENT=G.ID_EXPERIMENT|;
	my $projQuery2=qq |SELECT FILE_FORMAT,E.NAME,S.NAME,A.NAME FROM EXPERIMENT E,SAMPLE S,ANALYSIS A
							WHERE A.ID_ANALYSIS=? AND S.ID_SAMPLE=A.ID_SAMPLE AND E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SPOT IS NULL|;
	push @sthAnaList,($dbh->prepare($projQuery1),$dbh->prepare($projQuery2)); # free sample then gel

	my $pepValidStrg=($virtualData)? '' : 'AND P.VALID_STATUS > 0';
	my $protVisStrg=($hiddenRule==2)? 'AND VISIBILITY>=1' : '';
	my $protPepSpecifStrg=($pepSpecificity=~/unique/)? 'AND PEP_SPECIFICITY=100' : '';
	my $isSpecificStrg=($pepSpecificity eq 'unique')? 'AND IS_SPECIFIC=1' : '';
	my $sthSites=$dbh->prepare("SELECT PPA.ID_PROTEIN,AP.VISIBILITY,AP.CONF_LEVEL,P.ID_PEPTIDE,P.PEP_SEQ,ABS(PEP_BEG),ABS(PEP_END),GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',COALESCE(PM.REF_POS_STRING,'') ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),P.DATA,P.SPEC_COUNT,P.SCORE
								FROM PEPTIDE_MODIFICATION PM
								INNER JOIN PEPTIDE P ON P.ID_PEPTIDE=PM.ID_PEPTIDE
								INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE $isSpecificStrg
								INNER JOIN ANALYSIS_PROTEIN AP ON AP.ID_PROTEIN=PPA.ID_PROTEIN AND AP.ID_ANALYSIS=PPA.ID_ANALYSIS $protPepSpecifStrg
								WHERE P.ID_ANALYSIS=? $pepValidStrg $protVisStrg
								GROUP BY ABS(PEP_BEG),P.ID_PEPTIDE,PPA.ID_PROTEIN");

	####<Scanning analyses>####
	my (@groupLabels,%proteinSites,%groupDistrib,%siteDistrib,%peptideIDs,%bestVisibility,%distinctAna,@numAnaInGr,%anaType,%proteinSeq,%sitesInProtein,%siteContext,%unlocalizedPep,%excludedProteins);
	my $newLine=($exportType)? "\n" : '<BR>';
	foreach my $g (0..$#analysisGroups) {
		$numAnaInGr[$g]=scalar @{$analysisGroups[$g]};
		my (%localDistrib,%bestSpecCount,%peptideSites);
		my $numAna=0;
		foreach my $codeID (@{$analysisGroups[$g]}) {
			if ($numGroups) {
				$numAna++;
				&updateProgress('Processing Group #'.($g+1)." / Analysis $numAna/$numAnaInGr[$g]..."); # $numAna can be > $numAnaInGr if >1 lists in gr!
			}
			else {
				&updateProgress('Processing Analysis '.($g+1)."/$numAnaGroups..."); # $numAna can be > $numAnaInGr if >1 lists in gr!
			}
			my ($type,$anaID)=split(':',$codeID); # Can only be 'A'
			my @typeInfo;
			foreach my $sthAna (@sthAnaList) {
				$sthAna->execute($anaID);
				($anaType{$anaID},@typeInfo)=$sthAna->fetchrow_array;
				last if $typeInfo[0]; # 2 loops in case of project or experiment
			}
			if ($numGroups) {
				unless ($groupLabels[$g][0]) {$groupLabels[$g][0]='Group #'.($g+1);} # groupPos is irrelevent here
				$groupLabels[$g][1].=$newLine if $groupLabels[$g][1];
				$groupLabels[$g][1].='-'.join(' > ',@typeInfo);
			}
			else {
				@{$groupLabels[$g]}=(&promsMod::shortenName($typeInfo[-1],20),join(' > ',@typeInfo));
			}

			###<Modif sites
			$sthSites->execute($anaID);
			my $siteCount=0;
			#while (my ($protID,$protVis,$pepID,$pepSeq,$pepBeg,$posStrg,$refPosStrg,$pepData,$specCount)=$sthSites->fetchrow_array) {#}
			while (my ($protID,$protVis,$protConf,$pepID,$pepSeq,$pepBeg,$pepEnd,$varModCode,$refModInfo,$pepData,$specCount,$score)=$sthSites->fetchrow_array) { #,$charge
				unless ($pepBeg) { # cannot locate site on protein
					$unlocalizedPep{$pepSeq}=1;
					$excludedProteins{$protID}=1; # temp: will be deleted later if prot is matched by localized pepSeq
					next;
				}
				next if ($protVis==0 && $hiddenRule); # hiddenRule=1 (cannot be 2 if vis=0)
				next if ($catFilter && ($catFilterRule eq 'restrict' && !$catProteins{$catFilter}{$protID}) || ($catFilterRule eq 'exclude' && $catProteins{$catFilter}{$protID}));
				next if (!$varModCode || $varModCode !~ /(^|&)$selModifID:([^:&]+)/); # must contains vmod of interest
				my $posStrg=$2;
				#<Generate coverted seq
				if ($exportType) {
					my $protLength=length($proteinSeq{$protID});
					if ($proteinSeq{$protID}) {
						$proteinSeq{$protID}.='X' x ($pepEnd-$protLength) if $protLength < $pepEnd;
					}
					else {$proteinSeq{$protID}='X' x $pepEnd;}
					substr($proteinSeq{$protID},$pepBeg-1,$pepEnd-$pepBeg+1,$pepSeq); # works on index
				}
				#<Site probability
				my (@prsData,%probMQ);
				if ($projectVarMods{$selModifID}[0]=~/Phospho/ && $pepData && $pepData =~ /PRS=([^##]+)/) {
					@prsData=($1);
					my ($status,$proba,$positions)=split(';',$prsData[0]);
					push @prsData,$proba,$pepBeg-1,$pepSeq;
				}
				if ($refModInfo && $refModInfo=~/(^|&)$selModifID:[^&#]*##PRB_MQ=([^&#]+)/) {
					my $modProbStrg=$2;
					while ($modProbStrg =~ /([^,:]+):([^,]+)/g) {
						$probMQ{$1}=$2;
					}
				}
				#my ($x,$refPosStrg)=($refModInfo=~/(^|&)$selModifID:.*#PROB_MQ=([^:&]+)/);
				my @sequence=split(//,$pepSeq);
				my $siteOK=0;
				foreach my $pos (split(/\./,$posStrg)) {
					my $siteRes='';
					my $protPos=0;
					if ($pos=~/;/) { # ambiguity
						next;
					}
					elsif ($pos=~/\D/) { # possible = - * +
						next if $pos=~/[=\*]/; # irrelevant for "Any-term" modifs
						next if ($targetedRes && $pos ne $targetedRes);
						#$siteRes=($pos eq '-')? 'Protein N-term' : ($pos eq '=')? 'N-term' : ($pos eq '+')? 'Protein C-term' : 'C-term'; # '*'
						#$siteRes=($pos eq '-')? 'Protein N-term' : 'Protein C-term';
						$siteRes=$pos;
						$protPos=$pos;
					}
					else {
						my $posIdx=$pos-1;
						$siteRes=$sequence[$posIdx];
						next if ($targetedRes && $siteRes ne $targetedRes);
						$protPos=$pos+$pepBeg-1;
					}
					next unless $siteRes;
					my $site="$protID-$siteRes$protPos";
					next if ($exportType && $exportType eq 'sel' && !$selectedProteins{$protID} && !$selectedProteins{$site});
					$siteOK=1;
					@{$proteinSites{$site}}=($protID,$siteRes,$protPos,{},{}) unless $proteinSites{$site}; # only for 1st time seen
					if (!$groupDistrib{$site} || !$groupDistrib{$site}{$g} || $groupDistrib{$site}{$g}[0] < $protVis) { # best vis in gr
						$groupDistrib{$site}{$g}[0]=$protVis;
						$groupDistrib{$site}{$g}[1]=$protConf;
					}
					if ($exportType) {
						$sitesInProtein{$protID}{$site}=$protPos;
					}
					$localDistrib{$site}{$anaID}++; # presence/absence in each ana in group (1 per ana) + occurence per ana in group
					$groupDistrib{$site}{$g}[3]++; # total occurence in gr (multiple per ana)
#my $pepIon="$pepSeq:$varModCode:$charge";
					my $seqVmod="$pepSeq:$varModCode"; # spectral count is charge-independant in myProMS. TODO: Check if this is the proper way to to do spectral count
					$bestSpecCount{$seqVmod}{$anaID}=$specCount if ($specCount && (!$bestSpecCount{$seqVmod} || !$bestSpecCount{$seqVmod}{$anaID} || $bestSpecCount{$seqVmod}{$anaID} < $specCount)); # keys best spec count/ana
					$peptideSites{$seqVmod}{$site}=1;
					push @{$peptideIDs{$site}{$g}},$pepID;
					#>Best PRS Proba
					if ($prsData[0]) {
						if (!$proteinSites{$site}[3]{$g} || $prsData[1] > $proteinSites{$site}[3]{$g}[1]) { # keep best PRS
							$proteinSites{$site}[3]{$g}=\@prsData;
						}
					}
					#>Best MQ Proba
					if ($probMQ{$pos}) {
						if (!$proteinSites{$site}[4]{$g} || $probMQ{$pos} > $proteinSites{$site}[4]{$g}) { # keep best PRS
							$proteinSites{$site}[4]{$g}=$probMQ{$pos};
						}
					}
					elsif ($anaType{$anaID} eq 'MAXQUANT.DIR' && $score) {$proteinSites{$site}[4]{$g}=1;} # proba 100%
				}
				next unless $siteOK;
				$bestVisibility{$protID}=$protVis if (!$bestVisibility{$protID} || $bestVisibility{$protID}<$protVis);
				$distinctAna{$protID}{$anaID}=1; # in case same ana in different groups

#if ($exportType) { # generate coverted seq
#	my $protLength=length($proteinSeq{$protID});
#	if ($proteinSeq{$protID}) {
#		$proteinSeq{$protID}.='X' x ($pepEnd-$protLength) if $protLength < $pepEnd;
#	}
#	else {$proteinSeq{$protID}='X' x $pepEnd;}
#	substr($proteinSeq{$protID},$pepBeg-1,$pepEnd-$pepBeg+1,$pepSeq); # works on index
#}

				$siteCount++;
				if ($siteCount > 5000) {
					&updateProgress('.','+=');
					$siteCount=0;
				}
			}
		}
		#<Update spectral count
		foreach my $seqVmod (keys %bestSpecCount) {
			my $totSpecCount=0;
			foreach my $specCount (values %{$bestSpecCount{$seqVmod}}) {$totSpecCount+=$specCount;} # sum accros all ana in group
			$totSpecCount=1 unless $totSpecCount;
			foreach my $site (keys %{$peptideSites{$seqVmod}}) {$groupDistrib{$site}{$g}[4]+=$totSpecCount;}
		}
		#<Update group occurence
		my $idx=($siteMeasure =~ /^spec/)? 4 : 2;
		foreach my $site (keys %localDistrib) {
			$groupDistrib{$site}{$g}[2]=1*(sprintf "%.1f",100 * (scalar keys %{$localDistrib{$site}})/$numAnaInGr[$g]);
			$siteDistrib{$site}[0]+=$groupDistrib{$site}{$g}[$idx]; # total occurence in dataset
		}
	}
	foreach my $sth (@sthAnaList) {$sth->finish;}
	$sthSites->finish;

	###<Setting undef values to 0 (required for sortModifSite) & computing Relative spec Count if siteMeasure=specRel
	foreach my $site (keys %groupDistrib) {
		my $maxSpec=0;
		foreach my $g (0..$#analysisGroups) {
			@{$groupDistrib{$site}{$g}}=(0,0,0,0,0) unless $groupDistrib{$site}{$g};
			$maxSpec=$groupDistrib{$site}{$g}[4] if $maxSpec < $groupDistrib{$site}{$g}[4];
		}
		if ($siteMeasure eq 'specRel') {
			foreach my $g (0..$#analysisGroups) {
				$groupDistrib{$site}{$g}[4]=1*(sprintf "%.1f",100 * $groupDistrib{$site}{$g}[4]/$maxSpec);
			}
		}
	}

	if ($sortOrder eq 'peptide') {
		foreach my $site (keys %siteDistrib) {
			my $aveOcc=$siteDistrib{$site}[0]/$numAnaGroups;
			$siteDistrib{$site}[1]=0;
			foreach my $g (0..$#analysisGroups) {
				$siteDistrib{$site}[2]+=abs($groupDistrib{$site}{$g}[3]-$aveOcc); # sum of deltas to average
			}
		}
	}
	####<Fetching protein info>####
	&updateProgress('Fetching protein annotation...'); # unless $exportType;

	my (%protAnnotation,%masterProteins);
	my %sthProt=(
				 SEQ=>$dbh->prepare("SELECT PROT_SEQ,ALIAS,PROT_DES,MW,ORGANISM,ID_MASTER_PROTEIN FROM PROTEIN WHERE ID_PROTEIN=?"),
				 NO_SEQ=>$dbh->prepare("SELECT NULL,ALIAS,PROT_DES,MW,ORGANISM,ID_MASTER_PROTEIN FROM PROTEIN WHERE ID_PROTEIN=?")
	);
	my ($geneNameID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
	my $sthMPG=$dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN=? AND ID_IDENTIFIER=$geneNameID ORDER BY RANK");
	my $sthMPS=$dbh->prepare("SELECT PROT_SEQ FROM MASTER_PROTEIN WHERE ID_MASTER_PROTEIN=?");

	my $count=0;
	foreach my $protID (keys %bestVisibility) {
		$count++;
		if ($count >= 1000) {
			&updateProgress('.','+='); # unless $exportType;
			$count=0;
		}
		delete $excludedProteins{$protID}; # remove protein from this list if it is part of it
		my @badCoverage;
		my $seqType='NO_SEQ'; # default
		if ($exportType) { # check if protSeq is needed
			my $currProtLength=length($proteinSeq{$protID});
			foreach my $site (keys %{$sitesInProtein{$protID}}) {
				my $sitePos=$sitesInProtein{$protID}{$site};
				my ($extraNt,$extraCt)=('','');
				my ($contBegIdx,$contLength)=(0,$CONTEXT_SIZE);
				if ($sitePos eq '-') { # N-term modif
					$contBegIdx=0;
				}
				elsif ($sitePos eq '+') { # C-term modif
					$contBegIdx=$currProtLength-$CONTEXT_SIZE-1;
				}
				else {
					if ($sitePos < $MID_CONTEXT+1) { # before beg of prot seq
						$extraNt='_' x ($MID_CONTEXT+1-$sitePos);
						$contBegIdx=0;
						$contLength-=($MID_CONTEXT+1-$sitePos);
					}
					else {$contBegIdx=$sitePos-($MID_CONTEXT+1);}
					if ($sitePos+$MID_CONTEXT > $currProtLength) { # after end of prot seq
						$extraCt='_' x ($sitePos+$MID_CONTEXT-$currProtLength);
						$contLength-=($sitePos+$MID_CONTEXT-$currProtLength);
					}
				}
				my $context=substr($proteinSeq{$protID},$contBegIdx,$contLength);
				if ($context=~/X/ || $extraCt) { # incomplete coverage
					push @badCoverage,$site;
					$seqType='SEQ';
				}
				@{$siteContext{$site}}=($extraNt,$context,$extraCt,$contBegIdx,$contLength);
			}
		}
		$sthProt{$seqType}->execute($protID);
		(my $protSeq,@{$protAnnotation{$protID}})=$sthProt{$seqType}->fetchrow_array;
		$protAnnotation{$protID}[2]=($protAnnotation{$protID}[2])? sprintf "%.1f",$protAnnotation{$protID}[2]/1000 : 0; # MW
		if ($protAnnotation{$protID}[4]) {
			my $masterProtID=$protAnnotation{$protID}[4];
			unless ($masterProteins{$masterProtID}) {
				@{$masterProteins{$masterProtID}}=();
				$sthMPG->execute($masterProtID);
				while (my ($gene)=$sthMPG->fetchrow_array) {
					push @{$masterProteins{$masterProtID}},$gene;
				}
			}
			if ($seqType eq 'SEQ') {
				next if (!$protSeq || $protSeq eq '-');
				if ($protSeq eq '+') {
					$sthMPS->execute($masterProtID);
					($protSeq)=$sthMPS->fetchrow_array;
				}
			}
		}
		#<Reextract surounding sequence from real protSeq
		foreach my $site (@badCoverage) { # only if exportType
			if ($siteContext{$site}[2]) { # extraCt => reextract from real protSeq
				my $sitePos=$sitesInProtein{$protID}{$site};
				my $protLength=length($protSeq);
				next if $sitePos > $protLength; # in case wrong prot seq recorded!!!
				$siteContext{$site}[2]=''; # reset extraCt
				$siteContext{$site}[4]=$CONTEXT_SIZE; # reset context length
				if ($sitePos+6 > $protLength) { # after end of prot seq
					$siteContext{$site}[2]='_' x ($sitePos+$MID_CONTEXT-$protLength);
					$siteContext{$site}[4]-=(length($siteContext{$site}[0])+length($siteContext{$site}[2]));
				}
			}
			$siteContext{$site}[1]=substr($protSeq,$siteContext{$site}[3],$siteContext{$site}[4]); #$siteContext{$site}[4]);
		}
	}
	$sthProt{SEQ}->finish;
	$sthProt{NO_SEQ}->finish;
	$sthMPG->finish;
	$sthMPS->finish;

	$dbh->disconnect;


	#####################
	####<Export mode>####
	#####################
	if ($exportType) {
		&updateProgress('Writing data to worksheet...');
		my $xlsFileName="compare_sites_$jobDir.xls";
		my $workbook=Spreadsheet::WriteExcel->new("$promsPath{tmp}/scratch/export/$jobDir/$xlsFileName");
		eval { # Perl 5.8 compatibility
			$workbook->set_properties(title=>"Compare Project Items",
									  author=>'myProMS server',
									  comments=>'Automatically generated with Perl and Spreadsheet::WriteExcel'
									  );
		};
		$workbook->set_custom_color(40,189,215,255); # light light color #BDD7FF (new myProMS colors V3.5+)
		$workbook->set_custom_color(41,128,179,255); # dark blue color  #80B3FF
		my %itemFormat=(
				title =>			$workbook->add_format(align=>'center',size=>18,bold=>1,border=>1),
				header =>			$workbook->add_format(align=>'center',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
				mergeRowHeader =>	$workbook->add_format(align=>'center',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
				mergeColHeader =>	$workbook->add_format(align=>'center',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
				identVis =>			$workbook->add_format(align=>'left',valign=>'top',size=>10,border=>1,bold=>1),
				identHid =>			$workbook->add_format(align=>'left',valign=>'top',size=>10,border=>1),
				identVisBC =>			$workbook->add_format(align=>'left',valign=>'top',size=>10,border=>1,bold=>1,color=>'grey'),
				identHidBC =>			$workbook->add_format(align=>'left',valign=>'top',size=>10,border=>1,color=>'grey'),
				identVisVP =>			$workbook->add_format(align=>'left',valign=>'top',size=>10,border=>1,bold=>1,italic=>1),
				identHidVP =>			$workbook->add_format(align=>'left',valign=>'top',size=>10,border=>1,italic=>1),
#mergeIdent =>		$workbook->add_format(align=>'left',valign=>'top',size=>10,border=>1),
				identBC =>			$workbook->add_format(align=>'left',size=>10,color=>'grey',border=>1),
#mergeIdentBC =>		$workbook->add_format(align=>'left',valign=>'top',size=>10,color=>'grey',border=>1),
				identVP =>			$workbook->add_format(align=>'left',size=>10,italic=>1,border=>1),
#mergeIdentVP =>		$workbook->add_format(align=>'left',valign=>'top',size=>10,italic=>1,border=>1),
				text =>				$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>0,border=>1),
				textWrap =>				$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>1,border=>1),
				mergeText =>		$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>0,border=>1),
				mergeColText =>		$workbook->add_format(align=>'left',size=>10,text_wrap=>0,border=>1),
				numberVis =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,bold=>1),
				numberHid =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1),
				numberVisBC =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,bold=>1,color=>'grey'),
				numberHidBC =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,color=>'grey'),
				numberVisVP =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,bold=>1,italic=>1),
				numberHidVP =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,italic=>1),
				number =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1),
				mergeNumber =>		$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1),
				number1d =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',num_format=>'0.0',border=>1),
#mergeNumber1d =>	$workbook->add_format(size=>10,align=>'center',valign=>'top',num_format=>'0.0',border=>1)
				);

		my $worksheet=$workbook->add_worksheet('Results');

		my $grColSpan=($numGroups)? 3 : 2;
		my $itemColSpan=$numAnaGroups*$grColSpan;
		my $lastColIndex=$itemColSpan+5;

		####<Start printing>####
		my $xlsRow=0;
		my $xlsCol=0;
		$worksheet->set_column(0,0,30); # identifier col length

		##<Title
		$worksheet->merge_range($xlsRow,0,$xlsRow,$lastColIndex,"$projectVarMods{$selModifID}[0]-Sites Comparison ($timeStamp2)",$itemFormat{'title'});
		##<List filter
		my $catFilterStrg;
		if ($catFilter) {
			if ($catFilterRule eq 'exclude') {$catFilterStrg="Proteins not in Theme '$filterClassName' > List '$filterCatName'";}
			else {$catFilterStrg="Proteins restricted to Theme '$filterClassName' > List '$filterCatName'";}
		}
		else {$catFilterStrg='No List restriction';}
		$worksheet->merge_range(++$xlsRow,0,$xlsRow,$lastColIndex,$catFilterStrg,$itemFormat{'mergeColText'});
		##<Modif Residue rules
		my $targetResStrg=($targetedRes)? "Target residue restricted to $targetedRes" : "No target residue restriction";
		$worksheet->merge_range(++$xlsRow,0,$xlsRow,$lastColIndex,$targetResStrg,$itemFormat{'mergeColText'});
		##<Excluded peptides
		my $numBadPep=scalar keys %unlocalizedPep;
		if ($numBadPep) {
			my $numExcludProt=scalar keys %excludedProteins;
			$worksheet->merge_range(++$xlsRow,0,$xlsRow,$lastColIndex,"WARNING: $numBadPep peptide sequence(s) could not be correctly aligned. $numExcludProt proteins were skipped. Modification sites number is underestimated.",$itemFormat{'mergeColText'});
		}
		##<Table header
		$xlsRow++;
		my $proteinStrg=($exportType eq 'all')? 'All site' : 'Selected sites';
		$worksheet->merge_range($xlsRow,0,$xlsRow+2,0,$proteinStrg,$itemFormat{'mergeRowHeader'});
		$worksheet->merge_range($xlsRow,1,$xlsRow+2,1,'Gene',$itemFormat{'mergeRowHeader'});
		$worksheet->merge_range($xlsRow,2,$xlsRow+2,2,'Context',$itemFormat{'mergeRowHeader'});
		my $siteStrg=($numGroups)? 'Sites distribution in groups' : 'Sites distribution in Analyses';
		$xlsCol=$itemColSpan+2;
		$worksheet->merge_range($xlsRow,3,$xlsRow,$xlsCol,$siteStrg,$itemFormat{'mergeColHeader'});
		$worksheet->merge_range($xlsRow,++$xlsCol,$xlsRow+2,$xlsCol,'MW (kDa)',$itemFormat{'mergeRowHeader'});
		$worksheet->set_column(++$xlsCol,$xlsCol,80); # col length
		$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow+2,$xlsCol,'Description',$itemFormat{'mergeRowHeader'});
		$worksheet->set_column(++$xlsCol,$xlsCol,30); # col length
		$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow+2,$xlsCol,'Species',$itemFormat{'mergeRowHeader'});
		#<Groups composition / analyses
		$xlsRow++;
		$xlsCol=3;
		my $maxItemInGr=1;
		foreach my $g (0..$#analysisGroups) {
			my $numItemInGr=scalar @{$analysisGroups[$g]};
			$maxItemInGr=$numItemInGr if $numItemInGr > $maxItemInGr;
		}
		$worksheet->set_row($xlsRow,15*(1+$maxItemInGr)) if $maxItemInGr > 1; # row height
		foreach my $g (0..$#analysisGroups) {
			my $labelStrg=($numGroups)? "$groupLabels[$g][0]:\n$groupLabels[$g][1]" : $groupLabels[$g][1];
			my $endCol=$xlsCol+$grColSpan-1;
			$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow,$endCol,decode_utf8($labelStrg),$itemFormat{'mergeColHeader'});
			$xlsCol=$endCol+1;
		}
		#<Occurence,Spectral count & site probability
		$xlsRow++;
		$xlsCol=2;
		foreach my $g (0..$#analysisGroups) {
			$worksheet->write_string($xlsRow,++$xlsCol,'Occ. %',$itemFormat{'header'}) if $numGroups;
			$worksheet->write_string($xlsRow,++$xlsCol,'Sp. Count',$itemFormat{'header'});
			$worksheet->write_string($xlsRow,++$xlsCol,'Prob. %',$itemFormat{'header'});
		}

		###<Looping through sites>###
		my $siteCount=0;
		foreach my $site (sort{&sortModifSite($sortOrder,\%proteinSites,\%protAnnotation,\%groupDistrib,\%siteDistrib)} keys %proteinSites) {
			my ($protID,$siteRes,$protPos,$refPRS,$refProb)=@{$proteinSites{$site}};
			#next if ($exportType eq 'sel' && !$selectedProteins{$protID} && !$selectedProteins{$site});
#print "$proteinSeq{$protID}<BR>\n$site: '",join('',@{$siteContext{$site}}[0..2]),"'<BR>\n";
			$xlsCol=0;
			##<Protein-Site
			my $visTag=($bestVisibility{$protID})? 'Vis' : 'Hid';
			my $dispProtSite=($siteRes eq '-')? "$projectVarMods{$selModifID}[1]~$protAnnotation{$protID}[0]" : ($siteRes eq '+')? "$protAnnotation{$protID}[0]~$projectVarMods{$selModifID}[1]" : "$protAnnotation{$protID}[0]-$siteRes$protPos";
			$worksheet->write_string(++$xlsRow,$xlsCol,$dispProtSite,$itemFormat{'ident'.$visTag}); #.$confTag
			##<Gene
			my $geneName=($protAnnotation{$protID}[4] && $masterProteins{$protAnnotation{$protID}[4]}[0])? $masterProteins{$protAnnotation{$protID}[4]}[0] : '-';
			$worksheet->write_string($xlsRow,++$xlsCol,$geneName,$itemFormat{'ident'.$visTag});
			##<Context
			my $context=join('',@{$siteContext{$site}}[0..2]);
			if ($sitesInProtein{$protID}{$site} eq '-') {}
			elsif ($sitesInProtein{$protID}{$site} eq '+') {}
			else {substr($context,$MID_CONTEXT,1,lc(substr($context,$MID_CONTEXT,1)));} # set middle AA to lower case
			$worksheet->write_string($xlsRow,++$xlsCol,$context,$itemFormat{'ident'.$visTag});
			##<Analyses/groups
			foreach my $g (0..$#analysisGroups) {
				if ($groupDistrib{$site}{$g}[2]) { # occ%
					my $gVis=($groupDistrib{$site}{$g}[0])? 'Vis' : 'Hid';
					my $gConf=($groupDistrib{$site}{$g}[1]==0)? 'VP' : ($groupDistrib{$site}{$g}[1]==1)? 'BC' : '';
					#<Occ%
					$worksheet->write_number($xlsRow,++$xlsCol,$groupDistrib{$site}{$g}[2],$itemFormat{'number'.$gVis.$gConf}) if $numGroups;
					#<Sp. Count
					$worksheet->write_number($xlsRow,++$xlsCol,$groupDistrib{$site}{$g}[4],$itemFormat{'number'.$gVis.$gConf});
					#<Probability
					my $probStrg='-';
					if ($refPRS->{$g}) {
						my ($status,$proba,$positions) = split(/;/,$refPRS->{$g}[0]);
						$probStrg='PRS:'.(1*(sprintf("%.1f",$proba)));
					}
					elsif (defined $refProb->{$g}) {$probStrg='MQ:'.($refProb->{$g}*100);}
					$worksheet->write_string($xlsRow,++$xlsCol,$probStrg,$itemFormat{'ident'.$gVis.$gConf});
				}
				else {
					++$xlsCol;
					$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow,$xlsCol+$grColSpan-1,'-',$itemFormat{'mergeNumber'});
					$xlsCol+=($grColSpan-1);
				}
			}
			##<MW
			if ($protAnnotation{$protID}[2]) {
				$worksheet->write_number($xlsRow,++$xlsCol,$protAnnotation{$protID}[2],$itemFormat{'number1d'});
			}
			else {
				$worksheet->write_string($xlsRow,++$xlsCol,'-',$itemFormat{'text'});
			}
			##<Description
			$worksheet->write_string($xlsRow,++$xlsCol,$protAnnotation{$protID}[1],$itemFormat{'text'});
			##<Species
			$worksheet->write_string($xlsRow,++$xlsCol,$protAnnotation{$protID}[3],$itemFormat{'text'});

			$siteCount++;
			if ($siteCount > 5000) {
				&updateProgress('.','+=');
				$siteCount=0;
			}
		}
		$workbook->close();
		&updateProgress(' Done.','+=');
		print qq
|<SCRIPT LANGUAGE="JavaScript">
document.getElementById('waitDIV').style.display='none';
</SCRIPT>
<CENTER>
<BR><BR><BR>
<INPUT type="button" class="title2" value="Download Worksheet" onclick="window.location='$promsPath{tmp_html}/scratch/export/$jobDir/$xlsFileName'"/>
</CENTER>
</BODY>
</HTML>
|;
		exit;
	} # end of export

	my $numAllSites=scalar keys %groupDistrib;
	if ($numAllSites==0) {
		print qq
|<SCRIPT type="text/javascript">document.getElementById('waitDIV').style.display='none';</SCRIPT>
<BR><BR><FONT class="title2">No sites matched!</FONT>
</CENTER>
</BODY>
</HTML>
|;
		exit;
	}

	####<Computing Venn Diagram>####
	my $numSharedByAll=0;
	my %vennDiagram;
	if ($numAnaGroups <= 5) {
		foreach my $site (keys %groupDistrib) {
			my $set='';
			foreach my $g (0..$#groupLabels) {$set.=chr(65+$g) if $groupDistrib{$site}{$g}[1];}
			$vennDiagram{$set}++;
		}
	}
	else {
		####<Global statistics>####
		foreach my $site (keys %groupDistrib) {
			my $okAllGr=1;
			foreach my $g (0..$#groupLabels) {
				unless ($groupDistrib{$site}{$g}[1]) {
					$okAllGr=0;
					last;
				}
			}
			$numSharedByAll+=$okAllGr;
		}
	}
	print qq
|<SCRIPT type="text/javascript">
document.getElementById('waitDIV').style.display='none';
function whenWindowLoaded() {
|;
	if ($numAnaGroups <= 5) {
		print qq
|	document.getElementById('masterVennDIV').style.display='';
	var groupLabelsObj={|;
		foreach my $g (0..$#groupLabels) {
			print chr(65+$g),":'$groupLabels[$g][0]'";
			print ',' if $g < $#analysisGroups;
		}
		print "};\n\tvar VD=new VennDiagram({div:'vennDIV',size:300,clickAction:updateSiteList,groupLabels:groupLabelsObj,setValues:{";
		my $numSets=scalar keys %vennDiagram;
		my $count=0;
		foreach my $set (sort keys %vennDiagram) {
			print "$set:$vennDiagram{$set}";
			$count++;
			print ',' if $count < $numSets;
		}
		print "}});\n";
	}
	my $onloadStrg.=($okAutocheck)? 'autoCheck(0); hideShowUnchecked(0);' : '';
	print qq
|	$onloadStrg
}
function updateSiteList(set) { // called when a set is selected in Venn diagram
	/* Clear all selections */
	resetSelection();
	/* Auto selection */
	if (set=='all') {return;}
	else if (set.match(':full')) {
		var setInfo=set.split(':');
		var a=set.charCodeAt(0)-65; // ascii code
		document.getElementById('numPep_'+a).selectedIndex=1; // 'Found'
	}
	else {
		var labelCode=['A','B','C','D','E'];
		for (var a=0; a < $numAnaGroups; a++) {
			document.getElementById('numPep_'+a).selectedIndex=(set.match(labelCode[a]))? 1 : 2; // 'Found' : 'Not found'
		}
	}
	autoCheck(0);
	hideShowUnchecked(0);
}
function resetSelection() {
	var myForm=document.protView;
	document.getElementById('logic').selectedIndex=0; // logic = 'and'
	updateComparisonStrategy('and');
	uncheckAll(myForm.chkProt);
	//hideShowUnchecked(0);
	for (var a=0; a < $numAnaGroups; a++) {
		document.getElementById('numPep_'+a).selectedIndex=0; // 'Any'
	}
}
</SCRIPT>
<DIV id="masterVennDIV" style="text-align:center;display:none">
	<INPUT type="button" value="Export as image" style="font-weight:bold;font-size:10px" onclick="exportSVGtoImg('vennDIV','VennDiagram','./exportSVG.cgi')"/>
	<DIV id="vennDIV" style="text-align:center;padding:3px;"></DIV>
</DIV>
|;

	####<Printing Store form>####
	&printStoreForm;
	my $numBadPep=scalar keys %unlocalizedPep;
	if ($numBadPep) {
		my $numExcludProt=scalar keys %excludedProteins;
		print "<FONT class=\"title3\" color=\"#DD0000\">WARNING: $numBadPep peptide sequence(s) could not be correctly aligned. $numExcludProt proteins were skipped. Modification sites number is underestimated.</FONT><BR>\n";
	}
	####<Listing sites comparison>####
	my $itemStrg=($numGroups)? 'groups' : 'Analyses';
	print "<FONT class=\"title2\">",scalar keys %groupDistrib," sites found. $numSharedByAll sites shared by all $itemStrg</FONT>\n" if $numAnaGroups > 5;
	print qq
|<FORM name="protView" method="post" style="margin:0">
<!--parameters for Export -->
<INPUT type="hidden" name="ACT" value="list">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="parentItem" value="$parentItem">
<INPUT type="hidden" name="compType" value="$compType">
<INPUT type="hidden" name="numGroups" value="$numGroups">
<INPUT type="hidden" name="catFilter" value="$catFilter">
<INPUT type="hidden" name="catFilterRule" value="$catFilterRule">
<INPUT type="hidden" name="pepRule" value="$peptideRule">
<INPUT type="hidden" name="pepThreshold" value="$pepThreshold">
<INPUT type="hidden" name="noHidden" value="$noHidden">
<INPUT type="hidden" name="hiddenRule" value="$hiddenRule">
<INPUT type="hidden" name="comparisonID" value="$comparisonID">
<INPUT type="hidden" name="virtualData" value="$virtualData">
<INPUT type="hidden" name="autocheck" value="">
<INPUT type="hidden" name="hideShowStatus" value="$hideShowStatus">
<INPUT type="hidden" name="sort" value="$sortOrder">
<INPUT type="hidden" name="anaList" value="$anaList">
<INPUT type="hidden" name="exportType" value="">
<!-- No need for peptide modification filters -->
<!-- Modif site-specific parameters -->
<INPUT type="hidden" name="siteRes" value="$targetedRes">
<INPUT type="hidden" name="siteMeas" value="$siteMeasure">
<!--parameters for Peptide Distribution -->
<INPUT type="hidden" name="what" value="checkbox">
<INPUT type="hidden" name="item" value="$item">
<INPUT type="hidden" name="id_item" value="$itemID">
<INPUT type="hidden" name="largWindow">
<INPUT type="hidden" name="graphView">

<TABLE border=0 cellspacing=0 cellpadding=2><TR><TD><DIV id="saveProtDIV" style="display:none;"></DIV></TD></TR></TABLE>

<TABLE border=0 cellspacing=0 cellpadding=2>
<TR><TH class="rbBorder" bgcolor=$color2 class="font11" rowspan=2 colspan=2 valign=top nowrap>&nbsp;<FONT class="title3">Auto-check :</FONT><BR>(<SPAN id=\"numProtChecked\">$numProtChecked</SPAN> checked)</TH>
|;
	my $freqStrg='';
	foreach my $g (0..$#analysisGroups) {
		print "<TH class=\"rbBorder\" bgcolor=$color2 nowrap><FONT onmouseover=\"popup('<B>$groupLabels[$g][1]</B>')\" onmouseout=\"popout()\">&nbsp;$groupLabels[$g][0]</FONT>&nbsp;</TH>\n";
		my $gNum=$g+1;
		$freqStrg.="<OPTION value=\"$gNum\"";
		$freqStrg.=' selected' if ($autocheckParams{'FREQ_ANAGR'} && $autocheckParams{'FREQ_ANAGR'}==$gNum);
		$freqStrg.=">$gNum</OPTION>";
	}
	my $disSaveCompStrg=($projectStatus > 0 || $projectAccess eq 'guest')? 'disabled' : ''; # project was ended/archived
	my $disSaveProtStrg=($projectStatus > 0 || $projectAccess eq 'guest')? 'disabled' : ''; # project was ended/archived
	my ($selOrStatus,$selFreqStatus)=($filterLogic eq 'or')? (' selected','') : ($filterLogic eq 'frequency')? ('',' selected') : ('','');
	my $thTableColSpan=($numAnaGroups==2 && !$listComparison)? 4 : 2;
	print qq
|<TH class="rbBorder" rowspan=2 colspan=$thTableColSpan height=100% valign=top><TABLE border=0 cellpadding=0 cellspacing=0 height=100%><TR>
	<TH class="rBorder" bgcolor=$color2 nowrap valign="top">
&nbsp;Frequency:&nbsp;<BR>&nbsp;<SELECT name="freqAnaGr" id="freqAnaGr" class="font11" disabled>$freqStrg</SELECT>/$numAnaGroups&nbsp;
	</TH>
	<TH class="rBorder" bgcolor=$color2 nowrap>
&nbsp;Strategy<SUP onmouseover="popup('Select<BR>&nbsp;&nbsp;&nbsp;-A <B>logic</B> between Group filters: <B>And</B> (intersection) or <B>Or</B> (union)<BR>&nbsp;&nbsp;&nbsp;or<BR>&nbsp;&nbsp;&nbsp;-a global <B>frequency</B> of detection')\" onmouseout=\"popout()\">*</SUP>:
<SELECT name="logic" id="logic" class="font11" onchange="updateComparisonStrategy(this.value)"><OPTION value="and">And</OPTION><OPTION value="or"$selOrStatus>Or</OPTION><OPTION value="frequency"$selFreqStatus>Frequency</OPTION></SELECT>&nbsp;
<BR>&nbsp;<INPUT type="button" value="Check Matching" class="font11" style="width:110px" onclick="autoCheck(1)">
<INPUT type="button" value="Uncheck All" class="font11" style="width:85px" onclick="uncheckAll(document.protView.chkProt)">
<INPUT type="button" name="hideShowButton" value="$hideShowStrg" class="font11" style="width:110px" onclick="hideShowUnchecked(1)">&nbsp;
	</TH>
	<TH nowrap valign=top><INPUT type="button" value="Save Comparison..." style="width:150px" onclick="showSaveCompForm('show')" $disSaveCompStrg><BR>
<INPUT type="button" id="saveFormBUTTON" value="Save Selected..." style="width:150px" onclick="ajaxManageSaveProteins('getThemes');" $disSaveProtStrg></TH>
	<TH nowrap valign=top><INPUT type="button" value="Export Selected" onclick="exportSites('sel')" style="width:120px"><BR><INPUT type="button" value="Export All" onclick="exportSites('all')" style="width:120px"></TD>
	<!--<TH nowrap valign=top><INPUT type="button" value="Peptide Distribution" style="width:150px" onclick="graphicalView()"></TH>-->
</TR></TABLE></TH></TR>
<TR>
|;
	foreach my $g (0..$#analysisGroups) {
		my $aGr=$g+1;
		#<Found / Not found
		my $selInfPep=($autocheckParams{'PEP'} && $autocheckParams{'PEP'}{$aGr}[0] eq '<=')? ' selected' : '';
		my $selFound=($autocheckParams{'PEP'} && $autocheckParams{'PEP'}{$aGr}[1] >= 1)? ' selected' : '';
		my $selNotFound=($selInfPep && $autocheckParams{'PEP'}{$aGr}[1]==0)? ' selected' : '';
		print qq
|<TH class="rbBorder" bgcolor=$color2 nowrap>&nbsp;<SELECT name="numPep_$g" id="numPep_$g" class="font11"><OPTION value="0">Any</OPTION><OPTION value="1"$selFound>Found</OPTION><OPTION value="-1">Not found</OPTION></SELECT>&nbsp;</TH>
|;
	}
	print qq
|</TR>
<TR><TH style="font-size:4px">&nbsp;</TH></TR>
|;

	my $proteinText="&nbsp;".(scalar keys %groupDistrib)." sites&nbsp;<BR>&nbsp;".(scalar keys %protAnnotation)." proteins&nbsp;";
	my $massText='MW <FONT class="small">(kDa)</FONT>';
	my $organismText='Species';
	my $occurenceText=($siteMeasure eq 'specAbs')? 'Spectral count' : ($siteMeasure eq 'specRel')? 'Spectral count (%)' : 'Occurence (%)';
	foreach my $g (0..$#analysisGroups) {
		$groupLabels[$g][0]="<FONT color=#DD0000>$groupLabels[$g][0]</FONT>" if $sortOrder eq 'group:'.($g+1);
	}
	if ($sortOrder eq 'protein') {
		$proteinText="<FONT color=#DD0000>$proteinText</FONT>";
	}
	elsif ($sortOrder eq 'mass') {
		$massText="<FONT color=#DD0000>$massText</FONT>";
	}
	elsif ($sortOrder eq 'organism') {
		$organismText="<FONT color=#DD0000>$organismText</FONT>";
	}
	elsif ($sortOrder eq 'peptide'){
		$occurenceText="<FONT color=#DD0000>$occurenceText</FONT>";
	}
	print qq
|<TR class="row_0">
<TH class="rbBorder" rowspan=2><A href="javascript:selectSort('protein')" onmouseover="popup('Click to sort sites by <B>ascending proteins name</B>.')" onmouseout="popout()">$proteinText</A>&nbsp;<!--<INPUT type="button" id="manChkBUT" value="Update Checked" class="font11" style="width:110px;\$manualChkVis" onclick="update2ndVennDiagram()" disabled>&nbsp;--></TH>
<TH class="rbBorder" rowspan=2>&nbsp;Gene&nbsp;</TH>
<TH class="rbBorder" colspan=$numAnaGroups><A href="javascript:selectSort('peptide')" onmouseover="popup('Click to sort sites by <B>best distribution</B>.')" onmouseout="popout()">$occurenceText</TH>
<TH class="rbBorder" width=75 rowspan=2><A href="javascript:selectSort('mass')" onmouseover="popup('Click to sort sites by <B>decreasing protein mass</B>.')" onmouseout="popout()">$massText</A></TH>
<TH class="bBorder" width=650 rowspan=2>Description - <A href="javascript:selectSort('organism')" onmouseover="popup('Click to sort sites by <B>ascending proteins species</B>.')" onmouseout="popout()">$organismText</A></TH>
</TR>
<TR class="row_0">
|;
	foreach my $g (0..$#analysisGroups) {
		my $numAnaStrg=($numGroups)? scalar @{$analysisGroups[$g]}.' items:<BR>' : '';
		print "<TH class=\"rbBorder\">&nbsp;<A href=\"javascript:selectSort('group:",$g+1,"')\" onmouseover=\"popup('<B>$numAnaStrg$groupLabels[$g][1]</B><BR>Click to sort proteins by <B>decreasing $occurenceText</B>.')\" onmouseout=\"popout()\">$groupLabels[$g][0]</A>&nbsp;</TH>\n";
	}
	print "</TR>\n";

	my $dataIdx=($siteMeasure eq 'occ')? 2 : 4;
	my $numSites=0;
	foreach my $site (sort{&sortModifSite($sortOrder,\%proteinSites,\%protAnnotation,\%groupDistrib,\%siteDistrib)} keys %proteinSites) {
		my ($protID,$siteRes,$protPos,$refPRS,$refProb)=@{$proteinSites{$site}};
		my $chkStatus=($selectedProteins{$protID} || $selectedProteins{$site})? ' checked' : '';
		my $displayStatusStrg=($chkStatus)? '' : " style=\"display:$hideShowStatus\"";
		my $trClass='row_'.(++$numSites % 2);
		my $chkBoxValue=join(':',keys %{$distinctAna{$protID}},$site);
		my $distAnaStrg=join(',',keys %{$distinctAna{$protID}});
		my ($protTag,$classStrg)=($bestVisibility{$protID})? ('TH','') : ('TD',' class="hiddenProt"'); # $classStrg is important for proper recording to custom Lists!!!

		#<Occurence (freq)/Spectral count computation
		my $dataString='';
		my @distrib;
		my $pepIdStrg='';
		foreach my $g (0..$#analysisGroups) {
			my ($freqTag,$occ)=($groupDistrib{$site}{$g}[0])? ('TH',1) : ('TD',0); # protein best visibility in gr
			my $prsStrg=($occ && $refPRS->{$g})? '&nbsp;'.&phosphoRS::printIcon($refPRS->{$g}[0],{format=>'text',posShift=>$refPRS->{$g}[2],pepSeq=>$refPRS->{$g}[3]}) : '';
			$prsStrg.='&nbsp;'.&promsMod::MaxQuantProbIcon($refProb->{$g},{popupText=>'MaxQuant probablity='.($refProb->{$g}*100).'%'}) if defined $refProb->{$g};
			push @distrib,$occ;
			$pepIdStrg.=',' if $g > 0;
			#$pepIdStrg.="$g:";
			$pepIdStrg.=join(':',@{$peptideIDs{$site}{$g}}) if $peptideIDs{$site}{$g};
			$dataString.="<TD align=\"center\">&nbsp;<FONT class=\"$freqTag\">$groupDistrib{$site}{$g}[$dataIdx]</FONT>$prsStrg&nbsp;";
			$dataString.="<INPUT type=\"hidden\" name=\"pepDistrib\" value=\"".(join(':',@distrib))."\"/>" if $g==$#analysisGroups;
			$dataString.="</TD>";
		}

		my $dispProtSite=($siteRes eq '-')? "<FONT class=\"$projectVarMods{$selModifID}[2]\">$projectVarMods{$selModifID}[1]~</FONT>$protAnnotation{$protID}[0]" : ($siteRes eq '+')? "$protAnnotation{$protID}[0]<FONT class=\"$projectVarMods{$selModifID}[2]\">~$projectVarMods{$selModifID}[1]</FONT>" : "$protAnnotation{$protID}[0]<FONT class=\"$projectVarMods{$selModifID}[2]\">-$siteRes$protPos</FONT>";
		print qq
|<TR id="prot_$site" class="list $trClass" valign="top"$displayStatusStrg>
<TD nowrap><INPUT type="checkbox" name="chkProt"$classStrg value="$chkBoxValue" onclick="manualCheck(this.checked)"$chkStatus><A class="$protTag" href="javascript:sequenceView($protID,'$distAnaStrg')">$dispProtSite</A><IMG id="IMG_$site" src="$promsPath{images}/plus.gif" align="top" onclick="ajaxGetModifSiteDetails(event,'$site','$pepIdStrg')"></TD>
|;
		#<Gene strg
		my $masterProtID=$protAnnotation{$protID}[4];
		if ($masterProtID) {
			if (scalar @{$masterProteins{$masterProtID}} > 1) {
				print "<TH>&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-".(join('<BR>&nbsp;&nbsp;-',@{$masterProteins{$masterProtID}}[1..$#{$masterProteins{$masterProtID}}]))."</B>')\" onmouseout=\"popout()\">$masterProteins{$masterProtID}[0]</A>&nbsp;</TH>";
			}
			elsif ($masterProteins{$masterProtID}[0]) {print "<TH>&nbsp;$masterProteins{$masterProtID}[0]&nbsp;</TH>";}
			else {print "<TH>&nbsp;-&nbsp;</TH>";} # no gene mapped
		}
		else {print "<TH>&nbsp;-&nbsp;</TH>";} # no mapping at all

		#<Occurence (freq)/Spectral count display
		print $dataString;

		#<MW
		my $mass=($protAnnotation{$protID}[2]==0)? '-' : $protAnnotation{$protID}[2];
		print "<TD class=\"$protTag\" align=\"right\">$mass&nbsp;</TD>\n";
		#>Description - Species
		print "<TD width=600>$protAnnotation{$protID}[1] <FONT class=\"org\">$protAnnotation{$protID}[3]</FONT></TD>\n";
		print "</TR>\n";
	}
	print qq
|</TABLE>
<B><I>End of list.</I></B>
</FORM>
<BR><BR>
<DIV id="siteModifDetailsDIV" style="z-index:200;background-color:#FFFFFF;border:solid 3px #999999;padding:5px;position:absolute;display:none;"></DIV> <!--filter:alpha(opacity=80);opacity:0.8;-->
<DIV id="divDescription" class="clDescriptionCont"></DIV>
<SCRIPT type="text/javascript">
setPopup();
whenWindowLoaded();
</SCRIPT>
</BODY>
</HTML>
|;

	exit;
}

##########################################################################
#############<HEAD for 1 vs 1 Peptide/Protein comparisons>################
##########################################################################
sub display1v1ComparisonHead {
	if ($exportType) {
		my $timeStamp1=strftime("%Y%m%d %H-%M",localtime);
		print header(-type=>"application/vnd.ms-excel",-attachment=>"Compare items_$timeStamp1.xls");
	}
	else {
		my $fontSize=($numAnaGroups>=40)? '8px' : ($numAnaGroups>=25)? '10px' : ($numAnaGroups>=15)? '11px' : '12px';
		print header(-'content-encoding'=>'no',-charset=>'utf-8');
		warningsToBrowser(1);
		print qq
|<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.rowHeader{background-color:$color2; height:50px; text-align:right;}
.colHeader{background-color:$color2; min-width:75px;}
#compTable TH{font-size:$fontSize;}
</STYLE>
<SCRIPT type="text/javascript">
function preparePopup(item1,item2,ratio1,ratio2) {
	var popContent="<TABLE cellspacing=0><TR><TD nowrap><B>Common to "+item1+" and "+item2+":</B><BR>&nbsp;&bull;<B>"+ratio1+"% of "+item1+"</B><BR>&nbsp;&bull;<B>"+ratio2+"% of "+item1+" U "+item2+"</B></TD></TR></TABLE>";
	popup(popContent);
}
|;
		&promsMod::popupInfo();
		print qq
|</SCRIPT>
</HEAD>
<BODY style="margin:2px;" background="$promsPath{images}/bgProMS.gif">
<CENTER>
<DIV id="waitDIV">
<BR><BR><BR><BR><BR><FONT class="title3">Fetching data...</FONT><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">
<BR><BR><BR><FONT class="title3">Status: <SPAN id="waitSPAN"><SPAN>...</FONT>
</DIV>
|;
	}
}

##########################################################################
#############<BODY for 1 vs 1 Peptide/Protein comparisons>################
##########################################################################
sub display1v1ComparisonBody {
	my ($elementType,$refGroupLabels,$refElementDistrib,$refGroupDistrib)=@_;

	####<Computing Venn Diagram>####
	my %vennDiagram;
	if ($numAnaGroups <= 5) {
		foreach my $elem (keys %{$refGroupDistrib}) {
			my $set='';
			foreach my $g (sort{$a<=>$b} keys %{$refGroupDistrib->{$elem}}) {$set.=chr(65+$g);}
			$vennDiagram{$set}++;
		}
	}

	#######################
	####<Starting BODY>####
	#######################
	if ($numAnaGroups<=5) {
		print qq
|<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/vennDiagram.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
|;
	}
	print qq
|<SCRIPT type="text/javascript">
document.getElementById('waitDIV').style.display='none';
function whenWindowLoaded() {
|;
	if ($numAnaGroups <= 5) {
		print "\tvar groupLabelsObj={";
		foreach my $g (0..$#{$refGroupLabels}) {
			print chr(65+$g),":'$refGroupLabels->[$g][0]'";
			print ',' if $g < $#analysisGroups;
		}
		print "};\n\tvar VD=new VennDiagram({div:'vennDIV',size:300,groupLabels:groupLabelsObj,setValues:{";
		my $numSets=scalar keys %vennDiagram;
		my $count=0;
		foreach my $set (sort keys %vennDiagram) {
			print "$set:$vennDiagram{$set}";
			$count++;
			print ',' if $count < $numSets;
		}
		print "}});\n";
	}
	print qq
|}
</SCRIPT>
|;
	####<Printing protein list form (1v1_prot only)>####
	if ($elementType eq 'protein') {
		print qq
|<FORM name="listForm" method="post" target='listFrame' style="margin:0">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="parentItem" value="$parentItem">
<INPUT type="hidden" name="compType" value="full_prot">
<INPUT type="hidden" name="ACT" value="list">
<INPUT type="hidden" name="call" value="1v1_prot">
<INPUT type="hidden" name="numGroups" value="">
<INPUT type="hidden" name="groupPos" value="">
<INPUT type="hidden" name="catFilter" value="$catFilter">
<INPUT type="hidden" name="catFilterRule" value="$catFilterRule">
<INPUT type="hidden" name="anaList" value="">
<INPUT type="hidden" name="sort" value="$sortOrder">
<INPUT type="hidden" name="pepRule" value="$peptideRule">
<INPUT type="hidden" name="pepThreshold" value="$pepThreshold">
<INPUT type="hidden" name="virtualData" value="$virtualData">
<INPUT type="hidden" name="noHidden" value="$noHidden">
<INPUT type="hidden" name="hiddenRule" value="$hiddenRule">
<INPUT type="hidden" name="noMissCut" value="$noMissCut">
|;
		print "<INPUT type=\"hidden\" name=\"modifList\" value=\"",join(',',@modifList),"\">\n";
		if ($modifFilters{active}) {
			foreach my $filter ('exclude','restrict','ignore') {
				foreach my $modID (sort{$a<=>$b} keys %{$modifFilters{$filter}}) {print "<INPUT type=\"hidden\" name=\"filterMod_$modID\" value=\"$filter\">\n";}
			}
		}
		print qq
|<INPUT type="hidden" name="restrictModifLogic" value="$restrictModifLogic">
<INPUT type="hidden" name="delocPhospho" value="$delocPhospho">
<INPUT type="hidden" name="autocheck" value="">
	<!--INPUT type="hidden" name="protID" value="-"-->
</FORM>
<SCRIPT type="text/javascript">
var anaJoin=($numGroups)? ',' : '.';
var anaGroups=[];
|;
		foreach my $gr (0..$#analysisGroups) {
			print "anaGroups[$gr]='",join('.',@{$analysisGroups[$gr]}),"';\n";
		}
		print qq
|function listProteins(groups) { //(numGroups,groupPos,anaList,autocheckStrg)
	var setsUsed=[],autocheckStrg='and#';
	if (groups.length==2) { // shared proteins between 2 groups
		setsUsed=groups;
		autocheckStrg+='1:>=:1;2:>=:1';
	}
	else { // length=1 Proteins unique to a group
		for (var i=0; i<=$#analysisGroups; i++) {
			setsUsed.push(i);
			if (i>0) {autocheckStrg+=';';}
			autocheckStrg+=(i==groups[0])? (i+1)+':>=:1' : (i+1)+':<=:0';
		}
	}
	var anaListStrg='',groupPosStrg='';
	for (var i=0; i<setsUsed.length; i++) {
		if (anaListStrg) {
			anaListStrg+=anaJoin; // , or .
			groupPosStrg+=',';
		}
		anaListStrg+=anaGroups[setsUsed[i]];
		groupPosStrg+=setsUsed[i]+1; // Pos not Index
	}
	myForm=document.listForm;
	myForm.numGroups.value=($numGroups)? groups.length : 0;
	myForm.anaList.value=anaListStrg;
	if ($numGroups) {myForm.groupPos.value=groupPosStrg;} else {myForm.groupPos.value='';}  // Only for groups: string to use gr position in 1v1_prot compare (eg. 4 vs 5)
	myForm.autocheck.value=autocheckStrg;
	myForm.submit();
}
</SCRIPT>
|;
	}
	else {
		print qq
|<SCRIPT type="text/javascript">
function listProteins() {} // do nothing if peptide comparison
</SCRIPT>
|;
	}

	####<Printing Store form>####
	&printStoreForm;

	####<Global statistics>####
	my $numSharedByAll=0;
	foreach my $elem (keys %{$refGroupDistrib}) {
		$numSharedByAll++ if scalar keys %{$refGroupDistrib->{$elem}}==$numAnaGroups;
	}

	my $itemStrg=($numGroups)? 'groups' : 'Analyses';
	my $elemStrg=($elementType eq 'peptide')? 'distinct peptides' : 'proteins';
	my $disSaveStrg=($projectStatus > 0 || $projectAccess eq 'guest')? 'disabled' : ''; # project was ended/archived
	print "<FONT class=\"title2\">",scalar keys %{$refGroupDistrib}," $elemStrg found. $numSharedByAll $elemStrg shared by all $itemStrg</FONT>\n";
	print "<INPUT type=\"button\" value=\"Save Comparison...\" style=\"width:150px\" onclick=\"showSaveCompForm('show')\" $disSaveStrg><BR>\n";

	####<Drawing peptide/protein comparison>####
	my (%numElementsInGroup,%elementShared,%elementUnique);
	my ($maxNumPep,$maxNumUnique)=(0,0);
	print "<BR>\n";
	print "<TABLE cellspacing=10><TR><TD>\n" if $numAnaGroups <= 5;
	print "<FONT class=\"font11\">(Click in green/red cells to display detailed comparison in main window)</FONT>\n" if ($elementType eq 'protein');
	print "<TABLE id=\"compTable\" border=0 cellpadding=5>\n</TR><TH></TH><TH class=\"colHeader\">Unique</TH>\n";
	foreach my $g (0..$#analysisGroups) {
		print "<TH class=\"colHeader\"><A href=\"javascript:void(null)\" onmouseover=\"popup('<B>$refGroupLabels->[$g][1]</B>')\" onmouseout=\"popout()\">$refGroupLabels->[$g][0]</A></TH>";
		$numElementsInGroup{$g}=scalar keys %{$refElementDistrib->{$g}};
		$maxNumPep=$numElementsInGroup{$g} if $numElementsInGroup{$g} > $maxNumPep;
		$elementUnique{$g}=0;
		foreach my $elem (keys %{$refElementDistrib->{$g}}) {
			$elementUnique{$g}++ if scalar keys %{$refGroupDistrib->{$elem}}==1;
		}
		$maxNumUnique=$elementUnique{$g} if $elementUnique{$g} > $maxNumUnique;
	}
	foreach my $g1 (0..$#analysisGroups) {
		print "<TR><TH class=\"rowHeader\" nowrap><A href=\"javascript:void(null)\" onmouseover=\"popup('<B>$refGroupLabels->[$g1][1]</B>')\" onmouseout=\"popout()\">$refGroupLabels->[$g1][0]</A></TH>";
		#my $ratioColorU=($maxNumUnique)? sprintf '%.0f',240*(1-($elementUnique{$g1}/$maxNumUnique)) : 240;
		#print "<TH style=\"background-color:rgb($ratioColorU,240,$ratioColorU)\">$elementUnique{$g1}</TH>";
		my $ratioColorU=($maxNumUnique)? sprintf '%.0f',50+205*(1-($elementUnique{$g1}/$maxNumUnique)) : 255;
		print "<TH style=\"background-color:rgb($ratioColorU,255,$ratioColorU)\"><A href=\"javascript:listProteins([$g1])\" onmouseover=\"popup('<B>Unique to $refGroupLabels->[$g1][0]</B>')\" onmouseout=\"popout()\">$elementUnique{$g1}</A></TH>";
		foreach my $g2 (0..$#analysisGroups) {
			if ($g2==$g1) {
				my $ratioColor=($maxNumPep)? sprintf '%.0f',150+105*(1-($numElementsInGroup{$g1}/$maxNumPep)) : 255;
				print "<TH style=\"background-color:rgb($ratioColor,$ratioColor,$ratioColor)\"><A href=\"javascript:void(null)\" onmouseover=\"popup('<B>$refGroupLabels->[$g1][0]</B>')\" onmouseout=\"popout()\">$numElementsInGroup{$g1}</A></TH>";
			}
			else { # compute intersection & union
				unless (defined $elementShared{"$g1:$g2"}) {
					my $inter=0;
					foreach my $elem (keys %{$refElementDistrib->{$g1}}) {
						$inter++ if $refElementDistrib->{$g2}{$elem};
					}
					$elementShared{"$g1:$g2"}=$elementShared{"$g2:$g1"}=$inter;
				}
				my ($ratio0,$ratio1,$ratio2);
				if ($numElementsInGroup{$g1}) {
					$ratio0=$elementShared{"$g1:$g2"}/$numElementsInGroup{$g1};
					$ratio1=1 * sprintf '%.1f',100*$ratio0;
					$ratio2=1 * sprintf '%.1f',100*($elementShared{"$g1:$g2"}/($numElementsInGroup{$g1}+$numElementsInGroup{$g2}-$elementShared{"$g1:$g2"}));
				}
				else {$ratio0=$ratio1=$ratio2=0;}
				#my $ratioColorS=sprintf '%.0f',230*(1-$ratio0);
				#print "<TH style=\"background-color:rgb(230,$ratioColorS,$ratioColorS)\" onmouseover=\"preparePopup('$refGroupLabels->[$g1][0]','$refGroupLabels->[$g2][0]',$ratio1,$ratio2)\" onmouseout=\"popout()\">",$elementShared{"$g1:$g2"},"</TH>";
				my $ratioColorS=sprintf '%.0f',100+155*(1-$ratio0);
				print "<TH style=\"background-color:rgb(255,$ratioColorS,$ratioColorS)\"><A href=\"javascript:listProteins([$g1,$g2])\" onmouseover=\"preparePopup('$refGroupLabels->[$g1][0]','$refGroupLabels->[$g2][0]',$ratio1,$ratio2)\" onmouseout=\"popout()\">",$elementShared{"$g1:$g2"},"</A></TH>";
			}
		}
		print "</TR>\n";
	}
	print "</TABLE>\n"; # heatmap table
	if ($numAnaGroups <= 5) {
		print qq
|</TD>
<TH><INPUT type="button" value="Export as image" style="font-weight:bold;font-size:10px" onclick="exportSVGtoImg('vennDIV','VennDiagram','./exportSVG.cgi')"/>
<DIV id="vennDIV" style="text-align:center;padding:3px;"></DIV>
</TH>
</TR>
</TABLE>
|;
	}
	print qq
|</CENTER>
<BR><BR>
<DIV id="protDetailsDIV" style="z-index:200;background-color:#FFFFFF;border:solid 3px #999999;padding:5px;position:absolute;display:none;"></DIV> <!--filter:alpha(opacity=80);opacity:0.8;-->
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">
setPopup();
whenWindowLoaded();
</SCRIPT>
</BODY>
</HTML>
|;

	exit;
}

sub updateProteinData { # global: $peptideRule,$numGroups,$type,   %proteinList,%bestAnalysis,%listInGroup,%bestAnalysis,%peptideAnaDistrib,%protAnalyses,%processedProtAna
	my ($g,$codeIndex,$type,$anaID,$anaIndex,$protID,$vis,$confLevel,$score,$numPep,$refBestScore,$refGroupBestVis,$refItemAnaScanIdx)=@_;

	if (!$proteinList{$protID} || !$proteinList{$protID}[$g]) { # 1st time protein is seen
		@{$proteinList{$protID}[$g]}=($numPep,$vis,$confLevel);
		if ($peptideRule=~/ba/) { # ba or banr (best Ana or no groups)$bestScore{$protID}=$score
			$refBestScore->{$protID}=$score;
			$bestAnalysis{$protID}{$g}=$anaID;
		}
	}
	elsif (($numGroups || $listInGroup{$g}) && $peptideRule=~/ba/) { # ba or banr (best Ana)
		#if (($vis>=1 || $vis==$proteinList{$protID}[$g][1]) && ($numPep>$proteinList{$protID}[$g][0] || ($numPep==$proteinList{$protID}[$g][0] && $score>$refBestScore->{$protID}))) { #}
		if ((!$noHidden || $vis>=1 || $vis==$proteinList{$protID}[$g][1]) && ($numPep > $proteinList{$protID}[$g][0] || ($numPep==$proteinList{$protID}[$g][0] && $score > $refBestScore->{$protID}))) { # PP 25/10/16: Ignore hidden => ignore visibility
			@{$proteinList{$protID}[$g]}=($numPep,$vis,$confLevel);
			$refBestScore->{$protID}=$score;
			$bestAnalysis{$protID}{$g}=$anaID;
		}
	}
	else { # Group-level => combine all ana data
		$proteinList{$protID}[$g][0]+=$numPep;
		$proteinList{$protID}[$g][1]=$vis if $proteinList{$protID}[$g][1] < $vis;
		$proteinList{$protID}[$g][2]=$confLevel if $proteinList{$protID}[$g][2] < $confLevel;
	}
	if (!defined($refGroupBestVis->{$protID}) || $refGroupBestVis->{$protID}[1]<$vis || ($refGroupBestVis->{$protID}[1]==$vis && $refGroupBestVis->{$protID}[2]<$numPep)) { # best vis is used even if !noHidden (PP 25/10/16)
		@{$refGroupBestVis->{$protID}}=($anaID,$vis,$numPep);
		#$bestAnalysis{$protID}{$g}=$anaID;
	}

	if ($type eq 'A') {
		$peptideAnaDistrib{$protID}[$g][$codeIndex]=$protAnalyses{$protID}{$anaID}[0] if ($numGroups || $listInGroup{$g});
#print "ANA1 $protID:$g:$codeIndex:$codeIndex => $peptideAnaDistrib{$protID}[$g][$codeIndex]<BR>\n";
	}
	else { # C
		$peptideAnaDistrib{$protID}[$g][$codeIndex][$anaIndex]=$protAnalyses{$protID}{$anaID}[0]; # 1 more dimension!!!
#print "CAT1 $protID:$g:$codeIndex:$anaIndex => $peptideAnaDistrib{$protID}[$g][$codeIndex][$anaIndex]<BR>\n";
		$refItemAnaScanIdx->{$protID}=$codeIndex;
	}
	$processedProtAna{$protID}{$anaID}=1;
}

sub storeComparison {

	####<Connecting to the database>####
	my $dbh=&promsConfig::dbConnect;

	my $pepParamStrg="noHidden=$noHidden;hiddenRule=$hiddenRule;pepRule=$peptideRule;pepThreshold=$pepThreshold;virtualData=$virtualData;noMissCut=$noMissCut;delocPhospho=$delocPhospho;pepSpecificity=$pepSpecificity";
	if ($modifFilters{active}) {
		my $modFilterStrg='';
		foreach my $filter ('exclude','restrict','ignore') {
			if (scalar keys %{$modifFilters{$filter}}) {
				$modFilterStrg.=',' if $modFilterStrg;
				$modFilterStrg.=$filter.':'.join('.',(sort{$a<=>$b} keys %{$modifFilters{$filter}}));
				$modFilterStrg.=':'.$restrictModifLogic if $filter eq 'restrict';
			}
		}
		$pepParamStrg.=';modifFilters='.$modFilterStrg;
	}
	my ($catFilterID,$catExclusion);
	if ($catFilter==0) {$catFilterID=$catExclusion='NULL';}
	else {
		$catFilterID=$catFilter;
		$catExclusion=($catFilterRule eq 'exclude')? 1 : 'NULL'; # or 0?
	}

	if ($comparisonID) { # updating an existing comparison
		#my $compComments=(param('compComments'))? $dbh->quote(param('compComments')) : 'NULL';
		my $updateAutochkStrg=($autochkStrg eq '=')? '' : ",AUTOCHECK_PARAMS='$autochkStrg'"; # to prevent 1_vs_1 from overwriting full_prot selection
		my $compComments=$dbh->quote(param('compComments'));
		$dbh->do("UPDATE COMPARISON SET ID_CATEGORY=$catFilterID,CAT_EXCLUSION=$catExclusion,COMMENTS=$compComments,NUM_GROUPS=$numGroups,PEPTIDE_PARAMS='$pepParamStrg',SORT_ORDER='$sortOrder'$updateAutochkStrg,UPDATE_DATE=NOW(),UPDATE_USER='$userID' WHERE ID_COMPARISON=$comparisonID");
		$dbh->do("DELETE FROM ANA_COMPARISON WHERE ID_COMPARISON=$comparisonID");
		$dbh->do("DELETE FROM CAT_COMPARISON WHERE ID_COMPARISON=$comparisonID");
	}
	else {
		($comparisonID)=$dbh->selectrow_array("SELECT MAX(ID_COMPARISON) FROM COMPARISON");
		$comparisonID=0 unless $comparisonID;
		$comparisonID++;
		my $compName=$dbh->quote(param('newCompName'));
		#my $compComments=(param('newCompComments'))? $dbh->quote(param('newCompComments')) : 'NULL';
		my $compComments=$dbh->quote(param('newCompComments'));
		#my ($insAutochkStrg1,$insAutochkStrg2)=($autochkStrg eq '=')? ('','') : (',AUTOCHECK_PARAMS',",'$autochkStrg'");
		if (!$autochkStrg || $autochkStrg eq '=') { # generate one
			$autochkStrg='and#';
			foreach my $g (0..$#analysisGroups) {
				$autochkStrg.=';' if $g > 0;
				$autochkStrg.=($g+1).':>=:0';
			}
		}
		$dbh->do("INSERT INTO COMPARISON (ID_COMPARISON,ID_PROJECT,ID_CATEGORY,CAT_EXCLUSION,NAME,COMMENTS,NUM_GROUPS,PEPTIDE_PARAMS,SORT_ORDER,AUTOCHECK_PARAMS,UPDATE_DATE,UPDATE_USER) VALUES ($comparisonID,$projectID,$catFilterID,$catExclusion,$compName,$compComments,$numGroups,'$pepParamStrg','$sortOrder','$autochkStrg',NOW(),'$userID')");
	}
	my %sthType;
	$sthType{'A'}=$dbh->prepare("INSERT INTO ANA_COMPARISON (ID_COMPARISON,ID_ANALYSIS,COMP_GROUP,ANA_POS) VALUES ($comparisonID,?,?,?)");
	$sthType{'C'}=$dbh->prepare("INSERT INTO CAT_COMPARISON (ID_COMPARISON,ID_CATEGORY,COMP_GROUP,CAT_POS) VALUES ($comparisonID,?,?,?)");
	if ($numGroups) {
		foreach my $g (0..$#analysisGroups) {
			my $itemPos=0;
			foreach my $codeID (@{$analysisGroups[$g]}) {
				my ($type,$typeID)=split(':',$codeID); # A or C
				$sthType{$type}->execute($typeID,($g+1),++$itemPos);
			}
		}
	}
	else {
		foreach my $g (0..$#analysisGroups) { # 1 item/fake group
			my $codeID=$analysisGroups[$g][0];
			my ($type,$typeID)=split(':',$codeID); # A or C
			$sthType{$type}->execute($typeID,0,($g+1));
		}
	}
	$sthType{'A'}->finish;
	$sthType{'C'}->finish;
	$dbh->commit;
	$dbh->disconnect;

	####<Starting HTML>####
	print header(-charset=>'utf-8'); warningsToBrowser(1);
	my $parentFrame=($compType eq '1v1_prot')? 'opener' : 'parent.selItemsFrame';
	print qq
|<HEAD>
<SCRIPT type="text/javascript">
$parentFrame.location="./compareAnalyses.cgi?ACT=selAna&id_project=$projectID&parentItem=$parentItem&comparisonID=$comparisonID&compType=$compType";
</SCRIPT>
</HEAD>
</HTML>
|;
	exit;
}

sub sortProt {
	my $sortOrder=$_[0];
	# Alias asc.
	if ($sortOrder eq 'protein') {return lc($protAnnotation{$a}[0]) cmp lc($protAnnotation{$b}[0])}
	# Peptide in a given group/analysis desc. (anaPep > Alias)
	if ($sortOrder=~/group:(\d+)/) {my $grIdx=$1-1; return $proteinList{$b}[$grIdx][0]<=>$proteinList{$a}[$grIdx][0] || lc($protAnnotation{$a}[0]) cmp lc($protAnnotation{$b}[0])}
	# Distribution desc. (num gr/ana desc > delta ave asc > num peptides desc > alias)
	elsif ($sortOrder eq 'peptide') {$protDistribution{$b}[0]<=>$protDistribution{$a}[0] || $protDistribution{$a}[1]<=>$protDistribution{$b}[1] || $protDistribution{$b}[2]<=>$protDistribution{$a}[2] || lc($protAnnotation{$a}[0]) cmp lc($protAnnotation{$b}[0])}
	# Protein mass desc (mass > alias)
	elsif ($sortOrder eq 'mass') {$protAnnotation{$b}[2]<=>$protAnnotation{$a}[2] || lc($protAnnotation{$a}[0]) cmp lc($protAnnotation{$b}[0])}
	# Organism asc. (organism > alias)
	elsif ($sortOrder eq 'organism') {lc($protAnnotation{$a}[3]) cmp lc($protAnnotation{$b}[3]) || lc($protAnnotation{$a}[0]) cmp lc($protAnnotation{$b}[0])}
	# absolute Delta (absDelta > relDelta > num Peptides g#2 > alias)
	elsif ($sortOrder eq 'absDelta') {$absoluteDelta{$b}<=>$absoluteDelta{$a} || $relativeDelta{$b}<=>$relativeDelta{$a} || $proteinList{$b}[1][0]<=>$proteinList{$a}[1][0] || lc($protAnnotation{$a}[0]) cmp lc($protAnnotation{$b}[0])}
	# relative Delta (relDelta > absDelta > num Peptides g#2 > alias)
	elsif ($sortOrder eq 'relDelta') {$relativeDelta{$b}<=>$relativeDelta{$a} || $absoluteDelta{$b}<=>$absoluteDelta{$a} || $proteinList{$b}[1][0]<=>$proteinList{$a}[1][0] || lc($protAnnotation{$a}[0]) cmp lc($protAnnotation{$b}[0])}
}

sub sortModifSite {
	my ($sortOrder,$refProteinSites,$refProtAnnotation,$refGroupDistrib,$refSiteDistrib)=@_;
	# Alias asc > site pos asc
	if ($sortOrder eq 'protein') {return lc($refProtAnnotation->{$refProteinSites->{$a}[0]}[0]) cmp lc($refProtAnnotation->{$refProteinSites->{$b}[0]}[0]) || $refProteinSites->{$a}[2]<=>$refProteinSites->{$b}[2]}
	# Peptide in a given group/analysis desc. (anaPep > Alias > site pos asc)
	if ($sortOrder=~/group:(\d+)/) {
		my $grIdx=$1-1;
		my $idx=($siteMeasure=~/^spec/)? 4 : 2;
		return $refGroupDistrib->{$b}{$grIdx}[$idx]<=>$refGroupDistrib->{$a}{$grIdx}[$idx] || lc($refProtAnnotation->{$refProteinSites->{$a}[0]}[0]) cmp lc($refProtAnnotation->{$refProteinSites->{$b}[0]}[0]) || $refProteinSites->{$a}[2]<=>$refProteinSites->{$b}[2]}

	# Distribution desc. (num gr/ana desc > delta ave asc > alias > site pos asc)
	elsif ($sortOrder eq 'peptide') {
		return $refSiteDistrib->{$b}[0]<=>$refSiteDistrib->{$a}[0] || $refSiteDistrib->{$a}[1]<=>$refSiteDistrib->{$b}[1] || lc($refProtAnnotation->{$refProteinSites->{$a}[0]}[0]) cmp lc($refProtAnnotation->{$refProteinSites->{$b}[0]}[0]) || $refProteinSites->{$a}[2]<=>$refProteinSites->{$b}[2]
	}

	# Protein mass desc (mass > alias > site pos asc)
	elsif ($sortOrder eq 'mass') {return $refProtAnnotation->{$refProteinSites->{$b}[0]}[2]<=>$refProtAnnotation->{$refProteinSites->{$a}[0]}[2] || lc($refProtAnnotation->{$refProteinSites->{$a}[0]}[0]) cmp lc($refProtAnnotation->{$refProteinSites->{$b}[0]}[0]) || $refProteinSites->{$a}[2]<=>$refProteinSites->{$b}[2]}
	# Organism asc. (organism > alias > site pos asc)
	elsif ($sortOrder eq 'organism') {return lc($refProtAnnotation->{$refProteinSites->{$a}[0]}[3]) cmp lc($refProtAnnotation->{$refProteinSites->{$b}[0]}[3]) || lc($refProtAnnotation->{$refProteinSites->{$a}[0]}[0]) cmp lc($refProtAnnotation->{$refProteinSites->{$b}[0]}[0]) || $refProteinSites->{$a}[2]<=>$refProteinSites->{$b}[2]}
}

sub getVarModCodes { # also in listProteins.cgi
	my ($refProtMods,$varModStrg,$peptide)=@_;
	$varModStrg=~s/^\s\+\s//;
	foreach my $varModPos (split(/ \+ /,$varModStrg)) {
		my ($varModCode,$resStrg)=&promsMod::convertVarModString($varModPos,1);
		next unless $varModCode;
		if ($projectVarMods{$varModCode}) { # must be in project list
			$refProtMods->{$varModCode}{'ALL'}++;
			$refProtMods->{$varModCode}{'NR'}{$peptide}=1 if $peptideRule=~/nr/; # prevents INTRA-analysis redundancy
		}
	}
}

##################################
##############<AJAX>##############
##################################
sub ajaxGetAnalysesList {
	my ($parentItem,$parentID)=split(':',param('branchID'));

	my @queryList;
	if ($parentItem eq 'EXPERIMENT') {
		my $expQS=qq |SELECT CONCAT('A:',A.ID_ANALYSIS),S.NAME,A.NAME
								FROM SAMPLE S,ANALYSIS A
								WHERE S.ID_EXPERIMENT=$parentID AND S.ID_SAMPLE=A.ID_SAMPLE AND S.ID_SPOT IS NULL AND A.VALID_STATUS>=1
								ORDER BY S.DISPLAY_POS ASC,A.DISPLAY_POS ASC|;
		my $expQG=qq |SELECT CONCAT('A:',A.ID_ANALYSIS),G.NAME,SP.NAME,A.NAME
								FROM GEL2D G,SPOT SP,SAMPLE S,ANALYSIS A
								WHERE G.ID_EXPERIMENT=$parentID AND G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT
								AND S.ID_SAMPLE=A.ID_SAMPLE AND A.VALID_STATUS>=1
								ORDER BY G.DISPLAY_POS ASC,SP.NAME ASC,A.DISPLAY_POS ASC|;
		push @queryList,($expQS,$expQG);
	}
	elsif ($parentItem eq 'GEL2D') {
		my $gelQuery=qq |SELECT CONCAT('A:',A.ID_ANALYSIS),SP.NAME,A.NAME
								FROM SPOT SP,SAMPLE S,ANALYSIS A
								WHERE SP.ID_GEL2D=$parentID AND SP.ID_SPOT=S.ID_SPOT
								AND S.ID_SAMPLE=A.ID_SAMPLE AND A.VALID_STATUS>=1
								ORDER BY SP.NAME ASC,A.DISPLAY_POS ASC|;
		push @queryList,$gelQuery;
	}
	#elsif ($parentItem eq 'SPOT') {
	#	push @queryList,"SELECT CONCAT('A:',A.ID_ANALYSIS),A.NAME FROM SAMPLE S,ANALYSIS A WHERE S.ID_SPOT=$parentID AND S.ID_SAMPLE=A.ID_SAMPLE AND A.VALID_STATUS>=1 ORDER BY A.DISPLAY_POS ASC";
	#}
	elsif ($parentItem eq 'SAMPLE') {
		push @queryList,"SELECT CONCAT('A:',ID_ANALYSIS),NAME FROM ANALYSIS WHERE ID_SAMPLE=$parentID AND VALID_STATUS>=1 ORDER BY DISPLAY_POS ASC";
	}
	elsif ($parentItem eq 'CLASSIFICATION') {
		push @queryList,"SELECT CONCAT('C:',ID_CATEGORY),NAME FROM CATEGORY WHERE ID_CLASSIFICATION=$parentID ORDER BY DISPLAY_POS ASC";
	}

	####<Starting HTML
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);

	##<Connecting to the database
	my $dbh=&promsConfig::dbConnect;

	foreach my $query (@queryList) {
		my $sthAna=$dbh->prepare($query);
		$sthAna->execute;
		while (my ($codeID,@anaInfo)=$sthAna->fetchrow_array) {
			for (my $i=0; $i<$#anaInfo-1; $i++) {
				$anaInfo[$i]=&promsMod::shortenName($anaInfo[$i],11);
			}
			print "$codeID#:#",join (' > ',@anaInfo),"\n";
		}
		$sthAna->finish;
	}

	$dbh->disconnect;

	exit;
}

sub ajaxRestrictProteinList {
	my $projectID=param('projectID');
	my $catFilter=param('catFilter') || 0;

	my $dbh=&promsConfig::dbConnect;
	my $sthL=$dbh->prepare("SELECT T.ID_CLASSIFICATION,T.NAME,L.ID_CATEGORY,L.NAME,L.DISPLAY_POS FROM CLASSIFICATION T,CATEGORY L WHERE T.ID_CLASSIFICATION=L.ID_CLASSIFICATION AND T.ID_PROJECT=$projectID");
	$sthL->execute;
	my (%savedLists,%themeInfo);
	while (my ($themeID,$themeName,$listID,$listName,$listPos)=$sthL->fetchrow_array) {
		$themeInfo{$themeID}=$themeName;
		@{$savedLists{$themeID}{$listID}}=($listName,$listPos);
	}
	$sthL->finish;
	$dbh->disconnect;

	###<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<SELECT name="catFilter" class="formSelect" style="width:270px" onchange="checkCatFilterRule(this.value);setComparisonAsModified(2)"><OPTION value="">All proteins</OPTION>
|;
	foreach my $themeID (sort{lc($themeInfo{$a}) cmp lc($themeInfo{$b})} keys %themeInfo) {
		print "<OPTGROUP label=\"$themeInfo{$themeID}\">\n";
		foreach my $listID (sort{lc($savedLists{$themeID}{$a}[1]) cmp lc($savedLists{$themeID}{$b}[1])} keys %{$savedLists{$themeID}}) {
			print "<OPTION value=\"$listID\"";
			print ' selected' if $listID==$catFilter; # in case a list was already selected
			print ">$savedLists{$themeID}{$listID}[0]</OPTION>\n";
		}
		print "</OPTGROUP>\n";
	}
	print "</SELECT>\n";

	exit;
}

sub ajaxUpdateCompGroupAnalyses {
	my $compID=param('compID');

	####<Starting HTML
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);

	##<Connecting to the database
	my $dbh=&promsConfig::dbConnect;

	my $sthAC=$dbh->prepare("SELECT ANA_POS,COMP_GROUP,ID_ANALYSIS FROM ANA_COMPARISON WHERE ID_COMPARISON=$compID");
	my $query1=qq |SELECT E.NAME,G.NAME,SP.NAME,A.NAME
							FROM EXPERIMENT E,GEL2D G,SPOT SP,SAMPLE S,ANALYSIS A
							WHERE E.ID_EXPERIMENT=G.ID_EXPERIMENT AND G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT
							AND S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS=?|;
	my $query2=qq |SELECT E.NAME,S.NAME,A.NAME
							FROM EXPERIMENT E,SAMPLE S,ANALYSIS A
							WHERE E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SAMPLE=A.ID_SAMPLE AND S.ID_SPOT IS NULL AND A.ID_ANALYSIS=?|;
	my @sthList=($dbh->prepare($query1),$dbh->prepare($query2));

	my $sthLC=$dbh->prepare("SELECT CAT_POS,COMP_GROUP,CC.ID_COMPARISON,CL.NAME,CA.NAME FROM CAT_COMPARISON CC,CATEGORY CA,CLASSIFICATION CL WHERE CC.ID_CATEGORY=CA.ID_CATEGORY AND CA.ID_CLASSIFICATION=CL.ID_CLASSIFICATION AND ID_COMPARISON=$compID");

	my %itemList;

	##<Analyses>##
	my %analysisInfo;
	$sthAC->execute;
	while (my ($anaPos,$group,$anaID)=$sthAC->fetchrow_array) {
#foreach my $refGA ([1,2],[1,3],[2,2],[2,1]) {
#my ($group,$anaID)=@{$refGA};
		$group=0 unless $group;
		unless ($analysisInfo{$anaID}) {
			foreach my $sthAna (@sthList) {
				$sthAna->execute($anaID);
				my @anaInfo=$sthAna->fetchrow_array;
				if ($anaInfo[1]) {
					for (my $i=0; $i<=$#anaInfo; $i++) {
						$anaInfo[$i]=&promsMod::shortenName($anaInfo[$i],21);
					}
					$analysisInfo{$anaID}=join(' > ',@anaInfo);
					last;
				}
			}
		}
		$itemList{$group}[--$anaPos]="$group#:#A:$anaID#:#$analysisInfo{$anaID}\n";
	}
	$sthAC->finish;
	$sthList[0]->finish;
	$sthList[1]->finish;

	##<Lists>##
	$sthLC->execute;
	while (my ($listPos,$group,$listID,$themeName,$listName)=$sthLC->fetchrow_array) {
		$themeName=&promsMod::shortenName($themeName,21);
		$listName=&promsMod::shortenName($listName,21);
		$itemList{$group}[--$listPos]="$group#:#C:$listID#:#$themeName > $listName\n";
	}
	$dbh->disconnect;

	foreach my $group (sort{$a<=>$b} %itemList) {
		foreach my $itemData (@{$itemList{$group}}) {
			print $itemData;
		}
	}

	exit;
}

sub ajaxGetProteinDetails {
#print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
	my $projectID=param('id_project');
	my $proteinID=param('protID');
	#my $selVarMod=param('VMOD');
	my $hiddenRule=param('hiddenRule');
	my $peptideRule=param('pepRule');
	my $pepThreshold=param('pepThreshold'); # not used
	my $virtualData=param('virtualData'); # XIC-based identifications
	my $pepSpecificity=param('pepSpecificity');
	my $anaList=param('anaList');
	my @parentItemAna=split(',',param('anaContext'));
	my $bestAnaList=(param('bestAna'))? param('bestAna') : '';
	my $numGroups=param('numGroups');
	my $groupPosStrg=param('groupPos') || ''; # only if parent full_prot called from 1v1_prot

	my (@analysisGroups,@groupPos);
	if ($numGroups || $anaList=~/,/) { # real groups or List comparison
		foreach my $grDataStrg (split(/,/,$anaList)) {
			push @analysisGroups,[split(/\./,$grDataStrg)];
		}
		if ($groupPosStrg) {
			foreach my $pos (split(',',$groupPosStrg)) {
				push @groupPos,$pos;
			}
		}
	}
	else { # fake groups (1 ana / gr)
		foreach my $codeID (split(/\./,$anaList)) {
			push @analysisGroups,[$codeID];
		}
	}
	my @bestAnaList=split(/:/,$bestAnaList);
	my $numAnaGroups=scalar @analysisGroups;

	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;

	###<Project PTMs
	#my ($projectPtmString)=$dbh->selectrow_array("SELECT RELEVANT_PTMS FROM PROJECT WHERE ID_PROJECT=$projectID");
	#my %projectVarMods;
	#if ($projectPtmString) {
	#	foreach my $varMod (split(';',$projectPtmString)) {$projectVarMods{$varMod}=$allPostTransModifs{$varMod};}
	#}
	my $sthGetPM=$dbh->prepare("SELECT ID_MODIFICATION FROM PROJECT_MODIFICATION WHERE ID_PROJECT=$projectID");
	$sthGetPM->execute;
	my %projectVarMods;
	while (my ($modID)=$sthGetPM->fetchrow_array) {
		$projectVarMods{$modID}=$allPostTransModifs{$modID};
	}
	$sthGetPM->finish;

	#<Queries
	my ($alias,$protDes,$organism,$protLength)=$dbh->selectrow_array("SELECT ALIAS,PROT_DES,ORGANISM,PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=$proteinID");
	$protLength='?' unless $protLength;
	my @sthAnaList;
	my $projQuery1=qq |SELECT E.NAME,G.NAME,SP.NAME,A.NAME FROM EXPERIMENT E,GEL2D G,SPOT SP,SAMPLE S,ANALYSIS A
							WHERE A.ID_ANALYSIS=? AND S.ID_SAMPLE=A.ID_SAMPLE AND SP.ID_SPOT=S.ID_SPOT AND G.ID_GEL2D=SP.ID_GEL2D AND E.ID_EXPERIMENT=G.ID_EXPERIMENT|;
	my $projQuery2=qq |SELECT E.NAME,S.NAME,A.NAME FROM EXPERIMENT E,SAMPLE S,ANALYSIS A
							WHERE A.ID_ANALYSIS=? AND S.ID_SAMPLE=A.ID_SAMPLE AND E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SPOT IS NULL|;
	push @sthAnaList,($dbh->prepare($projQuery1),$dbh->prepare($projQuery2)); # free sample then gel
	my $sthCatInfo=$dbh->prepare("SELECT CLASSIFICATION.NAME,CATEGORY.NAME FROM CLASSIFICATION,CATEGORY WHERE CLASSIFICATION.ID_CLASSIFICATION=CATEGORY.ID_CLASSIFICATION AND ID_CATEGORY=?");
	my $sthAP=$dbh->prepare("SELECT VISIBILITY FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=$proteinID AND ID_ANALYSIS=?");
	my $sthCP=$dbh->prepare("SELECT 1 FROM CATEGORY_PROTEIN WHERE ID_PROTEIN=$proteinID AND ID_CATEGORY=?");
	my $pepBegStrg=($virtualData)? '' : 'AND VALID_STATUS > 0'; #'AND PEP_BEG > 0';
	my $isSpecificStrg=($pepSpecificity eq 'unique')? 'AND IS_SPECIFIC=1' : '';
	my $sthPep=$dbh->prepare("SELECT P.ID_PEPTIDE,PEP_SEQ,ABS(PEP_BEG),ABS(PEP_END),VALID_STATUS,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&')
											FROM PEPTIDE P
											LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
											INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE $isSpecificStrg
											WHERE ID_PROTEIN=$proteinID AND P.ID_ANALYSIS=? $pepBegStrg GROUP BY ABS(PEP_BEG),P.ID_PEPTIDE"); # ABS(PEP_BEG), in case multiple occurence in protein

	my (@bestProtVisibility,%peptideList,%peptideOccurence,%peptideByPos,%peptideGroups,%allAnaList);
	foreach my $g (0..$#analysisGroups) {
		my $codeIndex=-1;
		my $contextUsed=0;
		foreach my $codeID (@{$analysisGroups[$g]}) {
			$codeIndex++;
			my ($type,$typeID)=split(':',$codeID); # A or C
			my (@typeInfo,@localAnaList);
			if ($type eq 'A') { # Analysis
				foreach my $sthAna (@sthAnaList) {
					$sthAna->execute($typeID);
					@typeInfo=$sthAna->fetchrow_array;
					last if $typeInfo[0]; # 2 loops in case of project or experiment
				}
				@localAnaList=($typeID);
			}
			else { # C -> List
				$sthCatInfo->execute($typeID); # catID in fact!
				@typeInfo=$sthCatInfo->fetchrow_array;
				@localAnaList=@parentItemAna;
			}
			if ($numGroups) {
				unless ($groupLabels[$g][0]) {$groupLabels[$g][0]=($groupPos[$g])? 'Group #'.$groupPos[$g] : 'Group #'.($g+1);}
				$groupLabels[$g][1].='<BR>' if $groupLabels[$g][1];
				$groupLabels[$g][1].='-'.join(' > ',@typeInfo);
			}
			else {
				@{$groupLabels[$g]}=(&promsMod::shortenName($typeInfo[-1],10),join(' > ',@typeInfo));
			}

			if ($type eq 'C') {
				$sthCP->execute($typeID);
				my ($protInCat)=$sthCP->fetchrow_array;
				next unless $protInCat; # protein is not in List
				next if $contextUsed; # context Analyses already scaned
				$contextUsed=1;
			}
my (%modCode2varModStrg,%varMods);
			foreach my $anaID (@localAnaList) {
				next if ($bestAnaList && $anaID != $bestAnaList[$g]); # no need to find best Ana again
				$sthAP->execute($anaID);
				my ($protVis)=$sthAP->fetchrow_array;
				#next if (!defined($vis) || ($vis==0 && $hiddenRule)); # !defined($vis) because protein might not be in analysis
				next if !defined($protVis); # !defined($vis) because protein might not be in analysis
$allAnaList{$anaID}=1;
				@bestProtVisibility=($protVis,$anaID) if (!$bestProtVisibility[0] || $bestProtVisibility[0] < $protVis);
				my $usedVis=($protVis)? 'VIS' : 'HIDDEN';

				####> PTMs
				$sthPep->execute($anaID);
				my %pep1stPosAna;

				#while (my ($pepID,$pepSeq,$pepBeg,$pepEnd,$valStat)=$sthPep->fetchrow_array) {
				#	my $varModStrg=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$anaID,$pepSeq,undef,$pepBeg-1); # $pepBeg-1: shift mod pos to pos in protein
				#	my $pepEntity="$pepSeq$varModStrg";
				#	unless ($peptideOccurence{$pepEntity}) { # first analysis where pep+vMod found
				#		$peptideList{$pepEntity}{'SEQ'}=$pepSeq; # join('',@pepAA);
				#		$varModStrg=~s/^ \+ //;
				#		$peptideList{$pepEntity}{'VMOD'}=($varModStrg)? $varModStrg : '-'; # %varModData;
				#	}
				#	unless ($pep1stPosAna{$pepID}) { # pep can be found more than once in same protein & analysis (repeated sequence)
				#		if ($peptideRule=~/nr/) {$peptideGroups{$pepEntity}{$g}{$usedVis}=1;}
				#		else {$peptideGroups{$pepEntity}{$g}{$usedVis}++;} # occurence in group
				#		$pep1stPosAna{$pepID}=$pepBeg;
				#		$peptideGroups{$pepEntity}{$g}{'STATUS'}=$valStat if (!$peptideGroups{$pepEntity}{$g}{'STATUS'} || $valStat > $peptideGroups{$pepEntity}{$g}{'STATUS'});
				#	}
				#	$peptideOccurence{$pepEntity}{$anaID}=1;
				#	if (scalar keys %{$peptideOccurence{$pepEntity}}==1) { # 1st ana
				#		$peptideList{$pepEntity}{'POS'}{$pepBeg}=$pepEnd;
				#		$peptideByPos{$pep1stPosAna{$pepID}}{$pepEntity}=1;
				#	}
				#}

				while (my ($pepID,$pepSeq,$pepBeg,$pepEnd,$valStat,$modCode)=$sthPep->fetchrow_array) {
					my $varModStrg;
					if ($modCode) {
						$varModStrg=$modCode2varModStrg{$modCode} || &promsMod::decodeVarMod($dbh,$pepSeq,$modCode,\%varMods,$pepBeg-1);
						$modCode2varModStrg{$modCode}=$varModStrg;
					}
					else {$modCode=$varModStrg='';}
					#my $pepEntity="$pepSeq$varModStrg";
					my $pepEntity="$pepSeq$modCode";
					unless ($peptideOccurence{$pepEntity}) { # first analysis where pep+vMod found
						$peptideList{$pepEntity}{'SEQ'}=$pepSeq; # join('',@pepAA);
						#$varModStrg=~s/^ \+ //;
						$peptideList{$pepEntity}{'VMOD'}=($varModStrg)? $varModStrg : '-'; # %varModData;
					}
					unless ($pep1stPosAna{$pepID}) { # pep can be found more than once in same protein & analysis (repeated sequence)
						if ($peptideRule=~/nr/) {$peptideGroups{$pepEntity}{$g}{$usedVis}=1;}
						else {$peptideGroups{$pepEntity}{$g}{$usedVis}++;} # occurence in group
						$pep1stPosAna{$pepID}=$pepBeg;
$peptideGroups{$pepEntity}{$g}{'STATUS'}=$valStat if (!$peptideGroups{$pepEntity}{$g}{'STATUS'} || $valStat > $peptideGroups{$pepEntity}{$g}{'STATUS'});
					}
					$peptideOccurence{$pepEntity}{$anaID}=1;
					if (scalar keys %{$peptideOccurence{$pepEntity}}==1) { # 1st ana
						$peptideList{$pepEntity}{'POS'}{$pepBeg}=$pepEnd;
						$peptideByPos{$pep1stPosAna{$pepID}}{$pepEntity}=1;
					}
				}

			}
		}
	}
	$sthAnaList[0]->finish; $sthAnaList[1]->finish;
	$sthCatInfo->finish;
	$sthAP->finish;
	$sthCP->finish;
	$sthPep->finish;

	$dbh->disconnect;

	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);
	my ($color1,$color2)=&promsConfig::getRowColors;
	my $protClass=($bestProtVisibility[0])? 'TH' : 'TD';
	my $modWidth=80*$numAnaGroups;
	my $countRulesStrg;
	if ($numGroups) { # grouped Analyses, options: all anr ba banr
		$countRulesStrg=($peptideRule eq 'all')? 'All peptides in Group' : ($peptideRule eq 'anr')? 'Distinct peptides in Group' : ($peptideRule eq 'ba')? 'All peptides in best Analysis' : ($peptideRule eq 'banr')? 'Distinct peptides in best Analysis' : '';
	}
	else { # no groups, options: ba banr
		$countRulesStrg=($peptideRule eq 'ba')? 'All peptides in Analysis' : 'Distinct peptides in Analysis';
	}
	$countRulesStrg.=($hiddenRule==2)? ' only if protein is visible' : ' even if protein is hidden';
	$countRulesStrg.=' [<I>XIC-based extra peptides included</I>]' if $virtualData;
	my $anaStrg=join(',',keys %allAnaList);
	print qq
|<A class="$protClass" href="javascript:sequenceView($proteinID,'$anaStrg')">$alias</A>: $protDes - <FONT class="org">$organism</FONT> ($protLength aa.)<BR>
<B>Count rule:</B> $countRulesStrg.
<TABLE border=0 cellspacing=0>
<TR bgcolor=$color2>
	<TH class="rbBorder">&nbsp;Start&nbsp;</TH><TH class="rbBorder">&nbsp;Peptides&nbsp;</TH><TH class="rbBorder">&nbsp;PTMs&nbsp;</TH>
|;
	foreach my $g (0..$#analysisGroups) {
		my $numAnaStrg=($numGroups)? scalar @{$analysisGroups[$g]}.' items:<BR>' : '';
		my $thClass=($g==$#analysisGroups)? 'bBorder' : 'rbBorder';
		print "\t<TH class=\"$thClass\" nowrap>&nbsp;<FONT class=\"visibleLinkElement\" onmouseover=\"popup('<B>$numAnaStrg$groupLabels[$g][1]</B>')\" onmouseout=\"popout()\">$groupLabels[$g][0]</FONT>&nbsp;</TH>\n";
	}
	print "</TR>\n";
	my $existBadPTMs=0;
	my $bgColor=$color1;
	my $prevPepBegStrg='';
	foreach my $pepBeg (sort{$a<=>$b} keys %peptideByPos) {
		foreach my $pepEntity (sort{$peptideList{$a}{'POS'}{$pepBeg}<=>$peptideList{$b}{'POS'}{$pepBeg} || lc($peptideList{$a}{'VMOD'}) cmp lc($peptideList{$b}{'VMOD'})} keys %{$peptideByPos{$pepBeg}}) {
			if ($hiddenRule==2) { # do not count hidden => Check visible in at least 1 group
				my $okPep=0;
				foreach my $g (0..$#analysisGroups) {
					if ($peptideGroups{$pepEntity}{$g} && $peptideGroups{$pepEntity}{$g}{'VIS'}) {
						$okPep=1;
						last;
					}
				}
				next unless $okPep;
			}
			print "<TR bgcolor=\"$bgColor\" class=\"list\" valign=\"top\">";
			#<Start pos
			my $posStrg='';
			foreach my $pos (sort{$a<=>$b} keys %{$peptideList{$pepEntity}{'POS'}}) { # multiple start pos
				$posStrg.=',' if $posStrg;
				$posStrg.=$pos;
			}
			if ($posStrg ne $prevPepBegStrg) {
				print "<TH align=right>&nbsp;$posStrg&nbsp;</TH>\n";
				$prevPepBegStrg=$posStrg;
			}
			else {print "<TH align=right></TH>\n";}

			#<Sequence
			print "<TD>&nbsp;<FONT class=\"seq\">$peptideList{$pepEntity}{SEQ}</FONT></TD>\n";
			#<PTMs
			#print "<TH align=left>&nbsp;";
			#if (scalar keys %{$peptideList{$pepEntity}{'VMOD'}}) { # $peptideList{$pepEntity}{'VMOD'} always defined
			#	my $vModCount=0;
			#	foreach my $varModCode (sort keys %{$peptideList{$pepEntity}{'VMOD'}}) {
			#		my $varModLabel=($allPostTransModifs{$varModCode})? $allPostTransModifs{$varModCode}[1] : $varModCode;
			#		my $varModClass=($projectVarMods{$varModCode})? $projectVarMods{$varModCode}[2] : 'badPTM';
			#		$existBadPTMs=1 if $varModClass eq 'badPTM';
			#		my $absPosStrg='';
			#		foreach my $absPos (sort{$a<=>$b} keys %{$peptideList{$pepEntity}{'VMOD'}{$varModCode}}) {
			#			$absPosStrg.=',' if $absPosStrg;
			#			$absPosStrg.=($absPos==0)? '(N-ter)' : ($absPos==999999)? '(C-ter)' : $absPos;
			#		}
			#		print '+' if $vModCount;
			#		print "<FONT class=\"$varModClass\">$varModLabel<FONT  class=\"font11\">$absPosStrg</FONT></FONT>";
			#		$vModCount++;
			#	}
			#}
			#else {print '-';}
			#print "&nbsp;</TH>\n";
			print "<TD>&nbsp;$peptideList{$pepEntity}{VMOD}&nbsp;</TD>\n";
			#<Analyses groups
			foreach my $g (0..$#analysisGroups) {
				my $distValueStrg;
				if ($peptideGroups{$pepEntity}{$g}) {
					my ($vTag1,$vTag2)=($peptideGroups{$pepEntity}{$g}{'STATUS'})? ('','') : ('<I>','</I>'); # virtual pep
					if ($peptideGroups{$pepEntity}{$g}{'VIS'}) {
						$distValueStrg="<B>$vTag1$peptideGroups{$pepEntity}{$g}{VIS}$vTag2</B>";
					}
					if ($hiddenRule==2) { # do not count hidden
						$distValueStrg='-' unless $distValueStrg;
					}
					else {
						if ($peptideGroups{$pepEntity}{$g}{'HIDDEN'}) {
							$distValueStrg.='+' if $distValueStrg;
							$distValueStrg.="$vTag1$peptideGroups{$pepEntity}{$g}{HIDDEN}$vTag2";
						}
					}
				}
				else {$distValueStrg='-';}
				print "<TD align=\"center\">$distValueStrg</TD>\n";
			}
			print "</TR>\n";
			$bgColor=($bgColor eq $color1)? $color2 : $color1;
		}
	}
	if ($existBadPTMs) {
		my $colSpan=3+$numGroups;
		print "<TR><TD colspan=$colSpan><I><FONT class=\"badPTM\"><B>PTMs</B> not selected in Project.</FONT></I></TD></TR>\n";
	}
	print qq
|</TABLE>
&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="unselectProtein()">
|;
	exit;
}


sub ajaxGetModificationSiteDetails {
	my $selModifID=param('modID');
	my $site=param('site');
	my $numGroups=param('numGroups');
	my $pepIdStrg=param('pepIDs');
	my $anaList=param('anaList');
	my ($proteinID,$siteRes,$sitePos)=($site=~/^(\d+)-(.)(\d*)/);
my %convertPos2Text=('-'=>'Protein N-term','='=>'Any N-term','+'=>'Protein C-term','*'=>'Any C-term');

#print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;

	##<Items hierarchy
	my @sthAnaList;
	my $projQuery1=qq |SELECT E.NAME,G.NAME,SP.NAME,A.NAME FROM EXPERIMENT E,GEL2D G,SPOT SP,SAMPLE S,ANALYSIS A
							WHERE A.ID_ANALYSIS=? AND S.ID_SAMPLE=A.ID_SAMPLE AND SP.ID_SPOT=S.ID_SPOT AND G.ID_GEL2D=SP.ID_GEL2D AND E.ID_EXPERIMENT=G.ID_EXPERIMENT|;
	my $projQuery2=qq |SELECT E.NAME,S.NAME,A.NAME FROM EXPERIMENT E,SAMPLE S,ANALYSIS A
							WHERE A.ID_ANALYSIS=? AND S.ID_SAMPLE=A.ID_SAMPLE AND E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SPOT IS NULL|;
	push @sthAnaList,($dbh->prepare($projQuery1),$dbh->prepare($projQuery2)); # free sample then gel
	if ($numGroups) { # real groups
		foreach my $grDataStrg (split(',',$anaList)) {
			push @analysisGroups,[split(/\./,$grDataStrg)]; # :
		}
	}
	else { # fake groups (1 ana / gr)
		foreach my $codeID (split(/\./,$anaList)) { # :
			push @analysisGroups,[$codeID];
		}
	}
	my (@groupLabels,%allAnaList);
	foreach my $g (0..$#analysisGroups) {
		foreach my $codeID (@{$analysisGroups[$g]}) {
			my ($type,$anaID)=split(':',$codeID); # Can only be 'A'
			$allAnaList{$anaID}=1;
			my @typeInfo;
			foreach my $sthAna (@sthAnaList) {
				$sthAna->execute($anaID);
				@typeInfo=$sthAna->fetchrow_array;
				last if $typeInfo[0]; # 2 loops in case of project or experiment
			}
			if ($numGroups) {
				unless ($groupLabels[$g][0]) {$groupLabels[$g][0]='Group #'.($g+1);} # groupPos is irrelevent here
				$groupLabels[$g][1].='<BR>' if $groupLabels[$g][1];
				$groupLabels[$g][1].='-'.join(' > ',@typeInfo);
			}
			else {
				@{$groupLabels[$g]}=(&promsMod::shortenName($typeInfo[-1],20),join(' > ',@typeInfo));
			}
		}
	}
	foreach my $sth (@sthAnaList) {$sth->finish;}

	##<Protein info
	my ($alias,$protDes,$organism,$protLength)=$dbh->selectrow_array("SELECT ALIAS,PROT_DES,ORGANISM,PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=$proteinID");
	$protLength='?' unless $protLength;

	##<Modif info
	my ($modifDispCode,$modifACC)=$dbh->selectrow_array("SELECT DISPLAY_CODE,UNIMOD_ACC FROM MODIFICATION WHERE ID_MODIFICATION=$selModifID");
	my $isPhospho=($modifACC==21)? 1 : 0;

	##<Peptide data
	my $sthSite=$dbh->prepare("SELECT P.ID_ANALYSIS,PEP_SEQ,CHARGE,GROUP_CONCAT(DISTINCT ABS(PEP_BEG) ORDER BY ABS(PEP_BEG) ASC),GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',COALESCE(PM.REF_POS_STRING,'') ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),SCORE,DATA,SPEC_COUNT,QUERY_NUM,PEP_RANK
								FROM PEPTIDE_MODIFICATION PM
								INNER JOIN PEPTIDE P ON P.ID_PEPTIDE=PM.ID_PEPTIDE
								INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
								WHERE ID_PROTEIN=$proteinID AND P.ID_PEPTIDE=?
								GROUP BY P.ID_PEPTIDE");
	my $sthAna=$dbh->prepare("SELECT FILE_FORMAT FROM ANALYSIS WHERE ID_ANALYSIS=?");

	my (%anaType,%peptides,%posInSequence,%varModString,%varMods);
	my $g=-1;
	foreach my $grData (split(',',$pepIdStrg)) {
		$g++;
		next unless $grData; # no sites in gr
		foreach my $pepID (split(':',$grData)) {
			$sthSite->execute($pepID);
			my ($anaID,$pepSeq,$charge,$pepBegStrg,$varModCode,$refModInfo,$score,$pepData,$specCount,$qNum,$rank)=$sthSite->fetchrow_array;
			unless ($anaType{$anaID}) {
				$sthAna->execute($anaID);
				($anaType{$anaID})=$sthAna->fetchrow_array;
			}
			my ($pepBeg)=$pepBegStrg=~/^(\d+)/; # in case multiple matches in protein
			unless ($posInSequence{$pepSeq}) {
				$posInSequence{$pepSeq}[0]=$pepBeg;
				$posInSequence{$pepSeq}[1]=$pepBegStrg;
				#$posInSequence{$pepSeq}[2]=($siteRes=~/[-=]/)? 0 : ($siteRes=~/[\+\*]/)? -1 : $sitePos-$pepBeg+1;
			}
			unless ($varModString{$varModCode}) {
				$varModString{$varModCode}=&promsMod::decodeVarMod($dbh,$pepSeq,$varModCode,\%varMods,$pepBeg-1);
			}
			#<Site probability
			my ($prsStrg,$probStrg)=('','');
			my @posProb;
			if ($isPhospho && $pepData && $pepData =~ /PRS=([^##]+)/) { # $projectVarMods{$selModifID}[0]=~/Phospho/ &&
				$prsStrg=$1;
			}
			# MaxQuant proba
			if ($anaType{$anaID} eq 'MAXQUANT.DIR') {
				if ($refModInfo && $refModInfo=~/(^|&)$selModifID:[^&#]*##PRB_MQ=([^&#]+)/) {
					my $modProbStrg=$2;
					foreach my $posProb (split(',',$modProbStrg)) {
						my ($pos,$prob)=split(':',$posProb);
						$pos+=($pepBeg-1) if $pos=~/\d/; # convert to protein pos
						$probStrg=$prob if $pos eq $sitePos;
						push @posProb,[$pos,$prob];
					}
				}
				elsif ($score) {
					$probStrg=1;
					my ($x,$modPosStrg)=$varModCode=~/(^|&)$selModifID:([^&]+)/;
					foreach my $pos (split(/\./,$modPosStrg)) {
						$pos+=($pepBeg-1) if $pos=~/\d/; # convert to protein pos
						push @posProb,[$pos,1]; # set to 100%
					}
				}
			}
			push @{$peptides{$pepSeq}{$varModCode}{$charge}{$g}},[$score,$prsStrg,$probStrg,\@posProb,$specCount,$pepID,$qNum,$rank];
		}
	}
	$sthSite->finish;
	$sthAna->finish;

	$dbh->disconnect;

	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);
	my ($color1,$color2)=&promsConfig::getRowColors;
	my $protClass='TH'; # ($bestProtVisibility[0])? 'TH' : 'TD';

	my $anaStrg=join(',',keys %allAnaList);
	my ($NtermStrg,$siteStrg,$CtermStrg)=('','','');
	if ($siteRes eq '-') {$NtermStrg="<FONT class=\"modifSite mod_$selModifID\">$modifDispCode~</FONT>";}
	elsif ($siteRes eq '+') {$CtermStrg="<FONT class=\"modifSite mod_$selModifID\">~$modifDispCode</FONT>";}
	else {$siteStrg="<FONT class=\"modifSite mod_$selModifID\">-$siteRes$sitePos</FONT>";}
	print qq
|$NtermStrg<A class="$protClass" href="javascript:sequenceView($proteinID,'$anaStrg')">$alias</A>$siteStrg$CtermStrg: $protDes - <FONT class="org">$organism</FONT> ($protLength aa.)<BR>
<TABLE border=0 cellspacing=0>
<TR bgcolor="$color2">
<TH class="rbBorder" rowspan=2>&nbsp;Start&nbsp;</TH><TH class="rbBorder" rowspan=2>&nbsp;Charge&nbsp;</TH><TH class="rbBorder" rowspan=2>&nbsp;Sequence&nbsp;</TH><TH class="rbBorder" rowspan=2>&nbsp;Modifications&nbsp;</TH>
|;
	foreach my $g (0..$#analysisGroups) {
		my $numAnaStrg=($numGroups)? scalar @{$analysisGroups[$g]}.' items:<BR>' : '';
		my $class=($g==$#analysisGroups)? 'bBorder' : 'rbBorder';
		print "<TH class=\"$class\" colspan=3 nowrap>&nbsp;<FONT onmouseover=\"popup('<B>$numAnaStrg$groupLabels[$g][1]</B>')\" onmouseout=\"popout()\">$groupLabels[$g][0]</FONT>&nbsp;</TH>\n";
	}
	print "</TR>\n<TR bgcolor=\"$color2\">";
	foreach my $g (0..$#analysisGroups) {
		my $class=($g==$#analysisGroups)? 'bBorder' : 'rbBorder';
		print "<TH class=\"rbBorder\">&nbsp;Sc.&nbsp;</TH><TH class=\"rbBorder\">&nbsp;Pb.&nbsp;</TH><TH class=\"$class\">&nbsp;Sp.&nbsp;</TH>\n";
	}
	my $numCols=($numGroups*3)+4;
	my $bgColor=$color1;
	my $minPepBeg=0;
	my $prevPepBegStrg='';
	#my $prevCharge=0;
	my $firstVmod=1;
	foreach my $pepSeq (sort{$posInSequence{$a}[0]<=>$posInSequence{$b}[0] || length($a)<=>length($b)} keys %posInSequence) {
		my ($pepBeg,$pepBegStrg)=@{$posInSequence{$pepSeq}};
		my $displayedSeq=$pepSeq;
		if ($sitePos) {
			for (substr $displayedSeq,$sitePos-$pepBeg,1) {$_="</FONT><FONT class=\"modifSite mod_$selModifID\">$_</FONT><FONT>";} # only 1 loop!
			$displayedSeq='<FONT class="seq">'.$displayedSeq.'</FONT>';
		}
		else {
			$displayedSeq=$NtermStrg.'<FONT class="seq">'.$pepSeq.'</FONT>'.$CtermStrg;
		}

		my $aaShiftStrg='';
		if ($minPepBeg) {
			$aaShiftStrg='<FONT class="seq" style="visibility:hidden">'.('X'x($pepBeg-$minPepBeg)).'</FONT>';
		}
		else {$minPepBeg=$posInSequence{$pepSeq}[0];}
		foreach my $varModCode (sort{length($a)<=>length($b)} keys %{$peptides{$pepSeq}}) {
			if ($firstVmod) {$firstVmod=0;}
			else {print "<TR><TD colspan=$numCols></TD></TR>\n";} # empty thin white row
			foreach my $charge (sort{$a<=>$b} keys %{$peptides{$pepSeq}{$varModCode}}) {
				my $maxPepIdx=0;
				foreach my $g (0..$#analysisGroups) {
					next unless $peptides{$pepSeq}{$varModCode}{$charge}{$g};
					my $lastPepIdx=$#{$peptides{$pepSeq}{$varModCode}{$charge}{$g}};
					$maxPepIdx=$lastPepIdx if $maxPepIdx < $lastPepIdx;
				}
				foreach my $pepIdx (0..$maxPepIdx) {
					print "<TR bgcolor=\"$bgColor\" class=\"list\">";
					if ($pepBegStrg ne $prevPepBegStrg) {
						print "<TH align=\"right\">&nbsp;$pepBegStrg&nbsp;</TH>";
						$prevPepBegStrg=$pepBegStrg;
					}
					else {print "<TH></TH>";}
					#if ($charge != $prevCharge) {
						print "<TD align=\"center\">&nbsp;$charge<SUP>+</SUP>&nbsp;</TD>";
						#$prevCharge=$charge;
					#}
					#else {print "<TD></TD>";}
					print "<TD>&nbsp;$aaShiftStrg$displayedSeq&nbsp;</TD><TD nowrap>&nbsp;$varModString{$varModCode}&nbsp;</TD>\n";
					foreach my $g (0..$#analysisGroups) {
						if ($peptides{$pepSeq}{$varModCode}{$charge}{$g} && $peptides{$pepSeq}{$varModCode}{$charge}{$g}[$pepIdx]) {
							my ($score,$prsStrg,$prob,$refPosProb,$specCount,$pepID,$qNum,$rank)=@{$peptides{$pepSeq}{$varModCode}{$charge}{$g}[$pepIdx]};
							if ($score) {
								$score=1*(sprintf "%.2f",$score);
								print "<TD align=\"center\">&nbsp;<A id=\"pep_$pepID\" href=\"javascript:drawSpectrum('pep_$pepID','pep_${pepID}_${qNum}_$rank')\" onmouseover=\"popup('Click to display fragmentation spectrum')\" onmouseout=\"popout()\">$score</A>&nbsp;</TD>";
							}
							else {print "<TD align=\"center\">?</TD>";}
							my $probStrg='?';
							if ($prsStrg) {$probStrg=&phosphoRS::printIcon($prsStrg,{format=>'text',posShift=>$pepBeg-1,pepSeq=>$pepSeq});}
							if ($prob) {
								#$probStrg=&promsMod::MaxQuantProbIcon($prob,{popupText=>'MaxQuant probablity='.($prob*100).'%'});
								my $mqProbStrg='<B>MaxQuant probabilities:</B>';
								foreach my $refPos (@{$refPosProb}) {
									my $pos=$convertPos2Text{$refPos->[0]} || $refPos->[0];
									my $posStrg=($refPos->[0] eq $sitePos)? '<B>'.$pos.'</B>' : $pos;
									$mqProbStrg.="<BR>&nbsp;-$posStrg:".&promsMod::MaxQuantProbIcon($refPos->[1],{text=>'&nbsp;'.($refPos->[1]*100).'%&nbsp;',inPopup=>1});
								}
								$probStrg=&promsMod::MaxQuantProbIcon($prob,{popupText=>$mqProbStrg});
							}
							$specCount='?' unless $specCount;
							print "<TD align=\"center\">&nbsp;$probStrg&nbsp;</TD><TD align=\"right\">&nbsp;$specCount&nbsp;</TD>\n";
						}
						else {
							my $strg=($peptides{$pepSeq}{$varModCode}{$charge}{$g})? '' : '-'; # less peptides in gr VS no peptides in gr
							print "<TD colspan=3 align=\"center\">$strg</TD>\n";
						}

					}
					#print "<TD>&nbsp;$aaShiftStrg$displayedSeq&nbsp;</TD><TD nowrap>&nbsp;$varModString{$varModCode}&nbsp;</TD></TR>\n";
					print "</TR>\n";
					$bgColor=($bgColor eq $color1)? $color2 : $color1;
				}
			}
		}
	}
	print qq
|</TABLE>
&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="unselectSite()">
|;

	exit;
}

####>Revision history<####
# 2.3.1 Minor JS bug fix occuring when no project-relevant PTMs are found (PP 14/11/18)
# 2.3.0 Added 'Frequency' filter (PP 23/10/18)
# 2.2.2 Handles project status=-1 [no auto-end validation] (PP 07/06/18)
# 2.2.1 Minor bug correction for Phospho-sites Venn Diagram (GA 23/05/18)
# 2.2.0 Added peptide specificity filter (PP 18/05/18)
# 2.1.2 Uses &promsMod::cleanDirectory instead of &promsMod::cleanExportDirectory (PP 05/02/18)
# 2.1.1 Warns if unaligned peptide sequences (PP 06/03/17)
# 2.1.0 New modification site comparison option (PP 22/02/17)
# 2.0.8 Floating header removed for better full protein list formatting (PP 20/12/16)
# 2.0.7 Now visibility is ignored when hidden rule is not selected (PP 25/10/16)
# 2.0.6 Added modifications filtering for '1 vs 1 peptide' & venn diagram (PP 18/08/16)
# 2.0.5 Proteins from '1 vs 1 protein' comparison can be listed (PP 23/03/16)
# 2.0.4 Fix minor bug in save comparison (GA 22/03/16)
# 2.0.3 Fix missing ghost peptide filter since v2.0.0 performance improvement (PP 25/02/16)
# 2.0.2 Bug fix in comparison recording (PP 02/02/16)
# 2.0.1 Changed default compare mode to 'full_prot' (PP 01/02/16)
# 2.0.0 Added "1 vs 1 peptide" & "1 vs 1 protein" comparisons & improved performance (PP 03/12/15)
# 1.7.3 Added process progression DIV to prevent server timeout (PP 23/04/15)
# 1.7.2 Ajax-based classification storage & minor bug fix in &getProtDetails (PP 25/08/14)
# 1.7.1 Fix forgotten virtualData param in storeForm (PP 29/04/14)
# 1.7.0 Added 2nd Venn diagram for selected proteins, peptide threshold & list exclusion (PP 25/04/14)
# 1.6.2 Fix UTF-8 issues between DB and Excel (PP 18/04/14)
# 1.6.1 Fix bug cell merge in export when 2 groups with no project relevant PTMS (PP 11/03/14)
# 1.6.0 Allow custom lists comparison & true Excel export & SVG export (PP 14/02/13)
# 1.5.0 Bug fix in peptide+modification query (PP 30/09/13)
# 1.4.9 Record XIC identifications option & optimisation in var mod query (PP 27/09/13)
# 1.4.8 Minor change (PP 24/09/13)
# 1.4.7 Virtual protein management & fix for PTM in protein details popup (PP 25/07/13)
# 1.4.6 Link to compareItems.cgi no longer active: compareItems.cgi is now obsolete (PP 31/05/13)
# 1.4.5 Update with getVariableModifications from promsMod<BR>Modification of ajaxGetProteinDetails (GA 15/05/13)
# 1.4.4 Minor modification following &convertVarModString update in promsMod (GA 17/04/13)
# 1.4.3 Fix for bad peptide count in case "Groups + Distinct in Best Analysis" (PP 01/02/13)
# 1.4.2 Minor change in display (PP 08/01/13)
# 1.4.1 Project status management & user lists renaming (PP 18/09/12)
# 1.4.0 Added Venn diagram & more space for protein name (PP 17/07/12)
# 1.3.9 Change in control panel locking text (PP 07/05/12)
# 1.3.8 Correction in varMods Class (PP 13/07/11)
# 1.3.7 full peptide data for protein details (PP 21/06/11)
# 1.3.6 Autoscrolling header + store proteins bug correction (PP 02/05/11)
# 1.3.5 Add PTM distribution details for selected PTM & protein (PP 28/04/11)
# 1.3.4 bug undefined bestAnalysis (PP 26/04/11)
# 1.3.3 minor bug fix (PP 22/04/11)
# 1.3.2 add PTM display & auto-selection (PP 30/03/11)
# 1.3.1 add Save/Use comparison options (PP 09/03/11)
# 1.3.0 add Classification/Category filter (PP 25/02/11)
# 1.2.9 add absolute delta & export abs/rel delta (PP 25/02/11)
# 1.2.8 adapting to new protein list management options (PP 22/02/11)
# 1.2.7 constant MAX_NUM_GROUPS added (PP 17/12/10)
