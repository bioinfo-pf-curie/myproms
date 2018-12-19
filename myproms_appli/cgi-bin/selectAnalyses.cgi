#!/usr/local/bin/perl -w

################################################################################
# selectAnalyses.cgi      2.6.7                                                #
# Authors: P. Poullet, G. arras, F. Yvon, S. Liva (Institut Curie)             #
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
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;

#print header; warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my %msTypeName=&promsConfig::getMsType;
$msTypeName{'MIS'}="MS/MS"; #redef, to keep space

###############################
####>Recovering parameters<####
###############################
my $projectID=param('id_project'); # defined for duplicate, appFilter & remFilter
my $branchID=param('ID');
my $goType=(param('goType'))? param('goType'):'ana';
my ($item,$itemID)=split(':',$branchID);
my $callType = param('callType') || 'report'; #can be report,list,delete,duplicate,appFilter,remFilter,clearSel,remFlag
$callType='report' if $callType eq 'send'; # just in case
my $srcType=(param('srcType'))? param('srcType') : ($callType eq 'goQuantiAnalysis')? 'design' : 'ana';
my $showParam = param('showParam') || 0;
$item=lc($item);
my $ITEM=uc($item);

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

my $titleString= ($callType eq 'list')? 'List of ' : 'Select ';
my $itemName;
if ($ITEM eq 'GO_ANALYSIS'){
	($itemName)=$dbh->selectrow_array("SELECT NAME FROM $ITEM WHERE ID_GOANALYSIS=$itemID");
}
else {
	($itemName)=$dbh->selectrow_array("SELECT NAME FROM $ITEM WHERE ID_$ITEM=$itemID");
}

###> Option Start Q. GO Analysis, Design or Analysis can be selected
if ($callType eq 'goQuantiAnalysis') {
	my ($isSelectedAna,$isSelectedDes)=($srcType eq 'ana')? (' selected','') : ('',' selected');
	$titleString.="<SELECT name='srcType' class='title' onchange=\"reloadFrame(this.value);\"><OPTION value='design' $isSelectedDes>Design</OPTION><OPTION value='ana' $isSelectedAna>Analysis</OPTION></SELECT> in ".&promsMod::getItemType($ITEM)." <FONT color=#DD0000>$itemName</FONT>";
}
elsif ($ITEM eq 'ANALYSIS') {
	$titleString=($callType eq 'peptidePTM')? "Analysis <FONT color=#DD0000>$itemName</FONT>" : "Selected Analysis <FONT color=#DD0000>$itemName</FONT>";
}
else {
	$titleString.='Analyses in '.&promsMod::getItemType($ITEM)." <FONT color=#DD0000>$itemName</FONT>";
}


########################################################
####>Recovering list of experiment, sample analysis<####
########################################################
my (%listDataBank,@itemAnalyses,@itemDesigns,%okDelete,%okRemFilter,%okActLowScores,%anaProteins,%listParam,@notReportableAna,%modifications,$numSelModifBoxes); #%listQuantif,

####>PhosphoRS<####
my (%prsParam,%prsRunning);

####>Recovering DBs name<####
my $sthAD = $dbh->prepare("SELECT D.ID_DATABANK,NAME FROM ANALYSIS_DATABANK AD,DATABANK D WHERE AD.ID_DATABANK=D.ID_DATABANK AND AD.ID_ANALYSIS=?");


####>Recovering quantitation name<####
#my $sthQuanti = $dbh->prepare("SELECT ID_QUANTIF_METHOD,NAME FROM QUANTIF_METHOD");
#$sthQuanti->execute;
#while (my ($idQuanti,$quantiName)= $sthQuanti->fetchrow_array) {
#	$listQuantif{$idQuanti}=$quantiName;
#}
#$sthQuanti->finish;

####>Generic queries<####
my @sthItem;
my $baseFieldString='ID_ANALYSIS,VALID_STATUS,MS_TYPE,DATA_FILE,FILE_FORMAT,WIFF_FILE,TAXONOMY,MAX_RANK,MIN_SCORE,INSTRUMENT,0,LABELING,LOWER_SCORES'; # 0 ->ref to multiDB !!!WARNING: @projHierarchy added afterwards!!!
if ($ITEM eq 'ANALYSIS') {
	$sthItem[0]=$dbh->prepare("SELECT $baseFieldString,NAME FROM ANALYSIS WHERE ID_ANALYSIS=$itemID ORDER BY DISPLAY_POS ASC");
}
elsif ($ITEM eq 'SAMPLE') {
	$sthItem[0]=$dbh->prepare("SELECT $baseFieldString,NAME FROM ANALYSIS WHERE ID_SAMPLE=$itemID ORDER BY DISPLAY_POS ASC");
}
elsif ($ITEM eq 'SPOT') {
	$sthItem[0]=$dbh->prepare("SELECT $baseFieldString,ANALYSIS.NAME FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_SPOT=$itemID ORDER BY ANALYSIS.DISPLAY_POS ASC");
}
elsif ($ITEM eq 'GEL2D') {
	$sthItem[0]=$dbh->prepare("SELECT $baseFieldString,'SPOT',SPOT.NAME,ANALYSIS.NAME FROM ANALYSIS,SAMPLE,SPOT WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT=SPOT.ID_SPOT AND ID_GEL2D=$itemID ORDER BY SPOT.NAME ASC, ANALYSIS.DISPLAY_POS ASC");
}
elsif ($ITEM eq 'EXPERIMENT') {
	$sthItem[0]=$dbh->prepare("SELECT $baseFieldString,'GEL2D',GEL2D.NAME,'SPOT',SPOT.NAME,ANALYSIS.NAME FROM ANALYSIS,SAMPLE,SPOT,GEL2D WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT=SPOT.ID_SPOT AND SPOT.ID_GEL2D=GEL2D.ID_GEL2D AND GEL2D.ID_EXPERIMENT=$itemID ORDER BY GEL2D.DISPLAY_POS ASC, SPOT.NAME ASC, ANALYSIS.DISPLAY_POS ASC");
	$sthItem[1]=$dbh->prepare("SELECT $baseFieldString,'SAMPLE',SAMPLE.NAME,ANALYSIS.NAME FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT IS NULL AND ID_EXPERIMENT=$itemID ORDER BY SAMPLE.DISPLAY_POS ASC, ANALYSIS.DISPLAY_POS ASC");
	if ($srcType eq 'design') {
		my $sthgetDesign=$dbh->prepare("SELECT ID_DESIGN,NAME,DES FROM DESIGN WHERE ID_EXPERIMENT=$itemID");
		$sthgetDesign->execute;
		while (my (@designData) = $sthgetDesign->fetchrow_array) {
			push @itemDesigns, \@designData;
		}
		$sthgetDesign->finish;
	}
}
elsif ($ITEM eq 'GO_ANALYSIS'){
	$sthItem[0]=$dbh->prepare("SELECT ANALYSIS.$baseFieldString,'GEL2D',GEL2D.NAME,'SPOT',SPOT.NAME,ANALYSIS.NAME FROM ANALYSIS,SAMPLE,SPOT,GEL2D,GOANA_ANALYSIS WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT=SPOT.ID_SPOT AND SPOT.ID_GEL2D=GEL2D.ID_GEL2D AND ANALYSIS.ID_ANALYSIS=GOANA_ANALYSIS.ID_ANALYSIS AND GOANA_ANALYSIS.ID_GOANALYSIS=$itemID ORDER BY GEL2D.DISPLAY_POS ASC, SPOT.NAME ASC, ANALYSIS.DISPLAY_POS ASC");
	$sthItem[1]=$dbh->prepare("SELECT ANALYSIS.$baseFieldString,'SAMPLE',SAMPLE.NAME,ANALYSIS.NAME FROM ANALYSIS,SAMPLE,GOANA_ANALYSIS WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ANALYSIS.ID_ANALYSIS=GOANA_ANALYSIS.ID_ANALYSIS AND ID_GOANALYSIS=$itemID AND SAMPLE.ID_SPOT IS NULL ORDER BY SAMPLE.DISPLAY_POS ASC, ANALYSIS.DISPLAY_POS ASC");
}
elsif ($ITEM eq 'PATHWAY_ANALYSIS') {
	$sthItem[0]=$dbh->prepare("SELECT A.$baseFieldString,'GEL2D',GEL2D.NAME,'SPOT',SPOT.NAME,A.NAME FROM ANALYSIS A,SAMPLE,SPOT,GEL2D,PATHWAYANA_ANALYSIS PA WHERE A.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT=SPOT.ID_SPOT AND SPOT.ID_GEL2D=GEL2D.ID_GEL2D AND A.ID_ANALYSIS=PA.ID_ANALYSIS AND PA.ID_GOANALYSIS=$itemID ORDER BY GEL2D.DISPLAY_POS ASC,SPOT.NAME ASC,A.DISPLAY_POS ASC");
	$sthItem[1]=$dbh->prepare("SELECT A.$baseFieldString,A.NAME FROM ANALYSIS A, PATHWAYANA_ANALYSIS PA WHERE PA.ID_ANALYSIS=A.ID_ANALYSIS AND PA.ID_PATHWAY_ANALYSIS=$itemID");
}
my $sthVVP=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE VISIBILITY>0 AND ID_ANALYSIS=?");
my $sthAVP=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=?");
my $sthSTP=$dbh->prepare ("SELECT COUNT(*) FROM PROTEIN_VALIDATION WHERE SEL_STATUS>=1 AND ID_ANALYSIS=? AND IDENTIFIER NOT LIKE 'DECOY_%'");
my $sthATP=$dbh->prepare ("SELECT COUNT(*) FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=? AND IDENTIFIER NOT LIKE 'DECOY_%'");


####>Call type-specific queries<####
my %sthCTS; # Call Type Specific
if ($callType eq 'remFilter') {
	$sthCTS{RF}=$dbh->prepare("SELECT ID_PROT_VALID FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=? AND SEL_STATUS=-3 LIMIT 0,1");
}
elsif ($callType eq 'report') {
	$sthCTS{LQ}=$dbh->prepare("SELECT 1 FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.FOCUS='protein' AND AQ.ID_ANALYSIS=? LIMIT 0,1"); # labeled Quanti
	$sthCTS{LFQ}=$dbh->prepare("SELECT 1 FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND QUANTIF_ANNOT NOT LIKE '%::SOFTWARE=PD%' AND ID_ANALYSIS=? LIMIT 0,1"); # labeled Quanti
}
elsif ($callType eq 'delete') {
	#>Comparison
	#$sthCTS{'COMP'}=$dbh->prepare("SELECT COUNT(ID_COMPARISON) FROM ANA_COMPARISON WHERE ID_ANALYSIS=?");
	#>Quantification
	#$sthCTS{'LABELING'}=$dbh->prepare("SELECT LABELING FROM ANALYSIS WHERE ID_ANALYSIS=?");
	$sthCTS{'NODESIGN_Q'}=$dbh->prepare("SELECT 1 FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.FOCUS='protein' AND Q.ID_DESIGN IS NULL AND AQ.ID_ANALYSIS=? LIMIT 0,1");
	$sthCTS{'DESIGN_Q'}=$dbh->prepare("SELECT 1 FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_DESIGN IS NOT NULL AND AQ.ID_ANALYSIS=? LIMIT 0,1");
	$sthCTS{'XIC_Q'}=$dbh->prepare("SELECT 1 FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q,QUANTIFICATION_METHOD QM WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND Q.FOCUS='peptide' AND QM.CODE='XIC' AND AQ.ID_ANALYSIS=? LIMIT 0,1");
	#>GO analysis
	$sthCTS{'GO'}=$dbh->prepare("SELECT 1 FROM GOANA_ANALYSIS WHERE ID_ANALYSIS=? LIMIT 0,1");
	$sthCTS{'PATHWAY'}=$dbh->prepare("SELECT 1 FROM PATHWAYANA_ANALYSIS WHERE ID_ANALYSIS=? LIMIT 0,1");
}
elsif ($callType eq 'peptidePTM') {
	$sthCTS{MOD_ID}=$dbh->prepare("SELECT ID_MODIFICATION,SPECIFICITY FROM ANALYSIS_MODIFICATION WHERE ID_ANALYSIS=? AND MODIF_TYPE='V'");
}
#my %sthLS;
#if ($callType eq 'lowScores') {
#	my $allRank=10;
#	foreach my $rank (1..$allRank) {
#		$sthLS{$rank}=$dbh->prepare("SELECT 1 FROM QUERY_VALIDATION WHERE INFO_PEP$rank LIKE '%SEL=-1%' AND ID_ANALYSIS=? LIMIT 0,1");
#	}
#}

foreach my $sthI (@sthItem) {
	$sthI->execute;
	while (my @anaData=$sthI->fetchrow_array) {
		my $anaID=$anaData[0];
		my $labeling=$anaData[11];
		$sthAD->execute($anaID);
		my @dbUsed;
		while (my ($dbID,$dbName)=$sthAD->fetchrow_array) {
			push @dbUsed,$dbID;
			$listDataBank{$dbID}=$dbName;
		}
		$anaData[10]=\@dbUsed;
		push @itemAnalyses,\@anaData;
		if ($callType eq 'delete') {
			$okDelete{$anaID}=1;
			#>Comparison (project-wide) <-- WARNING: Potential empty comparison after deletion!!!
			#$sthCTS{'COMP'}->execute($anaID);
			#my ($numComp)=$sthCTS{'COMP'}->fetchrow_array;
			#if ($numComp) {$okDelete{$anaID}=0;}
			# else { #}
			#>Check for children
			foreach my $child ('XIC_Q','NODESIGN_Q','DESIGN_Q','GO','PATHWAY') {
				$sthCTS{$child}->execute($anaID);
				my ($hasChild)=$sthCTS{$child}->fetchrow_array;
				if ($hasChild) {
					$okDelete{$anaID}=0;
					last;
				}
			}
		}
		if ($anaData[1]>=1) { # valid proteins
			$sthVVP->execute($anaID);
			push @{$anaProteins{$anaID}},$sthVVP->fetchrow_array;
			$sthAVP->execute($anaID);
			push @{$anaProteins{$anaID}},$sthAVP->fetchrow_array;

			#>peptidePTM (PTM distrib)
			if ($callType eq 'peptidePTM') {
				$sthCTS{MOD_ID}->execute($anaID);
				while (my ($modID,$specif)=$sthCTS{MOD_ID}->fetchrow_array) {
					if ($specif) {
						foreach my $res (split(',',$specif)) {
							$modifications{$modID}{RES}{$res}=1;
						}
					}
					else {$modifications{$modID}{RES}{'_'}=1;} # observed cases in BD (why?!!!)
				}
			}
		}
		else {# non-validated proteins
			$sthSTP->execute($anaID);
			push @{$anaProteins{$anaID}},$sthSTP->fetchrow_array;
			$sthATP->execute($anaID);
			push @{$anaProteins{$anaID}},$sthATP->fetchrow_array;
		}
		if ($anaData[1]<=1) {
			if ($callType eq 'remFilter') { # exist Filter?
				$sthCTS{RF}->execute($anaID);
				($okRemFilter{$anaID})=$sthCTS{RF}->fetchrow_array;
			}
			#elsif ($callType eq 'lowScores') { # exist lower Scores?
			#	#foreach my $rank (sort{$a<=>$b} keys %sthLS) {
			#	#	$sthLS{$rank}->execute($anaID);
			#	#	($okActLowScores{$anaID})=$sthLS{$rank}->fetchrow_array;
			#	#	last if $okActLowScores{$anaID};
			#	#}
			#	$okActLowScores{$anaID}=($anaData[12] && $anaData[12]==2)? 0 : 1;
			#}
		}
		%{$listParam{$anaID}}=&promsMod::getSearchParam($dbh,$anaID) if ($showParam==1 || $callType eq 'phosphoRS');

		#>Report: check for quantifications
		if ($callType eq 'report') {
			my ($validStatus,$labeling)=($anaData[1],$anaData[11]);
			if ($validStatus==1) {
				my $sthHasQ=($labeling)? $sthCTS{LQ} : $sthCTS{LFQ};
				$sthHasQ->execute($anaID);
				my ($numQuanti)=$sthHasQ->fetchrow_array;
				push @notReportableAna,$anaID if $numQuanti;
			}
		}

		#>PhosphoRS
		if ($callType eq 'phosphoRS') {
			my $prsParamFile = "$promsPath{valid}/ana_$anaID/PRSparam_ana_$anaID.txt"; # new naming
			unless (-e $prsParamFile) { # check old naming
				my $dataFile = $anaData[3];
				$dataFile =~ s/\.\w{3}$//;
				$dataFile=quotemeta($dataFile);
				$prsParamFile = "$promsPath{valid}/ana_$anaID/PRSparam_$dataFile.txt";
			}
			$prsRunning{$anaID}=0; # always defined
			if (-e $prsParamFile) {
				open PRSparam, $prsParamFile;
				while(<PRSparam>){
					chomp;
					my @line = split /:/;
					$prsParam{$anaID}{$line[0]} = $line[1];
				}
				close PRSparam;
			}
			elsif (-e "$promsPath{tmp}/phosphoRS" && -e "$promsPath{tmp}/phosphoRS/current") {
				my ($errorFile)=(glob "$promsPath{tmp}/phosphoRS/current/$anaID\_*_error.txt")[0];
				if ($errorFile) {$prsRunning{$anaID}=(-s $errorFile)? -1 : 1;} # -1 => error in file
			}
		}

	}
	$sthI->finish;
}
$sthAD->finish;
foreach my $sthI (@sthItem) {$sthI->finish;}
$sthVVP->finish;
$sthAVP->finish;
$sthSTP->finish;
$sthATP->finish;

foreach my $sth (values %sthCTS) {$sth->finish;}


my %quantiData;
my $quantiDataJS = '';
if ($callType eq 'goQuantiAnalysis') {
	my ($ratioMethID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='PROT_RATIO_PEP'");
	my $sthQuanti = ($srcType eq 'ana') ?
			$dbh->prepare("SELECT Q.ID_QUANTIFICATION, Q.NAME, Q.QUANTIF_ANNOT
				      FROM QUANTIFICATION Q, ANA_QUANTIFICATION
				      WHERE Q.ID_QUANTIFICATION=ANA_QUANTIFICATION.ID_QUANTIFICATION
					  AND ID_QUANTIFICATION_METHOD=$ratioMethID
				      AND Q.ID_DESIGN IS NULL
				      AND ID_ANALYSIS=?"):
			$dbh->prepare("SELECT Q.ID_QUANTIFICATION, Q.NAME, Q.QUANTIF_ANNOT
				      FROM QUANTIFICATION Q
				      WHERE ID_QUANTIFICATION_METHOD=$ratioMethID
				      AND ID_DESIGN=?");

	$quantiDataJS .= "var quantiData = {};\n";
	if ($srcType eq 'ana') {
		foreach my $refAnaData (@itemAnalyses){
			my $anaID = $refAnaData->[0];

			# Skip if not SILAC or iTRAQ
			next unless ($refAnaData->[11] && ($refAnaData->[11] =~ /^(SILAC|iTRAQ)/));

			$sthQuanti->execute($anaID);
			while(my ($idQuanti, $quantiName, $annot) = $sthQuanti->fetchrow_array){
				$quantiData{$anaID}{$idQuanti}{'name'} = $quantiName;
				#my ($allRatioString) = ($annot =~ /RATIOS=([^:]+)/);
				my ($allRatioString) = ($annot =~ /::RATIOS=([^:]+)/); # Be careful: MIN_INF_RATIOS comes first in annot string so leave double colons in pattern matching
				next unless $allRatioString;
				my ($states) = ($annot =~ /::STATES=([^:]+)/);
				my @states = split /;/, $states;
				$quantiDataJS .= "quantiData[$idQuanti] = new Array();\n";
				my $rank = 0;
				while($allRatioString =~ /(\d+)\/(\d+)/g){
					$rank++;
					my ($ratioNum, $ratioDenom) = ($1, $2);
					my ($ratioNumName) = ($states[$ratioNum-1] =~ /\d+,([^,]+),\d+/);
					my ($ratioDenomName) = ($states[$ratioDenom-1] =~ /\d+,([^,]+),\d+/);
					$quantiData{$anaID}{$idQuanti}{'ratio'}{$rank} = "$ratioNumName/$ratioDenomName";
					$quantiDataJS .= "quantiData[$idQuanti][$rank-1] = '$ratioNumName/$ratioDenomName';\n";
				}
			}
		}
	}
	else {

		#LABEL=FREE::ALTER=two.sided::BIAS_CORRECTION=TRUE::
		#CONFIDENCE_LEVEL=0.95::
		#FDR_ALPHA=5::
		#FDR_CONTROL=TRUE::
		#FDR_METHOD=FDR-BH::
		#INV_RATIO=FALSE::
		#MIN_INF_RATIOS=0::
		#P_VALUE_OUTLIER=0.05::
		#PEPTIDES=unique;0;0;all;all::
		#QUANTIFICATION_METHOD=Ratio
		#::RATIOS=#117/#116::
		#STATES=6,#127:#2044:#1575:0.#132:#2044:#1576:0.#131:#2044:#1577:0.#130:#2044:#1578:0.#129:#2044:#1579:0.#128:#2044:#1580:0,#116;3,#133:#2045:#1581:0.#135:#2045:#1582:0.#134:#2045:#1583:0,#117
		#::BIAS_COEFF=0.901;0.794;0.602;0.595;0.904;0.901;1.753;1.555;1.759

		foreach my $refAnaData (@itemDesigns){
			my $designID = $refAnaData->[0];
			$sthQuanti->execute($designID);
			while(my ($idQuanti, $quantiName, $annot) = $sthQuanti->fetchrow_array){
				$quantiData{$designID}{$idQuanti}{'name'} = $quantiName;
				my ($allRatioString) = ($annot =~ /::RATIOS=([^:]+)/); # Be careful: MIN_INF_RATIOS comes first in annot string so leave double colons in pattern matching
				next unless $allRatioString;
				my ($states) = ($annot =~ /::STATES=([^:]+)/);
				my @states = split /;/, $states;
				$quantiDataJS .= "quantiData[$idQuanti] = [];\n";
				my $rank = 0;
				#while($allRatioString =~ /#(\d+)\/#(\d+)/g){
				while($allRatioString =~ /#(\d+)(%?\d*)\/#(\d+)(%?\d*)/g){
					$rank++;
					my ($ratioNum, $superNum, $ratioDenom , $superDenom) = ($1, $2, $3, $4);
					my ($ratioNumName)=$dbh->selectrow_array("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=$ratioNum");
					my ($ratioDenomName)=$dbh->selectrow_array("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=$ratioDenom");
					my $normTagNum=($superNum)? "°" :'';
					my $normTagDenom=($superDenom)? "°" :'';
					$quantiData{$designID}{$idQuanti}{'ratio'}{$rank} = "$ratioNumName$normTagNum/$ratioDenomName$normTagDenom";
					$quantiDataJS .= "quantiData[$idQuanti][$rank-1] = '$ratioNumName/$ratioDenomName';\n";
				}
			}
		}
	}
	$sthQuanti->finish;
}
elsif ($callType eq 'peptidePTM') {
	$numSelModifBoxes=(scalar keys %modifications < 4)? scalar keys %modifications : 4;
	my $sthMod=$dbh->prepare("SELECT PSI_MS_NAME,INTERIM_NAME,DES,SYNONYMES FROM MODIFICATION WHERE ID_MODIFICATION=?");
	foreach my $modID (keys %modifications) {
		$sthMod->execute($modID);
		my ($psi,$interim,$des,$syn)=$sthMod->fetchrow_array;
		if ($syn) {
			$syn=~s/^##//;
			$syn=~s/##\Z//;
			$syn=join('/',(split('##',$syn)));
		}
		$modifications{$modID}{NAME}=($psi)? $psi : ($interim)? $interim : ($des)? $des : $syn;
	}
	$sthMod->finish;
}
elsif ($callType eq 'list' && param('quanti')) {
	my $sthGOQuanti = $dbh->prepare("SELECT Q.NAME, GQ.RATIO, Q.QUANTIF_ANNOT
					FROM QUANTIFICATION Q, GOANA_QUANTIFICATION GQ, GO_ANALYSIS G, GOANA_ANALYSIS GA
					WHERE GA.ID_ANALYSIS=?
					AND GA.ID_GOANALYSIS=G.ID_GOANALYSIS
					AND G.ID_GOANALYSIS=GQ.ID_GOANALYSIS
					AND GQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION");
	foreach my $refAnaData (@itemAnalyses){
		my $anaID = $refAnaData->[0];
		$sthGOQuanti->execute($anaID);
		while(my ( $quantiName, $ratio, $annot) = $sthGOQuanti->fetchrow_array){
			my ($allRatioString) = ($annot =~ /RATIOS=([^:]+)/);
			my ($states) = ($annot =~ /STATES=([^:]+)/);
			my @states = split /;/, $states;
			my $rank = 0;
			while($allRatioString =~ /(\d+)\/(\d+)/g){
				$rank++;
				if($rank == $ratio){
					my ($ratioNum, $ratioDenom) = ($1, $2);
					my ($ratioNumName) = ($states[$ratioNum-1] =~ /\d+,([^,]+),\d+/);
					my ($ratioDenomName) = ($states[$ratioDenom-1] =~ /\d+,([^,]+),\d+/);
					$quantiData{$anaID}{'Ratio'} = "$ratioNumName/$ratioDenomName";
					last;
				}
			}
			$quantiData{$anaID}{'Quantification'} = $quantiName;
		}
	}
	$sthGOQuanti->finish;
}

$dbh->disconnect;

#######################
####>Starting HTML<####
#######################
my ($lightColor,$darkColor)=&promsConfig::getRowColors;

print header(-'content-encoding'=>'no',-charset=>'UTF-8');
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Select Analysis</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TD.center {text-align:center} //{font-weight:bold;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
$quantiDataJS
function cancelAction() {
	//top.promsFrame.selectedAction='summary'; // set default action to 'summary'
	top.promsFrame.optionFrame.selectOption();
}
|;
if ($callType eq 'list') {
	my $addParam = '';
	if(param('quanti')){
		$addParam = "+'&quanti=1'";
	}
	print qq
|function changeViewParam(expMode) {
	window.location='./selectAnalyses.cgi?callType=list&ID=$branchID&showParam='+expMode$addParam;
}
|;
}
else {
	if ($callType eq 'report') {
		print "var notReportableAna=[",join(',',@notReportableAna),"];\n";
		print qq
|function checkAnalysesForQuantif(reportCode) {
	if (reportCode==0 \|\| reportCode==2) { // report or report + end validation
		var okAlert=false;
		for (var i=0; i < notReportableAna.length; i++) {
			var chkbox=document.getElementById('chk_'+notReportableAna[i]);
			if (!okAlert && chkbox.checked) {
				alert('WARNING: Some Analyses are associated with quantification data.\\nThese data must be deleted before new Reporting.');
				okAlert=true;
			}
			chkbox.disabled=true;
		}
	}
	else { // reactivate all disabled checkboxes
		var anaBox=document.selAnaForm.anaList;
		if (!anaBox) return;
		if (anaBox.length) { // more than 1 checkboxes
			for (i=0; i < anaBox.length; i++){
				anaBox[i].disabled=false;
			}
		}
		else {anaBox.disabled=false;} // 1 checkbox only
	}
}
|;
	}
	print qq
|function testCheckbox() {
	var anaBox=document.selAnaForm.anaList;
	if (!anaBox) return 0; // no selectable analyses
	var selected=0;
	if (anaBox.length) { // more than 1 checkboxes
		for (i=0; i < anaBox.length; i++){
			if (anaBox[i].checked==true && anaBox[i].disabled==false) {selected=1; break;}
		}
	}
	else if (anaBox.checked==true && anaBox.disabled==false){selected=1;}
	return selected;
}
function checkall(checkStatus){
	var anaBox=document.selAnaForm.anaList;
	if (!anaBox) return; // no selectable analyses
	if (anaBox.length) { // more than 1 checkboxes
		for (i=0; i < anaBox.length; i++){
			if(anaBox[i].disabled == false){
				anaBox[i].checked=checkStatus;
			}
		}
	}
	else {anaBox.checked=checkStatus;} // Only 1 checkbox
}
|;
	if ($callType eq 'phosphoRS') {
		print qq
|function checkForm(myForm) {
	if (testCheckbox()==0) {alert('No Analyses selected !'); return false;}

	if(myForm.probThreshold.value > 100 \|\| myForm.probThreshold.value <= 0){
		alert('Enter a valid threshold');
		return false;
	}

	if (myForm.overwrite.checked==true) {
		if(!confirm("PhosphoRS will overwrite all changes made for position in manual validation. Are you sure?")){
			return false;
		}
	}
	return true;
}
|;
	}
	elsif ($callType eq 'goQuantiAnalysis') {
		my ($listName,$itemName)=($srcType eq 'ana')? ('anaList','Analysis') : ('desList','Design');
		print qq|
function reloadFrame(srcType){
	window.location='./selectAnalyses.cgi?ID=$branchID&id_project=$projectID&callType=$callType&srcType='+srcType;
}
function loadRatios(anaOrDesID,quantiID){
	var selectRatio = document.getElementById('ratio'+anaOrDesID);

	selectRatio.options.length = 0;

	if(quantiData[quantiID].length > 1){
		selectRatio.add(new Option('-= Select =-', 0));
	}
	for(var i=0;i<quantiData[quantiID].length;i++){
		selectRatio.add(new Option(quantiData[quantiID][i], i+1));
	}
	selectRatio.selectedIndex = 0;
}
function checkForm(myForm) {
	if (document.selAnaForm.$listName.length) {
	    for(var i=0;i<document.selAnaForm.$listName.length;i++){
			if(document.selAnaForm.${listName}[i].checked){
		        selectQuantiAndRatio(document.selAnaForm.${listName}[i].value);
		        return true;
			}
	    }
	}
	else {
	    if (document.selAnaForm.$listName.checked) {
		selectQuantiAndRatio(document.selAnaForm.$listName.value);
		return true;
	    }
	}
	alert('ERROR: Select 1 $itemName');
	return false;
}
function selectQuantiAndRatio(id) {
        document.getElementById('quanti').value = document.getElementById('quanti'+id).value;
        document.getElementById('ratio').value = document.getElementById('ratio'+id).value;
}
|;
	}
	elsif ($callType eq 'peptidePTM') {
		print qq
|function updateModifFiltering(type) {
	var modFilter=document.getElementById('modFilter');
	if (type) {
		modFilter.disabled=false;
	}
	else {
		modFilter.disabled=true;
		modFilter.options.selectedIndex=0;
		updateModifSelectionList();
	}
	document.getElementById('siteCountSPAN').style.display=(type=='restrict')? '' : 'none';
}
var prevModID=0;
function updateModifSelectionList(modID) {
	if (prevModID) { // enable previously disabled modif
		for (var s=1; s<=$numSelModifBoxes; s++) {
			var curSelBox=document.getElementById('modifID_'+s);
			for (var i=0; i< curSelBox.options.length; i++) {
				if (curSelBox.options[i].value==prevModID) {
					curSelBox.options[i].disabled=false;
					break;
				}
			}
		}
	}
	if (modID) {
		for (var s=1; s<=$numSelModifBoxes; s++) {
			var curSelBox=document.getElementById('modifID_'+s);
			for (var i=0; i< curSelBox.options.length; i++) {
				if (curSelBox.options[i].value==modID) {
					curSelBox.options[i].disabled=true;
					if (curSelBox.options[i].selected) {
						curSelBox.options[i].selected=false;
						delete selModifInBox[s][modID];
					}
					break;
				}
			}
		}
	}
	prevModID=modID;
}
var selModifInBox={},numSelInBox={};
for (var i=1; i<=$numSelModifBoxes; i++) {selModifInBox[i]={};numSelInBox[i]=0;}
function updateModifSelection(setNumber) {
	var curSelBox=document.getElementById('modifID_'+setNumber);
	//Fetching selected option
	var curSel={},numSel=0;
	for (var i=0; i< curSelBox.options.length; i++) {
		if (curSelBox.options[i].selected) {
			curSel[curSelBox.options[i].value]=i; // record each value and its index
			numSel++;
		}
	}
	//Comparing with previous list
	if (numSel == numSelInBox[setNumber]) { // switch in option selected (possible if only 1 option is selected)
		// remove previous selected option
/* Commented to allow selection of same modif in different sets (PP 17/03/16)
*		for (var selVal in selModifInBox[setNumber]) { // only 1 selVal
*			updateOtherBoxes(setNumber,selModifInBox[setNumber][selVal],false);
*		}
*/
		selModifInBox[setNumber]={};
		// add new option
		for (var selVal in curSel) {
			selModifInBox[setNumber][selVal]=curSel[selVal];
//			updateOtherBoxes(setNumber,curSel[selVal],true); // disable option in other selects
		}
	}
	else if (numSel > numSelInBox[setNumber]) { // new option selected
		for (var selVal in curSel) {
			if (selModifInBox[setNumber][selVal] != undefined) continue; // can be 0
			selModifInBox[setNumber][selVal]=curSel[selVal];
			numSelInBox[setNumber]++;
//			updateOtherBoxes(setNumber,curSel[selVal],true); // disable option in other selects
			break;
		}
	}
	else { // option deselected
		for (var selVal in selModifInBox[setNumber]) {
			if (curSel[selVal] != undefined) continue; // can be 0
			numSelInBox[setNumber]--;
//			updateOtherBoxes(setNumber,selModifInBox[setNumber][selVal],false); // activate option in other selects
			delete selModifInBox[setNumber][selVal];
			break;
		}
	}
	//Updating target residues
	var residueSet=document.selAnaForm['residues_'+setNumber];
	if (residueSet.length) {
		if (numSelInBox[setNumber]) {
			for (var i=0; i<residueSet.length; i++) {
				var modIDs=residueSet[i].value.split(':'); // 1st index is residue code! (required for submit)
				var matchedMod=false;
				for (var j=1; j<modIDs.length; j++) {
					if (selModifInBox[setNumber][modIDs[j]] != undefined) { // can be 0
						matchedMod=true;
						break;
					}
				}
				residueSet[i].checked=(matchedMod)? true : false;
				residueSet[i].disabled=(matchedMod)? false : true;
				document.getElementById('residuesDIV_'+setNumber+':'+i).style.display=(matchedMod)? '' : 'none';
			}
		}
		else {
			for (var i=0; i<residueSet.length; i++) {
				residueSet[i].checked=false;
				residueSet[i].disabled=true;
				document.getElementById('residuesDIV_'+setNumber+':'+i).style.display='none';
			}
		}
	}
	else {
		if (numSelInBox[setNumber]) {
			residueSet.checked=true;
			residueSet.disabled=false;
			document.getElementById('residuesDIV_'+setNumber+':0').style.display='';
		}
		else {
			residueSet.checked=false;
			residueSet.disabled=true;
			document.getElementById('residuesDIV_'+setNumber+':0').style.display='none';
		}
	}
}
/* Commented to allow selection of same modif in different sets (PP 17/03/16)
*function updateOtherBoxes(setNumber,valIdx,disabStatus) {
*	for (var i=1; i<=$numSelModifBoxes; i++) {
*		if (i==setNumber) continue;
*		var curBox=document.getElementById('modifID_'+i).options[valIdx].disabled=disabStatus;
*	}
*}
*/
function showHideTargetResidues(chkStatus) {
	document.getElementById('setInterSPAN').style.display=(chkStatus)? '' : 'none';
	for (var i=1; i<=$numSelModifBoxes; i++) {
		document.getElementById('modifResDIV_'+i).style.display=(chkStatus)? '' : 'none';
	}
}
function checkForm(myForm) {
	if (myForm.modFilterType.value && !myForm.modFilter.value) {
		alert('Modification filter not selected!');
		return false;
	}
	var okModifs=false;
	for (var i=1; i<=$numSelModifBoxes; i++) {
		if (document.getElementById('modifID_'+i).value) {
			okModifs=true;
			break;
		}
	}
	if (!okModifs) {alert('ERROR: No modification selected!'); return false;}
	if (testCheckbox()==0) {alert('ERROR: No Analyses selected!'); return false;}
	return true;
}
|;
	}
	else {
		print qq
|function checkForm(myForm) {
	if (testCheckbox()==0) {alert('No Analyses selected!'); return false;}
	if (myForm.endVal && myForm.endVal.value==-1) {alert('ERROR: No action selected!'); return false;} // endVal only if callType=report
	if (myForm.clearMode && myForm.clearMode.value==-1) {alert('ERROR: No action selected!'); return false;} // clearMode only if callType=clearSel
	if (myForm.templateFile && !myForm.templateFile.value) {alert('ERROR: No template file provided!'); return false;} // only if callType=impElution
	return true;
}
|;
	}
}
print qq
|</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" >
<CENTER>
<FONT class="title">$titleString</FONT><BR><BR>
|;
if ($callType eq 'goQuantiAnalysis' && $srcType eq 'design' && $#itemDesigns < 0) {
	print "<BR><BR>No design available (switch select to option <B>Analyses</B>).<BR>\n";
	print "</BODY>\n</HTML>\n";
	exit;
}

my $formAction=($callType eq 'impElution')? 'importElution.cgi' :
	($callType eq 'report')? 'send2Biologist.cgi' :
	($callType eq 'delete')? 'deleteProjectItem.cgi' :
	($callType eq 'lowScores')? 'autoSelect.cgi' :
	($callType=~/Filter/ || $callType eq 'duplicate')? 'filterAnalysis.cgi' :
	($callType eq 'clearSel')? 'autoSelect.cgi' :
	($callType eq 'remFlag')? 'filterPeptideAnalysis.cgi' :
	($callType eq 'phosphoRS')? 'analysePhospho.cgi' :
	($callType eq 'goQuantiAnalysis')? 'startGOQuantiAnalysis.cgi' :
	($callType eq 'peptidePTM')? 'showPepModifDistribution.cgi' :
	'selectAnalysis.cgi';
my $hrColspan;
if ($callType eq 'list') {
	$hrColspan=1;
	my ($buttonAction,$expMode)=($showParam)? ('Hide',0) : ('Show',1);
	print qq
|<TABLE border=0 cellspacing=0 cellpadding=0>
<TR><TD colspan=8><INPUT type="button" value="$buttonAction Search Parameters" onclick="changeViewParam($expMode)"></TD></TR>
<TR bgcolor="$darkColor">
|;
}
else {
	$hrColspan=2;
	print "<FORM name=\"selAnaForm\" action=\"./$formAction\" method=\"post\" onsubmit=\"return (checkForm(this));\" enctype=\"multipart/form-data\">\n";
	if ($callType eq 'delete') {print "<INPUT type=\"hidden\" name=\"form_ID\" value=\"$branchID\">\n";}
	elsif ($callType=~/Filter/) {
		print "<INPUT type=\"hidden\" name=\"id_project\" value=\"$projectID\">\n";
		if ($callType eq 'remFilter') {print "<INPUT type=\"hidden\" name=\"remFilter\" value=\"1\">\n";}
	}
	elsif ($callType eq 'impElution') {
		print qq
|<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="save" value="1"/>
|;
	}
	elsif ($callType eq 'duplicate') {
		print qq
|<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="duplicate" value="1">
<INPUT type="hidden" name="apply" value="1">
|;
	}
	elsif ($callType eq 'lowScores') {
		print qq
|<INPUT type="hidden" name="ACT" value="$callType">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="ITEM" value="$ITEM">
<INPUT type="hidden" name="ID" value="$itemID">
|;
	}
	elsif($callType eq 'remFlag'){
		print qq
|<INPUT type="hidden" name="apply" value="Remove Flags">
|;
	}
	elsif($callType eq 'phosphoRS'){
		print qq
|	<TABLE border=0 align="center" cellspacing=0 cellpadding=5>
<TR bgcolor="$darkColor">
	<TH colspan=2><FONT class="title2">PhosphoRS Analysis Rules</FONT></TH>
</TR>
<TR bgcolor="$lightColor">
	<TH width="10px"></TH>
	<TD>Position modification threshold:&nbsp;<INPUT type="text" value="50" name="probThreshold" size=3>%<BR>
	Activation type: <SELECT name="activationType">
					<OPTION selected>CID</OPTION>
					<OPTION>ETD</OPTION>
					<OPTION>HCD</OPTION>
				</SELECT><BR>
	Mass Deviation: <INPUT type="text" value="0.5" name="massTolerance" size=3>Da<BR>
	<INPUT type="checkbox" name="overwrite" id="overwrite" value="1"> Check to overwrite manual phosphorylation validation of positions
	</TD>
</TR>
</TABLE>
<BR>
|;
	}
	elsif($callType eq 'goQuantiAnalysis') {
		print qq
|<INPUT type="hidden" name="quanti" id="quanti" value="">
<INPUT type="hidden" name="ratio" id="ratio" value="">
|;
	}
	elsif ($callType eq 'peptidePTM') {
		print qq
|<INPUT type="hidden" name="ID" value="$branchID">
<INPUT type="hidden" name="numSelBoxes" value="$numSelModifBoxes"/>
<TABLE bgcolor="$darkColor">
<TR><TH align="right" valign="top">&nbsp;Peptide selection :</TH><TH bgcolor="$lightColor" align="left">
&nbsp;<INPUT type="checkbox" name="virtualPep" value="1" checked> Include peptides identified from XIC extractions (if any)&nbsp;<BR>
&nbsp;<INPUT type="checkbox" name="quantifPep" value="1"> Restrict to peptides with quantification data&nbsp;<BR>
&nbsp;<INPUT type="checkbox" name="residuePep" value="1" onclick="showHideTargetResidues(this.checked)"> Restrict to modifiable peptides (select target residues)&nbsp;<BR>
<SPAN id="setInterSPAN" style="display:none">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="setInter" value="1"> Combine required target residues across sets&nbsp;<BR></SPAN>
&nbsp;&nbsp;Modification filter: <SELECT name="modFilterType" onchange="updateModifFiltering(this.value)"><OPTION value="">-= No filtering =-</OPTION><OPTION value="restrict">Restrict to</OPTION><OPTION value="exclude">Exclude</OPTION></SELECT>
<SELECT name="modFilter" id="modFilter" onchange="updateModifSelectionList(this.value)" disabled><OPTION value="">-= None =-</OPTION>
|;
		foreach my $modID (sort{lc($modifications{$a}{NAME}) cmp lc($modifications{$b}{NAME})} keys %modifications) {
			print "<OPTION value=\"$modID\">$modifications{$modID}{NAME}</OPTION>\n";
		}
		print qq
|</SELECT><BR>
<SPAN id="siteCountSPAN" style="display:none">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="sitesFromRestrict" value="1"> Count  sites from this modification&nbsp;<BR></SPAN>
</TH></TR>
<TR><TH align="right" valign="top">&nbsp;Modification selection :</TH><TD bgcolor="$lightColor" valign="top"><TABLE cellpadding=0></TR>
|;
		my %convertResidue=('-'=>'Protein N-terminal','='=>'N-terminal','+'=>'Protein C-terminal','*'=>'C-terminal','_'=>'Unspecified');
		my %residue2Modifs;
		foreach my $modID (keys %modifications) {
			foreach my $res (keys %{$modifications{$modID}{RES}}) {
				$residue2Modifs{$res}{$modID}=1;
			}
		}
		my $selBoxHeight=(scalar keys %modifications <= 2)? '40px' : '75px';
		foreach my $i (1..$numSelModifBoxes) {
			print "<TH bgcolor=\"$lightColor\" nowrap valign=\"top\">&nbsp;Set #$i<BR><SELECT multiple name=\"modifID_$i\" id=\"modifID_$i\" style=\"height:$selBoxHeight; font-weight:bold\" onchange=\"updateModifSelection($i)\">\n";
			foreach my $modID (sort{lc($modifications{$a}{NAME}) cmp lc($modifications{$b}{NAME})} keys %modifications) {
				print "<OPTION value=\"$modID\">$modifications{$modID}{NAME}</OPTION>\n";
			}
			print "</SELECT><DIV id=\"modifResDIV_$i\" style=\"text-align:left;display:none\">";
			my $idx=0;
			foreach my $res (sort keys %residue2Modifs) {
				my $dispRes=($convertResidue{$res})? $convertResidue{$res} : $res;
				print "<DIV id=\"residuesDIV_$i:$idx\" style=\"display:none\">&nbsp;&nbsp;<INPUT type=\"checkbox\" name=\"residues_$i\" value=\"$res:",join(':',keys %{$residue2Modifs{$res}}),"\" disabled>$dispRes</DIV>\n";
				$idx++;
			}
			print "</DIV></TH>\n";
		}
		print qq
|</TR></TABLE></TD></TR>
<TR><TH align="right" valign="top">&nbsp;Protein distribution :</TH><TH bgcolor="$lightColor" align="left">&nbsp;<INPUT type="checkbox" name="protMulti" value="1"> A protein can be "with" and "without" selected modification(s) at the same time&nbsp;</TH></TR>
</TABLE><BR>
|;
	}
	print qq
|<TABLE border=0 cellspacing=0 cellpadding=0>
<TR bgcolor="$darkColor">
	<TH class="rbBorder">
|;
	unless($callType eq 'goQuantiAnalysis') {
		print qq|<INPUT type="checkbox" onclick="checkall(this.checked)">|;
	}
	else {print "&nbsp;";}
}
if ($srcType eq 'ana') {
	print qq
|	</TH>
	<TH class="rbBorder" nowrap colspan=2>&nbsp;Analysis name&nbsp;</TH>
	<TH class="rbBorder">&nbsp;MS type&nbsp;<BR>& nile</TH>
	<TH class="rbBorder">&nbsp;Instrument&nbsp;</TH>
	<TH class="rbBorder">Search file<BR>&nbsp;& Engine&nbsp;</TH>
	<TH class="rbBorder">Databank(s)<BR>&nbsp;Taxonomy&nbsp;</TH>
	<TH class="rbBorder" nowrap>&nbsp;Min. score&nbsp;<BR>&nbsp;Max. rank&nbsp;</TH>
	<TH class="rbBorder">&nbsp;Labeling&nbsp;</TH>
|;
}
else {
	print qq
|	</TH>
	<TH class="rbBorder" nowrap colspan=2>Design Name</TH>
|;
}
if ($callType eq 'phosphoRS') {
	print qq
|	<TH class="rbBorder">&nbsp;Selected&nbsp;<BR>&nbsp;proteins&nbsp;</TH>
	<TH class="bBorder">&nbsp;PhosphoRS&nbsp;<BR>&nbsp;status&nbsp;</TH>
|;
}
elsif ($callType eq 'goQuantiAnalysis') {
	print "<TH class=\"rbBorder\">&nbsp;Selected&nbsp;<BR>&nbsp;proteins&nbsp;</TH>\n" unless $srcType eq 'design';
	print qq
|	<TH class="rbBorder">&nbsp;Select Quantification&nbsp;</TH>
	<TH class="rbBorder">&nbsp;Ratio&nbsp;</TH>
|;
}
elsif($callType eq 'list' and param('quanti')){
	print qq
|	<TH class="rbBorder">&nbsp;Selected&nbsp;<BR>&nbsp;proteins&nbsp;</TH>
	<TH class="rbBorder">&nbsp;Selected Quantification&nbsp;</TH>
	<TH class="rbBorder">&nbsp;Selected&nbsp;<BR>&nbsp;Ratio&nbsp;</TH>
|;
}
else {
	print "<TH class=\"bBorder\">&nbsp;Selected&nbsp;<BR>&nbsp;proteins&nbsp;</TH>";
}
my %itemIcones=&promsConfig::getItemIcones;
#my $bgColor=($ITEM eq 'SAMPLE' || $ITEM eq 'SPOT')? $lightColor : $darkColor;
my $bgColor=$darkColor;
my %prevItemName;
my $disabSubmit=' disabled';
if ($srcType eq 'ana') {
	foreach my $refAnaData (@itemAnalyses) {
		my ($anaID,$valStat,$msType,$dataFile,$fileFormat,$wiffFile,$taxonomy,$maxRank,$minScore,$instrument,$refDbUsed,$labeling,$lowerScores,@projHierarchy)=@{$refAnaData};
		if ($taxonomy) {
			$taxonomy=~s/\(.*\)//;
			$taxonomy='('.$taxonomy.')';
		}
		else {$taxonomy='';}
		$msType=$msTypeName{$msType}; # global
		$fileFormat=~s/\..*//;
		$instrument='-' unless $instrument;
		$maxRank='-' unless $maxRank;
		$minScore='-' unless defined $minScore;
		#my $quantifMethod=($quantifID)? $listQuantif{$quantifID} : 'None';
		$labeling='None' unless $labeling;
		$labeling=~s/ .+//; # iTRAQ 4plex -> iTRAQ
		##>Row color
		my $fatherIt=$projHierarchy[-3];
		if (($fatherIt && (!$prevItemName{$fatherIt} || $prevItemName{$fatherIt} ne $projHierarchy[-2])) || $ITEM =~ /SAMPLE|SPOT/) { # keep color if same analysis parent item (SAMPLE or SPOT)
			$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
		}
		elsif ($ITEM eq 'EXPERIMENT' || $ITEM eq 'GEL2D') {
			my $colspan=($callType eq 'goQuantiAnalysis')? 10 :($callType eq 'phosphoRS')? 9: 8;
			print "<TR bgcolor=$bgColor><TD colspan=$hrColspan></TD><TD colspan=$colspan><HR width=98%></TD></TR>\n";
		}
		print "<TR valign=middle bgcolor=$bgColor>\n";
		##>Checkbox
		if ($callType ne 'list') {
			my $disabStrg = '';
			if($callType eq 'phosphoRS'){
				$disabStrg = ' disabled' unless ($listParam{$anaID}{'g:Variable modifications'} && $listParam{$anaID}{'g:Variable modifications'} =~ /Phospho/ && !$prsRunning{$anaID});
			}
			elsif ($callType eq 'goQuantiAnalysis'){
				$disabStrg = ' disabled' unless ($quantiData{$anaID});
			}

			#my $boxStr = ($callType ne 'delete' && ($valStat>1 || ($callType eq 'remFilter' && !$okRemFilter{$anaID}) || ($callType eq 'lowScores' && !$okActLowScores{$anaID}) || ($callType eq 'report' && $valStat==-1) || ($callType eq 'impElution' && $instrument ne 'MALDI-TOF-TOF')))? '-' : "<INPUTtype=\"checkbox\" name=\"anaList\" value=\"$anaID\">";
			my $boxType = ($callType eq 'goQuantiAnalysis')? 'radio' : 'checkbox';
			my $checkStr=($ITEM eq 'ANALYSIS')? 'checked' : '';
			my $boxStr = (($callType eq 'delete' && !$okDelete{$anaID}) ||
							($callType eq 'peptidePTM' && $valStat <= 0) ||
							($callType ne 'delete' && $callType ne 'peptidePTM' &&
								($valStat>1 || ($callType eq 'remFilter' && !$okRemFilter{$anaID})
									#|| ($callType eq 'lowScores' && !$okActLowScores{$anaID})
									|| ($callType eq 'lowScores' && $lowerScores && $lowerScores==2)
									|| ($callType eq 'report' && $valStat==-1)
									|| ($callType eq 'impElution' && $instrument ne 'MALDI-TOF-TOF')
								)
							)
						)? '-' : "<INPUT type=\"$boxType\" name=\"anaList\" id=\"chk_$anaID\" value=\"$anaID\"$disabStrg$checkStr>";
			print "\t<TH valign=middle>$boxStr</TH>\n";
			$disabSubmit='' if $boxStr ne '-' and $callType ne 'phosphoRS'; # at least 1 selectable analysis => activate Submit
		}
		##>Parents
		my $parentStrg='';
		for (my $i=0;$i<=$#projHierarchy-2;$i+=2) { # stops before ana name
			my $IT=$projHierarchy[$i];
			my $itName=$projHierarchy[$i+1];
			if ($prevItemName{$IT} && $prevItemName{$IT} eq $itName) {
				$parentStrg.="<FONT color=\"$bgColor\">$itName&nbsp;&nbsp;&gt;&nbsp;&nbsp;</FONT>";
			}
			else {
				$parentStrg.="$itName&nbsp;&nbsp;&gt;&nbsp;&nbsp;";
				$prevItemName{$projHierarchy[$i]}=$itName;
				for (my $j=$i+2;$j<$#projHierarchy-1;$j+=2) {$prevItemName{$projHierarchy[$j]}='';}
			}
		}
		##>Analysis
		my $anaCode=($valStat==-1)? 'analysis:no_scan' : ($valStat==0)? 'analysis:no_val' : ($valStat==1)? 'analysis:part_val' : 'analysis:val';
		my $dbStrg;
		foreach my $dbID (@{$refDbUsed}) {
			$dbStrg.='&bull;' if $dbStrg;
			$dbStrg.=$listDataBank{$dbID};
		}
		print qq
|	<TH nowrap align=left valign=middle>&nbsp;&nbsp;$parentStrg</TH>
	<TH nowrap align=left valign=middle><IMG src="$promsPath{images}/$itemIcones{$anaCode}">&nbsp;$projHierarchy[-1]&nbsp;</TH>
	<TD class="center">$msType<BR>&nbsp;$wiffFile&nbsp;</TD>
	<TD class="center">&nbsp;$instrument&nbsp;</TD>
	<TD class="center">&nbsp;$dataFile&nbsp;<BR>&nbsp;$fileFormat&nbsp;</TD>
	<TD class="center">&nbsp;$dbStrg&nbsp;<BR>&nbsp;$taxonomy&nbsp;</TD>
	<TD class="center">&nbsp;$minScore&nbsp;<BR>&nbsp;$maxRank&nbsp;</TD>
	<TD class="center">&nbsp;$labeling&nbsp;</TD>
|;
		##>Selected proteins
		if ($valStat>=1) {print "<TD class=\"center\">$anaProteins{$anaID}[0] ($anaProteins{$anaID}[1])</TD>\n";}
		else {print "<TD class=\"center\">$anaProteins{$anaID}[0] / $anaProteins{$anaID}[1]</TD>\n";}

		##>phosphoRS
		if ($callType eq 'phosphoRS') {
			my $pRSstatus;
			if($valStat>1){ # unavailable for phosphoRS
				$pRSstatus = "&nbsp;-&nbsp;";
			}
			elsif ($listParam{$anaID}{'g:Variable modifications'} !~ /Phospho/){
				$pRSstatus = "&nbsp;No phosphorylation&nbsp;";
			}
			else {
				$disabSubmit = '';
				if ($prsParam{$anaID}) {
					$pRSstatus = "<TABLE border=0><TD>";
					while(my ($name,$value) = each %{$prsParam{$anaID}}){
						$pRSstatus .= "&nbsp;$name: $value&nbsp;\n<BR>";
					}
					$pRSstatus =~ s/<BR>$//;
					$pRSstatus .= "</TD><TD align=\"center\" valign=\"center\">&nbsp;<INPUT type=\"button\" value=\"Revert\" class=\"font11\" onclick=\"window.location='$promsPath{cgi}/analysePhospho.cgi?deleteID=$anaID&branchID=$branchID';\">&nbsp;</TD></TABLE>";
				}
				elsif ($prsRunning{$anaID}==1) {
					$pRSstatus = '&nbsp;<FONT style="color:#00A000;font-weight:bold">On-going...</FONT>&nbsp;';
				}
				elsif ($prsRunning{$anaID}==-1) {
					$pRSstatus = '&nbsp;<FONT style="color:#DD0000;font-weight:bold">Error detected</FONT>&nbsp;';
				}
				else {
					$pRSstatus = '&nbsp;Not performed&nbsp;';
				}
			}
			print qq|<TD class="center">$pRSstatus</TD>|;
		}
		##> GO Quanti
		elsif ($callType eq 'goQuantiAnalysis') {
			if(scalar keys %{$quantiData{$anaID}}){
				print "<TD>&nbsp;<SELECT name=\"quanti$anaID\" id=\"quanti$anaID\" onchange=\"loadRatios($anaID, this.options[this.selectedIndex].value);\"\>\n";
				my $selectedQuantiID;
				foreach my $quantiID (keys %{$quantiData{$anaID}}){
					my $selectedString = '';
					unless($selectedQuantiID){
						$selectedQuantiID = $quantiID;
						$selectedString = 'selected';
					}
					print "<OPTION value=\"$quantiID\" $selectedString>$quantiData{$anaID}{$quantiID}{'name'}</OPTION>\n";
				}
				print "</SELECT>&nbsp;</TD>\n";
				print "<TD><SELECT name=\"ratio$anaID\" id=\"ratio$anaID\">\n";
				foreach my $rank (sort{$a<=>$b} keys %{$quantiData{$anaID}{$selectedQuantiID}{'ratio'}}){
					print "<OPTION value=\"$rank\">$quantiData{$anaID}{$selectedQuantiID}{'ratio'}{$rank}</OPTION>\n";
				}
				print "</SELECT>&nbsp;</TD>\n";
			}
			else {
				print qq|<TD class="center">&nbsp;-&nbsp;</TD><TD class="center">&nbsp;-&nbsp;</TD>|;
			}

		}
		elsif ($callType eq 'list' and param('quanti')) {
			print "<TD>&nbsp;$quantiData{$anaID}{'Quantification'}&nbsp;</TD>\n";
			print "<TD>&nbsp;$quantiData{$anaID}{'Ratio'}&nbsp;</TD>\n";
		}
		print "</TR>\n";
		if ($showParam==1) {
			print qq
|<TR valign=top bgcolor=$bgColor>
	<TD colspan=3>
	</TD><TD colspan=6 align=right><FIELDSET><LEGEND><B>Search parameters:</B></LEGEND><TABLE border=0 cellpadding=0 width=100%>
|;
			foreach my $param (sort {$a cmp $b} keys %{$listParam{$anaID}}) {
				(my $trueParam=$param)=~s/^\w://; # remove sort tag
				if ($trueParam eq 'Databank') {
					#>Fetching number of databanks searched
					my $lastIndex=0;
					foreach my $subParam (keys %{$listParam{$anaID}{$param}}) {
						$lastIndex=$#{$listParam{$anaID}{$param}{$subParam}};
						last;
					}
					$trueParam='Multiple databanks' if $lastIndex >= 1;
					print "<TR valign=top><TH align=right nowrap>$trueParam:</TH><TD width=100%>";
					foreach my $i (0..$lastIndex) {
						my $first=1;
						foreach my $subParam (sort {$a cmp $b} keys %{$listParam{$anaID}{$param}}) {
							(my $trueSubParam=$subParam)=~s/^\w://; # remove sort tag
							if ($first) {print "&nbsp;&bull;";} else {print "&nbsp;&nbsp;&nbsp;";}
							print "<B>$trueSubParam:</B>&nbsp;$listParam{$anaID}{$param}{$subParam}[$i]<BR>";
							$first=0;
						}
					}
				}
				else {print "<TR valign=top><TH align=right nowrap>$trueParam:</TH><TD width=100%>&nbsp;$listParam{$anaID}{$param}</TD></TR>\n";}
			}
			print qq
	|	</TABLE></FIELDSET></TD>
	</TR>
	|;
		}
	}
}
else {
	$disabSubmit='';
	foreach my $refDesignData (@itemDesigns) {
		my ($designID,$name,$des)=@{$refDesignData};

		print "<TR valign=middle bgcolor=$bgColor>\n";
		###> CheckBox
		my $disabStrg = ($quantiData{$designID})? '':'disabled';
		my $boxType = 'radio' ;
		my $boxStr = "<INPUT type=\"$boxType\" name=\"desList\" id=\"chk_$designID\" value=\"$designID\" $disabStrg>";
		print "\t<TH valign=middle>$boxStr</TH>\n";

		print "	<TH nowrap align=left valign=middle>&nbsp;<IMG src=\"$promsPath{images}/$itemIcones{'design'}\">&nbsp;$name&nbsp;</TH><TD></TD>\n";
		##> GO Quanti
		if(scalar keys %{$quantiData{$designID}}){
			print "<TD>&nbsp;<SELECT name=\"quanti$designID\" id=\"quanti$designID\" onchange=\"loadRatios($designID, this.options[this.selectedIndex].value);\"\>\n";
			my $selectedQuantiID;
			foreach my $quantiID (keys %{$quantiData{$designID}}){
				my $selectedString = '';
				unless($selectedQuantiID){
					$selectedQuantiID = $quantiID;
					$selectedString = 'selected';
				}
				print "<OPTION value=\"$quantiID\" $selectedString>$quantiData{$designID}{$quantiID}{'name'}</OPTION>\n";
			}
			print "</SELECT>&nbsp;</TD>\n";
			print "<TD><SELECT name=\"ratio$designID\" id=\"ratio$designID\">\n";
			foreach my $rank (sort{$a<=>$b} keys %{$quantiData{$designID}{$selectedQuantiID}{'ratio'}}){
				print "<OPTION value=\"$rank\">$quantiData{$designID}{$selectedQuantiID}{'ratio'}{$rank}</OPTION>\n";
			}
			print "</SELECT>&nbsp;</TD>\n";
		}
		else {
			print qq|<TD class="center">&nbsp;-&nbsp;</TD><TD class="center">&nbsp;-&nbsp;</TD>|;
		}
		print "</TR>\n";
	}
}


if ($callType eq 'report') {
	print qq
|<TR>
	<TD colspan=10><SELECT name="endVal" class="title3" onchange="checkAnalysesForQuantif(this.value)"$disabSubmit>
		<OPTION value=-1>-=Select an action=-</OPTION>
		<OPTION value=0>Report validation</OPTION>
		<OPTION value=2>Report and End validation</OPTION>
		<OPTION value=1>End validation</OPTION>
	</SELECT><INPUT type="submit" name ="Submit" value="Proceed" class="title3"$disabSubmit>&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="cancelAction();"></TD>
</TR>
|;
}
elsif ($callType eq 'clearSel') {
	print qq
|<TR>
	<TD colspan=10><SELECT name="clearMode" class="title3"$disabSubmit>
		<OPTION value=-1>-=Select an action=-</OPTION>
		<OPTION value='clearAuto'>Clear auto-selection</OPTION>
		<OPTION value='clearAll'>Clear all selections</OPTION>
	</SELECT><INPUT type="submit" name="Submit" value="Proceed" class="title3"$disabSubmit>&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="cancelAction();"></TD>
</TR>
|;
}
else {
	if ($callType eq 'impElution') {
		print "<TR bgcolor=\"$darkColor\"><TH colspan=3 align=left valign=top nowrap>&nbsp;&nbsp;Plate template file :</TH><TD colspan=7 nowrap><INPUT type=\"file\" name=\"templateFile\" size=80/><BR><SMALL><B><I>(XML format)</I></B></SMALL></TD></TR>\n";
	}
	print "<TR><TD colspan=10>";
	if ($callType ne 'list') {
		my $submitName=($callType eq 'delete')? 'Delete Analyses' :
			($callType eq 'duplicate')? 'Duplicate Analyses' :
			($callType eq 'appFilter')? 'Go to Filter Settings' :
			($callType eq 'remFilter')? 'Remove Filter' :
			($callType eq 'lowScores')? 'Activate Peptides' :
			($callType eq 'impElution')? 'Import Elution Data' :
			($callType eq 'remFlag')? "Remove Flags" :
			($callType eq 'phosphoRS')? 'Start PhosphoRS' :
			($callType eq 'goQuantiAnalysis')? 'Load Protein Set' :
			($callType eq 'peptidePTM')? 'Show Distribution' :
			'';
		print "<INPUT type=\"submit\" name=\"Submit\" value=\"$submitName\" class=\"title3\"$disabSubmit>&nbsp;&nbsp;";
	}
	print "<INPUT type=\"button\" value=\"Cancel\" onclick=\"cancelAction();\"></TD></TR>\n";
}
print "</TABLE>\n";
print "</FORM>\n" if $callType ne 'list';
print "</BODY>\n</HTML>\n";

####>Revision history<####
# 2.6.7 Add check for on-going phosphoRS analysis (PP 08/11/18)
# 2.6.6 Add branchID to allow phosphoRS parallelization (GA 02/11/18)
# 2.6.5 Add checkbox to overwrite query manually modified (GA 31/01/18)
# 2.6.4 [FIX] bug when 'design' is default value for parameter 'goType' & changed parameter/variable 'goType' to 'srcType' (PP 23/01/18)
# 2.6.3 Restrict GoQuantification quantifs to PROT_RATIO_PEP & minor display change (PP 15/01/18)
# 2.6.2 Compatible with PRS new (xxx_ana_&lt;anaID&gt;.xxx) and old parameter file naming (PP 20/04/17)
# 2.6.1 Added protein and site distribution options in modification distribution form (PP 23/03/17)
# 2.6.0 Added 1 option in modification distribution form (PP 21/04/16)
# 2.5.9 Comparison of same modif on different residues is now possible (PP 17/03/16)
# 2.5.8 Modification filter option for peptide modifications distribution (PP 12/01/16)
# 2.5.7 Allow new report on Proteome Discoverer label-free quantification (PP 22/04/15)
# 2.5.6 Minor modif for $sthCTS{'PATHWAY'} query (GA 27/02/15)
# 2.5.5 Add Analysis selection for peptide modifications distribution (PP 03/12/14)
# 2.5.4 Add Pathway check to Analysis deletablity (PP 13/11/14)
# 2.5.3 compatible with pathway analysis (Sl 21/10/14)
# 2.5.2 Minor modification for Super Ratio Q.GO selection (GA 17/07/14)
# 2.5.1 Add design for Start Q. GO Analysis (GA 14/05/14)
# 2.5.0 Minor modif for option goQuantiAnalysis in pattern matching to retrieve RATIOS (GA 30/04/14)
# 2.4.9 Excludes for decoy proteins from counts (PP 24/03/14)
# 2.4.8 Now authorizing PhosphoRS analyses on PARAGON files (FY 30/09/13)
# 2.4.7 Preventing PARAGON-analysis selection for PhosphoRS (FY 24/09/13)
# 2.4.6 Minor bug fix (PP 09/07/13)
# 2.4.5 Better detection of analysis deletability (PP 05/07/13)
# 2.4.4 Fix display bug in search parameters since multi-databank management (PP 04/07/13)
# 2.4.3 Uses LOWER_SCORES for lower-scoring peptides activation (PP 24/04/13)
# 2.4.2 Filter on iTRAQ and SILAC labeled analyses when starting a GO Q. analysis (FY 29/03/13)
# 2.4.1 Add GO Quanti Analysis dataset selection form and item list (FY 26/11/12)
# 2.4.0 Multi-databank search (PP 18/12/12)
# 2.3.9 Add notification if unselectable analysis for PhosphoRS +<BR> Uninitialized value warning fix (FY 13/04/12)
# 2.3.8 Add phosphoRS starting form (FY 01/03/12)
# 2.3.7 Prevents analysis deletion if 2ndary quantification data (PP 24/11/11)
# 2.3.6 No reporting if exist quantification data (PP 31/10/11)
# 2.3.5 Labeling replaces Quantif. Method (PP 02/08/11)
# 2.3.4 Add GO Analysis item (FY 07/07/11)
# 2.3.3 Minor correction for 'list' action (PP 21/04/11)
# 2.3.2 Add Comparison & Condition checks for analysis deletion (PP 14/03/2011)<BR>Add clear selection and clear flags cases (FY 22/03/2011)
