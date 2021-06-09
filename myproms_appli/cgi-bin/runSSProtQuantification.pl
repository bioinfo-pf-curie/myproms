#!/usr/local/bin/perl -w

################################################################################
# runSSProtQuantification.pl       1.8.0                                       #
# Component of site myProMS Web Server                                         #
# Authors: P. Poullet (Institut Curie)                                         #
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

### myProMS software versions ($MYPROMS_XIC_VERSION) for quantification (based on R scripts & data preprocessing in runXIC...)
### myProMS SSPA versions:
# 1.1: log2 rescaled XIC
# 2.0: Down-shift of log XIC
# 2.1: New "All sets" category for protein/site with adj. p-value > threshold p-val and deltaPC < threshold PC

$| = 1;
use promsConfig;
use promsMod;
use promsQuantif;
use strict;
#use XML::SAX;
use MIME::Base64;
use POSIX qw(strftime);
use File::Copy::Recursive qw(dirmove);
use File::Path;

# exit;
#######################
####>Configuration<####
#######################
my $MYPROMS_SSPA_VERSION='2.1';
#my $REF_XIC_MEDIAN_VALUE=10; # No longer used
my %promsPath=&promsConfig::getServerInfo('no_user');

###############################
####>Recovering parameters<####
###############################
my ($quantifID,$quantifDate)=@ARGV;
my $quantifDir="$promsPath{tmp}/quantification/$quantifDate";
unless (-e $quantifDir) {
	die "ERROR: Job directory '$quantifDate' not found!";
}
my %quantifParameters=&promsQuantif::getQuantificationParameters("$quantifDir/quantif_info.txt");

my $dbh=&promsConfig::dbConnect('no_user');
$dbh->do("UPDATE QUANTIFICATION SET STATUS=0,QUANTIF_ANNOT=REPLACE(QUANTIF_ANNOT,'SOFTWARE=myProMS','SOFTWARE=myProMS;$MYPROMS_SSPA_VERSION') WHERE ID_QUANTIFICATION=$quantifID"); # running
$dbh->commit;

my $runDir="$quantifDir/quanti_$quantifID";
my $dataDir="$runDir/data";
my $resultDir="$runDir/results";
my $graphDir="$resultDir/graph";
#make_path($dataDir,$graphDir,{verbose=>0,mode=>0775}); # $runDir,$resultDir will be created automatically
mkdir $runDir;
mkdir $dataDir;
mkdir $resultDir;
mkdir $graphDir;

# open (DEBUG,">$promsPath{tmp}/quantification/debug_SSPA.txt"); # DEBUG
my $fileStat="$quantifDir/status_$quantifID.out";
open(FILESTAT,">$fileStat") || die "Error while opening $fileStat";
print FILESTAT "Started ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
print FILESTAT "1/4 Fetching data\n";
close FILESTAT;

my $projectID=&promsMod::getProjectID($dbh,$quantifID,'QUANTIFICATION');

#my ($quantifAnnot,$quantifiedModifID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_MODIFICATION FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID");
my ($quantifAnnot,$quantifiedModifID,$quantifModStrg)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,Q.ID_MODIFICATION,GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ',')
																					 FROM QUANTIFICATION Q
																					 LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
																					 WHERE Q.ID_QUANTIFICATION=$quantifID
																					 GROUP BY Q.ID_QUANTIFICATION");
my ($isModifQuantif,$trueModifQuantif,$isFreeResQuantif)=(0,0,0);
my (%modifQuantifRank,%rankQuantifModif,@modifQuantifs);
my $matchQuantifModStrg='';
if ($quantifiedModifID || $quantifModStrg) {
	$isModifQuantif=1;
	$quantifModStrg=$quantifiedModifID if $quantifiedModifID;
	my $modifRank=0;
	foreach my $modID (split(',',$quantifModStrg)) {
		$modifRank++;
		$modifQuantifRank{$modID}=$modifRank;
		$rankQuantifModif{$modifRank}=$modID;
		push @modifQuantifs,$modID; # ordered list of modifs
		if ($modID > 0) {$trueModifQuantif=1;} else {$isFreeResQuantif=1;} # -1 for fake "FreeResidue"
	}
	$quantifiedModifID=$modifQuantifs[0] if $modifRank==1; # single PTM: back compatibility in case QUANTIFICATION.ID_MODIFICATION no longer used
	$matchQuantifModStrg=join('|',@modifQuantifs); # Needed later for regexp on pep varmod
}
my $keepSeqContext=($isModifQuantif && $quantifParameters{'DB'}{'KEEP_SEQ_CONTEXT'} && $quantifParameters{'DB'}{'KEEP_SEQ_CONTEXT'}[0])? 1 : 0; #keep sequence context of PTM (=> Do not combine site from fully cleaved & missed cleaved peptides)

###>Protein enriched by PTMs?
my $isProteinByPTM=0;
my (%proteinPTMcontexts,%peptideContextMatch,%proteinSequence);
if ($quantifParameters{'DB'}{'PROTEIN_PTM'}) { # Restrict peptides to specific PTMs & context
	$isProteinByPTM=1;
	&promsQuantif::preparePTMcontexts(\%proteinPTMcontexts,$quantifParameters{'DB'}{'PROTEIN_PTM'});
}

$quantifAnnot=~s/::/\$/g;
my ($labeling)=($quantifAnnot=~/LABEL=([^\$]+)/);
$labeling=uc($labeling);
my ($stateStrg)=($quantifAnnot=~/STATES=([^\$]+)/);
$stateStrg=~s/#//g; # remove all id tags


################################
####>Get states composition<####
################################
my ($pepSpecificity,$pepMissedCleavage,$ptmFilter,$pepFocus,$skipRescued)=@{$quantifParameters{'DB'}{'PEPTIDES'}};
my (%quantifDesign,%obs2Ana,%peptideQuantifs,%fragQuantifs,%obs2targetPos,$algoType); # SSPA
my (%ana2Obs,%labelModifs,%anaLabelModifs);

my $sthALM=$dbh->prepare("SELECT M.ID_MODIFICATION FROM ANALYSIS_MODIFICATION AM,MODIFICATION M WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND AM.ID_ANALYSIS=? AND M.IS_LABEL=1");
my $sthObsMod=$dbh->prepare("SELECT M.ID_MODIFICATION FROM OBSERVATION O LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION WHERE O.ID_OBSERVATION=?"); # returns NULL if no modif!!!
#my $sthObsBioSamp=$dbh->prepare("SELECT ID_BIOSAMPLE FROM OBSERVATION WHERE ID_OBSERVATION=?");
#my $sthObsMsSamp=$dbh->prepare("SELECT ID_SAMPLE FROM ANALYSIS A WHERE ID_ANALYSIS=?");

my $numStates=0;
foreach my $state (split(/;/,$stateStrg)) { # 'state1;state2;---;stateN',  stateX='nbBioRep,bioRep1.bioRep2.---.bioRepN,#condID', bioRepX='techRep1=techRep2=---=techRepN', techRepX='frac1+frac2+---+fracN'
	my ($numBioRep,$quantiObsIDs,$condID)=split(/,/,$state);
#next if $condID==54; # SILAC ref
	$numStates++;
	@{$quantifDesign{$numStates}}=();
	#my $bioRepCount=0;
	foreach my $bioReplicate (split(/\./,$quantiObsIDs)) {
		#$bioRepCount++;
#last if $bioRepCount==2;
		my @bioRep;
		#my %nameList;
		my $numBioRepObs=0;
		foreach my $techReplicate (split(/&/,$bioReplicate)) {
			my @fractionObs;
#my $numObs=0;
			foreach my $fraction (split(/\+/,$techReplicate)) {
#last if ++$numObs==3;
				$numBioRepObs++;
				my ($obsID,$parQuantifID,$anaID,$targetPos)=split(/:/,$fraction); # $parQuantifID can be 0 for SSPA
				$obs2Ana{$obsID}=$anaID;
				$obs2targetPos{$obsID}=$targetPos;
				push @fractionObs,$obsID;
				if ($pepFocus eq 'xic') {
					unless ($algoType) { # same for all peptide quantifs
						($algoType)=$dbh->selectrow_array("SELECT QM.CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD QM WHERE Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$parQuantifID");
					}
					if ($algoType=~/^(TDA|DIA)$/) { # MS2 XIC
						$fragQuantifs{$parQuantifID}{$anaID}=1; #$obsID;
					}
					else { # MS1 XIC
						$peptideQuantifs{$parQuantifID}{$targetPos}=1; #$obsID;
					}
				}
				#$sthObsBioSamp->execute($obsID);
				#my ($bsID)=$sthObsBioSamp->fetchrow_array;
				#$nameList{BIOSAMPLE}{$bsID}++ if $bsID;
				#$sthObsMsSamp->execute($anaID);
				#my ($sampID)=$sthObsMsSamp->fetchrow_array;
				#$nameList{MSSAMPLE}{$sampID}++;
				%{$labelModifs{$obsID}}=();
				push @{$ana2Obs{$anaID}},$obsID;
				if ($labeling eq 'SILAC') { # SILAC QUANTI
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
				elsif ($labeling =~ /ITRAQ|TMT/) {
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
				elsif ($labeling eq 'FREE') {$anaLabelModifs{$anaID}{0}=1; $labelModifs{$obsID}{0}=1;}
				else {die "ERROR: Unrecognized labeling ('$labeling')!";}
			} # end of fraction
			push @bioRep,\@fractionObs;

		} # end of tech replicate (sample)

		##>Guess most suitable name for bioRep
		#if ($nameList{BIOSAMPLE}) { # Bio-sample
		#	my $bsID=(sort{$nameList{BIOSAMPLE}{$b} <=> $nameList{BIOSAMPLE}{$a} || $a<=>$b} keys %{$nameList{BIOSAMPLE}})[0];
		#	$bioRep{ID}='b'.$bsID if $nameList{BIOSAMPLE}{$bsID}==$numBioRepObs;
		#}
		#unless ($bioRep{ID}) { # MS Sample
		#	my $sampID=(sort{$nameList{MSSAMPLE}{$b} <=> $nameList{MSSAMPLE}{$a} || $a <=> $b} keys %{$nameList{MSSAMPLE}})[0];
		#	$bioRep{ID}='s'.$sampID if $nameList{MSSAMPLE}{$sampID}==$numBioRepObs;
		#	$bioRep{ID}='c'.$bioRepCount unless $bioRep{ID};
		#}
		#$bioRep{LABEL}="State$statePos.BioRep$bioRepCount";
		push @{$quantifDesign{$numStates}},\@bioRep;

	} # end of bio-replicate

} # end of cond/state

$sthALM->finish;
$sthObsMod->finish;
#$sthObsBioSamp->finish;
#$sthObsMsSamp->finish;

###>Handling Modification Quantif phosphoRS/MaxQuant threshold & PTM position ambiguity settings
my ($ptmProbSoftCode,$ptmProbThreshold,$ambiguousModifPos,%quantifResMatchStrg,%qQuantifResMatchStrg)=('',0,0,'',''); # default ($qQuantifResMatchStrg PP 2017/06/22)
# my (%targetableRes,$matchProtCterm,$matchProtNterm); # for recreated PTMs only
if ($isModifQuantif) {
	if ($isFreeResQuantif) { # assumes single-modif quantif!!!!!!!
		# my $context=''; # for now
		# my $modID=$modifQuantifs[0]; # -1 (fake "freeResidus")
		# foreach my $res (@{$quantifParameters{'DB'}{'FREE_RESIDUES'}}) {
		# 	push @{$targetableRes{$modID}},[$res,$context];
		# 	$quantifResMatchStrg{$modID}.=$res;
		# 	$matchProtCterm=1 if ($res eq '+' || $context eq '+'); # for at least 1 modifID
		# 	$matchProtNterm=1 if ($res eq '-' || $context eq '-');
		# }
		&promsQuantif::preparePTMcontexts(\%proteinPTMcontexts,$quantifParameters{'DB'}{'SITE_CONTEXT'});
	}
	if ($trueModifQuantif) { # true PTMs to quantify
		if ($quantifParameters{'DB'}{'PTM_POS'}[0]=~/(\w+):(\d+)/) { # PRS Phospho or MQ any PTM
			$ptmProbSoftCode=$1; # PRS|MQ|SPC|PTMRS
			$ptmProbThreshold=$2; $ptmProbThreshold/=100 if $ptmProbSoftCode =~ /^(MQ|SPC|PTMRS)$/; # 0-1
			$ambiguousModifPos=1 if $quantifParameters{'DB'}{'PTM_POS'}[1] eq 'ambiguous';
		}
		elsif ($quantifParameters{'DB'}{'PTM_POS'}[0] eq 'ambiguous') { # other PTMs valid/ambiguous
			$ambiguousModifPos=1;
		}
	
		#>Find list of targeted residues: needed even if !$ambiguousModifPos
		my $sthAMR=$dbh->prepare("SELECT ID_MODIFICATION,SPECIFICITY FROM ANALYSIS_MODIFICATION WHERE ID_MODIFICATION IN ($quantifModStrg) AND ID_ANALYSIS=? AND MODIF_TYPE='V'");
		my %resList;
		foreach my $anaID (keys %ana2Obs) {
			$sthAMR->execute($anaID);
			my ($modID,$specifStrg)=$sthAMR->fetchrow_array;
			foreach my $res (split(',',$specifStrg)) {$resList{$modID}{$res}=1;}
		}
		$sthAMR->finish;
		my $sthMS=$dbh->prepare('SELECT SPECIFICITY FROM MODIFICATION WHERE ID_MODIFICATION=?');
		foreach my $modID (@modifQuantifs) {
			if ($resList{$modID}) {
				$quantifResMatchStrg{$modID}=join('',sort keys %{$resList{$modID}}); # eg 'STY'
			}
			else { # Fallback in case not specified during search(es)
				$sthMS->execute($modID);
				($quantifResMatchStrg{$modID})=$sthMS->fetchrow_array;
				$quantifResMatchStrg{$modID}=~s/,//g;
			}
		}
	}
	foreach my $modID (keys %quantifResMatchStrg) {
		$qQuantifResMatchStrg{$modID}=quotemeta($quantifResMatchStrg{$modID}); # +* -> \+\* (PP 2017/06/22)
	}
}


#####################################
####>Protein lists for filtering<####
#####################################
my ($protSelectionType,%selectExcludeProteins);
my $sthList=$dbh->prepare("SELECT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=? ORDER BY ID_PROTEIN");

###>Selection/exclusion list
if ($quantifParameters{'DB'}{'PROTEINS'}) {
    ($protSelectionType,my $listID)=@{$quantifParameters{'DB'}{'PROTEINS'}};
    $listID=~s/#//; # remove id flag
    $sthList->execute($listID);
    while (my ($protID)=$sthList->fetchrow_array) {
        $selectExcludeProteins{$protID}=1;
    }
}

###>Contaminants
my %contaminants;
if ($quantifParameters{'DB'}{'EXCL_CONTAMINANTS'} && $quantifParameters{'DB'}{'EXCL_CONTAMINANTS'}[0]) {
	my $sthCont=$dbh->prepare("SELECT DISTINCT ID_PROTEIN FROM ANALYSIS_PROTEIN AP
								INNER JOIN ANALYSIS_DATABANK AD ON AP.ID_ANALYSIS=AD.ID_ANALYSIS
								INNER JOIN DATABANK D ON AD.ID_DATABANK=D.ID_DATABANK
								WHERE AP.ID_ANALYSIS=? AND D.IS_CRAP=1 AND AP.DB_RANK LIKE CONCAT('%',AD.DB_RANK,'%')"); # LIKE because 2 can be in 2 or 12!!!
	foreach my $anaID (keys %ana2Obs) {
		$sthCont->execute($anaID);
		while (my ($protID)=$sthCont->fetchrow_array) {$contaminants{$protID}=1;}
	}
	$sthCont->finish;
}


###################################################
####>Fetching peptide xic values (feature=xic)<####
###################################################
my %pepQuantifValues;
if ($pepFocus eq 'xic') {
	if ($algoType=~/^(TDA|DIA)$/) { # MS2 XIC
		foreach my $ms2Quantif (keys %fragQuantifs) {
			foreach my $anaID (keys %{$fragQuantifs{$ms2Quantif}}) {
				open (FRAG_QUANTIF,"$promsPath{quantification}/project_$projectID/quanti_$ms2Quantif/swath_ana_$anaID.txt") || die $!;
				while (<FRAG_QUANTIF>) {
					next if $.==1;
					s/\s\Z//g; # chomp is not enough. More hidden Windows character
					my ($pepID,$fragMZ,$fragCharge,$fragType,$fragRes,$fragArea,$fragRT)=split(/!/,$_);
					next if $fragArea=~/^(N\/A|None)$/;
					next if $fragArea <= 0;
					$pepQuantifValues{"$pepID:0"}+=$fragArea; # Sum all fragments targetPos=0!!!!
				}
				close FRAG_QUANTIF;
			}
		}
	}
	else { # MS1 XIC
		my $sthPepQP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION Q,QUANTIFICATION_PARAMETER QP WHERE Q.ID_QUANTIFICATION_METHOD=QP.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=? AND (CODE LIKE '%_INTENSITY' OR CODE LIKE '%_AREA')");
		foreach my $xicQuantif (keys %peptideQuantifs) {
			my %usedParamIDs;
			$sthPepQP->execute($xicQuantif);
			while (my ($paramID)=$sthPepQP->fetchrow_array) {
				$usedParamIDs{$paramID}=1;
			}
			my $signalParamID;
			foreach my $targetPos (keys %{$peptideQuantifs{$xicQuantif}}) {
				#my $obsID=$obs2targetPos{$targetPos};
				my $pepQuantifFile=($labeling eq 'FREE')? 'peptide_quantification.txt' : "peptide_quantification_$targetPos.txt";
				open (PEP_QUANTIF,"$promsPath{quantification}/project_$projectID/quanti_$xicQuantif/$pepQuantifFile") || die $!;
				while(<PEP_QUANTIF>) {
					next if $.==1;
					chomp;
					my ($paramID,$pepID,$quantifValue)=split(/\t/,$_);
					#next if $paramID != $peptideAreaParamID{$xicQuantif};
					next if $quantifValue <= 0; # Isobaric bug
					unless ($signalParamID) {
						next unless $usedParamIDs{$paramID};
						$signalParamID=$paramID;
					}
					next if $paramID != $signalParamID;
					$pepQuantifValues{"$pepID:$targetPos"}=$quantifValue;
				}
				close PEP_QUANTIF;
			}
		}
		$sthPepQP->finish;
	}
}

#################################################################
####>Fetching and filtering peptides associated with samples<####
#################################################################
my $useHiddenProt=$quantifParameters{'DB'}{'HIDDEN_PROT'}[0]; # count protein where hidden if visible once
my %visibleProteins;
my $visibilityThres=1;
if ($useHiddenProt) {
	my $sthVisProt=$dbh->prepare("SELECT DISTINCT ID_PROTEIN FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS IN (".join(',',keys %ana2Obs).") AND VISIBILITY >= 1");
	$sthVisProt->execute;
	while (my($protID)=$sthVisProt->fetchrow_array) {$visibleProteins{$protID}=1;}
	$sthVisProt->finish;
	$visibilityThres=0;
}
# Retrieve all visible & hidden proteins but keeps only those visible in at least 1 analysis
my $chargeStrg=($pepFocus=~/sp_count|_ion|xic/)? 'CHARGE' : '0'; # sp_count|all_ion|all_pep|dist_ion|dist_pep|dist_seq (before: peptide or pepSeq)
my $peptideQuery=qq|SELECT GROUP_CONCAT(DISTINCT PPA.ID_PROTEIN,':',ABS(PPA.PEP_BEG),':',AP.VISIBILITY ORDER BY AP.VISIBILITY DESC SEPARATOR '&'),P.ID_PEPTIDE,PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),CHARGE,SPEC_COUNT,MISS_CUT
		FROM PEPTIDE P
		LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
		INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
		INNER JOIN ANALYSIS_PROTEIN AP ON PPA.ID_PROTEIN=AP.ID_PROTEIN
		WHERE P.ID_ANALYSIS=AP.ID_ANALYSIS AND AP.ID_ANALYSIS=? AND AP.VISIBILITY >= $visibilityThres|; # No need for group_concat on PEP_BEG since only 1 position is usable
	$peptideQuery.=($pepSpecificity eq 'unique')? ' AND PPA.IS_SPECIFIC=1' : ($pepSpecificity eq 'unique_shared')? ' AND AP.PEP_SPECIFICITY=100' : ''; # Filter at DB level
	$peptideQuery.=' AND MISS_CUT=0' if $pepMissedCleavage==0; # 0: exclude only missed-cut peptides, 1: allow, -1: also exclude overlapping fully cut
	$peptideQuery.=' AND P.VALID_STATUS > 0' if $skipRescued; # No ghost/MBR/rescued peptides
	$peptideQuery.=' GROUP BY P.ID_PEPTIDE'; # PPA.ID_PROTEIN,
#print DEBUG "\nQUERY=$peptideQuery\n";
my $sthGetPep=$dbh->prepare($peptideQuery);

$ptmFilter=~s/#//g; # remove id flag
my ($ptmAllowed,@selectedPTMs)=split(/[:,]/,$ptmFilter); # @selectedPTMs defined only if $ptmAllowed >= 2
my %allowedPtmID;
foreach my $modID (@selectedPTMs) {$allowedPtmID{$modID}=1;}
foreach my $modID (@modifQuantifs) {
	$allowedPtmID{$modID}=1; # to allow modif quantification
}

##>Pre-processing in case $pepMissedCleavage = -1
my %missedCutExcludedSeq; # only if $pepMissedCleavage=-1 (exclude missCut AND overlapping cut peptides)
if ($pepMissedCleavage==-1) {
	my (%missedCutPeptides,%beg2pepSeq,%seenData);
	foreach my $statePos (keys %quantifDesign) {
		foreach my $refBioRep (@{$quantifDesign{$statePos}}) {
			foreach my $refTechRep (@{$refBioRep}) {
				foreach my $obsID (@{$refTechRep}) { # 1 obs = 1 fraction
					my $anaID=$obs2Ana{$obsID};
					my (%missingSpCount,%obsIonCount);
					if (!$seenData{$anaID}) {
						$sthGetPep->execute($anaID);
						while (my ($protStrg,$pepID,$pepSeq,$varModStrg,$charge,$spCount,$misscut)=$sthGetPep->fetchrow_array) {
							my $pepLength;
							if ($misscut) {
								$missedCutExcludedSeq{$pepSeq}=1;
								$pepLength=length($pepSeq);
							}
							my @protList;
							foreach my $protData (split('&',$protStrg)) {
								my ($protID,$pepBeg,$visibility)=split(':',$protData);
								last unless defined $visibility; # truncated string if peptide matches too many proteins (response max length=1024 chars)
								next if ($visibility==0 && !$visibleProteins{$protID}); # $visibility is 0 => must be $useHiddenProt=1 
								if ($misscut) {
									$missedCutPeptides{$protID}{$pepBeg}=$pepLength if (!$missedCutPeptides{$protID} || !$missedCutPeptides{$protID}{$pepBeg} || $missedCutPeptides{$protID}{$pepBeg} < $pepLength);
								}
								else {$beg2pepSeq{$protID}{$pepBeg}=$pepSeq;} # fully cut pep seq (only 1 seq per beg)
							}
						}
						$seenData{$anaID}=1;
					}
				}
			}
		}
	}
	##>Filter fully cut peptides included in missed-cut. Match is CROSS-ANALYSIS!
	foreach my $protID (keys %missedCutPeptides) {
		next unless $beg2pepSeq{$protID}; # no fully cut peptides for protID
		foreach my $missBeg (keys %{$missedCutPeptides{$protID}}) {
			my $missEnd=$missBeg + $missedCutPeptides{$protID}{$missBeg} - 1;
			foreach my $pepBeg (keys %{$beg2pepSeq{$protID}}) {
				if ($pepBeg >= $missBeg && $pepBeg <= $missEnd) { # match! beg overlaps with missed cut seq
					my $pepSeq=$beg2pepSeq{$protID}{$pepBeg};
					$missedCutExcludedSeq{$pepSeq}=1;
					delete $beg2pepSeq{$protID}{$pepSeq};
					delete $beg2pepSeq{$protID} unless scalar keys %{$beg2pepSeq{$protID}};
				}
			}
		}
	}
}


##>Data file
my %recreatedVarMods; # to prevent re-creating same varMod multiple times
# my $sthProtL=$dbh->prepare('SELECT PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=?');; # Only if sites to be recreated at protein C-term
my $sthProtSeq=$dbh->prepare('SELECT P.PROT_SEQ,MP.PROT_SEQ FROM PROTEIN P LEFT JOIN MASTER_PROTEIN MP ON P.ID_MASTER_PROTEIN=MP.ID_MASTER_PROTEIN WHERE P.ID_PROTEIN=?');
# my %proteinLength; # Only if sites to be recreated at protein C-term
my (@bioRepOrder,%anaData,%proteinData,%bestVisibility,%excludedSeq);
my (%missedCutPeptides,%beg2pepSeq); # only if $pepMissedCleavage=-1 (exclude missCut AND overlapping cut peptides)
my (%ionAllXicValues,%tmpProteinData); # for pepFocus=xic ,%stateBioReps
foreach my $statePos (sort{$a<=>$b} keys %quantifDesign) {
	my $bioRepPos=0;
	foreach my $refBioRep (@{$quantifDesign{$statePos}}) {
		$bioRepPos++;
#$stateBioReps{$statePos}=$bioRepPos; # for xic only
		my $bioRepLabel="State$statePos.BioRep$bioRepPos";
		push @bioRepOrder,$bioRepLabel;

		####>Data
		my %techRepData;
		my $numTechRep=0;
		foreach my $refTechRep (@{$refBioRep}) {
			$numTechRep++;
			my (%fracSpCount,%fracXic,%protModProt);
			foreach my $obsID (@{$refTechRep}) { # 1 obs = 1 fraction
				%{$fracSpCount{$obsID}}=(); # defined so no need to checked later
				%{$fracXic{$obsID}}=(); # defined so no need to checked later
				my $anaID=$obs2Ana{$obsID};
				my $targetPos=$obs2targetPos{$obsID};
				my (%missingSpCount,%obsIonCount);
				#my @pepIdList;
				if (!$anaData{$anaID}) {
					$sthGetPep->execute($anaID);
					PEP:while (my ($protStrg,$pepID,$pepSeq,$varModStrg,$charge,$spCount,$misscut)=$sthGetPep->fetchrow_array) {
						next if $excludedSeq{$pepSeq}; # skip any future seq
						next if ($pepMissedCleavage==-1 && $missedCutExcludedSeq{$pepSeq});
						next if ($pepFocus eq 'xic' && !$pepQuantifValues{"$pepID:$targetPos"});
			
						$varModStrg='' unless $varModStrg;
						my (@protList,%peptideBeg,%protVarModStrg);
						foreach my $protData (split('&',$protStrg)) {
							my ($protID,$pepBeg,$visibility)=split(':',$protData);
							last unless defined $visibility; # truncated string if peptide matches too many proteins (response max length=1024 chars)
							next if ($visibility==0 && !$visibleProteins{$protID}); # $visibility is 0 => must be $useHiddenProt=1

							#>Protein filtering
							next if $contaminants{$protID}; # only if option was checked
							if ($protSelectionType) {
								if ($protSelectionType eq 'exclude') {next if $selectExcludeProteins{$protID};}
								else {next unless $selectExcludeProteins{$protID};} # restrict
							}

							push @protList,$protID;
							$peptideBeg{$protID}=$pepBeg; # local to a specific pepSeq
							#$bestVisibility{$protID}=$visibility if (!$bestVisibility{$protID} || $bestVisibility{$protID} < $visibility);
						}
						next unless scalar @protList; # protein(s) did not pass filter(s)

						#my $trueVarModStrg=$varModStrg;
						my @validProtList;
						if ($isFreeResQuantif) {
							# if ($matchProtNterm || $matchProtCterm) { # if create protein N/C-term PTM, same $varModStrg can be different for each protein
							# 	foreach my $protID (@protList) {
							# 		if ($matchProtCterm && !defined $proteinLength{$protID}) {
							# 			$sthProtL->executed($protID);
							# 			($proteinLength{$protID})=$sthProtL->fetchrow_array;
							# 		}
							# 		my $protVmod=&promsQuantif::modifyFreeResidues($pepSeq,$varModStrg,$peptideBeg{$protID},$proteinLength{$protID},\%targetableRes,\%recreatedVarMods);
							# 		next if (!$protVmod || $protVmod !~ /(^|&)($matchQuantifModStrg):/ || $protVmod=~/\d:(\.|&|$)/); # next $protID!!! missing mod pos data!
							# 		push @validProtList,$protID;
							# 		$protVarModStrg{$protID}=$protVmod; # local to a specific pepSeq
							# 	}
							# }
							# else { # same recreated vMod for all proteins
							# 	my $protID0=$protList[0];
							# 	if ($matchProtCterm && !defined $proteinLength{$protID0}) {
							# 		$sthProtL->executed($protID0);
							# 		($proteinLength{$protID0})=$sthProtL->fetchrow_array;
							# 	}
							# 	my $protVmod=&promsQuantif::modifyFreeResidues($pepSeq,$varModStrg,$peptideBeg{$protID0},0,\%targetableRes,\%recreatedVarMods);
							# 	next if (!$protVmod || $protVmod !~ /(^|&)($matchQuantifModStrg):/ || $protVmod=~/\d:(\.|&|$)/); # missing mod pos data!
							# 	foreach my $protID (@protList) {
							# 		push @validProtList,$protID;
							# 		$protVarModStrg{$protID}=$protVmod; # local to a specific pepSeq
							# 	}
							# }
							foreach my $protID (@protList) {
								my $protVmod=&promsQuantif::modifyFreeResidues($pepSeq,$varModStrg,$peptideBeg{$protID},$protID,$sthProtSeq,\%proteinPTMcontexts,\%recreatedVarMods,\%proteinSequence);									push @validProtList,$protID;
								next if $protVmod eq $varModStrg; # no site creation
								push @validProtList,$protID;
								$protVarModStrg{$protID}=$protVmod; # local to a specific pepSeq
							}
						}
						else {
							next if (($isModifQuantif && (!$varModStrg || $varModStrg !~ /(^|&)($matchQuantifModStrg):/)) || ($varModStrg && $varModStrg=~/\d:(\.|&|$)/)); # missing mod pos data!
							foreach my $protID (@protList) {
								if ($isProteinByPTM) {
									my $okPeptide=&promsQuantif::checkPTMcontext($pepSeq,$varModStrg,$peptideBeg{$protID},$protID,$sthProtSeq,\%proteinPTMcontexts,\%peptideContextMatch,\%proteinSequence);
									next unless $okPeptide;
								}
								push @validProtList,$protID;	
								$protVarModStrg{$protID}=$varModStrg; # local to a specific pepSeq
							}
						}
						next unless scalar @validProtList;
						@protList=@validProtList; # exclude non-matching proteins

						$spCount=0 unless $spCount; # in case only ghost (not defined), set to 1 if feature is not sp_count
						
						#>Processing var mods & matching label channel (checking allowed & removing label if labeled quantif)
						my @varMods=($varModStrg)? split(/&/,$varModStrg) : (); # Not modified by FREE_RESIDUES
						#my $varModCode='';
						##if ($labeling eq 'SILAC') {
						my $matchedChannel=0;
						foreach my $anaObsID (@{$ana2Obs{$anaID}}) {
							my $matchedObs=0;
							if ($labeling eq 'FREE') {
								$matchedObs=1;
							}
							else {
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
							}
							if ($matchedObs) {
								my $hasExcludedPTM=0;
								
								foreach my $protID (@protList) {
									my $varModCode='';
									my @protVarMods=($protVarModStrg{$protID})? split(/&/,$protVarModStrg{$protID}) : ();
									foreach my $vMod (@protVarMods) {
										my ($modID)=$vMod=~/^(-*\d+)/; # -1 for Free residues quantif
# unless (defined $modID) {
# 	print DEBUG "vMod='$vMod', protVarModStrg{$protID}='$protVarModStrg{$protID}', varModStrg='$varModStrg', pepID=$pepID\n";
# }
										next if $labelModifs{$anaObsID}{$modID}; # skip labeling modifs
										next if $proteinPTMcontexts{$modID}; # skip context PTMs
										if ($ptmAllowed != 1 && !$allowedPtmID{$modID}) {
											$hasExcludedPTM=1;
											last;
										}
										$varModCode.='&' if $varModCode;
										$varModCode.=$vMod;
									}
									if ($hasExcludedPTM) { # this modification is not allowed
										if ($ptmAllowed <= -1) { # -1 or -2 unmodified peptide not allowed if exists a modified version (unmodif pep may have been recorded before detection => clean below)
											foreach my $protID2 (@protList) {$excludedSeq{$pepSeq}{$protID2}=1;}
										}
										next PEP; # skip all other obs for this ana
									}
									$protVarModStrg{$protID}=$varModCode;		
								}
								
								if ($anaObsID == $obsID) {$matchedChannel=1;}
								else { # stored to be recalled when proper obsID is wanted
									push @{$anaData{$anaID}{$anaObsID}},[\@protList,\%peptideBeg,$pepID,$pepSeq,\%protVarModStrg,$charge,$spCount];
								}
								last; # obs in Ana: current anaObs was matched by peptide: no need to scan others
							}
						}

						###>Recording peptides
						if ($matchedChannel) {
							#my $seqVarMod="$pepSeq:$varModCode";
							foreach my $protID (@protList) {
								my $seqVarMod="$pepSeq:$protVarModStrg{$protID}";
								if (!$protModProt{$protID} || !$protModProt{$protID}{$seqVarMod}) {
									$protModProt{$protID}{$seqVarMod}=($isModifQuantif)? &createModProtID($protID,$pepSeq,$protVarModStrg{$protID},$peptideBeg{$protID},\@modifQuantifs) : $protID;
								}
								$protVarModStrg{$protID}='' if $pepFocus eq 'dist_seq';
							}
							#$varModCode='' if $pepFocus eq 'dist_seq';
							#my $ionKey="$pepSeq:$varModCode:$charge";
							if ($pepFocus eq 'sp_count') { # record best spCount/ion in each fraction
								foreach my $protID (@protList) {
									my $ionKey="$pepSeq:$protVarModStrg{$protID}:$charge";
									my $modProtID=$protModProt{$protID}{"$pepSeq:$protVarModStrg{$protID}"};
									#$techRepData{$protID}{$numTechRep}{"$pepSeq:$varModCode:$charge"}=$spCount if (!$techRepData{$protID} || !$techRepData{$protID}{$numTechRep} || !$techRepData{$protID}{$numTechRep}{"$pepSeq:$varModCode:$charge"} || $techRepData{$protID}{$numTechRep}{"$pepSeq:$varModCode:$charge"} < $spCount);
									$fracSpCount{$obsID}{$modProtID}{$ionKey}=$spCount if (!$fracSpCount{$obsID}{$modProtID} || !defined $fracSpCount{$obsID}{$modProtID}{$ionKey} || $fracSpCount{$obsID}{$modProtID}{$ionKey} < $spCount);
									$missingSpCount{$modProtID}{$ionKey}=1 if $spCount==0;
									$obsIonCount{$modProtID}{$ionKey}++;
								}
							}
							elsif ($pepFocus eq 'xic') {
								foreach my $protID (@protList) {
									my $modProtID=$protModProt{$protID}{"$pepSeq:$protVarModStrg{$protID}"};
									$fracXic{$obsID}{$modProtID}{"$pepSeq:$protVarModStrg{$protID}:$charge"}=$pepQuantifValues{"$pepID:$targetPos"};
									#$fracXic{$obsID}{$protID}{"$pepSeq:$protVarModStrg{$protID}:$charge"}=log($pepQuantifValues{"$pepID:$targetPos"})/$log2; #$logValue;
								}
#$techRepData{$protID}{$numTechRep}{"$pepSeq:$protVarModStrg{$protID}:$charge"}+=$pepQuantifValues{"$pepID:$obsID"}; # sum fractions						
#push @{$ionAllXicValues{"$pepSeq:$protVarModStrg{$protID}:$charge"}},$pepQuantifValues{"$pepID:$obsID"}; # records all pep XICs associated with this ion (needed later for xic rescaling)
							}
							else { # all|distinct ion|pep|seq
								foreach my $protID (@protList) {
									my $modProtID=$protModProt{$protID}{"$pepSeq:$protVarModStrg{$protID}"};
									$techRepData{$modProtID}{$numTechRep}{"$pepSeq:$protVarModStrg{$protID}:$charge"}++;
								}
#@pepIdList,$pepID if $pepFocus eq 'xic';
							}
						}
					}
				}
				else { # analysis data already fetched & assigned to an Obs =>  $anaData{$anaID}{$obsID} defined
					foreach my $refPep (@{$anaData{$anaID}{$obsID}}) {
						my ($refProtList,$refPeptideBeg,$pepID,$pepSeq,$refProtVarModStrg,$charge,$spCount)=@{$refPep};
						#my $seqVarMod="$pepSeq:$varModCode";
						foreach my $protID (@{$refProtList}) {
							my $seqVarMod="$pepSeq:".$refProtVarModStrg->{$protID};
							if (!$protModProt{$protID} || !$protModProt{$protID}{$seqVarMod}) {
								$protModProt{$protID}{$seqVarMod}=($isModifQuantif)? &createModProtID($protID,$pepSeq,$refProtVarModStrg->{$protID},$refPeptideBeg->{$protID}{$pepSeq},\@modifQuantifs) : $protID;
							}
							$refProtVarModStrg->{$protID}='' if $pepFocus eq 'dist_seq';
						}
						#$varModCode='' if $pepFocus eq 'dist_seq';
						#my $ionKey="$pepSeq:$varModCode:$charge";
						if ($pepFocus eq 'sp_count') { # record best spCount for each ion in each fraction (cannot be directly summed because of multiple occurences of same ion with same spCount)
							foreach my $protID (@{$refProtList}) {
								my $ionKey="$pepSeq:".$refProtVarModStrg->{$protID}.":$charge";
								my $modProtID=$protModProt{$protID}{"$pepSeq:".$refProtVarModStrg->{$protID}};
								$fracSpCount{$obsID}{$modProtID}{$ionKey}=$spCount if (!$fracSpCount{$obsID}{$modProtID} || !defined $fracSpCount{$obsID}{$modProtID}{$ionKey} || $fracSpCount{$obsID}{$modProtID}{$ionKey} < $spCount);
								$missingSpCount{$modProtID}{$ionKey}=1 if $spCount==0;
								$obsIonCount{$modProtID}{$ionKey}++;
							}
						}
						elsif ($pepFocus eq 'xic') {
							foreach my $protID (@{$refProtList}) {
								my $modProtID=$protModProt{$protID}{"$pepSeq:".$refProtVarModStrg->{$protID}};
								$fracXic{$obsID}{$modProtID}{"$pepSeq:".$refProtVarModStrg->{$protID}.":$charge"}=$pepQuantifValues{"$pepID:$targetPos"};
							}
						}
						else { # all|distinct ion|pep|seq
							foreach my $protID (@{$refProtList}) {
								my $modProtID=$protModProt{$protID}{"$pepSeq:".$refProtVarModStrg->{$protID}};
								$techRepData{$modProtID}{$numTechRep}{"$pepSeq:".$refProtVarModStrg->{$protID}.":$charge"}++;
							}
						}
					}
				}
				##>Replace sp_count with all_ions where not enough sp_count
				if ($pepFocus eq 'sp_count') {
					foreach my $modProtID (keys %missingSpCount) {
						foreach my $ionKey (keys %{$missingSpCount{$modProtID}}) {
							$fracSpCount{$obsID}{$modProtID}{$ionKey}=$obsIonCount{$modProtID}{$ionKey} if $fracSpCount{$obsID}{$modProtID}{$ionKey} < $obsIonCount{$modProtID}{$ionKey};
						}
					}
				}
			} # End of obs/fraction

			###>Sum each ion's spCount across fractions
			if ($pepFocus eq 'sp_count') {
				foreach my $obsID (keys %fracSpCount) {
					foreach my $modProtID (keys %{$fracSpCount{$obsID}}) {
						foreach my $ionKey (keys %{$fracSpCount{$obsID}{$modProtID}}) {
							$techRepData{$modProtID}{$numTechRep}{$ionKey}+=$fracSpCount{$obsID}{$modProtID}{$ionKey};
						}
					}
				}
			}
			elsif ($pepFocus eq 'xic') {
				my %techRepXic;
				foreach my $obsID (keys %fracXic) { # sum fractions if any (combine all obs in techRep)
					foreach my $modProtID (keys %{$fracXic{$obsID}}) {
						foreach my $ionKey (keys %{$fracXic{$obsID}{$modProtID}}) {
							$techRepData{$modProtID}{$numTechRep}{$ionKey}+=$fracXic{$obsID}{$modProtID}{$ionKey};
							$techRepXic{$ionKey}+=$fracXic{$obsID}{$modProtID}{$ionKey};
						}
					}
				}
				foreach my $ionKey (keys %techRepXic) {
					push @{$ionAllXicValues{$ionKey}},$techRepXic{$ionKey}; # rescale must occur on XIC from summed fractions
				}
			}

			###>Clean excluded pepSeq ($ptmAllowed <= -1)
			foreach my $pepSeq (keys %excludedSeq) {
				my $firstProt=1;
				foreach my $protID (keys %{$excludedSeq{$pepSeq}}) {
					foreach my $seqVarMod (grep{/^$pepSeq:/} keys %{$protModProt{$protID}}) {
						my $modProtID=$protModProt{$protID}{$seqVarMod};
						#my @peps=grep{/^$pepSeq:/} keys %{$techRepData{$protID}{$numTechRep}};
						foreach my $ionKey (grep{/^$pepSeq:/} keys %{$techRepData{$modProtID}{$numTechRep}}) {
							delete $techRepData{$modProtID}{$numTechRep}{$ionKey};
							delete $techRepData{$modProtID}{$numTechRep} unless scalar keys %{$techRepData{$modProtID}{$numTechRep}};
							delete $techRepData{$modProtID} unless scalar keys %{$techRepData{$modProtID}};
							delete $ionAllXicValues{$ionKey} if ($pepFocus eq 'xic' && $firstProt);
						}
						$firstProt=0;
					}
				}
				%{$excludedSeq{$pepSeq}}=(); # reset list of matching prot to 0 before next techRep
			}

		} # End of techRep

		###>Counting peptides (aggregate techReps -> mean of pepCount)
		if ($pepFocus=~/^(sp|all)_/) { # sp_count | all_(ion|pep)
			foreach my $modProtID (keys %techRepData) {
				foreach my $techRepPos (keys %{$techRepData{$modProtID}}) {
					foreach my $count (values %{$techRepData{$modProtID}{$techRepPos}}) {
						$proteinData{$modProtID}{$bioRepLabel}+=$count;
					}
				}
				$proteinData{$modProtID}{$bioRepLabel}/=$numTechRep;
			}
		}
		elsif ($pepFocus eq 'xic') { # Average techReps here! ### No longer true => "Divide by num techRep here & compute mean for bioRep later by just summing"
			foreach my $modProtID (keys %techRepData) {
				foreach my $techRepPos (keys %{$techRepData{$modProtID}}) {
					foreach my $ionKey (keys %{$techRepData{$modProtID}{$techRepPos}}) {
						#$techRepData{$modProtID}{$techRepPos}{$ionKey}/=$numTechRep;
						$tmpProteinData{$bioRepLabel}{$modProtID}{$ionKey}+=$techRepData{$modProtID}{$techRepPos}{$ionKey};
					}
				}
				foreach my $ionKey (keys %{$tmpProteinData{$bioRepLabel}{$modProtID}}) {
					$tmpProteinData{$bioRepLabel}{$modProtID}{$ionKey}/=$numTechRep;
				}
			}
			#$tmpProteinData{$bioRepLabel}=\%techRepData;
		}
		else { # ion or peptide or pepSeq (distinct +/- charge)
			foreach my $modProtID (keys %techRepData) {
				foreach my $techRepPos (keys %{$techRepData{$modProtID}}) {
					$proteinData{$modProtID}{$bioRepLabel}+=scalar keys %{$techRepData{$modProtID}{$techRepPos}};
				}
				$proteinData{$modProtID}{$bioRepLabel}/=$numTechRep;
			}
		}
	} # End of bioRep

} # End of State
# $sthProtL->finish;
$sthProtSeq->finish;
$sthGetPep->finish;

####>Log-transformation & Geometric sum of XICs<####
if ($pepFocus eq 'xic') {
	undef %pepQuantifValues; # no longer needed

	my $logMod=($isModifQuantif)? log(2) : log(10); # <======= IMPORTANT CHANGE!!! (PP 03/12/20)

##>Computing min threshold with 0.1% smallest values for down shift of logged values
my @allXicValues;
foreach my $ionKey (keys %ionAllXicValues) {
	push @allXicValues,@{$ionAllXicValues{$ionKey}};
}
undef %ionAllXicValues;
my $numMinValues=int(0.5 + 0.001 * scalar @allXicValues); $numMinValues=1 if $numMinValues==0;
my @sortedMinValues=(sort{$a<=>$b} @allXicValues)[0..($numMinValues-1)];
my $logThreshold=0;
foreach my $val (@sortedMinValues) {$logThreshold+=log($val)/$logMod;}
$logThreshold=($logThreshold/$numMinValues)-1; # geometric mean - 1

	##>Correcting all values
	foreach my $bioRepLabel (keys %tmpProteinData) {
		if ($isModifQuantif) { # log-transformed max ion!!!! & down shift
			foreach my $modProtID (keys %{$tmpProteinData{$bioRepLabel}}) {
				my $maxXic=(sort{$b<=>$a} values %{$tmpProteinData{$bioRepLabel}{$modProtID}})[0];
				my $corrXic=(log($maxXic)/$logMod)-$logThreshold;
				$corrXic=1 if $corrXic < 1;
				$proteinData{$modProtID}{$bioRepLabel}+=$corrXic;
			}
		}
		else { # Whole protein: Geometric sum of all down shifted ions
			foreach my $protID (keys %{$tmpProteinData{$bioRepLabel}}) {
				foreach my $ionKey (keys %{$tmpProteinData{$bioRepLabel}{$protID}}) {
					my $corrXic=(log($tmpProteinData{$bioRepLabel}{$protID}{$ionKey})/$logMod)-$logThreshold;
					$corrXic=1 if $corrXic < 1;
					$proteinData{$protID}{$bioRepLabel}+=$corrXic;
				}
			}
		}
	}
}
####>Log-transformation, Rescaling & Summing XICs (for feature=xic)<####
if ($pepFocus eq 'xic_v2') {
	undef %pepQuantifValues; # no longer needed
	
	###>Computing lower outlier threshold based on IQR (boxplot)<###  <---- No rescaling applied. Uncomment block below to reapply (PP 14/10/20)
	my $logMod=($isModifQuantif)? log(2) : log(10);
	my @log2XicValues;
	foreach my $ionKey (keys %ionAllXicValues) {
		push @log2XicValues,map {log($_)/$logMod} @{$ionAllXicValues{$ionKey}};
	}
	undef %ionAllXicValues;
	my @sortedLog2XicValues=sort{$a<=>$b} @log2XicValues;
	my $xicLength=scalar @sortedLog2XicValues;
	my $midLength=0.5*$xicLength; # 1/2
	my $medianValue;
	if ($xicLength % 2) { # odd number
		$medianValue=$sortedLog2XicValues[$midLength]; # int => index
	}
	else { # even number: take mean of boundary values
		$midLength=int($midLength);
		$medianValue=($sortedLog2XicValues[$midLength-1] + $sortedLog2XicValues[$midLength]) / 2;
	}
	my $boxMinRk=0.25*$xicLength; # 1/4
	my $boxMaxRk=0.75*$xicLength; # 3/4
	my $lowerQuartile=(int($boxMinRk)==$boxMinRk)? $sortedLog2XicValues[$boxMinRk] : ($sortedLog2XicValues[int($boxMinRk)] + $sortedLog2XicValues[int($boxMinRk)+1]) / 2;
	my $upperQuartile=(int($boxMaxRk)==$boxMaxRk)? $sortedLog2XicValues[$boxMaxRk] : ($sortedLog2XicValues[int($boxMaxRk)] + $sortedLog2XicValues[int($boxMaxRk)+1]) / 2;
	my $minLogXicValue=$lowerQuartile - (($upperQuartile-$lowerQuartile) * 1.5); # wisker = 1.5 box size
	$minLogXicValue-=1; # smallest non-outlier value ~ 1
	
	###>Rescaling & summing XICs for each protein<###
	foreach my $bioRepLabel (keys %tmpProteinData) {
		my %techMean;
		my $refTechRepData=$tmpProteinData{$bioRepLabel};
		foreach my $modProtID (keys %{$refTechRepData}) {
			foreach my $techRepPos (keys %{$refTechRepData->{$modProtID}}) {
				foreach my $ionKey (keys %{$refTechRepData->{$modProtID}{$techRepPos}}) {
					#$techMean{$modProtID}+=log($refTechRepData->{$modProtID}{$techRepPos}{$ionKey})/$logMod; # geometric mean
					$techMean{$modProtID}+=log($refTechRepData->{$modProtID}{$techRepPos}{$ionKey})/$logMod; # Mean here be
				}
			}
			$techMean{$modProtID}/=scalar keys %{$refTechRepData->{$modProtID}}; # geometric mean
		}
		foreach my $modProtID (keys %techMean) {
			#$proteinData{$modProtID}{$bioRepLabel}=sprintf '%.0f',$techMean{$modProtID}; #int(0.5 + 10*$techMean{$modProtID})/10; #1/10 precision
#my $scaledValue=sprintf '%.0f',(log($techMean{$modProtID})/$logMod)-$minLogXicValue;
			my $scaledValue=$techMean{$modProtID}-$minLogXicValue;
			#$proteinData{$modProtID}{$bioRepLabel}=($scaledValue < 0)? 0 : $scaledValue;
			$proteinData{$modProtID}{$bioRepLabel}=($scaledValue < 1)? 1 : $scaledValue; # Smallest values = 1
		}
	}
	$dbh->do("UPDATE QUANTIFICATION SET QUANTIF_ANNOT=CONCAT(QUANTIF_ANNOT,'::LOG2_XIC_THRESHOLD=$minLogXicValue') WHERE ID_QUANTIFICATION=$quantifID");
	$dbh->commit;
}
if ($pepFocus eq 'xic_v1') { # WARNING: compatibility with data variables must be verified
	undef %pepQuantifValues; # no longer needed
	
	###>Computing scale for each ion<###
	#>Set median value to $REF_XIC_MEDIAN_VALUE
	##my %ionScale;
#my $minStateBioRep=(sort{$a<=>$b} values %stateBioReps)[0];
#$minStateBioRep/=2; # half of smallest set

##my $minLogXicValue=my $min2Value=my $min1Value=999999;
my @allXicValues;
	foreach my $ionKey (keys %ionAllXicValues) {
push @allXicValues,@{$ionAllXicValues{$ionKey}};
		##my @sortedXICs=sort{$a<=>$b} @{$ionAllXicValues{$ionKey}};
		##delete $ionAllXicValues{$ionKey};
		##my $xicLength=scalar @sortedXICs;
		##my $midLength=$xicLength/2;
		##my $medianValue;
		##if ($xicLength % 2) { # odd number
		##	$medianValue=$sortedXICs[int($midLength)]; # int => index
		##}
		##else { # even number: take mean of boundary values
		##	$medianValue=( $sortedXICs[$midLength-1] + $sortedXICs[$midLength] ) / 2;
		##}
		####$ionScale{$ionKey}=($medianValue)? $REF_XIC_MEDIAN_VALUE/$medianValue : 1;
		##if ($medianValue) {
		##	$min1Value=$medianValue if ($min1Value > $medianValue);
		##	$min2Value=$medianValue if ($xicLength >= 2 && $min2Value > $medianValue); # at least 3 values
		##	$minLogXicValue=$medianValue if ($xicLength >= 3 && $minLogXicValue > $medianValue); # at least 3 values
		##}
	}
	undef %ionAllXicValues;
	##if ($minLogXicValue==999999) {
	##	$minLogXicValue=($min2Value < 999999)? $min2Value : $min1Value;
	##}
	
	###>Compute mean of 0.1% (or 100 max) smallest values of entire set
my $logMod=($isModifQuantif)? log(2) : log(10);
my $numXicValues=scalar @allXicValues;
my $minPc=int(0.5 + $numXicValues / 1000);
$minPc=100 if $minPc > 100; $minPc=1 if $minPc < 1;
my @minXicValues=(sort{$a<=>$b} @allXicValues)[0..($minPc-1)];
my $minLogXicValue=0;
foreach my $val (@minXicValues) {$minLogXicValue+=$val;}
$minLogXicValue/=$minPc;
#$minLogXicValue-=1; # => smallest median rescaled to 1
$minLogXicValue=(log($minLogXicValue)/$logMod) - 1; # => smallest prot value (1pep) will be rescaled to ~1
	
	###>Rescaling & summing XICs for each protein<###
	foreach my $bioRepLabel (keys %tmpProteinData) {
		my %techMean;
		my $refTechRepData=$tmpProteinData{$bioRepLabel};
		foreach my $protID (keys %{$refTechRepData}) {
			foreach my $techRepPos (keys %{$refTechRepData->{$protID}}) {
				foreach my $ionKey (keys %{$refTechRepData->{$protID}{$techRepPos}}) {
					#$techMean{$protID}+=($refTechRepData->{$protID}{$techRepPos}{$ionKey} * $ionScale{$ionKey});
#$techMean{$protID}+=$refTechRepData->{$protID}{$techRepPos}{$ionKey};
$techMean{$protID}+=log($refTechRepData->{$protID}{$techRepPos}{$ionKey})/$logMod;
				}
			}
			$techMean{$protID}/=scalar keys %{$refTechRepData->{$protID}};
#$techMean{$protID}-=$minLogXicValue;
		}
		foreach my $protID (keys %techMean) {
			#$proteinData{$protID}{$bioRepLabel}=sprintf '%.0f',$techMean{$protID}; #int(0.5 + 10*$techMean{$protID})/10; #1/10 precision
#my $scaledValue=sprintf '%.0f',(log($techMean{$protID})/$logMod)-$minLogXicValue;
my $scaledValue=$techMean{$protID}-$minLogXicValue;
			$proteinData{$protID}{$bioRepLabel}=($scaledValue < 0)? 0 : $scaledValue;
		}
	}
}

####################################
####>Printing data to table.txt<####
####################################
open(FILESTAT,">>$fileStat") || die "Error while opening $fileStat";
print FILESTAT "2/4 Printing data to file\n";
close FILESTAT;

open(DATA,">$dataDir/table.txt"); # Name of data table
##>Header
print DATA 'Protein';
foreach my $bioRepLabel (@bioRepOrder) {print DATA "\t$bioRepLabel";}
##>Data
if (scalar keys %proteinData) {
	foreach my $protID (sort {&promsMod::sortSmart($a,$b)} keys %proteinData) {
		#next if $bestVisibility{$protID}==0;
		print DATA "\n".$protID;
		foreach my $bioRepLabel (@bioRepOrder) {
			if ($proteinData{$protID}{$bioRepLabel} && $proteinData{$protID}{$bioRepLabel} > 0) {print DATA "\t".(sprintf '%.1f',$proteinData{$protID}{$bioRepLabel});}
			else {print DATA "\t0";}
		}
	}
}
else {warn "ERROR: No proteins quantified!\n";}
close DATA;

&checkForErrors;

$dbh->disconnect;
# close DEBUG;
#die "Data file created!"; # DEBUG



								############################################
			#############################> RUNNING R SCRIPTS <##########################
								############################################
open(FILESTAT,">>$fileStat");
print FILESTAT "3/4 Running SSP Analysis\n";
close FILESTAT;

my %cluster=&promsConfig::getClusterInfo; #('debian') # default is 'centos'
my $pathR=($cluster{'on'})? $cluster{'path'}{'R'} : $promsPath{'R'};

open(R_SCRIPT,">$runDir/analysisCounting.R");
print R_SCRIPT qq
|
###################################################################################
# Launcher for quantification R scripts (provides path for sourcing dependencies) #
###################################################################################

filepath="$promsPath{R_scripts}/"
source(paste(filepath,"analysisCounting.R",sep=""))
|;
close R_SCRIPT;

my $RcommandString="export LANG=en_US.UTF-8; cd $runDir; $pathR/R CMD BATCH --no-save --no-restore analysisCounting.R";

### ALREADY ON CLUSTER (PP 17/09/18)!!! ###
###if ($cluster{'on'}) { ###>Run job on cluster
###	my $numPepValues=(scalar keys %proteinData) * (scalar @bioRepOrder);
###	my $maxHours=int(0.5+($numPepValues/25000)); $maxHours=2 if $maxHours < 2; $maxHours=48 if $maxHours > 48; # 1 M lines -> 40 h
###	my $maxMem=int(1.5 + 1E-6 * $numPepValues);
###	$maxMem.='Gb';
###	my %jobParameters=(
###		commandFile=>"$runDir/runSSProtQuantification.sh",
###		maxMem=>$maxMem,
###		numCPUs=>1,
###		maxHours=>$maxHours,
###		jobName=>"myProMS_SSPA_$quantifID"
###	);
###	my ($pbsError,$pbsErrorFile)=$cluster{'runJob'}->($runDir,$RcommandString,\%jobParameters);
###	if ($pbsError) { # move PBS error message to job error file
###		system "cat $pbsErrorFile >> $promsPath{tmp}/quantification/current/$quantifID\_$quantifDate\_error.txt";
###	}
###}
###else { ###>Run job on Web server
	system $RcommandString;
###}
sleep 3;

&checkForErrors;

$dbh=&promsConfig::dbConnect('no_user'); # reconnect

####>ERROR Management<####
my $RoutStrg=`tail -3 $runDir/analysisCounting.Rout`;
unless ($RoutStrg=~/proc\.time\(\)/) {
	$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID");
	$dbh->commit;
	$dbh->disconnect;
	$RoutStrg=`tail -20 $runDir/analysisCounting.Rout`;
	my $RerrorStrg="R script has generated an error!";
	my $inError=0;
	foreach my $line (split(/\n/,$RoutStrg)) {
		next if (!$inError && $line !~ /^Error in/); # skip everything before "Error in..."
		$inError=1;
		$RerrorStrg.="\n$line";
	}
	die $RerrorStrg;
}

#die "R quantification is finished!!!"; # DEBUG


##################################################
####>Parsing R results & importing data in DB<####
##################################################
open(FILESTAT,">>$fileStat");
print FILESTAT "4/4 Parsing results\n";
close FILESTAT;

####>Fetching list of quantification parameters<####
my %quantifParamIDs;
my $sthQP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,P.CODE FROM QUANTIFICATION_PARAMETER P,QUANTIFICATION_METHOD M WHERE P.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND M.ID_QUANTIFICATION_METHOD=(SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID)");
$sthQP->execute;
while (my ($paramID,$code)=$sthQP->fetchrow_array) {
    $quantifParamIDs{$code}=$paramID;
}
$sthQP->finish;

####>Fetching list of sets used<####
my ($bestDeltaThreshold,$logicThreshold,$pValueThreshold)=@{$quantifParameters{'DB'}{'SPECIF_THRESHOLDS'}}; # default 50,0.05; # hard coded for now
my %usedSets;
my ($allSetsKey,$recAllSetsKey,$allSetsDistKey)=('','','');
foreach my $pos (1..$numStates) {
	$usedSets{"State$pos"}=1; # Always keep single-state sets (required for recording mean state values)
	if ($pos > 1) {
		$recAllSetsKey.='+';
		$allSetsDistKey.=' + ';
	}
	$recAllSetsKey.=$pos;
	$allSetsKey.='State'.$pos;
	$allSetsDistKey.='State'.$pos;
}
open(RES,"$resultDir/outTable.txt"); # reading outTable.txt for the 1st time
while(<RES>) {
	next if $.==1;
	chomp;
	s/"//g;

	my ($modProtID,$setStrg,$effect,$pVal,@stateMeans)=split(/\t/,$_); #"Protein"	"couple"	"effect"	"pvalue"	"State1.mean"	...	"StateN.mean"	"bestDelta"	"model"	"pvalue_adjusted"
	my $pValAdj=pop @stateMeans;
	pop @stateMeans; # remove "model" (algo v2+)
	my $deltaPc=pop @stateMeans; # since analysisCounting v1.1.1
	if ($deltaPc eq 'NA' || ($logicThreshold eq 'or' && $pValAdj > $pValueThreshold && $deltaPc < $bestDeltaThreshold) || ($logicThreshold eq 'and' && ($pValAdj > $pValueThreshold || $deltaPc < $bestDeltaThreshold))) {
		$usedSets{$allSetsKey}=1;
		next;
	}
	# my ($protID,$setStrg)=split(/\t/,$_);
	my ($testStrg,$refStrg)=split(/ \/ /,$setStrg);
	(my $setKey=$testStrg)=~s/ \+ //g;
	$usedSets{$setKey}=1;
}
close RES;

####>Computing all possible sets & positions of used<####
my @sets;
my (@flags) = (0) x $numStates;
$flags[0] = 1; # skip empty subset
my $stateIdx;
while (1) {
	my @subset = ();
	for ($stateIdx=0; $stateIdx<$numStates; $stateIdx++) {
		push @subset, $stateIdx+1 if $flags[$stateIdx] == 1; # record state positions
	}
	last if scalar @subset == $numStates; # exclude complete subset (size=numStates) and exit while loop
	push @sets,\@subset;
	for ($stateIdx=0; $stateIdx<$numStates; $stateIdx++) {
		if($flags[$stateIdx] == 0) {
			$flags[$stateIdx] = 1;
			last;
		}
		$flags[$stateIdx] = 0;
	}
}
my (%setPosition,@recordedSets);
my $setPos=0;
foreach my $refSubset (sort{scalar @{$a} <=> scalar @{$b} || $a->[0]<=>$b->[0]} @sets) { # size then 1st state in subset (should be enough)
	my $setKey='State'.join('State',@{$refSubset});
	next unless $usedSets{$setKey};
	$setPosition{$setKey}=++$setPos;
	my $recSetKey=join('+',@{$refSubset});
	push @recordedSets,$recSetKey;
}
if ($usedSets{$allSetsKey}) {
	$setPosition{$allSetsKey}=++$setPos;
	push @recordedSets,$recAllSetsKey;
}
#my ($setStrg)=($quantifAnnot=~/SETS=([^\$]+)/);
#my %setPosition;
#my $tgtPos=0;
#foreach my $set (split(/;/,$setStrg)) {
#	my $setKey='';
#	foreach my $pos (split('\+',$set)) {
#		$setKey.='State'.$pos;
#	}
#	$setPosition{$setKey}=++$tgtPos;
#}

my $dbhLite=&promsQuantif::dbCreateProteinQuantification($quantifID,$runDir); # SQLite
my $sthInsProt=$dbhLite->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_PROTEIN,SITE_CODE,ID_QUANTIF_PARAMETER,QUANTIF_VALUE,TARGET_POS) VALUES (?,?,?,?,?)"); # PK autoincrement / No site SSPA for now

my %setDistribution;
open(RES,"$resultDir/outTable.txt"); # reading outTable.txt a 2nd time
while(<RES>) {
	next if $.==1;
	chomp;
	s/"//g;
	my ($modProtID,$setStrg,$effect,$pVal,@stateMeans)=split(/\t/,$_); #"Protein"	"couple"	"effect"	"pvalue"	"State1.mean"	...	"StateN.mean"	"bestDelta"	"model"	"pvalue_adjusted"
	my ($protID,$modStrg)=split(/-/,$modProtID);
	my $pValAdj=pop @stateMeans;
	pop @stateMeans; # remove "model" (algo v2+)
	my $deltaPc=pop @stateMeans; # since analysisCounting v1.1.1
	#next if ($deltaPc eq 'NA' || $deltaPc==0); # Skip: no information there
	my ($setKey,$distribSet);
	my ($testStrg,$refStrg)=split(/ \/ /,$setStrg);
	if ($deltaPc eq 'NA' || ($pValAdj > $pValueThreshold && $deltaPc < $bestDeltaThreshold)) {
		$setKey=$allSetsKey;
		$distribSet=$allSetsDistKey;
	}
	else {
		($setKey=$testStrg)=~s/ \+ //g;
		$distribSet=$setStrg;
	}
	# my ($bestState,$refState); #,$testState
	my $bestState;
	foreach my $state (split(/ \+ /,$testStrg)) {
		$bestState=$state if (!$bestState || $stateMeans[$setPosition{$state}-1] > $stateMeans[$setPosition{$bestState}-1]);
		#$testState=$state if (!$testState || $stateMeans[$setPosition{$state}-1] < $stateMeans[$setPosition{$testState}-1]);
	}
	my $bestMean=$stateMeans[$setPosition{$bestState}-1];
next if $bestMean==0; # skip protein (just to be safe when xic)
	$setDistribution{$distribSet}++;
	# foreach my $state (split(/ \+ /,$refStrg)) {
	# 	$refState=$state if (!$refState || $stateMeans[$setPosition{$state}-1] > $stateMeans[$setPosition{$refState}-1]);
	# }

	#>States mean
	foreach my $i (0..$#stateMeans) {
		$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'MEAN'},$stateMeans[$i],$setPosition{'State'.($i+1)}); # states are ordered from table.txt
	}
	#>Effect
	$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'EFFECT'},$effect,$setPosition{$setKey});
	#>Best delta (%) 
	#my $testMean=$stateMeans[$setPosition{$testState}-1];
	#my $refMean=$stateMeans[$setPosition{$refState}-1];
	#my $deltaPc=100*($testMean-$refMean)/$bestMean;
	$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'DELTA_PC'},$deltaPc,$setPosition{$setKey});
	#>p-values
	$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'PVAL'},$pVal,$setPosition{$setKey}) unless $pVal==1;
	$sthInsProt->execute($protID,$modStrg,$quantifParamIDs{'PVAL_ADJ'},$pValAdj,$setPosition{$setKey}) unless $pValAdj==1;
}
close RES;

$sthInsProt->finish;
$dbhLite->commit;
$dbhLite->disconnect;

###>Set distribution
open(DIST,">$resultDir/distribution.txt");
print DIST "Set\tCount\n";
foreach my $setStrg (sort{length($a)<=>length($b) || &promsMod::sortSmart($a,$b)} keys %setDistribution) {
	print DIST "$setStrg\t$setDistribution{$setStrg}\n";
}
close DIST;

&checkForErrors;

###> Quantification is finished.
my $status=($quantifParameters{'DB'}{'VISIBILITY'})? $quantifParameters{'DB'}{'VISIBILITY'}[0] : 1;
my $setString=join(';',@recordedSets);
$dbh->do("UPDATE QUANTIFICATION SET STATUS=$status,QUANTIF_ANNOT=CONCAT(QUANTIF_ANNOT,'::SETS=$setString') WHERE ID_QUANTIFICATION=$quantifID");
$dbh->commit;

$dbh->disconnect;
#exit; # DEBUG!!!

open(FILESTAT,">>$fileStat");
print FILESTAT "Ended ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
close(FILESTAT);

####>Move all files outside tmp directory
unlink "$runDir/analysisCounting.R";

mkdir "$promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
dirmove($runDir,"$promsPath{quantification}/project_$projectID/quanti_$quantifID");

sleep 2;



sub createModProtID { # GLOBALS: $quantifiedModifID, $keepSeqContext, %modifQuantifRank, $ambiguousModifPos, %quantifResMatchStrg
	my ($protID,$pepSeq,$vmod,$pepSeqBeg,$refModifQuantifs)=@_;

	my $pepSeqContext=($keepSeqContext)? '['.$pepSeqBeg.'.'.length($pepSeq).']' : '';
	my @seq=split(//,$pepSeq);
	
	my $nonAmbigModProtID=$protID; # default
	my $ambigModProtID=$protID; # default
	my $seqIsAmbig=0;
	
	my $firstModQ=1;
	foreach my $modID (@{$refModifQuantifs}) {
		my ($sep,$posStrg)=($vmod=~/(^|&)$modID:([\d+\.\-\=\+\*]+)/); # added N/C-term
		next unless $posStrg; # in multi-modif quantif context: varmod does not necessarly contains all quantified modif
		#>Non-ambiguous pos
		my @pos=(); # reset
		foreach my $pepPos (split(/\./,$posStrg)) {
			if ($pepPos=~/\d+/) {push @pos,$seq[$pepPos-1].($pepSeqBeg+$pepPos-1);}
			#elsif ($pepPos=~/[=\*]/) {next;} # skip peptide N/C-term: cannot be in protein TODO: add to table.txt for data normalisation but do not quantify
			elsif ($pepPos eq '=') {push @pos,'n'.$pepSeqBeg;} # peptide N-term: kept for normalization purpose
			elsif ($pepPos eq '-') {push @pos,'n0';} # protein N-term: kept for normalization purpose
			elsif ($pepPos eq '*') {push @pos,'c'.($pepSeqBeg+$#seq);} # peptide C-term: kept for normalization purpose
			elsif ($pepPos eq '+') {push @pos,'c0';} # protein C-term: kept for normalization purpose
		}
		#>THIS IS WHERE THE PROTEIN & MODIFICATION SITE(S) ARE CODED!!!!!!
		# Single-modif quantif:
		#	-Non-ambiguous: "ProtID-Res1Pos1[.Res2Pos2...]" eg "1234-S56.T78"
		#	-Ambiguous: "ProtID-StartRes~EndRes:NumModifs/NumAcceptors" eg "1234-56~78:2/5"
		# Mutli-modif quantif:
		#	-Non-ambiguous: "ProtID-ModifXRank#Res1Pos1[.Res2Pos2...][&ModifYRank#Res3Pos3[.Res4Pos4...]]" eg "1234-1#S56.T78&2#C61" sorted by increasing rank!!!
		#	-Ambiguous: "ProtID-ModifXRank#StartRes~EndRes:NumModifs/NumAcceptors" eg "1234-1#56~78:2/5&2#C61" (PP: Not allowed from start form as of 30/06/20)
		my $nonAmbigMod=($quantifiedModifID)? join('.',@pos) : $modifQuantifRank{$modID}.'#'.join('.',@pos); # Add '<modif rank>#' if multi-modif quantif
		$nonAmbigModProtID.=($firstModQ)? '-' : '&';
		$nonAmbigModProtID.=$nonAmbigMod;
		
		# ==> Look for ambiguity
		if ($ambiguousModifPos) {
			$ambigModProtID.=($firstModQ)? '-' : '&';
			$ambigModProtID.=$modifQuantifRank{$modID}.'#' unless $quantifiedModifID;
			#my $numModifRes = () = $posStrg =~ /(\d+)/g;
			my $numModifRes=scalar @pos;
			my @targetResIdx;
			while ($pepSeq =~ /[$qQuantifResMatchStrg{$modID}]/g) {push @targetResIdx,$-[0];} # position of match (1st=0)
			my $numTargetRes=scalar @targetResIdx;
			my ($NtermMatch,$CtermMatch)=('','');
			if ($quantifResMatchStrg{$modID}=~/=/) {$NtermMatch='n'.$pepSeqBeg;} # peptide N-term
			elsif ($quantifResMatchStrg{$modID}=~/\*/) {$CtermMatch='c'.($pepSeqBeg+$#seq);} # peptide C-term
			elsif ($quantifResMatchStrg{$modID}=~/-/) {$NtermMatch='n0' if $pepSeqBeg==1;} # protein N-term
			#elsif ($quantifResMatchStrg=~/\+/) {$CtermMatch='c0' if $pepSeqBeg+length($pepSeq)-1==<protein length>;} # TODO: fetch prot length
			my $numAlltargets=$numTargetRes;
			$numAlltargets++ if $NtermMatch;
			$numAlltargets++ if $CtermMatch;
			if ($numModifRes < $numAlltargets) { # Ambiguity! protID-[modifPos=]<1stPos>~<lastPos>:<numModif>/<numSites>
				$seqIsAmbig=1;
				#$ambigModProtID=$protID.'-';
				if ($NtermMatch) {$ambigModProtID.=$NtermMatch;} # N-term
				else {$ambigModProtID.=($pepSeqBeg+$targetResIdx[0]);} # internal pos
				$ambigModProtID.='~';
				if ($CtermMatch) {$ambigModProtID.=$CtermMatch;} # C-term
				else {$ambigModProtID.=($pepSeqBeg+$targetResIdx[-1]);} # internal pos
				$ambigModProtID.=':'.$numModifRes.'/'.$numAlltargets;
			}
			else {$ambigModProtID.=$nonAmbigMod;}
		}
		$firstModQ=0;
	}
	my $modProtID=($seqIsAmbig)? $ambigModProtID : $nonAmbigModProtID;

	return $modProtID.$pepSeqContext;
}



sub checkForErrors {
	#my ($errorFile) = glob "$promsPath{tmp}/quantification/current/*_$quantifDate\_error.txt";
	my $errorFile = "$promsPath{tmp}/quantification/current/$quantifID\_$quantifDate\_error.txt";
	if ($errorFile && -s $errorFile) {
		$dbh->do("UPDATE QUANTIFICATION SET STATUS=-2 WHERE ID_QUANTIFICATION=$quantifID"); # Failed
		$dbh->commit;
		$dbh->disconnect;
		die "Aborting quantification due to errors.";
	}
}


####>Revision history<####
# 1.8.0 [UPDATE] New "None" category for proteins/sites not specific to any (set of) states (PP 04/05/21)
# 1.7.2 [UPDATE] SSPA code version management moved here (PP 11/04/21) 
# 1.7.1 [BUGFIX] Fix no data in PTM-enriched proteins when all PTMs are excluded (PP 15/03/21)
# 1.7.0 [CHANGE] Compatible with PTM-enriched protein quantification & added forgotten protein filtering (PP 09/03/21)
# 1.6.2 [CHANGE] Compatible with Free residues & down shift of log values for xic (PP 03/12/20)
# 1.6.1 [BUGFIX] In Label modification filtering (PP 13/10/20)
# 1.6.0 [FEATURE] Compatible with sites (PP 13/10/20)
# 1.5.0 [UPDATE] Uses SQLite database to store quantification data & read best delta from outTable.txt (PP 12/08/20)
# 1.4.0 [ENHANCEMENT] New log2-based rescaling for xic feature & set distribution file (PP 29/03/20)
# 1.3.4 [BUGFIX] In DIA-based xic feature leading to empty table.txt (PP 25/03/20)
# 1.3.3 [FEATURE] Compatible with visible/hidden proteins, peptide xic as feature and only used sets are recorded (PP 28/01/20)
# 1.3.2 [FEATURE] Compatible with MBR-peptide filter (PP 22/01/20)
# 1.3.1 [BUGFIX] commented forgotten exit & [ENHANCEMENT] compatible with quantification visibility status (PP 19/12/19)
# 1.3.0 [FEATURE] Added optional co-exclusion of fully cut peptides included in missed-cut ones & No cluster call for R since caller script is already running on cluster (PP 19/11/19)
# 1.2.4 Do not remove status log file to match new monitoring system behavior (VS 10/10/19)
# 1.2.3 Uses single $cluster{runJob}->() command (PP 11/07/19)
# 1.2.2 Uses quanti_$quantifID instead of quantif_$quantifID (PP 18/05/18)
# 1.2.1 [FIX] Bug due to truncated protein:visibility SQL response string if &gt; 1024 characters (PP 18/01/18)
# 1.2.0 Now relies on &promsConfig::getClusterInfo (PP 30/11/17)
# 1.1.0 Handle multiple peptide counting methods including spectral count (PP 02/08/17)
# 1.0.0 Started (PP 04/08/16)
