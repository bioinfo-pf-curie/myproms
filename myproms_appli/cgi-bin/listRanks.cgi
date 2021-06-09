#!/usr/local/bin/perl -w

################################################################################
# listRanks.cgi             2.3.11                                             #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Lists the peptides from a protein or query being validated                   #
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
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use strict;
use phosphoRS;
use URI::Escape;#
#use String::Unquotemeta;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my $allRank=10;
my $DELTA_MASS=0.25; # delta mass allowed when matching a peptide in multiple analyses (Mexp => for charge = +1)
my $DELTA_RANGE=2*$DELTA_MASS;
my %convertPos2Text=('-'=>'Protein N-term','='=>'Any N-term','+'=>'Protein C-term','*'=>'Any C-term');

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

####################
####>Parameters<####
####################
my $analysisID=param('ANAID');
#my $msType=param('MSTYPE');
my ($minScore,$maxRank,$fileFormat,$msType)=$dbh->selectrow_array("SELECT MIN_SCORE,MAX_RANK,FILE_FORMAT,MS_TYPE FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
#$msType='MIS' unless $msType eq 'PMF';
unless ($minScore) {
	$minScore=&promsConfig::getMinScore($fileFormat);
}
$maxRank=($msType ne 'MIS')? 10 : ($maxRank)? $maxRank : &promsConfig::getMaxRank; # <--- Maximum number of ranks retrieved (Absolute max = 10)

my $item=param('ITEM'); # query, protein or myList
my $itemID=param('ID'); # ID_PROT_VALID or ID_QUERY
my ($page,$maxPage)=split('-',param('PAGE'));
my $rangeStrg=(param('RANGE'))? param('RANGE') : 0; # range of proteins (match groups) or queries (id1-idn). 0 for myList
my $selectedSort=param('SORT'); # sort for proteins/queries!
my $maxRankDisplay=(param('maxDisp'))? param('maxDisp') : $maxRank;
my $maxRankSearch=(param('maxSrch'))? param('maxSrch') : $maxRank;
my $multiAna=param('multiAna') || 0;
my $showFiltered=param('SHOWFILT');
&validateMatchGroup if param('valMG');
&processForm if param('save'); # lower case means form submission !
&changeProteinStatus if defined(param('STATUS')); # item = protein, status can be 0
my $projectID=&promsMod::getProjectID($dbh,$analysisID,'analysis');
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
# my $disableLink=($projectAccess=~/power/)? 'disabled' : ''; # all others have full access
# my $disableForm=$disableLink; # Link:interpretations, Form: validations
my $disableForm=($projectAccess =~ /(bioinfo|mass|manag|super)/)? '' : 'disabled'; # have full access
my $disableLink=$disableForm; # disableForm can be changed later if protein selSatus=-2
my $disableExclude=($disableForm && $projectAccess !~ /power/)? 'disabled' : '';
my $rankSort=(param('sortType'))? param('sortType') : 'score'; # sort for peptides!
$rankSort='query' if ($msType eq 'PMF' && $rankSort eq 'score');
my $activatedPeptide=(param('actPep'))? param('actPep') : ''; # make peptide selectable and return to main
my $DeltaType = (param('Dlt'))?param('Dlt'):'ppm';
my $showPPM=($DeltaType eq 'ppm')?1:0;
my $ptmFilter=(param('varValue'))? param('varValue') : 'Indifferent';
(my $fixRes)=($ptmFilter =~ s/^(.):://g) ? $1 : '';
my $oldPtmFilter=$ptmFilter;
my ($pepRangeMin,$pepRangeMax); # Protein range for listing peptides
if ($item eq 'protein' || $item eq 'myList') {
	($pepRangeMin,$pepRangeMax)=param('pepRange')? split(/:/,param('pepRange')) : (1,0); # prot length not yet known
}

#############################################
####>Retrieval of variable modifications<####
#############################################
my @modifications=('Indifferent','None','Any');
#my %infoSubmit=promsMod::getSearchParam($dbh,$analysisID);
#push @modifications,sort(split(/,/,$infoSubmit{'g:Variable modifications'})) if $infoSubmit{'g:Variable modifications'};
my %ptmList;
#my $sthgetAMV=$dbh->prepare("SELECT AM.ID_MODIFICATION,AM.SPECIFICITY,PSI_MS_NAME,INTERIM_NAME,SYNONYMES FROM ANALYSIS_MODIFICATION AM, MODIFICATION M WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND MODIF_TYPE='V' AND ID_ANALYSIS=$analysisID");
my $sthgetAMV=$dbh->prepare("SELECT AM.ID_MODIFICATION,AM.SPECIFICITY,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,MODIF_TYPE FROM ANALYSIS_MODIFICATION AM, MODIFICATION M WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND ID_ANALYSIS=$analysisID");
$sthgetAMV->execute;
while (my ($modID,$specificity,$psiName,$interName,$altNames,$modifType)=$sthgetAMV->fetchrow_array) { # Get only ANALYSIS_MODIFICATION Variable (not interested in Fix ones)
	#my $modName=($psiName)?$psiName:($interName)?$interName:"mod_$modID";
	$altNames=~ s/##/,/g if $altNames;
	$altNames=~ s/^,//; $altNames=~ s/,\Z//; # remove starting & trailing "," if any
	my $modName=($psiName)?$psiName:($interName)?$interName:($altNames)? $altNames :"mod_$modID";
	next if $modName eq 'MappingL';
	$ptmList{$modID}{'SPECIFICITY'}=$specificity;
	$ptmList{$modID}{'NAME'}=$modName;
	$ptmList{$modID}{'MODIF_TYPE'}=$modifType;
}
$sthgetAMV->finish;

############################################
####>Fetching data for selected protein<####
############################################
my ($identifier,$dbRank,$weight,$protDes,$protOrg,$score,$confLevel,$numDatabanks);
my $selStatus=0; # initialize to something in case item=query
my %pepList; # keys will be queryIDs (protein view) or queryNum (query view)
my %boundaryAA; # stores list of flanking Nter/Cter residues for protein matches
my %pepMatchInfo; # store matching info so as to apply the range filter afterwards
if ($item eq 'protein' || $item eq 'myList') {
	($numDatabanks)=$dbh->selectrow_array("SELECT COUNT(*) FROM ANALYSIS_DATABANK WHERE ID_ANALYSIS=$analysisID");
	($identifier,$dbRank,$weight,$protDes,$protOrg,$selStatus,$score,$confLevel,my $protLength)=$dbh->selectrow_array("SELECT IDENTIFIER,DB_RANK,MW,PROT_DES,ORGANISM,SEL_STATUS,SCORE,CONF_LEVEL,PROT_LENGTH FROM PROTEIN_VALIDATION WHERE ID_PROT_VALID=$itemID");
	unless ($identifier) { # no match => protein was deleted from DB (due to user rank update)
		print header(-charset=>'UTF-8'),"\n";
		print qq
|<HTML>
<BODY bgcolor="white" background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR><BR>
<FONT class="title2">Currently selected protein must have deleted. Select a new one from list.</FONT>
</CENTER>
</BODY>
</HTML>
|;
		$dbh->disconnect;
		exit;
	}

	###>Matching peptides<###
 	my $sth=$dbh->prepare("SELECT QUERY_NUM,PEP_RANK,MATCH_INFO FROM RANK_PROTEIN_MATCH WHERE IDENTIFIER='$identifier' AND ID_ANALYSIS=$analysisID");
	$sth->execute || die $sth->errstr;
	my $refData=$sth->fetchall_arrayref; # reference to an array
	$sth->finish;

	###>Peptide range filter<###
	$protLength=0 unless $protLength;
	$pepRangeMin = 1 if $pepRangeMin < 1;
	if ($protLength) {
		$pepRangeMax = $protLength if $pepRangeMax > $protLength;
		($pepRangeMin,$pepRangeMax) = (1,$protLength) if $pepRangeMin > $pepRangeMax;
	}

	###>Storing matchInfo in %pepList<###
	my $maxBeg=0;
	foreach my $refMatch (@{$refData}) { # ascending identifier
		my ($queryNum,$rank,$matchInfo)=@{$refMatch};

		next if $rank > $maxRankDisplay; # rank display filter

		##>Checking if peptide is in peptide Range<##
		if ($pepRangeMax) {
			my $isInRange=0;
			foreach my $match (split(/:/,$matchInfo)) {
				my ($beg,$aa1,$aa2)=split(/,/,$match);
				if ($beg >= $pepRangeMin && $beg <= $pepRangeMax) { # only beg must be in range!
					$isInRange=1;
					last;
				}
			}
			next unless $isInRange;
		}
		else {
			foreach my $match (split(/:/,$matchInfo)) {
				my ($beg,$aa1,$aa2)=split(/,/,$match);
				$maxBeg=$beg if $maxBeg<$beg;
			}
		}

		$pepMatchInfo{"$queryNum:$rank"}=$matchInfo;
		push @{$pepList{$queryNum}},$rank; # <= queryNum != queryID !!!
		my $firstMatch=(split(/:/,$matchInfo))[0];
		my ($beg,$aa1,$aa2)=split(/,/,$firstMatch);
		if ($aa1) {$aa1.='.';} else {$aa1='';}
		if ($aa2) {$aa2=".$aa2";} else {$aa2='';}
		@{$boundaryAA{"$queryNum:$rank"}}=&promsMod::chkDef($aa1,$aa2);
	}
	$pepRangeMax=$maxBeg unless $pepRangeMax;
	$disableForm='disabled' if $selStatus<=-2;
}
else { # item = query
	$identifier=''; # needs to be defined for javascript
	foreach my $rank (1..$maxRankDisplay) { #(my $rank=1;$rank<=$maxRank;$rank++) {
		push @{$pepList{$itemID}},$rank; # <= queryID != queryNum !!!
		@{$boundaryAA{"$itemID:$rank"}}=('','');
	}
}

#######################
####>Starting HTML<####
#######################
print header(-'content-encoding'=>'no',-charset=>'UTF-8'),"\n";
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Rank List</TITLE>
<LINK rel=stylesheet href="$promsPath{html}/promsStyle.css" type=text/css>
<STYLE type="text/css">
TD {font-weight:bold;text-align:center;}
TD.left {text-align:left;}
TD.right {text-align:right;}
</STYLE>
<SCRIPT language="JavaScript">
|;
&promsMod::popupInfo();
print qq
|</SCRIPT>
</HEAD>
<BODY topmargin=2px bottommargin=0px leftmargin=5px background="$promsPath{images}/bgProMS.gif">
<DIV id="waitDiv" style="position:absolute;top:40%;width:100%">
<TABLE align=center>
<TR><TH>Fetching data. Please, wait<SPAN id="waitSpan">...</SPAN></TH></TR>
<TR><TH><IMG src="$promsPath{images}/scrollbarGreen.gif"/></TH></TR>
</TABLE>
</DIV>
|;

################################################
####>Fetching data for query+peptide(s) set<####
################################################
my (%queryNum,%queryStatus,%queryCharge,%queryPtmProb,%pepInfo,%massExp,%massObs,%elutionTime,%varModNames);
my %rankPos; # idem as %pepList but keys are always queryIDs

##>Activating peptide?<##
my ($actQID,$actQNum,$actRank)=($activatedPeptide)? split(/_/,$activatedPeptide) : (0,0,0);

my $dbQuery='SELECT QV.ID_QUERY,QUERY_NUM,VALID_STATUS,MASS_DATA,CHARGE,ELUTION_TIME,GROUP_CONCAT(CONCAT(QM.ID_MODIFICATION, "=", SUBSTRING_INDEX(QM.REF_POS_STRING, "##", -1)) SEPARATOR "&&")';
foreach my $rank (1..$maxRankDisplay) {$dbQuery.=",INFO_PEP$rank";} # rank starts at 1 not 0 anymore!
$dbQuery.=' FROM QUERY_VALIDATION QV LEFT JOIN QUERY_MODIFICATION QM ON QM.ID_QUERY=QV.ID_QUERY WHERE ';
if ($item eq 'query') {$dbQuery.='QV.ID_QUERY=?';}
else {$dbQuery.="QUERY_NUM=? AND ID_ANALYSIS=$analysisID";}
$dbQuery.=" GROUP BY QV.ID_QUERY";
my $sthSelQ=$dbh->prepare("$dbQuery");
foreach my $query (keys %pepList) { # qID or qNUM
	$sthSelQ->execute($query);
	my ($queryID,$qNum,$valid,$massData,$charge,$elutTime,$ptmProb,@pepInfoList)=$sthSelQ->fetchrow_array;
	$queryNum{$queryID}=$qNum;
	$queryStatus{$queryID}=$valid;
	$queryCharge{$queryID}=$charge;
	
	if ($ptmProb && $ptmProb =~ /PRB_/) { # non-PRS position probability
		$ptmProb =~ s/PRB_[^=]+=//;
		foreach my $modifInfos (split(/\&\&/, $ptmProb)) {
			my ($modifID, $ptms) = split(/=/, $modifInfos);
			next unless($modifID && $ptms);
			if(!defined($varModNames{$modifID})) { # For PTM not retrieved from promsMod::getVariableModifications (i.e. not valid dispcode and dispcolor)
				my ($psiMsName,$interimName,$altNames)=$dbh->selectrow_array("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES FROM MODIFICATION WHERE ID_MODIFICATION=$modifID");
				my $name=($psiMsName)?$psiMsName:($interimName)?$interimName:($altNames)?$altNames:'';
				$varModNames{$modifID}=$name;
			}
			$queryPtmProb{$queryID}{$modifID} = $ptms;
		}
	}
	($massExp{$queryID})=($massData=~/EXP=(\d+\.*\d*)/);
	($massObs{$queryID})=($massData=~/OBS=(\d+\.*\d*)/);
	if ($elutTime) {
		if ($elutTime=~/^sp/) {$elutionTime{$queryID}=$elutTime;}
		else {
			$elutTime=($elutTime=~/et(\d+\.\d+);/)? $1 : $elutTime;
			$elutTime=~s/^et//; # remove et tag if any
			$elutionTime{$queryID}="$elutTime min.";
		}
	}
	else {$elutionTime{$queryID}='?';}
	my %curPepInfo;
	my $i=0;
	foreach my $rank (1..$maxRankDisplay) {
last unless $pepInfoList[$i]; # no more interpretations
		$curPepInfo{$rank}=$pepInfoList[$i];
		$i++;
	}
	foreach my $rank (@{$pepList{$query}}) {
		push @{$rankPos{$queryID}},$rank; #(=@{$pepList{$query}}) keeps track of ranks (may not be continuous if protein view)
		##>Activating peptide
		if ($actQNum==$qNum && $actRank==$rank) {
			$queryStatus{$queryID}=-1 if $queryStatus{$queryID}==-3; # cannot be -4
			$curPepInfo{$actRank}=~s/SEL=-1/SEL=0/;
			$dbh->do("UPDATE QUERY_VALIDATION SET VALID_STATUS=$queryStatus{$queryID},INFO_PEP$rank='$curPepInfo{$rank}' WHERE ID_QUERY=$actQID");
			#>Updating matched proteins
			my $sthMax=$dbh->prepare("SELECT ID_PROT_VALID,MAX_MATCH,MAX_SCORE FROM PROTEIN_VALIDATION PV,RANK_PROTEIN_MATCH RPM WHERE PV.IDENTIFIER=RPM.IDENTIFIER AND PV.ID_ANALYSIS=$analysisID AND RPM.ID_ANALYSIS=$analysisID AND QUERY_NUM=$actQNum AND PEP_RANK=$actRank");
			my $sthUpMax=$dbh->prepare("UPDATE PROTEIN_VALIDATION SET MAX_MATCH=?,MAX_SCORE=? WHERE ID_PROT_VALID=?");
			$sthMax->execute;
			while (my ($protValidID,$maxMatch,$maxScore)=$sthMax->fetchrow_array) {
				$maxMatch++; # match frequency is not considered
				my ($actScore)=($curPepInfo{$actRank}=~/SC=(\-?\d+\.?\d*)/);
				$maxScore+=$actScore;
				$sthUpMax->execute($maxMatch,$maxScore,$protValidID);
			}
			$sthMax->finish;
			$sthUpMax->finish;
			$dbh->do("UPDATE ANALYSIS SET LOWER_SCORES=1 WHERE ID_ANALYSIS=$analysisID AND LOWER_SCORES !=2");

			promsMod::updateAnalysisHistory($dbh, $analysisID, undef, 'lowP_m');

			$dbh->commit;
		}
		push @{$pepInfo{$queryID}},$curPepInfo{$rank};
	}
}
$sthSelQ->finish;
&updateWaitBox;

###############################################################################
####>Finding number of other analyses matching peptides + Recording Scores<####
###############################################################################
###>Fetching list of sequences in current Analysis
my (%peptideList,%sequenceList,%varModList,%seqMassObs,%scoreList,%absDeltaList,%commentList,%existBetterPep);
#my (%absPpmList) ;
my $noSelect;
foreach my $queryID (keys %queryNum) {
	foreach my $p (0..$#{$pepInfo{$queryID}}) { # p is index not rank
		last unless $pepInfo{$queryID}[$p];
		my ($sequence)=($pepInfo{$queryID}[$p]=~/SEQ=(\w+)/);
		$sequenceList{$sequence}=1; # list of bare sequences
		($scoreList{"$queryID:$p"})=($pepInfo{$queryID}[$p]=~/SC=(\-?\d+\.?\d*)/);
		my ($deltaDa)=($pepInfo{$queryID}[$p]=~/DELT=([^,]+)/);
		$absDeltaList{"$queryID:$p"}=abs($deltaDa);
		#if ($showPPM) {
		#	my ($massCalc)=($pepInfo{$queryID}[$p]=~/CALC=(\d+\.*\d*)/);
		#	$absPpmList{"$queryID:$p"}=1000000*abs($deltaDa/$massCalc);
		#}
		($commentList{"$queryID:$p"})=($pepInfo{$queryID}[$p]=~/COM=(.+)/); # COM must be at the end!
		$commentList{"$queryID:$p"}='' unless $commentList{"$queryID:$p"};
		my $varMod=&promsMod::toStringVariableModifications($dbh,'QUERY',$queryID,$analysisID,$sequence,$rankPos{$queryID}[$p]); # compute for all peptides!
		#$varMod =~ s/Phospho\s\(.+?(\d+)/Phospho \($1/;
		$varModList{"$queryID:$p"}=$varMod;
		my ($sel)=($pepInfo{$queryID}[$p]=~/SEL=(-?\d)/);
		if ($sel==-1) { # better score exists
			$noSelect=1; # exist not selectable ranks
			next;
		}
		#keeping best score, unless undefined
		if (!defined($peptideList{"$sequence$varMod"}) ) {
			$peptideList{"$sequence$varMod"}="$queryID:$p";
		}
		elsif ($scoreList{"$queryID:$p"} > $scoreList{$peptideList{"$sequence$varMod"}}) { #this one in the current best
			$existBetterPep{$peptideList{"$sequence$varMod"}}="$queryID:$p" ; #link to better
			$peptideList{"$sequence$varMod"}="$queryID:$p" ;	#reattribute
		}
		else {
			$existBetterPep{"$queryID:$p"}=$peptideList{"$sequence$varMod"}; #link to current best
		}
	}
}
&updateWaitBox;

###>Fetching list of sequences in all other MS/MS Analyses
my (%analysisMatch,%extendAnalysisMatch);
if ($multiAna && $msType ne 'PMF') { # not if PMF
	&updateWaitBox('<BR>Scanning all analyses for matching interpretations...');
	my $sthAllAna=$dbh->prepare("SELECT ID_ANALYSIS,VALID_STATUS,MAX_RANK FROM ANALYSIS,SAMPLE,EXPERIMENT WHERE MS_TYPE !='PMF' AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_EXPERIMENT=EXPERIMENT.ID_EXPERIMENT AND ID_PROJECT=$projectID");
	$sthAllAna->execute();
	my $refAnaList=$sthAllAna->fetchall_arrayref();

	##>Query for non validated analyses
	my %sthNoValPep;
	foreach my $Rank (1..$maxRankSearch) {
		my $infoString='ID_QUERY,MASS_DATA,CHARGE';
		foreach my $rank (1..$Rank) {$infoString.=",INFO_PEP$rank";}
		$sthNoValPep{$Rank}=$dbh->prepare("SELECT $infoString FROM QUERY_VALIDATION WHERE ID_ANALYSIS=? AND VALID_STATUS >=-2 AND QUERY_NUM > 0");
	}

	##>Query for validated analyses
	my $sthValPep=$dbh->prepare("SELECT ID_PEPTIDE,MR_EXP,CHARGE FROM PEPTIDE WHERE ID_ANALYSIS=? AND PEP_SEQ=? AND SCORE>0"); #,VAR_MOD   SCORE>0 in case PMF within mixed search (SQ) or Ghosts
	foreach my $refAna (@{$refAnaList}) {
		my ($anaID,$validStatus,$otherAnaMaxrank)=@{$refAna};
		next if $anaID==$analysisID;
		$otherAnaMaxrank = $maxRankSearch if (!$otherAnaMaxrank || $otherAnaMaxrank > $maxRankSearch);

		##>Non (partially) validated Analyses
		if ($validStatus<=1) {
			$sthNoValPep{$otherAnaMaxrank}->execute($anaID);
			while (my ($qID,$massData,$charge,@infoData)=$sthNoValPep{$otherAnaMaxrank}->fetchrow_array) {
				my $rank=0;
				#my ($massOBS)=($massData=~/OBS=(\d+\.*\d*)/);
				my ($massEXP)=($massData=~/EXP=(\d+\.*\d*)/);

				#>loop on all peptides displayed
				my @massOKPeptides;
				foreach my $pep (keys %peptideList) { # pep=seq+vModbis
					my ($queryID,$p)=split(/:/,$peptideList{$pep});
					push @massOKPeptides,$pep if abs($massEXP-$massExp{$queryID}) < $DELTA_RANGE; # out of acceptable delta range
				}
				foreach my $info (@infoData) {
					$rank++;
					last unless $info;
					next if $info=~/SEL=-1/; # better score exists !!!!!!!!!!!!!!!!!!! may be should be included too
					next if $info=~/SC=0,/; # no score => PMF within mixed search (SQ)
					my ($seq)=($info=~/SEQ=(\w+)/);
					next unless $sequenceList{$seq}; # seq is not in displayed list
					#my $vMod=&promsMod::toStringVariableModifications($dbh,'QUERY',$qID,$anaID,$seq,$rank);
					#$vMod=~s/Phospho\s\(.+?(\d+)/Phospho \($1/;
					my $pep2=$seq.&promsMod::toStringVariableModifications($dbh,'QUERY',$qID,$anaID,$seq,$rank); # vMod

					#>loop on all mass OK peptides displayed
					foreach my $pep (@massOKPeptides) { # pep=seq+vModbis
						next unless $pep eq $pep2;
						my ($queryID,$p)=split(/:/,$peptideList{$pep});
						push @{$analysisMatch{$peptideList{$pep}}},$validStatus."_$qID"."_$rank" if ($charge && $charge==$queryCharge{$queryID});
						push @{$extendAnalysisMatch{$peptideList{$pep}}},$validStatus."_$qID"."_$rank";# any charge state (if not elsif => sum)
					}
				}
			}
		}

		##>Validated Analyses
		else {
			foreach my $pep (keys %peptideList) {
				my ($queryID,$p)=split(/:/,$peptideList{$pep});
				my ($seq,$vMod)=($pep=~/(\w+)(.*)/); # pep=seq + vMod
				$sthValPep->execute($anaID,$seq);
				while (my ($pepID,$massEXP,$charge)=$sthValPep->fetchrow_array) {
					next unless abs($massEXP-$massExp{$queryID}) < $DELTA_RANGE;
					my $varMod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$anaID,$seq);
					#$varMod=~s/Phospho\s\(.+?(\d+)/Phospho \($1/; # bis
					next unless $varMod eq $vMod; # positions must be the same (PP 24/04/12)
					push @{$analysisMatch{$peptideList{$pep}}},"2_$pepID" if $charge==$queryCharge{$queryID};
					push @{$extendAnalysisMatch{$peptideList{$pep}}},"2_$pepID"; # any charge state (if not elsif => sum)
				}
			}
		}
	}
	$sthAllAna->finish;
	for my $Rank (1..$maxRankSearch) {
		$sthNoValPep{$Rank}->finish;
	}
	$sthValPep->finish;
# 	$sthValPepNull->finish;
}

#Looking for commun sequence in current analysis (in case lower score selectable)
foreach my $pep2 (keys %scoreList) {
	my ($queryID,$p)=split(/:/,$pep2); # p is index!
	next unless $pepInfo{$queryID}[$p];
	my ($sel)=($pepInfo{$queryID}[$p]=~/SEL=(-*\d)/);
	next if $sel==-1;
	my ($seq)=($pepInfo{$queryID}[$p]=~/SEQ=(\w+)/);
	next unless $seq;
	$seq.=$varModList{"$queryID:$p"};
# Added 15/12/2010 in 2.1.2 by ppoullet (idem as old version) -->
	#my $vMod=&promsMod::toStringVariableModifications($dbh,'QUERY',$queryID,$analysisID,$seq,$rankPos{$queryID}[$p]);
	#$seq.=$vMod;
	##$varMod=join(' + ',sort (split(/ \+ /,$varMod))); # not needed: same analysis
	##$vMod =~ s/Phospho\s\(.+?(\d+)/Phospho \($1/; # not needed: same analysis
# <--- end of addition

	#Old version
	#my ($vMod)=($info=~/VMOD=([^,]+)/);
	#$vMod='' unless $vMod;
	#my $varModbis = $vMod;
	#$varModbis =~ s/Phospho\s\(.+?(\d+)/Phospho \($1/;
	#if (defined($peptideList{"$seq$varModbis"}) && $peptideList{"$seq$varModbis"} ne $pep2) { # a sequence+mod is matchedin this analyses
	#	my ($qID,$rank)=split(/:/,$peptideList{"$seq$varModbis"});
	#	$p++ ;
	#	if ($massObs{$qID}>$massObs{$queryID}-1 && $massObs{$qID}<$massObs{$queryID}+1) { #same charge state
	#		push @{$analysisMatch{$peptideList{"$seq$varModbis"}}},$sel."_$queryID"."_$p";
	#	}
	#	push @{$extendAnalysisMatch{$peptideList{"$seq$varModbis"}}},$sel."_$queryID"."_$p";
	#}
	if (defined($peptideList{$seq}) && $peptideList{$seq} ne $pep2) { # a sequence is matched in this analyses
		my ($qID,$p0)=split(/:/,$peptideList{$seq});
		#$p++;
		#if (abs(($massObs{$qID}*$queryCharge{$qID})-($massObs{$queryID}*$queryCharge{$queryID})) < $DELTA_RANGE) { #same charge state, delta allowed =+/-DELTA_MASS
		#	push @{$analysisMatch{$peptideList{"$seq"}}},$sel."_$queryID"."_$p";
		#}
		#push @{$extendAnalysisMatch{$peptideList{"$seq"}}},$sel."_$queryID"."_$p";

		if (abs($massExp{$qID}-$massExp{$queryID}) < $DELTA_RANGE) { # delta allowed =+/-DELTA_MASS
			$p++; # from index to rank
			push @{$analysisMatch{$peptideList{"$seq"}}},$sel."_$queryID"."_$p" if $queryCharge{$qID}==$queryCharge{$queryID}; # same charge state
			push @{$extendAnalysisMatch{$peptideList{"$seq"}}},$sel."_$queryID"."_$p";
		}
	}
}
&updateWaitBox;

########################################################################
####>Finding number of non-filtered proteins matching each peptides<####
########################################################################
my %numSelProt;
if (!$showFiltered) {
	my %identList;
	my $sthPV=$dbh->prepare("SELECT IDENTIFIER FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$analysisID AND IDENTIFIER NOT LIKE 'DECOY_%' AND SEL_STATUS > -3");
	$sthPV->execute;
	while (my ($identifier)=$sthPV->fetchrow_array) {$identList{$identifier}=1;}
	$sthPV->finish;

	my $sthMP=$dbh->prepare("SELECT IDENTIFIER FROM RANK_PROTEIN_MATCH WHERE QUERY_NUM=? AND PEP_RANK=? AND ID_ANALYSIS=$analysisID");
	foreach my $queryID (keys %queryNum) {
		foreach my $p (0..$#{$pepInfo{$queryID}}) {
			next unless $pepInfo{$queryID}[$p];
			my ($sel)=($pepInfo{$queryID}[$p]=~/SEL=(-*\d)/);
			next if ($sel==-1 && $item ne 'query');
			$sthMP->execute($queryNum{$queryID},$rankPos{$queryID}[$p]);
			my $numMatches=0;
			while (my ($identifier)=$sthMP->fetchrow_array) {$numMatches++ if $identList{$identifier};}
			$numSelProt{"$queryID:$p"}=$numMatches;
		}
	}
	$sthMP->finish;
}

&updateWaitBox;

#######################################################################################
####>Assigning a top score to a given sequence+varMod and grouping interpretations<####
#######################################################################################
my (%scoreHierarchy,%topPeptide);
my $existPepData=0;
#if($ptmFilter ne 'None' && $ptmFilter ne 'Indifferent' && $ptmFilter ne 'Any'){
#	$ptmFilter=~s/\)//;
#	$ptmFilter=~s/\(/\\(/;
#}
my $qPTMFilter=quotemeta($ptmFilter);
foreach my $pep (sort{&sortPeptides("$msType:default")} keys %scoreList) {
	$existPepData++;
	my ($queryID,$p)=split(/:/,$pep);
	my ($seq)=($pepInfo{$queryID}[$p]=~/SEQ=(\w+)/);
	my $vMod=$varModList{$pep}; #&promsMod::toStringVariableModifications($dbh,'QUERY',$queryID,$analysisID,$seq,$rankPos{$queryID}[$p]);
	##########################################################################################
	####>Print only the queries with the PTM that interest the massist in validation mode<####
	##########################################################################################
	if ($item ne 'query') {
		if ($ptmFilter eq 'Any'){
			next unless ($vMod || ($fixRes && $seq =~ /$fixRes/));
		}elsif ($ptmFilter eq 'None'){
			next if($vMod || ($fixRes && $seq =~ /$fixRes/));
		}elsif ($ptmFilter ne 'Indifferent'){
			if ($fixRes) { # fix modification
				next unless $seq =~ /$fixRes/;
			}
			elsif(!$vMod || ($vMod && $vMod !~ m/$qPTMFilter/)){
				next;
			}
		}
	}
	################################################################
	####>Print only the queries in the range chosen by the user<####
	################################################################
	#my $isInRange;
	##TODO: find the beginning and the end of the peptides matching the sequence so as to print it in listRanks.cgi
	#foreach my $match (split(/:/,$pepMatchInfo{"$queryNum{$queryID}:$rankPos{$queryID}[$p]"})){
	#	my ($beg,$aa1,$aa2)=split(/,/,$match);
	#	if ($beg >= $pepRangeMin && $beg+length($seq)-1 <= $pepRangeMax) {
	#		$isInRange=1;
	#		last;
	#	}
	#}
	#next unless $isInRange;
	$vMod='' unless $vMod;
	my $seqMod="$seq$vMod";
	if (!$topPeptide{$seqMod}) {
		@{$scoreHierarchy{$pep}}=();
		$topPeptide{$seqMod}=$pep;
	}
	elsif ($scoreList{$pep} > $scoreList{$topPeptide{$seqMod}}) {
		#>Copy old top contents to new top
		@{$scoreHierarchy{$pep}}=@{$scoreHierarchy{$topPeptide{$seqMod}}};
		push @{$scoreHierarchy{$pep}},$topPeptide{$seqMod};
		#>Delete previous top
		@{$scoreHierarchy{$topPeptide{$seqMod}}}=();
		delete $scoreHierarchy{$topPeptide{$seqMod}};
		#>Record new top
		$topPeptide{$seqMod}=$pep;
	}
	else {push @{$scoreHierarchy{$topPeptide{$seqMod}}},$pep;}
}
$dbh->disconnect;

###############################
####>Setting display order<####
###############################
my @displayOrder;
my $usedMsType=($msType eq 'PMF')? 'PMF' : 'MIS';
foreach my $topPep (sort{&sortPeptides("$usedMsType:$rankSort")} keys %scoreHierarchy) { # user-selected order
	push @displayOrder,$topPep;
	foreach my $lowPep (sort{&sortPeptides("$usedMsType:$rankSort")} @{$scoreHierarchy{$topPep}}) { # score order
		push @displayOrder,$lowPep;
	}
}

#########################################
####>Finding the 1st unselected rank<####
#########################################
my ($selectedRankID,$defaultRankID);
my $tableData; # flag for data to be displayed in table (in case query with no interpretations) ... <- Not necessary because of min score threshold
foreach my $pep (@displayOrder) {
	my ($queryID,$p)=split(/:/,$pep);
	next unless $pepInfo{$queryID}[$p]; # skip rank0 (used to be 'last')
	$tableData=1;
	my $rank=$rankPos{$queryID}[$p];
	my $itemTag=($disableLink || $scoreList{$pep}==0  || $fileFormat eq 'TDM.PEP.XML' || $fileFormat eq 'TDM.XML:NoSpec')? 'prot' : 'seq';
	$defaultRankID=join '_',($itemTag,$queryID,$queryNum{$queryID},$rank) unless $defaultRankID; # 1st pep that matches something
	last if $disableLink;
	if ($pepInfo{$queryID}[$p]=~/SEL=0/) { # peptide is not selected
		$selectedRankID=join '_',($itemTag,$queryID,$queryNum{$queryID},$rank);
		last;
	}
}
if ($existPepData) {$tableData=1};#To distinguish the case between no interpretation and no PTMs found
$selectedRankID=$defaultRankID if !$selectedRankID;
$selectedRankID='N/A' if !$defaultRankID; # query with no interpretation


#########################
####>Restarting HTML<####
#########################
my $protString='';
my $pepRangeString='';
if ($item ne 'query') {
	$protString="&PROT=$identifier";
	$pepRangeString="&pepRange=$pepRangeMin:$pepRangeMax";
}
my $protectedIden=uri_escape($identifier);
print qq
|<SCRIPT language="JavaScript">
parent.listRankSort ='$rankSort';
parent.listRankDlt ='$DeltaType';
parent.varValue = '$oldPtmFilter';
//parent.maxDisp = '$maxRankDisplay';
//parent.maxSrch = '$maxRankSearch';
var pepRangeString="$pepRangeString";
function drawSpectrum(rankId) {
	selectObject(rankId);
	parent.spectrumFrame.location="$promsPath{cgi}/drawSpectrum.cgi?IDENT=$protectedIden&RID="+rankId;
}
function editVarModification(vmodId,seqLength) {
	//highlight var mod
	if (selectedVarModId) {
		document.getElementById(selectedVarModId).style.color='#000000';
	}
	document.getElementById(vmodId).style.color='#DD0000'; //style.background = '#DD0000';
	selectedVarModId=vmodId;
	//popup window
	var winLocation="$promsPath{cgi}/editVarModification.cgi?CALL=listRanks&ID="+vmodId;
	var winWidth=400+seqLength*25;
	if (!varModWindow \|\| varModWindow.closed) {
		varModWindow=window.open(winLocation,'VarModWindow','width='+winWidth+',height=300,location=no,resizable=yes,scrollbars=yes');
		varModWindow.moveTo(250,300);
	}
	else {
		varModWindow.resizeTo(winWidth,250);
		varModWindow.location=winLocation;
	}
	varModWindow.focus();
}
function matchProteins(matchId,sequence) {
	selectObject(matchId);
	parent.spectrumFrame.location="$promsPath{cgi}/listMatchedProteins.cgi?ANAID=$analysisID&MSTYPE=$msType&ITEM=rank&SHOWFILT=$showFiltered&DIS_EXCLUDE=$disableExclude&MATCH="+matchId+"&SEQ="+sequence+"$protString";
}
function matchAnalyses(matchId,matchInfo,matchType) {
	var extraID=(matchType)? ':'+matchType : '';
	selectObject(matchId+extraID);
	parent.spectrumFrame.location="$promsPath{cgi}/listMatchedAnalyses.cgi?MSTYPE=$msType&START=1&DISABFORM=$disableLink&MATCH="+matchId+"&INFO="+matchInfo+"&ChargSts="+matchType;
}
function matchPattern(matchId) {
	selectObject(matchId);
	parent.spectrumFrame.location="$promsPath{cgi}/listMatchedProteins.cgi?ANAID=$analysisID&MSTYPE=$msType&ITEM=list&PROT=$identifier&SHOWFILT=$showFiltered&DIS_EXCLUDE=$disableExclude";
}
function selectObject(newId) {
	var oldObject=document.getElementById(selectedId);
	if (oldObject) {oldObject.style.color='#000000';}
	var newObject=document.getElementById(newId);
//	if (newObject) {
		newObject.style.color='#DD0000'; //style.background = '#DD0000';
		selectedId=newId;
//	}
}
function sequenceView(analysisId,proteinId,idType,ms){
	var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana=$analysisID&call=validation&id_prot=$itemID&id_type=notValid&msdata=0";
	top.openProtWindow(winLocation);
}
function changeProtStatus(selStatus) {
	var proceed;
	if (selStatus==-3) {
		if (confirm('Protein has been removed from list with a filter!\\nDo you really want to include it ?')) {proceed=1;}
	} else {proceed=1;}
	if (proceed) {
		window.location="./listRanks.cgi?ITEM=$item&ID=$itemID&ANAID=$analysisID&MSTYPE=$msType&IDENTIF=$identifier&SHOWFILT=$showFiltered&PAGE=$page-$maxPage&RANGE=$rangeStrg&STATUS="+selStatus+"&SORT="+parent.selectedSort+"&sortType="+parent.listRankSort+"&Dlt="+parent.listRankDlt+"&varValue="+parent.varValue+"$pepRangeString&multiAna=$multiAna";
	}
}
function checkPeptide(box1,box2,pepSel) {
	activateForm();
	if (box1.checked) {box2.checked=false;}
	else if (pepSel != 0 && pepSel != 1) {box2.checked=true;} // if already selected or rejected: cannot go back to neutral => uncheck box1=check box2
}
function checkAllPeptides(myForm,type) {
	var allFormInputs=myForm.getElementsByTagName('INPUT');
	for (var i=0; i<allFormInputs.length; i++) {
		if (allFormInputs[i].type=='checkbox') {
			var pepData=allFormInputs[i].name.match(/sel_(\\d+_\\d+_\\d+)/);
			if (pepData != null) { //[0]=whole matched string
				var pepInfoStrg=document.getElementById('pepInfo_'+pepData[1]).value;
				var selStatus=pepInfoStrg.match(/SEL=(-*\\d)/);
				if (type=='select') {
					if (allFormInputs[i].checked==true) continue;
					allFormInputs[i].checked=true;
					checkPeptide(allFormInputs[i],allFormInputs[i+1],selStatus[1]);
				}
				else {
					if (allFormInputs[i+1].checked==true) continue;
					allFormInputs[i+1].checked=true;
					checkPeptide(allFormInputs[i+1],allFormInputs[i],selStatus[1]);
				}
			}
		}
	}
}
function editPepComments(queryId,queryNum,rank) {
	var pepCom='pepCom_'+queryId+'_'+queryNum+'_'+rank;
	var comments=document.getElementById(pepCom).value;
	var newComments=prompt('Type comments for peptide selection :',comments);
	if (newComments==null) {newComments='';}
	document.getElementById(pepCom).value=newComments;
}
function popComments(queryId,queryNum,rank) {
	var pepCom='pepCom_'+queryId+'_'+queryNum+'_'+rank;
	var comments=document.getElementById(pepCom).value;
	if (comments) {popup(comments);}
	else {popup('No comments.');}
}
function activateForm() {
	if (activeForm) {return;}
	document.rankForm.save.style.fontWeight='bold';
	document.rankForm.save.style.color='#000000';
	activeForm=1;
}
function getSort() {
	document.rankForm.SORT.value=parent.selectedSort;
}
function setRowVisibility() {
	var img=document.getElementsByName('plusMinus')[0];
//	var rowList=document.getElementsByName('noSelect');
	var rowList=document.getElementsByTagName('table');
//alert(rowList.length);
	if (rowDisplay=="none") { // not selectable rows are hidden => show
		img.src='$promsPath{images}/minus1.gif';
		rowDisplay='block';
	}
	else { // not selectable rows are shown => hide
		img.src='$promsPath{images}/plus.gif';
		rowDisplay='none';
	}
	for (var i=0;i<rowList.length;i++) {
		if (rowList[i].getAttribute('id')) {
			rowList[i].style.display=rowDisplay;
		}
	}
}
function findPageAndUpdate(queryId) {
	var page;
	for (var i=0; i<parent.queryPages.length; i++) {
		var pageRange=parent.queryPages[i].split(/-/);
		if (queryId>=pageRange[0] && queryId<=pageRange[1]) {
			page=i+1;
			break;
		}
	}
	parent.updateValidation(queryId,0,0,0,0,'query',page,1);
}
function validMatchGroup(action) {
	var okToGo=false;
	if (action=='validate') {
		if (confirm('Select all interpretions with score >= '+document.rankForm.MGscore.value+' for current Match Group?')) {okToGo=true;}
	}
	else if (confirm('Clear all selections for current Match Group?')) {okToGo=true;}
	if (okToGo==true) {
		window.location="./listRanks.cgi?ANAID=$analysisID&MSTYPE=$msType&ITEM=$item&ID=$itemID&SHOWFILT=$showFiltered&PAGE=$page-$maxPage&RANGE=$rangeStrg&valMG="+action+"&MGscore="+document.rankForm.MGscore.value+"&SORT="+parent.selectedSort+"&sortType="+parent.listRankSort+"&Dlt="+parent.listRankDlt+"&varValue="+parent.varValue+"$pepRangeString&multiAna=$multiAna";
	}
}
function clearPMFGroup() {
	window.location="./listRanks.cgi?ANAID=$analysisID&MSTYPE=$msType&ITEM=$item&ID=$itemID&SHOWFILT=$showFiltered&PAGE=$page-$maxPage&RANGE=$rangeStrg&valMG=clear&SORT="+parent.selectedSort+"&sortType="+parent.listRankSort+"&Dlt="+parent.listRankDlt+"&varValue="+parent.varValue+"$pepRangeString&multiAna=$multiAna";
}
function selectSort(rankSort) {
	if (rankSort != parent.listRankSort) {
		window.location="./listRanks.cgi?ANAID=$analysisID&MSTYPE=$msType&ITEM=$item&SORT=$selectedSort&ID=$itemID&SHOWFILT=$showFiltered&PAGE=$page-$maxPage&RANGE=$rangeStrg&sortType="+rankSort+"&Dlt="+parent.listRankDlt+"&varValue="+parent.varValue+"$pepRangeString&multiAna=$multiAna";
	}
}
function multiAnaScan(chkStatus) {
	parent.multiAnaScan=(chkStatus)? 1 : 0;
	window.location="./listRanks.cgi?ANAID=$analysisID&MSTYPE=$msType&ITEM=$item&SORT=$selectedSort&ID=$itemID&SHOWFILT=$showFiltered&PAGE=$page-$maxPage&RANGE=$rangeStrg&sortType="+parent.listRankSort+"&Dlt="+parent.listRankDlt+"$pepRangeString&multiAna="+parent.multiAnaScan;
}
function selectShow(showValue) {
	if (showValue != parent.listRankDlt) {
		window.location="./listRanks.cgi?ANAID=$analysisID&MSTYPE=$msType&ITEM=$item&SORT=$selectedSort&ID=$itemID&SHOWFILT=$showFiltered&PAGE=$page-$maxPage&RANGE=$rangeStrg&sortType="+parent.listRankSort+"&Dlt="+showValue+"$pepRangeString&multiAna=$multiAna";
	}
}
function activatePeptide(rankID) {
	var pepData=rankID.split('_');
	if (confirm('Make peptide corresponding to query#'+pepData[1]+' rank#'+pepData[2]+' selectable?')) {
		if (!activeForm \|\| confirm('Last changes made to selection will not be saved! Proceed anyway?')) {
			window.location="./listRanks.cgi?ANAID=$analysisID&MSTYPE=$msType&ITEM=$item&SORT=$selectedSort&ID=$itemID&SHOWFILT=$showFiltered&PAGE=$page-$maxPage&RANGE=$rangeStrg&actPep="+rankID+"&sortType="+parent.listRankSort+"&Dlt="+parent.listRankDlt+"&maxSrch=$maxRankSearch&maxDisp=$maxRankDisplay"+"&varValue="+parent.varValue+"$pepRangeString&multiAna=$multiAna";
		}
	}
}
function changeMaxDisp(myItem) {
	parent.maxDisp = myItem.value;
	window.location="./listRanks.cgi?ANAID=$analysisID&MSTYPE=$msType&ITEM=$item&SORT=$selectedSort&ID=$itemID&SHOWFILT=$showFiltered&PAGE=$page-$maxPage&RANGE=$rangeStrg&sortType=$rankSort&Dlt=$DeltaType&maxSrch=$maxRankSearch&maxDisp="+myItem.value+"&varValue="+parent.varValue+"$pepRangeString&multiAna=$multiAna";
}
function changeMaxSearch (myItem) {
	parent.maxSrch=myItem.value;
	window.location="./listRanks.cgi?ANAID=$analysisID&MSTYPE=$msType&ITEM=$item&SORT=$selectedSort&ID=$itemID&SHOWFILT=$showFiltered&PAGE=$page-$maxPage&RANGE=$rangeStrg&sortType=$rankSort&Dlt=$DeltaType&maxDisp=$maxRankDisplay&maxSrch="+myItem.value+"&varValue="+parent.varValue+"$pepRangeString&multiAna=$multiAna";
}
function selectVarModView (varValue)  {
	if (varValue != parent.varValue) {
		parent.varValue=varValue;
		window.location="./listRanks.cgi?ANAID=$analysisID&MSTYPE=$msType&ITEM=$item&SORT=$selectedSort&ID=$itemID&SHOWFILT=$showFiltered&PAGE=$page-$maxPage&RANGE=$rangeStrg&sortType="+parent.listRankSort+"&Dlt="+parent.listRankDlt+"&maxSrch=$maxRankSearch&maxDisp=$maxRankDisplay&varValue="+parent.varValue+"$pepRangeString&multiAna=$multiAna";
	}
}
function applyPepRange () {
	var pepRangeMin=document.rankForm.rangeMin.value;
	var pepRangeMax=document.rankForm.rangeMax.value;
	window.location="./listRanks.cgi?ANAID=$analysisID&MSTYPE=$msType&ITEM=$item&SORT=$selectedSort&ID=$itemID&SHOWFILT=$showFiltered&PAGE=$page-$maxPage&RANGE=$rangeStrg&sortType="+parent.listRankSort+"&Dlt="+parent.listRankDlt+"&maxSrch=$maxRankSearch&maxDisp=$maxRankDisplay&varValue="+parent.varValue+"&multiAna=$multiAna&pepRange="+pepRangeMin+":"+pepRangeMax;
}
document.getElementById('waitDiv').style.display='none';
var selectedId='$selectedRankID';
var selectedVarModId;
var varModWindow; // window used to edit variable modifications
var activeForm=0;
var rowDisplay='none';
</SCRIPT>
<FORM name="rankForm" method="post"> <!--onSubmit="return(checkForm());"-->
<INPUT type="hidden" name="ANAID" value="$analysisID" />
<INPUT type="hidden" name="MSTYPE" value="$msType" />
<INPUT type="hidden" name="ITEM" value="$item" />
<INPUT type="hidden" name="ID" value="$itemID" />
<INPUT type="hidden" name="SORT" value="" />
<INPUT type="hidden" name="PAGE" value="$page-$maxPage" />
<INPUT type="hidden" name="RANGE" value="$rangeStrg" />
|;
# print "<BODY topmargin=2px bottommargin=0px leftmargin=5px onLoad=\"$onLoadString\" background=\"$promsPath{images}/bgProMS.gif\">\n"; # does not work if put in HEAD
# print "<FORM name=\"rankForm\" method=\"post\"> # onSubmit="return(checkForm());"
# print hidden(-name=>'ANAID',-default=>"$analysisID"),"\n";
# print hidden(-name=>'MSTYPE',-default=>"$msType"),"\n";
# print hidden(-name=>'ITEM',-default=>"$item"),"\n";
# print hidden(-name=>'ID',-default=>"$itemID"),"\n";
# print hidden(-name=>'SORT',-default=>""),"\n";
# print hidden(-name=>'PAGE',-default=>"$page-$maxPage"),"\n";
# print hidden(-name=>'RANGE',-default=>"$rangeStrg"),"\n";
# print hidden(-name=>'modifName',-default=>''),"\n";
foreach my $queryID (sort{$a<=>$b} keys %queryNum) {
	print "<INPUT type=\"hidden\" name=\"status_$queryID\" value=\"$queryStatus{$queryID}\" />\n";
	foreach my $p (0..$#{$pepInfo{$queryID}}) {
		next unless $pepInfo{$queryID}[$p]; # when no rank0
# 		next if $pepInfo{$queryID}[$p]=~/SEL=-1/;
		my $rank=$rankPos{$queryID}[$p];
		my $rankID="$queryID"."_$queryNum{$queryID}"."_$rank"; # item -> query !!!!!!!!!!!!
		print "<INPUT type=\"hidden\" name=\"pepInfo_$rankID\" id=\"pepInfo_$rankID\" value=\"$pepInfo{$queryID}[$p]\" />\n";
		next if $pepInfo{$queryID}[$p]=~/SEL=-1/; # better score exists
		print "<INPUT type=\"hidden\" name=\"pepCom_$rankID\" id=\"pepCom_$rankID\" value=\"$commentList{\"$queryID:$p\"}\" />\n";
	}
}
if ($item eq 'query') {
	print "<FONT class=\"title2\">Interpretations of <FONT color=#DD0000>query $queryNum{$itemID}</FONT> :</FONT>";
	printf "<FONT class=\"title3\">&nbsp;&nbsp;Mr(exp): %.2f Da &nbsp;&nbsp; Mr(obs): %.2f Da.",$massExp{$itemID},$massObs{$itemID};
	print "&nbsp;&nbsp;Elution: $elutionTime{$itemID}</FONT>\n";
}
else { # protein or myList
	my $identifierStrg=($numDatabanks>1)? $dbRank.'::'.$identifier : $identifier;
	print "<FONT class=\"title2\">Peptides matching protein <FONT color=#DD0000>$identifierStrg</FONT></FONT>";
	if ($weight==0) {print '<FONT class="title3">&nbsp;&nbsp;Mass: - Da';}
	else {printf "<FONT class=\"title3\">&nbsp;&nbsp;Mass: %.0f Da",$weight;}
	my $scoreLabel=($msType eq 'PMF')? 'Max. Score' : 'Score';
	printf "&nbsp;&nbsp;&nbsp;$scoreLabel: %.1f",$score if $selStatus>-2;
	print "</FONT>\n";
	print "&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" value=\"More info\" onclick=\"sequenceView()\">\n";
	my ($action,$text);
	if ($selStatus<=-2) { # (excluded)
		print "<BR><FONT style=\"font-weight:bold;text-decoration:underline;\">Description:</FONT> $protDes. <FONT class=\"org\">$protOrg</FONT>\n";
		my $word=($selStatus==-2)? 'excluded' : 'filtered';
		print "<BR><FONT class=\"title3\" color=#DD0000>This protein has been $word from the list of selectable proteins.</FONT>\n";
		print "&nbsp;&nbsp;&nbsp;<INPUT type='button' value=\"Include in list\" onclick=\"changeProtStatus($selStatus)\" $disableExclude><BR>\n";
# 		print end_html; # do not list matches
# 		exit;
	}
	else {
		print "&nbsp;<INPUT type=\"button\" value=\"Exclude from list\" onclick=\"changeProtStatus($selStatus)\" $disableExclude><BR>\n";
		print "<FONT style=\"font-weight:bold;text-decoration:underline;\">Description:</FONT> $protDes. <FONT class=\"org\">$protOrg</FONT>\n";
	}
}

if ($item eq 'query' && !$tableData) { # query has no interpretations
		print qq
|<BR><FONT class="title2">This query has no interpretations or a score below threshold.</FONT><BR>
<SCRIPT type="text/javascript">parent.spectrumFrame.location='$promsPath{html}/nothing.html'</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

my ($color1,$color2)=&promsConfig::getRowColors;
my $bgColor=$color2;
my $display=($item ne 'query')? 'none' : 'block';

print qq
|<TABLE width=1105 cellspacing=0 cellpadding=0>
<TR class="header" bgcolor="$color2">
	<TH width=40>Auto</TH>
	<TH width=30><IMG class="button" src="$promsPath{images}/good.gif" onclick="checkAllPeptides(document.rankForm,'select')" onmouseover="popup('Click to <B>select</B> all peptides.')" onmouseout="popout()"></TH>
	<TH width=30><IMG class="button" src="$promsPath{images}/bad.gif" onclick="checkAllPeptides(document.rankForm,'reject')" onmouseover="popup('Click to <B>reject</B> all peptides.')" onmouseout="popout()"></TH>
|;
my %sortColor;
foreach my $sortMode ('query','score','delta','ppm','charge','sequence') {
	$sortColor{$sortMode}=($sortMode eq $rankSort)? '#0000BB' : 'black';
}
if ($item ne 'query') {
	print "\t<TH width=70><A href=\"javascript:selectSort('query')\" style=\"font-style:italic;color:$sortColor{query}\" onmouseover=\"popup('Click to sort peptides by <B>ascending query number</B>.')\" onmouseout=\"popout()\">Query</A></TH>\n";
}
print "\t<TH width=45><ACRONYM onmouseover=\"popup('Click on peptides <B>rank</B> to add <B>comments</B> regarding validation.')\" onmouseout=\"popout()\"><I>Rank</I></ACRONYM></TH>\n";

if ($msType eq 'PMF') {print "\t<TH width=45>Score</TH>\n";}
else {
	print "\t<TH width=55><A href=\"javascript:selectSort('score')\" style=\"font-style:italic;color:$sortColor{score}\" onmouseover=\"popup('Click to sort peptides by <B>descending score</B>.')\" onmouseout=\"popout()\">Score</A></TH>\n";
}

if ($item eq 'query') {
	print "\t<TH width=65>Mr(calc)</TH>\n"; # new
}
else {
	print "\t<TH width=65>Mr(obs)</TH>\n"; #Observed
}
# if ($showPPM) {
	# print "\t<TH width=45><A href=\"javascript:selectSort('ppm')\" style=\"font-style:italic;color:$sortColor{ppm}\" onmouseover=\"popup('Click to sort peptides by <B>ascending mass delta (ppm)</B>.')\" onmouseout=\"popout()\">PPM</A></TH>\n";
# }
# else {
	# print "\t<TH width=45><A href=\"javascript:selectSort('delta')\" style=\"font-style:italic;color:$sortColor{delta}\" onmouseover=\"popup('Click to sort peptides by <B>ascending absolute mass delta (Da)</B>.')\" onmouseout=\"popout()\">Delta</A></TH>\n";
# }
unless ($showPPM) {
	print "\t<TH width=65><A href=\"javascript:selectShow('ppm')\" style=\"font-style:italic;color:$sortColor{ppm}\" onmouseover=\"popup('Click to show delta in <B>relative</B> mode (ppm).')\" onmouseout=\"popout()\">Delta</A></TH>\n";
}
else {
	print "\t<TH width=65><A href=\"javascript:selectShow('Da')\" style=\"font-style:italic;color:$sortColor{delta}\" onmouseover=\"popup('Click to show delta in <B>absolute</B> mode (Da).')\" onmouseout=\"popout()\">ppm</A></TH>\n";
}
print "\t<TH width=40>Miss</TH>\n";
my $checkMultiAna=($multiAna)? 'checked' : '';
if ($item eq 'query') {
	print "\t<TH width=50>Prot.</TH>\n";
	print "\t<TH width=80><INPUT type=\"checkbox\" name=\"multiAna\" value=\"1\" onclick=\"multiAnaScan(this.checked)\" $checkMultiAna/>&nbsp;<ACRONYM onmouseover=\"popup('Click on numbers to list other analyses containing the <B>same peptide</B>.')\" onmouseout=\"popout()\"><I>Ana.</I></ACRONYM></TH>\n";
	print "\t<TH width=70><A href=\"javascript:selectSort('charge')\" style=\"font-style:italic;color:$sortColor{charge}\" onmouseover=\"popup('Click to sort peptides by <B>charge</B>.')\" onmouseout=\"popout()\">Charge</A></TH>\n";
	print "\t<TH align=left>&nbsp&nbsp;"; # width=470
	print "<A href=\"javascript:selectSort('sequence')\" style=\"font-style:italic;color:$sortColor{sequence}\" onmouseover=\"popup('Click to sort peptides by <B>sequence</B>.')\" onmouseout=\"popout()\">Peptide sequence</A></TH>\n";
	print "</TR></TABLE>\n";
}
else { # protein, myList
	print "\t<TH width=50><A id=\"match_list\" href=\"javascript:matchPattern('match_list')\" onmouseover=\"popup('Click to view proteins with similar match patterns than <B>$identifier</B>.')\" onmouseout=\"popout()\"><I>Prot.</I></A></TH>\n";
	print "\t<TH width=80><INPUT type=\"checkbox\" name=\"multiAna\" value=\"1\" onclick=\"multiAnaScan(this.checked)\" $checkMultiAna/>&nbsp;<ACRONYM onmouseover=\"popup('Click on numbers to list other analyses containing the <B>same peptide</B>.')\" onmouseout=\"popout()\"><I>Ana.</I></ACRONYM></TH>\n";
	print "\t<TH width=70><A href=\"javascript:selectSort('charge')\" style=\"font-style:italic;color:$sortColor{charge}\" onmouseover=\"popup('Click to sort peptides by <B>charge</B>.')\" onmouseout=\"popout()\">Charge</A></TH>\n";
	print "\t<TH align=left valign=bottom>&nbsp;"; # width=470
	print "<IMG name=\"plusMinus\" src=\"$promsPath{images}/plus.gif\" onclick=\"javascript:setRowVisibility()\" onmouseover=\"popup('Click to show/hide <B>lower scoring</B> interpretations.')\" onmouseout=\"popout()\" width=16 height=22 align=top>" if $noSelect;

	print "<A href=\"javascript:selectSort('sequence')\" style=\"font-style:italic;color:$sortColor{sequence}\" onmouseover=\"popup('Click to sort peptides by <B>sequence</B>.')\" onmouseout=\"popout()\">Peptide sequence</A> with PTMs:<SELECT name=\"ptmFilter\" class=\"title3\" onchange=\"selectVarModView(this.value)\">\n";
	foreach my $varMod (@modifications) {
		my $selMod=($varMod eq $oldPtmFilter)? 'selected' : '';
		print "<OPTION value=\"$varMod\" $selMod>$varMod</OPTION>";
	}
	foreach my $modID (sort{lc($ptmList{$a}{'NAME'}) cmp lc($ptmList{$b}{'NAME'})} keys %ptmList) {
		my $selMod=($ptmList{$modID}{'NAME'} eq $ptmFilter)? 'selected' : '';
		my $extraName=($ptmList{$modID}{'MODIF_TYPE'} eq 'F') ? "$ptmList{$modID}{'SPECIFICITY'}::" : '';
		print "<OPTION value=\"$extraName$ptmList{$modID}{'NAME'}\" $selMod>$ptmList{$modID}{'NAME'}</OPTION>";
	}
	print "</SELECT></TD>\n</TR>\n";
	print "</TR></TABLE>\n";
}

#########################################
####>Looping through interpretations<####
#########################################
my $rankIdList='';
my $numSelectable;
my $prevPeptide='';
my %sequencePosStrg; # position of sequence on protein (not if item is 'query')
$dbh=&promsConfig::dbConnect;
foreach my $pep (@displayOrder) {
	my ($queryID,$p)=split(/:/,$pep);
	my $rank=$rankPos{$queryID}[$p];
	my $info=$pepInfo{$queryID}[$p];
	my $charge=$queryCharge{$queryID};
	if (!$info) { # only if query view
		print "<TABLE><TR><TH align=left>&nbsp&nbsp;<I>No more interpretations</I></TH></TR></TABLE>\n";
		last;
	}
	my ($select)=($info=~/SEL=(-?\d)/);
	my ($searchRank)=($info=~/SRK=(\d+)/); $searchRank=$rank unless $searchRank; # mixed decoy DB
	my ($missCut)=($info=~/MIS=(\d)/);
	my ($sequence)=($info=~/SEQ=(\w+)/);
	my ($numProt)=($info=~/MATCH=(\d+)/);
	my ($deltaDa)=($info=~/DELT=([^,]+)/); # change in 03/12/13 order to get all for Paragon MATCH that could be written like this DELT=-2.18348e-005
	my ($massCalc)=($info=~/CALC=([^,]+)/);
	my $varMod=$varModList{$pep}; #&promsMod::toStringVariableModifications($dbh,'QUERY',$queryID,$analysisID,$sequence,$rank); #
	my ($sub)=($info=~/SUBST=([^,]+)/);
	if ($sub) {
		my %susbtitutions; # Could exist several substitutions on the same peptide !
		foreach my $subst (split(' \+',$sub)) {
			next unless $subst;
			if ($subst =~ /(\w)->(\w) \((\d+)\)/){
			    my ($origin,$substituted,$pos)=($1,$2,$3);
			    $susbtitutions{$pos}="($1&rarr;$2)";
			}
		}
		foreach my $pos (sort {$b <=> $a} keys %susbtitutions) { # Modify the sequence in descending order so as to keep the AA position integrity (because string replacement size is bigger thant one letter... )
 		    substr($sequence , $pos-1 , 1 , $susbtitutions{$pos} );
		}
	}
	my ($filter)=($info=~/FLT=(-?\d)/); $filter=0 unless $filter;
	my $qValityStrg;
	if ($rank==1) {
		my ($PEP,$qValue)=($info=~/PEP=(.+),QVAL=([^,]+)/);
		if ($PEP) {
			$qValityStrg=($qValue==0)? sprintf "<B>q-value:</B> $qValue<BR><B>Post. Err. Prob.:</B> %1.2e",$PEP : sprintf "<B>q-value:</B> %1.2e<BR><B>Post. Err. Prob.:</B> %1.2e",$qValue,$PEP;
		}
	}
	#chomp($varMod);
	my $sequence2=quotemeta($sequence);# protection of metacharacter like & or () that can be put in the sequence (substitution)
	if (length($sequence)>47) {
		my ($subSeq,@newSeq);
		do {
			$subSeq=substr($sequence2,0,46);
			if ($subSeq) {
				push @newSeq,$subSeq;
				$sequence2=~s/$subSeq//;
			}
		} until !$subSeq;
		$sequence2=join('\<BR>',@newSeq);
	}
	#$sequence2=unquotemeta($sequence2);
	$sequence2=~s/\\//g;
	$sequence2=~s/[^<]([X,B,Z])/<FONT style="font-size=15px">$1<\/FONT>/g; # ambigous residue(s) shown bigger. [^<] prevents \<BR> (if any) to be recognized by match.

	my $rankID="$queryID"."_$queryNum{$queryID}"."_$rank";
	$rankIdList.="$rankID:" unless $select==-1; # better score exists

	##>Checking Selection status
	my ($autoString,$selString,$rejString);
	if ($select==-1) { # no checkbox
		if ($item eq 'query' || "$sequence$varMod" ne $prevPeptide) { # seqVarMod may be diff from prevPep if rank filter is set ("parent" pep did not pass but hidden one did)
			$bgColor=($bgColor eq $color1)? $color2 : $color1;
		}
		print "<TABLE id=\"noSelect_$rankID\" class=\"list\" width=1103 bgcolor=\"$bgColor\" cellspacing=0 cellpadding=0 style=\"display:$display\"><TR>\n"; # width=$width
		$autoString='';
		$selString=($disableForm)? '-' : "<IMG src=\"$promsPath{images}/plus.gif\" onclick=\"javascript:activatePeptide('$rankID')\" onmouseover=\"popup('Click to make peptide <B>selectable</B>.')\" onmouseout=\"popout()\" width=16 height=22 align=top>";
		$rejString='-';
	}
	else { # checkbox
		if ("$sequence$varMod" ne $prevPeptide) {
			$bgColor=($bgColor eq $color1)? $color2 : $color1;
		}
		print "<TABLE class=\"list\" width=1103 bgcolor=\"$bgColor\" cellspacing=0 cellpadding=0><TR>\n"; # width=$width
		my $autoImg=($select==1)? 'good.gif' : ($select==-2)? 'bad.gif' : 'blank.gif';
		my $checkSel=($select==2)? 'checked' : ''; # 1
		my $checkRej=($select==-3)? 'checked' : ''; # -2
		$autoString="<IMG src=\"$promsPath{images}/$autoImg\">";
		$selString="<INPUT type=\"checkbox\" name=\"sel_$rankID\" $checkSel onclick=\"checkPeptide(document.rankForm.sel_$rankID,document.rankForm.rej_$rankID,$select)\" $disableForm/>";
		$rejString="<INPUT type=\"checkbox\" name=\"rej_$rankID\" $checkRej onclick=\"checkPeptide(document.rankForm.rej_$rankID,document.rankForm.sel_$rankID,$select)\" $disableForm/>";
		$numSelectable++;
	}
	$prevPeptide="$sequence$varMod";

	##>Checkboxes
	print "\t<TH width=42>$autoString</TH>\n";
	print "\t<TH width=32>$selString</TH>\n"; # select 52
	print "\t<TH width=32>$rejString</TH>\n"; # reject

	##>Query number
	if ($item ne 'query') {
		my $imgString="&nbsp;<IMG src=\"$promsPath{images}";
		#if ($queryStatus{$queryID}==-1) {
		#	if ($filter==0) {$imgString.='/lightGray1.gif';}
		#	elsif  ($filter==-1) {$imgString.='/lightOrange1.gif';}
		#	elsif  ($filter==1) {$imgString.='/lightYellow1.gif';}
		#}
		if ($queryStatus{$queryID}==-1) {$imgString.='/lightGray1.gif';}
		elsif ($queryStatus{$queryID}==0 || $queryStatus{$queryID}==-3) {$imgString.='/lightRed1.gif';} # -3: worse scores
		elsif ($queryStatus{$queryID}==-2) {$imgString.='/lightGrayRed1.gif';} # reject + not verif
		else {$imgString.='/lightGreen1.gif';} # >=1
		$imgString.='" hspace=0 border=0 height=11 width=11>&nbsp;';
		print "\t<TD width=72 class=\"left\">$imgString";
		print "<A href=\"javascript:findPageAndUpdate($queryID)\" onmouseover=\"popup('<B>Elution:</B> $elutionTime{$queryID}<BR>Click to view all interpretations of <B>query $queryNum{$queryID}</B>.')\" onmouseout=\"popout()\">" unless $queryStatus{$queryID}==-4;
		print $queryNum{$queryID};
		print '</A>' unless $queryStatus{$queryID}==-4;
		print "</TD>\n";
	}

	##>Rank
	my $srkStrg=($searchRank==$rank)? '' : "<SMALL>:$searchRank</SMALL>";
	if ($select==-1) {print "\t<TD width=47>$rank$srkStrg</TD>\n";} # better score exists
	else {
		my $flagImgString=($filter==1)? "<IMG src=\"$promsPath{images}/lightYellow1.gif\">" : ($filter==-1)? "<IMG src=\"$promsPath{images}/lightOrange1.gif\">" : '';
		print "\t<TD width=47 class=\"right\">$flagImgString <A href=\"javascript:editPepComments($queryID,$queryNum{$queryID},$rank)\" onmouseover=\"popComments($queryID,$queryNum{$queryID},$rank)\" onmouseout=\"popout()\">$rank$srkStrg</A>&nbsp;&nbsp;&nbsp;&nbsp;";
# 		if ($filter) {
# 			my $image=($filter==-1)? '/lightOrange1.gif' : '/lightYellow1.gif';
# 			print "&nbsp<IMG src=\"$promsPath{images}/$image\" hspace=0 border=0 height=11 width=11>";
# 		}
		print "</TD>\n";
	}

	##>Score
	my $scoreString;
	if ($scoreList{$pep}==0) {$scoreString='-';} # PMF or PMF in mixed search
	else { # MS/MS
		$scoreString=($scoreHierarchy{$pep} || $select > -1)? sprintf "%.1f",$scoreList{$pep} : sprintf "(%.1f)",$scoreList{$pep};
	}
	if ($qValityStrg) {
		print "\t<TD width=57><A href=\"javascript:void(null)\" onmouseover=\"popup('$qValityStrg')\" onmouseout=\"popout()\">$scoreString</A></TD>\n";
	}
	else {print "\t<TD width=57>$scoreString</TD>\n";}

	##>Masses
	if ($item eq 'query') {
		print "\t<TD width=67><A href=\"javascript:void(null)\" onmouseover=\"popup('<B>Mr(calc):</B> $massCalc')\" onmouseout=\"popout()\">";
		printf "%.2f</A></TD>\n",$massCalc;
	}
	else {
		print "\t<TD width=67><A href=\"javascript:void(null)\" onmouseover=\"popup('<B>Mr(obs):</B> $massObs{$queryID}<BR><B>Mr(exp):</B> $massExp{$queryID}<BR><B>Mr(calc):</B> $massCalc<BR><B>Charge:</B> $queryCharge{$queryID}+')\" onmouseout=\"popout()\">";
		printf "%.2f</A></TD>\n",$massObs{$queryID};
	}

	##>Delta
	my $deltaPPM=sprintf "%.3f",1000000*($deltaDa/$massCalc);
	$deltaDa=sprintf  "%.1e",$deltaDa;
	if ($showPPM) {
		#>PPM
		print "\t<TD width=67><A href=\"javascript:void(null)\" onmouseover=\"popup('<B>Delta:</B> $deltaDa Da ($deltaPPM ppm)')\" onmouseout=\"popout()\">";
		printf "%.1f</A></TD>\n",$deltaPPM;
	}
	else {
		#>Dalton
		print "\t<TD width=67><A href=\"javascript:void(null)\" onmouseover=\"popup('<B>Delta:</B> $deltaDa Da ($deltaPPM ppm)')\" onmouseout=\"popout()\">";
		print "$deltaDa</A></TD>\n";
	}

	##>Missed cuts
	if ($select==-1 && $item ne 'query') {
		print "\t<TD width=42>-</TD>\n";
	}
	elsif ($existBetterPep{$pep}) {
		print "\t<TH width=42>*</TH>\n";
	}
	else {
		print "\t<TD width=42>$missCut</TD>\n";
	}

	##>Proteins
	if ($select==-1 && $item ne 'query') {
		print "\t<TD width=52>-</TD>\n";
	}
	elsif ($existBetterPep{$pep}) {
		print "\t<TH width=52>*</TH>\n";
	}
	else {
		print "\t<TD width=52>";
		if ($showFiltered) {
			print "<A id=\"prot_$rankID\" href=\"javascript:matchProteins('prot_$rankID','$sequence')\" onmouseover=\"popup('Click to view all proteins matched by <B>$sequence</B>.')\" onmouseout=\"popout()\">&nbsp;$numProt&nbsp</A>";
		}
		else { # hide filtered proteins
			my $plusString=($numSelProt{$pep}<$numProt)? "<SUP>+</SUP>" : '&nbsp'; # &nbsp($numProt)
			print "<A id=\"prot_$rankID\" href=\"javascript:matchProteins('prot_$rankID','$sequence')\" onmouseover=\"popup('Click to view all proteins matched by <B>$sequence</B>.')\" onmouseout=\"popout()\">&nbsp;$numSelProt{$pep}$plusString</A>";
		}
		print "</TD>\n";
	}

	##>Analyses
	print "\t<TH width=82>&nbsp;";
	if (!$multiAna || $select==-1 || $scoreList{$pep}==0) { # PMF
		print '-';
	}
	elsif ($existBetterPep{$pep}) {
		print '*';
	}
	else {
		my $numMatchAna=1; # add current MS/MS analysis
		my $extendedNumMatchAna=1 ;
		$numMatchAna+=scalar (@{$analysisMatch{$pep}}) if $analysisMatch{$pep};
		$extendedNumMatchAna+=scalar (@{$extendAnalysisMatch{$pep}}) if $extendAnalysisMatch{$pep};
		if ($numMatchAna>1) { # && !$disableLink) {
			my $matchInfo=join(':',@{$analysisMatch{$pep}});
			print "<A id=\"ana_$rankID\" href=\"javascript:matchAnalyses('ana_$rankID','$matchInfo','')\" onmouseover=\"popup('Click to view this peptide in other MS/MS Analyses.')\" onmouseout=\"popout()\">$numMatchAna</A>";
		}
		else {print $numMatchAna;}
		if ($extendedNumMatchAna>$numMatchAna) {
			print '/';
			if ($extendedNumMatchAna>1) { # && !$disableLink) {
				my $matchInfo=join(':',@{$extendAnalysisMatch{$pep}});
				print "<A id=\"ana_$rankID:all\" href=\"javascript:matchAnalyses('ana_$rankID','$matchInfo','all')\" onmouseover=\"popup('Click to view this peptide in other MS/MS Analyses with all charge states.')\" onmouseout=\"popout()\">$extendedNumMatchAna</A>";
			}
# 			else {print "$numMatchAna";} # cannot be
		}
	}
	print "&nbsp;</TD>\n";

	##>Charge
	print "\t<TD width=70 align='center'>$charge+</TD>\n";
	
	##>Sequence
	my $ref=($item ne 'query')? "$queryNum{$queryID}:$rank" : "$queryID:$rank"; #$pep;
	print "\t<TD class=\"left\">&nbsp;&nbsp;$boundaryAA{$ref}[0]"; # width=470
	if ($item eq 'query') {$sequencePosStrg{$sequence}='No further information';}
	elsif ($item ne 'query' && !$sequencePosStrg{$sequence}) {#Print the position information
		$sequencePosStrg{$sequence}='';
		foreach my $match (split(/:/,$pepMatchInfo{"$queryNum{$queryID}:$rank"})) {
			my($beg,$aa1,$aa2)=split(/,/,$match);
			my $end=$beg+length($sequence)-1;
			$sequencePosStrg{$sequence}.=",$beg-$end";
		}
		$sequencePosStrg{$sequence}=(substr $sequencePosStrg{$sequence}, 1);
		$sequencePosStrg{$sequence}="<B>Position:</B> $sequencePosStrg{$sequence}";
	}
	if ($scoreList{$pep}==0 || $fileFormat eq 'TDM.PEP.XML' || $fileFormat eq 'TDM.XML:NoSpec') { # PMF
		#print "$sequence2";
		print "<A href=\"javascript:void(null)\" onmouseover=\"popup('$sequencePosStrg{$sequence}.')\" onmouseout=\"popout()\">$sequence2</A>";
	}
	else {
		print "<A id=\"seq_$rankID\" href=\"javascript:drawSpectrum('seq_$rankID')\" onmouseover=\"popup('Click to view fragmentation spectrum of <B>$sequence</B><BR>$sequencePosStrg{$sequence}.')\" onmouseout=\"popout()\">$sequence2</A>";
	}
	print "$boundaryAA{$ref}[1]";
	if ($varMod) {
		#print "<FONT style=\"font-size:10px\">$varMod</FONT>";
		$varMod=~s/([:.])-1/$1?/g;
		
		# MaxQuant/Spectronaut probabilities #
		my $ptmProbStrg='';
		if($queryPtmProb{$queryID}) {
			my $fileFormat = $dbh->selectrow_array("SELECT NAME,VALID_STATUS,DATA_FILE,FILE_FORMAT FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
			my $softwarePTM = ($fileFormat eq 'MAXQUANT.DIR') ? 'MaxQuant' : ($fileFormat =~ /SPECTRONAUT\.XLS/) ? 'Spectronaut' : 'PtmRS';
			$ptmProbStrg="<B>$softwarePTM probabilities:</B>";
			
			foreach my $modID (sort{lc($varModNames{$a}) cmp lc($varModNames{$b})} keys %{$queryPtmProb{$queryID}}) {
				$ptmProbStrg.="<BR>&bull;$varModNames{$modID}:";
				foreach my $ptmInfo (split(/,/, $queryPtmProb{$queryID}{$modID})) {
					my ($ptmPos, $ptmProb, $ptmBestProb) = split(/:/, $ptmInfo);
					my $posStrg=$convertPos2Text{$ptmPos} || $ptmPos;
					my $ptmBestProbStr = ($ptmBestProb) ? '; Query Grp: '.($ptmBestProb*100).'%&nbsp;' : '';
					$ptmProbStrg.="<BR>&nbsp;&nbsp;-$ptmPos:".&promsMod::PTMProbIcon(($ptmBestProb) ? $ptmBestProb : $ptmProb,{text=>'&nbsp;'.($ptmProb*100).'%&nbsp;'.$ptmBestProbStr,inPopup=>1});
				}
			}
			$ptmProbStrg.="<br/><br/>";
		}
		
		print "<A id=\"vmod_$rankID\" style=\"font-size:10px\" href=\"javascript:editVarModification('vmod_$rankID',",length($sequence),")\" onmouseover=\"popup('${ptmProbStrg}Click to edit modifications <B>position</B>.')\" onmouseout=\"popout()\">$varMod</A>";
		my ($PRS) = ($info =~ /PRS=([^,]+)/);
		if ($PRS){
			&phosphoRS::printIcon($PRS);
		}
	}
	print "</TD>\n</TR></TABLE>\n";
}
$dbh->disconnect;
if (scalar @displayOrder==0){ # no matching peptides after PTMs/rank filtering
	print "<TABLE><TR><TH align=left>&nbsp;&nbsp;<FONT style=\"color:#DD0000;font-style:italic\">No peptides matching selected PTM/Rank filter.</FONT></TH></TR></TABLE>\n";
}

print hidden(-name=>'rankIdList',-default=>"$rankIdList"),"\n";
print "<TABLE><TR><TH align=left>\n";
if ($item ne 'query') {
	my ($selGood,$selPoor); # 1=good 2=poor
	$confLevel=2 if !$confLevel;
	if ($confLevel==1) {$selGood=''; $selPoor='selected';}
	else {$selGood='selected'; $selPoor='';}
	print "&nbsp;&nbsp;Confidence level :<SELECT name='conf_level' style=\"font-weight:bold;\" $disableForm><OPTION value=2 $selGood>Good</OPTION><OPTION value=1 $selPoor>Poor</OPTION></SELECT>\n";
	print hidden(-name=>'oldConf_level',-default=>"$confLevel"),"\n";
}
if ($numSelectable) {
	print "<INPUT type=\"submit\" name=\"save\" value=\"Save selection / Next >>\" style=\"width:190px;\" onclick=\"getSort()\" $disableForm/>\n";
}
if ($item ne 'query' && $msType ne 'PMF') {
# 	print "&nbsp&nbsp&nbsp;or&nbsp&nbsp&nbsp&nbsp<INPUT type=\"button\" value=\"Validate Match Group\" style=\"width:160px;\" onclick=\"validMatchGroup('validate');\" $disableForm/> for interpretations &#8805<INPUT type=\"text\" name=\"MGscore\" value=\"$minScore\" size=\"2\" $disableForm/>\n";
	print "&nbsp;&nbsp;&nbsp;or&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" value=\"Select interpretations\" style=\"width:160px;\" onclick=\"validMatchGroup('validate');\" $disableForm/> with score &#8805<INPUT type=\"text\" name=\"MGscore\" value=\"$minScore\" size=\"2\" $disableForm/>\n";
	print "&nbsp;&nbsp;&nbsp;or&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" value=\"Clear Match Group\" style=\"width:140px;\" onclick=\"validMatchGroup('clear');\" $disableForm/>\n";
}
elsif ($msType eq 'PMF') {
	print "&nbsp;&nbsp;&nbsp;or&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" value=\"Clear Match Group\" style=\"width:140px;\" onclick=\"clearPMFGroup();\" $disableForm/>\n";

}
my $rankHandleString='';
if ($msType ne 'PMF'){
	$rankHandleString = qq
|<TABLE cellpadding=0>
	<TR bgcolor="$color2">
		<TD>&nbsp;Max rank for display <SELECT name="maxDisp" onchange='changeMaxDisp(this);'>|;
	foreach my $rank (1..$allRank) {
		my $isSelected=($rank==$maxRankDisplay)? 'selected' : '';
		$rankHandleString .= "<OPTION $isSelected value=\"$rank\">$rank</OPTION>";
	}
	$rankHandleString .= qq
|		</SELECT>&nbsp;</TD>
		<TD nowrap>&nbsp;Max rank for searching in other analyses <SELECT name="maxSrch" onchange='changeMaxSearch(this);'>|;
	foreach my $rank (1..$allRank) {
		my $isSelected=($rank==$maxRankSearch)? 'selected' : '';
		$rankHandleString .= "<OPTION $isSelected value=\"$rank\">$rank</OPTION>";
	}
	$rankHandleString .= '</SELECT>&nbsp;</TD>';
	if ($item ne 'query') {
		$rankHandleString .= qq
|		<TD nowrap>&nbsp;List only peptides from residues <INPUT type="text" name="rangeMin" value="$pepRangeMin" size="4" $disableForm/> to <INPUT type="text" name="rangeMax" value="$pepRangeMax" size="4" $disableForm/>
		<INPUT type="button" value="Apply" onclick="applyPepRange();" $disableForm/>&nbsp;</TD>|;
	}
	print "</TABLE>\n";
}

my $onLoadString='';
# if ($selectedRankID=~/prot_(\d+)_\d+_0/) { #} rank=0
if (scalar @displayOrder) {
	my @rankData=split(/_/,$selectedRankID);
	if ($rankData[0] eq 'prot') { # rank=0 or access=~/power/
	# 	my $queryID=$1;
		my ($sequence)=($pepInfo{$rankData[1]}[0]=~/SEQ=(\w+)/);
		$onLoadString.="matchProteins('$selectedRankID','$sequence')";
	}
	else { # normal case
		$onLoadString.="drawSpectrum('$selectedRankID')";
	}
}
else {
	$onLoadString.="parent.spectrumFrame.location='$promsPath{html}/nothing.html';";
}
print qq
|</TH></TR>
</TABLE>
$rankHandleString
</FORM>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">
$onLoadString;
setPopup()
</SCRIPT>
</BODY>
</HTML>
|;


##################################<<< SUBROUTINES >>>###########################################

#######################################################
####<<<Routines for sorting interpretations list>>>####
#######################################################
sub sortPeptides {
	my ($queryIDA,$pA)=split(/:/,$a);
	my ($queryIDB,$pB)=split(/:/,$b);
	
	#<PMF
	if ($_[0] eq 'PMF:default') {$absDeltaList{$a}<=>$absDeltaList{$b} || $a cmp $b}
	elsif ($_[0] eq 'PMF:delta') {$absDeltaList{$a}<=>$absDeltaList{$b} || $a cmp $b}
	#elsif ($_[0] eq 'PMF:ppm') {$absPpmList{$a}<=>$absPpmList{$b} || $a cmp $b}
	elsif ($_[0] eq 'PMF:query') {
		$queryNum{$queryIDA}<=>$queryNum{$queryIDB} || $absDeltaList{$a}<=>$absDeltaList{$b} || $a cmp $b
	}
	elsif ($_[0] eq 'PMF:charge') {
		$queryCharge{$queryIDB}<=>$queryCharge{$queryIDA} || $scoreList{$b}<=>$scoreList{$a} || $a cmp $b
	}
	elsif ($_[0] eq 'PMF:sequence') {
		my ($sequenceA)=($pepInfo{$queryIDA}[$pA]=~/SEQ=(\w+)/); my ($sequenceB)=($pepInfo{$queryIDB}[$pB]=~/SEQ=(\w+)/);
		$sequenceA cmp $sequenceB || $absDeltaList{$a}<=>$absDeltaList{$b} || $a cmp $b
	}
	
	#<MIS
	elsif ($_[0] eq 'MIS:default') {$scoreList{$b}<=>$scoreList{$a} || $a cmp $b}
	elsif ($_[0] eq 'MIS:score') {$scoreList{$b}<=>$scoreList{$a} || $a cmp $b}
	#elsif ($_[0] eq 'MIS:delta') {$absDeltaList{$a}<=>$absDeltaList{$b} || $a cmp $b}
	#elsif ($_[0] eq 'MIS:ppm') {$absPpmList{$a}<=>$absPpmList{$b} || $a cmp $b}
	elsif ($_[0] eq 'MIS:query'){
		$queryNum{$queryIDA}<=>$queryNum{$queryIDB} || $scoreList{$b}<=>$scoreList{$a} || $a cmp $b
	}
	elsif ($_[0] eq 'MIS:charge') {
		$queryCharge{$queryIDB}<=>$queryCharge{$queryIDA} || $scoreList{$b}<=>$scoreList{$a} || $a cmp $b
	}
	elsif ($_[0] eq 'MIS:sequence') {
		my ($sequenceA)=($pepInfo{$queryIDA}[$pA]=~/SEQ=(\w+)/); my ($sequenceB)=($pepInfo{$queryIDB}[$pB]=~/SEQ=(\w+)/);
		$sequenceA cmp $sequenceB || $scoreList{$b}<=>$scoreList{$a} || $a cmp $b
	}
}

#########################
####<<<processForm>>>####
#########################
sub processForm {
#print header(-'content-encoding'=>'no',-charset=>'UTF-8'),"\n"; warningsToBrowser(1); # DEBUG
#print "STARTING...<BR>\n";
	my (%newVerifQueryList,%validProtList,%notValidProtList);	# parameters to be transmitted to update validation menu
	my $totalNewPeptides=0;										# parameters to be transmitted to update validation menu
	my %queryStatus; # must be hash because same queryID more than once (eg if query view)!
	my %listIdentifiers; # list of identifiers affected
	my @rankIdList=split(/:/,param('rankIdList')); # list of selectable peptides in form (:qID_qNum_rank)
	my %allRankInfo; # contains info for all peptides already extracted from DB (prevents accessing DB more than once for same data)
	my $pepSelForm; # number of selected peptides in form
	my $formConfLevel=param('conf_level'); # defined only if item is a protein
	my $changeConfLevel;
	if ($formConfLevel) {
		$changeConfLevel=($formConfLevel eq param('oldConf_level'))? 0 : 1;
	}

	####<List of DB queries>####
	my $sthSEL1=$dbh->prepare("SELECT IDENTIFIER FROM RANK_PROTEIN_MATCH WHERE QUERY_NUM=? AND PEP_RANK=? AND ID_ANALYSIS=$analysisID");
	my %sthUpRank;
	my $infoString='';
	foreach my $rank (1..$maxRankDisplay) {
		$sthUpRank{$rank}=$dbh->prepare("UPDATE QUERY_VALIDATION SET INFO_PEP$rank=?,VALID_STATUS=? WHERE ID_QUERY=?");
		$infoString.="INFO_PEP$rank";
		$infoString.=',' if $rank<$maxRankDisplay;
	}
	my $sthSelInfo=$dbh->prepare("SELECT $infoString FROM QUERY_VALIDATION WHERE ID_QUERY=?");

	####<Looping through all matching queries>####
	foreach my $rankID (@rankIdList) {
		my ($queryID,$queryNum,$rank)=split(/_/,$rankID);

		my $rankInfo=param("pepInfo_$rankID");
		my ($oldSelect)=($rankInfo=~/SEL=(-*\d)/); # cannot be -1 but can be -2,-3

		$queryStatus{$queryID}=param("status_$queryID") unless defined($queryStatus{$queryID}); # defined because can be 0
#! 		$queryStatus{$queryID}=0 if $queryStatus{$queryID}==-1; # no previous validations (cannot be -3 [-2!] -> SEL=-1)
		#                      ^ not 0 here

		my $select;
		my $updateQuery=0;
		my $computeStatus=0;
		my $newPeptide=0;
		if (param("sel_$rankID")) { # select box is checked
			$select=2;
			if ($oldSelect < 2) { # change in query selection status (cannot be -1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!***********
				$updateQuery=($oldSelect==1)? 1 : 2; # 1: update rank only, 2: update all
				if ($queryStatus{$queryID}==-1 || $queryStatus{$queryID}==-2) {$queryStatus{$queryID}=1;}
				else {$queryStatus{$queryID}++;} # was >=0 before
				$newPeptide=1;
			}
			$pepSelForm++; # needed for confLevel
		}
		elsif (param("rej_$rankID")) { # reject box is checked
			$select=-3;
			if ($oldSelect > -3) { # change in selection status (cannot be -1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!**********
				$updateQuery=($oldSelect==-2)? 1 : 2; # 1: update rank only, 2: update all
				if ($queryStatus{$queryID}==-1) {$queryStatus{$queryID}=-2;} # 1st reject
				elsif ($oldSelect >= 1) {
					$computeStatus=1;
					# $queryStatus{$queryID}--;
					$newPeptide=-1;
				}
				else {$computeStatus=1;} # queryStatus was -2 (previous reject)
			}
		} # else => neither the select box nor the reject box is checked (never-verified rank)

		###<Comments>###
		my ($oldPepComments)=($rankInfo=~/COM=(.+)/); $oldPepComments='' unless $oldPepComments;
		my $newPepComments=param("pepCom_$rankID"); $newPepComments='' unless $newPepComments;
		$rankInfo=~s/COM=.*//; # if $oldPepComments; # clearing previous comments
		$rankInfo.="COM=$newPepComments" if $newPepComments;
		$updateQuery=1 if $newPepComments ne $oldPepComments;
#print "$rankID: Old=$oldPepComments New=$newPepComments Info=$rankInfo<BR>\n";

		if ($computeStatus) { # Select must be performed to know new valid_status
			$sthSelInfo->execute($queryID);
			my @pepInfoList=$sthSelInfo->fetchrow_array;
			my %curPepInfo;
			my $i=0;
			foreach my $rk (1..$maxRankDisplay) {
				$curPepInfo{$rk}=$pepInfoList[$i];
				$i++;
			}
			my ($totRank,$selRank,$rejRank)=(0,0,0);
			foreach my $rk (1..$maxRankDisplay) {
				next unless $curPepInfo{$rk};
				my $sel;
				if ($rk==$rank) {$sel=$select;} # use new valid_status, not old one
				else {($sel)=($curPepInfo{$rk}=~/SEL=(-*\d)/);}
				$totRank++ if $sel != -1; # better score exists
				if ($sel>=1) {$selRank++;}
				elsif ($sel<=-2) {$rejRank++;}
			}
			if ($selRank) {$queryStatus{$queryID}=$selRank;}
			elsif ($rejRank==$totRank) {$queryStatus{$queryID}=0;}
			else {$queryStatus{$queryID}=-2;} # rejected + not verified
		}

		###<Updating table QUERY_VALIDATION>###
		if ($updateQuery) {
			$newVerifQueryList{$queryID}=1 if param("status_$queryID")==-1; # A new query is being verified, hash prevents duplicates
			$rankInfo=~s/SEL=-*\d/SEL=$select/; # starting SEL cannot be -1
# 				my $requete="UPDATE QUERY_VALIDATION SET INFO_PEP$rank='$rankInfo',VALID_STATUS=$queryStatus{$queryID} WHERE ID_QUERY=$queryID";
# 			 	$dbh->do($requete) || die $dbh->errstr;
			$sthUpRank{$rank}->execute($rankInfo,$queryStatus{$queryID},$queryID) || die $dbh->errstr;
		}

		###<Updating table PROTEIN_VALIDATION>###
		if ($updateQuery==2 || $changeConfLevel) { # protein update necessary
			###<Finding identifier for matched proteins from RANK_PROTEIN_MATCH table
			$sthSEL1->execute($queryNum,$rank) || die $sthSEL1->errstr;
			while (my ($identifier) = $sthSEL1->fetchrow_array) {
				$listIdentifiers{$identifier}=1;
			}
			$totalNewPeptides+=$newPeptide; # 0 if !updateQuery
		}

		$allRankInfo{"$queryNum:$rank"}=$rankInfo; # no need for later retreval from DB
	}
	$sthSEL1->finish;
	$sthSelInfo->finish;
	foreach my $rank (keys %sthUpRank) {$sthUpRank{$rank}->finish;}

	####<Protein validation>####
	my $newValidProteins=&updateProtValidation(\%listIdentifiers,\%allRankInfo,$formConfLevel,$pepSelForm);

	####<Finding next unselected item>####
	my ($newItemID,$newPage);
	if ($item eq 'query') {

		###<Searching in current page 1st>###
		my ($sortedItem,$order)=split(/:/,$selectedSort);
		my ($lowQueryID,$upQueryID)=split(/-/,$rangeStrg);
		my ($firstID,$firstNotVerifID);
		my $sthQ1=$dbh->prepare("SELECT ID_QUERY,VALID_STATUS FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID AND ID_QUERY>=$lowQueryID AND ID_QUERY<=$upQueryID ORDER BY $sortedItem $order,ID_QUERY ASC");
		$sthQ1->execute || die $sthQ1->errstr;
		while (my($queryID,$validStatus)=$sthQ1->fetchrow_array) {
			$firstID=$queryID if (!$firstID && $maxPage==1);
			if ($validStatus==-1) {
				$newItemID=$queryID;
				last;
			}
		}
		$sthQ1->finish;
		if ($newItemID) {$newPage=$page;}
		elsif ($maxPage==1) {
			$newItemID=$firstID;
			$newPage=1;
		}
		else {
			###<Searching in all pages>###
			my $sthQ2=$dbh->prepare("SELECT ID_QUERY,VALID_STATUS FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID ORDER BY ID_QUERY ASC");
			$sthQ2->execute || die $sthQ2->errstr;
			my $refQueryData=$sthQ2->fetchall_arrayref; # reference to an array
			$sthQ2->finish;

			my $queryPageSize=&promsConfig::getQueryPageSize;
			my $queryCount=0;
			my $pageCount=1;
			foreach my $refQuery (@{$refQueryData}) { # ascending query id/num
				$queryCount++;
				$firstID=$refQuery->[0] unless $firstID;
				if ($refQuery->[1]==-1 && !$newItemID) {
					$newItemID=$refQuery->[0];
					$newPage=$pageCount;
					last;
				}
				if ($queryCount >= $queryPageSize) {
					$queryCount=0;
					$pageCount++;
				}
			}
			unless ($newItemID) { # if all checked: select first one
				$newItemID=$firstID;
				$newPage=1;
			}
		}
	}
	elsif ($item eq 'protein') {
		($newItemID,$newPage)=&getNextProtID;
	}
	else {($newItemID,$newPage)=&getMyNextProtID;} # myList

	$dbh->commit;
	&promsMod::updateAnalysisHistory($dbh,$analysisID,'',"manual");
	$dbh->disconnect;

	my $newVerifQueries=scalar keys %newVerifQueryList;
	my $newValidDecoy=($newVerifQueries || $totalNewPeptides || $newValidProteins)? -2 : 0;

	##############
	####<HTML>####
	##############
	print header(-charset=>'UTF-8'); warningsToBrowser(1);
	my $updateRank=($item eq 'myList')? 2 : 1;
	print qq
|<HTML><HEAD>
<SCRIPT LANGUAGE="JavaScript">
if (top.protWindow && !top.protWindow.closed) {top.protWindow.location.reload(true);} // in case changes in validation affect protein being displayed in protWindow
parent.updateValidation($newItemID,$newVerifQueries,$totalNewPeptides,$newValidDecoy,$newValidProteins,'$item',$newPage,$updateRank);
</SCRIPT>
</HEAD></HTML>
|;
	exit;
}


########################################
####<<<Validate/Clear Match Group>>>####
######################################## Globals: $analysisID,$itemID,$item
sub validateMatchGroup {
	my $action=param('valMG');
	my $minMGscore=param('MGscore');
	my $totalNewPeptides=0;
	my $newVerifQueries=0;
	my $newValidProteins=0;

	####<List of proteins in Match Group & matching peptides>####
	my %listPeptides;
	my ($matchGroup)=$dbh->selectrow_array("SELECT MATCH_GROUP FROM PROTEIN_VALIDATION WHERE ID_PROT_VALID=$itemID");
	my $sthSelMG=$dbh->prepare("SELECT IDENTIFIER FROM PROTEIN_VALIDATION WHERE MATCH_GROUP=$matchGroup AND ID_ANALYSIS=$analysisID");
	my $sthSelQR=$dbh->prepare("SELECT QUERY_NUM,PEP_RANK FROM RANK_PROTEIN_MATCH WHERE IDENTIFIER=? AND ID_ANALYSIS=$analysisID");
	$sthSelMG->execute;
	while(my ($identifier)=$sthSelMG->fetchrow_array) {
		$sthSelQR->execute($identifier);
		while (my($queryNum,$rank)=$sthSelQR->fetchrow_array) {
			$listPeptides{$queryNum}{$rank}=1;
		}
	}
	$sthSelMG->finish;
	$sthSelQR->finish;

	####<Peptide validation>####
	my %allRankInfo;
	my (%sthSelRank,%sthUpRank);
	foreach my $rank (1..$allRank) {
		$sthSelRank{$rank}=$dbh->prepare("SELECT ID_QUERY,VALID_STATUS,INFO_PEP$rank FROM QUERY_VALIDATION WHERE QUERY_NUM=? AND ID_ANALYSIS=$analysisID");
		$sthUpRank{$rank}=$dbh->prepare("UPDATE QUERY_VALIDATION SET INFO_PEP$rank=?,VALID_STATUS=? WHERE ID_QUERY=?");
	}
	if ($action eq 'validate') {
		QUERY:foreach my $queryNum (keys %listPeptides) {
			foreach my $rank (sort{$a<=>$b} keys %{$listPeptides{$queryNum}}) {
				$sthSelRank{$rank}->execute($queryNum);
				my ($queryID,$validStatus,$rankInfo)=$sthSelRank{$rank}->fetchrow_array;
				my ($oldSelect)=($rankInfo=~/SEL=(-*\d)/);
				#unless ($rankInfo=~/SEL=[01]/) { #} check only never validated (no overwrite)
				if ($oldSelect==-3 || $oldSelect==-1 || $oldSelect==2) { # -3: do not overwrite manual rejection
					$allRankInfo{"$queryNum:$rank"}=$rankInfo; # no need for later retrieval from DB
					next;
				}
				my ($pepScore)=($rankInfo=~/SC=(\d+\.*\d*)/);
				if ($pepScore >= $minMGscore) {
					$rankInfo=~s/SEL=$oldSelect/SEL=2/;
					$allRankInfo{"$queryNum:$rank"}=$rankInfo; # no need for later retreval from DB
					my $queryStatus=($validStatus==-1)? 1 : ($validStatus>=0 && $oldSelect<=0)? $validStatus+1 : $validStatus;
					$sthUpRank{$rank}->execute($rankInfo,$queryStatus,$queryID);
					$totalNewPeptides++ if $oldSelect<=0;
					$newVerifQueries++ if ($validStatus==-1 || $validStatus==0);
					next QUERY; # select only 1 rank/query (best)
				}
			}
		}
	}
	else { # action = clear
		foreach my $queryNum (keys %listPeptides) {
			my $selPep=0;
			foreach my $rank (sort{$a<=>$b} keys %{$listPeptides{$queryNum}}) {
				$sthSelRank{$rank}->execute($queryNum);
				my ($queryID,$validStatus,$rankInfo)=$sthSelRank{$rank}->fetchrow_array;
				my ($oldSelect)=($rankInfo=~/SEL=(-*\d)/);
				if ($oldSelect <= -2 || $oldSelect >= 1) {
					my $queryStatus=-1;
					$rankInfo=~s/SEL=$oldSelect/SEL=0/;
					if ($oldSelect>=1) {
						$totalNewPeptides--;
						$selPep=1;
					}
					$sthUpRank{$rank}->execute($rankInfo,$queryStatus,$queryID);
				}
				$allRankInfo{"$queryNum:$rank"}=$rankInfo; # no need for later retrieval from DB
			}
			$newVerifQueries-- if $selPep;
		}
	}
	foreach my $rank (1..$allRank) {
		$sthSelRank{$rank}->finish;
		$sthUpRank{$rank}->finish;
	}

	####<List of all matched proteins>####
	my %listIdentifiers;
	my $sthSelProt=$dbh->prepare("SELECT IDENTIFIER FROM RANK_PROTEIN_MATCH WHERE QUERY_NUM=? AND PEP_RANK=? AND ID_ANALYSIS=$analysisID");
	foreach my $queryNum (keys %listPeptides) {
		foreach my $rank (keys %{$listPeptides{$queryNum}}) {
			$sthSelProt->execute($queryNum,$rank);
			while (my($identifier)=$sthSelProt->fetchrow_array) {
				$listIdentifiers{$identifier}=1;
			}
		}
	}
	$sthSelProt->finish;

	####<Protein validation>####
	$newValidProteins=&updateProtValidation(\%listIdentifiers,\%allRankInfo,0,0);

	####<Finding next selected item>####
	my ($newItemID,$newPage);
	if ($action eq 'validate') {
		if ($item eq 'protein') {($newItemID,$newPage)=&getNextProtID;}
		else {($newItemID,$newPage)=&getMyNextProtID;} # myList
	}
	else {
		$newItemID=$itemID;
		$newPage=$page;
	}

	$dbh->commit;
	$dbh->disconnect;

	my $newValidDecoy=($newVerifQueries || $totalNewPeptides || $newValidProteins)? -2 : 0;

	##############
	####<HTML>####
	##############
	my $updateRank=($item eq 'myList')? 2 : 1;
	print header(-charset=>'UTF-8');
	print qq
|<HTML><HEAD>
<SCRIPT LANGUAGE="JavaScript">
if (top.protWindow && !top.protWindow.closed) {top.protWindow.location.reload(true);} // in case changes in validation affect protein being displayed in protWindow
parent.updateValidation($newItemID,$newVerifQueries,$totalNewPeptides,$newValidDecoy,$newValidProteins,'$item',$newPage,$updateRank);
</SCRIPT>
</HEAD></HTML>
|;
	exit;
}

#######################################
####<<<Update protein validation>>>#### subroutine common to processForm & validateMatchGroup
####################################### Globals : $analysisID,$dbh
sub updateProtValidation {
	my ($refListIdentifiers,$refAllRankInfo,$formConfLevel,$pepSelForm)=@_;

	####<Updating info for all matched proteins>####
	my $sthSEL2=$dbh->prepare("SELECT ID_PROT_VALID,SEL_STATUS,MAX_MATCH,CONF_LEVEL FROM PROTEIN_VALIDATION WHERE IDENTIFIER=? AND ID_ANALYSIS=$analysisID");
	my $sthSEL3=$dbh->prepare("SELECT QUERY_NUM,PEP_RANK,MATCH_MULTI FROM RANK_PROTEIN_MATCH WHERE IDENTIFIER=? AND ID_ANALYSIS=$analysisID");
	my $scoreString=($msType eq 'PMF')? '' : ',SCORE=?';
	my $sthUP2=$dbh->prepare("UPDATE PROTEIN_VALIDATION SET NUM_MATCH=?,SEL_STATUS=?$scoreString,CONF_LEVEL=? WHERE ID_PROT_VALID=?");
	my %sthSelRank;
	foreach my $rank (1..$allRank) {
		$sthSelRank{$rank}=$dbh->prepare("SELECT INFO_PEP$rank FROM QUERY_VALIDATION WHERE QUERY_NUM=? AND ID_ANALYSIS=$analysisID");
	}

	my $newValidProt=0; # number of new selected proteins
	foreach my $identifier (keys %{$refListIdentifiers}) {

		###<Fetching protein data from table PROTEIN_VALIDATION>###
		$sthSEL2->execute($identifier) || die $sthSEL2->errstr;
		my ($protValidID,$oldSelStatus,$maxMatch,$confLevel)=$sthSEL2->fetchrow_array;
		next if $oldSelStatus <= -2; # protein is excluded/filtered from list => do not update selStatus

		###<Fetching info on all matching peptides>###
		my $selStatus=$oldSelStatus; # default
		my $numMatch=0; # reset to 0
		my $protScore=0; # reset to 0
		my $rejectedMatch=0;
		$sthSEL3->execute($identifier) || die $sthSEL3->errstr;
		while (my($queryNum,$rank,$matchFreq)=$sthSEL3->fetchrow_array) {
			my $qrKey="$queryNum:$rank";
			unless ($refAllRankInfo->{$qrKey}) { # rankInfo not already extracted from DB
# 				($refAllRankInfo->{$qrKey})=$dbh->selectrow_array("SELECT INFO_PEP$rank FROM QUERY_VALIDATION WHERE QUERY_NUM=$queryNum AND ID_ANALYSIS=$analysisID");
				next unless ($sthSelRank{$rank}) ; #mauvaise gestion bug, mais evite plantage FP ;
				$sthSelRank{$rank}->execute($queryNum);
				($refAllRankInfo->{$qrKey})=$sthSelRank{$rank}->fetchrow_array;
			}
			if ($refAllRankInfo->{$qrKey} =~ /SEL=[12]/) {
				$numMatch++;
				my ($pepScore)=($refAllRankInfo->{$qrKey}=~/SC=(\d+\.*\d*)/); # Score corresponding to peptide
				$protScore+=($matchFreq*$pepScore);
			}
			elsif ($refAllRankInfo->{$qrKey} =~ /SEL=-[23]/) {$rejectedMatch++;}
		}

		if ($numMatch) { # at least 1 peptide selected
			$selStatus=($numMatch+$rejectedMatch==$maxMatch)? 2 : 1;
		}
		else { # no peptides selected
			$selStatus=($rejectedMatch==$maxMatch)? 0 : -1;
		}

		###<Updating number of selected proteins for main validation window>###
		if ($selStatus > 0 && $oldSelStatus <= 0) {$newValidProt++;}
		elsif ($selStatus <= 0 && $oldSelStatus > 0) {$newValidProt--;}

		###<Updating protein confidence level>###
		$confLevel=2 unless $confLevel;
		if ($formConfLevel) { # test even if not changed because new proteins can be matched by new checked peptides
			if ($pepSelForm==$numMatch || $confLevel < $formConfLevel) { # all matching peptides are in form || better conf_level in form
				$confLevel=$formConfLevel;
			}# else: confLevel=2(Good) => do not change
		}

		###<Updating protein entries in table PROTEIN_VALIDATION>###
		if ($msType eq 'PMF') {$sthUP2->execute($numMatch,$selStatus,$confLevel,$protValidID) || die $sthUP2->errstr;}
		else {$sthUP2->execute($numMatch,$selStatus,$protScore,$confLevel,$protValidID) || die $sthUP2->errstr;}
	}

	$sthSEL2->finish;
	$sthSEL3->finish;
	$sthUP2->finish;
	foreach my $rank (1..$allRank) {$sthSelRank{$rank}->finish;}

	return $newValidProt;

}


###################################
####<<<Change Protein Status>>>#### Changes protein sel_status: to -2 (exclude) or to 0/1 (include)
################################### Globals: $analysisID,$itemID,$item
sub changeProteinStatus {

# print header(-'content-encoding'=>'no',-charset=>'UTF-8'),"DEBUG changeProteinStatus:<BR>\n"; warningsToBrowser(1); # DEBUG

	my $identifier=param('IDENTIF');
	my $oldStatus=param('STATUS');

	my ($selStatus,$newValidProteins,$newItemID,$newPage);

	####<Excluding protein>####
	if ($oldStatus > -2) {
		$selStatus=-2;
		$newValidProteins=($oldStatus>0)? -1 : 0;

		###<Updating table with new status>###
		my $scoreString=($msType eq 'PMF')? '' : ',SCORE=0';
		$dbh->do("UPDATE PROTEIN_VALIDATION SET SEL_STATUS=$selStatus,NUM_MATCH=0 $scoreString WHERE ID_PROT_VALID=$itemID") || die $dbh->errstr;

		###<Finding next selected protein ID>###
		if ($item eq 'protein') {($newItemID,$newPage)=&getNextProtID;}
		else {($newItemID,$newPage)=&getMyNextProtID;} # myList
	}

	####<Including protein>####
	else {
		###<Fetching matching peptides and their selection status>###
		my ($numMatch,$protScore,$rejectedMatch)=(0,0,0);
		my $maxMatch;
		my $sthSEL1=$dbh->prepare("SELECT QUERY_NUM,PEP_RANK,MATCH_MULTI FROM RANK_PROTEIN_MATCH WHERE IDENTIFIER='$identifier' AND ID_ANALYSIS=$analysisID");
		$sthSEL1->execute;
		while (my ($queryNum,$rank,$matchFreq) = $sthSEL1->fetchrow_array) {
			my ($rankInfo)=$dbh->selectrow_array("SELECT INFO_PEP$rank FROM QUERY_VALIDATION WHERE QUERY_NUM=$queryNum AND ID_ANALYSIS=$analysisID");
			if ($rankInfo =~ /SEL=[12]/) {
				$numMatch++;
				my ($pepScore)=($rankInfo=~/SC=(\d+\.*\d*)/); # Score corresponding to peptide
				$protScore+=($matchFreq*$pepScore);
			}
			elsif ($rankInfo =~ /SEL=-[23]/) {$rejectedMatch++;}
			$maxMatch++ unless $rankInfo =~ /SEL=-1/;
		}
		if ($numMatch) { # at least 1 peptide selected
			$selStatus=($numMatch+$rejectedMatch==$maxMatch)? 2 : 1;
		}
		else { # no peptides selected
			$selStatus=($rejectedMatch==$maxMatch)? 0 : -1;
		}

		###<Updating table with new status & score>###
		my $scoreString=($msType eq 'PMF')? '' : ",SCORE=$protScore";
		$dbh->do("UPDATE PROTEIN_VALIDATION SET NUM_MATCH=$numMatch,SEL_STATUS=$selStatus $scoreString,CONF_LEVEL=2 WHERE ID_PROT_VALID=$itemID") || die $dbh->errstr;

		$newItemID=$itemID;
		$newPage=$page;
		$newValidProteins=($selStatus>0)? 1 : 0;
	}

	&promsMod::updateAnalysisHistory($dbh,param('ANAID'),'','incexc');

	$dbh->commit;
	$dbh->disconnect;

	my $newDecoyPeptides=($newValidProteins)? -2 : 0;

	##############
	####<HTML>####
	##############
	my $updateRank=($item eq 'myList')? 2 : 1;
	print header(-charset=>'UTF-8');
	print qq
|<HTML>
<HEAD>
<SCRIPT LANGUAGE="JavaScript">
if (top.protWindow && !top.protWindow.closed && top.protWindow.selectedIdentifier=='$identifier') {top.protWindow.location.reload(true);} // refreshing change in protein status being displayed in protWindow
|;
	print "parent.selectableProteins++;\n" if $selStatus > -2;
	print "parent.updateValidation($newItemID,0,0,$newDecoyPeptides,$newValidProteins,'$item',$newPage,$updateRank);\n";
	print "</SCRIPT>\n</HEAD>\n</HTML>\n";
	exit;
}

############################################ subroutine common to processForm & validateMatchGroup
####<<<Get ID of next protein in list>>>#### Finds the next protein ID to be selected by listTempProteins.cgi
############################################ Globals: $analysisID,$itemID,$selectedSort
sub getNextProtID {
	my ($newProtID,$newPage);

	####<Query 1>####
	my @sortRules=split(/,/,$selectedSort);
	my $sortString='';
	foreach my $rule (@sortRules) {
		my ($sortedItem,$order)=split(/:/,$rule);
		$sortString.="$sortedItem $order,";
	}
	$sortString.='ID_PROT_VALID ASC';
	my $sth1=$dbh->prepare ("SELECT ID_PROT_VALID,SEL_STATUS,MATCH_GROUP FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$analysisID AND IDENTIFIER NOT LIKE 'DECOY_%' ORDER BY $sortString");
	$sth1->execute || die $sth1->errstr;
	my $refData1=$sth1->fetchall_arrayref; # reference to an array
	$sth1->finish;

	###################################
	####<Searching in current page>####
	###################################
	####<Converting match group string to list>####
	my (%matchGroups,%currentMGr);
	my @tmpList=split(',',$rangeStrg);
	foreach my $subStrg (@tmpList) {
		if ($subStrg=~/-/) {
			my ($start,$end)=split(/-/,$subStrg);
			foreach my $mg ($start..$end) {$currentMGr{$mg}=1};
		}
		else {$currentMGr{$subStrg}=1;}
	}
	####<Scanning>####
	my ($firstProtID,$reachedItem);
	foreach my $refProt (@{$refData1}) {
		$matchGroups{$refProt->[2]}++; # needed for all-page search
		next unless $currentMGr{$refProt->[2]};
		$firstProtID=$refProt->[0] if (!$firstProtID && ($showFiltered || (!$showFiltered && $refProt->[1] > -3)));
		if ($refProt->[0]==$itemID) {$reachedItem=1; next;} # skip everything before $itemID
		next unless $reachedItem;
		if ($refProt->[1]==-1) {
			$newProtID=$refProt->[0];
			last;
		}
	}
	if ($newProtID) {return ($newProtID,$page);}
	elsif ($maxPage==1) {return ($firstProtID,1);}

	################################
	####<Searching in all pages>#### if no  unverified protein found in current page
	################################
	my $pageSize=&promsConfig::getProteinPageSize;

	####<Query 2>####
	my $sth2=$dbh->prepare ("SELECT ID_PROT_VALID,SEL_STATUS,MATCH_GROUP FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$analysisID AND IDENTIFIER NOT LIKE 'DECOY_%' ORDER BY MAX_SCORE DESC,ID_PROT_VALID ASC");
	$sth2->execute || die $sth2->errstr;
	my $refData2=$sth2->fetchall_arrayref; # reference to an array
	$sth2->finish;

	####<Splitting proteins in multiple sets>####
	my $protCount=0;
	my $pageIndex=0;
	my @pageMGr;
	my %usedMGr;
	foreach my $refProt (@{$refData2}) {
		next if $usedMGr{$refProt->[2]};
		$usedMGr{$refProt->[2]}=1;
		$protCount+=$matchGroups{$refProt->[2]};
		$pageMGr[$pageIndex]{$refProt->[2]}=1;
		if ($protCount >= $pageSize) {
			$protCount=0;
			$pageIndex++;
		}
	}

	###>Scanning protein sets for 1st selectable protein<###
	$firstProtID=undef;
	my ($readFirst,$goodPage);
	my $lastPgIdx=$page-2; $lastPgIdx=0 if $lastPgIdx < 0;
	PAGE:foreach my $pgIndex ($page..$#pageMGr,0..$lastPgIdx) { # current index=page-1! (end of list, then beg, skip current page)
		$readFirst=1 if $pgIndex==0; # 1st page has been read
		foreach my $refProt (@{$refData2}) {
			if ($readFirst && !$firstProtID && ($showFiltered || (!$showFiltered && $refProt->[1] > -3))) {
				$firstProtID=$refProt->[0];
				$goodPage=$pgIndex+1;
			}
			next unless $pageMGr[$pgIndex]{$refProt->[2]}; # restrict to selected MGr
			if ($refProt->[1]==-1 && !$newProtID) {
				$newProtID=$refProt->[0];
				$newPage=$pgIndex+1;
				last PAGE;
			}
		}
	}
	unless ($newProtID) { # if all checked: select first one
		$newProtID=$firstProtID;
		$newPage=$goodPage;
	}
	return ($newProtID,$newPage);
}

sub getMyNextProtID {
	my ($newProtID,$newPage);

	####<Fetching private list of proteins>####
	my %myProteinList;
	my %sortKeyPos=(
		'IDENTIFIER'=>0,
		'SEL_STATUS'=>1,
		'NUM_MATCH'=>2,
		'MAX_MATCH'=>3,
		'MATCH_GROUP'=>4,
		'SCORE'=>5,
		'MAX_SCORE'=>6
	);
	my $sthD=$dbh->prepare("SELECT IDENTIFIER,SEL_STATUS,NUM_MATCH,MAX_MATCH,MATCH_GROUP,SCORE,MAX_SCORE FROM PROTEIN_VALIDATION WHERE ID_PROT_VALID=?");
	my $myProtFile="$promsPath{valid}/ana_$analysisID/my_proteins.txt";
	open (MYPROT,$myProtFile);
	while (<MYPROT>) {
		my ($protID,$page)=($_=~/^\S+\t(\d+)\t(\d+)/);
		next unless $protID; # not in analysis
		$sthD->execute($protID);
		my @data=$sthD->fetchrow_array;
		next if (!$showFiltered && $data[1]<=-3);
		@{$myProteinList{$protID}}=(@data,$page);
	}
	close MYPROT;

	$sthD->finish;

	my @sortRules=split(',',$selectedSort);
	my ($sortKey,$sortOrder)=split(':',$sortRules[0]);

	my ($firstProtID); # 1st backup selection (1st unverified protein)
	foreach my $protID (sort{&sortProt(\%myProteinList,$sortKeyPos{$sortKey},$sortOrder)} keys %myProteinList) {
		#<data for 1st protein in list
		unless ($firstProtID) {
			$firstProtID=$protID;
		}
		#<data for 1st unverified protein in list
		if ($myProteinList{$protID}[1]==-1 && !$newProtID) {
			$newProtID=$protID;
			last;
		}
	}
	$newProtID=$firstProtID unless $newProtID; # default selected protein is not in list
	$newPage=$myProteinList{$newProtID}[-1];

	return ($newProtID,$newPage);
}

############## Protein sort subroutines ###############
sub sortProt {
	my ($refProtList,$sortPos,$sortOrder)=@_;
	# Match group. (MG asc> MaxScore desc > Identifier asc)
	if ($sortPos==4) {$refProtList->{$a}[4]<=>$refProtList->{$b}[4] || $refProtList->{$b}[6]<=>$refProtList->{$a}[6] || lc($refProtList->{$a}[0]) cmp lc($refProtList->{$b}[0])}
	# Distribution asc. (Item asc > Identifier asc)
	elsif ($sortOrder eq 'ASC') {$refProtList->{$a}[$sortPos]<=>$refProtList->{$b}[$sortPos] || lc($refProtList->{$a}[0]) cmp lc($refProtList->{$b}[0])}
	# Distribution desc. (Item desc > Identifier asc)
	else {$refProtList->{$b}[$sortPos]<=>$refProtList->{$a}[$sortPos] || lc($refProtList->{$a}[0]) cmp lc($refProtList->{$b}[0])}
}

sub updateWaitBox {
	my $text=$_[0] || '.';
	print "<SCRIPT language=\"JavaScript\">document.getElementById('waitSpan').innerHTML+='$text';</SCRIPT>\n";
}
sub printPRSIcon{
	my $color = shift;
	print "<FONT style=\"background-color:$color\">PRS</FONT>";
}

####>Revision history<####
# 2.3.11 [BUGFIX] Fixed non-PRS PTM probability detection (PP 22/04/21)
# 2.3.10 [BUGFIX] Fixed forgotten test for undefined $ptmProb (PP 01/10/20)
# 2.3.9 [BUGFIX] Fixed displaying of all peptides instead of modified ones (VS 09/09/20)
# 2.3.8 [ENHANCEMENT] Added display of PtmRS site localization probabilities (VS 21/08/20)
# 2.3.7 [ENHANCEMENT] Added display of ion/query charge (VS 30/03/20)
# 2.3.6 Update PTM filter so as to keep fix modifications (GA 17/08/17)
# 2.3.5 Minor modif to disable spectrum visualization for TDM.PEP.XML -> like PMF (MLP 16/12/16)
# 2.3.4 Minor change to put ppm by default for $DeltaType (GA 13/04/16)
# 2.3.3 Minor change for ALTNAMES of modifications (GA 06/03/14)
# 2.3.2 Minor change for DELT pattern grouping in order to get Paragon data (GA 03/12/13)
# 2.3.1 Speed/code optimization & optional multi-analysis scan (PP 26/09/13)
# 2.3.0 Calling updateAnalysisHistory if user activate a low-score peptide (FY 17/09/13)
# 2.2.9 Removing VAR_MOD from script (GA 31/05/13)
# 2.2.8 Updates Analysis lower scoring peptide activation flag (PP 29/03/13)
# 2.2.7 Add Substitution visualization (GA 21/03/13)
# 2.2.6 Multi-databank flag on proteins (PP 11/12/12)
# 2.2.5 Minor modification: uri_escape for identifiers that can contain '#' (GA 28/11/12)
# 2.2.4 Displays Qvality PEP & q-value in a score popup (PP 07/08/12)
# 2.2.3 Fixed bug row color not changed when a "parent" peptide is filtered out but not its hidden "child" (PP 30/05/12)
# 2.2.2 Added varMod comparison to multi-analysis check: same positions required (PP 24/04/12)
# 2.2.1 Set delta mass allowed for multi-analysis to +/-0.25 on mass exp (PP 13/04/12)
# 2.2.0 Display PRS status on phosphopeptides (FY 27/03/12)
# 2.1.9 Manager has full access validation data (PP 01/10/11)
# 2.1.8 header(-'content-encoding'=>'no') (PP 01/07/11)
# 2.1.7 Compatibility with mixed (SQ) search (PP 16/04/2011)
# 2.1.6 Minor changes for bigger popup varMod edition (PP 11/04/2011)
# 2.1.5 Added minor changes (2.1.3m) made on production server (PP 11/04/2011)
# 2.1.4 Adding updating analysis history (FY 02/2011)
# 2.1.3 Minor update to print well the elution time (GA 07/01/11)<BR>See 2.6.5 modification in storeAnalyses.cgi
# 2.1.2 Correct bug (since 2.0.8) check matched Ana for activated lower-scoring peptides<BR>Correct all charges count in matched Ana (PP 16/12/2010)
