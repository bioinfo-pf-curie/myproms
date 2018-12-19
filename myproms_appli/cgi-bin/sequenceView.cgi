#!/usr/local/bin/perl -w

################################################################################
# sequenceView.cgi     3.2.9	                                               #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Displays detailed information on a protein:                                  #
#  sequence, peptides, analyses where found, ...                               #
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
use promsConfig;
use promsMod;
use strict;
use phosphoRS;
use List::Util qw(max);
use LWP::UserAgent;

#print header(-'content-encoding'=>'no'); warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my %msTypeNames=&promsConfig::getMsType;
my $userID=$ENV{'REMOTE_USER'};
my ($color1,$color2)=&promsConfig::getRowColors;
my %convertPos2Text=('-'=>'Protein N-term','='=>'Any N-term','+'=>'Protein C-term','*'=>'Any C-term');

####################
#    Parameters    #
####################
my $proteinID=&promsMod::cleanNumericalParameters(param('id_prot'));
if (param('AJAX')) {
	if (param('AJAX') eq 'connProt') {&ajaxUpdateConnectedProteins;}
	elsif (param('AJAX') eq 'vmodPep') {&ajaxDisplayVarModPeptides;}
	elsif (param('AJAX') eq 'searchInteract') {&ajaxSearchInteractors;}
	elsif (param('AJAX') eq 'quantifList') {&ajaxGetQuantificationList;}
	exit;
}
#my $frame=(param('frame'))? param('frame') : 'top';
my $idType=(param('id_type'))? param('id_type') : 'valid'; # (or notValid) type of type provided for proteinID
#my @analysesID = split /,/, param('id_ana');
#my $analysisID = $analysesID[0];
#my $simGr=param('simGr'); # script is called from Clustering mode
my $showMSdata=param('msdata');
#my $call=(param('call'))? param('call') : 'list'; # script is called from valid protein 'list' or from 'validation' mode
my $showAllAna=param('showAllAna') || 0;

####>Connecting to database<####
my $dbh=&promsConfig::dbConnect;

####>General info on project<####
my (@tempAnalysesID)=(param('id_ana'))? split(/,/,param('id_ana')) : $dbh->selectrow_array("SELECT ID_ANALYSIS FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=$proteinID LIMIT 0,1") ; # if no analysis was sent, take the first ID_ANALYSIS returned by this SELECT query
@tempAnalysesID=&promsMod::cleanNumericalParameters(@tempAnalysesID);
my $strgAnaList=join(",",@tempAnalysesID);
my $itemType=&promsMod::getItemType('ANALYSIS');
my (%nameItem,%validStatus,%dataFile,%fileFormat);
my $bestProtScore=0;
my $bestConfidence;
my @analysesID;
if (scalar @tempAnalysesID > 1) { # only validated analyses
	my $sthInfo=$dbh->prepare("SELECT NAME,VALID_STATUS,DATA_FILE,FILE_FORMAT,SCORE,CONF_LEVEL FROM ANALYSIS A,ANALYSIS_PROTEIN AP WHERE A.ID_ANALYSIS=? AND AP.ID_ANALYSIS=? AND ID_PROTEIN=$proteinID");
	foreach my $anaID (@tempAnalysesID) { # protein may not be in all ana listed...
		$sthInfo->execute($anaID,$anaID);
		my @info=$sthInfo->fetchrow_array;
		if ($info[0]) { # protein is indeed found in ana
			push @analysesID,$anaID;
			($nameItem{$anaID},$validStatus{$anaID},$dataFile{$anaID},$fileFormat{$anaID},my $protSc,my $protConf)=@info;
			$bestProtScore=$protSc if $protSc > $bestProtScore;
			$bestConfidence=$protConf if (!defined $bestConfidence || $protConf > $bestConfidence);
		}
	}
	$sthInfo->finish;
}
else { # val or not val analysis ($bestProtScore & $bestConfidence initialized later)
	@analysesID=($tempAnalysesID[0]);
	my $sthInfo = $dbh->prepare("SELECT NAME,VALID_STATUS,DATA_FILE,FILE_FORMAT FROM ANALYSIS WHERE ID_ANALYSIS=$analysesID[0]");
	$sthInfo->execute;
	($nameItem{$analysesID[0]},$validStatus{$analysesID[0]},$dataFile{$analysesID[0]},$fileFormat{$analysesID[0]})=$sthInfo->fetchrow_array;
	$sthInfo->finish;
}

my $projectID=&promsMod::getProjectID($dbh,$analysesID[0],'analysis');
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
my ($projectStatus)=$dbh->selectrow_array("SELECT STATUS FROM PROJECT WHERE ID_PROJECT=$projectID");
$projectStatus=0 unless $projectStatus;

####>Fetching protein status in Project<####
my ($currProtStatus,$validProtID,$notValidProtID,$masterProtID,%isoforms);
my ($protSeq,$protLength,$desProt,$identifier,$massProt,$alias,$comments,$nameOrganism,$date,$updateUser);

$currProtStatus = 0;

if ($idType eq 'valid') { # proteinID provided is from PROTEIN table
	$validProtID=$proteinID;
	$notValidProtID=0;
	$currProtStatus=1;
}
else { # proteinID provided is from PROTEIN_VALIDATION table
	$notValidProtID=$proteinID;
	($validProtID)=$dbh->selectrow_array("SELECT ID_PROTEIN FROM PROTEIN WHERE ID_PROJECT=$projectID AND IDENTIFIER=(SELECT IDENTIFIER FROM PROTEIN_VALIDATION WHERE ID_PROT_VALID=$notValidProtID)");
	#if ($validStatus==1 && $validProtID) { # partial validation
	#	($currProtStatus)=$dbh->selectrow_array("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=$analysisID AND ID_PROTEIN=$validProtID");
	#}
	#else { # never validated
	#	$currProtStatus=0;
	#}
	$validProtID=0 unless $validProtID;
}

my %proteinStatus;
my $sthProtCount = $dbh->prepare("SELECT 1 FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=? AND ID_PROTEIN=?");
foreach my $anaID (@analysesID){
	my $protStatus;
	if ($validStatus{$anaID}>=1 && $validProtID) { # partial validation
		$sthProtCount->execute($anaID, $validProtID);
		($protStatus) = $sthProtCount->fetchrow_array;
		$currProtStatus = 1;
		($bestProtScore,$bestConfidence)=$dbh->selectrow_array("SELECT SCORE,CONF_LEVEL FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=$validProtID AND ID_ANALYSIS=$anaID") unless defined $bestConfidence; # not yet set if single validated analysis
	}
	else {
		$protStatus = 0;
	}
	$proteinStatus{$anaID} = $protStatus;
}
$sthProtCount->finish;
####>Fetching protein data<####
if ($validProtID) { # Protein has been validated at least once
	my $protQuery="SELECT PROT_SEQ,PROT_LENGTH,PROT_DES,IDENTIFIER,MW,ALIAS,COMMENTS,ORGANISM,UPDATE_DATE,UPDATE_USER,ID_MASTER_PROTEIN FROM PROTEIN WHERE ID_PROTEIN=$validProtID"; #ID_CLUSTER,
	($protSeq,$protLength,$desProt,$identifier,$massProt,$alias,$comments,$nameOrganism,$date,$updateUser,$masterProtID)=$dbh->selectrow_array($protQuery); #$clusterID,
	$protSeq='' unless $protSeq;
	$protSeq=~s/\s+//g; # in case of \n
	if (length($protSeq) <= 1 && $masterProtID) { # should be '+' but in case unexpected problem ('-','' or NULL)
		($protSeq)=$dbh->selectrow_array("SELECT PROT_SEQ FROM MASTER_PROTEIN WHERE ID_MASTER_PROTEIN=$masterProtID");
	}
	$comments=&promsMod::HTMLcompatible($comments);

	###>Looking for isoforms (same UniProt entry)<###
	if ($masterProtID) {
		my $sthIso=$dbh->prepare("SELECT ID_PROTEIN,IDENTIFIER,ALIAS FROM PROTEIN WHERE ID_PROJECT=$projectID AND ID_MASTER_PROTEIN=$masterProtID");
		$sthIso->execute;
		while (my ($protID,$ident,$alias)=$sthIso->fetchrow_array) {
			@{$isoforms{$protID}}=($ident,$alias);
		}
		$sthIso->finish;
	}
}
else { # protein never validated
	my $protQuery="SELECT IDENTIFIER,PROT_LENGTH,PROT_DES,MW,ORGANISM FROM PROTEIN_VALIDATION WHERE ID_PROT_VALID=$notValidProtID";
	($identifier,$protLength,$desProt,$massProt,$nameOrganism)=$dbh->selectrow_array($protQuery);
	$alias=$identifier;
	$protSeq=&getProtSequence($identifier,$analysesID[0]);
}
$protLength=0 unless $protLength;

####>Var Mod info<####
my (%projectVarMods,%varModList,%varModPosOnProtein,%varModPosPeptides,%varModNames);
my $varModInProject=0;
my %allPostTransModifs = &promsMod::getVariableModifications($dbh);
#if ($ptmString) {
#	foreach my $varMod (split(';',$ptmString)) {$projectVarMods{$varMod}=$allPostTransModifs{$varMod};}
#}
#my $sthPepMod;
my $sthGetPM=$dbh->prepare("SELECT ID_MODIFICATION FROM PROJECT_MODIFICATION WHERE ID_PROJECT=$projectID");
$sthGetPM->execute;
while (my ($modID)=$sthGetPM->fetchrow_array) {
	$projectVarMods{$modID}=$allPostTransModifs{$modID};
}
$sthGetPM->finish;


####>Fetching matching peptide data<####

###>Peptide data for validated analyses
my $ghostFilterStrg=($bestConfidence && !$showMSdata)? 'AND VALID_STATUS>=1' : ''; # ghosts (ALWAYS display ghost peptides IF protein is virtual in all Analyses!)
my $MSdataStrg=($showMSdata)? ',MISS_CUT,MR_EXP,MR_CALC,MR_OBS,MR_DELTA,CHARGE' : '';
my $sthValPep=$dbh->prepare(qq
|SELECT GROUP_CONCAT(PEP_BEG SEPARATOR ':'),GROUP_CONCAT(PEP_END SEPARATOR ':'),P.ID_PEPTIDE,GROUP_CONCAT(FLANKING_AA SEPARATOR ':'),GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',COALESCE(PM.REF_POS_STRING,'') ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),DATA,PEP_SEQ,QUERY_NUM,PEP_RANK,0,VALID_STATUS,COMMENTS,SCORE$MSdataStrg,ELUTION_TIME,SUBST
	FROM PEPTIDE P
	LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
	JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
	WHERE P.ID_ANALYSIS=? AND ID_PROTEIN=$validProtID $ghostFilterStrg GROUP BY P.ID_PEPTIDE ORDER BY ABS(PEP_BEG),ABS(PEP_END),P.ID_PEPTIDE|
);
my (%peptideInfo,%pep_occurrence,%phosphoRsData,%maxQuantProb);
my (@pep_beg,@pep_end,@id_peptide,@flanking_AA,%pep_specificity);
my (%pepAna,%pepFile);

foreach my $analysisID (@analysesID){

	my $validStatus = $validStatus{$analysisID};
	my $dataFile = $dataFile{$analysisID};
	my $fileFormat = $fileFormat{$analysisID};

	if ($proteinStatus{$analysisID}) { # data for validated protein

		###>Peptide data for selected analysis
		$sthValPep->execute($analysisID);
		while (my ($begStrg,$endStrg,$pepID,$flAAStrg,$varModCode,$refModInfo,$pepData,@data)=$sthValPep->fetchrow_array) {
			my @begs=split(':',$begStrg); # in case multiple matches in protein sequence
			my @ends=split(':',$endStrg);
			my @flAAs=split(':',$flAAStrg);
			foreach my $i (0..$#begs) {
				push @pep_beg,abs($begs[$i]);
				push @pep_end,abs($ends[$i]);
				push @flanking_AA,$flAAs[$i];
				push @id_peptide,$pepID; # same idPep could be found more than once in array (multiple matches)
				$pep_occurrence{$pepID}{abs($begs[$i])}=abs($ends[$i]);
			}
			$pepAna{$pepID} = $analysisID;
			$data[1]=0 unless $data[1]; # queryNum (undef for MaxQuant SILAC peptides)
			$data[2]=0 unless $data[2]; # rank (undef for MaxQuant SILAC peptides)
			$peptideInfo{$pepID}=\@data;
			#$peptideInfo{$pepID}[3]=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$analysisID,$peptideInfo{$pepID}[0]); # pepSeq
			#$peptideInfo{$pepID}[3]=~s/([:.])-1/$1?/g; # if unknown position, -1 has to be removed
			$peptideInfo{$pepID}[3]=($varModCode)? ' + '.&promsMod::decodeVarMod($dbh,$data[0],$varModCode,\%varModNames) : ''; # $data[0]=pepSeq

			if ($varModCode) {
				# PhosphoRS #
				$phosphoRsData{$pepID} = $1 if ($pepData && $pepData =~ /PRS=([^##]+)/);

				# MaxQuant probalities #
				if ($fileFormat{$analysisID} eq 'MAXQUANT.DIR') {
					while ($varModCode=~/(\d+):([^&]+)/g) {
						my ($modID,$modPosStrg)=($1,$2);
						if ($refModInfo && $refModInfo=~/(^|&)$modID:[^&#]*##PRB_MQ=([^&#]+)/) {
							my $modProbStrg=$2;
							while ($modProbStrg =~ /([^,:]+):([^,]+)/g) {
								push @{$maxQuantProb{$pepID}{$modID}},[$1,$2];
							}
						}
						elsif ($data[6]) { # score
							foreach my $pos (split(/\./,$modPosStrg)) {
								push @{$maxQuantProb{$pepID}{$modID}},[$pos,1]; # set to 100%
							}
						}
					}
				}
			}
		}
		next if $#id_peptide <0; # @id_peptide is empty and make query $sthSpec fail !
		my $sthSpec=$dbh->prepare("SELECT ID_PEPTIDE,COUNT(DISTINCT ID_PROTEIN) FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_ANALYSIS=? AND ID_PEPTIDE IN (".(join(',',@id_peptide)).") GROUP BY ID_PEPTIDE");
		$sthSpec->execute($analysisID);
		while (my ($pID,$numMatchedProt)=$sthSpec->fetchrow_array) {
			$pep_specificity{$pID}=1*(sprintf '%.2f',100/$numMatchedProt);
		}
		$sthSpec->finish;

		#$sthPepMod=$dbh->prepare("SELECT ID_MODIFICATION,POS_STRING FROM PEPTIDE_MODIFICATION WHERE ID_PEPTIDE=?");
		if ($validStatus==2) {
			if (defined($dataFile)) {
				$dataFile='msms.txt' if $fileFormat eq 'MAXQUANT.DIR';
				$pepFile{$analysisID}="$promsPath{peptide}/proj_$projectID/ana_$analysisID/$dataFile"; # new data file structure (04/03/13)
				unless (-e $pepFile{$analysisID}) {
					if ($fileFormat eq 'MASCOT.DAT') {$pepFile{$analysisID}=~s/\.dat/_min\.dat/;} # "old" minimal file
					elsif ($fileFormat=~/\.XML/ && $fileFormat ne 'PARAGON.XML') {$pepFile{$analysisID}=~s/\.xml/\.pgf/;}
				}
			}
		}
		else {$pepFile{$analysisID}="$promsPath{valid}/ana_$analysisID/$dataFile";}
	}
	else { # data for non-validated protein !!! Only 1 ana selected !!!
		my $selStatus;
		($bestProtScore,$bestConfidence,$selStatus)=$dbh->selectrow_array("SELECT SCORE,CONF_LEVEL,SEL_STATUS FROM PROTEIN_VALIDATION WHERE ID_PROT_VALID=$notValidProtID");
		$bestConfidence=2 unless $bestConfidence;
		if ($selStatus>=1) {
			my $sthMatch=$dbh->prepare("SELECT QUERY_NUM,PEP_RANK,MATCH_INFO FROM RANK_PROTEIN_MATCH WHERE ID_ANALYSIS=$analysisID AND IDENTIFIER='$identifier'");
			#$sthPepMod=$dbh->prepare("SELECT ID_MODIFICATION,POS_STRING FROM QUERY_MODIFICATION WHERE ID_QUERY=?");
			my (%sthPep,%begList,%flkList);
			foreach my $r (1..10) {$sthPep{$r}=$dbh->prepare("SELECT ID_QUERY,INFO_PEP$r,MASS_DATA,CHARGE,ELUTION_TIME FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID AND QUERY_NUM=?");}
			my $sthSpec=$dbh->prepare("SELECT COUNT(DISTINCT RPM.IDENTIFIER) FROM RANK_PROTEIN_MATCH RPM,PROTEIN_VALIDATION PV WHERE RPM.IDENTIFIER=PV.IDENTIFIER AND RPM.ID_ANALYSIS=$analysisID AND PV.SEL_STATUS>=1 AND RPM.QUERY_NUM=? AND RPM.PEP_RANK=?");
			$sthMatch->execute;
			while (my ($qNum,$rank,$matchInfo)=$sthMatch->fetchrow_array) {
				$sthPep{$rank}->execute($qNum);
				my ($qID,$pepInfo,$massData,$charge,$elutionTime)=$sthPep{$rank}->fetchrow_array;
				my ($sel)=($pepInfo=~/SEL=(-?\d)/);
				next unless $sel>=1;
				my $rankID="$qID"."_$qNum"."_$rank";
				foreach my $match (split(':',$matchInfo)) { # beg,Nter,Cter:...:...
					#my ($beg)=($match=~/^(\d+)/);
					my ($beg,$flkNter,$flkCter)=split(',',$match);
					push @{$begList{$beg}},$rankID;
					$flkNter='' unless $flkNter;
					$flkCter='' unless $flkCter;
					$flkList{$beg}="$flkNter$flkCter";
				}
				$sthSpec->execute($qNum,$rank);
				my ($numMatchedProt)=$sthSpec->fetchrow_array;
				$pep_specificity{$rankID}=1*(sprintf '%.2f',100/$numMatchedProt); # %
				my ($pepSeq)=($pepInfo=~/SEQ=(\w+)/);
				#my ($varMod)=($pepInfo=~/VMOD=([^,]+)/); $varMod='' unless $varMod;
				my ($sub)=($pepInfo=~/SUBST=([^,]+)/); $sub='' unless $sub;
				my ($pepComments)=($pepInfo=~/COM=(.+)/);
				my ($pepSc)=($pepInfo=~/SC=(\-?\d+\.?\d*)/);
				#@{$peptideInfo{$rankID}}=($pepSeq,$qNum,$rank,$varMod,$sel,$pepComments,$pepSc);
				@{$peptideInfo{$rankID}}=($pepSeq,$qNum,$rank,$sel,$pepComments,$pepSc);
				if ($showMSdata) {
					my ($charge)=($pepInfo=~/MIS=(\d)/);
					my ($missCut)=($pepInfo=~/MIS=(\d)/);
					my ($mrExp)=($massData=~/EXP=(\d+\.*\d*)/);
					my ($mrCalc)=($pepInfo=~/CALC=(\d+\.*\d*)/);
					my ($mrObs)=($massData=~/OBS=(\d+\.*\d*)/);
					my ($mrDelta)=($pepInfo=~/DELT=(-*\d+\.*\d*)/);
					#my $deltaPPM=sprintf "%.3f",1000000*($mrDelta/$mrCalc);
					push @{$peptideInfo{$rankID}},($missCut,$mrExp,$mrCalc,$mrObs,$mrDelta,$charge,$elutionTime,$sub);
				}
				my ($PRS) = ($pepInfo =~ /PRS=([^,]+)/);
				$phosphoRsData{$rankID} = $PRS if $PRS;
			}
			foreach my $beg (sort{$a<=>$b} keys %begList) {
				foreach my $rankID (sort{length($peptideInfo{$a}[0])<=>length($peptideInfo{$b}[0])} @{$begList{$beg}}) {
					my ($qNum,$rank)=split(':',$rankID);
					push @pep_beg,$beg;
					push @pep_end,$beg+length($peptideInfo{$rankID}[0])-1;
					push @id_peptide,$rankID;
					$pepAna{$rankID} = $analysisID;
					push @flanking_AA,$flkList{$beg};
					$pep_occurrence{$rankID}{$beg}=1;
				}
			}
			$sthMatch->finish;
			foreach my $r (keys %sthPep) {$sthPep{$r}->finish;}
			$sthSpec->finish;
		}

		if ($selStatus<=-2) {
			my $word=($selStatus==-2)? 'excluded' : 'filtered';
			$comments="<FONT class=\"title3\" color='#DD0000'>This protein has been $word from the list of selectable proteins.</FONT>";
		}
		else {$comments="<FONT class=\"title3\" color='#DD0000'>This protein is still undergoing validation.</FONT>";}

		$pepFile{$analysisID}="$promsPath{valid}/ana_$analysisID/$dataFile";
	}

}

$sthValPep->finish;

####>Processing protein sequence<####
my $sthPepMod=($idType eq 'valid')? $dbh->prepare("SELECT ID_MODIFICATION,POS_STRING FROM PEPTIDE_MODIFICATION WHERE ID_PEPTIDE=?") : $dbh->prepare("SELECT ID_MODIFICATION,POS_STRING FROM QUERY_MODIFICATION WHERE ID_QUERY=?");
my (%boundaryStatus,@boundaryList,%unmatchedPeptides);
my $pepCoverage=0;
my $covProtLength=$protLength;
my $seqEndMatched=0;
my $seqMismatch=($protLength)? 0 : 1;
if (scalar @pep_beg) { # not if virtual protein
	$boundaryStatus{0}=0; # first residue is set to 0
	foreach my $i (0..$#pep_beg) {
		if ($pep_beg[$i]==0) { # MaxQuant with no protein sequence
			$unmatchedPeptides{$id_peptide[$i]}=1;
			next;
		}
		$boundaryStatus{$pep_beg[$i]-1}++; # -1 -> index not pos!
		$boundaryStatus{$pep_end[$i]-1}--; # -1 -> index not pos!
		$covProtLength=$pep_end[$i] if $pep_end[$i] > $covProtLength; # length might not be right due to protein isoforms
		$seqEndMatched=1 if (!$seqEndMatched && $flanking_AA[$i] && $flanking_AA[$i]=~/.-/); # peptide matches end of sequence
		unless ($seqMismatch) {
			$seqMismatch=1 if (!$seqMismatch && ($pep_end[$i]>$protLength || $peptideInfo{$id_peptide[$i]}[0] ne substr($protSeq,$pep_beg[$i]-1,$pep_end[$i]-$pep_beg[$i]+1)));
		}
	}
	##>Peptide match & coverage
	my $prevPos=0;
	my $matchBeg=0;
	foreach my $pos (sort {$a<=>$b} keys %boundaryStatus) {
		push @boundaryList,$pos;
		next if $pos==0;
		$boundaryStatus{$pos}+=$boundaryStatus{$prevPos}; # summing up all status
		if ($matchBeg==0 && $boundaryStatus{$pos}>=1) { # match begins
			$matchBeg=$pos;
		}
		elsif ($matchBeg>0 && $boundaryStatus{$pos}==0) { # match ends
			$pepCoverage+=($pos-$matchBeg+1);
			$matchBeg=0;
		}
		$prevPos=$pos;
	}
	$pepCoverage=($covProtLength)? sprintf "%.1f",($pepCoverage*100)/$covProtLength : 0; # if $protLength>0;
	$pepCoverage*=(!$seqEndMatched && $covProtLength > $protLength)? -1 : 1; # 25.0 -> 25   -1: flag for protLength problem

	if ($seqMismatch && $covProtLength) {
		$protSeq='X' x $covProtLength;
		$protSeq.='...' unless $seqEndMatched;
		foreach my $i (0..$#pep_beg) {
#print ">$peptideInfo{$id_peptide[$i]}[4] [$pep_beg[$i]-$pep_end[$i]]: ",$pep_beg[$i]-1,"-",$pep_end[$i]-$pep_beg[$i]+1," '$peptideInfo{$id_peptide[$i]}[0]' (",length($peptideInfo{$id_peptide[$i]}[0])," aa)<BR>\n";
			substr($protSeq,$pep_beg[$i]-1,$pep_end[$i]-$pep_beg[$i]+1)=$peptideInfo{$id_peptide[$i]}[0] if $pep_beg[$i]>0;
		}
	}

	###>PTMs
	#foreach my $i (0..$#id_peptide) {
	#	my ($pepSeq,$varModString)=@{$peptideInfo{$id_peptide[$i]}}[0,3];
	#	my @residueList = split(//,$pepSeq);
	#	$varModString=~s/^\s\+\s//; # string starts with ' + ';
	#	foreach my $varMod (split(/ \+ /,$varModString)) {
	#		my ($varModCode,$resid,$positions)=&promsMod::convertVarModString($varMod,1);
	#		if ($resid eq '-' || $resid eq '=' ) {
	#			$varModList{$varModCode}{$resid}{0}++;
	#			push @{$varModPosOnProtein{0}},$varModCode;
	#			push @{$varModPosPeptides{0}{$varModCode}},$id_peptide[$i];
	#		}
	#		elsif ($resid eq '+' || $resid eq '*' ) {
	#			$varModList{$varModCode}{$resid}{999999}++;
	#			push @{$varModPosOnProtein{999999}},$varModCode;
	#			push @{$varModPosPeptides{999999}{$varModCode}},$id_peptide[$i];
	#		}
	#		else {
	#			foreach my $pos (split(/\./,$positions)) {
	#				#--- Part modified on  07/06/12
	#				#if ($pos==-1) { # unknown position => find 1st matching residue on seq
	#				#	$resid=(join('|',split(//,$resid)));
	#				#	$pos=pos($pepSeq) if $pepSeq=~/$resid/g;
	#				#}
	#				#my $residue = $residueList[$pos-1];
	#				#my $posOnProt=$pep_beg[$i]+$pos-1;
	#				#$varModList{$varModCode}{$residue}{$posOnProt}++;
	#				#$varModPosOnProtein{$posOnProt}=$varModCode;
	#				#push @{$varModPosPeptides{$posOnProt}},$id_peptide[$i];
	#
	#				if ($pos==-1) { # unknown position => find all matching residues on seq
	#					my $tagPos=0;
	#					$resid=(join('|',split(//,$resid)));
	#					foreach my $pepResidu (@residueList) {
	#						$tagPos++;
	#						if ($resid =~ /$pepResidu/ ) {
	#							my $posOnProt=$pep_beg[$i]+$tagPos-1;
	#							$varModList{$varModCode}{$pepResidu}{$posOnProt}++;
	#							$varModPosOnProtein{$posOnProt}{$varModCode}=-1;# tag to know if the modification is ambiguous
	#							push @{$varModPosPeptides{$posOnProt}{$varModCode}},$id_peptide[$i];
	#						}
	#					}
	#				}
	#				else {
	#					my $residue = $residueList[$pos-1];
	#					my $posOnProt=$pep_beg[$i]+$pos-1;
	#					$varModList{$varModCode}{$residue}{$posOnProt}++;
	#					$varModPosOnProtein{$posOnProt}{$varModCode}=1;
	#					push @{$varModPosPeptides{$posOnProt}{$varModCode}},$id_peptide[$i];
	#				}
	#			}
	#		}
	#	}
	#}

	##>PTMs
	foreach my $i (0..$#id_peptide) {
next if $pep_beg[$i]==0; # skip if beg is not known
		my $pepSeq=$peptideInfo{$id_peptide[$i]}[0];
		my @residueList = split(//,$pepSeq);
		$sthPepMod->execute($id_peptide[$i]);
		while (my ($modID,$posString) = $sthPepMod->fetchrow_array) {
			if (!defined($varModNames{$modID})) { # For PTM not retrieved from promsMod::getVariableModifications (i.e. not valid dispcode and dispcolor)
				my ($psiMsName,$interimName,$altNames)=$dbh->selectrow_array("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES FROM MODIFICATION WHERE ID_MODIFICATION=$modID");
				my $name=($psiMsName)?$psiMsName:($interimName)?$interimName:($altNames)?$altNames:'';
				$varModNames{$modID}=$name;
			}
			foreach my $pos (split(/\./,$posString)) {
				###> Update $pos for N/C-term to avoid an error for if($pos==1) -> Argument "=" isn't numeric in numeric eq (==)
				# Moved on 18/11/13
				#$pos=0 if ($pos eq '-' || $pos eq '=' );
				#$pos=999999 if ($pos eq '+' || $pos eq '*' );
				#if ($pos==-1) { # unknown position => find all matching residues on seq
				if ($pos =~ /(-(\d+|=|-|\+|\*);{0,1})+/) {
					#my $tagPos=0;
					#my ($resid)=$dbh->selectrow_array("SELECT SPECIFICITY FROM MODIFICATION WHERE ID_MODIFICATION=$modID");
					#$resid=(join('|',split(/,/,$resid)));
					#foreach my $pepResidu (@residueList) {
					#	$tagPos++;
					#	if ($resid =~ /$pepResidu/ ) {
					#		my $posOnProt=$pep_beg[$i]+$tagPos-1;
					#		$varModList{$modID}{$pepResidu}{$posOnProt}++;
					#		$varModPosOnProtein{$posOnProt}{$modID}=-1;# tag to know if the modification is ambiguous
					#		push @{$varModPosPeptides{$posOnProt}{$modID}},$id_peptide[$i];
					#	}
					#}
					foreach my $negPos (split (/;/,$pos)) {
						$negPos=0 if ($negPos eq '--' || $negPos eq '-=' );
						$negPos=999999 if ($negPos eq '-+' || $negPos eq '-*' );
						my $residue = $residueList[abs($negPos)-1];
						my $posOnProt=$pep_beg[$i]+abs($negPos)-1;
						$varModList{$modID}{$residue}{$posOnProt}-- unless $varModList{$modID}{$residue}{$posOnProt};
						$varModPosOnProtein{$posOnProt}{$modID}=-1;# tag to know if the modification is ambiguous
						push @{$varModPosPeptides{$posOnProt}{$modID}},$id_peptide[$i];
					}
				}
				else {
					#my $residue = $residueList[$pos-1];
					#my $posOnProt=$pep_beg[$i]+$pos-1;
					# Change on 18/11/13 because several problems were reported in 'Detailed sequence coverage' for Protein N-term and N-term
					$pos=999999 if ($pos eq '+' || $pos eq '*' );
					my $residue= ($pos eq '-' || $pos eq '=') ? $residueList[0] : $residueList[$pos-1];
					my $posOnProt=($pos eq '-') ? 0 : ($pos eq '=') ? $pep_beg[$i] : $pep_beg[$i]+$pos-1;
					$varModList{$modID}{$residue}{$posOnProt}++;
					$varModPosOnProtein{$posOnProt}{$modID}=1;
					push @{$varModPosPeptides{$posOnProt}{$modID}},$id_peptide[$i];
				}
			}
		}
	}
}
$sthPepMod->finish;

##BEGIN
####>Finding info on Experiments, Samples and Analyses where protein is found<####
my (%analysisInfo,%analysisParents,%protVarMods);
my $parentsQuery1=qq |SELECT 2,EXPERIMENT.NAME,EXPERIMENT.DISPLAY_POS,SAMPLE.NAME,SAMPLE.DISPLAY_POS,'A' FROM EXPERIMENT,SAMPLE,ANALYSIS
						WHERE ID_ANALYSIS=? AND SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND EXPERIMENT.ID_EXPERIMENT=SAMPLE.ID_EXPERIMENT AND SAMPLE.ID_SPOT IS NULL|;
my $parentsQuery2=qq |SELECT 1,EXPERIMENT.NAME,EXPERIMENT.DISPLAY_POS,GEL2D.NAME,GEL2D.DISPLAY_POS,SPOT.NAME FROM EXPERIMENT,GEL2D,SPOT,SAMPLE,ANALYSIS
						WHERE ID_ANALYSIS=? AND SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND SPOT.ID_SPOT=SAMPLE.ID_SPOT AND GEL2D.ID_GEL2D=SPOT.ID_GEL2D AND EXPERIMENT.ID_EXPERIMENT=GEL2D.ID_EXPERIMENT|;
my @sthAnaParents=($dbh->prepare($parentsQuery1),$dbh->prepare($parentsQuery2)); # free sample then gel

###>Analyses where protein is validated<###
#my $numConnProt=0;
if ($validProtID) {
	my $sthVal;
	if ($showAllAna) {
		$sthVal=$dbh->prepare("SELECT A.ID_ANALYSIS,NAME,DISPLAY_POS,MS_TYPE,VALID_STATUS,CONF_LEVEL,VISIBILITY,NUM_PEP,MATCH_GROUP,PEP_SPECIFICITY,PEP_COVERAGE,ID_SAMPLE FROM ANALYSIS A,ANALYSIS_PROTEIN AP WHERE A.ID_ANALYSIS=AP.ID_ANALYSIS AND AP.ID_PROTEIN=$validProtID");
	}
	else {
		my $strgID = join(',',@analysesID);
		$sthVal=$dbh->prepare("SELECT A.ID_ANALYSIS,NAME,DISPLAY_POS,MS_TYPE,VALID_STATUS,CONF_LEVEL,VISIBILITY,NUM_PEP,MATCH_GROUP,PEP_SPECIFICITY,PEP_COVERAGE,ID_SAMPLE FROM ANALYSIS A,ANALYSIS_PROTEIN AP WHERE A.ID_ANALYSIS=AP.ID_ANALYSIS AND A.ID_ANALYSIS IN ($strgID) AND AP.ID_PROTEIN=$validProtID");
	}
	my $sthGr1=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE MATCH_GROUP=? AND ID_ANALYSIS=?");
	my $sthGr2=$dbh->prepare("SELECT MAX(NUM_PEP) FROM ANALYSIS_PROTEIN WHERE MATCH_GROUP=? AND ID_ANALYSIS=? AND ID_PROTEIN != $validProtID");
	my $sthGr3=$dbh->prepare("SELECT ALIAS,NUM_PEP FROM ANALYSIS_PROTEIN,PROTEIN WHERE ANALYSIS_PROTEIN.ID_PROTEIN=PROTEIN.ID_PROTEIN AND VISIBILITY=2 AND MATCH_GROUP=? AND ID_ANALYSIS=?");
	$sthVal->execute;
	my $sthPepVMod=$dbh->prepare("SELECT DISTINCT ID_MODIFICATION FROM PEPTIDE_PROTEIN_ATTRIB PPA,PEPTIDE_MODIFICATION PM WHERE PPA.ID_PEPTIDE=PM.ID_PEPTIDE AND PPA.ID_ANALYSIS=? AND ID_PROTEIN=$validProtID"); # AND VALID_STATUS > 0 include ghost
	my $sthGhostPep=$dbh->prepare("SELECT COUNT(*) FROM (SELECT PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&')
									FROM PEPTIDE P
									LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
									JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
									WHERE P.ID_ANALYSIS=? AND ID_PROTEIN=$validProtID AND VALID_STATUS=0 GROUP BY P.ID_PEPTIDE) AS PEP"); # should use DISTINCT or Not?

	while (my ($anaID,$anaName,$anaPos,$msType,$anaValid,$conf,$vis,$numPep,$matchGr,$pepSpec,$pepCov,$sampID)=$sthVal->fetchrow_array) {
		$pepCov='?' unless $pepCov;
		foreach my $sthPar (@sthAnaParents) { # 2 parent paths are possible
			$sthPar->execute($anaID);
			@{$analysisParents{$anaID}}=$sthPar->fetchrow_array;
			last if $analysisParents{$anaID}[0];
		}
		my ($groupBest,$groupClass);
		$sthGr1->execute($matchGr,$anaID);
		my ($groupSize)=$sthGr1->fetchrow_array;
		if ($groupSize>1) {
			if ($vis==2) {
				$groupBest=$alias;
				$sthGr2->execute($matchGr,$anaID);
				my ($maxPep)=$sthGr2->fetchrow_array;
				$groupClass=($numPep > $maxPep)? 'matchGroupOK' : ($numPep==$maxPep)? 'matchGroupEQ' : 'matchGroupBAD';
			}
			else { # protein is not best of its group
				$sthGr3->execute($matchGr,$anaID);
				($groupBest,my $numPepBest)=$sthGr3->fetchrow_array;
				if ($vis==1) {
					$groupClass=($numPep < $numPepBest)? 'matchGroupNO' : ($numPep==$numPepBest)? 'matchGroupEQ' : ($numPep==$numPepBest)? 'matchGroupEQ' : 'matchGroupBAD';
				}
				else { # hidden
					$groupClass=($numPep < $numPepBest)? 'matchGroupOK' : ($numPep==$numPepBest)? 'matchGroupEQ' : 'matchGroupBAD';
				}
			}
		}
		else {
			$groupBest='none';
			$groupClass='matchGroupOK';
		}
		if ($conf) { # Not a ghost protein => look for ghost peptides
			$sthGhostPep->execute($anaID);
			my ($numGhostPep)=$sthGhostPep->fetchrow_array;
			$numPep.="+<FONT class=\"virtualProt\" onmouseover=\"popup('Peptides retrieved from XIC extraction')\" onmouseout=\"popout()\">$numGhostPep</FONT>" if $numGhostPep;
		}
		@{$analysisInfo{$anaID}}=($anaName,$anaPos,$conf,$vis,$matchGr,$groupBest,$msType,$anaValid,1,$numPep,$pepSpec,$pepCov,$groupClass); # [8]=1: prot is validated in ana

		###>> Get PTM of all peptides for this protein in this $anaID
		$sthPepVMod->execute($anaID);
		while (my ($modID)=$sthPepVMod->fetchrow_array) {
			$protVarMods{$anaID}{$modID}=1 if $projectVarMods{$modID};
		}
	}
	$sthVal->finish;
	$sthGr1->finish;
	$sthGr2->finish;
	$sthGr3->finish;
	$sthPepVMod->finish;
	$sthGhostPep->finish;

	###>All connected proteins in project
	#my %connectedProteins;
	#&getNumberConnectedProteins($dbh,$validProtID,$projectID,\%connectedProteins);
	#$numConnProt=(scalar keys %connectedProteins)-1;
}

###>Analyses where protein is not validated<###
if ($showAllAna) {
	my $noValQuery=qq|SELECT A.ID_ANALYSIS,A.NAME,A.DISPLAY_POS,A.MS_TYPE,A.VALID_STATUS,PV.ID_PROT_VALID,PV.NUM_MATCH,PV.MAX_MATCH,A.ID_SAMPLE
						FROM PROTEIN_VALIDATION PV,ANALYSIS A,SAMPLE S,EXPERIMENT E
						WHERE PV.IDENTIFIER='$identifier'
							AND PV.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_EXPERIMENT=E.ID_EXPERIMENT
							AND E.ID_PROJECT=$projectID|;
	my $sthNoValP=$dbh->prepare($noValQuery);
	$sthNoValP->execute;
	while (my ($anaID,$anaName,$anaPos,$msType,$anaValid,$protID,$numPep,$maxPep,$sampID)=$sthNoValP->fetchrow_array) {
		next if $analysisInfo{$anaID}; # ana is partially validated and prot is validated
		foreach my $sthPar (@sthAnaParents) { # 2 parent pathes are possible
			$sthPar->execute($anaID);
			@{$analysisParents{$anaID}}=$sthPar->fetchrow_array;
			last if $analysisParents{$anaID}[0];
		}
		@{$analysisInfo{$anaID}}=($anaName,$anaPos,2,1,0,'',$msType,$anaValid,0,"$numPep/$maxPep",'-','-',$protID);  # [8]=0: prot is not validated in ana
	}
	#foreach my $sthPar (@sthAnaParents) {$sthPar->finish;}
	$sthNoValP->finish;
}
foreach my $sthPar (@sthAnaParents) {
	$sthPar->finish;
}


####>Quantifications<#### finding list of designs in project
my $sthDsg=$dbh->prepare("SELECT D.ID_DESIGN,E.DISPLAY_POS,E.NAME,D.NAME FROM EXPERIMENT E,DESIGN D WHERE E.ID_EXPERIMENT=D.ID_EXPERIMENT AND E.ID_PROJECT=$projectID");
my %designs;
$sthDsg->execute;
while (my ($designID,@desInfo)=$sthDsg->fetchrow_array) {
	@{$designs{$designID}}=@desInfo;
}
$sthDsg->finish;

####>Gene info<####
my %geneInfo;
#if ($validProtID) {
#	my $sthCR=$dbh->prepare("SELECT GENE.ID_GENE,GENE_NAME,GENE_SYMBOLS FROM GENE,GENE_PROTEIN WHERE GENE.ID_GENE=GENE_PROTEIN.ID_GENE AND ID_PROTEIN=$validProtID");
#	$sthCR->execute;
#	while (my ($geneID,$geneName,$geneSymbStrg)=$sthCR->fetchrow_array) {
#		$geneSymbStrg='' unless $geneSymbStrg;
#		@{$geneInfo{$geneID}}=($geneSymbStrg,$geneName);
#	}
#	$sthCR->finish;
#}
if ($masterProtID) {
	my $sthGN=$dbh->prepare("SELECT MI.VALUE FROM MASTERPROT_IDENTIFIER MI,IDENTIFIER I WHERE ID_MASTER_PROTEIN=$masterProtID AND MI.ID_IDENTIFIER=I.ID_IDENTIFIER AND I.CODE='GN' ORDER BY RANK");
	$sthGN->execute;
	while (my ($geneName)=$sthGN->fetchrow_array) {
		push @{$geneInfo{'GN'}},$geneName;
	}
	$sthGN->finish;

	@{$geneInfo{'OS'}}=$dbh->selectrow_array("SELECT SCIENTIFIC_NAME,COMMON_NAME FROM MASTER_PROTEIN MP,SPECIES S WHERE ID_MASTER_PROTEIN=$masterProtID AND MP.ID_SPECIES=S.ID_SPECIES");
}

####>Alias info<####
my %isoformProteins;
if ($masterProtID) {
	#my $sthAlias=$dbh->prepare("SELECT P.ID_PROTEIN,ID_ANALYSIS,ALIAS,MAX(VISIBILITY),MW,PROT_DES,ORGANISM FROM PROTEIN P,ANALYSIS_PROTEIN AP WHERE P.ID_PROTEIN=AP.ID_PROTEIN AND ID_MASTER_PROTEIN=$masterProtID AND ID_PROJECT=$projectID ORDER BY GROUP BY P.ID_PROTEIN");
	my $sthAlias=$dbh->prepare("SELECT ID_PROTEIN,ID_ANALYSIS,ALIAS,MAX(CONF_LEVEL),MAX(VISIBILITY),MW,PROT_DES,ORGANISM FROM (SELECT P.ID_PROTEIN,ID_ANALYSIS,ALIAS,CONF_LEVEL,VISIBILITY,MW,PROT_DES,ORGANISM FROM PROTEIN P,ANALYSIS_PROTEIN AP WHERE P.ID_PROTEIN=AP.ID_PROTEIN AND ID_MASTER_PROTEIN=$masterProtID AND ID_PROJECT=$projectID ORDER BY VISIBILITY DESC) AS TMP_RESULTS GROUP BY ID_PROTEIN");
	$sthAlias->execute;
	while (my ($protID,@protData)=$sthAlias->fetchrow_array) {
		#next if $protID==$validProtID;
		@{$isoformProteins{$protID}}=@protData;
	}
	$sthAlias->finish;
}
####>Cross references<####
my (%crossReference,%crossRefInfo,@uniprotItem);
if ($masterProtID) {
	my $sthXref=$dbh->prepare("SELECT I.ID_IDENTIFIER,MI.RANK,MI.VALUE,I.CODE,I.NAME,I.RESRC_URL,I.RESRC_NAME FROM MASTERPROT_IDENTIFIER MI,IDENTIFIER I WHERE ID_MASTER_PROTEIN=$masterProtID AND MI.ID_IDENTIFIER=I.ID_IDENTIFIER ORDER BY I.ID_IDENTIFIER,RANK");
	$sthXref->execute;
	while (my ($identID,$rank,$identValue,$identCode,$identName,$identURL,$resrcName)=$sthXref->fetchrow_array) {
		next unless $identURL;
		#push @{$crossReference{$identID}},$uniprotID if ($identCode eq 'AC' && $rank==1);
		push @{$crossReference{$identID}},$identValue;
		push @uniprotItem,$identValue if $identCode eq 'AC';
		@{$crossRefInfo{$identID}}=($identCode,$identName,$identURL,$resrcName);
	}
	$sthXref->finish;
}
#>Isoform AC for ProtVista
my ($isoformAC,$disabProtVista)=($identifier=~/(\||^)([A-Z]\d\w{3}\d\.\d+)(\||$)/)? ($2,'') : ($uniprotItem[0])? ($uniprotItem[0],'') : ('-','disabled');
#>GI info
if ($identifier=~/(gi\|\d+)/) { # adds gi link
	my $identValue=$1;
	my ($identID,@giData)=$dbh->selectrow_array("SELECT ID_IDENTIFIER,CODE,NAME,RESRC_URL,RESRC_NAME FROM IDENTIFIER WHERE CODE='GI'");
	push @{$crossReference{$identID}},$identValue;
	@{$crossRefInfo{$identID}}=@giData;
}

$dbh->disconnect;

###########################
#      Starting HTML      #
###########################
my ($protVistaScript,$protVistaLink)=($disabProtVista)? ('','') : (
	"<SCRIPT src=\"http://ebi-uniprot.github.io/CDN/protvista/protvista.js\"></SCRIPT>\n",
	"<LINK href=\"http://ebi-uniprot.github.io/CDN/protvista/css/main.css\" rel=\"stylesheet\"/>\n"
);
print header(-'content-encoding'=>'no',-charset=>'UTF-8');
warningsToBrowser(1); # writes PERL warnings as HTML comments <!-- --> must be activated after header
print qq
|<HEAD>
<TITLE>Protein Entry: $alias</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
$protVistaLink
<STYLE type="text/css">
.selectedAna {color:#DD0000;}
.hidden {visibility:hidden;}
.matchGroupOK {background-color:#50DD50;}
.matchGroupEQ {background-color:#DDDD00;}
.matchGroupNO {background-color:#FF8A30;}
.matchGroupBAD {background-color:#FF6060;}
.overlapRes0 {color:black;}
.overlapRes1 {font-weight:bold;color:red;}
.overlapRes2 {font-weight:bold;color:blue;}
|;
foreach my $varMod (sort{$a cmp $b} keys %projectVarMods) {
	print ".$projectVarMods{$varMod}[2] \{font-weight:bold; $projectVarMods{$varMod}[3]\}\n";
	$varModInProject++;
}
print qq
|.modRes {text-decoration:overline} /*text-transform:lowercase;font-style:italic;*/
.popup {background-color:#FFFFFF;border:solid 3px #999999;padding:5px;box-shadow:10px 10px 20px #808080;position:absolute;display:none;} //z-index:999; <- conflits with mouseover popup
</STYLE>
$protVistaScript
<SCRIPT LANGUAGE="JavaScript">
|;
&promsMod::popupInfo();
my $varModCall=($currProtStatus)? 'sequenceView' : 'listRanks';
my $anaIDstring = join ',', @analysesID;
print qq
|function updateSectionVisibility() {
	if (!opener.top.protSectionVisibility) {opener.top.protSectionVisibility={sequence:true}} // safety
	var visSections=opener.top.protSectionVisibility;
	//console.log(visSections);
	for (var sectionName in visSections) {
		//console.log(sectionName);
		if (visSections[sectionName]==true && !document.getElementById(sectionName+'BUTTON').disabled) {
			var optParam=(sectionName=='protVista')? '$isoformAC' : null;
			showHideSection(sectionName,optParam);
		}
	}
}
var protVistaDisplayed=0;
function showHideSection(sectionName,optParam) {
	var vis,butValue;
	var visSections=opener.top.protSectionVisibility;
	var sectionDiv=document.getElementById(sectionName+'DIV');
	if (sectionDiv.style.display=='') { // hide
		vis='none';
		butValue='Show';
		visSections[sectionName]=false;
	}
	else { // show
		vis='';
		butValue='Hide';
		visSections[sectionName]=true;
	}
	sectionDiv.style.display=vis;
	document.getElementById(sectionName+'BUTTON').value=butValue;
	if (sectionName=='analysis') {
		document.getElementById('anaButtons').style.display=vis;
		//document.getElementById('pepDistrib').style.display=vis;
	}
	else if (sectionName=='protVista' && protVistaDisplayed==0) {
		protVistaDisplayed=1;
		var pvDiv = document.getElementById('protVistaTrueDIV');
        var ProtVista = require('ProtVista');
        var instance = new ProtVista({
            el: pvDiv,
            uniprotacc: optParam
        });
	}
}
function extractSequence(seq) {
	var covTable=document.getElementById('covTABLE');
	document.getElementById('extractionTEXT').value=seq;
	var extrDiv=document.getElementById('extractionDIV');
	var divX=getDomOffset(covTable,'offsetLeft') + 20;
    var divY=getDomOffset(covTable,'offsetTop');
	extrDiv.style.left = divX + 'px';
	extrDiv.style.top = divY + 'px';
	extrDiv.style.display='block';
}
function inNavigationMode(action) {
	if (!opener \|\| !opener.promsFrame \|\| !opener.promsFrame.resultFrame) {
		alert(action+' is available only in Project Navigation mode.');
		return false;
	}
	return true;
}
function sequenceView(analysisId,proteinId,idType,ms) {
    if (idType == undefined) {
		window.location="$promsPath{cgi}/sequenceView.cgi?id_ana="+analysisId+"&id_prot="+proteinId+"&msdata="+opener.top.showMSdata+"&showAllAna=$showAllAna";
	}
	else{
		if (ms==-1 \|\| ms==2) {opener.top.protSectionVisibility.peptide=true;} // make sure peptide section is visible after reload
		opener.top.showMSdata=(ms > 0)? 1 : 0; // 2 or -1 if just changed 1/0 otherwised
		window.location="$promsPath{cgi}/sequenceView.cgi?id_ana="+analysisId+"&id_prot="+proteinId+"&id_type="+idType+"&msdata="+opener.top.showMSdata+"&showAllAna=$showAllAna";
	}

}
function editProt() {
	window.location="$promsPath{cgi}/editProtein.cgi?id_ana=$anaIDstring&id_prot=$validProtID&id_project=$projectID";
}
function editMatchGroup(analysisId,group) {
	if (!inNavigationMode('Match Group Edition')) {return;}
	var selBranch=opener.promsFrame.selectedBranchID.split(":");
	opener.top.promsFrame.resultFrame.location="$promsPath{cgi}/editMatchGroup.cgi?id_item="+selBranch[1]+"&item="+selBranch[0]+"&id_analysis="+analysisId+"&group="+group;
}
function graphicalView(field){
	if (testCheckbox(field,'Analysis')==0) {return;}
	document.anaList.largWindow.value=(document.body)? document.body.clientWidth : window.innerWidth;
	document.anaList.graphView.value=opener.graphView;
	document.anaList.action="./graphicalView.cgi";
	document.anaList.submit();
}
function updateAnalysesList(showAll) {
	document.anaList.showAllAna.value=showAll;
	document.anaList.submit();
}
function graphicalView2(field) {
	if (testCheckbox(field,'Protein')==0) {return;}
	document.protList.largWindow.value=(document.body)? document.body.clientWidth : window.innerWidth;
	document.protList.graphView.value=opener.graphView;
	document.protList.submit();
}
function drawSpectrum(pepId,pepInfo) {
	if (pepId != selectedPep \|\| spectWin.closed) {
		selectPeptide(pepId);
		var call=($currProtStatus==1)? 'pep' : 'seq';
		var paramString="RID="+pepInfo+"&CALL="+call;
		spectWin=window.open("$promsPath{cgi}/drawSpectrum.cgi?"+paramString,'SpectrumWindow','width=950,height=950,location=no,resizable=yes,scrollbars=yes');
	}
	spectWin.focus();
}
function closeSpectWin() {
	if (spectWin != null && !spectWin.closed) {
		spectWin.close();
	}
}
function selectPeptide(newPep) {
	if (selectedPep && document.getElementById(selectedPep)) { // selectedPep could be in a closed varMod popup
		document.getElementById(selectedPep).style.color='#000000';
	}
	document.getElementById(newPep).style.color='#DD0000'; //style.background = '#DD0000';
	selectedPep=newPep;
}
function editVarModification(vmodId,seqLength) {
	//highlight var mod
	if (selectedVarModId) {
		document.getElementById(selectedVarModId).style.color='#000000';
	}
	document.getElementById(vmodId).style.color='#DD0000'; //style.background = '#DD0000';
	selectedVarModId=vmodId;
	//popup window
	var winLocation="$promsPath{cgi}/editVarModification.cgi?CALL=$varModCall&ID="+vmodId;
	var winWidth=500+seqLength*25;
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
function startValidation(anaID) {
	opener.promsFrame.location="$promsPath{cgi}/startValidation.cgi?ID="+anaID+"&SEL=$identifier";
}
function testCheckbox(field,label) {
	var checked=0;
	if (field.length) {
		for (i=0; i < field.length; i++) {
			if (field[i].checked==true) {checked++;}
		}
	}
	else {
		if (field.checked==true) {checked++;}
	}
	if (checked==0) {
		alert('You must select at least 1 '+label);
	}
	return checked;
}
function check(field) {
	if (checkflag==false) {
		if (field.length) {
			for (i=0; i < field.length; i++) {
				field[i].checked=true;
			}
		}
		else {field.checked=true;}
		checkflag=true;
	}
	else {
		if (field.length) {
			for (i=0; i < field.length; i++) {
				field[i].checked=false;
			}
		}
		else {field.checked=false;}
		checkflag=false;
	}
}

var checkflag=false;
var spectWin;
var selectedPep;
var selectedIdentifier='$identifier'; // used remotely by listRanks if changes in protein status
var varModWindow; // window used to edit variable modifications
var selectedVarModId;

// AJAX --->
var XHR=null;
var selectedImg;
function ajaxGetProteinsMatchGroup(anaID) {
	if (selectedImg) {selectedImg.src='$promsPath{images}/plus.gif';}
	var grDIV=document.getElementById('groupDIV');
	var linkImg=document.getElementById('img_$proteinID:'+anaID);
	if (grDIV.style.display == 'none' \|\| selectedImg.id != linkImg.id) {
		linkImg.src='$promsPath{images}/minus1.gif';
		grDIV.style.display='block';
	}
	else {
		linkImg.src='$promsPath{images}/plus.gif';
		grDIV.style.display='none';
		return;
	}
	selectedImg=linkImg;
	//wait image
	grDIV.innerHTML='<BR><TABLE><TR><TD>&nbsp;&nbsp;&nbsp;</TD><TD><IMG src="$promsPath{images}/scrollbarGreen.gif"></TD><TD>&nbsp;&nbsp;&nbsp;</TD></TR></TABLE><BR>';
	var divX=getDomOffset(linkImg,'offsetLeft') - 700;
    if (divX < 0) divX=0;
    var divY=getDomOffset(linkImg,'offsetTop') + linkImg.offsetHeight;
	grDIV.style.left = divX + 'px';
	grDIV.style.top = divY + 'px';

	//If XHR object already exists, the request is canceled & the object is deleted
	if(XHR && XHR.readyState != 0){
		XHR.abort();
		delete XHR;
	}
	//Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/listConflicts.cgi?AJAX=MGr&anaID="+anaID+"&protID=$proteinID",true);
    XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			grDIV.innerHTML=(XHR.responseText);
		}
	}
	XHR.send(null);
}
function ajaxGetVarModPeptideInfo(varModCode,pos,peptideStrg,cumulate) {
	var vModA=document.getElementById('vmod_'+pos); // <A> DOM object
	var pepDIV=document.getElementById('vmodPeptideDIV');
	pepDIV.style.display='block';

	if(!cumulate){
		cumulate = 0;
	}

	//wait image
	pepDIV.innerHTML='<BR><TABLE><TR><TD>&nbsp;&nbsp;&nbsp;</TD><TD><IMG src="$promsPath{images}/scrollbarGreen.gif"></TD><TD>&nbsp;&nbsp;&nbsp;</TD></TR></TABLE><BR>';
	var divX=getDomOffset(vModA,'offsetLeft'); // - Math.round(pepDIV.offsetWidth/2) + 5
    if (divX < 0) divX=0;
    var divY=getDomOffset(vModA,'offsetTop') + vModA.offsetHeight + 15;
	pepDIV.style.left = divX + 'px';
	pepDIV.style.top = divY + 'px';

	//If XHR object already exists, the request is canceled & the object is deleted
	if(XHR && XHR.readyState != 0){
		XHR.abort();
		delete XHR;
	}
	//Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/sequenceView.cgi?AJAX=vmodPep&projectID=$projectID&id_prot=$proteinID&id_type=$idType&vModCode="+varModCode+"&position="+pos+"&pepStrg="+peptideStrg+"&cumulate="+cumulate,true);
    XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			pepDIV.innerHTML=(XHR.responseText);
		}
	}
	XHR.send(null);
}
function ajaxGetConnectedProteins(proteinID) {
	//check if query already performed
	if (queryWasPerformed) {
		document.getElementById('connProtDiv').style.display='';
		document.getElementById('showProt').style.display='none';
		document.getElementById('hideProt').style.display='';
		return;
	}

	//wait image
	protDiv=document.getElementById('connProtDiv');
	protDiv.innerHTML='&nbsp;&nbsp;&nbsp;<IMG src="$promsPath{images}/scrollbarGreen.gif">';
	protDiv.style.display='block';

	//If XHR object already exists, the request is canceled & the object is deleted
	if(XHR && XHR.readyState != 0){
		XHR.abort();
		delete XHR;
	}
	//Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/sequenceView.cgi?AJAX=connProt&projectID=$projectID&id_prot="+proteinID,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			listConnectedProteins(XHR.responseText);
		}
	}
	XHR.send(null);
}
function ajaxGetInteractors(uniAC) {
	//check if query already performed
	if (interactSearchWasPerformed) {
		document.getElementById('interactDiv').style.display='';
		document.getElementById('showInteract').style.display='none';
		document.getElementById('hideInteract').style.display='';
		return;
	}
	//wait image
	var interDiv=document.getElementById('interactDiv');
	interDiv.innerHTML='&nbsp;&nbsp;&nbsp;<IMG src="$promsPath{images}/scrollbarGreen.gif">';
	interDiv.style.display='';
	//If XHR object already exists, the request is canceled & the object is deleted
	if(XHR && XHR.readyState != 0){
		XHR.abort();
		delete XHR;
	}
	//Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/searchInteractors.cgi?ACT=ajaxSearchInteractors&uniAC="+uniAC+"&ID=$strgAnaList&ITEM=ANALYSIS&id_project=$projectID",true);
		XHR.onreadystatechange=function() {
			if (XHR.readyState==4 && XHR.responseText) {
				listInteractors(XHR.responseText);
			}
	}
	XHR.send(null);
}

function ajaxGetQuantificationsList(designBranch){
	var quantifDiv=document.getElementById('quantifListDiv');
	if (!designBranch) { // nothing selected
		quantifDiv.innerHTML='';
		showHideSection('quantif');
		return;
	}
	//wait image
	quantifDiv.innerHTML='&nbsp;&nbsp;&nbsp;<IMG src="$promsPath{images}/scrollbarGreen.gif">';
	quantifDiv.style.display='';
	//If XHR object already exists, the request is canceled & the object is deleted
	if(XHR && XHR.readyState != 0){
		XHR.abort();
		delete XHR;
	}
	//Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}

	XHR.open("GET","$promsPath{cgi}/sequenceView.cgi?AJAX=quantifList&id_prot=$proteinID&designBranch="+designBranch,true);
		XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			quantifDiv.innerHTML=XHR.responseText;
		}
	}
	XHR.send(null);
}


function getXMLHTTP(){
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
function listInteractors(text) {
	document.getElementById('interactDiv').innerHTML=text;
	document.getElementById('showInteract').style.display='none';
	document.getElementById('hideInteract').style.display='';
	interactSearchWasPerformed=true;
}
function listConnectedProteins(text) {
	protDiv=document.getElementById('connProtDiv');
	protDiv.innerHTML=text;
	document.getElementById('showProt').style.display='none';
	document.getElementById('hideProt').style.display='';
	queryWasPerformed=true;
}
function hideConnectedProteins() {
	document.getElementById('connProtDiv').style.display='none';
	document.getElementById('showProt').style.display='';
	document.getElementById('hideProt').style.display='none';
}
function hideInteractors() {
	document.getElementById('interactDiv').style.display='none';
	document.getElementById('showInteract').style.display='';
	document.getElementById('hideInteract').style.display='none';
}
var queryWasPerformed=false;
var interactSearchWasPerformed=false;
function getDomOffset( Obj, Prop ) {
	var iVal = 0;
	while (Obj && Obj.tagName != 'BODY') {
		eval('iVal += Obj.' + Prop + ';');
		Obj = Obj.offsetParent;
	}
	return iVal;
}
function mergePeptideData() {
	var anaIDs = new Array();
	for (var i=0;i<document.anaList.check_ana.length;i++) {
		if (document.anaList.check_ana[i].checked) {
			anaIDs.push(document.anaList.check_ana[i].value);
		}
	}
	if (anaIDs.length > 0) {
		window.location = "$promsPath{cgi}/sequenceView.cgi?id_ana=" + anaIDs.join(',') + "&id_prot=$proteinID&id_type=$idType&msdata=$showMSdata";
	}
}
window.onload=updateSectionVisibility;
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif"> <!--onload="updateSectionVisibility()" onunload="closeSpectWin()" -->
<CENTER>
<FONT class="title">Protein Entry: <FONT color="#DD0000">$alias</FONT></FONT><BR>|;
my $title = (scalar @analysesID == 1)? qq
|Peptide matches for $msTypeNames{$analysisInfo{$analysesID[0]}[6]} ($msTypeNames{"s$analysisInfo{$analysesID[0]}[6]"}) Analysis <FONT color="#DD0000">$analysisInfo{$analysesID[0]}[0]</FONT>.|:
"Peptide matches for selected Analyses.";
print qq|
<FONT class="title2">$title</FONT><BR>
</CENTER>
|;

#############################
####>Protein information<####
#############################
print "<BR><TABLE border=0><TR><TD bgcolor=$color2>\n";
print "<TABLE border=0>\n";
print "<TR><TH align=right valign=top>Name in myProMS :</TH><TD bgcolor=\"$color1\"> ";
my $fWeight=(scalar @analysesID == 1 && $analysisInfo{$analysesID[0]}[3])? 'bold' : 'normal';
my $fColor=($bestConfidence == 1)? 'gray' : 'black';
print "<FONT style=\"font-weight:$fWeight;color:$fColor;\">$alias</FONT>";
print " (Protein is hidden in $msTypeNames{\"s$analysisInfo{$analysesID[0]}[6]\"} Analysis <FONT style=\"font-weight:bold;color:#DD0000;\">$analysisInfo{$analysesID[0]}[0]</FONT>)" unless scalar @analysesID > 1 || $analysisInfo{$analysesID[0]}[3];
print "</TD></TR>\n";
print "<TR><TH align=right>Original identifier :</TH><TD bgcolor=\"$color1\">$identifier&nbsp;</TD></TR>\n";
print "<TR><TH align=right valign=top>Description :</TH><TD bgcolor=\"$color1\">$desProt&nbsp;</TD></TR>\n";
my ($massString,$sizeString)=($protLength)? ("$massProt Da","&nbsp;($protLength aa)") : ('Unknown','');
#print "<TR><TH align=right>Size :</TH><TD bgcolor=$color1>$sizeString</TD></TR>\n";
print "<TR><TH align=right nowrap>&nbsp;Nominal mass (Mr) :</TH><TD bgcolor=\"$color1\">$massString$sizeString&nbsp;</TD></TR>\n";
if (scalar keys %isoforms > 1) {
	print "<TR><TH align=right nowrap>&nbsp;Isoforms :</TH><TD bgcolor=\"$color1\"><SELECT name=\"isoforms\" onchange=\"sequenceView('$strgAnaList',this.value,'$idType',$showMSdata)\">";
	foreach my $protID (sort{&promsMod::sortSmart(lc($isoforms{$a}[1]),lc($isoforms{$b}[1]))} keys %isoforms) {
		my $selectStrg=($protID==$proteinID)? ' selected' : '';
		my $aliasStrg=$isoforms{$protID}[1];
		$aliasStrg.=" ($isoforms{$protID}[0])" if $isoforms{$protID}[0] ne $isoforms{$protID}[1];
		print "<OPTION value=\"$protID\"$selectStrg>$aliasStrg</OPTION>\n";
	}
	print "</SELECT></TD></TR>\n";
}

if ($geneInfo{'GN'}) {
	print "<TR><TH align=right>Gene name :</TH><TD bgcolor=\"$color1\"><B>$geneInfo{GN}[0]</B>";
	if (scalar @{$geneInfo{'GN'}} > 1) {
		print ' (';
		foreach my $i (1..$#{$geneInfo{'GN'}}) {
			print ',' if $i>1;
			print $geneInfo{GN}[$i];
		}
		print ')';
	}
	print "&nbsp;</TD></TR>\n";
}
my $speciesStrg=($geneInfo{'OS'})? "<I>$geneInfo{OS}[0]</I> ($geneInfo{OS}[1])" : "<I>$nameOrganism</I>"; # Modification on 28/11/12 -> Some user defined proteins
print "<TR><TH align=right>Species :</TH><TD bgcolor=$color1>$speciesStrg</TD></TR>\n";
if ($comments) {print "<TR><TH align=right valign=top>Comments :</TH><TD bgcolor=$color1>$comments</TD></TR>\n";}
print "<TR><TH align=right>Last modified :</TH><TD bgcolor=$color1>";
if ($date) {print "$date <B>by</B> $updateUser&nbsp;&nbsp;";} else {print "Never&nbsp;&nbsp;";}
if ($validProtID && $projectAccess ne 'guest' && $projectStatus <= 0) {
	###>Edit button
	print "&nbsp;<INPUT type=\"button\" value=\"Edit Protein\" onclick=\"editProt()\">";
}
print "</TD></TR>\n";
print "</TABLE></TD>\n</TR></TABLE>\n";
if ($bestConfidence == 0) {print "<FONT class=\"title3\"><FONT color='#DD0000'>Warning:</FONT>This protein was identified based on ion profile comparison with other analyses.</FONT><BR>\n";}
elsif ($bestConfidence == 1) {print "<FONT class=\"title3\"><FONT color='#DD0000'>Warning:</FONT>Protein with bad confidence level.</FONT><BR>\n";}


#######################################
####>Aliases of protein in Project<####
#######################################
my $disabAlias=(scalar keys %isoformProteins > 1)? '' : 'disabled=true';
print qq
|<BR>
<TABLE width="100%" bgcolor="$color2"><TR>
<TD><INPUT type="button" class="font11" id="aliasBUTTON" value="Show" onclick="showHideSection('alias')" style="width:55px" $disabAlias/></TD>
<TH bgcolor="$color1" align="left" width="100%" nowrap>&nbsp;Synonyms or isoforms of $alias in Project&nbsp;</TH>
</TR></TABLE>
<DIV id="aliasDIV" style="display:none;margin-left:20px">
<BR style="font-size:5px">
|;
if (scalar keys %isoformProteins > 1) {
	print qq
|<TABLE border=0 cellspacing=0>
<TR bgcolor="$color2">
<TH class="rbBorder" align=left>&nbsp;&nbsp;Protein&nbsp;</TH>
<TH class="rbBorder">&nbsp;MW<FONT style="font-size:9px;"> kDa</FONT>&nbsp;</TH>
<TH class="bBorder" width=600>&nbsp;Description - Species&nbsp;</TH>
</TR>
|;
	my $bgColor=$color1;
	foreach my $protID (sort{lc($isoformProteins{$a}[1]) cmp lc($isoformProteins{$b}[1])} keys %isoformProteins) {
		my ($anaID,$alias,$bestConf,$bestVis,$mw,$des,$org)=@{$isoformProteins{$protID}};
		$mw=($mw)? sprintf "%.1f",$mw/1000 : 0;
		my ($pClass,$pPopup)=&promsMod::getProteinClass($bestConf,$bestVis);
		my $colorStrg=($protID==$validProtID)? ' style="color:#DD0000"' : '';
		print "<TR bgcolor=\"$bgColor\"><TD>&nbsp;<A$colorStrg class=\"$pClass\" href=\"javascript:sequenceView($anaID,$protID,'valid',$showMSdata);\" onmouseover=\"popup('$pPopup')\" onmouseout=\"popout()\">$alias</A>&nbsp;</TD><TD class=\"$pClass\" align=right>&nbsp;$mw&nbsp;</TD><TD>&nbsp;$des <FONT class=\"org\">$org</FONT>&nbsp;</TD></TR>\n";
		$bgColor=($bgColor eq $color1)? $color2 : $color1;
	}

	print "</TABLE>\n<BR>";
}
print "</DIV>\n";

#########################################
####>Analyses where protein is found<####
#########################################
print qq
|<BR style="font-size:5px">
<TABLE width="100%" bgcolor="$color2"><TR>
<TD><INPUT type="button" class="font11" id="analysisBUTTON" value="Show" onclick="showHideSection('analysis')" style="width:55px"/></TD>
<TD bgcolor="$color1" width="100%"><TABLE cellspacing=0 cellpadding=0><TR><TH align="left" nowrap>&nbsp;List of Analyses where $alias is found&nbsp;</TH><TD nowrap><DIV id="anaButtons" style="display:none">
|;
if ($validProtID) {
	print qq
|<INPUT type="button" class="font11" value="Check / Unckeck All" onclick="check(document.anaList.check_ana)"> <INPUT type="button" value="Peptide Distribution" style="font-weight:bold" onclick="graphicalView(document.anaList.check_ana)">|;
}
if ($idType eq 'valid') {
	print qq| <INPUT type="button" value="Cumulate Analyses" style="font-weight:bold" onclick="mergePeptideData();">|;
}

my $mgPopupStrg=qq
|<TABLE><TR><TH colspan=\\'2\\' nowrap>Confidence in visibility assignment</TH></TR><TR><TD class=\\'matchGroupOK\\'>&nbsp;&nbsp;</TD><TD>Reliable</TD></TR><TR><TD class=\\'matchGroupEQ\\'>&nbsp;&nbsp;</TD><TD>Potential ambiguity</TD></TR><TR><TD class=\\'matchGroupNO\\'>&nbsp;&nbsp;</TD><TD>Non-optimal</TD></TR><TR><TD class=\\'matchGroupBAD\\'>&nbsp;&nbsp;</TD><TD>Contradicted by peptide distribution</TD></TR></TABLE>|;
print qq
|</DIV></TD></TR></TABLE></TD></TR></TABLE>
<DIV id="analysisDIV" style="display:none;margin-left:20px">
<FORM name="anaList" method="post" target="ProteinWindow">
<INPUT type="hidden" name="what" value="protein">
<INPUT type="hidden" name="id_item" value="$anaIDstring">
<INPUT type="hidden" name="item" value="project">
<INPUT type="hidden" name="id_prot" value="$validProtID">
<INPUT type="hidden" name="largWindow">
<INPUT type="hidden" name="graphView">
<INPUT type="hidden" name="msdata" value="$showMSdata">
<INPUT type="hidden" name="id_ana" value="$anaIDstring">
<INPUT type="hidden" name="showAllAna" value="$showAllAna">
<BR style="font-size:5px">
<TABLE border=0 cellspacing=0>
<TR bgcolor="$color2">
<TH class="rbBorder" align=left>&nbsp;&nbsp;Parents in Project&nbsp;&nbsp;</TH>
<TH class="rbBorder" align=left>&nbsp;&nbsp;Analysis&nbsp;&nbsp</TH>
<TH class="rbBorder">&nbsp;&nbsp;Peptides&nbsp;&nbsp;</TH>
<TH class="rbBorder" nowrap>&nbsp;&nbsp;Spec.<FONT style="font-size: 9px;"> (%)</FONT>&nbsp;&nbsp;</TH>
<TH class="rbBorder" nowrap>&nbsp;&nbsp;Cov.<FONT style="font-size: 9px;"> (%)</FONT>&nbsp;&nbsp;</TH>
<TH class="bBorder" nowrap>&nbsp;&nbsp;&nbsp;<A href="void(null)" onmouseover="popup('$mgPopupStrg')" onmouseout="popout()">Match group<SUP>*</SUP></A>&nbsp;&nbsp;&nbsp</TH>
<TH class="bBorder">&nbsp</TH>
</TR>
|;
my $bgColor=$color2;
my %itemIcones=&promsConfig::getItemIcones;
my $numValAnalyses=0;
my ($prevExpName,$prevGel2dName,$prevSpotName,$prevSampName)=('','','','');
my $anaNamesString='';
foreach my $anaID (sort{$analysisParents{$a}[0]<=>$analysisParents{$b}[0] || $analysisParents{$a}[2]<=>$analysisParents{$b}[2] || $analysisParents{$a}[4]<=>$analysisParents{$b}[4] || lc($analysisParents{$a}[5]) cmp lc($analysisParents{$b}[5]) || $analysisInfo{$a}[1]<=>$analysisInfo{$b}[1]} keys %analysisInfo) {
	##>Parents
	my $parentStrg;
	if ($analysisParents{$anaID}[0]==1) {# Gel
		#>Exp
		if ($analysisParents{$anaID}[1] eq $prevExpName) { # same Exp
			$parentStrg="<FONT color=\"$bgColor\">$analysisParents{$anaID}[1]&nbsp;&nbsp;&gt;&nbsp;&nbsp;";
			#>Gel
			if ($analysisParents{$anaID}[3] eq $prevGel2dName) {
				$parentStrg.="$analysisParents{$anaID}[3]&nbsp;&nbsp;&gt;&nbsp;&nbsp;";
				#>Spot
				if ($analysisParents{$anaID}[5] eq $prevSpotName) {$parentStrg.="$analysisParents{$anaID}[5]&nbsp;&nbsp;&gt;&nbsp;&nbsp;</FONT>";}
				else {$parentStrg.="</FONT>$analysisParents{$anaID}[5]&nbsp;&nbsp;&gt;&nbsp;&nbsp;";}
			}
			else {$parentStrg.="</FONT>$analysisParents{$anaID}[3]&nbsp;&nbsp;&gt;&nbsp;&nbsp;$analysisParents{$anaID}[5]&nbsp;&nbsp;&gt;&nbsp;&nbsp;";}
		}
		else { # diff Exp
			$parentStrg="$analysisParents{$anaID}[1]&nbsp;&nbsp;&gt;&nbsp;&nbsp;$analysisParents{$anaID}[3]&nbsp;&nbsp;&gt;&nbsp;&nbsp;$analysisParents{$anaID}[5]&nbsp;&nbsp;&gt;&nbsp;&nbsp;";
			$bgColor=($bgColor eq $color1)? $color2 : $color1;
		}
		$prevExpName=$analysisParents{$anaID}[1];
		$prevGel2dName=$analysisParents{$anaID}[3];
		$prevSpotName=$analysisParents{$anaID}[5];
	}
	else {# No Gel
		#>Exp
		if ($analysisParents{$anaID}[1] eq $prevExpName) { # same Exp
			$parentStrg="<FONT class=\"hidden\">$analysisParents{$anaID}[1]&nbsp;&nbsp;&gt;&nbsp;&nbsp;";
			#>Samp
			if ($analysisParents{$anaID}[3] eq $prevSampName) {$parentStrg.="$analysisParents{$anaID}[3]&nbsp;&nbsp;&gt;&nbsp;&nbsp;</FONT>";}
			else {$parentStrg.="</FONT>$analysisParents{$anaID}[3]&nbsp;&nbsp;&gt;&nbsp;&nbsp;";}
		}
		else { # diff Exp
			$parentStrg="$analysisParents{$anaID}[1]&nbsp;&nbsp;&gt;&nbsp;&nbsp;$analysisParents{$anaID}[3]&nbsp;&nbsp;&gt;&nbsp;&nbsp;";
			$bgColor=($bgColor eq $color1)? $color2 : $color1;
		}
		$prevExpName=$analysisParents{$anaID}[1];
		$prevSampName=$analysisParents{$anaID}[3];
	}
	print "<TR bgcolor=\"$bgColor\" class=\"list\"><TH nowrap align=left valign=top>&nbsp;&nbsp;$parentStrg</TH>";

	##>Analysis
	#print '<TD nowrap>';
	my $anaTypeCode=$msTypeNames{"s$analysisInfo{$anaID}[6]"};
	my $valStat=$analysisInfo{$anaID}[7];
	my $anaCode=($valStat==-1)? 'analysis:no_scan' : ($valStat==0)? 'analysis:no_val' : ($valStat==1)? 'analysis:part_val' : 'analysis:val';
	#print "<IMG src=\"$promsPath{images}/$itemIcones{$anaCode}\">&nbsp;";
	my $iconStrg="<IMG src=\"$promsPath{images}/$itemIcones{$anaCode}\">&nbsp;";
	my $ptmExtraString='';
	if ($protVarMods{$anaID}) {$ptmExtraString=&convertVarMods(\%{$protVarMods{$anaID}});}
	my ($protID,$typeID)=($analysisInfo{$anaID}[8])? ($validProtID,'valid') : ($analysisInfo{$anaID}[12],'notValid');
	#if (grep {$_ == $anaID} @analysesID) {
	#	if ($analysisInfo{$anaID}[3] > 0) {print "<B><FONT color=#DD0000>$analysisInfo{$anaID}[0]</FONT>$ptmExtraString</B>\n";} # protein is visible
	#	else {print "<FONT color=#DD0000>$analysisInfo{$anaID}[0]</FONT>$ptmExtraString\n";}
	#	$anaNamesString .= ', ' if $anaNamesString;
	#	$anaNamesString .= $analysisInfo{$anaID}[0];
	#}
	#else {
	#	my ($pClass,$pPopup)=&promsMod::getProteinClass($analysisInfo{$anaID}[2],$analysisInfo{$anaID}[3]); # confidence,visibility
	#	my ($protID,$typeID)=($analysisInfo{$anaID}[8])? ($validProtID,'valid') : ($analysisInfo{$anaID}[12],'notValid');
	#	print "<A class=\"$pClass\" href=\"javascript:sequenceView($anaID,$protID,'$typeID',$showMSdata);\" onmouseover=\"popup('$pPopup')\" onmouseout=\"popout()\">$analysisInfo{$anaID}[0]$ptmExtraString</A>";
	#}
	#print "&nbsp;&nbsp;</TD>\n";
	my $anaClassStrg='';
	if (grep {$_ == $anaID} @analysesID) { # in list of cummlated ana
		$anaClassStrg.='selectedAna'; # if grep {$_ == $anaID} @analysesID;
		$anaNamesString .= ', ' if $anaNamesString;
		$anaNamesString .= $analysisInfo{$anaID}[0];
	}
	$anaClassStrg.=($analysisInfo{$anaID}[3] > 0)? ' visibleProt' : ' hiddenProt'; # protein is visible or hidden
	if (grep {$_ == $anaID} @analysesID) { # selected analysis
		$anaClassStrg.=' virtualProt' if $analysisInfo{$anaID}[2]==0; # confidence (selection is not compatible with &promsMod::getProteinClass for virtual prot)
		if (scalar @analysesID > 1) {
			print "<TD nowrap>$iconStrg<A class=\"$anaClassStrg\" href=\"javascript:sequenceView($anaID,$protID,'$typeID',$showMSdata);\">$analysisInfo{$anaID}[0]$ptmExtraString</A>"; # keep single ana link
		}
		else {
			print "<TD nowrap class=\"$anaClassStrg\">$iconStrg$analysisInfo{$anaID}[0]$ptmExtraString"; # no single ana link
		}
	}
	else {
		my ($pClass,$pPopup)=&promsMod::getProteinClass($analysisInfo{$anaID}[2],$analysisInfo{$anaID}[3]); # confidence,visibility
		print "<TD nowrap>$iconStrg<A class=\"$pClass\" href=\"javascript:sequenceView($anaID,$protID,'$typeID',$showMSdata);\" onmouseover=\"popup('$pPopup')\" onmouseout=\"popout()\">$analysisInfo{$anaID}[0]$ptmExtraString</A>";
	}
	print "&nbsp;&nbsp;</TD>\n";

	##>Peptides
	print "<TH>$analysisInfo{$anaID}[9]";
	if ($valStat<2 && $analysisInfo{$anaID}[2] && $projectAccess=~/bioinfo|mass|manag|super/) { # can enter in validation mode (not for virtual prot!)
		print "&nbsp;<A href=\"javascript:startValidation($anaID);\" onmouseover=\"popup('Click to access all identification data of <B>$identifier</B> in this Analysis')\" onmouseout=\"popout()\"><IMG class=\"button\" src=\"$promsPath{images}/good.gif\"></A>";
	}
	print "</TH>\n";

	##>Specificity
	if ($valStat>=1) {
		if (defined $analysisInfo{$anaID}[10]) {print "<TH>$analysisInfo{$anaID}[10]</TH>";} else {print "<TH>?</TH>";}
	}
	else {print "<TH>-</TH>";}

	##>Coverage
	if ($valStat>=1) {
		my $pepCov=($analysisInfo{$anaID}[11]=~/^-([\d\.]+)/)? "&lt;$1" : $analysisInfo{$anaID}[11];
		print "<TH>$pepCov</TH>";
	}
	else {print "<TH>-</TH>";}

	##>Match group
	if (!$analysisInfo{$anaID}[5]) { # match group
		print "<TH>-</TH>\n";
	}
	elsif ($analysisInfo{$anaID}[5] eq 'none') {
		print "<TH class=\"$analysisInfo{$anaID}[12]\">&nbsp;<A onmouseover=\"popup('Protein does not belong to a <B>Match group</B> in $anaTypeCode Analysis <B>$analysisInfo{$anaID}[0]</B>')\" onmouseout=\"popout()\">none</A></TH>\n";
	}
	else {
		if ($projectAccess eq 'guest') {
			print "<TH class=\"$analysisInfo{$anaID}[12]\">&nbsp;$analysisInfo{$anaID}[5]</TH>\n";
		}
		else {
			print "<TH class=\"$analysisInfo{$anaID}[12]\">&nbsp;<A href=\"javascript:editMatchGroup($anaID,$analysisInfo{$anaID}[4])\" onmouseover=\"popup('Click to <B>edit</B> Match group in $anaTypeCode Analysis <B>$analysisInfo{$anaID}[0]</B>')\" onmouseout=\"popout()\">$analysisInfo{$anaID}[5]</A>";
			print "<IMG border=0 align=top id=\"img_$proteinID:$anaID\" src=\"$promsPath{images}/plus.gif\" onclick=\"ajaxGetProteinsMatchGroup($anaID)\" onmouseover=\"popup('Click to <B>display</B> Match group in $anaTypeCode Analysis <B>$analysisInfo{$anaID}[0]</B>')\" onmouseout=\"popout()\"></TH>\n";
		}
	}
	if ($analysisInfo{$anaID}[8]) { # protein is validated in Ana
		my $checkString=(grep {$_ == $anaID} @analysesID)? 'checked' : '';
		print "<TH><INPUT type=\"checkbox\" name=\"check_ana\" value=\"$anaID\" $checkString></TH>\n";
		$numValAnalyses++;
	}
	else {print "<TH></TH>\n";}
	print "</TR>\n"; # end of row
}

$bgColor=($bgColor eq $color1)? $color2 : $color1;
print "<TR bgcolor=\"$bgColor\"><TH align=\"left\" colspan=\"7\">";
if ($showAllAna) {
	print "&nbsp;End of list&nbsp;<INPUT type=\"button\" value=\"Display partial list\" onclick=\"updateAnalysesList(0)\">";
}
else {
	print qq "<FONT color=\"#DD0000\">&nbsp;... Warning: Partial list</FONT>&nbsp;<INPUT type=\"button\" value=\"Display full list\" onclick=\"updateAnalysesList(1)\">";
}
print qq
|</TH></TR>
</TABLE>
</FORM>
</DIV>
|;

###############################################
####>Quantification where protein is found<####
###############################################
print qq
|<BR style="font-size:5px">
<TABLE width="100%" bgcolor="$color2"><TR>
<TD><INPUT type="button" id="quantifBUTTON" class="font11" value="Show" onclick="showHideSection('quantif')" style="width:55px" /></TD>
<TH bgcolor="$color1" align="left" width="100%" nowrap>&nbsp;List of Quantifications where $alias is found&nbsp;</TH>
</TR></TABLE>
<DIV id="quantifDIV" style="display:none;margin-left:20px">
<BR style="font-size:5px">
<TABLE><TR>
<TD valign="top"><TABLE><TR>
	<TH bgcolor="$color2" nowrap>&nbsp;&nbsp;Select:<SELECT name="design" onchange="ajaxGetQuantificationsList(this.value)"><OPTION value="">-= Select =-</OPTION>
	<OPTGROUP label="-=Design-based quantifications=-">
|;
my $prevExpPos=0;
foreach my $designID (sort{$designs{$a}[0]<=>$designs{$b}[0] || lc($designs{$a}[2]) cmp lc($designs{$b}[2])} keys %designs) {
	my ($expPos,$expName,$desName)=@{$designs{$designID}};
	if ($prevExpPos != $expPos) {
		print "</OPTGROUP>\n" if $prevExpPos;
		print "<OPTGROUP label=\"Exp. $expName: \">\n";
		$prevExpPos=$expPos;
	}
	print "<OPTION value=\"design:$designID\">$desName</OPTION>\n";
}
if ($prevExpPos) {print "</OPTGROUP>\n";}
else {print "<OPTION disabled>None</OPTION>\n";}

print qq
|</OPTGROUP>
<OPTGROUP label="-=Internal-Analysis quantifications=-">
<OPTION value="analysis:$strgAnaList">List from selected Analyses</OPTION>
<OPTION value="analysis:0">List from all Analyses</OPTION>
</OPTGROUP>
</SELECT></TD></TR></TABLE></TD>
<TD><DIV id="quantifListDiv"></DIV></TD>
</TR></TABLE>
</DIV>
|;

############################
####>Connected proteins<####
############################
#my $protStrg=($numConnProt<=1)? 'protein' : 'proteins';
my $disabConnProt=($currProtStatus)? '' : 'disabled=true';
print qq
|<BR style="font-size:5px">
<TABLE width="100%" bgcolor="$color2"><TR>
<TD><INPUT type="button" id="showProt" class="font11" value="Show" onclick="ajaxGetConnectedProteins($validProtID)" style="width:55px" $disabConnProt/>
	<INPUT type="button" id="hideProt" class="font11" value="Hide" style="width:55px; display:none" onclick="hideConnectedProteins()"/>
</TD>
<TH bgcolor="$color1" align="left" width="100%" nowrap>&nbsp;Check if $alias shares peptides with other proteins in Project&nbsp;</TH>
</TR></TABLE>
<DIV id="connProtDiv" style="display:none;margin-left:20px"></DIV>
|;


##################
####>PTM list<####
##################
my $disabPTM=(scalar @pep_beg)? '' : 'disabled=true'; # not for virtual proteins
#if (scalar @pep_beg) { # not for virtual proteins
print qq
|<BR style="font-size:5px">
<TABLE width="100%" bgcolor="$color2"><TR>
<TD><INPUT type="button" class="font11" id="ptmBUTTON" value="Show" onclick="showHideSection('ptm')" style="width:55px" $disabPTM/></TD>
<TH bgcolor="$color1" align="left" width="100%" nowrap>&nbsp;Post-translational modifications relevant to Project&nbsp;</TH>
</TR></TABLE>
<DIV id="ptmDIV" style="display:none;margin-left:20px">
<BR style="font-size:5px">
|;
if ($varModInProject) {
	$bgColor = $color2;
	print qq|<TABLE border=0 cellspacing=0>
	<TR bgcolor="$bgColor">
	<TH class="rbBorder">&nbsp;Modification type&nbsp;</TH><TH class=\"rbBorder\">&nbsp;AA&nbsp;</TH><TH class="bBorder">&nbsp;Positions&nbsp;</TH></TR>|;
	my $matchedPTM=0;
	foreach my $varModCode (sort keys %varModList){
		next unless $projectVarMods{$varModCode};
		$matchedPTM++;
		my ($varModName,$varModTagStrg)=($allPostTransModifs{$varModCode})? ($allPostTransModifs{$varModCode}[0]," [<FONT class=\"$allPostTransModifs{$varModCode}[2]\">$allPostTransModifs{$varModCode}[1]</FONT>]") : ($varModCode,'');
		foreach my $residue (sort keys %{$varModList{$varModCode}}){
			my $max = max values %{$varModList{$varModCode}{$residue}};
			next unless $max>0; # All positions for this residue are ambiguous -> Do not need to print it
			my @posList;
			$bgColor = ($bgColor eq $color2) ? $color1 : $color2;
			print "<TR bgcolor=\"$bgColor\"><TH align=\"right\">&nbsp;$varModName$varModTagStrg:&nbsp;</TH><TH>&nbsp;$residue&nbsp;</TH>";
			foreach my $pos (sort {$a <=> $b} keys %{$varModList{$varModCode}{$residue}}){
				next unless $varModPosOnProtein{$pos}{$varModCode} != -1;
				my $countMatch = $varModList{$varModCode}{$residue}{$pos};
				if($pos ==0){ $pos = 'N-ter' } elsif($pos == 999999){ $pos = 'C-ter' };
				if($countMatch >= 2){
					$pos.="<SUP>$countMatch</SUP>";
				}
				push @posList, $pos;
			}
			print "<TH align=\"left\">&nbsp;";
			print join ",", @posList if(scalar @posList);
			print "&nbsp;</TH></TR>\n";
		}
	}
	print "<TR bgcolor=\"$color1\"><TH align=\"left\" colspan=3>&nbsp;No PTMs matched&nbsp;</TH></TR>\n" unless $matchedPTM;
	print "</TABLE>\n";
}
else {print "None.<BR>\n";}
print "<BR></DIV>\n";
#}

#######################################
####>Sequence with peptide matches<####
#######################################
my $anaText=(scalar @analysesID > 1)? 'selected Analyses' : "Analysis <FONT color=\"#DD0000\">$anaNamesString</FONT>";
print qq
|<BR style="font-size:5px">
<TABLE width="100%" bgcolor="$color2"><TR>
<TD><INPUT type="button" class="font11" id="sequenceBUTTON" value="Show" onclick="showHideSection('sequence')" style="width:55px"/></TD>
<TH bgcolor="$color1" align="left" width="100%" nowrap>&nbsp;Detailed sequence coverage in $anaText&nbsp;</TH>
</TR></TABLE>
<DIV id="sequenceDIV" style="display:none;margin-left:20px;">
<BR style="font-size:5px">
|;

if ($seqMismatch) {
	print "<B><FONT color=\"#DD0000\">Warning:</FONT> ";
	if ($protLength) {print 'Peptide/protein sequence mismatch detected (isoform?).';}
	else {print 'Protein sequence not recorded.'}
	if ($covProtLength) {print ' Peptide-based sequence is displayed.</B><BR>';}
	elsif ($protSeq) {print ' Sequence from canonical entry is displayed.</B><BR>';} # masterProt
}
my $coveredSeq='';
#if ($protLength) {
my $textBoxHeight=int((length($protSeq)/60) + 0.5);
print qq
|(Matching peptides are shown in <B><FONT color="red">bold red</FONT></B>, overlapping peptides in <B><FONT color="blue">bold blue</FONT></B>).
&nbsp;<INPUT type="button" value="Extract sequence" onclick="extractSequence(protSequence)" class="font11"/><INPUT type="button" value="Extract covered sequence" onclick="extractSequence(coveredSequence)" class="font11"/>
<DIV id="extractionDIV" class="popup">
<INPUT type="button" class="font11" value="Close" onclick="document.getElementById('extractionDIV').style.display='none'"/><BR>
<TEXTAREA id="extractionTEXT" name="extractionTEXT" rows=$textBoxHeight cols=60></TEXTAREA>
</DIV>
<TABLE id="covTABLE" border=0 cellpadding=0 cellspacing=0>
|;
my $index=-1;
#my $aaColor='black';
my $aaClass='overlapRes0';
my @protResidues=split(//,$protSeq);
my $lineLastPosIdx;
my $cumulate = (scalar @analysesID > 1)? 1:0;
foreach my $posIdx (0..$#protResidues) {
	my $resPos=$posIdx+1;
	my $newLine=0;
	if ($posIdx/60==int($posIdx/60)) {# new line
		$newLine=1;
		$lineLastPosIdx=$posIdx+59; $lineLastPosIdx=$#protResidues if $lineLastPosIdx > $#protResidues;
		#print "</B></FONT>" if $aaColor ne 'black';
		#print "</FONT>";
		#print "</TD></TR>\n" if $posIdx>0;
		#>N-ter varMod
		my $NterModStrg='';
		my $colSpan=1;
		if ($varModPosOnProtein{0}) {
			my $varModTagStrg='';
			if ($posIdx==0) {
				my @allVarModOnPos=keys %{$varModPosOnProtein{0}};
				if(keys %{$varModPosOnProtein{0}} == 1){
					my $varModCode=$allVarModOnPos[0];
					my $dispVarMod=($projectVarMods{$varModCode})? $allPostTransModifs{$varModCode}[1] : '*';
					$varModTagStrg=($allPostTransModifs{$varModCode})? "&nbsp;<FONT class=\"$allPostTransModifs{$varModCode}[2]\" onmouseover=\"popup('<B>$allPostTransModifs{$varModCode}[0] (Protein N-term)</B>')\" onmouseout=\"popout()\">$dispVarMod</FONT>~" : "&nbsp;<FONT onmouseover=\"popup('<B>$varModCode (Protein N-term)</B>')\" onmouseout=\"popout()\">$dispVarMod</FONT>~";
				}else{# if more than one modification -> print &
					$varModTagStrg="&nbsp;<FONT onmouseover=\"popup('<B> Multiple PTM on Protein N-term</B>')\" onmouseout=\"popout()\">&</FONT>~" ;
				}
			}
			$NterModStrg="<TD align=right class=\"seq\">$varModTagStrg</TD>";
			$colSpan=2;
		}
		#>VarMod line (Checking if varMod within the next 60 res)
		my $lastModPos=0;
		foreach my $modPos (sort{$a<=>$b} keys %varModPosOnProtein) {
			if ($modPos >= $resPos && $modPos < $resPos+60) {$lastModPos=$modPos;}
			elsif ($modPos >= $resPos+60) {last;}
		}
		print "<TR><TD colspan=$colSpan></TD><TD class=\"seq\"><BR>";
		my $curPos=$resPos; # pos on seq
		my $disPos=$curPos; # pos on displayed seq (with spaces)
		while ($curPos <= $lastModPos) { # no loop if no varMod on seq line
			print '&nbsp;' if ($curPos>$resPos && ($disPos-1)/10==int(($disPos-1)/10)); # space in seq => shift +1 space
			if ($varModPosOnProtein{$curPos}) {
				my $varModCode;
				my @allVarModOnPos=keys %{$varModPosOnProtein{$curPos}};
				if(keys %{$varModPosOnProtein{$curPos}} == 1){$varModCode=$allVarModOnPos[0];}
				else{$varModCode = 'multi';}
				my $dispVarMod='*'; # default
				my $jsActionStrg='javascript:void(null)'; # default
				my $popupInfoStrg=''; # default
				if ($projectVarMods{$varModCode}) {
					$dispVarMod=$allPostTransModifs{$varModCode}[1] or die Dumper(%allPostTransModifs);
					if ($currProtStatus==1) { # validated prot only
						$jsActionStrg="javascript:ajaxGetVarModPeptideInfo('$varModCode',$curPos,'".join(',',@{$varModPosPeptides{$curPos}{$varModCode}})."',$cumulate)"; # list matching peptides id (string)
						$popupInfoStrg='<BR>Click for more details.';
					}
				}
				elsif ($varModCode eq 'multi'){
					$dispVarMod='&';
					$varModCode='';
					$popupInfoStrg='Multiple modifications<BR>Click for more details.';
					my ($ajaxVarModCodes,$ajaxPeptides)=('','');
					foreach my $specVarMod (sort{$a cmp $b} keys %{$varModPosOnProtein{$curPos}}) {
						next unless $varModPosPeptides{$curPos}{$specVarMod};
						if(!$ajaxVarModCodes){
							$ajaxVarModCodes="$specVarMod";
							$ajaxPeptides=join(',',@{$varModPosPeptides{$curPos}{$specVarMod}});
						}
						else{
							$ajaxVarModCodes.=",,$specVarMod";
							$ajaxPeptides.=",,".join(',',@{$varModPosPeptides{$curPos}{$specVarMod}});
						}
					}
					$jsActionStrg="javascript:ajaxGetVarModPeptideInfo('$ajaxVarModCodes',$curPos,'$ajaxPeptides',$cumulate)"; # list matching peptides id (string)
				}
				#>check if enough space to display varMod
				my $vModExtraChars=length($dispVarMod)-1;
				my $numSpaces=0;
				if ($vModExtraChars) { # varMod is more than 1 letter-long
					my $numVisChars=1;
					for (my $i=$curPos+1; $i<=$curPos+$vModExtraChars; $i++) {
						if (($i-1)/10==int(($i-1)/10)) { # space in seq => 1 more char allowed
							$numVisChars++;
							$numSpaces++;
						}
						last if ($varModPosOnProtein{$i} && $i<=$lineLastPosIdx+1); # next varMod found current line
						$numVisChars++;
					}
					$dispVarMod=~s/^(.{$numVisChars}).*/$1/; # take only letter(s) that can fit
				}
				my $modResIdx=$curPos-1;
				#my $varModTagStrg=($varModPosOnProtein{$curPos}{$varModCode} && $varModPosOnProtein{$curPos}{$varModCode}<0)? "<A id=\"vmod_$curPos\" class=\"$allPostTransModifs{$varModCode}[2]\" href=\"$jsActionStrg\" onmouseover=\"popup('<B>$allPostTransModifs{$varModCode}[0] -> Ambiguous position for modification<BR>Click for more details.</B>')\" onmouseout=\"popout()\">-</A>" :($allPostTransModifs{$varModCode})? "<A id=\"vmod_$curPos\" class=\"$allPostTransModifs{$varModCode}[2]\" href=\"$jsActionStrg\" onmouseover=\"popup('<B>$allPostTransModifs{$varModCode}[0] ($protResidues[$modResIdx]:$curPos)$popupInfoStrg</B>')\" onmouseout=\"popout()\">$dispVarMod</A>" : "<A id=\"vmod_$curPos\" href=\"$jsActionStrg\" onmouseover=\"popup('<B>$varModCode$popupInfoStrg</B>')\" onmouseout=\"popout()\">$dispVarMod</A>";
				my $varModTagStrg=($varModPosOnProtein{$curPos}{$varModCode} && $varModPosOnProtein{$curPos}{$varModCode}<0 && $allPostTransModifs{$varModCode})? "<A id=\"vmod_$curPos\" class=\"$allPostTransModifs{$varModCode}[2]\" href=\"$jsActionStrg\" onmouseover=\"popup('<B>$allPostTransModifs{$varModCode}[0] -> Ambiguous position for modification<BR>Click for more details.</B>')\" onmouseout=\"popout()\">?</A>" :($allPostTransModifs{$varModCode})? "<A id=\"vmod_$curPos\" class=\"$allPostTransModifs{$varModCode}[2]\" href=\"$jsActionStrg\" onmouseover=\"popup('<B>$allPostTransModifs{$varModCode}[0] ($protResidues[$modResIdx]:$curPos)$popupInfoStrg</B>')\" onmouseout=\"popout()\">$dispVarMod</A>" : ($varModNames{$varModCode})? "<A id=\"vmod_$curPos\" href=\"$jsActionStrg\" onmouseover=\"popup('<B>$varModNames{$varModCode}  ($protResidues[$modResIdx]:$curPos)$popupInfoStrg</B>')\" onmouseout=\"popout()\">$dispVarMod</A>" :  "<A id=\"vmod_$curPos\" href=\"$jsActionStrg\" onmouseover=\"popup('<B>$popupInfoStrg</B>')\" onmouseout=\"popout()\">$dispVarMod</A>";
				print $varModTagStrg;
				$curPos+=$vModExtraChars-$numSpaces;
				$disPos+=$vModExtraChars;
			}
			else {print '&nbsp;'}
			$curPos++;
			$disPos++;
			#print '!' if ($curPos>$resPos && ($curPos-1)/10==int(($curPos-1)/10)); # space in seq => shift +1 space
		}
		print "</TD></TR>\n";
		#>Sequence line
		print "<TR><TD align=right>$resPos&nbsp;</TD>$NterModStrg<TD class=\"seq\">";
		#print "<FONT color=$aaColor><B>" if $aaColor ne 'black';
		#print "<FONT class=\"$aaClass\">";
	}
	elsif ($posIdx>0 && $posIdx/10==int($posIdx/10)) { # space
		print '&nbsp;';
	}

	my @allVarModOnPos=keys %{$varModPosOnProtein{$resPos}};
	my $residueStrg=($allVarModOnPos[0])? "<FONT class=\"modRes\">$protResidues[$posIdx]</FONT>" : $protResidues[$posIdx];

	##>Boundary residue
	if (defined($boundaryStatus{$posIdx})) {
		#my $printResidue;
		if ($posIdx==0) { # first residue
			$aaClass=($boundaryStatus{$posIdx} < 2)? "overlapRes$boundaryStatus{$posIdx}" : 'overlapRes2';
			print "<FONT class=\"$aaClass\">$residueStrg";
			$coveredSeq.=($aaClass ne 'overlapRes0')? $protResidues[$posIdx] : 'X';
		}
		elsif ($boundaryStatus{$posIdx}==1) { # entering simple match (or rarely 1 aa overlap!)
			if ($boundaryStatus{$boundaryList[$index]}>=2) { # leaving an overlap
				#$aaColor='red';
				#$printResidue="$protResidues[$posIdx]</B></FONT><FONT color=$aaColor><B>" ; # leaving an overlap
				$aaClass='overlapRes1';
				print "$residueStrg</FONT>";
				print "<FONT class=\"$aaClass\">" unless $posIdx==$lineLastPosIdx;
			}
			elsif ($boundaryStatus{$boundaryList[$index]}==0) { # leaving unmatched sequence
				#$aaColor='red';
				#$printResidue="<FONT color=$aaColor><B>$protResidues[$posIdx]";
				$aaClass='overlapRes1';
				print '</FONT>' unless $newLine;
				print "<FONT class=\"$aaClass\">$residueStrg";
			}
			else { # 1 aa overlap ($boundaryStatus{$boundaryList[$index]}=1)!!!
				#$aaColor='red';
				#$printResidue="</B></FONT><FONT color=blue><B>$protResidues[$posIdx]</B></FONT><FONT color=$aaColor><B>";
				$aaClass='overlapRes1';
				print '</FONT>' unless $newLine;
				print "<FONT class=\"overlapRes2\">$residueStrg</FONT>";
				print "<FONT class=\"$aaClass\">" unless $posIdx==$lineLastPosIdx;
			}
			$coveredSeq.=$protResidues[$posIdx];
		}
		elsif ($boundaryStatus{$posIdx}>=2) { # entering overlap
			if ($boundaryStatus{$boundaryList[$index]}>=2) { # leaving/entering other overlap level
				print "<FONT class=\"$aaClass\">" if $newLine;
				#$printResidue=$protResidues[$posIdx];
				print $residueStrg;
			}
			elsif ($boundaryStatus{$boundaryList[$index]}==1) { # leaving a simple match
				#$aaColor='blue';
				#$printResidue="</B></FONT><FONT color=$aaColor><B>$protResidues[$posIdx]";
				$aaClass='overlapRes2';
				print '</FONT>' unless $newLine;
				print "<FONT class=\"$aaClass\">$residueStrg";
			}
			else { # leaving unmatched sequence
				#$aaColor='blue';
				#$printResidue="<FONT color=$aaColor><B>$protResidues[$posIdx]";
				$aaClass='overlapRes2';
				print '</FONT>' unless $newLine;
				print "<FONT class=\"$aaClass\">$residueStrg";
			}
			$coveredSeq.=$protResidues[$posIdx];
		}
		else { # entering an unmatched sequence (next residue) ($boundaryStatus{$posIdx}=0)
			#$printResidue="$protResidues[$posIdx]</B></FONT>";
			#$aaColor='black'; # applies to next residue
			$aaClass='overlapRes0'; # applies to next residue
			print "$residueStrg</FONT><FONT class=\"$aaClass\">";
			$coveredSeq.=$protResidues[$posIdx];
		}
		#print "$printResidue";
		$index++;
	}
	##>Normal residue
	else { # Not a boundary residue
		print "<FONT class=\"$aaClass\">" if $newLine;
		print $residueStrg;
		#$coveredSeq.=($aaColor ne 'black')? $protResidues[$posIdx] : 'X';
		$coveredSeq.=($aaClass ne 'overlapRes0')? $protResidues[$posIdx] : 'X';
		#print "</FONT></TD></TR>\n" if $posIdx==$lineLastPosIdx;
		#next;
	}

	##>End of line
	if ($posIdx==$lineLastPosIdx) {
		print '</FONT>';

		##>Cter modification
		if ($posIdx==$#protResidues && $varModPosOnProtein{999999}) {
			my @allVarModOnPos=keys %{$varModPosOnProtein{999999}};
			my $varModCode=$allVarModOnPos[0];
			my $dispVarMod=($projectVarMods{$varModCode})? $allPostTransModifs{$varModCode}[1] : '*';
			my $varModTagStrg=($allPostTransModifs{$varModCode})? "~<FONT class=\"$allPostTransModifs{$varModCode}[2]\" onmouseover=\"popup('<B>$allPostTransModifs{$varModCode}[0] (Protein C-term)</B>')\" onmouseout=\"popout()\">$dispVarMod</FONT>" : "~<FONT onmouseover=\"popup('<B>$varModCode (Protein C-term)</B>')\" onmouseout=\"popout()\">$dispVarMod</FONT>";
			print $varModTagStrg;
		}

		print "</TD></TR>\n";
	}
}
my $bestSpecificity;
foreach my $anaID (@analysesID){
	if ($analysisInfo{$anaID}[10] && (!$bestSpecificity || $analysisInfo{$anaID}[10] > $bestSpecificity)) {
		$bestSpecificity = $analysisInfo{$anaID}[10];
	}
}
my $specificityStrg=($bestSpecificity)? "$bestSpecificity %" : 'Not recorded';
$bestProtScore='N/A' unless $bestConfidence;
print qq
|</TABLE><BR>
<SCRIPT LANGUAGE="JavaScript">
var protSequence='$protSeq';
var coveredSequence='$coveredSeq';
</SCRIPT>|;
if (scalar @analysesID ==1) {
	print "<B>Protein score:</B> $bestProtScore<BR>\n";
}
else {
	print "<B>Best protein score:</B> $bestProtScore<BR>\n";
}
print qq
|<B>Best peptide specificity:</B> $specificityStrg<BR>
<B>Peptide coverage:</B>
|;
if ($pepCoverage) {
	if ($pepCoverage < 0) {print ' &lt;',abs($pepCoverage)," %<BR>\n";} # seq error flag
	else {print " $pepCoverage %<BR>\n";}
}
if (scalar keys %unmatchedPeptides) {
	if ($pepCoverage) {print " &gt; $pepCoverage %";} else {print ' N/A';}
	print " (<FONT style=\"font-weight:bold;color:#DD0000;\">",scalar keys %unmatchedPeptides," peptide(s) with unknown position</FONT>).<BR>\n";
}
elsif (!$pepCoverage) {print " N/A (<FONT style=\"font-weight:bold;color:#DD0000;\">No peptides selected</FONT>).<BR>\n";}
#else {
#	print qq
#|<B>Unknown sequence.</B><BR>
#<B>Score:</B> $bestProtScore<BR>
#<B>Peptide specificity:</B> $analysisInfo{$analysisID}[10] %<BR>
#<B>Peptide coverage:</B> Unknown<BR>
#|;
#}
print "<BR></DIV>\n";

#####################################
####>List of peptides identified<####
#####################################
my $disabPeptide=(scalar @pep_beg)? '' : 'disabled=true'; # not for virtual proteins
print qq
|<BR style="font-size:5px">
<TABLE width="100%" bgcolor="$color2"><TR>
<TD><INPUT type="button" class="font11" id="peptideBUTTON" value="Show" onclick="showHideSection('peptide')" style="width:55px" $disabPeptide/></TD>
<TD bgcolor="$color1" width="100%"><TABLE cellspacing=0 cellpadding=0><TR>
<TH align="left" nowrap>&nbsp;Peptide list&nbsp;</TH>
|;
if (scalar @pep_beg) {
	my ($msData,$image)=($showMSdata)? (-1,'minus1.gif') : (2,'plus.gif');
	print qq|<TD><IMG src="$promsPath{images}/$image" onclick="sequenceView('$anaIDstring',$proteinID,'$idType',$msData)" onmouseover="popup('<B>Click to show/hide MS data & hidden peptides if any</B>')" onmouseout="popout()"></TD>|;
}
print qq
|</TR></TABLE></TD>
</TR></TABLE>
<DIV id="peptideDIV" style="display:none;margin-left:20px">
<BR style="font-size:5px">
|;
if (scalar @pep_beg) { # not for virtual proteins
	my %usedPeptide;
	my $onclickStrg=(scalar @analysesID == 1 && $validStatus{$analysesID[0]} < 2 && $projectAccess=~/bioinfo|mass|manag|super/)? " onclick=\"startValidation($analysesID[0])\"" : '';
	print qq
|<TABLE border=0 cellspacing=0>
<TR bgcolor="$color2">
<TH class="rbBorder" width=25><I>#</I></TH>
|;
print qq
|<TH class="rbBorder" width=20><ACRONYM onmouseover=\"popup('<B>Peptide validation method:</B><BR>&nbsp;<B>?</B>: Unknown<BR>&nbsp;<B>A</B>: Automated<BR>&nbsp;<B>M</B>: Manual')\" onmouseout=\"popout()\"><IMG src="$promsPath{images}/good.gif"$onclickStrg></ACRONYM></TH>
<TH class="rbBorder">&nbsp;Start&nbsp;</TH>
<TH class="rbBorder">&nbsp;End&nbsp;</TH>
<TH class="rbBorder" nowrap>&nbsp;Spec.<FONT style="font-size:9px;"> %</FONT>&nbsp;</TH>
|;
	if ($showMSdata) {
		print qq
|<TH class="rbBorder">&nbsp;Query&nbsp;</TH>
<TH class="rbBorder" nowrap>&nbsp;RT&nbsp;(min)&nbsp;</FONT></TH>
<TH class="rbBorder">&nbsp;<ACRONYM onmouseover="popup('Mouseover peptides <B>rank</B> to display <B>comments</B> regarding validation.')" onmouseout="popout()">Rank</ACRONYM>&nbsp;</TH>
<TH class="rbBorder">&nbsp;Score&nbsp;</TH>
<TH class="rbBorder">&nbsp;Charge&nbsp;</TH>
<TH class="rbBorder">&nbsp;Mr(Obs)&nbsp;</TH>
<TH class="rbBorder">&nbsp;Mr(Exp)&nbsp;</TH>
<TH class="rbBorder">&nbsp;Mr(calc)&nbsp;</TH>
<TH class="rbBorder">&nbsp;Delta&nbsp;</TH>
<TH class="rbBorder">&nbsp;Miss&nbsp;</TH>
|;
	}

	if (scalar @analysesID > 1) {
		print qq
|<TH class="rbBorder" nowrap align=left>&nbsp;Sequence</TH>
<TH class="bBorder">&nbsp;Analysis&nbsp;</TH>
|;
	} else {
		print "<TH class=\"bBorder\" nowrap align=left>&nbsp;Sequence</TH>\n";
	}
	print "</TR>\n";

	$bgColor=$color1;
	my $pepCount=0;
	foreach my $i (sort {$pep_beg[$a] <=> $pep_beg[$b] || $pep_end[$a] <=> $pep_end[$b] || $id_peptide[$a] <=> $id_peptide[$b]} 0..$#id_peptide) {
		next if $usedPeptide{$id_peptide[$i]}; # peptide is already listed due to multiple occurrences
		$pepCount++;
		my $validMethod=(!defined($peptideInfo{$id_peptide[$i]}[4]))? '*' : ($peptideInfo{$id_peptide[$i]}[4]==0)? '-' : ($peptideInfo{$id_peptide[$i]}[4]==1)? 'A' : 'M'; # select: undef,0,1,2
		my ($beg,$end)=($pep_beg[$i])? ($pep_beg[$i],$pep_end[$i]) : ('?','?');
		print qq
|<TR class="list" bgcolor="$bgColor">
<TH class="rBorder">$pepCount</TH>
<TH>$validMethod</TH>
<TH>$beg</TH>
<TH>$end</TH>
<TH>$pep_specificity{$id_peptide[$i]}</TH>
|;
		my $pepComments=($peptideInfo{$id_peptide[$i]}[4])? '' : '<B>Recovered peptide:</B> Retrieved by ion-profile comparison with other Analyses during quantification step).<BR>';
		$pepComments.=($peptideInfo{$id_peptide[$i]}[5])? quotemeta($peptideInfo{$id_peptide[$i]}[5]) : ($pepComments)? '' : 'No comments.';
		my ($charge,$sub,$elutTimeStrg)=(0,'','');
		if ($showMSdata) {
			if ($peptideInfo{$id_peptide[$i]}[1]) {print "<TH>$peptideInfo{$id_peptide[$i]}[1]</TH>\n";} else {print "<TH>-</TH>\n";}  # query (+/- ghost)
			($charge,$elutTimeStrg,$sub)=@{$peptideInfo{$id_peptide[$i]}}[12..14];
			#$elutTimeStrg=$peptideInfo{$id_peptide[$i]}[13];
			my $elutTime;
			if ($elutTimeStrg) {
				#unless ($elutTimeStrg=~/^sp/) {
				#	$elutTime=($elutTimeStrg=~/et(\d+\.\d+);/)? $1 : $elutTimeStrg;
				#	$elutTime=~s/^et//; # remove et tag if any
				#	#$elutTime.=' <SMALL>min</SMALL>';
				#}
				#$elutTime=($elutTimeStrg=~/^sc(\d+);*/)? "<FONT onmouseover=\"popup('Scan Number')\" onmouseout=\"popout()\">$1</FONT>" : $elutTime;
				if ($elutTimeStrg=~/et(\d+\.*\d*)/) {$elutTime=$1;}
				elsif ($elutTimeStrg=~/sc(\d+)/) {$elutTime="<FONT onmouseover=\"popup('Scan Number')\" onmouseout=\"popout()\">$1</FONT>";}
				elsif ($elutTimeStrg=~/sp(\d+)/) {$elutTime="<FONT onmouseover=\"popup('Spectrum Number')\" onmouseout=\"popout()\">$1</FONT>";}
				else {$elutTime=$elutTimeStrg;}
			}
			else {$elutTime='?';}
			print "<TH nowrap>$elutTime</TH>\n"; # elution time
			my $rank=(!$peptideInfo{$id_peptide[$i]}[2])? '-' : $peptideInfo{$id_peptide[$i]}[2];
			print "<TH><ACRONYM onmouseover=\"popup('$pepComments')\" onmouseout=\"popout()\">$rank</ACRONYM></TH>\n"; # rank + comments
			if ($peptideInfo{$id_peptide[$i]}[6]) {printf "\t<TH>%.1f</TH>\n",$peptideInfo{$id_peptide[$i]}[6];} # score
			else {print "\t<TH>-</TH>\n";}
			print "<TH>$charge<SUP>+</SUP></TH>";
			if (!$peptideInfo{$id_peptide[$i]}[2]) { # ghost
				if ($peptideInfo{$id_peptide[$i]}[10]) {printf "\t<TH>%.2f</TH>\n",$peptideInfo{$id_peptide[$i]}[10];} else {print "\t<TH>-</TH>\n";} # obs
				if ($peptideInfo{$id_peptide[$i]}[8]) {printf "\t<TH>%.2f</TH>\n",$peptideInfo{$id_peptide[$i]}[8];} else {print "\t<TH>-</TH>\n";} # exp
				print "\t<TH>-</TH>\n"; # calc
				print "\t<TH>-</TH>\n"; # delta
			}
			else {
				my $mrCalc=($peptideInfo{$id_peptide[$i]}[9])? sprintf "%.2f",$peptideInfo{$id_peptide[$i]}[9] :'-';
				my $mrObs=($peptideInfo{$id_peptide[$i]}[10])? sprintf "%.2f",$peptideInfo{$id_peptide[$i]}[10] :'-';
				my $mrExp=($peptideInfo{$id_peptide[$i]}[8])? sprintf "%.2f",$peptideInfo{$id_peptide[$i]}[8] :'-';
				print "\t<TH>$mrObs</TH>\n"; # obs
				print "\t<TH>$mrExp</TH>\n"; # exp
				print "\t<TH>$mrCalc</TH>\n"; # calc
				if ($peptideInfo{$id_peptide[$i]}[9] && $peptideInfo{$id_peptide[$i]}[11]) {
                    my $deltaPPM=sprintf "%.3f",1000000*($peptideInfo{$id_peptide[$i]}[11]/$peptideInfo{$id_peptide[$i]}[9]);
					print "\t<TH><A href=\"javascript:void(null)\" onmouseover=\"popup('<B>Delta:</B> $peptideInfo{$id_peptide[$i]}[11] Da ($deltaPPM ppm)')\" onmouseout=\"popout()\">";
					printf "%.2f</A></TH>\n",$deltaPPM; # deltaPPM
                }
                else {print "\t<TH>-</TH>\n";}
			}
			print "<TH>$peptideInfo{$id_peptide[$i]}[7]</TH>\n"; # missed cleavage
		}
		else {
			($elutTimeStrg,$sub)=@{$peptideInfo{$id_peptide[$i]}}[7..8]; # $elutTimeStrg needed to check availability of MAXQUANT virtual peptides fragmentation spectra
		}
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
			    substr($peptideInfo{$id_peptide[$i]}[0] , $pos-1 , 1 , $susbtitutions{$pos} );
			}
		}
		$peptideInfo{$id_peptide[$i]}[0]=~s/[^<]([X,B,Z])/<FONT style="font-size:15px">$1<\/FONT>/g; # ambigous residue(s) shown bigger.
		print "<TH nowrap align=left>&nbsp;"; #<FONT style=\"font-family:monospace;font-size:12pt\">
		my ($flkAA_Nter,$flkAA_Cter);
		if ($flanking_AA[$i]) {
			($flkAA_Nter,$flkAA_Cter)=split(//,$flanking_AA[$i]);
			print "$flkAA_Nter.";
		}
		#if ($analysisInfo{$analysisID}[6] eq 'MIS' && -e $pepFile) {#} # && ($projectAccess eq 'bioinfo' || $projectAccess eq 'mass' || $projectAccess =~ /super/)

		if ($pepFile{$pepAna{$id_peptide[$i]}} && -e $pepFile{$pepAna{$id_peptide[$i]}} && ($peptideInfo{$id_peptide[$i]}[4] || ($fileFormat{$pepAna{$id_peptide[$i]}} eq 'MAXQUANT.DIR' && $elutTimeStrg=~/sc\d+/))) { # not if ghost peptide (except some MAXQUANT)
			my $pepInfo;
			if ($currProtStatus) { # validated protein
				$pepInfo="pep_$id_peptide[$i]_$peptideInfo{$id_peptide[$i]}[1]_$peptideInfo{$id_peptide[$i]}[2]"; # pep_pepID_qNum_rank
			}
			else {$pepInfo="seq_$id_peptide[$i]";} # non-validated protein
			if (($fileFormat{$analysesID[0]} !~ /TDM.PEP.XML|TDM.XML:NoSpec/ && $peptideInfo{$id_peptide[$i]}[6]) || $fileFormat{$analysesID[0]} =~ /SKYLINE\.CSV/ || $fileFormat{$analysesID[0]} =~ /OPENSWATH\.TSV/ ||$fileFormat{$analysesID[0]} =~ /SPECTRONAUT\.XLS/ ) { # no score => no MSMS data (MaxQuant)
				print "<A id=\"$id_peptide[$i]\" href=\"javascript:drawSpectrum('$id_peptide[$i]','$pepInfo')\" onmouseover=\"popup('$pepComments<BR>Click to display <B>Fragmentation Spectrum</B>.')\" onmouseout=\"popout()\">"; # unless $fileFormat{$analysesID[0]}=~/TDM.PEP.XML|TDM.XML:NoSpec/;
			}
			else {
				my $class2=($peptideInfo{$id_peptide[$i]}[4])? '' : ' bi';
				print "<A class=\"no_link$class2\" onmouseover=\"popup('$pepComments')\" onmouseout=\"popout()\">";
			}
		}
		else {
			my $class2=($peptideInfo{$id_peptide[$i]}[4])? '' : ' bi';
			print "<A class=\"no_link$class2\" onmouseover=\"popup('$pepComments')\" onmouseout=\"popout()\">";
		}
		print "$peptideInfo{$id_peptide[$i]}[0]</A>";
		print ".$flkAA_Cter" if $flkAA_Cter;
		if ($peptideInfo{$id_peptide[$i]}[3]) { # var mods
			$peptideInfo{$id_peptide[$i]}[3]=~s/([:.])-1/$1?/g; # if unknown position, -1 has to be removed
			# MaxQuant probabilities #
			my $mqProbStrg='';
			if ($maxQuantProb{$id_peptide[$i]}) {
				$mqProbStrg='<B>MaxQuant probabilities:</B>';
				foreach my $modID (sort{lc($varModNames{$a}) cmp lc($varModNames{$b})} keys %{$maxQuantProb{$id_peptide[$i]}}) {
					$mqProbStrg.="<BR>&bull;$varModNames{$modID}:";
					foreach my $refPos (@{$maxQuantProb{$id_peptide[$i]}{$modID}}) {
						my $posStrg=$convertPos2Text{$refPos->[0]} || $refPos->[0];
						$mqProbStrg.="<BR>&nbsp;&nbsp;-$posStrg:".&promsMod::MaxQuantProbIcon($refPos->[1],{text=>'&nbsp;'.($refPos->[1]*100).'%&nbsp;',inPopup=>1});
					}
				}
			}
			if ($projectAccess=~/bioinfo|mass|manag|super/) {
				$mqProbStrg.='<BR>' if $mqProbStrg;
				print "<A id=\"vmod_$id_peptide[$i]\" style=\"font-size:8pt\" href=\"javascript:editVarModification('vmod_$id_peptide[$i]',",length($peptideInfo{$id_peptide[$i]}[0]),")\" onmouseover=\"popup('${mqProbStrg}Click to edit modifications <B>position</B>.')\" onmouseout=\"popout()\">$peptideInfo{$id_peptide[$i]}[3]</A>";
			}
			else {
				print "<FONT style=\"font-size:8pt\"";
				print " onmouseover=\"popup('$mqProbStrg')\" onmouseout=\"popout()\"" if $mqProbStrg;
				print "> $peptideInfo{$id_peptide[$i]}[3]</FONT>"; # var mods
			}
			# PhosphoRS #
			if ($phosphoRsData{$id_peptide[$i]}) {
				&phosphoRS::printIcon($phosphoRsData{$id_peptide[$i]});
			}
		}
		print "&nbsp;</TH>\n"; #</FONT>
		if (scalar @analysesID > 1) {
			print "<TH align=left>&nbsp;$analysisInfo{$pepAna{$id_peptide[$i]}}[0]&nbsp;</TH>\n";
		}
		print "</TR>\n";

		###>Dealing with multiple occurrences<###
		my $firstBeg=1;
		foreach my $beg (sort{$a<=>$b} keys %{$pep_occurrence{$id_peptide[$i]}}) {
			if ($firstBeg) {$firstBeg=0; next;}
			my $colSpan=(scalar @analysesID > 1)? 3 : 2;
			$colSpan+=10 if $showMSdata;
			print qq
|<TR class="list" bgcolor="$bgColor">
<TH class="rBorder">&nbsp</TH>
<TH></TH>
<TH>$beg</TH>
<TH>$pep_occurrence{$id_peptide[$i]}{$beg}</TH>
<TH colspan=$colSpan></TH><TR>
|; # do not repeat Specificity
		}
		$usedPeptide{$id_peptide[$i]}=1;
		$bgColor=($bgColor eq $color1)? $color2 : $color1;
	}
	print "</TABLE>\n<BR>";
}
print "</DIV>\n";

################################
####>Search for Interactors<####
################################
my ($uniprotAC,$disabInteractors) = ($uniprotItem[0])? ($uniprotItem[0],'') : ('-','disabled');
print qq
|<BR style="font-size:5px">
<TABLE width="100%" bgcolor="$color2"><TR>
<TD><INPUT type="button" id="showInteract" class="font11" value="Show" onclick="ajaxGetInteractors('$uniprotAC')" style="width:55px" $disabInteractors/>
	<INPUT type="button" id="hideInteract" class="font11" value="Hide" style="width:55px; display:none" onclick="hideInteractors()"/>
</TD>
<TH bgcolor="$color1" align="left" width="100%" nowrap>&nbsp;Proteins interacting with $alias in selected Analyses&nbsp;</TH>
</TR></TABLE>
<DIV id="interactDiv" style="display:none;margin-left:20px"></DIV>
|;

########################
####>ProtVista view<####
########################
my $isoformStrg=($isoformAC ne $alias && !$disabProtVista)? " ($isoformAC)" : '';
print qq
|<BR style="font-size:5px">
<TABLE width="100%" bgcolor="$color2"><TR>
<TD><INPUT type="button" id="protVistaBUTTON" class="font11" value="Show" onclick="showHideSection('protVista','$isoformAC')" style="width:55px" $disabProtVista/></TD>
<TH bgcolor="$color1" align="left" width="100%" nowrap>&nbsp;General features of $alias$isoformStrg &nbsp;</TH>
</TR></TABLE>
<DIV id="protVistaDIV" style="display:none;margin-left:20px"><BR style="font-size:5px"><DIV id="protVistaTrueDIV"></DIV></DIV>
|;

########################################
####>Cross-databanks external links<####
########################################
my $disabLink=($masterProtID)? '' : 'disabled=true';
print qq
|<BR style="font-size:5px">
<TABLE width="100%" bgcolor="$color2"><TR>
<TD><INPUT type="button" class="font11" id="linkBUTTON" value="Show" onclick="showHideSection('link')" style="width:55px" $disabLink/></TD>
<TH bgcolor="$color1" align="left" width="100%" nowrap>&nbsp;Links to external resources for $alias&nbsp;</TH>
</TR></TABLE>
<DIV id="linkDIV" style="display:none;margin-left:20px">
<BR style="font-size:5px">
<TABLE border=0 cellspacing=0>
<TR bgcolor="$color2">
<TH class="rbBorder" colspan=2>Links</TH>
<TH class="bBorder" width=500>Resource description</TH>
</TR>
|;
$bgColor=$color1;
foreach my $identID (sort{lc($crossRefInfo{$a}[1]) cmp lc($crossRefInfo{$b}[1])} keys %crossRefInfo) {
	print "<TR bgcolor=\"$bgColor\">\n<TH align=left valign=top nowrap>&nbsp;$crossRefInfo{$identID}[1]:</TH>\n<TH align=left>";
	my $rank=1;
	foreach my $identValue (@{$crossReference{$identID}}) {
		print ',' unless $rank % 2;
		if ($rank > 1 && $crossRefInfo{$identID}[0] eq 'AC') {
			print "&nbsp;$identValue";
		}
		else {
			my $identURL=$crossRefInfo{$identID}[2];
			$identURL=~s/XXXX/$identValue/;
			if ($crossRefInfo{$identID}[0] eq 'RefSeq') {
				my $extraURLstrg=($identValue=~/^NP_/)? 'entrez/viewer.fcgi?db=protein&id=' : 'nuccore/';
				$identURL=~s/VVVV/$extraURLstrg/;
			}
			elsif ($crossRefInfo{$identID}[0] eq 'UniGene') {
				my ($org,$val)=split('\.',$identValue);
				$identURL=~s/YYYY/$org/;
				$identURL=~s/ZZZZ/$val/;
			}
			print "&nbsp;<A href=\"$identURL\" target=\"_blank\">$identValue</A>";
		}
		print "<BR>" unless $rank % 2;
		$rank++;
	}
	print "</TH>\n<TD valign=top>&nbsp;$crossRefInfo{$identID}[3]</TD>\n</TR>\n";
	$bgColor=($bgColor eq $color1)? $color2 : $color1;
}
print qq
|</TABLE>
<FONT class="font11"><B><I>No more matches.</I></B></FONT>
</DIV><BR><BR>
<DIV id="groupDIV" class="popup" style="display:none"></DIV>
<DIV id="vmodPeptideDIV" class="popup"></DIV>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup()</SCRIPT>

</BODY>
</HTML>
|;

######################################<<< SUBROUTINES >>>#########################################

########################################
####<Convert VarMods hash to string>####
########################################
sub convertVarMods {
	my ($refVarMods)=@_;
	my $varModText='';
	foreach my $varMod (sort{$a cmp $b} keys %{$refVarMods}) {
		$varModText.='+' unless $varModText;
		$varModText.="<FONT class=\"$projectVarMods{$varMod}[2]\">$projectVarMods{$varMod}[1]</FONT>";
	}
	#return '+<FONT class="Acetyl">A</FONT><FONT class="Carbamidomethyl">C</FONT><FONT class="Methyl">M</FONT><FONT class="Oxidation">O</FONT><FONT class="Phospho">P</FONT>';
	return $varModText;
}

sub getProtSequence {
	my ($identifier,$analysisID)=@_;
	my $protSeq='';

	###>Extracting protein sequence from file<###
	my $fastaFile="$promsPath{valid}/ana_$analysisID/analysis.fasta";
	if (-e $fastaFile) {
		my $match;
		(my $matchIdentifier=$identifier)=~s/\|/\\\|/g; # converting | into \|
		open (FAS,$fastaFile);
		while (<FAS>) {
			if (/^>/) {
				last if $match;
				if (/^>$matchIdentifier\s/) {$match=1; next;}
			}
			next unless $match;
			$protSeq.=$_;
		}
		$protSeq=~s/\s+//g; # in case of \n
		close FAS;
	}
	return $protSeq;
}

sub ajaxUpdateConnectedProteins {
	my $projectID=&promsMod::cleanNumericalParameters(param('projectID')); # only if ajax call

	###<Connecting to database>###
	my $dbh=&promsConfig::dbConnect;

	my %connnectedProteins;
	&getNumberConnectedProteins($dbh,$proteinID,$projectID,\%connnectedProteins);

	my $sthProtInfo=$dbh->prepare("SELECT ALIAS,MW,PROT_DES,ORGANISM FROM PROTEIN WHERE ID_PROTEIN=?");
	my $sthProtAna=$dbh->prepare("SELECT ID_ANALYSIS,VISIBILITY,CONF_LEVEL FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=? ORDER BY VISIBILITY DESC,NUM_PEP DESC,SCORE DESC,CONF_LEVEL DESC");
	foreach my $protID (keys %connnectedProteins) {
		$sthProtInfo->execute($protID);
		@{$connnectedProteins{$protID}}=$sthProtInfo->fetchrow_array;
		$sthProtAna->execute($protID);
		my ($bestAnaVis,$bestAnaConfLevel)=(0,0);
		my @anaList;
		while (my ($anaID,$visibility,$confLevel)=$sthProtAna->fetchrow_array) {
			push @anaList,$anaID;
			unless ($bestAnaVis) {
				$bestAnaVis=$visibility;
				$bestAnaConfLevel=$confLevel;
			}
		}
		push @{$connnectedProteins{$protID}},($bestAnaVis,$bestAnaConfLevel,\@anaList);

	}
	$sthProtInfo->finish;
	$sthProtAna->finish;
	$dbh->disconnect;

	###<HTML>###
	print header(-charset=>'UTF-8');
	warningsToBrowser(1);

	# Ref protein
	my ($refProtAlias,$refMass,$refProtDes,$refSpecies,$refBestAnaVis,$refBestAnaConf,$refRefAnaList)=@{$connnectedProteins{$proteinID}};
	my ($refClass,$refPopup)=&promsMod::getProteinClass($refBestAnaConf,$refBestAnaVis);
	my $chkRefValueStrg=join(':',@{$refRefAnaList}).':'.$proteinID;
	if ($refMass) {$refMass=sprintf "%.1f",$refMass/1000;} else {$refMass='-';}
	print qq
|<FORM name="protList" action="./graphicalView.cgi" method="post" target="ProteinWindow">
<INPUT type="hidden" name="what" value="checkbox">
<INPUT type="hidden" name="id_item" value="$projectID">
<INPUT type="hidden" name="item" value="project">
<INPUT type="hidden" name="largWindow">
<INPUT type="hidden" name="graphView">
<INPUT type="button" class="font11" value="Check / Unckeck All" onclick="check(document.protList.chkProt)"><INPUT type="button" value="Peptide Distribution" onclick="graphicalView2(document.protList.chkProt)"/>
<TABLE border=0 cellspacing=0>
<TR bgcolor=$color2><TH class="rbBorder">Protein</TH><TH class="rbBorder" nowrap>&nbsp;MW<FONT style="font-size:9px;"> kDa</FONT>&nbsp;</TH><TH class="bBorder" width=550>Description - Species</TH></TR>
<TR class="list" bgcolor="$color1" valign="top">
<TD nowrap><INPUT type="checkbox" name="chkProt" value="$chkRefValueStrg"/><FONT class="$refClass" color="#DD0000">$refProtAlias</FONT></TD>
<TD class="$refClass" align="center">$refMass</TD><TD>&nbsp;$refProtDes <FONT class="org">$refSpecies</FONT></TD>
</TR>
|;
	my $bgColor=$color2;
	foreach my $protID (sort{lc($connnectedProteins{$a}[0]) cmp lc($connnectedProteins{$b}[0])} keys %connnectedProteins) {
		next if $protID==$proteinID;
		my ($protAlias,$mass,$protDes,$species,$bestAnaVis,$bestAnaConf,$refAnaList)=@{$connnectedProteins{$protID}};
		print "<TR class=\"list\" bgcolor=\"$bgColor\" valign=\"top\">\n";
		# Protein & analyses
		my $chkValueStrg=join(':',@{$refAnaList}).':'.$protID;
		my ($pClass,$pPopup)=&promsMod::getProteinClass($bestAnaConf,$bestAnaVis);
		print "<TD nowrap><INPUT type=\"checkbox\" name=\"chkProt\" value=\"$chkValueStrg\"/><A class=\"$pClass\" href=\"javascript:sequenceView($refAnaList->[0],$protID,'valid',opener.top.showMSdata)\" onmouseover=\"popup('$pPopup')\" onmouseout=\"popout()\">$protAlias</A></TD>";
		# MW
		if ($mass) {$mass=sprintf "%.1f",$mass/1000;}
		else {$mass='-';}
		print "<TD class=\"$pClass\" align=\"center\">$mass</TD>";
		# Desc Species
		print "<TD>&nbsp;$protDes <FONT class=\"org\">$species</FONT></TD></TR>\n";
		$bgColor=($bgColor eq $color1)? $color2 : $color1;
	}
	print qq
|</TABLE>
</FORM>
|;
}

sub getNumberConnectedProteins { # ghost peptides are also used!!!
	my ($dbh,$validProtID,$projectID,$refConnProt)=@_;
	my $sthProtPep=$dbh->prepare("SELECT ID_PEPTIDE FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_PROTEIN=$validProtID");
	my $sthConnProt=$dbh->prepare("SELECT ID_PROTEIN FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_PEPTIDE=?");
	#my $sthProtPep=$dbh->prepare("SELECT DISTINCT PEP_SEQ FROM PEPTIDE PEP,PEPTIDE_PROTEIN_ATTRIB PPA WHERE PPA.ID_PEPTIDE=PEP.ID_PEPTIDE AND ID_PROTEIN=$validProtID");
	#my $sthConnProt=$dbh->prepare("SELECT DISTINCT PPA.ID_PROTEIN FROM PEPTIDE PEP,PEPTIDE_PROTEIN_ATTRIB PPA, PROTEIN PRO WHERE PPA.ID_PEPTIDE=PEP.ID_PEPTIDE AND PRO.ID_PROTEIN=PPA.ID_PROTEIN AND ID_PROJECT=$projectID AND PEP_SEQ=?");
	$sthProtPep->execute;
	while (my ($pepSeq)=$sthProtPep->fetchrow_array) {
		$sthConnProt->execute($pepSeq);
		while (my ($protID)=$sthConnProt->fetchrow_array) {
			@{$refConnProt->{$protID}}=(); # includes $validProtID!!!
		}
	}
	$sthConnProt->finish;
	$sthProtPep->finish;
}

sub ajaxDisplayVarModPeptides {
	my $projectID=&promsMod::cleanNumericalParameters(param('projectID')); # only if ajax call
	my $varModCode=param('vModCode');
	my $selVarModPos=param('position');
	my $vModPeptides=param('pepStrg');
	my $cumulate=param('cumulate');

	###> In case more info is checked on validation mode, proteinID has to be updated otherwise, no information will be retrieved by AJAX (added on 19/11/13)
	if (param('id_type') eq 'notValid') {
		my $dbh=&promsConfig::dbConnect;
		($proteinID)=$dbh->selectrow_array("SELECT ID_PROTEIN FROM PROTEIN WHERE ID_PROJECT=$projectID AND IDENTIFIER=(SELECT IDENTIFIER FROM PROTEIN_VALIDATION WHERE ID_PROT_VALID=$proteinID)");
		$dbh->disconnect;
	}

	##<HTML CODE>###
	print header(-charset=>'UTF-8');
	warningsToBrowser(1);
	print qq
|<INPUT type="button" class="font11" value="Close" onclick="document.getElementById('vmodPeptideDIV').style.display='none'"/>
|;
	if ($varModCode =~ /,,/) {# multiple modification
		my @varModCodeTable=split(',,',$varModCode);
		my @vModPeptidesTable=split(',,',$vModPeptides);
		for my $i (0..$#varModCodeTable){
			&ajaxDisplayVarModPeptides2($projectID,$varModCodeTable[$i],$selVarModPos,"$vModPeptidesTable[$i]",$cumulate);
		}
	}else{
		&ajaxDisplayVarModPeptides2($projectID,$varModCode,$selVarModPos,$vModPeptides,$cumulate);
	}
}

sub ajaxDisplayVarModPeptides2 {
	my ($projectID,$varModID,$selVarModPos,$vModPeps,$cumulate)=@_;
	my @vModPeptides=split(',',$vModPeps);

	###<Database queries>###
	my $dbh=&promsConfig::dbConnect;
	my ($psiMsName,$interimName,$altNames)=$dbh->selectrow_array("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES FROM MODIFICATION WHERE ID_MODIFICATION=$varModID");
	my $vModName=($psiMsName)?$psiMsName:($interimName)?$interimName:($altNames)?$altNames:'';
	my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
	my $projectAccess=${$userInfo[2]}{$projectID};
	my (%peptideData,$residue,%varModName);

	#my $sthPep=$dbh->prepare("SELECT ABS(PEP_BEG),ABS(PEP_END),PEP_SEQ,0,SCORE,QUERY_NUM,PEP_RANK,DATA,P.COMMENTS,A.NAME
	#			FROM PEPTIDE P,PEPTIDE_PROTEIN_ATTRIB PA, ANALYSIS A
	#			WHERE P.ID_PEPTIDE=? AND PA.ID_PEPTIDE=? AND PA.ID_PROTEIN=$proteinID
	#			AND P.ID_ANALYSIS=A.ID_ANALYSIS AND $selVarModPos >= ABS(PEP_BEG) AND $selVarModPos <= ABS(PEP_END)"); # $selVarModPos constraint in case multiple occurence of peptide

		my $sthPep=$dbh->prepare(qq
|SELECT ABS(PEP_BEG),ABS(PEP_END),
	PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),
	GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',COALESCE(PM.REF_POS_STRING,'') ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),
	SCORE,QUERY_NUM,PEP_RANK,DATA,P.COMMENTS,A.ID_ANALYSIS,A.NAME,A.FILE_FORMAT,A.DATA_FILE,A.VALID_STATUS
	FROM PEPTIDE P
	INNER JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
	INNER JOIN PEPTIDE_PROTEIN_ATTRIB PA ON P.ID_PEPTIDE=PA.ID_PEPTIDE AND PA.ID_PROTEIN=$proteinID
	INNER JOIN ANALYSIS A ON P.ID_ANALYSIS=A.ID_ANALYSIS
	WHERE P.ID_PEPTIDE=? AND $selVarModPos >= ABS(PEP_BEG) AND $selVarModPos <= ABS(PEP_END)
	GROUP BY P.ID_PEPTIDE
|);
# INNER JOIN on PEPTIDE_MODIFICATION because at least 1 PTM
# $selVarModPos constraint in case multiple occurence of peptide

	foreach my $pepID (@vModPeptides) {
		$sthPep->execute($pepID);
		@{$peptideData{$pepID}}=$sthPep->fetchrow_array;
		$peptideData{$pepID}[5]=0 unless $peptideData{$pepID}[5]; # score
		## Add VMOD in $peptideData{$pepID}[3]
		#my ($anaID)=$dbh->selectrow_array("SELECT ID_ANALYSIS FROM PEPTIDE WHERE ID_PEPTIDE=$pepID");
		#$peptideData{$pepID}[3]=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$anaID,$peptideData{$pepID}[2]);
		#$peptideData{$pepID}[3]=~s/^ *\+ *//; # vmod
		#$peptideData{$pepID}[3]=~s/([:.])-1/$1?/g; # if unknown position, -1 has to be removed
		unless ($residue) {
			my @seqArray=split('',$peptideData{$pepID}[2]); # seq
			$residue=$seqArray[$selVarModPos-$peptideData{$pepID}[0]];
		}
	}
	$sthPep->finish;

	print qq
|&nbsp;&nbsp;<FONT class="title3">$vModName on $residue<SUB>$selVarModPos</SUB></FONT>
<TABLE border=0 cellspacing=0>
<TR bgcolor="$color2">|;
	print qq
|<TH class="rbBorder">&nbsp;Start&nbsp;</TH><TH class="rbBorder">&nbsp;End&nbsp;</TH><TH class="rbBorder">Sequence</TH><TH class="rbBorder">PTMs</TH>
|;
	if ($cumulate) {
		print "<TH class=\"rbBorder\">&nbsp;Score&nbsp;</TH><TH class=\"bBorder\">&nbsp;Analysis&nbsp;</TH>";
	} else {
		print "<TH class=\"bBorder\">&nbsp;Score&nbsp;</TH>";
	}
	print "</TR>\n";

	my $bgColor=$color1;
	my $minStart;
	foreach my $pepID (sort{$peptideData{$a}[0]<=>$peptideData{$b}[0] || $peptideData{$a}[1]<=>$peptideData{$b}[1] || $peptideData{$b}[5] <=> $peptideData{$a}[5] || lc($peptideData{$a}[11]) cmp lc($peptideData{$b}[11]) } keys %peptideData) {
		my ($start,$end,$seq,$vModCode,$refModInfo,$score,$qNum,$rank,$data,$pepComments,$anaID,$anaName,$fileFormat,$dataFile,$validStatus)=@{$peptideData{$pepID}};
		my $vMod=&promsMod::decodeVarMod($dbh,$seq,$vModCode,\%varModName); # if $vModCode;
		$qNum=0 unless $qNum;
		$minStart=$start if !defined $minStart;
		my $spaceStrg='&nbsp;' x ($start-$minStart);
		my $vModSeqIdx=$selVarModPos-$start;
		my @seqArray=split('',$seq);
		$seqArray[$vModSeqIdx]='<FONT color="#DD0000">'.$seqArray[$vModSeqIdx].'</FONT>';
		my $dispSeq=join('',@seqArray);
		my $dispVarMod=$vMod;

		my ($PRS) = $data =~ /PRS=([^##]+)/ if $data;
		$pepComments=($pepComments)? quotemeta($pepComments) : 'No comments.';
		$dataFile='msms.txt' if $fileFormat eq 'MAXQUANT.DIR';
		my $pepFile=($validStatus==2)? "$promsPath{peptide}/proj_$projectID/ana_$anaID/$dataFile" : "$promsPath{valid}/ana_$anaID/$dataFile";
		print qq
|<TR bgcolor="$bgColor">
<TH align=right>$start&nbsp;</TH><TH align=right>$end&nbsp;</TH><TH class="seq" nowrap align=left>$spaceStrg|;
		if ($score && -e $pepFile) {
			my $pepInfo='pep_'.$pepID.'_'.$qNum.'_'.$rank;
			print "<A id=\"pep2_$pepID\" href=\"javascript:drawSpectrum('pep2_$pepID','$pepInfo')\" onmouseover=\"popup('$pepComments<BR>Click to display <B>Fragmentation Spectrum</B>.')\" onmouseout=\"popout()\">$dispSeq</A>";
		}
		else {print "<A href=\"javascript:void(null)\" onmouseover=\"popup('$pepComments')\" onmouseout=\"popout()\">$dispSeq</A>";}
		print "&nbsp;</TH><TH nowrap align=left>";
		if ($projectAccess=~/bioinfo|mass|manag|super/) {
			print "<A id=\"vmod2_$pepID\" style=\"font-size:8pt\" href=\"javascript:editVarModification('vmod2_$pepID',",length($seq),")\" onmouseover=\"popup('Click to edit modifications <B>position</B>.')\" onmouseout=\"popout()\">$dispVarMod</A>";
		}
		else {
			print "<FONT style=\"font-size:8pt\">&nbsp;$dispVarMod</FONT>"; # var mods
		}
		#<PhosphoRS
		if ($PRS) {
			&phosphoRS::printIcon($PRS);
		}
		#<MaxQuant
		if ($fileFormat eq 'MAXQUANT.DIR') {
			my @maxQuantProb;
			if ($refModInfo && $refModInfo=~/(^|&)$varModID:[^&#]*##PRB_MQ=([^&#]+)/) {
				my $modProbStrg=$2;
				while ($modProbStrg =~ /([^,:]+):([^,]+)/g) {
					push @maxQuantProb,[$1,$2];
				}
			}
			elsif ($score) {
				my ($b,$modPosStrg)=$vModCode=~/(^|&)$varModID:([^&]+)/;
				foreach my $pos (split(/\./,$modPosStrg)) {
					push @maxQuantProb,[$pos,1]; # set to 100%
				}
			}
			foreach my $refPos (@maxQuantProb) {
				print &promsMod::MaxQuantProbIcon($refPos->[1],{popupText=>"MaxQuant probability for pos. <B>$refPos->[0]</B>: <B>".($refPos->[1]*100).'%</B>'});
			}
		}
		if ($score) {$score=sprintf "%.2f",$score;}
		else {$score='-';} # PMF or ghost
		print "&nbsp;</TH><TH>$score</TH>\n";

		if ($cumulate) {
			print "<TH align=center>&nbsp;$anaName&nbsp;</TH>"
		}

		print "</TR>\n";

		$bgColor=($bgColor eq $color1)? $color2 : $color1;
	}
	print "</TABLE>";

	$dbh->disconnect;
}

sub ajaxGetQuantificationList {
#print header(-'content-encoding'=>'no',-charset=>'UTF-8');
	require promsQuantif;
	my %proteinQuantifFamilies=&promsQuantif::getProteinQuantifFamilies;
	my $dbh=&promsConfig::dbConnect;
	my ($item,$itemIdStrg)=split(':',param('designBranch'));
	my @itemID=&promsMod::cleanNumericalParameters(split(',',$itemIdStrg));

	my (%quantifList,%anaQuantifs,%modifications);
	my (%quantifValues,%quantifInfo,%quantifHierarchy,%quantifPos);
	my @sthQuantifs;
	if ($item eq 'design') { # Ratio/MQ quantifs
		$sthQuantifs[0]=$dbh->prepare("SELECT DISTINCT QM.CODE,Q.ID_QUANTIFICATION,PQ.TARGET_POS,Q.ID_MODIFICATION,PQ.ID_QUANTIF_PARAMETER
										FROM QUANTIFICATION Q,PROTEIN_QUANTIFICATION PQ,QUANTIFICATION_METHOD QM
										WHERE Q.ID_QUANTIFICATION=PQ.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD
										AND Q.ID_DESIGN=$itemID[0] AND PQ.ID_PROTEIN=$proteinID AND TARGET_POS IS NOT NULL");
	}
	else { # internal-analysis
		my $analysisStrg;
		if ($itemID[0]==0) {
			my $sthAna=$dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=$proteinID");
			my @analysisList;
			$sthAna->execute;
			while (my ($anaID)=$sthAna->fetchrow_array) {
				push @analysisList,$anaID;
			}
			$analysisStrg=join(',',@analysisList);
		}
		else {$analysisStrg=join(',',@itemID);}

		#<Ratio-quantifs
		$sthQuantifs[0]=$dbh->prepare("SELECT DISTINCT QM.CODE,Q.ID_QUANTIFICATION,PQ.TARGET_POS,Q.ID_MODIFICATION
										FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q,PROTEIN_QUANTIFICATION PQ,QUANTIFICATION_METHOD QM
										WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION=PQ.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD
										AND Q.ID_DESIGN IS NULL AND PQ.ID_PROTEIN=$proteinID AND TARGET_POS IS NOT NULL AND AQ.ID_ANALYSIS IN ($analysisStrg)");
		#<Non-ratio quantifs (SIN,EMPAI)
		$sthQuantifs[1]=$dbh->prepare("SELECT DISTINCT QM.CODE,Q.ID_QUANTIFICATION,1,Q.ID_MODIFICATION
										FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q,PROTEIN_QUANTIFICATION PQ,QUANTIFICATION_METHOD QM
										WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION=PQ.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD
										AND Q.ID_DESIGN IS NULL AND PQ.ID_PROTEIN=$proteinID AND QM.CODE IN ('SIN','EMPAI') AND AQ.ID_ANALYSIS IN ($analysisStrg)");
	}

	###<Fetching quantifs>###
	my $ratioParamID;
	foreach my $sthQuantif (@sthQuantifs) {
		$sthQuantif->execute;
		while (my ($code,$quantifID,$targetPos,$modifID,$quantifParamID)=$sthQuantif->fetchrow_array) { # anaID only for internal quantifs
			if ($code eq 'PROT_RATIO_PEP') {
				unless ($ratioParamID) {
					($ratioParamID)=$dbh->selectrow_array("SELECT QP.ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER QP,QUANTIFICATION_METHOD QM WHERE QP.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND QM.CODE='PROT_RATIO_PEP' AND QP.CODE='RATIO'");
				}
				next if $quantifParamID != $ratioParamID; # this prevents non-ratio target_pos to be included (state mean, ...)
			}
			my $quantifFamily=($code=~/_RATIO_|TNPQ/)? 'RATIO' : $code;
			push @{$quantifList{$quantifFamily}},$quantifID."_".$targetPos;
			@{$modifications{$modifID}}=() if $modifID;
		}
		$sthQuantif->finish;
	}

	###<Fetching ratio-quantification data>###
	if ($quantifList{'RATIO'}) { # RATIO
		my %parameters=(QUANTIF_FAMILY=>'RATIO',VIEW=>'',NUM_PEP_CODE=>'NUM_PEP_USED',QUANTIF_LIST=>$quantifList{'RATIO'}); # VIEW & NUM_PEP_CODE only for RATIO
		&promsQuantif::fetchQuantificationData($dbh,\%parameters,\%quantifInfo,\%quantifValues,undef,{$proteinID=>1});
		foreach my $quantif (@{$quantifList{'RATIO'}}) {
			my ($quantifID,$ratioPos)=split('_',$quantif);
			my $ratioIdx=$ratioPos-1;
			my ($testCondID,$refCondID)=split(/\//,$quantifInfo{$quantifID}[1]->{'RATIOS'}[$ratioIdx]); # id flags removed by &promsQuantif::extractQuantificationParameters
			my $superRatioTag=($testCondID=~/%/)? '°' : '';
			$testCondID=~s/%\d+//;
			$refCondID=~s/%\d+//;
			my $ratioName=$quantifInfo{$quantifID}[2]->{$testCondID}{'NAME'}."$superRatioTag/".$quantifInfo{$quantifID}[2]->{$refCondID}{'NAME'}.$superRatioTag;
			my $quantifPath;
			if ($quantifInfo{$quantifID}[3]->[-2]{'ITEM'} eq 'DESIGN') {
				$quantifPath=$quantifInfo{$quantifID}[0];
				$quantifPos{$quantif}=0;
			}
			else { # ANALYSIS
				$quantifPos{$quantif}='1';
				foreach my $refParent (@{$quantifInfo{$quantifID}[3]}) {
					next if $refParent->{'ITEM'} eq 'PROJECT';
					$quantifPath.=' > ' if $quantifPath;
					$quantifPath.=$refParent->{'NAME'};
					$quantifPos{$quantif}.=sprintf "%03d",$refParent->{'POS'};
				}
			}
			@{$quantifHierarchy{$quantif}}=($quantifPath,$ratioName,$ratioPos);
		}
	}
	if ($quantifList{'MQ'}) {
		my %parameters=(QUANTIF_FAMILY=>'MQ',VIEW=>'',MEASURE=>['MQ_INT','MQ_IBAQ','MQ_LFQ'],NUM_PEP_CODE=>'PEPTIDES',QUANTIF_LIST=>$quantifList{'MQ'});
		&promsQuantif::fetchQuantificationData($dbh,\%parameters,\%quantifInfo,\%quantifValues,undef,{$proteinID=>1});
		foreach my $quantif (@{$quantifList{'MQ'}}) {
			my ($quantifID,$ratioPos)=split('_',$quantif);
			my $quantifPath;
			$quantifPos{$quantif}='1';
			foreach my $refParent (@{$quantifInfo{$quantifID}[3]}) {
				next if $refParent->{'ITEM'} eq 'PROJECT';
				last if $refParent->{'ITEM'} eq 'QUANTIFICATION'; # skip quantif Name
				$quantifPath.=' > ' if $quantifPath;
				$quantifPath.=$refParent->{'NAME'};
				$quantifPos{$quantif}.=sprintf "%03d",$refParent->{'POS'} if $refParent->{'ITEM'} eq 'EXPERIMENT';
			}
			$quantifPos{$quantif}.=sprintf "%05d",$quantifID; # samp,ana,quantif generated together
			@{$quantifHierarchy{$quantif}}=($quantifPath,$quantifInfo{$quantifID}[0]);
		}
	}
	foreach my $quantifFamily ('EMPAI','SIN') { # all measures are fetched
		if ($quantifList{$quantifFamily}) {
			my %parameters=(QUANTIF_FAMILY=>$quantifFamily,VIEW=>'',NUM_PEP_CODE=>'',QUANTIF_LIST=>$quantifList{$quantifFamily});
			&promsQuantif::fetchQuantificationData($dbh,\%parameters,\%quantifInfo,\%quantifValues,undef,{$proteinID=>1});
			foreach my $quantif (@{$quantifList{$quantifFamily}}) {
				my ($quantifID,$ratioPos)=split('_',$quantif);
				my $quantifPath;
				$quantifPos{$quantif}='1';
				foreach my $refParent (@{$quantifInfo{$quantifID}[3]}) {
					next if $refParent->{'ITEM'} eq 'PROJECT';
					last if $refParent->{'ITEM'} eq 'QUANTIFICATION'; # skip quantif Name
					$quantifPath.=' > ' if $quantifPath;
					$quantifPath.=$refParent->{'NAME'};
					$quantifPos{$quantif}.=sprintf "%03d",$refParent->{'POS'}; # 1 -> 001 (eg EXP1+SAMP1+ANA1=1001001001001 < EXP1+SAMP2+ANA1=1001001002001)
				}
				@{$quantifHierarchy{$quantif}}=($quantifPath,$quantifInfo{$quantifID}[0]);
			}
		}
	}

	###<Modification quantifications>###
	if (scalar keys %modifications) {
		my $sthModif=$dbh->prepare("SELECT DISPLAY_CODE,DISPLAY_COLOR FROM MODIFICATION WHERE ID_MODIFICATION=?");
		foreach my $modifID (keys %modifications) {
			$sthModif->execute($modifID);
			@{$modifications{$modifID}}=$sthModif->fetchrow_array;
		}
		$sthModif->finish;
	}

	$dbh->disconnect;

	###<Display results>###
	print header(-charset=>'UTF-8');
	warningsToBrowser(1);

	foreach my $quantifFamily (sort keys %quantifList) { # RATIO | MQ | EMPAI | SIN
		print "<FONT class=\"title3\">&bull;$proteinQuantifFamilies{NAME}{$quantifFamily}:<FONT>\n";
		if ($quantifFamily eq 'RATIO') {
			print qq
|<TABLE border=0 cellspacing=0>
<TR bgcolor="$color2">
<TH class="rbBorder" align="left">&nbsp;Quantification&nbsp;</TH>
<TH class="rbBorder" nowrap>&nbsp;Modification&nbsp;</TH>
<TH class="rbBorder">&nbsp;Ratio&nbsp;</TH>
<TH class="rbBorder">&nbsp;p-value&nbsp;</TH>
<TH class="bBorder">&nbsp;Peptides&nbsp;</TH>
</TR>
|;
			$bgColor=$color1;
			my $prevQuantifName='';
			foreach my $quantif (sort{lc($quantifHierarchy{$a}[0]) cmp lc($quantifHierarchy{$b}[0]) || $quantifHierarchy{$a}[2]<=>$quantifHierarchy{$b}[2] || &promsMod::sortSmart($a,$b)} @{$quantifList{$quantifFamily}}) {
				my $quantifNameStrg=($prevQuantifName eq $quantifHierarchy{$quantif}[0])? "<FONT style=\"visibility:hidden\">$quantifHierarchy{$quantif}[0]&nbsp;:&nbsp;</FONT>" : "$quantifHierarchy{$quantif}[0]&nbsp;:&nbsp;";
				$prevQuantifName=$quantifHierarchy{$quantif}[0];
				next if $quantifHierarchy{$quantif}[1] eq '/'; # Some quantitation went wrong and there are no quantitation ratios computed !
				print qq
|<TR bgcolor="$bgColor" class="list">
<TH nowrap align="left" valign=top>&nbsp;$quantifNameStrg$quantifHierarchy{$quantif}[1]&nbsp;</TH><TH align="left">|;
				my ($quantifID,$targetPos)=split('_',$quantif);
				my $modID=$quantifInfo{$quantifID}[4];
				my $modStrg;
				if ($modID) {
					my ($dispCode, $dispColor)=@{$modifications{$modID}};
					$modStrg="&nbsp;<FONT color=\"$dispColor\">$dispCode</FONT>";
					foreach my $protMod (keys %{$quantifValues{$quantif}}) {
						my ($protID,$modif)=split('-',$protMod);
						$modif='?' unless $modif;
						$modif=~s/^\[$modID\]//;
						&promsQuantif::decodeModificationSite($modif);
						print "$modStrg-$modif&nbsp;<br>";
					}
				}
				else {print "&nbsp;-&nbsp;"}
				print "</TH>\n<TH>";
				# ratio
				foreach my $protMod (keys %{$quantifValues{$quantif}}) {
					my $ratio=$quantifValues{$quantif}{$protMod}{'RATIO'};
					$ratio=($ratio == 1000)? "&infin;" : ($ratio == 0.001)? "1/&infin;" : ($ratio <1)? sprintf('1/%.2f',1/$ratio) : sprintf('%.2f',$ratio);
					print "&nbsp;$ratio&nbsp;<br>";
				}
				print "</TH>\n<TH>";
				# p-value
				foreach my $protMod (keys %{$quantifValues{$quantif}}) {
					my $pValue=(!defined($quantifValues{$quantif}{$protMod}{'P_VALUE'}))? '-' : ($quantifValues{$quantif}{$protMod}{'P_VALUE'}==1 || $quantifValues{$quantif}{$protMod}{'P_VALUE'}==0)? $quantifValues{$quantif}{$protMod}{'P_VALUE'} : sprintf('%.2e',$quantifValues{$quantif}{$protMod}{'P_VALUE'});
					print "&nbsp;$pValue&nbsp;<br>";
				}
				print "</TH>\n<TH>";
				# peptides
				foreach my $protMod (keys %{$quantifValues{$quantif}}) {
					print "&nbsp;$quantifValues{$quantif}{$protMod}{NUM_PEP_USED}&nbsp;<BR>\n";
				}
				print "</TH></TR>\n";
				$bgColor=($bgColor eq $color1)? $color2 : $color1;
			}
			print "</TABLE>\n<BR>";
		}

		else { # MQ,EMPAI,SIN
			print qq
|<TABLE border=0 cellspacing=0>
<TR bgcolor="$color2">
<TH class="rbBorder" align="left">&nbsp;Quantification&nbsp;</TH>
|;
			if ($quantifFamily eq 'MQ') {print "<TH class=\"rbBorder\">&nbsp;Intensity&nbsp;</TH><TH class=\"rbBorder\">&nbsp;iBAQ&nbsp;</TH><TH class=\"rbBorder\">&nbsp;LFQ&nbsp;</TH>";}
			elsif ($quantifFamily eq 'EMPAI') {print "<TH class=\"rbBorder\">&nbsp;emPAI&nbsp;</TH><TH class=\"rbBorder\">&nbsp;emPAI (Mol %)&nbsp;</TH><TH class=\"rbBorder\">&nbsp;emPAI (Mr %)&nbsp;</TH>";}
			elsif ($quantifFamily eq 'SIN') {print "<TH class=\"rbBorder\">&nbsp;SI<SUB>N</SUB>&nbsp;</TH>";}
			print "<TH class=\"bBorder\">&nbsp;Peptides&nbsp;</TH></TR>\n";
			
			$bgColor=$color1;
			my $prevParentsName='';
			foreach my $quantif (sort{$quantifPos{$a}<=>$quantifPos{$b} || lc($quantifHierarchy{$a}[0]) cmp lc($quantifHierarchy{$b}[0]) || lc($quantifHierarchy{$a}[1]) cmp lc($quantifHierarchy{$b}[1])} @{$quantifList{$quantifFamily}}) {
				my $quantifNameStrg=($prevParentsName eq $quantifHierarchy{$quantif}[0])? "<FONT style=\"visibility:hidden\">$quantifHierarchy{$quantif}[0]&nbsp;>&nbsp;</FONT>" : "$quantifHierarchy{$quantif}[0]&nbsp;>&nbsp;";
				$prevParentsName=$quantifHierarchy{$quantif}[0];
				my $quantifValueStrg='';
				if ($quantifFamily eq 'MQ') {
					my $intensity=$quantifValues{$quantif}{$proteinID}{MQ_INT} || '-';
					my $iBaq=$quantifValues{$quantif}{$proteinID}{MQ_IBAQ} || '-';
					my $lfq=$quantifValues{$quantif}{$proteinID}{MQ_LFQ} || '-';
					$quantifValueStrg="<TH align=\"right\">&nbsp;$intensity&nbsp;</TH><TH align=\"right\">&nbsp;$iBaq&nbsp;</TH><TH align=\"right\">&nbsp;$lfq&nbsp;</TH>";
				}
				elsif ($quantifFamily eq 'EMPAI') {
					my $empai=(defined($quantifValues{$quantif}{$proteinID}{EMPAI}))? $quantifValues{$quantif}{$proteinID}{EMPAI} : '-';
					my $empaiMol=(defined($quantifValues{$quantif}{$proteinID}{EMPAI_MOL}))? $quantifValues{$quantif}{$proteinID}{EMPAI_MOL} : '-';
					my $empaiMr=(defined($quantifValues{$quantif}{$proteinID}{EMPAI_MR}))? $quantifValues{$quantif}{$proteinID}{EMPAI_MR} : '-';
					$quantifValueStrg="<TH align=\"right\">&nbsp;$empai&nbsp;</TH><TH align=\"right\">&nbsp;$empaiMol&nbsp;</TH><TH align=\"right\">&nbsp;$empaiMr&nbsp;</TH>";
				}
				elsif ($quantifFamily eq 'SIN') {
					my $sin=$quantifValues{$quantif}{$proteinID}{SIN_SIN} || '-';
					$quantifValueStrg="<TH align=\"right\">&nbsp;$sin&nbsp;</TH>";
				}
				#my $quantifValue=sprintf '%.2e',$quantifValues{$quantifID};
				print qq
|<TR bgcolor="$bgColor" class="list">
<TH nowrap align="left" valign=top>&nbsp;$quantifNameStrg$quantifHierarchy{$quantif}[1]&nbsp;</TH>$quantifValueStrg<TH>&nbsp;$quantifValues{$quantif}{$proteinID}{PEPTIDES}&nbsp;</TH></TR>
|;
				$bgColor=($bgColor eq $color1)? $color2 : $color1;
			}
			print "</TABLE>\n<BR>";
		}
	}
	print "<FONT class=\"title3\" color=\"#DD0000\">No quantifications found!</FONT>\n" unless scalar keys %quantifList;
	exit;
}

####>Revision history<####
# 3.2.9 Improved support for display of MQ, EMPAI and SIN quantifications (PP 02/11/18)
# 3.2.8 Handles project status=-1 [no auto-end validation] (PP 07/06/18)
# 3.2.7 Improved check for spectrum availability & added MaxQuant position probability in &ajaxDisplayVarModPeptides2 (PP 22/03/18)
# 3.2.6 [Fix] Bug in AJAX &ajaxDisplayVarModPeptides2 to point to good peptide position when multiple occurences (PP 07/03/18)
# 3.2.5 Display peptide fragmentation spectrum for spectronaut quantification (MLP 22/01/18)
# 3.2.4 Fix minor bug to remove [modifID] in modification site quantification list (PP 11/12/17)
# 3.2.3 Minor modif to avoid error when clicking on a hidden protein in Match group column (GA 03/10/17)
# 3.2.2 Minor modif to display fragmentation spectrum for openswath data (MLP 01/09/17)
# 3.2.1 Minor modif to display fragmentation spectrum for skyline data (MLP 23/06/17)
# 3.2.0 ProtVista features & other minor changes (PP 03/05/17)
# 3.1.4 Minor modificiation to prevent Protein ratio display (GA 27/03/17)
# 3.1.3 Displays MaxQuant PTM position probabilities & ghost peptides count when available (PP 15/02/17)
# 3.1.2 Change to match use of mqpar.xml file in ANALYSIS.DATA_FILE field for MaxQuant & code cleaning (PP 08/02/17)
# 3.1.1 Add a next in a loop to prevent query $sthSpec to fail when there is no valid peptides for a validated protein<BR>This issue is linked with MassChroQ retrieved ghost-peptides (GA 19/01/17)
# 3.1.0 Minor modif to desable fragmentation spectrum vizualisation for TMD.PEP.XML and TDM.XML:NoSpec (MLP 16/12/16)
# 3.0.9 Bug fix in peptide valid status query (PP 14/12/16)
# 3.0.8 Bug fix in protein sequence display when peptide positions are missing (PP 08/12/16)
# 3.0.7 Compatibility with MaxQuant & partial code optimisation of validated peptides data retrieval (PP 29/11/16)
# 3.0.6 Minor modifications for mr obs, mr calc and delta (MLP 20/06/16)
# 3.0.5 Minor CSS adjustment (PP 05/04/16)
# 3.0.4 Added &lt;TITLE&gt; tag (PP 14/03/16)
# 3.0.3 Displays list of candidate isoforms (PP 05/02/16)
# 3.0.2 Extends protein quantification data to SIN and emPAI (PP 10/12/15)
# 3.0.1 Bug fix in match group call & update in Analyses listing (PP 30/10/15)
# 3.0.0 Fix minor bug for MaxQuant case (GA 30/10/15)
# 2.9.9 add protein quantification list view, display only first and all analysis (SL 03/09/15)
# 2.9.8 Interactor search disabled if no UniProt AC mapped (PP 19/03/15)
# 2.9.7 Fixes print of undef query/rank for ghost peptides in peptide table (PP 12/01/15)
# 2.9.6 Fix for ghost peptides in &ajaxDisplayVarModPeptides2 (PP 23/12/14)
# 2.9.5 Minor display chnage (PP 04/12/14)
# 2.9.4 Minor display changes (PP 29/10/14)
# 2.9.3 Add interactions list (SL 06/10/14)
# 2.9.2 Added background color code for confidence in visibility assignment (PP 18/09/14)
# 2.9.1 Updated because of parameter change in graphicalView.cgi (PP 25/08/14)
# 2.9.0 Changed analyses listing text when cumulating multiple analyses (PP 10/07/14)
# 2.8.9 Add popup for delta (GA 02/07/14)
# 2.8.8 Minor bug fix on undef $flanking_AA[$i] (PP 06/06/14)
# 2.8.7 Added AJAX call to display selected match group in popup div (PP 22/04/14)
# 2.8.6 Better check on protein sequence (PP 18/03/14)
# 2.8.5 Minor display bug fix (PP 06/04/14)
# 2.8.4 Two minor change in ajaxDisplayVarModPeptides and for Protein N-term modifications (GA 19/11/13)
# 2.8.3 Restrict ambiguous PTMs + Bug correction for Protein N-term modifications (GA 17/09/13)
# 2.8.2 No link to validation mode for virtual proteins (PP 06/09/13)
# 2.8.1 Minor modification for ambigous position (GA 05/09/13)
# 2.8.0 Minor change in CSS (PP 03/09/13)
# 2.7.9 Fix multiple abnormal behaviors with cummulated analyses or virtual protein (PP 27/08/13)
# 2.7.8 Minor warning correction for PTMs when it occurs on N-term/C-term positions (GA 22/08/13)
# 2.7.7 Minor changes in ajaxDisplayVarModPeptides and PTMs information (GA 05/07/13)
# 2.7.6 Modification for PTM viewing<BR>Minor bug correction in El. time column for Peptide List viewing -> print Scan Number (GA 14/06/13)
# 2.7.5 Remove VAR_MOD from script (GA 23/05/13)
# 2.7.4 Fix analysis checkbox selection bug (FY 17/05/13) +<BR> Displaying analysis name in var mod view for heterogeneous var mods at same position
# 2.7.3 Modifications for &getVariableModifications, now in promsMod<BR>Add query to get ID_ANALYSIS if param('id_ana') is empty (GA 14/05/13)
# 2.7.2 Add protein cumulating option (FY 26/04/13)
# 2.7.1 Fix coverage bug due to extra SCORE field in Ana/Prot query (PP 24/04/13)
# 2.7.0 Minor modification following &convertVarModString update in promsMod (GA 18/04/13)
# 2.6.9 Add substitution visualization (GA 21/03/13)
# 2.6.8 Handles GI accession as dbXlink (PP 12/03/13)
# 2.6.7 Handles PROT_SEQ='+' & new search file path for fully validated analysis (PP 04/03/13)
# 2.6.6 No frames & handles virtual proteins (PP 05/02/13)
# 2.6.5 Keep popup with PTM commentary information active even if problem with spectrum link in "Peptide List" (GA 29/01/13)
# 2.6.4 Protein isoform management through MASTER_PROTEIN (PP 24/01/13)
# 2.6.3 Manages peptide/protein sequence incompatibility (PP 24/01/13)
# 2.6.2 Update for PARAGON searches (GA 11/12/12)
# 2.6.1 Minor modification: for user defined proteins, $geneInfo{'OS'} could not be defined... (GA 28/11/12)
# 2.6.0 External links use new identifier mapping (PP 15/11/12)
# 2.5.8 Project status management (PP 20/09/12)
# 2.5.7 Update of PICR report parsing (PP 29/06/12)
# 2.5.6 Bug correction (PP 13/06/12)
# 2.5.5 Aesthetic modification for var-mod printing on sequence and pop-out in case of ambiguity and/or multiple modifications (GA 13/06/12)
# 2.5.4 Section organization & varMod popup info (PP 23/05/12)
# 2.5.3 Change PEPTIDE.DATA string separator ',' to '##' (FY 26/04/12)
# 2.5.2 Add PhosphoRS notifications on peptide list (FY 03/04/12)
# 2.5.1 Fixes bug when highlight overlap >=2 at beg of seq line (PP 27/03/12)
# 2.5.0 Handles positioning of PTMs with unknown position (PP 26/03/12)
# 2.4.9 Bug fix for peptide overlap starting at pos 1 (PP 07/03/12)
# 2.4.8 New varMod display on sequence (PP 02/03/12)
# 2.4.7 Manager has full access validation data (PP 01/10/11)
# 2.4.6 handles unknown (not listed in promsConfig.pm) varMods (PP 01/08/11)
# 2.4.5 Handling ghost peptides: valid_status=0 & bug from SQ search (PP 13/07/11)
# 2.4.4 SQ search & full datafile conservation managements (PP 01/06/11)
# 2.4.3 Correcting bug in peptide connection check when protein is not validated (PP 13/05/11)
# 2.4.2 Adding PTM view (FY 27/04/11)
# 2.4.1 Minor bug correction sorting PICR accessions (PP 15/04/11)
# 2.4.0 adapting to new protein list management options (PP 22/02/2011)
# 2.3.9 Minor update to print well the elution time (GA 07/01/11)<BR>See 2.6.5 modification in storeAnalyses.cgi
