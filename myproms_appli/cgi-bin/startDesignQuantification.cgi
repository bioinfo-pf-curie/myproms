#!/usr/local/bin/perl -w

################################################################################
# startDesignQuantification.cgi      1.5.1                                     #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
################################################################################
#----------------------------------CeCILL License-------------------------------
# This file is part of myProMS
#
# Copyright Institut Curie 2019
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
use POSIX qw(strftime); # to get the time
use File::Copy qw(copy move);
use File::Copy::Recursive qw(dirmove dircopy);
use File::Path qw(rmtree); # remove_tree
use File::Basename;  # To easily parse file names
use promsConfig;
use promsMod;
use promsQuantif;

#print header(-charset=>'UTF-8'); warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my %xicSoftware=('PD'=>'Proteome Discoverer','MCQ'=>'MassChroQ','MAS'=>'Mascot','PAR'=>'Paragon','PKV'=>'PeakView','MQ'=>'MaxQuant','OS'=>'OpenSwath','SKY'=>'Skyline','SPC'=>'Spectronaut','?'=>'Unknown');
#my @freeResList=('Protein Nter','A','C'..'I','K'..'N','P'..'T','V','W','Y','Protein Cter');
#my $MAX_CHANNELS=3;  # 3 max for SILAC
my $updateFrameString="";
my $maxLabPerChannel=10; # max num for isotope extraction...
my %normalizationNames=&promsQuantif::getQuantifNormalizationName;
my $MAX_PROT_MANUAL_PEP=10;

###############################
####>Recovering parameters<####
###############################
my $designID=&promsMod::cleanNumericalParameters(param('ID'));
if (param('AJAX')) {
	if (param('AJAX') eq 'fetchRefQuantif') {&ajaxFetchReferenceQuantifications;}
	elsif (param('AJAX') eq 'fetchRefProtPep') {&ajaxFetchRefProtPeptides;}
	exit;
}

############################
####>Form was submitted<####
############################
if (param('launch')) {
	&launchQuantifications;
	exit;
}

#######################################
####>Resume failed quantification<####
#######################################
if (param('RESUME')) {
	&resumeQuantification;
	exit;
}

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;
my $projectID=&promsMod::getProjectID($dbh,$designID,'design');
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
#my $projectAccess=${$userInfo[2]}{$projectID};
my $bioProjectAccessStrg=''; # for template selection
if ($userInfo[1]=~/^(bio|manag)$/) { # for template selection
	@userInfo=&promsMod::getUserInfo($dbh,$userID); # scan all accessible projects
	$bioProjectAccessStrg=join(',',keys %{$userInfo[2]});
}
my ($designName)=$dbh->selectrow_array("SELECT NAME FROM DESIGN WHERE ID_DESIGN=$designID");
my $titleString="Start Quantification from Design <FONT color=#DD0000>$designName</FONT>";
my $oldAlgoAccessStrg=''; #($ENV{HTTP_HOST}=~/curie\.fr/)? "<BR>\n<INPUT type=\"button\" class=\"title3\" value=\"Use old alogrithms\" onclick=\"window.location='./selAna4Quantification_OA.cgi?ID=design:$designID&quantifType=DESIGN'\"/>" : '';

####>Project PTMs<####
my $sthGetPM=$dbh->prepare("SELECT PM.ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DISPLAY_CODE,DISPLAY_COLOR FROM PROJECT_MODIFICATION PM,MODIFICATION M WHERE PM.ID_MODIFICATION=M.ID_MODIFICATION AND PM.ID_PROJECT=$projectID");
$sthGetPM->execute;
my %projectVarMods;
while (my ($modID,$psiName,$interName,$synName,$code,$color)=$sthGetPM->fetchrow_array) {
	@{$projectVarMods{$modID}}=($psiName || $interName,$code,$color);
	unless ($projectVarMods{$modID}[0]) {
		$synName=~s/^##//; $synName=~s/##$//; $synName=~s/##/, /;
		$projectVarMods{$modID}[0]=$synName;
	}
}
$sthGetPM->finish;

my $numPTMs=scalar keys %projectVarMods;
my ($phosphoID,$allowFreeResidues)=(0,0); # ,$freeResSpecif
my $sthPFR=$dbh->prepare('SELECT ID_MODIFICATION,UNIMOD_ACC FROM MODIFICATION WHERE UNIMOD_ACC IN (-1,21)');
$sthPFR->execute;
while (my ($ptmID,$unimod)=$sthPFR->fetchrow_array) {
	if ($unimod==21) { # Phospho
		$phosphoID=$ptmID if $projectVarMods{$ptmID}; # enabled only if used in project;
	}
	else { # Free residues (fake PMT)
		$allowFreeResidues=1;
		# $freeResSpecif=$specif || '';
		$numPTMs++;
	}
}
$sthPFR->finish;

########################################################
####>Recovering list of experiment, sample analysis<####
########################################################
my (%listDataBank,@itemAnalyses,%anaProteins,%listParam,%anaLabelMods,%modifications,%anaLabeling,@referenceDesign); #,%refQuantifications %okRemFilter,%okActLowScores,%listQuantif,
####>Labeled quantifs<####
#my (%anaPeptideQuantifs,%quantifChannels,$numChannels,$maxReplicates); # SILAC, iTRAQ internal quanti
my (%categoryList,%quantifMethods,%quantifTemplate,%projectTemplate,%defaultAnaFormSelValues);
foreach my $quantMethod ('NONE','XIC','SILAC','ITRAQ','TMT','MQ','DIA','TDA') {$quantifMethods{$quantMethod}=0;}
my (%designInfo,%anaObs); # For TnPQ and Prop Pep Ratio (quantifications associated to a DESIGN)
my $nbCond=0;
my $numProtInDesign=0;
my ($selCondStrg1,$selCondStrg2)=("<OPTION value=\"Select\">-= Select =-</OPTION>","<OPTION value=\"Select\">-= Select =-</OPTION>");
my $varGlobSelCond="var stateCondVal=[];\n"; # JS array to record Index of selected cond for each State
my $jsXicSoftStrg; # for DESIGN only
my ($ptmProbSoft,$ptmProbSoftCode,$isMultiPtmProb)=('','',0); #('PhosphoRS','PRS'); # default

my $sthAE=$dbh->prepare("SELECT E.NAME,E.ID_EXPCONDITION,BS.NAME,O.ID_ANALYSIS,O.TARGET_POS,OE.FRACTION_GROUP,OE.TECH_REP_GROUP,O.ID_OBSERVATION,GROUP_CONCAT(M.ID_MODIFICATION ORDER BY M.ID_MODIFICATION SEPARATOR ':')
							FROM OBSERVATION O
							LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION
							LEFT JOIN BIOSAMPLE BS ON O.ID_BIOSAMPLE=BS.ID_BIOSAMPLE
							JOIN OBS_EXPCONDITION OE ON O.ID_OBSERVATION=OE.ID_OBSERVATION
							JOIN EXPCONDITION E ON OE.ID_EXPCONDITION=E.ID_EXPCONDITION
							WHERE E.ID_DESIGN=$designID GROUP BY O.ID_OBSERVATION,E.ID_EXPCONDITION ORDER BY E.DISPLAY_POS,E.ID_EXPCONDITION");
my $expCondID=0;
my %anaList;

$sthAE->execute;
while (my ($expConditionName,$expConditionID,$bioSampName,$anaID,$targetPos,$fracGroup,$techRepGroup,$obsID,$modCode) = $sthAE->fetchrow_array) {

	$targetPos=0 unless $targetPos; # just to be safe
	if ($expCondID ne $expConditionID) {
		$selCondStrg1.="<OPTION id=\"$expConditionName\" value=\"$expConditionID\">$expConditionName</OPTION>";
		$selCondStrg2.="<OPTION id=\"$expConditionName\" value=\"$expConditionID\" disabled>$expConditionName</OPTION>";
		$varGlobSelCond.="stateCondVal[$nbCond]=0;\n";
		$nbCond++;
		$expCondID=$expConditionID;
		$designInfo{'EXPCONDITION'}{$expConditionID}{'NAME'}=$expConditionName;
		$designInfo{'EXPCONDITION'}{$expConditionID}{'DISPLAY_POS'}=$nbCond;
	}
	my $obsCode=($modCode)? "$anaID:$targetPos:$modCode" : "$anaID:$targetPos";
	$designInfo{'EXPCONDITION'}{$expConditionID}{'ANALYSIS'}{$anaID}=1;
	$designInfo{'EXPCONDITION'}{$expConditionID}{'OBSERVATION'}{$obsCode}=1;
	$designInfo{'OBSERVATION'}{$obsCode}{'ID_EXPCONDITION'}=$expConditionID;
	$designInfo{'OBSERVATION'}{$obsCode}{'BIOSAMPLE'}=$bioSampName || 'Unknown Biosample';
	push @{$anaObs{$expConditionID}{$anaID}},[$obsID,$obsCode,$fracGroup,$techRepGroup];
	$anaList{$anaID}=1;
}
$sthAE->finish;

###> Get ANALYSIS Information
my $anaString='('.join(',',keys %anaList).')';
my @sthGetAnaInfo=( # 2d-gel or sample
	$dbh->prepare("SELECT ID_ANALYSIS,A.NAME,CONCAT(G.NAME,'&nbsp;>&nbsp;',SP.NAME,'&nbsp;>&nbsp;',A.NAME) FROM ANALYSIS A,SAMPLE S,SPOT SP,GEL2D G,EXPERIMENT E,DESIGN D WHERE D.ID_EXPERIMENT=E.ID_EXPERIMENT AND E.ID_EXPERIMENT=G.ID_EXPERIMENT AND G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT AND S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS IN $anaString ORDER BY G.DISPLAY_POS,SP.NAME,S.DISPLAY_POS,A.DISPLAY_POS"),
	$dbh->prepare("SELECT ID_ANALYSIS,A.NAME,CONCAT(S.NAME,'&nbsp;>&nbsp;',A.NAME) FROM ANALYSIS A,SAMPLE S,EXPERIMENT E,DESIGN D WHERE D.ID_EXPERIMENT=E.ID_EXPERIMENT AND E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SAMPLE=A.ID_SAMPLE AND ID_SPOT IS NULL AND A.ID_ANALYSIS IN $anaString ORDER BY S.DISPLAY_POS,A.DISPLAY_POS")
);
#my $sthObs=$dbh->prepare("SELECT ID_ANALYSIS,VALID_STATUS,A.NAME FROM ANALYSIS A,ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND FOCUS='peptide' AND A.ID_ANALYSIS IN $anaString ORDER BY A.DISPLAY_POS,Q.NAME");
my $anaPos=0;
my ($firstAnalysis,$lastAnalysis)=(0,0);
foreach my $sthAnaInfo (@sthGetAnaInfo) {
	$sthAnaInfo->execute;
	while (my ($anaID,$anaName,$hierarchy) = $sthAnaInfo->fetchrow_array)  {
		$designInfo{'ANALYSIS'}{$anaID}{'HIERARCHY'}=$hierarchy;
		$designInfo{'ANALYSIS'}{$anaID}{'NAME'}=$anaName;
		$designInfo{'ANALYSIS'}{$anaID}{'POS'}=++$anaPos;
		#$designInfo{'ANALYSIS'}{$anaID}{'ANA.DISPLAY_POS'}=$anaDPos;
		#$designInfo{'ANALYSIS'}{$anaID}{'SAMP.DISPLAY_POS'}=$sampDPos;
		$firstAnalysis=$anaID unless $firstAnalysis;
		$lastAnalysis=$anaID;
	}
	$sthAnaInfo->finish;
}

###> Get QUANTIFICATION ANALYSIS
#my ($xicIDMethod)=$dbh->selectrow_array('SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE=\'XIC\'');

#my $sthQInfo=$dbh->prepare("SELECT Q.ID_QUANTIFICATION,NAME,QUANTIF_ANNOT,STATUS,ID_ANALYSIS,UPDATE_DATE FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND FOCUS='peptide' AND ID_QUANTIFICATION_METHOD=$xicIDMethod AND ID_ANALYSIS IN $anaString");
my $sthQInfo=$dbh->prepare("SELECT Q.ID_QUANTIFICATION,Q.NAME,QUANTIF_ANNOT,STATUS,ID_ANALYSIS,UPDATE_DATE,QM.CODE,QM.DES FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ,QUANTIFICATION_METHOD QM WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND FOCUS='peptide' AND STATUS >=1 AND CODE != 'SIN' AND ID_ANALYSIS IN $anaString");
my $sthIsobar=$dbh->prepare("SELECT DISTINCT M.ID_MODIFICATION FROM MODIFICATION M,ANALYSIS_MODIFICATION AM WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND IS_LABEL=1 AND ID_ANALYSIS=? LIMIT 1"); #  AND AM.MODIF_TYPE='V'
my $sthAna=$dbh->prepare("SELECT FILE_FORMAT FROM ANALYSIS WHERE ID_ANALYSIS=?");
my $hasPhospho=($phosphoID > 0)? 0 : -1;
$sthQInfo->execute;
while (my ($quantiID,$quantiName,$quantifAnnot,$qStatus,$anaID,$date,$methodCode,$methodDes) = $sthQInfo->fetchrow_array) {
#$quantiName=~s/3plex.+\| /2plex (/;
	next unless $quantifAnnot;
#print "$quantiID,/$quantiName,/$qStatus,/$anaID,/'$methodCode'<BR>\n";
	unless ($hasPhospho) {
		my ($okPhospho)=$dbh->selectrow_array("SELECT 1 FROM ANALYSIS_MODIFICATION WHERE ID_ANALYSIS=$anaID AND ID_MODIFICATION=$phosphoID");
		$hasPhospho=$okPhospho || -1;
	}
	$quantifMethods{$methodCode}=1;
	my ($xicSoftCode,$xicSoftVersion);
	$quantifAnnot=~s/::SOFTWARE=(\w+);?(\d|\.)*//; # remove software info for back compatibility
	if ($1) {
		$xicSoftCode=$1;
		$xicSoftVersion=$2;
	}
	elsif ($quantifAnnot=~/EXTRACTION_ALGO=/) {$xicSoftCode='MCQ';}
	else {
		$sthAna->execute($anaID);
		my ($fileFormat)=$sthAna->fetchrow_array;
		$xicSoftCode=($fileFormat=~/\.PDM\Z/)? 'PD' : ($fileFormat=~/^MASCOT/)? 'MAS' : ($fileFormat=~/PARAGON/)? 'PAR' : '?';
	}
	if ($xicSoftCode eq 'MQ') {
		($ptmProbSoft,$ptmProbSoftCode,$isMultiPtmProb)=('MaxQuant','MQ',1);
	}
	else {
		unless ($ptmProbSoftCode) {
			my ($isPRS)=($hasPhospho > 0)? $dbh->selectrow_array("SELECT 1 FROM PEPTIDE WHERE ID_ANALYSIS=$anaID AND DATA LIKE '%PRS=%' LIMIT 1") : (0);
			if ($isPRS) {
				($ptmProbSoft,$ptmProbSoftCode)=('PhosphoRS','PRS');
			}
			else {
				my ($modifID,$probStrg)=$dbh->selectrow_array("SELECT ID_MODIFICATION,REF_POS_STRING FROM PEPTIDE_MODIFICATION PM INNER JOIN PEPTIDE P ON P.ID_PEPTIDE=PM.ID_PEPTIDE WHERE ID_ANALYSIS=$anaID AND REF_POS_STRING LIKE '%PRB_%' LIMIT 1");
				if ($probStrg && $probStrg=~/PRB_(\w+)=/) {
					$ptmProbSoftCode=$1;
					$ptmProbSoft=($ptmProbSoftCode eq 'PTMRS')? 'ptmRS' : $xicSoftware{$ptmProbSoftCode} || $ptmProbSoftCode;
					if ($modifID != $phosphoID) {$isMultiPtmProb=1;} # Assume multi-modifs
					else {
						#my $sthMulti=$dbh->prepare("SELECT DISTINCT ID_MODIFICATION FROM PEPTIDE_MODIFICATION PM INNER JOIN PEPTIDE P ON P.ID_PEPTIDE=PM.ID_PEPTIDE WHERE ID_ANALYSIS=$anaID AND REF_POS_STRING LIKE '%PRB_PTMRS%'");
						#$sthMulti->execute;
						#my $numModif=$sthMulti->rows;
						#$isMultiPtmProb=($numModif > 1)? 1 : 0;
						my ($nonPhosphoModifID)=$dbh->selectrow_array("SELECT ID_MODIFICATION FROM PEPTIDE_MODIFICATION PM INNER JOIN PEPTIDE P ON P.ID_PEPTIDE=PM.ID_PEPTIDE WHERE ID_ANALYSIS=$anaID AND ID_MODIFICATION != $phosphoID AND REF_POS_STRING LIKE '%PRB_$ptmProbSoftCode%' LIMIT 1");
						$isMultiPtmProb=($nonPhosphoModifID)? 1 : 0;
					}
				}
			}
			$ptmProbSoftCode='-' unless $ptmProbSoftCode;
		}
		if ($xicSoftCode eq 'SPC') {
			($xicSoftVersion)=$quantifAnnot=~/SOFTWARE_VERSION=([^:]+)/ unless $xicSoftVersion;
		}
	}
	my $labelType=($quantifAnnot=~/LABEL=([^:;]+)/)? $1 : 'FREE';
	$labelType=uc($labelType);
	if ($labelType=~/FREE|NONE/) {
		$labelType='FREE';
		my $targetPos=0;
		if ($designInfo{'OBSERVATION'}{"$anaID:$targetPos"} && $designInfo{'OBSERVATION'}{"$anaID:$targetPos"}{'ID_EXPCONDITION'}) {
			$designInfo{'OBSERVATION'}{"$anaID:$targetPos"}{'NAME'}='Label-free';
			#$designInfo{'QUANTIFICATION'}{$quantiID}{TARGET}{0}='Label-free';
			$designInfo{'OBSERVATION'}{"$anaID:$targetPos"}{'ID_QUANTIFICATION'}{$quantiID}=$targetPos; # true channel pos used in quanti
			#$designInfo{'EXPCONDITION'}{ $designInfo{'OBSERVATION'}{"$anaID:$targetPos"}{'ID_EXPCONDITION'} }{'ID_QUANTIFICATION'}{$quantiID}=$targetPos; # true channel pos used in quanti
		}
		else {next;} # skip quantif without matching(label) observation
	}
	#if ($quantifAnnot=~/EXTRACTION_ALGO=/) { # data generated by MassChroQ XIC: LABEL=FREE,SILAC
	#}
	elsif ($labelType=~/SILAC|ITRAQ|TMT/) { # data imported from search file
		my ($labelTypeStrg,@labelInfo)=split('::',$quantifAnnot);
		#my ($labelType)=($labelTypeStrg=~/LABEL=(.+)/);
		my $isobarModID;
		foreach my $infoStrg (@labelInfo) {
			last if $infoStrg !~ /;/; # MassChroQ parameters
			my ($targetPos,$chanName,$labelStrg)=split(';',$infoStrg);
			my ($obsCode,$targetName);
			if ($labelType eq 'SILAC') {
				$obsCode="$anaID:-1"; # SILAC: no longer relies on target pos in OBSERVATION (PP 02/14)
				$targetName="$chanName ";
				my $first=1;
				#my @mods;
				my %mods; # hash because same mod can be used on different aa (mod is repeated)
				foreach my $modStrg (split ('@',$labelStrg)) { # multiple mods for same label channel
					my @labelData=split('#',$modStrg);
					if ($first) {$first=0;}
					else {$targetName.='+';}
					#if (!$labelData[1] || $labelData[1] =~ /No label/i) { #} #TODO: Fix Guillaume's MCQ XIC annotation!!!!!!!!
					if ($modStrg =~ /No label/i) {
						$targetName.="[No label]";
					}
					else {
						$targetName.="[$labelData[1]&raquo;$labelData[2]]"; #&raquo; &middot; &rarr;
						my $modID=($labelData[4] && $labelData[4]>0)? $labelData[4] : &promsMod::getModificationIDfromString($dbh,$labelData[1]); # $labelData[4] if MassChroQ XIC (-1 if no label!!!)
						$modID=0 unless $modID;
						#push @mods,$modID;
						$mods{$modID}=1;
					}
				}
				if (scalar keys %mods) { # label channel
					$obsCode.=':'.join(':',sort{$a<=>$b} keys %mods); # order is needed for later compare with %expCondition
				}
#print "Q=$quantiID, OBS='$obsCode' ($targetPos) '$labelStrg'<BR>\n";
			}
			elsif ($labelType=~/ITRAQ|TMT/) {
				$obsCode="$anaID:$targetPos"; # iTRAQ relies on target pos in OBSERVATION
				$targetName=$chanName;
				unless ($isobarModID) { # only once
					$sthIsobar->execute($anaID);
					($isobarModID)=$sthIsobar->fetchrow_array;
				}
				$obsCode.=":$isobarModID" if $isobarModID;
			}
			if ($designInfo{'OBSERVATION'}{$obsCode} && $designInfo{'OBSERVATION'}{$obsCode}{'ID_EXPCONDITION'}) {
				$designInfo{'OBSERVATION'}{$obsCode}{'NAME'}=$targetName;
				$designInfo{'OBSERVATION'}{$obsCode}{'ID_QUANTIFICATION'}{$quantiID}=$targetPos; # true channel pos used in quanti
			}
			else {next;} # skip quantif without matching(label) observation
		}
	}
	$designInfo{'LABEL'}{$labelType}=1;
	$designInfo{'QUANTIFICATION'}{$quantiID}{'NAME'}=$quantiName.' ['.$methodCode.']';
	$designInfo{'QUANTIFICATION'}{$quantiID}{'METHOD'}=$methodCode;
	$designInfo{'QUANTIFICATION'}{$quantiID}{'ANNOT'}=$quantifAnnot;
	$designInfo{'QUANTIFICATION'}{$quantiID}{'STATUS'}=$qStatus;
	$designInfo{'QUANTIFICATION'}{$quantiID}{'DATE'}=$date;
	$designInfo{'QUANTIFICATION'}{$quantiID}{'LABEL'}=$labelType;
	$designInfo{'QUANTIFICATION'}{$quantiID}{'XIC_SOFTWARE'}=$xicSoftCode;
	$designInfo{'QUANTIFICATION'}{$quantiID}{'ANALYSIS'}{$anaID}=1;
	$designInfo{'ANALYSIS'}{$anaID}{'LABEL'}{$labelType}{$quantiID}=1;
	$designInfo{'ANALYSIS'}{$anaID}{'HAS_QUANTIF'}=1;
	my %quantifParameters;
	foreach my $parameter (split(/::/,$quantifAnnot) ) {
		next unless $parameter;
		my ($parameterName,$parameterValue)=split(/=/,$parameter);
		$quantifParameters{$parameterName}=$parameterValue;
	}
	my $popupInfo="<U><B>$methodCode:</B></U><BR>$methodDes<BR><U><B>Quantification parameters:</B></U>";
	if ($methodCode eq 'XIC') {
		$popupInfo.="<BR><B>Software:</B> $xicSoftware{$xicSoftCode}";
		$popupInfo.=" v.$xicSoftVersion" if $xicSoftVersion;
		if ($xicSoftCode eq 'PD') {
			foreach my $paramName (sort{lc($a) cmp lc($b)} keys %quantifParameters) {
				next if $paramName eq 'LABEL';
				$popupInfo.="<BR><B>$paramName:</B> $quantifParameters{$paramName}";
			}
		}
		elsif ($xicSoftCode eq 'MQ') { # MaxQuant
			# Nothing for now
		}
		elsif ($xicSoftCode eq 'SKY') { # Skyline
			# Nothing for now
		}
		elsif ($xicSoftCode eq 'SPC') { # Spectronaut
			$popupInfo.="<BR><B>Library:</B> $quantifParameters{LIBRARY_NAME}";
		}
		elsif ($xicSoftCode eq 'MCQ') { # assume MCQ
			$popupInfo.="<BR><B>Extraction type:</B> $quantifParameters{RAWDATA_ACQUISITION} (for mzXML)";
			$popupInfo.="<BR><B>Type of XIC:</B> $quantifParameters{XIC_EXTRACTION_TYPE}";
			$popupInfo.='<BR><B>Alignment algo.:</B> ';
			if ($quantifParameters{'EXTRACTION_ALGO'} eq 'OBI') {
				$popupInfo.='OBI-Warp';
				$popupInfo.="<BR><B>Align from</B> $quantifParameters{MZ_ALIGN_RANGE_MIN} <B>to</B> $quantifParameters{MZ_ALIGN_RANGE_MAX}";
			}
			else { # ms2
				$popupInfo.='ms2';
				$popupInfo.="<BR><B>Tendency:</B> $quantifParameters{MS2_TENDENCY}, <B>MS2 smoothing:</B> $quantifParameters{MS2_SMOUTHING}, <B>MS1 smoothing:</B> $quantifParameters{MS1_SMOUTHING}";
			}
			unless ($designInfo{ANALYSIS}{$quantifParameters{REFERENCE}}{NAME}) { # Reference may not be part of current design!
				($designInfo{ANALYSIS}{$quantifParameters{REFERENCE}}{NAME})=$dbh->selectrow_array("SELECT NAME FROM ANALYSIS WHERE ID_ANALYSIS=$quantifParameters{REFERENCE}");
			}
			$popupInfo.="<BR><B>Reference:</B> $designInfo{ANALYSIS}{$quantifParameters{REFERENCE}}{NAME}";
		}
	}
	elsif ($methodCode eq 'DIA') {
		if ($xicSoftCode eq 'PKV') {
			my $modExclStrg=($quantifParameters{EXCLUDE_MODIFIED_PEPTIDES} eq 'NOT CHECKED' || $quantifParameters{EXCLUDE_MODIFIED_PEPTIDES}==0)? ' not' : '';
			$popupInfo.="<BR><B>Modified peptides are$modExclStrg excluded</B>";
			$popupInfo.="<BR><B>Num. transitions per peptide:</B> $quantifParameters{NB_TRANSITION_PER_PEPTIDE}";
			$popupInfo.="<BR><B>XIC extraction window:</B> $quantifParameters{XIC_EXTRACTION_WINDOW} min.";
			$quantifParameters{XIC_EXTRACTION_WINDOW}=~s/_DA/ Da/i;
			$quantifParameters{XIC_EXTRACTION_WINDOW}=~s/_PPM/ ppm/i;
			$popupInfo.="<BR><B>XIC width:</B> $quantifParameters{XIC_EXTRACTION_WINDOW}";
		}
		elsif ($xicSoftCode eq 'SPC') { # Spectronaut
			$popupInfo.="<BR><B>Software:</B> $xicSoftware{$xicSoftCode}";
			$popupInfo.=" v. $xicSoftVersion" if $xicSoftVersion;
			$popupInfo.="<BR><B>Library:</B> $quantifParameters{LIBRARY_NAME}";
		}
	}
	else {
		$popupInfo.="<BR>Provided by search results file";
	}
	$designInfo{'QUANTIFICATION'}{$quantiID}{'ANNOT_STRING'}=$popupInfo;
}
$sthQInfo->finish;
$sthIsobar->finish;
$sthAna->finish;

unless ($designInfo{'LABEL'}) { # no quantif at all
	$designInfo{'LABEL'}{'FREE'}=1;
	$quantifMethods{'NONE'}=1;
}

###>JS list of quantif=>XIC software
$jsXicSoftStrg="var quantiXicSoftware={\n";
my $firstSoft=1;
foreach my $quantiID (sort{$a<=>$b} keys %{$designInfo{'QUANTIFICATION'}}) {
	$jsXicSoftStrg.=",\n" unless $firstSoft;
	$jsXicSoftStrg.="\t$quantiID:'$designInfo{QUANTIFICATION}{$quantiID}{XIC_SOFTWARE}'";
	$firstSoft=0;
}
$jsXicSoftStrg.="\n};\n";

###>Reference quantif management for intra-protein normalization (only used for modif quantifs)
my $sthRefE=$dbh->prepare("SELECT DISTINCT D.ID_DESIGN,CONCAT(E.NAME,' > ',D.NAME),E.DISPLAY_POS,D.NAME FROM EXPERIMENT E
							INNER JOIN DESIGN D ON E.ID_EXPERIMENT=D.ID_EXPERIMENT
							INNER JOIN QUANTIFICATION Q ON D.ID_DESIGN=Q.ID_DESIGN
							INNER JOIN QUANTIFICATION_METHOD QM ON Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD
							WHERE E.ID_PROJECT=$projectID AND Q.ID_MODIFICATION IS NULL AND Q.FOCUS='protein' AND STATUS=1 AND QM.CODE='PROT_RATIO_PEP' ORDER BY E.DISPLAY_POS,D.NAME");
$sthRefE->execute;
while (my($desID,$desFullName)=$sthRefE->fetchrow_array) {
	push @referenceDesign,[$desID,$desFullName];
}
$sthRefE->finish;

###>Check if manual peptide selection is possible (<= $MAX_PROT_MANUAL_PEP protein in design) (Too slow at design-level => check 1st & last analysis)
#my $sthNumProt=$dbh->prepare("SELECT DISTINCT(ID_PROTEIN) FROM EXPCONDITION E
# 									INNER JOIN OBS_EXPCONDITION OE ON E.ID_EXPCONDITION=OE.ID_EXPCONDITION
# 									INNER JOIN OBSERVATION O ON OE.ID_OBSERVATION=O.ID_OBSERVATION
# 									INNER JOIN ANALYSIS_PROTEIN AP ON O.ID_ANALYSIS=AP.ID_ANALYSIS
# 									WHERE ID_DESIGN=? AND ABS(AP.VISIBILITY)=2 LIMIT ".($MAX_PROT_MANUAL_PEP+1)); # +1 to check if >$MAX_PROT_MANUAL_PEP
# $sthNumProt->execute($designID);
my $sthNumProt=$dbh->prepare("SELECT DISTINCT(ID_PROTEIN) FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS IN ($firstAnalysis,$lastAnalysis)");
$sthNumProt->execute;
$numProtInDesign=$sthNumProt->rows;
$sthNumProt->finish;

###>List of non-label modification used
my $sthMod=$dbh->prepare("SELECT DISTINCT M.ID_MODIFICATION,M.PSI_MS_NAME,M.INTERIM_NAME,M.DES,M.SYNONYMES FROM MODIFICATION M
							INNER JOIN ANALYSIS_MODIFICATION AM ON M.ID_MODIFICATION=AM.ID_MODIFICATION
							INNER JOIN OBSERVATION O ON AM.ID_ANALYSIS=O.ID_ANALYSIS
							INNER JOIN OBS_EXPCONDITION OC ON O.ID_OBSERVATION=OC.ID_OBSERVATION
							INNER JOIN EXPCONDITION C ON OC.ID_EXPCONDITION=C.ID_EXPCONDITION
							WHERE ID_DESIGN=? AND IS_LABEL != 1 AND AM.MODIF_TYPE='V'"
						);
$sthMod->execute($designID);
while (my ($modID,$psiName,$interName,$des,$synName,$specifStrg)=$sthMod->fetchrow_array) {
	$modifications{$modID}=$psiName || $interName || $des;
	unless ($modifications{$modID}) {
		$synName=~s/^##//; $synName=~s/##$//; $synName=~s/##/, /;
		$modifications{$modID}=$synName;
	}
}
$sthMod->finish;

my %unusedProjVarMods;
my %terminiCodes=('-'=>'Protein N-term','='=>'Any N-term','+'=>'Protein C-term','*'=>'Any C-term');
my $sthGetSpecif=$dbh->prepare("SELECT SPECIFICITY FROM MODIFICATION WHERE ID_MODIFICATION=?");
foreach my $modID (keys %projectVarMods) {
	unless ($modifications{$modID}) {
		$sthGetSpecif->execute($modID);
		my ($specifStrg)=$sthGetSpecif->fetchrow_array;
		foreach my $resCode (split(',',$specifStrg)) {
			next unless $resCode; # problem eg 'C1,C2,,C3,...'
			my ($res,$context)=split(';',$resCode);
			next if length($res) != 1; # unexpected problem
			my $dispRes;
			if ($res=~/\w/) {
				$dispRes=$res;
				$dispRes.=" ($terminiCodes{$context})" if $context;
			}
			else {
				$dispRes=$terminiCodes{$res};
			}
			$resCode=~s/;/,/; # ';' is not compatible with QUANTIF_ANNOT data structure
			push @{$unusedProjVarMods{$modID}},[$resCode,$dispRes];
		}
	}
}
$sthGetSpecif->finish;

####>Custom lists (Categories)<####
my $sthCat=$dbh->prepare("SELECT CL.ID_CLASSIFICATION,CL.NAME,ID_CATEGORY,CA.NAME,CA.LIST_TYPE FROM CLASSIFICATION CL,CATEGORY CA WHERE CL.ID_CLASSIFICATION=CA.ID_CLASSIFICATION AND ID_PROJECT=$projectID ORDER BY CL.NAME ASC,CA.DISPLAY_POS ASC");
$sthCat->execute;
while (my ($classID,$className,$catID,$catName,$listType)=$sthCat->fetchrow_array) {
	$listType='PROT' unless $listType;
next unless $listType eq 'PROT'; # restrict to proteins due to id matching problem in R if sites vs sites
	$categoryList{$classID}{'CL_NAME'}=$className unless $categoryList{$classID};
	push @{$categoryList{$classID}{'CAT'}},[$catID,$catName,$listType];
}
$sthCat->finish;

###> Get TEMPLATES Information
#my $sthCatTemplate=$dbh->prepare("SELECT CL.ID_CLASSIFICATION,CL.NAME,ID_CATEGORY,CA.NAME,CA.LIST_TYPE FROM CLASSIFICATION CL,CATEGORY CA WHERE CL.ID_CLASSIFICATION=CA.ID_CLASSIFICATION AND ID_PROJECT=$projectID AND ID_CATEGORY=? ORDER BY CL.NAME ASC,CA.DISPLAY_POS ASC");
#my $sthQTemplateInfo=$dbh->prepare("SELECT ID_MODIFICATION,ID_QUANTIFICATION,NAME,QUANTIF_ANNOT FROM QUANTIFICATION WHERE QUANTIF_ANNOT LIKE '%::TEMPLATE=1%'");
my $sthQTemplateInfo=($userInfo[1] eq 'bio')?
	$dbh->prepare("SELECT E.ID_PROJECT,ID_MODIFICATION,ID_QUANTIFICATION,Q.NAME,QUANTIF_ANNOT FROM QUANTIFICATION Q,DESIGN D,EXPERIMENT E WHERE Q.ID_DESIGN=D.ID_DESIGN AND D.ID_EXPERIMENT=E.ID_EXPERIMENT AND E.ID_PROJECT IN ($bioProjectAccessStrg) AND QUANTIF_ANNOT LIKE '%::TEMPLATE=1%'")
	: $dbh->prepare("SELECT E.ID_PROJECT,ID_MODIFICATION,ID_QUANTIFICATION,Q.NAME,QUANTIF_ANNOT FROM QUANTIFICATION Q,DESIGN D,EXPERIMENT E WHERE Q.ID_DESIGN=D.ID_DESIGN AND D.ID_EXPERIMENT=E.ID_EXPERIMENT AND QUANTIF_ANNOT LIKE '%::TEMPLATE=1%'");
my $sthProjName=$dbh->prepare('SELECT NAME FROM PROJECT WHERE ID_PROJECT=?');
$sthQTemplateInfo->execute;
while(my ($projID,$modID,$quantiID,$quantiName,$quantifAnnot)=$sthQTemplateInfo->fetchrow_array){
	my $labelType=($quantifAnnot=~/LABEL=([^:;]+)/)? $1 : 'FREE';
	next unless $designInfo{'LABEL'}{$labelType}; # Restrict to same labeling(s)!!!!

	#>Remove useless data
	$quantifAnnot=~s/::/&/g;
	$quantifAnnot=~s/&RATIOS=[^&]+//;
	$quantifAnnot=~s/&STATES=[^&]+//;
	$quantifAnnot=~s/&/::/g;

	next unless (($quantifAnnot !~ /SOFTWARE/ && $quantifAnnot !~ /ALGO_TYPE/ ) || ($quantifAnnot =~ /SOFTWARE=DIA/ && $quantifMethods{'DIA'}) || ($quantifAnnot =~ /SOFTWARE=SKY/ && $quantifMethods{'TDA'}) || ($quantifAnnot !~ /SOFTWARE=DIA/ && ($quantifMethods{'XIC'} || $quantifMethods{'SILAC'} || $quantifMethods{'ITRAQ'} || $quantifMethods{'TMT'})));
##next if ($quantifAnnot=~/ALGO_TYPE=PEP_RATIO/ && $quantifMethods{'XIC'}); # PP: not good for label

	##> check project modifications compatibility
	next if ($modID && !$projectVarMods{$modID});

	##> check protein list compatibility
	my $listID=($quantifAnnot=~/PROTEINS=\w+;#(\d)/)? $1: undef;
	#if ($listID){
	#	$sthCatTemplate->execute($listID);
	#	next if $sthCatTemplate->rows==0;
	#}
	next if ($listID && $projID != $projectID); # cannot use a list from another project

	##> check ana modifications compatibility
	my $peptideInfo=($quantifAnnot=~/PEPTIDES=([^A-Z]+)/)? $1: '';
	my @peptideInfos=split(/;/,$peptideInfo);
	if ($peptideInfos[2]=~/#/){
		$peptideInfos[2]=~s/\d://;
		$peptideInfos[2]=~s/\D//g;
		my @modIDs=split(//,$peptideInfos[2]);
		my $notMatch=1;
		foreach my $pepModID (@modIDs){
			unless ($modifications{$pepModID}){
				$notMatch=0;
				last;
			}
		}
		next unless $notMatch;
	}
	$quantifTemplate{$labelType}{$projID}{$quantiID}{'NAME'}=$quantiName;
	$quantifTemplate{$labelType}{$projID}{$quantiID}{'QUANTIFANNOT'}=$quantifAnnot;
	$quantifTemplate{$labelType}{$projID}{$quantiID}{'MODID'}=$modID;

	unless ($projectTemplate{$projID}) {
		$sthProjName->execute($projID);
		($projectTemplate{$projID})=$sthProjName->fetchrow_array;
	}
}
$sthQTemplateInfo->finish;
$sthProjName->finish;
#}

$dbh->disconnect;

#######################
####>Starting HTML<####
#######################
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
print header(-charset=>'utf-8'); # because of ° in source selection
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Start Design Quantification</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TD.center {text-align:center}
INPUT.keyword {font-weight:bold; font-size:11px;}
.bold {font-weight:bold;}
</STYLE>
<SCRIPT type="text/javascript">
// Dictionnary of help poupup texts (avoid using '"')
const helpText= {
	trackMS2:	'Also track number of truely identified peptides (by MS/MS)<BR>(unlike those recovered by MBWR: Match Between or Within Runs).',
	seqContext:
`Data from the same site on differently cleaved region will <B>not</B> be aggregated.<BR>
For instance, <B>Phosphorylation on Ser53</B> in<BR>
&nbsp;<B><SUB>51</SUB>LAS*PELER<SUB>58</SUB></B> and<BR>
&nbsp;<B><SUB>51</SUB>LAS*PELERAMLK<SUB>62</SUB></B> will be computed separately.`,
	scopePTM:
`-<B>Whole proteins:</B> Selected peptides will be used to quantify proteins<BR>
-<B>Whole proteins by PTM:</B> Peptides matching selected PTM(s) and context(s) will be used to quantify proteins<BR>
-<B>PTM sites:</B> PTM sites will be quantified instead whole proteins,<BR>
-<B>Free residues:</B> Selected residues will be quantified unless modified by a PTM</B>`,
	// freeResidues:	'Selected protein residues will be quantified unless modified by a PTM.',
	contextPTM:
`Restrict to peptides with PTMs matching this <B>context/motif</B> (Leave empty for no restriction).<BR>
&nbsp;&bull;Use <B>lowercase</B> to identify modified residue(s) and <B>uppercase</B> for other residue(s).<BR>
&nbsp;&bull;Use "<B>X</B>" or "<B>x</B>" for any residue.<BR>
&nbsp;&bull;Use "<B>/</B>" to separate multiple possible residues at same position.<BR>
&nbsp;&bull;Use "<B>{num}</B>" (or "<B>{min,max}</B>") for (variable) number of repeated residues.<BR>
&nbsp;&bull;Use "<B>-</B>" as position separator for lisibility (optional).<BR>
&nbsp;&bull;Use "<B>&gt;</B>" for protein N-term and "<B>&lt;</B>" for protein C-term.<BR>
&nbsp;&bull;Use " " (blank space/tab) to separate multiple contexts for the same PTM.<BR>
Examples: "<B>s/t</B>", "<B>&gt;x</B>", "<B>n-X-S/T</B>" (same as "<B>nXS/T</B>"), "<B>R-P/L-X{1,2}-e/d-A/G{3}-P</B>".`,
	delocalized:
`If "<B>delocalized</B>" is selected, sites with position confidence lower than threshold will be <B>delocalized</B> on peptide sequence:<BR>
eg. for <B><SUB>35</SUB>NAP<U>S</U>*GG<U>S</U>DR<SUB>43</SUB></B>:<BR>
&nbsp;&nbsp;-normal case: "<B>S38</B>"<BR>
&nbsp;&nbsp;-delocalized: "<B>38~41:1/2</B>" (site is between residues 38 and 41, 1 PTM for 2 acceptors).<BR>
<B>WARNING:</B> A site will be delocalized across the <B>entire</B> dataset if it does not pass the threshold in <B>one</B> sample.`,
	protLevelCorr:	'Corrects for variations due to whole protein fold-change.<BR>Not compatible with <B>Bias correction</B> with manual selection of peptides.',
	specificity:
`<B>For each protein:</B><BR>
&nbsp;&bull;<B>Proteotypic:</B> Only peptides matching a <B>single</B> protein are used.<BR>
&nbsp;&bull;<B>Proteotypic & others:</B> All peptides are used if at least 1 proteotypic exists.<BR>
&nbsp;&bull;<B>All:</B> All available peptides are used.`,
	modifications:	'Peptides carrying <B>allowed</B> modifications will contribute to quantification.',
	MGSharedPep:
`<B>Assign to best:</B> Peptides will be assigned to <B>only 1</B> (best) associated Match Group (based on proteotypicity, number of peptides and occurence in dataset)<BR>
-<B>Share:</B> Peptides will be used by <B>all</B> associated Match Groups<BR>
-<B>Exclude:</B> Peptides will <B>not</B> be used for quantification.`,
	includeMBWR:	'Include data from peptides identified by <B>Match Between or Within Run</B>. Not compatible with <B>Peptide/spectrum matches</B>.', // SSPA only
	protSelection:	'Restrict quantification to a subset of proteins.',
	exclContaminants:'Applies only to proteins identified from a contaminant databank.',
	visHidProt:		'A protein can be <B>visible</B> in some samples but <B>hidden</B> in others due to <B>Match group</B> creation rules.<BR>Ignoring hidden proteins can potentially result in miss-leading <B>present/absent</B> SSPA status.',
	normSharedPep:	'Only data from peptide ions found in <B>all</B> states will be considered for bias estimate.<BR>Not recommanded when numerous states are used.',
	biasCorrLevel:	'Bias correction can be performed at <B>different levels</B> for <B>Abundance</B> quantification</B> (Except for <B>Absolute</B> quantification).',
	avoidInfRatio:	'Compute a ratio if at least 1 peptide ion is shared across states regardless of the proportion of non-shared ions.',
	sspaThresholds:	'Thresholds for assigning a protein/site to a specific (set of) State(s)',
	
	// absolute quantification help
	pqiMethod: 'The Protein Quantification Index (PQI) used to convert peptide intensities to protein intensities and estimate the absolute quantities.',
	pqiNormalization: 'How the normalization of PQI values will be performed.<BR>If all the samples are similar, Globally will allow PQI values to be comparable across all samples.<BR>However, it is better to normalize Per Condition if the conditions are quite different from one another.',
	outUnits: 'What kind of values do you want as output ?<BR>Composition only (%) or composition + amounts.',
	absQuantifModel: 'The model to use to convert PQI to absolute quantities.<BR>Choose Direct Proportionality if you do not have any known peptides/proteins to calibrate the model,<BR>otherwise the Calibration Curve should be more accurate.',
	quantityFile:
`A comma-separated values file containing the known quantities in each sample whether it is the total protein quantity or the quantity of your calibration proteins.<BR>
<B>Two columns</B> are always required:<BR>
-'<B>Analysis</B>' with the name of the analyses, which must match exactly those in MyProMS<BR>
-'<B>Quantity</B>' for the known quantities.<BR>
If your samples are labeled (SILAC, TMT, etc.), you must add a '<B>Channel</B>' column and indicate the quantities per analysis/channel pairs.<BR>
Finally, if you use the calibration model, a '<B>Protein</B>' column is also required to indicate the protein names of the calibration proteins. In that case, protein names must match those displayed in myProMS and the quantities must be given by protein, for each analysis or analysis/channel pair.`,
	providedQuantities: 'The unit of measurement in which the quantities of the .csv file are provided. Required to do the right conversions.',
	peptideTopX: 'A positive integer value of the top X peptides to consider to compute the PQI.',
	peptideStrictness: 'Whether we should only consider proteins with the minimal peptide number given by Peptide TopX ("strict") or compute TopX on all proteins even if they have less peptides ("loose").'
};
|;
&promsMod::popupInfo();
###>Templates in JS
print "var quantifTemplate={";
my $first=1;
foreach my $label (keys %quantifTemplate) {
	foreach my $projID (keys %{$quantifTemplate{$label}}) {
		foreach my $quantiID (keys %{$quantifTemplate{$label}{$projID}}) {
			if ($first) {$first=0;} else {print ',';}
			print "\n\t$quantiID:['$quantifTemplate{$label}{$projID}{$quantiID}{'NAME'}','$quantifTemplate{$label}{$projID}{$quantiID}{QUANTIFANNOT}'";
			print ",'$quantifTemplate{$label}{$projID}{$quantiID}{MODID}'" if $quantifTemplate{$label}{$projID}{$quantiID}{'MODID'};
			print "]";
		}
	}
}
print "\n};\n";
###>Unmatched PTMs in JS
print "var unusedVarMods=[];\n";
if (scalar keys %unusedProjVarMods) {
	print "unusedVarMods=[",join(',',keys %unusedProjVarMods),"];\n";
}

print qq
|window.onload=function() { // needed for design quantif because of auto-selection if only 1 labeling available
	updateLabeling(document.selAnaForm.labeling.value);
	updateAlgoTypeSelection(document.selAnaForm.algoType.value);
	// recoverDefaultFormParam(); TEMPORARY DISABLED
};

var defaultAnaFormValues=[];
function recoverDefaultFormParam() {
	var params=document.getElementsByClassName('template');
	var i=0; // important
	for (i=0; i<params.length; i++) {
		if (params[i].type=='select-one') {
			defaultAnaFormValues[i]=[params[i].name,params[i].id,params[i].selectedIndex,params[i].type];
		}
		else if (params[i].type=='checkbox') {
			if (document.getElementsByName(params[i].name).length>1){
				defaultAnaFormValues[i]=[params[i].name,params[i].id,params[i].checked,params[i].type,'',params[i].value];
			}
			else {
				defaultAnaFormValues[i]=[params[i].name,params[i].id,params[i].checked,params[i].type];
			}
		}
		else if (params[i].tagName == 'DIV' \|\| params[i].tagName == 'SPAN'){
			if (params[i].style.visibility){
				defaultAnaFormValues[i]=[params[i].tagName,params[i].id,params[i].style.visibility,'style.visibility'];
			}
			else {
				defaultAnaFormValues[i]=[params[i].tagName,params[i].id,params[i].style.display,'style.display'];
			}
		}
		else {
			defaultAnaFormValues[i]=[params[i].name,params[i].id,params[i].value,params[i].type];
		}
		if (params[i].disabled) {
			defaultAnaFormValues[i][4]=params[i].disabled;
		}
	}
	/*var allDIV=document.getElementsByTagName("div");
	for (let j=0; j<allDIV.length; j++) {
		if (allDIV[j].id) {
			if (allDIV[j].style.visibility) {
				defaultAnaFormValues[i]=['DIV',allDIV[j].id,allDIV[j].style.visibility,'style.visibility'];
			}
			else {
				defaultAnaFormValues[i]=['DIV',allDIV[j].id,allDIV[j].style.display,'style.display'];
			}
			i++;
		}
	}
	var allSPAN=document.getElementsByTagName("span");
	for (let j=0; j<allSPAN.length; j++) {
		if (allSPAN[j].id) {
			if (allSPAN[j].style.visibility) {
				defaultAnaFormValues[i]=['SPAN',allSPAN[j].id,allSPAN[j].style.visibility,'style.visibility'];
			}
			else {
				defaultAnaFormValues[i]=['SPAN',allSPAN[j].id,allSPAN[j].style.display,'style.display'];
			}
			i++;
		}
	}*/
}

function useTemplate(quantiID,skipLabeling) {
	/*** Restore default parameters ***/
	for (let i=0; i<defaultAnaFormValues.length; i++) {
		if (defaultAnaFormValues[i][3] == 'hidden' \|\| defaultAnaFormValues[i][3] == 'button') continue;
		var paramIdent='';
		if (defaultAnaFormValues[i][1] != '') {
			paramIdent=document.getElementById(defaultAnaFormValues[i][1]);
		}
		else {
			paramIdent=document.getElementsByName(defaultAnaFormValues[i][0])[0];
		}

		if (defaultAnaFormValues[i][3] == 'select-one'){
			if (paramIdent.name != 'labeling' \|\| !skipLabeling) {
				paramIdent.selectedIndex=defaultAnaFormValues[i][2];
			}
		}
		else if (defaultAnaFormValues[i][3] == 'checkbox') {
			if (defaultAnaFormValues[i][5]) {
				var selCheckbox=document.getElementsByName(defaultAnaFormValues[i][0]);
				for (let m=0; m<selCheckbox.length; m++) {
					if (selCheckbox[m].value == defaultAnaFormValues[i][5]) {
						selCheckbox[m].checked=defaultAnaFormValues[i][2];
					}
				}
			}
			else {
				paramIdent.checked=defaultAnaFormValues[i][2];
			}
		}
		else if (defaultAnaFormValues[i][3] == 'style.visibility') {
			paramIdent.style.visibility=defaultAnaFormValues[i][2];
		}
		else if (defaultAnaFormValues[i][3] == 'style.display') {
			paramIdent.style.display=defaultAnaFormValues[i][2];
			if (defaultAnaFormValues[i][1] == 'refCommentSPAN') {
				paramIdent.innerHTML='';
			}
		}
		else{
			paramIdent.value=defaultAnaFormValues[i][2];
		}
		if (defaultAnaFormValues[i][4]){
			paramIdent.disabled=defaultAnaFormValues[i][4];
		}
		if (defaultAnaFormValues[i][0] == 'minInfRatios'){
			paramIdent.disabled=false;
		}
	}

	/*** Extract parameters from template ***/
	if (quantiID == 'select') return;
	var quantifParam=quantifTemplate[quantiID][1].split('::');
	var software; // undefined or array
	var peptopN='';
	var singleRef=2;
	var ratio='';
	var pepCharge='';
	var pepSource='';
	var algoType='';
	var topN='';
	var useListBias='';
	var idListBias='';
	var normalizationMethod='';
	var biasCorrection='';
	var ambiguousPos=0;
	var ptmID=(quantifTemplate[quantiID][2])? quantifTemplate[quantiID][2] : '';
	var createSiteRes; // undefined or array
	var badProbAct;

	// Template : select template name
	var selTemplate=document.getElementById('template');
	for(let i=1; i<selTemplate.options.length;i++){
		if (selTemplate.options[i].value == quantiID){
			selTemplate.selectedIndex=i;
		}
	}

	for (let p=0; p<quantifParam.length; p++){
		var param=quantifParam[p].split('=');

		// Strategy
		if (param[0] == 'SOFTWARE') {
			software=param[1].split(';');
		}
		else if (param[0] == 'SINGLE_REF') {
			singleRef=param[1];
		}
		else if (param[0] == 'RATIO_TYPE') {
			ratio=param[1];
		}
		else if(param[0] == 'ALGO_TYPE') {
			algoType=param[1];
		}
		else if (param[0] == 'LABEL') {
			var selLabel=document.getElementsByName('labeling');
			for (let i=0; i<selLabel[0].options.length; i++) {
				if (selLabel[0].options[i].value==param[1]) {
					selLabel[0].selectedIndex=i;
					break;
				}
			}
			updateLabeling(param[1]);
		}
		else if (param[0] == 'TOP_N') {
			topN=param[1];
		}
		else if (param[0] == 'CREATE_PTM') {
			createSiteRes=param[1].split(';');
		}

		// peptide selection
		else if (param[0] == 'PEPTIDES') {
			var peptideInfo=param[1].split(';');
			var selPepSpecif=document.getElementsByName('pepSpecifity');
			for (let i=0; i<selPepSpecif[0].options.length; i++) {
				if (selPepSpecif[0].options[i].value==peptideInfo[0]) {
					selPepSpecif[0].selectedIndex=i;
					break;
				}
			}
			var selPepCleav=document.getElementById('pepCleavage');
			for (let i=0; i<selPepCleav.options.length; i++) {
				if (selPepCleav.options[i].value==peptideInfo[1]) {
					selPepCleav.selectedIndex=i;
					break;
				}
			}

			var selPepPTM=document.getElementsByName('pepPTM');
			var selPTM=peptideInfo[2].split(':');
			for (let i=0; i<selPepPTM[0].options.length; i++) {
				if (selPepPTM[0].options[i].value==selPTM[0]) {
					selPepPTM[0].selectedIndex=i;
					break;
				}
			}
			updateCustomPtmDisplay(selPTM[0]);

			if (selPTM[0]==2) {
				var selCustomPTM=document.getElementsByName('customPTM');
				var mod=selPTM[1].split(/[,#]+/);
				for (let i=0; i<selCustomPTM.length; i++) {
					for (let j=1; j<mod.length ; j++){
						if (selCustomPTM[i].value==mod[j]) {
							selCustomPTM[i].checked=true;
						}
					}
				}
			}
			else if (selPTM[0]==-2) {
				var selCustomPTM=document.getElementsByName('customPTM');
				var mod=selPTM[1].split(/[,#]+/);
				var arrayCheckBox=[];
				for (let i=0; i<selCustomPTM.length; i++){
					arrayCheckBox[i]=selCustomPTM[i].value;
				}
				for (let j=1; j<mod.length ; j++) {
					for (let i=0; i<arrayCheckBox.length ; i++) {
						if (arrayCheckBox[i] == mod[j]) {
							arrayCheckBox.splice(i,1);
						}
					}
				}
				for (let i=0; i<selCustomPTM.length; i++) {
					for (let j=0; j<arrayCheckBox.length ; j++) {
						if (selCustomPTM[i].value==arrayCheckBox[j]) {
							selCustomPTM[i].checked=true;
						}
					}
				}
			}
			if (peptideInfo[3]) {
				pepCharge=peptideInfo[3];
			}
			if (peptideInfo[4]){
				pepSource=peptideInfo[4];
			}
		}
		// Target
		else if (param[0] == 'PTM_POS') {
			if (param[1].match(';')) {
				var ptmPosInfo=param[1].split(';');
				badProbAct=ptmPosInfo[1];
			}
			else if (param[1] == 'ambiguous') {
				ambiguousPos=1;
			}
		}
		// protein selection
		else if (param[0] == 'PROTEINS') {
			var selectType=param[1].split(';');
			if (selectType[0] == 'restrict'){
				document.getElementsByName('protSelType')[0].selectedIndex=1;
			}
			else if (selectType[0] == 'exclude'){
				document.getElementsByName('protSelType')[0].selectedIndex=0;
			}
			var selProtSelect=document.getElementsByName('protSelection');
			for (let i=1; i<selProtSelect[0].options.length; i++) {
				if (selProtSelect[0].options[i].value==selectType[1].split('#')[1]) {
					selProtSelect[0].selectedIndex=i;
					break;
				}
			}
		}
		// Quantification settings
		else if (param[0] == 'BIAS_CORRECTION') {
			if (param[1].match(';')){
				var biasInfo=param[1].split(';');
				biasCorrection=biasInfo[0];
				idListBias=biasInfo[1].substr(1,1);
				useListBias=biasInfo[2];
			}
			else {
				biasCorrection=param[1];
			}
		}
		else if (param[0] == 'NORMALIZATION_METHOD') {
			normalizationMethod=param[1];
		}
		else if (param[0] == 'MIN_INF_RATIOS' && param[1] == 1) {
			document.getElementsByName('minInfRatios')[0].checked=true;
		}
		else if (param[0] == 'FDR_CONTROL'){
			var selFDR=document.getElementsByName('fdrControl');
			for (let i=0; i<selFDR[0].options.length; i++) {
				if (selFDR[0].options[i].value==param[1]) {
					selFDR[0].selectedIndex=i;
					break;
				}
			}
		}
/* Too risky for auto-selected by template (PP)
		else if (param[0] == 'RESIDUAL_VAR'){
			var selResidualVar=document.getElementsByName('residualVar');
			for (let i=0; i<selResidualVar[0].options.length; i++) {
				if (selResidualVar[0].options[i].value==param[1]) {
					selResidualVar[0].selectedIndex=i;
					break;
				}
			}
			document.getElementById('resVarSPAN').style.visibility=(param[1]=='biological')? 'visible' : 'hidden';
		}
*/
	}

	/*** Apply template parameters ***/
	// algorithm selection
	var selAlgo=document.getElementById('algoType');
	if (!algoType){
		algoType=(software && software[0].match('DIA'))? 'MSstats' : 'SSPA';
	}
	var algoTypeStrg='';
	if (algoType == 'SSPA') {
		algoTypeStrg=algoType;
	}
	else {
		var ref=(singleRef == 1) ? 'singleRef' : 'multiRef';
		algoTypeStrg=algoType+':'+ratio+':'+ref;
	}
	for (let i=1; i<selAlgo.options.length; i++) {
		if (selAlgo.options[i].value==algoTypeStrg) {
			selAlgo.selectedIndex=i;
			break;
		}
	}
	updateAlgoTypeSelection(algoTypeStrg);
	if (algoType == 'PEP_INTENSITY' \|\| algoType == 'PEP_RATIO') {
		var selPepCharge=document.getElementsByName('pepCharge');
		for (let i=1; i<selPepCharge[0].options.length; i++) {
			if (selPepCharge[0].options[i].value==pepCharge) {
				selPepCharge[0].selectedIndex=i;
				break;
			}
		}
		updatePeptideSource(pepCharge);
		if (algoType == 'PEP_INTENSITY'){
			topN=(topN)? topN : '0';
			var selTopN=document.getElementsByName('topN');
			for (let i=1; i<selTopN[0].options.length; i++) {
				if (selTopN[0].options[i].value==topN) {
					selTopN[0].selectedIndex=i;
					break;
				}
			}
		}
	}
	else if (algoType == 'SSPA') {
		var selPepFocus=document.getElementsByName('pepFocus');
		for (let i=1; i<selPepFocus[0].options.length; i++) {
			if (selPepFocus[0].options[i].value==pepCharge) {
				selPepFocus[0].selectedIndex=i;
				break;
			}
		}
	}
	// Bias correction
	if (biasCorrection == 'TRUE') {
		var selBiasCorrect=document.getElementById('biasCorrect');
		for (let i=1; i<selBiasCorrect.options.length; i++) {
			if (selBiasCorrect.options[i].value==normalizationMethod) {
				selBiasCorrect.selectedIndex=i;
				break;
			}
		}
		updateBias(normalizationMethod);
		if (idListBias) {
			var selBiasRefProt=document.getElementsByName('refProt');
			for (let i=1; i<selBiasRefProt[0].options.length; i++) {
				if (selBiasRefProt[0].options[i].value==idListBias) {
					selBiasRefProt[0].selectedIndex=i;
					break;
				}
			}
		}
		if (useListBias) {
			document.getElementsByName('refProtSelType')[0].selectedIndex=(useListBias == 'not') ? '1' : '0';
		}
	}

	// Target
	if (ptmID) {
		var selPtmID=(createSiteRes)? -ptmID : ptmID;
		var selPtmScope=document.selAnaForm.ptmScope;
		for (let i=1; i<selPtmScope.options.length; i++) {
			if (selPtmScope.options[i].value==selPtmID) {
				selPtmScope.selectedIndex=i;
				break;
			}
		}
		// updateRefQuantification(selPtmID);
		updateTopN();
		
		if (createSiteRes) {
			// Apply target residue selection
			if (document.getElementById('targetResSPAN_'+ptmID)) {
				var tgtResParam=document.getElementsByName('targetRes_'+ptmID);
				if (tgtResParam.length) {
					for (let i=0; i<tgtResParam.length; i++) {
						for (let j=0; j<createSiteRes.length; j++) {
							if (tgtResParam[i].value==createSiteRes[j]) {tgtResParam[i].checked=true;}
						}
					}
				}
				else {
					for (let j=0; j<createSiteRes.length; j++) {
						if (tgtResParam.value==createSiteRes[j]) {tgtResParam.checked=true;}
					}
				}
				// document.getElementById('recreateSPAN').style.display='';
				document.getElementById('targetResSPAN_'+ptmID).style.display='';
			}
			else {
				alert('WARNING: Template PTM site was not matched.');
			}
		}
		else {
			if (ambiguousPos) {
				document.getElementsByName('ambiguousPos')[0].checked=true;
			}
			else if (badProbAct) {
				document.getElmentsByName('badProb')[0].selectedIndex=(badProbAct == 'ambiguous')? 1 : 0;
			}
		}
	}
}
function clearTemplate() {
	var templateSel=document.getElementById('template');
	if (templateSel.value != 'select') {
		templateSel.selectedIndex=0;
		useTemplate('select',true);
	}
}
function updateQuantifScope(scope) {
	updateTopN();
	const myForm=document.selAnaForm;
	if (scope.match('protein')) { // protein,proteinPTM
		// Clear PTM filter spans
		document.getElementById('seqContextSPAN').style.display='none';
		// var ptmProbSpan=document.getElementById('ptmProbSPAN'); // optional
		// if (ptmProbSpan) ptmProbSpan.style.display='none';
		// Clear reference quantif
		document.getElementById('refQuantifDIV').style.display='none';
		myForm.refQuantifType.selectedIndex=0;
		ajaxFetchReferenceQuantif(0);
		// Clear any free residue selection
		document.getElementById('createPtmDIV').style.display='none';
		// Update LFQ min # ratio
		myForm.minRatioLFQ.value=2;
		// Disable manual peptide selection for Bias correction
		if (myForm.refProt.value=='manual') { // manual peptide => reset to "All proteins"
			myForm.refProt.options[1].disabled=true;
			updateBiasCorrectionSet('all');
		}
		if (scope==='protein') {
			// Clear any multi-PTMs selection
			document.getElementById('multiPtmDIV').style.display='none';
/*
			var multiPtmChkBoxes=document.selAnaForm.multiPtmScope;
			if (multiPtmChkBoxes.length) {
				for (let i=0; i<multiPtmChkBoxes.length; i++) {
					multiPtmChkBoxes[i].checked=false;
				}
			}
			else {
				multiPtmChkBoxes.checked=false;
			}
*/
		}
	}
	if (scope.match('PTM')) { // multiPTM,create_multiPTM,proteinPTM
		if (scope==='create_multiPTM') {
			document.getElementById('multiPtmDIV').style.display='none';
			document.getElementById('createPtmDIV').style.display='';
		}
		else { // multiPTM,proteinPTM
			document.getElementById('multiPtmDIV').style.display='';
			document.getElementById('createPtmDIV').style.display='none';
			// updateMultiPtmSelection();
		}
		if (scope.match('multiPTM')) { // multiPTM,create_multiPTM
			document.getElementById('seqContextSPAN').style.display='';
			document.getElementById('refQuantifDIV').style.display=(document.getElementById('algoType').value.match('^(MSstats\|SSPA\|Abundance)'))? 'none' : ''; // TDA?
			// Update LFQ min # ratio
			myForm.minRatioLFQ.value=1;
			// Enable manual peptide selection for Bias correction
			myForm.refProt.options[1].disabled=(myForm.algoType.value.match('(PEP_INTENSITY\|DIA\|TDA):SimpleRatio') && myForm.refQuantifType.value*1 == 0 && $numProtInDesign <= $MAX_PROT_MANUAL_PEP)? false: true;
		}
	}
	
	// PTM position probability, context motif, other allowed PTMs, manual reference quantif
	updateMultiPtmSelection();

	// MG shared peptides
	updateMGSharedPeptides(undefined,undefined,scope);

	// Bias correction options text
	[myForm.biasCorrLevel.options[1].text,myForm.biasCorrLevel.options[2].text]=(scope.match('protein'))? ['Protein','Peptide & protein'] : ['Site','Peptide & site'];
}
function updateMultiPtmSelection() {
	// Fetch list of selected PTM-sites
	var selPtmList={},
		numSelMatched=0,
		multiPTMs=document.selAnaForm.multiPtmScope,
		ptmScope=document.selAnaForm.ptmScope.value;
	if (!multiPTMs) {return;} // no project PTMs declared
	if (multiPTMs.length) {
		for (let i=0; i<multiPTMs.length; i++) {
			var ptmID=Math.abs(multiPTMs[i].value*1);
			if (multiPTMs[i].checked) {
				selPtmList[ptmID]=multiPTMs[i].value;
				numSelMatched++;
			}
			document.getElementById('contextSPAN_'+ptmID).style.display=(multiPTMs[i].checked && ptmScope==='proteinPTM')? '' : 'none'; // Context only for proteinPTM (for now)
		}
	}
	else {
		var ptmID=Math.abs(multiPTMs.value*1);
		if (multiPTMs.checked) {
			selPtmList[ptmID]=multiPTMs.value;
			numSelMatched=1;
		}
		document.getElementById('contextSPAN_'+ptmID).style.display=(multiPTMs.checked && ptmScope==='proteinPTM')? '' : 'none'; // Context only for proteinPTM (for now)
	}
	// Re-define custom PTMs disable status
	var cusPTMs=document.selAnaForm.customPTM;
	if (cusPTMs.length) {
		for (let i=0; i<cusPTMs.length; i++) {
			cusPTMs[i].disabled=((ptmScope==='proteinPTM' \|\| ptmScope==='multiPTM') && selPtmList[cusPTMs[i].value])? true : false;
		}
	}
	else {cusPTMs.disabled=((ptmScope==='proteinPTM' \|\| ptmScope==='multiPTM') && selPtmList[cusPTMs.value])? true : false;}

	// Phospho (or any PTM if MaxQuant) prob threshold
	var ptmProbSpan=document.getElementById('ptmProbSPAN'),
		// ambigChk=document.selAnaForm.ambiguousPos;
		badProb=document.selAnaForm.badProb,
	    phosphoChecked=(selPtmList['$phosphoID'])? true : false;
	if (ptmProbSpan) {
		if (ptmScope==='multiPTM' && ((phosphoChecked && '$ptmProbSoftCode') \|\| $isMultiPtmProb) && !document.selAnaForm.algoType.value.match('^TDA') ) { // pos prob filtering not yet handled in TDA
			ptmProbSpan.style.display='';
			if (badProb) { // MaxQuant or Phospho: Ambiguity allowed ONLY for single-PTM quantif
				if (numSelMatched > 1) {
					badProb.selectedIndex=0;
					badProb.options[1].disabled=true;
				}
				else {badProb.options[1].disabled=false;}
			}
		}
		else {ptmProbSpan.style.display='none';}
	}
	
	// Update manual peptide selection (if displayed)
	if (ptmScope.match('multiPTM')) {
		if (document.selAnaForm.refQuantifType.value==-1 && document.selAnaForm.manualRefPep.value=='manual') { // if (document.getElementById('manualRefPepDIV').style.display != 'none')
			ajaxFetchRefProtPeptides('manual','protCorr');
		}
		else if (document.selAnaForm.refProt.value=='manual') { // if (document.getElementById('manualBiasPepDIV').style.display != 'none')
			ajaxFetchRefProtPeptides('manual','biasCorr');
		}
	}
}
function updateTopN() {
	const myForm=document.selAnaForm;
	const algoType=myForm.algoType.value;
	if (algoType.match('PEP_INTENSITY\|DIA\|TDA') && myForm.ptmScope.value.match('protein')) {
		myForm.topN.disabled=false;
		myForm.matchingPep.disabled=false;
		document.getElementById('topSPAN').style.display='';
	}
	else {
		myForm.topN.disabled=true;
		myForm.matchingPep.disabled=true;
		document.getElementById('topSPAN').style.display='none';
	}
	document.getElementById('trackTrueSPAN').style.display=(algoType.match('PEP_INTENSITY'))? 'block' : 'none'; // :SimpleRatio block for newline w/o <BR>
}
function updateMultiQuantifOption(algoType) {
	var myForm=document.selAnaForm;
	if (!algoType) algoType=myForm.algoType.value;
	var numSelStates=0;
	for (let i=1; i <= $nbCond; i++) {
		if (!myForm['state_'+i].disabled && myForm['state_'+i].value != 'Select') numSelStates++;
	}
	var [multiQVis,multiQDisab]=(numSelStates >= 3 && algoType.match(':SimpleRatio:'))? ['',false] : ['none',true]; // :singleRef
	document.getElementById('multiQuantifSPAN').style.display=multiQVis;
	var chkMultiQ=document.getElementById('multiQuantifCHK');
	chkMultiQ.disabled=multiQDisab;
	
	/* Protein FC correction for Sites */
	if (myForm.ptmScope.value==='multiPTM') {
		document.getElementById('refQuantifDIV').style.display=(algoType.match('^(MSstats\|DIA\|TDA\|Abundance\|SSPA)'))? 'none' : '';
	}
}
function ajaxFetchReferenceQuantif(designID) {
	var myForm=document.selAnaForm;
	if (myForm.labeling.value=="") {
		alert('Select a labeling strategy first!');
		myForm.refQuantifType.selectedIndex=0;
		return;
	}

	// Reset prot-level corr. manual peptide selection
	myForm.manualRefPep.selectedIndex=0;
	document.getElementById('manualRefPepDIV').innerHTML='';
	document.getElementById('manualRefPepDIV').style.display='none';
	document.getElementById('manualRefPepSPAN').style.display=(designID >= 0)? 'none' : ''; 
	var algoSel=document.getElementById('algoType');
	// Disable 3+ states for PEP_RATIO (until handled by R scripts)
	//if (algoSel.value.match('PEP_RATIO:')) {
	//	var disab3States=(designID==0)? false : true;
	//	for (let i=3; i<=$nbCond; i++) {myForm['state_'+i].disabled=disab3States;}
	//}
	// Disable/Enable 'All vs All' algos
	var selectedIdx=algoSel.options.selectedIndex;
	for (let i=1; i<algoSel.options.length; i++) { // skip -= Select =-
		if (algoSel.options[i].value.match('MSstats:')) continue;
		if (designID==0) { // enable option
			if (algoSel.options[i].value=='PEP_INTENSITY:SimpleRatio:multiRef') {algoSel.options[i].disabled=false;}
			else if (algoSel.options[i].value=='PEP_RATIO:SuperRatio:multiRef' && ($quantifMethods{SILAC} \|\| $quantifMethods{ITRAQ} \|\| $quantifMethods{TMT})) {algoSel.options[i].disabled=false;}
			else if (algoSel.options[i].value=='TDA:SimpleRatio:multiRef' && $quantifMethods{TDA}) {algoSel.options[i].disabled=false;}
		}
		else if (algoSel.options[i].value.match(':multiRef')) { // disable option
			algoSel.options[i].disabled=true;
		}
	}
	if (designID != 0 && algoSel.options[selectedIdx].disabled) { // selected option was disabled => use 1 just before
		var ratioStrg=(algoSel.options[selectedIdx].value=='PEP_RATIO:SuperRatio:multiRef')? '"Super ratios"' : '"All vs All" ratios';
		alert('WARNING: '+ratioStrg+' are not possible with "Protein-level fold change correction". Switching to "All vs State1"');
		algoSel.options.selectedIndex=selectedIdx-1;
		updateAlgoTypeSelection(algoSel.value);
	}

	// Disable/Enable manual pept. selection for Exp. bias correction
	if (designID != 0) {
		myForm.refProt.options[1].disabled=true; // 'manual'
		updateBiasCorrectionSet('all','protCorr');
	}
	else {
		updateBias(myForm.biasCorrect.value);
	}
	
	// Hide refQuantif selection
	var refQuantSpan=document.getElementById('refQuantifSPAN');
	if (designID <= 0) {
		refQuantSpan.innerHTML='<INPUT type="hidden" name="refQuantif" value="-1"/>';
		return;
	}
		
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/startDesignQuantification.cgi?AJAX=fetchRefQuantif&ID="+designID+"&labeling="+document.selAnaForm.labeling.value+"&algoType="+document.selAnaForm.algoType.value,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			refQuantSpan.innerHTML='<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Select a reference dataset from:'+XHR.responseText+' <FONT style="color:#DD0000">(WARNING: States order must match with current dataset)</FONT>';
		}
	}
	XHR.send(null);
}
function applyNewRefQuantifSet(quantifSet) {
	if (document.selAnaForm.manualRefPep.selectedIndex==0) {return;} // "Current rules"
	if (!quantifSet) { // -= Select =-
		document.selAnaForm.manualRefPep.selectedIndex=0;
	}
	ajaxFetchRefProtPeptides(document.selAnaForm.manualRefPep.value);
}
function ajaxFetchRefProtPeptides(selVal,call='protCorr') { // call=protCorr or biasCorr
	var manualRefPepDiv=(call=='protCorr')? document.getElementById('manualRefPepDIV') : document.getElementById('manualBiasPepDIV');
	var otherRefPepDiv=(call=='protCorr')? document.getElementById('manualBiasPepDIV') : document.getElementById('manualRefPepDIV');
	otherRefPepDiv.style.display='none';
	otherRefPepDiv.innerHTML='';
	if (selVal != 'manual') { // selVal=='0' \|\| (call=='biasCorr' && selVal > 0)
		manualRefPepDiv.style.display='none';
		manualRefPepDiv.innerHTML='';
		document.selAnaForm.normSharedPep.disabled=false;
		return;
	}
	// Check reference dataset
	if (call=='protCorr' && !document.selAnaForm.refQuantif.value) { // "-= Select =-" <- A design was selected but not a quantif ratio
		alert('Select a reference dataset first!');
		document.selAnaForm.manualRefPep.selectedIndex=0;
		return;
	}
	// Check states (at least 2 must be selected)
	document.selAnaForm.normSharedPep.disabled=true;
	var anaList=[];
	if (document.selAnaForm.refQuantif.value==-1) { // current dataset
		var numSelStates=0;
		for (let sPos=1; sPos <= $nbCond; sPos++) {
			var okState=0;
			if (document.selAnaForm['state_'+sPos].value != 'Select' && !document.selAnaForm['state_'+sPos].disabled) {
				var condID=document.selAnaForm['state_'+sPos].value;
				var pepQuant=document.selAnaForm['condObs:'+condID];
				if (pepQuant.length) {
					for (let i=0; i < pepQuant.length; i++) {
						if (pepQuant[i].checked && document.selAnaForm['anaXICQuanti:'+pepQuant[i].value].value) {
							var [obsID,pepQuantifID,anaID,targetPos,fracGroup]=document.selAnaForm['anaXICQuanti:'+pepQuant[i].value].value.split(':');
							anaList.push(anaID);
							okState=1;
						}
					}
				}
				else {
					if (pepQuant.checked && document.selAnaForm['anaXICQuanti:'+pepQuant.value].value) {
						var [obsID,pepQuantifID,anaID,targetPos,fracGroup]=document.selAnaForm['anaXICQuanti:'+pepQuant.value].value.split(':');
						anaList.push(anaID);
						okState=1;
					}
				}
				if (okState==1) {numSelStates++;} // make sure 1+ pep quantif(s) selected 
			}
		}
		if (numSelStates < 2) {
			alert('Select the States to be compared first!');
			document.selAnaForm.manualRefPep.selectedIndex=0; // protCorr
			document.selAnaForm.refProt.selectedIndex=0; // biasCorr
			return;
		}
	}
	//Check quantified modifs//
	var quantifModifList=[];
	if (document.selAnaForm.ptmScope.value==='multiPTM') {
		//var multiPtmChkBoxes=document.selAnaForm.multiPtmScope;
		//for (let i=0; i<multiPtmChkBoxes.length; i++) {
		//	if (multiPtmChkBoxes[i].checked) quantifModifList.push(multiPtmChkBoxes[i].value);
		//}
		var multiScope=document.selAnaForm.multiPtmScope;
		if (!multiScope.length) {multiScope=[multiScope.value];} // force to array if single element
		multiScope.forEach(function(chkBox) {
			if (chkBox.checked) {quantifModifList.push(chkBox.value);}
		});
	}
	else {quantifModifList=(call=='protCorr')? [document.selAnaForm.refQuantif.value] : [-1];}
/* Test below has an hugly side effect after full multi-modif quantif deselection
	if (call=='protCorr' && quantifModifList.length==0) {
		alert('Select PTM(s) to be quantified!');
		document.selAnaForm.manualRefPep.selectedIndex=0;
		return;
	}
*/
	
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/startDesignQuantification.cgi?AJAX=fetchRefProtPep&ID=$designID&refQuantif="+document.selAnaForm.refQuantif.value+"&modID="+(quantifModifList.join(':'))+"&labeling="+document.selAnaForm.labeling.value+"&selAna="+anaList.join(':'),true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			manualRefPepDiv.innerHTML=XHR.responseText;
			manualRefPepDiv.style.display='';
		}
	}
	XHR.send(null);
}
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
// Remember what was selected before for each state
$varGlobSelCond
function updateMGSharedPeptides(pepSpecif=document.getElementById('pepSpecifity').value,algoType=document.getElementById('algoType').value,ptmScope=document.getElementById('ptmScope').value) {
	const [dispStatus,disabStatus1,disabStatus2]=(pepSpecif === 'unique' \|\| !ptmScope.match('protein'))? ['none',true,true] : (algoType.match('Abundance'))? ['',false,false] : ['',true,false];
	document.getElementById('MGSharedPepDIV').style.display=dispStatus;
	const sharedPepSel=document.getElementById('MGSharedPep');
	sharedPepSel.options[1].disabled=disabStatus1; // Share with all
	sharedPepSel.options[2].disabled=disabStatus2; // Exclude
	if (sharedPepSel.options[sharedPepSel.selectedIndex].disabled) sharedPepSel.selectedIndex=0; // falls back to dafault (Best MG)
}
function updateCustomPtmDisplay(filter) {
	document.getElementById('customPtmDIV').style.display=(Math.abs(filter) <= 1)? 'none' : '';
}
function updatePeptideSource(selCharge) {
	var pepSrc=document.getElementById('pepSource');
	if (!pepSrc) return;
	if (selCharge=='best') {
		pepSrc.selectedIndex=1; // best
		pepSrc.disabled=true;
	}
	else {
		pepSrc.selectedIndex=0; // all
		pepSrc.disabled=false;
	}
}
function updateBias(bias) {
	const myForm=document.selAnaForm;
	const algoType=myForm.algoType.value;
	if (algoType.match('MSstats')) {
		document.getElementById('biasProtDIV').style.visibility=(bias && bias.match('REF_PROT\|globalStandards'))? 'visible' : 'hidden';
		myForm.normSharedPep.disabled=true; // shared peptides
		updateBiasCorrectionSet('all');		
	}
	else {
		if (bias && !bias.match('^none\|quantile')) {
			document.getElementById('biasProtDIV').style.visibility='visible';
			myForm.normSharedPep.disabled=false; // shared peptides
			// Manual pept. for bias correction (only possible for LF SimpleRatio in myProMS)
			if (myForm.ptmScope.value==='multiPTM' && algoType.match('(PEP_INTENSITY\|DIA\|TDA):SimpleRatio') && myForm.refQuantifType.value*1 == 0) {
				myForm.refProt.options[1].disabled=($numProtInDesign <= $MAX_PROT_MANUAL_PEP)? false : true; // Incompatible with protein-level FC correction
			}
			else {
				if (myForm.refProt.value=='manual') {
					myForm.refProt.selectedIndex=0;
				}
				myForm.refProt.options[1].disabled=true;
			}
//console.log('updateBias',bias,myForm.refProt.options[1].disabled);//
		}
		else { // *Select*/none/quantile
			document.getElementById('biasProtDIV').style.visibility='hidden';
			myForm.normSharedPep.disabled=true; // shared peptides
			updateBiasCorrectionSet('all');
		}
	}
	document.getElementById('biasCorrLevelSPAN').style.display=(bias && bias !== 'none.none' && algoType.match('Abundance'))? '' : 'none'; // Update bias correction level
}
function updateBiasCorrectionSet(setCode,call='biasCorr') { // 'all':all prot / 'shared':shared proteins / 'manual':manual pep selection / <listID>:custom list  (BEFORE 0:all prot / -1:manual pep selection / <listID>:custom list)
	// "(Do not) Use" switch
	var selType=document.selAnaForm.refProtSelType;
	if (!setCode.match(/^\\d/)) {
		if (setCode != document.selAnaForm.refProt.value) { // updateBiasCorrectionSet also called from other JS functions
			document.selAnaForm.refProt.selectedIndex=0; // can only be 'all' (not 'shared*' or 'manual')
		}
		selType.selectedIndex=0;
		selType.options[1].disabled=(document.selAnaForm.algoType.value.match('(PEP_INTENSITY\|DIA\|TDA):SimpleRatio') && $numProtInDesign <= $MAX_PROT_MANUAL_PEP)? false: true;
	}
	else { // Custom list
		selType.options[1].disabled=false;
	}
//console.log('updateBiasCorrectionSet',setCode,document.selAnaForm.refProt.options[1].disabled);//
	if (call=='protCorr' \|\| document.selAnaForm.refQuantifType.value == 0) { // Do not allow if protein-level correction
		ajaxFetchRefProtPeptides(setCode,'biasCorr'); // Manual peptide selection for bias correction
	}
}
function updateBiasCorrectionLevel() { // called when an abundance measure is (un)checked to update available correction levels
	const measList=document.selAnaForm.abundMeasures,
		biasCorrLevelSel=document.selAnaForm.biasCorrLevel;
	let disabProtLevel=true;
	for (let i=0; i<measList.length; i++) {
		if (measList[i].value==='ABSOLUTE') continue; // absolute quant is not considered
		if (measList[i].checked) {
			disabProtLevel=false;
			break;
		}
	}
	biasCorrLevelSel.options[1].disabled=disabProtLevel;
	biasCorrLevelSel.options[2].disabled=disabProtLevel;
	if (disabProtLevel) {biasCorrLevelSel.selectedIndex=0;}
}
function autoSelectStates(selMode) {
	var myForm=document.selAnaForm;
	var selectedStates=[];
	if (selMode && selMode=='range') {
		var rangeStrg=myForm.stateRange.value;
		var tmpStrg=rangeStrg.replace(/[\\d,-]+/g,'');
		if (tmpStrg.length) {
			alert("Error in range definition (only digits,',' and '-' allowed)!");
			return;
		}
		var selectedStates=[], rangePieces=rangeStrg.split(',');
		for (let i=0; i < rangePieces.length; i++) {
			let [start,end]=rangePieces[i].split('-');
			if (end===undefined) {if (start*1 <= $nbCond) {selectedStates.push(start*1);}}
			else {
				for (let j=start*1; j <= end*1; j++) {
					if (j > $nbCond) {break;}
					selectedStates.push(j);
				}
			}
		}
	}
	else { // use all states
		for (let c=1; c <= $nbCond; c++) {selectedStates.push(c);} // skip 0 -> '-= Select =-'
	}

	var nextUnusedIndex=0;
	for (let c=1; c <= $nbCond; c++) {
		let selCond=myForm['state_'+c];
		if (selCond.selectedIndex > 0) {continue;} // manually selected: do not overwrite
		for (let i=nextUnusedIndex; i < selectedStates.length; i++) {
			if (!selCond.options[selectedStates[i]].disabled) { // not already used by another State
				selCond.selectedIndex=selectedStates[i];
				updateStateSelection(selCond,c);
				nextUnusedIndex=i+1;
				break;
			}
		}
	}
}
function clearAllStates() {
	var myForm=document.selAnaForm;
	for (let c=1; c <= $nbCond; c++) {
		let selCond=myForm['state_'+c];
		selCond.selectedIndex=0;
		updateStateSelection(selCond,c);
	}
}
function updateStateSelection(selCond,pos) { // pos=State pos
	var stateSelected=pos-1;
	var dispDiv=document.getElementById('fillDiv');
	var selCondIdx=selCond.selectedIndex;
	var oldCondIdx=stateCondVal[stateSelected]; // !=0 if another Cond was previously selected for this State
	// Update selectable Conditions
	for (let i = 1; i <= $nbCond; i++) {
		if (i == pos) continue;
		var stateSEL=document.selAnaForm['state_'+i];
		if (selCondIdx != 0) {stateSEL.options[selCondIdx].disabled=true;} // Remove this newly selected Cond from all other States options
		if (oldCondIdx != 0) {stateSEL.options[oldCondIdx].disabled=false;} // new Cond/Select replaces an old one => Free old Cond for all other States options
	}
	// Update Select of Condition display
	dispDiv.options[selCondIdx].disabled=false; // enable corresponding option
	dispDiv.selectedIndex=selCondIdx; // select corresponding option
	if (oldCondIdx != 0) { // new Cond/Select replaces an old one => disable display of old Cond
		dispDiv.options[oldCondIdx].disabled=true;
	}
	// Update stateCondVal
	stateCondVal[stateSelected]=selCondIdx;
	showDivCond(selCond.value); //condID
	// Update multi-quantif option
	//updateMultiQuantifOption(); //????????????????
}
function showDivCond(condID){
	var stateDiv=document.getElementsByName('condDiv');
	if (condID != 'Select') {
		for (let i=0; i < stateDiv.length; i++) {
			if (stateDiv[i].id == condID) {
				stateDiv[i].style.display='inline-block';
			}
			else {
				stateDiv[i].style.display='none';
			}
		}
	}
	else {
		for (let i=0; i < stateDiv.length; i++) {
			stateDiv[i].style.display='none';
		}
	}
}
function updateObsCheckStatus(chkSetName,myChkbox) {
	var chkboxSet=document.getElementsByName(chkSetName);
	var reachedChk=false;
	for (var i=0; i<chkboxSet.length; i++) {
		if (reachedChk) {
			chkboxSet[i].checked=myChkbox.checked;
		}
		else if (chkboxSet[i].value==myChkbox.value) {
			reachedChk=true;
		}
	}
}
function synchronizePropagation(propagSel) {
	var allPropagSels=document.getElementsByName("propagQuantScope");
	for (let i=0; i<allPropagSels.length; i++) {
		if (allPropagSels[i].id != propagSel.id) {allPropagSels[i].selectedIndex = propagSel.selectedIndex;}
	}
}
/* quantiXicSoftware: To allow quanti selection propagation if more than 1 XIC algo */
$jsXicSoftStrg
function propagateQuantifSelection(selQuanti,force) {
	if (!selQuanti.value) return; // -= Select =-
	var valInfo=selQuanti.value.split(':'), // obsID:quantiID:...
		quantifName=selQuanti.options[selQuanti.selectedIndex].text,
		selState=selQuanti.dataset.state,
		xicSoft=quantiXicSoftware[valInfo[1]],
		propagScope=(force)? 'all' : document.getElementById('propagQuantScope:'+selState).value; // none,state*,all
	if (propagScope=='none') {return;}
	if (propagScope=='all' && !force && !confirm('WARNING: This XIC quantification will be propagated to ALL other state(s) where available. Proceed?')) {return;}
	var allSelects=document.getElementsByTagName("select");
	for (let i=0; i<allSelects.length; i++) {
		if (allSelects[i].name.match('anaXICQuanti:') && allSelects[i].name != selQuanti.name) {

			if (propagScope=='state' && allSelects[i].dataset.state != selState) continue;
			//propagScope='all' --->

			let matchLevel={ID:null,NAME:null,SOFT:[]};
			for (let j=1; j<allSelects[i].options.length; j++) { // skip Select
				var idInfo=allSelects[i].options[j].id.split(':'); // opt:<obsID>:<quantiID>
				//Quantif id comparison
				if (idInfo[2]==valInfo[1]) {
					matchLevel.ID=j;
					break; // no need to check other match options
				}
				//XIC Software comparison
				if (quantiXicSoftware[idInfo[2]]==xicSoft) {
					matchLevel.SOFT.push(j);
					//Quantif name comparison (only if XIC soft match)
					if (allSelects[i].options[j].text==quantifName) {
						matchLevel.NAME=j;
						break; // XIC soft & name match
					}
				}
			}
			if (matchLevel.ID) {allSelects[i].selectedIndex=matchLevel.ID;}
			else if (matchLevel.NAME) {allSelects[i].selectedIndex=matchLevel.NAME;}
			else if (matchLevel.SOFT[0]) {allSelects[i].selectedIndex=matchLevel.SOFT[0];}
			else {allSelects[i].selectedIndex=0;} // insure coherent dataset
		}
	}
}
|;

print "var labelTypes=['",join("','",keys %{$designInfo{'LABEL'}}),"'];\n";
print qq
|var algoNorms= {
	PEP_RATIO:[
			['$normalizationNames{"median.none"}','median.none'],
			['$normalizationNames{"median.scale"}*','median.scale'],
			['$normalizationNames{"none.scale"}','none.scale'],
			['$normalizationNames{"quantile"}','quantile']
	],
	PEP_INTENSITY:[
			['$normalizationNames{"median.none"}','median.none'],
			['$normalizationNames{"median.scale"}*','median.scale'],
			['$normalizationNames{"none.scale"}','none.scale'],
			['$normalizationNames{"quantile"}','quantile']
	],
	MSstats:[
			['Equalize medians*','equalizeMedians'],
			['Quantile','quantile'],
			['Global standards','globalStandards']
	]
};
algoNorms.TDA=algoNorms.PEP_INTENSITY; //myProMS algo
algoNorms.DIA=algoNorms.PEP_INTENSITY; //myProMS algo

function updateAlgoTypeSelection(algoType) {
	const myForm=document.selAnaForm;

	/* Update quantif name %...% span & States display */
	let nameSpanDisp='none',
	    referenceText='';
	if (algoType.match('S.+Ratio')) { // myProMS v2+ or MSstats
		nameSpanDisp='';
		referenceText=($nbCond==2)? '&nbsp;&lt;&lt;&lt;&nbsp;Reference state' : (algoType.match(':multiRef'))? '&nbsp;Each state will be used as reference for all following states' : '&nbsp;&lt;&lt;&lt;&nbsp;Common reference';
	}
	//else if (algoType != 'select') { // old algos
	//	referenceText=(algoType=='SSPA')? '&nbsp;: At least 3 replicates required per state' : '&nbsp;&lt;&lt;&lt;&nbsp;Reference state'; //'&nbsp;Each state will be used as reference for all following states'; //
	//}
	document.getElementById('nameSPAN').style.display=nameSpanDisp;
	document.getElementById('refCommentSPAN').innerHTML=referenceText;
	
	/* Update TopN */
	updateTopN();

	/* Multi-quantif option */
	updateMultiQuantifOption(algoType);
	
	/* Update peptide selection & PTM quantif options */
	let dispRatioPepDiv='none',dispSspaPepDiv='none',disabAnyPep=true;
	if (algoType=='SSPA') {
		dispSspaPepDiv='';
		myForm.pepRescued.disabled=(myForm.pepFocus.value==='sp_count')? true : none;
/* SSPA for SITE enabled PP 05/10/20
		if (myForm.ptmScope.selectedIndex > 1) { // Whole protein (by PTM)
			myForm.ptmScope.selectedIndex=0;
			updateQuantifScope(myForm.ptmScope.value);
		}
		for (let i=2; i<myForm.ptmScope.options.length; i++) { // disable PTM quantif
			myForm.ptmScope.options[i].disabled=true;
		}
*/
	}
	else {
		if (!algoType.match('select\|MSstats')) {
			dispRatioPepDiv='';
			if (algoType.match('Abundance\|PEP_INTENSITY\|DIA\|TDA')) {disabAnyPep=false;}
		}
		for (let i=1; i<myForm.ptmScope.options.length; i++) { // enable (if exist PTM) PTM-based quantif
			myForm.ptmScope.options[i].disabled=($numPTMs)? false : true;
		}
	}
	document.getElementById('ratioPepDIV').style.display=dispRatioPepDiv;
	document.getElementById('sspaPepSPAN').style.display=dispSspaPepDiv;
	document.getElementById('rescuedPepDIV').style.display=dispSspaPepDiv;
	document.getElementById('hiddenProtSPAN').style.display=dispSspaPepDiv;
	//if (disabAnyPep) {myForm.matchingPep.selectedIndex=1;} // matching pep
	updateMGSharedPeptides(undefined,algoType);

	if (algoType.match('Abundance')) {
		myForm.matchingPep.selectedIndex=0; // "Any"
		myForm.topN.selectedIndex=5; // "N"
		document.getElementById('abundMeasuresSPAN').style.display='block'; // block for newline w/o <BR>
	}
	else {
		myForm.matchingPep.selectedIndex=1; // "matching"
		myForm.topN.selectedIndex=2; // "3"
		document.getElementById('abundMeasuresSPAN').style.display='none';
	}
	myForm.matchingPep.options[0].disabled=disabAnyPep;
//updateRefQuantification(myForm.ptmScope.value); // because of SSPA & MSstats constraints
//updateMultiPtmSelection(); // because of SSPA & MSstats constraints

	/* Update bias correction options */
	var normSelOpt=myForm.biasCorrect.options;
	normSelOpt.length=2;
	//var trueAlgo=(algoType.match('PEP_'))? 'SxxxRatio' : (algoType.match('MSstats'))? 'MSstats' : (algoType.match(':Ratio'))? 'RatioTnPQ' : algoType;
	var [trueAlgo,ratioType,ratioPairs]=algoType.split(':');
	if (algoNorms[trueAlgo]) {
		myForm.biasCorrect.disabled=false;
		//var labelType=myForm.labeling.value;
		//for (var i=0; i<algoNorms[trueAlgo][labelType].length; i++) { //}
		//	normSelOpt[i+2]=new Option(algoNorms[trueAlgo][labelType][i][0],algoNorms[trueAlgo][labelType][i][1]);
		//}
		for (let i=0; i<algoNorms[trueAlgo].length; i++) { //}
			normSelOpt[i+2]=new Option(algoNorms[trueAlgo][i][0],algoNorms[trueAlgo][i][1]);
		}
	}
	else { // SSPA, select
		myForm.biasCorrect.disabled=true;
	}
	updateBias();

	/* Update remaining quantification options */
	document.getElementById('biasProtDIV').style.visibility='hidden'; // hide norm prot list in case visible
	var disabResVar=true, visResVarSpan='hidden';
	if (algoType.match('PEP_\|TDA\|DIA') && !algoType.match('Abundance')) {
		disabResVar=false;
		if (myForm.residualVar.value=='biological') {visResVarSpan='visible';}
	}
	myForm.residualVar.disabled=disabResVar;
	document.getElementById('resVarSPAN').style.visibility=visResVarSpan;
	myForm.fdrControl.disabled=(algoType.match('select\|Abundance'))? true : false;
	if (algoType==='SSPA') { // force to Benjamini-Hochberg (fdr) as only option
		for (let i=0; i < myForm.fdrControl.options.length; i++) {
			if (myForm.fdrControl.options[i].value=='fdr') {myForm.fdrControl.selectedIndex=i;}
			else {myForm.fdrControl.options[i].disabled=true;}
		}
	}
	else {
		for (let i=0; i < myForm.fdrControl.options.length; i++) {myForm.fdrControl.options[i].disabled=false;}
	}
	[document.getElementById('minInfRatiosDIV').style.display,myForm.minInfRatios.disabled]=(algoType.match('Ratio') && !algoType.match('MSstats'))? ['',false] : ['none',true];
	
	[document.getElementById('sspaSettingsDIV').style.display,document.getElementById('noSspaSettingsDIV').style.display]=(algoType==='SSPA')? ['','none'] : ['none',''];

	/* Check for matching quantif type data (if multiple types available) */
	var compatibleMethods={};
	if (algoType.match('MSstats\|DIA')) {compatibleMethods.DIA=1;}
	else if (algoType.match('TDA')) {compatibleMethods.TDA=1;}
	else if (algoType.match('PEP_RATIO')) {
		compatibleMethods={
			SILAC: 1,
			ITRAQ: 1,
			TMT: 1
		};
	}
	else if (algoType.match('PEP_INTENSITY')) {
		compatibleMethods={
			XIC: 1,
			SILAC: 1,
			ITRAQ: 1,
			TMT: 1
		};
	}
	else { // SSPA or Abundance
		compatibleMethods={
			XIC: 1,
			SILAC: 1,
			ITRAQ: 1,
			TMT: 1,
			DIA: 1
		};
	}
	var labelType=myForm.labeling.value;
	var optList=(labelType)? document.getElementsByName(labelType) : [];
	var firstGoodOpt, propagate=false;
	for (let i=0; i<optList.length; i++) {
		//if (optList[i].dataset.quantif_method==goodMethod) { //}
		if (compatibleMethods[optList[i].dataset.quantif_method]) {
			if (!firstGoodOpt) firstGoodOpt=optList[i];
			optList[i].disabled=false;
		}
		else {
			if (optList[i].selected) propagate=true; // at least 1 selected option was disabled
			optList[i].disabled=true;
		}
	}
	if (propagate && firstGoodOpt) {
		var idInfo=firstGoodOpt.id.split(':');
		var selQuanti=document.getElementById('anaXICQuanti:'+idInfo[1]);
		firstGoodOpt.selected=true;
		propagateQuantifSelection(selQuanti,true); // true -> propagate selection to all
	}
}
function updateLabeling(labelType) { // FREE, SILAC, ITRAQ, TMT
	/* Update list of compatible algos */
	var algoSelOpt=document.getElementById('algoType').options;
	algoSelOpt[0].selected=true; // -= Select =-
	for (var i=1; i<algoSelOpt.length; i++) {
		algoSelOpt[i].disabled=true;
	}
	if (labelType) {
		if ($quantifMethods{XIC} \|\| $quantifMethods{SILAC} \|\| $quantifMethods{ITRAQ} \|\| $quantifMethods{TMT}) {
			for (let i=1; i<algoSelOpt.length; i++) {
				if (algoSelOpt[i].value.match('MSstats\|DIA\|TDA')) continue; // disable MSstats
				if (algoSelOpt[i].value.match('PEP_RATIO') && labelType == 'FREE') continue; // disable PEP_RATIO algo
				algoSelOpt[i].disabled=false; // enable anything else
			}
		}
		if ($quantifMethods{DIA}) {
			//if ($nbCond > 2) {algoSelOpt[2].disabled=false;} // enable 'All vs All > MSstats (SWATH)'
			//algoSelOpt[5].disabled=false; // enable 'Common ref > MSstats (SWATH)'
			//algoSelOpt[8].disabled=false; // enable 'Peptide count > SSP Analysis'
			for (let i=1; i<algoSelOpt.length; i++) {
				if (algoSelOpt[i].value.match('MSstats\|DIA')) {algoSelOpt[i].disabled=false;} // enable MSstats
				else if (algoSelOpt[i].value=='SSPA') {algoSelOpt[i].disabled=false;} // enable SSPA
			}
		}
		if ($quantifMethods{TDA}) {
			for (let i=1; i<algoSelOpt.length; i++) {
				if (algoSelOpt[i].value.match('TDA:')) {algoSelOpt[i].disabled=false;} // enable TDA
				//else if (algoSelOpt[i].value=='SSPA') {algoSelOpt[i].disabled=false;} // enable SSPA
			}
		}
		if ($quantifMethods{NONE}) {
			for (let i=1; i<algoSelOpt.length; i++) {
				if (algoSelOpt[i].value=='SSPA') {algoSelOpt[i].disabled=false;} // enable SSPA
			}
		}
/*
		if (labelType == 'FREE') {
			//algoSelOpt[1].disabled=true; // disable 'Common ref > Simple ratio'
			//algoSelOpt[4].disabled=true; // disable 'Common ref > Super ratio'
		}
		else {
			algoSelOpt[1].disabled=true; // disable 'All vs All > Simple ratio' (not handled by new ML scripts)
			if (labelType != 'SILAC' \|\| $nbCond <= 2) {
				algoSelOpt[4].disabled=true; // disable 'Common ref > Super ratio'
			}
		}
*/
	}

	/* Update TopN */
	updateTopN();

	/* Update other options */
	document.getElementById('pepCleavage').selectedIndex=(labelType == 'SILAC')? 0 : 1; // Allow missed cleavage for SILAC
	document.getElementById('biasCorrect').options.length=2; // reset list of normalizations (except 'Select' and 'None')
	/* Select matching datasets & unselect all others */
	var optList=(labelType)? document.getElementsByName(labelType) : [];
	if (optList.length) {
		var idInfo=optList[0].id.split(':'); // Update 1st matching only (opt:obsID)
		var selQuanti=document.getElementById('anaXICQuanti:'+idInfo[1]);
/*
		for (let i=1; i<selQuanti.options.length; i++) { // skip '-= Select =-'
			if (selQuanti.options[i].value==optList[0].value) {
				selQuanti.selectedIndex=i;
				break;
			}
		}
*/
		optList[0].selected=true;
		propagateQuantifSelection(selQuanti,true); // true -> propagate selection to all
	}
	/* Unselect all datasets (Not really necessary) */
	else { // -= Select =- or labelType with no quantif (just to be safe)
		var allSelects=document.getElementsByTagName("select");
		for (let i=0; i<allSelects.length; i++) {
			if (allSelects[i].name.match('anaXICQuanti:')) {
				allSelects[i].selectedIndex=0;
			}
		}
	}
}
var popupTexts={}; // set below by observations XIC quantifications
function preparePopup(selectValue) {
	if (selectValue && popupTexts[selectValue]) {popup(popupTexts[selectValue]);}
}
function sequenceView(id_protein,anaIdStrg) { // Callable during manual reference peptides selection
	var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+anaIdStrg+"&id_prot="+id_protein+"&msdata="+top.showMSdata;
	top.openProtWindow(winLocation);
}
function cancelAction() {
	//top.promsFrame.selectedAction='summary'; // set default action to 'summary'
	top.promsFrame.optionFrame.selectOption();
}
function checkPtmSelection(myForm,multiPtmScope) {
	var ptmID='freeRes'; // default
	if (multiPtmScope) {
		if (!multiPtmScope.checked) {return 0;}
		ptmID=multiPtmScope.value;
	} 
	if (myForm.ptmScope.value.match(/proteinPTM\|create_multiPTM/)) {
		/* Check context */
		const contextStrg0=myForm['context_'+ptmID].value;
		if (contextStrg0) {
			// Remove leading/trailing spaces/- if any
			let contextStrg=contextStrg0.replace(/^[\\s-]+/,'');
			contextStrg=contextStrg.replace(/[\\s-]+\$/,'');
			// Evaluate all contexts
			let contexts=contextStrg.split(/\\s+/);
			let errorMessage;
			for (let i=0; i<contexts.length; i++) {
				contexts[i]=contexts[i].replace(/^-+/,''); // clean again
				contexts[i]=contexts[i].replace(/-+\$/,''); // clean again
				if (!contexts[i].length) continue;
				let matchRes=[];
				if (contexts[i].match(/[^-\\w\\{\\},\\/><]/) \|\| contexts[i].match(/[BJOUZ]/i)) { // All possible caracters
					errorMessage='Illegal character detected';
				}
				else if (contexts[i].toUpperCase()===contexts[i] \|\| contexts[i].match(/[A-Z][\\/]+[a-z]/) \|\| contexts[i].match(/[a-z][\\/]+[A-Z]/)) { // Mixture of upper and lowercases sites
					errorMessage='Site residue(s) can only be in lowercase';
				}
				else if (contexts[i].match(/.>\|<./)) {
					errorMessage='Wrong protein N/C-term declaration';
				}
				else if (contexts[i].match(/[a-z]-*[a-z]/)) { // other residues in uppercase
					errorMessage='Residue(s) around site must be in uppercase';
				}
				else if (contexts[i].match(/[a-z][^\\/]+[a-z]/)) { // Too many sites
					errorMessage='Only 1 position allowed for site residue(s)';
				}
				else if (contexts[i].match(/[a-z]\\{/)) { // Bad or missplaced {"min,max"}
					errorMessage='Range cannot apply to site';
				}
				else if (contexts[i].match(/(^\|[^\\w])\\{/) \|\| contexts[i].match(/\\{(\\D\|\$)/) \|\| contexts[i].match(/\\{\\}/) \|\| contexts[i].match(/(^\|\\D)\\}/)) { // check range boundaries
					errorMessage='Wrong range declaration';
				}
				else if (contexts[i].match(/[^\\d],/) \|\| contexts[i].match(/,[^\\d]/) ) { // check range ',' position
					errorMessage='Wrong range declaration';
				}
				else if ( contexts[i].match(/[^\\{\\d,]\\d/) \|\| contexts[i].match(/\\d[^\\d,\\}]/) ) { // check range number(s) position
					errorMessage='Wrong range declaration';
				}
				else if (matchRes=[...contexts[i].matchAll(/\\{(\\d+),(\\d+)\\}/g)]) { // scan all {min,max} to make sure min < max
					for (let j=0; j<matchRes.length; j++) {
						if (matchRes[j][1] >= matchRes[j][2]) {
							errorMessage='Wrong range declaration';
							break;
						}
					}
				}
				if (errorMessage) {
					alert('ERROR: '+errorMessage+' in motif "'+contexts[i]+'"');
					return -1; // implicite break
				}
			}
			myForm['context_'+ptmID].value=contexts.join(' '); // Update with cleaned context
		}
	}
	/* PTM ProbSoft settings */
	else if (myForm.ptmScope.value==='multiPTM' && ($isMultiPtmProb \|\| ptmID==$phosphoID) && (myForm.okProb.value < 0 \|\| myForm.okProb.value > 100)) {
		alert('ERROR: $ptmProbSoft probability must be between 0 and 100%');
		return -1;
	}
	return 1;
}
function checkForm(myForm) {
	/* Quantif name */
	if (!myForm.quantifName.value) {alert('ERROR: Missing name for quantification!'); return false;}
	
	/* PTM sites & residues (proteinPTM,multiPTM,create_multiPTM) */
	if (myForm.ptmScope.value.match('PTM')) {
		if (myForm.ptmScope.value==='create_multiPTM') { // Free residues
			if (!myForm.context_freeRes.value) {
				alert('ERROR: Enter at least 1 free residue to be quantified !');
				return false;
			}
			var ptmStatus=checkPtmSelection(myForm);
			if (ptmStatus==-1) {return false;}
		}
		else { // proteinPTM,multiPTM
			if (myForm.multiPtmScope.length) {
				for (let i=0; i<myForm.multiPtmScope.length; i++) {
					var ptmStatus=checkPtmSelection(myForm,myForm.multiPtmScope[i]);
					if (ptmStatus==1) {okPTM=true;}
					else if (ptmStatus==-1) {return false;}
				}
			}
			else {
				var ptmStatus=checkPtmSelection(myForm,myForm.multiPtmScope);
				if (ptmStatus==1) {okPTM=true;}
				else if (ptmStatus==-1) {return false;}
			}
			if (!okPTM) {
				alert('ERROR: Select at least 1 PTM for quantification !');
				return false;
			}
		}
	}

	/* PTM filter */
	if (Math.abs(myForm.pepPTM.value)==2) {
		var okPTM=false;
		if (myForm.customPTM.length) {
			for (let i=0; i<myForm.customPTM.length; i++) {
				if (myForm.customPTM[i].checked) {
					okPTM=true;
					break;
				}
			}
		}
		else {okPTM=myForm.customPTM.checked;}
		if (!okPTM) {
			alert('ERROR: No peptide modifications selected!');
			return false;
		}
	}

	/* Labeling */
	if (!myForm.labeling.value) {
		alert('ERROR: No labeling method selected.');
		return false;
	}
	var isLabeled=(myForm.labeling.value=='FREE')? false : true;

	/* Chosen Algorithm */
	if (myForm.algoType.value.match('select')) {
		alert('ERROR: Select an algorithm!');
		return false;
	}
	
	/* Abundance metrics */
	if (myForm.algoType.value.match('Abundance')) {
		let okMeas=false;
		for (let i=0; i<myForm.abundMeasures.length; i++) {
			if (myForm.abundMeasures[i].checked) {
				okMeas=true;
				break;
			}
		}
		if (!okMeas) {
			alert('ERROR: No abundance metric selected!');
			return false;
		}

		/* Absolute quantif */
		if (document.getElementById('absCheckbox').checked) {
			if (myForm.pqi.value == "") {
				alert('ERROR: Please choose a Protein Quantification Index !');
				return false;
			}
			if (myForm.abs_model.value=='calibration') {
				if (myForm.concentration_file.value == "") {
					alert('ERROR: Calibration model selected but no calibration file provided !');
					return false;
				} else {
					var ext = myForm.concentration_file.value.split('.').pop();
					if (ext != "csv" && ext != "txt") {
						alert('ERROR: Please provide a .csv file (or .txt with .csv formatting) as the calibration file.');
						return false;
					}
				}
			} else {  // proportionality model
				if ( (myForm.out_units.value == 'amounts') && (myForm.total_conc_file.value == "") ) {
					alert('ERROR: To get amounts you must provide a file with initial total quantities in your samples. Otherwise only the composition can be computed.');
					return false;
				} else if (myForm.out_units.value == 'amounts') {
					var ext = myForm.total_conc_file.value.split('.').pop();
					if (ext != "csv" && ext != "txt") {
						alert('ERROR: Please provide a .csv file (or .txt with .csv formatting) as the total quantity file.');
						return false;
					}
				}
			}
			if ((myForm.concentration_file.value != "") \|\| (myForm.total_conc_file.value != "")) {
				if (myForm.units.value == "") {
					alert('ERROR: You did not specify the unit of the provided quantities !');
					return false;
				}
			}
		}
	}
	else if (myForm.algoType.value.match('SSPA')) {
		let pValue=parseFloat(myForm.sspaPvalue.value);
		if (pValue <= 0 \|\| pValue > 1) {
			alert('ERROR: p-value for State-specificity threshold is out of range ]0-1] !');
			return false;
		}
		myForm.sspaPvalue.value=pValue;
	}

	/* States & Observations selection */
	var numCondition=0;
	var refStateAnaList={}; // for labeled quantif only
	for (let i=1; i <= $nbCond; i++) {
		var condID=myForm['state_'+i].value;
		if (condID=='Select') {
			if (i==1) {
				alert("ERROR: Reference state (State #1) must be set!");
				return false;
			}
			continue;
		}
		numCondition++;
		var obsCondOK=false;
		var checkMatchedAna=(!isLabeled \|\| myForm.algoType.value.match('SSPA\|PEP_INTENSITY'))? false : true;
		var matchedAna=!checkMatchedAna;
		var condObsBox=myForm['condObs:'+condID];
		if (condObsBox.length) {
			for (let j=0; j < condObsBox.length; j++) {
				if (i > 1 && obsCondOK && matchedAna) break;
				if (condObsBox[j].checked && myForm['anaXICQuanti:'+condObsBox[j].value].value) {
					obsCondOK=true;
					if (checkMatchedAna) {
						var obsData=myForm['anaXICQuanti:'+condObsBox[j].value].value.split(':');
						if (i==1) {refStateAnaList[obsData[2]]=1; matchedAna=true;} // records anaID
						else if (refStateAnaList[obsData[2]]) {matchedAna=true;} // at least 1 match is enough
					}
				}
			}
		}
		else if (condObsBox.checked && myForm['anaXICQuanti:'+condObsBox.value].value) {
			obsCondOK=true;
			if (checkMatchedAna) {
				var obsData=myForm['anaXICQuanti:'+condObsBox.value].value.split(':');
				if (i==1) {refStateAnaList[obsData[2]]=1;matchedAna=true;} // records anaID
				else if (refStateAnaList[obsData[2]]) {matchedAna=true;} // at least 1 match is enough
			}
		}
		if (!obsCondOK) {
			alert('ERROR: No observation selected for State #'+i+'!');
			return false;
		}
		if (checkMatchedAna && !matchedAna) {
			alert('ERROR: Reference (State #1) is not suitable for selected State #'+i+'! Check your design or state selection.');
			return false;
		}
	}
	if (numCondition < 2) {
		alert('ERROR: At least 2 distinct states have to be chosen!');
		return false;
	}
	
	/* Bias correction */
	if (!myForm.biasCorrect.value && myForm.algoType.value != 'SSPA') {
		alert('ERROR: No bias correction selected.');
		return false;
	}
	
	/* Multi-quantification */
	if (myForm.algoType.value.match(':singleRef') && myForm.multiQuantif.checked) {
		if (numCondition >= 3 && !myForm.quantifName.value.match(/%(TEST\|#)%/i) && !confirm('WARNING: All quantifications have same name. Proceed anyway?')) {
			return false;
		}
		if (!confirm((numCondition-1)+' quantifications will be launched. Proceed?')) {
			return false;
		}
	}
	return true;
}

function updatePQISelection(pqi) {
	var top_params = ['PEP_TOPX', 'PEP_STRICT', 'moreAbsDefaults', 'lessAbsDefaults'];
	var top_labels = ['PEP_TOPX_LABEL', 'PEP_STRICT_LABEL'];

	if (pqi == "aLFQ_TOP") {
		top_params.forEach(function(item, index, array) {
			document.getElementById(item).style.display = '';
			document.getElementById(item).disabled = false;
		});
		top_labels.forEach(function(item, index, array) {
			document.getElementById(item).style.display = '';
		});
		updateAbsParams('less');
	} else {
		top_params.forEach(function(item, index, array) {
			document.getElementById(item).disabled = true;
			document.getElementById(item).style.display = 'none';
		});
		top_labels.forEach(function(item, index, array) {
			document.getElementById(item).style.display = 'none';
		});
	}
}

function updateQtyFile() {
	var model = document.getElementById('ABS_MODEL').value;
	var output = document.getElementById('OUT_UNITS').value;

	if (model == "calibration") {
		document.getElementById('anchor_conc_label').style.display = '';
		document.getElementById('concentration_file').disabled = false;
		document.getElementById('concentration_file').style.display = '';
		document.getElementById('tpc_label').style.display = 'none';
		document.getElementById('total_conc_file').disabled = true;
		document.getElementById('total_conc_file').style.display = 'none';
		document.getElementById('total_conc_file').value = "";
		updateFileUnits(document.getElementById('concentration_file').value);
	} else if (output == 'amounts') {  // model is proportionality
		document.getElementById('anchor_conc_label').style.display = 'none';
		document.getElementById('concentration_file').disabled = true;
		document.getElementById('concentration_file').style.display = 'none';
		document.getElementById('concentration_file').value = "";
		document.getElementById('tpc_label').style.display = '';
		document.getElementById('total_conc_file').disabled = false;
		document.getElementById('total_conc_file').style.display = '';
		updateFileUnits(document.getElementById('total_conc_file').value);
	} else {  // model is proportionality and only composition is asked for as output
		document.getElementById('anchor_conc_label').style.display = 'none';
		document.getElementById('concentration_file').disabled = true;
		document.getElementById('concentration_file').style.display = 'none';
		document.getElementById('concentration_file').value = "";
		document.getElementById('tpc_label').style.display = 'non';
		document.getElementById('total_conc_file').disabled = true;
		document.getElementById('total_conc_file').style.display = 'none';
		document.getElementById('total_conc_file').value = "";
		updateFileUnits("");
	}
}

function updateAbsParams(act) {
	if (act=='more') {
		document.getElementById('moreAbsDefaults').style.display = 'none';
		document.getElementById('lessAbsDefaults').style.display = '';
		document.getElementById('absParamsDIV').style.display = '';
	}
	else {
		document.getElementById('moreAbsDefaults').style.display = '';
		document.getElementById('lessAbsDefaults').style.display = 'none';
		document.getElementById('absParamsDIV').style.display = 'none';
	}
}

function updateFileUnits(file) {
	if (file == "") {
		document.getElementById('units_label').style.display = 'none';
		document.getElementById('UNITS').style.display = 'none';
		document.getElementById('UNITS').disabled = true;
	} else {
		document.getElementById('units_label').style.display = '';
		document.getElementById('UNITS').style.display = '';
		document.getElementById('UNITS').disabled = false;
	}
}

</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" >
<CENTER>
<FONT class="title">$titleString</FONT>
$oldAlgoAccessStrg
<BR>
<FORM name="selAnaForm" method="post" onsubmit="return(checkForm(this));" enctype="multipart/form-data">
<INPUT type="hidden" name="ID" value="$designID">

<BR>
|;

##################################################
####>Displaying quantification-specific forms<####
##################################################
#my $disabSuperRatio=($projectAccess eq 'bioinfo')? '' : ' disabled'; # temporary
print qq
|<INPUT type="hidden" name="nbMaxCond" value="$nbCond">
<INPUT type="hidden" name="phosphoID" value="$phosphoID">
<INPUT type="hidden" name="ptmProbSoftCode" value="$ptmProbSoftCode">
<TABLE bgcolor="$darkColor">
<TR><TH class="title3" align=right valign=top>Name :</TH><TD bgcolor="$lightColor"><INPUT type="text" name="quantifName" value="" placeholder="Name of quantification" class="title3" style="width:400px"/>
<SPAN id="nameSPAN" style="display:none">
<BR>&nbsp;<B>Use: <INPUT type="button" class="keyword" value="%TEST%" onclick="document.selAnaForm.quantifName.value+=this.value"/>, <INPUT type="button" class="keyword" value="%REF%" onclick="document.selAnaForm.quantifName.value+=this.value"/> or  <INPUT type="button" class="keyword" value="%#%" onclick="document.selAnaForm.quantifName.value+=this.value"/> for test-state name, reference-state name and quantification rank respectively.</B>
</SPAN>
</TD></TR>
|;
if ($userInfo[1]=~/mass|bioinfo/) { # <<<<------------TEMPLATING HIDDEN (display:none) UNTIL MADE COMPATIBLE WITH NEW OPTIONS
	print qq
|<TR style="display:none"><TH class="title3" align=right valign=top>Apply template :</TH><TD bgcolor="$lightColor" class="title3" nowrap><SELECT id="template" name="template" class="title3" onchange="useTemplate(this.value)">
|;
	if (scalar keys %quantifTemplate) {
		print "<OPTION value=\"select\">-= Select =-</OPTION>\n";
		my ($labelStrg)=(scalar keys %{$designInfo{'LABEL'}} == 1)? keys %{$designInfo{'LABEL'}} : undef;
		if ($labelStrg) {
			foreach my $projID (sort {lc($projectTemplate{$a}) cmp lc($projectTemplate{$b})} keys %projectTemplate) {
				if ($userInfo[1] ne 'bio' || scalar keys %projectTemplate > 1 || $projID != $projectID) { # always displayed for bioinfo|mass
					print "<OPTION value=\"project\" disabled style=\"color:black;\">&nbsp;&nbsp;&bull;$projectTemplate{$projID}:</OPTION>\n";
				}
				foreach my $quantiID (sort {lc($quantifTemplate{$labelStrg}{$projID}{$a}{'NAME'}) cmp lc($quantifTemplate{$labelStrg}{$projID}{$b}{'NAME'})} keys %{$quantifTemplate{$labelStrg}{$projID}}) {
					print "<OPTION value=\"$quantiID\">$quantifTemplate{$labelStrg}{$projID}{$quantiID}{NAME}</OPTION>\n";
				}
			}
		}
		else {
			foreach my $label (sort {$a cmp $b} keys %quantifTemplate) {
				print "<OPTGROUP label=\"$label\">";
				foreach my $projID (sort {lc($projectTemplate{$a}) cmp lc($projectTemplate{$b})} keys %projectTemplate) {
					if ($userInfo[1] ne 'bio' || scalar keys %projectTemplate > 1 || $projID != $projectID) { # always displayed for bioinfo|mass
						print "<OPTION value=\"project\" disabled style=\"color:black;\">&nbsp;&nbsp;&bull;$projectTemplate{$projID}:</OPTION>\n";
					}
					foreach my $quantiID (sort {lc($quantifTemplate{$label}{$projID}{$a}{'NAME'}) cmp lc($quantifTemplate{$label}{$projID}{$b}{'NAME'})} keys %{$quantifTemplate{$label}{$projID}}) {
						print "<OPTION value=\"$quantiID\">$quantifTemplate{$label}{$projID}{$quantiID}{NAME}</OPTION>\n";
					}
				}
				print "</OPTGROUP>";
			}
		}
	}
	else{
		print "<OPTION value=\"select\">***No template recorded***</OPTION>";
	}
	print "</SELECT></TD></TR>\n";
}
print qq
|<TR><TH align=right class="title3" nowrap valign="top">&nbsp;Strategy :</TH><TD bgcolor="$lightColor">
<TABLE cellpadding=0 cellspacing=0><TR><TD class="title3" nowrap valign="top">
&nbsp;Labeling:<SELECT name="labeling" class="title3 template" onchange="clearTemplate();updateLabeling(this.value)"><OPTION value="">-= Select =-</OPTION>
|;
my $selLabelStrg=(scalar keys %{$designInfo{'LABEL'}} == 1)? 'selected' : '';
foreach my $labelType (sort keys %{$designInfo{'LABEL'}}) {
	print "<OPTION value=\"$labelType\" $selLabelStrg";
	#print ' selected' if scalar keys %{$designInfo{'LABEL'}}==1;
	if ($labelType eq 'FREE') {print ">Label-free</OPTION>";}
	elsif ($labelType eq 'ITRAQ') {print ">iTRAQ</OPTION>";}
	else {print ">$labelType</OPTION>";}
}
#my ($disabSites,$siteMessage)=($numPTMs)? ('','PTM sites') : ('disabled','**No project PTMs**');
# my $freeResText=($allowFreeResidues)? ' / Free residues' : '';
my $disabFreeRes=($allowFreeResidues)? 'disabled' : '';
print qq
|</SELECT>&nbsp;&nbsp;&nbsp;&nbsp;
</TD><TD class="title3" nowrap valign="top">
Method:<SELECT id="algoType" name="algoType" class="title3 template" onchange="updateAlgoTypeSelection(this.value)">
<OPTION value="select">-= Select =-</OPTION>
<OPTGROUP label="Peptide ratio (labeled data):">
<OPTION value="PEP_RATIO:SimpleRatio:singleRef" disabled>All vs State1</OPTION>
<OPTION value="PEP_RATIO:SuperRatio:multiRef" disabled>Super ratio</OPTION>
</OPTGROUP>
<OPTGROUP label="Peptide intensity:">
<OPTION value="PEP_INTENSITY:SimpleRatio:singleRef" disabled>All vs State1</OPTION>
<OPTION value="PEP_INTENSITY:SimpleRatio:multiRef" disabled>All vs All</OPTION>
<OPTION value="PEP_INTENSITY:Abundance" disabled>Abundance in States</OPTION>
</OPTGROUP>
<OPTGROUP label="DIA/SWATH:">
<OPTION value="MSstats:SimpleRatio:singleRef" disabled>All vs State1 (MSstats)</OPTION>
<OPTION value="MSstats:SimpleRatio:multiRef" disabled>All vs All (MSstats)</OPTION>
<OPTION value="DIA:SimpleRatio:singleRef" disabled>All vs State1</OPTION>
<OPTION value="DIA:SimpleRatio:multiRef" disabled>All vs All</OPTION>
<OPTION value="DIA:Abundance" disabled>Abundance in States</OPTION>
</OPTGROUP>
<OPTGROUP label="TDA (PRM,SRM,MRM):">
<OPTION value="TDA:SimpleRatio:singleRef" disabled>All vs State1</OPTION>
<OPTION value="TDA:SimpleRatio:multiRef" disabled>All vs All</OPTION>
</OPTGROUP>
<OPTGROUP label="Feature counts:">
<OPTION value="SSPA" disabled>SSP Analysis</OPTION>
</OPTGROUP>
</SELECT>&nbsp;&nbsp;&nbsp;&nbsp;
<SPAN id="abundMeasuresSPAN" style="font-size:12px; display:none">Abundance metrics:<BR>
&nbsp;&nbsp;<LABEL><INPUT type="checkbox" name="abundMeasures" class="template" value="SUM_INT" onclick="updateBiasCorrectionLevel()" checked/> Sum</LABEL><BR>
&nbsp;&nbsp;<LABEL><INPUT type="checkbox" name="abundMeasures" class="template" value="MEAN_INT" onclick="updateBiasCorrectionLevel()"/> Mean</LABEL><BR>
&nbsp;&nbsp;<LABEL><INPUT type="checkbox" name="abundMeasures" class="template" value="GEO_MEAN_INT" onclick="updateBiasCorrectionLevel()"/> Geometric mean</LABEL><BR>
&nbsp;&nbsp;<LABEL><INPUT type="checkbox" name="abundMeasures" class="template" value="MEDIAN_INT" onclick="updateBiasCorrectionLevel()"/> Median</LABEL><BR>
&nbsp;&nbsp;<LABEL><INPUT type="checkbox" name="abundMeasures" class="template" value="MY_LFQ" onclick="updateBiasCorrectionLevel(); document.getElementById('lfqRatioDIV').style.display=(this.checked)? '' : 'none'" checked/> LFQ:</LABEL>
<DIV id="lfqRatioDIV" style="font-weight:normal;margin-left:50px">
- Use &ge;<INPUT type="number" name="minRatioLFQ" class="template" value=2 min=1 max=5 style="width:35px;height:20px"/> peptide ratios<BR>
-<LABEL><INPUT type="checkbox" name="stabRatioLFQ" value=1 checked/>Stabilize large ratios</LABEL>
</DIV><DIV>
&nbsp;&nbsp;<LABEL><INPUT type="checkbox" id="absCheckbox" name="abundMeasures" class="template" value="ABSOLUTE" onclick="document.getElementById('AbsoluteDIV').style.display=(this.checked)? '' : 'none'"/> Absolute quantification:</LABEL>
	<DIV id="AbsoluteDIV" style="font-weight:normal;margin-left:50px;display:none">
		<LABEL id="pqi_label" for="PQI">- Protein Quantification Index<SUP onmouseover="popup(helpText.pqiMethod)" onmouseout="popout()">?</SUP></LABEL>
		<SELECT id="PQI" name="pqi" onchange="updatePQISelection(this.value)">
			<OPTION value="" selected>-= Select =-</OPTION>
			<OPTION value="aLFQ_iBAQ">iBAQ</OPTION>
			<OPTION value="aLFQ_TOP">Top</OPTION>
		</SELECT>
		<BR>
		<LABEL id="pqi_norm_label" for="PQI_NORM">- PQI normalization<SUP onmouseover="popup(helpText.pqiNormalization)" onmouseout="popout()">?</SUP></LABEL>
		<SELECT id="PQI_NORM" name="pqi_normalization">
			<OPTION value="global" selected>Globally</OPTION>
			<OPTION value="per_state">Per Condition</OPTION>
			<OPTION value="none">None</OPTION>
		</SELECT>
		<BR>
		<LABEL id="out_label" for="OUT_UNITS">- Desired type of values<SUP onmouseover="popup(helpText.outUnits)" onmouseout="popout()">?</SUP></LABEL>
		<SELECT id="OUT_UNITS" name="out_units" onchange="updateQtyFile()">
			<OPTION value="composition" selected>Composition only (mol%, mass%)</OPTION>
			<OPTION value="amounts">Composition + amounts (mol, mass etc.)</OPTION>
		</SELECT>
		<BR>
		<LABEL id="abs_model_label" for="ABS_MODEL">- Quantification model<SUP onmouseover="popup(helpText.absQuantifModel)" onmouseout="popout()">?</SUP></LABEL>
		<SELECT id="ABS_MODEL" name="abs_model" onchange="updateQtyFile()">
			<OPTION value="proportionality" selected>Direct proportionality</OPTION>
			<OPTION value="calibration">Calibration curve</OPTION>
		</SELECT>
		<DIV>
		<LABEL id="tpc_label" for="total_conc_file" style="display:none">- Total protein quantity of each sample (.csv file, optionnal for composition)<SUP onmouseover="popup(helpText.quantityFile)" onmouseout="popout()">?</SUP></LABEL>
		<LABEL id="anchor_conc_label" for="concentration_file" style="display:none">- Quantity of anchor proteins (.csv file)<SUP onmouseover="popup(helpText.quantityFile)" onmouseout="popout()">?</SUP></LABEL>
		</DIV><DIV>
		<INPUT type="file" id="total_conc_file" name="total_conc_file" accept=".csv" onchange="updateFileUnits(this.value)" style="display:none" disabled>
		<INPUT type="file" id="concentration_file" name="concentration_file" accept=".csv" onchange="updateFileUnits(this.value)" style="display:none" disabled>
		</DIV><DIV>
		<LABEL id="units_label" for="UNITS" style="display:none">- Unit of the provided quantities<SUP onmouseover="popup(helpText.providedQuantities)" onmouseout="popout()">?</SUP></LABEL>
		<SELECT id="UNITS" name="units" style="display:none" disabled>
			<OPTION value="" selected>-= Select =-</OPTION>
			<OPTGROUP label="Number">
				<OPTION value="mol">mol</OPTION>
				<OPTION value="mmol">mmol</OPTION>
				<OPTION value="copy_number">Copy Number</OPTION>
			</OPTGROUP>
			<OPTGROUP label="Concentration">
				<OPTION value="conc_mol_l">mol/L or mmol/mL</OPTION>
				<OPTION value="conc_mmol_l">mmol/L</OPTION>
			</OPTGROUP>
			<OPTGROUP label="Mass">
				<OPTION value="mass_g">g</OPTION>
				<OPTION value="mass_mg">mg</OPTION>
				<OPTION value="mass_ug">μg</OPTION>
			</OPTGROUP>
			<OPTGROUP label="Mass Concentration">
				<OPTION value="mass_conc_g_l">g/L or mg/mL or μg/μL</OPTION>
				<OPTION value="mass_conc_mg_l">mg/L or μg/mL</OPTION>
				<OPTION value="mass_conc_ug_l">μg/L</OPTION>
			</OPTGROUP>
		</SELECT>
		</DIV><DIV>
		<INPUT type="button" id="moreAbsDefaults" class="font11" value="More parameters" onclick="updateAbsParams('more')" style="display:none"/>
		<INPUT type="button" id="lessAbsDefaults" class="font11" value="Less parameters" style="display:none" onclick="updateAbsParams('less')"/>
		</DIV><DIV id="absParamsDIV" style="display:none">
			<TABLE cellpadding=0 cellspacing=0>
				<TR>
					<TD align=right nowrap><LABEL id="PEP_TOPX_LABEL" for="PEP_TOPX" style="display:none">Peptide TopX<SUP onmouseover="popup(helpText.peptideTopX)" onmouseout="popout()">?</SUP>:</LABEL></TD>
					<TD nowrap><INPUT type=number id="PEP_TOPX" name="pep_topx" value="3" min=1 step=1 style="width:40px;display:none" disabled></TD>
				</TR>
				<TR>
					<TD align=right nowrap><LABEL id="PEP_STRICT_LABEL" for="PEP_STRICT" style="display:none">Peptide strictness<SUP onmouseover="popup(helpText.peptideStrictness)" onmouseout="popout()">?</SUP>:</LABEL></TD>
					<TD nowrap>
						<SELECT id="PEP_STRICT" name="pep_strictness" style="display:none" disabled>
							<OPTION value="strict">Strict</OPTION>
							<OPTION value="loose" selected>Loose</OPTION>
						</SELECT>
					</TD>
				</TR>
			</TABLE>
		</DIV>
	</DIV>
</DIV>
</SPAN>
</TD>
<TD class="title3" nowrap valign="top">
<SPAN id="topSPAN" class="template" style="display:none">Use:<SELECT name="matchingPep" class="title3 template" disabled><OPTION value="0" disabled>any</OPTION><OPTION value="1" selected>matching</OPTION></SELECT>
&nbsp;top<SELECT name="topN" class="title3 template" disabled><OPTION value="1">1</OPTION><OPTION value="2">2</OPTION><OPTION value="3" selected>3</OPTION><OPTION value="5">5</OPTION><OPTION value="10">10</OPTION><OPTION value="0">N</OPTION></SELECT>&nbsp;peptides</SPAN>
<SPAN id="trackTrueSPAN" style="font-size:12px; display:none"><LABEL><INPUT type="checkbox" name="trackMsmsPep" value="1" checked/>Track truly identified peptides</LABEL><SUP onmouseover="popup(helpText.trackMS2)" onmouseout="popout()">?</SUP></SPAN>
<SPAN id="sspaPepSPAN" style="display:none" class="template">&nbsp;&nbsp;&nbsp;&nbsp;Use:<SELECT name="pepFocus" class="title3 template" onchange="const mbwrCkb=document.selAnaForm.pepRescued; mbwrCkb.checked=(this.value==='xic')? true : false; mbwrCkb.disabled=(this.value==='sp_count')? true : false;"><OPTION value="xic">Peptide abundance</OPTION><OPTION value="sp_count" selected>Peptide/spectrum matches</OPTION><OPTION value="all_ion">All peptide ions</OPTION><OPTION value="all_pep">All peptides</OPTION><OPTION value="dist_ion">Distinct peptide ions</OPTION><OPTION value="dist_pep">Distinct peptides</OPTION><OPTION value="dist_seq">Distinct peptide sequences</OPTION></SELECT></SPAN></TD></TR></TABLE>
</TD></TR>
<TR><TH class="title3" align=right valign=top>Focus<SUP onmouseover="popup(helpText.scopePTM)" onmouseout="popout()">?</SUP> :</TH><TH align=left bgcolor="$lightColor" nowrap><SELECT name="ptmScope" id="ptmScope" class="title3 template" onchange="updateQuantifScope(this.value)">
<OPTION value="protein">Whole proteins</OPTION>
<OPTION value="proteinPTM">Whole proteins by PTM</OPTION>
<OPTION value="multiPTM">PTM sites</OPTION>
<OPTION value="create_multiPTM" $disabFreeRes>Free residues</OPTION>
</SELECT>
<SPAN id="seqContextSPAN" style="display:none" class="bold template">&nbsp;<LABEL><INPUT type="checkbox" name="keepSeqContext" value="1" class="template">Keep sequence context</LABEL><SUP onmouseover="popup(helpText.seqContext)" onmouseout="popout()">?</SUP></SPAN>
|;

## UL with all PTMs
print "<DIV id=\"multiPtmDIV\" style=\"display:none;max-height:200px;overflow:auto\"><TABLE cellspacing=0>\n"; # class=\"ulList\"
if ($numPTMs > scalar keys %unusedProjVarMods) {
	# print "<TR><TD class=\"bold\">&nbsp;PTM sites:</TD></TR>\n" if scalar keys %unusedProjVarMods;
	foreach my $modID (sort{lc($projectVarMods{$a}[0]) cmp lc($projectVarMods{$b}[0])} keys %projectVarMods) {
		next if $unusedProjVarMods{$modID};
		print qq
|<TR><TD>&nbsp;&nbsp;<LABEL><INPUT type="checkbox" name="multiPtmScope" value="$modID" onchange="updateMultiPtmSelection()">$projectVarMods{$modID}[0] [<FONT style="font-weight:bold;color:#$projectVarMods{$modID}[2]">$projectVarMods{$modID}[1]</FONT>]-sites</LABEL>
<SPAN id="contextSPAN_$modID" style="margin-left:5px;display:none">Context motif(s)<SUP onmouseover="popup(helpText.contextPTM)" onmouseout="popout()">?</SUP>:<INPUT type="text" name="context_$modID" style="width:150px" value=""/></SPAN>
|;
		if ($modID==$phosphoID && $ptmProbSoftCode && !$isMultiPtmProb) {
			print qq
|<SPAN id="ptmProbSPAN" style="margin-left:5px;display:none" class="template">Positions with $ptmProbSoft probability&ge;<INPUT type="number" name="okProb" class="template" min="0" max="100" value="75" size="4">% are confirmed. The others are <SELECT name="badProb" class="template"><OPTION value="exclude" selected>excluded</OPTION><OPTION value="ambiguous">delocalized</OPTION></SELECT><SUP onmouseover="popup(helpText.delocalized)" onmouseout="popout()">?</SUP></SPAN>
|;
		}
		print "</TD></TR>\n";
	}
}
print "</TABLE></DIV>\n";
if ($isMultiPtmProb) { # Position probabilities available for all modifs
	print qq
|<SPAN id="ptmProbSPAN" class="bold" style="display:none" class="bold template">&nbsp;Positions with $ptmProbSoft probability&ge;<INPUT type="number" name="okProb" class="template" min="0" max="100" value="75" size="4">% are confirmed. The others are <SELECT name="badProb" class="bold template"><OPTION value="exclude" selected>excluded</OPTION><OPTION value="ambiguous">delocalized</OPTION></SELECT><SUP onmouseover="popup(helpText.delocalized)" onmouseout="popout()">?</SUP></SPAN>
|;
}

##>Free residues<##
if ($allowFreeResidues) {
	print qq
|<DIV id="createPtmDIV" style="display:none" class="bold">
&nbsp;Residue(s) context motif(s)<SUP onmouseover="popup(helpText.contextPTM)" onmouseout="popout()">?</SUP>:<INPUT type="text" name="context_freeRes" style="width:200px" value=""/>
</DIV>
|;
# 	print qq
# |<DIV id="createPtmDIV" style="display:none;max-height:200px;overflow:auto"><TABLE cellspacing=0>
# <TR><TD><B>&nbsp;Add new residue:</B><SELECT onchange="if (this.value) {document.getElementById('freeResSPAN_'+this.value).style.display='';}">
# <OPTION value="">-= Select =-</OPTION>
# |;
# 	my $resIdx=-1;
# 	foreach my $res (@freeResList) {
# 		$resIdx++;
# 		my $resValue=($res=~/Nter/)? '-' : ($res=~/Cter/)? '+' : $res;
# 		my $qResValue=quotemeta($resValue);
# 		next if $freeResSpecif=~/$qResValue/;
# 		print "<OPTION value=\"$resIdx\">$res</OPTION>\n";
# 	}
# 	print qq
# |</SELECT></TD></TR>
# <TR><TD><INPUT type="hidden" name="freeResSpecif" value="$freeResSpecif"/>
# |;
# 	$resIdx=-1;
# 	foreach my $res (@freeResList) {
# 		$resIdx++;
# 		my $resValue=($res=~/Nter/)? '-' : ($res=~/Cter/)? '+' : $res;
# 		my $qResValue=quotemeta($resValue);
# 		my $displayStrg=($freeResSpecif=~/$qResValue/)? '' : 'style="display:none"';
# 		print qq |<SPAN id="freeResSPAN_$resIdx" $displayStrg>&nbsp;&nbsp;<LABEL><INPUT type="checkbox" name="freeRes" class="template" value="$resValue">$res</LABEL><BR></SPAN>|;
# 	}
# 	print qq
# |</TD></TR>
# </TABLE></DIV>
# |;
}

print qq
|<DIV id="refQuantifDIV" style="display:none" class="bold">
&nbsp;Protein-level fold change correction<SUP onmouseover="popup(helpText.protLevelCorr)" onmouseout="popout()">?</SUP>:<SELECT name="refQuantifType" class="bold" onchange="ajaxFetchReferenceQuantif(this.value)"><OPTION value="0">None</OPTION><OPTION value="-1">Use current dataset</OPTION><OPTGROUP label="Use dataset from:">
|;
foreach my $refDes (@referenceDesign) {
	print "<OPTION value=\"$refDes->[0]\">$refDes->[1]</OPTION>\n";
}
print "<OPTION value=\"0\" disabled>** None found **</OPTION>\n" unless scalar @referenceDesign;
print qq
|</SELECT>
<SPAN id="refQuantifSPAN"><INPUT type="hidden" name="refQuantif" value="-1"/></SPAN>&nbsp;
</DIV>
</TH></TR>
<TR><TH class="title3" align=right valign=top>States :</TH><TD bgcolor=$lightColor>
<INPUT type="button" value="Use all states" onclick="autoSelectStates('all')">&nbsp;<B>or use state range:</B><INPUT type="text" name="stateRange" style="width:80px" value="" placeholder="e.g. 1,5-8"><INPUT type="button" value="Apply" onclick="autoSelectStates('range')">&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Clear all states" onclick="clearAllStates()">
<DIV style="max-height:150px;overflow:auto"><TABLE>
|;
# print the states that would be chosen for ratio computation
foreach my $y (1..$nbCond) {
	print "<TR><TH align=right nowrap>#$y:</TH><TD><SELECT name=\"state_$y\" onchange=\"updateStateSelection(this,$y)\">$selCondStrg1</SELECT>";
	print "<SPAN id=\"refCommentSPAN\" class=\"template\" style=\"color:#DD0000;font-weight:bold\"></SPAN>" if $y==1;
	print "</TD></TR>\n";
}

print qq
|</TABLE></DIV><SPAN id="multiQuantifSPAN" class="template" style="display:none"><INPUT type="checkbox" name="multiQuantif" class="template" id="multiQuantifCHK" value=1 disabled><LABEL style="color:#DD0000;font-weight:bold" for="multiQuantifCHK">Perform a separate quantification for each ratio</LABEL></SPAN>
</TD></TR>\n
|;
#	if ($quantifType eq 'DESIGN') {
#		print qq
#|<TR><TH align=right class="title3" nowrap>&nbsp;Labeling method :</TH><TD bgcolor=$lightColor nowrap><SELECT name="labeling" class="title3" onchange="updateLabeling(this.value)"><OPTION value="">-= Select =-</OPTION>|;
#		foreach my $labelType (sort keys %{$designInfo{'LABEL'}}) {
#			print "<OPTION value=\"$labelType\"";
#			print ' selected' if scalar keys %{$designInfo{'LABEL'}}==1;
#			if ($labelType eq 'FREE') {print ">Label-free</OPTION>";}
#			elsif ($labelType eq 'ITRAQ') {print ">iTRAQ</OPTION>";}
#			else {print ">$labelType</OPTION>";}
#		}
#		print "</SELECT></TD></TR>\n";
#	}
print qq
|<TR><TH align=right nowrap valign=top>&nbsp;Peptide selection :</TH><TD bgcolor="$lightColor" nowrap>
	<TABLE cellspacing=0 cellpadding=0><TR>
	<TD valign="top" nowrap>&nbsp;<B>Specificity<SUP onmouseover="popup(helpText.specificity)" onmouseout="popout()">?</SUP>:</B><SELECT name="pepSpecifity" id="pepSpecifity" class="template" onchange="updateMGSharedPeptides(this.value)"><OPTION value="unique">Proteotypic</OPTION><OPTION value="unique_shared">Proteotypic & others</OPTION><OPTION value="all" selected>All</OPTION></SELECT>
	&nbsp;&nbsp;&nbsp;<B>Missed cleav.:<B><SELECT name="pepCleavage" id="pepCleavage" class="template"><OPTION value="1">Allowed</OPTION><OPTION value="0" selected>Not allowed</OPTION><OPTION value="-1">Also exclude cleaved counterparts</OPTION></SELECT>
	&nbsp;&nbsp;&nbsp;<B>Modifications<SUP onmouseover="popup(helpText.modifications)" onmouseout="popout()">?</SUP>:<B>
	</TD>
	<TD valign="top" nowrap rowspan=2><SELECT name="pepPTM" onchange="updateCustomPtmDisplay(this.value)" class="template"><OPTION value="1">Allow all</OPTION><OPTION value="2">Allow selected</OPTION>
	<OPTION value="0" selected>Exclude all</OPTION><OPTION value="-1">Exclude all & unmodified forms</OPTION>
	<OPTION value="-2">Exclude selected & unmodified forms</OPTION>
	</SELECT><BR>
	<DIV id="customPtmDIV" class="template" style="display:none">|;
foreach my $modID (sort{lc($modifications{$a}) cmp lc($modifications{$b})} keys %modifications) {
	$modifications{$modID}=&promsMod::resize($modifications{$modID},40);
	print "<LABEL><INPUT type=\"checkbox\" class=\"template\" name=\"customPTM\" value=\"$modID\"/> $modifications{$modID}</LABEL><BR>\n";
}
print "<FONT style=\"color:#DD0000\">There are no modifications involved!</FONT>\n" unless scalar keys %modifications;
print "<INPUT type=\"hidden\" name=\"listPTM\" value=\"",join(',',sort{$a<=>$b} keys %modifications),"\">\n";
my $disabManualStrg=($numProtInDesign > $MAX_PROT_MANUAL_PEP)? 'disabled' : '';
print qq 
|</DIV></TD>
	<TD valign="top" nowrap rowspan=2><DIV id="ratioPepDIV" class="template">&nbsp;&nbsp;&nbsp;<B>Charges:<B><SELECT name="pepCharge" class="template" onchange="updatePeptideSource(this.value)"><OPTION value="all">All</OPTION><OPTION value="best">Best signal</OPTION></SELECT>
&nbsp;&nbsp;&nbsp;<B>Sources:<B><SELECT name="pepSource" class="template" id="pepSource"><OPTION value="all">All</OPTION><OPTION value="best">Best signal</OPTION></SELECT>
</DIV>
</TD></TR>
<TR><TD colspan=2>
	<DIV id="MGSharedPepDIV">&nbsp;<B>Peptides shared between Match Groups<SUP onmouseover="popup(helpText.MGSharedPep)" onmouseout="popout()">?</SUP>:</B><SELECT name="MGSharedPep" id="MGSharedPep"><OPTION value="best">Assign to best</OPTION><OPTION value="share">Share</OPTION><OPTION value="exclude">Exclude</OPTION></SELECT></DIV>
	<DIV id="rescuedPepDIV" style="display:none"><LABEL><INPUT type="checkbox" name="pepRescued" value="1"/><B>Include data from MBWR-rescued peptides<SUP onmouseover="popup(helpText.includeMBWR)" onmouseout="popout()">?</SUP></B></LABEL></DIV>
</TD></TR>
</TABLE>
<SPAN id="manualRefPepSPAN" style="display:none">
&nbsp;<B>Protein-level fold change correction dataset:</B><SELECT name="manualRefPep" onchange="ajaxFetchRefProtPeptides(this.value,'protCorr')"><OPTION value="0">Current rules</OPTION><OPTION value="manual" $disabManualStrg>Manual selection</OPTION></SELECT>
</SPAN>
<DIV id="manualRefPepDIV" style="background-color:white;padding:2px;max-height:500px;overflow:auto;display:none"></DIV>
</TD></TR>
<TR><TH align=right nowrap>Protein selection<SUP onmouseover="popup(helpText.protSelection)" onmouseout="popout()">?</SUP> :</TH><TD bgcolor=$lightColor nowrap><SELECT name="protSelType" class="template"><OPTION value="exclude">Exclude</OPTION><OPTION value="restrict">Restrict to</OPTION></SELECT>&nbsp;<B>proteins from List:</B><SELECT name="protSelection" class="template"><OPTION value="">-= Select =-</OPTION>
|;
foreach my $classID (sort{lc($categoryList{$a}{'CL_NAME'}) cmp lc($categoryList{$b}{'CL_NAME'})} keys %categoryList) {
	print "<OPTGROUP label=\"$categoryList{$classID}{CL_NAME}:\">\n";
	foreach my $refCat (@{$categoryList{$classID}{'CAT'}}) {
		my $typeStrg=($refCat->[2] eq 'SITE')? ' [Sites]' : '';
		print '<OPTION value="',$refCat->[0],'">',$refCat->[1],"$typeStrg</OPTION>\n";
	}
	print "</OPTGROUP>\n";
}
print qq
|</SELECT>
&nbsp;&nbsp;&nbsp;<LABEL><INPUT type="checkbox" name="exclContaminants" value="1"/><B>Exclude contaminants<SUP onmouseover="popup(helpText.exclContaminants)" onmouseout="popout()">?</SUP></B></LABEL>
<SPAN id="hiddenProtSPAN" style="display:none">&nbsp;&nbsp;&nbsp;<LABEL><INPUT type="checkbox" name="protHidden" value="1"/><B>Include proteins where hidden if also found visible<SUP onmouseover="popup(helpText.visHidProt)" onmouseout="popout()">?</SUP></B></LABEL></SPAN>
&nbsp;</TD></TR>
<TR><TH align=right nowrap valign=top>&nbsp;Quantification settings :
	</TH><TD bgcolor=$lightColor nowrap valign=top>
<DIV id="noSspaSettingsDIV">
	<TABLE cellpadding=0 cellspacing=0><TR>
	<TH align=right valign="top">&nbsp;-Bias correction:</TH>
	<TD><SELECT name="biasCorrect" class="template" id="biasCorrect" onchange="updateBias(this.value)"><OPTION value="">-= Select =-</OPTION><OPTION value="none.none">None</OPTION></SELECT>&nbsp;</TD>
	<TD valign="top" nowrap><DIV id="biasProtDIV" class="template" style="visibility:hidden">&nbsp;<B><SELECT name="refProtSelType" class="template"><OPTION value="use">Use</OPTION><OPTION value="not" disabled>Exclude</OPTION></SELECT>
data from</B> <SELECT name="refProt" class="template" onchange="updateBiasCorrectionSet(this.value)"><OPTION value="all">All proteins</OPTION><OPTION value="manual" disabled>Manual peptide selection (TDA only)</OPTION>
|;
foreach my $classID (sort{lc($categoryList{$a}{'CL_NAME'}) cmp lc($categoryList{$b}{'CL_NAME'})} keys %categoryList) {
	print "<OPTGROUP label=\"$categoryList{$classID}{CL_NAME}:\">\n";
	foreach my $refCat (@{$categoryList{$classID}{'CAT'}}) {
		my $typeStrg=($refCat->[2] eq 'SITE')? ' [Sites]' : '';
		print '<OPTION value="',$refCat->[0],'">',$refCat->[1],"$typeStrg</OPTION>\n";
	}
	print "</OPTGROUP>\n";
}
# Accepted p-value correction methods:
# Algo=(Super/Simple)Ratio (new): "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none" (fdr <=> FDR-BH)
# Algo=Ratio|TnPQ (old): "FDR-BH" (Benjamini-Hochberg), "FDR-ABH" (Benjamini-Hochberg (Adaptative)), "FWER-BONF" (Bonferroni), "Qvalue" (Qvalue (Storey et al.))
print qq
|	</SELECT> <B>(<LABEL><INPUT type="checkbox" name="normSharedPep" value="1">use only shared peptides</LABEL><SUP onmouseover="popup(helpText.normSharedPep)" onmouseout="popout()">?</SUP>) for normalization.</B>&nbsp;</DIV></TD>
	</TR>
	<TR>
		<TD colspan=3><SPAN id="biasCorrLevelSPAN" style="margin-left:40px; display:none"><B>Perform correction at<SUP onmouseover="popup(helpText.biasCorrLevel)" onmouseout="popout()">?</SUP>:</B><SELECT name="biasCorrLevel" id="biasCorrTgt"><OPTION value="ion">Peptide</OPTION><OPTION value="prot">Protein</OPTION><OPTION value="ion_prot">Peptide & protein</OPTION></SELECT><B>-level</B></SPAN></TD>
	</TR>
	</TABLE>
<DIV id="manualBiasPepDIV" style="background-color:white;padding:2px;max-height:500px;overflow:auto;display:none"></DIV>
	&nbsp;<B>-Biological replicates:</B><SELECT name="residualVar" class="template" onchange="document.getElementById('resVarSPAN').style.visibility=(this.value=='biological')? 'visible' : 'hidden'"><OPTION value="auto" selected>Auto-detect*</OPTION><OPTION value="technical">No</OPTION><OPTION value="biological">Yes</OPTION></SELECT>
<SPAN id="resVarSPAN" style="visibility:hidden" class="template"><FONT style="font-weight:bold;color:#DD0000"> Requires at least 2 biological replicates in each state</FONT></SPAN><BR>
	<DIV id="minInfRatiosDIV" style="display:none">&nbsp;<B>-</B><LABEL><INPUT type="checkbox" name="minInfRatios" class="template" value="1"><B>Avoid infinite ratios when possible.<SUP onmouseover="popup(helpText.avoidInfRatio)" onmouseout="popout()">?</SUP></B></LABEL></DIV>
</DIV>
<DIV id="sspaSettingsDIV" style="display:none">
&nbsp;<B>-State specificity thresholds<SUP onmouseover="popup(helpText.sspaThresholds)" onmouseout="popout()">?</SUP>:
&nbsp;Best delta &ge;</B><INPUT name="sspaBestDelta" type="number" value="50" min="0" max="100" style="width:45px"/><B>%</B>&nbsp;
<SELECT name="sspaLogic"><OPTION value="or">or</OPTION><OPTION value="and">and</OPTION></SELECT>
&nbsp;<B>p-value &le;</B><INPUT name="sspaPvalue" type="text" value="0.05" style="width:80px"/>
</DIV>
	&nbsp;<B>-p-value correction:</B><SELECT name="fdrControl" class="template"><OPTION value="none">None</OPTION><OPTION value="fdr" selected>Benjamini-Hochberg*</OPTION><OPTION value="bonferroni">Bonferroni</OPTION></SELECT><BR>
	&nbsp;<SMALL><SUP>*</SUP>Recommanded options</SMALL>
	</TD></TR>
|;
if ($userInfo[1] =~ /mass|bioinfo/) {
	print qq
|<TR><TH align=right>Visibility :</TH><TD bgcolor="$lightColor"><SELECT name="quantifVisibility"><OPTION value="1">Public</OPTION><OPTION value="2">Hide from collaborators</OPTION>|;
	print '<OPTION value="3">For bioinformaticians only</OPTION>' if $userInfo[1] eq 'bioinfo';
	print "</SELECT></TD></TR>\n";
}
print "</TABLE>\n";

my %itemIcones=&promsConfig::getItemIcones;
my $bgColor=$darkColor;
my %prevItemName;
my $disabSubmit=' disabled';
print qq
|<BR>
<INPUT type="submit" name="launch" value="Launch Quantification" class="title3">&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="cancelAction();">
<BR><BR><FONT class="title2">Observations associated with State
<SELECT name="fillDiv" id="fillDiv" class="title3" onchange="showDivCond(this.value)">$selCondStrg2</SELECT>: </FONT>
<BR>
|;
###> Print the Design
foreach my $expConditionID (sort{ $designInfo{'EXPCONDITION'}{$a}{'DISPLAY_POS'} <=> $designInfo{'EXPCONDITION'}{$b}{'DISPLAY_POS'} } keys %{$designInfo{'EXPCONDITION'}}) {
	#print "<DIV id=\"quantiParam\" name=\"$expConditionID\" style=\"display:none\">\n";
	print "<DIV id=\"$expConditionID\" name=\"condDiv\" style=\"max-height:400px;overflow-y:scroll;display:none\">\n";
	print "<TABLE border=0 cellspacing=0>\n";
	###> Header
	print qq
|<TR bgcolor=$darkColor>
<TH class="rbBorder title3" colspan=2 valign="middle">&nbsp;Observations&nbsp;</TH>
<TH class="bBorder" colspan=2>&nbsp;Peptide quantifications&nbsp;<BR>&nbsp;Extend selection to:<SELECT name="propagQuantScope" id="propagQuantScope:$expConditionID" onchange="synchronizePropagation(this)"><OPTION value="none">None</OPTION><OPTION value="state" selected>This state</OPTION><OPTION value="all">All states</OPTION></SELECT>&nbsp;</TH>
</TR>
|;
	###> Each Rows.
	my $bgColor=$darkColor;
	my $prevFracGroup=0;
	foreach my $anaID (sort{$designInfo{'ANALYSIS'}{$a}{'POS'} <=> $designInfo{'ANALYSIS'}{$b}{'POS'}} keys %{$designInfo{'EXPCONDITION'}{$expConditionID}{'ANALYSIS'}}) {
		foreach my $refObsData (sort{$a->[0]<=>$b->[0]} @{$anaObs{$expConditionID}{$anaID}}) {
			my ($obsID,$obsCode,$fracGroup,$techRepGroup)=@{$refObsData};
			$fracGroup=0 unless $fracGroup;
			if ($fracGroup==0 || ($fracGroup != $prevFracGroup)) {
				$bgColor=($bgColor eq $darkColor)? $lightColor : $darkColor;
				$prevFracGroup=$fracGroup;
			}
			$techRepGroup=0 unless $techRepGroup;
			#my $obsCode="$anaID:$targetPos";
			my $fracGroupStrg=($fracGroup)? "[Sample. #$fracGroup]" : '';
			my $techRepStrg=($techRepGroup)? "[Bio. repl. # $techRepGroup]" : '';
			my $labelStrg=($designInfo{'ANALYSIS'}{$anaID}{'HAS_QUANTIF'})? $designInfo{'OBSERVATION'}{$obsCode}{'NAME'} : 'Label-free';
			print "<TR bgcolor=$bgColor class=\"list\"><TH align=\"left\" nowrap>&nbsp;($designInfo{OBSERVATION}{$obsCode}{BIOSAMPLE})&nbsp;$designInfo{ANALYSIS}{$anaID}{HIERARCHY}&nbsp;>&nbsp;$labelStrg&nbsp;</TH><TH align=\"left\" nowrap>&nbsp;$fracGroupStrg&nbsp;$techRepStrg&nbsp;</TH>";
			print qq
|<TD><INPUT type="checkbox" name="condObs:$expConditionID" value="$obsID" onclick="updateObsCheckStatus('condObs:$expConditionID',this)" checked ></TD>
|;
			if ($designInfo{'ANALYSIS'}{$anaID}{'HAS_QUANTIF'}) {
				print qq
|<TD><SELECT name="anaXICQuanti:$obsID" id="anaXICQuanti:$obsID" data-state="$expConditionID" onchange="propagateQuantifSelection(this)" onmouseover="preparePopup(this.value)" onmouseout="popout()"><OPTION value="">-= Select =-</OPTION>
|;
				#my $firstOpt=1;
				foreach my $labelType (sort{$designInfo{'ANALYSIS'}{$anaID}{'LABEL'}{$a} cmp $designInfo{'ANALYSIS'}{$anaID}{'LABEL'}{$b}} keys %{$designInfo{'ANALYSIS'}{$anaID}{'LABEL'}}) {
					if ($labelType eq 'FREE') {print "<OPTGROUP label=\"Label-free:\">";} else {print "<OPTGROUP label=\"$labelType:\">";}
					foreach my $quantiID (sort{$a<=>$b} keys %{$designInfo{'ANALYSIS'}{$anaID}{'LABEL'}{$labelType}}) {
#print "<BR>*Q=$quantiID:  ($obsCode) '$designInfo{'OBSERVATION'}{$obsCode}{'ID_QUANTIFICATION'}{$quantiID}'<BR>\n";
						if (defined $designInfo{'OBSERVATION'}{$obsCode}{'ID_QUANTIFICATION'}{$quantiID}) {
							my $trueTargetPos=$designInfo{'OBSERVATION'}{$obsCode}{'ID_QUANTIFICATION'}{$quantiID};
							#my $optValue="$obsID:$quantiID:$anaID:$trueTargetPos:$fracGroup:$techRepGroup"; # true channel pos used in pep quantif
							#$optValue.=":$fracGroup" if $fracGroup; # fracGroup added for readility by human
							#my $selStrg=($firstOpt)? ' selected' : '';
							print qq
|<OPTION name="$labelType" id="opt:$obsID:$quantiID" value="$obsID:$quantiID:$anaID:$trueTargetPos:$fracGroup:$techRepGroup" data-quantif_method="$designInfo{QUANTIFICATION}{$quantiID}{METHOD}">$designInfo{QUANTIFICATION}{$quantiID}{NAME}</OPTION>
<SCRIPT type="text/javascript">popupTexts['$obsID:$quantiID:$anaID:$trueTargetPos:$fracGroup:$techRepGroup']='$designInfo{QUANTIFICATION}{$quantiID}{ANNOT_STRING}';</SCRIPT>
|;
							#$firstOpt=0;
						}
					}
					print "</OPTGROUP>\n";
				}
				print "</SELECT></TD>\n";
			}
			else {print "<TH align=\"left\" width=100% colspan=2>&nbsp;None found  (SSPA only!)&nbsp;<INPUT type=\"hidden\" name=\"anaXICQuanti:$obsID\" value=\"$obsID:0:$anaID:0:$fracGroup:$techRepGroup\"></TH>";} # no pepQuantif, no trueTargetPos
			print "</TR>\n";
		}
	}
	print "</TABLE>\n</DIV>\n";
}
print qq
|</FORM>
</CENTER>
<BR><BR>
<DIV id="divDescription" class="clDescriptionCont"></DIV>
<SCRIPT type="text/javascript">setPopup();</SCRIPT>
</BODY>
<HTML>
|;



#########################<<< SUBROUTINE >>>######################
sub launchQuantifications {

	####<Starting HTML>####
	print header(-'content-encoding'=>'no');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Launching Quantifications</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER><FONT class="title1">Launching Quantification(s)</FONT></CENTER>
<BR><BR><BR>
|;

	my @quantItemList=param('anaList'); # not defined for DESIGN

	my ($projectID,$experimentID); # defined later for Design. needed for nav frame refresh

	my @selectedConds; # for Design only
	my $usedCondNum=0;
	my $dbh=&promsConfig::dbConnect;
	($experimentID)=$dbh->selectrow_array("SELECT ID_EXPERIMENT FROM DESIGN WHERE ID_DESIGN=$designID");
	$projectID=&promsMod::getProjectID($dbh,$experimentID,'experiment');
	#$dbh->disconnect;

	
	my $labeling=param('labeling');
	my ($algoType,$ratioType,$refType)=split(':',param('algoType')); # $ratioType & $refType undef for SSPA
	$ratioType='' unless $ratioType;
	$refType='' unless $refType;
	my $pepCharge=param('pepCharge') || 'all';
	my $pepSrc=($pepCharge eq 'best')? 'best' : (param('pepSource'))? param('pepSource') : 'all'; # only for SILAC
	my $minInfRatios=param('minInfRatios') || 0;
	my $ptmScope=param('ptmScope') || 'protein';
	my (@ptmScopeModifIDs,%ptmScopeHash);
	if ($ptmScope =~ /PTM/) { # only for DB: protein by PTM / site / free-res quantification
		if ($ptmScope eq 'create_multiPTM') {@ptmScopeModifIDs=(-1);} # modID=-1 for free residue (0 not valid in DB!!!!)
		else {@ptmScopeModifIDs=sort{abs($a)<=>abs($b)} param('multiPtmScope');} # sort is critical: Modif_rank is based on increasing modifID
		%ptmScopeHash = map {$_=>1} @ptmScopeModifIDs;
	}
	my $phosphoID=param('phosphoID');
	
	foreach my $nbCond (1..param('nbMaxCond')){ # Information to compute the ratios
		my $condID=param("state_$nbCond");
		push @selectedConds,$condID if ($condID && $condID ne 'Select');
	}
	$usedCondNum=scalar @selectedConds; # for Design only

	###<%(REF|TEST)%>###
	my $quantifBaseName=param('quantifName');
	my %condNames;
	if ($quantifBaseName=~/%(REF[^%]*|TEST)%/) {
		my $condIDs=join(',',@selectedConds);
		my $sthConds=$dbh->prepare("SELECT ID_EXPCONDITION,NAME FROM EXPCONDITION WHERE ID_EXPCONDITION IN ($condIDs)");
		$sthConds->execute;
		while (my ($condID,$condName)=$sthConds->fetchrow_array) {
			$condNames{$condID}=$condName;
		}
		$sthConds->finish;
	}
	$dbh->disconnect;
	
	###<Multi-quantif>###
	my $numJobs=0; # default
	my %job2Conds;
	if (param('multiQuantif') && $ratioType ne 'SuperRatio') { # only for 'All vs State1' (SuperRatio just to be safe)
		#$numJobs=($refType eq 'singleRef')? $usedCondNum-1 : $usedCondNum*($usedCondNum-1)/2;
		foreach my $i (0..$#selectedConds-1) {
			foreach my $j ($i+1..$#selectedConds) {
				$numJobs++;
				@{$job2Conds{$numJobs}}=($selectedConds[$i],$selectedConds[$j]); # REF,TEST
			}
			last if $refType eq 'singleRef';
		}
	}
	else {$numJobs=1;}
	#}

	mkdir "$promsPath{tmp}/quantification" unless -e "$promsPath{tmp}/quantification";
	my $currentQuantifPath="$promsPath{tmp}/quantification/current";
	mkdir $currentQuantifPath unless -e $currentQuantifPath;

	#my $quantifBaseName=param('quantifName');
	my $refConfName; # for multi-Q only

	my $quantifDate=strftime("%Y%m%d%H%M%S",localtime);
	while (-e "$promsPath{tmp}/quantification/$quantifDate") { # to prevent multi-user launch collision
		sleep 2;
		$quantifDate=strftime("%Y%m%d%H%M%S",localtime);
	}

	print '<BR><FONT class="title3">Launching quantification process';
	if ($numJobs > 1) {print "es [Master job=#$quantifDate]:\n";}
	else {print " [job #$quantifDate]...";}

	my $fullQuantifType=($algoType=~/^(MSstats|SSPA|TDA|DIA)$/)? 'DESIGN:'.$algoType : 'DESIGN';
	my @jobList;
	foreach my $jobPos (1..$numJobs) {

		sleep 1 if $jobPos > 1; # wait 1 sec between jobs
		my $jobDir=$quantifDate;
		$jobDir.=".$jobPos" if $numJobs > 1;

		my $quantifPath="$promsPath{tmp}/quantification/$jobDir";
		mkdir $quantifPath || die "ERROR detected: $!";

		print "<BR>&nbsp;&nbsp;-$jobPos/$numJobs" if $numJobs > 1;

		my $quantifName=$quantifBaseName;
		my @states;
		if ($numJobs==1) {
			@states=@selectedConds;
			$quantifName=~s/%REF[^%]*%/$condNames{$states[0]}/gi;
			if (scalar @states == 2) {
				$quantifName=~s/%TEST%/$condNames{$states[1]}/gi;
			}
			else {
				$quantifName=~s/%TEST%/All-test-states/gi;
			}
			$quantifName=~s/%#%/$jobPos/g;
		}
		else { # multi-quantif
			#@states=($selectedConds[0],$selectedConds[$jobPos]);
			@states=@{$job2Conds{$jobPos}};			
			#if ($quantifName=~/%(TEST|REF[^%]*)%/i) {
			#	my $dbh=&promsConfig::dbConnect;
			#	my $sthCond=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
			#	if ($quantifName=~/%REF[^%]*%/i) {
			#		unless ($refConfName) { # do only once
			#			$sthCond->execute($states[0]);
			#			($refConfName)=$sthCond->fetchrow_array;
			#		}
			#		$quantifName=~s/%REF[^%]*%/$refConfName/gi;
			#	}
			#	if ($quantifName=~/%TEST%/i) {
			#		$sthCond->execute($states[1]);
			#		my ($testCondName)=$sthCond->fetchrow_array;
			#		$quantifName=~s/%TEST%/$testCondName/gi;
			#	}
			#	$sthCond->finish;
			#	$dbh->disconnect;
			#}
			$quantifName=~s/%REF[^%]*%/$condNames{$states[0]}/gi;
			$quantifName=~s/%TEST%/$condNames{$states[1]}/gi;
			$quantifName=~s/%#%/$jobPos/g;
		}
		my $usedStateNum=scalar @states;

		###<jobs info & data >###
		open (INFO,">$quantifPath/quantif_info.txt"); # item valueR valueDB
		print INFO "USER=$userID\n";
		print INFO "TYPE=DESIGN\n";
#		my $labeling=param('labeling');
#		##($algoType,my $refType)=split(':',param('algoType')); # $refType defined only for SimpleRatio & MSstats (multiRef or singleRef)
#		##my $ratioType=($algoType=~/S\w+Ratio/)? $algoType : ($algoType eq 'MSstats')? 'SimpleRatio' : ($algoType eq 'SSPA')? 'SSPA' : 'Ratio'; # ratioType is Ratio also for TnPQ!!!
#		my ($algoType,$ratioType,$refType)=split(':',param('algoType')); # $ratioType & $refType undef for SSPA
#		$ratioType='' unless $ratioType;
#		$refType='' unless $refType;
##$algoType=$ratioType='SuperRatio';
##$refType='singleRef';
#		#my $numChannels=param('numChannels'); # not for DESIGN
#		#my $maxReplicates=param('maxReplicates'); # not for DESIGN
#		my $pepCharge=param('pepCharge') || 'all';
#		my $pepSrc=($pepCharge eq 'best')? 'best' : (param('pepSource'))? param('pepSource') : 'all'; # only for SILAC
#		my $minInfRatios=param('minInfRatios') || 0;
		print INFO "PARAMETERS:\n";
		print INFO "VISIBILITY\t\t",param('quantifVisibility'),"\n" if param('quantifVisibility');
		print INFO "LABEL\t\t$labeling\n"; # only for DB not R (same as $quantifType for internal quantif)
		print INFO "QUANTIF_NAME\t\t$quantifName\n"; # only for DB
		if ($ptmScope =~ /protein/) {
			if ($algoType =~ /^(PEP_INTENSITY|DIA|TDA)/) {
				print INFO "TOP_N\t\t",param('topN'),"\n" if param('topN'); # only for DB not R (topN is undef if PTM quantif)
				print INFO "MATCH_PEP\t\t",param('matchingPep'),"\n" if param('matchingPep'); # only for DB not R
			}
			if ($ptmScope eq 'proteinPTM') {
				print INFO "PROTEIN_PTM\t";
				foreach my $modID (@ptmScopeModifIDs) {
					print INFO "\t#$modID";
					if (param("context_$modID")) {
						my @contexts=split(/\s+/,param("context_$modID"));
						print INFO '&'.join('&',@contexts);
					}
				}
				print INFO "\n";
			}
		}
		else { # site, free res
			#@ptmScopeModifIDs=(sort{abs($a)<=>abs($b)} param('multiPtmScope')); # sort is critical: Modif_rank is based on increasing modifID
			print INFO "PTM_SCOPE\t\t",join("\t",@ptmScopeModifIDs),"\n";
			print INFO "KEEP_SEQ_CONTEXT\t\t",(param('keepSeqContext') || 0),"\n";
			my ($createPTMStrg,$hasMatchedPTMs,$hasPhosphoID)=('',0,0);
			#my $phosphoID=param('phosphoID');
			#my $modifRank=0;
			foreach my $modID (@ptmScopeModifIDs) {
				#$modifRank++;
				if ($modID == -1) { # Free residue
					my @freeResContexts=split(/\s+/,param('context_freeRes'));
					print INFO "SITE_CONTEXT\t\t#$modID&",join('&',@freeResContexts),"\n";
					# print INFO "FREE_RESIDUES\t\t",join("\t",param('freeRes')),"\n";
					# my $freeResSpecif=param('freeResSpecif');
					# my %newResidues;
					# foreach my $resValue (param('freeRes')) {
					# 	my $qResValue=quotemeta($resValue);
					# 	$newResidues{$resValue}=1 if $freeResSpecif !~ /$qResValue/;
					# }
					# if (scalar %newResidues) { # New residue(s) used!
					# 	my @usedResidues;
					# 	foreach my $res (@freeResList) {
					# 		my $resValue=($res=~/Nter/)? '-' : ($res=~/Cter/)? '+' : $res;
					# 		my $qResValue=quotemeta($resValue);
					# 		push @usedResidues,$resValue if ($freeResSpecif =~ /$qResValue/ || $newResidues{$resValue});
					# 	}
					# 	my $newSpecif=join(',',@usedResidues);
					# 	my $dbh=&promsConfig::dbConnect;
					# 	$dbh->do("UPDATE MODIFICATION SET SPECIFICITY='$newSpecif' WHERE ID_MODIFICATION=-1"); # update list of usable free residues
					# 	$dbh->commit;
					# 	$dbh->disconnect;
					# }
				}
				else { # normal PTMs
					$hasMatchedPTMs++;
					$hasPhosphoID=1 if $modID==$phosphoID; # phospho has "Other" is not considered
				}
			}
			if ($hasMatchedPTMs) { # normal PTMs
				print INFO "PTM_POS\t\t";
				my $ptmProbSoftCode=param('ptmProbSoftCode');
				if (($hasPhosphoID || $ptmProbSoftCode eq 'MQ') && $algoType !~ /^TDA/) { # phospho-quantif w PRS or anyPTM-quantif with MaxQuant (Not for DIA:MSstats)
					print INFO "$ptmProbSoftCode:",param('okProb'),"\t",param('badProb'),"\n"; # PRS or MQ
				}
				else { # non-phospho PTM quantif
					# my $ambigPos=(param('ambiguousPos') && $hasMatchedPTMs==1 && !$createPTMStrg)? 'ambiguous' : 'valid';
					# print INFO "$ambigPos\n";
					print INFO "valid\n"; # for now, only MQ has any PTM prob (will change with PTMRS)
				}
				#my $ambigPos=param('ambiguousPos') || 0;
				#print INFO "AMBIGUOUS_QPTM_POS\t\t$ambigPos\n";
			}
			if (param('refQuantifType')) {
				print INFO "INTRA_PROT_REF\t\t";
				if (param('refQuantifType')==-1) {print INFO "-1\n";} # -1 (current dataset)
				elsif (param('refQuantif')) {
					print INFO param('refQuantif'),"\n"; # quantifID_ratioPos:refCondID:testCondID
				}
			}
		}
		my $ptmFilter;
		if (abs(param('pepPTM'))==2) { # new custom selection options
			if (param('pepPTM')==2) { # allow selected
				$ptmFilter='2:#'.join(',#',param('customPTM'));
			}
			else { # -2 exclude selected & sequence -> convert into list of selected!!!
				my (%notAllowed,@allowed);
				foreach my $ptmID (param('customPTM')) {$notAllowed{$ptmID}=1;}
				foreach my $ptmID (split(',',param('listPTM'))) {push @allowed,$ptmID if (!$notAllowed{$ptmID} && !$ptmScopeHash{$ptmID});} # Quantified PTM(s) not listed in allowed PTMs
				$ptmFilter='-2:#'.join(',#',@allowed);
			}
		}
		else {$ptmFilter=param('pepPTM');} # old options

		if ($algoType eq 'SSPA') {
			print INFO "PEPTIDES\t\t",param('pepSpecifity'),"\t",param('pepCleavage'),"\t$ptmFilter\t",param('pepFocus'),"\t";
			print INFO (param('pepRescued') && param('pepFocus') ne 'sp_count')? "0\n" : "1\n"; # flag for skipping MBWR data 0=>MBWR, 1=>No MBWR (back compatibility) (always "1" if sp_count)
			print INFO "HIDDEN_PROT\t\t",(param('protHidden') || 0),"\n";
		}
		else {
			print INFO "PEPTIDES\t\t",param('pepSpecifity'),"\t",param('pepCleavage'),"\t$ptmFilter\t$pepCharge\t$pepSrc","\n"; # only for DB ***Leave pepCharge & pepSrc for SWATH for compatibility with showProtQuantif***
			print INFO "MG_SHARED_PEP\t\t",param('MGSharedPep'),"\n" if ($ratioType eq 'Abundance' && $ptmScope eq 'protein'); # only for protein abundance (for now)
			if ($ptmScope ne 'protein' && param('refQuantifType') && param('refQuantifType')==-1) {
				print INFO "PEPTIDES_REF\t\t";
				if (param('manualRefPep')) {
					print INFO join("\t",('manual',param('manualPep'))),"\n"; # manual\tprotID:pep1ID,pep2ID,...,pepNID
				}
				else {
					print INFO "current\n"; # same as PEPTIDES
				}
			}
		}
		#<Protein selection (since v1.7.2)
		if (param('protSelection')) {
			print INFO "PROTEINS\t\t",param('protSelType'),"\t#",param('protSelection'),"\n";
		}
		print INFO "EXCL_CONTAMINANTS\t\t",(param('exclContaminants') || 0),"\n";
		print INFO "ALGO_TYPE\t\t$algoType\n"; # if $ratioType ne 'Ratio'; # needed for DB (not needed if Ratio) / Overwritten by launchQuantif if PEP_xxx incompatibilty
		if ($algoType ne 'SSPA') {
			my $ratioTypeCode=($ratioType eq 'Abundance')? 'None' : $ratioType;
			print INFO "RATIO_TYPE\t\t$ratioTypeCode\n";
			print INFO "NUM_TRUE_USED\t\t1\n" if ($algoType eq 'PEP_INTENSITY' && param('trackMsmsPep'));
			if ($ratioType=~/S\w+Ratio/) {
				print INFO "SINGLE_REF\t\t1\n" if ($refType && $refType eq 'singleRef');
			}
			elsif ($ratioType eq 'Abundance') { # <====== TEMP (parameter name to be changed)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				print INFO "\tnormalization.only\tyes\n";
#print INFO "\tresidual.variability\ttechnical\n"; # useless but still required
#print INFO "\tpAdj.method\tfdr\n"; # useless but still required
				my $abunMeasStrg=join("\t",param('abundMeasures'));
				foreach my $meas (param('abundMeasures')) {
					if ($meas eq 'MY_LFQ') {
						#my $minPepRatioLFQ=($ptmScope eq 'multiPTM')? 1 : 2; # Hard-coded but could become a user-defined parameter
						#print INFO "NUM_LFQ_RATIOS\t\t",param('minRatioLFQ'),"\n";
						#print INFO "STAB_LARGE_RATIOS\t\t",param('stabRatioLFQ'),"\n" if param('stabRatioLFQ');
						my $lfqParamStrg=','.param('minRatioLFQ').',';
						$lfqParamStrg.=(param('stabRatioLFQ'))? param('stabRatioLFQ') : 0;
						$abunMeasStrg=~s/MY_LFQ/MY_LFQ$lfqParamStrg/;
					}
					elsif ($meas eq 'ABSOLUTE') {
						my $absPath = "$quantifPath/absolute";
						mkdir $absPath;
						my $aLFQInfoFile = "$absPath/aLFQ_info.txt";
						my $confLevel = param('conf_level') || 0.95;
						my $pqi      = param('pqi');
						my $outType  = param('out_units');
						my $absModel = param('abs_model');
						my $pqiNorm  = param('pqi_normalization');
						my $qtyUnit  = param('units');
						my $absParamStrg = "$pqi,$absModel,$pqiNorm,$outType";

						my ($peptideTopX, $peptideStrictness);
						if ($pqi eq 'aLFQ_TOP') {
							$peptideTopX = &promsMod::cleanNumericalParameters(param('pep_topx'));
							$peptideStrictness = &promsMod::cleanParameters(param('pep_strictness'));
							$absParamStrg .= ",$peptideTopX,$peptideStrictness";
						}
						$absParamStrg .= ",$qtyUnit" if ($qtyUnit);
						my @absMeasures = ($absParamStrg);
						if ($outType eq "composition") {
							push @absMeasures, ("MOL_PERCENT", "MASS_PERCENT");
						}
						else {
							# Double check that a file with quantities and their unit are given to compute amounts
							if ((!param('total_conc_file') && !param('concentration_file')) || !$qtyUnit) {
								warn("There is no file with protein quantities, cannot compute amounts !");
								exit;
							}
							if ($qtyUnit =~ /conc/) {  # Is the quantity provided a concentration ?
								push @absMeasures, ("MOL_PERCENT", "MOL_CONC", "MASS_PERCENT", "MASS_CONC");
							} else {  # Or is it an actual amount ?
								push @absMeasures, ("MOL_PERCENT", "MOL", "MASS_PERCENT", "MASS");
							}
						}
						my $absMeasAnnot = join(';', @absMeasures);
						$abunMeasStrg =~ s/ABSOLUTE/$absMeasAnnot/;

						&processAbsQuantifInfo($aLFQInfoFile, $algoType, $designID);
						print INFO "CONFIDENCE_LEVEL\t\t$confLevel\n";
					}
				}
				print INFO "ABUND_MEASURES\t\t$abunMeasStrg\n";
			}
			#<Bias correction
			my $normMethod=param('biasCorrect');
			my $biasCorrect=($normMethod eq 'none.none')? 'FALSE' : 'TRUE';
			$biasCorrect.=','.param('biasCorrLevel') if ($biasCorrect eq 'TRUE' && $ratioType eq 'Abundance'); # Specify bias correction level(s). ONLY for ABUNDANCE
			my $okSharedPep=0;
			print INFO "BIAS_CORRECTION\t\t$biasCorrect"; # for DB only
			if (param('refProt') ne 'manual' && ($normMethod=~/REF_PROT|globalStandards/ || ($algoType=~/^(PEP_|TDA|DIA)/ && $normMethod !~/quantile/i))) { # all algos except SSPA
				# nothing written for 'all'
				if (param('refProt')=~/^\d+$/) {print INFO "\t#",param('refProt'),"\t",param('refProtSelType');} # listID & use/exclude
				if (param('normSharedPep') && $algoType ne 'MSstats') {
					print INFO "\tsharedPep";
					$okSharedPep=1;
				}
			}
			print INFO "\n";
			print INFO "\tcommon.peptides\tyes\n" if $okSharedPep; # R params (use only peptides shared across all states for normalization)
			if ($ratioType=~/S\w+Ratio|Abundance/) { # S*Ratio & MSstats
				my ($usedMethod)=($normMethod eq 'none.none' && $algoType eq 'MSstats')? 'FALSE' : $normMethod;
				$usedMethod=~s/scale/none/ if ($usedMethod=~/scale/ && param('refProt') eq 'manual' && scalar param('manualPep') == 1); # no scale if only 1 peptide ion for normalization
				print INFO "NORMALIZATION_METHOD\tnormalization.method\t$usedMethod\n";

				if ($ratioType=~/S\w+Ratio/) { # since SuperRatio
					print INFO "PEPTIDES_NORM\t\t",join("\t",('manual',param('manualPep'))),"\n" if (param('refProt') eq 'manual' && param('manualPep')); # manual\tprotID:pep1ID,pep2ID,...,pepNID
					#<SWATH with MSstats (list of ratios to be computed: C2/C1 C3/C1 ... Cn/C1 C3/C2 ... Cn/C2 ...)
					if ($algoType eq 'MSstats') {
						#>Ratios computed
						print INFO "\tcontrasts.matrix\t"; # R only
						foreach my $i (0..$#states-1) {
							foreach my $j ($i+1..$#states) {
								print INFO ";" if $j > 1;
								print INFO 'State'.($j+1).'/State'.($i+1); # positions not indexes
							}
							last if $refType eq 'singleRef';
						}
						print INFO "\n";
						#>CPU to use
						my %clusterInfo=&promsConfig::getClusterInfo;
						#my $clusters=($clusterInfo{'on'} && $clusterInfo{'maxCPUs'})? $clusterInfo{'maxCPUs'} : 1;
						my $clusters=($clusterInfo{'on'})? 8 : 1; # to match launchQuantification.pl
						print INFO "\tclusters\t$clusters"; # R only
						print INFO "\n";
					}
					else {
						#print INFO "\tsavegraph\tTRUE\n"; # only in R (for selected proteins) ?
						#print INFO "\tdisplaygraph\tTRUE\n"; # only in R
						#my $labelingR=($labeling eq 'FREE')? 'LabelFree' : ($usedStateNum==2)? 'SILAC' : ($ratioType eq 'SuperRatio')? 'SUPERSILAC' : 'multiSILAC'; # ***multiSILAC does not work if only 2 states! (25/02/15)***
						#my $labelingR=($labeling eq 'FREE')? 'LabelFree' : ($usedStateNum==2 && $labeling ne 'ITRAQ')? 'SILAC' : ($ratioType eq 'SuperRatio')? 'SUPERSILAC' : 'multiSILAC'; # ***multiSILAC does not work if only 2 states! (25/02/15)***
					}
				}
			}
		}
		print INFO "ID_DESIGN\t\t$designID\n";
		print INFO "STATES\t\t",join("\t",@states),"\n";

		#foreach my $quantiIDAnaIDCondID (@quantItemList){# One XIC-quanti could contain multi ANA quantified
		#	my @condInfo=split(/:/,$quantiIDAnaIDCondID);
		#	print INFO "QUANTITOCOND\t\t".param("anaXICQuanti:$quantiIDAnaIDCondID")."\n" if $states{$condInfo[1]};
		#}

		# QUANTITOCOND  (obsID quantifID fracGroup)
		#print INFO "ANALYSES/OBSERVATIONS & PARENT QUANTIFICATIONS:\n";

		my $residualVar=param('residualVar'); # defined only for $algoType=~/S\w+Ratio/
		my $checkBioRep=($algoType=~/^(PEP_|TDA|DIA)/ && $ratioType ne 'Abundance' && $residualVar eq 'auto')? 1 : 0;
		my $numBioRepOK=0;
		my $isState1=1;
		my $okSuperRatio=1; # flag for compatibility with SuperRatio
		foreach my $condID (@states) {
			my $filterQuantifs=0;
			my %usedTestParQuantifs; # multi-quantif with labeling only
			if (param('multiQuantif') && $labeling ne 'FREE' && $isState1) { # Filter for unmatched obs
				$filterQuantifs=1; # multi-quantif with labeling only
				foreach my $obsID (param("condObs:$states[1]")) {
					my ($obsID,$parQuantifID,$anaID,$targetPos,$fracGroup,$bioRep)=split(':',param("anaXICQuanti:$obsID"));
					$usedTestParQuantifs{"$parQuantifID:$anaID"}=1;
				}
			}
			#my ($bioRepCount,$prevTechRepGroup,$prevFracGroup)=(0,0,0);
			my $bioRepCount=0;
			my (%usedFracGroups,%usedTechRepGroups);
			my @condObsData;
			foreach my $obsID (param("condObs:$condID")) {
				if ($filterQuantifs) { # multi-quantif with labeling only
					foreach my $quantifCode (keys %usedTestParQuantifs) {
						if (param("anaXICQuanti:$obsID")=~/^\d+:$quantifCode:/) {
							push @condObsData,param("anaXICQuanti:$obsID");
							last;
						}
					}
				}
				else {push @condObsData,param("anaXICQuanti:$obsID");}
				if ($checkBioRep) {
					my ($fracGroup,$techRepGroup)=(split(':',param("anaXICQuanti:$obsID")))[4,5];
					#$bioRepCount++ if ((!$fracGroup || $fracGroup != $prevFracGroup) && (!$techRepGroup || $techRepGroup != $prevTechRepGroup));
					#$prevFracGroup=$fracGroup if $fracGroup;
					#$prevTechRepGroup=$techRepGroup if $techRepGroup;
					$bioRepCount++ if ((!$fracGroup || !$usedFracGroups{$fracGroup}) && (!$techRepGroup || !$usedTechRepGroups{$techRepGroup}));
					$usedFracGroups{$fracGroup}=1 if $fracGroup;
					$usedTechRepGroups{$techRepGroup}=1 if $techRepGroup;
				}
			}
			print INFO "QUANTITOCOND_$condID\t\t",join("\t",@condObsData),"\n";
			$numBioRepOK++ if $bioRepCount >= 2; # at least 2 bioRep / state (for $residualVar eq 'auto' only)

			##>Checking SuperRatio compatibility (only 1 bioRep/anaID for State1)
			if ($ratioType eq 'SuperRatio' && $usedStateNum > 2 && $isState1) {
				my %anaBioRep;
				foreach my $obsData (@condObsData) {
					my ($obsID,$parQuantifID,$anaID,$targetPos,$fracGroup,$bioRep)=split(/:/,$obsData);
					$anaBioRep{$anaID}{$bioRep}=1;
					if (scalar keys %{$anaBioRep{$anaID}} > 1) { # More than 1 bioRep/anaID => change SUPERRATIO to LABELFREE
						$okSuperRatio=0;
						last;
					}
				}
			}
			$isState1=0;
		}
		if ($algoType=~/^(PEP_|TDA|DIA)/) { # Not Ratio,TnPQ,SSPA,MSstats
			#my $Rdesign=($labeling eq 'FREE')? 'LABELFREE' : ($ratioType eq 'SuperRatio')? 'SUPERRATIO' : 'LABELED'; # ***New algo (17/08/17)***
			my $Rdesign=($algoType eq 'PEP_RATIO' && $okSuperRatio)? 'PEP_RATIO' : 'PEP_INTENSITY'; # will overwrite ALGO_TYPE above if necessary (in launchQuantifications.cgi)
			print INFO "\tdesign\t$Rdesign\n"; # only in R

			if ($ratioType ne 'Abundance') {
				if ($checkBioRep) {
					$residualVar=($numBioRepOK==$usedStateNum)? 'biological' : 'technical';
				}
				print INFO "RESIDUAL_VAR\tresidual.variability\t$residualVar\n";
				print INFO "FDR_CONTROL\tpAdj.method\t",param('fdrControl'),"\n";
			}
		}
		elsif ($algoType eq 'SSPA') {
			print INFO "SPECIF_THRESHOLDS\t\t",param('sspaBestDelta'),"\t",param('sspaLogic'),"\t",param('sspaPvalue'),"\n";
			print INFO "FDR_CONTROL\t\tfdr\n"; # hard-coded in analysisCounting.R for now
		}
		print INFO "MIN_INF_RATIOS\t\t$minInfRatios\n" if ($algoType !~ /MSstats|SSPA/ && $ratioType ne 'Abundance'); # if $usedStateNum==2; # +/- infinite ratio switch (DB only)

		open(FLAG,">$currentQuantifPath/$jobDir\_request.flag"); # flag file used by watchQuantif
		print FLAG "#";
		close FLAG;

		close INFO;

		#my $fullQuantifType=($algoType=~/^(MSstats|SSPA|TDA|DIA)$/)? 'DESIGN:'.$algoType : 'DESIGN';
		if ($numJobs==1) {
			###<Forking to launch quantifications>###
			my $childPid = fork;
			unless ($childPid) { # child here
				#>Disconnecting from server
				open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
				open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
				#open STDERR, ">>$promsPath{logs}/launchQuantification.log";
				system "./launchQuantifications.pl single $ENV{REMOTE_USER} $jobDir $fullQuantifType";
#system "./launchQuantifications.pl single $ENV{REMOTE_USER} $jobDir $fullQuantifType 2> $promsPath{tmp}/quantification/$jobDir\_errorQuantif.txt";
				exit;
			}
		}
		else {
			push @jobList,"$jobDir#$fullQuantifType";
		}
	} # end of jobs loop
	print '<BR><BR>' if $numJobs > 1;
	print " Done.</FONT><BR>\n";

	####>Multi-job: fork & call launchQuantification.pl in multi-job mode<####
	if ($numJobs > 1) {
		my $childPid = fork;
		unless ($childPid) { # child here
			#>Disconnecting from server
			open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
			open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
			system "./launchQuantifications.pl multi $ENV{REMOTE_USER} ".join(',',@jobList);
#system "./launchQuantifications.pl multi $ENV{REMOTE_USER} ".join(',',@jobList)." 2> $promsPath{tmp}/quantification/errorMultiQuantif.txt";
			exit;
		}
	}

	##<Refresh display>##
	#my $fullQuantifType = ($algoType=~/^(MSstats|SSPA|TDA)$/)? 'DESIGN:'.$algoType : 'DESIGN';
	my $jobType = ($ratioType eq 'Abundance')? "DESIGN:$algoType:Abund" : "DESIGN:$algoType";
	&refreshHTML($projectID,$experimentID,"design:$designID",$jobType);
	exit;
}

sub resumeQuantification { # Global: $designID, %promsPath, $userID
	my $quantifID=&promsMod::cleanNumericalParameters(param('RESUME'));
	my $projectID=&promsMod::cleanNumericalParameters(param('PROJ_ID'));
	my $dbh=&promsConfig::dbConnect;
	my ($name,$quantifAnnot)=$dbh->selectrow_array("SELECT NAME,QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID");

	####<Starting HTML>####
	print header(-'content-encoding'=>'no');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Resume Quantifications</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER><FONT class="title1">Resuming Quantification "<FONT color='#DD0000'>$name</FONT>"</FONT></CENTER>
<BR><BR><BR>
|;

	####<Transfer tmp data files>####
	my $currentQuantifPath="$promsPath{tmp}/quantification/current";
	#<Find old job info
	my ($oldQuantifDate,$jobFeatures)=$dbh->selectrow_array("SELECT ID_JOB,FEATURES FROM JOB_HISTORY WHERE FEATURES LIKE '\%ID_QUANTIFICATION=$quantifID;\%'");
	my $oldQuantifPath="$promsPath{tmp}/quantification/$oldQuantifDate";
	my ($jobType)=$jobFeatures=~/TYPE=([^;]+)/;
	# my $jobStatFile=`ls -l $currentQuantifPath/$quantifID\_* | head -1`;
	# my ($oldQuantifDate)=$jobStatFile=~/$quantifID\_(.+)_/;

	
	#<Create new job dir
	my $quantifDate=strftime("%Y%m%d%H%M%S",localtime);
	my $quantifPath="$promsPath{tmp}/quantification/$quantifDate";
	mkdir $quantifPath || die "ERROR detected: $!";
	copy("$oldQuantifPath/quantif_info.txt","$quantifPath/quantif_info.txt");
	dirmove("$oldQuantifPath/quanti_$quantifID","$quantifPath/quanti_$quantifID");
	#system "echo $quantifID > $quantifPath/resume.txt"; # flag for launchQuantification.pl

	#<Clean old files & info
	$dbh->do("DELETE FROM JOB_HISTORY WHERE ID_JOB='$oldQuantifDate'");
	rmtree($oldQuantifPath);
	unlink glob("$currentQuantifPath/$quantifID\_*") if glob("$currentQuantifPath/$quantifID\_*");

	#<Set quantif as not launched
	$dbh->do("UPDATE QUANTIFICATION SET STATUS=-1,UPDATE_USER='$userID' WHERE ID_QUANTIFICATION=$quantifID");
	$dbh->commit;

	my ($experimentID)=$dbh->selectrow_array("SELECT ID_EXPERIMENT FROM DESIGN WHERE ID_DESIGN=$designID");
	$dbh->disconnect;


	print "<BR><FONT class=\"title3\">Resuming quantification process [job #$quantifDate]...";

	###<Forking to launch quantifications>###
	my $fullQuantifType=$jobType;
	$fullQuantifType=~s/:Abund//;
	$fullQuantifType=~s/:PEP_INTENSITY//;  # DESIGN(:(DIA|TDA|MSstats|SSPA)) non-abundance not resumable as of 12/01/21
	my $childPid = fork;
	unless ($childPid) { # child here
		#>Disconnecting from server
		open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
		open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
		#open STDERR, ">>$promsPath{logs}/launchQuantification.log";
		system "./launchQuantifications.pl single $ENV{REMOTE_USER} $quantifDate $fullQuantifType $quantifID";
#system "./launchQuantifications.pl single $ENV{REMOTE_USER} $quantifDate $fullQuantifType $quantifID 2> $promsPath{tmp}/quantification/$quantifID\_errorQuantif.txt";
		exit;
	}

	####<Refresh display>####
	&refreshHTML($projectID,$experimentID,"quantification:$quantifID",$jobType);
	exit;
}

sub refreshHTML {
	my ($projectID,$experimentID,$branchID,$jobType)=@_;
	print "<BR><FONT class=\"title3\"><BR>This page will refresh itself in a few seconds.</FONT>\n";

	###>Calling watch popup window<###
	sleep 3;
	print qq
|<SCRIPT type="text/javascript">
var monitorJobsWin=window.open("$promsPath{cgi}/monitorJobs.cgi?filterType=Quantification [$jobType]&filterDateNumber=1&filterDateType=DAY&filterStatus=Queued&filterStatus=Running&filterProject=$projectID",'monitorJobsWindow','width=1200,height=500,scrollbars=yes,resizable=yes');
monitorJobsWin.focus();
</SCRIPT>
|;
#exit; # DEBUG!!!!!
	sleep 5;
	print qq
|<SCRIPT type="text/javascript">
// top.promsFrame.selectedAction='summary';
parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=$branchID&ACT=experiment&EXPERIMENT=EXPERIMENT:$experimentID&ISNAVFRAME=0&VIEW=quanti";
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}


sub ajaxFetchReferenceQuantifications {
	my $labeling=&promsMod::cleanParameters(param('labeling'));
	my $algoType=(split(':',param('algoType')))[0];
	$algoType=&promsMod::cleanParameters($algoType);
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
	
	my $dbh=&promsConfig::dbConnect;
	my $sthQuantif=$dbh->prepare("SELECT Q.ID_QUANTIFICATION,Q.NAME,COUNT(ID_EXPCONDITION),GROUP_CONCAT(DISTINCT MQ.ID_MODIFICATION) AS MULTI_MOD FROM QUANTIFICATION Q
									LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
									INNER JOIN QUANTIFICATION_METHOD QM ON Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND QM.CODE='PROT_RATIO_PEP'
									INNER JOIN EXPCONDITION_QUANTIF EQ ON Q.ID_QUANTIFICATION=EQ.ID_QUANTIFICATION
									WHERE ID_DESIGN=? AND Q.ID_MODIFICATION IS NULL AND FOCUS='protein' AND STATUS=1
											AND QUANTIF_ANNOT LIKE 'LABEL=$labeling:%' AND QUANTIF_ANNOT LIKE '%ALGO_TYPE=$algoType:%'
									GROUP BY Q.ID_QUANTIFICATION HAVING MULTI_MOD IS NULL
									ORDER BY Q.NAME");
	$sthQuantif->execute($designID);
	my $quantifOptStrg='';
	my $numQuantif=0;
	while (my ($quantifID,$quantifName,$numStates)=$sthQuantif->fetchrow_array) {
		$numQuantif++;
		$quantifOptStrg.="<OPTION value=\"$quantifID\">$quantifName [$numStates states]</OPTION>\n";
	}
	print "<SELECT name=\"refQuantif\" onchange=\"applyNewRefQuantifSet(this.value)\">";
	if ($numQuantif) {print "<OPTION value=\"\">-= Select =-</OPTION>\n$quantifOptStrg";}
	else {print "<OPTION value=\"\">** None found **</OPTION>\n";}
	print "</SELECT>\n";
	
	$sthQuantif->finish;

	$dbh->disconnect;
	exit;
}

sub ajaxFetchRefProtPeptides {
	my $refQuantifID=&promsMod::cleanNumericalParameters(param('refQuantif'));
	#my $modifID=&promsMod::cleanNumericalParameters(abs(param('modID'))); # abs in case site re-creation. TODO: handle site re-creation here
	my $quantifiedModifStrg=param('modID') || 'X';
	$quantifiedModifStrg=~s/-//g; # abs in case site re-creation. TODO: handle site re-creation here
	$quantifiedModifStrg=~s/:/\|/g; # for matching later
	my $labeling=param('labeling');
	my $selAnaStrg=param('selAna'); # only if $refQuantif is -1
	
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
	my ($lightColor,$darkColor)=&promsConfig::getRowColors;
	
	my $dbh=&promsConfig::dbConnect;
	my $anaStrg;
	if ($refQuantifID eq '-1') { # current design context
		my %selAna;
		foreach my $anaID (split(':',$selAnaStrg)) {$selAna{$anaID}=1;} # remove duplicates
		$anaStrg=join(',',keys %selAna);
	}
	else { # NOT USED!!! manual peptide selection only allowed in current design context
		#my ($quantID,$ratioPos,$testCondID,$refCondID)=split(/[_:]/,$refQuantif);
		#my $sthRefAna=$dbh->prepare("SELECT DISTINCT(O.ID_ANALYSIS) FROM OBS_EXPCONDITION OE
		#							INNER JOIN OBSERVATION O ON OE.ID_OBSERVATION=O.ID_OBSERVATION
		#							INNER JOIN EXPCONDITION_QUANTIF EQ ON EQ.ID_EXPCONDITION=OE.ID_EXPCONDITION
		#							WHERE OE.ID_EXPCONDITION IN (?,?) AND EQ.ID_QUANTIFICATION=?");
		#$sthRefAna->execute($testCondID,$refCondID,$quantID);
		#my @anaList;
		#while (my ($anaID)=$sthRefAna->fetchrow_array) {push @anaList,$anaID;}
		#$sthRefAna->finish;
		#$anaStrg=join(',',@anaList);
		($anaStrg)=$dbh->selectrow_array("SELECT GROUP_CONCAT(ID_ANALYSIS SEPARATOR ',') FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$refQuantifID GROUP BY ID_QUANTIFICATION");	
	}
	
	#my ($anaStrg)=$dbh->selectrow_array("SELECT GROUP_CONCAT(DISTINCT O.ID_ANALYSIS),ID_DESIGN FROM EXPCONDITION E
	#										INNER JOIN OBS_EXPCONDITION OE ON E.ID_EXPCONDITION=OE.ID_EXPCONDITION
	#										INNER JOIN OBSERVATION O ON OE.ID_OBSERVATION=O.ID_OBSERVATION
	#										WHERE ID_DESIGN=$designID GROUP BY ID_DESIGN");
	my $sthTestProt=$dbh->prepare("SELECT DISTINCT P.ID_PROTEIN,ALIAS,PROT_DES,PROT_LENGTH,ORGANISM FROM PROTEIN P
										INNER JOIN ANALYSIS_PROTEIN AP ON P.ID_PROTEIN=AP.ID_PROTEIN
										WHERE ID_ANALYSIS IN ($anaStrg)
                                        AND ABS(AP.VISIBILITY)=2 LIMIT $MAX_PROT_MANUAL_PEP");
	my $sthRefProtPep=$dbh->prepare("SELECT P.ID_PEPTIDE,ABS(PEP_BEG),PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&') AS VMOD_CODE,CHARGE,PPA.IS_SPECIFIC,MISS_CUT
									FROM PEPTIDE P
									LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
									INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
									WHERE P.ID_ANALYSIS=PPA.ID_ANALYSIS AND P.ID_ANALYSIS IN ($anaStrg) AND PPA.ID_PROTEIN=? GROUP BY P.ID_PEPTIDE ORDER BY VMOD_CODE DESC"); # order by to fetch modified pep 1srt => needed to filter any pepSeq found modifed on another peptide
	
	my $labelModStrg;
	if ($labeling ne 'FREE') {
		my @labelModifs;
		my $sthALM=$dbh->prepare("SELECT DISTINCT M.ID_MODIFICATION FROM ANALYSIS_MODIFICATION AM,MODIFICATION M WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND AM.ID_ANALYSIS IN ($anaStrg) AND M.IS_LABEL=1");
		$sthALM->execute;
		while (my ($modID)=$sthALM->fetchrow_array) {push @labelModifs,$modID;}
		$sthALM->finish;
		$labelModStrg=join('|',@labelModifs);
	}
	my (%protInfo,%protPeptides,%pepInfo,%varMod,%varModName,%skipPepSeq,%skipPepIons);
	$sthTestProt->execute;
	while (my ($protID,$alias,$protDes,$protLength,$organism)=$sthTestProt->fetchrow_array) {
		@{$protInfo{$protID}}=($alias,$protDes,$protLength,$organism);
		$sthRefProtPep->execute($protID);
		while (my ($pepID,$pepBeg,$pepSeq,$modCode,$charge,$isSpecif,$missCut)=$sthRefProtPep->fetchrow_array) {
			if ($skipPepSeq{$pepSeq}) {
				my $vMod=$modCode || '-';
				$skipPepIons{$protID}{$pepSeq.'_'.$vMod.'_'.$charge}=1;
				next;
			}
			if ($modCode) {
				if ($modCode=~/(^|&)($quantifiedModifStrg):/) {
					$skipPepSeq{$pepSeq}=1;
					$skipPepIons{$protID}{$pepSeq.'_'.$modCode.'_'.$charge}=1;
					next;
				}
				if ($labeling ne 'FREE') {
					$modCode=~s/(^|&)($labelModStrg):[\d\.]+//g;
					$modCode=~s/^&//;
				}
				$varMod{$modCode}=' + '.&promsMod::decodeVarMod($dbh,$pepSeq,$modCode,\%varModName) unless $varMod{$modCode};
			}
			else {
				$modCode='';
				$varMod{''}='';
			}
			$isSpecif=($isSpecif)? 'Yes' : 'No';
			$missCut='No' unless $missCut;
			push @{$protPeptides{$protID}{$pepBeg}{$pepSeq}{$modCode}{$charge}},$pepID;
			@{$pepInfo{$protID}{$pepBeg}{$pepSeq}{$modCode}}=($isSpecif,$missCut);
		}
	}
	$sthTestProt->finish;
	$sthRefProtPep->finish;
	$dbh->disconnect;
	
	print "<TABLE border=0 cellspacing=0 cellpadding=2 width=100%>\n";
	my $numProt=scalar keys %protPeptides;
	my $protCount=0;
	foreach my $protID (sort{$protInfo{$a}[0] cmp $protInfo{$b}[0]} keys %protPeptides) {
		$protCount++;
		my ($alias,$protDes,$protLength,$organism)=@{$protInfo{$protID}};
		print qq
|<TR bgcolor="$darkColor"><TD class="bBorder" colspan=8><TABLE>
	<TR><TH valign=top>#$protCount.<A href="javascript:sequenceView($protID,'$anaStrg')">$alias</A>:</TH><TD bgcolor="$lightColor" width=100%>$protDes <FONT class="org">$organism</FONT> ($protLength aa)</TD></TR>
	</TABLE></TD></TR>
<TR><TH width=20>&nbsp;&nbsp;</TH><TH bgcolor="$darkColor" class="rbBorder" width=25>#</TH>
<TH bgcolor="$darkColor" class="rbBorder" width=50>&nbsp;Start&nbsp;</TH>
<TH bgcolor="$darkColor" class="rbBorder">Peptide&nbsp;&nbsp;</TH>
<TH bgcolor="$darkColor" class="rbBorder" width=50>&nbsp;Charge&nbsp;</TH>
<TH bgcolor="$darkColor" class="rbBorder" width=100>&nbsp;Proteotypic&nbsp;</TH>
<TH bgcolor="$darkColor" class="bBorder" width=60>&nbsp;Miss.&nbsp;cut&nbsp;</TH>
</TR>
|;
		my $bgColor=$lightColor;
		my $numPep=0;
		foreach my $pepBeg (sort{$a<=>$b} keys %{$protPeptides{$protID}}) {
			foreach my $pepSeq (sort{length($a)<=>length($b)} keys %{$protPeptides{$protID}{$pepBeg}}) {
				my ($disabPepSeq,$popupStrg)=($skipPepSeq{$pepSeq})? ('disabled','onmouseover="popup(\'<B>Disabled: This peptide sequence is also found with a PTM site modification.</B>\')" onmouseout="popout()"') : ('','');
				foreach my $modCode (sort{length($a)<=>length($b)} keys %{$protPeptides{$protID}{$pepBeg}{$pepSeq}}) {
					foreach my $charge (sort{$a<=>$b} keys %{$protPeptides{$protID}{$pepBeg}{$pepSeq}{$modCode}}) {
						$numPep++;
						my ($isSpecif,$missCut)=@{$pepInfo{$protID}{$pepBeg}{$pepSeq}{$modCode}};
						my $pepIdStrg=$protID.':'.join(',',@{$protPeptides{$protID}{$pepBeg}{$pepSeq}{$modCode}{$charge}});
						print qq
|<TR>
	<TD></TD>
	<TD bgcolor="$bgColor" class="rBorder" align=right>&nbsp;$numPep&nbsp;</TD>
	<TD bgcolor="$bgColor" align=right>$pepBeg</TD>
	<TH bgcolor="$bgColor" class="font11" align=left nowrap><LABEL><INPUT type="checkbox" name="manualPep" value="$pepIdStrg"$disabPepSeq $popupStrg>$pepSeq$varMod{$modCode}</LABEL>&nbsp;</TH>
	<TD bgcolor="$bgColor" align=center>$charge<SUP>+</SUP></TD>
	<TD bgcolor="$bgColor" align=center>$isSpecif</TD>
	<TD bgcolor="$bgColor" align=center>$missCut</TD>
	<TD></TD>
</TR>
|;
						$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
					}
				}
			}
		}
		if ($skipPepIons{$protID}) {
			print "<TR bgcolor=\"$darkColor\"><TH colspan=8 align=\"left\">&nbsp;",(scalar keys %{$skipPepIons{$protID}})," peptide ions were excluded</TD></TR>\n";
		}
		print "<TR><TD colspan=8>&nbsp;</TD></TR>\n" if $protCount < $numProt;
	}
	print "</TABLE>\n";
	
	exit;
}


sub processAbsQuantifInfo {
	my ($infoFile, $algoType, $designID)=@_;
	my ($fileName, $dirName) = fileparse($infoFile, qr/\.[^.]*/);

	my ($aLFQInputType, $pqi, $absQuantifModel, $qtyFile, $pqiNorm, $qtyUnit, $peptideTopX, $peptideStrictness);
	my ($hasQuantity, $isMass, $isConcentration, $confLevel);

	if ($algoType eq 'PEP_INTENSITY') {
		$aLFQInputType = "openmslfq";  # resultsPep.txt will be converted to an OpenMS like format for aLFQ
	} elsif ($algoType eq 'DIA') {
		$aLFQInputType = "openswath";
	} else {
		die "Problem recognizing the type of MS search (DDA/DIA) for absolute quantification : expected PEP_INTENSITY or DIA, got $algoType";
	}
	$confLevel       = param('conf_level') || 0.95;
	$pqi 			 = (param('pqi') eq 'aLFQ_TOP') ? 'top' : (param('pqi') eq 'aLFQ_iBAQ') ? 'iBAQ' : '';
	$absQuantifModel = param('abs_model');
	$qtyFile 		 = ($absQuantifModel eq 'proportionality') ? &promsMod::cleanParameters(param('total_conc_file')) : &promsMod::cleanParameters(param('concentration_file'));
	$pqiNorm 		 = param('pqi_normalization');
	$qtyUnit 		 = param('units');

	if ($pqi eq 'top') {
		$peptideTopX = &promsMod::cleanNumericalParameters(param('pep_topx'));
		$peptideStrictness = &promsMod::cleanParameters(param('pep_strictness'));
	}

	### Moving uploaded file (if provided)
	if ($qtyFile) {
		$hasQuantity = "TRUE";
		my $tmpFile = tmpFileName($qtyFile);
		my $userQtyFile = "${dirName}user_prot_quantities.csv";
		move($tmpFile, $userQtyFile);

		if ($qtyUnit =~ /mass/) {
			$isMass = "TRUE";
		} else {
			$isMass = "FALSE";
		}
		if ($qtyUnit =~ /conc/) {
			$isConcentration = "TRUE";
		} else {
			$isConcentration = "FALSE";
		}
	} else {
		$hasQuantity = "FALSE";
		$isMass = "FALSE";
		$isConcentration = "FALSE";
	}

	open(aLFQ_INFO, ">$infoFile");
	print aLFQ_INFO "input_type\t$aLFQInputType\n";
	print aLFQ_INFO "model\t$absQuantifModel\n";
	print aLFQ_INFO "pqi_normalization\t$pqiNorm\n";
	print aLFQ_INFO "has_quantity\t$hasQuantity\n";
	print aLFQ_INFO "is_mass\t$isMass\n";
	print aLFQ_INFO "is_concentration\t$isConcentration\n";
	print aLFQ_INFO "quantity_units\t$qtyUnit\n" if ($qtyUnit);
	print aLFQ_INFO "peptide_method\t$pqi\n";
	print aLFQ_INFO "peptide_topx\t$peptideTopX\n" if (defined($peptideTopX));
	print aLFQ_INFO "peptide_strictness\t$peptideStrictness\n" if (defined($peptideStrictness));
	print aLFQ_INFO "conf_level\t$confLevel\n";
	print aLFQ_INFO "consensus_peptides\tFALSE\n";
	print aLFQ_INFO "consensus_proteins\tFALSE\n";

	close aLFQ_INFO;
}


# TODO: handle PTM site re-creation in &ajaxFetchRefProtPeptides
####> Revision history
# 1.5.1 [FEATURE] Added SSPA specificity threholds (PP 30/04/21)
# 1.5.0 [FEATURE] PTM-enriched protein quantification and free residues, both with site context (PP 09/02/21)
# 1.4.5 [MINOR] Removed forgotten debug mode (PP 08/02/21)
# 1.4.4 [BUGFIX] Restrict data source to completed (STATUS >= 1) peptide quantification (PP 02/02/21)
# 1.4.3 [MERGE] master into ppoullet (PP 02/02/21)
# 1.4.2 [FEATURE] Multi-level bias correction & MG shared peptides rule & contanimants exclusion & resume failed Abundance quantification (PP 01/02/21)
# 1.4.1 [ENHANCEMENT] Add confidence level to quantifAnnot for absolute quantification (VL 14/01/21)
# 1.4.0 [FEATURE] Free residue quantification management & disabled MBWR for SSPA with Spectra counting (PP 23/11/20)
# 1.3.20b [ENHANCEMENT] Extend categories to abundance quantifications for job monitoring (VL 19/11/20)
# 1.3.20 [BUGFIX] Minor fix in JS for display of manual peptide selection option in protein-level correction context (PP 14/11/20)
# 1.3.19 [ENHANCEMENT] Add absolute quantification with aLFQ package (VL 13/10/20)
# 1.3.18 [FEATURE] Enabled SSPA for Sites (PP 05/10/20)
# 1.3.17 [MINOR] Added project selection when opening monitor jobs windows (VS 02/09/20)
# 1.3.16 [FEATURE] Compatibility with Spectronaut & PTMRS PTM position probabilities and structured help popup text (PP 26/08/20)
# 1.3.15 [BUGFIX] Fix and optimized broken TDA manual peptide selection (PP 20/08/20)
# 1.3.14 [UPDATE] Change in "shared peptides" parameter handling & remove "shared proteins" (PP 17/08/20)
# 1.3.13 [FEATURE] New normalization filter option "Shared peptides" for myProMS algo (PP 14/08/20)
# 1.3.12 [UX] Minor change in size of PTM selection block & minor JS bug fixes (PP 12/08/20)
# 1.3.11 [ENHANCEMENT] Optimized SQL query for checking for possibility of manual peptide selection & other minor uninitialized values fixes (PP 11/08/20)
# 1.3.10 [BUGFIX] Fix minor HTML bug preventing SSPA peptide options to be displayed (PP 26/06/20)
# 1.3.9 [FEATURE] Added "Shared proteins" option for bias correction and "Stabilize large ratios" option for LFQ (PP 05/06/20)
# 1.3.8 [FEATURE] Abundance measures are now optional & added option to track MBWR peptides in LF ratio quantif (PP 20/05/20)
# 1.3.7 [BUGFIX] Fix no requirement for analysis matching for Abundance of labeled data & changed some default settings (PP 03/04/20)
# 1.3.6 [BUGFIX] Fix hidden prot parameter being written with wrong value in quantif_info.txt (PP 24/03/20)
# 1.3.5 [ENHANCEMENT] Force BH fdr as only option and include FDR_CONTROL as parameter for SSPA & hide templating (PP 25/02/20)
# 1.3.4 [ENHANCEMENT] Added NUM_LFQ_RATIOS as hard-coded parameter (PP 05/02/20)
# 1.3.3 [BUGFIX] Activate forgotten "Any peptide" for DIA:Abundance and biais correction for DIA by myProMS (PP 29/01/20)
# 1.3.2 [FEATURE] Added use of visible/hidden proteins and peptide xic option for SSPA (PP 28/01/20)
# 1.3.1 [FEATURE] Added MBR-peptide filter option for SSPA (PP 22/01/20)
# 1.3.0 [FEATURE] Added DIA and Abundance quantifications with myProMS (PP 31/12/19) 
# 1.2.4 [FEATURE] Added option to manually select peptides for bias correction of (single!) protein sites quantification (26/12/19)
# 1.2.3 [CHANGES] Use new job monitoring window opening parameters (VS 18/11/19)
# 1.2.2 [FEATURE] Added option to extend missed cleavage exclusion to overlapping peptides & [UX] Added letter code and color to modification-site name (PP 15/11/19)
# 1.2.1 [BUGFIX] Removed duplicate declaration of some variables (PP 12/11/19)
# 1.2.0 [FEATURE] 3+ States allowed with protein-level normalization in label-free & multi-job option for "All vs All" (PP 31/10/19)
# 1.1.1 [MODIF] Switch from watchQuantification to monitorJobs script (VS 21/10/19)
# 1.1.0 [FEATURE] Added multi-modif and TDA by myProMS algo (PP 26/08/19)
# 1.0.1 [FIX] Bug due to undefined JS variable unusedVarMods (PP 08/04/19)
# 1.0.0 Forked from selAna4Quantification.cgi v2.2.0 to become design-specific & handle PTM-site creation (PP 04/02/19)

