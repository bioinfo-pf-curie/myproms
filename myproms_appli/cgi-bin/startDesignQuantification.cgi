#!/usr/local/bin/perl -w

################################################################################
# startDesignQuantification.cgi      1.3.4                                     #
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
use promsConfig;
use promsMod;
use promsQuantif;

#print header(-charset=>'UTF-8'); warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my %xicSoftware=('PD'=>'Proteome Discoverer','MCQ'=>'MassChroQ','MAS'=>'Mascot','PAR'=>'Paragon','PKV'=>'PeakView','MQ'=>'MaxQuant','OS'=>'OpenSwath','SKY'=>'Skyline','?'=>'Unknown');
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

my ($phosphoID)=$dbh->selectrow_array('SELECT ID_MODIFICATION FROM MODIFICATION WHERE UNIMOD_ACC=21');
$phosphoID=-1 unless $projectVarMods{$phosphoID}; # disable if not used in project
my $numPTMs=scalar keys %projectVarMods;

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
my ($ptmProbSoft,$ptmProbSoftCode)=('PhosphoRS','PRS'); # default

my $sthAE=$dbh->prepare("SELECT E.NAME,E.ID_EXPCONDITION,BS.NAME,O.ID_ANALYSIS,O.TARGET_POS,OE.FRACTION_GROUP,OE.TECH_REP_GROUP,O.ID_OBSERVATION,GROUP_CONCAT(M.ID_MODIFICATION ORDER BY M.ID_MODIFICATION SEPARATOR ':')
							FROM OBSERVATION O
							LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION
							LEFT JOIN BIOSAMPLE BS ON O.ID_BIOSAMPLE=BS.ID_BIOSAMPLE
							JOIN OBS_EXPCONDITION OE ON O.ID_OBSERVATION=OE.ID_OBSERVATION
							JOIN EXPCONDITION E ON OE.ID_EXPCONDITION=E.ID_EXPCONDITION
							WHERE E.ID_DESIGN=$designID GROUP BY O.ID_OBSERVATION ORDER BY E.DISPLAY_POS,E.ID_EXPCONDITION");
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
#my $sthAnaInfo=$dbh->prepare("SELECT ID_ANALYSIS,A.NAME,VALID_STATUS,A.DISPLAY_POS,S.DISPLAY_POS FROM ANALYSIS A, SAMPLE S WHERE S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS IN $anaString");
#$sthAnaInfo->execute;
my @sthGetAnaInfo=( # 2d-gel or sample
	#$dbh->prepare("SELECT ID_ANALYSIS,VALID_STATUS,'gel2d',G.NAME,'spot',SP.NAME,'analysis',A.NAME FROM ANALYSIS A,SAMPLE S,SPOT SP,GEL2D G,EXPERIMENT E,DESIGN D WHERE D.ID_EXPERIMENT=E.ID_EXPERIMENT AND E.ID_EXPERIMENT=G.ID_EXPERIMENT AND G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT AND S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS IN $anaString ORDER BY G.DISPLAY_POS,SP.NAME,S.DISPLAY_POS,A.DISPLAY_POS"),
	#$dbh->prepare("SELECT ID_ANALYSIS,VALID_STATUS,'sample',S.NAME,'analysis',A.NAME FROM ANALYSIS A,SAMPLE S,EXPERIMENT E,DESIGN D WHERE D.ID_EXPERIMENT=E.ID_EXPERIMENT AND E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SAMPLE=A.ID_SAMPLE AND ID_SPOT IS NULL AND A.ID_ANALYSIS IN $anaString ORDER BY S.DISPLAY_POS,A.DISPLAY_POS")

	$dbh->prepare("SELECT ID_ANALYSIS,A.NAME,CONCAT(G.NAME,'&nbsp;>&nbsp;',SP.NAME,'&nbsp;>&nbsp;',A.NAME) FROM ANALYSIS A,SAMPLE S,SPOT SP,GEL2D G,EXPERIMENT E,DESIGN D WHERE D.ID_EXPERIMENT=E.ID_EXPERIMENT AND E.ID_EXPERIMENT=G.ID_EXPERIMENT AND G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT AND S.ID_SAMPLE=A.ID_SAMPLE AND A.ID_ANALYSIS IN $anaString ORDER BY G.DISPLAY_POS,SP.NAME,S.DISPLAY_POS,A.DISPLAY_POS"),
	$dbh->prepare("SELECT ID_ANALYSIS,A.NAME,CONCAT(S.NAME,'&nbsp;>&nbsp;',A.NAME) FROM ANALYSIS A,SAMPLE S,EXPERIMENT E,DESIGN D WHERE D.ID_EXPERIMENT=E.ID_EXPERIMENT AND E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SAMPLE=A.ID_SAMPLE AND ID_SPOT IS NULL AND A.ID_ANALYSIS IN $anaString ORDER BY S.DISPLAY_POS,A.DISPLAY_POS")
);

#my $sthObs=$dbh->prepare("SELECT ID_ANALYSIS,VALID_STATUS,A.NAME FROM ANALYSIS A,ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND FOCUS='peptide' AND A.ID_ANALYSIS IN $anaString ORDER BY A.DISPLAY_POS,Q.NAME");

my $anaPos=0;
foreach my $sthAnaInfo (@sthGetAnaInfo) {
	$sthAnaInfo->execute;
	while (my ($anaID,$anaName,$hierarchy) = $sthAnaInfo->fetchrow_array)  {
		$designInfo{'ANALYSIS'}{$anaID}{'HIERARCHY'}=$hierarchy;
		$designInfo{'ANALYSIS'}{$anaID}{'NAME'}=$anaName;
		$designInfo{'ANALYSIS'}{$anaID}{'POS'}=++$anaPos;
		#$designInfo{'ANALYSIS'}{$anaID}{'ANA.DISPLAY_POS'}=$anaDPos;
		#$designInfo{'ANALYSIS'}{$anaID}{'SAMP.DISPLAY_POS'}=$sampDPos;
	}
	$sthAnaInfo->finish;
}

###> Get QUANTIFICATION ANALYSIS
#my ($xicIDMethod)=$dbh->selectrow_array('SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE=\'XIC\'');

#my $sthQInfo=$dbh->prepare("SELECT Q.ID_QUANTIFICATION,NAME,QUANTIF_ANNOT,STATUS,ID_ANALYSIS,UPDATE_DATE FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND FOCUS='peptide' AND ID_QUANTIFICATION_METHOD=$xicIDMethod AND ID_ANALYSIS IN $anaString");
my $sthQInfo=$dbh->prepare("SELECT Q.ID_QUANTIFICATION,Q.NAME,QUANTIF_ANNOT,STATUS,ID_ANALYSIS,UPDATE_DATE,QM.CODE,QM.DES FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ,QUANTIFICATION_METHOD QM WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND FOCUS='peptide' AND CODE != 'SIN' AND ID_ANALYSIS IN $anaString");
my $sthIsobar=$dbh->prepare("SELECT DISTINCT M.ID_MODIFICATION FROM MODIFICATION M,ANALYSIS_MODIFICATION AM WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND IS_LABEL=1 AND ID_ANALYSIS=? LIMIT 1"); #  AND AM.MODIF_TYPE='V'
my $sthAna=$dbh->prepare("SELECT FILE_FORMAT FROM ANALYSIS WHERE ID_ANALYSIS=?");
$sthQInfo->execute;
while (my ($quantiID,$quantiName,$quantifAnnot,$qStatus,$anaID,$date,$methodCode,$methodDes) = $sthQInfo->fetchrow_array) {
#$quantiName=~s/3plex.+\| /2plex (/;
	next unless $quantifAnnot;
#print "$quantiID,/$quantiName,/$qStatus,/$anaID,/'$methodCode'<BR>\n";
	$quantifMethods{$methodCode}=1;
	my ($xicSoftCode,$xicSoftVersion);
	$quantifAnnot=~s/::SOFTWARE=(\w+);?(\d|\.)*//; # remove software info for back compatibility
	if ($1) {
		$xicSoftCode=$1;
		$xicSoftVersion=$2 if $2;
	}
	elsif ($quantifAnnot=~/EXTRACTION_ALGO=/) {$xicSoftCode='MCQ';}
	else {
		$sthAna->execute($anaID);
		my ($fileFormat)=$sthAna->fetchrow_array;
		$xicSoftCode=($fileFormat=~/\.PDM\Z/)? 'PD' : ($fileFormat=~/^MASCOT/)? 'MAS' : ($fileFormat=~/PARAGON/)? 'PAR' : '?';
	}
	if ($xicSoftCode eq 'MQ') {
		$ptmProbSoft='MaxQuant';
		$ptmProbSoftCode='MQ';
	}
	my $labelType=($quantifAnnot=~/LABEL=([^:;]+)/)? $1 : 'FREE';
	$labelType=uc($labelType);
	if ($labelType=~/FREE|NONE/) {
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
		$popupInfo.=" v. $xicSoftVersion" if $xicSoftVersion;
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
			$popupInfo.="<BR><B>Reference:</B> $designInfo{ANALYSIS}{$quantifParameters{REFERENCE}}{NAME}";
		}
	}
	elsif ($methodCode eq 'DIA' && $xicSoftCode eq 'PKV') {
		my $modExclStrg=($quantifParameters{EXCLUDE_MODIFIED_PEPTIDES} eq 'NOT CHECKED' || $quantifParameters{EXCLUDE_MODIFIED_PEPTIDES}==0)? ' not' : '';
		$popupInfo.="<BR><B>Modified peptides are$modExclStrg excluded</B>";
		$popupInfo.="<BR><B>Num. transitions per peptide:</B> $quantifParameters{NB_TRANSITION_PER_PEPTIDE}";
		$popupInfo.="<BR><B>XIC extraction window:</B> $quantifParameters{XIC_EXTRACTION_WINDOW} min.";
		$quantifParameters{XIC_EXTRACTION_WINDOW}=~s/_DA/ Da/i;
		$quantifParameters{XIC_EXTRACTION_WINDOW}=~s/_PPM/ ppm/i;
		$popupInfo.="<BR><B>XIC width:</B> $quantifParameters{XIC_EXTRACTION_WINDOW}";
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
my $sthRefE=$dbh->prepare("SELECT DISTINCT D.ID_DESIGN,CONCAT(E.NAME,' > ',D.NAME) FROM EXPERIMENT E
							INNER JOIN DESIGN D ON E.ID_EXPERIMENT=D.ID_EXPERIMENT
							INNER JOIN QUANTIFICATION Q ON D.ID_DESIGN=Q.ID_DESIGN
							INNER JOIN QUANTIFICATION_METHOD QM ON Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD
							WHERE E.ID_PROJECT=$projectID AND Q.ID_MODIFICATION IS NULL AND Q.FOCUS='protein' AND STATUS=1 AND QM.CODE='PROT_RATIO_PEP' ORDER BY E.DISPLAY_POS,D.NAME");
$sthRefE->execute;
while (my($desID,$desFullName)=$sthRefE->fetchrow_array) {
	push @referenceDesign,[$desID,$desFullName];
}
$sthRefE->finish;

###>Check if manual peptide selection is possible (<= $MAX_PROT_MANUAL_PEP protein in design)
my $sthNumProt=$dbh->prepare("SELECT DISTINCT(ID_PROTEIN) FROM EXPCONDITION E
									INNER JOIN OBS_EXPCONDITION OE ON E.ID_EXPCONDITION=OE.ID_EXPCONDITION
									INNER JOIN OBSERVATION O ON OE.ID_OBSERVATION=O.ID_OBSERVATION
									INNER JOIN ANALYSIS_PROTEIN AP ON O.ID_ANALYSIS=AP.ID_ANALYSIS
									WHERE ID_DESIGN=? AND ABS(AP.VISIBILITY)=2 LIMIT ".($MAX_PROT_MANUAL_PEP+1)); # +1 to check if >$MAX_PROT_MANUAL_PEP
$sthNumProt->execute($designID);
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
			my ($res,$context)=split(';',$resCode);
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
</STYLE>
<SCRIPT type="text/javascript">
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
	recoverDefaultFormParam();
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
	//Bias correction
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
		//updateRefQuantification(selPtmID);
		updateTopN();
		
		if (createSiteRes) {
			//Apply target residue selection
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
				//document.getElementById('recreateSPAN').style.display='';
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
function updateQuantifFocus(scope) {
	//updateMultiPtmSelection();
	updateTopN();
	var myForm=document.selAnaForm;
	if (scope=='protein') {
		// Clear PTM filter spans
		document.getElementById('ptmSPAN').style.display='none';
		var phosphoSpan=document.getElementById('phosphoSPAN'); // optional
		if (phosphoSpan) phosphoSpan.style.display='none';
		// Clear reference quantif
		document.getElementById('refQuantifDIV').style.display='none';
		myForm.refQuantifType.selectedIndex=0;
		ajaxFetchReferenceQuantif(0);
		// Clear any multi-PTMs selection
		document.getElementById('multiPtmUL').style.display='none';
		var multiPtmChkBoxes=document.selAnaForm.multiPtmScope;
		if (multiPtmChkBoxes.length) {
			for (let i=0; i<multiPtmChkBoxes.length; i++) {
				multiPtmChkBoxes[i].checked=false;
				clearTargetRes(multiPtmChkBoxes[i].value);
			}
		}
		else {
			multiPtmChkBoxes.checked=false;
			clearTargetRes(multiPtmChkBoxes.value);
		}
		// Disable manual peptide selection for Bias correction
		if (myForm.refProt.selectedIndex==1) { // manual peptide
			myForm.refProt.selectedIndex=0;
			myForm.refProt.options[1].disabled=true;
			updateBiasCorrectionSet(0);
		}
	}
	else { // (ptm/phospho)SPANs "on" handled updateMultiPtmSelection()
		document.getElementById('refQuantifDIV').style.display=(document.getElementById('algoType').value.match('SSPA\|Abundance'))? 'none' : '';
		document.getElementById('multiPtmUL').style.display='';
updateMultiPtmSelection();
		// Enable manual peptide selection for Bias correction
		myForm.refProt.options[1].disabled=(myForm.algoType.value.match('(PEP_INTENSITY\|DIA\|TDA):SimpleRatio') && $numProtInDesign > $MAX_PROT_MANUAL_PEP)? false: true;
	}
	//updateMultiPtmSelection();
}
function clearTargetRes(ptmID) {
	if (ptmID >= 0) return; // matched PTM: nothing to unselect
	ptmID*=-1;
	var targetRes=document.getElementsByName('targetRes_'+ptmID);
	if (targetRes.length) {
		for (let i=0; i<targetRes.length; i++) {
			targetRes[i].checked=false;
			targetRes[i].disabled=false;
		}
	}
	else {
		targetRes.checked=false;
		targetRes.disabled=false;
	}
}
function updateMultiPtmSelection() { // former updateRefQuantification(ptmID)
	//Fetch list of selected PTM-sites
	var selPtmList={}, numSelMatched=0, numSelOther=0, multiPTMs=document.selAnaForm.multiPtmScope;
	if (!multiPTMs) {return;} // no project PTMs declared
	if (multiPTMs.length) {
		for (let i=0; i<multiPTMs.length; i++) {
			var ptmID=Math.abs(multiPTMs[i].value*1);
			var dispTargRes='none';
			if (multiPTMs[i].checked) {
				selPtmList[ptmID]=multiPTMs[i].value;
				if (multiPTMs[i].value*1 > 0) {numSelMatched++;}
				else { // Other PTMs
					numSelOther++;
					dispTargRes='';
				}
			}
			if (multiPTMs[i].value*1 < 0) {
				document.getElementById('targetResSPAN_'+ptmID).style.display=dispTargRes;
			}
		}
	}
	else {
		var ptmID=Math.abs(multiPTMs.value*1);
		var dispTargRes='none';
		if (multiPTMs.checked) {
			selPtmList[ptmID]=multiPTMs.value;
			if (multiPTMs.value*1 > 0) {numSelMatched++;}
			else {
				numSelOther++;
				dispTargRes='';
			}
		}
		if (multiPTMs.value*1 < 0) {
			document.getElementById('targetResSPAN_'+ptmID).style.display=dispTargRes;
		}
	}
	//Re-define custom PTMs disable status
	var cusPTMs=document.selAnaForm.customPTM;
	if (cusPTMs.length) {
		for (let i=0; i<cusPTMs.length; i++) {
			cusPTMs[i].disabled=(selPtmList[cusPTMs[i].value])? true : false;
		}
	}
	else {cusPTMs.disabled=(selPtmList[cusPTMs.value])? true : false;}

	//Other
	var phosphoSpan=document.getElementById('phosphoSPAN'),
		ptmSpan=document.getElementById('ptmSPAN'),
		ambigChk=document.selAnaForm.ambiguousPos;
		badProb=document.selAnaForm.badProb,
	    phosphoChecked=(selPtmList['$phosphoID'] && selPtmList['$phosphoID'] > 0)? true : false;
	if ((phosphoChecked \|\| '$ptmProbSoftCode'=='MQ') && !document.selAnaForm.algoType.value.match('MSstats')) {
		if (phosphoSpan) phosphoSpan.style.display='';
		ptmSpan.style.display='none';
		ambigChk.checked=false;
	}
	else {
		if (phosphoSpan) phosphoSpan.style.display='none';
		if (numSelMatched==1 && numSelOther==0) { // Possible only for single-modif quantif
			ptmSpan.style.display='';
		}
		else {
			ptmSpan.style.display='none';
			ambigChk.checked=false;
		}
	}
	if (badProb) { // MaxQuant or Phospho: Ambiguity allowed ONLY for single-PTM quantif
		if (numSelMatched + numSelOther > 1) {
			badProb.selectedIndex=0;
			badProb.options[1].disabled=true;
		}
		else {badProb.options[1].disabled=false;}
	}
	//Update manual peptide selection (if displayed)
	if (document.selAnaForm.refQuantifType.value==-1 && document.selAnaForm.manualRefPep.value==1) { // if (document.getElementById('manualRefPepDIV').style.display != 'none')
		ajaxFetchRefProtPeptides(1,'protCorr');
	}
	else if (document.selAnaForm.refProt.value==-1) { // if (document.getElementById('manualBiasPepDIV').style.display != 'none')
		ajaxFetchRefProtPeptides(-1,'biasCorr');
	}
}
function updateTargetResidues(ptmID,myCkbRes) { // (dis/en)able checked residue for all other modifs
	if (unusedVarMods.length <= 1) return; // only 1 modif
	for (let i=0; i<unusedVarMods.length; i++) {
		if (unusedVarMods[i]==ptmID) continue;
		var targetResCkb=document.getElementsByName('targetRes_'+unusedVarMods[i]);
		for (let j=0; j<targetResCkb.length; j++) {
			if (targetResCkb[j].value==myCkbRes.value) {
				targetResCkb[j].disabled=(myCkbRes.checked)? true : false;
				break;
			}
		}
	}
}
function updateTopN() {
	var myForm=document.selAnaForm;
	var algoType=myForm.algoType.value;
	if (myForm.labeling.value == 'FREE' && algoType.match('PEP_INTENSITY\|DIA\|TDA') && myForm.ptmScope.value == 'protein') {
		myForm.topN.disabled=false;
		myForm.matchingPep.disabled=false;
		document.getElementById('topSPAN').style.display='';
	}
	else {
		myForm.topN.disabled=true;
		myForm.matchingPep.disabled=true;
		document.getElementById('topSPAN').style.display='none';
	}
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
	if (myForm.ptmScope.value=='multiPTM') {
		document.getElementById('refQuantifDIV').style.display=(algoType.match('MSstats\|Abundance\|SSPA'))? 'none' : '';
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
		//if (myForm.refProt.selectedIndex==1) {
		//	myForm.refProt.selectedIndex=0;
		//}
		myForm.refProt.options[1].disabled=true;
		updateBiasCorrectionSet(0);
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
function ajaxFetchRefProtPeptides(selVal,call='protCorr') {
	var manualRefPepDiv=(call=='protCorr')? document.getElementById('manualRefPepDIV') : document.getElementById('manualBiasPepDIV');
	if (selVal=='0' \|\| (call=='biasCorr' && selVal > 0)) {
		manualRefPepDiv.style.display='none';
		manualRefPepDiv.innerHTML='';
		return;
	}
	// Check reference dataset
	if (call=='protCorr' && !document.selAnaForm.refQuantif.value) { // "-= Select =-" <- A design was selected but not a quantif ratio
		alert('Select a reference dataset first!');
		document.selAnaForm.manualRefPep.selectedIndex=0;
		return;
	}
	// Check states (at least 2 must be selected)
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
	//Check quantified modifs
	var quantifModifList=[];
	if (document.selAnaForm.ptmScope.value=='multiPTM') {
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
	var myForm=document.selAnaForm;
	if (document.getElementById('algoType').value.match('^Ratio\|MSstats')) {
		document.getElementById('biasProtDIV').style.visibility=(bias.match('REF_PROT\|globalStandards'))? 'visible' : 'hidden';
		//myForm.refProt.selectedIndex=0;
		updateBiasCorrectionSet(0);		
	}
	else {
		if (bias && !bias.match('^none\|quantile')) {
			document.getElementById('biasProtDIV').style.visibility='visible';
			// Manual pept. for bias correction (only possible for LF SimpleRatio in myProMS)
			if (myForm.ptmScope.value=='multiPTM' && myForm.algoType.value.match('(PEP_INTENSITY\|DIA\|TDA):SimpleRatio') && myForm.refQuantifType.value*1 == 0) {
				myForm.refProt.options[1].disabled=(myForm.refQuantifType.value*1==0 && $numProtInDesign > $MAX_PROT_MANUAL_PEP)? false : true; // Incompatible with protein-level FC correction
			}
			else {
				if (myForm.refProt.selectedIndex==1) {
					myForm.refProt.selectedIndex=0;
				}
				myForm.refProt.options[1].disabled=true;
			}				
		}
		else {
			document.getElementById('biasProtDIV').style.visibility='hidden';
			//myForm.refProt.selectedIndex=0;
			updateBiasCorrectionSet(0);
		}
	}
}
function updateBiasCorrectionSet(setCode) { // 0:all prot / -1:manual pep selection / <listID>:custom list
	// "(Do not) Use" switch
	var selType=document.selAnaForm.refProtSelType;
	if (setCode <= 0) {
		if (setCode == 0) {document.selAnaForm.refProt.selectedIndex=0;} // updateBiasCorrectionSet also called from other JS functions
		selType.selectedIndex==0;
		selType.options[1].disabled=true;
		//if (setCode == -1) { // Manual peptide selection
		//	ajaxFetchRefProtPeptides(setCode,'biasCorr');
		//}
	}
	else { // Custom list
		selType.options[1].disabled=false;
	}
	ajaxFetchRefProtPeptides(setCode,'biasCorr'); // Manual peptide selection
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
	if (propagScope=='all' && !force && !confirm('WARNING: This peptide quantification will be propagated to ALL other state(s) were available. Proceed?')) {return;}
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
	var myForm=document.selAnaForm;

	/* Update quantif name %...% span & States display */
	var nameSpanDisp='none',
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
	var dispRatioPepDiv='none',dispSspaPepDiv='none',disabAnyPep=true;
	if (algoType=='SSPA') {
		dispSspaPepDiv='';
// DISABLE SITES FOR NOW ------->
		if (myForm.ptmScope.selectedIndex != 0) { // Whole protein
			myForm.ptmScope.selectedIndex=0;
			updateQuantifFocus(myForm.ptmScope.value);
		}
		for (let i=1; i<myForm.ptmScope.options.length; i++) { // disable PTM quantif
			myForm.ptmScope.options[i].disabled=true;
		}
	}
	else {
		if (!algoType.match('select\|MSstats')) {
			dispRatioPepDiv='';
			if (algoType.match('Abundance')) {disabAnyPep=false;}
		}
		for (let i=1; i<myForm.ptmScope.options.length; i++) { // enable PTM quantif
			myForm.ptmScope.options[i].disabled=($numPTMs)? false : true;
		}
	}
	document.getElementById('ratioPepDIV').style.display=dispRatioPepDiv;
	document.getElementById('sspaPepSPAN').style.display=dispSspaPepDiv;
	document.getElementById('rescuedPepDIV').style.display=dispSspaPepDiv;
	document.getElementById('hiddenProtSPAN').style.display=dispSspaPepDiv;
	if (disabAnyPep) {myForm.matchingPep.selectedIndex=1;} // matching pep
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
		for (var i=0; i<algoNorms[trueAlgo].length; i++) { //}
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
	myForm.minInfRatios.disabled=(algoType.match('MSstats\|TDA\|SSPA\|Abundance'))? true : false;
	
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
	if (!multiPtmScope.checked) {return 0;}
	var ptmValue=multiPtmScope.value*1, ptmID=Math.abs(ptmValue);
	if (ptmValue > 0) { // normal cases (Matched)
		/* PhosphoProbSoft settings */
		if (('$ptmProbSoftCode'=='MQ' \|\| ptmID==$phosphoID) && (myForm.okProb.value < 0 \|\| myForm.okProb.value > 100)) {
			alert('ERROR: $ptmProbSoft probability must be between 0 and 100%');
			return -1;
		}
	}
	else {
		/* Recreated sites */
		var tgtResParam=myForm['targetRes_'+ptmID];
		var okRes=false;
		if (tgtResParam.length) {
			for (let j=0; j<tgtResParam.length; j++) {
				if (tgtResParam[j].checked && !tgtResParam[j].disabled) {
					okRes=true;
					break;
				}
			}
		}
		else if (!tgtResParam.disabled) {okRes=tgtResParam.checked;}
		if (!okRes) {
			alert('ERROR: Select at least 1 target residue for '+multiPtmScope.dataset.name);
			return -1;
		}
	}
	return 1;
}
function checkForm(myForm) {
	/* Quantif name */
	if (!myForm.quantifName.value) {alert('ERROR: Missing name for quantification!'); return false;}
	
	/* PTM sites */
	if (myForm.ptmScope.value=='multiPTM') {
		var numSelected=0;
		if (myForm.multiPtmScope.length) {
			for (let i=0; i<myForm.multiPtmScope.length; i++) {
				var ptmStatus=checkPtmSelection(myForm,myForm.multiPtmScope[i]);
				if (ptmStatus==1) {numSelected++;}
				else if (ptmStatus==-1) {return false;}
			}
		}
		else {
			var ptmStatus=checkPtmSelection(myForm,myForm.multiPtmScope);
			if (ptmStatus==1) {numSelected++;}
			else if (ptmStatus==-1) {return false;}
		}
		if (numSelected==0) {alert('ERROR: Select at least 1 PTM to be quantifed !'); return false;}
	}
			
	/* PTM filter */
	if (Math.abs(myForm.pepPTM.value)==2) {
		var okPTM=false;
		if (myForm.customPTM.length) {
			for (var i=0; i<myForm.customPTM.length; i++) {
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
	if ( myForm.algoType.value.match('select') ) {
		alert('ERROR: Select an algorithm!');
		return false;
	}

	/* States & Observations selection */
	var numCondition=0;
	var refStateAnaList={}; // for labeled quantif only
	for (var i=1; i <= $nbCond; i++) {
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
		var matchedAna=(isLabeled)? false : true;
		var condObsBox=myForm['condObs:'+condID];
		if (condObsBox.length) {
			for (var j=0; j < condObsBox.length; j++) {
				if (i > 1 && obsCondOK && matchedAna) break;
				if (condObsBox[j].checked && myForm['anaXICQuanti:'+condObsBox[j].value].value) {
					obsCondOK=true;
					if (isLabeled) {
						var obsData=myForm['anaXICQuanti:'+condObsBox[j].value].value.split(':');
						if (i==1) {refStateAnaList[obsData[2]]=1; matchedAna=true;} // records anaID
						else if (refStateAnaList[obsData[2]]) {matchedAna=true;} // at least 1 match is enough
					}
				}
			}
		}
		else if (condObsBox.checked && myForm['anaXICQuanti:'+condObsBox.value].value) {
			obsCondOK=true;
			if (isLabeled) {
				var obsData=myForm['anaXICQuanti:'+condObsBox.value].value.split(':');
				if (i==1) {refStateAnaList[obsData[2]]=1;matchedAna=true;} // records anaID
				else if (refStateAnaList[obsData[2]]) {matchedAna=true;} // at least 1 match is enough
			}
		}
		if (!obsCondOK) {
			alert('ERROR: No observation selected for State #'+i+'!');
			return false;
		}
		if (!matchedAna && myForm.algoType.value != 'SSPA') {
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
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" >
<CENTER>
<FONT class="title">$titleString</FONT>
$oldAlgoAccessStrg
<BR>
<FORM name="selAnaForm" method="post" onsubmit="return(checkForm(this));">
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
if ($userInfo[1]=~/mass|bioinfo/) {
	print qq
|<TR><TH class="title3" align=right valign=top>Apply template :</TH><TD bgcolor="$lightColor" class="title3" nowrap><SELECT id="template" name="template" class="title3" onchange="useTemplate(this.value)">
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
}
print qq
|</SELECT></TD></TR>
<TR><TH align=right class="title3" nowrap>&nbsp;Strategy :</TH><TD bgcolor="$lightColor" class="title3" nowrap>&nbsp;Labeling:<SELECT name="labeling" class="title3 template" onchange="clearTemplate();updateLabeling(this.value)"><OPTION value="">-= Select =-</OPTION>
|;
my $selLabelStrg=(scalar keys %{$designInfo{'LABEL'}} == 1)? 'selected' : '';
foreach my $labelType (sort keys %{$designInfo{'LABEL'}}) {
	print "<OPTION value=\"$labelType\" $selLabelStrg";
	#print ' selected' if scalar keys %{$designInfo{'LABEL'}}==1;
	if ($labelType eq 'FREE') {print ">Label-free</OPTION>";}
	elsif ($labelType eq 'ITRAQ') {print ">iTRAQ</OPTION>";}
	else {print ">$labelType</OPTION>";}
}
my ($disabSites,$siteMessage)=($numPTMs)? ('','PTM sites') : ('disabled','**No project PTMs**');
print qq
|</SELECT>
&nbsp;&nbsp;&nbsp;&nbsp;Method:<SELECT id="algoType" name="algoType" class="title3 template" onchange="updateAlgoTypeSelection(this.value)">
<OPTION value="select">-= Select =-</OPTION>
<OPTGROUP label="Peptide ratio:">
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
</SELECT>

<SPAN id="topSPAN" class="template" style="display:none">&nbsp;&nbsp;&nbsp;&nbsp;Use:<SELECT name="matchingPep" class="title3 template" disabled><OPTION value="0" disabled>any</OPTION><OPTION value="1" selected>matching</OPTION></SELECT>&nbsp;top<SELECT name="topN" class="title3 template" disabled><OPTION value="1">1</OPTION><OPTION value="2">2</OPTION><OPTION value="3" selected>3</OPTION><OPTION value="5">5</OPTION><OPTION value="10">10</OPTION><OPTION value="0">N</OPTION></SELECT>&nbsp;peptides</SPAN>
<SPAN id="sspaPepSPAN" style="display:none" class="template">&nbsp;&nbsp;&nbsp;&nbsp;Use:<SELECT name="pepFocus" class="title3 template"><OPTION value="xic">Peptide abundance</OPTION><OPTION value="sp_count" selected>Peptide/spectrum matches</OPTION><OPTION value="all_ion">All peptide ions</OPTION><OPTION value="all_pep">All peptides</OPTION><OPTION value="dist_ion">Distinct peptide ions</OPTION><OPTION value="dist_pep">Distinct peptides</OPTION><OPTION value="dist_seq">Distinct peptide sequences</OPTION></SELECT></SPAN>
</TD></TR>
<TR><TH class="title3" align=right valign=top>Focus :</TH><TH align=left bgcolor="$lightColor" nowrap><SELECT name="ptmScope" class="title3 template" onchange="updateQuantifFocus(this.value)">
<OPTION value="protein">Whole proteins</OPTION>
<OPTION value="multiPTM" $disabSites>$siteMessage</OPTION>
</SELECT>
<SPAN id="ptmSPAN" style="display:none" class="template">&nbsp;&nbsp;<INPUT type="checkbox" name="ambiguousPos" value="1" class="template">Positions of modification sites are not reliable</SPAN>
|;
if ($ptmProbSoftCode eq 'MQ') { # Position probabilities available for all modifs
	print qq
|<SPAN id="phosphoSPAN" style="display:none" class="template">&nbsp;&nbsp;Positions with $ptmProbSoft probability&ge;<INPUT type="text" name="okProb" class="title3 template" value="75" size="2">% are confirmed. The others are <SELECT name="badProb" class="title3 template"><OPTION value="exclude" selected>excluded</OPTION><OPTION value="ambiguous">delocalized</OPTION></SELECT></SPAN>
|;
}
#my $numPtmMatched=0;
#foreach my $modID (sort{lc($projectVarMods{$a}[0]) cmp lc($projectVarMods{$b}[0])} keys %projectVarMods) {
#	next if $unusedProjVarMods{$modID};
#	$numPtmMatched++;
#	print "<OPTGROUP label=\"Matched sites:\">\n" if $numPtmMatched==1;
#	print "<OPTION value=\"$modID\">$projectVarMods{$modID}[0]-sites</OPTION>\n";
#}
#print "</OPTGROUP>\n" if $numPtmMatched;
#if (scalar keys %unusedProjVarMods) {
#	print "<OPTGROUP label=\"Other sites:\">\n";
#	foreach my $modID (sort{lc($projectVarMods{$a}[0]) cmp lc($projectVarMods{$b}[0])} keys %unusedProjVarMods) {
#		print "<OPTION value=\"-$modID\">$projectVarMods{$modID}[0]-sites</OPTION>\n";
#	}
#	print "</OPTGROUP>\n";
#}
## UL with all PTMs
print "<TABLE id=\"multiPtmUL\" style=\"display:none\" cellspacing=0>\n"; # class=\"ulList\"
if ($numPTMs > scalar keys %unusedProjVarMods) {
	print "<TR><TD class=\"title3\">Matched sites:</TD></TR>\n" if scalar keys %unusedProjVarMods;
	foreach my $modID (sort{lc($projectVarMods{$a}[0]) cmp lc($projectVarMods{$b}[0])} keys %projectVarMods) {
		next if $unusedProjVarMods{$modID};
		print "<TR><TD>&nbsp;&nbsp;<INPUT type=\"checkbox\" name=\"multiPtmScope\" id=\"multiPtmScope_$modID\" value=\"$modID\" onchange=\"updateMultiPtmSelection()\"><LABEL for=\"multiPtmScope_$modID\">$projectVarMods{$modID}[0] [<FONT style=\"font-weight:bold;color:#$projectVarMods{$modID}[2]\">$projectVarMods{$modID}[1]</FONT>\]-sites</LABEL>";
		if ($modID==$phosphoID && $ptmProbSoftCode ne 'MQ') { # Position probabilities available only for Phospho (PRS)
			print qq
|<SPAN id="phosphoSPAN" style="display:none" class="template">: Positions with $ptmProbSoft probability&ge;<INPUT type="text" name="okProb" class="template" value="75" size="2">% are confirmed. The others are <SELECT name="badProb" class="template"><OPTION value="exclude" selected>excluded</OPTION><OPTION value="ambiguous">delocalized</OPTION></SELECT></SPAN>
|;
		}
		print "</TD></TR>\n";
	}
}
if (scalar keys %unusedProjVarMods) {
	print "<TR><TD class=\"title3\">Other sites:</TD></TR>\n";
	foreach my $modID (sort{lc($projectVarMods{$a}[0]) cmp lc($projectVarMods{$b}[0])} keys %unusedProjVarMods) {
		#print "<OPTION value=\"-$modID\">$projectVarMods{$modID}[0]-sites</OPTION>\n";
		print qq
|<TR><TD>&nbsp;&nbsp;<INPUT type="checkbox" name="multiPtmScope" id="multiPtmScope_$modID" value="-$modID" data-name="$projectVarMods{$modID}[0]-sites" onchange="updateMultiPtmSelection()"><LABEL for="multiPtmScope_$modID">$projectVarMods{$modID}[0] [<FONT style="font-weight:bold;color:#$projectVarMods{$modID}[2]">$projectVarMods{$modID}[1]</FONT>]-sites&nbsp;</LABEL>
<SPAN id="targetResSPAN_$modID" class="template">targetting free residues:
|;
		my $resIdx=-1;
		foreach my $resInfo (@{$unusedProjVarMods{$modID}}) {
			$resIdx++;
			print qq|<INPUT type="checkbox" name="targetRes_$modID" class="template" value="$resInfo->[0]" onchange="updateTargetResidues($modID,this)">$resInfo->[1]|;
			print ',' if $resIdx < $#{$unusedProjVarMods{$modID}};
		}
		print qq
|&nbsp;</SPAN>
</TD></TR>
|;
	}
}
print "</TABLE>\n";

#<SPAN id="recreateSPAN" style="display:none" class="template">&nbsp;&nbsp;Sites will be <FONT color=#DD0000>recreated</FONT> on selected residues:|;
#foreach my $modID (keys %unusedProjVarMods) {
#	print "<SPAN id=\"targetResSPAN_$modID\" style=\"display:none\" class=\"template\">&nbsp;";
#	my $resIdx=-1;
#	foreach my $resInfo (@{$unusedProjVarMods{$modID}}) {
#		$resIdx++;
#		print "<INPUT type=\"checkbox\" name=\"targetRes_$modID\" class=\"template\" value=\"$resInfo->[0]\">$resInfo->[1]";
#		print ',' if $resIdx < $#{$unusedProjVarMods{$modID}};
#	}
#	print "&nbsp;</SPAN>";
#}	
#print qq
#|</SPAN>
print qq
|<DIV id="refQuantifDIV" style="display:none">
<FONT class="title3">Protein-level fold change correction:</FONT><SELECT name="refQuantifType" class="title3" onchange="ajaxFetchReferenceQuantif(this.value)"><OPTION value="0">None</OPTION><OPTION value="-1">Use current dataset</OPTION><OPTGROUP label="Use dataset from:">
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
	<TD valign="top" nowrap>&nbsp;<B>Specificity<SUP onmouseover="popup('<B>For each protein:</B><BR>&nbsp;&bull;<B>Proteotypic:</B> Only proteotypic peptides are used.<BR>&nbsp;&bull;<B>Proteo.+shared:</B> all peptides are used if at least 1 proteotypic exists.<BR>&nbsp;&bull;<B>All:</B> all available peptides are used.')" onmouseout="popout()">?</SUP>:</B><SELECT name="pepSpecifity" class="template"><OPTION value="unique">Proteotypic</OPTION><OPTION value="unique_shared">Proteotypic + shared</OPTION><OPTION value="all" selected>All</OPTION></SELECT>
	&nbsp;&nbsp;&nbsp;<B>Missed cleav.:<B><SELECT name="pepCleavage" id="pepCleavage" class="template"><OPTION value="1">Allowed</OPTION><OPTION value="0" selected>Not allowed</OPTION><OPTION value="-1">Also exclude cleaved counterparts</OPTION></SELECT>
	&nbsp;&nbsp;&nbsp;<B>Modifications:<B>
	<DIV id="rescuedPepDIV" style="display:none"><LABEL><INPUT type="checkbox" name="pepRescued" value="1"/><B>Include data from MBR-rescued peptides</B></LABEL></DIV>
	</TD>
	<TD valign="top" nowrap><SELECT name="pepPTM" onchange="updateCustomPtmDisplay(this.value)" class="template"><OPTION value="1">Allow all</OPTION><OPTION value="2">Allow selected</OPTION>
	<OPTION value="0" selected>Exclude all</OPTION><OPTION value="-1">Exclude all & unmodified forms</OPTION>
	<OPTION value="-2">Exclude selected & unmodified forms</OPTION>
	</SELECT><BR>
	<DIV id="customPtmDIV" class="template" style="display:none">|;
foreach my $modID (sort{lc($modifications{$a}) cmp lc($modifications{$b})} keys %modifications) {
	$modifications{$modID}=&promsMod::resize($modifications{$modID},40);
	print "<INPUT type=\"checkbox\" class=\"template\" name=\"customPTM\" value=\"$modID\"/> $modifications{$modID}<BR>\n";
}
print "<FONT style=\"color:#DD0000\">There are no modifications involved!</FONT>\n" unless scalar keys %modifications;
print "<INPUT type=\"hidden\" name=\"listPTM\" value=\"",join(',',sort{$a<=>$b} keys %modifications),"\">\n";
my $disabManualStrg=($numProtInDesign > $MAX_PROT_MANUAL_PEP)? 'disabled' : '';
print qq |</DIV></TD>
	<TD valign="top" nowrap><DIV id="ratioPepDIV" class="template">&nbsp;&nbsp;&nbsp;<B>Charges:<B><SELECT name="pepCharge" class="template" onchange="updatePeptideSource(this.value)"><OPTION value="all">All</OPTION><OPTION value="best">Best signal</OPTION></SELECT>
&nbsp;&nbsp;&nbsp;<B><SUP>°</SUP>Sources:<B><SELECT name="pepSource" class="template" id="pepSource"><OPTION value="all">All</OPTION><OPTION value="best">Best signal</OPTION></SELECT>
</DIV>
</TD></TR></TABLE>
<SPAN id="manualRefPepSPAN" style="display:none">
&nbsp;<B>Protein-level fold change correction dataset:</B><SELECT name="manualRefPep" onchange="ajaxFetchRefProtPeptides(this.value,'protCorr')"><OPTION value="0">Current rules</OPTION><OPTION value="1" $disabManualStrg>Manual selection</OPTION></SELECT>
</SPAN>
<DIV id="manualRefPepDIV" style="background-color:white;padding:2px;max-height:500px;overflow:auto;display:none"></DIV>
</TD></TR>
<TR><TH align=right nowrap>Protein selection :</TH><TD bgcolor=$lightColor nowrap><SELECT name="protSelType" class="template"><OPTION value="exclude">Exclude</OPTION><OPTION value="restrict">Restrict to</OPTION></SELECT>&nbsp;<B>proteins from List:</B><SELECT name="protSelection" class="template"><OPTION value="">-= Select =-</OPTION>
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
<SPAN id="hiddenProtSPAN" style="display:none">&nbsp;&nbsp;&nbsp;<LABEL><INPUT type="checkbox" name="protHidden" value="1"/><B>Include proteins where hidden if also found visible</B></LABEL></SPAN>
&nbsp;</TD></TR>
<TR><TH align=right nowrap valign=top>&nbsp;Quantification settings :
	</TH><TD bgcolor=$lightColor nowrap valign=top>
	<TABLE border=0 cellpadding=0 cellspacing=0><TR>
	<TH align=right valign="top">&nbsp;-Bias correction:</TH>
	<TD><SELECT name="biasCorrect" class="template" id="biasCorrect" onchange="updateBias(this.value)"><OPTION value="">-= Select =-</OPTION><OPTION value="none.none">None</OPTION></SELECT>&nbsp;</TD>
	<TD valign="top" nowrap><DIV id="biasProtDIV" class="template" style="visibility:hidden">&nbsp;<B><SELECT name="refProtSelType" class="template"><OPTION value="use">Use</OPTION><OPTION value="not" disabled>Do not use</OPTION></SELECT>
data from</B> <SELECT name="refProt" class="template" onchange="updateBiasCorrectionSet(this.value*1)"><OPTION value="0">All proteins</OPTION><OPTION value="-1" disabled>Manual peptide selection</OPTION>
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
|	</SELECT> <B>for normalization.</B>&nbsp;</DIV></TD>
	</TR></TABLE>
<DIV id="manualBiasPepDIV" style="background-color:white;padding:2px;max-height:500px;overflow:auto;display:none"></DIV>
	&nbsp;<B>-Biological replicates:</B><SELECT name="residualVar" class="template" onchange="document.getElementById('resVarSPAN').style.visibility=(this.value=='biological')? 'visible' : 'hidden'"><OPTION value="auto" selected>Auto-detect*</OPTION><OPTION value="technical">No</OPTION><OPTION value="biological">Yes</OPTION></SELECT>
<SPAN id="resVarSPAN" style="visibility:hidden" class="template"><FONT style="font-weight:bold;color:#DD0000"> Requires at least 2 biological replicates in each state</FONT></SPAN><BR>
	&nbsp;<B>-p-value correction:</B><SELECT name="fdrControl" class="template"><OPTION value="none">None</OPTION><OPTION value="fdr" selected>Benjamini-Hochberg*</OPTION><OPTION value="bonferroni">Bonferroni</OPTION></SELECT><BR>
	&nbsp;<B>-</B><INPUT type="checkbox" name="minInfRatios" class="template" value="1"><B>Avoid infinite ratios when possible.</B><BR>
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
	my $ptmScope=param('ptmScope');
	my @ptmScopeModifIDs;
	if ($ptmScope eq 'multiPTM') { # only for DB: PTM quantification
		@ptmScopeModifIDs=(sort{abs($a)<=>abs($b)} param('multiPtmScope')); # sort is critical: Modif_rank is based on increasing modifID
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
	my $currentQuantifDir="$promsPath{tmp}/quantification/current";
	mkdir $currentQuantifDir unless -e $currentQuantifDir;

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

		my $quantifDir="$promsPath{tmp}/quantification/$jobDir";
		mkdir $quantifDir || die "ERROR detected: $!";

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
		open (INFO,">$quantifDir/quantif_info.txt"); # item valueR valueDB
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
		print INFO "TOP_N\t\t",param('topN'),"\n" if param('topN'); # only for DB not R (topN is undef if PTM quantif)
		print INFO "MATCH_PEP\t\t",param('matchingPep'),"\n" if param('matchingPep'); # only for DB not R
		#my $ptmScope=param('ptmScope');
		#my @ptmScopeModifIDs;
		if ($ptmScope eq 'multiPTM') { # only for DB: PTM quantification
			#@ptmScopeModifIDs=(sort{abs($a)<=>abs($b)} param('multiPtmScope')); # sort is critical: Modif_rank is based on increasing modifID
			print INFO "PTM_SCOPE\t\t",join("\t",@ptmScopeModifIDs),"\n";
			my ($createPTMStrg,$hasMatchedPTMs,$hasPhosphoID)=('',0,0);
			#my $phosphoID=param('phosphoID');
			my $modifRank=0;
			foreach my $modID (@ptmScopeModifIDs) {
				$modifRank++;
				if ($modID > 0) {
					$hasMatchedPTMs++;
					$hasPhosphoID=1 if $modID==$phosphoID; # phospho has "Other" is not considered
				}
				else {
					$modID*=-1;
					$createPTMStrg="CREATE_PTM\t" unless $createPTMStrg;
					$createPTMStrg.="\t$modifRank\:".join(".",param("targetRes_$modID"));
				}
			}
			if ($createPTMStrg) {
				print INFO "$createPTMStrg\n";
			}
			if ($hasMatchedPTMs) { # normal PTMs
				print INFO "PTM_POS\t\t";
				my $ptmProbSoftCode=param('ptmProbSoftCode');
				if (($hasPhosphoID || $ptmProbSoftCode eq 'MQ') && $algoType ne 'MSstats') { # phospho-quantif w PRS or anyPTM-quantif with MaxQuant (Not for DIA:MSstats)
					print INFO "$ptmProbSoftCode:",param('okProb'),"\t",param('badProb'),"\n"; # PRS or MQ
				}
				else { # non-phospho PTM quantif
					my $ambigPos=(param('ambiguousPos') && $hasMatchedPTMs==1 && !$createPTMStrg)? 'ambiguous' : 'valid';
					print INFO "$ambigPos\n";
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
				foreach my $ptmID (split(',',param('listPTM'))) {push @allowed,$ptmID if (!$notAllowed{$ptmID} && $ptmID != param('ptmScope'));} # Quantified PTM not listed here
				$ptmFilter='-2:#'.join(',#',@allowed);
			}
		}
		else {$ptmFilter=param('pepPTM');} # old options

		if ($algoType eq 'SSPA') {
			print INFO "PEPTIDES\t\t",param('pepSpecifity'),"\t",param('pepCleavage'),"\t$ptmFilter\t",param('pepFocus'),"\t";
			print INFO (param('pepRescued'))? "0\n" : "1\n"; # flag for skipping MBR data 0=>MBR, 1=>No MBR (back compatibility)
			print INFO "HIDDEN_PROT\t\t",(param('pepSpecifity') || 0),"\n";
		}
		else {
			print INFO "PEPTIDES\t\t",param('pepSpecifity'),"\t",param('pepCleavage'),"\t$ptmFilter\t$pepCharge\t$pepSrc","\n"; # only for DB ***Leave pepCharge & pepSrc for SWATH for compatibility with showProtQuantif***
			
			if ($ptmScope && param('refQuantifType') && param('refQuantifType')==-1) {
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
		print INFO "ALGO_TYPE\t\t$algoType\n"; # if $ratioType ne 'Ratio'; # needed for DB (not needed if Ratio) / Overwritten by launchQuantif if PEP_xxx incompatibilty
		if ($algoType ne 'SSPA') {
			my $ratioTypeCode=($ratioType eq 'Abundance')? 'None' : $ratioType;
			print INFO "RATIO_TYPE\t\t$ratioTypeCode\n";
			if ($ratioType=~/S\w+Ratio/) {
				print INFO "SINGLE_REF\t\t1\n" if ($refType && $refType eq 'singleRef');
			}
			elsif ($ratioType eq 'Abundance') { # <====== TEMP (parameter name to be changed)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				print INFO "\tnormalization.only\tyes\n";
				print INFO "\tresidual.variability\ttechnical\n"; # useless but still required
				print INFO "\tpAdj.method\tfdr\n"; # useless but still required
			}
			if ($ratioType eq 'Abundance') {
				my $minPepRatioLFQ=($ptmScope eq 'multiPTM')? 1 : 2; # Hard-coded but could become a user-defined parameter
				print INFO "NUM_LFQ_RATIOS\t\t$minPepRatioLFQ\n";
			}
			#<Bias correction
			my $normMethod=param('biasCorrect');
			my $biasCorrect=($normMethod eq 'none.none')? 'FALSE' : 'TRUE';
			print INFO "BIAS_CORRECTION\t\t$biasCorrect"; # for DB only
			print INFO "\t#",param('refProt'),"\t",param('refProtSelType') if (param('refProt') > 0 && ($normMethod=~/REF_PROT|globalStandards/ || ($algoType=~/^(PEP_|TDA|DIA)/ && $normMethod !~/quantile/i))); # all algos except SSPA
			print INFO "\n";
			if ($ratioType=~/S\w+Ratio|Abundance/) { # S*Ratio & MSstats
				my ($usedMethod)=($normMethod eq 'none.none' && $algoType eq 'MSstats')? 'FALSE' : $normMethod;
				$usedMethod=~s/scale/none/ if ($usedMethod=~/scale/ && param('refProt')==-1 && scalar param('manualPep') == 1); # no scale if only 1 peptide ion for normalization
				print INFO "NORMALIZATION_METHOD\tnormalization.method\t$usedMethod\n";

				if ($ratioType=~/S\w+Ratio/) { # since SuperRatio
					print INFO "PEPTIDES_NORM\t\t",join("\t",('manual',param('manualPep'))),"\n" if (param('refProt')==-1 && param('manualPep')); # manual\tprotID:pep1ID,pep2ID,...,pepNID
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
						my $clusters=($clusterInfo{'on'} && $clusterInfo{'maxCPUs'})? $clusterInfo{'maxCPUs'} : 1;
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
		print INFO "MIN_INF_RATIOS\t\t$minInfRatios\n" if ($algoType !~ /MSstats|SSPA/ && $ratioType ne 'Abundance'); # if $usedStateNum==2; # +/- infinite ratio switch (DB only)

		open(FLAG,">$currentQuantifDir/$jobDir\_request.flag"); # flag file used by watchQuantif
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
	print "<BR><FONT class=\"title3\"><BR>This page will refresh itself in a few seconds.</FONT>\n";

	###>Calling watch popup window<###
	sleep 3;
	#my $fullQuantifType = ($algoType=~/^(MSstats|SSPA|TDA)$/)? 'DESIGN:'.$algoType : 'DESIGN';
	print qq
|<SCRIPT type="text/javascript">
var monitorJobsWin=window.open("$promsPath{cgi}/monitorJobs.cgi?filterType=Quantification [$fullQuantifType]&filterDateNumber=1&filterDateType=DAY&filterStatus=Queued&filterStatus=Running",'monitorJobsWindow','width=1200,height=500,scrollbars=yes,resizable=yes');
monitorJobsWin.focus();
</SCRIPT>
|;
#exit; # DEBUG!!!!!
	sleep 5;
	print qq
|<SCRIPT LANGUAGE="JavaScript">
//top.promsFrame.selectedAction='summary';
parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=design:$designID&ACT=experiment&EXPERIMENT=EXPERIMENT:$experimentID&ISNAVFRAME=0&VIEW=quanti";
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
	
	my $sthTestProt=$dbh->prepare("SELECT DISTINCT AP.ID_PROTEIN,ALIAS,PROT_DES,PROT_LENGTH,ORGANISM FROM EXPCONDITION E
										INNER JOIN OBS_EXPCONDITION OE ON E.ID_EXPCONDITION=OE.ID_EXPCONDITION
										INNER JOIN OBSERVATION O ON OE.ID_OBSERVATION=O.ID_OBSERVATION
										INNER JOIN ANALYSIS_PROTEIN AP ON O.ID_ANALYSIS=AP.ID_ANALYSIS
										INNER JOIN PROTEIN P ON AP.ID_PROTEIN=P.ID_PROTEIN
										WHERE ID_DESIGN=? AND ABS(AP.VISIBILITY)=2 GROUP BY AP.ID_PROTEIN LIMIT $MAX_PROT_MANUAL_PEP"); # LIMIT to be safe
	
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
	$sthTestProt->execute($designID);
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
<TR><TH width=20>&nbsp;&nbsp;</TH><TH bgcolor="$darkColor" class="rbBorder"" width=25>#</TH>
<TH bgcolor="$darkColor" class="rbBorder" width=50>&nbsp;Start&nbsp;</TH>
<TH bgcolor="$darkColor" class="rbBorder">Peptide&nbsp;&nbsp;</TH>
<TH bgcolor="$darkColor" class="rbBorder" width=50>&nbsp;Charge&nbsp;</TH>
<TH bgcolor="$darkColor" class="rbBorder" width=100>&nbsp;Proteotypic&nbsp;</TH>
<TH bgcolor="$darkColor" class="bBorder" width=60>&nbsp;Miss.&nbsp;cut&nbsp;</TH>
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
	<TH bgcolor="$bgColor" class="font11" align=left nowrap><INPUT type="checkbox" name="manualPep" value="$pepIdStrg"$disabPepSeq $popupStrg >$pepSeq$varMod{$modCode}&nbsp;</TH>
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


# TODO: handle PTM site re-creation in &ajaxFetchRefProtPeptides
####> Revision history
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

