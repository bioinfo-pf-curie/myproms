#!/usr/local/bin/perl -w

################################################################################
# selAna4Quantification.cgi      2.0.0                                         #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
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
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use promsQuantif;

#print header(-charset=>'UTF-8'); warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my %msTypeName=&promsConfig::getMsType;
my $userID=$ENV{'REMOTE_USER'};
$msTypeName{'MIS'}="MS/MS"; #redef, to keep space
my %quantifProcesses=('XICMCQ'=>'Ext. ion chrom.','EMPAI'=>'emPAI','SIN'=>'SI<SUB>N</SUB>','SILAC'=>'SILAC','ITRAQ'=>'iTRAQ','TMT'=>'TMT','DESIGN'=>'Protein-Ratio'); #,'MCQ'=>'Ext. ion chrom.''TnPQ or Peptide ratio'
my %xicSoftware=('PD'=>'Proteome Discoverer','MCQ'=>'MassChroQ','MAS'=>'Mascot','PAR'=>'Paragon','PKV'=>'PeakView','MQ'=>'MaxQuant','OS'=>'OpenSwath','?'=>'Unknown');
#my $MAX_CHANNELS=3;  # 3 max for SILAC
my $updateFrameString="";
my $maxLabPerChannel=10; # max num for isotope extraction...
my %normalizationNames=&promsQuantif::getQuantifNormalizationName;

###############################
####>Recovering parameters<####
###############################
my $quantifType=uc(param('quantifType')); # can be XICMCQ,EMPAI,SIN,DESIGN
my $branchID=param('ID');
my ($item,$itemID)=split(':',$branchID);
$item=lc($item);
my $ITEM=uc($item);

############################
####>form was submitted<####
############################
if (param('launch')) {
	&launchQuantifications;
	exit;
}

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;
my $projectID=&promsMod::getProjectID($dbh,$itemID,$item);
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
#my $projectAccess=${$userInfo[2]}{$projectID};
my $bioProjectAccessStrg=''; # for template selection
if ($quantifType eq 'DESIGN' && $userInfo[1]=~/^(bio|manag)$/) { # for template selection
	@userInfo=&promsMod::getUserInfo($dbh,$userID); # scan all accessible projects
	$bioProjectAccessStrg=join(',',keys %{$userInfo[2]});
}
my ($itemName)=$dbh->selectrow_array("SELECT NAME FROM $ITEM WHERE ID_$ITEM=$itemID");
my $selectionStrg=($quantifType eq 'DESIGN')? 'Observations' : 'Analyses';
my $titleString="Select $selectionStrg in ".&promsMod::getItemType($ITEM)." <FONT color=#DD0000>$itemName</FONT><BR>for $quantifProcesses{$quantifType} Quantification";
my $oldAlgoAccessStrg=($ENV{HTTP_HOST}=~/curie\.fr/ && $ITEM eq 'DESIGN')? "<BR>\n<INPUT type=\"button\" class=\"title3\" value=\"Use old alogrithms\" onclick=\"window.location='./selAna4Quantification_OA.cgi?ID=$branchID&quantifType=$quantifType'\"/>" : '';

####>Project PTMs<####
my $sthGetPM=$dbh->prepare("SELECT PM.ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES FROM PROJECT_MODIFICATION PM,MODIFICATION M WHERE PM.ID_MODIFICATION=M.ID_MODIFICATION AND PM.ID_PROJECT=$projectID");
$sthGetPM->execute;
my %projectVarMods;
while (my ($modID,$psiName,$interName,$synName)=$sthGetPM->fetchrow_array) {
	$projectVarMods{$modID}=$psiName || $interName;
	unless ($projectVarMods{$modID}) {
		$synName=~s/^##//; $synName=~s/##$//; $synName=~s/##/, /;
		$projectVarMods{$modID}=$synName;
	}
}
$sthGetPM->finish;

my ($phosphoID)=$dbh->selectrow_array('SELECT ID_MODIFICATION FROM MODIFICATION WHERE UNIMOD_ACC=21');
$phosphoID=-1 unless $projectVarMods{$phosphoID}; # disable if not used in project


########################################################
####>Recovering list of experiment, sample analysis<####
########################################################
my (%listDataBank,@itemAnalyses,%anaProteins,%listParam,%anaLabelMods,%modifications,%anaLabeling); #,%refQuantifications %okRemFilter,%okActLowScores,%listQuantif,
if ($quantifType ne 'DESIGN') {

	####>Recovering DBs name<####
	my $sthAD = $dbh->prepare("SELECT D.ID_DATABANK,NAME FROM ANALYSIS_DATABANK AD,DATABANK D WHERE AD.ID_DATABANK=D.ID_DATABANK AND AD.ID_ANALYSIS=?");

	####>Recovering quantitation name<####
	#my $sthQuanti = $dbh->prepare("SELECT ID_QUANTIF_METHOD,NAME FROM QUANTIF_METHOD");
	#$sthQuanti->execute;
	#while (my ($idQuanti,$quantiName)= $sthQuanti->fetchrow_array) {
	#	$listQuantif{$idQuanti}=$quantiName;
	#}
	#$sthQuanti->finish;

	my @sthItem;
	my $baseFieldString='ID_ANALYSIS,VALID_STATUS,MS_TYPE,DATA_FILE,FILE_FORMAT,WIFF_FILE,TAXONOMY,MAX_RANK,MIN_SCORE,INSTRUMENT,0,LABELING';
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
	}
	my $sthVVP=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE VISIBILITY>0 AND ID_ANALYSIS=?");
	my $sthAVP=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=?");
	my $sthSTP=$dbh->prepare("SELECT COUNT(*) FROM PROTEIN_VALIDATION WHERE SEL_STATUS>=1 AND ID_ANALYSIS=?");
	my $sthATP=$dbh->prepare("SELECT COUNT(*) FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=?");
	my $sthAM=$dbh->prepare("SELECT AM.ID_MODIFICATION,AM.SPECIFICITY,IS_LABEL,M.VALID_STATUS,M.MONO_MASS,M.PSI_MS_NAME,M.INTERIM_NAME,M.DES,M.SYNONYMES FROM ANALYSIS_MODIFICATION AM,MODIFICATION M WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND ID_ANALYSIS=? AND AM.MODIF_TYPE='V'");

	foreach my $sth (@sthItem) {
		$sth->execute;
		while (my @anaData=$sth->fetchrow_array) {
			my $anaID=$anaData[0];
			#>Databank(s)
			$sthAD->execute($anaID);
			my @dbUsed;
			while (my ($dbID,$dbName)=$sthAD->fetchrow_array) {
				push @dbUsed,$dbID;
				$listDataBank{$dbID}=$dbName;
			}
			$anaData[10]=\@dbUsed;
			push @itemAnalyses,\@anaData;
			if ($anaData[1]>=1) { # valid proteins
				$sthVVP->execute($anaID);
				push @{$anaProteins{$anaID}},$sthVVP->fetchrow_array;
				$sthAVP->execute($anaID);
				push @{$anaProteins{$anaID}},$sthAVP->fetchrow_array;
			}
			else {# non-validated proteins
				$sthSTP->execute($anaID);
				push @{$anaProteins{$anaID}},$sthSTP->fetchrow_array;
				$sthATP->execute($anaID);
				push @{$anaProteins{$anaID}},$sthATP->fetchrow_array;
			}
			$sthAM->execute($anaID);
			while (my ($modID,$anaSpec,$isLabel,$isValid,$monoMass,$psiName,$interName,$des,$synName)=$sthAM->fetchrow_array) {
				if ($isLabel) {
					@{$anaLabelMods{$modID}}=($anaSpec,$psiName,$interName,$isValid,$monoMass);
					$anaLabeling{$anaID}=(($psiName && $psiName=~/Label:/) || ($des && $des=~/Label:|SILAC/i) || ($synName && $synName=~/Label:|SILAC/i))? 'SILAC' : 'OTHER';
				}
				else {
					$modifications{$modID}=$psiName || $interName || $des;
					unless ($modifications{$modID}) {
						$synName=~s/^##//; $synName=~s/##$//; $synName=~s/##/, /;
						$modifications{$modID}=$synName;
					}
					$anaLabeling{$anaID}='FREE';
				}
			}
		}
		$sth->finish;
	}
	$sthAD->finish;
	foreach my $sth (@sthItem) {$sth->finish;}
	$sthVVP->finish;
	$sthAVP->finish;
	$sthSTP->finish;
	$sthATP->finish;
	$sthAM->finish;
}

####>Labeled quantifs<####
#my (%anaPeptideQuantifs,%quantifChannels,$numChannels,$maxReplicates); # SILAC, iTRAQ internal quanti
my (%categoryList,%quantifMethods,%quantifTemplate,%projectTemplate,%defaultAnaFormSelValues);
foreach my $quantMethod ('XIC','SILAC','ITRAQ','TMT','MQ','DIA') {$quantifMethods{$quantMethod}=0;}
my (%designInfo,%anaObs); # For TnPQ and Prop Pep Ratio (quantifications associated to a DESIGN)
my $nbCond=0;
my ($selCondStrg1,$selCondStrg2)=("<OPTION value=\"Select\">-= Select =-</OPTION>","<OPTION value=\"Select\">-= Select =-</OPTION>");
my $varGlobSelCond="var stateCondVal=[];\n"; # JS array to record Index of selected cond for each State
my $jsXicSoftStrg; # for DESIGN only
if ($quantifType eq 'DESIGN') {
	#($numChannels,$maxReplicates,$nbCond)=(0,0,0);

	###> Get OBSERVATION related to this design
	#my $sthAE=$dbh->prepare("SELECT NAME,EC.ID_EXPCONDITION,ID_ANALYSIS,TARGET_POS,FRACTION_GROUP,ID_OBSERVATION FROM EXPCONDITION EC,OBSERVATION O WHERE EC.ID_EXPCONDITION=O.ID_EXPCONDITION AND EC.ID_DESIGN=$itemID ORDER BY NAME ASC");
	#my $sthAE=$dbh->prepare("SELECT NAME,E.ID_EXPCONDITION,ID_ANALYSIS,TARGET_POS,FRACTION_GROUP,O.ID_OBSERVATION,GROUP_CONCAT(ID_MODIFICATION ORDER BY ID_MODIFICATION SEPARATOR ':') FROM OBSERVATION O LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION JOIN EXPCONDITION E ON O.ID_EXPCONDITION=E.ID_EXPCONDITION WHERE ID_DESIGN=$itemID GROUP BY O.ID_OBSERVATION ORDER BY ID_EXPCONDITION");
	my $sthAE=$dbh->prepare("SELECT E.NAME,E.ID_EXPCONDITION,BS.NAME,O.ID_ANALYSIS,O.TARGET_POS,OE.FRACTION_GROUP,OE.TECH_REP_GROUP,O.ID_OBSERVATION,GROUP_CONCAT(M.ID_MODIFICATION ORDER BY M.ID_MODIFICATION SEPARATOR ':')
								FROM OBSERVATION O
								LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION
								LEFT JOIN BIOSAMPLE BS ON O.ID_BIOSAMPLE=BS.ID_BIOSAMPLE
								JOIN OBS_EXPCONDITION OE ON O.ID_OBSERVATION=OE.ID_OBSERVATION
								JOIN EXPCONDITION E ON OE.ID_EXPCONDITION=E.ID_EXPCONDITION
								WHERE E.ID_DESIGN=$itemID GROUP BY O.ID_OBSERVATION ORDER BY E.DISPLAY_POS,E.ID_EXPCONDITION");
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
	my $sthQInfo=$dbh->prepare("SELECT Q.ID_QUANTIFICATION,Q.NAME,QUANTIF_ANNOT,STATUS,ID_ANALYSIS,UPDATE_DATE,QM.CODE FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ,QUANTIFICATION_METHOD QM WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND FOCUS='peptide' AND CODE != 'SIN' AND ID_ANALYSIS IN $anaString");
	my $sthIsobar=$dbh->prepare("SELECT DISTINCT M.ID_MODIFICATION FROM MODIFICATION M,ANALYSIS_MODIFICATION AM WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND IS_LABEL=1 AND ID_ANALYSIS=? LIMIT 0,1"); #  AND AM.MODIF_TYPE='V'
	my $sthAna=$dbh->prepare("SELECT FILE_FORMAT FROM ANALYSIS WHERE ID_ANALYSIS=?");
	$sthQInfo->execute;
	while (my ($quantiID,$quantiName,$quantifAnnot,$qStatus,$anaID,$date,$methodCode) = $sthQInfo->fetchrow_array) {
#$quantiName=~s/3plex.+\| /2plex (/;
		next unless $quantifAnnot;
#print "$labelType,/$quantiID,/$quantiName,/$qStatus,/$anaID,/$methodCode<BR>\n";
		$quantifMethods{$methodCode}=1;
		my $xicSoftCode;
		$quantifAnnot=~s/::SOFTWARE=(\w+);?(\d|\.)*//; # remove software info for back compatibility
		if ($1) {$xicSoftCode=$1;}
		elsif ($quantifAnnot=~/EXTRACTION_ALGO=/) {$xicSoftCode='MCQ';}
		else {
			$sthAna->execute($anaID);
			my ($fileFormat)=$sthAna->fetchrow_array;
			$xicSoftCode=($fileFormat=~/\.PDM\Z/)? 'PD' : ($fileFormat=~/^MASCOT/)? 'MAS' : ($fileFormat=~/PARAGON/)? 'PAR' : '?';
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
		$designInfo{'QUANTIFICATION'}{$quantiID}{'NAME'}=$quantiName;
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
		my $popupInfo='<U><B>Quantification parameters:</B></U>';
		if ($methodCode eq 'XIC') {
			$popupInfo.="<BR><B>Software:</B> $xicSoftware{$xicSoftCode}";
			if ($xicSoftCode eq 'PD') {
				foreach my $paramName (sort{lc($a) cmp lc($b)} keys %quantifParameters) {
					next if $paramName eq 'LABEL';
					$popupInfo.="<BR><B>$paramName:</B> $quantifParameters{$paramName}";
				}
			}
			elsif($xicSoftCode eq 'MQ') { # MaxQuant
				# Nothing for now
			}
			else { # assume MCQ
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

	###>JS list of quantif=>XIC software
	$jsXicSoftStrg="var quantiXicSoftware={\n";
	my $firstSoft=1;
	foreach my $quantiID (sort{$a<=>$b} keys %{$designInfo{'QUANTIFICATION'}}) {
		$jsXicSoftStrg.=",\n" unless $firstSoft;
		$jsXicSoftStrg.="\t$quantiID:'$designInfo{QUANTIFICATION}{$quantiID}{XIC_SOFTWARE}'";
		$firstSoft=0;
	}
	$jsXicSoftStrg.="\n};\n";


	###>List of potential reference quantifs in case mod quantif
	#my $sthRefQ=$dbh->prepare("SELECT ID_QUANTIFICATION,NAME FROM QUANTIFICATION WHERE ID_DESIGN=$itemID AND ID_MODIFICATION IS NULL");
	#$sthRefQ->execute;
	#while (my($qID,$qName)=$sthRefQ->fetchrow_array) {
	#	$refQuantifications{$qID}=$qName;
	#}
	#$sthRefQ->finish;

	###>List of non-label modification used
	my $sthMod=$dbh->prepare("SELECT DISTINCT M.ID_MODIFICATION,M.PSI_MS_NAME,M.INTERIM_NAME,M.DES,M.SYNONYMES FROM MODIFICATION M
								INNER JOIN ANALYSIS_MODIFICATION AM ON M.ID_MODIFICATION=AM.ID_MODIFICATION
								INNER JOIN OBSERVATION O ON AM.ID_ANALYSIS=O.ID_ANALYSIS
								INNER JOIN OBS_EXPCONDITION OC ON O.ID_OBSERVATION=OC.ID_OBSERVATION
								INNER JOIN EXPCONDITION C ON OC.ID_EXPCONDITION=C.ID_EXPCONDITION
								WHERE ID_DESIGN=$itemID AND IS_LABEL != 1 AND AM.MODIF_TYPE='V'"
							);
	$sthMod->execute;
	while (my ($modID,$psiName,$interName,$des,$synName)=$sthMod->fetchrow_array) {
		$modifications{$modID}=$psiName || $interName || $des;
		unless ($modifications{$modID}) {
			$synName=~s/^##//; $synName=~s/##$//; $synName=~s/##/, /;
			$modifications{$modID}=$synName;
		}
	}
	$sthMod->finish;

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

		next unless (($quantifAnnot!~/SOFTWARE/ && $quantifAnnot!~/ALGO_TYPE/ ) || ($quantifAnnot=~/SOFTWARE=DIA/ && $quantifMethods{'DIA'}) || ($quantifAnnot!~/SOFTWARE=DIA/ && ($quantifMethods{'XIC'} || $quantifMethods{'SILAC'} || $quantifMethods{'ITRAQ'} || $quantifMethods{'TMT'})));
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
}

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
<TITLE>Select Analysis</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TD.center {text-align:center}
</STYLE>
<SCRIPT type="text/javascript">
|;
&promsMod::popupInfo();
if ($quantifType eq 'XICMCQ') {
	###> For label purporse MCQXIC extraction
	print qq
|//For label purporse MCQXIC extraction
var refValue=-1;
var analabels=[];
analabels[0]=[-1,'Light',0.0000,''];
|;
	#print "var analabels=new Array();\n";
	#print "analabels[0]=new Array();\n";
	#print "analabels[0][0]=-1;\n";
	#print "analabels[0][1]=\"No Label\";\n";
	#print "analabels[0][2]=0.0000;\n";
	#print "analabels[0][3]=\"\";\n";
	my $labelIndx=0;
	foreach my $modID (keys %anaLabelMods) {
		$labelIndx++;
		my ($anaSpec,$psiMsName,$interimName,$isValid,$monoMass)=@{$anaLabelMods{$modID}};
		print "analabels[$labelIndx]=[$modID,'$interimName',$monoMass,'$anaSpec'];\n";
		#print "analabels[$labelIndx]=new Array();\n";
		#print "analabels[$labelIndx][0]=$modID;\n";
		#print "analabels[$labelIndx][1]=\"$interimName\";\n";
		#print "analabels[$labelIndx][2]=$monoMass;\n";
		#print "analabels[$labelIndx][3]=\"$anaSpec\";\n";
	}
	print qq
|function updateLabels(modID,channelNum,labelNum){
	for(var i=0; i< analabels.length; i++) {
		if(analabels[i][0] == modID) {
			var name=document.getElementById('lab'+labelNum+':name_channel'+channelNum);
			name.value=analabels[i][1];
			// var dm=document.getElementById('lab'+labelNum+':dm_channel'+channelNum);
			// dm.value=analabels[i][2];
			var sp=document.getElementById('lab'+labelNum+':sp_channel'+channelNum);
			sp.value=analabels[i][3];
		}
	}
}
function updateChargeState(xicType) {
	var checkBoxField=document.getElementById('allChargeStates');
	if (xicType == 'max'){checkBoxField.disabled=false;}
	else{checkBoxField.disabled=true;}
}
function updateTextState(name,target) {
	var textField=document.getElementById(name);
	if (target == 'ANY'){textField.disabled=false;}
	else{textField.disabled=true;}
}
function addQuantificationLabel(channelNum,labelNum,action) {
	if (action == 'show') {
		if(labelNum < $maxLabPerChannel) {
			var newLabel=labelNum+1;
			// Show-Hide button
			document.getElementById('show:'+channelNum+':'+labelNum).style.display='none';
			document.getElementById('hide:'+channelNum+':'+labelNum).style.display='none';
			document.getElementById('hide:'+channelNum+':'+newLabel).style.display='';
			if (newLabel == $maxLabPerChannel) {
				document.getElementById('show:'+channelNum+':'+newLabel).style.display='none';
			}
			// Show new label fieldset
			document.getElementById('fs:'+channelNum+':'+newLabel).style.display='';
		}
	}
	else {
		var newLabel=labelNum-1;
		document.getElementById('fs:'+channelNum+':'+labelNum).style.display='none';
		document.getElementById('show:'+channelNum+':'+newLabel).style.display='';
		if (newLabel > 1) {document.getElementById('hide:'+channelNum+':'+newLabel).style.display='';}
	}
}
|;
}
elsif ($quantifType eq 'DESIGN') {
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

	print qq
|window.onload=function() { // needed for design quantif because of auto-selection if only 1 labeling available
	updateLabeling(document.selAnaForm.labeling.value);
	updateAlgoTypeSelection(document.selAnaForm.algoType.value);
	recoverDefaultFormParam();
}

var defaultAnaFormValues=[];
function recoverDefaultFormParam(){
	var params=document.getElementsByClassName('template');
	var i=0;
	for (i=0; i<params.length ; i++){
		if (params[i].type=='select-one'){
			defaultAnaFormValues[i]=[params[i].name,params[i].id,params[i].selectedIndex,params[i].type];
		}
		else if (params[i].type=='checkbox'){
			if (document.getElementsByName(params[i].name).length>1){
				defaultAnaFormValues[i]=[params[i].name,params[i].id,params[i].checked,params[i].type,'',params[i].value];
			}
			else{
				defaultAnaFormValues[i]=[params[i].name,params[i].id,params[i].checked,params[i].type];
			}
		}
		else if (params[i].tagName == 'DIV' \|\| params[i].tagName == 'SPAN'){
			if (params[i].style.visibility){
				defaultAnaFormValues[i]=[params[i].tagName,params[i].id,params[i].style.visibility,'style.visibility'];
			}
			else{
				defaultAnaFormValues[i]=[params[i].tagName,params[i].id,params[i].style.display,'style.display'];
			}
		}
		else{
			defaultAnaFormValues[i]=[params[i].name,params[i].id,params[i].value,params[i].type];
		}
		if (params[i].disabled){
			defaultAnaFormValues[i][4]=params[i].disabled;
		}
	}
	/*var allDIV=document.getElementsByTagName("div");
	for (let j=0; j<allDIV.length ; j++){
		if (allDIV[j].id){
			if (allDIV[j].style.visibility){
				defaultAnaFormValues[i]=['DIV',allDIV[j].id,allDIV[j].style.visibility,'style.visibility'];
			}
			else{
				defaultAnaFormValues[i]=['DIV',allDIV[j].id,allDIV[j].style.display,'style.display'];
			}
			i++;
		}
	}
	var allSPAN=document.getElementsByTagName("span");
	for (let j=0; j<allSPAN.length ; j++){
		if (allSPAN[j].id){
			if (allSPAN[j].style.visibility){
				defaultAnaFormValues[i]=['SPAN',allSPAN[j].id,allSPAN[j].style.visibility,'style.visibility'];
			}
			else{
				defaultAnaFormValues[i]=['SPAN',allSPAN[j].id,allSPAN[j].style.display,'style.display'];
			}
			i++;
		}
	}*/
}


function useTemplate(quantiID,skipLabeling){
	// restore default parameters
	for (let i=0; i<defaultAnaFormValues.length; i++){
		if (defaultAnaFormValues[i][3] == 'hidden' \|\| defaultAnaFormValues[i][3] == 'button') continue;
		var paramIdent='';
		if (defaultAnaFormValues[i][1]!=''){
			paramIdent=document.getElementById(defaultAnaFormValues[i][1]);
		}
		else{
			paramIdent=document.getElementsByName(defaultAnaFormValues[i][0])[0];
		}

		if (defaultAnaFormValues[i][3] == 'select-one'){
			if (paramIdent.name != 'labeling' \|\| !skipLabeling) {
				paramIdent.selectedIndex=defaultAnaFormValues[i][2];
			}
		}
		else if (defaultAnaFormValues[i][3] == 'checkbox'){
			if (defaultAnaFormValues[i][5]){
				var selCheckbox=document.getElementsByName(defaultAnaFormValues[i][0]);
				for (let m=0; m<selCheckbox.length; m++){
					if (selCheckbox[m].value == defaultAnaFormValues[i][5]){
						selCheckbox[m].checked=defaultAnaFormValues[i][2]
					}
				}
			}
			else{
				paramIdent.checked=defaultAnaFormValues[i][2];
			}
		}
		else if (defaultAnaFormValues[i][3] == 'style.visibility'){
			paramIdent.style.visibility=defaultAnaFormValues[i][2];
		}
		else if (defaultAnaFormValues[i][3] == 'style.display'){
			paramIdent.style.display=defaultAnaFormValues[i][2];
			if (defaultAnaFormValues[i][1] == 'refCommentSpan'){
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


	// display template parameters
	if (quantiID == 'select') return;
	var quantifParam=quantifTemplate[quantiID][1].split('::');
	var software='';
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
	var prsPos;

	// Template : select template name
	var selTemplate=document.getElementById('template');
	for(var i=1; i<selTemplate.options.length;i++){
		if (selTemplate.options[i].value == quantiID){
			selTemplate.selectedIndex=i;
		}
	}

	for (var i=0; i<quantifParam.length; i++){
		var param=quantifParam[i].split('=');

		// Strategy
		if (param[0] ==  'SOFTWARE'){
			software=param[1].split(';');
		}
		if (param[0] ==  'SINGLE_REF'){
			singleRef=param[1];
		}
		if (param[0] ==  'RATIO_TYPE'){
			ratio=param[1];
		}
		if(param[0] ==  'ALGO_TYPE'){
			algoType=param[1];
		}
		if (param[0] == 'LABEL'){
			var selLabel=document.getElementsByName('labeling');
			for (let i=0; i<selLabel[0].options.length; i++) {
				if (selLabel[0].options[i].value==param[1]) {
					selLabel[0].selectedIndex=i;
					break;
				}
			}
			updateLabeling(param[1]);
		}
		if (param[0] == 'TOP_N'){
			topN=param[1];
		}

		// peptide selection
		if (param[0] == 'PEPTIDES'){
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

			if(selPTM[0]==2){
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
			else if (selPTM[0]==-2){
				var selCustomPTM=document.getElementsByName('customPTM');
				var mod=selPTM[1].split(/[,#]+/);
				var arrayCheckBox=[];
				for (let i=0; i<selCustomPTM.length; i++){
					arrayCheckBox[i]=selCustomPTM[i].value;
				}
				for (let j=1; j<mod.length ; j++){
					for (let i=0; i<arrayCheckBox.length ; i++){
						if (arrayCheckBox[i] == mod[j]) {
							arrayCheckBox.splice(i,1);
						}
					}
				}
				for (let i=0; i<selCustomPTM.length; i++) {
					for (let j=0; j<arrayCheckBox.length ; j++){
						if (selCustomPTM[i].value==arrayCheckBox[j]) {
							selCustomPTM[i].checked=true;
						}
					}
				}
			}
			if (peptideInfo[3]){
				pepCharge=peptideInfo[3];
			}
			if (peptideInfo[4]){
				pepSource=peptideInfo[4];
			}
		}

		// Target
		if (param[0] == 'PTM_POS'){
			if (param[1].match(';')){
				var ptmPosInfo=param[1].split(';');
				prsPos=ptmPosInfo[1];
			}
			else if (param[1] == 'ambiguous'){
				ambiguousPos=1;
			}
		}


		// protein selection
		if (param[0] == 'PROTEINS'){
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
		if (param[0] == 'BIAS_CORRECTION'){
			if (param[1].match(';')){
				var biasInfo=param[1].split(';');
				biasCorrection=biasInfo[0];
				idListBias=biasInfo[1].substr(1,1);
				useListBias=biasInfo[2];
			}
			else{
				biasCorrection=param[1];
			}
		}
		if (param[0] == 'NORMALIZATION_METHOD'){
			normalizationMethod=param[1];
		}
		if (param[0] == 'MIN_INF_RATIOS' && param[1] == 1){
			document.getElementsByName('minInfRatios')[0].checked=true;
		}
		if (param[0] == 'FDR_CONTROL'){
			var selFDR=document.getElementsByName('fdrControl');
			for (let i=0; i<selFDR[0].options.length; i++) {
				if (selFDR[0].options[i].value==param[1]) {
					selFDR[0].selectedIndex=i;
					break;
				}
			}
		}
/* Too risky for auto-selected by template (PP)
		if (param[0] == 'RESIDUAL_VAR'){
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

	// algorithm selection
	var selAlgo=document.getElementById('algoType');
	if (!algoType){
		algoType=(software && software[0].match('DIA'))? 'MSstats' : 'SSPA';
	}
	var algoTypeStrg='';
	if (algoType == 'SSPA'){
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
	if (algoType == 'PEP_INTENSITY' \|\| algoType == 'PEP_RATIO'){
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
	else if (algoType == 'SSPA'){
		var selPepFocus=document.getElementsByName('pepFocus');
		for (let i=1; i<selPepFocus[0].options.length; i++) {
			if (selPepFocus[0].options[i].value==pepCharge) {
				selPepFocus[0].selectedIndex=i;
				break;
			}
		}
	}

	//Bias correction
	if (biasCorrection == 'TRUE'){
		var selBiasCorrect=document.getElementById('biasCorrect');
		for (let i=1; i<selBiasCorrect.options.length; i++) {
			if (selBiasCorrect.options[i].value==normalizationMethod) {
				selBiasCorrect.selectedIndex=i;
				break;
			}
		}
		updateBias(normalizationMethod);
		if (idListBias){
			var selBiasrefProt=document.getElementsByName('refProt');
			for (let i=1; i<selBiasrefProt[0].options.length; i++) {
				if (selBiasrefProt[0].options[i].value==idListBias) {
					selBiasrefProt[0].selectedIndex=i;
					break;
				}
			}
		}

		if (useListBias){
			document.getElementsByName('refProtSelType')[0].selectedIndex=(useListBias == 'not') ? '1' : '0';
		}
	}

	// Target
	if (ptmID){
		var selPtmScope=document.getElementsByName('ptmScope');
		for (let i=1; i<selPtmScope[0].options.length; i++) {
			if (selPtmScope[0].options[i].value==ptmID) {
				selPtmScope[0].selectedIndex=i;
				break;
			}
		}
		updateRefQuantification(ptmID);
		updateTopN();
		if (ambiguousPos){
			document.getElementsByName('ambiguousPos')[0].checked=true;
		}
		else if (prsPos){
			document.getElmentsByName('badPRS')[0].selectedIndex=(prsPos == 'ambiguous')? 1 : 0;
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

function updateRefQuantification(ptmID) {
	var modQuantVis;
	var cusPTMs=document.selAnaForm.customPTM;
	//Clear any previous disabling
	if (cusPTMs.length) {
		for (var i=0; i<cusPTMs.length; i++) {
			cusPTMs[i].disabled=false;
		}
	}
	else {cusPTMs.disabled=false;}
	//Update with new selection
	if (ptmID==0) { // do nothing
		modQuantVis='hidden';
	}
	else {
		modQuantVis='visible';
		if (cusPTMs.length) { // disable selected quantified PTM
			for (var i=0; i<cusPTMs.length; i++) {
				cusPTMs[i].disabled=(cusPTMs[i].value==ptmID)? true : false;
			}
		}
		else {cusPTMs.disabled=true;} // must be =ptmID
	}
	document.getElementById('phosphoSPAN').style.display='none';
	document.getElementById('ptmSPAN').style.display='none';
	if (ptmID != 0) { // cannot be SSPA
		if (ptmID==$phosphoID && !document.selAnaForm.algoType.value.match('MSstats')) {
			document.getElementById('phosphoSPAN').style.display='';
		}
		else {
			document.getElementById('ptmSPAN').style.display='';
		}
	}
}
function updateTopN() {
	var myForm=document.selAnaForm;
	var algoType=myForm.algoType.value;
	if (myForm.labeling.value == 'FREE' && algoType.match('PEP_INTENSITY\|TnPQ') && myForm.ptmScope.value*1 == 0) {
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
	if (document.getElementById('algoType').value.match('^Ratio\|MSstats')) {
		document.getElementById('refProtDIV').style.visibility=(bias.match('REF_PROT\|globalStandards'))? 'visible' : 'hidden';
	}
	else {
		document.getElementById('refProtDIV').style.visibility=(bias && !bias.match('^none\|quantile'))? 'visible' : 'hidden';
	}
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
		propagScope=(force)? 'all' : document.getElementById('propagQuantScope:'+selState).value; // all,state,none
	if (propagScope=='none') {return;}
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
}

if ($quantifType eq 'XICMCQ') {
	print qq
|function updateRefAna(reference){
   var anaBox=document.selAnaForm.anaList;
   // Uncheck last reference analysis
   if (refValue > 0 && anaBox.length) {
		for (var i=0; i < anaBox.length; i++){
		    if(anaBox[i].value==refValue){
			    anaBox[i].checked=false;
		    }
		}
   }
   // Check reference analysis
   if(anaBox.length) {// more than 1 checkboxes
		for (var i=0; i < anaBox.length; i++){
		    if(anaBox[i].value==reference.value){
			    anaBox[i].checked=true;
			    refValue=anaBox[i].value;
		    }
        }
   }
   else{
		anaBox.checked=true;
   }

}
function updateSettings(act) {
	if (act=='more') {
		document.getElementById('moreSettings').style.display='none';
		document.getElementById('lessSettings').style.display='';
		document.getElementById('advancedSetDIV').style.display='';
	}
	else {
		document.getElementById('moreSettings').style.display='';
		document.getElementById('lessSettings').style.display='none';
		document.getElementById('advancedSetDIV').style.display='none';
	}
}
function showMZRange(){
	var chButton=document.getElementById('allChargeStates');
	if(chButton.checked){
		document.getElementById('massrangeDIV').style.display='';
	}
	else{
		document.getElementById('massrangeDIV').style.display='none';
	}
}
function showParameters(alignType) {
	if (alignType=='OBI') {
		document.getElementById('paramObiwarp').style.display='';
		document.getElementById('paramMs2').style.display='none';
	}
	else if (alignType=='MS2'){
		document.getElementById('paramObiwarp').style.display='none';
		document.getElementById('paramMs2').style.display='';
	}
}
function updateLabeling(labeling) {
	if (labeling=='SILAC'){
		document.getElementById('paramLabel').style.display='';
		document.getElementById('paramAlign').style.display='none';
	}
	else if (labeling=='FREE'){
		document.getElementById('paramLabel').style.display='none';
		document.getElementById('paramAlign').style.display='';
	}
	//Filtering Reference Analysis
	var refAnaSelect=document.getElementById('refAna_alignment');
	for (let i=1; i < refAnaSelect.options.length; i++) { //[0]= '-= Select =-'
		refAnaSelect.options[i].disabled=(refAnaSelect.options[i].getAttribute('data-labeling')==labeling)? false : true;
	}
	if (refAnaSelect.options[refAnaSelect.selectedIndex].disabled) {refAnaSelect.selectedIndex=0;} //deselect if disabled
	//Filtering Checkable analyses
	var anaBox=document.selAnaForm.anaList;
	if (anaBox.length) {
		for (let i=0; i < anaBox.length; i++) {
			var disStatus=(anaBox[i].getAttribute('data-labeling')==labeling)? false : true;
			anaBox[i].disabled=disStatus;
			document.getElementById('file_'+anaBox[i].value).disabled=disStatus;
		}
	}
	else {
		var disStatus=(anaBox.getAttribute('data-labeling')==labeling)? false : true;
		anaBox.disabled=disStatus;
		document.getElementById('file_'+anaBox.value).disabled=disStatus;
	}
}
function getRadioVal() {
  var rads = document.getElementsByName('extractL');

  for(var rad in rads) {
    if(rads[rad].checked)
      return rads[rad].value;
  }

  return null;
}
|;
}
elsif ($quantifType eq 'DESIGN') {
	print "var labelTypes=['",join("','",keys %{$designInfo{'LABEL'}}),"'];\n";
	print qq
|var algoNorms= {
	RatioTnPQ:{SILAC:[['Scale normalization*','SCALE'],['Reference protein(s)','REF_PROT']]}, // text:value pairs
	SxxxRatio:{
				SILAC:[
						['$normalizationNames{"median.none"}','median.none'],
						['$normalizationNames{"median.scale"}*','median.scale'],
						['$normalizationNames{"loess.none"}','loess.none'],
						['$normalizationNames{"loess.scale"}','loess.scale'],
						['$normalizationNames{"quantile"}','quantile']
				],
				FREE:[
						['$normalizationNames{"median.none"}','median.none'],
						['$normalizationNames{"median.scale"}*','median.scale'],
						['$normalizationNames{"quantile"}','quantile']
				]
	},
	MSstats:{FREE:[['Equalize medians*','equalizeMedians'],['Quantile','quantile'],['Global standards','globalStandards']]}
};
algoNorms.RatioTnPQ.FREE=algoNorms.RatioTnPQ.SILAC; algoNorms.RatioTnPQ.ITRAQ=algoNorms.RatioTnPQ.SILAC; algoNorms.RatioTnPQ.TMT=algoNorms.RatioTnPQ.SILAC;
//algoNorms.TnPQ=algoNorms.Ratio;
//algoNorms.SuperRatio=algoNorms.SimpleRatio;
algoNorms.SxxxRatio.ITRAQ=algoNorms.SxxxRatio.SILAC; algoNorms.SxxxRatio.TMT=algoNorms.SxxxRatio.SILAC;
algoNorms.MSstats.SILAC=algoNorms.MSstats.FREE;

function updateAlgoTypeSelection(algoType) {
	var myForm=document.selAnaForm;

	/* Update States display */
	var referenceText='';
	if (algoType=='select') {referenceText='';}
	if (algoType.match('S.+Ratio')) { // myProMS v2+ or MSstats
		referenceText=($nbCond==2)? '&nbsp;&lt;&lt;&lt;&nbsp;Reference state' : (algoType.match(':multiRef'))? '&nbsp;Each state will be used as reference for all following states' : '&nbsp;&lt;&lt;&lt;&nbsp;Common reference';
	}
	else if (algoType != 'select') { // old algos
		referenceText=(algoType=='SSPA')? '&nbsp;: At least 3 replicates required per state' : '&nbsp;&lt;&lt;&lt;&nbsp;Reference state'; //'&nbsp;Each state will be used as reference for all following states'; //
	}
	document.getElementById('refCommentSpan').innerHTML=referenceText;
	//Restrict old algo to 2 states
	var disabStatus=(algoType.match(':Ratio'))? true : false; // Ratio or TnPQ (old algos)
	for (let i=3; i <= $nbCond; i++) {
		myForm['state_'+i].disabled=disabStatus;
	}

	/* Update TopN */
	updateTopN();

	/* Multi-quantif option */
	var [multiQVis,multiQDisab]=($nbCond >= 3 && algoType.match(':singleRef'))? ['block',false] : ['none',true];
	document.getElementById('multiQuantifSPAN').style.display=multiQVis;
	var chkMultiQ=document.getElementById('multiQuantifCHK');
	chkMultiQ.disabled=multiQDisab;
	document.getElementById('multiQuantifDIV').style.display=(multiQVis=='block' && chkMultiQ.checked)? 'block' : 'none';

	/* Update peptide selection & PTM quantif options */
	var dispRatioPepDiv='none',dispSspaPepDiv='none';
	if (algoType=='SSPA') {
		dispSspaPepDiv='';
		if (myForm.ptmScope.value != 0) {
			myForm.ptmScope.selectedIndex=0;
			for (let i=1; i<myForm.ptmScope.options.length; i++) { // disable PTM quantif
				myForm.ptmScope.options[i].disabled=true;
			}
		}
	}
	else {
		if (!algoType.match('select\|MSstats')) {dispRatioPepDiv='';}
		for (let i=1; i<myForm.ptmScope.options.length; i++) { // enable PTM quantif
			myForm.ptmScope.options[i].disabled=false;
		}
	}
	document.getElementById('ratioPepDIV').style.display=dispRatioPepDiv;
	document.getElementById('sspaPepSPAN').style.display=dispSspaPepDiv;
	updateRefQuantification(myForm.ptmScope.value); // because of SSPA & MSstats constaints

	/* Update normalization options */
	var normSelOpt=myForm.biasCorrect.options;
	normSelOpt.length=2;
//	if (algoType=='select') {return;}
	//var trueAlgo=(algoType.match('SimpleRatio'))? 'SimpleRatio' : (algoType.match('MSstats'))? 'MSstats' : algoType;
	var trueAlgo=(algoType.match('PEP_'))? 'SxxxRatio' : (algoType.match('MSstats'))? 'MSstats' : (algoType.match(':Ratio'))? 'RatioTnPQ' : algoType;
	if (algoNorms[trueAlgo]) {
		myForm.biasCorrect.disabled=false;
		var labelType=myForm.labeling.value;
		for (var i=0; i<algoNorms[trueAlgo][labelType].length; i++) {
			normSelOpt[i+2]=new Option(algoNorms[trueAlgo][labelType][i][0],algoNorms[trueAlgo][labelType][i][1]);
		}
	}
	else { // SSPA, select
		myForm.biasCorrect.disabled=true;
	}

	/* Update remaining quantification settings */
	document.getElementById('refProtDIV').style.visibility='hidden'; // hide norm prot list in case visible
	var disabResVar=true, visResVarSpan='hidden';
	if (algoType.match('PEP_')) {
		disabResVar=false;
		if (myForm.residualVar.value=='biological') {visResVarSpan='visible';}
	}
	myForm.residualVar.disabled=disabResVar;
	document.getElementById('resVarSPAN').style.visibility=visResVarSpan;
	myForm.fdrControl.disabled=(algoType=='select')? true : false; // Now active for all algos (old algos too if only FDR & BH)
	myForm.minInfRatios.disabled=(algoType.match('MSstats\|SSPA'))? true : false; // No control on INF ratio before MSstats & SSPA
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
				//if (i==2 \|\| i==5) continue; // disable 'MSstats (SWATH)'
				//if (i==4 && labelType != 'SILAC') continue; // disable Super ratio

				if (algoSelOpt[i].value.match('MSstats')) continue; // disable MSstats
				if (algoSelOpt[i].value.match('PEP_RATIO') && labelType == 'FREE') continue; // disable PEP_RATIO algo

				algoSelOpt[i].disabled=false; // enable anything else
			}
		}
		if ($quantifMethods{DIA}) {
			//if ($nbCond > 2) {algoSelOpt[2].disabled=false;} // enable 'All vs All > MSstats (SWATH)'
			//algoSelOpt[5].disabled=false; // enable 'Common ref > MSstats (SWATH)'
			//algoSelOpt[8].disabled=false; // enable 'Peptide count > SSP Analysis'
			for (let i=1; i<algoSelOpt.length; i++) {
				if (algoSelOpt[i].value.match('MSstats')) {algoSelOpt[i].disabled=false;} // enable MSstats
				else if (algoSelOpt[i].value=='SSPA') {algoSelOpt[i].disabled=false;} // enable SSPA
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
		for (let i=1; i<selQuanti.options.length; i++) { // skip '-= Select =-'
			if (selQuanti.options[i].value==optList[0].value) {
				selQuanti.selectedIndex=i;
				break;
			}
		}
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
|;
}

if ($quantifType ne 'DESIGN') {
	print qq
|function testCheckbox() {
	var anaBox=document.selAnaForm.anaList;
	if (!anaBox) return 0; // no selectable analyses
	var selected=0;
	if (anaBox.length) { // more than 1 checkboxes
		for (i=0; i < anaBox.length; i++){
			if (!anaBox[i].disabled && anaBox[i].checked==true) {
				if (!testQuantifSelection(anaBox[i].value)) {
					selected=-1;
					break;
				}
				selected=1;
			}
		}
	}
	else if (!anaBox.disabled && anaBox.checked==true){
		selected=(!testQuantifSelection(anaBox.value))? -1 : 1;
	}
	return selected;
}
function checkall(checkStatus){
	var anaBox=document.selAnaForm.anaList;
	if (!anaBox) return; // no selectable analyses
	if (anaBox.length) { // more than 1 checkboxes
		for (i=0; i < anaBox.length; i++){
			anaBox[i].checked=checkStatus;
			updateQuantifSelection(anaBox[i]);
		}
	}
	else { // Only 1 checkbox
		anaBox.checked=checkStatus;
		updateQuantifSelection(anaBox);
	}
}
function testQuantifSelection(anaID) {return true;}
function updateQuantifSelection(chkbox) {}
|;
}
print qq
|function cancelAction() {
	//top.promsFrame.selectedAction='summary'; // set default action to 'summary'
	top.promsFrame.optionFrame.selectOption();
}

function checkForm(myForm) {
|;
if ($quantifType ne 'DESIGN') {
	print qq
|	if (testCheckbox()==0) {alert('ERROR: No Analyses selected!'); return false;}
	else if (testCheckbox()==-1) {alert('ERROR: Peptide quantification not selected for at least 1 Analysis!'); return false;} // only for SILAC,ITRAQ
|;
}
if ($quantifType eq 'XICMCQ') {
	print qq
|	var anaBox=myForm.anaList;
	var missingFile=false;
	if (anaBox.length) { // more than 1 checkbox
		for (let i=0; i < anaBox.length; i++){
			if (!anaBox[i].disabled && anaBox[i].checked==true && !document.getElementById('file_'+anaBox[i].value).value.match('mzX*ML')) {
				missingFile=true;
				break;
			}
		}
	}
	else if (!anaBox.disabled && anaBox.checked==true && !document.getElementById('file_'+anaBox.value).value.match('mzX*ML')) {
		missingFile=true;
	}
	if (missingFile) {
		alert('ERROR: Missing or not-mzXML/mzML file detected.');
		return false;
	}
	if (!myForm.refAna_alignment.value && getRadioVal()=='NO'){
		alert('ERROR: No reference selected for alignment.');
		return false;
	}
	if (!myForm.refAna_alignment.value){
		alert('ERROR: No reference selected.');
		return false;
	}
|;
}
elsif ($quantifType eq 'DESIGN') {
	print qq
|	//name
	if (!myForm.quantifName.value) {alert('ERROR: Missing name for quantification!'); return false;}
	//PhosphoRS settings
	if (myForm.ptmScope.value==$phosphoID && (myForm.okPRS.value<0 \|\| myForm.okPRS.value>100)) {
		alert('ERROR: PhosphoRS probability must be between 0 and 100%');
		return false;
	}
	//PTM filter
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

	// Labeling
	if (!myForm.labeling.value) {
		alert('ERROR: No labeling method selected.');
		return false;
	}
	var isLabeled=(myForm.labeling.value=='FREE')? false : true;

	// Chosen Algorithm
	if ( myForm.algoType.value.match('select') ) {
		alert('ERROR: Select an algorithm!');
		return false;
	}

	// States & Observations selection
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
	// Bias correction
	if (!myForm.biasCorrect.value && myForm.algoType.value != 'SSPA') {
		alert('ERROR: No bias correction selected.');
		return false;
	}
	//Multi-quantification
	if (myForm.algoType.value.match(':singleRef') && myForm.multiQuantif.checked) {
		if (numCondition >= 3 && !myForm.quantifName.value.match(/%(TEST\|#)%/i) && !confirm('WARNING: All quantifications have same name. Proceed anyway?')) {
			return false;
		}
		if (!confirm((numCondition-1)+' quantifications will be launched. Proceed?')) {
			return false;
		}
	}
|;
}
print qq
|	return true;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" >
<CENTER>
<FONT class="title">$titleString</FONT>
$oldAlgoAccessStrg
<BR>
<FORM name="selAnaForm" method="post" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="ID" value="$branchID">
<INPUT type="hidden" name="quantifType" value="$quantifType">
<BR>
|;

##################################################
####>Displaying quantification-specific forms<####
##################################################
if ($quantifType eq 'DESIGN') { # Design-based quanti method
	#my $disabSuperRatio=($projectAccess eq 'bioinfo')? '' : ' disabled'; # temporary
	print qq
|<INPUT type="hidden" name="nbMaxCond" value="$nbCond">
<INPUT type="hidden" name="phosphoID" value="$phosphoID">
<TABLE bgcolor="$darkColor">
<TR><TH class="title3" align=right valign=top>Name :</TH><TD bgcolor="$lightColor"><INPUT type="text" name="quantifName" value="" placeholder="Name of quantification" class="title3" style="width:400px"/>
	<DIV id="multiQuantifDIV" style="display:none">&nbsp;<B>Use: "<FONT color="#DD0000">%TEST%</FONT>", "<FONT color="#DD0000">%REF%</FONT>" and  "<FONT color="#DD0000">%#%</FONT>" for test-state name, reference-state name and quantification rank respectively.</B></DIV>
</TD></TR>

	|;
	if ($userInfo[1] eq 'mass' || $userInfo[1] eq 'bioinfo'){
print qq |
<TR><TH class="title3" align=right valign=top>Use parameters :</TH><TD bgcolor="$lightColor" class="title3" nowrap>&nbsp;Template:<SELECT id="template" name="template" class="title3" onchange="useTemplate(this.value)">
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
print qq |
</SELECT></TD></TR>
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
	print qq
|</SELECT>
&nbsp;&nbsp;&nbsp;&nbsp;Algorithm:<SELECT id="algoType" name="algoType" class="title3 template" onchange="updateAlgoTypeSelection(this.value)">
<OPTION value="select">-= Select =-</OPTION>
<OPTGROUP label="Peptide ratio:">
<OPTION value="PEP_RATIO:SimpleRatio:singleRef" disabled>All vs State1</OPTION>
<OPTION value="PEP_RATIO:SuperRatio:multiRef" disabled>Super ratio</OPTION>
</OPTGROUP>
<OPTGROUP label="Peptide intensity:">
<OPTION value="PEP_INTENSITY:SimpleRatio:singleRef" disabled>All vs State1</OPTION>
<OPTION value="PEP_INTENSITY:SimpleRatio:multiRef" disabled>All vs All</OPTION>
</OPTGROUP>
<OPTGROUP label="DIA (MSstats):">
<OPTION value="MSstats:SimpleRatio:singleRef" disabled>All vs State1</OPTION>
<OPTION value="MSstats:SimpleRatio:multiRef" disabled>All vs All</OPTION>
<!--
<OPTGROUP label="2-state (old algos):">
<OPTION value="Ratio:Ratio:singleRef" disabled>Simple ratio (old)</OPTION>
<OPTION value="TnPQ:Ratio:singleRef" disabled>TnPQ</OPTION>
</OPTGROUP>
-->
<OPTGROUP label="Peptide count:">
<OPTION value="SSPA" disabled>SSP Analysis</OPTION>
</OPTGROUP>
</SELECT>

<SPAN id="topSPAN" class="template" style="display:none">&nbsp;&nbsp;&nbsp;&nbsp;Use:<SELECT name="matchingPep" class="title3 template" disabled><OPTION value="0" disabled>any</OPTION><OPTION value="1" selected>matching</OPTION></SELECT>&nbsp;top<SELECT name="topN" class="title3 template" disabled><OPTION value="1">1</OPTION><OPTION value="2">2</OPTION><OPTION value="3" selected>3</OPTION><OPTION value="0">N</OPTION></SELECT>&nbsp;peptides</SPAN>
<SPAN id="sspaPepSPAN" style="display:none" class="template">&nbsp;&nbsp;&nbsp;&nbsp;Use:<SELECT name="pepFocus" class="title3 template"><OPTION value="sp_count">Peptide/spectrum matches</OPTION><OPTION value="all_ion">All peptide ions</OPTION><OPTION value="all_pep">All peptides</OPTION><OPTION value="dist_ion">Distinct peptide ions</OPTION><OPTION value="dist_pep">Distinct peptides</OPTION><OPTION value="dist_seq">Distinct peptide sequences</OPTION></SELECT></SPAN>
</TD></TR>
<TR><TH class="title3" align=right valign=top>Target :</TH><TD bgcolor="$lightColor" class="title3" nowrap><SELECT name="ptmScope" class="title3 template" onchange="updateRefQuantification(this.value); updateTopN()">
<OPTION value="0">Whole proteins</OPTION>
|;
	foreach my $modID (sort{lc($projectVarMods{$a}) cmp lc($projectVarMods{$b})} keys %projectVarMods) {
		print "<OPTION value=\"$modID\">$projectVarMods{$modID}-proteins</OPTION>\n";
	}
	print qq
|</SELECT>
<SPAN id="phosphoSPAN" style="display:none" class="template">&nbsp;&nbsp;Positions with PhosphoRS probability&ge;<INPUT type="text" name="okPRS" class="title3 template" value="51" size="2">% are confirmed. The others are <SELECT name="badPRS" class="title3 template"><OPTION value="exclude">excluded</OPTION><OPTION value="ambiguous" selected>delocalized</OPTION></SELECT></SPAN>
<SPAN id="ptmSPAN" style="display:none" class="template">&nbsp;&nbsp;<INPUT type="checkbox" name="ambiguousPos" value="1" class="template">Positions of modification sites are not reliable</SPAN>
|;
#<!-- Hide reference quantification option for now (PP 17/12/14)
#&nbsp;&nbsp;&nbsp;&nbsp;Reference quantification:<SELECT name="refQuantifID" class="title3"><OPTION value="0">-= Select =-</OPTION>
#|;
#	foreach my $qID (sort{lc($refQuantifications{$a}) cmp lc($refQuantifications{$b})} keys %refQuantifications) {
#		print "<OPTION value=\"$qID\">$refQuantifications{$qID}</OPTION>\n";
#	}
#	print qq
#|</SELECT>&nbsp;&nbsp;
#-->
#<SPAN>
	print qq
|</TD></TR>
<TR><TH class="title3" align=right valign=top>States :</TH><TD bgcolor=$lightColor>
|;
	if ($nbCond > 5) {
		print qq
|<INPUT type="button" value="Use all states" onclick="autoSelectStates('all')">&nbsp;<B>or use state range:</B><INPUT type="text" name="stateRange" style="width:80px" value="" placeholder="e.g. 1,5-8"><INPUT type="button" value="Apply" onclick="autoSelectStates('range')">&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Clear all states" onclick="clearAllStates()">
|;
	}
	print qq
|<DIV style="max-height:150px;overflow:auto"><TABLE>
|;
	# print the states that would be chosen for ratio computation
	foreach my $y (1..$nbCond) {
		print "<TR><TH align=right nowrap>#$y:</TH><TD><SELECT name=\"state_$y\" onchange=\"updateStateSelection(this,$y)\">$selCondStrg1</SELECT>";
		print "<SPAN id=\"refCommentSpan\" class=\"template\" style=\"color:#DD0000;font-weight:bold\"></SPAN>" if $y==1;
		print "</TD></TR>\n";
	}

	print qq
|</TABLE></DIV><SPAN id="multiQuantifSPAN" class="template" style="display:none"><LABEL style="color:#DD0000;font-weight:bold"><INPUT type="checkbox" name="multiQuantif" class="template" id="multiQuantifCHK" value=1 onchange="document.getElementById('multiQuantifDIV').style.display=(this.checked)? 'block' : 'none';" disabled>Perform a separate quantification for each ratio</LABEL></SPAN>
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
|<TR><TH align=right nowrap valign=top>&nbsp;Peptide selection :</TH><TD bgcolor=$lightColor nowrap>
	<TABLE cellspacing=0 cellpadding=0><TR>
	<TD valign="top" nowrap>&nbsp;<B>Specificity:</B><SELECT name="pepSpecifity" class="template"><OPTION value="unique">Proteotypic</OPTION><OPTION value="unique_shared">Proteotypic + shared</OPTION><OPTION value="all" selected>All</OPTION></SELECT>
	&nbsp;&nbsp;&nbsp;<B>Missed cleav.:<B><SELECT name="pepCleavage" id="pepCleavage" class="template"><OPTION value="1">Allowed</OPTION><OPTION value="0" selected>Not allowed</OPTION></SELECT>
	&nbsp;&nbsp;&nbsp;<B>Modifications:<B>
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
	print qq |</DIV></TD>
	<TD valign="top" nowrap><DIV id="ratioPepDIV" class="template">&nbsp;&nbsp;&nbsp;<B>Charges:<B><SELECT name="pepCharge" class="template" onchange="updatePeptideSource(this.value)"><OPTION value="all">All</OPTION><OPTION value="best">Best signal</OPTION></SELECT>
&nbsp;&nbsp;&nbsp;<B><SUP>°</SUP>Sources:<B><SELECT name="pepSource" class="template" id="pepSource"><OPTION value="all">All</OPTION><OPTION value="best">Best signal</OPTION></SELECT>
</DIV>
</TD></TR></TABLE>
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
|</SELECT></TD></TR>
<TR><TH align=right nowrap valign=top>&nbsp;Quantification settings :
	</TH><TD bgcolor=$lightColor nowrap valign=top>
	<TABLE border=0 cellpadding=0 cellspacing=0><TR>
	<TH align=right valign="top">&nbsp;-Bias correction:</TH>
	<TD><SELECT name="biasCorrect" class="template" id="biasCorrect" onchange="updateBias(this.value)"><OPTION value="">-= Select =-</OPTION><OPTION value="none.none">None</OPTION></SELECT>&nbsp;</TD>
	<TD valign="top" nowrap><DIV id="refProtDIV" class="template" style="visibility:hidden">&nbsp;<B><SELECT name="refProtSelType" class="template"><OPTION value="use">Use</OPTION><OPTION value="not" disabled>Do not use</OPTION></SELECT>
data from</B> <SELECT name="refProt" class="template" onchange="var selType=document.selAnaForm.refProtSelType; if (this.value==0) {selType.selectedIndex=0; selType.options[1].disabled=true;} else {selType.options[1].disabled=false;}"><OPTION value="0">All proteins</OPTION>
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
	&nbsp;<B>-Biological replicates:</B><SELECT name="residualVar" class="template" onchange="document.getElementById('resVarSPAN').style.visibility=(this.value=='biological')? 'visible' : 'hidden'"><OPTION value="auto" selected>Auto-detect*</OPTION><OPTION value="technical">No</OPTION><OPTION value="biological">Yes</OPTION></SELECT>
<SPAN id="resVarSPAN" style="visibility:hidden" class="template"><FONT style="font-weight:bold;color:#DD0000"> Requires at least 2 biological replicates in each state</FONT></SPAN><BR>
	&nbsp;<B>-p-value correction:</B><SELECT name="fdrControl" class="template"><OPTION value="none">None</OPTION><OPTION value="fdr" selected>Benjamini-Hochberg*</OPTION><OPTION value="bonferroni">Bonferroni</OPTION></SELECT><BR>
	&nbsp;<B>-</B><INPUT type="checkbox" name="minInfRatios" class="template" value="1"><B>Avoid infinite ratios when possible.</B><BR>
	&nbsp;<SMALL><SUP>*</SUP>Recommanded options</SMALL>
	</TD></TR>
</TABLE>
|;
}
elsif ($quantifType eq 'XICMCQ') {
	print qq
|<TABLE bgcolor=$darkColor>
<TR><TH class="title2" align=right>Name :</TH><TD bgcolor=$lightColor><INPUT type="text" name="quantifName" value="$quantifProcesses{$quantifType} extraction" class="title3" style="width:400px"/></TD></TR>
<TR>
<TH align=right nowrap valign=top>Raw-data settings :</TH>
<TD bgcolor=$lightColor nowrap>&nbsp;<B>Extraction type:</B>
<SELECT name="rawdata_extraction">
<OPTION value="centroid">Centroid</OPTION>
<OPTION value="profile">Profile</OPTION>
</SELECT>
<FONT class="font11">&nbsp;(for mzXML)</FONT>
</TD>
</TR>
<TR>
<TH align=right nowrap valign=top>Isotope labeling :</TH>
<TD bgcolor=$lightColor nowrap>&nbsp;<B><SELECT name="extractL" value="YES" onchange="updateLabeling(this.value)"><OPTION value=\"FREE\">NONE</OPTION><OPTION value=\"SILAC\">SILAC</OPTION></SELECT></B>
|;
	###> Label informations
	my $optMods="<OPTION value=\"\">-= Select =-</OPTION>\n";
	$optMods.="<OPTION value=\"-1\"> Light / +0.0000 Da</OPTION>\n";
	foreach my $modID (keys %anaLabelMods) {
		my ($anaSpec,$psiMsName,$interimName,$isValid,$monoMass)=@{$anaLabelMods{$modID}};
		my $isDisabled=($isValid)? "" : "disabled";
		my $monoMassStg=($monoMass>0)?  "+$monoMass Da" : "$monoMass Da";
		$optMods.="<OPTION value=\"$modID\" $isDisabled>$psiMsName / $monoMassStg</OPTION>\n";
	}
	#print "<TH align=right nowrap valign=top></TH><TD bgcolor=$lightColor nowrap><B>\n";
	print "<DIV id=\"paramLabel\" style=\"display:none\"><B>\n";
	#my @colors=("#80ACFF","#5195FF","#428CFF");
	foreach my $channel (1..3) {
		print "<BR>\n" if $channel > 1;
		print "<FIELDSET style=\"border:3px groove threedface;\"><LEGEND>Name for channel #$channel:<INPUT type=\"text\" id=\"name_channel$channel\" name=\"name_channel$channel\" value=\"\" style=\"width:250px\"></LEGEND>\n";
		#my $color=$colors[$channel-1];
		foreach my $label (1..$maxLabPerChannel) { # In triplex SILAC,channel can be composed of several labels: Heavy=Arg10 and Lys8 ; Medium=
			my $state=($label==1)? '' : 'none';
			print qq
|<DIV id="fs:$channel:$label" style="display:$state">
<FIELDSET>
<LEGEND>Quantification Label:</LEGEND>
<TABLE width=100% bgcolor=$darkColor>
<TR><TH align=right nowrap valign=middle>Label Name:</TH><TD><INPUT type="text" id="lab$label:name_channel$channel" name="lab$label:name_channel$channel" value="" style="width:250px"></TD></TR>
<TR><TH align=right nowrap valign=middle>Modification target:</TH><TD><SELECT id="lab$label:tg_channel$channel" name="lab$label:tg_channel$channel" onchange="updateTextState('lab$label:sp_channel$channel',this.value)"><OPTION value="ANY">Side chain</OPTION><OPTION value="Nter">N-Ter</OPTION><OPTION value="Cter">C-Ter</OPTION></SELECT></TD></TR>
<TR><TH align=right nowrap valign=middle>Modification:</TH><TD><SELECT id="lab$label:channel$channel" name="lab$label:channel$channel" onchange="updateLabels(this.value,$channel,$label)">$optMods</SELECT> on <INPUT type="text" id="lab$label:sp_channel$channel" name="lab$label:sp_channel$channel" value="" size="2"></TD></TR>
</TABLE>
</FIELDSET>
<INPUT id="show:$channel:$label" type="button" value="Add quantification label" onclick="addQuantificationLabel($channel,$label,'show')" style="font-size:11px"/><INPUT id="hide:$channel:$label" type="button" value="Remove quantification label" onclick="addQuantificationLabel($channel,$label,'hide')" style="display:none; font-size:11px;"/>
</DIV>
|;
		}
		#print "<TR><TD nowrap>&nbsp;Quan Channel $channel:<SELECT id=\"channel$channel\" name=\"channel$channel\" onchange=\"updateLabels(this.value,$channel)\">$optMods</SELECT>\n";
		#print "<";
		#print "&nbsp;Channel name:<INPUT type=\"text\" id=\"name_channel$channel\" name=\"name_channel$channel\" value=\"\" size=\"6\">";
		#print "&nbsp;Delta-Mass: <INPUT type=\"text\" id=\"dm_channel$channel\" name=\"dm_channel$channel\" value=\"\" size=\"6\">";
		#print "&nbsp;Specificity: <INPUT type=\"text\" id=\"sp_channel$channel\" name=\"sp_channel$channel\" value=\"\" size=\"2\"></TD></TR>\n";
		print "\n</FIELDSET>\n";
	}
	print "</B></DIV>";
	print "</TD></TR>\n";
	print qq
|<DIV id="paramAlign">
<TR>
<TH align=right nowrap valign=top>Alignment settings :</TH><TD bgcolor=$lightColor><TABLE cellpadding=0>
	<TR><TD nowrap>&nbsp;<B>Alignment algorithm:</B>
	<SELECT name="alignment_method" onchange="showParameters(this.value)">
	<OPTION value="MS2">ms2</OPTION>
	<OPTION value="OBI">OBI-Warp</OPTION>
	</SELECT>
	&nbsp;&nbsp;&nbsp;<B>Reference:</B>
	<SELECT name="refAna_alignment" id="refAna_alignment" required">
	<OPTION value="">-= Select =-</OPTION>
|;
	foreach my $refAnaData (@itemAnalyses) {##> Option of reference for alignment
		my ($anaID,$valStat,$msType,$dataFile,$fileFormat,$wiffFile,$taxonomy,$maxRank,$minScore,$instrument,$refDbUsed,$labelStrg,@projHierarchy)=@{$refAnaData};
		my $disabStrg=($anaLabeling{$anaID} eq 'FREE')? '' : ' disabled';
		print "	<OPTION value=\"$anaID\" data-labeling=\"$anaLabeling{$anaID}\" onclick=\"updateRefAna(this);\"$disabStrg>$projHierarchy[-1]</OPTION>\n" unless $valStat<1;
	}
	print qq
|	</SELECT>
	</TD>
	</TR>
	<TD bgcolor=$lightColor nowrap valign=top>
		<DIV id="paramObiwarp" style="display:none">
		<TABLE cellpadding=0 cellspacing=0>
		<TR>
		<TD nowrap>&nbsp;<B>Align from <INPUT type="text" name="mz_start" value="400" size="2"> to <INPUT type="text" name="mz_stop" value="1200" size="3"> m/z window</B></TD>
		</TR>
		</TABLE>
		</DIV>
		<DIV id="paramMs2">
		<TABLE cellpadding=0 cellspacing=0>
		<TR><TD nowrap>&nbsp;<B>Tendency <INPUT type="text" name="ms2_tendency" value="10" size="2"></TR>
		<TR><TD nowrap>&nbsp;<INPUT type="checkbox" name="ms2Smooth" value="1" checked><B>MS2 smoothing <INPUT type="text" name="ms2_smoothing" value="5" size="1"></TD></TR>
		<TR><TD nowrap>&nbsp;<INPUT type="checkbox" name="ms1Smooth" value="1" checked><B>MS1 smoothing <INPUT type="text" name="ms1_smoothing" value="3" size="1"></B></TD>
		</TR>
		</TABLE>
		</DIV>
	</TD>
</TABLE>
</TD>
</TR>
</DIV>
<TR>
<TH align=right nowrap valign=top>Peptide selection :</TH>
<TD bgcolor=$lightColor nowrap>
	<TABLE cellpadding=0 cellspacing=0>
	<TD bgcolor=$lightColor nowrap><INPUT type="checkbox" name="allChargeStates" id="allChargeStates" value="1" onchange="showMZRange()"><B>Extract all charge states of the peptides (even if no MS/MS exists for it)</B>&nbsp;
	</TABLE>
	<DIV id="massrangeDIV" style="display:none">
	<TABLE cellpadding=0 cellspacing=0>
	<TR>
	<TD nowrap>&nbsp;<B>Range (m/z) <INPUT type="text" name="mzRange_start" value="400" size="2"> - <INPUT type="text" name="mzRange_stop" value="1200" size="3"></B></TD>
	</TR>
	</TABLE>
	</DIV>
</TD>
</TR>
<TR>
<TH align=right nowrap valign=top>Extract XIC traces :</TH>
<TD bgcolor=$lightColor nowrap>&nbsp;<B>No<INPUT type="radio" name="traces" id="traces" value="0" checked>Yes<INPUT type="radio" name="traces" id="traces" value="1"></B>
</TD>
</TR>
<TR>
<TH align=right nowrap valign=top>&nbsp;Quantification settings :</TH><TD bgcolor=$lightColor nowrap valign=top>
	&nbsp;<B>Type of XIC:</B><SELECT name="XIC_type" onchange="updateChargeState(this.value)"><OPTION value="sum">TIC XIC</OPTION><OPTION value="max">BasePeak XIC</OPTION></SELECT>
	<BR><INPUT type="button" id="moreSettings" class="font11" value="More settings" onclick="updateSettings('more')"/><INPUT type="button" id="lessSettings" class="font11" value="Less settings" style="display:none" onclick="updateSettings('less')"/>
	<DIV id="advancedSetDIV" style="display:none">
	&nbsp;&bull;<B><U>Advanced settings</U>:</B>
	<TABLE cellpadding=0 cellspacing=0>
	<TR>
	<TD nowrap>&nbsp;<B>Size of mass tolerance window for XIC:</B> min=<INPUT type="text" name="mztol_min" value="5" size="1"> max=<INPUT type="text" name="mztol_max" value="5" size="1"></B><SELECT name="XIC_range"><OPTION value="ppm">ppm</OPTION><OPTION value="mz">mz</OPTION></SELECT></TD>
	</TR>
	<TR><TD nowrap>&nbsp;<B>Peak matching for XIC:</B>
	<SELECT name="XIC_rt">
	<OPTION value="post_matching">Post matching mode</OPTION>
	<OPTION value="real_or_mean">Real or mean mode</OPTION>
	<OPTION value="mean">Mean mode</OPTION>
	</SELECT>
	</TD>
	</TR>
	<TR>
	<TD nowrap>&nbsp;<B>Detection threshold between <INPUT type="text" name="dt_start" value="30000" size="4"> to <INPUT type="text" name="dt_stop" value="50000" size="4"></B></TD>
	</TR>
	<TR>
	<TD>&nbsp;&nbsp;<B><U>XIC filtering</U>:</B></TD>
	</TR>
	<TR><TD bgcolor=$lightColor nowrap><INPUT type="checkbox" name="asFilter" value="1" checked>&nbsp;<B>Anti-Spike:</B><INPUT type="text" name="anti_spike" value="5" size="1"></TD><TR>
	<TR><TD bgcolor=$lightColor nowrap><INPUT type="checkbox" name="hmFilter" value="1">&nbsp;<B>Half-Mediane:</B> min=<INPUT type="text" name="med_min" value="5" size="1"> max=<INPUT type="text" name="med_max" value="40" size="1"></TD><TR>
	<TR><TD bgcolor=$lightColor nowrap><INPUT type="checkbox" name="smFilter" value="1">&nbsp;<B>Smoothing:</B><INPUT type="text" name="smooth_val" value="3" size="1"></TD><TR>
	</TABLE>
	</DIV>
	</TD>
</TR>
</TABLE>
<BR>
|;
}

if ($quantifType ne 'DESIGN') {
	print qq
|<BR><TABLE border=0 cellspacing=0 cellpadding=0>
<TR bgcolor="$darkColor">
	<TH class="rbBorder"><INPUT type="checkbox" onclick="checkall(this.checked)"></TH>
	<TH class="rbBorder" nowrap colspan=2>&nbsp;Analysis&nbsp;</TH>
|;
	print "<TH class=\"rbBorder\">&nbsp;Peptide quantification&nbsp;</TH>\n";# if $quantifType=~/SILAC|ITRAQ|TMT/;
	print qq
|	<TH class="rbBorder">&nbsp;MS type&nbsp;<BR>& File</TH>
	<TH class="rbBorder">&nbsp;Instrument&nbsp;</TH>
	<TH class="rbBorder">&nbsp;Search file<BR>&nbsp;& Engine&nbsp;</TH>
	<TH class="rbBorder">&nbsp;Databank(s)&nbsp;<BR>&nbsp;Taxonomy&nbsp;</TH>
	<TH class="rbBorder" nowrap>&nbsp;Min. score&nbsp;<BR>&nbsp;Max. rank&nbsp;</TH>
	<TH class="bBorder">&nbsp;Validated&nbsp;<BR>&nbsp;proteins&nbsp;</TH>
</TR>
|;
}

my %itemIcones=&promsConfig::getItemIcones;
my $bgColor=($ITEM eq 'SAMPLE' || $ITEM eq 'SPOT')? $lightColor : $darkColor;
my %prevItemName;
my $disabSubmit=' disabled';
if ($quantifType ne 'DESIGN') {
	my %filesList;
	my $mzXMLPath="$promsPath{tmp}/upload/project_$projectID";

	if (-e $mzXMLPath){
		opendir (DIR, $mzXMLPath) || print "ERROR: Unable to read '$mzXMLPath' !<BR>\n";
		while (defined (my $currentmzXMLFile = readdir (DIR))) {
			next unless ( $currentmzXMLFile =~ /.+\.mzXML\Z/ || $currentmzXMLFile =~ /.+\.mzML\Z/ );
			$filesList{$currentmzXMLFile}=0;
		}
		closedir DIR;
	}

	foreach my $refAnaData (@itemAnalyses) {
		my ($anaID,$valStat,$msType,$dataFile,$fileFormat,$wiffFile,$taxonomy,$maxRank,$minScore,$instrument,$refDbUsed,$labelStrg,@projHierarchy)=@{$refAnaData};
		$taxonomy='Unknown' unless $taxonomy;
		$taxonomy=~s/\(.*\)//;
		#my $okQuantifAna=($valStat>=1 && $msType ne 'PMF' && (($quantifType eq 'EMPAI' && ($fileFormat eq 'MASCOT.DAT' || $fileFormat eq 'MASCOT.PDM')) || ($quantifType eq 'SIN' && $fileFormat ne 'PHENYX.XML') || $quantifType eq 'XIC'))? 1 : 0;
		my $okQuantifAna=0;
		if ($valStat>=1 && $msType ne 'PMF') {
			if (($quantifType eq 'EMPAI' && ($fileFormat eq 'MASCOT.DAT' || $fileFormat eq 'MASCOT.PDM'))) {$okQuantifAna=1;}
			elsif ($quantifType eq 'SIN' && $fileFormat eq 'MASCOT.DAT') {$okQuantifAna=1;}
			#elsif ($quantifType eq 'XIC') {$okQuantifAna=1;}
			elsif ($quantifType eq 'XICMCQ') {$okQuantifAna=1;}
			elsif ($labelStrg && $labelStrg=~/::$quantifType/) { # SILAC, iTRAQ
				#foreach my $refQuantif (@{$anaPeptideQuantifs{$anaID}}) {
				#	if ($refQuantif->[3]==1) {
				#		$okQuantifAna=1;
				#		last;
				#	}
				#}
			}
		}

		$msType=$msTypeName{$msType}; # global
		$fileFormat=~s/\..*//;
		$instrument='-' unless $instrument;
		$maxRank='-' unless $maxRank;
		$minScore='-' unless defined $minScore;
		$disabSubmit='' if $okQuantifAna; # at least 1 selectable analysis => activate Submit
		##>Row color
		my $fatherIt=$projHierarchy[-3];
		if ($fatherIt && (!$prevItemName{$fatherIt} || $prevItemName{$fatherIt} ne $projHierarchy[-2])) { # keep color if same analysis parent item (SAMPLE or SPOT)
			$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
		}
		elsif ($ITEM eq 'EXPERIMENT' || $ITEM eq 'GEL2D') {print "<TR bgcolor=$bgColor><TD colspan=2></TD><TD colspan=8><HR width=98%></TD></TR>\n";}
		print "<TR valign=middle bgcolor=$bgColor>\n";
		##>Checkbox
		#my $boxStr='-';
		#if ($okQuantifAna) {
		#	#$boxStr="<INPUT type=\"checkbox\" name=\"anaList\" value=\"$anaID";
		#	#$boxStr.=".$anaPeptideQuantifs{$anaID}[0]" if $quantifType=~/SILAC|ITRAQ|TMT/;
		#	#$boxStr.="\"\>";
		#	$boxStr="<INPUT type=\"checkbox\" name=\"anaList\" value=\"$anaID\"/>";
		#}
		#print "\t<TH valign=middle>$boxStr</TH>\n";

		##>Checkbox
		if ($okQuantifAna) {
			my $disabStrg=($quantifType eq 'XICMCQ' && $anaLabeling{$anaID} ne 'FREE')? ' disabled' : '';
			print "\t<TH valign=middle><INPUT type=\"checkbox\" name=\"anaList\" value=\"$anaID\" data-labeling=\"$anaLabeling{$anaID}\"$disabStrg/></TH>\n";
		}
		else {print "\t<TH valign=middle>-</TH>\n";}
		#>Parents & Analysis
		my $parentStrg='';
		for (my $i=0;$i<=$#projHierarchy-2;$i+=2) { # stops before ana name
			my $IT=$projHierarchy[$i];
			my $itName=$projHierarchy[$i+1];
			if ($prevItemName{$IT} && $prevItemName{$IT} eq $itName) {
				#$parentStrg.="<FONT style=\"visibility:hidden\">&nbsp;&nbsp;$itName&nbsp;&nbsp;&gt;&nbsp;&nbsp;</FONT>";
			}
			else {
				$parentStrg.="&nbsp;&nbsp;$itName&nbsp;&nbsp;&gt;&nbsp;&nbsp;";
				$prevItemName{$projHierarchy[$i]}=$itName;
				for (my $j=$i+2;$j<$#projHierarchy-1;$j+=2) {$prevItemName{$projHierarchy[$j]}='';}
			}
		}
		my $anaCode=($valStat==-1)? 'analysis:no_scan' : ($valStat==0)? 'analysis:no_val' : ($valStat==1)? 'analysis:part_val' : 'analysis:val';
		print qq
|	<TH nowrap align=left valign=middle>$parentStrg</TH>
	<TH nowrap align=left valign=middle><IMG src="$promsPath{images}/$itemIcones{$anaCode}">&nbsp;$projHierarchy[-1]&nbsp;</TH>
|;

		##>Labeling method
		print "<TD>";
		if (!$labelStrg) {print '&nbsp;None&nbsp;';}
		#elsif ($quantifType=~/SILAC|ITRAQ|TMT/ && $valStat>=1) { # No peptide quantification for non-validated analysis
		#	print "<SELECT id=\"quantif_$anaID\" name=\"quantif_$anaID\" style=\"width:250px\" disabled>\n";
		#	print "<OPTION value=\"\">-= Select =-</OPTION>\n" if scalar @{$anaPeptideQuantifs{$anaID}} > 1;
		#	foreach my $refQuantif (@{$anaPeptideQuantifs{$anaID}}) {
		#		my ($quantifID,$quantifName,$quantifInfoStrg,$status)=@{$refQuantif};
		#		print "<OPTION value=\"$quantifID\" onmouseover=\"popup('<B>Channels:</B><BR>$quantifInfoStrg')\" onmouseout=\"popout()\"";
		#		print ' disabled' if $status != 1;
		#		print ">$quantifName</OPTION>\n";
		#	}
		#	print "</SELECT>\n";
		#	#my $labelingName=(split('::',$labelStrg))[0];
		#	#print "&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<B>Name:</B> $anaPeptideQuantifs{$anaID}[1]<BR><B>Channels:</B><BR>$anaPeptideQuantifs{$anaID}[2]')\" onmouseout=\"popout()\">$labelingName</A>&nbsp;";
		#}
		else {print $labelStrg;}
		print "</TD>\n";

		my $dbStrg;
		foreach my $dbID (@{$refDbUsed}) {
			$dbStrg.='&bull;' if $dbStrg;
			$dbStrg.=$listDataBank{$dbID};
		}
		print qq
|	<TD class="center">$msType<BR>&nbsp;$wiffFile&nbsp;</TD>
	<TD class="center">&nbsp;$instrument&nbsp;</TD>
	<TD class="center">&nbsp;$dataFile&nbsp;<BR>&nbsp;$fileFormat&nbsp;</TD>
	<TD class="center">&nbsp;$dbStrg&nbsp;<BR>&nbsp;$taxonomy&nbsp;</TD>
	<TD class="center">&nbsp;$minScore&nbsp;<BR>&nbsp;$maxRank&nbsp;</TD>
|;
		##>Validated proteins
		if ($valStat>=1) {print "<TD class=\"center\">$anaProteins{$anaID}[0] ($anaProteins{$anaID}[1])</TD>\n";}
		else {print "<TD class=\"center\">$anaProteins{$anaID}[0] / $anaProteins{$anaID}[1]</TD>\n";}
		print "</TR>\n";
		if ($quantifType eq 'XICMCQ' && $okQuantifAna) {
			###> Pre-selection of mzXML file
			my $listmzXMLFiles='<OPTION value=\"\">-= Select =-</OPTION>';
			my $selected;
			my $msName='';
			($msName=$wiffFile)=~s/\.[^\.]+//;
			$msName=~ s/\s+\Z//;
			foreach my $dataFile (sort {lc($a) cmp lc($b)} keys %filesList) {
				#$bgColor=($bgColor eq $darkColor)? $lightColor : $darkColor;
				$selected=($dataFile=~/$msName/)? ' selected': '';
				$listmzXMLFiles.="<OPTION value=\"$dataFile\"$selected>$dataFile</OPTION>\n";
			}
			$listmzXMLFiles.="</SELECT>";
			my $disabStrg=($anaLabeling{$anaID} eq 'FREE')? '' : ' disabled';
			#print "<TR bgcolor=$bgColor><TD colspan=2></TD><TD colspan=8><B>mzXML file:<INPUT type=\"file\" name=\"file_$anaID\" id=\"file_$anaID\" value=\"\" size=80></TD></TR>\n";
			print "<TR bgcolor=$bgColor><TD colspan=2></TD><TD colspan=8><B>mzXML file:<SELECT name=\"file_$anaID\" id=\"file_$anaID\"$disabStrg>$listmzXMLFiles</TD></TR>\n";
			if ($ITEM ne 'EXPERIMENT' && $ITEM ne 'GEL2D') {
				$bgColor=($bgColor eq $darkColor)? $lightColor : $darkColor;
			}
		}
	}

	#my $xicString=($quantifType eq 'XIC')? "<TR><TD colspan=10><B>Criteria for XIC extraction:</B><BR>Mass-Window=<INPUT name=\"windowMass\" value=\"0.1\" size=\"5\"> Da.<BR>Time-window=<INPUT name=\"windowTime\" value=\"90\" size=\"2\"> seconds.<BR>Peptide Rank=<INPUT name=\"pepRank\" value=\"2\" size=\"1\">.</TD></TR>":'';
	#print "<TR bgcolor=$bgColor><TD colspan=10><B>Criteria used for XIC extraction:</B> <INPUT name=\"windowMass\" value=\"0.1\" size=\"5\"> Da - Time-window=<INPUT name=\"windowTime\" value=\"90\" size=\"5\"> seconds. - Peptide Rank=<INPUT name=\"pepRank\" value=\"2\" size=\"3\"></TD>";
	print qq
|<TR><TD colspan=10><INPUT type="submit" name="launch" value="Launch Quantification" class="title3"$disabSubmit>&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="cancelAction();"></TD></TR>
</TABLE>
</FORM>
</CENTER>
<DIV id="divDescription" class="clDescriptionCont"></DIV>
<SCRIPT type="text/javascript">setPopup();</SCRIPT>
</BODY>
<HTML>
|;
}
else {# DESIGN Method

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
<TH class="bBorder" colspan=2>&nbsp;Peptide quantifications&nbsp;<BR>&nbsp;Extend selection to:<SELECT name="propagQuantScope" id="propagQuantScope:$expConditionID" onchange="synchronizePropagation(this)"><OPTION value="all">All states</OPTION><OPTION value="state">This state</OPTION><OPTION value="none">None</OPTION></SELECT>&nbsp;</TH>
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
				print "<TR bgcolor=$bgColor class=\"list\"><TH align=\"left\" nowrap>&nbsp;($designInfo{OBSERVATION}{$obsCode}{BIOSAMPLE})&nbsp;$designInfo{ANALYSIS}{$anaID}{HIERARCHY}&nbsp;>&nbsp;$designInfo{OBSERVATION}{$obsCode}{NAME}&nbsp;</TH><TH align=\"left\" nowrap>&nbsp;$fracGroupStrg&nbsp;$techRepStrg&nbsp;</TH>";
				print qq
|<TD><INPUT type="checkbox" name="condObs:$expConditionID" value="$obsID" onclick="updateObsCheckStatus('condObs:$expConditionID',this)" checked ></TD>
|;
				if ($designInfo{'ANALYSIS'}{$anaID}{'HAS_QUANTIF'}) {
					print qq
|<TD><SELECT name="anaXICQuanti:$obsID" id="anaXICQuanti:$obsID" data-state="$expConditionID" onchange="propagateQuantifSelection(this)"><OPTION value="">-= Select =-</OPTION>
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
								print "<OPTION name=\"$labelType\" id=\"opt:$obsID:$quantiID\" value=\"$obsID:$quantiID:$anaID:$trueTargetPos:$fracGroup:$techRepGroup\" onmouseover=\"popup('$designInfo{QUANTIFICATION}{$quantiID}{ANNOT_STRING}')\" onmouseout=\"popout()\">$designInfo{QUANTIFICATION}{$quantiID}{NAME}</OPTION>\n"; #$selStrg
								#$firstOpt=0;
							}
						}
						print "</OPTGROUP>\n";
					}
					print "</SELECT></TD>\n";
				}
				else {print "<TH align=\"left\" bgcolor=$bgColor colspan=2>&nbsp;* None found  (SSPA Only!) *&nbsp;<INPUT type=\"hidden\" name=\"anaXICQuanti:$obsID\" value=\"$obsID:0:$anaID:0:$fracGroup:$techRepGroup\"></TH>";} # no pepQuantif, no trueTargetPos
				print "</TR>\n";
			}
		}
		print "</TABLE>\n</DIV>\n";
	}

#	$dbh=&promsConfig::dbConnect;
#	my $bgColor=$lightColor;
#	foreach my $quantiIDAnaID (sort{ $anaXIC{$a}{'QUANTI_NAME'} cmp $anaXIC{$b}{'QUANTI_NAME'} || $anaXIC{$a}{'ANA_NAME'} cmp $anaXIC{$b}{'ANA_NAME'} } keys %anaXIC) {# Order by analysis name
#		my ($statusString,$selected)=($anaXIC{$quantiIDAnaID}{'STATUS'}==1)? ('Finished',''):('Not finished','disabled');
#		my $xicParamUsedString='';
#		if ($anaXIC{$quantiIDAnaID}{'ANNOT'} =~ /^elutionWindow/){ # Former XIC (in-house) extraction -> Not MCQ
#			my ($eluWindowTag,$massWindowTag,$maxRankTag)=split(/;/,$anaXIC{$quantiIDAnaID}{'ANNOT'});
#			my ($eluWindowStrg,$eluWindowVal)=split(/=/,$eluWindowTag);
#			my ($massWindowStrg,$massWindowVal)=split(/=/,$massWindowTag);
#			my ($maxRankStrg,$maxRankVal)=split(/=/,$maxRankTag);
#			$xicParamUsedString=qq
#|<B>Mass-window(Da):</B> $massWindowVal<BR>
#<B>Time-window(sec.):</B> $eluWindowVal<BR>
#<B>Max-rank used:</B> $maxRankVal
#|;
#		}
#		else{ # MassChroQ parameters annotation
#			my %quantifParameters;
#			foreach my $parameter (split(/::/,$anaXIC{$quantiIDAnaID}{'ANNOT'}) ){
#				next unless $parameter;
#				my ($parameterName,$parameterValue)=split(/=/,$parameter);
#				$quantifParameters{$parameterName}=$parameterValue;
#			}
#			if ($quantifParameters{'REFERENCE'} ){
#				my ($refName)=$dbh->selectrow_array("SELECT NAME FROM ANALYSIS WHERE ID_ANALYSIS=$quantifParameters{'REFERENCE'}");
#				my $extraAlgoString=($quantifParameters{'EXTRACTION_ALGO'} eq 'OBI')? 'OBI-Warp' : 'ms2';
#				$xicParamUsedString=qq
#|<B>Type of XIC:</B> $quantifParameters{XIC_EXTRACTION_TYPE}<BR>
#<B>Reference:</B> $refName<BR>
#<B>Alignment algo.:</B> $extraAlgoString<BR>
#<B>Align from</B> $quantifParameters{MZ_ALIGN_RANGE_MIN} <B>to</B> $quantifParameters{MZ_ALIGN_RANGE_MAX}
#|;
#			}
#			else{$xicParamUsedString='-';}
#		}
#		print qq
#|<TR valign=middle bgcolor=$bgColor>
#<TD class="center"><INPUT type="checkbox" name="anaList" value="$quantiIDAnaID" $selected></TD>
#<TH>&nbsp;$anaXIC{$quantiIDAnaID}{'ANA_NAME'}&nbsp;</TH>
#<TD class="center"><SELECT id="Cond:$quantiIDAnaID" name="Cond:$quantiIDAnaID" onchange="">$selCondStrg</TD>
#<TD class="center">&nbsp;$anaXIC{$quantiIDAnaID}{'QUANTI_NAME'}&nbsp;</TD>
#<TD class="center">&nbsp;$anaXIC{$quantiIDAnaID}{'MZXML'}&nbsp;</TD>
#<TD nowrap>$xicParamUsedString</TD>
#<TD class="center">&nbsp;$statusString&nbsp;</TD>
#</TR>
#|;
#		$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
#	}
#	$dbh->disconnect;
#	print qq
#|<TR><TD colspan=5><INPUT type="submit" name="launch" value="Launch Quantification" class="title3">&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="cancelAction();"></TD></TR>
#</TABLE>
#</FORM>
#</CENTER>
	print qq
|</FORM>
</CENTER>
<BR><BR>
<DIV id="divDescription" class="clDescriptionCont"></DIV>
<SCRIPT type="text/javascript">setPopup();</SCRIPT>
</BODY>
<HTML>
|;
}


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
	my $numJobs=1; # default

	my @selectedConds; # for Design only
	my $usedCondNum=0;
	if ($quantifType eq 'DESIGN') {

		my $dbh=&promsConfig::dbConnect;
		($experimentID)=$dbh->selectrow_array("SELECT ID_EXPERIMENT FROM DESIGN WHERE ID_DESIGN=$itemID");
		$projectID=&promsMod::getProjectID($dbh,$experimentID,'experiment');
		$dbh->disconnect;

		foreach my $nbCond (1..param('nbMaxCond')){ # Information to compute the ratios
			my $condID=param("state_$nbCond");
			push @selectedConds,$condID if ($condID && $condID ne 'Select');
		}
		$usedCondNum=scalar @selectedConds; # for Design only

		###<Multi-quantif>###
		if (param('multiQuantif')) { # only for 'All vs State1'
			$numJobs=$usedCondNum-1;
		}
	}
	elsif ($quantifType=~/SIN|EMPAI/) { # not for XICMCQ
		$numJobs=scalar @quantItemList;
	}

	mkdir "$promsPath{tmp}/quantification" unless -e "$promsPath{tmp}/quantification";
	my $currentQuantifDir="$promsPath{tmp}/quantification/current";
	mkdir $currentQuantifDir unless -e $currentQuantifDir;

	my $quantifBaseName=param('quantifName');
	my $refConfName; # for multi-Q only

	my $quantifDate=strftime("%Y%m%d%H%M%S",localtime);
	while (-e "$promsPath{tmp}/quantification/$quantifDate") { # to prevent multi-user launch collision
		sleep 2;
		$quantifDate=strftime("%Y%m%d%H%M%S",localtime);
	}

	print '<BR><FONT class="title3">Launching quantification process';
	if ($numJobs > 1) {print "es [Master job=#$quantifDate]:\n";}
	else {print "[job #$quantifDate]...";}

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
		if ($numJobs==1) {@states=@selectedConds;} # normal case
		elsif ($quantifType eq 'DESIGN') { # multi-quantif "All vs State1" splitted into 1 quantif/ratio
			@states=($selectedConds[0],$selectedConds[$jobPos]);
			if ($quantifName=~/%(TEST|REF[^%]*)%/i) {
				my $dbh=&promsConfig::dbConnect;
				my $sthCond=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
				if ($quantifName=~/%REF[^%]*%/i) {
					unless ($refConfName) { # do only once
						$sthCond->execute($states[0]);
						($refConfName)=$sthCond->fetchrow_array;
					}
					$quantifName=~s/%REF[^%]*%/$refConfName/gi;
				}
				if ($quantifName=~/%TEST%/i) {
					$sthCond->execute($states[1]);
					my ($testCondName)=$sthCond->fetchrow_array;
					$quantifName=~s/%TEST%/$testCondName/gi;
				}
				$sthCond->finish;
				$dbh->disconnect;
				$quantifName=~s/%#%/$jobPos/g;
			}
		}
		my $usedStateNum=scalar @states;

		###<jobs info & data >###
		open (INFO,">$quantifDir/quantif_info.txt"); # item valueR valueDB
		print INFO "USER=$userID\n";
		if ($quantifType eq 'XICMCQ'){print INFO "TYPE=XIC\n";} # MCQ -> Alternative to XIC extraction
		else {print INFO "TYPE=$quantifType\n";}
		my $algoType;
		if ($quantifType eq 'DESIGN') {
			my $labeling=param('labeling');
			##($algoType,my $refType)=split(':',param('algoType')); # $refType defined only for SimpleRatio & MSstats (multiRef or singleRef)
			##my $ratioType=($algoType=~/S\w+Ratio/)? $algoType : ($algoType eq 'MSstats')? 'SimpleRatio' : ($algoType eq 'SSPA')? 'SSPA' : 'Ratio'; # ratioType is Ratio also for TnPQ!!!
			($algoType,my $ratioType,my $refType)=split(':',param('algoType')); # $ratioType & $refType undef for SSPA
#$algoType=$ratioType='SuperRatio';
#$refType='singleRef';
			#my $numChannels=param('numChannels'); # not for DESIGN
			#my $maxReplicates=param('maxReplicates'); # not for DESIGN
			my $pepCharge=param('pepCharge') || 'all';
			my $pepSrc=($pepCharge eq 'best')? 'best' : (param('pepSource'))? param('pepSource') : 'all'; # only for SILAC
			my $minInfRatios=param('minInfRatios') || 0;
			print INFO "PARAMETERS:\n";
			print INFO "LABEL\t\t$labeling\n"; # only for DB not R (same as $quantifType for internal quantif)
			print INFO "QUANTIF_NAME\t\t$quantifName\n"; # only for DB
			print INFO "TOP_N\t\t",param('topN'),"\n" if param('topN'); # only for DB not R (topN is undef if PTM quantif)
			print INFO "MATCH_PEP\t\t",param('matchingPep'),"\n" if param('matchingPep'); # only for DB not R
			my $ptmScope=param('ptmScope');
			if ($ptmScope) { # only for DB: PTM quantification
				print INFO "PTM_SCOPE\t\t$ptmScope\n";
				print INFO "PTM_POS\t\t";
				my $phosphoID=param('phosphoID');
				if ($ptmScope==$phosphoID && $algoType ne 'MSstats') { # phospho quantif
					print INFO 'PRS:',param('okPRS'),"\t",param('badPRS'),"\n"; # No PRS for DIA
				}
				else { # other PTM quantif
					my $ambigPos=(param('ambiguousPos'))? 'ambiguous' : 'valid';
					print INFO "$ambigPos\n";
				}
				#my $ambigPos=param('ambiguousPos') || 0;
				#print INFO "AMBIGUOUS_QPTM_POS\t\t$ambigPos\n";
				#print INFO "REF_QUANTIF\t\t",param('refQuantifID'),"\n" if param('refQuantifID');
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
				print INFO "PEPTIDES\t\t",param('pepSpecifity'),"\t",param('pepCleavage'),"\t$ptmFilter\t",param('pepFocus'),"\n";
			}
			else {
				print INFO "PEPTIDES\t\t",param('pepSpecifity'),"\t",param('pepCleavage'),"\t$ptmFilter\t$pepCharge\t$pepSrc","\n"; # only for DB ***Leave pepCharge & pepSrc for SWATH for compatibility with showProtQuantif***
			}
			#<Protein selection (since v1.7.2)
			if (param('protSelection')) {
				print INFO "PROTEINS\t\t",param('protSelType'),"\t#",param('protSelection'),"\n";
			}
			print INFO "ALGO_TYPE\t\t$algoType\n"; # if $ratioType ne 'Ratio'; # needed for DB (not needed if Ratio) / Overwritten by launchQuantif if PEP_xxx incompatibilty
			if ($algoType ne 'SSPA') {
				print INFO "RATIO_TYPE\t\t$ratioType\n";
				if ($ratioType=~/S\w+Ratio/) { # S*Ratio & DIA
					print INFO "SINGLE_REF\t\t1\n" if ($refType && $refType eq 'singleRef');
				}
				my $normMethod=param('biasCorrect');
				my $biasCorrect=($normMethod eq 'none.none')? 'FALSE' : 'TRUE';
				print INFO "BIAS_CORRECTION\t\t$biasCorrect"; # for DB only
				print INFO "\t#",param('refProt'),"\t",param('refProtSelType') if (param('refProt') && ($normMethod=~/REF_PROT|globalStandards/ || ($algoType=~/PEP_/ && $normMethod !~/quantile/i))); # all algos except SSPA
				print INFO "\n";
				if ($ratioType=~/S\w+Ratio/) { # S*Ratio & MSstats
					my ($usedMethod)=($normMethod eq 'none.none' && $algoType eq 'MSstats')? 'FALSE' : $normMethod;
					print INFO "NORMALIZATION_METHOD\tnormalization.method\t$usedMethod\n";
				}
			}
			if ($ratioType=~/S\w+Ratio/) { # since SuperRatio
				#>SWATH with MSstats (list of ratios to be computed: C2/C1 C3/C1 ... Cn/C1 C3/C2 ... Cn/C2 ...)
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
			print INFO "ID_DESIGN\t\t$itemID\n";
			print INFO "STATES\t\t",join("\t",@states),"\n";

			#foreach my $quantiIDAnaIDCondID (@quantItemList){# One XIC-quanti could contain multi ANA quantified
			#	my @condInfo=split(/:/,$quantiIDAnaIDCondID);
			#	print INFO "QUANTITOCOND\t\t".param("anaXICQuanti:$quantiIDAnaIDCondID")."\n" if $states{$condInfo[1]};
			#}

			# QUANTITOCOND  (obsID quantifID fracGroup)
			#print INFO "ANALYSES/OBSERVATIONS & PARENT QUANTIFICATIONS:\n";

			my $residualVar=param('residualVar'); # defined only for $algoType=~/S\w+Ratio/
			my $checkBioRep=($algoType=~/PEP_/ && $residualVar eq 'auto')? 1 : 0;
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
			if ($algoType=~/PEP_/) { # Not Ratio,TnPQ,SSPA,MSstats
				#my $Rdesign=($labeling eq 'FREE')? 'LABELFREE' : ($ratioType eq 'SuperRatio')? 'SUPERRATIO' : 'LABELED'; # ***New algo (17/08/17)***
				my $Rdesign=($algoType eq 'PEP_RATIO' && $okSuperRatio)? 'PEP_RATIO' : 'PEP_INTENSITY'; # will overwrite ALGO_TYPE above if necessary (in launchQuantifications.cgi)
				print INFO "\tdesign\t$Rdesign\n"; # only in R

				if ($checkBioRep) {
					$residualVar=($numBioRepOK==$usedStateNum)? 'biological' : 'technical';
				}
				print INFO "RESIDUAL_VAR\tresidual.variability\t$residualVar\n";
				print INFO "FDR_CONTROL\tpAdj.method\t",param('fdrControl'),"\n";
			}
			print INFO "MIN_INF_RATIOS\t\t$minInfRatios\n" if $algoType !~ /MSstats|SSPA/; # if $usedStateNum==2; # +/- infinite ratio switch (DB only)

		}

		elsif ($quantifType eq 'XICMCQ') { # MassChroQ extraction
			my $raw=param('rawdata_extraction'); # profile or centroid
			my $algo=param('alignment_method'); # OBI-Warp or in-house algo
			my $anaReference=param('refAna_alignment'); # Reference analysis
			my $mzAnaRangeMin=param('mz_start'); # range for alignment -OBIwarp
			my $mzAnaRangeMax=param('mz_stop'); # range for alignment -OBIwarp
			my $mzRangeMin=param('mzRange_start'); # range for allChargeState computation
			my $mzRangeMax=param('mzRange_stop'); # range for allChargeState computation
			my $tendency=param('ms2_tendency'); # ms2 parameter
			my $smouthMS2=param('ms2_smoothing'); # ms2 parameter
			my $smouthMS1=param('ms1_smoothing'); # ms2 parameter
			my $xicType=param('XIC_type'); # sum or max
			my $xicRange=param('XIC_range');
			my $mzTolMin=param('mztol_min');
			my $mzTolMax=param('mztol_max');
			my $xicValType=param('XIC_rt'); # real_or_mean, mean or post_matching
			my $detectionThresholdMin=param('dt_start');
			my $detectionThresholdMax=param('dt_stop');
			my $antiSpike=param('anti_spike'); # anti-spike value
			my $medMin=param('med_min'); # half-median min
			my $medMax=param('med_max'); # half-median max
			my $smoothVal=param('smooth_val'); # smoothing value
			my $allChargeStates=param('allChargeStates')?param('allChargeStates') : 0;
			my $extractTraces=param('traces')?param('traces'):0;

			print INFO "PARAMETERS:\n";
			print INFO "QUANTIF_NAME\t\t",param('quantifName'),"\n";
			print INFO "RAWDATA_ACQUISITION\t\t$raw\n";
			print INFO "EXTRACTION_ALGO\t\t$algo\n";
			print INFO "REFERENCE\t\t$anaReference\n";
			if ($algo eq 'OBI') {
				print INFO "MZ_ALIGN_RANGE_MIN\t\t$mzAnaRangeMin\n";
				print INFO "MZ_ALIGN_RANGE_MAX\t\t$mzAnaRangeMax\n";
			}
			else{
				print INFO "MS2_TENDENCY\t\t$tendency\n";
				if (param('ms2Smooth')){
					print INFO "MS2_SMOUTHING\t\t$smouthMS2\n";
				}
				if (param('ms1Smooth')){
					print INFO "MS1_SMOUTHING\t\t$smouthMS1\n";
				}
			}

			print INFO "XIC_EXTRACTION_TYPE\t\t$xicType\n";
			print INFO "XIC_RANGE\t\t$xicRange\n";
			print INFO "MZTOL_MIN\t\t$mzTolMin\n";
			print INFO "MZTOL_MAX\t\t$mzTolMax\n";
			print INFO "XIC_VAL\t\t$xicValType\n";
			print INFO "DT_START\t\t$detectionThresholdMin\n";
			print INFO "DT_STOP\t\t$detectionThresholdMax\n";
			### Only one filter applies... Modified on 05/09/2014
			if (param('asFilter')) {
				print INFO "ANTISPIKE\t\t$antiSpike\n";
			}
			if (param('smFilter')){
				print INFO "SMOOTH\t\t$smoothVal\n";
			}
			if (param('hmFilter')){
				print INFO "MED_MIN\t\t$medMin\n";
				print INFO "MED_MAX\t\t$medMax\n";
			}
			print INFO "ALLCHARGESTATES\t\t$allChargeStates\n";
			if( $allChargeStates ) {
				print INFO "MZ_RANGE_MIN\t\t$mzRangeMin\n";
				print INFO "MZ_RANGE_MAX\t\t$mzRangeMax\n";
			}
			print INFO "TRACES\t\t$extractTraces\n";

			if (param('extractL') eq 'SILAC') {
				###> Get the labeled modification so as to write QUANTIF_ANNOT according to SILAC
				my %anaLabelMods;
				my $dbh=&promsConfig::dbConnect;
				my $sthMods=$dbh->prepare("SELECT ID_MODIFICATION,PSI_MS_NAME,MONO_MASS FROM MODIFICATION WHERE VALID_STATUS=1 AND IS_LABEL=1");
				$sthMods->execute;
				while (my ($modID,$psiMsName,$monoMass)=$sthMods->fetchrow_array) {
					@{$anaLabelMods{$modID}}=($psiMsName,$monoMass);
				}
				$sthMods->finish;
				$dbh->disconnect;

				my @channelsToExtract;
				foreach my $channel (1..3) {
					my @labelsToExtract;
					foreach my $label (1..$maxLabPerChannel) {
						if (param("lab$label:name_channel$channel") && param("lab$label:name_channel$channel") ne "") {
							# exemple: Lys8#Label:13C(6)15N(2)#K#Lys+8 -> $labelModifAlias#$labelModifName#$modifRes#$searchModifName
							my $modID=param("lab$label:channel$channel");
							my ($labName,$psiMsName,$deltaMass) = ($modID > 0 )? (param("lab$label:name_channel$channel"),@{$anaLabelMods{$modID}}) : ('None','No label','0.00000');
							my $modifRes=(param("lab$label:tg_channel$channel") ne 'ANY') ? param("lab$label:tg_channel$channel") : param("lab$label:sp_channel$channel");
							push @labelsToExtract, "$labName#$psiMsName#$modifRes##$modID"; #
						}
					}

					if (param("name_channel$channel") && param("name_channel$channel") ne "") {
						my $labelString=join("@",@labelsToExtract);
						push @channelsToExtract, "$channel;".param("name_channel$channel").";$labelString" if $labelString;
					}
				}
				my $channelString=join("::",@channelsToExtract);
				print INFO "CHANNELS\t\t$channelString\n";
			}

			foreach my $anaID (@quantItemList){
				print INFO "MZXML\t\t$anaID:".param("file_$anaID")."\n";
			}
			###print INFO "ANALYSES:\n";
		}
		else { # SIN,EMPAI,...?
			print INFO "PARAMETERS:\n";
			my $qRootName=($quantifType eq 'SIN')? 'PEP_SI_' : $quantifType;
			print INFO "QUANTIF_NAME\t\t$qRootName\n";
			print INFO "ANALYSES:\n";
		}

		if ($quantifType =~ /SIN|EMPAI/) { # =SIN|EMPAI (no more SILAC|ITRAQ|TMT). For DESIGN-based wait flag created by launchquantification, no analysis are done
			#foreach my $quantItemID (@quantItemList) { # anaID or anaID.parentQuantifID.ratioNum for SILAC & iTRAQ
				my $anaID=$quantItemList[$jobPos-1];
				#$quantItemID.='.'.param("quantif_$quantItemID") if $quantifType=~/SILAC|ITRAQ|TMT/;
				print INFO "$anaID\n";
				open(FLAG,">$currentQuantifDir/$anaID\_$jobDir\_wait.flag"); # flag file used by watchQuantif
				print FLAG "#";
				close FLAG;
			#}
		}
		elsif ($quantifType =~ /DESIGN/) {
			open(FLAG,">$currentQuantifDir/$jobDir\_request.flag"); # flag file used by watchQuantif
			print FLAG "#";
			close FLAG;
		}

		close INFO;

		my $fullQuantifType=($quantifType eq 'DESIGN' && $algoType=~/MSstats|SSPA/)? $quantifType.':'.$algoType : $quantifType;
		my $quantItemStrg=join(':',@quantItemList);
		$quantItemStrg='' unless $quantItemStrg;
		if ($numJobs==1) {
			###<Forking to launch quantifications>###
			my $childPid = fork;
			unless ($childPid) { # child here
				#>Disconnecting from server
				open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
				open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
				#open STDERR, ">>$promsPath{logs}/launchQuantification.log";
				system "./launchQuantifications.pl single $ENV{REMOTE_USER} $jobDir $fullQuantifType $quantItemStrg";
	#system "./launchQuantifications.pl single $ENV{REMOTE_USER} $jobDir $fullQuantifType $quantItemStrg 2> $promsPath{tmp}/quantification/$jobDir\_errorQuantif.txt";
				exit;
			}
		}
		else {
			push @jobList,"$jobDir#$fullQuantifType#$quantItemStrg";
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
	print qq
|<SCRIPT type="text/javascript">
var watchQuantifWin=window.open("$promsPath{cgi}/watchQuantifications.cgi",'WatchQuantifWindow','width=1000,height=500');
watchQuantifWin.focus();
//parent.optionFrame.watchQuantifications(); // starts watchQuantificationw.cgi
</SCRIPT>
|;
#exit; # DEBUG!!!!!
	sleep 5;
	print qq
|<SCRIPT LANGUAGE="JavaScript">
//top.promsFrame.selectedAction='summary';
|;
	if ($quantifType eq 'DESIGN') {
		#print "$updateFrameString;\n";
		print "parent.itemFrame.location=\"$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=design:$itemID&ACT=experiment&EXPERIMENT=EXPERIMENT:$experimentID&ISNAVFRAME=0&VIEW=quanti\";\n";
	}
	else {
		print "parent.optionFrame.selectOption(parent.optionFrame.document.getElementById('summary'));\n"; #// refresh optionFrame with summary option
	}
	print qq
|</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

####> Revision history
# 2.0.0 Full support for multi-job launches (PP 17/09/18)
# 1.9.9 Added more contraints on peptide/PTM quantif/PhosphoRS options based on algo selection (PP 19/07/18)
# 1.9.8 Improved auto-propagation of peptide quantification selection (PP 01/06/18)
# 1.9.7 Form & template are reset by labeling change, skip template-based biological replicates selection & template's project is displayed (PP 09/05/18)
# 1.9.6 Add a class ('template') to all parameters that will be modified by the template selection (MLP 07/05/18)
# 1.9.5 Add a function to restore default parameters (use for templates selection) (MLP 26/04/18)
# 1.9.4 [Fix] minor bug in R parameter value for MSstats when no normalization (PP 23/04/18)
# 1.9.3 More modif to select a template (MLP 18/04/18)
# 1.9.2 Modif to select a template (MLP 16/04/18)
# 1.9.1 Minor change to allow proper frame refreshing (PP 06/04/18)
# 1.9.0 Multi-quantif launch support & auto-selection of all states (PP 28/03/18)
# 1.8.7 Declare NORMALIZATION_METHOD parameter even if no normalization (PP 13/03/18)
# 1.8.6 [Fix] javascript bugs for MassChroQ XML file selection check & new checks on labeling coherence (PP 13/02/18)
# 1.8.5 Minor modif (MLP 24/01/18)
# 1.8.4 Added 'clusters' parameter for MSstats (PP 24/01/18)
# 1.8.3 Minor modifications for XIC extraction with MassChroQ (GA 14/12/17)
# 1.8.2 Modification of checked options for XIC extraction with MassChroQ (GA 13/12/17)
# 1.8.1 Better detection of bio/tech/replicates from design (PP 04/12/17)
# 1.8.0 Change in parameters to match upgraded ratio quantification R scripts (PP 27/11/17)
# 1.7.8 New option for peptide selection for SSPA, fix bug to disable Super ratio for Label-free & commented substitution 3plex/2plex (PP 20/04/17)
# 1.7.7 Compatible with TMT labeling (PP 18/01/17)
# 1.7.6 Minor modif: change "required" (SELECT "refAna_alignment") by checkForm verification (MLP 09/01/17)
# 1.7.5 Add "required" for SELECT "refAna_alignment" (MLP 09/01/17)
# 1.7.4 Minor modification of masschroq parameters and color used for labelling (GA 10/10/16)
# 1.7.3 Minor modification (GA 24/08/16)
# 1.7.2 Added Protein selection/exclusion option, MSstats global standard & "Distinct peptide sequences" option for SSPA & (PP 19/08/16)
# 1.7.1 Uses states position (not indexes) in contrasts.matrix for SWATH (PP 04/08/16)
# 1.7.0 Adding Swath (MSstats) & SSPA (PP 02/08/16)
# 1.6.6 Minor modification (GA 12/05/16)
# 1.6.5 Removed all p-value correction methods except BH & hide topN selection on Label-Free modif quantification (PP 26/04/16)
# 1.6.4 Minor bug fix in normalization method selection for label-free (PP 24/03/16)
# 1.6.3 Minor display change to better distinguish Samples (PP 08/03/16)
# 1.6.2 Added matching peptides option for label-free quantification (PP 29/02/16)
# 1.6.1 Fixed uninitialized $nbCond (PP 16/02/16)
# 1.6.0 Handles complex design,"Simple ratios"=multiSILAC algorithm for SILAC triplex (PP 29/01/16)
# 1.5.2 Modification of XIC extraction form (GA 05/09/14)
# 1.5.1 New Observation management: OBS_EXPCONDITION table & check on ref/test states compatibility for labeled quantif (PP 24/07/14)
# 1.5.0 Super Ratio option & Bug fix in condition detection is some cases (PP 21/05/14)
# 1.4.1 Modification of XIC extraction form (GA 09/04/14)
# 1.4.0 Add xic-traces option for masschroq extractions (GA 31/03/14)
# 1.3.9 Modification of MassChroQ parameters (GA 28/03/14)
# 1.3.8 Minor modification of parameters for MassChroQ extraction (GA 25/03/14)
# 1.3.7 display always header peptide quantification (SL 25/03/14)
# 1.3.6 Syntax bug fix (PP 19/03/14)
# 1.3.5 Uses mkdir instead of make_path (PP 10/03/14)
# 1.3.4 Change syntax $x=$y // $z to old ()? one (PP 10/03/14)
# 1.3.3 Peptide quantif selection & improved state selection for internal quantif (PP 26/03/14)
# 1.3.2 Various Fixes & updates <BR> Fix SILAC No Label channel annotation for MassChroQ XIC (PP 20/02/14)
# 1.3.1 Form update for label XIC extraction (GA 11/12/13)
# 1.3.0 system command removal (PP 08/11/13)
# 1.2.9 Bug fix in infinite ratios switch for DESIGN (PP 23/10/13)
# 1.2.8 Added infinite ratios switch parameter (PP 21/10/13)
# 1.2.7 iTRAQ -> ITRAQ for QUANTIF_ANNOT field (PP 07/10/13)
# 1.2.6 Minor display bug fix (PP 30/09/13)
# 1.2.5 MCQ -> XICMCQ & launch prepared like DESIGN (PP 09/09/13)
# 1.2.4 'TNPQ'->'DESIGN' & multi-labeling management (PP 29/08/13)
# 1.2.3 Minor change (PP 12/07/13)
# 1.2.2 Added pool management (PP 12/07/13)
# 1.2.1 Fix bugs popup MS2 & condition display (PP 05/07/13)
# 1.2.0 Add MS2 parameters for MassChroQ alingment (GA 04/07/13)
# 1.1.9 Better label detection & extended design ... (PP 04/07/13)
# 1.1.8 Modification for QUANTITATION launch based on new ANA_EXPCONDITION table and log file addition launchQuantification.log (GA 08/03/13)
# 1.1.7 Added more peptide filters: missed cleavage and PTM/unmodified peptide pair (PP 11/01/13)
# 1.1.6 Multi-databank search, bug fix in SILAC/iTRAQ quantif_info.txt file & code cleaning (PP 08/01/13)
# 1.1.5 Minor modification : remove 3 replicates option for Pep-Ratio + add checkbox for MCQ extraction (GA 05/11/12)
# 1.1.4 Modification for XIC extraction (quantif-name) and XIC parameters + Modification in TNPQ parameters (GA 07/11/12)
# 1.1.3 Minor modification ACT=experiment to make it homogene with openProject.cgi (GA 18/10/12)
# 1.1.2 Add MassChroQ extraction option -> similar to XIC (GA 04/09/12)
# 1.1.1 Check on SILAC/iTRAQ peptide quantif status before allowing protein quantif (PP 07/06/12)
# 1.1.0 Minor changes & bug fixes (PP 30/04/12)
# 1.0.9 Minor bug correction -> remove printError call (GA 18/04/12)
# 1.0.8 minor bug correction for perl/javascript: 'eq' vs '==' (GA 04/04/12)
# 1.0.7 merge 1.0.6GA & 1.0.5PP (PP 02/04/12)
# 1.0.6 Updates for label-free based quantification like TNPQ, XIC, etc (GA 16/01/12)
# 1.0.5 Updates for label-based internal quantifications (PP 26/08/11)
# 1.0.4 More Quantif methods (GA ...)
# 1.0.3 Fix file detection bug for T3PQ
