#!/usr/local/bin/perl -w

################################################################################
# manageDesignCondition.cgi    1.3.4                                           #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Create, edit & store Design, Experimental States & associate observations    #
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
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use strict;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};

#print header(-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
#############################
####>Fetching parameters<####
#############################
my $action=param('ACT');
my $rowID=param('ID'); # if action=add => rowID=parentID
my $item=uc(param('ITEM')); # item being processed (child if 'add')
my $parentItem=param('PARENT'); # only for add
my $parentID=param('PARENT_ID'); #-- action=store|update
my $projectID=param('PROJECT_ID');
my $restrictBranch=param('RESTRICT_BRANCH') || ''; # undef sometimes
my ($trueItem,$trueItemID)=($item,$rowID);
my $designID;
my $MAX_NUM_STATES=100; # Maximal number of states that can be drawn on a frame.
#my $numExpCond; # needed when deleting/adding an ExpState for proper optionFrame refreshing

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;
#my ($userStatus)=$dbh->selectrow_array("SELECT USER_STATUS FROM USER_LIST WHERE ID_USER='$userID'");
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=$userInfo[2]->{$projectID};
if ($action eq 'add') {
	#$projectID=&promsMod::getProjectID($dbh,$rowID,$parentItem) unless $projectID;
	$parentID=$rowID;
}

############################
####>Fetching item info<####
############################
my $itemID;
my $itemType=&promsMod::getItemType($item);
#my @itemInfo;
my %itemIcones=&promsConfig::getItemIcones;
my (@nameText,@valueText); # Globals

##############################
####>Applying the effects<####
##############################
my ($name,$description);
my $newItemID;
my $expStateStrg='';
my $startNewStates=0;
my $duplicateWarningStrg='';

if ($action =~ /edit/ || $action eq 'summary'){

	$itemID=$rowID;

	@nameText=('Name','Description');
	$itemID=param('ITEM_ID') unless $itemID;
	if ($action =~ /edit/ && $item eq 'EXPCONDITION') { # WARNING:  $itemID & $item are redefined!!!
		($itemID)=$dbh->selectrow_array("SELECT ID_DESIGN FROM EXPCONDITION WHERE ID_$item=$itemID");
		$item='DESIGN';
	}

	my ($designID)=($item eq 'DESIGN')? $itemID : ($item eq 'EXPCONDITION')? $dbh->selectrow_array("SELECT ID_DESIGN FROM EXPCONDITION WHERE ID_EXPCONDITION=$itemID") : -1 ;

	($name,$description,my $experimentID)=$dbh->selectrow_array("SELECT NAME,DES,ID_EXPERIMENT FROM DESIGN WHERE ID_DESIGN=$designID");
	$restrictBranch='experiment:'.$experimentID unless $restrictBranch;

	##>Fetching list Experiment children to restrict data set selection [Commented ---> Not finished]
	my $restrictStrg=qq
|Restrict display to data in <SELECT name="RESTRICT_BRANCH" class="title3" style="color:#DD0000" onchange="updateDisplayRestriction(this.value)">
	<OPTION value="experiment:$experimentID">Whole Experiment</OPTION>
|;
	my $sthExpChild=$dbh->prepare("SELECT ID_SAMPLE,NAME FROM SAMPLE WHERE ID_EXPERIMENT=$experimentID ORDER BY DISPLAY_POS");
	$sthExpChild->execute;
	while (my ($sampID,$sampName)=$sthExpChild->fetchrow_array) {
		my $branchID="sample:$sampID";
		$restrictStrg.="<OPTION value=\"$branchID\"";
		$restrictStrg.=' selected' if $branchID eq $restrictBranch;
		$restrictStrg.=">$sampName</OPTION>\n";
	}
	$restrictStrg.="</SELECT>\n";
	$sthExpChild->finish;

	my $sthExpCond=$dbh->prepare("SELECT ID_EXPCONDITION,NAME,DES FROM EXPCONDITION WHERE ID_DESIGN=$designID");
	$sthExpCond->execute;
	my $sthgetQ1=$dbh->prepare("SELECT COUNT(*) FROM EXPCONDITION_QUANTIF WHERE ID_EXPCONDITION=?");
	my ($light,$dark)=&promsConfig::getRowColors;
	$expStateStrg="\n<BR><BR><FONT class=\"title2\"><FONT color=#DD0000>States</FONT> associated with the design<BR></FONT><FONT class=\"title3\">$restrictStrg</FONT><BR><BR>";
	#$expStateStrg="\n<BR><BR><FONT class=\"title2\"><FONT color=#DD0000>States</FONT> associated with the design</FONT><BR><BR>";

	my (@sthGetAnaInfo,$maxFracGroup,$maxTechRepGroup);
	my $observationStrg;
	if ($restrictBranch =~ /experiment:/) {
		@sthGetAnaInfo=( # 2d-gel or sample
			$dbh->prepare("SELECT ID_ANALYSIS,'gel2d',G.NAME,'spot',SP.NAME,'analysis',A.NAME FROM ANALYSIS A,SAMPLE S,SPOT SP,GEL2D G,EXPERIMENT E,DESIGN D WHERE D.ID_EXPERIMENT=E.ID_EXPERIMENT AND E.ID_EXPERIMENT=G.ID_EXPERIMENT AND G.ID_GEL2D=SP.ID_GEL2D AND SP.ID_SPOT=S.ID_SPOT AND S.ID_SAMPLE=A.ID_SAMPLE AND VALID_STATUS > 0 AND ID_DESIGN=$designID ORDER BY G.DISPLAY_POS,SP.NAME,S.DISPLAY_POS,A.DISPLAY_POS"), # VALID_STATUS,
			$dbh->prepare("SELECT ID_ANALYSIS,'sample',S.NAME,'analysis',A.NAME FROM ANALYSIS A,SAMPLE S,EXPERIMENT E,DESIGN D WHERE D.ID_EXPERIMENT=E.ID_EXPERIMENT AND E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SAMPLE=A.ID_SAMPLE AND ID_SPOT IS NULL AND VALID_STATUS > 0 AND ID_DESIGN=$designID ORDER BY S.DISPLAY_POS,A.DISPLAY_POS") # VALID_STATUS,
		);
		$observationStrg='Observations :';
	}
	elsif ($restrictBranch =~ /sample:(\d+)/) {
		my $sampID=$1;
		@sthGetAnaInfo=(
			$dbh->prepare("SELECT ID_ANALYSIS,'sample',S.NAME,'analysis',A.NAME FROM ANALYSIS A,SAMPLE S WHERE S.ID_SAMPLE=$sampID AND A.ID_SAMPLE=$sampID AND ID_SPOT IS NULL AND VALID_STATUS > 0 ORDER BY A.DISPLAY_POS")
		);
		$observationStrg='<FONT color="#DD0000">Observations :</FONT>';
		if ($action eq 'editAE') {
			($maxFracGroup)=$dbh->selectrow_array("SELECT MAX(FRACTION_GROUP) FROM OBS_EXPCONDITION WHERE ID_EXPCONDITION=$trueItemID");
			($maxTechRepGroup)=$dbh->selectrow_array("SELECT MAX(TECH_REP_GROUP) FROM OBS_EXPCONDITION WHERE ID_EXPCONDITION=$trueItemID");
		}
	}
	# elsif Gel restriction?
	$maxFracGroup=0 unless $maxFracGroup;
	$maxTechRepGroup=0 unless $maxTechRepGroup;

	my $sthLA=$dbh->prepare("SELECT DISTINCT M.ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES FROM MODIFICATION M,ANALYSIS_MODIFICATION AM WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND IS_LABEL=1 AND ID_ANALYSIS=? ORDER BY M.ID_MODIFICATION");
	my $sthLQ=$dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND FOCUS='peptide' AND ID_ANALYSIS=? LIMIT 0,1");
	#my $sthObs=$dbh->prepare("SELECT O.ID_EXPCONDITION,TARGET_POS,FRACTION_GROUP,GROUP_CONCAT(ID_MODIFICATION ORDER BY ID_MODIFICATION SEPARATOR ',') FROM OBSERVATION O LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION JOIN EXPCONDITION E ON O.ID_EXPCONDITION=E.ID_EXPCONDITION WHERE ID_DESIGN=$designID AND ID_ANALYSIS=? GROUP BY O.ID_OBSERVATION");
	my $sthObs=$dbh->prepare("SELECT OE.ID_EXPCONDITION,O.TARGET_POS,OE.ID_OBSERVATION,OE.FRACTION_GROUP,OE.TECH_REP_GROUP,GROUP_CONCAT(ID_MODIFICATION ORDER BY ID_MODIFICATION SEPARATOR ',')
								FROM OBSERVATION O
								LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION
								JOIN OBS_EXPCONDITION OE ON O.ID_OBSERVATION=OE.ID_OBSERVATION
								JOIN EXPCONDITION E ON OE.ID_EXPCONDITION=E.ID_EXPCONDITION
								WHERE E.ID_DESIGN=$designID AND O.ID_ANALYSIS=? GROUP BY O.ID_OBSERVATION"
							);

	###> Query to get the quantitations associated to the design to prevent to break an association between an analysis and an expstate linked by a quantification.
	my $sthGetQ2=$dbh->prepare("SELECT COUNT(ID_QUANTIFICATION) FROM EXPCONDITION_QUANTIF WHERE ID_EXPCONDITION=? AND ID_QUANTIFICATION IN
									(SELECT ANA_QUANTIFICATION.ID_QUANTIFICATION
									 FROM QUANTIFICATION,ANA_QUANTIFICATION
									 WHERE QUANTIFICATION.ID_QUANTIFICATION=ANA_QUANTIFICATION.ID_QUANTIFICATION
									 AND ID_ANALYSIS=? AND ID_DESIGN=$designID
									)"
								);

	###> Get Analysis information in order to fill OBSERVATION, OBS_EXPCONDITION & OBS_MODIFICATION tables
	my (%hierarchy,%position,%expState,%stateChildren,%fractionGroup,%numFractionsInGroup,%numAnaCondQuantif,%technicalRepGroup,%numTechRepsInGroup,%usedTechRepGroups); # %validStatus
	my (%obs,%pool,%numPools,%techRepGr,%noQuantifAna); # for State summary
	my $i=0;
	foreach my $sthGAI (@sthGetAnaInfo) {
		$sthGAI->execute;
		while (my ($anaID,@anaHierarchy) = $sthGAI->fetchrow_array) { #$valStatus,
			#$validStatus{$anaID}=$valStatus;

			$sthLQ->execute($anaID);
			my ($quantifAnnot)=$sthLQ->fetchrow_array;

			##>A peptide quantification associated with Analysis
			# IMPORTANT NOTE (31/01/14 PP): TARGET_POS has meaning only for iTRAQ-like quantif
			# Otherwise: 0 for label-free and -1 for all channel pos of SILAC-like (true match given by OBS_MODIFCATION)
			if ($quantifAnnot) {
				$quantifAnnot=~s/::SOFTWARE=[^:]+//; # remove software & version info for back compatibility
				$quantifAnnot=~s/::CORRECTION=[^:]+//; # remove correction data if any
				my ($labelTypeStrg,@labelInfo)=split('::',$quantifAnnot);
				my ($labelType)=($labelTypeStrg =~ /^LABEL=(.+)/);
				if ($labelType =~ /FREE|NONE/) { # no labeling in analysis
					@{$hierarchy{"$anaID,0"}}=(@anaHierarchy,'label','Label-free');
					$position{"$anaID,0"}=++$i;
				}
				else {
					my $isobarModID;
					foreach my $infoStrg (@labelInfo) {
						my ($targetPos,$chanName,$labelStrg)=split(';',$infoStrg);
						my $obsCode; #="$anaID,$targetPos";
						my $targetName;
						if ($labelType=~/SILAC/i) {
							$obsCode="$anaID,-1";
							$targetName="$chanName ";
							$targetName.='[No label]' unless $labelStrg; # poorly annotated unlabeled channel
							my $first=1;
							#my @mods;
							my %mods; # hash because same mod can be used on different aa (mod is repeated)
							foreach my $modStrg (split ('@',$labelStrg)) { # multiple mods for same label channel
								my @labelData=split('#',$modStrg);
								if ($first) {$first=0;}
								else {$targetName.='+';}
								if ($labelData[1] =~ /No label/i) {
									$targetName.="[$labelData[1]]";
								}
								else {
									$targetName.="[$labelData[1]&raquo;$labelData[2]]"; #&raquo; &middot; &rarr;
									my $modID=($labelData[4])? $labelData[4] : &promsMod::getModificationIDfromString($dbh,$labelData[1]); # $labelData[4] if MassChroQ XIC
									$modID=0 unless $modID;
									#$obsCode.=",$modID";
									#push @mods,$modID;
									$mods{$modID}=1;
								}
							}
							if (scalar keys %mods) { # label channel
								$obsCode.=','.join(',',sort{$a<=>$b} keys %mods); # order is needed for later compare with %expState
							}
						}
						elsif ($labelType=~/ITRAQ|TMT/i) {
							$obsCode="$anaID,$targetPos";
							$targetName=$chanName;
							unless ($isobarModID) { # only once
								$sthLA->execute($anaID);
								($isobarModID)=$sthLA->fetchrow_array;
							}
							$obsCode.=",$isobarModID" if $isobarModID;
						}
						@{$hierarchy{$obsCode}}=(@anaHierarchy,'label',$labelType,'target',$targetName);
						$position{$obsCode}=++$i;
					}
				}
			}
			##>No peptide quantification associated with Analysis => guess from ANALYSIS_MODIFICATION table
			else {
				$noQuantifAna{$anaID}=1; # only used if action is 'editAE'
				$sthLA->execute($anaID);
				my $isLabeled=0;
				my $obsCode;
				my $channelPos=0;
				while (my ($modID,$psiName,$intName,$synName)=$sthLA->fetchrow_array) {
					my $labelName=($psiName)? $psiName : ($intName)? $intName : (split(/##/,$synName))[1];
					if ($labelName=~/ITRAQ|TMT/i) {
						$isLabeled=1;
						my %targetInfo;
						my $targetPos=0;
						my %infoSubmit=&promsMod::getSearchParam($dbh,$anaID);
						if ($infoSubmit{'q:Quantification'}) {
							my $first=1;
							foreach my $targetName (split(':',$infoSubmit{'q:Quantification'})) {
								if ($first) { # skip quantif Name
									$first=0;
									next;
								}
								$targetInfo{++$targetPos}=$targetName;
							}
						}
						else { # backup solution
							foreach my $name ($psiName,$intName,$synName) {
								next unless $name;
								if ($name=~/(\d+)\s*plex/i) {
									my $plexNum=$1;
									if ($labelName=~/ITRAQ/i) {
										foreach my $targetPos (1..$plexNum) {
											$targetInfo{$targetPos}=112+$targetPos;
											$targetInfo{$targetPos}=121 if $targetInfo{$targetPos}==120; # replaces 120 with 121 (120 does not exist...?)
										}
									}
									else { # TMT
										foreach my $targetPos (1..$plexNum) {
											$targetInfo{$targetPos}=125+$targetPos;
										}
									}
									last;
								}
							}
						}
						foreach my $targetPos (sort{$a<=>$b} keys %targetInfo) {
							@{$hierarchy{"$anaID,$targetPos,$modID"}}=(@anaHierarchy,'label',$labelName,'target',$targetInfo{$targetPos});
							$position{"$anaID,$targetPos,$modID"}=++$i;
						}
						next;
					}
					unless ($isLabeled) { # 1st loop: create a non-labeled channel
						$channelPos++;
						@{$hierarchy{"$anaID,$channelPos"}}=(@anaHierarchy,'label','No label');
						$position{"$anaID,$channelPos"}=++$i;
						$isLabeled=1;
					}
					# !!!!!!!!!! TODO: Not compatible with single state defined by multiple modifications !!!!!!!
					$channelPos++;
					@{$hierarchy{"$anaID,$channelPos,$modID"}}=(@anaHierarchy,'label',$labelName);
					$position{"$anaID,$channelPos,$modID"}=++$i;
				}
				unless ($isLabeled) { # no labeling in analysis
					@{$hierarchy{"$anaID,0"}}=(@anaHierarchy,'label','Label-free');
					$position{"$anaID,0"}=++$i;
				}
			}

			# Get whether or not this analysis is already associated to an EXPCONDITION & compute ExpState content summay!
			$sthObs->execute($anaID);
			while (my ($expStateID,$targetPos,$obsID,$fracGroup,$techRepGroup,$modCode)=$sthObs->fetchrow_array) {
				$targetPos=-1 if ($quantifAnnot && $quantifAnnot=~/SILAC/); # back compatibility for Design before 31/01/14
				my $obsCode="$anaID,$targetPos";
				$obsCode.=",$modCode" if $modCode;
				$expState{$obsCode}=$expStateID; # if $expStateID;
				$fractionGroup{$obsCode}=$fracGroup;
				$numFractionsInGroup{$expStateID}{$fracGroup}++ if $fracGroup;
				$technicalRepGroup{$obsCode}=$techRepGroup;
				$numTechRepsInGroup{$expStateID}{$techRepGroup}++ if $techRepGroup;
				$usedTechRepGroups{$expStateID}{$techRepGroup}=1 if $techRepGroup;
				$sthGetQ2->execute($expStateID,$anaID);
				($numAnaCondQuantif{$obsCode})=$sthGetQ2->fetchrow_array;
				$stateChildren{$expStateID}+=$numAnaCondQuantif{$obsCode};
				if ($action eq 'editAE' && $trueItem eq 'EXPCONDITION' && $trueItemID == $expStateID) {
					$maxFracGroup=$fracGroup if ($fracGroup && $fracGroup > $maxFracGroup);
					$maxTechRepGroup=$techRepGroup if ($techRepGroup && $techRepGroup > $maxTechRepGroup);
				}
				#>Summary
				$obs{$expStateID}++; # $obsID always defined & unique
				my $usedFracGr=$fracGroup || 'O'.$obsID;
				$pool{$expStateID}{$usedFracGr}++;
				$numPools{$expStateID}++ if $pool{$expStateID}{$usedFracGr}==2;
				my $usedTechRepGr=$techRepGroup || 'F'.$usedFracGr;
				$techRepGr{$expStateID}{$usedTechRepGr}++;
			}
		}
		$sthGAI->finish;
	}
	$maxFracGroup+=20; $maxTechRepGroup+=20; # !!! +20 is arbitrary !!!
	$sthLA->finish;
	$sthLQ->finish;
	$sthGetQ2->finish;
	$sthObs->finish;


	my $stateCount=0;
	while (my ($expStateID,$namec,$desc) = $sthExpCond->fetchrow_array) {
		$stateCount++;
		###> For each state, get the number of Quantitations to adapt javascript alert when deletion
		$sthgetQ1->execute($expStateID);
		my ($numQuanti)=$sthgetQ1->fetchrow_array;
		$desc=&promsMod::HTMLcompatible($desc);
		$expStateStrg.="\n<DIV id=\"state$stateCount\">\n";
		###> Start form
		if ($trueItem eq 'EXPCONDITION' && $trueItemID == $expStateID){
			if ($action eq 'edit') {
				#($parentID)=$dbh->selectrow_array("SELECT ID_DESIGN FROM EXPCONDITION WHERE ID_EXPCONDITION=$expStateID");
				$expStateStrg.=qq
|<FORM name="addItemForm" method="post">
<INPUT type="hidden" name="ACT" value="update">
<INPUT type="hidden" name="ITEM" value="EXPCONDITION">
<INPUT type="hidden" name="ITEM_ID" value="$expStateID">
<INPUT type="hidden" name="PARENT" value="DESIGN">
<INPUT type="hidden" name="PARENT_ID" value="$designID">
<INPUT type="hidden" name="PROJECT_ID" value="$projectID">
<INPUT type="hidden" name="RESTRICT_BRANCH" value="$restrictBranch">
|;
			}
			elsif ($action eq 'editAE') {
				$expStateStrg.=qq
|<FORM name="editObsExpForm" method="post">
<INPUT type="hidden" name="ACT" value="store">
<INPUT type="hidden" name="ITEM" value="EXPCONDITION">
<INPUT type="hidden" name="PARENT_ID" value="$itemID">
<INPUT type="hidden" name="PARENT" value="DESIGN">
<INPUT type="hidden" name="PROJECT_ID" value="$projectID">
<INPUT type="hidden" name="ID" value="$expStateID">
<INPUT type="hidden" name="RESTRICT_BRANCH" value="$restrictBranch">
<INPUT type="hidden" name="OBSEXPCONDVALUES" value="">
</FORM>
|;
			}
		}
		#>Summary string
		my $numObs=$obs{$expStateID} || 0;
		my $fractSamp=$numPools{$expStateID} || 0;
		my $fullSamp=($pool{$expStateID})? (scalar keys %{$pool{$expStateID}}) - $fractSamp : 0;
		my $numBioRep=($techRepGr{$expStateID})? scalar keys %{$techRepGr{$expStateID}} : 0;
		my $summaryStrg='';
		if ($numBioRep > 1) {
			$summaryStrg.="$numBioRep biological replicates";
			#$summaryStrg.='s' if $numBioRep > 1;
		}
		if ($fractSamp) {
			$summaryStrg.=' / ' if $summaryStrg;
			$summaryStrg.="$fractSamp fractionnated sample";
			$summaryStrg.='s' if $fractSamp > 1;
		}
		if ($fullSamp) {
			$summaryStrg.=' / ' if $summaryStrg;
			$summaryStrg.="$fullSamp sample";
			$summaryStrg.='s' if $fullSamp > 1;
		}
		$summaryStrg.=' / ' if $summaryStrg;
		$summaryStrg.=($numObs)? "$numObs observation" : 'Empty';
		$summaryStrg.='s' if $numObs > 1;


		$expStateStrg.="<BR><TABLE border=0 width=1200 id=\"$expStateID\" cellpadding=2 cellspacing=2 bgcolor=$dark>\n";

		if ($action eq 'edit' && $trueItem eq 'EXPCONDITION' && $trueItemID == $expStateID) {
			$expStateStrg.=qq
|<TR><TH align=left valign=top class="title2" nowrap>&nbsp;State #$stateCount :</TH><TH align=left>$summaryStrg</TH></TR>
<TR><TH align=right valign=top width=130 nowrap>Name :</TH><TD valign=left bgcolor=$light><INPUT type="text" name="name" size="50" value="$namec"></TD></TR>
<TR><TH align=right valign=top nowrap>Description :</TH><TD valign=left bgcolor=$light><TEXTAREA name="des" rows="2" cols="65">$desc</TEXTAREA></TD></TR>
<TR><TD align=left></TD><TD><INPUT type="submit" value="Save">&nbsp;<INPUT type="button" value="Cancel" onclick="cancelAction()"></TD></TR>
|;
		}
		else {
			$expStateStrg.=qq
|<TR><TH align=left valign=top class="title2" nowrap>&nbsp;State #$stateCount :</TH><TH align=left>$summaryStrg</TH></TR>
<TR><TH align=right valign=top width=130 nowrap>Name :</TH><TD valign=left bgcolor=$light>&nbsp;$namec</TD></TR>
<TR><TH align=right valign=top nowrap>Description :</TH><TD valign=left bgcolor=$light>&nbsp;$desc</TD></TR>
<TR><TH align=right valign=top nowrap>$observationStrg</TH><TD valign=left bgcolor=$light>
|;
			my $maxHierarchyIndex=0;
			my $maxNumGroups=scalar keys %position;
			foreach my $refH (values %hierarchy) {$maxHierarchyIndex=$#{$refH} if $maxHierarchyIndex < $#{$refH};}
			my $obsString='';
			my @prevHierarchy=();
			my $chkIdx=-1;
			my %checkedObsCodes; # num obscode in cond if "Manage Observations"
			if ($action eq 'editAE' && $trueItem eq 'EXPCONDITION' && $trueItemID == $expStateID) {
				foreach my $obsCode (keys %position) {
					$checkedObsCodes{$obsCode}=1 if ($expState{$obsCode} && $expState{$obsCode}==$expStateID);
				}
			}
			my $bgColor=$dark;
my %fracGr2bioSamp;
			foreach my $obsCode (sort{$position{$a}<=>$position{$b}} keys %position) {
				my ($anaID,$targetPos,@mods)=split(',',$obsCode);
				if (($action eq 'editAE' && $trueItem eq 'EXPCONDITION' && $trueItemID == $expStateID && !$expState{$obsCode}) || ($expState{$obsCode} && $expState{$obsCode}==$expStateID)) {
					#my $anaCode=($validStatus{$anaID}==-1)? 'analysis:no_scan' : ($validStatus{$anaID}==0)? 'analysis:no_val' : ($validStatus{$anaID}==1)? 'analysis:part_val' : 'analysis:val';
					my ($disabStrg,$visStrg);
					if ($expState{$obsCode} && $expState{$obsCode} != $expStateID) { # Cannot associate a already associated analysis
						$disabStrg=' disabled';
						$visStrg=' style="visibility:hidden"';
					}
					else {
						$disabStrg='';
						$visStrg='';
					}
					my $checkStrg=($checkedObsCodes{$obsCode})? ' checked' : '';
					#if ($expState{$obsCode} && $expState{$obsCode} == $expStateID) { #}
					if ($checkedObsCodes{$obsCode}) {
						$checkStrg=' checked';
						$disabStrg=' disabled' if $numAnaCondQuantif{"$expStateID,$anaID"};
					}
					if (!$prevHierarchy[0] || $prevHierarchy[1] ne $hierarchy{$obsCode}[1]) {
						$bgColor=($bgColor eq $dark)? $light : $dark;
					}
					$obsString.="<TR class=\"list\" bgcolor=\"$bgColor\">";
					for (my $i=0;$i<$#{$hierarchy{$obsCode}}; $i+=2) {
						$obsString.="<TD valign='bottom' nowrap>";
						my $hideItem=1;
						for (my $j=0; $j<=$i; $j+=2) {
							if (!$prevHierarchy[$j] || $prevHierarchy[$j+1] ne $hierarchy{$obsCode}[$j+1]) {
								$hideItem=0;
								last;
							}
						}
						unless ($hideItem) {
							#my $itemCode=($hierarchy{$obsCode}[$i] eq 'analysis')? $anaCode : ($hierarchy{$obsCode}[$i] eq 'label' && $hierarchy{$obsCode}[$i+1]=~/No label|free/i)? 'no_label' : $hierarchy{$obsCode}[$i];
							#$obsString.="<IMG src=\"$promsPath{images}/$itemIcones{$itemCode}\">&nbsp;" if $itemIcones{$itemCode};
							$obsString.="<B>".$hierarchy{$obsCode}[$i+1];
							$obsString.='<SUP style="color:#DD0000">*</SUP>' if ($noQuantifAna{$anaID} && $i==$#{$hierarchy{$obsCode}}-3); # flag analysis
							$obsString.="&nbsp;>" if $i < $#{$hierarchy{$obsCode}}-1;
							$obsString.="</B>";
						}
						$obsString.="</TD>";
						if ($i==$#{$hierarchy{$obsCode}}-1) {
							for (my $j=$i+2; $j<$maxHierarchyIndex; $j+=2) {$obsString.="<TD></TD>";}
							$obsString.="<TD valign='bottom' nowrap><INPUT type=\"checkbox\" name=\"obsList\" value=\"$obsCode\" onclick=\"extendObservationCheck(".++$chkIdx.")\" $checkStrg$visStrg$disabStrg>&nbsp;</TD>" if ($action eq 'editAE' && $trueItem eq 'EXPCONDITION' && $trueItemID == $expStateID);
						}
					}
					if ($maxNumGroups > 1) {
						##>Pool groups & technical replicates
						if ($action eq 'editAE' && $trueItem eq 'EXPCONDITION' && $trueItemID == $expStateID) {
							my ($disabStrg,$visStrg);
							if ($expState{$obsCode} && $expState{$obsCode} != $expStateID) { # Cannot associate a already associated analysis
								$disabStrg=' disabled';
								$visStrg=' style="visibility:hidden"';
							}
							else {
								$disabStrg=($numAnaCondQuantif{"$expStateID,$anaID"})? ' disabled' : '';
								$visStrg='';
							}
							#>Pool groups
							$obsString.="<TD><SELECT name=\"fractionGroup\" onchange=\"extendFractionSelection($chkIdx,this.options.selectedIndex)\" $visStrg$disabStrg><OPTION value=\"-1\">No fractions</OPTION>";
							foreach my $f (1..$maxFracGroup) { # $maxNumGroups
								my $selStrg=($fractionGroup{$obsCode} && $fractionGroup{$obsCode}==$f)? ' selected' : '';
								#$obsString.="<OPTION value=\"$f\"$selStrg>Sample #$f</OPTION>";
								$obsString.="<OPTION value=\"$f\"$selStrg>Fraction of Sample #$f</OPTION>";
							}
							$obsString.="</SELECT></TD>";
							#>Technical replicates
							$obsString.="<TD><SELECT name=\"techRepGroup\" onchange=\"extendTechRepSelection($chkIdx,this.options.selectedIndex)\" $visStrg$disabStrg><OPTION value=\"-1\">No tech. replicates</OPTION>";
							foreach my $f (1..$maxTechRepGroup) { # $maxNumGroups
								my $selStrg=($technicalRepGroup{$obsCode} && $technicalRepGroup{$obsCode}==$f)? ' selected' : '';
								#$obsString.="<OPTION value=\"$f\"$selStrg>Bio. repl #$f</OPTION>";
								$obsString.="<OPTION value=\"$f\"$selStrg>Tech. repl. of BioSample #$f</OPTION>";
							}
							$obsString.="</SELECT></TD>";
						}
						else {
							#$obsString.=($fractionGroup{$obsCode})? "<TD nowrap>&nbsp;&nbsp;[Sample #$fractionGroup{$obsCode}]</TD>" : "<TD></TD>";
							#$obsString.=($technicalRepGroup{$obsCode})? "<TD nowrap>&nbsp;&nbsp;[Bio. repl. #$technicalRepGroup{$obsCode}]</TD>" : "<TD></TD>";
							#>Fractions declaration
							if ($fractionGroup{$obsCode}) {
								my $fractionStrg=($numFractionsInGroup{$expStateID}{ $fractionGroup{$obsCode} } > 1)? 'Fraction of Sample' : 'Sample';
								$obsString.="<TD nowrap>&nbsp;&nbsp;[$fractionStrg #$fractionGroup{$obsCode}]</TD>";
							}
							else {$obsString.="<TD></TD>";}
							#>Tech Reps declaration
							if ($technicalRepGroup{$obsCode}) {
								my $techRepStrg=($numTechRepsInGroup{$expStateID}{ $technicalRepGroup{$obsCode} } > 1)? 'Tech. repl. of BioSample' : 'BioSample';
								$obsString.="<TD nowrap>&nbsp;&nbsp;[$techRepStrg #$technicalRepGroup{$obsCode}]</TD>";
								$fracGr2bioSamp{ $fractionGroup{$obsCode} }{ $technicalRepGroup{$obsCode} }=1 if $fractionGroup{$obsCode};
							}
							#else {$obsString.="<TD></TD>";}
							else {
								if ($numBioRep > 1 && ($usedTechRepGroups{$expStateID} || !$numFractionsInGroup{$expStateID})) { # Auto-generate BioSample number based on fraction & techRep context
									my $bioSampNumStrg='';
									if ($fractionGroup{$obsCode} && $fracGr2bioSamp{ $fractionGroup{$obsCode} }) {$bioSampNumStrg='&nbsp;#'.$fracGr2bioSamp{ $fractionGroup{$obsCode} };}
									else {
										my $newNum=1;
										if ($usedTechRepGroups{$expStateID}) {
											$newNum=(sort{$b<=>$a} keys %{$usedTechRepGroups{$expStateID}})[0];
											$newNum++;
											$fracGr2bioSamp{ $fractionGroup{$obsCode} }=$newNum if $fractionGroup{$obsCode};
										}
										$usedTechRepGroups{$expStateID}{$newNum}=1;
										$bioSampNumStrg='&nbsp;#'.$newNum;
									}
									$obsString.="<TD>&nbsp;&nbsp;[BioSample$bioSampNumStrg]</TD>";
								}
								else {$obsString.="<TD></TD>";}
							}
						}
					}
					$obsString.="<TD width=100%></TD></TR>\n";
					@prevHierarchy=@{$hierarchy{$obsCode}};
				}
			}
			my ($autExtStrg,$maxHeight)=('','300px'); # default
			if ($action eq 'editAE' && $trueItem eq 'EXPCONDITION' && $trueItemID == $expStateID && $obsString) {
				#$autExtStrg="<INPUT type=\"checkbox\" name=\"autoExtend\" value=\"1\" checked onchange=\"autoExtFlag=this.checked;\"> <B>Auto-extend selection</B>\n";
				$autExtStrg="<B>Select by set of:</B><SELECT name=\"autoExtend\" id=\"autoExtend\">\n<OPTION value=\"1000\">All</OPTION>\n<OPTION value=\"1\" selected>1</OPTION>\n"; # onchange="autoExtFlag=this.value*1"
				foreach my $i (2..20) {$autExtStrg.="<OPTION value=\"$i\">$i</OPTION>";}
				$autExtStrg.="<OPTION value=\"1000\">All</OPTION></SELECT>\n";
				$autExtStrg.='&nbsp;&nbsp;&nbsp;<FONT style="font-weight:bold;color:#DD0000"><SUP>*</SUP>WARNING: Missing peptide quantification(s)! Labeling channels listed may not be accurate.</FONT>' if scalar keys %noQuantifAna;
				$maxHeight='1000px';
			}
			$expStateStrg.=($obsString)? "$autExtStrg<DIV style=\"max-height:$maxHeight;overflow:auto\"><TABLE border=0 cellspacing=0>\n$obsString</TABLE></DIV>" : '&nbsp;None.';
			$expStateStrg.="</TD></TR>\n";

			if ($action eq 'editAE' && $trueItem eq 'EXPCONDITION' && $trueItemID == $expStateID){
				$expStateStrg.="<TR><TD align=left width=110></TD><TD align=left>";
				$expStateStrg.="<INPUT type=\"button\" value=\"Save associations\" onclick=\"addObsExpstate()\"><INPUT type=\"button\" value=\"Cancel\" onclick=\"cancelAction()\"></TD></TR>\n";
			}
			else {
				#my $disabDelStrg=($userStatus =~/mass|bioinfo|manag/ && !$stateChildren{$expStateID})? '' : ' disabled';
				#my $disabObsStrg=($userStatus =~/mass|bioinfo|manag/)? '' : ' disabled';
				my $disabDelStrg=($projectAccess eq 'guest' || $stateChildren{$expStateID})? ' disabled' : '';
				my $disabModifStrg=($projectAccess eq 'guest')? ' disabled' : '';
				$expStateStrg.=qq
|<TR><TD align=left width=110></TD><TD><INPUT type="button" value="Edit" onclick="editExpstate($expStateID)"$disabModifStrg>
<INPUT type="button" value="Delete" onclick="deleteExpstate($expStateID,$numQuanti)"$disabDelStrg>
<INPUT type="button" value="Manage Observations" onclick="updateAnaExpstate($expStateID)"$disabModifStrg></TD></TR>
|;
			}
		}
		$expStateStrg.="</TABLE>\n";

		if ($action eq 'edit' && $trueItem eq 'EXPCONDITION' && $trueItemID == $expStateID) {$expStateStrg.="</FORM>\n";}
		$expStateStrg.="</DIV>\n";

	}
	$sthExpCond->finish;
	$sthgetQ1->finish;
	#$sthGetQ2->finish;

	if ($action ne 'edit') {
$expStateStrg.= qq
|<FORM name="expCondForm" method="POST">
<INPUT type="hidden" name="ACT" value="store">
<INPUT type="hidden" name="ITEM" value="EXPCONDITION">
<INPUT type="hidden" name="PARENT_ID" value="$itemID">
<INPUT type="hidden" name="PARENT" value="DESIGN">
<INPUT type="hidden" name="PROJECT_ID" value="$projectID">
<INPUT type="hidden" name="RESTRICT_BRANCH" value="$restrictBranch">
<INPUT type="hidden" name="EXPCONDVALUES" value="">
</FORM>
|;
		$startNewStates=++$stateCount;
		while ($stateCount <= $MAX_NUM_STATES) {
			$expStateStrg.=qq
|<DIV id="state$stateCount" style="display:none">
<TABLE border=0 width=1200 cellpadding=2 cellspacing=2 bgcolor=$dark>
<TR><TH align=left valign=top class="title2" nowrap colspan=2>&nbsp;State #$stateCount</TH></TR>
<TR><TH align=right valign=top width=130 nowrap>Name :</TH><TD align=left bgcolor=$light><INPUT type="text" id="name$stateCount" name="name" size="55"></TD></TR>
<TR><TH align=right width=110 valign=top>Description :</TH><TD align=left bgcolor=$light><TEXTAREA name="description" id="description$stateCount" rows="2" cols="65"></TEXTAREA></TD></TR>
<TR><TD align=left width=110><TD valign=left><INPUT type="button" value="Cancel" onclick="showHideStates('hide',$stateCount)"></TD></TR>
</TABLE><BR>
</DIV>
|;
			$stateCount++;
		}
		my $disabAddStrg=($projectAccess eq 'guest')? ' disabled' : '';
		$expStateStrg.=qq
|<INPUT type="button" name="addExpState" value="Add new State" onclick="showHideStates('show',-1)"$disabAddStrg>
<INPUT type="button" id="addAllExpState" value="Save all States" onclick="addExpstate()" disabled>
|;
	}
	if ($action =~ /edit/ && $trueItem eq 'EXPCONDITION') {
		$expStateStrg.=qq
|<SCRIPT language="Javascript">
document.getElementById('$rowID').scrollIntoView(true);
</SCRIPT>
|;
	}

	if ($item eq 'DESIGN') {
		($parentID)=$dbh->selectrow_array("SELECT ID_EXPERIMENT FROM DESIGN WHERE ID_DESIGN=$itemID");
	}
	elsif ($item eq 'EXPCONDITION') {
		#($parentID)=$dbh->selectrow_array("SELECT ID_DESIGN FROM EXPCONDITION WHERE ID_EXPCONDITION=$itemID");
		$parentID=$designID;
	}
}
elsif ($action eq 'add') {
	$itemID=$rowID;
	@nameText=('Name','Description');
	#@itemInfo=&promsMod::getItemInfo($dbh,$parentItem,$rowID,1);
	#push @itemInfo,{'ITEM'=>$item,'TYPE'=>$itemType,'ID'=>undef,'NAME'=>''};
}
elsif ($action eq 'store') { # Store the new DESIGN or EXPCONDITION and show directly the summary
	if ($item eq 'DESIGN') {
		$name=param('name');
		$description=param('des');
		($newItemID)=$dbh->selectrow_array("SELECT MAX(ID_DESIGN) FROM DESIGN");
		my $sthAdd=$dbh->prepare("INSERT INTO DESIGN (ID_DESIGN,ID_EXPERIMENT,NAME,DES,UPDATE_DATE,UPDATE_USER) VALUES (?,?,?,?,NOW(),?)");
		$sthAdd->execute(++$newItemID,$parentID,$name,$description,$userID) || die $dbh->errstr();
		$designID=$newItemID;
		$itemID=$newItemID;
	}
	else {# one or more expstate added to the design
		$itemID=$rowID;
		if (defined(param('OBSEXPCONDVALUES'))){ # Association of analysis (could be an empty '')
			my (%obsChecked,%fractionGroup,%technicalRepGroup,%obsInDB);
			#$numExpCond=1; # any non-0 value is OK

			###> 1st Get checked observations
			my $obsList=param('OBSEXPCONDVALUES');
			foreach my $obsData (split(/:/,$obsList)) {
				my ($obsCode,$fracGroup,$techRepGroup)=split(/\./,$obsData);
#print "CHK_OBS='$obsCode' ($fracGroup)<BR>\n";
				#my ($anaID,$targetPos,@mods)=split(/,/,$obsCode);
				#@{$obsChecked{"$anaID,$targetPos"}}=($fracGroup,\@mods);
				@{$obsChecked{$obsCode}}=($fracGroup,$techRepGroup);
			}
			###> 2nd Get observations in DB associated to that state
			#my $sthGetObs=$dbh->prepare("SELECT O.ID_OBSERVATION,ID_ANALYSIS,TARGET_POS,GROUP_CONCAT(ID_MODIFICATION ORDER BY ID_MODIFICATION SEPARATOR ','),FRACTION_GROUP FROM OBSERVATION O LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION WHERE ID_EXPCONDITION=$itemID GROUP BY O.ID_OBSERVATION");
			my $sthGetObs=$dbh->prepare("SELECT O.ID_OBSERVATION,O.ID_BIOSAMPLE,O.ID_ANALYSIS,O.TARGET_POS,GROUP_CONCAT(ID_MODIFICATION ORDER BY ID_MODIFICATION SEPARATOR ','),OE.FRACTION_GROUP,OE.TECH_REP_GROUP
											FROM OBSERVATION O
											LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION
											JOIN OBS_EXPCONDITION OE ON O.ID_OBSERVATION=OE.ID_OBSERVATION
											WHERE OE.ID_EXPCONDITION=$itemID GROUP BY O.ID_OBSERVATION");
			my $sthDelOEC=$dbh->prepare("DELETE FROM OBS_EXPCONDITION WHERE ID_OBSERVATION=? AND ID_EXPCONDITION = $itemID");
			my $sthOmultiE=$dbh->prepare("SELECT COUNT(*) FROM OBS_EXPCONDITION WHERE ID_OBSERVATION=?");
			my $sthDelOM=$dbh->prepare("DELETE FROM OBS_MODIFICATION WHERE ID_OBSERVATION=?");
			my $sthDelObs=$dbh->prepare("DELETE FROM OBSERVATION WHERE ID_OBSERVATION=?");
			my ($restrictDisplay,%anaRestrictList);
			if ($restrictBranch =~ /sample:(\d+)/) {
				my $sampID=$1;
				$restrictDisplay=1;
				my $sthARL=$dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS A,SAMPLE S WHERE S.ID_SAMPLE=$sampID AND A.ID_SAMPLE=$sampID AND ID_SPOT IS NULL AND VALID_STATUS > 0");
				$sthARL->execute;
				while (my ($anaID)=$sthARL->fetchrow_array) {$anaRestrictList{$anaID}=1;}
				$sthARL->finish;
			}
			$sthGetObs->execute;
			while (my ($obsID,$bioSampID,$anaID,$targetPos,$modCode,$fracGroup,$techRepGroup)=$sthGetObs->fetchrow_array) {
				next if ($restrictDisplay && !$anaRestrictList{$anaID}); # skip all obsExp that where not displayed
				my $obsCode="$anaID,$targetPos";
				$obsCode.=",$modCode" if $modCode;
				$fracGroup=-1 unless $fracGroup;
				$techRepGroup=-1 unless $techRepGroup;
#print "DB_OBS='$obsCode' ($fracGroup)<BR>\n";
				@{$obsInDB{$obsCode}}=($obsID,$fracGroup,$techRepGroup);
				unless ($obsChecked{$obsCode}) { # an Observation was unchecked
					$sthDelOEC->execute($obsID); # => delete OBSERVATION/EXPCONDITION link
					unless ($bioSampID) { # if not linked to a biosample
						$sthOmultiE->execute($obsID); # check if used by other ExpCond (in other designs!)
						my ($numExpC)=$sthOmultiE->fetchrow_array;
						if ($numExpC==0) { # if not used => delete Obs
							$sthDelOM->execute($obsID);
							$sthDelObs->execute($obsID);
						}
					}
				}
			}
			$sthGetObs->finish;
			$sthOmultiE->finish;
			$sthDelOEC->finish;
			$sthDelOM->finish;
			$sthDelObs->finish;

			###> 3rd Add or update
			#my $sthInsObs=$dbh->prepare("INSERT INTO OBSERVATION (ID_EXPCONDITION,ID_ANALYSIS,TARGET_POS,FRACTION_GROUP) VALUES ($itemID,?,?,?)");
			my $sthInsObs=$dbh->prepare("INSERT INTO OBSERVATION (ID_ANALYSIS,TARGET_POS) VALUES (?,?)");
			#my $sthUpObs=$dbh->prepare("UPDATE OBSERVATION SET FRACTION_GROUP=? WHERE ID_OBSERVATION=?"); # only pool group can be updated
			my $sthUpObsExp=$dbh->prepare("UPDATE OBS_EXPCONDITION SET FRACTION_GROUP=?,TECH_REP_GROUP=? WHERE ID_OBSERVATION=? AND ID_EXPCONDITION=$itemID"); # only pool group can be updated
			my $sthSiObs0=$dbh->prepare("SELECT O.ID_OBSERVATION FROM OBSERVATION O LEFT OUTER JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION WHERE O.ID_ANALYSIS=? AND M.ID_MODIFICATION IS NULL LIMIT 0,1");
			my $sthSiObs1=$dbh->prepare("SELECT O.ID_OBSERVATION FROM OBSERVATION O,OBS_MODIFICATION M WHERE O.ID_OBSERVATION=M.ID_OBSERVATION AND ID_ANALYSIS=? AND M.ID_MODIFICATION=? LIMIT 0,1"); # SILAC Obs
			my $sthOthObs=$dbh->prepare("SELECT ID_OBSERVATION FROM OBSERVATION WHERE ID_ANALYSIS=? AND TARGET_POS=? LIMIT 0,1"); # Other labels (Label-free, iTRAQ)
			my $sthInsObsExp=$dbh->prepare("INSERT INTO OBS_EXPCONDITION (ID_EXPCONDITION,ID_OBSERVATION,FRACTION_GROUP,TECH_REP_GROUP) VALUES ($itemID,?,?,?)");
			my $sthInsOM=$dbh->prepare("INSERT INTO OBS_MODIFICATION (ID_OBSERVATION,ID_MODIFICATION) VALUES (?,?)");
			foreach my $obsCode (keys %obsChecked) {
				my $fracGroup=($obsChecked{$obsCode}[0]==-1)? undef : $obsChecked{$obsCode}[0];
				my $techRepGroup=($obsChecked{$obsCode}[1]==-1)? undef : $obsChecked{$obsCode}[1];
				if ($obsInDB{$obsCode}) {
					$sthUpObsExp->execute($fracGroup,$techRepGroup,$obsInDB{$obsCode}[0]) if ($obsInDB{$obsCode}[1] != $obsChecked{$obsCode}[0] || $obsInDB{$obsCode}[2] != $obsChecked{$obsCode}[1]); # fraction/tech rep group was changed
				}
				else {
					my ($anaID,$targetPos,@mods)=split(/,/,$obsCode);

					##>Check if Obs already defined but not linked to ExpCond
					my $obsID;
					if ($targetPos==-1) { # SILAC-type labeling
						if ($mods[0]) { # labeled channel
							$sthSiObs1->execute($anaID,$mods[0]); # 1st mod in list should be enough to find the right Obs
							($obsID)=$sthSiObs1->fetchrow_array;
						}
						else { # un-labeled channel => find Obs with no link to Mod
							$sthSiObs0->execute($anaID);
							($obsID)=$sthSiObs0->fetchrow_array;
						}
					}
					else { # Label-free, iTRAQ
						$sthOthObs->execute($anaID,$targetPos);
						($obsID)=$sthOthObs->fetchrow_array;
					}
					unless ($obsID) { # Obs not found => create new one
						$sthInsObs->execute($anaID,$targetPos);
						$obsID=$dbh->last_insert_id(undef,undef,'OBSERVATION','ID_OBSERVATION');
						foreach my $modID (@mods) { #@{$obsChecked{$obsCode}[1]}
							$sthInsOM->execute($obsID,$modID);
						}
					}
					$sthInsObsExp->execute($obsID,$fracGroup,$techRepGroup); # link Obs to ExpCond
				}
			}
			$sthUpObsExp->finish;
			$sthSiObs0->finish;
			$sthSiObs1->finish;
			$sthOthObs->finish;
			$sthInsObs->finish;
			$sthInsOM->finish;
			$sthInsObsExp->finish;
		}
		else {
			my $expCondVal=param('EXPCONDVALUES');
			($newItemID)=$dbh->selectrow_array("SELECT MAX(ID_EXPCONDITION) FROM EXPCONDITION");
			#($numExpCond)=$dbh->selectrow_array("SELECT COUNT(*) FROM EXPCONDITION WHERE ID_DESIGN=$parentID");
			my $sthInsC=$dbh->prepare("INSERT INTO EXPCONDITION (ID_EXPCONDITION,ID_DESIGN,NAME,DES) VALUES (?,?,?,?)");
			foreach my $condInfo (split(/##/, $expCondVal)) {
				($name,$description)=split(/::/,$condInfo);
				if ($name) {
					$sthInsC->execute(++$newItemID,$parentID,$name,$description) || die $dbh->errstr();
				}
			}
			$sthInsC->finish;
		}
		$designID=$parentID;
	}
	$dbh->commit;

	#@itemInfo=&promsMod::getItemInfo($dbh,$item,$newItemID,1);

}
elsif ($action eq 'delete'){#DELETE
	$itemID=$rowID;
	if ($item eq 'DESIGN') {
		###>Get information before removing the item
		($parentID)=$dbh->selectrow_array("SELECT ID_EXPERIMENT FROM DESIGN WHERE ID_DESIGN=$itemID");
		$parentItem='EXPERIMENT';

		#!!! Design is deletable only if no ExpStates
		###> Delete EXPCONDITION ASSOCIATED TO the DESIGN ABOUT TO BE REMOVED
		#my $sthDeleteAnaExpCond=$dbh->prepare("SELECT ID_EXPCONDITION FROM EXPCONDITION WHERE ID_DESIGN=$itemID");
		#$sthDeleteAnaExpCond->execute;
		#while ( (my $idExpCond)=$sthDeleteAnaExpCond->fetchrow_array) {
		#	&deleteExpstate($dbh,$idExpCond);
		#}
		#$sthDeleteAnaExpCond->finish;

		###> Delete the Design
		$dbh->do("DELETE FROM DESIGN WHERE ID_DESIGN=$itemID") || die $dbh->errstr;

	}
	elsif ($item eq 'EXPCONDITION') {
		###>Get information before removing the item
		($parentID)=$dbh->selectrow_array("SELECT ID_DESIGN FROM EXPCONDITION WHERE ID_EXPCONDITION=$itemID");
		#($numExpCond)=$dbh->selectrow_array("SELECT COUNT(*) FROM EXPCONDITION WHERE ID_DESIGN=$parentID");
		$parentItem='DESIGN';
		#&deleteExpstate($dbh,$itemID);

		###>Delete observation associations
		#$dbh->do("DELETE OBS_MODIFICATION FROM OBS_MODIFICATION INNER JOIN OBSERVATION ON OBS_MODIFICATION.ID_OBSERVATION=OBSERVATION.ID_OBSERVATION WHERE OBSERVATION.ID_EXPCONDITION=$itemID");
		#$dbh->do("DELETE FROM OBSERVATION WHERE ID_EXPCONDITION=$itemID") || die $dbh->errstr;
		my @deletableObs;
		my $sthObsExp=$dbh->prepare("SELECT O.ID_OBSERVATION FROM OBS_EXPCONDITION OE,OBSERVATION O WHERE OE.ID_OBSERVATION=O.ID_OBSERVATION AND O.ID_BIOSAMPLE IS NULL AND OE.ID_EXPCONDITION=$itemID");
		my $sthOmultiE=$dbh->prepare("SELECT COUNT(*) FROM OBS_EXPCONDITION WHERE ID_OBSERVATION=? AND ID_EXPCONDITION != $itemID");
		$sthObsExp->execute;
		while (my ($obsID)=$sthObsExp->fetchrow_array) {
			$sthOmultiE->execute($obsID);
			my ($numExpC)=$sthOmultiE->fetchrow_array;
			push @deletableObs,$obsID if $numExpC==0;
		}
		$sthObsExp->finish;
		$sthOmultiE->finish;

		$dbh->do("DELETE FROM OBS_EXPCONDITION WHERE ID_EXPCONDITION=$itemID") || die $dbh->errstr;
		my $sthDelOMod=$dbh->prepare("DELETE FROM OBS_MODIFICATION WHERE ID_OBSERVATION=?");
		my $sthDelObs=$dbh->prepare("DELETE FROM OBSERVATION WHERE ID_OBSERVATION=?");
		foreach my $obsID (@deletableObs) {
			$sthDelOMod->execute;
			$sthDelObs->execute;
		}
		$sthDelOMod->finish;
		$sthDelObs->finish;

		###>If Quantitation associated to EXPCONDITION, it is removed!
		##my $sthGetCond=$dbh->prepare("SELECT ID_QUANTIFICATION FROM EXPCONDITION_QUANTIF WHERE ID_EXPCONDITION=$itemID");
		##$sthGetCond->execute;
		##while (my ($quantifID)=$sthGetCond->fetchrow_array) {
		##	&deleteQuantification($dbh,$quantifID);
		##}
		##$sthGetCond->finish;
		$dbh->do("DELETE FROM EXPCONDITION WHERE ID_EXPCONDITION=$itemID") || die $dbh->errstr;
	}
	$dbh->commit;

}
elsif ($action eq 'update'){###> Save the update from EDIT option
	$itemID=param('ITEM_ID');
	$name=param('name');
	$description=param('des');
	my $sthUp=$dbh->prepare("UPDATE $item SET NAME=?,DES=? WHERE ID_$item=$itemID");
	$sthUp->execute($name,$description);
	$sthUp->finish;
	$dbh->commit;
}
elsif ($action eq 'duplicate') { # item is DESIGN

	###>Duplicate Design<###

	##>Fetch target Experiment data & potential Observations<##
	my (%destAnaNameID,%destObsModif,%destObsData);
	my $sthDestData=$dbh->prepare(qq|SELECT A.ID_ANALYSIS,A.NAME,O.ID_OBSERVATION,O.ID_BIOSAMPLE,O.TARGET_POS,M.ID_MODIFICATION FROM ANALYSIS A
										LEFT JOIN OBSERVATION O ON A.ID_ANALYSIS=O.ID_ANALYSIS
										LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION
										INNER JOIN SAMPLE S ON A.ID_SAMPLE=S.ID_SAMPLE
										WHERE S.ID_EXPERIMENT=?|
								);
	$sthDestData->execute($parentID);
	while (my ($anaID,$anaName,$obsID,$bioSampID,$targetPos,$modID)=$sthDestData->fetchrow_array) {
		$destAnaNameID{$anaName}=$anaID;
		if ($obsID) { # existing Obs for Analysis
			@{$destObsData{$anaID}{$obsID}}=($targetPos,$bioSampID); # bioSampID can be undef
			push @{$destObsModif{$obsID}},$modID if $modID;
		}
	}
	$sthDestData->finish;

	my $sthSrcExp=$dbh->prepare("SELECT NAME,DISPLAY_POS FROM EXPERIMENT WHERE ID_EXPERIMENT=?");
	my $sthSrcDes=$dbh->prepare("SELECT NAME,DES FROM DESIGN WHERE ID_DESIGN=?");
	my $sthAdd=$dbh->prepare("INSERT INTO DESIGN (ID_DESIGN,ID_EXPERIMENT,NAME,DES,UPDATE_DATE,UPDATE_USER) VALUES (?,?,?,?,NOW(),?)");
	my ($newCondID)=$dbh->selectrow_array("SELECT MAX(ID_EXPCONDITION) FROM EXPCONDITION");
	my $sthSrcC=$dbh->prepare("SELECT ID_EXPCONDITION,NAME,DES,DISPLAY_POS FROM EXPCONDITION WHERE ID_DESIGN=? ORDER BY ID_EXPCONDITION");
	my $sthInsC=$dbh->prepare("INSERT INTO EXPCONDITION (ID_EXPCONDITION,ID_DESIGN,NAME,DES,DISPLAY_POS) VALUES (?,?,?,?,?)");
	my $sthSrcData=$dbh->prepare(qq|SELECT A.ID_ANALYSIS,A.NAME,O.ID_OBSERVATION,O.ID_BIOSAMPLE,O.TARGET_POS,M.ID_MODIFICATION,OC.ID_EXPCONDITION,OC.FRACTION_GROUP,OC.TECH_REP_GROUP FROM ANALYSIS A
										INNER JOIN OBSERVATION O ON A.ID_ANALYSIS=O.ID_ANALYSIS
										LEFT JOIN OBS_MODIFICATION M ON O.ID_OBSERVATION=M.ID_OBSERVATION
										INNER JOIN OBS_EXPCONDITION OC ON O.ID_OBSERVATION=OC.ID_OBSERVATION
										INNER JOIN EXPCONDITION C ON OC.ID_EXPCONDITION=C.ID_EXPCONDITION
										WHERE C.ID_DESIGN=?|
								);
	my $sthInsObs=$dbh->prepare("INSERT INTO OBSERVATION (ID_ANALYSIS,TARGET_POS,ID_BIOSAMPLE) VALUES (?,?,?)");
	my $sthInsOM=$dbh->prepare("INSERT INTO OBS_MODIFICATION (ID_OBSERVATION,ID_MODIFICATION) VALUES (?,?)");
	my $sthInsObsExp=$dbh->prepare("INSERT INTO OBS_EXPCONDITION (ID_EXPCONDITION,ID_OBSERVATION,FRACTION_GROUP,TECH_REP_GROUP) VALUES (?,?,?,?)");

	($newItemID)=$dbh->selectrow_array("SELECT MAX(ID_DESIGN) FROM DESIGN");
	my $curDesignID=$newItemID;
	$designID=$itemID=++$newItemID;

	my (%numAnaUnmatched,%numObsUnmatched,%experimentInfo);
	foreach my $srcDesign (param('SRC_DESIGN')) {
		my ($srcExpID,$srcDesignID)=split(':',$srcDesign);
		$sthSrcDes->execute($srcDesignID);
		my ($curDesName,$curDescription)=$sthSrcDes->fetchrow_array;
		$sthAdd->execute(++$curDesignID,$parentID,$curDesName,$curDescription,$userID) || die $dbh->errstr();
		if ($curDesignID==$designID) {
			$name=$curDesName;
			$description=$curDescription;
		}

		###>Duplicate Conditions<###
		my %srcCond2DestCond;
		$sthSrcC->execute($srcDesignID);
		while (my ($srcCondID,$condName,$condDes,$condPos)=$sthSrcC->fetchrow_array) {
			$sthInsC->execute(++$newCondID,$curDesignID,$condName,$condDes,$condPos);
			$srcCond2DestCond{$srcCondID}=$newCondID;
		}

		###>Fetching source Observations<###
		my (%srcAnaName,%srcObsModif,%srcObsData);
		$sthSrcData->execute($srcDesignID);
		while (my ($anaID,$anaName,$obsID,$bioSampID,$targetPos,$modID,$expCondID,$fracGroup,$techRepGroup)=$sthSrcData->fetchrow_array) {
			$srcAnaName{$anaName}=$anaID;
			@{$srcObsData{$anaID}{$obsID}}=($targetPos,$bioSampID,$expCondID,$fracGroup,$techRepGroup); # bioSampID can be undef
			push @{$srcObsModif{$obsID}},$modID if $modID;
		}

		###>Match Analyses based on name<###
		my %srcAna2DestAna;
		foreach my $anaName (keys %destAnaNameID) {
			next unless $srcAnaName{$anaName}; # names must match!
			my $srcAnaID=$srcAnaName{$anaName};
			$srcAna2DestAna{$srcAnaID}=$destAnaNameID{$anaName};
		}

		###>Duplicate Observations<###
		foreach my $srcAnaID (sort{$a<=>$b} keys %srcObsData) {
			unless ($srcAna2DestAna{$srcAnaID}) {
				$numAnaUnmatched{$srcExpID}{$srcDesignID}++;
				$numObsUnmatched{$srcExpID}{$srcDesignID}+=scalar keys %{$srcObsData{$srcAnaID}};
				$experimentInfo{$srcExpID}{DESIGN}{$srcDesignID}=$curDesName;
				unless ($experimentInfo{$srcExpID}{INFO}) {
					$sthSrcExp->execute($srcExpID);
					(@{$experimentInfo{$srcExpID}{INFO}})=$sthSrcExp->fetchrow_array;
				}
				next;
			}
			my $destAnaID=$srcAna2DestAna{$srcAnaID};
			foreach my $srcObsID (keys %{$srcObsData{$srcAnaID}}) {
				my ($srcTargetPos,$srcBioSampID,$srcCondID,$fracGroup,$techRepGroup)=@{$srcObsData{$srcAnaID}{$srcObsID}};
				my $srcModifKey=($srcObsModif{$srcObsID})? join(':',sort{$a<=>$b} @{$srcObsModif{$srcObsID}}) : '-';
				my $obsMatchID=0;
				if ($destObsData{$destAnaID}) { # obs exist(s) in dest Ana => look for obs match
					foreach my $destObsID (keys %{$destObsData{$destAnaID}}) {
						my ($destTargetPos,$destBioSampID)=@{$destObsData{$destAnaID}{$destObsID}};
						my $destModifKey=($destObsModif{$destObsID})? join(':',sort{$a<=>$b} @{$destObsModif{$destObsID}}) : '-';
						if ($destTargetPos==$srcTargetPos && $destModifKey eq $srcModifKey) { # Match!!!
							$obsMatchID=$destObsID;
							last;
						}
					}
				}
				##>Create new Observation
				unless ($obsMatchID) {
					$sthInsObs->execute($destAnaID,$srcTargetPos,$srcBioSampID);
					$obsMatchID=$dbh->last_insert_id(undef,undef,'OBSERVATION','ID_OBSERVATION');
					@{$destObsData{$destAnaID}{$obsMatchID}}=($srcTargetPos,$srcBioSampID); # Record new Observation in case of another src design duplication
					if ($srcObsModif{$srcObsID}) {
						foreach my $modID (@{$srcObsModif{$srcObsID}}) {
							$sthInsOM->execute($obsMatchID,$modID);
							push @{$destObsModif{$obsMatchID}},$modID; # Record new Observation in case of another src design duplication
						}
					}
				}
				##>Connect (new/matched) Obs to new ExpCondition
				$sthInsObsExp->execute($srcCond2DestCond{$srcCondID},$obsMatchID,$fracGroup,$techRepGroup);
			}
		}
	}
	$sthSrcExp->finish;
	$sthSrcDes->finish;
	$sthAdd->finish;
	$sthSrcC->finish;
	$sthInsC->finish;
	$sthSrcData->finish;
	$sthInsObs->finish;
	$sthInsOM->finish;
	$sthInsObsExp->finish;
#$dbh->rollback;
	$dbh->commit;
	if (scalar keys %numObsUnmatched) {
		$duplicateWarningStrg='<FONT class="title2" color="#DD0000">WARNING:</FONT><FONT class="title3">';
		foreach my $expID (sort{$experimentInfo{$a}{INFO}[1]<=>$experimentInfo{$b}{INFO}[1]} keys %experimentInfo) {
			foreach my $desID (sort{lc($experimentInfo{$expID}{DESIGN}{$a}) cmp lc($experimentInfo{$expID}{DESIGN}{$b})} keys %{$experimentInfo{$expID}{DESIGN}}) {
				$duplicateWarningStrg.='<BR>&nbsp;&nbsp;&nbsp;-';
				if ($numAnaUnmatched{$expID}{$desID}) {
					$duplicateWarningStrg.=($numAnaUnmatched{$expID}{$desID}==1)? "$numAnaUnmatched{$expID}{$desID} Analysis and " : "$numAnaUnmatched{$expID}{$desID} Analyses and ";
				}
				$duplicateWarningStrg.="$numObsUnmatched{$expID}{$desID} Observation";
				$duplicateWarningStrg.='s' if $numObsUnmatched{$expID}{$desID} > 1;
				$duplicateWarningStrg.=" were not matched in $experimentInfo{$expID}{INFO}[0] &gt; $experimentInfo{$expID}{DESIGN}{$desID}";
			}
		}
		$duplicateWarningStrg.='</FONT>';
	}
}

($name,$description)=&promsMod::chkDef($name,$description);

if ($action eq 'summary') {
	push @valueText,'&nbsp;'.&promsMod::HTMLcompatible($name);
	push @valueText,'&nbsp;'.&promsMod::HTMLcompatible($description);
}
elsif ($action =~ /edit/ || $action eq 'add') {
	if ($item ne param('ITEM')) {# Edition of an expstate -> don't need to put an INPUT in that part of the HTML page
		push @valueText,'&nbsp;'.&promsMod::HTMLcompatible($name);
		push @valueText,'&nbsp;'.&promsMod::HTMLcompatible($description);
	}
	else {# Edition of a design
		push @valueText,"<INPUT type='text' name='name' value='$name' size='50'>";
		push @valueText,"<TEXTAREA name='des' rows='2' cols='65'>$description</TEXTAREA>";
	}
}
my @sourceDesigns;
if ($item eq 'DESIGN' && $action ne 'store' && $action ne 'duplicate') {
	$designID=$itemID;

	if ($action eq 'add') {
		my $sthSrcDes=$dbh->prepare("SELECT E.ID_EXPERIMENT,E.NAME,D.ID_DESIGN,D.NAME FROM EXPERIMENT E,DESIGN D WHERE E.ID_EXPERIMENT=D.ID_EXPERIMENT AND E.ID_PROJECT=? ORDER BY E.DISPLAY_POS,D.NAME");
		my $prevExpID=0;
		$sthSrcDes->execute($projectID);
		while (my ($expID,$expName,$desID,$desName)=$sthSrcDes->fetchrow_array) {
			next if $expID==$rowID; # skip parent Exp
			if ($expID != $prevExpID) {
				push @sourceDesigns,[$expID,$expName,[]];
				$prevExpID=$expID;
			}
			push @{$sourceDesigns[$#sourceDesigns]->[2]},[$desID,$desName];
		}
	}
}
elsif ($item eq 'EXPCONDITION' && $action ne 'store') {
	($designID)=$dbh->selectrow_array("SELECT ID_DESIGN FROM EXPCONDITION WHERE ID_EXPCONDITION=$itemID");
}
$dbh->disconnect;


################
####> HTML <####
################
my $titleString;
print header(-charset=>'utf-8'); warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
|;
$designID=0 unless $designID;

if ($action=~/store|update|duplicate/) {
	if ($item eq 'DESIGN') {
		print qq
|function refreshFrames() {
	top.promsFrame.selectedAction='summary';
	parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID&DESIGN=DESIGN:$designID&branchID=design:$designID&EXPERIMENT=$parentItem:$parentID&ACT=experiment&VIEW=quanti"; //&ISNAVFRAME=0
}
|;
		if ($duplicateWarningStrg) {
			print qq
|</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Design duplication</FONT>
<BR><BR><BR><BR>
$duplicateWarningStrg
<BR><BR>
<INPUT type="button" value=" OK " class="title3" onclick="refreshFrames()"/>
<BR><BR>
</CENTER>
</BODY>
</HTML>
|;
		}
		else {
			print qq
|refreshFrames();
</SCRIPT>
</HEAD>
</HTML>
|;
		}
	}
	else {
		#if ($action eq 'store' && $numExpCond==0) {
			print qq
|parent.optionFrame.location="$promsPath{cgi}/selectOptionQuanti.cgi?branchID=$parentItem:$parentID";
</SCRIPT>
</HEAD>
</HTML>
|;
		#}
		#else {print "parent.optionFrame.selectOption();";}
	}

	exit;
}
elsif ($action eq 'delete') {
	if ($item eq 'DESIGN') {
		print "parent.itemFrame.location=\"$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=$parentItem:$parentID&EXPERIMENT=$parentItem:$parentID&ACT=experiment&VIEW=quanti\";"; #&ISNAVFRAME=1
	}
	else { # ExpState deletion
		#if ($numExpCond<=1) {
			print "parent.optionFrame.location=\"$promsPath{cgi}/selectOptionQuanti.cgi?branchID=$parentItem:$parentID\";";
		#}
		#else {print "parent.optionFrame.selectOption();";}
	}
	print qq
|</SCRIPT>
</HEAD>
</HTML>
|;
	exit;
}
$titleString=($action eq 'add')? 'New Design' : "Design <FONT color=#DD0000>$name</FONT>";
print qq
|function showHideStates(action,position) {
	if (action == 'show'){
		for (var i = $startNewStates ; i <= $MAX_NUM_STATES ; i++ ) {
			if(document.getElementById('state'+i).style.display=='none' ){
				document.getElementById('state'+i).style.display='';
				document.getElementById('addAllExpState').disabled=false;
				break;
			}
		}
	}
	else {// 'hide a certain state'
		document.getElementById('state'+position).style.display='none';
		var nbNone=0;
		for (var i = $startNewStates ; i <= $MAX_NUM_STATES ; i++ ) {
			if(document.getElementById('state'+i).style.display=='none' ){
				nbNone++;
			}
		}
		if(nbNone == $MAX_NUM_STATES){
			document.getElementById('addAllExpState').disabled=true;
		}
	}
}
function addExpstate(){
	var expCondValues='';
	for (var i = $startNewStates ; i <= $MAX_NUM_STATES ; i++ ) {
		var condName=document.getElementById('name'+i);
		if (condName.value != ''){
			var condDesc=document.getElementById('description'+i);
			if (expCondValues) expCondValues+='##';
			expCondValues+=condName.value+'::'+condDesc.value;
		}
	}
	if (expCondValues == ''){
		alert('Please, fill the name of one state !');
		return;
	}
	// alert(expCondValues);
	document.expCondForm.EXPCONDVALUES.value=expCondValues;
	document.expCondForm.submit();
}
function deleteExpstate(myExpCond,numQuanti) {
	var confString;
	if (numQuanti) {
		confString='WARNING: This State has already been linked to a Quantification.\\nDeletion will remove all associated Quantifications.\\nProceed ? ';
	}
	else {
		confString='Delete this state ?';
	}

	if (confirm(confString)) {
		top.promsFrame.resultFrame.location="$promsPath{cgi}/manageDesignCondition.cgi?ACT=delete&ITEM=EXPCONDITION&ID="+myExpCond+"&PROJECT_ID=$projectID&PARENT=DESIGN&RESTRICT_BRANCH=$restrictBranch"; //&PARENT_ID=$itemID rowID
	}
}
function editExpstate(myExpCond) {
	window.location="$promsPath{cgi}/manageDesignCondition.cgi?ACT=edit&ITEM=EXPCONDITION&ID="+myExpCond+"&PROJECT_ID=$projectID&PARENT=DESIGN&RESTRICT_BRANCH=$restrictBranch";
}
function updateAnaExpstate(myExpCond) {
	window.location="$promsPath{cgi}/manageDesignCondition.cgi?ACT=editAE&ITEM=EXPCONDITION&ID="+myExpCond+"&PROJECT_ID=$projectID&PARENT=DESIGN&RESTRICT_BRANCH=$restrictBranch";
}
function updateDisplayRestriction(restrictBranch) {
	window.location="$promsPath{cgi}/manageDesignCondition.cgi?ACT=$action&ITEM=$trueItem&ID=$trueItemID&PROJECT_ID=$projectID&PARENT=DESIGN&RESTRICT_BRANCH="+restrictBranch;
}
function addObsExpstate() {
	var obsExpCondValues='';
	var checkBoxInfo=document.getElementsByName('obsList');
	var fractGroup=document.getElementsByName('fractionGroup');
	var techRepGroup=document.getElementsByName('techRepGroup');
	if (checkBoxInfo.length) {
		for (var i = 0 ; i < checkBoxInfo.length ; i++ ) {
			if (checkBoxInfo[i].checked){
				if (obsExpCondValues) obsExpCondValues+=':';
				obsExpCondValues+=checkBoxInfo[i].value+'.'+fractGroup[i].value+'.'+techRepGroup[i].value;
			}
		}
	}
	else if (checkBoxInfo.checked) {obsExpCondValues=checkBoxInfo.value+'.'+fractGroup.value+'.'+techRepGroup.value;} // single observation
	document.editObsExpForm.OBSEXPCONDVALUES.value=obsExpCondValues;
	document.editObsExpForm.submit();
}
function extendObservationCheck(chkIdx) {
	var autoExtCount=(document.getElementById('autoExtend').value * 1) - 1;
	if (autoExtCount==0) return;
	var checkBoxInfo=document.getElementsByName('obsList');
	if (!checkBoxInfo.length) {return;}
	var chkStatus=checkBoxInfo[chkIdx].checked;
	var labelCode=checkBoxInfo[chkIdx].value.replace(/^\\d+,/,''); // remove pep quantID, keep channel & modifID(s)
	for (let i=chkIdx+1; i < checkBoxInfo.length; i++) {
		if (checkBoxInfo[i].value.replace(/^\\d+,/,'') == labelCode) { // label channels & modif must match
			checkBoxInfo[i].checked=chkStatus;
			autoExtCount--;
			if (autoExtCount==0) break;
		}
	}
/* Old code. Replaced by shorter one above
	var refObs=checkBoxInfo[chkIdx].value.split(',');
	var refChannel=refObs[1];
	var refModIdStrg='';
	if (refObs.length > 2) { // labeled channel
		refObs.shift(); refObs.shift(); // keep only modId list
		refModIdStrg=refObs.join(',');
	}
	for (let i=chkIdx+1; i < checkBoxInfo.length; i++) {
		var obs=checkBoxInfo[i].value.split(',');
		if (obs[1]==refChannel) {
			if (refChannel==0) {checkBoxInfo[i].checked=chkStatus; autoExtCount--;} // unlabeled channel
			else { // -1 (SILAC), 1..n (iTRAQ/TMT)
				obs.shift(); obs.shift(); // keep only modId list
				if (refModIdStrg==obs.join(',')) {checkBoxInfo[i].checked=chkStatus; autoExtCount--;}
			}
			if (autoExtCount==0) break;
		}
	}
*/

}
function extendFractionSelection(chkIdx,selIdx) {
	var autoExtCount=(document.getElementById('autoExtend').value * 1) - 1;
	var checkBoxInfo=document.getElementsByName('obsList');
	if (!checkBoxInfo.length) {return;}
	var labelCode=checkBoxInfo[chkIdx].value.replace(/^\\d+,/,''); // remove pep quantID, keep channel & modifID(s)
	var fractionInfo=document.getElementsByName('fractionGroup');
	for (let i=chkIdx+1; i < checkBoxInfo.length; i++) {
		if (checkBoxInfo[i].checked && checkBoxInfo[i].value.replace(/^\\d+,/,'') == labelCode) { // label channels & modif must match
			fractionInfo[i].options.selectedIndex=selIdx;
			autoExtCount--;
			if (autoExtCount==0) break;
		}
	}
}
function extendTechRepSelection(chkIdx,selIdx) {
	var autoExtCount=(document.getElementById('autoExtend').value * 1) - 1;
	if (autoExtCount==0) return;
	var checkBoxInfo=document.getElementsByName('obsList');
	if (autoExtCount==0) return;
	var labelCode=checkBoxInfo[chkIdx].value.replace(/^\\d+,/,''); // remove pep quantID, keep channel & modifID(s)
	var techRepInfo=document.getElementsByName('techRepGroup');
	for (let i=chkIdx+1; i < checkBoxInfo.length; i++) {
		if (checkBoxInfo[i].checked && checkBoxInfo[i].value.replace(/^\\d+,/,'') == labelCode) { // label channels & modif must match
			techRepInfo[i].options.selectedIndex=selIdx;
			autoExtCount--;
			if (autoExtCount==0) break;
		}
	}
}
function cancelAction() {
	top.promsFrame.selectedAction='summary';
	top.promsFrame.optionFrame.selectOption();
}
function checkAllSourceDesigns(desIdx) {
	var desSrc=document.addItemForm.SRC_DESIGN;
	if (desSrc.length) {
		var chkStatus=(desSrc[desIdx].checked)? false : true;
		var [refExpID,refDesID]=desSrc[desIdx].value.split(':');
		for (let i=desIdx; i<desSrc.length; i++) {
			let [expID,desID]=desSrc[i].value.split(':');
			if (expID==refExpID) {desSrc[i].checked=chkStatus;}
			else {break;}
		}
	}
	else {desSrc.checked=(desSrc.checked)? false : true;}
}
function submitDuplicateDesign() {
	var desSrc=document.addItemForm.SRC_DESIGN;
	var okChecked=false;
	if (desSrc.length) {
		for (let i=0; i<desSrc.length; i++) {
			if (desSrc[i].checked) {
				okChecked=true;
				break;
			}
		}
	}
	else if (desSrc.checked) {okChecked=true;}
	if (!okChecked) {
		alert('ERROR: Select at least 1 Design to be duplicated.');
		return;
	}
	document.addItemForm.ACT.value='duplicate';
	document.addItemForm.submit();
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">$titleString</FONT><BR><BR>
|;

######################
####>Writing form<####
######################
$action = 'store' unless $action ne 'add';
$action = 'update' if $action eq 'edit';

if ($action ne 'summary' && $item eq $trueItem) {
	print qq
|<FORM name="addItemForm" method="post">
<INPUT type="hidden" name="ACT" value="$action">
<INPUT type="hidden" name="ITEM" value="$item">
<INPUT type="hidden" name="PARENT" value="$parentItem">
<INPUT type="hidden" name="PARENT_ID" value="$parentID">
<INPUT type="hidden" name="PROJECT_ID" value="$projectID">
<INPUT type="hidden" name="RESTRICT_BRANCH" value="$restrictBranch">
|;
	print "<INPUT type=\"hidden\" name=\"ITEM_ID\" value=\"$itemID\">\n" if ($action eq 'edit' || $action eq 'update');

	#if ($action eq 'edit') {
	#	print "<INPUT type=\"hidden\" name=\"PARENT_ID\" value=\"$itemInfo[-2]{ID}\">\n"; # unless ($item eq 'PROJECT' || $action eq 'edit');#
	#	print "<INPUT type=\"hidden\" name=\"PARENT\" value=\"$itemInfo[-2]{ITEM}\">\n";
	#}
	#print "\n";

}

####>Printing data in table<####
my ($light,$dark)=&promsConfig::getRowColors;
print qq
|<TABLE border=0 width=800>
<TR><TD><TABLE border=0 cellpadding=2 bgcolor=$dark width=100%>
|;
for my $i (0..$#nameText) {
	$valueText[$i]="" unless $valueText[$i];
	print "<TR><TH align=right valign=top width=130>$nameText[$i] :</TH>";
	print "<TD bgcolor=$light>$valueText[$i]</TD></TR>\n";
}

if ($action ne 'summary' && $item eq param('ITEM')) {

	print "<TR><TD colspan=2 align=center>\n";
	###>SUBMIT button
	print '<INPUT type="submit" name="save" value=" Save " >',"\n"; # style="color: #A0A0A0; background-color: #FF9999"
	###>CLEAR button
	my $resetString=($action eq 'add' || $action eq 'store')? '  Clear  ' : ' Cancel changes ';
	print "&nbsp;&nbsp;&nbsp;<INPUT type='reset' value='$resetString' >";
	###>CANCEL button
	print "&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" value=\"Cancel\" onclick=\"cancelAction()\" >\n";

	print qq
|</TD></TR>
</TABLE>
</TD></TR>
|;
	if ($item eq 'DESIGN' && $action eq 'store' && scalar @sourceDesigns) {
		print qq
|<TR><TH><BR>OR<BR><BR></TH></TR>
<TR><TD><TABLE border=0 cellpadding=2 bgcolor=$dark width=100%>
	<TR>
		<TH align=right valign=top width=130>Duplicate from :</TH>
		<TD bgcolor=$light>Select 1 or more designs from the list below:
<UL style="max-height:250px;overflow:auto;margin:0;padding:5px 0px 15px 0px;">
|;
		my $desIndex=0;
		foreach my $refExp (@sourceDesigns) {
			print "<LI style=\"background-color:$dark;list-style-type:none;font-weight:bold;padding:3px 5px 3px 5px;\">Experiment: $refExp->[1] <INPUT type=\"button\" onclick=\"checkAllSourceDesigns($desIndex)\" value=\"Check/Uncheck all\"/></LI>\n";
			foreach my $refDes (@{$refExp->[2]}) {
				print "<LI style=\"list-style-type:none;padding:2px 2px 2px 10px;\"><INPUT type=\"checkbox\" name=\"SRC_DESIGN\" value=\"$refExp->[0]:$refDes->[0]\">$refDes->[1]</LI>\n";
				$desIndex++;
			}
		}
		print qq
|</UL>
<B>Warning:</B> Only Observations from Analyses with matching names will be transferred.</TD></TR>
<TR><TH colspan=2><INPUT type="button" value=" GO " onclick="submitDuplicateDesign()"/></TD></TR>
</TABLE></TD></TR>
|;
	}
	print "</TABLE>\n";
	print "</FORM>\n";
}
else {print "</TABLE></TD></TR>\n</TABLE>\n";}

print "$expStateStrg</CENTER>\n<BR><BR></BODY>\n</HTML>\n";


#############################<<< SUBROUTINES >>>##############################

#########################################################################################
###> Delete a quantification, its child/children and every information related to it <###
#########################################################################################
###sub deleteQuantification {
###
###	my ($dbh,$id)=@_;
###
###	###>1- Check if there are any children related to that QUANTIFICATION
###	my $sthGetChildren=$dbh->prepare("SELECT ID_QUANTIFICATION FROM PARENT_QUANTIFICATION WHERE ID_PARENT_QUANTIFICATION=$id");
###	$sthGetChildren->execute;
###	while (my ($childID)=$sthGetChildren->fetchrow_array) {
###		&deleteQuantification($dbh,$childID);
###		$dbh->do("DELETE FROM PARENT_QUANTIFICATION WHERE ID_PARENT_QUANTIFICATION=$id AND ID_QUANTIFICATION=$childID") || die $dbh->errstr();
###	}
###
###	###>2- Delete PEPTIDE_SET, PEPTIDE AND PROTEIN information : TODO check the compatibility with PEPTIDE_SET and QUANTIFICATION once it has been created
###	#$dbh->do("DELETE FROM PEPSET_QUANTIFICATION WHERE ID_QUANTIFICATION=$id") || die $dbh->errstr();
###	$dbh->do("DELETE FROM PEPTIDE_QUANTIFICATION WHERE ID_QUANTIFICATION=$id") || die $dbh->errstr();
###	$dbh->do("DELETE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=$id") || die $dbh->errstr();
###
###	###>3- Delete ANA_CONDITION information
###	$dbh->do("DELETE FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$id") || die $dbh->errstr();
###
###	###>4- Delete EXPCONDITION_QUANTIF information
###	$dbh->do("DELETE FROM EXPCONDITION_QUANTIF WHERE ID_QUANTIFICATION=$id") || die $dbh->errstr();
###
###	###>5- Delete the QUANTIFICATION itself
###	$dbh->do("DELETE FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$id") || die $dbh->errstr();
###
###	$dbh->commit;
###}


####>Revision history<####
# 1.3.4 Fraction and tehcRep auto-extension compatible with 3+ label channels such as iTRAQ/TMT (PP 24/03/19)
# 1.3.3 Minor change to account for possible CORRECTION tag in peptide quantif QUANTIF_ANNOT field (PP 22/03/19)
# 1.3.1m2 Improved fraction/tech. rep./biosample declaration (PP 13/11/18)
# 1.3.1m [Fix] minor bug in SOFTWARE string removal from QUANTIF_ANNOT field value (PP 13/09/18)
# 1.3.1 Multiple Designs duplication at once (PP 14/03/18)
# 1.3.0 Added option to duplicate entire Design from another Experiment (PP 18/01/18)
# 1.2.4 Max number of conditions set to 100 (PP 16/01/18)
# 1.2.3 Non-guest users can modify design (PP 09/11/17)
# 1.2.2 Added warning on analysis with no peptide quantification in case of inaccurate labeling channel listing in "Manage Observations" (PP 11/05/17)
# 1.2.1 Commented code allowing child quantification deletion upon State deletion since deletion of a State is allowed only if no child quantification (PP 10/01/17)
# 1.2.0 Replace Condition with State & and state summary (PP 06/05/16)
# 1.1.8 Display can be restricted to a Sample in case of large datasets (PP 10/03/16)
# 1.1.7 Observations management is now within a Form sent in POST (PP 08/02/16)
# 1.1.6 Technical replicate management & better ergonomy (PP 09/12/15)
# 1.1.5 New Observation management: OBS_EXPCONDITION table & Obs unicity (PP 11/07/14)
# 1.1.4 Auto extension of observation and poolGroup selection (PP 12/05/14)
# 1.1.3 Already used observations are no longer displayed in association form (PP 23/04/14)
# 1.1.2 Minor display improvement in observation list (PP 17/04/14)
# 1.1.1 Modification of OBSERVATION.TARGET_POS behavior because of potential conflict if multiple peptide quantifications (PP 31/01/14)
# 1.1.0 Bug fix for multi-mod label channel (PP 30/09/13)
# 1.0.9 PEPTIDE_SET no longer used (PP 16/09/13)
# 1.0.8 Removed quantification icons (PP 03/09/13)
# 1.0.7 Table OBS_MODIFICATION is now updated for label-quantifs (PP 16/07/13)
# 1.0.6 Added pooling management (PP 12/07/13)
# 1.0.5 Compatible with 2d Gels and labeled Analysis: TARGET_POS in table OBSERVATION (PP 04/07/13)
# 1.0.4 Modification for ANA_EXPCONDITION table (GA 08/03/13)
# 1.0.3 Replace parent.navFrame.view into top.experimentView to be homogene with openProject.cgi (GA 18/10/12)
# 1.0.2 Removal of redundant code and better addition of multiple expconditions (GA 18/04/12)
# 1.0.1 Updating the display of conditions (GA 02/04/12)
# 1.0.0 New script to handle the opening, the selection and the edition of DESIGN and EXPCONDITION (GA 03/03/2011)
