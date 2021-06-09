#!/usr/local/bin/perl -w

#############################################################################
# displayValidationHistory.cgi       1.1.6                                  #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                   #
# Contact: myproms@curie.fr                                                 #
#############################################################################
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
use POSIX qw(strftime); # to get the time
use Spreadsheet::WriteExcel;
use List::Util qw(max);
use strict;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID = $ENV{'REMOTE_USER'};
my %valTypeTitle = ('quali'=>'Qualitative auto-selection',
		'manual'=>'Manual selection&nbsp;&nbsp;&nbsp;&nbsp;',
		'comp_s'=>'Comparative auto-selection',
		'comp_f'=>'Comparative auto-selection (flagging)',
		'r_flag'=>'Flag removal',
		'filter'=>'Protein filtering',
		'r_filt'=>'Protein filter removal',
		'incexc'=>'Manual protein exclusion',
		'untrac'=>'Untraced selection&nbsp;&nbsp;&nbsp;',
		'prs'=>'PhosphoRS analysis&nbsp;&nbsp;&nbsp;&nbsp;',
		'lowP_m'=>'Lower-scoring peptide activation',
		'lowP_a'=>'Lower-scoring peptide activation',
		'varMod'=>'Change of variable-modification position'
		);

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

##########################
my $analysisID = param('ID');
my $projectID = $dbh->selectrow_array("SELECT EXPERIMENT.ID_PROJECT
				      FROM EXPERIMENT,SAMPLE,ANALYSIS
				      WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_EXPERIMENT=EXPERIMENT.ID_EXPERIMENT AND ANALYSIS.ID_ANALYSIS=$analysisID");
my @userInfos = &promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess = $userInfos[2]->{$projectID};
my $projectFullAccess=($projectAccess =~ /bioinfo|mass|manag|super/)? 1 : 0;
#my $usedUserStatus = (param('BIO')) ? 'bio' : $userInfos[1];
my $view = (param('view')) ? param('view') : ($projectFullAccess)? 'mass' : 'bio';
$view='bio' unless $projectFullAccess; # in case direct URL typing by biologist
my $format = (param('format'))? lc(param('format')) : 'html';
my ($analysisName,$msType,$fileFormat,$decoy,$minScore,$maxRank) = $dbh->selectrow_array("SELECT NAME,MS_TYPE,FILE_FORMAT,DECOY,MIN_SCORE,MAX_RANK FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");

###############################
####>Fetching history data<####
###############################
# Checking display conditions : 0 => validated but not yet reported, 1=> Validated and reported , -1 => Reported but cleared by massist, still visible by biologists until next report
# Massists and super users see validation history before report, Biologists see only reported validations
my $excludedStatus=($view eq 'mass')? -1 : 0;
my $sthGetHistory = $dbh->prepare("SELECT ID_VAL_HISTORY,STEP,STATUS,VAL_TYPE,PARAM_STRG,QUERY_VAL_STRG,PEP_VAL_STRG,PROT_VAL_STRG,START_DATE,VALID_USER
				FROM VALIDATION_HISTORY
				WHERE ID_ANALYSIS=$analysisID AND STATUS != $excludedStatus  AND STEP>0 ORDER BY STEP");
$sthGetHistory->execute;
my %history;
while( my ($id,$step,$status,$valType,$paramStrg,$queryValStrg,$pepValStrg,$protValStrg,$date,$validUser) = $sthGetHistory->fetchrow_array){
	#next if (($view eq 'mass' && $status==-1) || ($view eq 'bio' && $status==0);
	$history{$id}{'step'} = $step;
	$history{$id}{'status'} = $status;
	$history{$id}{'valType'} = $valType;
	$history{$id}{'paramStrg'} = $paramStrg;
	$history{$id}{'queryValStrg'} = $queryValStrg;
	$history{$id}{'pepValStrg'} = $pepValStrg;
	$history{$id}{'protValStrg'} = $protValStrg;
	$history{$id}{'date'} = ($date)? &promsMod::formatDate($date) : 'Unknown' ;
	$history{$id}{'validUser'} = ($validUser)? $validUser : 'Unknown' ;
}
$sthGetHistory->finish;

# Unordered actions (e.g PRS) #
my $sthUnordered = $dbh->prepare("SELECT ID_VAL_HISTORY,STEP,STATUS,VAL_TYPE,PARAM_STRG,QUERY_VAL_STRG,PEP_VAL_STRG,PROT_VAL_STRG,START_DATE,VALID_USER
				FROM VALIDATION_HISTORY
				WHERE ID_ANALYSIS=$analysisID AND STATUS != $excludedStatus  AND STEP<0");
$sthUnordered->execute;
my %unorderedHist;
while( my ($id,$step,$status,$valType,$paramStrg,$queryValStrg,$pepValStrg,$protValStrg,$date,$validUser) = $sthUnordered->fetchrow_array){
	$unorderedHist{$id}{'status'} = $status;
	$unorderedHist{$id}{'valType'} = $valType;
	$unorderedHist{$id}{'paramStrg'} = $paramStrg;
	$unorderedHist{$id}{'queryValStrg'} = $queryValStrg;
	$unorderedHist{$id}{'pepValStrg'} = $pepValStrg;
	$unorderedHist{$id}{'protValStrg'} = $protValStrg;
	$unorderedHist{$id}{'date'} = ($date)? &promsMod::formatDate($date) : 'Unknown' ;
	$unorderedHist{$id}{'validUser'} = ($validUser)? $validUser : 'Unknown' ;
}
$sthUnordered->finish;

my $sthCountQM = $dbh->prepare("SELECT COUNT(DISTINCT(QV.ID_QUERY)) FROM QUERY_VALIDATION QV, QUERY_MODIFICATION QM
						       WHERE QV.ID_QUERY=QM.ID_QUERY
						       AND QM.REF_POS_STRING IS NOT NULL
						       AND QM.REF_POS_STRING!=QM.POS_STRING
						       AND ID_ANALYSIS=?");
$sthCountQM->execute($analysisID);
my ($countQueryVMod) = $sthCountQM->fetchrow_array;
$sthCountQM->finish;

my $sthCountPM = $dbh->prepare("SELECT COUNT(DISTINCT(P.ID_PEPTIDE)) FROM PEPTIDE P, PEPTIDE_MODIFICATION PM
						       WHERE P.ID_PEPTIDE=PM.ID_PEPTIDE
						       AND REF_POS_STRING IS NOT NULL AND REF_POS_STRING NOT LIKE '##PRB_%'
							   AND ((REF_POS_STRING LIKE '%##PRB_%' AND SUBSTRING(REF_POS_STRING,1,LOCATE('#',REF_POS_STRING)-1) != POS_STRING) OR (REF_POS_STRING NOT LIKE '%##PRB_%' AND REF_POS_STRING != POS_STRING))
						       AND ID_ANALYSIS=?"); #AND REF_POS_STRING!=POS_STRING
$sthCountPM->execute($analysisID);
my ($countPeptideVMod) = $sthCountPM->fetchrow_array;
$sthCountPM->finish;

$dbh->disconnect;


##################
# Stringing data #
##################
my $prefix = ($view eq 'mass')? 'V' : 'R' ;
my $manualValidation = 0; # Boolean to check for manual validation and disable FDR display
my @stepHistory;
ID :foreach my $id (sort{ $history{$a}{'step'} <=> $history{$b}{'step'} } keys %history ){
#	print qq
#|<TR bgcolor=$bgColor><TD colspan=2 valign=top>&nbsp<FONT class="title2">$history{$id}{step}. $valTypeTitle{$history{$id}{valType}}</FONT></TD><TD>
#<B>User:</B> $history{$id}{validUser}&nbsp;&nbsp;&nbsp;<B>Date:</B> $history{$id}{date}&nbsp;</TD></TR>
#|;
	my @strgList = ();
	### Qualitative autoselect ###
	if($history{$id}{'valType'} eq 'quali' ){

		my ($selGoodInt) = ($history{$id}{'paramStrg'} =~ /selGoodInt:check:(\d*);/);
		my ($rejBadInt) = ($history{$id}{'paramStrg'} =~ /rejBadInt:check:(\d*);/);
		if($selGoodInt || $rejBadInt){
                        my $optionStrg = '';
			if($selGoodInt){
				$optionStrg .= ($rejBadInt)? "<B>Select</B> and <B>reject</B> " : "<B>Select</B> ";
			} elsif ($rejBadInt && !$selGoodInt){
				$optionStrg .= "<B>Reject</B> ";
			}
			$optionStrg .= "interpretations following these criteria:";
			push @strgList, $optionStrg;
			my ($minInt2Sc1,$minInt2Sc2,$minInt2Sc3,$minInt3Sc1,$minInt3Sc2,$minInt3Sc3,$minInt4Sc1,$minInt4Sc2,$minInt4Sc3);
			my ($minIntSc1,$minIntSc2,$minIntSc3);
			if($fileFormat =~ /SEQUEST/ ){
				($minInt2Sc1) = ($history{$id}{'paramStrg'} =~ /minInt2Sc1:text:([\w\.]*);/); if($minInt2Sc1){ push @strgList, "&nbsp&nbsp&nbsp - Minimum score if 1 peptide 2+/protein : $minInt2Sc1"};
				($minInt2Sc2) = ($history{$id}{'paramStrg'} =~ /minInt2Sc2:text:([\w\.]*);/); if($minInt2Sc2){ push @strgList, "&nbsp&nbsp&nbsp - Minimum score if 2 peptides 2+/protein : $minInt2Sc2"};
				($minInt2Sc3) = ($history{$id}{'paramStrg'} =~ /minInt2Sc3:text:([\w\.]*);/); if($minInt2Sc3){ push @strgList, "&nbsp&nbsp&nbsp - Minimum score if 3 peptides 2+/protein : $minInt2Sc3"};
				($minInt3Sc1) = ($history{$id}{'paramStrg'} =~ /minInt3Sc1:text:([\w\.]*);/); if($minInt3Sc1){ push @strgList, "&nbsp&nbsp&nbsp - Minimum score if 1 peptide 3+/protein : $minInt3Sc1"};
				($minInt3Sc2) = ($history{$id}{'paramStrg'} =~ /minInt3Sc2:text:([\w\.]*);/); if($minInt3Sc2){ push @strgList, "&nbsp&nbsp&nbsp - Minimum score if 2 peptides 3+/protein : $minInt3Sc2"};
				($minInt3Sc3) = ($history{$id}{'paramStrg'} =~ /minInt3Sc3:text:([\w\.]*);/); if($minInt3Sc3){ push @strgList, "&nbsp&nbsp&nbsp - Minimum score if 3 peptides 3+/protein : $minInt3Sc3"};
				($minInt4Sc1) = ($history{$id}{'paramStrg'} =~ /minInt4Sc1:text:([\w\.]*);/); if($minInt4Sc1){ push @strgList, "&nbsp&nbsp&nbsp - Minimum score if 1 peptide 4+/protein : $minInt4Sc1"};
				($minInt4Sc2) = ($history{$id}{'paramStrg'} =~ /minInt4Sc2:text:([\w\.]*);/); if($minInt4Sc2){ push @strgList, "&nbsp&nbsp&nbsp - Minimum score if 2 peptides 4+/protein : $minInt4Sc2"};
				($minInt4Sc3) = ($history{$id}{'paramStrg'} =~ /minInt4Sc3:text:([\w\.]*);/); if($minInt4Sc3){ push @strgList, "&nbsp&nbsp&nbsp - Minimum score if 3 peptides 4+/protein : $minInt4Sc3"};
			} else {
				($minIntSc1) = ($history{$id}{'paramStrg'} =~ /minIntSc1:text:([\w\.]*);/); if($minIntSc1){ push @strgList, "&nbsp&nbsp&nbsp - Minimum score if 1 peptide/protein : $minIntSc1"};
				($minIntSc2) = ($history{$id}{'paramStrg'} =~ /minIntSc2:text:([\w\.]*);/); if($minIntSc2){ push @strgList, "&nbsp&nbsp&nbsp - Minimum score if 2 peptides/protein : $minIntSc2"};
				($minIntSc3) = ($history{$id}{'paramStrg'} =~ /minIntSc3:text:([\w\.]*);/); if($minIntSc3){ push @strgList, "&nbsp&nbsp&nbsp - Minimum score if &#8805 3 peptides/protein : $minIntSc3"};
			}
			my ($noReject) = ($history{$id}{'paramStrg'} =~ /noReject:check:(\d*);/); if($noReject){ push @strgList, "&nbsp&nbsp&nbsp - Exclude already rejected interpretation from count."};
			my ($minSize) = ($history{$id}{'paramStrg'} =~ /minSize:text:(\w*);/);
			my ($maxSize) = ($history{$id}{'paramStrg'} =~ /maxSize:text:(\w*);/);
			if($minSize || $maxSize){
				my $sizeStrg = '&nbsp&nbsp&nbsp - ';
				$sizeStrg .= ($minSize)? "$minSize aa &#8804 " : '';
				$sizeStrg .= "size";
				$sizeStrg .= ($maxSize)? " &#8804 $maxSize aa." : '.';
				push @strgList, $sizeStrg;
			}
			my ($minDelta) = ($history{$id}{'paramStrg'} =~ /minDelta:text:(\w*);/); if($minDelta){ push @strgList, "&nbsp&nbsp&nbsp - Mass error &#8804 $minDelta."};
                        if($msType eq 'PMF'){
                                my ($newMaxRankMS1) = ($history{$id}{'paramStrg'} =~ /newMaxRankMS1:select:(\d*);/);
                                if($newMaxRankMS1){ push @strgList, "&nbsp&nbsp&nbsp - Rank &#8804 $newMaxRankMS1."};
                        } elsif($msType eq 'MIS'){
                                my ($newMaxRankMS2) = ($history{$id}{'paramStrg'} =~ /newMaxRankMS2:select:(\d*);/);
                                if($newMaxRankMS2){ push @strgList, "&nbsp&nbsp&nbsp - Rank &#8804 $newMaxRankMS2."};
                        } elsif($msType eq 'SQ'){
                                my ($newMaxRankMS1) = ($history{$id}{'paramStrg'} =~ /newMaxRankMS1:select:(\d*);/);
                                my ($newMaxRankMS2) = ($history{$id}{'paramStrg'} =~ /newMaxRankMS2:select:(\d*);/);
                                if($newMaxRankMS1 && $newMaxRankMS2){
                                        push @strgList, "&nbsp&nbsp&nbsp - Rank &#8804 $newMaxRankMS1 for PMF.<BR>\n&nbsp&nbsp&nbsp - Rank &#8804 $newMaxRankMS2 for MS/MS."
                                }
                        }
			my ($flaggedUp) = ($history{$id}{'paramStrg'} =~ /flaggedUp:check:(\d*);/); if($flaggedUp){ push @strgList, "&nbsp&nbsp&nbsp - Interpretation was <IMG src=\"$promsPath{images}/lightYellow1.gif\" hspace=0 border=0 height=11 width=11> flagged."};
			my ($flaggedDown) = ($history{$id}{'paramStrg'} =~ /flaggedDown:check:(\d*);/); if($flaggedDown){ push @strgList, "&nbsp&nbsp&nbsp - Interpretation was not <IMG src=\"$promsPath{images}/lightOrange1.gif\" hspace=0 border=0 height=11 width=11> flagged."};
		}
		my ($oneInt) = ($history{$id}{'paramStrg'} =~ /oneInt:check:(\d*);/); if($oneInt){ push @strgList, "<B>Select</B> only one interpretation per query (best matching rank)"};
		my ($oneSeqType) = ($history{$id}{'paramStrg'} =~ /oneSeqType:text:(\w*);/);
		my ($oneSeq) = ($history{$id}{'paramStrg'} =~ /oneSeq:check:(\d*);/); if($oneSeq) { push @strgList, "<B>Select</B> only one query per $oneSeqType (best score)"};
		my ($overPep) = ($history{$id}{'paramStrg'} =~ /overPep:check:(\d*);/); if($overPep){ push @strgList, "<B>Overwrite</B> previous selections/rejections."};
		my ($selProt) = ($history{$id}{'paramStrg'} =~ /selProt:check:(\d*);/); if($selProt){ push @strgList, "<B>Select</B> only proteins meeting these criteria:"};
		my ($minPep) = ($history{$id}{'paramStrg'} =~ /minPep:text:(\w*);/); if($selProt && $minPep){ push @strgList, "&nbsp&nbsp&nbsp - containing at least $minPep peptides."};
		my ($minProtSc) = ($history{$id}{'paramStrg'} =~ /minProtSc:text:([\w\.]*);/); if($selProt && $minProtSc){ push @strgList, "&nbsp&nbsp&nbsp - with score &#8805 $minProtSc."};
		my ($minCov) = ($history{$id}{'paramStrg'} =~ /minCov:text:(\w*);/); if($selProt && $minCov){ push @strgList, "&nbsp&nbsp&nbsp - with peptide coverage &#8805 $minCov%."};
		my ($bestMatch) = ($history{$id}{'paramStrg'} =~ /bestMatch:check:(\d*);/); if($bestMatch){ push @strgList, "<B>Keep</B> only best matching protein(s) for each match group."};
		my ($overProt) = ($history{$id}{'paramStrg'} =~ /overProt:check:(\d*);/); if($overProt){ push @strgList, "<B>Overwrite</B> previous exclusions."};

		#print "<TR bgcolor=$bgColor><TD width=30></TD><TD>";
		#foreach my $string(@strgList){
		#	print "$string<BR>\n";
		#}
		#print "</TD>";
	}
	### Comparative selection ###
	elsif($history{$id}{'valType'} =~ /comp/){ # comp_s or comp_f
		#my @strgList ;
		my ($reference) = ($history{$id}{'paramStrg'} =~ /COMP:(\w+.?\w+);/ ); push @strgList, "Reference: <B>$reference</B>";
		my ($actionPep) = ($history{$id}{'paramStrg'} =~ /actionPep:radio:(\w*);/);
		if($actionPep eq 'selRej'){
			my ($selectPep) = ($history{$id}{'paramStrg'} =~ /selectPep:check:(\d*);/); if($selectPep){ push @strgList, "<B>Select</B> peptides matching the option below."}
			my ($rejectPep) = ($history{$id}{'paramStrg'} =~ /rejectPep:check:(\d*);/); if($rejectPep){ push @strgList, "<B>Reject</B> peptides that do not match the option below."}
			my ($anaAction) = ($history{$id}{'paramStrg'} =~ /anaAction:radio:(\w*);/);
			if(($selectPep || $rejectPep) && $anaAction eq "inter"){
				push @strgList, "&nbsp&nbsp&nbsp - common to Reference and $analysisName.";
			} elsif(($selectPep || $rejectPep) && $anaAction eq "diff"){
				push @strgList, "&nbsp&nbsp&nbsp - unique to $analysisName.";
			}
		} elsif($actionPep eq 'flag'){
			my ($FlagUp) = ($history{$id}{'paramStrg'} =~ /FlagUp:check:(\d*);/); if($FlagUp){push @strgList, "<B>Flag</B> (<IMG src=\"$promsPath{images}/lightYellow1.gif\" hspace=0 border=0 height=11 width=11>) peptides matching the option below."}
			my ($FlagDown) = ($history{$id}{'paramStrg'} =~ /FlagDown:check:(\d*);/); if($FlagDown){push @strgList, "<B>Flag</B> (<IMG src=\"$promsPath{images}/lightOrange1.gif\" hspace=0 border=0 height=11 width=11>) peptides that do not match the option below."}
			my ($FltAction) = ($history{$id}{'paramStrg'} =~ /FltAction:radio:(\w*);/);
			if(($FlagUp || $FlagDown) && $FltAction eq 'inter'){
				push @strgList, "&nbsp&nbsp&nbsp - common to Reference and $analysisName.";
			} elsif(($FlagUp || $FlagDown) && $FltAction eq 'diff'){
				push @strgList, "&nbsp&nbsp&nbsp - unique to $analysisName.";
			}
		}
		my ($selMaxRank) = ($history{$id}{'paramStrg'} =~ /selMaxRank:select:(\d*);/); if($selMaxRank){ push @strgList, "<B>Compare</B> only peptides with rank &#8804 $selMaxRank."};
		my ($validPept) = ($history{$id}{'paramStrg'} =~ /validPept:check:(\d*);/); if($validPept){ push @strgList, "<B>Use</B> only peptides already selected in Reference."};
		my ($rejMinValPept) = ($history{$id}{'paramStrg'} =~ /rejMinValPept:check:(\d*);/);
		my ($smPepSc) = ($history{$id}{'paramStrg'} =~ /smPepSc:text:([\d\.]*);/); if($rejMinValPept && $smPepSc){ push @strgList, "<B>Ignore</B> peptides in $analysisName with score below $smPepSc."}
		my ($ignRefMinScore) = ($history{$id}{'paramStrg'} =~ /ignRefMinScore:check:(\d*);/);
		my ($refMinScore) = ($history{$id}{'paramStrg'} =~ /refMinScore:text:([\w\.]*);/); if($refMinScore && $ignRefMinScore){ push @strgList, "<B>Ignore</B> peptides in Reference with score below $refMinScore."}
		my ($discardEltTime) = ($history{$id}{'paramStrg'} =~ /discardEltTime:check:(\d*);/);
		my ($maxEltTime) = ($history{$id}{'paramStrg'} =~ /maxEltTime:text:(\d*);/); if($discardEltTime && $maxEltTime){ push @strgList, "<B>Ignore</B> peptides with elution time differing by more than $maxEltTime."}

		#print "<TR bgcolor=$bgColor><TD width=30></TD><TD>";
		#foreach my $string(@strgList){
		#	print "$string<BR>\n";
		#}
		#print "</TD>";
	}
	### Removing flags ###
	elsif ($history{$id}{'valType'} eq 'r_flag'){ # || $history{$id}{'valType'} eq 'r_filt'
		#next ID; ## only removing flags title will be displayed
		@strgList = ();
	}
	### Manual selection ###
	elsif ($history{$id}{'valType'} eq 'manual'){
		$manualValidation = 1;
		#print "<TR bgcolor=$bgColor><TD width=30></TD><TD></TD>";
		@strgList = ();
	}
	### Protein filtering ###
	elsif ($history{$id}{'valType'} eq 'filter'){
		$manualValidation = 1;
		#my @strgList;
		# If using other analyses as filter
		my ($anaFileList) = ($history{$id}{'paramStrg'} =~ /anaFileList:([^;]+);/);
		if($anaFileList){
			my @filenames = split /,/, $anaFileList;
			my $anaActionString = '';
			my ($anaAction) = ($history{$id}{'paramStrg'} =~ /anaAction:radio:(\w+);/);
			if(scalar @filenames != 1){
				$anaActionString.= "<TABLE cellpadding=0 cellspacing=0 border=0><TR>";
				$anaActionString.= ($anaAction eq 'diff')? "<TD valign='top'><B>Keep</B> proteins unique to $analysisName comparing to these analyses: &nbsp</TD><TD>" : ($anaAction eq 'inter')? "<TD valign='top'>Keep proteins common to $analysisName and these analyses: &nbsp</TD><TD>" : '';
				foreach my $filename (@filenames){
						$anaActionString.="<B>$filename</B>\n<BR>"
				}
				$anaActionString =~ s/<BR>\Z//;
				$anaActionString.= "</TD></TABLE>";
			} else {
				$anaActionString.= ($anaAction eq 'diff')? "<B>Keep</B> proteins unique to $analysisName comparing to <B>$filenames[0]</B>" : ($anaAction eq 'inter')? "<B>Keep</B> proteins common to $analysisName and <B>$filenames[0]</B>" : '';
			}
			push @strgList, $anaActionString;
			my ($rejProt) = ($history{$id}{'paramStrg'} =~ /rejProt:check:(\w*);/);
			my ($smPepSc) = ($history{$id}{'paramStrg'} =~ /smPepScore:text:(\d+);/);
			if($rejProt || $smPepSc ){
				my $nonValidStrg = "<B>If</B> non-validated Analyses are used as filter, ";
				if($rejProt){
					$nonValidStrg.= "ignore already rejected proteins";
					$nonValidStrg.= ($smPepSc)? " and " : ".";
				}
				if($smPepSc){
					$nonValidStrg.= "ignore proteins matched by a single peptide with score &#8804 $smPepSc (MS/MS analyses only)."
				}
				push @strgList, $nonValidStrg;
			}
			my ($groupAnaAction) = ($history{$id}{'paramStrg'} =~ /groupAnaAction:check:(\w*);/);
			if($groupAnaAction){
				push @strgList, "<B>Extend</B> filter to match group.";
			}
		}
		# If using organisms
		my ($orgList) = ($history{$id}{'paramStrg'} =~ /orgList:list:([^;]+);/);
		if($orgList){
			my @orgNames = split /,/, $orgList;
			my $orgString = '';
			my ($orgAction) = ($history{$id}{'paramStrg'} =~ /orgAction:radio:(\w+);/);
			if(scalar @orgNames != 1){
				$orgString.= "<TABLE cellpadding=0 cellspacing=0 border=0><TR>";
				$orgString.= ($orgAction eq 'keep')? "<TD valign='top'><B>Keep</B> proteins from these species: &nbsp</TD><TD>" : ($orgAction eq 'exclude')? "<TD valign='top'><B>Exclude</B> proteins from these species: &nbsp</TD><TD>" : '';
				foreach my $orgname (@orgNames){
						$orgString.="<B>$orgname</B>\n<BR>";
				}
				$orgString =~ s/<BR>\Z//;
				$orgString.= "</TD></TABLE>";
			} else {
				$orgString.= ($orgAction eq 'keep')? "<B>Keep</B> proteins from <B>$orgNames[0]</B>." : ($orgAction eq 'exclude')? "<B>Exclude</B> proteins from <B>$orgNames[0]</B>." : '';
			}
			push @strgList, $orgString;
		}
		# If using protein file
		my ($protFile) = ($history{$id}{'paramStrg'} =~ /protFile:file:([^;]+);/);
		if($protFile){
			my $protFileString = '<B>Keep</B> proteins ';
			my ($fileAction) = ($history{$id}{'paramStrg'} =~ /fileAction:radio:(\w+);/ );
			my ($groupFileAction) = ($history{$id}{'paramStrg'} =~ /groupFileAction:check:(\w+);/ );
			if($fileAction eq 'diff'){
				$protFileString.="unique to $analysisName comparing to <b>$protFile</b>";
			} elsif ($fileAction eq 'inter'){
				$protFileString.="common to $analysisName and <b>$protFile</b>";
			}
			if($groupFileAction){
				$protFileString.=" and extend filter to match group"
			}
			$protFileString.=".";
			push @strgList, $protFileString;
		}
		#print "<TR bgcolor=$bgColor><TD width=30></TD><TD>";
		#foreach my $string(@strgList){
		#	print "$string";
		#	print "<BR>" if ($string !~ /<\/TABLE>\Z/ ); # no <br> if string is a table
		#	print "\n";
		#}
		#print "</TD>";
	}
	elsif ($history{$id}{'valType'} eq 'incexc' || $history{$id}{'valType'} eq 'r_filt'){
		$manualValidation = 1;
		#print "<TR bgcolor=$bgColor><TD width=30></TD><TD></TD>";
	}
        elsif ($history{$id}{'valType'} eq 'untrac'){
                my ($manVal) = ( $history{$id}{'paramStrg'} =~ /manualValidation=(\d);/ );
                $manualValidation = 1 if ($manVal);
                #print "<TR bgcolor=$bgColor><TD width=30></TD><TD></TD>";
       }# } elsif($history{$id}{'valType'} eq 'prs'){
#		my ($massTolerance) = ($history{$id}{'paramStrg'} =~ /massTolerance:([^;]+)/);
#		my ($activationType) = ($history{$id}{'paramStrg'} =~ /activationType:([^;]+)/);
#		my ($threshold) = ($history{$id}{'paramStrg'} =~ /threshold:([^;]+)/);
#
#		print qq
#|<TR bgcolor=$bgColor><TD width=30></TD><TD>
#<B>Threshold</B>: $threshold%<BR>
#<B>Activation type</B>: $activationType<BR>
#<B>Mass tolerance</B>: $massTolerance
#|;
#	}

	### Displaying queries/peptides/proteins stats ###
	my @additionalInfos;
	#my $prefix = ($usedUserStatus eq 'bioinfo' || $usedUserStatus eq 'massist' || ($usedUserStatus eq 'bio' && $projectAccess =~ /super/) )? 'V' : 'R' ;
	if ($history{$id}{'valType'} eq 'comp_f'){ # valType is comp_f (flag)
		my ($FlagUp) = ($history{$id}{'pepValStrg'} =~ /$prefix\_FlagUp\=(\d+);/);
		my ($FlagDown) = ($history{$id}{'pepValStrg'} =~ /$prefix\_FlagDown\=(\d+);/);
#		print qq
#|<TD valign=top align=right><TABLE border=0 cellpadding=0 width=220><TD nowrap>
#$FlagUp peptides <IMG src="$promsPath{images}/lightYellow1.gif" hspace=0 border=0 height=11 width=11> flagged.<BR>
#$FlagDown peptides <IMG src="$promsPath{images}/lightOrange1.gif" hspace=0 border=0 height=11 width=11> flagged.
#</TD></TABLE></TD></TR>
#|;
		push @additionalInfos, qq|$FlagUp peptides <IMG src="$promsPath{images}/lightYellow1.gif" hspace=0 border=0 height=11 width=11> flagged.|;
		push @additionalInfos, qq|$FlagDown peptides <IMG src="$promsPath{images}/lightOrange1.gif" hspace=0 border=0 height=11 width=11> flagged.|;
	}
	elsif($history{$id}{'valType'} eq 'filter'){
		my ($protFiltered) = ($history{$id}{'protValStrg'} =~ /$prefix\_Filt\=(\d+);/);
		my ($protAll) = ($history{$id}{'protValStrg'} =~ /$prefix\_All\=(\d+);/);
#		print qq
#|<TD valign=top align=right><TABLE border=0 cellpadding=0 width=220><TD nowrap>
#$protFiltered/$protAll filtered proteins.
#</TD></TABLE></TD></TR>
#|;
		push @additionalInfos, "$protFiltered/$protAll filtered proteins.";
	}
	elsif($history{$id}{'valType'} eq 'r_filt'){
		my ($protValid) = ($history{$id}{'protValStrg'} =~ /$prefix\_Valid\=(\d+);/);
		my ($protAll) = ($history{$id}{'protValStrg'} =~ /$prefix\_All\=(\d+);/);
#		print qq
#|<TD valign=top align=right><TABLE border=0 cellpadding=0 width=220><TD nowrap>
#$protValid/$protAll selected proteins.
#</TD></TABLE></TD></TR>
#|;
		push @additionalInfos, "$protValid/$protAll selected proteins.";
	}
	elsif($history{$id}{'valType'} eq 'incexc'){
		my ($protExcluded) = ($history{$id}{'protValStrg'} =~ /$prefix\_Exclu\=(\d+);/);
		my ($protAll) = ($history{$id}{'protValStrg'} =~ /$prefix\_All\=(\d+);/);
#		print qq
#|<TD valign=top align=right><TABLE border=0 cellpadding=0 width=220><TD nowrap>
#$protExcluded/$protAll excluded proteins.
#</TD></TABLE></TD></TR>
#|;
		push @additionalInfos, "$protExcluded/$protAll excluded proteins.";
	}
	elsif($history{$id}{'valType'} =~ /^lowP/){
		#print "<TR bgcolor=$bgColor><TD width=30></TD><TD colspan=2>";
		if ($history{$id}{'valType'} eq 'lowP_a') {
			push @strgList, "All lower-scoring peptides activated."
		}
		elsif($history{$id}{'valType'} eq 'lowP_m'){
			push @strgList, "Some lower-scoring peptides activated."
		}
		#print "</TD></TR>\n";
	}
	else {
		my ($queryString,$pepString,$protString,$FDRstring) = ('','','','');
		if($history{$id}{'queryValStrg'}){
			my ($queryVerif) = ($history{$id}{'queryValStrg'} =~ /$prefix\_Verif\=(\d+);/);
			my ($queryAll) = ($history{$id}{'queryValStrg'} =~ /$prefix\_All\=(\d+);/);
			$queryString.= "$queryVerif/$queryAll queries verified.";
		}
		if($history{$id}{'pepValStrg'}){
			my ($pepVal) = ($history{$id}{'pepValStrg'} =~ /$prefix\_Valid\=(\d+);/);
			my ($pepDecoy) = ($history{$id}{'pepValStrg'} =~ /$prefix\_Decoy\=(\d+);/);
			$pepString.= "$pepVal peptides validated.<BR>$pepDecoy decoy peptides validated.";
			#my $FDR = ($pepDecoy)? 100*$pepDecoy/$pepVal : 0;
			#my $FDRrounded = sprintf("%.2f", $FDR);
			#$FDRstring.= (!$decoy) ? '' : ($manualValidation) ? '<B>FDR:</B> N/A (Manual validation)' : "<B>FDR:</B> $FDRrounded% ($pepDecoy decoy peptides)";
		}
		if($history{$id}{'protValStrg'}){
			my ($protValid) = ($history{$id}{'protValStrg'} =~ /$prefix\_Valid\=(\d+);/);
			my ($protAll) = ($history{$id}{'protValStrg'} =~ /$prefix\_All\=(\d+);/);
			$protString.= "$protValid/$protAll proteins selected.";
		}
#		print qq
#|<TD valign=top align=right><TABLE border=0 cellpadding=0 width=220><TD nowrap>
#$queryString
#$pepString
#$protString
#$FDRstring
#|;
		#push @additionalInfos, ($queryString,$pepString,$protString,$FDRstring);
        push @additionalInfos, ($queryString,$pepString,$protString);
                if($history{$id}{'valType'} eq 'untrac'){
                        my ($protAll) = ($history{$id}{'protValStrg'} =~ /$prefix\_All\=(\d+);/);
                        my ($protExcluded) = ($history{$id}{'protValStrg'} =~ /$prefix\_Exclu\=(\d+);/);
                        my ($FlagUp) = ($history{$id}{'pepValStrg'} =~ /$prefix\_FlagUp\=(\d+);/);
                        my ($FlagDown) = ($history{$id}{'pepValStrg'} =~ /$prefix\_FlagDown\=(\d+);/);
                        my ($protFiltered) = ($history{$id}{'protValStrg'} =~ /$prefix\_Filt\=(\d+);/);
                        #print qq|<BR>$protExcluded/$protAll excluded proteins<BR>
                        #$protFiltered/$protAll filtered proteins<BR>
                        #$FlagUp peptides <IMG src="$promsPath{images}/lightYellow1.gif" hspace=0 border=0 height=11 width=11> flagged.<BR>
                        #$FlagDown peptides <IMG src="$promsPath{images}/lightOrange1.gif" hspace=0 border=0 height=11 width=11> flagged.
                        #|;
			push @additionalInfos, ("$protExcluded/$protAll excluded proteins.",
						"$protFiltered/$protAll filtered proteins.",
						"$FlagUp peptides <IMG src=\"$promsPath{images}/lightYellow1.gif\" hspace=0 border=0 height=11 width=11> flagged.",
						"$FlagDown peptides <IMG src=\"$promsPath{images}/lightOrange1.gif\" hspace=0 border=0 height=11 width=11> flagged.");
                }
#print qq|</TD></TABLE></TD></TR>
#|;
	}

	push @stepHistory, {title => $valTypeTitle{$history{$id}{'valType'}},
			  step => $history{$id}{step},
			  user => $history{$id}{validUser},
			  date => $history{$id}{date},
			  mainInfos => \@strgList,
			  additionalInfos => \@additionalInfos};

}

##############
### Others ###
##############
my @importActions; # displayed at the beginning of the history
my @others; # displayed at the end of the history in "Other actions" section

## Import parameters ##
my $maxRankString = "<B>Max. rank</B>: $maxRank";
my $minScoreString = "<B>Min. score</B>: $minScore";
my ($decoyMethod,$maxFDR)=split(',',$decoy);
my $fdrString = "<B>False discovery</B>:";
if ($maxFDR) {
	$maxFDR =~ s/FDR=//;
	$minScoreString .= ' [FDR-based]';
	my ($FDRvalue,$FDRalgo)=split(':',$maxFDR);
	$FDRalgo=(!$FDRalgo || $FDRalgo eq 'qvality')? 'Qvality' : ($FDRalgo eq 'DTcount')? 'DT count' : 'precomputed FDR';
	$fdrString .= "<BR>&nbsp;&nbsp;- Targeted FDR: $FDRvalue% ($FDRalgo)";
	$minScoreString = "<B>Min. score</B>: N/A" if $FDRalgo =~ /precomputed/;
}
else {$minScoreString.=' [User-defined]';}

$fdrString.='<BR>&nbsp;&nbsp;- Decoy method: ';
if ($decoyMethod eq 'INT:SEARCH') {$fdrString.='Automated decoy search.';}
elsif ($decoyMethod eq 'INT:BANK') {$fdrString.='Mixed databank with valid and decoy sequences.';}
else { # External decoy search
	my ($decoyFile)=($decoyMethod=~/EXT:(.+)/);
	$fdrString.="External search file ($decoyFile).";
}
my @importInfos = ($maxRankString,$minScoreString,$fdrString);

push @importActions, {title => 'Import parameters', mainInfos => \@importInfos};

## History in DB ##
foreach my $id (keys %unorderedHist){
	my @mainInfos = ();
	my @additionalInfos = ();

	# PhosphoRS #
	if($unorderedHist{$id}{'valType'} eq 'prs'){
		my ($massTolerance) = ($unorderedHist{$id}{'paramStrg'} =~ /massTolerance:([^;]+)/);
		my ($activationType) = ($unorderedHist{$id}{'paramStrg'} =~ /activationType:([^;]+)/);
		my ($threshold) = ($unorderedHist{$id}{'paramStrg'} =~ /threshold:([^;]+)/);

		@mainInfos = ("<B>Threshold</B>: $threshold%",
			      "<B>Activation type</B>: $activationType",
			      "<B>Mass tolerance</B>: $massTolerance");

		my ($total) = ($unorderedHist{$id}{'queryValStrg'} =~ /total:(\d+)/);
		my ($confirmed) = ($unorderedHist{$id}{'queryValStrg'} =~ /confirmed:(\d+)/);
		my ($changed) = ($unorderedHist{$id}{'queryValStrg'} =~ /changed:(\d+)/);
		my ($flagged) = ($unorderedHist{$id}{'queryValStrg'} =~ /flagged:(\d+)/);
		my ($unchanged) = ($unorderedHist{$id}{'queryValStrg'} =~ /unchanged:(\d+)/);

		@additionalInfos = ("<B>Phosphopeptides</B>:",
				    "$confirmed/$total confirmed",
				    "$changed/$total changed",
				    "$flagged/$total different but < threshold",
				    "$unchanged/$total confirmed but < threshold");
	}

	push @others, {title => $valTypeTitle{$unorderedHist{$id}{valType}},
		       user => $unorderedHist{$id}{validUser},
		       date => $unorderedHist{$id}{date},
		       mainInfos => \@mainInfos,
		       additionalInfos => \@additionalInfos};
}

## Query VarMods ##
if ($countQueryVMod && $view eq 'mass') {
	my $mainInfo = getVModMainInfo($countQueryVMod);
	push @others, {title => $valTypeTitle{'varMod'}, mainInfos => [$mainInfo]};
}

# Peptide VarMods #
if ($countPeptideVMod && $view eq 'bio') {
	my $mainInfo = getVModMainInfo($countPeptideVMod);
	push @others, {title => $valTypeTitle{'varMod'}, mainInfos => [$mainInfo]};
}


my ($light,$dark)=&promsConfig::getRowColors;
if ($format eq 'xls') {
	generateSpreadsheet();
}

############
###>HTML<###
############
print header(-charset=>'utf-8');
warningsToBrowser(1);
my $bgColor = $light;
print qq|
<HTML>
<HEAD>
<TITLE>$analysisName Validation History</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css"></HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Analysis <FONT color="red">$analysisName</FONT></FONT>|;
if ($projectFullAccess){
	my $selBioStrg=($view eq 'bio')? ' selected' : '';
	print qq
| <SELECT class="title" onchange="window.location='./displayValidationHistory.cgi?ID=$analysisID&view='+this.value">
	<OPTION value="mass">On-going</OPTION>
	<OPTION value="bio"$selBioStrg>Reported</OPTION>
</SELECT>
|;
}
print qq
|<FONT class="title"> Validation History</FONT>
<BR><BR>
<INPUT type="button" value="Export" class="title3" onclick="window.location='$promsPath{cgi}/displayValidationHistory.cgi?ID=$analysisID&format=xls';">
<BR><BR>
|;

###################################
# Generating Import-Actions Table #
###################################
if (scalar @importActions) {
	print generateTable(\@importActions);
	print "<BR><BR>\n";
}


#################################
# Generating Step History Table #
#################################
if (scalar @stepHistory) {
	print generateTable(\@stepHistory);
	print "<BR><BR>\n";
}


###############################
# Generating Other HTML table #
###############################
if (scalar @others) {
	print qq
|<FONT class="title">Other actions</FONT>
<BR><BR>
|;
	print generateTable(\@others);
}


print "<FONT class=\"title2\">No history recorded.</FONT>" unless (scalar @importActions) or (scalar @stepHistory) or (scalar @others);  # can be visible by massists on biologist view
print qq
|<BR><INPUT type="button" value=" Close Window " onclick="window.close()"/>
</CENTER>
</BODY>
</HTML>
|;

sub getVModMainInfo{
	my $count = shift;

	my $mainInfo = "$count peptide";
	if ($count > 1){
		$mainInfo .= 's have';
	} else {
		$mainInfo .= ' has';
	}
	$mainInfo .= " been modified.";

	return $mainInfo;
}

sub generateTable{
	my $refData = shift;

	my $html;

	$html .= "<TABLE align=\"center\" border=0 cellspacing=0 width=\"750px\">\n";
	foreach my $rowData (@$refData){
		$bgColor = ($bgColor eq $dark)?$light:$dark;
		$html .= generateRow(%{$rowData}, bgColor => $bgColor);
	}
	$html .= "</TABLE>\n";

	return $html;
}


sub generateRow{
	my (%data) = @_;

	my $bgColor = $data{bgColor};
	my $step = (defined($data{step}))? $data{step} : 0;
	my $title = $data{title};
	my $user = $data{user};
	my $date = $data{date};
	my @mainInfos = ($data{mainInfos})? @{$data{mainInfos}}:();
	my @additionalInfos = ($data{additionalInfos})? @{$data{additionalInfos}}:();

	my $titleString = ($step)?"$step. $title":$title;
	my $titleColspan = ($user && $date)?2:3;

	#--- Header ---#
	my $html = qq
|<TR bgcolor="$bgColor"><TD colspan=$titleColspan valign=top>&nbsp;<FONT class="title2">$titleString</FONT>&nbsp;&nbsp;</TD>
|;
	if ($user && $date) {
		$html .= qq
|<TD><B>User:</B> $user&nbsp;&nbsp;&nbsp;<B>Date:</B> $date&nbsp;</TD>
|;
	}
	$html .= "</TR>\n";

	#--- Body ---#
	# Left cell #
	$html .= "<TR bgcolor=\"$bgColor\"><TD width=30></TD><TD>";
	$html .= join "<BR>\n", @mainInfos;
	$html .= "</TD>";
	# Right cell #
	$html .= "<TD valign=\"top\" align=\"right\">";
	if (scalar @additionalInfos) {
		$html .= "<TABLE border=0 cellpadding=0 width=220><TD nowrap>";
		$html .= join "<BR>\n", @additionalInfos;
		$html .= "</TD></TABLE>";
	}
	$html .= "</TD></TR>\n";

	return $html;
}

sub generateSpreadsheet{

	print header(-type=>"application/vnd.ms-excel",-attachment=>"$analysisName - Validation History.xls");
	my $workbook = Spreadsheet::WriteExcel->new("-");
	my $worksheet = $workbook->add_worksheet();

	#--- Formatting ---#
	my %format;
	my $darkColor = $workbook->set_custom_color(40, $dark);
	my $lightColor = $workbook->set_custom_color(41, $light);
	# Main title that contains analysis name
	$format{mainHeader} = $workbook->add_format(align => 'center', bg_color => $darkColor);
	$format{mainHeader}->set_bold();
	# Step title
	$format{header} = $workbook->add_format(align => 'left', bg_color => $lightColor);
	$format{header}->set_bold();
	# User and date column headers
	$format{subHeader} = $workbook->add_format(align => 'center', bg_color => $lightColor);
	$format{subHeader}->set_bold();
	# Cells that contain user name and date
	$format{vTop} = $workbook->add_format(align => 'center', valign => 'top');
	# Generic centered cells
	$format{center} = $workbook->add_format(align => 'center');

	my $spacerHeight = 2;

	my $NCOL = 4;
	$worksheet->merge_range(0,0,0,$NCOL-1,"$analysisName Validation History",$format{mainHeader});
	$worksheet->set_column(0,0,70);
	$worksheet->set_column(1,1,40);
	$worksheet->set_column(2,2,20);
	$worksheet->set_column(3,3,20);

	my $r=1;
	# Step history
	foreach my $stepData (@importActions,undef,@stepHistory,undef,@others){ # undef is used to make a spacer

		if (defined $stepData) {
			my $mergeStepTitleRow = 1;
			$mergeStepTitleRow++ unless $stepData->{user};
			$mergeStepTitleRow++ unless $stepData->{date};
			my $userCol = ($stepData->{date})?2:3;
			my $dateCol = 3;

			my $stepTitle = $stepData->{title};
			$stepTitle = $stepData->{step}.'. '.$stepTitle if $stepData->{step};
			$worksheet->merge_range($r,0,$r,$mergeStepTitleRow,$stepTitle,$format{header});
			$worksheet->write_string($r,$userCol,'User',$format{subHeader}) if $stepData->{user};
			$worksheet->write_string($r,$dateCol,'Date',$format{subHeader}) if $stepData->{date};
			$r++;

			my @mainInfos = ($stepData->{mainInfos})? map {cleanString($_)} @{$stepData->{mainInfos}}:();
			my @additionalInfos = ($stepData->{additionalInfos})? map {cleanString($_)} @{$stepData->{additionalInfos}}:();

			my $nrow = max(scalar @mainInfos,scalar @additionalInfos,1);

			my $mainRowSize = $mergeStepTitleRow;
			$mainRowSize++ unless scalar @additionalInfos;
			if (scalar @mainInfos) {
				for(my $i=0;$i<=$#mainInfos;$i++){
					if ($mainRowSize == 1) {
						$worksheet->write_string($r+$i,0,$mainInfos[$i]);
					}
					else {
						$worksheet->merge_range($r+$i,0,$r+$i,$mainRowSize-1,$mainInfos[$i],$workbook->add_format());
					}
				}
			}
			if (scalar @additionalInfos) {
				$worksheet->write_col($r,1,\@additionalInfos);
			}
			if ($stepData->{user} || $stepData->{date}) {
				if ($nrow > 1) {
					$worksheet->merge_range($r,$userCol,$r+$nrow-1,2,$stepData->{user},$format{vTop}) if $stepData->{user};
					$worksheet->merge_range($r,$dateCol,$r+$nrow-1,3,$stepData->{date},$format{vTop}) if $stepData->{date};
				}
				else {
					$worksheet->write_string($r,$userCol,$stepData->{user},$format{center}) if $stepData->{user};
					$worksheet->write_string($r,$dateCol,$stepData->{date},$format{center}) if $stepData->{date};
				}
			}
			$r += $nrow;
		}
		else {
			# spacer
			foreach (1..$spacerHeight){
				$worksheet->merge_range($r+$_-1,0,$r+$_-1,$NCOL-1,"",$workbook->add_format());
			}
			$r+=$spacerHeight;
		}
	}

	$workbook->close();

	exit;
}

sub cleanString{
	# Remove all HTML elements from a string in order to display it properly in XLS format
	my $string = shift;

	$string =~ s/<IMG[^>]+Yellow[^>]+>/up/;
	$string =~ s/<IMG[^>]+Orange[^>]+>/down/;
	$string =~ s/<\/?[^>]+>//g;
	$string =~ s/&#8804/<=/g;
	$string =~ s/&#8805/>=/g;
	$string =~ s/&nbsp;?/ /g;
	$string =~ s/\s+\Z//;

	return $string;
}


####>Revision history<####
# 1.1.6 [ENHANCEMENT] Added displaying of best query per ion selection parameter (VS 30/03/20)
# 1.1.5 Compatible with "##PRB_MQ=xxx" in PEPTIDE_MODIFICATION.REF_POS_STRING for MaxQuant (PP 14/02/17)
# 1.1.4 Minor bug correction in Manual selection display and -charset=>'utf-8' (GA 02/09/14)
# 1.1.3 Handles precomputed FDR by Percolator in PD (PP 14/05/14)
# 1.1.2 Handles FDR algo info (PP 14/04/14)
# 1.1.1 Change syntax $x=$y // $z to old ()? one (PP 10/03/14)
# 1.1.0 New cases,code improvements and XLS export (FY 19/09/13)
# 1.0.7 Add phosphoRS analysis (FY 30/03/12)
# 1.0.6 Sorting correction & massist display correction (FY 17/11/11)
# 1.0.5 Manager has full access validation data (PP 01/10/11)
# 1.0.4 SQ ranks management (FY 01/06/11)
# 1.0.3 Minor bug & add untraced selection category (FY 26/04/11)
# 1.0.2 Minor bug (FY 18/04/11)
# 1.0.1 Minor changes in display & code formatting (PP 12/04/11)
# 1.0.0 New script to display validation history of an analysis (FY 16/03/2011)
