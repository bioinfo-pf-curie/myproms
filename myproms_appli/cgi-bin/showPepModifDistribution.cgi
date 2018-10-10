#!/usr/local/bin/perl -w

################################################################################
# showPepModifDistribution.cgi          1.3.1                                  #
# Authors: P. Poullet (Institut Curie)                                         #
# Contact: myproms@curie.fr                                                    #
# Displays peptide modification(s) distribution across a set of Analyses       #
# Called by selectAnalyses.cgi                                                 #
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

#print header;warningsToBrowser(1);
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my %itemIcones=&promsConfig::getItemIcones;
my %terminalModifCodes=('-'=>'Prot. N-term','='=>'N-term','+'=>'Prot. C-term','*'=>'C-term');

###################
####>Arguments<####
###################
my $branchID=lc(param('ID')); # just to be safe
my ($item,$itemID)=split(':',$branchID);
my @selAnalyses=param('anaList');
my $virtualPep=param('virtualPep');
my $quantifPep=param('quantifPep');
my $residuePep=param('residuePep');
my $setIntersection=($residuePep && param('setInter'))? 1 : 0;
my $modFilterType=param('modFilterType') || '';
my $modFilterID=($modFilterType && param('modFilter'))? param('modFilter') : 0;
my $sitesFromRestrict=($modFilterType eq 'restrict' && $modFilterID && param('sitesFromRestrict'))? 1 : 0;
my $protMulti=param('protMulti') || 0; # a protein can be with and without sel modif(s) at the same time
my $totAnalyses=scalar @selAnalyses;
my (%selModifications,%restrictResidues,%restrictTermini,%restrictTerminiMod,%dispRestricRes,$dispRestrictResStrg);
my $usedBoxPos=0;
foreach my $i (1..param('numSelBoxes')) {
	if (param("modifID_$i")) {
		@{$selModifications{++$usedBoxPos}}=param("modifID_$i");
		if ($residuePep) {
			foreach my $resData (param("residues_$i")) {
				my ($res,@modIDs)=split(':',$resData);
				#my ($res,@modIDs)=split(':',$resData{$i});
				if ($res=~/[-+]/) { # Protein N/C-term (no need to filter if peptide termini =*: true for all peptides)
					$restrictTermini{$res}=1;
					foreach my $modID (@modIDs) {$restrictTerminiMod{$modID}=$res;}
				}
				elsif ($res=~/\w/) { # && $setIntersection
					$restrictResidues{$res}=1;
				}
				$dispRestricRes{$usedBoxPos}{$res}=$terminalModifCodes{$res} || $res;
			}
			if ($setIntersection) {
				$dispRestrictResStrg.=' & ' if $dispRestrictResStrg;
				$dispRestrictResStrg.= '['.join(',',sort values %{$dispRestricRes{$usedBoxPos}}).']';
			}
		}
	}
}

my $restrictResStrg=join('',keys %restrictResidues);
#print "*RESTRICT_RES=$restrictResStrg<BR>\n";


#######################
####>Starting HTML<####
#######################
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my $bgColor=$lightColor;

print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);
#<SCRIPT src="$promsPath{html}/js/Raphael/g.raphael-min.js"></SCRIPT>
#<SCRIPT src="$promsPath{html}/js/Raphael/g.pie-min.js"></SCRIPT>
print qq
|<HTML>
<HEAD>
<TITLE>PTM Distribution</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE>
.barCell {
	background-color: #F3F3F3;
	border-style: solid;
	border-width: 1px;
	border-color: #000;
	width: 400px;
}
</STYLE>
<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/barPlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/vennDiagram.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/genericPlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
<SCRIPT type="text/javascript">
function changeSequenceCountType(countType) {
	var visC='',visP='none';
	if (countType=='percent') {
		visC='none',visP='';
	}
	document.getElementById('vennSeqDIV_C').style.display=visC;
	document.getElementById('vennSeqDIV_P').style.display=visP;
}
</SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
|;

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

####>Fetching Analyses info<####
my @sthItem;
if ($item eq 'experiment') {
	$sthItem[0]=$dbh->prepare("SELECT A.NAME,A.VALID_STATUS,CONCAT(G.NAME,' > ',SP.NAME) FROM ANALYSIS A,SAMPLE S,SPOT SP,GEL2D G WHERE A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_SPOT=SP.ID_SPOT AND SP.ID_GEL2D=G.ID_GEL2D AND A.ID_ANALYSIS=?");
	$sthItem[1]=$dbh->prepare("SELECT A.NAME,A.VALID_STATUS,S.NAME FROM ANALYSIS A,SAMPLE S WHERE A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_SPOT IS NULL AND A.ID_ANALYSIS=?");
}
elsif ($item eq 'gel2d') {
	$sthItem[0]=$dbh->prepare("SELECT A.NAME,A.VALID_STATUS,SP.NAME FROM ANALYSIS A,SAMPLE S,SPOT SP WHERE A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_SPOT=SP.ID_SPOT AND A.ID_ANALYSIS=?");
}
elsif ($item=~/sample|spot|analysis/) {
	$sthItem[0]=$dbh->prepare("SELECT NAME,VALID_STATUS FROM ANALYSIS WHERE ID_ANALYSIS=?");
}
print qq
|<INPUT type="button" value=" <<< Back " onclick="window.location='$promsPath{cgi}/selectAnalyses.cgi?ID=$branchID&callType=peptidePTM'"/>
<CENTER>
|;
my %analysisName;
if ($item eq 'analysis') {
	$sthItem[0]->execute($selAnalyses[0]);
	my ($anaName,$valStat)=$sthItem[0]->fetchrow_array;
	print qq
|<FONT class="title" style="text-align:center">Peptide Modifications Distribution of Analysis <FONT color="#DD0000">$anaName</FONT></FONT>
<BR><BR>
|;
	$sthItem[0]->finish;
}
else {
	my $anaStrg=($totAnalyses > 1)? 'Analyses' : 'Analysis';
	my $numCols=int($totAnalyses/5);
	$numCols++ if $totAnalyses % 5;
	print qq
|<FONT class="title" style="text-align:center">Peptide Modifications Distribution</FONT>
<BR><BR>
<TABLE bgcolor="$darkColor">
<TR><TH class="title2" colspan="$numCols">&nbsp;Selected $anaStrg&nbsp;</TH></TR>
<TR bgcolor="$lightColor"><TH align="left" valign=\"top\">
|;
	my $count=0; my $totCount=0;
	foreach my $anaID (@selAnalyses) {
		foreach my $sth (@sthItem) {
			$sth->execute($anaID);
			my ($anaName,$valStat,$parentsName)=$sth->fetchrow_array;
			if ($anaName) {
				my $anaCode=($valStat==-1)? 'analysis:no_scan' : ($valStat==0)? 'analysis:no_val' : ($valStat==1)? 'analysis:part_val' : 'analysis:val';
				print '&bull;';
				print $parentsName.' > ' if $parentsName;
				print "<IMG src=\"$promsPath{images}/$itemIcones{$anaCode}\">&nbsp;$anaName&nbsp;<BR>";
				$analysisName{$anaID}=$anaName;
				last;
			}
		}
		$count++;
		$totCount++;
		if ($count==5) {
			print "</TH>";
			print "<TH align=\"left\" valign=\"top\">" if $totCount < $totAnalyses;
			$count=0;
		}
	}
	print "</TH>" if $totAnalyses % 5;
	print "</TR></TABLE><BR>\n";
	foreach my $sth (@sthItem) {$sth->finish;}
}

my $sthMod=$dbh->prepare("SELECT PSI_MS_NAME,INTERIM_NAME,DES,SYNONYMES FROM MODIFICATION WHERE ID_MODIFICATION=?");
my $filterModName;
print "<TABLE bgcolor=\"$darkColor\"><TR><TH valign=\"top\" class=\"title2\" nowrap>&nbsp;Parameters :</TH><TH bgcolor=\"$lightColor\" align=\"left\">";
if ($virtualPep || $quantifPep || $residuePep) {
	print "-Peptides identified from XIC extractions are included.&nbsp;<BR>\n" if $virtualPep;
	print "-Only peptides with quantification data are considered.&nbsp;<BR>\n" if $quantifPep;
	if ($dispRestrictResStrg) {
		print "-Only peptides modifiable on $dispRestrictResStrg are considered.&nbsp;<BR>\n" if $dispRestrictResStrg=~/\w/; # no restriction for -,=,+,*
	}
	elsif ($residuePep) {
		print "-All peptides modifiable on any specified residues are considered.&nbsp;<BR>\n";
	}
	if ($modFilterID) {
		$sthMod->execute($modFilterID);
		my ($psi,$interim,$des,$syn)=$sthMod->fetchrow_array;
		$filterModName=($psi)? $psi : ($interim)? $interim : ($des)? $des : $syn;
		my $withStrg=($modFilterType eq 'restrict')? 'with' : 'without';
		print "-Only peptides $withStrg <FONT color=\"#DD0000\">$filterModName</FONT> are considered.&nbsp;<BR>\n";
		print "-Modification sites counted are <FONT color=\"#DD0000\">$filterModName</FONT>-sites.&nbsp;<BR>\n" if $sitesFromRestrict;
	}
	if ($protMulti) {
		print "-A protein can be \"with\" and \"without\" selected modification(s) at the same time.&nbsp;<BR>\n";
	}
	else {
		print "-A protein cannot be \"with\" and \"without\" selected modification(s) at the same time.&nbsp;<BR>\n";
	}
}
else {print 'None';}
print qq
|</TH></TR></TABLE><BR>
<DIV id="waitDIV">
<BR><BR><BR><BR><BR><FONT class="title3">Fetching data...</FONT><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">
<BR><BR><BR><FONT class="title3">Status:<SPAN id="waitSPAN">...<SPAN></FONT>
</DIV>
|;

####>Fetching modification info<####
my (%modifications,%modifGroup);
#my $labelModUsed=0;
$modifications{0}{0}=($residuePep)? 'W/o selected modified residues' : 'W/o selected modifications'; # Control
$modifGroup{0}='A';
foreach my $setPos (keys %selModifications) {
	foreach my $modID (@{$selModifications{$setPos}}) {
		$sthMod->execute($modID);
		my ($psi,$interim,$des,$syn)=$sthMod->fetchrow_array;
		if ($syn) {
			$syn=~s/^##//;
			$syn=~s/##\Z//;
			$syn=join('/',(split('##',$syn)));
		}
		else {$syn='modification #'.$modID;}
		$modifications{$setPos}{$modID}=($psi)? $psi : ($interim)? $interim : ($des)? $des : $syn;
		if ($residuePep) {
			$modifications{$setPos}{$modID}.='['.join('/',sort values %{$dispRestricRes{$setPos}}).']';
		}
		#$labelModUsed=1 if $modifications{$setPos}{$modID}=~/^Label:/; # SILAC
	}
	$modifGroup{$setPos}=chr($setPos + 65); # B,C,D,...
}
$sthMod->finish;

####>Fetching all peptides<####
my $virtualPepQuery=($virtualPep)? '' : 'AND P.VALID_STATUS > 0'; #'AND QUERY_NUM IS NOT NULL';
#my $quantifPepQuery=($quantifPep)? 'INNER JOIN PEPTIDE_QUANTIFICATION PQ ON P.ID_PEPTIDE=PQ.ID_PEPTIDE' : '';
my ($flankingResQuery,$protTermQuery)=('','');
if ($residuePep) {
	if ($setIntersection) {
		my $protTermPattern=($restrictTermini{'-'} && $restrictTermini{'+'})? '%-%' : ($restrictTermini{'-'})? '-%' : ($restrictTermini{'+'})? '%-' : '';
		$protTermQuery="AND PPA.FLANKING_AA LIKE '$protTermPattern'" if $protTermPattern;
	}
	if (scalar keys %restrictTerminiMod) {
		$flankingResQuery=",GROUP_CONCAT(DISTINCT PPA.FLANKING_AA SEPARATOR ':')";
	}
}
my $sthPep=$dbh->prepare("SELECT P.ID_PEPTIDE,PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING SEPARATOR '&'),CHARGE,GROUP_CONCAT(DISTINCT AP.ID_PROTEIN,':',ABS(PEP_BEG) SEPARATOR '&') $flankingResQuery
							FROM PEPTIDE P
							LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
							INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE $protTermQuery
							INNER JOIN ANALYSIS_PROTEIN AP ON PPA.ID_PROTEIN=AP.ID_PROTEIN AND VISIBILITY >= 1 AND P.ID_ANALYSIS=AP.ID_ANALYSIS
							WHERE P.ID_ANALYSIS=? $virtualPepQuery GROUP BY P.ID_PEPTIDE"); #$quantifPepQuery

my ($quantifParentDir,%anaPepQuantif,%quantifFileRead,%pepIsQuantified);
if ($quantifPep) { # restrict to quantified peptides
	my $projectID=&promsMod::getProjectID($dbh,$itemID,$item);
	$quantifParentDir="$promsPath{quantification}/project_$projectID";
	my $sthAnaQuant=$dbh->prepare("SELECT Q.ID_QUANTIFICATION FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND FOCUS='peptide' AND ID_ANALYSIS=?");
	foreach my $anaID (@selAnalyses) {
		$sthAnaQuant->execute($anaID);
		while (my ($pepQuantifID)=$sthAnaQuant->fetchrow_array) {
			push @{$anaPepQuantif{$anaID}},$pepQuantifID;
		}
	}
	$sthAnaQuant->finish;
}
my (%modifPeptides,%anaModifPeptides,%modifChargePeptides,%groupDistribution,%groupSize,%vennGroups,%distinctProteins,%distinctSites,%allProteins,%unlocalizedPep);
my $curAnalysis=0;
foreach my $anaID (@selAnalyses) {
	$curAnalysis++;
	print qq
|<SCRIPT type="text/javascript">document.getElementById('waitSPAN').innerHTML=' Processing Analysis $curAnalysis/$totAnalyses...';</SCRIPT>
|;

	##>Fetch list of quantified peptides from quantification file(s)
	if ($anaPepQuantif{$anaID}) {
		foreach my $pepQuantifID (@{$anaPepQuantif{$anaID}}) {
			next if $quantifFileRead{$pepQuantifID};
			foreach my $quantFile (glob ("$quantifParentDir/quanti_$pepQuantifID/peptide_quantification*")) {
				open (QUANT,$quantFile);
				while (<QUANT>) {
					next if $.==1;
					my $pepID=(split(/\t/,$_))[1];
					$pepIsQuantified{$pepID}=1;
				}
				close QUANT;
			}
			$quantifFileRead{$pepQuantifID}=1;
		}
	}

	$sthPep->execute($anaID);
	while (my ($pepID,$pepSeq,$modCode,$charge,$protMatchStrg,$flankingAAStrg)=$sthPep->fetchrow_array) { # $flankingAAStrg not always called!!!
		next if ($restrictResStrg && $pepSeq !~ /[$restrictResStrg]/); # Needed to clean up group {0} only if $setIntersection
		if ($quantifPep) {
			if (!$pepIsQuantified{$pepID}) {next;}
			else {delete $pepIsQuantified{$pepID};} # no longer needed
		}
		$flankingAAStrg=($flankingAAStrg)? ':'.$flankingAAStrg.':' : '';
		my %pepProtBeg;
		foreach my $protData (split(/&/,$protMatchStrg)) {
			my ($protID,$pepBeg)=split(/:/,$protData);
			$pepProtBeg{$protID}=$pepBeg; # multi-pos on same prot not handled! => counted as 1.
			$unlocalizedPep{$pepSeq}=1 unless $pepBeg;
		}
		my %matchedMods;
		if ($modCode) {
			$modCode='&'.$modCode;
			next if (($modFilterType eq 'exclude' && $modCode=~/&$modFilterID:/) || ($modFilterType eq 'restrict' && $modCode !~ /&$modFilterID:/));
			my %matchedSets; # only if setIntersection
			foreach my $setPos (keys %selModifications) {
				#>Checking if at least 1 modif of set is on selected residue(s)
				if ($residuePep) {
					my $okModPos=0;
					MODIF:foreach my $modID (@{$selModifications{$setPos}}) {
						#next if (!$setIntersection && $restrictTerminiMod{$modID} && $flankingAAStrg !~ /$restrictTerminiMod{$modID}/); # - or + for protein N/C-term
						#if (!$setIntersection && $restrictTerminiMod{$modID}) {
						#	next if (($restrictTerminiMod{$modID} eq '-' && $flankingAAStrg !~ /:-/) || ($restrictTerminiMod{$modID} eq '+' && $flankingAAStrg !~ /-:/)); # - or + for protein N/C-term
						#}
						if ($modCode=~/&$modID:([^&]+)/) { # /:$modID:/
							my @modPos=split(/\./,$1);
							next unless $dispRestricRes{$setPos};
							foreach my $pos (@modPos) {
								my $modRes;
								if ($pos=~/[=*]/) {$modRes=$pos;} # any peptide terminus
								elsif ($pos=~/[-+]/) { # Prot N/C-term
									next unless $restrictTerminiMod{$modID};
									next if (($restrictTerminiMod{$modID} eq '-' && $flankingAAStrg !~ /:-/) || ($restrictTerminiMod{$modID} eq '+' && $flankingAAStrg !~ /-:/)); # - or + for protein N/C-term
									$modRes=$pos;
								} # N/C-term modif
								else { # internal modif
									$pos--;
									($modRes)=($pepSeq=~/.{$pos}(.)/);
								}
								if ($dispRestricRes{$setPos}{$modRes}) { # MATCH!
									$okModPos=1;
									last MODIF;
								}
							}
						}
					}
					next unless $okModPos; # next setPos
					$matchedSets{$setPos}=1; # if $setIntersection;
				}
				else {$matchedSets{$setPos}=1;} # actual filtering occurs in bloc below
			}
			my $numMatchedSets=scalar keys %matchedSets;
			#next if ($setIntersection && $numMatchedSets < $usedBoxPos); # next peptide
#next ION if ($numMatchedSets < $usedBoxPos); # next peptide ion
			##>Encoding entities
			if (($setIntersection && $numMatchedSets == $usedBoxPos) || (!$setIntersection && $numMatchedSets >= 1)) {
				#foreach my $setPos (keys %selModifications) {#}
				foreach my $setPos (keys %matchedSets) {
					foreach my $modID (@{$selModifications{$setPos}}) {
						if ($modCode=~/&$modID:([^&]+)/) {
							$modCode=~/&$modFilterID:([^&]+)/ if $sitesFromRestrict; # overwrites with $modFilterID
							my @modPos=split(/\./,$1);
							foreach my $protID (keys %pepProtBeg) {
								#my $sitePosCode=$protID;
								foreach my $pos (@modPos) {
									next if $pos=~/;/; # skip with position ambiguity
									my $sitePosCode=($pos=~/[=-]/)? $protID.':N' : ($pos=~/[+*]/)? $protID.':C' : $protID.':'.($pepProtBeg{$protID}+$pos-1); # $protID.':'.$pepSeq.':'.$pos; #
									$distinctSites{$setPos}{$sitePosCode}=1; # Counts separately multiple sites on same pep!
								}
								#$distinctSites{$setPos}{$sitePosCode}=1; # Counts multiple sites on same pep as 1!
								$distinctProteins{$setPos}{$protID}=1;
								$allProteins{$protID}=1;
							}
							$matchedMods{$setPos}=1;
							#$modCode=~s/:$modID://; # remove mod info from mod list
							#$modCode=~s/&$modID:[^&]+//; # remove mod info from mod list
						}
					}
				}
			}
		}
		else { # no PTM
			next if $modFilterType eq 'restrict';
			if ($residuePep) { # check if pepSeq is modifiable
				my $okSeq=0;
				SET:foreach my $setPos (keys %dispRestricRes) {
					foreach my $modRes (keys %{$dispRestricRes{$setPos}}) {
						if ($pepSeq=~/$modRes/) {
							$okSeq=1;
							last SET;
						}
					}
				}
				next unless $okSeq;
			}
			$modCode='';
		}
		if (scalar keys %matchedMods) { # Modif set
			foreach my $setPos (keys %matchedMods) {
				$modifPeptides{$setPos}{$pepSeq.$modCode}=1;
				$modifChargePeptides{$setPos}{$pepSeq.$modCode.$charge}=1;
				$anaModifPeptides{$setPos}{$anaID}{$pepSeq.$modCode}=1;
				$groupDistribution{SEQ}{$pepSeq}{$modifGroup{$setPos}}=1;
				foreach my $protID (keys %pepProtBeg) {
					$groupDistribution{PROT}{$protID}{$modifGroup{$setPos}}=1;
				}
			}
		}
		else { # W/o modif set
			delete $unlocalizedPep{$pepSeq} unless $sitesFromRestrict; # bad localization is irrelevant for "w/o modif" group except if counting sites from restrict modif
			foreach my $protID (keys %pepProtBeg) {
				$allProteins{$protID}=1;
				$distinctProteins{0}{$protID}=1;
				$groupDistribution{PROT}{$protID}{'A'}=1;
			}
			next if ($residuePep && !$setIntersection && (($restrictTermini{'-'} && $flankingAAStrg !~ /:-/) || ($restrictTermini{'+'} && $flankingAAStrg !~ /-:/))); # - or + for protein N/C-term (already filtered if $setIntersection)
			$modifPeptides{0}{$pepSeq.$modCode}=1;
			$modifChargePeptides{0}{$pepSeq.$modCode.$charge}=1;
			$anaModifPeptides{0}{$anaID}{$pepSeq.$modCode}=1;
			$groupDistribution{SEQ}{$pepSeq}{'A'}=1; # $pepSeq.$modCode
			#foreach my $protID (keys %pepProtBeg) {
			#	$allProteins{$protID}=1;
			#	$distinctProteins{0}{$protID}=1;
			#	$groupDistribution{PROT}{$protID}{'A'}=1;
			#}
			if ($sitesFromRestrict && $modCode && $modCode=~/&$modFilterID:([^&]+)/) { # !!!!!! TEST !!!!!!!
				my @modPos=split(/\./,$1);
				foreach my $protID (keys %pepProtBeg) {
					foreach my $pos (@modPos) {
						next if $pos=~/;/; # skip with position ambiguity
						my $sitePosCode=($pos=~/[=-]/)? $protID.':N' : ($pos=~/[+*]/)? $protID.':C' : $protID.':'.($pepProtBeg{$protID}+$pos-1); # $protID.':'.$pepSeq.':'.$pos; #
						$distinctSites{0}{$sitePosCode}=1; # Counts separately multiple sites on same pep!
					}
				}
			}
		}
	}

}
$sthPep->finish;
$dbh->disconnect;

####>Venn diagram for sequences<####
foreach my $seq (keys %{$groupDistribution{SEQ}}) {
	my $group=join('',(sort{$a cmp $b} keys %{$groupDistribution{SEQ}{$seq}}));
	$vennGroups{SEQ}{$group}++;
	foreach my $modGr (keys %{$groupDistribution{SEQ}{$seq}}) {
		$groupSize{SEQ}{$modGr}+=$groupDistribution{SEQ}{$seq}{$modGr};
	}
}

####>Venn diagram for proteins & "w/o modif" group cleaning<####
foreach my $protID (keys %{$groupDistribution{PROT}}) {
	if (!$protMulti && $groupDistribution{PROT}{$protID}{'A'} && scalar keys %{$groupDistribution{PROT}{$protID}} > 1) { # clean protein distribution (remove from "w/o modif" if belongs to a modif group)
		delete $groupDistribution{PROT}{$protID}{'A'};																		# if $labelModUsed => assume label-channel comparison so same prot can be w & w/o modif
		delete $distinctProteins{0}{$protID};
	}
	my $group=join('',(sort{$a cmp $b} keys %{$groupDistribution{PROT}{$protID}}));
	$vennGroups{PROT}{$group}++;
}
####>Filling missing values (if any)<####
my $totalProteins=scalar keys %allProteins;
my ($totalPSM,$totalPep,$totalSites);
foreach my $setPos (keys %modifications) {
	%{$modifPeptides{$setPos}}=() unless $modifPeptides{$setPos};
	%{$modifChargePeptides{$setPos}}=() unless $modifChargePeptides{$setPos};
	foreach my $anaID (@selAnalyses) {
		%{$anaModifPeptides{$setPos}{$anaID}}=() if (!$anaModifPeptides{$setPos} || !$anaModifPeptides{$setPos}{$anaID});
	}
	$groupSize{SEQ}{$modifGroup{$setPos}}=0 if !$groupSize{SEQ}{$modifGroup{$setPos}};
	$totalPSM+=scalar keys %{$modifChargePeptides{$setPos}};
	$totalPep+=scalar keys %{$modifPeptides{$setPos}};
	$totalSites+=scalar keys %{$distinctSites{$setPos}};
}

my $groupLabelStrg='';
foreach my $setPos (sort{$a<=>$b} keys %modifications) {
	$groupLabelStrg.=',' if $groupLabelStrg;
	$groupLabelStrg.=$modifGroup{$setPos}.":'".join(' + ',sort{lc($a) cmp lc($b)} values %{$modifications{$setPos}})."'";
}
#------------------ Peptide ions HTML
print qq
|<SCRIPT type="text/javascript">
document.getElementById('waitDIV').style.display='none';
</SCRIPT>
<TABLE border=0 width=1000><TR>
	<TD valign="top" style="height:100%;width:100%"><TABLE bgcolor="$darkColor" style="height:100%;width:100%">
		<TR><TH class="title2" bgcolor="$darkColor">Peptide ions</TH></TR>
		<TR><TH class="barCell"><DIV id="barPsmDIV"></DIV></TH></TR>
		<TR><TH align="left" bgcolor="$lightColor" class="font11" nowrap valign=top style="height:100%;width:100%">
|;
my $numSets=scalar keys %modifications;
my (%legendsArray,%valuesArray);
foreach my $setPos (sort{$a<=>$b} keys %modifications) {
	my $numPSM=scalar keys %{$modifChargePeptides{$setPos}};
	my $pcPSMStrg=($numSets==2)? ' ('.(sprintf "%.1f",100*$numPSM/$totalPSM).' %)' : '';
	my $setName=join(' + ',sort{lc($a) cmp lc($b)} values %{$modifications{$setPos}});
	print "&bull;$setName: $numPSM$pcPSMStrg<BR>\n";
	$legendsArray{PSM}{$setPos}=$setName;
	$valuesArray{PSM}{$setPos}=$numPSM;
	$totalPSM+=$numPSM;
}
#------------------ Peptides (with modif)
print qq
|</TH></TR>
	</TABLE></TD>
	<TD valign="top" style="height:100%;width:100%"><TABLE bgcolor="$darkColor" style="height:100%;width:100%">
		<TR><TH class="title2" bgcolor="$darkColor">Peptides <FONT class="title3">(with modifications)</FONT></TH></TR>
		<TR><TH class="barCell"><DIV id="barPepDIV"></DIV></TH></TR>
		<TR><TH align="left" bgcolor="$lightColor" class="font11" nowrap valign=top style="height:100%;width:100%">
|;
foreach my $setPos (sort{$a<=>$b} keys %modifications) {
	my $numPep=scalar keys %{$modifPeptides{$setPos}};
	my $pcPepStrg=($numSets==2)? ' ('.(sprintf "%.1f",100*$numPep/$totalPep).' %)' : '';
	my $setName=join(' + ',sort{lc($a) cmp lc($b)} values %{$modifications{$setPos}});
	print "&bull;$setName: $numPep$pcPepStrg<BR>\n";
	$legendsArray{PEP}{$setPos}=$setName;
	$valuesArray{PEP}{$setPos}=$numPep;
}
#------------------ Venn sequences HTML
print qq
|</TH></TR>
	</TABLE></TD>
	<TD valign="top" style="height:100%;width:100%"><TABLE bgcolor="$darkColor" style="height:100%;width:100%">
		<TR><TH bgcolor="$darkColor"><FONT class="title2">Peptide sequences</FONT> <SELECT name="countType" class="title3" onchange="changeSequenceCountType(this.value)"><OPTION value="counts" selected>counts</OPTION><OPTION value="percent">%</OPTION></SELECT></TH></TR>
		<TR><TD><DIV id="vennSeqDIV_C"></DIV><DIV id="vennSeqDIV_P"></DIV></TD></TR>
		<TR><TH align="left" bgcolor="$lightColor" class="font11" nowrap valign=top style="height:100%;width:100%">
|;
foreach my $setPos (sort{$a<=>$b} keys %modifications) {
	print "&bull;",join(' + ',sort{lc($a) cmp lc($b)} values %{$modifications{$setPos}}),": $groupSize{SEQ}{$modifGroup{$setPos}}<BR>\n";
}
#------------------ Sites HTML
my $siteTitle=($sitesFromRestrict)? "$filterModName-sites" : 'Modification sites';
print qq
|</TH></TR>
	</TABLE></TD>
</TR><TR>
	<TD valign="top" style="height:100%;width:100%"><TABLE bgcolor="$darkColor" style="height:100%;width:100%">
		<TR><TH class="title2" bgcolor="$darkColor">$siteTitle</TH></TR>
		<TR><TH class="barCell"><DIV id="barSiteDIV"></DIV></TH></TR>
		<TR><TH align="left" bgcolor="$lightColor" class="font11" nowrap valign=top style="height:100%;width:100%">
|;
foreach my $setPos (sort{$a<=>$b} keys %modifications) {
	#next if $setPos == 0; # no sites for unmodified set
	my $setName=join(' + ',sort{lc($a) cmp lc($b)} values %{$modifications{$setPos}});
	my $numSites=scalar keys %{$distinctSites{$setPos}};
	my $pcSiteStrg=($numSets > 2)? ' ('.(sprintf "%.1f",100*$numSites/$totalSites).' %)' : '';
	print "&bull;$setName: $numSites$pcSiteStrg<BR>\n";
	$legendsArray{SITE}{$setPos}=$setName;
	$valuesArray{SITE}{$setPos}=$numSites;
}
#foreach my $pepSeq (sort keys %unlocalizedPep) {print "$pepSeq<BR>\n";}
my $numBadPep=scalar keys %unlocalizedPep;
if ($numBadPep) {
	print "<FONT color=\"#DD0000\">&bull;WARNING:<BR>&nbsp;$numBadPep peptide sequence(s) could not be correctly aligned.<BR>&nbsp;Modification sites number may be underestimated.</FONT>\n";
}
#------------------ Proteins HTML
print qq
|</TH></TR>
	</TABLE></TD>
	<TD valign="top" style="height:100%;width:100%"><TABLE bgcolor="$darkColor" style="height:100%;width:100%">
		<TR><TH class="title2" bgcolor="$darkColor">Proteins</TH></TR>
		<TR><TH class="barCell"><DIV id="barProtDIV"></DIV></TH></TR>
		<TR><TH align="left" bgcolor="$lightColor" class="font11" nowrap valign=top style="height:100%;width:100%">
|;
foreach my $setPos (sort{$a<=>$b} keys %modifications) {
	my $numProt=scalar keys %{$distinctProteins{$setPos}};
	my $pcProtStrg=' ('.(sprintf "%.1f",100*$numProt/$totalProteins).' %)';
	my $setName=join(' + ',sort{lc($a) cmp lc($b)} values %{$modifications{$setPos}});
	print "&bull;$setName: $numProt$pcProtStrg<BR>\n";
	$legendsArray{PROT}{$setPos}=$setName;
	$valuesArray{PROT}{$setPos}=$numProt;
}
#------------------ Venn proteins HTML
print qq
|</TH></TR>
	</TABLE></TD>
		<TD valign="top" style="height:100%;width:100%"><TABLE bgcolor="$darkColor" style="height:100%;width:100%">
		<TR><TH bgcolor="$darkColor"><FONT class="title2">Protein distribution</FONT></TH></TR>
		<TR><TD><DIV id="vennProtDIV"></DIV></TD></TR>
		<TR><TH align="left" bgcolor="$lightColor" class="font11" nowrap valign=top style="height:100%;width:100%"></TH></TR>
	</TABLE></TD>
</TR>
</TABLE>
|;
my %setColors=(0=>'#AA0000',1=>'#0000FF',2=>'#C0C000',3=>'#00FF00',4=>'#555');
if (scalar @selAnalyses > 1) {
	print qq
|<BR><BR><TABLE bgcolor="$darkColor"><TR><TH class="title2">Analysis-specific distribution
<DIV id="linePlotDIV" style="background-color:#FFFFFF;padding:2px;"></DIV>
</TH></TR></TABLE>
|;
}
print qq
|<SCRIPT type="text/javascript">
var psmData={};
|;
#------------------ Ions JS
foreach my $setPos (sort{$a<=>$b} keys %modifications) {
	print qq
|psmData['$setPos']={};
psmData['$setPos'].label = '$legendsArray{PSM}{$setPos}';
psmData['$setPos'].value = $valuesArray{PSM}{$setPos};
psmData['$setPos'].color = '$setColors{$setPos}';
|;
}
my $drawWidth=300;
my $padX=50;
my $barWidth=int(($drawWidth-2*$padX)*0.9/$numSets);
print qq
|var psmBarPlot = new barPlot({height:320,
								width:$drawWidth,
								paddingX:$padX,
								barWidth:$barWidth
								});
psmBarPlot.draw(psmData,'barPsmDIV');

var pepData={};
|;
#------------------ Peptides JS
foreach my $setPos (sort{$a<=>$b} keys %modifications) {
	print qq
|pepData['$setPos']={};
pepData['$setPos'].label = '$legendsArray{PEP}{$setPos}';
pepData['$setPos'].value = $valuesArray{PEP}{$setPos};
pepData['$setPos'].color = '$setColors{$setPos}';
|;
}
print qq
|var pepBarPlot = new barPlot({height:320,
								width:$drawWidth,
								paddingX:$padX,
								barWidth:$barWidth
								});
pepBarPlot.draw(pepData,'barPepDIV');

var siteData={};
|;
#------------------Sites JS
foreach my $setPos (sort{$a<=>$b} keys %modifications) {
	#next if $setPos==0;
	print qq
|siteData['$setPos']={};
siteData['$setPos'].label = '$legendsArray{SITE}{$setPos}';
siteData['$setPos'].value = $valuesArray{SITE}{$setPos};
siteData['$setPos'].color = '$setColors{$setPos}';
|;
}
print qq
|var siteBarPlot = new barPlot({height:320,
								width:$drawWidth,
								paddingX:$padX,
								barWidth:$barWidth
								});
siteBarPlot.draw(siteData,'barSiteDIV');

var protData={};
|; #------------------ Proteins
foreach my $setPos (sort{$a<=>$b} keys %modifications) {
	print qq
|protData['$setPos']={};
protData['$setPos'].label = '$legendsArray{PROT}{$setPos}';
protData['$setPos'].value = $valuesArray{PROT}{$setPos};
protData['$setPos'].color = '$setColors{$setPos}';
|;
}
print qq
|var protBarPlot = new barPlot({height:320,
								width:$drawWidth,
								paddingX:$padX,
								barWidth:$barWidth
								});
protBarPlot.draw(protData,'barProtDIV');
|;
#------------------ Venn sequences
my $totalCount=0;
foreach my $countType ('C','P') {
	my $divID='vennSeqDIV_'.$countType;
	print "\tnew VennDiagram({div:'$divID',size:330,legendPosition:'bottom',groupLabels:{$groupLabelStrg},setValues:{";
	my $firstGroup=1;
	foreach my $group (sort{length($a)<=>length($b) || $a cmp $b} keys %{$vennGroups{SEQ}}) {
		if ($firstGroup) {$firstGroup=0;} else {print ',';}
		if ($countType eq 'C') {
			print "$group:$vennGroups{SEQ}{$group}";
			$totalCount+=$vennGroups{SEQ}{$group};
		}
		else {
			my $pc=sprintf "%.2f",100*$vennGroups{SEQ}{$group}/$totalCount;
			$pc=~s/0+\Z//;
			print "$group:$pc";
		}
	}
	print "}});\n";
}
print qq
|	document.getElementById('vennSeqDIV_P').style.display='none'; // must stay visible at drawing

	new VennDiagram({div:'vennProtDIV',size:330,legendPosition:'bottom',groupLabels:{$groupLabelStrg},setValues:{
|;#------------------ Venn proteins
my $firstGroup=1;
foreach my $group (sort{length($a)<=>length($b) || $a cmp $b} keys %{$vennGroups{PROT}}) {
	if ($firstGroup) {$firstGroup=0;} else {print ',';}
	print "$group:$vennGroups{PROT}{$group}";
}
print "}});\n";

#------------ Analysis peptide count
if (scalar @selAnalyses > 1) {
	my $width=scalar @selAnalyses * 50;
	if ($width < 200) {$width=200;}
	elsif ($width > 1000) {$width=1000;}
	print qq
|	var GP=new genericPlot({
		div:'linePlotDIV',name:'GP',
		width:$width,height:300,
		axisX:{title:'Analysis position'},axisY:{title:'# Peptides'},
		zoomable:true,
		axisClosure:true
    });
|;
	my $setIdx=0;
	my %anaTotal;
	foreach my $setPos (sort{$modifGroup{$a} cmp $modifGroup{$b}} keys %modifications) {
		my $setName=join(' + ',sort{lc($a) cmp lc($b)} values %{$modifications{$setPos}});
		print "\tGP.addDataSet($setIdx,{name:'$setName',axisX:'X',axisY:'Y',color:'$setColors{$setPos}',line:{type:'line'}});\n";
		print "\tGP.addDataAsString($setIdx,'";
		# loop through all analyses
		my $anaPos=0;
		foreach my $anaID (@selAnalyses) {
			print ';' if $anaPos > 0;
			print "$analysisName{$anaID},$anaID,",++$anaPos,",",(scalar keys %{$anaModifPeptides{$setPos}{$anaID}});
			$anaTotal{$anaID}+=scalar keys %{$anaModifPeptides{$setPos}{$anaID}};
		}
		print "');\n";
		$setIdx++;
	}
	#Total
	print "\tGP.addDataSet($setIdx,{name:'Total',axisX:'X',axisY:'Y',color:'#000',line:{type:'line'}});\n";
	print "\tGP.addDataAsString($setIdx,'";
	my $anaPos=0;
	foreach my $anaID (@selAnalyses) {
		print ';' if $anaPos > 0;
		print "$analysisName{$anaID},$anaID,",++$anaPos,",$anaTotal{$anaID}";
	}
	print "');\n";
	print "\tGP.draw();\n";
}

print qq
|</SCRIPT>
</CENTER>
<BR><BR>
</BODY>
</HTML>
|;

####>Revision history<####
# 1.3.1 [Fix] Bug in an incomplete query $sthAnaQuant (PP 28/08/18)
# 1.3.0 Checks for quantified peptides in file(s) not DB (PP 11/05/18)
# 1.2.1 Protein can belong (or not) to both "with" and "w/o" category. Sites can be computed from the restriction modification (PP 23/03/17)
# 1.2.0 Added sites and proteins plots (PP 02/03/17)
# 1.1.1 Plot category renaming to remove ambiguity (PP 10/08/16)
# 1.1.0 Allows +/-setIntersection & bug fix for Protein N/C-term modifications (PP 21/04/16)
# 1.0.6 Comparison of same modif on different residues is now possible (PP 17/03/16)
# 1.0.5 Modification filter option for peptide modifications distribution & barplots (PP 18/01/16)
# 1.0.4 Improved query performance (PP 10/12/15)
# 1.0.3 Added forgotten body background image (PP 08/09/15)
# 1.0.2 Added progression status (PP 11/05/15)
# 1.0.1 Added modif sites and protein counts (PP 07/01/15)
# 1.0.0 Created by PP (PP 03/12/14)
