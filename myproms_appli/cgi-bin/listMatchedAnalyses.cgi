#!/usr/local/bin/perl -w

################################################################################
# listMatchedAnalyses.cgi      1.6.7                                           #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Lists query info of all analyses matched by a given peptide                  #
# during a validation process                                                  #
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

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

####################
####>Parameters<####
####################
&createFrames if param('START');
# print header;warningsToBrowser(1); # DEBUG
&processForm if param('save');
my $selectedItem=param('MATCH');
my $disableForm=(param('DISABFORM'))? 'disabled' : '';
my $sortType = (param('sortType'))? param('sortType') : 'score';
my $ChargSts = (param('ChargSts'))? param('ChargSts') : '';
my $matchInfo=param('INFO');
my ($color1,$color2)=&promsConfig::getRowColors;
my @selData=($selectedItem=~/seq_(\d+)_\d+_(\d+)/); # ($queryID,$queryNum,$rank)
my @matchInfo=($selData[0]."_$selData[1]"); # queryID_rank
push @matchInfo,split(/:/,param('INFO')); # list of matching peptides (validated: 1_pepID, non val: 0_qID_rNum) from other analyses


#########################################################################
####>Fetching data for query+peptide(s) sets in all matched analyses<####
#########################################################################
my (%anaStatus,%analysisID,%queryStatus,%queryCharge,%pepInfo,%score,%massExp,%massObs,%massDelta,%massDeltaPPM,%elutionTime);
my $sthValPep=$dbh->prepare("SELECT ID_ANALYSIS,QUERY_NUM,PEP_RANK,SCORE,MR_OBS,CHARGE,MR_EXP,MR_CALC,MR_DELTA,PEP_SEQ,ELUTION_TIME,VALID_STATUS FROM PEPTIDE WHERE ID_PEPTIDE=?");
my $refAna=1;
my ($sequence,$uniqCharge,$varMod);
foreach my $info (@matchInfo) {
	my @pepData=split('_',$info);
	##>Non validated Analyses
	if ($refAna || $pepData[0]<=1) { # refAna always <= 1
		my ($qID,$rNum);
		if ($refAna) {($qID,$rNum)=@pepData;} # qID_rank
		else {($qID,$rNum)=@pepData[1,2];} # valStat_qID_rank
		my $sthQuery=$dbh->prepare("SELECT ID_ANALYSIS,QUERY_NUM,VALID_STATUS,MASS_DATA,CHARGE,ELUTION_TIME,INFO_PEP$rNum FROM QUERY_VALIDATION WHERE ID_QUERY=$qID");
		$sthQuery->execute;
		my ($anaID,$qNum,$qStat,$mData,$charge,$eTime,$pInfo)=$sthQuery->fetchrow_array;
		my $rankID="$qID"."_$qNum"."_$rNum";
		if ($refAna) {
			($anaStatus{$rankID})=$dbh->selectrow_array("SELECT VALID_STATUS FROM ANALYSIS WHERE ID_ANALYSIS=$anaID");
			$refAna=0;
		}
		else {$anaStatus{$rankID}=$pepData[0];} # -1,0,1
		$analysisID{$rankID}=$anaID;
		$queryStatus{$rankID}=$qStat;
		$pepInfo{$rankID}=$pInfo;
		($score{$rankID})=($pInfo=~/SC=(-?\d+\.*\d*)/); # required for sorting list
		($massExp{$rankID})=($mData=~/EXP=(\d+\.*\d*)/);
		($massObs{$rankID})=($mData=~/OBS=(\d+\.*\d*)/);
		($massDelta{$rankID})=($pepInfo{$rankID}=~/DELT=(-?\d+\.*\d*)/);
		my ($massCalc)=($pepInfo{$rankID}=~/CALC=(\d+\.*\d*)/);
		$massDeltaPPM{$rankID}=sprintf "%.3f",1000000*($massDelta{$rankID}/$massCalc);
		$queryCharge{$rankID}=$charge;
		if ($eTime) {
			if ($eTime=~/et(\d+\.\d+);/) {
				$elutionTime{$rankID}=$1;#New format for elution time
			}
			else {
				if ($eTime=~/sp\d+;sc\d+;/) { # scan data: not handled
					$elutionTime{$rankID}='?';
				}
				else {
					$elutionTime{$rankID}=$eTime; # old format
				}
			}
		}
		$elutionTime{$rankID}='?' unless $elutionTime{$rankID};

		$sthQuery->finish;
		unless ($sequence) { # do only once
			($sequence)=($pInfo=~/SEQ=(\w+)/);
			#$charge=int(0.5 + $massExp{$rankID}/$massObs{$rankID});
			#($varMod)=($pInfo=~/VMOD=([^,]+)/);
			$varMod=&promsMod::toStringVariableModifications($dbh,'QUERY',$qID,$anaID,$sequence,$rNum);
		}
		$uniqCharge=$charge;
	}
	##>Validated Analyses
	else {
		my $pepID=$pepData[1];

		$sthValPep->execute($pepID);
		my ($anaID,$qNum,$rNum,$score,$mObs,$charge,$mExp,$mCalc,$mDelta,$seq,$eTime,$vStatus)=$sthValPep->fetchrow_array;
		my $rankID="$pepID"."_$qNum"."_$rNum";
		$vStatus=0 unless $vStatus;
		$anaStatus{$rankID}=2+$vStatus; # 2=unknown, 3=auto, 4=manual, 1 used by SEL for peptide in not-validated analysis (old:$pepData[0]; # 2)
		$analysisID{$rankID}=$anaID;
		$score{$rankID}=$score; # required for sorting list
		$massExp{$rankID}=$mExp;
		$massObs{$rankID}=$mObs;
		$massDelta{$rankID}=$mDelta;
		$massDeltaPPM{$rankID}=sprintf "%.3f",1000000*($mDelta/$mCalc);
		$queryCharge{$rankID}=$charge;
		if($eTime=~/et(\d+\.\d+);/){
			$elutionTime{$rankID}=$1;#New format for elution time
		}else{
			if($eTime=~/sp\d+;sc\d+;/){
				$elutionTime{$rankID}='?';
			}else{
				$elutionTime{$rankID}=($eTime)? $eTime : '?';#Before changing the format of elution time
			}
		}
		unless ($sequence) { # do only once
			$sequence=$seq;
			#$charge=int(0.5 + $mExp/$mObs);
			$varMod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$anaID,$seq);
		}
		$uniqCharge=$charge;
	}
}
$sthValPep->finish;
$varMod='' unless $varMod;

################################################
####>Fetching data for all matched analyses<####
################################################
my (%analysisName,%analysisParents,%analysisInfo);
my $selAna=$dbh->prepare("SELECT NAME,VALID_STATUS,DATA_FILE,WIFF_FILE,TAXONOMY,FILE_FORMAT FROM ANALYSIS WHERE ID_ANALYSIS=?");
foreach my $anaID (values %analysisID) {
	$selAna->execute($anaID);
	my ($validStatus,$dfile,$wfile,$tax, $fileFormat);
	($analysisName{$anaID},$validStatus,$dfile,$wfile,$tax, $fileFormat)=$selAna->fetchrow_array();
	$tax='Unknown' unless $tax;
	chomp($wfile) if $wfile; # \n at end before storeAnalysis.cgi 3.0.7
	my @itemInfo=&promsMod::getItemInfo($dbh,'ANALYSIS',$anaID);
	$analysisParents{$anaID}='';
	foreach my $i (1..$#itemInfo-1) {
		$analysisParents{$anaID}.="$itemInfo[$i]{NAME} > ";
	}
	#$analysisParents{$anaID}="$refItemInfo->[5] > $refItemInfo->[3] > ";
	$analysisInfo{$anaID}="<B>Search file:</B> $dfile<BR><B>MS file:</B> $wfile<BR><B>Taxonomy:</B> $tax";
	###>Moving data file to tmp directory (/tmp/(my)proms/validation already exists)
	#unless (-e "$promsPath{tmp}/validation/ana_$anaID" || $validStatus==2) {
	#	system "mkdir $promsPath{tmp}/validation/ana_$anaID";
	#	system "cp $promsPath{valid}/ana_$anaID/*.dat $promsPath{tmp}/validation/ana_$anaID/." if $fileFormat eq 'MASCOT.DAT';
	#	system "cp $promsPath{valid}/ana_$anaID/*.pgf $promsPath{tmp}/validation/ana_$anaID/." if $fileFormat eq 'PHENYX.XML';
	#}
}
$selAna->finish;
$dbh->disconnect;

#######################
####>Starting HTML<####
#######################
print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Matched Analyses</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TD {font-weight:bold;text-align:center;}
TD.left {text-align:left;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
| ;
&promsMod::popupInfo();
print "top.promsFrame.listMatchedAnalysesSort ='$sortType'; \n" if $sortType;
print qq
|function drawSpectrum(rankId,call) {
	selectObject(rankId);
	top.promsFrame.spectrumFrame.anaSpecFrame.location="$promsPath{cgi}/drawSpectrum.cgi?RID="+rankId+"&CALL="+call;
}
function selectObject(newId) {
	var oldObject=document.getElementById(selectedId);
	oldObject.style.color='#000000';
	var newObject=document.getElementById(newId);
//	if (newObject) {
		newObject.style.color='#DD0000'; //style.background = '#DD0000';
		selectedId=newId;
//	}
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
function activateForm() {
	if (activeForm) {return;}
	document.anaForm.save.style.fontWeight='bold';
	document.anaForm.save.style.color='#000000';
	activeForm=1;
}
function selectSort(type) {
	window.location="./listMatchedAnalyses.cgi?MATCH=$selectedItem&INFO=$matchInfo&DISABFORM=$disableForm&sortType="+type;
}
var selectedId='$selectedItem';
var activeForm=0;
</SCRIPT>
</HEAD>
<BODY topmargin=0 onload="drawSpectrum('$selectedItem','ana')" background="$promsPath{images}/bgProMS.gif">
<FORM name="anaForm" method="post">
<INPUT type="hidden" name="MATCH" value="$selectedItem" />
|;

my $rankIdList='';
foreach my $rankID (keys %score) {
	next if $anaStatus{$rankID}>=2; # validated analysis
	print "<INPUT type=\"hidden\" name=\"status_$rankID\" value=\"$queryStatus{$rankID}\"/>\n";
	print "<INPUT type=\"hidden\" id=\"pepInfo_$rankID\" name=\"pepInfo_$rankID\" value=\"$pepInfo{$rankID}\"/>\n";
	print "<INPUT type=\"hidden\" name=\"anaId_$rankID\" value=\"$analysisID{$rankID}\"/>\n";
	$rankIdList.="$rankID:";
}
print "<INPUT type=\"hidden\" name=\"rankIdList\" value=\"$rankIdList\"/>\n" if $rankIdList;

my %sortColor;
foreach my $sortMode ('score','delta','ppm','name','masse') {
	$sortColor{$sortMode}=($sortMode eq $sortType)? '#0000BB' : 'black';
}
my $chargeString = ($ChargSts eq 'all')? '' : "($uniqCharge+)";
my $massHeaderString = ($ChargSts ne 'all')? 'MR(obs)' : "<A href=\"javascript:selectSort('masse')\" style=\"font-style:italic;color:$sortColor{masse}\" onmouseover=\"popup('Click to sort analyses by <B>decreasing observed masse</B>.')\" onmouseout=\"popout()\">MR(obs)</A>";
print qq
|<CENTER>
<FONT class="title2">Data relative to <FONT color="#DD0000">$sequence$varMod $chargeString</FONT> in all MS/MS Analyses</FONT><BR>

<TABLE border=0 cellpadding=0 cellspacing=0>
<TR bgcolor="$color2">
	<TH class="rbBorder" width=40>Auto</TH>
	<TH class="rbBorder" width=30><IMG class="button" src="$promsPath{images}/good.gif" onclick="checkAllPeptides(document.anaForm,'select')" onmouseover="popup('Click to <B>select</B> all peptides.')" onmouseout="popout()"></TH>
	<TH class="rbBorder" width=30><IMG class="button" src="$promsPath{images}/bad.gif" onclick="checkAllPeptides(document.anaForm,'reject')" onmouseover="popup('Click to <B>reject</B> all peptides.')" onmouseout="popout()"></TH>
	<TH class="rbBorder" width=65>Query</TH>
	<TH class="rbBorder" width=50>Rank</TH>
	<TH class="rbBorder" width=55><A href="javascript:selectSort('score')" style="font-style:italic;color:$sortColor{score}" onmouseover="popup('Click to sort analyses by <B>descending score</B>.')" onmouseout="popout()">Score</A></TH>
	<TH class="rbBorder" width=75>$massHeaderString</TH>
	<TH class="rbBorder" width=70>Mr(exp)</TH>
	<TH class="rbBorder" width=50><A href="javascript:selectSort('delta')" style="font-style:italic;color:$sortColor{delta}" onmouseover="popup('Click to sort analyses by <B>ascending absolute mass delta (Da)</B>.')" onmouseout="popout()">Delta</A></TH>
	<TH class="rbBorder" width=50><A href="javascript:selectSort('ppm')" style="font-style:italic;color:$sortColor{ppm}" onmouseover="popup('Click to sort analyses by <B>ascending relative mass delta (ppm)</B>.')" onmouseout="popout()">ppm</A></TH>
	<TH class="rbBorder" width=60>Elution</TH>
	<TH class="bBorder" width=450 align=left>&nbsp&nbsp;<A href="javascript:selectSort('name')" style="font-style:italic;color:$sortColor{name}" onmouseover="popup('Click to sort analyses by <B>item name</B>.')" onmouseout="popout()">MS/MS Analysis</A></TH>
</TR>
|;
my $bgColor=$color1;
foreach my $rankID (sort{&sortAnalyses($sortType)} keys %analysisID) {
	print "<TR class=\"list\" bgcolor=\"$bgColor\">\n";

	##>Selec/Reject
	my ($autoString,$selectString,$rejectString);
	if ($anaStatus{$rankID}<=1) { # non validated
		my ($select)=($pepInfo{$rankID}=~/SEL=(-*\d)/); # cannot be -1
		my $autoImg=($select==1)? 'good.gif' : ($select==-2)? 'bad.gif' : 'blank.gif';
		my $checkSel=($select==2)? 'checked' : ''; # 1
		my $checkRej=($select==-3)? 'checked' : ''; # -2
		$autoString="<IMG src=\"$promsPath{images}/$autoImg\">";
		$selectString="<INPUT type=\"checkbox\" name=\"sel_$rankID\" $checkSel onClick=\"checkPeptide(document.anaForm.sel_$rankID,document.anaForm.rej_$rankID,$select)\" $disableForm/>";
		$rejectString="<INPUT type=\"checkbox\" name=\"rej_$rankID\" $checkRej onClick=\"checkPeptide(document.anaForm.rej_$rankID,document.anaForm.sel_$rankID,$select)\" $disableForm/>";
	}
	else { # validated
		if ($anaStatus{$rankID}==2) {$autoString=$selectString='<FONT color="#00A000">?</FONT>';} # valid method not defined
		elsif ($anaStatus{$rankID}==3) { # auto
			$autoString="<IMG src=\"$promsPath{images}/good.gif\">";
			$selectString='';
		}
		else { # 4 manual
			$autoString='';
			$selectString="<IMG src=\"$promsPath{images}/good.gif\">";
		}
		$rejectString='';
	}
	print "\t<TH>$autoString</TH>\n"; #
	print "\t<TH>$selectString</TH>\n"; # select
	print "\t<TH>$rejectString</TH>\n"; # reject

	##>Query + rank
	my ($queryNum,$rank)=($rankID=~/\d+_(\d+)_(\d+)/);
	my $imgString="&nbsp<IMG src=\"$promsPath{images}";
	if ($anaStatus{$rankID}<=1) {
		$imgString.=($queryStatus{$rankID}==-1)? '/lightGray1.gif' : ($queryStatus{$rankID}==0 || $queryStatus{$rankID}==-2)? '/lightRed1.gif' : '/lightGreen1.gif';
	}
	else {$imgString.='/lightGreen1.gif';}
	$imgString.='" hspace=0 border=0 height=11 width=11>&nbsp;';
	print "\t<TD class=\"left\">$imgString$queryNum</TD>\n";
	print "\t<TD>$rank</TD>\n";
# <A href=\"javascript:void(null)\" onmouseover=\"popup('')\" onmouseout=\"popout()\">
	##>Score
	printf "\t<TD><A href=\"javascript:void(null)\" onmouseover=\"popup('$score{$rankID}')\" onmouseout=\"popout()\">%.1f</A></TD>\n",$score{$rankID};

	##>Masses
	printf "\t<TD><A href=\"javascript:void(null)\" onmouseover=\"popup('<B>MR(obs):</B> $massObs{$rankID}<BR><B>Charge:</B> $queryCharge{$rankID}+')\" onmouseout=\"popout()\">%.2f</A></TD>\n",$massObs{$rankID};
	printf "\t<TD><A href=\"javascript:void(null)\" onmouseover=\"popup('$massExp{$rankID}')\" onmouseout=\"popout()\">%.2f</A></TD>\n",$massExp{$rankID};
	printf "\t<TD><A href=\"javascript:void(null)\" onmouseover=\"popup('<B>Delta:</B> $massDelta{$rankID} Da')\" onmouseout=\"popout()\">%.2f</A></TD>\n",$massDelta{$rankID};
	printf "\t<TD><A href=\"javascript:void(null)\" onmouseover=\"popup('<B>Delta:</B> $massDeltaPPM{$rankID} ppm')\" onmouseout=\"popout()\">%.1f</A></TD>\n",$massDeltaPPM{$rankID};

	##>Elution time
	my $eTime=($elutionTime{$rankID}=~/(.+) to /)? "$1+" : $elutionTime{$rankID};
	print "\t<TD><A href=\"javascript:void(null)\" onmouseover=\"popup('$elutionTime{$rankID} min')\" onmouseout=\"popout()\">$eTime</A></TD>\n";  #affichage partiels temps

	##>MS/MS Analysis
	my $anaStatImg=($anaStatus{$rankID}==-1)? 'lightGray1.gif' : ($anaStatus{$rankID}==0)? 'lightRed1.gif' : ($anaStatus{$rankID}==1)? 'lightYellow1.gif' : 'lightGreen1.gif';
	my $call=($anaStatus{$rankID}<=1)? 'ana' : 'val';
	print qq
|	<TD class="left"><IMG src="$promsPath{images}/$anaStatImg" hspace=0 border=0 height=11 width=11>&nbsp;$analysisParents{$analysisID{$rankID}}
	<A id="seq_$rankID" href="javascript:drawSpectrum('seq_$rankID','$call')" onmouseover="popup('$analysisInfo{$analysisID{$rankID}}')" onmouseout="popout()">$analysisName{$analysisID{$rankID}}</A>
</TD>
</TR>
|;
	$bgColor=($bgColor eq $color1)? $color2 : $color1;
}

print qq
|<TR><TH align=left colspan=9>
<INPUT type="submit" name="save" value="Save selection" style="width:125px;" $disableForm/>
</TH></TR>
</TABLE></CENTER></FORM>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup()</SCRIPT>
</HTML>
|;
exit;



##################################################################################################
##################################<<< Form Processing >>>#########################################
##################################################################################################
sub processForm {

	my $refRankID=param('MATCH');
	$refRankID=~s/seq_//;
	my $refAnalysisID=param("anaId_$refRankID");
	my @rankIdList=split(/:/,param('rankIdList'));
	my $maxRank=&promsConfig::getMaxRank; # <--- Maximum number of ranks retrieved (Absolute max = 10)
	my $newRefPeptide=0; # counter for updating number of validated peptides in summary window
	my $newRefQuery=0; # counter for updating number of verified queries in summary window
	my $newRefProteins=0; # counter for updating number of validated proteins in summary window
	my %allRankInfo; # contains info for all peptides already extracted from DB (prevents accessing DB more than once for same data)
	my %listIdentifiers; # list of identifiers/analysis affected

	####<Updating query & peptide data>####
	my $sthSEL1=$dbh->prepare("SELECT IDENTIFIER FROM RANK_PROTEIN_MATCH WHERE QUERY_NUM=? AND PEP_RANK=? AND ID_ANALYSIS=?");
	my %sthUpRank;
	my $infoString='';
	foreach my $rank (1..$maxRank) {
		$sthUpRank{$rank}=$dbh->prepare("UPDATE QUERY_VALIDATION SET INFO_PEP$rank=?,VALID_STATUS=? WHERE ID_QUERY=?");
		$infoString.="INFO_PEP$rank";
		$infoString.=',' if $rank<$maxRank;
	}
	my $sthSelInfo=$dbh->prepare("SELECT $infoString FROM QUERY_VALIDATION WHERE ID_QUERY=?");

	foreach my $rankID (@rankIdList) {
		my $select;
		my $updateQuery;
		my ($queryID,$queryNum,$rank)=split(/_/,$rankID);

		###<Checking if query & peptide need update>###
		my $analysisID=param("anaId_$rankID");
		my $queryStatus=param("status_$rankID");
		my $rankInfo=param("pepInfo_$rankID");

		$queryStatus=0 if $queryStatus==-1; # no previous validations (cannot be -2 => SEL=-1) 0=Starting value => final only if query is updated

		my ($oldSelect)=($rankInfo=~/SEL=(-*\d)/);  # cannot be -1 but can be -2
		my $computeStatus=0;
		if (param("sel_$rankID")) { # select box is checked
			$select=2;
			if ($oldSelect < 2) { # change in selection status (cannot be -1)
				$updateQuery=($oldSelect==1)? 1 : 2; # 1: update rank only, 2: update all
				#$queryStatus++;
				if ($queryStatus==-1 || $queryStatus==-2) {$queryStatus=1;}
				else {$queryStatus++;} # was >=0 before
				$newRefPeptide=1 if $rankID eq $refRankID;
			}
		}
		elsif (param("rej_$rankID")) { # reject box is checked
			$select=-3;
			if ($oldSelect > -3) { # change in selection status (cannot be -1)
				$updateQuery=($oldSelect==-2)? 1 : 2; # 1: update rank only, 2: update all
				if ($queryStatus==-1) {$queryStatus=-2;} # 1st reject
				elsif ($oldSelect >= 1) {
					$computeStatus=1;
					#$queryStatus--;
					$newRefPeptide=-1 if $rankID eq $refRankID;
				}
				else {$computeStatus=1;} # queryStatus was -2 (previous reject)
			}
		} # else => peptide is neither selected nor rejected (never-verified peptide)

		if ($computeStatus) { # Select must be performed to know new valid_status
			$sthSelInfo->execute($queryID);
			my @pepInfoList=$sthSelInfo->fetchrow_array;
			my %curPepInfo;
			my $i=0;
			foreach my $rk (1..$maxRank) {
				$curPepInfo{$rk}=$pepInfoList[$i];
				$i++;
			}
			my ($totRank,$selRank,$rejRank)=(0,0,0);
			foreach my $rk (1..$maxRank) {
				next unless $curPepInfo{$rk};
				my $sel;
				if ($rk==$rank) {$sel=$select;} # use new valid_status, not old one
				else {($sel)=($curPepInfo{$rk}=~/SEL=(-*\d)/);}
				$totRank++ if $sel != -1; # better score exists
				if ($sel>=1) {$selRank++;}
				elsif ($sel<=-2) {$rejRank++;}
			}
			if ($selRank) {$queryStatus=$selRank;}
			elsif ($rejRank==$totRank) {$queryStatus=0;}
			else {$queryStatus=-2;} # rejected + not verified
		}

		###<Updating table QUERY_VALIDATION>###
		if ($updateQuery) { # update necessary
			$newRefQuery=1 if (param("status_$rankID")==-1 && $rankID eq $refRankID); # A new query is being verified in current analysis
			$rankInfo=~s/SEL=-*\d/SEL=$select/; # starting SEL cannot be -1
			#$dbh->do("UPDATE QUERY_VALIDATION SET INFO_PEP$rank='$rankInfo',VALID_STATUS=$queryStatus WHERE ID_QUERY=$queryID") || die $dbh->errstr;
			$sthUpRank{$rank}->execute($rankInfo,$queryStatus,$queryID) || die $dbh->errstr;

		 	###<Fetching list of identifiers from RANK_PROTEIN_MATCH table>###
			if ($updateQuery==2) {
				$sthSEL1->execute($queryNum,$rank,$analysisID) || die $sthSEL1->errstr;
				while (my ($identifier) = $sthSEL1->fetchrow_array) {
					$listIdentifiers{$analysisID}{$identifier}=1;
				}
			}
		}
		$allRankInfo{$analysisID}{"$queryNum:$rank"}=$rankInfo; # no need for later retrieval from DB
	}
	$sthSEL1->finish;
	$sthSelInfo->finish;
	foreach my $rank (keys %sthUpRank) {$sthUpRank{$rank}->finish;}

	####<Updating info for all matched proteins>####
	my $sthSEL2=$dbh->prepare("SELECT ID_PROT_VALID,SEL_STATUS,MAX_MATCH FROM PROTEIN_VALIDATION WHERE IDENTIFIER=? AND ID_ANALYSIS=?");
	my $sthSEL3=$dbh->prepare("SELECT QUERY_NUM,PEP_RANK,MATCH_MULTI FROM RANK_PROTEIN_MATCH WHERE IDENTIFIER=? AND ID_ANALYSIS=?");
	my $sthUP2=$dbh->prepare("UPDATE PROTEIN_VALIDATION SET NUM_MATCH=?,SEL_STATUS=?,SCORE=?,CONF_LEVEL=2 WHERE ID_PROT_VALID=?");
	my %sthSelRk;
	foreach my $rank (1..$maxRank) {
		$sthSelRk{$rank}=$dbh->prepare("SELECT INFO_PEP$rank FROM QUERY_VALIDATION WHERE QUERY_NUM=? AND ID_ANALYSIS=?");
	}

	foreach my $analysisID (keys %listIdentifiers) {
		foreach my $identifier (keys %{$listIdentifiers{$analysisID}}) {

			###<Fetching protein data from table PROTEIN_VALIDATION>###
			$sthSEL2->execute($identifier,$analysisID) || die $sthSEL2->errstr;
			my ($protValidID,$oldSelStatus,$maxMatch)=$sthSEL2->fetchrow_array;

			next if $oldSelStatus <= -2; # protein is excluded/filtered from list => do not update selStatus

			###<Fetching info on all matching peptides>###
			my ($numMatch,$protScore,$rejectedMatch)=(0,0,0);
			my $selStatus;
			$sthSEL3->execute($identifier,$analysisID) || die $sthSEL3->errstr;
			while (my($queryNum,$rank,$matchFreq)=$sthSEL3->fetchrow_array) {
				my $qrKey="$queryNum:$rank";
				unless ($allRankInfo{$analysisID}{$qrKey}) { # rankInfo not already extracted from DB
					#($allRankInfo{$analysisID}{$qrKey})=$dbh->selectrow_array("SELECT INFO_PEP$rank FROM QUERY_VALIDATION WHERE QUERY_NUM=$queryNum AND ID_ANALYSIS=$analysisID");
					$sthSelRk{$rank}->execute($queryNum,$analysisID);
					($allRankInfo{$analysisID}{$qrKey})=$sthSelRk{$rank}->fetchrow_array;
				}
				if ($allRankInfo{$analysisID}{$qrKey} =~ /SEL=[12]/) {
					$numMatch++;
					my ($pepScore)=($allRankInfo{$analysisID}{$qrKey}=~/SC=(\d+\.*\d*)/); # Score corresponding to peptide
					$protScore+=($matchFreq*$pepScore);
				}
				elsif ($allRankInfo{$analysisID}{$qrKey} =~ /SEL=-[23]/) {$rejectedMatch++;}
			}

			if ($numMatch) { # at least 1 peptide selected
				$selStatus=($numMatch+$rejectedMatch==$maxMatch)? 2 : 1;
			}
			else { # no peptides selected
				$selStatus=($rejectedMatch==$maxMatch)? 0 : -1;
			}

			###<Updating number of selected proteins for main validation window>###
			if ($analysisID == $refAnalysisID) {
				if ($selStatus > 0 && $oldSelStatus <= 0) {$newRefProteins++;}
				elsif ($selStatus <= 0 && $oldSelStatus > 0) {$newRefProteins--;}
			}

			###<Updating protein entries in table PROTEIN_VALIDATION>###
			$sthUP2->execute($numMatch,$selStatus,$protScore,$protValidID) || die $sthUP2->errstr;
		} # end of (identifier,matchFreq) loop
	} # end of $rankID loop

	$sthSEL2->finish;
	$sthSEL3->finish;
	$sthUP2->finish;
	foreach my $rank (keys %sthSelRk) {$sthSelRk{$rank}->finish;}
	$dbh->commit;
	$dbh->disconnect;

	my $newValidDecoy=($newRefQuery || $newRefPeptide || $newRefProteins)? -2 : 0;

	##############
	####<HTML>####
	##############
	print header(-charset=>'utf-8');
	print qq
|<HTML><HEAD>
<SCRIPT LANGUAGE="JavaScript">
// !!! top.promsFrame or parent.parent because of subframes !!!
if (top.promsFrame.selectedView=='protein') {
	top.promsFrame.updateValidation(top.promsFrame.selectedProtId,$newRefQuery,$newRefPeptide,$newValidDecoy,$newRefProteins,top.promsFrame.selectedView,top.promsFrame.selProteinPage,1);
}
else { // view=query
	top.promsFrame.updateValidation(top.promsFrame.selectedQueryId,$newRefQuery,$newRefPeptide,$newValidDecoy,$newRefProteins,top.promsFrame.selectedView,top.promsFrame.selQueryPage,1);
}
</SCRIPT>
</HEAD></HTML>
|;
	exit;
}


#########################################################################################
##################################<<< FRAMES >>>#########################################
#########################################################################################
sub createFrames {
	####<Parameters>####
	my $selectedItem=param('MATCH');
	$selectedItem=~s/ana/seq/; # same as in listRanks.cgi
	my $matchInfo=param('INFO');
	my $disableForm=(param('DISABFORM'))? param('DISABFORM') : '';
	my ($ChargSts) = (param('ChargSts'))?param('ChargSts'):'';

	####<Starting HTML>####
	print header(-charset=>'utf-8');
		# warningsToBrowser(1);
	print qq
|<HTML>

<HEAD>
<SCRIPT LANGUAGE="JavaScript">
function loadFrame () {
top.promsFrame.spectrumFrame.matchAnaFrame.location="$promsPath{cgi}/listMatchedAnalyses.cgi?MATCH=$selectedItem&ChargSts=$ChargSts&INFO=$matchInfo&DISABFORM=$disableForm&sortType="+top.promsFrame.listMatchedAnalysesSort;
};
</SCRIPT>
</HEAD>

<FRAMESET rows="130,*" onload ="loadFrame ()" >
	<FRAME name="matchAnaFrame"  >
	<FRAME name="anaSpecFrame">
</FRAMESET>

</HTML>
|;
	exit;
}

#########################################################################################
##################################<<< Sort Analyses >>>##################################
#########################################################################################

sub sortAnalyses {
	if ($_[0] eq 'name') {lc($analysisParents{$analysisID{$a}}) cmp lc($analysisParents{$analysisID{$b}}) || lc($analysisName{$analysisID{$a}}) cmp lc($analysisName{$analysisID{$b}}) || $score{$b}<=>$score{$a}}
	elsif ($_[0] eq 'score'){$score{$b}<=>$score{$a} || lc($analysisName{$analysisID{$a}}) cmp lc($analysisName{$analysisID{$b}})}
	elsif ($_[0] eq 'masse'){$massObs{$b}<=>$massObs{$a} || lc($analysisName{$analysisID{$a}}) cmp lc($analysisName{$analysisID{$b}}) || $score{$b}<=>$score{$a}}
	elsif ($_[0] eq 'delta'){abs($massDelta{$a})<=>abs($massDelta{$b}) || lc($analysisName{$analysisID{$a}}) cmp lc($analysisName{$analysisID{$b}}) || $score{$b}<=>$score{$a}}
	elsif ($_[0] eq 'ppm'){abs($massDeltaPPM{$a})<=>abs($massDeltaPPM{$b}) || lc($analysisName{$analysisID{$a}}) cmp lc($analysisName{$analysisID{$b}}) || $score{$b}<=>$score{$a}}
}

####>Revision history<####
# 1.6.7 Check for undefined Analysis taxonomy (PP 04/03/16)
# 1.6.6 Remove VAR_MOD and VMOD from script (GA 23/05/13)
# 1.6.5 chomp on wiff file name because of remaining \n (PP 11/03/13)
# 1.6.4 Better handling of undefined elution time (PP 13/04/12)
# 1.6.3 Minor revision to print correctly the elution time (garras 01/02/2011)
# 1.6.2 Add charge state to MRobs popup (ppoullet 16/12/2010)
