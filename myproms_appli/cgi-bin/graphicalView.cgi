#!/usr/local/bin/perl -w

################################################################################
# graphicalView.cgi     2.2.9                                                  #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Draws a graphical view of selected proteins peptides distribution            #
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

#print header; warningsToBrowser(1); #DEBUG
#######################
####>Configuration<####
#######################
my $DEF_PROT_LENGTH=1000; # used if no protein length can be computed
my %promsPath=&promsConfig::getServerInfo;

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

###########################
####>Parameters 1 of 2<####
###########################
my $what=param('what');
my $itemID=param('id_item');
my $item=lc(param('item')); # if what=protein: id_item=ref analysis BUT item=project!
my $largWindow=param('largWindow');
#>peptide filter
my $specifFilter=(param('specifFilter'))? param('specifFilter') : 'All';
my $ptmFilter=(param('ptmFilter'))? param('ptmFilter') : 'Indifferent';
my $speciesFilter=(param('speciesFilter'))? param('speciesFilter') : 'All';
#>View
my $view=param('graphView'); #'CumComp'; # 'CumExp'; # 'StdComp'; # 'CumComp'; #'StdExp';
my $hideEmptyProt=(param('hideEmptyProt'))? 1 : 0;
my $includeFixedMods=(param('inclFixMods'))? 1 : 0;
#>Other options
my $chkIgnore=(param('varModBox'))? 'checked' : '';
my %listProteins;
my ($refAnaID,$refGroup)=(0,0); # for what=group,protein
my $projectAccess='';

####>$what = 'checkbox'<####
if ($what eq 'checkbox') {
	if (param('newView')) { # page reloaded due to change in display view
		my @dataList=split(/ /,param('listProt'));
		foreach my $data (@dataList) {
			my ($protID,$anaID)=split(/:/,$data);
			$listProteins{$protID}{$anaID}=1;
		}
	}
	else { # 1st call
		###>Fetching data from list of checkbox fields<### param('listMode') no longer required
		foreach my $info (param('chkProt')) {
			my @infoData=split(':',$info); # anaID1:...:anaIDn:protID
			for (my $i=0; $i<$#infoData; $i++) {
				$listProteins{$infoData[-1]}{$infoData[$i]}=1;
			}
		}
	}
}

####>$what = 'group'<####
elsif ($what eq 'group') {
	my $userID=$ENV{'REMOTE_USER'};
	my $projectID=&promsMod::getProjectID($dbh,$itemID,$item);
	my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
	$projectAccess=${$userInfo[2]}{$projectID};
	($refAnaID,$refGroup)=split(/:/,param('groupInfo'));
	my $sth=$dbh->prepare("SELECT ID_PROTEIN FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=$refAnaID AND MATCH_GROUP=$refGroup");
	$sth->execute();
	while (my ($protID)=$sth->fetchrow_array) {
		$listProteins{$protID}{$refAnaID}=1;
	}
	$sth->finish;
}
####>$what = 'protein'<####
elsif ($what eq 'protein') { # item set to project
	my $protID=param('id_prot');
	my @anaList=(param('newView'))? split(/ /,param('listAna')) : param('check_ana'); #all selected analyses where protein is found
	foreach my $anaID (@anaList) {
		$listProteins{$protID}{$anaID}=1;
	}
	$refAnaID=$itemID;
}

##############
####>Main<####
##############
my $ITEM=uc($item);
my $nameItem;
if ($what ne 'protein') {
	($nameItem)=$dbh->selectrow_array("SELECT NAME FROM $ITEM WHERE ID_$ITEM=$itemID");
}

####>Preparing all queries<####
my $sthProt1=$dbh->prepare("SELECT PROT_LENGTH,IDENTIFIER,PROT_DES,ALIAS,ORGANISM FROM PROTEIN WHERE ID_PROTEIN=?");
my $sthPep1=$dbh->prepare("SELECT DISTINCT ID_PEPTIDE FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_ANALYSIS=? AND ID_PROTEIN=?");
#my $sthPep2=$dbh->prepare("SELECT PEP_SEQ,PEP_LENGTH,SCORE,VALID_STATUS FROM PEPTIDE WHERE ID_PEPTIDE=? AND VALID_STATUS>=1");
my $sthPep2=$dbh->prepare("SELECT PEP_SEQ,CHARGE,PEP_LENGTH,SCORE,VALID_STATUS FROM PEPTIDE WHERE ID_PEPTIDE=?");
my $sthPep3=$dbh->prepare("SELECT ABS(PEP_BEG) FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_ANALYSIS=? AND ID_PROTEIN=? AND ID_PEPTIDE=?");
my $sthPep4=$dbh->prepare("SELECT COUNT(DISTINCT ID_PROTEIN) FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_ANALYSIS=? AND ID_PEPTIDE=?");
my $sthProt2=$dbh->prepare("SELECT SCORE,VISIBILITY,CONF_LEVEL FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=? AND ID_PROTEIN=?");
my $sthAM=$dbh->prepare("SELECT AM.ID_MODIFICATION,MODIF_TYPE,AM.SPECIFICITY,PSI_MS_NAME,INTERIM_NAME FROM ANALYSIS_MODIFICATION AM, MODIFICATION M WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND ID_ANALYSIS=?");
my $sthPM=$dbh->prepare("SELECT DISTINCT PM.ID_MODIFICATION,PM.MODIF_TYPE,SPECIFICITY,PSI_MS_NAME,INTERIM_NAME FROM MODIFICATION M, PEPTIDE_MODIFICATION PM, PEPTIDE P WHERE M.ID_MODIFICATION=PM.ID_MODIFICATION AND PM.ID_PEPTIDE=P.ID_PEPTIDE AND P.ID_ANALYSIS=?");
my $sthPepMod=$dbh->prepare("SELECT DISTINCT ID_MODIFICATION FROM PEPTIDE_MODIFICATION WHERE ID_PEPTIDE=?");
my @itemQueries; # Find item hierarchy in project
if ($item eq 'project') {
	$itemQueries[0]="SELECT EXPERIMENT.NAME,SAMPLE.NAME,ANALYSIS.NAME FROM EXPERIMENT,SAMPLE,ANALYSIS WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_EXPERIMENT=EXPERIMENT.ID_EXPERIMENT AND ID_ANALYSIS=?";
}
elsif ($item eq 'experiment') {
	$itemQueries[0]="SELECT SAMPLE.NAME,ANALYSIS.NAME FROM SAMPLE,ANALYSIS WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT IS NULL AND ID_ANALYSIS=?";
	$itemQueries[1]="SELECT GEL2D.NAME,SPOT.NAME,ANALYSIS.NAME FROM GEL2D,SPOT,SAMPLE,ANALYSIS WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT=SPOT.ID_SPOT AND SPOT.ID_GEL2D=GEL2D.ID_GEL2D AND ID_ANALYSIS=?";
}
elsif ($item eq 'gel2d') {
	$itemQueries[0]="SELECT SPOT.NAME,ANALYSIS.NAME FROM SPOT,SAMPLE,ANALYSIS WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT=SPOT.ID_SPOT AND ID_ANALYSIS=?";
}
else { # sample, spot or analysis
	$itemQueries[0]="SELECT NAME FROM ANALYSIS WHERE ID_ANALYSIS=?";
}
my @sthItems;
foreach my $itemQuery (@itemQueries) {
	push @sthItems,$dbh->prepare($itemQuery);
}

####>Executing all queries<#####
my (%protLength,%identifier,%desProt,%alias,%organism,%peptideMatch,%listPos,%confidence,%peptideSpec,%numIdentifications,%speciesList);
my (%score,%visibility,%ptmList,%allFixedMods,%trueVarMods);
#my %infoSubmit=promsMod::getSearchParam($dbh,$analysisID);
#push @modifications,sort(split(/,/,$infoSubmit{'g:Variable modifications'}));

my %analysisString;
my $matchGrAlias;
foreach my $protID (keys %listProteins) {
	$sthProt1->execute($protID);
	($protLength{$protID},$identifier{$protID},$desProt{$protID},$alias{$protID},$organism{$protID})=$sthProt1->fetchrow_array;
	$speciesList{$organism{$protID}}{$protID}=1;

	foreach my $anaID (keys %{$listProteins{$protID}}) {

		##>Fixed modifications<##
		my %fixedMods;
		if ($includeFixedMods) { # Add fixed mods to list of PTMs
			#my %infoSubmit=promsMod::getSearchParam($dbh,$anaID);
			#if ($infoSubmit{'f:Fixed modifications'}) {
			#	foreach my $mod (split(/,/,$infoSubmit{'f:Fixed modifications'})) {
			#		my ($type)=($mod=~/\(([^:\)]+)/);
			#		$type='N-term' if $type eq 'Any N-Terminus';
			#		$mod=~s/\s*\(.+//;
			#		my $fixedModName="$mod ($type)";
			#		$ptmList{$fixedModName}=1;
			#		$fixedMods{$fixedModName}=1;
			#		$allFixedMods{$fixedModName}=1; # records all fixedMods
			#	}
			#}
			$sthAM->execute($anaID);
			while (my ($modID,$modifType,$specificity,$psiName,$interName)=$sthAM->fetchrow_array) {
				next unless $modifType eq 'F';
				my $modName=($psiName)?$psiName:($interName)?$interName:"mod_$modID";
				$ptmList{$modID}{'SPECIFICITY'}=$specificity;
				$ptmList{$modID}{'NAME'}=$modName;
				$fixedMods{$modID}=1;
				$allFixedMods{$modID}=1;
			}
		}
		$sthPM->execute($anaID);
		while (my ($modID,$modifType,$specificity,$psiName,$interName)=$sthPM->fetchrow_array) {
			my $modName=($psiName)?$psiName:($interName)?$interName:"mod_$modID";
			$ptmList{$modID}{'SPECIFICITY'}=$specificity;
			$ptmList{$modID}{'NAME'}=$modName;
		}
		my $cumAna=($view=~/Std/ || $what eq 'group')? $anaID : 0;

		%{$peptideMatch{$protID}{$cumAna}}=() unless $peptideMatch{$protID}{$cumAna}; # required for display of unmatched proteins
		%{$listPos{$protID}{$cumAna}}=() unless $listPos{$protID}{$cumAna}; # required for display of unmatched proteins

		if ($speciesFilter eq 'All' || $organism{$protID} eq $speciesFilter) { # species filter
			$sthPep1->execute($anaID,$protID);
			while (my ($pepID)=$sthPep1->fetchrow_array) {
				$sthPep2->execute($pepID);
				my ($seq,$charge,$length,$score,$valStatus)=$sthPep2->fetchrow_array;
				#next unless $valStatus>=1; # skip ghost peptides
				$score=0 unless $score;
				my $varMod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$anaID,$seq);
				#>Adding fixed Mods
				foreach my $fixModID (keys %fixedMods) {
					if ($ptmList{$fixModID}{'SPECIFICITY'} =~ /=|-|\+|\*/ ) {$varMod.=" + $ptmList{$fixModID}{'NAME'}";} # end modifications
					elsif ($ptmList{$fixModID}{'SPECIFICITY'}) { # specific aa modifications
						foreach my $aa (split(/\./,$ptmList{$fixModID}{'SPECIFICITY'})) {
							$varMod.=" + $ptmList{$fixModID}{'NAME'}" if $seq=~/$aa/;
						}
					}
				}
				if ($varMod) {
					my $okVarMod=($ptmFilter eq 'Indifferent' || $ptmFilter eq 'Any')? 1 : 0;
					#my @modifs=split(/\s\+\s/,$varMod);
					#shift @modifs;
					$sthPepMod->execute($pepID);
					while (my ($modID) = $sthPepMod->fetchrow_array) {
						$okVarMod=1 if $ptmList{$modID}{'NAME'} eq $ptmFilter;
						$trueVarMods{$modID}=1 unless $fixedMods{$modID}; # records true varMods
					}
					#foreach my $mod (@modifs) {
					#	my ($type)=($mod=~/\(([^:\)]+)/);
					#	$type='N-term' if $type eq 'Any N-Terminus';
					#	$mod=~s/\s*\(.+//;
					#	my $varModName="$mod ($type)";
					#	$okVarMod=1 if $varModName eq $ptmFilter;
					#	$ptmList{$varModName}=1;
					#	$trueVarMods{$varModName}=1 unless $fixedMods{$varModName}; # records true varMods
					#}
					next if (!$okVarMod || $ptmFilter eq 'None');
				}
				elsif ($ptmFilter ne 'Indifferent' && $ptmFilter ne 'None') {next;} # filter is set on a specific PTM

				@{$peptideMatch{$protID}{$cumAna}{$pepID}}=($seq,$charge,$length,$varMod,$score);
				$numIdentifications{"$protID:$cumAna"}=1 if $specifFilter eq 'Shared';

				$sthPep3->execute($anaID,$protID,$pepID);
				while (my ($pepBeg)=$sthPep3->fetchrow_array) {
					$listPos{$protID}{$cumAna}{$seq}{$pepBeg}++;
					push @{$peptideSpec{$seq}{"$protID:$cumAna"}{$pepID}},$pepBeg if $specifFilter ne 'All'; # recording specificity (hash because of possible duplication due to PMT or lower scoring)
				}
			}
		}
		$sthProt2->execute($anaID,$protID);
		my ($score,$vis,$conf)=$sthProt2->fetchrow_array;
		$score{$protID}{$cumAna}=$score if (!$score{$protID}{$cumAna} || $score{$protID}{$cumAna}<$score);
		$visibility{$protID}{$cumAna}=$vis if (!$visibility{$protID}{$cumAna} || $visibility{$protID}{$cumAna}<$vis);
		$confidence{$protID}{$cumAna}=$conf if (!$confidence{$protID}{$cumAna} || $confidence{$protID}{$cumAna}<$conf);

		if ($view=~/Std/) {
			my @names;
			foreach my $sthItem (@sthItems) {
				$sthItem->execute($anaID);
				@names=$sthItem->fetchrow_array;
				last if scalar @names; # experiment to analysis through sample or gel
			}
			$analysisString{$anaID}=join(' > ',@names);
		}
		else {
			$analysisString{$cumAna}='Cumulated MS/MS Analyses';
		}

	}
	if ($what eq 'protein') {$nameItem=$alias{$protID};}
	elsif ($what eq 'group' && $visibility{$protID}{$refAnaID}==2) {
		$nameItem="<FONT color='#DD0000'>Match Group of $alias{$protID}</FONT>";
	}
}
$sthProt1->finish;
$sthPep1->finish;
$sthPep2->finish;
$sthPep3->finish;
$sthProt2->finish;
$sthAM->finish;
$sthPM->finish;
$sthPepMod->finish;
foreach my $sthItem (@sthItems) {$sthItem->finish;}

$dbh->disconnect;

####>Finding longest protein sequence for later scale calculation<#### (moved before filtering in case no peptides left)
my $maxProtLength=0;
my %usedLength;
foreach my $protID (keys %listProteins){
	if ($protLength{$protID} == 0){
		my $maxMatchEnd=0;
		foreach my $anaID (keys %{$peptideMatch{$protID}}) {
			foreach my $pepID (keys %{$peptideMatch{$protID}{$anaID}}) {
				#my ($pepSeq,$pepLength)=@{$peptideMatch{$protID}{$anaID}{$pepID}}[0,1];
				my ($pepSeq,$pepLength)=@{$peptideMatch{$protID}{$anaID}{$pepID}}[0,2];
				foreach my $matchPos (sort{$b<=>$a} keys %{$listPos{$protID}{$anaID}{$pepSeq}}){
					my $posEnd=$matchPos+$pepLength-1;
					$maxMatchEnd=$posEnd if $maxMatchEnd<$posEnd;
					last; # count only biggest matchPos
				}
			}
		}
		$usedLength{$protID}=$maxMatchEnd;
	}
	else {$usedLength{$protID}=$protLength{$protID};}
	$maxProtLength=$usedLength{$protID} if $maxProtLength < $usedLength{$protID};
}
$maxProtLength=$DEF_PROT_LENGTH unless $maxProtLength; # no matching peptides left

####>Filtering on specificity<####
if ($specifFilter ne 'All') {
	my $numIdent=scalar keys %numIdentifications;
	foreach my $pepSeq (keys %peptideSpec) {
		if (($specifFilter eq 'Unique' && scalar keys %{$peptideSpec{$pepSeq}} > 1) || ($specifFilter eq 'Shared' && scalar keys %{$peptideSpec{$pepSeq}} < $numIdent)) {
			###>Removing peptide matches
			foreach my $protAna (keys %{$peptideSpec{$pepSeq}}) {
				my ($protID,$anaID)=split(':',$protAna);
				foreach my $pepID (keys %{$peptideSpec{$pepSeq}{$protAna}}) {
					my $seq=$peptideMatch{$protID}{$anaID}{$pepID}[0];
					foreach my $beg (@{$peptideSpec{$seq}{$protAna}{$pepID}}) {
						$listPos{$protID}{$anaID}{$seq}{$beg}--;
						delete $listPos{$protID}{$anaID}{$seq}{$beg} unless $listPos{$protID}{$anaID}{$seq}{$beg};
						delete $listPos{$protID}{$anaID}{$seq} unless scalar keys %{$listPos{$protID}{$anaID}{$seq}};
					}
					delete $peptideMatch{$protID}{$anaID}{$pepID};
				}
			}
		}
	}
}

####>Calculating peptide coverage<####
my %pepCoverage;
my $cumCoverage;
my %compMatches;
foreach my $protID (keys %listPos) {
	my %cumStatus;
	foreach my $anaID (keys %{$listPos{$protID}}) {
#next unless scalar keys %{$peptideMatch{$protID}{$anaID}}; # filter effect
		foreach my $pepSeq (keys %{$listPos{$protID}{$anaID}}){
			my $pepLength=length($pepSeq);
			foreach my $beg (keys %{$listPos{$protID}{$anaID}{$pepSeq}}){
				my $end=$beg+$pepLength-1;
				$cumStatus{$beg}+=$listPos{$protID}{$anaID}{$pepSeq}{$beg};
				$cumStatus{$end}-=$listPos{$protID}{$anaID}{$pepSeq}{$beg};
				$compMatches{$protID}{$anaID}{$beg}+=$listPos{$protID}{$anaID}{$pepSeq}{$beg};
				$compMatches{$protID}{$anaID}{$end}-=$listPos{$protID}{$anaID}{$pepSeq}{$beg};
			}
		}
		$pepCoverage{$protID}{$anaID}=&promsMod::getCoverage($protLength{$protID},\%{$compMatches{$protID}{$anaID}}) if $protLength{$protID};
	}
	$cumCoverage=&promsMod::getCoverage($protLength{$protID},\%cumStatus) if ($protLength{$protID} && $what eq 'protein');
}


#######################
####>Starting HTML<####
#######################
my ($color1,$color2)=&promsConfig::getRowColors;
print header(-charset=>'utf-8');
warningsToBrowser(1); # writes PERL warnings as HTML comments <!-- --> must be activated after header
print qq
|<HEAD>
<LINK rel="icon" type="image/png" href="$promsPath{images}/myProMS_icon.png">
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
A {text-decoration:none}
ACRONYM {cursor: help;}
</STYLE>
<TITLE>Peptide Distribution</TITLE>
<SCRIPT LANGUAGE="JavaScript">
var groupList=new Array;
var oldSequence;
var oldVarMod;
|;
##>Array with list of fixedMods NOT in varMods
print "var onlyFixedMods=[";
my $countFixedMods=0;
foreach my $modID (sort{lc($ptmList{$a}{'NAME'}) cmp lc($ptmList{$b}{'NAME'})} keys %allFixedMods) {
	next if $trueVarMods{$modID};
	$countFixedMods++;
	print ',' if $countFixedMods>1;
	print "\"$ptmList{$modID}{'NAME'}\"";
}
print "];\n";
&promsMod::popupInfo(250,'rel',0,'rel',20);
print qq
|function viewProt(idProt,idAna){
	window.location="./sequenceView.cgi?id_ana="+idAna+"&id_prot="+idProt+"&msdata="+opener.showMSdata;
}
function showPeptides(sequence,varMod) {
	if (!sequence) {return;}
	//Set all matching sequences to normal
	if (oldSequence) {
		var pepImg=document.getElementsByTagName("IMG");
		for (i=0; i<pepImg.length; i++) {
			//var imgSrc=(pepImg[i].name.match(' \+ '))? 'prot_FFFF00.gif' : 'peptide.gif';
			var pepSeq=pepImg[i].name.split(' + ');
			var imgSrc=(pepSeq.length > 1)? 'prot_FFFF00.gif' : 'peptide.gif';
			if (pepSeq[0]==oldSequence) {
				pepImg[i].src='$promsPath{images}/'+imgSrc;
			}
		}
	}
	//Highlight selected peptides
	if (document.viewForm.varModBox.checked) { // ignore var mods
		var pepImg=document.getElementsByTagName("IMG");
		for (i=0; i<pepImg.length; i++) {
			var pepSeq=pepImg[i].name.split(' + ');
			if (pepSeq[0]==sequence) {
				pepImg[i].src='$promsPath{images}/peptideBlue.gif';
			}
		}
	}
	else { // include var mods
		var pepImg=document.getElementsByName(sequence+varMod);
		for (i=0; i<pepImg.length; i++) {
			pepImg[i].src='$promsPath{images}/peptideBlue.gif';
		}
	}
	oldSequence=sequence;
	oldVarMod=varMod;
}
function hideShowEmptyProteins(hideEmpty) {
	var protAnaDispStatus=(hideEmpty)? 'none' : 'block';
	var allDivList=document.getElementsByTagName('DIV'); // getElementsByName does NOT work with IE!!!
	for (var i=0; i<allDivList.length; i++) {
		if (allDivList[i].getAttribute('name')=='protAnaDiv') { // allDivList[i].name does NOT work with Firefox!!!
			if (allDivList[i].className=='empty') {allDivList[i].style.display=protAnaDispStatus;}
		}
		else if (allDivList[i].getAttribute('name')=='protDiv') { // Only in checkbox or protein mode (what)
			var protDispStatus='none'; // default
			if (hideEmpty) {
				var childList=allDivList[i].childNodes;
				for (var j=0; j<childList.length; j++) {
					if (childList[j].tagName=='DIV' && childList[j].className=='matched') { // 1 visible child is enough
						protDispStatus='block';
						break;
					}
				}
			}
			else {protDispStatus='block';}
			allDivList[i].style.display=protDispStatus;
		}
	}
}
function includeFixedModifications(chkStatus) {
	if (!chkStatus) {
		var myForm=document.viewForm;
		for (var i=0; i<onlyFixedMods.length; i++) {
			if (onlyFixedMods[i]==myForm.ptmFilter.value) { // a only fixedMod is selected
				myForm.ptmFilter.selectedIndex=0; // set to Indifferent
				break;
			}
		}
	}
	changeView();
}
function changeView() {
	var myForm=document.viewForm;
	opener.graphView=myForm.graphView.value;
	myForm.largWindow.value=(document.body)? window.document.body.clientWidth : window.innerWidth;
	myForm.submit();
}
</SCRIPT>
</HEAD>
<BODY bgcolor='white' background='$promsPath{images}/bgProMS.gif'>
<CENTER><FONT style="font-size:22px;font-weight:bold;">Peptide Distribution - $nameItem</FONT></CENTER>
<INPUT type="button" value="<< Back" onclick="history.back()">
<FORM name="viewForm" method="post">
<INPUT type="hidden" name="newView" value="1">
<INPUT type="hidden" name="what" value="$what">
<INPUT type="hidden" name="item" value="$item">
<INPUT type="hidden" name="id_item" value="$itemID">
<INPUT type="hidden" name="largWindow">
|;
if ($what eq 'checkbox') { # From a protein list
	print "<INPUT type=\"hidden\" name=\"listProt\" value=\"";
	foreach my $protID (keys %listProteins) {
		foreach my $anaID (keys %{$listProteins{$protID}}) {
			print "$protID:$anaID ";
		}
	}
	print "\">\n";
}
elsif ($what eq 'group') { # From a match group
	print "<INPUT type=\"hidden\" name=\"groupInfo\" value=\"",param('groupInfo'),"\">\n";
}
else { # From sequenceView ($what eq 'protein')
	print "<INPUT type=\"hidden\" name=\"id_prot\" value=\"",param('id_prot'),"\">\n";
	my @listAna=(param('newView'))? param('listAna') : param('check_ana');
	print "<INPUT type=\"hidden\" name=\"listAna\" value=\"@listAna\">\n";
}

print "<TABLE border=0 align=center><TR><TD bgcolor=$color2>\n";
print "<TABLE border=0>\n";

###>View<###
my $viewString='<B>Standard:</B> Each MS/MS Analysis with corresponding peptide matches is displayed separatly.<BR>';
$viewString.='<B>Cumulated:</B> All MS/MS Analyses with corresponding peptide matches are cumulated.<HR width=80%>';
$viewString.='<B>Expanded:</B> Matches for each peptide are displayed on separate lines.<BR>';
$viewString.='<B>Compact:</B> Peptide matches are projected on protein.';
print "<TR><TH align=right nowrap><FONT class=\"title2\"><ACRONYM onmouseover=\"popup('$viewString')\" onmouseout=\"popout()\">View<SUP><SMALL>*</SMALL></SUP></ACRONYM> :</FONT></TH>\n";
print "<TD colspan=3 bgcolor=$color1 nowrap><SELECT name=\"graphView\" style=\"font-weight:bold;font-size:16\" onchange=\"changeView()\">\n";
my $selecString=($view eq 'StdExp')? 'selected' : '';
print "<OPTION value=\"StdExp\" $selecString>Standard Expanded</OPTION>\n";
$selecString=($view eq 'StdComp')? 'selected' : '';
print "<OPTION value=\"StdComp\" $selecString>Standard Compact</OPTION>\n";
$selecString=($view eq 'CumExp')? 'selected' : '';
print "<OPTION value=\"CumExp\" $selecString>Cumulated Expanded</OPTION>\n";
$selecString=($view eq 'CumComp')? 'selected' : '';
print "<OPTION value=\"CumComp\" $selecString>Cumulated Compact</OPTION>\n";
my $chkUnMatchStatus=($hideEmptyProt)? 'checked' : '';
my $chkInclFixStatus=($includeFixedMods)? 'checked' : '';
print qq
|</SELECT>
&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="hideEmptyProt" value="1" onclick="hideShowEmptyProteins(this.checked)" $chkUnMatchStatus><FONT class="title3">Hide unmatched proteins</FONT>
&nbsp;&nbsp;&nbsp;<INPUT type="checkbox" name="inclFixMods" value="1" onclick="includeFixedModifications(this.checked)" $chkInclFixStatus><FONT class="title3">Include fixed modifications</FONT>&nbsp;&nbsp;
</TD></TR>
|;

###>Filters<###
print qq
|<TR><TH align=right nowrap><FONT class="title2">Filters :</FONT></TH>
<TD bgcolor=$color1 nowrap><FONT class="title3">Species:</FONT><SELECT name="speciesFilter" class="title3" onchange="changeView()"><OPTION value="All">All</OPTION>
|;
foreach my $species (sort{lc($a) cmp lc($b)} keys %speciesList) {
	my $selSpec=($species eq $speciesFilter)? 'selected' : '';
	print "<OPTION value=\"$species\" $selSpec>$species</OPTION>";
}
my ($disUniqStr,$disSharStr)=($specifFilter eq 'Unique')? (' selected','') : ($specifFilter eq 'Shared')? ('',' selected') : ('','');
print qq
|</SELECT></TD>
<TD bgcolor=$color1 nowrap><FONT class="title3">&nbsp;Peptides in matched set:</FONT><SELECT name="specifFilter" class="title3" onchange="changeView()"><OPTION value="All">All</OPTION><OPTION value="Unique"$disUniqStr>Unique</OPTION><OPTION value="Shared"$disSharStr>Shared</OPTION></SELECT></TD>
<TD bgcolor=$color1 nowrap><FONT class="title3">&nbsp;PTMs:</FONT><SELECT name="ptmFilter" class="title3" onchange="changeView()">
|;
my @modifications=('Indifferent','None','Any');
foreach my $modName (@modifications) {
	my $selMod=($modName eq $ptmFilter)? 'selected' : '';
	print "<OPTION value=\"$modName\" $selMod>$modName</OPTION>";
}
@modifications=sort{lc($ptmList{$a}{'NAME'}) cmp lc($ptmList{$b}{'NAME'})} keys %ptmList;
foreach my $modID (@modifications) {
	my $selMod=($ptmList{$modID}{'NAME'} eq $ptmFilter)? 'selected' : '';
	print "<OPTION value=\"$ptmList{$modID}{'NAME'}\" $selMod>$ptmList{$modID}{'NAME'}</OPTION>";
}
print qq
|</SELECT></TD>\n</TR>
</TABLE></TD></TR></TABLE>
|; #&nbsp;<INPUT type="button" style="font-weight:bold" value="Proceed" onclick="changeView()"/>

$largWindow=0.95 * $largWindow;
if ($what eq 'checkbox' || $what eq 'protein') {
	if ($view=~/Exp/) {&pepLegend;} else {print "<BR>\n";}
	foreach my $protID (sort {lc($identifier{$a}) cmp lc($identifier{$b}) || $a<=>$b} keys %listProteins) {
		my $lengthString;
		if(!$protLength{$protID}){$lengthString='unknown size';}
		else{$lengthString="$protLength{$protID} aa";}
		my $protVisStatus=($hideEmptyProt)? 'none' : 'block';
		#>First pass
		my %numPepAna;
		foreach my $anaID (keys %{$listPos{$protID}}) {
			$numPepAna{$anaID}=scalar keys %{$peptideMatch{$protID}{$anaID}};
			$protVisStatus='block' if $numPepAna{$anaID};
		}
		print "<DIV name=\"protDiv\" style=\"display:$protVisStatus\"><B>$alias{$protID}</B> : $desProt{$protID} ($lengthString) <I><U>$organism{$protID}</U></I><BR>\n";
		#>2nd pass
		foreach my $anaID (sort{$visibility{$protID}{$b}<=>$visibility{$protID}{$a} || scalar (keys %{$listPos{$protID}{$b}})<=>scalar (keys %{$listPos{$protID}{$a}}) || $score{$protID}{$b}<=>$score{$protID}{$a}} keys %{$listPos{$protID}}) {
			#my $className=($peptideMatch{$protID} && $peptideMatch{$protID}{$anaID})? 'matched' : 'empty';
			my $className=($numPepAna{$anaID})? 'matched' : 'empty';
			my $visStatus=(!$hideEmptyProt || $className eq 'matched')? 'block' : 'none';
			print "<DIV name=\"protAnaDiv\" class=\"$className\" style=\"display:$visStatus\">\n";
			my ($class,$titleString)=&promsMod::getProteinClass($confidence{$protID}{$anaID},$visibility{$protID}{$anaID});
			if ($view=~/Std/) {
				print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<A class=\"$class\" href=javascript:viewProt($protID,$anaID) onmouseover=\"popup('$titleString')\" onmouseout=\"popout()\">$analysisString{$anaID}</A>&nbsp(score: $score{$protID}{$anaID}";
			}
			else {
				print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<A class=\"$class\">$analysisString{$anaID}</A>&nbsp(best score: $score{$protID}{$anaID}";
			}
			if ($pepCoverage{$protID}{$anaID}) {
				if ($view=~/Std/) {printf ", coverage: %.1f",$pepCoverage{$protID}{$anaID};}
				else {printf ", cumulated coverage: %.1f",$pepCoverage{$protID}{$anaID};}
				print ' %';
			}
			print ")<BR>\n";
			if ($view=~/Exp/) {
				&drawExpMatches($protID,$anaID);
			}
			else {
				&drawCompMatches($protID,$anaID);
			}
			print "</DIV>\n";
		}
		print "</DIV>\n";
	}
}
elsif ($what eq 'group') {
	if ($view=~/Exp/) {&pepLegend;} else {print "<BR>\n";}
	foreach my $protID (sort{$visibility{$b}{$refAnaID}<=>$visibility{$a}{$refAnaID} || scalar (keys %{$listPos{$b}{$refAnaID}})<=>scalar (keys %{$listPos{$a}{$refAnaID}}) || $score{$b}{$refAnaID}<=>$score{$a}{$refAnaID} || lc($identifier{$a}) cmp lc($identifier{$b})} keys %listProteins) {
		#my $className=($peptideMatch{$protID} && $peptideMatch{$protID}{$refAnaID})? 'matched' : 'empty';
		my $className=(scalar keys %{$peptideMatch{$protID}{$refAnaID}})? 'matched' : 'empty';
		my $visStatus=(!$hideEmptyProt || $className eq 'matched')? 'block' : 'none';
		print "<DIV name=\"protAnaDiv\" class=\"$className\" style=\"display:$visStatus\">\n";
		my $lengthString;
		if (!$protLength{$protID}) {$lengthString='unknown size';}
		else {$lengthString="$protLength{$protID} aa";}
		my ($class,$titleString)=&promsMod::getProteinClass($confidence{$protID}{$refAnaID},$visibility{$protID}{$refAnaID});
		print "<A class=\"$class\" href=javascript:viewProt($protID,$refAnaID) onmouseover=\"popup('$titleString')\" onmouseout=\"popout()\">$alias{$protID}</A>&nbsp(score: $score{$protID}{$refAnaID}";
		if ($pepCoverage{$protID}{$refAnaID}) {
			printf ", coverage: %.1f",abs($pepCoverage{$protID}{$refAnaID});
			print ' %';
		}
		print "): $desProt{$protID} ($lengthString) <I><U>$organism{$protID}</U></I><BR>\n";
		if ($view=~/Exp/) {
			&drawExpMatches($protID,$refAnaID);
		}
		else {
			&drawCompMatches($protID,$refAnaID);
		}
		print "</DIV>\n";
	}
}

if ($what eq 'protein' && $view=~/Std/) {
	print '<B>Cumulated peptide coverage:';
	if ($cumCoverage) {printf " %.1f",$cumCoverage; print ' %';} else {print ' Unknown.';}
	print "</B><BR>\n";
}

print "<INPUT type=\"button\" value=\"<< Back\" onClick=\"history.back()\">\n" if scalar keys %listProteins > 12;

print qq
|</FORM>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup()</SCRIPT>

</BODY>
</HTML>
|;



#################################
####<<<Drawing subroutines>>>####
#################################
sub pepLegend {
	print qq
|<BR><B>
Click on a peptide to highlight all its occurrences
(<INPUT type="checkbox" name="varModBox" onclick="showPeptides(oldSequence,oldVarMod)" $chkIgnore> Ignore PTMs)
</B>
<TABLE><TR>
<TH align=left><IMG src="$promsPath{images}/peptide.gif" width=20 height=7 border=0>: default peptide&nbsp&nbsp&nbsp</TH>
<TH align=left><IMG src="$promsPath{images}/prot_FFFF00.gif" width=20 height=7 border=0>: modified peptide&nbsp&nbsp&nbsp</TH>
<TH align=left><IMG src="$promsPath{images}/peptideBlue.gif" width=20 height=7 border=0>: highlighted peptide</TH>
</TR></TABLE><BR>
|;
}
sub drawExpMatches {
	my ($protID,$anaID)=@_;
	unless ($usedLength{$protID}) { # no prot length & no matching peptide after filtering!
		&drawUnknownSize($DEF_PROT_LENGTH);
		print "<BR>\n";
		return;
	}
	my $newProtLength=int(0.5+($largWindow*($usedLength{$protID}/$maxProtLength)));
	print "<IMG src=\"$promsPath{images}/protein.gif\" width=$newProtLength height=7 border=0 onmouseover=\"popup('<B>$identifier{$protID}</B>')\" onmouseout=\"popout()\">";
	&drawUnknownSize(0) unless $protLength{$protID};
	print "<BR>\n";
	#foreach my $pepID (sort{$peptideMatch{$protID}{$anaID}{$b}[1]<=>$peptideMatch{$protID}{$anaID}{$a}[1] || $peptideMatch{$protID}{$anaID}{$a}[0] cmp $peptideMatch{$protID}{$anaID}{$b}[0]  || $peptideMatch{$protID}{$anaID}{$a}[2] cmp $peptideMatch{$protID}{$anaID}{$b}[2] || $peptideMatch{$protID}{$anaID}{$b}[3]<=>$peptideMatch{$protID}{$anaID}{$a}[3]} keys %{$peptideMatch{$protID}{$anaID}}) { # length,seq,varMod,score
	foreach my $pepID (sort{$peptideMatch{$protID}{$anaID}{$b}[2]<=>$peptideMatch{$protID}{$anaID}{$a}[2] || $peptideMatch{$protID}{$anaID}{$a}[0] cmp $peptideMatch{$protID}{$anaID}{$b}[0]  || $peptideMatch{$protID}{$anaID}{$a}[3] cmp $peptideMatch{$protID}{$anaID}{$b}[3] || $peptideMatch{$protID}{$anaID}{$b}[4]<=>$peptideMatch{$protID}{$anaID}{$a}[4]} keys %{$peptideMatch{$protID}{$anaID}}) { # length,seq,varMod,score
		my ($pepSeq,$pepCharge,$pepLength,$varMod,$pepScore)=@{$peptideMatch{$protID}{$anaID}{$pepID}};
		$pepScore='N/A (hidden peptide)' unless $pepScore;
		my $length=int(0.5+($pepLength*$newProtLength)/$usedLength{$protID});
		my %varModPos;
		if($varMod){
			my $varModString = $varMod;
			$varModString =~s/^\s\+\s// ; # string starts with ' + ';
			my @pepSeq = split //, $pepSeq;
			foreach my $var (split / \+ /, $varModString){
				my ($varModCode,$resid,$positions) = &promsMod::convertVarModString($var,0);
				if($resid eq '-'){
					push @{$varModPos{$varModCode}}, "Protein N-term";
				} elsif ($resid eq '+'){
					push @{$varModPos{$varModCode}}, "Protein C-term";
				} elsif ($resid eq '='){
					push @{$varModPos{$varModCode}}, "N-term";
				} elsif ($resid eq '*'){
					push @{$varModPos{$varModCode}}, "C-term";
				} else {
					foreach my $pos (split /\./, $positions){
						next if $pos eq '?'; # do not consider ambigous positions
						my $residue = $pepSeq[$pos-1]; #fetch corresponding residue on peptide sequence
						push @{$varModPos{$varModCode}}, "$residue:$pos";
					}
				}
			}
		}
		foreach my $pepBeg (sort{$a<=>$b} keys %{$listPos{$protID}{$anaID}{$pepSeq}}){
			my $space=int(0.5+(($pepBeg-1)*$newProtLength)/$usedLength{$protID});
			print "<IMG src=\"$promsPath{images}/space.gif\" width=$space height=7 border=0>" if $space;
			my $pepEnd=$pepBeg+$pepLength-1;
			my $infoString="<B>Sequence:</B> $pepSeq<BR>";
			if($varMod){
				$infoString.= "<B>PTMs:</B> ";
				foreach my $var (keys %varModPos){
					$infoString .= "$var ";
					foreach my $resPos (@{$varModPos{$var}}){
						my ($pos) = ($resPos =~ /^\w\:(\d+)/);
						if($pos){
							my $protPos = $pos + $pepBeg - 1;
							$infoString .= "$resPos (Prot:$protPos), ";
						} else {
							$infoString .= "$resPos, ";
						}
					}
					$infoString =~ s/\, \Z/<BR> \+ /;
				}
				$infoString =~ s/ \+ \Z//;
			}

			$infoString.= "<B>Charge:</B> $pepCharge<SUP>+</SUP><BR><B>Score:</B> $pepScore<BR><B>Length:</B> $pepLength aa<BR><B>Match position:</B> residues $pepBeg to $pepEnd";
			$infoString.="<BR><B>Occurence:</B> x$listPos{$protID}{$anaID}{$pepSeq}{$pepBeg}" if $view=~/Cum/;
			my $imgString=($varMod)? 'prot_FFFF00.gif' : 'peptide.gif';
			print "<A href=\"javascript:showPeptides('$pepSeq','$varMod')\" onmouseover=\"popup('$infoString')\" onmouseout=\"popout()\"><IMG name=\"$pepSeq$varMod\" src=\"$promsPath{images}/$imgString\" width=$length height=7 border=0></A>";
			print "<BR>\n";
		}
	}
	print "<BR>\n";
}

sub drawCompMatches {
	my ($protID,$anaID)=@_;
	unless ($usedLength{$protID}) { # no prot length & no matching peptide after filtering!
		&drawUnknownSize($DEF_PROT_LENGTH);
		print "<BR>\n";
		return;
	}

	my $scale=$largWindow/$maxProtLength;

	###<Summing match boundaries status>###
	my %boundarySum;
	my $prevRes=1;
	$boundarySum{$prevRes}=0;
	foreach my $res (sort{$a<=>$b} keys %{$compMatches{$protID}{$anaID}}) {
		$boundarySum{$res}+=$boundarySum{$prevRes}+$compMatches{$protID}{$anaID}{$res};
		$prevRes=$res;
	}

	###<Mapping colors to match status>###
	$prevRes=1;
	my $prevStatus=0;
	my $status=0;
	my $size;
	my %image=(0=>'protein.gif', 1=>'peptide.gif', 2=>'peptideBlue.gif');
	foreach my $res (sort{$a<=>$b} keys %boundarySum) {
		my $newStatus=($boundarySum{$res}==0)? 0 : ($boundarySum{$res}==1)? 1 : 2;
		next if ($newStatus && $newStatus==$status); # no need to draw
		my $beg=($status>$prevStatus)? $prevRes : $prevRes+1;
		my $end=($status>$newStatus)? $res : $res-1;
		$size=int(0.5+($end-$beg+1)*$scale);
		my $infoString=($status==0)? '<B>No match</B>' : ($status==1)? '<B>Single match:</B><BR>' : "<B>$boundarySum{$prevRes} overlapping matches:</B><BR>";
		$infoString.="From residues $beg to $end" if $status;
		print "<IMG src=\"$promsPath{images}/$image{$status}\" width=$size height=7 border=0 onmouseover=\"popup('$infoString')\" onmouseout=\"popout()\">" if $size;
		$prevStatus=$status;
		$status=$newStatus;
		$prevRes=$res;
	}
	$size=int(0.5+($usedLength{$protID}-$prevRes)*$scale);
	print "<IMG src=\"$promsPath{images}/protein.gif\" width=$size height=7 border=0 onmouseover=\"popup('<B>No match</B>')\" onmouseout=\"popout()\">" if $size;

	###<Unknow protein size>###
	&drawUnknownSize(0) unless $protLength{$protID};
	print "<BR><BR>\n";
}

sub drawUnknownSize {
	my $numTicks=($_[0])? int(($_[0]/5)+0.5) : 3;
	print "<ACRONYM onmouseover=\"popup('<B>Protein size is not known.</B>')\" onmouseout=\"popout()\">";
	foreach (1..$numTicks){
		print "<IMG src=\"$promsPath{images}/space.gif\" width=2 height=7 border=0>";
		print "<IMG src=\"$promsPath{images}/protein.gif\" width=2 height=7 border=0>";
	}
	print "</ACRONYM>";
}

####>Revision history<####
# 2.2.9 Updated due to changes in protein checkbox in calling scripts (PP 25/08/14)
# 2.2.8 Modification for ambigous positions (GA 05/09/13)
# 2.2.7 Include ghosts peptides in visualization (PP 28/08/13)
# 2.2.6 Remove VAR_MOD from script (GA 29/05/13)
# 2.2.5 Minor modification following &convertVarModString update in promsMod (GA 17/04/13)
# 2.2.4 Uses &promsMod::getProteinClass instead of getLinkClass (PP 06/02/13)
# 2.2.3 Handling ghost peptides: valid_status=0 (PP 15/07/11)
# 2.2.2 Adding PTM position on protein in peptide info popup (FY 02/05/11)
# 2.2.1 Correction bug /0 for max prot length calculation when no peptides left after filtering (PP 15/04/11)
# 2.2.0 adapting to new protein list management options + bug correction (PP 22/02/2011)
