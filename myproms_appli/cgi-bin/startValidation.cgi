#!/usr/local/bin/perl -w

################################################################################
# startValidation.cgi       2.5.0                                              #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Generates the HTML frames for the validation process                         #
# Generates the validation menu through JavaScript                             #
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

#print header; warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
#my $baseUrl="http://$ENV{HTTP_HOST}$promsPath{cgi}"; # Full path required for Macintosh compatibility
my $userID=$ENV{'REMOTE_USER'};
my $queryPageSize=&promsConfig::getQueryPageSize;
my $protPageSize=&promsConfig::getProteinPageSize;
my ($lightColor,$darkColor)=&promsConfig::getRowColors;

########################
####>>>Parameters<<<####
########################
my $analysisID=param('ID');
my $alertTax=param('alertTax');
my $selIdentifier=param('SEL');
my $PTMFilter=(param('PTMFilter'))? param('PTMFilter') : 'Indifferent';
$PTMFilter=~ s/^.:://; # fix modification
#my $okReport=(param('okReport'))? param('okReport') : 0;

##########################
####>Connecting to db<####
##########################
my $dbh=&promsConfig::dbConnect;

#############################################
####>Retrieval of variable modifications<####
#############################################
my @modifications=('Indifferent','None','Any');
#my %infoSubmit=&promsMod::getSearchParam($dbh,$analysisID);
#push @modifications,sort(split(/,/,$infoSubmit{'g:Variable modifications'})) if $infoSubmit{'g:Variable modifications'};
my %ptmList;
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


################################
####>Fetching Analysis Info<####
################################
my @itemInfo=&promsMod::getItemInfo($dbh,'analysis',$analysisID);
my $analysisName=$itemInfo[-1]{'NAME'}; $analysisName=&promsMod::shortenName($analysisName,20);
my $decoy=($itemInfo[-1]{'DECOY'})? 1 : 0;
my $projectID=$itemInfo[0]{'ID'};
my $navUrlString="$promsPath{cgi}/openProject.cgi?ID=$projectID";
$navUrlString.=($itemInfo[2]{'ITEM'} eq 'SAMPLE')? "&branchID=analysis:$analysisID&ACT=open" : "&branchID=gel2d:$itemInfo[2]{ID}&itemBranchID=analysis:$analysisID&ACT=open";
my ($validStatus,$msType,$dataFile,$taxonomy,$fileType,$instrument,$anaMaxRank)=$dbh->selectrow_array("SELECT VALID_STATUS,MS_TYPE,DATA_FILE,TAXONOMY,FILE_FORMAT,INSTRUMENT,MAX_RANK FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
$taxonomy='Unknown' unless $taxonomy;
$instrument='ESI-QUAD-TOF' unless $instrument;
$instrument =~ s/\s+\Z//;
my ($usrParam,$refFragRules,$deltaParent,$deltaFrag,$tolFrag,$nrLevel) = &promsMod::getInstrParam($dbh,$instrument,$userID);
my $dbOrganism=1;
if ($taxonomy=~/All entries/) {
	#($dbOrganism)=$dbh->selectrow_array("SELECT ORGANISM FROM DATABANK WHERE ID_DATABANK=$databankID"); # empty if multiple species
	my $sthOrg=$dbh->prepare("SELECT ORGANISM FROM ANALYSIS_DATABANK AD,DATABANK D WHERE ID_ANALYSIS=$analysisID AND AD.ID_DATABANK=D.ID_DATABANK");
	$sthOrg->execute;
	my @orgList;
	while (my ($org)=$sthOrg->fetchrow_array) {
		unless ($org) {
			$dbOrganism=0;
			last;
		}
	}
	$sthOrg->finish;
}
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
my $disableString=($projectAccess =~ /bioinfo|mass|manag|super/)? '' : 'disabled'; # have full access
my ($disabledReport,$hasQuantif);
if ($validStatus==1 && !$disableString) {
	#if ($labeling) {
	#	($hasQuantif)=$dbh->selectrow_array("SELECT COUNT(Q.ID_QUANTIFICATION) FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.FOCUS='protein' AND AQ.ID_ANALYSIS=$analysisID"); # labeled Quanti
	#}
	#else {
	#	($hasQuantif)=$dbh->selectrow_array("SELECT COUNT(*) FROM ANA_QUANTIFICATION WHERE ID_ANALYSIS=$analysisID");
	#}

	#$disabledReport=($okReport)? '' : 'disabled';
	#>Check for reported quantif data
	($hasQuantif)=$dbh->selectrow_array("SELECT 1 FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q,QUANTIFICATION_METHOD QM WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND QM.CODE !='SILAC' AND QM.CODE !='ITRAQ' AND QUANTIF_ANNOT NOT LIKE '%::SOFTWARE=PD%' AND AQ.ID_ANALYSIS=$analysisID LIMIT 0,1");
	$disabledReport=($hasQuantif)? 'disabled' : '';
}
else {
	$disabledReport=($validStatus==-1)? 'disabled' : $disableString;
}
my $disabledEnd=($validStatus==-1)? 'disabled' : $disableString;
my $aliasAccess=&promsMod::getProfileAlias($projectAccess);
#my $maxRank=($msType eq 'PMF')? 10 : &promsConfig::getMaxRank;
#$anaMaxRank=$maxRank unless ($anaMaxRank) && ($msType ne 'PMF') ;
$anaMaxRank=($msType ne 'MIS')? 10 : ($anaMaxRank)? $anaMaxRank: &promsConfig::getMaxRank;
my $maxRank=$anaMaxRank;

####>Creating the filterPTMString
my $PTMFilterString="";
foreach my $varMod (@modifications) {
	my $selMod=($varMod eq $PTMFilter)? 'selected' : '';
	$PTMFilterString.="<OPTION $selMod value=\"$varMod\">$varMod</OPTION>";
}
foreach my $modID (sort{lc($ptmList{$a}{'NAME'}) cmp lc($ptmList{$b}{'NAME'})} keys %ptmList) {
	my $selMod=($ptmList{$modID}{'NAME'} eq $PTMFilter)? 'selected' : '';
	my $extraName=($ptmList{$modID}{'MODIF_TYPE'} eq 'F') ? "$ptmList{$modID}{'SPECIFICITY'}::" : '';
	$PTMFilterString.="<OPTION value=\"$extraName$ptmList{$modID}{'NAME'}\" $selMod>$ptmList{$modID}{'NAME'}</OPTION>";
}

####>Copying data file to temp validation directory<####
#unless ($msType eq 'PMF') {
#	system ("mkdir $promsPath{tmp}") unless -e "$promsPath{tmp}";
#	system ("mkdir $promsPath{tmp}/validation") unless -e "$promsPath{tmp}/validation";
#	system ("mkdir $promsPath{tmp}/validation/ana_$analysisID") unless -e "$promsPath{tmp}/validation/ana_$analysisID";
#	if ($fileType eq "MASCOT.DAT") {
#		system ("cp $promsPath{valid}/ana_$analysisID/*.dat $promsPath{tmp}/validation/ana_$analysisID/.");
#	}
#	elsif ($fileType eq "PHENYX.XML" || $fileType eq "MASCOT.XML") {
#		system ("cp $promsPath{valid}/ana_$analysisID/*.pgf $promsPath{tmp}/validation/ana_$analysisID/.");
#		system ("cp $promsPath{valid}/ana_$analysisID/*.htm $promsPath{tmp}/validation/ana_$analysisID/."); # in case bad extension for xml file
#	}
#}

#############################
####>Fetching query data<####
#############################
my $numVerifQueries=0;
my $numValidPeptides=0;
my $numValidDecoyPeptides=-1;
my $firstQueryID;
my $selectedQueryID;
my $pepInfoStrg='';
if ($decoy) {
	foreach my $rank (1..$anaMaxRank) {$pepInfoStrg.=",INFO_PEP$rank";}
}
my $sthQ=$dbh->prepare("SELECT ID_QUERY,VALID_STATUS$pepInfoStrg FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID AND QUERY_NUM>0 ORDER BY ID_QUERY ASC"); #  same as QUERY_NUM ASC
$sthQ->execute || die $sthQ->errstr;
my $refQueryData=$sthQ->fetchall_arrayref; # reference to an array
$sthQ->finish;

my @queryPages;
my $queryCount;
my $queryPage;
my $queryPgIndex=0;
my $lowQueryID;

my $numQueries=scalar @{$refQueryData};
foreach my $refQuery (@{$refQueryData}) { # ascending query id/num
	$numVerifQueries++ if $refQuery->[1] != -1;
	#$numValidPeptides+=$refQuery->[1] if $refQuery->[1] > 0;
	$queryCount++;
	$firstQueryID=$refQuery->[0] unless $firstQueryID;
	$lowQueryID=$refQuery->[0] unless $lowQueryID;
	if ($refQuery->[1]==-1 && !$selectedQueryID) {
		$selectedQueryID=$refQuery->[0];
		$queryPage=$queryPgIndex+1;
	}
	if ($queryCount >= $queryPageSize) {
		$queryPages[$queryPgIndex]="$lowQueryID-".$refQuery->[0];
		$lowQueryID=undef;
		$queryCount=0;
		$queryPgIndex++;
	}
	# Validated pep count && manual validation check
	if ($numValidDecoyPeptides != -2) {
		foreach my $p (2..$#{$refQuery}) {
			last unless $refQuery->[$p]; # empty field
			if ($refQuery->[$p]=~/SEL=1/) {
				my ($specCount)=($refQuery->[$p]=~/SPC=(\d+)/); # top + lower-scoring peptides
				$specCount=1 unless $specCount; # old data
				#$numValidPeptides+=$specCount; # see editProjectItem.cgi where this number is stored in a FDR string...
				$numValidPeptides++;
			}
			elsif ($refQuery->[$p]=~/SEL=(2|-3)/) { # manual selection or rejection
				$numValidDecoyPeptides=-2;
				last;
			}
		}
	}
}
$queryPages[$queryPgIndex]="$lowQueryID-".$refQueryData->[$numQueries-1]->[0] if $queryCount;
unless ($selectedQueryID) { # if all checked: select first one
	$selectedQueryID=($firstQueryID)? $firstQueryID : 0;
	$queryPage=1;
}

####>Decoy<####
#my $numValidDecoyPeptides=-1;
if ($decoy) {
	###>Checking for manual validation
	#QUERY:foreach my $refQuery (@{$refQueryData}) {
	#	foreach  my $p (2..$#{$refQuery}) {
	#		last unless $refQuery->[$p]; # empty field
	#		if ($refQuery->[$p]=~/SEL=2/ || $refQuery->[$p]=~/SEL=-3/) { # maual selection or rejection
	#			$numValidDecoyPeptides=-2;
	#			last QUERY;
	#		}
	#	}
	#}
	if ($numValidDecoyPeptides>-2) { # no manual validation detected
		$numValidDecoyPeptides=0;
		my $sthDQ=$dbh->prepare("SELECT VALID_STATUS$pepInfoStrg FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID AND QUERY_NUM<0 AND VALID_STATUS>0");
		$sthDQ->execute;
		while (my ($valStatus,@pepData)=$sthDQ->fetchrow_array) {
			#$numValidDecoyPeptides+=$valStatus;
			foreach my $pepInfo (@pepData) {
				last unless $pepInfo; # empty field
				if ($pepInfo=~/SEL=1/) {
					my ($specCount)=($pepInfo=~/SPC=(\d+)/); # top + lower-scoring peptides
					$specCount=1 unless $specCount; # old data
					#$numValidDecoyPeptides+=$specCount; # includes lower-scoring
					$numValidDecoyPeptides++;
				}
			}
		}
		$sthDQ->finish;
	}
}

################################
####>Fetching proteins data<####
################################
my $numValidProteins=0;
my $filter=0;
my $showFiltered=0;
my $selectedProtID=0;
my $selectableProteins=0;
# my ($selectableProteins)=$dbh->selectrow_array("SELECT count(*) FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$analysisID AND SEL_STATUS>-3");
my $identStrg=($selIdentifier)? ',IDENTIFIER' : '';
my $sthP=$dbh->prepare("SELECT ID_PROT_VALID,SEL_STATUS,MATCH_GROUP$identStrg FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$analysisID AND IDENTIFIER NOT LIKE 'DECOY_%' ORDER BY MAX_MATCH DESC,MAX_SCORE DESC,ID_PROT_VALID ASC");
$sthP->execute || die $sthP->errstr;
my $refProtData=$sthP->fetchall_arrayref; # reference to an array
$sthP->finish;
$dbh->disconnect;

####>General protein info<####
my %matchGroups;
my $selectedMGr=0;
foreach my $refProt (@{$refProtData}) {
	$matchGroups{$refProt->[2]}++;
	$selectableProteins++ if $refProt->[1] > -3;
	$numValidProteins++ if $refProt->[1] >= 1;
	$filter=1 if $refProt->[1]==-3;
	if ($selIdentifier && $refProt->[3] eq $selIdentifier) {
		$selectedProtID=$refProt->[0];
		$selectedMGr=$refProt->[2];
		$showFiltered=1 if $refProt->[1] <= -2; # excluded prot will be visible
	}
}

####>Splitting proteins in multiple sets<####
my @proteinPages;
my $protPage; # 1...n
my $protCount=0;
my $protPgIndex=0; # 0..n-1
my @pageMGr;
my %usedMGr;
#foreach my $refProt (@{$refProtData}) {
#	next if $usedMGr{$refProt->[2]};
#	$usedMGr{$refProt->[2]}=1;
#	$protPage=$protPgIndex+1 if $refProt->[2]==$selectedMGr;
#	$protCount+=$matchGroups{$refProt->[2]};
#	$pageMGr[$protPgIndex]{$refProt->[2]}=1;
#	if ($protCount >= $protPageSize) {
#		$proteinPages[$protPgIndex]=&groupsToString($pageMGr[$protPgIndex]);
#		$protCount=0;
#		$protPgIndex++;
#	}
#}
foreach my $matchGr (sort{$a<=>$b} keys %matchGroups) {
	$protPage=$protPgIndex+1 if $matchGr==$selectedMGr;
	$pageMGr[$protPgIndex]{$matchGr}=1;
	$protCount+=$matchGroups{$matchGr};
	if ($protCount >= $protPageSize) {
		$proteinPages[$protPgIndex]=&groupsToString($pageMGr[$protPgIndex]);
		$protCount=0;
		$protPgIndex++;
	}
}
$proteinPages[$protPgIndex]=&groupsToString($pageMGr[$protPgIndex]) if $protCount;

###>Scanning protein pages for 1st selectable protein<### (no pre-selected identifier provided)
unless ($protPage) {
	my ($firstProtID,$firstProtPage);
	PAGE:foreach my $pgIndex (0..$#pageMGr) {
		foreach my $refProt (@{$refProtData}) {
			if (!$firstProtID && $refProt->[1] > -3) {
				$firstProtID=$refProt->[0];
				$firstProtPage=$pgIndex+1;
			}
			next unless $pageMGr[$pgIndex]{$refProt->[2]}; # restrict to selected MGr
			if ($refProt->[1]==-1 && !$selectedProtID) {
				$selectedProtID=$refProt->[0];
				$protPage=$pgIndex+1;
				last PAGE;
			}
		}
	}
	unless ($selectedProtID) { # if all checked: select first one
		if ($firstProtID) {
			$selectedProtID=$firstProtID;
			$protPage=$firstProtPage;
		}
		else { # !!!NO PROTEIN AT ALL!!!
			$selectedProtID=0;
			$protPage=1;
		}
	}
}

#######################
####>Starting HTML<####
#######################
print header(-'content-encoding'=>'no',-charset=>'UTF-8');
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Main Query Validation Window</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
<!--  to hide script contents from old browsers
|;
####>Storing proteins/query sets in parent frame (top)<####
print qq
|// Page management variables
var proteinPages=new Array();
var selProteinPage=$protPage;
var queryPages=new Array();
var selQueryPage=$queryPage;
var ptmFilter='$PTMFilter';
var ptmInitPage=1;
|;
foreach my $i (0..$#proteinPages) {
	print "proteinPages[$i]='$proteinPages[$i]';\n";
}
print "proteinPages[0]='0';\n" if $#proteinPages==-1;
foreach my $i (0..$#queryPages) {
	print "queryPages[$i]='$queryPages[$i]';\n";
}
print "queryPages[0]='0';\n" if $#queryPages==-1;
print qq
|function changeFilter (filterBox) {
	if (filterBox.checked) {showFiltered=1;} else {showFiltered=0;}
	var itemId;
	if (selectedView=='query') {updateValidation(selectedQueryId,0,0,0,0,selectedView,selQueryPage,1);}
	else {updateValidation(selectedProtId,0,0,0,0,selectedView,selProteinPage,1);} // protein, myList
}
function updateValidation(itemId,newVerifQuery,newValidPeptide,newValidDecoyPeptide,newValidProteins,view,page,updateRank) {
	numVerifQueries+=newVerifQuery;
	numValidPeptides+=newValidPeptide;
	if (newValidDecoyPeptide<0) {numValidDecoyPeptides=-2;} // manual validation has occured
	else if (numValidDecoyPeptides>=0) {
		numValidDecoyPeptides+=newValidDecoyPeptide;
	}
	numValidProteins+=newValidProteins;
	selectedView=view;
	if (selectedView=='query') {
		selectedQueryId=itemId;
		selQueryPage=page;
	}
	else { // protein, myList
		selectedProtId=itemId;
		selProteinPage=page;
	}
	writeMenu(itemId,updateRank); // 1 => update rankFrame
}
function writeMenu(selItem,updateRank) {
	var selQuery=' '; var selProt=' '; var selMyList=' ';
	var newQueryFrameURL;
	if (selectedView=='query') {
		selQuery=' selected ';
		selectedSort=selQuerySort;
		var pageStrg=selQueryPage+'-'+queryPages.length;
		newQueryFrameURL="$promsPath{cgi}/listQueries.cgi?ID=$analysisID&MSTYPE=$msType&SEL="+selItem+"&SORT="+selectedSort+"&UPRANK="+updateRank+"&PAGE="+pageStrg+"&QR="+queryPages[selQueryPage-1];
	}
	else {
		selectedSort=selProtSort;
		var pageStrg=selProteinPage+'-'+proteinPages.length;
		if (selectedView=='protein') {
			selProt=' selected ';
			newQueryFrameURL="$promsPath{cgi}/listTempProteins.cgi?ID=$analysisID&MSTYPE=$msType&SEL="+selItem+"&SORT="+selectedSort+"&UPRANK="+updateRank+"&SHOWFILT="+showFiltered+"&PAGE="+pageStrg+"&MG="+proteinPages[selProteinPage-1]+"&DISPLAYPTM="+displayPTM+"&PTMFILTER="+ptmFilter;
		}
		else { // myList
			selMyList=' selected ';
			newQueryFrameURL="$promsPath{cgi}/listMyTempProteins.cgi?ACT=list&ID=$analysisID&MSTYPE=$msType&SEL="+selItem+"&SORT="+selectedSort+"&UPRANK="+updateRank+"&SHOWFILT="+showFiltered+"&PAGE="+pageStrg;
		}
	}

	//Updating queryFrame (1st because takes longer to be loaded than MenuFrame)
	if (newQueryFrameURL==curQueryFrameURL) {
		queryFrame.window.document.open(); // clear queryFrame for proper update
		queryFrame.window.document.close();
	}
	else {curQueryFrameURL=newQueryFrameURL;}
	queryFrame.location=newQueryFrameURL;

	//Updating menuFrame
	var doc=menuFrame.window.document;
	var queryString=(numVerifQueries>1)? 'queries' : 'query';
	var protString=(numValidProteins>1)? 'proteins' : 'protein';
	var pepString=(numValidPeptides>1)? 'peptides' : 'peptide';
	var FDRstring;
	if (numValidDecoyPeptides==-2) {FDRstring='N/A <FONT style="font-size:10px;">(Manual validation)</FONT>';}
	else if (numValidDecoyPeptides==-1) {FDRstring='N/A <FONT style="font-size:10px;">(No decoy data)</FONT>';}
	else {
		if (numValidPeptides>0) { // check for division by 0
			var FDRpc=100*numValidDecoyPeptides/numValidPeptides;
			FDRpc=Math.round(FDRpc*100)/100;
		}
		else {FDRpc=0;}
		var decoyPepString=(numValidDecoyPeptides>1)? 'peptides' : 'peptide';
		FDRstring=FDRpc+'% <FONT style="font-size:10px;">('+numValidDecoyPeptides+' decoy '+decoyPepString+')</FONT>';
	}
	doc.open();
	doc.writeln('<HTML><HEAD>');
	doc.writeln('<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">');
|;
my $disabRestoreStrg=($validStatus==1)? $disableString : 'disabled';
print qq
|	doc.writeln('</HEAD><BODY topmargin="0" leftmargin="0"><CENTER>');
	doc.writeln('<TABLE cellpadding=2 cellspacing=2>');
	doc.writeln('<TR><TH align=left nowrap colspan=2><FONT style="font-size:15px;color:#DD0000">$analysisName</FONT><FONT style="font-size:10px;">&nbsp;($dataFile)</FONT></TH></TR>');

	doc.writeln('<TR><TH colspan=2 nowrap>');
/*
*	doc.write('<BUTTON style="width:100px;height:35px;" onclick="parent.autoSelection()" $disableString><FONT style="line-height:13px;">Automated<BR>Selection</FONT></BUTTON>');
*	doc.write('<BUTTON style="width:100px;height:35px;" onclick="parent.clearSelection()" $disableString><FONT style="line-height:13px;">Clear<BR>Selection</BUTTON><BR>');
*	doc.write('<BUTTON style="width:100px;height:35px;" onclick="parent.restoreSelection()" $disabRestoreStrg><FONT style="line-height:13px;font-weight:bold;color:#00AA00;">Restore<BR>Selection</FONT></BUTTON>');
*	doc.write('<BUTTON style="width:100px;height:35px;" onclick="parent.location='+urlString+'"><FONT style="line-height:13px;font-weight:bold;">Exit<BR>Validation</FONT></BUTTON><BR>');
*	doc.write('<BUTTON style="width:100px;height:35px;" onclick="parent.send2Biologist(0);" $disabledReport><FONT style="line-height:13px;font-weight:bold;color:#0000DD;">Report<BR>Validation</FONT></BUTTON>');
*	doc.write('<BUTTON style="width:100px;height:35px;" onclick="parent.endValidation();" $disabledEnd><FONT style="line-height:13px;font-weight:bold;color:#DD0000;">End<BR>Validation</FONT></BUTTON>');
*/
	doc.write('<BUTTON style="width:180px;height:22px" onclick="parent.restoreSelection()" $disabRestoreStrg><FONT style="line-height:13px;font-weight:bold">Restore selection</FONT></BUTTON><BR>');
	doc.write('<BUTTON style="width:180px;height:22px" onclick="parent.send2Biologist(0);" $disabledReport><FONT style="line-height:13px;font-weight:bold;">Report validation</FONT></BUTTON><BR>');
	doc.write('<BUTTON style="width:180px;height:22px" onclick="parent.location='+urlString+'"><FONT style="line-height:13px;font-weight:bold;">Exit validation</FONT></BUTTON>');
	doc.writeln('</TH></TR>');

	doc.writeln('<TR><TH colspan=2 nowrap><TABLE cellpadding=0 cellspacing=0>');
		doc.writeln('<TR><TH nowrap align=right>'+numVerifQueries+'/$numQueries</TH><TH align=left nowrap>&nbsp;verified '+queryString+'.</TH></TR>');
		doc.writeln('<TR><TH nowrap align=right>'+numValidPeptides+'</TH><TH align=left nowrap>&nbsp;validated '+pepString+'.</TH></TR>');
		doc.writeln('<TR><TH nowrap align=right>'+numValidProteins+'/'+selectableProteins+'</TH><TH align=left nowrap>&nbsp;selected '+protString+'.</TH></TR>');
		doc.writeln('<TR><TH nowrap align=left colspan=2>FDR: '+FDRstring+'</TH></TR>');
	doc.writeln('</TABLE></TD></TR>');
	doc.writeln('</TABLE>');
	doc.writeln('<FONT style="font-size:10px;"><HR width=90%><FONT>');
	doc.writeln('<TABLE cellpadding=0 cellspacing=0>');
	doc.writeln('<TR><TH nowrap align=right>View:</TH>');
	doc.writeln('<TD><SELECT name="selectView" id="selectView" style="font-weight:bold;" onchange="parent.setView(this.value)">');
	doc.writeln('<OPTION'+selQuery+'value="query">ALL QUERIES</OPTION>');
	doc.writeln('<OPTION'+selProt+'value="protein">ALL PROTEINS</OPTION>');
	doc.writeln('<OPTION'+selMyList+'value="myList">MY PROTEINS</OPTION>');
	doc.writeln('</SELECT></TD></TR>');
	doc.writeln('<TR><TH nowrap align=right>Sort by:</TH><TD><SELECT name="selectSort" style="font-weight:bold;font-size:12px" onchange="parent.setSort(this)">');
	var i;
	var selString;
	if (selectedView=='query') {
		for (i=0;i<querySortItem.length;i++) {
			if (querySortItem[i]==selectedSort) {selString=' selected ';}
			else {selString=' ';}
			doc.writeln('<OPTION'+selString+'value="'+querySortItem[i]+'">'+querySortText[i]+'</OPTION>');
		}
	}
	else { //protein or myList
		for (i=0;i<protSortItem.length;i++) {
			if (protSortItem[i]==selectedSort) {selString=' selected ';}
			else {selString=' ';}
			doc.writeln('<OPTION'+selString+'value="'+protSortItem[i]+'">'+protSortText[i]+'</OPTION>');
		}
	}
	doc.writeln('</SELECT></TD></TR>');


	doc.writeln('<TR><TH nowrap align=right>PTM filter:</TH>');
	doc.writeln('<TD><SELECT name="PTMFilter" style="width:138px;font-weight:bold;font-size:12px" onchange="parent.updatePTMFilter(this.value)">');
	doc.writeln('$PTMFilterString');
	doc.writeln('</SELECT></TD></TR>');

	if ($filter==1) {
		doc.writeln('<TR><TH nowrap align=left colspan=2>&nbsp<INPUT type="checkbox" name="checkFilter"');
		if (showFiltered) {doc.writeln(' checked ');}
		doc.writeln('onclick="parent.changeFilter(this)"/>&nbsp;Show filtered proteins.</TD></TR>');
	}

	doc.writeln('</TABLE></CENTER></BODY></HTML>');
	doc.close();
}
function setView(newView) { // query,protein,myList
	var prevView=selectedView;
	selectedView=newView;
	if (selectedView=='query') {
		writeMenu(selectedQueryId,1);
	}
	else { // protein, myList
		var updateRank=(prevView=='myList' && queryFrame.selectedProtId)? 0 : (prevView=='query' && selectedView=='myList')? 2 : 1; // 2: forces listMyTempProteins.cgi to update
		writeMenu(selectedProtId,updateRank);
	}
}
function setSort(selectFormSort) {
	selectedSort=selectFormSort.value;
	if (selectedView=='query') {
		selQuerySort=selectFormSort.value;
		var pageStrg=selQueryPage+'-'+queryPages.length;
		queryFrame.location="$promsPath{cgi}/listQueries.cgi?ID=$analysisID&MSTYPE=$msType&SEL="+selectedQueryId+"&SORT="+selectedSort+"&PAGE="+pageStrg+"&QR="+queryPages[selQueryPage-1];;
	}
	else if (selectedView=='protein') {
		selProtSort=selectFormSort.value;
		var pageStrg=selProteinPage+'-'+proteinPages.length;
		queryFrame.location="$promsPath{cgi}/listTempProteins.cgi?ID=$analysisID&MSTYPE=$msType&SEL="+selectedProtId+"&SORT="+selectedSort+"&SHOWFILT="+showFiltered+"&PAGE="+pageStrg+"&MG="+proteinPages[selProteinPage-1]+"&DISPLAYPTM="+displayPTM+"&PTMFILTER="+ptmFilter;
	}
	else { // myList
		selProtSort=selectFormSort.value;
		queryFrame.location="$promsPath{cgi}/listMyTempProteins.cgi?ID=$analysisID&MSTYPE=$msType&SEL="+selectedProtId+"&SORT="+selectedSort+"&SHOWFILT="+showFiltered;
	}
}
function endValidation() {
	if (!confirm("This will terminate Validation.   \\n\\n         Proceed ?")) {return;}
	var sendData=0;
	if (confirm("Send data to biologists first ?")) {
		sendData=send2Biologist(1);
		if (sendData==-1) {return;} // cancelled from send2Biologist()
	}
	else if (!confirm("End Validation anyway ?")) {return;}
	var endVal=(sendData)? 2 : 1;
	window.location="$promsPath{cgi}/send2Biologist.cgi?ID=$analysisID&MSTYPE=$msType&endVal="+endVal;
}
function send2Biologist(endVal) {
	if (numValidProteins==0 && !confirm("There are no proteins selected for this Analysis !   \\n\\n                          Send anyway ?")) {
		if (endVal && !confirm("End Validation anyway ?")) {return -1;}
		else {return 0;}
	}
	var proceed=1;
	var notVerified = $numQueries - numVerifQueries;
	var queryText=(notVerified==1)? ' query has' : ' queries have';
	if (notVerified > 0 && !confirm(notVerified + queryText + " not been verified.        \\n\\n            Send anyway ?")) {
		if (endVal && !confirm("End Validation anyway ?")) {return -1;}
		else {proceed=0;}
	}
	if (endVal) {return proceed;}
	else if (proceed==1) {window.location="$promsPath{cgi}/send2Biologist.cgi?ID=$analysisID&MSTYPE=$msType&endVal=0";}
}
function autoSelection() {
	rankFrame.location="$promsPath{cgi}/autoSelect.cgi?id_project=$projectID&ITEM=ANALYSIS&ID=$analysisID&MSTYPE=$msType&ACT=select";
}
function clearSelection() {
	rankFrame.location="$promsPath{cgi}/autoSelect.cgi?id_project=$projectID&ITEM=ANALYSIS&ID=$analysisID&MSTYPE=$msType&ACT=clearSelection";
}
function restoreSelection() {
	if (confirm('Clear current selection and restore previously transferred one ?')) {
		rankFrame.location="$promsPath{cgi}/autoSelect.cgi?id_project=$projectID&ITEM=ANALYSIS&ID=$analysisID&MSTYPE=$msType&ACT=restore";
	}
}
function updatePTMFilter(newVal){
	if (ptmFilter != newVal){
		ptmFilter = newVal;
		ptmInitPage=selProteinPage;
		if (selectedView=='protein') {
			var pageStrg=selProteinPage+'-'+proteinPages.length;
			queryFrame.location="$promsPath{cgi}/listTempProteins.cgi?ID=$analysisID&MSTYPE=$msType&SHOWFILT="+showFiltered+"&SEL="+selectedProtId+"&SORT="+selectedSort+"&PAGE="+pageStrg+"&MG="+proteinPages[selProteinPage-1]+"&DISPLAYPTM="+displayPTM+"&PTMFILTER="+ptmFilter;
		}
	}
}
|;
if ($alertTax && $taxonomy=~/All entries/ && !$dbOrganism) {
	print "alert(\"WARNING: The taxonomy for this Analysis is '$taxonomy'.\\nPeptides might match several orthologs in different species.\");\n";
}
if ($hasQuantif) {
	print "alert(\"WARNING: This Analysis is associated with quantification data.\\nThese data must be deleted before new Reporting.\");\n";
}
my $menuRow=($filter)? 270 : 250; # 300: 280; #
my $rkFrameSize=($msType ne 'PMF')? 170 : 300;
print qq
|var urlString="'$navUrlString'";
var curQueryFrameURL='';
var maxDisp=$anaMaxRank;
var maxSrch=$maxRank;
top.spectrumTolerance=$tolFrag;
//top.spectrumTolUnit='$deltaFrag'; NOT USED
top.spectrumMinInt=$nrLevel;
var listRankSort='score';
var listRankDlt='$deltaParent';
var varValue='Indifferent';//use in listRanks.cgi
var multiAnaScan=0;//use in listRanks.cgi
var listMatchedAnalysesSort='name';
var numVerifQueries=$numVerifQueries;
var numValidPeptides=$numValidPeptides;
var numValidDecoyPeptides=$numValidDecoyPeptides;
var numValidProteins=$numValidProteins;
var selectableProteins=$selectableProteins;
var selectedQueryId=$selectedQueryID;
var selectedProtId=$selectedProtID;
var showFiltered=$showFiltered;
var displayPTM=0;
selectedView='protein'; // if changed to 'query', change selectedSort and <FRAMESET onload...> accordingly !!!
var querySortItem=['MAX_SCORE:DESC','MAX_SCORE:ASC','QUERY_NUM:ASC','QUERY_NUM:DESC','VALID_STATUS:DESC','VALID_STATUS:ASC'];
var querySortText=['max. score desc.','max. score asc.','query number asc.','query number desc.','verified query first','verified query last'];
var protSortItem=['MATCH_GROUP:ASC,MAX_SCORE:DESC','MAX_SCORE:DESC','MAX_SCORE:ASC','SCORE:DESC','SCORE:ASC','MAX_MATCH:DESC','MAX_MATCH:ASC','NUM_MATCH:DESC','NUM_MATCH:ASC','SEL_STATUS:DESC','SEL_STATUS:ASC','IDENTIFIER:ASC'];
var protSortText=['match group','max. score desc.','max. score asc','score des.','score asc.','max. matches desc','max. matches asc.','valid matches desc.','valid matches asc.','selected first','selected last','identifier'];
var selQuerySort='QUERY_NUM:ASC';
var selProtSort='MATCH_GROUP:ASC,MAX_SCORE:DESC';
var selectedSort=selProtSort; // if changed to selQuerySort: change selectedView and <FRAMESET onload...> accordingly !!!
//peptideView.cgi variables
//var peptideViewTolerance = 0.4 ;
// end hiding contents from old browsers  -->
</SCRIPT>
</HEAD>

<FRAMESET cols="270,*" bordercolor="$darkColor" onload="writeMenu($selectedProtID,1)"><!--$selectedQueryID-->
	<FRAMESET rows="$menuRow,*" border=0>
		<FRAME name="menuFrame" scrolling=no>
		<FRAME name="queryFrame">
	</FRAMESET>
	<FRAMESET rows="$rkFrameSize,*" bordercolor="$darkColor">
		<FRAME name="rankFrame">
		<FRAME name="spectrumFrame">
	</FRAMESET>
</FRAMESET>

</HTML>
|;
############################<< SUBROUTINES >>###############################


###############################################
####<Converting match group list to string>####
###############################################
sub groupsToString {
	my $refMGr=$_[0];
	my ($mgString,$startGr,$prevGr,$currGr);
	my $firstLoop=1;
	foreach my $matchGr (sort{$a<=>$b} keys %{$refMGr}) {
		$currGr=$matchGr;
		if ($firstLoop) {
			$mgString="$matchGr";
			$startGr=$prevGr=$matchGr;
			$firstLoop=0;
			next;
		}
		if ($matchGr==$prevGr+1 && $startGr==$prevGr) {
			$mgString.='-';
		}
		elsif ($matchGr > $prevGr+1) {
			$mgString.="$prevGr" unless $prevGr==$startGr;
			$mgString.=",$matchGr";
			$startGr=$matchGr;
		}
		$prevGr=$matchGr;
	}
	$mgString.="$currGr" if $startGr < $currGr;
	return $mgString;
}

####>Revision history<####
# 2.5.0 Update PTM filter so as to do the same as listRanks.cgi (GA 17/05/17)
# 2.4.9 Fix typo on DISPLAYPTM argument (PP 16/06/16)
# 2.4.8 Change in main action buttons & added FRAMESET bordercolor (PP 24/03/16)
# 2.4.7 Change $numValidPeptides and $numValidDecoyPeptides count to avoid over-estimation of peptide number (GA 06/08/15)
# 2.4.6 Allow new report on Proteome Discoverer label-free quantification (PP 22/04/15)
# 2.4.5 Does not rely on selectOption okReport param (PP 06/03/15)
# 2.4.4 Check for uninitialiazed $taxonomy (PP 18/12/14)
# 2.4.3 Corrected decoy/target count (PP 07/05/14)
# 2.4.2 Added multiAnaScan as listRanks global parameter (PP 26/09/13)
# 2.4.1 Minor bug fix when no variable modification (PP 29/04/13)
# 2.4.0 Multi-databank search (PP 18/12/12)
# 2.3.9 Larger frame for protein list (PP 11/12/12)
# 2.3.8 Simplified check for allowing reporting (PP 07/05/12)
# 2.3.7 No reporting if exist quantification data (PP 31/10/11)
# 2.3.6 Manager has full access validation data (PP 01/10/11)
# 2.3.5 Compatibility with mixed (SQ) search (PP 16/04/2011)
# 2.3.4 No protein bug correction (FY 29/03/11)
# 2.3.3 Add DISPLAYPTM in listTempProteins.cgi call (GA xx/01/11)<BR>See 2.1.0 of listTempProteins.cgi
