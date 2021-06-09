#!/usr/local/bin/perl -w

################################################################################
# editBranch.cgi                     2.4.2                                     #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Prepares deletion of any item and its contents -> deleteProjectItem.cgi      #
# Moves any item (except project) and its contents                             #
# Ends or archives project                                                     #
# Not used to specifically delete analyses                                     #
# Also ends, archives or re-opens a project                                    #
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
use POSIX qw(strftime); # to get the time
use Archive::Tar; # compression requires IO::Zlib!!!
use File::Path qw(rmtree); #remove_tree
use File::Copy;
use File::Copy::Recursive qw(dirmove);
use File::Find qw(find finddepth);
no warnings 'File::Find';

#print header; warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my $date=strftime("%Y-%m-%d %H:%M:%S", localtime);

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

###################
####>Arguments<####
###################
my $item=&promsMod::cleanParameters(param('ITEM'));
my $itemID=&promsMod::cleanNumericalParameters(param('ID'));
my $action=param('ACT');
if ($action=~/utoEndVal/) {&setAutoEndValidation; exit;} # for project only
elsif ($action eq 'continue') {&continueProject; exit;} # for project only
elsif ($action=~/^(end|archive)$/) {&endOrArchiveProject; exit;} # for project only
elsif ($action eq 'error') {&writeError; exit;} # if ajax call to send2Biologist.cgi fails
#my $delete=(param('GO'));
my $move=(param('move'))? 1 : 0; # defined after form submission
my $newParentBranchID=(param('newBranchID'))? param('newBranchID') : 0;
my $duplClassif=(param('duplClassif'))?1: 0;

my @itemInfo=&promsMod::getItemInfo($dbh,$item,$itemID);
my $itemName=$itemInfo[-1]{NAME};
my $itemType=$itemInfo[-1]{'TYPE'};
my $projectID=$itemInfo[0]{'ID'};
my $parentID=$itemInfo[-2]{'ID'}; # -2 not possible for project!
my $parentItem=$itemInfo[-2]{'ITEM'};
my $parentType=($item eq 'ANALYSIS')? 'Sample or Spot' : $itemInfo[-2]{'TYPE'};
my $actionString=($action eq 'delete')? 'Deleting' : 'Moving';

#######################
####>Starting HTML<####
#######################
print header(-charset=>'UTF-8');
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Edit Branch</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE>
</STYLE>
<SCRIPT type="text/javascript">
|;
if ($action eq 'move') {
	print qq
|function checkForm(myForm) {
	if (!myForm.newBranchID.value) { // no Parent selected
		alert('You must select a new $parentType for this $itemType.');
		return false;
	}
	return true;
}
|;
}
print qq
|function cancelAction() {
	top.promsFrame.selectedAction='summary'; // set default action to 'summary'
	top.promsFrame.optionFrame.selectOption();
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">$actionString $itemType <FONT color="#DD0000">$itemName</FONT> and its Contents</FONT>
<BR><BR><BR>
|;


#########################################
####>Storing info on selected Branch<####
#########################################
my $sthG2D=$dbh->prepare("SELECT ID_GEL2D FROM GEL2D WHERE ID_EXPERIMENT=?");
my $sthSpot=$dbh->prepare("SELECT ID_SPOT FROM SPOT WHERE ID_GEL2D=?");
my $sthSamp1=$dbh->prepare("SELECT ID_SAMPLE FROM SAMPLE WHERE ID_EXPERIMENT=? AND ID_SPOT IS NULL");
my $sthSamp2=$dbh->prepare("SELECT ID_SAMPLE FROM SAMPLE WHERE ID_SPOT=?");
my $sthAna=$dbh->prepare("SELECT ID_ANALYSIS,VALID_STATUS FROM ANALYSIS WHERE ID_SAMPLE=?");
my $sthProt=$dbh->prepare("SELECT ID_PROTEIN FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=?");
my $sthQuantif=$dbh->prepare("SELECT A.ID_QUANTIFICATION,Q.ID_DESIGN FROM ANA_QUANTIFICATION A,QUANTIFICATION Q WHERE A.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND A.ID_ANALYSIS=?");
my $sthBioSampObs=$dbh->prepare("SELECT 1 FROM OBSERVATION WHERE ID_ANALYSIS=? AND ID_BIOSAMPLE IS NOT NULL");

my (%itemList,%movedValidAnalyses,%proteinList,%quantifList);
my ($numNotImp,$numUnValid,$numPartValid,$numValid,$numGOAna,$existDesignQuantif,$numExplorAna,$numPathwayAna,$numGseaAna,$linkedBioSamp);
if ($item eq 'EXPERIMENT') {
	@{$itemList{'EXPERIMENT'}}=($itemID);
	&getGel2dList($itemID);
	&getSampList($sthSamp1,$itemID);
	if (scalar keys %proteinList) {
		($numGOAna)=$dbh->selectrow_array("SELECT COUNT(*) FROM GO_ANALYSIS WHERE ID_EXPERIMENT=$itemID");
		($numExplorAna)=$dbh->selectrow_array("SELECT COUNT(*) FROM EXPLORANALYSIS WHERE ID_EXPERIMENT=$itemID");
		($numPathwayAna)=$dbh->selectrow_array("SELECT COUNT(*) FROM PATHWAY_ANALYSIS WHERE ID_EXPERIMENT=$itemID AND ANALYSIS_TYPE='PATHWAY'");
		($numGseaAna)=$dbh->selectrow_array("SELECT COUNT(*) FROM PATHWAY_ANALYSIS WHERE ID_EXPERIMENT=$itemID AND ANALYSIS_TYPE='GSEA'");
	}
}
elsif ($item eq 'GEL2D') {
	@{$itemList{'GEL2D'}}=($itemID);
	&getSpotList($itemID);
}
elsif ($item eq 'SPOT') { # for delete only
	@{$itemList{'SPOT'}}=($itemID);
	my ($sampID)=$dbh->selectrow_array("SELECT ID_SAMPLE  FROM SAMPLE WHERE ID_SPOT=$itemID");
	@{$itemList{'SAMPLE'}}=($sampID);
	&getAnaList($sampID);
}
elsif ($item eq 'SAMPLE') {
	@{$itemList{'SAMPLE'}}=($itemID);
	&getAnaList($itemID);
}
elsif ($item eq 'ANALYSIS') {
	@{$itemList{'ANALYSIS'}}=($itemID);
	my ($valStatus)=$dbh->selectrow_array("SELECT VALID_STATUS FROM ANALYSIS WHERE ID_ANALYSIS=$itemID");
	if ($valStatus>=1) {
		&getProtList($itemID);
		$movedValidAnalyses{$itemID}=$valStatus;
	}
}
$sthG2D->finish;
$sthSpot->finish;
$sthSamp1->finish;
$sthSamp2->finish;
$sthAna->finish;
$sthProt->finish;
$sthQuantif->finish;
$sthBioSampObs->finish;

#my %anaInfo; # global, used for deletion
#my ($numSamp,$numValid,$numPartValid,$numUnValid,$numNotImp,$numValProt)=(0,0,0,0,0);
#my $sthSelSamp=$dbh->prepare("SELECT ID_SAMPLE FROM SAMPLE WHERE ID_EXPERIMENT=?");
#my $sthSelAna;
#if ($item eq 'ANALYSIS'){
#	$sthSelAna=$dbh->prepare("SELECT VALID_STATUS,DATA_FILE FROM ANALYSIS WHERE ID_ANALYSIS=?");
#}
#else {
#	$sthSelAna=$dbh->prepare("SELECT ID_ANALYSIS,VALID_STATUS,DATA_FILE FROM ANALYSIS WHERE ID_SAMPLE=?");
#}
#if ($item eq 'EXPERIMENT') {
#	$sthSelSamp->execute($itemID);
#	while (my ($sampID)=$sthSelSamp->fetchrow_array) {
#		%{$anaInfo{$sampID}}=(); # in case no analyses in sample
#		$sthSelAna->execute($sampID);
#		while (my ($anaID,$validStatus,$dataFile)=$sthSelAna->fetchrow_array) {
#			@{$anaInfo{$sampID}{$anaID}}=($validStatus,$dataFile);
#			if ($validStatus==2) {$numValid++;} elsif ($validStatus==1) {$numPartValid++;} elsif ($validStatus==0) {$numUnValid++;} else {$numNotImp++;}
#		}
#		$numSamp++;
#	}
#	if ($numValid || $numPartValid) {
#		($numValProt)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM ANALYSIS_PROTEIN,ANALYSIS,SAMPLE WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_EXPERIMENT=$itemID");
#	}
#}
#elsif ($item eq 'SAMPLE') {
#	$sthSelAna->execute($itemID);
#	while (my ($anaID,$validStatus,$dataFile)=$sthSelAna->fetchrow_array) {
#		@{$anaInfo{$itemID}{$anaID}}=($validStatus,$dataFile);
#		if ($validStatus==2) {$numValid++;} elsif ($validStatus==1) {$numPartValid++;} elsif ($validStatus==0) {$numUnValid++;} else {$numNotImp++;}
#	}
#	$numSamp=1;
#	if ($numValid || $numPartValid) {
#		($numValProt)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM ANALYSIS_PROTEIN,ANALYSIS WHERE ANALYSIS_PROTEIN.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ID_SAMPLE=$itemID");
#	}
#}
#else { # ANALYSIS
#	$sthSelAna->execute($itemID);
#	my ($validStatus,$dataFile)=$sthSelAna->fetchrow_array;
#	@{$anaInfo{0}{$itemID}}=($validStatus,$dataFile);
#	if ($validStatus>=1) {
#		if ($validStatus==2) {$numValid++;} else {$numPartValid++;}
#		($numValProt)=$dbh->selectrow_array("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=$itemID");
#	} elsif ($validStatus==0) {$numUnValid++;} else {$numNotImp++;}
#}
#$sthSelSamp->finish;
#$sthSelAna->finish;


####>Listing items to be affected<####
if ($action eq 'delete') {
	print "<TABLE>\n<TR><TD class=\"title2\" colspan=3><FONT color=\"#DD0000\">WARNING!<BR>The following items will be permanently deleted:</FONT></TD></TR>\n";
}
else {
	print "<TABLE>\n<TR><TD class=\"title2\" colspan=3>This action will affect:</TD></TR>\n";
}
if ($item eq 'EXPERIMENT') {
	print "<TR><TD width=20></TD><TD class=\"title3\" align=right>1</TD><TD class=\"title3\"> Experiment</TD></TR>\n";
}
foreach my $it ('GEL2D','SPOT','SAMPLE') {
	next unless $itemList{$it};
	my $numIt=scalar @{$itemList{$it}};
	my $itemStrg=($numIt>1)? 's' : '';
	print "<TR><TD width=20></TD><TD class=\"title3\" align=right>$numIt</TD><TD class=\"title3\"> ",&promsMod::getItemType($it),"$itemStrg</TD></TR>\n";
}
if ($numNotImp) {
	my $anaStrg=($numNotImp>1)? 'es' : 'is';
	print "<TR><TD width=20></TD><TD class=\"title3\" align=right>$numNotImp</TD><TD class=\"title3\"> non-imported MS/MS Analys$anaStrg</TD></TR>\n";
}
if ($numUnValid) {
	my $anaStrg=($numUnValid>1)? 'es' : 'is';
	print "<TR><TD width=20></TD><TD class=\"title3\" align=right>$numUnValid</TD><TD class=\"title3\"> non-validated MS/MS Analys$anaStrg</TD></TR>\n";
}
if ($numPartValid) {
	my $anaStrg=($numPartValid>1)? 'es' : 'is';
	print "<TR><TD width=20></TD><TD class=\"title3\" align=right>$numPartValid</TD><TD class=\"title3\"> partially validated MS/MS Analys$anaStrg</TD></TR>\n";
}
if ($numValid) {
	my $anaStrg=($numValid>1)? 'es' : 'is';
	print "<TR><TD width=20></TD><TD class='title3' align=right>$numValid</TD><TD class=\"title3\"> validated MS/MS Analys$anaStrg</TD></TR>\n";
}
my $numValProt=scalar keys %proteinList;
if ($numValProt) {
	my $protStrg=($numValProt>1)? 'Proteins' : 'Protein';
	print "<TR><TD width=20></TD><TD class=\"title3\" align=right>$numValProt</TD><TD class=\"title3\"> validated $protStrg</TD></TR>\n";
	#>Quantifs
	my $numQuantif=scalar keys %quantifList;
	if ($numQuantif) {
		my $quantifStrg=($numQuantif>1)? 'Quantifications' : 'Quantification';
		print "<TR><TD width=20></TD><TD class=\"title3\" align=right>$numQuantif</TD><TD class=\"title3\"> $quantifStrg</TD></TR>\n";
	}
	#>Explor ana
	if ($numExplorAna) {
		my $explorStrg=($numExplorAna>1)? 'Analyses' : 'Analysis';
		print "<TR><TD width=20></TD><TD class=\"title3\" align=right>$numExplorAna</TD><TD class=\"title3\"> Exploratory $explorStrg</TD></TR>\n";
	}
	#>GO
	if ($numGOAna) {
		my $goStrg=($numGOAna>1)? 'Analyses' : 'Analysis';
		print "<TR><TD width=20></TD><TD class=\"title3\" align=right>$numGOAna</TD><TD class=\"title3\"> Gene Ontlogy $goStrg</TD></TR>\n";
	}
	if ($numPathwayAna) {
		my $pathwayStrg=($numPathwayAna>1)? 'Analyses' : 'Analysis';
		print "<TR><TD width=20></TD><TD class=\"title3\" align=right>$numPathwayAna</TD><TD class=\"title3\"> Pathway $pathwayStrg</TD></TR>\n";
	}
	if ($numGseaAna) {
		my $gseaStrg=($numGseaAna>1)? 'Analyses' : 'Analysis';
		print "<TR><TD width=20></TD><TD class=\"title3\" align=right>$numGseaAna</TD><TD class=\"title3\"> Gene Set Enrichment $gseaStrg</TD></TR>\n";
	}
}
print "</TABLE><BR><BR>\n";


##################
####>Deletion<####
##################
if ($action eq 'delete') {

	$dbh->disconnect;
	my $disabStrg=($item ne 'EXPERIMENT' && ($numGOAna || $existDesignQuantif || $numExplorAna || $numPathwayAna || $numGseaAna || $linkedBioSamp))? 'disabled' : '';
	print "<FORM name=\"selAnaForm\" action=\"./deleteProjectItem.cgi\" method=\"post\">\n";
	print "<INPUT type=\"hidden\" name=\"ITEM\" value=\"$item\">\n";
	print "<INPUT type=\"hidden\" name=\"ID\" value=\"$itemID\">\n";

	#foreach my $it ('EXPERIMENT','GEL2D','SPOT','SAMPLE','ANALYSIS') {
	#	next unless defined($itemList{$it});
	#	my $inputName=($it eq 'EXPERIMENT')? 'expList' : ($it eq 'GEL2D')? 'gel2dList' : ($it eq 'SPOT')? 'spotList' : ($it eq 'SAMPLE')? 'sampleList' : 'anaList';
	#	foreach my $itID (@{$itemList{$it}}) {
	#		print "<INPUT type=\"hidden\" name=\"$inputName\" value=\"$itID\">\n";
	#	}
	#}

	print qq
|<TABLE cellpadding=4><TR>
<TD bgcolor=#DD0000><INPUT type="submit" value="Delete" $disabStrg></TD>
<TD>&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="cancelAction()"></TD>
</TR></TABLE>
</FORM>
|;
	if ($disabStrg) {
		print "<BR><BR><FONT class=\"title3\" color=\"#DD0000\">Cannot delete $itemType involved in design-Quantification and/or GO/Exploratory/Pathway/GSEA analyses";
		#print " or linked to Biological Samples" if $item eq 'EXPERIMENT';
		print "!</FONT>\n";
	}
	print qq
|</CENTER>
</BODY>
</HTML>
|;
	exit;
}

#######################
####>Moving Branch<####
#######################
else { # 'move'
	####>Move Form<####
	print "<FORM name=\"moveForm\" method=\"POST\" onsubmit=\"return checkForm(this);\">\n";
	print "<INPUT type=\"hidden\" name=\"ITEM\" value=\"$item\" />\n";
	print "<INPUT type=\"hidden\" name=\"ID\" value=\"$itemID\" />\n";
	print "<INPUT type=\"hidden\" name=\"ACT\" value=\"$action\" />\n";
	my $disabStrg=($item ne 'ANALYSIS' && ($numGOAna || $existDesignQuantif || $numExplorAna || $numPathwayAna || $numGseaAna || ($item ne 'EXPERIMENT' && $linkedBioSamp)))? 'disabled' : '';
	print "<FONT class=\"title2\">Move $itemType to :<SELECT name=\"newBranchID\" style=\"font-size:16px;font-weight:bold;\" $disabStrg>\n";
	print "<OPTION value=\"\">-= Select a destination $parentType =-</OPTION>\n";

	if ($item eq 'EXPERIMENT') {
		my $sthName=$dbh->prepare("SELECT NAME,STATUS FROM PROJECT WHERE ID_PROJECT=?");
		my @userInfo=&promsMod::getUserInfo($dbh,$userID); # $userInfo[2] is ref to access info for all projects
		my %editProjects;
		foreach my $projID (keys %{$userInfo[2]}) {
			if (${$userInfo[2]}{$projID} ne 'guest') {
				next if $projID==$projectID;
				$sthName->execute($projID);
				my ($name,$status)=$sthName->fetchrow_array;
				$status=0 unless $status;
				$editProjects{$projID}=$name if $status <= 0; # project must be still on-going
			}
		}
		$sthName->finish;
		foreach my $projID (sort{lc($editProjects{$a}) cmp lc($editProjects{$b})} keys %editProjects) {
			my $selectString=($newParentBranchID eq "PROJECT:$projID")? ' selected' : '';
			print "<OPTION value=\"PROJECT:$projID\"$selectString>&nbsp;$editProjects{$projID}</OPTION>\n";
		}
	}
	else {
		my @sthPar;
		my $restrict2ExpID=0;
		if ($item eq 'GEL2D' || $item eq 'SAMPLE') {
			$sthPar[0]=$dbh->prepare("SELECT 'EXPERIMENT',ID_EXPERIMENT,NAME FROM EXPERIMENT WHERE ID_PROJECT=$projectID ORDER BY DISPLAY_POS ASC");
		}
		#elsif ($item eq 'SPOT') {
		#	$sthPar[0]=$dbh->prepare("SELECT 'GEL2D',ID_GEL2D,EXPERIMENT.NAME,GEL2D.NAME FROM EXPERIMENT,GEL2D WHERE EXPERIMENT.ID_EXPERIMENT=GEL2D.ID_EXPERIMENT AND ID_PROJECT=$projectID ORDER BY EXPERIMENT.DISPLAY_POS ASC,GEL2D.DISPLAY_POS ASC");
		#}
		elsif ($item eq 'ANALYSIS') {
			if ($existDesignQuantif) {
				my @anaInfo=&promsMod::getItemInfo($dbh,'ANALYSIS',$itemID);
				$restrict2ExpID=$anaInfo[1]{ID}; # parent experiment ID
			}
			$sthPar[0]=$dbh->prepare("SELECT 'SPOT',SPOT.ID_SPOT,CONCAT_WS(' > ',EXPERIMENT.NAME,GEL2D.NAME,SPOT.NAME),EXPERIMENT.ID_EXPERIMENT FROM EXPERIMENT,GEL2D,SPOT,SAMPLE
										WHERE EXPERIMENT.ID_EXPERIMENT=GEL2D.ID_EXPERIMENT AND GEL2D.ID_GEL2D=SPOT.ID_GEL2D
											AND SPOT.ID_SPOT=SAMPLE.ID_SPOT AND ID_PROJECT=$projectID
										ORDER BY EXPERIMENT.DISPLAY_POS ASC,GEL2D.DISPLAY_POS ASC,SPOT.NAME ASC"
									);
			$sthPar[1]=$dbh->prepare("SELECT 'SAMPLE',ID_SAMPLE,CONCAT_WS(' > ',EXPERIMENT.NAME,SAMPLE.NAME),EXPERIMENT.ID_EXPERIMENT FROM EXPERIMENT,SAMPLE WHERE EXPERIMENT.ID_EXPERIMENT=SAMPLE.ID_EXPERIMENT AND ID_SPOT IS NULL AND ID_PROJECT=$projectID ORDER BY EXPERIMENT.DISPLAY_POS ASC,SAMPLE.DISPLAY_POS ASC");
		}
		foreach my $sth (@sthPar) {
			$sth->execute;
			while (my ($parItem,$parID,$hierarchy,$expID)=$sth->fetchrow_array) { # $expID only for analyses
				next if ($parItem eq $itemInfo[-2]{'ITEM'} && $parID==$itemInfo[-2]{'ID'}); # skip current parent
				next if ($item eq 'ANALYSIS' && $existDesignQuantif && $expID != $restrict2ExpID); # restrict to same experiment
				my $selectString=($newParentBranchID eq "$parItem:$parID")? ' selected' : '';
				print "<OPTION value=\"$parItem:$parID\"$selectString>$hierarchy</OPTION>\n";
			}
			$sth->finish;
		}
	}


	#else { # SAMPLE or ANALYSIS
	#	my $sthExp=$dbh->prepare("SELECT ID_EXPERIMENT,NAME FROM EXPERIMENT WHERE ID_PROJECT=$projectID ORDER BY DISPLAY_POS ASC");
	#	if ($item eq 'SAMPLE') {
	#		$sthExp->execute();
	#		while (my($expID,$expName)=$sthExp->fetchrow_array) {
	#			next if $expID==$parentID;
	#			my $selectString=($newParentBranchID==$expID)? 'selected' : '';
	#			print "<OPTION value=\"$expID\" $selectString>&nbsp;$expName</OPTION>\n";
	#		}
	#	}
	#	elsif ($item eq 'ANALYSIS') {
	#		my $sthSamp=$dbh->prepare("SELECT ID_SAMPLE,NAME FROM SAMPLE WHERE ID_EXPERIMENT=? ORDER BY DISPLAY_POS ASC");
	#		$sthExp->execute();
	#		while (my($expID,$expName)=$sthExp->fetchrow_array) {
	#			print "<OPTION value=\"0\">&nbsp;+ $expName:</OPTION>\n";
	#			$sthSamp->execute($expID);
	#			while (my($sampID,$sampName)=$sthSamp->fetchrow_array) {
	#				next if $sampID==$parentID;
	#				my $selectString=($newParentBranchID==$sampID)? 'selected' : '';
	#				print "<OPTION value=\"$sampID\" $selectString>&nbsp&nbsp&nbsp&nbsp;- $sampName</OPTION>\n";
	#			}
	#		}
	#		$sthSamp->finish;
	#	}
	#	$sthExp->finish;
	#}
	print "</SELECT>\n<BR><TABLE cellpadding=4><TR>\n<TD class='title3' colspan=2>";
	if ($item eq 'EXPERIMENT') {print "<INPUT type=checkbox name=\"duplClassif\" $disabStrg/>Duplicate custom Lists and Themes";}
	else {print "&nbsp;";}
	print qq
|</TD></TR><TR>
<TD bgcolor=#DD0000 width=170><INPUT type="submit" name="move" value="Move $itemType" style="width:170px;" $disabStrg/></TD>
<TD>&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="cancelAction()"/></TD>
</TR></TABLE>
</FORM>
|;

	if ($disabStrg) {
		print "<BR><BR><FONT class=\"title3\" color=\"#DD0000\">Cannot move $itemType involved in design-Quantification and/or GO/Exploratory/Pathway/GSEA analyses";
		print " or linked to Biological Samples" if $item eq 'EXPERIMENT';
		print "!</FONT>\n";
	}
	unless ($move) {
		print "</CENTER></BODY>\n</HTML>\n";
		$dbh->disconnect;
		exit;
	}
	else {
		####>Effectively moving Branch<####
		print "</CENTER><BR><BR>\n";
		print "<FONT class=\"title3\">Moving $itemType and all its contents...";
		my ($newParentItem,$newParentID)=split(':',$newParentBranchID);
# exit;
		###>Updating item itself<####
		if ($newParentItem eq 'SPOT') {
			my ($maxDisplayPos,$sampleID)=$dbh->selectrow_array("SELECT MAX(DISPLAY_POS),ID_SAMPLE FROM ANALYSIS WHERE ID_SAMPLE=(SELECT ID_SAMPLE FROM SAMPLE WHERE ID_SPOT=$newParentID) GROUP BY ID_SAMPLE");
			$maxDisplayPos++;
			$dbh->do("UPDATE ANALYSIS SET ID_SAMPLE=$sampleID,DISPLAY_POS=$maxDisplayPos,UPDATE_USER='$userID',UPDATE_DATE='$date' WHERE ID_ANALYSIS=$itemID");
			&updatePosition('SAMPLE',$sampleID); # updating former seeblings
		}
		else {
			my ($maxDisplayPos)=$dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM $item WHERE ID_$newParentItem=$newParentID");
			$maxDisplayPos++;
			$dbh->do("UPDATE $item SET ID_$newParentItem=$newParentID,DISPLAY_POS=$maxDisplayPos,UPDATE_USER='$userID',UPDATE_DATE='$date' WHERE ID_$item=$itemID");
			&updatePosition($parentItem,$parentID); # updating former seeblings
		}

		###>Upadting samples linked to spots<####
		if ($item eq 'GEL2D') {
			my $sthUpSa = $dbh->prepare("UPDATE SAMPLE SET ID_EXPERIMENT=? WHERE ID_SAMPLE=?");
			my $sthSeSa = $dbh->prepare("SELECT ID_SAMPLE FROM SAMPLE, SPOT WHERE SPOT.ID_GEL2D=? AND SAMPLE.ID_SPOT=SPOT.ID_SPOT");
			$sthSeSa->execute($itemID);
			while (my ($sampleID) = $sthSeSa->fetchrow_array) {
				$sthUpSa->execute($newParentID, $sampleID);
			}
			$sthSeSa->finish;
			$sthUpSa->finish;
		}


		###>Updating validated proteins if any and new Project<###
		if ($item eq 'EXPERIMENT' && $numValProt) { # new Project => update proteins

			##>Fetching new Project's proteins<##
			my %newProjProteins;
			my $sthNewProt=$dbh->prepare("SELECT ID_PROTEIN,IDENTIFIER FROM PROTEIN WHERE ID_PROJECT=$newParentID");
			$sthNewProt->execute;
			while (my ($protID,$identifier)=$sthNewProt->fetchrow_array) {
				$newProjProteins{$identifier}=$protID;
			}
			$sthNewProt->finish;
			print '.';

			##>Finding whether proteins are unique to selected Experiment<##
			my (%movedProtAnalyses,%validIdentifiers,%uniqProtAna,%notUniqProteins);
			my $sthAnaProt=$dbh->prepare("SELECT PROTEIN.ID_PROTEIN,IDENTIFIER FROM PROTEIN,ANALYSIS_PROTEIN WHERE PROTEIN.ID_PROTEIN=ANALYSIS_PROTEIN.ID_PROTEIN AND ID_ANALYSIS=?");
			foreach my $anaID (keys %movedValidAnalyses) {
				$sthAnaProt->execute($anaID);
				while (my ($protID,$identifier)=$sthAnaProt->fetchrow_array) {
					$validIdentifiers{$protID}=$identifier;
				}
			}
			$sthAnaProt->finish;
			my $sthProtAna=$dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=?");
			foreach my $protID (keys %validIdentifiers) {
				$sthProtAna->execute($protID);
				my $numAna;
				while (my ($anaID)=$sthProtAna->fetchrow_array) {
					$numAna++;
					push @{$movedProtAnalyses{$protID}},$anaID if $movedValidAnalyses{$anaID};
				}
				if (scalar @{$movedProtAnalyses{$protID}}==$numAna) {
					$uniqProtAna{$protID}=1;
				}
				else {
					$notUniqProteins{$protID}=1;
				}
			}
			$sthProtAna->finish;
			print '.';

			####>Finding classification
			my (%protList,%classList,%catList,%catType,%catLink,%existingClassif);
			if ($duplClassif==1) {
				#my $sthClassifProt=$dbh->prepare("SELECT CATEGORY.ID_CATEGORY,CATEGORY.ID_CLASSIFICATION FROM CATEGORY WHERE CATEGORY.ID_CATEGORY IN (SELECT CATEGORY_PROTEIN.ID_CATEGORY FROM CATEGORY_PROTEIN WHERE CATEGORY_PROTEIN.ID_PROTEIN=?)");
				my $sthClassifProt=$dbh->prepare("SELECT C.ID_CATEGORY,C.ID_CLASSIFICATION,GROUP_CONCAT(CP.ID_CATEGORY_PROTEIN) FROM CATEGORY C,CATEGORY_PROTEIN CP WHERE C.ID_CATEGORY=CP.ID_CATEGORY AND CP.ID_PROTEIN=? GROUP BY C.ID_CATEGORY");
				foreach my $protID (keys %validIdentifiers) {
					$sthClassifProt->execute($protID);
					while (my ($catID,$classID,$listType,$cpIdList)=$sthClassifProt->fetchrow_array) {
						@{$protList{$protID}{$catID}}=split(',',$cpIdList);
						$classList{$classID}=0;
						$catList{$classID}{$catID}=0;
					}
				}
				$sthClassifProt->finish;

				##>Duplicating classification and category
				#Finding existing name
				my $sthClassifName=$dbh->prepare ("SELECT NAME,ID_CLASSIFICATION FROM CLASSIFICATION WHERE ID_PROJECT=$newParentID");
				$sthClassifName->execute;
				while (my ($className, $classID)=$sthClassifName->fetchrow_array) {
					$existingClassif{$className}=$classID;
				}
				$sthClassifName->finish;

				#Prepraring DB access
				my $sthSelectClass=$dbh->prepare("SELECT NAME,DES,COMMENTS FROM CLASSIFICATION WHERE ID_CLASSIFICATION=?");
				my $sthInsertClass=$dbh->prepare("INSERT INTO CLASSIFICATION (ID_CLASSIFICATION,ID_PROJECT,NAME,DES,COMMENTS,UPDATE_DATE,UPDATE_USER) values (?,$newParentID,?,?,?,NOW(),'$userID')");
				my ($newClassID)=$dbh->selectrow_array("SELECT MAX(ID_CLASSIFICATION) FROM CLASSIFICATION");
				my $sthSelectCat=$dbh->prepare("SELECT NAME,DES,COMMENTS,LIST_TYPE,DISPLAY_POS FROM CATEGORY WHERE ID_CATEGORY=?");
				my $sthInsertCat=$dbh->prepare("INSERT INTO CATEGORY (ID_CATEGORY,ID_CLASSIFICATION,NAME,DES,COMMENTS,LIST_TYPE,DISPLAY_POS,UPDATE_DATE,UPDATE_USER) values (?,?,?,?,?,?,NOW(),'$userID')");
				my ($newCatID)=$dbh->selectrow_array("SELECT MAX(ID_CATEGORY) FROM CATEGORY");

				foreach my $classID (keys %classList) {
					$newClassID++;
					$classList{$classID}=$newClassID;
					$sthSelectClass->execute($classID);
					my ($name,$des,$comment)=$sthSelectClass->fetchrow_array;
					if (defined($existingClassif{$name})) {
						if ($name=~/\s#(\d+)$/) {
							my $number=$1;
							my $newNumber=$number+1;
							$name =~ s/(.*\s#)$number/$1$newNumber/;
						}
						else {$name.=" #2";}
					}
					$sthInsertClass->execute($newClassID,$name,$des,$comment);
					foreach my $catID (keys %{$catList{$classID}}) {
						$newCatID++;
						$catList{$classID}{$catID}=$newCatID;
						$catLink{$catID}=$newCatID;
						$sthSelectCat->execute($catID);
						my ($name,$des,$comment,$listType,$display_pos)=$sthSelectCat->fetchrow_array;
						$sthInsertCat->execute($newCatID,$newClassID,$name,$des,$comment,$listType,$display_pos);
						$catType{$catID}=$listType || 'PROT';
					}
				}
				$sthSelectClass->finish;
				$sthInsertClass->finish;
				$sthSelectCat->finish;
				$sthInsertCat->finish;
			}

			##>Peptides and analysis/protein attributions update queries<##
			my $sthUpPepAtt=$dbh->prepare("UPDATE PEPTIDE_PROTEIN_ATTRIB SET ID_PROTEIN=? WHERE ID_PROTEIN=? AND ID_ANALYSIS=?");
			my $sthUpAnaProt=$dbh->prepare("UPDATE ANALYSIS_PROTEIN SET ID_PROTEIN=? WHERE ID_PROTEIN=? AND ID_ANALYSIS=?");

			##>Moving unique proteins and updating peptide and Analysis attributions<##
			my $sthModSite=$dbh->prepare("SELECT ID_MODIFICATION,RESIDUE,POSITION FROM MODIFICATION_SITE WHERE ID_CATEGORY_PROTEIN=?");
			my $sthMovProt=$dbh->prepare("UPDATE PROTEIN SET ID_PROJECT=$newParentID WHERE ID_PROTEIN=?");
			my $sthDelModSite=$dbh->prepare("DELETE MS FROM MODIFICATION_SITE MS INNER JOIN CATEGORY_PROTEIN CP ON MS.ID_CATEGORY_PROTEIN=CP.ID_CATEGORY_PROTEIN WHERE CP.ID_PROTEIN=?");
			my $sthDelCatProt=$dbh->prepare("DELETE FROM CATEGORY_PROTEIN WHERE ID_PROTEIN=?");
			my $sthDelProt=$dbh->prepare("DELETE FROM PROTEIN WHERE ID_PROTEIN=?");
			my $sthInsCatProt=$dbh->prepare("INSERT INTO CATEGORY_PROTEIN (ID_CATEGORY,ID_PROTEIN) VALUES (?,?)");
			my $sthInsModSite=$dbh->prepare("INSERT INTO MODIFICATION_SITE (ID_CATEGORY_PROTEIN,ID_CATEGORY,ID_MODIFICATION,RESIDUE,SITE) VALUES (?,?,?,?,?)");

			foreach my $protID (keys %uniqProtAna) {
				my %modifSites;
				if ($duplClassif==1) {
					#>Store Site info if Modif site list
					foreach my $catID (keys %{$protList{$protID}}) {
						next if $catType{$protID} eq 'PROT';
						foreach my $catProtID (@{$protList{$protID}{$catID}}) {
							$sthModSite->execute($catProtID);
							while (my ($modID,$res,$pos)=$sthModSite->fetchrow_array) {
								push @{$modifSites{$catID}{$catProtID}},[$modID,$res,$pos];
							}
						}
					}
				}
				$sthDelModSite->execute($protID) ; # remove protein from any Site
				$sthDelCatProt->execute($protID) ; # remove protein from any Category
				my $identifier=$validIdentifiers{$protID};
				if (defined($newProjProteins{$identifier})) { # protein already exists in new Project=>keep new protID
					foreach my $anaID (@{$movedProtAnalyses{$protID}}) {
						$sthUpPepAtt->execute($newProjProteins{$identifier},$protID,$anaID); # update PepAtt
						$sthUpAnaProt->execute($newProjProteins{$identifier},$protID,$anaID); # update AnaProt
					}
					$sthDelProt->execute($protID); # delete protein from old Project
					if ($duplClassif==1)  {
						foreach my $oldCatID (keys %{$protList{$protID}}) {
							if ($catType{$oldCatID} eq 'PROT') {
								$sthInsCatProt->execute($catLink{$oldCatID},$newProjProteins{$identifier});
							}
							else { # Modif-prot list
								foreach my $oldCatProtID (keys %{$modifSites{$oldCatID}}) {
									$sthInsCatProt->execute($catLink{$oldCatID},$newProjProteins{$identifier});
									my ($newCatProtID)=$dbh->last_insert_id(undef,undef,'CATEGORY_PROTEIN','ID_CATEGORY_PROTEIN');
									foreach my $refSite (@{$modifSites{$oldCatID}{$oldCatProtID}}) {
										$sthInsModSite->execute($newCatProtID,$catLink{$oldCatID},@{$refSite});
									}
								}
							}
						}
					}
				}
				else { # move protein (no change to PepAtt and AnaProt)
					$sthMovProt->execute($protID);
					if ($duplClassif==1) {
						foreach my $oldCatID (keys %{$protList{$protID}}) {
							if ($catType{$oldCatID} eq 'PROT') {
								$sthInsCatProt->execute($catLink{$oldCatID},$protID);
							}
							else { # Modif-prot list
								foreach my $oldCatProtID (keys %{$modifSites{$oldCatID}}) {
									$sthInsCatProt->execute($catLink{$oldCatID},$protID);
									my ($newCatProtID)=$dbh->last_insert_id(undef,undef,'CATEGORY_PROTEIN','ID_CATEGORY_PROTEIN');
									foreach my $refSite (@{$modifSites{$oldCatID}{$oldCatProtID}}) {
										$sthInsModSite->execute($newCatProtID,$catLink{$oldCatID},@{$refSite});
									}
								}
							}
						}
					}
				}
			}
			$sthMovProt->finish;
			$sthDelProt->finish;
			$sthDelCatProt->finish;
			$sthDelModSite->finish;
			print '.';

			##>Duplicating non-unique proteins and updating Peptide and Analysis attributions<##
			my $sthSelProt=$dbh->prepare("SELECT ID_MASTER_PROTEIN,IDENTIFIER,ALIAS,PROT_DES,ORGANISM,MW,PROT_SEQ,PROT_LENGTH,COMMENTS,UPDATE_USER,UPDATE_DATE FROM PROTEIN WHERE ID_PROTEIN=?");
			my $sthDupProt=$dbh->prepare("INSERT INTO PROTEIN (ID_PROTEIN,ID_PROJECT,ID_MASTER_PROTEIN,IDENTIFIER,ALIAS,PROT_DES,ORGANISM,MW,PROT_SEQ,PROT_LENGTH,COMMENTS,UPDATE_USER,UPDATE_DATE) VALUES (?,$newParentID,?,?,?,?,?,?,?,?,?,?,?)");
			my ($maxProtID)=$dbh->selectrow_array("SELECT MAX(ID_PROTEIN) FROM PROTEIN");

			foreach my $protID (keys %notUniqProteins) {
				my %modifSites;
				if ($duplClassif==1) {
					#>Store Site info if Modif site list
					foreach my $catID (keys %{$protList{$protID}}) {
						next if $catType{$protID} eq 'PROT';
						foreach my $catProtID (@{$protList{$protID}{$catID}}) {
							$sthModSite->execute($catProtID);
							while (my ($modID,$res,$pos)=$sthModSite->fetchrow_array) {
								push @{$modifSites{$catID}{$catProtID}},[$modID,$res,$pos];
							}
						}
					}
				}
				my $identifier=$validIdentifiers{$protID};
				if (defined($newProjProteins{$identifier})) { # protein already exists in new Project=>keep new protID
					foreach my $anaID (@{$movedProtAnalyses{$protID}}) {
						$sthUpPepAtt->execute($newProjProteins{$identifier},$protID,$anaID); # update PepAtt
						$sthUpAnaProt->execute($newProjProteins{$identifier},$protID,$anaID); # update AnaProt
					}
					if ($duplClassif==1) {
						foreach my $oldCatID (keys(%{$protList{$protID}})){
							if ($catType{$oldCatID} eq 'PROT') {
								$sthInsCatProt->execute($catLink{$oldCatID},$newProjProteins{$identifier});
							}
							else { # Modif-prot list
								foreach my $oldCatProtID (keys %{$modifSites{$oldCatID}}) {
									$sthInsCatProt->execute($catLink{$oldCatID},$newProjProteins{$identifier});
									my ($newCatProtID)=$dbh->last_insert_id(undef,undef,'CATEGORY_PROTEIN','ID_CATEGORY_PROTEIN');
									foreach my $refSite (@{$modifSites{$oldCatID}{$oldCatProtID}}) {
										$sthInsModSite->execute($newCatProtID,$catLink{$oldCatID},@{$refSite});
									}
								}
							}
						}
					}
				}
				else { # duplicate protein (create new entry with new protID)
					$sthSelProt->execute($protID);
					my @protData=$sthSelProt->fetchrow_array;
					$sthDupProt->execute(++$maxProtID,@protData);
					foreach my $anaID (@{$movedProtAnalyses{$protID}}) {
						$sthUpPepAtt->execute($maxProtID,$protID,$anaID); # update PepAtt
						$sthUpAnaProt->execute($maxProtID,$protID,$anaID); # update AnaProt
					}
					if ($duplClassif==1) {
						foreach my $oldCatID (keys(%{$protList{$protID}})){
							if ($catType{$oldCatID} eq 'PROT') {
								$sthInsCatProt->execute($catLink{$oldCatID},$maxProtID);
							}
							else { # Modif-prot list
								foreach my $oldCatProtID (keys %{$modifSites{$oldCatID}}) {
									$sthInsCatProt->execute($catLink{$oldCatID},$maxProtID);
									my ($newCatProtID)=$dbh->last_insert_id(undef,undef,'CATEGORY_PROTEIN','ID_CATEGORY_PROTEIN');
									foreach my $refSite (@{$modifSites{$oldCatID}{$oldCatProtID}}) {
										$sthInsModSite->execute($newCatProtID,$catLink{$oldCatID},@{$refSite});
									}
								}
							}
						}
					}
				}
			}
			$sthSelProt->finish;
			$sthDupProt->finish;
			$sthUpPepAtt->finish;
			$sthUpAnaProt->finish;
			$sthModSite->finish;
			$sthInsCatProt->finish;
			$sthInsModSite->finish;
			print '.';

			##>Making new project directory if not exist<##
			mkdir "$promsPath{peptide}/proj_$newParentID" unless -e "$promsPath{peptide}/proj_$newParentID";

			##>Moving peptide .dat files to new project<##
			foreach my $anaID (keys %movedValidAnalyses) {
				dirmove("$promsPath{peptide}/proj_$projectID/ana_$anaID","$promsPath{peptide}/proj_$newParentID/ana_$anaID") if $movedValidAnalyses{$anaID} == 2; # new data file structure
			}
		}

		####>Moving gel files<####
		if ($item eq 'EXPERIMENT') {
			my @gelIDs;
			my $sthGel = $dbh->prepare("SELECT ID_GEL2D FROM GEL2D WHERE ID_EXPERIMENT=?");
			$sthGel->execute($itemID);
			while (my ($gelID) = $sthGel->fetchrow_array) {
				push @gelIDs, $gelID;
			}
			$sthGel->finish;
			mkdir "$promsPath{gel_unix}/project_$newParentID" unless -e "$promsPath{gel_unix}/project_$newParentID";
			foreach my $gelID (@gelIDs){
				move("$promsPath{gel_unix}/project_$projectID/gel2D_$gelID.jpg","$promsPath{gel_unix}/project_$newParentID/gel2D_$gelID.jpg");
			}
		}


		$dbh->commit;

		print " Done.</FONT><BR>\n";

		####>Updating Frames<####
		#exit; #debug
		sleep 2;
		my $selBranchID=lc($item).":$itemID";
		my $URLstring;
		if ($item eq 'EXPERIMENT') {
			$URLstring="parent.location=\"$promsPath{cgi}/openProject.cgi?ID=$newParentID&branchID=$selBranchID&ACT=open\";";
		}
		elsif ($item eq 'ANALYSIS') { # cannot be SPOT or Associated SAMPLE
			my @newParentInfo=&promsMod::getItemInfo($dbh,$newParentItem,$newParentID);
			if ($newParentInfo[2]{'ITEM'} eq 'GEL2D') {
				$URLstring="parent.navFrame.location=\"$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=gel2d:$newParentInfo[2]{'ID'}&itemBranchID=$selBranchID&ACT=nav\";";
			}
			else {
				$URLstring="parent.navFrame.location=\"$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=$selBranchID&ACT=nav\";";
			}
		}
		else {
			$URLstring="parent.navFrame.location=\"$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=$selBranchID&ACT=nav\";";
		}
		$dbh->disconnect;
		print qq
|<SCRIPT type="text/javascript">
top.promsFrame.selectedAction='summary';
$URLstring
</SCRIPT>
</BODY>
</HTML>
|;
		exit;
	}
}


###############################
####>Item List subroutines<####
###############################
sub getGel2dList {
	my ($expID)=@_;
	$sthG2D->execute($expID);
	while (my ($gelID)=$sthG2D->fetchrow_array) {
		push @{$itemList{'GEL2D'}},$gelID;
		&getSpotList($gelID);
	}
	$sthG2D->finish;
}
sub getSpotList {
	my ($gelID)=@_;
	$sthSpot->execute($gelID);
	while (my ($spotID)=$sthSpot->fetchrow_array) {
		push @{$itemList{'SPOT'}},$spotID;
		&getSampList($sthSamp2,$spotID);
	}
	$sthSpot->finish;
}
sub getSampList {
	my ($sthSamp,$parentID)=@_;
	$sthSamp->execute($parentID);
	while (my ($sampID)=$sthSamp->fetchrow_array) {
		push @{$itemList{'SAMPLE'}},$sampID;
		&getAnaList($sampID);
	}
	$sthSamp->finish;
}
sub getAnaList {
	my ($sampID)=@_;
	$sthAna->execute($sampID);
	while (my ($anaID,$valStatus)=$sthAna->fetchrow_array) {
		push @{$itemList{'ANALYSIS'}},$anaID;
		if ($valStatus>=1) {
			&getProtList($anaID);
			$movedValidAnalyses{$anaID}=$valStatus;
		}
		if ($valStatus==-1) {$numNotImp++;} elsif ($valStatus==0) {$numUnValid++;} elsif ($valStatus==1) {$numPartValid++;} else {$numValid++;}
		unless ($linkedBioSamp) {
			$sthBioSampObs->execute($anaID);
			($linkedBioSamp)=$sthBioSampObs->fetchrow_array;
		}
	}
}
sub getProtList {
	my ($anaID)=@_;
	$sthProt->execute($anaID);
	while (my ($protID)=$sthProt->fetchrow_array) {
		$proteinList{$protID}=1;
	}
	# Quantif
	$sthQuantif->execute($anaID);
	while (my ($quantifID,$designID)=$sthQuantif->fetchrow_array) {
		$quantifList{$quantifID}=1;
		$existDesignQuantif=1 if $designID;
	}
}

#####################################################
####<Update display position of former seeblings>####
#####################################################
sub updatePosition {
	my ($parItem,$parID)=@_;
	my $sthSel=$dbh->prepare("SELECT ID_$item FROM $item WHERE ID_$parItem=$parID ORDER BY DISPLAY_POS ASC");
	my $sthUp=$dbh->prepare("UPDATE $item SET DISPLAY_POS=? WHERE ID_$item=?");
	$sthSel->execute;
	my $displayPos=0;
	while (my ($itID)=$sthSel->fetchrow_array) {
		$displayPos++;
		$sthUp->execute($displayPos,$itID);
	}
	$sthSel->finish;
	$sthUp->finish;
}

#####################
####<End Project>####
#####################
sub endOrArchiveProject {

	my $titleString=($action eq 'end')? 'Ending Project' : 'Archiving Project';
	print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>End or Archive Project</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT type="text/javascript">
// AJAX --->
function processAnalyses(anaList) {
	//Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if (!XHR) {return false;}
	XHR.open("GET","./send2Biologist.cgi?call=ajax&endVal=1&anaList="+anaList,true);
    XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText && !XHR.responseText.match('__OK__')) {
			sendProcessError();
		}
	}
	XHR.send(null);
}
function sendProcessError() {
	//Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if (!XHR) {return false;}
	XHR.open("GET","./editBranch.cgi?ACT=error",true);
    XHR.onreadystatechange=function() {
		if (XHR.readyState==4) {
			//nothing to do
		}
	}
	XHR.send(null);
}
var XHR=null;
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
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<IMG src="$promsPath{images}/engrenage.gif">
<BR>
<FONT class="title">$titleString</FONT>
</CENTER>
<BR><BR>
|;
	my ($status)=$dbh->selectrow_array("SELECT STATUS FROM PROJECT WHERE ID_PROJECT=$itemID");
	$status=0 unless $status;

	####<Ending all on-going analyses (AJAX call to send2Biologist.cgi)>####
	if ($status <= 0) { # On-going project

		###<Fetching list of not fully validated analyses>###
		my $sthAna=$dbh->prepare("SELECT A.ID_ANALYSIS FROM EXPERIMENT E,SAMPLE S,ANALYSIS A WHERE E.ID_PROJECT=$itemID AND E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SAMPLE=A.ID_SAMPLE AND A.VALID_STATUS < 2");
		$sthAna->execute;
		my @analysisList;
		while (my ($anaID)=$sthAna->fetchrow_array) {push @analysisList,$anaID;}
		$sthAna->finish;

		if (scalar @analysisList) {
			$dbh->disconnect; # safer
			unlink glob "$promsPath{tmp}/ajaxAnaReport_$userID.*" if glob "$promsPath{tmp}/ajaxAnaReport_$userID.*"; # just in case
			print "<FONT class=\"title2\">Ending Analysis validation...";

			for (my $i=0; $i<=$#analysisList; $i+=5) {
				my $lastIdx=$i+4; $lastIdx=$#analysisList if $lastIdx>$#analysisList;
				my $subAnaList=join(':',@analysisList[$i..$lastIdx]);
				print qq |<SCRIPT type="text/javascript">processAnalyses('$subAnaList');</SCRIPT>|;
				my $processEnd="$promsPath{tmp}/ajaxAnaReport_$userID.end";
				my $processError="$promsPath{tmp}/ajaxAnaReport_$userID.error";
				my $procResult='';
				my $count=0;
				while ($procResult eq '') {
					sleep 1;
					$count++;
					$procResult=(-e $processEnd)? 'OK' : (-e $processError || $count > 500)? 'ERROR' : '';
					print "<!--*-->\n" unless $count % 30; # keeps connection alive ~ every 30 sec
				}

				####<Processing response>####
				if ($procResult eq 'OK') {
					unlink $processEnd;
					print '.';
				}
				else {
					unlink $processError if -e $processError;
					print qq
|</FONT><BR><BR>
<FONT class="title3">ERROR: Unexpected response from myProMS Server<BR>
Operation aborted.</FONT>
&nbsp;&nbsp;<INPUT type="button" value=" OK " onclick="parent.optionFrame.selectOption(parent.optionFrame.document.getElementById('summary'))"/>
</BODY>
</HTML>
|;
					exit;
				}
			}

			print " Done.</FONT><BR>\n";

			$dbh=&promsConfig::dbConnect; # reconnect
		}

		###<Cleaning Upload directory (mzXML files)>###
		#remove_tree("$promsPath{tmp}/upload/project_$itemID") if -e "$promsPath{tmp}/upload/project_$itemID";
		rmtree("$promsPath{tmp}/upload/project_$itemID") if -e "$promsPath{tmp}/upload/project_$itemID";
	}

	####<Compressing all project's data files>####
	if ($action eq 'archive') {
		print "<FONT class=\"title2\">Archiving data files...</FONT><BR>\n";

		###<Peptide data>###
		print "<FONT class=\"title3\">&nbsp;&nbsp;-Peptides...";
		&archiveDirectory("$promsPath{peptide}","proj_$itemID");
		print " Done.</FONT><BR>\n";

		###<2D Gels data>###
		print "<FONT class=\"title3\">&nbsp;&nbsp;-Gels...";
		&archiveDirectory("$promsPath{gel_unix}","project_$itemID");
		print " Done.</FONT><BR>\n";

		###<Quantification data>###
		print "<FONT class=\"title3\">&nbsp;&nbsp;-Quantifications...";
		&archiveDirectory("$promsPath{quantification}","project_$itemID");
		print " Done.</FONT><BR>\n";

		###<GO data>###
		print "<FONT class=\"title3\">&nbsp;&nbsp;-Gene Ontology analyses...";
		&archiveDirectory("$promsPath{go_unix}","project_$itemID");
		print " Done.</FONT><BR>\n";

		###<Exploratory Analyses>###
		print "<FONT class=\"title3\">&nbsp;&nbsp;-Exploratory analyses...";
		&archiveDirectory("$promsPath{explorAna}","project_$itemID");
		print " Done.</FONT><BR>\n";

		###<Pathway Analyses>###
		print "<FONT class=\"title3\">&nbsp;&nbsp;-Pathway analyses...";
		&archiveDirectory("$promsPath{pathAna}","project_$itemID");
		print " Done.</FONT><BR>\n";

		###<GSEA Analyses>###
		print "<FONT class=\"title3\">&nbsp;&nbsp;-GSEA analyses...";
		&archiveDirectory("$promsPath{gsea}","project_$itemID");
		print " Done.</FONT><BR>\n";

	}

	####<Updating Project Status>####
	my $date=strftime("%Y-%m-%d %H:%M:%S",localtime);
	my $newStatus=($action eq 'end')? 1 : 2;
	$dbh->do("UPDATE PROJECT SET STATUS=$newStatus,UPDATE_USER='$userID',UPDATE_DATE='$date' WHERE ID_PROJECT=$itemID");

	$dbh->commit;
	$dbh->disconnect;

	if ($action eq 'end') {
		print qq
|<BR><BR><FONT class="title2">Project was successfully ended.</FONT>
&nbsp;&nbsp;<INPUT type="button" value=" OK " onclick="top.promsFrame.location='$promsPath{cgi}/openProject.cgi?&ACT=open&ID=$itemID&branchID=project:$itemID'"/>
</BODY>
</HTML>
|;
	}
	else {
		print qq
|<BR><BR><FONT class="title2">Project was successfully archived.</FONT>
&nbsp;&nbsp;<INPUT type="button" value=" OK " onclick="top.promsFrame.location='$promsPath{cgi}/selectProject.cgi'"/>
</BODY>
</HTML>
|;
	}
	exit;
}

##############################################
####<Continue an Ended or Archived Project>####
##############################################
sub continueProject {
	my ($status)=$dbh->selectrow_array("SELECT STATUS FROM PROJECT WHERE ID_PROJECT=$itemID");
	$status=0 unless $status;
	my $title=($status==2)? 'Restoring Project' : 'Continuing Project';
	print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>$title</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<IMG src="$promsPath{images}/engrenage.gif">
<BR>
<FONT class="title">$title</FONT>
</CENTER>
<BR><BR>
|;

	####<Un-tarring all project's data files>####
	if ($status==2) {
		print "<FONT class=\"title2\">Un-archiving data files...</FONT><BR>\n";

		###<Peptide data>###
		print "<FONT class=\"title3\">&nbsp;&nbsp;-Peptides...";
		&unArchiveFile("$promsPath{peptide}","proj_$itemID.tar");
		print " Done.</FONT><BR>\n";

		###<2D Gels data>###
		print "<FONT class=\"title3\">&nbsp;&nbsp;-Quantifications...";
		&unArchiveFile("$promsPath{gel_unix}","project_$itemID.tar");
		print " Done.</FONT><BR>\n";

		###<Quantification data>###
		print "<FONT class=\"title3\">&nbsp;&nbsp;-Quantifications...";
		&unArchiveFile("$promsPath{quantification}","project_$itemID.tar");
		print " Done.</FONT><BR>\n";

		###<GO data>###
		print "<FONT class=\"title3\">&nbsp;&nbsp;-Gene Ontology analyses...";
		&unArchiveFile("$promsPath{go_unix}","project_$itemID.tar");
		print " Done.</FONT><BR>\n";

		###<Exploratory Analyses>###
		print "<FONT class=\"title3\">&nbsp;&nbsp;-Exploratory analyses...";
		&unArchiveFile("$promsPath{explorAna}","project_$itemID.tar");
		print " Done.</FONT><BR>\n";

		###<Pathway Analyses>###
		print "<FONT class=\"title3\">&nbsp;&nbsp;-Pathway analyses...";
		&unArchiveFile("$promsPath{pathAna}","project_$itemID.tar");
		print " Done.</FONT><BR>\n";

		###<GSEA Analyses>###
		print "<FONT class=\"title3\">&nbsp;&nbsp;-GSEA analyses...";
		&unArchiveFile("$promsPath{gsea}","project_$itemID.tar");
		print " Done.</FONT><BR>\n";
	}

	####<Updating Project Status>####
	my $date=strftime("%Y-%m-%d %H:%M:%S",localtime);
	my ($newStatus,$actionStrg)=($status==2)? (1,'un-archived') : (0,'re-opened'); # $status is 2 or 1 only
	$dbh->do("UPDATE PROJECT SET STATUS=$newStatus,UPDATE_USER='$userID',UPDATE_DATE='$date' WHERE ID_PROJECT=$itemID");

	$dbh->commit;
	$dbh->disconnect;

	print qq
|<BR><BR><FONT class="title2">Project was successfully $actionStrg.</FONT>
&nbsp;&nbsp;<INPUT type="button" value=" OK " onclick="top.promsFrame.location='$promsPath{cgi}/openProject.cgi?&ACT=open&ID=$itemID&branchID=project:$itemID'"/>
</BODY>
</HTML>
|;
	exit;
}


#####################################################################
####<(Un)Set project sensitivity to auto-end analysis validation>####
#####################################################################
sub setAutoEndValidation {
	my $newStatus=($action eq 'noAutoEndVal')? -1 : 0;
	$dbh->do("UPDATE PROJECT SET STATUS=$newStatus WHERE ID_PROJECT=$itemID");
	$dbh->commit;
	$dbh->disconnect;
	print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT type="text/javascript">
	top.promsFrame.selectedAction='summary';
	top.promsFrame.optionFrame.selectOption();
</SCRIPT>
</HEAD>
</HTML>
|;
	exit;
}


sub archiveDirectory { # with all subdirectory
	my ($parentDir,$selectedDir)=@_;
	return unless -e "$parentDir/$selectedDir";
	chdir $parentDir;
	my @files;
	finddepth(sub {
		return if ($_ eq '.' || $_ eq '..');
		push @files, $File::Find::name;
	},$selectedDir);
	my $tar=Archive::Tar->new();
	for (my $i=0; $i<=$#files; $i+=10) {
		my $lastIdx=$i+9; $lastIdx=$#files if $lastIdx>$#files;
		$tar->add_files(@files[$i..$lastIdx]) or die "<B>**Error**: Archiving failed:</B> $!<BR>\n";
		print '.';
	}
	$tar->write("$selectedDir.tar",COMPRESS_GZIP) or die "<B>**Error**: Archiving failed:</B> $!<BR>\n";

	#remove_tree($selectedDir); # delete original dir
	rmtree($selectedDir); # delete original dir
}

sub unArchiveFile {
	my ($parentDir,$archive)=@_;
	return unless -e "$parentDir/$archive";
	chdir $parentDir;
	my $tar=Archive::Tar->new;
	$tar->read($archive) or die "<B>**Error**: Un-archiving failed:</B> $!<BR>\n";
    $tar->extract() or die "<B>**Error**: Un-archiving failed:</B> $!<BR>\n";

	unlink("$parentDir/$archive");
}

sub writeError {
	print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);
	print "<HTML><BODY>\n";
	open (P_ERR,">$promsPath{tmp}/ajaxAnaReport_$userID.error");
	print P_ERR "__ERROR__\n";
	close P_ERR;
	print "</BODY></HTML>\n";
	exit;
}

####>Revision history<####
# 2.4.2 Add consideration of GSEA analyses and adapt query on PATHWAY_ANALYSIS table (VL 18/11/20)
# 2.4.1 Handles project status=-1 [no auto-end validation] & &setAutoEndValidation (PP 07/06/18)
# 2.4.0 Only uses ITEM & ID parameters for deletion by deleteProjectItem.cgi (PP 02/03/18)
# 2.3.3 Compatible with MODIFICATION_SITE table used for list of modification sites (PP 28/07/17)
# 2.3.2 Prevents deletion of 2D gel/Sample/Analysis if involved in design-quantification (PP 08/11/16)
# 2.3.1 Pathway analysis deletion & archiving<BR>TODO: Make "Move Item" compatible with explor/funct analyses & bioSamples (PP 03/12/14)
# 2.3.0 Added Exploratory analysis folder to project archiving (PP 21/08/14)
# 2.2.9 Uses rmtree instead of remove_tree (PP 10/03/14)
# 2.2.8 Uses Ajax instead of LWP when calling send2Biologist.cgi (PP 08/03/14)
# 2.2.7 system command removal (PP 08/11/13)
# 2.2.6 New attempt to improve URL definition for LWP call (PP 23/10/13)
# 2.2.5 Improve URL definition for LWP call of send2Biologist.cgi (PP 21/10/13)
# 2.2.4 Moving 2D gels while moving EXPERIMENT item (FY 13/09/13)
# 2.2.3 Archiving 2D gels (FY 03/09/13)
# 2.2.2 Move restriction & Project status management: End/Archive/Restore/Continue & new valid file path (PP 04/03/13)
# 2.2.1 full data file management for validated analyses (PP 07/06/11)
