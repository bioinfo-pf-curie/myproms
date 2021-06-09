#!/usr/local/bin/perl -w

################################################################################
# deleteProjectItem.cgi                  2.6.9                                 #
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
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use strict ;
use File::Path qw(rmtree); # remove_tree
use File::Copy qw(move); # needed in &promsMod::removeFromMultiAna
use promsConfig;
use promsMod;
use goAnalysis;
use promsQuantif;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;

#print header(-'content-encoding'=>'no'); warningsToBrowser(1); # DEBUG
##################
####>Argument<####
##################
my ($item,$itemID);
my (@anaList,@sampleList,@spotList,@gel2dList,@expList);
my $projectID;

if (param('ITEM')) { # delete a single EMPTY! Item (except Spot w/ asso Sample)
	$item=lc(param('ITEM'));
	$itemID=param('ID') ;
	if ($item eq 'analysis') {push @anaList,$itemID;}
	elsif ($item eq 'sample') {push @sampleList,$itemID;}
	elsif ($item eq 'spot') {push @spotList,$itemID;}
	elsif ($item eq 'gel2d') {push @gel2dList,$itemID;}
	elsif ($item eq 'experiment') {push @expList,$itemID;}
	elsif ($item eq 'project') {$projectID=$itemID;}
}
else { # Delete an item AND its contents
	if (param('anaList')) {
		@anaList=param('anaList');
		$item='analysis';
		$itemID=$anaList[0];
	}
	if (param('sampleList')) {
		@sampleList=param('sampleList');
		$item='sample';
		$itemID=$sampleList[0];
	}
	if (param('spotList')) {
		@spotList=param('spotList');
		$item='spot';
		$itemID=$spotList[0];
	}
	if (param('gel2dList')) {
		@gel2dList=param('gel2dList');
		$item='gel2d';
		$itemID=$gel2dList[0];
	}
	if (param('expList')) {
		@expList=param('expList');
		$item='experiment';
		$itemID=$expList[0];
	}
}

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;
my ($parentItem,$parentID);
my @itemInfo=&promsMod::getItemInfo($dbh,$item,$itemID); # only for 1st item if multiple

#######################
####>Starting HTML<####
#######################
print header(-'content-encoding'=>'no'); warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Deleting Items</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<BR><BR>
<FONT class="title2">Deleting...|;

if ($item ne 'project')  {
	my $parentIndex=($item eq 'analysis' && $itemInfo[2]{'ITEM'} eq 'GEL2D')? -3 : -2;
	$parentItem=$itemInfo[$parentIndex]{'ITEM'};
	$parentID=$itemInfo[$parentIndex]{'ID'};
	$projectID=$itemInfo[0]{'ID'};
}

####>Analyses<####
if (scalar @anaList) {
	&deleteAnalyses($dbh,$projectID,\@anaList);
}

####>Samples<####
if (scalar @sampleList) {
	&deleteSamples($dbh,$projectID,\@sampleList);
}

####>Spots<####
if (scalar @spotList) {
	&deleteSpots($dbh,$projectID,\@spotList);
}

####>2D-gels<####
if (scalar @gel2dList) {
	&delete2DGels($dbh,$projectID,\@gel2dList);
}

####>Experiments<####
if (scalar @expList) {
	&deleteExperiments($dbh,$projectID,\@expList);
}

####>Projects<####
if ($item eq 'project') { # only 1 project
	my $sthPE=$dbh->prepare("SELECT ID_EXPERIMENT FROM EXPERIMENT WHERE ID_PROJECT=$projectID");
	my @expList;
	$sthPE->execute;
	while (my ($expID)=$sthPE->fetchrow_array) {
		push @expList,$expID;
	}
	$sthPE->finish;
	&deleteExperiments($dbh,$projectID,\@expList);

	###>Deleting Comparisons (before CATEGORY!)<###
	$dbh->do("DELETE FROM CAT_COMPARISON WHERE ID_COMPARISON IN (SELECT ID_COMPARISON FROM COMPARISON WHERE ID_PROJECT=$projectID)") || die $dbh->errstr;
	$dbh->do("DELETE FROM COMPARISON WHERE ID_PROJECT=$projectID") || die $dbh->errstr; # ANA_COMPARISON deleted at analysis level

	###>Deleting Classifications<###
	$dbh->do("DELETE FROM CATEGORY WHERE ID_CLASSIFICATION IN (SELECT ID_CLASSIFICATION FROM CLASSIFICATION WHERE ID_PROJECT=$projectID)") || die $dbh->errstr;
	$dbh->do("DELETE FROM CLASSIFICATION WHERE ID_PROJECT=$projectID") || die $dbh->errstr;

	###>Deleting relevant modifications<###
	$dbh->do("DELETE FROM PROJECT_MODIFICATION WHERE ID_PROJECT=$projectID") || die $dbh->errstr;

	###>Deleting project access<###
	$dbh->do("DELETE FROM PROJECT_ACCESS WHERE ID_PROJECT=$projectID") || die $dbh->errstr;

	###>Deleting BioSamples<###
	my $sthPBS=$dbh->prepare("SELECT ID_BIOSAMPLE FROM PROJECT_BIOSAMPLE WHERE ID_PROJECT=$projectID");
	my $sthUsedBS=$dbh->prepare("SELECT 1 FROM PROJECT_BIOSAMPLE WHERE ID_BIOSAMPLE=? AND ID_PROJECT != $projectID LIMIT 0,1");
	my @delBioSamples;
	$sthPBS->execute;
	while (my($bioSampleID)=$sthPBS->fetchrow_array) {
		$sthUsedBS->execute($bioSampleID);
		my ($usedBioSamp)=$sthUsedBS->fetchrow_array;
		push @delBioSamples,$bioSampleID unless $usedBioSamp;
	}
	$sthPBS->finish;
	$sthUsedBS->finish;

	$dbh->do("DELETE FROM PROJECT_BIOSAMPLE WHERE ID_PROJECT=$projectID");
	my $sthDelBSProp=$dbh->prepare("DELETE FROM BIOSAMPLE_PROPERTY WHERE ID_BIOSAMPLE=?");
	my $sthDelBS=$dbh->prepare("DELETE FROM BIOSAMPLE WHERE ID_BIOSAMPLE=?");
	foreach my $bioSampleID (@delBioSamples) {
		$sthDelBSProp->execute($bioSampleID);
		$sthDelBS->execute($bioSampleID); # Linked Observations have been deleted at analysis level
	}
	$sthDelBSProp->finish;
	$sthDelBS->finish;

	###>Deleting Project-specific BioSample Properties<###
	my $sthPProp=$dbh->prepare("SELECT ID_PROPERTY FROM PROJECT_PROPERTY WHERE ID_PROJECT=$projectID");
	my $sthUsedProp=$dbh->prepare("SELECT 1 FROM PROJECT_PROPERTY WHERE ID_PROPERTY=? AND ID_PROJECT != $projectID LIMIT 0,1");
	my @delProperties;
	$sthPProp->execute;
	while (my($propID)=$sthPProp->fetchrow_array) {
		$sthUsedProp->execute($propID);
		my ($usedProp)=$sthUsedProp->fetchrow_array;
		push @delProperties,$propID unless $usedProp;
	}
	$sthPProp->finish;
	$sthUsedProp->finish;

	$dbh->do("DELETE FROM PROJECT_PROPERTY WHERE ID_PROJECT=$projectID");
	my $sthDelProp=$dbh->prepare("DELETE FROM PROPERTY WHERE ID_PROPERTY=?");
	foreach my $propID (@delProperties) {
		$sthDelProp->execute($propID);
	}
	$sthDelProp->finish;
	
	###> Deleting related METADATA
	deleteMetadata($dbh, "PROJECT", $projectID);

	###>Deleting entry in PROJECT table<###
	$dbh->do("DELETE FROM PROJECT WHERE ID_PROJECT=$projectID") || die $dbh->errstr;
	
	###>Deleting project data on disk<###
	my @pathList=("$promsPath{explorAna}/project_$projectID",
				  "$promsPath{gel_unix}/project_$projectID",
				  "$promsPath{go_unix}/project_$projectID",
				  "$promsPath{peptide}/proj_$projectID",
				  "$promsPath{quantification}/project_$projectID",
				  "$promsPath{tmp}/upload/project_$projectID",
				 );
	foreach my $path (@pathList) {
		#remove_tree($path) if -e $path;
		rmtree($path) if -e $path;
	}
	print '.';
	#$dbh->commit;
}

$dbh->commit;
$dbh->disconnect;

print " Done.</FONT>\n";
sleep 2;
#exit ; #debug

##########################
####>Reloading frames<####
##########################
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
if ($item eq 'project') {
	print "parent.location=\"$promsPath{cgi}/selectProject.cgi\";\n";
}
else {
	print "top.promsFrame.selectedAction='summary';\n";
	my $selBranchID=lc($parentItem).":$parentID";
	#my $selBranchID=($itemInfo[-2]{'ITEM'} eq 'SAMPLE' && $itemInfo[2]{'ITEM'} eq 'GEL2D')? "spot:$itemInfo[3]{ID}" : "sample:$itemInfo[-2]{ID}";

	if ($parentItem eq 'SPOT') {
		print "parent.itemFrame.location=\"$promsPath{cgi}/openProject.cgi?ID=$projectID&GEL=gel2d:$itemInfo[2]{ID}&branchID=$selBranchID&ACT=gel\";\n";
	}
	elsif ($parentItem eq 'GEL2D') { # (re)select parent Gel
		print "parent.navFrame.selectItem(parent.navFrame.getItemIndex('$selBranchID'));\n";
	}
	#>Refresh navFrame
	else {
		print "parent.navFrame.location=\"$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=$selBranchID&ACT=nav&VIEW=\"+parent.navFrame.view;\n";
	}

	#if ($item ne 'gel2d' && $itemInfo[2] && $itemInfo[2]{'ITEM'} eq 'GEL2D') { # itemFrame
	#	print "parent.itemFrame.location=\"./openProject.cgi?ID=$projectID&GEL=gel2d:$itemInfo[2]{ID}&branchID=$selBranchID&ACT=gel\";\n";
	#}
	#elsif ($existGels || $item ne 'gel2d') { #navFrame
	#	print "parent.navFrame.location=\"./openProject.cgi?ID=$projectID&branchID=$selBranchID&ACT=nav\";\n";
	#}
	#else { # remove itemFrame
	#	print "top.promsFrame.location=\"./openProject.cgi?ID=$projectID&branchID=$selBranchID\";\n";
	#}
}
print qq
|</SCRIPT>
</BODY>
</HTML>
|;
exit;

######################################################################################

###########################
####<<< SUBROUTINES >>>####
###########################
sub deleteAnalyses {
	my ($dbh,$projectID,$refAnaList)=@_;
	my %modSamples;
	my $sthDelAD=$dbh->prepare("DELETE FROM ANALYSIS_DATABANK WHERE ID_ANALYSIS=?");
	my $sthDelAM=$dbh->prepare("DELETE FROM ANALYSIS_MODIFICATION WHERE ID_ANALYSIS=?");
	my $sthDelASL=$dbh->prepare("DELETE FROM ANALYSIS_SWATH_LIB WHERE ID_ANALYSIS=?");
	my $sthSelASL=$dbh->prepare("SELECT ID_SWATH_LIB FROM ANALYSIS_SWATH_LIB WHERE ID_ANALYSIS=?");
	my $sthS=$dbh->prepare("SELECT ID_SAMPLE FROM ANALYSIS WHERE ID_ANALYSIS=?");
	foreach my $anaID (@{$refAnaList}) {
		$sthDelAD->execute($anaID);
		$sthDelAM->execute($anaID);
		$sthS->execute($anaID);
		$sthSelASL->execute($anaID);
		my $swath=$sthSelASL->fetchrow_array;
		$sthDelASL->execute($anaID) if $swath;
		my ($sampID)=$sthS->fetchrow_array;
		$modSamples{$sampID}=1;
		my ($validStatus,$fileFormat,$dataFile)=$dbh->selectrow_array("SELECT VALID_STATUS,FILE_FORMAT,DATA_FILE FROM ANALYSIS WHERE ID_ANALYSIS=$anaID;");
		if ($validStatus<=1) {
			###<Cleanup RANK_PROTEIN_MATCH
			$dbh->do("DELETE FROM RANK_PROTEIN_MATCH WHERE ID_ANALYSIS=$anaID") || die $dbh->errstr;
			###<Cleanup QUERY_MODIFICATION
			$dbh->do("DELETE QUERY_MODIFICATION FROM QUERY_MODIFICATION INNER JOIN QUERY_VALIDATION ON QUERY_MODIFICATION.ID_QUERY=QUERY_VALIDATION.ID_QUERY WHERE QUERY_VALIDATION.ID_ANALYSIS=$anaID") || die $dbh->errstr;
			###<Cleanup QUERY_VALIDATION
			$dbh->do("DELETE FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$anaID") || die $dbh->errstr;
			###<Cleanup PROTEIN_VALIDATION
			$dbh->do("DELETE FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$anaID") || die $dbh->errstr;

			###<Deleting temporary data
			#remove_tree("$promsPath{valid}/ana_$anaID");
			rmtree("$promsPath{valid}/ana_$anaID");
			if ($fileFormat =~ /\.PDM/) { # Proteome Discoverer
				(my $anaFile=$dataFile) =~ s/\_\d+\.*\d*\.pdm/\.ana/g;
				&promsMod::removeFromMultiAna($anaID,$projectID,$anaFile);
			}
		}
		&deleteValidData($dbh,$projectID,$anaID,$fileFormat,$dataFile) if $validStatus>=1; # delete validated data

		#### Cleanup validation history of current analysis
		$dbh->do("DELETE FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$anaID") || die $dbh->errstr;

		####>Cleanup ITEM table
		$dbh->do("DELETE FROM ANALYSIS WHERE ID_ANALYSIS=$anaID") || die $dbh->errstr;
		print '.';
		$dbh->commit;
	}
	$sthDelAD->finish;
	$sthDelAM->finish;
	$sthS->finish;
	
	####> Update PROJECT_MODIFICATION data
	#&updateProjectModification($dbh,$projectID);
	foreach my $sampID (keys %modSamples) {
		&updateBrothersPosition('SAMPLE',$sampID,'ANALYSIS');
	}
}

#sub updateProjectModification {
#	my ($dbh,$projectID)=@_;
#	my %modInAM=();
#	my $sthGetAMinP=$dbh->prepare("SELECT DISTINCT(ID_MODIFICATION) FROM EXPERIMENT E, SAMPLE S, ANALYSIS A , ANALYSIS_MODIFICATION AM WHERE E.ID_EXPERIMENT=S.ID_EXPERIMENT AND S.ID_SAMPLE=A.ID_SAMPLE AND AM.ID_ANALYSIS=A.ID_ANALYSIS AND E.ID_PROJECT=$projectID");
#	$sthGetAMinP->execute;
#	while (my ($modificationID) = $sthGetAMinP->fetchrow_array  ) {
#		$modInAM{$modificationID}=1;
#	}
#	$sthGetAMinP->finish;
#	my $sthGetPM=$dbh->prepare("SELECT ID_MODIFICATION FROM PROJECT_MODIFICATION WHERE ID_PROJECT=$projectID");
#	my $sthDelPM=$dbh->prepare("DELETE FROM PROJECT_MODIFICATION WHERE ID_PROJECT=$projectID AND ID_MODIFICATION=?");
#	$sthGetPM->execute;
#	while (my ($modificationID) = $sthGetPM->fetchrow_array  ) {
#		if (!defined($modInAM{$modificationID})) {
#			$sthDelPM->execute($modificationID);
#		}
#	}
#	$sthGetPM->finish;
#	$sthDelPM->finish;
#	$dbh->commit;
#}

sub deleteMetadata {
	my ($dbh,$item,$itemID)=@_;
	my $sthSelME=$dbh->prepare("SELECT ID_META_ANNOTATION FROM META_ANNOTATION WHERE ID_$item=?");
	my $sthDelAI=$dbh->prepare("DELETE FROM ANNOTATION_ITEM WHERE ID_META_ANNOTATION=?");
	my $sthDelME=$dbh->prepare("DELETE FROM META_ANNOTATION WHERE ID_$item=?");
	
	my $metaPath = "$promsPath{metadata}/proj_$projectID/";
	if ($item eq 'EXPERIMENT') {
		$metaPath .= "exp_$itemID/";
	} elsif($item eq 'SAMPLE') {
		my $sthSE=$dbh->prepare("SELECT ID_EXPERIMENT FROM SAMPLE WHERE ID_SAMPLE=?");
		$sthSE->execute($itemID);
		my ($expID)=$sthSE->fetchrow_array;
		$sthSE->finish;
		
		$metaPath .= "exp_$expID/samp_$itemID/";
	}
	
	my @metaIDs;
	rmtree($metaPath) if(-e $metaPath);
	$sthSelME->execute($itemID);
	while (my ($metaID)=$sthSelME->fetchrow_array) {
		push(@metaIDs, $metaID);
	}
	
	if(@metaIDs) {
		foreach my $metaID (@metaIDs) {
			$sthDelAI->execute($metaID);
		}
		$sthDelME->execute($itemID);
	}
	
	$sthSelME->finish;
	$sthDelME->finish;
	$sthDelAI->finish;
}

sub deleteSamples {
	my ($dbh,$projectID,$refSampleList)=@_;
	my %modExperiments;
	my $sthSE=$dbh->prepare("SELECT ID_SPOT,ID_EXPERIMENT FROM SAMPLE WHERE ID_SAMPLE=?");
	my $sthSA=$dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS WHERE ID_SAMPLE=?");
	my $sthDel=$dbh->prepare("DELETE FROM SAMPLE WHERE ID_SAMPLE=?");
	
	foreach my $sampID (@{$refSampleList}) {
		$sthSE->execute($sampID);
		my ($spotID,$expID)=$sthSE->fetchrow_array;
		$modExperiments{$expID}=1 unless $spotID;
		$sthSA->execute($sampID);
		
		deleteMetadata($dbh, "SAMPLE", $sampID);
		
		# Remove analysis
		my @sampAna;
		while (my ($anaID)=$sthSA->fetchrow_array) {push @sampAna,$anaID;}
		&deleteAnalyses($dbh,$projectID,\@sampAna);
		
		# Remove sample
		$sthDel->execute($sampID);
		
		print '.';
		$dbh->commit;
		
	}
	$sthSE->finish;
	$sthSA->finish;
	$sthDel->finish;
	
	foreach my $expID (keys %modExperiments) {
		&updateBrothersPosition('EXPERIMENT',$expID,'SAMPLE');
	}
	$dbh->commit;
}

sub deleteSpots {
	my ($dbh,$projectID,$refSpotList)=@_;
	my $sthSS=$dbh->prepare("SELECT ID_SAMPLE FROM SAMPLE WHERE ID_SPOT=?");
	my @sampleList;
	foreach my $spotID (@{$refSpotList}) {
		$sthSS->execute($spotID);
		while (my ($sampID)=$sthSS->fetchrow_array) {push @sampleList,$sampID;}
	}
	$sthSS->finish;
	&deleteSamples($dbh,$projectID,\@sampleList);

	my $sthDel=$dbh->prepare("DELETE FROM SPOT WHERE ID_SPOT=?");
	foreach my $spotID (@{$refSpotList}) {
		$sthDel->execute($spotID);
	}
	$sthDel->finish;
}


sub delete2DGels {
	my ($dbh,$projectID,$refGelList)=@_;
	my $sthGE=$dbh->prepare("SELECT ID_EXPERIMENT FROM GEL2D WHERE ID_GEL2D=?");
	my $sthGS=$dbh->prepare("SELECT ID_SPOT FROM SPOT WHERE ID_GEL2D=?");
	my $sthDel=$dbh->prepare("DELETE FROM GEL2D WHERE ID_GEL2D=?");
	my (@spotList,%modExperiments);
	foreach my $gelID (@{$refGelList}) {
		$sthGE->execute($gelID);
		my ($expID)=$sthGE->fetchrow_array;
		$modExperiments{$expID}=1;
		$sthGS->execute($gelID);
		while (my ($spotID)=$sthGS->fetchrow_array) {push @spotList,$spotID;}
		&deleteSpots($dbh,$projectID,\@spotList);
		$sthDel->execute($gelID);
		unlink "$promsPath{gel_unix}/project_$projectID/gel2D_$gelID.jpg";
		print '.';
		$dbh->commit;
	}
	$sthGE->finish;
	$sthGS->finish;
	$sthDel->finish;
	foreach my $expID (keys %modExperiments) {
		&updateBrothersPosition('EXPERIMENT',$expID,'SAMPLE');
	}
	$dbh->commit;
}

sub deleteExperiments {
	my ($dbh,$projectID,$refExpList)=@_;
	#<Exploratory Analysis
	my $sthEA=$dbh->prepare("SELECT ID_EXPLORANALYSIS FROM EXPLORANALYSIS WHERE ID_EXPERIMENT=? ORDER BY ID_EXPLORANALYSIS DESC");
	my $sthDelAS=$dbh->prepare("DELETE FROM ANNOTATIONSET WHERE ID_EXPLORANALYSIS=?");
	my $sthDelEAQ=$dbh->prepare("DELETE FROM EXPLORANA_QUANTIF WHERE ID_EXPLORANALYSIS=?");
	my $sthDelPEA=$dbh->prepare("DELETE FROM PARENT_EXPLORANALYSIS WHERE ID_EXPLORANALYSIS=?"); # for motif analysis
	my $sthDelEAA=$dbh->prepare("DELETE FROM EXPLORANA_ANA WHERE ID_EXPLORANALYSIS=?"); # for motif analysis
	my$sthDelEA=$dbh->prepare("DELETE FROM EXPLORANALYSIS WHERE ID_EXPLORANALYSIS=?");
	#<Pathway Analysis
	my $sthPA=$dbh->prepare("SELECT ID_PATHWAY_ANALYSIS FROM PATHWAY_ANALYSIS WHERE ID_EXPERIMENT=? AND ANALYSIS_TYPE='PATHWAY'");
	my $sthDelPAQ=$dbh->prepare("DELETE FROM PATHWAYANA_QUANTIFICATION WHERE ID_PATHWAY_ANALYSIS=?");
	my $sthDelPAA=$dbh->prepare("DELETE FROM PATHWAYANA_ANALYSIS WHERE ID_PATHWAY_ANALYSIS=?");
	my $sthDelPA=$dbh->prepare("DELETE FROM PATHWAY_ANALYSIS WHERE ID_PATHWAY_ANALYSIS=?");
	#<GSEA
	my $sthGSEA=$dbh->prepare("SELECT ID_PATHWAY_ANALYSIS FROM PATHWAY_ANALYSIS WHERE ID_EXPERIMENT=? AND ANALYSIS_TYPE='GSEA'");
	my $sthDelGSEAQ=$dbh->prepare("DELETE FROM PATHWAYANA_QUANTIFICATION WHERE ID_PATHWAY_ANALYSIS=?");
	my $sthDelGSEA=$dbh->prepare("DELETE FROM PATHWAY_ANALYSIS WHERE ID_PATHWAY_ANALYSIS=?");
	#<GO Analysis
	my $sthGo=$dbh->prepare("SELECT ID_GOANALYSIS FROM GO_ANALYSIS WHERE ID_EXPERIMENT=?");
	#<Design-based Quantification
	my $sthQD=$dbh->prepare("SELECT ID_QUANTIFICATION FROM QUANTIFICATION WHERE ID_DESIGN=?");
	#<Design, Condition & obsCondition (OBSERVATION,OBS_MODIFICATION & internal-ana quantifs deleted at analysis level)
	my $sthED=$dbh->prepare("SELECT ID_DESIGN FROM DESIGN WHERE ID_EXPERIMENT=?");
	my $sthDelOC=$dbh->prepare("DELETE OEC FROM OBS_EXPCONDITION OEC INNER JOIN EXPCONDITION EC ON OEC.ID_EXPCONDITION=EC.ID_EXPCONDITION WHERE ID_DESIGN=?");
	my $sthDelC=$dbh->prepare("DELETE FROM EXPCONDITION WHERE ID_DESIGN=?");
	my $sthDelD=$dbh->prepare("DELETE FROM DESIGN WHERE ID_DESIGN=?");
	# User lock experiment
	my $sthDelULE=$dbh->prepare("DELETE FROM USER_EXPERIMENT_LOCK WHERE ID_EXPERIMENT=?");
	#<Exp items,Experiment
	my $sthES=$dbh->prepare("SELECT ID_SAMPLE FROM SAMPLE WHERE ID_EXPERIMENT=?");
	my $sthEG=$dbh->prepare("SELECT ID_GEL2D FROM GEL2D WHERE ID_EXPERIMENT=?");
	my $sthDelE=$dbh->prepare("DELETE FROM EXPERIMENT WHERE ID_EXPERIMENT=?");

	foreach my $expID (@{$refExpList}) {
		##<Metadata
		deleteMetadata($dbh, "EXPERIMENT", $expID);
		
		#<User lock
		$sthDelULE->execute($expID);
		
		print '/';
		##<Exploratory analyses & annotationsets
		$sthEA->execute($expID);
		while (my ($eaID)=$sthEA->fetchrow_array) {
			$sthDelAS->execute($eaID);
			$sthDelEAQ->execute($eaID);
			$sthDelPEA->execute($eaID);
			$sthDelEAA->execute($eaID);
			$sthDelEA->execute($eaID);
			rmtree("$promsPath{explorAna}/project_$projectID/$eaID");
			print '.';
		}

		##<Pathway analyses
		$sthPA->execute($expID);
		while (my ($paID)=$sthPA->fetchrow_array) {
			$sthDelPAQ->execute($paID);
			$sthDelPAA->execute($paID);
			$sthDelPA->execute($paID);
			rmtree("$promsPath{pathAna}/project_$projectID/$paID");
			print '.';
		}

		##<GSEA analyses
		$sthGSEA->execute($expID);
		while (my ($gseaID) = $sthGSEA->fetchrow_array) {
			$sthDelGSEAQ->execute($gseaID);
			$sthDelGSEA->execute($gseaID);
			rmtree("$promsPath{gsea}/project_$projectID/gsea_$gseaID");
			print '.';
		}

		##<GO analyses
		$sthGo->execute($expID);
		while (my ($goID)=$sthGo->fetchrow_array) {
			&goAnalysis::deleteGOAnalysis($dbh,$projectID,$goID);
			print '.';
		}

		##<Designs, obsCondition, conditions & design-based quantifications
		print '/';
		$sthED->execute($expID);
		while (my ($desID)=$sthED->fetchrow_array) {
			$sthQD->execute($desID);
			while (my ($qID)=$sthQD->fetchrow_array) {
				&promsQuantif::deleteQuantification($dbh,$projectID,$qID,{VERBOSE=>1});
				print '.';
			}
			$sthDelOC->execute($desID);
			#$sthDelCQ->execute($desID);
			$sthDelC->execute($desID);
			$sthDelD->execute($desID);
			print '.';
		}

		##<Sample & children
		print '/';
		my @expSamp;
		$sthES->execute($expID);
		while (my ($sampID)=$sthES->fetchrow_array) {push @expSamp,$sampID;}
		&deleteSamples($dbh,$projectID,\@expSamp);

		##<2DGels & children
		my @exp2dGels;
		$sthEG->execute($expID);
		while (my ($gelID)=$sthEG->fetchrow_array) {push @exp2dGels,$gelID;}
		&delete2DGels($dbh,$projectID,\@exp2dGels) if scalar @exp2dGels;

		##>Experiment
		$sthDelE->execute($expID);
		$dbh->commit;
	}
	$sthEA->finish;
	$sthDelEAA->finish;
	$sthDelPEA->finish;
	$sthDelAS->finish;
	$sthDelEA->finish;
	$sthDelEAQ->finish;
	$sthPA->finish;
	$sthDelPAQ->finish;
	$sthDelPAA->finish;
	$sthDelPA->finish;
	$sthGSEA->finish;
	$sthDelGSEAQ->finish;
	$sthDelGSEA->finish;
	$sthGo->finish;
	$sthQD->finish;
	$sthED->finish;
	$sthDelOC->finish;
	$sthDelC->finish;
	$sthDelD->finish;
	$sthES->finish;
	$sthEG->finish;
	$sthDelE->finish;
	
	&updateBrothersPosition('PROJECT',$projectID,'EXPERIMENT');
	print '.';
}

###############################################
####<Deleting validated data from analysis>####
###############################################
# !!!IMPORTANT: Analysis deletablity is controlled by selectAnalysis.cgi (callType=delete) & selectOption.cgi
# Deletion should NOT be allowed at non-Experiment level:
#	-if quantif, except for internal labeled peptide quantif (SILAC from MSF)
#	-if used for GO & Pathway Analyses
#	-OK if Comparison or Observation: handled by deleteValidData (WARNING: possible empty comparison after deletion)
sub deleteValidData {
	my ($dbh,$projectID,$anaID,$fileFormat,$dataFile)=@_;
	my %promsPath=&promsConfig::getServerInfo; # <- promsConfig declared in main script

	####<Deleting internal analysis quantification data>####
	my $sthQ=$dbh->prepare("SELECT Q.ID_QUANTIFICATION FROM ANA_QUANTIFICATION AQ INNER JOIN QUANTIFICATION Q ON AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION WHERE ID_DESIGN IS NULL AND ID_ANALYSIS=$anaID AND FOCUS=? ORDER BY Q.ID_QUANTIFICATION DESC");
	#my $sthAnaCount=$dbh->prepare("SELECT COUNT(ID_ANALYSIS) FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=?");
	#my $sthDelQPep=$dbh->prepare("DELETE PQ FROM PEPTIDE_QUANTIFICATION PQ INNER JOIN PEPTIDE P ON PQ.ID_PEPTIDE=P.ID_PEPTIDE WHERE ID_QUANTIFICATION=? AND ID_ANALYSIS=$anaID");
	#my $sthDelQAna=$dbh->prepare("DELETE FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=? AND ID_ANALYSIS=$anaID");
	#my $sthDelQRT=$dbh->prepare("DELETE FROM QUANTIF_REFRT WHERE ID_QUANTIFICATION=?");

	#my %childQuantifList;
	###my @quantifList;
	foreach my $focus ('protein','peptide') { # peptide: quantif file will be deleted even if multi-ana quantif (eg. MassChroQ,MaxQuant) <- assumes deletion is trigered by Experiment-wide deletion
		$sthQ->execute($focus);
		while (my ($quantifID)=$sthQ->fetchrow_array) { # All deletion steps handled by &promsQuantif::deleteQuantification
			#if ($focus eq 'peptide') {
			#	$sthAnaCount->execute($quantifID);
			#	my ($numAna)=$sthAnaCount->fetchrow_array;
			#	if ($numAna > 1) { # Other analyses involved
			#		$sthDelQPep->execute($quantifID);
			#		$sthDelQAna->execute($quantifID);
			#		print '.';
			#		next;
			#	}
			#	$sthDelQRT->execute($quantifID);
			#}
			&promsQuantif::deleteQuantification($dbh,$projectID,$quantifID); # protein: child-quantifs should have been deleted
			print '.';
		}
	}
	$sthQ->finish;
	#$sthDelQPep->finish;
	#$sthDelQAna->finish;
	#$sthDelQRT->finish;

	####<Deleting reference RT link>####
	$dbh->do("DELETE FROM ANALYSIS_REFRT WHERE ID_ANALYSIS=$anaID");

	####<Observations>#### OBS_EXPCONDITION,OBS_MODIFICATION,OBSERVATION
	$dbh->do("DELETE OEC FROM OBS_EXPCONDITION OEC INNER JOIN OBSERVATION O ON OEC.ID_OBSERVATION=O.ID_OBSERVATION WHERE O.ID_ANALYSIS=$anaID"); # just to be safe
	$dbh->do("DELETE OM FROM OBS_MODIFICATION OM INNER JOIN OBSERVATION O ON OM.ID_OBSERVATION=O.ID_OBSERVATION WHERE O.ID_ANALYSIS=$anaID");
	$dbh->do("DELETE FROM OBSERVATION WHERE ID_ANALYSIS=$anaID");

	####<Deleting Peptide modification link>####
	$dbh->do("DELETE PM FROM PEPTIDE_MODIFICATION PM INNER JOIN PEPTIDE P ON PM.ID_PEPTIDE=P.ID_PEPTIDE WHERE P.ID_ANALYSIS=$anaID") || die $dbh->errstr;

	####<Deleting Analysis comparison link>####
	$dbh->do("DELETE FROM ANA_COMPARISON WHERE ID_ANALYSIS=$anaID"); # called from parent item because cannot delete specific ANA if it has children

	####<Deleting GO Analysis link>#### Allowed only at experiment level!!!
	$dbh->do("DELETE FROM GOANA_ANALYSIS WHERE ID_ANALYSIS=$anaID");

	####<Deleting peptides from both peptide tables>####
	$dbh->do("DELETE FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_ANALYSIS=$anaID");
	$dbh->do("DELETE FROM PEPTIDE WHERE ID_ANALYSIS=$anaID");

	###<Deleting peptide file (new procedure)>###
	#remove_tree("$promsPath{peptide}/proj_$projectID/ana_$anaID");
	if ($fileFormat =~ /\.PDM/) { # Proteome Discoverer
		(my $msfFile=$dataFile)=~s/_[0-9]*\.*[0-9]*.pdm/\.msf/;
		if (!-l "$promsPath{peptide}/proj_$projectID/ana_$anaID/$msfFile") { # if it is not a symlink
			my ($moveOK,$filePath)=(0,'');
			foreach my $directory (glob("$promsPath{peptide}/proj_$projectID/*")) {
				next if $directory eq "$promsPath{peptide}/proj_$projectID/ana_$anaID";
				if (-l "$directory/$msfFile") {
					if ($moveOK) {
						unlink("$directory/$msfFile");
						symlink($filePath,"$directory/$msfFile");
					}
					else{
						unlink("$directory/$msfFile");
						move("$promsPath{peptide}/proj_$projectID/ana_$anaID/$msfFile","$directory/$msfFile");
						$moveOK=1;
						$filePath="$directory/$msfFile";
					}
				}
			}
		}
	}
	rmtree("$promsPath{peptide}/proj_$projectID/ana_$anaID");

	####<Deleting proteins>####
	my %modifMasterProteins; # records list of affected master proteins
	###<Queries>###
	my $sthAnaProt=$dbh->prepare("SELECT PROTEIN.ID_PROTEIN,ID_MASTER_PROTEIN FROM ANALYSIS_PROTEIN,PROTEIN WHERE PROTEIN.ID_PROTEIN=ANALYSIS_PROTEIN.ID_PROTEIN AND ID_ANALYSIS=$anaID");
	my $sthNumAna=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=?");
	my $sthDelProt1=$dbh->prepare("DELETE FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=$anaID AND ID_PROTEIN=?");
	my $sthDelProt2=$dbh->prepare("DELETE MS FROM MODIFICATION_SITE MS INNER JOIN CATEGORY_PROTEIN CP ON MS.ID_CATEGORY_PROTEIN=CP.ID_CATEGORY_PROTEIN AND CP.ID_PROTEIN=?");
	my $sthDelProt3=$dbh->prepare("DELETE FROM CATEGORY_PROTEIN WHERE ID_PROTEIN=?");
	my $sthDelProt4=$dbh->prepare("DELETE FROM PROTEIN WHERE ID_PROTEIN=?");

	###<Fetching list of proteins in analysis>###
	$sthAnaProt->execute;
	while (my($protID,$masterProtID)=$sthAnaProt->fetchrow_array) { # reference to an array;

		##<Checking if protein is unique to analysis>##
		$sthNumAna->execute($protID);
		my ($numAna)=$sthNumAna->fetchrow_array;

		##<Deleting protein>##
		$sthDelProt1->execute($protID); # ANALYSIS_PROTEIN
		if ($numAna==1) { # protein is unique to analysis
			$modifMasterProteins{$masterProtID}=1 if $masterProtID;
			$sthDelProt2->execute($protID); # MODIFICATION_SITE
			$sthDelProt3->execute($protID); # CATEGORY_PROTEIN
			$sthDelProt4->execute($protID); # PROTEIN
		}
	}
	$sthAnaProt->finish;
	$sthNumAna->finish;
	$sthDelProt1->finish;
	$sthDelProt2->finish;
	$sthDelProt3->finish;
	$sthDelProt4->finish;

	##<Deleting unused master proteins
	&promsMod::deleteUnusedMasterProteins($dbh,\%modifMasterProteins);
}

sub updateBrothersPosition {
	my ($parItem,$parID,$childItem) =@_;
	#my $item =($parItem eq 'SAMPLE')? 'ANALYSIS' : ($parItem eq 'EXPERIMENT')?'SAMPLE':($parItem eq 'PROJECT')?'EXPERIMENT':die('illegal call from update_position');
	#my @items=&promsMod::getItemChild($parItem);
	#my $item=$items[0];
	my $spotStrg=($childItem eq 'SAMPLE' && $parItem eq 'EXPERIMENT')? ' ID_SPOT IS NULL AND' : '';
	my $sthSel=$dbh->prepare("SELECT ID_$childItem FROM $childItem WHERE$spotStrg ID_$parItem=$parID ORDER BY DISPLAY_POS ASC");
	my $sthUp=$dbh->prepare("UPDATE $childItem SET DISPLAY_POS=? WHERE ID_$childItem=?");
	$sthSel->execute;
	my $displayPos=0;
	while (my ($itID)=$sthSel->fetchrow_array) {
		$sthUp->execute(++$displayPos,$itID);
	}
	$sthSel->finish;
	$sthUp->finish;
}

####>Revision history<####
# 2.6.9 [FEATURE] Add deletion of GSEA with experiment (VL 18/11/20)
# 2.6.8 [ENHANCEMENT] Added verbose=1 to quantification deletion in Experiment-wide deletion context to prevent server timeout (PP 22/07/20)
# 2.6.7 [BUGFIX] Delete user_lock experiment on project item deletion (VS 09/04/20)
# 2.6.6 [BUGFIX] Delete metadata on project item deletion (VS 16/12/19)
# 2.6.5 [ENHANCEMENT] Simplify metadata deletion (VS 15/11/19)
# 2.6.4 [FIX] Properly delete a file related to a metadata when it is deleted/its project item is deleted (VS 05/06/19)
# 2.6.3 Change project data path for metadata path (VS 11/06/19)
# 2.6.2 Delete Metadata on item deletion (VS 05/06/19)
# 2.6.1 [Fix] bugs in Experiment global deletion (PP 27/06/18)
# 2.6.0 Peptide quantification deletion steps moved to &promsQuanti::deleteQuantification (PP 11/05/18)
# 2.5.3 Compatible with MODIFICATION_SITE table used for list of modification sites (PP 28/07/17)
# 2.5.2 Deletes only analysis-specific data from peptide quantif if multi-(ana/sample) quantif (PP 03/01/17)
# 2.5.1 Update deleteAnalyses for Swath (delete ANALYSIS_SWATH_LIB) (MLP 14/06/16)
# 2.5.0 Handles refrence RT link (PP 21/03/16)
# 2.4.9 Update deleteValidData for split-files (GA 30/03/15)
# 2.4.8 Change for MSF split-mode and .ana update & minor bug fix (GA,PP 12/03/15)
# 2.4.7 Minor syntax error fix (PP 20/02/15)
# 2.4.6 Change in quantification deletion & handles Pathway Analyses (PP 13/11/14)
# 2.4.5 Updated for OBS_EXPCONDITION, Exploratory Analyses and BioSamples (PP 08/08/14)
# 2.4.4 Uses rmtree instead of remove_tree (PP 10/03/14)
# 2.4.3 Update for List comparison (PP 21/01/14)
# 2.4.2 Added File::Copy for &promsMod::removeFromMultiAna (PP 08/11/13)
# 2.4.1 Fix forgotten quanti directory deletion (PP 11/10/13)
# 2.4.0 Updated gel path (project_x subdir) (FY 03/09/13)
# 2.3.9 Debug for quantification deletion (PP 05/07/13)
# 2.3.8 Updated for OBSERVATION & OBS_MODIFICATION tables (PP 02/07/13)
# 2.3.7 Deletes relevant PTMs upon project deletion (PP 28/06/13)
# 2.3.6 Minor bug correction in deleteExperiments (GA 07/06/13)
# 2.3.5 Add PEPTIDE_MODIFICATION deletion in &deleteValidData (GA 02/05/13)
# 2.3.4 Commented PROJECT_MODIFICATION management in &deleteAnalyses (PP 26/04/13)
# 2.3.3 Update analysis based on new modification tables (GA 16/04/13)
# 2.3.2 Remove deleting in obsolete ANA_QUANTIF_PROT table (FY 09/04/13)
# 2.3.1 &promsMod::removeFromMultiAna not called for fully validated analysis (PP 20/03/13)
# 2.3.0 Major update to handle deletion of new analysis and experiment-dependent items (PP 13/03/13)
# 2.2.5 Merge 2.2.GA & 2.2.3b (PP 02/04/12)
# 2.2.3b Text update during deletion (PP 28/03/12)
# 2.2.3 Deleting validation history before deleting analyses in database (FY 22/03/11)
