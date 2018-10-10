#!/usr/local/bin/perl -w

################################################################################
# storeSpectrum.cgi     1.4.4                                                  #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Stores a reference spectrum for a peptide sequence.                          #
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
use strict ;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use XML::Simple; # needed by &promsMod::extractData

#print header; warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my $date=strftime("%Y-%m-%d %H:%M:%S", localtime);

#############################
####>Connecting to ProMS<####
#############################
my $dbh=&promsConfig::dbConnect;

####################
####>Parameters<####
####################
my $action=param('ACT');
my $refSpecID=param('RID'); # RID <= 0 (-ID_SPECTRUM)
&deleteReference if $action eq 'delete';
my $objectID=param('ID');
my $rank=param('RANK');
my $call=param('CALL');
my $fileType= param('TYPE');


#############################
####>Fetching query info<####
#############################

my ($projectID,$analysisID,$validStatus,$fileName,$dataFile,$msfFile,$queryNum,$searchRank,$sequence,$massObs,$massExp,$massCalc,$score,$delta,$missCut,$sub,$comments,$instrument,$elutionTime);
if ($call eq 'pep' || $call eq 'val') { # validated spectrum, called from sequenceView if 'pep' or listMatchedAnalyses.cgi if 'val'
	($analysisID,$queryNum,$elutionTime,$searchRank,$sequence,$massObs,$massExp,$massCalc,$score,$delta,$missCut)=$dbh->selectrow_array("SELECT ID_ANALYSIS,QUERY_NUM,ELUTION_TIME,SEARCH_RANK,PEP_SEQ,MR_OBS,MR_EXP,MR_CALC,SCORE,MR_DELTA,MISS_CUT FROM PEPTIDE WHERE ID_PEPTIDE=$objectID");
	$searchRank=$rank unless $searchRank; # mixed decoy DB

	($validStatus,$fileName)=$dbh->selectrow_array("SELECT VALID_STATUS,DATA_FILE FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
	$fileName =~s/\.xml/\.pgf/  if $fileType=~/(PHENYX|MASCOT|PARAGON)\.XML/; #replace file to peaklist file in case of phenyx
	$projectID=&promsMod::getProjectID($dbh,$analysisID,'analysis');
	#my $pepFileRoot=sprintf "P%06d",$analysisID;
	#$dataFile="$promsPath{peptide}/proj_$projectID/$pepFileRoot.dat" if $fileType eq 'MASCOT.DAT';
	#$dataFile="$promsPath{peptide}/proj_$projectID/$pepFileRoot.pgf" if ($fileType eq 'PHENYX.XML' ||  $fileType eq 'MASCOT.XML');
	#unless (-e $dataFile) {
	#	$dataFile="$promsPath{peptide}/proj_$projectID/$analysisID"."_$fileName"; # new data file conservation procedure (01/06/11)
	#	$dataFile=~s/\.xml/\.pgf/ if $fileType=~/\.XML/;
	#}
	if ($validStatus==2) {
		$dataFile="$promsPath{peptide}/proj_$projectID/ana_$analysisID/$fileName";
		($msfFile=$dataFile) =~ s/\_\d+\.pdm/\.msf/ if $fileType=~/\.PDM/;
		$dataFile=~s/\.dat/_min\.dat/ unless -e $dataFile; # assume minimal dat file
	}
	else {
		$dataFile="$promsPath{valid}/ana_$analysisID/$fileName";
		if ($fileType=~/\.PDM/) {
			(my $msfFileName=$fileName) =~ s/\_\d+\.pdm/\.msf/;
			$msfFile="$promsPath{valid}/multi_ana/proj_$projectID/$msfFileName";
		}
	}
}
elsif ($call eq 'rank' || $call eq 'ana' || $call eq 'seq') { # call=rank or ana (called during validation)
	($analysisID,$queryNum,$elutionTime,$fileName,my $massData,my $pepInfo,$analysisID)=$dbh->selectrow_array("SELECT A.ID_ANALYSIS,Q.QUERY_NUM,Q.ELUTION_TIME,A.DATA_FILE,Q.MASS_DATA,Q.INFO_PEP$rank,Q.ID_ANALYSIS FROM QUERY_VALIDATION Q,ANALYSIS A WHERE Q.ID_ANALYSIS=A.ID_ANALYSIS AND ID_QUERY=$objectID");
	$fileName =~s/\.xml/\.pgf/  if ($fileType eq 'PHENYX.XML' || $fileType eq 'MASCOT.XML') ; #replace file to pealist file in case of phenyx
	$dataFile="$promsPath{valid}/ana_$analysisID/$fileName";
	if ($fileType=~/\.PDM/) {
		$projectID=&promsMod::getProjectID($dbh,$analysisID,'analysis');
		(my $msfFileName=$fileName) =~ s/\_\d+\.pdm/\.msf/;
		$msfFile="$promsPath{valid}/multi_ana/proj_$projectID/$msfFileName";
	}
	($massExp)=($massData=~/EXP=(\d+\.*\d*)/);
	($massObs)=($massData=~/OBS=(\d+\.*\d*)/);
	($searchRank)=($pepInfo=~/SRK=(\d+)/); $searchRank=$rank unless $searchRank; # mixed decoy DB
	($sequence)=($pepInfo=~/SEQ=(\w+)/);
	($score)=($pepInfo=~/SC=(-?\d+\.*\d*)/);
	($delta)=($pepInfo=~/DELT=(-*\d+\.*\d*)/);
	($massCalc)=($pepInfo=~/CALC=(\d+\.*\d*)/);
	($missCut)=($pepInfo=~/MIS=(\d)/);
	($sub)=($pepInfo=~/SUBST=([^,]+)/);
	if ($pepInfo=~ /COM=(.+)/) {
		$comments = $1;
	};
}

else {  # Problem
	print header(-'content-encoding'=>'no',-charset=>'utf-8'); # start_html,"\n"; # start_html required to force update
	warningsToBrowser(1) ;
	print ("<H3> Unknown call command \"$call\" <H3> <BR>");
	die;
}
($instrument) = $dbh->selectrow_array("SELECT INSTRUMENT FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");

###################################
####>Fetching (new) spectrumID<####
###################################
my $spectrumID;
if ($refSpecID==0) { # no reference spectrum => create new primary key
	($spectrumID)=$dbh->selectrow_array("SELECT MAX(ID_SPECTRUM) FROM SPECTRUM");
	$spectrumID++;
	$dbh->do("INSERT INTO SPECTRUM (ID_SPECTRUM) VALUES ($spectrumID)") || die $dbh->errstr; # reserving/protecting ID
}
else { # old reference spectrum exists => get its pk and overwrite
	$spectrumID=$refSpecID;
	$dbh->do("DELETE FROM SPECTRUM_MODIFICATION WHERE ID_SPECTRUM=$spectrumID") || die $dbh->errstr; # erase previous PTMs
}
my $zeroString='';
my $zeroNum=6-length($spectrumID);
foreach my $i (1..$zeroNum) {$zeroString.='0';}
#my $specFileName=($fileType eq 'MASCOT.DAT')? "S$zeroString$spectrumID.dat" : ($fileType eq 'PHENYX.XML' || $fileType eq 'MASCOT.XML' )? "S$zeroString$spectrumID.pgf" : "toto.dat";
my $specFileName=($fileType eq 'MASCOT.DAT')? "S$zeroString$spectrumID.dat" : ($fileType eq 'PHENYX.XML' || $fileType eq 'MASCOT.XML' || $fileType eq 'PARAGON.XML' )? "S$zeroString$spectrumID.pgf" : ($fileType =~ /\.PDM/ )? "S$zeroString$spectrumID.pdm" : "S$zeroString$spectrumID.xxx";
my $specFile="$promsPath{spectrum}/$specFileName";


#######################################
####>Extracting data from dataFile<####
#######################################
my %sectionField=('parameters'=>'all','masses'=>'all',"query$queryNum"=>'all');
$sectionField{'summary'}{"qmass$queryNum="}=1;
$sectionField{'summary'}{"qexp$queryNum="}=1;
$sectionField{'summary'}{"qmatch$queryNum="}=1;
$sectionField{'peptides'}{"q$queryNum"."_p$searchRank="}=1;
if ($action eq 'store') {
	if ($fileType =~ /\.PDM/) {
		my $extSpectrumID;
		if ($validStatus>=1) { # >=1 because $objectID is ID_PEPTIDE!!!
			my ($pepData)=$dbh->selectrow_array("SELECT DATA FROM PEPTIDE WHERE ID_PEPTIDE=$objectID");
			($extSpectrumID)=($pepData=~/EXT_SPECTRUMID=(\d+)/);
		}
		else {($extSpectrumID)=$dbh->selectrow_array("SELECT EXT_SPECTRUMID FROM QUERY_VALIDATION WHERE ID_QUERY=$objectID");}
		my ($projectID)=&promsMod::getProjectID($dbh,$analysisID,'analysis') unless $projectID; # should be defined already
		$sectionField{'projectid'}=$projectID;
		$sectionField{'extspectrumid'}=$extSpectrumID;
	}elsif ($fileType eq 'PARAGON.XML') {
		my ($elutionTime)=$dbh->selectrow_array("SELECT ELUTION_TIME FROM QUERY_VALIDATION WHERE ID_QUERY=$objectID");
		my ($extSpectrumID)=($elutionTime=~/sp([\d+|\.]*);/ );
		$sectionField{'extspectrumid'}=$extSpectrumID;
		$sectionField{'searchRank'}=$searchRank;
		$sectionField{'queryNum'}=$queryNum;
	}
}
&promsMod::extractData(\%sectionField,$dataFile,$specFile,$fileType,$msfFile);


########################################
####>Storing data in table SPECTRUM<####
########################################
my $sthUpRS=$dbh->prepare(qq
|UPDATE SPECTRUM SET
	QUERY_NUM=$queryNum,PEP_RANK=$rank,SEARCH_RANK=$searchRank,PEP_SEQ='$sequence',SPEC_FILE='$specFileName',DATA_FILE='$fileName',SCORE=$score,MISS_CUT=$missCut,MR_OBS=$massObs,MR_EXP=$massExp,MR_CALC=$massCalc,MR_DELTA=$delta,UPDATE_USER='$userID',UPDATE_DATE='$date',
	ELUTION_TIME=?,SUBST=?,COMMENTS=?,INSTRUMENT=?
	WHERE ID_SPECTRUM=$spectrumID|);
#$dbh->do("UPDATE SPECTRUM SET QUERY_NUM=$queryNum,ELUTION_TIME='$elutionTime',PEP_RANK=$rank,PEP_SEQ='$sequence',SPEC_FILE='$specFileName',DATA_FILE='$fileName',SCORE=$score,MISS_CUT=$missCut,MR_OBS=$massObs,MR_EXP=$massExp,MR_CALC=$massCalc,MR_DELTA=$delta,UPDATE_USER='$userID',UPDATE_DATE='$date' WHERE ID_SPECTRUM=$spectrumID") || die $dbh->errstr;

$sthUpRS->execute($elutionTime,$sub,$comments,$instrument); # possible undef values
###> New VMOD handling
my $sthGetVMods= ($call eq 'pep' || $call eq 'val')?$dbh->prepare("SELECT PM.ID_MODIFICATION,AM.MODIF_TYPE,AM.SPECIFICITY,POS_STRING FROM PEPTIDE_MODIFICATION PM, MODIFICATION M, ANALYSIS_MODIFICATION AM WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND PM.ID_MODIFICATION=M.ID_MODIFICATION AND ID_PEPTIDE=$objectID AND ID_ANALYSIS=$analysisID"):
					$dbh->prepare("SELECT QM.ID_MODIFICATION,AM.MODIF_TYPE,AM.SPECIFICITY,POS_STRING FROM QUERY_MODIFICATION QM, MODIFICATION M, ANALYSIS_MODIFICATION AM WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND QM.ID_MODIFICATION=M.ID_MODIFICATION AND ID_QUERY=$objectID AND PEP_RANK=$rank AND ID_ANALYSIS=$analysisID");

$sthGetVMods->execute;

my $sthInsMods=$dbh->prepare("INSERT INTO SPECTRUM_MODIFICATION (ID_MODIFICATION,ID_SPECTRUM,MODIF_TYPE,SPECIFICITY,POS_STRING) VALUES (?,?,?,?,?) ");
while (my ($modID,$modifType,$specificity,$posString) = $sthGetVMods->fetchrow_array ) {
	$sthInsMods->execute($modID,$spectrumID,$modifType,$specificity,$posString);
}
$sthGetVMods->finish;

my $sthGetFMods=$dbh->prepare("SELECT ID_MODIFICATION,SPECIFICITY FROM ANALYSIS_MODIFICATION WHERE ID_ANALYSIS=$analysisID AND MODIF_TYPE='F'");
$sthGetFMods->execute;
while (my ($modID,$specificity) = $sthGetFMods->fetchrow_array ) {
	$sthInsMods->execute($modID,$spectrumID,'F',$specificity,undef);
}
$sthGetFMods->finish;

$sthInsMods->finish;


#if ($varModsString) {
#	$dbh->do("UPDATE SPECTRUM SET VAR_MOD='$varModsString' WHERE ID_SPECTRUM=$spectrumID") || die $dbh->errstr;
#}
#if ($comments) {
#	$comments=$dbh->quote ($comments) ;
#	$dbh->do("UPDATE SPECTRUM SET COMMENTS=$comments WHERE ID_SPECTRUM=$spectrumID") || die $dbh->errstr;
#}
#if ($instrument) {
#	$instrument=$dbh->quote ($instrument) ;
#	$dbh->do("UPDATE SPECTRUM SET INSTRUMENT=$instrument WHERE ID_SPECTRUM=$spectrumID") || die $dbh->errstr;
#}

$dbh->commit;
$dbh->disconnect;
#exit; #debug

#######################
####>Starting HTML<####
#######################
my $rankID=($call eq 'pep')? 'pep_' : 'seq_';
$rankID.="$objectID"."_$queryNum"."_$rank";
print redirect("./drawSpectrum.cgi?RID=$rankID&CALL=$call");


########################################<<<SUBROUTINE>>>#################################

#######################################
####<Deleting a reference spectrum<####
#######################################
sub deleteReference {

	##<Deleting spectrum file
	my $zeroString='';
	my $zeroNum=6-length($refSpecID);
	for my $i (1..$zeroNum) {$zeroString.='0';}
	my $specFileName="S$zeroString$refSpecID.dat";
	unlink "$promsPath{spectrum}/$specFileName";

	##<Deleting entry in table SPECTRUM
	$dbh->do("DELETE FROM SPECTRUM_MODIFICATION WHERE ID_SPECTRUM=$refSpecID");
	$dbh->do("DELETE FROM SPECTRUM WHERE ID_SPECTRUM=$refSpecID");

	$dbh->commit;
	$dbh->disconnect;

	print header(-'content-encoding'=>'no',-charset=>'utf-8'); # start_html,"\n"; # start_html required to force update
	warningsToBrowser(1);
	print qq
|<HTML>
<BODY>
<SCRIPT LANGUAGE="JavaScript">
//	alert('delete: '+parent.location.href);
	parent.location.replace(parent.location.href); // reload() creates crazy loop in refSpecFrame!!!
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

####>Revision history<####
# 1.4.4 system command removal & added 'use XML::Simple' (PP 13/11/13)
# 1.4.3 Remove VAR_MOD from script and store new modifications (GA 23/05/13)
# 1.4.2 Update for substitution information (GA 21/03/13)
# 1.4.1 Corrects msf file extraction from fully validated analyses (PP 18/03/13)
# 1.4.0 Uses new data file path for valid analysis (PP 06/03/13)
# 1.3.9 Make this script work for PDM files (GA 13/11/12)
# 1.3.8 full datafile conservation management<BR>& no longer uses DATA_FILE from QUERY_VALIDATION table (PP 14/06/11)
