#!/usr/local/bin/perl -w

################################################################################
# mapProteinIdentifiers.pl        1.5.0                                        #
# Authors: P. Poullet, G. Arras & F. Yvon (Institut Curie)                     #
# Contact: myproms@curie.fr                                                    #
# Maps protein identifiers over the internet.                                  #
# Called by send2Biologist.cgi.                                                #
# Runs in background.                                                          #
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
use POSIX qw(strftime); # to get the time
use XML::Simple;
use File::Path qw(mkpath); # make_path
use LWP::UserAgent;
use promsConfig;
use promsMod;
#exit;
#############################
####>Fetching parameters<####
#############################
my ($userID,$itemIdString,$dbType,$update)=@ARGV; ## $itemIdString = $projectID in case of force update ($update='force' ONLY called project-wide!!!)
$userID=&promsMod::cleanParameters($userID);
$dbType='AUTO' unless $dbType;
$update="" unless $update;

########################
####>Immediate fork<####
########################
#unless ($update) { # $update only for update to myPROMS v.3.0 (called by updateIdentifierMapping.pl)
if (!$update  || $update eq 'force') {
	my $childPid = fork;
	if ($childPid) { # parent here
		#sleep 2;
		exit; # returns control to caller script
	}
}

##############>CHILD FROM NOW ON -------------------------------------->

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo($userID);
#mkdir $promsPath{'tmp'} unless -e $promsPath{'tmp'};
mkdir $promsPath{'logs'} unless -e $promsPath{'logs'};

####>Disconnecting from parent's STDs<####
my $year=(!$update || $update eq 'force')? strftime("%Y",localtime) : 'update';
my $logFile="$promsPath{logs}/mapping_$year.log";
open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
open STDERR, ">>$logFile";
open STDOUT, ">>$logFile";
$| = 1;

#########################################
####>Waiting for previous job to end<####
#########################################
my $jobFlagFile="$promsPath{logs}/mapping.job";
my $timeStamp;
if (-e $jobFlagFile) {
	open(PREV_JOB,$jobFlagFile);
	my @jobInfo=<PREV_JOB>;
	close PREV_JOB;
	$timeStamp=(split(/\s/,$jobInfo[0]))[0];
	unless ($timeStamp) {
		unlink($jobFlagFile);
	}
}
my $prevLogSize=(stat($logFile))[7];
sleep 5;
while (-e $jobFlagFile) {
	my $newSize=(stat($logFile))[7];
	if ($newSize > $prevLogSize) { # file size
		$prevLogSize=$newSize;
		$timeStamp=time;
	}
	elsif (time-$timeStamp > 900) { # more than 15 min w/o change => assume previous job has failed!
		unlink($jobFlagFile);
		last;
	}
	sleep 10;
}

####>Starting new job<####
my $jobStartTime=strftime("%Y%m%d%H%M%S",localtime);
my $itemStrg=($update eq 'force')? 'proj' : 'ana';
print "\n>JOB START $jobStartTime user=$userID [$dbType] $itemStrg=$itemIdString\n";


##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect($userID);

###>New job flag file<###
my $projectID=($update eq 'force')? &promsMod::cleanNumericalParameters($itemIdString) : &promsMod::getProjectID($dbh,(split(',',$itemIdString))[0],'ANALYSIS');
open(JOB,">$jobFlagFile");
print JOB time." $jobStartTime $userID $projectID";
close JOB;


##############################
####>Preloading some data<####
##############################
my (%species,%identifierTypes);
my $sthSp=$dbh->prepare("SELECT ID_SPECIES,TAXONID FROM SPECIES");
$sthSp->execute;
while (my($spID,$taxID)=$sthSp->fetchrow_array) {
	$species{$taxID}=$spID;
}
$sthSp->finish;

my $sthIT=$dbh->prepare("SELECT ID_IDENTIFIER,CODE FROM IDENTIFIER");
$sthIT->execute;
while (my($identID,$dbCode)=$sthIT->fetchrow_array) {
	$identifierTypes{$dbCode}=$identID;
}
$sthIT->finish;


my @analysisList;

my (%forceRemovedProteins,%forceUnusedMasters,%forceUpdateProteins,%backupUniprotAC);
if ($update eq 'force') {

	###>Delete previous mapping older than 1 month<### (PP 15/03/18, 20/09/18)
	print "+Removing link to mapping older than 1 month...";
	#my $sthMapProt=$dbh->prepare("SELECT ID_PROTEIN,ID_MASTER_PROTEIN,PROT_SEQ FROM PROTEIN WHERE ID_PROJECT=? AND ID_MASTER_PROTEIN IS NOT NULL");
	my $sthMapProt=$dbh->prepare("SELECT P.ID_PROTEIN,P.ID_MASTER_PROTEIN,P.PROT_SEQ FROM PROTEIN P,MASTER_PROTEIN MP WHERE P.ID_MASTER_PROTEIN=MP.ID_MASTER_PROTEIN AND ID_PROJECT=$projectID AND MP.UPDATE_DATE-DATE_SUB(NOW(),INTERVAL 1 MONTH) < 0");
	#my $sthUpProt1=$dbh->prepare("UPDATE PROTEIN SET ID_MASTER_PROTEIN=NULL,PROT_SEQ='-',PROT_LENGTH=0 WHERE ID_PROTEIN=?"); # extend to des & species?
	my $sthUpProt1=$dbh->prepare("UPDATE PROTEIN P INNER JOIN MASTER_PROTEIN M ON P.ID_MASTER_PROTEIN=M.ID_MASTER_PROTEIN SET P.PROT_SEQ=M.PROT_SEQ,P.ID_MASTER_PROTEIN=NULL WHERE P.ID_PROTEIN=?"); # reinject prot_seq (in case unable to re-map)
	my $sthUpProt2=$dbh->prepare("UPDATE PROTEIN SET ID_MASTER_PROTEIN=NULL WHERE ID_PROTEIN=?"); # extend to des & species?
	my $protCount=0;
	$sthMapProt->execute;
	while (my ($protID,$masterProtID,$protSeq)=$sthMapProt->fetchrow_array) {
		$forceRemovedProteins{$protID}=$masterProtID; # keep track in case unable to re-map
		$forceUnusedMasters{$masterProtID}=1; # keep track in case unable to re-map
		if ($protSeq eq '+') {$sthUpProt1->execute($protID);}
		else {$sthUpProt2->execute($protID);}
		$protCount++;
		print '.' unless $protCount % 1000;
	}
	$sthMapProt->finish;
	$sthUpProt1->finish;
	$sthUpProt2->finish;
	my ($totUnmapped)=$dbh->selectrow_array("SELECT COUNT(*) FROM PROTEIN WHERE ID_PROJECT=$projectID AND ID_MASTER_PROTEIN IS NULL");

	###>Retrieve uniprotID & backup uniprot ACC from master proteins<###
	my $sthForceID=$dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN=? AND ID_IDENTIFIER=$identifierTypes{ID}"); # force update only
	my $sthForceAC=$dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN=? AND ID_IDENTIFIER=$identifierTypes{AC} AND RANK=1"); # force update only
	foreach my $masterProtID (keys %forceUnusedMasters) {
		$sthForceID->execute($masterProtID);
		my ($uniID)=$sthForceID->fetchrow_array;
		$forceUpdateProteins{$uniID}=$masterProtID;
		$sthForceAC->execute($masterProtID);
		my ($uniAC)=$sthForceAC->fetchrow_array; # backup ACC in case ID is no longer valid
		$backupUniprotAC{$uniAC}=$uniID;
	}
	$sthForceID->finish;
	$sthForceAC->finish;
	print " Done.\n-",(scalar keys %forceUpdateProteins)," proteins (matching ",(scalar keys %forceUnusedMasters)," UniProt IDs) were unmapped\n";
	print "-$totUnmapped proteins without mapping in project to be updated.\n";

	###>Retrieve list of analyses in project <###
	my $sthSelAnalysisByProject=$dbh->prepare("SELECT DISTINCT(AP.ID_ANALYSIS) FROM PROTEIN P, ANALYSIS_PROTEIN AP WHERE P.ID_PROJECT=$projectID AND P.ID_PROTEIN=AP.ID_PROTEIN");
	$sthSelAnalysisByProject->execute;
	while (my ($anaID)=$sthSelAnalysisByProject->fetchrow_array) {
		push @analysisList,$anaID;
	}
	$sthSelAnalysisByProject->finish;
}
else {
	@analysisList=&promsMod::cleanNumericalParameters(split(',',$itemIdString));
}

###########################################
####>Fetching identifiers to be mapped<####
###########################################
my ($projectIdentMapID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM PROJECT WHERE ID_PROJECT=$projectID AND ID_IDENTIFIER != (SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GI')"); # skip GI_ACCESSION mapping!
#my $sthDB=$dbh->prepare("SELECT IDENTIFIER_TYPE FROM ANALYSIS A,DATABANK D WHERE A.ID_ANALYSIS=? AND A.ID_DATABANK=D.ID_DATABANK");
my $sthDB=$dbh->prepare("SELECT IDENTIFIER_TYPE,DB_RANK FROM ANALYSIS_DATABANK AD,DATABANK D WHERE AD.ID_ANALYSIS=? AND AD.ID_DATABANK=D.ID_DATABANK");
#my $sthID=$dbh->prepare("SELECT P.ID_PROTEIN,IDENTIFIER,DB_RANK,PROT_LENGTH,ORGANISM FROM ANALYSIS_PROTEIN A,PROTEIN P WHERE ID_ANALYSIS=? AND P.ID_PROTEIN=A.ID_PROTEIN AND ID_MASTER_PROTEIN IS NULL");
#my $strgSQL = ($update eq 'force')? "" : "AND (ID_MASTER_PROTEIN IS NULL OR PROT_SEQ='-')";
my $sthID=$dbh->prepare("SELECT P.ID_PROTEIN,IDENTIFIER,DB_RANK,PROT_LENGTH,ORGANISM,PROT_SEQ,ID_MASTER_PROTEIN FROM ANALYSIS_PROTEIN A,PROTEIN P WHERE ID_ANALYSIS=? AND P.ID_PROTEIN=A.ID_PROTEIN AND (ID_MASTER_PROTEIN IS NULL OR PROT_SEQ='-')"); # $strgSQL
my $sthMP=$dbh->prepare("SELECT ID_MASTER_PROTEIN FROM PROTEIN WHERE IDENTIFIER=? AND ID_MASTER_PROTEIN IS NOT NULL LIMIT 0,1"); # get 1st
my $sthUpMP=$dbh->prepare("UPDATE PROTEIN SET ID_MASTER_PROTEIN=?,ALIAS=? WHERE ID_PROTEIN=?");
my $sthOrg=$dbh->prepare("SELECT P.ID_PROTEIN FROM PROTEIN P,ANALYSIS_PROTEIN A WHERE ID_ANALYSIS=? AND ID_MASTER_PROTEIN IS NOT NULL AND P.ID_PROTEIN=A.ID_PROTEIN AND ORGANISM='unknown organism'");
my $sthAC=($projectIdentMapID)? $dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER WHERE ID_IDENTIFIER=$projectIdentMapID AND ID_MASTER_PROTEIN=? AND RANK=1 LIMIT 0,1") : undef;

my (%proteinList,%identifier2Master,%missingProteins,%noOrganismProteins,%directProteinAnnotation,%master2AliasConv);
my $numDirectMapped=0;
foreach my $anaID (@analysisList) {
	my %dbRankIdentType;
	if ($dbType eq 'AUTO') {
		$sthDB->execute($anaID);
		while (my ($identType,$dbRank)=$sthDB->fetchrow_array) {
			$identType='UNKNOWN' unless $identType;
			$dbRankIdentType{$dbRank}=$identType;
#print $identType."#".$dbRank."\n<br>";
		}
	}
	$sthID->execute($anaID);
	while (my ($protID,$identifier,$dbRank,$protLength,$organism,$protSeq,$masterProtID0)=$sthID->fetchrow_array) {
		my $identType=($dbType eq 'AUTO')? $dbRankIdentType{$dbRank} : $dbType;
#print "$protID,$identifier,$dbRank,$protLength,$organism,$protSeq,master=$masterProtID0\n<br>";
		if ($masterProtID0 && $update ne 'force') { # already mapped
			$identifier2Master{$identifier}=$masterProtID0;
		}
		else {
			#>Check if this identifier is not already mapped for previous analysis or another project
			if ($update ne 'force'){
				unless ($identifier2Master{$identifier}) {
					$sthMP->execute($identifier);
					my ($masterProtID)=$sthMP->fetchrow_array;
					if ($masterProtID) {
						$identifier2Master{$identifier}=$masterProtID;
						if ($projectIdentMapID && !$master2AliasConv{$masterProtID}) {
							$sthAC->execute($masterProtID);
							($master2AliasConv{$masterProtID})=$sthAC->fetchrow_array;
						}
					}
				}
				if ($identifier2Master{$identifier}) {
					my $alias=($master2AliasConv{$identifier2Master{$identifier}})? $master2AliasConv{$identifier2Master{$identifier}} : $identifier;
					$sthUpMP->execute($identifier2Master{$identifier},$alias,$protID);
					$numDirectMapped++;
				}
			}
			#else { # to be mapped
			if ($update eq 'force' || !$identifier2Master{$identifier}){
				if ($identType=~/^UPS_(\w+)/) { # UPS Standards
					my $id=$1;
					my ($modIdentType,$modIdentifier);
					if ($id=~/ACCESSION|ALL/) {
						$modIdentType='UNIPROT_ACCESSION';
						($modIdentifier)=($identifier=~/^([\w\.\-]+)ups/);
					}
					else {
						$modIdentType='UNIPROT_ID';
						($modIdentifier)=($identifier=~/(\w+)[\.\-]*_UPS\Z/);
					}
					$proteinList{$modIdentType}{$protID}=$modIdentifier;
				}
				elsif ($identType eq 'UNIPROT_ID') {
					my ($modIdentifier)=($identifier=~/^(\w+)/); # exclude myProMS-specific isoform tag if any (ACC with iso tag are OK)
					$proteinList{$identType}{$protID}=$modIdentifier;
				}
				else {$proteinList{$identType}{$protID}=$identifier;}
#print "protID=$protID##$identType##$identifier\n<br>";
			}
		}
		if (!$protLength || $protSeq eq '-') { #  all annotations missing or just prot seq
			if ($identType eq 'GI_ACCESSION' || $identType eq 'NCBI_ALL') {
				push @{$directProteinAnnotation{'GI_ACCESSION'}{'ANA'}{$protID}},$anaID;
				my ($modIdentifier)=($identifier=~/gi\|(\d+)/);
				if ($modIdentifier) {$directProteinAnnotation{'GI_ACCESSION'}{'ID'}{$modIdentifier}=$protID;}
				else {push @{$missingProteins{$protID}},$anaID;} # record as missing annotation in table PROTEIN
			}
			elsif ($identType=~/UNIPROT_(ACCESSION|ALL)/ && $identifier=~/-/) { # Isoform!!!
				push @{$directProteinAnnotation{'UNIPROT_ACCESSION'}{'ANA'}{$protID}},$anaID;
				#my ($modIdentifier)=($identifier=~/([A-Z]\w{5}-\d+)/); # eg. P31947-2 (?? why not /([A-Z]\w{5}-*\d*)/ ?? PP 03/03/17)
				my ($modIdentifier)=($identifier=~/([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}-\d+)/); # eg. P31947-2
				if ($modIdentifier) {$directProteinAnnotation{'UNIPROT_ACCESSION'}{'ID'}{$modIdentifier}=$protID;}
				else {push @{$missingProteins{$protID}},$anaID;} # record as missing annotation in table PROTEIN
			}
			elsif ($identType=~/UPS_(ACCESSION|ALL)/) { # UPS Standards
				#my $id=$1;
				#my ($modIdentType,$modIdentifier);
				#if ($id=~/ACCESSION|ALL/) {
				#	my $modIdentType='UNIPROT_ACCESSION';
					my ($modIdentifier)=($identifier=~/^([\w\.\-]+)ups/);
				#}
				#else {
				#	$modIdentType='UNIPROT_ID';
				#	($modIdentifier)=($identifier=~/(\w+)_UPS\Z/);
				#}
				push @{$directProteinAnnotation{'UNIPROT_ACCESSION'}{'ANA'}{$protID}},$anaID;
				if ($modIdentifier) {$directProteinAnnotation{'UNIPROT_ACCESSION'}{'ID'}{$modIdentifier}=$protID;}
			}
			else {
				push @{$missingProteins{$protID}},$anaID; # record as missing annotation in table PROTEIN
			}
		}
		elsif (!$organism || $organism=~/unknown/i) {
			$noOrganismProteins{$protID}=1;
		}
	}
	#>More proteins with missing organism (master protein is known) # databank parsing error correction
	$sthOrg->execute($anaID);
	while (my ($protID)=$sthOrg->fetchrow_array) {
		$noOrganismProteins{$protID}=1;
	}
	$dbh->commit;
}

$sthDB->finish;
$sthID->finish;
$sthMP->finish;
$sthUpMP->finish;
$sthOrg->finish;
$sthAC->finish if $projectIdentMapID;

$dbh->disconnect;
#exit;

###########################################################################
####>Direct retrieval of missing annotations from dedicated ressources<####
###########################################################################
&directProteinAnnotation(\%directProteinAnnotation) if scalar keys %directProteinAnnotation;

my $numProt2Map=0;
foreach my $identType (keys %proteinList) {
	$numProt2Map+=scalar keys %{$proteinList{$identType}};
}
print "+$numProt2Map/",($numProt2Map+$numDirectMapped)," proteins to annotate ($numDirectMapped proteins directly annotated).\n";


############################################
####>Mapping identifiers to uniprot IDs<####
############################################
my %uniprotIdList;
my $numProtSkipped=0;
foreach my $identType (keys %proteinList) {
	my $okIdentType=0;
	if ($identType eq 'UNIPROT_ID') {
		$okIdentType=&processUniprotIdentifier($proteinList{$identType},\%uniprotIdList);
#print "ident=$okIdentType ## size prot = ",scalar keys %uniprotIdList,"\n<br>";
	}
	elsif ($identType eq 'UNIPROT_ACCESSION') {
		$okIdentType=&processUniprotAccession($proteinList{$identType},\%uniprotIdList);
	}
	elsif ($identType eq 'UNIPROT_ALL') {
		$okIdentType=&processUniprotAccession($proteinList{$identType},\%uniprotIdList);
	}
	elsif ($identType eq 'IPI_ACCESSION') {
		my %identifier2ProtID;
		foreach my $protID (keys %{$proteinList{$identType}}) {
			my ($modIdentifier)=($proteinList{$identType}{$protID}=~/(IPI\d+)/); # ignore version number IPIXXXX.n -> IPIXXXX
			next unless $modIdentifier;
			push @{$identifier2ProtID{$modIdentifier}},$protID;
			$okIdentType=1;
		}
		&uniprotMapping(\%identifier2ProtID,'P_IPI',\%uniprotIdList) if $okIdentType;
	}
	elsif ($identType eq 'GI_ACCESSION' || $identType eq 'NCBI_ALL') {
		$okIdentType=&convertGI($proteinList{$identType},\%uniprotIdList);
	}
	elsif ($identType eq 'FLYBASE_ID') {
		$okIdentType=&convertFlybase($proteinList{$identType},\%uniprotIdList);
	}
	unless ($okIdentType) { # Unknown or bad identifier type => guess
		my $firstIdentifier=(values %{$proteinList{$identType}})[0];
		if ($firstIdentifier=~/[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}/) { # [A-Z]\d\w{3}\d
			$okIdentType=&processUniprotAccession($proteinList{$identType},\%uniprotIdList);
		}
		if (!$okIdentType && $firstIdentifier=~/\w+_\w+/) {
			$okIdentType=&processUniprotIdentifier($proteinList{$identType},\%uniprotIdList);
		}
		if (!$okIdentType && $firstIdentifier=~/gi\|\d+/) {
			$okIdentType=&convertGI($proteinList{$identType},\%uniprotIdList);
		}
		if (!$okIdentType && $firstIdentifier=~/^FB(pp|gn|tr)/) {
			$okIdentType=&convertFlybase($proteinList{$identType},\%uniprotIdList);
		}
		if (!$okIdentType) {
			my $numProt=scalar keys %{$proteinList{$identType}};
			$numProtSkipped+=$numProt;
			print "*Warning: Unknown or inappropriate identifier type found: '$identType' for identifiers like '$firstIdentifier' ($numProt proteins will not be mapped).\n";
		}
	}
}
print "-$numProtSkipped proteins excluded from mapping due to unrecognized identifier type.\n" if $numProtSkipped;
$numProt2Map-=$numProtSkipped;

$dbh=&promsConfig::dbConnect('no_user'); # reconnect

unless (scalar keys %uniprotIdList) {
	print "-None of the protein identifiers could be mapped to UniProt IDs.\n";
	#$dbh=&promsConfig::dbConnect('no_user'); # reconnect
	&updateProteinAnnotationFromMaster($dbh,\%missingProteins) if scalar keys %missingProteins;
	&updateProteinOrganism($dbh,\%noOrganismProteins) if scalar keys %noOrganismProteins;
	#&checkForMissingSequences($dbh,\@analysisList);
	$dbh->disconnect;
	my $jobEndTime=strftime("%Y%m%d%H%M%S",localtime);
	print "<JOB END $jobEndTime\n";
	unlink $jobFlagFile;
	exit;
}


##########################
####>Multi-query hash<####
##########################
my %sthQueries=(
	'INS_MP'=>$dbh->prepare('INSERT INTO MASTER_PROTEIN (ID_SPECIES,PROT_DES,PROT_LENGTH,MW,PROT_SEQ,UPDATE_DATE) VALUES (?,?,?,?,?,NOW())'), # autoincrement
	'INS_MI'=>$dbh->prepare('INSERT INTO MASTERPROT_IDENTIFIER (ID_IDENTIFIER,ID_MASTER_PROTEIN,RANK,VALUE) VALUES (?,?,?,?)'),
	'INS_SP'=>$dbh->prepare('INSERT INTO SPECIES (ID_SPECIES,COMMON_NAME,SCIENTIFIC_NAME,TAXONID) VALUES (?,?,?,?)'),
	'UP_PROT'=>$dbh->prepare("UPDATE PROTEIN SET ID_MASTER_PROTEIN=? WHERE ID_PROTEIN=? AND ID_MASTER_PROTEIN IS NULL"), # AND ID_MASTER_PROTEIN IS NULL <- just to be safe
	'UP_PSQ'=>$dbh->prepare("UPDATE PROTEIN SET PROT_SEQ='+' WHERE ID_PROTEIN=? AND (PROT_SEQ='-' OR PROT_SEQ IS NULL OR PROT_SEQ=?)") # in case no protSeq or protSeq = masterSeq
);
if ($projectIdentMapID) {
	$sthQueries{'SEL_PIM'}=$dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER WHERE ID_IDENTIFIER=$projectIdentMapID AND ID_MASTER_PROTEIN=? AND RANK=1 LIMIT 0,1");
	$sthQueries{'UP_PROT2'}=$dbh->prepare("UPDATE PROTEIN SET ID_MASTER_PROTEIN=?,ALIAS=? WHERE ID_PROTEIN=? AND ID_MASTER_PROTEIN IS NULL"); # AND ID_MASTER_PROTEIN IS NULL <- just to be safe
}
if ($update eq 'force') {
	$sthQueries{'SEL_PIM2'}=$dbh->prepare("SELECT ID_MASTER_PROTEIN FROM MASTERPROT_IDENTIFIER WHERE ID_IDENTIFIER=$identifierTypes{ID} AND VALUE=? LIMIT 0,1");
	$sthQueries{'SEL_PSEQ'}=$dbh->prepare("SELECT ID_PROTEIN,PROT_SEQ FROM PROTEIN WHERE ID_PROJECT=$projectID AND PROT_LENGTH=?");
	$sthQueries{'SEL_MI'}=$dbh->prepare("SELECT I.CODE,GROUP_CONCAT(MI.VALUE ORDER BY MI.RANK SEPARATOR ':') FROM MASTERPROT_IDENTIFIER MI,IDENTIFIER I WHERE MI.ID_MASTER_PROTEIN=? AND MI.ID_IDENTIFIER=I.ID_IDENTIFIER GROUP BY MI.ID_IDENTIFIER");
	$sthQueries{'SEL_MPA'}=$dbh->prepare("SELECT PROT_LENGTH,MW,PROT_SEQ FROM MASTER_PROTEIN WHERE ID_MASTER_PROTEIN=?");
	#$sthQueries{'UP_OLD_SEQ'}=$dbh->prepare("UPDATE PROTEIN P INNER JOIN MASTER_PROTEIN M ON P.ID_MASTER_PROTEIN=M.ID_MASTER_PROTEIN SET P.PROT_SEQ=M.PROT_SEQ WHERE P.ID_MASTER_PROTEIN=? AND P.PROT_SEQ='+'"); # reinject original prot_seq
	$sthQueries{'UP_MP'}=$dbh->prepare('UPDATE MASTER_PROTEIN SET ID_SPECIES=?,PROT_DES=?,PROT_LENGTH=?,MW=?,PROT_SEQ=?,UPDATE_DATE=NOW() WHERE ID_MASTER_PROTEIN=?');
	$sthQueries{'DEL_MI'}=$dbh->prepare("DELETE FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN=?");
	$sthQueries{'UP_MPD'}=$dbh->prepare("UPDATE MASTER_PROTEIN SET UPDATE_DATE=NOW() WHERE ID_MASTER_PROTEIN=?");
}

#####################################################
####>Checking for already referenced UniProt IDs<####
#####################################################
my %directUpdateProteins;
if ($update ne 'force') {
	print "+Checking for already referenced UniProt IDs...";
	##my ($fromSQL,$condSQL) = ($update eq "force")? (",MASTER_PROTEIN MP"," AND MI.ID_MASTER_PROTEIN=MP.ID_MASTER_PROTEIN AND UPDATE_DATE <= DATE_SUB(NOW(),INTERVAL 15 DAY)") : ("","");
	##my $sthUID=$dbh->prepare("SELECT MI.ID_MASTER_PROTEIN FROM MASTERPROT_IDENTIFIER MI $fromSQL WHERE MI.ID_IDENTIFIER=$identifierTypes{ID} AND MI.VALUE=? $condSQL LIMIT 0,1 ");
	my $sthUID=$dbh->prepare("SELECT ID_MASTER_PROTEIN FROM MASTERPROT_IDENTIFIER WHERE ID_IDENTIFIER=$identifierTypes{ID} AND VALUE=? LIMIT 0,1");
	my $count=0;
	my $numMappableProt=0;
	my $numUniId2Annot=0;
	my $numKnownUniIds=0;
	my $numUpdatedProt=0;
	foreach my $uniID (keys %uniprotIdList) {
		$count++;
		print '.' unless $count % 100;

		$sthUID->execute($uniID);
		my ($masterProtID)=$sthUID->fetchrow_array;
		unless ($masterProtID) { # Unknown UniProt ID
			$numUniId2Annot++;
			$numMappableProt+=scalar @{$uniprotIdList{$uniID}};
			#delete $uniprotIdList{$uniID} if ($update eq "force");
			next;
		}
	#@{$directUpdateProteins{$masterProtID}}=@{$uniprotIdList{$uniID}}; # DB update performed later on
	#$numUpdatedProt+=scalar @{$directUpdateProteins{$masterProtID}};
		##if ($update eq 'force' && $dateDiff < 0) { # mapping to master prot is > 1 month old => fetch annotations again
		##	$sthForceAC->execute($masterProtID);
		##	my ($uniAC)=$sthForceAC->fetchrow_array; # backup ACC in case ID is no longer valid
		##	$forceUpdateProteins{$uniID}=$masterProtID;
		##	$backupUniprotAC{$uniAC}=$uniID;
		##	$numForceUpdatedMaster++;
		##}
		##else { # mapping to master prot is up to date => no need to fetch annotations
	@{$directUpdateProteins{$masterProtID}}=@{$uniprotIdList{$uniID}}; # DB update performed later on
	$numUpdatedProt+=scalar @{$directUpdateProteins{$masterProtID}};
			@{$uniprotIdList{$uniID}}=();
			delete $uniprotIdList{$uniID};
		##}
		$numKnownUniIds++;
	}
	$sthUID->finish;

	###>Performing direct protein update<### (if any)
	print '/';
	$count=0;
	foreach my $masterProtID (keys %directUpdateProteins) {
		my $convAlias;
		if ($sthQueries{'SEL_PIM'}) { # project-level identifier conversion
			$sthQueries{'SEL_PIM'}->execute($masterProtID);
			($convAlias)=$sthQueries{'SEL_PIM'}->fetchrow_array;
		}
		foreach my $protID (@{$directUpdateProteins{$masterProtID}}) {
			if ($convAlias) { # alias update
				$sthQueries{'UP_PROT2'}->execute($masterProtID,$convAlias,$protID);
			}
			else { # no alias update
				$sthQueries{'UP_PROT'}->execute($masterProtID,$protID);
			}
			if ($update eq 'force') { # remove from tracking
				delete $forceRemovedProteins{$protID};
			}
		}
		$count++;
		print '.' unless $count % 100;
	}
	$dbh->commit;

	print " Done\n";
	print "-$numUpdatedProt proteins matched $numKnownUniIds UniProt IDs already referenced.\n" if $numUpdatedProt;
	my $numLostProt=$numProt2Map-$numMappableProt-$numUpdatedProt;
	print "-$numLostProt proteins could not be mapped to UniProt IDs.\n" if $numLostProt;
	print "+$numMappableProt proteins to annotate (with $numUniId2Annot UniProt IDs).\n";

}

############################################
####>Fetching annotations via UniProtKB<####
############################################
print "+Retrieving annotations from UniProtKB...";
#my $baseURL='http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-noSession+-ascii+-f+id%20acc%20des%20gen%20org%20txi%20drx%20mow%20seq';
my $baseURL='http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&format=uniprot&style=raw';
my @identifierList=keys %uniprotIdList;
my $numAnnotUniIds=0;
my $numAnnotProt=0;
my $idString;
my $try = 0;
SPLICE: while(my @subIdentList=splice(@identifierList,0,25)) {
	$try++;
	#$idString = join ('|',@subIdentList);
	#my $fullURL="$baseURL+[UNIPROT-id:$idString]";
	$idString = join (',',@subIdentList);
#print "\n>>STRG=$idString\n"; #=====================================================================
	my $fullURL="$baseURL&id=$idString";
	my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
	$agent->timeout(360);
	if (!$promsPath{'http_proxy'}) {$agent->env_proxy;}
	elsif (lc($promsPath{'http_proxy'}) eq 'no') {$agent->no_proxy;}
	else {$agent->proxy(['http','ftp'],$promsPath{'http_proxy'});}
	my $response = $agent->get($fullURL);
	while (my $wait = $response->header('Retry-After')) {
		#print "Waiting ($wait)...\n";
		sleep $wait;
		$response = $agent->get($response->base);
	}
	if ($response->is_success) { #importAnnotationsFromSRS
#print ">>RESP:\n",$response->content,"\n";
		my ($newAnnotUniIds,$newAnnotProt)=&importAnnotationsFromUniprotKB($dbh,\%sthQueries,$response->content,\%uniprotIdList,\%species,\%identifierTypes,\%forceUpdateProteins,\%forceRemovedProteins,\%backupUniprotAC); # commit in subroutine
		$numAnnotUniIds+=$newAnnotUniIds;
		$numAnnotProt+=$newAnnotProt;
		$try = 0;
	}
	else {
		if ($try <= 5) {
			sleep 3; # required by uniprot server
			redo SPLICE;
		}
		else {
			print "***ERROR: Got from UniprotKB: ",$response->status_line," for ",$response->request->uri,"\n"; #chomp()
		}
	}
	print '.';
	sleep 1; # required by server
#last;
}
print " Done.\n-$numAnnotProt proteins ($numAnnotUniIds UniProt IDs) annotated.\n";


###############################################################################################
####>Re-linking unmapped proteins to old their mapping (if any) & deleting unused mappings<####
###############################################################################################
if ($update eq 'force') {
	foreach my $protID (keys %forceRemovedProteins) {
		$sthQueries{'UP_PROT'}->execute($forceRemovedProteins{$protID},$protID);
	}
	print "-",(scalar keys %forceRemovedProteins)," proteins with obselete mapping could not be updated.\n";# if scalar keys %forceRemovedProteins;
	#my $numMasterDeleted=&promsMod::deleteUnusedMasterProteins($dbh); #,\%forceUnusedMasters
	#print "-$numMasterDeleted master proteins deleted.\n"; # if $numMasterDeleted;
}

foreach my $sth (values %sthQueries) {$sth->finish;}


##############################################
####>Updating missing protein annotations<####
##############################################
&updateProteinAnnotationFromMaster($dbh,\%missingProteins) if scalar keys %missingProteins;
&updateProteinOrganism($dbh,\%noOrganismProteins) if scalar keys %noOrganismProteins;
#&checkForMissingSequences($dbh,\@analysisList);
&checkForUnlocalizedPeptides($dbh,\@analysisList);
$dbh->commit;

###############################
####>Cleaning & correction<####
###############################
&correctDuplicateMaster($dbh); # just to be safe (updates proteins but does not delete unused duplicates: Done by &promsMod::deleteUnusedMasterProteins)
$dbh->commit;
my $numMasterDeleted=&promsMod::deleteUnusedMasterProteins($dbh);
print "-$numMasterDeleted master proteins deleted.\n";
$dbh->commit;

$dbh->disconnect;
my $jobEndTime=strftime("%Y%m%d%H%M%S",localtime);
print "<JOB END $jobEndTime\n";
unlink $jobFlagFile;


########################################<<< SUBROUTINES >>>############################################

sub processUniprotIdentifier {
	my ($refProtList,$refUniprotIdList)=@_;
	my %identifier2ProtID;
	my $okIdentType=0;
	my $mapACC=0;
	foreach my $protID (keys %{$refProtList}) {
		my ($identifier)=($refProtList->{$protID}=~/(\w+_\w+)/);
		next unless $identifier; # true UniProt ID is already in list
		if ($identifier=~/^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})_/) { # unstable ID (<AC>_SPECIES) [A-Z]\d\w{3}\d
			push @{$identifier2ProtID{$1}},$protID;
			$mapACC=1;
			$okIdentType=1;
		}
		else { # stable ID
			push @{$refUniprotIdList->{$identifier}},$protID;
			$okIdentType=1;
		}
	}
	&uniprotMapping(\%identifier2ProtID,'ACC',$refUniprotIdList) if $mapACC;
	return $okIdentType;
}
sub processUniprotAccession {
	my ($refProtList,$refUniprotIdList)=@_;
	my %identifier2ProtID;
	my $okIdentType=0;
	foreach my $protID (keys %{$refProtList}) {
		my ($identifier)=($refProtList->{$protID}=~/([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})/); # skip version & isoform number   [A-Z]\d\w{3}\d
		next unless $identifier;
		$okIdentType=1;
		push @{$identifier2ProtID{$identifier}},$protID;
	}
	&uniprotMapping(\%identifier2ProtID,'ACC',$refUniprotIdList) if $okIdentType;
	return $okIdentType;
}

sub convertGI {
	my ($refGiList,$refUniprotIdList)=@_;
#print '*1 ',scalar keys %{$refGiList}," *\n";
	my (@unmappedGIs,%giNum2ProtID,%mapIdentifiers);
	foreach my $protID (keys %{$refGiList}) {
		my ($giNum)=($refGiList->{$protID}=~/gi\|(\d+)/);
		next unless $giNum;
		push @unmappedGIs,$giNum;
		$giNum2ProtID{$giNum}=$protID; # same project so 1 to 1 match???
		$mapIdentifiers{$giNum}=undef; # undef values will indicate unmapped gi identifiers
	}
	my $numTotalGi=scalar (keys %giNum2ProtID);
	return 0 unless $numTotalGi;

	print "\tConverting $numTotalGi gi identifiers to UniProt IDs...";

	####<Trying direct mapping on Uniprot server>####
	my $giString;
    # Splitting request for each 50 proteins #
    my $try = 0;
	my $hasUnexpResp=0;
    SPLICE: while(my @subGiList = splice(@unmappedGIs, 0, 50)) {
        $try++;
#print "TRY #$try...\n";
		$giString = join (' ',@subGiList);
		##my $params = {
		##	from => 'P_GI',
		##	to => 'ID',
		##	format => 'tab',
		##	query => $giString
		##};
        my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
        $agent->timeout(360);
		if (!$promsPath{'http_proxy'}) {$agent->env_proxy;}
		elsif (lc($promsPath{'http_proxy'}) eq 'no') {$agent->no_proxy;}
		else {$agent->proxy(['http','ftp'],$promsPath{'http_proxy'});}
        ##push @{$agent->requests_redirectable}, 'POST';
        ##my $response = $agent->post("http://www.uniprot.org/mapping/", $params);
		my $response = $agent->get("http://www.uniprot.org/mapping?from=P_GI&to=ID&format=tab&query=$giString");
        while (my $wait = $response->header('Retry-After')) {
			#print "Waiting ($wait)...\n";
			sleep $wait;
			$response = $agent->get($response->base);
        }
        if ($response->is_success) {
#print ">RESPONSE:\n",$response->content,"\n";
            my @responseContent = split("\n",$response->content);
            my $respHeader=shift @responseContent; # skip header
            foreach my $mapping (@responseContent){
                my @identData = split "\t", $mapping;
                if ($#identData != 1) {
					sleep 1; # required by uniprot server
                    if ($try <= 5) {
						sleep 3; # required by uniprot server
                        redo SPLICE;
                    }
					else {
						print "***ERROR: Unexpected response content from UniProt mapping service:\n";
						print "----------------\n----Parameters used:\nfrom=>'P_GI', to=>'ID', format=>'tab', query=>'$giString'\n----Response recieved:\n$respHeader\n",join("\n",@responseContent[0..3]),"\n...\n----------------\n" if $hasUnexpResp==0;
						$hasUnexpResp++;
						sleep 3; # required by uniprot server
						next SPLICE;
                    }
                }
				unless ($mapIdentifiers{$identData[0]}) {# use first match only (in case more than 1)
					@{$mapIdentifiers{$identData[0]}}=('spID',$identData[1]);
#print "gi|$identData[0] => $identData[1]\n";
				}
            }
            $try = 0;
        }
		else {
            if ($try <= 5) {
				sleep 3; # required by uniprot server
                redo SPLICE;
            }
			else {
                print "***ERROR: Got from Uniprot: ",$response->status_line," for ",$response->request->uri,"\n";
            }
        }
		print '.';
		sleep 3; # required by server
    }

	####<Trying mapping on NCBI server for missing Uniprot matches>####
	@unmappedGIs=(); # reset
	foreach my $giNum (keys %mapIdentifiers) {
		push @unmappedGIs,$giNum unless $mapIdentifiers{$giNum};
    }
	my $numUnmapped=scalar @unmappedGIs;
	print " ($numUnmapped gi unmapped after 1st pass)...";

#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id=13278382,56405010,126032329
#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=13278382,56405010,126032329&rettype=docsum&retmode=text

	my $numNewMapped=0;
	$try = 0;
	SPLICE2: while(my @subGiList = splice(@unmappedGIs,0,100)) {
		$try++;
		$giString = join (',',@subGiList);
#print "TRY #$try ($giString)...\n";
		my $params = {
			db=>'protein',
			rettype=>'docsum',
			retmode=>'text',
			id=>$giString
		};
		my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
		$agent->timeout(360);
		if (!$promsPath{'http_proxy'}) {$agent->env_proxy;}
		elsif (lc($promsPath{'http_proxy'}) eq 'no') {$agent->no_proxy;}
		else {$agent->proxy(['http','ftp'],$promsPath{'http_proxy'});}
		push @{$agent->requests_redirectable}, 'POST';
		my $response = $agent->post("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",$params);
		while (my $wait = $response->header('Retry-After')) {
			#print "Waiting ($wait)...\n";
			sleep $wait;
			$response = $agent->get($response->base);
		}
		if ($response->is_success) {
			my @responseContent = split("\n",$response->content);
#			foreach my $line (@responseContent) {
#				if ($line=~/gi\|(\d+)/) {
#					my $giNum=$1;
#					next if $mapIdentifiers{$giNum}; # already seen in a previous line
#					my ($dbCode,$dbIdent);
#					if ($line=~/\|sp\|\w+\|(w+_\w+)/) {
#						$dbIdent=$1;
#						if ($dbIdent=~/([A-Z]\d\w{3}\d)_/) { # unstable ID (<AC>_SPECIES)
#							$dbIdent=$1;
#							$dbCode='spAC';
#						}
#						else {$dbCode='spID';} # stable ID
#					}
#					elsif ($line=~/\|(ref|gb|emb|dbj)\|([\w]+)/) { # version number (XXXXX.n) is ignored
#						$dbCode=$1;
#						$dbIdent=$2;
#					}
#					if ($dbIdent) {
#						@{$mapIdentifiers{$giNum}}=($dbCode,$dbIdent) ;
#						$numNewMapped++;
##print "-$line: => ($dbCode) $dbIdent\n";
#					}
#				}
#			}

			my ($giNum,$dbCode,$dbIdent);
			foreach my $line (@responseContent) {
				if ($line=~/<Id>(\d+)<\/Id>/) {
					$giNum=$1;
					$dbCode=$dbIdent=undef;
				}
				elsif ($line=~/<Item Name="Extra" .+>gi\|\d+\|(.+)<\/Item>/) {
					my $identStrg=$1;
					if ($identStrg=~/^sp\|(\w+)\|(w+_\w+)/) {
						my $uniAC=$1;
						$dbIdent=$2;
						if ($dbIdent=~/[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}_/) { # unstable ID (<AC>_SPECIES)  [A-Z]\d\w{3}\d_
							$dbIdent=$uniAC;
							$dbCode='spAC';
						}
						else {$dbCode='spID';} # stable ID
					}
					elsif ($identStrg=~/^(ref|gb|emb|dbj)\|([\w]+)/) { # version number (XXXXX.n) is ignored
						$dbCode=$1;
						$dbIdent=$2;
					}
				}
				elsif ($line=~/<Item Name="ReplacedBy" .+>(.+)<\/Item>/) {
					$dbIdent=$1;
					$dbIdent=~s/\.\d+//;
					$dbCode=($dbIdent=~/.P_\d/)? 'ref' : 'spAC';
				}
				elsif ($line=~/<\/DocSum>/ && $dbIdent) { # end of entry
					@{$mapIdentifiers{$giNum}}=($dbCode,$dbIdent);
					$numNewMapped++;
#print "$giNum\t$dbCode\t$dbIdent\n";
				}
			}
			$try = 0;
		}
		else {
			if ($try <= 5) {
				sleep 3; # required by uniprot server
				redo SPLICE2;
			}
			else {
				print "***ERROR: Got from NCBI: ",$response->status_line," for ",$response->request->uri,"\n";
			}
		}
		print '.';
		sleep 3; # required by server
	}
	print " Done (",$numUnmapped-$numNewMapped," gi unmapped after 2nd pass).\n";

	#my %srsDbCodes=('sp'=>'UNIPROT-id','gb'=>'EMBLCDS','emb'=>'EMBLCDS','dbj'=>'EMBLCDS','ref'=>'refseqp-ID');
	my %uniprotIdentCodes=('spID'=>'UNIPROT_ID','spAC'=>'ACC','gb'=>'EMBL','emb'=>'EMBL','dbj'=>'EMBL','ref'=>'P_REFSEQ_AC');
	my %newIdentifierList;
	my %unknownBdCodes;
	foreach my $giNum (keys %mapIdentifiers) {
		next unless $mapIdentifiers{$giNum};
		my ($bdCode,$newIdent)=@{$mapIdentifiers{$giNum}};
		unless ($uniprotIdentCodes{$bdCode}) {
			$unknownBdCodes{$bdCode}++;
			next;
		}
		if ($uniprotIdentCodes{$bdCode} eq 'UNIPROT_ID') { # matched at 1st pass
			push @{$refUniprotIdList->{$newIdent}},$giNum2ProtID{$giNum};
		}
		else { # prepare for 2nd pass
			push @{$newIdentifierList{$uniprotIdentCodes{$bdCode}}{$newIdent}},$giNum2ProtID{$giNum}; # push @{...{db}{newIdent}},protID
		}
	}
	if (scalar keys %unknownBdCodes) {
		my $numIdent=0;
		foreach my $bdCode (sort keys %unknownBdCodes) {$numIdent+=$unknownBdCodes{$bdCode};}
		print "*Warning: Unknown databank code(s) found in GI some entries: ",join(', ',keys %unknownBdCodes)," ($numIdent GI will not be mapped).\n";
	}

	####<Uniprot mapping>####
	foreach my $identType (keys %newIdentifierList) {
		#next if $identType eq 'UNIPROT_ID'; # no need
		&uniprotMapping($newIdentifierList{$identType},$identType,$refUniprotIdList);
	}
	return 1;
}

sub convertFlybase {
	my ($refProtList,$refUniprotIdList)=@_;

	my (%FBgene2ProtID,%FBidentifier2ProtID,@FB2convert);
	foreach my $protID (keys %{$refProtList}) {
		if ($refProtList->{$protID}=~/^FBgn/) {push @{$FBgene2ProtID{$refProtList->{$protID}}},$protID;}
		elsif ($refProtList->{$protID}=~/^FB/) {
			push @FB2convert,$refProtList->{$protID};
			$FBidentifier2ProtID{$refProtList->{$protID}}=$protID;
		}
	}
	my $num2convert=scalar @FB2convert;
	my $numTotalFB=scalar (keys %FBgene2ProtID) + $num2convert;
	return 0 unless $numTotalFB;

	####<STEP 1: Converting Flybase peptides (or transcripts) to Flybase genes>####
	print "\tConverting $num2convert/$numTotalFB Flybase identifiers to Flybase genes...";
	my $idString;
	# Splitting request for each 50 proteins #
    my $try = 0;
    SPLICE: while(my @subIdList = splice(@FB2convert, 0, 100)) {
        $try++;
#print "TRY #$try...\n";
		$idString = join (' ',@subIdList);
        my $params = {
			mode => 'convert',
			convert => 'fbgn',
			ids => $idString
        };
        my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
        $agent->timeout(360);
        if (!$promsPath{'http_proxy'}) {$agent->env_proxy;}
		elsif (lc($promsPath{'http_proxy'}) eq 'no') {$agent->no_proxy;}
		else {$agent->proxy(['http','ftp'],$promsPath{'http_proxy'});}
        push @{$agent->requests_redirectable}, 'POST';
        my $response = $agent->post("http://flybase.org/cgi-bin/export2batch.pl", $params);
        while (my $wait = $response->header('Retry-After')) {
			#print "Waiting ($wait)...\n";
			sleep $wait;
			$response = $agent->get($response->base);
        }
        if ($response->is_success) {
#print ">RESPONSE:\n",$response->content,"\n";
			(my $htmlCode=$response->content)=~s/<TR>/\n<TR>/; # make sure each table row is on separate line
			foreach my $line (split(/\n/, $htmlCode)) {
				if ($line=~/(FBgn\d{7})/) {
					my $FBgene=$1;
					my ($identifier)=($line=~/<td [^>]+>&nbsp;(FB\w{9})/);
#print ">$identifier => $FBgene\n";
					push @{$FBgene2ProtID{$FBgene}},$FBidentifier2ProtID{$identifier};
				}
			}
			$try = 0;
        }
		else {
            if ($try <= 5) {
				sleep 3; # required by uniprot server
                redo SPLICE;
            }
			else {
                print "***ERROR: Got from Flybase: ",$response->status_line," for ",$response->request->uri,"\n";
            }
        }
		sleep 3; # required by server
		print '.';
	}
	print " Done (",scalar (keys %FBgene2ProtID)," genes found).\n";

	####<STEP 2: Converting Flybase genes to uniprot IDs>####
	&uniprotMapping(\%FBgene2ProtID,'FLYBASE_ID',\%uniprotIdList) if scalar keys %FBgene2ProtID;
}

sub uniprotMapping {
	my ($refIdentList,$identType,$refUniprotIdList)=@_;
#print "IDENT_TYPE: $identType\n";
	my @identifierList=keys %{$refIdentList};
	my $numIdent2Map=scalar @identifierList;
	print "\tMapping $numIdent2Map $identType identifiers to UniProt IDs...";
#print "\n$identType: @identifierList\n";

	my %mapIdentifiers;
	my $identString;
    # Splitting request for each 50 proteins #
	my $SET_SIZE=50;
    my $try = 0;
	my $hasUnexpResp=0;
    SPLICE: while(my @subIdentList = splice(@identifierList,0,$SET_SIZE)) {
        $try++;
#print "TRY #$try...\n";
		$identString = join (' ',@subIdentList);
#print "$identString\n";
		##my $params = {
		##	from => $identType,
		##	to => 'ID',
		##	format => 'tab',
		##	query => $identString
		##};
        my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
        $agent->timeout(360);
        if (!$promsPath{'http_proxy'}) {$agent->env_proxy;}
		elsif (lc($promsPath{'http_proxy'}) eq 'no') {$agent->no_proxy;}
		else {$agent->proxy(['http','ftp'],$promsPath{'http_proxy'});}
        ##push @{$agent->requests_redirectable}, 'POST';
        ##my $response = $agent->post("http://www.uniprot.org/mapping/", $params);
		my $response = $agent->get("http://www.uniprot.org/mapping?from=$identType&to=ID&format=tab&query=$identString");
        while (my $wait = $response->header('Retry-After')) {
			#print "Waiting ($wait)...\n";
			sleep $wait;
			$response = $agent->get($response->base);
        }
        if ($response->is_success) {
#print ">RESPONSE:\n",$response->content,"\n";
            my @responseContent = split("\n",$response->content);
            my $respHeader=shift @responseContent; # skip header
            foreach my $mapping (@responseContent) {
                my @identData = split "\t", $mapping;
                if ($#identData != 1) {
					sleep 1; # required by uniprot server
                    if ($try <= 5) {
						sleep 3; # required by server
                        redo SPLICE;
                    }
					else {
						print "***ERROR: Unexpected response content from UniProt mapping service:\n";
						print "----------------\n----Parameters used:\nfrom=>'$identType', to=>'ID', format=>'tab', query=>'$identString'\n----Response recieved:\n$respHeader\n",join("\n",@responseContent[0..3]),"\n...\n----------------\n" if $hasUnexpResp==0;
						$hasUnexpResp++;
						sleep 3; # required by server
						next SPLICE;
                    }
                }
				unless ($mapIdentifiers{$identData[0]}) {# use first match only (in case more than 1)
					$mapIdentifiers{$identData[0]}=$identData[1];
#print "$identData[0] => $identData[1]\n";
				}
            }
            $try = 0;
        }
		else {
            if ($try <= 5) {
				sleep 3; # required by uniprot server
                redo SPLICE;
            }
			else {
                print "***ERROR: Got from UniProt: ",$response->status_line," for ",$response->request->uri,"\n";
            }
        }
		sleep 3; # required by server
		print '.';
    }
	my $numMapped=0;
	foreach my $oldIdent (keys %mapIdentifiers) {
		push @{$refUniprotIdList->{$mapIdentifiers{$oldIdent}}},@{$refIdentList->{$oldIdent}}; # push @{...{db}{newIdent}},protID list
		$numMapped++;
	}
	print " Done ($numMapped identifiers mapped).\n";
}


sub importAnnotationsFromUniprotKB { #importAnnotationsFromSRS
	my ($dbh,$refQueries,$response,$refUniprotIdList,$refSpecies,$refIdentifierTypes,$refForceUpdateProteins,$refForceRemovedProteins,$refBackupUniprotAC)=@_;

	###<Parsing response text>###
	my (%identifierAnnot,%identifierCrossRef,%speciesMatched);

	my ($currentIdent,$description,$spString,$taxonID,$sequenceSection, $geneSection);
	$response=~s/[\.;]\s+(AC|DE|GN|OS|OX|DR|SQ)(\s)/;\n$1$2/g;
#print "==============================\n$response\n==============================\n";
	foreach my $line (split(/\n/,$response)) {
		#print "line1=$line\n<br>";
		if ($line=~/<HTML>/) { # No match at all or unexpected problem with query
			last;
			return (0,0);
		}
		if ($line=~/^\/\//) {
			#<Species
			unless ($speciesMatched{$taxonID}) {
				#my ($spName,$spScName);
				#($spName)=($spString=~/\((.+)\)\Z/);
				#if ($spName) { # not always there
				#	while ($spName=~/^[^\(]+\)/g) {
				#		$spName=~s/^[^\(]+\(//;
				#	}
				#	my $qspName=quotemeta($spName);
				#	($spScName)=($spString=~/(.+)\s\($qspName\)\Z/);
				#}
				#else {
				#	$spName=$spScName=$spString;
				#}
				my ($spScName,$spName)=($spString=~/^([^(]+) \(([^)]+)/);
				unless ($spName) { # not always there
					$spName=$spScName=$spString;
				}
				@{$speciesMatched{$taxonID}}=($spName,$spScName);
			}
			$spString=undef;
			if ($description) {
				$description=~s/;\Z//;
				$identifierAnnot{$currentIdent}{'DE'}=$description;
			}
			$sequenceSection=undef;
			next;
		}
		elsif ($sequenceSection) {
			$line=~s/\s+//g;
			$identifierAnnot{$currentIdent}{'SQ'}.=$line;
			next;
		}
		next unless ($line=~/^\s*([A-Z]{2})\s+/);
		my $tag=$1;
		if ($tag eq 'ID') {
			($currentIdent)=($line=~/ID\s+(\w+)/);
#print ">$currentIdent:\n";
			%{$identifierAnnot{$currentIdent}}=();
			%{$identifierCrossRef{$currentIdent}}=();
			$description=$spString=$taxonID=undef;
			next;
		}
		elsif ($tag eq 'AC') {
			$line=~/^\s*AC\s+(.+)\s*/;
			push @{$identifierCrossRef{$currentIdent}{'AC'}}, split(';\s*',$1);
#print "AC=@{$identifierCrossRef{$currentIdent}{AC}}\n";
			next;
		}
		elsif ($tag eq 'DE') {
			$description.=' ' if $description;
			$line=~/^\s*DE\s+(.+)\s*/;
			$description.=$1;
			next;
		}
		elsif ($tag eq 'GN') {
#print "--$line\n";
			#my @tempGenes;
			#if ($line=~/\sName=([^;]+)/) {
			#	push @tempGenes,$1;
			#}
			#if ($line=~/\sSynonyms=([^;]+)/) {
			#	push @tempGenes,split(', ',$1);
			#}
			#if ($line=~/\sOrderedLocusNames=([^;]+)/) {
			#	push @tempGenes,split(', ',$1);
			#}
			#if ($line=~/\sORFNames=([^;]+)/) { # ???
			#	push @tempGenes,split(', ',$1);
			#}
			#foreach my $gene (@tempGenes) {
			#	$gene=~s/ \{.+//; # 'Git1 {ECO:0000313|Ensembl:ENSMUSP00000098375,...' -> 'Git1'
			#	$gene=~s/,\Z//;
			#	push @{$identifierCrossRef{$currentIdent}{'GN'}},$gene;
			#}
#print "Genes=@{$identifierCrossRef{$currentIdent}{GN}}\n";

			#####Change parsing due to multi line in gene name
			#####Now concat of line and parsing in OS part
			chomp($line);
			$line=~s/GN//;
			$geneSection.=$line;
			next;
		}
		elsif ($tag eq 'OS') { # multi-line possible
			##parsing gene string
			if ($geneSection) {
				my @tempGenes=();
#print "gene1=$geneSection => $currentIdent<br>\n";
				$geneSection=~s/ \{[^\{]+\}//g; #CG15877 {ECO:0000313|EMBL:AAF47641.1,   ECO:0000313|FlyBase:FBgn0035337},   Dmel_CG15877 {ECO:0000313|EMBL:AAF47641.1}; -> CG15877 Dmel_CG15877
#print "gene2=$geneSection<br>\n";
				foreach my $gnName (split(/\;/, $geneSection)) {
				    if ($gnName=~/\sName=([^;]+)/) {
						push @tempGenes,$1;
				    }
				    if ($gnName=~/\sSynonyms=([^;]+)/) {
						push @tempGenes,split(', ',$1);
				    }
				    if ($gnName=~/\sOrderedLocusNames=([^;]+)/) {
						push @tempGenes,split(', ',$1);
				    }
				    if ($gnName=~/\sORFNames=([^;]+)/) { # ???
						push @tempGenes,split(', ',$1);
				    }
				}
				$geneSection="";
				foreach my $gene (@tempGenes) {
					#$gene=~s/ \{[^\{]+\}//g; #CG15877 {ECO:0000313|EMBL:AAF47641.1,   ECO:0000313|FlyBase:FBgn0035337},   Dmel_CG15877 {ECO:0000313|EMBL:AAF47641.1}; -> CG15877 Dmel_CG15877
					#$gene=~s/ \{.+//; # 'Git1 {ECO:0000313|Ensembl:ENSMUSP00000098375,...' -> 'Git1'
					$gene=~s/,\Z//;
					$gene=~s/\s+//;
					#print $gene."\n";
					push @{$identifierCrossRef{$currentIdent}{'GN'}},$gene;
				}
#print "GN=@{$identifierCrossRef{$currentIdent}{GN}}\n";
			}

			$line=~s/^OS\s+//;
			$line=~s/;\Z//;
			$spString.=' ' if $spString;
			$spString.=$line;
			#my ($spScName,$spName)=($line=~/OS\s+([^\(]+)\s\(([^\)]+)\)/);
			#($spScName,$spName)=($line=~/OS\s+(.+)\s\((.+)\);/);
			next;
		}
		elsif ($tag eq 'OX') {
			($taxonID)=($line=~/NCBI_TaxID=(\d+)/);
#print "TaxonID=$taxonID\n";
			$identifierAnnot{$currentIdent}{'OX'}=$taxonID;
			next;
		}
		elsif ($tag eq 'SQ') { # 2 SQ lines
			if ($line=~/(\d+) MW;/) {
				($identifierAnnot{$currentIdent}{'MW'})=$1;
				($identifierAnnot{$currentIdent}{'AA'})=($line=~/(\d+) AA;/);
			}
			$sequenceSection=1;
			next;
		}
		elsif ($tag eq 'DR') {
			$line=~s/^\s*DR\s+//;
			if ($line=~/^(IPI|PIR|RefSeq|UniGene|Ensembl|GeneID|KEGG|IntAct|STRING|HPA|MIM|MGI|FlyBase|SGD|WormBase);/) {
#print "   '$line' => ";
				my $dbCode=$1;
				$line=~s/;\Z//; # . at end of lines replace by ;
#print "'$line' -> $dbX:";
				my @data=split(/;\s/,$line);
				#if ($dbCode=~/PIR|RefSeq|Ensembl/) {
				if ($dbCode=~/IntAct|STRING/) {
					@{$identifierCrossRef{$currentIdent}{$dbCode}}=($data[1]);
				}
				elsif ($dbCode=~/WormBase/) { # example : 'DR   WormBase; Y48G9A.3; CE48240; WBGene00021697; gcn-1.' > only the three first ids can be linked to the database
					my @dbCodes=('W_CDS','W_P','W_G');
					for (my $i = 1 ; $i < 4 ; $i++){
						@{$identifierCrossRef{$currentIdent}{$dbCodes[$i-1]}}=($data[$i]) if $dbCodes[$i-1] && $data[$i];
					}
				}
				else {
					foreach my $identX (@data[1..$#data]) {
						next if ($dbCode eq "FlyBase" && $identX !~ /FBgn/);
						next if $identX eq '-';
#print "*$identX";
						push @{$identifierCrossRef{$currentIdent}{$dbCode}},$identX;
					}
				}
				#}
				#else {
				#	@{$identifierCrossRef{$currentIdent}{$dbCode}}=($data[1]);
				#}
#print "\n";
			}
		}
	}
#print ">SPECIES:\n";
#foreach my $taxonID (keys %speciesMatched) {
#	print "$taxonID: @{$speciesMatched{$taxonID}}\n";
#}
#print ">MAPPING:\n";
#foreach my $uniprotID (keys %identifierAnnot) {
#	print "+$uniprotID\n";
#	print " OX: $identifierAnnot{$uniprotID}{'OX'}\n DE: $identifierAnnot{$uniprotID}{'DE'}\n AA: $identifierAnnot{$uniprotID}{'AA'}\n MW: $identifierAnnot{$uniprotID}{'MW'}\n SEQ: $identifierAnnot{$uniprotID}{'SQ'}\n";
#	foreach my $dbCode (keys %{$identifierCrossRef{$uniprotID}}) {
#		print " $dbCode: @{$identifierCrossRef{$uniprotID}{$dbCode}}\n";
#	}
#}

	###<Importing new species in DB>###
	my ($maxSpeciesID)=$dbh->selectrow_array('SELECT MAX(ID_SPECIES) FROM SPECIES');
	foreach my $taxonID (keys %speciesMatched) {
		next if $refSpecies->{$taxonID}; # already stored
		$refQueries->{'INS_SP'}->execute(++$maxSpeciesID,@{$speciesMatched{$taxonID}},$taxonID);
		$refSpecies->{$taxonID}=$maxSpeciesID;
	}
	#$dbh->commit;

	###<Importing mapping in DB>###
	my $newAnnotProt=0;
	my %masterProteins;
	my %matchedUniprotIDs;
	foreach my $uniprotID (keys %identifierAnnot) {
#print "\n>>$uniprotID:\n"; # if $uniprotID eq 'CLU_HUMAN';

		###<Force update
		my (%oldIdentifierXref, %oldIdentifierAnnot);
		my ($forceUpMasterID,$isDifferentXref,$isDifferentAnnot,$obsUniprotID,$dbMasterID,$usedMasterID) = (0,0,0,0,0,0); # default
		my $usedUniprotID=$uniprotID; # default
		my $directUpdateUniID=0; # force only
		my $masterProtID;
		if ($update eq 'force') {
			if ($refForceUpdateProteins->{$uniprotID}) {$forceUpMasterID=$refForceUpdateProteins->{$uniprotID};} # uniprotID in DB is still valid
			else { # uniprotID in DB is no longer valid => try to match with backup uniprotACC
				foreach my $uniAC (@{$identifierCrossRef{$uniprotID}{'AC'}}) {
					if ($refBackupUniprotAC->{$uniAC}) { # ACC match!!!
						$obsUniprotID=$refBackupUniprotAC->{$uniAC}; # set uniprotID in DB as obsolete
						$forceUpMasterID=$refForceUpdateProteins->{$obsUniprotID}; # keep same ID in DB but Annot & Xref will be replaced (unless new uniprotID already in DB)
						##<Check if ACC-based mapped uniprotID is already in DB
						$refQueries->{'SEL_PIM2'}->execute($uniprotID);
						($dbMasterID)=$refQueries->{'SEL_PIM2'}->fetchrow_array;
						unless ($dbMasterID) { # new uniprotID is not yet in DB => set forceUpMasterID to be updated
							$dbMasterID=0;
							$isDifferentAnnot=1;
							$isDifferentXref=1;
						} # annot & Xref will be checked later if dbMasterID
						last; # ACC
					}
				}
			}
			unless ($forceUpMasterID) {
				if ($refUniprotIdList->{$uniprotID}) { # Try direct uniID match
#print " Direct Match!\n"; # if $uniprotID eq 'CLU_HUMAN';
					$directUpdateUniID=1; # Assume there are no matching master prot because $forceUpMasterID would be !0
					$isDifferentAnnot=1;
					$isDifferentXref=1;
				}
				else { # Try protein sequence comparison
#print " Check PROT_SEQ ($identifierAnnot{$uniprotID}{AA} res)\n"; # if $uniprotID eq 'CLU_HUMAN';
					$refQueries->{'SEL_PSEQ'}->execute($identifierAnnot{$uniprotID}{'AA'});
					PSEQ:while (my ($protID,$protSeq)=$refQueries->{'SEL_PSEQ'}->fetchrow_array) {
						foreach my $obsUniID (keys %{$refUniprotIdList}) {
							next if $matchedUniprotIDs{$obsUniID};
							foreach my $protID2 (@{$refUniprotIdList->{$obsUniID}}) {
								if ($protID2==$protID) { #protein with same sequence!!!
									$obsUniprotID=$obsUniID;
									$forceUpMasterID=$refForceUpdateProteins->{$obsUniprotID} || 0; # 0 for obsolete uniID with no previous mapping in DB
									##<Check if prot seq-based mapped uniprotID is already in DB
									$refQueries->{'SEL_PIM2'}->execute($uniprotID);
									($dbMasterID)=$refQueries->{'SEL_PIM2'}->fetchrow_array;
									unless ($dbMasterID) { # new uniprotID is not yet in DB => set forceUpMasterID to be updated
										$dbMasterID=0;
										$isDifferentAnnot=1;
										$isDifferentXref=1;
									} # annot & Xref will be checked later if dbMasterID
#print " MATCH PROT=$protID OBS=$obsUniprotID DB=$dbMasterID\n"; # if $uniprotID eq 'CLU_HUMAN';
									last PSEQ;
								}
							}
						}
					}
				}
			}

			next if (!$forceUpMasterID && !$obsUniprotID && !$directUpdateUniID); # could not match this uniprotID to anything => skip (unlikely event)

			$usedUniprotID=$obsUniprotID || $uniprotID;
			next unless $usedUniprotID; # just to safe

			#<Update tracking
			$matchedUniprotIDs{$usedUniprotID}=1;
			foreach my $protID (@{$refUniprotIdList->{$usedUniprotID}}) {
				delete $refForceRemovedProteins->{$protID}; # remove from tracking
			}

			$usedMasterID=$dbMasterID || $forceUpMasterID; # still =0 if $directUpdateUniID

			##<Check Xref
			unless ($isDifferentXref) {
				$refQueries->{'SEL_MI'}->execute($usedMasterID); # $forceUpMasterID
				while (my ($code, $strgValue)=$refQueries->{'SEL_MI'}->fetchrow_array) {
					$oldIdentifierXref{$code}=$strgValue if $code ne 'ID';
				}
				if (scalar(keys %oldIdentifierXref) == scalar (keys %{$identifierCrossRef{$uniprotID}})) { # same number of Xref => check values
					foreach my $code (keys %oldIdentifierXref) {
						if ($identifierCrossRef{$uniprotID}{$code}) {
							if ($oldIdentifierXref{$code} eq join(":",@{$identifierCrossRef{$uniprotID}{$code}})) {
								next;
							}
							else {
								$isDifferentXref=1;
								last;
							}
						}
						else {
							$isDifferentXref=1;
							last;
						}
					}
				}
				else {
					$isDifferentXref=1;
				}
			}
			$refQueries->{'DEL_MI'}->execute($usedMasterID) if $isDifferentXref; # Delete previous annotations


			##<Check annot
			unless ($isDifferentAnnot) {
				$refQueries->{'SEL_MPA'}->execute($usedMasterID); # $forceUpMasterID
				($oldIdentifierAnnot{'AA'},$oldIdentifierAnnot{'MW'},$oldIdentifierAnnot{'SQ'})=$refQueries->{'SEL_MPA'}->fetchrow_array;
				foreach my $code (keys %oldIdentifierAnnot) {
					if ($identifierAnnot{$uniprotID}{$code}) {
						if ($oldIdentifierAnnot{$code} eq $identifierAnnot{$uniprotID}{$code}) {
							next;
						}
						else {
							$isDifferentAnnot=1;
							last;
						}
					}
					else { # unlikely!!!
						$isDifferentAnnot=1;
					}
				}
			}
			if ($isDifferentAnnot) {
				$refQueries->{'UP_MP'}->execute($refSpecies->{$identifierAnnot{$uniprotID}{'OX'}},$identifierAnnot{$uniprotID}{'DE'},$identifierAnnot{$uniprotID}{'AA'},$identifierAnnot{$uniprotID}{'MW'},$identifierAnnot{$uniprotID}{'SQ'},$usedMasterID); # $forceUpMasterID
			}
			$refQueries->{'UP_MPD'}->execute($usedMasterID) if $usedMasterID; # set UPDATE_DATE as NOW() no matter the changes
			$masterProtID=$usedMasterID; # $forceUpMasterID;
		} # end of update eq 'force'

		unless ($masterProtID) { # normal case OR update='force' without previous mapping
			##<Create master prot in DB
			$refQueries->{'INS_MP'}->execute($refSpecies->{$identifierAnnot{$uniprotID}{'OX'}},$identifierAnnot{$uniprotID}{'DE'},$identifierAnnot{$uniprotID}{'AA'},$identifierAnnot{$uniprotID}{'MW'},$identifierAnnot{$uniprotID}{'SQ'});
			$masterProtID=$dbh->last_insert_id(undef,undef,'MASTER_PROTEIN','ID_MASTER_PROTEIN');
		}
		$masterProteins{$masterProtID}=1;

		##<Insert (up to date) Xref values
		if ((!$forceUpMasterID && !$dbMasterID) || $isDifferentXref) {
			$refQueries->{'INS_MI'}->execute($refIdentifierTypes->{'ID'},$masterProtID,1,$uniprotID); # UniProtID
			foreach my $dbCode (keys %{$identifierCrossRef{$uniprotID}}) {
				my $rank=0;
				foreach my $identX (@{$identifierCrossRef{$uniprotID}{$dbCode}}) {
					$refQueries->{'INS_MI'}->execute($refIdentifierTypes->{$dbCode},$masterProtID,++$rank,$identX);
				}
			}
		}

		##<(Re)-Link protein(s) to master prot (Mandatory for new & update mappings!!!)
		my $convAlias;
		if ($refQueries->{'SEL_PIM'}) { # project-level identifier conversion
			$refQueries->{'SEL_PIM'}->execute($masterProtID);
			($convAlias)=$refQueries->{'SEL_PIM'}->fetchrow_array;
		}
		foreach my $protID (@{$refUniprotIdList->{$usedUniprotID}}) {
			if ($convAlias) { # alias update
				$refQueries->{'UP_PROT2'}->execute($masterProtID,$convAlias,$protID);
			}
			else { # no alias update
				$refQueries->{'UP_PROT'}->execute($masterProtID,$protID);
			}
			$refQueries->{'UP_PSQ'}->execute($protID,$identifierAnnot{$uniprotID}{'SQ'}); # PROT_SEQ -> '+'
			$newAnnotProt++;
		}
	}
	$dbh->commit;

	#return (scalar keys %identifierAnnot,$newAnnotProt);
	return (scalar keys %masterProteins,$newAnnotProt);
}

####>Direct annotation retrieval from internet resources (might not match Master Protein annotation)<####
sub directProteinAnnotation {
	my ($refDirectProtAnnot)=@_;
	my $numProt2Annot=my $numAnnotProt=0;
	print "+Updating missing annotations...";

	foreach my $identType (keys %{$refDirectProtAnnot}) { # GI or UNIPROT_ACCESSION
		my ($baseURL,$params,$identKey,$joinKey);
		if ($identType eq 'GI_ACCESSION') {
			$baseURL='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi';
			$params = {
				db=>'protein',
				rettype=>'fasta',
				retmode=>'xml',
			};
			$identKey='id';
			$joinKey=',';
		}
		elsif ($identType eq 'UNIPROT_ACCESSION') {
			$baseURL='http://uniprot.org/batch/?';
			$params = {
				format=>'fasta'
			};
			$identKey='query';
			$joinKey=' ';
		}

		####<Fetching annotations>####
		my %dataAnnot;
		my @missingList=keys %{$refDirectProtAnnot->{$identType}{'ID'}};
		my $identString;
		my $try=0;
		SPLICE: while (my @subList = splice(@missingList,0,100)) {
			$try++;
			$identString = join ($joinKey,@subList);
			$params->{$identKey}=$identString; # update identifier list

			my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
			$agent->timeout(360);
			if (!$promsPath{'http_proxy'}) {$agent->env_proxy;}
			elsif (lc($promsPath{'http_proxy'}) eq 'no') {$agent->no_proxy;}
			else {$agent->proxy(['http','ftp'],$promsPath{'http_proxy'});}
			push @{$agent->requests_redirectable}, 'POST';
			my $response = $agent->post($baseURL,$params);
			while (my $wait = $response->header('Retry-After')) {
				#print "Waiting ($wait)...\n";
				sleep $wait;
				$response = $agent->get($response->base);
			}
			if ($response->is_success) {
				if ($identType eq 'GI_ACCESSION') { # GI
					my ($protIDxml,$organism)=('','');
					my $xml = new XML::Simple(forcearray=>['TSeq']);
					my $xmlData = $xml->XMLin($response->content);
					foreach my $protXML (@{$xmlData->{'TSeq'}}){
						$protIDxml = $protXML->{'TSeq_gi'};
						$dataAnnot{$protIDxml}{'des'} = $protXML->{'TSeq_defline'};
						$dataAnnot{$protIDxml}{'seq'} = $protXML->{'TSeq_sequence'};
						$dataAnnot{$protIDxml}{'length'} = $protXML->{'TSeq_length'};
						if ($protXML->{'TSeq_orgname'}) {
							$organism = $protXML->{'TSeq_orgname'};
							$dataAnnot{$protIDxml}{'des'} =~ s/ \[$organism\]//; #removing organism from description
						}
						$dataAnnot{$protIDxml}{'organism'} = $organism;
					}
				}
				elsif ($identType eq 'UNIPROT_ACCESSION') { # Uniprot ACC for Isoforms
					my ($protACC,$seq);
					foreach my $line (split(/\n/,$response->content)) {
						if ($line=~/^>/) {
							if ($protACC) {
								$dataAnnot{$protACC}{'seq'} = $seq;
								$dataAnnot{$protACC}{'length'} = length($seq);
							}
							($protACC)=($line=~/^>sp\|([^\|]+)/);
							($dataAnnot{$protACC}{'organism'})=($line=~/ OS=(.+)/);
							$dataAnnot{$protACC}{'organism'}=~s/ \w{2}=.+//;
							$line=~s/ \w{2}=.+//;
							($dataAnnot{$protACC}{'des'})=($line=~/^\S+ (.+)/);
							$seq='';
						}
						else {$seq.=$line;}
					}
					if ($seq) {
						$dataAnnot{$protACC}{'seq'} = $seq;
						$dataAnnot{$protACC}{'length'} = length($seq);
					}
				}
			}
			else {
				if ($try <= 5) {
					sleep 3; # required by uniprot server
					redo SPLICE;
				}
				else {
					print "***ERROR: Got from $baseURL: ",$response->status_line," for ",$response->request->uri,"\n";
				}
			}
			sleep 3; # required by server
			print '.';
		}

		####<Storing data in DB>####
		my $dbh=&promsConfig::dbConnect('no_user');
		my $sthPP=$dbh->prepare('SELECT PEP_BEG,PEP_END,FLANKING_AA FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_ANALYSIS=? AND ID_PROTEIN=? AND PEP_BEG > 0');
		my $sthUpP=$dbh->prepare("UPDATE PROTEIN SET PROT_DES=?,PROT_SEQ=?,PROT_LENGTH=?,MW=?,ORGANISM=?,UPDATE_DATE=NOW(),UPDATE_USER='$userID' WHERE ID_PROTEIN=?");
		my $sthUpPC=$dbh->prepare('UPDATE ANALYSIS_PROTEIN SET PEP_COVERAGE=? WHERE ID_ANALYSIS=? AND ID_PROTEIN=?');
		my %massAAave=&promsConfig::getMassAAave; # Average mass, needed for protein mass calculation
		my %massATave=&promsConfig::getMassATave; # Average mass, needed for protein mass calculation

		foreach my $mIdent (keys %dataAnnot) {
			my %countAA;
			foreach my $aa (split(//,$dataAnnot{$mIdent}->{'seq'})) {
				$countAA{$aa}++;
			}
			my $protMass=0;
			foreach my $aa (keys %countAA) {
				$protMass+=($massAAave{$aa}*$countAA{$aa}) if $massAAave{$aa}; # some characters are not amino-acids
			}
			$protMass+=($massATave{H}+$massATave{H}+$massATave{O}) if $protMass; # H(Nter) + OH(Cter)
			$protMass=sprintf "%.2f",$protMass; # no need for more precision
			my $protID=$refDirectProtAnnot->{$identType}{'ID'}{$mIdent};
			$sthUpP->execute($dataAnnot{$mIdent}->{'des'},$dataAnnot{$mIdent}->{'seq'},$dataAnnot{$mIdent}->{'length'},$protMass,$dataAnnot{$mIdent}->{'organism'},$protID);

			#<Updating peptide coverage
			&updatePeptideCoverage($sthPP,$sthUpPC,$protID,$dataAnnot{$mIdent}->{'length'},$refDirectProtAnnot->{$identType}{'ANA'}{$protID});
		}
		$sthPP->finish;
		$sthUpP->finish;
		$sthUpPC->finish;
		$dbh->commit;
		$dbh->disconnect;
		$numProt2Annot+=scalar keys %{$refDirectProtAnnot->{$identType}{'ID'}};
		$numAnnotProt+=scalar keys %dataAnnot;
	}
	print " Done ($numAnnotProt/$numProt2Annot proteins updated).\n";
}

####>Updates protein annotation based on MASTER_PROTEIN table (PROTEIN.PROT_SEQ set to '+')<#####
sub updateProteinAnnotationFromMaster {
	my ($dbh,$refMissingProteins)=@_;
	print "+Updating proteins with missing annotation...";
	my $sthMP=$dbh->prepare("SELECT ID_MASTER_PROTEIN FROM PROTEIN WHERE ID_PROTEIN=?");
	my $sthMPA=$dbh->prepare("SELECT PROT_DES,PROT_LENGTH,MW,SCIENTIFIC_NAME FROM MASTER_PROTEIN M,SPECIES S WHERE M.ID_SPECIES=S.ID_SPECIES AND ID_MASTER_PROTEIN=?");
	my $sthPP=$dbh->prepare('SELECT PEP_BEG,PEP_END,FLANKING_AA FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_ANALYSIS=? AND ID_PROTEIN=? AND PEP_BEG > 0');
	my $sthUpP=$dbh->prepare("UPDATE PROTEIN SET PROT_DES=?,PROT_SEQ='+',PROT_LENGTH=?,MW=?,ORGANISM=?,UPDATE_DATE=NOW(),UPDATE_USER='$userID' WHERE ID_PROTEIN=?");
	my $sthUpPC=$dbh->prepare('UPDATE ANALYSIS_PROTEIN SET PEP_COVERAGE=? WHERE ID_ANALYSIS=? AND ID_PROTEIN=?');

	my %master2MissingProteins;
	foreach my $protID (keys %{$refMissingProteins}) {
		$sthMP->execute($protID);
		my ($masterProtID)=$sthMP->fetchrow_array;
		push @{$master2MissingProteins{$masterProtID}},$protID if $masterProtID;
	}
	$sthMP->finish;
	my $numProtAnnot=0;
	foreach my $masterProtID (keys %master2MissingProteins) {
		$sthMPA->execute($masterProtID);
		my ($des,$length,$mw,$organism)=$sthMPA->fetchrow_array;

		foreach my $protID (@{$master2MissingProteins{$masterProtID}}) {
			$sthUpP->execute($des,$length,$mw,$organism,$protID);

			##<Computing peptide coverage
			&updatePeptideCoverage($sthPP,$sthUpPC,$protID,$length,$refMissingProteins->{$protID});
			$numProtAnnot++;
		}
	}
	$sthMPA->finish;
	$sthPP->finish;
	$sthUpP->finish;
	$sthUpPC->finish;
	$dbh->commit;
	print " Done ($numProtAnnot/",scalar keys %{$refMissingProteins}," proteins updated).\n";
}

sub updatePeptideCoverage {
	my ($sthPP,$sthUpPC,$protID,$length,$refAnaList)=@_;
#print "\n>$protID:\n";
	##<Computing peptide coverage
	foreach my $anaID (@{$refAnaList}) {
		my $covLength=$length;
		my %boundaryStatus;
		my $endMatched=0;
		$sthPP->execute($anaID,$protID);
		while (my ($pepBeg,$pepEnd,$flkAA)=$sthPP->fetchrow_array) {
			$boundaryStatus{$pepBeg}++;
			$boundaryStatus{$pepEnd}--;
			$covLength=$pepEnd if $pepEnd > $covLength; # length might not be right due to protein isoforms
			$endMatched=1 if (!$endMatched && $flkAA && $flkAA=~/.-/); # peptide matches end of sequence ($flkAA can be undef for very old analyses)
		}
		my $coverage=0;
		my $hasPeptide=0;
		my $boundaryNter=0;
		foreach my $boundary (sort{$a<=>$b} keys %boundaryStatus) {
			if ($hasPeptide==0) { # start of peptide region (cannot become <0!)
				$boundaryNter=$boundary;
			}
			$hasPeptide+=$boundaryStatus{$boundary};
			if ($hasPeptide==0) { # end of peptide region (should be true for last boundary too)
				$coverage+=($boundary-$boundaryNter)+1;
			}
		}
		my $pepCoverage=sprintf "%.1f",($coverage*100)/$covLength;
		$pepCoverage*=(!$endMatched && $covLength > $length)? -1 : 1; # 25.0 -> 25   -1: flag for protLength problem
		$sthUpPC->execute($pepCoverage,$anaID,$protID) if scalar keys %boundaryStatus; # to prevent  0% in case no peptides with beg > 0
#print "$anaID: $pepCoverage\n";
	}
}

sub updateProteinOrganism {
	my ($dbh,$refNoOrganismProteins)=@_;
	print "+Updating proteins with missing organism...";
	my $sthOS=$dbh->prepare("SELECT SCIENTIFIC_NAME FROM PROTEIN P,MASTER_PROTEIN M,SPECIES S WHERE ID_PROTEIN=? AND P.ID_MASTER_PROTEIN=M.ID_MASTER_PROTEIN AND M.ID_SPECIES=S.ID_SPECIES");
	my $sthUpOrg=$dbh->prepare("UPDATE PROTEIN SET ORGANISM=? WHERE ID_PROTEIN=?");
	my $numProtOrg=0;
	foreach my $protID (keys %{$refNoOrganismProteins}) {
		$sthOS->execute($protID);
		my ($organism)=$sthOS->fetchrow_array;
		next unless $organism;
		$sthUpOrg->execute($organism,$protID);
		$numProtOrg++;
	}
	$sthOS->finish;
	$sthUpOrg->finish;
	$dbh->commit;
	print " Done ($numProtOrg/",scalar keys %{$refNoOrganismProteins}," proteins updated).\n";
}


#sub checkForMissingSequences { # Fixing missing protein sequence despite master protein. All other fields are ok due to (re)mapping
#	my ($dbh,$refAnaList)=@_;
#	print "+Checking for missing protein sequences...";
#	my $sthSQ=$dbh->prepare("SELECT P.ID_PROTEIN FROM ANALYSIS_PROTEIN A,PROTEIN P WHERE A.ID_PROTEIN=P.ID_PROTEIN AND ID_MASTER_PROTEIN IS NOT NULL AND PROT_SEQ='-' AND ID_ANALYSIS=?");
#	my $sthUpSQ=$dbh->prepare("UPDATE PROTEIN SET PROT_SEQ='+' WHERE ID_PROTEIN=?");
#	my $missingSeq=0;
#	foreach my $anaID (@{$refAnaList}) {
#		my $missingAna=0;
#		$sthSQ->execute($anaID);
#		while (my ($protID)=$sthSQ->fetchrow_array) {
#			$sthUpSQ->execute($protID);
#			$missingAna++;
#		}
#		$dbh->commit if $missingAna;
#		$missingSeq+=$missingAna;
#	}
#	$sthUpSQ->finish;
#	$sthSQ->finish;
#	print " Done ($missingSeq proteins updated).\n";
#}

sub checkForUnlocalizedPeptides {
	my ($dbh,$refAnaList)=@_;
	print "+Checking for unlocalized peptides...";
	my $sthBad=$dbh->prepare("SELECT P.ID_PROTEIN,PE.ID_PEPTIDE,PEP_SEQ,IS_SPECIFIC,VALID_STATUS FROM PEPTIDE_PROTEIN_ATTRIB PPA,PROTEIN P,PEPTIDE PE WHERE PPA.ID_PROTEIN=P.ID_PROTEIN AND PPA.ID_PEPTIDE=PE.ID_PEPTIDE AND PEP_BEG=0 AND PROT_SEQ != '-' AND PROT_SEQ IS NOT NULL AND PPA.ID_ANALYSIS=? ORDER BY P.ID_PROTEIN,PE.ID_PEPTIDE");
	my $sthPS=$dbh->prepare("SELECT PROT_SEQ,ID_MASTER_PROTEIN FROM PROTEIN WHERE ID_PROTEIN=?");
	my $sthMP=$dbh->prepare("SELECT PROT_SEQ FROM MASTER_PROTEIN WHERE ID_MASTER_PROTEIN=?");
	my $sthDelPPA=$dbh->prepare("DELETE FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_PROTEIN=? AND ID_PEPTIDE=?");
	my $sthUpPPA=$dbh->prepare("INSERT INTO PEPTIDE_PROTEIN_ATTRIB (ID_ANALYSIS,ID_PROTEIN,ID_PEPTIDE,PEP_BEG,PEP_END,FLANKING_AA,IS_SPECIFIC) VALUES (?,?,?,?,?,?,?)");

	my (%peptideIDs,%peptideSpecif,%pepID2anaID,%validStatus,%proteinData);
	foreach my $anaID (@{$refAnaList}) {
		$sthBad->execute($anaID);
		while (my ($protID,$pepID,$pepSeq,$isSpecific,$valStatus)=$sthBad->fetchrow_array) {
			unless ($proteinData{$protID}) {
				$sthPS->execute($protID);
				@{$proteinData{$protID}}=$sthPS->fetchrow_array;
				if ($proteinData{$protID}[1] && $proteinData{$protID}[0] eq '+') {
					$sthMP->execute($proteinData{$protID}[1]);
					($proteinData{$protID}[0])=$sthMP->fetchrow_array;
				}
			}
			push @{$peptideIDs{$protID}{$pepSeq}},$pepID;
			$peptideSpecif{$protID}{$pepSeq}=$isSpecific if (!$peptideSpecif{$protID} || !$peptideSpecif{$protID}{$pepSeq}); # keep best is_specif if any
			$pepID2anaID{$pepID}=$anaID;
			$validStatus{$pepID}=$valStatus || 0;
		}
	}
	$sthBad->finish;
	$sthPS->finish;
	$sthMP->finish;

	my (%numPepUpdated,%numProtMatched,%anaMatchedHash);
	foreach my $protID (keys %peptideIDs) {
		my $protSeq=$proteinData{$protID}[0];
		my %numMatches;
		foreach my $pepSeq (keys %{$peptideIDs{$protID}}) {
			while ($protSeq=~/(\w?)$pepSeq(\w?)/g) { # potential multiple matches
				my ($resBeg,$resEnd)=($1,$2);
				my $startPos=$-[0] + 1;
				my $endPos=$+[0];
				if ($resBeg) {$startPos++;}
				else {$resBeg='-';}
				if ($resEnd) {$endPos--;}
				else {$resEnd='-';}
				@{$numMatches{$pepSeq}{$startPos}}=($resBeg,$resEnd,$endPos);
			}
		}
		foreach my $pepSeq (keys %numMatches) {
			foreach my $pepID (@{$peptideIDs{$protID}{$pepSeq}}) {
				my $factor=($validStatus{$pepID})? 1 : -1;
				$sthDelPPA->execute($protID,$pepID);
				foreach my $startPos (keys %{$numMatches{$pepSeq}}) {
					my ($resBeg,$resEnd,$endPos)=@{$numMatches{$pepSeq}{$startPos}};
					$sthUpPPA->execute($pepID2anaID{$pepID},$protID,$pepID,$factor*$startPos,$factor*$endPos,$resBeg.$resEnd,$peptideSpecif{$pepSeq} || undef);
				}
				$numPepUpdated{$pepID}=1;
				$numProtMatched{$protID}=$validStatus{$pepID} unless $numProtMatched{$protID}; # never goes back to 0 if already 1
				$anaMatchedHash{$protID}{$pepID2anaID{$pepID}}=1;
			}
		}
	}
	$sthDelPPA->finish;
	$sthUpPPA->finish;

	#<Update peptide coverage
	if (scalar keys %numPepUpdated) {
		my $sthPP=$dbh->prepare('SELECT PEP_BEG,PEP_END,FLANKING_AA FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_ANALYSIS=? AND ID_PROTEIN=? AND PEP_BEG > 0');
		my $sthUpPC=$dbh->prepare('UPDATE ANALYSIS_PROTEIN SET PEP_COVERAGE=? WHERE ID_ANALYSIS=? AND ID_PROTEIN=?');
		foreach my $protID (keys %numProtMatched) {
			next if $numProtMatched{$protID}==0; # only virtual peptides were updated
			my @anaMatched=(keys %{$anaMatchedHash{$protID}});
			&updatePeptideCoverage($sthPP,$sthUpPC,$protID,length($proteinData{$protID}[0]),\@anaMatched);
		}
		$sthPP->finish;
		$sthUpPC->finish;
	}

	$dbh->commit;

	print " Done (",scalar keys %numPepUpdated," peptides updated, matching ",scalar keys %numProtMatched," proteins).\n";
}



sub correctDuplicateMaster {
	my ($dbh)=@_;
	print "+Checking for duplicate master proteins...";
	my ($uniIdID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='ID'");
	my %duplicateMasters;
	my $sthDM=$dbh->prepare("SELECT GROUP_CONCAT(ID_MASTER_PROTEIN ORDER BY ID_MASTER_PROTEIN DESC SEPARATOR ',') FROM MASTERPROT_IDENTIFIER WHERE ID_IDENTIFIER=$uniIdID GROUP BY VALUE HAVING COUNT(*) > 1");
	my $sthUpProt=$dbh->prepare("UPDATE PROTEIN SET ID_MASTER_PROTEIN=? WHERE ID_MASTER_PROTEIN=?");
	$sthDM->execute;
	my $numDuplicates=0;
	while (my ($masterIdStrg)=$sthDM->fetchrow_array) {
		$numDuplicates++;
		my @idList=split(',',$masterIdStrg);
		my $usedMasterID=shift @idList;
		foreach my $dupMasterID (@idList) {
			$sthUpProt->execute($usedMasterID,$dupMasterID);
			$duplicateMasters{$dupMasterID}=1;
		}
		print '.';
	}
	$sthDM->finish;
	$sthUpProt->finish;
	#my $numMasterDeleted=($numDuplicates)? &promsMod::deleteUnusedMasterProteins($dbh,\%duplicateMasters) : 0;
	#$dbh->commit;
	print " Done.\n-$numDuplicates duplicate master proteins found.\n";
}

####>Revision history<####
# 1.5.0 Major improvement in project-wide annotation update (PP 25/09/18)
# 1.4.6 Uses get instead of post for www.uniprot.org/mapping service (PP 18/09/18)
# 1.4.5 Uniprot response content is written to log file in case unexpected format (PP 02/07/18)
# 1.4.4 Changed a forgetten UniProt ACC old format to new one (PP 13/04/18)
# 1.4.3 Clears all mapping data prior to (re)mapping when update='force' & compatibility with new UniProt ACC format (PP 19/03/18)
# 1.4.2 Minor modification of &processUniprotAccession (GA 14/03/18)
# 1.4.1 Add wormbase as external resource in importAnnotationsFromUniprotKB (GA 30/08/17)
# 1.4.0 Now updates peptide position on protein when missing (PEP_BEG=0) (PP 03/03/17)
# 1.3.1 Minor change to display project id in log file during project remapping (PP 26/06/16)
# 1.3.0 Minor change to allow mapping of UNIPROT_ID with isoform tag (PP 08/02/16)
# 1.2.9 Fix bug preventing mapping in non-force case (PP 23/04/15)
# 1.2.8 writing in mapping_locatime log if update eq force (SL 19/03/15)
# 1.2.7 Fix bug undef $condSQL (PP 19/03/15)
# 1.2.6 add force update to update annotations (SL 26/01/15)
# 1.2.5 Remove optional ' {....}' after gene name in GN (PP 28/10/14)
# 1.2.4 Bug fix: UNIPROT_ACC -> UNIPROT_ACCESSION (PP 18/04/14)
# 1.2.3 Support for UPS standards mapping (PP 21/03/14)
# 1.2.2 Scan for missing protein sequence ('-') even if master protein is defined (PP 18/03/14)
# 1.2.1 Fix bug no sequence extraction from UniProtKB (PP 12/03/14)
# 1.2.0 Uses UniprotKB dbfetch instead of SRS & mkdir instead of make_path (PP 11/03/14)
# 1.1.1 Extends proxy settings to FTP protocol (PP 20/11/13)
# 1.1.0 Minor change (PP 13/11/13)
# 1.0.9 Better proxy declaration for LWP (PP 02/07/13)
# 1.0.8 Handles Flybase FB(gn|pp|tr) mapping (PP 22/04/13)
# 1.0.7 Tests if $flkAA & $organism are defined (PP 16/04/13)
# 1.0.6 Skips GI_ACCESSION mapping (PP 12/03/13)
# 1.0.5 Bug fix for missing FLANKING_AA in one of $sthPP queries (PP 01/03/13)
# 1.0.4 Handles Uniprot isoform number in ACC & added _&lt;year&gt; tag to mapping log file (PP 22/02/13)
# 1.0.3 Handles NCBI_ALL as identifier type (PP 28/01/13)
# 1.0.2 Flag on peptide coverage if peptide position &gt; protein length (PP 23/01/13)
# 1.0.1 Multi-databank searches (PP 12/121/12)
# 1.0.0 Replaces wgetMapIdentifiers.cgi & wgetMissingProteins.cgi (PP 22/11/12)
