#!/usr/local/bin/perl -w

################################################################################
# storeProjectItem.cgi         3.0.4                                           #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Stores/updates project items                                                 #
# Called by editProjectItem.cgi form submission                                #
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
#use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use strict;
use File::Copy;


#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
#my $date=strftime("%Y-%m-%d %H:%M:%S",localtime);

##################################
####>Watching a databank scan<#### For distant Mascot server. Called by promsMod::getProtInfo in verbose mode
##################################
if (param('WATCH')) {
	my $analysisID=param('ITEM_ID');
	my $numProteins=param('MAXPROT');
	my $dbSize=param('DBSIZE');
	my $anaParentDir=param('DIR');
	my $fastaFile="$anaParentDir/ana_$analysisID/analysis.fasta";
	my $endFile="$anaParentDir/ana_$analysisID/end.txt";
	print header(-'content-encoding'=>'no'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY topmargin='0'>
|;
	my $numEntries=0;
	my $prevRatio=0;
	print "<FONT class=\"title3\">0%";
	while (!-e $fastaFile) {sleep 1;} # wait for fasta if not created yet
	my $sleepTime=1+($dbSize/250000); $sleepTime=5 if $sleepTime>5;
	while (!-e $endFile && $prevRatio<100) {
		$numEntries=`grep -c '>' $fastaFile`;
		my $newRatio=10*(int(0.5+(10*$numEntries/$numProteins)));
		for (my $pc=$prevRatio+10;$pc<=$newRatio;$pc+=10) {print "...$pc%";}
		$prevRatio=$newRatio;
		sleep $sleepTime;
	}
	for (my $pc=$prevRatio+10;$pc<=100;$pc+=10) {print "...$pc%";}
	unlink $endFile;
	print " Done.</FONT>\n</BODY>\n</HTML>\n";
	exit;
}



##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

#print header; warningsToBrowser(1); # DEBUG
######################################
####>Updating Match Groups status<####
######################################
if (defined(param('MG'))) {
	my $newStatus=param('MG'); # 0/1
	my $analysisID=param('ITEM_ID');
	$dbh->do("UPDATE ANALYSIS SET VERIFIED_MG=$newStatus,UPDATE_DATE=NOW(),UPDATE_USER='$userID' WHERE ID_ANALYSIS=$analysisID");
	$dbh->commit;
	$dbh->disconnect;
	print redirect("./editProjectItem.cgi?ACT=summary&ITEM=ANALYSIS&ID=$analysisID");
	exit;
}


#######################
####>Starting HTML<####
#######################
print header(-'content-encoding'=>'no'); # start_html,"\n"; # start_html required to force update
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
|;

#############################
####>Fetching parameters<####
#############################
my $action=param('ACT');
my $item=uc(param('ITEM'));
my $parentItem=param('PARENT'); # add (not PROJECT)
my $parentID=param('PARENT_ID'); # add (not PROJECT)
my $itemID;
my @colName; my @colValue;
#my $maxRank=&promsConfig::getMaxRank; # Maximum number of interpretations allowed per query
#my $minScore=&promsConfig::getMinScore; # Minimum score allowed for interpretations
# my %massAAave=&promsConfig::getMassAAave; # Average mass, needed for protein mass calculation
# my %massATave=&promsConfig::getMassATave; # Average mass, needed for protein mass calculation
my @percent=(10,20,30,40,50,60,70,80,90,100); # needed to monitor process progression


##################################
####>Starting data processing<####
##################################
push @colName,('UPDATE_DATE','UPDATE_USER');
push @colValue,('NOW()',"'$userID'");

my $scanDB=(param('scanDB') && param('scanDB') eq 'now')? 1 : 0;
#if (param('scanDB') && param('scanDB') eq 'now') {
#	$scanDB=1;
#	print "<CENTER>\n<FONT class=\"title\">Importing ";
#	print "Analysis Data and " if $action eq 'add';
#	print "Protein Annotations</FONT><BR>\n";
#	print "<FONT class=\"title3\">You can add a new Analysis while the current one is being imported.</FONT><BR><BR>\n";
#	print "<IMG src=\"$promsPath{images}/engrenage.gif\">\n</CENTER><BR><BR>\n";
#}

if ($action eq 'add') {
	###>Fetching a new ID for item
	($itemID)=$dbh->selectrow_array("SELECT MAX(ID_$item) FROM $item");
	$itemID++;
	if ($item eq 'PROJECT') {
		push @colName,"ID_$item";
		push @colValue,$itemID;
	}
	push @colName,'START_DATE';
	push @colValue,'NOW()';
	if ($item ne 'PROJECT') {
		if ($item eq 'GEL2D') {
			my $imageFileName=(split(/[\\\/]/,param('image_file')))[-1];
			push @colName,'IMAGE_FILE';
			push @colValue,$dbh->quote($imageFileName);
			##>Copying image file to gelPath/gel2D_ID>###
			if (!-d $promsPath{gel_unix}) {
				mkdir $promsPath{gel_unix};
			}
			my $projectID = promsMod::getProjectID($dbh, $parentID, $parentItem);
			if (!-d "$promsPath{gel_unix}/project_$projectID") {
				mkdir "$promsPath{gel_unix}/project_$projectID";
			}

			my $gelFile="$promsPath{gel_unix}/project_$projectID/gel2D_$itemID.jpg";
			my $tmpfile = tmpFileName(upload('image_file')); # name of temp file being uploaded
			move($tmpfile,$gelFile);
			chmod 0664, $gelFile;
		}
		push @colName,"ID_$parentItem"; #&promsMod::getItemFkName($item);
		push @colValue,$parentID;
	}
}
else { # edit or reval
	$itemID=param('ITEM_ID');
}

my $itemStartingID=$itemID; # records the 1st ID (if addition of multiple entries)

my $name; # Global
if ($action ne 'reval') {
	$name=param('name');
	if ($action eq 'edit' || ($action eq 'add' && $item eq 'PROJECT')) {
		push @colName,'NAME';
		push @colValue,$dbh->quote($name);
		if ($item eq 'ANALYSIS') { # edit
			push @colName,'LAB_CODE';
			push @colValue,$dbh->quote(param('labCode'));
		}
	}
	if ($item eq 'EXPERIMENT') {
		my $prefSpeciesID=param('prefSpecies') || 'NULL';
		push @colName,'ID_SPECIES';
		push @colValue,$prefSpeciesID;
	}
	push @colName,('DES','COMMENTS');
	push @colValue,( $dbh->quote(param('des')), $dbh->quote(param('comments')) );
}
my ($projectID,$workgroup);
if ($item eq 'PROJECT') {
	push @colName,'PROT_VISIBILITY';
	push @colValue,param('protVis');
	push @colName,'ID_IDENTIFIER';
	my $identConv=(param('mapping'))? param('mapping') : undef;
	push @colValue,$dbh->quote($identConv);
	#push @colName,'RELEVANT_PTMS';
	#my $ptmString=join(';',param('relevantPTMs'));
	#$ptmString=undef unless $ptmString;
	#push @colValue,$dbh->quote($ptmString);
	push @colName,'OWNER';
	my $owner=param('owner');
	push @colValue,$dbh->quote($owner);
	push @colName,'WORK_GROUP';
	$workgroup=param('workgroup');
	$workgroup=param('newWorkgroup') if ($workgroup && $workgroup eq '#new#');
	push @colValue,$dbh->quote($workgroup);
	$projectID=$itemID;
}
else {
	$projectID=param('PROJECT_ID');
}

#my ($databankID); # Globals (used for analysis by at least 2 subs)
#if ($item eq 'ANALYSIS') { # edit or reval
#	#$msFileName=(split(/\\/,param('data_file')))[-1]; # works for add and edit
#	$databankID=param('databankID');
#	if (param('scanDB')) { # now or later (add or edit)
#		push @colName,'VALID_STATUS';
#		my $validStatus=($action eq 'reval')? 1 : ($scanDB)? 0 : -1; # 1:Partial validation, -1:Protein data not imported, 0:Protein data imported
#		push @colValue,$validStatus;
#	}
#	#if ($action eq 'edit') {
#	#	my $newQuantiID=(param('changeQuanti'))? param('changeQuanti') : 'NULL';
#	#	push @colName,'ID_QUANTIF_METHOD';
#	#	push @colValue,$newQuantiID;
#	#}
#	#$fileFormat=param('file_format');
#}


################################
####>Generating SQL Queries<####
################################
if ($action eq 'add') {
	my $colNameString=join(",",@colName);
	my $colValueString=join(",",@colValue);
	if ($item eq 'PROJECT') {
		$dbh->do("INSERT INTO $item ($colNameString) VALUES ($colValueString)") || die $dbh->errstr();
		###>Fetching user info
		my @userInfo=&promsMod::getUserInfo($dbh,$userID);
		###>Updating PROJECT_ACCESS table
		#if ($userInfo[1] eq 'bio') { # user is a biologist
		#	my $accessID=(scalar @{$userInfo[7]})? 32 : 30; # if (mascotIDs) => super_admin else => admin
		#	$dbh->do("INSERT INTO PROJECT_ACCESS (ID_USER,ID_PROJECT,ID_PROFILE,UPDATE_DATE,UPDATE_USER) VALUES ('$userID',$projectID,$accessID,'$date','$userID')") || die $dbh->errstr(); # 30=administrator
		#}
	}
	else {
		###>Generating a new DISPLAY_POS default value<###
		my ($displayPos) = $dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM $item WHERE ID_$parentItem=$parentID");
		$displayPos++;

		###>Processing labels (EXPERIMENT or SAMPLE)<###
		my %labelList;
		if (param('labels')) {
			my $labelString=param('labels');
			$labelString=~s/[^\d,-]//g;
			my @tempList=split(',',$labelString);
			foreach my $label (@tempList) {
				my ($beg,$end)=split('-',$label);
				$labelList{$beg}=1;
				if ($end) {
					for (my $l=$beg+1;$l<=$end;$l++) {
						$labelList{$l}=1;
					}
				}
			}
		}

		if (scalar keys %labelList) { # valid label(s)
			$name=&promsMod::resize($name,40); # instead of 50 just to safe
			foreach my $label (sort{$a<=>$b} keys %labelList) {
				my $labelName="$name$label";
				$labelName=$dbh->quote($labelName);
# print "INSERT INTO $item (ID_$item,NAME,DISPLAY_POS,$colNameString) VALUES ($itemID,$labelName,$displayPos,$colValueString)<BR>\n";
				$dbh->do("INSERT INTO $item (ID_$item,NAME,DISPLAY_POS,$colNameString) VALUES ($itemID,$labelName,$displayPos,$colValueString)") || die $dbh->errstr();
				$itemID++;
				$displayPos++;
			}
		}
		else { # no valid labels => only 1 entry
			$name=$dbh->quote($name);
#print "INSERT INTO $item (ID_$item,NAME,DISPLAY_POS,$colNameString) VALUES ($itemID,$name,$displayPos,$colValueString)<BR>\n";
			$dbh->do("INSERT INTO $item (ID_$item,NAME,DISPLAY_POS,$colNameString) VALUES ($itemID,$name,$displayPos,$colValueString)") || die $dbh->errstr();
		}
		#$dbh->commit if $scanDB; # id_analysis is protected as soon as possible!
	}
}
else { # action = edit or reval

	###>Identifier mapping<###
	if ($item eq 'PROJECT') {
		$workgroup='' unless $workgroup;
		my ($oldWorkgroup,$oldMapping)=$dbh->selectrow_array("SELECT WORK_GROUP,ID_IDENTIFIER FROM PROJECT WHERE ID_PROJECT=$itemID");
		$oldWorkgroup='' unless $oldWorkgroup;
		$oldMapping='' unless $oldMapping;
		#>Delete new workgroup managers project-based access (no longer required)
		if ($workgroup ne $oldWorkgroup) { # change in workgroup
			my $wgSelectStrg=($workgroup)? "WORK_GROUP='$workgroup'" : 'WORK_GROUP IS NULL';
			my $sthWGM=$dbh->prepare("SELECT ID_USER FROM USER_LIST WHERE USER_STATUS='manag' AND $wgSelectStrg");
			my $sthDelPA=$dbh->prepare("DELETE FROM PROJECT_ACCESS WHERE ID_PROJECT=$itemID AND ID_USER=?");
			$sthWGM->execute;
			while (my ($uID)=$sthWGM->fetchrow_array) {
				$sthDelPA->execute($uID);
			}
			$sthWGM->finish;
			$sthDelPA->finish;
		}
		my $newMapping=(param('mapping'))? param('mapping') : 0;
		if ($newMapping ne $oldMapping) {
			##>Convert identifiers to Gene Symbols
			if ($newMapping) {
				my $sthUpP=$dbh->prepare("UPDATE PROTEIN SET ALIAS=? WHERE ID_PROTEIN=?");
				my ($giIdentID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GI'");
				if ($newMapping==$giIdentID) { # special case: info is extracted from default identifier
					my $sthDI=$dbh->prepare("SELECT ID_PROTEIN,IDENTIFIER FROM PROTEIN WHERE ID_PROJECT=$itemID");
					$sthDI->execute;
					while (my ($protID,$identifier)=$sthDI->fetchrow_array) {
						my ($alias)=($identifier=~/(gi\|\d+)/)? $1 : $identifier;
						$sthUpP->execute($alias,$protID);
					}
					$sthDI->finish;
				}
				else { # info in master protein
					my $sthMI=$dbh->prepare("SELECT VALUE,ID_PROTEIN FROM PROTEIN P,MASTERPROT_IDENTIFIER MI WHERE ID_PROJECT=$itemID AND ID_IDENTIFIER=$newMapping AND P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND RANK=1");
					$sthMI->execute;
					while (my ($mapIdent,$protID)=$sthMI->fetchrow_array) {
						$sthUpP->execute($mapIdent,$protID);
					}
					$sthMI->finish;
				}
				$sthUpP->finish;
			}
			##>Restore default identifiers
			else {
				$dbh->do("UPDATE PROTEIN SET ALIAS=IDENTIFIER WHERE ID_PROJECT=$itemID"); # most frequent case

				##>Check if identifier is modified for alias
				#my %identModString=&promsConfig::getIdent2AliasModString;
				#my $sthAna=$dbh->prepare("SELECT ID_ANALYSIS FROM EXPERIMENT,SAMPLE,ANALYSIS WHERE ID_PROJECT=$itemID AND EXPERIMENT.ID_EXPERIMENT=SAMPLE.ID_EXPERIMENT AND SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE");
				#my $sthIT=$dbh->prepare("SELECT IDENTIFIER_TYPE FROM ANALYSIS,DATABANK WHERE DATABANK.ID_DATABANK=ANALYSIS.ID_DATABANK AND ID_ANALYSIS=?");
				#my $sthAP=$dbh->prepare("SELECT ID_PROTEIN,IDENTIFIER FROM PROTEIN WHERE ID_PROTEIN IN (SELECT ID_PROTEIN FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=?)");
				#$sthAna->execute;
				#while (my ($anaID)=$sthAna->fetchrow_array) {
				#	$sthIT->execute($anaID);
				#	my ($identifierType)=$sthIT->fetchrow_array;
				#	next unless $identModString{$identifierType}; # alias=identifier
				#	$sthAP->execute($anaID);
				#	while (my ($protID,$identifier)=$sthAP->fetchrow_array) {
				#		my $alias=&promsMod::modifyIdentifier($identifier,$identModString{$identifierType});
				#		$sthUpP->execute($alias,$protID);
				#	}
				#}
				#$sthAna->finish;
				#$sthIT->finish;
				#$sthAP->finish;
			}
		}
	}
	elsif ($item eq 'SPOT') {
		my $isoPoint=(param('isoPoint')=~/^[\d+|\.]+\Z/)? param('isoPoint') : 'NULL'; # number or NULL
		my $molWeight=(param('molWeight')=~/^[\d+|\.]+\Z/)? param('molWeight') : 'NULL'; # number or NULL
		my $intensity=(param('intensity')=~/^[\d+|\.]+\Z/)? param('intensity') : 'NULL'; # number or NULL
		my $externalID=(param('externalID')=~/.+/)? $dbh->quote(param('externalID')) : 'NULL'; # any or NULL
		push @colName,('ISOELECTRIC_POINT','MOLECULAR_WEIGHT','INTENSITY','EXTERNAL_ID');
		push @colValue,($isoPoint,$molWeight,$intensity,$externalID);
	}
	my $updateString;
	for (my $c=0;$c<=$#colName;$c++) { # looping through all conditions
		$updateString.="$colName[$c]=$colValue[$c]";
		$updateString.=',' if $c<$#colName;
	}

	$dbh->do("UPDATE $item SET $updateString WHERE ID_$item=$itemID") || die $dbh->errstr();

	###>Update Spot sample<###
	if ($item eq 'SPOT') {
		my $newSampleID=param('assoSampleID');
		my $oldSampleID=param('oldSampleID');
		if ($newSampleID != $oldSampleID) {
			my $sthUpSamp=$dbh->prepare("UPDATE SAMPLE SET ID_SPOT=?,DISPLAY_POS=? WHERE ID_SAMPLE=?");
			if ($oldSampleID || $newSampleID>0) {
				my $refSampleID=($oldSampleID)? $oldSampleID : $newSampleID;
				my $sthDP=$dbh->prepare("SELECT ID_SAMPLE FROM SAMPLE WHERE ID_EXPERIMENT=(SELECT ID_EXPERIMENT FROM SAMPLE WHERE ID_SAMPLE=$refSampleID) AND ID_SPOT IS NULL ORDER BY DISPLAY_POS ASC");
				my $sthUpDP=$dbh->prepare("UPDATE SAMPLE SET DISPLAY_POS=? WHERE ID_SAMPLE=?");
				my $maxDisPos=0; #$dhb->selectrow_array("SELECT MAX(DISPLAY_POS) FROM SAMPLE WHERE ID_SPOT=$itemID");
				$sthDP->execute;
				while (my ($brotherID)=$sthDP->fetchrow_array) { # never ==oldSampleID
					next if $brotherID==$newSampleID; # skip new associated sample
					$sthUpDP->execute(++$maxDisPos,$brotherID); # update brother display
				}
				$sthDP->finish;
				$sthUpDP->finish;
				if ($oldSampleID) {
					my ($numAna)=$dbh->selectrow_array("SELECT COUNT(*) FROM ANALYSIS WHERE ID_SAMPLE=$oldSampleID");
					if ($numAna) {$sthUpSamp->execute(undef,++$maxDisPos,$oldSampleID);} # dissociate old sample from spot & put last in list
					else {$dbh->do("DELETE FROM SAMPLE WHERE ID_SAMPLE=$oldSampleID");} # delete empty sample
				}
			}
			if ($newSampleID==-1) { # new sample created on the fly
				my $name=$dbh->quote(param('newSampleName'));
				my ($maxSampleID)=$dbh->selectrow_array("SELECT MAX(ID_SAMPLE) FROM SAMPLE");
				$maxSampleID++;
				my ($expID)=$dbh->selectrow_array("SELECT ID_EXPERIMENT FROM GEL2D WHERE ID_GEL2D=(SELECT ID_GEL2D FROM SPOT WHERE ID_SPOT=$itemID)");
				$dbh->do("INSERT INTO SAMPLE (ID_SAMPLE,ID_EXPERIMENT,ID_SPOT,NAME,START_DATE,UPDATE_DATE,UPDATE_USER) VALUES ($maxSampleID,$expID,$itemID,$name,NOW(),NOW(),'$userID')");
			}
			elsif ($newSampleID==0) { # no associated sample
				my ($maxDisPos)=$dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM SAMPLE WHERE ID_EXPERIMENT=(SELECT ID_EXPERIMENT FROM SAMPLE WHERE ID_SAMPLE=$oldSampleID) AND ID_SPOT IS NULL");
				$sthUpSamp->execute(undef,++$maxDisPos,$newSampleID);
			}
			else {$sthUpSamp->execute($itemID,undef,$newSampleID);} # associate new sample (if exists) with spot
			$sthUpSamp->finish;
		}
	}

	###>Updating displayPos<###
	my ($newDisplayPos,$oldDisplayPos)=(param('displayPos'),param('oldDisplayPos'));
	if ($newDisplayPos && $newDisplayPos != $oldDisplayPos) { # edit + EXP or SAMP or GEL or ANA
		my @itemParents=&promsMod::getItemParent($item);
		my $sthDP=$dbh->prepare("SELECT ID_$item,DISPLAY_POS FROM $item WHERE ID_$itemParents[0]=(SELECT ID_$itemParents[0] FROM $item WHERE ID_$item=$itemID)");
		my $sthUpDP=$dbh->prepare("UPDATE $item SET DISPLAY_POS=? WHERE ID_$item=?");
		$sthDP->execute;
		while (my ($brotherID,$dispPos)=$sthDP->fetchrow_array) {
			next unless $dispPos; # Sample is associated with a Spot
			if ($brotherID==$itemID) {
				$sthUpDP->execute($newDisplayPos,$brotherID);
			}
			elsif ($newDisplayPos > $oldDisplayPos) {
				$sthUpDP->execute($dispPos-1,$brotherID) if ($dispPos > $oldDisplayPos && $dispPos <= $newDisplayPos);
			}
			else { # new < old
				$sthUpDP->execute($dispPos+1,$brotherID) if ($dispPos >= $newDisplayPos && $dispPos < $oldDisplayPos);
			}
		}
		$sthDP->finish;
		$sthUpDP->finish;
	}

	###>Update Proteins visibility (only for Projects)<###
	my $protVisibility=param('protVis');
	my $oldProtVisibility=param('oldProtVis');
	if ($item eq 'PROJECT' && (param('upVis') || $protVisibility != $oldProtVisibility)) {
		&promsMod::applyVisibilityRule($dbh,$projectID); # update classifications as well
	}
}

###>Project relevant PTMs<###
if ($item eq 'PROJECT' && $action=~/add|edit/) {
	###> Update PROJECT_MODIFICATION table
	$dbh->do("DELETE FROM PROJECT_MODIFICATION WHERE ID_PROJECT=$itemID") if $action eq 'edit';
	foreach my $modID (param('relevantPTMs')) {
		$modID=int($modID);
		$dbh->do("INSERT INTO PROJECT_MODIFICATION (ID_PROJECT,ID_MODIFICATION) VALUES($itemID,$modID)");
	}
}


########################################
####>Scaning databank (MS analysis)<####
########################################
#my (%protList,%protDes,%protMW,%protOrg,%protLength); # Globals (used for analysis by at least 2 subs)
#my (%massExp, %massObs,%queryInfo,%maxQueryScore,%numValid,%rankProtMatch,%matchList,%elutionTime,%matchGroup,%protSelStatus,%maxProtScore); #useless
#my ($queryID,$protValID); # Globals (used for analysis by at least 2 subs)

if ($scanDB) { # Edit + Import Now
	system "./scanDatabank.pl $userID $itemID"; # background scan
	print "<BR><CENTER><FONT class=\"title3\">Protein annotations are being imported in background.</FONT></CENTER>\n";

	#my $sthProt=$dbh->prepare("SELECT IDENTIFIER,ID_PROT_VALID FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$itemID AND IDENTIFIER NOT LIKE 'DECOY_%'");
	#$sthProt->execute;
	#while (my($identifier,$protValidID)=$sthProt->fetchrow_array) {
	#	$protList{$identifier}=$protValidID;
	#}
	#$sthProt->finish;
	#
	#####>Extracting data from databank file<####   ($dbh,$databankID, $itemID,$type, $protDesRef, $protMWRef,$protOrgRef, $protLengthRef, %maxProtMatch)
	##&promsMod::getProtInfo_old($dbh,$databankID,$itemID,'verbose',\%protDes,\%protMW,\%protOrg,\%protLength,\%protList);
	#&promsMod::getProtInfo('verbose',$dbh,$databankID,[$itemID],\%protDes,\%protMW,\%protOrg,\%protLength,\%protList);
	#
	#####>Updating PROTEIN_VALIDATION table<####
	#print "<FONT class=\"title2\">Updating database:<BR>\n";
	#my $sthUp=$dbh->prepare("UPDATE PROTEIN_VALIDATION SET MW=?,PROT_LENGTH=?,PROT_DES=?,ORGANISM=? WHERE ID_PROT_VALID=?");
	#my $numEntry=scalar keys (%protList);
	#my @limitValue;
	#foreach my $pc (@percent) {push @limitValue,int(0.5+($numEntry*$pc/100));}
	#my $index=0;
	#my $counter1=0;
	#my $counter2=0;
	#my $maxCounter2=int(0.5+($numEntry/100));
	#print "<FONT class=\"title3\">0%";
	#foreach my $identifier (keys %protList) {
	#	my $des=&promsMod::resize($protDes{$identifier},250); # max length allowed in table
	#	my $organism=&promsMod::resize($protOrg{$identifier},100); # max length allowed in table
	#	$sthUp->execute($protMW{$identifier},$protLength{$identifier},$des,$organism,$protList{$identifier}) || die $sthUp->errstr();
	#	$counter1++;
	#	if ($counter1>=$limitValue[$index]) {
	#		print "$percent[$index]%"; # keeping connection alive
	#		$index++;
	#	}
	#	$counter2++;
	#	if ($counter2==$maxCounter2) {print '.'; $counter2=0;}
	#}
	#$sthUp->finish;
	#print "</FONT><BR>Done.</FONT><BR><BR>\n";
}

###>Checking which frame to update<###
$item=lc($item); #!!!! lower case now!!!
my $gelID;
if ($item eq 'spot') {
	($gelID)=$dbh->selectrow_array("SELECT ID_GEL2D FROM SPOT WHERE ID_SPOT=$itemStartingID");
}
elsif ($item eq 'sample') {
	($gelID)=$dbh->selectrow_array("SELECT ID_GEL2D FROM SAMPLE,SPOT WHERE SAMPLE.ID_SPOT=SPOT.ID_SPOT AND ID_SAMPLE=$itemStartingID");
}
elsif ($item eq 'analysis') {
	($gelID)=$dbh->selectrow_array("SELECT ID_GEL2D FROM ANALYSIS,SAMPLE,SPOT WHERE SPOT.ID_SPOT=SAMPLE.ID_SPOT AND SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND ID_ANALYSIS=$itemStartingID");
}
$gelID=0 unless $gelID;
my $targetFrame=($gelID && $item ne 'spot')? 'itemFrame' : 'navFrame';

$dbh->commit;
$dbh->disconnect;
# print "<H3>All tables loaded!</H3>\n";
sleep 2 if $scanDB;
#exit;  #debug
print qq
|<SCRIPT LANGUAGE="JavaScript">
top.promsFrame.selectedAction='summary';
if ('$item'=='project' \|\| ('$item'=='gel2d' && !parent.itemFrame)) {
	top.promsFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=$item:$itemStartingID&ACT=open"; // top.promsFrame but ! parent in case new project!
}
else if ('$targetFrame'=='navFrame') {
	if ('$item'=='spot') { // reload navFrame in case modified sample association
		parent.navFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=gel2d:$gelID&itemBranchID=$item:$itemID&ACT=nav&VIEW="+parent.navFrame.view;
	}
	else {
		parent.navFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=$item:$itemStartingID&ACT=nav&VIEW="+parent.navFrame.view;
	}
}
else { //itemFrame
	parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID&GEL=gel2d:$gelID&branchID=$item:$itemStartingID&ACT=gel";
}
</SCRIPT>
</BODY>
</HTML>
|;

####>Revision history<####
# 3.0.4 Minor bug fix for Add Experiment with Preferred species (PP 04/08/17)
# 3.0.3 Restores previously commented management of parameter WATCH called by promsMod::getProtInfo (PP 04/01/17)
# 3.0.2 Bug fix for Experiment when preferred species is undef (PP 15/03/16)
# 3.0.1 Added Preferred species selection for Experiment (PP 01/03/16)
# 3.0.0 Lab Code for Analysis (PP 24/07/14)
# 2.9.9 Syntax error fix (PP 18/11/13)
# 2.9.8 system command removal (PP 12/11/13)
# 2.9.7 Creating "gels" directory if it does not exist, and uploading gels into a project subdir (FY 02/09/13)
# 2.9.6 Fix bug creation project when relevant PTMs selected (PP 28/06/13)
# 2.9.5 Minor change for relevant PTMs update -> save Project relevant PTMs (GA 06/05/13)
# 2.9.4 Code cleaning related to databank scan (PP 29/04/13)
# 2.9.3 Bug fix due to deleted ANALYIS.ID_QUANTIF_METHOD (PP 17/04/13)
# 2.9.2 Handles GI_ACCESSION for identifier conversion (PP 12/03/13)
# 2.9.1 New identifier mapping procedure (PP 20/11/12)
# 2.9.0 Delete manager's project-based access when Project's workgroup matches manager's one<BR>& new call &getProtInfo (PP 13/08/12)
# 2.8.9 Minor update -> new parameter (VIEW) given to openProject.cgi for QUANTI/GO mode (GA 08/12/11)
# 2.8.8 header(-'content-encoding'=>'no') (PP 02/11/11)
# 2.8.7 Workgroup management & cancels 2.8.6  (PP 20/09/11)
# 2.8.6 handles project creation by super biologist (PP 22/08/11)
