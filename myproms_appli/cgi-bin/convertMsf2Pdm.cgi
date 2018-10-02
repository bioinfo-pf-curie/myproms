#!/usr/local/bin/perl -w

################################################################################
# convertMsf2Pdm.cgi     1.2.4                                                 #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Converts Proteome Discoverer .msf file into a Mascot-like .pdm file          #
# Called by storeAnalyses.cgi                                                  #
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
$| = 1; # only for STDIN
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use IO::Handle; # to flush file buffer
use promsConfig;
use promsMod;
use strict;
use XML::Simple;


#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'} || 'unknown';

#############################
####>Fetching parameters<####
#############################
my $processingNodeNumber=param('node'); # 4; #
my $msfFile=param('file'); # 'SILAC-EF1j-Top10CID-RIPA-10B-uni.msf'; #
my $filePath=param('path') || ''; # '/bioinfo/projects_dev/myproms/ppoullet/tmp/batch/ppoullet'; #
my $fileID=param('fileID')?param('fileID') : '' ; # '108'; # list of rawfiles in FileInfos SQLite table

$filePath=~s/\/\Z//; # removes ending '/' if exists
$filePath="/$filePath" if $filePath !~ /^\//; # adds '/' to path
if (param('error')) {
	print header(-'content-encoding'=>'no',-charset=>'utf-8'); # start_html,"\n"; # start_html required to force update
	warningsToBrowser(1);
	my $processError;
	if ($fileID) {
		($processError=$msfFile)=~s/\.msf\Z/_$processingNodeNumber\.$fileID\.$userID\.error/;
	}else{
		($processError=$msfFile)=~s/\.msf\Z/_$processingNodeNumber\.$userID\.error/;
	}
	open (ERROR,">$promsPath{tmp}/$processError");
	print ERROR '__ERROR__';
	close ERROR;
	print "<HTML></HTML>\n";
	exit;
}
my $selMinScore=param('minScore') || 0; # 'default'; #
my $percolatorThres=(param('percolThr'))? param('percolThr')/100 : undef; # 1% -> 0.01
my $databankID=param('databankID'); # 17; #


#######################
####>Starting HTML<####
#######################
print header(-'content-encoding'=>'no',-charset=>'utf-8'); # start_html,"\n"; # start_html required to force update
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Convert MSF to PDM</TITLE>
<SCRIPT language="javascript">
/*
################################################################################
################################################################################
Required so that IE can start streaming the output
################################################################################
################################################################################
*/
</SCRIPT>
</HEAD>
<BODY>
|;

##############################
####>Connect to SQLite DB<####
##############################
my $dbsqlite = DBI->connect("dbi:SQLite:$filePath/$msfFile", "", "", {PrintError => 1,RaiseError => 1});

###################################################
####>Parsing Sequest .msf file (sqlite format)<####
###################################################
my ($proteomeDiscovererVersion)=$dbsqlite->selectrow_array("SELECT SoftwareVersion FROM SchemaInfo ORDER BY rowid ASC LIMIT 1"); # last record
$proteomeDiscovererVersion=~s/^(\d\.\d+).*/$1/; # x.x.x.xx -> x.x (numerical value)

my $refSearches=&getAllSearches;

#(my $pdmFile="$filePath/$msfFile")=~s/\.msf\Z/_$processingNodeNumber\.pdm/;
my $usedFileName;
if ($fileID) {
	($usedFileName=$msfFile)=~s/\.msf\Z/_$processingNodeNumber\.$fileID\.pdm/;
}else{
	($usedFileName=$msfFile)=~s/\.msf\Z/_$processingNodeNumber\.pdm/;
}

my $pdmFile="$promsPath{tmp}/batch/$userID/$usedFileName"; # pdm generated in user's directory no matter where msf is

my (%vmodsInUnimod,%queryOrder,%spectrumInfo,%modifications,%masses); # Globals (used for analysis by at least 2 subs)  @queryToSpectrum,@queryRT,@queryCharge,%protSelStatus


####<Open PDM file>####
open (FILE,">$pdmFile");
FILE->autoflush; # <- IO::Handle

my $boundary='--gc0p4Jq0M2Yt08jU534c0p';


#########################################
####> Specific functions for PD 2.2 <####
#########################################
if ($proteomeDiscovererVersion >= 2.2) {
		my ($dbFile,$decoy,$minScorePep)=&printParametersMSF2_2(\%modifications,\%masses);
		#my $addSpInfoStg=($fileID)?" AND StudyFileId=$fileID":" "; # If a merge was done
		my $addSpInfoStg=($fileID)?" WHERE SpectrumFileID=$fileID":" "; # If a merge was done
		my $sthSpInfo = $dbsqlite->prepare("SELECT SpectrumID, Mass, Charge, MassOverCharge FROM MSnSpectrumInfo $addSpInfoStg ORDER BY Mass ASC, Charge ASC");# Change the order to be like mascot
		$sthSpInfo->execute;
		my $queryNum=0;
		while (my ($spectrumID,$mass,$charge,$massPeak)= $sthSpInfo->fetchrow_array) {
			$queryOrder{++$queryNum}=$spectrumID;
			@{$spectrumInfo{$spectrumID}}=(
				$mass-$masses{'Hydrogen'}+$masses{'Electron'},	# 0 MR
				$charge,										# 1 CHARGE
				$massPeak										# 2 MASSPEAK
			)
		}
		$sthSpInfo->finish;
		&printHeaderMSF2_2($dbFile,$queryNum);
		&printSummaryMSF($decoy,\%queryOrder,\%spectrumInfo);
		&printPeptidesMSF2_2($decoy,$minScorePep,\%queryOrder,\%spectrumInfo,\%modifications,\%masses);
		close FILE;
		$dbsqlite->disconnect;
		my $processEnd;
		if ($fileID) {
			($processEnd=$msfFile)=~s/\.msf\Z/_$processingNodeNumber\.$fileID\.$userID\.end/;
		}else{
			($processEnd=$msfFile)=~s/\.msf\Z/_$processingNodeNumber\.$userID\.end/;
		}
		
		open (P_END,">$promsPath{tmp}/$processEnd");
		print P_END "__END__\n";
		close P_END;
		
		print "<!--__OK__--></BODY>\n</HTML>\n";
		
		exit;
}

####<1st step: Print header section & load hashes>####
my ($dbFile,$decoy,$minScorePep)=&printParametersMSF(\%modifications,\%masses);

####<Extraction of all spectra & sort into queries>####
#my $sthQO = $dbsqlite->prepare("SELECT SpectrumID FROM SpectrumHeaders, MassPeaks, Spectra WHERE Spectra.UniqueSpectrumID=SpectrumHeaders.UniqueSpectrumID AND SpectrumHeaders.MassPeakID=MassPeaks.MassPeakID ORDER BY SpectrumHeaders.Mass asc");
#my $sthQO = $dbsqlite->prepare("SELECT SpectrumID FROM SpectrumHeaders ORDER BY Charge ASC, Mass ASC");
#my $sthSpInfo = $dbsqlite->prepare("SELECT SpectrumID, SpectrumHeaders.Mass, MassPeaks.Charge, MassPeaks.Mass FROM SpectrumHeaders, MassPeaks WHERE SpectrumHeaders.MassPeakID=MassPeaks.MassPeakID ORDER BY MassPeaks.Charge ASC, MassPeaks.Mass ASC");
#my $sthSpInfo = $dbsqlite->prepare("SELECT SpectrumID, SpectrumHeaders.Mass, SpectrumHeaders.Charge, MassPeaks.Mass FROM SpectrumHeaders, MassPeaks WHERE SpectrumHeaders.MassPeakID=MassPeaks.MassPeakID ORDER BY MassPeaks.Charge ASC, MassPeaks.Mass ASC");
my $addSpInfoStg=($fileID)?" AND FileID=$fileID":" "; # If a merge was done
my $sthSpInfo = $dbsqlite->prepare("SELECT SpectrumID, SpectrumHeaders.Mass, SpectrumHeaders.Charge, MassPeaks.Mass FROM SpectrumHeaders, MassPeaks WHERE SpectrumHeaders.MassPeakID=MassPeaks.MassPeakID$addSpInfoStg ORDER BY SpectrumHeaders.Mass ASC, MassPeaks.Charge ASC");# Change the order to be like mascot
$sthSpInfo->execute;
my $queryNum=0;
while (my ($spectrumID,$mass,$charge,$massPeak)= $sthSpInfo->fetchrow_array) {
	$queryOrder{++$queryNum}=$spectrumID;
	@{$spectrumInfo{$spectrumID}}=(
		$mass-$masses{'Hydrogen'}+$masses{'Electron'},	# 0 MR
		$charge,										# 1 CHARGE
		$massPeak										# 2 MASSPEAK
	)
}
$sthSpInfo->finish;

####<2nd step: Print header section>####
&printHeaderMSF($dbFile,$queryNum);

####<3rd step: Print summary section>####
&printSummaryMSF($decoy,\%queryOrder,\%spectrumInfo);

####<4th step: Print peptides section & fasta file>####
&printPeptidesMSF($decoy,$minScorePep,\%queryOrder,\%spectrumInfo,\%modifications,\%masses);

close FILE;

$dbsqlite->disconnect;

#close OUT;
my $processEnd;
if ($fileID) {
	($processEnd=$msfFile)=~s/\.msf\Z/_$processingNodeNumber\.$fileID\.$userID\.end/;
}else{
	($processEnd=$msfFile)=~s/\.msf\Z/_$processingNodeNumber\.$userID\.end/;
}

open (P_END,">$promsPath{tmp}/$processEnd");
print P_END "__END__\n";
close P_END;

print "<!--__OK__--></BODY>\n</HTML>\n";

exit;

####################<<< SUBROUTINES >>>#########################

#############################################################
####<Fetch search params for all Processing Node Numbers>####
#############################################################
sub getAllSearches { # GLOBALS: $dbsqlite
	my %searches;
	#my ($rawfilenumber,$processingNodeNumber);
	#my $typeInfo="";
#!!	my ($rawfilenumber) = $dbsqlite->selectrow_array("SELECT ProcessingNodeNumber FROM ProcessingNodes WHERE FriendlyName='Spectrum Files';"); #commented by PP

	#my $sth = $dbsqlite->prepare("SELECT ProcessingNodeNumber, ProcessingNodeID, NodeComment FROM ProcessingNodes WHERE NodeName='SequestNode' OR NodeName='Mascot';");
	#my $sth = $dbsqlite->prepare("SELECT ProcessingNodeNumber, ProcessingNodeID, FriendlyName FROM ProcessingNodes WHERE (NodeName='SequestNode' OR NodeName='Mascot')");# Change on 25/04/12 for PD 1.2 and 1.3, NodeComment was empty in PD1.3 version

	###> For PD2.0 -> ProcessingNodes table does not exist anymore and this information is kept in Workflows table as an XML file
	###> For PD1.4 -> new name in the database in ProcessingNodes: NodeName='IseNode' & FriendlyName='Sequest HT' for SEQUEST searches
	my $sthGetParams;
	if ($proteomeDiscovererVersion >= 2.0) {
		my ($wfXML) = $dbsqlite->selectrow_array("SELECT WorkflowXML FROM Workflows");
		my $xml = new XML::Simple();
		my $xmlData = $xml->XMLin($wfXML);
		my ($parenProcessingNumber,$processingNodeName,$processingNodeNumber);
		my $prsNodeID=0;
		foreach my $processingNodeParameters (sort{$a->{ParentProcessingNodeNumber} cmp $b->{ParentProcessingNodeNumber}} @{$xmlData->{WorkflowTree}->{WorkflowNode}}) {
			next unless $processingNodeParameters->{FriendlyName} eq 'Mascot' or uc($processingNodeParameters->{FriendlyName}) =~ 'SEQUEST';
			$prsNodeID++;
			$processingNodeNumber=$processingNodeParameters->{ProcessingNodeNumber};
			$searches{$processingNodeNumber}{'NODEID'}=$prsNodeID;
			$searches{$processingNodeNumber}{'DESCRIPTOR'}=$processingNodeParameters->{FriendlyName};
			$searches{$processingNodeNumber}{'NODEGUID'}=$processingNodeParameters->{Guid};
		}
		foreach my $processingNodeParameters (sort{$a->{ParentProcessingNodeNumber} cmp $b->{ParentProcessingNodeNumber}} @{$xmlData->{WorkflowTree}->{WorkflowNode}}) {
			next unless $processingNodeParameters->{FriendlyName} eq 'Percolator' or uc($processingNodeParameters->{FriendlyName}) =~ /PHOSPHORS/;
			$prsNodeID++;
			my ($tagName)=($processingNodeParameters->{FriendlyName} eq 'Percolator')? 'PERCOLATOR':'PHOSPHORS';
			foreach my $parNode (split(';',$processingNodeParameters->{ParentProcessingNodeNumber})) {
				if ($searches{$parNode}) {
					$searches{$parNode}{"${tagName}_NODE"}=$processingNodeParameters->{ProcessingNodeNumber};
					$searches{$parNode}{"${tagName}_NODEID"}=$prsNodeID;
				}
			}
		}
		return (\%searches);
	}
	elsif ($proteomeDiscovererVersion >= 1.4 && $proteomeDiscovererVersion < 2.0) {
		$sthGetParams = $dbsqlite->prepare("SELECT ProcessingNodeNumber, ProcessingNodeID, FriendlyName FROM ProcessingNodes WHERE (UPPER(FriendlyName) LIKE 'SEQUEST%' OR NodeName='Mascot')");
	}
	else {
		$sthGetParams = $dbsqlite->prepare("SELECT ProcessingNodeNumber, ProcessingNodeID, FriendlyName FROM ProcessingNodes WHERE (NodeName='SequestNode' OR NodeName='Mascot')");
	}
	#my $sthPercol = $dbsqlite->prepare("SELECT ProcessingNodeNumber,ProcessingNodeID,ProcessingNodeParentNumber FROM ProcessingNodes WHERE FriendlyName='Percolator'");
	my $sthPercolPRS = $dbsqlite->prepare("SELECT ProcessingNodeNumber,ProcessingNodeID,ProcessingNodeParentNumber,FriendlyName FROM ProcessingNodes WHERE FriendlyName LIKE 'Percolator%' OR FriendlyName LIKE '%phosphoRS%'");
	$sthGetParams->execute;
	while (my ($procNodeNumber,$processingNodeID,$descriptor) = $sthGetParams->fetchrow_array){
		$searches{$procNodeNumber}{'NODEID'}=$processingNodeID;
		$searches{$procNodeNumber}{'DESCRIPTOR'}=$descriptor;
		$sthPercolPRS->execute;
		while (my ($percolprsNodeNumber,$percolprsNodeID,$parentNodeStrg,$friendlyName)=$sthPercolPRS->fetchrow_array) {
			my ($tagName)=($friendlyName =~ /^Percolator/)? 'PERCOLATOR':'PHOSPHORS'; # Percolator[ (64Bit)]
			foreach my $parNode (split(';',$parentNodeStrg)) {
				if ($parNode==$procNodeNumber) {
					$searches{$procNodeNumber}{"${tagName}_NODE"}=$percolprsNodeNumber;
					$searches{$procNodeNumber}{"${tagName}_NODEID"}=$percolprsNodeID;
				}
			}
		}
	}
	$sthGetParams->finish;
	$sthPercolPRS->finish;

	return (\%searches);
}

##############################################
####<Print parameters section (.pdm file)>####
##############################################
### Be careful!!! In mascot header file, if there is a fixed modification on an amino acid,
### the mass is updated so as to incorporate this fixed modification in the header
### For example: a carbamidomethyl on cysteine residue add 57.021464 to the mass
### of the cysteine which is 103.01. Then, it will be written that the mass of
### cystein is 160.030 which is true and false. True because that is the mass that MASCOT
### used for that amino acid. False because it is not the true mass of this amino acid.
### In the dat file produced by this script, the mass written is the mass of each amino-acid without
### any fixe modification.
sub printParametersMSF { # GLOBALS: $dbsqlite, $processingNodeNumber, $refSearches, $proteomeDiscovererVersion, FILE, $boundary
	my ($refModifications,$refMasses)=@_;
	print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>-Parameters...";

	my ($file,$db,$tol,$tolu,$itol,$itolu,$cle,$pfa,$decoy,$taxonomy,$instrument)=('','','','','','','','','','','');

	#my ($rawfilenumber) = $dbsqlite->selectrow_array("SELECT ProcessingNodeNumber FROM ProcessingNodes WHERE FriendlyName='Spectrum Files'");
	my $workFlowTable='WorkflowInfo';
	if ($proteomeDiscovererVersion >= 2.0) {
		$workFlowTable='Workflows';
		my ($wfXML) = $dbsqlite->selectrow_array("SELECT WorkflowXML FROM Workflows");
		my $xml = new XML::Simple();
		my $xmlData = $xml->XMLin($wfXML);
		foreach my $processingNodeParameters (sort{$a->{ParentProcessingNodeNumber} cmp $b->{ParentProcessingNodeNumber}} @{$xmlData->{WorkflowTree}->{WorkflowNode}}) {
			#next unless ($processingNodeParameters->{ProcessingNodeName} eq 'SpectrumFilesNode'  || $processingNodeParameters->{ProcessingNodeNumber} == $processingNodeNumber);
			if ($processingNodeParameters->{ProcessingNodeName} eq 'SpectrumFilesNode') {
				my @msfFiles=split(/;/,$processingNodeParameters->{ProcessingNodeParameters}->{ProcessingNodeParameter}->{content});
				my $addQuery=($fileID)?" AND FileID=$fileID":'';
				foreach my $mFile (@msfFiles){
					my @name=split(/\\/,$mFile);
					my $fName=pop(@name);
					my $fPath=join("\\",@name);
					$file.="File Name: $fName\n";
					$file.="File Path: \\$fPath\n";
				}
			}
			next unless ($processingNodeParameters->{ProcessingNodeNumber} == $processingNodeNumber);
			foreach my $processingNodeParameter (@{$processingNodeParameters->{ProcessingNodeParameters}->{ProcessingNodeParameter}}){
				if ($processingNodeParameter->{Name} eq 'Database'){
					$db=$processingNodeParameter->{content};
				}
				elsif ($processingNodeParameter->{Name} eq 'FastaDatabase'&& uc($refSearches->{$processingNodeNumber}{'DESCRIPTOR'})=~/SEQUEST/){
					$db=$processingNodeParameter->{content};
				}
				elsif($processingNodeParameter->{Name} eq 'PeptideTolerance'){
					($tol,$tolu)=split(/\s/,$processingNodeParameter->{content});
				}
				elsif($processingNodeParameter->{Name} eq 'FragmentTolerance'){
					($itol,$itolu)=split(/\s/,$processingNodeParameter->{content})
				}
				elsif($processingNodeParameter->{IntendedPurpose} eq 'CleavageReagent'){
					my $enzymeXML=new XML::Simple();
					my $enzymeDesc=$enzymeXML->XMLin($processingNodeParameter->{content}); # <Enzyme Version="1" Name="Trypsin" CleavageSites="KR" CleavageInhibitors="P" Offset="1" CleavageSpecificity="SpecificAtBothEnds" />
					$cle=$enzymeDesc->{Name};
				}
				elsif ($processingNodeParameter->{Name} eq 'Taxonomy'){
					$taxonomy=$processingNodeParameter->{content};
				}
				elsif($processingNodeParameter->{Name} =~ /(Max)?MissedCleavages/){
					$pfa=$processingNodeParameter->{content};
				}
				elsif($processingNodeParameter->{Name} eq 'UseDecoyDatabase' || $processingNodeParameter->{Name} eq 'DecoySearch'){
					$decoy=($processingNodeParameter->{IsValueSet} eq 'True')? 1 : 0;
				}
				elsif($processingNodeParameter->{Name} eq 'Instrument'){
					$instrument=$processingNodeParameter->{content};
				}
			}
			$refSearches->{$processingNodeNumber}{'DB_NAME'}=$db;
		}
	}
	else{
		my $sthParam1 = $dbsqlite->prepare("SELECT ParameterName, ValueDisplayString, ParameterValue FROM ProcessingNodeParameters WHERE ProcessingNodeNumber=$processingNodeNumber OR ProcessingNodeNumber=(SELECT ProcessingNodeNumber FROM ProcessingNodes WHERE FriendlyName='Spectrum Files')");
		$sthParam1->execute;
	    # In Sequest searches, true is not the taxonomy information
		while (my ($paramN,$valueDisplayStg,$paramV) = $sthParam1->fetchrow_array) {#Get all the parameters from the file
			if ($paramN eq 'SpectrumFileNames') {# Raw original filename(s) - if it is a merge, the filenames are separated by ';'
				#$file=$paramV; # changed on 24/07/12
				my @msfFiles=split(/;/,$valueDisplayStg);
				my $addQuery=($fileID)?" AND FileID=$fileID":'';
				foreach my $mFile (@msfFiles){
					my ($fTime)=$dbsqlite->selectrow_array("SELECT FileTime FROM FileInfos WHERE FileName='$mFile'$addQuery");
					next unless $fTime;
					my @name=split(/\\/,$mFile);
					my $fName=pop(@name);
					my $fPath=join("\\",@name);
					$file.="File Name: $fName\n";
					$file.="File Path: \\$fPath\n";
					$file.="File Time: $fTime\n\n";
				}
			}
			elsif ($paramN eq 'Database') {
				$db=$valueDisplayStg;
			}
			elsif ($paramN eq 'FastaDatabase' && uc($refSearches->{$processingNodeNumber}{'DESCRIPTOR'})=~/SEQUEST/) {
				$db=$valueDisplayStg;
				$refSearches->{$processingNodeNumber}{DB_ID}=$paramV;
			}
			elsif ($paramN eq 'PeptideTolerance') {
				($tol,$tolu)=split(/\s/,$valueDisplayStg);
			}
			elsif ($paramN eq 'FragmentTolerance') {
				($itol,$itolu)=split(/\s/,$valueDisplayStg);
			}
			elsif ($paramN eq 'Enzyme') {
				$cle=$valueDisplayStg;
			}
			elsif ($paramN eq 'Taxonomy') {
				$taxonomy=$valueDisplayStg;
			}
			elsif ($paramN eq 'MissedCleavages' || $paramN eq 'MaxMissedCleavages') {
				$pfa=$valueDisplayStg;
			}
			elsif ($paramN eq 'UseDecoyDatabase' || $paramN eq 'DecoySearch') {
				$decoy=($valueDisplayStg eq 'True')? 1 : 0;
			}
			elsif ($paramN eq 'Instrument') {
				$instrument=$valueDisplayStg;
			}
		}
		$sthParam1->finish;
	}

	my ($isSequest,$minScorePep,$com); #to know if the dat that is created is generated from an MSF file
	if (uc($refSearches->{$processingNodeNumber}{'DESCRIPTOR'})=~/SEQUEST/){
		&getSequestMonoMass($refMasses);
		$isSequest="SEQUEST";
		if ($decoy) {$minScorePep=0;} # Assumes Qvality-like FDR
		elsif ($selMinScore=~/def/i) {$minScorePep=&promsConfig::getMinScore('SEQUEST.PDM');} # default
		else {$minScorePep=$selMinScore;}
		($com)=$dbsqlite->selectrow_array("SELECT WorkflowName FROM $workFlowTable WHERE rowid=1"); # get the name of the project Workflow
	}
	else{ # MASCOT
		&getMascotMonoMass($refMasses);
		$isSequest="MASCOT";
		if ($decoy) {$minScorePep=0;} # Assumes FDR-based filtering
		elsif ($selMinScore=~/def/i) {$minScorePep=&promsConfig::getMinScore('MASCOT.PDM');} # default
		else {$minScorePep=$selMinScore;}
		my ($message) = $dbsqlite->selectrow_array("SELECT Message FROM WorkflowMessages WHERE Message LIKE \"Received Mascot result file%\"");
		#my @result=split(/[\/\\]/,$message);# $message="Mascot result on server (filename=../data/20100803/F008722.dat)";
		#$com=$result[-1];# $com="F008722.dat)";
		#$com=substr($com,0,length($com)-1);# $com="F008722.dat";
		($com)=($message=~/([^\/\\]+)\)\Z/); # "Mascot result on server (filename=../data/20100803/F008722.dat) -> F008722.dat
	}

	#########################################################################################################
	###>Get all the lines corresponding to modifications (static and dynamic) in ProcessingNodeParameters<###
	#########################################################################################################
	#my $sthParam2 = $dbsqlite->prepare("SELECT ParameterName,ValueDisplayString FROM ProcessingNodeParameters WHERE ProcessingNodeID=$refSearches->{$processingNodeNumber}{'NODEID'} AND ProcessingNodeNumber=$processingNodeNumber AND ParameterName LIKE '%Mod%' AND ParameterName NOT LIKE 'Search%'");
	#my $sthParam2 = $dbsqlite->prepare("SELECT ParameterName,ValueDisplayString FROM ProcessingNodeParameters WHERE ProcessingNodeID=$refSearches->{$processingNodeNumber}{'NODEID'} AND ProcessingNodeNumber=$processingNodeNumber AND ParameterName LIKE ? ESCAPE ? AND ParameterName NOT LIKE ?");
	#$sthParam2->execute('%Mod\_%',"\\",'Search%');
	my $sthDelMass = $dbsqlite->prepare("SELECT DeltaMass,AminoAcidModificationID FROM AminoAcidModifications WHERE ModificationName=? LIMIT 0,1"); # multiple rows are possible
	my $sthParam2;

	my $sthNLInfo;
	if ($proteomeDiscovererVersion >= 1.3){
		$sthNLInfo=$dbsqlite->prepare("SELECT Name,AANL.NeutralLossID,AANL.MonoisotopicMass FROM AminoAcids AA, AminoAcidModificationsNeutralLosses AANL, AminoAcidModificationsAminoAcidsNL AAMAANL WHERE AA.AminoAcidID=AAMAANL.AminoAcidID AND AAMAANL.NeutralLossID=AANL.NeutralLossID AND AminoAcidModificationID=? AND AA.OneLetterCode=? ORDER BY AANL.NeutralLossID ASC");
		#$sthParam2=$dbsqlite->prepare("SELECT ParameterName,ValueDisplayString FROM ProcessingNodeParameters WHERE ProcessingNodeID=$refSearches->{$processingNodeNumber}{NODEID} AND ProcessingNodeNumber=$processingNodeNumber AND ParameterName LIKE ? AND ValueDisplayString LIKE ?");
		# Retrieve Static_1 but not Static_X that is not a real fix-modification
		if ($proteomeDiscovererVersion < 2.0){
			$sthParam2=$dbsqlite->prepare("SELECT ParameterName,ValueDisplayString FROM ProcessingNodeParameters WHERE ProcessingNodeID=$refSearches->{$processingNodeNumber}{NODEID} AND ProcessingNodeNumber=$processingNodeNumber AND (ParameterName LIKE ? OR (ParameterName LIKE '%Static#_%' ESCAPE '#' AND ParameterName NOT LIKE '%Static#_X' ESCAPE '#')) AND ValueDisplayString LIKE ?");
			$sthParam2->execute('%Mod%','%+%Da%');
		}
	}
	elsif ($proteomeDiscovererVersion == 1.2) {
		$sthParam2=$dbsqlite->prepare("SELECT ParameterName,ValueDisplayString FROM ProcessingNodeParameters WHERE ProcessingNodeID=$refSearches->{$processingNodeNumber}{NODEID} AND ProcessingNodeNumber=$processingNodeNumber AND ParameterName LIKE ? AND ParameterName NOT LIKE ?");
		$sthParam2->execute('%Mod%','Search%');
	}
	else{
		$sthParam2=$dbsqlite->prepare("SELECT ParameterName,ValueDisplayString FROM ProcessingNodeParameters WHERE ProcessingNodeID=$refSearches->{$processingNodeNumber}{NODEID} AND ProcessingNodeNumber=$processingNodeNumber AND ParameterName LIKE ?");
		$sthParam2->execute('%Mod%');
	}
	#Fill the static and dynamic modification informations
	my ($mods,$itmods,$modsString,$itmodsString)=('','','','');
	my ($i,$j)=(1,1);
	if ($proteomeDiscovererVersion >= 2.0){
		my ($wfXML) = $dbsqlite->selectrow_array("SELECT WorkflowXML FROM Workflows");
		my $xml = new XML::Simple();
		my $xmlData = $xml->XMLin($wfXML);
		foreach my $processingNodeParameters (@{$xmlData->{WorkflowTree}->{WorkflowNode}}) {
			next unless $processingNodeParameters->{ProcessingNodeNumber} == $processingNodeNumber;
			foreach my $processingNodeParameter (@{$processingNodeParameters->{ProcessingNodeParameters}->{ProcessingNodeParameter}}){
				if(($processingNodeParameter->{IntendedPurpose} =~ /Static(Terminal)?Modification/) && $processingNodeParameter->{IsValueSet} eq "True"){
					my $statmodXML=new XML::Simple();
					my $statmodDesc=$statmodXML->XMLin($processingNodeParameter->{content}); #<Modification Version="2" AminoAcids="X" Name="MappingL" Abbreviation="L" ID="1167" UnimodAccession="-1" DeltaMass="113.08406" DeltaAverageMass="113.15890" IsSubstitution="True" LeavingGroup="" Substitution="C6H11NO" PositionType="Any" />
					next if $statmodDesc->{Name} eq "MascotXValue" || $statmodDesc->{Name} eq "MappingL";
					my $delta=$statmodDesc->{DeltaMass};
					my $res=($statmodDesc->{AminoAcids})? $statmodDesc->{AminoAcids} : $statmodDesc->{Terminus};
					if ($statmodDesc->{Terminus}) { # converts "Any_N_Terminus" to standard "Any N-Terminus" (PP 01/07/16)
						$res=~s/Any_/Any /i;
						$res=~s/_Term/-Term/i;
					}
					my $valueD="$statmodDesc->{Name} ($res)";
					$mods.=",$valueD";
					$modsString.="FixedMod$i=$delta,$valueD\nFixedModResidues$i=$res\n";
					$i++;
				}
				elsif(($processingNodeParameter->{IntendedPurpose} =~ /Dynamic(Terminal)?Modification/) && $processingNodeParameter->{IsValueSet} eq "True"){
					my $dynmodXML=new XML::Simple();
					my $dynmodDesc=$dynmodXML->XMLin($processingNodeParameter->{content}); #<Modification Version="2" AminoAcids="C" Name="Carbamidomethyl" Abbreviation="Carbamidomethyl" ID="8" UnimodAccession="4" DeltaMass="57.02146" DeltaAverageMass="57.05130" IsSubstitution="False" LeavingGroup="" Substitution="H(3) C(2) N O" PositionType="Any" />
					my $res=($dynmodDesc->{AminoAcids})? $dynmodDesc->{AminoAcids} : $dynmodDesc->{Terminus};
					if ($dynmodDesc->{Terminus}) { # converts "Any_N_Terminus" to standard "Any N-Terminus" (PP 01/07/16)
						$res=~s/Any_/Any /i;
						$res=~s/_Term/-Term/i;
					}
					my $valueD="$dynmodDesc->{Name} ($res)";
					my $paramN=$processingNodeParameter->{Name};
					$refModifications->{$paramN}{'NAME'}=$valueD;
					$refModifications->{$paramN}{'SUBST'}=$dynmodDesc->{Name}; # for terminal modification
					$refModifications->{$paramN}{'DELTA'}=$dynmodDesc->{DeltaMass};
					$refModifications->{$paramN}{'ID'}=$dynmodDesc->{ID};
					$refModifications->{$paramN}{'VALUE'}=$j;
					$valueD=~ s/, //g;
					$itmods.=",$valueD";
					$itmodsString.="delta$j=$dynmodDesc->{DeltaMass},$valueD\n";
					my $addNLString='';
					if ($dynmodDesc->{AminoAcids}) {
						my %nl;# nl <-> Neutral-Losses
						foreach my $aa (split(//,$dynmodDesc->{AminoAcids})) {
							$sthNLInfo->execute($dynmodDesc->{ID},$aa);
							while (my ($nlName,$nlID,$nlMass)=$sthNLInfo->fetchrow_array){
								$nl{$nlID}=$nlMass;
							}
						}
						foreach my $nlID (sort{$a<=>$b} keys %nl) {# Print neutral-loss like Mascot (v2.3) DAT files.
							if ($addNLString) {
								$addNLString.="NeutralLoss${j}_master=$nl{$nlID}\n";
							}else{
								$addNLString.="NeutralLoss$j=$nl{$nlID}\n";
							}
						}
						###> For Phospho(Y) var-mod, must write specifically this line.
						if (!$addNLString && $dynmodDesc->{Name} =~ /Phospho/ && $dynmodDesc->{AminoAcids} eq 'Y') {
							$addNLString.="NeutralLoss$j=0.000000\n";
						}#code
					}
					$itmodsString.=$addNLString;
					$refModifications->{$paramN}{'PRINTED_NAME'}=$valueD;
					$j++;
				}
			}
		}
	}
	else{
		while (my ($paramN,$valueD) = $sthParam2->fetchrow_array) {
			next unless $valueD;
			my $oldValueD=$valueD;
			# Substitution to transform the new ProteomeDiscover value Oxidation / +15.995 Da (M)
			# into a more readable value Oxidation (M)
			#$valueD =~ s/[0-9]*//g;
			#$valueD =~ s/\.//g;
			#$valueD =~ s/ \///g;
			#$valueD =~ s/\+//g;
			#$valueD =~ s/ Da //g;
			$valueD=~s/\/.+\(/\(/;
			#my @modif= split(/\(/,$valueD);
			#$modif[0]=~ s/ //g;
			#$refModifications->{$paramN}{'SUBST'}=$modif[0];
			(my $subValueD=$valueD)=~s/ \([^\(]+\)\Z//; # Oxidation (M) -> Oxidation
			#my $sth2 = $dbsqlite->prepare("SELECT DeltaMass FROM AminoAcidModifications WHERE ModificationName LIKE \"$modif[0]\";");
			$sthDelMass->execute($subValueD);
			my ($delta,$aamodifID)=$sthDelMass->fetchrow_array;
			if ($paramN =~ /^Stat/) { # fix mod
				$valueD=~ s/, //g;
				$mods.=",$valueD";
				#my @res= split(/\(/,$valueD);
				#$res[1]=~ s/\)//g;
				my ($res)=$valueD=~/\(([^\(]+)\)\Z/; # Acetyl (Any N-Terminus) -> Any N-Terminus
				$modsString.="FixedMod$i=$delta,$valueD\nFixedModResidues$i=$res";
				$i++;
			}
			elsif ($paramN =~ /^Dyn/) { # var mod
				$refModifications->{$paramN}{'NAME'}=$valueD;
				$refModifications->{$paramN}{'SUBST'}=$subValueD; # for terminal modification
				$refModifications->{$paramN}{'DELTA'}=$delta;
				$refModifications->{$paramN}{'ID'}=$aamodifID;
				$refModifications->{$paramN}{'VALUE'}=$j;
				$valueD=~ s/, //g;
				$itmods.=",$valueD";
				$itmodsString.="delta$j=$delta,$valueD\n";
				if ($proteomeDiscovererVersion >= 1.3){# Version has to be 1.3 at least to get the NL information
					my $addNLString='';
					my ($aas)=($oldValueD =~ /\((\w+)\)/)?$1:'';# There might be more than one AA
					if($aas){
						my %nl;# nl <-> Neutral-Losses
						foreach my $aa (split(//,$aas)) {
							$sthNLInfo->execute($aamodifID,$aa);
							while (my ($nlName,$nlID,$nlMass)=$sthNLInfo->fetchrow_array){
								$nl{$nlID}=$nlMass;
							}
						}
						foreach my $nlID (sort{$a<=>$b} keys %nl) {# Print neutral-loss like Mascot (v2.3) DAT files.
							if ($addNLString) {
								$addNLString.="NeutralLoss${j}_master=$nl{$nlID}\n";
							}else{
								$addNLString.="NeutralLoss$j=$nl{$nlID}\n";
							}
						}
						###> For Phospho(Y) var-mod, must write specifically this line.
						if (!$addNLString && $subValueD =~ /Phospho/ && $aas eq 'Y') {
							$addNLString.="NeutralLoss$j=0.000000\n";
						}
					}
					$itmodsString.=$addNLString;
				}
				$refModifications->{$paramN}{'PRINTED_NAME'}=$valueD;# value printed in the dat file
				$j++;
			}
			else{
				# Ecrire un fichier error_log pour signifier que le cas n'a pas été retenu
			}
		}
		$sthParam2->finish;
		$sthDelMass->finish;
	}

	if ($proteomeDiscovererVersion >= 1.3){$sthNLInfo->finish;}
	$mods = ($mods) ? (substr $mods, 1) : ''; # remove starting ','
	$itmods = ($itmods) ? (substr $itmods, 1) : ''; # remove starting ','

	##############################################################
	###>Find Quantification if any in ProcessingNodeParameters<###
	##############################################################
	my ($quantifName,$paramValue,$processMatched,$processCount)=('','',0,0);
	if ($proteomeDiscovererVersion >= 2.0 ) {
		my ($wfXML) = $dbsqlite->selectrow_array("SELECT WorkflowXML FROM Workflows");
		my $xml = new XML::Simple();
		my $xmlData = $xml->XMLin($wfXML);
		### PrecursorIonsQuantifierNode -> SILAC
		### ReporterIonQuantifierNode -> TMT
		foreach my $processingNodeParameters (@{$xmlData->{WorkflowTree}->{WorkflowNode}}) {
			next unless ($processingNodeParameters->{ProcessingNodeName} eq 'PrecursorIonsQuantifierNode' || $processingNodeParameters->{ProcessingNodeName} eq 'ReporterIonQuantifierNode');
			foreach my $processingNodeParameter (@{$processingNodeParameters->{ProcessingNodeParameters}->{ProcessingNodeParameter}}){
				if($processingNodeParameter->{IntendedPurpose} eq "QuantificationMethod"){
					$paramValue=$processingNodeParameter->{content};
					($quantifName)=($paramValue=~/^<ProcessingMethod name=\"(.*)\" version/);
				}
			}
			if ($paramValue) {
				my $parentNodeStrg=$processingNodeParameters->{ParentProcessingNodeNumber};
				$processMatched=0;
				foreach my $parNodeNum (split(/;/,$parentNodeStrg)) {
					if ($parNodeNum==$processingNodeNumber) {
						$processMatched=1;
						last;
					}
				}
			}
		}
	}
	else{

		my $sthQN=$dbsqlite->prepare("SELECT ProcessingNodeParentNumber,ValueDisplayString,ParameterValue FROM ProcessingNodeParameters,ProcessingNodes WHERE ProcessingNodeParameters.ProcessingNodeID=ProcessingNodes.ProcessingNodeID AND ParameterName LIKE 'Quanti%ationMethod'"); # Quantification or Quantitation
		#$sthQN->execute;
		#QNODE:while ((my $parentNodeStrg,$quantifName,$paramValue)=$sthQN->fetchrow_array) {
		#	foreach my $parNodeNum (split(/;/,$parentNodeStrg)) {
		#		last QNODE if $parNodeNum==$processingNodeNumber;
		#	}
		#}
		#$sthQN->finish;
		#$quantifName='' unless $quantifName;# if there is no quantification in the MSF file, an error is printed in the HTML file: "warning: Use of uninitialized value $quantifName in concatenation (.)"

		my $sthPN=$dbsqlite->prepare("SELECT ProcessingNodeParentNumber FROM ProcessingNodes WHERE ProcessingNodeNumber=?");
		$sthQN->execute;
		while (my ($parentNodeStrg,$qName,$parValue)=$sthQN->fetchrow_array) {
			$processCount++;
			$quantifName=$qName;
			$paramValue=$parValue;
			$processMatched=&checkSearchNode($sthPN,$processingNodeNumber,$parentNodeStrg);
			last if $processMatched;
		}
		$sthQN->finish;
		$sthPN->finish;
	}

	if ($processMatched==0 && $processCount > 1) { # multiple quantifs that cannot be linked to the $processingNodeNumber used => skip quantif info
		$quantifName=$paramValue='';
		print "[*Unable to match quantification process to search process: Skipping quantification info*]";
	}
	my $isPercolator=($refSearches->{$processingNodeNumber}{PERCOLATOR_NODE})? 'Percolator' : '';
	my $isPhosphoRS=($refSearches->{$processingNodeNumber}{PHOSPHORS_NODE})? 1 : 0;

	#######################
	###>Writing to file<###
	#######################
#LICENSE=Licensed to: Institut Curie, Lab. de Spectrom. de Mass Proteomique (No. 2G-49264 Loew), (4 processors).
	print FILE qq
|MIME-Version: 1.0 (Generated by Mascot version 1.0)
Content-Type: multipart/mixed; boundary=gc0p4Jq0M2Yt08jU534c0p

$boundary
Content-Type: application/x-Mascot; name="parameters"

LICENSE=
MP=
NM=
COM=$com
IATOL=
IA2TOL=
IASTOL=
IBTOL=
IB2TOL=
IBSTOL=
IYTOL=
IY2TOL=
IYSTOL=
SEG=
SEGT=
SEGTU=
LTOL=
TOL=$tol
TOLU=$tolu
ITH=
ITOL=$itol
ITOLU=$itolu
PFA=$pfa
DB=$db
MODS=$mods
MASS=Monoisotopic
CLE=$cle
FILE=$file
PEAK=
QUE=
TWO=
SEARCH=MIS
USERNAME=
USEREMAIL=
CHARGE=2+ and 3+
INTERMEDIATE=
REPORT=AUTO
OVERVIEW=
FORMAT=Mascot generic
FORMVER=1.01
FRAG=
IT_MODS=$itmods
USER00=
USER01=
USER02=
USER03=
USER04=
USER05=
USER06=
USER07=
USER08=
USER09=
USER10=
USER11=
USER12=
PRECURSOR=
TAXONOMY=$taxonomy
ACCESSION=
REPTYPE=
SUBCLUSTER=
ICAT=
INSTRUMENT=$instrument
ERRORTOLERANT=
FRAMES=
CUTOUT=
USERID=0
QUANTITATION=$quantifName
DECOY=$decoy
SEARCH_ALGO=$isSequest
FDR_ALGO=$isPercolator
PHOSPHORS=$isPhosphoRS
PEP_ISOTOPE_ERROR=
RULES=1,2,8,9,10,13,14,15
INTERNALS=0.0,700.0
$boundary
Content-Type: application/x-Mascot; name="masses"

A=$refMasses->{'A'}
B=$refMasses->{'B'}
C=$refMasses->{'C'}
D=$refMasses->{'D'}
E=$refMasses->{'E'}
F=$refMasses->{'F'}
G=$refMasses->{'G'}
H=$refMasses->{'H'}
I=$refMasses->{'I'}
J=$refMasses->{'J'}
K=$refMasses->{'K'}
L=$refMasses->{'L'}
M=$refMasses->{'M'}
N=$refMasses->{'N'}
O=$refMasses->{'O'}
P=$refMasses->{'P'}
Q=$refMasses->{'Q'}
R=$refMasses->{'R'}
S=$refMasses->{'S'}
T=$refMasses->{'T'}
U=$refMasses->{'U'}
V=$refMasses->{'V'}
W=$refMasses->{'W'}
X=$refMasses->{'X'}
Y=$refMasses->{'Y'}
Z=$refMasses->{'Z'}
Hydrogen=$refMasses->{'Hydrogen'}
Carbon=$refMasses->{'Carbon'}
Nitrogen=$refMasses->{'Nitrogen'}
Oxygen=$refMasses->{'Oxygen'}
Electron=$refMasses->{'Electron'}
C_term=$refMasses->{'C_term'}
N_term=$refMasses->{'N_term'}
$itmodsString$modsString
$boundary
|;

	if ($paramValue) {# Print the XML information of the quantitation method used by PD for SILAC !
		###> Add a synonym to this modification...
		$paramValue=~s/""/"/g;
		my $xml = new XML::Simple();
		my $xmlData = $xml->XMLin($paramValue);
		my $labelXLM=$xmlData->{MethodPart}{QuanChannels}{MethodPart};

		my $dbh=&promsConfig::dbConnect;
		foreach my $labelName (keys %{$labelXLM}){
			if ($labelXLM->{$labelName}{MethodPart}{MethodPart}{MethodPart}) { # No label or single label
				my $labelModifAlias=$labelXLM->{$labelName}{MethodPart}{MethodPart}{name};
				if ($labelXLM->{$labelName}{MethodPart}{MethodPart}{MethodPart}{MethodPart}{SideChainModification}{Parameter}) {
					my $labelParamStrg=$labelXLM->{$labelName}{MethodPart}{MethodPart}{MethodPart}{MethodPart}{SideChainModification}{Parameter}{content};
					my ($labelModifName)=$labelParamStrg=~/ Name="([^"]+)/;
					my ($modifRes)=$labelParamStrg=~/ AminoAcids="([^"]+)/;
					my $modificationID=&promsMod::getModificationIDfromString($dbh,$labelModifName,$modifRes,\%vmodsInUnimod);
					#&addSynonyms($dbh,$modificationID,$labelModifAlias);
					#&addSynonyms($dbh,$modificationID,$quantifName);
				}
			}
			else {
				foreach my $labelModifAlias (sort{lc($a) cmp lc($b)} keys %{$labelXLM->{$labelName}{MethodPart}{MethodPart}}) {
					my $labelParamStrg=$labelXLM->{$labelName}{MethodPart}{MethodPart}{$labelModifAlias}{MethodPart}{MethodPart}{SideChainModification}{Parameter}{content};
					my ($labelModifName)=$labelParamStrg=~/ Name="([^"]+)/;
					my ($modifRes)=$labelParamStrg=~/ AminoAcids="([^"]+)/;
					my $modificationID=&promsMod::getModificationIDfromString($dbh,$labelModifName,$modifRes,\%vmodsInUnimod);
					#&addSynonyms($dbh,$modificationID,$labelModifAlias);
					#&addSynonyms($dbh,$modificationID,$quantifName);
				}
			}
		}
		$dbh->disconnect;

		print FILE qq
|Content-Type: application/x-Mascot; name="quantitation"

$paramValue
$boundary
|;
	}
	print " Done.</B><BR>\n";
	return ($db,$decoy,$minScorePep);
}

##########################################
####<Print header section (.pdm file)>####
##########################################
sub printHeaderMSF { # GLOBALS: $dbsqlite, $processingNodeNumber, FILE
	print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>-Header...";
    my ($db,$nbquery)=@_;
    my ($wftable)=($proteomeDiscovererVersion < 2.0)?'WorkflowInfo':'Workflows';
    my ($date)=$dbsqlite->selectrow_array("SELECT WorkflowStartDate FROM $wftable");
    my ($dateCode,$startTime)=split(/\s/,$date);
	$startTime=~s/\.\d+//; # 12:34:05.625 -> 12:34:05
    my ($dbpath)=($refSearches->{$processingNodeNumber}{DB_ID})? $dbsqlite->selectrow_array("SELECT FileName FROM FastaFiles WHERE FastaFileID=$refSearches->{$processingNodeNumber}{DB_ID}"): ($refSearches->{$processingNodeNumber}{'DB_NAME'})? $refSearches->{$processingNodeNumber}{'DB_NAME'} : '';
    my ($friendlyName)=$refSearches->{$processingNodeNumber}{DESCRIPTOR};
    #my ($friendlyName)=$dbsqlite->selectrow_array("SELECT FriendlyName FROM ProcessingNodes WHERE ProcessingNodeNumber=$processingNodeNumber");

    my ($dbName,$mascotVersion,$dbNbSeq,$dbNbSeqAT,$nbQueries)=("","","","","");
    if (uc($friendlyName)=~/MASCOT/ && $proteomeDiscovererVersion >= 1.3){# In PD 1.3+, some other information is available in CustomDataProcessingNodes
		my ($dbInfo)=$dbsqlite->selectrow_array("SELECT FieldValue FROM CustomDataProcessingNodes WHERE ProcessingNodeNumber=$processingNodeNumber AND FieldID=(SELECT FieldID FROM CustomDataFIelds WHERE DisplayName='Fasta database information')");
		foreach my $field (split("\n",$dbInfo)) {
			my ($title,$value)=split(": ",$field);
			if($title eq 'Version'){ $mascotVersion=$value;}
			elsif($title eq 'Number of sequences'){ $dbNbSeq=$value;}
			elsif($title eq 'Number of sequences after taxonomy'){ $dbNbSeqAT=$value;}
			elsif($title eq 'FASTA db'){ $db=$value;}
			#elsif($title eq 'Number of queries'){ $nbquery=$value;}
		}
    }

	print FILE qq
|Content-Type: application/x-Mascot; name="header"

sequences=$dbNbSeq
sequences_after_tax=$dbNbSeqAT
residues=
distribution=
exec_time=
date=
time=$startTime
queries=$nbquery
max_hits=
version=$mascotVersion
fastafile=$dbpath
release=$db
taskid=
$boundary
|;
	print " Done.</B><BR>\n";
}

###################################################################
####<Print Summary: the queries one by one with the ions found>####
###################################################################
sub printSummaryMSF { # GLOBALS: $dbsqlite, FILE
	print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>-Summary...";
	my ($decoy,$refQueries,$refSpectra)=@_;

	my @summaryTypes=('summary');
	push @summaryTypes,'decoy_summary' if $decoy;

	foreach my $summary (@summaryTypes) { # Target (and Decoy)
		print FILE "Content-Type: application/x-Mascot; name=\"$summary\"\n\n";
		foreach my $queryNum (sort{$a<=>$b} keys %{$refQueries}) {
			my $spectrumID=$refQueries->{$queryNum};
			print FILE "qmass$queryNum=$refSpectra->{$spectrumID}[0]\n";
			print FILE "qexp$queryNum=$refSpectra->{$spectrumID}[2],$refSpectra->{$spectrumID}[1]+\n";
		}
		print FILE "num_hits=0\n";
		print FILE "$boundary\n";
		print '.';
	}
	print " Done.</B><BR>\n";
}

###################################################################
####<Print Summary: the queries one by one with the ions found>####
###################################################################
sub printPeptidesMSF { # GLOBALS: $pdmFile, $dbsqlite, $processingNodeNumber, $refSearches, FILE, $databankID, $userID
	my ($decoy,$minScorePep,$refQueries,$refSpectra,$refModifications,$refMasses)=@_;

	########################################
	####<Get protein info from MSF file>####
	########################################
#print '1>';
	print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>-Proteins...";
	my %proteinsInfo;
	my $dbh=&promsConfig::dbConnect;
	my ($parseRules,$identType,$defIdentType)=$dbh->selectrow_array("SELECT PARSE_RULES,IDENTIFIER_TYPE,DEF_IDENT_TYPE FROM DATABANK,DATABANK_TYPE WHERE DATABANK.ID_DBTYPE=DATABANK_TYPE.ID_DBTYPE AND ID_DATABANK=$databankID"); # $databankID global. Only 1 db-search allowed in MSF
	$defIdentType='UNKNOWN' unless $defIdentType;
	$dbh->disconnect;
	$identType=$defIdentType unless $identType;
	my @rules=split(',:,',$parseRules);
	my ($idRule)=($rules[0]=~/ID=(.+)/); #<<<<<<<<<<<<<< Synchronise with promsMod::getProtInfo >>>>>>
	my $tempDbFile;
	if ($fileID) {
		($tempDbFile ="$promsPath{tmp}/$msfFile")=~s/\.msf\Z/_$processingNodeNumber\.$fileID\.$userID.fasta/;
	}else{
		($tempDbFile ="$promsPath{tmp}/$msfFile")=~s/\.msf\Z/_$processingNodeNumber\.$userID.fasta/; # same dir as msf file
	}
	open (FASTA,">$tempDbFile");
#	print "";
	my $sthProt = ($proteomeDiscovererVersion >= 2.0) ? $dbsqlite->prepare("SELECT Proteins.ProteinID, Sequence, Description FROM Proteins, ProteinAnnotations WHERE Proteins.ProteinID=ProteinAnnotations.ProteinID") : $dbsqlite->prepare("SELECT Proteins.ProteinID, Sequence, Description FROM Proteins, ProteinAnnotations, FastaFilesProteinAnnotations WHERE Proteins.ProteinID=ProteinAnnotations.ProteinID AND ProteinAnnotations.ProteinAnnotationID=FastaFilesProteinAnnotations.ProteinAnnotationID AND FastaFileID=(SELECT ParameterValue FROM ProcessingNodeParameters WHERE ParameterName='FastaDatabase' AND ProcessingNodeNumber=$processingNodeNumber)");
	$sthProt->execute;
	while (my ($proteinID,$sequence,$description)= $sthProt->fetchrow_array) {
		print FASTA "$description\n$sequence\n";
		@{$proteinsInfo{$proteinID}}=(
			[], 		# 0 identifiers
			$sequence	# 1
		);
        $description=~s/^>//;
		foreach my $entry (split "\001",$description) { #\001 means SOH<-> Start of header
			my ($identifier)=($entry=~/$idRule/);
			if ($identifier && $identType eq 'UNIPROT_ID') { #check for isoforms in 1st keyword before |
				$entry=~s/^sp\|//;
				if ($entry=~/^[^|]+-(\d+)\|/) {
					$identifier.="-$1";
				}
			}
			else {
				($identifier)=($description=~/^(\S+)/) unless $identifier; # fall back to fasta full identifier
			}
			push @{$proteinsInfo{$proteinID}[0]},$identifier;
		}
	}
	$sthProt->finish;
	close FASTA;
	print " Done.</B><BR>\n";

	########################################################
	####<Get peptide info: Looping from target to decoy>####
	########################################################
	#>Get Scoring type used
	my $queryScore=($proteomeDiscovererVersion >= 2.0)? "SELECT MIN(ScoreID) FROM ProcessingNodeScores WHERE NodeGuid='$refSearches->{$processingNodeNumber}{NODEGUID}'" : "SELECT MIN(ScoreID) FROM ProcessingNodeScores WHERE ProcessingNodeID=$refSearches->{$processingNodeNumber}{NODEID}";
	my ($scoreID) = $dbsqlite->selectrow_array($queryScore);
	my $activeModStrg=(uc($refSearches->{$processingNodeNumber}{'DESCRIPTOR'})=~/SEQUEST/)? 'AND isActive=1' : ''; # for proper modification selection
	#my ($proteomeDiscovererVersion)=$dbsqlite->selectrow_array("SELECT SoftwareVersion FROM SchemaInfo");
	#$proteomeDiscovererVersion=~s/^(\d\.\d+).*/$1/; # x.x.x.xx -> x.x (numerical value)
	my $addMissCutInfo=($proteomeDiscovererVersion >= 1.3)? 'MissedCleavages' : "''"; # Since PD 1.3, MissedCleavage information is found in Peptides and PeptidesDecoy

	my @peptideTypes=(['Target','','peptides']);
	push @peptideTypes,['Decoy','_decoy','decoy_peptides'] if $decoy;


	my %percolatorMeasures;
	my $displayTag=($proteomeDiscovererVersion >= 2.0)?'Percolator ':'';
	if ($percolatorThres && $refSearches->{$processingNodeNumber}{'PERCOLATOR_NODE'}) {
		#my $sthPercol=$dbsqlite->prepare("SELECT DisplayName,FieldID FROM CustomDataFields WHERE SourceNodeNumber=$refSearches->{$processingNodeNumber}{PERCOLATOR_NODE} AND (DisplayName='q-Value' OR DisplayName='PEP')");
		my $sthPercol=$dbsqlite->prepare("SELECT DisplayName,FieldID FROM CustomDataFields WHERE SourceNodeNumber=$refSearches->{$processingNodeNumber}{PERCOLATOR_NODE} AND TargetNodeNumber=$processingNodeNumber AND (DisplayName='${displayTag}q-Value' OR DisplayName='${displayTag}PEP')");
		$sthPercol->execute;
		while (my($displayName,$fieldID)=$sthPercol->fetchrow_array) {
			push @{$percolatorMeasures{$displayName}},$fieldID; # 2 entries: 1st for target peptides and 2nd for decoy
		}
		$sthPercol->finish;
	}

	my %phosphoRSMeasures;
	if ($refSearches->{$processingNodeNumber}{'PHOSPHORS_NODE'}) {
		my $sthPRScol=$dbsqlite->prepare("SELECT DisplayName,FieldID FROM CustomDataFields WHERE SourceNodeNumber=$refSearches->{$processingNodeNumber}{PHOSPHORS_NODE} AND TargetNodeNumber=$processingNodeNumber AND (DisplayName LIKE 'phosphoRS%')");
		# 3 Informations provided in CustomDataFields :
		# - phosphoRS Binomial Peptide Score
		# - phosphoRS Isoform Probability
		# - phosphoRS Site Probabilities
		$sthPRScol->execute;
		while (my($displayName,$fieldID)=$sthPRScol->fetchrow_array) {
			push @{$phosphoRSMeasures{"$displayName"}},$fieldID; # 3 entries: 1st Binomial Peptide Score, 2nd Isoform Probability and 3rd Site Probabilities
		}
		$sthPRScol->finish;
	}

	foreach my $refType (@peptideTypes) { # Target (and Decoy)
		my ($typeName,$tableFlag,$sectionName)=@{$refType};

		print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>-$typeName peptides...";
		my (%peptidesInfo,%spectrumPeptides); #,%peptidesDecoyInfo,%spectrumDecoyPeptides

		######################################################################
		###<Get all the information from Peptides table in SQLite datafile>###
		######################################################################
#print "2 $typeName>";
		###<Percolator filtering>###
		my $sthPepInfo;
		#if ($percolatorThres && $refSearches->{$processingNodeNumber}{'PERCOLATOR_NODE'}) {
		my ($qValID,$PEP_ID)=($typeName eq 'Target')? ($percolatorMeasures{"${displayTag}q-Value"}[0],$percolatorMeasures{"${displayTag}PEP"}[0]) : ($percolatorMeasures{"${displayTag}q-Value"}[1],$percolatorMeasures{"${displayTag}PEP"}[1]);
		if ($percolatorThres && defined($qValID) && defined($PEP_ID)) {
			$sthPepInfo = $dbsqlite->prepare("SELECT P.PeptideID,ScoreValue,SpectrumID,MatchedIonsCount,Sequence,$addMissCutInfo,GROUP_CONCAT(FieldValue,':') FROM PeptideScores$tableFlag PS,Peptides$tableFlag P,CustomDataPeptides$tableFlag CD WHERE P.ProcessingNodeNumber=$processingNodeNumber AND PS.PeptideID=P.PeptideID AND P.PeptideID=CD.PeptideID AND ScoreID=$scoreID AND FieldID IN ($qValID,$PEP_ID) GROUP BY CD.PeptideID");
		}
		else {
			$sthPepInfo = $dbsqlite->prepare("SELECT P.PeptideID,ScoreValue,SpectrumID,MatchedIonsCount,Sequence,$addMissCutInfo FROM PeptideScores$tableFlag PS,Peptides$tableFlag P WHERE P.ProcessingNodeNumber=$processingNodeNumber AND PS.PeptideID=P.PeptideID AND ScoreID=$scoreID AND ScoreValue>$minScorePep");
		}
		#my $step=5000; #1000; # int($numPeptides/10);
		my $count=0;
		$sthPepInfo->execute;
		while (my($peptideID,$score,$spectrumID,$matchedIonsCount,$sequence,$missCut,$percolData) = $sthPepInfo->fetchrow_array) {
			$count++;
			if ($count>=10000) { # 2500
				$count=0;
				print '.';
			}
			next unless $spectrumInfo{$spectrumID};
			if ($percolData) {
				my ($qVal,$PEP)=split(':',$percolData);
				next if $qVal > $percolatorThres;
			}

			@{$peptidesInfo{$peptideID}}=(
				$score,							# 0
				$sequence,						# 1
				$matchedIonsCount,				# 2
				0.000000,						# 3
				$missCut,						# 4 only defined if $proteomeDiscovererVersion >= 1.3;
				"0" x (length($sequence)+2),	# 5 String that explain the modifications on the peptides
				[]								# 6 matching proteins
			);
			push @{$peptidesInfo{$peptideID}},$percolData if $percolData; # 7 Percolator data
			push @{$spectrumPeptides{$spectrumID}},$peptideID; # 1 spectrum <-> several peptides
		}
		$sthPepInfo->finish;
		print '/';

#print "3 $typeName>";
		if ($typeName eq 'Target' || $proteomeDiscovererVersion >= 1.3) { # Before v1.3: no PeptidesProteins_decoy table
			my $sthPepProt = $dbsqlite->prepare("SELECT PeptideID, ProteinID FROM PeptidesProteins$tableFlag WHERE ProcessingNodeNumber=$processingNodeNumber;");
			$sthPepProt->execute;
			while( my($peptideID,$proteinID) = $sthPepProt->fetchrow_array){
				next unless $peptidesInfo{$peptideID};
				push @{$peptidesInfo{$peptideID}[6]},$proteinID;
			}
			$sthPepProt->finish;
			print '.';
		}

		#######################################################################
		###<Update the mass value of the peptide according to modifications>###
		#######################################################################
		##>Regular modifications
#print "4 $typeName>";
		#my $sthRegMod = $dbsqlite->prepare("SELECT Position, PeptideID, DeltaMass, ModificationName FROM PeptidesAminoAcidModifications, AminoAcidModifications WHERE AminoAcidModifications.AminoAcidModificationID=PeptidesAminoAcidModifications.AminoAcidModificationID AND PeptidesAminoAcidModifications.ProcessingNodeNumber=$processingNodeNumber AND isActive=1 ORDER BY PeptideID, Position ASC");#Normal modifications
		my $queryModif="SELECT Position, PeptideID, DeltaMass, ModificationName FROM PeptidesAminoAcidModifications$tableFlag, AminoAcidModifications WHERE AminoAcidModifications.AminoAcidModificationID=PeptidesAminoAcidModifications$tableFlag.AminoAcidModificationID AND PeptidesAminoAcidModifications$tableFlag.ProcessingNodeNumber=$processingNodeNumber $activeModStrg ORDER BY PeptideID, Position ASC";
		my $sthRegMod = $dbsqlite->prepare($queryModif);#Normal modifications
		$sthRegMod->execute;
		$count=0;
		while (my($position,$peptideID,$deltaMass,$modificationName) = $sthRegMod->fetchrow_array) { #,$isActive
			$count++;
			if ($count>=20000) {
				$count=0;
				print '.';
			}
			#if($isActive == 1){
				next unless $peptidesInfo{$peptideID};
				$modificationName = quotemeta(uc($modificationName)); # quote because match below
				$peptidesInfo{$peptideID}[3]+=$deltaMass;
				my @sequence = split(//,$peptidesInfo{$peptideID}[1]);
				foreach my $modif (keys %{$refModifications}) {
					#if (uc($refModifications->{$modif}{'NAME'}) =~ /$modificationName \(\w*$sequence[$position]\w*\)/){
					#if (uc($refModifications->{$modif}{'PRINTED_NAME'}) =~ /^$modificationName \(\w*$sequence[$position]\w*\)/){
					if (uc($refModifications->{$modif}{'PRINTED_NAME'}) =~ /^$modificationName \(\w*$sequence[$position]\w*\)/){ # if ' is forgotten, then, for Methyl, Dimethyl, Trimethyl modification, it does not find the correct one!
						substr($peptidesInfo{$peptideID}[5] , $position+1 , 1 , $refModifications->{$modif}{'VALUE'} );
					}
				}
			#}
		}
		$sthRegMod->finish;

		##>Terminal modifications
#print "5 $typeName>";
		#my $sthTermMod = $dbsqlite->prepare("SELECT PeptideID, DeltaMass, ModificationName FROM PeptidesTerminalModifications, AminoAcidModifications WHERE PeptidesTerminalModifications.TerminalModificationID=AminoAcidModifications.AminoAcidModificationID AND ProcessingNodeNumber=$processingNodeNumber AND isActive=1");#Terminal modifications
		my $queryTermMod="SELECT PeptideID, DeltaMass, ModificationName FROM PeptidesTerminalModifications$tableFlag, AminoAcidModifications WHERE PeptidesTerminalModifications$tableFlag.TerminalModificationID=AminoAcidModifications.AminoAcidModificationID AND ProcessingNodeNumber=$processingNodeNumber $activeModStrg";
		my $sthTermMod = $dbsqlite->prepare($queryTermMod);#Terminal modifications
		$sthTermMod->execute;
		while( my($peptideID,$deltaMass,$modificationName) = $sthTermMod->fetchrow_array) { #,$isActive
			#if($isActive == 1){
				next unless $peptidesInfo{$peptideID};
				$peptidesInfo{$peptideID}[3]+=$deltaMass;
				foreach my $modif (keys %{$refModifications}) {
					# Remove on 03/10/12
					#next if $refModifications->{$modif}{'PRINTED_NAME'} !~ /[NC]-term/; # terminal Modif => skip regular ones [eg. Acetyl(Protein N-term) vs Acetyl (K)]
					next if $refModifications->{$modif}{'PRINTED_NAME'} !~ /[NC][-_]term/i; # added by PP 10/01/13, modified 01/07/16
					if ($modificationName =~ /$refModifications->{$modif}{'SUBST'}/) {
						substr($peptidesInfo{$peptideID}[5] , 0 , 1 , $refModifications->{$modif}{'VALUE'});
					}
				}
			#}
		}
		$sthTermMod->finish;
		print '/'; # 3

		##> PhosphoRS information
		if ($phosphoRSMeasures{'phosphoRS Site Probabilities'}[0]) {
			my $sthPRS=$dbsqlite->prepare("SELECT PeptideID,FieldValue FROM CustomDataPeptides$tableFlag WHERE FieldID=$phosphoRSMeasures{'phosphoRS Site Probabilities'}[0]");
			$sthPRS->execute;
			while( my($peptideID,$prsData) = $sthPRS->fetchrow_array) {
				next unless $peptidesInfo{$peptideID};
				$peptidesInfo{$peptideID}[8]=$prsData;
			}
			$sthPRS->finish;
		}

		#####################################################
		###<Get all the peptides associated to this query>###
		#####################################################
#print "6 $typeName>";
		#print FILE "$boundary\n"; #<--- Boundary already written at end of previous section
		print FILE "Content-Type: application/x-Mascot; name=\"$sectionName\"\n\n";
		$count=0;
		foreach my $queryNum (sort{$a<=>$b} keys %{$refQueries}) {
			$count++;
			if ($count>=2000) {
				$count=0;
				print '.';
			}
			my $spectrumID=$refQueries->{$queryNum};
			#my $peptideID=$spectrumPeptides{$query};#PeptideID is an array which contains all the peptideIDs related to a spectrum
			if (!defined $spectrumPeptides{$spectrumID}){
				print FILE "q$queryNum",'_p1=-1',"\n";
				next;
			}
			my $i=0;
			#foreach my $pep (@{$peptideID}){ #}
			foreach my $peptideID (sort{$peptidesInfo{$b}[0]<=>$peptidesInfo{$a}[0] || $peptidesInfo{$a}[1] cmp $peptidesInfo{$b}[1]} @{$spectrumPeptides{$spectrumID}}) { # score - sequence
				my ($miscleavages,$proteinString,$termsString)=(0,'','');
				if($proteomeDiscovererVersion >= 1.3) {
					$miscleavages=$peptidesInfo{$peptideID}[4];
				}
				else {
					$miscleavages=($peptidesInfo{$peptideID}[1] =~ tr/K^R//);
					$miscleavages-=1 if $peptidesInfo{$peptideID}[1] =~/R$|K$/; # Prevent to forget to count a miss-cut in a C-Terminal peptide (not necessarily a K or a R that ends the sequence)
					#my $miscleavages=($peptidesInfo{$peptideID}[1] =~ tr/K^R//)-1;
					#$miscleavages=0 if $miscleavages==-1; #it is the case when it is C-Terminal type peptide (not necessarily a K or a R that ends the sequence)
				}

				if ($typeName eq 'Target') { #< Peptides ####################################
					foreach my $prot (@{$peptidesInfo{$peptideID}[6]}) {
						my $beg=index($proteinsInfo{$prot}[1],$peptidesInfo{$peptideID}[1])+1;
						my $end=$beg+length($peptidesInfo{$peptideID}[1])-1;
						my ($resbeg,$resend)=("","");
						if ($beg <= 1) {$resbeg="-";}
						else {$resbeg=substr($proteinsInfo{$prot}[1],$beg-2,1);}
						if ($end >= length($proteinsInfo{$prot}[1])) {$resend="-";}
						else {$resend=substr($proteinsInfo{$prot}[1],$end,1);}

						#Add all identifiers corresponding to the same protein
						foreach my $ident (@{$proteinsInfo{$prot}[0]}){
							$proteinString="$proteinString,\"$ident\":0:$beg:$end:0";
							$termsString="$termsString:$resbeg,$resend";
						}
					}
				}
				else { #< Decoy ###################################
					if ($proteomeDiscovererVersion >= 1.3) { # && scalar @{$peptidesInfo{$peptideID}[6]}
						foreach my $prot (@{$peptidesInfo{$peptideID}[6]}) {
							my $end=length($peptidesInfo{$peptideID}[1]);
							#Add all identifiers corresponding to the same protein
							foreach my $ident (@{$proteinsInfo{$prot}[0]}){
								$proteinString="$proteinString,\"$ident\":0:1:$end:0";
								$termsString="$termsString:-,-";
							}
						}
					}
					else {
						# Before v1.3: no PeptidesProteins_decoy table => all decoy peptides are linked to the same fake decoy protein in pdm
						# Decoy searches: there is not the same information in a mascot(.dat) and a proteomediscover(.msf) file
						# In MASCOT.DAT => there is the sequence and a link to a fake gi accession number
						# In SEQUEST.MSF => there is just a peptide and its score
						my $end=length($peptidesInfo{$peptideID}[1]);
						$proteinString=",\"no_identifier\":0:1:$end:0";
						$termsString=":-,-";
					}
				}

				# For Mascot research. Sometimes, data are not imported in the MSF file (there is a peptide with no protein related to it)
				if ($termsString) {
					$termsString = substr $termsString, 1; # removes starting ','
					$proteinString = substr $proteinString, 1; # removes starting ':'
					my $mass=&calculatePeptideMass($peptideID,$refMasses,\%peptidesInfo);
					my $delta=$refSpectra->{$spectrumID}[0]-$mass;
					$i++;
					print FILE "q$queryNum","_p$i=$miscleavages,$mass,";
					printf FILE "\%.6f",$delta;# Avoid to print e-05 for very small values...
					print FILE ",$peptidesInfo{$peptideID}[2],$peptidesInfo{$peptideID}[1],,$peptidesInfo{$peptideID}[5],$peptidesInfo{$peptideID}[0],,0,0;$proteinString\n";
					print FILE "q$queryNum","_p$i","_terms=$termsString\n";

					if ($peptidesInfo{$peptideID}[7]) { # Percolator data
						print FILE "q$queryNum","_p$i","_percolator=$peptidesInfo{$peptideID}[7]\n";
					}
					if ($peptidesInfo{$peptideID}[8]) { # PhosphoRS data
						print FILE "q$queryNum","_p$i","_phosphors=$peptidesInfo{$peptideID}[8]\n";
					}
				}
			}
			if ($i == 0) {
				print FILE "q$queryNum","_p1=-1\n";
			}
		}

		undef %peptidesInfo;
		undef %spectrumPeptides;

		print " Done.</B><BR>\n";

		print FILE "$boundary\n";
	}

}

#############################################################################
####<Walk up the Processing Nodes tree to looking for parent search node>####
#############################################################################
sub checkSearchNode { # same in send2Biologist.cgi
	my ($sthPN,$searchNodeNum,$parentNodeStrg)=@_;
	my $match=0;
	foreach my $parNodeNum (split(/;/,$parentNodeStrg)) {
		if ($parNodeNum==$searchNodeNum) {
			$match=1;
			last;
		}
	}
	unless ($match) {
		QNODE:foreach my $parNodeNum (split(/;/,$parentNodeStrg)) {
			$sthPN->execute($parNodeNum);
			while (my ($grandParentNodeStrg)=$sthPN->fetchrow_array) {
				$match=&checkSearchNode($sthPN,$searchNodeNum,$grandParentNodeStrg);
				last QNODE if $match;
			}
		}
	}
	return $match;
}

#################################################################
####<Compute the mass of a peptide according to his sequence>####
#################################################################
sub calculatePeptideMass {
	my ($peptideID,$refMasses,$refPeptidesInfo)=@_;
	my $mass=0;
	foreach my $aa (split(//,$refPeptidesInfo->{$peptideID}[1])) {
		$mass+=$refMasses->{$aa};
	}
	$mass+=$refMasses->{'C_term'}+$refMasses->{'N_term'};
	$mass+=$refPeptidesInfo->{$peptideID}[3]; # update the mass with the corresponding modifications (variable & fix!!!)
	return $mass;
}


# SEQUEST Monoisotopic Masses
sub getSequestMonoMass {
    my $masses=$_[0];
    $masses->{'A'}=71.03711375;
    $masses->{'B'}=114.534935165000001;
    $masses->{'C'}=103.00918445;
    $masses->{'D'}=115.02694295;
    $masses->{'E'}=129.04259301;
    $masses->{'F'}=147.06841387;
    $masses->{'G'}=57.02146369;
    $masses->{'H'}=137.05891180999998;
    $masses->{'I'}=113.08406393;
    $masses->{'J'}=113.08406393;
    $masses->{'K'}=128.09496295999998;
    $masses->{'L'}=113.08406393;
    $masses->{'M'}=131.04048457;
    $masses->{'N'}=114.04292738;
    $masses->{'O'}=237.14772677;
    $masses->{'P'}=97.05276381;
    $masses->{'Q'}=128.05857744;
    $masses->{'R'}=156.10111095999997;
    $masses->{'S'}=87.03202835;
    $masses->{'T'}=101.04767841;
    $masses->{'U'}=150.95363555;
    $masses->{'V'}=99.06841387;
    $masses->{'W'}=186.0793129;
    $masses->{'X'}=0.000000;
    $masses->{'Y'}=163.06332847000002;
    $masses->{'Z'}=128.550585225;
    $masses->{'C_term'}=17.00273963;
    $masses->{'N_term'}=1.00782503;
    $masses->{'Hydrogen'}=1.007825;
    $masses->{'Carbon'}=12.000000;
    $masses->{'Nitrogen'}=14.003074;
    $masses->{'Oxygen'}=15.994915;
    $masses->{'Electron'}=0.000549;
}

# Mascot Monoisotopic Masses
sub getMascotMonoMass {
    my $masses=$_[0];
    $masses->{'A'}=71.037114;
    $masses->{'B'}=114.534940;
    $masses->{'C'}=103.009185;
    $masses->{'D'}=115.026943;
    $masses->{'E'}=129.042593;
    $masses->{'F'}=147.068414;
    $masses->{'G'}=57.021464;
    $masses->{'H'}=137.058912;
    $masses->{'I'}=113.084064;
    $masses->{'J'}=0.000000;
    $masses->{'K'}=128.094963;
    $masses->{'L'}=113.084064;
    $masses->{'M'}=131.040485;
    $masses->{'N'}=114.042927;
    $masses->{'O'}=0.000000;
    $masses->{'P'}=97.052764;
    $masses->{'Q'}=128.058578;
    $masses->{'R'}=156.101111;
    $masses->{'S'}=87.032028;
    $masses->{'T'}=101.047679;
    $masses->{'U'}=150.953630;
    $masses->{'V'}=99.068414;
    $masses->{'W'}=186.079313;
    $masses->{'X'}=111.000000;
    $masses->{'Y'}=163.063329;
    $masses->{'Z'}=128.550590;
    $masses->{'C_term'}=17.002740;
    $masses->{'N_term'}=1.007825;
    $masses->{'Hydrogen'}=1.007825;
    $masses->{'Carbon'}=12.000000;
    $masses->{'Nitrogen'}=14.003074;
    $masses->{'Oxygen'}=15.994915;
    $masses->{'Electron'}=0.000549;
}

############################################################
####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!####
####< Same functions for new MSF files for PD 2.2 only >####
####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!####
############################################################

##########################################
####<Print header section (.pdm file)>####
##########################################
sub printHeaderMSF2_2 { # GLOBALS: $dbsqlite, $processingNodeNumber, FILE
	print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>-Header...";
    my ($db,$nbquery)=@_;
    my ($date)=$dbsqlite->selectrow_array("SELECT WorkflowStartDate FROM Workflows");
    my ($dateCode,$startTime)=split(/\s/,$date);
	$startTime=~s/\.\d+//; # 12:34:05.625 -> 12:34:05
    my ($dbpath)=($refSearches->{$processingNodeNumber}{'DB_NAME'})? $refSearches->{$processingNodeNumber}{'DB_NAME'} : '';
    my ($friendlyName)=$refSearches->{$processingNodeNumber}{DESCRIPTOR};
    #my ($friendlyName)=$dbsqlite->selectrow_array("SELECT FriendlyName FROM ProcessingNodes WHERE ProcessingNodeNumber=$processingNodeNumber");

    my ($dbName,$mascotVersion,$dbNbSeq,$dbNbSeqAT,$nbQueries)=("","","","","");
#    if (uc($friendlyName)=~/MASCOT/ && $proteomeDiscovererVersion >= 1.3){# In PD 1.3+, some other information is available in CustomDataProcessingNodes
#		my ($dbInfo)=$dbsqlite->selectrow_array("SELECT FieldValue FROM CustomDataProcessingNodes WHERE ProcessingNodeNumber=$processingNodeNumber AND FieldID=(SELECT FieldID FROM CustomDataFIelds WHERE DisplayName='Fasta database information')");
#		foreach my $field (split("\n",$dbInfo)) {
#			my ($title,$value)=split(": ",$field);
#			if($title eq 'Version'){ $mascotVersion=$value;}
#			elsif($title eq 'Number of sequences'){ $dbNbSeq=$value;}
#			elsif($title eq 'Number of sequences after taxonomy'){ $dbNbSeqAT=$value;}
#			elsif($title eq 'FASTA db'){ $db=$value;}
#			#elsif($title eq 'Number of queries'){ $nbquery=$value;}
#		}
#    }

	print FILE qq
|Content-Type: application/x-Mascot; name="header"

sequences=$dbNbSeq
sequences_after_tax=$dbNbSeqAT
residues=
distribution=
exec_time=
date=
time=$startTime
queries=$nbquery
max_hits=
version=$mascotVersion
fastafile=$dbpath
release=$db
taskid=
$boundary
|;
	print " Done.</B><BR>\n";
}

sub printParametersMSF2_2 { # GLOBALS: $dbsqlite, $processingNodeNumber, $refSearches, $proteomeDiscovererVersion, FILE, $boundary
	my ($refModifications,$refMasses)=@_;
	print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>-Parameters...";

	my ($file,$db,$tol,$tolu,$itol,$itolu,$cle,$pfa,$decoy,$taxonomy,$instrument)=('','','','','','','','','','','');

	#my ($rawfilenumber) = $dbsqlite->selectrow_array("SELECT ProcessingNodeNumber FROM ProcessingNodes WHERE FriendlyName='Spectrum Files'");
	my ($wfXML) = $dbsqlite->selectrow_array("SELECT WorkflowXML FROM Workflows");
    my $xml = new XML::Simple();
	my $xmlData = $xml->XMLin($wfXML);
	foreach my $processingNodeParameters (sort{$a->{ParentProcessingNodeNumber} cmp $b->{ParentProcessingNodeNumber}} @{$xmlData->{WorkflowTree}->{WorkflowNode}}) {
		if ($processingNodeParameters->{ProcessingNodeName} eq 'SpectrumFilesNode') {
				my @msfFiles=split(/;/,$processingNodeParameters->{ProcessingNodeParameters}->{ProcessingNodeParameter}->{content});
				my $addQuery=($fileID)?" AND FileID=$fileID":'';
				foreach my $mFile (@msfFiles){
					my @name=split(/\\/,$mFile);
					my $fName=pop(@name);
					my $fPath=join("\\",@name);
					$file.="File Name: $fName\n";
					$file.="File Path: \\$fPath\n";
				}
		}
		next unless ($processingNodeParameters->{ProcessingNodeNumber} == $processingNodeNumber);
		foreach my $processingNodeParameter (@{$processingNodeParameters->{ProcessingNodeParameters}->{ProcessingNodeParameter}}){
				if ($processingNodeParameter->{Name} eq 'Database'){
					$db=$processingNodeParameter->{content};
				}
				elsif ($processingNodeParameter->{Name} eq 'FastaDatabase'&& uc($refSearches->{$processingNodeNumber}{'DESCRIPTOR'})=~/SEQUEST/){
					$db=$processingNodeParameter->{content};
				}
				elsif($processingNodeParameter->{Name} eq 'PeptideTolerance'){
					($tol,$tolu)=split(/\s/,$processingNodeParameter->{content});
				}
				elsif($processingNodeParameter->{Name} eq 'FragmentTolerance'){
					($itol,$itolu)=split(/\s/,$processingNodeParameter->{content});
				}
				elsif($processingNodeParameter->{IntendedPurpose} eq 'CleavageReagent'){
					my $enzymeXML=new XML::Simple();
					my $enzymeDesc=$enzymeXML->XMLin($processingNodeParameter->{content}); # <Enzyme Version="1" Name="Trypsin" CleavageSites="KR" CleavageInhibitors="P" Offset="1" CleavageSpecificity="SpecificAtBothEnds" />
					$cle=$enzymeDesc->{Name};
				}
				elsif ($processingNodeParameter->{Name} eq 'Taxonomy'){
					$taxonomy=$processingNodeParameter->{content};
				}
				elsif($processingNodeParameter->{Name} =~ /(Max)?MissedCleavages/){
					$pfa=$processingNodeParameter->{content};
				}
				elsif($processingNodeParameter->{Name} eq 'UseDecoyDatabase' || $processingNodeParameter->{Name} eq 'DecoySearch'){
					$decoy=($processingNodeParameter->{IsValueSet} eq 'True')? 1 : 0;
				}
				elsif($processingNodeParameter->{Name} eq 'Instrument'){
					$instrument=$processingNodeParameter->{content};
				}
		}
	    $refSearches->{$processingNodeNumber}{'DB_NAME'}=$db;
    }

	my ($isSequest,$minScorePep,$com); #to know if the dat that is created is generated from an MSF file
	if (uc($refSearches->{$processingNodeNumber}{'DESCRIPTOR'})=~/SEQUEST/){
		&getSequestMonoMass($refMasses);
		$isSequest="SEQUEST";
		if ($decoy) {$minScorePep=0;} # Assumes Qvality-like FDR
		elsif ($selMinScore=~/def/i) {$minScorePep=&promsConfig::getMinScore('SEQUEST.PDM');} # default
		else {$minScorePep=$selMinScore;}
		($com)=$dbsqlite->selectrow_array("SELECT WorkflowName FROM Workflows WHERE rowid=1"); # get the name of the project Workflow
	}
	else{ # MASCOT
		&getMascotMonoMass($refMasses);
		$isSequest="MASCOT";
		if ($decoy) {$minScorePep=0;} # Assumes FDR-based filtering
		elsif ($selMinScore=~/def/i) {$minScorePep=&promsConfig::getMinScore('MASCOT.PDM');} # default
		else {$minScorePep=$selMinScore;}
		my ($message) = $dbsqlite->selectrow_array("SELECT Message FROM WorkflowMessages WHERE Message LIKE \"Received Mascot result file%\"");
		#my @result=split(/[\/\\]/,$message);# $message="Mascot result on server (filename=../data/20100803/F008722.dat)";
		#$com=$result[-1];# $com="F008722.dat)";
		#$com=substr($com,0,length($com)-1);# $com="F008722.dat";
		($com)=($message=~/([^\/\\]+)\)\Z/); # "Mascot result on server (filename=../data/20100803/F008722.dat) -> F008722.dat
	}

	#########################################################################################################
	###>Get all the lines corresponding to modifications (static and dynamic) in ProcessingNodeParameters<###
	#########################################################################################################
	my $sthDelMass=$dbsqlite->prepare("SELECT DeltaMonoisotopicMass,ModificationID FROM FoundModifications WHERE Name=? LIMIT 0,1"); # multiple rows are possible
	my $sthNLInfo=$dbsqlite->prepare("SELECT NeutralLossName,NeutralLossID,NeutralLossMass FROM FoundModificationsAminoAcids FMAA, AminoAcids AA WHERE AA.AminoAcidID=FMAA.AminoAcidsAminoAcidID AND FoundModificationsModificationID=? AND AA.OneLetterCode=? ORDER BY NeutralLossID ASC");

	#Fill the static and dynamic modification informations
	my ($mods,$itmods,$modsString,$itmodsString)=('','','','');
	my ($i,$j)=(1,1);
	foreach my $processingNodeParameters (@{$xmlData->{WorkflowTree}->{WorkflowNode}}) {
		next unless $processingNodeParameters->{ProcessingNodeNumber} == $processingNodeNumber;
		foreach my $processingNodeParameter (@{$processingNodeParameters->{ProcessingNodeParameters}->{ProcessingNodeParameter}}){
				if(($processingNodeParameter->{IntendedPurpose} =~ /Static(Terminal)?Modification/) && $processingNodeParameter->{IsValueSet} eq "True"){
					my $statmodXML=new XML::Simple();
					my $statmodDesc=$statmodXML->XMLin($processingNodeParameter->{content}); #<Modification Version="2" AminoAcids="X" Name="MappingL" Abbreviation="L" ID="1167" UnimodAccession="-1" DeltaMass="113.08406" DeltaAverageMass="113.15890" IsSubstitution="True" LeavingGroup="" Substitution="C6H11NO" PositionType="Any" />
					next if $statmodDesc->{Name} eq "MascotXValue" || $statmodDesc->{Name} eq "MappingL";
					my $delta=$statmodDesc->{DeltaMass};
					my $res=($statmodDesc->{AminoAcids})? $statmodDesc->{AminoAcids} : $statmodDesc->{Terminus};
					if ($statmodDesc->{Terminus}) { # converts "Any_N_Terminus" to standard "Any N-Terminus" (PP 01/07/16)
						$res=~s/Any_/Any /i;
						$res=~s/_Term/-Term/i;
					}
					my $valueD="$statmodDesc->{Name} ($res)";
					$mods.=",$valueD";
					$modsString.="FixedMod$i=$delta,$valueD\nFixedModResidues$i=$res\n";
					$i++;
				}
				elsif(($processingNodeParameter->{IntendedPurpose} =~ /Dynamic(Terminal)?Modification/) && $processingNodeParameter->{IsValueSet} eq "True"){
					my $dynmodXML=new XML::Simple();
					my $dynmodDesc=$dynmodXML->XMLin($processingNodeParameter->{content}); #<Modification Version="2" AminoAcids="C" Name="Carbamidomethyl" Abbreviation="Carbamidomethyl" ID="8" UnimodAccession="4" DeltaMass="57.02146" DeltaAverageMass="57.05130" IsSubstitution="False" LeavingGroup="" Substitution="H(3) C(2) N O" PositionType="Any" />
					my $res=($dynmodDesc->{AminoAcids})? $dynmodDesc->{AminoAcids} : $dynmodDesc->{Terminus};
					if ($dynmodDesc->{Terminus}) { # converts "Any_N_Terminus" to standard "Any N-Terminus" (PP 01/07/16)
						$res=~s/Any_/Any /i;
						$res=~s/_Term/-Term/i;
					}
					my $valueD="$dynmodDesc->{Name} ($res)";
					my $paramN=$processingNodeParameter->{Name};
					$refModifications->{$paramN}{'NAME'}=$valueD;
					$refModifications->{$paramN}{'SUBST'}=$dynmodDesc->{Name}; # for terminal modification
					$refModifications->{$paramN}{'DELTA'}=$dynmodDesc->{DeltaMass};
					$refModifications->{$paramN}{'ID'}=$dynmodDesc->{ID};
					$refModifications->{$paramN}{'VALUE'}=$j;
					$valueD=~ s/, //g;
					$itmods.=",$valueD";
					$itmodsString.="delta$j=$dynmodDesc->{DeltaMass},$valueD\n";
					my $addNLString='';
					if ($dynmodDesc->{AminoAcids}) {
						my %nl;# nl <-> Neutral-Losses
						foreach my $aa (split(//,$dynmodDesc->{AminoAcids})) {
							$sthNLInfo->execute($dynmodDesc->{ID},$aa);
							while (my ($nlName,$nlID,$nlMass)=$sthNLInfo->fetchrow_array){
								$nl{$nlID}=$nlMass;
							}
						}
						foreach my $nlID (sort{$a<=>$b} keys %nl) {# Print neutral-loss like Mascot (v2.3) DAT files.
							next unless $nl{$nlID};
							if ($addNLString) {
								$addNLString.="NeutralLoss${j}_master=$nl{$nlID}\n";
							}else{
								$addNLString.="NeutralLoss$j=$nl{$nlID}\n";
							}
						}
						###> For Phospho(Y) var-mod, must write specifically this line.
						if (!$addNLString && $dynmodDesc->{Name} =~ /Phospho/ && $dynmodDesc->{AminoAcids} eq 'Y') {
							$addNLString.="NeutralLoss$j=0.000000\n";
						}#code
					}
					$itmodsString.=$addNLString;
					$refModifications->{$paramN}{'PRINTED_NAME'}=$valueD;
					$j++;
				}
		}
	}

	$mods = ($mods) ? (substr $mods, 1) : ''; # remove starting ','
	$itmods = ($itmods) ? (substr $itmods, 1) : ''; # remove starting ','

	##############################################################
	###>Find Quantification if any in ProcessingNodeParameters<###
	##############################################################
	my ($quantifName,$paramValue,$processMatched,$processCount)=('','',0,0);
	### PrecursorIonsQuantifierNode -> SILAC
	### ReporterIonQuantifierNode -> TMT
	foreach my $processingNodeParameters (@{$xmlData->{WorkflowTree}->{WorkflowNode}}) {
		next unless ($processingNodeParameters->{ProcessingNodeName} =~ /PrecursorIonsQuantifierNode|ReporterIonQuantifierNode|MinoraFeatureCreatorNode/);
		foreach my $processingNodeParameter (@{$processingNodeParameters->{ProcessingNodeParameters}->{ProcessingNodeParameter}}){
				if($processingNodeParameter->{IntendedPurpose} eq "QuantificationMethod"){
					$paramValue=$processingNodeParameter->{content};
				    next unless $paramValue;
					($quantifName)=($paramValue=~/^<ProcessingMethod name=\"(.*)\" version/);
				}
		}
		if ($paramValue) {
				my $parentNodeStrg=$processingNodeParameters->{ParentProcessingNodeNumber};
				$processMatched=0;
				foreach my $parNodeNum (split(/;/,$parentNodeStrg)) {
					if ($parNodeNum==$processingNodeNumber) {
						$processMatched=1;
						last;
					}
				}
		}
	}

	if ($processMatched==0 && $processCount > 1) { # multiple quantifs that cannot be linked to the $processingNodeNumber used => skip quantif info
		$quantifName=$paramValue='';
		print "[*Unable to match quantification process to search process: Skipping quantification info*]";
	}
	my $isPercolator=($refSearches->{$processingNodeNumber}{PERCOLATOR_NODE})? 'Percolator' : '';
	my $isPhosphoRS=($refSearches->{$processingNodeNumber}{PHOSPHORS_NODE})? 1 : 0;

	#######################
	###>Writing to file<###
	#######################
#LICENSE=Licensed to: Institut Curie, Lab. de Spectrom. de Mass Proteomique (No. 2G-49264 Loew), (4 processors).
	print FILE qq
|MIME-Version: 1.0 (Generated by Mascot version 1.0)
Content-Type: multipart/mixed; boundary=gc0p4Jq0M2Yt08jU534c0p

$boundary
Content-Type: application/x-Mascot; name="parameters"

LICENSE=
MP=
NM=
COM=$com
IATOL=
IA2TOL=
IASTOL=
IBTOL=
IB2TOL=
IBSTOL=
IYTOL=
IY2TOL=
IYSTOL=
SEG=
SEGT=
SEGTU=
LTOL=
TOL=$tol
TOLU=$tolu
ITH=
ITOL=$itol
ITOLU=$itolu
PFA=$pfa
DB=$db
MODS=$mods
MASS=Monoisotopic
CLE=$cle
FILE=$file
PEAK=
QUE=
TWO=
SEARCH=MIS
USERNAME=
USEREMAIL=
CHARGE=2+ and 3+
INTERMEDIATE=
REPORT=AUTO
OVERVIEW=
FORMAT=Mascot generic
FORMVER=1.01
FRAG=
IT_MODS=$itmods
USER00=
USER01=
USER02=
USER03=
USER04=
USER05=
USER06=
USER07=
USER08=
USER09=
USER10=
USER11=
USER12=
PRECURSOR=
TAXONOMY=$taxonomy
ACCESSION=
REPTYPE=
SUBCLUSTER=
ICAT=
INSTRUMENT=$instrument
ERRORTOLERANT=
FRAMES=
CUTOUT=
USERID=0
QUANTITATION=$quantifName
DECOY=$decoy
SEARCH_ALGO=$isSequest
FDR_ALGO=$isPercolator
PHOSPHORS=$isPhosphoRS
PEP_ISOTOPE_ERROR=
RULES=1,2,8,9,10,13,14,15
INTERNALS=0.0,700.0
$boundary
Content-Type: application/x-Mascot; name="masses"

A=$refMasses->{'A'}
B=$refMasses->{'B'}
C=$refMasses->{'C'}
D=$refMasses->{'D'}
E=$refMasses->{'E'}
F=$refMasses->{'F'}
G=$refMasses->{'G'}
H=$refMasses->{'H'}
I=$refMasses->{'I'}
J=$refMasses->{'J'}
K=$refMasses->{'K'}
L=$refMasses->{'L'}
M=$refMasses->{'M'}
N=$refMasses->{'N'}
O=$refMasses->{'O'}
P=$refMasses->{'P'}
Q=$refMasses->{'Q'}
R=$refMasses->{'R'}
S=$refMasses->{'S'}
T=$refMasses->{'T'}
U=$refMasses->{'U'}
V=$refMasses->{'V'}
W=$refMasses->{'W'}
X=$refMasses->{'X'}
Y=$refMasses->{'Y'}
Z=$refMasses->{'Z'}
Hydrogen=$refMasses->{'Hydrogen'}
Carbon=$refMasses->{'Carbon'}
Nitrogen=$refMasses->{'Nitrogen'}
Oxygen=$refMasses->{'Oxygen'}
Electron=$refMasses->{'Electron'}
C_term=$refMasses->{'C_term'}
N_term=$refMasses->{'N_term'}
$itmodsString$modsString
$boundary
|;

	if ($paramValue) {# Print the XML information of the quantitation method used by PD for SILAC !
		###> Add a synonym to this modification...
		$paramValue=~s/""/"/g;
		my $xml = new XML::Simple();
		my $xmlData = $xml->XMLin($paramValue);
		my $labelXLM=$xmlData->{MethodPart}{QuanChannels}{MethodPart};

		my $dbh=&promsConfig::dbConnect;
		foreach my $labelName (keys %{$labelXLM}){
			if ($labelXLM->{$labelName}{MethodPart}{MethodPart}{MethodPart}) { # No label or single label
				my $labelModifAlias=$labelXLM->{$labelName}{MethodPart}{MethodPart}{name};
				if ($labelXLM->{$labelName}{MethodPart}{MethodPart}{MethodPart}{MethodPart}{SideChainModification}{Parameter}) {
					my $labelParamStrg=$labelXLM->{$labelName}{MethodPart}{MethodPart}{MethodPart}{MethodPart}{SideChainModification}{Parameter}{content};
					my ($labelModifName)=$labelParamStrg=~/ Name="([^"]+)/;
					my ($modifRes)=$labelParamStrg=~/ AminoAcids="([^"]+)/;
					my $modificationID=&promsMod::getModificationIDfromString($dbh,$labelModifName,$modifRes,\%vmodsInUnimod);
					#&addSynonyms($dbh,$modificationID,$labelModifAlias);
					#&addSynonyms($dbh,$modificationID,$quantifName);
				}
			}
			else {
				foreach my $labelModifAlias (sort{lc($a) cmp lc($b)} keys %{$labelXLM->{$labelName}{MethodPart}{MethodPart}}) {
					my $labelParamStrg=$labelXLM->{$labelName}{MethodPart}{MethodPart}{$labelModifAlias}{MethodPart}{MethodPart}{SideChainModification}{Parameter}{content};
					my ($labelModifName)=$labelParamStrg=~/ Name="([^"]+)/;
					my ($modifRes)=$labelParamStrg=~/ AminoAcids="([^"]+)/;
					my $modificationID=&promsMod::getModificationIDfromString($dbh,$labelModifName,$modifRes,\%vmodsInUnimod);
					#&addSynonyms($dbh,$modificationID,$labelModifAlias);
					#&addSynonyms($dbh,$modificationID,$quantifName);
				}
			}
		}
		$dbh->disconnect;

		print FILE qq
|Content-Type: application/x-Mascot; name="quantitation"

$paramValue
$boundary
|;
	}
	print " Done.</B><BR>\n";
	return ($db,$decoy,$minScorePep);
}


###################################################################
####<Print Summary: the queries one by one with the ions found>####
###################################################################
sub printPeptidesMSF2_2 { # GLOBALS: $pdmFile, $dbsqlite, $processingNodeNumber, $refSearches, FILE, $databankID, $userID
	my ($decoy,$minScorePep,$refQueries,$refSpectra,$refModifications,$refMasses)=@_;

	########################################
	####<Get protein info from MSF file>####
	########################################
#print '1>';
	print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>-Proteins...";
	my %proteinsInfo;
	my $dbh=&promsConfig::dbConnect;
	my ($parseRules,$identType,$defIdentType)=$dbh->selectrow_array("SELECT PARSE_RULES,IDENTIFIER_TYPE,DEF_IDENT_TYPE FROM DATABANK,DATABANK_TYPE WHERE DATABANK.ID_DBTYPE=DATABANK_TYPE.ID_DBTYPE AND ID_DATABANK=$databankID"); # $databankID global. Only 1 db-search allowed in MSF
	$defIdentType='UNKNOWN' unless $defIdentType;
	$dbh->disconnect;
	$identType=$defIdentType unless $identType;
	my @rules=split(',:,',$parseRules);
	my ($idRule)=($rules[0]=~/ID=(.+)/); #<<<<<<<<<<<<<< Synchronise with promsMod::getProtInfo >>>>>>
	my $tempDbFile;
	if ($fileID) {
		($tempDbFile ="$promsPath{tmp}/$msfFile")=~s/\.msf\Z/_$processingNodeNumber\.$fileID\.$userID.fasta/;
	}else{
		($tempDbFile ="$promsPath{tmp}/$msfFile")=~s/\.msf\Z/_$processingNodeNumber\.$userID.fasta/; # same dir as msf file
	}
	open (FASTA,">$tempDbFile");
#	print "";
	my $sthProt = $dbsqlite->prepare("SELECT UniqueSequenceID, Sequence, FastaTitleLines, Accession FROM TargetProteins");
	$sthProt->execute;
	while (my ($proteinID,$sequence,$description,$accession)= $sthProt->fetchrow_array) {
		print FASTA "$description\n$sequence\n";
		@{$proteinsInfo{$proteinID}}=(
			[$accession], 		# 0 identifiers
			$sequence	# 1
		);
	}
	$sthProt->finish;
	close FASTA;
	print " Done.</B><BR>\n";

	########################################################
	####<Get peptide info: Looping from target to decoy>####
	########################################################
	#>Get Scoring type used
	my ($scoreID,$scoreName) = $dbsqlite->selectrow_array("SELECT MIN(ScoreID),ScoreName FROM ProcessingNodeScores WHERE NodeGuid='$refSearches->{$processingNodeNumber}{NODEGUID}'");
	my @peptideTypes=(['Target','','peptides']);
	push @peptideTypes,['Decoy','_decoy','decoy_peptides'] if $decoy;

	#my %phosphoRSMeasures;
	#if ($refSearches->{$processingNodeNumber}{'PHOSPHORS_NODE'}) {
	#	my $sthPRScol=$dbsqlite->prepare("SELECT DisplayName,FieldID FROM CustomDataFields WHERE SourceNodeNumber=$refSearches->{$processingNodeNumber}{PHOSPHORS_NODE} AND TargetNodeNumber=$processingNodeNumber AND (DisplayName LIKE 'phosphoRS%')");
	#	# 3 Informations provided in CustomDataFields :
	#	# - phosphoRS Binomial Peptide Score
	#	# - phosphoRS Isoform Probability
	#	# - phosphoRS Site Probabilities
	#	$sthPRScol->execute;
	#	while (my($displayName,$fieldID)=$sthPRScol->fetchrow_array) {
	#		push @{$phosphoRSMeasures{"$displayName"}},$fieldID; # 3 entries: 1st Binomial Peptide Score, 2nd Isoform Probability and 3rd Site Probabilities
	#	}
	#	$sthPRScol->finish;
	#}
    my $workFlowID;
	foreach my $refType (@peptideTypes) { # Target (and Decoy)
		my ($typeName,$tableFlag,$sectionName)=@{$refType};

		print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>-$typeName peptides...";
		my (%peptidesInfo,%spectrumPeptides); #,%peptidesDecoyInfo,%spectrumDecoyPeptides

		######################################################################
		###<Get all the information from Peptides table in SQLite datafile>###
		######################################################################
#print "2 $typeName>";
		###<Percolator filtering>###
		my $sthPepInfo=$dbsqlite->prepare("SELECT PeptideID,$scoreName,MSnSpectrumInfoSpectrumID,MatchedIonsCount,Sequence,MissedCleavages,PercolatorqValue,PercolatorPEP FROM ${typeName}Psms TP, ${typeName}PsmsMSnSpectrumInfo TPMSI WHERE TP.WorkflowID=TPMSI.${typeName}PsmsWorkflowID AND TP.PeptideID=TPMSI.${typeName}PsmsPeptideID");
		#my $step=5000; #1000; # int($numPeptides/10);
		my $count=0;
		$sthPepInfo->execute;
		while (my($peptideID,$score,$spectrumID,$matchedIonsCount,$sequence,$missCut,$qVal,$PEP) = $sthPepInfo->fetchrow_array) {
			$count++;
			if ($count>=10000) { # 2500
				$count=0;
				print '.';
			}
			next unless $spectrumInfo{$spectrumID};
			next if $qVal && $qVal > $percolatorThres;

			@{$peptidesInfo{$peptideID}}=(
				$score,							# 0
				$sequence,						# 1
				$matchedIonsCount,				# 2
				0.000000,						# 3
				$missCut,						# 4 only defined if $proteomeDiscovererVersion >= 1.3;
				"0" x (length($sequence)+2),	# 5 String that explain the modifications on the peptides
				[]								# 6 matching proteins
			);
			push @{$peptidesInfo{$peptideID}},"$qVal:$PEP" if $qVal ; # 7 Percolator data
			push @{$spectrumPeptides{$spectrumID}},$peptideID; # 1 spectrum <-> several peptides
		}
		$sthPepInfo->finish;
		print '/';

#print "3 $typeName>";
		if ($typeName eq 'Target') { # Before v1.3: no PeptidesProteins_decoy table
		    ($workFlowID)=$dbsqlite->selectrow_array("SELECT WorkflowID FROM Workflows");
			my $sthPepProt = $dbsqlite->prepare("SELECT TargetPsmsPeptideID, TargetProteinsUniqueSequenceID FROM TargetProteinsTargetPsms WHERE TargetPsmsWorkflowID=$workFlowID");
			$sthPepProt->execute;
			while( my($peptideID,$proteinID) = $sthPepProt->fetchrow_array){
				next unless $peptidesInfo{$peptideID};
				push @{$peptidesInfo{$peptideID}[6]},$proteinID;
			}
			$sthPepProt->finish;
			print '.';
		}

		#######################################################################
		###<Update the mass value of the peptide according to modifications>###
		#######################################################################
		##>Regular modifications
#print "4 $typeName>";
		###> !!! WARNING !!! Position is not stored the same in PD 2_2 and lower versions
		my $queryModif="SELECT Position, ${typeName}PsmsPeptideID, DeltaMonoisotopicMass, Name, UnimodAccession FROM ${typeName}PsmsFoundModifications PFM, FoundModifications FM WHERE PFM.FoundModificationsModificationID=FM.ModificationID AND ${typeName}PsmsWorkflowID=$workFlowID ORDER BY ${typeName}PsmsPeptideID, Position ASC";
		my $sthRegMod = $dbsqlite->prepare($queryModif);#Normal modifications
		$sthRegMod->execute;
		$count=0;
		while (my($position,$peptideID,$deltaMass,$modificationName,$unimodID) = $sthRegMod->fetchrow_array) { #,$isActive
			$position--;
			$count++;
			if ($count>=20000) {
				$count=0;
				print '.';
			}
			next unless $peptidesInfo{$peptideID};
			$modificationName = quotemeta(uc($modificationName)); # quote because match below
			$peptidesInfo{$peptideID}[3]+=$deltaMass;
			my @sequence = split(//,$peptidesInfo{$peptideID}[1]);
			foreach my $modif (keys %{$refModifications}) {
				#if (uc($refModifications->{$modif}{'NAME'}) =~ /$modificationName \(\w*$sequence[$position]\w*\)/){
				#if (uc($refModifications->{$modif}{'PRINTED_NAME'}) =~ /^$modificationName \(\w*$sequence[$position]\w*\)/){
				if (uc($refModifications->{$modif}{'PRINTED_NAME'}) =~ /^$modificationName \(\w*$sequence[$position]\w*\)/){ # if ' is forgotten, then, for Methyl, Dimethyl, Trimethyl modification, it does not find the correct one!
					substr($peptidesInfo{$peptideID}[5] , $position+1 , 1 , $refModifications->{$modif}{'VALUE'} );
				}
			}
		}
		$sthRegMod->finish;
		print '/'; # 3

		###> PhosphoRS information
		#if ($phosphoRSMeasures{'phosphoRS Site Probabilities'}[0]) {
		#	my $sthPRS=$dbsqlite->prepare("SELECT PeptideID,FieldValue FROM CustomDataPeptides$tableFlag WHERE FieldID=$phosphoRSMeasures{'phosphoRS Site Probabilities'}[0]");
		#	$sthPRS->execute;
		#	while( my($peptideID,$prsData) = $sthPRS->fetchrow_array) {
		#		next unless $peptidesInfo{$peptideID};
		#		$peptidesInfo{$peptideID}[8]=$prsData;
		#	}
		#	$sthPRS->finish;
		#}

		#####################################################
		###<Get all the peptides associated to this query>###
		#####################################################
#print "6 $typeName>";
		#print FILE "$boundary\n"; #<--- Boundary already written at end of previous section
		print FILE "Content-Type: application/x-Mascot; name=\"$sectionName\"\n\n";
		$count=0;
		foreach my $queryNum (sort{$a<=>$b} keys %{$refQueries}) {
			$count++;
			if ($count>=2000) {
				$count=0;
				print '.';
			}
			my $spectrumID=$refQueries->{$queryNum};
			#my $peptideID=$spectrumPeptides{$query};#PeptideID is an array which contains all the peptideIDs related to a spectrum
			if (!defined $spectrumPeptides{$spectrumID}){
				print FILE "q$queryNum",'_p1=-1',"\n";
				next;
			}
			my $i=0;
			#foreach my $pep (@{$peptideID}){ #}
			foreach my $peptideID (sort{$peptidesInfo{$b}[0]<=>$peptidesInfo{$a}[0] || $peptidesInfo{$a}[1] cmp $peptidesInfo{$b}[1]} @{$spectrumPeptides{$spectrumID}}) { # score - sequence
				my ($miscleavages,$proteinString,$termsString)=(0,'','');
				if($proteomeDiscovererVersion >= 1.3) {
					$miscleavages=$peptidesInfo{$peptideID}[4];
				}
				else {
					$miscleavages=($peptidesInfo{$peptideID}[1] =~ tr/K^R//);
					$miscleavages-=1 if $peptidesInfo{$peptideID}[1] =~/R$|K$/; # Prevent to forget to count a miss-cut in a C-Terminal peptide (not necessarily a K or a R that ends the sequence)
					#my $miscleavages=($peptidesInfo{$peptideID}[1] =~ tr/K^R//)-1;
					#$miscleavages=0 if $miscleavages==-1; #it is the case when it is C-Terminal type peptide (not necessarily a K or a R that ends the sequence)
				}

				if ($typeName eq 'Target') { #< Peptides ####################################
					foreach my $prot (@{$peptidesInfo{$peptideID}[6]}) {
						my $beg=index($proteinsInfo{$prot}[1],$peptidesInfo{$peptideID}[1])+1;
						my $end=$beg+length($peptidesInfo{$peptideID}[1])-1;
						my ($resbeg,$resend)=("","");
						if ($beg <= 1) {$resbeg="-";}
						else {$resbeg=substr($proteinsInfo{$prot}[1],$beg-2,1);}
						if ($end >= length($proteinsInfo{$prot}[1])) {$resend="-";}
						else {$resend=substr($proteinsInfo{$prot}[1],$end,1);}

						#Add all identifiers corresponding to the same protein
						foreach my $ident (@{$proteinsInfo{$prot}[0]}){
							$proteinString="$proteinString,\"$ident\":0:$beg:$end:0";
							$termsString="$termsString:$resbeg,$resend";
						}
					}
				}
				else { #< Decoy ###################################
					if ($proteomeDiscovererVersion >= 1.3) { # && scalar @{$peptidesInfo{$peptideID}[6]}
						foreach my $prot (@{$peptidesInfo{$peptideID}[6]}) {
							my $end=length($peptidesInfo{$peptideID}[1]);
							#Add all identifiers corresponding to the same protein
							foreach my $ident (@{$proteinsInfo{$prot}[0]}){
								$proteinString="$proteinString,\"$ident\":0:1:$end:0";
								$termsString="$termsString:-,-";
							}
						}
					}
					else {
						# Before v1.3: no PeptidesProteins_decoy table => all decoy peptides are linked to the same fake decoy protein in pdm
						# Decoy searches: there is not the same information in a mascot(.dat) and a proteomediscover(.msf) file
						# In MASCOT.DAT => there is the sequence and a link to a fake gi accession number
						# In SEQUEST.MSF => there is just a peptide and its score
						my $end=length($peptidesInfo{$peptideID}[1]);
						$proteinString=",\"no_identifier\":0:1:$end:0";
						$termsString=":-,-";
					}
				}

				# For Mascot research. Sometimes, data are not imported in the MSF file (there is a peptide with no protein related to it)
				if ($termsString) {
					$termsString = substr $termsString, 1; # removes starting ','
					$proteinString = substr $proteinString, 1; # removes starting ':'
					my $mass=&calculatePeptideMass($peptideID,$refMasses,\%peptidesInfo);
					my $delta=$refSpectra->{$spectrumID}[0]-$mass;
					$i++;
					print FILE "q$queryNum","_p$i=$miscleavages,$mass,";
					printf FILE "\%.6f",$delta;# Avoid to print e-05 for very small values...
					print FILE ",$peptidesInfo{$peptideID}[2],$peptidesInfo{$peptideID}[1],,$peptidesInfo{$peptideID}[5],$peptidesInfo{$peptideID}[0],,0,0;$proteinString\n";
					print FILE "q$queryNum","_p$i","_terms=$termsString\n";

					if ($peptidesInfo{$peptideID}[7]) { # Percolator data
						print FILE "q$queryNum","_p$i","_percolator=$peptidesInfo{$peptideID}[7]\n";
					}
					if ($peptidesInfo{$peptideID}[8]) { # PhosphoRS data
						print FILE "q$queryNum","_p$i","_phosphors=$peptidesInfo{$peptideID}[8]\n";
					}
				}
			}
			if ($i == 0) {
				print FILE "q$queryNum","_p1=-1\n";
			}
		}

		undef %peptidesInfo;
		undef %spectrumPeptides;

		print " Done.</B><BR>\n";

		print FILE "$boundary\n";
	}

}

####>Revision history<####
# 1.2.4 Minor modification for split mode file (GA 08/12/17)
# 1.2.3 Update SQLite queries to PD 2.2 version (GA 21/08/17)
# 1.2.2 Minor modification in PDM creation : add a \n in FixedMod (GA 05/05/17)
# 1.2.1 Minor modification to include TMT in PDM for 2.1 versions of Proteome Discoverer<BR>TODO: Check for lower versions (GA 31/03/17)
# 1.2.0 Update to convert PD 2.0+ "Any_N_Terminus" into "Any N-Terminus"<BR>TODO: Match peptide modifications on delta mass instead of name (PP 01/07/16)
# 1.1.9 Detects and moves UniProt isoform tag from ACC to ID if ID is selected (PP 05/02/15)
# 1.1.8 Checks for defined fasta file in &printHeaderMSF (PP 30/11/15)
# 1.1.7 Minor modification for Sequest pattern matching (GA 10/07/15)
# 1.1.6 Compatible with Percolator 64Bit and PD 2.0 (PP 16/04/15 + GA 25/06/15)
# 1.1.5b Check SQLite queries for 2.0 PD MSF (GA 12/06/15)
# 1.1.5a Compatible with Percolator 64Bit (PP 16/04/15)
# 1.1.4 Falls back to full protein identifier if parse rule fails (PP 23/03/15)
# 1.1.3 Add split-like conversion of MSF for each rawFiles of merge searches (GA 03/03/15)
# 1.1.2 Add PhosphoRS in "peptide" section (GA 07/01/15)
# 1.1.1 Minor modif to get Sequest HT or SEQUEST searches from PD1.4 (GA 28/11/14)
# 1.1.0 Bug fix to prevent double addition of fix mod mass to peptide mass (PP 13/10/14)
# 1.0.9 Bug fix to skip printing of internal Percolator data for non-matching peptides (PP 17/09/14)
# 1.0.8 Change query for fix modification otherwise this information is not saved in ANALYSIS_MODIFICATION (GA 01/09/14)
# 1.0.7 Minor modification in %percolatorMeasures retrieval (GA 22/05/14)
# 1.0.6 Fix bug for variable mods being initialized before /^Dyn/ check (PP 16/05/14)
# 1.0.5 Adds data from internal Percolator analysis (PP 12/05/14)
# 1.0.4 Restores former handling of fake decoy protein for version &le; 1.3 & improved match between quantification and search nodes<BR>pdm always generated in user's directory (PP 15/04/14)
# 1.0.3 decoy not imported if proteodiscover version == 1.2  (SL 27/04/14)
# 1.0.2 minor display bug fix (PP 05/04/14)
# 1.0.1 Better error handling (PP 22/01/14)
# 1.0.0 Extracted from storeAnalyses (PP 17/01/14)