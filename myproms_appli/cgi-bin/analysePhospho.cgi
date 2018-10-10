#!/usr/local/bin/perl -w

#############################################################################
# analysePhospho.cgi         1.1.3                                          #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                   #
# Contact: myproms@curie.fr                                                 #
# Script processing PhosphoRS analyses started by selectAnalyses.cgi        #
#############################################################################
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
use strict;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use phosphoRS;
#use IO::Uncompress::Unzip qw(unzip $UnzipError); # needed for &promsMod::extractSpectrumMSF
use XML::Simple; # needed for &promsMod::extractSpectrumMSF

my %promsPath=&promsConfig::getServerInfo;
my ($infoPepStrg, @infoPepList);
@infoPepList = map { "INFO_PEP$_" } (1..10);
$infoPepStrg = join(', ',@infoPepList);

#---------------------#
# Fetching parameters #
#---------------------#
my $probThreshold = param('probThreshold');
my $massTolerance = param('massTolerance');
my $activationType = param('activationType');
my $overWrite=(param('overwrite'))?param('overwrite'):0;
my @analysisList = param('anaList');

my $dbh=&promsConfig::dbConnect;
my ($phosphoModID)=$dbh->selectrow_array("SELECT ID_MODIFICATION FROM MODIFICATION WHERE UNIMOD_ACC=21");

#------------------#
# Only delete case #
#------------------#
if(param('deleteID')){
    my $branchID = param('branchID');
    print header(-'content-encoding' => 'no', -'charset' => 'UTF-8');
    print qq
|<HEAD>
<TITLE>Deleting PhosphoRS Analysis</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY>
|;
	&deletePRS($dbh,param('deleteID'));
	$dbh->disconnect;

	print qq
|<SCRIPT language="Javascript">
window.location="$promsPath{cgi}/selectAnalyses.cgi?ID=$branchID&callType=phosphoRS";
</SCRIPT>
</BODY>
|;
	exit;
}

$dbh->disconnect;
#------#
# HTML #
#------#
print header(-'content-encoding' => 'no', -'charset' => 'UTF-8');
warningsToBrowser(1);

print qq
|<HEAD>
<TITLE>Performing Phosphorylation Analyses</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT language="Javascript">
function newPhosphoRS() {
	var branchID=top.promsFrame.navFrame.getSelectedBranchID();
	window.location="./selectAnalyses.cgi?ID="+branchID+"&callType=phosphoRS";
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER><FONT class="title">Performing Phosphorylation Site Analysis</FONT></CENTER>
<BR>
<BR>
|;

#------------------#
# Looping analyses #
#------------------#
foreach my $anaID (@analysisList){

    # Fetching analysis data #
    $dbh=&promsConfig::dbConnect;
    my $startTime = time;

	my ($specificity)=$dbh->selectrow_array("SELECT SPECIFICITY FROM ANALYSIS_MODIFICATION WHERE ID_MODIFICATION=$phosphoModID AND ID_ANALYSIS=$anaID");
	$specificity=~s/,//g; # S,T,Y -> STY

    my ($countUnchanged,$countFlagged,$countChanged,$countConfirmed,$countTotal) = (0,0,0,0,0);
    my ($anaName,$dataFile,$fileFormat,$maxRank) = $dbh->selectrow_array("SELECT NAME,DATA_FILE,FILE_FORMAT,MAX_RANK FROM ANALYSIS WHERE ID_ANALYSIS=$anaID");
    $dbh->disconnect;
	my $anadir = "$promsPath{valid}/ana_$anaID";
    my $fullDataFileName = "$anadir/$dataFile";

    # Printing parameters file #
    #my @anaFilePath = split /\//, $fullDataFileName;
    my $fileName = $dataFile;
    $fileName =~ s/\.\w{3}$//;
    #my $anadir = join '/', @anaFilePath;
	my $paramFile = "$promsPath{valid}/ana_$anaID/PRSparam_ana_$anaID.txt";

    my $childPid=fork;
    if ($childPid == 0) { # child here
		$dbh=&promsConfig::dbConnect;
		print "<FONT class=\"title2\">+ PhosphoRS analysis of $anaName:</FONT></BR>\n";

		# Removing previous storing parameters file if exists #
		#my ($dataFileCutName) = ($dataFile =~ /^(.+)\.\w{3}$/);
		#my $paramFile = "$promsPath{valid}/ana_$anaID/PRSparam_$dataFileCutName.txt";
		&deletePRS($dbh,$anaID) if -e $paramFile;

		# Starting PhosphoRS #

		my $phosphoRS = new phosphoRS(AnaID => $anaID, File => $fullDataFileName, FileFormat => $fileFormat, MassTolerance => $massTolerance, ActivationType => $activationType);

		my %noAmbiguityQueries;
		if ($fullDataFileName =~ /([^\/]+)\.pdm$/){
			print "<FONT class=\"title3\">&nbsp;&nbsp;- Fetching spectrum list...";
			my $cutFileName = $1;
			$cutFileName =~ s/_\d$//;
			my $projectID = &promsMod::getProjectID($dbh, $anaID, 'ANALYSIS');
			my $dir = "$promsPath{valid}/multi_ana/proj_$projectID";
			my $msfFileName = "$dir/$cutFileName.msf";
			die "File: $dir/$cutFileName.msf not found!" unless -e $msfFileName;

			#my @infoPepList = map { "INFO_PEP$_ LIKE \"%Phospho%\"" } (1..10);
			#my $pepInfoQuery = join ' OR ', @infoPepList;
			#my $sthSpecId = $dbh->prepare("SELECT EXT_SPECTRUMID, QUERY_NUM FROM QUERY_VALIDATION WHERE VALID_STATUS>=? AND ID_ANALYSIS=? AND QUERY_NUM>0 AND ($pepInfoQuery)");
			my @spectrumList;
			my %sthPep;
			foreach my $i (1..$maxRank){
				$sthPep{$i} = $dbh->prepare("SELECT INFO_PEP$i FROM QUERY_VALIDATION WHERE ID_QUERY=?");
			}
			#my $sthSpecId = $dbh->prepare("SELECT EXT_SPECTRUMID,QUERY_NUM,V.ID_QUERY,GROUP_CONCAT(CONCAT(PEP_RANK,':',POS_STRING) SEPARATOR ',') FROM QUERY_VALIDATION V,QUERY_MODIFICATION M WHERE V.ID_QUERY=M.ID_QUERY AND VALID_STATUS>=-3 AND QUERY_NUM>0 AND ID_MODIFICATION=$phosphoModID AND ID_ANALYSIS=$anaID GROUP BY QUERY_NUM");
			my ($addQuery)=($overWrite)?'':' M.REF_POS_STRING IS NULL AND';###> Avoid to take into account modified queries
			my $sthSpecId = $dbh->prepare("SELECT EXT_SPECTRUMID,QUERY_NUM,V.ID_QUERY,GROUP_CONCAT(CONCAT(PEP_RANK,':',POS_STRING) SEPARATOR ',') FROM QUERY_VALIDATION V,QUERY_MODIFICATION M WHERE V.ID_QUERY=M.ID_QUERY AND$addQuery VALID_STATUS>=-3 AND QUERY_NUM>0 AND ID_MODIFICATION=$phosphoModID AND ID_ANALYSIS=$anaID GROUP BY QUERY_NUM");
			$sthSpecId->execute;
			while(my ($extSpectrumID,$queryNum,$queryID,$phosPosStrg) = $sthSpecId->fetchrow_array){
				#>Filetring for obvious unambiguous positions
				my $ambiguity=0;
				foreach my $pepPosData (split(',',$phosPosStrg)) {
					my ($rank,$posStrg)=split(':',$pepPosData);
					$sthPep{$rank}->execute($queryID);
					my ($infoPep)=$sthPep{$rank}->fetchrow_array;
					my ($seq) = ($infoPep=~ /SEQ=([^,]+),/);
					my @acceptors=$seq=~/[$specificity]/g;
					my @phosphoPos=split(/\./,$posStrg);
					if (scalar @acceptors > scalar @phosphoPos) { # more acceptors than modifs -> phosphoRS
						$ambiguity=1;
						last;
					}
				}
				if ($ambiguity) {
					push @spectrumList,[$queryNum,$extSpectrumID];
				}
				else { # no need for phosphRS
					$noAmbiguityQueries{$queryID}=$phosPosStrg;
				}
			}
			print " Done (",scalar keys %noAmbiguityQueries," unambiguous spectra will be skipped).</FONT><BR>\n";
			$sthSpecId->finish;
			foreach my $sth (values %sthPep){
				$sth->finish;
			}
			$dbh->disconnect;

			my $numSpectra=scalar @spectrumList;
			print "<FONT class=\"title3\">&nbsp;&nbsp;- Extracting $numSpectra spectra from MSF file: </FONT>";
			print "<B>0%" if $numSpectra >= 1000;
			my $dbsqlite = DBI->connect("dbi:SQLite:$msfFileName","","",{PrintError=>1,RaiseError=>1});
			my $c = 0;
			my $pc=0.1; # 1/10 fraction of 1
			foreach my $refSpectrum (@spectrumList) {
				$c++;
				my ($queryNum,$extSpectrumID)=@{$refSpectrum};
				my ($parentFile,@spectrumInfos) = &promsMod::extractSpectrumMSF($dbsqlite, $extSpectrumID); #, $dir, 'no_user');
				my $refIons = $spectrumInfos[0];
				my $charge = $spectrumInfos[4];
				my $peaks = '';

				foreach my $value (@{$refIons}){
					$peaks .= $value->{X} . ':' . $value->{Y} . ',';
				}
				$peaks =~ s/,$//;

				$phosphoRS->addSpectrumInfo(Peaks => $peaks, Charge => $charge, Activation => $activationType, Query => $queryNum) or die "Cannot add spectrum info for query $queryNum";
				print '<!--.-->'; # invisible (to keep connection)
				print '.' if (($c % 50) == 0); # 50 => 20 points between x%
				if ($numSpectra >= 1000 && $c/$numSpectra >= $pc) {
					print 100*$pc,'%';
					$pc+=0.1;
				}
			}
			$dbsqlite->disconnect;
			print "</B><FONT class=\"title3\"> Done.</FONT><BR>\n";
		}

		# Running PhosphoRS
		print "<FONT class=\"title3\">&nbsp;&nbsp;- Running PhosphoRS...";
		my $code = $phosphoRS->startAnalysis;
		#if($code != 0){
		#	print "<BR><FONT class=\"title2\" color=\"red\">Failed, code:$code</FONT><BR><BR>";
		#	$phosphoRS->cleanFiles;
		#	next;
		#}
		print " Done.</FONT><BR>\n";

		# Fetching queries from DB #
		print "<FONT class=\"title3\">&nbsp;&nbsp;- Updating phosphorylation sites position...";
		$dbh = &promsConfig::dbConnect;
		my %sthUpQ;
		foreach my $i (1..$maxRank){
			$sthUpQ{$i} = $dbh->prepare("UPDATE QUERY_VALIDATION SET INFO_PEP$i=? WHERE ID_QUERY=?");
		}
		my $sthUpQM1 = $dbh->prepare("UPDATE QUERY_MODIFICATION SET REF_POS_STRING=POS_STRING WHERE ID_MODIFICATION=$phosphoModID AND ID_QUERY=? AND PEP_RANK=? AND REF_POS_STRING IS NULL");
		my $sthUpQM2 = $dbh->prepare("UPDATE QUERY_MODIFICATION SET POS_STRING=? WHERE ID_MODIFICATION=$phosphoModID AND ID_QUERY=? AND PEP_RANK=?");
		my $sthQuery = $dbh->prepare("SELECT ID_QUERY,QUERY_NUM,$infoPepStrg FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$anaID AND QUERY_NUM>0 AND VALID_STATUS >=-3");
		my $sthPhosphoPosStg = $dbh->prepare("SELECT POS_STRING,REF_POS_STRING FROM QUERY_MODIFICATION WHERE ID_QUERY=? AND PEP_RANK=? AND ID_MODIFICATION=$phosphoModID");
		$sthQuery->execute;
		while (my ($idQuery,$queryNum,@infoPep) = $sthQuery->fetchrow_array) {
			if ($noAmbiguityQueries{$idQuery}) {
				foreach my $pepPosData (split(',',$noAmbiguityQueries{$idQuery})) {
					my ($rank,$posStrg)=split(':',$pepPosData);
					my $infoPepStrg = $infoPep[$rank-1];
					$infoPepStrg .= "PRS=3;100;,";
					$sthUpQ{$rank}->execute($infoPepStrg,$idQuery);
					$countConfirmed++;
					$countTotal++;
				}
				next;
			}
			# Browsing each info pep #
			my $i=0;
			while( my $infoPepStrg = shift @infoPep){
				$i++;
				$sthPhosphoPosStg->execute($idQuery,$i);
				my ($positionStrg,$refpos)=$sthPhosphoPosStg->fetchrow_array;
				$positionStrg=$refpos if ($refpos && $positionStrg =~ /-\d/);# one site is ambigous
				if($infoPepStrg =~ /VMOD=([^,]+),/){# keep VMOD for RMOD
					my $vmodStrg = $1;
					my $rmodStrg;
					if ($infoPepStrg !~ /RMOD=/) {
						$rmodStrg = $vmodStrg;
					}
					#my $positionStrg;
					# Fetching all phosphorylation positions #
					my %positions;
					#while($vmodStrg =~ /Phospho \(([A-Z]+):([^\)]+)\)/g){
					#	$positionStrg .= ($positionStrg)? ".$2" : $2;
					#	$positions{$1} = $2;
					#}
					$positions{'STY'}=$positionStrg;
					next unless $positionStrg;
					my ($seq) = ($infoPepStrg=~ /SEQ=([^,]+),/);
					#next unless ($seq =~ /[STY].*[STY]/);makePhosphoRS.sh

					$countTotal++;

					# Fetching best isoform from PhosphoRS #
					#print "Query: $queryNum , i: $i , Positions: $positionStrg<BR>\n";
					if(my $isoformRef = $phosphoRS->getIsoforms($queryNum,$i)){
						my @isoforms = @{$isoformRef};
						my $bestIsoformRef = (sort {$b->[0] <=> $a->[0]} @isoforms )[0];
						my ($proba,$positionsRef) = @{$bestIsoformRef};
						$proba *= 100;
						my $bestPositions = join('.',@{$positionsRef});
						$infoPepStrg =~ s/PRS=[^,]+,//g;

						# Updating info pep data strings #
						if($bestPositions ne $positionStrg){
							if($proba >= $probThreshold){
								# Change Mascot positions to phosphoRS positions
								my ($sequence) = ($infoPepStrg =~ /SEQ=([^,]+)/);
								my @sequence = split(//,$sequence);
								my %newPositions;
								foreach my $position (@{$positionsRef}){
									my $aa = $sequence[$position-1];
									#my $phosphoType = ($aa eq 'S' or $aa eq 'T')? 'ST' : $aa; # assuming that S & T phosphorylations are always regrouped
									my $phosphoType = ($aa=~/[$specificity]/)? $specificity : $aa; # uses search specif (PP 19/12/14)
									push @{$newPositions{$phosphoType}}, $position;
								}
								$vmodStrg =~ s/ \+ Phospho[^\)]+\)//g;
								foreach my $phosphoType (keys %newPositions){
									my $posStrg = join '.', @{$newPositions{$phosphoType}};
									$vmodStrg .= " + Phospho ($phosphoType:$posStrg)";
								}
								$infoPepStrg =~ s/VMOD=[^,]+/VMOD=$vmodStrg/;
								# Store old positions
								$infoPepStrg .= "RMOD=$rmodStrg," if $rmodStrg;
								$sthUpQM1->execute($idQuery,$i); # Record original position
								$sthUpQM2->execute($bestPositions,$idQuery,$i);
								my $aaModPos;
								foreach my $aa ( keys %positions ){ # to keep phospho aa type (lost by PhosphoRS)
									$aaModPos .= "[$aa]";
									$aaModPos .= $positions{$aa};
								}
								$infoPepStrg .= "PRS=2;$proba;$aaModPos,"; # (status ; PRS score ; previous Mascot positions)
								$countChanged++;
							}
							else {
								# PRS != Mascot, but unchanged because PRS proba < threshold
								$infoPepStrg .= "PRS=1;$proba;$bestPositions,"; # (status ; PRS score ; PRS best positions)
								$countFlagged++;
							}
						}
						else {
							# PRS = Mascot
							if($proba >= $probThreshold){
								$infoPepStrg .= "PRS=3;$proba;,"; # same isoform for PRS and Mascot
								$countConfirmed++;
							}
							else {
								$infoPepStrg .= "PRS=0;$proba;,"; # same isoform but PRS score < threshold
								$countUnchanged++;
							}
						}
					}
					else {
					   $infoPepStrg .= "PRS=4;;,";
					}
					$sthUpQ{$i}->execute($infoPepStrg,$idQuery);
				}
			}
		}
		foreach my $sth (values %sthUpQ){
			$sth->finish;
		}
	    $sthUpQM1->finish;
	    $sthUpQM2->finish;
		$sthPhosphoPosStg->finish;
	    $dbh->commit;
	    $dbh->disconnect;
	    open (INFILE, ">$anadir/PRS_${anaID}_done.txt");
	    print INFILE "$countConfirmed,$countUnchanged,$countChanged,$countFlagged,$countTotal\n";
	    close(INFILE);
	    exit;
    }
    else {
		while (! -e "$anadir/PRS_${anaID}_done.txt") {
			sleep 60;
			print "."
		}
		open (INFILE, "$anadir/PRS_${anaID}_done.txt");
		my $line=<INFILE>;
		chomp($line);
		($countConfirmed,$countUnchanged,$countChanged,$countFlagged,$countTotal)=split(/,/,$line);
		close(INFILE);
		unlink("$anadir/PRS_${anaID}_done.txt");

		## Printing parameters file #
		#my @anaFilePath = split /\//, $fullDataFileName;
		#my $fileName = pop @anaFilePath;
		#$fileName =~ s/\.\w{3}$//;
		#my $dir = join '/', @anaFilePath;
		#open PRSPARAM, ">$anadir/PRSparam_$fileName.txt";
		open (PRSPARAM, ">$paramFile");
		print PRSPARAM qq
|Threshold:$probThreshold%
Mass Tolerance:$massTolerance Da
Activation Type:$activationType|;
		close PRSPARAM;
		print " Done.</FONT><BR>\n";

		# Summary display #
		my $processTime = time - $startTime;
		my $processTimeStrg = int($processTime/60).' min '.($processTime % 60).' sec';
		print qq
|<FONT class="title3">&nbsp;&nbsp;- Analysis completed in $processTimeStrg:<BR>
<FONT class="title3">
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$countConfirmed/$countTotal phosphosites confirmed by PhosphoRS (probability > $probThreshold%)<BR>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$countUnchanged/$countTotal phosphosites confirmed by PhosphoRS (probability < $probThreshold%)<BR>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$countChanged/$countTotal phosphosites were changed<BR>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$countFlagged/$countTotal phosphosites were differently determinated by PhosphoRS but unchanged (probability < $probThreshold%)<BR>
</FONT><BR>
|;
		$dbh = &promsConfig::dbConnect;
		# Updating history #
		my $paramStrg = "threshold:$probThreshold;massTolerance:$massTolerance;activationType:$activationType;";
		&promsMod::updateAnalysisHistory($dbh, $anaID, $paramStrg, 'prs');
		$dbh->commit;
		$dbh->disconnect;
    }
}

print qq
|<CENTER><INPUT type="button" value="New PhosphoRS Analysis" onclick="newPhosphoRS();"></CENTER>
</BODY>
</HTML>
|;

sub deletePRS{
    my ($dbh,$anaID) = @_;

	print "<FONT class=\"title3\">&nbsp;&nbsp;- Deleting previous PhosphoRS analysis...";
    ## Restore DB entries ##
	my ($maxRank) = $dbh->selectrow_array("SELECT MAX_RANK FROM ANALYSIS WHERE ID_ANALYSIS=$anaID");
	my @infoPepList = map { "INFO_PEP$_" } (1..$maxRank);
	my $infoPepStrg = join(', ',@infoPepList);

	# table QUERY_MODIFICATION #
	$dbh->do("UPDATE QUERY_MODIFICATION M INNER JOIN QUERY_VALIDATION V ON M.ID_QUERY=V.ID_QUERY AND V.ID_ANALYSIS=$anaID AND M.ID_MODIFICATION=$phosphoModID AND M.REF_POS_STRING IS NOT NULL SET M.POS_STRING=M.REF_POS_STRING");
	$dbh->do("UPDATE QUERY_MODIFICATION M INNER JOIN QUERY_VALIDATION V ON M.ID_QUERY=V.ID_QUERY SET M.REF_POS_STRING=NULL WHERE V.ID_ANALYSIS=$anaID AND M.ID_MODIFICATION=$phosphoModID");

	# table QUERY_VALIDATION #
    my $sthInfoPep = $dbh->prepare("SELECT ID_QUERY,$infoPepStrg FROM QUERY_VALIDATION WHERE ID_ANALYSIS=? AND QUERY_NUM>0");
    my %sthUpInfoPep;
    for(my $i=1; $i<=$maxRank; $i++){
		$sthUpInfoPep{$i} = $dbh->prepare("UPDATE QUERY_VALIDATION SET INFO_PEP$i=? WHERE ID_QUERY=?");
    }
    $sthInfoPep->execute($anaID);
    while(my ($idQuery,@infoPep) = $sthInfoPep->fetchrow_array){
		my $i = 0;
		while(my $infoPep = shift @infoPep){
			$i++;
			if($infoPep =~ /VMOD=([^,]+Phospho[^,]+),/){
				my $vmodStrg = $1;
				if($infoPep =~ /PRS=(\d);/){
					if($1 == 2){
						# Restore old positions
						my ($aaModPos) = ($infoPep =~ /PRS=\d;[^;]+;([^,]+)/);
						my $restoredPhosphoStrg = '';
						while($aaModPos =~ /\[(\w+)\]([^\[,]+)/g){
							$restoredPhosphoStrg .= " + Phospho ($1:$2)";
						}
						$vmodStrg =~ s/ \+ Phospho \([^\)]+\)//g;
						$vmodStrg .= "$restoredPhosphoStrg";

						$infoPep =~ s/VMOD=[^,]+/VMOD=$vmodStrg/;
						$infoPep =~ s/RMOD=[^,]+,//; # Deleting PRS data
					}
					$infoPep =~ s/PRS=[^,]+,//; # Deleting PRS data
					$sthUpInfoPep{$i}->execute($infoPep,$idQuery);
				}
			}
		}
    }
    $sthInfoPep->finish;
    foreach my $i (keys %sthUpInfoPep){
		$sthUpInfoPep{$i}->finish;
    }

    # Update history #
    $dbh->do("DELETE FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$anaID AND VAL_TYPE='prs' AND STATUS<1");
    $dbh->do("UPDATE VALIDATION_HISTORY SET STATUS=-1 WHERE ID_ANALYSIS=$anaID AND VAL_TYPE='prs' AND STATUS=1");
    $dbh->commit;

    # Delete PRS files #
    #my ($dataFileName) = $dbh->selectrow_array("SELECT DATA_FILE FROM ANALYSIS WHERE ID_ANALYSIS=$anaID");
    #$dataFileName =~ s/\.\w{3}//;
    my $dir = "$promsPath{valid}/ana_$anaID";
    #unlink "$dir/PRS_$dataFileName.xml" or warn "Cannot delete PRS file: $!";
    #unlink "$dir/PRSparam_$dataFileName.txt" or warn "Cannot delete PRS file: $!";
	#unlink "$dir/PRS_ana_$anaID.xml" or warn "Cannot delete PRS output file: $!";
	#unlink "$dir/PRSparam_ana_$anaID.txt" or warn "Cannot delete PRS parameter file: $!";
	unlink glob "$dir/PRSparam_*.txt"; # old & new param naming
	unlink glob "$dir/PRS_*"; # old & new result.xml naming + status.txt
	unlink "$dir/PRSerrors.txt" if -e "$dir/PRSerrors.txt";
	print " Done.</FONT><BR>\n";
}

####>Revision history<####
# 1.1.3 Minor bug fix in call for a new PhosphoRS analysis (PP 15/06/18)
# 1.1.2 Modification to avoid positionStg stored in INFO_PEP (GA 31/01/18)
# 1.1.1 Uses xxx_ana_&lt;anaID&gt;.xxx instead of data file name for PRS parameter and output files (PP 21/03/17)
# 1.1.0 Add a fork for big xml files that take too long with XMLin in phosphoRS.pm (GA 11/03/16)
# 1.0.9 Minor change (PP 26/08/15)
# 1.0.8 Compatibility update du to new output of &promsMod::extractSpectrumMSF v3.3.9 (PP 18/02/15)
# 1.0.7 Full compatibily with MSF files, updates table QUERY_MODIFICATION & speed improvement (PP 23/12/14)
# 1.0.6 Added 'use XML::Simple' (PP 13/11/13)
# 1.0.5 Added file format to phosphoRS parameters (FY 30/09/13)
# 1.0.4 Using RMOD tag in INFO_PEP string to store previous positions (FY 11/04/13)
# 1.0.3 New status (4) if no phospho-isoform for a given peptide (FY 06/11/12)
# 1.0.2 Management of any kind of phosphorylation (FY 26/07/12)
# 1.0.1 Management of Proteome Discoverer files +<BR>Add phosphopeptides filtering (FY 24/05/12)
# 1.0.0 New script processing PhosphoRS analyses (FY 28/03/12)