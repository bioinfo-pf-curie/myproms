#!/usr/local/bin/perl -w

################################################################################
# runemPAIQuantification.pl   1.1.9                                            #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Allows quantification of proteins by retrieving the emPAI of the proteins    #
# that were computed by Mascot (need to configure Mascot web-server)           #
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
use promsConfig;
use promsMod;
use POSIX qw(strftime);
use strict;
use LWP::UserAgent;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo('no_user');
my %mascotServers=&promsConfig::getMascotServers;
my $mascotServer;
foreach my $server (keys %mascotServers){
	$mascotServer=$server;
	last; #take first one in list (hopefully only 1 mascot server declared)
}

###############################
####>Recovering parameters<####
###############################
my ($quantiID,$quantifDate,$quantItemID)=@ARGV;
my ($analysisID,$qID)=split(/\./,$quantItemID);
my $quantifDir="$promsPath{tmp}/quantification/$quantifDate";
my $fileStat= "$quantifDir/status_$quantItemID.out";
open(FILESTAT,">$fileStat") || die "Error while opening $fileStat";
#print FILESTAT "#@ARGV#\n";
#print strftime("%H:%M:%S %d/%m/%Y",localtime);
print FILESTAT "Started ",strftime("%H:%M:%S %d/%m/%Y",localtime),"\n";
###########################################################
####>Check if the analysis was done with Mascot Server<####
###########################################################
my (%emPAI,%protIDs,%mr);
my $dbh=&promsConfig::dbConnect('no_user');
my ($fileFormat,$datafile)=$dbh->selectrow_array("SELECT FILE_FORMAT,DATA_FILE FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
my ($projectID)=&promsMod::getProjectID($dbh,$analysisID,'analysis');
$dbh->disconnect;

if ($fileFormat =~ /MASCOT/) {
    my ($datFile,$oldName)=("$promsPath{'valid'}/ana_$analysisID/$datafile",'');
    if (!-e "$promsPath{'valid'}/ana_$analysisID/$datafile") {# Validation has been finished and the dat file has been removed... (anterior version)
		my $pepFileExt=($fileFormat=~/\.PDM/)? 'pdm' : 'dat';
		my $oldPepFileName=sprintf "P%06d.$pepFileExt",$analysisID;
		$datFile="$promsPath{peptide}/proj_$projectID/$oldPepFileName";
    }
    open (FILE, $datFile) || (print FILESTAT "Unable to open $datFile file\n");
    my ($section,$dateCode)=('','');
    while (my $line=<FILE>) {
        $section=$1 if $line=~/name="(\w+)"/;
        last if ($section eq 'summary');
        next if ($section && $section ne 'parameters' && $section ne 'header');
        $line=~s/\s*\Z//; # delete \n and trailing spaces
        if ($fileFormat =~ /\.PDM/ && $line =~ /^COM=(.+)/) {
			$oldName=$datafile;
			$datafile=$1;#the name of the original DATFILE is put at this point
        }
		elsif ($line=~ /^date=(\d+)/ && $fileFormat !~ /\.PDM/){
			###> Only for mascot dat. Conversion from seconds to a datecode!
			###> localtime vs gmtime : should be the same af mascot server so if mascot is gmtime, you should use gmtime 
			#$dateCode=strftime("%Y%m%d",localtime($1));
			$dateCode=strftime("%Y%m%d",gmtime($1));
        }
    }
    print FILESTAT "oldname=$oldName datafile=$datafile datecode$dateCode datFile=$datFile $fileFormat\n";
    close FILE;
    ###>Need to get the proper dateCode
    if($fileFormat=~ /\.PDM/){
		$oldName =~ s/\_[0-9]*\.pdm/\.ana/g;
		(my $msfFile = $oldName) =~ s/\.ana/\.msf/;
		my $dbsqlite = DBI->connect( "dbi:SQLite:$promsPath{valid}/multi_ana/proj_$projectID/$msfFile", "", "", {PrintError => 1,RaiseError => 1});
		###> Be careful, the commented code won't work for DAT files generated overnight (the MSF was created before the production of the DAT therefore the datecode is incorrect)
		#my ($proteomeDiscovererVersion)=$dbsqlite->selectrow_array("SELECT SoftwareVersion FROM SchemaInfo ORDER BY rowid ASC LIMIT 1"); # last record
		#$proteomeDiscovererVersion=~s/^(\d\.\d+).*/$1/; # x.x.x.xx -> x.x (numerical value)
		#my ($wftable)=($proteomeDiscovererVersion < 2.0)?'WorkflowInfo':'Workflows';
		#my ($date)=$dbsqlite->selectrow_array("SELECT WorkflowStartDate FROM $wftable");
		#my ($dateMSF,$startTime)=split(/\s/,$date);
		#($dateCode=$dateMSF) =~ s/-//g;
		my ($message)=$dbsqlite->selectrow_array("SELECT Message FROM WorkflowMessages WHERE Message LIKE 'Received Mascot result file %'");
		$dbsqlite->disconnect;
		if ($message =~ /filename=..\/data\/(\d+)\//) {	$dateCode=$1; }
    }

    ###>Get the identifiers of the validated proteins:
    print FILESTAT "1/3 Getting protein data from myProMS database\n";
    close FILESTAT;
    $dbh=&promsConfig::dbConnect('no_user');
    my ($numDatabanks)=$dbh->selectrow_array("SELECT COUNT(*) FROM ANALYSIS_DATABANK WHERE ID_ANALYSIS=$analysisID");
    #my $sthProtValid=$dbh->prepare("SELECT PROTEIN.ID_PROTEIN,IDENTIFIER FROM PROTEIN,ANALYSIS_PROTEIN WHERE ANALYSIS_PROTEIN.ID_PROTEIN=PROTEIN.ID_PROTEIN AND ID_ANALYSIS=$analysisID");
    #my $sthProtValid=$dbh->prepare("SELECT PROTEIN.ID_PROTEIN,IDENTIFIER,DB_RANK FROM PROTEIN,ANALYSIS_PROTEIN WHERE ANALYSIS_PROTEIN.ID_PROTEIN=PROTEIN.ID_PROTEIN AND ID_ANALYSIS=$analysisID");
    my $sthProtValid=$dbh->prepare("SELECT PROTEIN.ID_PROTEIN,ID_MASTER_PROTEIN,IDENTIFIER,DB_RANK,MW FROM PROTEIN,ANALYSIS_PROTEIN WHERE ANALYSIS_PROTEIN.ID_PROTEIN=PROTEIN.ID_PROTEIN AND ID_ANALYSIS=$analysisID");
    my $sthGetWeight=$dbh->prepare("SELECT MW FROM MASTER_PROTEIN WHERE ID_MASTER_PROTEIN=?");
    $sthProtValid->execute;
    while (my ($idProtein,$masterProtID,$identifier,$dbRank,$mw) = $sthProtValid->fetchrow_array) {
	    my $identifierStrg=($numDatabanks>1)? $dbRank.'::'.$identifier : $identifier;
	    $protIDs{$identifierStrg}=$idProtein;
	    if ($mw) {
		    $mr{$identifierStrg}=$mw;
	    }elsif($masterProtID){
		    $sthGetWeight->execute($masterProtID);
		    ($mr{$identifierStrg})=$sthGetWeight->fetchrow_array;
	    }else{
		    $mr{$identifierStrg}=0;
	    }
    }
    $sthProtValid->finish;
	$sthGetWeight->finish;
    $dbh->disconnect;

    open(FILESTAT,">>$fileStat");
    #print FILESTAT "2/3 Getting emPAI from Mascot web-server ($nbRetProts emPAI retrieved among $nbProts)\n";
    print FILESTAT "2/3 Getting emPAI from Mascot web-server\n";
    close FILESTAT;
    my ($sumemPAI,$sumemPAImr)=&getMascotemPAI($datafile,$dateCode,$analysisID,\%emPAI,\%mr);

    ###>Create new Quantitation in myProMS database.
    open(FILESTAT,">>$fileStat");
    print FILESTAT "3/3 Storing emPAI into myProMS database\n";
    close FILESTAT;

    $dbh=&promsConfig::dbConnect('no_user');

    my %empaiPARAMS;
    my $sthGetParamIDs=$dbh->prepare("SELECT CODE,ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE LIKE 'EMPAI%'");#In QUANTIFICATION_PARAMETER table -> emPAI value
    $sthGetParamIDs->execute();
    while (my ($code,$idQP) = $sthGetParamIDs->fetchrow_array) {
	    $empaiPARAMS{$code}=$idQP;
    }
    
    my $sthProtQuantif=$dbh->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_PROTEIN,ID_QUANTIFICATION,ID_QUANTIF_PARAMETER,QUANTIF_VALUE) VALUES (?,$quantiID,?,?)");
    foreach my $protID (keys(%emPAI)){
	    if ($protIDs{$protID} && $emPAI{$protID}{$analysisID}){ #Check if the protein was validated for this Analysis
		    $sthProtQuantif->execute($protIDs{$protID},$empaiPARAMS{'EMPAI'},$emPAI{$protID}{$analysisID});
		    $sthProtQuantif->execute($protIDs{$protID},$empaiPARAMS{'EMPAI_MOL'},$emPAI{$protID}{$analysisID}*100.00/$sumemPAI) unless !$empaiPARAMS{'EMPAI_MOL'} or $sumemPAI == 0;
		    $sthProtQuantif->execute($protIDs{$protID},$empaiPARAMS{'EMPAI_MR'},$emPAI{$protID}{$analysisID}*$mr{$protID}*100.00/$sumemPAImr) unless !$empaiPARAMS{'EMPAI_MR'} or $sumemPAImr == 0;
	    }
    }
    $sthProtQuantif->finish;

    $dbh->do("UPDATE QUANTIFICATION SET STATUS=1 WHERE ID_QUANTIFICATION=$quantiID") || die $dbh->errstr();

    $dbh->commit;
    $dbh->disconnect;

    open(FILESTAT,">>$fileStat");
    print FILESTAT "Quantification Ended";
    close(FILESTAT);

    sleep 2;
}
else {
    open(FILESTAT,">>$fileStat");
    print FILESTAT "Quantification Ended\n";
    print FILESTAT "WARNING : emPAI quantification is only possible to searches performed on Mascot...";
    close(FILESTAT);
}

##############################################################################################
######> Get Mascot emPAI using the myproms4emPAI.pl script in cgi <########################
##############################################################################################
sub getMascotemPAI {

	my($datFile,$datCalendar,$analysisID,$refemPAI,$refmr)=@_;
	my ($sumemPAI,$sumemPAImr)=(0.00,0.00);
    
	### > Before 2.5.1
	#my $params = {
	#	file => "../data/$datCalendar/$datFile",
	#	do_export => 1,
	#	prot_hit_num => 1,
	#	prot_acc => 1,
	#	pep_query => 1,
	#	pep_rank => 1,
	#	pep_isunique => 1,
	#	pep_exp_mz => 1,
	#	export_format => 'CSV',
	#	_sigthreshold => 0.05,
	#	report => 'AUTO',
	#	_server_mudpit_switch => 0.000000001,
	#	show_same_sets => 1,
	#	_showsubsets => 1,
	#	protein_master => 1,
	#	prot_empai => 1,
	#	peptide_master => 1,
	#	pep_seq => 1
	#};
	my $params = {
		file => "../data/$datCalendar/$datFile",
		_minpeplen => 5,
		_showsubsets => 1,
		_server_mudpit_switch => 0.000000001,
		_sigthreshold => 0.05,
		export_format => 'CSV',
		do_export => 1,
		group_family => 1,
		pep_exp_mz => 1,
		pep_isbold => 1,
		pep_isunique => 1,
		pep_query => 1,
		pep_rank => 1,
		peptide_master => 1,
		prot_acc => 1,
		prot_empai => 1,
		prot_hit_num => 1,
		protein_master => 1,
		sessionid => 'all_secdisabledsession',
		use_homology => 1
	};
    
	my $agent = LWP::UserAgent->new;
	$agent->timeout(360);
	push @{$agent->requests_redirectable}, 'POST';
	if ($mascotServers{$mascotServer}{proxy}) { # proxy settings
		if ($mascotServers{$mascotServer}{proxy} eq 'no') {$agent->no_proxy($mascotServers{$mascotServer}{url});}
		else {$agent->proxy('http', $mascotServers{$mascotServer}{proxy});}
	} else {$agent->env_proxy;}
	#my $response = $agent->post("$mascotServers{$mascotServer}{url}/cgi/myproms4emPAI_v2.pl",$params);
	my $response = $agent->post("$mascotServers{$mascotServer}{url}/cgi/myproms4emPAI.pl",$params);
    
	while (my $wait = $response->header('Retry-After')) {
		#print "Waiting ($wait)...\n";
		sleep $wait;
		$response = $agent->get($response->base);
	}
    
	###> 3rd: read the response of the request (if success)
	if ($response->is_success) {
		#print $response->content;
		my @lines=split(/\n/,$response->content);
		foreach my $line ( @lines ) {
			next unless $line =~ /emPAI/;
			#my @info = split(/,/, $line);# prot_hit_num, prot_acc, pep_query, pep_rank, pep_isbold, pep_isunique, pep_exp_mz, pep_res_before , pep_seq , pep_res_after , emPAI
			my @info = split(/,/, $line); #prot_hit_num,prot_family_member,prot_acc,pep_query,pep_rank,pep_isbold,pep_isunique,pep_exp_mz
			$info[2] =~ s/"//g;
			$refemPAI->{$info[2]}{$analysisID}=$info[9];
			if ($refmr->{$info[2]}) { # Some proteins found in DAT file could have been removed in validation-step. They should not be used for normalization
				$sumemPAI+=$info[9];
				$sumemPAImr+=$info[9]*$refmr->{$info[2]};
			}
		}
	}
	else{
		print $response->status_line, "<BR>\n";
	}
	return ($sumemPAI,$sumemPAImr);
}

###> Obsolete -> 24/12/2012
##############################################################################################
######> Alternative to the precedent thing                                                <###
######> It gets the emPAI (when it is available) from the valid list of proteins directly <###
##############################################################################################
###sub getMascotemPAI {
###
###    my($datFile,$datCalendar,$identifier,$analysisID,$maxNbProt,$refemPAI)=@_;
###
###    ###> 1st : writing the information sent for wget
###    my $postData="..%2Fdata%2F$datCalendar%2F$datFile;";
###    $postData.="_ignoreionsscorebelow=0;_proteinfamilyswitch=0;_sigthreshold=0.05;_sortunassigned=scoredown;";
###    $postData.="mimetype=application%2Fvnd.matrixscience.reportnode%2Bxml;";
###    if ($identifier eq 'getAll'){
###	### pr.show=quantitation option allows to show in Protein Family mode the quantitation tab with all the proteins
###	$postData.="pr.show=quantitation;pr.ftype=accession;pr.per_page=0;report=$maxNbProt;sub=tc%3Arf";
###    }else{
###	$identifier =~ s/\|/%7C/g;
###	$postData.="pr.find=$identifier;pr.ftype=accession;pr.per_page=0;report=$maxNbProt;sub=tc%3Arf";
###    }
###
###    ###> 2nd : send the request to the agent
###    my $ua = LWP::UserAgent->new;
###    $ua->timeout(360);
###    my $res = $ua->get("$mascotServer/cgi/master_results_2.pl?file=$postData");
###
###    ###> Block to prevent from problems with the server
###    while (my $wait = $res->header('Retry-After')) {
###            #print "Waiting ($wait)...\n";
###            sleep $wait;
###            $res = $ua->request($res->base);
###    }
###
###    ###> 3rd: read the response of the request (if success)
###    if ($res->is_success) {
###	#print $res->content;
###	my ($onGi,$onForm,$myId)=(0,0);
###	my ($prot);
###	my @lines=split(/\n/,$res->content);
###	foreach my $line ( @lines ) {
###	    if ($identifier eq 'getAll') {# For merge raw files, the second part did not work (no emPAI information given)
###		if($onForm==1){
###		    if ($line =~/<td class="empai">([^<]+)<\/td>/) {
###			$refemPAI->{$prot}{$analysisID}=$1;
###		    }elsif($line =~/<\/tr>/){
###			$onForm=0;
###		    }
###		}elsif($line =~ /<td class="accession"/){
###		    $onForm=1;
###		    if($line=~/>([^>]+)<\/a><\/td>/){# Works for SwissProt and NCBI accession number
###			$prot=$1;
###		    }
###		}
###	    }else {
###		if($onForm==1){
###		    if($line =~/title="Detailed view of ([^>]+)">/){# Works for SwissProt and NCBI accession number
###			$onGi=1;
###			$prot= $1;#TODO: Consider other accession formats!!!
###			#print " $prot";
###		    }elsif($onGi=1 && $line =~/<td class="empai">([^<]+)<\/td>/){
###			$refemPAI->{$prot}{$analysisID}=$1;
###			#print " empai=$1 <BR>\n";
###		    }elsif($line =~/<\/tr>/){
###			$onGi=0;
###		    }elsif($line=~/<\/form>/){
###			last;
###		    }
###		}else{
###		    if($line =~/<th class="empai" scope="col"/ || $line =~/<th scope="col" class="empai"/){
###			$onForm=1;#the emPAI value is supposed to be on this form of the file
###		    }
###		}
###	    }
###	}
###    }else{
###	print $res->status_line, "<BR>\n";
###    }
###}

####>Revision history<####
# 1.1.9 [MODIFICATION] Small modifs on log printing (VS 09/10/19) 
# 1.1.8 Uses new &promsConfig::getMascotServers function (PP 25/06/19)
# 1.1.7 Minor fix to prevent DBI warning (PP 30/10/18)
# 1.1.6 Add a test for emPAI normalization and change datecode retrievement for Mascot DAT searches -> use gmtime and not localtime (GA 24/02/16)
# 1.1.5 Minor change script version -> myproms4emPAI_v2.pl to myproms4emPAI.pl (GA 11/01/16)
# 1.1.4 Update for new Mascot web-server 2.5.1 version + change datecode retrieval for PDM searches (GA 17/07/15)
# 1.1.3 Minor modification for proteins without weight information -> get the master protein one (GA 10/06/14)
# 1.1.2 Include computing of molar and weight percentage normalization (GA 28/05/14)
# 1.1.1 Change args of perl file (GA 16/04/14)
# 1.1.0 Fix multi-databank bug (GA 08/04/14)
# 1.0.9 Rename export_da2_ter.pl into myproms4emPAI.pl (GA 19/09/13)
# 1.0.8 Remove ID_ANALYSIS value from insert in PROTEIN_QUANTIFICATION & use autoincrement (PP 12/09/13)
# 1.0.7 Better proxy declaration for LWP (PP 02/07/13)
# 1.0.6 Minor modification to add subsets in the CSV file (GA 21/06/13)
# 1.0.5 Modification of the script to use LWP from Mascot web-server (GA 24/12/12)
# 1.0.4 Minor modification: replace wget by lwp command (GA 12/10/12)
# 1.0.3 Modification of query-order sent to Mascot (GA 11/10/12)
# 1.0.2 Modification of query sent to Mascot so as to work for merge files (GA 04/09/12)
# 1.0.1 Minor update: add ID_ANALYSIS in the insert of PROTEIN_QUANTIFICATION (GA 07/07/2011)
# 1.0.0 New script to compute emPAI quantification