#!/usr/local/bin/perl -w

################################################################################
# runPathwayAnalysis.pl       1.0.0                                            #
# Authors: P. Poullet, S.Liva (Institut Curie)                                 #
# Contact: myproms@curie.fr                                                    #
# Launches Pathway Analysisis                                                  #
# called by startPathwayAnalysis.cgi                                           #
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
use File::Path qw(make_path remove_tree);
use File::Copy qw(move);
use promsConfig;
use promsMod;
use LWP::UserAgent;
use LWP::Simple qw(getstore);
use JSON;
my %promsPath=&promsConfig::getServerInfo;
my $dbh = &promsConfig::dbConnect;
my %protIDs;

my ($pathwayID, $catID)=@ARGV;
my $projectID=&promsMod::getProjectID($dbh,$pathwayID,'PATHWAY_ANALYSIS');
my $projectPath="$promsPath{pathAna}/project_$projectID";

my $tmpPathwayPath="$promsPath{tmp}/pathway_analysis";
my $tmpPathwayPathID="$tmpPathwayPath/$pathwayID";

make_path("$projectPath", { verbose => 0, mode => 0755}) unless -e "$projectPath";
make_path("$tmpPathwayPath", { verbose => 0, mode => 0755}) unless -e "$tmpPathwayPath";
make_path("$tmpPathwayPathID", { verbose => 0, mode => 0755}) unless -e "$tmpPathwayPathID";

my $codeIdent="AC";
my ($identifierID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='$codeIdent'");
my ($strgFromSQL, $strgWhereSQL, $strgCondSQL)=($catID)? ("CATEGORY_PROTEIN C,","C.ID_CATEGORY=? and C.ID_PROTEIN=P.ID_PROTEIN and","") : ("PATHWAYANA_ANALYSIS PA,ANALYSIS_PROTEIN AP,","PA.ID_ANALYSIS = AP.ID_ANALYSIS and AP.ID_PROTEIN = P.ID_PROTEIN and", "and PA.ID_PATHWAY_ANALYSIS=? and AP.VISIBILITY>=1");
my $sthProt=$dbh->prepare("SELECT P.ID_PROTEIN, MI.VALUE  FROM $strgFromSQL PROTEIN P, MASTER_PROTEIN MP, MASTERPROT_IDENTIFIER MI
			  where $strgWhereSQL
			  P.ID_MASTER_PROTEIN = MP.ID_MASTER_PROTEIN and
			  MP.ID_MASTER_PROTEIN = MI.ID_MASTER_PROTEIN and
			  MI.ID_IDENTIFIER=? and MI.RANK=1 $strgCondSQL");
($catID)? $sthProt->execute($catID, $identifierID) : $sthProt->execute($identifierID, $pathwayID);
while(my ($protID, $uniprot) = $sthProt->fetchrow_array){
	$protIDs{$uniprot}{$protID}=1;
}
$sthProt->finish;
my $listRequest = join("\n", keys %protIDs);
#print $listRequest."<br>";
my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
$agent->timeout(10);
$agent->env_proxy;
$agent->show_progress(1);
my $resource="UNIPROT";

##Give a list of protein to reactome web service
my $responseFormPost = $agent->post('http://www.reactome.org:80/AnalysisService/identifiers/form?pageSize=200&page=1&sortBy=NAME&order=ASC&resource=UNIPROT',
				Content_Type => 'form-data',
				Content => [
					file => [
					    undef,
					    'token.txt',
					    'Content_Type' => 'text/plain',
					    'Content' => $listRequest
					],
				],
			);
die $responseFormPost->status_line if !$responseFormPost->is_success;
my $contentPost = $responseFormPost->content;
my $decodedPost = decode_json($contentPost);
my $token = $decodedPost->{'summary'}->{'token'};

my $urlFound='http://www.reactome.org:80/AnalysisService/download/'.$token.'/entities/found/UNIPROT/uniprotFound.csv';
my $urlNotFound='http://www.reactome.org:80/AnalysisService/download/'.$token.'/entities/notfound/uniprotNotFound.csv';
my $urlPathway='http://www.reactome.org:80/AnalysisService/download/'.$token.'/pathways/TOTAL/pathway.csv';

##Get file with uniprot found in original list
my $downloadFound = $agent->get('http://www.reactome.org:80/AnalysisService/download/'.$token.'/entities/found/UNIPROT/uniprotFound.csv');
die $downloadFound->status_line if !$downloadFound->is_success;
my $foundSave = $tmpPathwayPathID."/uniprotFound.csv";
getstore($urlFound, $foundSave);

##Get file with uniprot not found in original list
my $downloadNotFound = $agent->get('http://www.reactome.org:80/AnalysisService/download/'.$token.'/entities/notfound/uniprotNotFound.csv');
die $downloadNotFound->status_line if !$downloadNotFound->is_success;
my $notFoundSave = $tmpPathwayPathID."/uniprotNotFound.csv";
getstore($urlNotFound, $notFoundSave);

##Get file with all information for pathways
my $downloadPathway = $agent->get('http://www.reactome.org:80/AnalysisService/download/'.$token.'/pathways/UNIPROT/pathway.csv');
die $downloadPathway->status_line if !$downloadPathway->is_success;
my $pathwaySave = $tmpPathwayPathID."/pathway.csv";
getstore($urlPathway, $pathwaySave);

if (-e $pathwaySave && -e $foundSave && -e $notFoundSave) {
	my $sthUpdatePathAna=$dbh->do("UPDATE PATHWAY_ANALYSIS SET STATUS=1 where ID_PATHWAY_ANALYSIS=$pathwayID");
	$dbh->commit;
	move($tmpPathwayPathID,"$projectPath/$pathwayID");
}
else {
	my $sthUpdatePathAna=$dbh->do("UPDATE PATHWAY_ANALYSIS SET STATUS=-1 where ID_PATHWAY_ANALYSIS=$pathwayID");
	$dbh->commit;
	remove_tree($tmpPathwayPathID);
}
$dbh->disconnect;

####>Revision history<####
# 1.0.0  new script, run reactome web service, replace runAndDisplayPathwayAnalysis (SL 19/11/14)