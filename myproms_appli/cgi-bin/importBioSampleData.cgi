#!/usr/local/bin/perl -w

################################################################################
# importBioSampleData.cgi             1.0.0                                    #
# Authors: P. Poullet, G. Arras, S.Liva (Institut Curie)              	       #
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
#use POSIX qw(strftime); # to get the time
use LWP::UserAgent;
use URI::Escape;
use promsConfig;
use promsMod;
use strict;

#print header; warningsToBrowser(1);#DEBUG

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};

#############################
####>Fetching Parameters<####
#############################
my $projectID = &promsMod::cleanNumericalParameters(param('projectID'));
my $context=param('context') || 'biosample'; # or property

#######################
####>Starting HTML<####
#######################
my ($lightColor,$darkColor)=&promsConfig::getRowColors;print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Biological Sample Data Import</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
|;

#########################
####>Submission form<####
#########################
unless (param('submitted')) {
	
	##>All Species
	my $title;
	my %allSpecies;
	my $speciesOptionStrg='';
	if ($context eq 'biosample') {
		$title='Importing Biological Samples and Annotations from file';
		my $dbh=&promsConfig::dbConnect;
		my $selSpecies = $dbh->prepare("SELECT ID_SPECIES, COMMON_NAME, SCIENTIFIC_NAME FROM SPECIES WHERE IS_REFERENCE = 1 ORDER BY COMMON_NAME");
		$selSpecies->execute;
		while (my($speciesID,$commonName,$scientifName) = $selSpecies -> fetchrow_array) {
			$speciesOptionStrg.=qq|<OPTION value="$speciesID">$commonName ($scientifName)</OPTION>|;
		}
		$selSpecies->finish;
		$dbh->disconnect;
	}
	else {
		$title='Importing Biological Sample Annotations from file';
	}
	
	my $bgColor=$lightColor;
	print qq
|<SCRIPT type="text/javascript">
function checkForm(myForm) {
	if (!myForm.dataFile.value) {
		alert('ERROR: No data file selected');
		return false;
	}
	return true;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" >
<CENTER>
<BR><FONT class="title">$title</FONT><BR><BR>
<FORM name="importForm" method="post" enctype="multipart/form-data" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="projectID" value="$projectID" />
<INPUT type="hidden" name="context" value="$context" />
<INPUT type="hidden" name="submitted" value="1" />
<TABLE border=0  bgcolor="$darkColor">
<TR><TH align="right" valign="top" width=100 nowrap>Select file :</TH><TD bgcolor="$lightColor"><INPUT type="file" name="dataFile"><BR>
<FONT class="font11"><B>Format:</B> tab-separated list of <B>Biological samples</B> and any number of associated <B>Properties</B> and/or <B>Treaments</B>.<BR>
<B>Header</B> required. Biological samples in <B>first</B> column. Prefix <B>Treatments</B> name with '<B>T:</B>' (eg. <B>T:DrugX</B>).
|;
	if ($context eq 'biosample') {
		print qq
|<BR><B>Optional</B> column: Sample description (header: <B>Description</B>)</FONT></TD></TR>
<TR><TH align="right" valign="top" nowrap>Species :</TH><TD bgcolor="$lightColor"><SELECT name="species" onchange="document.getElementById('speciesSPAN').style.visibility=(this.value=='in_file')? 'visible' : 'hidden'"><OPTION value="in-file">In file</OPTION>$speciesOptionStrg</SELECT>
	<SPAN id="speciesSPAN"><B>Format:</B> Column name must be <B>Species</B>. Use <B>taxonID</B> or <B>Scientific name</B>.</SPAN>&nbsp;</TD></TR>
|;
	}
	else {
		print "</FONT></TD></TR>\n";
	}
	print qq
|<TR><TH align="right" valign="top" nowrap>Option :</TH><TD bgcolor="$lightColor"><LABEL><INPUT type="checkbox" name="overwrite" value=1>Overwrite pre-existing conflicting annotation values.</LABEL></TD></TR>
<TR><TH colspan=2><INPUT type="submit" class="title3" value=" Import "></TH></TR>
</TABLE>
</FORM>
</CENTER>
</BODY>
</HTML>
|;
	exit;
}

##############################	
#####>Form was submitted<#####	
##############################
print qq
|<SCRIPT type="text/javascript">
function updateNavFrame() {
	top.promsFrame.selectedAction = 'properties';
	parent.itemFrame.location.reload(); // no impact on optionFrame
	parent.optionFrame.selectOption();
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" >
<CENTER>
<BR><FONT class="title">Importing Biological Sample data from file</FONT><BR><BR>
</CENTER>
|;

####>Parameters<####
my $dataFile = tmpFileName(upload('dataFile'));
my $refSpeciesID=($context eq 'property')? -1 : (param('species') eq 'in-file')? 0 : param('species');
my $overwrite=param('overwrite') || 0;
my $error='';

####>Parsing file<####
print "&bull; <FONT class=\"title3\">Reading data file...";
my (%sampleData,%samplePos,@annotationNames,%annotationType,%annotAllValues,%speciesList,%annot2BioSamples);
my $fileHead=`head -2 $dataFile`;
my $newLine='';
$newLine.='CR' if $fileHead=~/\r/;
$newLine.='LF' if $fileHead=~/\n/;
local $/=($newLine eq 'CR')? "\r" : ($newLine eq 'CRLF')? "\r\n" : "\n"; # Perl new line separator 
my ($speciesIdx,$descripIdx);
open(DATA,$dataFile) || die "Unable to open $dataFile";
while (my $line=<DATA>) {
	$line=~s/^ +//; # Clean leading space(s) if any
	$line=~s/\s+$//; # chomp is not enough <- Windows
	my ($bsName,@values)=split(/ *\t */,$line); # remove starting/trailing spaces
	if ($.==1) { # 1st line of the file
		@annotationNames=@values;
		foreach my $idx (0..$#values) {
			if (lc($values[$idx]) eq 'species') {
				$annotationType{$annotationNames[$idx]}='species';
				$speciesIdx=$idx;
			}
			elsif (lc($values[$idx]) eq 'description') {
				$annotationType{$annotationNames[$idx]}='description';
				$descripIdx=$idx;
			}
			else {
				my $annType='O'; # property (Other)
				if ($annotationNames[$idx]=~/^T:/) {
					$annotationNames[$idx]=~s/^T://;
					$annType='T';
				}
				$annotationType{$annotationNames[$idx]}=$annType;
			}
		}
		if (!$refSpeciesID && !defined $speciesIdx) {
			$error='No species column found';
			last;
		}
		next;
	}
	next unless $bsName;
	@{$sampleData{$bsName}}=(); # in case no annotations at all
	$samplePos{$bsName}=$.;
	foreach my $idx (0..$#values) {
		last if !defined $annotationNames[$idx]; # in case no header for value (should not happen unless badly-formatted file)
		if (!defined $values[$idx] || length($values[$idx])==0 || uc($values[$idx])=~/^N\/*A$/) {$values[$idx]=undef;}
		elsif ($annotationType{$annotationNames[$idx]}=~/^(O|T)$/) {
			$annotAllValues{$annotationNames[$idx]}{$values[$idx]}=1;
			$annot2BioSamples{$annotationNames[$idx]}{$bsName}=1;
		}
		push @{$sampleData{$bsName}},$values[$idx];
	}
	if (!$refSpeciesID && !$values[$speciesIdx]) {
		$error="No species declared for Biological sample '$bsName' at line $.";
		last;
	}
	$speciesList{$values[$speciesIdx]}=-1 if defined $speciesIdx; # list of distinct species
}
close DATA;
unlink $dataFile;
&exitOnError($error) if $error;

if (scalar keys %sampleData == 0) {
	&exitOnError('No Biological samples detected! Make sure your file is properly formatted (especially regarding compatible newline character)');
}
my $speciesInfo=($context eq 'property')? '' : ($refSpeciesID)? " (1 species)" : ' ('.(scalar keys %speciesList).' species)';
print " Done: ",scalar keys %sampleData," Biological samples$speciesInfo and ",scalar keys %annotAllValues," Properties found.</FONT><BR>\n";

####>Connecting to DB<####
my $dbh=&promsConfig::dbConnect;

####>Fetching list of existing Biological samples & matching properties<####
my (%projectBioSamples,%projBioSampProp);
my $sthPBS=$dbh->prepare("SELECT BS.NAME,BS.ID_BIOSAMPLE,BS.DES,BS.ID_SPECIES FROM BIOSAMPLE BS INNER JOIN PROJECT_BIOSAMPLE PB ON PB.ID_BIOSAMPLE=BS.ID_BIOSAMPLE WHERE PB.ID_PROJECT=$projectID");
my $sthBSProp=$dbh->prepare("SELECT ID_PROPERTY,PROPERTY_VALUE FROM BIOSAMPLE_PROPERTY WHERE ID_BIOSAMPLE=?");
$sthPBS->execute;
while (my ($bsName,@bsData)=$sthPBS->fetchrow_array) {
	next unless $sampleData{$bsName};
	@{$projectBioSamples{$bsName}}=@bsData;
	$sthBSProp->execute($bsData[0]); # id
	while (my ($propID,$propValue)=$sthBSProp->fetchrow_array) {
		$projBioSampProp{$bsName}{$propID}=$propValue;
	}
}
$sthPBS->finish;
$sthBSProp->finish;

####>Comparing lists of BioSamples<####
if ($context eq 'property') { # No new BioSamples allowed
	my @extraBioSamples;
	foreach my $bsName (keys %sampleData) {
		push @extraBioSamples,$bsName unless $projectBioSamples{$bsName};
	}
	my $numExtraBS=scalar @extraBioSamples;
	if ($numExtraBS) { #  => delete from list
		foreach my $bsName (@extraBioSamples) {
			delete $sampleData{$bsName};
		}
		foreach my $name (keys %annot2BioSamples) {
			foreach my $bsName (@extraBioSamples) {
				delete $annot2BioSamples{$name}{$bsName} if $annot2BioSamples{$name}{$bsName};
			}
			unless (scalar keys %{$annot2BioSamples{$name}}) { # delete unused annot
				delete $annotAllValues{$name};
				delete $annotationType{$name};
			}
		}
		my $bsStrg=($numExtraBS==1)? 'Biological sample was' : "$numExtraBS Biological samples were";
		print "&bull; <FONT class=\"title3\">The following $bsStrg not found in Project and will be ignored:<BR>";
		foreach my $bsName (sort{$samplePos{$a}<=>$samplePos{$b}} @extraBioSamples) {
			print "&nbsp;&nbsp;&nbsp;-$bsName<BR>\n";
		}
		print "</FONT><BR>\n";
	}
}

####>Mapping species to ID_SPECIES<####
if (!$refSpeciesID && scalar keys %speciesList) {
	my $baseURL='http://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon';
	my $sthSPSN=$dbh->prepare("SELECT ID_SPECIES FROM SPECIES WHERE SCIENTIFIC_NAME=? LIMIT 1");
	my $sthSPTX=$dbh->prepare("SELECT ID_SPECIES FROM SPECIES WHERE TAXONID=? LIMIT 1");
	my $sthInsSP=$dbh->prepare("INSERT INTO SPECIES (COMMON_NAME,SCIENTIFIC_NAME,TAXONID) VALUES (?,?,?)");
	my %uniqueSpecies;
	foreach my $spValue (keys %speciesList) {
		my $sthSP=($spValue=~/\D/)? $sthSPSN : $sthSPTX; # scientific name or taxonid
		$sthSP->execute($spValue);
		my ($speciesID)=$sthSP->fetchrow_array;
		if ($speciesID) {$speciesList{$spValue}=$speciesID;} # matched!
		else { # Fetch from EBI taxonomy DB
			my %speciesData;
			my $fullURL=$baseURL;
			my $isTaxonID=0;
			if ($spValue=~/\D/) { # scientific name
				$speciesData{'SCIENTIFIC_NAME'}=$spValue;
				$fullURL.='/scientific-name/'.uri_escape($spValue);		
			}
			else { # assume taxonid
				$speciesData{'TAXONID'}=$spValue;
				$fullURL.="/tax-id/$spValue";
				$isTaxonID=1;
			}
			my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr',ssl_opts => { verify_hostname => 0 });
			$agent->timeout(360);
			if (!$promsPath{'http_proxy'}) {$agent->env_proxy;}
			elsif (lc($promsPath{'http_proxy'}) eq 'no') {$agent->no_proxy;}
			else {$agent->proxy(['http','ftp'],$promsPath{'http_proxy'});}
			my $try=1;
			my $error='';
			while (!$error && $try) {
				my $response = $agent->get($fullURL);
				while (my $wait = $response->header('Retry-After')) {
					sleep $wait;
					$response = $agent->get($response->base);
				}
				if ($response->content=~/No results/) {
					my $infoStrg=($isTaxonID)? "Taxonid '$speciesData{TAXONID}'" : "Scientific name '$speciesData{SCIENTIFIC_NAME}'";
					$error="$infoStrg does not match any species";
					$try=0;
				}
				elsif ($response->is_success) {
					foreach my $line (split(/\n/,$response->content)) {
						if ($line=~/"taxId":\s+"(\d+)"/) {$speciesData{'TAXONID'}=$1 unless $speciesData{'TAXONID'};} # unless to keep user/1st in list
						elsif ($line=~/"scientificName":\s+"([^"]+)"/) {$speciesData{'SCIENTIFIC_NAME'}=$1 unless $speciesData{'SCIENTIFIC_NAME'};}
						elsif ($line=~/"commonName":\s+"([^"]+)"/) {$speciesData{'COMMON_NAME'}=$1;}
					}
					$try = 0;
					unless ($speciesData{'COMMON_NAME'}) {
						if (($isTaxonID && $speciesData{'SCIENTIFIC_NAME'}) || (!$isTaxonID && $speciesData{'TAXONID'})) { # common name does ot always exist
							$speciesData{'COMMON_NAME'}=$speciesData{'SCIENTIFIC_NAME'};
						}
						else {
							my $infoStrg=($isTaxonID)? "Taxonid '$speciesData{TAXONID}'" : "Scientific name '$speciesData{SCIENTIFIC_NAME}'";
							$error="$infoStrg returned unexpected results";
						}
					}
				}
				else {
					if ($try <= 5) {
						sleep 3; # required by server
						$try++;
					}
					else {
						$error="Got from EBI: ".$response->status_line." for ".$response->request->uri; #chomp()
					}
				}
			}
			if ($error) {
				$sthSPSN->finish;
				$sthSPTX->finish;
				$sthInsSP->finish;
				$dbh->rollback;
				$dbh->disconnect;
				&exitOnError($error);
			}
			else {
				if ($uniqueSpecies{$speciesData{'TAXONID'}}) { # same species!!! do not duplicate
					$speciesList{$spValue}=$uniqueSpecies{$speciesData{'TAXONID'}};
				}
				else { # 1rst time seen
					$sthInsSP->execute($speciesData{'COMMON_NAME'},$speciesData{'SCIENTIFIC_NAME'},$speciesData{'TAXONID'});
					$speciesList{$spValue}=$dbh->last_insert_id(undef,undef,'SPECIES','ID_SPECIES');
					$uniqueSpecies{$speciesData{'TAXONID'}}=$speciesList{$spValue};
				}
			}
		}
	}
	$sthSPSN->finish;
	$sthSPTX->finish;
	$sthInsSP->finish;
	print "&bull; <FONT class=\"title3\">",scalar keys %uniqueSpecies," new Species were added to myProMS.</FONT><BR>\n" if scalar keys %uniqueSpecies;
}

####>Checking stored properties values & updating if necessary<####
print "&bull; <FONT class=\"title3\">Processing Annotation data...";
my $sthProp=$dbh->prepare("SELECT P.ID_PROPERTY,P.NAME,P.PROPERTY_TYPE,P.POSSIBLE_VALUES FROM PROPERTY P,PROJECT_PROPERTY PP WHERE PP.ID_PROPERTY=P.ID_PROPERTY AND PP.ID_PROJECT=$projectID");
my $sthUpNewVal=$dbh->prepare("UPDATE PROPERTY SET POSSIBLE_VALUES=? WHERE ID_PROPERTY=?");
$sthProp->execute;
my %storedProperties;
my $numPropUpdated=0;
while (my ($propID,$name,$type,$possValueStrg)=$sthProp->fetchrow_array) {
	next unless defined $annotAllValues{$name};
	$storedProperties{$name}=$propID;
	my $numNewValues=0;
	next unless $possValueStrg; # nothing to update
	my %propValues;
	foreach my $value (split(':#:',$possValueStrg)) {$propValues{$value}=1;}
	my $numPrevValues=scalar keys %propValues;
	foreach my $value (keys %{$annotAllValues{$name}}) {$propValues{$value}=1;}
	my $numAllValues=scalar keys %propValues;
	next unless $numAllValues > $numPrevValues;
	if ($numAllValues <= 10) {
		$possValueStrg=join(':#:',sort{&promsMod::sortSmart($a,$b)} keys %propValues);
	}
	else {$possValueStrg=undef;} # too many possible values => do not record
	$sthUpNewVal->execute($possValueStrg,$propID);
	$numPropUpdated++;
}
$sthProp->finish;
$sthUpNewVal->finish;

####>Storing new properties<####
my $sthInsProp=$dbh->prepare("INSERT INTO PROPERTY (NAME,PROPERTY_TYPE,USE_IN_ANALYSIS,IS_VERIFIED,POSSIBLE_VALUES) VALUES (?,?,1,0,?)");
my $sthInsPP=$dbh->prepare("INSERT INTO PROJECT_PROPERTY (ID_PROJECT,ID_PROPERTY) VALUES (?,?)");
my $numPropAdded=0;
foreach my $name (@annotationNames) {
	next if (!$annotationType{$name} || $annotationType{$name} !~ /^(O|T)$/);
	next if $storedProperties{$name}; # already recorded
	#next unless $annotationType{$name}=~/^(O|T)$/;
	my @newValues;
	my $numNewValues=scalar keys %{$annotAllValues{$name}};
	my $possValueStrg=($numNewValues <= 10)? join(':#:',sort{&promsMod::sortSmart($a,$b)} keys %{$annotAllValues{$name}}) : undef;
	$sthInsProp->execute($name,$annotationType{$name},$possValueStrg);
	$storedProperties{$name}=$dbh->last_insert_id(undef,undef,'PROPERTY','ID_PROPERTY');
	$sthInsPP->execute($projectID,$storedProperties{$name});
	$numPropAdded++;
}
$sthInsProp->finish;
$sthInsPP->finish;
print " Done:";
print " $numPropUpdated annotations updated." if $numPropUpdated;
print " $numPropAdded annotations added." if $numPropAdded;
print " No update necessary." if (!$numPropUpdated && !$numPropAdded);
print "</FONT><BR>\n";

####>Storing/Updating Biological samples<####
print "&bull; <FONT class=\"title3\">Processing Biological sample data...";
my $sthInsBS = $dbh->prepare("INSERT INTO BIOSAMPLE (ID_SPECIES,NAME,DES,RECORD_DATE,UPDATE_DATE,UPDATE_USER) VALUES (?,?,?,NOW(),NOW(),?)");
my $sthUpBS = $dbh->prepare("UPDATE BIOSAMPLE SET ID_SPECIES=?,DES=?,UPDATE_DATE=NOW(),UPDATE_USER=? WHERE ID_BIOSAMPLE=?");
my $sthInsProjB = $dbh->prepare("INSERT INTO PROJECT_BIOSAMPLE(ID_PROJECT,ID_BIOSAMPLE) VALUES (?,?)");
my $sthBSPRK = $dbh->prepare("SELECT MAX(PROPERTY_RANK) FROM BIOSAMPLE_PROPERTY WHERE ID_BIOSAMPLE=?");
my $sthInsBSP = $dbh->prepare("INSERT INTO BIOSAMPLE_PROPERTY(ID_BIOSAMPLE,ID_PROPERTY,PROPERTY_RANK,PROPERTY_VALUE) VALUES (?,?,?,?)");
my $sthUpBSP = $dbh->prepare("UPDATE BIOSAMPLE_PROPERTY SET PROPERTY_VALUE=? WHERE ID_BIOSAMPLE=? AND ID_PROPERTY=?");
my $sthDelBSP = $dbh->prepare("DELETE FROM BIOSAMPLE_PROPERTY WHERE ID_BIOSAMPLE=? AND ID_PROPERTY=?");
my (%bsUpdated,%bsAdded);
foreach my $bsName (sort{$samplePos{$a}<=>$samplePos{$b}} keys %sampleData) {
	next if ($context eq 'property' && !$projectBioSamples{$bsName}); # no new BioSample added in this context
	my $des=(defined $descripIdx && defined $sampleData{$bsName}[$descripIdx])? $sampleData{$bsName}[$descripIdx] : undef;
	my $speciesID; # stays undef if context is property (not needed)
	if ($refSpeciesID > 0) {$speciesID=$refSpeciesID;}
	elsif (defined $speciesIdx) {
		my $spValue=$sampleData{$bsName}[$speciesIdx];
		$speciesID=$speciesList{$spValue};
	}
	my $bioSampID;
	if ($projectBioSamples{$bsName}) { # biosample already recorded
		$bioSampID=$projectBioSamples{$bsName}[0];
		if ($overwrite && $context eq 'biosample') { # BioSample is NOT overwritten if context is property
			my $okOverwrite=0;
			$okOverwrite=1 if ( ($des && (!$projectBioSamples{$bsName}[1] || $des ne $projectBioSamples{$bsName}[1])) || (!$des && $projectBioSamples{$bsName}[1]) );
			$okOverwrite=1 if $speciesID != $projectBioSamples{$bsName}[2];
			if ($okOverwrite) {
				$des=$projectBioSamples{$bsName}[1] unless $des;
				$speciesID=$projectBioSamples{$bsName}[2] unless $speciesID;
				$sthUpBS->execute($speciesID,$des,$userID,$bioSampID);
				$bsUpdated{$bioSampID}=1;
			}
		}
	}
	else {
		$sthInsBS->execute($speciesID,$bsName,$des,$userID);
		$bioSampID=$dbh->last_insert_id(undef,undef,'BIOSAMPLE','ID_BIOSAMPLE');
		$sthInsProjB->execute($projectID,$bioSampID);
		$bsAdded{$bioSampID}=1;
	}
	#>Sample Properties
	$sthBSPRK->execute($bioSampID);
	my ($rank)=$sthBSPRK->fetchrow_array;
	$rank=0 unless $rank;

	foreach my $idx (0..$#{$sampleData{$bsName}}) {
		next if ((defined $speciesIdx && $idx==$speciesIdx) || (defined $descripIdx && $idx==$descripIdx));
		my $propID=$storedProperties{ $annotationNames[$idx] };
		if ($projBioSampProp{$bsName} && $projBioSampProp{$bsName}{$propID}) {
			if ($overwrite && (!defined $sampleData{$bsName}[$idx] || $sampleData{$bsName}[$idx] ne $projBioSampProp{$bsName}{$propID})) {
				if (defined $sampleData{$bsName}[$idx]) {
					$sthUpBSP->execute($sampleData{$bsName}[$idx],$bioSampID,$propID);
				}
				else { # no value for property => delete BS_P entry
					$sthDelBSP->execute($bioSampID,$propID); # WARNING: introduces a gap in ranks
				}
				$bsUpdated{$bioSampID}=1;
			}
		}
		elsif (defined $sampleData{$bsName}[$idx]) {
			$rank++;
			$sthInsBSP->execute($bioSampID,$propID,$rank,$sampleData{$bsName}[$idx]);
			$bsUpdated{$bioSampID}=1 unless $bsAdded{$bioSampID};
		}
	}
}
$sthInsBS->finish;
$sthUpBS->finish;
$sthInsProjB->finish;
$sthBSPRK->finish;
$sthInsBSP->finish;
$sthUpBSP->finish;
$sthDelBSP->finish;
print " Done:";
my $numBsUpdated=scalar keys %bsUpdated;
my $numBsAdded=scalar keys %bsAdded;
print " $numBsUpdated Biological Samples updated." if $numBsUpdated;
print " $numBsAdded Biological Samples added." if $numBsAdded;
print " No update necessary." if (!$numBsUpdated && !$numBsAdded);
print "</FONT><BR>\n";

$dbh->commit;
$dbh->disconnect;

print qq
|<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="button" class="title2" value=" Done " onclick="updateNavFrame()"/>
</BODY>
</HTML>
|;

sub exitOnError {
	my ($errorStrg)=@_;
	return unless $errorStrg;
	print qq
|<BR><BR><FONT class="title2" color="#DD0000">ERROR: $errorStrg!</FONT>
<BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="button" class="title2" value=" OK " onclick="parent.optionFrame.selectOption();"/>
<BR><BR>
</BODY>
</HTML>
|;
	exit;
}

####>Revision history<####
# 1.0.0 [FEATURE] Batch import of biological sample and annotations from file (PP 16/03/20)
