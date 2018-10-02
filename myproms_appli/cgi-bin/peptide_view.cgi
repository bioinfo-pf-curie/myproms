#!/usr/local/bin/perl -w

################################################################################
# peptide_view.cgi              3.2.1                                          #
# Authors: P. Poullet, G. Arras, F. Yvon, M. Le Picard (Institut Curie)        #
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
use strict ;
use CGI ':standard';
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use XML::SAX;
use XML::Twig;
use XML::Simple; # used by &promsMod::extractSpectrumMSF
#use GD;

#print header; warningsToBrowser(1); #debug
#######################################
####>Configuration & Connection DB<####
#######################################
my %promsPath=&promsConfig::getServerInfo;
my %msTypeNames=&promsConfig::getMsType;
my %massValueAA=&promsConfig::getMassAAmono;
my %massValueAtom =&promsConfig::getMassATMono;
my %fragmentDef=&promsConfig::getFragmentDef; #value of fragmentation
my %fragmentClassif=&promsConfig::getFragmentClassif;

my $userID=$ENV{'REMOTE_USER'};
my $dbh=&promsConfig::dbConnect;
my ($color1,$color2)=&promsConfig::getRowColors;
my $maxRef=&promsConfig::getMaxRefSpectra - 1;
my ($higlightColor,$whiteColor,$blackColor,$colorA,$colorB)=&promsConfig::getSpectrumColors;

##############################
####>Recovering parameter<####
##############################
my $DatFile = param('file');
unless (-e $DatFile) {
	print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Missing Spectrum File</TITLE>
<LINK rel=stylesheet href="$promsPath{html}/promsStyle.css" type=text/css>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<FONT class="title2" color="#DD0000"><BR><BR>ERROR: Spectrum data file was not found!</FONT>
<BR><BR><FONT class="title3">Check your Analysis import procedure or contact your Mass Spectrometry facility.</FONT>
|;
	exit;
}
my $sequence=param('SEQUENCE') || '';
my $massObs=param('massObs') || 0;
my $massExp=param('massExp') || 0;
my $query_number=param('query') || 0;
my $rank=param('hit') || 0;
my $objectID=&promsMod::cleanNumericalParameters(param('ID')) || 0;
my $projectID=&promsMod::cleanNumericalParameters(param('PROJECT_ID')); # needed for spectrum extraction from msf file
my $refSpectrum=param('REF') || 0;
my ($refPos,$numRef)=(param('REF_POS'))? split(/_/,param('REF_POS')): (0,0); # eg: 02 -> 1st in 3
my $call=param('CALL') || '0';
my $disableStore=param('disStore') || '';
my $refScore=param('REF_SC');
my $upUser=param('UP_USER') || '';
my $fileType = param('TYPE');
if ($fileType eq "REFSPEC") {
	$fileType="MASCOT.DAT" if $DatFile =~ /\.dat/;
	$fileType="PHENYX.XML" if $DatFile =~ /\.pgf/;
	$fileType="sptxt" if $DatFile =~ /\.sptxt/;
	if ($DatFile =~ /\.pdm/) {
		$fileType="MASCOT.PDM";
		my $onParameters=0;
		open (FILE, $DatFile) || die "Unable to open $DatFile";
		while (my $line=<FILE>) {
			if ($line =~ /name=\"parameters\"/) {
				$onParameters=1;
			}
			if ($onParameters && $line =~ /SEARCH_ALGO=SEQUEST/){
				$fileType="SEQUEST.PDM";
			}
		}
		close FILE;
	}
	if ($fileType eq "PHENYX.XML") { # paragon and phenyx have the same extension -> pgf / (paragon|phenyx) generic file...
		my ($isParagon,$onMasses)=(0,0);
		open (FILE, $DatFile) || die "Unable to open $DatFile";
		while (my $line=<FILE>) {
			chomp($line);
			if ($line =~ /Paragon/) {
				$fileType="PARAGON.XML";
				last;
			}
		}
		close FILE;
	}
}
my $minValue=param('minValue') if param('minValue');
my $maxValue=param('maxValue') if param('maxValue');
my $mzIntervall = 1 if (defined $minValue && defined $maxValue && $minValue<$maxValue);
my $defMinValue=99999; my $defMaxValue=0;
my $frameWidth = (param('width') && param('width')=~/(\d+)/) ? $1 : "1050";#IE dont send a corrrect value, return undefined
my $spectrumWidth=(param('spWidth') && param('spWidth')=~/(\d+)/) ? $1 : $frameWidth-50;
my $spectrumHeight=(param('spHeight') && param('spHeight')=~/(\d+)/) ? $1 : "220";
#my ($spectrumID)=($fileType eq "MASCOT.DAT")? ($DatFile=~/.*S(\d+)\.dat/) : ($fileType eq "PHENYX.XML" || $fileType eq "MASCOT.XML")? ($DatFile=~/.*S(\d+)\.pgf/) : "" if $refSpectrum > 0;
my ($spectrumID)=($fileType eq "MASCOT.DAT")? ($DatFile=~/.*S(\d+)\.dat/) : ($fileType eq "PHENYX.XML" || $fileType eq "MASCOT.XML" || $fileType eq "PARAGON.XML")? ($DatFile=~/.*S(\d+)\.pgf/) : ($fileType =~ /\.PDM/)? ($DatFile=~/.*S(\d+)\.pdm/)  : ($fileType =~ /sptxt/) ? ($DatFile=~/SwLib_(\d+)\/.*\.sptxt/) : "" if $refSpectrum > 0;

my $useInternalFragment = 0;
my $maxFragmentMass = 900;
my $warningString='';
my $searchRank='';# For Paragon
my $checkRefSpectrum=1;
my %massValueFixedMods;   #   %massValueFixedMods{residue}  = delta mass
my %massValueMods;
my %massValueFragMod = ('PHOS'=>97.976896, 'OxiM'=>64);  # trouver valeur exactes
my (%allowModif,%allowVARModif);

########################################
####>Recovering information from DB<####
########################################
my ($analysisID,$scan,$charge,$peptideSequence,$validStatus,$titleSequence,$varModsString,$sub,$comments,$instrument,$elutionTime,$sthVmodMass,$spectrumCode,$wiffFile); # global variable set at this section
my $extSpectrumID;
my $libID;
if ($call eq 'lib') {
	$elutionTime=param('irt');
	($peptideSequence,$charge,$varModsString)=split(/_/,$sequence);
	($sequence=$peptideSequence)=~s/n*\[\w*\]//g;
	$titleSequence=$sequence;
	$libID=&promsMod::cleanNumericalParameters(param('libID'));
	($instrument)=$dbh->selectrow_array ("SELECT INSTRUMENT FROM SWATH_LIB WHERE ID_SWATH_LIB=$libID");
	$instrument='ESI-QUAD-TOF' unless $instrument;
	#$varModsString=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$objectID,$analysisID,$sequence);
	$varModsString=~s/@/\+/g;
	$varModsString=($varModsString)? $varModsString : '';
	$sthVmodMass=$dbh->prepare("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES,SPECIFICITY,MONO_MASS FROM MODIFICATION");
	$sthVmodMass->execute;
	while (my ($psiMsName,$interimName,$altNames,$specificity,$monomass) = $sthVmodMass->fetchrow_array ) {
		$altNames=~ s/##/,/g if $altNames;
		$altNames=~ s/^,//; $altNames=~ s/,\Z//; # remove starting & trailing "," if any
		my $name=($psiMsName)? $psiMsName : ($interimName)? $interimName : ($altNames)? $altNames : '';
		$massValueMods{$name}=$monomass;
	}
	$sthVmodMass->finish;
}
else {
	if ($refSpectrum <= 0) {   ####>data spectrum<####
		if ($call eq 'rank' || $call eq 'ana' || $call eq 'seq' ) {
			($analysisID,$charge,$elutionTime,my $massData,my $infoPep)=$dbh->selectrow_array ("SELECT ID_ANALYSIS,CHARGE,ELUTION_TIME,MASS_DATA,INFO_PEP$rank FROM QUERY_VALIDATION WHERE ID_QUERY=$objectID");
			#($validStatus,$instrument,$wiffFile)=$dbh->selectrow_array ("SELECT VALID_STATUS,INSTRUMENT,WIFF_FILE FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
			#my ($instrument,$massData,$infoPep)=$dbh->selectrow_array ("SELECT INSTRUMENT,MASS_DATA,INFO_PEP$rank FROM ANALYSIS A,QUERY_VALIDATION QV WHERE A.ID_ANALYSIS=QV.ID_ANALYSIS AND ID_QUERY=$objectID");
			($sequence) = ($infoPep=~ /SEQ=(\w+)/);
			($massObs) = ($massData=~/OBS=(\d+\.*\d*)/);
			($massExp) = ($massData=~/EXP=(\d+\.*\d*)/);
			($searchRank) = ($infoPep=~/SRK=(\d+)/)? $1 : '';
			$varModsString=&promsMod::toStringVariableModifications($dbh,'QUERY',$objectID,$analysisID,$sequence,$rank);
			if ($infoPep=~ /COM=([^,]*),?/) {
				$comments = $1;
				chomp ($comments);
			}
			$titleSequence=$sequence;
			($sub)=($infoPep=~/SUBST=([^,]+)/);
			$sthVmodMass=$dbh->prepare("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES,M.SPECIFICITY,MONO_MASS FROM ANALYSIS_MODIFICATION AM,MODIFICATION M WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND ID_ANALYSIS=$analysisID");
			$spectrumCode="Q$objectID"."R$rank";
		}
		elsif ($call eq 'pep' || $call eq 'val') {
			($analysisID,$charge,$elutionTime,$sequence,$massObs,$massExp,$sub,$comments)=$dbh->selectrow_array("SELECT ID_ANALYSIS,CHARGE,ELUTION_TIME,PEP_SEQ,MR_OBS,MR_EXP,SUBST,COMMENTS FROM PEPTIDE WHERE ID_PEPTIDE=$objectID");
			$varModsString=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$objectID,$analysisID,$sequence);
			$titleSequence=$sequence;
			$sthVmodMass=$dbh->prepare("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES,M.SPECIFICITY,MONO_MASS FROM ANALYSIS_MODIFICATION AM,MODIFICATION M WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND ID_ANALYSIS=$analysisID");
			$spectrumCode="P$objectID";
		}
		else { #error !	for debug
			die "<H3>Unknow call parameter $call</H3><BR>";
		}
		($validStatus,$instrument,$wiffFile)=$dbh->selectrow_array("SELECT VALID_STATUS,INSTRUMENT,WIFF_FILE FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
	}
	else {      ####>reference spectrum<####
		if ($fileType eq "sptxt"){
			($analysisID,$charge,$sequence,$sub,$comments)=$dbh->selectrow_array("SELECT ID_ANALYSIS,CHARGE,PEP_SEQ,SUBST,COMMENTS FROM PEPTIDE WHERE ID_PEPTIDE=$query_number");
			$varModsString=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$query_number,$analysisID,$sequence);
			$titleSequence=$sequence;
			$sthVmodMass=$dbh->prepare("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES,M.SPECIFICITY,MONO_MASS FROM SWATH_LIB_MODIFICATION SLM,MODIFICATION M WHERE SLM.ID_MODIFICATION=M.ID_MODIFICATION AND ID_SWATH_LIB=$spectrumID");
			$spectrumCode="S$spectrumID";

		}else{
			($elutionTime,$massObs,$massExp,$sequence,$sub,$comments,$instrument)=$dbh->selectrow_array("SELECT ELUTION_TIME,MR_OBS,MR_EXP,PEP_SEQ,SUBST,COMMENTS,INSTRUMENT FROM SPECTRUM WHERE ID_SPECTRUM=$spectrumID");
			$charge=sprintf "%1D",($massExp/$massObs)+.5; # $massObs & $massExp must be defined!
			$varModsString=&promsMod::toStringVariableModifications($dbh,'SPECTRUM',$spectrumID,$analysisID,$sequence);
			$titleSequence=$sequence;
			$sthVmodMass=$dbh->prepare("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES,M.SPECIFICITY,MONO_MASS FROM SPECTRUM_MODIFICATION SM,MODIFICATION M WHERE SM.ID_MODIFICATION=M.ID_MODIFICATION AND ID_SPECTRUM=$spectrumID");
			$spectrumCode="S$spectrumID";
		}
	}
	$instrument='ESI-QUAD-TOF' unless $instrument;
	$instrument =~ s/\s+\Z//;

	$sthVmodMass->execute;
	while (my ($psiMsName,$interimName,$altNames,$specificity,$monomass) = $sthVmodMass->fetchrow_array ) {
		$altNames=~ s/##/,/g if $altNames;
		$altNames=~ s/^,//; $altNames=~ s/,\Z//; # remove starting & trailing "," if any
		my $name=($psiMsName)? $psiMsName : ($interimName)? $interimName : ($altNames)? $altNames : '';
		$massValueMods{$name}=$monomass;
		#foreach my $specString (split(/,/,$specificity)) {
		#	my $aa=$specString;
		#	$aa='N-term' if $specString =~ /=/;
		#	$aa='Protein N-term' if $specString =~ /-/;
		#	$aa='C-term' if $specString =~ /\*/;
		#	$aa='Protein C-term' if $specString =~ /\+/;
		#	#$massValueMods{"$name ($aa)"}=$monomass;
		#
		#}
	}
	$sthVmodMass->finish;

	if ($sub) { # SUBSTITUTION INFORMATION
		my %susbtitutions; # Could exist several substitutions on the same peptide !
		foreach my $subst (split(' \+',$sub)) {
			next unless $subst;
			if ($subst =~ /(\w)->(\w) \((\d+)\)/){
				my ($origin,$substituted,$pos)=($1,$2,$3);
				$susbtitutions{$pos}="($1&rarr;$2)";
				substr($sequence , $pos-1 , 1 , $substituted );
			}
		}
		foreach my $pos (sort {$b <=> $a} keys %susbtitutions) { # Modify the sequence in descending order so as to keep the AA position integrity (because string replacement size is bigger thant one letter... )
			substr($titleSequence , $pos-1 , 1 , $susbtitutions{$pos} );
		}
	}

	####>number of the spectrum if it is a query from a MSF file
	#my $extSpectrumID;
	if ($fileType=~/\.PDM/ && param('TYPE') ne 'REFSPEC') {
		if ($call eq 'rank' || $call eq 'ana' || $call eq 'seq'){
			($extSpectrumID)=$dbh->selectrow_array("SELECT EXT_SPECTRUMID FROM QUERY_VALIDATION WHERE ID_QUERY=$objectID");
		}
		elsif ($call eq 'pep' || $call eq 'val') {
			#($extSpectrumID)=$dbh->selectrow_array("SELECT EXT_SPECTRUMID FROM QUERY_VALIDATION,PEPTIDE WHERE QUERY_VALIDATION.ID_ANALYSIS=PEPTIDE.ID_ANALYSIS AND ID_PEPTIDE=$objectID AND PEPTIDE.QUERY_NUM=QUERY_VALIDATION.QUERY_NUM");
			my ($data)=$dbh->selectrow_array("SELECT DATA FROM PEPTIDE WHERE ID_PEPTIDE=$objectID");
			($extSpectrumID)=($data=~/EXT_SPECTRUMID=(\d+)/) if $data;
		}
	}
	elsif ($fileType eq 'MAXQUANT.DIR') { # MAXQUANT
		$elutionTime=~s/;sc(\d+);*//;
		$extSpectrumID=$1; # scan number
	}
	elsif($fileType=~/^TDM\./){
		($scan)=($elutionTime=~/sc(\d+);et.+/);
		($elutionTime)=($elutionTime=~/sc\d+;et(.+)/);
	}

	if ($elutionTime) { # For PARAGON with RT problems -> et-1.00;sp1.1.1.17512.1.1;
		if ($elutionTime=~/et(-?\d+\.*\d*);/) {
			$extSpectrumID=$elutionTime unless $extSpectrumID;
			$elutionTime=$1;
			$elutionTime='?' if $elutionTime eq '-1.00';
		}
		else {
			if($elutionTime =~/sp\d+;sc\d+;/){$elutionTime = '?';}
			#else{$elutionTime="$elutionTime min.";}
		}
		$elutionTime=~s/^(et|rt)//;
	}
	else {$elutionTime='?';}

}
$massExp=($massObs-1.007825032)*$charge unless $massExp || $fileType eq "sptxt";
my $isRef = ($refSpectrum <= 0)? "false" : "true";


####>Recovering info for Instrument settings
my ($usrParam,$refFragRules,$deltaParent,$deltaFrag,$tolFrag,$nrLevel) = &promsMod::getInstrParam ($dbh,$instrument,$userID); #$usrParam = -1 : instrument not defined; 0 : defined for all user, 1 user defined
my %fragmentationRules=%{$refFragRules};

####>user Preference
my $userPref =$dbh->selectrow_array("SELECT USER_PREF FROM USER_LIST WHERE ID_USER='$userID'");
$userPref='' unless $userPref;
my $useSpecApp=($userPref=~/specApp=1/)? 1 : 0; # using interactive SVG graph
my $isVert = ($userPref=~/vertLab=1/)? "true" : "false";

####>deconnecting DB
$dbh->disconnect;

#settings  $intPercent and $matchInterval, by defualt values from instrument or user préférence for THIS spectrum
my $intPercent=(defined(param('miniInt')) && param('miniInt')>0 && param('miniInt')<=1)? param('miniInt') : $nrLevel;
my $matchInterval=(defined(param('tolerance')) && param('tolerance')>0)? param('tolerance') : $tolFrag;

##############################################################################
####>Recovering experimental ionic serie and other information from files<####
##############################################################################

my %expIonTable; #  $expIonTable{mz}[]  ; [0] intensite, [1] type de fragmentaion si match
my $query_desc;  #description du query titre, elution etc ...
my $fixedModsString;
my %FixMods;  # modification Fixe   [0] description    [1]  deltamass    [2] residue where modif occurs
if ($fileType=~/\.(DAT|PDM)/) { #use .dat type file for params

	my $onQuery=0;
	my $onParameters=0;
	my $onMasses=0;
	my $localIonSerie;
	my $parentFile='';
	open (FILE, $DatFile) || die "Unable to open $DatFile";
	#my $useAsDat=($fileType=~/\.PDM/ && $call ne 'pep' && param('TYPE') ne 'REFSPEC')? 0 : 1;
	my $useAsDat=($fileType=~/\.PDM/ && param('TYPE') ne 'REFSPEC')? 0 : 1; # pdm spectra are extracted in Sxxxxxx.pdm file
	unless ($useAsDat) {
		my @splitDatFile = split(/\//,$DatFile);
		my $msfFile = $splitDatFile[-1];
		$msfFile=~ s/\_[0-9]*\.*[0-9]*\.pdm/\.msf/g;
		#my $dbsqlite = DBI->connect("dbi:SQLite:$promsPath{valid}/multi_ana/proj_$projectID/$msfFile", "", "", {PrintError => 1,RaiseError => 1});
		#my ($zipSpectrumData)=$dbsqlite->selectrow_array("SELECT Spectrum FROM SpectrumHeaders, MassPeaks, Spectra WHERE Spectra.UniqueSpectrumID=SpectrumHeaders.UniqueSpectrumID AND SpectrumHeaders.MassPeakID=MassPeaks.MassPeakID AND SpectrumHeaders.SpectrumID=$extSpectrumID;");
		#my $zipFile="$promsPath{tmp}/spectrum_$userID.zip"; #"/tmp/myproms/batch/garras/spectrum.zip";
		#open (OUTPUT, ">$zipFile") || die "Error: Cannot write $zipFile file!\n";
		#binmode OUTPUT;
		#print OUTPUT $zipSpectrumData;
		#close(OUTPUT);
		#my $xmlString=`unzip -p $zipFile`;#copy in the stdout the xml file that is used directly in perl
		#system "rm $zipFile";
		#my $xml = new XML::Simple();
		#my $xmlData = $xml->XMLin($xmlString);
		#foreach my $value (@{$xmlData->{'PeakCentroids'}{'Peak'}}) { #copy all the ions for this specific query
		#	my ($mz,$intens)=($value->{'X'},$value->{'Y'});
		#	$expIonTable{$mz}[0]=sprintf "%.4f",$intens;
		#}
		#my $retTimeSec=$xmlData->{'Header'}{'SpectrumIdentifiers'}{'SpectrumIdentifier'}{'RetentionTime'}*60;
		#$query_desc="title=Spectrum$xmlData->{'Header'}{'SpectrumID'}%20scans%3a$xmlData->{'Header'}{'SpectrumIdentifiers'}{'SpectrumIdentifier'}{'ScanNumber'}%2c";
		#$query_desc =~ s/%([\dA-Fa-f][\dA-Fa-f])/pack ("C", hex($1))/eg;
		#$query_desc.=" Elution: $retTimeSec sec";
		#$query_desc.=sprintf "(%.2f min).",$xmlData->{'Header'}{'SpectrumIdentifiers'}{'SpectrumIdentifier'}{'RetentionTime'};
		my $msfFullFile=($validStatus==2)? "$promsPath{peptide}/proj_$projectID/ana_$analysisID/$msfFile" : "$promsPath{valid}/multi_ana/proj_$projectID/$msfFile";
		($parentFile,my ($refIons,$retTimeMin,$spectrumIDXML,$scanNumber))=&promsMod::extractSpectrumMSF($msfFullFile,$extSpectrumID);
		foreach my $value (@{$refIons}) { # copy all the ions for this specific query
			my ($mz,$intens)=($value->{'X'},$value->{'Y'});
			$expIonTable{$mz}[0]=sprintf "%.4f",$intens;
		}
		my $retTimeSec=$retTimeMin*60;
		$query_desc="title=Spectrum$spectrumIDXML%20scans%3a$scanNumber%2c";
		$query_desc =~ s/%([\dA-Fa-f][\dA-Fa-f])/pack ("C", hex($1))/eg;
		$query_desc.=" Elution: $retTimeSec sec";
		$query_desc.=sprintf "(%.2f min).",$retTimeMin;
	}
	while (my $line=<FILE>) {
		if ($line =~ /name=\"parameters\"/) {
			$onParameters=1;
		}
		if ($onParameters==1 && $line =~ /FILE=(.+)/ && !$parentFile) {
			$parentFile=(split(/[\/\\]/,$1))[-1];
			chomp $parentFile;
		}
		elsif ($line =~ /name=\"masses\"/) {
			$onMasses=1;
		}
		elsif ($line =~ /name=\"query$query_number\"/) {
			$onQuery=1;
		}
		elsif ($line =~ /gc0p4Jq0M2Yt08jU534c0p/) {
			$onQuery=0;
			$onParameters=0;
			$onMasses=0;
		}
		elsif ($onMasses==1 && $line =~ /FixedMod(\d+)=(.+)\,(.+)/) {
			$FixMods{$1}[0]= $3; #$3 description
			$FixMods{$1}[1]= $2; #$2 delta mass
		}
		elsif ($onMasses==1 && $line =~ /FixedModResidues(\d+)=(\w+)/) {
			$FixMods{$1}[2]= $2; #residue who fix mods occurs
		}
		elsif ($onMasses==1 && $line =~ /delta\d+=(.+)\,(.+\))/) {
			my ($mass,$modifName)=($1,$2);
			$modifName=~s/ \([^(]+\)$//; # removes residu info
			$massValueMods{$modifName} = $mass unless $massValueMods{$modifName}; # Do not overwrite mass from MODIFICATION table <- Unimod
		}
		elsif ($onQuery==1 && $line =~ /^title=(.+)/ && $useAsDat) {
			$query_desc=$1;
			$query_desc =~ s/%([\dA-Fa-f][\dA-Fa-f])/pack ("C", hex($1))/eg;
		}
		elsif ($onQuery==1 && $line =~ /^rtinseconds=(\d+)/ && $useAsDat) { # Orbitrap
			my $retTimeSec=$1;
			my $retTimeMin=sprintf "%.2f",$retTimeSec/60;
			$query_desc.=" Elution: $retTimeMin min. ($retTimeSec sec.)";
		}
		elsif ($onQuery==1 && $line =~ /Ions1=(.+)/ && $useAsDat) {
			$localIonSerie = $1 ;
			chomp ($localIonSerie) ;
			foreach my $Ion (split(/,/,$localIonSerie)) {
				my ($mz,$intens)=split(/:/,$Ion);
				$expIonTable{$mz}[0]=sprintf "%.4f",$intens;
			}
			last;
		}
	}
	close FILE;
	$query_desc="File: $parentFile, $query_desc" unless $query_desc=~/File/;

	####>creating% massValueFixedMods hash<####
	foreach my $numMod (keys (%FixMods)) {
		$massValueFixedMods{$FixMods{$numMod}[2]} = $FixMods{$numMod}[1] ;
	}
}


elsif ($fileType eq "PHENYX.XML" || $fileType eq "MASCOT.XML") {#use .pgf file
	my $onQuery=0;
	my $onMasses=0;
	my $onAAMasses=0;
	my $onParameters=0;
	my $localIonSerie;
	my $parentFile='';

	my $pgfFile = $DatFile;
	$pgfFile =~ s/\.xml/\.pgf/;
	open (FILE, $pgfFile) || die "Unable to open $pgfFile";
	while (my $line=<FILE>) {
		if ($line =~ /name=\"masses\"/) {
			$onMasses=1;
		}
		elsif ($line =~ /name=\"AAmass\"/) {
			$onAAMasses=1;
		}
		elsif ($line =~ /name=\"parameters\"/) {
			$onParameters=1;
		}
		if ($onParameters==1 && $line =~ /FILE=(.+)/) {
			$parentFile=$1;
			chomp $parentFile;
		}
		elsif ($line =~ /name=\"query$query_number\"/) {
			$onQuery=1;
		}
		elsif ($line =~ /gc0p4Jq0M2Yt08jU534c0p/) {
			$onQuery=0;
			$onMasses=0;
			$onAAMasses=0;
			$onParameters=0;
		}
		elsif ($fileType eq "MASCOT.XML" && $onParameters==1 && $line =~/^MODS=(.*)/) {
			$fixedModsString = $1;
		}
		elsif ($onMasses==1 && $line =~/delta\d+=(.+)\,(.+)/) {
			my ($mass,$modifName)=($1,$2);
			$modifName=~s/ \([^(]+\)$//; # removes residu info
			$massValueMods{$modifName} = $mass unless $massValueMods{$modifName}; # Do not overwrite mass from MODIFICATION table <- Unimod
		}
		elsif ($onMasses==1 && $line =~/Fixed=(.+)\,(.+)\,(.+)/) {		#Fixed=$massModif,$modifName,$residue
			$massValueFixedMods{$3}=$1 ;
 			# %FixMods ;  # modification Fixe   [0] description    [1]  deltamass    [2] residue where modif occurs
			$FixMods{$2}[0] =$2;
			$FixMods{$2}[1] =$1;
			$FixMods{$2}[2] =$3;
		}
		elsif ($fileType eq "MASCOT.XML" && $onAAMasses==1 && $line=~/(\w+)=(\d+.?\d*)/) {
			$massValueAA{$1} = $2;
		}
		elsif ($onQuery==1 && $line=~/title=(.+)/) {
		$query_desc = $1;
		}
		elsif ($onQuery==1 && $line =~ /Ions1=(.+)/) {
			$localIonSerie = $1;
			chomp ($localIonSerie);
			my @expIonArray  = split(/,/,$localIonSerie);
			foreach my $Ion (@expIonArray ) {
				my ($mz,$intens) = split(/:/,$Ion);
				$expIonTable{$mz}[0] = sprintf("%.4f",$intens);
			}
		}
	}
	close FILE;
	$query_desc="File: $parentFile, $query_desc" unless $query_desc=~/File/;
}


elsif ($fileType eq "PARAGON.XML") {
	###> 1st: read the PGF file with variable modification and search descriptions
	my $vmodFile=$DatFile;
	$vmodFile=~s/\.xml/\.pgf/g;
	open (DESC, $vmodFile) || die "Unable to open $vmodFile";
	my ($onMasses,$onSearch,$onQuery,$localIonSerie)=(0,0,0,'');
	my %search;
	while (my $line=<DESC>) {
		chomp($line);
		if ($line eq '--gc0p4Jq0M2Yt08jU534c0p' ){
			($onMasses,$onSearch,$onQuery)=(0,0,0);
			next;
		}
		elsif ($line eq 'Content-Type: Paragon; name="masses"'){
			$onMasses=1;
			next;
		}
		elsif ($line eq 'Content-Type: Paragon; name="search"'){
			$onSearch=1;
			next;
		}
		elsif ($line =~ /name=\"query$query_number\"/) {
			$onQuery=1;
		}
		if ($onQuery && $line =~ /title=(.+)/) {
			$query_desc="File: $1";
			chomp $query_desc;
		}
		elsif ($onQuery==1 && $line =~ /Ions1=(.+)/) {
			my $localIonSerie = $1;
			chomp ($localIonSerie);
			my @expIonArray  = split(/,/,$localIonSerie);
			foreach my $Ion (@expIonArray ) {
				my ($mz,$intens) = split(/:/,$Ion);
				$expIonTable{$mz}[0] = sprintf("%.4f",$intens);
			}
		}
		if ($onMasses && $line =~/delta\d+=(.+),(.+)/){
			my ($mass,$modifName)=($1,$2);
			$modifName=~s/ \([^)]+\)$//; # removes residu info
			$massValueMods{$modifName} = $mass unless $massValueMods{$modifName}; # Do not overwrite mass from MODIFICATION table <- Unimod
		}
		if ($onSearch && $line !~ 'Content' ){
			my ($searchID,$rawfileName)=split(/,/,$line);
			$search{$searchID}=$rawfileName;
		}
	}
	close DESC;

	###> 2nd: read the spectrum
	open (FILE, $DatFile) || die "Unable to open $DatFile";
	if($extSpectrumID=~/sp([\d+|\.]*);/){# Get the spectrum xml:id
		$extSpectrumID=$1;
	}
	if(!$localIonSerie) {
		my $onMSMSpeaks=0;
		my $nbMatch=0;
		while (my $line=<FILE>) {
			chomp($line);
			$onQuery=1 if ($line =~ /<SPECTRUM .* xml:id="$extSpectrumID/ ) ;
			next unless $onQuery;
			if ($onMSMSpeaks) {
				my $peakInfo=$line;
				$peakInfo =~ s/<\!\[CDATA\[//g;
				$peakInfo =~ s/\]\]>//g;
				my ($mz,$chargeState,$intens)=split(/\t/,$peakInfo);
				$expIonTable{$mz}[0]=sprintf "%.4f",$intens;
			}
			$onMSMSpeaks=1 if ($line =~ /<MSMSPEAKS/);
			if ($line =~ /<MATCH/){
				$nbMatch++;
				if ($nbMatch == $searchRank) {
					my ($searchID)=($line =~ /<MATCH .* searches=\"(.+)\" seq=/);
					$query_desc="File: $search{$searchID} (spectrum number in ProteinPilot software: $extSpectrumID) ";
				}
			}
			last if ($line=~ /]]>/);
		}
		close FILE;
	}
}

elsif ($fileType eq 'MAXQUANT.DIR') {
	my %msmsColNum;
	open (FILE, $DatFile) || die "Unable to open $DatFile";
	while (<FILE>) {
		my @parameters=split(/ *\t */,$_); # remove starting/trailing spaces
        if ($.==1) {
			my $ncol=0;
			foreach my $colName (@parameters) {
				$msmsColNum{$colName}=$ncol;
				$ncol++;
			}
			next;
        }
		if ($parameters[$msmsColNum{'Scan number'}] && $parameters[$msmsColNum{'Scan number'}]==$extSpectrumID) {
			#my $charge=$parameters[$msmsColNum{'Charge'}];
			my @intensities=split(';',$parameters[$msmsColNum{'Intensities'}]);
			my @masses=split(';',$parameters[$msmsColNum{'Masses'}]);
			foreach my $i (0..$#masses) {
				$expIonTable{$masses[$i]}[0] = sprintf("%.4f",$intensities[$i]);
			}
			$query_desc="Scan number: $extSpectrumID; Elution: $parameters[$msmsColNum{'Retention time'}] min. (calibrated: $elutionTime min.); Raw file: $parameters[$msmsColNum{'Raw file'}]";
			last;
		}
	}
	close FILE;
}


elsif ($fileType eq 'sptxt') { # SWATH library
	open (FILE, $DatFile) || die "Unable to open $DatFile";
	my $onSeq=0;
	my @file=split('/',$DatFile);
	$query_desc="File: $file[-1] \t";
	unless ($call eq 'lib'){
		if ($varModsString){
			my $dbh=&promsConfig::dbConnect;
			my $sthModification =$dbh->prepare("SELECT M.ID_MODIFICATION,POS_STRING,MONO_MASS FROM MODIFICATION M,PEPTIDE_MODIFICATION IM WHERE M.ID_MODIFICATION=IM.ID_MODIFICATION AND ID_PEPTIDE=$query_number ORDER BY POS_STRING");
			$sthModification->execute;
			my @aaSeqPep=split(//,$sequence);
			while (my ($modID,$pos,$mass)=$sthModification->fetchrow_array){
				my @modPos=split(/\./,$pos);
				foreach my $position (@modPos){
					$aaSeqPep[$position-1].='['.sprintf("%0.f",$massValueAA{$aaSeqPep[$position-1]}+$mass).']';
				}
			}
			my $pepSeq=join('',@aaSeqPep);
			$peptideSequence=$pepSeq;
			$dbh->disconnect;
		}
		else{
			$peptideSequence=$sequence;
		}
	}
	my $peptideSequenceQuote=quotemeta($peptideSequence);
	if ($peptideSequenceQuote){
		while (my $line=<FILE>) {
			if ($line=~/^Name:/ && $line=~/$peptideSequenceQuote\/$charge/) {
				$onSeq=1;
			}
			elsif ($line=~/^Name:/ && $line!~/$peptideSequenceQuote\/$charge/) {
				$onSeq=0;
			}
			if ($onSeq && $line=~/PrecursorMZ: (\d+\.?\d*)/){
				$massObs=$1;
			}
			if ($onSeq && $line=~/iRT=\d+\.?\d*,(\d+\.?\d*),\d+\.?\d*/) {
				$elutionTime=sprintf("%0.2f",$1);
			}
			if (substr($line,0,2)=~/\d+/ && $onSeq==1) {
				my @massValues=split(/\t/,$line);
				my $mz=$massValues[0];
				my $intens=$massValues[1];
				$expIonTable{$mz}[0] = sprintf("%.4f",$intens);
			}
		}
		$massExp=($massObs-1.007825032)*$charge unless $massExp;

		unless (%expIonTable){
			$checkRefSpectrum=0;
			my $matchColor=$colorB;
			print header(-charset=>'utf-8'); warningsToBrowser(1);
			print "<BR><CENTER><FONT class=\"title3\">No reference MS/MS fragmentation of <FONT color=\"$matchColor\">$peptideSequence</FONT> $varModsString found.</FONT><BR></CENTER>";
		}
	}
}
elsif ($fileType eq 'SWATH.PKV' || $fileType eq 'SKYLINE.CSV' || $fileType eq 'OPENSWATH.TSV' || $fileType eq 'SPECTRONAUT.XLS') {
	open (FILE, $DatFile) || die "Unable to open $DatFile";
	while (my $line=<FILE>) {
		next if $.==1; # skip headers
		my ($peptideID,$mz,$fragCharge,$ionType,$residue,$area)=split(/!/,$line);
		next if (!$area || $area !~ /[\d\.]/i);
		if ($peptideID==$objectID) {
            $expIonTable{$mz}[0]=sprintf("%.4f",$area);
        }
    }
	$query_desc="Elution: $elutionTime min.; Raw file: $wiffFile";
}
elsif ($fileType=~/^TDM\./) {
	open(IN,"<$DatFile") or die ("open : $!");
	my (@mz,@intensity);
	my ($matchSpectreMZ,$matchSpectreI,$matchMZ,$matchIntensity)=(0,0,0,0);
	while (<IN>) {
        if($_=~/$scan\.spectrum/ && $_=~/GAML:Xdata/){
			$matchSpectreMZ=1;
		}elsif($_=~/\/GAML:Xdata/){
			$matchSpectreMZ=0;
		}elsif($_=~/$scan\.spectrum/ && $_=~/GAML:Ydata/){
			$matchSpectreI=1;
		}elsif($_=~/\/GAML:Ydata/){
			$matchSpectreI=0;
		}elsif($matchSpectreMZ && $_=~/<GAML:values/){
			$matchMZ=1;
		}elsif($matchSpectreI && $_=~/<GAML:values/){
			$matchIntensity=1;
		}elsif($matchSpectreI && $_=~/<\/GAML:values/){
			$matchIntensity=0;
		}elsif($matchSpectreMZ && $_=~/<\/GAML:values/){
			$matchMZ=0;
		}
		if ($matchMZ && $_!~/GAML/) {
			my @mass=split(/\s/,$_);
			foreach my $mass (@mass){
				push @mz,$mass;
			}
		}
		if ($matchIntensity && $_!~/GAML/) {
			my @intens=split(/\s/,$_);
			foreach my $intensity (@intens){
				push @intensity,$intensity;
			}
        }
		last if (@intensity && @mz && $matchIntensity==0 && $matchMZ==0);
    }
    close IN;

	my $i=0;
	foreach my $intensity (@intensity){
		next unless $intensity;
		$expIonTable{$mz[$i]}[0]=$intensity;
		$i++;
	}


	####### with Twig (with handlers) #####
	#my $twig=new XML::Twig(Twig_handlers => {'group/group/GAML:trace' => \&group});
	#$twig->parsefile($DatFile);

	##### with SAX #####
	#my $parser = XML::SAX::ParserFactory->parser(Handler => TDMHandler->new($DatFile,$scan,\%expIonTable) );
	#$parser->parse_uri($DatFile);


	##### with Twig (without handlers) #####
#	my $twig=new XML::Twig;
#	$twig->parsefile($DatFile);
#	my $root=$twig->root;
#	foreach my $group ($root->children){
#		next unless $group->att('id') && $group->att('id')==$scan;
#		foreach my $group2 ($group->children){
#			next unless $group2->att('label')=~/fragment ion mass spectrum/;
#            foreach my $gaml_trace ($group2->children){
#				next unless $gaml_trace->att('label')=~/$scan.spectrum/;
#				my @mz;
#				foreach my $mz (split(/\s/,$gaml_trace->field('GAML:Xdata'))){
#					next unless $mz;
#					push @mz,$mz;
#				}
#				my $i=0;
#				foreach my $intensity (split(/\s/,$gaml_trace->field('GAML:Ydata'))){
#					next unless $intensity;
#					#print
#					$expIonTable{$mz[$i]}[0]=$intensity;
#					$i++;
#				}
#			}
#		}
#	}

	##### with XML Simple #####
	#my $xml=XML::Simple->new(KeepRoot=>1);
	#my $xmlData=$xml->XMLin($DatFile);
	#foreach my $group (@{$xmlData->{'bioml'}->{'group'}}){
	#	if (defined $group->{'id'} && $group->{'id'}=~/$scan/) {
	#		foreach my $test2 (@{$group->{'group'}}){
	#			if ($test2->{'label'}=~/fragment ion mass spectrum/) {
	#				my @intensity;
	#				foreach my $intensity (split(/\s/,$test2->{'GAML:trace'}->{'GAML:Ydata'}->{'GAML:values'}->{'content'})){
	#					next if $intensity eq '';
	#					push @intensity,$intensity;
	#				}
	#				my $i=0;
	#				foreach my $mz (split(/\s/,$test2->{'GAML:trace'}->{'GAML:Xdata'}->{'GAML:values'}->{'content'})){
	#					next if $mz eq '';
	#					$expIonTable{$mz}[0]=$intensity[$i];
	#					$i++;
	#				}
	#			}
	#		}
	#	}
	#}
	$query_desc="Scan number: $scan; Elution: $elutionTime min.; Raw file: $wiffFile";
}

###################################
####>recovering  variable Mods<####
###################################
my %varModsTable;
if (defined ($varModsString)) {
	$varModsString =~s/^\s\+\s// ;
	foreach my $modif (split(/\s\+\s/, $varModsString)) {
		if ($modif=~ /(.+)\s\(([^:]+):([\d\.]+)\)/) {
			my $modifName = $1;
			my $modiCompl = $2;
			my @modifsPos = split(/\./,$3);
			foreach my $modifPos (@modifsPos) {
				# Changed on 04/09/2013
				#my $vmodCompString=($fileType eq "PARAGON.XML")? "$modifName":"$modifName ($modiCompl)";
				#if (!defined($massValueMods{"$vmodCompString"})){$warningString .= "Warning: Unhandled $vmodCompString <BR> \n"; next;}
				#$varModsTable{$modifPos} = "$vmodCompString";
				if (!defined($massValueMods{"$modifName"})){$warningString .= "Warning: Unhandled $modifName! <BR> \n"; next;}
				$varModsTable{$modifPos} = "$modifName";
			}
		}
		#elsif ($modif=~ /(\w+)\s\(N\-term\)/) {
		elsif ($modif=~ /(.+)\s\((.*N-term.*)\)/i) { # eg. (N-term), (Protein N-term), ...
			#$varModsTable{0} = "$1 ($2)";
			$varModsTable{0} = $1;
		}
		elsif ($modif ne "") {
			$warningString .= "Warning: Unhandled $modif!<BR>\n"; #error with this modification, theorical calculation will not be good !!!
		}
	}
}

####################
####>processing<####
####################

################################
####>fragmentation settings<####
################################

my %ionSerieN_term;
foreach my $fragment (@{$fragmentClassif{"N_term"}}) {
	if ($fragmentationRules{$fragment} && $fragmentationRules{$fragment}>0) {
		$ionSerieN_term{$fragment} = $fragmentDef{$fragment};
		foreach my $neutral (@{$fragmentClassif{"neutral_loss"}}) {#° or *
			if ($fragmentationRules{$fragment.$neutral} && $fragmentationRules{$fragment.$neutral}>0) {
				$ionSerieN_term{$fragment.$neutral} = $fragmentDef{$fragment} + $fragmentDef{$neutral};
			}
		}
	}
}
my %ionSerieC_term;
foreach my $fragment (@{$fragmentClassif{"C_term"}}) {
	if ($fragmentationRules{$fragment} && $fragmentationRules{$fragment}>0) {
		$ionSerieC_term{$fragment} = $fragmentDef{$fragment};
		foreach my $neutral (@{$fragmentClassif{"neutral_loss"}}) { #° or *
			if ($fragmentationRules{$fragment.$neutral} && $fragmentationRules{$fragment.$neutral}>0) {
				$ionSerieC_term{$fragment.$neutral} = $fragmentDef{$fragment} + $fragmentDef{$neutral};
			}
		}
	}
}
my %ionSerieIntern;
foreach my $fragment (@{$fragmentClassif{"intern"}}) {
	if ($fragmentationRules{$fragment} && $fragmentationRules{$fragment}>0) {
		$useInternalFragment = 1;
		$ionSerieIntern{$fragment} = $fragmentDef{$fragment};
		foreach my $neutral (@{$fragmentClassif{"neutral_loss"}}) { #° or *
			if ($fragmentationRules{$fragment.$neutral} && $fragmentationRules{$fragment.$neutral}>0) {
				$ionSerieIntern{$fragment.$neutral} =  $fragmentDef{$fragment} + $fragmentDef{$neutral};
			}
		}
	}
}

my %chargeSerie ;
$chargeSerie{'+'} = 1 if ($fragmentationRules{'+'} && $fragmentationRules{'+'}>0);
$chargeSerie{'2+'} = 2 if ($fragmentationRules{'2+'} && $fragmentationRules{'2+'}>0);
if ($fragmentationRules{'3+'} && $fragmentationRules{'3+'}>0){
	# Displaying all charges <= parent ion charge
	foreach(3..$charge){
		$chargeSerie{"$_+"} = $_;
	}
}


###########################
####>fetching Sequence<####
###########################
my @arraySequence ;
my $allowFragOxiM = 0;
my $allowFragPhospho = 0;
$arraySequence[0] ="N_term";
for (my $i=0; $i<length($sequence); $i++) {
	my $localAA =substr($sequence,$i,1);
	$arraySequence[$i+1]=$localAA;
##Uncommented by PP 26/06/14
	$allowFragOxiM = 1 if (defined($varModsTable{$i+1}) && $varModsTable{$i+1} eq "Oxidation (M)");
	#$allowFragPhospho = 1 if defined ($varModsTable{$i+1}) && ($varModsTable{$i+1} =~ /Phospho\s\(ST/ );
	$allowFragPhospho = 1 if (defined($varModsTable{$i+1}) && $varModsTable{$i+1} =~ /Phospho/ && $arraySequence[$i+1] =~ /S|T/);
	$allowVARModif{"N_term"}{"PHOS"}[$i+1] = $allowFragPhospho;
	$allowVARModif{"N_term"}{"OxiM"}[$i+1] = $allowFragOxiM;
}

push (@arraySequence, "C_term");
$allowFragOxiM = 0;
$allowFragPhospho = 0;
for (my $i=length($sequence)-1; $i>=0 ; $i--) {
	$allowFragOxiM = 1 if (defined($varModsTable{$i+1}) && $varModsTable{$i+1} eq "Oxidation (M)");
	#$allowFragPhospho = 1 if (defined($varModsTable{$i+1}) && $varModsTable{$i+1} =~ /Phospho\s\(ST/);
	$allowFragPhospho = 1 if (defined($varModsTable{$i+1}) && $varModsTable{$i+1} =~ /Phospho/ && $arraySequence[$i+1] =~ /S|T/); # Due to 04/09/2013, it had to be changed
	$allowVARModif{"C_term"}{"PHOS"}[$i+1] = $allowFragPhospho;
	$allowVARModif{"C_term"}{"OxiM"}[$i+1] = $allowFragOxiM;
}

###########################################
####>fetching N-terminal fragmentation<####
###########################################
my %massTable ; #  $massTable{ordre}{type de fragmentation}{charge}[]   ; [0] th mass, [2] infomatch [3] exp intensite [4] exp moz
my $allowLossH2O = 0;
my $allowLossNH3 = 0;
my $peptideStartMass=0;
for (my $i=1 ;$i<$#arraySequence ;$i++) {
	$peptideStartMass += $massValueAA{$arraySequence[$i]} ;
	$peptideStartMass  += $massValueMods{$varModsTable{$i}} if defined($varModsTable{$i}) ;
	$peptideStartMass  += $massValueFixedMods{$arraySequence[$i]} if defined $massValueFixedMods{$arraySequence[$i]};
	if ($i == 1){ # N-term
		$peptideStartMass += $massValueAA{$arraySequence[0]} ;
		$peptideStartMass  += $massValueMods{$varModsTable{0}} if defined($varModsTable{0}) ;
		$peptideStartMass  += $massValueFixedMods{$arraySequence[0]} if defined $massValueFixedMods{$arraySequence[0]};
	}
	if ($i == $#arraySequence-1) { # C-term
		#$peptideStartMass += $massValueAA{$arraySequence[$i+1]} ;
		$peptideStartMass  += $massValueMods{$varModsTable{$i+1}} if defined($varModsTable{$i+1}) ;
		$peptideStartMass  += $massValueFixedMods{$arraySequence[$i+1]} if defined $massValueFixedMods{$arraySequence[$i+1]};
	}
	#next if ($arraySequence[$i] eq "C_term" || $arraySequence[$i] eq "N_term");
	$allowLossH2O = 1 if ($allowLossH2O == 0 && $arraySequence[$i] =~ /[STED]/);
	$allowLossNH3 = 1 if ($allowLossNH3 == 0 && $arraySequence[$i] =~ /[NQRK]/);
	foreach my $fragment (keys (%ionSerieN_term)) {
		next if (($allowLossH2O == 0 && $fragment=~/°/) || ($allowLossNH3 == 0 && $fragment=~/\*/));
		foreach my $charge (keys %chargeSerie) {
			$massTable{$i}{$fragment}{$charge}[0] = ($peptideStartMass + $ionSerieN_term{$fragment} + $chargeSerie{$charge}*$massValueAtom{"H"})/$chargeSerie{$charge};
			foreach my $varMod (keys (%{$allowVARModif{"N_term"}})) {
				if ($allowVARModif{"N_term"}{$varMod}[$i] == 1) {
					$massTable{$i}{"$fragment-$varMod"}{$charge}[0] = $massTable{$i}{$fragment}{$charge}[0] - $massValueFragMod{$varMod}/$chargeSerie{$charge};
				}
			}
		}
	}
	if ($i == $#arraySequence-1) { # C-term
		$peptideStartMass += $massValueAA{$arraySequence[$i+1]} ;
		#$peptideStartMass  += $massValueMods{$varModsTable{$i+1}} if defined($varModsTable{$i+1}) ;
		#$peptideStartMass  += $massValueFixedMods{$arraySequence[$i+1]} if defined $massValueFixedMods{$arraySequence[$i+1]};
	}
}
my $peptideMass = $peptideStartMass;

###########################################
####>fetching C-terminal fragmentation<####
###########################################
my $peptideEndmass =0;
$allowLossH2O = 0;
$allowLossNH3 = 0;
for (my $i=$#arraySequence-1 ;$i>0 ;$i--)  {
	my $Endpos = (@arraySequence-1)-$i;
	$peptideEndmass += $massValueAA{$arraySequence[$i]} ;
	$peptideEndmass += $massValueMods{$varModsTable{$i}} if defined($varModsTable{$i}) ;
	$peptideEndmass += $massValueFixedMods{$arraySequence[$i]} if defined($massValueFixedMods{$arraySequence[$i]});
	if ($i==$#arraySequence-1) { # C-term
		$peptideEndmass += $massValueAA{$arraySequence[$i+1]};
		$peptideEndmass += $massValueMods{$varModsTable{$i+1}} if defined($varModsTable{$i+1});
		$peptideEndmass += $massValueFixedMods{$arraySequence[$i+1]} if defined($massValueFixedMods{$arraySequence[$i+1]});
	}
	if ($i==1) { # N-term
		#$peptideEndmass += $massValueAA{$arraySequence[0]};
		$peptideEndmass += $massValueMods{$varModsTable{0}} if defined($varModsTable{0}) ;
		$peptideEndmass += $massValueFixedMods{$arraySequence[0]} if defined($massValueFixedMods{$arraySequence[0]});
	}
	#next if ($arraySequence[$i] eq "C_term") || ($arraySequence[$i] eq "N_term") ;
	$allowLossH2O = 1 if ($allowLossH2O == 0 && $arraySequence[$i] =~ /[STED]/);
	$allowLossNH3 = 1 if ($allowLossNH3 == 0 && $arraySequence[$i] =~ /[NQRK]/);
	foreach my $fragment (keys %ionSerieC_term) {
		next if (($allowLossH2O == 0 && $fragment=~/°/) || ($allowLossNH3 == 0 && $fragment=~/\*/));
		foreach my $charge (keys %chargeSerie) {
			$massTable{$Endpos}{$fragment}{$charge}[0] = ($peptideEndmass+ $ionSerieC_term{$fragment} + $chargeSerie{$charge}*$massValueAtom{"H"})/$chargeSerie{$charge};
			foreach my $varMod (keys %{$allowVARModif{"C_term"}}) {
				if ($allowVARModif{"C_term"}{$varMod}[$i] == 1) {
					$massTable{$Endpos}{"$fragment-$varMod"}{$charge}[0] = $massTable{$Endpos}{$fragment}{$charge}[0] - $massValueFragMod{$varMod}/$chargeSerie{$charge};
				}
			}
		}
	}
	if ($i==1) { # N-term
		$peptideEndmass += $massValueAA{$arraySequence[0]};
		#$peptideEndmass += $massValueMods{$varModsTable{0}} if defined($varModsTable{0}) ;
		#$peptideEndmass += $massValueFixedMods{$arraySequence[0]} if defined($massValueFixedMods{$arraySequence[0]});
	}
}
#for (my $i=@arraySequence-1 ;$i>=0 ;$i--)  {
#	my $Endpos = (@arraySequence-1)-$i ;
#	print '<BR>$Endpos=',$Endpos;
#	$peptideEndmass += $massValueAA{$arraySequence[$i]} ;
#	$peptideEndmass += $massValueMods{$varModsTable{$i}} if defined($varModsTable{$i}) ;
#	$peptideEndmass += $massValueFixedMods{$arraySequence[$i]} if defined $massValueFixedMods{$arraySequence[$i]} ;
#	print '<BR>$massValueAA{$arraySequence[$i]} = ',$massValueAA{$arraySequence[$i]};
#	print '<BR>$varModsTable{$i}=',$varModsTable{$i};
#	print '<BR>$massValueMods{$varModsTable{$i}}=',$massValueMods{$varModsTable{$i}};
#	print '<BR>$massValueFixedMods{$arraySequence[$i]}=',$massValueFixedMods{$arraySequence[$i]};
#	print '<BR>$peptideEndMass=',$peptideEndmass," i=",$i," seq=",$arraySequence[$i],'<BR><BR>';
#	next if ($arraySequence[$i] eq "C_term") || ($arraySequence[$i] eq "N_term") ;
#	$allowLossH2O = 1 if ($allowLossH2O == 0) && ($arraySequence[$i] =~ /[STED]/) ;
#	$allowLossNH3 = 1 if ($allowLossNH3 == 0) && ($arraySequence[$i] =~ /[NQRK]/) ;
#	foreach my $fragment (keys (%ionSerieC_term)) {
#		next if ( (($allowLossH2O == 0)&&($fragment=~/°/)) || (($allowLossNH3 == 0)&&($fragment=~/\*/)) );
#		foreach my $charge (keys (%chargeSerie)) {
#			$massTable{$Endpos}{$fragment}{$charge}[0] = ($peptideEndmass+ $ionSerieC_term{$fragment} + $chargeSerie{$charge}*$massValueAtom{"H"})/$chargeSerie{$charge};
#			print '<BR>massTable',$massTable{$Endpos}{$fragment}{$charge}[0] ;
#			foreach my $varMod (keys (%{$allowVARModif{"C_term"}})) {
#				if ($allowVARModif{"C_term"}{$varMod}[$i] == 1 ) {
#				$massTable{$Endpos}{"$fragment-$varMod"}{$charge}[0] = $massTable{$Endpos}{$fragment}{$charge}[0] - $massValueFragMod{$varMod}/$chargeSerie{$charge} }
#			}
#		}
#	}
#}

#########################################
####>Fetching internal fragmentation<####
#########################################
my %massInternTable ;
if ($useInternalFragment ==1) {
	for (my $j=2; $j<@arraySequence; $j++) {
		my $fragmentStartMass=0;
		my $fragmentSequence='';
		for (my $i=$j; $i<@arraySequence; $i++) {
			$fragmentStartMass += $massValueAA{$arraySequence[$i]} ;
			$fragmentStartMass  += $massValueMods{$varModsTable{$i}} if defined($varModsTable{$i}) ;
			$fragmentStartMass  += $massValueFixedMods{$arraySequence[$i]} if defined($massValueFixedMods{$arraySequence[$i]});
			$fragmentSequence .= $arraySequence[$i];
			last if $fragmentStartMass>$maxFragmentMass;
			next if ($arraySequence[$i] eq "C_term" || $arraySequence[$i] eq "N_term");
			next if defined($massInternTable{$fragmentSequence}) ;
			next if $i==$j;
			foreach my $fragSerie (keys %ionSerieIntern) {
				foreach my $charge (keys %chargeSerie) {
					$massInternTable{$fragmentSequence}{$fragSerie}{$charge}[0] = ($fragmentStartMass + $ionSerieIntern{$fragSerie} + $chargeSerie{$charge}*$massValueAtom{"H"})/$chargeSerie{$charge};
				}
			}
		}
	}
}

##########################################
####>Looking for match in ionic serie<####
##########################################

####>Setting min intensity for match<####
my $higtIntensity=0;
foreach my $mz (keys %expIonTable) {
	$higtIntensity = $expIonTable{$mz}[0] if $expIonTable{$mz}[0]>$higtIntensity;
}
my $minLevel2Match = $intPercent*$higtIntensity;

#########creating @limitedExpIon
my @limitedExpIon;
foreach my $mz (keys %expIonTable) {
	if ($expIonTable{$mz}[0]>$minLevel2Match) {
		push @limitedExpIon,$mz;
	}
}

####>Looking for matched fragmentation in ionic serie<####
foreach my $number (keys %massTable) {
	foreach my $modif (keys %{$massTable{$number}}) {
		foreach my $charge (keys %{$massTable{$number}{$modif}}) {
			foreach my $mz (@limitedExpIon) {
				if (abs($massTable{$number}{$modif}{$charge}[0]-$mz)<$matchInterval)  {
					$massTable{$number}{$modif}{$charge}[2] = 1;  # %massTable [2] infomatch
					if (!$massTable{$number}{$modif}{$charge}[3] || $massTable{$number}{$modif}{$charge}[3]<$expIonTable{$mz}[0]) {
						$massTable{$number}{$modif}{$charge}[3]=$expIonTable{$mz}[0];
						$massTable{$number}{$modif}{$charge}[4]=$mz;
					}
					#my $labelString = "$modif$number";
					#$labelString .= "$charge" if $charge>1 ;
					#$expIonTable{$mz}[1].= "$labelString," ;
				}
			}
		}
	}
}
####>Looking for matched fragmentation in  Internal Fragmentation ionic serie<####
if ($useInternalFragment==1) { #$massInternTable{$fragmentSequence}{'yb'}{$charge}[0]
	foreach my $fragmentSequence (keys %massInternTable) {
		foreach my $frag_type (keys %{$massInternTable{$fragmentSequence}}) {
			foreach my $charge (keys %{$massInternTable{$fragmentSequence}{$frag_type}}) {
				foreach my $mz (@limitedExpIon) {
					if (abs($massInternTable{$fragmentSequence}{$frag_type}{$charge}[0]-$mz)<$matchInterval) {
						$massInternTable{$fragmentSequence}{$frag_type}{$charge}[2] = 1;  # %massTable [2] infomatch
						if (!$massInternTable{$fragmentSequence}{$frag_type}{$charge}[3] || $massInternTable{$fragmentSequence}{$frag_type}{$charge}[3]<$expIonTable{$mz}[0]) {
							$massInternTable{$fragmentSequence}{$frag_type}{$charge}[3]=$expIonTable{$mz}[0];
							$massInternTable{$fragmentSequence}{$frag_type}{$charge}[4]=$mz;
						}
					}
				}
			}
		}
	}
}

##########>Fetching %expIonTable,  from %massTable with priority<##############
foreach my $number (keys %massTable) {
	foreach my $modif (keys %{$massTable{$number}}) {
		my $parentModif = ($modif=~/(.*)-/)? $1 : $modif;
		foreach my $charge (keys %{$massTable{$number}{$modif}}) {
			next unless $massTable{$number}{$modif}{$charge}[4];
			#my $labelString = "$modif$number";
			my $labelString = "$modif($number)"; # PP 16/04/14
			$labelString .= "++" if $charge eq "2+";
			#$labelString = " " if $fragmentationRules{$parentModif} !=2;
			$labelString = " " if ($fragmentationRules{$parentModif} != 2 || $fragmentationRules{$charge} != 2); # modified by PP

			if (!$expIonTable{$massTable{$number}{$modif}{$charge}[4]}[1]) {
				$expIonTable{$massTable{$number}{$modif}{$charge}[4]}[1] = "$labelString";
			}
			elsif ($modif eq "y" || $modif eq "b") { #priority to y and b, for lsmp lab !
				if ($expIonTable{$massTable{$number}{$modif}{$charge}[4]}[1] =~ /y|b\d+/)  {
					$expIonTable{$massTable{$number}{$modif}{$charge}[4]}[1] .= ",$labelString";
					$expIonTable{$massTable{$number}{$modif}{$charge}[4]}[1] =~ s/(.+),(y\d.*)/$2,$1/; #place y in first position
				}
				else {
					$expIonTable{$massTable{$number}{$modif}{$charge}[4]}[1] = "$labelString";
				}
			}
			elsif ($expIonTable{$massTable{$number}{$modif}{$charge}[4]}[1] =~ /y|b\d+/) { #priority to y and b, for lsmp lab !
				next;
			}

			#else {
			#	next;
			#}
		}
	}
}


##########>Fetching %expIonTable,  from %massInternTable <##############
if ($useInternalFragment==1) { #$massInternTable{$fragmentSequence}{'yb'}{$charge}[0]
	foreach my $fragmentSequence (keys %massInternTable) {
		foreach my $frag_type (keys %{$massInternTable{$fragmentSequence}}) {
			foreach my $charge (keys %{$massInternTable{$fragmentSequence}{$frag_type}}) {
				next unless $massInternTable{$fragmentSequence}{$frag_type}{$charge}[4];
				if (!$expIonTable{$massInternTable{$fragmentSequence}{$frag_type}{$charge}[4]}[1]) {
					$expIonTable{$massInternTable{$fragmentSequence}{$frag_type}{$charge}[4]}[1] = ($fragmentationRules{$frag_type}==2)? $fragmentSequence : " ";
				}
			}
		}
	}
}


####>Looking for matched molecular pic in  ionic serie<####
foreach my $mz (@limitedExpIon) {
	if ($massObs > $mz-$matchInterval && $massObs < $mz+$matchInterval) {
		$expIonTable{$mz}[1]="M";
	}
}

###################################################
####>Creating files to by send to msms_gif.cgi<####
###################################################
my $file;

if ($useSpecApp==0) {
	if ($call eq 'lib'){
        unlink glob "$promsPath{tmp}/lib_$userID*.pgf" if glob "$promsPath{tmp}/lib_$userID*.pgf"; # cleaning previous files
		$file="$promsPath{tmp}/lib_$userID"."_$query_number"."_$rank.pgf";
    }
    else{
		if (defined($objectID)) { #analytic spectrum
			unlink glob "$promsPath{tmp}/qry_$userID*.pgf" if glob "$promsPath{tmp}/qry_$userID*.pgf"; # cleaning previous files
			$file="$promsPath{tmp}/qry_$userID"."_$query_number"."_$rank.pgf";
		}
		else { #refSpectrum
			unlink glob "$promsPath{tmp}/ref_$userID*.pgf" if glob "$promsPath{tmp}/ref_$userID*.pgf"; # cleaning previous files
			$file="$promsPath{tmp}/ref_$userID"."_$query_number"."_$rank.pgf";
		}
	}
	open (DESC,">$file");
	foreach my $mz (keys %expIonTable) {
		$defMinValue=$mz if $mz <= $defMinValue;
		$defMaxValue=$mz if $mz >= $defMaxValue;
		if (!defined($mzIntervall) || $mz>$minValue && $mz<$maxValue) {
			my $ionSerie = "\n ionserie=$mz:$expIonTable{$mz}[0]";
			$ionSerie .= ":$expIonTable{$mz}[1]" if defined $expIonTable{$mz}[1];
			print DESC $ionSerie;
		}
	}

	close DESC;
}


######################################
######>Strings for SVG graph<#########
######################################
my $IonString = "";
if ($useSpecApp==1) {
	foreach my $mz (keys (%expIonTable)) {
		$defMinValue=$mz if $mz <= $defMinValue;
		$defMaxValue=$mz if $mz >= $defMaxValue;
		if (!defined($mzIntervall) || $mz>$minValue && $mz<$maxValue) {
			$IonString.= sprintf "%.2f:%.2f",$mz,$expIonTable{$mz}[0];
			if (defined($expIonTable{$mz}[1])) {
				my $ionLabel = $expIonTable{$mz}[1];
				$ionLabel=~ s/-OxiM\((\d+)\)/\($1\)-64/g; #replacing neutral loss Ox M by -64
				$ionLabel=~ s/-PHOS\((\d+)\)/\($1\)-98/g; #replacing neutral loss Phospho  by -98
				#$ionLabel=~ s/(,?)([^,]+?)\+\+/$1\($2\)\+\+/g; #setting (  ) if ++  # COMMENTED PP 16/04/14
				$IonString .= "=$ionLabel";
			}
			$IonString .=";";
		}
	}
}
$defMinValue=int($defMinValue);
$defMaxValue=int($defMaxValue+1);

######################
#######> HTML <#######
######################
if ($checkRefSpectrum){
	print header(-charset=>'utf-8'); warningsToBrowser(1);
}
#head and javascript
print qq
|<HEAD>
<TITLE>Peptide View</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
|;
if ($useSpecApp) {
	my $unifYStrg=(abs($refSpectrum)==1)? "\n						onViewChange:parent.unifyViews," : '';
	print qq
|<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/spectrumPlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
<SCRIPT type="text/javascript">
var SP;
window.onload=function() {
    SP=new spectrumPlot({div:'spDIV',
						width:$spectrumWidth-100,height:$spectrumHeight-40,
						tolerance:$matchInterval,
						isRef:$isRef,
						isVert:$isVert,$unifYStrg
						ionSeries:'$IonString'
						});
|;
	if (abs($refSpectrum)==1) { # There is a reference
		print qq
|	var okRange=false;
	var waitCount=0;
	(function waitLoop() {
		setTimeout(function() {
			waitCount++;
			okRange=parent.unifyRange(SP);
			if (!okRange && waitCount<10) {waitLoop();}
			else {SP.draw();}
		},100)
	})();
|;
	}
	else {print "	SP.draw();\n"}
	print qq
|}
</SCRIPT>
|;
}
print qq
|<SCRIPT LANGUAGE="JavaScript">
var refTop=(top.promsFrame)? top : (top.opener.top.spectrumTolerance)? top.opener.top : top.opener.top.opener.top; // going to main myProMS window through Validation mode or protein quantification or sequence View 'ProteinWindow'!!!
refTop.spectrumTolerance = $matchInterval;
refTop.spectrumMinInt = $intPercent;
|;
#if ($refSpectrum <= 0 && $call ne "pep" && $call ne "seq") {
#	print qq
#|top.promsFrame.spectrumTolerance = $matchInterval;
#top.promsFrame.spectrumMinInt = $intPercent;
#|;
#}
#elsif ($refSpectrum > 0)  {
#	print qq
#|top.promsFrame.refSpecTolerance = $matchInterval;
#top.promsFrame.refSpecMinInt = $intPercent;
#|;
#}

if ($refSpectrum <= 0 && ($call eq "pep" || $call eq "val" )) {
	print qq
|function editComments() {
	var commentEditWinbows=window.open("$promsPath{cgi}/editSpectrumComments.cgi?EditComm=pep&ID=$objectID","","toolbar=no, location=no, directories=no, status=no, menubar=no, scrollbars=1, resizable=yes, top=400, left=300, width=500, height=180");
}
|;
}
elsif ($refSpectrum <= 0 && ($call eq "rank" || $call eq 'ana' || $call eq 'seq')) {
	print qq
|function editComments() {
	var commentEditWinbows=window.open("$promsPath{cgi}/editSpectrumComments.cgi?EditComm=rank&ID=$objectID&hit=$rank","","toolbar=no, location=no, directories=no, status=no, menubar=no, scrollbars=1, resizable=yes, top=400, left=300, width=500, height=180");
}
|;
}

my  $frame=($refSpectrum)? 'parent' : 'window';
if ($refSpectrum <= 0) {
	print qq
|function storeSpectrum(action) {
	// default (if no reference spectrum)
	var refSpecID=0;
	var rank=$rank;
|;
	if ($refSpectrum < 0) {
		print qq
|	// overwrite current reference spectrum
	if (action=='over') {
		var refPos=parent.currentPos;
		var refRank=refPos+1;
		if (!confirm ('Overwrite reference spectrum #'+refRank+' ?')) {
			return;
		}
		// Fetching old reference spectrum info from parent frame array
		refSpecID=parent.spectra[refPos][0];
		rank=parent.spectra[refPos][3];
	}
|;
	} # end of $refSpectrum < 0
	print qq
|$frame.location="$promsPath{cgi}/storeSpectrum.cgi?ACT=store&ID=$objectID&RANK="+rank+"&RID="+refSpecID+"&CALL=$call&TYPE=$fileType";
}
|;
}

if ($refSpectrum > 0) {
	my $specRank=$refPos+1;
	print qq
|	function deleteReference() {
	if (confirm("Delete reference spectrum #$specRank ?")) {
		var refSpecID=parent.spectra[$refPos][0];
		window.location="./storeSpectrum.cgi?ACT=delete&RID="+refSpecID;
	}
}
function editComments() {
	var commentEditWinbows=window.open("$promsPath{cgi}/editSpectrumComments.cgi?EditComm=ref&ID=$spectrumID","","toolbar=no, location=no, directories=no, status=no, menubar=no, scrollbars=1, resizable=yes, top=400, left=300, width=500, height=180");
}
|;
}
print qq
|function defaultRange(myForm) {
	myForm.minValue.value=$defMinValue;
	myForm.maxValue.value=$defMaxValue;
}
function checkForm(myForm) {
	if (!myForm.minValue.value \|\| myForm.minValue.value.match(/\\D/)) {alert("Enter a valid min. value for spectrum range"); return false;}
	else if (!myForm.maxValue.value \|\| myForm.maxValue.value.match(/\\D/)) {alert("Enter a valid max. value for spectrum range"); return false;}
	else if (myForm.maxValue.value*1 <= myForm.minValue.value*1) {alert("Enter valid interval for spectrum range"); return false;}

	if (!myForm.spHeight.value \|\| myForm.spHeight.value.match(/\\D/)) {alert("Enter a valid height for spectrum"); return false;}
	if (!myForm.spWidth.value \|\| myForm.spWidth.value.match(/\\D/)) {alert("Enter a valid width for spectrum"); return false;}

	if (!myForm.tolerance.value)  {alert("Enter tolerance value"); return false;}
	if (!myForm.miniInt.value)  {alert("Enter min. intensity value"); return false;}
	return true;
}
|;
&promsMod::popupInfo();
print qq
|</SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
|;
print "<FONT class=\"title2\" style=\"color:#DD0000\">$warningString</FONT>\n" if $warningString;

##########################
####>printing Result <####
##########################

####>Head<####
if ($refSpectrum <= 0) { # query spectrum
	my $matchColor=$colorA;
	print "<CENTER><FONT class=\"title2\">MS/MS Fragmentation of <FONT color=\"$matchColor\">$titleSequence</FONT></FONT>";
	print "<FONT class=\"title3\"> + $varModsString</FONT>" if $varModsString;
	print "<FONT class=\"title2\"> in Selected MS/MS Analysis</FONT>" if $call eq 'ana';
	print "</CENTER>\n";
	print "<BR><B>Comments:</B> ". &promsMod::HTMLcompatible($comments) ."\n" if $comments;
	#my $disabNew=((defined($numRef) && $numRef>=$maxRef ) || $disableStore || $fileType=~/\.PDM/)? 'disabled' : ''; # Modif 13/11/12
	#my $disabOver=($refSpectrum==0 || $disableStore || $fileType=~/\.PDM/)? 'disabled' : ''; # Modif 13/11/12
	my $disabNew=((defined($numRef) && $numRef>=$maxRef ) || $disableStore)? 'disabled' : '';
	my $disabOver=($refSpectrum==0 || $disableStore)? 'disabled' : '';
	print qq|<BR><TABLE border=0 cellpadding=0>|;
	if ($call ne 'lib' && $fileType ne 'MAXQUANT.DIR') {print "<TH align=\"left\" nowrap>Save spectrum as <INPUT type=\"button\" value=\"New reference\" onclick=\"storeSpectrum('new')\" $disabNew/> <B>or</B> <INPUT type=\"button\" value=\"Overwrite current\" onclick=\"storeSpectrum('over')\" $disabOver></TH><TD>&nbsp;&nbsp;</TD>";}
	print qq |<TD nowrap><INPUT type="button" value="Edit Comments" onclick="editComments()" |;if ($call eq 'lib'){print $disabOver} else{print $disableStore} print qq |/></TD>|; #if $call eq "pep" (on last TD);
}
elsif ($checkRefSpectrum) { # Drawing a Reference Spectrum
	my $matchColor=$colorB;
	print "<CENTER><FONT class=\"title3\">Reference MS/MS Fragmentation of <FONT color=\"$matchColor\">$titleSequence</FONT>";
	print " + $varModsString" if $varModsString;
	printf "&nbsp;(Score: %.1f)",$refScore if $refScore;
	print "&nbsp;&nbsp;&nbsp;Reference:";
	for (my $s=0;$s<=$maxRef;$s++) {
		if ($s==$refPos) {
			print "&nbsp;&nbsp;",$s+1,"</FONT>";
		}
		else {
			my $disab=($s > $numRef)? 'disabled' : '';
			print "&nbsp;&nbsp;<INPUT type=\"button\" value=\"",$s+1,"\" style=\"width:25px;font-weight:bold\" onclick=\"parent.loadRefSpectrum($s)\" $disab/>";
		}
	}
	print "</FONT></CENTER>\n";
	print "<B> (Spectrum recorded by '$upUser')</B>\n" if $upUser;
	print "<BR><B>Comments:</B> ". &promsMod::HTMLcompatible($comments) ."\n" if $comments;
	print "<BR>";
	print qq
|<TABLE border=0 cellpadding=0><TD nowrap><INPUT type="button" value="Delete reference" onclick="deleteReference()" $disableStore/>
<INPUT type="button" value="Edit Comments" onclick="editComments()" $disableStore/></TD>
| unless $fileType eq 'sptxt';
}
print qq
|<TD>&nbsp;&nbsp;</TD>
<TH align="left" nowrap>Display:<SELECT onchange="SP.toBlackAndWhite(this.value*1)"><OPTION value="0">Colored</OPTION><OPTION value="1">Black & White</OPTION></SELECT>
<INPUT type="button" value="Export as PNG" onclick="exportSVGtoImg('spDIV','Spectrum_$spectrumCode','./exportSVG.cgi','png')"/><INPUT type="button" value="Export as SVG" onclick="exportSVGtoImg('spDIV','Spectrum_$spectrumCode','./exportSVG.cgi','svg')"/></TH>
| if $useSpecApp && $checkRefSpectrum;
print "</TR></TABLE>\n";

####>Graphical view of ionic serie<####
my $spectreType =($refSpectrum <= 0)? "query" : "ref";
my $intervallValue =(defined($mzIntervall) && $mzIntervall==1)? $minValue."-".$maxValue:"all";
if ($useSpecApp==0) {
	print qq
|<TABLE><TR><TD>
	<IMG SRC="./msms_gif.cgi?file=$file&width=$spectrumWidth&height=$spectrumHeight&TYPE=$spectreType&scall=$intervallValue" WIDTH=$spectrumWidth HEIGHT=$spectrumHeight align=center BORDER=2 >
	</TD></TR></TABLE>
|;
}
elsif ($checkRefSpectrum) { # codebase="$promsPath{'html'}/java"
	print "<DIV id=\"spDIV\"></DIV>\n";
#	my $appletWidth =$frameWidth-50;
#	print qq
#|<TABLE border=1 cellspacing=0 cellpadding=0><TR><TD>
#	<APPLET code="SpView.class" codebase="$promsPath{java_spectrum}" WIDTH=$appletWidth HEIGHT=220 >
#	<PARAM NAME="isRef", VALUE="$isRef">
#	<PARAM NAME="isVert", VALUE="$isVert">
#	<PARAM NAME="ionSerie", VALUE="$IonString">
#	<BR><CENTER>
#	<FONT class="title3" style="color:#DD0000">There is a problem with Java on your browser.<BR>&nbsp;Fix it or disable interactive spectrum usage in your profile.&nbsp;</FONT>
#	</CENTER><BR>
#</APPLET></TD></TR></TABLE>
#|;
}

####>Analysis summary<####
if($checkRefSpectrum){
	print "$query_desc\n<BR>";
	print "Elution: $elutionTime min.&nbsp;&nbsp;&nbsp;&nbsp;\n" unless $query_desc=~/Elution/;
	printf "Mr(calc): %.5f; Mr(exp): %.5f; Mr(obs): %.5f; charge %1D+<BR>\n",$peptideMass,$massExp,$massObs,$charge;
	####>Instrument<####
	print "Acquired on <B>$instrument<B><BR>\n" if $instrument; #$refSpectrum>0

	####>Fixed modif summary<####
	foreach my $numMod (keys %FixMods) {
		print "Fixed Modification $FixMods{$numMod}[0]: $FixMods{$numMod}[1]<BR>\n";
	}
	if ($fileType eq "MASCOT.XML" && defined($fixedModsString))	{
		print "Fixed Modification: $fixedModsString<BR>\n";
	}
	####>var modif summary<####
	foreach my $modifPos (sort {$a<=>$b} keys (%varModsTable)) {
		print "$varModsTable{$modifPos} at position $arraySequence[$modifPos]";
		print " $modifPos" if $arraySequence[$modifPos] !~/_term/;
		print "<BR>\n";
	}


	####>Table header<####
	print "<BR><TABLE border=1 cellspacing=0 cellpadding=2><TR bgcolor=\"$color2\">\n";
	foreach my $modif (sort{if ($a =~ /b/ && $b !~ /b/) {return -1;} elsif ($b =~ /b/ && $a !~ /b/) {return 1;} else {$a cmp $b}} keys %ionSerieN_term) {
		foreach my $charge (sort{$a cmp $b} keys %chargeSerie) {
			print "<TH>$modif<SUP>$charge</SUP></TH>";
		}
	}
	print "<TH>#</TH><TH>aminoacid</TH><TH>#</TH>";
	foreach my $modif (sort{$a cmp $b} keys %ionSerieC_term) {
		foreach my $charge (sort{$a cmp $b} keys %chargeSerie) {
			print "<TH>$modif<SUP>$charge</SUP></TH>";
		}
	}
	print "</TR>\n";



	####>ion sequence<####
	for (my $i=0 ;$i<@arraySequence-1 ;$i++) {
		next if ($arraySequence[$i] eq "C_term" || $arraySequence[$i] eq "N_term");
		print "<TR align=center>\n";
		foreach my $modif (sort{if ($a =~ /b/ && $b !~ /b/) {return -1;} elsif ($b =~ /b/ && $a !~ /b/) {return 1 ; } else { $a cmp $b ;}} keys %ionSerieN_term) {
			foreach my $charge (sort{$a cmp $b} keys %chargeSerie) {
				my $localColor=(defined($massTable{$i}{$modif}{$charge}[2]))? $higlightColor : $whiteColor;
				#if ($massTable{$i}{$modif}{$charge}[0]) {printf "<TD bgcolor=$localColor>%.2f</TD>",$massTable{$i}{$modif}{$charge}[0];}
				if ($massTable{$i}{$modif}{$charge}[0]) {
					if ($useSpecApp==1) {
						print "<TD bgcolor=\"$localColor\"><A href=\"javascript:void(null)\" onmouseover=\"if (SP && SP.ready) {SP.highlightMZ('on',$massTable{$i}{$modif}{$charge}[0])}\" onmouseout=\"if (SP && SP.ready) {SP.highlightMZ('off',$massTable{$i}{$modif}{$charge}[0])}\">";
						printf "%.2f</A></TD>",$massTable{$i}{$modif}{$charge}[0];
					}
					else {printf "<TD bgcolor=$localColor>%.2f</TD>",$massTable{$i}{$modif}{$charge}[0];}
				}
				else {print "<TD></TD>";}
			}
		}
		my $endPos = (@arraySequence-1)- $i;

		my $aaSymbol = $arraySequence[$i];
		if ($varModsTable{$i} || ($i == 1 && $varModsTable{0})) {
			my $popupText = '';
			if ($i == 1 && $varModsTable{0}) {
				$popupText .= "$varModsTable{0}<BR>";
			}
			if ($varModsTable{$i}) {
				$popupText .= $varModsTable{$i};
			}

			$aaSymbol .= "<A onmouseover=\"popup('$popupText');\" onmouseout=\"popout();\"><FONT color=FF00FF >*</FONT></A>";
		}

		print "<TD><B>$i</B> </TD><TD><FONT color=FF00FF ><B>$aaSymbol</B></FONT></TD><TD><B>$endPos</B></TD>";
		foreach my $modif (sort{$a cmp $b} keys %ionSerieC_term) {
			foreach my $charge (sort{$a cmp $b} keys %chargeSerie) {
				my $localColor=(defined($massTable{$endPos}{$modif}{$charge}[2]))? $higlightColor : $whiteColor;
				#if ($massTable{$endPos}{$modif}{$charge}[0] == 1131.62 ) { printf "hello"; }
				#if ($massTable{$endPos}{$modif}{$charge}[0]) {printf "<TD bgcolor=$localColor>%.2f</TD>",$massTable{$endPos}{$modif}{$charge}[0];}
				if ($massTable{$endPos}{$modif}{$charge}[0]) {
					if ($useSpecApp==1) {
						print "<TD bgcolor=\"$localColor\"><A href=\"javascript:void(null)\" onmouseover=\"if (SP && SP.ready) {SP.highlightMZ('on',$massTable{$endPos}{$modif}{$charge}[0])}\" onmouseout=\"if (SP && SP.ready) {SP.highlightMZ('off',$massTable{$endPos}{$modif}{$charge}[0])}\">";
						printf "%.2f</A></TD>",$massTable{$endPos}{$modif}{$charge}[0];
					}
					else {printf "<TD bgcolor=$localColor>%.2f</TD>",$massTable{$endPos}{$modif}{$charge}[0];}
				}
				else {print "<TD></TD>";}
			}
		}
		print "</TR>\n";
	}

	print "</TABLE><BR>\n";
}
#####################################################
###############Fragment Table########################
#####################################################
if ($useInternalFragment ==1 && $checkRefSpectrum) {
	my $column_num = 4;
	my $header_motif = "<TH>Sequence</TH>";
	foreach my $fragment_name (sort {$a cmp $b} keys %ionSerieIntern) {
		$header_motif .= "<TH>$fragment_name</TH>";
	}
	my $full_string;
	for (my $i=1; $i<=$column_num; $i++) {
		$full_string .= $header_motif."\n";
	}
	print qq
|<BR>
<TABLE BORDER=1 CELLSPACING=0 CELLPADDING=2>
<TR bgcolor="$color2">
$full_string
</TR><TR>
|;
	my $count=0;
	foreach my $fragmentSequence (sort {$a cmp $b} keys %massInternTable) {
		print "<TH>$fragmentSequence</TH>";
		foreach my $fragment_name (sort {$a cmp $b} keys %ionSerieIntern) {
			my $localColor=(defined($massInternTable{$fragmentSequence}{$fragment_name}{"+"}[2]))? $higlightColor : $whiteColor;
			printf "<TD bgcolor = $localColor>%.2f</TD>",$massInternTable{$fragmentSequence}{$fragment_name}{"+"}[0];
		}
		$count++;
		if ($count ==$column_num) {
			$count=0;
			print "</TR><TR>";
		}
	}
	print qq
|</TR></TABLE>
<BR>
|;
}

####>Parameter Form<####
my $projectString=($projectID)?"<INPUT type=hidden name='PROJECT_ID' value='$projectID'>":''; # No $projectID for reference spectrum
print qq
|<FORM name="rangeForm" method="post" onsubmit="return checkForm(this)">
<INPUT type=hidden name='file' value="$DatFile">
<INPUT type= hidden name='query' value="$query_number">
<INPUT type="hidden" name="hit" value="$rank">
<INPUT type=hidden name='REF' value="$refSpectrum">
<INPUT type=hidden name='TYPE' value="$fileType">
<INPUT type=hidden name='width' value="$frameWidth">
$projectString
|;

my $ref_Pos = $refPos if defined $refPos;
$ref_Pos .="_".$numRef if (defined($ref_Pos) && defined($numRef));
print "<INPUT type=hidden name='REF_POS' value='$ref_Pos'> \n" if defined($ref_Pos);
print "<INPUT type=hidden name='ID' value='$objectID'> \n"  if defined($objectID);
print "<INPUT type=hidden name='CALL' value='$call'> \n"  if $call;
if ($call eq 'lib') {
	print "<INPUT type=hidden name='irt' value='$elutionTime'> \n"  if $elutionTime;
	print "<INPUT type=hidden name='massObs' value='$massObs'> \n"  if $massObs;
	print "<INPUT type=hidden name='massExp' value='$massExp'> \n"  if $massExp;
	print "<INPUT type=hidden name='SEQUENCE' value='$sequence'> \n"  if $sequence;
	print "<INPUT type=hidden name='libID' value='$libID'> \n"  if $libID;
}
print "<INPUT type=hidden name= 'disStore' value='$disableStore'> \n" if $disableStore;
print "<INPUT type=hidden name= 'REF_SC' value='$refScore'> \n" if $refScore;
print "<INPUT type=hidden name= 'UP_USER' value='$upUser'> \n" if $upUser;

my $usedMinValue=($minValue)? $minValue : $defMinValue;
my $usedMaxValue=($maxValue)? $maxValue : $defMaxValue;
print qq
|<FONT class="title3">Parameters for drawing spectrum :</FONT>
<TABLE border=0 bgcolor="$color2">
<TR><TH align="right">Range :</TH><TD bgcolor="$color1"><B>Min value:</B><INPUT type="text" name="minValue" value="$usedMinValue" size=6>&nbsp;&nbsp;&nbsp;&nbsp;<B>Max value:</B><INPUT type="text" name="maxValue" value="$usedMaxValue" size=6>
<INPUT type="button" value="Default range" onclick="defaultRange(document.rangeForm)"></TD></TR>
<TR><TH align="right">Size :</TH><TD bgcolor="$color1"><B>Height:</B><INPUT type="text" name="spHeight" value="$spectrumHeight" size=6> px&nbsp;&nbsp;&nbsp;&nbsp;<B>Width:</B><INPUT type="text" name="spWidth" value="$spectrumWidth" size=8> px</TD></TR>
<TR><TH align="right" valign="top">&nbsp;Peak matching :</TH>
<TD nowrap bgcolor="$color1"><B>Tolerance:</B><INPUT type=text name=tolerance value=$matchInterval size=4>Da<BR>
<B>Minimal intensity:</B><INPUT type="text" name="miniInt" value="$intPercent" size=8> (fraction of highest peak)
</TD></TR>
<TR><TH colspan=2><INPUT type="submit" name="submitParam" value="Apply"></TH></TR>
</TABLE>
</FORM>
<DIV id="divDescription" class="clDescriptionCont">
	<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">
	setPopup();
</SCRIPT>
</BODY>
</HTML>
| if $checkRefSpectrum;

sub group {
	my ($twig,$gaml_trace)=@_;
	if ($gaml_trace->att('label')=~/$scan.spectrum/) {
        my @mz;
		foreach my $mz (split(/\s/,$gaml_trace->field('GAML:Xdata'))){
			next unless $mz;
			push @mz,$mz;
		}
		my $i=0;
		foreach my $intensity (split(/\s/,$gaml_trace->field('GAML:Ydata'))){
			next unless $intensity;
			$expIonTable{$mz[$i]}[0]=$intensity;
			$i++;
		}
    }
	#if ($group->att('id') && $group->att('id')==$scan){
	#	foreach my $group2 ($group->children){
	#		next unless $group2->att('label')=~/fragment ion mass spectrum/;
	#		foreach my $gaml_trace ($group2->children){
	#			next unless $gaml_trace->att('label')=~/$scan.spectrum/;
	#			my @mz;
	#			foreach my $mz (split(/\s/,$gaml_trace->field('GAML:Xdata'))){
	#				next unless $mz;
	#				push @mz,$mz;
	#			}
	#			my $i=0;
	#			foreach my $intensity (split(/\s/,$gaml_trace->field('GAML:Ydata'))){
	#				next unless $intensity;
	#				print $intensity,"<BR>";
	#				$expIonTable{$mz[$i]}[0]=$intensity;
	#				$i++;
	#			}
	#		}
	#	}
	#}
	$twig->purge;
	return;
}


package TDMHandler; {
	my (@mz,$refExpIonTable);
	my $i=0;
	sub new {
		(my $type,my $file,my $scan,$refExpIonTable)=@_;
		my $self=bless({},$type);
		return $self;
	}

	sub start_document{
		my ($self)=@_;
	}

	sub end_document{
		my $self=$_[0];
	}

	sub start_element{
		my ($self,$element)=@_;
		$self->{'Element'}=$element->{'Name'};
		if ($element->{'Name'}=~/group/ && !$self->{'isGroup1'}) {
			$self->{'isGroup1'}=1;
			if ($element->{'Attributes'}->{'{}id'}->{'Value'} && $element->{'Attributes'}->{'{}id'}->{'Value'}==$scan) {
                $self->{'isPep'}=1;
            }
        }elsif($element->{'Name'}=~/group/ && $self->{'isGroup1'}){
			$self->{'isGroup2'}=1;
		}elsif($self->{'isPep'}){
			if($element->{'Name'}=~/GAML:Xdata/){
				if ($element->{'Attributes'}->{'{}label'}->{'Value'}=~/$scan.spectrum/) {
					$self->{'isMass'}=1;
				}else{$self->{'isMass'}=0;}

			}elsif($element->{'Name'}=~/GAML:Ydata/){
				if ($element->{'Attributes'}->{'{}label'}->{'Value'}=~/$scan.spectrum/) {
					$self->{'isIntensity'}=1;
				}else{$self->{'isIntensity'}=0;}
			}
			elsif($element->{'Name'}=~/GAML:values/){
				if($self->{'isIntensity'}){
					$self->{'isIntensityValue'}=1;
				}
				elsif($self->{'isMass'}){
					$self->{'isMassValue'}=1;
				}
			}
		}
	}

	sub end_element{
		my ($self,$element)=@_;
		if ($element->{'Name'}=~/group/ && !$self->{'isGroup2'}) {
			$self->{'isPep'}=0;
			$self->{'isGroup1'}=0;
			$i=0;
		}elsif($element->{'Name'}=~/group/ && $self->{'isGroup2'}){
			$self->{'isGroup2'}=0;
		}elsif($element->{'Name'}=~/GAML:values/){
			#next if $self->{'isTable'};
			$self->{'isIntensityValue'}=0;
			$self->{'isMassValue'}=0;
		}
		$self->{'Element'}='';
	}

	sub characters {
        my ($self, $element) = @_;
		if ($self->{'isMassValue'}) {
			foreach my $mz (split(/\s/,$element->{'Data'})){
				next if $mz eq '';
				push @mz,$mz;
			}
        }elsif($self->{'isIntensityValue'}){
			foreach my $intensity (split(/\s/,$element->{'Data'})){
				next if $intensity eq '';
				$refExpIonTable->{$mz[$i]}[0]=$intensity;
				$i++;
			}
			$self->{'isTable'}=1;
		}
    }
}1;

####>Revision history<####
# 3.2.1 Cleaner exit on missing spectrum data file (PP 23/03/18)
# 3.2.0 Minor modif to allow Spectronaut spectrum drawing (MLP 23/01/18)
# 3.1.10 Minor modif to allow DIA reference spectrum drawing (MLP 19/12/17)
# 3.1.9 Minor modif to allow OpenSwath spectrum drawing (MLP 01/09/17)
# 3.1.8 Modification to extract spectrum informations for X! Tandem files (MLP 16/02/17)
# 3.1.7 Fix uninitialized variables for SWATH.PKV (PP 20/01/17)
# 3.1.6 Update to display spectrum for X! Tandem (MLP 16/12/16)
# 3.1.5 Compatible with MaxQuant except Reference spectrum (PP 29/11/16)
# 3.1.4 Bug fix in non-Unimod modification mass detection in DAT file (PP 02/11/16)
# 3.1.3 Minor bug correction to get informations from sptxt files (MLP 18/07/2016)
# 3.1.2 Minor correction (MLP 16/06/2016)
# 3.1.1 Minor bug correction to pass $libID in "Parameters for drawing spectrum" (MLP 02/05/2016)
# 3.1.0 Modification to adapt peptide_view for Swath Libraries (MLP 11/04/2016)
# 3.0.7 Update to allow call from protein quantification raw peptide data (PP 23/09/15)
# 3.0.6 Minor modif due to split-mode that changes $msfFile (GA 27/02/15)
# 3.0.5 5 digit precision on peptide masses in summary info (PP 17/02/15)
# 3.0.4 Get parentFile from promsMod::extractSpectrumMSF (GA 28/01/15)
# 3.0.3 Black&white + resize + SVG export options for spectrum & '-PHOS/OxiM(x)' changed to '(x)-98/64' + activated for N-term (PP 26/06/14)
# 3.0.2 Interactive spectrum is exportabled as image (PP 18/06/14)
# 3.0.1 Change is ion label y7 -> y(7) (PP 16/04/14)
# 3.0.0 JAVA applet replaced by SVG. Requires spectrumPlot.js (PP 31/03/14)
# 2.8.8 Minor modif due to modif of listRanks for altnames of modifications (GA 06/03/14)
# 2.8.7 Minor modif in $allowFragPhospho test due to minor modif 2.8.4 (GA 27/11/13)
# 2.8.6 Added 'use XML::Simple' (PP 13/11/13)
# 2.8.5 Minor code rearrangement (FY 25/09/13)
# 2.8.4 Same bug correction than 2.8.2 that was not considering all the possibilities (GA 04/09/13)
# 2.8.3 Minor bug correction to not use $projectID in form for reference spectrum (GA 02/09/13)
# 2.8.2 Minor bug correction to make $massValueMods up-to-date with mods specificity (GA 23/08/13)
# 2.8.1 Fix missing PROJECT_ID parameter in range form (PP 09/07/13)
# 2.8.0 Remove VAR_MOD from script (GA 22/05/13)
# 2.7.9 Conversion to UTF-8 (PP 22/04/13)
# 2.7.8 Displaying an asterisk with popup for var mods on each AA in fragmentation table (FY 17/04/13)
# 2.7.7 Displaying Paragon Substitution in validation mode! (GA 21/03/13)
# 2.7.6 Displaying >= 3+ charges in table and spectrum (FY 19/03/13)
# 2.7.5 Uses new file path for fully validated analyses (PP 07/03/13)
# 2.7.4 Update peptide spectrum for proteome-discoverer: get $extSpectrumID for pep call (GA 28/11/12)
# 2.7.3 Minor modif to save reference spectrum for PDM files (GA 13/11/12)
# 2.7.2 Modification to get wiff file name information for PARAGON xml merged files (GA 02/10/12)
# 2.7.1 Modification to extract spectrum information for PARAGON xml files (GA 18/09/12)
# 2.7.0 Fixed missing comma and parsing rules for COM attribute in INFO_PEP (FY 20/07/12)
# 2.6.9 Skip java activation test & checks on uninitialized values (PP 18/04/12)
# 2.6.8 Fix deletion warnings when tmp files do not exist (FY 03/04/12)
# 2.6.7 Improved Reg Exp for N-term PMTs detection (PP 08/04/11)
# 2.6.6 Minor update to print well the elution time (GA 07/01/11)<BR>See 2.6.5 modification in storeAnalyses.cgi
