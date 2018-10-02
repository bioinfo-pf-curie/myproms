#!/usr/local/bin/perl -w

################################################################################
# importMaxquant.cgi       1.3.3                                               #
# Component of site myProMS Web Server                                         #
# Authors: P. Poullet, G. Arras, S. Liva (Institut Curie)                      #
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

use strict;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use POSIX qw(strftime);  # Core module
#use IO::Uncompress::Gunzip qw(gunzip);
use XML::Simple;
use File::Path qw(rmtree); # remove_tree
use File::Copy qw(copy move); # Core module
use File::Spec::Functions qw(splitpath); # Core module
use promsConfig;
use promsMod;
#exit;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my ($color1,$color2)=&promsConfig::getRowColors;
my $MAX_DB=3;

my $experimentID=&promsMod::cleanNumericalParameters(param('ID'));
my $projectID=&promsMod::cleanNumericalParameters(param('id_project'));
my $action=param('ACT') || 'form';

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

print header(-'content-encoding'=>'no',-charset=>'utf-8'); warningsToBrowser(1);

if ($action eq 'form') {

	my (%DBlist,%isCrapDB);

	my $sthDB=$dbh->prepare("SELECT D.ID_DATABANK,D.NAME,VERSION_NAME,FASTA_FILE,DT.NAME,IS_CRAP FROM DATABANK D,DATABANK_TYPE DT WHERE D.ID_DBTYPE=DT.ID_DBTYPE AND USE_STATUS='yes'");
	$sthDB->execute;
	while (my ($dbID,$name,$version,$fastaFile,$dbankType,$isCrap)=$sthDB->fetchrow_array) {
		my $dbSource;
		if ($fastaFile=~/:/) {
			my ($mascotServer,$dbankDir,$fileName)=split(':',$fastaFile);
			$dbSource=$mascotServer;
			$version=$dbankDir;
		}
		else {
			if (!-e "$promsPath{banks}/db_$dbID/$fastaFile") {next;}
			$dbSource='Local';
		}
		$DBlist{$dbSource}{$dbID}=$name;
		$DBlist{$dbSource}{$dbID}.=" ($version)" if $version;
		$DBlist{$dbSource}{$dbID}.=" [$dbankType]";
		$isCrapDB{$dbID}=$isCrap;
	}
	$sthDB->finish;
	my $databaseString="<OPTION selected value=\"\">-= Select =-</OPTION>\n";
	my $databaseContString="<OPTION selected value=\"\">*Contaminants not searched*</OPTION>\n";

	foreach my $dbSource (sort{lc($a) cmp lc($b)} keys %DBlist) {
		$databaseString.="<OPTGROUP label=\"$dbSource:\">\n";
		foreach my $dbID (sort{lc($DBlist{$dbSource}{$a}) cmp lc($DBlist{$dbSource}{$b})} keys %{$DBlist{$dbSource}}) {
			$databaseString.="<OPTION value=\"$dbID\">$DBlist{$dbSource}{$dbID}</OPTION>\n";
			$databaseContString.="<OPTION value=\"$dbID\">$DBlist{$dbSource}{$dbID}</OPTION>\n" if $isCrapDB{$dbID};
		}
		$databaseString.="</OPTGROUP>\n";
		$databaseContString.="</OPTGROUP>\n";
	}

	$dbh->disconnect;

	print qq
|<HTML>
<HEAD>
<TITLE>Select MaxQuant files</TITLE>
<LINK rel=stylesheet href="$promsPath{html}/promsStyle.css" type=text/css>
<SCRIPT type="text/javascript">
|;
	&promsMod::browseDirectory_JavaScript if $promsPath{'shared'};
	print qq
|function updateDatabankSelection(dbNum,optIdx) {
	if (dbNum==$MAX_DB) return;
	var dbNext=dbNum+1;
	if (optIdx) { // selection
		document.getElementById('dbDIV_'+dbNext).style.display='';
		document.getElementById('db_'+dbNext).disabled=false;
	}
	else {
		for (var d=dbNext; d<=$MAX_DB; d++) {
			document.getElementById('dbDIV_'+d).style.display='none';
			document.getElementById('db_'+d).disabled=true;
			document.getElementById('db_'+d).options.selectedIndex=0;
		}
	}
}
function cancelAction() {
	window.location="./processAnalyses.cgi?ID=$experimentID";
}
function checkFileForm(importForm) {
	if (document.getElementById('formSubmit').disabled==true) return false; // in case "Enter" key is pressed after 1st submission
//	if (importForm.mqpar.files.length == 0) {
//		alert('ERROR: Select the XML parameter file for your MaxQuant analysis.');
//		return false;
//	}
	var missingFiles=false;
	if (importForm.uplArch.files.length == 0) { // nothing in upload option
		missingFiles=true;
		if (importForm.sharedDirFiles) { // shared dir available
			if (importForm.sharedDirFiles.length) { // multiple files found
				var numFiles=0;
				for (let i=0; i<importForm.sharedDirFiles.length; i++) {
					if (importForm.sharedDirFiles[i].checked) {
//console.log(importForm.sharedDirFiles[i].value);
						if (importForm.sharedDirFiles[i].value.match(/\\.(zip\|gz)\$/)) {
							missingFiles=false;
							break;
						}
						else if (importForm.sharedDirFiles[i].value.match(/(evidence\|peptides)\\.txt\$/)) {
//console.log('OK',i);
							if (++numFiles==2) {
								missingFiles=false;
								break;
							}
						}
					}
				}
			}
			else if (importForm.sharedDirFiles.checked && importForm.sharedDirFiles.value.match(/\\.(zip\|gz)\$/)) {missingFiles=false;} // single file found
		}
	}
	if (missingFiles) {
		alert('ERROR: Missing MaxQuant data file(s).');
		return false;
	}
	if (!importForm.databank1.value) {
		alert('ERROR: Select at least 1 sequence databank to be used !');
		return false;
	}
	if (!importForm.mgType.value) {
		alert('ERROR: Select a rule for Match groups generation !');
		return false;
	}
	document.getElementById('formSubmit').disabled=true;
	document.getElementById('waitDIV').style.display='block';
	return true;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Select MaxQuant file to be imported</FONT>
<BR><BR>
<FORM name="mqListForm" action="./importMaxquant.cgi" method="post" onsubmit="return checkFileForm(this);" enctype="multipart/form-data">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="ID" value="$experimentID">
<INPUT type="hidden" name="ACT" value="upload">
<TABLE>
<TR><TD colspan=2>
	<TABLE bgcolor=$color2 style="min-width:800px">
	<TR>
		<TH align=right width=250 valign=top>MaxQuant parameter file :</TH>
		<TD bgcolor=$color1 nowrap><INPUT type="file" name="mqpar"><BR>(The mqpar.xml file used for your MaxQuant analysis)</TD>
	</TR>
	<TR>
		<TH align=right valign=top>MaxQuant data files* :</TH>
		<TD bgcolor=$color1 nowrap>
		<TABLE cellpadding=0>
			<TR><TH align="right">Upload archive file**:</TH><TD><INPUT type="file" name="uplArch"></TD></TR>
|;
	if ($promsPath{'shared'}) {
		print qq
|			<TR><TH align="right">or</TH><TH></TH></TR>
			<TR><TH align="right" valign=top nowrap>Use shared directory:</TH><TD><DIV id="sharedDirDIV" style="width:600px;max-height:300px;overflow:auto;border:2px solid $color2">
|;
		&promsMod::browseDirectory_getFiles($promsPath{'shared'},'sharedDirFiles',{fileMatch=>qr/\.(txt|gz|zip)$/i}); # List of handled file extensions
		print qq
|		</DIV></TD>
	</TR>
|;
	}
	print qq
|		</TABLE>
		</TD>
	</TR>
	<TR>
		<TH align=right valign="top" nowrap>&nbsp;Databank(s) used for search :</TH>
		<TD bgcolor=$color1 nowrap>|;
	foreach my $d (1..$MAX_DB) {
		my ($dispStrg,$disabStrg)=($d==1)? ('','') : ('style="display:none"',' disabled');
		print qq
|<DIV id="dbDIV_$d" $dispStrg><B>#$d:</B><SELECT id="db_$d" name="databank$d" style="width:550px" onchange="updateDatabankSelection($d,this.options.selectedIndex)">$databaseString</SELECT></DIV>
|;
	}
	print qq
|</TD>
	</TR>
	<TR>
		<TH align=right valign="top">Contaminants :</TH>
		<TD bgcolor=$color1 nowrap><INPUT type="checkbox" name="excludeCON" value="1" onclick="document.getElementById('dbDIV_con').style.display=(this.checked)? 'none' : ''"><B>Exclude from protein list<BR><DIV id="dbDIV_con">&nbsp;<B>Databank:</B><SELECT name="databankCON" style="max-width:600px">$databaseContString</SELECT></DIV></TD>
	</TR>
	<TR>
		<TH align=right valign="top">Quantification :</TH>
		<TD bgcolor=$color1 nowrap><INPUT type="checkbox" name="protQuantif"  value="1" checked><B>Import protein quantification data</B></TD>
	</TR>
	<TR>
		<TH align=right>Match groups rule :</TH>
		<TD bgcolor=$color1 nowrap><SELECT name="mgType"><OPTION value="">-= Select =-</OPTION><OPTION value="MaxQuant">MaxQuant</OPTION><OPTION value="myProMS">myProMS</OPTION></SELECT></TD>
	</TR>
	<TR><TH colspan=2>
		<INPUT type="submit" name="save" id="formSubmit" value=" Proceed ">
		&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="cancelAction();"></TH>
	</TR>
</TABLE></TD></TR>
<TR>
	<TD valign=top nowrap>
	<DIV style="float:left">*</DIV><DIV style="float:left;height:100px">MaxQuant files from the 'combined/txt' directory:
	<UL style="margin:0">
		<LI>evidence.txt</LI>
		<LI>msms.txt <I>(Optional. For MS/MS spectrum display)</I></LI>
		<LI>peptides.txt</LI><LI>proteinGroups.txt <I>(Required for protein quantification only)</I></LI>
		<LI>parameters.txt <I>(Required only if mqpar.xml is not available)</I></LI>
		<LI>summary.txt <I>(Required only if mqpar.xml is not available)</I></LI>
	</UL>
	</DIV>
	<DIV>&nbsp;**zip or gz archive&nbsp;</DIV>
	</TD>
</TR>
</TABLE>
</FORM>
<DIV id="waitDIV" style="display:none">
<BR><BR>
<FONT class="title2">Uploading file. Please wait...</FONT>
</DIV>
</CENTER>
</BODY>
</HTML>
|;
	exit;
}

### Import Archive file
	print qq
|<HTML>
<HEAD>
<TITLE>Import MaxQuant files</TITLE>
<LINK rel=stylesheet href="$promsPath{html}/promsStyle.css" type=text/css>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER><FONT class="title1">Importing MaxQuant data</FONT></CENTER>
<BR>
|;


my $evidenceFile='evidence.txt'; # Contains all peptide identified and peptide quantitation
my $peptideFile='peptides.txt'; # Contains all peptide identified and peptide quantitation
my $proteinGroupsFile='proteinGroups.txt'; # Contains all quantitation information at protein level (LFQ, iBAQ, Intensity,...)
my $msmsFile='msms.txt'; # Contains m/z and intensities information
my $summaryFile='summary.txt'; # Contains some redundant information found in mqpar.xml (varmod, label type, experimental design)
my $parametersFile='parameters.txt'; # Contains some redundant information found in mqpar.xml (fixmod, peptide parameters, software version)

#mkdir "$promsPath{tmp}/MaxQuant" unless -e "$promsPath{tmp}/MaxQuant";
&promsMod::cleanDirectory("$promsPath{tmp}/MaxQuant",'1d'); # [1 day]  also created if does not exists yet
my $tmpFilesDir="$promsPath{tmp}/MaxQuant/".(strftime("%Y%m%d%H%M%S",localtime)).'_'.$userID;
rmtree $tmpFilesDir if -e $tmpFilesDir; # clean directory
mkdir $tmpFilesDir;

my $importPepQuantif=1;
my $importProtQuantif=param('protQuantif');

my ($mqparFile);
if (param('mqpar')) {
		$mqparFile=(split(/[\\\/]/,param('mqpar')))[-1]; # XML-formatted parameters & design (replaces parameters.txt & summary.txt)
move(tmpFileName(upload('mqpar')),"$tmpFilesDir/$mqparFile");
}

my $newFile;
my $numFiles=0;
if (param('uplArch')) {
	#my $uplFile=(split(/[\\\/]/,param('uplArch')))[-1]; # Windows or Linux
	my (undef,$path,$uplFile)=splitpath(param('uplArch'));

	$newFile="$tmpFilesDir/$uplFile";
#my $tmpfile = tmpFileName(upload("uplArch"));
	move(tmpFileName(upload("uplArch")),$newFile); # name of temp file being uploaded
	$numFiles=1;
}
else { # from shared directory
	foreach my $sharedFile (param('sharedDirFiles')) {
		my (undef,$path,$fileName)=splitpath($sharedFile);
		$newFile="$tmpFilesDir/$fileName";
		move("$promsPath{shared}/$sharedFile",$newFile);
		if ($fileName =~/\.(gz|zip)\Z/) { # keep only archive if extra files
			$numFiles=1;
			last;
		}
		else {$numFiles++;}
	}
}

###>Inflating file
if ($numFiles==1) { # Must be an archive
	if ($newFile =~/\.(gz|zip)\Z/) {
		print "<BR><FONT class='title3'>Extracting files...";
		my $type=$1;
		if ($type eq 'gz') {
			if ($newFile =~ /\.tar\.gz\Z/) {
				system "tar -zxf $newFile -C $tmpFilesDir";
			}
			elsif ($newFile =~ /\.gz\Z/) {
				system "gunzip -c $newFile > $tmpFilesDir";
			}
		}
		elsif ($type eq 'zip') {
			#system "unzip -q -d $tmpFilesDir $newFile";
			print "</FONT><BR>\n";
			&promsMod::unzipArchive($newFile,$tmpFilesDir,{mode=>'verbose',txtBefore=>'<B>&nbsp;&nbsp;-',txtAfter=>"</B><BR>\n"});
			print "<FONT class='title3'>";
		}
		print " Done.</FONT><BR>\n";

		###>Deleting zip file
		unlink $newFile;
	}
	else {
		print "<BR><FONT class='title3' color=\"#DD0000\">ERROR: The archive type is not recognized (Use zip or gz only).</FONT><BR>\n";
		rmtree $tmpFilesDir;
		exit;
	}
}

###>Moving all txt files to top work directory (Only 1st-level sub-directories are checked!) <- in case extracted from archive with extra hierarchy
opendir (DIR,$tmpFilesDir);
while (my $item = readdir (DIR)) {
	next if $item=~/^\.+$/; #  skip '.' & '..' directories
	if (-d "$tmpFilesDir/$item") { # directory
		foreach my $txtFile (glob("$tmpFilesDir/$item/*.txt")) {
			my (undef, $path, $fileName) = splitpath($txtFile);
			move($txtFile,"$tmpFilesDir/$fileName");
		}
		rmtree "$tmpFilesDir/$item";
	}
}
close DIR;

my @requiredFiles=($evidenceFile,$peptideFile);
push @requiredFiles,$mqparFile if $mqparFile;
push @requiredFiles,$proteinGroupsFile if $importProtQuantif;
push @requiredFiles,($summaryFile,$parametersFile) if !$mqparFile;
foreach my $file (@requiredFiles) {
	if (!-e "$tmpFilesDir/$file") {
		$dbh->disconnect;
		print "<BR><FONT class='title3' color=\"#DD0000\">ERROR: File $file is missing!</FONT><BR>\n</BODY>\n</HTML>\n";
		exit;
	}
}

###>Link to a myProMS database in order to be consistent
my (@databankIDs,%idParseRules); # allows to used sub identifier instead of whole one (eg search with UNIPROT_ALL imported as UNIPROT_ID/ACC)
my $sthPR=$dbh->prepare("SELECT PARSE_RULES,IS_CRAP FROM DATABANK D,DATABANK_TYPE DT WHERE D.ID_DBTYPE=DT.ID_DBTYPE AND ID_DATABANK=?");
foreach my $d (1..$MAX_DB) {
	next unless param("databank$d");
	my $dbID=&promsMod::cleanNumericalParameters(param("databank$d"));
	push @databankIDs,$dbID;
	$sthPR->execute($dbID);
	my ($parseRules,$isCrap)=$sthPR->fetchrow_array;
	my @rules=split(',:,',$parseRules);
	($idParseRules{$dbID})=($rules[0]=~/ID=(.+)/);
}
$sthPR->finish;

my $contaminantDB=0;
my $excludeCON=param('excludeCON') || 0;
if (!$excludeCON && param('databankCON')) {
	$contaminantDB=&promsMod::cleanNumericalParameters(param('databankCON'));
	push @databankIDs,$contaminantDB;
}

###>Rule to compute match groups
my $matchGroupType=param('mgType') || 'myProMS';

###>Parameter group
my $paramGrIdx=param('paramGr') || 0;

#my @summaryColumns=('Raw file','Experiment','Enzyme','Enzyme mode','Variable modifications','Max. missed cleavages','Labels0','Labels1','Labels2');
my $dataFile=($mqparFile)?$mqparFile:'No mqpar.xml file';
my $sthInsS=$dbh->prepare("INSERT INTO SAMPLE (ID_SAMPLE,ID_EXPERIMENT,NAME,DISPLAY_POS,START_DATE,UPDATE_DATE,UPDATE_USER) VALUES (?,$experimentID,?,?,NOW(),NOW(),?)") || die $dbh->errstr;
my $sthInsA=$dbh->prepare("INSERT INTO ANALYSIS (ID_SAMPLE,NAME,START_DATE,DATA_FILE,VALID_STATUS,VALID_USER,LOWER_SCORES,FIRST_VALID_DATE,VALID_DATE,LABELING,WIFF_FILE,DECOY,DISPLAY_POS,FILE_FORMAT,MS_TYPE,MAX_RANK,MIN_SCORE,UPDATE_DATE,UPDATE_USER) VALUES (?,?,NOW(),'$dataFile',2,?,0,NOW(),NOW(),?,?,?,?,'MAXQUANT.DIR','MIS',1,0,NOW(),?)");
my $sthInsAM=$dbh->prepare("INSERT INTO ANALYSIS_MODIFICATION (ID_ANALYSIS,ID_MODIFICATION,MODIF_TYPE,SPECIFICITY) VALUES (?,?,?,?)");
my $sthLabMod=$dbh->prepare("UPDATE MODIFICATION SET IS_LABEL=1 WHERE ID_MODIFICATION=?");
my $sthInsAD=$dbh->prepare("INSERT INTO ANALYSIS_DATABANK (ID_ANALYSIS,ID_DATABANK,DB_RANK) VALUES (?,?,?)");
my $sthInsO=$dbh->prepare("INSERT INTO OBSERVATION (ID_ANALYSIS,TARGET_POS) VALUES (?,?)");
my $sthInsOM=$dbh->prepare("INSERT INTO OBS_MODIFICATION (ID_OBSERVATION,ID_MODIFICATION) VALUES(?,?)");
my ($maxSampleID) = $dbh->selectrow_array("SELECT MAX(ID_SAMPLE) FROM SAMPLE");
my ($displayPosSamp) = $dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM SAMPLE WHERE ID_EXPERIMENT=$experimentID");


#my ($quantifMethodID,$areaParamID);
my $sthInsQ=$dbh->prepare("INSERT INTO QUANTIFICATION (ID_QUANTIFICATION_METHOD,ID_DESIGN,NAME,FOCUS,STATUS,QUANTIF_ANNOT,UPDATE_USER,UPDATE_DATE) VALUES (?,?,?,?,1,?,'$userID',NOW())");
my $sthInsAQ=$dbh->prepare("INSERT INTO ANA_QUANTIFICATION (ID_QUANTIFICATION,ID_ANALYSIS) VALUES (?,?)");
my (%mqRepIon) = &getModIDforReporterIon;


my (%anaInfo,%anaNames,%vmodsInUnimod,%modifName2ID,%designInfo);
my (%anaVmods,%labels);
my (%sampIDs,%ana2Observation,%samplePos);


####> Global parameters & design used for MaxQuant search <####
my (@fixedMods,@varMods,%isobarModifInfo);
my (%labelInfos,%ana2pepQuantification);
my ($labeling,$labelingPlex);
my ($fdrPc,$pepUsed,$peptideFilterStrg,$versionStrg,$xmlParams);
my (@designRawFiles,@designExperiments,@designFractions,@designLabels);

if ($mqparFile) {
###>mqpar.xml<###
print "<BR><FONT class='title3'>Importing experimental design from mqpar file:<BR>\n";
my $xml = new XML::Simple();
		$xmlParams = $xml->XMLin("$tmpFilesDir/$mqparFile",ForceArray=>['parameterGroup','string','short','int'],SuppressEmpty=>undef);
		$fdrPc=$xmlParams->{peptideFdr} * 100;
#>--Fixed modifs--<#
if ($xmlParams->{fixedModifications}{string}) { # 'string' attribute missing if no modif at all!
	foreach my $modifStrg (@{$xmlParams->{fixedModifications}{string}}) {
		next unless $modifStrg;
		my ($varModName,$specificity)=&promsMod::convertVarModString($modifStrg);
		my $modID=&promsMod::getModificationIDfromString($dbh,$varModName,$specificity,\%vmodsInUnimod); # %vmodsInUnimod: prevents re-parsing of unimods_table file
		push @fixedMods,[$modID,$specificity];
		$modifName2ID{$varModName}=$modID;
	}
}
#>--Variable modifs--<#
if ($xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{variableModifications}{string}) {
	foreach my $modifStrg (@{$xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{variableModifications}{string}}) {
		next unless $modifStrg;
		my ($varModName,$specificity)=&promsMod::convertVarModString($modifStrg);
		my $modID=($modifName2ID{$varModName})? $modifName2ID{$varModName} : &promsMod::getModificationIDfromString($dbh,$varModName,$specificity,\%vmodsInUnimod);
		push @varMods,[$modID,$specificity,$modifStrg];
		$modifName2ID{$varModName}=$modID;
	}
}
#>--Peptides used for quantif--<#
		$pepUsed=($xmlParams->{quantMode}==1)? 'razor' : ($xmlParams->{quantMode}==2)? 'unique' : 'all'; # 0: all, 1: unique+razor, 2: unique
		$peptideFilterStrg="$pepUsed;1;"; # peptide used;missedCut;
if ($xmlParams->{restrictMods}{string}[0]) { # PTM used
	$peptideFilterStrg.='-' if $xmlParams->{useCounterparts} eq 'false'; # exclude unmodified matching peptides
	#>Check if all varMods are listed
	my @matchedVmods;
	foreach my $modifStrg (@{$xmlParams->{restrictMods}{string}}) {
		foreach my $refVarMod (@varMods) {
			if ($modifStrg eq $refVarMod->[2]) {
				push @matchedVmods,$modifStrg;
				last;
			}
		}
	}
	if (scalar @matchedVmods == scalar @varMods) {$peptideFilterStrg.='1';} # all allowed
	else {
		$peptideFilterStrg.='2:';
		my $firstPTM=1;
		foreach my $modifStrg (@matchedVmods) {
			if ($firstPTM) {$firstPTM=0;}
			else {$peptideFilterStrg.=',';}
			my ($varModName,$specificity)=&promsMod::convertVarModString($modifStrg);
			$peptideFilterStrg.='#'.$modifName2ID{$varModName};
		}
	}
}
else {$peptideFilterStrg.='0';} # no PTM allowed
$peptideFilterStrg.=';all;all'; # charge;source used
		@designRawFiles=@{$xmlParams->{filePaths}{string}};
		@designExperiments=@{$xmlParams->{experiments}{string}};
		@designFractions=@{$xmlParams->{fractions}{short}};
		@designLabels=@{$xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{labelMods}{string}};

}
else{
		print "<BR><FONT class='title3'>Importing experimental design from summary.txt and parameters.txt files:<BR>\n";

		open (SUMMARY,"$tmpFilesDir/$summaryFile") || die "Unable to open $tmpFilesDir/$summaryFile";
		my %summaryColNum;
		while (my $line=<SUMMARY>) {
				chomp($line);
				my @parameters=split(/\t/,$line);
				if ($parameters[0] eq 'Raw file') {
						my $ncol=0;
						foreach my $colName (@parameters) {
								$summaryColNum{$colName}=$ncol;
								$ncol++;
						}
				}
				elsif ($parameters[0] eq 'Total' ||  defined($designInfo{'Experiment'}{$parameters[0]}) ){
						# Experiment information is not kept so skip to the end !
						last;
				}
				else{
						if ($parameters[$summaryColNum{'LC-MS run type'}] =~ /Reporter ion MS/) {
								#code
								print "<BR><FONT class='title3' color=\"#DD0000\">WARNING: for Reporter ion MS2 or MS3, you need to provide the mqpar.xml file.</FONT><BR>\n";
								exit;
						}
						###> Add sample according to experiment rawname
						push @designRawFiles,$parameters[$summaryColNum{'Raw file'}];
						push @designExperiments,$parameters[$summaryColNum{'Experiment'}] if $summaryColNum{'Experiment'};
						push @designFractions,$parameters[$summaryColNum{'Fraction'}] if $summaryColNum{'Fraction'};
						if ($#designLabels <0) {
								if ( $summaryColNum{'Labels0'} ){
										$parameters[$summaryColNum{'Labels0'}]='' unless $parameters[$summaryColNum{'Labels0'}];
								}
								push @designLabels,$parameters[$summaryColNum{'Labels0'}] if $summaryColNum{'Labels0'};
								push @designLabels,$parameters[$summaryColNum{'Labels1'}] if $summaryColNum{'Labels1'};
								push @designLabels,$parameters[$summaryColNum{'Labels2'}] if $summaryColNum{'Labels2'};
								if ($summaryColNum{'Variable modifications'} && $parameters[$summaryColNum{'Variable modifications'}]) {
										foreach my $modifStrg (split(/;/,$parameters[$summaryColNum{'Variable modifications'}])) {
												next unless $modifStrg;
												my ($varModName,$specificity)=&promsMod::convertVarModString($modifStrg);
												my $modID=($modifName2ID{$varModName})? $modifName2ID{$varModName} : &promsMod::getModificationIDfromString($dbh,$varModName,$specificity,\%vmodsInUnimod);
												push @varMods,[$modID,$specificity,$modifStrg];
												$modifName2ID{$varModName}=$modID;
										}
								}
						}
						my ($sampleName)=($summaryColNum{'Experiment'}) ? $parameters[$summaryColNum{'Experiment'}] : "No_experiments_defined_in_maxquant";
						%{$designInfo{'Experiment'}{$sampleName}}=();
				}
		}
		close SUMMARY;
		
		my $addPepString='';
		open (PARAMFILE,"$tmpFilesDir/$parametersFile") || die "Unable to open $tmpFilesDir/$parametersFile";
		while (my $line=<PARAMFILE>) {
				chomp($line);
				$line=~s/\s+$//; # trailing spaces if any
				my @parameters=split(/\t/,$line);
				if ($parameters[0] eq 'Version') {
						##>Maxquant XIC software & version
						$parameters[1]=~s/\s+$//; # trailing spaces if any
						$versionStrg=";$parameters[1]";
				}
				elsif($parameters[0] eq 'PSM FDR') {
						$fdrPc=$parameters[1] * 100;
				}
				elsif($parameters[0] eq 'Fixed modifications' && $parameters[1] =~ /\w+/ ) {
						my $fixStg=$parameters[1];
						$fixStg=~ s/\s+$//g; # Remove last whitespace: 'Carbamidomethyl (C) ' becomes 'Carbamidomethyl (C)'
						foreach my $modifStrg (split(/;/,$fixStg)){
								next unless $modifStrg;
								my ($varModName,$specificity)=&promsMod::convertVarModString($modifStrg);
								my $modID=&promsMod::getModificationIDfromString($dbh,$varModName,$specificity,\%vmodsInUnimod); # %vmodsInUnimod: prevents re-parsing of unimods_table file
								push @fixedMods,[$modID,$specificity];
								$modifName2ID{$varModName}=$modID;
						}
				}
				elsif($parameters[0] eq 'Modifications included in protein quantification') {
						$peptideFilterStrg=";1;"; # peptide used;missedCut;
						#>Check if all varMods are listed
						my @matchedVmods;
						foreach my $modifStrg (split(/;/,$parameters[1])) {
								foreach my $refVarMod (@varMods) {
										if ($modifStrg eq $refVarMod->[2]) {
												push @matchedVmods,$modifStrg;
												last;
										}
								}
								if (scalar @matchedVmods == scalar @varMods) {$addPepString.='1';} # all allowed
								else {
										$addPepString.='2:';
										my $firstPTM=1;
										foreach my $modifStrg (@matchedVmods) {
												if ($firstPTM) {$firstPTM=0;}
												else {$addPepString.=',';}
												my ($varModName,$specificity)=&promsMod::convertVarModString($modifStrg);
												$addPepString.='#'.$modifName2ID{$varModName};
										}
								}
						}

				}
				elsif($parameters[0] eq 'Peptides used for protein quantification') {
						$pepUsed=($parameters[1]=~ /Razor/)? 'razor' : ($parameters[1]=~ /Unique/)? 'unique' : 'all';
				}
				elsif($parameters[0] eq 'Discard unmodified counterpart peptides') {
						$peptideFilterStrg=$peptideFilterStrg."-" if $parameters[1] =~ /True/;
						$peptideFilterStrg=$pepUsed.$peptideFilterStrg.$addPepString;
				}
		}
		close PARAMFILE;
		if (!$peptideFilterStrg){ $peptideFilterStrg="$pepUsed;1;0"; }
		$peptideFilterStrg.=';all;all'; # charge;source used
		print " Done.</FONT><BR>\n";


}

##>Labeling<##
my (%labelModifSpecificity,%isotopeLabelInfo,%labelName,@labelList,%label2targetPos);
my ($pepQuantifName,$labelingName,$pepQuantifAnnot,$quantifMethodID,$areaParamID);


#>------SILAC-----<# !!! Actually any Isotopic labeling !!!
#if (($xmlParams && $xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{multiplicity} > 1)||$#designLabels>1) { # $xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{labelMods}{string}[1]
if ($#designLabels>0) {
	$labelingPlex=($xmlParams)?$xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{multiplicity}:1+$#designLabels; #scalar @{$xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{labelMods}{string}};
	$labeling='SILAC';
	$labelingName='SILAC '.$labelingPlex.'plex';
	$pepQuantifName='SILAC made by MaxQuant';
	print "&nbsp;-$labelingName detected<BR>\n";

	$pepQuantifAnnot='LABEL=SILAC';
	%labelName=('L'=>'Light','M'=>'Medium','H'=>'Heavy');
	my (@colNames,@labNames);
	if ($labelingPlex == 2) { # possible ambiguity L/H or L/M or M/H? => read evidence.txt
		#@colNames=('L','H');
		#@labNames=('Light','Heavy');
		my $headEvStrg=`head -1 $tmpFilesDir/$evidenceFile`;
		chomp($headEvStrg);
		my @colHeaders=split(/\t/,$headEvStrg);
		my $targetPos=0;
		foreach my $colLabel ('L','M','H') {
			foreach my $header (@colHeaders) {
				if ($header eq "Intensity $colLabel") {
					push @colNames,$colLabel;
					push @labNames,$labelName{$colLabel};
					$label2targetPos{$colLabel}=++$targetPos;
					last;
				}
			}
		}
	}
	else {
		@colNames=('L','M','H');
		@labNames=('Light','Medium','Heavy');
		%label2targetPos=('L'=>1,'M'=>2,'H'=>3);
	}
	@labelList=@colNames;

	foreach my $labelIdx (0..$#designLabels) {
		#my $label=($labelIdx==0)? 'NA' :
		my @modifications=($labelIdx==0)? ('None') : split(/;/,$designLabels[$labelIdx]);
		my $labelCol=shift @colNames;
		my $labName=shift @labNames;
		my $targetPos=$labelIdx+1;
		$pepQuantifAnnot.="::$targetPos;$labName;";
		my @residueMods;
		foreach my $repIon (@modifications) {
			if ($repIon eq 'None') {
				push @residueMods,'None#No label##';
			}
			else {
				my $specificity=($repIon =~ /Lys/)?'K':($repIon =~ /Arg/)?'R':($repIon =~ /Nter/) ? '=' : ($repIon =~ /Ile7/)? 'I': ($repIon =~ /Leu7/) ? 'L' : ($repIon =~ /ICAT/) ? 'C' : ($repIon =~ /18O/) ? '*' : '';
				my ($unimodID)=$mqRepIon{$repIon}{'UNIMOD_ACC'};
				my ($modID,$psiName,$interimName)=$dbh->selectrow_array("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME FROM MODIFICATION WHERE UNIMOD_ACC=$mqRepIon{$repIon}{UNIMOD_ACC}");
				unless ($modID) { # not yet in DB
					$modID=&promsMod::getModificationIDfromString($dbh,$mqRepIon{$repIon}{'INTERIM_NAME'},$specificity,\%vmodsInUnimod);
					$sthLabMod->execute($modID); # set as label modif (just to be safe)
					($psiName,$interimName)=$dbh->selectrow_array("SELECT PSI_MS_NAME,INTERIM_NAME FROM MODIFICATION WHERE ID_MODIFICATION=$modID");
				}
				$psiName=$interimName unless $psiName;
				$psiName='' unless $psiName;
				push @residueMods,"$repIon#$psiName#$specificity#$psiName#$modID";
				$isotopeLabelInfo{$labelCol}{$modID}{$specificity}=1;
				$labelModifSpecificity{$modID}{$specificity}=1;
			}
		}
		$pepQuantifAnnot.=join('@',@residueMods);
	}
	($quantifMethodID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='SILAC'"); # SILAC
	($areaParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE='ISO_AREA' AND ID_QUANTIFICATION_METHOD=$quantifMethodID");
}
#>------iTRAQ or TMT-----<#
elsif ($xmlParams && $xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{isobaricLabels}) { # iTRAQ, TMT
	my %isobarTag;
	my %aaThree2OneLetter=&promsMod::getAAThree2OneLetter;
	foreach my $modifTgt (@{$xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx]{isobaricLabels}{string}}) {
		#my ($modif,$res,$tag)=($modifTgt=~/(.+)-([^-\d]+)(\d+)$/);
		###> Examples:
		###> TMT: TMT10plex-Nter126C,  TMT10plex-Lys127C, etc.
		###> iTRAQ: iTRAQ8plex-Nter113, iTRAQ8plex-Lys114, etc.
		my ($modif,$res,$tag)=($modifTgt=~/(.+)-([^-\d]+)(\d+\w?)$/);
		$isobarModifInfo{NAME}=$modif;
		$res=($res eq 'Nter')? '=' : ($res eq 'Cter')? '*' : $aaThree2OneLetter{$res};
		$isobarModifInfo{RES}{$res}=1;
		$isobarTag{$tag}=1;
	}
	if ($isobarModifInfo{NAME}=~/iTRAQ/i) {
		$labeling='ITRAQ';
		$pepQuantifName='iTRAQ made by MaxQuant';
	}
	elsif ($isobarModifInfo{NAME}=~/TMT/i) {
		$labeling='TMT';
		$pepQuantifName='TMT made by MaxQuant';
	}
	else { # Warning: Should not happen!
		$labeling='ISOBAR';
		$pepQuantifName='Isobaric made by MaxQuant';
	}
	($labelingPlex)=($isobarModifInfo{NAME}=~/(\d+)plex/i);
	$labelingName=$isobarModifInfo{NAME};
	print "&nbsp;-$labelingName labeling detected<BR>\n";

	$pepQuantifAnnot='LABEL='.$labeling;
	my $tagPos=0;
	@{$isobarModifInfo{TAGS}}=();
	foreach my $tag (sort{$a cmp $b} keys %isobarTag) {
		push @{$isobarModifInfo{TAGS}},$tag; # index!!!
		$tagPos++;
		$pepQuantifAnnot.="::$tagPos;$tag;";
	}
	my $specificity=join(',',sort keys %{$isobarModifInfo{RES}});
	$isobarModifInfo{ID}=&promsMod::getModificationIDfromString($dbh,$isobarModifInfo{NAME},$specificity,\%vmodsInUnimod);
	$sthLabMod->execute($isobarModifInfo{ID}); # set as label modif (just to be safe)
	$labelModifSpecificity{$isobarModifInfo{ID}}=$isobarModifInfo{RES};
	($quantifMethodID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='$labeling'");
	($areaParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE='REP_AREA' AND ID_QUANTIFICATION_METHOD=$quantifMethodID");
	($isobarModifInfo{MONO_MASS})=$dbh->selectrow_array("SELECT MONO_MASS FROM MODIFICATION WHERE ID_MODIFICATION=$isobarModifInfo{ID}");

	foreach my $repIdx (0..($labelingPlex-1)) {
		$label2targetPos{$repIdx}=$repIdx+1;
		$labelName{$repIdx}=$isobarModifInfo{TAGS}[$repIdx];
		push @labelList,$repIdx;
	}
}
#>------label-free-----<#
else {
	$labeling='FREE';
	$labelingPlex=0;
	$pepQuantifName='Label-free made by MaxQuant';
	$labelingName=undef;
	print "&nbsp;-No labeling<BR>\n";

	$pepQuantifAnnot='LABEL=FREE';
	($quantifMethodID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='XIC'"); # label-free
	($areaParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE='XIC_AREA' AND ID_QUANTIFICATION_METHOD=$quantifMethodID");
}
$sthLabMod->finish;

##>Maxquant XIC software & version
$versionStrg=';'.$xmlParams->{maxQuantVersion} if (!$versionStrg && $xmlParams && $xmlParams->{maxQuantVersion});
$pepQuantifAnnot.='::SOFTWARE=MQ'.$versionStrg;

##>Insert peptide QUANTIFICATION
$sthInsQ->execute($quantifMethodID,undef,$pepQuantifName,'peptide',$pepQuantifAnnot);
my $peptideQuantifID=$dbh->last_insert_id(undef,undef,'QUANTIFICATION','ID_QUANTIFICATION');

##>MaxQuant analysis design<##
my %displayPosAna;
#my $firstSampID;
my $expCondPos=0;
foreach my $anaIdx (0..$#designRawFiles) {
	if ($xmlParams && $xmlParams->{paramGroupIndices}{int}[$anaIdx] != $paramGrIdx) {
		print "<BR><FONT class='title3' color=\"#DD0000\">ERROR: Import of more than 1 group is not supported (see paramGroupIndices in mqpar.xml).</FONT><BR>\n";
		rmtree $tmpFilesDir;
		exit;
	}

	my ($sampleName,$realSampleName)=($#designExperiments>0 && $designExperiments[$anaIdx])? ($designExperiments[$anaIdx],' '.$designExperiments[$anaIdx]) : ("No experiment",'');

	my $fraction;
	$fraction=($#designFractions>0 && $designFractions[$anaIdx])? $designFractions[$anaIdx] : 0;
	$fraction=0 if $fraction==32767; # 32767 if no fraction

	unless ($sampIDs{$sampleName}) {
		$maxSampleID++;
		$displayPosSamp++;
		$sthInsS->execute($maxSampleID,$sampleName,$displayPosSamp,$userID) || die $dbh->errstr;
		$sampIDs{$sampleName}=$maxSampleID;
		#$firstSampID=$maxSampleID unless $firstSampID;
		$displayPosAna{$sampleName}=0; # sample-specific needed in case mixture of fractionnated & non-fractionnated samples (sample files may be mixted)
	}
	$displayPosAna{$sampleName}++;
	my $usedDisplayPos=$fraction || $displayPosAna{$sampleName}; # WARNING: There should not be mixture of fractions and no fraction for the same sample!!!

	my ($rawFile)=($designRawFiles[$anaIdx]=~/([^\\\/]+)$/);
	(my $anaName=$rawFile)=~s/\.[^\.]+$//;

	$sthInsA->execute($sampIDs{$sampleName},$anaName,$userID,$labelingName,$rawFile,"INT:SEARCH,FDR=$fdrPc:precomputed",$usedDisplayPos,$userID) || die $dbh->errstr;
	my $analysisID=$dbh->last_insert_id(undef,undef,'ANALYSIS','ID_ANALYSIS');

	$anaNames{$analysisID}=$anaName;
	$anaInfo{$anaName}{'ID_ANALYSIS'}=$analysisID;

	##>data directories
	mkdir "$promsPath{peptide}/proj_$projectID" unless -e "$promsPath{peptide}/proj_$projectID";
	mkdir "$promsPath{peptide}/proj_$projectID/ana_$analysisID";
	copy("$tmpFilesDir/$mqparFile","$promsPath{peptide}/proj_$projectID/ana_$analysisID/$mqparFile") if $mqparFile;

	#>ANALYSIS_DATABANK
	my $dbRank=1;
	foreach my $dbID (@databankIDs) {
		$sthInsAD->execute($analysisID,$dbID,$dbRank);
		$dbRank++;
	}

	#>ANALYSIS_MODIFICATION
	#-Fixed
	foreach my $refMod (@fixedMods) { # [modID,specif]
		$sthInsAM->execute($analysisID,$refMod->[0],'F',$refMod->[1]);
	}
	#-Variable
	foreach my $refMod (@varMods) { # [vModID,Specif,vModStrg]
		$sthInsAM->execute($analysisID,$refMod->[0],'V',$refMod->[1]);
		$anaVmods{$analysisID}{$refMod->[2]}=$refMod->[0];
	}
	#-Labeling
	foreach my $modID (keys %labelModifSpecificity) {
		$sthInsAM->execute($analysisID,$modID,'V',join(',',sort keys %{$labelModifSpecificity{$modID}}));
	}

	#>ANALYSIS_QUANTIFICATION
	$sthInsAQ->execute($peptideQuantifID,$analysisID);
	$ana2pepQuantification{$analysisID}=$peptideQuantifID;

#next unless $importProtQuantif;

	#>------SILAC design------<#
	if ($labeling eq 'SILAC') {
		my @colNames=('L','H');
		#$labels{$analysisID}{'Labels0'}='NA'; # add label-free channel
		#$labels{$analysisID}{'Labels1'}=$parameters[$summaryColNum{'Labels1'}];
		if ($labelingPlex==3) {
			@colNames=('L','M','H');
			#$labels{$analysisID}{'Labels2'}=$parameters[$summaryColNum{'Labels2'}];
		}
		foreach my $colName (@colNames) {
			$sthInsO->execute($analysisID,-1); # $targetPosObs
			my $obsID=$dbh->last_insert_id(undef,undef,'OBSERVATION','ID_OBSERVATION');
			$ana2Observation{$analysisID}{$colName}=$obsID;
			if ($isotopeLabelInfo{$colName}) {
				foreach my $modID (keys %{$isotopeLabelInfo{$colName}}) {
					$sthInsOM->execute($obsID,$modID);
				}
			}
			my $expCondName="$colName$realSampleName";
			$designInfo{'Expcondition'}{$expCondName}{'ANALYSIS'}{$analysisID}=1;
			$designInfo{'ExpconditionOriginalName'}{$sampleName}{$expCondName}=1;
			$designInfo{'Fraction'}{$expCondName}{$fraction}{$analysisID}=1 if $fraction;
			$designInfo{'Expcondition2Label'}{$expCondName}=$colName; # Keep a memory of which label was used in order to link Observation and modification in the end
			$designInfo{'Expcondition'}{$expCondName}{'POS'}=++$expCondPos unless $designInfo{'Expcondition'}{$expCondName}{'POS'}; # only once/sample
		}
	}
	#>------Isobaric design------<#
	elsif ($isobarModifInfo{ID}) { # $labeling=~/ITRAQ|TMT/
		foreach my $repIdx (0..($labelingPlex-1)) {
			$sthInsO->execute($analysisID,$repIdx+1); # $targetPosObs
			my $obsID=$dbh->last_insert_id(undef,undef,'OBSERVATION','ID_OBSERVATION');
			$ana2Observation{$analysisID}{$repIdx}=$obsID;
			$sthInsOM->execute($obsID,$isobarModifInfo{ID});
			#my $expCondName="Reporter intensity $repIdx$realSampleName";
			my $expCondName="$repIdx$realSampleName";
			$designInfo{'Expcondition'}{$expCondName}{'ANALYSIS'}{$analysisID}=1;
			$designInfo{'ExpconditionOriginalName'}{$sampleName}{$expCondName}=1;
			$designInfo{'Fraction'}{$expCondName}{$fraction}{$analysisID}=1 if $fraction;
			$designInfo{'Expcondition2Label'}{$expCondName}=$repIdx; #"Reporter intensity $repIdx"; # Keep a memory of which label was used in order to link Observation and modification in the end
			$designInfo{'Expcondition'}{$expCondName}{'POS'}=++$expCondPos unless $designInfo{'Expcondition'}{$expCondName}{'POS'}; # only once/sample
		}
	}
	elsif ($labeling eq 'FREE') { # Label-free
		$sthInsO->execute($analysisID,0); # $targetPosObs
		$ana2Observation{$analysisID}{'NONE'}=$dbh->last_insert_id(undef, undef, 'OBSERVATION', 'ID_OBSERVATION');
		$designInfo{'Expcondition'}{$sampleName}{'ANALYSIS'}{$analysisID}=1;
		$designInfo{'ExpconditionOriginalName'}{$sampleName}{$sampleName}=1;
		$designInfo{'Fraction'}{$sampleName}{$fraction}{$analysisID}=1 if $fraction;
		$designInfo{'Expcondition'}{$sampleName}{'POS'}=++$expCondPos unless $designInfo{'Expcondition'}{$sampleName}{'POS'}; # only once/sample
#$designInfo{'Expcondition2Label'}{$ampleName}='NONE';
	}
	$designInfo{'Experiment'}{$sampleName}{$analysisID}=1;
}
$sthInsAM->finish;
$sthInsOM->finish;
$sthInsS->finish;
$sthInsA->finish;
$sthInsAD->finish;
$sthInsO->finish;

my $numSamples=scalar keys %sampIDs;
my $sampStrg=($numSamples==1)? 'Sample' : 'Samples';
my $numAnalyses=scalar keys %anaNames;
my $anaStrg=($numAnalyses==1)? 'Analysis' : 'Analyses';
print "&nbsp;-$numSamples $sampStrg and $numAnalyses $anaStrg created.<BR>\n";


###############################################
####>Check mqpar & txt files compatibility<####
###############################################
my $peptideHeadStrg=`head -1 $tmpFilesDir/$peptideFile`;
chomp($peptideHeadStrg);
my ($numExpPeptide,$numExpMatched)=(0,0);
my $numExpParam=($numSamples==1 && (keys %sampIDs)[0] eq 'No experiment')? 0 : $numSamples;
foreach my $colHeader (split(/\t/,$peptideHeadStrg)) {
	$colHeader=~s/^\s+//; $colHeader=~s/\s+$//; # remove starting & trailing spaces if any
	if ($colHeader=~/^Experiment (.+)/) {
		my $expName=$1;
		$numExpMatched++ if $sampIDs{$expName};
		$numExpPeptide++;
	}
}
if ($numExpMatched < $numExpParam || $numExpPeptide > $numExpParam) {
	$dbh->rollback;
	$dbh->disconnect;
	rmtree $tmpFilesDir;
	print "<BR><FONT class='title3' color=\"#DD0000\">ERROR: Parameter file '$mqparFile' does not match data files!</FONT><BR>\n</BODY>\n</HTML>\n";
	exit;
#print "<BR><FONT class='title3' color=\"#DD0000\">ERROR: Parameter file '$mqparFile' does not match data files!</FONT><BR>\n";
}

########################
####>myproMS Design<####
########################
#print "<BR><FONT class='title3'>Generating protein quantification design(s)...";
my $importDate=strftime("%Y/%m/%d %H:%M:%S",localtime);
my ($designID)=$dbh->selectrow_array("SELECT MAX(ID_DESIGN) FROM DESIGN");
$designID++;
$dbh->do("INSERT INTO DESIGN (ID_DESIGN,ID_EXPERIMENT,NAME,DES,UPDATE_DATE,UPDATE_USER) VALUES ($designID,$experimentID,'MaxQuant [$importDate]','Automatically generated during MaxQuant data import',NOW(),'$userID')");
my ($maxExpCondID)=$dbh->selectrow_array("SELECT MAX(ID_EXPCONDITION) FROM EXPCONDITION");
my $sthInsEC=$dbh->prepare("INSERT INTO EXPCONDITION (ID_EXPCONDITION,ID_DESIGN,NAME) VALUES (?,?,?)");
my $sthInsOE=$dbh->prepare("INSERT INTO OBS_EXPCONDITION (ID_EXPCONDITION,ID_OBSERVATION,FRACTION_GROUP) VALUES(?,?,?)");

#my (%label2targetPos,%labelName,@labelList);
#if ($labeling eq 'SILAC') {
#	if ($labelingPlex==2) {
#		%label2targetPos=('L'=>1,'H'=>2);
#		@labelList=('L','H');
#	}
#	else {
#		%label2targetPos=('L'=>1,'M'=>2,'H'=>3);
#		@labelList=('L','M','H');
#	}
#	%labelName=('L'=>'Light','M'=>'Medium','H'=>'Heavy');
#
#}
#elsif ($isobarModifInfo{ID}) { # $labeling=~/ITRAQ|TMT/
#	foreach my $repIdx (0..($labelingPlex-1)) {
#		$label2targetPos{$repIdx}=$repIdx+1;
#		$labelName{$repIdx}=$isobarModifInfo{TAGS}[$repIdx];
#		push @labelList,$repIdx;
#	}
#}
##else {%label2targetPos=('NONE'=>0);}

my $numStates=0;
if ($labeling ne 'FREE') {
	#>EXPCONDITION
	foreach my $label (@labelList) {
		$maxExpCondID++;
		$sthInsEC->execute($maxExpCondID,$designID,$labelName{$label});
		$designInfo{'LABEL_EXPCOND'}{$label}{'ID'}=$maxExpCondID;
		$numStates++;
	}
}

my %labelReplicCount;
foreach my $sampleName (sort{$sampIDs{$a}<=>$sampIDs{$b}} keys %{$designInfo{'ExpconditionOriginalName'}}) {
	foreach my $expCondName (sort{$designInfo{'Expcondition'}{$a}{'POS'}<=>$designInfo{'Expcondition'}{$b}{'POS'}} keys %{$designInfo{'ExpconditionOriginalName'}{$sampleName}}) { # %{$designInfo{'Expcondition'}}
		my ($labType,$targetPos);
		if ($labeling eq 'FREE') { # Label-free -> each MQ experiment is a state (EXPCONDITION)
			$labType='NONE';
			$targetPos=0;
			$maxExpCondID++;
			$designInfo{'Expcondition'}{$expCondName}{'ID'}=$maxExpCondID; # needed for EXPCONDITION_QUANTIF
			$sthInsEC->execute($maxExpCondID,$designID,$expCondName);
			$numStates++;
		}
		else {
			$labType=$designInfo{'Expcondition2Label'}{$expCondName};
			$targetPos=$label2targetPos{$labType};
			$designInfo{'Expcondition'}{$expCondName}{'ID'}=$designInfo{'LABEL_EXPCOND'}{$labType}{'ID'}; # needed for EXPCONDITION_QUANTIF
		}
		#$targetPos=($expCondName=~/^([LMH]) /)? $label2targetPos{$1} : 0;
		#my ($labType)=($designInfo{'Expcondition2Label'}{$expCondName})? $designInfo{'Expcondition2Label'}{$expCondName} : 'NONE';
		my %fractions;
		if ($designInfo{'Fraction'} && $designInfo{'Fraction'}{$expCondName}) {
			my %replicCount;
			foreach my $fraction (sort{$a<=>$b} keys %{$designInfo{'Fraction'}{$expCondName}}){
				#my $techRep=0;
				foreach my $analysisID (sort{$a<=>$b} keys %{$designInfo{'Fraction'}{$expCondName}{$fraction}}){ # only 1 possible!? (PP 21/11/16)
					#$techRep++;
					$replicCount{$fraction}++;
					$labelReplicCount{$labType}{$fraction}++;
#next unless $ana2pepQuantification{$analysisID}; # No peptide quantif data for this analysis -> skip (PP 08/11/16)
					if ($labeling eq 'FREE') {
						$sthInsOE->execute($maxExpCondID,$ana2Observation{$analysisID}{$labType},$replicCount{$fraction});
					}
					else {
						$sthInsOE->execute($designInfo{'LABEL_EXPCOND'}{$labType}{'ID'},$ana2Observation{$analysisID}{$labType},$labelReplicCount{$labType}{$fraction}); # TODO: check $labelReplicCount{}{}
					}
					push @{$fractions{$replicCount{$fraction}}},"#$ana2Observation{$analysisID}{$labType}:#$ana2pepQuantification{$analysisID}:#$analysisID:$targetPos";
				}
			}
		}
		else {
			my $replicNum=0;
			foreach my $analysisID (keys %{$designInfo{'Expcondition'}{$expCondName}{'ANALYSIS'}}) {
				$replicNum++;
				$labelReplicCount{$labType}++;
#next unless $ana2pepQuantification{$analysisID}; # No peptide quantif data for this analysis -> skip (PP 08/11/16)
				if ($labeling eq 'FREE') {
					$sthInsOE->execute($maxExpCondID,$ana2Observation{$analysisID}{$labType},undef);
				}
				else {
					$sthInsOE->execute($designInfo{'LABEL_EXPCOND'}{$labType}{'ID'},$ana2Observation{$analysisID}{$labType},$labelReplicCount{$labType}); # TODO: check $labelReplicCount{}
				}
				push @{$fractions{$replicNum}},"#$ana2Observation{$analysisID}{$labType}:#$ana2pepQuantification{$analysisID}:#$analysisID:$targetPos";
			}
		}
		my @replicates;
		foreach my $replic (sort{$a<=>$b} keys %fractions) {
			push @replicates,join('+',@{$fractions{$replic}});
		}
		#$designInfo{'ALL_STATES'}.=';' if $designInfo{'STATES'};
		#$designInfo{'ALL_STATES'}.=(scalar @replicates).','.join('.',@replicates).',#'.$maxExpCondID;
		#push @{$designInfo{'STATES'}{$sampleName}{$expCondName}},(scalar @replicates).','.join('.',@replicates).',#'.$maxExpCondID;
		push @{$designInfo{'STATES'}{$sampleName}},(scalar @replicates).','.join('.',@replicates).',#'.$designInfo{'Expcondition'}{$expCondName}{'ID'};
		push @{$designInfo{'ALL_STATES'}{$targetPos}},@replicates;

	}
}
$sthInsEC->finish;
$sthInsOE->finish;

my $stateStrg=($numStates==1)? 'State' : 'States';
print "&nbsp;-1 Design and $numStates $stateStrg generated.<BR>\n";


###################################
####>Protein Quantification(s)<####
###################################
my (%quantifParamIDs,%requiredParams,@protQuantifList,%recordedParams);
my ($MQMethID,$ratioQuantMethID);
my $numProtQuantif=0;
my  $correctedStrg=($labeling eq 'TMT')? 'not corrected' : 'corrected'; # for isobaric only. *******WARNING: Could be controlled by a parameter rather than by isobaric label type!!!********
if ($importProtQuantif) {

	####>Fetching list of quantification parameters<####
	($MQMethID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='MQ'");
	my $sthQP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,P.CODE FROM QUANTIFICATION_PARAMETER P,QUANTIFICATION_METHOD M WHERE P.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND M.ID_QUANTIFICATION_METHOD=?");
	$sthQP->execute($MQMethID);
	while (my ($paramID,$code)=$sthQP->fetchrow_array) {
		$quantifParamIDs{$code}=$paramID;
	}

	if ($labeling eq 'SILAC') {
		($ratioQuantMethID)=$dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='PROT_RATIO_PEP'");
		$sthQP->execute($ratioQuantMethID);
		while (my ($paramID,$code)=$sthQP->fetchrow_array) {
			$code="RATIO:$code" if $code=~/PEPTIDES|RAZ_UNI_PEP|UNIQUE_PEP/; # to distinguish from same codes used in MQ method
			$quantifParamIDs{$code}=$paramID;
		}
	}
	$sthQP->finish;


	my %quantifParamCodes=(
		'MQ_INT'=>'Intensity','MQ_IBAQ'=>'iBAQ','MQ_LFQ'=>'LFQ intensity','MQ_SC'=>'MS/MS Count',
		'ISOBAR:MQ_INT'=>"Reporter intensity $correctedStrg",
		'RATIO_RAW'=>'Ratio #%#','RATIO'=>'Ratio #%# normalized','RATIO_VAR'=>'Ratio #%# variability [%]',
		'PEPTIDES'=>'Peptides','RAZ_UNI_PEP'=>'Razor + unique peptides','UNIQUE_PEP'=>'Unique peptides',
		'NUM_PEP_USED'=>'Ratio #%# count'
	);
	#my @noRatioParamCodes=('MQ_INT','MQ_IBAQ','MQ_LFQ','MQ_SC','PEPTIDES','RAZ_UNI_PEP','UNIQUE_PEP');
	#my @ratioParamCodes=('RATIO_RAW','RATIO_NORM','RATIO_VAR','PEPTIDES','RAZ_UNI_PEP','UNIQUE_PEP','NUM_PEP_USED');
	my @stateIntParamCodes=($isobarModifInfo{ID})? ('ISOBAR:MQ_INT') : ('MQ_INT','MQ_IBAQ','MQ_LFQ','MQ_SC');
	my @ratioParamCodes=('RATIO_RAW','RATIO','RATIO_VAR','NUM_PEP_USED');
	my @allStatesParamCodes=('PEPTIDES','RAZ_UNI_PEP','UNIQUE_PEP');

	my $sthInsEQ=$dbh->prepare("INSERT INTO EXPCONDITION_QUANTIF (ID_EXPCONDITION,ID_QUANTIFICATION) VALUES (?,?)");
	my $sthInsParQ=$dbh->prepare("INSERT IGNORE INTO PARENT_QUANTIFICATION (ID_PARENT_QUANTIFICATION,ID_QUANTIFICATION) VALUES (?,?)");

	####
	####>For each sample: Create a non-ratio quantification (Intensity,iBAQ,LFQ,SpecCount) containing all (non-)label channels<####
	####
	foreach my $sampleName (sort{$sampIDs{$a}<=>$sampIDs{$b}} keys %{$designInfo{'ExpconditionOriginalName'}}) { # eg. label '123','456'
		my $sampleNameStrg=($sampleName eq 'No experiment')? '' : " $sampleName";
		$numProtQuantif++;
		my $quantifAnnot="LABEL=$labeling";
		$quantifAnnot.='::SOFTWARE=MQ'.$versionStrg.'::PEPTIDES='.$peptideFilterStrg.'::RATIO_TYPE=None::STATES='.join(';',@{$designInfo{'STATES'}{$sampleName}});
		my $quantifName=($numSamples==1)? $sampleName : "$sampleName [Intensities]";
		$sthInsQ->execute($MQMethID,$designID,$quantifName,'protein',$quantifAnnot) || die $sthInsQ->errstr();
		my $quantifID=$dbh->last_insert_id(undef,undef,'QUANTIFICATION','ID_QUANTIFICATION');
		push @protQuantifList,$quantifID;
		$sthInsParQ->execute($peptideQuantifID,$quantifID) || die $sthInsParQ->errstr();
		my %linkedAna;
		my $targetPos=0;
		foreach my $expCondName (sort{$designInfo{'Expcondition'}{$a}{'POS'}<=>$designInfo{'Expcondition'}{$b}{'POS'}} keys %{$designInfo{'ExpconditionOriginalName'}{$sampleName}} ) { # label: ' L 123','H 123',... label-free: '$sampleName'
			$targetPos++;
			my $expCondNameStrg=($expCondName eq 'No experiment')? '' : " $expCondName";
			#$quantifAnnot='LABEL=FREE::SOFTWARE=MQ'.$versionStrg.'::RATIO_TYPE=None::STATE='.$designInfo{'Expcondition'}{$expCondName}{'STATE'};
			#$sthinsQ->execute($MQMethID,$designID,$expCondName,'protein',$quantifAnnot) || die $sthinsQ->errstr();
			#my $quantifID=$dbh->last_insert_id(undef,undef,'QUANTIFICATION','ID_QUANTIFICATION');
			$sthInsEQ->execute($designInfo{'Expcondition'}{$expCondName}{'ID'},$quantifID) || die $sthInsEQ->errstr();

			if ($designInfo{'Fraction'} && $designInfo{'Fraction'}{$expCondName}) {
				foreach my $fraction (keys %{$designInfo{'Fraction'}{$expCondName}}){
					foreach my $analysisID (keys %{$designInfo{'Fraction'}{$expCondName}{$fraction}}) {
#next unless $ana2pepQuantification{$analysisID}; # No peptide quantif data for this analysis -> skip (PP 08/11/16)
						next if $linkedAna{$analysisID}; # do only once
						$sthInsAQ->execute($quantifID,$analysisID) || die $sthInsAQ->errstr();
						#$sthInsParQ->execute($ana2pepQuantification{$analysisID},$quantifID) || die $sthInsParQ->errstr();
						$linkedAna{$analysisID}=1;
					}
				}
			}
			else {
				#foreach my $analysisID (keys %{$designInfo{'Experiment'}{$expCondName}}){
				foreach my $analysisID (keys %{$designInfo{'Expcondition'}{$expCondName}{'ANALYSIS'}}) {
#next unless $ana2pepQuantification{$analysisID}; # No peptide quantif data for this analysis -> skip (PP 08/11/16)
					next if $linkedAna{$analysisID}; # do only once
					$sthInsAQ->execute($quantifID,$analysisID) || die $sthInsAQ->errstr();
					#$sthInsParQ->execute($ana2pepQuantification{$analysisID},$quantifID) || die $sthInsParQ->errstr();
					$linkedAna{$analysisID}=1;
				}
			}

			###>Make parameters state-specific: channel-spec for label ('Intensity L 123'), sample-spec for label-free ('Intensity')
			foreach my $paramCode (@stateIntParamCodes) {
				#my $usedParamName=($isobarModifInfo{ID} && $paramCode eq 'MQ_INT')? "Reporter intensity $correctedStrg" : $quantifParamCodes{$paramCode};
				(my $trueParamCode=$paramCode)=~s/ISOBAR://; # eg. ISOBAR:CODE -> CODE
				push @{$recordedParams{"$quantifParamCodes{$paramCode}$expCondNameStrg"}},[$quantifID,$trueParamCode,$targetPos]; # push because same param used for different quantifs (in case of ratios)
			}
			#$requiredParams{$quantifID}{$targetPos}=($isobarModifInfo{ID})? "Reporter intensity $correctedStrg".$expCondNameStrg : 'Intensity'.$expCondNameStrg;
			$requiredParams{$quantifID}{$targetPos}=($isobarModifInfo{ID})? $quantifParamCodes{'MQ_INT'}.$sampleNameStrg : $quantifParamCodes{'MQ_INT'}.$expCondNameStrg;
		}

		##>Make peptide params sample-specific: These peptide params not state-specific (global) for sample-level quantif
		foreach my $paramCode (@allStatesParamCodes) {
			push @{$recordedParams{"$quantifParamCodes{$paramCode}$sampleNameStrg"}},[$quantifID,$paramCode,0]; # push because same param used for different quantifs (in case of ratios)
		}

	}

	####
	####>Add a set of quantitations to keep labeling information (Ratio value)<####
	####
	if ($labeling eq 'SILAC') { # Add Ratio (No ratios for iTRAQ!!!)
		#my @colNames=($labelingPlex == 2)? ('L','H') : ('L','M','H');

		my $baseQuantifAnnot='LABEL='.$labeling.'::SOFTWARE=MQ'.$versionStrg.'::BIAS_CORRECTION=TRUE::NORMALIZATION_METHOD=global.normalization.median::PEPTIDES='.$peptideFilterStrg.'::RATIO_TYPE=SimpleRatio';

		my $globalQuantiID;
		if ($numSamples > 1) { # -> Global ratio(s) != sample-specific ratio(s)
			my $quantifAnnot=$baseQuantifAnnot;
			#>RATIOS
			my $ratioStrg='';
			foreach my $refIdx (0..$#labelList-1) {
				foreach my $testIdx (1..$#labelList) {
					next if $testIdx==$refIdx;
					$ratioStrg.=';' if $ratioStrg;
					$ratioStrg.='#'.$designInfo{'LABEL_EXPCOND'}{$labelList[$testIdx]}{'ID'}.'/#'.$designInfo{'LABEL_EXPCOND'}{$labelList[$refIdx]}{'ID'};
				}
			}
			$quantifAnnot.='::RATIOS='.$ratioStrg;
			#>STATES
			my $labelStateStrg='';
			#$globalQuantiID=$quantifID+$numSamples+1; # use highest DB ID
			$numProtQuantif++;
			foreach my $targetPos (sort{$a<=>$b} keys %{$designInfo{'ALL_STATES'}}) {
				$labelStateStrg.=';' if $labelStateStrg;
				my $label=$labelList[$targetPos-1];
				$labelStateStrg.=(scalar @{$designInfo{'ALL_STATES'}{$targetPos}}).','.join('.',@{$designInfo{'ALL_STATES'}{$targetPos}}).',#'.$designInfo{'LABEL_EXPCOND'}{$label}{'ID'};
			}
			$quantifAnnot.='::STATES='.$labelStateStrg;
			$sthInsQ->execute($ratioQuantMethID,$designID,"Global ratios [$importDate]",'protein',$quantifAnnot) || die $sthInsQ->errstr();
			$globalQuantiID=$dbh->last_insert_id(undef,undef,'QUANTIFICATION','ID_QUANTIFICATION');
			#>EXPCONDITION_QUANTIF
			foreach my $label (keys %{$designInfo{'LABEL_EXPCOND'}}) {
				$sthInsEQ->execute($designInfo{'LABEL_EXPCOND'}{$label}{'ID'},$globalQuantiID) || die $sthInsEQ->errstr();
			}
			$sthInsParQ->execute($peptideQuantifID,$globalQuantiID) || die $sthInsParQ->errstr();
			push @protQuantifList,$globalQuantiID;
		}

		my @allStatesAnnot;
		foreach my $sampleName (sort{$sampIDs{$a}<=>$sampIDs{$b}} keys %{$designInfo{'ExpconditionOriginalName'}}) {
			my $sampleNameStrg=($sampleName eq 'No experiment')? '' : " $sampleName";

			my $quantifAnnot=$baseQuantifAnnot;
			#>RATIOS
			my $ratioStrg='';
			foreach my $refIdx (0..$#labelList-1) {
				my $refExpCondName=$labelList[$refIdx].$sampleNameStrg;
				foreach my $testIdx (1..$#labelList) {
					next if $testIdx==$refIdx;
					my $testExpCondName=$labelList[$testIdx].$sampleNameStrg;
					$ratioStrg.=';' if $ratioStrg;
					$ratioStrg.='#'.$designInfo{'Expcondition'}{$testExpCondName}{'ID'}.'/#'.$designInfo{'Expcondition'}{$refExpCondName}{'ID'};
				}
			}
			$quantifAnnot.='::RATIOS='.$ratioStrg;
			#>STATES
			$quantifAnnot.='::STATES='.join(';',@{$designInfo{'STATES'}{$sampleName}});

			$numProtQuantif++;
			$sthInsQ->execute($ratioQuantMethID,$designID,"$sampleName [Ratios]",'protein',$quantifAnnot) || die $sthInsQ->errstr();
			my $quantifID=$dbh->last_insert_id(undef,undef,'QUANTIFICATION','ID_QUANTIFICATION');
			$sthInsParQ->execute($peptideQuantifID,$quantifID) || die $sthInsParQ->errstr();

			push @protQuantifList,$quantifID;

			my %linkedAna;
			foreach my $expCondName ( keys  %{$designInfo{'ExpconditionOriginalName'}{$sampleName}} ) {
				$sthInsEQ->execute($designInfo{'Expcondition'}{$expCondName}{'ID'},$quantifID) || die $sthInsEQ->errstr();
				if ($designInfo{'Fraction'} && $designInfo{'Fraction'}{$expCondName}) {
					foreach my $fraction (keys %{$designInfo{'Fraction'}{$expCondName}}){
						foreach my $analysisID (keys %{$designInfo{'Fraction'}{$expCondName}{$fraction}}){
	#next unless $ana2pepQuantification{$analysisID}; # No peptide quantif data for this analysis -> skip (PP 08/11/16)
							next if $linkedAna{$analysisID}; # do only once
							$sthInsAQ->execute($quantifID,$analysisID) || die $sthInsAQ->errstr();
							#$sthInsParQ->execute($ana2pepQuantification{$analysisID},$quantifID) || die $sthInsParQ->errstr();
							if ($globalQuantiID) {
								$sthInsAQ->execute($globalQuantiID,$analysisID) || die $sthInsAQ->errstr();
								#$sthInsParQ->execute($ana2pepQuantification{$analysisID},$globalQuantiID) || die $sthInsParQ->errstr();
							}
							$linkedAna{$analysisID}=1;
						}
					}
				}
				else {
					#foreach my $analysisID (keys %{$designInfo{'Experiment'}{$expCondName}}){
					foreach my $analysisID (keys %{$designInfo{'Expcondition'}{$expCondName}{'ANALYSIS'}}) {
	#next unless $ana2pepQuantification{$analysisID}; # No peptide quantif data for this analysis -> skip (PP 08/11/16)
						next if $linkedAna{$analysisID}; # do only once
						$sthInsAQ->execute($quantifID,$analysisID) || die $sthInsAQ->errstr();
						$sthInsParQ->execute($ana2pepQuantification{$analysisID},$quantifID) || die $sthInsParQ->errstr();
						if ($globalQuantiID) {
							$sthInsAQ->execute($globalQuantiID,$analysisID) || die $sthInsAQ->errstr();
							#$sthInsParQ->execute($ana2pepQuantification{$analysisID},$globalQuantiID) || die $sthInsParQ->errstr();
						}
						$linkedAna{$analysisID}=1;
					}
				}
			}

			###> Make parameters ratio-specific (eg H/L)
			my $targetPos=0;
			foreach my $refIdx (0..$#labelList-1) {
				foreach my $testIdx (1..$#labelList) {
					next if $testIdx==$refIdx;
					$targetPos++;
					my $ratioName=$labelList[$testIdx].'/'.$labelList[$refIdx];
					foreach my $paramCode (@ratioParamCodes) {
						my $paramText="$quantifParamCodes{$paramCode}$sampleNameStrg";
						$paramText=~s/#%#/$ratioName/;
						push @{$recordedParams{$paramText}},[$quantifID,$paramCode,$targetPos]; # push because same param used for different quantifs (in case of ratios)
					}
					$requiredParams{$quantifID}{$targetPos}="Ratio $ratioName normalized";
				}
			}
			##>These peptide counts are not channel-specific
			foreach my $paramCode (@allStatesParamCodes) { # "RATIO:$paramCode" to distinguish PROT_RATIO_PEP parameters from MQ parameters
				push @{$recordedParams{"$quantifParamCodes{$paramCode}$sampleNameStrg"}},[$quantifID,"RATIO:$paramCode",0]; # push because same param used for different quantifs (in case of ratios)
			}
		}

		####>Global ratio quantif
		if ($globalQuantiID) {
			my $targetPos=0;
			foreach my $refIdx (0..$#labelList-1) {
				foreach my $testIdx (1..$#labelList) {
					next if $testIdx==$refIdx;
					$targetPos++;
					my $ratioName=$labelList[$testIdx].'/'.$labelList[$refIdx];
					foreach my $paramCode (@ratioParamCodes) {
						my $paramText=$quantifParamCodes{$paramCode};
						$paramText=~s/#%#/$ratioName/;
						push @{$recordedParams{$paramText}},[$globalQuantiID,$paramCode,$targetPos]; # push because same param used for different quantifs (in case of ratios)
					}
					$requiredParams{$globalQuantiID}{$targetPos}="Ratio $ratioName normalized";
				}
			}
			##>These peptide counts are not channel-specific
			foreach my $paramCode (@allStatesParamCodes) {
				push @{$recordedParams{$quantifParamCodes{$paramCode}}},[$globalQuantiID,"RATIO:$paramCode",0]; # push because same param used for different quantifs (in case of ratios)
			}
		}
	}
	$sthInsEQ->finish;
	$sthInsParQ->finish;

}
$sthInsAQ->finish;
$sthInsQ->finish;
$dbh->commit;
#$dbh->rollback; # DEBUG

print "&nbsp-1 peptide and $numProtQuantif protein Quantifications imported.<BR>\n";
print "Done.</FONT><BR>\n";


#$dbh->disconnect; exit; # DEBUG

########################
####>msmms.txt file<####
########################
my $counter=0;
my %msmsMissedCut; # missed cleavage data missing in evidence.txt before MaxQuant v~1.5 but exists in msms.txt
if (-e "$tmpFilesDir/$msmsFile") {
	print "<BR><FONT class='title3'>Reading msms file...";
	open (MSMS,"$tmpFilesDir/$msmsFile") || die "Unable to open $tmpFilesDir/$msmsFile";
	my $firstmsmsLine=<MSMS>;
	my @parameters=split(/ *\t */,$firstmsmsLine);
	my %msmsColNum;
	my $ncol=0;
	foreach my $colName (@parameters) {
		$msmsColNum{$colName}=$ncol;
		$ncol++;
	}
	###> Copy the msmsFile in all peptide data directories
	foreach my $rawName (keys %anaInfo) {
		my $anaID=$anaInfo{$rawName}{'ID_ANALYSIS'};
		next unless $anaID;
		open($anaInfo{$rawName}{'INFILE'},'>',"$promsPath{peptide}/proj_$projectID/ana_$anaID/msms.txt") or die;
		my $fhavalue=$anaInfo{$rawName}{'INFILE'};
		print $fhavalue $firstmsmsLine;
	}
	while (<MSMS>) {
		$counter++;
		my @infos=split(/ *\t */,$_); # remove starting/trailing spaces
		#$infos[0]=~s/\s+$//; # remove trailing spaces
		my $fhavalue=$anaInfo{$infos[0]}{'INFILE'};
		$msmsMissedCut{$infos[$msmsColNum{'Sequence'}]}=$infos[$msmsColNum{'Missed cleavages'}];
		print $fhavalue $_;
		if ($counter > 10000) {
			print '.';
			$counter=0;
		}
	}
	close MSMS;
	foreach my $rawName (keys %anaInfo) {
		my $fhavalue=$anaInfo{$rawName}{'INFILE'};
		next unless $fhavalue;
		close($fhavalue);
		delete($anaInfo{$rawName}{'INFILE'})
	}
	print " Done.</FONT><BR>\n";
}


###########################
####>evidence.txt file<####
###########################
###> TODO : check DECOY + FDR + MAX_RANK for ANALYSIS
print "<BR><FONT class='title3'>Reading evidence file...";
open (EVIDENCE,"$tmpFilesDir/$evidenceFile")  || die "Unable to open $tmpFilesDir/$evidenceFile";
#my @evidenceColumns=('Sequence','Length','Modifications','Modified sequence','Missed cleavages','Proteins','Leading proteins','Leading razor protein',
#                     'Raw file','Charge','m/z','Mass','Mass Error [ppm]','Mass Error [Da]','Retention time','Calibrated retention time',
#                     'PEP','MS/MS Count','MS/MS Scan Number','Score','Intensity','Peptide ID');
my (%evidenceColNum,%pepInfo,%queryNum,%pepVmod);
my (%evidence2peptide,%maxProtMatch,%actualSeq2seq,%protDbRank,%matchList,%bestScore,%proteinRank); #,%bestScorebad
my (%seq2posString,%peptideXIC,%isGhost,%featureCount);
my %raw2UsedIdentifiers;
my $missedCutData=1;
$counter=0;
while (my $line=<EVIDENCE>) {
	$counter++;
	chomp($line);
	my @parameters=split(/ *\t */,$line); # remove starting/trailing spaces
	if ($.==1) { # 1st line of the file
		my $ncol=0;
		foreach my $colName (@parameters) {
			$evidenceColNum{uc($colName)}=$ncol; # Upper case for all column headers!!!
			$ncol++;
		}
		$missedCutData=0 unless $evidenceColNum{uc('Missed cleavages')};
		next;
	}
	next unless $anaInfo{$parameters[$evidenceColNum{uc('Raw file')}]}; # parameter group filtering??????????????
	my $analysisID=$anaInfo{$parameters[$evidenceColNum{uc('Raw file')}]}{'ID_ANALYSIS'};
	$queryNum{$analysisID}++;
	###> $parameters[$evidenceColNum{'Proteins'}] can be empty if only one REV__ protein is matched !
	my $actualSequence=$parameters[$evidenceColNum{uc('Modified sequence')}];
	my ($matchGood,$matchBad,$matchCON)=(0,0,0);
	my $proteinList=($parameters[$evidenceColNum{uc('Proteins')}])? $parameters[$evidenceColNum{uc('Proteins')}] : $parameters[$evidenceColNum{uc('Leading proteins')}];
#my $protRank=1;
	foreach my $rawIdentifier (split(/;/,$proteinList)) {
		my $identifier=$rawIdentifier; # default
		$raw2UsedIdentifiers{$rawIdentifier}=$rawIdentifier; # default
		if ($identifier =~ /REV__/) {
			#$identifier =~ s/REV__/DECOY_/g;
			$matchBad=1;
		}
		elsif ($excludeCON && $identifier=~/CON__/) {
			$matchCON=1;
		}
		else {
			$matchGood=1;
			my $tmpIdenfitier;
			foreach my $dbID (@databankIDs) {
				($tmpIdenfitier)=($rawIdentifier=~/$idParseRules{$dbID}/);
				if ($tmpIdenfitier) {
					$identifier=$raw2UsedIdentifiers{$rawIdentifier}=$tmpIdenfitier;
					last;
				}
			}

#$protDbRank{$identifier}=($identifier=~/CON__/ && $contaminantDB)? 2 : 1;
#$matchList{$protDbRank{$identifier}}{$identifier}{$analysisID}=1;
			$matchList{$identifier}{$analysisID}=1;
			$maxProtMatch{$analysisID}{$identifier}{$actualSequence}++;

#$proteinRank{$identifier}=$protRank if (!$proteinRank{$identifier} || $proteinRank{$identifier} > $protRank); # needed to select best MQ prot during MG re-creation (pseudo MG! Better to read proteinGroups.txt)
			unless ($proteinRank{$identifier}) {
				$proteinRank{$identifier}=($rawIdentifier eq $parameters[$evidenceColNum{uc('Leading Razor Protein')}])? 1 : 2; # needed to select best MQ prot during MG re-creation
			}
		}
#$protRank++;
#$protDbRank{$identifier}=($identifier=~/CON__/ && $contaminantDB)? 2 : 1;
#$matchList{$protDbRank{$identifier}}{$identifier}{$analysisID}=1;
#$maxProtMatch{$analysisID}{$identifier}{$actualSequence}++;
	}
	if ($matchBad) {
		$featureCount{'DECOY'}{$analysisID}{$actualSequence}++; # for FDR
#$bestScorebad{$analysisID}{$actualSequence}=$score if ($score && $bestScorebad{$analysisID} || !$bestScorebad{$analysisID}{$actualSequence});
#$bestScorebad{$analysisID}{$actualSequence}=$score if (!defined($bestScorebad{$analysisID}{$actualSequence}) || ($bestScorebad{$analysisID}{$actualSequence} && $score>$bestScorebad{$analysisID}{$actualSequence}));
# NOT NEEDED # @{$pepInfo{$analysisID}{"-$queryNum{$analysisID}"}}=($actualSequence,$parameters[$evidenceColNum{uc('Sequence')}],$parameters[$evidenceColNum{uc('Length')}],"-$queryNum{$analysisID}",$score,$missedCut,$massExp,$parameters[$evidenceColNum{uc('Mass')}],$parameters[$evidenceColNum{uc('m/z')}],$massErrorDa,$charge,$etString,$validStatus,"PEP=$pep",$parameters[$evidenceColNum{uc('MS/MS Count')}]);
# NOT NEEDED # @{$evidence2peptide{$parameters[$evidenceColNum{uc('id')}]}}=($analysisID,"-$queryNum{$analysisID},$actualSequence") unless $evidence2peptide{$actualSequence}{$parameters[$evidenceColNum{uc('id')}]};
#print "Bad '$actualSequence' $analysisID -$queryNum{$analysisID} score$bestScorebad{$analysisID}{$actualSequence}<BR>\n";
	}
	elsif ($matchCON && !$matchGood) {
		$featureCount{'TARGET'}{$analysisID}{$actualSequence}++; # for FDR
	}
	next unless $matchGood;
	#if ($matchGood) {
	$featureCount{'TARGET'}{$analysisID}{$actualSequence}++; # for FDR

	#>Missed Cleavages
	my $missedCut=0;
	if ($missedCutData) {$missedCut=$parameters[$evidenceColNum{uc('Missed cleavages')}];}
	elsif ($msmsMissedCut{$parameters[$evidenceColNum{uc('Sequence')}]}) {$missedCut=$msmsMissedCut{$parameters[$evidenceColNum{uc('Sequence')}]};}
	#>Scores
	my ($score,$pepDataStrg);
	if ($parameters[$evidenceColNum{uc('Score')}] =~ /\d/) {
		$score=$parameters[$evidenceColNum{uc('Score')}]*1;
		$pepDataStrg='DELTA_SC='.($parameters[$evidenceColNum{uc('Delta score')}]*1);
		$pepDataStrg.='##PEP='.$parameters[$evidenceColNum{'PEP'}] if $parameters[$evidenceColNum{'PEP'}] =~ /\d/;
		$pepDataStrg.='##'; # for MQ_EVI below
	}
	else {$score=0;}
	$pepDataStrg.='MQ_EVI='.$parameters[$evidenceColNum{uc('id')}];
	my $charge=$parameters[$evidenceColNum{uc('Charge')}];
	#$queryNum{$analysisID}++;
	my $massExp=($parameters[$evidenceColNum{uc('m/z')}]-1.007825032)*$charge; # same way computed in storeAnalyses.cgi for Paragon
	my $massErrorDa=($evidenceColNum{uc('Mass Error [Da]')})? $parameters[$evidenceColNum{uc('Mass Error [Da]')}] : ($parameters[$evidenceColNum{uc('Mass Error [ppm]')}])? $parameters[$evidenceColNum{uc('Mass Error [ppm]')}]*$parameters[$evidenceColNum{uc('Mass')}]/1000000 : undef;
	#>Check for isobaric modif (Not listed as variable modif!!!!)
	if ($isobarModifInfo{ID} && $massErrorDa > 10) {
		my $numModif=int(0.5+($massErrorDa*$charge/$isobarModifInfo{MONO_MASS}));
		if ($numModif >= 1) {
			my ($Nterm,$protNterm,@positions);
			my $massShift=($isobarModifInfo{MONO_MASS}/$charge)+1.007825032;
			foreach my $res (sort keys %{$isobarModifInfo{RES}}) {
				if ($res eq '=') {
					$Nterm=1;
					$numModif--;
					$massErrorDa-=$massShift;
				}
				elsif ($res eq '-') {
					$protNterm=1;
					$numModif--;
					$massErrorDa-=$massShift;
				}
				else {
					while ($parameters[$evidenceColNum{uc('Sequence')}]=~/$res/g && $numModif > 0) {
						push @positions,$-[0]+1;
						$massErrorDa-=$massShift;
						$numModif--;
					}
				}
				last if $numModif==0;
			}
			@positions=sort{$a<=>$b} @positions;
			if ($Nterm) {unshift @positions,'=';} # 1rst element
			elsif ($protNterm) {unshift @positions,'-';} # assumes N-term and Protein N-term are incompatible
			$isobarModifInfo{ANA}{$analysisID}{$queryNum{$analysisID}}=join('.', @positions);
		}
	}

	my $etString=($parameters[$evidenceColNum{uc('MS/MS Scan Number')}] =~ /\d/)? "et$parameters[$evidenceColNum{'CALIBRATED RETENTION TIME'}];sc$parameters[$evidenceColNum{'MS/MS SCAN NUMBER'}];":"et$parameters[$evidenceColNum{'CALIBRATED RETENTION TIME'}];";
	$actualSeq2seq{$actualSequence}=$parameters[$evidenceColNum{uc('Sequence')}];
	#if ($matchGood) {#}
	my $validStatus=0;
	if (!$bestScore{$analysisID} || !defined($bestScore{$analysisID}{$actualSequence})) {
		$validStatus=1 if $score;
		$bestScore{$analysisID}{$actualSequence}=$score; # ranked by descending scores in file
	}
	$isGhost{$analysisID}{$queryNum{$analysisID}}=1 unless $validStatus;
	#$validStatus=1 if ($score && !$bestScore{$analysisID}{$actualSequence}); # if there is no score, then, it is a ghost peptide
	@{$pepInfo{$analysisID}{$queryNum{$analysisID}}}=($actualSequence,$parameters[$evidenceColNum{uc('Sequence')}],$parameters[$evidenceColNum{uc('Length')}],$queryNum{$analysisID},$score,$missedCut,$massExp,$parameters[$evidenceColNum{uc('Mass')}],$parameters[$evidenceColNum{uc('m/z')}],$massErrorDa,$charge,$etString,$validStatus,$pepDataStrg); # ,$parameters[$evidenceColNum{uc('MS/MS Count')}]
	@{$evidence2peptide{$parameters[$evidenceColNum{uc('id')}]}}=($analysisID,$queryNum{$analysisID},$actualSequence);
	if ($labeling eq 'SILAC') {
		foreach my $colLabel ('L','M','H') {
			$peptideXIC{$analysisID}{$queryNum{$analysisID}}{'XIC'}{$colLabel}=$parameters[$evidenceColNum{uc("Intensity $colLabel")}] if $evidenceColNum{uc("Intensity $colLabel")};
		}
	}
	elsif ($isobarModifInfo{ID}) { # $labeling=~/ITRAQ|TMT/
		foreach my $repIdx (0..($labelingPlex-1)) {
			$peptideXIC{$analysisID}{$queryNum{$analysisID}}{'XIC'}{$repIdx}=$parameters[$evidenceColNum{uc("Reporter intensity $correctedStrg $repIdx")}] if ($evidenceColNum{uc("Reporter intensity $correctedStrg $repIdx")});
		}
	}
	elsif ($labeling eq 'FREE') {
		$peptideXIC{$analysisID}{$queryNum{$analysisID}}{'XIC'}{'NONE'}=$parameters[$evidenceColNum{uc('Intensity')}] if $parameters[$evidenceColNum{uc('Intensity')}];
	}
#$featureCount{'TARGET'}{$analysisID}{$actualSequence}++;
	#print "Good '$actualSequence' $analysisID $queryNum{$analysisID} $parameters[$evidenceColNum{'id'}] score=$bestScore{$analysisID}{$actualSequence}<BR>\n";
	#}
	foreach my $vmod (keys %{$anaVmods{$analysisID}}) {
		if ($parameters[$evidenceColNum{uc($vmod)}] && $evidenceColNum{uc($vmod)} >= 1) {
			@{$pepVmod{$analysisID}{$queryNum{$analysisID}}{$vmod}}=(undef); # default
			#next if $seq2posString{$actualSequence}{$vmod}; # If already computed, no need to search for position...
			my $posString;
			if ($evidenceColNum{uc("$vmod Probabilities")}) {
				my $probSequence=$parameters[$evidenceColNum{uc("$vmod Probabilities")}]; # "C(0.837)DPRLGKYMAC(0.16)C(0.003)LLYR"
				my @probValues= $probSequence =~ /\(([^\)]+)\)/g; # @probValues = ((0.837) , (0.16) , (0.003))
				my @slices=split(/\([^\)]+\)/,$probSequence);
				my $pos=0;
				my (%probMod,@positions);
				my $probString='##PRB_MQ=';
				my $keepProb=0;
				#for (my $i=0; $i < $#slices; $i++) {
				for (my $i=0; $i <= $#slices; $i++) { # Change on 2017/12/06 because of GlyGly modification with probability on last residue !!!
					last if $i > $#probValues;# Be careful, 2 cases for GlyGly: 'AAAAAAALQAK(1)SDEK(1)AAVAGK' or 'AAAAAAALQAK(0.714)SDEK(0.286)'
					$pos+=length($slices[$i]);
					$probMod{$pos}=$probValues[$i];
					$probString.=',' if $i > 0;
					$probString.=$pos.':'.$probValues[$i];
					$keepProb=1 if ($probValues[$i] > 0 && $probValues[$i] < 1); # no need to keep prob with only 100% & 0%
				}
				$pepVmod{$analysisID}{$queryNum{$analysisID}}{$vmod}[0]=$probString if $keepProb;
				next if $seq2posString{$actualSequence}{$vmod}; # If already computed, no need to search for position...
				my $nbMod=$parameters[$evidenceColNum{uc($vmod)}];
				foreach my $pos (sort {$probMod{$b} <=> $probMod{$a} || $a<=>$b} keys %probMod) {
					last if $nbMod == 0;
					$nbMod--;
					push @positions,$pos;
				}
				$posString=join('.', sort{$a<=>$b} @positions);
			}
			else {
				if ($vmod =~ /C-TERM/i) {
					$posString=($vmod =~ /PROTEIN/i) ? '+' : '*';
				}
				if ($vmod =~ /N-TERM/i) {
					$posString=($vmod =~ /PROTEIN/i) ? '-' : '=';
				}
			}
			$seq2posString{$actualSequence}{$vmod}=$posString;
		}
	}
	if ($counter>100000) {
		print '.';
		$counter=0;
	}
}
close EVIDENCE;

##>Update FDR info
my $sthUpAna=$dbh->prepare("UPDATE ANALYSIS SET FALSE_DISCOVERY=? WHERE ID_ANALYSIS=?");
foreach my $analysisID (keys %maxProtMatch) {
	my %pepCount=(TARGET=>0,DECOY=>0);
	my %specCount=(TARGET=>0,DECOY=>0);
	foreach my $type (keys %pepCount) {
		next if (!$featureCount{$type} || !$featureCount{$type}{$analysisID});
		$pepCount{$type}=scalar keys %{$featureCount{$type}{$analysisID}};
		foreach my $actualSequence (keys %{$featureCount{$type}{$analysisID}}) {$specCount{$type}+=$featureCount{$type}{$analysisID}{$actualSequence};}
		$specCount{$type}-=$pepCount{$type};
	}
	$sthUpAna->execute("$pepCount{DECOY}:$specCount{DECOY}:$pepCount{TARGET}:$specCount{TARGET}",$analysisID); # numDecoy:(specDecoy-numDecoy):numTarget:(spectarget-numTarget)
}
$sthUpAna->finish;

print " Done.</FONT><BR>\n";

print "<BR><FONT class='title3'>Preprocessing protein list...";
my ($protVisibility,$projectIdentMapID)=$dbh->selectrow_array("SELECT PROT_VISIBILITY,ID_IDENTIFIER FROM PROJECT WHERE ID_PROJECT=$projectID");
$protVisibility=0 unless $protVisibility;
$projectIdentMapID=0 unless $projectIdentMapID;
my ($giIdentID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GI'");
my (%projectProt,%bestProtVis); #,%proteinLength,%incompleteProtein,%noDescriptionProtein);
#my $sthP=$dbh->prepare("SELECT ID_PROTEIN,IDENTIFIER,PROT_LENGTH,PROT_DES FROM PROTEIN WHERE ID_PROJECT=$projectID");
my $sthP=($protVisibility)?
	$dbh->prepare("SELECT P.ID_PROTEIN,IDENTIFIER,MAX(VISIBILITY) FROM PROTEIN P,ANALYSIS_PROTEIN AP WHERE P.ID_PROTEIN=AP.ID_PROTEIN AND ID_PROJECT=$projectID GROUP BY P.ID_PROTEIN")
	: $dbh->prepare("SELECT ID_PROTEIN,IDENTIFIER FROM PROTEIN WHERE ID_PROJECT=$projectID"); # no need for best vis of $protVisibility==0
$sthP->execute || die $sthP->errstr;
$counter=0;
while (my ($protID,$identifier,$maxVis)=$sthP->fetchrow_array) {
	next unless $identifier; # identifier is null (previous ID protection not removed)
	$counter++;
	$projectProt{$identifier}=$protID;
	$bestProtVis{$identifier}=$maxVis // 0;
	if ($counter>100000) {
		print '.';
		$counter=0;
	}
}
$sthP->finish;
print '/.';
$dbh->do("DELETE FROM PROTEIN WHERE ID_PROJECT=$projectID AND IDENTIFIER IS NULL") || $dbh->errstr; # in case of error in previous process (rare!)
#my $insProt=$dbh->prepare("INSERT INTO PROTEIN (ID_PROTEIN,ID_PROJECT) VALUES (?,$projectID)") || $dbh->errstr;
##>Counting new proteins
my %newProteins;
$counter=0;
foreach my $identifier (keys %matchList) {
	foreach my $analysisID (keys %{$matchList{$identifier}}) {
		##if ($projectProt{$identifier}) { # new proteins added to %projectProt
		##	$matchList{$identifier}{$analysisID}=$projectProt{$identifier};
		##}
		##else {
		##	$newProteins{$identifier}=1;
		##	$bestProtVis{$identifier}=0;
		##	#$matchList{$identifier}{$analysisID}=1;
		##}
		unless ($projectProt{$identifier}) {
			$newProteins{$identifier}=1;
			$bestProtVis{$identifier}=0;
		}
		$counter++;
		if ($counter>100000) {
			print '.';
			$counter=0;
		}
	}
}
##my ($maxProteinID)=$dbh->selectrow_array("SELECT MAX(ID_PROTEIN) FROM PROTEIN");
##my $protectProtID = $maxProteinID + 2 + scalar keys %newProteins; # +1 should be enough
##$dbh->do("INSERT INTO PROTEIN (ID_PROTEIN,ID_PROJECT) VALUES ($protectProtID,$projectID)") || $dbh->errstr; # ID protection
##$dbh->commit;

##foreach my $identifier (keys %matchList) {
##	foreach my $analysisID (keys %{$matchList{$identifier}}) {
##		if ($projectProt{$identifier}) { # new proteins added to %projectProt
##			$matchList{$identifier}{$analysisID}=$projectProt{$identifier};
##		}
##		else {
##			$matchList{$identifier}{$analysisID}=++$maxProteinID;
##			$projectProt{$identifier}=$maxProteinID;
##			$newProteins{$identifier}=$maxProteinID;
##			$bestProtVis{$identifier}=0;
##		}
##		$counter++;
##		if ($counter>100000) {
##			print '.';
##			$counter=0;
##		}
##	}
##}

print " Done.</FONT><BR>\n";

print "<BR><FONT class='title3'>Scanning databank(s):</FONT><BR>\n";
my (%contMatchList,%protDes,%protMW,%protOrg,%protLength);
my (@anaIDs)=sort{$a<=>$b} keys %bestScore;
my $prevRefMatchList;
my $dbRank=0;
foreach my $dbID (@databankIDs) {
	$dbRank++;
	my $prefix;
	my $refMatchList={};
	if ($dbRank==1) {
		if ($contaminantDB) { # separate proteins from contaminants
			my $contDbRank=scalar @databankIDs; # last in list
			foreach my $identifier (keys %matchList) {
				if ($identifier=~/^CON__/) {
					$contMatchList{$identifier}=$matchList{$identifier};
					$protDbRank{$identifier}=$contDbRank;
				}
				else {$refMatchList->{$identifier}=$matchList{$identifier};} # add to list to be scanned
			}
		}
		else {$refMatchList=\%matchList;}
	}
	else {
		foreach my $identifier (keys %{$prevRefMatchList}) {
			if ($protLength{$identifier}) {$protDbRank{$identifier}=$dbRank-1;} # matched during previous db scan
			elsif ($dbID==$contaminantDB) {$protDbRank{$identifier}=1;} # set to 1 because no more db scan
			else {$refMatchList->{$identifier}=$matchList{$identifier};} # add to list to be scanned
		}
		if ($dbID==$contaminantDB) {
			$refMatchList=\%contMatchList;
			$prefix='CON__';
		}
	}
#print '.';
	#my $prefix=($dbID==$contaminantDB)? 'CON__' : undef;
	#&promsMod::getProtInfo('silent',$dbh,$dbID,\@anaIDs,\%protDes,\%protMW,\%protOrg,\%protLength,\%{$matchList{$dbRank}},undef,$prefix);
	&promsMod::getProtInfo('verbose',$dbh,$dbID,\@anaIDs,\%protDes,\%protMW,\%protOrg,\%protLength,$refMatchList,undef,$prefix) if scalar keys %{$refMatchList};
	$prevRefMatchList=$refMatchList;
#print '.';
}
unless ($contaminantDB) {
	foreach my $identifier (keys %{$prevRefMatchList}) {
		$protDbRank{$identifier}=1 unless $protDbRank{$identifier};
	}
}
print "<FONT class='title3'>Done.</FONT><BR>\n";


###> Compute Match Info based on analysis.fasta files
print "<BR><FONT class='title3'>Processing peptide/protein match data...";
my (%numMatches,%sequenceList);
$counter=0;
foreach my $analysisID (@anaIDs) {
	if (-e "$promsPath{peptide}/proj_$projectID/ana_$analysisID/analysis.fasta") {
		my $identifier; #,$des,$org,$mw,$length;
		my $sequence='';
		open (FAS,"$promsPath{peptide}/proj_$projectID/ana_$analysisID/analysis.fasta") || die "Unable to open $promsPath{peptide}/proj_$projectID/ana_$analysisID/analysis.fasta";
		while(<FAS>) {
			$counter++;
			if (/^>(\S+)/) {
				my $newIdentifier=$1;
				if ($sequence && $maxProtMatch{$analysisID}{$identifier}) { # test
					$sequenceList{$identifier}=$sequence;
					foreach my $actualSequence (keys %{$maxProtMatch{$analysisID}{$identifier}}) {
						my $pepSeq=$actualSeq2seq{$actualSequence};
						#my ($seqBefore,$resEnd)=($sequence=~/(.*)$pepSeq(\w?)/);
						if (!$numMatches{$analysisID} || !$numMatches{$analysisID}{$identifier} || !$numMatches{$analysisID}{$identifier}{$pepSeq}) {
							while ($sequence=~/(\w?)$pepSeq(\w?)/g) { # potential multiple matches
								my ($resBeg,$resEnd)=($1,$2);
								my $startPos=$-[0] + 1;
								my $endPos=$+[0];
								if ($resBeg) {$startPos++;}
								else {$resBeg='-';}
								if ($resEnd) {$endPos--;}
								else {$resEnd='-';}
				#@{$matchInfos{$analysisID}{$identifier}{$actualSequence}{$startPos}}=($resBeg,$resEnd,$endPos);
								@{$numMatches{$analysisID}{$identifier}{$pepSeq}{$startPos}}=($resBeg,$resEnd,$endPos);
							}
						}
					}
				}
				$identifier=$newIdentifier;
				$sequence='' if defined($identifier);
            }
			else {
				#chomp;
				$_=~s/\W+//g; # chomp not enough? (Windows!)
				$sequence.=$_;
			}
			if ($counter>10000) {
				print '.';
				$counter=0;
			}
		}
		close FAS;

		if ($sequence && $maxProtMatch{$analysisID}{$identifier}) {
			$sequenceList{$identifier}=$sequence;
			foreach my $actualSequence (keys %{$maxProtMatch{$analysisID}{$identifier}}) {
				my $pepSeq=$actualSeq2seq{$actualSequence};
				if (!$numMatches{$analysisID} || !$numMatches{$analysisID}{$identifier} || !$numMatches{$analysisID}{$identifier}{$pepSeq}) {
					while ($sequence=~/(\w?)$pepSeq(\w?)/g) { # potential multiple matches
						my ($resBeg,$resEnd)=($1,$2);
						my $startPos=$-[0] + 1;
						my $endPos=$+[0];
						if ($resBeg) {$startPos++;}
						else {$resBeg='-';}
						if ($resEnd) {$endPos--;}
						else {$resEnd='-';}
						#@{$matchInfos{$analysisID}{$identifier}{$actualSequence}{$startPos}}=($resBeg,$resEnd,$endPos);
						@{$numMatches{$analysisID}{$identifier}{$pepSeq}{$startPos}}=($resBeg,$resEnd,$endPos);
					}
				}
			}
		}
		unlink "$promsPath{peptide}/proj_$projectID/ana_$analysisID/analysis.fasta"; # no longer needed
    }
}
print " Done.</FONT><BR>\n";
#die "DEBUG";

###> Read peptide file to get context and be able to write MATCH_MULTI and MATCH_INFO pos,AA-1,AA+1
###> 2nd step : create a Design, Condition, Observation, etc.
print "<BR><FONT class='title3'>Reading peptide file...";
my (%peptideColNum,%pepProteins);
open (PEPTIDE,"$tmpFilesDir/$peptideFile") || die "Unable to open $tmpFilesDir/$peptideFile";
$counter=0;
while (my $line=<PEPTIDE>) {
	$counter++;
	chomp($line);
	my @parameters=split(/ *\t */,$line); # remove starting/trailing spaces
	if ($.==1) { # 1st line of the file
		my $ncol=0;
		foreach my $colName (@parameters) {
			$peptideColNum{$colName}=$ncol;
			$ncol++;
		}
		next;
	}
	#'Proteins'
	#'Evidence IDs' --> Get the evidence that gave this peptide
	my $rawRazorProtein=$parameters[$peptideColNum{'Leading razor protein'}];
	$rawRazorProtein =~ s/REV__/DECOY_/;
	my $razorProtein = $raw2UsedIdentifiers{$rawRazorProtein} || $rawRazorProtein;
	foreach my $rawIdentifier ( split(/;/,$parameters[$peptideColNum{'Proteins'}]) ) { # Be careful : CON__ prefix stands for contaminants and REV__ prefix stands for decoy/reverse db
		my $matchGood=($rawIdentifier !~ /REV__/ && ($rawIdentifier !~ /CON__/ || !$excludeCON))? 1 : 0;
		#if ($matchingProt =~ /REV__/ ) {
			$rawIdentifier =~ s/REV__/DECOY_/;
		#}
		my $identifier=$raw2UsedIdentifiers{$rawIdentifier} || $rawIdentifier;
		foreach my $evID (split(/;/,$parameters[$peptideColNum{'Evidence IDs'}])) {
			next unless $evidence2peptide{$evID}; # matches a decoy or excluded contaminant (or skipped parameter group data????????????????????)
			my ($analysisID,$queryNum,$actualSequence)=@{$evidence2peptide{$evID}};
			my $pepSeq=$actualSeq2seq{$actualSequence};
			$pepProteins{$analysisID}{$pepSeq}{$identifier}=1 if $matchGood; # for peptide/protein specificity

			#if (!$numMatches{$analysisID} || !$numMatches{$analysisID}{$identifier} || !$numMatches{$analysisID}{$identifier}{$pepSeq}) {
			#	if ($identifier eq $razorProtein && $peptideColNum{'Start position'}) { # info available in file from v~1.5 only for $razorProtein
			#		my $startPos=$parameters[$peptideColNum{'Start position'}];
			#		my $endPos=$startPos + length($pepSeq)-1;
			#		@{$numMatches{$analysisID}{$identifier}{$pepSeq}{$startPos}}=($parameters[$peptideColNum{'Amino acid before'}],$parameters[$peptideColNum{'Amino acid after'}],$endPos);
			#	}
			#	else {@{$numMatches{$analysisID}{$identifier}{$pepSeq}{0}}=('X','X',0);}
			#}
if ($identifier eq $razorProtein && $peptideColNum{'Start position'}) { # info available in file from v~1.5 only for $razorProtein
	%{$numMatches{$analysisID}{$identifier}{$pepSeq}}=(); # reset in case initialized with fasta file
	my $startPos=$parameters[$peptideColNum{'Start position'}];
	my $endPos=($peptideColNum{'End position'})? $parameters[$peptideColNum{'End position'}] : $startPos + length($pepSeq)-1;
	@{$numMatches{$analysisID}{$identifier}{$pepSeq}{$startPos}}=($parameters[$peptideColNum{'Amino acid before'}],$parameters[$peptideColNum{'Amino acid after'}],$endPos);
}
elsif (!$numMatches{$analysisID} || !$numMatches{$analysisID}{$identifier} || !$numMatches{$analysisID}{$identifier}{$pepSeq}) {
	@{$numMatches{$analysisID}{$identifier}{$pepSeq}{0}}=('X','X',0); # non-razor not mapped in fasta file
}
		}
	}
	if ($counter>100000) {
		print '.';
		$counter=0;
	}
}
close PEPTIDE;
print " Done.</FONT><BR>\n";

print "<BR><FONT class='title3'>Recording peptide data..."; # SLOW!!!!!!!!!!!!!
###> Store Peptide Information
my $sthInsPep=$dbh->prepare("INSERT INTO PEPTIDE (ID_ANALYSIS,PEP_SEQ,PEP_LENGTH,QUERY_NUM,PEP_RANK,SEARCH_RANK,SCORE,MISS_CUT,MR_EXP,MR_CALC,MR_OBS,MR_DELTA,CHARGE,ELUTION_TIME,VALID_STATUS,DATA,SPEC_COUNT) VALUES (?,?,?,?,1,1,?,?,?,?,?,?,?,?,?,?,?)");
my $sthInsPM=$dbh->prepare("INSERT INTO PEPTIDE_MODIFICATION (ID_MODIFICATION,ID_PEPTIDE,MODIF_TYPE,POS_STRING,REF_POS_STRING) VALUES (?,?,'V',?,?)");
my (%pepIDs,%specificSequences);
$counter=0;
foreach my $analysisID (sort{$a<=>$b} keys %pepInfo) {
	foreach my $qNum (sort{$a<=>$b} keys %{$pepInfo{$analysisID}}) {
#next if $queryNum < 0 ; # decoy matches are not kept in myProMS DB
		if ($counter>10000) {
			print '.';
			$counter=0;
		}
		$counter++;
		my ($actualSequence,@data)=@{$pepInfo{$analysisID}{$qNum}};
		##>Ghost?
		#$data[2]=0 if ($isGhost{$analysisID} && $isGhost{$analysisID}{$qNum}); # queryNum
		if ($isGhost{$analysisID} && $isGhost{$analysisID}{$qNum}) {
			next unless $peptideXIC{$analysisID}{$qNum}; # ghost without quantif data
			$data[2]=0;
		}
		##>Peptide is specific ?
		my ($isSpecific)=($pepProteins{$analysisID}{$actualSeq2seq{$actualSequence}} && scalar keys %{$pepProteins{$analysisID}{$actualSeq2seq{$actualSequence}}}==1)? 1:0;
		$specificSequences{$analysisID}{$data[0]}=1 if $isSpecific;
		#print "$peptideID,$analysisID,@data<BR>\n";
		#next unless ($data[3] && $data[3] == $bestScore{$analysisID}{$actualSequence}); # Do not consider lower scoring peptides
		$sthInsPep->execute($analysisID,@data,$featureCount{'TARGET'}{$analysisID}{$actualSequence}) || die $sthInsPep->errstr();
		my $peptideID = $dbh->last_insert_id(undef,undef,'PEPTIDE','ID_PEPTIDE');
		#@{$pepIDs{$analysisID}{$actualSequence}{$peptideID}}=($data[0],$data[11]); # keep peqSeq and validStatus for further coverage computing
		$pepIDs{$analysisID}{$actualSequence}{$peptideID}=$data[11]; # validStatus
		$peptideXIC{$analysisID}{$qNum}{'ID_PEPTIDE'}=$peptideID;
		foreach my $vmod (keys %{$pepVmod{$analysisID}{$qNum}}) {
			$sthInsPM->execute($anaVmods{$analysisID}{$vmod},$peptideID,$seq2posString{$actualSequence}{$vmod},$pepVmod{$analysisID}{$qNum}{$vmod}[0]);
		}
		if ($isobarModifInfo{ANA} && $isobarModifInfo{ANA}{$analysisID} && $isobarModifInfo{ANA}{$analysisID}{$qNum}) { # reconstructed variable isobaric modif positions
			$sthInsPM->execute($isobarModifInfo{ID},$peptideID,$isobarModifInfo{ANA}{$analysisID}{$qNum},undef);
		}
	}
}
print " Done.</FONT><BR>\n";
$sthInsPep->finish;
$sthInsPM->finish;


print "<BR><FONT class='title3'>Recording protein data:</FONT><BR>\n";
my (%maxProtScore,%ppa);
my $sthInsProt=$dbh->prepare("INSERT INTO PROTEIN (ID_PROJECT,IDENTIFIER,ALIAS,PROT_DES,ORGANISM,MW,PROT_SEQ,PROT_LENGTH) VALUES ($projectID,?,?,?,?,?,?,?)");
my $insAP=$dbh->prepare("INSERT INTO ANALYSIS_PROTEIN (ID_ANALYSIS,ID_PROTEIN,SCORE,CONF_LEVEL,DB_RANK,NUM_PEP,NUM_MATCH,PEP_COVERAGE,MATCH_GROUP,PEP_SPECIFICITY,VISIBILITY) VALUES (?,?,?,?,?,?,?,?,?,?,?)");
my $sthTop=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=? AND VISIBILITY=2");
my $sthAtt=$dbh->prepare("INSERT INTO PEPTIDE_PROTEIN_ATTRIB (ID_PEPTIDE,ID_PROTEIN,ID_ANALYSIS,PEP_BEG,PEP_END,FLANKING_AA,IS_SPECIFIC) VALUES (?,?,?,?,?,?,?)");

##>MQ match groups
my (%MQmatchGroup,%MQvisibility);
if ($matchGroupType eq 'MaxQuant') {
	#&createMaxQuantMatchGroups($protVisibility,\%MQmatchGroup,\%MQvisibility,\%projectProt,\%bestProtVis); # same for all analyses
	&createMaxQuantMatchGroups($protVisibility,\%MQmatchGroup,\%MQvisibility,\%matchList,\%bestProtVis,\%raw2UsedIdentifiers); # same for all analyses
}
my @newAnaMapping;
my %protTopMG; # prot is in top of a match group in project
my $totalAna=scalar keys %maxProtMatch;
my $anaCounter=0;
foreach my $analysisID (sort{$a<=>$b} keys %maxProtMatch) {
	$anaCounter++;
	print "&nbsp;-$anaNames{$analysisID} ($anaCounter/$totalAna)...";

	####>Fetching starting ID and protecting ID space for table PROTEIN
	##my ($proteinID)=$dbh->selectrow_array("SELECT MAX(ID_PROTEIN) FROM PROTEIN");
	##$proteinID=0 unless $proteinID;
	##my $maxProteinID=$proteinID + scalar (keys %{$maxProtMatch{$analysisID}}) + 1;
	##$dbh->do("INSERT INTO PROTEIN (ID_PROTEIN,ID_PROJECT) VALUES ($maxProteinID,$projectID)") || $dbh->errstr;
	##$dbh->commit;

	foreach my $identifier (keys %{$maxProtMatch{$analysisID}}) {
		#my $refBestScore=($identifier=~/DECOY_/)? \%bestScorebad : \%bestScore; # no longer decoy in %maxProtMatch (PP 15/11/16)
		foreach my $actualSequence (keys %{$maxProtMatch{$analysisID}{$identifier}}) {
			#$maxProtScore{$analysisID}{$identifier}+=$refBestScore->{$analysisID}{$actualSequence};
			$maxProtScore{$analysisID}{$identifier}+=$bestScore{$analysisID}{$actualSequence};
		}
	}

	###################################
	####>Updating ANALYSIS_PROTEIN<####
	###################################
	####>Computing PEP_SPECIFICITY
	my %pepSpecificity;#$pepProteins{$actualSequence}{$identifier}
	foreach my $pepSeq (keys %{$pepProteins{$analysisID}}) {
		my $specificity=sprintf "%.2f",100/(scalar keys %{$pepProteins{$analysisID}{$pepSeq}});
		$specificity*=1; # 100.00 -> 100
		foreach my $identifier (keys %{$pepProteins{$analysisID}{$pepSeq}}) {
			$pepSpecificity{$identifier}=$specificity if (!$pepSpecificity{$identifier} || $pepSpecificity{$identifier}<$specificity);
		}
	}
	print '.';

	####>Finding match groups
	print '/.';
	my ($refmatchGroup,$refvisibility)=($matchGroupType eq 'MaxQuant')? (\%MQmatchGroup,\%MQvisibility) : &createMatchGroups($analysisID,$protVisibility,\%bestProtVis,\%maxProtScore,\%maxProtMatch,\%proteinRank);

	####>Storing data in DB
	print '/';
	my $newProtein=0;
	$counter=0;
	my (%boundaryStatus,%maxCovLength,%seqEndMatched);
	foreach my $identifier (keys %{$maxProtMatch{$analysisID}}) {
#next if $identifier =~ /DECOY_/;
		$counter++;
		my $des=&promsMod::resize($protDes{$identifier},250); # max length allowed in table
		my $organism=&promsMod::resize($protOrg{$identifier},100); # max length allowed in table
		my $score=($maxProtScore{$analysisID}{$identifier})?$maxProtScore{$analysisID}{$identifier}:0;
		if ($newProteins{$identifier}) { # protein is new to project=> update values !!!
			my $alias=($projectIdentMapID==$giIdentID && $identifier=~/(gi\|\d+)/)? $1 : $identifier; # "GI" restriction on alias
			$sthInsProt->execute($identifier,$alias,$des,$organism,$protMW{$identifier},$sequenceList{$identifier},$protLength{$identifier}) || die $sthInsProt->errstr;
my $proteinID=$dbh->last_insert_id(undef,undef,'PROTEIN','ID_PROTEIN');
$projectProt{$identifier}=$proteinID; # update list of project prot IDs
#foreach my $anaID (keys %{$matchList{$identifier}}) { # update for all matching analyses
#	$matchList{$identifier}{$anaID}=$proteinID;
#}
			$newProtein=1;
			delete $newProteins{$identifier};
		}
		###
		###
		###
		#my $numMatch=0;
		foreach my $actualSequence (keys %{$maxProtMatch{$analysisID}{$identifier}}){
			#my %usedBeg; # moved locally to sequence to allow '0' value when no beg info (PP 10/11/16)
			###my ($beg,$flAA_Nter,$flAA_Cter)=split(/,/,$matchInfos{$analysisID}{$identifier}{$actualSequence});
			####next if $usedBeg{$beg};
			###my $end;
			my $isSpecific=($specificSequences{$analysisID}{$actualSeq2seq{$actualSequence}})? 1 : undef;
			foreach my $beg (sort{$a<=>$b} keys %{$numMatches{$analysisID}{$identifier}{$actualSeq2seq{$actualSequence}}}) {
				my ($flAA_Nter,$flAA_Cter,$end)=@{$numMatches{$analysisID}{$identifier}{$actualSeq2seq{$actualSequence}}{$beg}};
				foreach my $peptideID (keys %{$pepIDs{$analysisID}{$actualSequence}}){
					#my ($pepSeq,$validStatus)=@{$pepIDs{$analysisID}{$actualSequence}{$peptideID}};
					my $validStatus=$pepIDs{$analysisID}{$actualSequence}{$peptideID};
					#$end=$beg+length($pepSeq)-1; # Length recorded before $pepIDs{$analysisID}{$actualSequence}{$peptideID}
					#my $isSpecific=($specificSequences{$analysisID}{$pepSeq})? 1 : undef;
					if ($validStatus==0) { # for Ghost Peptides,
						$sthAtt->execute($peptideID,$projectProt{$identifier},$analysisID,-abs($beg),-abs($end),$flAA_Nter.$flAA_Cter,$isSpecific) || die $sthAtt->errstr;
					}
					else {
						$sthAtt->execute($peptideID,$projectProt{$identifier},$analysisID,$beg,$end,$flAA_Nter.$flAA_Cter,$isSpecific) || die $sthAtt->errstr;
					}
					@{$ppa{$peptideID}{$projectProt{$identifier}}}=($beg,$end,"$flAA_Nter$flAA_Cter");
				}
				$boundaryStatus{$analysisID}{$identifier}{$beg}++;
				$boundaryStatus{$analysisID}{$identifier}{$end}--;
				$maxCovLength{$analysisID}{$identifier}=$end if (!$maxCovLength{$analysisID}{$identifier} || $maxCovLength{$analysisID}{$identifier} < $end);
				$flAA_Nter='' unless $flAA_Nter;
				$flAA_Cter='' unless $flAA_Cter;
				$seqEndMatched{$analysisID}{$identifier}=1 if $flAA_Cter eq '-'; # peptide matches end of sequence
				#$numMatch{$analysisID}{$identifier}++;
			}
			#$numMatch+=scalar keys %{$matchInfos{$analysisID}{$identifier}{$actualSequence}};
		}
		###>Computing number of peptides & matches
		my ($numPep,$numMatch)=(0,0);
		my %usedPepSeq;
		foreach my $actualSequence (keys %{$maxProtMatch{$analysisID}{$identifier}}) {
			next unless $bestScore{$analysisID}{$actualSequence}; # =0 if only ghost
			$numPep++;
			my $pepSeq=$actualSeq2seq{$actualSequence};
			next if $usedPepSeq{$pepSeq}; # Actual number of matches on protein. Can be < numPep because of PTM...
			$numMatch+=scalar keys %{$numMatches{$analysisID}{$identifier}{$pepSeq}};
			$usedPepSeq{$pepSeq}=1;
		}

		###>Computing PEP_COVERAGE
		my $coverage=0;
		my $hasPeptide=0;
		my $boundaryNter=0;
		my $pepCoverage;
		if ($protLength{$identifier}) {
			foreach my $boundary (sort{$a<=>$b} keys %{$boundaryStatus{$analysisID}{$identifier}}) {
				if ($hasPeptide==0) { # start of peptide region (cannot become <0!)
					$boundaryNter=$boundary;
				}
				$hasPeptide+=$boundaryStatus{$analysisID}{$identifier}{$boundary};
				if ($hasPeptide==0) { # end of peptide region (should be true for last boundary too)
					$coverage+=($boundary-$boundaryNter)+1;
				}
			}
			my $usedProtLength=($maxCovLength{$analysisID}{$identifier} <= $protLength{$identifier})? $protLength{$identifier} : ($seqEndMatched{$analysisID}{$identifier})? $maxCovLength{$analysisID}{$identifier} : -1*$maxCovLength{$analysisID}{$identifier}; # -: flag for protLength problem
			$pepCoverage=sprintf "%.1f",(100*$coverage)/$usedProtLength;
			$pepCoverage*=1; # 25.0 -> 25
		}

		###
		###
		###
		#my $numPep=scalar keys %{$maxProtMatch{$analysisID}{$identifier}};
		#print "'$identifier' NUMPEP=$numPep $protDbRank{$identifier},$protMW{$identifier},$protLength{$identifier},$des,$organism,$score,$maxProtScore{$analysisID}{$identifier},$refmatchGroup->{$identifier},$projectProt{$identifier}<BR>\n";
		my $confLevel = ($score)? 2 : 0;
		#print "$analysisID,$projectProt{$identifier},$score,$protDbRank{$identifier},$confLevel,",scalar keys %{$matchInfos{$analysisID}{$identifier}},",",scalar keys %{$maxProtMatch{$analysisID}{$identifier}},",$pepCoverage,$refmatchGroup->{$identifier},$pepSpecificity{$identifier},$refvisibility->{$identifier},<BR>\n";
		$insAP->execute($analysisID,$projectProt{$identifier},$score,$confLevel,$protDbRank{$identifier},$numPep,$numMatch,$pepCoverage,$refmatchGroup->{$identifier},$pepSpecificity{$identifier},$refvisibility->{$identifier}) || die $insAP->errstr();
		$protTopMG{$identifier}=1 if (defined($refvisibility->{$identifier}) && $refvisibility->{$identifier}==2);
		if ($counter > 300) {
			print ".";
			$counter=0;
		}

	}

	push @newAnaMapping,$analysisID if $newProtein;

	####>Removing protecting IDs
	##$dbh->do("DELETE FROM PROTEIN WHERE ID_PROTEIN=$maxProteinID") || $dbh->errstr;
	##$dbh->commit;

	print " Done.<BR>\n";

}
$sthInsProt->finish;
$insAP->finish;
$sthTop->finish;
#$dbh->do("DELETE FROM PROTEIN WHERE ID_PROTEIN=$protectProtID") || $dbh->errstr; # in case of error in previous process (rare!)
$dbh->commit;
print "<FONT class=\"title3\">&nbsp;Done.</FONT><BR>\n";

#$dbh->disconnect; exit; # DEBUG


#######################################
####> Peptide quantification data <####
#######################################
#$importPepQuantif=0; # DEBUG
if ($importPepQuantif) {

print "<BR><FONT class='title3'>Recording peptide XIC data:</FONT><BR>\n";

#my $sthInsPepQ=$dbh->prepare("INSERT INTO PEPTIDE_QUANTIFICATION (ID_QUANTIFICATION,ID_QUANTIF_PARAMETER,ID_PEPTIDE,QUANTIF_VALUE,TARGET_POS) VALUES (?,?,?,?,?)");
mkdir "$promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
my %pepQuantHandle;
foreach my $pepQuantifID (values %ana2pepQuantification) {
	my $quantifDir="$promsPath{quantification}/project_$projectID/quanti_$pepQuantifID";
	mkdir $quantifDir;
	if ($labeling eq 'FREE') {
		open($pepQuantHandle{$pepQuantifID}{0},">$quantifDir/peptide_quantification.txt") || die $!;
		print {$pepQuantHandle{$pepQuantifID}{0}} "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n";
	}
	else {
		foreach my $colLabel (@labelList) {
			my $targetPos=$label2targetPos{$colLabel};
			open($pepQuantHandle{$pepQuantifID}{$targetPos},">$quantifDir/peptide_quantification_$targetPos.txt") || die $!;
			print {$pepQuantHandle{$pepQuantifID}{$targetPos}} "ID_QUANTIF_PARAMETER\tID_PEPTIDE\tQUANTIF_VALUE\n";
		}
	}
}

$anaCounter=0;
###> SILAC<###
if ($labeling eq 'SILAC') {
	my $sthAddVP=$dbh->prepare("INSERT INTO PEPTIDE (ID_ANALYSIS,PEP_SEQ,PEP_LENGTH,MISS_CUT,CHARGE,ELUTION_TIME,DATA,VALID_STATUS,QUERY_NUM,PEP_RANK) VALUES (?,?,?,?,?,?,?,0,0,0)");
	my $sthAddVPM=$dbh->prepare("INSERT INTO PEPTIDE_MODIFICATION (ID_PEPTIDE,ID_MODIFICATION,MODIF_TYPE,POS_STRING,REF_POS_STRING) VALUES (?,?,'V',?,?)");
	my $sthUpPep=$dbh->prepare("UPDATE PEPTIDE SET DATA=CONCAT(DATA,?) WHERE ID_PEPTIDE=?");
	foreach my $analysisID (sort{$a<=>$b} keys %peptideXIC) {
		$anaCounter++;
		$counter=0;
		print "&nbsp;-$anaNames{$analysisID} ($anaCounter/$totalAna)...";
		my $qsetNum=0;
		foreach my $queryNum (sort{$a<=>$b} keys %{$peptideXIC{$analysisID}}) {
			next unless $peptideXIC{$analysisID}{$queryNum}{'XIC'}; # no quantif data all
			$qsetNum++;
			foreach my $colLabel (@labelList) {
				$counter++;
				if ($counter > 2000) {
					print '.';
					$counter=0;
				}
				my $targetPos=$label2targetPos{$colLabel};
				if ($peptideXIC{$analysisID}{$queryNum}{'XIC'}{$colLabel}) {
					#if ($labelInfos{$analysisID}{$colLabel}{'MODIFICATION'}) { #}
					if ($isotopeLabelInfo{$colLabel}) {
						my $pepID=$peptideXIC{$analysisID}{$queryNum}{'ID_PEPTIDE'};
						my ($actualSequence,@data)=@{$pepInfo{$analysisID}{$queryNum}};# @data=(pepSeq,Length,$anaID,$score,$missCut,$massExp,$Mass,$mz,$massError,$charge,$et,$vamidStatus,$data,$specCount)
						my $pepSeq=$data[0];
						###> Create a Ghost Peptide
						my $isSpecific=($specificSequences{$analysisID}{$pepSeq})? 1 : undef;
						$sthAddVP->execute($analysisID,$pepSeq,$data[1],$data[4],$data[9],$data[10],"QSET=$qsetNum");
						my $vpPepID=$dbh->last_insert_id(undef,undef,'PEPTIDE','ID_PEPTIDE');
						###> Add modifications
						my $posString='';
						###> Label ones
						foreach my $modID (keys %{$isotopeLabelInfo{$colLabel}}) {
							if ($isotopeLabelInfo{$colLabel}{$modID}{'*'}) {
								$posString='*';
							}
							elsif ($isotopeLabelInfo{$colLabel}{$modID}{'='}) {
								$posString='=';
							}
							else {
								my @aas=split(//,$pepSeq);
								my @pos;
								for (my $i = 0 ; $i <=$#aas ; $i++) {
									push @pos,$i+1 if $isotopeLabelInfo{$colLabel}{$modID}{$aas[$i]};
								}
								$posString=join('.',@pos);
							}
							$sthAddVPM->execute($vpPepID,$modID,$posString,undef);
						}
						###> Regular ones (vmod like oxidation, acetylation,...)
						foreach my $vmod (keys %{$pepVmod{$analysisID}{$queryNum}}) {
							$sthAddVPM->execute($vpPepID,$anaVmods{$analysisID}{$vmod},$seq2posString{$actualSequence}{$vmod},$pepVmod{$analysisID}{$queryNum}{$vmod}[0]); # use same site probability as validated peptide (good, bad?)
						}
						###> Create a link to PEPTIDE_PROTEIN_ATTRIB with this ghost peptide
						foreach my $protID (keys %{$ppa{$peptideXIC{$analysisID}{$queryNum}{'ID_PEPTIDE'}}}) {
							my ($beg,$end,$matchInfo)=@{$ppa{$peptideXIC{$analysisID}{$queryNum}{'ID_PEPTIDE'}}{$protID}};
							$sthAtt->execute($vpPepID,$protID,$analysisID,-abs($beg),-abs($end),$matchInfo,$isSpecific);
						}
						###> Add quantification for this specific Ghost peptide
#$sthInsPepQ->execute($ana2pepQuantification{$analysisID},$areaParamID,$vpPepID,$peptideXIC{$analysisID}{$queryNum}{'XIC'}{$colLabel},$targetPos);
print {$pepQuantHandle{$ana2pepQuantification{$analysisID}}{$targetPos}} "$areaParamID\t$vpPepID\t$peptideXIC{$analysisID}{$queryNum}{XIC}{$colLabel}\n";
						#print "$quantiSILAC{$analysisID},$areaParamID,$vpPepID,$peptideXIC{$analysisID}{$queryNum}{\"XIC $colLabel\"},$labelInfos{$analysisID}{$colLabel}{'TARGET_POS'}<BR>\n";
					}
					else { # No modification ID was found so this is Light label
						#print "$quantiSILAC{$analysisID},$areaParamID,$peptideXIC{$analysisID}{$queryNum}{'ID_PEPTIDE'},$peptideXIC{$analysisID}{$queryNum}{\"XIC $colLabel\"},$labelInfos{$analysisID}{$colLabel}{'TARGET_POS'}<BR>";
#$sthInsPepQ->execute($ana2pepQuantification{$analysisID},$areaParamID,$peptideXIC{$analysisID}{$queryNum}{'ID_PEPTIDE'},$peptideXIC{$analysisID}{$queryNum}{'XIC'}{$colLabel},$targetPos);
print {$pepQuantHandle{$ana2pepQuantification{$analysisID}}{$targetPos}} "$areaParamID\t$peptideXIC{$analysisID}{$queryNum}{ID_PEPTIDE}\t$peptideXIC{$analysisID}{$queryNum}{XIC}{$colLabel}\n";
						#my ($data)=$dbh->selectrow_array("SELECT DATA FROM PEPTIDE WHERE ID_PEPTIDE=$peptideXIC{$analysisID}{$queryNum}{'ID_PEPTIDE'}");
						#$data.="##QSET=$qsetNum";
						$sthUpPep->execute("##QSET=$qsetNum",$peptideXIC{$analysisID}{$queryNum}{'ID_PEPTIDE'});
					}
				}
			}
		}
		$dbh->commit;
		print " Done.<BR>\n";
	}
	$sthAddVP->finish;
	$sthAddVPM->finish;
	$sthUpPep->finish;
}
###>iTRAQ or TMT <###
elsif ($isobarModifInfo{ID}) { # $labeling=~/ITRAQ|TMT/
	foreach my $analysisID (sort{$a<=>$b} keys %peptideXIC) {
		$anaCounter++;
		$counter=0;
		print "&nbsp;-$anaNames{$analysisID} ($anaCounter/$totalAna)...";
		foreach my $queryNum (sort{$a<=>$b} keys %{$peptideXIC{$analysisID}}) {
			next unless $peptideXIC{$analysisID}{$queryNum}{'XIC'}; # no quantif data all
			foreach my $repIdx (@labelList) { #(sort{$a<=>$b} keys %{$peptideXIC{$analysisID}{$queryNum}{XIC}})
				next unless $peptideXIC{$analysisID}{$queryNum}{'XIC'}{$repIdx}; # skip 0/undef values
				$counter++;
				if ($counter > 2000) {
					print '.';
					$counter=0;
				}
#$sthInsPepQ->execute($ana2pepQuantification{$analysisID},$areaParamID,$peptideXIC{$analysisID}{$queryNum}{'ID_PEPTIDE'},$peptideXIC{$analysisID}{$queryNum}{'XIC'}{$repIdx},$repIdx+1);
my $targetPos=$repIdx+1;
print {$pepQuantHandle{$ana2pepQuantification{$analysisID}}{$targetPos}} "$areaParamID\t$peptideXIC{$analysisID}{$queryNum}{ID_PEPTIDE}\t$peptideXIC{$analysisID}{$queryNum}{XIC}{$repIdx}\n";
			}
		}
#$dbh->commit;
		print " Done.<BR>\n";
	}
}
###>Label-free<###
elsif ($labeling eq 'FREE') {
	foreach my $analysisID (sort{$a<=>$b} keys %peptideXIC) {
		$anaCounter++;
		$counter=0;
		print "&nbsp;-$anaNames{$analysisID} ($anaCounter/$totalAna)...";
		foreach my $queryNum (sort{$a<=>$b} keys %{$peptideXIC{$analysisID}}) {
			$counter++;
			if ($counter > 2000) {
				print '.';
				$counter=0;
			}
#$sthInsPepQ->execute($ana2pepQuantification{$analysisID},$areaParamID,$peptideXIC{$analysisID}{$queryNum}{'ID_PEPTIDE'},$peptideXIC{$analysisID}{$queryNum}{'XIC'}{'NONE'},undef) if $peptideXIC{$analysisID}{$queryNum}{'XIC'};
print {$pepQuantHandle{$ana2pepQuantification{$analysisID}}{0}} "$areaParamID\t$peptideXIC{$analysisID}{$queryNum}{ID_PEPTIDE}\t$peptideXIC{$analysisID}{$queryNum}{XIC}{NONE}\n" if $peptideXIC{$analysisID}{$queryNum}{'XIC'};
		}
		$dbh->commit;
		print " Done.<BR>\n";
	}
}
#$sthInsPepQ->finish;
foreach my $pepQuantifID (keys %pepQuantHandle) {
	foreach my $targetPos (keys %{$pepQuantHandle{$pepQuantifID}}) {close $pepQuantHandle{$pepQuantifID}{$targetPos};}
}

print "<FONT class=\"title3\">&nbsp;Done.</FONT><BR>\n";

} # END of $importPepQuantif

$sthAtt->finish;


#######################################
####> Protein quantification data <####
#######################################
if ($importProtQuantif) {

	print "<BR><FONT class='title3'>Importing protein quantification data from proteinGroups file...";
	my $pepTotStg=($pepUsed eq 'razor')? 'Peptide counts (razor+unique)' : ($pepUsed eq 'unique')? 'Peptide counts (unique)' : 'Peptide counts (all)';
	my $sthInsProtQ=$dbh->prepare("INSERT INTO PROTEIN_QUANTIFICATION (ID_PROTEIN,ID_QUANTIFICATION,ID_QUANTIF_PARAMETER,QUANTIF_VALUE,TARGET_POS) VALUES (?,?,?,?,?)");
	open (PROTEIN,"$tmpFilesDir/$proteinGroupsFile") || die "Unable to open $tmpFilesDir/$proteinGroupsFile";

	my %proteinColNum;
	$counter=0;
	while (my $line=<PROTEIN>) {
		chomp($line);
		my @parameters=split(/ *\t */,$line); # remove starting/trailing spaces
		if ($.==1) { # 1st line of the file
			my $ncol=0;
			foreach my $colName (@parameters) {
				$proteinColNum{$colName}=$ncol;
				$ncol++;
			}
			next;
		}
		###> In proteinGroups.txt, you find in the column the distribution of peptides found for each entry
		###> The idea is to keep only quantification values made for the same set of peptides.
		###> ex : Protein IDs column => O43790;CON__O43790;CON__Q14533;Q14533
		###>      Peptide counts (all) column => 3;3;3;2
		###>      Then, only protein O43790, CON__O43790 and CON__Q14533 will be kept in MQ quantitation
		###> PP 24/11/16: Restrict rule above to proteins at top of a match group in at least 1 analysis
		my @pepTotal=split(/;/,$parameters[$proteinColNum{$pepTotStg}]);
		my @identifiers; #=split(/;/,$parameters[$proteinColNum{'Majority protein IDs'}]); # 'Protein IDs'
		foreach my $rawIdentifier (split(/;/,$parameters[$proteinColNum{'Majority protein IDs'}])) {
			push @identifiers,$raw2UsedIdentifiers{$rawIdentifier} || $rawIdentifier;
		}
		next unless $projectProt{$identifiers[0]}; # Skip "dirty" lines... (PP 10/11/16)
		my %selProtFromMQgroup;
		$selProtFromMQgroup{$projectProt{$identifiers[0]}}=$pepTotal[0] if $protTopMG{$identifiers[0]}; # must be top of a MG in >= 1 MQ ana (bestVis=2) Should be always true for this identifier # if $projectProt{$identifiers[0]};
		for (my $i=1; $i<=$#identifiers; $i++){
			last if $pepTotal[$i] < $pepTotal[0]; # in case of very large groups, peptide counts are not recorded for all member
			$selProtFromMQgroup{$projectProt{$identifiers[$i]}}=$pepTotal[$i] if $protTopMG{$identifiers[$i]}; # must be top of a MG in >= 1 MQ ana (bestVis=2)
		}
		foreach my $paramText (keys %recordedParams) {  # %{$designInfo{'Quantification'}}
			next unless $proteinColNum{$paramText};
			my $ncol = $proteinColNum{$paramText};
next unless $ncol>=0; # ?
			$parameters[$ncol]=~s/;.+//; # Global 'Peptide counts xxx' column ('X;Y;Z;..') is used when no design. Only 1st value is kept (PP 08/11/16)
			next if (!$parameters[$ncol] || $parameters[$ncol] eq 'NaN' || $parameters[$ncol] == 0); # avoid to store 0 values in myProMS db
			foreach my $refQuantif (@{$recordedParams{$paramText}}) {
				my ($qID,$qCode,$targetPos)=@{$refQuantif}; # @{$designInfo{'Quantification'}{$colName}};
				next unless $quantifParamIDs{$qCode};
				my $paramOK=0;
				if ($targetPos==0) { # required quantif value of at least 1 channel must be OK
					foreach my $tgPos (keys %{$requiredParams{$qID}}) {
						if ($parameters[$proteinColNum{$requiredParams{$qID}{$tgPos}}] && $parameters[$proteinColNum{$requiredParams{$qID}{$tgPos}}] ne 'NaN') {
							$paramOK=1;
							last;
						}
					}
					next unless $paramOK;
				}
				elsif ($parameters[$proteinColNum{$requiredParams{$qID}{$targetPos}}] && $parameters[$proteinColNum{$requiredParams{$qID}{$targetPos}}] ne 'NaN') { # do not record any values if required one fails
					$paramOK=1;
				}
				next unless $paramOK;
				my $usedTargetPos=$targetPos || undef; # 0 -> undef
				foreach my $protID (keys %selProtFromMQgroup) {
					$sthInsProtQ->execute($protID,$qID,$quantifParamIDs{$qCode},$parameters[$ncol],$usedTargetPos);
				}
			}
		}

		$counter++;
		if ($counter > 250) {
			print '.';
			$counter=0;
		}
	}
	close PROTEIN;
	$sthInsProtQ->finish;
	$dbh->commit;
	print " Done.</FONT><BR>\n";

}

$dbh->disconnect;

###>Move all files in peptide quantitation folder & link to protein quantification folder(s) !
mkdir "$promsPath{quantification}/project_$projectID" unless -e "$promsPath{quantification}/project_$projectID";
foreach my $quantifID ($peptideQuantifID,@protQuantifList) {
	mkdir "$promsPath{quantification}/project_$projectID/quanti_$quantifID"; # unless -e "$promsPath{quantification}/project_$projectID/quanti_$quantifID";
	foreach my $file (@requiredFiles) {
		if ($quantifID==$peptideQuantifID) {
			move("$tmpFilesDir/$file","$promsPath{quantification}/project_$projectID/quanti_$quantifID/$file");
#copy("$tmpFilesDir/$file","$promsPath{quantification}/project_$projectID/quanti_$quantifID/$file"); # DEBUG
		}
		else {
			symlink("$promsPath{quantification}/project_$projectID/quanti_$peptideQuantifID/$file","$promsPath{quantification}/project_$projectID/quanti_$quantifID/$file");
		}
	}
	### OR symlink whole directory (not safe for deletion)
}
sleep 2;
rmtree $tmpFilesDir;

#print '**DONE**'; exit; #DEBUG !!!

################################
####>New identifier Mapping<####
################################
if (scalar @newAnaMapping) { # true if new valid proteins in project
	my $anaStrg=join(',',@newAnaMapping);
	system "./mapProteinIdentifiers.pl $userID $anaStrg";
}
print "<BR><BR><FONT class=\"title2\">Import was successful! Updating display...</FONT><BR><BR>\n";
sleep 2;
#exit; # DEBUG
print qq
|<SCRIPT type="text/javascript">
top.promsFrame.selectedAction='summary';
top.experimentView='quanti';
parent.navFrame.location="$promsPath{cgi}/openProject.cgi?ID=$projectID&branchID=experinment:$experimentID&ACT=nav";
</SCRIPT>
</BODY>
</HTML>
|;

################################
####<Generates Match Groups>####
################################
sub createMatchGroups {
	my ($analysisID,$protVisibility,$refBestProtVis,$refMaxScore,$refMaxProtMatch,$refProteinRank)=@_;
	my (%matchGroup,%visibility);
	my $numGroup=0;
#print "$analysisID\n";
	my @sortedIdentifiers=(sort{scalar (keys %{$refMaxProtMatch->{$analysisID}{$b}})<=>scalar (keys %{$refMaxProtMatch->{$analysisID}{$a}}) || $refMaxScore->{$analysisID}{$b}<=>$refMaxScore->{$analysisID}{$a} || $refProteinRank->{$a}<=>$refProteinRank->{$b} || $a cmp $b} keys %{$refMaxProtMatch->{$analysisID}});
	my $counter=0;
	foreach my $i (0..$#sortedIdentifiers) {
		next if $matchGroup{$sortedIdentifiers[$i]}; # already assigned to a match group
		$matchGroup{$sortedIdentifiers[$i]}=++$numGroup;
        $visibility{$sortedIdentifiers[$i]}=2; # Reference protein
		$refBestProtVis->{$sortedIdentifiers[$i]}=2; # update bestVis
#next; # SKIP grouping!!!
		foreach my $j ($i+1..$#sortedIdentifiers) {
           next if $matchGroup{$sortedIdentifiers[$j]}; # already assigned to a match group
			##<Comparing peptide contents of identifier#1 and identifier#2>## All peptides must match!
			my $matchOK=1;
			foreach my $seq (keys %{$refMaxProtMatch->{$analysisID}{$sortedIdentifiers[$j]}}) {
				$counter++;
				if ($counter >= 100000) {
					print '.';
					$counter=0;
				}
 				if (!$refMaxProtMatch->{$analysisID}{$sortedIdentifiers[$i]}{$seq}) {
					delete $refMaxProtMatch->{$analysisID}{$sortedIdentifiers[$i]}{$seq}; # to be safe
					$matchOK=0;
					last;
				}
			}
			if ($matchOK) {
				$matchGroup{$sortedIdentifiers[$j]}=$matchGroup{$sortedIdentifiers[$i]};
				$visibility{$sortedIdentifiers[$j]}=($protVisibility && defined($refBestProtVis->{$sortedIdentifiers[$j]}) && ($refBestProtVis->{$sortedIdentifiers[$j]}==2 || ($protVisibility==2 && $refBestProtVis->{$sortedIdentifiers[$j]})))? 1 : 0;
				$refBestProtVis->{$sortedIdentifiers[$j]}=$visibility{$sortedIdentifiers[$j]} if (!defined($refBestProtVis->{$sortedIdentifiers[$j]}) || $visibility{$sortedIdentifiers[$j]} > $refBestProtVis->{$sortedIdentifiers[$j]}); # update bestVis
			}
		}
	}
	return(\%matchGroup,\%visibility);
}
sub createMaxQuantMatchGroups {
	my ($protVisibility,$refMQmatchGroup,$refMQvisibility,$refmatchList,$refBestProtVis,$refRaw2UsedIdentifiers)=@_;
	print "&nbsp;-Applying MaxQuant match group rules...";
	open (PROTEIN,"$tmpFilesDir/$proteinGroupsFile") || die "Unable to open $tmpFilesDir/$proteinGroupsFile";
	my $counter=0;
	my $mgNum=0;
	my ($allProtColumnIdx,$topProtColumnIdx)=(0,1);
	while (my $line=<PROTEIN>) {
		my @parameters=split(/ *\t */,$line); # remove starting/trailing spaces
		if ($.==1) { # 1st line of the file
			my $ncol=0;
			foreach my $colName (@parameters) {
				if ($colName eq 'Proteins IDs') {$allProtColumnIdx=$ncol;}
				elsif ($colName eq 'Majority protein IDs') {$topProtColumnIdx=$ncol;}
				$ncol++;
			}
			next;
		}
		#my @topIdentifiers=split(/;/,$parameters[$topProtColumnIdx]);
my $topIdentifier;
		foreach my $rawIdentifier (split(/;/,$parameters[$topProtColumnIdx])) {
			next if $rawIdentifier =~/REV__/;
			#$rawIdentifier =~ s/REV__/DECOY_/; # shouldn't be useful since DECOY are not in %{$refmatchList}
			my $identifier=$refRaw2UsedIdentifiers->{$rawIdentifier} || $rawIdentifier;
if ($refmatchList->{$identifier}) {
	$topIdentifier=$identifier;
	last;
}
		}
		next unless $topIdentifier; # Skip "dirty" lines... (PP 10/11/16)
		#next unless $refmatchList->{$topIdentifiers[0]}; # Skip "dirty" lines... (PP 10/11/16)
		$mgNum++;
		foreach my $rawIdentifier (split(/;/,$parameters[$allProtColumnIdx])) {
			my $identifier=$refRaw2UsedIdentifiers->{$rawIdentifier} || $rawIdentifier;
			$refMQmatchGroup->{$identifier}=$mgNum;
			$refMQvisibility->{$identifier}=($identifier eq $topIdentifier)? 2 : ($protVisibility && defined($refBestProtVis->{$identifier}) && ($refBestProtVis->{$identifier}==2 || ($protVisibility==2 && $refBestProtVis->{$identifier})))? 1 : 0;
		}
		$counter++;
		if ($counter >= 250) {
			print '.';
			$counter=0;
		}
	}
	close PROTEIN;
	print " Done.<BR>\n";
}

###########################################################
####<Table of Modifications for report ion in MaxQuant>####
###########################################################
sub getModIDforReporterIon { # Modification with corresponding Unimod ID as checked on Unimod website 09 of November 2015
	my %mqRepIon=();
	$mqRepIon{'18O'}{'UNIMOD_ACC'}=193;
	$mqRepIon{'Arg10'}{'UNIMOD_ACC'}=267;
	$mqRepIon{'Arg6'}{'UNIMOD_ACC'}=188;
	$mqRepIon{'DimethLys0'}{'UNIMOD_ACC'}=36;
	$mqRepIon{'DimethLys2'}{'UNIMOD_ACC'}=199;
	$mqRepIon{'DimethLys4'}{'UNIMOD_ACC'}=510;
	$mqRepIon{'DimethLys6'}{'UNIMOD_ACC'}=1291;
	$mqRepIon{'DimethLys8'}{'UNIMOD_ACC'}=330;
	$mqRepIon{'DimethNter0'}{'UNIMOD_ACC'}=36;
	$mqRepIon{'DimethNter2'}{'UNIMOD_ACC'}=199;
	$mqRepIon{'DimethNter4'}{'UNIMOD_ACC'}=510;
	$mqRepIon{'DimethNter6'}{'UNIMOD_ACC'}=1291;
	$mqRepIon{'DimethNter8'}{'UNIMOD_ACC'}=330;
	$mqRepIon{'ICAT-0'}{'UNIMOD_ACC'}=105;
	$mqRepIon{'ICAT-9'}{'UNIMOD_ACC'}=106;
	$mqRepIon{'ICPL-Lys0'}{'UNIMOD_ACC'}=365;
	$mqRepIon{'ICPL-Lys10'}{'UNIMOD_ACC'}=687;
	$mqRepIon{'ICPL-Lys4'}{'UNIMOD_ACC'}=364;
	$mqRepIon{'ICPL-Lys6'}{'UNIMOD_ACC'}=866;
	$mqRepIon{'ICPL-Nter0'}{'UNIMOD_ACC'}=365;
	$mqRepIon{'ICPL-Nter10'}{'UNIMOD_ACC'}=687;
	$mqRepIon{'ICPL-Nter4'}{'UNIMOD_ACC'}=364;
	$mqRepIon{'ICPL-Nter6'}{'UNIMOD_ACC'}=866;
	$mqRepIon{'Ile7'}{'UNIMOD_ACC'}=695;
	$mqRepIon{'Leu7'}{'UNIMOD_ACC'}=695;
	$mqRepIon{'Lys4'}{'UNIMOD_ACC'}=481;
	$mqRepIon{'Lys6'}{'UNIMOD_ACC'}=188;
	$mqRepIon{'Lys8'}{'UNIMOD_ACC'}=259;
	$mqRepIon{'mTRAQ-Lys0'}{'UNIMOD_ACC'}=888;
	$mqRepIon{'mTRAQ-Lys4'}{'UNIMOD_ACC'}=889;
	$mqRepIon{'mTRAQ-Lys8'}{'UNIMOD_ACC'}=1302;
	$mqRepIon{'mTRAQ-Nter0'}{'UNIMOD_ACC'}=888;
	$mqRepIon{'mTRAQ-Nter4'}{'UNIMOD_ACC'}=889;
	$mqRepIon{'mTRAQ-Nter8'}{'UNIMOD_ACC'}=1302;

	$mqRepIon{'18O'}{'INTERIM_NAME'}='double_O18';
	$mqRepIon{'Arg10'}{'INTERIM_NAME'}='13C6-15N4';
	$mqRepIon{'Arg6'}{'INTERIM_NAME'}='13C6';
	$mqRepIon{'DimethLys0'}{'INTERIM_NAME'}='di-Methylation ';
	$mqRepIon{'DimethLys2'}{'INTERIM_NAME'}='CHD2';
	$mqRepIon{'DimethLys4'}{'INTERIM_NAME'}='C13HD2';
	$mqRepIon{'DimethLys6'}{'INTERIM_NAME'}='Dimethyl:2H(6)';
	$mqRepIon{'DimethLys8'}{'INTERIM_NAME'}='Dimethyl:2H(6)13C(2)';
	$mqRepIon{'DimethNter0'}{'INTERIM_NAME'}='di-Methylation';
	$mqRepIon{'DimethNter2'}{'INTERIM_NAME'}='CHD2';
	$mqRepIon{'DimethNter4'}{'INTERIM_NAME'}='C13HD2';
	$mqRepIon{'DimethNter6'}{'INTERIM_NAME'}='Dimethyl:2H(6)';
	$mqRepIon{'DimethNter8'}{'INTERIM_NAME'}='Dimethyl:2H(6)13C(2)';
	$mqRepIon{'ICAT-0'}{'INTERIM_NAME'}='ICAT_light';
	$mqRepIon{'ICAT-9'}{'INTERIM_NAME'}='ICAT_heavy';
	$mqRepIon{'ICPL-Lys0'}{'INTERIM_NAME'}='ICPL_light';
	$mqRepIon{'ICPL-Lys10'}{'INTERIM_NAME'}='ICPL_medium';
	$mqRepIon{'ICPL-Lys4'}{'INTERIM_NAME'}='ICPL_heavy';
	$mqRepIon{'ICPL-Lys6'}{'INTERIM_NAME'}='ICPL:13C(6)2H(4)';
	$mqRepIon{'ICPL-Nter0'}{'INTERIM_NAME'}='ICPL_light';
	$mqRepIon{'ICPL-Nter10'}{'INTERIM_NAME'}='ICPL_medium';
	$mqRepIon{'ICPL-Nter4'}{'INTERIM_NAME'}='ICPL_heavy';
	$mqRepIon{'ICPL-Nter6'}{'INTERIM_NAME'}='ICPL:13C(6)2H(4)';
	$mqRepIon{'Ile7'}{'INTERIM_NAME'}='Label:13C(6)15N(1)';
	$mqRepIon{'Leu7'}{'INTERIM_NAME'}='Label:13C(6)15N(1)';
	$mqRepIon{'Lys4'}{'INTERIM_NAME'}='Lys4';
	$mqRepIon{'Lys6'}{'INTERIM_NAME'}='13C6';
	$mqRepIon{'Lys8'}{'INTERIM_NAME'}='13C6-15N2 ';
	$mqRepIon{'mTRAQ-Lys0'}{'INTERIM_NAME'}='mTRAQ';
	$mqRepIon{'mTRAQ-Lys4'}{'INTERIM_NAME'}='mTRAQ:13C(3)15N(1)';
	$mqRepIon{'mTRAQ-Lys8'}{'INTERIM_NAME'}='mTRAQ:13C(6)15N(2)';
	$mqRepIon{'mTRAQ-Nter0'}{'INTERIM_NAME'}='mTRAQ';
	$mqRepIon{'mTRAQ-Nter4'}{'INTERIM_NAME'}='mTRAQ:13C(3)15N(1)';
	$mqRepIon{'mTRAQ-Nter8'}{'INTERIM_NAME'}='mTRAQ:13C(6)15N(2)';
	return(%mqRepIon);
}


####>Revision history<####
# 1.3.3 Minor change to revert code block for modifications allowed for quantification to v1.3.1 (PP 01/10/18)
# 1.3.2 Add extra files for maxquant import other than TMT or iTRAQ: summary.txt and parameters.txt (GA 09/08/18)
# 1.3.1 [Fix] bug peptide specificity computation now excludes decoy and optionally contaminent matches (PP 15/06/18)
# 1.3.0 Peptide quantification data now written to file $promsPath{quantification}/project_$projectID/quanti_$quantifID/peptide_quantification(_$targetPos).txt (PP 09/05/18)
# 1.2.1 Minor modification in myProMS import (GA 10/04/18)<BR> TODO : test the import of many groups
# 1.2.0 Updated for auto-increment in tables ANALYSIS,PROTEIN,QUANTIFICATION & added compatibility with shared directory (PP 23/03/18)
# 1.1.5 Modification in <EVIDENCE> parsing because of GlyGly modification (GA 06/12/17)
# 1.1.4 Added MaxQuant version string in QUANTIF_ANNOT field (PP 13/09/17)
# 1.1.3 Matches Analysis display pos with fraction number if used (PP 09/08/17)
# 1.1.2 Bug fix in peptide start position for non-razor proteins & records identification evidence id in PEPTIDE.DATA (PP 08/03/17)
# 1.1.1 Records modification position probability & bug fix in SILAC ghost peptides & optional import of protein quantifications (PP 17/02/17)
# 1.1.0 Major revision in import mechanism and data storage design. Not compatible with older imports (PP 20/01/17)
# 1.0.6 Add UPDATE_DATE and UPDATE_USER in INSERT queries for ANALYSIS and QUANTIFICATION tables (GA 19/10/2016)
# 1.0.5 Change to deal with autoincrement in PEPTIDE table (GA 19/10/2016)
# 1.0.4 Minor modification (GA 05/10/2016)
# 1.0.3 No import for quantitative data if values are articially put to zero (GA 12/09/2016)
# 1.0.2 Minor change for MaxQuant query. MAXQUANT code becomes MQ (GA 27/07/2016)
# 1.0.1 Add labeled import (SILAC) & update design information (GA 09/05/2016)<BR> TODO : tests with iTRAQ data
# 1.0.0 Import for MaxQuant quantitated folder<BR> TODO: check for replicate<->fraction depedencies and update Observation according to technical replicates (GA 24/09/15)<BR>TODO Import as well Light and Heavy information
