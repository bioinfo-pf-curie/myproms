#!/usr/local/bin/perl -w

################################################################################
# drawSpectrum.cgi                  1.8.6                                      #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Checks if a reference spectrum exists                                        #
# before launching peptide_view.cgi (1 or 2 frames)                            #
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
use promsConfig;
use promsMod;
use strict;

# print header; warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};

####################
####>Parameters<####
####################
my ($objectID,$queryNum,$rank)=(split(/_/,param('RID')))[1..3]; # objectID=queryID or pepID
my $call=(param('CALL'))? param('CALL') : 'rank';
#my $identifier=param('IDENT'); # required by on-line Mascot

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

###############################################
####>Fetching Peptide Sequence + Data File<####
###############################################
my ($projectID,$analysisID,$dataFile,$sequence,$varMod,$massObs,$comments,$msType);
my $fileFormat;
if ($call eq 'pep' ||  $call eq 'val') { # called from sequenceView or validated analysis called by listMatchedAnalysis
	($analysisID,$sequence,$massObs,$comments)=$dbh->selectrow_array("SELECT ID_ANALYSIS,PEP_SEQ,MR_OBS,COMMENTS FROM PEPTIDE WHERE ID_PEPTIDE=$objectID");
	(my $validStatus,$fileFormat,my $dataFileName,$msType)=$dbh->selectrow_array("SELECT VALID_STATUS,FILE_FORMAT,DATA_FILE,MS_TYPE FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
	$projectID=&promsMod::getProjectID($dbh,$analysisID,'analysis');
	$varMod=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$objectID,$analysisID,$sequence);

	#my $pepFileRoot=sprintf "P%06d",$analysisID;
	#$dataFile="$promsPath{peptide}/proj_$projectID/$pepFileRoot.";
	#$dataFile.=($fileFormat=~/\.XML/)? 'pgf' : ($fileFormat=~/\.PDM/)? 'pdm' : 'dat';
	#unless (-e $dataFile) {
	#	$dataFile="$promsPath{peptide}/proj_$projectID/$analysisID"."_$dataFileName"; # new data file conservation procedure (01/06/11)
	#	$dataFile=~s/\.xml/\.pgf/ if ($fileFormat=~/\.XML/ && $fileFormat ne 'PARAGON.XML');
	#}
	$dataFileName = "swath_ana_$analysisID.txt" if($fileFormat eq 'SPECTRONAUT.XLS');
	$dataFileName = 'msms.txt' if($fileFormat eq 'MAXQUANT.DIR');
	$dataFile=($validStatus==2)? "$promsPath{peptide}/proj_$projectID/ana_$analysisID/$dataFileName" : "$promsPath{valid}/ana_$analysisID/$dataFileName"; # new data file structure (04/03/13)
	unless (-e $dataFile) {
		if ($fileFormat eq 'MASCOT.DAT') {$dataFile=~s/\.dat/_min\.dat/;} # "old" minimal file
		elsif ($fileFormat=~/\.XML/ && $fileFormat ne 'PARAGON.XML') {$dataFile=~s/\.xml/\.pgf/;}
	}
	$comments='' unless $comments;
}
elsif ($call eq 'rank' || $call eq 'ana' || $call eq 'seq') { # call=rank or ana (called during validation) or seq (called from sequence view of validated)
	my ($fileName,$pepInfo,$massData);
	($analysisID,$fileName,$msType,$pepInfo,$massData)=$dbh->selectrow_array("SELECT Q.ID_ANALYSIS,A.DATA_FILE,A.MS_TYPE,Q.INFO_PEP$rank,Q.MASS_DATA FROM QUERY_VALIDATION Q,ANALYSIS A WHERE Q.ID_ANALYSIS=A.ID_ANALYSIS AND ID_QUERY=$objectID");
	$projectID=&promsMod::getProjectID($dbh,$analysisID,'analysis');
	$dataFile="$promsPath{valid}/ana_$analysisID/$fileName";
	($sequence)=($pepInfo=~/SEQ=(\w+)/);
	($varMod)=&promsMod::toStringVariableModifications($dbh,'QUERY',$objectID,$analysisID,$sequence,$rank);
	($massObs)=($massData=~/OBS=(\d+\.*\d*)/);
	($comments)=($pepInfo=~/COM=(.+)/);
	$comments='' unless $comments;
	$fileFormat=$dbh->selectrow_array("SELECT FILE_FORMAT FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
}
else {
	print header; warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Peptide Fragmentation Spectrum</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY><FONT class="title3">Unknown call parameters: '$call'</FONT><BR></BODY>
</HTML>
|;
	exit;
}

####>Variable modification(s) w/o position<####
if ($varMod && $varMod=~/-1/) {
	$varMod=~s/-1/<FONT color='#DD0000'>?<\/FONT>/g;
	my $commentsString="";
	if($comments){ $commentsString="<BR><BR><FONT class=\"title2\">Comments : </FONT><FONT class=\"title3\">$comments</FONT>";}
	print header; warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Peptide Fragmentation Spectrum</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY>
<CENTER>
<BR><BR><BR>
<FONT class="title2">No fragmentation available for $sequence</FONT><FONT class="title3">$varMod</FONT>
$commentsString
<BR><BR>
<FONT class="title3">(<FONT color='#DD0000'>Position of at least 1 modification is unknown!</FONT>)</FONT>
</CENTER>
</BODY>
</HTML>
|;
	exit;
}

my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};
my $disableStrg=($projectAccess =~ /bioinfo|mass|manag/)? '' : '&disStore=disabled'; # not for super_bio
$dataFile=~ s/#/%23/g;
my $urlString="$promsPath{cgi}/peptide_view.cgi?file=$dataFile&query=$queryNum&hit=$rank&ID=$objectID&CALL=$call&TYPE=$fileFormat&PROJECT_ID=$projectID$disableStrg";
#if ($fileFormat=~/\.DAT/ || $fileFormat=~/\.PDM/) { # Mascot or Proteome Discover
#	my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
#	my $projectAccess=${$userInfo[2]}{$projectID};
#	$disableStrg=($projectAccess eq 'bioinfo' || $projectAccess eq 'mass')? '' : '&disStore=disabled';
#	$urlString="$promsPath{cgi}/peptide_view.cgi?file=$dataFile&query=$queryNum&hit=$rank&ID=$objectID&CALL=$call$disableStrg&TYPE=$fileFormat&PROJECT_ID=$projectID";
#}
#elsif ($fileFormat eq 'MASCOT.XML') {
#	my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
#	my $projectAccess=${$userInfo[2]}{$projectID};
#	$disableStrg=($projectAccess eq 'bioinfo' || $projectAccess eq 'mass')? '' : '&disStore=disabled';
#	$urlString="$promsPath{cgi}/peptide_view.cgi?file=$dataFile&query=$queryNum&hit=$rank&ID=$objectID&CALL=$call$disableStrg&TYPE=$fileFormat";
#}
#elsif ($fileFormat eq 'PHENYX.XML') {
#	my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
#	my $projectAccess=${$userInfo[2]}{$projectID};
#	$disableStrg=($projectAccess eq 'bioinfo' || $projectAccess eq 'mass')? '' : '&disStore=disabled';
#	$urlString="$promsPath{cgi}/peptide_view.cgi?file=$dataFile&query=$queryNum&hit=$rank&ID=$objectID&CALL=$call$disableStrg&TYPE=$fileFormat";
#}

#############################################################################
####>Checking if Reference Spectrum (with correct charge/massObs) exists<####
#############################################################################
my @refSpectra;
my ($refSpecID,$refSpecFile,$refQueryNum,$refRank,$refScore,$refUpUser);
my $sthRefSpec=$dbh->prepare("SELECT ID_SPECTRUM,MR_OBS,SPEC_FILE,QUERY_NUM,PEP_RANK,SCORE,UPDATE_USER FROM SPECTRUM WHERE PEP_SEQ='$sequence' ORDER BY SCORE DESC");
$sthRefSpec->execute();
while (my ($specID,$refMassObs,$fileName,$qNum,$rank,$score,$upUser)=$sthRefSpec->fetchrow_array) {
	if (abs($refMassObs-$massObs)<2) { # +/-2 daltons
		$upUser='' unless $upUser;
		unless ($refSpecID) {
			$refSpecID=$specID;
			$refSpecFile="$promsPath{spectrum}/$fileName";
			($refQueryNum,$refRank,$refScore)=($qNum,$rank,$score);
			$refUpUser=$upUser;
# 			last;
		}
		push @refSpectra,[$specID,$refSpecFile,$qNum,$rank,$score,$upUser];
	}
}
$sthRefSpec->finish;


if ($msType eq 'DIA') {
	####> FOR DIA $refSpecID = library ID
	my $sthRefSpec=$dbh->prepare("SELECT ASL.ID_SWATH_LIB, NAME FROM ANALYSIS_SWATH_LIB ASL, SWATH_LIB SL WHERE ID_ANALYSIS=$analysisID AND ASL.ID_SWATH_LIB=SL.ID_SWATH_LIB");
	$sthRefSpec->execute();
	$rank = 1;
	while (my ($libraryID, $libraryName, $upUser)=$sthRefSpec->fetchrow_array) {
		$upUser='' unless $upUser;
		my $specFile = (-e "$promsPath{data}/swath_lib/SwLib_$libraryID/$libraryName.sptxt") ? "$promsPath{data}/swath_lib/SwLib_$libraryID/$libraryName.sptxt" : "$promsPath{data}/swath_lib/SwLib_$libraryID/$libraryName.tsv";
		unless ($refSpecID) {
			$refSpecID=$libraryID;
			$refSpecFile= $specFile;
			($refQueryNum,$refRank,$refScore)=($objectID,$rank,'0') if(!$refQueryNum);
			$refUpUser=$upUser;
		}
		push @refSpectra,[$refSpecID, $specFile,$objectID,$rank++,'',$upUser];
	}
}
$dbh->disconnect;


#######################
####>Starting HTML<####
#######################
print header;
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Peptide Fragmentation Spectrum</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
var refTop=(top.promsFrame)? top : (top.opener.top.spectrumTolerance)? top.opener.top : top.opener.top.opener.top; // going to main myProMS window through Validation mode or protein quantification or sequence View 'ProteinWindow'!!!
|;
my $globalParamSpect='&tolerance="+refTop.spectrumTolerance+"&miniInt="+refTop.spectrumMinInt+"';
#my $globalParamSpect = ($call eq 'rank' || $call eq 'ana')? '&tolerance="+top.spectrumTolerance+"&miniInt="+top.spectrumMinInt':'"' ;
##my $globalParamRef = (($call ne "pep") && ($call ne "seq") )?'&tolerance="+top.refSpecTolerance+"&miniInt="+top.refSpecMinInt':'"' ;
if ($refSpecID) { # There are reference spectra for this peptide

	##>Storing all sprectra info in JS 2D array<##
	print "var spectra=new Array();\n";
	my $numSpec=-1;
	foreach my $refSpec (@refSpectra) {
		$numSpec++;
		print "spectra[$numSpec]=new Array();\n";
		foreach my $specInfo (@{$refSpec}) {
			print "spectra[$numSpec].push('$specInfo');\n";
		}
	}

	print qq
|var currentPos=0; // default starting reference spectrum position
function loadRefSpectrum(specPos) {
	refSpecFrame.location="$promsPath{cgi}/peptide_view.cgi?REF=1&file="+spectra[specPos][1]+"&query="+spectra[specPos][2]+"&hit="+spectra[specPos][3]+"&REF_SC="+spectra[specPos][4]+"&TYPE=REFSPEC&UP_USER="+spectra[specPos][5]+"&REF_POS="+specPos+"_$numSpec$disableStrg&width="+window.innerWidth;
	currentPos=specPos; // update reference spectrum position
}
function loadQueryFrames() {
	querySpecFrame.location="$urlString&REF=-1&REF_POS=0_$numSpec$disableStrg$globalParamSpect&width="+window.innerWidth ;
	refSpecFrame.location="$promsPath{cgi}/peptide_view.cgi?REF=1&file=$refSpecFile&query=$refQueryNum&hit=$refRank&TYPE=REFSPEC&REF_SC=$refScore&UP_USER=$refUpUser&REF_POS=0_$numSpec$disableStrg$globalParamSpect&width="+window.innerWidth;
}

var refSP,querySP;
var minRefVal,maxRefVal,minQueryVal,maxQueryVal;
function unifyRange(SP) {
	if (SP.isReference) {
		refSP=SP;
		minRefVal=SP.minValueX;
		maxRefVal=SP.maxValueX;
	}
	else {
		querySP=SP;
		minQueryVal=SP.minValueX;
		maxQueryVal=SP.maxValueX;
	}
//console.log('RANGE: ('+SP.isReference+') '+SP.unifyView);
	if (minRefVal != null && maxRefVal != null && minQueryVal != null && maxQueryVal != null) {
		refSP.minValueX=querySP.minValueX=Math.min(minRefVal,minQueryVal);
		refSP.maxValueX=querySP.maxValueX=Math.max(maxRefVal,maxQueryVal);
		return true;
	}
	else {return false;}
}
function unifyViews(SP) {
	var fromSP,toSP,toFrame;
	if (SP.isReference) {fromSP=refSP; toSP=querySP; toFrame=querySpecFrame;}
	else {fromSP=querySP; toSP=refSP; toFrame=refSpecFrame;}
//console.log('VIEW: '+fromSP.unifyView+' -> '+toSP.unifyView);
	for (var prop in fromSP.chartSettings.current.X) { // copy X axis properties
		toSP.chartSettings.current.X[prop]=fromSP.chartSettings.current.X[prop];
	}
//console.log(toFrame.SP.isReference);
	toSP.unifyView=false; // prevents infinite unify loop!!!
	toFrame.clearChart(toSP);
	toFrame.plotChart(toSP);
	toSP.unifyView=true;
}
</SCRIPT>
</HEAD>
<FRAMESET rows="50%,*" onload="loadQueryFrames()">
	<FRAME name="refSpecFrame">
	<FRAME name="querySpecFrame">
</FRAMESET>
</HTML>
|;
}
else { # no reference spectra
	print qq
|window.location="$urlString&REF=0$globalParamSpect&width="+window.innerWidth;
</SCRIPT>
</HEAD></HTML>
|;
}

####>Revision history<####
# 1.8.6 [FEATURE] Handles Spectronaut library spectrum visualization (VS 06/06/20)
# 1.8.5 [Fix]Minor bug for undefined $msType when called from validation mode (PP 30/04/18)
# 1.8.4 Minor modif to allow DIA reference spectrum drawing (MLP 19/12/17)
# 1.8.3 Change to match use of mqpar.xml file in ANALYSIS.DATA_FILE field for MaxQuant (PP 02/02/17)
# 1.8.2 Compatible with MaxQuant msms.txt file (PP 28/11/16)
# 1.8.1 Update to allow call from protein quantification raw peptide data (PP 23/09/15)
# 1.8.0 Added JS functions to manage SVG spectra (PP 26/03/14)
# 1.7.1 Remove VAR_MOD from script (GA 22/05/13)
# 1.7.0 New GPL license (PP 12/03/13)
# 1.6.9 Handles new search file path for fully validated analysis (PP 04/03/13)
# 1.6.8 Minor modification for PARAGON searches (GA 11/12/12)
# 1.6.7 Minor change to protect new PDM indexed names (the '#' is not understood in URL and launch an Internal Servor Error)<BR>See storeAnalysis modification 2.8.2 (GA 24/07/12)
# 1.6.6 Minor change to show comments to make it visible not only in validation mode but as well in protein view (CALL=pep) (GA 04/06/12)
# 1.6.5 Change printing to show comments (GA 24/05/12)
# 1.6.4 Skip java activation test (PP 13/04/12)
# 1.6.3 Manager has full access validation data (PP 01/10/11)
# 1.6.2 bug unit value on $varMod (PP 05/09/11)
# 1.6.1 no longer uses DATA_FILE from QUERY_VALIDATION table (PP 14/06/11)
# 1.6.0 SQ & full datafile conservation managements (PP 01/06/11)
