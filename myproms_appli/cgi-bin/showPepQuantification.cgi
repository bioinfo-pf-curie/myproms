#!/usr/local/bin/perl -w

#############################################################################
# showPepQuantification.cgi        1.9.10                                   #
# Authors: P. Poullet, G. Arras, M. Le Picard, V. Sabatet (Institut Curie)  #
# Contact: myproms@curie.fr                                                 #
# Displays peptide quantification data                                      #
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

#TODO: Fix problem with export: JS links in Excel!!!
$|=1;
use warnings;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use utf8; # Tells Perl that characters are UTF-8. Necessary for Excel export to work with UTF-8 characters
use Spreadsheet::WriteExcel;
use File::Path qw(rmtree); # remove_tree
use File::Copy;
#exit;

#print header; warningsToBrowser(1); #DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
# Limit ratios for +/-infinity switch in volcano plot
my $MAX_INF_RATIO=1000;
my $MIN_INF_RATIO=1/$MAX_INF_RATIO; 
my %xicSoftware=('PD'=>'Proteome Discoverer','MCQ'=>'MassChroQ','MAS'=>'Mascot','PAR'=>'Paragon','?'=>'Unknown','PKV'=>'PeakView','MQ'=>'MaxQuant','SKY'=>'Skyline','OS'=>'OpenSwath');

####################
####>Parameters<####
####################
my $call=param('CALL') || 'quanti'; # ana or quanti
my $selQuantifID=(param('id_quantif'))? param('id_quantif') : 0; # 0=> call=ana
my $action=(param('ACT'))? param('ACT') : ($selQuantifID)? 'view' : 'select';
my $analysisID=(param('id_ana'))? param('id_ana') : 0; # 0=> call=quanti

####>Connect to the database
my $dbh=&promsConfig::dbConnect;

###> Export quanti-values settings
#my ($workbook,$worksheet,$xlsRow,$rowColor1,$rowColor2,$xlsName,$xlsCol,$maxProtColLength,$maxPepColLength,%format,%itemFormat);
my ($workbook,%itemFormat);
my $doExport=param('export')? 1 : 0;

#my $act='display';# parameter pass to displayQuantification() javascript method
my $maxNumProtRT=10; # Maximum number of proteins to print RT values for XIC extractions


################
####>Main 1<####
################
if ($selQuantifID eq 'XIC') {&showLog10Plot($analysisID); exit;}

my ($quantifType,$quantifAnnot,$quantifMethDesc,$selQuantifName)=$dbh->selectrow_array("SELECT M.CODE,QUANTIF_ANNOT,M.DES,Q.NAME FROM QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$selQuantifID");
my ($numAnaUsed)=$dbh->selectrow_array("SELECT COUNT(*) FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
my ($projectID)=&promsMod::getProjectID($dbh,$selQuantifID,'quantification');

my (%quantifParamInfo, %anaQuantifList, %quantificationTypes, %methodList, %condToState, %sourceFiles);

my (@anaOrder, %peptideInAna, %listAna); # Stores analysis data
my (%selMGTopProt, %trueMGTopProt, %maxProtMatch, %idtoidentifier, %protInfo, %posBeg, %pepProtPos); # Stores proteins data
my (%pepAll, %pepInfo, %peptideQuant, %quantifData, %peptideFragments, %peptideSets, %labeledPeptideSets, %nbQuantifPeptides, %peptideMrobs, %peptideRT, %peptideScore, %pepDataSource); # Stores peptides data

my $title=($numAnaUsed==1)? 'Single-Analysis Peptide Quantification' : 'Multi-Analysis Quantification';
my $isTrace=0; # ($quantifAnnot =~ /TRACES=1/)? 1 : 0; # To know if some XIC_TRACES were computed ### DISABLED until call to PEPTIDE_QUANTIFICATION table is replaced (PP 27/06/18) ###
my $xicSoftCode;
my $xicSoftVersionStrg=' (Unknown version)';
#$quantifAnnot=~s/::SOFTWARE=([^:|;]+);?\d*\.*\d*//; # remove software info for back compatibility
$quantifAnnot=~s/::SOFTWARE=([^:|;]+);?([\d*\.?]*)//; # change on 2017/11/03 to deal with software versions with multiple dots like 2.1.0.81

if ($1) {
	$xicSoftCode=$1;
	$xicSoftVersionStrg=" ($2)" if $2;
} elsif ($quantifAnnot=~/EXTRACTION_ALGO=/) {
	$xicSoftCode='MCQ';
} else {
	my ($fileFormat)=$dbh->selectrow_array("SELECT FILE_FORMAT FROM ANALYSIS A,ANA_QUANTIFICATION Q WHERE A.ID_ANALYSIS=Q.ID_ANALYSIS AND ID_QUANTIFICATION=$selQuantifID");
	$xicSoftCode=($fileFormat=~/\.PDM\Z/)? 'PD' : ($fileFormat=~/^MASCOT/)? 'MAS' : ($fileFormat=~/PARAGON/)? 'PAR' : '?';
}
my %dataCorrection;
if ($quantifAnnot=~/::CORRECTION=([^:]+)/) {
	my ($isLocal,$coefStrg)=split(';',$1);
	$isLocal=~s/#//;
	$dataCorrection{COEF}=$coefStrg;
	if ($isLocal) {
		$dataCorrection{SRC}=1;
		($dataCorrection{NAME},$dataCorrection{PRODUCT},$dataCorrection{LOT})=$dbh->selectrow_array("SELECT IC.NAME,PRODUCT_NUMBER,LOT_NUMBER FROM ISOTOPIC_CORRECTION IC,QUANTIFICATION Q WHERE Q.ID_PRODUCT=IC.ID_PRODUCT AND ID_QUANTIFICATION=$selQuantifID");
		#$dataCorrection{COMPANY}='ThermoFisher';
	}
	else {$dataCorrection{SRC}=-1;}
}
else {$dataCorrection{SRC}=0;}

###>Fetching quantification parameters
my $sthQP2=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,P.CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD M,QUANTIFICATION_PARAMETER P WHERE Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND M.ID_QUANTIFICATION_METHOD=P.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$selQuantifID");
$sthQP2->execute;
while (my ($paramID,$paramCode)=$sthQP2->fetchrow_array) {
	$quantifParamInfo{"$paramCode"}=$paramID;
}
$sthQP2->finish;

if ($action eq 'delete') {&deleteQuantification; exit;} # only for label-based 2ndary quantifications
elsif ($action eq 'ajaxXicTrace') {&ajaxXicTrace($selQuantifID); exit;}
elsif ($action eq 'ajaxShowPepTDAGraph'){&ajaxShowPepTDAGraph; exit;}
elsif ($action eq 'ajaxShowPepTDAList'){&ajaxShowPepTDAList; exit;}
elsif ($action eq 'edit') { ###>Edition submission<###
	my $qName=$dbh->quote(param('name'));
	$dbh->do("UPDATE QUANTIFICATION SET NAME=$qName WHERE ID_QUANTIFICATION=$selQuantifID");
	$dbh->commit;
	$action='display';
}

if ($call eq 'ana') {
	my $sthAQ=$dbh->prepare("SELECT Q.ID_QUANTIFICATION,Q.NAME,Q.FOCUS,STATUS,M.ID_QUANTIFICATION_METHOD,M.NAME,M.CODE,Q.QUANTIF_ANNOT FROM ANA_QUANTIFICATION A,QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE ID_ANALYSIS=$analysisID AND A.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_DESIGN IS NULL AND Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD");
	$sthAQ->execute;
	while (my ($quantifID,$name,$focus,$status,$methodID,$methodName,$methodCode,$qAnnot)=$sthAQ->fetchrow_array) {
		@{$anaQuantifList{$quantifID}}=($name,$focus,$status,$methodID,$qAnnot);
		@{$methodList{$methodID}}=($methodName,$methodCode);
		my $typeName=($focus eq 'peptide' && $methodCode=~/XIC|SIN|SILAC|ITRAQ|TMT|TDA/)? 'Peptide quantification'
			: ($focus eq 'peptide' && $methodCode eq 'DIA')? 'SWATH-MS/DIA'
			: ($focus eq 'peptide')? 'Peptide ratio'
			: ($focus eq 'protein' && $methodCode=~/PROT_RATIO/)? 'Protein ratio' # should never be TNPQ (call=ana)
			: ($focus eq 'protein')? 'Protein quantification' # EMPAI
			: 'Unknown';
		push @{$quantificationTypes{$typeName}},$quantifID;
		#push @{$quantificationTypes{'Peptide quantification'}},$quantifID if ($designID && $methodCode=~/PROT_RATIO|TNPQ/);
	}
	$sthAQ->finish;
}
else { #quanti
	my ($methodID,$methodName,$methodCode)=$dbh->selectrow_array("SELECT M.ID_QUANTIFICATION_METHOD,M.NAME,M.CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=$selQuantifID");
	@{$methodList{$methodID}}=($methodName,$methodCode);
	($analysisID,@{$anaQuantifList{$selQuantifID}})=$dbh->selectrow_array("SELECT A.ID_ANALYSIS,Q.NAME,Q.FOCUS,STATUS,M.ID_QUANTIFICATION_METHOD FROM ANA_QUANTIFICATION A,QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE A.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=$selQuantifID LIMIT 0,1");
	# $analysisID = first matching analysis
}

&getAnalysisDataSources($analysisID, $projectID, \%sourceFiles);

$dbh->disconnect if $action eq 'select'; #'display'; # $selQuantifID;

#######################
####>Starting HTML<####
#######################

if ($doExport) {
	# Information for the web-browser otherwise, an error will occur
	my $timeStamp1=strftime("%Y%m%d %H-%M",localtime);

 	$workbook=Spreadsheet::WriteExcel->new("-");
	$workbook->set_properties(title=>"Peptide quantification data",
							author=>'myProMS server',
							comments=>'Automatically generated with Perl and Spreadsheet::WriteExcel'
							);
	$workbook->set_custom_color(40,189,215,255); # light light color  #BDD7FF (new myProMS colors V3.5+)
	#$workbook->set_custom_color(40,224,224,255); # light light color #E0E0FF
	%itemFormat=(
			title =>			$workbook->add_format(align=>'center',size=>18,bold=>1,border=>1),
			header =>			$workbook->add_format(align=>'center',size=>10,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			headerR =>			$workbook->add_format(align=>'right',valign=>'top',size=>10,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			mergeRowHeader =>	$workbook->add_format(align=>'center',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			mergeColHeader =>	$workbook->add_format(align=>'center',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			mergeColHeaderL =>	$workbook->add_format(align=>'left',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			text =>				$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>0,border=>1),
			textC =>			$workbook->add_format(align=>'center',size=>10,valign=>'top',text_wrap=>0,border=>1),
			textWrap =>			$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>1,border=>1),
			#mergeText =>		$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>0,border=>1),
			mergeColText =>		$workbook->add_format(align=>'left',size=>10,text_wrap=>1,border=>1),
			#numberVis =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,bold=>1),
			#numberHid =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1),
			#numberVisBC =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,bold=>1,color=>'grey'),
			#numberHidBC =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,color=>'grey'),
			#numberVisVP =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,bold=>1,italic=>1),
			#numberHidVP =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1,italic=>1),
			number =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1),
			#mergeNumber =>		$workbook->add_format(size=>10,align=>'center',valign=>'top',border=>1),
			number1d =>			$workbook->add_format(size=>10,align=>'center',valign=>'top',num_format=>'0.0',border=>1),
#mergeNumber1d =>	$workbook->add_format(size=>10,align=>'center',valign=>'top',num_format=>'0.0',border=>1)
			);

	print header(-type=>"application/vnd.ms-excel",-attachment=>"Peptide_quantification_$timeStamp1.xls");

}
else {
	print header(-'content-encoding'=>'no',-charset=>'utf-8',-cache_control=>"no-cache") ;
	warningsToBrowser(1); # writes PERL warnings as HTML comments <!-- --> must be activated after header
	print qq |
	<HTML>
		<HEAD>
			<TITLE>Peptide Quantification Results</TITLE>
			<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
			<STYLE type="text/css">
				.TD {font-weight:normal;}
				.TH {font-weight:bold;}
				.LINK {cursor:pointer;}
				.popup {z-index:999;background-color:#FFFFFF;border:solid 3px #999999;padding:5px;position:absolute;display:none;box-shadow:10px 10px 20px #808080;}
			</STYLE>
			<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
			<SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
			<SCRIPT src="$promsPath{html}/js/local/heatMap.js"></SCRIPT>
			<SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
			<SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
			<SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
			
			<SCRIPT LANGUAGE="JavaScript">
|;
	&promsMod::popupInfo();
	if ($call eq 'ana') {
		print qq |
			function displayQuantification(quantifInfoStrg,action) {
				var quantifInfo=quantifInfoStrg.split(':'); // focus:quantifID
				var focus=quantifInfo[0], quantifID=quantifInfo[1];
				if (quantifID==0) {action='select';}
				if (action=='delete' && !confirm('Delete selected quantification?')) {return;}
				if (focus=='peptide') {
					window.location="$promsPath{cgi}/showPepQuantification.cgi?CALL=$call&id_ana=$analysisID&id_quantif="+quantifID+"&ACT="+action;
				}
				else {
					window.location="$promsPath{cgi}/showProtQuantification.cgi?CALL=$call&id_ana=$analysisID&id_quantif="+quantifID+"&ACT="+action+"&view=graph";
				}
			}
		|;
			}
		print qq|
			function sequenceView(id_protein,anaIdStrg){
				var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+anaIdStrg+"&id_prot="+id_protein+"&msdata="+top.showMSdata;
				top.openProtWindow(winLocation);
			}
			
			
			function updateSettings(act) {
				if (act=='more') {
					document.getElementById('moreSettings').style.display='none';
					document.getElementById('lessSettings').style.display='block';
					document.getElementById('advancedSetDIV').style.display='block';
				}
				else {
					document.getElementById('moreSettings').style.display='block';
					document.getElementById('lessSettings').style.display='none';
					document.getElementById('advancedSetDIV').style.display='none';
				}
			}
			
			function exportQuanti() {
				window.location="$promsPath{cgi}/showPepQuantification.cgi?CALL=$call&id_ana=$analysisID&id_quantif=$selQuantifID&ACT=$action&export=1";
			}
			
			var XHR=null;
			var isNav = (navigator.appName.indexOf("Netscape") !=-1);
			function ajaxDrawXICData(e,pepID) {
				var displayDiv=document.getElementById('displayDIV');
			
				var divX,divY;
				divX = (isNav)? e.pageX : event.clientX + document.body.scrollLeft; //divX-=5;
				divY = (isNav)? e.pageY : event.clientY + document.body.scrollTop; //divY+=10;
				displayDiv.style.left = divX + 'px';
				displayDiv.style.top = divY + 'px';
				displayDiv.style.display='block';
			
				//wait image
				var infoDIV=document.getElementById('infoDIV');
				infoDIV.innerHTML="<BR><BR><IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR>";
			
				//If XHR object already exists, the request is canceled & the object is deleted
				if(XHR && XHR.readyState != 0){
					XHR.abort();
					delete XHR;
				}
			
				//Creation of the XMLHTTPRequest object
				XHR = getXMLHTTP();
				if (!XHR) {
					return false;
				}
				var URL = "$promsPath{cgi}/showPepQuantification.cgi?CALL=$call&id_quantif=$selQuantifID&ACT=ajaxXicTrace&id_peptide="+pepID;
			
				XHR.open("GET",URL,true);
				XHR.onreadystatechange=function() {
					if (XHR.readyState==4 && XHR.responseText){
						infoDIV.innerHTML=XHR.responseText;
						reloadImage('localPicture'); // Need to reload so as to avoid to print cache images
						reloadImage('globalPicture');
					}
				}
				XHR.send(null);
			}
			
			function ajaxShowGraphicPepTDA (protID) {
				if(document.getElementById('TDAButtonGraph'+protID).value=='Show graphic') {
					document.getElementById('listTDA'+protID).style.display='none';
					document.getElementById('graphPepTDA'+protID).style.display='';
					document.getElementById('TDAButtonGraph'+protID).value='Hide graphic';
					document.getElementById('TDAButtonList'+protID).value='Show list';
					
					if (protID && document.getElementById('graphPepTDA' + protID).innerHTML == ""){
						document.getElementById('graphPepTDA'+protID).innerHTML="<BR><BR><IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR>";
						var XHR = getXMLHTTP();
						XHR.open("GET", "$promsPath{cgi}/showPepQuantification.cgi?id_quantif=$selQuantifID&ACT=ajaxShowPepTDAGraph&id_protein="+protID , true);
						XHR.onreadystatechange=function() {
							if (XHR.readyState==4 && XHR.responseText) {
								eval(XHR.responseText);
							}
						}
						XHR.send(null);
					}
				}else if(document.getElementById('TDAButtonGraph'+protID).value=='Hide graphic'){
					document.getElementById('graphPepTDA'+protID).style.display='none';
					document.getElementById('TDAButtonGraph'+protID).value='Show graphic';
				}
			}
			
			function ajaxShowListPepTDA(protID, anaID = ''){
				if(document.getElementById('TDAButtonList'+protID).value=='Show list'){
					document.getElementById('TDAButtonGraph'+protID).value='Show graphic';
					document.getElementById('graphPepTDA'+protID).style.display='none';
					document.getElementById('listTDA'+protID).style.display='';
					document.getElementById('TDAButtonList'+protID).value='Hide list';
					
					if (protID && document.getElementById('listTDA'+protID).innerHTML == ""){
						document.getElementById('listTDA'+protID).innerHTML="<BR><BR><IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR>";
						var XHR = getXMLHTTP();
						XHR.open("GET", "$promsPath{cgi}/showPepQuantification.cgi?id_quantif=$selQuantifID&ACT=ajaxShowPepTDAList&id_protein="+protID+"&id_analysis="+anaID, true);
						XHR.onreadystatechange=function() {
							if (XHR.readyState==4 && XHR.responseText) {
								document.getElementById('listTDA'+protID).innerHTML=XHR.responseText;
							}
						}
						XHR.send(null);
					}
				}else if(document.getElementById('TDAButtonList'+protID).value=='Hide list'){
					document.getElementById('listTDA'+protID).style.display='none';
					document.getElementById('TDAButtonList'+protID).value='Show list';
				}
			}
			
			function ajaxShowFragTDA(e,anaPepId){
				var displayDiv=document.getElementById('displayDIV');
				displayDiv.style.display='none';
				displayDiv.style.width='auto';
				var infoDiv=document.getElementById('infoDIV');
				infoDiv.innerHTML="<BR><BR><IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR>";
			
				var left=e.clientX;
				var top=e.clientY;
				if(document.body.scrollTop){
					top+=document.body.scrollTop
				}
				displayDiv.style.left=left+'px';
				displayDiv.style.top=top+'px';
				displayDiv.style.display='block';
			
			
				if (anaPepId){
					var XHR = getXMLHTTP();
					XHR.open("GET", "$promsPath{cgi}/showPepQuantification.cgi?id_quantif=$selQuantifID&ACT=ajaxShowFragTDA&anaPepId="+anaPepId, true);
					XHR.onreadystatechange=function() {
						if (XHR.readyState==4 && XHR.responseText) {
							var codeParts=XHR.responseText.split('#==========#');
							infoDiv.innerHTML=codeParts[0];
							eval(codeParts[1]);
							var newWidth=Math.round(displayDiv.offsetWidth/10);
							displayDiv.style.width=displayDiv.offsetWidth + newWidth + 'px';
						}
					}
					XHR.send(null);
				}
			}
			
			function getXMLHTTP(){
				var xhr=null;
				if(window.XMLHttpRequest) {// Firefox & others
					xhr = new XMLHttpRequest();
				}
				else if(window.ActiveXObject){ // Internet Explorer
					try {
					  xhr = new ActiveXObject("Msxml2.XMLHTTP");
					} catch (e) {
						try {
							xhr = new ActiveXObject("Microsoft.XMLHTTP");
						} catch (e1) {
							xhr = null;
						}
					}
				}
				else { // XMLHttpRequest not supported by browser
					alert("Your browser does not support XMLHTTPRequest objects...");
				}
				return xhr;
			}
			
			var spectWin;
			function drawSpectrum(pepId,pepInfo) {
				var paramString="RID="+pepInfo+"&CALL=pep";
				spectWin=window.open("$promsPath{cgi}/drawSpectrum.cgi?"+paramString,'SpectrumWindow','width=950,height=950,location=no,resizable=yes,scrollbars=yes');
				spectWin.focus();
			}
			
			function reloadImage(imgID){
				var reloadImage=document.getElementById(imgID);
				reloadImage.src = reloadImage.src + "?" + new Date().getTime();
			}
		</SCRIPT>	
	</HEAD>
	<BODY background='$promsPath{images}/bgProMS.gif'> <!-- Do not add onload code here! Use window.onload instead. -->
		<CENTER>
			<FONT class="title">$title</FONT>
			<DIV id="displayDIV" class="popup"> <!--filter:alpha(opacity=80);opacity:0.8;-->
				<DIV id="infoDIV"></DIV>
			</DIV><BR/>
|;
}


#print "*** ACTION=$action ***<BR>\n";

################
####>Main 1<####
################
#if (!$designID) {
#	print qq
#|<TABLE bgcolor="$darkColor">
#<TR><TD class="title2" align=right>Select:</TD><TD><SELECT class="title2" onchange="displayQuantification(this.value,'display','$view')"><OPTION value="">-= Select =-</OPTION>
#|;
#}
#else {
#	print qq
#|<TABLE bgcolor="$darkColor">
#<TR><TD class="title2" align=right>Select:</TD><TD><SELECT class="title2" onchange="displayQuantification($selQuantifID,this.value,'graph')"><OPTION value="">-= Select =-</OPTION>
#|;
#}

if ($call eq 'ana' && !$doExport) {
	print qq
|<TABLE bgcolor="$darkColor">
<TR><TD class="title2" align=right>Select:</TD><TD><SELECT class="title2" onchange="displayQuantification(this.value,'display')"><OPTION value="peptide:0">-= Select =-</OPTION>
| ;
	foreach my $typeName (sort keys %quantificationTypes) {
		print "<OPTGROUP label=\"$typeName:\">\n";
		foreach my $quantifID (sort{lc($anaQuantifList{$a}[0]) cmp lc($anaQuantifList{$b}[0])} @{$quantificationTypes{$typeName}}) {
			print "<OPTION value=\"$anaQuantifList{$quantifID}[1]:$quantifID\"";
			my $errorFlag='';
			if ($anaQuantifList{$quantifID}[2] == -1 || $anaQuantifList{$quantifID}[2] == 0) {
				print ' disabled'; # just created or running
			}
			elsif ($anaQuantifList{$quantifID}[2] == -2) { # Error
				$errorFlag='***';
			}
			print ' selected' if $quantifID==$selQuantifID;
			my $methodName=$methodList{$anaQuantifList{$quantifID}[3]}[0];
			my $qAnnot=$anaQuantifList{$quantifID}[4];
			my $softwareExt='';
			if ($methodList{$anaQuantifList{$quantifID}[3]}[1] eq 'XIC'){ # Need to distinguish between MassChroQ, MaxQuant and Proteome-Discoverer
				$softwareExt=($qAnnot =~ /SOFTWARE=MCQ/) ? '-MCQ' : ($qAnnot =~ /SOFTWARE=MQ/) ? '-MQ' : ($qAnnot =~ /SOFTWARE=PD/) ? '-PD' : '';
			}
			print ">$errorFlag$anaQuantifList{$quantifID}[0] [$methodName]$softwareExt$errorFlag</OPTION>\n";
		}
		print "</OPTGROUP>\n";
	}
	my $notDeletable=1;
	if ($selQuantifID) {
		my ($methodCode)=$dbh->selectrow_array("SELECT CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=$selQuantifID");
		#if (($methodCode =~ /SWATH/ && $xicSoftCode eq 'PKV') || $xicSoftCode eq 'SKY' || $xicSoftCode eq 'OS') {#}	# SWATH PKV is never deletable
		if (($methodCode=~/SWATH/ && $xicSoftCode eq 'PKV') || $xicSoftCode=~/^(SKY|OS|MQ)$/) {
			$notDeletable=1;
		}
		elsif ($methodCode !~ /SILAC|ITRAQ|TMT/ && $xicSoftCode ne 'PD') { # SILAC & ITRAQ & label-free PD are never deletable
			($notDeletable)=$dbh->selectrow_array("SELECT 1 FROM PARENT_QUANTIFICATION WHERE ID_PARENT_QUANTIFICATION=$selQuantifID LIMIT 0,1");
		}
	}
	my $disabDeleteStrg=($notDeletable)? ' disabled' : '';

	print qq
|</SELECT></TD>
<TD><INPUT type="button" value="Edit" onclick="document.getElementById('editDIV').style.display='block'"/></TD>
<TD bgcolor="#DD0000"><INPUT type="button" value="Delete" style="color:#DD0000" onclick="displayQuantification('peptide:$selQuantifID','delete')"$disabDeleteStrg/></TD>
</TR></TABLE><BR>
|;


	if ($action eq 'select') {
		print "</CENTER>\n</BODY>\n</HTML>\n";
		exit;
	}

	####>Hidden edition form<####
	print qq
|<FORM name="editForm" method="POST">
<INPUT type="hidden" name="ACT" value="edit">
<INPUT type="hidden" name="CALL" value="ana">
<INPUT type="hidden" name="id_ana" value="$analysisID">
<INPUT type="hidden" name="id_quantif" value="$selQuantifID">
<DIV id="editDIV" style="display:none">
<TABLE bgcolor=$darkColor border=0>
<TR><TH align=right>&nbsp;Name :</TH><TD bgcolor=$lightColor><INPUT type="text" name="name" value="$anaQuantifList{$selQuantifID}[0]" size="50"></TD></TR>
<TR><TH colspan=2><INPUT type="submit" value=" Save ">&nbsp;&nbsp;&nbsp;<INPUT type="button" value=" Cancel " onclick="document.getElementById('editDIV').style.display='none'"/></TR>
</TABLE>
</DIV>
</FORM>
|;
}

###################################################
####>Fetching data for selected Quantification<####
###################################################
print qq
|<DIV id="waitDIV">
<BR><BR><BR><BR><BR><FONT class="title3">Fetching data...</FONT><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">
<BR><BR><FONT class="title3">Status:&nbsp;<SPAN id="waitSPAN"><SPAN>...</FONT>
</DIV>
| unless $doExport;

my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
my ($labelType)=($labelStrg)? ($labelStrg=~/LABEL=(.+)/) : ('FREE');
$labelType=uc($labelType); # iTRAQ -> ITRAQ


#################################
####>LABEL FREE/DESIGN BLOCK<####
#################################
#if ($quantifType =~ /(PROT_RATIO|TNPQ|XIC|SIN|EMPAI)/){ #}
if ($labelType eq 'FREE') {
	#####
	# In label-free experiments, several analysis are used, therefore, a random ID_ANALYSIS
	# is set to let the javascript function works
	#####
	#if (!$analysisID) {
	#	($analysisID)=$dbh->selectrow_array("SELECT ANALYSIS.ID_ANALYSIS FROM ANALYSIS,ANA_QUANTIFICATION WHERE ANA_QUANTIFICATION.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ID_QUANTIFICATION IN (SELECT ID_PARENT_QUANTIFICATION FROM PARENT_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID) ORDER BY NAME ASC LIMIT 1");
	#}

	my ($viewProtanaID, $count, $totAnalyses, $numAna);
	my $multiAnaQuanti = ($call eq 'quanti') ? 1 : 0; # 'quanti' => Multi Analysis Quantification || 'ana' => Analysis level quantification
	my $shouldLoadData = ($quantifType =~ /^(SWATH|TDA|DIA)$/ || $xicSoftCode eq 'SKY') ? 0 : 1;
	my $nbTransition   = 1;
	
	printWaitingMsg("Fetching data...", 0);
	
	if ($multiAnaQuanti) { # Fetch all analysis and corresponding peptides from the quantification
		&fetchQuantiAnalysis($selQuantifID);
	} else { # Internal Quantifications
		@anaOrder = ($analysisID);
	}
	
	## Fetching Analysis
	$totAnalyses = scalar @anaOrder;
	$numAna = 0;
	foreach my $anaID (@anaOrder) {
		my $anaName = $listAna{$anaID};
		$viewProtanaID=$anaID unless $viewProtanaID;# Arbitrary chose an ID_ANALYSIS to print the protein sequence!
		$numAna++;

		printWaitingMsg("Fetching data: Analysis $numAna/$totAnalyses...", 0);
		
		### Fetch all proteins for current analysis
		fetchAnaProteins($anaID, $shouldLoadData); 
		
		### Fetch peptides from all proteins
		if($shouldLoadData) {
			printWaitingMsg('/', 1);
			fetchPeptidesFromFile($anaID);
		}
	}
	
	if(!$multiAnaQuanti) {
		## Fetching peptide info & grouping isoforms (+/- labeled)
		my $sthPI=$dbh->prepare("SELECT P.ID_PEPTIDE,PEP_SEQ,GROUP_CONCAT(PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),MR_OBS,ELUTION_TIME,CHARGE,GROUP_CONCAT(DISTINCT(ABS(PEP_BEG)) ORDER BY ABS(PEP_BEG) SEPARATOR ','),SCORE,DATA
									FROM PEPTIDE P
									LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
									INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
									WHERE P.ID_ANALYSIS=$analysisID AND ID_PROTEIN=? GROUP BY P.ID_PEPTIDE ORDER BY PEP_SEQ,CHARGE");

		foreach my $matchGroup (keys %trueMGTopProt) {
			$nbQuantifPeptides{$trueMGTopProt{$matchGroup}} = 0;
			my $identifier = $idtoidentifier{$trueMGTopProt{$matchGroup}};
			
			$sthPI->execute($trueMGTopProt{$matchGroup}); #$protID
			while (my ($pepID,$pepSeq,$modCode,$mrObs,$rtSc,$charge,$beg,$score,$pepData)=$sthPI->fetchrow_array) {
				$pepData = '' unless $pepData;
				my ($srcRk) = ($pepData =~ /SOURCE_RANK=(\d+)/); $srcRk = 0 unless $srcRk; # NOT ACTUALLY USED!!!
				$pepDataSource{$pepID} = $sourceFiles{$analysisID}{$srcRk}; # use for export only
				
				#my $varModStrg=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$analysisID,$pepSeq);
				my $varModStrg = ($modCode) ? ' + '.&promsMod::decodeVarMod($dbh, $pepSeq, $modCode) : '';
				
				my @beg = split(/,/,$beg);
				@{$posBeg{"$pepSeq$varModStrg"}} = (\@beg, $mrObs);
				$pepProtPos{$trueMGTopProt{$matchGroup}}{$pepID} = $beg; # prot ID ; pepID ; beg

				$peptideMrobs{$pepID} = ($mrObs) ? $mrObs : '-';
				$peptideScore{$pepID} = ($score) ? $score : 0; # ghost peptides don't have scores

				#(my $rt)=($rtSc=~m/sc/) ? ($rtSc=~/sc.+;et(.+)/) : $rtSc;
				my $rt=(!$rtSc) ? 0 : ($rtSc=~/et([\d\.-]+)/) ? $1 : ($rtSc=~/^([\d\.-]+)/)? $1 : 0; #$rtSc; # force to 0 because unrecognized (PP 18/04/18)
				$peptideRT{$pepID} = $rt; # can be empty
				
				push @{$peptideSets{$pepID}}, [$pepID, $charge, $beg[0], $pepSeq, $varModStrg];
				push @{$labeledPeptideSets{$matchGroup}{"$pepSeq$varModStrg"}{$charge}}, $pepID;
			
				if ($maxProtMatch{$identifier}) {
					$maxProtMatch{$identifier}{"$pepSeq$varModStrg"} += 1;
				}
				else {
					$maxProtMatch{$identifier}{"$pepSeq$varModStrg"} = 1;
				}
				
				$nbQuantifPeptides{$trueMGTopProt{$matchGroup}}++;
			}
		}
		$sthPI->finish;
	}
	
	###>Displaying summary<###
	print qq |
		<SCRIPT LANGUAGE="Javascript">document.getElementById('waitDIV').style.display='none';</SCRIPT>
	| unless $doExport;
	&displayPeptideSummary($labelType,$xicSoftCode);
	
	if ($doExport) {
		if (!$shouldLoadData) {
			fetchAnaProteins($analysisID, 1) if ($multiAnaQuanti);
			$nbTransition = &fetchTDAPeptides($analysisID);
		}
		
		my %exportParameters=('refQuantifData'=>\%quantifData,
							'refPepDataSource'=>\%pepDataSource,
							'refProtInfo'=>\%protInfo,
							'refPosBeg'=>\%posBeg,
							'reflabeledPeptideSets'=>\%labeledPeptideSets,
							'refPepMrObs'=>\%peptideMrobs,
							'refPeptideScore'=>\%peptideScore,
							'refPeptideSets'=>\%peptideSets,
							#'refPeptideRT'=>\%peptideRT,
							'refPeptideFragments'=>\%peptideFragments,
							'nbTransition'=>$nbTransition);
		
		if($multiAnaQuanti) {
			$exportParameters{'refProtMatch'}  = \%maxProtMatch;
			$exportParameters{'refPepAll'}     = \%pepAll;
			$exportParameters{'refListAna'}    = \%listAna;
			$exportParameters{'refAnaOrder'}   = \@anaOrder;
		} else {
			$exportParameters{'refProtMatch'}  = \%selMGTopProt;
			$exportParameters{'refPepProtPos'} = \%pepProtPos,
		}
		
		&exportProteinList(\%exportParameters);
		$workbook->close();
		exit;
	}

	###>Displaying data<###
	my $colName = ($quantifType eq 'XIC') ? "XIC" : "SI";
	my $paramCode = ($quantifType eq 'XIC') ? "XIC_AREA" : "SIN_SI";
	my $printMCQ_RT = ($xicSoftCode eq 'MCQ' && scalar keys %selMGTopProt <= $maxNumProtRT) ? 1 : 0;

	foreach my $matchGroup (sort{$a<=>$b} keys %selMGTopProt) { # ALIAS,PROT_DES,MW,PROT_LENGTH,ORGANISM //,NUM_PEP,NUM_MATCH,SCORE,CONF_LEVEL,PEP_COVERAGE,PEP_SPECIFICITY
		my ($protID, $alias, $des, $mw, $length, $org) = @{$selMGTopProt{$matchGroup}};
		my $colSpan = 10; # last col is empty for longer prot title
		
		## Displaying protein header
		if (!$shouldLoadData) { # TDA/DIA/Skyline XIC
			my $graphButtonDisplay = ($totAnalyses > 1) ? "" : "none";
			my $showListArgs = "$protID";
			$showListArgs .= ", ".$anaOrder[0] if (!$multiAnaQuanti || scalar @anaOrder == 1);
			
			print qq |
				<TABLE border=0 cellspacing=0 cellpadding=2 width="100%">
					<TR bgcolor="$darkColor">
						<TH valign=top><A href="javascript:sequenceView($protID,$analysisID)">$alias</A>:</TH><TD bgcolor="$lightColor" width=100%>$des <FONT class="org">$org</FONT></TD>
						<TD><INPUT type=\"button\" id=\"TDAButtonGraph$protID\" value=\"Show graphic\" onclick=\"ajaxShowGraphicPepTDA($protID)\" style='display: $graphButtonDisplay'></TD>
						<TD><INPUT type=\"button\" id=\"TDAButtonList$protID\" value=\"Show list\" onclick=\"ajaxShowListPepTDA($showListArgs)\"></TD>
					</TR>
				</TABLE>
				
				<DIV id=\"listTDA$protID\" style=\"display:none\"></DIV>
				<DIV id=\"graphPepTDA$protID\" style=\"display:none\"></DIV>
			|;

		} else {
			print qq |
				<TABLE border=0 cellspacing=0 cellpadding=2 width="100%">
					<TR bgcolor="$darkColor">
						<TD class="bBorder" colspan=$colSpan>
							<TABLE>
								<TR>
									<TH valign=top><A href="javascript:sequenceView($protID,$analysisID)">$alias</A>:</TH>
									<TD bgcolor="$lightColor" width=100%>$des <FONT class="org">$org</FONT> ($length aa)</TD>
								</TR>
							</TABLE>
						</TD>
					</TR>
			|;
			
			## Displaying fragments name as header
			if ($multiAnaQuanti && ($printMCQ_RT || $xicSoftCode eq 'MQ')) {
				print "<TR><TH colspan=6>&nbsp;&nbsp;</TH>\n";
				
				foreach my $anaID (@anaOrder) {
					### Compute fragment columns length based on their amount
					$colSpan = ($xicSoftCode eq 'MQ') ? 2 : 3; # 3: MCQ, PD; 2: MQ
					print "<TH bgcolor=\"$darkColor\" class=\"rbBorder\" colspan=$colSpan>&nbsp;$listAna{$anaID}&nbsp;</TH>\n";
				}
				
				print "</TR>\n";
			}
			
			print qq |
					<TR><TH>&nbsp;&nbsp;</TH><TH bgcolor="$darkColor" class="rbBorder">#</TH>
					<TH bgcolor="$darkColor" class="rbBorder" >&nbsp;Peptides&nbsp;</TH>
					<TH bgcolor="$darkColor" class="rbBorder">&nbsp;Mr(Obs)&nbsp;</TH>
					<TH bgcolor="$darkColor" class="rbBorder">&nbsp;Start&nbsp;</TH>
					<TH bgcolor="$darkColor" class="rbBorder">&nbsp;Charge&nbsp;</TH>
			|;
			print "<TH bgcolor=\"$darkColor\" class=\"rbBorder\">&nbsp;Scores&nbsp;</TH>\n" if (!$multiAnaQuanti);
			
			## Displaying quantification type-specific header(s)
			foreach my $anaID (@anaOrder) {
				if ($printMCQ_RT) {
					print "<TH bgcolor=\"$darkColor\" class=\"rbBorder\">&nbsp;RT<sub>beg</sub>&nbsp;</TH>\n";
					print "<TH bgcolor=\"$darkColor\" class=\"rbBorder\">&nbsp;RT<sub>end</sub>&nbsp;</TH>\n";
				} elsif($xicSoftCode eq 'MQ') {
					print "<TH bgcolor=\"$darkColor\" class=\"rbBorder\">&nbsp;RT<SUB>min</SUB></TH>";
					print "<TH bgcolor=\"$darkColor\" class=\"rbBorder\">&nbsp;$colName&nbsp;</TH>";
				} elsif ($multiAnaQuanti) { # Print Analysis Name as columns
					print "<TH bgcolor=\"$darkColor\" class=\"rbBorder\">&nbsp;$listAna{$anaID}&nbsp;</TH>\n";
				} else {
					print "<TH bgcolor=\"$darkColor\" class=\"rbBorder\">&nbsp;$colName&nbsp;</TH>\n";
				}
			}
			print qq |
						<TD width=50%></TD>
					</TR>
			|;
		
			if($multiAnaQuanti) {
				&displayQuantificationValues($protID, $colSpan);
			
				print qq |
						<TR><TD colspan=$colSpan>&nbsp;</TD></TR>
					</TABLE>
				|;
			} else {
				my $bgColor=$lightColor;
				my $pepCount=0;
				
				foreach my $seqVarMod (sort{$posBeg{$a}[0][0]<=>$posBeg{$b}[0][0] || lc($a) cmp lc($b) || $a cmp $b} keys %{$maxProtMatch{$idtoidentifier{$protID}}} ){
					foreach my $charge (sort keys %{$labeledPeptideSets{$matchGroup}{$seqVarMod}}) {
						foreach my $pepID (@{$labeledPeptideSets{$matchGroup}{$seqVarMod}{$charge}}) {
							$pepCount++;
							my @pepInfo=@{$peptideSets{$pepID}[0]};
							my $score=($peptideScore{$pepID})? (sprintf "%.2f",$peptideScore{$pepID})*1 : '-';
							my $mrObs=($peptideMrobs{$pepID})? (sprintf "%.2f",$peptideMrobs{$pepID})*1 : '-';
							print qq
	|<TR bgcolor="$bgColor" class="list"><TD bgcolor="#FFFFFF"></TD><TD class="rBorder" align=right>&nbsp;$pepCount&nbsp;</TD><TH class="font11" align=left nowrap>$seqVarMod&nbsp;</TH><TD align=center>$mrObs</TD><TD align=center >$pepProtPos{$protID}{$pepID}</TD><TD align=center>$charge<SUP>+</SUP></TD>
	<TD align=center>&nbsp;$score&nbsp;</TD>
	|;
	#						print qq
	#|<TR><TD></TD><TD bgcolor="$bgColor" class="rBorder" align=right>&nbsp;$pepCount&nbsp;</TD><TH bgcolor="$bgColor" class="font11" align=left nowrap>$seqVarMod&nbsp;</TH><TD bgcolor="$bgColor" >$mrObs</TD><TD bgcolor="$bgColor" align=center >$pepInfo[2]</TD><TD bgcolor="$bgColor" align=center>$charge<SUP>+</SUP></TD>
	#<TD bgcolor="$bgColor" align=center>&nbsp;$score&nbsp;</TD>
	#|;
							#if ($xicSoftCode eq 'PKV' || $xicSoftCode eq 'SKY' || $xicSoftCode eq 'OS') { #}
							if ($xicSoftCode eq 'MQ') {
								my $rt=($peptideRT{$pepID})? (sprintf "%.2f",$peptideRT{$pepID})*1 : '-';
								print "<TD align=center>&nbsp;$rt&nbsp;</TD>";
							}
							if ($peptideFragments{$pepID}) {
								my $fragMZList;
								my $nbFragment=scalar keys %{$peptideFragments{$pepID}};
								foreach my $fragID (sort {"".$peptideFragments{$pepID}{$a}[5] cmp "".$peptideFragments{$pepID}{$b}[5] || "".$peptideFragments{$pepID}{$a}[6] cmp "".$peptideFragments{$pepID}{$b}[6]} keys %{$peptideFragments{$pepID}}){
									my ($fragMZ,$fragCharge,$fragRT,$fragType,$fragArea)=@{$peptideFragments{$pepID}{$fragID}};
									if($fragRT){print "<TD align=center><FONT onmouseover=\"popup('<B>M/Z : </B>$fragMZ')\" onmouseout=\"popout()\">&nbsp;$fragType&nbsp;$fragArea\@$fragRT&nbsp;</FONT></TD>";}
									else{print "<TD align=center><FONT onmouseover=\"popup('<B>M/Z : </B>$fragMZ')\" onmouseout=\"popout()\">&nbsp;$fragType&nbsp;$fragArea&nbsp;</FONT></TD>";}
								}
								if ($nbTransition != $nbFragment) {
									while ($nbFragment != $nbTransition) {
										print "<TD align=center>&nbsp;-&nbsp;</TD>";
										$nbFragment++;
									}
								}
							}
							else {
								if ($printMCQ_RT && defined($quantifData{$pepID}{$quantifParamInfo{'RT_BEGIN'}}[1])) {
									print "<TD align=center>&nbsp;".&formatRTinmin($quantifData{$pepID}{$quantifParamInfo{'RT_BEGIN'}}[1])."&nbsp;</TD>";
									print "<TD align=center>&nbsp;".&formatRTinmin($quantifData{$pepID}{$quantifParamInfo{'RT_END'}}[1])."&nbsp;</TD>";
								}
								elsif ($printMCQ_RT) {
									print "<TD align=center>&nbsp;-&nbsp;</TD>";
									print "<TD align=center>&nbsp;-&nbsp;</TD>";
								}
								
								my $qData= $quantifData{$pepID}{$quantifParamInfo{$paramCode}}[1]*1;
								my $qDataLink=($isTrace)? "<A href=\"javascript:void(null)\" onclick=\"ajaxDrawXICData(event, $pepID)\"/>$qData</A>" : $qData;
								print "<TD align=right>&nbsp;$qDataLink&nbsp;</TD>\n";
							}
							print "</TR>\n";
							$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
						}
					}
				}
				print "<TR><TD colspan=$colSpan>&nbsp;</TD></TR></TABLE>\n";
			}
		}
		print "<br/>\n\n";
	}
	
	print "End of list.\n\n";

}
#########################
####>LABELED QUANTIF<####
#########################
else {

	#######################################
	####>Analysis-level quantification<####
	#######################################
	if ($call eq 'ana' || $numAnaUsed==1) {
		printWaitingMsg("Fetching data...", 0);

		#($analysisID)=$dbh->selectrow_array("SELECT ID_ANALYSIS FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID") unless $analysisID; # not defined if call=quanti
		#my ($quantifName)=($call eq 'ana')? $dbh->selectrow_array("SELECT NAME,QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID") : '';

		####>Fetching protein info<####
		my (%selMGTopProt,%trueMGTopProt,%posBeg,%protInfo,%idtoidentifier,%protMG); #,%maxPepMG
		my $sthMG=$dbh->prepare("SELECT IDENTIFIER,P.ID_PROTEIN,MATCH_GROUP,NUM_PEP,ALIAS,PROT_DES,MW,PROT_LENGTH,ORGANISM FROM ANALYSIS_PROTEIN A,PROTEIN P WHERE A.ID_PROTEIN=P.ID_PROTEIN AND ID_ANALYSIS=$analysisID AND VISIBILITY=2"); #,NUM_PEP,NUM_MATCH,SCORE,CONF_LEVEL,PEP_COVERAGE,PEP_SPECIFICITY
		$sthMG->execute;
		my $count=0;
		while (my ($identifier,$protID,$matchGroup,$numPep,@protInfo)=$sthMG->fetchrow_array) { #$vis,
			#if ($vis==2) { # could be manually modified
				@{$selMGTopProt{$matchGroup}}=($protID,@protInfo);
				$protMG{$protID}=$matchGroup;
			#}
			#if (!$maxPepMG{$matchGroup} || $maxPepMG{$matchGroup}<$numPep) { # true best prot be able to fetch all MG peptides
				$trueMGTopProt{$matchGroup}=$protID;
			#	$maxPepMG{$matchGroup}=$numPep;
			#}
			$idtoidentifier{$protID}=$identifier;
			@{$protInfo{$identifier}}=($protID,@protInfo);
			$count++;
			if ($count==1000) {
				$count=0;
				printWaitingMsg(".", 1);
			}
		}
		$sthMG->finish;

		printWaitingMsg("/", 1);

		my (%labelingInfo,%sumValues);

		####>SILAC<####
		if ($labelType eq 'SILAC') {
			my (@channelList,@labelModList);
			#my $maxChanNum=0;
			#my $noLabelChannel=0;
			foreach my $infoStrg (@labelInfo) {
				my ($chanNum,$chanName,$labelStrg)=split(';',$infoStrg); # '1;Light;No label;'   '2;Heavy;Label:13C(6);K'
				last if $chanNum !~ /^\d/;
				push @channelList,$chanNum;
				$labelingInfo{$chanNum}{NAME}=$chanName;
				$sumValues{$chanNum}=0;
				if ($labelStrg=~/^None#/) {
					push @{$labelingInfo{$chanNum}{'PARAMS'}},['None','No label','','',''];
					#$noLabelChannel=$chanNum;
				}
				else {
					foreach my $label (split(/\@/,$labelStrg)) {
						my ($labelModifAlias,$labelModifName,$modifRes,$searchMod)=split(/#/,$label);
						$searchMod=~s/( \([A-Z]+)\)\Z/$1/; # 'Labelxxx (K)' -> 'Labelxxx (K' : To allow AA position-independent match
						push @{$labelingInfo{$chanNum}{'PARAMS'}},[$labelModifAlias,$labelModifName,$modifRes,$searchMod,quotemeta($searchMod)];
						push @labelModList,(quotemeta($labelModifName),quotemeta($searchMod)); # record both for compatibility pre/post table MODIFICATION
					}
				}
				#$maxChanNum=$chanNum if $maxChanNum<$chanNum;
			}

			###>Fetching quantification data<###
			#my $sthQP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,ID_PEPTIDE,QUANTIF_VALUE,TARGET_POS FROM PEPTIDE_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID"); # only XIC_AREA is recorded
			#$sthQP->execute;
			$count=0;
			#while (my ($paramID,$pepID,$paramValue,$chanNum)=$sthQP->fetchrow_array) { #}
			foreach my $chanNum (keys %labelingInfo) {
				open (QUANTI,"$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/peptide_quantification_$chanNum.txt");
				while (<QUANTI>) {
					next if $.==1;
					chomp;
					my ($paramID,$pepID,$paramValue)=split(/\t/,$_);
					@{$quantifData{$pepID}{$paramID}}=($chanNum,$paramValue);
					$sumValues{$chanNum}+=$paramValue if $paramValue;
					$count++;
					if ($count==5000) {
						$count=0;
						&printWaitingMsg(".", 1);
					}
				}
				close QUANTI;
			}
			#$sthQP->finish;
#print "</CENTER><TABLE border=1>\n";
			##>Fetching peptide info & grouping isoforms (+/- labeled)
			my (%peptideScore,%sequenceBeg,%labeledPeptideSets,%peptideSets,%peptideMrobs,%pepAll,%pepDataSource,%quantSetVarMods);
			my $extraString=($xicSoftCode eq 'MCQ')? "MCQSET_$selQuantifID=" : "QSET="; # If several quantitation were performed (PD or MassChroQ), some peptides would not be retieved
			my $paramCode='ISO_AREA';
			#my $sthPI=$dbh->prepare("SELECT P.ID_PEPTIDE,PEP_SEQ,MR_OBS,CHARGE,ABS(PEP_BEG),SCORE,DATA FROM PEPTIDE_PROTEIN_ATTRIB PPA,PEPTIDE P WHERE PPA.ID_PEPTIDE=P.ID_PEPTIDE AND P.ID_ANALYSIS=$analysisID AND ID_PROTEIN=? AND DATA LIKE '%$extraString%' ORDER BY PEP_SEQ,CHARGE"); #
			my $sthPI=$dbh->prepare("SELECT ID_PROTEIN,P.ID_PEPTIDE,PEP_SEQ,GROUP_CONCAT(PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),MR_OBS,CHARGE,ABS(PEP_BEG),SCORE,DATA
											FROM PEPTIDE P
											LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
											INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
											WHERE P.ID_ANALYSIS=$analysisID GROUP BY P.ID_PEPTIDE ORDER BY PEP_SEQ,CHARGE");
			$sthPI->execute;
			while (my ($protID,$pepID,$pepSeq,$modCode,$mrObs,$charge,$beg,$score,$pepData)=$sthPI->fetchrow_array) {
				next unless $protMG{$protID};
				next unless $quantifData{$pepID}{$quantifParamInfo{$paramCode}}; # # No quantif data ***qSetStrg may still be defined (incomplete validated set => no quanti data at all)***
				my $matchGroup=$protMG{$protID};
				$pepData='' unless $pepData;
				my ($srcRk)=($pepData=~/SOURCE_RANK=(\d+)/); $srcRk=0 unless $srcRk;
				$pepDataSource{$pepID}=$sourceFiles{$analysisID}{$srcRk};
				my ($quantSetID)=$pepData=~/$extraString(\d+)/;
				#my $varModStrg=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$analysisID,$pepSeq);
				my $varModStrg=($modCode)? ' + '.&promsMod::decodeVarMod($dbh, $pepSeq, $modCode) : '';
#print "**$pepID ($quantSetID: $pepData): $pepSeq//$varModStrg<BR>\n" if $pepSeq eq 'LMVPLLK'; # == 219154; # if (qSetStrg && !$quantifData{$pepID});
#print "<TR><TD>$pepID</TD><TD>$pepSeq</TD><TD>$varModStrg</TD><TD>$charge+</TD><TD>";
				$peptideScore{$pepID}=($score)? $score : 0; # ghost peptides don't have scores
				$peptideMrobs{$pepID}=$mrObs;
				foreach my $qLabelMod (@labelModList) { # remove label mods from vmod string
					$varModStrg=~s/ \+ $qLabelMod \([^\(]+\)(\Z| \+)/$1/;
				}
				@{$quantSetVarMods{$quantSetID}}=("$pepSeq$varModStrg",$matchGroup); # unlabeled reference peptide seqVmod
#print "$pepSeq$varModStrg<BR>\n" if $pepSeq eq 'LMVPLLK'; # $pepID == 219154;
				#if ($quantifData{$pepID}[0] == $noLabelChannel) { # chanNum
				#	$quantSetVarMods{$quantSetID}="$pepSeq$varModStrg";
				#}

				push @{$peptideSets{$quantSetID}},[$pepID,$charge,$beg,$pepSeq,$varModStrg];

				my @beg=split(/,/,$beg);
				@{$posBeg{"$pepSeq$varModStrg"}}=(\@beg,$mrObs);

				$count++;
				if ($count==5000) {
					$count=0;
					printWaitingMsg(".", 1);
				}
			}
			$sthPI->finish;

			foreach my $quantSetID (keys %quantSetVarMods) {
				my ($seqVarMod,$matchGroup)=@{$quantSetVarMods{$quantSetID}};
				foreach my $refPep (@{$peptideSets{$quantSetID}}) {
					my ($pepID,$charge,$beg,$pepSeq,$varModStrg)=@{$refPep};
					$sequenceBeg{$seqVarMod}=$beg if (!$sequenceBeg{$seqVarMod} || $sequenceBeg{$seqVarMod} > $beg); # in case of sequence repetition
					my $chanNum=$quantifData{$pepID}{$quantifParamInfo{$paramCode}}[0];
					#my $srcRk=$quantifData{$pepID}{$quantifParamInfo{$paramCode}}[2];
					push @{$labeledPeptideSets{$matchGroup}{$seqVarMod}{$charge}{$pepDataSource{$pepID}}{$chanNum}},$pepID; # {chanNum} !!!multiple instances of same peptide possible!!!
					push @{$pepAll{$pepID}},[$pepID,$charge,$beg,$pepSeq,$varModStrg];
					$count++;
					if ($count==5000) {
						$count=0;
						printWaitingMsg(".", 1);
					}
				}
			}

			###>Displaying summary<###
			print qq
|<SCRIPT type="text/javascript">document.getElementById('waitDIV').style.display='none'</SCRIPT>
| unless $doExport;
			#print "<FONT class=\"title2\">$quantifName</FONT><BR>\n" if $call eq 'ana';
			&displayPeptideSummary($labelType,$xicSoftCode,\%labelingInfo,\%sumValues);
			if ($doExport) {
				my %exportParameters=('refQuantifData'=>\%quantifData,
								'refPepDataSource'=>\%pepDataSource,
								'refProtInfo'=>\%protInfo,
								'refProtMatch'=>\%selMGTopProt,
								'refPosBeg'=>\%posBeg,
								'reflabeledPeptideSets'=>\%labeledPeptideSets,
								'refPepMrObs'=>\%peptideMrobs,
								'refPeptideScore'=>\%peptideScore,
								'refPeptideSets'=>\%peptideSets,
								'refChannelList'=>\@channelList,
								'refLabelingInfo'=>\%labelingInfo,
								'refPepAll'=>\%pepAll);
				&exportProteinList(\%exportParameters);
				$workbook->close();
				exit;
			}

			###>Displaying data<###
			my $numDataSrc=scalar keys %{$sourceFiles{$analysisID}}; #dataSources;
			#my $colSpan=7 + $maxChanNum;
			my $colSpan=7 + scalar @channelList;
			my $lastChannel=$channelList[-1];
			$colSpan++ if $numDataSrc > 1;
			print "<TABLE border=0 cellspacing=0 cellpadding=2>\n";

			foreach my $matchGroup (sort{$a<=>$b} keys %selMGTopProt) { #ALIAS,PROT_DES,MW,PROT_LENGTH,ORGANISM,NUM_PEP,NUM_MATCH,SCORE,CONF_LEVEL,PEP_COVERAGE,PEP_SPECIFICITY
				my ($protID,$alias,$des,$mw,$length,$org)=@{$selMGTopProt{$matchGroup}};
				print qq
|<TR bgcolor="$darkColor"><TD class="bBorder" colspan=$colSpan><TABLE>
	<TR><TH valign=top><A href="javascript:sequenceView($protID,$analysisID)">$alias</A>:</TH><TD bgcolor="$lightColor" width=100%>$des <FONT class="org">$org</FONT> ($length aa)</TD></TR>
	</TABLE></TD></TR>
<TR><TH>&nbsp;&nbsp;</TH><TH bgcolor="$darkColor" class="rbBorder">#</TH>
<TH bgcolor="$darkColor" class="rbBorder">Peptide sets&nbsp;&nbsp;</TH>
<TH bgcolor="$darkColor" class="rbBorder">&nbsp;Start&nbsp;</TH>
<TH bgcolor="$darkColor" class="rbBorder">&nbsp;Charge&nbsp;</TH>
<TH bgcolor=\"$darkColor\" class=\"rbBorder\">&nbsp;Scores&nbsp;</TH>
|;
				foreach my $chanNum (@channelList) { #1..$maxChanNum
					my $tdClass=($chanNum==$lastChannel && $numDataSrc <= 1)? 'bBorder' : 'rbBorder'; # $maxChanNum
					print "<TH bgcolor=\"$darkColor\" class=\"$tdClass\">&nbsp;$labelingInfo{$chanNum}{NAME}&nbsp;</TH>\n";
				}
				print "<TH bgcolor=\"$darkColor\" class=\"bBorder\">&nbsp;Source&nbsp;</TH>" if $numDataSrc > 1;
				print "<TD width=50%></TD></TR>\n";
				my $bgColor=$lightColor;
				my $numPep=0;
				#foreach my $seqVarMod (sort{$sequenceBeg{$a}<=>$sequenceBeg{$b} || $a cmp $b} keys %{$peptideSets{$matchGroup}}) { #}
				foreach my $seqVarMod (sort{lc($a) cmp lc($b) || $a cmp $b} keys %{$labeledPeptideSets{$matchGroup}}) {
					foreach my $charge (sort keys %{$labeledPeptideSets{$matchGroup}{$seqVarMod}}) {
						#if ($labeledPeptideSets{$matchGroup} && $labeledPeptideSets{$matchGroup}{$seqVarMod} && $labeledPeptideSets{$matchGroup}{$seqVarMod}{$charge}) {
							foreach my $dataSrc (sort keys %{$labeledPeptideSets{$matchGroup}{$seqVarMod}{$charge}}) {
								my (@scores,%quantifValues);
								my $numQuantChannels=0;
								foreach my $chanNum (@channelList) { #1..$maxChanNum
									my $sc='-';
									if ($labeledPeptideSets{$matchGroup}{$seqVarMod}{$charge}{$dataSrc}{$chanNum}) {
										my $pepID=(sort{$peptideScore{$b}<=>$peptideScore{$a}} @{$labeledPeptideSets{$matchGroup}{$seqVarMod}{$charge}{$dataSrc}{$chanNum}})[0];
										$sc=($peptideScore{$pepID})? sprintf "%.2f",$peptideScore{$pepID} : '-'; # ghost peptides don't have scores
										if ($quantifData{$pepID}{$quantifParamInfo{$paramCode}}[1]) {
											$quantifValues{$chanNum}=sprintf "%.1f",$quantifData{$pepID}{$quantifParamInfo{$paramCode}}[1];
											$numQuantChannels++;
										}
									}
									push @scores,$sc;
									#$quantifValues{$chanNum}=$qVal;
								}
#next unless $numQuantChannels==scalar @channelList; # $maxChanNum;
								$numPep++;
								my $startPos=$sequenceBeg{$seqVarMod} || '-';
								my $scoreStrg=join('/',@scores);
								print "<TR><TD></TD><TD bgcolor=\"$bgColor\" class=\"rBorder\" align=right>&nbsp;$numPep&nbsp;</TD><TH bgcolor=\"$bgColor\" class=\"font11\" align=left nowrap>$seqVarMod&nbsp;</TH><TD bgcolor=\"$bgColor\" align=center>$startPos</TD><TD bgcolor=\"$bgColor\" align=center>$charge<SUP>+</SUP></TD>";
								print "<TD bgcolor=\"$bgColor\" align=center>&nbsp;$scoreStrg&nbsp;</TD>";
								foreach my $chanNum (@channelList) { # 1..$maxChanNum
									my $value=(defined $quantifValues{$chanNum})? sprintf "%.1f",$quantifValues{$chanNum} : '-';
									print "<TD bgcolor=\"$bgColor\" align=center>&nbsp;$value&nbsp;</TD>";
								}
								print "<TD bgcolor=\"$bgColor\" align=center nowrap>&nbsp;$dataSrc&nbsp;</TD>" if $numDataSrc > 1;
								print "<TD></TD></TR>\n";
								$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
							}
						#}
						#else { # no quantification data at all
						#	$numPep++;
						#	print "<TR><TD></TD><TD bgcolor=\"$bgColor\" class=\"rBorder\" align=right>&nbsp;$numPep&nbsp;</TD><TH bgcolor=\"$bgColor\" class=\"font11\" align=left nowrap>$seqVarMod&nbsp;</TH><TD bgcolor=\"$bgColor\" align=center>$sequenceBeg{$seqVarMod}</TD><TD bgcolor=\"$bgColor\" align=center>$charge<SUP>+</SUP></TD>";
						#	my $pepID=(sort{$peptideScore{$b}<=>$peptideScore{$a}} @{$peptideSets{$matchGroup}{$seqVarMod}{$charge}})[0];
						#	print "<TD bgcolor=\"$bgColor\" align=center>&nbsp;$peptideScore{$pepID}&nbsp;</TD>";
						#	foreach my $chanNum (1..$maxChanNum) {
						#		print "<TD bgcolor=\"$bgColor\" align=center>&nbsp;-&nbsp;</TD>";
						#	}
						#	print "<TD bgcolor=\"$bgColor\" align=center>-</TD>" if $numDataSrc > 1;
						#	print "<TD></TD></TR>\n";
						#	$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
						#}
					}
				}
				print "<TR><TD colspan=$colSpan>&nbsp;</TD></TR>\n";
#last;
			}
			print "<TR><TD colspan=$colSpan><B>End of list.</B></TD></TR>\n</TABLE>\n";
		}

		####>iTRAQ or TMT<####
		elsif ($labelType=~/ITRAQ|TMT/) {
			my $maxReporterPos=0;
			foreach my $infoStrg (@labelInfo) {
				my ($repPos,$reporter,$monoMass)=split(';',$infoStrg);
				last if $repPos !~ /^\d/;
				@{$labelingInfo{$repPos}}=($reporter,$monoMass);
				$sumValues{$repPos}=0;
				$maxReporterPos=$repPos if $maxReporterPos < $repPos;
			}

			##>Fetching quantification data
			my (%quantifValues,%labeledPeptideSets,%peptideMrobs,%peptideScore,%peptideSets);
			my ($repValue,$signalParamID); # can be REP_INTENSITY or REP_AREA
			my ($intensityParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER P,QUANTIFICATION_METHOD M WHERE P.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND P.CODE='REP_INTENSITY' AND M.CODE='$labelType'");
			my ($areaParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER P,QUANTIFICATION_METHOD M WHERE P.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND P.CODE='REP_AREA' AND M.CODE='$labelType'");
			my %paramID2Code=($intensityParamID=>'REP_INTENSITY',$areaParamID=>'REP_AREA');			
			$count=0;
			foreach my $chanNum (keys %labelingInfo) {
				open (QUANTI,"$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/peptide_quantification_$chanNum.txt");
				while (<QUANTI>) {
					next if $.==1;
					chomp;
					my ($paramID,$pepID,$paramValue)=split(/\t/,$_);
					unless ($signalParamID) {
						next unless $paramID2Code{$paramID};
						$signalParamID=$paramID;
						$repValue=$paramID2Code{$paramID};
					}
					next unless $paramID == $signalParamID; # For TMT, for each reporter, there are two values for each peptide (REP_INTENSITY AND REP_MASS) 
					push @{$quantifValues{$pepID}{$repValue}{$chanNum}},($paramValue);
					$sumValues{$chanNum}+=$paramValue if $paramValue;
					if ($count==5000) {
						$count=0;
						&printWaitingMsg(".", 1);
					}
				}
				close QUANTI;
			}

			##>Fetching peptide info & matchgroups
			my (%peptideData, %pepDataSource);
			#my $sthPI=$dbh->prepare("SELECT P.ID_PEPTIDE,PEP_SEQ,MR_OBS,CHARGE,PEP_BEG,SCORE,DATA FROM PEPTIDE_PROTEIN_ATTRIB PPA,PEPTIDE P WHERE PPA.ID_PEPTIDE=P.ID_PEPTIDE AND P.ID_ANALYSIS=$analysisID AND ID_PROTEIN=?");
			my $sthPI=$dbh->prepare("SELECT P.ID_PEPTIDE,PEP_SEQ,GROUP_CONCAT(PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),MR_OBS,CHARGE,ABS(PEP_BEG),SCORE,DATA
											FROM PEPTIDE P
											LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
											INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
											WHERE P.ID_ANALYSIS=$analysisID AND ID_PROTEIN=? GROUP BY P.ID_PEPTIDE ORDER BY PEP_SEQ,CHARGE");
			foreach my $matchGroup (keys %trueMGTopProt) {
				$sthPI->execute($trueMGTopProt{$matchGroup}); #$protID
				while (my ($pepID,$pepSeq,$modCode,$mrObs,$charge,$beg,$score,$pepData)=$sthPI->fetchrow_array) {
					next unless $score; # ghost peptides don't have scores
					$pepData='' unless $pepData;
					my ($srcRk)=($pepData=~/SOURCE_RANK=(\d+)/); $srcRk=0 unless $srcRk;
					$pepDataSource{$pepID}=$sourceFiles{$analysisID}{$srcRk}; # used for export only
					#my $varModStrg=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$analysisID,$pepSeq);
					my $varModStrg=($modCode)? ' + '.&promsMod::decodeVarMod($dbh, $pepSeq, $modCode) : '';
					if ($peptideData{$matchGroup} && $peptideData{$matchGroup}{$pepID}) { # possible repeats in protein sequence
						$peptideData{$matchGroup}{$pepID}[1]=$beg if $peptideData{$matchGroup}{$pepID}[1] > $beg;
					}
					else {@{$peptideData{$matchGroup}{$pepID}}=("$pepSeq$varModStrg",$beg,$charge,$score);} # first time seen
					push @{$labeledPeptideSets{$matchGroup}{"$pepSeq$varModStrg"}{$charge}},$pepID;
					$peptideMrobs{$pepID}=$mrObs;
					$peptideScore{$pepID}=$score;
					push @{$peptideSets{$pepID}},[$pepID,$charge,$beg,$pepSeq,$varModStrg];
					push @{$labeledPeptideSets{$matchGroup}{"$pepSeq$varModStrg"}{$charge}},$pepID;
					my @beg=split(/,/,$beg);
					@{$posBeg{"$pepSeq$varModStrg"}}=(\@beg,$mrObs);
				}
			}
			$sthPI->finish;

			##>Displaying summary
			print qq
|<SCRIPT LANGUAGE="Javascript">document.getElementById('waitDIV').style.display='none'</SCRIPT>
| unless $doExport;
			&displayPeptideSummary($labelType,$xicSoftCode,\%labelingInfo,\%sumValues);
			if ($doExport) {
				my %exportParameters=('refQuantifData'=>\%quantifValues,
								'refPepDataSource'=>\%pepDataSource,
								'refProtInfo'=>\%protInfo,
								'refProtMatch'=>\%selMGTopProt,
								'refPosBeg'=>\%posBeg,
								'reflabeledPeptideSets'=>\%labeledPeptideSets,
								'refPepMrObs'=>\%peptideMrobs,
								'refPeptideScore'=>\%peptideScore,
								'refPeptideSets'=>\%peptideSets,
								'refPeptideData'=>\%peptideData,
								'refLabelingInfo'=>\%labelingInfo);
				&exportProteinList(\%exportParameters);
				$workbook->close();
				exit;
			}

			##>Displaying data
			my $colSpan=7 + $maxReporterPos; # last col is empty for longer prot title
			print "<TABLE border=0 cellspacing=0 cellpadding=2>\n";

			foreach my $matchGroup (sort{$a<=>$b} keys %selMGTopProt) { #ALIAS,PROT_DES,MW,PROT_LENGTH,ORGANISM //,NUM_PEP,NUM_MATCH,SCORE,CONF_LEVEL,PEP_COVERAGE,PEP_SPECIFICITY
				my ($protID,$alias,$des,$mw,$length,$org)=@{$selMGTopProt{$matchGroup}};
				print qq
|<TR bgcolor="$darkColor"><TD class="bBorder" colspan=$colSpan><TABLE>
	<TR><TH valign=top><A href="javascript:sequenceView($protID,$analysisID)">$alias</A>:</TH><TD bgcolor="$lightColor" width=100%>$des <FONT class="org">$org</FONT> ($length aa)</TD></TR>
	</TABLE></TD></TR>
<TR><TH>&nbsp;&nbsp;</TH><TH bgcolor="$darkColor" class="rbBorder">#</TH>
<TH bgcolor="$darkColor" class="rbBorder">Peptide&nbsp;&nbsp;</TH>
<TH bgcolor="$darkColor" class="rbBorder">&nbsp;Start&nbsp;</TH>
<TH bgcolor="$darkColor" class="rbBorder">&nbsp;Charge&nbsp;</TH>
<TH bgcolor="$darkColor" class="rbBorder">&nbsp;Scores&nbsp;</TH>
|;
				foreach my $repPos (1..$maxReporterPos) {
					my $tdClass=($repPos==$maxReporterPos)? 'bBorder' : 'rbBorder';
					print "<TH bgcolor=\"$darkColor\" class=\"$tdClass\">&nbsp;$labelingInfo{$repPos}[0]&nbsp;</TH>\n";
				}
				print "<TD width=50%></TD></TR>\n";
				my $bgColor=$lightColor;
				my $numPep=0;
				foreach my $pepID (sort{$peptideData{$matchGroup}{$a}[1]<=>$peptideData{$matchGroup}{$b}[1] || $peptideData{$matchGroup}{$a}[0] cmp $peptideData{$matchGroup}{$b}[0] || $peptideData{$matchGroup}{$a}[2]<=>$peptideData{$matchGroup}{$b}[2] || $peptideData{$matchGroup}{$b}[3]<=>$peptideData{$matchGroup}{$a}[3]} keys %{$peptideData{$matchGroup}}) {
					$numPep++;
					my $startPos=$peptideData{$matchGroup}{$pepID}[1] || '-';
					print qq
|<TR>
	<TD></TD>
	<TD bgcolor="$bgColor" class="rBorder" align=right>&nbsp;$numPep&nbsp;</TD>
	<TH bgcolor="$bgColor" class="font11" align=left nowrap>$peptideData{$matchGroup}{$pepID}[0]&nbsp;</TH>
	<TD bgcolor="$bgColor" align=center>$startPos</TD>
	<TD bgcolor="$bgColor" align=center>$peptideData{$matchGroup}{$pepID}[2]<SUP>+</SUP></TD>
	<TD bgcolor="$bgColor" align=center>$peptideData{$matchGroup}{$pepID}[3]</TD>
|;
					if ($quantifValues{$pepID}) {
						foreach my $repPos (1..$maxReporterPos) {
							my $value=($quantifValues{$pepID}{$repValue}{$repPos}[0])? sprintf "%.1f",$quantifValues{$pepID}{$repValue}{$repPos}[0] : '-';
							print "<TD bgcolor=\"$bgColor\" align=center>&nbsp;$value&nbsp;</TD>";
						}
					}
					else {
						foreach my $repPos (1..$maxReporterPos) {
							print "<TD bgcolor=\"$bgColor\" align=center>&nbsp;-&nbsp;</TD>";
						}
					}
					print "<TD></TD></TR>\n";
					$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
				}
				print "<TR><TD colspan=$colSpan>&nbsp;</TD></TR>\n";
				#last  if $matchGroup >= 4;
			}
			print "<TR><TD colspan=$colSpan><B>End of list.</B></TD></TR>\n</TABLE>\n";
		}
#	print qq
#|<DIV id="graphicDiv" class="popup"> <!--filter:alpha(opacity=80);opacity:0.8;-->
#	<DIV id="dashboardDiv"><TABLE border=0>
#		<TR><TD align=right><INPUT type="button" onclick="document.getElementById('graphicDiv').style.display='none'" value=" Close "/></TD><TH><DIV id="titleDiv"></DIV></TH></TR>
#		<TR><TD align=right><DIV id="sliderDiv"></DIV></TD><TD><DIV id="chartDiv"></DIV></TD></TR>
#	</TABLE></DIV>
#</DIV>
#|;
	}

	############################################
	####>MultiAnalysis-level quantification<####
	############################################
	else { # call=quanti
		# TODO (for Maxquant)
		# Implement multi-analysis display similar to Label-free
		# ...
		print qq
|<SCRIPT LANGUAGE="Javascript">document.getElementById('waitDIV').style.display='none';</SCRIPT>
<BR><FONT class="title3">Please, select one analysis from the navigation window and then click on 'Internal Quantification'.</FONT>
| unless $doExport;
	}
}

$dbh->disconnect;

print qq |
		<DIV id="displayDiv" class="popup"> <!--filter:alpha(opacity=80);opacity:0.8;-->
			<DIV id="infoDiv"></DIV>
		</DIV>
		
		<DIV id="divDescription" class="clDescriptionCont"></DIV>
		
		<SCRIPT LANGUAGE="javascript">
			setPopup();
		</SCRIPT>
	</BODY>
</HTML>
|;

#############################################<<<SUBROUTINES>>>###########################################

sub printWaitingMsg {
	my ($msg, $append) = @_;
	my $appendHTML = ($append) ? "+=" : "=";
	print qq |
		<SCRIPT LANGUAGE="JavaScript">document.getElementById('waitSPAN').innerHTML $appendHTML '$msg';</SCRIPT>
	| unless $doExport;	
}

sub fetchQuantiAnalysis {
	my ($quantiID, $anaID) = @_;
	
	my $sthAnaStr = "SELECT A.ID_ANALYSIS, ID_QUANTIFICATION, NAME
				  FROM ANALYSIS A, ANA_QUANTIFICATION AQ
				  WHERE A.ID_ANALYSIS=AQ.ID_ANALYSIS AND ID_QUANTIFICATION=$quantiID";
	$sthAnaStr .= " AND A.ID_ANALYSIS=$anaID" if $anaID;
	$sthAnaStr .= " ORDER BY NAME ASC";

	my $sthAna  = $dbh->prepare($sthAnaStr);
	my $sthAPep = $dbh->prepare("SELECT ID_PEPTIDE
								 FROM PEPTIDE
								 WHERE ID_ANALYSIS=?");
	
	$sthAna->execute;
	while (my ($anaID,$quantiID,$anaName) = $sthAna->fetchrow_array) {
		#$dataFile=~s/_\d+\.*\d*\.pdm/\.msf/; # in case PD
		$listAna{$anaID} = $anaName;
		push @anaOrder, $anaID;
		&getAnalysisDataSources($anaID, $projectID, \%sourceFiles);
		$sthAPep->execute($anaID);
		if ($quantifType eq 'XIC') {
			while (my ($pepID)=$sthAPep->fetchrow_array) {
				push @{$peptideInAna{$anaID}}, $pepID;
			}
		}
	}
	$sthAna->finish;
	$sthAPep->finish;
}


sub fetchAnaProteins {
	my ($anaID, $shouldLoadData, $protID) = @_;
	my $count = 0;
	my $sthPIStr = "SELECT IDENTIFIER, PR.ID_PROTEIN, MATCH_GROUP, NUM_PEP, ALIAS, PROT_DES, MW, PROT_LENGTH, ORGANISM, NUM_MATCH, SCORE, CONF_LEVEL, PEP_COVERAGE, PEP_SPECIFICITY
					FROM ANALYSIS_PROTEIN A, PROTEIN PR
					WHERE A.ID_PROTEIN=PR.ID_PROTEIN AND ID_ANALYSIS=? AND VISIBILITY=2";
	$sthPIStr .= "  AND PR.ID_PROTEIN=$protID" if $protID;
	
	my $sthPPIStr = "SELECT IDENTIFIER, PR.ID_PROTEIN, MATCH_GROUP, NUM_PEP, ALIAS, PROT_DES, MW, PROT_LENGTH, ORGANISM, NUM_MATCH, A.SCORE, CONF_LEVEL, PEP_COVERAGE, PEP_SPECIFICITY, 
					 P.ID_PEPTIDE,PEP_SEQ,GROUP_CONCAT(PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),MR_OBS,ELUTION_TIME,CHARGE,GROUP_CONCAT(DISTINCT(ABS(PEP_BEG)) ORDER BY ABS(PEP_BEG) SEPARATOR ','),P.SCORE,DATA
					 FROM PEPTIDE P
					 LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
					 INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE
					 INNER JOIN PROTEIN PR ON PR.ID_PROTEIN=PPA.ID_PROTEIN 
					 INNER JOIN ANALYSIS_PROTEIN A ON A.ID_PROTEIN=PR.ID_PROTEIN
					 WHERE P.ID_ANALYSIS=? AND VISIBILITY=2";
	$sthPPIStr .= "  AND PR.ID_PROTEIN=$protID" if $protID;
	$sthPPIStr .= "  GROUP BY P.ID_PEPTIDE ORDER BY PEP_SEQ,CHARGE";
	
	my $sthPPI = $dbh->prepare($sthPPIStr);
	my $sthPI = $dbh->prepare($sthPIStr);
	
	if(!$shouldLoadData) {
		$sthPI->execute($anaID);
		while (my ($identifier, $protID, $matchGroup, $numPep, $alias, $protDes, $protMW, $protLength, $org, $numMatch, $protScore, $confLevel, $pepCoverage, $pepSpecificity)=$sthPI->fetchrow_array) {
			my @protInfo = ($alias, $protDes, $protMW, $protLength, $org, $numPep, $numMatch, $protScore, $confLevel, $pepCoverage, $pepSpecificity);
			@{$selMGTopProt{$matchGroup}} = ($protID, @protInfo); # could be manually modified
			@{$protInfo{$identifier}} = ($protID, @protInfo);
			$idtoidentifier{$protID} = $identifier;
			$trueMGTopProt{$matchGroup} = $protID;
		}
		$sthPI->finish;
	} else {
		$sthPPI->execute($anaID);
		while (my ($identifier, $protID, $matchGroup, $numPep, $alias, $protDes, $protMW, $protLength, $org, $numMatch, $protScore, $confLevel, $pepCoverage, $pepSpecificity, $pepID, $pepSeq, $modCode, $mrObs, $rtSc, $charge, $beg, $score, $pepData)=$sthPPI->fetchrow_array) { # $vis,
			my @protInfo = ($alias, $protDes, $protMW, $protLength, $org, $numPep, $numMatch, $protScore, $confLevel, $pepCoverage, $pepSpecificity);
			@{$selMGTopProt{$matchGroup}} = ($protID, @protInfo); # could be manually modified
			@{$protInfo{$identifier}} = ($protID, @protInfo);
			$idtoidentifier{$protID} = $identifier;
			
			### Fetching current protein peptides
			$pepData = '' unless $pepData;
			my $srcRk = ($pepData=~/SOURCE_RANK=(\d+)/) ? $1 : 0;
			$pepDataSource{$pepID} = $sourceFiles{$anaID}{$srcRk};
			
			#my $varModStrg=&promsMod::toStringVariableModifications($dbh,'PEPTIDE',$pepID,$anaID,$pepSeq);
			my $varModStrg = ($modCode) ? ' + '.&promsMod::decodeVarMod($dbh, $pepSeq, $modCode) : '';
			@{$pepInfo{$pepID}} = ($pepSeq, $varModStrg, $charge);
			
			my @beg = split(/,/,$beg);
			@{$posBeg{"$pepSeq$varModStrg"}} = (\@beg, $mrObs);
			
			$peptideMrobs{$pepID} = $mrObs;
			$peptideScore{$pepID} = ($score) ? $score : 0; # ghost peptides don't have scores
			
			#(my $rt)=($rtSc=~m/sc/) ? ($rtSc=~/sc.+;et(.+)/) : $rtSc;
			my $rt = (!$rtSc) ? 0 : ($rtSc =~ /et([\d\.-]+)/) ? $1 : ($rtSc =~ /^([\d\.-]+)/) ? $1 : 0; #$rtSc; # force to 0 because unrecognized (PP 18/04/18)
			$peptideRT{$pepID} = $rt;
			
			push @{$peptideSets{$pepID}}, [$pepID, $charge, $beg[0], $pepSeq, $varModStrg];
			push @{$labeledPeptideSets{$matchGroup}{"$pepSeq$varModStrg"}{$charge}}, $pepID;
			
			if ($maxProtMatch{$identifier}) {
				$maxProtMatch{$identifier}{"$pepSeq$varModStrg"} += 1;
			} else {
				$maxProtMatch{$identifier}{"$pepSeq$varModStrg"} = 1;
			}
			
			$count++;
			if ($count == 5000) {
				$count = 0;
				printWaitingMsg('.', 1);
			}
		}
		$sthPPI->finish;
	}
}



sub fetchPeptidesFromFile {
	my($anaID) = @_;
	my $chanNum = undef;
	my $count = 0;
	
	### Fetch peptide data
	open(QUANTI,"$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/peptide_quantification.txt");
	while (<QUANTI>) {
		next if $.==1;
		chomp;
		
		my ($paramID, $pepID, $paramValue) = split(/\t/, $_);
		
		if(!$peptideInAna{$anaID}) { # Internal Quantification
			@{$quantifData{$pepID}{$paramID}} = ($chanNum, $paramValue);
		} else {
			$peptideQuant{$pepID}{$paramID} = $paramValue;
		}
		
		$count++;
		if ($count == 10000) {
			$count = 0;
			printWaitingMsg('.', 1);
		}
	}
	close QUANTI;
	
	if($peptideInAna{$anaID}) { # Multi-analysis quantification
		$count = 0;
		foreach my $pepID (@{$peptideInAna{$anaID}}) {
			next unless $peptideQuant{$pepID};
			next unless $pepInfo{$pepID}; # Ghost peptide assigned to prot with vis<2 in current analysis (happens with multi-ana MassChroQ)
			
			my ($pepSeq, $varModStrg, $charge) = @{$pepInfo{$pepID}};
			push @{$pepAll{"$pepSeq$varModStrg"}{$charge}{$anaID}}, $pepID;
			
			foreach my $paramID (keys %{$peptideQuant{$pepID}}) {
				@{$quantifData{$pepID}{$paramID}} = ($chanNum, $peptideQuant{$pepID}{$paramID}); #$paramValue
			}
			
			$count++;
			if ($count == 10000) {
				$count=0;
				printWaitingMsg('.', 1)
			}
			delete $peptideQuant{$pepID}; # no longer needed
		}
	}
}

sub fetchTDAPeptides {
	my ($anaID, $protID) = @_;
	my ($line, $fragID, $nbTransition);
	my $swathFile = $promsPath{"quantification"}."/project_$projectID/quanti_$selQuantifID/swath_ana_$anaID.txt";
	my $count = 0;
	
	# Has fragments
	if (-e $swathFile) {
		
		### Fetch peptides data
		my $sthQP=$dbh->prepare("SELECT ID_PEPTIDE FROM PEPTIDE P,ANA_QUANTIFICATION AQ WHERE ID_QUANTIFICATION=$selQuantifID AND P.ID_ANALYSIS=AQ.ID_ANALYSIS AND P.ID_ANALYSIS=$anaID");
		$sthQP->execute;
		while (my ($pepID)=$sthQP->fetchrow_array) {
			next unless exists($pepInfo{$pepID}); ### Avoid missing peptides between sthQP and sthPI queries
			my ($pepSeq, $varModStrg, $charge) = @{$pepInfo{$pepID}};

			push @{$pepAll{"$pepSeq$varModStrg"}{$charge}{$anaID}}, $pepID;
			$quantifData{$pepID} = 1;
			
			$count++;
			if ($count == 10000) {
				$count = 0;
				printWaitingMsg('.', 1);
			}
		}
		$sthQP->finish;
		
		### Fetching peptide fragments data
		my $nbCol = 0;
		open(FRAGMENTFILE, "<", $swathFile) or die ("open: $!");
		while (<FRAGMENTFILE>) {
			next if ($_=~m/ID_PEPTIDE/);
			chomp;
			
			my ($pepID,$fragMZ,$fragCharge,$fragIonType,$fragRes,$fragArea,$rt) = split(/!/,$_);
			$fragMZ = ($fragMZ) ? (sprintf "%0.2f", $fragMZ)*1 : '-';
			$fragArea = ($fragArea=~/\d+/) ? (sprintf "%0.2f", $fragArea)*1 : '-';
			
			my $fragType = "$fragIonType$fragRes";
			if ($peptideFragments{$pepID}) {#$fragID = number of fragment per peptide
				my @sortHash = sort {$b <=> $a} keys %{$peptideFragments{$pepID}};
				$fragID = $sortHash[0]+1;
			} else {
				$fragID = 1;
			}
			
			@{$peptideFragments{$pepID}{$fragID}} = ($fragMZ, $fragCharge, $rt, $fragType, $fragArea, $fragIonType, $fragRes);
			$nbCol = $fragID if $nbCol < $fragID;
		}
		close FRAGMENTFILE;
		
		### Recover number of fragment per petide
		my $sthQuantiAnnot=$dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
		$sthQuantiAnnot->execute;
		my $quantifAnnot = $sthQuantiAnnot->fetchrow_array;
		foreach my $annot (split(/::/,$quantifAnnot)){
			if ($annot =~ /NB_TRANSITION_PER_PEPTIDE=(\d+)/) {
				$nbTransition = $1;
			}
		}
		$nbTransition = $nbCol unless $nbTransition;
	} else {
		fetchPeptidesFromFile($anaID);
	}
	
	return ($nbTransition) ? $nbTransition : ($fragID) ? $fragID : 1;
}

sub displayQuantificationValues {
	my ($protID, $nbTransition) = @_;
	my $bgColor = $lightColor;
	my $numPep = 0;
	my $printMCQ_RT = ($xicSoftCode eq 'MCQ' && scalar keys %maxProtMatch <= $maxNumProtRT) ? 1 : 0;
	$nbTransition = 1 if(!$nbTransition);
	
	## Displaying peptide quantification values
	foreach my $seqVarMod (sort{$posBeg{$a}[0][0]<=>$posBeg{$b}[0][0] || lc($a) cmp lc($b) || $a cmp $b} keys %{$maxProtMatch{$idtoidentifier{$protID}}} ){
		foreach my $charge (sort{$a <=> $b} keys %{$pepAll{$seqVarMod}}) {
			
			my ($refBeg,$mrObs) = @{$posBeg{$seqVarMod}};
			my $beg = join(',',@{$refBeg});
			$beg = '-' unless $beg;
			
			my $maxLocalPepIdx = 0;
			foreach my $anaID (@anaOrder) {
				next unless $pepAll{$seqVarMod}{$charge}{$anaID};
				$maxLocalPepIdx = $#{$pepAll{$seqVarMod}{$charge}{$anaID}} if $maxLocalPepIdx < $#{$pepAll{$seqVarMod}{$charge}{$anaID}};
			}

			foreach my $pepIdx (0..$maxLocalPepIdx) { # occurence of SAME peptide ion in current analysis!!!??? (PP 20/04/18)
				$numPep++;
				my @xicValues = (); my $pepID;
				my $anaPepID;
				foreach my $anaID (@anaOrder) {
					if ($pepAll{$seqVarMod}{$charge}{$anaID} && $pepAll{$seqVarMod}{$charge}{$anaID}[$pepIdx]) {
						$pepID = $pepAll{$seqVarMod}{$charge}{$anaID}[$pepIdx]; # if defined($peptideMrobs{$pepAll{$seqVarMod}{$charge}{$anaID}[$pepIdx]});
						
						if ($xicSoftCode =~ /^(MCQ|PD)$/) { 
							my ($rtb,$rte) = (&formatRTinmin($quantifData{$pepID}{$quantifParamInfo{'RT_BEGIN'}}[1]),&formatRTinmin($quantifData{$pepID}{$quantifParamInfo{'RT_END'}}[1]));
							push @xicValues, "<TD align=center>&nbsp;$rtb&nbsp;</TD><TD align=center>&nbsp;$rte&nbsp;</TD>" if ($printMCQ_RT && $rtb);
							push @xicValues, "<TD align=center>&nbsp;-&nbsp;</TD><TD  align=center>&nbsp;-&nbsp;</TD>" if ($printMCQ_RT && !$rtb);
							my $qDataLink = ($isTrace) ? "<A href=\"javascript:void(null)\" onclick=\"ajaxDrawXICData(event,$pepID)\"/>$quantifData{$pepID}{$quantifParamInfo{'XIC_AREA'}}[1]</A>" : $quantifData{$pepID}{$quantifParamInfo{'XIC_AREA'}}[1];
							push @xicValues, "<TD align=right>&nbsp;$qDataLink&nbsp;</TD>";
						}
						elsif ($xicSoftCode eq 'MQ') {
							my $rt = (sprintf "%.2f",$peptideRT{$pepID})*1;
							push @xicValues, "<TD align=center>&nbsp;$rt&nbsp;</TD>";
							push @xicValues, "<TD align=right>&nbsp;$quantifData{$pepID}{$quantifParamInfo{XIC_AREA}}[1]&nbsp;</TD>";
						}
						elsif($xicSoftCode eq 'SKY') { # TDA
							my $rt = (sprintf "%.2f",$peptideRT{$pepID})*1;
							push @xicValues, "<TD align=center>&nbsp;$rt&nbsp;</TD>";
							
							if(!$peptideFragments{$pepID}) { # MS1
								push @xicValues, "<TD align=right>&nbsp;$quantifData{$pepID}{$quantifParamInfo{XIC_AREA}}[1]&nbsp;</TD>";
							}
							$anaPepID .= $anaID.'_'.$pepID.'@';
						} else { # DIA
							my $rt = (sprintf "%.2f",$peptideRT{$pepID})*1;
							push @xicValues, "<TD align=center>&nbsp;$rt&nbsp;</TD>";
							$anaPepID .= $anaID.'_'.$pepID.'@';
						}
						
						if ($peptideFragments{$pepID}) { # MS 2
							my ($fragMZList, $nbFrag);
							my $nbFragment = scalar keys %{$peptideFragments{$pepID}};
							foreach my $fragID (sort {$peptideFragments{$pepID}{$a}[5] cmp $peptideFragments{$pepID}{$b}[5] || $peptideFragments{$pepID}{$a}[6] cmp $peptideFragments{$pepID}{$b}[6]}keys %{$peptideFragments{$pepID}}){
								my ($fragMZ, $fragCharge, $fragRT, $fragType, $fragArea) = @{$peptideFragments{$pepID}{$fragID}};
								
								my $cellHTMLContent = "<TD align=right><FONT onmouseover=\"popup('<B>M/Z : </B>$fragMZ')\" onmouseout=\"popout()\">&nbsp;$fragType&nbsp;$fragArea";
								$cellHTMLContent .= "\@$fragRT" if ($fragRT);
								$cellHTMLContent .= "&nbsp;</FONT></TD>";
								push @xicValues, $cellHTMLContent;
								
								$nbFrag++;
							}
							for (my $i=$nbFrag; $i<$nbTransition; $i++) {
								push @xicValues, "<TD align=right>&nbsp;-&nbsp;</TD>";
							}
						}
					} else {
						# No RT
						push @xicValues, "<TD align=center>&nbsp;-&nbsp;</TD><TD align=center>&nbsp;-&nbsp;</TD>" if $printMCQ_RT;
						push @xicValues, "<TD align=center>&nbsp;-&nbsp;</TD>" if $xicSoftCode eq 'MQ';
						
						# No XIC
						push @xicValues, "<TD align=right>&nbsp;-&nbsp;</TD>";
						
						if ($xicSoftCode=~/^(PKV|SKY|OS)$/ && %peptideFragments) {
							my $i=0;
							while ($i!=$nbTransition){
								push @xicValues, "<TD align=right>&nbsp;-&nbsp;</TD>";
								$i++;
							}
						}
					}
				}
				
				$mrObs=($pepID && $peptideMrobs{$pepID})? (sprintf "%0.2f",$peptideMrobs{$pepID})*1 : "-";
				my $qSeqVarMod = quotemeta($seqVarMod);
				#if($xicSoftCode eq 'SKY' || $xicSoftCode eq 'PKV' || $xicSoftCode eq 'OS'){#}
				if ($xicSoftCode=~/^(PKV|SKY|OS)$/ && %peptideFragments) { # DIA/TDA MS2 values
					print "<TR class=\"list\" bgcolor=\"$bgColor\"><TD bgcolor=\"#FFFFFF\"></TD><TD class=\"rBorder\" align=right>&nbsp;$numPep&nbsp;</TD><TH class=\"font11\" align=left nowrap ><FONT onmouseover=\"popup('Click to display <B>peptide fragment raw data</B>.')\" onmouseout=\"popout()\" ><A href=\"javascript:void(null)\" onclick=\"ajaxShowFragTDA(event,'$anaPepID')\">$seqVarMod&nbsp;</A></FONT></TH><TD align=center>$mrObs</TD><TD align=center>$beg</TD><TD align=center>$charge<SUP>+</SUP></TD>";
				} else{
					print "<TR class=\"list\" bgcolor=\"$bgColor\"><TD bgcolor=\"#FFFFFF\"></TD><TD class=\"rBorder\" align=right>&nbsp;$numPep&nbsp;</TD><TH class=\"font11\" align=left nowrap>$seqVarMod&nbsp;</TH><TD>$mrObs</TD><TD align=center>$beg</TD><TD align=center>$charge<SUP>+</SUP></TD>";
				}
				print join("\n",@xicValues),"</TR>\n";

				#if ($quantifType eq 'SWATH') {
					$bgColor = ($bgColor eq $lightColor) ? $darkColor : $lightColor;
				#}
			}
		}
	}
}

sub displayPeptideSummary {
	my ($labelType,$xicSoftCode,$refInfo,$refSum)=@_;
	my @posList=($refInfo)? sort{$a<=>$b} keys %{$refInfo} : ();
	my $numExtraCol=scalar @posList;
	$numExtraCol=2 if $numExtraCol<=1;
	#my ($lightColor,$darkColor)=&promsConfig::getRowColors;

	####################################
	####<Exporting summary to Excel>####
	####################################
	if ($doExport) {
		my $worksheet1=$workbook->add_worksheet('Summary');
		my $xlsRow=0;
		if ($labelType !~ /ITRAQ|TMT/) {
			$worksheet1->set_column(0,0,50);
			$worksheet1->set_column(1,1,40);
			$worksheet1->set_column(2,2,30);
		}
		else{
			$worksheet1->set_column(0,0,30);
		}
		$worksheet1->merge_range($xlsRow,0,$xlsRow,$numExtraCol,'Peptide Quantification',$itemFormat{'title'});

		##<Label-free>##
		if ($labelType eq 'FREE') {
			##<Parameters>##
			$xlsRow++;
			$worksheet1->merge_range(++$xlsRow,0,$xlsRow,2,'Quantification parameters',$itemFormat{'mergeColHeader'});
			$worksheet1->write_string(++$xlsRow,0,'Software :',$itemFormat{'headerR'});
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$xicSoftware{$xicSoftCode}.$xicSoftVersionStrg,$itemFormat{'mergeColText'});
			$worksheet1->write_string(++$xlsRow,0,'XIC type :',$itemFormat{'headerR'});
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$quantifMethDesc,$itemFormat{'mergeColText'});

			my %quantifParameters;
			foreach my $parameter (split(/::/,$quantifAnnot) ){
				next unless $parameter;
				my ($parameterName,$parameterValue)=split(/=/,$parameter);
				next if $parameterName eq 'LABEL';
				if ($quantifType eq 'SWATH' || $xicSoftCode eq 'SKY') {
                    $parameterName=~s/_/ /g;
					$parameterValue=~s/_/ /g;
					$parameterName=lc($parameterName);
					$parameterName=ucfirst($parameterName);
                }
				if ($parameterValue!~/N\/A/) {
					if ($parameterName eq 'Exclude modified peptides') {
						if ($parameterValue==1) {
							$quantifParameters{$parameterName}='On';
						}
					}
					else{$quantifParameters{$parameterName}=$parameterValue;}
                }
			}
			$worksheet1->write_string(++$xlsRow,0,'XIC quantification Name :',$itemFormat{'headerR'});
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$selQuantifName,$itemFormat{'mergeColText'});
			#<MassChroQ>#
			if ($xicSoftCode eq 'MCQ') {
				$worksheet1->write_string(++$xlsRow, 0, "Raw-data settings :",$itemFormat{'headerR'});
				$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$quantifParameters{'RAWDATA_ACQUISITION'},$itemFormat{'mergeColText'});
				#<Alignment parameters
				if ( $quantifParameters{'REFERENCE'} ){
					my ($extraAlgoString,$extraAlgoParams)=($quantifParameters{'EXTRACTION_ALGO'} eq 'OBI')? ('OBI-Warp',"Align between $quantifParameters{'MZ_ALIGN_RANGE_MIN'} to $quantifParameters{'MZ_ALIGN_RANGE_MAX'} m/z window ") : ('ms2',"Tendency: $quantifParameters{MS2_TENDENCY}\n•Smoothing: $quantifParameters{MS2_SMOUTHING} (MS/MS) and - $quantifParameters{MS1_SMOUTHING} (MS)");
					my ($refName)=$dbh->selectrow_array("SELECT NAME FROM ANALYSIS WHERE ID_ANALYSIS=$quantifParameters{'REFERENCE'}");
					$worksheet1->set_row(++$xlsRow,13*3);
					$worksheet1->write_string($xlsRow, 0, "Alignment settings :",$itemFormat{'headerR'});
					$worksheet1->merge_range($xlsRow,1,$xlsRow,2,"•Alignment algorithm: $extraAlgoString\n•Reference: $refName\n•$extraAlgoParams",$itemFormat{'mergeColText'});
				}

				my $chargeStatesStrg=($quantifParameters{'ALLCHARGESTATE'})? 'All charge states extracted' : 'Validated charge states extracted';
				$worksheet1->write_string(++$xlsRow, 0, "Charge states : ",$itemFormat{'headerR'});
				$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$chargeStatesStrg,$itemFormat{'mergeColText'});
				#<Advanced settings of quantification
				$worksheet1->set_row(++$xlsRow,13*5);
				$worksheet1->write_string($xlsRow, 0, "Quantification settings :",$itemFormat{'headerR'});
				my $xicFilteringStg="\n";
				$xicFilteringStg.="•Anti-Spike: $quantifParameters{'ANTISPIKE'} " if $quantifParameters{'ANTISPIKE'};
				$xicFilteringStg.="- Half-Mediane:  min=$quantifParameters{'MED_MIN'} max=$quantifParameters{'MED_MAX'} " if $quantifParameters{'MED_MIN'};
				$xicFilteringStg.="- Smoothing: $quantifParameters{'SMOOTH'}" if $quantifParameters{'SMOOTH'};
				$worksheet1->merge_range($xlsRow,1,$xlsRow,2,"•Type of XIC: $quantifParameters{'XIC_EXTRACTION_TYPE'}\n•Mass tolerance for XIC extraction: min=$quantifParameters{'MZTOL_MIN'} max=$quantifParameters{'MZTOL_MAX'}\n•Retention-time of XIC: $quantifParameters{'XIC_VAL'}\n•Detection threshold between $quantifParameters{'DT_START'} to $quantifParameters{'DT_STOP'}$xicFilteringStg",$itemFormat{'mergeColText'});
			}
			#<Proteome Discoverer>#
			else {
				foreach my $param (sort{lc($a) cmp lc($b)} keys %quantifParameters) {
					$worksheet1->write_string(++$xlsRow,0,"$param :",$itemFormat{'headerR'});
					$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$quantifParameters{$param},$itemFormat{'mergeColText'});
				}
			}
		}
		##<Isotope summary>##
		else {# label
			$xlsRow++;
			$worksheet1->write_string(++$xlsRow,0,'Software :',$itemFormat{'headerR'});
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$xicSoftware{$xicSoftCode}.$xicSoftVersionStrg,$itemFormat{'mergeColText'});
			$worksheet1->merge_range(++$xlsRow,0,$xlsRow,$numExtraCol,'Isotope extraction parameters',$itemFormat{'mergeColHeader'});
			$worksheet1->write_string(++$xlsRow,0,'Label :',$itemFormat{'headerR'});
			my $dispLabelType=($labelType eq 'ITRAQ')? 'iTRAQ' : $labelType;
			$worksheet1->merge_range($xlsRow,1,$xlsRow,$numExtraCol,$dispLabelType,$itemFormat{'mergeColText'});
			#$worksheet1->write_string(++$xlsRow,0,'Channel :',$itemFormat{'headerR'});
			if ($labelType eq 'SILAC') {
				$worksheet1->write_string(++$xlsRow,0,'Channel :',$itemFormat{'headerR'});
				my $posNum=1;
				foreach my $pos (@posList) {
					$worksheet1->write_string($xlsRow,$posNum,$pos,$itemFormat{'textC'});
					$posNum++;
				}
			}
			$worksheet1->write_string(++$xlsRow,0,'Signal name :',$itemFormat{'headerR'});
			my $posNum=1;
			foreach my $pos (@posList) {
				my $condName=($labelType eq 'SILAC')? $refInfo->{$pos}{NAME} : $refInfo->{$pos}[0];
				$worksheet1->write_string($xlsRow,$posNum,$condName,$itemFormat{'textC'});
				$posNum++;
			}
			if ($labelType eq 'SILAC') {
				my ($posNum,$numLines)=(1,1);
				foreach my $pos (@posList) {
					foreach my $refParam (@{$refInfo->{$pos}{PARAMS}}) {
						$numLines = scalar @{$refInfo->{$pos}{PARAMS}} if scalar @{$refInfo->{$pos}{PARAMS}} > $numLines;
					}
				}
				$worksheet1->set_row(++$xlsRow,13*$numLines);
				$worksheet1->write_string($xlsRow,0,'Signal name :',$itemFormat{'headerR'});
				foreach my $pos (@posList) {
					my $count=0;
					my $isoString="";
					foreach my $refParam (@{$refInfo->{$pos}{PARAMS}}) {
						$count++;
						$isoString.="$refParam->[0]";
						$isoString.=": $refParam->[1]" if $refParam->[0] ne 'None';
						$isoString.=' ['.$refParam->[3].']' if ($refParam->[3] && $refParam->[3] ne $refParam->[1]);
						$isoString.=' ('.$refParam->[2].')' if $refParam->[2];
						$isoString.="\n" if $count < scalar @{$refInfo->{$pos}{PARAMS}};
					}
					$worksheet1->write_string($xlsRow,$posNum,$isoString,$itemFormat{'textC'});
					$posNum++;
				}
				$worksheet1->write_string(++$xlsRow,0,'Total signal :',$itemFormat{'headerR'});
				$posNum=1;
				foreach my $pos (@posList) {
					$worksheet1->write_string($xlsRow,$posNum,$refSum->{$pos},$itemFormat{'textC'});
					$posNum++;
				}
			}
		}
		#$workbook->close();
		#exit;
	}

	####<Quantification summary>####
	if ($labelType ne 'FREE' && !$doExport) { # label
		my $dispLabelType=($labelType eq 'ITRAQ')? 'iTRAQ' : $labelType;
		my $colSpan=scalar @posList;
		print qq
|<TABLE bgcolor="$darkColor">
<TR bgcolor="$lightColor"><TH align=right bgcolor="$darkColor" class="title3">&nbsp;Label :</TH><TD class="title3" colspan="$colSpan">&nbsp;$dispLabelType&nbsp;&nbsp;&nbsp;<INPUT type="button" name="export" value="Export data" onclick="exportQuanti();"/>&nbsp;</TD></TR>
<TR><TH align=right bgcolor="$darkColor">Software :</TH><TD bgcolor="$lightColor" colspan="$colSpan">&nbsp;<B>$xicSoftware{$xicSoftCode}$xicSoftVersionStrg</B>&nbsp;</TD></TR>
|;
		if ($dataCorrection{SRC}) {
			print "<TR><TH align=right bgcolor=\"$darkColor\">Correction :</TH><TD bgcolor=\"$lightColor\" colspan=\"$colSpan\">&nbsp;";
			if ($dataCorrection{SRC}>=1) {print "<B>$dataCorrection{NAME} (Product #$dataCorrection{PRODUCT}, Lot #$dataCorrection{LOT})</B>";}
			else {print "<B>Unknown product</B>";}
			print "&nbsp;</TD></TR>\n";
		}
		if ($labelType eq 'SILAC') {
			print "<TR bgcolor=\"$lightColor\"><TH align=right bgcolor=\"$darkColor\">&nbsp;Channel :</TH>";
			foreach my $pos (@posList) {
				print "<TH>&nbsp;$pos&nbsp;</TH>";
			}
			print "</TR>\n";
		}
		print "<TR bgcolor=\"$lightColor\"><TH align=right bgcolor=\"$darkColor\">&nbsp;Signal name :</TH>";
		foreach my $pos (@posList) {
			my $condName=($labelType eq 'SILAC')? $refInfo->{$pos}{NAME} : $refInfo->{$pos}[0];
			print "<TH>&nbsp;$condName&nbsp;</TH>";
		}
		print "</TR>\n";
		if ($labelType eq 'SILAC') {
			print "<TR bgcolor=\"$lightColor\"><TH align=right bgcolor=\"$darkColor\">&nbsp;Isotope(s) :</TH>";
			foreach my $pos (@posList) {
				print "<TH>&nbsp;";
				my $count=0;
				foreach my $refParam (@{$refInfo->{$pos}{PARAMS}}) {
					$count++;
					print "$refParam->[0]";
					print ": $refParam->[1]" if $refParam->[0] ne 'None';
					print ' [',$refParam->[3],']' if ($refParam->[3] && $refParam->[3] ne $refParam->[1]);
					print ' (',$refParam->[2],')' if $refParam->[2];
					print '<BR>' if $count < scalar @{$refInfo->{$pos}{PARAMS}};
				}
				print '&nbsp;</TH>';
			}
			print "</TR>\n";
		}
		print "<TR bgcolor=\"$lightColor\"><TH align=right bgcolor=\"$darkColor\">&nbsp;Total signal :</TH>";
		foreach my $pos (@posList) {
			printf "<TH>&nbsp;%.2e&nbsp;</TH>",$refSum->{$pos};
		}
		print "</TR>\n</TABLE>\n<BR>\n</CENTER>\n";
	}
	if ($quantifAnnot =~ /EXTRACTION_ALGO/) { # XIC extraction for example -> display XIC parameters of extraction
		# example: '::ANTISPIKE=5::DT_START=10::DT_STOP=5000::EXTRACTION_ALGO=OBI::MED_MAX=20::MED_MIN=5::MZ_ALIGN_RANGE_MAX=1200::MZ_ALIGN_RANGE_MIN=400::MZTOL_MAX=0.3::MZTOL_MIN=0.3::QUANTIF_NAME=Ext. ion chrom. extraction::RAWDATA_ACQUISITION=profile::REFERENCE=1311::SMOOTH=3::XIC_EXTRACTION_TYPE=sum::XIC_VAL=real_or_mean'
		my @quantiList;
		#if ($action eq 'xicView') {
		#	my ($realFocus)=$dbh->selectrow_array("SELECT FOCUS FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
		#	if ($realFocus ne 'peptide') {
		#		my $sthParentQ=$dbh->prepare("SELECT ID_PARENT_QUANTIFICATION FROM PARENT_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
		#		$sthParentQ->execute;
		#		while (my ($idParentQ)= $sthParentQ->fetchrow_array) {
		#			push @quantiList,$idParentQ;
		#		}
		#	}
		#	else {
		#		push @quantiList,$selQuantifID;
		#	}
		#}
		#else{
			push @quantiList,$selQuantifID;
		#}

		foreach my $quantiXIC (@quantiList) {
			#print "<TABLE bgcolor=\"$darkColor\">\n";
			my ($quantifName,$quantifAnnotLocal)=$dbh->selectrow_array("SELECT NAME,QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantiXIC");
			my %quantifParameters;
			foreach my $parameter (split(/::/,$quantifAnnotLocal) ){
				next unless $parameter;
				my ($parameterName,$parameterValue)=split(/=/,$parameter);
				next unless $parameterValue;
				$quantifParameters{$parameterName}=$parameterValue;
			}
			if (!$doExport) {
				print "<BR>\n" if $call eq 'quanti';
				print qq
|<TABLE bgcolor="$darkColor" align=center>
<TR>
<TH align=right nowrap valign=top>XIC quantification Name :</TH>
<TH bgcolor="$lightColor" nowrap align=left>&nbsp;$quantifName&nbsp;|;
				print "&nbsp;&nbsp;<INPUT type=\"button\" name=\"export\" value=\"Export data\" onclick=\"exportQuanti();\"/>&nbsp;" if $labelType eq 'FREE';
				print qq
|</TH></TR>
<TR><TH align=right nowrap>Software :</TH><TD bgcolor="$lightColor">&nbsp;$xicSoftware{$xicSoftCode}$xicSoftVersionStrg&nbsp;</TD></TR>
<TR>
<TH align=right nowrap valign=top>Raw-data settings :</TH>
<TD bgcolor="$lightColor" nowrap>&nbsp;<B>Extraction type:</B>
$quantifParameters{RAWDATA_ACQUISITION}
<FONT class="font11">&nbsp;(for mzXML)</FONT>
</TD>
</TR>
|;
			}
			if ( $quantifParameters{'REFERENCE'} ){
				my $extraAlgoString=($quantifParameters{'EXTRACTION_ALGO'} eq 'OBI')? 'OBI-Warp' : 'ms2';
				my ($refName)=$dbh->selectrow_array("SELECT NAME FROM ANALYSIS WHERE ID_ANALYSIS=$quantifParameters{'REFERENCE'}");
				my $alignParams="<TD nowrap align=\"right\">&nbsp;";
				if($quantifParameters{'EXTRACTION_ALGO'} eq 'MS2'){
					$alignParams.="<B>Tendency:</B> $quantifParameters{MS2_TENDENCY} " if $quantifParameters{MS2_TENDENCY};
					$alignParams.="<B> - Smoothing:</B> $quantifParameters{MS2_SMOUTHING} <B>(MS/MS) " if $quantifParameters{MS2_SMOUTHING};
					$alignParams.=" $quantifParameters{MS1_SMOUTHING} <B>(MS)</B>" if $quantifParameters{MS1_SMOUTHING};
					$alignParams.="</TD>";
				}
				else{
					$alignParams="<TD nowrap>&nbsp;<B>Align between</B> $quantifParameters{'MZ_ALIGN_RANGE_MIN'} <B>to</B> $quantifParameters{'MZ_ALIGN_RANGE_MAX'} <B>m/z window</B></TD>";
				}
				if (!$doExport) {
					print qq
|<TH align=right nowrap valign=top>Alignment settings :</TH><TD bgcolor="$lightColor"><TABLE>
	<TR><TD nowrap align=\"right\">&nbsp;<B>Alignment algorithm:</B>
	$extraAlgoString
	&nbsp;&nbsp;&nbsp;<B>Reference:</B>
	$refName
	</TD>
	</TR>
	<TR>
	$alignParams
	</TR>
</TABLE>
</TD>
</TR>
|;
				}

			}
			my $chargeStatesStrg=($quantifParameters{'ALLCHARGESTATE'})? 'All charge states extracted' : 'Validated charge states extracted';

			if (!$doExport) {
				my $xicTypeStg=($quantifParameters{'XIC_EXTRACTION_TYPE'} eq 'sum')? 'TIC XIC':'BasePeak XIC';
				my $peakMatchingStg=($quantifParameters{'XIC_VAL'} eq 'post_matching')? 'Post matching mode': ($quantifParameters{'XIC_VAL'} eq 'mean')? 'Mean mode' : 'Real or mean mode';
				my $xicRange=($quantifParameters{'XIC_RANGE'} && $quantifParameters{'XIC_RANGE'} eq 'ppm')? 'ppm': 'mz';# Before 28/03/14, every extraction was made on mz tol...
				my $xicFiltering='';
				$xicFiltering.="&nbsp;<B>Anti-Spike:</B>&nbsp;$quantifParameters{'ANTISPIKE'}\n" if $quantifParameters{'ANTISPIKE'};
				$xicFiltering.="&nbsp;&nbsp;&nbsp;<B>Half-Mediane:</B> min=$quantifParameters{'MED_MIN'} max=$quantifParameters{'MED_MAX'}\n" if $quantifParameters{'MED_MAX'};
				$xicFiltering.="&nbsp;&nbsp;&nbsp;<B>Smoothing:</B>&nbsp;$quantifParameters{'SMOOTH'}" if $quantifParameters{'SMOOTH'};
				$xicFiltering='&nbsp;<B>NONE</B>' unless $xicFiltering;
				print qq
|<TR>
<TH align=right nowrap valign=top>Charge states :</TH><TD bgcolor=$lightColor>&nbsp;$chargeStatesStrg</TD></TR>
<TR>
<TH align=right nowrap valign=top>Quantification settings :<BR><INPUT type="button" id="moreSettings" class="font11" value="More settings" onclick="updateSettings('more')"/><INPUT type="button" id="lessSettings" class="font11" value="Less settings" style="display:none" onclick="updateSettings('less')"/>
</TH><TD bgcolor=$lightColor nowrap valign=top>
	<TABLE cellpadding=0 cellspacing=0>
	<TR><TD nowrap>&nbsp;<B>Type of XIC:</B>
	$xicTypeStg
	</TR>
	</TABLE>
	<DIV id="advancedSetDIV" style="display:none">
	&nbsp;&bull;<B><U>Advanced settings:</U></B>
	<TABLE cellpadding=0 cellspacing=0>
	<TR>
	<TD nowrap>&nbsp;<B>Size of mass tolerance window for XIC:</B> min=$quantifParameters{'MZTOL_MIN'} max=$quantifParameters{'MZTOL_MAX'} $xicRange</TD>
	</TR>
	<TR><TD nowrap>&nbsp;<B>Peak matching for XIC:</B>
	$peakMatchingStg
	</TD>
	</TR>
	<TR>
	<TD nowrap>&nbsp;<B>Detection threshold between </B>$quantifParameters{'DT_START'}<B> to </B>$quantifParameters{'DT_STOP'}</TD>
	</TR>
	<TR><TD>&nbsp;&nbsp;<B><U>XIC filtering</U>:</B></TD></TR>
	<TR>
	<TD bgcolor=$lightColor nowrap>$xicFiltering
	</TD>
	</TR>
</TABLE></TR>
</DIV>
</TR>
</TABLE>
<BR>
</CENTER>
|;
			}
			#print "</TR>\n</TABLE>\n<BR>\n</CENTER>\n";
		}
	}
	#elsif ($labelType eq 'FREE' && ($xicSoftCode eq 'PD' || $xicSoftCode eq 'PKV' || $xicSoftCode eq 'SKY' || $xicSoftCode=~/OS/)) { #}
	elsif ($labelType eq 'FREE' && $xicSoftCode=~/^(PKV|SKY|OS|MQ)$/) {
		my %quantifParameters;
		foreach my $parameter (split(/::/,$quantifAnnot)){
			my ($parameterName,$parameterValue)=split(/=/,$parameter);
			next if $parameterName eq 'LABEL' || !$parameterValue;
			if ($xicSoftCode=~/^(OS|PKV|SKY)$/) {
				$parameterName=~s/_/ /g;
				$parameterValue=~s/_/ /g;
				$parameterName=lc($parameterName);
				$parameterName=ucfirst($parameterName);
			}
			if ($parameterValue!~/N\/A/) {
				if ($parameterName eq 'Exclude modified peptides') {
                    if ($parameterValue==1) {
                        $quantifParameters{$parameterName}='On';
                    }
                }
                else{$quantifParameters{$parameterName}=$parameterValue;}
			}
		}
		if (!$doExport) {
			print "<BR>\n" if $call eq 'quanti';
			print qq
|<TABLE bgcolor="$darkColor" align=center>
<TR>
<TH align=right nowrap valign=top>XIC quantification Name :</TH>
<TH bgcolor="$lightColor" nowrap align=left>&nbsp;$selQuantifName&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" name=\"export\" value=\"Export data\" onclick=\"exportQuanti();\"/>&nbsp;</TH>
</TR>
<TR><TH align=right nowrap>Labeling :</TH><TD bgcolor="$lightColor">&nbsp;None&nbsp;</TD></TR>
<TR><TH align=right nowrap>Software :</TH><TD bgcolor="$lightColor">&nbsp;$xicSoftware{$xicSoftCode}$xicSoftVersionStrg&nbsp;</TD></TR>
<TR><TH align=right nowrap valign="top">Parameters :</TH><TD bgcolor="$lightColor">
|;
			if (scalar keys %quantifParameters) {
				print "<TABLE cellpadding=0>\n";
				foreach my $param (sort{lc($a) cmp lc($b)} keys %quantifParameters) {
					print "<TR><TD nowrap align=\"right\">&nbsp;<B>$param:</B></TD><TD>&nbsp;$quantifParameters{$param}&nbsp;</TD></TR>\n";
				}
				print "</TABLE>";
			}
			else {print "&nbsp;None recorded.";}
			print qq
|</TD></TR>
</TABLE>
<BR>
|;
		}
	}
}


sub deleteQuantification {
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1); # writes PERL warnings as HTML comments <!-- --> must be activated after header
	print qq
|<HEAD>
<TITLE>Deleting Quantification</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<BODY background='$promsPath{images}/bgProMS.gif'>
<FONT class="title2"><BR><BR><BR>Deleting quantification...|;

	####<Connect to the database
	my $dbh=&promsConfig::dbConnect;

	my ($projID)=&promsMod::getProjectID($dbh,$selQuantifID,'quantification');
	#my ($status)=$dbh->selectrow_array("SELECT STATUS FROM QUANTIFICATION WHERE ID_QUANTIFCATION=$selQuantifID");

	#<Ghost peptides
	my $sthgetAnaQ=$dbh->prepare("SELECT ID_ANALYSIS FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	#my $sthgetVP=$dbh->prepare("SELECT PQ.ID_PEPTIDE FROM PEPTIDE P,PEPTIDE_QUANTIFICATION PQ WHERE PQ.ID_PEPTIDE=P.ID_PEPTIDE AND ID_QUANTIFICATION=$selQuantifID AND ID_ANALYSIS=? AND VALID_STATUS=0 AND SCORE IS NULL");# If a single extraction is deleted, just the ghost peptides of this analysis will be suppressed
	my $sthgetVP=$dbh->prepare("SELECT ID_PEPTIDE FROM PEPTIDE WHERE ID_ANALYSIS=? AND VALID_STATUS=0 AND SCORE IS NULL");# If a single extraction is deleted, just the ghost peptides of this analysis will be suppressed
	my $sthgetQuanA=$dbh->prepare("SELECT Q.ID_QUANTIFICATION FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND FOCUS='peptide' AND ID_ANALYSIS=?");
	#my $sthgetVPQ=$dbh->prepare("SELECT COUNT(*) FROM PEPTIDE_QUANTIFICATION WHERE ID_PEPTIDE=?");
	my $sthGetProtID=$dbh->prepare("SELECT ID_PROTEIN FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_PEPTIDE=?");
	my $sthGetPPA=$dbh->prepare("SELECT COUNT(*) FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_PROTEIN=? AND ID_ANALYSIS=? AND PEP_BEG>0");
	my $sthDelPPA1=$dbh->prepare("DELETE FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_PEPTIDE=?");
	my $sthDelPPA2=$dbh->prepare("DELETE FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_PROTEIN=? AND ID_ANALYSIS=?");
	#my $sthDelPepQ=$dbh->prepare("DELETE FROM PEPTIDE_QUANTIFICATION WHERE ID_PEPTIDE=?");
	my $sthDelPepM=$dbh->prepare("DELETE FROM PEPTIDE_MODIFICATION WHERE ID_PEPTIDE=?");
	my $sthDelPep=$dbh->prepare("DELETE FROM PEPTIDE WHERE ID_PEPTIDE=?");
	my $sthDelAP=$dbh->prepare("DELETE FROM ANALYSIS_PROTEIN WHERE CONF_LEVEL=0 AND ID_PROTEIN=? AND ID_ANALYSIS=?");
	my $sthGetPData=$dbh->prepare("SELECT ID_PEPTIDE,DATA FROM PEPTIDE WHERE ID_ANALYSIS=? AND DATA LIKE '%##MCQSET_$selQuantifID=%'");
	my $sthUpPData=$dbh->prepare("UPDATE PEPTIDE SET DATA=? WHERE ID_PEPTIDE=?");

	my ($numPepQ,$numPPA);

	my (%virtualPep,%quantifList);
	$sthgetAnaQ->execute;
	while (my ($anaID) = $sthgetAnaQ->fetchrow_array) {
		$sthgetVP->execute($anaID);
		my $hasVirtualPep=0;
		while (my ($vpPepID)=$sthgetVP->fetchrow_array) {
			$virtualPep{$vpPepID}=$anaID;
			$hasVirtualPep=1;
		}
		if ($hasVirtualPep) {
			$sthgetQuanA->execute($anaID);
			while (my ($quantID) = $sthgetQuanA->fetchrow_array) {
				$quantifList{$quantID}=1;
			}
		}
	}

	##<Fetch list of quantified virtual peptides
	my %peptideQuant;
	foreach my $quantID (keys %quantifList) {
		foreach my $quantFile (glob ("$promsPath{quantification}/project_$projID/quanti_$quantID/peptide_quantification*")) {
			open(QUANTI,$quantFile);
			while (<QUANTI>) {
				next if $.==1;
				my $pepID=(split(/\t/,$_))[0];
				next unless $virtualPep{$pepID};
				$peptideQuant{$virtualPep{$pepID}}{$pepID}{$quantID}=1; # {anaID}{pepID}{quantID}
			}
			close QUANTI;
		}
	}

	foreach my $anaID (sort{$b<=>$a} keys %peptideQuant) {
		foreach my $vpPepID (keys %{$peptideQuant{$anaID}}) {
			next if scalar keys %{$peptideQuant{$anaID}{$vpPepID}} > 1;
			next unless $peptideQuant{$anaID}{$vpPepID}{$selQuantifID};
			#$sthgetVPQ->execute($vpPepID);
			#($numPepQ)=$sthgetVPQ->fetchrow_array;
			#next if $numPepQ>1;# This ghost peptide is still associated to another quantification than this one, do not delete it
			# This ghost peptide was associated to that quantification only -> need to be deleted !
			$sthGetProtID->execute($vpPepID);# Get all the proteins associated to that peptide (local to 1 analysis)
			$sthDelPPA1->execute($vpPepID);# Delete all PPA associated to this peptide in this analysis
			#$sthDelPepQ->execute($vpPepID); # Quantification ## Line to be removed after deletion of PEPTIDE_QUANTIFICATION table
			$sthDelPepM->execute($vpPepID); # Modification
			$sthDelPep->execute($vpPepID);# Delete ghost-peptide
			while (my ($protID)=$sthGetProtID->fetchrow_array) {
				$sthGetPPA->execute($protID,$anaID);
				($numPPA)=$sthGetPPA->fetchrow_array;
				# Delete Virtual Proteins
				if ($numPPA==0) {
					$sthDelPPA2->execute($protID,$anaID);
					$sthDelAP->execute($protID,$anaID);
				}
			}
		}
		$sthGetPData->execute($anaID);
		while (my ($pepID,$data)=$sthGetPData->fetchrow_array) { # Update DATA info in PEPTIDE table if quantification of isotope extraction was performed with MassChroQ
			$data=~s/##MCQSET_$selQuantifID=\d+//;
			$sthUpPData->execute($data,$pepID);
		}
	}

	$sthgetVP->finish;
	#$sthgetVPQ->finish;
	$sthGetProtID->finish;
	$sthGetPPA->finish;
	$sthDelPPA1->finish;
	$sthDelPPA2->finish;
	#$sthDelPepQ->finish;
	$sthDelPepM->finish;
	$sthDelPep->finish;
	$sthDelAP->finish;
	$sthGetPData->finish;
	$sthUpPData->finish;

	#<Peptides
	#$dbh->do("DELETE FROM PEPTIDE_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID"); ## Line to be removed after deletion of PEPTIDE_QUANTIFICATION table
	print '.';

	#<Proteins
	$dbh->do("DELETE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	print '.';

	#<Parent & analysis
	$dbh->do("DELETE FROM PARENT_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	$dbh->do("DELETE FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	$dbh->do("DELETE FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	print '.';

	if (-e "$promsPath{quantification}/project_$projID") {
		#remove_tree("$promsPath{quantification}/project_$projID/quanti_$selQuantifID") if -e "$promsPath{quantification}/project_$projID/quanti_$selQuantifID";
		rmtree("$promsPath{quantification}/project_$projID/quanti_$selQuantifID") if -e "$promsPath{quantification}/project_$projID/quanti_$selQuantifID";
		my @remaindDirs = glob "$promsPath{quantification}/project_$projID/quanti_*";
		#remove_tree("$promsPath{quantification}/project_$projID") unless scalar @remaindDirs;
		rmtree("$promsPath{quantification}/project_$projID") unless scalar @remaindDirs;
	}

	# Remove log file of this quantification
	unlink "$promsPath{logs}/quanti_$selQuantifID.log" if -e "$promsPath{logs}/quanti_$selQuantifID.log";

	print " Done</FONT>\n";

	###<Rebuild match groups
	my $numAna=scalar keys %peptideQuant;
	print "<FONT class=\"title2\"><BR><BR>Rebuilding match groups for $numAna Analyses:</FONT>\n";
	my %projectData;
	my %options=(VERBOSE=>1,COMMIT=>0,USE_GHOST=>1);
	my $anaCount=0;
	foreach my $anaID (sort{$b<=>$a} keys %peptideQuant) {
		$anaCount++;
		print "<FONT class=\"title2\"><BR>&nbsp;+Analysis #$anaCount:</FONT><BR>\n";
		&promsMod::updateMatchGroups($dbh,$anaID,\%projectData,\%options);
	}

	$dbh->commit;
	$dbh->disconnect;

	print "<FONT class=\"title2\"><BR>Done.</FONT>\n";

	sleep 2;

	print qq
|<SCRIPT LANGUAGE="JavaScript">
parent.optionFrame.location.reload(); // will update report button if only primary peptide quantif remains
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

sub exportProteinList {
	#my ($refQuantifData,$refPepDataSource,$refProtInfo,$refProtMatch,$reflabeledPeptideSets,$refPepMrObs,$refPeptideScore,$refPeptideSets,$refPeptideFragments,$refPepProtPos,$nbTransition,@extraParams)=@_;
	my ($refExportParameters)=@_;
	my $refQuantifData=$refExportParameters->{refQuantifData};
	my $refPepDataSource=$refExportParameters->{refPepDataSource};
	my $refProtInfo=$refExportParameters->{refProtInfo};
	my $refProtMatch=$refExportParameters->{refProtMatch};
	my $refPosBeg=$refExportParameters->{refPosBeg};
	my $reflabeledPeptideSets=$refExportParameters->{reflabeledPeptideSets};
	my $refPepMrObs=$refExportParameters->{refPepMrObs};
	my $refPeptideScore=$refExportParameters->{refPeptideScore};
	my $refPeptideSets=$refExportParameters->{refPeptideSets};
	my $refPeptideFragments=$refExportParameters->{refPeptideFragments} || undef;
	my $refPepProtPos=$refExportParameters->{refPepProtPos} || undef;
	my $nbTransition=$refExportParameters->{nbTransition} || undef;
	my $refPepAll=$refExportParameters->{refPepAll} || undef;
	my $refLabelingInfo=$refExportParameters->{refLabelingInfo} || undef;
	my $refListAna=$refExportParameters->{refListAna} || undef;
	my $refAnaOrder=$refExportParameters->{refAnaOrder} || undef;
	my $refChannelList=$refExportParameters->{refChannelList} || undef;
	my $refPeptideData=$refExportParameters->{refPeptideData} || undef;
	#my $refPeptideRT=$refExportParameters->{refPeptideRT} || undef;

	# my @extraParams; # Removed on March 2017
	# @extraParams : XIC single -> undef
	# @extraParams : XIC multiple/SWATH -> %pepAll,%listAna,@anaOrder
	# @extraParams : SILAC -> @channelList,%labelingInfo,%pepAll
	# @extraParams : iTRAQ -> %peptideData,%labelingInfo
	$nbTransition=($nbTransition)? $nbTransition : '';
	my $printMCQ_RT=($xicSoftCode eq 'MCQ' && scalar keys %{$refProtMatch} <= $maxNumProtRT)? 1 : 0;

	####<Start printing>####
	my $worksheet2=$workbook->add_worksheet('Results');
	my $xlsRow=0;
	my $xlsCol=0;

	####<Header>####
	$xlsCol=0;
	###> Protein name + desc
	if ($quantifType =~/^(SWATH|TDA|DIA)$/) {
		$xlsRow++;
        $worksheet2->set_column($xlsCol,$xlsCol,20); # col length
		$worksheet2->merge_range(0,$xlsCol,1,$xlsCol,'Identifier',$itemFormat{'mergeRowHeader'}) ;
		$worksheet2->set_column(++$xlsCol,$xlsCol,30); # col length
		$worksheet2->merge_range(0,$xlsCol,1,$xlsCol,'Sequence',$itemFormat{'mergeRowHeader'});
		$worksheet2->set_column(++$xlsCol,$xlsCol,30); # col length
		$worksheet2->merge_range(0,$xlsCol,1,$xlsCol,'PTMs',$itemFormat{'mergeRowHeader'});
		$worksheet2->set_column(++$xlsCol,$xlsCol,30); # col length
		$worksheet2->merge_range(0,$xlsCol,1,$xlsCol,'Mr(Calc)',$itemFormat{'mergeRowHeader'});
		$worksheet2->merge_range(0,++$xlsCol,1,$xlsCol,'Start',$itemFormat{'mergeRowHeader'});
		$worksheet2->merge_range(0,++$xlsCol,1,$xlsCol,'Charge',$itemFormat{'mergeRowHeader'});
		$worksheet2->merge_range(0,++$xlsCol,1,$xlsCol,'Scores',$itemFormat{'mergeRowHeader'});
    }
    else {
		$worksheet2->set_column($xlsCol,$xlsCol,20); # col length
		$worksheet2->write_string($xlsRow,$xlsCol,"Identifier",$itemFormat{'header'});
		$worksheet2->set_column(++$xlsCol,$xlsCol,30); # col length
		$worksheet2->write_string($xlsRow,$xlsCol,"Sequence",$itemFormat{'header'});
		$worksheet2->set_column(++$xlsCol,$xlsCol,30); # col length
		$worksheet2->write_string($xlsRow,$xlsCol,"PTMs",$itemFormat{'header'});
		$worksheet2->write_string($xlsRow,++$xlsCol,"Mr(Obs)",$itemFormat{'header'}) unless $labelType eq 'SILAC';
		$worksheet2->write_string($xlsRow,++$xlsCol,"Start",$itemFormat{'header'});
		$worksheet2->write_string($xlsRow,++$xlsCol,"Charge",$itemFormat{'header'});
		$worksheet2->write_string($xlsRow,++$xlsCol,"Scores",$itemFormat{'header'});
	}

	if ($labelType eq 'FREE') {
		if ($call eq 'quanti') {
			foreach my $anaID (@{$refAnaOrder}) {
				$worksheet2->write_string($xlsRow,++$xlsCol,"RTbeg",$itemFormat{'header'}) if $printMCQ_RT;
				$worksheet2->write_string($xlsRow,++$xlsCol,"RTend",$itemFormat{'header'}) if $printMCQ_RT;
				$worksheet2->write_string($xlsRow,++$xlsCol,$refListAna->{$anaID}[1],$itemFormat{'header'}) unless $quantifType=~/^(SWATH|TDA|DIA)$/;
			}
		}
		else { # call = analysis
			$worksheet2->write_string($xlsRow,++$xlsCol,"RTbeg",$itemFormat{'header'}) if $printMCQ_RT;
			$worksheet2->write_string($xlsRow,++$xlsCol,"RTend",$itemFormat{'header'}) if $printMCQ_RT;
			$worksheet2->write_string($xlsRow,++$xlsCol,"XIC",$itemFormat{'header'}) unless $quantifType=~/^(SWATH|TDA|DIA)$/;
		}
	}
	elsif ($labelType eq 'SILAC') {
		foreach my $chanNum (@{$refChannelList}) {
			$worksheet2->write_string($xlsRow,++$xlsCol,$refLabelingInfo->{$chanNum}{NAME},$itemFormat{'header'});
		}
	}
	elsif ($labelType=~/ITRAQ|TMT/) {
		my $maxReporterPos=scalar keys %{$refLabelingInfo};
		foreach my $repPos (1..$maxReporterPos){
			$worksheet2->write_string($xlsRow,++$xlsCol,$refLabelingInfo->{$repPos}[0],$itemFormat{'header'});
		}
	}
	if ($quantifType !~ /^(SWATH|TDA|DIA)$/){
		$worksheet2->set_column(++$xlsCol,$xlsCol,30); # col length
		$worksheet2->write_string($xlsRow,$xlsCol,"Source",$itemFormat{'header'});
		$worksheet2->write_string($xlsRow,++$xlsCol,"MW (kDa)",$itemFormat{'header'});
		$worksheet2->set_column(++$xlsCol,$xlsCol,80); # col length
		$worksheet2->write_string($xlsRow,$xlsCol,"Description",$itemFormat{'header'});
		$worksheet2->set_column(++$xlsCol,$xlsCol,30); # col length
		$worksheet2->write_string($xlsRow,$xlsCol,"Species",$itemFormat{'header'});
	}
	else {
		$worksheet2->set_column(++$xlsCol,$xlsCol,30); # col length
		$worksheet2->merge_range(0,$xlsCol,1,$xlsCol,'Source',$itemFormat{'mergeRowHeader'});
		$worksheet2->merge_range(0,++$xlsCol,1,$xlsCol,'MW (Da)',$itemFormat{'mergeRowHeader'});
		$worksheet2->set_column(++$xlsCol,$xlsCol,80); # col length
		$worksheet2->merge_range(0,$xlsCol,1,$xlsCol,'Description',$itemFormat{'mergeRowHeader'});
		$worksheet2->set_column(++$xlsCol,$xlsCol,30); # col length
		$worksheet2->merge_range(0,$xlsCol,1,$xlsCol,'Species',$itemFormat{'mergeRowHeader'});
		$xlsRow--;
		for (my $i=1;$i<=$nbTransition;$i++) {
			$worksheet2->set_column(++$xlsCol,$xlsCol,80); # col length
			$worksheet2->merge_range($xlsRow,$xlsCol,$xlsRow,2+$xlsCol,"Fragment $i",$itemFormat{'mergeRowHeader'});
			$worksheet2->set_column($xlsCol,$xlsCol,15); # col length
			$worksheet2->write_string(++$xlsRow,$xlsCol,"Area",$itemFormat{'header'});
			$worksheet2->write_string($xlsRow,++$xlsCol,"M/Z",$itemFormat{'header'});
			$worksheet2->write_string($xlsRow,++$xlsCol,"Type",$itemFormat{'header'});
			$xlsRow--;
		}
		$xlsRow++;
	}

	my ($refProteins,$refPeptides,$refCharges);
	if ($labelType eq 'FREE' && $call eq 'quanti') {
		($refProteins,$refPeptides,$refCharges)=($refProtInfo,$refProtMatch,$refPepAll);
	}
	else {
		($refProteins,$refPeptides,$refCharges)=($refProtMatch,$reflabeledPeptideSets,$reflabeledPeptideSets);
	}

	foreach my $matchGroup (sort{$refProteins->{$b}[3] <=> $refProteins->{$a}[3] || $a cmp $b} keys %{$refProteins}) { #ALIAS,PROT_DES,MW,PROT_LENGTH,ORGANISM,NUM_PEP,NUM_MATCH,SCORE,CONF_LEVEL,PEP_COVERAGE,PEP_SPECIFICITY
		my ($protID,$alias,$des,$mw,$length,$org)=@{$refProteins->{$matchGroup}};
		#next if ($refnbQuantifPeptides->{$reftrueMGTopProt->{$matchGroup}}==0);

		if ($labelType=~/ITRAQ|TMT/) {
			my ($reportValue)=($labelType=~/ITRAQ/)?'REP_AREA':'REP_INTENSITY';
			my $maxReporterPos=scalar keys %{$refLabelingInfo};
			foreach my $pepID (sort{$refPeptideData->{$matchGroup}{$a}[1]<=>$refPeptideData->{$matchGroup}{$b}[1] || $refPeptideData->{$matchGroup}{$a}[0] cmp $refPeptideData->{$matchGroup}{$b}[0] || $refPeptideData->{$matchGroup}{$a}[2]<=>$refPeptideData->{$matchGroup}{$b}[2] || $refPeptideData->{$matchGroup}{$b}[3]<=>$refPeptideData->{$matchGroup}{$a}[3]} keys %{$refPeptideData->{$matchGroup}}) {
				$xlsCol=0;
				my @pepInfo=@{$refPeptideSets->{$pepID}[0]};# $pepID,$charge,$beg,$pepSeq,$varModStrg
				#my $dataSrc;
				$worksheet2->write_string(++$xlsRow,$xlsCol,$alias,$itemFormat{'text'});
				$worksheet2->write_string($xlsRow,++$xlsCol,$pepInfo[3],$itemFormat{'text'});		# $pepSeq
				$worksheet2->write_string($xlsRow,++$xlsCol,$pepInfo[4],$itemFormat{'text'});		# $varModStrg
				if($refPepMrObs->{$pepID} eq "-"){
					$worksheet2->write_string($xlsRow,++$xlsCol,$refPepMrObs->{$pepID},$itemFormat{'text'});
				}
				else{
					$worksheet2->write_number($xlsRow,++$xlsCol,$refPepMrObs->{$pepID},$itemFormat{'number'});
				}
				$worksheet2->write_number($xlsRow,++$xlsCol,$pepInfo[2],$itemFormat{'number'});		# $beg
				$worksheet2->write_string($xlsRow,++$xlsCol,"$pepInfo[1]+",$itemFormat{'textC'});		# $charge
				my $score=($refPeptideScore->{$pepID})?$refPeptideScore->{$pepID}:'-';
				$worksheet2->write_string($xlsRow,++$xlsCol,$score,$itemFormat{'text'});
				if ($refQuantifData->{$pepID}{$reportValue}) {
					foreach my $repPos (1..$maxReporterPos) {
						if($refQuantifData->{$pepID}{$reportValue}{$repPos}){
							$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifData->{$pepID}{$reportValue}{$repPos}[0],$itemFormat{'number'});
							#$dataSrc=$refQuantifData->{$pepID}{'REP_AREA'}{$repPos}[1] unless $dataSrc;
						}
						else {
							$worksheet2->write_string($xlsRow,++$xlsCol,'-',$itemFormat{'text'});
						}
					}
				}
				else{
					foreach my $repPos (1..$maxReporterPos) {
						$worksheet2->write_string($xlsRow,++$xlsCol,'-',$itemFormat{'text'});
					}
				}
				$worksheet2->write_string($xlsRow,++$xlsCol,$refPepDataSource->{$pepID},$itemFormat{'text'});
				$worksheet2->write_number($xlsRow,++$xlsCol,$mw,$itemFormat{'number'});
				$worksheet2->write_string($xlsRow,++$xlsCol,"$des",$itemFormat{'text'});
				$worksheet2->write_string($xlsRow,++$xlsCol,$org,$itemFormat{'text'});
			}
		}
		else {
			foreach my $seqVarMod (sort{$refPosBeg->{$a}[0][0]<=>$refPosBeg->{$b}[0][0] || lc($a) cmp lc($b) || $a cmp $b} keys %{$refPeptides->{$matchGroup}} ){

				my @refCharges=($labelType eq 'FREE' && $call eq 'quanti')? keys %{$refPepAll->{$seqVarMod}} : keys %{$reflabeledPeptideSets->{$matchGroup}{$seqVarMod}};

				foreach my $charge (sort @refCharges) {

# TODO: $extraParams[0]->{$seqVarMod}{$charge}{$anaID} is now an array ref => update following code to match (PP 08/02/17)
# TODO: Check @extraParams compatibility with all quantif/labeling types (PP 08/02/17)

					$xlsCol=0;
					my ($pepID,$score,$source,@pepInfo,@currentQuantiValues,@currentSourceValues);
					my $paramCode;
					if ($call eq 'quanti' && $labelType eq 'FREE') { # XIC : multi ana
						my (@scores,$currentPepID);
						foreach my $anaID (@{$refAnaOrder}) {
							my $sc='-';
							if ($refPepAll->{$seqVarMod}{$charge}{$anaID}) {
								$currentPepID=$refPepAll->{$seqVarMod}{$charge}{$anaID}[0]; #??? possible multiple occurence not handled != list display
								if ($printMCQ_RT) {
									if (defined($refQuantifData->{$currentPepID}{$quantifParamInfo{'RT_BEGIN'}}[1])) {
										push @currentQuantiValues,(&formatRTinmin($refQuantifData->{$currentPepID}{$quantifParamInfo{'RT_BEGIN'}}[1]),&formatRTinmin($refQuantifData->{$currentPepID}{$quantifParamInfo{'RT_END'}}[1]));
									}
									else {
										push @currentQuantiValues,("-","-");
									}
								}
								if ($xicSoftCode ne 'PKV' && $xicSoftCode ne 'OS') {
									push @currentQuantiValues,$refQuantifData->{$currentPepID}{$quantifParamInfo{'XIC_AREA'}}[1];
								}
								push @currentSourceValues,$refPepDataSource->{$currentPepID} || "-";
								$pepID=$currentPepID unless $pepID;
								$pepID=$currentPepID if $refPepMrObs->{$pepID};
								$sc=($refPeptideScore->{$pepID})?$refPeptideScore->{$pepID} : '-';
							}
							else {
								push @currentQuantiValues,"-";
								push @currentSourceValues,"-";
							}

							push @scores,$sc;
						}
						$score=join('/',@scores);
						$source=join('/',@currentSourceValues);
						$paramCode='XIC_AREA';
					}
					elsif ($call eq 'ana' && $labelType eq 'FREE') {# XIC: single ana
						$pepID=$reflabeledPeptideSets->{$matchGroup}{$seqVarMod}{$charge}[0];
						if ($printMCQ_RT) {
							if (defined($refQuantifData->{$pepID}{$quantifParamInfo{'RT_BEGIN'}}[1])) {
								push @currentQuantiValues,(&formatRTinmin($refQuantifData->{$pepID}{$quantifParamInfo{'RT_BEGIN'}}[1]),&formatRTinmin($refQuantifData->{$pepID}{$quantifParamInfo{'RT_END'}}[1]));
							}
							else {
								push @currentQuantiValues,("-","-");
							}
						}
						#@pepInfo=@{$refPeptideSets->{$pepID}[0]};# $pepID,$charge,$beg,$pepSeq,$varModStrg
						if ($xicSoftCode ne 'PKV' &&  $xicSoftCode ne 'OS') {
							push @currentQuantiValues,$refQuantifData->{$pepID}{$quantifParamInfo{'XIC_AREA'}}[1];
						}
						#else{
						#	push @currentQuantiValues,
						#}
						$score=($refPeptideScore->{$pepID})?$refPeptideScore->{$pepID}:'-';
						$source=$refPepDataSource->{$pepID} || '-';
						$paramCode='XIC_AREA';
					}
					elsif ($labelType eq 'SILAC') {
						my (@scores);
						foreach my $dataSrc (sort keys %{$reflabeledPeptideSets->{$matchGroup}{$seqVarMod}{$charge}}) {
							foreach my $chanNum (@{$refChannelList}) { #1..$maxChanNum
								my $sc='-';
								if ($reflabeledPeptideSets->{$matchGroup}{$seqVarMod}{$charge}{$dataSrc}{$chanNum}) {
									my $currentPepID=(sort{$refPeptideScore->{$b}<=>$refPeptideScore->{$a}} @{$reflabeledPeptideSets->{$matchGroup}{$seqVarMod}{$charge}{$dataSrc}{$chanNum}})[0];
									$sc=($refPeptideScore->{$currentPepID})? $refPeptideScore->{$currentPepID} : '-'; # ghost peptides don't have scores
									if ($refQuantifData->{$currentPepID}{$quantifParamInfo{'ISO_AREA'}}[1]) {
										push @currentQuantiValues,$refQuantifData->{$currentPepID}{$quantifParamInfo{'ISO_AREA'}}[1];
										$pepID=$currentPepID;
									}
								}else{
									push @currentQuantiValues,'-';
								}
								push @scores,$sc;
							}
							#$score=join('/',@scores);
							$score="-";
							@pepInfo=@{$refPepAll->{$pepID}[0]};# $pepID,$charge,$beg,$pepSeq,$varModStrg
							$worksheet2->write_string(++$xlsRow,$xlsCol,$alias,$itemFormat{'text'});
							$worksheet2->write_string($xlsRow,++$xlsCol,$pepInfo[3],$itemFormat{'text'});
							$worksheet2->write_string($xlsRow,++$xlsCol,$pepInfo[4],$itemFormat{'text'});
							$worksheet2->write_number($xlsRow,++$xlsCol,$pepInfo[2],$itemFormat{'number'});
							$worksheet2->write_string($xlsRow,++$xlsCol,"$charge+",$itemFormat{'text'});
							$worksheet2->write_string($xlsRow,++$xlsCol,$score,$itemFormat{'text'});
							foreach my $qdata (@currentQuantiValues) {
								if ($qdata eq '-') {
									$worksheet2->write_string($xlsRow,++$xlsCol,$qdata,$itemFormat{'number'});
								}
								else{
									$worksheet2->write_number($xlsRow,++$xlsCol,$qdata,$itemFormat{'number'});
								}
							}
							$worksheet2->write_string($xlsRow,++$xlsCol,$dataSrc,$itemFormat{'text'});
							$worksheet2->write_number($xlsRow,++$xlsCol,$mw,$itemFormat{'number'});
							$worksheet2->write_string($xlsRow,++$xlsCol,"$des",$itemFormat{'text'});
							$worksheet2->write_string($xlsRow,++$xlsCol,$org,$itemFormat{'text'});
						}
						$paramCode='ISO_AREA';
						next;
					}
					#else{# iTRAQ
					#	$pepID=(sort{$refPeptideScore->{$b}<=>$refPeptideScore->{$a}} @{$reflabeledPeptideSets->{$matchGroup}{$seqVarMod}{$charge}})[0];
					#	push @currentQuantiValues,$refQuantifData->{$pepID}[1];
					#}
					@pepInfo=@{$refPeptideSets->{$pepID}[0]};# $pepID,$charge,$beg,$pepSeq,$varModStrg
					$worksheet2->write_string(++$xlsRow,$xlsCol,$alias,$itemFormat{'text'}); # prot
					$worksheet2->write_string($xlsRow,++$xlsCol,$pepInfo[3],$itemFormat{'text'}); # seq
					$worksheet2->write_string($xlsRow,++$xlsCol,$pepInfo[4],$itemFormat{'text'}); # vmod
					#$worksheet2->write_string($xlsRow,++$xlsCol,$pepInfo[3],$itemFormat{'text'});
					#$worksheet2->write_string($xlsRow,++$xlsCol,$pepInfo[4],$itemFormat{'text'});
					if ($refPepMrObs->{$pepID} eq "-") { # MrObs
						$worksheet2->write_string($xlsRow,++$xlsCol,$refPepMrObs->{$pepID},$itemFormat{'text'});
					}
					else {
						$worksheet2->write_number($xlsRow,++$xlsCol,$refPepMrObs->{$pepID},$itemFormat{'number'});
					}

					$worksheet2->write_number($xlsRow,++$xlsCol,$pepInfo[2],$itemFormat{'number'}); # prot start


					#$worksheet2->write_number($xlsRow,++$xlsCol,$pepInfo[2],$itemFormat{'number'});
					$worksheet2->write_string($xlsRow,++$xlsCol,"$charge+",$itemFormat{'text'});
					$worksheet2->write_string($xlsRow,++$xlsCol,$score,$itemFormat{'text'});
					if ($quantifType !~ /^(SWATH|TDA|DIA)$/) {
						foreach my $qdata (@currentQuantiValues) {
							if ($qdata eq '-') {
								$worksheet2->write_string($xlsRow,++$xlsCol,$qdata,$itemFormat{'text'});
							}
							else{
								$worksheet2->write_number($xlsRow,++$xlsCol,$qdata,$itemFormat{'number'});
							}
						}
					}
					$worksheet2->write_string($xlsRow,++$xlsCol,$source,$itemFormat{'text'});
					$worksheet2->write_number($xlsRow,++$xlsCol,$mw,$itemFormat{'number'});
					$worksheet2->write_string($xlsRow,++$xlsCol,"$des",$itemFormat{'text'});
					$worksheet2->write_string($xlsRow,++$xlsCol,$org,$itemFormat{'text'});
					if ($quantifType=~/^(SWATH|TDA|DIA)$/) {
						my ($i,$fragMWList,$fragTypeList);
						foreach my $fragment (keys %{$refPeptideFragments->{$pepID}}){
							my ($fragMW,$fragCharge,$fragRT,$fragType,$fragArea,$fragIonType,$fragRes)=@{${$refPeptideFragments}{$pepID}{$fragment}};
							#$fragMWList=($fragMWList) ? ($fragMWList.="\\$fragMW") : $fragMW;
							#$fragTypeList=($fragTypeList) ? ($fragTypeList.="\\$fragType") : $fragType;
							$worksheet2->write_number($xlsRow,++$xlsCol,$fragArea,$itemFormat{'number'});
							$worksheet2->write_number($xlsRow,++$xlsCol,$fragMW,$itemFormat{'number'});
							$worksheet2->write_string($xlsRow,++$xlsCol,$fragType,$itemFormat{'text'});
						}
						#$worksheet2->write_string($xlsRow,++$xlsCol,$fragMWList,$itemFormat{'number'});
						#$worksheet2->write_string($xlsRow,++$xlsCol,$fragTypeList,$itemFormat{'text'});
					}
				}
			}
		}
	}
}

sub getAnalysisDataSources {
	my ($anaID,$projectID,$refSourceFiles)=@_;
	return if $refSourceFiles->{$anaID}; # already retrieved
	my ($validStatus,$dataFile)=$dbh->selectrow_array("SELECT VALID_STATUS,DATA_FILE FROM ANALYSIS WHERE ID_ANALYSIS=$anaID");
	my $anaDir=($validStatus==1)? "$promsPath{valid}/ana_$analysisID" : "$promsPath{peptide}/proj_$projectID/ana_$anaID";

	$refSourceFiles->{$anaID}{0}=$dataFile;
	if (-e "$anaDir/mergedFiles.txt") {
		open(MERGE,"$anaDir/mergedFiles.txt") or die $!;
		while (<MERGE>) {
			chomp;
			my ($rank,$fileName,$fileID)=split(/\t/,$_);
			$refSourceFiles->{$anaID}{$rank}=$fileName;
		}
		close MERGE;
	}
}

sub showLog10Plot {
	my ($anaID)=@_;

	#######################
	####>Starting HTML<####
	#######################
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1); # writes PERL warnings as HTML comments <!-- --> must be activated after header
	print qq
|<HTML>
<HEAD>
<TITLE>List of Quantifications</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.TD {font-weight:normal;}
.TH {font-weight:bold;}
.LINK {cursor:pointer;}
.popup {z-index:999;background-color:#FFFFFF;border:solid 3px #999999;padding:5px;position:absolute;display:none;box-shadow:10px 10px 20px #808080;}
</STYLE>
<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/genericPlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/peptidePlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
</SCRIPT>
</HEAD>
<DIV id="waitDiv"><CENTER><BR><BR><BR><BR><BR><FONT class="title3">Fetching data. Please wait...</FONT><BR><IMG src="$promsPath{images}/scrollbarGreen.gif"></CENTER></DIV>
|;

	################################################
	####>Retrieve peptide quantification values<####
	################################################
	my $dbh=&promsConfig::dbConnect;

	my (%pepSequence);
	my $sthPA=$dbh->prepare("SELECT P.ID_PEPTIDE,PEP_SEQ,GROUP_CONCAT(DISTINCT PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),CHARGE FROM PEPTIDE P LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE WHERE P.ID_ANALYSIS=? GROUP BY P.ID_PEPTIDE");
	$sthPA->execute($anaID);
	while (my ($pepID,$pepSeq,$modCode,$charge)=$sthPA->fetchrow_array) {
		#next unless $usedPeptides{$pepID};
		my $modStrg=($modCode)? '+'.&promsMod::decodeVarMod($dbh, $pepSeq, $modCode) : '';
		$pepSequence{$pepID}=$pepSeq.$modStrg.' '.$charge.'+';
	}
	$sthPA->finish;

	my %allQuantifParams;
	my $sthQPar=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,CODE FROM QUANTIFICATION_PARAMETER");
	$sthQPar->execute;
	while (my ($paramID,$paramCode)=$sthQPar->fetchrow_array) {
		$allQuantifParams{$paramID}=$paramCode;
	}
	$sthQPar->finish;

	my ($projectID)=&promsMod::getProjectID($dbh,$anaID,'analysis');

	my $sthPQ=$dbh->prepare("SELECT Q.ID_QUANTIFICATION,NAME,QUANTIF_ANNOT FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE Q.ID_QUANTIFICATION=AQ.ID_QUANTIFICATION AND FOCUS='peptide' AND ID_ANALYSIS=?");
	$sthPQ->execute($anaID);
	my (%pepQData,%quantiNames,@quantiOrder,%channelNames); #,%usedPeptides
	while (my ($quantiID,$name,$quantifAnnot)=$sthPQ->fetchrow_array) {
		$quantiNames{$quantiID}=$name;
		push @quantiOrder,$quantiID;
		my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
		my ($labelType)=($labelStrg)? ($labelStrg=~/LABEL=(.+)/) : ('FREE');
		$labelType=uc($labelType);
		my %pepQuantiFiles;
		if ($labelType eq 'FREE') {
			$pepQuantiFiles{0}='peptide_quantification.txt';
			$channelNames{0}{$quantiID}='';
		}
		elsif ($labelType eq 'SILAC') {
			foreach my $infoStrg (@labelInfo) {
				my ($chanNum,$chanName,$labelStrg)=split(';',$infoStrg); # '1;Light;No label;'   '2;Heavy;Label:13C(6);K'
				last if $chanNum !~ /^\d/;
				$pepQuantiFiles{$chanNum}="peptide_quantification_$chanNum.txt";
				$channelNames{$chanNum}{$quantiID}=$chanName;
			}
		}
		elsif ($labelType=~/ITRAQ|TMT/) {
			foreach my $infoStrg (@labelInfo) {
				my ($repPos,$reporter,$monoMass)=split(';',$infoStrg);
				last if $repPos !~ /^\d/;
				$pepQuantiFiles{$repPos}="peptide_quantification_$repPos.txt";
				$channelNames{$repPos}{$quantiID}=$reporter;
			}
		}
		foreach my $targetPos (keys %pepQuantiFiles) {
			open(QUANT,"$promsPath{quantification}/project_$projectID/quanti_$quantiID/$pepQuantiFiles{$targetPos}");
			while(<QUANT>) {
				next if $.==1;
				chomp;
				my ($paramID,$pepID,$qValue)=split(/\t/,$_);
				next if $allQuantifParams{$paramID} !~ /_(INTENSITY|AREA)/;
				next unless $pepSequence{$pepID}; # restrict to peptide in selected analysis (in case multi-ana quantif)
				$pepQData{$targetPos}{$pepID}{$quantiID}=$qValue;
				#$usedPeptides{$pepID}=1;
			}
			close QUANT;
		}
	}
	$sthPQ->finish;

	$dbh->disconnect;


	print qq
|<SCRIPT LANGUAGE="JavaScript">
var hColors=['#E18B6B','#95B9C7','#7E2217','#9A9A9A','#8AFB17','#FBB917','#F660AB'];
var GP;
window.onload=function() {
	GP=new genericPlot({div:'mainGraphDIV',name:'GP',width:600,height:600,
						axisX:{title:'Log10(dataset 1)'}, // ,zoomable:true
						axisY:{title:'Log10(dataset 2)'}, // ,zoomable:true
						sameScale: true,
						axisClosure:true,
						zoomable:true,
						pointLabel:myPointLabel,
						}); //,pointOnList:listSelected
	//GP.addThreshold({axis:'X',label:'Ratio=1',value:0,color:'#555',pattern:'-'});
	//GP.addThreshold({axis:'Y',label:'Ratio=1',value:0,color:'#555',pattern:'-'});
|;


	###>Quantif sets<###
	my $quantifIdx=0;
	foreach my $quantifIdx1 (0..$#quantiOrder-1) {
		foreach my $quantifIdx2 (($quantifIdx1+1)..$#quantiOrder) {
			foreach my $targetPos (sort{$a<=>$b} keys %channelNames) { # SAME LABELING is expected between all quanti!!!!
				my $channelStrg1=($targetPos)? ' ['.$channelNames{$targetPos}{$quantiOrder[$quantifIdx1]}.']' : '';
				my $channelStrg2=($targetPos)? ' ['.$channelNames{$targetPos}{$quantiOrder[$quantifIdx2]}.']' : '';
				print "\tGP.addDataSet($quantifIdx,{name:'$quantiNames{$quantiOrder[$quantifIdx1]}$channelStrg1 vs $quantiNames{$quantiOrder[$quantifIdx2]}$channelStrg2',axisX:'X',axisY:'Y',color:hColors[$quantifIdx]});\n";
				$quantifIdx++;
			}
		}
	}
	###>Quantif data<###$selectedQuantifications[$quantifIdx1]
	$quantifIdx=0;
	my $log10=log(10);
	foreach my $quantifIdx1 (0..$#quantiOrder-1) {
		foreach my $quantifIdx2 (($quantifIdx1+1)..$#quantiOrder) {
			foreach my $targetPos (sort{$a<=>$b} keys %channelNames) {
				my $count=0;
				foreach my $pepID (keys %{$pepQData{$targetPos}}) {
					next unless $pepQData{$targetPos}{$pepID}{$quantiOrder[$quantifIdx1]} && $pepQData{$targetPos}{$pepID}{$quantiOrder[$quantifIdx2]};
					if ($count==0) {
						print "\tGP.addDataAsString($quantifIdx,'";
					}
					else {print ";";}
					$count++;
					print "$pepSequence{$pepID},,";
					printf "%.4f",log($pepQData{$targetPos}{$pepID}{$quantiOrder[$quantifIdx1]})/$log10;
					print ',';
					printf "%.4f",log($pepQData{$targetPos}{$pepID}{$quantiOrder[$quantifIdx2]})/$log10;
					print ',1';
					if ($count==100) {
						print "');\n";
						$count=0;
					}
				}
				if ($count > 0) {
					print "');\n";
				}
				$quantifIdx++;
			}
		}
	}
	print qq
|	GP.draw();
	document.getElementById('waitDiv').style.display='none';
}
function myPointLabel(dp,type) {
    return dp.label+'\\nX='+Math.round(10**dp.x)+'\\nY='+Math.round(10**dp.y); //dp.dataSet.params.name+'\\n'+
}
</SCRIPT>
<DIV id="mainGraphDIV" style="float:left;background-color:#FFFFFF;border:solid 3px $darkColor;padding:2px;"></DIV>
<DIV id="displayDIV" class="popup"> <!--filter:alpha(opacity=80);opacity:0.8;-->
	<DIV id="infoDIV"></DIV>
</DIV>
</BODY>
</HTML>
|;

}

sub formatRTinmin {
	my ($rtNumber)=@_;
	return $rtNumber unless $rtNumber;
	$rtNumber=sprintf '%.2f', $rtNumber / 60;
	return($rtNumber);
}

sub ajaxXicTrace {
	my ($quantiID)=@_;
	my $pepID=int param("id_peptide");

	unlink "$promsPath{tmp}/${userID}_globalXIC.png" if -e "$promsPath{tmp}/${userID}_globalXIC.png";
	unlink "$promsPath{tmp}/${userID}_globalXIC.png" if -e "$promsPath{tmp}/${userID}_globalXIC.png";

	#-------------------------------------------------------------------------------#
	# Read TraceFile in directory and get the rt trace file of the selected peptide #
	#-------------------------------------------------------------------------------#
	my $dbh=&promsConfig::dbConnect;
	my $projectID=promsMod::getProjectID($dbh,$quantiID,'QUANTIFICATION');
	my ($pepSeq,$mrObs,$charge)=$dbh->selectrow_array("SELECT PEP_SEQ,MR_OBS,CHARGE FROM PEPTIDE WHERE ID_PEPTIDE=$pepID");
	$mrObs=($mrObs)? "Mr(Obs)=$mrObs" : "";
	my ($rtBegin)=$dbh->selectrow_array("SELECT QUANTIF_VALUE FROM PEPTIDE_QUANTIFICATION PQ, QUANTIFICATION_PARAMETER QP WHERE PQ.ID_QUANTIF_PARAMETER=QP.ID_QUANTIF_PARAMETER AND ID_QUANTIFICATION=$quantiID AND ID_PEPTIDE=$pepID AND CODE='RT_BEGIN'");
	my ($rtEnd)=$dbh->selectrow_array("SELECT QUANTIF_VALUE FROM PEPTIDE_QUANTIFICATION PQ, QUANTIFICATION_PARAMETER QP WHERE PQ.ID_QUANTIF_PARAMETER=QP.ID_QUANTIF_PARAMETER AND ID_QUANTIFICATION=$quantiID AND ID_PEPTIDE=$pepID AND CODE='RT_END'");
	$dbh->disconnect;
	open(PEPTRACE,"$promsPath{quantification}/project_$projectID/quanti_$quantiID/pepID_to_traces.txt");
	my $traceFile;
	while (my $line=<PEPTRACE>) {
		chomp($line);
		my ($id,$rtFile)=split(/\t/,$line);
		if ($id eq $pepID) {
			$traceFile="$promsPath{quantification}/project_$projectID/quanti_$quantiID/all_xics_traces/$rtFile";
			last;
		}
	}
	close(PEPTRACE);

	#----------#
	# R script #
	#----------#
	my $userID = $ENV{'REMOTE_USER'};
	my $tmpDir = strftime("%Y%m%d%H%M%S",localtime).$userID;
	mkdir "$promsPath{tmp}/$tmpDir";
	my $rScript = "$promsPath{tmp}/$tmpDir/plot.R";

	open Rscript, ">$rScript";
	print Rscript qq
|traceFile <- read.csv("$traceFile",sep="\t",skip=1)
#par(mfrow=c(ncol(traceFile)-3,1))
#for(i in 2:(ncol(traceFile)-2))
#{
#	plot(traceFile[,1],traceFile[,i],type='l',xlab="Retention time",ylab=colnames(traceFile)[i])
#	points(traceFile[,1],traceFile[,ncol(traceFile)],pch=23,col=554)
#}
###> All chromatogram intensities
png(filename="$promsPath{tmp}/${userID}_globalXIC.png", width=1000, height=500, units = "px")
plot(traceFile[,1]/60,traceFile[,2],type="l",xlab="Retention time",ylab="Intensity")
points(traceFile[,1]/60,traceFile[,7],pch=23,col=554)
dev.off()

###> XIC intensities
png(filename="$promsPath{tmp}/${userID}_localXIC.png", width=1000, height=500, units = "px")
rtindBeg <- max(which(traceFile[,1]<$rtBegin))
rtindEnd <- max(which(traceFile[,1]<$rtEnd))
plot(traceFile[rtindBeg:rtindEnd,1]/60,traceFile[rtindBeg:rtindEnd,2],type="l",xlab="Retention time",ylab="Intensity")
dev.off()

|;
	close Rscript;
	system "cd $promsPath{tmp}/$tmpDir; $promsPath{R}/R CMD BATCH --no-save --no-restore $rScript";

	unlink "$promsPath{tmp}/$tmpDir/plot.R";
	unlink "$promsPath{tmp}/$tmpDir/plot.Rout";
	rmtree "$promsPath{tmp}/$tmpDir";

	move("$promsPath{tmp}/$userID\_globalXIC.png", "$promsPath{quantification}/project_$projectID/quanti_$quantiID");
	move("$promsPath{tmp}/$userID\_localXIC.png", "$promsPath{quantification}/project_$projectID/quanti_$quantiID");

	#------#
	# HTML #
	#------#
	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'utf-8',-cache_control=>"no-cache");
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Display Extracted Ion Chromatogram (XIC)</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.popup {z-index:999;background-color:#FFFFFF;border:solid 3px #999999;padding:5px;position:absolute;display:none;box-shadow:10px 10px 20px #808080;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript"></SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<CENTER>
<FONT class="title">Local and global trace files for $pepSeq [z=$charge] $mrObs</FONT>
&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'"/>
<BR><BR>
<IMG src="$promsPath{html}/data/quantification_data/project_$projectID/quanti_$quantiID/${userID}_localXIC.png" id="localPicture"/>
<IMG src="$promsPath{html}/data/quantification_data/project_$projectID/quanti_$quantiID/${userID}_globalXIC.png" id="globalPicture"/>
</SCRIPT>
</CENTER>
</BODY>
</HTML>
|;


}

sub ajaxShowPepTDAList {
	my $nbTransition = 1;
	my $protID = param('id_protein');
	my $analysisID = param('id_analysis');
	
	print header(-type=>'text/plain',-charset=>'utf-8',-cache_control=>"no-cache");
	warningsToBrowser(1);
	
	if ($analysisID) { # Fetch all analysis and corresponding peptides from the quantification
		&fetchQuantiAnalysis($selQuantifID, $analysisID);
	} else { # Internal Quantifications
		&fetchQuantiAnalysis($selQuantifID);
	}
	
	### Fetch all proteins for current analysis
	foreach my $anaID (@anaOrder) {
		&fetchAnaProteins($anaID, 1, $protID);
		$nbTransition = &fetchTDAPeptides($anaID, $protID);
	}
	
	# Displaying protein data
	print qq |
			<br/>
			<TABLE border=0 cellspacing=0 cellpadding=2>
	|;
	
	if (%peptideFragments || !$analysisID) { # Is MS2 or Multi Analysis
		print "<TR><TH colspan=6>&nbsp;&nbsp;</TH>\n";
		
		foreach my $anaID (@anaOrder) {
			### Compute fragment columns length based on their amount
			my $colSpan = $nbTransition+1;
			print "<TH bgcolor=\"$darkColor\" class=\"rbBorder\" colspan=$colSpan>&nbsp;$listAna{$anaID}&nbsp;</TH>\n";
		}
		
		print "</TR>\n";
	}
	
	# Displaying default table header 
	print qq |
		<TR>
			<TH>&nbsp;&nbsp;</TH><TH bgcolor="$darkColor" class="rbBorder">#</TH>
			<TH bgcolor="$darkColor" class="rbBorder">Peptides&nbsp;&nbsp;</TH>
			<TH bgcolor="$darkColor" class="rbBorder">&nbsp;Mr(Obs)&nbsp;</TH>
			<TH bgcolor="$darkColor" class="rbBorder">&nbsp;Start&nbsp;</TH>
			<TH bgcolor="$darkColor" class="rbBorder">&nbsp;Charge&nbsp;</TH>
	|;
	
	## Displaying quantification type-specific header(s)
	foreach my $anaID (@anaOrder) {
		if(%peptideFragments) {
			print "<TH bgcolor=\"$darkColor\" class=\"rbBorder\">&nbsp;RT<SUB>mean</SUB></TH>\n";
			for (my $i=1; $i<$nbTransition+1; $i++) {
				print "<TH bgcolor=\"$darkColor\" class=\"rbBorder\">&nbsp;Fragment&nbsp;$i&nbsp;</TH>";
			}
		} else {
			print "<TH bgcolor=\"$darkColor\" class=\"rbBorder\">&nbsp;RT&nbsp;</TH>\n";
			print "<TH bgcolor=\"$darkColor\" class=\"rbBorder\">&nbsp;XIC&nbsp;</TH>\n";
		}	
	}
	print qq |
			<TH width=50%></TH>
		</TR>
	|;
	
	# Displaying quantification values
	&displayQuantificationValues($protID, $nbTransition);
				
	print qq |
			<TR><TD colspan=$nbTransition>&nbsp;</TD></TR>
		</TABLE>
	|;
}

sub ajaxShowPepTDAGraph {
	## Store results data as JS string 
	my ($jsAnaColumnsStr, $jsAddRowStr, $jsPepGroupsStr, $jsPeptideIdStr) = "";
	my $isMS2 = 0; # Specify if it should look for peptide or fragments data
	
	## Connecting to the database
	my $dbh=&promsConfig::dbConnect;
	
	my $protID=param('id_protein');
	my ($protAlias,$protLength)=$dbh->selectrow_array("SELECT ALIAS,PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=$protID");
	
	my $projectID=&promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION');
	
	
	### Fetching analysis info
	my $sthRun=$dbh->prepare("SELECT AQ.ID_ANALYSIS, A.NAME, Q.NAME FROM ANA_QUANTIFICATION AQ, ANALYSIS A, QUANTIFICATION Q WHERE AQ.ID_QUANTIFICATION=$selQuantifID AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND A.ID_ANALYSIS=AQ.ID_ANALYSIS");
	$sthRun->execute;
	my %runInfo;
	while(my ($anaID,$anaName,$quantiName)=$sthRun->fetchrow_array){
		@{$runInfo{$anaID}}=($anaName,$quantiName);
	}
	$sthRun->finish;
	my $anaString=join(',',keys %runInfo);

	
	## Setting columns: runs and States
	my $anaCount = 0; # index
	$jsAnaColumnsStr .= "HMF_$protID.setColumns([";
	foreach my $anaID (sort {$runInfo{$a}[0] cmp $runInfo{$b}[0]} keys %runInfo) {
		$anaCount++;
		$jsAnaColumnsStr .= "['$runInfo{$anaID}[0]',$anaID,'$runInfo{$anaID}[1]']";
		$jsAnaColumnsStr .= ',' if $anaCount < scalar keys %runInfo;
		
		# Check for fragments data
		my $swathFile = $promsPath{"quantification"}."/project_$projectID/quanti_$selQuantifID/swath_ana_$anaID.txt";
		$isMS2 = 1 if(-e $swathFile && !$isMS2);
	}
	$jsAnaColumnsStr .= "]);\n";
	
	## Fetching peptides info
	my (%peptideInfo, %peptideBeg, %varModText, %peptideID, %peptideDataGost);
	my $sthPI=$dbh->prepare("SELECT PEP_SEQ,GROUP_CONCAT(PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),CHARGE,ABS(PEP_BEG),P.ID_ANALYSIS,P.VALID_STATUS,P.ID_PEPTIDE,QUERY_NUM,PEP_RANK,SCORE
							FROM PEPTIDE P
							LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
							INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE AND P.ID_ANALYSIS IN ($anaString)
							WHERE ID_PROTEIN=? GROUP BY P.ID_PEPTIDE");
	$sthPI->execute($protID);
	while (my ($pepSeq,$varModCode,$charge,$beg,$anaID,$validStatus,$pepID,$queryNum,$pepRank,$score)=$sthPI->fetchrow_array) {
		my $vModStrg;
		if ($varModCode) {
			$vModStrg='&'.$varModCode;
		} else {
			$varModCode = $vModStrg = '';
		}
		
		$peptideInfo{$pepSeq}{$varModCode}{$charge}{$anaID} = $pepID;
		$peptideBeg{$pepSeq} = $beg;
		$peptideID{$pepID} = {
			"seq"    	 => $pepSeq,
			"charge" 	 => $charge,
			"anaID"	 	 => $anaID,
			"varMod" 	 => $varModCode,
			"queryNum"	 => $queryNum,
			"pepRank"	 => $pepRank,
			"score"		 => $score,
			"ghost"		 => ($validStatus == 0) ? 1 : 0,
		};
		$varModText{$varModCode} = &promsMod::decodeVarMod($dbh, $pepSeq, $varModCode) if $varModCode;
	}
	$sthPI->finish;
	$dbh->disconnect;

	## SHOULD USE PEPTIDE FRAGMENTS DATA
	if($isMS2) {
		## Fetching fragments info
		my %peptideData;
		foreach my $anaID (keys %runInfo){
			my $swathFile = $promsPath{"quantification"}."/project_$projectID/quanti_$selQuantifID/swath_ana_$anaID.txt";
			open(IN, $swathFile);
			while(<IN>){
				next if $. == 1;
				my ($pepID,$fragMZ,$fragCharge,$fragType,$fragRes,$area,$fragRT) = split('!', $_);
				$fragRT = ($fragRT) ? $fragRT : '';
				chomp $fragRT;
				chomp $area;
				if($peptideID{$pepID}) {
					# pepKey = seq&varMod_charge if varMod, otherwise seq_charge
					my $pepKey= ($peptideID{$pepID}{"varMod"})? $peptideID{$pepID}{"seq"}.'&'.$peptideID{$pepID}{"varMod"}.'_'.$peptideID{$pepID}{"charge"} : $peptideID{$pepID}{"seq"}.'_'.$peptideID{$pepID}{"charge"};
					if($peptideID{$pepID}{"ghost"}) {
						@{$peptideDataGost{$pepKey}{$fragType.$fragRes.'_'.$fragCharge}{$anaID}}=($area,$fragRT,$fragMZ);
					} else {
						@{$peptideData{$pepKey}{$fragType.$fragRes.'_'.$fragCharge}{$anaID}}=($area,$fragRT,$fragMZ);
					}
				}
			}
			close IN;
		}
		
		### Adding rows: peptide fragments
		my (%peptideGroups, %fragmentRT, %fragmentGostRT);
		my $pepCount = 0; # already used for List view
		my $fragCount = 0;
		$jsPeptideIdStr ="var peptideId_$protID={\n";
		
		foreach my $pepSeq (sort{$peptideBeg{$a}<=>$peptideBeg{$b} || length($a)<=>length($b)} keys %peptideInfo) {
			foreach my $varMod (sort keys %{$peptideInfo{$pepSeq}}) {
				foreach my $charge (sort{$a<=>$b} keys %{$peptideInfo{$pepSeq}{$varMod}}) {
					my $lastFragNb = $fragCount;
					my $nbFrag;
					my $pepKey=($varMod)? $pepSeq.'&'.$varMod.'_'.$charge : $pepSeq.'_'.$charge; # pepKey = seq&varMod.charge if varMod, otherwise seq.charge
					$pepCount++;
					
					@{$peptideGroups{$pepCount}{'PEP'}} = ($pepSeq, $varMod, $charge);
					$peptideGroups{$pepCount}{'FRAG'}[0] = $fragCount; # index unless $peptideGroups{$pepCount}{'FRAG'};
					
					$jsPeptideIdStr .= ",\n" if $pepCount > 1;
					$jsPeptideIdStr .= "\t'$pepCount':{";

					# Add "regular" peptides
					foreach my $fragCharge (sort{&promsMod::sortSmart($a,$b)} keys %{$peptideData{$pepKey}}) {
						my ($frag, $ch) = split('_', $fragCharge);
						$peptideGroups{$pepCount}{'FRAG'}[1] = $fragCount;
						
						$jsAddRowStr .= "HMF_$protID.addRow(['$frag','$pepCount:$fragCount','Charge: $ch+'],[";
		
						# Lopping through columns
						$anaCount = 0;
						foreach my $anaID (sort {$runInfo{$a}[0] cmp $runInfo{$b}[0]} keys %runInfo) {
							push @{$fragmentRT{$pepKey}{$anaID}},$peptideData{$pepKey}{$fragCharge}{$anaID}[1];
						}
						$nbFrag++; $fragCount++;
						
						foreach my $anaID (sort {$runInfo{$a}[0] cmp $runInfo{$b}[0]} keys %runInfo) {
							my $pepID = $peptideInfo{$pepSeq}{$varMod}{$charge}{$anaID};
							my $area = ($peptideData{$pepKey}{$fragCharge}{$anaID}[0]) ? ($peptideData{$pepKey}{$fragCharge}{$anaID}[0] eq "None") ? '': $peptideData{$pepKey}{$fragCharge}{$anaID}[0] : '';
							$anaCount++;
							
							$jsAddRowStr .= $area;
							$jsAddRowStr .= ',' if $anaCount < scalar keys %runInfo;
							
							if($nbFrag == scalar keys %{$peptideData{$pepKey}}) {
								if($peptideInfo{$pepSeq}{$varMod}{$charge}{$anaID}) {
									$jsPeptideIdStr .= "'$anaID':['pep_$pepID\_".$peptideID{$pepID}{"queryNum"}."_".$peptideID{$pepID}{"pepRank"}."',".$peptideID{$pepID}{"score"}.",[".(join(',', @{$fragmentRT{$pepKey}{$anaID}}))."],$lastFragNb]";
								} else {
									$jsPeptideIdStr .= "'$anaID':null";
								}
								$jsPeptideIdStr.="," if $anaCount < scalar keys %runInfo;
							}
						}
						$jsAddRowStr .= "]);\n";
					}
					
					# Add ghost peptides
					foreach my $fragCharge (sort{&promsMod::sortSmart($a,$b)} keys %{$peptideDataGost{$pepKey}}) {
						my ($frag,$ch) = split('_',$fragCharge);
						$peptideGroups{$pepCount}{'FRAG'}[1] = $fragCount;
						
						$jsAddRowStr .= "HMF_$protID.addRow(['$frag','$pepCount:$fragCount','Charge: $ch+'],[";
	
						# Lopping through columns
						$anaCount = 0;
						foreach my $anaID (sort {$runInfo{$a}[0] cmp $runInfo{$b}[0]} keys %runInfo) {
							push @{$fragmentGostRT{$pepKey}{$anaID}},$peptideDataGost{$pepKey}{$fragCharge}{$anaID}[1];
						}
						$nbFrag++; $fragCount++;
						
						foreach my $anaID (sort {$runInfo{$a}[0] cmp $runInfo{$b}[0]} keys %runInfo) {
							my $pepID = $peptideInfo{$pepSeq}{$varMod}{$charge}{$anaID};
							my $area = ($peptideDataGost{$pepKey}{$fragCharge}{$anaID}[0] && $peptideDataGost{$pepKey}{$fragCharge}{$anaID}[0] eq "None") ? '': ($peptideDataGost{$pepKey}{$fragCharge}{$anaID}[0])? $peptideDataGost{$pepKey}{$fragCharge}{$anaID}[0] : '';
							$anaCount++;
							
							$jsAddRowStr .= $area;
							$jsAddRowStr .= ',' if $anaCount < scalar keys %runInfo;
							
							if($nbFrag == scalar keys %{$peptideDataGost{$pepKey}}) {
								if($peptideInfo{$pepSeq}{$varMod}{$charge}{$anaID}) {
									$jsPeptideIdStr .= "'$anaID':['pep_$pepID\_".$peptideID{$pepID}{"queryNum"}."_".$peptideID{$pepID}{"pepRank"}."',".$peptideID{$pepID}{"score"}.",[".(join(',', @{$fragmentGostRT{$pepKey}{$anaID}}))."],$lastFragNb]";
								} else{
									$jsPeptideIdStr .= "'$anaID':null";
								}
								$jsPeptideIdStr .= "," if $anaCount < scalar keys %runInfo;
							}
						}
						$jsAddRowStr .= "]);\n";
					}
					$jsPeptideIdStr .= "}";
				}
			}
		}
		$jsPeptideIdStr .= "\n};\n";
		
		### Defining groups
		$jsPepGroupsStr .= "HMF_$protID.defineGroups('row',[";
		foreach my $pepCount (sort{$a<=>$b} keys %peptideGroups) {
			$jsPepGroupsStr .= ',' if $pepCount > 1;
			my ($pepSeq, $varMod, $charge) = @{$peptideGroups{$pepCount}{'PEP'}};
			my $pepSeqVarMod = ($varMod) ? $pepSeq.'+'.$varModText{$varMod} : $pepSeq;
			my $pepEnd = $peptideBeg{$pepSeq} + length($pepSeq)-1;
			$jsPepGroupsStr .= "['$pepSeqVarMod',$pepCount,'Charge: $charge+\\nPosition: $peptideBeg{$pepSeq}-$pepEnd',$peptideGroups{$pepCount}{FRAG}[0],$peptideGroups{$pepCount}{FRAG}[1]]";
		}
		$jsPepGroupsStr .= "]);";
	
	## SHOULD USE PEPTIDE DATA
	} else {
		my $peptideFile = $promsPath{"quantification"}."/project_$projectID/quanti_$selQuantifID/peptide_quantification.txt";
		exit unless(-e $peptideFile);
		
		## Fetching peptide quantification info
		open (QUANTI, $peptideFile);
		while (<QUANTI>) {
			next if $.==1;
			chomp;
			my ($paramID,$pepID,$pepQuanti)=split(/\t/,$_);
			next if(!$peptideID{$pepID});
			$peptideID{$pepID}{"quanti"} = $pepQuanti;
		}
		close QUANTI;
	
		## Adding rows: peptide fragments
		my $pepCount = 1; # already used for List view
		$jsPeptideIdStr = "var peptideId_$protID={\n";
		foreach my $pepSeq (sort{$peptideBeg{$a}<=>$peptideBeg{$b} || length($a)<=>length($b)} keys %peptideInfo) {
			foreach my $varMod (sort keys %{$peptideInfo{$pepSeq}}) {
				foreach my $charge (sort{$a<=>$b} keys %{$peptideInfo{$pepSeq}{$varMod}}) {
					$jsPeptideIdStr.= ",\n" if $pepCount > 1;
					$jsPeptideIdStr.= "\t'$pepCount':{";
					
					# One peptide = one row
					my $pepSeqVarMod=($varMod)? $pepSeq.'+'.$varModText{$varMod} : $pepSeq;
					$jsAddRowStr .= "HMF_$protID.addRow(['$pepSeqVarMod','$pepCount:0','Charge: $charge+\\nPosition: $peptideBeg{$pepSeq}-".($peptideBeg{$pepSeq}+length($pepSeq)-1)."'],[";
					
					my $anaCount=0; # index
					foreach my $anaID (sort {$runInfo{$a}[0] cmp $runInfo{$b}[0]} keys %runInfo) {
						$anaCount++;
						my $pepID = $peptideInfo{$pepSeq}{$varMod}{$charge}{$anaID};
						$jsAddRowStr .= $peptideID{$pepID}{"quanti"};
						$jsAddRowStr .= ',' if $anaCount < scalar keys %runInfo;
						
						$jsPeptideIdStr .= "'$anaID':['pep_$pepID\_".$peptideID{$pepID}{"queryNum"}."_".$peptideID{$pepID}{"pepRank"}."',";
						$jsPeptideIdStr .= $peptideID{$pepID}{"score"}.",";
						$jsPeptideIdStr .= "[],0]";
						$jsPeptideIdStr .= "," if $anaCount < scalar keys %runInfo;
					}
					$jsPeptideIdStr .= "}";
					$jsAddRowStr .= "]);\n";
					$pepCount++;
				}
			}
		}
		$jsPeptideIdStr .= "\n};\n";
	}
	
	if($jsAnaColumnsStr && $jsAddRowStr) {
		printHeatMap($protID, $isMS2, $jsAnaColumnsStr, $jsAddRowStr, $jsPepGroupsStr, $jsPeptideIdStr);
	}
}

sub printHeatMap {
	my ($protID, $isMS2, $jsAnaColumnsStr, $jsAddRowStr, $jsPepGroupsStr, $jsPeptideIdStr) = @_;
	my $rowEntitiesTxt = ($isMS2) ? "Fragments" : "Peptides";
	my $flagTextTxt = ($isMS2) ? "flagText:'excluded fragments'," : "";
	
	print header(-type=>'text/plain',-charset=>'utf-8',-cache_control=>"no-cache");
	
	# Draw HeatMap
	print qq |
var HMF_$protID=new heatMap({
	div:'graphPepTDA$protID',
	cellHeight:15,
	moveLabel:true,
	//editable:true,
	editableItems:{scope:true,type:true,color:true},
	entities:{row:'$rowEntitiesTxt',column:'Run'},
	$flagTextTxt
	cellOnMouseOver: function(value,hCell) {
		var valueStrg;
		if (value==null) {valueStrg='* no value *';}
		else {valueStrg='Intensity: '+(Math.round(10*value)/10);}
		if (hCell) {
			var rowIdData=rowList_$protID\[hCell.rowIdx].id.split(':'), colId=columnList_$protID\[hCell.columnIdx].id;
			if (peptideId_$protID\[rowIdData[0]][colId]) {
				//alert (peptideId_$protID\[rowIdData[0]][colId]);
				//alert(rowIdData);
				
				if($isMS2) {
					valueStrg+='\\nPeptide score: '+(Math.round(100*peptideId_$protID\[rowIdData[0]][colId][1])/100);
					valueStrg+='\\nFragment RT: '+(peptideId_$protID\[rowIdData[0]][colId][2][hCell.rowIdx-peptideId_$protID\[rowIdData[0]][colId][3]]);
				}
			}
			else {valueStrg+='\\n* no peptide identified *';}
		}
		return valueStrg;
	},
	singleCellSelection:true,
	cellOnClick:function(rowId,colId,status) {
		if (!status \|\| !$isMS2) return;
		var rowIdData=rowId.split(':');
		if (peptideId_$protID\[rowIdData[0]][colId]) {drawSpectrum(null,peptideId_$protID\[rowIdData[0]][colId][0]);}
		else {alert('Peptide was not identified in this run');}
	},
	/*rowGroupOnClick:function(groupId){
		var anaPepId='';
		for(var i=0; i<columnList_$protID.length;i++){
			if(peptideId_$protID\[groupId][columnList_$protID\[i].id]){
				var pepIdData=peptideId_$protID\[groupId][columnList_$protID\[i].id][0].split('_');
				anaPepId+=(columnList_$protID\[i].id)+'_'+(pepIdData[1])+"@";
			}
		}
		if(anaPepId){ajaxShowFragTDA(anaPepId);}
	},*/
	normalization:{scope:'row',reference:'user',limitValues:{ref:0,min:0}},
	exportAsImage:['Export as image','Peptide_fragments_heat_map','$promsPath{cgi}/exportSVG.cgi']
});

$jsAnaColumnsStr
$jsAddRowStr
$jsPepGroupsStr
HMF_$protID.draw();

var rowList_$protID=HMF_$protID.getRowList();
var columnList_$protID=HMF_$protID.getColumnList();
$jsPeptideIdStr
|;
}


####>Revision history<####
# 1.9.10 Displays isobaric data correction info if performed by myProMS (PP 24/01/19)
# 1.9.9 Added RT for DIA & TDA (MS1/MS2) (VS 27/11/18)
# 1.9.8 Optimized data loading (VS 26/11/18)
# 1.9.7 Added TDA MS1 Heat Map (VS 21/11/18)
# 1.9.6 Added TDA (Skyline) parameters handling (VS 20/11/18)
# 1.9.5 Modified swath files paths (VS 08/11/2018)
# 1.9.4 Minor modifications (VS 18/10/18)
# 1.9.3 Minor change for TMT display (GA 09/10/18)
# 1.9.2 Minor modifications (GA 26/09/18)
# 1.9.1 Removed usage of PEPTIDE_QUANTIFICATION in &showLog10Plot (PP 27/06/18)
# 1.9.0 Reads peptide quantification data from file & rebuilds match groups after deletion of a peptide quantification (PP 06/06/18)
# 1.8.7 Minor bug fix in display for MQ vs DIA (PP 30/05/18)
# 1.8.6 Bug correction on data export (MLP 24/04/18)
# 1.8.5 Better parsing of elution time string & compatible with MaxQuant (18/04/18)<BR>TODO: check export for all quantif types.
# 1.8.4 Minor modif on heat map display for OpenSwath (MLP 18/12/17)
# 1.8.3 Minor modification of XIC filtering display (GA 13/12/17)
# 1.8.2 Minor correction (MLP 06/12/17)
# 1.8.1 Minor bug corrections (GA 03/11/17)
# 1.8.0 Adapt script for DIA (OpenSwath) (MLP 01/09/17)
# 1.7.9 Adapt script for TDA (PRM, SRM, MRM) (MLP 27/07/17)
# 1.7.8 Adapt script for TMT (GA 03/04/17)
# 1.7.7 Bug fix to display XIC for lower-scoring peptide in quanti call mode<BR> Fix export for all quantif/label/call types (PP/GA 08/03/17)
# 1.7.6 Minor modif for $beg and $sthPI -> for multiple peptide match on a protein (MLP 08/02/17)
# 1.7.5 Compatible with TMT labeling & minor change in peptide MCQSET removal following peptide quantification deletion (PP 18/01/17)
# 1.7.4 Fix missed 2nd useless retention time formating for non-PeakView data (PP 19/12/16)
# 1.7.3 Minor modif for missing fragments on peakview data display (MLP 01/12/2016)
# 1.7.2 Bug fix in useless retention time formating for non-PeakView data (PP 15/11/16)
# 1.7.1 Minor update to sort %maxProtMatch by protein's identifier in data display (MLP 28/07/2016)
# 1.7.0 Updated for SWATH (PeakView quantification) (MLP 09/06/2016)
# 1.6.1 Minor modification to show all MS/MS XICs of a same instance peptide ie seqVmodCharge (GA 29/04/16)
# 1.6.0 Speed optimization & displayed value rounding for scores & XICs (PP 18/02/16)
# 1.5.6 Minor changes for handling undefined $pepData (PP 09/02/16)
# 1.5.5 Updated for new peptide data source management. TODO: Check export of Label-free PD & display APEX_RT (PP 11/08/15)
# 1.5.4 Bug fix in iTRAQ data display & compatibility with PD label-free (PP 30/04/15)
# 1.5.3 Minor bug correction in Mr(Obs)<BR>Add RTbeg & RTend values for XIC values<BR>Add R code for XIC trace plots (GA 14/10/14)
# 1.5.2 Minor modification of MassChroQ parameter extraction summary (GA 31/03/14)
# 1.5.1 Add basic log10 option to compare XIC extractions for dev purposes<BR>TODO: change color definition (GA 26/03/14)
# 1.5.0 Fix export for labeled quantif + data error message for unlabeled (GA 14/03/14)
# 1.4.9 Uses rmtree instead of remove_tree (PP 10/03/14)
# 1.4.8 Change syntax $x=$y // $z to old ()? one (PP 10/03/14)Uses rmtree instead of remove_tree (PP 10/03/14)
# 1.4.7 Minor modif for isotope MCQ extraction (GA 29/01/14)
# 1.4.6 Minor update is source file display & system command removal (PP 08/11/13)
# 1.4.5 iTRAQ -> ITRAQ & minor change in total signal display (PP 07/10/13)
# 1.4.4 Removing obsolete javascript code (FY 04/10/13)
# 1.4.3 Test for defined channel quantif value (PP 30/09/13)
# 1.4.2 Minor update by commenting 'next' line in SILAC to print peptides with not all channels available (GA 23/09/13)
# 1.4.1 PEPTIDE_SET no longer used (PP 16/09/13)
# 1.4.0 Generated from listAnaQuantifications.cgi 1.3.3 (PP 02/09/13)
# 1.3.3 Minor modification in XLS formating for XIC-export (GA 20/08/13)
# 1.3.2 Add export option for XIC-quantification results (GA 14/08/13)
# 1.3.1 Another improvement of new design display (PP ../07/13)
# 1.3.0 Improved display of new design options (PP 11/07/13)
# 1.2.9 Minor change (PP 08/07/13)
# 1.2.8 Fix missing parameters for ms2 alignment during MCQ extraction (GA 05/07/13)
# 1.2.7 extended design management (PP 04/07/13)
# 1.2.6 Better label detection & fix bug ghost peptide modification (PP 04/07/13)
# 1.2.5 Remove VAR_MOD from script (GA 28/05/13)
# 1.2.4 Modifs for unified naming of files between label and unlabel quantifs (PP 16/04/13)
# 1.2.3 Restrict GO decoration selection to qualitative GO analyses (FY 04/04/13)
# 1.2.2 ABS(PEP_BEG) on Label-quantif because of negative values for ghost-peptides (PP 03/04/13)
# 1.2.1 Minor modification in &deleteQuantification -> delete virtual proteins (GA 08/03/13)
# 1.2.0 Minor modification -> ABS(PEP_BEG) because of negative values for ghost-peptides (GA 08/02/13)
# 1.1.9 Ratio selection & compatibility with new peptide filters (PP 15/01/13)
# 1.1.8 Minor display bugs correction (PP 09/01/13)
# 1.1.7 Minor adding column Mr(obs) for XIC printing (GA 07/12/12)
# 1.1.6 Correction of bug in ajaxListSelectedProteins that was not working properly + print XIC (GA 23/11/12)
# 1.1.5 Aesthetic modifications for 'Peptides used' column (GA 16/11/12)
# 1.1.4 Add quantification parameters for XIC extractions + add deletion of ghost peptides associated to a quantification (GA 22/10/12)
# 1.1.3 Compatible with new data path .../results/project_ID/GoID (PP 13/09/12)
# 1.1.2 Minor bug fixes & improvements (PP 30/07/12)
# 1.1.1 Added custom protein list filtering on display (PP 18/07/12)
# 1.1.0 Find all GO analyses linked to Experiment & true number of filtered proteins (PP 08/06/12)
# 1.0.9 Added save proteins to List. Still not fully stable (PP 07/06/12)
# 1.0.8 More code factorisation & minor display changes: Not fully stable (PP 23/04/12)
# 1.0.7 Bug correction & code factorization (PP 19/04/12)
# 1.0.6 Merge 1.0.5GA & 1.0.4c (PP 06/04/12)
# 1.0.4c PEPTIDE.QUANTI_DATA change to DATA (PP 05/04/12)
# 1.0.4b minor change in peptide selection display (PP 27/03/12)
# 1.0.4 added graphical display of peptide for labeled quanti (PP 27/02/12)
# 1.0.3 Added selected sources (merged MSF) & graphical displays (PP 02/02/12)
# 1.0.2 Protein-level internal analysis quantification (PP 19/09/11)
# 1.0.1 No longer uses CHANNEL as quantif parameter (PP 29/07/11)
# 1.0.0 Created (PP 25/07/11)
