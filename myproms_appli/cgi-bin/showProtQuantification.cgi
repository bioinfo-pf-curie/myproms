#!/usr/local/bin/perl -w

################################################################################
# showProtQuantification.cgi     2.4.5                                         #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Displays protein quantification data                                         #
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
$|=1;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use promsQuantif;
use phosphoRS;
use strict;
use Spreadsheet::WriteExcel;
use utf8; # Tells Perl that characters are UTF8. Necessary for Excel export to work with UTF-8 characters ...
use Encode qw(encode_utf8); # ... Encode needed if hard UTF8 chars in code!!!
use File::Path qw(rmtree); # remove_tree

#print header(-charset=>'utf-8'); warningsToBrowser(1); #DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my $projectAccess; # must be global
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my %labelingName=('FREE'=>'Label-free','SILAC'=>'SILAC','ITRAQ'=>'iTRAQ','TMT'=>'TMT');
#my $MAX_INF_RATIO_DB=1000; # true infinite value in database
#my $MIN_INF_RATIO_DB=0.001;
# Limit ratios for +/-infinity switch in volcano plot
my $MAX_INF_RATIO=64;
my $MIN_INF_RATIO=1/$MAX_INF_RATIO;
my $log2=log(2);
my %normalizationNames=&promsQuantif::getQuantifNormalizationName;
my %pValueCorrections=('none'=>'None','fdr'=>'Benjamini-Hochberg (FDR)','bonferroni'=>'Bonferroni'); # for algo v3

####################
####>Parameters<####
####################
my $call=param('CALL') || 'quanti'; # ana or quanti
my $action=param('ACT') || 'select';
my $analysisID=param('id_ana') || 0; # 0=> call=quanti
my $selQuantifID=param('id_quantif') || 0; # 0=> call=ana
my $view=param('view') || 'graph'; # 'list';
my $foldChangeType=param('foldChgType') || 'abs'; $foldChangeType='abs' if $view eq 'graph';
my $dispFoldChange=param('foldChange') || 2; $dispFoldChange=1 if ($dispFoldChange=~/[^\d\.]/ || $dispFoldChange < 1);
my $foldChangeRaw=param('foldChgRaw') || 0; # MaxQuant ratio only (0->NORM or 1->RAW)
my $ratioParamCode=($foldChangeRaw)? 'RATIO_RAW' : 'RATIO'; # raw ratio only for Software MaxQuant
#my $dispPvalue=(param('pValue'))? param('pValue') : 0.05; # undef at start
my $dispPvalue=param('pValue') || 0.05; $dispPvalue=1 if ($dispPvalue=~/[^\d\.e-]/ || $dispPvalue > 1);
#my $dispStdDev=param('stdDev') || 0; $dispStdDev=0 if $dispStdDev < 0; # fraction
#my $dispStdDevFrac=$dispStdDev/100;
my $dispCV=param('coefVar') // 0; $dispCV=0 if $dispCV=~/[^\d\.]/; # includes negative numbers
my $dispCVFrac=$dispCV/100;
my $dispPepType=param('pepType') || 'all';
my $numPepCode=($dispPepType eq 'all')? 'NUM_PEP_USED' : ($dispPepType eq 'distinct')? 'DIST_PEP_USED' : $dispPepType;
my $dispNumPep=param('numPep') || 3; $dispNumPep=3 if $dispNumPep=~/\D/;
my $dispSort=param('sort') || 'peptide'; # for list only
my (%dispRatios,%dispStates,%condToState,%ratioCode,%quantifParamInfo,$designID,$ratioType); # globals also for some ajax calls (%dispRatios & %dispStates mutually exclusve)
if (param('dispRatios')) {foreach my $pos (param('dispRatios')) {$dispRatios{$pos}=1;}} # Positions of ratios displayed in results
elsif (param('dispStates')) {foreach my $pos (param('dispStates')) {$dispStates{$pos}=1;}} # Positions of states displayed in results
my $dispMeasure=param('dispMeasure') || 'MQ_INT'; # for MaxQuant non-ratio quantif only
$dispSort='peptide' if (param('dispRatios') && $dispSort=~/ratio_(-*\d+)/ && !$dispRatios{$1});
my $restrictListID=param('restrictList') || 0;
if ($action eq 'delete') {&deleteQuantification($analysisID,$selQuantifID); exit;} # only for label-based 2ndary quantifications
elsif ($action eq 'ajaxDispCorrel') {&ajaxDisplayCorrelationMatrix; exit;}
elsif ($action eq 'ajaxProtStat') {&ajaxProteinStatistics; exit;}
elsif ($action eq 'ajaxPepRatio') {&ajaxPeptideRatios; exit;}
elsif ($action eq 'ajaxPepMaxQuant') {&ajaxPeptideMaxQuant; exit;}
elsif ($action eq 'ajaxPepSwath') {&ajaxPeptideSwath; exit;}
elsif ($action eq 'ajaxListProt') {&ajaxListSelectedProteins; exit;}
elsif ($action eq 'ajaxListDecGraph') {&ajaxListDecorateGraph; exit;}
elsif ($action eq 'ajaxGoTerms') {&ajaxGetGoTerms; exit;}
elsif ($action eq 'ajaxTermProt') {&ajaxGetGoTermProteins; exit;}
elsif ($action eq 'ajaxSaveProt') {&ajaxManageSaveProteins; exit;}
elsif ($action eq 'ajaxRestrictList') {&ajaxRestrictProteinList; exit;}
elsif ($action eq 'manageTemplate'){&ajaxManageTemplate; exit;}

my $dispResults=param('displayRes') || 0; # display form was submitted

####>Connect to the database
my $dbh=&promsConfig::dbConnect;
my $projectID=($selQuantifID)? &promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION') : &promsMod::getProjectID($dbh,$analysisID,'ANALYSIS');
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
$projectAccess=${$userInfo[2]}{$projectID};

###>Edition submission (internal-analysis only)<###
if ($action eq 'edit') {
	my $qName=$dbh->quote(param('name'));
	$dbh->do("UPDATE QUANTIFICATION SET NAME=$qName WHERE ID_QUANTIFICATION=$selQuantifID");
	$dbh->commit;
	$action='display';
}


################
####>Main 1<####
################
my (%anaQuantifList,%quantificationTypes,@quantificationInfo,%methodList);
my ($quantifType,$quantifSoftware,$algoVersion,$isPepIntensity,$quantifMethDesc,$title,$quantifModID,$modifName);
my ($focus,$quantifAnnot,$quantifStatus,$quantifMethodID,$updateDate,$updateUser);
my $highlighMatchStrg='';
my ($workbook,%itemFormat); # globals for export to Excel

if ($call eq 'ana') {
	$title=($action eq 'export')? 'Protein Quantification' : 'List of Single-Analysis Quantifications';
	my $sthAQ=$dbh->prepare("SELECT Q.ID_QUANTIFICATION,Q.NAME,Q.FOCUS,STATUS,M.ID_QUANTIFICATION_METHOD,M.NAME,M.CODE,M.DES,Q.QUANTIF_ANNOT FROM ANA_QUANTIFICATION A,QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE ID_ANALYSIS=$analysisID AND A.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_DESIGN IS NULL AND Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD");
	$sthAQ->execute;
	while (my ($quantifID,$name,$focus,$status,$methodID,$methodName,$methodCode,$methodDesc,$qAnnot)=$sthAQ->fetchrow_array) {
		$name='No name' unless $name;
		@{$anaQuantifList{$quantifID}}=($name,$focus,$status,$methodID,$qAnnot);
		@{$methodList{$methodID}}=($methodName,$methodCode,$methodDesc);
		my $typeName=($focus eq 'peptide' && $methodCode=~/XIC|SIN|SILAC|ITRAQ|TMT/)? 'Peptide quantification'
			: ($focus eq 'peptide' && $methodCode eq 'DIA')? 'SWATH-MS/DIA'
			: ($focus eq 'peptide')? 'Peptide ratio'
			: ($focus eq 'protein' && $methodCode=~/PROT_RATIO/)? 'Protein ratio' # should never be TNPQ (call=ana)
			: ($focus eq 'protein')? 'Protein quantification' # EMPAI
			: 'Unknown';
		push @{$quantificationTypes{$typeName}},$quantifID;
		#push @{$quantificationTypes{'Peptide quantification'}},$quantifID if ($designID && $methodCode=~/PROT_RATIO|TNPQ/);
	}
	$sthAQ->finish;
	$dbh->disconnect if $action eq 'select'; #'display'; # $selQuantifID;

	if ($selQuantifID) {
		($quantifType,$quantifMethDesc)=@{$methodList{$anaQuantifList{$selQuantifID}[3]}}[1,2];
		($focus,$quantifAnnot,$quantifStatus,$quantifMethodID,$updateDate,$updateUser)=$dbh->selectrow_array("SELECT FOCUS,QUANTIF_ANNOT,STATUS,ID_QUANTIFICATION_METHOD,UPDATE_DATE,UPDATE_USER FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	}
}
else { # called from quanti frame
	($quantifType,$designID,$quantifMethDesc,$quantifModID,$modifName,$focus,$quantifAnnot,$quantifStatus,$quantifMethodID,$updateDate,$updateUser)=
	$dbh->selectrow_array("SELECT QM.CODE,ID_DESIGN,QM.DES,M.ID_MODIFICATION,PSI_MS_NAME,FOCUS,QUANTIF_ANNOT,STATUS,QM.ID_QUANTIFICATION_METHOD,UPDATE_DATE,UPDATE_USER FROM QUANTIFICATION Q
			LEFT JOIN MODIFICATION M ON Q.ID_MODIFICATION=M.ID_MODIFICATION
			INNER JOIN QUANTIFICATION_METHOD QM ON Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD
			WHERE ID_QUANTIFICATION=$selQuantifID");
	$action='proteinView' if ($action eq 'select' && $designID);
	#@quantificationInfo=$dbh->selectrow_array("SELECT Q.ID_QUANTIFICATION,Q.NAME,Q.FOCUS,STATUS,M.ID_QUANTIFICATION_METHOD,M.NAME,M.CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=$selQuantifID");
	$title=($designID)? 'Design-based Quantification of ' : 'Single-Analysis Quantification of';
	if ($quantifModID) {
		$title.="<FONT color=\"#DD0000\">$modifName</FONT>-Proteins";
		$highlighMatchStrg=",'^###-'"; # for matching modProtID with protID in volcanoPlot
	}
	else {$title.='Whole Proteins';}
	if ($projectAccess=~/bioinfo|mass/) {
		$title.='&nbsp;<SPAN id="templateSpan">';
		if ($quantifAnnot=~/::TEMPLATE=1/) {
			$title.='<SPAN class="template">Template&nbsp;&nbsp;<INPUT type="button" class="font11" value="Remove" onclick="manageTemplate(\'remove\')"></SPAN>';
		}
		elsif ($quantifAnnot=~/SOFTWARE=DIA/ || ($quantifAnnot=~/SOFTWARE=myProMS;(\d+)/ && $1 >= 3)) {
			$title.='<INPUT type="button" class="title3" value="Save as template" onclick="manageTemplate(\'add\')">';
		}
		$title.='</SAPN>';
	}
}
$quantifMethDesc=&promsMod::resize($quantifMethDesc,50);
#$ratioParamCode='RATIO' if $quantifType ne 'MQ';
$quantifModID=0 unless $quantifModID;

if ($action eq 'export') {
	#################################
	####>Prepare export to Excel<####
	#################################
	my $timeStamp1=strftime("%Y%m%d %H-%M",localtime);
	#my $timeStamp2=strftime("%Y-%m-%d %H:%M",localtime);

	$workbook=Spreadsheet::WriteExcel->new("-");
	$workbook->set_properties(title=>"Protein quantification data",
							  author=>'myProMS server',
							  comments=>'Automatically generated with Perl and Spreadsheet::WriteExcel'
							  );
	$workbook->set_custom_color(40,189,215,255); # light light color  #BDD7FF (new myProMS colors V3.5+)
	#$workbook->set_custom_color(41,128,179,255); # dark blue color  #80B3FF
	%itemFormat=(
			title =>			$workbook->add_format(align=>'center',size=>18,bold=>1,border=>1),
			header =>			$workbook->add_format(align=>'center',size=>10,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			headerR =>			$workbook->add_format(align=>'right',valign=>'top',size=>10,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			mergeRowHeader =>	$workbook->add_format(align=>'center',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			mergeColHeader =>	$workbook->add_format(align=>'center',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			mergeColHeaderL =>	$workbook->add_format(align=>'left',valign=>'vcenter',size=>12,text_wrap=>1,bold=>1,bg_color=>40,border=>1),
			text =>				$workbook->add_format(align=>'left',size=>10,valign=>'top',text_wrap=>0,border=>1),
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

	print header(-type=>"application/vnd.ms-excel",-attachment=>"Protein_quantification_$timeStamp1.xls");

}
else {
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
.template {font-size:15px;color:#F5F5F5;background-color:#00A000;border-radius:5px;padding:2px 2px 2px 10px;}
.TD {font-weight:normal;}
.TH {font-weight:bold;}
.LINK {cursor:pointer;}
.highlight{width:250px;}
.popup {z-index:100;background-color:#FFFFFF;border:solid 3px #999999;padding:5px;position:absolute;display:none;box-shadow:10px 10px 20px #808080;}
</STYLE>
<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/heatMap.js"></SCRIPT>
|;
	if ($action eq 'summary') {
		print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
		&promsMod::popupInfo;
	}
	else {
		print qq
|<SCRIPT src="$promsPath{html}/js/local/genericPlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/volcanoPlot2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/peptidePlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
<SCRIPT LANGUAGE="JavaScript">
|;
	&promsMod::popupInfo;
	print qq
|var view='$view';
function displayQuantification(quantifInfoStrg,action) {
	var quantifInfo=quantifInfoStrg.split(':'); // focus:quantifID
	var focus=quantifInfo[0], quantifID=quantifInfo[1];
	if (quantifID==0) {action='select';}
	if (action=='delete' && !confirm('Delete selected quantification?')) {return;}
	if (focus=='protein') {
		//window.location="$promsPath{cgi}/showProtQuantification.cgi?CALL=$call&id_ana=$analysisID&id_quantif="+quantifID+"&ACT="+action+"&view="+view+"&foldChange=$dispFoldChange&pValue=$dispPvalue&sort=$dispSort"; //&stdDev=\$dispStdDev
		window.location="$promsPath{cgi}/showProtQuantification.cgi?CALL=$call&id_ana=$analysisID&id_quantif="+quantifID+"&ACT="+action+"&view=$view"; //&foldChange=\$dispFoldChange&pValue=\$dispPvalue&sort=$dispSort"; //&stdDev=\$dispStdDev
	}
	else { // peptide
		window.location="$promsPath{cgi}/showPepQuantification.cgi?CALL=$call&id_ana=$analysisID&id_quantif="+quantifID+"&ACT="+action;
	}
}
function sequenceView(id_protein,anaIdStrg) {
	var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+anaIdStrg+"&id_prot="+id_protein+"&msdata="+top.showMSdata;
	top.openProtWindow(winLocation);
}
var spectWin,selectedPep;
function drawSpectrum(pepId,pepInfo) {
	//if (pepId != selectedPep \|\| (spectWin && spectWin.closed)) {
		if (selectedPep && document.getElementById(selectedPep)) { // selectedPep could have been remove by another AJAX call
			document.getElementById(selectedPep).style.color='#000000';
		}
		if (pepId) {document.getElementById(pepId).style.color='#DD0000';} // null if called from SVG;
		selectedPep=pepId;
		var paramString="RID="+pepInfo+"&CALL=pep";
		spectWin=window.open("$promsPath{cgi}/drawSpectrum.cgi?"+paramString,'SpectrumWindow','width=950,height=950,location=no,resizable=yes,scrollbars=yes');
	//}
	spectWin.focus();
}
|;
	}
	print qq
|</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif"> <!-- Do not add onload code here! Use window.onload instead. -->
<CENTER>
<FONT class="title">$title</FONT>
<BR>
|;


#print "*** ACTION=$action ***<BR>\n";

	######################################################
	####>Internal quantifications drop-down selection<####
	######################################################
	if ($call eq 'ana') {
		print qq
|<TABLE bgcolor="$darkColor">
<TR><TD class="title2" align=right>Select:</TD><TD><SELECT class="title2" onchange="displayQuantification(this.value,'display')"><OPTION value="protein:0">-= Select =-</OPTION>
|;
		my $prevFocus='';
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
			if ($typeName eq 'Peptide quantification' && scalar @{$quantificationTypes{$typeName}} > 1) {
				print "<OPTION value=\"peptide:XIC\">Compare multiple XIC extractions</OPTION>";
			}
			print "</OPTGROUP>\n";
		}
		my $notDeletable=1;
		if ($selQuantifID) {
			($notDeletable)=$dbh->selectrow_array("SELECT COUNT(*) FROM PARENT_QUANTIFICATION WHERE ID_PARENT_QUANTIFICATION=$selQuantifID");
		}
		my $disabDeleteStrg=($notDeletable || $designID)? ' disabled' : '';

		print qq
|</SELECT></TD>
<TD><INPUT type="button" value="Edit" onclick="document.getElementById('editDIV').style.display='block'"/></TD>
<TD bgcolor="#DD0000"><INPUT type="button" value="Delete" style="color:#DD0000" onclick="displayQuantification('protein:$selQuantifID','delete')"$disabDeleteStrg/></TD>
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
<INPUT type="hidden" name="view" value="$view">
<INPUT type="hidden" name="id_ana" value="$analysisID">
<INPUT type="hidden" name="id_quantif" value="$selQuantifID">
<INPUT type="hidden" name="foldChgType" value="$foldChangeType">
<INPUT type="hidden" name="foldChange" value="$dispFoldChange">
<!--INPUT type="hidden" name="stdDev" value="\$dispStdDev"-->
<INPUT type="hidden" name="coefVar" value="$dispCV">
<INPUT type="hidden" name="pValue" value="$dispPvalue">
<INPUT type="hidden" name="pepType" value="$dispPepType">
<INPUT type="hidden" name="numPep" value="$dispNumPep">
<INPUT type="hidden" name="sort" value="$dispSort">
<DIV id="editDIV" style="display:none">
<TABLE bgcolor=$darkColor border=0>
<TR><TH align=right>&nbsp;Name :</TH><TD bgcolor=$lightColor><INPUT type="text" name="name" value="$anaQuantifList{$selQuantifID}[0]" size="50"></TD></TR>
<TR><TH colspan=2><INPUT type="submit" value=" Save ">&nbsp;&nbsp;&nbsp;<INPUT type="button" value=" Cancel " onclick="document.getElementById('editDIV').style.display='none'"/></TR>
</TABLE>
</DIV>
</FORM>
|;
	}
	#print "<DIV id=\"waitDiv\"><BR><BR><BR><BR><BR><FONT class=\"title3\">Fetching data. Please wait...</FONT><BR><IMG src=\"$promsPath{images}/scrollbarGreen.gif\"></DIV>\n";
}

###################################################
####>Fetching data for selected Quantification<####
###################################################
#($quantifType)=$dbh->selectrow_array("SELECT M.CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$selQuantifID") unless $quantifType; # internal quanti (SILAC, iTRAQ)
#if ($quantifType !~ /(PROT_RATIO|TNPQ|XIC|SIN|EMPAI)/) {
#my ($focus,$quantifAnnot,$quantifStatus,$quantifMethodID,$updateDate,$updateUser)=$dbh->selectrow_array("SELECT FOCUS,QUANTIF_ANNOT,STATUS,ID_QUANTIFICATION_METHOD,UPDATE_DATE,UPDATE_USER FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
if ($action=~/summary|export/) {
	if ($updateUser) {
		my $sthUser=$dbh->prepare("SELECT USER_NAME FROM USER_LIST WHERE ID_USER=?");
		$sthUser->execute($updateUser);
		my ($userName)=$sthUser->fetchrow_array;
		$updateUser=$userName || $updateUser;
	}
	else {$updateUser='?';}
}
if ($quantifType eq 'XIC' && $call eq 'quanti') { # To print the XIC in a good way
	$focus='protein';
	$action='xicView';
}
my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
my ($labelType)=($labelStrg)? ($labelStrg=~/LABEL=(.+)/) : ('FREE');
$labelType='FREE' unless $labelType;


###>Number of proteins and isoforms<###
my ($numProteinsQuantified,$numIsoFormsQuantified);
if ($quantifType !~ /EMPAI|SIN/) { # |MQ
	my $paramCode=($quantifType eq 'MQ')? $dispMeasure : $ratioParamCode; # 'NUM_PEP_TOTAL';
	my ($paramID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER WHERE CODE='$paramCode' AND ID_QUANTIFICATION_METHOD=$quantifMethodID");
	($numProteinsQuantified)=$dbh->selectrow_array("SELECT COUNT(DISTINCT ID_PROTEIN) FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID AND ID_QUANTIF_PARAMETER=$paramID");
	($numIsoFormsQuantified)=($quantifModID)? $dbh->selectrow_array("SELECT COUNT(ID_PROT_QUANTIF) FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID AND ID_QUANTIF_PARAMETER=$paramID") : 0;
}
#print "*** QUANTIF_TYPE='$quantifType' // LABEL_TYPE='$labelType' // FOCUS='$focus' // VIEW='$view' ***<BR>\n";
my $protQuantifiedStrg="$numProteinsQuantified protein"; $protQuantifiedStrg.='s' if $numProteinsQuantified && $numProteinsQuantified > 1;
if ($quantifModID) {
	$protQuantifiedStrg.=", $numIsoFormsQuantified $modifName-form"; $protQuantifiedStrg.='s' if $numIsoFormsQuantified > 1;
}
$protQuantifiedStrg.=' quantified';


##################################
####>Protein restriction list<####
##################################
my %restrictProteins;
if ($restrictListID) {
	my $useSites;
	if ($quantifModID) {
		($useSites)=$dbh->selectrow_array("SELECT 1 FROM CATEGORY WHERE ID_CATEGORY=$restrictListID AND LIST_TYPE='SITE'");
	}
	if ($useSites) {
		&promsQuantif::fetchSiteList($dbh,$restrictListID,\%restrictProteins,$quantifModID);
	}
	else { # Normal proteins
		my $sthRP=$dbh->prepare("SELECT DISTINCT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$restrictListID");
		$sthRP->execute;
		while (my($protID)=$sthRP->fetchrow_array) {
			$restrictProteins{$protID}=1;
		}
		$sthRP->finish;
	}
}

###################################
####>LABEL FREE/DESIGN QUANTIF<####
###################################
#if ($quantifType =~ /(PROT_RATIO|TNPQ|XIC|SIN|EMPAI)/){ #}
my $minPvalue=1;
if ($labelType eq 'FREE' || $designID) {

	################################
	####>TNPQ or PROT_RATIO_PEP<####
	################################
	#if (($designID && $quantifType ne 'MQ') || ($quantifType eq 'MQ' && $labelType eq 'SILAC')) {#} # multi-ana quantif
	if ($designID) { # multi-ana quantif

		####>View TNPQ/RATIO (Protein data)<####
		my ($usedAnaID)=($analysisID)? ($analysisID) : $dbh->selectrow_array("SELECT ANALYSIS.ID_ANALYSIS FROM ANALYSIS,ANA_QUANTIFICATION WHERE ANA_QUANTIFICATION.ID_ANALYSIS=ANALYSIS.ID_ANALYSIS AND ID_QUANTIFICATION IN (SELECT ID_PARENT_QUANTIFICATION FROM PARENT_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID) ORDER BY NAME ASC LIMIT 0,1");# Select one of the ANALYSIS related to that TNPQ quantification so as to make sequenceView method work...
#my ($quantifAnnot,$quantifMethodID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_QUANTIFICATION_METHOD FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
		#my ($parentQuantifType)=$dbh->selectrow_array("SELECT CODE FROM QUANTIFICATION,PARENT_QUANTIFICATION,QUANTIFICATION_METHOD WHERE QUANTIFICATION.ID_QUANTIFICATION=PARENT_QUANTIFICATION.ID_PARENT_QUANTIFICATION AND QUANTIFICATION.ID_QUANTIFICATION_METHOD=QUANTIFICATION_METHOD.ID_QUANTIFICATION_METHOD AND PARENT_QUANTIFICATION.ID_QUANTIFICATION=$selQuantifID LIMIT 1");
#my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
		my %labelingInfo;
		#($labelingInfo{'LABELTYPE'})=($labelStrg)? $labelStrg=~/LABEL=(.+)/ : ('FREE');
		my $sthQP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,NAME,CODE FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=$quantifMethodID");
		$sthQP->execute;
		while (my ($paramID,$paramName,$paramCode)=$sthQP->fetchrow_array) {
			@{$quantifParamInfo{$paramCode}}=($paramID,$paramName);
		}
		$sthQP->finish;
		foreach my $infoStrg (@labelInfo) {
			my ($setting,$valueStrg)=split('=',$infoStrg);
			if ($setting eq 'RATIOS' && $designID) { # $parentQuantifType && $parentQuantifType eq 'XIC'
				$valueStrg=~s/\#//g; # for design experiments, RATIOS are displayed with '#' for condition IDs
			}
			if ($setting eq 'PEPTIDES') {$valueStrg=~s/\#//g;} # remove ID tag of selected modifications if any
			@{$labelingInfo{$setting}}=split(';',$valueStrg);
		}

		##>Non-quantified proteins (myProMS v3)
		if ($labelingInfo{'LOST_PROTEINS'}) {
			my $proteoStrg=($quantifModID)? "$modifName-form" : 'protein';
			if ($labelingInfo{'LOST_PROTEINS'}[0]) {
				$protQuantifiedStrg.=" / $labelingInfo{LOST_PROTEINS}[0] $proteoStrg";
				$protQuantifiedStrg.=($labelingInfo{'LOST_PROTEINS'}[0] == 1)? ' was ' : 's were ';
				$protQuantifiedStrg.=' not quantified';
			}
			else { # 0
				$proteoStrg.='s';
				$protQuantifiedStrg.=" / All $proteoStrg were quantified";
			}
		}

		$ratioType=($labelingInfo{'RATIO_TYPE'})? $labelingInfo{'RATIO_TYPE'}[0] : (!$labelingInfo{'RATIOS'})? 'None' : 'Ratio';
$numPepCode='PEPTIDES' if ($ratioType eq 'None' && $numPepCode=~/_PEP_USED/); # (NUM/DIST)_PEP_USED not available in no-ratio MaxQuant
		my $numStates=scalar @{$labelingInfo{'STATES'}};
		$algoVersion=($ratioType eq 'Ratio')? 1 : (!$labelingInfo{'SOFTWARE'})? 2 : ($labelingInfo{'SOFTWARE'}[0] eq 'myProMS')? $labelingInfo{'SOFTWARE'}[1] : 0; # myProMS version ONLY!
		$isPepIntensity=($labelingInfo{'MEAN_STATE'} && $labelingInfo{'MEAN_STATE'}[0]==1)? 1 : ($algoVersion==2 && $labelType eq 'FREE')? 1 : ($algoVersion==3 && $labelingInfo{'ALGO_TYPE'}[0] eq 'PEP_INTENSITY')? 1 : 0; # 0 for v1 Ratio, v2 Label & v3 PEP_RATIO
		#$quantifSoftware=($algoVersion > 0)? "myProMS v$algoVersion" : ($labelingInfo{'SOFTWARE'}[0] eq 'MQ')? 'MaxQuant' : $labelingInfo{'SOFTWARE'}[0];
		$quantifSoftware=($algoVersion > 0)? 'myProMS' : ($labelingInfo{'SOFTWARE'}[0] eq 'MQ')? 'MaxQuant' : $labelingInfo{'SOFTWARE'}[0];
		if ($algoVersion) { # myProMS version ONLY!
			$quantifMethDesc=($isPepIntensity)?	'Protein ratio based on peptide intensity' : 'Protein ratio based on peptide ratio'; # overwrites $quantifMethDesc
		}
		%ratioCode=&generateRatioCodes($labelType,$ratioType,$numStates,$isPepIntensity,$algoVersion) if $ratioType ne 'None';

		###>Gene ontology annotation
		my %goAnalyses;
		if ($view eq 'graph' && $dispResults) {
			my $sthGO=$dbh->prepare("SELECT ID_GOANALYSIS,NAME,ASPECT FROM GO_ANALYSIS WHERE ID_EXPERIMENT=(SELECT D.ID_EXPERIMENT FROM DESIGN D,QUANTIFICATION Q WHERE D.ID_DESIGN=Q.ID_DESIGN AND Q.ID_QUANTIFICATION=$selQuantifID) AND GOA_TYPE=\"graph\"");
			$sthGO->execute;
			while (my ($goID,$name,$aspectStrg)=$sthGO->fetchrow_array) {
				@{$goAnalyses{$goID}}=($name,$aspectStrg);
			}
		}

		###>Displaying summary
		#my $numRatios=scalar @{$labelingInfo{'RATIOS'}};
		my %stateInfo=&printQuantificationSummary($dbh,$usedAnaID,\%labelingInfo,\%goAnalyses);

		#### ---> Form here on action cannot be 'summary' --->
		print "<DIV id=\"waitDiv\"><BR><BR><BR><BR><BR><FONT class=\"title3\">Fetching data. Please wait...</FONT><BR><IMG src=\"$promsPath{images}/scrollbarGreen.gif\"></DIV>\n" if $action ne 'export';


		###############################################
		#####>Fetching protein quantification data<####
		###############################################
		my (%quantifValues,%proteinInfo,%usedIsoforms);
		if (scalar keys %restrictProteins) {
			foreach my $protID (keys %restrictProteins) {@{$proteinInfo{$protID}}=();} # in case prot in restrict are not quantified (=>correct num of proteins displayed)
		}
		###my $sthProtQ=$dbh->prepare("SELECT ID_PROTEIN,TARGET_POS,QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID AND ID_QUANTIF_PARAMETER=?");
		my $sthProtQ=$dbh->prepare("SELECT PQ.ID_PROTEIN,GROUP_CONCAT(RESIDUE,POSITION ORDER BY POSITION SEPARATOR '.'),TARGET_POS,QUANTIF_VALUE
										FROM PROTEIN_QUANTIFICATION PQ
										LEFT JOIN PROTQUANTIF_MODRES PQMR ON PQ.ID_PROT_QUANTIF=PQMR.ID_PROT_QUANTIF
										LEFT JOIN MODIFIED_RESIDUE MR ON PQMR.ID_MODIF_RES=MR.ID_MODIF_RES AND MR.ID_QUANTIFICATION=PQ.ID_QUANTIFICATION
										WHERE PQ.ID_QUANTIFICATION=$selQuantifID AND ID_QUANTIF_PARAMETER=? GROUP BY PQ.ID_PROT_QUANTIF");
		my $pvalueCode;
		my @params;
		if ($quantifSoftware eq 'MaxQuant') {
			$pvalueCode='';
			if ($ratioType eq 'None') { # Intensity quantifications
				#@params=($view eq 'graph')? ($dispMeasure) : ('MQ_LFQ','MQ_INT','MQ_IBAQ','MQ_SC');
				@params=($dispMeasure);
			}
			else {
				@params=($ratioParamCode,'RATIO_VAR'); # filtering can apply on RATIO_VAR in list view
			}
		}
		else {
			$pvalueCode=($ratioType=~/S\w+Ratio/ || $labelingInfo{'FDR_CONTROL'}[0] eq 'TRUE')? 'PVAL_ADJ' : 'PVAL';
			@params=('RATIO',$pvalueCode);
			push @params,'SD_GEO' if ($quantifType eq 'PROT_RATIO_PEP' && $quantifSoftware ne 'MaxQuant'); # no SD_GEO for TNPQ!!!
			push @params,('CONF_LOW','CONF_UP') if ($action eq 'export' && $ratioType eq 'Ratio');
			#push @params,'DIST_PEP_USED' if ($ratioType=~/S\w+Ratio/ && ($numPepCode eq 'DIST_PEP_USED' || $action eq 'export')); # for all ratios
		}
		foreach my $paramCode (@params) {
			$sthProtQ->execute($quantifParamInfo{$paramCode}[0]);
			while (my ($protID,$modResStrg,$targetPos,$qValue)=$sthProtQ->fetchrow_array) {
				#next unless $dispRatios{$targetPos};
				if ($ratioType eq 'None') { # no-ratio MaxQuant
					next unless $dispStates{$targetPos};
				}
				else {
					if ($dispRatios{-$targetPos}) { # reverse ratio management
						$targetPos=-$targetPos;
						$qValue=1/$qValue if $paramCode=~/^RATIO/;
					}
					elsif (!$dispRatios{$targetPos}) {next;}
				}
				my $protID0=$protID;
				if ($modResStrg) { # quantif of modification
					$protID.='-'.&promsQuantif::decodeModificationSite($modResStrg);
				}
				next if ($restrictListID && !$restrictProteins{$protID});
				@{$proteinInfo{$protID0}}=(); # actual number of quantified proteins (only list filtering)
				#$quantifValues{$protID}{$paramCode}{$targetPos}=($paramCode eq 'SD_GEO')? $qValue*100 : $qValue;
				$quantifValues{$protID}{$paramCode}{$targetPos}=$qValue;
				$minPvalue=$qValue if ($paramCode eq $pvalueCode && $qValue > 0 && $qValue < $minPvalue);
				$usedIsoforms{$protID}=1;
			}
		}
		$minPvalue/=10;
		my @pepParams;
		if ($ratioType=~/S\w+Ratio/) { # No peptide data for super ratios before 18/02/15: use those from linked primary ratios
# TODO: Simplify (remove linkedRatios) to make compatible with only quantif >18/02/15
			my %linkedRatios; # for SuperRatio only
if ($ratioType eq 'SuperRatio') {
			my $rPos=0;
			foreach my $ratioData (@{$labelingInfo{'RATIOS'}}) { # Find linked primary ratios
				$rPos++;
				my ($testStatePos,$refStatePos)=split(/\//,$ratioData);
				if ($testStatePos=~/%(\d+)/) { # ratio of ratio
					#next unless $dispRatios{$rPos};
					next if (!$dispRatios{$rPos} && !$dispRatios{-$rPos});
					my $linkedTestPos=$1;
					my ($linkedRefPos)=($refStatePos=~/%(\d+)/);
					$testStatePos=~s/%\d+//;
					$refStatePos=~s/%\d+//;
					push @{$linkedRatios{$linkedTestPos}},$rPos;
					push @{$linkedRatios{$linkedRefPos}},$rPos;
				}
			}
}
			#my @pepSupParams=('NUM_PEP_USED'); # NUM_PEP_USED always needed (point size in graph, displayed in list)
			#push @pepSupParams,'DIST_PEP_USED' if ($numPepCode eq 'DIST_PEP_USED' || ($quantifType ne 'MQ' && $action eq 'export')); # for all ratios
			my @pepSupParams; # channel-specific params
			@pepSupParams=($numPepCode) if $numPepCode=~/(NUM|DIST)_PEP_USED/;
			if ($view eq 'list' || $action eq 'export') {
				push @pepSupParams,'NUM_PEP_USED' if $numPepCode ne 'NUM_PEP_USED';
				push @pepSupParams,'DIST_PEP_USED' if ($numPepCode ne 'DIST_PEP_USED' && $quantifSoftware ne 'MaxQuant');
			}
			foreach my $paramCode (@pepSupParams) {
				$sthProtQ->execute($quantifParamInfo{$paramCode}[0]);
				while (my ($protID,$modResStrg,$targetPos,$qValue)=$sthProtQ->fetchrow_array) {
					next unless $proteinInfo{$protID}; # $linkedRatios{$targetPos} could bring additional protID NOT in $proteinInfo{$protID}
					my $ratioPos0=$targetPos; # reverse ratio management
					#next if (!$dispRatios{$targetPos} && !$linkedRatios{$targetPos});
					if ($ratioType eq 'None') { # no-ratio MaxQuant
						next unless $dispStates{$targetPos};
					}
					else {
						$targetPos=-$targetPos if $dispRatios{-$targetPos};
						next if (!$dispRatios{$targetPos} && !$linkedRatios{$ratioPos0});
						#next if ($restrictListID && !$restrictProteins{$protID});
						if ($modResStrg) { # quantif of modification
							$protID.='-'.&promsQuantif::decodeModificationSite($modResStrg);
							next unless $usedIsoforms{$protID}; # $linkedRatios{$targetPos} could bring additional isoforms
						}
					}
					$quantifValues{$protID}{$paramCode}{$targetPos}=$qValue if $dispRatios{$targetPos};
					if ($linkedRatios{$ratioPos0}) { # SuperRatio
						foreach my $rPos (@{$linkedRatios{$ratioPos0}}) {
							#$quantifValues{$protID}{$paramCode}{$rPos}+=$qValue;
							##$quantifValues{$protID}{$paramCode}{$rPos}=$qValue if ($quantifValues{$protID}{$paramCode} && $quantifValues{$protID}{$paramCode}{$rPos} && $quantifValues{$protID}{$paramCode}{$rPos} > $qValue); # keep smallest value
							$quantifValues{$protID}{$paramCode}{$rPos}=$qValue if ($dispRatios{$rPos} && (!$quantifValues{$protID}{$paramCode} || !$quantifValues{$protID}{$paramCode}{$rPos} || $quantifValues{$protID}{$paramCode}{$rPos} > $qValue)); # keep smallest value
							$quantifValues{$protID}{$paramCode}{-$rPos}=$qValue if ($dispRatios{-$rPos} && (!$quantifValues{$protID}{$paramCode} || !$quantifValues{$protID}{$paramCode}{-$rPos} || $quantifValues{$protID}{$paramCode}{-$rPos} > $qValue)); # keep smallest value
						}
					}
				}
			}
			push @pepParams,$numPepCode if $numPepCode !~ /(NUM|DIST)_PEP_USED/; # MaxQuant PEPTIDES, RAZ_UNI_PEP or UNIQUE_PEP
			if ($view eq 'list' || $action eq 'export') {
				if ($quantifSoftware eq 'MaxQuant') {
					foreach my $param ('PEPTIDES','RAZ_UNI_PEP','UNIQUE_PEP') {
						push @pepParams,$param if $param ne $numPepCode;
					}
				}
				else {@pepParams=('NUM_PEP_TOTAL');}
			}
#print "PEP='@pepParams', SUP='@pepSupParams' SEL=$numPepCode,$quantifParamInfo{$numPepCode}[0]<BR>\n";
		}
		else {
			if ($view eq 'graph') {@pepParams=($numPepCode);}
			else {
				if ($ratioType eq 'None') {
					@pepParams=('PEPTIDES','RAZ_UNI_PEP','UNIQUE_PEP')
				}
				else {
					@pepParams=($numPepCode eq 'DIST_PEP_USED' || ($quantifSoftware ne 'MaxQuant' && $action eq 'export'))? ('DIST_PEP_USED','NUM_PEP_USED','NUM_PEP_TOTAL') : ('NUM_PEP_USED','NUM_PEP_TOTAL');
				}
			}
		}

		my $refDispSets=($ratioType eq 'None')? \%dispStates : \%dispRatios;
		foreach my $paramCode (@pepParams) { # empty for SuperRatio in graph view
			$sthProtQ->execute($quantifParamInfo{$paramCode}[0]);
			while (my ($protID,$modResStrg,$targetPos,$qValue)=$sthProtQ->fetchrow_array) {
				next unless $proteinInfo{$protID}; # NUM_PEP_TOTAL has no ratioPos & could bring additional protID
				#next if ($restrictListID && !$restrictProteins{$protID});
				if ($modResStrg) { # quantif of modification
					$protID.='-'.&promsQuantif::decodeModificationSite($modResStrg);
					next unless $usedIsoforms{$protID};
				}
				#next if ($targetPos && !$dispRatios{$targetPos}); # no ratioPos for NUM_PEP_TOTAL
				if ($targetPos) { # no ratioPos for NUM_PEP_TOTAL or MQ intensities
					if ($ratioType eq 'None') { # no-ratio MaxQuant
						next unless $dispStates{$targetPos};
					}
					else {
						if ($dispRatios{-$targetPos}) {$targetPos=-$targetPos;} # reverse ratio management
						elsif (!$dispRatios{$targetPos}) {next;}
					}
					$quantifValues{$protID}{$paramCode}{$targetPos}=$qValue if $refDispSets->{$targetPos}; # NUM_PEP_TOTAL has no ratioPos & could bring additional isoforms
				}
				else { # undef for quantified proteins
					foreach my $tPos (keys %{$refDispSets}) { # extend to all displayed ratios ($targetPos not always defined for single-ratio quantif)
						$quantifValues{$protID}{$paramCode}{$tPos}=$qValue;
					}
				}
			}
		}

		$sthProtQ->finish;

		###>Fetching protein info
		#my %proteinInfo;
		#my $extraFieldStrg=($view eq 'graph')? ',MW,PROT_DES,ORGANISM,PROT_LENGTH' : ',MW,PROT_DES,ORGANISM,PROT_LENGTH';
		#my $extraFieldStrg=($view eq 'graph')? '' : ',MW,PROT_DES,ORGANISM';
		my $sthPI=$dbh->prepare("SELECT DISTINCT P.ID_PROTEIN,ALIAS,MW,PROT_DES,ORGANISM,PROT_LENGTH FROM PROTEIN_QUANTIFICATION Q,PROTEIN P WHERE Q.ID_PROTEIN=P.ID_PROTEIN AND ID_QUANTIFICATION=$selQuantifID");
		my ($GNidentID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
		my $sthGN=$dbh->prepare("SELECT MI.VALUE FROM PROTEIN P,MASTERPROT_IDENTIFIER MI WHERE P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND MI.ID_IDENTIFIER=$GNidentID AND P.ID_PROTEIN=? ORDER BY MI.RANK");
		#my $sthPA=$dbh->prepare("SELECT AP.ID_ANALYSIS FROM ANALYSIS_PROTEIN AP,ANA_QUANTIFICATION AQ WHERE AP.ID_ANALYSIS=AQ.ID_ANALYSIS AND ID_PROTEIN=? AND AQ.ID_QUANTIFICATION=$selQuantifID ORDER BY VISIBILITY DESC,NUM_PEP DESC LIMIT 0,1");# To be sure that the PROTEIN is in the ANALYSIS
		my ($allAnaStrg)=$dbh->selectrow_array("SELECT GROUP_CONCAT(ID_ANALYSIS) FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID GROUP BY ID_QUANTIFICATION");
		my $sthPA=$dbh->prepare("SELECT GROUP_CONCAT(ID_ANALYSIS) FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS IN ($allAnaStrg) AND ID_PROTEIN=? AND VISIBILITY >= 1");# To be sure that the PROTEIN is in the ANALYSIS

		$sthPI->execute;
		while (my ($protID,@info)=$sthPI->fetchrow_array) {
			next unless $proteinInfo{$protID};
			#next if ($restrictListID && !$restrictProteins{$protID});
			$info[0]=~s/ .*//; # clean alias from badly parsed characters (identifier conversion)
			my @geneList;
			if ($view eq 'list' || $action eq 'export') {
				$info[1]=sprintf "%.1f",$info[1]/1000;
				$sthGN->execute($protID);
				while (my ($gene)=$sthGN->fetchrow_array) {push @geneList,$gene;}
			}
			my $anaStrg;
			if ($view eq 'list') {
				$sthPA->execute($protID);
				($anaStrg)=$sthPA->fetchrow_array;
				$anaStrg=0 unless $anaStrg;
			}
			else {$anaStrg=$allAnaStrg;}
			@{$proteinInfo{$protID}}=(@info,\@geneList,$anaStrg); # [5] for genelist
			#$proteinInfo{$protID}[1]=sprintf "%.1f",$proteinInfo{$protID}[1]/1000 if $view eq 'list'; # MW
		}
		$sthPI->finish;
		$sthGN->finish;
		$sthPA->finish;

		$dbh->disconnect;

		####################
		####>Export XLS<####
		####################
		if ($action eq 'export') {
			&exportProteinList($labelType,$ratioType,\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,$pvalueCode);
			$workbook->close();
			exit;
		}
		#print qq |<DIV id="waitDiv"><BR><BR><BR><BR><BR><FONT class="title3">Fetching data. Please wait...</FONT><BR><IMG src="$promsPath{images}/scrollbarGreen.gif"></DIV><SCRIPT type="text/javascript">document.getElementById('waitDiv').style.display='none'</SCRIPT>|;
		print qq |<SCRIPT type="text/javascript">document.getElementById('waitDiv').style.display='none'</SCRIPT>|;
		#my $numRatios=scalar @{$labelingInfo{RATIOS}};

		########################
		####>Graphical view<####
		########################
		if ($view eq 'graph') {
			if ($quantifSoftware eq 'MaxQuant') {
				if ($ratioType eq 'None') {
					#&displayGenericPlot(\%quantifValues,\%proteinInfo,'MQ_INT','Intensity',$usedAnaID);
					&displayIntensityPlot(\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo);
				}
				else {
					# Ratios but no p-values!!!
					&displayRatioPlot(\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo);
				}
			}
			else {
				&displayVolcanoPlot(\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,$pvalueCode);
			}
		}

		###################
		####>List view<####
		###################
		else { # View type is list
			print qq
|<SCRIPT type="text/javascript">
window.onload=function() {
	ajaxUpdateRestrict();
}
</SCRIPT>
|;
			if ($quantifSoftware eq 'MaxQuant') {
				&printMaxQuantProteinList($labelType,\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo);
			}
			else {
				&printProteinRatioList($labelType,$ratioType,\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,$pvalueCode); #$labelingInfo{FDR_ALPHA}[0]/100,$dispStdDev
			}
		}
	}

	######################
	####>EMPAI or SIN<####
	######################
	elsif ($quantifType =~ /EMPAI|SIN/) { # |MQ

		print "<DIV id=\"waitDiv\"><BR><BR><BR><BR><BR><FONT class=\"title3\">Fetching data. Please wait...</FONT><BR><IMG src=\"$promsPath{images}/scrollbarGreen.gif\"></DIV>\n" if $action ne 'export';

		my @qpcodes=($quantifType eq 'EMPAI')? ('EMPAI','EMPAI_MOL','EMPAI_MR') : ('SIN_SIN'); #:($quantifType eq 'MQ')? ('MQ_LFQ','MQ_INT','MQ_IBAQ','MQ_SC')
		my ($usedAnaID)=($analysisID)? ($analysisID) : $dbh->selectrow_array("SELECT ID_ANALYSIS FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID"); # internal quantif
		#my $sthProtQ=$dbh->prepare("SELECT ID_PROTEIN,QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID AND QUANTIF_VALUE IS NOT NULL");
		my $sthProtQ=$dbh->prepare("SELECT ID_PROTEIN,QUANTIF_VALUE,QP.CODE FROM PROTEIN_QUANTIFICATION PQ, QUANTIFICATION_PARAMETER QP WHERE PQ.ID_QUANTIF_PARAMETER=QP.ID_QUANTIF_PARAMETER AND ID_QUANTIFICATION=$selQuantifID AND QUANTIF_VALUE IS NOT NULL");
		#my $sthProtI=$dbh->prepare("SELECT P.ID_PROTEIN,ALIAS,MW,PROT_DES,ORGANISM,NUM_PEP FROM PROTEIN_QUANTIFICATION Q,PROTEIN P WHERE Q.ID_PROTEIN=P.ID_PROTEIN AND ID_QUANTIFICATION=$selQuantifID");
		my $sthProtI=$dbh->prepare("SELECT P.ID_PROTEIN,ALIAS,NUM_PEP,MW,PROT_DES,ORGANISM FROM PROTEIN P,ANALYSIS_PROTEIN AP WHERE AP.ID_PROTEIN=P.ID_PROTEIN AND VISIBILITY>=1 AND ID_ANALYSIS=$usedAnaID");
		my ($GNidentID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
		my $sthGN=$dbh->prepare("SELECT MI.VALUE FROM PROTEIN P,MASTERPROT_IDENTIFIER MI WHERE P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND MI.ID_IDENTIFIER=$GNidentID AND P.ID_PROTEIN=? ORDER BY MI.RANK");

		$sthProtI->execute;
		my (%proteinInfo,%quantifValues); #,@protOrder
		while (my ($protID,@info)=$sthProtI->fetchrow_array) {
			next if ($restrictListID && !$restrictProteins{$protID});
			@{$proteinInfo{$protID}}=@info;
			$proteinInfo{$protID}[2]=sprintf "%.1f",$proteinInfo{$protID}[2]/1000;# if $view eq 'list'; # MW
			my @geneList;
			$sthGN->execute($protID);
			while (my ($gene)=$sthGN->fetchrow_array) {push @geneList,$gene;}
			push @{$proteinInfo{$protID}},\@geneList;
			foreach my $paramC (@qpcodes) { $quantifValues{$protID}{$paramC}=0; }# for later sort
		}
		$sthProtI->finish;
		$sthGN->finish;
		$sthProtQ->execute;
		while (my ($protID,$qValue,$code)=$sthProtQ->fetchrow_array) {
			next unless $proteinInfo{$protID}; # skip hidden & restricted protein even if they were quantified!
			#push @protOrder,$protID;
			($quantifValues{$protID}{$code})=$qValue;
		}
		$sthProtQ->finish;
		my $numQuantProteins=scalar keys %quantifValues;
		#&displayIntensityPlot(\%quantifValues,\%proteinInfo,'MQ_INT','Intensity',$usedAnaID) if ($quantifType eq 'MQ' && $action ne 'export');
		my $locationStrg='';
		if ($action eq 'export') {
			my @quantifInfo=&promsMod::getItemInfo($dbh,'QUANTIFICATION',$selQuantifID);
			foreach my $it (@quantifInfo) {
				$locationStrg.=' > ' if $locationStrg;
				$locationStrg.=$it->{'NAME'};
			}
		}

		$dbh->disconnect;

		my $numTotProteins=scalar keys %proteinInfo;
		my $proteinTitle="$numQuantProteins/$numTotProteins Proteins";
		my $colName=($quantifType eq 'EMPAI')? 'emPAI' : ($action eq 'export')? 'SIN' : 'SI<SUB>N</SUB>'; #:($quantifType eq 'MQ')? 'MaxQuant'  <BR>extracted  <BR>computed
		my $nbCol=($quantifType eq 'EMPAI')? 8 : 6;
		#my ($lightColor,$darkColor)=&promsConfig::getRowColors;

		if ($action eq 'export') {
			####<Start printing>####
			my $worksheet=$workbook->add_worksheet('Results');

			##>Table header<##
			my $xlsRow=0;
			my $xlsCol=0;
			$worksheet->set_column(0,0,30); # identifier col length
			$worksheet->merge_range($xlsRow,$xlsCol,$xlsRow,$xlsCol+$nbCol,'myProMS '.$colName.' quantification : '.$locationStrg,$itemFormat{'mergeColHeaderL'});
			$worksheet->write_string(++$xlsRow,$xlsCol,$proteinTitle,$itemFormat{'header'});
			$worksheet->write_string($xlsRow,++$xlsCol,'Gene & Synonyms',$itemFormat{'header'});
			if ($quantifType eq 'EMPAI'){
				$worksheet->write_string($xlsRow,++$xlsCol,'emPAI',$itemFormat{'header'});
				$worksheet->write_string($xlsRow,++$xlsCol,'emPAI (Mol %)',$itemFormat{'header'});
				$worksheet->write_string($xlsRow,++$xlsCol,'emPAI (Mr %)',$itemFormat{'header'});
			}
			#elsif ($quantifType eq 'MQ'){
			#	$worksheet->write_string($xlsRow,++$xlsCol,'LFQ',$itemFormat{'header'});
			#	$worksheet->write_string($xlsRow,++$xlsCol,'Intensity',$itemFormat{'header'});
			#	$worksheet->write_string($xlsRow,++$xlsCol,'iBAQ',$itemFormat{'header'});
			#	$worksheet->write_string($xlsRow,++$xlsCol,'Spectral Count',$itemFormat{'header'});
			#}
			else{
				$worksheet->write_string($xlsRow,++$xlsCol,$colName,$itemFormat{'header'});
			}
			$worksheet->write_string($xlsRow,++$xlsCol,'Peptides',$itemFormat{'header'});
			$worksheet->write_string($xlsRow,++$xlsCol,'MW (kDa)',$itemFormat{'header'});
			$worksheet->set_column(++$xlsCol,$xlsCol,80); # col length
			$worksheet->write_string($xlsRow,$xlsCol,'Description',$itemFormat{'header'});
			$worksheet->set_column(++$xlsCol,$xlsCol,30); # col length
			$worksheet->write_string($xlsRow,$xlsCol,'Species',$itemFormat{'header'});

			##>Protein list
			foreach my $protID (sort{$quantifValues{$b}{$qpcodes[0]}<=>$quantifValues{$a}{$qpcodes[0]}} keys %proteinInfo) { # @protOrder
				$xlsCol=0;
				#>Identifier
				$worksheet->write_string(++$xlsRow,$xlsCol,$proteinInfo{$protID}[0],$itemFormat{'text'});
				#>Gene
				$worksheet->write_string($xlsRow,++$xlsCol,join(',',@{$proteinInfo{$protID}[5]}),$itemFormat{'text'});
				#>Quantif value
				foreach my $paramC (@qpcodes) {
					if ($quantifValues{$protID}{$paramC}) {$worksheet->write_number($xlsRow,++$xlsCol,$quantifValues{$protID}{$paramC},$itemFormat{'number'});}
					else {$worksheet->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
				}
				#>Num peptides
				$worksheet->write_number($xlsRow,++$xlsCol,$proteinInfo{$protID}[1],$itemFormat{'number'});
				#>MW
				$worksheet->write_number($xlsRow,++$xlsCol,$proteinInfo{$protID}[2],$itemFormat{'number1d'});
				#>Des
				$worksheet->write_string($xlsRow,++$xlsCol,$proteinInfo{$protID}[3],$itemFormat{'textWrap'});
				#>Species
				$worksheet->write_string($xlsRow,++$xlsCol,$proteinInfo{$protID}[4],$itemFormat{'text'});
			}
			$workbook->close();
			exit;
		}
		else {
			print qq
|<SCRIPT type="text/javascript">document.getElementById('waitDiv').style.display='none'</SCRIPT>
<TABLE border=0 cellspacing=0 cellpadding=2>
<TR><TD colspan=5><INPUT type="button" value="Export data" onclick="window.location='./showProtQuantification.cgi?CALL=$call&id_ana=$analysisID&id_quantif=$selQuantifID&ACT=export'"/></TD></TR>
<TR bgcolor="$darkColor">
<TH class="rbBorder" align=left nowrap>&nbsp;$proteinTitle&nbsp;</TH><TH class="rbBorder">&nbsp;Gene&nbsp;</TH>
|;

#			if ($quantifType eq 'MQ') {
#				print qq
#|<TH class="rbBorder" nowrap>&nbsp;LFQ&nbsp;</TH>
#<TH class="rbBorder" nowrap>&nbsp;Instensity&nbsp;</TH>
#<TH class="rbBorder" nowrap>&nbsp;iBAQ&nbsp;</TH>
#<TH class="rbBorder" nowrap>&nbsp;Spectral Count&nbsp;</TH>
#|;
#			}
#			elsif
			if ($quantifType eq 'EMPAI') {
				print qq
|<TH class="rbBorder" nowrap>&nbsp;emPAI&nbsp;</TH>
<TH class="rbBorder" nowrap>&nbsp;emPAI (Mol %)&nbsp;</TH>
<TH class="rbBorder" nowrap>&nbsp;emPAI (Mr %)&nbsp;</TH>
|;
			}
			else {print "<TH class=\"rbBorder\" nowrap>&nbsp;$colName&nbsp;</TH>";}
			print qq
|<TH class="rbBorder" nowrap>&nbsp;Peptides&nbsp;</TH>
<TH class="rbBorder" nowrap>&nbsp;MW<SMALL> kDa</SMALL>&nbsp;</TH>
<TH class="bBorder" width=700>&nbsp;Description - Species&nbsp;</TH>
</TR>
|;
			my $bgColor=$lightColor;
			foreach my $protID (sort{$quantifValues{$b}{$qpcodes[0]}<=>$quantifValues{$a}{$qpcodes[0]}} keys %proteinInfo) { # @protOrder
				print "<TR valign=\"top\" bgcolor=\"$bgColor\">";
				print "<TD class=\"TH\" nowrap=\"nowrap\">&nbsp;<A href=\"javascript:sequenceView($protID,'$usedAnaID')\">$proteinInfo{$protID}[0]</A>&nbsp;</TH>"; # alias
				if (scalar @{$proteinInfo{$protID}[5]} > 1) { # gene
					print "<TH>&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-",join('<BR>&nbsp;&nbsp;-',@{$proteinInfo{$protID}[5]}[1..$#{$proteinInfo{$protID}[5]}]),"</B>')\" onmouseout=\"popout()\">",$proteinInfo{$protID}[5][0],"</A>&nbsp;</TH>";
				}
				elsif ($proteinInfo{$protID}[5][0]) {print "<TH>&nbsp;$proteinInfo{$protID}[5][0]&nbsp;</TH>";}
				else {print '<TH>&nbsp;-&nbsp;</TH>';} # no gene mapped
				foreach my $paramC (@qpcodes) {
					if ($quantifValues{$protID}{$paramC}) {
						if($quantifValues{$protID}{$paramC}<1000){
							printf "<TH>&nbsp;$quantifValues{$protID}{$paramC}&nbsp;</TH>";
						}
						else{
							printf "<TH>&nbsp;%.2e&nbsp;</TH>",$quantifValues{$protID}{$paramC};
						}
					}
					else {print "<TH>-</TH>"};
				}
				print "<TH>&nbsp;$proteinInfo{$protID}[1]&nbsp;</TH><TD class=\"TH\" align=right>$proteinInfo{$protID}[2]&nbsp;</TD><TD>$proteinInfo{$protID}[3] <U><I>$proteinInfo{$protID}[4]</I></U></TD>";
				print "</TR>\n";
				$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
			}
			print "<TR><TD colspan=5><B>End of list.</B></TD></TR>\n</TABLE>\n";
		}
	} # end of SIN/EMPAI HTML display
	&endHTML;
}

##################################
####>INTERNAL LABELED QUANTIF<#### Cannot be MaxQuant
##################################
else { #if ($focus eq 'protein') { #} protein

	if ($call eq 'quanti') { # possible since openProject v1.5.9
		($analysisID,@{$anaQuantifList{$selQuantifID}})=$dbh->selectrow_array("SELECT A.ID_ANALYSIS,Q.NAME,Q.FOCUS,STATUS,M.ID_QUANTIFICATION_METHOD FROM ANA_QUANTIFICATION A,QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE A.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND Q.ID_QUANTIFICATION=$selQuantifID");
	}
	my %labelingInfo;
	#my $quantifType=$methodList{$anaQuantifList{$selQuantifID}[3]}[1];
	my $sthQP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,NAME,CODE FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=$anaQuantifList{$selQuantifID}[3]");
	$sthQP->execute;
	while (my ($paramID,$paramName,$paramCode)=$sthQP->fetchrow_array) {
		@{$quantifParamInfo{$paramCode}}=($paramID,$paramName);
	}
	$sthQP->finish;

	########################
	####>PROT_RATIO_PEP<####
	########################
	if ($quantifType eq 'PROT_RATIO_PEP') { # || $quantifType eq 'TNPQ'
		foreach my $infoStrg (@labelInfo) {
			my ($setting,$valueStrg)=split('=',$infoStrg);
			@{$labelingInfo{$setting}}=split(';',$valueStrg);
		}
		$ratioType=($labelingInfo{'RATIO_TYPE'})? $labelingInfo{'RATIO_TYPE'}[0] : 'Ratio';
		$quantifSoftware='myProMS';
		$algoVersion=($ratioType eq 'Ratio')? 1 : 2;
		my $numStates=(scalar @{$labelingInfo{'STATES'}}); #-1;
		foreach my $ratioPos (1..$numStates*($numStates-1)/2) {$ratioCode{$ratioPos}=$ratioPos; $ratioCode{-$ratioPos}=-$ratioPos;}

		###>Gene ontology annotation
		my %goAnalyses;
		if ($view eq 'graph' && $dispResults) {
			my @anaParents=&promsMod::getItemInfo($dbh,'ANALYSIS',$analysisID);
			#my $sthGO=$dbh->prepare("SELECT GO.ID_GOANALYSIS,NAME,ASPECT FROM GO_ANALYSIS GO,GOANA_ANALYSIS GA WHERE GO.ID_GOANALYSIS=GA.ID_GOANALYSIS AND ID_ANALYSIS=$analysisID");
			my $sthGO=$dbh->prepare("SELECT ID_GOANALYSIS,NAME,ASPECT FROM GO_ANALYSIS WHERE ID_EXPERIMENT=$anaParents[1]{ID} AND GOA_TYPE=\"graph\"");
			$sthGO->execute;
			while (my ($goID,$name,$aspectStrg)=$sthGO->fetchrow_array) {
				@{$goAnalyses{$goID}}=($name,$aspectStrg);
			}
		}

		###>Displaying summary
		#my $numRatios=scalar @{$labelingInfo{'RATIOS'}};
		my %stateInfo=&printQuantificationSummary($dbh,$analysisID,\%labelingInfo,\%goAnalyses);

		###>Error in Quantification
		if ($action ne 'export') {
			if ($anaQuantifList{$selQuantifID}[2] == -1) {
				$dbh->disconnect;
				print qq
|<SCRIPT type="text/javascript">document.getElementById('waitDiv').style.display='none'</SCRIPT>
<BR><FONT class=\"title2\" color='#DD0000'>***An error has occured during quantification***</FONT>
</CENTER>
</BODY>
</HTML>
|;
				exit;
			}
			print "<DIV id=\"waitDiv\"><BR><BR><BR><BR><BR><FONT class=\"title3\">Fetching data. Please wait...</FONT><BR><IMG src=\"$promsPath{images}/scrollbarGreen.gif\"></DIV>\n";
		}

		###>Fetching protein quantification data
		my (%quantifValues,%proteinInfo,%usedIsoforms);
		if (scalar keys %restrictProteins) {
			foreach my $protID (keys %restrictProteins) {@{$proteinInfo{$protID}}=();} # in case prot in restrict are not quantified (=>correct num of proteins displayed)
		}
		#my $sthProtQ=$dbh->prepare("SELECT ID_PROTEIN,TARGET_POS,QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID AND ID_QUANTIF_PARAMETER=?");
		my $sthProtQ=$dbh->prepare("SELECT PQ.ID_PROTEIN,GROUP_CONCAT(RESIDUE,POSITION ORDER BY POSITION SEPARATOR '.'),TARGET_POS,QUANTIF_VALUE
										FROM PROTEIN_QUANTIFICATION PQ
										LEFT JOIN PROTQUANTIF_MODRES PQMR ON PQ.ID_PROT_QUANTIF=PQMR.ID_PROT_QUANTIF
										LEFT JOIN MODIFIED_RESIDUE MR ON PQMR.ID_MODIF_RES=MR.ID_MODIF_RES AND MR.ID_QUANTIFICATION=PQ.ID_QUANTIFICATION
										WHERE PQ.ID_QUANTIFICATION=$selQuantifID AND ID_QUANTIF_PARAMETER=? GROUP BY PQ.ID_PROT_QUANTIF");
		my $pvalueCode=($labelingInfo{'FDR_CONTROL'}[0] eq 'TRUE')? 'PVAL_ADJ' : 'PVAL';
		my @extraParams=($action eq 'export')? ('CONF_LOW','CONF_UP') : ();
		foreach my $paramCode ('RATIO','SD_GEO',$pvalueCode,@extraParams) {
			$sthProtQ->execute($quantifParamInfo{$paramCode}[0]);
			while (my ($protID,$modResStrg,$ratioPos,$qValue)=$sthProtQ->fetchrow_array) {
				next unless $dispRatios{$ratioPos};
				next if ($restrictListID && !$restrictProteins{$protID});
				@{$proteinInfo{$protID}}=(); # actual number of quantified proteins (only list filtering)
				if ($modResStrg) { # quantif of modification
					$protID.='-'.&promsQuantif::decodeModificationSite($modResStrg);
				}
				#$quantifValues{$protID}{$paramCode}{$ratioPos}=($paramCode eq 'SD_GEO')? $qValue*100 : $qValue;
				$quantifValues{$protID}{$paramCode}{$ratioPos}=$qValue;
				$minPvalue=$qValue if ($paramCode eq $pvalueCode && $qValue > 0 && $qValue < $minPvalue);
				$usedIsoforms{$protID}=1;
			}
		}
		$minPvalue/=10;

		foreach my $paramCode ('DIST_PEP_USED','NUM_PEP_USED','NUM_PEP_TOTAL') { # NO SUPER RATIO EXPECTED FOR INTERNAL QUANTIF!!!
			$sthProtQ->execute($quantifParamInfo{$paramCode}[0]);
			while (my ($protID,$modResStrg,$ratioPos,$qValue)=$sthProtQ->fetchrow_array) {
				next if ($ratioPos && !$dispRatios{$ratioPos}); # no ratioPos for NUM_PEP_TOTAL
				next unless $proteinInfo{$protID}; # NUM_PEP_TOTAL has no ratioPos & could bring additional protID
				#next if ($restrictListID && !$restrictProteins{$protID});
				if ($modResStrg) { # quantif of modification
					$protID.='-'.&promsQuantif::decodeModificationSite($modResStrg);
					next unless $usedIsoforms{$protID};
				}
				if ($ratioPos) { # defined for +/-inf ratios
					$quantifValues{$protID}{$paramCode}{$ratioPos}=$qValue if $dispRatios{$ratioPos};
				}
				else { # undef for quantified proteins
					foreach my $r (keys %dispRatios) { # extend to all displayed ratios
						$quantifValues{$protID}{$paramCode}{$r}=$qValue;
					}
				}
			}
		}
		$sthProtQ->finish;

		###>Fetching protein info
		#my %proteinInfo;
		#my $extraFieldStrg=($view eq 'graph' && $action ne 'export')? ',PROT_LENGTH' : ',MW,PROT_DES,ORGANISM';
		#my $sthProtI=$dbh->prepare("SELECT P.ID_PROTEIN,ALIAS$extraFieldStrg FROM ANALYSIS_PROTEIN A,PROTEIN P WHERE A.ID_PROTEIN=P.ID_PROTEIN AND VISIBILITY>=1 AND ID_ANALYSIS=$analysisID");
		my $sthProtI=$dbh->prepare("SELECT DISTINCT P.ID_PROTEIN,ALIAS,MW,PROT_DES,ORGANISM,PROT_LENGTH FROM PROTEIN_QUANTIFICATION PQ,PROTEIN P WHERE PQ.ID_PROTEIN=P.ID_PROTEIN AND ID_QUANTIFICATION=$selQuantifID");
		my ($GNidentID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
		my $sthGN=$dbh->prepare("SELECT MI.VALUE FROM PROTEIN P,MASTERPROT_IDENTIFIER MI WHERE P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND MI.ID_IDENTIFIER=$GNidentID AND P.ID_PROTEIN=? ORDER BY MI.RANK");
		$sthProtI->execute;
		while (my ($protID,@info)=$sthProtI->fetchrow_array) {
			next if ($restrictListID && !$restrictProteins{$protID});
			next unless $quantifValues{$protID}; # only quantified proteines
			my @geneList;
			if ($view eq 'list' || $action eq 'export') {
				$info[1]=sprintf "%.1f",$info[1]/1000;
				$sthGN->execute($protID);
				while (my ($gene)=$sthGN->fetchrow_array) {push @geneList,$gene;}
			}
			$info[0]=~s/[,;']/\./g; # Clean identifier: MaxQuant contaminant can be crappy!!!
			@{$proteinInfo{$protID}}=(@info,\@geneList,$analysisID); # to matched multi-ana quantif
			#$proteinInfo{$protID}[1]=sprintf "%.1f",$proteinInfo{$protID}[1]/1000 if ($view eq 'list' || $action eq 'export'); # MW
		}
		$sthProtI->finish;
		$sthGN->finish;

		$dbh->disconnect;

		if ($action eq 'export') {
			&exportProteinList($labelType,$ratioType,\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,$pvalueCode);
			$workbook->close();
			exit;
		}

		print qq |<SCRIPT type="text/javascript">document.getElementById('waitDiv').style.display='none'</SCRIPT>|;
		#my $numRatios=scalar @{$labelingInfo{RATIOS}};
		#my $numPepCode=($dispPepType eq 'all')? 'NUM_PEP_USED' : 'DIST_PEP_USED';

		########################
		####>Graphical view<####
		########################
		if ($view eq 'graph') {
			&displayVolcanoPlot(\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,$pvalueCode);
		}

		###################
		####>List view<####
		###################
		else {
			print qq
|<SCRIPT type="text/javascript">
window.onload=function() {
	ajaxUpdateRestrict();
}
</SCRIPT>
|;
			&printProteinRatioList($labelType,$ratioType,\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,$pvalueCode); #$labelingInfo{FDR_ALPHA}[0]/100,$dispStdDev
		}
	}
	&endHTML;
}


#############################################<<<SUBROUTINES>>>###########################################

sub endHTML {
	print qq
|</CENTER>
<DIV id="displayDIV" class="popup"> <!--filter:alpha(opacity=80);opacity:0.8;-->
	<DIV id="infoDIV"></DIV>
</DIV>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT LANGUAGE="javascript">
setPopup();
</SCRIPT>
</BODY>
</HTML>
|;
}


sub generateRatioCodes {
	my ($labelType,$ratioType,$numStates,$isPepIntensity,$algoVersion)=@_;
	my %ratioCode;
	if ($ratioType=~/S\w+Ratio/) {
		if ($isPepIntensity) { # $labelType eq 'FREE'
			my $ratioPos=0;
			foreach my $refStatePos (1..$numStates-1) {
				foreach my $testStatePos ($refStatePos+1..$numStates) {
					next if $testStatePos==$refStatePos;
					$ratioPos++;
					$ratioCode{$ratioPos}=($algoVersion==3)? "State$testStatePos.State$refStatePos" : "State$testStatePos-State$refStatePos";
					$ratioCode{-$ratioPos}='-'.$ratioCode{$ratioPos};
				}
			}
		}
		elsif ($algoVersion==3) {
			my $numPrimRatios=$numStates-1;
			foreach my $ratioPos (1..$numPrimRatios) {
				$ratioCode{$ratioPos}='State'.($ratioPos+1).'.State1';
				$ratioCode{-$ratioPos}='-'.$ratioCode{$ratioPos};
			}
			if ($ratioType eq 'SuperRatio') {
				my $ratioPos=$numPrimRatios;
				foreach my $x (1..$numPrimRatios-1) { # ratios of ratios
					foreach my $y ($x+1..$numPrimRatios){
						next if $y==$x;
						$ratioPos++;
						$ratioCode{$ratioPos}="State$y.State1.State$x.State1";
						$ratioCode{-$ratioPos}='-'.$ratioCode{$ratioPos};
					}
				}
			}
		}
		else {
			if ($ratioType eq 'SuperRatio') {
				my $numPrimRatios=$numStates-1;
				foreach my $ratioPos (1..$numPrimRatios) {
					$ratioCode{$ratioPos}=chr(64+$ratioPos);
					$ratioCode{-$ratioPos}='-'.$ratioCode{$ratioPos};
				}
				my $ratioPos=$numPrimRatios;
				foreach my $x (1..$numPrimRatios-1) { # ratios of ratios
					foreach my $y ($x+1..$numPrimRatios){
						next if $y==$x;
						$ratioPos++;
						$ratioCode{$ratioPos}=chr(64+$y).'-'.chr(64+$x);
						$ratioCode{-$ratioPos}='-'.$ratioCode{$ratioPos};
					}
				}
			}
			else {
				my $numPrimRatios=$numStates*($numStates-1)/2;
				foreach my $ratioPos (1..$numPrimRatios) {
					$ratioCode{$ratioPos}=chr(64+$ratioPos);
					$ratioCode{-$ratioPos}='-'.$ratioCode{$ratioPos};
				}
			}
		}
	}
	else {
		foreach my $ratioPos (1..$numStates*($numStates-1)/2) {$ratioCode{$ratioPos}=$ratioPos; $ratioCode{-$ratioPos}=-$ratioPos;}
	}
	return %ratioCode;
}

sub printQuantificationSummary { # GLOBALS: $selQuantifID, $view, $projectID, $designID, $ratioType, $lightColor, $darkColor, $quantifType, $quantifSoftware
	my ($dbh,$anaID,$refInfo,$refGoAnalyses)=@_;
	#my ($lightColor,$darkColor)=&promsConfig::getRowColors;
	my $resultDir="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/results";
	my $graphDirHtml="$promsPath{quanti_html}/project_$projectID/quanti_$selQuantifID/results/graph";
	my ($newLine,$startB,$endB)=($action eq 'export')? ("\n",'','') : ('<BR>','<B>','</B>');

	my $ratioTypeStrg=($ratioType eq 'Ratio')? 'Simple ratios '.$startB.'(old version)'.$endB : ($ratioType eq 'SuperRatio')? 'Super ratios' : ($ratioType eq 'None')? 'No ratio' : 'Simple ratios';
	#my $quantifSoftware=($refInfo->{'SOFTWARE'})? $refInfo->{'SOFTWARE'}[0] : 'myProMS';
	#$quantifSoftware='MaxQuant' if $quantifSoftware eq 'MQ';
	my $softVersionStrg=($algoVersion)? ' v'.$algoVersion : ($refInfo->{'SOFTWARE'} && $refInfo->{'SOFTWARE'}[1])? ' v'.$refInfo->{'SOFTWARE'}[1] : ''; # any Software

	my $topNString='';
	if ($labelType eq 'FREE' && !$quantifModID && $quantifSoftware !~ /SWATH|DIA|MaxQuant/) {
		$topNString=($refInfo->{'MATCH_PEP'} && $refInfo->{'MATCH_PEP'}[0])? ' - '.$startB.'Matching'.$endB : ' - '.$startB.'Any'.$endB;
		$topNString.=($refInfo->{'TOP_N'})? ' top'.$startB.$refInfo->{'TOP_N'}[0].$endB.' peptide' : ($quantifSoftware eq 'MaxQuant')? '' : ' top'.$startB.'N'.$endB.' peptide';
		$topNString.='s' if ($topNString && (!$refInfo->{'TOP_N'} || $refInfo->{'TOP_N'}[0]>1));
	}
	my $swathStrg=($quantifSoftware=~/SWATH/)? ' (SWATH)' : ($quantifSoftware=~/DIA/) ? ' (DIA)' : '';

	#<Conditions & replicates
	my (%stateInfo,%replicateName,%numTechReplicates,%bioRepLabel);
	#my $existReplicates=0;
	#my $supRatioReplicates=0;
	my $numStates=0;
	###> Get parent quantifType to distinguish from label/label-free quantification
	#my ($parentQuantifType)=$dbh->selectrow_array("SELECT CODE FROM QUANTIFICATION,PARENT_QUANTIFICATION,QUANTIFICATION_METHOD WHERE QUANTIFICATION.ID_QUANTIFICATION=PARENT_QUANTIFICATION.ID_PARENT_QUANTIFICATION AND QUANTIFICATION.ID_QUANTIFICATION_METHOD=QUANTIFICATION_METHOD.ID_QUANTIFICATION_METHOD AND PARENT_QUANTIFICATION.ID_QUANTIFICATION=$selQuantifID LIMIT 1");
	#if ($parentQuantifType && $parentQuantifType eq 'XIC'){#} # Be careful with Pep-Ratio from label vs. label-free experiments... STATES are not formated the same!!!
	if ($designID) {
		my $sthAnaName=$dbh->prepare("SELECT NAME FROM ANALYSIS WHERE ID_ANALYSIS=?");
		my $sthExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
		my $sthPepQuantif=$dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
		my $sthObsBioSamp=$dbh->prepare("SELECT ID_BIOSAMPLE FROM OBSERVATION WHERE ID_OBSERVATION=?");
		my $sthObsMsSamp=$dbh->prepare("SELECT ID_SAMPLE FROM ANALYSIS A WHERE ID_ANALYSIS=?");
		my $sthBsName=$dbh->prepare("SELECT NAME FROM BIOSAMPLE WHERE ID_BIOSAMPLE=?");
		my $sthMsName=$dbh->prepare("SELECT NAME FROM SAMPLE WHERE ID_SAMPLE=?");
		my %anaName;
		foreach my $stateData (@{$refInfo->{STATES}}) {
			$stateData=~s/#//g; # removes id tags
			$numStates++;
			my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$stateData);
			#if ($ratioType eq 'SuperRatio') { #=~/S\w+Ratio/
			#	if ($numStates==1) {$supRatioReplicates=$numBioRep;}
			#	elsif ($numBioRep==1) {$supRatioReplicates--;}
			#}
			#$expCondID=~ s/^.//; # To remove the first character (which is a '#' to show that it is an ID)
			$sthExpCondName->execute($expCondID);
			($stateInfo{$numStates}{'NAME'})=$sthExpCondName->fetchrow_array;
			$condToState{$expCondID}=$numStates;
			$stateInfo{$numStates}{'POPUP'}='';
			#my (%replicateInfo,%replicateStrg);
			my (%channelName,%bioRepHeader);
			my ($currBioRep,$allTechRepCount,$numPools,$numObs)=(0,0,0,0);
			foreach my $bioReplicate (split(/\./,$quantiObsIDs)) { # bio rep separator
				$currBioRep++;
				my (%fracAnaNames,%numPoolObs,%bioRepNames);
				my $numTechRep=0;
				my $numObsInBioRep=0;
				foreach my $techReplicate (split(/&/,$bioReplicate)) { # tech rep separator
					$numTechRep++;
					#my @fracAnaNames;
					$numPoolObs{$numTechRep}=0;
					foreach my $pooledObs (split(/\+/,$techReplicate)) {
						$numPoolObs{$numTechRep}++;
						my ($obsID,$parQuantifID,$anaID,$targetPos);
						my @obsData=split(/:/,$pooledObs);
						if (scalar @obsData==4) {($obsID,$parQuantifID,$anaID,$targetPos)=@obsData;} # design post 07/2013 -> obsID,parQuantifID,anaID,targetPos
						else { # older design -> parQuantifID,anaID
							($parQuantifID,$anaID)=@obsData;
							$obsID=0;
							$targetPos=0;
						}
						unless ($anaName{$anaID}) {
							$sthAnaName->execute($anaID);
							($anaName{$anaID})=$sthAnaName->fetchrow_array;
						}
						if ($obsID) {
							$sthObsBioSamp->execute($obsID);
							my ($bsID)=$sthObsBioSamp->fetchrow_array;
							$bioRepNames{BIOSAMPLE}{$bsID}++ if $bsID;
						}
						$sthObsMsSamp->execute($anaID);
						my ($sampID)=$sthObsMsSamp->fetchrow_array;
						$bioRepNames{MSSAMPLE}{$sampID}++;
						my $targetName='';
						if ($targetPos) {
							unless ($channelName{$parQuantifID}) {
								$sthPepQuantif->execute($parQuantifID);
								my ($qAnnot)=$sthPepQuantif->fetchrow_array;
								my @quantifInfo=split(/::/,$qAnnot);
								foreach my $channelInfo (@quantifInfo[1..$#quantifInfo]) {
									my ($chanPos,$chanName)=split(/;/,$channelInfo); # + other unsued elements
									$channelName{$parQuantifID}{$chanPos}=$chanName;
								}
							}
							$targetName=' > '.$channelName{$parQuantifID}{$targetPos};
						}
						push @{$fracAnaNames{$numTechRep}},$anaName{$anaID}.$targetName;
					}
					$numObsInBioRep+=$numPoolObs{$numTechRep};
				}
				$allTechRepCount+=$numTechRep;
				$replicateName{$numStates}{$currBioRep}=($action eq 'export')? $stateInfo{$numStates}{'NAME'} : $startB.$stateInfo{$numStates}{'NAME'}.$endB;
				$numTechReplicates{"$numStates:$currBioRep"}=$numTechRep;

				##<Guess most suitable name for bioRep
				if ($bioRepNames{BIOSAMPLE} && scalar keys %{$bioRepNames{BIOSAMPLE}}==1) { # Bio-sample
					my $bsID=(keys %{$bioRepNames{BIOSAMPLE}})[0];
					#if ($bioRepNames{BIOSAMPLE}{$bsID}==$numObsInBioRep) {
						$sthBsName->execute($bsID);
						($bioRepLabel{"$numStates:$currBioRep"})=$sthBsName->fetchrow_array;
					#}
				}
				if (!$bioRepLabel{"$numStates:$currBioRep"} && scalar keys %{$bioRepNames{MSSAMPLE}}==1) { # MS Sample
					my $sampID=(keys %{$bioRepNames{MSSAMPLE}})[0];
					#if ($bioRepNames{MSSAMPLE}{$sampID}==$numObsInBioRep) {
						$sthMsName->execute($sampID);
						($bioRepLabel{"$numStates:$currBioRep"})=$sthMsName->fetchrow_array;
					#}
				}
				if (!$bioRepLabel{"$numStates:$currBioRep"}) {
					$bioRepLabel{"$numStates:$currBioRep"}='Replicate #'.$currBioRep;
				}

				if ($numObsInBioRep > 1) { # 2 or more fractions and/or tech rep
					if ($numBioRep > 1) {
						$replicateName{$numStates}{$currBioRep}.="$startB:$endB ".$bioRepLabel{"$numStates:$currBioRep"};
						if ($action eq 'export') {$stateInfo{$numStates}{'POPUP'}.=$startB.$bioRepLabel{"$numStates:$currBioRep"}.$endB.$newLine;}
						else {@{$bioRepHeader{$currBioRep}}=($bioRepLabel{"$numStates:$currBioRep"},$numTechRep);}
					}
					foreach my $techRep (1..$numTechRep) {
						$numPools++;
						if ($numTechRep > 1) {
							$stateInfo{$numStates}{'POPUP'}.=($action eq 'export')? "•Tech. rep. #$techRep" : "<TD valign=top nowrap>&bull;<U>Tech. rep. #$techRep";
							$stateInfo{$numStates}{'POPUP'}.=" (Pool)" if $numPoolObs{$techRep} > 1;
							$stateInfo{$numStates}{'POPUP'}.=($action eq 'export')? ":$newLine" : "</U>:$newLine";
						}
						elsif ($numPoolObs{$techRep} > 1) {
							$stateInfo{$numStates}{'POPUP'}.=($action eq 'export')? "•Pool:$newLine" : "<TD valign=top nowrap>&bull;<U>Pool:</U>$newLine";
						}
						if ($numTechRep > 1 || $numPoolObs{$techRep} > 1) {
							$stateInfo{$numStates}{'POPUP'}.=($action eq 'export')?
															"  +".join("$newLine  +",@{$fracAnaNames{$techRep}}).$newLine
															: "&nbsp;&nbsp;+".join("$newLine&nbsp;&nbsp;+",@{$fracAnaNames{$techRep}})."<BR>";
						}
						else {
							$replicateName{$numStates}{$currBioRep}.="$startB:$endB $fracAnaNames{$techRep}[0]";
							$stateInfo{$numStates}{'POPUP'}.=($action eq 'export')? "•$fracAnaNames{$techRep}[0]$newLine" : "<TD valign=top nowrap>&bull;$fracAnaNames{$techRep}[0]<BR>";
						}
					}
				}
				else { # only 1 obs in replicate
					$replicateName{$numStates}{$currBioRep}.="$startB:$endB $fracAnaNames{1}[0]";
					$stateInfo{$numStates}{'POPUP'}.=($action eq 'export')? "•$fracAnaNames{1}[0]$newLine" : "<TD valign=top nowrap>&bull;$fracAnaNames{1}[0]<BR>";
				}
				$numObs+=$numObsInBioRep;
				#$stateString.=" $replicateName{$numStates}{$currBioRep},";
			}
			#chop($stateString);
			#$stateInfo{$numStates}{'NAME_LONG'}=$stateString;
			$stateInfo{$numStates}{'NAME_LONG'}=($action eq 'export')? ': ' : ':&nbsp;';
			if ($numBioRep>1) {
				#$existReplicates=1;
				$stateInfo{$numStates}{'REPLIC'}=($action eq 'export' || $ratioType=~/S\w+Ratio/ || $quantifSoftware eq 'MaxQuant')? '' : "&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Correlation\" onclick=\"displayReplicCorrelation(this,$numStates,'$stateInfo{$numStates}{NAME}')\">"; # not true replicates if SuperRatio
				$stateInfo{$numStates}{'NAME_LONG'}.=($ratioType eq 'SuperRatio' && $numStates==1 && $numBioRep==scalar @{$refInfo->{STATES}})? "$numBioRep sets / " : "$numBioRep bio. rep. / ";
			}
			else {$stateInfo{$numStates}{'REPLIC'}='';}
			if ($allTechRepCount > $numBioRep) {
				$stateInfo{$numStates}{'NAME_LONG'}.="$allTechRepCount tech. rep. / ";
			}
			elsif ($numPools) {
				$stateInfo{$numStates}{'NAME_LONG'}.=($numPools > 1)? "$numPools pools / " : "$numPools pool / ";
			}
			$stateInfo{$numStates}{'NAME_LONG'}.="$numObs observation";
			$stateInfo{$numStates}{'NAME_LONG'}.='s' if $numObs > 1;

			if (scalar keys %bioRepHeader) {
				my $popupStrg='<TABLE border=1 cellspacing=0><TR>';
				foreach my $currBioRep (sort{$a<=>$b} keys %bioRepHeader) {
					$popupStrg.="<TH colspan=$bioRepHeader{$currBioRep}[1]>$bioRepHeader{$currBioRep}[0]</TH>";
				}
				$popupStrg.='</TR><TR>';
				$stateInfo{$numStates}{'POPUP'}=$popupStrg.'</TR><TR>'.$stateInfo{$numStates}{'POPUP'}.'</TR></TABLE>';
			}
#print "#$numStates: $numBioRep, $supRatioReplicates<BR>\n";
		}
		$sthAnaName->finish;
		$sthExpCondName->finish;
		$sthPepQuantif->finish;
		$sthObsBioSamp->finish;
		$sthObsMsSamp->finish;
		$sthBsName->finish;
		$sthMsName->finish;
	}
	else { # no design
		foreach my $stateData (@{$refInfo->{STATES}}) {
			$numStates++;
			(my $numBioRep,$stateInfo{$numStates}{'NAME'},my $repPosStrg)=split(',',$stateData);
			my $currBioRep=0;
			foreach my $replicName (split(/\./,$stateInfo{$numStates}{'NAME'})) {
				$currBioRep++;
				$replicateName{$numStates}{$currBioRep}=$replicName;
				$numTechReplicates{"$numStates:$currBioRep"}=1;
			}
			$stateInfo{$numStates}{'NAME'}=~s/\./\+/g;
			if ($numBioRep>1) {
				$stateInfo{$numStates}{'REPLIC'}=($action eq 'export')? " ($numBioRep replicates)" : "&nbsp;($numBioRep replicates)&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Correlation\" onclick=\"displayReplicCorrelation(this,$numStates,'$stateInfo{$numStates}{NAME}')\">";
				#$existReplicates=1;
			}
			else {$stateInfo{$numStates}{'REPLIC'}='';}
		}
	}

	#<Setting default displayed ratios (0 for MaxQuant no-ratio)
	if ($ratioType ne 'None') {
		my $maxDefaultRatios=$numStates-1; # all ratios with 1st reference
		unless (scalar keys %dispRatios) {
			foreach my $pos (1..$maxDefaultRatios) {$dispRatios{$pos}=1;}
		}
	}

	#<Peptides (retro-compatibility in case no miss-cut param)
	my ($pepRangeStrg,$pepMissCutStrg,$pepPtmStrg,$pepChargeStrg,$pepSourceStrg);
	if ($refInfo->{PEPTIDES}) {
		if ($quantifSoftware eq 'MaxQuant') {
			$pepRangeStrg=($refInfo->{PEPTIDES}[0] eq 'unique')? 'Unique peptides only' : ($refInfo->{PEPTIDES}[0] eq 'razor')? 'Razor+unique peptides allowed' : 'All peptides allowed';
		}
		else {
			$pepRangeStrg=($refInfo->{PEPTIDES}[0] eq 'unique')? 'Proteotypic peptides only' : ($refInfo->{PEPTIDES}[0] eq 'unique_shared')? 'Shared peptides allowed if exist proteotypic' : 'Non-proteotypic peptides allowed';
		}
		$pepMissCutStrg='Missed cleavage allowed';
		my $pepIdx=0;
		if (scalar @{$refInfo->{PEPTIDES}}==5) { # Missed cleavage option exists
			$pepMissCutStrg='Missed cleavage not allowed' unless $refInfo->{PEPTIDES}[1]==1;
			$pepIdx++;
		}
		$pepIdx++;
		my ($ptmAllowed,@selectedPTMs)=split(/[:,]/,$refInfo->{PEPTIDES}[$pepIdx]);
		if (abs($ptmAllowed)<=1) { # old options
			$pepPtmStrg=($ptmAllowed==1)? 'All modifications allowed' : ($ptmAllowed==-1)? 'Modified and unmodified matching peptides not allowed' : 'No modifications allowed';
		}
		else { # new custom selection options (2 or -2)
			my @ptmNames;
			my $sthVM=$dbh->prepare('SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DES FROM MODIFICATION WHERE ID_MODIFICATION=?');
			foreach my $modID (@selectedPTMs) {
				$sthVM->execute($modID);
				my ($psiName,$interName,$synName,$des)=$sthVM->fetchrow_array;
				my $ptmName=$psiName || $interName || $des;
				unless ($ptmName) {
					$synName=~s/^##//; $synName=~s/##$//; $synName=~s/##/, /;
					$ptmName=$synName;
				}
				push @ptmNames,$ptmName;
			}
			$sthVM->finish;
			if ($action eq 'export') {$pepPtmStrg='Modifications allowed: '.join(', ',sort{lc($a) cmp lc($b)} @ptmNames);}
			else {
				$pepPtmStrg="<A href=\"javascript:void(null)\" onmouseover=\"popup('<B>Modifications allowed:</B><BR>&nbsp;&bull;".join('<BR>&nbsp;&bull;',sort{lc($a) cmp lc($b)} @ptmNames)."')\" onmouseout=\"popout()\">Some modifications allowed<SUP>*</SUP></A>";
			}
			$pepPtmStrg.=' (Other modified and unmodified forms not allowed)' if $ptmAllowed==-2;
		}
		#if ($quantifType !~ /(PROT_RATIO|TNPQ|XIC|SIN|EMPAI)/) { #} Modified in 29/10/12
		#if ($quantifType !~ /XIC|SIN|EMPAI/) {
			$pepIdx++;
			$pepChargeStrg=($refInfo->{PEPTIDES}[$pepIdx] eq 'best')? 'Best charge state' : 'All charge states';
			$pepIdx++;
			$pepSourceStrg=($refInfo->{PEPTIDES}[$pepIdx] && $refInfo->{PEPTIDES}[$pepIdx] eq 'best')? 'Best source' : 'All sources';
		#}
	}
	my %convPtmParam=('ambiguous'=>'delocalized','exclude'=>'excluded','valid'=>'not delocalized');
	my $ptmQuantifStrg=(!$refInfo->{PTM_POS})? '' : ($refInfo->{PTM_POS}[0]=~/PRS:(.+)/)? "<B>$modifName</B>-site positions are <B>confirmed</B> if PRS probability &ge; <B>$1%</B>, others are <B>$convPtmParam{$refInfo->{PTM_POS}[1]}</B>" : "<B>$modifName</B>-sites are $convPtmParam{$refInfo->{PTM_POS}[0]}";

	#<Bias correction
	my %fdrMethods=('FDR-BH'=>'Benjamini-Hochberg',
					'FDR-ABH'=>'Benjamini-Hochberg (Adaptative)',
					'FWER-BONF'=>'Bonferroni',
					'Qvalue'=>'Qvalue (Storey et al.)'
				   );
	my $biasCorrFlag=(!$refInfo->{'BIAS_CORRECTION'} || $refInfo->{'BIAS_CORRECTION'}[0] eq 'FALSE')? 'false' : 'true';
	my ($biasCorrectStrg,$biasCorrectStrg2)=('','');
	if ($action eq 'export') {
		if ($refInfo->{'BIAS_CORRECTION'}[0] eq 'FALSE') {$biasCorrectStrg='None';}
		else {
			if ($refInfo->{'BIAS_CORRECTION'}[1]) { # Prot list used for bias correction
				my ($listName)=$dbh->selectrow_array("SELECT NAME FROM CATEGORY WHERE ID_CATEGORY=$refInfo->{BIAS_CORRECTION}[1]");
				my $notStrg=($refInfo->{'BIAS_CORRECTION'}[2] && $refInfo->{'BIAS_CORRECTION'}[2] eq 'not')? ' not' : '';
				$biasCorrectStrg="Proteins$notStrg in list: $listName";
			}
			else {
				$biasCorrectStrg='Unknown'; # default
				if ($ratioType=~/S\w+Ratio/) {
					if ($refInfo->{'NORMALIZATION_METHOD'}[0]) {
						$biasCorrectStrg=($normalizationNames{$refInfo->{'NORMALIZATION_METHOD'}[0]})? $normalizationNames{$refInfo->{'NORMALIZATION_METHOD'}[0]} : $refInfo->{'NORMALIZATION_METHOD'}[0];
					}
					else {$biasCorrectStrg=$normalizationNames{'global.mad.normalization'};}
				}
				else {$biasCorrectStrg="Scale normalization";}
			}
			#$biasCorrectStrg.="</TD></TR>\n<TR><TD></TD><TD>";
			if ($refInfo->{'BIAS_COEFF'}) {
				my $coeffIdx=0;
				foreach my $statePos (sort{$a<=>$b} keys %replicateName) {
					foreach my $replicPos (sort{$a<=>$b} keys %{$replicateName{$statePos}}) {
						foreach my $numTechRep (1..$numTechReplicates{"$statePos:$replicPos"}) {
							$biasCorrectStrg2.=$newLine if $coeffIdx;
							$biasCorrectStrg2.=$replicateName{$statePos}{$replicPos};
							$biasCorrectStrg2.=' : Tech. Rep. #'.$numTechRep if $numTechReplicates{"$statePos:$replicPos"} > 1;
							$biasCorrectStrg2.=($refInfo->{'BIAS_COEFF'}[$coeffIdx])? ' / '.$refInfo->{'BIAS_COEFF'}[$coeffIdx] : ' / ?';
							$coeffIdx++;
						}
					}
				}
			}
		}
	}
	else {
		if ($refInfo->{'BIAS_CORRECTION'}) {
			my $biasCorrButtonStrg=($quantifStatus<=0 || $quantifSoftware eq 'MaxQuant')? '' : "<INPUT type=\"button\" class=\"font11\" value=\"Box plots\" onclick=\"displayNormalization(this,$biasCorrFlag)\"/>";
			if ($refInfo->{'BIAS_CORRECTION'}[0] eq 'FALSE') {$biasCorrectStrg="<TD valign=top>None</TD><TD>&nbsp;&nbsp;$biasCorrButtonStrg</TD></TR>\n";}
			else {
				#my $normalizationMethod=($ratioType=~/S\w+Ratio/ && $refInfo->{'NORMALIZATION_METHOD'}[0])? $normalizationNames{$refInfo->{'NORMALIZATION_METHOD'}[0]} : ($ratioType=~/S\w+Ratio/)? $normalizationNames{'global.mad.normalization'} : "Scale normalization";
				my $normalizationMethod='Unknown'; # default
				if ($ratioType=~/S\w+Ratio/) {
					if ($refInfo->{'NORMALIZATION_METHOD'}[0]) {
						$normalizationMethod=($normalizationNames{$refInfo->{'NORMALIZATION_METHOD'}[0]})? $normalizationNames{$refInfo->{'NORMALIZATION_METHOD'}[0]} : $refInfo->{'NORMALIZATION_METHOD'}[0];
					}
					else {$normalizationMethod=$normalizationNames{'global.mad.normalization'};}
				}
				else {$normalizationMethod="Scale normalization";}
				$biasCorrectStrg="<TD colspan=2 nowrap><B>$normalizationMethod</B>";

				if ($refInfo->{'BIAS_CORRECTION'}[1]) { # Prot list used/not used for bias correction
					$refInfo->{'BIAS_CORRECTION'}[1]=~s/#//;
					my ($listName)=$dbh->selectrow_array("SELECT CONCAT(CL.NAME,' > ',C.NAME) FROM CATEGORY C,CLASSIFICATION CL WHERE C.ID_CLASSIFICATION=CL.ID_CLASSIFICATION AND ID_CATEGORY=$refInfo->{BIAS_CORRECTION}[1]");
					my $notStrg=($refInfo->{'BIAS_CORRECTION'}[2] && $refInfo->{'BIAS_CORRECTION'}[2] eq 'not')? ' not' : '';
					$biasCorrectStrg.=" based on <B>set of proteins$notStrg</B> in '$listName'";
				}

				$biasCorrectStrg.="&nbsp;&nbsp;$biasCorrButtonStrg" if ($quantifStatus > 0 && $algoVersion < 3 && !$refInfo->{'BIAS_COEFF'}); # display only boxplot before normalization has occured

				$biasCorrectStrg.="</TD></TR>\n<TR><TD></TD><TD nowrap>";
				if ($quantifStatus > 0) {
					if ($refInfo->{'BIAS_COEFF'}) {
						my $valueStrg='';
						my $coeffIdx=0;
						foreach my $statePos (sort{$a<=>$b} keys %replicateName) {
							foreach my $replicPos (sort{$a<=>$b} keys %{$replicateName{$statePos}}) {
								foreach my $numTechRep (1..$numTechReplicates{"$statePos:$replicPos"}) {
									$valueStrg.='<BR>' if $coeffIdx;
									$valueStrg.='&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&bull;'.$replicateName{$statePos}{$replicPos};
									$valueStrg.=' : Tech. Rep. #'.$numTechRep if $numTechReplicates{"$statePos:$replicPos"} > 1;
									$valueStrg.=($refInfo->{'BIAS_COEFF'}[$coeffIdx])? ' / '.$refInfo->{'BIAS_COEFF'}[$coeffIdx] : ' / ?';
									$valueStrg.='&nbsp;';
									$coeffIdx++;
								}
							}
						}
						$valueStrg="<DIV style=\"max-height:100px;overflow-y:scroll;display:inline-block\">$valueStrg</DIV>" if $coeffIdx > 4;
						$biasCorrectStrg.="$valueStrg</TD><TD valign=\"middle\">&nbsp;&nbsp;&nbsp;$biasCorrButtonStrg</TD></TR>\n";
					}
					elsif ($algoVersion==3) { # fetch data from file
						open(BIAS,"$resultDir/allbias_coef.txt");
						#my @biasData=<BIAS>;
						my (%header2Idx,%biasType);
						my $lastStateIdx=-1;
						while (<BIAS>) {
							chomp;
							my ($rowName,@lineData)=split(/\s+/,$_);
							if ($.==1) {
								my $colIdx=-1;
								foreach my $colName (@lineData) {
									$colIdx++;
									$header2Idx{$colName}=$colIdx;
									$lastStateIdx++ if $colName=~/State/; # not yet in Experiment section
								}
							}
							@{$biasType{$rowName}}=@lineData; # includes header
						}
						close BIAS;
						my %biasData;
						foreach my $i (0..$lastStateIdx) {
							my ($exp,$statePos,$bioRepPos)=$biasType{'run'}[$i]=~/^([^_]+)_State(\d+).*_rep(\d+)/;
							foreach my $biasName ('coefAdd','coefMul') {
								my $biasValue=(!$biasType{$biasName})? 1 : ($biasName eq 'coefAdd')? (sprintf "%.3f",1/(2**($biasType{$biasName}[$i]+$biasType{$biasName}[$header2Idx{$exp}])))*1 : (sprintf "%.3f",($biasType{$biasName}[$i]*$biasType{$biasName}[$header2Idx{$exp}]))*1;
								push @{$biasData{"$statePos:$bioRepPos"}},$biasValue;
							}
						}
						my $valueStrg='';
						my $numAllReplics=0;
						foreach my $statePos (sort{$a<=>$b} keys %replicateName) {
							foreach my $replicPos (sort{$a<=>$b} keys %{$replicateName{$statePos}}) {
								$valueStrg.='<BR>' if $valueStrg;
								$numAllReplics++;
								if ($isPepIntensity) {
									$valueStrg.=$replicateName{$statePos}{$replicPos}.': 1/'.$biasData{"$statePos:$replicPos"}[0].' ~ '.$biasData{"$statePos:$replicPos"}[1];
								}
								else { # PEP_RATIO
									next if $statePos==1; # State1 not listed
									$valueStrg.=$replicateName{$statePos}{$replicPos}.' <B>/</B> '.$replicateName{1}{$replicPos}.': 1/'.$biasData{"$statePos:$replicPos"}[0].' ~ '.$biasData{"$statePos:$replicPos"}[1];
								}
							}
						}
						$valueStrg="<DIV style=\"max-height:100px;overflow-y:scroll;display:inline-block\">$valueStrg</DIV>" if $numAllReplics > 5;
						$biasCorrectStrg.="$valueStrg</TD><TD valign=\"middle\">&nbsp;&nbsp;&nbsp;$biasCorrButtonStrg</TD></TR>\n";
					}
				}
			}
		}
	}
	my ($fdrStrg,$pvalStrg,$alterStrg)=('','','');
	if ($ratioType eq 'Ratio') {
		($fdrStrg,$pvalStrg)=($refInfo->{'FDR_CONTROL'}[0] eq 'FALSE')? ('None','Not given') : ($refInfo->{'FDR_ALPHA'}[0].'% ('.$fdrMethods{$refInfo->{'FDR_METHOD'}[0]}.')' , $refInfo->{P_VALUE_OUTLIER}[0] );
		$alterStrg=(!$refInfo->{ALTER})? '' : ($refInfo->{ALTER}[0] eq 'two.sided')? 'Two sided' : ($refInfo->{ALTER}[0] eq 'less')? 'Lesser' : 'Greater';
	}

	####################################
	####<Exporting summary to Excel>####
	####################################
	if ($action eq 'export') {
		my $worksheet1=$workbook->add_worksheet('Summary');
		my $xlsRow=0;
		$worksheet1->set_column(0,0,50);
		$worksheet1->set_column(1,1,40);
		$worksheet1->set_column(2,2,30);
		$worksheet1->merge_range($xlsRow,0,$xlsRow,2,$title,$itemFormat{'title'});

		##<Source>##
		my @quantifInfo=&promsMod::getItemInfo($dbh,'QUANTIFICATION',$selQuantifID);
		$xlsRow++;
		$worksheet1->merge_range(++$xlsRow,0,$xlsRow,2,'Source',$itemFormat{'mergeColHeader'});
		#<Name
		$worksheet1->write_string(++$xlsRow,0,'Name :',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$quantifInfo[-1]{'NAME'},$itemFormat{'mergeColText'});
		#<Location
		$worksheet1->write_string(++$xlsRow,0,'Location in myProMS :',$itemFormat{'headerR'});
		my $locationStrg='';
		foreach my $it (@quantifInfo) {
			last if $it->{'ITEM'} eq 'QUANTIFICATION';
			$locationStrg.=' > ' if $locationStrg;
			$locationStrg.=$it->{'NAME'};
		}
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$locationStrg,$itemFormat{'mergeColText'});
		#<Design
		$worksheet1->write_string(++$xlsRow,0,'Design level :',$itemFormat{'headerR'});
		my $designStrg=($designID)? 'Experiment' : 'Analysis';
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$designStrg,$itemFormat{'mergeColText'});
		#<Peptide quantif
		$worksheet1->write_string(++$xlsRow,0,'Peptide XIC extraction engine :',$itemFormat{'headerR'});
		my ($isMCQ)=$dbh->selectrow_array("SELECT 1 FROM PARENT_QUANTIFICATION PQ,QUANTIFICATION P WHERE PQ.ID_QUANTIFICATION=$selQuantifID && PQ.ID_PARENT_QUANTIFICATION=P.ID_QUANTIFICATION AND P.FOCUS='peptide' AND P.QUANTIF_ANNOT LIKE '%::ALLCHARGESTATES=%' LIMIT 0,1");
		my $pepQuantStrg=($isMCQ)? 'MassChroQ' : 'Proteome Discoverer';
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$pepQuantStrg,$itemFormat{'mergeColText'});
		#<Labeling
		$worksheet1->write_string(++$xlsRow,0,'Labeling :',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$labelingName{$labelType}.$swathStrg,$itemFormat{'mergeColText'});

		##<Parameters>##
		$xlsRow++;
		$worksheet1->merge_range(++$xlsRow,0,$xlsRow,2,'Quantification parameters',$itemFormat{'mergeColHeader'});
		$worksheet1->write_string(++$xlsRow,0,'Quantification method :',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$quantifMethDesc.' ['.$ratioTypeStrg.$topNString.' - Software: '.$quantifSoftware.$softVersionStrg.']',$itemFormat{'mergeColText'});
		#<States
		foreach my $statePos (sort{$a<=>$b} keys %stateInfo) {
			$worksheet1->write_string(++$xlsRow,0,"State #$statePos :",$itemFormat{'headerR'});
			my $replicateNames=($stateInfo{$statePos}{'NAME_LONG'})? $stateInfo{$statePos}{'NAME_LONG'} : '';
			if ($stateInfo{$statePos}{'POPUP'}) {
				$worksheet1->write_string($xlsRow,1,"$stateInfo{$statePos}{NAME}$replicateNames",$itemFormat{'text'});
				$worksheet1->write_string($xlsRow,2,$stateInfo{$statePos}{'POPUP'},$itemFormat{'textWrap'});
			}
			else {$worksheet1->merge_range($xlsRow,1,$xlsRow,2,"$stateInfo{$statePos}{NAME}$replicateNames",$itemFormat{'mergeColText'});}
		}
		#<Peptide selection
		my ($pepString,$numPepLines)=($quantifType !~ /SIN|EMPAI/ && $quantifSoftware !~ /SWATH|DIA/ )? ("\n•$pepChargeStrg.\n•$pepSourceStrg.",5) : ('',3);
		if ($ptmQuantifStrg) {
			$numPepLines++;
			$ptmQuantifStrg="\n•$ptmQuantifStrg.";
		}
		$worksheet1->set_row(++$xlsRow,13*$numPepLines);
		$worksheet1->write_string($xlsRow,0,'Peptide selection :',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,"•$pepRangeStrg.\n•$pepMissCutStrg.\n•$pepPtmStrg.$pepString$ptmQuantifStrg",$itemFormat{'mergeColText'});
		#<Protein selection
		my $protSelectionStrg='All visible proteins'; # default
		if ($refInfo->{'PROTEINS'}) {
			$protSelectionStrg=($refInfo->{'PROTEINS'}[0] eq 'exclude')? 'Exclude proteins in List ' : 'Restrict to proteins in List ';
			(my $listID=$refInfo->{'PROTEINS'}[1])=~s/#//;
			my ($listStrg)=$dbh->selectrow_array("SELECT CONCAT(T.NAME,' > ',L.NAME) FROM CATEGORY L,CLASSIFICATION T WHERE L.ID_CLASSIFICATION=T.ID_CLASSIFICATION AND ID_CATEGORY=$listID");
			$protSelectionStrg.="$listStrg";
		}
		$worksheet1->write_string(++$xlsRow,0,'Protein selection :',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$protSelectionStrg.' ('.$protQuantifiedStrg.')',$itemFormat{'mergeColText'});
		#<Bias
		$worksheet1->write_string(++$xlsRow,0,'Bias correction :',$itemFormat{'headerR'});
		if ($biasCorrectStrg2) {
			$worksheet1->write_string($xlsRow,1,$biasCorrectStrg,$itemFormat{'text'});
			$worksheet1->write_string($xlsRow,2,$biasCorrectStrg2,$itemFormat{'textWrap'});
		}
		else {$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$biasCorrectStrg,$itemFormat{'mergeColText'});}
		#<Inf ratios
		$worksheet1->write_string(++$xlsRow,0,'Infinite ratios :',$itemFormat{'headerR'});
		my $infRatioString=($refInfo->{'MIN_INF_RATIOS'} && $refInfo->{'MIN_INF_RATIOS'}[0])? 'Avoided whenever possible' : 'Not avoided';
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$infRatioString,$itemFormat{'mergeColText'});

		#<Advanced settings
		if ($ratioType eq 'Ratio') {
			if ($refInfo->{'THRESHOLD_CV'}) { # replicates used TODO: Make sure THRESHOLD_CV is reported for design quantif
				$worksheet1->write_string(++$xlsRow,0,'Variation coefficient threshold between replicates :',$itemFormat{'headerR'});
				my $biasCorrStrg.=($refInfo->{'THRESHOLD_CV'}[0] eq 'FALSE')? 'Auto' : $refInfo->{'THRESHOLD_CV'}[1];
				$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$biasCorrStrg,$itemFormat{'mergeColText'});
			}
			$worksheet1->write_string(++$xlsRow,0,'FDR control :',$itemFormat{'headerR'});
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$fdrStrg,$itemFormat{'mergeColText'});
			$worksheet1->write_string(++$xlsRow,0,'p-value threshold for outlier detection :',$itemFormat{'headerR'});
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$pvalStrg,$itemFormat{'mergeColText'});
			if ($refInfo->{'ALTER'}) {
				$worksheet1->write_string(++$xlsRow,0,'Alternative hypothesis for comparison :',$itemFormat{'headerR'});
				$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$alterStrg,$itemFormat{'mergeColText'});
			}
			if ($refInfo->{'CONFIDENCE_LEVEL'}) {
				$worksheet1->write_string(++$xlsRow,0,'Confidence interval on protein abundance :',$itemFormat{'headerR'});
				$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$refInfo->{'CONFIDENCE_LEVEL'}[0],$itemFormat{'mergeColText'});
			}
		}

		##<Creation date & user
		$updateDate=&promsMod::formatDate($updateDate);
		$worksheet1->write_string(++$xlsRow,0,'Created :',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,"$updateDate by $updateUser",$itemFormat{'mergeColText'});

		##<Export settings>##
		$xlsRow++;
		$worksheet1->merge_range(++$xlsRow,0,$xlsRow,2,'Export settings',$itemFormat{'mergeColHeader'});
		$worksheet1->set_row(++$xlsRow,13 * scalar keys %dispRatios);
		$worksheet1->write_string($xlsRow,0,'Ratios exported :',$itemFormat{'headerR'});
		my $disRatiosStrg='';
		foreach my $rPos (sort{$a<=>$b} keys %dispRatios) {
			$disRatiosStrg.="\n" if $disRatiosStrg;
			my ($testStatePos,$refStatePos)=split(/\//,$refInfo->{RATIOS}[$rPos-1]);
			my $ratioTag='';
			if ($designID) {
				if ($testStatePos=~/%/) { # SuperRatio
					$ratioTag='°';
					$testStatePos=~s/%\d+//;
					$refStatePos=~s/%\d+//;
				}
				$testStatePos=$condToState{$testStatePos};
				$refStatePos=$condToState{$refStatePos};
			}
			$disRatiosStrg.="•$stateInfo{$testStatePos}{NAME}$ratioTag/$stateInfo{$refStatePos}{NAME}$ratioTag";
		}
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$disRatiosStrg,$itemFormat{'mergeColText'});
		my ($foldChgTypeStrg,$dispFcStrg)=($foldChangeType eq 'abs')? ('Absolute fold change ≥',"$dispFoldChange (∞=1000, 1/∞=0.001)") : ($foldChangeType eq 'up')? ('Fold change ≥',"$dispFoldChange (∞=1000)") : ('Fold change ≤',"1/$dispFoldChange (1/∞=0.001)");
		$worksheet1->write_string(++$xlsRow,0,$foldChgTypeStrg,$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,"$dispFcStrg",$itemFormat{'mergeColText'});
		$worksheet1->write_string(++$xlsRow,0,'p-value ≤',$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$dispPvalue,$itemFormat{'mergeColText'});
		if ($quantifType eq 'PROT_RATIO_PEP') {
			#$worksheet1->write_string(++$xlsRow,0,'Standard deviation ≤',$itemFormat{'headerR'});
			#$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$dispStdDev,$itemFormat{'mergeColText'});
			$worksheet1->write_string(++$xlsRow,0,'Coefficient of variation ≤',$itemFormat{'headerR'});
			my $cvStrg=($dispCV)? $dispCV.'%' : '*No filter*';
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$cvStrg,$itemFormat{'mergeColText'});
		}
		my $numPepStrg=($dispPepType eq 'distinct')? 'Distinct peptides used ≥' : 'Peptides used ≥';
		$worksheet1->write_string(++$xlsRow,0,$numPepStrg,$itemFormat{'headerR'});
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$dispNumPep,$itemFormat{'mergeColText'});
		$worksheet1->write_string(++$xlsRow,0,'Sort by :',$itemFormat{'headerR'});
		my $sortOption;
		if ($dispSort=~/ratio_(\d+)/) {
			my $rPos=$1;
			my ($testStatePos,$refStatePos)=split(/\//,$refInfo->{RATIOS}[$rPos-1]);
			my $ratioTag='';
			if ($designID) {
				if ($testStatePos=~/%/) { # SuperRatio
					$ratioTag=encode_utf8('°');
					$testStatePos=~s/%\d+//;
					$refStatePos=~s/%\d+//;
				}
				$testStatePos=$condToState{$testStatePos};
				$refStatePos=$condToState{$refStatePos};
			}
			$sortOption="Ratio: '$stateInfo{$testStatePos}{NAME}$ratioTag/$stateInfo{$refStatePos}{NAME}$ratioTag'";
		}
		if ($dispSort=~/state_(\d+)/) {
			$sortOption="$quantifParamInfo{$dispMeasure}[1]: $stateInfo{$1}{NAME}";
		}
		else {$sortOption=($dispSort eq 'identifier')? 'Identifiers' : ($dispSort eq 'peptide')? 'Peptides used' : ($dispSort eq 'mw')? 'Molecular weight' : '';}
		$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$sortOption,$itemFormat{'mergeColText'});
		if ($restrictListID) {
			my ($listStrg)=$dbh->selectrow_array("SELECT CONCAT(CL.NAME,' > ',C.NAME) FROM CATEGORY C,CLASSIFICATION CL WHERE C.ID_CLASSIFICATION=CL.ID_CLASSIFICATION AND ID_CATEGORY=$restrictListID");
			$worksheet1->write_string(++$xlsRow,0,'Restrict to proteins in List :',$itemFormat{'headerR'});
			$worksheet1->merge_range($xlsRow,1,$xlsRow,2,$listStrg,$itemFormat{'mergeColText'});
		}

		return %stateInfo;
	} # end of &printQuantificationSummary when action='export'


	####################################
	####<Displaying summary in HTML>####
	####################################
	my $correlImgFileName=($quantifSoftware eq 'myProMS')? 'Beforematrixcorrelation.jpeg' : 'AllProtein_Before_matrixcorrelation.jpeg'; # !!! Temp! Remove "AllProtein_" for SWATH!!!!
	#my $correlMatFileName=($quantifSoftware eq 'myProMS')? 'Beforematrixcorrelation.txt' : 'AllProtein_Before_matrixcorrelation.txt'; # !!! Temp! Remove "AllProtein_" for SWATH!!!!
	my $correlMatFileName=($algoVersion >= 2)? 'Beforematrixcorrelation.txt' : 'AllProtein_Before_matrixcorrelation.txt'; # !!! Temp! Remove "AllProtein_" for SWATH!!!!
	my $ajaxPepAction=($quantifSoftware=~/SWATH|DIA/)? 'ajaxPepSwath' : ($quantifSoftware eq 'MaxQuant')? 'ajaxPepMaxQuant' : 'ajaxPepRatio';
	my $dispMeasParamStrg=($quantifSoftware eq 'MaxQuant' && $ratioType eq 'None')? "&dispMeasure=$dispMeasure" : '';
	print qq
|<SCRIPT type="text/javascript">
// AJAX ------>
function getXMLHTTP() {
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

// <--- AJAX
function getDomOffset( Obj, Prop ) {
	var iVal = 0;
	while (Obj && Obj.tagName != 'BODY') {
		eval('iVal += Obj.' + Prop + ';');
		Obj = Obj.offsetParent;
	}
	return iVal;
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
function getElementPosition(e) {
	var left=0;
	var top=0;
	while (e.offsetParent != undefined && e.offsetParent != null) { //Adding parent item position
		left += e.offsetLeft + (e.clientLeft != null ? e.clientLeft : 0);
		top += e.offsetTop + (e.clientTop != null ? e.clientTop : 0);
		e = e.offsetParent;
	}
	return [top,left];
}
function manageTemplate(value){
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=manageTemplate&id_quantif=$selQuantifID&projectID=$projectID&TEMPLATE="+value,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			//document.getElementById('templateSpan').innerHTML='';
			document.getElementById('templateSpan').innerHTML=XHR.responseText;
			if (value == 'add') {
				alert("Succesfully saved as template.");
			}
			else {
				alert("Succesfully removed from templates.");
			}
		}
	}
	XHR.send(null);
}
|;
	if ($quantifSoftware ne 'MaxQuant') {
		print qq
|function displayNormalization(biasButton,hasCorrection) {
	var imageStrg="<TABLE cellpadding=0 cellspacing=0><TR>";
	if ('$quantifSoftware'.match('SWATH\|DIA')) {imageStrg+="<TH class=\\"title2\\">Fragment intensities";}
	else if (hasCorrection) {imageStrg+="<TH colspan=2 class=\\"title2\\">Bias correction";}
	else {imageStrg+="<TH class=\\"title2\\">Peptide intensities";}
	imageStrg+="&nbsp;&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></TH></TR>";
	if (hasCorrection) {imageStrg+="<TR><TH>Before</TH><TH>After</TH></TR>";}
	var pepBoxRepBefore,pepBoxRepAfter;
	if ('$ratioType'=='Ratio') {
		pepBoxRepBefore='pep_box_repl_before.png';
		pepBoxRepAfter='pep_box_repl_after.png';
	}
	else { // (Super/Simple)Ratio or SWATH
		pepBoxRepBefore='Beforeallpeptide.jpeg';
		pepBoxRepAfter='Afterallpeptide.jpeg';
	}
	imageStrg+="<TR><TD><IMG src='$graphDirHtml/"+pepBoxRepBefore+"' width=480 heigth=480></TD>";
	if (hasCorrection) {imageStrg+="<TD><IMG src='$graphDirHtml/"+pepBoxRepAfter+"' width=480 heigth=480></TD>";}
	imageStrg+="</TR>";
	var infoDiv=document.getElementById('infoDIV').innerHTML=imageStrg;
	var [top,left]=getElementPosition(biasButton);
	top+=30;
	left-=(hasCorrection)? 460 : 230;
//alert('window.innerWidth='+window.innerWidth+', document.body.clientWidth='+document.body.clientWidth+' ,document.documentElement.clientWidth='+document.documentElement.clientWidth);
	//var left=(document.body.clientWidth)? (document.body.clientWidth/2)-500 : (window.innerWidth)? (window.innerWidth/2)-500 : (document.documentElement.clientWidth)? (document.documentElement.clientWidth/2)-500 : ;
	if (left < 0) {left=0;}
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = left+'px';
	displayDiv.style.top = top+'px';
	displayDiv.style.display='block';
}
function displayReplicCorrelation(corrButton,statePos,replicName) {
	var [top,left]=getElementPosition(corrButton);
	top+=20;
	left-=250;
	if (left < 0) {left=0;}
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = left+'px';
	displayDiv.style.top = top+'px';
	displayDiv.style.display='block';
	var infoDiv=document.getElementById('infoDIV');
	var htmlStrg;
	if (statePos) {
		htmlStrg="<TABLE cellpadding=0 cellspacing=0><TR><TH colspan=2 class=\\"title2\\">Correlation between replicates for '"+replicName+"'&nbsp;&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></TH></TR>";
		htmlStrg+='<TR><TD><IMG src="$graphDirHtml/pep_corr_repl_cond'+statePos+'.png"></TD></TR></TABLE>';
		infoDiv.innerHTML=htmlStrg;
	}
	else { // S(uper/imple)Ratio
		htmlStrg="<TABLE cellpadding=0 cellspacing=0><TR><TH class=\\"title2\\" colspan=2>Correlation between all replicates&nbsp;&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></TH></TR>";
		htmlStrg+='<TR><TD></TD><TH class="darkBg">Before normalization</TH></TR>';
		htmlStrg+='<TR><TH class="darkBg">A<BR>f<BR>t<BR>e<BR>r<BR><BR>n<BR>o<BR>r<BR>m<BR>a<BR>l<BR>i<BR>z<BR>a<BR>t<BR>i<BR>o<BR>n</TH>';
		htmlStrg+='<TD><DIV id="globCorrDIV"><BR><BR><IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR></DIV></TD></TR></TABLE>';
		infoDiv.innerHTML=htmlStrg;
		//Creation of the XMLHTTPRequest object
		var XHR = getXMLHTTP();
		if (!XHR) {
			return false;
		}
		XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxDispCorrel&id_quantif=$selQuantifID&projectID=$projectID",true);
		XHR.onreadystatechange=function() {
			if (XHR.readyState==4 && XHR.responseText) {
				document.getElementById('globCorrDIV').innerHTML='';
				eval(XHR.responseText); // javascript code
			}
		}
		XHR.send(null);
	}
}
/*
function displayReplicCorrelation_old(corrButton,statePos,replicName) {
	var imageStrg;
	if (statePos) {
		imageStrg="<TABLE cellpadding=0 cellspacing=0><TR><TH colspan=2 class=\\"title2\\">Correlation between replicates for '"+replicName+"'&nbsp;&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></TH></TR>";
		imageStrg+='<TR><TD><IMG src="$graphDirHtml/pep_corr_repl_cond'+statePos+'.png"></TD></TR></TABLE>';
	}
	else { // S(uper/imple)Ratio
		imageStrg="<TABLE cellpadding=0 cellspacing=0><TR><TH class=\\"title2\\">Correlation between all replicates&nbsp;&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></TH></TR>";
		imageStrg+='<TR><TD><IMG src="$graphDirHtml/$correlImgFileName" width=500 heigth=500></TD></TR></TABLE>';
	}
	var infoDiv=document.getElementById('infoDIV').innerHTML=imageStrg;
	var [top,left]=getElementPosition(corrButton);
	top+=20;
	left-=250;
	if (left < 0) {left=0;}
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = left+'px';
	displayDiv.style.top = top+'px';
	displayDiv.style.display='block';
}
*/
function displayImageFromButton(buttonObj,title,image,imgWidth) {
	if (!imgWidth) {imgWidth=700;}
	var imageStrg="<TABLE cellpadding=0 cellspacing=0><TR><TH colspan=2 class=\\"title2\\">"+title+"&nbsp;&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></TH></TR>";
	imageStrg+="<TR><TD><IMG src='$graphDirHtml/"+image+"' width="+imgWidth+"></TD></TR></TABLE>";
	var infoDiv=document.getElementById('infoDIV').innerHTML=imageStrg;
	var [top,left]=getElementPosition(buttonObj);
	top+=20;
	left-=350;
	if (left < 0) {left=0;}
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = left+'px';
	displayDiv.style.top = top+'px';
	displayDiv.style.display='block';
}
function displayRatioStat(linkID,ratioCode,ratioName) {
	var ratioStrg='Ratio: '+ratioName;
	if (ratioCode.match(/^-/)) { // reverse ratio
		ratioCode=ratioCode.replace('-',''); // true ratioCode (replace w/o regExp -> 1st '-' only)
		var rName=ratioName.split('/');
		ratioName=rName[1]+'/'+rName[0]+' (Reversed)';
		ratioStrg='Reverse ratio: '+rName[1]+'/'+rName[0]
	}
	var ratioLink=document.getElementById(linkID);
	var imageStrg="<TABLE cellpadding=0 cellspacing=0><TR><TH colspan=2 class=\\"title2\\">"+ratioStrg+"'&nbsp;<INPUT type=\\"button\\" class=\\"font11\\" value=\\" Close \\" onclick=\\"document.getElementById('displayDIV').style.display='none'\\"></TH></TR>";

	if ('$ratioType'=='Ratio') {
		imageStrg+="<TR><TH>Boxplot</TH><TH>Distribution</TH></TR>";
		imageStrg+="<TR><TD><IMG src=\\"$graphDirHtml/prot_ratio_box_"+ratioCode+".png\\"></TD><TD><IMG src=\\"$graphDirHtml/prot_ratio_hist_"+ratioCode+".png\\"></TD></TR>";
	}
	else if ($algoVersion==3) {
		imageStrg+="<TR><TH colspan=2><BR>Fold Change distribution<BR><IMG src=\\"$graphDirHtml/distriblog2FC_"+ratioCode+".jpeg\\" width=600></TH></TR>";
		imageStrg+="<TR><TH colspan=2><BR>p-value distribution<BR><IMG src=\\"$graphDirHtml/distribPValue_"+ratioCode+".jpeg\\" width=600></TH></TR>";
	}
	else { // (Super/Simple)Ratio
		if ('$labelType'=='FREE') {
			imageStrg+="<TR><TH colspan=2>Distribution</TH></TR>";
			imageStrg+="<TR><TD colspan=2><IMG src=\\"$graphDirHtml/Graph_Box_Hist_ratioprot_"+ratioCode+".png\\"></TD></TR>";
		}
		else {
			if (ratioCode.match('-') \|\| '$refInfo->{BIAS_CORRECTION}[0]'=='FALSE') { // ratio of ratio or no normalization
				imageStrg+="<TR><TD colspan=2><IMG src=\\"$graphDirHtml/Graph_Box_Hist_ratioprot_"+ratioCode+".png\\"></TD></TR>";
			}
			else { // primary ratio with normalization
				//imageStrg+="<TR><TD><TABLE cellpadding=0 cellspacing=0><TR><TH>Before normalization</TH><TH>After normalization</TH></TR><TR><TD><IMG src=\\"$graphDirHtml/MAplotPep_"+ratioCode+"Before.png\\" width=240></TD><TD><IMG src=\\"$graphDirHtml/MAplotPep_"+ratioCode+"After.png\\" width=240></TD></TR></TABLE></TD><TD><IMG src=\\"$graphDirHtml/Graph_Box_Hist_ratioprot_"+ratioCode+".png\\"></TD></TR>";
				imageStrg+="<TR><TH>Before normalization<BR><IMG src=\\"$graphDirHtml/MAplotPep_"+ratioCode+"Before.png\\" width=320><BR>After normalization<BR><IMG src=\\"$graphDirHtml/MAplotPep_"+ratioCode+"After.png\\" width=320></TH><TD valign=\\"top\\"><IMG src=\\"$graphDirHtml/Graph_Box_Hist_ratioprot_"+ratioCode+".png\\"></TD></TR>";
			}
		}
	}
	imageStrg+="</TABLE>";
	var infoDiv=document.getElementById('infoDIV').innerHTML=imageStrg;
	var [top,left]=getElementPosition(ratioLink);
	top+=20;
	left-=500;
	if (left < 0) {left=0;}
	var displayDiv=document.getElementById('displayDIV');
	displayDiv.style.left = left+'px';
	displayDiv.style.top = top+'px';
	displayDiv.style.display='block';
}
|;
	}
	if ($action eq 'summary') {
		print qq
|</SCRIPT>
<BR>
|;
	}
	else {
		print qq
|function updateSortBy() {
	if ('$view'=='graph') {
		if (document.getElementById('listFromGraph') && document.getElementById('protListDIV').style.display != 'none') { // a protein list is displayed from graph
			ajaxListSelectedProteins(lastListFromGraph[0],lastListFromGraph[1]); // re-fetch list to update sort order
		}
	}
	else {document.displayForm.submit();} // view=list => resubmit form
}
function updatedListDisplay() {
	var foldChange=document.getElementById('foldChange').value;
	if (foldChange < 1) {
		alert('ERROR: Fold change must be >= 1!');
		return;
	}
/*
	var stdDev=document.getElementById('stdDev').value;
	if (stdDev < 1) {
		alert('ERROR: Standard deviation must be >= 1!');
		return;
	}
*/
	var pValue=document.getElementById('pValue').value;
	if (pValue < 0) {
		alert('ERROR: p-value must be >= 0!');
		return;
	}
	var sort=document.getElementById('sort').value;
	window.location="$promsPath{cgi}/showProtQuantification.cgi?CALL=$call&id_ana=$analysisID&id_quantif=$selQuantifID&ACT=$action&view=$view&foldChange="+foldChange+"&pValue="+pValue+"&sort="+sort; //+"&stdDev="+stdDev
}
var selGoID=0;
var goAspects={};
|;
		foreach my $goID (keys %{$refGoAnalyses}) {
			print "goAspects[$goID]='$refGoAnalyses->{$goID}[1]';\n";
		}
		print qq
|function updateGoAspects(goID) {
	var selAspect=document.getElementById('goAspect');
	selAspect.options.length=1;
	document.getElementById('goTermsDIV').innerHTML='';
	if (goID==0) {
		selAspect.style.visibility='hidden';
	}
	else {
		var i=1;
		if (goAspects[goID].match('P')) {selAspect.options[i]=new Option('Biological process',goID+',P'); i++;}
		if (goAspects[goID].match('C')) {selAspect.options[i]=new Option('Cellular component',goID+',C'); i++;}
		if (goAspects[goID].match('F')) {selAspect.options[i]=new Option('Molecular function',goID+',F');}
		selAspect.style.visibility='visible';
	}
}
function checkAllProteins(chkStatus) {
	var checkBoxList=document.protForm.chkProt;
	if (checkBoxList.length) {
		for (var i=0; i < checkBoxList.length; i++) {checkBoxList[i].checked=chkStatus;}
	}
	else {checkBoxList.checked=chkStatus;}
}
function exportProteins() {
	var defAction=document.displayForm.ACT.value;
	document.displayForm.ACT.value='export';
	document.displayForm.submit();
	document.displayForm.ACT.value=defAction;
}
function exportIsoforms() {
	//TODO: Check for selected proteins
	document.protForm.action='../sliva/createPhosphoSequence.cgi';
	document.protForm.submit();
	document.protForm.action='./showProtQuantification.cgi';
}
|;
		if ($view eq 'graph') {
			my $dispTargetPosStrg=($ratioType eq 'None')? join(',',sort{$a<=>$b} keys %dispStates) : join(',',sort{abs($a)<=>abs($b)} keys %dispRatios);
			#my $maxRatioAllowed=($ratioType eq 'Ratio')? scalar @{$refInfo->{RATIOS}} : (scalar @{$refInfo->{STATES}})-1;
			print qq
|//datasetIdx to targetPos
var dataset2TargetPos=[$dispTargetPosStrg];
// Color list for plot
var hColors=['#E18B6B','#95B9C7','#7E2217','#9A9A9A','#8AFB17','#FBB917','#F660AB','#000000','#4AA02C']; // reverse of volcanoPlot library
var hColorIdx=0;
var usedHighLights={},usedColorIdx=[];
for (let idx=0; idx<hColors.length; idx++) {usedColorIdx[idx]=0;}
function updateUsedColors(hlName,action) { // callback for chart highlighting
	if (action != 'delete') return;
	var idx=usedHighLights[hlName];
	delete usedHighLights[hlName];
	usedColorIdx[idx]--;
	if (usedColorIdx[idx]==0) {hColorIdx=idx;}
}
function proteinLink(dSetIdx,identifier,protID) {
	ajaxProteinData(null,protID,dataset2TargetPos[dSetIdx],'$ajaxPepAction');
}
// AJAX --->
var lastListFromGraph=[]; // needed to recall ajaxListSelectedProteins if "Sort by" is changed in graph view
function ajaxListSelectedProteins(selectedPoints,thresholds) {
	lastListFromGraph=[selectedPoints,thresholds];
	saveListGlobals.themesFetched=false; // used for themes management
	document.getElementById('displayDIV').style.display='none'; // in case popup is displayed
	// Adjusting protListDIV position to graph (graph can move due to summary option display)
	var graphDiv=document.getElementById('mainGraphDIV');
	var divX = getDomOffset(graphDiv,'offsetLeft') + graphDiv.offsetWidth + 2;
	var divY = getDomOffset(graphDiv,'offsetTop');
	var listDiv=document.getElementById('protListDIV');
	listDiv.style.left = divX + 'px';
	listDiv.style.top = divY + 'px';
	// Updating contents
	listDiv.innerHTML="<BR><BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<FONT class=\\"title3\\">Fetching data. Please wait...</FONT><BR>&nbsp;&nbsp;&nbsp;&nbsp;<IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR>";
	listDiv.style.display='block';
	//Parameters (extracting no-duplicates list of proteins)
	var noDupProteins=new Object();
	for (var gr in selectedPoints) {
		for (var i=0; i<selectedPoints[gr].length; i++) {
			noDupProteins[selectedPoints[gr][i]]=1;
		}
	}
	var paramStrg='ACT=ajaxListProt&id_ana=$anaID&id_quantif=$selQuantifID&view=$view&dispTargetPosAjax=$dispTargetPosStrg&pepType=$dispPepType&numPep=$dispNumPep&sort='+document.getElementById('sort').value;
	if (document.displayForm.dispMeasure) { // for MaxQuant intensities
		paramStrg+='&dispMeasure='+document.displayForm.dispMeasure.value+'&selProt=';
	}
	else { // Ratio
		paramStrg+='&foldChange=$dispFoldChange&pValue=$dispPvalue&coefVar=$dispCV&selProt='; //&stdDev=\$dispStdDev
	}
	var p1=true;
	for (var prot in noDupProteins) {
		if (!p1) {paramStrg+=',';}
		paramStrg+=prot;
		p1=false;
	}
//paramStrg+='&thresholds='+thresholds.join(','); // no longer used
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("POST","$promsPath{cgi}/showProtQuantification.cgi",true);
	//Send the proper header information along with the request
	XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
	XHR.setRequestHeader("Content-length", paramStrg.length);
	XHR.setRequestHeader("Connection", "close");
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			listDiv.innerHTML='<INPUT type="hidden" id="listFromGraph" value="1"/>'+XHR.responseText; // listFromGraph: flag
		}
	}
	XHR.send(paramStrg);
}
function ajaxListDecorateGraph(listID) {
	if (listID==0) return;
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxListDecGraph&modifID=$quantifModID&listID="+listID,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			var listData=XHR.responseText.split('::'); // list_name::sub_match(0/1)::prot1;prot2;...;protN
			var matchPattern=(listData[1]=='1')? '^###-' : null;
			addHighlighting(mainChart,listData[0],hColors[hColorIdx],{'-1':listData[2].split(';')},matchPattern);
			usedHighLights[listData[0]]=hColorIdx; usedColorIdx[hColorIdx]++;
			hColorIdx++;
			if (hColorIdx==hColors.length) hColorIdx=0;
		}
	}
	XHR.send(null);
}
function ajaxGoDecorateGraph(termIdStrg) {
	if (!termIdStrg) return;
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxTermProt&projectID=$projectID&goStrg="+termIdStrg,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			var termData=termIdStrg.split(',');
			addHighlighting(mainChart,termData[2],hColors[hColorIdx],{'-1':XHR.responseText.split(';')}$highlighMatchStrg);
			usedHighLights[listData[0]]=hColorIdx; usedColorIdx[hColorIdx]++;
			hColorIdx++;
			if (hColorIdx==hColors.length) hColorIdx=0;
		}
	}
	XHR.send(null);
}
|;
		}
		else {print "// AJAX --->\n";}
		#my $maxRatioAllowed=($ratioType eq 'Ratio')? scalar @{$refInfo->{RATIOS}} : (scalar @{$refInfo->{STATES}})-1;
		print qq
|var isNav = (navigator.appName.indexOf("Netscape") !=-1);
function ajaxProteinData(e,protID,ratioPos,action) {
	var displayDiv=document.getElementById('displayDIV');
	var divX,divY;
	if (e) { // called from protein list
		divX = (isNav)? e.pageX : event.clientX + document.body.scrollLeft; //divX-=5;
		divY = (isNav)? e.pageY : event.clientY + document.body.scrollTop; //divY+=10;
	}
	else {
		var graphDiv=document.getElementById('mainGraphDIV');
		divX = getDomOffset(graphDiv,'offsetLeft') + graphDiv.offsetWidth +2;
		divY = getDomOffset(graphDiv,'offsetTop');
	}
	displayDiv.style.left = divX + 'px';
	displayDiv.style.top = divY + 'px';
	displayDiv.style.display='block';

	var infoDiv=document.getElementById('infoDIV');
	infoDiv.innerHTML="<BR><BR><IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR>";

	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	//var extraParams=(action=='$ajaxPepAction')? '&pepType=$dispPepType' : '';
	//XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT="+action+"&id_ana=$anaID&id_quantif=$selQuantifID&view=$view&sort=$dispSort&foldChange=$dispFoldChange&pValue=$dispPvalue&coefVar=$dispCV&pepType=$dispPepType&numPep=$dispNumPep&id_prot="+protID+"&ratio="+ratioPos,true); //+extraParams    &stdDev=\$dispStdDev
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT="+action+"&id_ana=$anaID&id_quantif=$selQuantifID&pepType=$dispPepType&id_prot="+protID+"$dispMeasParamStrg&ratio="+ratioPos,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			if (action=='$ajaxPepAction') {
				var codeParts=XHR.responseText.split('#==========#');
				infoDiv.innerHTML=codeParts[0]; // HTML part
				eval(codeParts[1]); // javascript part
			}
			else {infoDiv.innerHTML=XHR.responseText;}
		}
	}
	XHR.send(null);
}
function ajaxUpdateGoTermList(goIdStrg) {
	var termsDiv=document.getElementById('goTermsDIV');
	if (!goIdStrg) {
		termsDiv.innerHTML='';
		return;
	}
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxGoTerms&projectID=$projectID&goStrg="+goIdStrg,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			termsDiv.innerHTML=XHR.responseText;
		}
	}
	XHR.send(null);
}
|;
		&promsMod::printAjaxManageSaveProteins($projectID,\%promsPath,'document.protForm.chkProt','ajaxUpdateRestrict');
		print qq
|function ajaxUpdateRestrict() {
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT=ajaxRestrictList&projectID=$projectID&restrictList=$restrictListID&noSelect=1",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			document.getElementById('restrictDIV').innerHTML='<SELECT name="restrictList" class="title3" onchange="document.displayForm.submit()"><OPTION value="">-= Select =-</OPTION>'+XHR.responseText;
			if ('$view'=='graph' && $dispResults==1) {
				document.getElementById('listDecoSPAN').innerHTML='<SELECT id="listDeco" onchange="ajaxListDecorateGraph(this.value,this.text)"><OPTION value="0">-= None =-</OPTION>'+XHR.responseText;
			}
		}
	}
	XHR.send(null);
}
</SCRIPT>
<BR>
<FORM name="displayForm" method="POST">
<INPUT type="hidden" name="ACT" value="$action">
<INPUT type="hidden" name="CALL" value="$call">
<INPUT type="hidden" name="id_ana" value="$anaID">
<INPUT type="hidden" name="id_quantif" value="$selQuantifID">
<INPUT type="hidden" name="displayRes" value="1">
|;
	}
	my $statusImage=($quantifStatus==-2)? 'lightRed1.gif' : ($quantifStatus==-1)? 'lightGray1.gif' : ($quantifStatus==0)? 'lightYellow1.gif' : 'lightGreen1.gif';
	#my $colspanStrg=($existReplicates)? ' colspan=2' : '';
	#my $colspanStrg=($existReplicates && $ratioSelection)? ' colspan=3' : ($existReplicates || $ratioSelection)? ' colspan=2' : '';
	my $colspanStrg=($ratioType eq 'None')? '' : ($quantifSoftware eq 'MaxQuant')? ' colspan=2' : ' colspan=3'; #($existReplicates)? ' colspan=3' : ' colspan=2';
	print qq
|<TABLE bgcolor=$darkColor border=0>
<TR><TH nowrap align=right>Labeling :</TH><TD bgcolor=$lightColor$colspanStrg valign=top>&nbsp;<B>$labelingName{$labelType}$swathStrg</B></TD></TR>
|;
	#my $maxDefaultRatios=$numStates-1; # all ratios with 1st reference
	#unless (scalar keys %dispRatios) { # setting default displayed ratios
	#	foreach my $pos (1..$maxDefaultRatios) {$dispRatios{$pos}=1;}
	#}
	my $colSpanCorrel=($quantifSoftware eq 'MaxQuant' || $quantifStatus > 0)? 1 : 2; # +/- correlation buttons
	foreach my $statePos (sort{$a<=>$b} keys %stateInfo) {
		#print "<TR><TH nowrap align=right>&nbsp;State #$statePos:</TH><TD bgcolor=$lightColor width=450><B>$stateInfo{$statePos}{NAME}</B>$stateInfo{$statePos}{REPLIC}</TD></TR>\n";
		my $replicateNames=($stateInfo{$statePos}{'NAME_LONG'})? $stateInfo{$statePos}{'NAME_LONG'} : '';
		print "<TR><TH nowrap align=right>&nbsp;State #$statePos :</TH><TD bgcolor=$lightColor colspan=$colSpanCorrel  style=\"min-width:700px\">&nbsp;";
		if ($ratioType eq 'None') {
			if ($action eq 'summary') {print "<IMG src=\"$promsPath{images}/$statusImage\"/>&nbsp;";}
			elsif (scalar keys %stateInfo > 1) {
				unless (scalar keys %dispStates) {%dispStates=(1=>1);} # set default
				my $selStatus=($dispStates{$statePos})? ' checked' : '';
				print "<INPUT type=\"checkbox\" name=\"dispStates\" value=\"$statePos\"$selStatus>&nbsp;";
			}
			else {print "<INPUT type=\"hidden\" name=\"dispStates\" value=\"$statePos\">&nbsp;"}
		}
		print "<A href=\"javascript:void(null)\" onmouseover=\"popup('$stateInfo{$statePos}{POPUP}')\" onmouseout=\"popout()\">" if $stateInfo{$statePos}{'POPUP'};
		print "<B>$stateInfo{$statePos}{NAME}$replicateNames</B>";
		print "</A>" if $stateInfo{$statePos}{'POPUP'};
		print $stateInfo{$statePos}{REPLIC};
		if ($quantifStatus > 0 && $ratioType=~/S\w+Ratio/ && $labelType eq 'FREE' && $quantifSoftware eq 'myProMS') {
			my $imageFile=($algoVersion==3)? "distriblog2FC_State$statePos.jpeg" : "Graph_Box_Hist_ratioprot_mean_State$statePos.png";
			print "&nbsp;<INPUT type=\"button\" class=\"font11\" value=\"Mean distrib.\" onclick=\"displayImageFromButton(this,'Protein Mean Distribution for State $stateInfo{$statePos}{NAME}','$imageFile')\">";
		}
		print "&nbsp;</TD>";
		if ($refInfo->{RATIOS}) {
			if ($statePos==1) { # && $existReplicates
				if ($quantifStatus > 0 && $quantifSoftware ne 'MaxQuant') {
					print qq |<TD bgcolor=$lightColor rowspan="$numStates" valign=middle>|;
					if ($ratioType=~/S\w+Ratio/) {
						#my $disabStrg=(-e "$resultDir/graph/$correlImgFileName")? '' : 'disabled';
						my $disabStrg=(-e "$resultDir/$correlMatFileName")? '' : 'disabled';
						print qq |&nbsp;<INPUT type="button" class="font11" value="Global correlation" style="width:120px" onclick="displayReplicCorrelation(this)" $disabStrg>&nbsp;<BR>|;
					}
					print qq |&nbsp;<INPUT type="button" class="font11" value="Global coeff. var." style="width:120px" onclick="displayImageFromButton(this,'Coefficient of variation between all replicates','pep_cv_hist.png')">&nbsp;| if -e "$resultDir/graph/pep_cv_hist.png"; # before myProMS v3! ($ratioType eq 'Ratio' || $supRatioReplicates); #displayPepVariation(this)
					print "</TD>\n";
				}
				if (scalar @{$refInfo->{RATIOS}} > 1) { # multiple ratios
					unless (scalar keys %dispRatios) {%dispRatios=(1=>1,2=>1,3=>1);} # set default
					print "<TH bgcolor=$lightColor rowspan=$numStates valign=top><TABLE bgcolor=$darkColor cellpadding=1><TR><TH bgcolor=$lightColor>";
					if ($action eq 'summary') {
						print "Ratios computed:</TH><TH>St.</TH>";
						foreach my $stPosT (2..$numStates) {print "<TH>#$stPosT</TH>"}
						print "</TR>\n";
						#my $numReferences=(($ratioType eq 'SimpleRatio' && $labelType ne 'FREE' && $quantifSoftware ne 'MaxQuant') || ($refInfo->{'SINGLE_REF'} && $refInfo->{'SINGLE_REF'}[0]==1))? 1 : $numStates-1;
my $numReferences=($refInfo->{'SINGLE_REF'} && $refInfo->{'SINGLE_REF'}[0]==1)? 1 : $numStates-1;
						my $ratioPos=0;
						foreach my $stPosR (1..$numReferences) {
							print "<TR>";
							if ($stPosR==1) {
								my $extraRefStrg.=($numReferences >= 4)? '<BR>e&nbsp;<BR>r&nbsp;<BR>e&nbsp;<BR>n&nbsp;<BR>c&nbsp;<BR>e&nbsp;<BR>s&nbsp;' : '';
								print "<TH bgcolor=\"$lightColor\" valign=\"middle\" rowspan=$numReferences align=right><FONT style=\"line-height:90%\">R&nbsp;<BR>e&nbsp;<BR>f&nbsp;$extraRefStrg</FONT></TH>";
							}
							print "<TH>#$stPosR</TH>";
							foreach my $stPosT (2..$numStates) {
								if ($stPosT<=$stPosR) {print "<TH bgcolor=$lightColor></TH>"; next;}
								$ratioPos++;
								my $ratioTag=($ratioType eq 'SuperRatio' && $ratioPos > $numReferences)? encode_utf8('°') : '';
								my $ratioName="$stateInfo{$stPosT}{NAME}$ratioTag/$stateInfo{$stPosR}{NAME}$ratioTag";
								print "<TH bgcolor=$lightColor onmouseover=\"popup('<B>$ratioName</B>')\" onmouseout=\"popout()\"><IMG id=\"menuRatio_$ratioPos\" src=\"$promsPath{images}/$statusImage\"";
								print " onclick=\"displayRatioStat('menuRatio_$ratioPos','$ratioCode{$ratioPos}','$ratioName')\"" if ($quantifSoftware eq 'myProMS' && $quantifStatus > 0);
								print "></TH>\n";
							}
							print "</TR>\n";
						}
					}
					else {
						print "Ratios displayed:</TH><TH>St.</TH>";
						foreach my $stPosT (1..$numStates) {print "<TH>#$stPosT</TH>"}
						print "</TR>\n";
						#my $singleRef=(($ratioType eq 'SimpleRatio' && $labelType ne 'FREE' && $quantifSoftware ne 'MaxQuant') || ($refInfo->{'SINGLE_REF'} && $refInfo->{'SINGLE_REF'}[0]==1))? 1 : 0;
my $singleRef=($refInfo->{'SINGLE_REF'} && $refInfo->{'SINGLE_REF'}[0]==1)? 1 : 0;
						my $ratioPos=0;
						my %revRatioPos;
						foreach my $stPosR (1..$numStates) {
							print "<TR>";
							if ($stPosR==1) {print "<TH bgcolor=$lightColor rowspan=$numStates align=right valign=top nowrap><FONT color=\"#00BB00\">&nbsp;Computed ratios</FONT>&nbsp;<BR><FONT color=\"#DD0000\">&nbsp;Reversed ratios</FONT>&nbsp;</TH>";}
							print "<TH>#$stPosR</TH>";
							foreach my $stPosT (1..$numStates) {
								if ($stPosT==$stPosR) {print "<TH bgcolor=$lightColor></TH>\n";}
								elsif ($stPosT > $stPosR) {
									$ratioPos++;
									$revRatioPos{"$stPosT:$stPosR"}=-$ratioPos;
									if ($singleRef && $stPosR > 1) {print "<TH></TH>\n";}
									else {
										my $ratioTag=($ratioType eq 'SuperRatio' && $stPosR > 1)? encode_utf8('°') : '';
										my $ratioName="$stateInfo{$stPosT}{NAME}$ratioTag/$stateInfo{$stPosR}{NAME}$ratioTag";
										print "<TH bgcolor=\"#00DD00\" onmouseover=\"popup('<B>$ratioName</B>')\" onmouseout=\"popout()\"><INPUT type=\"checkbox\" name=\"dispRatios\" id=\"$ratioPos\" value=\"$ratioPos\" onclick=\"if (this.checked) {document.getElementById('-$ratioPos').checked=false;}\"";
										print ' checked' if $dispRatios{$ratioPos};
										print "/>";
										print "<IMG id=\"menuRatio_$ratioPos\" src=\"$promsPath{images}/plus.gif\" onclick=\"displayRatioStat('menuRatio_$ratioPos','$ratioCode{$ratioPos}','$ratioName')\">" if $quantifSoftware eq 'myProMS';
										print "</TH>\n";
									}
								}
								else { # reverse ratio
									if ($singleRef && $stPosT > 1) {print "<TH></TH>\n";}
									else {
										my $revRatPos=$revRatioPos{"$stPosR:$stPosT"};
										my $ratioTag=($ratioType eq 'SuperRatio' && $stPosT > 1)? encode_utf8('°') : '';
										my $ratioName="$stateInfo{$stPosT}{NAME}$ratioTag/$stateInfo{$stPosR}{NAME}$ratioTag";
										print "<TH bgcolor=\"#DD0000\" onmouseover=\"popup('<B>$ratioName</B>')\" onmouseout=\"popout()\"><INPUT type=\"checkbox\" name=\"dispRatios\" id=\"$revRatPos\" value=\"$revRatPos\" onclick=\"if (this.checked) {document.getElementById((-1*$revRatPos)+'').checked=false;}\"";
										print ' checked' if $dispRatios{$revRatPos};
										print "/></TH>\n";
									}
								}
							}
							print "</TR>\n";
						}
					}
					print "</TABLE>\n";
				}
				else { # single ratio
					print "<TD bgcolor=$lightColor rowspan=2 class=\"title3\">";
					if ($action eq 'summary') {
						print "&nbsp;<IMG src=\"$promsPath{images}/$statusImage\">";
						if ($quantifSoftware eq 'myProMS') {
							print "<A id=\"menuRatio_1\" href=\"javascript:displayRatioStat('menuRatio_1','$ratioCode{1}','$stateInfo{2}{NAME}/$stateInfo{1}{NAME}')\">$stateInfo{2}{NAME}/$stateInfo{1}{NAME}</A>&nbsp;";
						}
						else {print "$stateInfo{2}{NAME}/$stateInfo{1}{NAME}&nbsp;";}
					}
					else {
						print "<TABLE cellspacing=0 cellpadding=0 ><TR><TD align=right class=\"title3\">&nbsp;<FONT color=\"#00BB00\">Computed ratio:</FONT></TD><TD class=\"title3\"><INPUT type=\"checkbox\" name=\"dispRatios\" id=\"1\" value=\"1\" onclick=\"if (this.checked) {document.getElementById('-1').checked=false;}\"";
						print ' checked' if $dispRatios{1};
						print "/>";
						if ($quantifSoftware eq 'myProMS') {
							print "<A id=\"menuRatio_1\" href=\"javascript:displayRatioStat('menuRatio_1','$ratioCode{1}','$stateInfo{2}{NAME}/$stateInfo{1}{NAME}')\">$stateInfo{2}{NAME}/$stateInfo{1}{NAME}</A>&nbsp;";
						}
						else {print "$stateInfo{2}{NAME}/$stateInfo{1}{NAME}&nbsp;";}
						print "</TD></TR><TR><TD align=right class=\"title3\">&nbsp;<FONT color=\"#DD0000\">Reversed ratio:</FONT></TD><TD class=\"title3\"><INPUT type=\"checkbox\" name=\"dispRatios\" id=\"-1\" value=\"-1\" onclick=\"if (this.checked) {document.getElementById('1').checked=false;}\"";
						print ' checked' if $dispRatios{-1};
						print "/>$stateInfo{1}{NAME}/$stateInfo{2}{NAME}&nbsp;</TD></TR></TABLE>\n";
					}
					print "</TD>";
				}
			}
		}
		print "</TR>\n";
	}

	#my $pepString=($quantifType !~ /XIC|SIN|EMPAI/)? "&nbsp;&nbsp;&bull;$pepChargeStrg&nbsp;&nbsp;&nbsp;&bull;$pepSourceStrg&nbsp;" : ''; # Change on 29/10/12
	my $pepString=($quantifSoftware eq 'myProMS')? "&nbsp;&nbsp;&bull;$pepChargeStrg&nbsp;&nbsp;&nbsp;&bull;$pepSourceStrg&nbsp;" : '';
	#my $pValueMeaningStrg=($algoVersion < 3)? '' : ($refInfo->{'RESIDUAL_VAR'}[0] eq 'biological')? ''
	my $infRatioString=($refInfo->{'MIN_INF_RATIOS'} && $refInfo->{'MIN_INF_RATIOS'}[0])? 'Avoided whenever possible' : 'Not avoided';
	$ptmQuantifStrg='<BR>&nbsp;&bull;'.$ptmQuantifStrg if $ptmQuantifStrg;
	my $protSelectionStrg='<B>All</B> visible proteins'; # default
	if ($refInfo->{'PROTEINS'}) {
		$protSelectionStrg=($refInfo->{'PROTEINS'}[0] eq 'exclude')? '<B>Exclude</B> proteins in List ' : '<B>Restrict</B> to proteins in List ';
		(my $listID=$refInfo->{'PROTEINS'}[1])=~s/#//;
		my ($listStrg)=$dbh->selectrow_array("SELECT CONCAT(T.NAME,' > ',L.NAME) FROM CATEGORY L,CLASSIFICATION T WHERE L.ID_CLASSIFICATION=T.ID_CLASSIFICATION AND ID_CATEGORY=$listID");
		$protSelectionStrg.="<B>$listStrg</B>";
	}
	print "<TR><TH align=right nowrap valign=top>&nbsp;Peptide selection :</TH><TD nowrap bgcolor=$lightColor$colspanStrg>&nbsp;";
	if ($refInfo->{PEPTIDES}) {print "&bull;$pepRangeStrg&nbsp;&nbsp;&nbsp;&bull;$pepMissCutStrg&nbsp;&nbsp;&nbsp;&bull;$pepPtmStrg&nbsp;$pepString&nbsp;$ptmQuantifStrg";}
	else {print 'All';}
	print qq
|</TD></TR>
<TR><TH align=right nowrap valign=top>&nbsp;Protein selection :</TH><TD nowrap bgcolor=$lightColor$colspanStrg>&nbsp;$protSelectionStrg&nbsp;(<B>$protQuantifiedStrg</B>)</TD></TR>
<TR><TH align=right nowrap valign=top>&nbsp;Quantification settings :|;
	if ($ratioType eq 'Ratio') {
		print qq
|<BR><INPUT type="button" id="moreSettings" class="font11" value="More settings" onclick="updateSettings('more')"/><INPUT type="button" id="lessSettings" class="font11" value="Less settings" style="display:none" onclick="updateSettings('less')"/>|;
	}
	#my $quantifSoftwareStrg=($quantifSoftware eq 'myProMS')? 'myProMS v'.$algoVersion : $quantifSoftware;
	print qq
|</TH><TD bgcolor=$lightColor$colspanStrg valign=top>
	&nbsp;-<B>Quantification method:</B> $quantifMethDesc [$ratioTypeStrg$topNString - Software: $quantifSoftware$softVersionStrg]<BR>
|;
	if ($algoVersion==3) {
		print "&nbsp;-<B>p-value</B> estimates <B>$refInfo->{'RESIDUAL_VAR'}[0]</B> variation [Correction: <B>",$pValueCorrections{$refInfo->{'FDR_CONTROL'}[0]},"</B>]<BR>\n";
	}
	if ($quantifSoftware ne 'MaxQuant') {print "&nbsp;-<B>Infinite ratios:</B> $infRatioString\n";}
	if ($refInfo->{RATIOS}) {
		print "<TABLE border=0 cellspacing=0><TR valign=top><TD nowrap>&nbsp;-<B>Bias correction:</B></TD>$biasCorrectStrg</TABLE>\n";
	}
	if ($ratioType eq 'Ratio') {
		print qq
|	<DIV id="advancedSetDIV" style="display:none">
	&nbsp;&bull;<B><U>Advanced settings:</U></B>
|;
		if ($refInfo->{THRESHOLD_CV}) { # replicates used TODO: Make sure THRESHOLD_CV is reported for design quantif
			my $biasCorrStrg='<B>Variation coefficient threshold between replicates:</B> ';
			$biasCorrStrg.=($refInfo->{THRESHOLD_CV}[0] eq 'FALSE')? 'Auto' : $refInfo->{THRESHOLD_CV}[1];
			print "<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-$biasCorrStrg\n";
		}
		print qq
|	<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-<B>FDR control:</B> $fdrStrg
	<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-<B>p-value threshold for outlier detection:</B> $pvalStrg
|;
		print "<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-<B>Alternative hypothesis for comparison:</B> $alterStrg\n" if $refInfo->{'ALTER'};
		print "<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-<B>Confidence interval on protein abundance:</B> $refInfo->{CONFIDENCE_LEVEL}[0]\n" if $refInfo->{'CONFIDENCE_LEVEL'};
		print "</DIV>\n";
	}
	print qq
|</TD></TR>
|;
	if ($action eq 'summary') {
		my $statusStrg=($quantifStatus==-2)? '<FONT color=#DD0000>Failed</FONT> (Click on "Monitor Quantifications" for more information)' : ($quantifStatus==-1)? 'Not launched yet' : ($quantifStatus==0)? 'On-going' : 'Finished';
		$updateDate=&promsMod::formatDate($updateDate);
		print qq
|<TR><TH align=right nowrap valign=top>&nbsp;Creation date :</TH><TD bgcolor=$lightColor$colspanStrg valign=top>&nbsp;$updateDate by <B>$updateUser</B></TD></TR>
<TR><TH align=right nowrap valign=top>&nbsp;Status :</TH><TH align="left" bgcolor=$lightColor$colspanStrg valign=top>&nbsp;$statusStrg</TH></TR>
</TABLE>
<BR>
|;
		$dbh->disconnect;

		if ($quantifStatus==-1 || $quantifStatus==0) { # on-going: reload every 5 sec
			print qq
|<SCRIPT type="text/javascript">
setTimeout(function(){parent.optionFrame.location.reload();},5000);
</SCRIPT>
|;
		}

		&endHTML;
		exit;
	}
	#my ($selGraphView,$selListView,$listDivVis)=($view eq 'graph')? (' selected','','none') : ('',' selected','block');
	my $trueFilterStyle='style="font-weight:bold;color:#DD0000"';
	my ($selListView,$foldChPvalStyle)=($view eq 'list')? (' selected',$trueFilterStyle) : ('','style="font-weight:bold"'); # 'graph' is default
	my ($selFcUp,$selFcDown)=($foldChangeType eq 'up')? (' selected','') : ($foldChangeType eq 'down')? ('',' selected') : ('','');
	$dispPvalue=($dispPvalue)? $dispPvalue : $refInfo->{FDR_ALPHA}[0]/100;
	#my ($selDistPepType,$selRazPepType,$selUniPepType)=($numPepCode eq 'DIST_PEP_USED')? (' selected','','') : ($dispPepType eq 'RAZ_UNI_PEP')? ('',' selected','') :  ($dispPepType eq 'UNIQUE_PEP')? ('','',' selected') : ('','','');
	my ($selDistPepType,$selRazPepType,$selUniPepType,$selRatioCountType)=($numPepCode eq 'DIST_PEP_USED')? (' selected','','','') : ($numPepCode eq 'RAZ_UNI_PEP')? ('',' selected','','') :  ($numPepCode eq 'UNIQUE_PEP')? ('','',' selected','') :  ($numPepCode eq 'NUM_PEP_USED')? ('','','',' selected') : ('','','','');
	$selDistPepType=($dispPepType eq 'distinct')? ' selected' : ''; # 'all' is default

	print qq
|<TR><TH align=right nowrap valign=top>&nbsp;<FONT class="title2">Display :</FONT></TH><TD bgcolor=$lightColor$colspanStrg><TABLE border=0 cellspacing=0 cellpadding=0>
	<TR><TD colspan=2 nowrap><FONT class="title3">&bull;View:</FONT><SELECT name="view" class="title3" onchange="document.displayForm.submit()"><OPTION value="graph">Graphical</OPTION><OPTION value="list"$selListView>List</OPTION></SELECT>
|;
	if ($ratioType eq 'None') {
		print "&nbsp;&nbsp;&nbsp;<FONT style=\"font-weight:bold\">&bull;Measure:</FONT><SELECT name=\"dispMeasure\" style=\"font-weight:bold\" onchange=\"document.displayForm.submit()\"><OPTION value=\"MQ_INT\">Intensity</OPTION>";
		my $sthMeas=$dbh->prepare("SELECT 1 FROM PROTEIN_QUANTIFICATION P,QUANTIFICATION_PARAMETER Q WHERE P.ID_QUANTIF_PARAMETER=Q.ID_QUANTIF_PARAMETER AND P.ID_QUANTIFICATION=$selQuantifID AND CODE=? LIMIT 0,1");
		foreach my $refCode (['MQ_IBAQ','iBAQ'],['MQ_LFQ','LFQ'],['MQ_SC','MS/MS count']) {
			$sthMeas->execute($refCode->[0]);
			my ($OK)=$sthMeas->fetchrow_array;
			my $selStrg=($OK && $refCode->[0] eq $dispMeasure)? 'selected' : '';
			my $disabStrg=($OK)? '' : 'disabled';
			print "<OPTION value=\"$refCode->[0]\" $disabStrg $selStrg>$refCode->[1]</OPTION>";
		}
		$sthMeas->finish;
		print "</SELECT>\n";
	}
	else {
		print "&nbsp;&nbsp;&nbsp;<FONT $foldChPvalStyle>&bull;";
		if ($quantifSoftware eq 'MaxQuant') {
			my $selRawStrg=($foldChangeRaw)? ' selected' : '';
			print "</FONT><SELECT name=\"foldChgRaw\" $foldChPvalStyle><OPTION value=\"0\">Norm.</OPTION><OPTION value=\"1\"$selRawStrg>Raw</OPTION></SELECT><FONT $foldChPvalStyle>";
		}
		print qq
|Ratio:</FONT><SELECT name="foldChgType" $foldChPvalStyle><OPTION value="abs">Up &ge; or Down &le; 1/</OPTION><OPTION value="up"$selFcUp>Up &ge;</OPTION><OPTION value="down"$selFcDown>Down &le; 1/</OPTION></SELECT><INPUT type="text" name="foldChange" value="$dispFoldChange" size=2 $foldChPvalStyle/>
|;
		if ($quantifSoftware eq 'MaxQuant') {print "<INPUT type=\"hidden\" name=\"pValue\" value=\"1\"/>";}
		else {print qq |&nbsp;&nbsp;&nbsp;<FONT $foldChPvalStyle>&bull;p-value &le;</FONT><INPUT type="text" name="pValue" value="$dispPvalue" size=5 $foldChPvalStyle/>|;}
		if ($quantifType eq 'PROT_RATIO_PEP') {
			if ($quantifSoftware eq 'MaxQuant') {print "<INPUT type=\"hidden\" name=\"coefVar\" value=\"0\"/>";}
			else {print qq |&nbsp;&nbsp;&nbsp;<FONT $trueFilterStyle onmouseover="popup('<B>Use 0% for no filtering</B>')" onmouseout="popout()" >&bull;Coeff. var.<SUP>*</SUP> &le;</FONT><INPUT type="text" name="coefVar" value="$dispCV" size=2 $trueFilterStyle/><FONT $trueFilterStyle>%</FONT>|;}
		}
		#else {} # TNPQ
	}
	# Peptide selection
	print qq |&nbsp;&nbsp;&nbsp;<FONT $trueFilterStyle>&bull;</FONT><SELECT name="pepType" $trueFilterStyle>|; # onchange="document.displayForm.submit()"
	if ($quantifSoftware eq 'MaxQuant') {
		my $ratioPepOptStrg=($ratioType eq 'None')? '' : "<OPTION value=\"NUM_PEP_USED\"$selRatioCountType>Ratio counts</OPTION>";
		print qq
|<OPTION value="PEPTIDES">All peptides</OPTION><OPTION value="RAZ_UNI_PEP"$selRazPepType>Razor+unique peptides</OPTION><OPTION value="UNIQUE_PEP"$selUniPepType>Unique peptides</OPTION>$ratioPepOptStrg
</SELECT> <FONT $trueFilterStyle>&ge;</FONT>
|;
	}
	else {
		print qq
|<OPTION value="all">All</OPTION><OPTION value="distinct"$selDistPepType>Dist.</OPTION>
</SELECT> <FONT $trueFilterStyle>peptides &ge;</FONT>
|;
	}
	print qq
|<INPUT type="text" name="numPep" value="$dispNumPep" size=2 $trueFilterStyle/>
&nbsp;&nbsp;&nbsp;<B>&bull;Sort by:</B><SELECT name="sort" id="sort" onchange="updateSortBy()" style="font-weight:bold">
|;
#<TR><TH align=right nowrap valign=top>&nbsp;<FONT class="title2">Visualization :</FONT></TH><TD bgcolor=$lightColor$colspanStrg>
#	<TABLE cellspacing=0 cellpadding=0><TR><TD valign=top><SELECT class="title2" onchange="displayQuantification($selQuantifID,'$action',this.value)"><OPTION value="graph"$selGraphView>Volcano plot</OPTION><OPTION value="list"$selListView>Full list</OPTION></SELECT></TD>
#	<TD><DIV id="listView" style="display:$listDivVis"><TABLE cellscaping=0 cellpadding=0><TR>
#		<TH align=left nowrap>&nbsp;&bull;<FONT class="title3">Flag:</FONT> Fold change &ge; <INPUT type="text" id="foldChange" value="$dispFoldChange" size=2/>
#<!--	and Std. deviation &le; <INPUT type="text" id="stdDev" value="$dispStdDev" size=3/></TH> -->
#	and p-value &le; <INPUT type="text" id="pValue" value="$dispPvalue" size=5/></TH>
#		<TH rowspan=2><INPUT type="button" value="Update" class="title3" onclick="updatedListDisplay()"/></TH>
#	</TR><TR>
#	<TH align=left>&nbsp;&bull;<FONT class="title3">Sort by:</FONT><SELECT id="sort" onchange="updatedListDisplay()">
#|;
	my @sortOptions=(['identifier','Identifiers']);
	if ($ratioType eq 'None') {
		if (scalar keys %stateInfo == 1) {
			push @sortOptions,['state_1','Measure'];
		}
		else {
			foreach my $statePos (sort{$a<=>$b} keys %stateInfo) {
				push @sortOptions,["state_$statePos","Measure: '$stateInfo{$statePos}{NAME}'"];
			}
		}
	}
	else {
		if (scalar @{$refInfo->{RATIOS}} == 1) {
			push @sortOptions,['ratio_1','Ratio'];
		}
		else {
			my $ratioPos=0;
			foreach my $ratioData (@{$refInfo->{RATIOS}}) {
				$ratioPos++;
				next unless $dispRatios{$ratioPos};
				my ($testStatePos,$refStatePos)=split(/\//,$ratioData);
				#if ($parentQuantifType && $parentQuantifType eq 'XIC'){#}
				my $ratioTag='';
				if ($designID) {
					if ($testStatePos=~/%/) {
						$ratioTag=encode_utf8('°');
						$testStatePos=~s/%\d+//;
						$refStatePos=~s/%\d+//;
					}
					$testStatePos=$condToState{$testStatePos};
					$refStatePos=$condToState{$refStatePos};
				}
				push @sortOptions,["ratio_$ratioPos","Ratio: '$stateInfo{$testStatePos}{NAME}$ratioTag/$stateInfo{$refStatePos}{NAME}$ratioTag'"];
			}
		}
	}
	push @sortOptions,(['peptide','Peptide filter'],['mw','Molecular weight']);
	foreach my $refOption (@sortOptions) {
		print "<OPTION value=\"$refOption->[0]\"";
		print ' selected' if $refOption->[0] eq $dispSort;
		print ">$refOption->[1]</OPTION>";
	}
	print qq
|</SELECT>&nbsp;&nbsp;&nbsp;</TD><TD rowspan=2><INPUT type="submit" value="Update" style="font-weight:bold"/></TD></TR>
	<TR><TD nowrap width=10%><FONT class="title3">&bull;Restrict to proteins in List:&nbsp;</FONT></TD><TD><DIV id="restrictDIV"><!--Custom list selection comes here with window.onload --></DIV></TD></TR>
</TABLE></TD></TR>
|;
##	if (scalar keys %{$refGoAnalyses}) {
##		print qq
##|<TR><TH align=right nowrap valign=top>&nbsp;<FONT class="title3">GO decoration :</FONT></TH><TD bgcolor=$lightColor$colspanStrg><TABLE border=0 cellspacing=0 cellpadding=0>
##	<TR><TD nowrap><SELECT id="goAna" class="title3" onchange="updateGoAspects(this.value)"><OPTION value="0">-= None =-</OPTION>
##|;
##		foreach my $goID (sort{lc($refGoAnalyses->{$a}[0]) cmp lc($refGoAnalyses->{$b}[0])} keys %{$refGoAnalyses}) {
##			print "<OPTION value=\"$goID\">$refGoAnalyses->{$goID}[0]</OPTION>";
##		}
##		print qq
##|</SELECT><FONT class="title3">:&nbsp;</FONT><SELECT id="goAspect" class="title3" style="visibility:hidden" onchange="ajaxUpdateGoTermList(this.value)"><OPTION value="">-= Select Aspect =-</OPTION></SELECT>&nbsp;</TD>
##<TD><DIV id="goTermsDIV"></DIV></TD></TR></TD></TABLE></TD></TR>
##|;
##	}

	if ($view eq 'graph' && $dispResults) {
		print qq
|<TR><TH align=right nowrap valign=top>&nbsp;<FONT class="title3">Chart highlighting :</FONT></TH><TD bgcolor=$lightColor$colspanStrg><TABLE border=0 cellspacing=0 cellpadding=0>
	<TR><TH align=right nowrap>&nbsp;Custom lists:</TH><TD nowrap colspan=3><SPAN id="listDecoSPAN"></SPAN></TD></TR>
|;
		if (scalar keys %{$refGoAnalyses}) {
			print qq
|	<TR><TH align=right nowrap>&nbsp;Gene Ontology:</TH><TD nowrap><SELECT id="goAna" onchange="updateGoAspects(this.value)"><OPTION value="0">-= None =-</OPTION>
|;
			foreach my $goID (sort{lc($refGoAnalyses->{$a}[0]) cmp lc($refGoAnalyses->{$b}[0])} keys %{$refGoAnalyses}) {
				print "<OPTION value=\"$goID\">$refGoAnalyses->{$goID}[0]</OPTION>";
			}
			print qq
|</SELECT><FONT class="title3">:&nbsp;</FONT><SELECT id="goAspect" style="visibility:hidden" onchange="ajaxUpdateGoTermList(this.value)"><OPTION value="">-= Select Aspect =-</OPTION></SELECT>&nbsp;</TD>
<TD><DIV id="goTermsDIV"></DIV></TD></TR>
|;
		}
		print "</TABLE></TD></TR>\n";
	}

	print qq
|</TABLE>
</FORM>
<BR>
|;
	unless ($dispResults) { # form was not submitted (1st time displayed)
		$dbh->disconnect;
		print qq
|<SCRIPT type="text/javascript">
window.onload=function() {
	ajaxUpdateRestrict();
}
</SCRIPT>
<FONT class="title1"><BR><BR><BR>Choose settings and click on "Update" to display data</FONT>
|;
		&endHTML;
		exit;
	}
	return %stateInfo;
} # end of &printQuantificationSummary


sub deleteQuantification {
	my ($anaID,$quantifID)=@_;
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
	my ($projID)=&promsMod::getProjectID($dbh,$anaID,'analysis');
	#my ($status)=$dbh->selectrow_array("SELECT STATUS FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID");

	#<Peptide sets
	#my $sthPSet=$dbh->prepare("SELECT DISTINCT ID_PEPTIDE_SET FROM PEPSET_QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID");
	#my $sthUsed=$dbh->prepare("SELECT COUNT(*) FROM PEPSET_QUANTIFICATION WHERE ID_QUANTIFICATION != $quantifID AND ID_PEPTIDE_SET=?");
	#my $sthDelPSetQuant=$dbh->prepare("DELETE FROM PEPSET_QUANTIFICATION WHERE ID_PEPTIDE_SET=?");
	#my $sthDelPepPSet=$dbh->prepare("DELETE FROM PEPSET_PEPTIDE WHERE ID_PEPTIDE_SET=?");
	#my $sthDelPSet=$dbh->prepare("DELETE FROM PEPTIDE_SET WHERE ID_PEPTIDE_SET=?");
	#$sthPSet->execute;
	#my $count=0;
	#while (my ($pepSetID)=$sthPSet->fetchrow_array) {
	#	$count++;
	#	if ($count==250) {
	#		$count=0;
	#		print '.';
	#	}
	#	$sthDelPSetQuant->execute($pepSetID);
	#	$sthUsed->execute($pepSetID);
	#	my ($isUsed)=$sthUsed->fetchrow_array;
	#	next if $isUsed; # peptideSet is used by another quantification
	#	$sthDelPepPSet->execute($pepSetID);
	#	$sthDelPSet->execute($pepSetID);
	#}
	#$sthPSet->finish;
	#$sthUsed->finish;
	#$sthDelPSetQuant->finish;
	#$sthDelPepPSet->finish;
	#$sthDelPSet->finish;

	#<Proteins
	$dbh->do("DELETE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID");
	print '.';

	#<Parent & analysis
	$dbh->do("DELETE FROM PARENT_QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID");
	$dbh->do("DELETE FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID");
	$dbh->do("DELETE FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID");
	print '.';
	$dbh->commit;
	$dbh->disconnect;

	if (-e "$promsPath{quantification}/project_$projID") {
		#remove_tree("$promsPath{quantification}/project_$projID/quanti_$quantifID") if -e "$promsPath{quantification}/project_$projID/quanti_$quantifID";
		rmtree("$promsPath{quantification}/project_$projID/quanti_$quantifID") if -e "$promsPath{quantification}/project_$projID/quanti_$quantifID";
		my @remaindDirs = glob "$promsPath{quantification}/project_$projID/quanti_*";
		#remove_tree("$promsPath{quantification}/project_$projID") unless scalar @remaindDirs;
		rmtree("$promsPath{quantification}/project_$projID") unless scalar @remaindDirs;
	}

	# Remove log file of this quantification
	unlink "$promsPath{logs}/quanti_$quantifID.log" if -e "$promsPath{logs}/quanti_$quantifID.log";

	print " Done</FONT>\n";
	sleep 2;

	print qq
|<SCRIPT type="text/javascript">
//window.location="$promsPath{cgi}/showProtQuantification.cgi?id_ana=$anaID";
parent.optionFrame.location.reload(); // will update report button if only primary peptide quantif remains
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}

sub sortProteins {
	my ($sort,$refProteinInfo,$refQuantifValues,$firstRatioUsed)=@_;
	my ($pID1,$modStrg1)=($a=~/^(\d+)(.*)/); $modStrg1='' unless $modStrg1;
	my ($pID2,$modStrg2)=($b=~/^(\d+)(.*)/); $modStrg2='' unless $modStrg2;
	if ($sort eq 'identifier') {lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))}
	elsif ($sort eq 'mw') {lc($refProteinInfo->{$pID2}[1]) <=> lc($refProteinInfo->{$pID1}[1]) || lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))}
	elsif ($sort =~ /ratio_(-*\d+)/) {
		my $ratioPos=$1;
		if ($refQuantifValues->{$a} && $refQuantifValues->{$b}) {
			#$refQuantifValues->{$b}{'RATIO'}{$ratioPos}<=>$refQuantifValues->{$a}{'RATIO'}{$ratioPos} || $refQuantifValues->{$b}{'NUM_PEP_USED'}{$ratioPos}<=>$refQuantifValues->{$a}{'NUM_PEP_USED'}{$ratioPos} || lc($refProteinInfo->{$a}[0]) cmp lc($refProteinInfo->{$b}[0])
			&sortCheckDefined($ratioParamCode,$refQuantifValues,$ratioPos,-1) || &sortCheckDefined('NUM_PEP_USED',$refQuantifValues,$ratioPos,-1) || &sortCheckDefined('NUM_PEP_TOTAL',$refQuantifValues,$ratioPos,-1) || lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))
		}
		elsif ($refQuantifValues->{$a}) {-1}
		elsif ($refQuantifValues->{$b}) {1}
		else {lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))}
	}
	elsif ($sort =~ /state_(\d+)/) {
		my $statePos=$1;
		if ($refQuantifValues->{$a} && $refQuantifValues->{$b}) {
			#$refQuantifValues->{$b}{'RATIO'}{$ratioPos}<=>$refQuantifValues->{$a}{'RATIO'}{$ratioPos} || $refQuantifValues->{$b}{'NUM_PEP_USED'}{$ratioPos}<=>$refQuantifValues->{$a}{'NUM_PEP_USED'}{$ratioPos} || lc($refProteinInfo->{$a}[0]) cmp lc($refProteinInfo->{$b}[0])
			&sortCheckDefined($dispMeasure,$refQuantifValues,$statePos,-1) || &sortCheckDefined('PEPTIDES',$refQuantifValues,$statePos,-1) || &sortCheckDefined('RAZ_UNI_PEP',$refQuantifValues,$statePos,-1) || &sortCheckDefined('UNIQUE_PEP',$refQuantifValues,$statePos,-1) || lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))
		}
		elsif ($refQuantifValues->{$a}) {-1}
		elsif ($refQuantifValues->{$b}) {1}
		else {lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))}
	}
	# defaults to peptide NUM_PEP_USED or DIST_PEP_USED
	else {
		if ($refQuantifValues->{$a} && $refQuantifValues->{$b}) {
			#$refQuantifValues->{$b}{$numPepCode}{1}<=>$refQuantifValues->{$a}{$numPepCode}{1} || $refQuantifValues->{$b}{'RATIO'}{1}<=>$refQuantifValues->{$a}{'RATIO'}{1} || lc($refProteinInfo->{$a}[0]) cmp lc($refProteinInfo->{$b}[0])
			&sortCheckDefined($numPepCode,$refQuantifValues,$firstRatioUsed,-1) || &sortCheckDefined($ratioParamCode,$refQuantifValues,$firstRatioUsed,-1) || lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))
		}
		elsif ($refQuantifValues->{$a}) {-1}
		elsif ($refQuantifValues->{$b}) {1}
		else {lc($refProteinInfo->{$pID1}[0]) cmp lc($refProteinInfo->{$pID2}[0]) || &promsMod::sortSmart(&promsMod::preparePtmString($modStrg1),&promsMod::preparePtmString($modStrg2))}
	}
}
sub sortCheckDefined {
	my ($qCode,$refQuantifValues,$ratioPos,$order)=@_;
	$order=1 if !$order;
	if (defined $refQuantifValues->{$a}{$qCode}{$ratioPos} && defined $refQuantifValues->{$b}{$qCode}{$ratioPos}) {
		if ($order==1) {$refQuantifValues->{$a}{$qCode}{$ratioPos}<=>$refQuantifValues->{$b}{$qCode}{$ratioPos}}
		else {$refQuantifValues->{$b}{$qCode}{$ratioPos}<=>$refQuantifValues->{$a}{$qCode}{$ratioPos}}
	}
	elsif (defined $refQuantifValues->{$a}{$qCode}{$ratioPos}) {1*$order}
	elsif (defined $refQuantifValues->{$b}{$qCode}{$ratioPos}) {-1*$order}
	else {0}
}
#sub preparePtmString { # reverse: S145.T234 -> 145S234T for better sort by sortSmart
#	my ($ptmStrg)=@_;
#	return '' unless $ptmStrg;
#	$ptmStrg=~s/^-//;
#	my $revStrg='';
#	foreach my $modPos (split(/\./,$ptmStrg)) {
#		my ($res,$pos)=($modPos=~/^(.)(\d+)/);
#		$revStrg.=$2.$1;
#	}
#	return $revStrg;
#}

sub displayVolcanoPlot {
	my ($refLabelingInfo,$refStateInfo,$refQuantifValues,$refProteinInfo,$pvalueCode)=@_;
	print qq
|<SCRIPT type="text/javascript">
var mainChart,VP;
window.onload=function() {
	ajaxUpdateRestrict();
	VP=new volcanoPlot({div:'mainGraphDIV',
						width:500,height:450,
						foldChange:$dispFoldChange,
						pValue:$dispPvalue,
						allowHighlight:true,
						updateHighlight:{callback:updateUsedColors,editable:false},
						pointOnclick:proteinLink,
						pointOnList:ajaxListSelectedProteins,
						exportAsImage:['Export as image','VolcanoPlot','./exportSVG.cgi']
						});
	mainChart=VP; // needed for highlighting
|;
	my $pepSizeStrg=($dispPepType eq 'all')? 'All pep. used' : 'Dist. pep. used';
	my $ratioIdx=0;
	my (%ratioPos2index,%supRatios);
	foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
		my ($testStatePos,$refStatePos)=($ratioPos > 0)? split(/\//,$refLabelingInfo->{'RATIOS'}[$ratioPos-1]) : reverse(split(/\//,$refLabelingInfo->{'RATIOS'}[abs($ratioPos)-1]));
		my $ratioTag='';
		if ($designID) {
			if ($testStatePos=~/%/) {
				$ratioTag=encode_utf8('°');
				$testStatePos=~s/%\d+//;
				$refStatePos=~s/%\d+//;
				$supRatios{$ratioPos}=1;
			}
			$testStatePos=$condToState{$testStatePos};
			$refStatePos=$condToState{$refStatePos};
		}
		print "\tVP.addDataSet($ratioIdx,{name:'$refStateInfo->{$testStatePos}{NAME}$ratioTag/$refStateInfo->{$refStatePos}{NAME}$ratioTag',sizeName:'$pepSizeStrg'});\n";
		$ratioPos2index{$ratioPos}=$ratioIdx;
		$ratioIdx++;
	}
	my (%dispProteins,%dispIsoforms);
	foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
		my $count=0;
		my $ratioIdx=$ratioPos2index{$ratioPos};
		foreach my $modProtID (keys %{$refQuantifValues}) {
			my ($protID,$modStrg)=($modProtID=~/^(\d+)(.*)/); $modStrg='' unless $modStrg;
			#next if (!$quantifValues{$modProtID}{$pvalueCode} || !$quantifValues{$modProtID}{$pvalueCode}{$ratioPos} || !$quantifValues{$modProtID}{'RATIO'} || !$quantifValues{$modProtID}{'RATIO'}{$ratioPos});
			next if (!$refQuantifValues->{$modProtID}{'RATIO'} || !$refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos});
			next if $refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos} < $dispNumPep;
			#next if ($dispConfInt && defined $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos} && 1.96*$refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos} > $dispConfInt);
			if ($dispCV && defined $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}) {
				if ($ratioType eq 'Ratio') {
					next if abs(log($refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos})/log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos})) > $dispCVFrac;
				}
				else {
					next if $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}/abs(log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos})/$log2) > $dispCVFrac;
				}
			}
			if ($count==0) {
				print "\tVP.addDataAsString($ratioIdx,'";
			}
			else {print ";";}
			$count++;
#$modStrg=~s/,/./g; # for compatibility with volcano-plot
#(my $usedProtID=$modProtID)=~s/,/./g; # for compatibility with volcano-plot
			print "$refProteinInfo->{$protID}[0]$modStrg,$modProtID,";
			#>Set missing p-values=1!!!
			if (!$refQuantifValues->{$modProtID}{$pvalueCode} || !defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}) {
				$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}=1;
			}
			elsif ($refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}==0) {
				$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}=$minPvalue;
			}
			#>Set +/- infinite ratios
			if ($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos} == 1000 || $refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos} > $MAX_INF_RATIO) {
				$refProteinInfo->{$protID}[4]=1000 unless $refProteinInfo->{$protID}[4]; # set prot length
				printf "+,%.2f",100*$refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos}/$refProteinInfo->{$protID}[4]; # num pep/100 res
			}
			elsif ($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos} == 0.001 || $refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos} < $MIN_INF_RATIO) {
				$refProteinInfo->{$protID}[4]=1000 unless $refProteinInfo->{$protID}[4]; # set prot length
				printf "-,%.2f",100*$refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos}/$refProteinInfo->{$protID}[4]; # num pep/100 res
			}
			else { # normal ratios within [$MIN_INF_RATIO-$MAX_INF_RATIO] boundaries
				#$quantifValues{$modProtID}{$pvalueCode}{$ratioPos}=1 if (!$quantifValues{$modProtID}{$pvalueCode} || !$quantifValues{$modProtID}{$pvalueCode}{$ratioPos});
				#if ($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos} >= $MAX_INF_RATIO) {$refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos}='+';}
				#elsif ($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos} <= $MIN_INF_RATIO) {$refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos}='-';}
				printf "$refQuantifValues->{$modProtID}{RATIO}{$ratioPos},%.2e",$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos};
			}
			print ",$refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos}";
			if ($count==100) {
				print "');\n";
				$count=0;
			}
			$dispProteins{$protID}=1;
			$dispIsoforms{$modProtID}=1 if $modStrg;
		}
		if ($count > 0) {
			print "');\n";
		}
	}
	my $numProtDispStrg=(scalar keys %dispProteins).' proteins';
	$numProtDispStrg.=' / '.(scalar keys %dispIsoforms).' isoforms' if scalar keys %dispIsoforms;
	$numProtDispStrg.=' displayed';
	print qq
|	VP.draw();
}
</SCRIPT>
</CENTER>
<TABLE cellspacing=0 cellpadding=0>
<TR><TD><FONT class="title3">$numProtDispStrg</FONT></TD></TR>
<TR><TD><DIV id="mainGraphDIV" style="float:left;background-color:#FFFFFF;border:solid 3px $darkColor;padding:2px;"></DIV></TD></TR>
</TABLE>
<DIV id="protListDIV" style="position:absolute;height:535;overflow:auto"></DIV>
|;
}

sub displayIntensityPlot {
	my ($refLabelingInfo,$refStateInfo,$refQuantifValues,$refProteinInfo)=@_; # globals: $dispMeasure,$numPepCode,%quantifParamInfo
	my $measureName=$quantifParamInfo{$dispMeasure}[1];
	my $numPepName=$quantifParamInfo{$numPepCode}[1];
	$numPepName.=' peptides' if $numPepCode ne 'PEPTIDES';
	my $isLogged=1;
	my $axisYtitle=($isLogged)? (1,"Log10($measureName)") : (0,$measureName);
	my $log10=log(10);
#print '**',join(':',sort{$a<=>$b} keys %dispStates),"**<BR>\n";
	print qq
|</CENTER>
<SCRIPT type="text/javascript">
function myPointLabel(dp,type) {
	var infoStrg=dp.label;
	if (type != 'min') {
		if (dp.dataSet.params.chart.dataSets.length > 1) {infoStrg+='\\nState: '+dp.dataSet.params.name;}
		var intensity=($isLogged==1)? (1*(Math.pow(10,dp.getY())).toFixed(1)) : dp.getY();
		infoStrg+='\\n$measureName='+intensity+'\\nrank='+dp.getX()+'\\n$numPepName='+dp.size;
	}
	return infoStrg;
}
var mainChart,IP;
window.onload=function() {
	ajaxUpdateRestrict();
	IP=new genericPlot({div:'mainGraphDIV',width:600,height:450,
						axisX:{title:'Rank'}, // ,zoomable:true
						axisY:{title:'$axisYtitle'}, // ,zoomable:true
						axisClosure:true,
						//hideDataSets:true, // hide dataset label if only 1 displayed!
						zoomable:true,
						pointLabel:myPointLabel,
						pointOnclick:proteinLink,
						pointOnList:ajaxListSelectedProteins,
						allowHighlight:true,
						updateHighlight:{callback:updateUsedColors,editable:false},
						exportAsImage:['Export as image','IntensityPlot','./exportSVG.cgi']
					});
	mainChart=IP; // needed for highlighting
|;
	my $dsetIdx=0;
	foreach my $statePos (sort{$a<=>$b} keys %dispStates) {
		print "\tIP.addDataSet($dsetIdx,{name:'$refStateInfo->{$statePos}{NAME}'});\n";
		$dsetIdx++;
	}
	my $protCount;
	$dsetIdx=0;
	my %dispProteins;
	foreach my $statePos (sort{$a<=>$b} keys %dispStates) {
		my $count=0;
		$protCount=0;
		foreach my $protID (sort{&sortCheckDefined($dispMeasure,$refQuantifValues,$statePos,1) || &sortCheckDefined($numPepCode,$refQuantifValues,$statePos,1)} keys %{$refQuantifValues}) {
			next if (!$refQuantifValues->{$protID}{$dispMeasure} || !$refQuantifValues->{$protID}{$dispMeasure}{$statePos});
			next if (!$refQuantifValues->{$protID}{$numPepCode} || !$refQuantifValues->{$protID}{$numPepCode}{$statePos} || $refQuantifValues->{$protID}{$numPepCode}{$statePos} < $dispNumPep); # Not in DB if ==0 (eg. UNIQUE_PEP)
			#next if ($dispConfInt && defined $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos} && 1.96*$refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos} > $dispConfInt);
			if ($count==0) {
				print "\tIP.addDataAsString($dsetIdx,'";
			}
			else {print ";";}
			$count++;
			$protCount++;
			$refProteinInfo->{$protID}[0]=~s/[,;']/\./g; # Clean MaxQuant crappy contaminant identifiers
			print "$refProteinInfo->{$protID}[0],$protID,$protCount,";
			if ($isLogged) {print log($refQuantifValues->{$protID}{$dispMeasure}{$statePos})/$log10;} else {print $refQuantifValues->{$protID}{$dispMeasure}{$statePos};}
			print ','.$refQuantifValues->{$protID}{$numPepCode}{$statePos};
			if ($count==100) {
				print "');\n";
				$count=0;
			}
			$dispProteins{$protID}=1;
		}
		if ($count > 0) {
			print "');\n";
		}
		$dsetIdx++;
	}
	my $numProtDisp=scalar keys %dispProteins;
	print qq
|	IP.draw();
}
</SCRIPT>
<TABLE cellspacing=0 cellpadding=0>
<TR><TD><FONT class="title3">$numProtDisp proteins displayed</FONT></TR></TD>
<TR><TD><DIV id="mainGraphDIV" style="float:left;background-color:#FFFFFF;border:solid 3px $darkColor;padding:2px;"></DIV></TR></TD>
</TABLE>
<DIV id="protListDIV" style="position:absolute;height:535;overflow:auto"></DIV>
|;
}

sub displayRatioPlot {
	my ($refLabelingInfo,$refStateInfo,$refQuantifValues,$refProteinInfo)=@_; # globals: $dispMeasure,$numPepCode,%quantifParamInfo
	my $numPepName=$quantifParamInfo{$numPepCode}[1];
	$numPepName.=' peptides' if $numPepCode ne 'PEPTIDES';
	my $revFoldChange=1/$dispFoldChange;
	my $log2=log(2);
	print qq
|</CENTER>
<SCRIPT type="text/javascript">
function myPointLabel(dp,type) {
	var infoStrg=dp.label;
	if (type != 'min') {
		if (dp.dataSet.params.chart.dataSets.length > 1) {infoStrg+='\\nState: '+dp.dataSet.params.name;}
		var ratio=(dp.getX() >= 0)? 1*((Math.pow(2,dp.getX())).toFixed(2)) : '1/'+(1*((1/(Math.pow(2,dp.getX()))).toFixed(2)));
		infoStrg+='\\nFC='+ratio+'\\nCoef. var.='+dp.getY()+'%\\n$numPepName='+dp.size;
	}
	return infoStrg;
}
function convertThreshold(axis,value) {
	if (axis=='X') {value=Math.log2(value);}
	return value;
}
var mainChart,RP;
window.onload=function() {
	ajaxUpdateRestrict();
	RP=new genericPlot({div:'mainGraphDIV',width:600,height:450,
						axisX:{title:'Log2(fold change)'}, // ,zoomable:true
						axisY:{title:'Coefficient variation (%)'}, // ,zoomable:true
						axisClosure:true,
						//hideDataSets:true,
						zoomable:true,
						editThreshold:true,
						convertValue:convertThreshold, // threshold lines only
						pointLabel:myPointLabel,
						pointOnclick:proteinLink,
						pointOnList:ajaxListSelectedProteins,
						allowHighlight:true,
						updateHighlight:{callback:updateUsedColors,editable:false},
						exportAsImage:['Export as image','RatioPlot','./exportSVG.cgi']
					});
	mainChart=RP; // needed for highlighting
	RP.addThreshold({axis:'Y',label:'CV thresold',value:$dispCV,color:'#FF0000',editable:true});
	RP.addThreshold({axis:'X',label:'min. fold change',value:$revFoldChange,color:'#00A000',keepAbove1:true,editable:true});
	RP.addThreshold({axis:'X',label:'max. fold change',value:$dispFoldChange,color:'#00A000',keepAbove1:true,editable:true});
|;
	my $ratioIdx=0;
	my %ratioPos2index;
	foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
		my ($testStateID,$refStateID)=($ratioPos > 0)? split(/\//,$refLabelingInfo->{'RATIOS'}[$ratioPos-1]) : reverse(split(/\//,$refLabelingInfo->{'RATIOS'}[abs($ratioPos)-1]));
		my $testStatePos=$condToState{$testStateID};
		my $refStatePos=$condToState{$refStateID};
		print "\tRP.addDataSet($ratioIdx,{name:'$refStateInfo->{$testStatePos}{NAME}/$refStateInfo->{$refStatePos}{NAME}'});\n";
		$ratioPos2index{$ratioPos}=$ratioIdx;
		$ratioIdx++;
	}
	my %dispProteins;
	foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
		my $count=0;
		foreach my $modProtID (keys %{$refQuantifValues}) {
			my ($protID,$modStrg)=($modProtID=~/^(\d+)(.*)/); $modStrg='' unless $modStrg;
next if (!$refQuantifValues->{$modProtID}{$ratioParamCode} || !$refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos});
			next if $refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos} < $dispNumPep;
#next if ($dispCV && defined $refQuantifValues->{$modProtID}{'RATIO_VAR'}{$ratioPos} && $refQuantifValues->{$modProtID}{'RATIO_VAR'}{$ratioPos} > $dispCV);
			if ($count==0) {
				print "\tRP.addDataAsString($ratioPos2index{$ratioPos},'";
			}
			else {print ";";}
			$count++;
			$refProteinInfo->{$protID}[0]=~s/[,;']/\./g; # Clean MaxQuant crappy contaminant identifiers
			my $coefVar=$refQuantifValues->{$modProtID}{RATIO_VAR}{$ratioPos} || 0; # not defined for < 3 peptides
			print "$refProteinInfo->{$protID}[0]$modStrg,$modProtID,",log($refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos})/$log2,",$coefVar,$refQuantifValues->{$protID}{$numPepCode}{$ratioPos}";
			if ($count==100) {
				print "');\n";
				$count=0;
			}
			$dispProteins{$protID}=1;
		}
		if ($count > 0) {
			print "');\n";
		}
	}
	my $numProtDisp=scalar keys %dispProteins;
	print qq
|	RP.draw();
}
</SCRIPT>
<TABLE cellspacing=0 cellpadding=0>
<TR><TD><FONT class="title3">$numProtDisp proteins displayed</FONT></TD></TR>
<TR><TD><DIV id="mainGraphDIV" style="float:left;background-color:#FFFFFF;border:solid 3px $darkColor;padding:2px;"></DIV></TD></TR>
</TABLE>
<DIV id="protListDIV" style="position:absolute;height:535;overflow:auto"></DIV>
|;
}

sub printProteinRatioList {
	my ($labelType,$ratioType,$refLabelingInfo,$refStateInfo,$refQuantifValues,$refProteinInfo,$pvalueCode)=@_; # $labelType not global for all calls
	my ($minFoldChange,$maxFoldChange,$maxPvalue)=(0.5,2,0.05); # default values for arrow flag
	my $invDispFoldChange=1/$dispFoldChange;
	my $quantifSoftware=($refLabelingInfo->{'SOFTWARE'})? $refLabelingInfo->{'SOFTWARE'}[0] : 'myProMS';
	my $numTotQuantifItems=scalar keys %{$refQuantifValues};
	my $numTotProteins=scalar keys %{$refProteinInfo};
	my $protString=($quantifModID)? "&nbsp;$numTotQuantifItems isoforms&nbsp;<BR>&nbsp;$numTotProteins proteins&nbsp;" : "&nbsp;$numTotProteins proteins&nbsp;";
	my $peptideTitle=($labelType eq 'SILAC')? 'Pept. sets' : 'Peptides';
	my $disabSave=($projectAccess eq 'guest')? ' disabled' : '';
	my $clearButtonStrg=($view eq 'list')? '' : '<INPUT type="button" value="Clear list" style="font-weight:bold" onclick="document.getElementById(\'protListDIV\').style.display=\'none\'">';
	my ($ratioColspan,$cvColStrg)=($quantifType eq 'TNPQ')? (3,'') : (4,'<TH class="rbBorder" nowrap>&nbsp;CV (%)&nbsp;</TH>'); # ratio/arrow/p-value/(+/-stdDev)
	#my $ratioType=($refLabelingInfo->{'RATIO_TYPE'})? $refLabelingInfo->{'RATIO_TYPE'}[0] : 'Ratio';
	my %supRatios;
	if ($ratioType eq 'SuperRatio') {
		foreach my $ratioPos (keys %dispRatios) {
			$supRatios{$ratioPos}=1 if $refLabelingInfo->{'RATIOS'}[abs($ratioPos)-1]=~/%/;
		}
	}
	my $numColumns=scalar (keys %dispRatios) * $ratioColspan + 5;
	print qq
|<FORM name="protForm" method="POST">
<TABLE border=0 cellspacing=0 cellpadding=2>
<TR><TD colspan=$numColumns><DIV id="saveProtDIV" style="display:none;"></DIV></TD></TR>
<TR><TD colspan=$numColumns>$clearButtonStrg
<INPUT type="button" value="Check all" onclick="checkAllProteins(true)"/>
<INPUT type="button" value="Uncheck all" onclick="checkAllProteins(false)"/>
<INPUT type="button" id="saveFormBUTTON" value="Save proteins..." onclick="ajaxManageSaveProteins('getThemes','PROT')"$disabSave/>
|;
	print "<INPUT type=\"button\" id=\"saveSiteFormBUTTON\" value=\"Save sites...\" onclick=\"ajaxManageSaveProteins('getThemes','SITE',$quantifModID)\"$disabSave/>\n" if $quantifModID;
	#print "<INPUT type=\"button\" value=\"Export data\" onclick=\"exportProteins()\"/>\n" if $action ne 'ajaxListProt';
print "<INPUT type=\"button\" value=\"Export data\" onclick=\"exportProteins()\"/>\n" if ($action ne 'ajaxListProt'); # && $algoVersion != 3); # !!!TEMP!!!
#print "<INPUT type=\"button\" value=\"Export isoform fragments\" onclick=\"exportIsoforms()\"/>\n" if ($quantifModID && $userID eq 'sliva'); # !!!!! TEMP !!!!
	print qq
|</TD></TR>
<TR bgcolor="$darkColor">
<TH class="rbBorder" align=left rowspan=2><DIV id="protCountDIV">$protString</DIV></TH><TH class="rbBorder" rowspan=2>Gene</TH>
|; #&nbsp;$proteinTitle&nbsp;
	foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
		my ($testStatePos,$refStatePos)=($ratioPos > 0)? split(/\//,$refLabelingInfo->{'RATIOS'}[$ratioPos-1]) : reverse(split(/\//,$refLabelingInfo->{'RATIOS'}[abs($ratioPos)-1]));
		my $ratioTag='';
		if ($designID) {
			if ($supRatios{$ratioPos}) {
				$ratioTag=encode_utf8('°');
				$testStatePos=~s/%\d+//;
				$refStatePos=~s/%\d+//;
			}
			$testStatePos=$condToState{$testStatePos};
			$refStatePos=$condToState{$refStatePos};
		}
		print "<TH class=\"rbBorder\" colspan=$ratioColspan nowrap>&nbsp;<A id=\"listRatio_$ratioPos\" href=\"javascript:displayRatioStat('listRatio_$ratioPos','$ratioCode{$ratioPos}','$refStateInfo->{$testStatePos}{NAME}$ratioTag/$refStateInfo->{$refStatePos}{NAME}$ratioTag')\">$refStateInfo->{$testStatePos}{NAME}$ratioTag/$refStateInfo->{$refStatePos}{NAME}$ratioTag</A>&nbsp;</TH>\n";
	}
	if ($ratioType eq 'Ratio') {
		my $popupInfo=($numPepCode eq 'DIST_PEP_USED')? " onmouseover=\"popup('<B>distinct/used/total</B>')\" onmouseout=\"popout()\"" : " onmouseover=\"popup('<B>used/total</B>')\" onmouseout=\"popout()\"" ;
		print "<TH class=\"rbBorder\" rowspan=2 nowrap$popupInfo>&nbsp;$peptideTitle&nbsp;<BR>used</TH>\n";
	}
	else { # Super/Simple Ratio
		print "<TH class=\"rbBorder\" rowspan=2 nowrap>All<BR>&nbsp;peptides&nbsp;</TH>\n";
	}
	print qq
|<TH class="rbBorder" rowspan=2 nowrap>&nbsp;MW<SMALL> kDa</SMALL>&nbsp;</TH><TH class="bBorder" rowspan=2 width=700 nowrap>&nbsp;Description - Species&nbsp;</TH>
</TR>
<TR bgcolor="$darkColor">
|;
	my $firstRatioUsed=0;
	foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
		$firstRatioUsed=$ratioPos unless $firstRatioUsed;
		print "<TH class=\"rbBorder\" colspan=2 nowrap>&nbsp;Ratio&nbsp;</TH><TH class=\"rbBorder\" nowrap>&nbsp;p-value&nbsp;</TH>";
		if ($ratioType=~/S\w+Ratio/) {
			my $popupInfo=($numPepCode eq 'DIST_PEP_USED')? " onmouseover=\"popup('<B>distinct/used</B>')\" onmouseout=\"popout()\"" : " onmouseover=\"popup('<B>used</B>')\" onmouseout=\"popout()\"" ;
			print "<TH class=\"rbBorder\" nowrap$popupInfo>&nbsp;$peptideTitle&nbsp;</TH>";
		}
		else {print $cvColStrg;}
	}
	print "</TR>\n";

	my $bgColor=$lightColor;
	my $ajaxPepAction=($quantifSoftware=~/SWATH|DIA/)? 'ajaxPepSwath' : ($quantifSoftware eq 'MaxQuant')? 'ajaxPepMaxQuant' : 'ajaxPepRatio';
	my $numDispItems=0;
	my %dispProteins;
	foreach my $modProtID (sort{&sortProteins($dispSort,$refProteinInfo,$refQuantifValues,$firstRatioUsed)} keys %{$refQuantifValues}) { #%{$refProteinInfo}
#last if $numDispItems >= 5;
		my ($protID,$modResStrg)=($modProtID=~/^(\d+)(.*)/); # eg. 12345-Y40,T152 -> 12345, -Y40,T152
		$modResStrg='' unless $modResStrg;
		#>Peptide filter
		my $refUsedRatioPos=0;
		if ($view eq 'list') {
			my $okPeptide=0;
			foreach my $ratioPos (keys %dispRatios) {
				next if (!$refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos} || $refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos} < $dispNumPep); # may not be defined for +/-inf ratios
				$okPeptide=1; # at least 1 ratio must pass the filter
				$refUsedRatioPos=$ratioPos;
				last;
			}
			next unless $okPeptide;
		}
		my $mass=(!$refProteinInfo->{$protID}[1])? '-' : $refProteinInfo->{$protID}[1];
		my $anaList = $refProteinInfo->{$protID}[-1];
		my $quantifDataStrg='';
		my $okFilters=($view eq 'list')? 0 : 1; # 1 in case ajax call from graphical view
		my %protCV;
		foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
			$refUsedRatioPos=$ratioPos unless $refUsedRatioPos;
			if (!$refQuantifValues->{$modProtID} || !$refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos}) {
				$quantifDataStrg.="<TH colspan=2>-</TH><TH>-</TH>";
				$quantifDataStrg.="<TH>-</TH>" if $quantifType ne 'TNPQ'; # Std dev or SuperRatio prim ratio
				next;
			}
			# Ratio
			my $ratio=$refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos};
			my ($absRatio,$arrowDir)=($ratio >= 1)? ($ratio,'up_') : (1/$ratio,'down_');
			my $ratioStrg=($absRatio == 1000)? '&infin;' : sprintf "%.2f",$absRatio; # '∞'
			$ratioStrg='1/'.$ratioStrg if ($ratio < 1 && $ratioStrg ne '1.00');
			# Arrow flag
			my $arrowStrg='';
			my $okRatio=($ratio <= $minFoldChange || $ratio >= $maxFoldChange)? 1 : 0;
			my $okPvalue=(defined $refQuantifValues->{$modProtID}{$pvalueCode} && defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} && $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} <= $maxPvalue)? 1 : 0;
			if ($okRatio || $okPvalue) {
				my $arrowColor=(!$okRatio || !$okPvalue)? 'gray' : ($ratio >= 1)? 'red' : 'green';
				$arrowStrg="<IMG class=\"LINK\" src=\"$promsPath{images}/$arrowDir$arrowColor.png\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'$ajaxPepAction')\">";
			}
			#if ($ratio <= $minFoldChange || $ratio >= $maxFoldChange) {
			#	my $arrowColor=(!$refQuantifValues->{$modProtID}{$pvalueCode} || !defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} || $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} > $maxPvalue)? 'gray' : ($ratio >= 1)? 'red' : 'green';
			#	$arrowStrg="<IMG class=\"LINK\" src=\"$promsPath{images}/$arrowDir$arrowColor.png\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'$ajaxPepAction')\">";
			#}
			$quantifDataStrg.="<TH nowrap align=right>&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'$ajaxPepAction')\"/>$ratioStrg</A></TH><TD>$arrowStrg&nbsp;</TD>";
			# p-value
			#my $pValueStrg=(!defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos})? '-' : ($refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}>=0.01)? sprintf "%.2f",$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} : sprintf "%.1e",$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos};

			if (defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}) {
				my $pValueStrg=($refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}>=0.01)? sprintf "%.2f",$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} : sprintf "%.1e",$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos};
				$quantifDataStrg.="<TH nowrap>&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'ajaxProtStat')\">$pValueStrg</A>&nbsp;</TH>";
			}
			else {$quantifDataStrg.="<TH nowrap>&nbsp;-&nbsp;</TH>";}
			# Std dev OR num Pep
			if ($ratioType eq 'Ratio') {
				if ($quantifType ne 'TNPQ') { # Std dev
					if ($refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}) {
						$protCV{$ratioPos}=($ratioType eq 'Ratio')? abs(log($refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos})/log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos})) : $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}/abs(log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos})/$log2);
						my $cvStrg=sprintf "%.1f",($protCV{$ratioPos}*100);
						$quantifDataStrg.="<TH nowrap>&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'ajaxProtStat')\">$cvStrg</A>&nbsp;</TH>";
					}
					else {$quantifDataStrg.="<TH nowrap>&nbsp;-&nbsp;</TH>";}
				}
			}
			else {
				$quantifDataStrg.="<TH nowrap>&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'$ajaxPepAction')\"/>";
				$quantifDataStrg.=($numPepCode eq 'DIST_PEP_USED')? "$refQuantifValues->{$modProtID}{DIST_PEP_USED}{$ratioPos}/$refQuantifValues->{$modProtID}{NUM_PEP_USED}{$ratioPos}" : $refQuantifValues->{$modProtID}{NUM_PEP_USED}{$ratioPos};
				$quantifDataStrg.="</A>&nbsp;</TH>";
			}
			unless ($okFilters) { # no ratioPos has yet passed filters
				if ($foldChangeType eq 'abs') {
					$okFilters=(($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} < 1 && $refQuantifValues->{$modProtID}{RATIO}{$ratioPos} <= $invDispFoldChange) || ($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos} >= 1 && $refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos} >= $dispFoldChange))? 1 : 0;
				}
				elsif ($foldChangeType eq 'up') {
					$okFilters=($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} >= 1 && $refQuantifValues->{$modProtID}{RATIO}{$ratioPos} >= $dispFoldChange)? 1 : 0;
				}
				else { # down
					$okFilters=($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} <= 1 && $refQuantifValues->{$modProtID}{RATIO}{$ratioPos} <= $invDispFoldChange)? 1 : 0;
				}
				$okFilters=($okFilters && ($dispPvalue>=1 || (defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} && $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} <= $dispPvalue)))? 1 : 0;
				#if ($okFilters && $dispConfInt && defined $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}) { # undefined always passes filter!
				#	$okFilters=(1.96*$refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos} <= $dispConfInt)? 1 : 0;
				#}
				if ($okFilters && $dispCV && defined $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}) { # undefined always passes filter!
					#$okFilters=0 if ($refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}/abs(log($refQuantifValues->{$modProtID}{RATIO}{$ratioPos}/$log2)) > $dispCVFrac);
					$okFilters=0 if $protCV{$ratioPos} > $dispCVFrac;
				}
			}
		}
#print "$modProtID<BR>\n" unless $okFilters;
		next unless $okFilters;
		$numDispItems++;
		$dispProteins{$protID}=1;
		print "<TR bgcolor=\"$bgColor\" valign=top><TD class=\"TH\" nowrap><INPUT type=\"checkbox\" name=\"chkProt\" value=\"$modProtID\"/><A href=\"javascript:sequenceView($protID,'$anaList')\">$refProteinInfo->{$protID}[0]$modResStrg</A>&nbsp;</TD>\n";
		# Gene(s)
		if (scalar @{$refProteinInfo->{$protID}[5]} > 1) {
			print "<TH>&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-",join('<BR>&nbsp;&nbsp;-',@{$refProteinInfo->{$protID}[5]}[1..$#{$refProteinInfo->{$protID}[5]}]),"</B>')\" onmouseout=\"popout()\">",$refProteinInfo->{$protID}[5][0],"</A>&nbsp;</TH>";
		}
		elsif ($refProteinInfo->{$protID}[5][0]) {print '<TH>&nbsp;',$refProteinInfo->{$protID}[5][0],'&nbsp;</TH>';}
		else {print '<TH>-</TH>';} # no gene mapped
		# quantif data
		print $quantifDataStrg;
		if ($refQuantifValues->{$modProtID} && $quantifType ne 'MQ') {# For SILAC quantifications made with MaxQuant, only redundant peptides are provided in proteinGroups.txt (column called 'Ratio H/L count')
			#my $distPepStrg=($refQuantifValues->{$modProtID}{'DIST_PEP_USED'}==$refQuantifValues->{$modProtID}{'NUM_PEP_USED'})? '' : $refQuantifValues->{$modProtID}{'DIST_PEP_USED'}.'/';
			#my $distPepUsedStrg=($refQuantifValues->{$modProtID}{DIST_PEP_USED})? $refQuantifValues->{$modProtID}{DIST_PEP_USED} : '-'; #
			my $distPepUsedStrg=($ratioType=~/S\w+Ratio/)? $refQuantifValues->{$modProtID}{NUM_PEP_TOTAL}{$refUsedRatioPos} : ($numPepCode eq 'DIST_PEP_USED')? "$refQuantifValues->{$modProtID}{DIST_PEP_USED}{$refUsedRatioPos}/$refQuantifValues->{$modProtID}{NUM_PEP_USED}{$refUsedRatioPos}/$refQuantifValues->{$modProtID}{NUM_PEP_TOTAL}{$refUsedRatioPos}": "$refQuantifValues->{$modProtID}{NUM_PEP_USED}{$refUsedRatioPos}/$refQuantifValues->{$modProtID}{NUM_PEP_TOTAL}{$refUsedRatioPos}"; #
			print "<TH>&nbsp;$distPepUsedStrg&nbsp;</TH>";
		}
		else {print "<TH>-</TH>"}
		# mw, des & species
		print "<TD class=\"TH\" align=right>$mass&nbsp;</TD><TD>$refProteinInfo->{$protID}[2] <U><I>$refProteinInfo->{$protID}[3]</I></U></TD>\n";
		print "</TR>\n";
		$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
	}
	$protString=($quantifModID)? "&nbsp;$numDispItems/$numTotQuantifItems isoforms<BR>&nbsp;".(scalar keys %dispProteins)."/$numTotProteins proteins&nbsp;" : "&nbsp;$numDispItems/$numTotProteins proteins&nbsp;";
	print qq
|<TR><TD colspan=$numColumns><B>End of list.</B>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<IMG src="$promsPath{images}/up_red.png"/><IMG src="$promsPath{images}/down_green.png"/>: Fold change &ge; 2 <B>and</B> p-value &le; 0.05&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<IMG src="$promsPath{images}/up_gray.png"/><IMG src="$promsPath{images}/down_gray.png"/>: Fold change &ge; 2 <B>or</B> p-value &le; 0.05</TD></TR>
</TABLE>
</FORM>
<SCRIPT LANGUAGE="javascript"><!-- Not executed if called from AJAX -->
document.getElementById('protCountDIV').innerHTML='$protString';
</SCRIPT>
|;
}


sub printMaxQuantProteinList {
	my ($labelType,$refLabelingInfo,$refStateInfo,$refQuantifValues,$refProteinInfo)=@_; # $labelType  not global for all calls
	my $measureName=$quantifParamInfo{$dispMeasure}[1]; # only for no-ratio
	my $invDispFoldChange=1/$dispFoldChange;
	my $numTotQuantifItems=scalar keys %{$refQuantifValues};
	my $numTotProteins=scalar keys %{$refProteinInfo};
	my $protString=($quantifModID)? "&nbsp;$numTotQuantifItems isoforms&nbsp;<BR>&nbsp;$numTotProteins proteins&nbsp;" : "&nbsp;$numTotProteins proteins&nbsp;";
	my $disabSave=($projectAccess eq 'guest')? ' disabled' : '';
	my $clearButtonStrg=($view eq 'list')? '' : '<INPUT type="button" value="Clear list" style="font-weight:bold" onclick="document.getElementById(\'protListDIV\').style.display=\'none\'">';
	my $numColumns=($ratioType eq 'None')? scalar (keys %dispStates) + 6 : scalar (keys %dispRatios) * 3 + 6;
	print qq
|<FORM name="protForm" method="POST">
<TABLE border=0 cellspacing=0 cellpadding=2>
<TR><TD colspan=$numColumns><DIV id="saveProtDIV" style="display:none;"></DIV></TD></TR>
<TR><TD colspan=$numColumns>$clearButtonStrg
<INPUT type="button" value="Check all" onclick="checkAllProteins(true)"/>
<INPUT type="button" value="Uncheck all" onclick="checkAllProteins(false)"/>
<INPUT type="button" id="saveFormBUTTON" value="Save proteins..." onclick="ajaxManageSaveProteins('getThemes','PROT')"$disabSave/>
|;
	print "<INPUT type=\"button\" id=\"saveSiteFormBUTTON\" value=\"Save sites...\" onclick=\"ajaxManageSaveProteins('getThemes','SITE',$quantifModID)\"$disabSave/>\n" if $quantifModID;
	print "<INPUT type=\"button\" value=\"Export data\" onclick=\"exportProteins()\" disabled/>\n" unless $action eq 'ajaxListProt';
	print qq
|</TD></TR>
<TR bgcolor="$darkColor">
<TH class="rbBorder" align=left rowspan=2><DIV id="protCountDIV">$protString</DIV></TH><TH class="rbBorder" rowspan=2>Gene</TH>
|;
	if ($ratioType eq 'None') {
		my $numStates=scalar keys %dispStates;
		my $stateStrg=($numStates > 1)? 'States' : 'State';
		print "<TH class=\"rbBorder\" colspan=$numStates>&nbsp$stateStrg&nbsp;</TH>";
	}
	else {
		foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
			my ($testStateID,$refStateID)=($ratioPos > 0)? split(/\//,$refLabelingInfo->{'RATIOS'}[$ratioPos-1]) : reverse(split(/\//,$refLabelingInfo->{'RATIOS'}[abs($ratioPos)-1]));
			print "<TH class=\"rbBorder\" colspan=2 nowrap>&nbsp;$refStateInfo->{$condToState{$testStateID}}{NAME}/$refStateInfo->{$condToState{$refStateID}}{NAME}&nbsp;</TH>\n";
		}
	}
	print qq
|<TH class="rbBorder" colspan=3 nowrap>&nbsp;Peptides&nbsp;</TH>
<TH class="rbBorder" rowspan=2 nowrap>&nbsp;MW<SMALL> kDa</SMALL>&nbsp;</TH><TH class="bBorder" rowspan=2 width=700 nowrap>&nbsp;Description - Species&nbsp;</TH>
</TR>
<TR bgcolor="$darkColor">
|;
	my $firstSetUsed;
	if ($ratioType eq 'None') {
		foreach my $statePos (sort{$a<=>$b} keys %dispStates) {
			$firstSetUsed=$statePos unless $firstSetUsed;
			print "<TH class=\"rbBorder\" nowrap>&nbsp;$refStateInfo->{$statePos}{NAME}&nbsp;</TH>";
		}
	}
	else {
		my $ratioTitle=($foldChangeRaw)? 'Raw ratio' : 'Norm. ratio';
		foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
			$firstSetUsed=$ratioPos unless $firstSetUsed;
			print "<TH class=\"rbBorder\" nowrap>&nbsp;$ratioTitle&nbsp;</TH><TH class=\"rbBorder\" nowrap>&nbsp;Ratio count&nbsp;</TH>\n";
		}
	}
	print qq
|<TH class="rbBorder"nowrap>&nbsp;Unique&nbsp;</TH><TH class="rbBorder" nowrap>&nbsp;Razor+uniq.&nbsp;</TH><TH class="rbBorder" nowrap>&nbsp;All&nbsp;</TH>
</TR>
|;
	my $refDispSets=($ratioType eq 'None')? \%dispStates : \%dispRatios;
	my $bgColor=$lightColor;
	my $numDispItems=0;
	my %dispProteins;
	foreach my $modProtID (sort{&sortProteins($dispSort,$refProteinInfo,$refQuantifValues,$firstSetUsed)} keys %{$refQuantifValues}) { #%{$refProteinInfo}
		my ($protID,$modResStrg)=($modProtID=~/^(\d+)(.*)/); # eg. 12345-Y40,T152 -> 12345, -Y40,T152
		$modResStrg='' unless $modResStrg;
		#>Peptide filter
		if ($view eq 'list') { # can also be called from graph view by ajax on subset of proteins => no filtering
			my $okPeptide=0;
			foreach my $setPos (keys %{$refDispSets}) {
				next if (!$refQuantifValues->{$modProtID}{$numPepCode} || !$refQuantifValues->{$modProtID}{$numPepCode}{$setPos} || $refQuantifValues->{$modProtID}{$numPepCode}{$setPos} < $dispNumPep); # may not be defined for +/-inf ratios
				$okPeptide=1; # at least 1 ratio must pass the filter
				last;
			}
			next unless $okPeptide;
		}
		my $mass=(!$refProteinInfo->{$protID}[1])? '-' : $refProteinInfo->{$protID}[1];
		my $anaList = $refProteinInfo->{$protID}[-1];
		my $quantifDataStrg='';
		if ($ratioType eq 'None') {
			foreach my $statePos (sort{$a<=>$b} keys %dispStates) {
				$quantifDataStrg.=($refQuantifValues->{$modProtID} && $refQuantifValues->{$modProtID}{$dispMeasure}{$statePos})? "<TH align=right>&nbsp;$refQuantifValues->{$modProtID}{$dispMeasure}{$statePos}&nbsp;</TH>" : '<TH>-</TH>';
			}
		}
		else {
			my $okFilters=($view eq 'list')? 0 : 1; # 1 in case ajax call from graphical view
			foreach my $ratioPos (sort{abs($a)<=>abs($b)} keys %dispRatios) {
				if (!$refQuantifValues->{$modProtID} || !$refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos}) {
					$quantifDataStrg.="<TH>-</TH><TH>-</TH>";
					next;
				}
				unless ($okFilters) { # no ratioPos has yet passed filters
					if ($foldChangeType eq 'abs') {
						$okFilters=(($refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} < 1 && $refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} <= $invDispFoldChange) || ($refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} >= 1 && $refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} >= $dispFoldChange))? 1 : 0;
					}
					elsif ($foldChangeType eq 'up') {
						$okFilters=($refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} >= 1 && $refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} >= $dispFoldChange)? 1 : 0;
					}
					else { # down
						$okFilters=($refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} <= 1 && $refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} <= $invDispFoldChange)? 1 : 0;
					}
					if ($okFilters && $dispCV && $refQuantifValues->{$modProtID}{'RATIO_VAR'} &&  defined $refQuantifValues->{$modProtID}{'RATIO_VAR'}{$ratioPos}) { # undefined always passes filter!
						$okFilters=0 if $refQuantifValues->{$modProtID}{'RATIO_VAR'}{$ratioPos} > $dispCV;
					}
					if ($okFilters && $refQuantifValues->{$modProtID}{$numPepCode} && $refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos}) {
						$okFilters=0 if $refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos} < $dispNumPep;
					}
				}
				my $ratioStrg=($refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} >= 1)? sprintf "%.2f",$refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos} : '1/'.sprintf "%.2f",1/$refQuantifValues->{$modProtID}{$ratioParamCode}{$ratioPos};
				#$quantifDataStrg.="<TH nowrap>&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'ajaxProtStat')\">$ratioStrg</A>&nbsp;</TH><TH>$refQuantifValues->{$modProtID}{NUM_PEP_USED}{$ratioPos}</TH>";
				$quantifDataStrg.="<TH nowrap>&nbsp;$ratioStrg&nbsp;</TH><TH>$refQuantifValues->{$modProtID}{NUM_PEP_USED}{$ratioPos}</TH>";
			}
			next unless $okFilters;
		}
		$numDispItems++;
		$dispProteins{$protID}=1;
		my $class=($anaList)? 'TH' : 'TD'; # quantified MaxQuant proteins can be hidden in myProMS
		print "<TR bgcolor=\"$bgColor\" valign=top><TD class=\"$class\" nowrap><INPUT type=\"checkbox\" name=\"chkProt\" value=\"$modProtID\"/><A href=\"javascript:sequenceView($protID,'$anaList')\">$refProteinInfo->{$protID}[0]$modResStrg</A>&nbsp;</TD>\n";
		# Gene(s)
		if (scalar @{$refProteinInfo->{$protID}[5]} > 1) {
			print "<TH>&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-",join('<BR>&nbsp;&nbsp;-',@{$refProteinInfo->{$protID}[5]}[1..$#{$refProteinInfo->{$protID}[5]}]),"</B>')\" onmouseout=\"popout()\">",$refProteinInfo->{$protID}[5][0],"</A>&nbsp;</TH>";
		}
		elsif ($refProteinInfo->{$protID}[5][0]) {print '<TH>&nbsp;',$refProteinInfo->{$protID}[5][0],'&nbsp;</TH>';}
		else {print '<TH>-</TH>';} # no gene mapped
		# quantif data
		print $quantifDataStrg;
		# peptides (all states)
		my $uniqStrg=$refQuantifValues->{$modProtID}{'UNIQUE_PEP'}{$firstSetUsed} || '-';
		my $razUniStrg=$refQuantifValues->{$modProtID}{'RAZ_UNI_PEP'}{$firstSetUsed} || '-';
		my $peptideStrg=$refQuantifValues->{$modProtID}{'PEPTIDES'}{$firstSetUsed} || '-';
		print "<TH>$uniqStrg</TH><TH>$razUniStrg</TH><TH>$peptideStrg</TH>";
		# mw, des & species
		print "<TD class=\"TH\" align=right>$mass&nbsp;</TD><TD>$refProteinInfo->{$protID}[2] <U><I>$refProteinInfo->{$protID}[3]</I></U></TD>\n";
		print "</TR>\n";
		$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
	}
	$protString=($quantifModID)? "&nbsp;$numDispItems/$numTotQuantifItems isoforms<BR>&nbsp;".(scalar keys %dispProteins)."/$numTotProteins proteins&nbsp;" : "&nbsp;$numDispItems/$numTotProteins proteins&nbsp;";
	print qq
|</TABLE>
</FORM>
<SCRIPT LANGUAGE="javascript"><!-- Not executed if called from AJAX -->
document.getElementById('protCountDIV').innerHTML='$protString';
</SCRIPT>
|;
}


sub exportProteinList { # Only for design or no-design labeled quantifs
	my ($labelType,$ratioType,$refLabelingInfo,$refStateInfo,$refQuantifValues,$refProteinInfo,$pvalueCode)=@_; #,$dispStdDev
	my $invDispFoldChange=1/$dispFoldChange;
	my $quantifSoftware=($refLabelingInfo->{'SOFTWARE'})? $refLabelingInfo->{'SOFTWARE'}[0] : 'myProMS';
	my $algoVersion=0; # myProMS only
	if ($quantifSoftware eq 'myProMS') {
		$algoVersion=($refLabelingInfo->{'SOFTWARE'} && $refLabelingInfo->{'SOFTWARE'}[1])? $refLabelingInfo->{'SOFTWARE'}[1] : ($ratioType eq 'Ratio')? 1 : 2;
	}
	my $numTotQuantifItems=scalar keys %{$refQuantifValues};
	my $numTotProteins=scalar keys %{$refProteinInfo};
	my $peptideTitle=($labelType eq 'SILAC')? 'Pept. sets' : 'Peptides';

	####<Start printing>####
	my $worksheet2=$workbook->add_worksheet('Results');
	my $xlsRow=0;
	my $xlsCol=0;

	###<Headers>###
	#<Identifiers
	$worksheet2->set_column(0,0,30); # identifier col length
	$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow+1,$xlsCol,'Gene & Synonyms',$itemFormat{'mergeRowHeader'});
	#$worksheet2->merge_range($xlsRow,$xlsCol,$xlsRow+1,$xlsCol,"Proteins",$itemFormat{'mergeRowHeader'});
	#<Ratios
	#my $ratioType=($refLabelingInfo->{'RATIO_TYPE'})? $refLabelingInfo->{'RATIO_TYPE'}[0] : 'Ratio';
	my (%ratioColDelta); # %supRatios,
	my $curRatioPos=0;
	my $firstRatioUsed=0;
	foreach my $ratioData (@{$refLabelingInfo->{RATIOS}}) {
		$curRatioPos++;
		next unless $dispRatios{$curRatioPos};
		$firstRatioUsed=$curRatioPos unless $firstRatioUsed;
		my ($testStatePos,$refStatePos)=split(/\//,$ratioData);
		my $ratioTag='';
		if ($designID) {
			$testStatePos=~s/#//;
			$refStatePos=~s/#//;
			if ($ratioData=~/%/) {
				#$supRatios{$curRatioPos}=1;
				$ratioTag='°';
				$testStatePos=~s/%\d+//;
				$refStatePos=~s/%\d+//;
			}
			$testStatePos=$condToState{$testStatePos};
			$refStatePos=$condToState{$refStatePos};
		}
		$ratioColDelta{$curRatioPos}=($quantifType eq 'TNPQ')? 4 : ($quantifSoftware eq 'myProMS')? 5 : 5;
		$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow,$xlsCol+$ratioColDelta{$curRatioPos},$refStateInfo->{$testStatePos}{'NAME'}.$ratioTag.'/'.$refStateInfo->{$refStatePos}{'NAME'}.$ratioTag,$itemFormat{'mergeColHeader'});
		$worksheet2->write_string($xlsRow+1,$xlsCol,'Ratio',$itemFormat{'header'});
		$worksheet2->write_comment($xlsRow+1,$xlsCol,'Ratio of 1000 or 0.001 means protein was not detected in one condition.');
		$worksheet2->write_string($xlsRow+1,++$xlsCol,'Log2',$itemFormat{'header'});
		#if ($ratioType eq 'Ratio') {
		#	$worksheet2->write_string($xlsRow+1,++$xlsCol,'Lower bound.',$itemFormat{'header'});
		#	$worksheet2->write_comment($xlsRow+1,$xlsCol,"Lower boundary of ratio confidence interval: There a ".($refLabelingInfo->{'CONFIDENCE_LEVEL'}[0]*100)."% certainty that true ratio is included between lower and upper boundaries.");
		#	$worksheet2->write_string($xlsRow+1,++$xlsCol,'Upper bound.',$itemFormat{'header'});
		#	$worksheet2->write_comment($xlsRow+1,$xlsCol,"Upper boundary of ratio confidence interval: There a ".($refLabelingInfo->{'CONFIDENCE_LEVEL'}[0]*100)."% certainty that true ratio is included between lower and upper boundaries.");
		#}
		$worksheet2->write_string($xlsRow+1,++$xlsCol,'p-value',$itemFormat{'header'});
		$worksheet2->write_string($xlsRow+1,++$xlsCol,'CV %',$itemFormat{'header'}) if $quantifSoftware eq 'myProMS';
		if ($ratioType=~/S\w+Ratio/) {
			$worksheet2->write_string($xlsRow+1,++$xlsCol,'Dist. pept. used',$itemFormat{'header'});
			$worksheet2->write_string($xlsRow+1,++$xlsCol,'Pept. used',$itemFormat{'header'});
		}
	}
	#<Peptides
	if ($ratioType eq 'Ratio') {
		$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow,$xlsCol+2,$peptideTitle,$itemFormat{'mergeColHeader'});
		$worksheet2->write_string($xlsRow+1,$xlsCol,'Distinct used',$itemFormat{'header'});
		$worksheet2->write_string($xlsRow+1,++$xlsCol,'Used',$itemFormat{'header'});
		$worksheet2->write_string($xlsRow+1,++$xlsCol,'Total',$itemFormat{'header'});
	}
	else {
		$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow+1,$xlsCol,"Total $peptideTitle",$itemFormat{'mergeRowHeader'});
	}
	#<MW, description & species
	$worksheet2->merge_range($xlsRow,++$xlsCol,$xlsRow+1,$xlsCol,'MW (kDa)',$itemFormat{'mergeRowHeader'});
	$worksheet2->set_column(++$xlsCol,$xlsCol,80); # col length
	$worksheet2->merge_range($xlsRow,$xlsCol,$xlsRow+1,$xlsCol,'Description',$itemFormat{'mergeRowHeader'});
	$worksheet2->set_column(++$xlsCol,$xlsCol,30); # col length
	$worksheet2->merge_range($xlsRow,$xlsCol,$xlsRow+1,$xlsCol,'Species',$itemFormat{'mergeRowHeader'});

	###<Looping through proteins>###
	$xlsRow++;
	my $numDispItems=0;
	my %dispProteins;
	foreach my $modProtID (sort{&sortProteins($dispSort,$refProteinInfo,$refQuantifValues,$firstRatioUsed)} keys %{$refQuantifValues}) { #%{$refProteinInfo}
#last if $numDispItems >= 5;
		my ($protID,$modResStrg)=($modProtID=~/^(\d+)(.*)/); # eg. 12345-Y40,T152 -> 12345, -Y40,T152
		$modResStrg='' unless $modResStrg;
		$xlsCol=0;
		my $refUsedRatioPos=0;
		#<Peptide filter
		my $okPeptide=0;
		foreach my $ratioPos (keys %dispRatios) {
			next if (!$refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos} || $refQuantifValues->{$modProtID}{$numPepCode}{$ratioPos} < $dispNumPep); # may not be defined for +/-inf ratios
			$okPeptide=1; # at least 1 ratio must pass the filter
			$refUsedRatioPos=$ratioPos;
			last;
		}
		next unless $okPeptide;

		#<Ratio/p-value/Std dev filters
		my $okFilters=0;
		foreach my $ratioPos (keys %dispRatios) {
			unless ($okFilters) { # no ratioPos has passed filters
				if (defined $refQuantifValues->{$modProtID}{RATIO}{$ratioPos}) {
					if ($foldChangeType eq 'abs') {
						$okFilters=(($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} < 1 && $refQuantifValues->{$modProtID}{RATIO}{$ratioPos} <= $invDispFoldChange) || ($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} >= 1 && $refQuantifValues->{$modProtID}{RATIO}{$ratioPos} >= $dispFoldChange))? 1 : 0;
					}
					elsif ($foldChangeType eq 'up') {
						$okFilters=($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} >= 1 && $refQuantifValues->{$modProtID}{RATIO}{$ratioPos} >= $dispFoldChange)? 1 : 0;
					}
					else { # down
						$okFilters=($refQuantifValues->{$modProtID}{RATIO}{$ratioPos} <= 1 && $refQuantifValues->{$modProtID}{RATIO}{$ratioPos} <= $invDispFoldChange)? 1 : 0;
					}
				}
				$okFilters=($okFilters && ($dispPvalue>=1 || (defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} && $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos} <= $dispPvalue)))? 1 : 0;
				if ($okFilters && $dispCV && defined $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}) { # undefined always passes filter!
					if ($ratioType eq 'Ratio') { # v1
						$okFilters=0 if abs(log($refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos})/log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos})) > $dispCVFrac;
					}
					else { # no CV filter for TNPQ or MaxQuant
						$okFilters=0 if $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}/abs(log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos})/$log2) > $dispCVFrac;
					}
				}
			}
		}
		next unless $okFilters;
		$numDispItems++;
		$dispProteins{$protID}=1;
		# Identifier
		$worksheet2->write_string(++$xlsRow,$xlsCol,$refProteinInfo->{$protID}[0].$modResStrg,$itemFormat{'text'});
		# Gene(s)
		$worksheet2->write_string($xlsRow,++$xlsCol,join(',',@{$refProteinInfo->{$protID}[5]}),$itemFormat{'text'});
		#<Ratios
		foreach my $ratioPos (sort{$a<=>$b} keys %dispRatios) {
			$refUsedRatioPos=$ratioPos unless $refUsedRatioPos;
			if (!$refQuantifValues->{$modProtID}{RATIO}{$ratioPos}) { # !$refQuantifValues->{$protID}{'RATIO'} ||
				foreach my $i (0..$ratioColDelta{$ratioPos}) {
					$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
				}
				next;
			}
			# Ratio
			$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{RATIO}{$ratioPos},$itemFormat{'number'});
			# log2 ratio
			my $log2Ratio=($refQuantifValues->{$modProtID}{RATIO}{$ratioPos}==1000)? 1000 : ($refQuantifValues->{$modProtID}{RATIO}{$ratioPos}==0.001)? -1000 : log($refQuantifValues->{$modProtID}{RATIO}{$ratioPos})/$log2;
			$worksheet2->write_number($xlsRow,++$xlsCol,$log2Ratio,$itemFormat{'number'});
			#if ($ratioType eq 'Ratio') {
			#	# low conf limit
			#	if ($refQuantifValues->{$protID}{'CONF_LOW'} && defined $refQuantifValues->{$modProtID}{'CONF_LOW'}{$ratioPos}) {
			#		$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'CONF_LOW'}{$ratioPos},$itemFormat{'number'});
			#	}
			#	else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
			#	# upper conf limit
			#	if ($refQuantifValues->{$protID}{'CONF_UP'} && defined $refQuantifValues->{$modProtID}{'CONF_UP'}{$ratioPos}) {
			#		$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'CONF_UP'}{$ratioPos},$itemFormat{'number'});
			#	}
			#	else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
			#}
			# p-value
			if (defined $refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos}) {
				$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{$pvalueCode}{$ratioPos},$itemFormat{'number'});
			}
			else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}

			if ($quantifSoftware eq 'myProMS') {
				# Std dev -> CV
				if ($refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}) {
					my $cv=($ratioType eq 'Ratio')? abs(log($refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos})/log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos})) : $refQuantifValues->{$modProtID}{'SD_GEO'}{$ratioPos}/abs(log($refQuantifValues->{$modProtID}{'RATIO'}{$ratioPos})/$log2);
					$worksheet2->write_number($xlsRow,++$xlsCol,100*$cv,$itemFormat{'number'});
				}
				else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
			}
			# SuperRatio peptides
			if ($ratioType=~/S\w+Ratio/) {
				# DIST_PEP_USED
				if ($refQuantifValues->{$modProtID}{'DIST_PEP_USED'} && $refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$ratioPos} && $refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$ratioPos}) {
					$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$ratioPos},$itemFormat{'number'});
				}
				else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
				# NUM_PEP_USED
				if ($refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$ratioPos} && $refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$ratioPos}) {
					$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$ratioPos},$itemFormat{'number'});
				}
				else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
			}
		}
		# Peptides
		if ($refQuantifValues->{$modProtID}) {
			if ($ratioType eq 'Ratio') {
				# DIST_PEP_USED
				if ($refQuantifValues->{$modProtID}{'DIST_PEP_USED'} && $refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$refUsedRatioPos} && $refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$refUsedRatioPos}) {
					$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'DIST_PEP_USED'}{$refUsedRatioPos},$itemFormat{'number'});
				}
				else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
				# NUM_PEP_USED
				if ($refQuantifValues->{$modProtID}{'NUM_PEP_USED'} && $refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$refUsedRatioPos}) {
					$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'NUM_PEP_USED'}{$refUsedRatioPos},$itemFormat{'number'});
				}
				else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
			}
			# NUM_PEP_TOTAL
			if ($refQuantifValues->{$modProtID}{'NUM_PEP_TOTAL'} && $refQuantifValues->{$modProtID}{'NUM_PEP_TOTAL'}{$refUsedRatioPos}) {
				$worksheet2->write_number($xlsRow,++$xlsCol,$refQuantifValues->{$modProtID}{'NUM_PEP_TOTAL'}{$refUsedRatioPos},$itemFormat{'number'});
			}
			else {$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});}
		}
		else {
			if ($ratioType eq 'Ratio') {
				$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
				$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
			}
			$worksheet2->write_blank($xlsRow,++$xlsCol,$itemFormat{'number'});
		}
		# MW, desc, species
		my $mass=(!$refProteinInfo->{$protID}[1])? 0 : $refProteinInfo->{$protID}[1];
		$worksheet2->write_number($xlsRow,++$xlsCol,$mass,$itemFormat{'number1d'});
		$worksheet2->write_string($xlsRow,++$xlsCol,$refProteinInfo->{$protID}[2],$itemFormat{'textWrap'});
		$worksheet2->write_string($xlsRow,++$xlsCol,$refProteinInfo->{$protID}[3],$itemFormat{'text'});
	}
	#<Identifier header (written last to get remaining number of proteins)
	#my $numTotQuantifItems=scalar keys %{$refQuantifValues};
	my $protString=($quantifModID)? "$numDispItems/$numTotQuantifItems isoforms\n".(scalar keys %dispProteins)."/$numTotProteins proteins" : "$numDispItems/$numTotProteins proteins";
	$worksheet2->merge_range(0,0,1,0,$protString,$itemFormat{'mergeRowHeader'});
}

sub ajaxListSelectedProteins {
	#print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
	my @selProteins=split(',',param('selProt'));
	#my ($minFoldChange,$maxFoldChange,$maxPvalue)=split(',',param('thresholds'));
	my @targetPosList=split(',',param('dispTargetPosAjax'));
	#foreach my $tgtPos (@targetPosList) {$dispRatios{$tgtPos}=1;} # %dispRatios is global but not defined yet

	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;
	my $projectID=&promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION');
	my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
	$projectAccess=${$userInfo[2]}{$projectID};


	###>Fetching protein quantification data
	my $sthQP2=$dbh->prepare("SELECT M.CODE,ID_QUANTIF_PARAMETER,P.NAME,P.CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD M,QUANTIFICATION_PARAMETER P WHERE Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND M.ID_QUANTIFICATION_METHOD=P.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$selQuantifID");
	$sthQP2->execute;
	while (my ($qType,$paramID,$paramName,$paramCode)=$sthQP2->fetchrow_array) {
		$quantifType=$qType;
		@{$quantifParamInfo{$paramCode}}=($paramID,$paramName);
	}
	$sthQP2->finish;

	####>Protein ratio (from peptides)<####
	my %labelingInfo;
	(my $quantifAnnot,$designID,$quantifModID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_DESIGN,ID_MODIFICATION FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
	my ($labelType)=($labelStrg)? ($labelStrg=~/LABEL=(.+)/) : ('FREE');
	$labelingInfo{'LABELTYPE'}=$labelType;
	my $ratioType=($labelingInfo{'RATIO_TYPE'})? $labelingInfo{'RATIO_TYPE'}[0] : (!$labelingInfo{'RATIOS'})? 'None' : 'Ratio';
	my $algoVersion=($ratioType eq 'Ratio')? 1 : (!$labelingInfo{'SOFTWARE'})? 2 : ($labelingInfo{'SOFTWARE'}[0] eq 'myProMS')? $labelingInfo{'SOFTWARE'}[1] : 0; # myProMS version ONLY!
	my $isPepIntensity=($labelingInfo{'MEAN_STATE'} && $labelingInfo{'MEAN_STATE'}[0]==1)? 1 : ($algoVersion==2 && $labelType eq 'FREE')? 1 : ($algoVersion==3 && $labelingInfo{'ALGO_TYPE'}[0] eq 'PEP_INTENSITY')? 1 : 0; # 0 for v1 Ratio, v2 Label & v3 PEP_RATIO


	#if ($quantifType eq 'PROT_RATIO_PEP') {
		foreach my $infoStrg (@labelInfo) {
			my ($setting,$valueStrg)=split('=',$infoStrg);
			#if ($setting eq 'RATIOS' && $designID) { # $parentQuantifType && $parentQuantifType eq 'XIC'
			#	$valueStrg=~s/#//g; # for design experiments, RATIOS are displayed with '#' for condition IDs
			#}
			$valueStrg=~s/#//g if $designID; # for design experiments, RATIOS are displayed with '#' for condition IDs
			@{$labelingInfo{$setting}}=split(';',$valueStrg);
		}
		$ratioType=($labelingInfo{'RATIO_TYPE'})? $labelingInfo{'RATIO_TYPE'}[0] : 'Ratio';
		my $numStates=scalar @{$labelingInfo{'STATES'}};
		%ratioCode=&generateRatioCodes($labelType,$ratioType,$numStates,$isPepIntensity,$algoVersion) if $ratioType ne 'None';
	#}
	my $quantifSoftware=($labelingInfo{SOFTWARE})? $labelingInfo{SOFTWARE}[0] : 'myProMS';
	$quantifSoftware='MaxQuant' if $quantifSoftware eq 'MQ';
	my %stateInfo;
	my $statePos=0;
	my $sthgetExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
	foreach my $stateData (@{$labelingInfo{STATES}}) {
		$statePos++;
		(my $numBioRep,$stateInfo{$statePos}{'NAME'},my $repPosStrg)=split(',',$stateData);
		$stateInfo{$statePos}{'NAME'}=~s/\./\+/g;
		if ($designID){
			my $expCondID=$repPosStrg;
			#$expCondID=~ s/#//;
			$condToState{$expCondID}=$statePos;
			$sthgetExpCondName->execute($expCondID);
			($stateInfo{$statePos}{'NAME'})=$sthgetExpCondName->fetchrow_array;
		}
	}
	$sthgetExpCondName->finish;

	my (%quantifValues,%proteinInfo,%codedModStrings);
	my $sthProtQ2=($selProteins[0]=~/\D/)? # modif quantif
		$dbh->prepare("SELECT TARGET_POS,QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION  P
						LEFT JOIN PROTQUANTIF_MODRES PQMR ON P.ID_PROT_QUANTIF=PQMR.ID_PROT_QUANTIF
						LEFT JOIN MODIFIED_RESIDUE MR ON PQMR.ID_MODIF_RES=MR.ID_MODIF_RES
						WHERE P.ID_QUANTIFICATION=$selQuantifID AND ID_QUANTIF_PARAMETER=? AND ID_PROTEIN=? GROUP BY P.ID_PROT_QUANTIF
						HAVING GROUP_CONCAT(RESIDUE,POSITION ORDER BY POSITION SEPARATOR '.')=?")
		: $dbh->prepare("SELECT TARGET_POS,QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID AND ID_QUANTIF_PARAMETER=? AND ID_PROTEIN=?");
	#my $sthProtI2=$dbh->prepare("SELECT ALIAS,MW,PROT_DES,ORGANISM FROM ANALYSIS_PROTEIN A,PROTEIN P WHERE A.ID_PROTEIN=P.ID_PROTEIN AND VISIBILITY>=1 AND ID_ANALYSIS=$analysisID AND P.ID_PROTEIN=?");
	my $sthProtI2=$dbh->prepare("SELECT ALIAS,MW,PROT_DES,ORGANISM,PROT_LENGTH FROM ANALYSIS_PROTEIN A,PROTEIN P WHERE A.ID_PROTEIN=P.ID_PROTEIN AND ID_ANALYSIS=? AND P.ID_PROTEIN=?"); # Modif on 23/11/2012 because for TNPQ / Prot Pep Ratio, the protein could not be found in a single analysis
	my ($GNidentID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
	my $sthGN=$dbh->prepare("SELECT MI.VALUE FROM PROTEIN P,MASTERPROT_IDENTIFIER MI WHERE P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND MI.ID_IDENTIFIER=$GNidentID AND P.ID_PROTEIN=? ORDER BY MI.RANK");
	# VISIBILITY>=1 removed because of virtual proteins. Why was it there in the 1st place?... (PP 12/07/13)
	my @params;
	my $pvalueCode='';
	if ($ratioType eq 'None') {
		foreach my $statePos (@targetPosList) {$dispStates{$statePos}=1;} # %dispStates is global but not defined yet
		@params=($dispMeasure);
	}
	else {
		foreach my $ratioPos (@targetPosList) {$dispRatios{$ratioPos}=1;} # %dispRatios is global but not defined yet
		@params=('RATIO');
		if ($quantifSoftware eq 'MaxQuant') {
			#push @params,'RATIO_VAR';
		}
		else {
			$pvalueCode=($ratioType=~/S\w+Ratio/ || $labelingInfo{'FDR_CONTROL'}[0] eq 'TRUE')? 'PVAL_ADJ' : 'PVAL';
			push @params,$pvalueCode;
			push @params,'SD_GEO' if ($quantifType eq 'PROT_RATIO_PEP' && $ratioType eq 'Ratio'); # no SD_GEO for TNPQ!!!
		}
	}
	#push @params,'DIST_PEP_USED' if ($ratioType=~/S\w+Ratio/ && $numPepCode eq 'DIST_PEP_USED'); # for all ratios
	foreach my $modProtID (@selProteins) {
		my ($protID,$modStrg)=($modProtID=~/^(\d+)-*(.*)/);
		$codedModStrings{$modStrg}=&dbEncodeModifPosition($modStrg) if ($modStrg && !$codedModStrings{$modStrg});
		foreach my $paramCode (@params) {
			if ($modStrg) {$sthProtQ2->execute($quantifParamInfo{$paramCode}[0],$protID,$codedModStrings{$modStrg});}
			else {$sthProtQ2->execute($quantifParamInfo{$paramCode}[0],$protID);}
			while (my ($targetPos,$qValue)=$sthProtQ2->fetchrow_array) {
				if ($ratioType eq 'None') {
					if ($dispStates{$targetPos}) {$quantifValues{$modProtID}{$paramCode}{$targetPos}=$qValue;}
				}
				else {
					if ($dispRatios{$targetPos}) {$quantifValues{$modProtID}{$paramCode}{$targetPos}=$qValue;} #($paramCode eq 'SD_GEO')? $qValue*100 :
					elsif ($dispRatios{-$targetPos}) { # reverse ratio management
						$qValue=1/$qValue if $paramCode=~/^RATIO/;
						$quantifValues{$modProtID}{$paramCode}{-$targetPos}=$qValue;
					}
				}
			}
		}
	}
	my @pepParams;
	my (%linkedRatios,%usedPrimRatios);
	if ($quantifType eq 'MQ') {
		@pepParams=('PEPTIDES','RAZ_UNI_PEP','UNIQUE_PEP');
	}
	elsif ($ratioType=~/S\w+Ratio/) { # Old quantifs have no peptide data for super ratios: use those from linked primary ratios
		if ($quantifSoftware eq 'MaxQuant') {
			@pepParams=('PEPTIDES','RAZ_UNI_PEP','UNIQUE_PEP','NUM_PEP_USED');
		}
		else {
			my $rPos=0;
			foreach my $ratioData (@{$labelingInfo{'RATIOS'}}) { # Find linked primary ratios
				$rPos++;
				my ($testStatePos,$refStatePos)=split(/\//,$ratioData);
				if ($testStatePos=~/%(\d+)/) { # ratio of ratio
					#next unless $dispRatios{$rPos};
					next if (!$dispRatios{$rPos} && !$dispRatios{-$rPos});
					my $linkedTestPos=$1;
					my ($linkedRefPos)=($refStatePos=~/%(\d+)/);
					$testStatePos=~s/%\d+//;
					$refStatePos=~s/%\d+//;
					@{$linkedRatios{$rPos}}=($linkedTestPos,$linkedRefPos);
					%{$usedPrimRatios{$linkedTestPos}}=();
					%{$usedPrimRatios{$linkedRefPos}}=();
				}
			}
			my @pepSupParams=('NUM_PEP_USED'); # NUM_PEP_USED always needed
			push @pepSupParams,'DIST_PEP_USED' if $numPepCode eq 'DIST_PEP_USED'; # for all ratios
			foreach my $modProtID (@selProteins) {
				my ($protID,$modStrg)=($modProtID=~/^(\d+)-*(.*)/);
				foreach my $paramCode (@pepSupParams) {
					if ($modStrg) {$sthProtQ2->execute($quantifParamInfo{$paramCode}[0],$protID,$codedModStrings{$modStrg});}
					else {$sthProtQ2->execute($quantifParamInfo{$paramCode}[0],$protID);}
					while (my ($ratioPos,$qValue)=$sthProtQ2->fetchrow_array) {
						$usedPrimRatios{$ratioPos}{$modProtID}{$paramCode}=$qValue if $usedPrimRatios{$ratioPos};
						#$quantifValues{$modProtID}{$paramCode}{$ratioPos}=$qValue if $dispRatios{$ratioPos};
						if ($dispRatios{$ratioPos}) {$quantifValues{$modProtID}{$paramCode}{$ratioPos}=$qValue;}
						elsif ($dispRatios{-$ratioPos}) {$quantifValues{$modProtID}{$paramCode}{-$ratioPos}=$qValue;} # rev ratio
					}
				}
			}
			#>Replace missing pep data with those from primary ratios
			RPOS:foreach my $ratioPos (keys %linkedRatios) {
				my $usedRatio=($dispRatios{$ratioPos})? $ratioPos : ($dispRatios{-$ratioPos})? -$ratioPos : 0;
				next unless $usedRatio;
				foreach my $modProtID (keys %quantifValues) {
					#if ($quantifValues{$modProtID}{'RATIO'}{$ratioPos}) { # to be displayed
					#	foreach my $paramCode (@pepSupParams) {
					#		last RPOS if $quantifValues{$modProtID}{$paramCode}{$ratioPos}; # numPep stored in DB for 2ndary ratio
					#		# Compute from primary ratios
					#		$quantifValues{$modProtID}{$paramCode}{$ratioPos}=($usedPrimRatios{$linkedRatios{$ratioPos}[0]}{$modProtID} && $usedPrimRatios{$linkedRatios{$ratioPos}[0]}{$modProtID}{$paramCode} <= $usedPrimRatios{$linkedRatios{$ratioPos}[1]}{$modProtID}{$paramCode})? $usedPrimRatios{$linkedRatios{$ratioPos}[0]}{$modProtID}{$paramCode} : $usedPrimRatios{$linkedRatios{$ratioPos}[1]}{$modProtID}{$paramCode};
					#	}
					#}
					foreach my $paramCode (@pepSupParams) {
						last RPOS if $quantifValues{$modProtID}{$paramCode}{$usedRatio}; # numPep stored in DB for 2ndary ratio
						# Compute from primary ratios
						$quantifValues{$modProtID}{$paramCode}{$usedRatio}=($usedPrimRatios{$linkedRatios{$ratioPos}[0]}{$modProtID} && $usedPrimRatios{$linkedRatios{$ratioPos}[0]}{$modProtID}{$paramCode} <= $usedPrimRatios{$linkedRatios{$ratioPos}[1]}{$modProtID}{$paramCode})? $usedPrimRatios{$linkedRatios{$ratioPos}[0]}{$modProtID}{$paramCode} : $usedPrimRatios{$linkedRatios{$ratioPos}[1]}{$modProtID}{$paramCode};
					}
				}
			}
			@pepParams=('NUM_PEP_TOTAL');
		}
	}
	else {
		@pepParams=($quantifType eq 'TNPQ')? ('NUM_PEP_USED','NUM_PEP_TOTAL') : ($numPepCode eq 'DIST_PEP_USED')? ('DIST_PEP_USED','NUM_PEP_USED','NUM_PEP_TOTAL') : ('NUM_PEP_USED','NUM_PEP_TOTAL');
	}

	foreach my $modProtID (@selProteins) {
		my ($protID,$modStrg)=($modProtID=~/^(\d+)-*(.*)/);
		foreach my $paramCode (@pepParams) {
			if ($modStrg) {$sthProtQ2->execute($quantifParamInfo{$paramCode}[0],$protID,$codedModStrings{$modStrg});}
			else {$sthProtQ2->execute($quantifParamInfo{$paramCode}[0],$protID);}
			while (my ($targetPos,$qValue)=$sthProtQ2->fetchrow_array) {
				if ($ratioType eq 'None') { # MaxQuant
					foreach my $tPos (keys %dispStates) { # global: no targetPos
						$quantifValues{$modProtID}{$paramCode}{$tPos}=$qValue;
					}
				}
				else {
					if ($targetPos) { # defined for +/-inf ratios
						$quantifValues{$modProtID}{$paramCode}{$targetPos}=$qValue if $dispRatios{$targetPos};
					}
					else { # undef for quantified proteins
						foreach my $rPos (keys %dispRatios) { # extend to all displayed ratios
							$quantifValues{$modProtID}{$paramCode}{$rPos}=$qValue;
						}
					}
				}
			}
		}

		###>Fetching protein info
		#if ($designID){ # $labelType eq 'FREE' && $quantifType =~ /PROT_RATIO|TNPQ/
		#	my ($anaID)=$dbh->selectrow_array("SELECT ID_ANALYSIS FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=$protID AND ID_ANALYSIS IN (SELECT ID_ANALYSIS FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID) ORDER BY NUM_PEP DESC LIMIT 1");# To be sure that the PROTEIN is in the ANALYSIS
		#	$sthProtI2->execute($anaID,$protID);
		#	@{$proteinInfo{$protID}}=$sthProtI2->fetchrow_array;
		#	$sthGN->execute($protID);
		#	my @geneList;
		#	while (my ($gene)=$sthGN->fetchrow_array) {push @geneList,$gene;}
		#	push @{$proteinInfo{$protID}},(\@geneList,$anaID);
		#}
		#else{
		#	$sthProtI2->execute($analysisID,$protID);
		#	@{$proteinInfo{$protID}}=$sthProtI2->fetchrow_array;
		#	$sthGN->execute($protID);
		#	my @geneList;
		#	while (my ($gene)=$sthGN->fetchrow_array) {push @geneList,$gene;}
		#	push @{$proteinInfo{$protID}},(\@geneList,$analysisID);
		#}
		my $usedAnaID;
		if ($designID){ # $labelType eq 'FREE' && $quantifType =~ /PROT_RATIO|TNPQ/
			($usedAnaID)=$dbh->selectrow_array("SELECT ID_ANALYSIS FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=$protID AND ID_ANALYSIS IN (SELECT ID_ANALYSIS FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID) ORDER BY NUM_PEP DESC LIMIT 1");# To be sure that the PROTEIN is in the ANALYSIS
		}
		else {$usedAnaID=$analysisID;}
		$sthProtI2->execute($usedAnaID,$protID);
		@{$proteinInfo{$protID}}=$sthProtI2->fetchrow_array;
		$sthGN->execute($protID);
		my @geneList;
		while (my ($gene)=$sthGN->fetchrow_array) {push @geneList,$gene;}
		push @{$proteinInfo{$protID}},(\@geneList,$usedAnaID);

		$proteinInfo{$protID}[1]=sprintf "%.1f",$proteinInfo{$protID}[1]/1000; # MW
	}
	$sthProtQ2->finish;
	$sthProtI2->finish;

	$dbh->disconnect;

	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);
	if ($quantifSoftware eq 'MaxQuant') {
		&printMaxQuantProteinList($labelType,\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo);
	}
	else {
		&printProteinRatioList($labelType,$ratioType,\%labelingInfo,\%stateInfo,\%quantifValues,\%proteinInfo,$pvalueCode); #,10
	}
	exit;
}

####################<<<ajaxPeptideRatios>>>#####################
sub ajaxPeptideRatios {
	# TODO: +/- infinite ratio: Add missing filters for peptides (charge, source)
	my $modProtID=param('id_prot');
	my ($protID,$modStrg)=($modProtID=~/^(\d+)-*(.*)/);
	my $selRatioPos=param('ratio');
	my $trueRatioPos=abs($selRatioPos);
	#my $selRatioIdx=$selRatioPos-1;
#print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;
	my $projectID=&promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION');
	my ($quantifAnnot,$designID,$quantifType,$quantifModID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_DESIGN,M.CODE,Q.ID_MODIFICATION FROM QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$selQuantifID");

	my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
	my ($labelType)=($labelStrg)? $labelStrg=~/LABEL=(.+)/ : ('FREE');
	$labelType=uc($labelType); # iTRAQ to ITRAQ
	my (%labelingInfo,%condToState,%condToObs);
	$labelingInfo{'LABELTYPE'}=$labelType;
	#if ($designID || $labelType =~ /SILAC|ITRAQ/) {
		foreach my $infoStrg (@labelInfo) {
			my ($setting,$valueStrg)=split('=',$infoStrg);
			$valueStrg=~s/#//g; # remove all ID tags
			@{$labelingInfo{$setting}}=split(';',$valueStrg);
		}
	#}
	my $ratioType=($labelingInfo{RATIO_TYPE})? $labelingInfo{RATIO_TYPE}[0] : 'Ratio';
	my $algoVersion=($ratioType eq 'Ratio')? 1 : ($labelingInfo{SOFTWARE} && $labelingInfo{SOFTWARE}[0] eq 'myProMS')? 3 : 2;
my $isPepIntensity=($algoVersion==2 && $labelType eq 'FREE')? 1 : ($algoVersion==3 && $labelingInfo{ALGO_TYPE}[0] eq 'PEP_INTENSITY')? 1 : 0; # 0 for v1 Ratio, v2 Label & v3 PEP_RATIO
	my $numRatios=scalar @{$labelingInfo{RATIOS}};
	my $numStates=scalar @{$labelingInfo{STATES}};
	my $idxShift=($quantifType eq 'PROT_RATIO_PEP')? 0 : 2;

	if ($modStrg && !$quantifModID) { # fall back to whole protID if called with a modProtID for a non-modif quantif (eg. called from log2 plot)
		$modProtID=$protID;
		$modStrg='';
	}
	my ($protAlias,$protLength)=$dbh->selectrow_array("SELECT ALIAS,PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=$protID");
	$protLength=0 unless $protLength;

	my $sthPR;
	if ($modStrg) { # modification quantif
		my $dbModStrg=&dbEncodeModifPosition($modStrg);
		$sthPR=$dbh->prepare("SELECT QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION P
											LEFT JOIN PROTQUANTIF_MODRES PQMR ON P.ID_PROT_QUANTIF=PQMR.ID_PROT_QUANTIF
											LEFT JOIN MODIFIED_RESIDUE MR ON PQMR.ID_MODIF_RES=MR.ID_MODIF_RES
											INNER JOIN QUANTIFICATION_PARAMETER Q ON P.ID_QUANTIF_PARAMETER=Q.ID_QUANTIF_PARAMETER AND Q.CODE=?
											WHERE P.ID_QUANTIFICATION=$selQuantifID AND ID_PROTEIN=$protID AND TARGET_POS=? GROUP BY P.ID_PROT_QUANTIF
											HAVING GROUP_CONCAT(RESIDUE,POSITION ORDER BY POSITION SEPARATOR '.')='$dbModStrg'");
	}
	else { # whole protein quantif
		$sthPR=$dbh->prepare("SELECT QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION P,QUANTIFICATION_PARAMETER Q WHERE P.ID_QUANTIF_PARAMETER=Q.ID_QUANTIF_PARAMETER AND Q.CODE=? AND ID_QUANTIFICATION=$selQuantifID AND ID_PROTEIN=$protID AND TARGET_POS=?");
	}


	##<Fetching labeling type
	#my ($quantifAnnot,$designID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_DESIGN FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	my ($pepQuantifIDs)=$dbh->selectrow_array("SELECT GROUP_CONCAT(ID_PARENT_QUANTIFICATION SEPARATOR ',') FROM PARENT_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID AND PAR_FUNCTION IS NULL GROUP BY ID_QUANTIFICATION");
	my ($apexRTParamID)=$dbh->selectrow_array("SELECT ID_QUANTIF_PARAMETER FROM QUANTIFICATION_PARAMETER QP,QUANTIFICATION Q WHERE QP.ID_QUANTIFICATION_METHOD=Q.ID_QUANTIFICATION_METHOD AND CODE='RT_APEX' AND ID_QUANTIFICATION IN ($pepQuantifIDs) LIMIT 0,1");

	my ($trueProtRatio,$trueTestStateName,$trueRefStateName,$is2ndaryRatio);
	my (%protRatio,%testStatePos,%refStatePos,%testStateName,%refStateName,%anaList,%protMean); #,%usedReplicates
my %replicDesign; # v3
	if ($designID) { # Design quanti ***statePos=CondID but converted later in true pos by %condToState***
		my ($testStateID,$refStateID)=($selRatioPos > 0)? split(/\//,$labelingInfo{RATIOS}[$trueRatioPos-1]) : reverse (split(/\//,$labelingInfo{RATIOS}[$trueRatioPos-1]));
		my ($testRatioPos,$refRatioPos);
		if ($testStateID=~/%/) {  # && $algoVersion < 3) { #} 2NDARY RATIO!!!
			$is2ndaryRatio=1;
			#$testStatePos=~s/%\d+//;
			#$refStatePos=~s/%\d+//;
			($testRatioPos)=($testStateID=~/%(\d+)/);
			($refRatioPos)=($refStateID=~/%(\d+)/);
			foreach my $ratioPos ($testRatioPos,$refRatioPos) {
				($testStatePos{$ratioPos},$refStatePos{$ratioPos})=split(/\//,$labelingInfo{RATIOS}[$ratioPos-1]); # Converted to real pos below
			}
			$sthPR->execute('RATIO',$trueRatioPos);
			($trueProtRatio)=$sthPR->fetchrow_array;
			$trueProtRatio=1/$trueProtRatio if $selRatioPos < 0; # reverse
		}
		else { # PRIMARY RATIO or FREE
			$testStatePos{$selRatioPos}=$testStateID; # Converted to real pos below
			$refStatePos{$selRatioPos}=$refStateID; # Converted to real pos below
		}
		my $statePos=0;
		foreach my $stateInfo (@{$labelingInfo{STATES}}) {
			$statePos++;
			my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$stateInfo);
my $curBioRep=0;
foreach my $bioReplicate (split(/\./,$quantiObsIDs)) {
	$curBioRep++;
	$replicDesign{$statePos}{$curBioRep}=($quantiObsIDs=~/&/)? 1 : 0;
}
			#$expCondID=~ s/^.//;
			$condToState{$expCondID}=$statePos; # absolute pos in all states list
			$condToObs{$expCondID}=$quantiObsIDs;
		}
		my $sthExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");

		foreach my $ratioPos (keys %testStatePos) {
			$sthPR->execute('RATIO',abs($ratioPos));
			my ($pRatio)=$sthPR->fetchrow_array;
			if ($is2ndaryRatio || $selRatioPos > 0) {$protRatio{$ratioPos}=$pRatio;}
			else {$protRatio{$ratioPos}=1/$pRatio;}
			my $testStateExpCondID=$testStatePos{$ratioPos};
			$sthExpCondName->execute($testStateExpCondID);
			($testStateName{$ratioPos})=$sthExpCondName->fetchrow_array;
			my $refStateExpCondID=$refStatePos{$ratioPos};
			$sthExpCondName->execute($refStateExpCondID);
			($refStateName{$ratioPos})=$sthExpCondName->fetchrow_array;
			$testStatePos{$ratioPos}=$condToState{$testStateExpCondID}; # absolute pos in all states list
			$refStatePos{$ratioPos}=$condToState{$refStateExpCondID}; # absolute pos in all states list
			#if ($protRatio{$ratioPos} >= $MIN_INF_RATIO_DB) {$usedReplicates{$ratioPos}=$condToObs{$refStateExpCondID};}
			#elsif ($protRatio{$ratioPos} <= $MAX_INF_RATIO_DB) {$usedReplicates{$ratioPos}=$condToObs{$testStateExpCondID};}
			if ($isPepIntensity) { # $labelType eq 'FREE' && $ratioType ne 'Ratio'
				$sthPR->execute('MEAN_STATE',$numRatios+$testStatePos{$ratioPos});
				($protMean{$testStatePos{$ratioPos}})=$sthPR->fetchrow_array;
				$sthPR->execute('MEAN_STATE',$numRatios+$refStatePos{$ratioPos});
				($protMean{$refStatePos{$ratioPos}})=$sthPR->fetchrow_array;
#print "TEST=$protMean{$testStatePos{$ratioPos}}, REF=$protMean{$refStatePos{$ratioPos}}<BR>\n";
			}
		}
		$sthExpCondName->finish;
		if ($is2ndaryRatio) { # SuperRatio
			$trueTestStateName=$testStateName{$testRatioPos}.encode_utf8('°');
			$trueRefStateName=$testStateName{$refRatioPos}.encode_utf8('°');
		}
		else {
			$trueTestStateName=$testStateName{$selRatioPos};
			$trueRefStateName=$refStateName{$selRatioPos};
		}
	}
	else { # internal quanti
		$sthPR->execute('RATIO',$trueRatioPos);
		($protRatio{$selRatioPos})=$sthPR->fetchrow_array;
		($testStatePos{$selRatioPos},$refStatePos{$selRatioPos})=($selRatioPos > 0)? split(/\//,$labelingInfo{RATIOS}[$trueRatioPos-1]) : reverse(split(/\//,$labelingInfo{RATIOS}[$trueRatioPos-1]));
		($testStateName{$selRatioPos})=$labelingInfo{STATES}[$testStatePos{$selRatioPos}-1]=~/,(.+),/; $testStateName{$selRatioPos}=~s/\./\+/g;
		($refStateName{$selRatioPos})=$labelingInfo{STATES}[$refStatePos{$selRatioPos}-1]=~/,(.+),/; $refStateName{$selRatioPos}=~s/\./\+/g;
#print "\$selRatioPos=$selRatioPos; \$testStatePos{\$selRatioPos}=$testStatePos{$selRatioPos}, \$testStateName{\$selRatioPos}=$testStateName{$selRatioPos}<BR>\n";
# Define $usedReplicates for infinite ratio !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		#if ($protRatio <= 0.001 || $protRatio >= 1000) { # infinite ratio
		#	my ($stateChannelStrg)=($protRatio <= 0.001)? $labelingInfo{STATES}[$refStatePos-1]=~/,([^,]+)\Z/ : $labelingInfo{STATES}[$testStatePos-1]=~/,([^,]+)\Z/;
		#	# get state data from parent peptide quanti
		#	# sthQuery....
		#	foreach my $channel (split(/\./,$stateChannelStrg)) {
		#		#...
		#		# execute($channel);
		#
		#	}
		#	# ...
		#}
		#($usedReplicates{$selRatioPos})=($protRatio{$selRatioPos} <= 0.001)? $labelingInfo{STATES}[$refStatePos{$selRatioPos}-1]=~/,([^,]+)\Z/ : $labelingInfo{STATES}[$testStatePos{$selRatioPos}-1]=~/,([^,]+)\Z/;
	}
	$sthPR->finish;

	$modStrg=($modStrg)? '-'.$modStrg : '';

	#my $sthPB=$dbh->prepare("SELECT ABS(PEP_BEG) FROM PEPTIDE_PROTEIN_ATTRIB WHERE ID_PROTEIN=$protID AND ID_PEPTIDE=? ORDER BY ABS(PEP_BEG) ASC");
	#my $pepQuery=($apexRTParamID)? "SELECT SCORE,DATA,ID_ANALYSIS,QUANTIF_VALUE,QUERY_NUM,PEP_RANK FROM PEPTIDE P
	#									LEFT JOIN PEPTIDE_QUANTIFICATION PQ ON P.ID_PEPTIDE=PQ.ID_PEPTIDE AND ID_QUANTIF_PARAMETER=$apexRTParamID AND ID_QUANTIFICATION IN ($pepQuantifIDs)
	#									WHERE P.ID_PEPTIDE=?"
	#								: "SELECT SCORE,DATA,ID_ANALYSIS,NULL,QUERY_NUM,PEP_RANK FROM PEPTIDE P,PEPTIDE_QUANTIFICATION PQ WHERE P.ID_PEPTIDE=PQ.ID_PEPTIDE AND P.ID_PEPTIDE=? AND ID_QUANTIFICATION IN ($pepQuantifIDs)";
	#my $sthPD=$dbh->prepare($pepQuery);

	my $dataDir="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/data";
	my $resultDir="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/results";


my (%protPeptideData,%peptideApexRT);
my $sthProtPep=$dbh->prepare(qq
|SELECT P.ID_PEPTIDE,GROUP_CONCAT(DISTINCT ABS(PEP_BEG) ORDER BY ABS(PEP_BEG) ASC SEPARATOR ','),SCORE,DATA,P.ID_ANALYSIS,QUERY_NUM,PEP_RANK
	FROM PEPTIDE P,PEPTIDE_PROTEIN_ATTRIB PPA,ANA_QUANTIFICATION AQ
	WHERE P.ID_PEPTIDE=PPA.ID_PEPTIDE AND PPA.ID_PROTEIN=? AND P.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION IN ($pepQuantifIDs) GROUP BY P.ID_PEPTIDE
|);
$sthProtPep->execute($protID);
while (my ($pepID,$begStrg,$score,$pepData,$anaID,$qNum,$rank)=$sthProtPep->fetchrow_array) {
	$score=0 unless $score;
	$pepData='' unless $pepData;
	@{$protPeptideData{$pepID}{BEG}}=split(',',$begStrg);
	@{$protPeptideData{$pepID}{INFO}}=($score,$pepData,$anaID,$qNum,$rank);
}
$sthProtPep->finish;
if ($apexRTParamID) {
	QUANTIF:foreach my $pepQuantifID (split(',',$pepQuantifIDs)) {
		foreach my $quantFile (glob ("$promsPath{quantification}/project_$projectID/quanti_$pepQuantifID/peptide_quantification*")) {
			open(QUANTI,$quantFile);
			while (<QUANTI>) {
				next if $.==1;
				my ($paramID,$pepID,$qValue)=split(/\t/,$_);
				$peptideApexRT{$pepID}=$qValue if ($paramID==$apexRTParamID && $protPeptideData{$pepID});
			}
			close QUANTI;
			last QUANTIF if scalar keys %peptideApexRT == scalar keys %protPeptideData;
		}
	}
}


	my (%peptideSets,%sequenceBeg,%dataSources);
my %dataSrcAlias; # v3
	my (%rawPeptideSets,%pepOccurence,%varModName); #,%sequenceBeg,%dataSources

	######################
	####<Normal ratio>####
	######################
	my %normalRatio;
my %usedAnalyses; # anaID used for selected ratio (needed for filtering table.txt for SuperRatio. No longer use Experiment) v3
my %usedStates; # (Super/Simple)Ratio (v2 v3)
	my $numNormalRatios=0;
	###my ($pepFile,$pepSetFactor)=($ratioType eq 'Ratio')? ("$resultDir/pep_mean_cv_ratio.txt",2) : ($labelType eq 'FREE')? ("$resultDir/Log2MeanPep.txt",1) : ("$resultDir/RatioPep.txt",2);
	#my $pepSetFactor=($ratioType eq 'Ratio' || $labelType ne 'FREE')? 2 : 1;
	my $pepFile=($algoVersion==3)? 'resultsPep.txt' : ($ratioType eq 'Ratio')? 'pep_mean_cv_ratio.txt' : ($labelType eq 'FREE')? 'Log2MeanPep.txt' : 'RatioPep.txt';
	$pepFile="$resultDir/$pepFile";
	foreach my $ratioPos (keys %protRatio) {
		$normalRatio{$ratioPos}=0;
		if ($protRatio{$ratioPos} && $protRatio{$ratioPos} != 0.001 && $protRatio{$ratioPos} != 1000) { # test on defined protRatio because of rare possible exclusion of all peptides by global CVout
			$normalRatio{$ratioPos}=1;
			$numNormalRatios++;
			if ($ratioType eq 'Ratio') { # v1
				my %varModName;
				my $prevProtID=0;
				open(PEP,$pepFile) || die "ERROR: $!";
				while(<PEP>) {
					next if $.==1;
					chomp;
					my ($pID,$pepSeq,$varMod,$charge,$pepIdKey,$protName,$protValid,@values)=split(/\t/,$_); # ratio	mean.cond1	mean.cond2	PERM.Repl_Light	PERM.Repl_Heavy	PERM_Out	Pep_Out_0.05
					my ($trueProtID)=($pID=~/^(\d+)/);
					last if ($prevProtID==$protID && $trueProtID != $protID); # sorted by true protID but not by modProtID
					next if $pID ne $modProtID;
					$prevProtID=$trueProtID;
					my $seqVarMod=$pepSeq;
					if ($varMod && $varMod ne 'NA') {
						if ($varMod=~/^&/) { # var mod code
							$varMod=&promsMod::decodeVarMod($dbh,$pepSeq,$varMod,\%varModName);
						}
						else { # var mod text
							$varMod=~s/^ \+ //;
						}
						$seqVarMod.=' + '.$varMod;
					}
					else {$varMod='';}
					unless ($sequenceBeg{$pepSeq}) {
						my ($pep1ID)=$pepIdKey=~/^(\d+)/;
						#$sthPB->execute($pep1ID);
						#while (my ($beg)=$sthPB->fetchrow_array) {
						#	push @{$sequenceBeg{$pepSeq}},$beg;
						#}
@{$sequenceBeg{$pepSeq}}=@{$protPeptideData{$pep1ID}{BEG}};
					}
					my ($dataSrc,$scoreStrg,$scoreMin,$scoreMax,$prsStrg,$prsMin,$prsMax,$rtStrg,$rtMin,$rtMax);
					$pepIdKey=~s/\s//g; # very old quantif with spaces around ';'
					foreach my $pepID (split(/[+;]/,$pepIdKey)) {
						#$sthPD->execute($pepID);
						#my ($score,$pepData,$anaID)=$sthPD->fetchrow_array;
						#$score=0 unless $score;
						#$pepData='' unless $pepData;
my ($score,$pepData,$anaID)=@{$protPeptideData{$pepID}{INFO}}[0..2];
						#$scoreStrg.='/' if $scoreStrg;
						#$scoreStrg.=($score)? $score : '-';
						unless ($dataSrc) {
							$dataSrc='#'.$anaID;
							$dataSrc.='.'.$1 if $pepData=~/SOURCE_RANK=(\d+)/;
						}
						if ($score) {
							$scoreMin=$score unless ($scoreMin && $score>$scoreMin);
							$scoreMax=$score unless ($scoreMax && $score<$scoreMax);
						}
						if ($quantifModID && $pepData=~/PRS=\d;([^;]+);/) {
							my $prs=$1;
							$prsMin=$prs unless ($prsMin && $prsMin>$prs);
							$prsMax=$prs unless ($prsMax && $prsMax<$prs);
						}
						#$dataSrc=$dSrc || '-';
						$anaList{$anaID}=1;
					}
					$scoreStrg='-';
					if ($scoreMin && $scoreMax) {
						if ($scoreMin == $scoreMax) {$scoreStrg=1*(sprintf '%.2f',$scoreMin);}
						else {$scoreStrg='['.(1*(sprintf '%.2f',$scoreMin)).'-'.(1*(sprintf '%.2f',$scoreMax)).']';}
					}
					if ($quantifModID) {
						$prsStrg='-';
						if ($prsMin && $prsMax) {
							if ($prsMin == $prsMax) {$prsStrg=1*(sprintf '%.1f',$prsMin);}
							else {$prsStrg='['.(1*(sprintf '%.1f',$prsMin)).'-'.(1*(sprintf '%.1f',$prsMax)).']';}
						}
					}
					$dataSrc='' unless $dataSrc;
					$dataSources{$dataSrc}=$dataSrc;
					@{$peptideSets{$ratioPos}{$pepIdKey}{INFO}}=($sequenceBeg{$pepSeq}[0],$pepSeq,$varMod,$charge,$scoreStrg,{NAME=>'',DATA=>[$dataSrc]});
					my $pepRatio=($selRatioPos > 0)? $values[$trueRatioPos-1] : 1/$values[$trueRatioPos-1];
					push @{$peptideSets{$ratioPos}{$pepIdKey}{DATA}},$pepRatio,
						$values[$numRatios+$testStatePos{$ratioPos}-1],$values[$numRatios+$refStatePos{$ratioPos}-1], # mean_up, mean_down
						$values[$numRatios+$numStates+$testStatePos{$ratioPos}-1],$values[$numRatios+$numStates+$refStatePos{$ratioPos}-1],$values[$numRatios+($numStates*2)], # cv_out_up, cv_out_down, cv_out
						$values[$numRatios+($numStates*2)+$trueRatioPos+$idxShift]; # pep_out_ratio
				}
				close PEP;
			}
			else { # (Super/Simple)Ratio
				my %selStateCode; # v3
				my $rPos=0;
				PR:foreach my $x (1..$numStates-1) {
					foreach my $y ($x+1..$numStates) {
						$rPos++;
						if ($rPos==$trueRatioPos) {
							if ($algoVersion==3) {
								if ($isPepIntensity) {
									$selStateCode{"State$x"}=$x; #$ratioPos;
									$selStateCode{"State$y"}=$y; #$ratioPos;
								}
								else {$selStateCode{"State$y.State$x"}=$ratioPos;}
							}
							$usedStates{"State$x"}=$usedStates{"State$y"}=1;
							last PR;
						}
					}
					last if $is2ndaryRatio; # no need to check $labelInfo{SINGLE_REF}
				}
				if ($is2ndaryRatio) {
					SR:foreach my $x (2..$numStates-1) {
						foreach my $y ($x+1..$numStates) {
							$rPos++;
							if ($rPos==$trueRatioPos) {
								if ($algoVersion==3) {
									$selStateCode{"State$x.State1"}=$x-1; # primary ratioPos
									$selStateCode{"State$y.State1"}=$y-1; # primary ratioPos
								}
								$usedStates{"State$x"}=$usedStates{"State$y"}=$usedStates{'State1'}=1;
								last SR;
							}
						}
					}
				}

				my (%peptideData, %seenExcluded); # simpleRatio with isPepIntensity
				open(PEP,$pepFile) || die "ERROR: $!";
				if ($algoVersion==3) {
					my (%pepHeader2Idx,%replicates);
					while(<PEP>) {
						#next if $.==1;
						chomp;
						#my ($condCode,$experiment,$bioRep,$techRep,$pID,$pepCode,$log2Value,$pepIdKey,$outlier)=split(/\t/,$_);
						#ProteinID	Peptide	Experiment	replicate	repTech	Condition	proteinName	PeptideId	proteinValidity	log2Measure	normProtein	out
						#my ($pID,$pepCode,$experiment,$bioRep,$techRep,$condCode,$protName,$pepIdKey,$protValid,$log2Value,$protNorm,$outlier)=split(/\t/,$_);
						#ProteinID	Peptide	Experiment	Condition	proteinValidity	log2Measure	replicate	repTech	proteinName	PeptideId	normProtein	out # v3b
						my @lineData=split(/\t/,$_);
						####>Header------------------
						if ($.==1) { # 1st line of the file
							my $colIdx=0;
							foreach my $colName (@lineData) {
								$pepHeader2Idx{$colName}=$colIdx;
								$colIdx++;
							}
							next;
						}
						####>Data------------------
						next if $lineData[$pepHeader2Idx{'ProteinID'}] ne $modProtID;
						next unless $selStateCode{$lineData[$pepHeader2Idx{'Condition'}]};
						my $rPos=$selStateCode{$lineData[$pepHeader2Idx{'Condition'}]}; # always primary ratio or meanProt!
						my $isOutlier=($lineData[$pepHeader2Idx{'out'}] eq 'NA')? 0 : 1;
						@{$peptideData{$rPos}{$lineData[$pepHeader2Idx{'Peptide'}]}{$lineData[$pepHeader2Idx{'replicate'}].$lineData[$pepHeader2Idx{'repTech'}]}}=($lineData[$pepHeader2Idx{'log2Measure'}],$lineData[$pepHeader2Idx{'PeptideId'}],$isOutlier);
my ($statePos)=$lineData[$pepHeader2Idx{'Condition'}]=~/^State(\d+)/;
my ($bioRepPos)=$lineData[$pepHeader2Idx{'replicate'}]=~/(\d+)/;
my ($techRepPos)=$lineData[$pepHeader2Idx{'repTech'}]=~/(\d+)/;
my $numBioRep=scalar keys %{$replicDesign{$statePos}};
$dataSrcAlias{$lineData[$pepHeader2Idx{'replicate'}].$lineData[$pepHeader2Idx{'repTech'}]}=($numBioRep==1 && $replicDesign{$statePos}{$bioRepPos}==0)? 'No replic.': ($numBioRep==1)? "tech.rep$techRepPos" : ($replicDesign{$statePos}{$bioRepPos}==0)? "bio.rep$bioRepPos" : "bio.rep$bioRepPos/tech.rep$techRepPos";
					}
				}
				else { # v2 No state recorded: only experiment!
					my $selExp=chr(64+abs($ratioPos)); # cannot be 2ndary ratio!!! OK for iTRAQ in results files
					while(<PEP>) {
						next if $.==1;
						chomp;
						my ($expState,$pID,$pepCode,$log2Value,$pepIdKey,$pepExcluded,$nbPep,$log2PepStrg)=split(/\t/,$_); # $log2PepStrg defined only for label-free since FunctionLimma.R v2.4.1
						##my ($trueProtID)=($pID=~/^(\d+)/);
						##last if ($prevProtID==$protID && $trueProtID != $protID); # sorted by true protID but maybe not by modProtID <- NEW: Not even true for LABEL-FREE!!!
						next if $pID ne $modProtID;
						##$prevProtID=$trueProtID;
						next if $log2Value eq 'NA'; # No $log2Value for excluded peptide!!!!!! TODO: record a value here
						my $statePos;
						my $usedPos;
						if ($labelType eq 'FREE') {
							($statePos)=($expState=~/(\d+)/);
							next unless $protMean{$statePos};
							$usedPos=$statePos;
						}
						else {
							if ($expState ne $selExp) {next;}
							$usedPos=$ratioPos;
						}
						@{$peptideData{$usedPos}{$pepCode}{ALL_REP}}=($log2Value,$pepIdKey,$pepExcluded);
					}
				}

				close PEP;

				foreach my $usedPos (keys %peptideData) {
					foreach my $pepCode (keys %{$peptideData{$usedPos}}) {

						my ($seqModCode,$charge)=split(/_/,$pepCode);
						my ($pepSeq,$varMod);
						if ($seqModCode=~/(\w+)&(.+)/) {
							$pepSeq=$1;
							$varMod=&promsMod::decodeVarMod($dbh,$pepSeq,$2,\%varModName);
						}
						else {
							$pepSeq=$seqModCode;
							$varMod='';
						}

						foreach my $repCode (keys %{$peptideData{$usedPos}{$pepCode}}) {
							my ($log2Value,$pepIdKey,$pepExcluded)=@{$peptideData{$usedPos}{$pepCode}{$repCode}};
							unless ($sequenceBeg{$pepSeq}) {
								my ($pep1ID)=$pepIdKey=~/^(\d+)/;
								#$sthPB->execute($pep1ID);
								#while (my ($beg)=$sthPB->fetchrow_array) {
								#	push @{$sequenceBeg{$pepSeq}},$beg;
								#}
@{$sequenceBeg{$pepSeq}}=@{$protPeptideData{$pep1ID}{BEG}};
							}
						#my @pepIdList=($pepIdKey);
						#my @pepLog2Mean=($log2Value);
						#if ($log2PepStrg) { # Add unagreggated peptides
						#	$log2PepStrg=~s/^_//;
						#	if ($pepIdKey=~/_/) { # or $nbPep > 1?
						#		push @pepIdList,split(/_/,$pepIdKey);
						#		push @pepLog2Mean,split(/_/,$log2PepStrg);
						#	}
						#}
						#my %dataSrcPepIndex; # only if $log2PepStrg is defined
							my (%dataSrc,$scoreStrg,$scoreMin,$scoreMax,$prsStrg,$prsMin,$prsMax,$rtStrg,$rtMin,$rtMax);

							foreach my $pepID (split(/[+_;\.=]/,$pepIdKey)) {
								#$sthPD->execute($pepID);
								#my ($score,$pepData,$anaID)=$sthPD->fetchrow_array;
								#$score=0 unless $score;
								#$pepData='' unless $pepData;
my ($score,$pepData,$anaID)=@{$protPeptideData{$pepID}{INFO}}[0..2];
								#$scoreStrg.='/' if $scoreStrg;
								#$scoreStrg.=($score)? $score : '-';
								if ($score) {
									$scoreMin=$score unless ($scoreMin && $score>$scoreMin);
									$scoreMax=$score unless ($scoreMax && $score<$scoreMax);
								}
								if ($quantifModID && $pepData=~/PRS=\d;([^;]+);/) {
									my $prs=$1;
									$prsMin=$prs unless ($prsMin && $prsMin>$prs);
									$prsMax=$prs unless ($prsMax && $prsMax<$prs);
								}
								#$dSrc='-' unless $dSrc;
								my $dSrc='#'.$anaID;
								$dSrc.='.'.$1 if $pepData=~/SOURCE_RANK=(\d+)/;
								$dataSrc{$dSrc}=$dSrc;
								$dataSources{$dSrc}=$dSrc;
								$anaList{$anaID}=1;
							}
							$scoreStrg='-';
							if ($scoreMin && $scoreMax) {
								if ($scoreMin == $scoreMax) {$scoreStrg=1*(sprintf '%.2f',$scoreMin);}
								else {$scoreStrg='['.(1*(sprintf '%.2f',$scoreMin)).'-'.(1*(sprintf '%.2f',$scoreMax)).']';}
							}
							if ($quantifModID) {
								$prsStrg='-';
								if ($prsMin && $prsMax) {
									if ($prsMin == $prsMax) {$prsStrg=1*(sprintf '%.1f',$prsMin);}
									else {$prsStrg='['.(1*(sprintf '%.1f',$prsMin)).'-'.(1*(sprintf '%.1f',$prsMax)).']';}
								}
							}
							my $pepUsed=($pepExcluded)? 0 : 1;
							my @srcList=keys %dataSrc;
							my %srcInfo=(NAME=>$dataSrcAlias{$repCode} || $repCode,DATA=>\@srcList);
							#if ($labelType eq 'FREE') { #}
							if ($isPepIntensity) {
								if ($pepExcluded && $algoVersion==2) { # v2: special case: same occurences are on different lines
									push @{$seenExcluded{$usedPos}{$pepCode}},[$pepIdKey,$sequenceBeg{$pepSeq}[0],$pepSeq,$varMod,$charge,$scoreStrg,\%srcInfo,2**$log2Value]; # (keys %dataSrc)[0]
								}
								else {
									@{$peptideSets{$usedPos}{$pepIdKey}{INFO}}=($sequenceBeg{$pepSeq}[0],$pepSeq,$varMod,$charge,$scoreStrg,\%srcInfo); # join('<BR>',sort keys %dataSrc)
									push @{$peptideSets{$usedPos}{$pepIdKey}{DATA}},(2**$log2Value,scalar (split(/[+_;\.=]/,$pepIdKey)),$pepUsed);
								}
							}
							else {
								@{$peptideSets{$usedPos}{$pepIdKey}{INFO}}=($sequenceBeg{$pepSeq}[0],$pepSeq,$varMod,$charge,$scoreStrg,\%srcInfo); # join('<BR>',sort keys %dataSrc)
								my $ratioVal=($ratioPos > 0)? 2**$log2Value : 2**(-$log2Value);
								push @{$peptideSets{$usedPos}{$pepIdKey}{DATA}},($ratioVal,(scalar (split(/[+_;\.=]/,$pepIdKey)))/2,$pepUsed);
							}
						}
					}
				}
				#$dataSources{'-'}='...'; # no datasource distinction for (Super/simple)Ratio
				if (scalar keys %seenExcluded) { # Label-free with SimpleRatio (algo v2 only)
					foreach my $statePos (keys %seenExcluded) {
						foreach my $pepCode (keys %{$seenExcluded{$statePos}}) {
							my (@pepID,$seqBeg,$pepSeq,$varMod,$charge,@scores,%dataSrc,$meanValue,$occ);
							foreach my $refPep (@{$seenExcluded{$statePos}{$pepCode}}) {
								push @pepID,$refPep->[0];
								($seqBeg,$pepSeq,$varMod,$charge)=@{$refPep}[1..4];
								push @scores,$refPep->[5] if $refPep->[5] ne '-';
								foreach my $src (keys %{$refPep->[6]}) {$dataSrc{$src}=1;}
								$meanValue+=$refPep->[7];
								$occ++;
							}
							my $pepIdKey=join('_',sort{$a<=>$b} @pepID);
							my $scoreStrg;
							my $numScores=scalar @scores;
							if ($numScores==0) {$scoreStrg='-';}
							elsif ($numScores==1) {$scoreStrg=$scores[0];}
							else {
								my ($scoreMin,$scoreMax)=(sort{$a<=>$b}@scores)[0,-1];
								$scoreStrg='['.(1*(sprintf '%.2f',$scoreMin)).'-'.(1*(sprintf '%.2f',$scoreMax)).']';
							}
							my @srcList=keys %dataSrc;
							@{$peptideSets{$statePos}{$pepIdKey}{INFO}}=($seqBeg,$pepSeq,$varMod,$charge,$scoreStrg,{NAME=>'',DATA=>\@srcList}); # join('<BR>',sort keys %dataSrc)
							push @{$peptideSets{$statePos}{$pepIdKey}{DATA}},($meanValue/$occ,$occ,0);
						}
					}
				}

			}
		}
	}

	##########################################
	####<Raw data (former Infinite ratio)>####
	##########################################
	if ($ratioType eq 'Ratio') {
		###<1: Assume data is in table.txt (new quantif procedure: runXICProtQuantification v2.2.3+)
		###open(DATA,"$dataDir/table.txt") || die "$dataDir<BR>ERROR: $!";
		###my $prevProtID=0;
		###while(<DATA>) {
		###	next if $.==1;
		###	chomp;
		###	my ($mpID,$pepSeq,$varMod,$charge,$pepIdKey,$identifier,$validity,@values)=split(/\t/,$_);
		###	my ($pID)=($mpID=~/^(\d+)/);
		###	last if ($pID != $protID && $prevProtID==$protID); # modProtID not always in continuous order (but protID is!)
		###	$prevProtID=$pID;
		###	next if $mpID ne $modProtID;
		###	if ($varMod) {
		###		if ($varMod=~/^&/) { # var mod code
		###			$varMod=&decodeVarMod($dbh,$pepSeq,$varMod,\%varModName);
		###		}
		###		else { # var mod text
		###			$varMod=~s/^ \+ //;
		###		}
		###		#$varMod=' + '.$varMod;
		###	}
		###	else {$varMod='';}
		###	unless ($sequenceBeg{$pepSeq}) {
		###		my ($pep1ID)=$pepIdKey=~/^(\d+)/;
		###		$sthPB->execute($pep1ID);
		###		while (my ($beg)=$sthPB->fetchrow_array) {
		###			push @{$sequenceBeg{$pepSeq}},$beg;
		###		}
		###	}
		###	my ($dataSrc,$scoreMin,$scoreMax,$prsMin,$prsMax,$rtMin,$rtMax);
		###	foreach my $pepID (split(/\s*;\s*/,$pepIdKey)) {
		###		$sthPD->execute($pepID);
		###		my ($score,$pepData,$anaID,$rtApex,$qNum,$rk,$tgPos)=$sthPD->fetchrow_array;
		###		$score=0 unless $score;
		###		$pepData='' unless $pepData;
		###		#$scoreStrg.='/' if $scoreStrg;
		###		#$scoreStrg.=($score)? $score : '-';
		###		if ($score) {
		###			$scoreMin=$score if (!$scoreMin || $scoreMin > $score);
		###			$scoreMax=$score if (!$scoreMax || $scoreMax < $score);
		###		}
		###		unless ($dataSrc) {
		###			$dataSrc='#'.$anaID;
		###			$dataSrc.='.'.$1 if $pepData=~/SOURCE_RANK=(\d+)/;
		###		}
		###		if ($quantifModID && $pepData=~/PRS=\d;([^;]+);/) {
		###			my $prs=$1;
		###			$prsMin=$prs if (!$prsMin || $prsMin > $prs);
		###			$prsMax=$prs if (!$prsMax || $prsMax < $prs);
		###		}
		###		if ($rtApex) {
		###			$rtMin=$rtApex if (!$rtMin || $rtMin > $rtApex);
		###			$rtMax=$rtApex if (!$rtMax || $rtMax < $rtApex);
		###		}
		###		#$dataSrc=$dSrc || '-'; # Not good for LABEL-FREE!!!! (!= for refState & testState)
		###		$anaList{$anaID}=1;
		###	}
		###	#$dataSrc='' unless $dataSrc;
		###	$dataSources{$dataSrc}=$dataSrc;
		###	my $uniquePep=$pepSeq.$varMod.'_'.$charge;
		###	foreach my $xic (@values) {
		###		$pepOccurence{$selRatioPos}{$uniquePep}++ if $xic ne 'NA';
		###	}
		###	my $pepCode=$pepSeq.':'.$varMod.'_'.$charge;
		###	$pepCode.="_$dataSrc" if $labelType ne 'FREE';
		###	@{$rawPeptideSets{$selRatioPos}{$pepCode}{INFO}}=($sequenceBeg{$pepSeq}[0],$pepSeq,$varMod,$charge,$uniquePep);
		###	push @{$rawPeptideSets{$selRatioPos}{$pepCode}{DATA}{$refStatePos{$selRatioPos}}},[$scoreMin,$values[$refStatePos{$selRatioPos}-1],$rtMax,$dataSrc,undef,undef,undef];
		###	push @{$rawPeptideSets{$selRatioPos}{$pepCode}{DATA}{$testStatePos{$selRatioPos}}},[$scoreMax,$values[$testStatePos{$selRatioPos}-1],$rtMin,$dataSrc,undef,undef,undef];
		###}


		my %anaTg2StatePos;
		if ($designID) {
			foreach my $expCondID (keys %condToObs) {
next if ($condToState{$expCondID} != $testStatePos{$selRatioPos} && $condToState{$expCondID} != $refStatePos{$selRatioPos});
				foreach my $replic (split(/\./,$condToObs{$expCondID})) {
					foreach my $fraction (split(/\+/,$replic)) { # there should be only 1 fraction in label-free...
						my ($obsID,$pepQuantID,$anaID,$tgPos)=split(':',$fraction);
						$tgPos=0 unless $tgPos;
						$anaTg2StatePos{"$anaID.$tgPos"}=$condToState{$expCondID}; # absolute pos in all states list
					}
				}
			}
		}
		else { # internal
			my ($anaID)=$dbh->selectrow_array("SELECT ID_ANALYSIS FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");

			my ($testChanNum)=$labelingInfo{STATES}[$testStatePos{$selRatioPos}-1]=~/,(\d+)$/;
			my ($refChanNum)=$labelingInfo{STATES}[$refStatePos{$selRatioPos}-1]=~/,(\d+)$/;
			$anaTg2StatePos{"$anaID.$testChanNum"}=$testStatePos{$selRatioPos};
			$anaTg2StatePos{"$anaID.$refChanNum"}=$refStatePos{$selRatioPos};
#			my ($parQuantiAnnot,$parQuantifID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,PQ.ID_PARENT_QUANTIFICATION FROM QUANTIFICATION Q,PARENT_QUANTIFICATION PQ WHERE Q.ID_QUANTIFICATION=PQ.ID_PARENT_QUANTIFICATION AND PQ.ID_QUANTIFICATION=$selQuantifID LIMIT 0,1"); # there should be only one
#			my ($labelStrg,@labelInfo)=split('::',$parQuantiAnnot);
#			foreach my $infoStrg (@labelInfo) {
#				last if $infoStrg !~ /^\d/; # entering MassChroQ params if from MassChroQ XIC
#				my ($chanNum,$chanName,$labelDataStrg)=split(';',$infoStrg);
#print "$chanNum,$chanName,$labelDataStrg<BR>\n";
#				#$anaTg2StatePos{"$anaID.$chanNum"}=$testPos;
#
#			}
			#my ($testPos,$refPos)=split(/\//,$labelingInfo{'RATIOS'}[$selRatioPos-1]);
			#$anaTg2StatePos{"$anaID.$testPos"}=$testPos;
			#$anaTg2StatePos{"$anaID.$refPos"}=$refPos;
		}

		###my $sthXic=$dbh->prepare("SELECT QUANTIF_VALUE,TARGET_POS FROM PEPTIDE_QUANTIFICATION WHERE ID_PEPTIDE=? AND ID_QUANTIFICATION IN ($pepQuantifIDs) LIMIT 0,1");
		my $sthTgPos=$dbh->prepare("SELECT TARGET_POS FROM PEPTIDE_QUANTIFICATION WHERE ID_PEPTIDE=? AND ID_QUANTIFICATION IN ($pepQuantifIDs)");

		my %targetPos2StatePos;
		open(DATA,"$dataDir/table.txt") || die "$dataDir<BR>ERROR: $!";
		my $prevProtID=0;
		while(<DATA>) {
			next if $.==1;
			chomp;
			my ($mpID,$pepSeq,$varMod,$charge,$pepIdKey,$identifier,$validity,@values)=split(/\t/,$_);
			my ($pID)=($mpID=~/^(\d+)/);
			last if ($pID != $protID && $prevProtID==$protID); # modProtID not always in continuous order (but protID is!)
			$prevProtID=$pID;
			next if $mpID ne $modProtID;
			if ($varMod) {
				if ($varMod=~/^&/) { # var mod code
					$varMod=&promsMod::decodeVarMod($dbh,$pepSeq,$varMod,\%varModName);
				}
				else { # var mod text
					$varMod=~s/^ \+ //;
				}
				#$varMod=' + '.$varMod;
			}
			else {$varMod='';}

			my (%multiPepData,%xicValue,%qSetList,$refAnaID);
			foreach my $statePos (values %anaTg2StatePos) {
				$xicValue{$statePos}=$values[$statePos-1];
			}
			my @peptideList=split(/\s*;\s*/,$pepIdKey);
			unless ($sequenceBeg{$pepSeq}) {
				#$sthPB->execute($peptideList[0]);
				#while (my ($beg)=$sthPB->fetchrow_array) {
				#	push @{$sequenceBeg{$pepSeq}},$beg;
				#}
@{$sequenceBeg{$pepSeq}}=@{$protPeptideData{$peptideList[0]}{BEG}};
			}
			my $uniquePep=$pepSeq.$varMod.'_'.$charge;
			#foreach my $xic (@values) {
			#	$pepOccurence{$selRatioPos}{$uniquePep}++ if $xic ne 'NA';
			#}
			my $basePepCode=$pepSeq.':'.$varMod.'_'.$charge;
			#my @targetList;

			MAIN_PEP:foreach my $pepIdStrg (@peptideList) { # No longer checks
				my $statePos;
				###my $xicValue=0;
				foreach my $pepID (split(/\+/,$pepIdStrg)) { # compatible with XIC summing
					#$sthPD->execute($pepID);
					#my ($score,$pepData,$anaID,$rtApex,$qNum,$rank)=$sthPD->fetchrow_array;
					#$score=0 unless $score;
					#$pepData='' unless $pepData;
my ($score,$pepData,$anaID,$qNum,$rank)=@{$protPeptideData{$pepID}{INFO}};
my $rtApex=$peptideApexRT{$pepID};
					###$sthXic->execute($pepID);
					###my ($xicVal,$tgPos)=$sthXic->fetchrow_array;
					###next unless $xicVal;
					###$tgPos=0 unless $tgPos;
					###next MAIN_PEP unless $anaTg2StatePos{"$anaID.$tgPos"}; # peptide does not belong to test nor ref conditions (multiple ratios in same quantif)
					###$statePos=$anaTg2StatePos{"$anaID.$tgPos"};
					my %matchedTgPos;
					$sthTgPos->execute($pepID);
					while (my ($tgPos)=$sthTgPos->fetchrow_array) { # multiple matches for same peptideID for iTRAQ. Only 1 for other (non-)labelling methods
						$tgPos=0 unless $tgPos;
						$matchedTgPos{$anaTg2StatePos{"$anaID.$tgPos"}}=1 if $anaTg2StatePos{"$anaID.$tgPos"};
					}
					next MAIN_PEP unless scalar keys %matchedTgPos; # peptide does not belong to test nor ref conditions (multiple ratios in same quantif)
					$anaList{$anaID}=1;
					$refAnaID=$anaID;
					###$xicValue{$statePos}+=$xicVal; # recompute sum
					my $dataSrc='#'.$anaID;
					$dataSrc.='.'.$1 if $pepData=~/SOURCE_RANK=(\d+)/;
					my $qSet;
					if ($labelType eq 'FREE') {$qSet='-';}
					else {
						$qSet=($pepData=~/QSET=([^#]+)/)? $1 : 0;
						$qSetList{$qSet}=1;
					}
					my $prsStrg=$1 if ($quantifModID && $pepData=~/PRS=([^#]+)/);
					#$pepCode.="_$dataSrc" if $labelType ne 'FREE';
					foreach my $statePos (keys %matchedTgPos) {
						@{$multiPepData{$statePos}{$qSet}}=($pepID,$score,$dataSrc,$rtApex,$prsStrg,$qNum,$rank) if (!$multiPepData{$statePos} || !$multiPepData{$statePos}{$qSet} || $multiPepData{$statePos}{$qSet}[1]<$score); # qset seen 1+ for same state!!! => prevents overwrite of best scoring pep in qset (may not be the best fix) PP 01/02/16
#print "+$statePos:$qSet: @{$multiPepData{$statePos}{$qSet}}<BR>\n";
					}
				}
				$pepOccurence{$selRatioPos}{$uniquePep}++;
			}
			my $quantifSet=($labelType eq 'FREE')? '-' : 'A'.$refAnaID.'_'.join('+',sort{$a<=>$b} keys %qSetList);
			my $pepCode=$basePepCode;
			@{$rawPeptideSets{$selRatioPos}{$pepCode}{INFO}}=($sequenceBeg{$pepSeq}[0],$pepSeq,$varMod,$charge,$uniquePep);
				#push @{$rawPeptideSets{$selRatioPos}{$pepCode}{DATA}{$qSet}{$ana2StatePos{$anaID}}},[$score,(sprintf "%.0f",$xicValue),$rtApex,$dataSrc,$prsStrg,$pepID,$qNum,$rank];
			foreach my $statePos (keys %multiPepData) {
				push @{$rawPeptideSets{$selRatioPos}{$pepCode}{DATA}{$quantifSet}{$statePos}},[$xicValue{$statePos},$multiPepData{$statePos}]; # push because of label-free (only 1 index for labeled because of dataSrc in pepCode)
			}
		}
		close DATA;
		###$sthXic->finish;
		$sthTgPos->finish;
	}
	else { # (Super/Simple)Ratio
		my (%statePos2RatioPos); # %selExp,
		foreach my $ratioPos (keys %protRatio) { # in case different ratios use the same (eg ref) state
			if ($isPepIntensity) { # $labelType eq 'FREE'
				foreach my $statePos (keys %protMean) {@{$statePos2RatioPos{$statePos}}=($ratioPos);}
			}
			else {
##$selExp{chr(64+$ratioPos)}=$ratioPos; # not for iTRAQ (exp='A' always)
				push @{$statePos2RatioPos{$testStatePos{$ratioPos}}},$ratioPos;
				push @{$statePos2RatioPos{$refStatePos{$ratioPos}}},$ratioPos;
			}
			%{$rawPeptideSets{$ratioPos}}=();
		}
		open(DATA,"$dataDir/table.txt") || die "$dataDir<BR>ERROR: $!";
		#my $selExp=chr(64+$selRatioPos); # cannot be 2ndary ratio!!!
		my $prevProtID=0;
		while(<DATA>) {
			next if $.==1;
			chomp;
			my ($mpID,$pepCode,$state,@data)=split(/\t/,$_);
#next unless $usedStates{$state}; # this state is not used for selected ratio
			my ($pID)=($mpID=~/^(\d+)/);
			last if ($pID != $protID && $prevProtID==$protID); # modProtID not always in continuous order (but protID is!)
			$prevProtID=$pID;
			next if $mpID ne $modProtID;

my $quantifSet=($labelType eq 'FREE')? '-' : $data[-5]; # count from end because +/- techRep
if ($ratioType eq 'SuperRatio' && $state eq 'State1') { # restrict State1 to selected ratio context (anaIDs)
	my ($anaID)=$quantifSet=~/^A(\d+)/;
	next unless $anaList{$anaID};
}
			my ($statePos)=($state=~/(\d+)/);
#if ($labelType eq 'FREE') {next unless $protMean{$statePos};} # state filter
#else {next unless $selExp{$data[-6]};} # Experiment: index from end because +/- techRep
	###next if ($labelType ne 'FREE' && $labelType ne 'ITRAQ' && $data[-6] ne 'NA' && !$selExp{$data[-6]}); # Experiment: index from end because +/- techRep
	###my $quantifSet=($labelType eq 'FREE')? '-' : $data[-5]; # count from end because +/- techRep
			my $pepIdKey=$data[-4]; # count from end because +/- techRep
			my ($pep1ID)=$pepIdKey=~/^(\d+)/; # count from end because +/- techRep
			#next unless ($exp eq 'NA' || $exp eq $selExp);
			#my $xicValue=sprintf "%.0f",$data[-1]; # count from end because +/- techRep
			my $xicValue=$data[-1]; # count from end because +/- techRep
			my ($seqModCode,$charge,$dataSrc)=split(/_/,$pepCode);
			$pepCode=~s/_[^_]+\Z//;# if $labelType eq 'FREE'; # remove dataSource for label-free quantif
			my $uniquePep=$seqModCode.'_'.$charge;

			foreach my $ratioPos (@{$statePos2RatioPos{$statePos}}) { # selected ratio(s) that use this state
				$pepOccurence{$ratioPos}{$uniquePep}++;
				unless ($rawPeptideSets{$ratioPos}{$pepCode}) { # 1st time seen
					$dataSrc='-' unless $dataSrc;
					my ($pepSeq,$varMod);
					if ($seqModCode=~/(\w+)&(.+)/) {
						$pepSeq=$1;
						$varMod=&promsMod::decodeVarMod($dbh,$pepSeq,$2,\%varModName);
					}
					else {
						$pepSeq=$seqModCode;
						$varMod='';
					}
					unless ($sequenceBeg{$pepSeq}) {
						#$sthPB->execute($pep1ID);
						#while (my ($beg)=$sthPB->fetchrow_array) {
						#	push @{$sequenceBeg{$pepSeq}},$beg;
						#}
@{$sequenceBeg{$pepSeq}}=@{$protPeptideData{$pep1ID}{BEG}};
					}
					@{$rawPeptideSets{$ratioPos}{$pepCode}{INFO}}=($sequenceBeg{$pepSeq}[0],$pepSeq,$varMod,$charge,$uniquePep);
				}
				$dataSources{$dataSrc}=$dataSrc; # old quantif with true data source file in table.xtx

				my %multiPepData;
				foreach my $pepIdStrg (split(/\+/,$pepIdKey)) {
					my ($pepID)=$pepIdStrg=~/^(\d+)/; # because of early quantifs: pepID1;peppID2;...
					#$sthPD->execute($pepID);
					#my ($score,$pepData,$anaID,$rtApex,$qNum,$rank)=$sthPD->fetchrow_array;
					#$score=0 unless $score;
					#$pepData='' unless $pepData;
my ($score,$pepData,$anaID,$qNum,$rank)=@{$protPeptideData{$pepID}{INFO}};
my $rtApex=$peptideApexRT{$pepID};
					$anaList{$anaID}=1;
					$dataSrc='#'.$anaID;
					$dataSrc.='.'.$1 if $pepData=~/SOURCE_RANK=(\d+)/;
					my $qSet=($labelType eq 'FREE')? '-' : ($pepData=~/QSET=([^#]+)/)? $anaID.':'.$1 : $anaID.':0'; # compatible with merged files?????
					my $prsStrg=$1 if ($quantifModID && $pepData=~/PRS=([^#]+)/);
					@{$multiPepData{$qSet}}=($pepID,$score,$dataSrc,$rtApex,$prsStrg,$qNum,$rank);
				}
				#push @{$rawPeptideSets{$ratioPos}{$pepCode}{DATA}{$quantifSet}{$statePos}},[$score,$xicValue,$rtApex,$dataSrc,$prsStrg,$pep1ID,$qNum,$rank]; # push because of label-free (only 1 index for labeled because of dataSrc in pepCode)
				push @{$rawPeptideSets{$ratioPos}{$pepCode}{DATA}{$quantifSet}{$statePos}},[$xicValue,\%multiPepData]; # push because of label-free (only 1 index for labeled because of dataSrc in pepCode)
			}
		}
		close DATA;
	}
	#$sthPB->finish;
	#$sthPD->finish;

	##<Converting dataSources>##
	my $numSources=scalar keys %dataSources;

	my %analysisSrc;
	my $sthDS=$dbh->prepare("SELECT VALID_STATUS,WIFF_FILE FROM ANALYSIS WHERE ID_ANALYSIS=?");
	foreach my $anaID (keys %anaList) {
		$sthDS->execute($anaID);
		my ($validStatus,$wiffFile)=$sthDS->fetchrow_array;
		my $anaDir=($validStatus==2)? "$promsPath{peptide}/proj_$projectID/ana_$anaID" : "$promsPath{valid}/ana_$anaID";
		if (-e "$anaDir/mergedFiles.txt") {
			open(MERGE,"$anaDir/mergedFiles.txt");
			while (<MERGE>) {
				my ($srcRk,$srcName,$srcID)=split(/\t/,$_);
				$analysisSrc{'#'.$anaID.'.'.$srcRk}=$wiffFile.'>'.$srcName;
			}
			close MERGE;
		}
		else {$analysisSrc{'#'.$anaID}=$wiffFile;}
	}
	$sthDS->finish;
	foreach my $srcKey (keys %analysisSrc) {
		$dataSources{$srcKey}=$analysisSrc{$srcKey} if $srcKey=~/^#/;
	}

	$dbh->disconnect;

	######################
	####<Starting HTML>###
	######################
	#my ($lightColor,$darkColor)=&promsConfig::getRowColors;
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);
	my ($biasStrg1,$biasStrg2)=($ratioType=~/S\w+Ratio/)? ('','') : ($labelingInfo{BIAS_CORRECTION}[0] eq 'TRUE')? ('<SUP>*</SUP>','<I>Intensities after data normalization</I>') : ('','<BR>');
	my $anaStrg=join(',',keys %anaList);

	if ($is2ndaryRatio) { # SuperRatio
		#my $numNormalRatios=0;
		#foreach my $ratioPos (keys %protRatio) {$numNormalRatios++ if $normalRatio{$ratioPos};}
		my $protRatioStrg;
		if ($trueProtRatio < 1) {$protRatioStrg=($numNormalRatios==2)? '1/'.sprintf "%.2f",1/$trueProtRatio : '1/&infin;';} #'∞'
		else {$protRatioStrg=($numNormalRatios==2)? sprintf "%.2f",$trueProtRatio : '&infin;';}
		print qq
|<TABLE><TR><TD nowrap>
<FONT class="title2">Peptide data for <A href="javascript:sequenceView($protID,'$anaStrg')">$protAlias$modStrg</A> in $trueTestStateName/$trueRefStateName (ratio=$protRatioStrg):</FONT>
</TD></TR></TABLE>
|;
	}
	my %targetList;
	my $lastPepStart=0;
	foreach my $ratioPos (sort{if ($selRatioPos > 0) {abs($b)<=>abs($a)} else {abs($a)<=>abs($b)}} keys %protRatio) { # test then ref
		my $protRatioStrg;
		if (!$protRatio{$ratioPos}) {$protRatioStrg='unknown <FONT class="title3">[All peptides excluded]</FONT>';} # In case global CVout excludes all peptides of this ratio
		elsif ($protRatio{$ratioPos} < 1) {$protRatioStrg=($normalRatio{$ratioPos})? '1/'.sprintf "%.2f",1/$protRatio{$ratioPos} : '1/&infin;';} #'∞'
		else {$protRatioStrg=($normalRatio{$ratioPos})? sprintf "%.2f",$protRatio{$ratioPos} : '&infin;';}
		if ($is2ndaryRatio) {
			print qq
|<HR width=100%>
<TABLE><TR><TD nowrap>
&nbsp;&nbsp;&nbsp;&bull;<FONT class="title2">Primary ratio: $testStateName{$ratioPos}/$refStateName{$ratioPos} (ratio=$protRatioStrg):</FONT>
|;
		}
		else {
			print qq
|<TABLE><TR><TD nowrap>
<FONT class="title2">Peptide data for <A href="javascript:sequenceView($protID,'$anaStrg')">$protAlias$modStrg</A> in $testStateName{$ratioPos}/$refStateName{$ratioPos} (ratio=$protRatioStrg):</FONT>
|;
		}
		print '&nbsp;&nbsp;<FONT class="title3">View:';
		my $rawDataVis='';
		if ($normalRatio{$ratioPos}) {
			$rawDataVis='none';
			print qq
|</FONT><SELECT class="title3" onchange="var divList=['pepPlot:$selQuantifID:$ratioPos','pepList:$selQuantifID:$ratioPos','rawData:$selQuantifID:$ratioPos']; for (var i=0; i<3; i++) {var vis=(divList[i]==this.value)? '' : 'none'; document.getElementById(divList[i]).style.display=vis;}">
<OPTION value="pepPlot:$selQuantifID:$ratioPos">Graphical</OPTION><OPTION value="pepList:$selQuantifID:$ratioPos">List</OPTION><OPTION value="rawData:$selQuantifID:$ratioPos">Raw data</OPTION>
</SELECT>
|;
		}
		else {print " Raw data</FONT>";}
		print qq
|&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'"/>
</TD></TR></TABLE>
<DIV id="rawData:$selQuantifID:$ratioPos" style="display:$rawDataVis">
|;
		my ($numColumns,$scColClass,$srcColStrg)=($isPepIntensity)? (11,'bBorder','') : ($numSources > 1)? (12,'rbBorder','<TH class="bBorder">&nbsp;Source&nbsp;</TH>') : (11,'bBorder',''); # $labelType eq 'FREE'
		#my $rtColStrg=($labelType eq 'FREE')? '' : "<TH class=\"$rtColClass\" nowrap>&nbsp;Ret. time&nbsp;</TH>";
		#if ($labelType eq 'ITRAQ') { # no need for RT column
		#	$rtColStrg='';
		#	$numColumns--;
		#}
		if ($ratioType eq 'Ratio' && !scalar keys %{$rawPeptideSets{$ratioPos}}) { # TODO: temp until update
			print "<BR>&nbsp;&nbsp;<FONT class=\"title3\">This quantification is too old for infinite ratios display.<BR>&nbsp;&nbsp;Delete it and run it again if necessary.<BR><BR><BR>\n";
		}
		else {
			my $scLabelStrg=($labelType=~/ITRAQ|TMT/)? 'Score' : 'Scores';
			print qq
|<TABLE cellspacing=0>
<TR bgcolor=$darkColor>
<TH class="rbBorder">#</TH><TH class="rbBorder">&nbsp;Pos.&nbsp;<TH class="rbBorder">Peptides</TH><TH class="rbBorder">&nbsp;Charge&nbsp;</TH><TH class="rbBorder">&nbsp;Occur.&nbsp;</TH><TH class="rbBorder" colspan=3>&nbsp;XIC values&nbsp;</TH><TH class="$scColClass" colspan=3>&nbsp;$scLabelStrg&nbsp;</TH>$srcColStrg
</TR>
|;
			my $bgColor=$lightColor;
			my $numPep=0;
			my $peptide;
			my $prevPeptide='';
			my $prevUniquePep='';
			my $prevBeg=0;
			foreach my $pepCode (sort{$rawPeptideSets{$ratioPos}{$a}{INFO}[0]<=>$rawPeptideSets{$ratioPos}{$b}{INFO}[0] || $rawPeptideSets{$ratioPos}{$a}{INFO}[1] cmp $rawPeptideSets{$ratioPos}{$b}{INFO}[1] || $rawPeptideSets{$ratioPos}{$a}{INFO}[2] cmp $rawPeptideSets{$ratioPos}{$b}{INFO}[2] || $rawPeptideSets{$ratioPos}{$a}{INFO}[3]<=>$rawPeptideSets{$ratioPos}{$b}{INFO}[3]} keys %{$rawPeptideSets{$ratioPos}}) {
				my $varMod=($rawPeptideSets{$ratioPos}{$pepCode}{INFO}[2])? ' + '.$rawPeptideSets{$ratioPos}{$pepCode}{INFO}[2] : '';
				my $uniquePep=$rawPeptideSets{$ratioPos}{$pepCode}{INFO}[4];
				my $numPepStrg='';
				$peptide=$rawPeptideSets{$ratioPos}{$pepCode}{INFO}[1].$varMod;
				foreach my $quantifSet (sort{&promsMod::sortSmart($rawPeptideSets{$ratioPos}{$pepCode}{DATA}{$a}{$testStatePos{$ratioPos}}[0][3],$rawPeptideSets{$ratioPos}{$pepCode}{DATA}{$b}{$testStatePos{$ratioPos}}[0][3])} keys %{$rawPeptideSets{$ratioPos}{$pepCode}{DATA}}) {
					if ($dispPepType eq 'distinct') {
						if ($peptide ne $prevPeptide) {
							$numPep++;
							$numPepStrg=$numPep;
						}
					}
					else {
						$numPep++;
						$numPepStrg=$numPep;
					}
					my $peptideStrg=($peptide ne $prevPeptide)? $rawPeptideSets{$ratioPos}{$pepCode}{INFO}[1].$varMod : '';
					my $begStrg=($prevBeg != $rawPeptideSets{$ratioPos}{$pepCode}{INFO}[0])? "$rawPeptideSets{$ratioPos}{$pepCode}{INFO}[0]&nbsp;" : '';
					my ($chargeStrg,$occStrg)=($prevUniquePep eq $uniquePep)? ('','') : ($rawPeptideSets{$ratioPos}{$pepCode}{INFO}[3].'<SUP>+</SUP>','x'.$pepOccurence{$ratioPos}{$uniquePep});
					print qq
|<TR bgcolor="$bgColor" valign="top">
	<TD class="rBorder" align="right">&nbsp;$numPepStrg&nbsp;</TD><TD align="right" valign="top">$begStrg</TD><TH class="font11" align="left" valign="top" nowrap>$peptideStrg&nbsp;</TH>
	<TH valign="top">$chargeStrg</TH><TD align="center" valign="top">&nbsp;$occStrg&nbsp;</TD>
|;
					#my ($sumRetTime,$numRetTime)=(0,0); # labeled
					my %numPeps=($testStatePos{$ratioPos}=>0,$refStatePos{$ratioPos}=>0); # FREE
					my %numPeps2=($testStatePos{$ratioPos}=>0,$refStatePos{$ratioPos}=>0);
					my (%peptID,%sc,%dtsrc,%rt,%prs,%qNum,%rank);
					foreach my $statePos ($testStatePos{$ratioPos},$refStatePos{$ratioPos}) {
						#my ($xicValue,$retTime)=($rawPeptideSets{$pepCode}{DATA}{$statePos} && $rawPeptideSets{$pepCode}{DATA}{$statePos}[1])? ($rawPeptideSets{$pepCode}{DATA}{$statePos}[1],$rawPeptideSets{$pepCode}{DATA}{$statePos}[2]) : ('-','-');
						my $xicStrg='';
						#my $prsStrg; # labeled
#print "*$statePos: ",scalar @{$rawPeptideSets{$ratioPos}{$pepCode}{DATA}{$statePos}},"<BR>\n";
						foreach my $refData (@{$rawPeptideSets{$ratioPos}{$pepCode}{DATA}{$quantifSet}{$statePos}}) {
							if ($isPepIntensity || $labelType eq 'FREE') { # $labelType eq 'FREE'
								$xicStrg.='<BR>' if $xicStrg;
					#my $rt=($refData->[2])? sprintf "%.2f",$refData->[2] : '-';
					#my $dtSrc=($refData->[0])? $dataSources{$refData->[0]} : '-';
								#$xicStrg.="&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<TABLE cellspacing=0><TR><TD align=right nowrap><U>Ret. time:</U></TD><TD>$rt<TD></TR><TR><TD align=right><U>Source:</U></TD><TD>$dtSrc</TD></TR></TABLE>')\" onmouseout=\"popout()\">$refData->[1]</A>&nbsp;";
								$xicStrg.=($refData->[0] && $refData->[0] ne 'NA')? sprintf "&nbsp;%1.2e&nbsp;",$refData->[0] : '-';
								$numPeps{$statePos}++;
							}
							else {
								#if ($refData->[2]) {
								#	$sumRetTime+=$refData->[2];
								#	$numRetTime++;
								#}
								#($dataSrc,$prsStrg)=@{$refData}[3..4];
								#$xicStrg=sprintf "&nbsp;%1.2e&nbsp;",$refData->[1] if $refData->[1];
								if ($refData->[0] && $refData->[0] ne 'NA') { # NA if algo Ratio
									$xicStrg=sprintf "&nbsp;%1.2e&nbsp;",$refData->[0];
									#if (!$refData->[6]) {$xicStrg="&nbsp;<I>$xicValue</I>";} # !query number -> ghost
									#else {
									#	$xicStrg="&nbsp;<A id=\"pep_$refData->[5]\" href=\"javascript:drawSpectrum('pep_$refData->[5]','pep_$refData->[5]_$refData->[6]_$refData->[7]')\" onmouseover=\"popup('Click to display <B>Fragmentation Spectrum</B>.')\" onmouseout=\"popout()\">$xicValue</A>";
									#}
									#if ($prsStrg) {$xicStrg.=&phosphoRS::printIcon($prsStrg,{format=>'string'});}
									#else {$xicStrg.='&nbsp;'}
								}
							}
							#>Multi-peptide data
							foreach my $qSet (keys %{$refData->[1]}) {
#print "*$qSet: @{$refData->[1]{$qSet}}<BR>\n";
								#($peptID{$qSet}{$statePos},$sc{$qSet}{$statePos},$rt{$qSet}{$statePos},$prs{$qSet}{$statePos},$qNum{$qSet}{$statePos},$rank{$qSet}{$statePos})=@{$refData->[2]{$qSet}};
								push @{$peptID{$qSet}{$statePos}},$refData->[1]{$qSet}[0]; # push because of new Label-Free (but only 1 qSet per dataScr)
								push @{$sc{$qSet}{$statePos}},$refData->[1]{$qSet}[1];
								push @{$dtsrc{$qSet}{$statePos}},$refData->[1]{$qSet}[2];
								push @{$rt{$qSet}{$statePos}},$refData->[1]{$qSet}[3];
								push @{$prs{$qSet}{$statePos}},$refData->[1]{$qSet}[4];
								push @{$qNum{$qSet}{$statePos}},$refData->[1]{$qSet}[5];
								push @{$rank{$qSet}{$statePos}},$refData->[1]{$qSet}[6];
								$numPeps2{$statePos}++;
							}
						}
						$xicStrg='&nbsp;-&nbsp;' unless $xicStrg;
						my $sepStrg;
						if ($isPepIntensity || $labelType eq 'FREE') { # $labelType eq 'FREE' for 'Ratio'
							my $numSep=(sort{$b<=>$a} values %numPeps)[0];
							$sepStrg='|<BR>' x $numSep;
						}
						else {$sepStrg='/';}
						#$retTimeStrg.=($retTime)? sprintf "%.2f",$retTime : '-';
						#$retTimeStrg.=($retTime)? &phosphoRS::printIcon($prsStrg,{format=>'string',text=>'&nbsp;'.(sprintf "%.2f",$retTime).'&nbsp;'}) : '-';
						#$retTimeStrg.=&phosphoRS::printIcon($prsStrg,{format=>'string'}) if $prsStrg;
						my $align;
						if ($statePos==$testStatePos{$ratioPos}) {
							$align='right';
							#$retTimeStrg.=($ratioType eq 'Ratio')? '-' : ' / '; # labeled
							#$retTimeStrg.='/';
						}
						else {
							$align='left';
							print "<TD>$sepStrg</TD>";
						}
						print "<TD align=\"$align\">$xicStrg</TD>";
					}

					#>Other columns
my $dataSrc='-';
my $align='right';
my $numSep2=(sort{$b<=>$a} values %numPeps2)[0];
foreach my $statePos ($testStatePos{$ratioPos},$refStatePos{$ratioPos}) {
	my $scoreStrg;
	my $plusStrg='';
	foreach my $qSet (sort keys %peptID) {
		if ($peptID{$qSet}{$statePos}) {
			foreach my $i (0..$#{$peptID{$qSet}{$statePos}}) { # i > 0 for Label-Free only!
				$scoreStrg.="<BR>" if $scoreStrg;
				($dataSrc,my $dtSrcStrg)=($dtsrc{$qSet}{$statePos}[$i])? ($dtsrc{$qSet}{$statePos}[$i],$dataSources{$dtsrc{$qSet}{$statePos}[$i]}) : ('-','-'); # changes if label-free (unique for labeled)
				my $dataPopStrg='';
				$dataPopStrg='<TABLE cellspacing=0>';
				if ($rt{$qSet}{$statePos}[$i]) {
					my $ret=1*(sprintf "%.2f",$rt{$qSet}{$statePos}[$i]);
					$dataPopStrg.="<TR><TD align=right nowrap><U>Ret. time:</U></TD><TD>$ret min.<TD></TR>";
				}
				$dataPopStrg.="<TR><TD align=right><U>Source:</U></TD><TD>$dtSrcStrg</TD></TR></TABLE>";
				my ($score,$linkStrg,$linkPopStrg)=('?','void(null)','');
				if ($sc{$qSet}{$statePos}[$i]) {
					$score=1 * (sprintf '%.1f',$sc{$qSet}{$statePos}[$i]);
					$linkStrg="drawSpectrum('pep_$peptID{$qSet}{$statePos}[$i]','pep_$peptID{$qSet}{$statePos}[$i]_$qNum{$qSet}{$statePos}[$i]_$rank{$qSet}{$statePos}[$i]')";
					$linkPopStrg='Click&nbsp;to&nbsp;display&nbsp;<B>Fragmentation&nbsp;Spectrum</B>';
				}
				$scoreStrg.="&nbsp;$plusStrg<A id=\"pep_$peptID{$qSet}{$statePos}[$i]\" href=\"javascript:$linkStrg\" onmouseover=\"popup('$dataPopStrg$linkPopStrg')\" onmouseout=\"popout()\">$score</A>";
				$scoreStrg.=($prs{$qSet}{$statePos}[$i])? &phosphoRS::printIcon($prs{$qSet}{$statePos}[$i],{format=>'string'}) : '&nbsp;';
			}
		}
		else {
			$scoreStrg.="<BR>" if $scoreStrg;
			$scoreStrg.='&nbsp;-&nbsp;';
		}
		$plusStrg='+';
	}
	if ($labelType=~/ITRAQ|TMT/) {
		print "<TD colspan=3 align=\"center\">$scoreStrg</TD>"; # same score for all states
		last;
	}
	else {
		print "<TD align=\"$align\">$scoreStrg</TD>";
		if ($statePos==$testStatePos{$ratioPos}) {
			my $sepStrg=($isPepIntensity || $labelType eq 'FREE')? '|<BR>' : '/<BR>'; # $labelType eq 'FREE' for Ratio
			print '<TD>',$sepStrg x $numSep2,'</TD>';
		}
		$align='left';
	}
}

					print "<TD>&nbsp;$dataSources{$dataSrc}&nbsp;</TD>" if (!$isPepIntensity && $numSources > 1); # $labelType ne 'FREE'
					print "</TR>\n";
					$prevBeg=$rawPeptideSets{$ratioPos}{$pepCode}{INFO}[0];
					$prevPeptide=$peptide;
					$prevUniquePep=$uniquePep;
					$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
				}
			}
			print "<TR bgcolor=\"$darkColor\"><TD colspan=$numColumns></TD></TR>\n" if $bgColor eq $darkColor;
			print "</TABLE><BR>\n";
		}
		print "</DIV>\n";

		####<Normal ratio>####
		if ($normalRatio{$ratioPos}) {
			my ($refClass,$sourceColStrg,$colspan)=($numSources > 1)? ('rbBorder','<TH class="bBorder">&nbsp;Source&nbsp;</TH>',8) : ('bBorder','',7);
			my ($labelFreeDivStrg,%titleStrg,$valueHeader);
			if ($isPepIntensity) { # $labelType eq 'FREE' && $ratioType ne 'Ratio'
				$valueHeader='Mean';
				foreach my $targetPos (sort{if ($selRatioPos > 0) {$b<=>$a} else {$a<=>$b}} keys %protMean) {
					push @{$targetList{$ratioPos}},$targetPos;
					my $stateName=($labelFreeDivStrg)? $trueRefStateName : $trueTestStateName;
					$titleStrg{$targetPos}="&nbsp;&nbsp;&nbsp;&bull;<FONT class=\"title2\">Mean values for $stateName:</FONT><BR>";
					$labelFreeDivStrg.="$titleStrg{$targetPos}<DIV id=\"pepPlot:$selQuantifID:$ratioPos:$targetPos\"></DIV>\n";
				}
			}
			else {
				$valueHeader='Ratio';
				@{$targetList{$ratioPos}}=($ratioPos);
				$labelFreeDivStrg=$titleStrg{$ratioPos}='';
			}
			print qq
|<DIV id="pepPlot:$selQuantifID:$ratioPos">$labelFreeDivStrg</DIV>
<DIV id="pepList:$selQuantifID:$ratioPos" style="display:none">
|;
			#foreach my $targetPos (@{$targetList{$ratioPos}}) { #}
			foreach my $tpIdx (0..$#{$targetList{$ratioPos}}) {
				my $targetPos=$targetList{$ratioPos}[$tpIdx];
				print qq
|$titleStrg{$targetPos}
<TABLE cellspacing=0>
<TR bgcolor=$darkColor>
<TH class="rbBorder">#</TH><TH class="rbBorder">&nbsp;Pos.&nbsp;<TH class="rbBorder">Peptides</TH><TH class="rbBorder">&nbsp;Charge&nbsp;</TH><TH class="rbBorder">&nbsp;$valueHeader&nbsp;</TH>
|;
				if ($ratioType eq 'Ratio') {$colspan++; print "<TH class=\"rbBorder\">&nbsp;$testStateName{$targetPos}$biasStrg1&nbsp;</TH><TH class=\"rbBorder\">&nbsp;$refStateName{$targetPos}$biasStrg1&nbsp;</TH>";}
				else {print "<TH class=\"rbBorder\">&nbsp;Occur.&nbsp;</TH>";}
				print "<TH class=\"$refClass\">&nbsp;Scores&nbsp;</TH>$sourceColStrg\n";
				print "</TR>\n";
				my $bgColor=$lightColor;
				my $numPep=0;
				my $prevPeptide='';
				my $prevBeg=0;
				my $maxBeg=0;
				foreach my $pepIdKey (sort{$peptideSets{$targetPos}{$a}{INFO}[0]<=>$peptideSets{$targetPos}{$b}{INFO}[0] || $peptideSets{$targetPos}{$a}{INFO}[1] cmp $peptideSets{$targetPos}{$b}{INFO}[1] || $peptideSets{$targetPos}{$a}{INFO}[2] cmp $peptideSets{$targetPos}{$b}{INFO}[2] || $peptideSets{$targetPos}{$a}{INFO}[3]<=>$peptideSets{$targetPos}{$b}{INFO}[3]} keys %{$peptideSets{$targetPos}}) { # beg,seq,vMod,charge
					$maxBeg=$sequenceBeg{$peptideSets{$targetPos}{$pepIdKey}{INFO}[1]}[-1];
					my $varMod=($peptideSets{$targetPos}{$pepIdKey}{INFO}[2])? " + $peptideSets{$targetPos}{$pepIdKey}{INFO}[2]" : '';
					my $numPepStrg='';
					if ($dispPepType eq 'distinct') {
						my $peptide="$peptideSets{$targetPos}{$pepIdKey}{INFO}[1]$varMod";
						if ($peptide ne $prevPeptide) {
							$numPep++;
							$numPepStrg=$numPep;
							$prevPeptide=$peptide;
						}
					}
					else {
						$numPep++;
						$numPepStrg=$numPep;
					}
					my $begStrg=($prevBeg != $peptideSets{$targetPos}{$pepIdKey}{INFO}[0])? "$peptideSets{$targetPos}{$pepIdKey}{INFO}[0]&nbsp;" : '';
					my $flagImgStrg=($peptideSets{$targetPos}{$pepIdKey}{DATA}[-1] eq 'NA')? '' : ($peptideSets{$targetPos}{$pepIdKey}{DATA}[-1]==1)? "<IMG src=\"$promsPath{images}/good.gif\">" : "<IMG src=\"$promsPath{images}/bad.gif\">";
					#my $ratio=($labelType eq 'FREE' && $ratioType eq 'SimpleRatio')? sprintf "%.0f",$peptideSets{$targetPos}{$pepIdKey}{DATA}[0] : ($peptideSets{$targetPos}{$pepIdKey}{DATA}[0]>=1)? sprintf "%.2f",$peptideSets{$targetPos}{$pepIdKey}{DATA}[0] : '1/'.sprintf "%.2f",1/$peptideSets{$targetPos}{$pepIdKey}{DATA}[0];
					my $ratio=($isPepIntensity)? sprintf "%.0f",$peptideSets{$targetPos}{$pepIdKey}{DATA}[0] : ($peptideSets{$targetPos}{$pepIdKey}{DATA}[0]>=1)? sprintf "%.2f",$peptideSets{$targetPos}{$pepIdKey}{DATA}[0] : '1/'.sprintf "%.2f",1/$peptideSets{$targetPos}{$pepIdKey}{DATA}[0];
					print qq
|<TR bgcolor="$bgColor">
<TD class="rBorder" align=right>&nbsp;$numPepStrg&nbsp;</TD><TD align=right>$begStrg</TD><TH class="font11" align=left nowrap>$flagImgStrg$peptideSets{$targetPos}{$pepIdKey}{INFO}[1]$varMod&nbsp;</TH>
<TH>$peptideSets{$targetPos}{$pepIdKey}{INFO}[3]<SUP>+</SUP></TH>
<TD align=center>&nbsp;$ratio&nbsp;</TD>
|;
					if ($ratioType eq 'Ratio') {
						my $testVal=sprintf "%.2e",$peptideSets{$targetPos}{$pepIdKey}{DATA}[1];
						my $refVal=sprintf "%.2e",$peptideSets{$targetPos}{$pepIdKey}{DATA}[2];
						print "<TD align=center>&nbsp;$testVal&nbsp;</TD><TD align=center>&nbsp;$refVal&nbsp;</TD>\n";
					}
					else {print "<TD align=center>&nbsp;x$peptideSets{$targetPos}{$pepIdKey}{DATA}[1]&nbsp;</TD>";} # (Super/Simple)Ratio: pep set occurence
					print "<TD align=center>&nbsp;$peptideSets{$targetPos}{$pepIdKey}{INFO}[4]&nbsp;</TD>\n";
					if ($numSources > 1) {
						my @srcFiles;
						foreach my $srcKey (@{$peptideSets{$targetPos}{$pepIdKey}{INFO}[5]{DATA}}) {push @srcFiles,$dataSources{$srcKey} || $srcKey;}
						if ($ratioType eq 'Ratio') {print "<TD>&nbsp;",join('<BR>',sort @srcFiles),"&nbsp;</TD>";} # $dataSources{$peptideSets{$targetPos}{$pepIdKey}{INFO}[5]}
						else {
							my $srcName=($algoVersion==3)? $peptideSets{$targetPos}{$pepIdKey}{INFO}[5]{NAME} : '...';
							print "<TD>&nbsp;<A href=\"javascript:void(null)\" onmouseover=\"popup('<U>Source(s):</U><BR>",join('<BR>',sort @srcFiles),"')\" onmouseout=\"popout()\">$srcName</A>&nbsp;</TD>";
						} # $peptideSets{$targetPos}{$pepIdKey}{INFO}[5]
					}
					print "</TR>\n";
					$prevBeg=$peptideSets{$targetPos}{$pepIdKey}{INFO}[0];
					$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
				}
				print "<TR bgcolor=\"$lightColor\"><TD colspan=\"$colspan\">&nbsp;<B>No peptide data available</B></TD></TR>\n" unless $numPep;
				print "<TR bgcolor=\"$darkColor\"><TD colspan=\"$colspan\"></TD></TR>\n" if $bgColor eq $darkColor;
				$lastPepStart=$maxBeg if $maxBeg > $lastPepStart;
				print "</TABLE>\n";
				print "<BR>\n" if $tpIdx < $#{$targetList{$ratioPos}};
			}
			print qq
|$biasStrg1$biasStrg2
</DIV><BR>
|;
		}
	}

	exit unless $numNormalRatios; # !!! CONDITIONAL EXIT !!!

	my $plotWidth=($protLength)? $protLength : $lastPepStart+50;
	$plotWidth=400 if $plotWidth < 400;
	$plotWidth=1000 if $plotWidth > 1000;
	#my $plotReferenceFlag=($labelType eq 'FREE' && $ratioType ne 'Ratio')? 'false' : 'true';

	###<javascipt for peptidePlot>###
	print "#==========#\n";
	foreach my $ratioPos (sort{if ($selRatioPos > 0) {abs($b)<=>abs($a)} else {abs($a)<=>abs($b)}} keys %protRatio) { # test then ref
		next unless $protRatio{$ratioPos}; # Undefined if all peptides were excluded by global CVout
		if ($normalRatio{$ratioPos}) {
			foreach my $targetPos (@{$targetList{$ratioPos}}) {
				my $trueTargetPos=abs($targetPos);
				my ($extraIdStrg,$protValueType,$protValue)=($isPepIntensity)? (':'.$targetPos,'protMean',$protMean{$targetPos}) : ('','protRatio',$protRatio{$ratioPos}); # $labelType eq 'FREE' && $ratioType ne 'Ratio'
				if (($isPepIntensity || $labelType eq 'FREE') && !defined($protValue)) { # no peptide at all for this state
					print qq
|document.getElementById('pepPlot:$selQuantifID:$ratioPos$extraIdStrg').innerHTML='&nbsp;&nbsp;&nbsp;<FONT class="title3">No peptide data available</FONT>'
|;
					next;
				}
				if ($isPepIntensity) { # $labelType eq 'FREE' && $ratioType ne 'Ratio'
					my $ratioPropStrg=($algoVersion==3)? ",'Source'" : '';
					my $pepColorStrg=($algoVersion==3)? "peptideColor: {map:'Source',type:'discrete'}" : '';
					print qq
|var PP_$trueTargetPos=new peptidePlot({div:'pepPlot:$selQuantifID:$ratioPos$extraIdStrg',width:$plotWidth,height:300,valueAxisLabel:'Log2(mean XIC)',valueLabel:'Mean XIC',
									protein: {name:'$protAlias$modStrg',length:$protLength,value:$protValue,relativeThresholds:[0.5,2]},
									peptideProperties: ['Score'$ratioPropStrg],
									convertValue: function(v) {return Math.log(v)/$log2},
									convertValueDisplayed: function(v) {return (v*1).toPrecision(3)},
									//convertValueDisplayed: function(v) {return Math.round(v*1)},
									$pepColorStrg
									});
|;
				}
				else {
					my $ratioPropStrg=($ratioType eq 'Ratio')? ",'Test','Source'" : ($algoVersion==3)? ",'Source'" : '';
					my $pepColorStrg=($ratioType eq 'Ratio' || $algoVersion==3)? "peptideColor: {map:'Source',type:'discrete'}" : '';
					print qq
|var PP_$trueTargetPos=new peptidePlot({div:'pepPlot:$selQuantifID:$ratioPos$extraIdStrg',width:$plotWidth,height:300,valueAxisLabel:'Log2(ratio)',valueLabel:'Ratio',
									protein: {name:'$protAlias$modStrg',length:$protLength,value:$protValue,relativeThresholds:[0.5,2]},
									reference: {label:'Reference: No fold change',value:1},
									peptideProperties: ['Score'$ratioPropStrg],
									convertValue: function(v) {return Math.log(v)/$log2},
									convertValueDisplayed: function(v) {return (v<1)? '1/'+(1/v).toPrecision(3) : (v*1).toPrecision(3)},
									$pepColorStrg
									});
|;
				}
				foreach my $pepIdKey (sort{$peptideSets{$targetPos}{$a}{INFO}[0]<=>$peptideSets{$targetPos}{$b}{INFO}[0] || $peptideSets{$targetPos}{$a}{INFO}[3]<=>$peptideSets{$targetPos}{$b}{INFO}[3]} keys %{$peptideSets{$targetPos}}) { # beg,charge
					my $vModStrg=($peptideSets{$targetPos}{$pepIdKey}{INFO}[2])? "'$peptideSets{$targetPos}{$pepIdKey}{INFO}[2]'" : 'null';
					my $excluded=($peptideSets{$targetPos}{$pepIdKey}{DATA}[-1] eq '0')? 1 : 0; # 'eq' because 'NA' possible
					print "PP_$trueTargetPos.addPeptide([[",join(',',@{$sequenceBeg{$peptideSets{$targetPos}{$pepIdKey}{INFO}[1]}}),"],'$peptideSets{$targetPos}{$pepIdKey}{INFO}[1]',$vModStrg,$peptideSets{$targetPos}{$pepIdKey}{INFO}[3],$peptideSets{$targetPos}{$pepIdKey}{DATA}[0],$excluded,'$peptideSets{$targetPos}{$pepIdKey}{INFO}[4]'"; # [beg],seq,mod,charge,ratio,exclu,score
					# source
					if ($algoVersion==3) {
						print ",'$peptideSets{$targetPos}{$pepIdKey}{INFO}[5]{NAME}'";
					}
					elsif ($ratioType eq 'Ratio') {
						my $testVal=sprintf "%.2e",$peptideSets{$targetPos}{$pepIdKey}{DATA}[1];
						my $refVal=sprintf "%.2e",$peptideSets{$targetPos}{$pepIdKey}{DATA}[2];
						print ",'$testVal/$refVal','",$dataSources{$peptideSets{$targetPos}{$pepIdKey}{INFO}[5]{DATA}->[0]},"'";
					}
					print "]);\n";
				}
				print "PP_$trueTargetPos.draw();\n";
			}
		}
	}
	exit;
}


####################<<<ajaxPeptideMaxQuant>>>#####################
sub ajaxPeptideMaxQuant {
	my $modProtID=param('id_prot');
	my ($protID,$modStrg)=($modProtID=~/^(\d+)-*(.*)/);
	my $dispMeasure=param('dispMeasure') || '';
	my $selTargetPos=param('ratio');
	my $trueTargetPos=abs($selTargetPos);

	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;
	my $projectID=&promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION');
	my ($quantifAnnot,$designID,$quantifType,$quantifModID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_DESIGN,M.CODE,Q.ID_MODIFICATION FROM QUANTIFICATION Q,QUANTIFICATION_METHOD M WHERE Q.ID_QUANTIFICATION_METHOD=M.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$selQuantifID");

	my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
	my ($labelType)=($labelStrg)? $labelStrg=~/LABEL=(.+)/ : ('FREE');
	$labelType=uc($labelType); # iTRAQ to ITRAQ
	my (%labelingInfo,%condToState,%condToObs);
	$labelingInfo{'LABELTYPE'}=$labelType;
	#if ($designID || $labelType =~ /SILAC|ITRAQ/) {
		foreach my $infoStrg (@labelInfo) {
			my ($setting,$valueStrg)=split('=',$infoStrg);
			$valueStrg=~s/#//g; # remove all ID tags
			@{$labelingInfo{$setting}}=split(';',$valueStrg);
		}
	#}
	my $ratioType=$labelingInfo{RATIO_TYPE}[0];

	if ($modStrg && !$quantifModID) { # fall back to whole protID if called with a modProtID for a non-modif quantif (eg. called from log2 plot)
		$modProtID=$protID;
		$modStrg='';
	}
	my ($protAlias,$protLength)=$dbh->selectrow_array("SELECT ALIAS,PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=$protID");
	$protLength=0 unless $protLength;

	my $paramCode=$dispMeasure || 'RATIO';
	my $sthPR;
	if ($modStrg) { # modification quantif
		my $dbModStrg=&dbEncodeModifPosition($modStrg);
		$sthPR=$dbh->prepare("SELECT QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION P
											LEFT JOIN PROTQUANTIF_MODRES PQMR ON P.ID_PROT_QUANTIF=PQMR.ID_PROT_QUANTIF
											LEFT JOIN MODIFIED_RESIDUE MR ON PQMR.ID_MODIF_RES=MR.ID_MODIF_RES
											INNER JOIN QUANTIFICATION_PARAMETER Q ON P.ID_QUANTIF_PARAMETER=Q.ID_QUANTIF_PARAMETER AND Q.CODE='$paramCode'
											WHERE P.ID_QUANTIFICATION=$selQuantifID AND ID_PROTEIN=$protID AND TARGET_POS=$trueTargetPos GROUP BY P.ID_PROT_QUANTIF
											HAVING GROUP_CONCAT(RESIDUE,POSITION ORDER BY POSITION SEPARATOR '.')='$dbModStrg'");
	}
	else { # whole protein quantif
		$sthPR=$dbh->prepare("SELECT QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION P,QUANTIFICATION_PARAMETER Q WHERE P.ID_QUANTIF_PARAMETER=Q.ID_QUANTIF_PARAMETER AND Q.CODE='$paramCode' AND ID_QUANTIFICATION=$selQuantifID AND ID_PROTEIN=$protID AND TARGET_POS=$trueTargetPos");
	}
	$sthPR->execute;
	my ($protValue)=$sthPR->fetchrow_array;
	$sthPR->finish;

	#########################
	###<Processing design>###
	#########################
	my ($quantifTargetName,$testStateID,$refStateID,$selStateID);
	my $sthExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
	if ($ratioType eq 'None') {
		(my $numBioRep,my $quantiObsIDs,$selStateID)=split(',',$labelingInfo{STATES}[$selTargetPos-1]);
		($quantifTargetName)=$sthExpCondName->execute($selStateID);
	}
	else {
		($testStateID,$refStateID)=($selTargetPos > 0)? split(/\//,$labelingInfo{RATIOS}[$trueTargetPos-1]) : reverse (split(/\//,$labelingInfo{RATIOS}[$trueTargetPos-1]));
		$sthExpCondName->execute($testStateID);
		my ($testStateName)=$sthExpCondName->fetchrow_array;
		$sthExpCondName->execute($refStateID);
		my ($refStateName)=$sthExpCondName->fetchrow_array;
		$quantifTargetName="$testStateName/$refStateName";
	}
	$sthExpCondName->finish;

	my (%anaList); #,%pepQuantifList
	foreach my $stateInfo (@{$labelingInfo{STATES}}) {
		my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$stateInfo);
		if ($ratioType eq 'None') {next if $expCondID != $selStateID;}
		else {next if $expCondID != $testStateID && $expCondID != $refStateID;}
		foreach my $bioRep (split(/\./,$quantiObsIDs)) {
			foreach my $techRep (split(/&/,$bioRep)) { # no fraction for SWATH (so far...)
				my ($obsID,$pepQuantID,$anaID,$tgPos)=split(':',$techRep);
				$anaList{$anaID}=1;
				#$pepQuantifList{$pepQuantID}=1;
			}
		}
	}
	my $anaString=join(',',keys %anaList);
	#my $pepQuantifString=join(',',keys %pepQuantifList);

	#my %runInfo;
	#my $sthRun=$dbh->prepare("SELECT A.ID_ANALYSIS,A.NAME,Q.NAME FROM ANALYSIS A,ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION IN ($pepQuantifString)");
	#$sthRun->execute;
	#while (my ($anaID,$anaName,$quantName)=$sthRun->fetchrow_array) {
	#	@{$runInfo{$anaID}}=($anaName,$quantName);
	#}
	#$sthRun->finish;

	$dbh->disconnect;

	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);

	my $protValueStrg;
	if ($ratioType eq 'None') {
		$protValueStrg=($dispMeasure eq 'MQ_INT')? 'Intensity' : ($dispMeasure eq 'MQ_IBAQ')? 'iBAQ'  : ($dispMeasure eq 'MQ_LFQ')? 'LFQ' : 'MS/MS count';
		$protValueStrg.='='.sprintf "%.0f",$protValue;
	}
	else {
		$protValueStrg='Ratio=';
		if ($protValue < 1) {$protValueStrg.='1/'.sprintf "%.2f",1/$protValue;}
		else {$protValueStrg.=sprintf "%.2f",$protValue;}
	}
	$modStrg=($modStrg)? '-'.$modStrg : '';

	print qq
|<TABLE><TR><TD nowrap><FONT class="title2">Peptide raw data for <A href="javascript:sequenceView($protID,'$anaString')">$protAlias$modStrg</A> in $quantifTargetName ($protValueStrg):</FONT>
&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'"/>
</TD></TR></TABLE>
<BR>
<BR><FONT class="title3">Not yet implemented...</FONT><BR><BR>
|;
	exit;
}

####################<<<ajaxPeptideSwath>>>#####################
sub ajaxPeptideSwath {
#print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1); # DEBUG
	my $modProtID=param('id_prot');
	my ($protID,$modStrg)=($modProtID=~/^(\d+)-*(.*)/);
	my $selRatioPos=param('ratio');
	my $trueRatioPos=abs($selRatioPos);
	#my $selRatioIdx=$selRatioPos-1;

	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;

	my $projectID=&promsMod::getProjectID($dbh,$selQuantifID,'QUANTIFICATION');
	my ($quantifAnnot,$quantifModID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_MODIFICATION FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
	my ($labelType)=($labelStrg)? $labelStrg=~/LABEL=(.+)/ : ('FREE');
	$labelType=uc($labelType); # iTRAQ to ITRAQ
	my (%labelingInfo,%condToState,%condToObs);
	$labelingInfo{'LABELTYPE'}=$labelType;
	foreach my $infoStrg (@labelInfo) {
		my ($setting,$valueStrg)=split('=',$infoStrg);
		$valueStrg=~s/#//g; # remove all ID tags
		@{$labelingInfo{$setting}}=split(';',$valueStrg);
	}
	#my $ratioType=$labelingInfo{RATIO_TYPE}[0]; # PROT_RATIO_PEP
	my $numRatios=scalar @{$labelingInfo{RATIOS}};
	my $numStates=scalar @{$labelingInfo{STATES}};

	if ($modStrg && !$quantifModID) { # fall back to whole protID if called with a modProtID for a non-modif quantif (eg. called from log2 plot)
		$modProtID=$protID;
		$modStrg='';
	}
	my ($protAlias,$protLength)=$dbh->selectrow_array("SELECT ALIAS,PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=$protID");
	$protLength=0 unless $protLength;

	my $trueProtRatio;
	if ($modStrg) { # modification quantif
		my $dbModStrg=&dbEncodeModifPosition($modStrg);
		($trueProtRatio)=$dbh->selectrow_array("SELECT QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION P
											LEFT JOIN PROTQUANTIF_MODRES PQMR ON P.ID_PROT_QUANTIF=PQMR.ID_PROT_QUANTIF
											LEFT JOIN MODIFIED_RESIDUE MR ON PQMR.ID_MODIF_RES=MR.ID_MODIF_RES
											INNER JOIN QUANTIFICATION_PARAMETER Q ON P.ID_QUANTIF_PARAMETER=Q.ID_QUANTIF_PARAMETER AND Q.CODE='RATIO'
											WHERE P.ID_QUANTIFICATION=$selQuantifID AND ID_PROTEIN=$protID AND TARGET_POS=$trueRatioPos GROUP BY P.ID_PROT_QUANTIF
											HAVING GROUP_CONCAT(RESIDUE,POSITION ORDER BY POSITION SEPARATOR '.')='$dbModStrg'");
	}
	else { # whole protein quantif
		($trueProtRatio)=$dbh->selectrow_array("SELECT QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION P,QUANTIFICATION_PARAMETER Q WHERE P.ID_QUANTIF_PARAMETER=Q.ID_QUANTIF_PARAMETER AND Q.CODE='RATIO' AND ID_QUANTIFICATION=$selQuantifID AND ID_PROTEIN=$protID AND TARGET_POS=$trueRatioPos");
	}
	my $protRatio=($selRatioPos < 0)? 1/$trueProtRatio : $trueProtRatio;


	#########################
	###<Processing design>###
	#########################
	my (%anaList,%pepQuantifList); #,%usedReplicates
	my ($testStateID,$refStateID)=($selRatioPos > 0)? split(/\//,$labelingInfo{RATIOS}[$trueRatioPos-1]) : reverse (split(/\//,$labelingInfo{RATIOS}[$trueRatioPos-1]));
	my $sthExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
	$sthExpCondName->execute($testStateID);
	my ($testStateName)=$sthExpCondName->fetchrow_array;
	$sthExpCondName->execute($refStateID);
	my ($refStateName)=$sthExpCondName->fetchrow_array;
	$sthExpCondName->finish;

	my ($testStatePos,$refStatePos);
	my $statePos=0;
	my $runCount=0;
	foreach my $stateInfo (@{$labelingInfo{STATES}}) {
		$statePos++;
		my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$stateInfo);
		if ($expCondID == $testStateID) {$testStatePos=$statePos;}
		elsif ($expCondID == $refStateID) {$refStatePos=$statePos;}
		foreach my $bioRep (split(/\./,$quantiObsIDs)) {
			foreach my $techRep (split(/&/,$bioRep)) { # no fraction for SWATH (so far...)
				$runCount++;
				next if ($expCondID != $testStateID && $expCondID != $refStateID);
				my ($obsID,$pepQuantID,$anaID,$tgPos)=split(':',$techRep);
				$anaList{$anaID}=$runCount;
				$pepQuantifList{$pepQuantID}=1;
			}
		}
	}
	my $anaString=join(',',keys %anaList);
	my $pepQuantifString=join(',',keys %pepQuantifList);

	my %runInfo;
	my $sthRun=$dbh->prepare("SELECT A.ID_ANALYSIS,A.NAME,Q.NAME FROM ANALYSIS A,ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION IN ($pepQuantifString)");
	$sthRun->execute;
	while (my ($anaID,$anaName,$quantName)=$sthRun->fetchrow_array) {
		@{$runInfo{$anaID}}=($anaName,$quantName);
	}
	$sthRun->finish;
#print "T=$testStatePos, R=$refStatePos<BR>\n";
	##>Processing ratio design
	my (%columnStrg,%design);
	#my @trueStatesPos=($selRatioPos > 0)? ($testStatePos,$refStatePos) : ($refStatePos,$testStatePos);
	foreach my $statePos ($testStatePos,$refStatePos) {
		my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$labelingInfo{STATES}[$statePos-1]);
		my $bioRepCount=0;
		foreach my $bioRep (split(/\./,$quantiObsIDs)) {
			$bioRepCount++;
			foreach my $techRep (split(/&/,$bioRep)) { # no fraction for SWATH (so far...)
				my ($obsID,$pepQuantID,$anaID,$tgPos)=split(':',$techRep);
				push @{$columnStrg{$statePos}},"['$runInfo{$anaID}[0]',$anaID,'$runInfo{$anaID}[0] > $runInfo{$anaID}[1]']";
				push @{$design{$statePos}},["State$statePos:BioRep$bioRepCount:$anaList{$anaID}",$anaID];
#print "State$statePos:BioRep$bioRepCount:$anaList{$anaID}<BR>\n";
			}
		}
	}


	####################
	###<Peptide data>###
	####################

	###<Loading peptide quantif data from table.txt & excluded.txt>###
	my $dataDir="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/data";
	my $resultDir="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/results";
	my %peptideData;
	foreach my $file ('table.txt','excluded.txt') {
		next unless -e "$dataDir/$file";
		my $excluded=($file eq 'excluded.txt')? 1 : 0;
		open(PEP,"$dataDir/$file") || die "ERROR: $!";
		my $protFound=0;
		while(<PEP>) {
			next if $.==1;
			if (!/^$modProtID\t/) {
				next unless $protFound;
				last; # file is ordered by prot
			}
			$protFound=1;
			chomp;
			my ($prID,$pepSeqVarMod,$charge,$fragIon,$fragIonCharge,$labelType,$state,$bioRep,$run,$intensity)=split(/\t/,$_);
			next if $state !~ /State($testStatePos|$refStatePos)/; # restrict to selected states
			next if $intensity eq 'NA';
			$bioRep=~s/$state.//; # '.' -> any single-character separator
			@{$peptideData{$pepSeqVarMod.'_'.$charge}{$fragIon.'_'.$fragIonCharge}{"$state:$bioRep:$run"}}=($intensity,$excluded);
#print "$_<BR>\n" if ($pepSeqVarMod =~ /^MIKPFFHSLSEK.+/ && $fragIon eq 'b5');
		}
		close PEP;
	}

	if (-e "$resultDir/CVoutPep.txt") {
		open(PEP,"$resultDir/CVoutPep.txt") || die "ERROR: $!";
		my $protFound=0;
		while(<PEP>) { # 75817	DSDQTTIAEDLMK_2_y5_1	State1	0.145431865367481
			next if $.==1;
			if (!/^$modProtID\t/) {
				next unless $protFound;
				last; # file is ordered by prot
			}
			$protFound=1;
			chomp;
			my ($prID,$pepCode,$state,$intensity)=split(/\t/,$_);
			my ($pepSeqVarMod,$charge,$fragIon,$fragIonCharge)=split('_',$pepCode);
			$pepSeqVarMod=~s/\*/:/g; # : converted to * because incompatible with MSstats
			#foreach my $bioRep (keys %{$peptideData{$pepSeqVarMod.'_'.$charge}{$state}}) { # bioReps & runs are not specified => assume all are excluded
			#	foreach my $run (keys %{$peptideData{$pepSeqVarMod.'_'.$charge}{$state}{$bioRep}}) {
			#		$peptideData{$pepSeqVarMod.'_'.$charge}{$fragIon.'_'.$fragIonCharge}{"$state:$bioRep:$run"}[1]=2 if $peptideData{$pepSeqVarMod.'_'.$charge}{$fragIon.'_'.$fragIonCharge}{"$state:$bioRep:$run"}; # set exclusion to level 2
			#	}
			#}
			foreach my $runKey (keys %{$peptideData{$pepSeqVarMod.'_'.$charge}{$fragIon.'_'.$fragIonCharge}}) {
				$peptideData{$pepSeqVarMod.'_'.$charge}{$fragIon.'_'.$fragIonCharge}{$runKey}[1]=2 if $runKey=~/^$state:/; # set exclusion to level 2
			}

		}
		close PEP;
	}

	###<Fetching peptides info from DB>###
	my (%peptideInfo,%peptideBeg,%varModName,%varModText);
	my $sthPI=$dbh->prepare("SELECT PEP_SEQ,GROUP_CONCAT(PM.ID_MODIFICATION,':',PM.POS_STRING ORDER BY PM.ID_MODIFICATION SEPARATOR '&'),CHARGE,ABS(PEP_BEG),P.ID_ANALYSIS,P.ID_PEPTIDE,QUERY_NUM,PEP_RANK,SCORE
							FROM PEPTIDE P
							LEFT JOIN PEPTIDE_MODIFICATION PM ON P.ID_PEPTIDE=PM.ID_PEPTIDE
							INNER JOIN PEPTIDE_PROTEIN_ATTRIB PPA ON P.ID_PEPTIDE=PPA.ID_PEPTIDE AND P.ID_ANALYSIS IN ($anaString)
							WHERE ID_PROTEIN=? AND PEP_RANK=1 GROUP BY P.ID_PEPTIDE");
	$sthPI->execute($protID);
	while (my ($pepSeq,$varModCode,$charge,$beg,$anaID,@info)=$sthPI->fetchrow_array) {
		my $vModStrg;
		if ($varModCode) {$vModStrg='&'.$varModCode;} else {$varModCode=$vModStrg='';}
		if ($peptideData{$pepSeq.$vModStrg.'_'.$charge}) { # match on seq+vMod+charge because ids are not recorded in table.txt)
			$peptideInfo{$pepSeq}{$varModCode}{$charge}{$anaID}=\@info; # pepID,qNum,rank,score
			$peptideBeg{$pepSeq}=$beg;
			$varModText{$varModCode}=&promsMod::decodeVarMod($dbh,$pepSeq,$varModCode,\%varModName) if $varModCode;
		}
	}
	$sthPI->finish;

	$dbh->disconnect;


	####<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);
	my $normalRatio=($protRatio==1000 || $protRatio==0.001)? 0 : 1; # inf or no inf
	my $protRatioStrg;
	if ($protRatio < 1) {$protRatioStrg=($normalRatio)? '1/'.sprintf "%.2f",1/$protRatio : '1/&infin;';} #'∞'
	else {$protRatioStrg=($normalRatio)? sprintf "%.2f",$protRatio : '&infin;';}
	$modStrg=($modStrg)? '-'.$modStrg : '';
	my $numTestRep=scalar @{$design{$testStatePos}};
	my $numRefRep=scalar @{$design{$refStatePos}};
	print qq
|<TABLE><TR><TD nowrap>
<FONT class="title2">Peptide fragments raw data for <A href="javascript:sequenceView($protID,'$anaString')">$protAlias$modStrg</A> in $testStateName/$refStateName (ratio=$protRatioStrg):</FONT>
&nbsp;&nbsp;<FONT class="title3">View:</FONT><SELECT class="title3" onchange="var divList=['pepHeat:$selQuantifID:$selRatioPos','pepList:$selQuantifID:$selRatioPos']; for (var i=0; i<2; i++) {var vis=(divList[i]==this.value)? '' : 'none'; document.getElementById(divList[i]).style.display=vis;}">
<OPTION value="pepHeat:$selQuantifID:$selRatioPos">Graphical</OPTION><OPTION value="pepList:$selQuantifID:$selRatioPos">List</OPTION>
</SELECT>
&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'"/>
</TD></TR></TABLE>
<DIV id="pepHeat:$selQuantifID:$selRatioPos"></DIV>
<DIV id="pepList:$selQuantifID:$selRatioPos" style="display:none">
<TABLE cellspacing=0>
<TR bgcolor="$darkColor">
<TH class="rbBorder">#</TH><TH class="rbBorder">&nbsp;Pos.&nbsp;<TH class="rbBorder">&nbsp;Peptides&nbsp;</TH><TH class="rbBorder">&nbsp;Charge&nbsp;</TH><TH class="rbBorder">&nbsp;Fragment&nbsp;</TH><TH class="rbBorder" colspan="$numTestRep">&nbsp;$testStateName&nbsp;</TH><TH class="bBorder" colspan="$numRefRep">&nbsp;$refStateName&nbsp;</TH></TR>
|;
	my $bgColor=$lightColor;
	my $pepCount=0;
	foreach my $pepSeq (sort{$peptideBeg{$a}<=>$peptideBeg{$b} || length($a)<=>length($b)} keys %peptideInfo) {
		foreach my $varMod (sort keys %{$peptideInfo{$pepSeq}}) {
			my $rootPepKey=($varMod)? $pepSeq.'&'.$varMod : $pepSeq;
			my $pepSeqVarMod=($varMod)? $pepSeq.'+'.$varModText{$varMod} : $pepSeq;
			foreach my $charge (sort{$a<=>$b} keys %{$peptideInfo{$pepSeq}{$varMod}}) {
				$pepCount++;
				my $pepKey=$rootPepKey.'_'.$charge;
				my $numFrag=scalar keys %{$peptideData{$pepKey}};
				print "<TR bgcolor=\"$bgColor\"><TD class=\"rBorder\" align=\"right\" valign=\"top\" rowspan=$numFrag>&nbsp;$pepCount&nbsp;</TD><TD align=\"right\" valign=\"top\" rowspan=$numFrag>&nbsp;$peptideBeg{$pepSeq}&nbsp;</TD><TH align=\"right\" valign=\"top\" nowrap rowspan=$numFrag>&nbsp;$pepSeqVarMod</TH><TH valign=\"top\"rowspan=$numFrag>&nbsp;$charge</SUP>+</SUP></TH>";
				my $firstFrag=1;
				foreach my $fragCharge (sort{&promsMod::sortSmart($a,$b)} keys %{$peptideData{$pepKey}}) {
					if ($firstFrag) {$firstFrag=0;}
					else {print "<TR bgcolor=\"$bgColor\">"}
					my ($frag,$ch)=split('_',$fragCharge);
					my $chStrg=($ch==1)? '+' : $ch.'+';
					print "<TH>&nbsp;$frag<SUP>$chStrg</SUP>&nbsp;</TH>";
					#>Looping through columns (runs)
					foreach my $statePos ($testStatePos,$refStatePos) {
						foreach my $refRun (@{$design{$statePos}}) {
							my ($runKey,$anaID)=@{$refRun};
							if ($peptideData{$pepKey}{$fragCharge}{$runKey}) {
								my $fragValue=($peptideData{$pepKey}{$fragCharge}{$runKey}[0]=~/None/)? 'None' : sprintf '%.0f',$peptideData{$pepKey}{$fragCharge}{$runKey}[0];
								print "<TD align=\"right\"><A href=\"javascript:drawSpectrum(null,'pep_",(join('_',@{$peptideInfo{$pepSeq}{$varMod}{$charge}{$anaID}}[0..2])),"')\">$fragValue&nbsp;</A>&nbsp;</TD>";
							}
							else {print "<TD align=\"right\">-&nbsp;</TD>";}
						}
					}
					print "</TR>\n";
				}
				$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
			}
		}
	}

	print qq
|</TABLE>
<BR>
</DIV>
<BR>
#==========#
var HMF=new heatMap({
	div:'pepHeat:$selQuantifID:$selRatioPos',
	//cellHeight:21,
	//editable:true,
	editableItems:{scope:true,type:true,color:true},
	entities:{row:'Fragments',column:'Run'},
	flagText:'excluded fragments',
	cellOnMouseOver: function(value,hCell) {
		var valueStrg;
		if (value==null) {valueStrg='* no value *';}
		else {valueStrg='Intensity: '+(Math.round(10*value)/10);}
		if (hCell) {
			var rowIdData=rowList[hCell.rowIdx].id.split(':'), colId=columnList[hCell.columnIdx].id;
			if (peptideId[rowIdData[0]][colId]) {valueStrg+='\\nPeptide score: '+(Math.round(100*peptideId[rowIdData[0]][colId][1])/100);}
			else {valueStrg+='\\n* no peptide identified *';}
		}
		return valueStrg;
	},
	singleCellSelection:true,
	cellOnClick:function(rowId,colId,status) {
		if (!status) return;
		var rowIdData=rowId.split(':');
		if (peptideId[rowIdData[0]][colId]) {drawSpectrum(null,peptideId[rowIdData[0]][colId][0]);}
		else {alert('Peptide was not identified in this run');}
	},
	normalization:{scope:'row',reference:'user',limitValues:{ref:0,min:0}},
	exportAsImage:['Export as image','Peptide_fragments_heat_map','$promsPath{cgi}/exportSVG.cgi']
});
|;

	##>Setting columns: runs and States
	my %stateGroups;
	my $anaCount=0; # index
	print 'HMF.setColumns([';
	foreach my $statePos ($testStatePos,$refStatePos) { # Test then Ref
		$stateGroups{$statePos}[0]=$anaCount; # index of 1st run in state
		print ',' if $statePos==$refStatePos;
		print join(',',@{$columnStrg{$statePos}});
		$anaCount+=scalar @{$columnStrg{$statePos}};
		$stateGroups{$statePos}[1]=$anaCount-1; # index of last run in state
	}
	print "])\n";

	##>Adding rows: peptide fragments
	my %peptideGroups;
	my $jsPeptideIdStrg="var peptideId={\n";
	$pepCount=0; # already used for List view
	my $fragCount=0;
	foreach my $pepSeq (sort{$peptideBeg{$a}<=>$peptideBeg{$b} || length($a)<=>length($b)} keys %peptideInfo) {
		foreach my $varMod (sort keys %{$peptideInfo{$pepSeq}}) {
			my $rootPepKey=($varMod)? $pepSeq.'&'.$varMod : $pepSeq;
			foreach my $charge (sort{$a<=>$b} keys %{$peptideInfo{$pepSeq}{$varMod}}) {
				$pepCount++;
				my $pepKey=$rootPepKey.'_'.$charge;
				@{$peptideGroups{$pepCount}{'PEP'}}=($pepSeq,$varMod,$charge);
				$peptideGroups{$pepCount}{'FRAG'}[0]=$fragCount; # index    unless $peptideGroups{$pepCount}{'FRAG'};
				my $firstFrag=1;
				foreach my $fragCharge (sort{&promsMod::sortSmart($a,$b)} keys %{$peptideData{$pepKey}}) {
					$peptideGroups{$pepCount}{'FRAG'}[1]=$fragCount;
					my ($frag,$ch)=split('_',$fragCharge);
					print "HMF.addRow(['$frag','$pepCount:$fragCount','Charge: $ch+'],[";
					if ($firstFrag) { # same peptideID for all fragments
						$jsPeptideIdStrg.=",\n" if $pepCount > 1;
						$jsPeptideIdStrg.="\t'$pepCount':{";
					}
					#>Looping through columns (runs)
					my @excluded;
					my $localFragCount=0; # index
					foreach my $statePos ($testStatePos,$refStatePos) {
						foreach my $refRun (@{$design{$statePos}}) {
							my ($runKey,$anaID)=@{$refRun};
							if ($firstFrag) { # generate colID 2 pep
								$jsPeptideIdStrg.=',' if $localFragCount;
								if ($peptideInfo{$pepSeq}{$varMod}{$charge}{$anaID}) {$jsPeptideIdStrg.="'$anaID':['pep_".(join('_',@{$peptideInfo{$pepSeq}{$varMod}{$charge}{$anaID}}[0..2]))."',$peptideInfo{$pepSeq}{$varMod}{$charge}{$anaID}[3]]";}
								else {$jsPeptideIdStrg.="'$anaID':null";}
							}
							print ',' if $localFragCount;

							if ($peptideData{$pepKey}{$fragCharge}{$runKey}) {
								if($peptideData{$pepKey}{$fragCharge}{$runKey}[0]=~/None/){print '';}
								else{print $peptideData{$pepKey}{$fragCharge}{$runKey}[0]}
								push @excluded,$localFragCount if $peptideData{$pepKey}{$fragCharge}{$runKey}[1];
							}
							$localFragCount++;
						}
					}
					print ']';
					if (scalar @excluded) {
						print ',[],[',join(',',@excluded),']';
					}
					print ");\n";
					$fragCount++;

					if ($firstFrag) {
						$jsPeptideIdStrg.='}';
						$firstFrag=0;
					}
				}
			}
		}
	}
	$jsPeptideIdStrg.="\n};\n";

	##>Defining groups
	print "HMF.defineGroups('column',[['$refStateName',1,'Reference state',$stateGroups{$refStatePos}[0],$stateGroups{$refStatePos}[1]],['$testStateName',2,'Test state',$stateGroups{$testStatePos}[0],$stateGroups{$testStatePos}[1]]]);\n";
	print "HMF.defineGroups('row',[";
	foreach my $pepCount (sort{$a<=>$b} keys %peptideGroups) {
		print ',' if $pepCount > 1;
		my ($pepSeq,$varMod,$charge)=@{$peptideGroups{$pepCount}{'PEP'}};
		my $pepSeqVarMod=($varMod)? $pepSeq.'+'.$varModText{$varMod} : $pepSeq;
		#my $rowID='pep_'.(join('_',@{$peptideInfo{$pepSeq}{$varMod}{$charge}}[0..2]));
		#print "['$pepSeqVarMod',$pepCount,'Charge: $charge+\\nScore: $peptideInfo{$pepSeq}{$varMod}{$charge}[3]\\nPosition: $peptideBeg{$pepSeq}-",($peptideBeg{$pepSeq}+length($pepSeq)-1),"',$peptideGroups{$pepCount}{FRAG}[0],$peptideGroups{$pepCount}{FRAG}[1]]";
		print "['$pepSeqVarMod',$pepCount,'Charge: $charge+\\nPosition: $peptideBeg{$pepSeq}-",($peptideBeg{$pepSeq}+length($pepSeq)-1),"',$peptideGroups{$pepCount}{FRAG}[0],$peptideGroups{$pepCount}{FRAG}[1]]";
	}
	print qq
|]);
HMF.draw();
var rowList=HMF.getRowList();
var columnList=HMF.getColumnList();
$jsPeptideIdStrg
|;

	exit;
}




sub ajaxProteinStatistics { # If multiple ratios were computed, ratio reversion (if any) applies only to selected ratio
	my $modProtID=param('id_prot');
	my $selRatioPos=param('ratio');
	my $trueRatioPos=abs($selRatioPos);
	my ($protID,$modStrg)=($modProtID=~/^(\d+)-*(.*)/);

	###<Connecting to the database>###
	my $dbh=&promsConfig::dbConnect;

	my ($protAlias)=$dbh->selectrow_array("SELECT ALIAS FROM PROTEIN WHERE ID_PROTEIN=$protID");
	my ($quantifAnnot,$designID)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT,ID_DESIGN FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$selQuantifID");
	my (%stateName,%condToState,@ratioInfo,$alterHyp,$adjStrg,$trueAlpha);
	my $ratioType='Ratio'; # default
	my ($quantifSoftware,$algoVersion)=('myProMS',0); # default
	foreach my $infoStrg (split('::',$quantifAnnot)) {
		my ($setting,$valueStrg)=split('=',$infoStrg);
		if ($setting eq 'SOFTWARE') {
			($quantifSoftware,$algoVersion)=split(';',$valueStrg);
		}
		elsif ($setting eq 'STATES') {
			my $statePos=0;
			if ($designID) {
				$valueStrg=~s/\#//g; # clean ID tags
				my $sthExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
				foreach my $stateData (split(';',$valueStrg)) {
					$statePos++;
					my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$stateData);
					$sthExpCondName->execute($expCondID);
					$condToState{$expCondID}=$statePos;
					($stateName{$statePos})=$sthExpCondName->fetchrow_array;
				}
			}
			else {
				foreach my $stateData (split(';',$valueStrg)) {
					$statePos++;
					(my $numBioRep,$stateName{$statePos},my $repPosStrg)=split(',',$stateData);
					$stateName{$statePos}=~s/\./\+/g;
				}
			}
		}
		elsif ($setting eq 'RATIOS') {
			$valueStrg=~s/#//g if $designID; # for label-free experiments, RATIOS are displayed with '#' to tell
			@ratioInfo=split(';',$valueStrg);
		}
		elsif ($setting eq 'FDR_CONTROL') {
			$adjStrg=($valueStrg eq 'TRUE')? '_ADJ' : '';
		}
		elsif ($setting eq 'FDR_ALPHA') { # not always defined
			$trueAlpha=$valueStrg/100;
		}
		elsif ($setting eq 'ALTER') {
			$alterHyp=$valueStrg;
		}
		elsif ($setting eq 'RATIO_TYPE') {
			$ratioType=$valueStrg;
			($quantifSoftware,$algoVersion)=('myProMS',1) if $ratioType eq 'Ratio';
		}
	}
	$algoVersion=2 if ($quantifSoftware eq 'myProMS' && !$algoVersion);
	$trueAlpha=0.05 unless $trueAlpha;
	$alterHyp='two.sided' unless $alterHyp;
	$adjStrg='_ADJ' if $ratioType=~/S\w+Ratio/; # overwrites FDR_CONTROL (just to be safe)
	my (%statData,%paramName);
	my $sthPS;
	if ($modStrg) { # modification quantif
		my $dbModStrg=&dbEncodeModifPosition($modStrg);
		$sthPS=$dbh->prepare("SELECT CODE,NAME,TARGET_POS,QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION P
						LEFT JOIN PROTQUANTIF_MODRES PQMR ON P.ID_PROT_QUANTIF=PQMR.ID_PROT_QUANTIF
						LEFT JOIN MODIFIED_RESIDUE MR ON PQMR.ID_MODIF_RES=MR.ID_MODIF_RES
						INNER JOIN QUANTIFICATION_PARAMETER Q ON P.ID_QUANTIF_PARAMETER=Q.ID_QUANTIF_PARAMETER
						WHERE P.ID_QUANTIFICATION=$selQuantifID AND ID_PROTEIN=$protID GROUP BY P.ID_PROT_QUANTIF
						HAVING GROUP_CONCAT(RESIDUE,POSITION ORDER BY POSITION SEPARATOR '.')='$dbModStrg'");
	}
	else {
		$sthPS=$dbh->prepare("SELECT CODE,NAME,TARGET_POS,QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION PQ,QUANTIFICATION_PARAMETER QP WHERE PQ.ID_QUANTIF_PARAMETER=QP.ID_QUANTIF_PARAMETER AND ID_QUANTIFICATION=$selQuantifID AND ID_PROTEIN=$protID"); # AND TARGET_POS=$ratioPos
	}
	$sthPS->execute;
	while (my ($parCode,$parName,$targetPos,$quantifValue)=$sthPS->fetchrow_array) {
		$targetPos=0 unless $targetPos; # multi-ratio data
		$statData{$parCode}{$targetPos}=$quantifValue;
		if ($parCode eq 'RATIO' && $targetPos==$trueRatioPos && $selRatioPos < 0) {
			$statData{$parCode}{$targetPos}=1/$statData{$parCode}{$targetPos};
		}
		if ($parCode eq 'SD_GEO') {$paramName{$parCode}='Coeff. variation';} # overwrites SD_GEO definition
		else {$paramName{$parCode}=$parName;}
	}
	$sthPS->finish;

	$dbh->disconnect;

	####<Starting HTML>###
	#my ($lightColor,$darkColor)=&promsConfig::getRowColors;
	print header(-type=>'text/plain',-charset=>'utf-8'); warningsToBrowser(1);
	my $numRatios=scalar @ratioInfo;
	my ($pepLabel,$distPepUsedStrg);
	if ($ratioType eq 'Ratio') {
		$pepLabel='Peptides used';
		$distPepUsedStrg=($numPepCode eq 'DIST_PEP_USED')? "$statData{DIST_PEP_USED}{0}/$statData{NUM_PEP_USED}{0}/$statData{NUM_PEP_TOTAL}{0}" : "$statData{NUM_PEP_USED}{0}/$statData{NUM_PEP_TOTAL}{0}";
	}
	else { # Super/Simple Ratio
		$pepLabel='All peptides';
		$distPepUsedStrg=$statData{NUM_PEP_TOTAL}{0};
	}
	$modStrg=($modStrg)? '-'.$modStrg : '';
	print qq
|<FONT class="title2">All statistics for $protAlias$modStrg:</FONT>&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type="button" class="font11" value=" Close " onclick="document.getElementById('displayDIV').style.display='none'"/>
<TABLE bgcolor=$darkColor>
<TR><TH nowrap align=right>$pepLabel :</TH><TH bgcolor=$lightColor align=left colspan=$numRatios>$distPepUsedStrg</TH></TR>
<TR><TH></TH>
|;
	my $ratioPos=0;
	my %linkedRatios;
	foreach my $ratioData (@ratioInfo) {
		$ratioPos++;
		my ($testStatePos,$refStatePos)=split(/\//,$ratioData);
		my $ratioTag='';
		if ($designID) {
			if ($testStatePos=~/%/) {
				$ratioTag=encode_utf8('°');
				my ($linkedTestPos)=$testStatePos=~/%(\d+)/;
				my ($linkedRefPos)=$refStatePos=~/%(\d+)/;
				@{$linkedRatios{$ratioPos}}=($linkedTestPos,$linkedRefPos); # needed for NUM_PEP_USED
				$testStatePos=~s/%\d+//;
				$refStatePos=~s/%\d+//;
			}
			$testStatePos=$condToState{$testStatePos};
			$refStatePos=$condToState{$refStatePos};
		}
		if ($ratioPos==$trueRatioPos) {
			if ($selRatioPos > 0) {print "<TH nowrap>&nbsp;<FONT color=#DD0000>$stateName{$testStatePos}$ratioTag/$stateName{$refStatePos}$ratioTag</FONT>&nbsp;</TH>";}
			else {print "<TH nowrap>&nbsp;<FONT color=#DD0000>$stateName{$refStatePos}$ratioTag/$stateName{$testStatePos}$ratioTag</FONT>&nbsp;</TH>";}
		}
		else {print "<TH nowrap>&nbsp;$stateName{$testStatePos}$ratioTag/$stateName{$refStatePos}$ratioTag&nbsp;</TH>";}
	}
	print "</TR>\n";
	my @quantifParams=('RATIO',"PVAL$adjStrg","NORM_PVAL$adjStrg",'SD_GEO');
	push @quantifParams,('CONF_LOW','CONF_UP') if $algoVersion==3;
	push @quantifParams,('DIST_PEP_USED','NUM_PEP_USED') if $ratioType=~/S\w+Ratio/;
	my %log2Ratio;
	foreach my $parCode (@quantifParams) {
		next if $parCode =~ /^CONF_(LOW|UP)/; # written with Coef Var
		next unless $paramName{$parCode}; # NORM_PVAL$adjStrg not always defined
		$paramName{$parCode} ='CV [Ratio 95% conf.]' if $parCode eq 'SD_GEO';
		print "<TR><TH nowrap align=right>&nbsp;$paramName{$parCode} :</TH>";
		foreach my $rPos (1..$ratioPos) {
			my $valueStrg='-';
			if (defined $statData{$parCode}{$rPos}) {
				if ($parCode eq 'RATIO') {
					my $absRatio=($statData{$parCode}{$rPos} >= 1)? $statData{$parCode}{$rPos} : 1/$statData{$parCode}{$rPos};
					$valueStrg=sprintf "%.2f",$absRatio;
					$valueStrg='1/'.$valueStrg if ($statData{$parCode}{$rPos} < 1 && $valueStrg > 1);
					$log2Ratio{$rPos}=log($statData{$parCode}{$rPos})/log(2);
					$valueStrg.=' (Log2: '.(sprintf "%.2f",$log2Ratio{$rPos}).')';
				}
				elsif ($parCode=~/PVAL/) {
					$valueStrg=($statData{$parCode}{$rPos})? sprintf "%.2e",$statData{$parCode}{$rPos} : '~0';
					if ($parCode eq "PVAL$adjStrg" && $statData{$parCode}{$rPos} <= $trueAlpha) {
						$valueStrg.="<IMG src=\"$promsPath{images}/good.gif\">";
					}
				}
				elsif ($parCode eq 'SD_GEO') {
					#$valueStrg=($ratioType eq 'Ratio')? sprintf "%.1f\%",100 * abs(log($statData{$parCode}{$rPos})/log($statData{RATIO}{$rPos})) : sprintf "%.1f",100*$statData{$parCode}{$rPos}/abs(log($statData{RATIO}{$rPos})/log(2));
					#$valueStrg.=($statData{'CONF_LOW'})? '% ['.(sprintf "%.2f",$statData{'CONF_LOW'}{$rPos}).' - '.(sprintf "%.2f",$statData{'CONF_UP'}{$rPos}).']' : '% (Conf. int.: &plusmn;'.(sprintf "%.2f",1.96*$statData{$parCode}{$rPos}).' Log2 unit)'; # +/- conf interval
					my $cvPc=($ratioType eq 'Ratio')? 100*abs(log($statData{$parCode}{$rPos})/log($statData{RATIO}{$rPos})) : 100*$statData{$parCode}{$rPos}/abs($log2Ratio{$rPos}); # v1 : v2,v3
					my ($ratioLow,$ratioUp)=($ratioType eq 'Ratio')? ($statData{'CONF_LOW'}{$rPos},$statData{'CONF_UP'}{$rPos}) : ($statData{'CONF_LOW'})? (2**$statData{'CONF_LOW'}{$rPos},2**$statData{'CONF_UP'}{$rPos}) : (2**($log2Ratio{$rPos} - 1.96*$statData{$parCode}{$rPos}),2**($log2Ratio{$rPos} + 1.96*$statData{$parCode}{$rPos})); # v1 : v3 : v2
					my $ratioLowStrg=($ratioLow >= 1)? sprintf "%.2f",$ratioLow : sprintf "1/%.2f",1/$ratioLow;
					my $ratioUpStrg=($ratioUp >= 1)? sprintf "%.2f",$ratioUp : sprintf "1/%.2f",1/$ratioUp;
					$valueStrg=(sprintf "%.1f",$cvPc)."% [$ratioLowStrg - $ratioUpStrg]"; #.' '.$statData{$parCode}{$rPos};
				}
				else {$valueStrg=$statData{$parCode}{$rPos};}
			}
			elsif ($linkedRatios{$rPos} && $parCode=~/(NUM|DIST)_PEP_USED/) {
				$valueStrg=($statData{$parCode}{$linkedRatios{$rPos}[0]} <= $statData{$parCode}{$linkedRatios{$rPos}[1]})? $statData{$parCode}{$linkedRatios{$rPos}[0]} : $statData{$parCode}{$linkedRatios{$rPos}[1]};
			}
			print "<TH bgcolor=$lightColor align=left nowrap>&nbsp;$valueStrg&nbsp;</TH>";
		}
		print "</TR>\n";
	}
	#if ($ratioType eq 'Ratio') { # old algo
	#	print "<TR><TH nowrap align=right>&nbsp;Confidence interval :</TH>";
	#	foreach my $rPos (1..$ratioPos) {
	#		my ($lowConf,$upConf);
	#		if ($alterHyp ne 'less') {
	#			$lowConf=(!defined $statData{'CONF_LOW'}{$rPos})? '?' : ($statData{'CONF_LOW'}{$rPos} >= 1)? sprintf "%.2f",$statData{'CONF_LOW'}{$rPos} : sprintf "1/%.2f",1/$statData{'CONF_LOW'}{$rPos};
	#			$lowConf=~s/0+\Z//;
	#		}
	#		if ($alterHyp ne 'greater') {
	#			$upConf=(!defined $statData{'CONF_UP'}{$rPos})? '?' : ($statData{'CONF_UP'}{$rPos} >= 1)? sprintf "%.2f",$statData{'CONF_UP'}{$rPos} : sprintf "1/%.2f",1/$statData{'CONF_UP'}{$rPos};
	#			$upConf=~s/0+\Z//;
	#		}
	#		my $confIntStrg=($alterHyp eq 'less')? "0-$upConf" : ($alterHyp eq 'greater')? "&ge;$lowConf" : "$lowConf-$upConf";
	#		print "<TH bgcolor=$lightColor align=left nowrap>&nbsp;$confIntStrg&nbsp;</TH>";
	#	}
	#	print "</TR>\n";
	#}
	print "</TABLE>\n";

	exit;
}

sub ajaxGetGoTerms {
	my ($goID,$aspect)=split(',',param('goStrg'));
	my $projectID=param('projectID');
	my $dbh=&promsConfig::dbConnect;
	my $selGOfromParentID=$dbh->prepare("SELECT ID_GOANALYSIS, PARAM_STRG from GO_ANALYSIS where ID_PARENT_GOANA=?");
	$selGOfromParentID->execute($goID);
	my %chilGoID;
	while (my($goAnaID, $paramStrg)=$selGOfromParentID->fetchrow_array) {
		my ($binStrg,$sizeStrg)=split(/;/,$paramStrg);
		my $binValue = (split(/=/,$binStrg))[1];
		$chilGoID{$goAnaID}=$binValue;
	}
	$selGOfromParentID->finish;
	my $onChange=param('onChange') || 'ajaxGoDecorateGraph';

	####<Starting HTML>###
	print header(-type=>'text/plain'); warningsToBrowser(1);
	print "<SELECT id=\"goTermsSELECT\" class=\"highlight\" onchange=\"$onChange(this.value,this.options[this.selectedIndex].text)\"><OPTION value=\"\">-= Select a GO term =-</OPTION>\n";

	if (!scalar(keys %chilGoID)) {
		my $resultFile="$promsPath{go_unix}/project_$projectID/$goID/results_$aspect.txt";
		open (TERMS,$resultFile);
		my $numTerms=0;
		while (<TERMS>) {
			if (/^(GO:\d+)\t([^\t]+)\t([^\t]+)/) {
				next if $3 > 1;
				my $pValue=sprintf "%.2e",$3;
				print "<OPTION value=\"$goID,$aspect,$1\">$2 [$1] ($pValue)</OPTION>\n";
				$numTerms++;
			}
		}
		close TERMS;
		unless ($numTerms) {print "<OPTION value='' disabled>No terms found</OPTION>\n";}
	}
	else {##there is a parent
		my @logRatioStrg=split(/;/,$dbh->selectrow_array("SELECT PARAM_STRG from GO_ANALYSIS where ID_GOANALYSIS=$goID"));
		my %ratioValue;
		foreach my $item (@logRatioStrg) {
			if ($item =~/logRatios/) {
				my ($value)=(split(/=/,$item))[1];
				my $count=0;
				my ($oldRatio, $strgRatioEnd);
				foreach my $ratio (split(/,/,$value)) {
					my $strgRatio="";
					if ($count == 0) {
						$strgRatio = "1/&infin;&nbsp;".&revertLog2($ratio);
					}
					elsif ($count==3) {
						$strgRatio = "$oldRatio&nbsp;".&revertLog2($ratio);
						$ratioValue{$count}=$strgRatio;
						$strgRatioEnd = &revertLog2($ratio)."&nbsp;&infin;";
						$ratioValue{$count+1}=$strgRatioEnd;
						last;
					}
					else {
						$strgRatio = "$oldRatio&nbsp;".&revertLog2($ratio);
					}
					$ratioValue{$count}=$strgRatio;
					$oldRatio=&revertLog2($ratio);
					$count++;
				}
			}

		}

		foreach my $goAnaID (sort{$chilGoID{$a} <=>$chilGoID{$b}} keys %chilGoID) {
			print "<OPTGROUP label=\"bin $ratioValue{$chilGoID{$goAnaID}} \">";
			my $resultFile="$promsPath{go_unix}/project_$projectID/$goAnaID/results_$aspect.txt";
			open (TERMS,$resultFile);
			my $numTerms=0;
			while (<TERMS>) {
				if (/^(GO:\d+)\t([^\t]+)\t([^\t]+)/) {
					next if $3 > 1;
					my $pValue=sprintf "%.2e",$3;
					print "<OPTION value=\"$goAnaID,$aspect,$1,$chilGoID{$goAnaID}:$ratioValue{$chilGoID{$goAnaID}}\">$2 [$1] ($pValue)</OPTION>\n";
					$numTerms++;
				}
			}
			close TERMS;
			print "</OPTGROUP>";
			unless ($numTerms) {print "<OPTION value='' disabled>No terms found</OPTION>\n";}
		}
	}
	print "</SELECT>\n";
	$dbh->disconnect;
}

sub ajaxGetGoTermProteins {
	my ($goID,$aspect,$termID)=split(',',param('goStrg'));
	my $projectID=param('projectID');

	####<Starting HTML>###
	print header(-type=>'text/plain',-charset => 'utf-8'); warningsToBrowser(1);
	my $resultFile="$promsPath{go_unix}/project_$projectID/$goID/results_$aspect.txt";
	my $protIdStrg;
	open (TERMS,$resultFile);
	while (<TERMS>) {
		if (/^$termID\t/) {
			$protIdStrg=(split(/\t/,$_))[5];
			last;
		}
	}
	close TERMS;
	print $protIdStrg;
}

sub ajaxListDecorateGraph {
    my $listID=param('listID');
	my $quantifModID=param('modifID');

    ###<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);

	###<Fetching list of proteins/sites>###
	my $dbh=&promsConfig::dbConnect;

    my $sthList=$dbh->prepare("SELECT NAME,LIST_TYPE FROM CATEGORY WHERE ID_CATEGORY=?");
	$sthList->execute($listID);
    my ($listName,$listType)=$sthList->fetchrow_array;
	$listName=&promsMod::shortenName($listName,30);
	$sthList->finish;
	my $subMatch=0;
	my $sthProt;
	if ($listType eq 'SITE' && $quantifModID > 0) { # compare 2 set of sites
		$sthProt=$dbh->prepare("SELECT ID_PROTEIN,ID_MODIFICATION,GROUP_CONCAT(RESIDUE,POSITION ORDER BY ID_MODIFICATION,POSITION SEPARATOR '.') FROM CATEGORY_PROTEIN CP,MODIFICATION_SITE MS
								WHERE CP.ID_CATEGORY_PROTEIN=MS.ID_CATEGORY_PROTEIN AND CP.ID_CATEGORY=? GROUP BY CP.ID_CATEGORY_PROTEIN");
	}
	else {
		$sthProt=$dbh->prepare("SELECT DISTINCT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=?");
		$subMatch=1 if $quantifModID; # match proteins on sites
	}
    print $listName.'::'.$subMatch.'::';
    $sthProt->execute($listID);
    my @protList;
    while (my ($protID,$modID,$modifCode)=$sthProt->fetchrow_array) {
		if ($modID) {
			push @protList,$protID.'-'.$modifCode if $modID==$quantifModID; # skip site if modif are different! => will not match anyway
		}
		else {push @protList,$protID;}
    }
    $sthProt->finish;
    $dbh->disconnect;
	if (scalar @protList) {print join(';',@protList);}
    else {print 0;}
    exit;
}


sub ajaxManageSaveProteins {
	my $job=param('job');
	my $listType=param('listType') || 'PROT';
	my $themeID=param('themeID') || 0;

	###<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);

	my $dbh=&promsConfig::dbConnect;

	###<Fetching list of available themes/classifications>###
	if ($job eq 'getThemes') {
		my $projectID=param('projID');
		my %themes;
		my $sthTh=$dbh->prepare("SELECT ID_CLASSIFICATION,NAME FROM CLASSIFICATION WHERE ID_PROJECT=$projectID");
		$sthTh->execute;
		while (my ($thID,$thName)=$sthTh->fetchrow_array) {
			$themes{$thID}=$thName;
		}
		$sthTh->finish;
		my $entities=($listType eq 'SITE')? 'sites' : 'proteins'; # used between SPAN because can be changed by JS ajaxManageSaveProteins() in promsMod
		print qq
|<TABLE border=0 bgcolor=$darkColor id="storeList"><TR>
<TH valign=top nowrap>&nbsp;<FONT style="font-size:15px;">Save checked <SPAN id="listTypeEntity1">$entities</SPAN> in</FONT></TH>
<TD valign=top>
<SELECT id="themeSELECT" style="width:250px;font-weight:bold;" onchange="ajaxManageSaveProteins(this.value,saveListGlobals.listType)">
|; # saveListGlobals.listType: global JS variable because of possible swtich between PROT and SITE without new AJAX call
		my $displayStatus;
		if (scalar keys %themes) { # exists at least 1 classification
			print "<OPTION selected value=\"getLists:0\">-=Select a Theme=-</OPTION>\n";
			foreach my $thID (sort{lc($themes{$a}) cmp lc($themes{$b})} keys %themes) {
				print "<OPTION value=\"getLists:$thID\">$themes{$thID}</OPTION>\n";
			}
			$displayStatus='none';
		}
		else {$displayStatus='block';}
		print qq
|<OPTION value="getLists:-1">-=Create a new Theme=-</OPTION>
</SELECT><BR>
<INPUT type="text" id="themeNameIN" value="" placeholder="New Theme name" style="width:250px;display:$displayStatus"/>
</TD>
<TD valign=top><DIV id="listSelectionDIV">|;
		my $disabSave='disabled';
		unless (scalar keys %themes) { # no Theme
			$disabSave='';
			print qq
|<SELECT id="listSelectionSELECT" style="width:250px;font-weight:bold" onchange="ajaxManageSaveProteins(this.value,'$listType')">
<OPTION value="selectList:-1:-1">-=Create a new List=-</OPTION>
</SELECT><BR>
<INPUT type="text" id="listNameIN" value="" placeholder="New List name" style="width:250px;"/>
|;
		}
		print qq
|</DIV></TD>
<TD nowrap valign=top><INPUT type="button" id="saveProtBUTTON" value=" Save " onclick="ajaxManageSaveProteins('saveProt')" $disabSave>&nbsp;<INPUT type="button" value="Cancel" onclick="ajaxManageSaveProteins('cancel')"></TD>
<TR><TD></TD>
<TH align=left colspan=3><INPUT type="radio" name="replace" id="replace0" value="0" checked>Add new <SPAN id="listTypeEntity2">$entities</SPAN> to List<BR>
<INPUT type="radio" name="replace" id="replace1" value="1">Replace List contents with new <SPAN id="listTypeEntity3">$entities</SPAN></TH>
</TR></TABLE>
|;
	}
	####<Fetching list of available lists/categories>####
	elsif ($job eq 'getLists') {
		my %lists;
		if ($themeID > 0) {
			my $sthL=$dbh->prepare("SELECT ID_CATEGORY,NAME,DISPLAY_POS FROM CATEGORY WHERE ID_CLASSIFICATION=$themeID AND LIST_TYPE='$listType'");
			my @sthNoMod=(
				$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,COMPARISON CO WHERE CO.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID AND LIST_TYPE='$listType'"), # COMPARISON
				$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,CAT_COMPARISON CC WHERE CC.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID AND LIST_TYPE='$listType'"), # CAT_COMPARISON
				$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,GO_ANALYSIS GA WHERE GA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID AND LIST_TYPE='$listType'"), # GO_ANALYSIS
				$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,PATHWAY_ANALYSIS PA WHERE PA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID AND LIST_TYPE='$listType'"), # PATHWAY_ANALYSIS
				$dbh->prepare("SELECT CA.ID_CATEGORY FROM CATEGORY CA,EXPLORANALYSIS EA WHERE EA.ID_CATEGORY=CA.ID_CATEGORY AND ID_CLASSIFICATION=$themeID AND LIST_TYPE='$listType'") # EXPLORANALYSIS
			);
			$sthL->execute;
			while (my ($listID,$listName,$displayPos)=$sthL->fetchrow_array) {
				@{$lists{$listID}}=($listName,$displayPos,0);
			}
			foreach my $sth (@sthNoMod) {
				$sth->execute;
				while (my ($listID)=$sth->fetchrow_array) {
					$lists{$listID}[2]=1 if $lists{$listID}; # update modif flag
				}
			}
			$sthL->finish;
			foreach my $sth (@sthNoMod) {$sth->finish;}
		}
		my $displayStatus;
		print "<SELECT id=\"listSelectionSELECT\" style=\"width:250px;font-weight:bold\" onchange=\"ajaxManageSaveProteins(this.value,'$listType')\">\n";
		if (scalar keys %lists) {
			print "<OPTION value=\"selectList:$themeID:0\">-=Select a List=-</OPTION>\n";
			foreach my $listID (sort{$lists{$a}[1] cmp $lists{$b}[1]} keys %lists) {
				print "<OPTION value=\"selectList:$themeID:$listID\"";
				print ' disabled' if $lists{$listID}[2]; # list is used and cannot be modified
				print ">$lists{$listID}[0]</OPTION>\n";
			}
			$displayStatus='none';
		}
		else {$displayStatus='block';}
		print qq
|<OPTION value="selectList:$themeID:-1">-=Create a new List=-</OPTION>
</SELECT>
<INPUT type="text" id="listNameIN" value="" placeholder="New List name" style="width:250px;display:$displayStatus"/>
|;
	}
	####<Saving proteins in List>####
	else {
		my $listID=param('listID');
		my $modifID=param('modifID');
		my $replace=(param('replace'))? param('replace') : 0;
		my ($themeName,$listName);
		if ($themeID==-1) { # Adding new theme
			my $projectID=param('projID');
			$themeName=param('themeName');
			($themeID)=$dbh->selectrow_array("SELECT MAX(ID_CLASSIFICATION) FROM CLASSIFICATION");
			$themeID++;
			my $sthInsCL=$dbh->prepare("INSERT INTO CLASSIFICATION (ID_CLASSIFICATION,NAME,ID_PROJECT) VALUES ($themeID,?,?)");
			$sthInsCL->execute($themeName,$projectID);
			$sthInsCL->finish;
			$listID=-1; # to be safe
		}
		else {
			$themeID=~s/\D+//g; # clean to prevent prepare
			($themeName)=$dbh->selectrow_array("SELECT NAME FROM CLASSIFICATION WHERE ID_CLASSIFICATION=$themeID");
		}
		my %oldList;
		my $checkOld=0;
		if ($listID==-1) { # Creating new list
			$replace=0; # no need to clean list
			$listName=param('listName');
			($listID)=$dbh->selectrow_array("SELECT MAX(ID_CATEGORY) FROM CATEGORY");
			$listID++;
			my ($displayPos)=$dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM CATEGORY WHERE ID_CLASSIFICATION=$themeID");
			$displayPos++;
			my $sthInsCA=$dbh->prepare("INSERT INTO CATEGORY (ID_CATEGORY,NAME,LIST_TYPE,DISPLAY_POS,ID_CLASSIFICATION,UPDATE_USER,UPDATE_DATE) VALUES ($listID,?,?,$displayPos,?,'$userID',NOW())");
			$sthInsCA->execute($listName,$listType,$themeID);
			$sthInsCA->finish;
		}
		else { # use existing list
			$listID=~s/\D+//g; # clean to prevent SQL injection
			($listName)=$dbh->selectrow_array("SELECT NAME FROM CATEGORY WHERE ID_CATEGORY=$listID");

			if ($replace) {
				###<Delete previous category contents
				$dbh->do("DELETE FROM MODIFICATION_SITE WHERE ID_CATEGORY=$listID");
				$dbh->do("DELETE FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$listID");
			}
			else { # Add to old contents
				###<Fetching existing protein/site contents to prevent duplicates
				if ($listType eq 'SITE') {
					#my $sthOld=$dbh->prepare("SELECT ID_PROTEIN,GROUP_CONCAT('[',ID_MODIFICATION,']',RESIDUE,POSITION ORDER BY ID_MODIFICATION,POSITION SEPARATOR '.') FROM CATEGORY_PROTEIN CP,MODIFICATION_SITE MS WHERE CP.ID_CATEGORY_PROTEIN=MS.ID_CATEGORY_PROTEIN AND CP.ID_CATEGORY=$listID GROUP BY CP.ID_CATEGORY_PROTEIN");
					#$sthOld->execute;
					#while (my ($protID,$modCode)=$sthOld->fetchrow_array) {
					#	if ($modCode=~/\]\d/) { # (occurence) => ambiguous
					#		$modCode=~s/(\[\d+\]\d+)\.(.+)/$2\.$1/; # Occurence is likely first when ORDER BY POSITION => send at the end. WARNING: Won't work for multi-modif site
					#	}
					#	$oldList{"$protID-$modCode"}=1;
					#}
					#$sthOld->finish;
					&promsQuantif::fetchSiteList($dbh,$listID,\%oldList); # no focusModID => format='[modID]modCode'
				}
				else {
					my $sthOld=$dbh->prepare("SELECT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$listID");
					$sthOld->execute;
					while (my ($protID)=$sthOld->fetchrow_array) {
						$oldList{$protID}=1;
					}
					$sthOld->finish;
				}
				$checkOld=1 if scalar keys %oldList;
			}
		}

		if ($listType eq 'PROT') {
			###<Removing duplicates (eg. if from modification quantif) and skipping already recorded prot
			my %newProtID;
			foreach my $modProtID (split(',',param('protID'))) {
				my ($protID,$modCode)=split('-',$modProtID);
				next if ($checkOld && $oldList{$protID}); # protein is already in category
				$newProtID{$protID}=1;
			}

			###<Storing protein list
			my $sthInsCP=$dbh->prepare("INSERT INTO CATEGORY_PROTEIN (ID_PROTEIN,ID_CATEGORY) VALUES (?,$listID)");
			foreach my $protID (keys %newProtID) {
				$sthInsCP->execute($protID);
				$oldList{$protID}=1; # In case multiple occurences of same protID in list (eg called from listProteins.cgi)
			}
			$sthInsCP->finish;
		}
		elsif ($listType eq 'SITE') {
			my $sthInsCP=$dbh->prepare("INSERT INTO CATEGORY_PROTEIN (ID_PROTEIN,ID_CATEGORY) VALUES (?,$listID)");
			my $sthInsMS=$dbh->prepare("INSERT INTO MODIFICATION_SITE (ID_CATEGORY,ID_MODIFICATION,ID_CATEGORY_PROTEIN,RESIDUE,POSITION) VALUES ($listID,$modifID,?,?,?)");
			foreach my $modProtID (split(',',param('protID'))) {
				my ($protID,$modCode)=split('-',$modProtID);
				next unless $modCode; # in case of mixte of sites and prot (eg from compareQuantifications.cgi)
				if ($checkOld) {
					next if $oldList{"$protID-[$modifID]$modCode"}; # add [modID] before modCode to match list format & check if already in list
				}
				#<Encode modCode for DB storage
				if ($modCode=~/~/) {
					$modCode=~s/^ProtNt~/=0./; # position was absent
					$modCode=~s/^PepNt(\d+)~/=$1./;
					$modCode=~s/^(\d+)~/-$1./; # normal res
					$modCode=~s/ProtCt:/\*99999./; # position was absent
					$modCode=~s/PepCt(\d+):/\*$1./;
					$modCode=~s/(\d+):/+$1./; # normal res
					$modCode=~s/\///;
				}
				else {
					$modCode=~s/ProtNt/n0/; # position was absent
					$modCode=~s/ProtCt/c99999/; # position was absent
					$modCode=~s/PepNt/n/;
					$modCode=~s/PepCt/c/;
				}
				$sthInsCP->execute($protID);
				my ($catProtID)=$dbh->last_insert_id(undef,undef,'CATEGORY_PROTEIN','ID_CATEGORY_PROTEIN');

				#<Separate individual site
				foreach my $site (split(/\./,$modCode)) {
					$site=~/(.)(\d+)/; # (site)(position)
					$sthInsMS->execute($catProtID,$1,$2);
				}

				$sthInsMS->finish;
				$sthInsCP->finish;
			}
		}

		###<Updating Theme
		$dbh->do("UPDATE CLASSIFICATION SET UPDATE_USER='$userID',UPDATE_DATE=NOW() WHERE ID_CLASSIFICATION=$themeID");
		$dbh->commit;

		print qq
|<B>Selected proteins were saved in List <FONT color=#DD0000>$listName</FONT> of Theme <FONT color=#DD0000>$themeName</FONT>.</B>
<!-- THEME_ID=$themeID THEME_NAME=$themeName### LIST_ID=$listID LIST_NAME=${listName}::: (parsed by JS when needed) -->
|;

	}

	$dbh->disconnect;
	exit;
}

sub ajaxRestrictProteinList {
	my $projectID=param('projectID');
	my $submitForm=param('submitForm') || '';
	my $noSelect=param('noSelect') || 0;

	my $dbh=&promsConfig::dbConnect;
	my $sthL=$dbh->prepare("SELECT T.ID_CLASSIFICATION,T.NAME,L.ID_CATEGORY,L.NAME,L.DISPLAY_POS,LIST_TYPE FROM CLASSIFICATION T,CATEGORY L WHERE T.ID_CLASSIFICATION=L.ID_CLASSIFICATION AND T.ID_PROJECT=$projectID");
	$sthL->execute;
	my (%savedLists,%themeInfo);
	while (my ($themeID,$themeName,$listID,$listName,$listPos,$type)=$sthL->fetchrow_array) {
		$type='PROT' unless $type;
		$themeInfo{$themeID}=$themeName;
		@{$savedLists{$themeID}{$listID}}=($listName,$listPos,$type);
	}
	$sthL->finish;
	$dbh->disconnect;

	###<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
	unless ($noSelect) { # otherwise <SELECT> and 1st "-= Selec =-" option => written by JS function
		print "<SELECT name=\"restrictList\" class=\"title3\"";
		print " onchange=\"document.$submitForm.submit()\"" if $submitForm;
		print ">\n<OPTION value=\"\">-= Select =-</OPTION>\n";
	}
	foreach my $themeID (sort{lc($themeInfo{$a}) cmp lc($themeInfo{$b})} keys %themeInfo) {
		print "<OPTGROUP label=\"$themeInfo{$themeID}\">\n";
		foreach my $listID (sort{lc($savedLists{$themeID}{$a}[1]) cmp lc($savedLists{$themeID}{$b}[1])} keys %{$savedLists{$themeID}}) {
			print "<OPTION value=\"$listID\"";
			print ' selected' if $listID==$restrictListID; # in case a list was already selected
			my $typeStrg=($savedLists{$themeID}{$listID}[2] eq 'SITE')? ' [Sites]' : '';
			print ">$savedLists{$themeID}{$listID}[0]$typeStrg</OPTION>\n";
		}
		print "</OPTGROUP>\n";
	}
	print "</SELECT>\n";
	exit;
}


sub ajaxDisplayCorrelationMatrix {
#print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1); # DEBUG
	my $projectID=param('projectID');

	my $dbh=&promsConfig::dbConnect;

	my $sthQuantif=$dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
	my $sthState=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
	my $sthAnaInfo=$dbh->prepare("SELECT A.NAME,S.NAME FROM ANALYSIS A,SAMPLE S WHERE A.ID_SAMPLE=S.ID_SAMPLE AND A.ID_ANALYSIS=?");

	$sthQuantif->execute($selQuantifID);
	my ($quantifAnnot)=$sthQuantif->fetchrow_array;
	$quantifAnnot=~s/::/\$/g; # easier to match
	my ($labelType)=($quantifAnnot=~/LABEL=([^\$]+)/);
	my ($quantifSoftware,$algoVersion)=($quantifAnnot=~/=([^=]+MSstats)/)? ($1,0) : ($quantifAnnot=~/=myProMS;(\d+)/)? ('myProMS',$1) : ('myProMS',2);
	my ($ratioType)=($quantifAnnot=~/RATIO_TYPE=([^\$]+)/); $ratioType='Ratio' unless $ratioType;
	my ($stateStrg)=($quantifAnnot=~/STATES=([^\$]+)/);
	$stateStrg=~s/#//g;

	my (%states,%labels,%labelData);
	my $minCorrelation=1;
	my $drawGroups=0;
	if ($quantifSoftware eq 'myProMS') {
		my $numStates=0;
		my (%labelsInfo,%sampChannelOcc);
		foreach my $stateData (split(';',$stateStrg)) {
			$numStates++;
			my ($numBioRep,$quantiObsIDs,$expCondID)=split(',',$stateData);
			$sthState->execute($expCondID);
			($states{$numStates})=$sthState->fetchrow_array;
			my (%parentSampleName,%anaName,%channelName);
			my $currBioRep=0;
			foreach my $bioReplicate (split(/\./,$quantiObsIDs)) { # bio rep separator
				$currBioRep++;
				my $multiTechRep=($bioReplicate=~/&/)? 1 : 0;
				my $numTechRep=0;
				foreach my $techReplicate (split(/&/,$bioReplicate)) { # tech rep separator
					$numTechRep++;
					#my $multiFrac=($techReplicate=~/\+/)? 1 : 0;
					my ($firstFrac)=($techReplicate=~/^([^\+]+)/);
					my @obsData=split(/:/,$firstFrac);
					my ($parQuantifID,$anaID,$targetPos);
					if (scalar @obsData==4) {($parQuantifID,$anaID,$targetPos)=@obsData[1..3];} # design post 07/2013 -> obsID,parQuantifID,anaID,targetPos
					else { # older design -> parQuantifID,anaID
						($parQuantifID,$anaID)=@obsData;
						$targetPos=0;
					}
					unless ($anaName{$anaID}) {
						$sthAnaInfo->execute($anaID);
						($anaName{$anaID},$parentSampleName{$anaID})=$sthAnaInfo->fetchrow_array;
					}
					my $targetName='';
					if ($targetPos) { # labeling
						unless ($channelName{$parQuantifID}) {
							$sthQuantif->execute($parQuantifID);
							my ($qAnnot)=$sthQuantif->fetchrow_array;
							my @quantifInfo=split(/::/,$qAnnot);
							foreach my $channelInfo (@quantifInfo[1..$#quantifInfo]) {
								my ($chanPos,$chanName)=split(/;/,$channelInfo); # + other unused elements
								$channelName{$parQuantifID}{$chanPos}=$chanName;
							}
						}
						$targetName=' > '.$channelName{$parQuantifID}{$targetPos};
					}
					@{$labelsInfo{"State$numStates-biorep$currBioRep-techRep$numTechRep"}}=($numStates,$parentSampleName{$anaID},$anaName{$anaID},$targetName);
#print "//State$numStates-biorep$currBioRep-techRep$numTechRep\n";
					$sampChannelOcc{$parentSampleName{$anaID}.$targetName}++;
				}
			}
		}

		$drawGroups=1 if $numStates < scalar keys %labelsInfo; # bio|tech replicates => draw label groups
		my( %labelsName,%groupName);
		foreach my $label (keys %labelsInfo) {
			$labelsName{$label}=($drawGroups)? $labelsInfo{$label}[1] : "[$states{$labelsInfo{$label}[0]}] $labelsInfo{$label}[1]";
			my $sampChan=$labelsInfo{$label}[1].$labelsInfo{$label}[3];
			$labelsName{$label}.=($sampChannelOcc{$sampChan}==1)? $labelsInfo{$label}[3] : ' > '.$labelsInfo{$label}[2].$labelsInfo{$label}[3];
#$labelsName{$label}.=" ($label)"; # DEBUG
		}

		my %correlStages=('BEFORE'=>'Beforematrixcorrelation.txt','AFTER'=>'Aftermatrixcorrelation.txt');
		foreach my $stage (keys %correlStages) {
			my $corrMatrixFile="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/results/$correlStages{$stage}";
			#my $corrMatrixFile="$promsPath{quantification}/project_$projectID/quanti_$selQuantifID/results/Afteroutmatrixcorrelation.txt";
			open(CORR,$corrMatrixFile) or die "ERROR: Cannot open file '$corrMatrixFile': $!";
			my (%convertLabels,@labelArray);
my @rawLabelList; # Temp
			while (<CORR>) {
				chomp;
				if ($.==1) {
					my %bioRepCountShift=(A=>0);
					my %numBioRepInState;
					my $prevExpCode='';
					foreach my $rawLabel (split(/\t/,$_)) {
						next unless $rawLabel; # top left cell is empty
push @rawLabelList,$rawLabel;
						my ($expCode,$stateNum,$bioRepNum,$techRepNum);
						if ($algoVersion==2) {
							($expCode,$stateNum,$bioRepNum,$techRepNum)=($rawLabel=~/^(\w*)-*State(\d+)-biorep(\d+)-techRep(\d+)/); # "\w-" only for SuperRatio
						}
						else { # v3: A_rep1_techRep1_State1
							#($stateNum,$expCode,$bioRepNum,$techRepNum)=($rawLabel=~/^[\D]+(\d+)_([^_]+)_[\D]+(\d+)_[\D]+(\d+)/); # v3b
							($stateNum,$expCode,$bioRepNum,$techRepNum)=($rawLabel=~/^State(\d+)_([^_]+)_rep(\d+)_techRep(\d+)/); # v3b
						}

						if ($ratioType eq 'SuperRatio' && $stateNum==1) { # Cummulate BioReps across all Exp for Ref SILAC (State1)
							my $nextCode=chr(ord($expCode)+1); # A -> B
							$bioRepCountShift{$nextCode}=$bioRepCountShift{$expCode} if $expCode ne $prevExpCode;
							$bioRepCountShift{$nextCode}++;
							$bioRepNum+=$bioRepCountShift{$expCode};
							$prevExpCode=$expCode;
						}
						my $label=$labelsName{"State$stateNum-biorep$bioRepNum-techRep$techRepNum"} || "*State$stateNum-biorep$bioRepNum-techRep$techRepNum*";
						@{$labels{$label}}=($stateNum,$bioRepNum,$techRepNum);
						$convertLabels{$rawLabel}=$label;
						push @labelArray,$label;
					}
					next;
				}
				my @data=split(/\t/,$_);
				my $rawLabel=(scalar @data > scalar keys %convertLabels)? shift @data : shift @rawLabelList; # 1st col is row label OR no raw label (assume same order than in header)
				foreach my $i (0..$#data) {
					if ($stage eq 'BEFORE') {next if $i < $.-1;} else {next if $i >= $.-1;} # TOP-RIGHT=BEFORE, BOTTOM-LEFT=AFTER
					$labelData{$convertLabels{$rawLabel}}{$labelArray[$i]}=$data[$i];
					$minCorrelation=$data[$i] if ($data[$i] ne 'NA' && $minCorrelation > $data[$i]);
				}
			}
			close CORR;
		}
	}
	elsif ($quantifSoftware =~ /MSstats/) {

	}

	$sthState->finish;
	$sthQuantif->finish;
	$sthAnaInfo->finish;

	$dbh->disconnect;

	###<Starting HTML>###
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
	my $numLabels=scalar keys %labels;
	#my $matSize=($numLabels < 5)? 100 : ($numLabels < 50)?  $numLabels*20 : 1000;
	#JS: width:$matSize,height:$matSize,
	my $cellSize=($numLabels == 2)? 100 : ($numLabels <= 4)? 50 : 25;
	print qq
|var GCM=new heatMap({
div:'globCorrDIV',cellWidth:$cellSize,cellHeight:$cellSize,
//editable:true,
colLabelOrientation:'vertical',
moveLabel:true,
//normalization:{scope:'row',reference:'z-score',colors:'GYR'}
normalization:{scope:'row',reference:'user',limitValues:{ref:0,max:1}} //$minCorrelation
});
|;
	my @sortedLabels=(sort{$labels{$a}[0]<=>$labels{$b}[0] || $labels{$a}[1]<=>$labels{$b}[1] ||$labels{$a}[2]<=>$labels{$b}[2]} keys %labels);

	#>Columns
	#print "GCM.setColumns([['",join("'],['",@sortedLabels),"']]);\n";
	print "GCM.setColumns([";
	my $firstColLabel=1;
	foreach my $colLabel (@sortedLabels) {
		print ',' unless $firstColLabel;
		print "['$colLabel',null,'$colLabel']";
		$firstColLabel=0;
	}
	print "]);\n";

	#>Rows
	foreach my $rowLabel (@sortedLabels) {
		print "GCM.addRow(['$rowLabel',null,'$rowLabel'],[";
		my $firstLabel=1;
		foreach my $colLabel (@sortedLabels) {
			print ',' unless $firstLabel;
			print $labelData{$rowLabel}{$colLabel} if $labelData{$rowLabel}{$colLabel} ne 'NA';
			$firstLabel=0;
		}
		print "]);\n";
	}
	if ($drawGroups) {
		my $groupStrg='[';
		my $prevStatePos=0;
		my $begGrIdx=0;
		foreach my $i (0..$#sortedLabels) {
			my $statePos=$labels{$sortedLabels[$i]}[0];
			if ($i > 0 && $statePos > $prevStatePos) {
				my $endGrIdx=$i-1;
				$groupStrg.=',' if $prevStatePos > 1;
				$groupStrg.="['$states{$prevStatePos}',$prevStatePos,,$begGrIdx,$endGrIdx]";
				$begGrIdx=$i;
			}
			$prevStatePos=$statePos;
		}
		$groupStrg.=",['$states{$prevStatePos}',$prevStatePos,,$begGrIdx,$#sortedLabels]]";
		print qq
|GCM.defineGroups('row',$groupStrg);
GCM.defineGroups('column',$groupStrg);
|;
	}
	print "GCM.draw();\n";
	exit;
}


sub ajaxManageTemplate {
	my $templateAction=param('TEMPLATE');
	my $dbh=&promsConfig::dbConnect;
	print header(-type=>'text/plain',-charset=>'UTF-8'); warningsToBrowser(1);
	if ($templateAction eq 'add') {
		$dbh->do("UPDATE QUANTIFICATION SET QUANTIF_ANNOT=CONCAT(QUANTIF_ANNOT,'::TEMPLATE=1') WHERE ID_QUANTIFICATION=$selQuantifID");
		#print "&nbsp;Yes&nbsp;&nbsp;<INPUT type=\"button\" id=\"remove\" value=\"Remove from templates\" onclick=\"manageTemplate(this.id)\">";
		print '<SPAN class="template">Template&nbsp;&nbsp;<INPUT type="button" class="font11" value="Remove" onclick="manageTemplate(\'remove\')"></SPAN>';
	}
	elsif ($templateAction eq 'remove') {
		$dbh->do("UPDATE QUANTIFICATION SET QUANTIF_ANNOT=REPLACE(QUANTIF_ANNOT,'::TEMPLATE=1','') WHERE ID_QUANTIFICATION=$selQuantifID");
		#print "&nbsp;No&nbsp;&nbsp;<INPUT type=\"button\" id=\"add\" value=\"Save as template\" onclick=\"manageTemplate(this.id)\">";
		print '<INPUT type="button" class="title3" value="Save as template" onclick="manageTemplate(\'add\')">';
	}
	$dbh->commit;
	$dbh->disconnect;
	exit;
}

sub dbEncodeModifPosition {
	my $modStrg=$_[0];
	my $dbModStrg='';
	if ($modStrg=~/~/) { # ambiguity format => convert to DB output format
		my %pos;
		my ($start)=($modStrg=~/^(\d+)/); $pos{$start}='-';
		my ($end)=($modStrg=~/~(\d+)/); $pos{$end}='+';
		my ($mods,$sites)=($modStrg=~/(\d+)\/(\d+)/); $pos{$sites}=$mods;
		foreach my $p (sort{$a<=>$b} keys %pos) {
			$dbModStrg.='.' if $dbModStrg;
			$dbModStrg.=$pos{$p}.$p;
		}
	}
	else {$dbModStrg=$modStrg;}
	return $dbModStrg;
}

sub revertLog2 {
	my $l2Val=shift;
	my $fc=exp($l2Val*log(2));
	my $formatFc=($fc > 1)? sprintf "%.2f",$fc : '1/'.sprintf "%.2f",1/$fc;
	$formatFc=~s/\.0+\Z//;
	$formatFc=1 if $formatFc eq '1/1';
	return $formatFc;
}


####>Revision history<####
# 2.4.5 [Fix] Bug in &ajaxRestrictProteinList (PP 11/09/18)
# 2.4.4 Chart highlighting now possible from custom lists (PP 07/09/18)
# 2.4.3 Clean protein alias from badly parsed characters in case identifier conversion (PP 09/07/18)
# 2.4.2 Updated quanti selection in ana call for compatibility with TMT and DIA. TODO: Remove usage of PEPTIDE_QUANTIFICATION table in peptide view for 'ratio' (PP 27/06/18)
# 2.4.1 [Fix] bug display of negative normalization proteins list & added = in pepID string split (PP 15/06/18)
# 2.4.0 Reads peptide quantif data from file(s) not DB (PP 11/05/18)
# 2.3.8 [Fix] Minor JS bug when displaying ratios distribution (PP 07/05/18)
# 2.3.7 [Fix] Proper image file selection for label-free v2 vs v3 algo (PP 16/04/18)
# 2.3.6 Modif to save a quanti as template (MLP 16/04/18)
# 2.3.5 Minor change is text display (PP 06/04/18)
# 2.3.4 Quantifications with myProMS algo v3 are now exportable in XLS<BR>TODO: Export MaxQuant (PP 16/02/18)
# 2.3.3 Minor update to handle NAs in algo v3 correlation matrix (PP 13/12/17)
# 2.3.2 Fix minor bug in &ajaxPeptideRatios for algo v2 %seenExcluded dataSrc (PP 11/12/17)
# 2.3.1 Uses color-codes for peptide source in peptide view for algo v3 (04/12/17)
# 2.3.0 Compatible with quantif algo myProMS v3 (PP 27/11/17)
# 2.2.3 Modifs to display DIA quantification from MLP (PP 16/11/17)
# 2.2.2 Now uses &promsQuantif::fetchSiteList, update in &ajaxManageSaveProteins function & auto-reload for on-going quantification (PP 11/08/17)
# 2.2.1 Filters sites if restriction list and quantif are sites (PP 07/08/17)
# 2.2.0 Compatible with MODIFICATION_SITE, improved BioRep naming & include GA 2.1.4b update (PP 03/08/17)
# 2.1.4 Updated Excel colors (PP 22/02/17)
# 2.1.3 Minor modification (MLP 26/01/17)
# 2.1.2 Improved handling of MaxQuant quantif (PP 20/01/17)
# 2.1.1 Fill missing format for 2 merge_range functions in XLS export (GA 20/12/16)
# 2.1.0 Major changes in MaxQuant data display<BR>TODO: Handle export (PP 13/12/16)
# 2.0.2 Minor modification for emPAI (GA 25/11/16)
# 2.0.1 Bug fix in SWATH fragment display for multi-ratio quantif (PP 26/10/16)
# 2.0.0 Reverse ratios handling & Coef. variation filtering & SWATH peptide display (PP 19/10/16)
# 1.8.4 Minor fix to keep +/-inf flag in log2 ratios during export (PP 10/08/16)
# 1.8.3 Uses SINGLE_REF flag for multi-ratio quantifs (PP 04/08/16)
# 1.8.2 Handles SWATH & now used by 'summary' call of protein ratio quantifs (PP 02/08/16)
# 1.8.1 Minor modification for MaxQuant (GA 01/08/16)
# 1.8.0 Bug fix for rare cases of missing ratio due to exclusion of all peptides by global CV outlier detection (PP ../05/16)
# 1.7.9 Display MQ quantifications in myProMS (intensity, iBAQ, LFQ, MS/MS count) (GA 09/05/16)
# 1.7.8 Compatibility with peptidePlot.js v1.1.0 & bug fix for raw peptide data display of iTRAQ & other minor updates (PP 28/04/16)
# 1.7.7 Added matching peptides option for label-free quantification (PP 25/02/16)
# 1.7.6 Minor bug fix in &ajaxListSelectedProteins for label-free (PP 19/02/16)
# 1.7.5 Minor change to prevent undefined quantif name (PP 11/02/16)
# 1.7.4 Bug fix in peptide display of internal quantif & +/-inf conversion limit for volcano only (PP 02/02/16)
# 1.7.3 Bug fix in $sthPR->execute for internal quantif (PP 30/01/16)
# 1.7.2 Compatible with quantifSet XIC summing (PP 19/01/16)
# 1.7.1 Compatible with old FC quantifs (PP 15/10/15)
# 1.7.0 Compatible with PTM quantification & called openProject v1.5.9 (PP 13/05/15)
# 1.6.5 display goTerm from child Go Analysis and pathway term from Pathway analysis, add revert2log function (SL 05/11/14)
# 1.6.4 Minor bug correction numStates for internal quantifs (PP 23/09/14)
# 1.6.3 Uses &promsMod::printAjaxManageSaveProteins function for List save management (PP 25/08/14)
# 1.6.2 Exploratory Analysis constraint on List deletion & fix in R image display (PP 25/08/14)
# 1.6.1 Update for EMPAI in export and printing (GA 26/06/14)
# 1.6.0 Support for SuperRatio, default view is now List & change "Category" to "List" (PP 19/06/14)
# 1.5.1 Default display is now 3 peptide/protein (PP 29/04/14)
# 1.5.0 More links to ratios statistics (PP 14/04/14)
# 1.4.9 Add log2 option to compare XIC extractions when call=analysis<BR>TODO: change color definition (GA 24/03/14)
# 1.4.8 Uses rmtree instead of remove_tree (PP 10/03/14)
# 1.4.7 ratio display management for internal quantif & graph view merge for design/internal quantif (PP 27/02/13)
# 1.4.6 Export to Excel (PP 20/02/14)
# 1.4.5 Export SVG & improvement of ajaxPepRatio action (PP 11/02/14)
# 1.4.4 system command removal (PP 12/11/13)
# 1.4.3 Displays minimum infinite ratios choice (PP 21/10/13)
# 1.4.2 Improved handling of infinite ratios (PP 09/10/13)
# 1.4.1 Minor modification to avoid to skip $labelingInfo information (GA 23/09/13)
# 1.4.0 Generated from listAnaQuantifications.cgi 1.3.3PP2 (PP 02/09/13)
# 1.3.3PP2 Volcano plot with missing p-values set to 1 & infinite ratios detection (PP 21/08/13)
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
