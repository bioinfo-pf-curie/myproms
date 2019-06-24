#!/usr/local/bin/perl -w

################################################################################
# selAna4Quantification.cgi      3.1.0                                         #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
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
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use promsQuantif;

#print header(-charset=>'UTF-8'); warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my %msTypeName=&promsConfig::getMsType;
my $userID=$ENV{'REMOTE_USER'};
$msTypeName{'MIS'}="MS/MS"; #redef, to keep space
my %quantifProcesses=('XICMCQ'=>'XIC extraction','EMPAI'=>'emPAI quantification','SIN'=>'SI<SUB>N</SUB> quantification','XICCORR'=>'Isobaric correction'); #,'SILAC'=>'SILAC','ITRAQ'=>'iTRAQ','TMT'=>'TMT','DESIGN'=>'Protein-Ratio'); #,'MCQ'=>'Ext. ion chrom.''TnPQ or Peptide ratio'
#my $updateFrameString="";
my $maxLabPerChannel=10; # max num for isotope extraction...

###############################
####>Recovering parameters<####
###############################
my $quantifType=uc(param('quantifType')); # can be XICMCQ,EMPAI,SIN,XICCORR
my $branchID=param('ID');
my ($item,$itemID)=split(':',$branchID);
$item=lc($item);
my $ITEM=uc($item);

############################
####>form was submitted<####
############################
if (param('launch')) {
	&launchQuantifications;
	exit;
}

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;
my $projectID=&promsMod::getProjectID($dbh,$itemID,$item);
my ($itemName)=$dbh->selectrow_array("SELECT NAME FROM $ITEM WHERE ID_$ITEM=$itemID");
my $selItems=($quantifType=~/XIC/)? 'Peptides Quantifications' : 'Analyses';
my $titleString="Select $selItems in ".&promsMod::getItemType($ITEM)." <FONT color=#DD0000>$itemName</FONT><BR>for $quantifProcesses{$quantifType}";

########################################################
####>Recovering list of experiment, sample analysis<####
########################################################
my (%listDataBank,@itemAnalyses,%anaProteins,%anaLabelMods,%anaLabeling,%anaPeptideQuantifs); #%modifications,

####>Recovering DBs name<####
my $sthAD = $dbh->prepare("SELECT D.ID_DATABANK,NAME FROM ANALYSIS_DATABANK AD,DATABANK D WHERE AD.ID_DATABANK=D.ID_DATABANK AND AD.ID_ANALYSIS=?");

####>Recovering quantitation name<####
#my $sthQuanti = $dbh->prepare("SELECT ID_QUANTIF_METHOD,NAME FROM QUANTIF_METHOD");
#$sthQuanti->execute;
#while (my ($idQuanti,$quantiName)= $sthQuanti->fetchrow_array) {
#	$listQuantif{$idQuanti}=$quantiName;
#}
#$sthQuanti->finish;

my @sthItem;
my $baseFieldString='ID_ANALYSIS,VALID_STATUS,MS_TYPE,DATA_FILE,FILE_FORMAT,WIFF_FILE,TAXONOMY,MAX_RANK,MIN_SCORE,INSTRUMENT,0,LABELING';
if ($ITEM eq 'ANALYSIS') {
	$sthItem[0]=$dbh->prepare("SELECT $baseFieldString,NAME FROM ANALYSIS WHERE ID_ANALYSIS=$itemID ORDER BY DISPLAY_POS ASC");
}
elsif ($ITEM eq 'SAMPLE') {
	$sthItem[0]=$dbh->prepare("SELECT $baseFieldString,NAME FROM ANALYSIS WHERE ID_SAMPLE=$itemID ORDER BY DISPLAY_POS ASC");
}
elsif ($ITEM eq 'SPOT') {
	$sthItem[0]=$dbh->prepare("SELECT $baseFieldString,ANALYSIS.NAME FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND ID_SPOT=$itemID ORDER BY ANALYSIS.DISPLAY_POS ASC");
}
elsif ($ITEM eq 'GEL2D') {
	$sthItem[0]=$dbh->prepare("SELECT $baseFieldString,'SPOT',SPOT.NAME,ANALYSIS.NAME FROM ANALYSIS,SAMPLE,SPOT WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT=SPOT.ID_SPOT AND ID_GEL2D=$itemID ORDER BY SPOT.NAME ASC, ANALYSIS.DISPLAY_POS ASC");
}
elsif ($ITEM eq 'EXPERIMENT') {
	$sthItem[0]=$dbh->prepare("SELECT $baseFieldString,'GEL2D',GEL2D.NAME,'SPOT',SPOT.NAME,ANALYSIS.NAME FROM ANALYSIS,SAMPLE,SPOT,GEL2D WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT=SPOT.ID_SPOT AND SPOT.ID_GEL2D=GEL2D.ID_GEL2D AND GEL2D.ID_EXPERIMENT=$itemID ORDER BY GEL2D.DISPLAY_POS ASC, SPOT.NAME ASC, ANALYSIS.DISPLAY_POS ASC");
	$sthItem[1]=$dbh->prepare("SELECT $baseFieldString,'SAMPLE',SAMPLE.NAME,ANALYSIS.NAME FROM ANALYSIS,SAMPLE WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND SAMPLE.ID_SPOT IS NULL AND ID_EXPERIMENT=$itemID ORDER BY SAMPLE.DISPLAY_POS ASC, ANALYSIS.DISPLAY_POS ASC");
}
my $sthVVP=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE VISIBILITY>0 AND ID_ANALYSIS=?");
my $sthAVP=$dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=?");
my $sthSTP=$dbh->prepare("SELECT COUNT(*) FROM PROTEIN_VALIDATION WHERE SEL_STATUS>=1 AND ID_ANALYSIS=?");
my $sthATP=$dbh->prepare("SELECT COUNT(*) FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=?");
my $sthAM=$dbh->prepare("SELECT AM.ID_MODIFICATION,AM.SPECIFICITY,IS_LABEL,M.VALID_STATUS,M.MONO_MASS,M.PSI_MS_NAME,M.INTERIM_NAME,M.DES,M.SYNONYMES FROM ANALYSIS_MODIFICATION AM,MODIFICATION M WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND ID_ANALYSIS=?"); # isobaric modif are FIXED!
my $sthPLQ=$dbh->prepare("SELECT Q.ID_QUANTIFICATION,Q.NAME,Q.QUANTIF_ANNOT FROM QUANTIFICATION Q,ANA_QUANTIFICATION A WHERE Q.ID_QUANTIFICATION=A.ID_QUANTIFICATION AND A.ID_ANALYSIS=? AND Q.FOCUS='peptide' AND Q.STATUS>=1 AND Q.ID_PRODUCT IS NULL");

foreach my $sth (@sthItem) {
	$sth->execute;
	while (my @anaData=$sth->fetchrow_array) {
		my $anaID=$anaData[0];
		$anaData[11]='' unless $anaData[11]; # LABELING
		#>Databank(s)
		$sthAD->execute($anaID);
		my @dbUsed;
		while (my ($dbID,$dbName)=$sthAD->fetchrow_array) {
			push @dbUsed,$dbID;
			$listDataBank{$dbID}=$dbName;
		}
		$anaData[10]=\@dbUsed;
		push @itemAnalyses,\@anaData;
		if ($anaData[1]>=1) { # valid proteins
			$sthVVP->execute($anaID);
			push @{$anaProteins{$anaID}},$sthVVP->fetchrow_array;
			$sthAVP->execute($anaID);
			push @{$anaProteins{$anaID}},$sthAVP->fetchrow_array;
		}
		else {# non-validated proteins
			$sthSTP->execute($anaID);
			push @{$anaProteins{$anaID}},$sthSTP->fetchrow_array;
			$sthATP->execute($anaID);
			push @{$anaProteins{$anaID}},$sthATP->fetchrow_array;
		}
		$sthAM->execute($anaID);
		$anaLabeling{$anaID}='FREE'; # default
		while (my ($modID,$anaSpec,$isLabel,$isValid,$monoMass,$psiName,$interName,$des,$synName)=$sthAM->fetchrow_array) {
			if ($isLabel) {
				$psiName='' unless $psiName; $interName='' unless $interName; $des='' unless $des; $synName='' unless $synName;
				@{$anaLabelMods{$modID}}=($anaSpec,$psiName,$interName,$isValid,$monoMass);
				my $allNames=$psiName.'##'.$interName.'##'.$des.'##'.$synName;
				$anaLabeling{$anaID}=($allNames=~/iTRAQ.+plex/i)? 'ITRAQ' : ($allNames=~/TMT.+plex/i)? 'TMT' : ($allNames=~/Label:|SILAC/i)? 'SILAC' : 'OTHER';
				last;
			}
			#else {
			#	##$modifications{$modID}=$psiName || $interName || $des;
			#	##unless ($modifications{$modID}) {
			#	##	$synName=~s/^##//; $synName=~s/##$//; $synName=~s/##/, /;
			#	##	$modifications{$modID}=$synName;
			#	##}
			#	$anaLabeling{$anaID}='FREE';
			#}
		}
		
		##>Peptide quantifications
		if ($quantifType eq 'XICCORR') { # && $anaData[11]=~/ITRAQ|TMT/i
			$sthPLQ->execute($anaID);
			while (my ($quantifID,$quantifName,$quantifAnnot)=$sthPLQ->fetchrow_array) {				
				next if (!$quantifAnnot || $quantifAnnot !~ /^LABEL=(ITRAQ|TMT)/i || $quantifAnnot=~/::CORRECTION=/); # A corrected quantif already exists
				push @{$anaPeptideQuantifs{$anaID}},[$quantifID,$quantifName];
			}
		}
	}
	$sth->finish;
}
$sthAD->finish;
foreach my $sth (@sthItem) {$sth->finish;}
$sthVVP->finish;
$sthAVP->finish;
$sthSTP->finish;
$sthATP->finish;
$sthAM->finish;
$sthPLQ->finish;

my %isotopicDistributions;
if ($quantifType eq 'XICCORR') {
	my $sthIsoDis=$dbh->prepare("SELECT ID_PRODUCT,LABEL_TYPE,NAME,PRODUCT_NUMBER,LOT_NUMBER FROM ISOTOPIC_CORRECTION WHERE USE_STATUS='yes'");
	$sthIsoDis->execute;
	while (my ($productID,$labelType,@info)=$sthIsoDis->fetchrow_array) {
		my ($label,$plex)=split(':',$labelType);
		@{$isotopicDistributions{$label}{$productID}}=(@info,$plex);
	}
	$sthIsoDis->finish;
}

$dbh->disconnect;

#######################
####>Starting HTML<####
#######################
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
print header(-charset=>'utf-8'); # because of ° in source selection
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Select Analyses for Non-Design Quantifications</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TD.center {text-align:center}
</STYLE>
<SCRIPT type="text/javascript">
|;
if ($quantifType eq 'XICMCQ') {
	###> For label purporse MCQXIC extraction
	print qq
|//For label purporse MCQXIC extraction
var refValue=-1;
var analabels=[];
analabels[0]=[-1,'Light',0.0000,''];
|;
	#print "var analabels=new Array();\n";
	#print "analabels[0]=new Array();\n";
	#print "analabels[0][0]=-1;\n";
	#print "analabels[0][1]=\"No Label\";\n";
	#print "analabels[0][2]=0.0000;\n";
	#print "analabels[0][3]=\"\";\n";
	my $labelIndx=0;
	foreach my $modID (keys %anaLabelMods) {
		$labelIndx++;
		my ($anaSpec,$psiMsName,$interimName,$isValid,$monoMass)=@{$anaLabelMods{$modID}};
		print "analabels[$labelIndx]=[$modID,'$interimName',$monoMass,'$anaSpec'];\n";
		#print "analabels[$labelIndx]=new Array();\n";
		#print "analabels[$labelIndx][0]=$modID;\n";
		#print "analabels[$labelIndx][1]=\"$interimName\";\n";
		#print "analabels[$labelIndx][2]=$monoMass;\n";
		#print "analabels[$labelIndx][3]=\"$anaSpec\";\n";
	}
	print qq
|function updateLabels(modID,channelNum,labelNum){
	for(var i=0; i< analabels.length; i++) {
		if(analabels[i][0] == modID) {
			var name=document.getElementById('lab'+labelNum+':name_channel'+channelNum);
			name.value=analabels[i][1];
			// var dm=document.getElementById('lab'+labelNum+':dm_channel'+channelNum);
			// dm.value=analabels[i][2];
			var sp=document.getElementById('lab'+labelNum+':sp_channel'+channelNum);
			sp.value=analabels[i][3];
		}
	}
}
function updateChargeState(xicType) {
	var checkBoxField=document.getElementById('allChargeStates');
	if (xicType == 'max'){checkBoxField.disabled=false;}
	else{checkBoxField.disabled=true;}
}
function updateTextState(name,target) {
	var textField=document.getElementById(name);
	if (target == 'ANY'){textField.disabled=false;}
	else{textField.disabled=true;}
}
function addQuantificationLabel(channelNum,labelNum,action) {
	if (action == 'show') {
		if(labelNum < $maxLabPerChannel) {
			var newLabel=labelNum+1;
			// Show-Hide button
			document.getElementById('show:'+channelNum+':'+labelNum).style.display='none';
			document.getElementById('hide:'+channelNum+':'+labelNum).style.display='none';
			document.getElementById('hide:'+channelNum+':'+newLabel).style.display='';
			if (newLabel == $maxLabPerChannel) {
				document.getElementById('show:'+channelNum+':'+newLabel).style.display='none';
			}
			// Show new label fieldset
			document.getElementById('fs:'+channelNum+':'+newLabel).style.display='';
		}
	}
	else {
		var newLabel=labelNum-1;
		document.getElementById('fs:'+channelNum+':'+labelNum).style.display='none';
		document.getElementById('show:'+channelNum+':'+newLabel).style.display='';
		if (newLabel > 1) {document.getElementById('hide:'+channelNum+':'+newLabel).style.display='';}
	}
}
function updateRefAna(reference){
   var anaBox=document.selAnaForm.anaList;
   // Uncheck last reference analysis
   if (refValue > 0 && anaBox.length) {
		for (var i=0; i < anaBox.length; i++){
		    if(anaBox[i].value==refValue){
			    anaBox[i].checked=false;
		    }
		}
   }
   // Check reference analysis
   if(anaBox.length) {// more than 1 checkboxes
		for (var i=0; i < anaBox.length; i++){
		    if(anaBox[i].value==reference.value){
			    anaBox[i].checked=true;
			    refValue=anaBox[i].value;
		    }
        }
   }
   else{
		anaBox.checked=true;
   }

}
function updateSettings(act) {
	if (act=='more') {
		document.getElementById('moreSettings').style.display='none';
		document.getElementById('lessSettings').style.display='';
		document.getElementById('advancedSetDIV').style.display='';
	}
	else {
		document.getElementById('moreSettings').style.display='';
		document.getElementById('lessSettings').style.display='none';
		document.getElementById('advancedSetDIV').style.display='none';
	}
}
function showMZRange(){
	var chButton=document.getElementById('allChargeStates');
	if(chButton.checked){
		document.getElementById('massrangeDIV').style.display='';
	}
	else{
		document.getElementById('massrangeDIV').style.display='none';
	}
}
function showParameters(alignType) {
	if (alignType=='OBI') {
		document.getElementById('paramObiwarp').style.display='';
		document.getElementById('paramMs2').style.display='none';
	}
	else if (alignType=='MS2'){
		document.getElementById('paramObiwarp').style.display='none';
		document.getElementById('paramMs2').style.display='';
	}
}
function updateLabeling(labeling) {
	if (labeling=='SILAC'){
		document.getElementById('paramLabel').style.display='';
		document.getElementById('paramAlign').style.display='none';
	}
	else if (labeling=='FREE'){
		document.getElementById('paramLabel').style.display='none';
		document.getElementById('paramAlign').style.display='';
	}
	//Filtering Reference Analysis
	var refAnaSelect=document.getElementById('refAna_alignment');
	for (let i=1; i < refAnaSelect.options.length; i++) { //[0]= '-= Select =-'
		refAnaSelect.options[i].disabled=(refAnaSelect.options[i].getAttribute('data-labeling')==labeling)? false : true;
	}
	if (refAnaSelect.options[refAnaSelect.selectedIndex].disabled) {refAnaSelect.selectedIndex=0;} //deselect if disabled
	//Filtering Checkable analyses
	var anaBox=document.selAnaForm.anaList;
	if (anaBox.length) {
		for (let i=0; i < anaBox.length; i++) {
			var disStatus=(anaBox[i].getAttribute('data-labeling')==labeling)? false : true;
			anaBox[i].disabled=disStatus;
			document.getElementById('file_'+anaBox[i].value).disabled=disStatus;
		}
	}
	else {
		var disStatus=(anaBox.getAttribute('data-labeling')==labeling)? false : true;
		anaBox.disabled=disStatus;
		document.getElementById('file_'+anaBox.value).disabled=disStatus;
	}
}
function getRadioVal() {
  var rads = document.getElementsByName('extractL');

  for(var rad in rads) {
    if(rads[rad].checked)
      return rads[rad].value;
  }

  return null;
}
|;
}
elsif ($quantifType eq 'XICCORR') {
	print qq
|window.onload=function() {updateCorrectionLabeling(document.selAnaForm.product);}
function updateCorrectionLabeling(prodSel) {
	var corrLabeling=prodSel.options[prodSel.selectedIndex].getAttribute('data-labeling');
	var anaBox=document.selAnaForm.anaList;
	if (!anaBox) return 0; // no selectable analyses
	
	if (anaBox.length) { // more than 1 checkboxes
		for (let i=0; i < anaBox.length; i++) {
			anaBox[i].disabled=(anaBox[i].getAttribute('data-labeling') == corrLabeling)? false : true;
		}
	}
	else {
		anaBox.disabled=(anaBox.getAttribute('data-labeling') == corrLabeling)? false : true;
	}
}
|;
}
print qq
|function testCheckbox() {
	var anaBox=document.selAnaForm.anaList;
	if (!anaBox) return 0; // no selectable analyses
	var selected=0;
	if (anaBox.length) { // more than 1 checkboxes
		for (i=0; i < anaBox.length; i++){
			if (!anaBox[i].disabled && anaBox[i].checked==true) {
				if (!testQuantifSelection(anaBox[i].value)) {
					selected=-1;
					break;
				}
				selected=1;
			}
		}
	}
	else if (!anaBox.disabled && anaBox.checked==true){
		selected=(!testQuantifSelection(anaBox.value))? -1 : 1;
	}
	return selected;
}
function checkall(checkStatus){
	var anaBox=document.selAnaForm.anaList;
	if (!anaBox) return; // no selectable analyses
	if (anaBox.length) { // more than 1 checkbox
		for (let i=0; i < anaBox.length; i++) {
			anaBox[i].checked=checkStatus;
			//updateQuantifSelection(anaBox[i]);
		}
	}
	else { // Only 1 checkbox
		anaBox.checked=checkStatus;
		//updateQuantifSelection(anaBox);
	}
}
function testQuantifSelection(anaID) {return true;}
//function updateQuantifSelection(chkbox) {}
function cancelAction() {
	//top.promsFrame.selectedAction='summary'; // set default action to 'summary'
	top.promsFrame.optionFrame.selectOption();
}
function checkForm(myForm) {
	if (testCheckbox()==0) {alert('ERROR: No Analyses selected!'); return false;}
	else if (testCheckbox()==-1) {alert('ERROR: Peptide quantification not selected for at least 1 Analysis!'); return false;} // only for SILAC,ITRAQ
|;
if ($quantifType eq 'XICMCQ') {
	print qq
|	var anaBox=myForm.anaList;
	var missingFile=false;
	if (anaBox.length) { // more than 1 checkbox
		for (let i=0; i < anaBox.length; i++){
			if (!anaBox[i].disabled && anaBox[i].checked==true && !document.getElementById('file_'+anaBox[i].value).value.match('mzX*ML')) {
				missingFile=true;
				break;
			}
		}
	}
	else if (!anaBox.disabled && anaBox.checked==true && !document.getElementById('file_'+anaBox.value).value.match('mzX*ML')) {
		missingFile=true;
	}
	if (missingFile) {
		alert('ERROR: Missing or not-mzXML/mzML file detected.');
		return false;
	}
	if (!myForm.refAna_alignment.value && getRadioVal()=='NO'){
		alert('ERROR: No reference selected for alignment.');
		return false;
	}
	if (!myForm.refAna_alignment.value){
		alert('ERROR: No reference selected.');
		return false;
	}
|;
}
elsif ($quantifType eq 'XICCORR') {
	print qq
|	if (!myForm.product.value){
		alert('ERROR: No Isotopic distribution selected.');
		return false;
	}
	if (!myForm.keep_old.checked && !confirm('WARNING: Original Peptide Quantifications will be deleted. Proceed?')) {
		return false;
	}
|;
}
print qq
|	return false;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif" >
<CENTER>
<FONT class="title">$titleString</FONT>
<BR>
<FORM name="selAnaForm" method="post" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="ID" value="$branchID">
<INPUT type="hidden" name="quantifType" value="$quantifType">
<BR>
|;

############################################
####>Displaying MassChroQ-specific form<####
############################################
if ($quantifType eq 'XICMCQ') {
	print qq
|<TABLE bgcolor=$darkColor>
<TR><TH class="title2" align=right>Name :</TH><TD bgcolor=$lightColor><INPUT type="text" name="quantifName" value="$quantifProcesses{$quantifType} extraction" class="title3" style="width:400px"/></TD></TR>
<TR>
<TH align=right nowrap valign=top>Raw-data settings :</TH>
<TD bgcolor=$lightColor nowrap>&nbsp;<B>Extraction type:</B>
<SELECT name="rawdata_extraction">
<OPTION value="centroid">Centroid</OPTION>
<OPTION value="profile">Profile</OPTION>
</SELECT>
<FONT class="font11">&nbsp;(for mzXML)</FONT>
</TD>
</TR>
<TR>
<TH align=right nowrap valign=top>Isotope labeling :</TH>
<TD bgcolor=$lightColor nowrap>&nbsp;<B><SELECT name="extractL" value="YES" onchange="updateLabeling(this.value)"><OPTION value=\"FREE\">NONE</OPTION><OPTION value=\"SILAC\">SILAC</OPTION></SELECT></B>
|;
	###> Label informations
	my $optMods="<OPTION value=\"\">-= Select =-</OPTION>\n";
	$optMods.="<OPTION value=\"-1\"> Light / +0.0000 Da</OPTION>\n";
	foreach my $modID (keys %anaLabelMods) {
		my ($anaSpec,$psiMsName,$interimName,$isValid,$monoMass)=@{$anaLabelMods{$modID}};
		my $isDisabled=($isValid)? "" : "disabled";
		my $monoMassStg=($monoMass>0)?  "+$monoMass Da" : "$monoMass Da";
		$optMods.="<OPTION value=\"$modID\" $isDisabled>$psiMsName / $monoMassStg</OPTION>\n";
	}
	#print "<TH align=right nowrap valign=top></TH><TD bgcolor=$lightColor nowrap><B>\n";
	print "<DIV id=\"paramLabel\" style=\"display:none\"><B>\n";
	#my @colors=("#80ACFF","#5195FF","#428CFF");
	foreach my $channel (1..3) {
		print "<BR>\n" if $channel > 1;
		print "<FIELDSET style=\"border:3px groove threedface;\"><LEGEND>Name for channel #$channel:<INPUT type=\"text\" id=\"name_channel$channel\" name=\"name_channel$channel\" value=\"\" style=\"width:250px\"></LEGEND>\n";
		#my $color=$colors[$channel-1];
		foreach my $label (1..$maxLabPerChannel) { # In triplex SILAC,channel can be composed of several labels: Heavy=Arg10 and Lys8 ; Medium=
			my $state=($label==1)? '' : 'none';
			print qq
|<DIV id="fs:$channel:$label" style="display:$state">
<FIELDSET>
<LEGEND>Quantification Label:</LEGEND>
<TABLE width=100% bgcolor=$darkColor>
<TR><TH align=right nowrap valign=middle>Label Name:</TH><TD><INPUT type="text" id="lab$label:name_channel$channel" name="lab$label:name_channel$channel" value="" style="width:250px"></TD></TR>
<TR><TH align=right nowrap valign=middle>Modification target:</TH><TD><SELECT id="lab$label:tg_channel$channel" name="lab$label:tg_channel$channel" onchange="updateTextState('lab$label:sp_channel$channel',this.value)"><OPTION value="ANY">Side chain</OPTION><OPTION value="Nter">N-Ter</OPTION><OPTION value="Cter">C-Ter</OPTION></SELECT></TD></TR>
<TR><TH align=right nowrap valign=middle>Modification:</TH><TD><SELECT id="lab$label:channel$channel" name="lab$label:channel$channel" onchange="updateLabels(this.value,$channel,$label)">$optMods</SELECT> on <INPUT type="text" id="lab$label:sp_channel$channel" name="lab$label:sp_channel$channel" value="" size="2"></TD></TR>
</TABLE>
</FIELDSET>
<INPUT id="show:$channel:$label" type="button" value="Add quantification label" onclick="addQuantificationLabel($channel,$label,'show')" style="font-size:11px"/><INPUT id="hide:$channel:$label" type="button" value="Remove quantification label" onclick="addQuantificationLabel($channel,$label,'hide')" style="display:none; font-size:11px;"/>
</DIV>
|;
		}
		#print "<TR><TD nowrap>&nbsp;Quan Channel $channel:<SELECT id=\"channel$channel\" name=\"channel$channel\" onchange=\"updateLabels(this.value,$channel)\">$optMods</SELECT>\n";
		#print "<";
		#print "&nbsp;Channel name:<INPUT type=\"text\" id=\"name_channel$channel\" name=\"name_channel$channel\" value=\"\" size=\"6\">";
		#print "&nbsp;Delta-Mass: <INPUT type=\"text\" id=\"dm_channel$channel\" name=\"dm_channel$channel\" value=\"\" size=\"6\">";
		#print "&nbsp;Specificity: <INPUT type=\"text\" id=\"sp_channel$channel\" name=\"sp_channel$channel\" value=\"\" size=\"2\"></TD></TR>\n";
		print "\n</FIELDSET>\n";
	}
	print "</B></DIV>";
	print "</TD></TR>\n";
	print qq
|<DIV id="paramAlign">
<TR>
<TH align=right nowrap valign=top>Alignment settings :</TH><TD bgcolor=$lightColor><TABLE cellpadding=0>
	<TR><TD nowrap>&nbsp;<B>Alignment algorithm:</B>
	<SELECT name="alignment_method" onchange="showParameters(this.value)">
	<OPTION value="MS2">ms2</OPTION>
	<OPTION value="OBI">OBI-Warp</OPTION>
	</SELECT>
	&nbsp;&nbsp;&nbsp;<B>Reference:</B>
	<SELECT name="refAna_alignment" id="refAna_alignment" required">
	<OPTION value="">-= Select =-</OPTION>
|;
	foreach my $refAnaData (@itemAnalyses) {##> Option of reference for alignment
		my ($anaID,$valStat,$msType,$dataFile,$fileFormat,$wiffFile,$taxonomy,$maxRank,$minScore,$instrument,$refDbUsed,$labelStrg,@projHierarchy)=@{$refAnaData};
		my $disabStrg=($anaLabeling{$anaID} eq 'FREE')? '' : ' disabled';
		print "	<OPTION value=\"$anaID\" data-labeling=\"$anaLabeling{$anaID}\" onclick=\"updateRefAna(this);\"$disabStrg>$projHierarchy[-1]</OPTION>\n" unless $valStat<1;
	}
	print qq
|	</SELECT>
	</TD>
	</TR>
	<TD bgcolor=$lightColor nowrap valign=top>
		<DIV id="paramObiwarp" style="display:none">
		<TABLE cellpadding=0 cellspacing=0>
		<TR>
		<TD nowrap>&nbsp;<B>Align from <INPUT type="text" name="mz_start" value="400" size="2"> to <INPUT type="text" name="mz_stop" value="1200" size="3"> m/z window</B></TD>
		</TR>
		</TABLE>
		</DIV>
		<DIV id="paramMs2">
		<TABLE cellpadding=0 cellspacing=0>
		<TR><TD nowrap>&nbsp;<B>Tendency <INPUT type="text" name="ms2_tendency" value="10" size="2"></TR>
		<TR><TD nowrap>&nbsp;<INPUT type="checkbox" name="ms2Smooth" value="1" checked><B>MS2 smoothing <INPUT type="text" name="ms2_smoothing" value="5" size="1"></TD></TR>
		<TR><TD nowrap>&nbsp;<INPUT type="checkbox" name="ms1Smooth" value="1" checked><B>MS1 smoothing <INPUT type="text" name="ms1_smoothing" value="3" size="1"></B></TD>
		</TR>
		</TABLE>
		</DIV>
	</TD>
</TABLE>
</TD>
</TR>
</DIV>
<TR>
<TH align=right nowrap valign=top>Peptide selection :</TH>
<TD bgcolor=$lightColor nowrap>
	<TABLE cellpadding=0 cellspacing=0>
	<TD bgcolor=$lightColor nowrap><INPUT type="checkbox" name="allChargeStates" id="allChargeStates" value="1" onchange="showMZRange()"><B>Extract all charge states of the peptides (even if no MS/MS exists for it)</B>&nbsp;
	</TABLE>
	<DIV id="massrangeDIV" style="display:none">
	<TABLE cellpadding=0 cellspacing=0>
	<TR>
	<TD nowrap>&nbsp;<B>Range (m/z) <INPUT type="text" name="mzRange_start" value="400" size="2"> - <INPUT type="text" name="mzRange_stop" value="1200" size="3"></B></TD>
	</TR>
	</TABLE>
	</DIV>
</TD>
</TR>
<TR>
<TH align=right nowrap valign=top>Extract XIC traces :</TH>
<TD bgcolor=$lightColor nowrap>&nbsp;<B>No<INPUT type="radio" name="traces" id="traces" value="0" checked>Yes<INPUT type="radio" name="traces" id="traces" value="1"></B>
</TD>
</TR>
<TR>
<TH align=right nowrap valign=top>&nbsp;Quantification settings :</TH><TD bgcolor=$lightColor nowrap valign=top>
	&nbsp;<B>Type of XIC:</B><SELECT name="XIC_type" onchange="updateChargeState(this.value)"><OPTION value="sum">TIC XIC</OPTION><OPTION value="max">BasePeak XIC</OPTION></SELECT>
	<BR><INPUT type="button" id="moreSettings" class="font11" value="More settings" onclick="updateSettings('more')"/><INPUT type="button" id="lessSettings" class="font11" value="Less settings" style="display:none" onclick="updateSettings('less')"/>
	<DIV id="advancedSetDIV" style="display:none">
	&nbsp;&bull;<B><U>Advanced settings</U>:</B>
	<TABLE cellpadding=0 cellspacing=0>
	<TR>
	<TD nowrap>&nbsp;<B>Size of mass tolerance window for XIC:</B> min=<INPUT type="text" name="mztol_min" value="5" size="1"> max=<INPUT type="text" name="mztol_max" value="5" size="1"></B><SELECT name="XIC_range"><OPTION value="ppm">ppm</OPTION><OPTION value="mz">mz</OPTION></SELECT></TD>
	</TR>
	<TR><TD nowrap>&nbsp;<B>Peak matching for XIC:</B>
	<SELECT name="XIC_rt">
	<OPTION value="post_matching">Post matching mode</OPTION>
	<OPTION value="real_or_mean">Real or mean mode</OPTION>
	<OPTION value="mean">Mean mode</OPTION>
	</SELECT>
	</TD>
	</TR>
	<TR>
	<TD nowrap>&nbsp;<B>Detection threshold between <INPUT type="text" name="dt_start" value="30000" size="4"> to <INPUT type="text" name="dt_stop" value="50000" size="4"></B></TD>
	</TR>
	<TR>
	<TD>&nbsp;&nbsp;<B><U>XIC filtering</U>:</B></TD>
	</TR>
	<TR><TD bgcolor=$lightColor nowrap><INPUT type="checkbox" name="asFilter" value="1" checked>&nbsp;<B>Anti-Spike:</B><INPUT type="text" name="anti_spike" value="5" size="1"></TD><TR>
	<TR><TD bgcolor=$lightColor nowrap><INPUT type="checkbox" name="hmFilter" value="1">&nbsp;<B>Half-Mediane:</B> min=<INPUT type="text" name="med_min" value="5" size="1"> max=<INPUT type="text" name="med_max" value="40" size="1"></TD><TR>
	<TR><TD bgcolor=$lightColor nowrap><INPUT type="checkbox" name="smFilter" value="1">&nbsp;<B>Smoothing:</B><INPUT type="text" name="smooth_val" value="3" size="1"></TD><TR>
	</TABLE>
	</DIV>
	</TD>
</TR>
</TABLE>
<BR>
|;
}

print qq
|<BR><TABLE border=0 cellspacing=0 cellpadding=0>
<TR bgcolor="$darkColor">
	<TH class="rbBorder"><INPUT type="checkbox" onclick="checkall(this.checked)"></TH>
	<TH class="rbBorder" nowrap colspan=2>&nbsp;Analysis&nbsp;</TH>
	<TH class="rbBorder">&nbsp;Peptide quantification&nbsp;</TH>
	<TH class="rbBorder">&nbsp;MS type&nbsp;<BR>& File</TH>
	<TH class="rbBorder">&nbsp;Instrument&nbsp;</TH>
	<TH class="rbBorder">&nbsp;Search file<BR>&nbsp;& Engine&nbsp;</TH>
	<TH class="rbBorder">&nbsp;Databank(s)&nbsp;<BR>&nbsp;Taxonomy&nbsp;</TH>
	<TH class="rbBorder" nowrap>&nbsp;Min. score&nbsp;<BR>&nbsp;Max. rank&nbsp;</TH>
	<TH class="bBorder">&nbsp;Validated&nbsp;<BR>&nbsp;proteins&nbsp;</TH>
</TR>
|;


my %itemIcones=&promsConfig::getItemIcones;
my $bgColor=($ITEM eq 'SAMPLE' || $ITEM eq 'SPOT')? $lightColor : $darkColor;
my %prevItemName;
my $disabSubmit=' disabled';

my %filesList;
my $mzXMLPath="$promsPath{tmp}/upload/project_$projectID";
if ($quantifType eq 'XICMCQ' && -e $mzXMLPath) {
	opendir (DIR, $mzXMLPath) || print "ERROR: Unable to read '$mzXMLPath' !<BR>\n";
	while (defined (my $currentmzXMLFile = readdir (DIR))) {
		next unless ( $currentmzXMLFile =~ /.+\.mzXML\Z/ || $currentmzXMLFile =~ /.+\.mzML\Z/ );
		$filesList{$currentmzXMLFile}=0;
	}
	closedir DIR;
}

foreach my $refAnaData (@itemAnalyses) {
	my ($anaID,$valStat,$msType,$dataFile,$fileFormat,$wiffFile,$taxonomy,$maxRank,$minScore,$instrument,$refDbUsed,$labelStrg,@projHierarchy)=@{$refAnaData};
	$taxonomy='Unknown' unless $taxonomy;
	$taxonomy=~s/\(.*\)//;
	#my $okQuantifAna=($valStat>=1 && $msType ne 'PMF' && (($quantifType eq 'EMPAI' && ($fileFormat eq 'MASCOT.DAT' || $fileFormat eq 'MASCOT.PDM')) || ($quantifType eq 'SIN' && $fileFormat ne 'PHENYX.XML') || $quantifType eq 'XIC'))? 1 : 0;
	my $okQuantifAna=0;
	if ($valStat>=1 && $msType ne 'PMF') {
		if (($quantifType eq 'EMPAI' && ($fileFormat eq 'MASCOT.DAT' || $fileFormat eq 'MASCOT.PDM'))) {$okQuantifAna=1;}
		elsif ($quantifType eq 'SIN' && $fileFormat eq 'MASCOT.DAT') {$okQuantifAna=1;}
		#elsif ($quantifType eq 'XIC') {$okQuantifAna=1;}
		elsif ($quantifType eq 'XICMCQ') {$okQuantifAna=1;}
		elsif ($quantifType eq 'XICCORR') {$okQuantifAna=1 if $anaPeptideQuantifs{$anaID};}
	}

	$msType=$msTypeName{$msType}; # global
	$fileFormat=~s/\..*//;
	$instrument='-' unless $instrument;
	$maxRank='-' unless $maxRank;
	$minScore='-' unless defined $minScore;
	$disabSubmit='' if $okQuantifAna; # at least 1 selectable analysis => activate Submit
	##>Row color
	my $fatherIt=$projHierarchy[-3];
	if ($fatherIt && (!$prevItemName{$fatherIt} || $prevItemName{$fatherIt} ne $projHierarchy[-2])) { # keep color if same analysis parent item (SAMPLE or SPOT)
		$bgColor=($bgColor eq $lightColor)? $darkColor : $lightColor;
	}
	elsif ($ITEM eq 'EXPERIMENT' || $ITEM eq 'GEL2D') {print "<TR bgcolor=$bgColor><TD colspan=2></TD><TD colspan=8><HR width=98%></TD></TR>\n";}
	print "<TR valign=middle bgcolor=$bgColor>\n";
	##>Checkbox
	if ($okQuantifAna) {
		my $disabStrg=(($quantifType eq 'XICMCQ' && $anaLabeling{$anaID} ne 'FREE') || ($quantifType eq 'XICCORR' && !$anaPeptideQuantifs{$anaID}))? ' disabled' : '';
		print "\t<TH valign=middle><INPUT type=\"checkbox\" name=\"anaList\" value=\"$anaID\" data-labeling=\"$anaLabeling{$anaID}\"$disabStrg/></TH>\n";
	}
	else {print "\t<TH valign=middle>-</TH>\n";}
	#>Parents & Analysis
	my $parentStrg='';
	for (my $i=0;$i<=$#projHierarchy-2;$i+=2) { # stops before ana name
		my $IT=$projHierarchy[$i];
		my $itName=$projHierarchy[$i+1];
		if ($prevItemName{$IT} && $prevItemName{$IT} eq $itName) {
			#$parentStrg.="<FONT style=\"visibility:hidden\">&nbsp;&nbsp;$itName&nbsp;&nbsp;&gt;&nbsp;&nbsp;</FONT>";
		}
		else {
			$parentStrg.="&nbsp;&nbsp;$itName&nbsp;&nbsp;&gt;&nbsp;&nbsp;";
			$prevItemName{$projHierarchy[$i]}=$itName;
			for (my $j=$i+2;$j<$#projHierarchy-1;$j+=2) {$prevItemName{$projHierarchy[$j]}='';}
		}
	}
	my $anaCode=($valStat==-1)? 'analysis:no_scan' : ($valStat==0)? 'analysis:no_val' : ($valStat==1)? 'analysis:part_val' : 'analysis:val';
	print qq
|	<TH nowrap align=left valign=middle>$parentStrg</TH>
	<TH nowrap align=left valign=middle><IMG src="$promsPath{images}/$itemIcones{$anaCode}">&nbsp;$projHierarchy[-1]&nbsp;</TH>
|;
	##>Labeling method
	print "<TD>";
	if (!$labelStrg) {print '&nbsp;None&nbsp;';}
	elsif ($quantifType eq 'XICCORR' && $anaPeptideQuantifs{$anaID}) { # No peptide quantification for non-validated analysis
		print "<SELECT id=\"quantif_$anaID\" name=\"quantif_$anaID\" style=\"width:250px\">\n";
		print "<OPTION value=\"\">-= Select =-</OPTION>\n" if scalar @{$anaPeptideQuantifs{$anaID}} > 1;
		foreach my $refQuantif (@{$anaPeptideQuantifs{$anaID}}) {
			my ($quantifID,$quantifName)=@{$refQuantif};
			print "<OPTION value=\"$quantifID\">$quantifName</OPTION>\n";
		}
		print "</SELECT>\n";
	}
	else {print $labelStrg;}
	print "</TD>\n";

	my $dbStrg;
	foreach my $dbID (@{$refDbUsed}) {
		$dbStrg.='&bull;' if $dbStrg;
		$dbStrg.=$listDataBank{$dbID};
	}
	print qq
|	<TD class="center">$msType<BR>&nbsp;$wiffFile&nbsp;</TD>
	<TD class="center">&nbsp;$instrument&nbsp;</TD>
	<TD class="center">&nbsp;$dataFile&nbsp;<BR>&nbsp;$fileFormat&nbsp;</TD>
	<TD class="center">&nbsp;$dbStrg&nbsp;<BR>&nbsp;$taxonomy&nbsp;</TD>
	<TD class="center">&nbsp;$minScore&nbsp;<BR>&nbsp;$maxRank&nbsp;</TD>
|;
	##>Validated proteins
	if ($valStat>=1) {print "<TD class=\"center\">$anaProteins{$anaID}[0] ($anaProteins{$anaID}[1])</TD>\n";}
	else {print "<TD class=\"center\">$anaProteins{$anaID}[0] / $anaProteins{$anaID}[1]</TD>\n";}
	print "</TR>\n";
	if ($quantifType eq 'XICMCQ' && $okQuantifAna) {
		###> Pre-selection of mzXML file
		my $listmzXMLFiles='<OPTION value=\"\">-= Select =-</OPTION>';
		my $selected;
		my $msName='';
		($msName=$wiffFile)=~s/\.[^\.]+//;
		$msName=~ s/\s+\Z//;
		foreach my $dataFile (sort {lc($a) cmp lc($b)} keys %filesList) {
			#$bgColor=($bgColor eq $darkColor)? $lightColor : $darkColor;
			$selected=($dataFile=~/$msName/)? ' selected': '';
			$listmzXMLFiles.="<OPTION value=\"$dataFile\"$selected>$dataFile</OPTION>\n";
		}
		$listmzXMLFiles.="</SELECT>";
		my $disabStrg=($anaLabeling{$anaID} eq 'FREE')? '' : ' disabled';
		#print "<TR bgcolor=$bgColor><TD colspan=2></TD><TD colspan=8><B>mzXML file:<INPUT type=\"file\" name=\"file_$anaID\" id=\"file_$anaID\" value=\"\" size=80></TD></TR>\n";
		print "<TR bgcolor=$bgColor><TD colspan=2></TD><TD colspan=8><B>mzXML file:<SELECT name=\"file_$anaID\" id=\"file_$anaID\"$disabStrg>$listmzXMLFiles</TD></TR>\n";
		if ($ITEM ne 'EXPERIMENT' && $ITEM ne 'GEL2D') {
			$bgColor=($bgColor eq $darkColor)? $lightColor : $darkColor;
		}
	}
}

#my $xicString=($quantifType eq 'XIC')? "<TR><TD colspan=10><B>Criteria for XIC extraction:</B><BR>Mass-Window=<INPUT name=\"windowMass\" value=\"0.1\" size=\"5\"> Da.<BR>Time-window=<INPUT name=\"windowTime\" value=\"90\" size=\"2\"> seconds.<BR>Peptide Rank=<INPUT name=\"pepRank\" value=\"2\" size=\"1\">.</TD></TR>":'';
#print "<TR bgcolor=$bgColor><TD colspan=10><B>Criteria used for XIC extraction:</B> <INPUT name=\"windowMass\" value=\"0.1\" size=\"5\"> Da - Time-window=<INPUT name=\"windowTime\" value=\"90\" size=\"5\"> seconds. - Peptide Rank=<INPUT name=\"pepRank\" value=\"2\" size=\"3\"></TD>";

if ($quantifType eq 'XICCORR') {
	print qq
|<TR><TD colspan=10 style="font-size:2px">&nbsp;</TD></TR>
<TR bgcolor="$darkColor"><TD colspan=10 style="font-size:2px">&nbsp;</TD></TR>
<TR bgcolor="$darkColor"><TD colspan=10 class="title3">&nbsp;Use isotopic distribution from:<SELECT name="product" onchange="updateCorrectionLabeling(this)">
|;
	if (scalar keys %isotopicDistributions > 1) {print "<OPTION value=\"\" data-labeling=\"-\">-= Select =-</OPTION>\n";}
	foreach my $label (sort {$a cmp $b} keys %isotopicDistributions) {
		print "<OPTGROUP label=\"$label:\">\n";
		foreach my $productID (sort {lc($isotopicDistributions{$label}{$a}[0]) cmp lc($isotopicDistributions{$label}{$b}[0])} keys %{$isotopicDistributions{$label}}) {
			print "<OPTION value=\"$productID\" data-labeling=\"$label\">$isotopicDistributions{$label}{$productID}[0] [$isotopicDistributions{$label}{$productID}[3]plex] (Product #$isotopicDistributions{$label}{$productID}[1], Lot #$isotopicDistributions{$label}{$productID}[2])</OPTION>\n";
		}
		print "</OPTGROUP>\n";
	}
	print qq
|</SELECT>&nbsp;</TH></TR>
<TR bgcolor="$darkColor"><TD colspan=10 style="font-size:2px">&nbsp;</TD></TR>
<TR bgcolor="$darkColor"><TD colspan=10 class="title3">&nbsp;<INPUT type="checkbox" name="keep_old" value="1" checked>Keep original peptide quantification(s)</TD></TR>
<TR bgcolor="$darkColor"><TD colspan=10 style="font-size:2px">&nbsp;</TD></TR>
|;
}

print qq
|<TR><TD colspan=10 style="font-size:2px">&nbsp;</TD>
<TR><TD colspan=10><INPUT type="submit" name="launch" value="Launch $quantifProcesses{$quantifType}" class="title3"$disabSubmit>&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="cancelAction();"></TD></TR>
</TABLE>
</FORM>
</CENTER>
</BODY>
<HTML>
|;

#########################<<< SUBROUTINE >>>######################
sub launchQuantifications {

	####<Starting HTML>####
	print header(-'content-encoding'=>'no');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>Launching Quantifications</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER><FONT class="title1">Launching Quantification(s)</FONT></CENTER>
<BR><BR><BR>
|;
	my @quantItemList=param('anaList');
	if ($quantifType eq 'XICCORR') {
		my %pepQuantifs;
		foreach my $anaID (@quantItemList) {
			my $pepQuantID=param("quantif_$anaID");
			next unless $pepQuantID;
			$pepQuantifs{$pepQuantID}=$anaID unless $pepQuantifs{$pepQuantID}; # remove duplicate pep quantif (just to be safe)
		}
		@quantItemList=();
		while (my ($pepQuantID,$anaID) = each %pepQuantifs) { # <- peptide quantif ids!!!
			push @quantItemList,$anaID.'.'.$pepQuantID;
		}
	}
	
	my $numJobs=1; # default
	if ($quantifType=~/SIN|EMPAI|XICCORR/) { # not for XICMCQ
		$numJobs=scalar @quantItemList;
	}

	mkdir "$promsPath{tmp}/quantification" unless -e "$promsPath{tmp}/quantification";
	my $currentQuantifDir="$promsPath{tmp}/quantification/current";
	mkdir $currentQuantifDir unless -e $currentQuantifDir;

	my $refConfName; # for multi-Q only

	my $quantifDate=strftime("%Y%m%d%H%M%S",localtime);
	while (-e "$promsPath{tmp}/quantification/$quantifDate") { # to prevent multi-user launch collision
		sleep 2;
		$quantifDate=strftime("%Y%m%d%H%M%S",localtime);
	}

	print '<BR><FONT class="title3">Launching quantification process';
	if ($numJobs > 1) {print "es [Master job=#$quantifDate]:\n";}
	else {print "[job #$quantifDate]...";}

	my @jobList;
	foreach my $jobPos (1..$numJobs) {

		sleep 1 if $jobPos > 1; # wait 1 sec between jobs
		my $jobDir=$quantifDate;
		$jobDir.=".$jobPos" if $numJobs > 1;

		my $quantifDir="$promsPath{tmp}/quantification/$jobDir";
		mkdir $quantifDir || die "ERROR detected: $!";

		print "<BR>&nbsp;&nbsp;-$jobPos/$numJobs" if $numJobs > 1;

		###<jobs info & data >###
		open (INFO,">$quantifDir/quantif_info.txt"); # item valueR valueDB
		print INFO "USER=$userID\n";
		if ($quantifType eq 'XICMCQ'){print INFO "TYPE=XIC\n";} # MCQ -> Alternative to XIC extraction
		else {print INFO "TYPE=$quantifType\n";}
		my $quantItemStrg='';
		if ($quantifType eq 'XICMCQ') { # MassChroQ extraction
			my $raw=param('rawdata_extraction'); # profile or centroid
			my $algo=param('alignment_method'); # OBI-Warp or in-house algo
			my $anaReference=param('refAna_alignment'); # Reference analysis
			my $mzAnaRangeMin=param('mz_start'); # range for alignment -OBIwarp
			my $mzAnaRangeMax=param('mz_stop'); # range for alignment -OBIwarp
			my $mzRangeMin=param('mzRange_start'); # range for allChargeState computation
			my $mzRangeMax=param('mzRange_stop'); # range for allChargeState computation
			my $tendency=param('ms2_tendency'); # ms2 parameter
			my $smouthMS2=param('ms2_smoothing'); # ms2 parameter
			my $smouthMS1=param('ms1_smoothing'); # ms2 parameter
			my $xicType=param('XIC_type'); # sum or max
			my $xicRange=param('XIC_range');
			my $mzTolMin=param('mztol_min');
			my $mzTolMax=param('mztol_max');
			my $xicValType=param('XIC_rt'); # real_or_mean, mean or post_matching
			my $detectionThresholdMin=param('dt_start');
			my $detectionThresholdMax=param('dt_stop');
			my $antiSpike=param('anti_spike'); # anti-spike value
			my $medMin=param('med_min'); # half-median min
			my $medMax=param('med_max'); # half-median max
			my $smoothVal=param('smooth_val'); # smoothing value
			my $allChargeStates=param('allChargeStates')?param('allChargeStates') : 0;
			my $extractTraces=param('traces')?param('traces'):0;

			print INFO "PARAMETERS:\n";
			print INFO "QUANTIF_NAME\t\t",param('quantifName'),"\n";
			print INFO "RAWDATA_ACQUISITION\t\t$raw\n";
			print INFO "EXTRACTION_ALGO\t\t$algo\n";
			print INFO "REFERENCE\t\t$anaReference\n";
			if ($algo eq 'OBI') {
				print INFO "MZ_ALIGN_RANGE_MIN\t\t$mzAnaRangeMin\n";
				print INFO "MZ_ALIGN_RANGE_MAX\t\t$mzAnaRangeMax\n";
			}
			else{
				print INFO "MS2_TENDENCY\t\t$tendency\n";
				if (param('ms2Smooth')){
					print INFO "MS2_SMOUTHING\t\t$smouthMS2\n";
				}
				if (param('ms1Smooth')){
					print INFO "MS1_SMOUTHING\t\t$smouthMS1\n";
				}
			}

			print INFO "XIC_EXTRACTION_TYPE\t\t$xicType\n";
			print INFO "XIC_RANGE\t\t$xicRange\n";
			print INFO "MZTOL_MIN\t\t$mzTolMin\n";
			print INFO "MZTOL_MAX\t\t$mzTolMax\n";
			print INFO "XIC_VAL\t\t$xicValType\n";
			print INFO "DT_START\t\t$detectionThresholdMin\n";
			print INFO "DT_STOP\t\t$detectionThresholdMax\n";
			### Only one filter applies... Modified on 05/09/2014
			if (param('asFilter')) {
				print INFO "ANTISPIKE\t\t$antiSpike\n";
			}
			if (param('smFilter')){
				print INFO "SMOOTH\t\t$smoothVal\n";
			}
			if (param('hmFilter')){
				print INFO "MED_MIN\t\t$medMin\n";
				print INFO "MED_MAX\t\t$medMax\n";
			}
			print INFO "ALLCHARGESTATES\t\t$allChargeStates\n";
			if( $allChargeStates ) {
				print INFO "MZ_RANGE_MIN\t\t$mzRangeMin\n";
				print INFO "MZ_RANGE_MAX\t\t$mzRangeMax\n";
			}
			print INFO "TRACES\t\t$extractTraces\n";

			if (param('extractL') eq 'SILAC') {
				###> Get the labeled modification so as to write QUANTIF_ANNOT according to SILAC
				my %anaLabelMods;
				my $dbh=&promsConfig::dbConnect;
				my $sthMods=$dbh->prepare("SELECT ID_MODIFICATION,PSI_MS_NAME,MONO_MASS FROM MODIFICATION WHERE VALID_STATUS=1 AND IS_LABEL=1");
				$sthMods->execute;
				while (my ($modID,$psiMsName,$monoMass)=$sthMods->fetchrow_array) {
					@{$anaLabelMods{$modID}}=($psiMsName,$monoMass);
				}
				$sthMods->finish;
				$dbh->disconnect;

				my @channelsToExtract;
				foreach my $channel (1..3) {
					my @labelsToExtract;
					foreach my $label (1..$maxLabPerChannel) {
						if (param("lab$label:name_channel$channel") && param("lab$label:name_channel$channel") ne "") {
							# exemple: Lys8#Label:13C(6)15N(2)#K#Lys+8 -> $labelModifAlias#$labelModifName#$modifRes#$searchModifName
							my $modID=param("lab$label:channel$channel");
							my ($labName,$psiMsName,$deltaMass) = ($modID > 0 )? (param("lab$label:name_channel$channel"),@{$anaLabelMods{$modID}}) : ('None','No label','0.00000');
							my $modifRes=(param("lab$label:tg_channel$channel") ne 'ANY') ? param("lab$label:tg_channel$channel") : param("lab$label:sp_channel$channel");
							push @labelsToExtract, "$labName#$psiMsName#$modifRes##$modID"; #
						}
					}

					if (param("name_channel$channel") && param("name_channel$channel") ne "") {
						my $labelString=join("@",@labelsToExtract);
						push @channelsToExtract, "$channel;".param("name_channel$channel").";$labelString" if $labelString;
					}
				}
				my $channelString=join("::",@channelsToExtract);
				print INFO "CHANNELS\t\t$channelString\n";
			}

			foreach my $anaID (@quantItemList){
				print INFO "MZXML\t\t$anaID:".param("file_$anaID")."\n";
			}
			###print INFO "ANALYSES:\n";
		}
		else {
			if ($quantifType eq 'XICCORR') {
				print INFO "ID_PRODUCT=",param('product'),"\n";
				my $keepOld=param('keep_old') || 0;
				print INFO "KEEP_OLD=$keepOld\n";
			}
			else { # SIN,EMPAI
				print INFO "PARAMETERS:\n";
				my $qRootName=($quantifType eq 'SIN')? 'PEP_SI_' : $quantifType;
				print INFO "QUANTIF_NAME\t\t$qRootName\n";
			}
			print INFO "ANALYSES:\n";
			my $quantItemID=$quantItemList[$jobPos-1];
			$quantItemStrg=$quantItemID;
			#$quantItemID.='.'.param("quantif_$quantItemID") if $quantifType=~/SILAC|ITRAQ|TMT/;
			print INFO "$quantItemID\n";
			open(FLAG,">$currentQuantifDir/$quantItemID\_$jobDir\_wait.flag"); # flag file used by watchQuantif
			print FLAG "#";
			close FLAG;
		}
		close INFO;

		#my $quantItemStrg=join(':',@quantItemList);
		#$quantItemStrg='' unless $quantItemStrg;
		if ($numJobs==1) {
			###<Forking to launch quantifications>###
			my $childPid = fork;
			unless ($childPid) { # child here
				#>Disconnecting from server
				open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
				open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
				#open STDERR, ">>$promsPath{logs}/launchQuantification.log";
				system "./launchQuantifications.pl single $ENV{REMOTE_USER} $jobDir $quantifType $quantItemStrg";
#system "./launchQuantifications.pl single $ENV{REMOTE_USER} $jobDir $fullQuantifType $quantItemStrg 2> $promsPath{tmp}/quantification/$jobDir\_errorQuantif.txt";
				exit;
			}
		}
		else {
			push @jobList,"$jobDir#$quantifType#$quantItemStrg";
		}
	} # end of jobs loop
	print '<BR><BR>' if $numJobs > 1;
	print " Done.</FONT><BR>\n";

	####>Multi-job: fork & call launchQuantification.pl in multi-job mode<####
	if ($numJobs > 1) {
		my $childPid = fork;
		unless ($childPid) { # child here
			#>Disconnecting from server
			open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
			open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
			system "./launchQuantifications.pl multi $ENV{REMOTE_USER} ".join(',',@jobList);
#system "./launchQuantifications.pl multi $ENV{REMOTE_USER} ".join(',',@jobList)." 2> $promsPath{tmp}/quantification/errorMultiQuantif.txt";
			exit;
		}
	}
	print "<BR><FONT class=\"title3\"><BR>This page will refresh itself in a few seconds.</FONT>\n";

	###>Calling watch popup window<###
	sleep 3;
	print qq
|<SCRIPT type="text/javascript">
var watchQuantifWin=window.open("$promsPath{cgi}/watchQuantifications.cgi",'WatchQuantifWindow','width=1200,height=500,scrollbars=yes,resizable=yes');
watchQuantifWin.focus();
//parent.optionFrame.watchQuantifications(); // starts watchQuantificationw.cgi
</SCRIPT>
|;
#exit; # DEBUG!!!!!
	sleep 5;
	print qq
|<SCRIPT LANGUAGE="JavaScript">
//top.promsFrame.selectedAction='summary';
parent.optionFrame.selectOption(parent.optionFrame.document.getElementById('summary')); // refresh optionFrame with summary option
</SCRIPT>
</BODY>
</HTML>
|;
	exit;
}


####> Revision history
# 3.1.0 Handles peptide quantifications selection for isobaric correction (PP 02/05/19)
# 3.0.0 Removed design-based quantifications management & compatible with Isobaric XIC correction (PP 19/02/19)
# 2.2.0 Compatible with PTM-quantif with MaxQuant probabilities for any PTM (PP 04/01/19)
# 2.1.0 Added protein-level normalization options for site quantification (PP 07/12/18)
# 2.0.1 [Fix] minor bug in %#% replacement (PP 05/12/18)
# 2.0.0 Full support for multi-job launches (PP 17/09/18)
# 1.9.9 Added more contraints on peptide/PTM quantif/PhosphoRS options based on algo selection (PP 19/07/18)
# 1.9.8 Improved auto-propagation of peptide quantification selection (PP 01/06/18)
# 1.9.7 Form & template are reset by labeling change, skip template-based biological replicates selection & template's project is displayed (PP 09/05/18)
# 1.9.6 Add a class ('template') to all parameters that will be modified by the template selection (MLP 07/05/18)
# 1.9.5 Add a function to restore default parameters (use for templates selection) (MLP 26/04/18)
# 1.9.4 [Fix] minor bug in R parameter value for MSstats when no normalization (PP 23/04/18)
# 1.9.3 More modif to select a template (MLP 18/04/18)
# 1.9.2 Modif to select a template (MLP 16/04/18)
# 1.9.1 Minor change to allow proper frame refreshing (PP 06/04/18)
# 1.9.0 Multi-quantif launch support & auto-selection of all states (PP 28/03/18)
# 1.8.7 Declare NORMALIZATION_METHOD parameter even if no normalization (PP 13/03/18)
# 1.8.6 [Fix] javascript bugs for MassChroQ XML file selection check & new checks on labeling coherence (PP 13/02/18)
# 1.8.5 Minor modif (MLP 24/01/18)
# 1.8.4 Added 'clusters' parameter for MSstats (PP 24/01/18)
# 1.8.3 Minor modifications for XIC extraction with MassChroQ (GA 14/12/17)
# 1.8.2 Modification of checked options for XIC extraction with MassChroQ (GA 13/12/17)
# 1.8.1 Better detection of bio/tech/replicates from design (PP 04/12/17)
# 1.8.0 Change in parameters to match upgraded ratio quantification R scripts (PP 27/11/17)
# 1.7.8 New option for peptide selection for SSPA, fix bug to disable Super ratio for Label-free & commented substitution 3plex/2plex (PP 20/04/17)
# 1.7.7 Compatible with TMT labeling (PP 18/01/17)
# 1.7.6 Minor modif: change "required" (SELECT "refAna_alignment") by checkForm verification (MLP 09/01/17)
# 1.7.5 Add "required" for SELECT "refAna_alignment" (MLP 09/01/17)
# 1.7.4 Minor modification of masschroq parameters and color used for labelling (GA 10/10/16)
# 1.7.3 Minor modification (GA 24/08/16)
# 1.7.2 Added Protein selection/exclusion option, MSstats global standard & "Distinct peptide sequences" option for SSPA & (PP 19/08/16)
# 1.7.1 Uses states position (not indexes) in contrasts.matrix for SWATH (PP 04/08/16)
# 1.7.0 Adding Swath (MSstats) & SSPA (PP 02/08/16)
# 1.6.6 Minor modification (GA 12/05/16)
# 1.6.5 Removed all p-value correction methods except BH & hide topN selection on Label-Free modif quantification (PP 26/04/16)
# 1.6.4 Minor bug fix in normalization method selection for label-free (PP 24/03/16)
# 1.6.3 Minor display change to better distinguish Samples (PP 08/03/16)
# 1.6.2 Added matching peptides option for label-free quantification (PP 29/02/16)
# 1.6.1 Fixed uninitialized $nbCond (PP 16/02/16)
# 1.6.0 Handles complex design,"Simple ratios"=multiSILAC algorithm for SILAC triplex (PP 29/01/16)
# 1.5.2 Modification of XIC extraction form (GA 05/09/14)
# 1.5.1 New Observation management: OBS_EXPCONDITION table & check on ref/test states compatibility for labeled quantif (PP 24/07/14)
# 1.5.0 Super Ratio option & Bug fix in condition detection is some cases (PP 21/05/14)
# 1.4.1 Modification of XIC extraction form (GA 09/04/14)
# 1.4.0 Add xic-traces option for masschroq extractions (GA 31/03/14)
# 1.3.9 Modification of MassChroQ parameters (GA 28/03/14)
# 1.3.8 Minor modification of parameters for MassChroQ extraction (GA 25/03/14)
# 1.3.7 display always header peptide quantification (SL 25/03/14)
# 1.3.6 Syntax bug fix (PP 19/03/14)
# 1.3.5 Uses mkdir instead of make_path (PP 10/03/14)
# 1.3.4 Change syntax $x=$y // $z to old ()? one (PP 10/03/14)
# 1.3.3 Peptide quantif selection & improved state selection for internal quantif (PP 26/03/14)
# 1.3.2 Various Fixes & updates <BR> Fix SILAC No Label channel annotation for MassChroQ XIC (PP 20/02/14)
# 1.3.1 Form update for label XIC extraction (GA 11/12/13)
# 1.3.0 system command removal (PP 08/11/13)
# 1.2.9 Bug fix in infinite ratios switch for DESIGN (PP 23/10/13)
# 1.2.8 Added infinite ratios switch parameter (PP 21/10/13)
# 1.2.7 iTRAQ -> ITRAQ for QUANTIF_ANNOT field (PP 07/10/13)
# 1.2.6 Minor display bug fix (PP 30/09/13)
# 1.2.5 MCQ -> XICMCQ & launch prepared like DESIGN (PP 09/09/13)
# 1.2.4 'TNPQ'->'DESIGN' & multi-labeling management (PP 29/08/13)
# 1.2.3 Minor change (PP 12/07/13)
# 1.2.2 Added pool management (PP 12/07/13)
# 1.2.1 Fix bugs popup MS2 & condition display (PP 05/07/13)
# 1.2.0 Add MS2 parameters for MassChroQ alingment (GA 04/07/13)
# 1.1.9 Better label detection & extended design ... (PP 04/07/13)
# 1.1.8 Modification for QUANTITATION launch based on new ANA_EXPCONDITION table and log file addition launchQuantification.log (GA 08/03/13)
# 1.1.7 Added more peptide filters: missed cleavage and PTM/unmodified peptide pair (PP 11/01/13)
# 1.1.6 Multi-databank search, bug fix in SILAC/iTRAQ quantif_info.txt file & code cleaning (PP 08/01/13)
# 1.1.5 Minor modification : remove 3 replicates option for Pep-Ratio + add checkbox for MCQ extraction (GA 05/11/12)
# 1.1.4 Modification for XIC extraction (quantif-name) and XIC parameters + Modification in TNPQ parameters (GA 07/11/12)
# 1.1.3 Minor modification ACT=experiment to make it homogene with openProject.cgi (GA 18/10/12)
# 1.1.2 Add MassChroQ extraction option -> similar to XIC (GA 04/09/12)
# 1.1.1 Check on SILAC/iTRAQ peptide quantif status before allowing protein quantif (PP 07/06/12)
# 1.1.0 Minor changes & bug fixes (PP 30/04/12)
# 1.0.9 Minor bug correction -> remove printError call (GA 18/04/12)
# 1.0.8 minor bug correction for perl/javascript: 'eq' vs '==' (GA 04/04/12)
# 1.0.7 merge 1.0.6GA & 1.0.5PP (PP 02/04/12)
# 1.0.6 Updates for label-free based quantification like TNPQ, XIC, etc (GA 16/01/12)
# 1.0.5 Updates for label-based internal quantifications (PP 26/08/11)
# 1.0.4 More Quantif methods (GA ...)
# 1.0.3 Fix file detection bug for T3PQ
