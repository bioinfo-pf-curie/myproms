#!/usr/local/bin/perl -w

################################################################################
# startGOAnalysis.cgi       1.1.5                                              #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Fetch and provide proteins and parameters for GO enrichment analysis         #
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
use strict;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use POSIX qw(strftime); # to get the time
use File::Copy;
use File::Copy::Recursive qw(dirmove);
use File::Path qw(rmtree);
use goAnalysis;
use promsOboParser;
use GO::AnnotationProvider::AnnotationParser;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
$CGITempFile::TMPDIRECTORY=$promsPath{'tmp'};

my $dbh=&promsConfig::dbConnect;
my $expID = param('id_exp');
my $projectID=&promsMod::getProjectID($dbh,$expID,'EXPERIMENT');

if(param('start')){
    &processForm;
    exit;
}
#my $sthGetProjectID = $dbh->prepare("SELECT ID_PROJECT FROM EXPERIMENT WHERE ID_EXPERIMENT=?");
#$sthGetProjectID->execute($expID);
#my ($projectID) = $sthGetProjectID->fetchrow_array;
#$sthGetProjectID->finish;

#############################
# Generating selection tree #
#############################
my %analysisClass=(-1=>'no_scan',0=>'no_val',1=>'part_val',2=>'val');
my %treeOptions;
$treeOptions{'BUTTON'}='font-size:9px;width:55px';
%{$treeOptions{'CHECKBOX'}}=('experiment'=>0,'gel2d'=>1,'spot'=>1,'sample'=>1,'analysis'=>1); # all items are potentially checkable
#my %giItem;

my $expName = $dbh->selectrow_array("SELECT NAME FROM EXPERIMENT WHERE ID_EXPERIMENT=$expID");
my $sthG2D=$dbh->prepare("SELECT ID_GEL2D,NAME FROM GEL2D WHERE ID_EXPERIMENT=? ORDER BY DISPLAY_POS ASC");
my $sthSpot=$dbh->prepare("SELECT SPOT.ID_SPOT,SPOT.NAME,ID_SAMPLE,SAMPLE.NAME FROM SPOT,SAMPLE WHERE SPOT.ID_SPOT=SAMPLE.ID_SPOT AND ID_GEL2D=? ORDER BY SPOT.NAME ASC");
my $sthSamp=$dbh->prepare("SELECT ID_SAMPLE,NAME FROM SAMPLE WHERE ID_EXPERIMENT=? AND ID_SPOT IS NULL ORDER BY DISPLAY_POS ASC");
my $sthAna=$dbh->prepare("SELECT ID_ANALYSIS,ANALYSIS.NAME,VALID_STATUS FROM ANALYSIS WHERE ID_SAMPLE=? ORDER BY DISPLAY_POS ASC");

my @experimentTree=(0,'experiment',$expID,'','',1,$expName,'');
## The experiment ##
$sthG2D->execute($expID);
while(my ($gelID,$gelName) = $sthG2D->fetchrow_array){
    my @gelTree=(1,'gel2d',$gelID,'','',0,$gelName,'');
    my $chkBoxGel=0;

    ## Gels ##
    $sthSpot->execute($gelID);
    while(my ($spotID,$spotName,$sampID,$sampName) = $sthSpot->fetchrow_array){
	my @spotTree=(2,'spot',$spotID,'','',0,$spotName,'');
	my $chkBoxSpot=0;

	## Samples ##
	if ($sampID) {
	    $sthAna->execute($sampID);
	    while (my ($anaID,$anaName,$validStatus)=$sthAna->fetchrow_array){
		push @spotTree,[3,'analysis',$anaID,'',$analysisClass{$validStatus},0,$anaName,''];
		if($validStatus >= 1){
		    $chkBoxSpot++
		} else {
		    $treeOptions{'DISCHECK'}{"analysis:$anaID"}=1;
		}
		#if($identifierType eq 'GI_ACCESSION'){
		#    $giItem{analysis}{$anaID} = 1;
		#    $giItem{sample}{$sampID} = 1;
		#    $giItem{gel}{$gelID} = 1;
		#}
	    }
	    push @gelTree,\@spotTree;
	    if ($chkBoxSpot > 0) {$chkBoxGel++;}
	    else {$treeOptions{'DISCHECK'}{"spot:$spotID"}=1;}
	}
    }
    push @experimentTree,\@gelTree;
    unless ($chkBoxGel > 0) {$treeOptions{'DISCHECK'}{"gel2d:$gelID"}=1;}
}
## Free samples ##
$sthSamp->execute($expID);
while (my ($sampID,$sampName)=$sthSamp->fetchrow_array) {
    my @sampleTree=(1,'sample',$sampID,'','',0,$sampName,'');
    my $chkBoxSamp=0;

    #>Analyses<#
    $sthAna->execute($sampID);
    while (my ($anaID,$anaName,$validStatus)=$sthAna->fetchrow_array) {
	push @sampleTree,[2,'analysis',$anaID,'',$analysisClass{$validStatus},0,$anaName,''];
	if($validStatus>=1){
	    $chkBoxSamp++;
	} else {
	    $treeOptions{'DISCHECK'}{"analysis:$anaID"}=1;
	}
	#if($identifierType eq 'GI_ACCESSION'){
	#	    $giItem{analysis}{$anaID} = 1;
	#	    $giItem{sample}{$sampID} = 1;
	#}
    }
    push @experimentTree,\@sampleTree;
    $treeOptions{'DISCHECK'}{"sample:$sampID"}=1 unless $chkBoxSamp>0;
}

$sthG2D->finish;
$sthSpot->finish;
$sthSamp->finish;
$sthAna->finish;

# Fetching ontology files #
my $disableSubmitStrg = '';
my $ontologyString;
my $sthOBO = $dbh->prepare("SELECT ID_ONTOLOGY, OBO_FILE, NAME FROM ONTOLOGY WHERE STATUS>=1");
$sthOBO->execute;
my $refOboList = $sthOBO->fetchall_arrayref;
if(scalar @{$refOboList}){
    $ontologyString =  "<SELECT name=\"ontology\">\n";
    foreach my $oboEntry (@{$refOboList}){
	my ($oboID,$oboFileName,$name) = @{$oboEntry};
	$ontologyString .= "<OPTION value=\"$oboID\">$name ($oboFileName)</OPTION>\n";
    }
    $ontologyString .= "</SELECT>";
} else {
    $ontologyString = "&nbsp;<FONT color=\"red\">No ontology file found in database</FONT>\n";
    $disableSubmitStrg = 'disabled';
}

# Fetching annotation files #
my $annotationString;
my $sthGOA = $dbh->prepare("SELECT ID_GOANNOTATION, NAME, DES, ID_SPECIES, ANNOT_FILE FROM GOANNOTATION WHERE STATUS=1");
$sthGOA->execute;
my $refGOAList = $sthGOA->fetchall_arrayref;
if(scalar @{$refGOAList}){
    $annotationString = "<SELECT name=\"annotation\">\n";
    foreach my $goaEntry (@{$refGOAList}){
		my ($id,$name,$description,$id_species,$annotFile) = @{$goaEntry};
		my ($speciesSciName,$speciesComName) = $dbh->selectrow_array("SELECT SCIENTIFIC_NAME, COMMON_NAME FROM SPECIES WHERE ID_SPECIES=$id_species");
		$annotationString .= "<OPTION value=\"$id\" onmouseover=\"popup('$description')\" onmouseout=\"popout()\">$name - $speciesSciName";
		if($speciesComName){
			$annotationString .= " ($speciesComName)";
		}
		$annotationString .= "</OPTION>\n";
    }
    $annotationString .= "</SELECT>";
}
else {
    $annotationString = "&nbsp;<FONT color=\"#DD0000\">No annotation file found in database</FONT>\n";
    $disableSubmitStrg = 'disabled';
}

###########################################################
# Fetching categories for background population selection #
###########################################################
my $sthClass = $dbh->prepare("SELECT ID_CLASSIFICATION,NAME FROM CLASSIFICATION WHERE ID_PROJECT=?");
my $sthCategory = $dbh->prepare("SELECT ID_CATEGORY,NAME FROM CATEGORY WHERE ID_CLASSIFICATION=?");
#my $populationStrg = "<SELECT name=\"population\" onchange=\"selectPopulation(this.selectedIndex)\" style=\"width:200px\"><OPTION value=\"0\">-= Select =-</OPTION>\n";
my $categoryStrg="<OPTION value=\"0\">-= Select =-</OPTION>\n";
$sthClass->execute($projectID);
while (my ($classID,$className) = $sthClass->fetchrow_array){
    $categoryStrg .= "<OPTGROUP label=\"Theme: $className\">\n";
    $sthCategory->execute($classID);
    while (my ($categoryID,$categoryName) = $sthCategory->fetchrow_array){
		$categoryStrg .= "<OPTION value=$categoryID>$categoryName</OPTION>\n";
    }
}
#$populationStrg .= "</SELECT>";
$sthClass->finish;
$sthCategory->finish;
$dbh->disconnect;
#my $protSetClassStrg = $populationStrg;
#$protSetClassStrg =~ s/name=\"population\"/name=\"protSet\"/;

###########################
###### Starting HTML ######
###########################
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
print header(-'content-encoding'=>'no',-charset=>'UTF-8');
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
	TD {font-size:13px;font-weight:bold;}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
|;
	&promsMod::popupInfo();
	&promsMod::writeJsTreeFunctions(\@experimentTree,\%treeOptions);

#print "var giItem = {};\n";
#foreach my $itemType (keys %giItem){
#    print " giItem['$itemType'] = {};\n";
#    foreach my $itemID (keys %{$giItem{$itemType}}){
#	print "  giItem['$itemType']['$itemID'] = 1;\n";
#    }
#}

print qq
|function checkForm(myForm){
    var criteriaValue = getRadioValue(myForm.criteria);
    var protSetRadioValue = getRadioValue(myForm.chProtSet);
    myForm.checkedItems.value=getCheckedItems();
    if(!myForm.anaName.value){
		alert('Type a name for this GO analysis.');
		return false;
    }
	else if (!myForm.aspect[0].checked && !myForm.aspect[1].checked && !myForm.aspect[2].checked){ // P,C,F
		alert('Select at least 1 domain for the analysis.');
		return false;
    }
	else if (criteriaValue == 'FDR' && (myForm.FDR.value > 100 \|\| myForm.FDR.value <= 0)){
		alert('Enter a valid FDR.')
		return false;
    }
	else if (criteriaValue == 'pval' && (myForm.pval.value > 1 \|\| myForm.pval.value <= 0)){
		alert('Enter a valid p-value threshold.')
		return false;
    }
	else if(!myForm.checkedItems.value && protSetRadioValue == 'analysis'){
		alert('Select at least 1 item for the analysis.');
		return false;
    }
	else if(myForm.protSet[myForm.protSet.selectedIndex].value==0 && protSetRadioValue == 'category') {
		alert('Select a category.');
		return false;
    }
/*
	else {
		if(giAreSelected(myForm.checkedItems.value)){
			return confirm('You have selected 1 or more analyses using GI identifiers. This will greatly increase the amount of unannotated proteins. UNIPROT identifiers are recommended. Continue ?');
		} else {
			return false;
		}
    }
*/
}
/*
function giAreSelected(chckItemString){
    var items = chckItemString.split("\+");
    for(var i=0;i<items.length;i++){
		var itemInfos = items[i].split(":");
		var itemType = itemInfos[0];
		var itemIds = itemInfos[1].split(',');
		for(var j=0;i<itemIds.length;j++){
			if(giItem[itemType][itemIds[j]]){
			return true;
			}
		}
    }
    return false;
}
*/
function selectPopulation(index){
    var myForm = document.enrichForm;
    var popValue = myForm.population.options[index].value;
    if (popValue == 0 && myForm.protNb.disabled == true){
		myForm.protNb.disabled = false;
    }
	else if (popValue != 0 && myForm.protNb.disabled == false){
		myForm.protNb.disabled = true;
    }
}
function getRadioValue(radio){
    var value;
    var len = radio.length;
    for(i = 0;i<len;i++){
		if(radio[i].checked){
			value = radio[i].value;
		}
    }
    return value;
}
function selectStatistics(){
    var myForm = document.enrichForm;
    var value = getRadioValue(myForm.criteria);
    if(value == 'FDR'){
		myForm.pval.disabled = true;
		myForm.bonferroni.disabled = true;
		myForm.FDR.disabled = false;
		myForm.FDRmethod.disabled = false;
    }
	else if (value == 'pval'){
		myForm.FDR.disabled = true;
		myForm.FDRmethod.disabled = true;
		myForm.pval.disabled = false;
		myForm.bonferroni.disabled = false;
    }
}
function checkMinPep(){
    var myForm = document.enrichForm;
    if(myForm.chMinPep.checked == true){
		myForm.minPep.disabled=false;
    }
	else {
		myForm.minPep.disabled=true;
    }
}
function actionOnSelect(){}
function actionOnCheck(){}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER><FONT class="title">Gene Ontology Enrichment Analysis</FONT><BR><BR><BR>

<TABLE border=0 cellpadding=0 cellspacing=0>
    <TR>
	<TD nowrap valign='top'>
	    <FORM name="enrichForm" method="post" onsubmit="return(checkForm(this))" enctype="multipart/form-data">
	    <INPUT type='hidden' name='id_exp' value="$expID" />
	    <INPUT type='hidden' name='checkedItems' value=""/>
	    <TABLE bgcolor=$darkColor border=0>
		<TR>
		    <TH valign='top' align='right'>&nbsp;Name :</TH><TD bgcolor=$lightColor><INPUT type='text' name='anaName' size='60'></TD>
		</TR>
		<TR>
		    <TH valign='top' align='right'>&nbsp;Description :</TH><TD bgcolor=$lightColor><TEXTAREA rows="2" cols="60" name="description"></TEXTAREA></TD>
		</TR>
		<TR>
		    <TH valign='top' align='right'>&nbsp;Ontology File :</TH><TD bgcolor=$lightColor>$ontologyString</TD>
		</TR>
		<TR>
		    <TH valign='top' align='right'>&nbsp;Annotation :</TH>
		    <TD  bgcolor=$lightColor>
			    $annotationString
		    </TD>
		</TR>
		<TR>
		    <TH valign='top' align='right'>&nbsp;Domain(s) :</TH>
		    <TD  bgcolor=$lightColor>
		    <INPUT type='checkbox' name='aspect' value='P' checked> Biological Process&nbsp;&nbsp;&nbsp;
		    <INPUT type='checkbox' name='aspect' value='C' checked> Cellular Component&nbsp;&nbsp;&nbsp;
		    <INPUT type='checkbox' name='aspect' value='F' checked> Molecular Function&nbsp;&nbsp;&nbsp;
		    </TD>
		</TR>
		<TR>
		    <TH valign='top' align='right'>&nbsp;&nbsp;Advanced Parameters :</TH><TD bgcolor=$lightColor>
			&nbsp;Estimated number of proteins in organism: <INPUT id='protNb' type='text' name='protNb' size=6> <I>(default is number of annotated proteins in organism)</I><BR>
			<TABLE border=0 cellspacing=0><TR><TD>&nbsp;Background population:</TD><TD><INPUT type="radio" name="chPopSource" value="category" checked> Select a List: <SELECT name="population" onchange="selectPopulation(this.selectedIndex)">$categoryStrg</SELECT></TD></TR>
			<TR><TD></TD><TD><INPUT type="radio" name="chPopSource" value="localFile"> Upload a local file: <INPUT type="file" name="localPopFile"></TD></TR></TABLE>
			&nbsp;Statistical settings:<BR>
			&nbsp;&nbsp;&nbsp;<INPUT type='radio' name='criteria' value='FDR' onClick="selectStatistics();" checked> Control FDR at <INPUT type='text' name='FDR' value='1' size=2>%
			with <SELECT name='FDRmethod'>
				    <OPTION value='BH' selected>Benjamini & Hochberg</OPTION>
				    <OPTION value='FDRsimulation'>Simulation</OPTION>
			</SELECT> method<BR>
			&nbsp;&nbsp;&nbsp;<INPUT type='radio' name='criteria' value='pval' onClick="selectStatistics();"> Use a p-value threshold: <INPUT type='text' name='pval' size=3 value='0.01' disabled>
			with <INPUT type='checkbox' name='bonferroni' value=1 checked disabled> Bonferroni correction<BR>
			&nbsp;<INPUT type='checkbox' name='nonSignificant'> Show non-significant terms in graphical view<BR>
			&nbsp;<INPUT type='checkbox' name='chMinPep' onclick="checkMinPep();"> Include only proteins containing at least <INPUT type='text' name='minPep' value='1' size=2 disabled> peptide(s)
		</TD></TR>
		<TR>
		    <TD align='center' colspan=2>
			<INPUT type='submit' name='start' value='Start Analysis' $disableSubmitStrg/>&nbsp;&nbsp;&nbsp;<INPUT type='reset' value='Clear' />
		    </TD>
		</TR>
	    </TABLE>
	</TD>
	<TD  valign='top'>
	    <TABLE width=100% border=0 cellspacing=0 bgcolor=$darkColor>
		<TR><TH nowrap>&nbsp;Select a protein set from:&nbsp;</TH></TR>
		<TR><TH align="left"><INPUT type="radio" name="chProtSet" value="category"> a List&nbsp;</TH></TR>
		<TR><TD bgcolor=$darkColor><SELECT name="protSet" style="width:200px">$categoryStrg</SELECT></TD></TR>
		<TR><TH bgcolor=$darkColor align="left"><INPUT type="radio" name="chProtSet" value="analysis" checked> Project items&nbsp;</TH></TR>
	    </TABLE>
	    </FORM>
	    	<DIV id="divDescription" class="clDescriptionCont">
	<!--Empty div-->
	</DIV>
	<SCRIPT type="text/javascript">setPopup();</SCRIPT>
|;
my $tableIndex=0;

&promsMod::printTree(\@experimentTree,\%treeOptions,\$tableIndex,1);

print qq|
</TABLE>
</CENTER>
</BODY>
</HTML>
|;
exit;

sub processForm{
    #my $dbh=&promsConfig::dbConnect;

    #######################
    # Fetching parameters #
    #######################
    #my %param; # contains arguments for termFinder.pl
    my $userID = $ENV{'REMOTE_USER'};
    #my $expID = param('expID');
    my $name = param('anaName');
    my $description = param('description');

    my $protNb = (param('protNb'))? param('protNb'): 0;

    my @warnings; # will be written in errorLog

    ##############
    # Start HTML #
    ##############
    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    print qq
    |<HTML>
    <HEAD>
    <LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
    </HEAD>
    <BODY background="$promsPath{images}/bgProMS.gif">
    <CENTER><BR><IMG src="$promsPath{images}/engrenage.gif"><BR><BR></CENTER>
    |;
    print "<FONT class=\"title3\">\n";
    print "Enrichment Analysis in Progress...<BR>\n";

    #-------------------#
    # Files and aspects #
    #-------------------#
    my $goaID = param('annotation');
    my $oboID = param('ontology');
    my $goaFile = "$promsPath{goa}/$goaID.goa";
    my $oboFile = "$promsPath{obo}/$oboID.obo";
    my @aspects;
    my %ontologies;

    print "Reading ontology...";
    foreach my $aspect (param('aspect')){
		push @aspects, $aspect;
		$ontologies{$aspect} = new promsOboParser(
												ontologyFile => $oboFile,
												aspect => $aspect
												);
    }
    print " OK<BR>\n";
    #$param{aspects} = \@aspects;
    #$param{ontologies} = \%ontologies;

    print "Reading annotation...";
    my $annotation = GO::AnnotationProvider::AnnotationParser->new(annotationFile=> $goaFile);
    print " OK<BR>\n";
    my ($identifierID,$identifierName) = $dbh->selectrow_array("SELECT IDENTIFIER.ID_IDENTIFIER,IDENTIFIER.NAME
							       FROM GOANNOTATION,IDENTIFIER
							       WHERE ID_GOANNOTATION=$goaID
							       AND IDENTIFIER.ID_IDENTIFIER=GOANNOTATION.ID_IDENTIFIER");
    $identifierID = 0 unless $identifierID;

    #---------------------------#
    # Fetching dataset proteins #
    #---------------------------#
	my $sthListProt = $dbh->prepare("SELECT DISTINCT P.ID_PROTEIN,P.IDENTIFIER
									FROM CATEGORY_PROTEIN CP,PROTEIN P,ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S
									WHERE CP.ID_PROTEIN=P.ID_PROTEIN
									AND P.ID_PROTEIN=AP.ID_PROTEIN
									AND AP.ID_ANALYSIS=A.ID_ANALYSIS
									AND A.ID_SAMPLE = S.ID_SAMPLE
									AND CP.ID_CATEGORY=?
									AND S.ID_EXPERIMENT=$expID"); # $expID to restrict list content to protein actually in Experiment
									#TODO: consider GEL path too

    # Min peptide value for each protein #
    my $minPep = (param('minPep'))? param('minPep'):1;

    my %protIDs; # ${uniprot} = DB_ID
    my @anaList;
    if(param('chProtSet') eq 'category'){
		print "Fetching protein set from Custom List...";
        my $catID = param('protSet');
        #my $sthProt = $dbh->prepare("SELECT DISTINCT P.ID_PROTEIN,P.IDENTIFIER
        #                        FROM CATEGORY_PROTEIN CP ,PROTEIN P,ANALYSIS A,SAMPLE S
        #                        WHERE CP.ID_PROTEIN=P.ID_PROTEIN
        #                        AND CP.ID_CATEGORY=$catID
        #                        AND A.ID_SAMPLE = S.ID_SAMPLE
        #                        AND S.ID_EXPERIMENT=$expID"); #TODO: consider GEL path too
        $sthListProt->execute($catID);
		while(my ($protID, $protName) = $sthListProt->fetchrow_array){
			if($identifierID){
			my @identifiers = goAnalysis::getProteinIds($dbh, $protID, $identifierID);
				if(scalar @identifiers){
					my $selectedID = goAnalysis::selectID(\@identifiers, $annotation);
					$protIDs{$selectedID} = $protID;
				} else {
					$protIDs{$protName} = $protID;
					push @warnings, "No $identifierName found for protein $protName, this identifier was used instead";
				}
			}
			else {
				$protIDs{$protName} = $protID;
			}
		}
        #$sthProt->finish;
		print " OK<BR>\n";
    }
    elsif (param('chProtSet') eq 'analysis'){
		print "Fetching protein set from Analyses...";
		# Fetching analyses IDs from checkedItems ids #
		my $totItemStrg = param('checkedItems');
		my %sthAnaList=(
			'gel2d'=>$dbh->prepare('SELECT ID_ANALYSIS FROM ANALYSIS,SAMPLE,SPOT
					   WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE
					   AND SAMPLE.ID_SPOT=SPOT.ID_SPOT
					   AND ID_GEL2D=?'),
			'spot'=>$dbh->prepare('SELECT ID_ANALYSIS FROM ANALYSIS,SAMPLE
					  WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE
					  AND ID_SPOT=?'),
			'sample'=>$dbh->prepare('SELECT ID_ANALYSIS FROM ANALYSIS
						WHERE ID_SAMPLE=?'),
			'analysis'=>$dbh->prepare('SELECT ID_ANALYSIS FROM ANALYSIS
						  WHERE ID_ANALYSIS=?') # after submit only
		);
		foreach my $itemInfo (split(/\+/,$totItemStrg)){ # ex: experiment:1,2,3+sample:1,3+analysis:2
			my ($chkItem,@listIDs)=split(/[:,]/,$itemInfo);
			foreach my $itID (@listIDs) {
				$sthAnaList{$chkItem}->execute($itID);
				while (my ($analysisID)=$sthAnaList{$chkItem}->fetchrow_array){
					push @anaList, $analysisID;
				}
			}
		}
		foreach my $sth (values %sthAnaList) {$sth->finish;}

		# Fetching proteins #
        my $sthProt = $dbh->prepare("SELECT PROTEIN.ID_PROTEIN,IDENTIFIER FROM PROTEIN,ANALYSIS_PROTEIN
                                WHERE PROTEIN.ID_PROTEIN=ANALYSIS_PROTEIN.ID_PROTEIN
                                AND ANALYSIS_PROTEIN.ID_ANALYSIS=?
                                AND VISIBILITY>=1
                                AND NUM_PEP>=?");
        foreach my $analysisID (@anaList){
            $sthProt->execute($analysisID, $minPep);
			while(my ($protID, $protName) = $sthProt->fetchrow_array){
				if($identifierID){
					my @identifiers = goAnalysis::getProteinIds($dbh, $protID, $identifierID);
					if(scalar @identifiers){
						my $selectedID = goAnalysis::selectID(\@identifiers, $annotation);
						$protIDs{$selectedID} = $protID;
					}
					else {
						$protIDs{$protName} = $protID;
						push @warnings, "No $identifierName found for protein $protName, this identifier was used instead";
					}
					#warn "$protID - $protName => $selectedID"; # debug
				}
				else {
					$protIDs{$protName} = $protID;
				}
			}
        }
        $sthProt->finish;
		print "OK<BR>\n";
    }

    #------------#
    # Population #
    #------------#
    my @population = ();
    if(param('chPopSource') eq 'category'){
		my $categoryID = param('population');
		print "Fetching background proteins from Custom List...";
		#my $sthProt = $dbh->prepare("SELECT CATEGORY_PROTEIN.ID_PROTEIN,IDENTIFIER FROM CATEGORY_PROTEIN,PROTEIN WHERE ID_CATEGORY=$categoryID
		#				AND CATEGORY_PROTEIN.ID_PROTEIN=PROTEIN.ID_PROTEIN");
		$sthListProt->execute($categoryID);
		while(my ($protID,$protName) = $sthListProt->fetchrow_array){
			if($identifierID){
				my @identifiers = goAnalysis::getProteinIds($dbh, $protID, $identifierID);
				if(scalar @identifiers){
					my $selectedID = goAnalysis::selectID(\@identifiers, $annotation);
					push @population, $selectedID;
				}
				else {
					push @population, $protName;
					push @warnings, "No $identifierName found for protein $protName, this identifier was used instead";
				}
				#warn "$protID - $protName => $selectedID"; # debug
			}
			else {
				push @population, $protName;
			}
		}
		print " OK<BR>\n";
    }
    elsif (param('chPopSource') eq 'localFile'){
		#$param{population} = "file:".tmpFileName(upload('localPopFile'));
		my $fileName = tmpFileName(upload('localPopFile'));
        print "# Fetching background proteins...\n";
        open(POPFILE,$fileName) or die $!;
        while(<POPFILE>){
            chomp;
            push @population, $_;
        }
        close(POPFILE);
        unlink($fileName);
        print "# OK<BR>\n";
    }
	$sthListProt->finish;


    #$param{population} = \@population if scalar @population;
    my $criteria = param('criteria');
    my $threshold = param($criteria);
    my $method;
    if ($criteria eq 'FDR'){
		$method = param('FDRmethod');
    }
	elsif ($criteria eq 'pval'){
		$method = (param('bonferroni') && param('bonferroni') == 1)? 'bonferroni': 'none';
    }

    my $drawGraph = (param('nonSignificant'))? 2 : 1;

    # Generating temp directory name for enrichment data files #
    my $fileName = strftime("%Y%m%d%H%M%S",localtime);

    # Checking if slim or complete obo #
    my ($status) = $dbh->selectrow_array("SELECT STATUS FROM ONTOLOGY WHERE ID_ONTOLOGY=$oboID");
    my $slim;
    if($status == 2){ # complete
		$slim = 0;
    }
	else { # slim
		#($slim) = $dbh->selectrow_array("SELECT ID_ONTOLOGY FROM ONTOLOGY WHERE STATUS=2 ORDER BY VERSION_DATE DESC LIMIT 0,1");
		$slim = 1;
    }

    ###########################
    ### Starting termFinder ###
    ###########################
    my $success = goAnalysis::termFinder(
			dbh => $dbh,
			ontologies => \%ontologies,
			annotation => $annotation,
			aspects => \@aspects,
			data => \%protIDs,
			protNb => $protNb,
			minPep => $minPep,
			criteria => $criteria,
			threshold => $threshold,
			method => $method,
			drawGraph => $drawGraph,
			fileName => $fileName,
			slim => $slim,
			population => \@population
			);


    ###################################################################
    ### Adding errors from this script to termFinder error log file ###
    ###################################################################
    if(scalar @warnings){
		my $errorLogFile = "$promsPath{'tmp'}/GO/$fileName/errorLog.txt";
		open ERR, ">>$errorLogFile";
		foreach (@warnings){
			print ERR "$_\n";
		}
		close ERR;
    }

    ##########################
    ### Saving in database ###
    ##########################
    if ($success) {
        my $paramStrg = "criteria=$criteria;threshold=$threshold;method=$method;drawGraph=$drawGraph;";
		if(param('chProtSet') eq 'category'){
			$paramStrg .= "protSetCategory=".param('protSet').';';
		}
		# Population #
		if(param('chPopSource') eq 'category'){
			$paramStrg .= "popCategory=".param('population').';';
		}
		elsif (param('chPopSource') eq 'localFile'){
			$paramStrg .= "popFile=".param('localPopFile').';';
		}
		# Min pep #
		if($minPep > 1){
			$paramStrg .= "minPep=$minPep;";
		}

		# Checking results for each aspect
		my $aspectStrg = join '',@aspects;
		foreach my $aspect (@aspects){
			unless(-e "$promsPath{'tmp'}/GO/$fileName/results_$aspect.txt"){
				$aspectStrg =~ s/$aspect//;
			}
		}

		my $sthGOASave = $dbh->prepare("INSERT INTO GO_ANALYSIS (ID_EXPERIMENT,ID_ONTOLOGY,ID_GOANNOTATION,NAME,DES,GOA_TYPE,ASPECT,PARAM_STRG,UPDATE_DATE,UPDATE_USER) VALUES (?,?,?,?,?,?,?,?,NOW(),?);");
		$sthGOASave->execute($expID,$oboID,$goaID,$name,$description,'graph',$aspectStrg,$paramStrg,$userID);
		$sthGOASave->finish;
		my $goAnaID = $dbh->last_insert_id(undef,undef,'GO_ANALYSIS','ID_GOANALYSIS');
		if(scalar(@anaList)){
			my $sthAnaSave=$dbh->prepare("INSERT INTO GOANA_ANALYSIS (ID_GOANALYSIS,ID_ANALYSIS) VALUES (?,?);");
			foreach my $analysisID (@anaList){
				$sthAnaSave->execute($goAnaID,$analysisID);
			}
			$sthAnaSave->finish;
		}
		my $errorText = '';
		# Check error log file #
		open(ERRLOG,"$promsPath{tmp}/GO/$fileName/errorLog.txt") or warn "Unable to read error file: $!";
		while(<ERRLOG>){
			chomp;
			$errorText .= "$_<BR>\n";
		}
		close(ERRLOG);
		# Moving temp files #
		unlink "$promsPath{tmp}/GO/$fileName/job.txt";
		unlink "$promsPath{tmp}/GO/$fileName/stdout.txt";
		mkdir $promsPath{'go_unix'} unless -d $promsPath{'go_unix'};
		mkdir "$promsPath{go_unix}/project_$projectID" unless -d "$promsPath{go_unix}/project_$projectID";
		#mkdir "$promsPath{go_unix}/project_$projectID/$goAnaID";
		#foreach(glob("$promsPath{tmp}/GO/$fileName/*")){
		#	move($_,"$promsPath{go_unix}/project_$projectID/$goAnaID/") || print STDERR "$_: $!";
		#}
		#rmtree("$promsPath{tmp}/GO/$fileName") || print STDERR $!;
		dirmove("$promsPath{tmp}/GO/$fileName","$promsPath{go_unix}/project_$projectID/$goAnaID");

		$dbh->commit;
		$dbh->disconnect;
		print "<BR>\nProcessing finished.</FONT>\n";
		sleep(3);
		print qq
|<SCRIPT language="JavaScript">
	parent.itemFrame.location='$promsPath{cgi}/openProject.cgi?ACT=experiment&VIEW=go&branchID=go_analysis:$goAnaID';
</SCRIPT>
</BODY>
</HTML>
|;
    }
	else { # failed process case
        $dbh->disconnect;
		#my $errorText = '';
		## Check error log file #
		#open(ERRLOG,"$promsPath{tmp}/GO/$fileName/errorLog.txt") or die "Unable to read error file: $!";
		#while(<ERRLOG>){
		#    chomp;
		#    $errorText .= "$_<BR>\n";
		#}
		#close(ERRLOG);
		# removing all temp files
		rmtree("$promsPath{'tmp'}/GO/$fileName") or warn "Cannot delete directory: $!";
        print qq
|<BR>
<FONT color="red">Processing failed !</FONT><BR>
</FONT>
<BR>
<CENTER>
<INPUT type="button" value="Go back" onclick="parent.resultFrame.location='$promsPath{cgi}/startGOAnalysis.cgi?id_exp=$expID';">
</CENTER>
|;
    }
}

####>Revision history<####
# 1.1.5 Now records user's login instead of user's name (PP 26/07/18)
# 1.1.4 SQL query optimization for protein in Lists (PP 20/08/15)
# 1.1.3 remove $user to tmp path (SL 17/12/14)
# 1.1.2 Changed scrollingbar image to engrenage (PP 07/11/14)
# 1.1.1 Use of use File::Copy::Recursive (PP 12/11/13)
# 1.1.0 Major bug fix in SQL request for category proteins (FY 09/10/13)
# 1.0.9 Minor bug and display correction (PP 28/05/13)
# 1.0.8 Fitting form processing with termFinder.pl conversion to goAnalysis.pm (FY 19/11/12)
# 1.0.7 Put arguments in the right order in getProjectID calling (FY 22/10/12)
# 1.0.6 Compatible with new data path .../results/project_ID/GoID (PP 13/09/12)
# 1.0.5 Add minimal peptide number criterion<BR>& Add alert message before performing test on GI_ACCESSION proteins<BR>& Change format of arguments to termFinder.pl<BR>& errorLog.txt is no longer displayed after a failed process since fatal errors are printed in JOB file (FY 19/07/12)
# 1.0.4 Send experience ID to termFinder.pl if dataset type is category<BR>& Printing a empty comment during process for each "NULL" key word in JOB file to restart server time out (FY 15/06/12)
# 1.0.3 Change in itemFrame refresh parameters (PP 16/04/12)
# 1.0.2 Fix uninitialized value warning +<BR>Use rmtree function to delete tmp files (FY 12/04/12)
# 1.0.1 Multiple minor improvements (FY 26/01/12)
# 1.0.0 New script for GO term enrichment analysis of a protein set (FY 10/06/11)
