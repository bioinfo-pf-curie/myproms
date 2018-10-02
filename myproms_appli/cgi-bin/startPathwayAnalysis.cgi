#!/usr/local/bin/perl -w

################################################################################
# startPathwayAnalysis.cgi       1.0.4                                         #
# Authors: P. Poullet, G. Arras, S.Liva (Institut Curie)                       #
# Contact: myproms@curie.fr                                                    #
# Fetch and provide proteins and parameters for Pathway enrichment analysis    #
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
use File::Path qw(make_path remove_tree);

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my $dbh=&promsConfig::dbConnect;
my $expID = param('id_exp');
my $projectID=&promsMod::getProjectID($dbh,$expID,'EXPERIMENT');

if(param('start')){
    &processForm;
    exit;
}

#############################
# Generating selection tree #
#############################
my %analysisClass=(-1=>'no_scan',0=>'no_val',1=>'part_val',2=>'val');
my %treeOptions;
$treeOptions{'BUTTON'}='font-size:9px;width:55px';
%{$treeOptions{'CHECKBOX'}}=('experiment'=>1,'gel2d'=>1,'spot'=>1,'sample'=>1,'analysis'=>1); # all items are potentially checkable
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
    }
    push @experimentTree,\@sampleTree;
    $treeOptions{'DISCHECK'}{"sample:$sampID"}=1 unless $chkBoxSamp>0;
}

$sthG2D->finish;
$sthSpot->finish;
$sthSamp->finish;
$sthAna->finish;

my $disableSubmitStrg = '';

###########################################################
# Fetching categories for background population selection #
###########################################################
my $sthClass = $dbh->prepare("SELECT ID_CLASSIFICATION,NAME FROM CLASSIFICATION WHERE ID_PROJECT=?");
my $sthCategory = $dbh->prepare("SELECT ID_CATEGORY,NAME FROM CATEGORY WHERE ID_CLASSIFICATION=?");
my $categoryStrg="<OPTION value=\"0\">-= Select =-</OPTION>\n";
$sthClass->execute($projectID);
while (my ($classID,$className) = $sthClass->fetchrow_array){
    $categoryStrg .= "<OPTGROUP label=\"Theme: $className\">\n";
    $sthCategory->execute($classID);
    while (my ($categoryID,$categoryName) = $sthCategory->fetchrow_array){
		$categoryStrg .= "<OPTION value=$categoryID>$categoryName</OPTION>\n";
    }
}
$sthClass->finish;
$sthCategory->finish;
$dbh->disconnect;

###########################
###### Starting HTML ######
###########################
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
print qq
|function checkForm(myForm){

    var protSetRadioValue = getRadioValue(myForm.chProtSet);
    myForm.checkedItems.value=getCheckedItems();
    if(!myForm.anaName.value){
	    alert('Type a name for this Pathway analysis.');
	    return false;
    }
    /*else if (myForm.FDR.value > 100 \|\| myForm.FDR.value <= 0){
        alert('Enter a valid FDR.');
        return false;
    }*/
    else if (myForm.pval.value > 1 \|\| myForm.pval.value <= 0){
        alert('Enter a valid p-value threshold.');
        return false;
    }
    else if(!myForm.checkedItems.value && protSetRadioValue == 'analysis'){
        alert('Select at least 1 item for the analysis.');
        return false;
    }
    else if(myForm.protSet[myForm.protSet.selectedIndex].value==0 && protSetRadioValue == 'category') {
        alert('Select a list.');
        return false;
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
function actionOnSelect(){}
function actionOnCheck(){}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER><FONT class="title">Pathway Enrichment Analysis</FONT><BR><BR><BR>
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
		    <TH valign='top' align='right'>&nbsp;&nbsp;Advanced Parameters :</TH><TD bgcolor=$lightColor>
			&nbsp;Statistical settings:<BR>
			<!--&nbsp;&nbsp;&nbsp;&bull;&nbsp;Control FDR at <INPUT type='text' name='FDR' value='1' size=2>%
			<BR>-->
			&nbsp;&nbsp;&nbsp;&bull;&nbsp;Use a p-value threshold: <INPUT type='text' name='pval' size=3 value='0.01'>
			<BR>
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

print qq
|</TABLE>
</CENTER>
</BODY>
</HTML>
|;
exit;

##insert database and refresh to runAndDisplayPathwayAnalysis.cgi
sub processForm {
    #######################
    # Fetching parameters #
    #######################
    my $userID = $ENV{'REMOTE_USER'};
    my $name = param('anaName');
    my $description = param('description');
    my $totItemStrg = param('checkedItems');

    ##INSERT PATHWAY_ANALYSIS
    my $sthInsertPathAna=$dbh->prepare("INSERT INTO PATHWAY_ANALYSIS(ID_EXPERIMENT, ID_CATEGORY, NAME, DES, PARAM_STRG, STATUS, RECORD_DATE, UPDATE_DATE, UPDATE_USER)
				      values (?, ?, ?, ?, ?, ?, NOW(), NOW(),?)");
    my $catID = param('protSet')? param('protSet') : undef;
    #my $FDRthreshold=param('FDR');
    my $FDRthreshold="1";
    my $PVALthreshold=param('pval');
    my $paramStrg="FDR=$FDRthreshold;PVAL=$PVALthreshold";
    $sthInsertPathAna->execute($expID, $catID, $name, $description, $paramStrg,'0', $userID);
    $sthInsertPathAna->finish;
    my $pathAnaID=$dbh->last_insert_id(undef,undef,'PATHWAY_ANALYSIS','ID_PATHWAY_ANALYSIS');
    if (param('chProtSet') eq 'category') {
		$dbh->commit;
    }
    elsif ($totItemStrg) {
		$sthInsertPathAna->finish;
		my %analysisList;
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
						  WHERE ID_ANALYSIS=?'), # after submit only
			'experiment'=>$dbh->prepare('SELECT ID_ANALYSIS FROM SAMPLE S, ANALYSIS A
						WHERE S.ID_SAMPLE=A.ID_SAMPLE
						AND S.ID_EXPERIMENT=?')
		);
		foreach my $itemInfo (split(/\+/,$totItemStrg)){ # ex: experiment:1,2,3+sample:1,3+analysis:2
			my ($chkItem,@listIDs)=split(/[:,]/,$itemInfo);
			foreach my $itID (@listIDs) {
				$sthAnaList{$chkItem}->execute($itID);
				while (my ($analysisID)=$sthAnaList{$chkItem}->fetchrow_array){
					$analysisList{$analysisID}=1;
				}
			}
		}
		##INSERT PATHWAYANA_ANALYSIS
		if (scalar(keys %analysisList)) {
			my $sthInsertPathAnaAna=$dbh->prepare("INSERT INTO PATHWAYANA_ANALYSIS(ID_ANALYSIS, ID_PATHWAY_ANALYSIS) values (?,?)");
			foreach my $anaID (sort{$a <=> $b} keys %analysisList) {
				$sthInsertPathAnaAna->execute($anaID, $pathAnaID);
			}
			$sthInsertPathAnaAna->finish;
		}
		$dbh->commit;
    }

    $dbh->disconnect;

     ###>Forking to launch Analysis<###
    if (defined $pathAnaID) {
	my $projectPath="$promsPath{pathAna}/project_$projectID";
        my $tmpPathwayPath="$promsPath{tmp}/pathway_analysis";
        my $tmpPathwayPathID="$tmpPathwayPath/$pathAnaID";

        make_path("$projectPath", { verbose => 0, mode => 0755}) unless -e "$projectPath";
        make_path("$tmpPathwayPath", { verbose => 0, mode => 0755}) unless -e "$tmpPathwayPath";
        make_path("$tmpPathwayPathID", { verbose => 0, mode => 0755}) unless -e "$tmpPathwayPathID";

	my $childPathwayAna= fork;
	unless ($childPathwayAna) { # child here
		#>Disconnecting from server
		open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
		open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
		open STDERR, '>/dev/null' or die "Can't open /dev/null: $!";
		system "./runPathwayAnalysis.pl $pathAnaID $catID  2> $tmpPathwayPathID/error.txt";
		exit;
	}
    }

    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    print qq
|<SCRIPT language="JavaScript">
	parent.itemFrame.location='$promsPath{cgi}/openProject.cgi?ACT=experiment&VIEW=functAna&branchID=pathwayanalysis:$pathAnaID';
</SCRIPT>
</BODY>
</HTML>
|;
}
####>Revision history<####
# 1.0.4 Comment checkForm FDR test to prevent launch without selection of analysis/category set (GA 01/03/17)
# 1.0.3 comment FDR statistics, put defaut p-value to 0.05  (SL 20/10/15)
# 1.0.2 add fork new script runPathwayAnalysis.pl, replace runAndDisplayPathwayAnalysis (SL 03/09/15)
# 1.0.1 Remove extra header in processForm function (PP 13/11/14)
# 1.0.0 New script for Pathway enrichment analysis of a protein set, insert database and launch runAndDisplayPathwayAnalysis.cgi (SL 16/10/14)
