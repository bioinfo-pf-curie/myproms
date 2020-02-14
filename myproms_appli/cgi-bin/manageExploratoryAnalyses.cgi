#!/usr/local/bin/perl -w

################################################################################
# manageExploratoryAnalyses.cgi       1.1.8                                    #
# Authors: P. Poullet, S.Liva (Institut Curie)                                 #
# Contact: myproms@curie.fr                                                    #
# Display PCA & clustering analysis                                            #
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
use promsQuantif;
use POSIX qw(strftime); # to get the time
use File::Copy;
use File::Path qw(rmtree);

#print header(-charset=>'utf-8');warningsToBrowser(1);#DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my %proteinQuantifFamilies=&promsQuantif::getProteinQuantifFamilies;
my %features=('protQuant'=>'Protein quantifications','pepQuant'=>'Peptide quantifications','pepCount'=>'Peptide count');
my %pepTypeDesc=('NUM_PEP_USED'=>'All','DIST_PEP_USED'=>'Distinct','RAZ_UNI_PEP'=>'Razor + unique','UNIQUE_PEP'=>'Unique','IDENT_PEP'=>'Identified');
my %anaTypeDesc=('cluster'=>'Clustering','clusterPEP'=>'Clustering','PCA'=>'PCA', 'PCAPEP'=>'PCA');
my $action = param('ACT')? param('ACT') : 'summary';
my $explorID = param('explorID');
my $experimentID = param('experimentID')? param('experimentID') : "";
my $ajax = param('AJAX')? param('AJAX') : '';
my $dbh = &promsConfig::dbConnect;
my $projectID=&promsMod::getProjectID($dbh,$experimentID,'EXPERIMENT');
my ($explorName,$anaType,$description,$status,$listID,$listExclusion,$filterStrg,$paramStrg) = $dbh -> selectrow_array("SELECT NAME,ANA_TYPE,DES,STATUS,ID_CATEGORY,CAT_EXCLUSION,FILTER_LIST,PARAM_LIST FROM EXPLORANALYSIS WHERE ID_EXPLORANALYSIS = $explorID");

if ($action eq "delete") {

    print header(-'content-encoding'=>'no',-charset=>'utf-8');
    warningsToBrowser(1); # DEBUG
    if ($anaType eq ('PCA') || $anaType eq ('cluster')) {
        #my $sthDeleteExplorQuantif = $dbh -> do("DELETE from EXPLORANA_QUANTIF where ID_EXPLORANALYSIS = $explorID");
        #my $sthDeleteAnnotSet = $dbh -> do("DELETE from ANNOTATIONSET where ID_EXPLORANALYSIS = $explorID");
        $dbh->do("DELETE from EXPLORANA_QUANTIF where ID_EXPLORANALYSIS = $explorID");
        $dbh->do("DELETE from ANNOTATIONSET where ID_EXPLORANALYSIS = $explorID");
    }
    else {
        #my $sthDeleteExplorAnaAna = $dbh->do("DELETE FROM EXPLORANA_ANA where ID_EXPLORANALYSIS = $explorID");
        $dbh->do("DELETE from ANNOTATIONSET where ID_EXPLORANALYSIS = $explorID");
        $dbh->do("DELETE FROM EXPLORANA_ANA where ID_EXPLORANALYSIS = $explorID");
    }
    #my $sthDeleteExplorAnalysis = $dbh -> do("DELETE from EXPLORANALYSIS where ID_EXPLORANALYSIS = $explorID");
    $dbh -> do("DELETE from EXPLORANALYSIS where ID_EXPLORANALYSIS = $explorID");
    $dbh -> commit;
    $dbh -> disconnect;
    my $pathToFile = (-e "$promsPath{explorAna}/project_$projectID/$explorID")? "$promsPath{explorAna}/project_$projectID/$explorID" : "$promsPath{tmp}/exploratory_analysis/$explorID";
    rmtree($pathToFile);
    print qq
|<HTML>
<HEAD>
<SCRIPT LANGUAGE="JavaScript">
top.promsFrame.selectedAction = 'summary';
parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=experiment:$experimentID&branchID=experiment:$experimentID&VIEW=explorAna";
</SCRIPT>
</HEAD>
</HTML>
|;
    exit;

}

if (param('submit')) {

    my $desc = param('desc')? param('desc') : undef;
    my $explorName = param('explorName');

    my $sthUpdateExplorAna = $dbh -> prepare("UPDATE EXPLORANALYSIS set NAME = ?, DES = ? where ID_EXPLORANALYSIS = ?") or die "Cannot prepare: " . $dbh->errstr();
    $sthUpdateExplorAna -> execute($explorName, $desc,$explorID) or die "Cannot execute: " . $sthUpdateExplorAna->errstr();
    $sthUpdateExplorAna -> finish;

    $dbh -> commit;
    $dbh -> disconnect;
    
    ####>Updating all frames<###
    print header(-'content-encoding'=>'no',-charset=>'utf-8');
    warningsToBrowser(1); # DEBUG
    print qq
|<HTML>
<HEAD>
<TITLE>Update All Frames</TITLE>
<SCRIPT LANGUAGE="JavaScript">
	top.promsFrame.selectedAction = 'summary';
	parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=experiment:$experimentID&branchID=exploranalysis:$explorID&VIEW=explorAna";
</SCRIPT>
</HEAD>
</HTML>
|;
    exit;
}

$description='' unless $description;
$filterStrg='' unless $filterStrg;
$paramStrg='' unless $paramStrg;
my ($isAnnot) = $dbh -> selectrow_array("SELECT COUNT(AN.ID_ANNOTATIONSET) from EXPLORANALYSIS EA,ANNOTATIONSET AN where EA.ID_EXPLORANALYSIS = $explorID and EA.ANA_TYPE = '$anaType' and EA.ID_EXPLORANALYSIS = AN.ID_EXPLORANALYSIS");

my $statusStrg = ($status == -2)? "<IMG src=$promsPath{images}/bad.gif>&nbsp;&nbsp;<FONT color=\"red\">***ERROR***</FONT>" : ($status == -1)? "<FONT color=\"orange\"><B>0n-going analysis</B></FONT>" : "<IMG src=$promsPath{images}/good.gif>&nbsp;<FONT color=\"green\">Finished</FONT>";
#my %convertItem = ("obs" => "Observation", "quanti" => "Quantification values", "prot" => "proteins", "peptide" => "Peptide number", "PV" => "pValue", "FC" => "foldChange");
my $annotLabel=($anaType eq 'PCA')? 'highlights' : 'annotations';

### START HTML
print header(-charset=>'utf-8');
warningsToBrowser(1);
#my $anaStrg = ($anaType eq 'cluster')? 'Clustering' : $anaType;
my $anaStrg=$anaTypeDesc{$anaType};
my $title = ($action eq 'summary')? "$anaStrg <FONT color='#DD0000'>$explorName</FONT>" : "Editing $anaStrg <FONT color='#DD0000'>$explorName</FONT>" ;
print qq
|<HTML>
<HEAD>
<TITLE>ExplorAnalysis</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">|;
if ($action eq 'edit') {
    print qq
|<SCRIPT langage="Javascript">
function checkForm(myForm) {
	if (!myForm.explorName.value) {
		alert('A name is expected for the analysis');
		return false;
	}
	return true;
}
|;
}
else {
    print qq
|<SCRIPT langage="Javascript">
function displayQuantiList() {
	var quantifDiv=document.getElementById('quantiList');
	var quantifBut=document.getElementById('buttonList')
	if (quantifDiv.style.display == 'none') {
		quantifDiv.style.display = '';
		quantifBut.value='Hide list';
		//document.getElementById('buttonList').innerHTML = '<INPUT type="button" value="Hide list" onclick="displayQuantiList(document.displayMenu,\\'hide\\')">';
	}
	else {
		quantifDiv.style.display = 'none';
		quantifBut.value='Show list';
		//document.getElementById('quantiList').style.display = 'none';
		//document.getElementById('buttonList').innerHTML = '<INPUT type="button" value="Display list" onclick="displayQuantiList(document.displayMenu,\\'show\\')">';
	}
}
|;
	print "setTimeout(function(){parent.optionFrame.location.reload();},5000);\n" if $status==-1;
}
print qq
|</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">$title</FONT><BR><BR>
|;
my ($light,$dark)=&promsConfig::getRowColors;

my (%quantifList,%dataParams,%dataFilters,%refQuantif,%analysisList,$focusStrg,$isNormalized);
if ($action eq 'summary') {
	my $modifID;
	
    if ($anaType eq 'PCA' || $anaType eq 'cluster') {
        my $sthDesQuantif = $dbh -> prepare("SELECT EQ.ID_QUANTIFICATION,EQ.TARGET_POS,CONCAT(D.NAME,' > ',Q.NAME),Q.QUANTIF_ANNOT,ID_MODIFICATION FROM EXPLORANA_QUANTIF EQ
                                                    INNER JOIN QUANTIFICATION Q ON EQ.ID_QUANTIFICATION = Q.ID_QUANTIFICATION
                                                    INNER JOIN DESIGN D ON Q.ID_DESIGN = D.ID_DESIGN
                                                    WHERE EQ.ID_EXPLORANALYSIS = $explorID");
        $sthDesQuantif->execute;
        while (my($quantifID,$targetPos,$quantifPath,$quantifAnnot,$modID) = $sthDesQuantif -> fetchrow_array) {
            $modifID=$modID;
            if ($targetPos) { # RATIO family
                my (%labelingInfo,%stateInfo);
                &promsQuantif::extractQuantificationParameters($dbh,$quantifAnnot,\%labelingInfo,\%stateInfo);
                if ($labelingInfo{'RATIOS'}) {
                    my ($testCondID,$refCondID)=split(/\//,$labelingInfo{'RATIOS'}[$targetPos-1]);
                    my $normTag='';
                    if ($testCondID=~/%/) { # Super ratio
                        $normTag='°';
                        $testCondID=~s/%\d+//;
                        $refCondID=~s/%\d+//;
                    }
                    $quantifList{$quantifID.'_'.$targetPos}="$quantifPath : $stateInfo{$testCondID}{NAME}$normTag/$stateInfo{$refCondID}{NAME}$normTag";
                }
                else { # MQ
					my ($numBioRep,$quantiObsIDs,$condID)=split(',',$labelingInfo{'STATES'}[$targetPos-1]);
                    $quantifList{$quantifID.'_'.$targetPos}="$quantifPath : $stateInfo{$condID}{NAME}";
                }
            }
        }
        $sthDesQuantif->finish;

        my $sthNonDesQuantif=$dbh->prepare("SELECT EQ.ID_QUANTIFICATION,EQ.TARGET_POS,CONCAT(S.NAME,' > ',A.NAME,' > ',Q.NAME),QUANTIF_ANNOT,ID_MODIFICATION FROM EXPLORANA_QUANTIF EQ
                                                    INNER JOIN QUANTIFICATION Q ON EQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND ID_DESIGN IS NULL
                                                    INNER JOIN ANA_QUANTIFICATION AQ ON Q.ID_QUANTIFICATION=AQ.ID_QUANTIFICATION
                                                    INNER JOIN ANALYSIS A ON AQ.ID_ANALYSIS=A.ID_ANALYSIS
                                                    INNER JOIN SAMPLE S ON A.ID_SAMPLE=S.ID_SAMPLE
                                                    WHERE EQ.ID_EXPLORANALYSIS = $explorID");
        $sthNonDesQuantif->execute;
        while (my($quantifID,$ratioPos,$quantifPath,$quantifAnnot,$modID) = $sthNonDesQuantif -> fetchrow_array) {
            $modifID=$modID;
            if ($ratioPos) { # RATIO family
                my (%labelingInfo,%stateInfo);
                &promsQuantif::extractQuantificationParameters($dbh,$quantifAnnot,\%labelingInfo,\%stateInfo);
                my ($testStatePos,$refStatePos)=split(/\//,$labelingInfo{'RATIOS'}[$ratioPos-1]);
                $quantifList{$quantifID.'_'.$ratioPos}="$quantifPath : $stateInfo{$testStatePos}{NAME}/$stateInfo{$refStatePos}{NAME}";
            }
            else {
                $quantifList{$quantifID.'_0'}=$quantifPath;
            }
        }
        $sthNonDesQuantif->finish;
    }
    else {
        my $sthExplorAnaAna = $dbh->prepare("SELECT EA.ID_ANALYSIS,EA.GROUP_POS,CONCAT(S.NAME,' > ',A.NAME) FROM EXPLORANA_ANA EA
                                            INNER JOIN ANALYSIS A ON EA.ID_ANALYSIS=A.ID_ANALYSIS
                                            INNER JOIN SAMPLE S ON A.ID_SAMPLE=S.ID_SAMPLE
                                            WHERE EA.ID_EXPLORANALYSIS=$explorID");
        $sthExplorAnaAna->execute;
        while (my($analysisID,$groupPos,$analysisPath) = $sthExplorAnaAna -> fetchrow_array) {
            $analysisList{$analysisID}=$analysisPath;
        }
    }
    
	foreach my $param (split("//", $paramStrg)) {
		next if !$param; # strg starts with //
		my ($item, $itemValue) = split("=", $param);
		$dataParams{$item} = $itemValue;
	}
	foreach my $filter (split("//", $filterStrg)) {
		next if !$filter; # strg starts with //
		my ($item, $itemValue) = split("=", $filter);
		$dataFilters{$item}=$itemValue;
	}

    if ($dataParams{'normalization'}) {
		foreach my $strg (split(/[;:]/, $dataParams{'normalization'})) {
			my ($testQuant, $refQuant)=split(/%/, $strg);
			$refQuantif{$testQuant}=$refQuant;
		}
		$focusStrg='Normalized ';
		$isNormalized=1;
    }
	if ($modifID) {
		my ($modifName)=$dbh->selectrow_array("SELECT PSI_MS_NAME FROM MODIFICATION WHERE ID_MODIFICATION=$modifID");
		$focusStrg.="$modifName-Proteins";
	}
	else {$focusStrg='Proteins';}

}
else { # edit
	print qq
|<FORM NAME="displayMenu" method="POST" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="explorID" value="$explorID" />
|;
}
print qq
|<TABLE border="0" cellspacing="2" cellpadding="2" bgcolor="$dark">
|;
my $titleHeader = ($anaType eq ('PCA') || $anaType eq ('PCAPEP')) ? "PCA name" : "Cluster name";
if ($action eq 'edit') {
	print qq|<INPUT type="hidden" name="experimentID" value="$experimentID" />|;
	print qq
|<TR><TH align=right nowrap>$titleHeader :</TH><TD bgcolor="$light"><INPUT type="text" name="explorName" size="50" maxlength="100" value="$explorName"/></TD></TR>
<TR><TH align=right valign="top" nowrap>&nbsp;Description :</TH><TD bgcolor="$light"><TEXTAREA name="desc" rows="2" cols="50">$description</TEXTAREA></TD>
<TR><TH colspan="2"><INPUT type="submit" name="submit" value="Update"></TD></TR>
|;
}
else {#summary
	$description=&promsMod::HTMLcompatible($description);
	my $measureStrg='';
	if ($dataParams{'quantCode'}) {
		foreach my $refMeas (@{$proteinQuantifFamilies{'MEASURES'}{$dataParams{quantFam}}}) {
			if ($refMeas->[0] eq $dataParams{'quantCode'}) {
				$measureStrg=" ($refMeas->[1])";
				last;
			}
		}
	}

	print qq
|<TR><TH align=right width=150 nowrap>$titleHeader :</TH><TD bgcolor="$light">$explorName</TD></TR>
<TR><TH align=right valign="top" nowrap>&nbsp;Description :</TH><TD bgcolor="$light">$description</TD></TR>
<TR><TH align=right valign="top" nowrap>&nbsp;Focus :</TH><TD bgcolor="$light">$focusStrg</TD></TR>
<TR><TH align=right valign="top" nowrap>&nbsp;Feature :</TH><TD bgcolor="$light">$features{$dataParams{feature}}&nbsp;:&nbsp;$dataParams{'quantFam'}|;
if ($anaType eq ('PCA') || $anaType eq ('cluster')) {
    print qq|&nbsp;[$proteinQuantifFamilies{NAME}{$dataParams{quantFam}}$measureStrg]|;
}
print qq|</TD></TR>|;
if ($anaType eq ('PCA') || $anaType eq ('cluster')) {
    print qq|<TR><TH align=right valign="top" nowrap>&nbsp;Data transform :</TH><TD nowrap bgcolor=$light>$dataParams{dataTrans}</TD></TR>|;
}
print qq|
<TR><TH align=right valign="top" nowrap>&nbsp;Data filtering :</TH><TD nowrap bgcolor=$light>|;
    if ($anaType eq ('PCA') || $anaType eq ('cluster')) {
        if ($dataParams{'quantFam'} eq 'RATIO') {
            $dataFilters{'FC'}=~s/:(\d+)//;
            my $fcOccur=$1 || 1; # back compatibility with older Analyses with no fold change occurence
            my $pvOccur;
            unless ($isNormalized) {
                $dataFilters{'PV'}=~s/:(\d+)//;
                $pvOccur=$1 || 1; # back compatibility with older Analyses with no p-value occurence
            }
            $dataFilters{'PEP'}='all:1' unless $dataFilters{'PEP'}; # back compatibility with older Analyses
            $dataFilters{'INF'}=34 if !defined $dataFilters{'INF'}; # back compatibility with older Analyses
            print "&nbsp;&bull;&nbsp;Abs. fold change &ge; <B>$dataFilters{FC}</B> in at least <B>$fcOccur</B> quantification(s)&nbsp;<BR>\n";
            if ($dataFilters{'INF'} < 0) {print "&nbsp;&bull;&nbsp;Infinite ratios are treated as <B>missing values</B>&nbsp;<BR>\n";}
            else {
                my ($infRatioStrg,$protStrg)=($dataFilters{'INF'}==100)? ('All','') : ($dataFilters{'INF'}==0)? ('No','') : ("$dataFilters{INF}%",' / protein');
                print "&nbsp;&bull;&nbsp;<B>$infRatioStrg</B> infinite ratios allowed$protStrg&nbsp;<BR>\n";
            }
            print "&nbsp;&bull;&nbsp;p-value &le; <B>$dataFilters{PV}</B> in at least <B>$pvOccur</B> quantification(s)&nbsp;<BR>\n" unless $isNormalized;
        }
        my ($pepType,$numPep)=split(':',$dataFilters{PEP});
        $pepType='NUM_PEP_USED' if $pepType eq 'all'; # compatibility with old explorana
        print "&nbsp;&bull;&nbsp;<B>$pepTypeDesc{$pepType}</B> peptides &ge; <B>$numPep</B><BR>\n";
    }
	
    $dataFilters{'NA'}=34 if !defined $dataFilters{'NA'}; # back compatibility with older Analyses
	my $missingValuesStrg=($dataFilters{'NA'}==100)? 'All' : ($dataFilters{'NA'}==0)? 'No' : ($dataFilters{'NA'}==-1)? 'Not authorized' : "$dataFilters{NA}%";
    my $focusStrg = ($anaType eq ('PCA') || $anaType eq ('cluster'))? "protein" : "peptide";
	print "&nbsp;&bull;&nbsp;<B>$missingValuesStrg</B> missing values allowed / $focusStrg&nbsp;<BR>\n";
    if (defined $dataFilters{'nbNA'}) {
        my $pcNA="";
        if ($anaType eq ('PCA') || $anaType eq ('cluster')) {
            $pcNA=sprintf "%.2f",100 * $dataFilters{'nbNA'}/($dataFilters{'nbAllProt'} * scalar (keys %quantifList));
        }
        else {
            $pcNA=sprintf "%.2f",100 * $dataFilters{'nbNA'}/($dataFilters{'nbAllProt'} * scalar (keys %analysisList));
        }
		print "&nbsp;&bull;&nbsp;<B>$pcNA%</B> missing values ($dataFilters{'nbNA'} values imputed)&nbsp;<BR>\n";
	}
	
    if ($anaType eq ('PCA') || $anaType eq ('cluster')) {
        if (defined $dataFilters{'nbProt'}) {
            my ($themeName,$listName) = $dbh -> selectrow_array("SELECT CL.NAME,CA.NAME FROM CATEGORY CA,CLASSIFICATION CL WHERE CA.ID_CLASSIFICATION=CL.ID_CLASSIFICATION AND ID_CATEGORY = $listID");
            if ($listExclusion == 1) {
                print "&nbsp;&bull;&nbsp;<B>$dataFilters{nbProt}</B> Proteins <B>excluded</B> from list <B>$themeName &gt; $listName</B>&nbsp;<BR>\n";
            }
            else {
                #my ($protNumber) = $dbh -> selectrow_array("SELECT COUNT(ID_PROTEIN) FROM CATEGORY_PROTEIN WHERE ID_CATEGORY = $listID");
                print "&nbsp;&bull;&nbsp;<B>Only $dataFilters{nbProt}</B> proteins in List <B>$themeName &gt; $listName</B>&nbsp;<BR>\n";
            }
        }
        if ($dataFilters{'nbAllProt'}) {
            print "&nbsp;&bull;&nbsp;<B>$dataFilters{nbAllProt}</B> proteins used&nbsp;<BR>\n";
        }
        #else {
        #	print "&nbsp;&bull;$convertItem{$itemValue}<BR>";
        #}
    }
	print "</TD></TR>\n";
    if (defined($dataParams{'groupMethod'})) {
        print qq
|<TR><TH align=right valign="top" bgcolor=$dark nowrap>&nbsp;Parameters :</TH><TD nowrap bgcolor=$light>
&nbsp;&bull;&nbsp;Grouping method: <b>$dataParams{'groupMethod'}<b> 
</TD></TR>|;
    }
    
	if ($anaType eq ('cluster') || $anaType eq ('clusterPEP')) {
		print qq
|<TR><TH align=right valign="top" bgcolor=$dark nowrap>&nbsp;Clustering parameters :</TH><TD nowrap bgcolor=$light>
&nbsp;&bull;&nbsp;Method: <B>$dataParams{method}</B><BR>
&nbsp;&bull;&nbsp;Metric: <B>$dataParams{metric}</B></TD></TR>
|;
	}
    
    my $strgList = ($anaType eq ('PCA') || $anaType eq ('cluster'))? "Quantifications" : "Analyses" ;
	print qq
|<TR><TH align=right bgcolor=$dark nowrap>Status :</TH><TD nowrap bgcolor=$light>&nbsp;&nbsp;$statusStrg</TD></TR>
<TR><TH align=right valign="top" bgcolor=$dark nowrap>$strgList used :</TH><TD nowrap bgcolor=$light>
	<INPUT id="buttonList" type="button" value="Show list" onclick="displayQuantiList(document.displayMenu,'show')">
	<DIV id="quantiList" style="display:none;">
		<TABLE cellpadding=0>
|;
	if (keys %refQuantif) {
		print "<TR><TH colspan=2 align=left><U>Normalization:</U></TH></TR>\n";
		my $count=0;
		foreach my $quantifID (sort{$a cmp $b} keys %quantifList)  {
			if ($refQuantif{$quantifID}) {
				$count++;
				print "<TR><TD align=right nowrap>&nbsp;#$count.</TD><TD nowrap>&nbsp;<B>$quantifList{$quantifID} - $quantifList{$refQuantif{$quantifID}}</B>&nbsp;</TD></TR>\n";
				delete $quantifList{$quantifID};
				delete $quantifList{$refQuantif{$quantifID}};
			}
		}
		print "<TR><TH colspan=2 align=left><U>Test:</U></TH></TR>\n";
	}
	if (keys %quantifList) {
		my $count=0;
		foreach my $quantifInfo (sort{$a cmp $b} values %quantifList) {
			$count++;
			print "<TR><TD align=right nowrap>&nbsp;#$count.</TD><TD nowrap>&nbsp;<B>$quantifInfo</B>&nbsp;</TD></TR>\n";
		}
	}
    
    if (keys %analysisList) {
        my $count=0;
		foreach my $analysisInfo (sort{$a cmp $b} values %analysisList) {
			$count++;
			print "<TR><TD align=right nowrap>&nbsp;#$count.</TD><TD nowrap>&nbsp;<B>$analysisInfo</B>&nbsp;</TD></TR>\n";
		}
    }
	print qq
|</TABLE></DIV>
</TD></TR>
|;
}
print qq
|</TABLE>
<BR>
|;
if ($action eq 'edit'){
	print "</FORM>\n";
}
elsif ($action eq 'summary' && $isAnnot) {
    print qq
|<BR>
<FONT class="title2">List of graphical $annotLabel used:</FONT>
<TABLE cellspacing=0>
|;
	my %convertGO = ("P" => "Biological process", "C" => "Cellular component", "F" => "Molecular function");
	my %annotSet;
	my $sthSelectAnnotation = $dbh -> prepare("SELECT NAME,RANK,ANNOT_TYPE,ANNOT_LIST FROM ANNOTATIONSET WHERE ID_EXPLORANALYSIS = $explorID");
	$sthSelectAnnotation -> execute;
	while (my($annotName, $annotRank, $annotType, $annotList) = $sthSelectAnnotation -> fetchrow_array) {
		@{$annotSet{$annotType}{$annotName}} = ($annotRank,$annotList);
	}
	$sthSelectAnnotation -> finish;
    
    my $strgQuantAna = ($anaType=~/PEP/)? "analyses" : "quantifications";
	foreach my $annotType (sort{$a cmp $b} keys %annotSet) {
		if ($annotType =~ /^property:/) { # PCA & cluster
			my $sthProp=$dbh->prepare("SELECT NAME FROM PROPERTY WHERE ID_PROPERTY=?");
			my $propLabel=($annotType eq 'property:O')? 'property' : 'treatment';
			print "<TR><TD>&nbsp;</TD></TR>\n<TR bgcolor=\"$dark\"><TH align=\"left\">&nbsp;Biosample $propLabel $annotLabel for $strgQuantAna:&nbsp;</TH></TR>\n";
			
            foreach my $annotName (sort{$annotSet{$annotType}{$a}[0] <=> $annotSet{$annotType}{$b}[0]} keys %{$annotSet{$annotType}}) {
                if ($anaType eq 'PCA') {
					my ($propID,$propertyValues)=split(':=:',$annotSet{$annotType}{$annotName}[1]);
					$propID=~s/#//;
					$sthProp->execute($propID);
					my ($propName)=$sthProp->fetchrow_array;
					if ($annotType eq 'property:O') {
						print "<TR bgcolor=\"$light\"><TD nowrap>&nbsp;$propName = \"$propertyValues\"&nbsp;</TD></TR>\n";
					}
					elsif ($annotType eq 'property:T') {
						print "<TR bgcolor=\"$light\"><TD nowrap>&nbsp;$propName = \"",&promsQuantif::convertTreatmentToText($propertyValues),"\"&nbsp;</TD></TR>\n";
					}
				}
				else { # cluster
					my ($propID)=($annotName=~/^#(\d+)/);
					$sthProp->execute($propID);
					my ($propName)=$sthProp->fetchrow_array;
					print "<TR bgcolor=\"$light\"><TD nowrap>&nbsp;$propName&nbsp;</TD></TR>\n";
				}
				print "<TR bgcolor=\"$dark\"><TD></TD></TR>\n";
			}
			$sthProp->finish;
		}
        elsif ($annotType =~ /^ms_sample:/) {
            print "<TR><TD>&nbsp;</TD></TR>\n<TR bgcolor=\"$dark\"><TH align=\"left\">&nbsp;Biosample other annotations for $strgQuantAna:&nbsp;</TH></TR>\n";
            print "<TR bgcolor=\"$light\"><TD nowrap>&nbsp;MS Sample&nbsp;</TD></TR>\n";
        }
		elsif ($annotType =~ /prot:GO/) { # PCA & cluster
			print "<TR><TD>&nbsp;</TD></TR>\n<TR><TH align=\"left\" bgcolor=\"$dark\">&nbsp;GO $annotLabel for proteins:&nbsp;</TH></TR>\n";
			my $sthGO=$dbh -> prepare("SELECT NAME from GO_ANALYSIS where ID_GOANALYSIS = ?");
			foreach my $annotName (sort{$annotSet{$annotType}{$a}[0] <=> $annotSet{$annotType}{$b}[0]} keys %{$annotSet{$annotType}}) {
				if ($anaType eq 'PCA') {
					my @annotTab = split(",", $annotName);
					my $nbItem = scalar(@annotTab);
					if ($nbItem == 3) {
						my ($goAnaID,$goAspect,$goID) = split(",", $annotName);
						$goAnaID=~s/^#//;
						$sthGO->execute($goAnaID);
						my ($analysisName) = $sthGO->fetchrow_array;
						print "<TR bgcolor=\"$light\"><TD nowrap>&nbsp;$analysisName > $convertGO{$goAspect} > $annotSet{$annotType}{$annotName}[1] [$goID]&nbsp;</TD></TR>\n";
					}
					else {
						my ($goAnaID,$goAspect,$goID,$binInfo) = split(",", $annotName);
						my $binValue=(split(":",$binInfo))[1];
						$goAnaID=~s/^#//;
						$sthGO->execute($goAnaID);
						my ($analysisName) = $sthGO->fetchrow_array;
						if (!$analysisName) {
							my ($parentGOAnaID) = $dbh->selectrow_array("SELECT ID_PARENT_GOANA FROM GO_ANALYSIS WHERE ID_GOANALYSIS=$goAnaID");
							$sthGO->execute($parentGOAnaID);
							($analysisName) = $sthGO->fetchrow_array;
						}

						print "<TR bgcolor=\"$light\"><TD nowrap>&nbsp;$analysisName > $convertGO{$goAspect} > $annotSet{$annotType}{$annotName}[1] [$binValue]&nbsp;</TD></TR>\n";
					}
				}
				else { # cluster
					my ($goAnaID,$goAspect,$goString) = split(':=:',$annotSet{$annotType}{$annotName}[1]);
					$goAnaID=~s/^#//;
					$sthGO->execute($goAnaID);
					my ($analysisName) = $sthGO->fetchrow_array;
					if (!$analysisName) {
						my ($parentGOAnaID) = $dbh->selectrow_array("SELECT ID_PARENT_GOANA FROM GO_ANALYSIS WHERE ID_GOANALYSIS=$goAnaID");
						$sthGO->execute($parentGOAnaID);
						($analysisName) = $sthGO->fetchrow_array;
					}
					$analysisName='?' unless $analysisName;
					print "<TR bgcolor=\"$light\"><TD nowrap>&nbsp;$annotName: $analysisName > $convertGO{$goAspect} >&nbsp;<BR>";
					foreach my $termStrg (split(':@:',$goString)) {
						my ($goID,$goTerm)=split('@',$termStrg);
						print "&nbsp;&nbsp;&nbsp;&nbsp;&bull;$goTerm [$goID]<BR>\n";
					}
				}
				print "<TR bgcolor=\"$dark\"><TD></TD></TR>\n";
			}
			$sthGO->finish;
		}
		elsif ($annotType =~ /prot:PA/) {
			print "<TR><TD>&nbsp;</TD></TR>\n<TR><TH align=\"left\" bgcolor=\"$dark\">&nbsp;Pathway $annotLabel for proteins:&nbsp;</TH></TR>\n";
			my $sthPA=$dbh->prepare("SELECT NAME FROM PATHWAY_ANALYSIS WHERE ID_PATHWAY_ANALYSIS=?");
			foreach my $annotName (sort{$annotSet{$annotType}{$a}[0] <=> $annotSet{$annotType}{$b}[0]} keys %{$annotSet{$annotType}}) {
				if ($anaType eq 'PCA') {
					my ($pathID,$reactNum)=split(",", $annotName);
					$pathID=~s/^#//;
					$sthPA->execute($pathID);
					my ($analysisName) = $sthPA->fetchrow_array;
					print "<TR bgcolor=\"$light\"><TD nowrap>&nbsp;$analysisName > $annotSet{$annotType}{$annotName}[1]&nbsp;<BR>";
				}
				else {##cluster
					my ($pathID, $pathString)=split(':=:',$annotSet{$annotType}{$annotName}[1]);
					$pathID=~s/^#//;
					$sthPA->execute($pathID);
					my ($analysisName) = $sthPA->fetchrow_array;
					print "<TR bgcolor=\"$light\"><TD nowrap>&nbsp;$annotName: $analysisName > &nbsp;<BR>";
					foreach my $termStrg (split(':@:',$pathString)) {
						my ($reactNum, $pathName)=split('@',$termStrg);
						print "&nbsp;&nbsp;&nbsp;&nbsp;&bull;$pathName<BR>\n";
					}
				}
			}
			$sthPA->finish;
		}
		elsif ($annotType eq 'prot:LIST') { #Cluster/PCA
			if ($anaType eq 'cluster') {
				print "<TR><TD>&nbsp;</TD></TR>\n<TR><TH align=\"left\" bgcolor=\"$dark\" colspan=\"2\">&nbsp;List of highlights for proteins:&nbsp;</TH></TR>\n";
				print "<TR><TH align=\"left\" bgcolor=\"$dark\">&nbsp;Name</TH><TH align=\"left\" bgcolor=\"$dark\">&nbsp;List</TH></TR>\n";
				foreach my $annotName (sort{$annotSet{$annotType}{$a}[0] <=> $annotSet{$annotType}{$b}[0]} keys %{$annotSet{$annotType}}) {
					$annotSet{$annotType}{$annotName}[1]=~s/#//g;
					my @list;
					foreach my $listID (split(':',$annotSet{$annotType}{$annotName}[1])) {
						my ($listName)=$dbh->selectrow_array("SELECT NAME FROM CATEGORY WHERE ID_CATEGORY=$listID");
						push @list,$listName;
					}
					print "<TR bgcolor=\"$light\"><TD nowrap>&nbsp;$annotName&nbsp;</TD><TD nowrap>&nbsp;".join("<br>&nbsp;",@list)."</TD></TR>\n";
					print "<TR bgcolor=\"$dark\"><TD colspan=\"2\"></TD></TR>\n";
				}
			}
			else {
				my $sthList=$dbh -> prepare("SELECT T.NAME,L.NAME FROM CATEGORY L,CLASSIFICATION T WHERE L.ID_CLASSIFICATION=T.ID_CLASSIFICATION AND L.ID_CATEGORY=?");
				print "<TR><TD>&nbsp;</TD></TR>\n<TR><TH align=\"left\" bgcolor=\"$dark\">&nbsp;List of highlights for proteins:&nbsp;</TH></TR>\n";
				foreach my $annotName (sort{$annotSet{$annotType}{$a}[0] <=> $annotSet{$annotType}{$b}[0]} keys %{$annotSet{$annotType}}) {
					my ($listID)=($annotName=~/^#(\d+)/);
					$sthList->execute($listID);
					my ($themeName,$listName)=$sthList->fetchrow_array;
					print "<TR bgcolor=\"$light\"><TD nowrap>&nbsp;$themeName > $listName&nbsp;</TD></TR>\n";
					print "<TR bgcolor=\"$dark\"><TD></TD></TR>\n";
				}
				$sthList->finish;
			}
		}
		elsif ($annotType eq 'prot:THEME') { # cluster
			my $sthTheme=$dbh -> prepare("SELECT NAME FROM CLASSIFICATION WHERE ID_CLASSIFICATION=?");
			print "<TR><TD>&nbsp;</TD></TR>\n<TR><TH align=\"left\" bgcolor=\"$dark\">&nbsp;List of annotations for proteins:&nbsp;</TH></TR>\n";
			foreach my $annotName (sort{$annotSet{$annotType}{$a}[0] <=> $annotSet{$annotType}{$b}[0]} keys %{$annotSet{$annotType}}) {
				my ($themeID)=($annotName=~/^#(\d+)/);
				next unless $themeID;
				$sthTheme->execute($themeID);
				my ($themeName)=$sthTheme->fetchrow_array;
				print "<TR bgcolor=\"$light\"><TD nowrap>&nbsp;All lists in theme: $themeName&nbsp;</TD></TR>\n";
				print "<TR bgcolor=\"$dark\"><TD></TD></TR>\n";
			}
			$sthTheme->finish;
		}
	}
    print "</TABLE>\n";
}
print qq
|</CENTER>
</BODY>
</HTML>
|;

$dbh -> disconnect;

####>Revision history<####
# 1.1.8 [ENHANCEMENT] Code update to match new behavior of &promsQuantif(v1.5.4)::extractQuantificationParameters for MaxQuant quantif (PP 21/11/19)
# 1.1.7 fix little bug (SL 31/01/18)
# 1.1.6 Compatible with peptide exploratory analysis (SL 30/11/17)
# 1.1.5 Compatible with MaxQuant non-ratio quantifs (PP 12/01/17)
# 1.1.4 Compatible with all non-design quantifs (PP 15/02/16)
# 1.1.3 Handles infinite ratios treated as missing values (PP 22/01/16)
# 1.1.2 Handles infinite ratio and missing value filters values (PP 09/11/15)
# 1.1.1 More checks for defined values in annotations (PP 26/10/15)
# 1.1.0 Displays new fold change occurence filter (PP 06/08/15)
# 1.0.9 display List from theme (SL 15/06/15)
# 1.0.8 Displays peptide filter (PP 09/04/15)
# 1.0.7 show list => display reference if exist (SL 24/03/15)
# 1.0.6 Auto-reload Summary every 5 sec if analysis is still on-going (PP 09/01/15)
# 1.0.5 correct go analysis display for parentGoAna and add display pathway annotation (SL 06/11/14)
# 1.0.4 add protein number (SL 10/09/14)
# 1.0.3 Update for proteins and quantifications annotations/highlights listing (PP 08/08/14)
# 1.0.2 PCA compatible with new prot:GO annotation format (PP 10/07/14)
# 1.0.1 Display quantif ratios for list of quantif used & other bug fixes (PP 27/06/14)
# 1.0.0 New script to display and store the results of PCA (SL 08/04/14)
