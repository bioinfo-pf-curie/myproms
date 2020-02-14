#!/usr/local/bin/perl -w

#############################################################################
# displayGOQuantiAnalysis.cgi    1.2.3                                      #
# Authors: P. Poullet, G. Arras ,F. Yvon & S. Liva (Institut Curie)         #
# Contact: myproms@curie.fr                                                 #
# Display heatmap of quantitative GO analyses                               #
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
$|=1;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use promsQuantif;
use strict;
use goAnalysis;
use File::Path qw(rmtree); # remove_tree

#print header(-charset=>'UTF-8'); warningsToBrowser(1); # DEBUG

my %promsPath = &promsConfig::getServerInfo;
my $dbh = &promsConfig::dbConnect;
my $userID = $ENV{'REMOTE_USER'};

my $goAnaID = param('ID') or die 'Missing ID';
my $sthMaster = $dbh->prepare("SELECT G.NAME, ID_PROJECT, ASPECT, PARAM_STRG FROM GO_ANALYSIS G, EXPERIMENT E
                              WHERE G.ID_EXPERIMENT=E.ID_EXPERIMENT
                              AND ID_GOANALYSIS=?");
$sthMaster->execute($goAnaID);
my ($goName,$projectID,$aspect,$paramStrg) = $sthMaster->fetchrow_array;
$sthMaster->finish;

my $sthQG = $dbh->prepare("SELECT Q.ID_QUANTIFICATION,RATIO,Q.ID_MODIFICATION,GROUP_CONCAT(MQ.ID_MODIFICATION ORDER BY MQ.MODIF_RANK SEPARATOR ','),QUANTIF_ANNOT,ID_DESIGN
			  FROM GOANA_QUANTIFICATION GQ
              INNER JOIN QUANTIFICATION Q ON Q.ID_QUANTIFICATION=GQ.ID_QUANTIFICATION
              LEFT JOIN MULTIMODIF_QUANTIFICATION MQ ON Q.ID_QUANTIFICATION=MQ.ID_QUANTIFICATION
			  WHERE ID_GOANALYSIS=? GROUP BY Q.ID_QUANTIFICATION");
$sthQG->execute($goAnaID);
my ($quantiID, $ratioPos,$quantifModID,$multiModifStrg,$annot,$designID) = $sthQG->fetchrow_array;
my $modifQuantifStrg=$quantifModID || $multiModifStrg || '';
$sthQG->finish;

my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};

if (param('AJAX')) {
    &getProtList;
    exit;
}

my ($criteria) = ($paramStrg =~ /criteria=([^;]+)/);
my ($threshold) = ($paramStrg =~ /threshold=([^;]+)/);
$threshold = 0.01 unless defined $threshold;
my ($logRatios) = ($paramStrg =~ /logRatios=([^;]+)/);
my @thrLogRatios = ($logRatios)? split(/,/,$logRatios) : ();
my %binData;

print header(-'content-encoding'=>'none',-charset=>'UTF-8');
warningsToBrowser(1);
my ($light,$dark)=&promsConfig::getRowColors;

my $log10=log(10);
my %goTerms;
my $sthChild = $dbh->prepare("SELECT ID_GOANALYSIS, PARAM_STRG FROM GO_ANALYSIS WHERE ID_PARENT_GOANA=?");
$sthChild->execute($goAnaID);
my $numLine=0;
while (my ($childID, $paramStrg) = $sthChild->fetchrow_array){
    my ($binIdx) = ($paramStrg =~ /bin=(\d)/);
    $binData{$binIdx}{'ID'} = $childID;
    ($binData{$binIdx}{'size'}) = ($paramStrg =~ /size=([^;]+)/);
    my ($countBinProt,$countPopProt);
	open RES,"$promsPath{go_unix}/project_$projectID/$childID/results_$aspect.txt" || die "Could not open result file: $!\n";
    while(<RES>){
		if (/^### (\d+)\s(\d+)/) { # 1st line in file
		    $countBinProt = $1;
		    $countPopProt = $2;
		    next;
	    }
		my ($goID,$goName,$pval,$numProt,$totNumProt,$protListStrg,$FDR) = split(/\t/,$_);
		next if (!$goID || $goID !~ '^GO:');
		$goName = &goAnalysis::formatGOName($goName);
	    $goTerms{$goID} = $goName;
		$binData{$binIdx}{'p-value'}{$goID}=-log($pval)/$log10;
		$binData{$binIdx}{'enrichment'}{$goID}=sprintf "%.2f",($numProt/$countBinProt)/($totNumProt/$countPopProt);
		$numLine++;
    }
    close RES;
}
$sthChild->finish;
$dbh->disconnect;

my $masterDir = "$promsPath{go_unix}/project_$projectID/$goAnaID";

###>Fetching p-values from merge.txt
#my $mergeFile = "$masterDir/merge.txt";
#my $numLine=0;
#open MERGE, "$mergeFile";
#while (my $line=<MERGE>) {
#    $numLine++;
#    chomp($line);
#    my ($goID,@bin)=split("\t", $line);
#    for (my $i=0; $i<=$#bin; $i++) {
#		$binData{$i}{'p-value'}{$goID}=-log($bin[$i])/$log10;
#    }
#}
#close MERGE;

if ($numLine < 1) { # || -z $mergeFile
        print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<CENTER>
<FONT class="title">Less than two significant terms in this Q.GO enrichment.<BR>No heatmap can be drawn with such a small set.</FONT>
</CENTER>
</BODY>
</HTML>
|;
    exit;
}

my @dendroElements;
my %orderedGoids;
my @order;
open ORDER, "$masterDir/goOrder.txt" or die $!;
while(my $lineOrder=<ORDER>){
    next if $. == 1;
    chomp($lineOrder);
    my ($index, $pos, $goID)=split("\t", $lineOrder);
    push @order, $goID;
    $orderedGoids{$pos}=$goID;
}
close ORDER;

open DENDRO, "$masterDir/goDendro.txt" or die $!;
while(<DENDRO>){
    next if $. == 1;
    chomp;
    my @line = split /\t/;
    $line[1]=~ s/-(\d+)/-$orderedGoids{$1}/;
    $line[2] =~ s/-(\d+)/-$orderedGoids{$1}/;
    push @dendroElements, join(',', @line);
}
close DENDRO;
my $dendroString = join ';', @dendroElements;

my @binNames;
my $numBins=scalar keys %binData;
my $midBinPos=int($numBins/2)+1;
my $prevValue;
foreach my $bin (sort {$a <=> $b} keys %binData) {
	my $binPos=$bin+1;
	my $binName;
	 if (defined $thrLogRatios[0]) {
		my $value=&revertLog2($thrLogRatios[$bin]) unless $binPos==$numBins;
		if ($binPos==1) {$binName="\]1/∞ - $value\]";} # ≤ ≥ ∞
		elsif ($binPos==$numBins) {$binName="\[$prevValue - ∞\[";}
		elsif ($binPos == $midBinPos) {$binName="\]$prevValue - $value\[";} # 3 for 1 to 5 bins
		elsif ($binPos < $midBinPos) {$binName="\]$prevValue - $value\]";}
		elsif ($binPos > $midBinPos) {$binName="\[$prevValue - $value\[";}
		else {$binName="bin$binPos";}
		$prevValue=$value unless $binPos==$numBins;
	}
	else {$binName="bin$binPos";} # older Q-GOA
	push @binNames, "$binName ($binData{$bin}{size}\%)";
}
$goName=~s/\W/_/g; # to be included in exported PNG image name
print qq
|<HTML>
<HEAD>
<TITLE>Display Heatmap</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.TD {font-weight:normal;}
.TH {font-weight:bold;}
.LINK {cursor:pointer;}
.row_0 {background-color:$dark;}
.row_1 {background-color:$light;}
.popup {z-index:100;background-color:#FFFFFF;border:solid 3px #999999;padding:5px;position:absolute;display:none;box-shadow:10px 10px 20px #808080;}
</STYLE>
<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/heatMap.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/peptidePlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
<SCRIPT LANGUAGE="JavaScript">
|;
&promsMod::popupInfo;
&promsMod::printAjaxManageSaveProteins($projectID,\%promsPath,'document.protForm.chkProt');
print qq
|var HM;
window.onload=function() {
    HM=new heatMap({
	div:'heatMap',
	name:'HM',
	editableItems:{scope:true,color:true},// type:false,
	//colLabelOrientation:'diagonal',
	entities:{row:'GO-term',column:'Bin'},
	rowOnClick:openGODescription,
	moveColumn:false,
	singleCellSelection:true,
	cellOnMouseOver:displayCellInfo,
	cellOnClick:getProtTable,
	normalization:{scope:'row',reference:'z-score',colors:'GBR'},
	exportAsImage:['Export as image','FC_GO_$goName','./exportSVG.cgi']
    });

    HM.setColumns([['$binNames[0]',$binData{0}{'ID'}],['$binNames[1]',$binData{1}{'ID'}],['$binNames[2]',$binData{2}{'ID'}],['$binNames[3]',$binData{3}{'ID'}],['$binNames[4]',$binData{4}{'ID'}]]);
|;

    foreach my $goID (@order) {
		#foreach my $pos (sort{$a <=> $b} keys %order) {
		#my $goid = $order{$pos};
		my $rowID = $goID;
		#$rowID =~ s/GO:0*//;
		my $rowString = ''; #join ',', map { $binData{$_}{'p-value'}{$goID} } (0..4);
		foreach my $bin (0..4) {
			$rowString.=',' if $bin > 0;
			$rowString.=$binData{$bin}{'p-value'}{$goID} if $binData{$bin}{'p-value'}{$goID};
		}
		# GO~term : ncRNA 3'-end processing contain -> do not use apostrophe for $goTerms{$goid} but quotation marks instead
        print "    HM.addRow([\"$goTerms{$goID}\",'$rowID',\"$goTerms{$goID}\"],[$rowString]);\n";
    }

    print qq
|	HM.addTree('row','$dendroString');
    HM.draw();
}
const enrichFactors=[|;
	my $curBin=0;
	foreach my $binIdx (sort{$a<=>$b} keys %binData) {
		print "\n\t[";
		my $firstGoID=1;
		foreach my $goID (@order) {
			print ',' unless $firstGoID;
			print $binData{$binIdx}{'enrichment'}{$goID} if $binData{$binIdx}{'enrichment'}{$goID}; # can be undef
			$firstGoID=0;
		}
		print ']';
		print ',' if $binIdx+1 < $numBins;
	}
	print qq
|
];
function displayCellInfo(hCellValue,hCell) {
	var valueStrg;
	if (hCell.value==null) {valueStrg='* no value *';}
	else {
		valueStrg='p-value: ' + (10**-hCellValue).toExponential(2);
		if (HM.normProcess.reference=='z-score' && !hCell.isExcluded) {
			if (hCell.zscore==null) {valueStrg+=' [no z-score]';}
			else {valueStrg+=' [z='+(Math.round(100*hCell.zscore)/100)+']';}
		}
		valueStrg+='\\nEnrichment factor: '+enrichFactors[hCell.columnIdx][hCell.rowIdx];
	}
	return valueStrg;
}
function getProtTable(rowID, colID, onStatus){
    var protTable=document.getElementById("protTable");
	if (onStatus==true) {
		protTable.style.display='';
		ajaxGetProteins(rowID, colID, protTable);
	}
	else {
		protTable.style.display='none';
	}
}
function checkAllProteins(chkStatus) {
	var checkBoxList=document.protForm.chkProt;
	if (checkBoxList.length) {
		for (var i=0; i < checkBoxList.length; i++) {checkBoxList[i].checked=chkStatus;}
	}
	else {checkBoxList.checked=chkStatus;}
}
function openGODescription(rowID){
    window.open('http://amigo.geneontology.org/amigo/term/'+rowID);
}
function sequenceView(id_protein,id_analysis){
	var winLocation="$promsPath{'cgi'}/sequenceView.cgi?id_ana="+id_analysis+"&id_prot="+id_protein+"&msdata="+top.showMSdata;
	top.openProtWindow(winLocation);
}
//var XHR=null;
function ajaxGetProteins(rowID,colID,protTable) {
	//wait image
	protTable.innerHTML='<BR><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">';

	//If XHR object already exists, the request is canceled & the object is deleted
	if(XHR && XHR.readyState != 0){
		XHR.abort();
		delete XHR;
	}

	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	var URL = "$promsPath{cgi}/displayGOQuantiAnalysis.cgi?AJAX=1&ID=$goAnaID&ASPECT=$aspect&ROW="+rowID+"&BIN="+colID;
	//var URL = "$promsPath{cgi}/showProtQuantification.cgi?id_quantif=$quantiID&ACT=ajaxListProt&sort=peptide&selProt="+protData[rowID][colID];

	XHR.open("GET",URL,true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText){
			protTable.innerHTML = XHR.responseText;
			protTable.scrollIntoView({block:"start",inline:"nearest",behavior:"smooth"});
		}
	}
	XHR.send(null);
}
var isNav = (navigator.appName.indexOf("Netscape") !=-1);
function ajaxProteinData(e,protID,ratioPos,action) {
	var displayDiv=document.getElementById('displayDIV');
	var divX,divY;
	if (e) { // called from protein list
		divX = (isNav)? e.pageX : event.clientX + document.body.scrollLeft; //divX-=5;
		divY = (isNav)? e.pageY : event.clientY + document.body.scrollTop; //divY+=10;
	}
	else {
		var graphDiv=document.getElementById('volcanoDIV');
		divX = getDomOffset(graphDiv,'offsetLeft') + graphDiv.offsetWidth +2;
		divY = getDomOffset(graphDiv,'offsetTop');
	}
	displayDiv.style.left = divX + 'px';
	displayDiv.style.top = divY + 'px';
	displayDiv.style.display='block';

	var infoDIV=document.getElementById('infoDIV');
	infoDIV.innerHTML="<BR><BR><IMG src=\\"$promsPath{images}/scrollbarGreen.gif\\"><BR><BR>";

	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("GET","$promsPath{cgi}/showProtQuantification.cgi?ACT="+action+"&id_ana=0&id_quantif=$quantiID&view=list&sort=peptide&foldChange=2&pValue=1&stdDev=10&pepType=all&numPep=1&id_prot="+protID+"&ratio="+ratioPos,true); //+extraParams

	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			if (action=='ajaxPepRatio') {
				var codeParts=XHR.responseText.split('#==========#');
				infoDIV.innerHTML=codeParts[0]; // HTML part
				eval(codeParts[1]); // javascript part
			}
			else {infoDIV.innerHTML=XHR.responseText;}
		}
	}
	XHR.send(null);
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
</SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<CENTER>
<FONT class="title">Fold Change-based Gene Ontology Enrichment</FONT><BR><BR>
<IMG src="$promsPath{html}/data/go/results/project_$projectID/$goAnaID/densityPlot.png"/>
<TABLE>
<TR><TH class="title2">Heatmap of -log10(p-value)</TH></TR>
<TR><TH align="left">Click on a cell to list proteins content. Click on a GO-term for detailed information.</TD></TR>
<TR><TD><DIV id="heatMap"></DIV></TD></TR>
</TABLE>
<DIV id="protTable"></DIV>
<DIV id="displayDIV" class="popup"> <!--filter:alpha(opacity=80);opacity:0.8;-->
	<DIV id="infoDIV"></DIV>
</DIV>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT LANGUAGE="javascript">
setPopup();
</SCRIPT>
</CENTER>
</BODY>
</HTML>
|;

sub revertLog2 {
	my $l2Val=shift;
	#my $fc=exp($l2Val*log(2));
	my $fc=2**$l2Val;
	my $formatFc=($fc > 1)? sprintf "%.2f",$fc : '1/'.sprintf "%.2f",1/$fc;
	$formatFc=~s/\.0+\Z//;
	$formatFc=1 if $formatFc eq '1/1';
	return $formatFc;
}

sub cleanFiles{
    my $dir = shift;
    #remove_tree($dir);
    rmtree($dir);
}


sub getProtList{
    my $goID=param('ROW');
	my $binID=param('BIN');

	my ($anaString)=$dbh->selectrow_array("SELECT GROUP_CONCAT(ID_ANALYSIS SEPARATOR ',') FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$quantiID GROUP BY ID_QUANTIFICATION");
    my $sthBin = $dbh->prepare("SELECT PARAM_STRG FROM GO_ANALYSIS WHERE ID_GOANALYSIS=?");
    $sthBin->execute($binID);
    my ($binParamStrg) = $sthBin->fetchrow_array;
    $sthBin->finish;

    my ($bin) = ($binParamStrg =~ /bin=([^;]+)/);
    $bin++;
    my ($binSize) = ($binParamStrg =~ /size=([^;]+)/);
    my $goName;

    # Fetching proteins IDs of corresponding GO ID in result file
    my $resultDirUnix = "$promsPath{go_unix}/project_$projectID/$binID";
    my ($enrichFactor,$pvalue,$countBinProt,$countPopProt,$countGOProt);
    open(RES,"$resultDirUnix/results_$aspect.txt") || die "Could not open result file";
	my %proteinsMatched;
    while(<RES>){
	    if(/^### (\d+)\s(\d+)/){ # 1st line in file
		    $countBinProt = $1;
		    $countPopProt = $2;
		    next;
	    }
	    my @info = split(/\t/, $_);
	    if ($info[0] eq $goID) {
		    if ($info[0] eq 'unannotated'){
			    $goName = 'Unannotated proteins' ;
		    } elsif ($info[0] eq 'discarded') {
			    $goName = 'Discarded proteins';
		    } else {
			    $goName = &goAnalysis::formatGOName($info[1]);
			    $enrichFactor = sprintf("%.2f",($info[3]/$countBinProt)/($info[4]/$countPopProt));
			    $pvalue = sprintf("%1.2e",$info[2]);
		    }
		    foreach my $protID (split(';',$info[5])){
				$proteinsMatched{$protID}=1;
			    $countGOProt++;
		    }
		    last;
	    }
    }
    close RES;
    
	my $isoformName='';
    #my  @quantifModifs;
	if ($modifQuantifStrg) {
		my @quantifModifs=split(',',$modifQuantifStrg);
		my $sthMod=$dbh->prepare("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DISPLAY_CODE,DISPLAY_COLOR FROM MODIFICATION WHERE ID_MODIFICATION=?");
		#my $modifRank=0;
        foreach my $modID (@quantifModifs) {
            #$modifRank++;
            #unless ($modificationInfo{$modID}) {
                $sthMod->execute($modID);
                my ($psiName,$interName,$synName,$displayCode,$displayColor)=$sthMod->fetchrow_array;
            #    my $modifName=$psiName || $interName || $synName;
            #    $modifName=~s/^##//; $modifName=~s/##.*$//;
            #    $modificationInfo{$modID}=[$modID,$displayCode,$displayColor,$modifName];
            #}
            $isoformName.='/' if $isoformName;
            $isoformName.="<FONT color='#$displayColor'>$displayCode</FONT>";
        }
        #$isoformName=$modifName.'-sites';
		$isoformName.='-sites';
	}
    
	my $numPepCode='NUM_PEP_USED';
	my $quantif=$quantiID.'_'.$ratioPos;
	my (%quantifValues,%quantifInfo,%proteinInfo,%dispModifSites);
	if ($countGOProt) {
		my %parameters=(QUANTIF_FAMILY=>'RATIO',VIEW=>'list',NUM_PEP_CODE=>$numPepCode,QUANTIF_LIST=>[$quantif]);
		my ($numPepThr) = ($paramStrg =~ /minPep=([^;]+)/);
		$parameters{STRICT_FILTER}{$numPepCode}=$numPepThr if ($numPepThr && $numPepThr > 1);
		my ($pValueThr) = ($paramStrg =~ /ratioPvalue=([^;]+)/);
		$parameters{STRICT_FILTER}{P_VALUE}=$pValueThr if ($pValueThr && $pValueThr > 1);
		&promsQuantif::fetchQuantificationData($dbh,\%parameters,\%quantifInfo,\%quantifValues,\%dispModifSites,\%proteinInfo,\%proteinsMatched);
	}

	$dbh->disconnect;

	##>Remove PTM-ratios outside bin range:  ]bin ratio range]
	if ($modifQuantifStrg) {
		my ($logRatios) = ($paramStrg =~ /logRatios=([^;]+)/);
		my @thrLogRatios = split(/,/,$logRatios);
		my $maxRatio=2**$thrLogRatios[$bin-1] if $bin < 5;
		my $minRatio=2**$thrLogRatios[$bin-2] if $bin > 1;
		foreach my $modProtID (keys %quantifValues) {
			if ($bin < 5) {
				#delete $quantifValues{$modProtID} if $quantifValues{$modProtID}{$quantif}{RATIO} > $maxRatio; # eg. for bin#2: delete if ratio >= $thrLogRatios[1]
				if ($quantifValues{$modProtID}{$quantif}{RATIO} > $maxRatio) {
					delete $quantifValues{$modProtID};
					next;
				}
			}
			if ($bin > 1) {
				delete $quantifValues{$modProtID} if $quantifValues{$modProtID}{$quantif}{RATIO} <= $minRatio;
			}
		}
	}


    ####>HTML<####
    print header(-type=>'text/plain',-charset => 'utf-8');
    warningsToBrowser(1);

    unless ($countGOProt) {
		print qq
|<BR>
<FONT class="title1">No protein from this term in bin $bin</FONT>
<BR><BR>
|;
		exit;
    }

    print qq
|<BR>
<FONT class="title"><a href="http://amigo.geneontology.org/amigo/term/$goID" target=_blank>$goName</a></FONT><BR>
<FONT class="title3">Bin $bin ($binSize%)</FONT><BR>
|;
	print qq|<FONT class="title3">Enrichment Factor: $enrichFactor, p-value=$pvalue</FONT>| unless ($goID eq 'unannotated') or ($goID eq 'discarded');
	my $disabSave=($projectAccess eq 'guest')? ' disabled' : '';
	print qq
|<BR><BR>
<FORM name="protForm" method="POST">
<TABLE border=0 cellspacing=0 cellpadding=2>
<TR><TD colspan=8><DIV id="saveProtDIV" style="display:none;"></DIV></TD></TR>
<TR><TD colspan=8>
<INPUT type="button" value="Check all" onclick="checkAllProteins(true)"/>
<INPUT type="button" value="Uncheck all" onclick="checkAllProteins(false)"/>
<INPUT type="button" id="saveFormBUTTON" value="Save proteins..." onclick="ajaxManageSaveProteins('getThemes','PROT')"$disabSave/>
|;
	print "<INPUT type=\"button\" id=\"saveSiteFormBUTTON\" value=\"Save sites...\" onclick=\"ajaxManageSaveProteins('getThemes','SITE','$modifQuantifStrg')\"$disabSave/>\n" if $modifQuantifStrg;
	print qq
|</TD></TR>
|;
	my $itemHeaderStrg='';
	if ($modifQuantifStrg) {
		my ($totIsoforms) = ($paramStrg =~ /numSiteSel=([^;]+)/);
		$itemHeaderStrg="&nbsp;".(scalar keys %quantifValues)."/$totIsoforms $isoformName&nbsp;<BR>";
	}
	$itemHeaderStrg.="&nbsp;$countGOProt/$countPopProt proteins&nbsp;";
    print qq
|<TR class="row_0">
<TH class="rbBorder" rowspan=2>$itemHeaderStrg</TH><TH class="rbBorder" rowspan=2>&nbsp;Gene&nbsp;</TH><TH class="rbBorder" colspan=4>Quantification data</TH>
<TH class="rbBorder" rowspan=2>&nbsp;MW&nbsp;<BR><SMALL>kDa</SMALL></TH><TH class=\"bBorder\" width=700 rowspan=2>Description - Species</TH></TR>
<TR class="row_0"><TH colspan=2 class="rbBorder">&nbsp;Ratio&nbsp;</TH><TH class=\"rbBorder\">&nbsp;p-value&nbsp;</TH><TH class="rbBorder">&nbsp;Pep.&nbsp;used&nbsp;</TH></TR>
|;
    my ($minFoldChange,$maxFoldChange,$maxPvalue)=(0.5,2,0.05);
    my $numDispItems=0;
    foreach my $modProtID (sort{$quantifValues{$b}{$quantif}{$numPepCode} <=> $quantifValues{$a}{$quantif}{$numPepCode}} keys %quantifValues) {
	    $numDispItems++;
		my ($protID,$modStrg)=($modProtID=~/^(\d+)(.*)/); # $modStrg='' unless $modStrg;
        $modStrg=($modStrg)? '-'.$dispModifSites{$modProtID} : '';
		my $trClass='list row_'.($numDispItems % 2);
		print "<TR class=\"$trClass\">\n";
		# Protein/Isoform
		print "<TD class=\"TH\"><INPUT type=\"checkbox\" name=\"chkProt\" value=\"$modProtID\"/><A href=\"javascript:sequenceView($protID,'$anaString')\">$proteinInfo{$protID}[0]$modStrg</A>&nbsp;</TD>\n";
	    # gene
		print "<TH>&nbsp;";
		if (scalar @{$proteinInfo{$protID}[5]} > 1) { # gene
			print "<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-",join('<BR>&nbsp;&nbsp;-',@{$proteinInfo{$protID}[5]}[1..$#{$proteinInfo{$protID}[5]}]),"</B>')\" onmouseout=\"popout()\">$proteinInfo{$protID}[5][0]</A>";
		}
		elsif ($proteinInfo{$protID}[5][0]) {print $proteinInfo{$protID}[5][0];}
		else {print '-';} # no gene mapped
		print '</TH>';
		# Ratio
	    my $ratio=$quantifValues{$modProtID}{$quantif}{RATIO};
	    my ($absRatio,$arrowDir)=($ratio >= 1)? ($ratio,'up_') : (1/$ratio,'down_');
	    my $ratioStrg=sprintf "%.2f",$absRatio; $ratioStrg*=1;
		if ($ratio>=1000) {$ratioStrg='&infin;';}
		elsif ($ratio<=0.001) {$ratioStrg='1/&infin;';}
	    elsif ($ratio < 1) {$ratioStrg='1/'.$ratioStrg;} # if ($ratio < 1 && $ratioStrg > 1);
	    # Arrow flag
	    my $arrowStrg='';
	    my $okRatio=($ratio <= $minFoldChange || $ratio >= $maxFoldChange)? 1 : 0;
	    my $okPvalue=($quantifValues{$modProtID}{$quantif}{P_VALUE} && $quantifValues{$modProtID}{$quantif}{P_VALUE} <= $maxPvalue)? 1 : 0;
	    if ($okRatio || $okPvalue) {
		    my $arrowColor=(!$okRatio || !$okPvalue)? 'gray' : ($ratio >= 1)? 'red' : 'green';
		    $arrowStrg="<IMG class=\"LINK\" src=\"$promsPath{images}/$arrowDir$arrowColor.png\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'ajaxPepRatio')\">";
	    }
	    print "<TH nowrap align=right valign=top>&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'ajaxPepRatio')\"/>$ratioStrg</A></TH><TD valign=top>$arrowStrg&nbsp;</TD>";
	    # p-value
	    my $pValueStrg=(!defined $quantifValues{$modProtID}{$quantif}{P_VALUE})? '-' : ($quantifValues{$modProtID}{$quantif}{P_VALUE} >= 0.01)? sprintf "%.2f",$quantifValues{$modProtID}{$quantif}{P_VALUE} : sprintf "%.1e",$quantifValues{$modProtID}{$quantif}{P_VALUE};
		print "<TH nowrap valign=top>&nbsp;<A href=\"javascript:void(null)\" onclick=\"ajaxProteinData(event,'$modProtID',$ratioPos,'ajaxProtStat')\">$pValueStrg</A>&nbsp;</TH>";
		# peptides
		print "<TH valign=top>&nbsp;$quantifValues{$modProtID}{$quantif}{$numPepCode}&nbsp;</TH>";
		# mass, description, species
		my $mass=(!$proteinInfo{$protID}[1])? '-' : $proteinInfo{$protID}[1];
	    print "<TD class=\"TH\" align=right>$mass&nbsp;</TD><TD>$proteinInfo{$protID}[2] <U><I>$proteinInfo{$protID}[3]</I></U></TD>";
		print "</TR>\n";
    }
    print qq
|</TABLE></FORM><BR><BR>
|;

}

####>Revision history
# 1.2.3 [BUGFIX] Fixed incomplete compatibility with multi-modif quantifications (PP 30/11/19)
# 1.2.2 [FEATURE] Compatible with multi-modif quantifications (PP 15/07/19)
# 1.2.1 [Fix] undef bin values when no term found for a bin (PP 11/12/18)
# 1.2.0 Handles PTM-based GO-quantification, merge.txt no longer used & updated list management (PP 31/08/18)
# 1.1.3 [Fix] Now uses /usr/local/bin/perl instead of /usr/bin/perl (PP 19/09/17)
# 1.1.2 Minor graphics improvement (PP 28/01/15)
# 1.1.1 replace $userID by $masterID to skip error in path ,modify name extractOrder to goOrder and extractDendro to goDendro , move R cluster script to R_script path and rename goQuantiCluster.R (SL 02/12/14)
# 1.0.9 Add " to path for cd and $rScript execution so as to avoid problems with VPN connexions<BR>TODO: keep heatmap pictures (GA 02/12/14)
# 1.0.8 Exit when the number of significant-terms is less than 2 to avoid R errors with hclust function (GA 09/07/14)
# 1.0.7 Minor SQL query change and DIV call (GA 14/05/14)
# 1.0.6 Minor bug corrections (GA 30/04/14)
# 1.0.5 Uses rmtree instead of remove_tree (PP 10/03/14)
# 1.0.4 Bin renaming & SVG export option (PP 11/02/14)
# 1.0.3 Requires chartLibrary2.js (PP 20/09/13)
# 1.0.2 Fixed bug in numPepUsed fetching (FY 13/09/13)
# 1.0.1 Support of any quantification method with protein ratio (FY 04/04/13)
# 1.0.0 New script to display Quantification GO analysis (FY 15/03/13)
