#!/usr/local/bin/perl -w
################################################################################
# heatMapMotifAnalysis.cgi       1.0.0                                         #
# Authors: P. Poullet, S.Liva  (Institut Curie)         					   #
# Contact: myproms@curie.fr                                                    #
# generate heatMap from motif enrichment analysis     				           #
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
use File::Path qw(make_path remove_tree);

#######################
####>Configuration<####
#######################
#print header(-charset=>'utf-8', -'content-encoding' => 'no'); ##DEBUG
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my $dbh=&promsConfig::dbConnect;
my $expID = param('id_exp');
my $projectID=&promsMod::getProjectID($dbh,$expID,'EXPERIMENT');

my %motifAnalysis;
my $sthSelAnalysis = $dbh->prepare("SELECT ID_EXPLORANALYSIS, NAME, PARAM_LIST, FILTER_LIST FROM EXPLORANALYSIS WHERE ID_EXPERIMENT=$expID and STATUS=1 and ANA_TYPE='MOTIF'");
$sthSelAnalysis->execute();
while(my ($explorAnaID, $name, $paramList,$filterList)=$sthSelAnalysis->fetchrow_array) {
	$motifAnalysis{$explorAnaID}=[$name, $paramList, $filterList];
}
$sthSelAnalysis->finish;


if (param('save')) {

	##save in database
	print header(-'content-encoding'=>'no',-charset=>'utf-8');
	warningsToBrowser(1);

	my $HMname = param('HMname');
	my $sthInsertExplo=$dbh->prepare("INSERT INTO EXPLORANALYSIS(ID_EXPLORANALYSIS,ID_EXPERIMENT,NAME,ANA_TYPE) VALUES (?,?,?,?)");
	my $sthInsertParentExplo=$dbh-> prepare("INSERT INTO PARENT_EXPLORANALYSIS(ID_EXPLORANALYSIS,ID_PARENT_EXPLORANALYSIS,DISPLAY_POS) VALUES (?,?,?)");

	####INSERT EXPLORANALYSIS TABLE
	my $exploID;
	if (defined($HMname)) {
		($exploID)=$dbh->selectrow_array("SELECT MAX(ID_EXPLORANALYSIS) FROM EXPLORANALYSIS");
		$exploID++;
		$sthInsertExplo->execute($exploID,$expID,$HMname,'HM_MOTIF');
	}

	##INSERT PARENT_EXPLORANALYSIS TABLE
	foreach my $motifIDstrg (split(/,/,param('motifIDlist'))) {
		my ($motifID,$pos)=split(/:/,$motifIDstrg);
		$sthInsertParentExplo -> execute($exploID,$motifID,$pos);
	}
	$sthInsertExplo->finish;
	$sthInsertParentExplo->finish;
	$dbh->commit;

	print qq
|<SCRIPT LANGUAGE="JavaScript">
top.promsFrame.selectedAction = 'summary';
parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=experiment:$expID&branchID=motifanalysis:$exploID&VIEW=functAna";
</SCRIPT>
</BODY>
</HTML>
|;
	$dbh->disconnect;
	exit;
}

if (param('submitted')) {

	my (%isMotif, %motifRatio, %motifScore, @files, @nameOrder, @nameInfo, @motifIDlist);
	my $pos=0;
	foreach my $motifID (param('heatMapMotif')) {
		$pos++;
		push @motifIDlist, $motifID.":".$pos;
		my $motifFilePath="$promsPath{explorAna}/project_$projectID/$motifID/";#motifResult.txt";
		push @files, $motifAnalysis{$motifID}[0];
		open(MOTIF,"$motifFilePath/motifResult.txt");
			while(my $line=<MOTIF>) {
			next if $.==1;
			chomp $line;
			my ($index, $motif, $score, $foregroundMatch, $foregroundSize, $backgroundMatch, $backgroundSize,$ratio)=split(/\t/, $line);
			$motif=~s/\./-/g;
			$motif=~s/"//g;
			$isMotif{$motif}=1;
			#$motifScore{$motifAnalysis{$motifID}[0]}{$motif}{'log10'}=log10($score);
			$motifRatio{$motifAnalysis{$motifID}[0]}{$motif}{'log10'}=log10($ratio);
			#$motifScore{$motifAnalysis{$motifID}[0]}{$motif}{'log2'}=log2($score);
			#$motifRatio{$motifAnalysis{$motifID}[0]}{$motif}{'log2'}=log2($ratio);
		}
		close(MOTIF);
		push @nameOrder, $motifAnalysis{$motifID}[0];
		push @nameInfo, "['$motifAnalysis{$motifID}[0]','$motifAnalysis{$motifID}[0]']";
	}

	my $strgName = join(",", sort{$a cmp $b} @nameInfo);
	my $strgListMotif = join(",", @motifIDlist);

####	my $pathScoreLog10 = "$promsPath{tmp}/motifX/Score_log10.txt";
	#my $pathRatioLog10 = "$promsPath{tmp}/motifX/Ratio_log10.txt";
####	my $pathScoreLog2 = "$promsPath{tmp}/motifX/Score_log2.txt";
####	my $pathRatioLog2 = "$promsPath{tmp}/motifX/Ratio_log2.txt";
####
####	open(MOTIF_SCORE_LOG10,">$pathScoreLog10");
	#open(MOTIF_RATIO_LOG10,">$pathRatioLog10");
####	open(MOTIF_SCORE_LOG2,">$pathScoreLog2");
####	open(MOTIF_RATIO_LOG2,">$pathRatioLog2");
####
####	print MOTIF_SCORE_LOG10 join("\t","Motif",sort{$a cmp $b} @files)."\n";
	#print MOTIF_RATIO_LOG10 join("\t","Motif",sort{$a cmp $b} @files)."\n";
####	print MOTIF_SCORE_LOG2 join("\t","Motif",sort{$a cmp $b} @files)."\n";
####	print MOTIF_RATIO_LOG2 join("\t","Motif",sort{$a cmp $b} @files)."\n";
####
	#foreach my $motif (keys %isMotif) {
####
####		print MOTIF_SCORE_LOG10 $motif;
		#print MOTIF_RATIO_LOG10 $motif;
####		print MOTIF_SCORE_LOG2 $motif;
####		print MOTIF_RATIO_LOG2 $motif;
####
		#foreach my $file (sort{$a cmp $b} @files ) {
####			my $scoreLog10 = ($motifScore{$file}{$motif}{'log10'})? $motifScore{$file}{$motif}{'log10'} : "0";
			#my $ratioLog10 = ($motifRatio{$file}{$motif}{'log10'})? $motifRatio{$file}{$motif}{'log10'} : "0";
####			my $scoreLog2 = ($motifScore{$file}{$motif}{'log2'})? $motifScore{$file}{$motif}{'log2'} : "0";
####			my $ratioLog2 = ($motifRatio{$file}{$motif}{'log2'})? $motifRatio{$file}{$motif}{'log2'} : "0";
####
####			#$score=~s/,/./g;
####
#        print MOTIF_SCORE_LOG10 "\t$scoreLog10";
			#print MOTIF_RATIO_LOG10 "\t$ratioLog10";
####			print MOTIF_SCORE_LOG2 "\t$scoreLog2";
####			print MOTIF_RATIO_LOG2 "\t$ratioLog2";
####
	#}
####	print MOTIF_SCORE_LOG10 "\n";
	#print MOTIF_RATIO_LOG10 "\n";
####	print MOTIF_SCORE_LOG2 "\n";
####	print MOTIF_RATIO_LOG2 "\n";
####
#}
####close(MOTIF_SCORE_LOG10);
####close(MOTIF_RATIO_LOG10);
####close(MOTIF_SCORE_LOG2);
####close(MOTIF_RATIO_LOG2);

#######################
####>Starting HTML<####
#######################
print header(-'content-encoding'=>'no',-charset=>'utf-8');
warningsToBrowser(1);

print qq
|<HTML>
<HEAD>
<TITLE>Heatmap Motif Analysis</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.TD {font-weight:normal;}
.TH {font-weight:bold;}
.row_0{background-color:$darkColor;}
.row_1{background-color:$lightColor;}
.annotation{width:300px;}
.popup {z-index:999;background-color:#FFFFFF;border:solid 3px #999999;padding:5px;position:absolute;display:none;box-shadow:10px 10px 20px #808080;}
</STYLE>
<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/chartLibrary2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/heatMap.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/peptidePlot.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/canvg.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/svgfix-0.2.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/other/rgbcolor.js"></SCRIPT>
<SCRIPT language="Javascript">|;
&promsMod::popupInfo();
print qq|
var HM;

function drawClustering() {
	HM = new heatMap({
		div:'heatMapDIV',
		editable:true,
		entities:{row:'Motif',column:'Name'},
		normalization:{scope:'all', limitValues:{ref:0}, reference:'user'},
		exportAsImage:['Export as image','Clustering','./exportSVG.cgi']
	     });
	HM.setColumns([$strgName]);
|;
	foreach my $motif (keys %isMotif) {
		my @values;
		foreach my $file (sort{$a cmp $b} @files ) {
			if ($motifRatio{$file}{$motif}{'log10'}) {
				push @values, $motifRatio{$file}{$motif}{'log10'};
			}
			else {
				push @values, 'null';
			}
		}
		print "\tHM.addRow(['$motif','$motif','$motif'],[",join(',',@values),"]";
		print ");\n";
	}

print qq|
	HM.draw();
}

function saveForm (myForm) {
	if (!myForm.HMname.value) {
		alert('Provide a name for Heatmap Motif Analysis');
		return false;
	}
	else {
		myForm.submit();
	}
}

</SCRIPT>
</HEAD>
<BODY>
<CENTER>
<FONT class="title">Display Heatmap Motif Enrichment Analysis</FONT><BR><BR>
<FORM name="heatMapMotifForm" method="POST" accept-charset="UTF-8">
<INPUT TYPE="hidden" name="id_exp" value="$expID">
<INPUT type="hidden" name="save" value="1"/>
<INPUT TYPE="hidden" name="motifIDlist" value=$strgListMotif>
<b>Heatmap name </b><INPUT TYPE="text" name="HMname" size="25" maxlength="100" value=""/>
<INPUT type="submit" name="save" value="Save" onclick="javascript:saveForm(document.heatMapMotifForm)"><br><br>
<DIV id="heatMapDIV" style="float:left;background-color:#FFFFFF;border:solid 3px $darkColor;padding:2px;"></DIV>
<DIV id="divDescription" class="clDescriptionCont"></DIV>
<SCRIPT type="text/javascript">
setPopup();

/********* DRAWING CLUSTERING *********/
drawClustering();
//document.getElementById('waitDIV').style.display='none';
//document.getElementById('resultDIV').style.visibility='visible';

</SCRIPT>
</FORM>
</BODY>
</HTML>
|;
	$dbh->disconnect;
	exit;

}

#check if exist motif enrichment results
my $pathToFile="$promsPath{explorAna}/project_$projectID";
my $disabButton=" ";
foreach my $exploID (keys %motifAnalysis) {
	opendir(my $dh, "$pathToFile/$exploID") || die "Can't opendir $pathToFile/$exploID: $!";
	my @logoMotif = grep { /^logoMotif/  } readdir($dh);
	if (!scalar(@logoMotif)) {$disabButton=" disabled";last;}
	close $dh;
}

###########################
###### Starting HTML ######
###########################
print header(-'content-encoding'=>'no',-charset=>'UTF-8');
warningsToBrowser(1);

print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">|;
&promsMod::popupInfo();

print qq|
var paramList = new Array();
var useMotif = {};

function addRemoveMotif(myForm, action) {

	var allMotif = myForm.allMotif;
	var heatMapMotif = myForm.heatMapMotif;
	var selected = false;
	var moved = false;
	paramList.length = 0;

	if (action == 'add') {
		for (var i = 0; i < allMotif.length; i++) {
			if (allMotif.options[i].selected) {
				useMotif[allMotif.options[i].value]=allMotif.options[i].value;
				selected = true;
				var paramTXT = allMotif.options[i].text;
				var paramID = allMotif.options[i].value;
				var paramConcat = paramID+"::"+paramTXT;
				paramList.push(paramConcat);
			}
		}
		if (selected) {
		    NEW_OPTION: for (var i = 0; i < paramList.length; i++) {
				    var listParam = paramList[i].split('\::');
				    for (var j = 0; j < heatMapMotif.options.length; j++) {
					if (listParam[0] == heatMapMotif.options[j].value) {
					    continue NEW_OPTION;
					}
				    }
				    heatMapMotif.options[heatMapMotif.options.length] = new Option(listParam[1],listParam[0]);
				}
		}
		else {
			alert('No param selected');
		}

		for (var i = 0; i < allMotif.length; i++) {
                    if (useMotif[allMotif.options[i].value]) {
                        allMotif.options[i].disabled = true;
                    }
                    else {
                        allMotif.options[i].disabled = false;
                    }
                    allMotif.options[i].selected = false;
                }


	}
	else {
		for (var i = 0; i < heatMapMotif.length; i++) {
		    if (heatMapMotif.options[i].selected) {
			selected = true;
		    }
		}
		if (selected) {
		    for (var i = 0; i < heatMapMotif.length; i++) {
			if (heatMapMotif.options[i].selected) {
			    delete useMotif[heatMapMotif.options[i].value];
			    var paramTXT = heatMapMotif.options[i].text;
			    var paramID = heatMapMotif.options[i].value;
			    heatMapMotif.removeChild(heatMapMotif.options[i]);
			    i--;
			}
		    }
		}
		else {
		    alert('No param selected');
		}
		for (var i = 0; i < allMotif.length; i++) {
		    if (useMotif[allMotif[i].value]) {
			allMotif.options[i].disabled = true;
		    }
		    else {
			allMotif.options[i].disabled = false;
		    }
		}
	}
}

function checkForm (myForm) {

	var heatMapMotif = myForm.heatMapMotif;
	var nbMot=0;
	for (var i = 0; i < heatMapMotif.length; i++) {
		heatMapMotif[i].selected=true;
		nbMot++;
	}

	if (nbMot >= 2) {
		myForm.submit();
	}
	else if (nbMot < 2) {
		alert('You have to choose at least 2 analyses for the heatmap ');
		return false;
	}
	else {
		alert('No motif analysis in the field');
		return false;
	}

}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER><FONT class="title">HeatMap Motif Enrichment Analysis</FONT><BR><BR><BR>
<FORM name="heatMapMotifForm" method="POST" accept-charset="UTF-8"  onsubmit="return(checkForm(this));">
<INPUT TYPE="hidden" name="id_exp" value="$expID">
<INPUT type="hidden" name="submitted" value="1"/>

<TABLE class="form" border=0 bgcolor="$darkColor">
    <TR>
        <TH class="gridH" width=300><FONT class="title3">Motif Enrichment Analysis </FONT></TH>
        <TH width=100>&nbsp</TH>
        <TH class="gridH" width=300><FONT class="title3">Motif Enrichment 4 heatMap</FONT></TH>
    </TR>
    <TR>
        <TD colspan=3>
            <TABLE>
                <TR>
                    <TD>
                        <SELECT multiple name="allMotif" style="width:400px;height:200px;font-weight:bold">|;
							foreach my $motifID (sort {$a cmp $b} keys %motifAnalysis) {
								print qq|<OPTION value=$motifID>@{$motifAnalysis{$motifID}}[0]</OPTION>\n|;
							}
							print qq|
                        </SELECT>
                    </TD>

                    <TD valign=middle>
						<INPUT type="button" name="add" value="> Add >" style="width:100px;font-weight:bold" onclick="addRemoveMotif(document.heatMapMotifForm,'add')"/><BR>
						<INPUT type="button" name="remove" value="< Remove <"  style="width:100px;font-weight:bold" onclick="addRemoveMotif(document.heatMapMotifForm,'remove')"/>
					</TD>
					<TD>
                        <SELECT multiple name="heatMapMotif" id="heatMapMotif" style="width:400px;height:200px;font-weight:bold">
                        </SELECT>
                    </TD>
                </TR>
            </TABLE>
        </TD>
    </TR>
    <TR>
	<TH colspan="3" class="center"><INPUT type="submit" value="Run HeapMap" onclick="javascript:checkForm(document.heatMapMotifForm)"$disabButton></TH>
    </TR>

</TABLE>

<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup();</SCRIPT>
</CENTER>
</BODY>
</HTML>|;

sub log2 {
    my $n = shift;
    return log($n)/log(2);
}

sub log10 {
    my $n = shift;
    return log($n)/log(10);
}
$dbh->disconnect;

####>Revision history<####
# 1.0.0 New script for creating, displaying and saving heatmap motif enrichment analysis (SL 04/07/17)