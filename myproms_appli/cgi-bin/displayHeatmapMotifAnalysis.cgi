#!/usr/local/bin/perl -w
################################################################################
# displayHeatmapMotifAnalysis.cgi       1.0.0                                  #
# Authors: P. Poullet, S.Liva (Institut Curie)         						   #
# Contact: myproms@curie.fr                                                    #
# display heatMap from motif enrichment analysis     				           #
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
my $explorID = param('ID');
my $projectID=param('PROJECT_ID');

my ($HMname)=$dbh->selectrow_array("SELECT NAME FROM EXPLORANALYSIS WHERE ID_EXPLORANALYSIS=$explorID");

my %motifAnalysis;
my $sthSelAnalysis = $dbh->prepare("SELECT ID_PARENT_EXPLORANALYSIS, DISPLAY_POS FROM PARENT_EXPLORANALYSIS WHERE ID_EXPLORANALYSIS=$explorID");
$sthSelAnalysis->execute();
while(my ($parentExplo, $displayPos)=$sthSelAnalysis->fetchrow_array) {
	my ($explorAnaID, $name)=$dbh->selectrow_array("SELECT ID_EXPLORANALYSIS, NAME FROM EXPLORANALYSIS WHERE ID_EXPLORANALYSIS=$parentExplo");
	$motifAnalysis{$explorAnaID}=[$name, $displayPos];
}
$sthSelAnalysis->finish;
$dbh->disconnect;

my (%isMotif, %motifRatio, %motifScore, @files, @nameOrder, @nameInfo, @motifIDlist);
foreach my $motifID (sort {$motifAnalysis{$a}[1]  <=> $motifAnalysis{$b}[1] } keys %motifAnalysis) {
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
		$motifRatio{$motifAnalysis{$motifID}[0]}{$motif}{'log10'}=log10($ratio);
	}
	close(MOTIF);
	push @nameOrder, $motifAnalysis{$motifID}[0];
	push @nameInfo, "['$motifAnalysis{$motifID}[0]','$motifAnalysis{$motifID}[0]']";
}

my $strgName = join(",", sort{$a cmp $b} @nameInfo);

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

</SCRIPT>
</HEAD>
<BODY>
<CENTER>
<FONT class="title">Display Heatmap Motif Enrichment Analysis for <font color="red">$HMname</FONT></FONT><BR><BR>

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

sub log10 {
    my $n = shift;
    return log($n)/log(10);
}

####>Revision history<####
# 1.0.0 displaying heatmap motif enrichment analysis (SL 17/07/17)