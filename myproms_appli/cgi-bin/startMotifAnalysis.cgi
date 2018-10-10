#!/usr/local/bin/perl -w
################################################################################
# startMotifAnalysis.cgi       1.0.2                                           #
# Authors: P. Poullet, S.Liva  (Institut Curie)         					   #
# Contact: myproms@curie.fr                                                    #
# Fetch and provide proteins and parameters for Motif enrichment analysis      #
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
my $ajax = param('AJAX') || "";
my $idModification=param('modif') || "";

if ($ajax eq 'displayListModif') {
	&ajaxDisplayListModif($idModification);
	exit;
}

if (param('submitted')) {

	print header(-charset=>'utf-8', -'content-encoding' => 'no');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR><BR><BR><BR><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">
<BR><FONT class="title3">Fetching data. Please wait...<BR>|;

	my $centralRes=param('centralRes');#central residu name
	my $width=param('width');#length of sequence
	my $indexPos=($width+1)/2;#central residue position
	my $occurence=param('occurence');
	my $motifName=param('MotifName');
	my $pvalCutoff=param('significance');
	my $bgValue=param('background');
	my $randomNbSeq=param('seq');
	my $pvalue=param('pValue');
	my $quantification= param('quantif') || "";#foreground or background selection
	my $foldChange=param('foldChange') || 2; $foldChange=1 if ($foldChange=~/[^\d\.]/ || $foldChange < 1);
	my $foldChangeType=param('foldChgType');
	my $chkInfRatio=param('infiniteRatio') || "";
	my $catID=param('listMenu') || undef;
	my $strgForeground=param('optionMenu');## quantification or list

	my $strgBackground=($bgValue eq "random")? "$bgValue:$randomNbSeq" : "$bgValue:$quantification";
	my $strgInfRatio=($chkInfRatio)? "Y" : "N";
	my $strgParamList=($strgForeground eq "quantification")? "datasource=quantif//" : "datasource=list//";
	$strgParamList.="centralRes=$centralRes//occurence=$occurence//width=$width//background=$strgBackground//significance=$pvalCutoff";
	my $strgFilterList=($strgForeground eq "quantification") ? "pValue=$pvalue//foldChange=$foldChange//foldChangeType=$foldChangeType//InfRatio=$strgInfRatio" : "" ;

	####INSERT database
	my $sthInsertMotif=$dbh->prepare("INSERT INTO EXPLORANALYSIS(ID_EXPLORANALYSIS,ID_CATEGORY,ID_EXPERIMENT,NAME,ANA_TYPE,PARAM_LIST,FILTER_LIST,CAT_EXCLUSION,STATUS,UPDATE_DATE,UPDATE_USER) VALUES (?,?,?,?,?,?,?,?,-1,NOW(),?)");
	my $sthInsertExplorQuantif=$dbh->prepare("INSERT INTO EXPLORANA_QUANTIF(ID_QUANTIFICATION, ID_EXPLORANALYSIS, TARGET_POS) VALUES (?,?,?)");
	my $motifID;
	if (defined $motifName) {
		($motifID)=$dbh->selectrow_array("SELECT MAX(ID_EXPLORANALYSIS) FROM EXPLORANALYSIS");
		$motifID++;
		$sthInsertMotif->execute($motifID,$catID,$expID,$motifName,'MOTIF',$strgParamList,$strgFilterList,undef,$userID);
		if ($strgForeground eq "quantification") {
			my ($quantifID,$targetPos)=split(/_/,$quantification);
			$sthInsertExplorQuantif->execute($quantifID,$motifID,$targetPos);
			$sthInsertExplorQuantif->finish;
		}
		$sthInsertMotif->finish;
		$dbh->commit;
	}

	my (@selectedQuantifications, %quantifValues, %quantifInfo, %proteinInfo, %proteinForeground, %proteinBackground, %proteinSequence, $minRatio, $maxRatio, $nbTotAA);
	my %compAA=("G"=>1,"A"=>1,"S"=>1,"P"=>1,"V"=>1,"T"=>1,"C"=>1,"L"=>1,"I"=>1,"N"=>1,"D"=>1,"Q"=>1,"K"=>1,"E"=>1,"M"=>1,"H"=>1,"F"=>1,"R"=>1,"Y"=>1,"W"=>1);#for random sequence
	my %isoforms=();

	if ($catID) {#list here
		&fetchProtMod($catID,$idModification,\%isoforms);##fetch and reconstruct isoforms from modification_site
		foreach my $protID (keys %isoforms) {
			my ($seq)=&fetchSequence($protID);## sequence
			$proteinSequence{$protID}=$seq;
			$proteinForeground{$protID}{$isoforms{$protID}}="";
			if ($bgValue eq "random") {##for future calculation of probability in proteome AA composition
				my @seqLength=split(//,$seq);
				for (my $i=0;$i<$#seqLength;$i++) {
					if ($compAA{$seqLength[$i]}) {
						$compAA{$seqLength[$i]}++;
						$nbTotAA++;
					}
				}
			}
		}
	}

	if ($bgValue eq "quanti" || ($quantification && $bgValue eq "random")) {##for foreground or background quantification
		push @selectedQuantifications, $quantification;
		my %parameters=(QUANTIF_FAMILY=>'RATIO',VIEW=>'list',NUM_PEP_CODE=>'NUM_PEP_USED',QUANTIF_LIST=>\@selectedQuantifications, VERBOSE=>'Y');
		&promsQuantif::fetchQuantificationData($dbh,\%parameters,\%quantifInfo,\%quantifValues,\%proteinInfo);
		$minRatio=$parameters{MIN_RATIO};
		$maxRatio=$parameters{MAX_RATIO};

		foreach my $protModID (sort{$a cmp $b} keys %quantifValues) {
			my ($protID,$mod)=split(/-/,$protModID);
			my ($protSeq)=&fetchSequence($protID);
			$proteinSequence{$protID}=$protSeq;
			if ($bgValue eq "quanti") {##Background quantification
				$proteinBackground{$protID}{$mod}="";
			}
			else {##RANDOM calculate AA frequency to build a random sequence
				my @seqLength=split(//,$protSeq);
				for (my $i=0;$i<$#seqLength;$i++) {
					if ($compAA{$seqLength[$i]}) {
						$compAA{$seqLength[$i]}++;
						$nbTotAA++;
					}
				}
			}

			if (!$catID) {##FILTERING PART
				if ($chkInfRatio) {next if ($quantifValues{$protModID}{$quantification}{'RATIO'} == 1000 || $quantifValues{$protModID}{$quantification}{'RATIO'} == 0.001);}
				my $protPvalue= ($quantifValues{$protModID}{$quantification}{'P_VALUE'})? $quantifValues{$protModID}{$quantification}{'P_VALUE'} : 1;
				next if ($protPvalue > $pvalue);
				if ($foldChangeType eq "abs") {next if ($quantifValues{$protModID}{$quantification}{'RATIO'} < $foldChange && $quantifValues{$protModID}{$quantification}{'RATIO'} > 1/$foldChange);}
				elsif ($foldChangeType eq "up") {next if ($quantifValues{$protModID}{$quantification}{'RATIO'} < $foldChange);}
				else {next if ($quantifValues{$protModID}{$quantification}{'RATIO'} > 1/$foldChange);}
				$proteinForeground{$protID}{$mod}=""; ##Foreground
			}
		}
	}

	###>Config<####
	mkdir "$promsPath{tmp}" unless -e "$promsPath{tmp}";
	mkdir "$promsPath{explorAna}/project_$projectID" unless -e "$promsPath{explorAna}/project_$projectID";
	mkdir "$promsPath{tmp}/motifX" unless -e "$promsPath{tmp}/motifX";
	mkdir "$promsPath{tmp}/motifX/$motifID" unless -e "$promsPath{tmp}/motifX/$motifID";

	if ($bgValue eq "random") {##calculate the probability of AA in our data to generate proteome random sequence.
		my $tmpMotifPath="$promsPath{tmp}/motifX/$motifID";
		open(COMP,">$tmpMotifPath/randomCompAA.txt");
		foreach my $char(sort{$a cmp $b} keys %compAA) {
			print COMP "$char\t".(($compAA{$char}/$nbTotAA)*100)."\n";
		}
		close(COMP);
	}

	##launch createMotifFile function, create files with sequence for foreground and background
	my $return=&createMotifFile(\%proteinBackground, \%proteinForeground, \%quantifInfo, \%proteinSequence, $width, $indexPos, $centralRes, $motifID, $bgValue,$catID);
	$dbh->disconnect;
#print qq|launchMotifEnrichmentAnalyses.pl $motifID $projectID $centralRes $occurence $pvalCutoff $bgValue $width $randomNbSeq|;
#	exit;
	if ($return) {
		my $childMotif = fork;
		unless ($childMotif) { # child here
			#>Disconnecting from server
			open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
			open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
			open STDERR, '>/dev/null' or die "Can't open /dev/null: $!";
			system "./launchMotifEnrichmentAnalyses.pl $motifID $projectID $centralRes $occurence $pvalCutoff $bgValue $width $randomNbSeq 2> $promsPath{tmp}/motifX/$motifID/error.txt";
			exit;
		}
	print qq
|<SCRIPT LANGUAGE="JavaScript">
top.promsFrame.selectedAction = 'summary';
parent.itemFrame.location="$promsPath{cgi}/openProject.cgi?ACT=experiment&EXPERIMENT=experiment:$expID&branchID=motifanalysis:$motifID&VIEW=functAna";
</SCRIPT>
</BODY>
</HTML>
|;
	}
	exit;

}

##To create select menu for modification and fetch all residues who belongs to modification
my (%allModif,%modifInfo);
my $sthSelModif=$dbh->prepare("SELECT M.ID_MODIFICATION,M.PSI_MS_NAME, M.INTERIM_NAME, M.SYNONYMES, M.SPECIFICITY FROM PROJECT_MODIFICATION PM, MODIFICATION M
								WHERE PM.ID_MODIFICATION=M.ID_MODIFICATION and ID_PROJECT=? ");
my $sthSpecificity=$dbh->prepare("select DISTINCT(SPECIFICITY) from SAMPLE S, ANALYSIS A , ANALYSIS_MODIFICATION AM
								 where S.ID_EXPERIMENT=? and S.ID_SAMPLE=A.ID_SAMPLE and A.ID_ANALYSIS = AM.ID_ANALYSIS and ID_MODIFICATION=?");
$sthSelModif->execute($projectID);
while(my ($modifID, $psiName, $interName, $synonymes, $specificity)=$sthSelModif->fetchrow_array) {
	my $strgName= ($psiName)? $psiName : ($interName)? $interName : (split(/##/,$synonymes))[0];
	my @modifRes;
	$sthSpecificity->execute($expID,$modifID);

	while (my ($spec)=$sthSpecificity->fetchrow_array) {
		foreach my $res (split/[,;]/,$spec) {
			next if ($res eq "=" || $res eq "-" || $res eq "\*");
			push @modifRes, $res;
		}
	}
	my $strgModifRes = join("",@modifRes);
	$modifInfo{$modifID}=[$strgName,$strgModifRes];
	$allModif{$modifID}=$strgName;
}
$sthSelModif->finish;
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
var modifInfo = {};|;
	foreach my $modID (sort{$a <=> $b} keys %modifInfo) {
		print "modifInfo[$modID]='$modifInfo{$modID}[1]';\n";
	}


print qq|
//AJAX
var XHR = null;
function getXMLHTTP() {
    var xhr = null;
    if (window.XMLHttpRequest) {// Firefox & others
        xhr = new XMLHttpRequest();
    }
    else if (window.ActiveXObject) { // Internet Explorer
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

function displayForegroundMenu (modID) {
	var quantifDiv=document.getElementById('displayQuantif');
	if (!modID) {
		document.getElementById('listOption').style.display= 'none';
		document.getElementById('listQuantif').style.display= 'none';
		document.getElementById('FG_FILTER').style.display= '';
		document.getElementById('optionMenu').selectedIndex = "";
		quantifDiv.innerHTML='<input type="radio" name="background" id="quant" value="quanti">&nbsp;quantification selected&nbsp;<br>'
		return;
	}
	document.getElementById('listOption').style.display ='block';

}


function ajaxDisplayData(item) {

	var quantifDiv=document.getElementById('displayQuantif');

	if (!item) {
		document.getElementById('listQuantif').style.display= 'none';
		document.getElementById('FG_FILTER').style.display= '';
		quantifDiv.innerHTML='<input type="radio" name="background" id="quant" value="quanti">&nbsp;quantification selected&nbsp;<br>'
		return;
	}
	else {
		document.getElementById('listQuantif').style.display= '';
	}

	var modID=document.getElementById('modif').value;
	createOptionResidues(modID);

	//Ajax
	//If XHR object already exists, the request is canceled & the object is deleted
	if (XHR && XHR.readyState != 0) {
		XHR.abort();
		delete XHR;
	}

	//Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if ( !XHR ) {
		return false;
	}
	var listDiv=document.getElementById('listQuantif');
	listDiv.innerHTML = '';

	if (item == "quantification") {
		if (document.getElementById('FG_FILTER').style.display == 'none') {document.getElementById('FG_FILTER').style.display='';}
		XHR.open("GET","./startExploratoryAnalysis.cgi?AJAX=displaySelect&MODIF="+modID+"&QUANTIF_FAM=RATIO&ID=$expID&CALL=motif",true);
		XHR.onreadystatechange=function() {
			if (XHR.readyState == 4 && XHR.responseText) {
					listDiv.innerHTML=XHR.responseText;
			}
		}
		XHR.send(null);
		quantifDiv.innerHTML='<input type="radio" name="background" id="quant" value="quanti">&nbsp;quantification selected&nbsp;<br>'
	}
	else {
		//if (document.getElementById('FG_FILTER').style.display = ''){document.getElementById('FG_FILTER').style.display='none';}
		document.getElementById('FG_FILTER').style.display='none';

		XHR.open("GET","./startMotifAnalysis.cgi?AJAX=displayListModif&MODIF="+modID+"&id_exp=$expID",false);
		XHR.onreadystatechange=function() {
			if (XHR.readyState == 4 && XHR.responseText) {
					listDiv.innerHTML=XHR.responseText;
			}
		}

		XHR.send(null);
		ajaxDisplayBackground(modID);

	}
}

function ajaxDisplayBackground (modID) {

	var quantifDiv=document.getElementById('displayQuantif');

	//Ajax
	//If XHR object already exists, the request is canceled & the object is deleted
	if (XHR && XHR.readyState != 0) {
		XHR.abort();
		delete XHR;
	}
	//Creation of the XMLHTTPRequest object
	XHR = getXMLHTTP();
	if ( !XHR ) {
		return false;
	}

	XHR.open("GET","./startExploratoryAnalysis.cgi?AJAX=displaySelect&MODIF="+modID+"&QUANTIF_FAM=RATIO&ID=$expID&CALL=motif",true);
	XHR.onreadystatechange=function() {
		if (XHR.readyState == 4 && XHR.responseText) {
				quantifDiv.innerHTML='<input type="radio" name="background" id="quant" value="quanti">'
				quantifDiv.innerHTML+=XHR.responseText;
		}
	}
	XHR.send(null);
}


function createOptionResidues(modID) {
	//console.log('modif='+modID);
	var centralResSEL = document.motifAnaForm.centralRes;
	centralResSEL.options.length=0;
	for (var i=0; i<modifInfo[modID].length; i++) {
		//console.log(modifInfo[modID][i]);
		centralResSEL.options[i]=new Option(modifInfo[modID][i],modifInfo[modID][i]);
	}

}

function checkForm(myForm) {

	if (!myForm.MotifName.value) {
		alert('Provide a name for Motif Analysis');
		return false;
	}

	if (myForm.optionMenu.value == 'quantification' && !myForm.quantif.value) {
		alert('Choose a foreground');
		return false;
	}

	if (myForm.optionMenu.value=='list') {
		if (myForm.quant.checked && myForm.quantif.value=="") {
			alert('Choose a background');
			return false;
		}
	}

	if (!myForm.rand.checked && !myForm.quant.checked) {
		alert('Choose a background');
		return false;
	}
}

</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER><FONT class="title">Motif Enrichment Analysis</FONT><BR><BR><BR>
<FORM name="motifAnaForm" method="POST" onsubmit="return(checkForm(this));">
<INPUT type="hidden" name="id_exp" value="$expID">
<TABLE border="0" bgcolor="$darkColor" >
	<TR>
		<TH align="right" nowrap>Motif Analysis Name :</TH>
		<TD bgcolor=$lightColor><INPUT type="text" name="MotifName" size="35" maxlength="50" value=""/>&nbsp;</TD>
	</TR>
	<TR>
		<TH align="right" nowrap>Modifications :</TH>
		<TD bgcolor=$lightColor>
		  <SELECT  name="modif" id="modif" onchange="displayForegroundMenu(this.value)">
			<OPTION value="">==Select==</OPTION>|;
			foreach my $modID (keys %allModif) {
				print qq|<OPTION value="$modID">$allModif{$modID}</OPTION>|;
			}
			print qq|
		  </SELECT>
		</TD>
	</TR>
	<TR >
		<TH align="right" nowrap>Foreground Selection :</TH>
		<TD bgcolor=$lightColor nowrap>
			<DIV id="listOption" style="display:none;float:left">
				<SELECT name="optionMenu" id="optionMenu" onchange="ajaxDisplayData(this.value)">
				<OPTION value="">-=Select=-</OPTION>
				<OPTION value="quantification">Quantification</OPTION>
				<OPTION value="list">List</OPTION>
				</SELECT>
			</DIV>
			<DIV id="listQuantif" style="overflow:hidden">
			</DIV>
		</TD>
	</TR>
	<TR id="FG_FILTER" >
		<TH valign="top" align="right" nowrap >Foreground filtering :</TH>
		<TD valign="top" nowrap bgcolor=$lightColor >
			Ratio&nbsp;
			<SELECT name="foldChgType">
				<OPTION value="abs">Up &ge; or Down &le; 1/</OPTION>
				<OPTION value="up">Up &ge;</OPTION>
				<OPTION value="down">Down &le; 1/</OPTION>
			</SELECT>

			<INPUT type="text" name="foldChange" value="2" size=2>
			&nbsp;<input type="checkbox" name="infiniteRatio" value="1">Exclude infinite ratios
			<br>
			p-value &le;&nbsp;&nbsp;<INPUT type="text" name="pValue" value="0.01" size="5">
		</TD>

	</TR>
	<TR>
		<TH align=right nowrap>Central residue :</TH>
		<TD bgcolor=$lightColor><SELECT NAME="centralRes"><!-- *Updated by JavaScript* --></SELECT></TD>
	</TR>
	<TR>
		<TH align=right nowrap>Width :</TH><TD bgcolor=$lightColor><INPUT type="number" name="width" min="1" step="2" value="11"></TD>
	</TR>
	<TR>
		<TH align=right nowrap>Min. Occurrence :</TH><TD bgcolor=$lightColor><INPUT type="number" name="occurence" min="5"  value="20"></TD>
	</TR>
	<TR>
		<TH align=right nowrap>Significance :</TH><TD bgcolor=$lightColor><INPUT type="text"   name="significance" value="0.000001" max="0.0005" ></TD>
	</TR>
	<TR>
		<TH valign=top nowrap>Background Selection :</TH>
		<TD align=left nowrap bgcolor=$lightColor>

		<DIV id="displayQuantif">
			<input type="radio" name="background" id="quant" value="quanti">&nbsp;quantification selected&nbsp;<br>
		</DIV>
		<input type="radio" name="background" id="rand" value="random">&nbsp;Random&nbsp;
		<input type="text" name="seq" value="100000">&nbsp;&nbsp;sequences
		</TD>
	</TR>
	<TR>
		<TD colspan="2" align="center"><INPUT type="submit" name="submitted" value="Run Motif Analysis"></TD>
	</TR>
</TABLE>
</FORM>
</CENTER>
</BODY>
|;

sub createMotifFile {

	my ($refBackProt, $refForeProt, $refInfo ,$refSequence, $width, $indexPos, $centralRes, $motifID, $backVal, $idCat)=@_;
	#print header(-'content-encoding'=>'no',-charset=>'UTF-8');
	#warningsToBrowser(1);
	my $tmpMotifPath="$promsPath{tmp}/motifX/$motifID";

	##create file with list of analysis_id for displayMotifAnalysis.cgi script
	open (ANA, ">$tmpMotifPath/analysisList.txt");
	if (!$idCat) {
		foreach my $quantifID (sort{ $refInfo->{$a}[5] <=> $refInfo->{$b}[5]} keys $refInfo) {
			print ANA  @{$refInfo->{$quantifID}}[5];
		}
	}
	else {
		my $sthSelAnaProt=$dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=? LIMIT 0,1");
		my @anaTab;
		foreach my $protID (keys $refForeProt) {
		$sthSelAnaProt->execute($protID);
			my ($anaID)=$sthSelAnaProt->fetchrow_array;
			push @anaTab, $anaID;
		}
		print ANA join(",",@anaTab);
	}
	close(ANA);

	if ($backVal eq "quanti") {#pour le background QUANTI
		open(MOTIF,">$tmpMotifPath/backgroundQuanti.txt");
		my %checkBackIsoforms;
		foreach my $protID (sort{$a <=> $b} keys $refBackProt) {
			my $sequence=$refSequence->{$protID};
			foreach my $modRes (keys %{$refBackProt->{$protID}}) {
				my @seqLength=split(//,$sequence);
				my $seqLength=scalar(@seqLength);
				if ($modRes=~/\~/) {
					my $position = (split(/:/,$modRes))[0];
					my ($itemStart, $itemEnd)=split(/~/,$position);
					my $lengthStrg=($itemEnd-$itemStart)+1;
					next if ( ($itemStart-1)+$lengthStrg > $seqLength);
					my $subSeq=substr($sequence,$itemStart-1,$lengthStrg);
					my $delta=0;
					foreach my $char (split//,$subSeq) {
						if ($char eq $centralRes) {
							my $pos=($itemStart)+$delta;
							if ($checkBackIsoforms{$protID}{$pos}{$char}) {
								$delta++;
								next;
							}
							if ( ($pos+$width)>$seqLength) {
								$delta++;
								next;
							}
							elsif (($pos-$width)<0) {
								$delta++;
								next;
							}
							else {
								my $ambSeq = substr($sequence,$pos-$indexPos,$width);
								print MOTIF "$ambSeq\n";
								$delta++;
								$checkBackIsoforms{$protID}{$pos}{$char}=1;
								next;
							}
						}
						$delta++;
					}
				}
				elsif ($modRes=~/\./) {
					foreach my $item (split(/\./,$modRes)) {
						my $aminoAcid = substr($item,0,1);
						my $pos = substr($item,1);
						next if (uc($aminoAcid) ne uc($centralRes));
						next if $checkBackIsoforms{$protID}{$pos}{$aminoAcid};
						if ( ($pos+$width)>$seqLength) {
							next;
						}
						elsif (($pos-$width)<0) {
							next;
						}
						else {
							my $subSeq = substr($sequence,$pos-$indexPos,$width);
							print MOTIF "$subSeq\n";
						}
						$checkBackIsoforms{$protID}{$pos}{$aminoAcid}=1;
					}
				}
				else {
					my $pos=substr($modRes,1);
					my $itemStart=$pos;
					my $itemEnd=$pos;
					my $aminoAcid = substr($modRes,0,1);
					next if (uc($aminoAcid) ne uc($centralRes));
					next if $checkBackIsoforms{$protID}{$pos}{$aminoAcid};
					if ( ($pos+$width)>$seqLength) {
						next;
					}
					elsif (($pos-$width)<0) {
						next;
					}
					else {
						my $subSeq = substr($sequence,$pos-$indexPos,$width);
						print MOTIF "$subSeq\n";
					}
					$checkBackIsoforms{$protID}{$pos}{$aminoAcid}=1;
				}
			}
		}
		close(MOTIF);
	}

	if (keys $refForeProt) {#FOREGROUND
		open(PROT_INFO,">$tmpMotifPath/foreground_protein.txt");
		open(MOTIF,">$tmpMotifPath/foreground_sequence.txt");
		my %checkIsoforms;
		foreach my $protID (sort{$a <=> $b} keys $refForeProt) {
			my $sequence=$refSequence->{$protID};
			foreach my $modRes (keys %{$refForeProt->{$protID}}) {
				my @seqLength=split(//,$sequence);
				my $seqLength=scalar(@seqLength);
				if ($modRes=~/\~/) {
					my $position = (split(/:/,$modRes))[0];
					my ($itemStart, $itemEnd)=split(/~/,$position);
					my $lengthStrg=($itemEnd-$itemStart)+1;
					next if ( ($itemStart-1)+$lengthStrg > $seqLength);
					my $subSeq=substr($sequence,$itemStart-1,$lengthStrg);
					my $delta=0;
					foreach my $char (split//,$subSeq) {
						if ($char eq $centralRes) {
							my $pos=($itemStart)+$delta;
							if ($checkIsoforms{$protID}{$pos}{$char}) {
								$delta++;
								next;
							}
							if ( ($pos+$width)>$seqLength) {
								$delta++;
								next;
							}
							elsif (($pos-$width)<0) {
								$delta++;
								next;
							}
							else {
								my $ambSeq = substr($sequence,$pos-$indexPos,$width);
								print PROT_INFO $protID."\t".$modRes."\t".$char."".$pos."\t".$ambSeq."\n";# if ($backVal ne "quanti");
								print MOTIF "$ambSeq\n";
								$delta++;
								$checkIsoforms{$protID}{$pos}{$char}=1;
								next;
							}
						}
						$delta++;
					}
				}
				elsif ($modRes=~/\./){
					foreach my $item (split(/\./,$modRes)) {
						my $aminoAcid = substr($item,0,1);
						my $pos = substr($item,1);
						next if (uc($aminoAcid) ne uc($centralRes));
						next if $checkIsoforms{$protID}{$pos}{$aminoAcid};
						if ( ($pos+$width)>$seqLength ) {
							next;
						}
						elsif (($pos-$width)<0) {
							next;
						}
						else {
							my $subSeq = substr($sequence,$pos-$indexPos,$width);
							print PROT_INFO $protID."\t".$modRes."\t".$item."\t".$subSeq."\n";# if ($backVal ne "quanti");
							print MOTIF "$subSeq\n";
						}
						$checkIsoforms{$protID}{$pos}{$aminoAcid}=1;
					}
				}
				else {
					my $pos=substr($modRes,1);
					my $itemStart=$pos;
					my $itemEnd=$pos;
					my $aminoAcid = substr($modRes,0,1);
					next if (uc($aminoAcid) ne uc($centralRes));
					next if $checkIsoforms{$protID}{$pos}{$aminoAcid};

					if ( ($pos+$width)>$seqLength) {
						next;
					}
					elsif (($pos-$width)<0) {
						next;
					}
					else {
						my $subSeq = substr($sequence,$pos-$indexPos,$width);
						print PROT_INFO $protID."\t".$modRes."\t".$subSeq."\n";
						print MOTIF "$subSeq\n";
					}
					$checkIsoforms{$protID}{$pos}{$aminoAcid}=1;
				}
			}
		}
		close(MOTIF);
		close(PROT_INFO)
	}
	return 1;
}

sub ajaxDisplayListModif {

	print header(-'content-encoding' => 'no',-charset => 'utf-8');
	warningsToBrowser(1);

	my ($modificationID) = @_;
	my (%modif,%checkModif)=();
	my $sthSelModifSite=$dbh->prepare("SELECT CL.ID_CLASSIFICATION,CL.NAME, CA.ID_CATEGORY, CA.NAME FROM CATEGORY CA, CLASSIFICATION CL
					   WHERE CL.ID_PROJECT=? and CL.ID_CLASSIFICATION=CA.ID_CLASSIFICATION
					   AND CA.LIST_TYPE ='SITE'");
	$sthSelModifSite->execute($projectID);
	while (my ($idClassification, $className, $idCategory, $catName)=$sthSelModifSite->fetchrow_array) {
		next if $checkModif{$idClassification}{$idCategory};
		push @{$modif{$idClassification."#".$className}},$idCategory."#".$catName;
		$checkModif{$idClassification}{$idCategory}=1;
	}
	print qq|<SELECT name="listMenu" id="listMenu">|;
	if (!scalar(keys %modif)) {
		print qq|<OPTION value="">No List</OPTION>|;

	}
	else {
		foreach my $classStrg (sort{$a cmp $b} keys %modif) {
			my ($classID, $className)=split(/#/,$classStrg);
			print qq|<OPTGROUP label="$className">|;
			foreach my $catStrg (@{$modif{$classStrg}}) {
				my ($catID, $catName)=split(/#/,$catStrg);
				print qq|<OPTION value="$catID">$catName</OPTION>|;
			}
		}

	}
	print qq|</SELECT>|;

}

sub fetchSequence {
	my ($protID)=@_;

	my $sthGetSequence=$dbh->prepare("SELECT ID_MASTER_PROTEIN, PROT_SEQ FROM PROTEIN WHERE ID_PROTEIN=? and ID_PROJECT=$projectID");
	my $sthGetMasterSeq=$dbh->prepare("SELECT PROT_SEQ FROM MASTER_PROTEIN WHERE ID_MASTER_PROTEIN=?");

	$sthGetSequence->execute($protID);
	my ($masterProtID, $protSeq)=$sthGetSequence->fetchrow_array;

	if ($protSeq eq "+") {
		$sthGetMasterSeq->execute($masterProtID);
		$protSeq=$sthGetMasterSeq->fetchrow_array;
	}

	$sthGetMasterSeq->finish;
	$sthGetSequence->finish;
	return $protSeq;
}

sub fetchProtMod {

	my ($catID,$modID,$refIsoforms)=@_;
	my $sthSelMod=$dbh->prepare("SELECT ID_PROTEIN, GROUP_CONCAT('[',ID_MODIFICATION,']',RESIDUE,POSITION ORDER BY ID_MODIFICATION,POSITION SEPARATOR '.') FROM CATEGORY_PROTEIN CP,MODIFICATION_SITE MS
  WHERE CP.ID_CATEGORY_PROTEIN=MS.ID_CATEGORY_PROTEIN AND CP.ID_CATEGORY=? GROUP BY CP.ID_CATEGORY_PROTEIN");
	$sthSelMod->execute($catID);
	while (my ($protID,$modStrg)=$sthSelMod->fetchrow_array) {
		if ($modStrg=~/^\[$modID\]/) { ##match the modification
			$modStrg=~s/\[$modID\]//g;
			my $strgIso="";
			if ($modStrg=~/[+-]/) { ##Ambigous isoforms
				my $nb=(split/\./,$modStrg)[0];
				my $nb1=(split//,$nb)[0];
				my $nb2=(split//,$nb)[1];
				my $start=(split/\./,$modStrg)[1];$start=~s/\-//;
				my $end=(split/\./,$modStrg)[2];$end=~s/\+//;
				$strgIso=$start."~".$end.":".$nb1."/".$nb2;
			}
			else { ##Other isoforms
				$strgIso=$modStrg;
			}
			$refIsoforms->{$protID}=$strgIso;
		}
	}
}


####>Revision history<####
# 1.0.2 change character -> residue (SL 14/03/18)
# 1.0.1 improve script, add 2 functions fetchSequence, fetchProtMod and available for protein list (24/07/17)
# 1.0.0 New script for motif enrichment analysis (SL 27/06/17)