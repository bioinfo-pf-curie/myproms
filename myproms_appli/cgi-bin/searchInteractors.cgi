#!/usr/local/bin/perl -w

################################################################################
# searchInteractors.cgi       1.1.0                                            #
# Authors: P. Poullet, S.Liva (Institut Curie)	                               #
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
use LWP::UserAgent;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
#######################
####>Argument / Processing form<####
#######################
my $itemID=param('ID');
my $item=param('ITEM');
my $projectID=param('id_project');
#my $searchString=param('searchString')? param('searchString') : '';
#my $species=param('species')? param('species') : '';
my $action=(param('ACT'))? param('ACT') : '';
my $uniAC=(param('uniAC'))? param('uniAC') : '';
my $refItem=(param('searchItem'))? param('searchItem') : $item;
my ($refItemID,$itemChkStrg,$listChkStrg,$projChkStrg)=($refItem eq 'PROJECT')? ($projectID,'','','checked') : ($refItem eq 'LIST')? (param('list'),'','checked','') : ($itemID,'checked','','');
my @itemIDList=split(",",$refItemID);
$itemChkStrg='checked' if $item eq 'PROJECT';
my $MAX_PROTEINS=20;

if ($action eq 'ajaxSearchProteins') {
    &ajaxSearchProteins;
    exit;
}

if ($action eq 'ajaxSearchInteractors') {
    &ajaxSearchInteractors;
    exit;
}

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;
my ($itemName)=$dbh->selectrow_array("SELECT NAME FROM $item WHERE ID_$item=$itemID");
my $itemType=&promsMod::getItemType($item);
my $listName; # if a list has been selected
my %speciesList;
my $sthSelSpecies=$dbh->prepare("SELECT TAXONID, COMMON_NAME, SCIENTIFIC_NAME FROM SPECIES WHERE IS_REFERENCE=1");
$sthSelSpecies->execute();
while (my($taxonID, $commonName, $scientifName)=$sthSelSpecies->fetchrow_array) {
    $speciesList{$taxonID}=$scientifName."&nbsp;(".$commonName.")&nbsp;";
}
$sthSelSpecies->finish;

my %classificationList=&promsMod::getListClass($dbh,$projectID);
my %categoryList;
my $sthCat=$dbh->prepare("SELECT ID_CATEGORY,NAME FROM CATEGORY WHERE ID_CLASSIFICATION=? ORDER BY DISPLAY_POS ASC");
foreach my $classID (keys %classificationList) {
    $sthCat->execute($classID);
    while (my ($catID,$catName)=$sthCat->fetchrow_array) {
	push @{$categoryList{$classID}},[$catID,$catName];
    }
}
$sthCat->finish;
$dbh->disconnect;

#####################################
####>Starting HTML / Search form<####
#####################################
print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Display, Search Proteins Interactions</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
function changeSearch(item) {
    if (item == 'prot') {
		window.location="./searchKwProtein.cgi?ITEM=$item&ID=$itemID&id_project=$projectID";
    }
}
function sequenceView(listAnaID, protID, idType, ms) {
    var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+listAnaID+"&id_prot="+protID+"&msdata="+ms+"&id_type="+idType;
    top.openProtWindow(winLocation);
}
function ajaxSearchProteins(myForm) {
	document.getElementById('titleDIV').innerHTML='';
    document.getElementById('resultDiv').style.display = 'none';
    document.getElementById('waitDIV').style.display = '';

    var searchString=myForm.searchString.value;
    var species=myForm.species.value;
	var extendSpecies=(myForm.extendSpecies.checked)? 1 : 0;

    if (!searchString) {
		alert('ERROR: No search string typed!');
		return;
    }
    if (!species) {
		alert('ERROR: No species selected!');
		return;
    }
    var params="ACT=ajaxSearchProteins&searchString="+searchString+"&species="+species+"&extendSpecies="+extendSpecies+"&id_project=$projectID&ITEM=$item&ID=$itemID";

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
    XHR.open("POST","./searchInteractors.cgi",true);
    XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=utf-8");
    XHR.setRequestHeader("Content-length", params.length);
    XHR.setRequestHeader("Connection", "close");
    XHR.onreadystatechange=function() {
	if (XHR.readyState==4 && XHR.responseText) {
	    var resp=XHR.responseText;
	    if (resp.match('singleUniprot')) {
		var uniInfo=XHR.responseText.split("=");
		ajaxSearchInteractors(myForm,uniInfo[1],uniInfo[2]);
	    }
	    else {
		document.getElementById('waitDIV').style.display = 'none';
		document.getElementById('resultDiv').style.display = 'block';
		document.getElementById('resultDiv').innerHTML=XHR.responseText;
	    }
	}
    }
    XHR.send(params);
}
function ajaxSearchInteractors(myForm, uniAC, uniID) {
    document.getElementById('waitDIV').style.display = 'block';
    document.getElementById('resultDiv').style.display = 'none';
    var params="ACT=ajaxSearchInteractors&uniAC="+uniAC+"&id_project=$projectID&ITEM=$item&ID=$itemID";

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
    XHR.open("POST","./searchInteractors.cgi",true);
    XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=utf-8");
    XHR.setRequestHeader("Content-length", params.length);
    XHR.setRequestHeader("Connection", "close");
    XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			var title='<FONT class="title">Search Interactors for <FONT color="green">'+uniID+' ('+uniAC+')</FONT></FONT>';
			document.getElementById('waitDIV').style.display='none';
			document.getElementById('titleDIV').innerHTML=title;
			document.getElementById('resultDiv').innerHTML=XHR.responseText;
			document.getElementById('resultDiv').style.display='';
		}
    }
    XHR.send(params);
}


var XHR=null;
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
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Search for </FONT>
<SELECT name="SearchType" onchange="changeSearch(this.value)" class="title">
<OPTION value="prot"> Proteins</OPTION>
<OPTION value="inter" selected> Interactors</OPTION>
</SELECT>
<BR><BR>
<FORM method="post" name="searchProtForm">
<INPUT type="hidden" name="id_project" value="$projectID">
<INPUT type="hidden" name="ID" value="$itemID">
<INPUT type="hidden" name="ITEM" value="$item">
<TABLE border=0 bgcolor="$darkColor" cellpadding=0 cellspacing=5>
<TR>
	<TH align=right valign=top width="120px" class="title2">Search in :</TH>
	<TD class="title2"><INPUT type="radio" name="searchItem" value="$item" $itemChkStrg>$itemType <FONT style="color:#DD0000">$itemName</FONT><BR>
	<INPUT type="radio" name="searchItem" value="LIST" $listChkStrg>List :<SELECT name="list" class="title2"><OPTION value="">-= Select =-</OPTION>
|;
foreach my $classID (sort{lc($classificationList{$a}[0]) cmp lc($classificationList{$b}[0])} keys %classificationList) {
	print "<OPTGROUP label=\"-Theme: $classificationList{$classID}[0]\">\n";
	foreach my $refCat (@{$categoryList{$classID}}) {
		print '<OPTION value="',$refCat->[0],'"';
		if ($refItem eq 'LIST' && $refCat->[0]==$refItemID) {
			print ' selected';
			$listName=$refCat->[1];
		}
		print '>',$refCat->[1],"</OPTION>\n";
	}
	print "</OPTGROUP>\n";
}
print qq|</SELECT><BR>\n|;
if ($item ne 'PROJECT') {
	print "<INPUT type=\"radio\" name=\"searchItem\" value=\"PROJECT\" $projChkStrg>Whole Project";
}
print qq|
</TR>
<TR>
    <TH align="right" valign="top" width="120px" class="title2">For :</TH>
    <TH align="left"><INPUT type="text" id="searchString"  value="" size="50" placeholder="Protein or gene identifier"/><BR>
	<FONT style="font-size:11px;font-weight:normal">Match is <B>not</B> case sensitive. Non-[a-z 0-9] are ignored.</FONT></TH>
</TR>
<TR>
    <TH align="right" class="title2" valign="top">Species :</TH>
    <TH align="left">
    <SELECT ID="species" NAME="species"><OPTION VALUE="">-=Select=-</OPTION>|;
    foreach my $taxonID (sort{&promsMod::sortSmart(lc($speciesList{$a}),lc($speciesList{$b}))} keys %speciesList) {
	print qq|<OPTION value="$taxonID:$speciesList{$taxonID}">$speciesList{$taxonID}</OPTION>|;
    }
    print qq
|    </SELECT><SPAN style="display:none"><BR>
<INPUT type="checkbox" name="extendSpecies" value="1" checked>Extend to all strains</SPAN><!-- hidden for simplication -->
    </TH>
</TR>
<TR><TD colspan=2 align=center><INPUT type="button" name="search" value="Search" onclick="ajaxSearchProteins(document.searchProtForm)" style="width:100px"/></TD></TR>
</TABLE>
</FORM><BR>
<DIV id="waitDIV" style="display:none">
<BR><BR><FONT class="title3">Fetching data. Please wait...</FONT><BR>&nbsp;&nbsp;&nbsp;&nbsp;<IMG src="$promsPath{images}/scrollbarGreen.gif"><BR><BR>
</DIV>
<DIV id="titleDIV"></DIV>
<DIV id="resultDiv"></DIV>
<BR><BR><BR>
</CENTER>
</BODY>
</HTML>
|;

sub ajaxSearchProteins {
    my $searchString=param('searchString');
    $searchString=~s/['",\.\*;\(\)\[\]]//g; # ignoring these characters
    my ($taxonID,$speciesInfo)=split(':',param('species'));
	my $speciesStrg=(param('extendSpecies') && param('extendSpecies')==1)? $speciesInfo : $taxonID;

    my $params = {
		query => "$searchString AND organism:$speciesStrg", # reviewed:yes name: and/or gene:
		sort=>'score',
		format=>'tab',
		columns=>'id,entry name,protein names,length,genes,organism', # All columns values: see "Customize" option at http://www.uniprot.org/uniprot
		limit=>$MAX_PROTEINS
	};

    my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
    $agent->timeout(10);
    $agent->env_proxy;
    push @{$agent->requests_redirectable}, 'POST';
    my $response = $agent->post("http://www.uniprot.org/uniprot/",$params);

    while (my $wait = $response->header('Retry-After')) {
		#print STDERR "Waiting ($wait)...\n";
		sleep $wait;
		$response = $agent->get($response->base);
    }

    print header(-type=>'text/plain',-charset=>'utf-8');
	my @resultLines;
    if ($response->is_success) {
		my $count=0;
		my (@uniAc,@uniId);
		@resultLines=split(/\n/,$response->content);
		foreach my $line (@resultLines) {
			++$count;
			next if $count==1;
			my ($uniAC,$uniID)=split(/\t/,$line);
			push @uniAc,$uniAC;
			push @uniId,$uniID;
		}

		if (!@uniAc) {
			print qq|<FONT class="title">No matches</FONT><BR>|;
			exit;
		}

		if ((scalar @uniAc) == 1) {
			print "singleUniprot=$uniAc[0]=$uniId[0]";
			exit;
		}


		print qq
|<FONT class="title2">Select a protein from the list to view its interactors</FONT><BR>
<BR style="font-size:5px">
<TABLE  border="0" bgcolor="$darkColor" cellspacing="0">
<TR>
	<TH class="rbBorder">&nbsp;Protein identifiers&nbsp;</TH>
	<TH class="rbBorder">Gene name(s)</TH>
	<TH class="rbBorder" style="width:700px">Description</TH>
	<TH class="rbBorder">&nbsp;Length (aa)&nbsp;</TH>
	<TH class="bBorder">&nbsp;Organism&nbsp;</TH>
</TR>
|;
		my %proteinList;
		my $lineCount=0;
		foreach my $line (@resultLines) {
			$lineCount++;
			next if $lineCount==1;
			my ($uniAC,$uniID,$pDes,$length,$gStrg,$organism)=split(/\t/,$line);
			my @gNames=split(' ',$gStrg);
			$gNames[0]="<B>$gNames[0]</B>" if @gNames;
			my $genes=join(', ',@gNames);
			@{$proteinList{$uniAC}}=($uniID,$pDes,$length,$genes,$organism);
		}
		my $strgColor=$lightColor;
		foreach my $uniAC (sort{$proteinList{$a}[0] cmp $proteinList{$b}[0]} keys %proteinList) {
			my ($uniID,$pDes,$length,$genes,$organism)=@{$proteinList{$uniAC}};
			print qq
|<TR bgcolor="$strgColor" class="list">
	<TH nowrap align=left valign="top">&nbsp;<A href="javascript:ajaxSearchInteractors(document.searchProtForm,'$uniAC','$uniID')"/>$uniID ($uniAC)</A>&nbsp;</TH>
	<TD align="center" valign="top">$genes</TD>
	<TD valign="top">$pDes</TD>
	<TD align="right" valign="top">&nbsp;$length&nbsp;</TD>
	<TD align="center" valign="top">&nbsp;$organism&nbsp;</TD>
</TR>
|;
			$strgColor = ($strgColor eq $lightColor)? $darkColor : $lightColor;
		}
	print qq
|</TABLE>|;
    }
    exit;
}

sub ajaxSearchInteractors {
    my (%concatInteract, %interactorsList, %responseList);
    my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
    $agent->timeout(10);
    $agent->env_proxy;

    $responseList{"imex"} = $agent->get('http://www.ebi.ac.uk/Tools/webservices/psicquic/imex/webservices/current/search/query/'.$uniAC);
    $responseList{"intact"} = $agent->get('http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/'.$uniAC);
    $responseList{"innatedb-all"} = $agent->get('https://psicquic.all.innatedb.com/webservices/current/search/query/'.$uniAC);
    $responseList{"irefindex"} = $agent->get('http://irefindex.vib.be/webservices/current/search/query/'.$uniAC);
    $responseList{"reactome-fi"} = $agent->get('http://www.ebi.ac.uk/Tools/webservices/psicquic/reactome-fi/webservices/current/search/query/'.$uniAC);
    $responseList{"uniprot"} = $agent->get('http://www.ebi.ac.uk/Tools/webservices/psicquic/uniprot/webservices/current/search/query/'.$uniAC);

    print header; warningsToBrowser(1);
    
    for my $responseDB (keys %responseList) {
        my $response = $responseList{$responseDB};
        
        my @resultLines;
        if ($response->is_success) {
            @resultLines = split("\n",$response->content);
            #print $response->decoded_content;  # or whatever
        }
        else {
            die $response->status_line;
        }
        
        for (my $i=0; $i<scalar(@resultLines); $i++) {
            my ($uniqInteract_A, $uniqInteract_B, $alterInteract_A, $alterInteract_B, $aliasesA, $aliasesB, $interactMethod, $firstAuthor, $publi, $taxonA, $taxonB, $interactType, $databases, $interactIdent, $score) = (split("\t", $resultLines[$i]));
            
            next if($responseDB eq 'irefindex' && $uniqInteract_A =~ /complex:/);
            
            my @interA=split(":", $uniqInteract_A);
            my @interB=split(":", $uniqInteract_B);
            shift(@interA);
            shift(@interB);
            my $interA=join(":",@interA);
            my $interB=join(":",@interB);
            $interA=~s/"//g;
            $interB=~s/"//g;
            
            my $geneAliasesA = ($responseDB eq 'reactome-fi') ? $alterInteract_A : $aliasesA; 
            my $geneAliasesB = ($responseDB eq 'reactome-fi') ? $alterInteract_B : $aliasesB; 
            
            my ($uniqInter, $uniqAliases) = ($interA ne $uniAC)? ($interA, $geneAliasesA) : ($interB, $geneAliasesB);
            my ($strgMeth)=($interactMethod=~ /\((.+)\)/ );
            my ($strgType)=($interactType=~ /\((.+)\)/ );
            $uniqAliases=~ s/"//g;
            
            # Compute uniProt aliases
            my @uniprotAliases = &getXrefByDbName($uniqAliases, "uniprotkb");
            if($responseDB eq 'innatedb-all') {
                $uniqInter = $uniprotAliases[0];
                @uniprotAliases = &getXrefByDbName($uniqAliases, "hgnc|mgi");
            }
            my $strgAliase=join("<br>\n", @uniprotAliases);
            
            # Compute Author reference
            if($responseDB eq 'irefindex') {
                $firstAuthor =~ s/[-_]/ /g;
                $firstAuthor = substr($firstAuthor, 0, -2);
                $firstAuthor =~ s/ (\d+)/. \($1\)/g;
            }
            $firstAuthor = "Article" if(!$firstAuthor);
            $firstAuthor = ucfirst($firstAuthor);
            my $pubmed=join(",", &getXrefByDbName($publi,"pubmed"));
            my $strgPubMed="<A href='http://www.ncbi.nlm.nih.gov/pubmed/?term=$pubmed' target='blank'>$firstAuthor</A>";#join(",",&getXrefByDbName($publi, "pubmed"));

            next if($concatInteract{$uniqInter} || $uniqInter eq $uniAC || $uniqInter eq '');
            push @{$concatInteract{$uniqInter}}, [$strgAliase, $strgMeth, $firstAuthor, $strgPubMed, $strgType, $databases];
            $interactorsList{$uniqInter}=1;
        }
    }
    
	my $referenceStrg='<FONT class="font11" style="font-weight:bold">Search performed with <A href="http://code.google.com/p/psicquic" target="_blank">PSICQUIC</A> (<A href="http://www.nature.com/nmeth/journal/v8/n7/full/nmeth.1637.html" target="_blank">Aranda, B. et al. Nature Methods 8, 2011</A>).</FONT>';
    if (!scalar keys %interactorsList) {
		print qq
|<BR><FONT class="title2">No interactors found!</FONT><BR>
$referenceStrg
|;
		exit;
    }

    #exit;
    my (%proteinList, @queriesList);
    my $dbh=&promsConfig::dbConnect;
    #my ($itemName)=$dbh->selectrow_array("SELECT NAME FROM $item WHERE ID_$item=$itemID");
    my $itemType = &promsMod::getItemType($item);
    my $existMatches = 0;
    my $idList = join(',', @itemIDList);
    
    if ($refItem eq 'PROJECT') {
	    push(@queriesList, "SELECT P.ID_PROTEIN,ID_ANALYSIS,AP.VISIBILITY FROM PROTEIN P,ANALYSIS_PROTEIN AP WHERE P.ID_PROTEIN=AP.ID_PROTEIN AND ID_PROJECT IN ($idList) GROUP BY AP.ID_PROTEIN ORDER BY NUM_PEP DESC,VISIBILITY DESC,CONF_LEVEL DESC"); # best ana 1st
    }
    elsif ($refItem eq 'EXPERIMENT') {
	    push(@queriesList, "SELECT ID_PROTEIN,AP.ID_ANALYSIS,AP.VISIBILITY FROM ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S WHERE AP.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_SPOT IS NULL AND ID_EXPERIMENT IN ($idList) GROUP BY AP.ID_PROTEIN ORDER BY NUM_PEP DESC,VISIBILITY DESC,CONF_LEVEL DESC");
	    push(@queriesList, "SELECT ID_PROTEIN,AP.ID_ANALYSIS,AP.VISIBILITY FROM ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S,SPOT SP,GEL2D G WHERE AP.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_SPOT=SP.ID_SPOT AND SP.ID_GEL2D=G.ID_GEL2D AND G.ID_EXPERIMENT IN ($idList) ORDER BY NUM_PEP DESC,VISIBILITY DESC,CONF_LEVEL DESC");
    }
    elsif ($refItem eq 'GEL2D') {
	    push(@queriesList, "SELECT ID_PROTEIN,AP.ID_ANALYSIS,AP.VISIBILITY FROM ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S,SPOT SP WHERE AP.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_SPOT=SP.ID_SPOT AND ID_GEL2D IN ($idList) GROUP BY AP.ID_PROTEIN ORDER BY NUM_PEP DESC,VISIBILITY DESC,CONF_LEVEL DESC");
    }
    elsif ($refItem eq 'SPOT') {
	    push(@queriesList, "SELECT ID_PROTEIN,AP.ID_ANALYSIS,AP.VISIBILITY FROM ANALYSIS_PROTEIN AP,ANALYSIS A,SAMPLE S WHERE AP.ID_ANALYSIS=A.ID_ANALYSIS AND A.ID_SAMPLE=S.ID_SAMPLE AND ID_SPOT IN ($idList) GROUP BY AP.ID_PROTEIN ORDER BY NUM_PEP DESC,VISIBILITY DESC,CONF_LEVEL DESC");
    }
    elsif ($refItem eq 'SAMPLE') {
	    push(@queriesList, "SELECT ID_PROTEIN,AP.ID_ANALYSIS,AP.VISIBILITY FROM ANALYSIS_PROTEIN AP,ANALYSIS A WHERE AP.ID_ANALYSIS=A.ID_ANALYSIS AND ID_SAMPLE IN ($idList) GROUP BY AP.ID_PROTEIN ORDER BY NUM_PEP DESC,VISIBILITY DESC,CONF_LEVEL DESC");
    }
    elsif ($refItem eq 'ANALYSIS') {
	    push(@queriesList, "SELECT ID_PROTEIN,ID_ANALYSIS,VISIBILITY FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS IN ($idList)");
    }
    else { #LIST
	    push(@queriesList, "SELECT ID_PROTEIN,-1,2 FROM CATEGORY_PROTEIN WHERE ID_CATEGORY IN ($idList)");
    }

    foreach my $query (@queriesList) {
        my $sth = $dbh->prepare($query);
        $sth->execute();
        while (my ($protID, $anaID, $visible) = $sth->fetchrow_array) {
            #print "prot:$protID<br>";
            @{$proteinList{$protID}} = ($anaID, $visible);
        }
		$sth->finish;
    }

    my %uniprotInfo;
    my $codeIdent="AC";
    my $allUniAC = "'".join("','", keys %interactorsList)."'";
    my $identifierID = $dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='$codeIdent'");
    my $sthExistProt=$dbh->prepare("SELECT MI.VALUE, P.ID_PROTEIN, P.ALIAS FROM MASTERPROT_IDENTIFIER MI INNER JOIN MASTER_PROTEIN MP ON MP.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN INNER JOIN PROTEIN P ON P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN WHERE MI.ID_IDENTIFIER=$identifierID AND MI.VALUE IN ($allUniAC)");
    $sthExistProt->execute();
    while (my ($uniAC, $protID, $alias) = $sthExistProt->fetchrow_array) {
        next unless $proteinList{$protID};
        $uniprotInfo{$uniAC}{$protID}=$alias;
        $interactorsList{$uniAC}=2;
        #print "uni:$uniAC##prot=$proteinID##alias=$alias<br>\n";
    }

    print qq
|<BR style="font-size:5px">
<TABLE  border="0" cellspacing="0">
<TR bgcolor="$darkColor">
	<TH class="rbBorder">&nbsp;Interactor&nbsp;</TH>
	<TH class="rbBorder" nowrap>&nbsp;Gene name(s)&nbsp;</TH>
	<TH class="rbBorder" nowrap>&nbsp;Interaction type&nbsp;</TH>
	<TH class="rbBorder" nowrap>&nbsp;Detection method&nbsp;</TH>
	<TH class="bBorder">&nbsp;Reference&nbsp;</TH>
</TR>
|;
    my $strgColor=$lightColor;
    foreach my $interact (sort{$interactorsList{$b} <=> $interactorsList{$a} || $a cmp $b} keys %concatInteract) {
		print qq|<TR bgcolor=$strgColor class="list"><TD valign="top">|;
		if ($uniprotInfo{$interact}) {
			my @matchedProtID;
			foreach my $protID (keys %{$uniprotInfo{$interact}}) {
				my $boldString;
				if ($proteinList{$protID}) {
					my $strgAnaID=($refItem eq 'ANALYSIS')? $itemID : $proteinList{$protID}[0];
					$boldString="<A href=\"javascript:sequenceView('$strgAnaID',$protID,'valid',0)\">$uniprotInfo{$interact}{$protID}</A>";
					$boldString="<B>$boldString</B>" if $proteinList{$protID}->[1] >= 1;
				}
				else {
					$boldString="<strike>$interact</strike>";
				}
				push @matchedProtID,"&nbsp;$boldString&nbsp;";
			}
			print join('<BR>',@matchedProtID);
		}
		else {
			print qq|&nbsp;<strike>$interact</strike>&nbsp;|;
		}
		print "</TD>";
		my (%methVal, %typeVal, %aliasVal, %pubVal);
		my (@interMethod,@interType, @interAlias, $strgAliase);
		foreach my $refInfo (@{$concatInteract{$interact}}) {
			my ($aliases, $method, $author, $publi, $type, $databases)=@{$refInfo};
			$methVal{$method}="&nbsp;$method&nbsp;";
			$typeVal{$type}="&nbsp;$type&nbsp;";
			$pubVal{$publi}="&nbsp;$publi&nbsp;";
			$strgAliase = $aliases;
		}
		my $strgMethod=join("<br>\n",values %methVal);
		my $strgType=join("<br>\n",values %typeVal);
		my $strgPubli=join("<br>\n",values %pubVal);
		print qq
|<TD nowrap>$strgAliase</TD>
<TD nowrap valign="top">$strgType</TD>
<TD nowrap valign="top">$strgMethod</TD>
<TD nowrap valign="top">$strgPubli</TD>
</TR>
|;
		$strgColor=($strgColor eq $lightColor)? $darkColor : $lightColor;
    }

print qq
|<TR><TD colspan=5>&nbsp;&nbsp;<STRIKE>Interactor not found in selected Project item</STRIKE>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$referenceStrg</TD></TR>
</TABLE>
|;
    $dbh->disconnect;
    exit;
}

sub getXrefByDbName { ##from EBI
    my ($xrefs, $dbname)=@_;
    my @listXref;
    for my $xref (split(/\|/, $xrefs)) {
	#print "1:$xref<br>\n";
	my ($db, $id,$txt)=split/[:\(\)]/, $xref;
	#print "2: $db, $id; $txt<br>\n";
	if ($db =~ /$dbname/) {
	    #print "3:$id<br>\n";
	    push @listXref, "$id";
	}
    }
    return @listXref;
}
####>Revision history<####
# 1.1.0 Add more services to retrieve Interactors (VS 04/06/19)
# 1.0.1 Minor display changes (PP 29/10/14)
# 1.0.0 new script to retrieve protein interactions through web services (SL 26/09/14)
