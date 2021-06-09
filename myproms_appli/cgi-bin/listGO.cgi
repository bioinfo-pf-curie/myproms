#!/usr/local/bin/perl -w

#############################################################################
# listGO.cgi           1.0.12                                                #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                   #
# Contact: myproms@curie.fr                                                 #
# Display all proteins of current item in GO terms                          #
# from specified ontology and depth, with barPlot and tables                #
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
# TODO: make a 3-tabs page for each aspect
# TODO: Use promsMod.pm function & AJAX call to showProtQuantification.cgi to manage protein record in Lists

$|=1;
use strict;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use promsOboParser;
use GO::AnnotationProvider::AnnotationParser;
use GO::Node;
use goAnalysis;

my $DEBUG = 0;

my $itemType = uc(param('ITEM'));# or do { print header; foreach (param){print $_,'=>',param($_),"<br>"} exit};
my $itemID = param('ID');

my %promsPath = promsConfig::getServerInfo;
my $dbh = promsConfig::dbConnect;

if(param('CATEGORY')){
    manageCategory($dbh);
    exit;
}

my $projectID = promsMod::getProjectID($dbh, $itemID, $itemType);
my ($projectStatus)=$dbh->selectrow_array("SELECT STATUS FROM PROJECT WHERE ID_PROJECT=$projectID");
$projectStatus=0 unless $projectStatus;

my $userID=$ENV{'REMOTE_USER'};
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=${$userInfo[2]}{$projectID};

if(param('AJAX')){
    displayGO($dbh);
    exit;
}

# Fetching OBO files #
my %obo;
my $sthOBO = $dbh->prepare("SELECT ID_ONTOLOGY,NAME FROM ONTOLOGY WHERE STATUS=2"); # only complete OBO
$sthOBO->execute;
while(my @infos = $sthOBO->fetchrow_array){
    $obo{$infos[0]} = $infos[1];
}
$sthOBO->finish;

# Fetching GOA files #
my %goa;
my $sthGOA = $dbh->prepare("SELECT ID_GOANNOTATION,NAME,SCIENTIFIC_NAME FROM GOANNOTATION,SPECIES WHERE SPECIES.ID_SPECIES=GOANNOTATION.ID_SPECIES AND STATUS > 0");
$sthGOA->execute;
while(my @infos = $sthGOA->fetchrow_array){
    $goa{$infos[0]} = "$infos[1] ($infos[2])";
}
$sthGOA->finish;

# Aspects #
my %aspect = ( P => 'Biological Process', C => 'Cellular Component', F => 'Molecular Function');

# Depth #
my %depth = map { $_ => $_ } (2..20);

# Fetching categories and classifications
my ($itemName)=$dbh->selectrow_array("SELECT NAME FROM $itemType WHERE ID_$itemType=$itemID");
my %classificationList=&promsMod::getListClass($dbh,$projectID);
my %categoryList;
#$listMode='child' if ($classificationID && !$classificationList{$classificationID}); # in unlikely case of deletion of selected classification
my $sthCat=$dbh->prepare("SELECT ID_CATEGORY,NAME,LIST_TYPE FROM CATEGORY WHERE ID_CLASSIFICATION=? ORDER BY DISPLAY_POS ASC");
foreach my $classID (keys %classificationList) {
	$sthCat->execute($classID);
	while (my ($catID,$catName,$type)=$sthCat->fetchrow_array) {
		$type='PROT' unless $type;
		push @{$categoryList{$classID}},[$catID,$catName,$type];
	}
}
$sthCat->finish;
$dbh->disconnect;

# Starting HTML #

my ($light,$dark)=&promsConfig::getRowColors;

    # Header #
print header(-'content-encoding'=>'no',-charset=>'UTF-8');
warningsToBrowser(1);
print qq|
<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT src="$promsPath{html}/js/local/barPlot.js"></SCRIPT>
<SCRIPT language="Javascript">
|;
&promsMod::popupInfo();

####<Add to/Remove from Category options>####
if ($projectAccess ne 'guest' && $projectStatus <= 0) { # skip if no categories in project
	my $classNum=0;
	print "\n//---Classifications---\nmenu=[];\nmenu[$classNum]=[];\n";
	foreach my $classID (sort{lc($classificationList{$a}[0]) cmp lc($classificationList{$b}[0])} keys %classificationList) {
		$classNum++;
		print "menu[$classNum]=[];\n";
		my $catNum=0;
		if ($categoryList{$classID}) { # at least 1 category
			#<Check if at least 1 protein Category
			my $classOK=0;
			foreach my $refCat (@{$categoryList{$classID}}) {
				if ($refCat->[2] eq 'PROT') {
					$classOK=1;
					last;
				}
			}
			if ($classOK) {
				print "menu[$classNum][$catNum]=new Option('-=Select a custom List=-','$classID:null');\n";
				foreach my $refCat (@{$categoryList{$classID}}) {
					$catNum++;
					my $siteStrg=($refCat->[2] eq 'SITE')? ' [Sites]' : '';
					print "menu[$classNum][$catNum]=new Option('$refCat->[1]$siteStrg','$classID:$refCat->[0]');\n";
				}
				$catNum++;
			}
		}
		print "menu[$classNum][$catNum]=new Option('-=Create a new List=-','$classID:-1');\n";
	}
}
print qq|
var XHR=null;
function changeCatList(numeroMenu){
    var catList = document.getElementById('catList');

    catList.options.length=0;
    for (i=0;i<menu[numeroMenu].length;i++){
	    catList.options[i]=new Option(menu[numeroMenu][i].text,menu[numeroMenu][i].value);
    }
    catList.selectedIndex=0;
    if (menu[numeroMenu].length==1) { // no category (only create option)
	    document.getElementById('newCat').style.display='block';
    }
}
function selectClassification(){
    var classList = document.getElementById('classList');
    var addRemove = document.getElementById('addRemove');
    var catList = document.getElementById('catList');

    var classValue=classList.options[classList.selectedIndex].value;
    if (classValue == -1) { // create new classification
	    if (addRemove.value=='remove') {
		    alert('Cannot remove proteins from non-existing Theme!');
		    classList.selectedIndex=0;
		    return;
	    }
	    document.getElementById('newClass').style.display='block';
	    document.getElementById('newCat').style.display='block';
	    catList.disabled=true;
    }
    else {
	    document.getElementById('newClass').style.display='none';
	    document.getElementById('newCat').style.display='none';
	    if (classValue == "null"){
		    catList.disabled=true;
	    }
	    else{
		    catList.disabled=false;
		    window.location.href = classValue;
	    }
    }
}
function checkCreateCat(optionValue){
    var addRemove = document.getElementById('addRemove');
    var catList = document.getElementById('catList');

    var IDs=optionValue.split(":");
    if (IDs[1] && IDs[1]==-1) {
	    if (addRemove.value=='remove') {
		    alert('Cannot remove proteins from non-existing List!');
		    catList.selectedIndex=0;
		    return;
	    }
	    document.getElementById('newCat').style.display='block';
    }
    else {
	    document.getElementById('newCat').style.display='none';
    }
}
function checkClassSelection(protOpt) {
    if (protOpt != 'remove') return;

    var classList = document.getElementById('classList');
    var catList = document.getElementById('catList');

    if (catList.value==-1) {
	    alert('Cannot remove proteins from non-existing List!');
	    catList.selectedIndex=0;
	    return;
    }
    else if (classList.value==-1) {
	    alert('Cannot remove proteins from non-existing Theme!');
	    classList.selectedIndex=0;
	    return;
    }
}
function addRemFromCategory() {
    var classList = document.getElementById('classList');
    var catList = document.getElementById('catList');
    var newClassName = document.getElementById('newClassName');
    var newCatName = document.getElementById('newCatName');
    var addRemove = document.getElementById('addRemove');

    if (classList.value=='null') { // no classification selected
	    alert('You must select a Theme and a List.');
	    return;
    }
    // create new classification
    if (classList.value==-1) {
	    if (!newClassName.value \|\| newClassName.value=='New Theme name') {
		    alert('You must provide a name for the new Theme.');
		    return;
	    }
	    if (!newCatName.value \|\| newCatName.value=='New List name') {
		    alert('You must provide a name for the new List.');
		    return;
	    }
    }
    // existing classification selected
    var IDs=catList.options[catList.selectedIndex].value.split(":");
    if (classList.value != -1 && IDs[1]=='null'){ // no category selected
	    alert('You must select a List.');
	    return;
    }
    if (IDs[1]==-1 && (!newCatName.value \|\| newCatName.value=='New List name')) { // create new category
	    alert('You must provide a name for the new List.');
	    return;
    }

    if(IDs[1]=='null'){
	IDs[1] = -1;
    }
    var arrayGoBox = document.getElementsByName('checkGO');
    var stringPOST = "CATEGORY="+IDs[1]+"&ACT="+addRemove.value;

    if(classList.value == -1){
	stringPOST += "&newClassName="+newClassName.value+"&projectID=$projectID";
    } else {
	stringPOST += "&CLASS="+IDs[0];
    }
    if(IDs[1]==-1 \|\| classList.value == -1){
	stringPOST += "&newCatName="+newCatName.value;
    }

    var noCheck = true;
    for(var i=0;i<arrayGoBox.length;i++){
	var arrayProtBox = document.getElementsByName('checkProt'+arrayGoBox[i].value);
	for(var j=0;j<arrayProtBox.length;j++){
	    if(arrayProtBox[j].checked){
		noCheck = false;
		stringPOST += '&protID='+arrayProtBox[j].value;
	    }
	}
    }
    if(noCheck){
	alert('You must select at least 1 protein.');
	return;
    }

    var addProtToListDiv = document.getElementById('addProtToListDiv');
    addProtToListDiv.innerHTML='<IMG src="$promsPath{images}/scrollbarGreen.gif">';

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
    var URL = "$promsPath{cgi}/listGO.cgi";
    XHR.open("POST",URL,true);
    XHR.setRequestHeader('Content-Type','application/x-www-form-urlencoded');

    XHR.onreadystatechange=function() {
	if (XHR.readyState==4 && XHR.responseText){
	    var response = XHR.responseText;
	    addProtToListDiv.innerHTML = response;
	}
    }
    XHR.send(stringPOST);
}
function sequenceView(id_protein,id_analysis){
    var winLocation="$promsPath{cgi}/sequenceView.cgi?id_ana="+id_analysis+"&id_prot="+id_protein+"&msdata="+top.showMSdata;
    top.openProtWindow(winLocation);
}

function toggleProtDomains(searchTerm) {
    var tables = document.getElementsByClassName('GOTable');
    
    // Search in table for specific protID or geneName
    for (var i=0; i < tables.length; i++) {
        var table = tables[i];
        var id = table.id.substring(5,);
        var domain = document.getElementById("domain"+id);
        var hideUnchecked = (document.getElementById('hideUnchecked')) ? document.getElementById('hideUnchecked').checked : false;
        var found = false;
        
        if(searchTerm !== undefined && searchTerm != '' && domain.dataset.name !== undefined && domain.dataset.name.toLowerCase().includes(searchTerm.toLowerCase())) {
            found = true;
        }
        
        var lines = tables[i].children[0].children;
        var countProt = 0;
        var countProtChecked = 0;
        for(var y=0; y < lines.length; y++) {
            var line = lines[y];
            var checkbox = line.children[0].children[0];
            line.style.display = "";
            
            if(checkbox && typeof line.dataset.protid != 'undefined' && typeof line.dataset.genename != 'undefined' && typeof line.dataset.genedesc != 'undefined') {
                if(searchTerm !== undefined) {
                    checkbox.checked = (searchTerm != '' && (line.dataset.protid.toLowerCase().includes(searchTerm.toLowerCase()) \|\| line.dataset.genename.toLowerCase().includes(searchTerm.toLowerCase()) \|\| line.dataset.genedesc.toLowerCase().includes(searchTerm.toLowerCase())));
                }
                
                if(checkbox.checked) {
                    found = true;
                    countProtChecked++;
                } else if(hideUnchecked) {
                    line.style.display = "none";
                }
                
                countProt++;
            }
        }
        
        if(searchTerm !== undefined) {
            updateProtCount(id, (searchTerm != '') ? countProtChecked : undefined);
        }
        
        if(!found && hideUnchecked) {
            domain.style.display = "none";
            expandOrCollapseTable(id, "none");
        } else if(domain.style.display == "none") {
            domain.style.display = "";
        }
        
        found = false;
    }
}


function checkProt(el, goId) {
    var table = document.getElementById('table' + goId);
    var nChecked = (!table.dataset.checkedprot \|\| Number.isNaN(table.dataset.checkedprot)) ? 0 : (el.checked) ? parseInt(table.dataset.checkedprot)+1 : parseInt(table.dataset.checkedprot)-1;
    updateProtCount(goId, nChecked);
}

function updateProtCount(goId, nChecked) {
    var table = document.getElementById('table' + goId);
    var lines = table.children[0].children;
    var totalProts = lines.length-2;
    
    var protCountersEl = document.getElementsByClassName("count"+goId);
    if(protCountersEl) {
        for(var y=0; y<protCountersEl.length; y++) {
            var counterEl = protCountersEl[y];
            if(counterEl) {
                counterEl.innerHTML = (nChecked !== undefined) ? nChecked + "/" + totalProts : totalProts;
                counterEl.innerHTML += (lines.length > 1) ? " proteins" : " protein";
            }
        }
    }
    
    table.dataset.checkedprot = (nChecked !== undefined && !Number.isNaN(nChecked)) ? nChecked : 0;
}

function toggleAllTables(display=undefined) {
    var domains = document.getElementsByClassName('GODomain');
    for (var i=0; i < domains.length; i++) {
        var domain = domains[i];
        if(domain.style.display != "none") {
            var id = domain.id.substring(6,);
            expandOrCollapseTable(id, display);
        }
    }
}

function expandOrCollapseTable(id, display=undefined){
    var table = document.getElementById("table"+id);
    var image = document.getElementById("button"+id);

    if((typeof display !== 'undefined' && display != "none") \|\| (typeof display === 'undefined' && table.style.display == "none")) {
        image.src = '$promsPath{images}/minus1.gif';
        table.style.display = "";
    } else {
        image.src = '$promsPath{images}/plus.gif';
        table.style.display = "none";
    }
}

function selectByGO(goBox) {
    var goid = goBox.value;

    var allChckBox = document.getElementsByName('checkProt'+goid);
    for(var i=0;i<allChckBox.length;i++){
        allChckBox[i].checked = goBox.checked;
    }
    
    updateProtCount(goid, (goBox.checked) ? allChckBox.length : undefined);
}

function ajaxShowGO(){
    var goTable = document.getElementById("goTable");
    var canvas = document.getElementById("canvas");
    var obo = document.getElementById("obo").value;
    var goa = document.getElementById("goa").value;
    var aspect = document.getElementById("aspect").value;
    var depth = document.getElementById("depth").value;
	var restrict = document.getElementById("restrict").value;
    var minCount = (document.getElementById("minCount").value)?document.getElementById("minCount").value:0;

    // Checking values
    if(obo == 0){
	alert('Select an ontology file.');
	return;
    }
    if(goa == 0){
	alert('Select an annotation file.');
	return;
    }

    goTable.innerHTML='<FONT class="title3">Computing GO summary...</FONT><BR><IMG src="$promsPath{images}/scrollbarGreen.gif">';
    canvas.innerHTML = ''; // cleaning previous barplot

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
    var URL = "$promsPath{cgi}/listGO.cgi?AJAX=1&OBO="+obo+"&GOA="+goa+"&aspect="+aspect+"&depth="+depth+"&restrict="+restrict+"&minCount="+minCount+"&ID=$itemID&ITEM=$itemType";
    XHR.open("GET",URL,true);
    XHR.onreadystatechange=function() {
        if (XHR.readyState==4 && XHR.responseText){
            var response = XHR.responseText.split("###DATA###");
            goTable.innerHTML = response[0];
            var table = {};
            var allData = response[1].split("\\n");
            for(var i=0;i<allData.length;i++){
                if(allData[i]){
                    var data = allData[i].split("\|\|");
                    table[data[0]] = {};
                    table[data[0]]['label'] = data[1];
                    table[data[0]]['value'] = data[2] * 1;
                    table[data[0]]['color'] = data[3];
                }
            }
            var plot = new barPlot({ 'maxWidth' : 1200,
                                   'minWidth' : 600,
                                   'onBarClick' : function(id){
                                            if(document.getElementById('domain'+id).style.display != "none") {
                                                if(document.getElementById('table'+id).style.display=="none"){
                                                    expandOrCollapseTable(id);
                                                }
                                                document.getElementById(id).scrollIntoView();
                                            }
                                        }});
            plot.draw(table, canvas);
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
function getWarnings(warningButton){
    var warningDiv=document.getElementById("warnings");
    if(warningButton.value == 'Show warnings'){
	warningDiv.style.display = 'block';
	warningButton.value = 'Hide warnings';
    } else if (warningButton.value == 'Hide warnings'){
	warningDiv.style.display = 'none';
	warningButton.value = 'Show warnings';
    }
}
</SCRIPT>
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'>
<CENTER>
<FONT class="title1">Gene Ontology Summary</FONT>
<BR><BR>
<TABLE border=0 cellspacing=0 cellpadding=2 bgcolor=$dark><TR><TH colspan=2 nowrap>
|;
    # OBO selectItem #
print "&nbsp;&nbsp;";
print &selectItem("Ontology","obo", \%obo,0,"Select an ontology file");
print "&nbsp;&nbsp;&nbsp;&nbsp;";
    # GOA selectItem #
print &selectItem("Annotation","goa", \%goa,0,"Select an annotation file");
print "&nbsp;&nbsp;&nbsp;&nbsp;";
    # Aspect selectItem #
print &selectItem("Domain","aspect", \%aspect,0);
print "&nbsp;&nbsp;&nbsp;&nbsp;";
    # Depth level #
print &selectItem("Depth","depth",\%depth,1);
print "&nbsp;&nbsp;</TD></TR>\n";
	# Restrict to a Custom list (former category) #
print "<TR><TD nowrap>&nbsp;&nbsp;<FONT class=\"title3\">Restrict to proteins in List :</FONT><SELECT name=\"restrict\" id=\"restrict\"><OPTION value=\"0\">-= No restriction =-</OPTION>\n";
foreach my $classID (sort{lc($classificationList{$a}[0]) cmp lc($classificationList{$b}[0])} keys %classificationList) {
	print "<OPTGROUP label=\"-Theme: $classificationList{$classID}[0]\">\n";
	foreach my $refCat (@{$categoryList{$classID}}) {
		my $siteStrg=($refCat->[2] eq 'SITE')? ' [Sites]' : '';
		print "<OPTION value=\"",$refCat->[0],"\">",$refCat->[1].$siteStrg,"</OPTION>\n";
	}
	print "</OPTGROUP>\n";
}
print "</SELECT>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;";
    # Minimal protein amount #
print "<FONT class=\"title3\">Keep terms with at least</FONT> <INPUT type=\"text\" name=\"minCount\" id=\"minCount\" size=2 value=\"5\"> <FONT class=\"title3\">proteins</FONT>&nbsp;&nbsp;</TD>";
    # Submit button #
print "<TD align=right><INPUT type=\"button\" class=\"title3\" value=\"Display\" onclick=\"ajaxShowGO();\">&nbsp;&nbsp;&nbsp;&nbsp;</TD></TR></TABLE><BR><BR>";
print "<DIV id=\"canvas\"></DIV>";
    # GO Div #
print "<DIV id=\"goTable\"></DIV>";

    # Footer #
print qq
|</CENTER>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">
setPopup();
</SCRIPT>
</BODY>
</HTML>
|;

sub displayGO {
    my $dbh = shift;
    my $oboID = param('OBO');
    my $goaID = param('GOA');
    my $aspect = param('aspect');
    my $depth = param('depth');
	my $catID = param('restrict');
    my $minCount = param('minCount');
    my $UNMAPPED_ID = 'unmapped';
    my $OTHERS_ID = 'others';
    my @warnings;

    # Fetching analyses #
    print header(-'content-encoding'=>'no',-charset=>'UTF-8') if $DEBUG;
    my $time = time;
    my @anaID;
    if($itemType ne 'ANALYSIS'){
        my $sth;
        if($itemType eq 'SAMPLE'){
            $sth = $dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS WHERE ID_SAMPLE=?");
        } elsif ($itemType eq 'EXPERIMENT'){
            $sth = $dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS,SAMPLE
                                WHERE SAMPLE.ID_EXPERIMENT=?
                                AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE");
        } elsif ($itemType eq 'SPOT'){
            $sth = $dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS,SAMPLE
                                WHERE SAMPLE.ID_SPOT=?
                                AND ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE");
        } elsif ($itemType eq 'GEL2D'){
            $sth = $dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS,SAMPLE,SPOT
                                WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE
                                AND SAMPLE.ID_SPOT=SPOT.ID_SPOT
                                AND SPOT.ID_GEL2D=?");
        } elsif ($itemType eq 'PROJECT'){
            $sth = $dbh->prepare("SELECT ID_ANALYSIS FROM ANALYSIS,SAMPLE,EXPERIMENT
                                WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE
                                AND SAMPLE.ID_EXPERIMENT=EXPERIMENT.ID_EXPERIMENT
                                AND EXPERIMENT.ID_PROJECT=?
                                AND EXPERIMENT.ID_EXPERIMENT NOT IN (
                                    SELECT E2.ID_EXPERIMENT
                                    FROM EXPERIMENT E2
                                    INNER JOIN USER_EXPERIMENT_LOCK EU ON EU.ID_EXPERIMENT=E2.ID_EXPERIMENT
                                    WHERE E2.ID_PROJECT=$itemID AND EU.ID_USER='$userID'
                                )");
        }
        $sth->execute($itemID);
        while(my ($anaID) = $sth->fetchrow_array){
            push @anaID, $anaID;
        }
        $sth->finish;
    }
	else {
        push @anaID, $itemID;
    }

	# Fetching proteins in restriction list #
	my %catProteins;
	if ($catID) {
		my $sthRes=$dbh->prepare("SELECT DISTINCT ID_PROTEIN FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=$catID");
		$sthRes->execute;
		while (my ($protID)=$sthRes->fetchrow_array) {$catProteins{$protID}=1;}
		$sthRes->finish;
	}

    # Checking which identifier to use #
    my $sthId = $dbh->prepare("SELECT ID_IDENTIFIER FROM GOANNOTATION WHERE ID_GOANNOTATION=?"); $sthId->execute($goaID);
    my ($identifierID) = $sthId->fetchrow_array;
    $sthId->finish;

    # Fetching proteins #
    my (%prot,%masterProteins,%notInList);
    my $sthProt = $dbh->prepare("SELECT PROTEIN.ID_PROTEIN,ALIAS,ID_MASTER_PROTEIN,MW,NUM_PEP,PROT_DES,ORGANISM,PEP_COVERAGE,PEP_SPECIFICITY FROM PROTEIN,ANALYSIS_PROTEIN
                                WHERE ID_ANALYSIS=? AND ANALYSIS_PROTEIN.ID_PROTEIN=PROTEIN.ID_PROTEIN AND VISIBILITY>=1");
    foreach my $anaID (@anaID){
        $sthProt->execute($anaID);
        while(my ($protID,$identifier,$masterProtID,$mw,$numPep,$desc,$species,$pepCoverage,$pepSpecificity) = $sthProt->fetchrow_array){
			if ($catID && !$catProteins{$protID}) {
				$notInList{$protID}=1;
				next;
			}
            unless($prot{$protID} && $prot{$protID}{numPep} > $numPep){
				if($identifierID){
					my @uniprotIDs = goAnalysis::getProteinIds($dbh, $protID, $identifierID);
					$prot{$protID}{identifiers} = \@uniprotIDs;
				}
                $prot{$protID}{identifier} = $identifier;

                $prot{$protID}{mw} = sprintf "%.1f", $mw/1000;
                $prot{$protID}{numPep} = $numPep;
                $prot{$protID}{desc} = $desc;
				$prot{$protID}{species} = $species;
                $prot{$protID}{anaID} = $anaID;
                $prot{$protID}{pepCov} = ($pepCoverage)? $pepCoverage : '?';
                $prot{$protID}{pepSpe} = $pepSpecificity;
				if ($masterProtID) {
					@{$masterProteins{$masterProtID}}=();
					$prot{$protID}{masterProt} = $masterProtID;
				}
            }
        }
    }
    $sthProt->finish;

	# Fetching master proteins info #
	my ($geneNameID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
	my $sthMP=$dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN=? AND ID_IDENTIFIER=$geneNameID ORDER BY IDENT_RANK");
	foreach my $masterProtID (keys %masterProteins) {
		$sthMP->execute($masterProtID);
		while (my ($gene)=$sthMP->fetchrow_array) {
			push @{$masterProteins{$masterProtID}},$gene;
		}
	}
	$sthMP->finish;
    print "Fetch proteins: ",time-$time," sec<br>" if $DEBUG;

    # Reading OBO and GOA #
    $time = time if $DEBUG;
    my $fullOboFileName = "$promsPath{'obo'}/$oboID.obo";
    my $ontology = new promsOboParser(ontologyFile => $fullOboFileName , aspect => $aspect);
    print "Read OBO: ",time-$time," sec<br>" if $DEBUG;

    $time=time if $DEBUG;
    my $fullGoaFileName = "$promsPath{'goa'}/$goaID.goa";
    my $annotation = new GO::AnnotationProvider::AnnotationParser(annotationFile => $fullGoaFileName);
    print "Read GOA: ",time-$time," sec<br>" if $DEBUG;

    # Mapping proteins to GO terms #
    $time=time;
    my %goToProt;
    foreach my $protID (keys %prot){

		my $chosenIdentifier;
		if(!$identifierID or !scalar @{$prot{$protID}{identifiers}}){
			$chosenIdentifier = $prot{$protID}{identifier};
		}
		else {
			$chosenIdentifier = goAnalysis::selectID($prot{$protID}{identifiers}, $annotation);
		}

        my $goIdListRef = $annotation->goIdsByName(name => $chosenIdentifier, aspect => $aspect);
        unless ($goIdListRef) {
            $goToProt{'unmapped'}{$protID}=1;
            next;
        }
        foreach my $goID (@$goIdListRef){
            my $node = $ontology->nodeFromId($goID) or do { push @warnings, "$goID (used to annotate $chosenIdentifier), was not found in selected ontology"; next };
            mapGOtoProt(\%goToProt, $protID, $node);
        }
    }
    print "Mapping to GO: ",time-$time," sec<br>" if $DEBUG;

    # Retrieving GO terms at specified depth #
    $time=time;
    #my %nodeList;
    #getGOatLevel($ontology->rootNode, $depth, \%nodeList);
    my @unsortedNodeList = $ontology->getNodesAtLevel($depth);
    print "Getting specified nodes: ",time-$time," sec<br>" if $DEBUG;
    $time=time;
    my @nodeList = sort { scalar keys %{$goToProt{$b->goid}} <=> scalar keys %{$goToProt{$a->goid}}} grep {$goToProt{$_->goid}} @unsortedNodeList;
    if($goToProt{$UNMAPPED_ID}){
        push @nodeList, new GO::Node( goid => $UNMAPPED_ID , term=> 'Unannotated');
    }
    print "Sorting nodes: ",time-$time," sec<br>" if $DEBUG;

    # Printing category form #
    my %classificationList=&promsMod::getListClass($dbh,$projectID);
    my ($color1, $color2) = promsConfig::getRowColors;

    $dbh->disconnect;
    print header(-'content-encoding'=>'no',-charset=>'UTF-8') unless $DEBUG;
    warningsToBrowser(1);

    if (scalar @warnings) {
		print qq
|<CENTER><INPUT type="button" value="Show warnings" class=\"font11\" onclick="getWarnings(this);">
<DIV id="warnings" style='display:none'>
|;
		foreach (@warnings){
			print "$_<BR><BR>\n";
		}
		print "</DIV></CENTER><BR>";
    }

    
	# Restriction #
	if ($catID) {
		my $numNotInItem=0;
		foreach my $protID (keys %catProteins) {
			$numNotInItem++ unless $prot{$protID};
		}
		print "<FONT style=\"font-weight:bold;color:#DD0000;\">";
		if ($numNotInItem) {
			my $s1Strg=($numNotInItem > 1)? 's' : '';
			print "$numNotInItem protein$s1Strg from List not found in selected Item &nbsp;&bull;&nbsp; ";
		}
		my $numNotInList=scalar keys %notInList;
		my $s2Strg=($numNotInList > 1)? 's' : '';
		print "$numNotInList protein$s2Strg excluded by List restriction.</FONT><BR><BR>\n";
	}

    print qq |
        <table border=0 cellpadding=3 cellspacing=0>
            <tr bgcolor=$color2><td colspan=4 align="center">
                <form onsubmit="toggleProtDomains(document.getElementById('searchProtInput').value); return false;">
                    <label style="font-weight:bold">Protein/gene name filter: <input type="text" id="searchProtInput" value="" /></label>
                    <input id="checkMatchingProtDomains" type="submit" value="Check matching" />
                    <input type="button" value="Uncheck all" onclick="toggleProtDomains('');" />
                    <input type="button" value="Expand/Collapse all" data-collapse="true" onclick="toggleAllTables((this.dataset.collapse == 'true') ? '' : 'none'); this.dataset.collapse = (this.dataset.collapse == 'true') ? false : true;" />
                    <label><input type="checkbox" id="hideUnchecked" onclick="toggleProtDomains(undefined)" />Hide unchecked</label>
                </form>
            </td></tr>
    |;
    
	# Save in List form #
	if ($projectAccess ne 'guest' && $projectStatus <= 0) {
        print qq |
            <tr bgcolor=$color2>
                <TD valign=top><SELECT name="addRemove" id="addRemove" style="font-size:14px;font-weight:bold;text-align:right" onchange="checkClassSelection(this.value)"><OPTION value="add">Add selected proteins to</OPTION><OPTION value="replace">Replace w/ selected proteins</OPTION><OPTION value="remove">Remove selected proteins from</OPTION></SELECT></TD>
                <TD valign=top>
                    <SELECT name="classList" id="classList" style="width:250px;font-weight:bold" onchange="selectClassification()">
        |;
        
		my $displayStatus;
		if (scalar keys %classificationList) {
			print "<OPTION selected value=\"null\">-=Select a Theme=-</OPTION>\n";
			my $classNum=0;
			foreach my $classID (sort{lc($classificationList{$a}[0]) cmp lc($classificationList{$b}[0])} keys %classificationList) {
				$classNum++;
				print "<OPTION value=\"Javascript:changeCatList($classNum)\">$classificationList{$classID}[0]</OPTION>\n";
			}
			$displayStatus='none';
		} else {$displayStatus='block';}
        
		print qq |
                        <OPTION value="-1">-=Create a new Theme=-</OPTION>
                    </SELECT><BR>
                    <TABLE cellpadding=0 cellspacing=0 id="newClass" style="display:$displayStatus"><TR><TD><INPUT type="text" name="newClassName" id="newClassName" value="New Theme name" style="width:250px"/></TD></TR></TABLE>
                </TD>
                <TD valign=top>
                    <SELECT name="catList" id="catList" style="width:250px;font-weight:bold" onchange="checkCreateCat(this.value)" disabled>
                        <OPTION value="null:null">-=Select a List=-</OPTION>
                    </SELECT><BR>
                    <TABLE cellpadding=0 cellspacing=0 id="newCat" style="display:$displayStatus"><TR><TD><INPUT type="text" name="newCatName" id="newCatName" value="New List name" style="width:250px"/></TD></TR></TABLE>
                </TD>
                <TD><INPUT type="button" style="font-size:14px" value="Save" onclick="javascript:addRemFromCategory()">&nbsp</TD>
            </TR>
        |;
	}
    
    print qq |
        </TABLE><BR>
        <DIV id="addProtToListDiv"></DIV>
    |;

    # Printing tables #
    my $dataString;
    my $barColor = $color2;
    while (my $node = shift @nodeList){
        my $goID = $node->goid;
        my $term = formatGOName($node->term);
        my $buttonId = "button$goID";
        my $tableId = "table$goID";
        my $domainId = "domain$goID";
        my $numProt = scalar keys %{$goToProt{$goID}};
        if($numProt < $minCount and $goID ne $UNMAPPED_ID and $goID ne $OTHERS_ID){
            splice @nodeList, -1, 0, new GO::Node( goid => $OTHERS_ID, term=> 'Others') unless $goToProt{$OTHERS_ID};
            @{$goToProt{$OTHERS_ID}}{keys %{$goToProt{$goID}}} = values %{$goToProt{$goID}};
            next;
        }
        my $termString = ($goID eq $UNMAPPED_ID)? $term : "<A href='http://amigo.geneontology.org/amigo/term/$goID' target='_blank'>$term</A>";

        # Header table (containing GO term) #
        my $proteinString = 'protein';
        $proteinString .= 's' if($numProt > 1);
        print qq
|<TABLE border=0 align="center" class="GODomain" id="$domainId" data-name="$term" style="width:1200px;background-color:$color2;margin-top: 10px;">
<TR><TH align=left bgcolor="$color1"><A id="$goID" href="javascript:expandOrCollapseTable('$goID');">
<IMG src="$promsPath{images}/plus.gif" style="position:relative; top:3px" id="$buttonId"/></A><span style="position:relative;bottom:3px">$termString - (<span class="count$goID">$numProt $proteinString</span>)</span></TH></TR></TABLE>
|;

        # Proteins #
        print qq
|<TABLE border=0 cellspacing=0 align="center" class="GOTable" data-checkedprot="0" id="$tableId" style="display:none; margin-top: 5px;"><TR bgcolor=$color2>
<TH class="rbBorder" width=250 align="left"><INPUT type="checkbox" name="checkGO" value="$goID" onclick="selectByGO(this)">&nbsp;<span class="count$goID">$numProt $proteinString</span>&nbsp;</TH>
<TH class="rbBorder" width=120>&nbsp;Gene name&nbsp;</TH>
<TH class="rbBorder" width=80>&nbsp;MW<FONT style="font-size:9px;"> kDa</FONT>&nbsp;</TH>
<TH class="rbBorder" width=80>&nbsp;Peptides&nbsp;</TH>
<TH class="rbBorder" width=80>&nbsp;Spec.<FONT style="font-size:9px;"> %</FONT>&nbsp;</TH>
<TH class="rbBorder" width=80>&nbsp;Cov.<FONT style="font-size:9px;"> %</FONT>&nbsp;</TH>
<TH class="bBorder" width=500>&nbsp;Description - Species&nbsp;</TH>
</TR>
|;
        my $bgColor = $color2;
        foreach my $protID (sort{ $prot{$b}{numPep} <=> $prot{$a}{numPep}} keys %{$goToProt{$node->goid}}){
            $bgColor = ($bgColor eq $color1)? $color2 : $color1;
			my ($geneStrg,$geneName1)=('-',''),; # default
            my $lastIdx=$#{$masterProteins{$prot{$protID}{masterProt}}};
			if ($lastIdx >= 0) {
				$geneName1=$masterProteins{$prot{$protID}{masterProt}}[0]; # cannot use 'shift' because array is used multiple time
				if ($lastIdx >= 1) {
					$geneStrg="<A href=\"javascript:void(null)\" onmouseover=\"popup('<B><U>Synonyms:</U><BR>&nbsp;&nbsp;-".join('<BR>&nbsp;&nbsp;-',@{$masterProteins{$prot{$protID}{masterProt}}}[1..$lastIdx])."</FONT>')\" onmouseout=\"popout()\">$geneName1</A>";
				}
				else {$geneStrg=$geneName1;}
			}
            print qq
|<TR bgcolor="$bgColor" data-protid="$prot{$protID}{identifier}" class="GOTableRow" id="TableRow$goID" data-bgcolor="$bgColor" data-genename="$geneName1" data-genedesc="$prot{$protID}{desc}">
<TH align="left" valign="top">
<INPUT type="checkbox" onchange="checkProt(this, '$goID');" name="checkProt$goID" value="$protID">&nbsp;<A href="javascript:sequenceView($protID,$prot{$protID}{anaID});">$prot{$protID}{identifier}</A>&nbsp;</TD>
<TH valign="top" align="left">&nbsp;$geneStrg&nbsp;</TD>
<TH align="right" valign="top">&nbsp;$prot{$protID}{mw}&nbsp;</TH>
<TH valign="top">&nbsp;$prot{$protID}{numPep}&nbsp;</TH>
<TH valign="top">&nbsp;$prot{$protID}{pepSpe}&nbsp;</TH>
<TH valign="top">&nbsp;$prot{$protID}{pepCov}&nbsp;</TH>
<TD>$prot{$protID}{desc} <FONT class="org">$prot{$protID}{species}</FONT></TD>
</TR>
|;
        }
        print "<TR><TD colspan=6>&nbsp;</TD></TR></TABLE>\n";

        # Data for bar plot
        $barColor = ($barColor eq $color2)? $color1: $color2;
        $dataString .= "$goID||$term||$numProt||$barColor\n";
    }

    print "###DATA###\n";
    print $dataString;
}

sub manageCategory {
    my $dbh = shift;
    my $categoryID = param('CATEGORY');
    my $action = param('ACT');
    my @protID = param('protID');

    my ($className,$catName);
    if($categoryID == -1){ # means there is a new classification and/or a new category to create
        my $classID;
        if(my $newClassName = param('newClassName')){
            my $projectID = param('projectID');
            ($classID) = $dbh->selectrow_array("SELECT MAX(ID_CLASSIFICATION) FROM CLASSIFICATION");
            $classID++;

            my $sthNewClass = $dbh->prepare("INSERT INTO CLASSIFICATION(ID_CLASSIFICATION, ID_PROJECT, NAME, UPDATE_DATE, UPDATE_USER)
                                            VALUES(?,?,?,NOW(),?)");
            $sthNewClass->execute($classID, $projectID, $newClassName, $userID);
            $sthNewClass->finish;

            $className = $newClassName;
        } else {
            $classID = param('CLASS');
            ($className) = $dbh->selectrow_array("SELECT NAME FROM CLASSIFICATION WHERE ID_CLASSIFICATION=$classID");
        }
        my $newCatName = param('newCatName');
        $categoryID = $dbh->selectrow_array("SELECT MAX(ID_CATEGORY) FROM CATEGORY");
        $categoryID++;

        my $sthDisplayPos = $dbh->prepare("SELECT MAX(DISPLAY_POS) FROM CATEGORY WHERE ID_CLASSIFICATION=?");
        $sthDisplayPos->execute($classID);
        my ($displayPos) = $sthDisplayPos->fetchrow_array();
        $sthDisplayPos->finish;
        $displayPos++;

        my $sthNewCat = $dbh->prepare("INSERT INTO CATEGORY(ID_CATEGORY, ID_CLASSIFICATION, NAME, DISPLAY_POS, UPDATE_DATE, UPDATE_USER)
                                      VALUES(?,?,?,?,NOW(),?)");
        $sthNewCat->execute($categoryID, $classID, $newCatName, $displayPos, $userID);
        $sthNewCat->finish;

        $catName = $newCatName;
    } else {
        ($catName) = $dbh->selectrow_array("SELECT NAME FROM CATEGORY WHERE ID_CATEGORY=$categoryID");
    }

    my $actionStrg;
    if($action eq 'add'){
        &addProteinsToCategory($dbh, $categoryID, \@protID);
        $actionStrg = 'were added to';
    }
    elsif($action eq 'replace'){
        &emptyCategory($dbh, $categoryID);
        &addProteinsToCategory($dbh, $categoryID, \@protID);
        $actionStrg = 'have replaced content of';
    }
    elsif($action eq 'remove'){
        &removeProteinsFromCategory($dbh, $categoryID, \@protID);
        $actionStrg = 'were removed from';
    }
    $dbh->commit;
    $dbh->disconnect;

    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    print "<B>Selected proteins $actionStrg Classification <FONT color=#DD0000>$className</FONT>, Category <FONT color=#DD0000>$catName</FONT>.</B><BR><BR>";
}

sub addProteinsToCategory { # assumes proteins not modified-proteins (used by listGO.cgi only!)
	my ($dbh,$categoryID,$protRefList) = @_;
	my $sthChk=$dbh->prepare("SELECT 1 FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=? AND ID_PROTEIN=?");
	my $sthInsCatProt = $dbh->prepare("INSERT INTO CATEGORY_PROTEIN(ID_CATEGORY,ID_PROTEIN) VALUES(?,?)");
	foreach my $protID (@{$protRefList}) {
		$sthChk->execute($categoryID,$protID);
		my ($exists)=$sthChk->fetchrow_array;
		$sthInsCatProt->execute($categoryID,$protID) unless $exists;
	}
	$sthInsCatProt->finish;
	$sthChk->finish;
	$dbh->commit;
}

sub removeProteinsFromCategory { # assumes proteins not modified-proteins (used by listGO.cgi only!)
	my ($dbh,$categoryID,$protRefList) = @_;
	my $sthDelCatProt = $dbh->prepare("DELETE FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=? AND ID_PROTEIN=?");
	foreach my $protID (@{$protRefList}) {
		$sthDelCatProt->execute($categoryID,$protID);
	}
	$sthDelCatProt->finish;
	$dbh->commit;
}

sub emptyCategory { # assumes proteins not modified-proteins (used by listGO.cgi only!)
	my ($dbh,$categoryID) = @_;
	my $sthDelCatProt = $dbh->prepare("DELETE FROM CATEGORY_PROTEIN WHERE ID_CATEGORY=?");
	$sthDelCatProt->execute($categoryID);
	$sthDelCatProt->finish;
	$dbh->commit;
}

sub selectItem {
    my ($name, $id, $hashRef, $sortNumbers, $noSelectionText) = @_;
    my $html = "<FONT class=\"title3\">$name :</FONT><SELECT name=\"$id\" id=\"$id\">\n";

    if ($noSelectionText) {
	$html .= "<OPTION value=0 selected>-= $noSelectionText =-</OPTION>\n";
    }

    foreach my $key (sort {if($sortNumbers){$hashRef->{$a} <=> $hashRef->{$b}} else {$hashRef->{$a} cmp $hashRef->{$b}}} keys %{$hashRef}){
        $html .= "<OPTION value=$key>".$hashRef->{$key}."</OPTION>\n";
    }
    $html .= "</SELECT>\n";
    return $html;
}

sub mapGOtoProt{
    my ($hashRef, $protID, $node) = @_;

    $hashRef->{$node->goid}{$protID} = 1;

    foreach my $parent ($node->parentNodes){
        mapGOtoProt($hashRef, $protID, $parent);
    }
}

#sub getGOatLevel{
#    my ($node, $level, $nodeHashRef) = @_;
#
#    if($level > 0 and !$node->isLeaf){
#        foreach my $child ($node->childNodes){
#            getGOatLevel($child, $level-1, $nodeHashRef);
#        }
#    } else {
#        $nodeHashRef->{$node->goid} = $node;
#    }
#}

sub formatGOName{
	my $goName = shift;

	# applying uppercase to term first letter, unless this letter must be lowercase (eg: mRNA, ncRNA...)
	unless($goName =~ /^\w{1,3}RNA/){
		$goName = ucfirst($goName);
	}

	return $goName;
}

####>Revision history<####
# 1.0.12 [BUGFIX] Fix go summary running from project context (VS 09/09/20)
# 1.0.11 [BUGFIX] Fix few bugs related to gene name displaying (PP 03/08/20)
# 1.0.10 [BUGFIX] Fix gene name displaying when a protein did not have one (VS 22/07/20)
# 1.0.9 [UPDATE] Changed RANK field to IDENT_RANK for compatibility with MySQL 8 (PP 04/03/20) 
# 1.0.8 [FEATURE] Add filters and search field to retrieve specific protein(s) (VS 09/11/19)
# 1.0.7 [FEATURE] Remove locked experiments from GO lists searches (VS 08/08/19)
# 1.0.6 Handles project status=-1 [no auto-end validation] (PP 07/06/18)
# 1.0.5 Compatible with MODIFICATION_SITE (PP 08/08/17)
# 1.0.4 Major bug fix with categories (FY 09/10/13)
# 1.0.3 Restricted to complete OBO files (PP 29/05/13)
# 1.0.2 Added Custom list restriction & gene names (PP 27/05/13)
# 1.0.1 No default selection for GOA and OBO file +<BR> Modifying form style (FY 25/04/13)
# 1.0.0 New script to display present GO terms in current item (FY 13/11/12)
