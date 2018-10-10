#!/usr/local/bin/perl -w

#############################################################################
# manageSpecies.cgi    1.0.9                                                #
# Authors: P. Poullet, G. Arras & F. Yvon (Institut Curie)                  #
# Contact: myproms@curie.fr                                                 #
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
$| = 1;
use strict;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;

####################
####>Parameters<####
####################
my $action = (param('ACT'))? param('ACT'): 'list';
my $nameType=(param('TYPE'))? param('TYPE') : 'C'; # Common name or Scientific name
my $startLetter=(param('LET'))? param('LET') : 'A';


##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

if (param('submit')) {
    &processForm;
    exit;
}

if ($action eq 'list') {
    my %typeList=('C'=>['COMMON_NAME','Common name'],'S'=>['SCIENTIFIC_NAME','Scientific name']);
    my $sthGetSpecies = ($startLetter eq '*')?
          $dbh->prepare("SELECT ID_SPECIES,ID_REF_SPECIES,COMMON_NAME,SCIENTIFIC_NAME,TAXONID,IS_REFERENCE FROM SPECIES WHERE $typeList{$nameType}[0] NOT REGEXP '^[A-Z]' ORDER BY $typeList{$nameType}[0]")
        : $dbh->prepare("SELECT ID_SPECIES,ID_REF_SPECIES,COMMON_NAME,SCIENTIFIC_NAME,TAXONID,IS_REFERENCE FROM SPECIES WHERE $typeList{$nameType}[0] LIKE '$startLetter%' ORDER BY $typeList{$nameType}[0]");
    $sthGetSpecies->execute;
    my $speciesRef = $sthGetSpecies->fetchall_arrayref();
    $sthGetSpecies->finish;
    my (%taxonRefSpecies,%speciesUsed);
    my $sthRefSpecies=$dbh->prepare("SELECT TAXONID FROM SPECIES WHERE ID_SPECIES=?");
    foreach my $entry (@{$speciesRef}) {
		$speciesUsed{$$entry[0]}=1 if $$entry[5];
        next unless $$entry[1]; # ID_REF_SPECIES
        $sthRefSpecies->execute($$entry[1]);
        ($taxonRefSpecies{$$entry[1]})=$sthRefSpecies->fetchrow_array;
    }
    $sthRefSpecies->finish;
    my @sthSpUsed=($dbh->prepare("SELECT DISTINCT S.ID_SPECIES FROM SPECIES S,IDENTIFIER I WHERE S.ID_SPECIES=I.ID_SPECIES"), # used by an identifier type
                   $dbh->prepare("SELECT DISTINCT S.ID_SPECIES FROM SPECIES S,MASTER_PROTEIN M WHERE S.ID_SPECIES=M.ID_SPECIES"), # used by a master protein
                   $dbh->prepare("SELECT DISTINCT S.ID_SPECIES FROM SPECIES S,GOANNOTATION G WHERE S.ID_SPECIES=G.ID_SPECIES"), # used in a GO analysis
                   $dbh->prepare("SELECT DISTINCT R.ID_SPECIES FROM SPECIES R,SPECIES S WHERE R.ID_SPECIES=S.ID_REF_SPECIES") # used as reference for another species
                  );
    foreach my $sth (@sthSpUsed) {
        $sth->execute;
        while (my ($spID)=$sth->fetchrow_array) {
            $speciesUsed{$spID}=1;
        }
        $sth->finish;
    }
    $dbh->disconnect;
    my ($light,$dark)=&promsConfig::getRowColors;
    my $bgColor=$light;

    print header;
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<TITLE>Species List</TITLE>
<SCRIPT language="Javascript">
function deleteSpecies(idSpecies,speciesName){
    if(confirm('Delete '+speciesName+' ?')){
        window.location='$promsPath{cgi}/manageSpecies.cgi?ACT=delete&TYPE=$nameType&LET=$startLetter&ID='+idSpecies;
    }
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<DIV style="float:top">
<BR>
<TABLE><TR><TH bgcolor="$dark">
<FONT class="title2">&nbsp;Go to:</FONT><SELECT class="title action" onchange="if (this.value) {window.location='$promsPath{cgi}/'+this.value;}">
	<OPTION value="">-= Select =-</OPTION>
	<OPTION value="promsMain.cgi">Main Window</OPTION>
	<OPTION value="selectProject.cgi">Project Selection</OPTION>
</SELECT>
</TH></TR></TABLE>
</DIV>
<CENTER>
<FONT class="title">List of Species</FONT>&nbsp;&nbsp;<INPUT type="button" class="title2" value="Add species" onclick="window.location='$promsPath{cgi}/manageSpecies.cgi?ACT=add&TYPE=$nameType&LET=$startLetter';">
<BR><BR>
<TABLE bgcolor=$dark>
|;
    foreach my $type (sort keys %typeList) {
        print "<TR><TH align=right class=\"title3\">&nbsp;$typeList{$type}[1]:</TH><TH bgcolor=$light>";
        foreach my $letter ('A'..'Z','*') {
            my $letterClass=($type eq $nameType && $letter eq $startLetter)? 'selected' : 'selectable';
            print "<A class=\"$letterClass\" href=\"$promsPath{cgi}/manageSpecies.cgi?ACT=list&TYPE=$type&LET=$letter\">&nbsp;$letter&nbsp;</A>";
            print "|" if $letter ne '*';
        }
        print "</TH></TR>\n";
    }
    print qq
|</TABLE>
<BR>
|;
    if (scalar(@{$speciesRef})) {
        print qq
|<TABLE border=0 cellspacing=0>
<TR bgcolor=$dark>
    <TH class="rbBorder title3" width=300>Common name</TH>
    <TH class="rbBorder title3" width=400>Scientific name</TH>
    <TH class="rbBorder title3" width=200>Taxon ID</TH>
    <TH class="bBorder title3">Action</TH>
</TR>
|;
        foreach my $entry (@{$speciesRef}) {
            my ($speciesID,$refSpeciesID,$comName,$sciName,$taxonId,$isReference)=@{$entry};
            $taxonId = '?' unless $taxonId;
            $taxonId .= " &#8658; <A href=\"$promsPath{cgi}/manageSpecies.cgi?ACT=edit&ID=$refSpeciesID&TYPE=$nameType&LET=$startLetter\">$taxonRefSpecies{$refSpeciesID}</A>" if $refSpeciesID;
            my $disabDel=($speciesUsed{$speciesID})? 'disabled' : '';
			my $isRefIconStrg=($isReference)? "<IMG src=\"$promsPath{images}/good.gif\">&nbsp;" : '';
            print qq
|<TR align="center" bgcolor=$bgColor>
    <TD class="title3" bgcolor=$bgColor>$isRefIconStrg$comName</TD>
    <TD><I>$sciName</I></TD>
    <TD>$taxonId</TD>
    <TD>
        <INPUT type="button" value="Edit" onclick="window.location='$promsPath{cgi}/manageSpecies.cgi?ACT=edit&ID=$speciesID&TYPE=$nameType&LET=$startLetter';">
        <INPUT type="button" value="Delete" onclick="deleteSpecies($speciesID,'$comName');" $disabDel>
    </TD>
</TR>
|;
            $bgColor = ($bgColor eq $light)?$dark:$light;
        }
        print "</TABLE>\n";
    }
    else {print "<BR><FONT class=\"title2\">No matching species found.</FONT><BR>\n";}
    print qq
|<BR><BR>
</CENTER>
</BODY>
</HTML>
|;
}
elsif ($action eq 'add' || $action eq 'edit'){
    my ($comName,$sciName,$taxonID,$isReference) = ('','','',0);
    my $idSpecies=(param('ID'))?param('ID'):undef;
    my ($refSpeciesID,@refSpeciesInfo,$hasSubSpecies);
    my $title;
    if ($action eq 'edit') {
        ($refSpeciesID,$comName,$sciName,$taxonID,$isReference) = $dbh->selectrow_array("SELECT ID_REF_SPECIES,COMMON_NAME,SCIENTIFIC_NAME,TAXONID,IS_REFERENCE FROM SPECIES WHERE ID_SPECIES=$idSpecies");
        if ($refSpeciesID) {
            @refSpeciesInfo=$dbh->selectrow_array("SELECT COMMON_NAME,SCIENTIFIC_NAME,TAXONID FROM SPECIES WHERE ID_SPECIES=$refSpeciesID");
        }
        ($hasSubSpecies)=$dbh->selectrow_array("SELECT COUNT(*) FROM SPECIES WHERE ID_REF_SPECIES=$idSpecies");
        $title = "Edit Species $comName";
    }
    else {$title = 'Add Species';}
    my @refSpecies;
    my $sthRefSpecies = $dbh->prepare("SELECT ID_SPECIES, COMMON_NAME, SCIENTIFIC_NAME FROM SPECIES WHERE IS_REFERENCE=1 ORDER BY COMMON_NAME");
    $sthRefSpecies->execute();
    while (my ($idSpecies,$cName,$sName) = $sthRefSpecies->fetchrow_array) {
	push @refSpecies, { ID => $idSpecies, commonName => $cName, scientificName => $sName};
    }
    $sthRefSpecies->finish;

    #my ($disabledGOAStrg,$checkedGOAStrg) = ($goaUrl)? ('','checked'): ('disabled','');
    my ($light,$dark)=&promsConfig::getRowColors;
    print header;
    warningsToBrowser(1);
    print qq
|<HEAD>
<SCRIPT language="Javascript">
function checkForm(myForm){
    if(!myForm.taxonID.value){
        alert('Enter a taxon ID.');
        return false;
    } else if (!myForm.comName.value){
        alert('Enter a common name.');
        return false;
    } else if (!myForm.sciName.value){
        alert('Enter a scientific name.');
        return false;
    } else {
        return true;
    }
}
function searchTaxonomy(myForm) {
    var searchText=(myForm.taxonID.value)? '?term='+myForm.taxonID.value : (myForm.comName.value)?'?term='+myForm.comName.value : (myForm.sciName.value)? '?term='+myForm.sciName.value : '';
    window.open('http://www.ncbi.nlm.nih.gov/taxonomy'+searchText,'NCBI_Tax');
}
function checkIsRef(box){
    var selectRefSp = document.getElementById('selectRefSpecies');
    if(box.checked){
	selectRefSp.value = 0;
	selectRefSp.disabled = true;
    } else {
	selectRefSp.disabled = false;
    }
}
</SCRIPT>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<TITLE>$title</TITLE>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR>
<FONT class="title">$title</FONT><BR>
<BR>
<FORM name="addSpeciesForm" method=POST onsubmit="return(checkForm(this));">
|;
    if ($action eq 'edit') {
        print qq
|<INPUT type="hidden" name="idSpecies" value="$idSpecies">
<INPUT type="hidden" name="TYPE" value="$nameType">
|;
    }
    my $checkRefString = ($isReference)? 'checked': '';
    my $checkRefComment;
    if ($hasSubSpecies) {
	$checkRefString .= ' disabled';
	$checkRefComment = "Used as reference for $hasSubSpecies other species";
    } else {
	$checkRefComment = "<I>Only reference species can be used for data processing (eg. GO analyses)</I>"
    }

    print qq
|<TABLE border=0 cellpadding=2 bgcolor=$dark>
<TR><TH align=right>Common name :</TH><TD bgcolor=$light><INPUT type="text" name="comName" value="$comName" size=40></TD></TR>
<TR><TH align=right>Scientific name :</TH><TD bgcolor=$light><INPUT type="text" name="sciName" value="$sciName" size=70></TD></TR>
<TR><TH align=right>Taxon ID :</TH><TD bgcolor=$light><INPUT type="text" name="taxonID" value="$taxonID" size=10>&nbsp;[<A href="javascript:searchTaxonomy(document.addSpeciesForm)">Search on NCBI Taxonomy</A>]</TD></TR>
<TR><TH align=right>Is reference :</TH><TD bgcolor=$light><INPUT type="checkbox" name="isReference" value=1 onclick="checkIsRef(this);" $checkRefString>&nbsp;$checkRefComment</TD></TR>
<TR><TH align=right valign=top>Reference species :</TH><TD bgcolor=$light nowrap>
|;
    my $refDisabledString = ($isReference)? ' disabled': '';
    print "<SELECT id=\"selectRefSpecies\" name=\"refSpecies\"$refDisabledString>\n<OPTION value=0 selected>No reference</OPTION>\n";
    foreach my $refSpecies (@refSpecies){
	my %species = %{$refSpecies};
	my $selectedString = ($refSpeciesID && $refSpeciesID==$species{ID})? 'selected':'';
	print "<OPTION value=$species{ID} $selectedString>$species{commonName} ($species{scientificName})</OPTION>\n";
    }
    print qq
|</SELECT>
</TD>
</TR>
<TR><TD colspan=2 align=center>
|;
    my $submitValue=($action eq 'add')? 'Add' : 'Save';
    print qq
|<INPUT type="submit" name="submit" value="$submitValue">
&nbsp;&nbsp;&nbsp;<INPUT type="reset" value="  Clear  ">
&nbsp;&nbsp;&nbsp;<INPUT type="button" value="Cancel" onclick="window.location='./manageSpecies.cgi?ACT=list&TYPE=$nameType&LET=$startLetter';">
</TD></TR>
</TABLE>
|;
}
elsif ($action eq 'delete') {
    my $speciesID = param('ID');
    $dbh->do("DELETE FROM SPECIES WHERE ID_SPECIES=$speciesID");
    $dbh->commit;
    $dbh->disconnect;

    print header;
    warningsToBrowser(1);
    print qq
|<HEAD>
</HEAD>
<BODY onload="window.location='./manageSpecies.cgi?ACT=list&TYPE=$nameType&LET=$startLetter';">
</BODY>
</HTML>
|;
}

# add or edit submit action
sub processForm{
    my $comName = param('comName');
    my $sciName = param('sciName');
    my $taxonID = param('taxonID');
    my $isReference = (param('isReference'))? param('isReference') : 0;
    my $idSpecies = param('idSpecies'); # undef if action='add'
    #my $goaURL = (param('goaURL'))? param('goaURL') : undef;
    #my $taxonRefSpecies=param('refTaxonID');
    my $refSpeciesID = (param('refSpecies'))? param('refSpecies') : 0;
#    if ($taxonRefSpecies) {
#        ($refSpeciesID)=$dbh->selectrow_array("SELECT ID_SPECIES FROM SPECIES WHERE TAXONID='$taxonRefSpecies'");
#        unless($refSpeciesID){
#            print header;
#            warningsToBrowser(1);
#            print qq
#|<HTML>
#<BODY background="$promsPath{images}/bgProMS.gif">
#<FONT color="red">Error: Reference species with Taxon ID '$taxonRefSpecies' not recorded!</FONT><BR>
#<INPUT type="button" value="Go back" onclick="javascript:history.back();">
#</BODY>
#</HTML>
#|;
#            exit;
#        }
#    }
    my $sthUpDB;
    if ($idSpecies) { #update
        $sthUpDB = $dbh->prepare("UPDATE SPECIES SET ID_REF_SPECIES=?,COMMON_NAME=?,SCIENTIFIC_NAME=?,TAXONID=?,IS_REFERENCE=? WHERE ID_SPECIES=$idSpecies");
    }
    else { #add
        # Check if taxon ID already exists
        my $sthTaxonExists = $dbh->prepare("SELECT ID_SPECIES FROM SPECIES WHERE TAXONID=?");
        $sthTaxonExists->execute($taxonID);
        if ($sthTaxonExists->fetchrow_array) {
            $sthTaxonExists->finish;
            $dbh->disconnect;

            print header;
            warningsToBrowser(1);
            print qq
|<HTML>
<BODY background="$promsPath{images}/bgProMS.gif">
<FONT color="red">Error: Taxon ID '$taxonID' already assigned!</FONT><BR>
<INPUT type="button" value="Go back" onclick="javascript:history.back();">
</BODY>
</HTML>
|;
            exit;
        }
        else {
            $sthUpDB = $dbh->prepare("INSERT INTO SPECIES (ID_REF_SPECIES,COMMON_NAME,SCIENTIFIC_NAME,TAXONID,IS_REFERENCE) VALUES (?,?,?,?,?)"); # PK autoincrement
        }
        $sthTaxonExists->finish;
    }
    my @params = ($refSpeciesID,$comName,$sciName,$taxonID,$isReference);
    $params[0] = undef unless $params[0]; # setting NULL in SQL
    $sthUpDB->execute(@params) || die $dbh->errstr;

    $sthUpDB->finish;
    $dbh->commit;
    $dbh->disconnect;

    print header;
    warningsToBrowser(1);
    ($startLetter)=($nameType eq 'C')? ($comName=~/^(.)/) : ($sciName=~/^(.)/);
    print qq
|<HTML>
<BODY onload="window.location='./manageSpecies.cgi?ACT=list&TYPE=$nameType&LET=$startLetter';">
</BODY>
</HTML>
|;
}

####>Revision history<####
# 1.0.9 Color on main "Go to" options (PP 25/03/16)
# 1.0.8 Change syntax $x=$y // $z to old ()? one (PP 10/03/14)
# 1.0.7 No deletion of reference species & added icon (PP 24/09/13)
# 1.0.6 IS_REFERENCE management (FY 16/09/13)
# 1.0.5 Removing GOA url fetching & fix DB column misnaming (FY 22/03/13)
# 1.0.4 Added reference species management (PP 31/01/13)
# 1.0.3 Added A to Z menu for species group selection: Requires promsStyle.css 2.0.7 (PP 04/01/13)
# 1.0.2 Called by promsMain.cgi && prevents deletion if species is used (PP 13/12/12)
# 1.0.1 Manage displaying if no species (FY 26/11/12)
# 1.0.0 New script to manage species (FY 13/11/12)