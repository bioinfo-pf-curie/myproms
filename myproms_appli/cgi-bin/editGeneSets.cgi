#!/usr/local/bin/perl -w

################################################################################
# editGeneSets.cgi         1.0.1                                               #
# Authors: V. Laigle (Institut Curie)                                          #
# Contact: myproms@curie.fr                                                    #
# Create, edit or delete Gene Sets databanks (used for GSEA)                   #
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

$| = 1;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use POSIX qw(strftime); # to get the time
use strict;
use warnings;
use File::Copy;
use File::Basename;
use promsConfig;
use promsMod;

#####################
### Configuration ###
#####################
my %promsPath = &promsConfig::getServerInfo;
my ($lightColor, $darkColor) = &promsConfig::getRowColors;
my $userID = $ENV{'REMOTE_USER'};
my $dbh = &promsConfig::dbConnect;
my $action = param('ACT');
my $geneSetsID = (param('ID'))? param('ID') : 0;  # 0 if ACT=add

if ($action eq 'add') {
    if (param('save')) {  # Form was sent
        &createGeneSetsDB;
    } else {
        &printAddGeneSetsForm;
    }
} elsif ($action eq 'edit') {
    if (param('save')) {  # Form was sent
        &editGeneSetsDB($geneSetsID);
    } else {
        &printEditGeneSetsForm($geneSetsID);
    }
} elsif ($action eq 'update') {
    if (param('save')) {  # Form was sent
        &updateGeneSetsDB($geneSetsID);
    } else {
        &printUpdateGeneSetsForm($geneSetsID);
    }
} elsif ($action eq 'restore') {
    if (param('save')) {  # Form was sent
        &restoreGeneSetsDb($geneSetsID);
    } else {
        &printRestoreGeneSetsForm($geneSetsID);
    }
} elsif ($action eq 'delete') {
    &deleteGeneSetsDB($geneSetsID);
}
$dbh->disconnect;
exit;


########################
### Adding Gene Sets ###
########################
sub printAddGeneSetsForm {  # Globals: %promsPath, $dbh, $userID, $action, $lightColor, $darkColor

    my $resetString = ($action eq 'add')? ' Clear ' : ' Clear Changes ';

    # Get all known species in DB (separate reference or not and sort alphabetically)
    my $speciesStrg = "<OPTION value=\"\">-= Select =-</OPTION>\n";
    $speciesStrg .= "<OPTION value=\"Custom\">No specific species (Custom or muti-species)</OPTION>\n";
    my $sthGetSpecies = $dbh->prepare("SELECT ID_SPECIES, COMMON_NAME, SCIENTIFIC_NAME FROM SPECIES WHERE IS_REFERENCE IS NOT NULL ORDER BY SCIENTIFIC_NAME, COMMON_NAME");
    $sthGetSpecies->execute;
    while (my ($speciesID, $spCName, $spSName) = $sthGetSpecies->fetchrow_array) {
        $speciesStrg .= "<OPTION value=\"$speciesID\">$spSName ($spCName)</OPTION>\n";
    }
    $sthGetSpecies->finish;

    # Get all identifier types in DB (sort alphabetically)
    my $identStrg = "<OPTION value=\"\">-= Select =-</OPTION>\n";
    $identStrg .= "<OPTION value=\"Custom\">Custom identifiers</OPTION>\n";
    my $sthGetIdent = $dbh->prepare("SELECT ID_IDENTIFIER, NAME FROM IDENTIFIER ORDER BY NAME");
    $sthGetIdent->execute;
    while (my ($identID, $identName) = $sthGetIdent->fetchrow_array) {
        $identStrg .= "<OPTION value=\"$identID\">$identName</OPTION>\n";
    }
    $sthGetIdent->finish;

    # Get all projects visible by user (group by owner as in project selection and sort alphabetically)
    my (%projects, $projectStrg);
    my ($userStatus, $adminProjects) = $dbh->selectrow_array(
        "SELECT UL.USER_STATUS, GROUP_CONCAT(PA.ID_PROJECT SEPARATOR ',')
        FROM USER_LIST UL
        LEFT JOIN PROJECT_ACCESS PA ON UL.ID_USER = PA.ID_USER
        LEFT JOIN USER_PROFILE UP ON PA.ID_PROFILE = UP.ID_PROFILE
        WHERE UL.ID_USER = '$userID' AND (UP.NAME IS NULL OR UP.NAME LIKE '\%administrator\%')
        GROUP BY UL.ID_USER"
    );

    my $projectQuery;
    if ($userStatus =~ /mass|bioinfo/) {
        $projectQuery = "SELECT ID_PROJECT, NAME, OWNER, WORK_GROUP FROM PROJECT";  # Select all projects for massists and bioinfo
    } else {
        $projectQuery = "SELECT ID_PROJECT, NAME, OWNER, WORK_GROUP FROM PROJECT WHERE ID_PROJECT IN ($adminProjects)";  # Select only projects where user is at least admin
    }

    my $sthGetProject = $dbh->prepare($projectQuery);
    $sthGetProject->execute;
    while (my ($projectID, $projectName, $projectOwner, $projectWG) = $sthGetProject->fetchrow_array) {
        $projectWG = "None" unless ($projectWG);
        $projects{$projectWG}{$projectOwner}{$projectID} = $projectName;
    }
    $sthGetProject->finish;

    my $indent = "&nbsp;&nbsp;&nbsp;&nbsp;";
    $projectStrg = "<OPTION value=\"\">-= Select =-</OPTION>\n";
    $projectStrg .= "<OPTION value=\"Public\">Public databank (No specific project)</OPTION>\n";
    foreach my $projOwner (sort{$a cmp $b} keys %{$projects{'None'}}) {  # All projects without work group first (usual)
        $projectStrg .= "<OPTGROUP label=\"$projOwner\"></OPTGROUP>\n";
        # $projectStrg .= "<OPTION value=\"owner_no_wg_$projOwner\">${indent}All projects from '$projOwner'</OPTION>\n";
        foreach my $projID (sort{$projects{'None'}{$projOwner}{$a} cmp $projects{'None'}{$projOwner}{$b}} keys %{$projects{'None'}{$projOwner}}) {
            $projectStrg .= "<OPTION value=\"$projID\">${indent}$projects{'None'}{$projOwner}{$projID}</OPTION>\n";
        }
    }

    foreach my $projWG (sort{$a cmp $b} keys %projects) {
        next if ($projWG eq "None");
        $projectStrg .= "<OPTGROUP label=\"WG: $projWG\"></OPTGROUP>\n";
        # $projectStrg .= "<OPTION value=\"wg_$projWG\">${indent}All projects from work group '$projWG'</OPTION>\n";
        foreach my $projOwner (sort{$a cmp $b} keys %{$projects{$projWG}}) {
            $projectStrg .= "<OPTGROUP label=\"${indent}$projOwner\"></OPTGROUP>\n";
            # $projectStrg .= "<OPTION value=\"owner_$projOwner\_wg_$projWG\">${indent}${indent}All projects from '$projOwner' in WG '$projWG'</OPTION>\n";
            foreach my $projID (sort{$projects{$projWG}{$projOwner}{$a} cmp $projects{$projWG}{$projOwner}{$b}} keys %{$projects{$projWG}{$projOwner}}) {
                $projectStrg .= "<OPTION value=\"$projID\">${indent}${indent}$projects{$projWG}{$projOwner}{$projID}</OPTION>\n";
            }
        }
    }


    ####################
    ### Display HTML ###
    ####################
    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<TITLE>Add New Gene Sets Databank</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT type="text/javascript">
// Dictionnary of help popup texts (avoid using '"')
const helpText= {
    gmtFile: 'A .gmt file containing your custom gene sets. Click to get to the Broad Institute wiki on GSEA and GMT files for explanations on how a .gmt file is formatted.',
    refSpecies: 'Only the reference species are available here. If the desired species is not present, ask us to set your species as reference.'
};
|;
    &promsMod::popupInfo();
    print qq
|function cancelAction() {
    window.location = './manageGeneSets.cgi?ACT=display';
}
function checkForm(myForm) {
    if (!myForm.gsName.value) {
        alert('Give a name to this Gene Sets Databank.');
        return false;
    } else if (!myForm.species.value) {
        alert('Select the species corresponding to this Gene Sets Databank. If no particular species is used, select this option.');
        return false;
    } else if (!myForm.identifier.value) {
        alert('Select a type of identifiers for this Gene Sets Databank. It can be a custom one.');
        return false;
    } else if (!myForm.gsProject.value) {
        alert('Select a project for this Gene Sets Databank or set it as public.');
        return false;
    } else if (!myForm.gmtFile.value) {
        alert('Provide a .gmt file for this Gene Sets Databank.');
        return false;
    } else if (!myForm.gsVersionDate.value) {
        alert('Provide the version date for this Gene Sets Databank.');
        return false;
    }
    return true;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Add new Gene Sets Databank</FONT>
<BR><BR>
<FORM name="newGSForm" method="post" onsubmit="return checkForm(this);" enctype="multipart/form-data">
<INPUT type="hidden" name="ACT" value="$action">
<TABLE border=0>
    <TR>
        <TD bgcolor=$darkColor>
            <TABLE border=0 cellpadding=2>
                <TR>
                    <TH align="right" valign="top" nowrap>Name of the Gene Sets Databank :</TH>
                    <TD bgcolor="$lightColor" nowrap>
                        <INPUT type="text" id="gsName" name="gsName" maxlength="50">
                    </TD>
                </TR>
                <TR>
                    <TH align="right" valign="top" nowrap>Description :</TH>
                    <TD bgcolor="$lightColor" nowrap>
                        <TEXTAREA id="gsDes" name="gsDes" rows="2" cols="60" maxlength="100"></TEXTAREA>
                    </TD>
                </TR>
                <TR>
                    <TH align="right" valign="top" nowrap>Species<SUP onmouseover="popup(helpText.refSpecies);" onmouseout="popout();">?</SUP>  :</TH>
                    <TD bgcolor="$lightColor" nowrap>
                        <SELECT id="species" name="species">
                            $speciesStrg
                        </SELECT>
                    </TD>
                </TR>
                <TR>
                    <TH align="right" valign="top" nowrap>Identifier type :</TH>
                    <TD bgcolor="$lightColor" nowrap>
                        <SELECT id="identifier" name="identifier">
                            $identStrg
                        </SELECT>
                    </TD>
                </TR>
                <TR>
                    <TH align="right" valign="top" nowrap>Link databank to a specific project :</TH>
                    <TD bgcolor="$lightColor" nowrap>
                        <SELECT id="gsProject" name="gsProject">
                            $projectStrg
                        </SELECT>
                    </TD>
                </TR>
                <TR>
                    <TH align="right" valign="top" nowrap>Gene Sets file (.gmt)<SUP onmouseover="popup(helpText.gmtFile);" onmouseout="popout();"><A href="https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29" target="_blank">?</A></SUP> :</TH>
                    <TD bgcolor="$lightColor" nowrap>
                        <INPUT type="file" name="gmtFile" id="gmtFile" accept=".gmt,.GMT">
                    </TD>
                </TR>
                <TR>
                    <TH align="right" valign="top" nowrap>Gene Sets Version Date</TH>
                    <TD bgcolor="$lightColor" nowrap>
                        <INPUT type="date" id="gsVersionDate" name="gsVersionDate">
                    </TD>
                </TR>
                <TR>
                    <TD colspan=2 align=center>
                        <INPUT type="submit" name="save" value=" Save ">
                        &nbsp;&nbsp;&nbsp;
                        <INPUT type="reset" value="$resetString">
                        &nbsp;&nbsp;&nbsp;
                        <INPUT type="button" value="Cancel" onclick="cancelAction()">
                    </TD>
                </TR>
            </TABLE>
        </TD>
    </TR>
</TABLE>
</FORM>
</CENTER>
<BR><BR>
<DIV id="divDescription" class="clDescriptionCont"></DIV>
<SCRIPT type="text/javascript">
    setPopup();
</SCRIPT>
</BODY>
</HTML>
|;

}


sub createGeneSetsDB {  # Globals: %promsPath, $dbh, $userID

    ####################
    ### Display HTML ###
    ####################
    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<TITLE>Manage Gene Sets Databank</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<FONT class="title">Adding new Gene Sets Databank...</FONT>
|;

    my $geneSetsID;
    my ($gsName, $gsDes, $speciesID, $identID, $projectID, $userGmtFile, $gsFileType, $numSets, $gsVersionDate, $gsStatus);
    $gsName         = param('gsName');
    $gsDes          = param('gsDes') || '-';
    $speciesID      = (param('species') eq "Custom") ? "NULL" : param('species');
    $identID        = (param('identifier') eq "Custom")? "NULL" : param('identifier');
    $projectID      = (param('gsProject') eq "Public")? "NULL" : param('gsProject');
    $userGmtFile    = (param('gmtFile'))? &promsMod::cleanParameters(param('gmtFile')) : '';
    $gsVersionDate  = param('gsVersionDate') || "CURDATE()";
    $gsStatus       = 1;  # Set Gene Sets in use automatically

    unless ($userGmtFile) {  # Double check, should not pass checkForm()
        $dbh->disconnect;
        die "There was no file provided for the gene sets !";
    }
    ($gsFileType) = $userGmtFile =~ /\.([^.]+)$/;
    $gsFileType = uc($gsFileType);
    my $tmpGmtFile = tmpFileName($userGmtFile);
    ($numSets) = (`wc -l $tmpGmtFile` =~ /^(\d+)\s/);

    my $sthInsNewGS = $dbh->prepare("INSERT INTO GENESETS (ID_SPECIES, ID_IDENTIFIER, NAME, DES, ID_PROJECT, GENESETS_FILE, GENESETS_TYPE, NUM_SETS, GENESETS_STATUS, VERSION_DATE, UPDATE_DATE, UPDATE_USER) VALUES ($speciesID, $identID, \"$gsName\", \"$gsDes\", $projectID, \"$tmpGmtFile\", \"$gsFileType\", $numSets, $gsStatus, \"$gsVersionDate\", NOW(), \"$userID\");");
    $sthInsNewGS->execute;
    $geneSetsID = $dbh->last_insert_id(undef, undef, 'GENESETS', 'ID_GENESETS');
    $sthInsNewGS->finish;

    # Create new gene sets DB folder with new ID, move gmt file in it and update DB
    my ($gmtFileDir, $gmtFilePath);
    if ($projectID eq "NULL") {
        $gmtFileDir = "$promsPath{'genesets'}/public/gsDB_$geneSetsID";
        $gmtFilePath = "public/gsDB_$geneSetsID/$userGmtFile";
        unless (-d "$gmtFileDir") {
            unless (-d "$promsPath{'genesets'}/public") {
                mkdir "$promsPath{'genesets'}" unless (-d "$promsPath{'genesets'}");
                mkdir "$promsPath{'genesets'}/public";
            }
            mkdir "$gmtFileDir";
        }
    } else {
        $gmtFileDir = "$promsPath{'genesets'}/private/project_$projectID/gsDB_$geneSetsID";
        $gmtFilePath = "private/project_$projectID/gsDB_$geneSetsID/$userGmtFile";
        unless (-d "$gmtFileDir") {
            unless (-d "$promsPath{'genesets'}/private/project_$projectID") {
                unless (-d "$promsPath{'genesets'}/private") {
                    mkdir "$promsPath{'genesets'}" unless (-d "$promsPath{'genesets'}");
                    mkdir "$promsPath{'genesets'}/private";
                }
                mkdir "$promsPath{'genesets'}/private/project_$projectID";
            }
            mkdir "$gmtFileDir";
        }
    }
    move($tmpGmtFile, "$promsPath{'genesets'}/$gmtFilePath");
    $dbh->do("UPDATE GENESETS SET GENESETS_FILE = \"$gmtFilePath\" WHERE ID_GENESETS = $geneSetsID") || die $dbh->errstr;
    $dbh->commit;

    print qq
|<BR><BR>
<FONT class="title3">
    Done.
    <BR>
    Creating archive for this version of the Gene Sets databank.
</FONT>
|;
    # Create archive for this Gene Sets DB version
    &createArchive($geneSetsID);

    print qq
|<BR><BR>
<FONT class="title3">
    Done.
    <BR>
    This page will refresh itself in a few seconds.
</FONT>
<BR>
|;
    sleep 5;
    print qq
|<SCRIPT type="text/javascript">
    window.location = "$promsPath{cgi}/manageGeneSets.cgi?ACT=display";
</SCRIPT>
</BODY>
</HTML>
|;

}


#########################
### Editing Gene Sets ###
#########################
sub printEditGeneSetsForm {  # Globals: %promsPath, $dbh, $action, $lightColor, $darkColor
    my ($geneSetsID) = @_;
    my ($gsName, $gsDes, $gsSpeciesID, $gsIdentID, $gsStatus) = $dbh->selectrow_array("SELECT NAME, DES, ID_SPECIES, ID_IDENTIFIER, GENESETS_STATUS FROM GENESETS WHERE ID_GENESETS = $geneSetsID");
    $gsSpeciesID = "Custom" unless ($gsSpeciesID);
    $gsIdentID   = "Custom" unless ($gsIdentID);

    my $selSpecies = ($gsSpeciesID eq "Custom")? " selected" : " ";
    my $selIdent   = ($gsIdentID eq "Custom")? " selected" : " ";

    my $resetString = ($action eq 'add')? ' Clear ' : ' Clear Changes ';

    # Get all known species in DB (separate reference or not and sort alphabetically)
    my $isRef = 1;
    my $speciesStrg = "<OPTION value=\"\">-= Select =-</OPTION>\n";
    $speciesStrg .= "<OPTION value=\"Custom\"$selSpecies>No specific species (Custom or muti-species)</OPTION>\n";
    $speciesStrg .= "<OPTGROUP label=\"Reference species\">\n";
    my $sthGetSpecies = $dbh->prepare("SELECT ID_SPECIES, COMMON_NAME, SCIENTIFIC_NAME, IS_REFERENCE FROM SPECIES ORDER BY IS_REFERENCE DESC, SCIENTIFIC_NAME, COMMON_NAME");
    $sthGetSpecies->execute;
    while (my ($speciesID, $spCName, $spSName, $spIsRef) = $sthGetSpecies->fetchrow_array) {
        $selSpecies = ($gsSpeciesID eq $speciesID)? " selected" : " ";  # Cannot set ==, gsSpeciesID may be a string
        if ($isRef && !$spIsRef) {
            $isRef = 0;
            $speciesStrg .= "</OPTGROUP>\n";
            $speciesStrg .= "<OPTGROUP label=\"Other species\">\n";
        }
        $speciesStrg .= "<OPTION value=\"$speciesID\"$selSpecies>$spSName ($spCName)</OPTION>\n";
    }
    $speciesStrg .= "</OPTGROUP>\n";
    $sthGetSpecies->finish;

    # Get all identifier types in DB (sort alphabetically)
    my $identStrg = "<OPTION value=\"\">-= Select =-</OPTION>\n";
    $identStrg .= "<OPTION value=\"Custom\"$selIdent>Custom identifiers</OPTION>\n";
    my $sthGetIdent = $dbh->prepare("SELECT ID_IDENTIFIER, NAME FROM IDENTIFIER ORDER BY NAME");
    $sthGetIdent->execute;
    while (my ($identID, $identName) = $sthGetIdent->fetchrow_array) {
        $selIdent = ($gsIdentID eq $identID)? " selected" : " ";  # Cannot set ==, gsIdentID may be a string
        $identStrg .= "<OPTION value=\"$identID\"$selIdent>$identName</OPTION>\n";
    }
    $sthGetIdent->finish;

    # Build status select string
    # my @statusArr = ([1, 'In Use'], [0, 'Not In Use'], [-1, 'Error']);
    # my $selStatus = "";
    my $statusStrg = "";
    # foreach my $stat (@statusArr) {
    #     $selStatus = ($stat->[0] == $gsStatus) ? " selected" : " ";
    #     $statusStrg .= "<OPTION value=\"$stat->[0]\"$selStatus>$stat->[1]</OPTION>\n";
    # }

    ####################
    ### Display HTML ###
    ####################
    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<TITLE>Edit Gene Sets Databank</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT type="text/javascript">
|;
    &promsMod::popupInfo();
    print qq
|function cancelAction() {
    window.location = './manageGeneSets.cgi?ACT=display';
}
function checkForm(myForm) {
    if (!myForm.gsName.value) {
        alert('Give a name to this Gene Sets Databank.');
        return false;
    } else if (!myForm.species.value) {
        alert('Select the species corresponding to this Gene Sets Databank. If no particular species is used, select this option.');
        return false;
    } else if (!myForm.identifier.value) {
        alert('Select a type of identifiers for this Gene Sets Databank. It can be a custom one.');
        return false;
    }
    return true;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Edit Gene Sets Databank <FONT color='#DD0000'>$gsName</FONT></FONT>
<BR><BR>
<FORM name="editGSForm" method="post" onsubmit="return checkForm(this);" enctype="multipart/form-data">
<INPUT type="hidden" name="ACT" value="$action">
<INPUT type="hidden" name="ID" value="$geneSetsID">
<TABLE border=0>
    <TR>
        <TD bgcolor=$darkColor>
            <TABLE border=0 cellpadding=2>
                <TR>
                    <TH align="right" valign="top" nowrap>Name of the Gene Sets Databank :</TH>
                    <TD bgcolor="$lightColor" nowrap>
                        <INPUT type="text" id="gsName" name="gsName" maxlength="50" value="$gsName">
                    </TD>
                </TR>
                <TR>
                    <TH align="right" valign="top" nowrap>Description :</TH>
                    <TD bgcolor="$lightColor" nowrap>
                        <TEXTAREA id="gsDes" name="gsDes" rows="2" cols="60" maxlength="100">
|;
    print &promsMod::HTMLcompatible($gsDes);
    print qq
|                       </TEXTAREA>
                    </TD>
                </TR>
                <TR>
                    <TH align="right" valign="top" nowrap>Species :</TH>
                    <TD bgcolor="$lightColor" nowrap>
                        <SELECT id="species" name="species">
                            $speciesStrg
                        </SELECT>
                    </TD>
                </TR>
                <TR>
                    <TH align="right" valign="top" nowrap>Identifier type :</TH>
                    <TD bgcolor="$lightColor" nowrap>
                        <SELECT id="identifier" name="identifier">
                            $identStrg
                        </SELECT>
                    </TD>
                </TR>
<!-- Status of the Gene Sets DB
                <TR>
                    <TH align="right" valign="top" nowrap>Status :</TH>
                    <TD bgcolor="$lightColor" nowrap>
                        <SELECT id="status" name="status">
                            $statusStrg
                        </SELECT>
                    </TD>
                </TR>
-->
                <TR>
                    <TD colspan=2 align=center>
                        <INPUT type="submit" name="save" value=" Save ">
                        &nbsp;&nbsp;&nbsp;
                        <INPUT type="reset" value="$resetString">
                        &nbsp;&nbsp;&nbsp;
                        <INPUT type="button" value="Cancel" onclick="cancelAction()">
                    </TD>
                </TR>
            </TABLE>
        </TD>
    </TR>
</TABLE>
</FORM>
</CENTER>
<BR><BR>
<DIV id="divDescription" class="clDescriptionCont"></DIV>
<SCRIPT type="text/javascript">
    setPopup();
</SCRIPT>
</BODY>
</HTML>
|;

}


sub editGeneSetsDB {  # Globals: %promsPath, $dbh, $userID
    my ($geneSetsID) = @_;

    ####################
    ### Display HTML ###
    ####################
    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<TITLE>Editing Gene Sets Databank</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<FONT class="title">Editing Gene Sets Databank...</FONT>
|;

    my $gsName    = param('gsName');
    my $gsDes     = param('gsDes') || '-';
    my $speciesID = (param('species') eq "Custom") ? "NULL" : param('species');
    my $identID   = (param('identifier') eq "Custom")? "NULL" : param('identifier');
    # my $status    = param('status');

    $dbh->do("UPDATE GENESETS SET NAME = \"$gsName\", DES = \"$gsDes\", ID_SPECIES = $speciesID, ID_IDENTIFIER = $identID, UPDATE_DATE = NOW(), UPDATE_USER = \"$userID\" WHERE ID_GENESETS = $geneSetsID") || die $dbh->errstr;
    $dbh->commit;

    print qq
|<BR><BR>
<FONT class="title3">
    Done.
    <BR>
    This page will refresh itself in a few seconds.
</FONT>
<BR>
|;
    sleep 5;
    print qq
|<SCRIPT type="text/javascript">
    window.location = "$promsPath{cgi}/manageGeneSets.cgi?ACT=display";
</SCRIPT>
</BODY>
</HTML>
|;

}


##########################
### Updating Gene Sets ###
##########################
sub printUpdateGeneSetsForm {  # Globals: %promsPath, $dbh, $userID, $action, $lightColor, $darkColor
    my ($geneSetsID) = @_;
    my ($gsName) = $dbh->selectrow_array("SELECT NAME FROM GENESETS WHERE ID_GENESETS = $geneSetsID");
    my $resetString = ($action eq 'add')? ' Clear ' : ' Clear Changes ';

    ####################
    ### Display HTML ###
    ####################
    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<TITLE>Update Gene Sets Databank</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT type="text/javascript">
// Dictionnary of help popup texts (avoid using '"')
const helpText= {
    gmtFile: 'A .gmt file containing your custom gene sets. Click to get to the Broad Institute wiki on GSEA and GMT files for explanations on how a .gmt file is formatted.'
};
|;
    &promsMod::popupInfo();
    print qq
|function cancelAction() {
    window.location = './manageGeneSets.cgi?ACT=display';
}
function checkForm(myForm) {
    if (!myForm.gmtFile.value) {
        alert('Provide a .gmt file for this Gene Sets Databank.');
        return false;
    } else if (!myForm.gsVersionDate.value) {
        alert('Provide the version date for this Gene Sets Databank.');
        return false;
    }
    return true;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Update Gene Sets Databank <FONT color='#DD0000'>$gsName</FONT></FONT>
<BR><BR>
<FORM name="updateGSForm" method="post" onsubmit="return checkForm(this);" enctype="multipart/form-data">
<INPUT type="hidden" name="ACT" value="$action">
<INPUT type="hidden" name="ID" value="$geneSetsID">
<TABLE border=0>
    <TR>
        <TD bgcolor=$darkColor>
            <TABLE border=0 cellpadding=2>
                <TR>
                    <TH align="right" valign="top" nowrap>Gene Sets file (.gmt)<SUP onmouseover="popup(helpText.gmtFile);" onmouseout="popout();"><A href="https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29" target="_blank">?</A></SUP> :</TH>
                    <TD bgcolor="$lightColor" nowrap>
                        <INPUT type="file" name="gmtFile" id="gmtFile" accept=".gmt,.GMT">
                    </TD>
                </TR>
                <TR>
                    <TH align="right" valign="top" nowrap>Gene Sets Version Date</TH>
                    <TD bgcolor="$lightColor" nowrap>
                        <INPUT type="date" id="gsVersionDate" name="gsVersionDate">
                    </TD>
                </TR>
                <TR>
                    <TD colspan=2 align=center>
                        <INPUT type="submit" name="save" value=" Save ">
                        &nbsp;&nbsp;&nbsp;
                        <INPUT type="reset" value="$resetString">
                        &nbsp;&nbsp;&nbsp;
                        <INPUT type="button" value="Cancel" onclick="cancelAction()">
                    </TD>
                </TR>
            </TABLE>
        </TD>
    </TR>
</TABLE>
</FORM>
</CENTER>
<BR><BR>
<DIV id="divDescription" class="clDescriptionCont"></DIV>
<SCRIPT type="text/javascript">
    setPopup();
</SCRIPT>
</BODY>
</HTML>
|;
}


sub updateGeneSetsDB {  # Globals: %promsPath, $dbh, $userID
    my ($geneSetsID) = @_;

    ####################
    ### Display HTML ###
    ####################
    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<TITLE>Updating Gene Sets Databank</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<FONT class="title">Updating Gene Sets Databank...</FONT>
|;

    my ($userGmtFile, $gsFileType, $numSets, $gsVersionDate, $gsStatus);
    $userGmtFile    = (param('gmtFile'))? &promsMod::cleanParameters(param('gmtFile')) : '';
    $gsVersionDate  = param('gsVersionDate') || "CURDATE()";
    $gsStatus       = 1;  # Set Gene Sets in use automatically

    unless ($userGmtFile) {  # Double check, should not pass checkForm()
        $dbh->disconnect;
        die "There was no file provided for the gene sets !";
    }
    ($gsFileType) = $userGmtFile =~ /\.([^.]+)$/;
    $gsFileType = uc($gsFileType);
    my $tmpGmtFile = tmpFileName($userGmtFile);
    ($numSets) = (`wc -l $tmpGmtFile` =~ /^(\d+)\s/);

    my ($previousGmt, $previousVersion, $otherVersions, $newGmt, $gmtPath);
    ($previousGmt, $previousVersion, $otherVersions) = $dbh->selectrow_array("SELECT GENESETS_FILE, VERSION_DATE, OTHER_VERSIONS FROM GENESETS WHERE ID_GENESETS = $geneSetsID");
    $gmtPath = dirname($previousGmt);
    $newGmt = "$gmtPath/$userGmtFile";
    $otherVersions = ($otherVersions) ? $otherVersions.';'.$previousVersion : "$previousVersion";

    # Check that archive exists for previous version and remove existing .gmt file if ok
    $previousVersion =~ s/-//g;
    my $previousArchive = "$promsPath{'genesets'}/$gmtPath/GS$geneSetsID\_$previousVersion.tar.gz";

    if (-s $previousArchive) {
        unlink "$promsPath{'genesets'}/$previousGmt" if (-f "$promsPath{'genesets'}/$previousGmt");
    } else {
        &createArchive($geneSetsID);
        if (-s $previousArchive) {  # Double check creation of archive
            unlink "$promsPath{'genesets'}/$previousGmt" if (-f "$promsPath{'genesets'}/$previousGmt");
        } else {
            $dbh->do("UPDATE GENESETS SET GENESETS_STATUS = -1, UPDATE_DATE = NOW(), UPDATE_USER = \"$userID\" WHERE ID_GENESETS = $geneSetsID") || die $dbh->errstr;  # 'Error' status
            $dbh->commit;
            $dbh->disconnect;
            die "Unable to create archive for file $previousGmt";
        }
    }

    # Do update of file and version in DB
    move($tmpGmtFile, "$promsPath{'genesets'}/$newGmt");
    $dbh->do("UPDATE GENESETS SET GENESETS_FILE = \"$newGmt\", GENESETS_TYPE = \"$gsFileType\", NUM_SETS = $numSets, GENESETS_STATUS = $gsStatus, VERSION_DATE = \"$gsVersionDate\", OTHER_VERSIONS = \"$otherVersions\", UPDATE_DATE = NOW(), UPDATE_USER = \"$userID\" WHERE ID_GENESETS = $geneSetsID") || die $dbh->errstr;
    $dbh->commit;

    print qq
|<BR><BR>
<FONT class="title3">
    Done.
    <BR>
    Creating archive for this version of the Gene Sets databank.
</FONT>
|;
    # Create archive for this Gene Sets DB version
    &createArchive($geneSetsID);

    print qq
|<BR><BR>
<FONT class="title3">
    Done.
    <BR>
    This page will refresh itself in a few seconds.
</FONT>
<BR>
|;
    sleep 5;
    print qq
|<SCRIPT type="text/javascript">
    window.location = "$promsPath{cgi}/manageGeneSets.cgi?ACT=display";
</SCRIPT>
</BODY>
</HTML>
|;

}


#########################################
### Restoring Other Gene Sets Version ###
#########################################
sub printRestoreGeneSetsForm {  # Globals: %promsPath, $dbh, $action, $lightColor, $darkColor
    my ($geneSetsID) = @_;
    my ($gsName, $currentVersion, $otherVersions) = $dbh->selectrow_array("SELECT NAME, VERSION_DATE, OTHER_VERSIONS FROM GENESETS WHERE ID_GENESETS = $geneSetsID");
    my $resetString = ($action eq 'add')? ' Clear ' : ' Clear Changes ';

    my $otherVersionsStrg = "<OPTION value=\"\">-= Select =-</OPTION>\n";
    foreach my $otherVersion (sort{$a cmp $b} split(';', $otherVersions)) {
        $otherVersionsStrg .= "<OPTION value=\"$otherVersion\">$otherVersion</OPTION>\n";
    }

    ####################
    ### Display HTML ###
    ####################
    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<TITLE>Restore Gene Sets Databank</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT type="text/javascript">
function cancelAction() {
    window.location = './manageGeneSets.cgi?ACT=display';
}
function checkForm(myForm) {
    if (!myForm.restoreVersion.value) {
        alert('Select a version of the Gene Sets Databank to restore !');
        return false;
    }
    return true;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Restore other version of Gene Sets Databank <FONT color='#DD0000'>$gsName</FONT></FONT>
<BR><BR>
<FORM name="restoreGSForm" method="post" onsubmit="return checkForm(this);" enctype="multipart/form-data">
<INPUT type="hidden" name="ACT" value="$action">
<INPUT type="hidden" name="ID" value="$geneSetsID">
<TABLE border=0>
    <TR>
        <TD bgcolor=$darkColor>
            <TABLE border=0 cellpadding=2>
                <TR>
                    <TH align="right" valign="top" nowrap>Current version date</TH>
                    <TD bgcolor="$lightColor" nowrap>&nbsp;$currentVersion</TD>
                </TR>
                <TR>
                    <TH align="right" valign="top" nowrap>Version to restore :</TH>
                    <TD bgcolor="$lightColor" nowrap>
                        <SELECT id="restoreVersion" name="restoreVersion">
                            $otherVersionsStrg
                        </SELECT>
                    </TD>
                </TR>
                <TR>
                    <TD colspan=2 align=center>
                        <INPUT type="submit" name="save" value=" Save ">
                        &nbsp;&nbsp;&nbsp;
                        <INPUT type="reset" value="$resetString">
                        &nbsp;&nbsp;&nbsp;
                        <INPUT type="button" value="Cancel" onclick="cancelAction()">
                    </TD>
                </TR>
            </TABLE>
        </TD>
    </TR>
</TABLE>
</FORM>
</CENTER>
<BR><BR>
<DIV id="divDescription" class="clDescriptionCont"></DIV>
<SCRIPT type="text/javascript">
    setPopup();
</SCRIPT>
</BODY>
</HTML>
|;
}


sub restoreGeneSetsDb {
    my ($geneSetsID) = @_;
    my $restoreVersion = param('restoreVersion');

    ####################
    ### Display HTML ###
    ####################
    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<TITLE>Restoring Gene Sets Databank</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<FONT class="title">Restoring Gene Sets Databank...</FONT>
|;

    # Check that archive exists for previous version and remove existing .gmt file if ok
    # Otherwise, create archive before removing file
    my ($previousGmt, $previousVersion, $otherVersions) = $dbh->selectrow_array("SELECT GENESETS_FILE, VERSION_DATE, OTHER_VERSIONS FROM GENESETS WHERE ID_GENESETS = $geneSetsID");
    my $gmtPath = dirname($previousGmt);
    $previousVersion =~ s/-//g;
    my $previousArchive = "$promsPath{'genesets'}/$gmtPath/GS$geneSetsID\_$previousVersion.tar.gz";

    if (-s $previousArchive) {
        unlink "$promsPath{'genesets'}/$previousGmt" if (-f "$promsPath{'genesets'}/$previousGmt");
    } else {
        &createArchive($geneSetsID);
        if (-s $previousArchive) {  # Double check creation of archive
            unlink "$promsPath{'genesets'}/$previousGmt" if (-f "$promsPath{'genesets'}/$previousGmt");
        } else {
            $dbh->do("UPDATE GENESETS SET GENESETS_STATUS = -1, UPDATE_DATE = NOW(), UPDATE_USER = \"$userID\" WHERE ID_GENESETS = $geneSetsID") || die $dbh->errstr;  # 'Error' status
            $dbh->commit;
            $dbh->disconnect;
            die "Unable to create archive for file $previousGmt";
        }
    }
    # Restore previous archive and extract .gmt file + update DB
    &restoreArchive($geneSetsID, $restoreVersion);

    print qq
|<BR><BR>
<FONT class="title3">
    Done.
</FONT>
<BR><BR>
<FONT class="title3">
    Done.
    <BR>
    This page will refresh itself in a few seconds.
</FONT>
<BR>
|;
    sleep 5;
    print qq
|<SCRIPT type="text/javascript">
    window.location = "$promsPath{cgi}/manageGeneSets.cgi?ACT=display";
</SCRIPT>
</BODY>
</HTML>
|;


}


######################################
### Deleting / Archiving Gene Sets ###
######################################
sub deleteGeneSetsDB {
    my ($geneSetsID) = @_;
    my ($gsName) = $dbh->selectrow_array("SELECT NAME FROM GENESETS WHERE ID_GENESETS = $geneSetsID");

    ####################
    ### Display HTML ###
    ####################
    print header(-'content-encoding'=>'no',-charset=>'UTF-8');
    warningsToBrowser(1);
    print qq
|<HTML>
<HEAD>
<TITLE>Deleting Gene Sets Databank</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<FONT class="title">Deleting Gene Sets Databank...</FONT>
|;

    # Check if db was used by a valid GSEA
    my %gseaUsingDB;
    my $isUsed = $dbh->selectrow_array("SELECT 1 FROM PATHWAY_ANALYSIS WHERE ANALYSIS_TYPE = 'GSEA' AND STATUS IN (0, 1) AND PARAM_STRG LIKE \"\%geneSetsDB=#$geneSetsID\%\" LIMIT 1");

    if ($isUsed) {  # Cannot delete the Gene Sets DB if used by GSEA
        print qq
|<BR><BR><BR>
<FONT class="title3">Cannot delete <FONT color='#DD0000'>$gsName</FONT></FONT> because it is linked to some valid Gene Set Enrichment Analyses that use this databank.
<BR><BR><BR>
<FONT class="title3">Archiving the databank and setting status <FONT color='#DD0000'>Not in use</FONT> instead</FONT>
|;
        # Check that archive exists for current version and remove existing .gmt file if ok
        my ($currentGmt, $currentVersion) = $dbh->selectrow_array("SELECT GENESETS_FILE, VERSION_DATE FROM GENESETS WHERE ID_GENESETS = $geneSetsID");
        $currentVersion =~ s/-//g;
        my $gmtPath = dirname($currentGmt);
        my $currentArchive = "$promsPath{'genesets'}/$gmtPath/GS$geneSetsID\_$currentVersion.tar.gz";

        if (-s $currentArchive) {
            unlink "$promsPath{'genesets'}/$currentGmt" if (-f "$promsPath{'genesets'}/$currentGmt");
        } else {
            &createArchive($geneSetsID);
            if (-s $currentArchive) {  # Double check creation of archive
                unlink "$promsPath{'genesets'}/$currentGmt" if (-f "$promsPath{'genesets'}/$currentGmt");
            } else {
                $dbh->do("UPDATE GENESETS SET GENESETS_STATUS = -1, UPDATE_DATE = NOW(), UPDATE_USER = \"$userID\" WHERE ID_GENESETS = $geneSetsID") || die $dbh->errstr;  # 'Error' status
                $dbh->commit;
                $dbh->disconnect;
                die "Unable to create archive for file $currentGmt";
            }
        }
        $dbh->do("UPDATE GENESETS SET GENESETS_STATUS = 0, UPDATE_DATE = NOW(), UPDATE_USER = \"$userID\" WHERE ID_GENESETS = $geneSetsID") || die $dbh->errstr;  # 'Not in use'

    } else {  # No GSEA launched with this Gene Sets DB
        my ($currentGmt, $currentVersion) = $dbh->selectrow_array("SELECT GENESETS_FILE, VERSION_DATE FROM GENESETS WHERE ID_GENESETS = $geneSetsID");
        my $gmtPath = dirname($currentGmt);
        $dbh->do("DELETE FROM GENESETS WHERE ID_GENESETS = $geneSetsID") || die $dbh->errstr;
        # Delete corresponding files and archives
        system "rm -r $promsPath{'genesets'}/$gmtPath" if -e "$promsPath{'genesets'}/$gmtPath";
        if ($gmtPath =~ /private/) {  # DB was specific to a project -> remove project dir if empty
            my $projDir = dirname($gmtPath);
            if (-d "$promsPath{'genesets'}/$projDir" && &isEmptyDir("$promsPath{'genesets'}/$projDir")) {
                system "rm -r $promsPath{'genesets'}/$projDir";
            }
        }
    }
    $dbh->commit;

    sleep 3;
    print qq
|<BR><BR><BR>
Done.
<BR><BR><BR>
Redirecting...
|;
    sleep 3;
    print qq
|<SCRIPT type="text/javascript">
    window.location = "$promsPath{cgi}/manageGeneSets.cgi?ACT=display";
</SCRIPT>
</BODY>
</HTML>
|;

}


####################
### Archive Subs ###
####################
sub createArchive {  # Globals: %promsPath, $dbh
    my ($geneSetsID) = @_;
    my ($currentGmt, $currentVersion) = $dbh->selectrow_array("SELECT GENESETS_FILE, VERSION_DATE FROM GENESETS WHERE ID_GENESETS = $geneSetsID");
    $currentVersion =~ s/-//g;

    my $gmtFile = basename($currentGmt);
    my $gmtPath = dirname($currentGmt);
    my $gsDir = "$promsPath{'genesets'}/$gmtPath";

    # Add .arch extension to the archived file name to recover it easily if needed
    my $archiveFile = "$gmtFile.arch";
    copy("$gsDir/$gmtFile", "$gsDir/$archiveFile") || die "Creation of a copy to archive failed: $!";

    my $archiveName = "GS$geneSetsID\_$currentVersion";
    system "cd $gsDir; tar -zcf $archiveName.tar.gz $archiveFile";
    unlink "$gsDir/$archiveFile";
}


sub restoreArchive {  # Globals: %promsPath, $dbh, $userID
    my ($geneSetsID, $newVersion) = @_;
    my $gsStatus = 1;  # Set Gene Sets in use automatically
    my $restoreVersion = $newVersion;
    $restoreVersion =~ s/-//g;

    my ($currentGmt, $currentVersion, $otherVersions) = $dbh->selectrow_array("SELECT GENESETS_FILE, VERSION_DATE, OTHER_VERSIONS FROM GENESETS WHERE ID_GENESETS = $geneSetsID");

    my $gmtPath = dirname($currentGmt);
    my $gsDir = "$promsPath{'genesets'}/$gmtPath";
    my $archiveName = "GS$geneSetsID\_$restoreVersion";
    system "cd $gsDir; tar -zxf $archiveName.tar.gz";

    # Retrieve extracted file and update DB
    opendir(DIR, "$gsDir");
    my @archFiles = grep(/\.arch$/, readdir(DIR));
    closedir(DIR);

    my $archFileNb = scalar @archFiles;
    unless ($archFileNb == 1) {  # We assume that there is no .arch file in $gsDir at first
        $dbh->do("UPDATE GENESETS SET GENESETS_STATUS = -1, UPDATE_DATE = NOW(), UPDATE_USER = \"$userID\" WHERE ID_GENESETS = $geneSetsID") || die $dbh->errstr;  # 'Error' status
        $dbh->commit;
        $dbh->disconnect;
        die "Cannot update DB after restoring archive. A single .arch file is expected in directory $gsDir, found $archFileNb.";
    }

    my $archiveFile = $gmtPath . "/" . $archFiles[0];
    (my $newGmt = $archiveFile) =~ s/\.arch$//;  # Remove artificial ext .arch introduced at archive creation
    copy("$promsPath{'genesets'}/$archiveFile", "$promsPath{'genesets'}/$newGmt") || die "Creation of a copy of archived file failed: $!";
    unlink "$promsPath{'genesets'}/$archiveFile";

    my ($gsFileType) = $newGmt =~ /\.([^.]+)$/;
    $gsFileType = uc($gsFileType);
    my ($numSets) = (`wc -l $promsPath{'genesets'}/$newGmt` =~ /^(\d+)\s/);

    my @allVersions = split(';', $otherVersions);
    push @allVersions, $currentVersion;
    @allVersions = sort @allVersions;
    my $newOtherVersions = join(';', grep(!/$newVersion/, @allVersions));

    $dbh->do("UPDATE GENESETS SET GENESETS_FILE = \"$newGmt\", GENESETS_TYPE = \"$gsFileType\", NUM_SETS = $numSets, GENESETS_STATUS = $gsStatus, VERSION_DATE = \"$newVersion\", OTHER_VERSIONS = \"$newOtherVersions\", UPDATE_DATE = NOW(), UPDATE_USER = \"$userID\" WHERE ID_GENESETS = $geneSetsID");
    $dbh->commit;
}


###################
### Utility sub ###
###################
sub isEmptyDir {  # Sub to check if a directory is empty
    my ($dirName) = @_;
    opendir(DIR, $dirName) || die "Error in opening dir $dirName\n";
    return scalar(grep { $_ ne "." && $_ ne ".." } readdir(DIR)) == 0;
}


####>Revision history<####
# 1.0.1 [MODIF] Removed possibility to create gene sets databank for non-reference species (VL 07/06/21)
# 1.0.0 Created (VL 10/11/20)
