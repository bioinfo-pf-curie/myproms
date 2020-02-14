#!/usr/local/bin/perl -w

################################################################################
# exportLibraryProteins.cgi         1.0.0                                      #
# Authors: V. Laigle (Institut Curie)                                          #
# Contact: myproms@curie.fr                                                    #
# Export a protein list from spectral library                                  #
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

$|=1;       # buffer flush (allows instant printing)
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use strict;
use promsConfig;
use promsMod;

my %promsPath=&promsConfig::getServerInfo;
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my $userID=$ENV{'REMOTE_USER'};

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

my $libID = param('ID');
my ($libraryName, $libVersion, $identifierType, $organism) = $dbh->selectrow_array("SELECT NAME, VERSION_NAME, IDENTIFIER_TYPE, ORGANISM FROM SWATH_LIB WHERE ID_SWATH_LIB='$libID' ") or die "Couldn't fetch library info: " . $dbh->errstr;
my $action=(param('submit'))? param('submit') : "";

if ($action) {
    my $format = &promsMod::cleanParameters(param("exportFormat"));

    my $workingDir      = "$promsPath{swath_lib}/SwLib_$libID";
    my $libFile         = "$workingDir/$libraryName.sptxt";
    my $downloadFile    = "prots_$libraryName.$format";

    my %proteins;  # To store all proteins from the library with their info
    my %proteinsNotInProms;  # In case some proteins are not yet recorded

    # Get all proteins from library file (including redundants)
    my @redundantProts;
    open (LIB_FILE, "<", $libFile) or die "Can't open library file $libFile : $!";
    while (my $line = <LIB_FILE>) {
        if ($line =~ /Protein=(\d+)\/(\S+)\s/ && $line !~ /Protein=(\d+)\/(?:sp|tr)\|Biognosys\|iRT/) {
            push @redundantProts, (split('/', $2));
        }
    }
    close LIB_FILE;

    # Prepare queries according to identifier type
    my $libPattern;  # Pattern of the protein name in the given library
    my $identifierCode;
    if ($identifierType=~/ID/) {
        $identifierCode = 'ID';
        $libPattern = qr/^(.*)$/;
    } elsif ($identifierType=~/ACCESSION/) {
        $identifierCode = 'AC';
        $libPattern = qr/^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$/;
    } elsif ($identifierType=~/ALL/) {
        $identifierCode = 'AC';
        $libPattern = qr/^(?:sp|tr)\|([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\|[\w\d]+$/;
    }

    my ($identifierID) = $dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='$identifierCode'") or die "Couldn't prepare statement: " . $dbh->errstr;
    my ($geneIdentifierID) = $dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'") or die "Couldn't prepare statement: " . $dbh->errstr;
    
    # Map species ID to organism name
    my %organisms;
    if ($organism) {
        my ($speciesID) = $dbh->selectrow_array("SELECT ID_SPECIES FROM SPECIES WHERE (SCIENTIFIC_NAME='$organism' OR COMMON_NAME='$organism') LIMIT 0,1") or die "Couldn't select species name: " . $dbh->errstr;
        $organisms{$speciesID} = $organism;
    } else {
        my $sthOrganism = $dbh->prepare("SELECT ID_SPECIES, SCIENTIFIC_NAME, COMMON_NAME FROM SPECIES") or die "Couldn't prepare statement: " . $dbh->errstr;
        $sthOrganism->execute;
        while (my ($speciesID, $scientificName, $commonName) = $sthOrganism->fetchrow_array) {
            $organisms{$speciesID} = ($scientificName)? $scientificName : $commonName;
        }
        $sthOrganism->finish;
    }

    # Remove redundants
    my %uniqueProts;
    foreach my $prot (@redundantProts) {
        my ($protIdentifier) = $prot =~ /$libPattern/;
        if ($protIdentifier) {
            if ($uniqueProts{$protIdentifier}) {
                next;
            } else {
                $uniqueProts{$protIdentifier} = {};
            }
        }
    }
    my $protString = join(',', map(qq('$_'), keys %uniqueProts));

    # Fetch proteins info from DB
    my $sthProtInfo = $dbh->prepare(
        "SELECT MI5.ID_MASTER_PROTEIN, GROUP_CONCAT(DISTINCT MI5.VALUE SEPARATOR ';'), 
        MP.ID_SPECIES, MP.PROT_DES, GROUP_CONCAT(DISTINCT MI6.VALUE SEPARATOR ';') 
        FROM MASTER_PROTEIN MP 
        INNER JOIN MASTERPROT_IDENTIFIER MI5 ON MP.ID_MASTER_PROTEIN=MI5.ID_MASTER_PROTEIN 
        INNER JOIN MASTERPROT_IDENTIFIER MI6 ON MI5.ID_MASTER_PROTEIN=MI6.ID_MASTER_PROTEIN 
        WHERE MI5.ID_IDENTIFIER=$identifierID AND MI5.VALUE IN ($protString) AND MI6.ID_IDENTIFIER=$geneIdentifierID 
        GROUP BY MI5.ID_MASTER_PROTEIN"
    );
    $sthProtInfo->execute();
    while (my ($masterProtID, $protIdentStrg, $speciesID, $desc, $geneNames) = $sthProtInfo->fetchrow_array) {
        if ($format eq 'csv') {
            $desc =~ s/,/./g;  # commas can be present in description, e.g. in molecules names, which breaks csv column
        }
        $proteins{$masterProtID} = {
            'protIdentifiers'   => $protIdentStrg,
            'organism'          => $organisms{$speciesID},
            'description'       => $desc,
            'geneNames'         => $geneNames
        }
    }
    $sthProtInfo->finish;
    $dbh->disconnect;

    # Check for proteins not present in MyProMS DB
    my @protsInProms;
    foreach my $protID (keys %proteins) {
        push @protsInProms, (split(';', $proteins{$protID}{'protIdentifiers'}));
    }
    
    foreach my $protIdentifier (keys %uniqueProts) {
        unless (grep (/$protIdentifier/, @protsInProms)) {
            $proteinsNotInProms{$protIdentifier} = 1;  # = {} to get external annotation
        }
    }

    # Exporting proteins info
    my $sep = ($format eq 'csv')? "," : "\t";  
	
    print header(-type=>"text/plain", -attachment=>$downloadFile);

    my $fileHeader = join($sep, ("Protein_$identifierCode", "Gene_names", "Organism", "Description"));
    print STDOUT "$fileHeader\n";
    foreach my $protID (sort {$a <=> $b} keys %proteins) {
        my $protLine = join($sep, ($proteins{$protID}->{'protIdentifiers'}, 
                                    $proteins{$protID}->{'geneNames'}, 
                                    $proteins{$protID}->{'organism'}, 
                                    $proteins{$protID}->{'description'}
                                    ));
        print STDOUT "$protLine\n"; 
    }
    foreach my $protIdentifier (keys %proteinsNotInProms) {
        print STDOUT "$protIdentifier\n";
    }
} else {
    #######################
    ####>Starting HTML<####
    #######################
    print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);

    # Form to select files
    print qq
|<html>
    <head>
        <title>Export Proteins</title>
        <link rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
        <style type="text/css">
            .popup {
                background-color : #FFFFFF;
                border : solid 3px #999999;
                padding : 5px;
                box-shadow : 10px 10px 20px #808080;
                position : absolute;
                display : none;
            } 
        </style>
        <script language="JavaScript">
|;

    &promsMod::popupInfo();
    print qq 
|       </script>
    </head>
|;
    ###> Print HTML form <###
    print qq 
|   <body background="$promsPath{images}/bgProMS.gif">
        <div style="float:top">
        <br>
        <table>
            <tr>
                <th bgcolor="$darkColor">
                <font class="title2">&nbsp;Go to:</FONT>
                <select class="title action" onchange="if (this.value) {window.location='$promsPath{cgi}/'+this.value;}">
                    <option value="">-= Select =-</option>
                    <option value="promsMain.cgi">Main Window</option>
                    <option value="selectProject.cgi">Project Selection</option>
                </select>
                </th>
            </tr>
        </table>
        </div>

        <center>
        
        <div id="form" style="display:" >
        <font class="title1">Export proteins from <font class="title1" color="#DD0000">$libraryName</font>
        </font>
        <br><br><br><br>
        <form method="POST" action="./exportLibraryProteins.cgi" name="export">
        <table bgcolor="$darkColor">
            <input type="hidden" name="ID" id="ID" value="$libID">
            <tr>
                <th align=right valign="top">Export format : </th>
                <td bgcolor="$lightColor" width=400px>
                    <select name="exportFormat" required>
                        <option value="">-= Select format =-</option>
                        <option value="csv">Comma-Separated Values (.csv)</option>
                        <option value="tsv">Tab-Separated Values (.tsv)</option>
                    </select><br>
                </td>
            </tr>
            <tr>
                <th colspan=2>
                <br>
                    <!-- SUBMIT button -->
                    <input type="submit" name="submit" value="Download">
                    <!-- CLEAR button -->
                    &nbsp &nbsp &nbsp
                    <input type="reset" value="Clear" />
                    <!-- CANCEL button -->
                    &nbsp &nbsp &nbsp
                    <input type="button" value="Cancel" onclick="window.location='./listSwathLibraries.cgi'">
                </th>
            </tr>
        </table>
        </form>
        </div>
        <br><br><br>
        <input type="button" class="buttonadd" value="Return to spectral libraries list." onclick="window.location='./listSwathLibraries.cgi'">
 
        <script type="text/javascript">
            setPopup()
        </script>
        
        </center>
    </body>
</html>
|;
}

$dbh->disconnect;

####>Revision history<####
# 1.0.0 Created to allow to export proteins from swath libraries (VL 16/09/2019)
