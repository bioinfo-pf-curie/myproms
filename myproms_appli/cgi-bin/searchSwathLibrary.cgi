#!/usr/local/bin/perl -w

################################################################################
# searchSwathLibrary.cgi         1.1.0                                         #
# Authors: M. Le Picard (Institut Curie)                                       #
# Contact: myproms@curie.fr                                                    #
#Search peptide in the libraries available in myProMS                          #
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
use IO::Handle;
use promsConfig;
use File::Copy;
use File::Basename;
use XML::Simple;
use promsMod;
use LWP::UserAgent;
use String::Util qw(trim);
use File::Compare;
use promsMod;

my %promsPath=&promsConfig::getServerInfo;
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my $bgColor=$lightColor;


##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;
my $userID=$ENV{'REMOTE_USER'};


#######################
####>Starting HTML<####
#######################
print header(-'content-encoding'=>'no',-charset=>'UTF-8'); warningsToBrowser(1);

#############################
####>Fetching parameters<####
#############################
##Library informations
my $libraryID=param('ID');
my ($selectLibraryName,$specie,$split,$identLibType)=$dbh->selectrow_array("SELECT NAME, ORGANISM, SPLIT, IDENTIFIER_TYPE FROM SWATH_LIB WHERE ID_SWATH_LIB='$libraryID' ") or die "Couldn't prepare statement: " . $dbh->errstr;

my %databank;
my $sthDB=$dbh->prepare("SELECT DS.ID_DATABANK, FASTA_FILE FROM DATABANK_SWATHLIB DS, DATABANK D WHERE ID_SWATH_LIB=? AND DS.ID_DATABANK=D.ID_DATABANK ") or die "Couldn't prepare statement: " . $dbh->errstr;
$sthDB->execute($libraryID);
while (my ($dbID,$fastaName)=$sthDB->fetchrow_array) {
    $databank{$dbID}=$fastaName;
}
$sthDB->finish;

my $specieLatinName;
if ($specie){$specieLatinName=$dbh->selectrow_array("SELECT SCIENTIFIC_NAME FROM SPECIES WHERE SCIENTIFIC_NAME='$specie' OR COMMON_NAME='$specie'");}
else{$specieLatinName=param('speciesList');}

##Library path
my $spcLibFormat = $dbh->selectrow_array("SELECT IF(PARAM_STRG LIKE '%Spectronaut%', 'SPC', 'TPP') FROM SWATH_LIB WHERE ID_SWATH_LIB=$libraryID");
my $searchFile = ($spcLibFormat eq 'TPP') ? "$promsPath{swath_lib}/SwLib_$libraryID/$selectLibraryName.sptxt" : "$promsPath{swath_lib}/SwLib_$libraryID/$selectLibraryName.tsv";


if (param('selectprot')){&searchLibrary(param('selectprot'),param('accessionList'),param('sequenceProt'));exit;}



################################
####>Form to search protein<####
################################

print qq
|<HTML>
<HEAD>
<TITLE>Search</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
.drawViewPept:hover {
    text-decoration: underline;
    cursor: pointer;
}

.row1{background-color:$lightColor;}
.row2{background-color:$darkColor;}
.overlapRes0 {color:black;}
.overlapRes1 {font-weight:bold;color:blue;}
.overlapRes2 {font-weight:bold;color:red;}
.padtab{padding-left:10px;padding-right:10px;}
.popup {background-color:#FFFFFF;border:solid 3px #999999;padding:5px;box-shadow:10px 10px 20px #808080;position:absolute;display:none;} 
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
    |;
    &promsMod::popupInfo();
    print qq |
    // AJAX --->
    function getXMLHTTP() {
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
        
    function ajaxSearchLibrary(selectprot,accessionList,sequenceProt){
        document.getElementById('waitDIV2').style.display='';
        document.getElementById('result').innerHTML="";
        var ID=document.getElementById('ID').value;
        //Creation of the XMLHTTPRequest object
        var XHR = getXMLHTTP();
        XHR.open("GET","searchSwathLibrary.cgi?selectprot="+selectprot+"&ID="+ID+"&accessionList="+accessionList+"&sequenceProt="+sequenceProt,true);
        XHR.onreadystatechange=function() {
            if (XHR.readyState==4 && XHR.responseText) {
               document.getElementById('result').innerHTML=XHR.responseText;
               document.getElementById('waitDIV2').style.display='none';
            }
        }
        
        XHR.send(null);
    }
    // <--- AJAX

    function showHideFunction(id) {
        if(document.getElementById(id+'DIV').style.display=="none"){
            document.getElementById(id+'BUTTON').value='Hide';
            document.getElementById(id+'DIV').style.display='';
        }
        else{
            document.getElementById(id+'BUTTON').value='Show';
            document.getElementById(id+'DIV').style.display='none';
        }
    }
    
    function protListView(pepnumber){
        var object=document.getElementById('DIVprotList'+pepnumber);
        if(object.style.display=='')
            object.style.display='none';
        else
            object.style.display='';
    }
    
    function pepView(pepId,sequence,irt,massObs,massExp,libID) {
        if (pepId != selectedPep \|\| spectWin.closed) {
            var call='lib';
            var typefile= ('$spcLibFormat' == 'TPP') ? 'sptxt' : 'tsv';
            var file='$searchFile';
            var irt=irt;
            selectPeptide(pepId);
            var paramString="file="+file+"&CALL="+call+"&TYPE="+typefile+"&SEQUENCE="+encodeURIComponent(sequence)+"&irt="+irt+"&massObs="+massObs+"&massExp="+massExp+"&libID="+libID;
            spectWin=window.open("$promsPath{cgi}/peptide_view.cgi?"+paramString,'SpectrumWindow','width=950,height=950,location=no,resizable=yes,scrollbars=yes');
        }
        spectWin.focus();
    }
    function selectPeptide(newPep) {
        if (selectedPep && document.getElementById(selectedPep)) { 
            document.getElementById(selectedPep).style.color='#000000';
        }
        document.getElementById(newPep).style.color='#DD0000'; //style.background = '#DD0000';
        selectedPep=newPep;
    }
    var selectedPep;
    

</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
    <DIV style="float:top">
        <BR>
        <TABLE><TR><TH bgcolor="$darkColor">
            <FONT class="title2">&nbsp;Go to:</FONT>
            <SELECT class="title action" onchange="if (this.value) {window.location='$promsPath{cgi}/'+this.value;}">
                <OPTION value="">-= Select =-</OPTION>
                <OPTION value="promsMain.cgi">Main Window</OPTION>
                <OPTION value="selectProject.cgi">Project Selection</OPTION>
            </SELECT>
        </TH></TR></TABLE>
    </DIV>
    <CENTER>
        <FONT class="title1">Search in <FONT class="title1" color="#DD0000">$selectLibraryName</FONT></FONT><BR><BR><BR><BR>
        <FORM method="POST" action="./searchSwathLibrary.cgi" name="parameters" enctype="multipart/form-data">
            <TABLE bgcolor="$darkColor"><TR><TH align=right valign=top>Entry : </TH><TD bgcolor="$lightColor"><TEXTAREA rows="2" cols="43" name='search' id='search' required>|; print param('search') if param('search'); print qq |</TEXTAREA></TD></TR>
            <INPUT type="hidden" name="ID" id="ID" value="$libraryID">|;
            if ($specie && $identLibType=~/UNIPROT|UniProt/){
                print "<TR><TH align=right valign=top>Species name : </TH><TD  bgcolor=\"$lightColor\"><INPUT type=\"checkbox\" name=\"specie\" checked><I>$specie</I></TD></TR>";
            }    
            elsif($identLibType=~/UNIPROT|UniProt/){    
                my $sthSpeciesList=$dbh->prepare("SELECT SCIENTIFIC_NAME FROM SPECIES") or die "Couldn't prepare statement: " . $dbh->errstr;
                    $sthSpeciesList->execute;
                    print "<TR><TH valign=top align=right>Species name : </TH><TD bgcolor=\"$lightColor\"><SELECT name=\"speciesList\" id=\"speciesList\"><option value=\"\">-= Select Species =-</option>";
                    
                    while (my($specieLatinName)=$sthSpeciesList->fetchrow_array){
                        print qq
                        |<option value="$specieLatinName"><I>$specieLatinName</I></option>|;
                    }
                    $sthSpeciesList->finish;
                print "</SELECT></TD></TR>";
            }
            

            print qq
            |<TR ><TH colspan=2><BR><INPUT type="submit" name="submit" value="Submit">
            <!-- CLEAR button -->
            &nbsp &nbsp &nbsp<INPUT type="reset" value="Clear" />
            <!-- CANCEL button -->
            &nbsp &nbsp &nbsp<INPUT type="button" value="Cancel" onclick="window.location='./listSwathLibraries.cgi'"></TH></TR>
            
            </TABLE>
        </FORM>
|;
        
        
###################################
####>List results from Uniprot<####
###################################
if (param('submit')){
    
    my $text=($identLibType=~/UNIPROT|UniProt/)? 'Uniprot' : $selectLibraryName ;
    print "<DIV id=\"waitDIV\"><BR><BR><BR><BR><BR><FONT class=\"title3\">Fetching data from $text ...</FONT><BR><BR><IMG src=\"$promsPath{images}/scrollbarGreen.gif\"><BR><BR></DIV>";
    
    
    $specieLatinName=~s/\(.*\)// if ($specieLatinName && $specieLatinName=~/\(.*\)/);
    
    my $checkSpecie=param('specie');
    my $trClass='row1';
    my @searchText=split('[;,\s]',param('search')) if param('search');
    @searchText=&promsMod::cleanParameters(@searchText);
    my $nbSearch=0;
    foreach my $search (@searchText){
        next if $search eq '';
        $search=trim($search);
        
        if ($identLibType=~/UNIPROT|UniProt/) {
            ###> Search in Uniprot
            my $params;
            if ($specieLatinName eq '' || $checkSpecie eq '' ) {
                $params = {
                    query => "$search",
                    sort=>'score',
                    format=>'tab',
                    columns=>'id,entry name,protein names,length,genes',
                    limit=>10
                }; 
            }
            else{
                $params = {
                    query => "$search organism:($specieLatinName)", # reviewed:yes name: and/or gene:}
                    sort=>'score',
                    format=>'tab',
                    columns=>'id,entry name,protein names,length,genes',
                    limit=>10
                };
            }   

            my $agent = LWP::UserAgent->new(); #env_proxy => 1
            $agent->proxy('http','http://www-cache.curie.fr:3128/');
            push @{$agent->requests_redirectable}, 'POST';
        
            my $responseSearchProt = $agent->post("http://www.uniprot.org/uniprot/",$params);
            
            while (my $wait = $responseSearchProt->header('Retry-After')) {
                sleep $wait;
                $responseSearchProt = $agent->get($responseSearchProt->base);
            }
            my $lineCount=0;
            
            my %uniprotResultList;
            ###> Post in a table all results
            if ($responseSearchProt->is_success) {
                ###>recovering all accession numbers 
                my (@tabUniAC,@result);  
                foreach my $line (split(/\n/,$responseSearchProt->content)) {
                    $lineCount++;
                    next if $lineCount==1;
                    @result=split(/\t/,$line);
                    push(@tabUniAC,$result[0]);
                }
                
               
                ###> recovering accession numbers and sequence of each protein found by uniprot
                my $acString = join (',',@tabUniAC);
                my $URL="http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&format=uniprot&style=raw&id=$acString";
                my $response = $agent->get($URL);
                while (my $wait = $response->header('Retry-After')) {
                    sleep $wait;
                    $response = $agent->get($response->base);
                }
                my (@uniACList,$firstAC,$seqProt);       #@uniACList = list of all AC of one protein
                my $aaNumber=1;
                $seqProt='';
                my %accessionSeqList;
                if ($response->is_success) {
                    foreach my $line (split(/\n/,$response->content)) {
                        ###> insert sequence and AC numbers in hash table 
                        if (substr($line,0,2) eq '//') {
                            @{$accessionSeqList{$firstAC}}=(@uniACList,"$seqProt");
                            @uniACList=();
                            $seqProt='';
                        }
                        elsif (substr($line,0,2) eq 'AC') {
                            if (scalar @uniACList ==0){
                                @uniACList=split(";",substr($line,3,-1));
                                $firstAC=trim($uniACList[0]);
                            }
                            else{push @uniACList,split(";",substr($line,3,-1));} 
                        }
                        elsif (substr($line,0,2) eq '  ') {
                            $line=trim($line);
                            $seqProt="$seqProt"."<BR>$line";
                        }
                    }
                }
                
                ###> Recovering results send by Uniprot
                if ($lineCount>=0) {
                    $lineCount=0;
                    foreach my $line (split(/\n/,$responseSearchProt->content)) {
                        $lineCount++;
                        next if $lineCount==1;
                        my ($uniAC,$uniID,$pNames,$length,$geneStrg)=split(/\t/,$line);
                        $uniAC=trim($uniAC);
                        my $searchMaj=uc($search);
                        my $searchMin=lc($search);
                        
                        if ("@tabUniAC" =~ /$search/ || "@tabUniAC" =~ /$searchMaj/ || "@tabUniAC" =~ /$searchMin/){
                            if ($uniAC eq $search || $uniAC eq $searchMaj || $uniAC eq $searchMin) {&listMatchProt($search,$uniAC,$uniID,$pNames,$length,$geneStrg,$lineCount,\%uniprotResultList,\%accessionSeqList);}
                        }
                        else{&listMatchProt($search,$uniAC,$uniID,$pNames,$length,$geneStrg,$lineCount,\%uniprotResultList,\%accessionSeqList);}
                    }
                }
            }
             ###> Display results 
            if (%uniprotResultList){
                print qq |
                <SCRIPT LANGAGE="JavaScript">document.getElementById('waitDIV').style.display='none'</SCRIPT>
                <DIV><BR><BR><CENTER><FONT class="title2">Results for "$search"</FONT><BR><BR>
                <FORM method="POST" action="./searchSwathLibrary.cgi" name="selectID" enctype="multipart/form-data">
                <TABLE width=45% cellspacing=0>
                <TR  bgcolor="$darkColor">
                <TH class="rbBorder" height=30px>&nbsp;Protein ID (AC)&nbsp;</TH>
                <TH class="rbBorder" height=30px>&nbsp;Gene Names&nbsp;</TH>
                <TH class="rbBorder" height=30px>&nbsp;Protein Names&nbsp;</TH>
                <TH class="rbBorder" height=30px>&nbsp;AA&nbsp;</TH>
                <TH class="bBorder"  height=30px>&nbsp;&nbsp;#&nbsp;Peptides&nbsp;&nbsp;</TH>
                </TR>
                |;
                
                foreach my $uniID (sort keys %uniprotResultList){
                    my ($uniAC,$genes,$protName,$length,$peptideNumber,$sequenceProt,$acList)=split('&', $uniprotResultList{$uniID});
                    print qq
                    |<TR class=\"$trClass\" >
                    <TD align="center" class="padtab" nowrap valign=center><A href="http://www.uniprot.org/uniprot/$uniAC" target="_blank"/>|; if ($peptideNumber==0){print "<S>$uniID ($uniAC)</S>";} else {print "$uniID ($uniAC)";} print qq | </A></TD>
                    <TD align="center" class="padtab" valign=top>|; if ($genes eq "") {print "N/A";} else {print $genes;} print qq |</TD>
                    <TD align="left" class="padtab" valign=top>$protName</TD>
                    <TD align="center" class="padtab" valign=center>$length</TD>
                    <TD align="center" class="padtab" valign=center>|; if ($peptideNumber==0){print "<B>-</B>";} else{print "<A href=\"javascript:ajaxSearchLibrary('$uniAC','$acList','$sequenceProt')\" /><B>$peptideNumber</B></A>"}print qq|</TD>
                    </TR>
                    |;
                    $trClass=($trClass eq 'row1')? 'row2' : 'row1';
                }
                print "</TABLE></FORM></CENTER></DIV>";
            }
            else {print "<SCRIPT LANGAGE=\"JavaScript\">document.getElementById('waitDIV').style.display='none'</SCRIPT><BR><BR><FONT class=\"title2\">No match for \"$search\". </FONT>";}
        }
        else{
            ## get peptide number associated to this protein
            my $peptideNumber=`grep -c $search $searchFile`;
            next unless $peptideNumber;
            chomp $peptideNumber;
            if ($peptideNumber) {
                ## recvering protein sequence in fasta file
                my $sequenceProt;
                foreach my $dbID (keys %databank){
                    open(DB,"<$promsPath{data}/banks/db_$dbID/$databank{$dbID}") or die ("open : $!");
                    my $match=0;
                    while (<DB>) {
                        if ($_=~/>/) {
                            if ($_=~/$search/) {$match=1;}
                            else{$match=0;}
                        }
                        else{
                            chomp $_ if $match;
                            $sequenceProt.=$_ if $match;
                        }
                    }
                    last if $sequenceProt;
                    close DB;
                }
                $nbSearch++;
                &searchLibrary($search,$search,$sequenceProt,$nbSearch);
            }
            else{print "<SCRIPT LANGAGE=\"JavaScript\">document.getElementById('waitDIV').style.display='none'</SCRIPT><BR><BR><FONT class=\"title2\">No match for \"$search\". </FONT>";}
            
        }         
    } 
}
print qq
|
<DIV id="waitDIV2" style="display:none"><BR><BR><FONT class="title3">Fetching data ...</FONT><BR><BR><IMG src="$promsPath{images}/scrollbarGreen.gif"><BR><BR></DIV>
<DIV id="result"></DIV>

<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">setPopup()</SCRIPT>
</CENTER>
</BODY>
</HTML>


|;



sub listMatchProt{
    my ($search,$uniAC,$uniID,$pNames,$length,$geneStrg,$lineCount,$refUniprotResultList,$refAccessionSeqList)=@_;
    next if $lineCount==1;
    
    my @protNames=split("[\[]",$pNames);
    my $protName=$protNames[0];
    my @geneName=split(' ',join('',split('[;] ',$geneStrg)));
    my $genes;
    if (scalar @geneName == 0) {
        $genes="";
    }
    else{
        $geneName[0]="<B>$geneName[0]</B>";
        $genes=join(', ',@geneName);
    }
    
    my $peptideNumber= ($spcLibFormat eq 'TPP') ? `grep -c $uniAC $searchFile` : `grep $uniAC $searchFile | cut -f4 | sort | uniq | wc -l`;
    chomp $peptideNumber;
   
    ####>Recovering AC numbers and sequence of the protein
    my $sequenceProt='';
    my $acList='';
    my @value;
    foreach my $value (values(@{$refAccessionSeqList->{$uniAC}})){
        $value=trim($value);
        push(@value,$value);
    }
    $sequenceProt=pop @value;
    $acList=join('@',@value);

    $refUniprotResultList->{$uniID}="$uniAC&$genes&$protName&$length&$peptideNumber&$sequenceProt&$acList";
}

###############################################
####>Fonction to search peptide in library<####
###############################################

sub searchLibrary{
    my ($selectProt,$accList,$seqProt,$nbSearch)=@_;
    $nbSearch=1 unless $nbSearch;
    
    ###> Linear sequence of the protein
    my $seqPeptaille=join("",split(/\s/,$seqProt));
    my @seqPeptaille=split(/<BR>/,$seqPeptaille);
    my $sequenceProtein=join('',@seqPeptaille); 
    
    my @accessList=split("@",$accList);
    my @searchIDProt=split('[\(]',$selectProt);
    my $searchProt=$searchIDProt[0];
    
    ###> Display peptide informations
    open (INFILE,"<","$searchFile") or die("open: $!");
    
    
    my $i=1;
    my $line;
    my ($peptide,$charge,$sequencePeptide,$masseTot,$masseParent,$protID,$specificity,$irt,$modification);
    my (@begPeptide,@endPeptide,%proteinList);
    my (%peptideInfo,%posBegPep);
    my $seqPB=0;
    if($spcLibFormat eq 'TPP') {
        while ($line=<INFILE>){
            if ($line=~/^Name:(.*)/){
                ($peptide,$charge,$sequencePeptide,$masseTot,$masseParent,$protID,$specificity,$irt,$modification)=('','','','','','','','','');
                my @pepInfo=split("[/\]",$1);
                ($peptide=$pepInfo[0])=~s/\s//g;
                #$peptide=~s/\s//g;
                #my @peptideModif=split(/n*\[\w*\]/,$peptide);
                #if ($peptide=~m/n*\[\w*\]/) {$sequencePeptide=join("",@peptideModif); }
                
                if ($peptide=~m/n*\[\d*\]/) {($sequencePeptide=$peptide)=~s/n|\[\d*\]//g;}
                else{$sequencePeptide=$peptide;}
                $charge=$pepInfo[1];
            }
            elsif (substr($line,0,3) eq "MW:") {$masseTot=substr($line,4,-1);}
            elsif (substr($line,0,12) eq "PrecursorMZ:") {$masseParent=substr($line,13,-1); }
            elsif (substr($line,0,8) eq "Comment:") {
                my @comment=split(' ',$line);
                foreach my $com (@comment){
                    if(substr($com,0,4) eq "iRT=") {
                        my @irtLine=split(",",substr($com,4));
                        $irt=$irtLine[0];
                        $irt=$irt*1;
                    }
                    elsif(substr($com,0,5) eq "Mods=") {
                        my %modList;
                        if (substr($com,0,6) eq "Mods=0") {$modification="-";}
                        else{
                            my $var=substr($com,5);
                            my @result=&promsMod::convertVarModStringSwath($var);
                            if (scalar @result ==1) {
                                my @resTab=split(/!/,$result[0]);
                                if ($resTab[2] eq "=") {$modification=$resTab[0]." ("."N-term".")";}
                                else{$modification=" @ ".$resTab[0]." (".$resTab[1].":".$resTab[2].")";}
                            }
                            else{
                                foreach my $res (@result){
                                    my @resTab=split(/!/,$res);
                                    if ($resTab[2] eq "=") {
                                        $modList{$resTab[0]}{"N-term"}{""}=1;
                                    }
                                    else{$modList{$resTab[0]}{"$resTab[1]:"}{$resTab[2]}=1;}
                                }
                                foreach my $varMod (keys(%modList)){
                                    foreach my $res (keys %{$modList{$varMod}}){
                                        if ($modification eq "") {
                                            $modification=" @ ".$varMod." ($res".join(".",keys %{$modList{$varMod}{$res}}).")";
                                        }
                                        else{
                                            $modification=$modification." @ ".$varMod." ($res".join(".",keys %{$modList{$varMod}{$res}}).")";
                                        }
                                    }
                                }
                            }
                        }    
                    }
                    elsif (substr($com,0,8) eq "Protein=") {
                        my (@protIDs,@protList);
                        if ($split==0) {
                            @protIDs=split(/\//,substr($com,8));
                            if ("@protIDs"=~/$searchProt/){$protID=$searchProt;}
                            else{
                                foreach my $ac (@accessList){
                                    if ("@protIDs"=~/$ac/) {$protID=$ac;}
                                }
                            }
                            $specificity=100/$protIDs[0];
                            shift(@protIDs);
                            foreach my $prot (@protIDs){
                                if ($prot=~/(tr|sp)\|(\w*)/) {
                                    push @protList,$2;
                                }
                                else{
                                    push @protList,$prot;
                                }
                            }
                            $specificity = sprintf("%0.1f", $specificity);
                            $specificity=$specificity*1;
                        }
                        elsif ($split==1){
                            @protIDs=split(/\//,substr($com,10));
                            my @numberIDs=split(/_/,$protIDs[0]);
                            if ("@protIDs"=~/$searchProt/) {$protID=$searchProt;}
                            else{
                                foreach my $ac (@protIDs){
                                    if ("@accessList"=~/$ac/ ){$protID=$ac;} 
                                }
                            }
                            $specificity=100/$numberIDs[2];
                            shift(@protIDs);
                            foreach my $prot (@protIDs){
                                if ($prot=~/(tr|sp)\|(\w*)\|/) {
                                    push @protList,$2;
                                }
                                else{push @protList,$prot;}
                            }
                            $specificity = sprintf("%0.1f", $specificity);
                            $specificity=$specificity*1;
                        }
                        if ($protID && "@accessList" =~ /$protID/) {
                            $proteinList{"$peptide\_$charge"}=\@protList;
                        }
                    }
                }
                if ($protID && "@accessList" =~ /$protID/) {
                    my $inversionLI=0;
                    my ($firstPosition,$lastPosition,$positionPeptide);
                    ##> match peptide sequence on protein sequence and recover peptide positions
                    if ($sequenceProtein=~m/$sequencePeptide/) {
                        while ($sequenceProtein=~m/$sequencePeptide/g){
                            $firstPosition=length($`)+1;
                            $lastPosition=$firstPosition+length($sequencePeptide)-1;
                            push @begPeptide,$firstPosition;
                            push @endPeptide,$lastPosition;
                            $positionPeptide=$positionPeptide."$firstPosition-$lastPosition";
                        }
                    }
                    else{ #for Leucine and Isoleucine inversion
                        my $sequencePeptideMod=$sequencePeptide;
                        $sequencePeptideMod=~s/[IL]/[IL]/g;
                        if ($sequenceProtein=~m/$sequencePeptideMod/ ) {
                            $firstPosition=length($`)+1;
                            $lastPosition=$firstPosition+length($sequencePeptide)-1;
                            $positionPeptide=$positionPeptide."$firstPosition-$lastPosition";
                            $inversionLI=1;
                            push @begPeptide,$firstPosition;
                            push @endPeptide,$lastPosition;
                        }
                    }
                    $posBegPep{$peptide}=$firstPosition;
                    push @{$peptideInfo{$peptide}{$charge}{"$modification\_$irt"}},$inversionLI;
                    push @{$peptideInfo{$peptide}{$charge}{"$modification\_$irt"}},$positionPeptide;
                    push @{$peptideInfo{$peptide}{$charge}{"$modification\_$irt"}},$masseParent;
                    push @{$peptideInfo{$peptide}{$charge}{"$modification\_$irt"}},$masseTot;
                    push @{$peptideInfo{$peptide}{$charge}{"$modification\_$irt"}},$specificity;
                    
                }
            }
        }
    } elsif($spcLibFormat eq 'SPC') {
        # Parse header to retrieve specific columns
        my $line=<INFILE>;
        my @headers = split(/[,;\t\s]/, $line);
        my (%colName2Index, %massMods);
        foreach my $i (0 .. $#headers) {
            $colName2Index{$headers[$i]} = $i;
        }
        
        my $accessListStrPattern = join('|', @accessList);
        while (($line=<INFILE>)) {
            my @lineContent = split(/[\t]/, $line);
            next unless($lineContent[$colName2Index{'ProteinGroups'}] =~ /($accessListStrPattern)/);
            
            my $protID = $1;
            my ($peptide,$charge,$sequencePeptide,$peptideInt,$massTot,$massParent,$specificity,$irt,$modification)=('','','','','','','','','');
            my ($firstPosition,$lastPosition,$positionPeptide);
            
            $sequencePeptide = $lineContent[$colName2Index{'StrippedPeptide'}];
            $peptideInt = $lineContent[$colName2Index{'IntModifiedPeptide'}];
            $peptideInt = substr($peptideInt, 1, -1) if($peptideInt);
            $peptide = $lineContent[$colName2Index{'ModifiedPeptide'}];
            $peptide = substr($peptide, 1, -1) if($peptide);
            $charge = $lineContent[$colName2Index{'PrecursorCharge'}];
            
            $irt = $lineContent[$colName2Index{'iRT'}];
            if($irt) {
                $irt =~ s/,/\./;
                $irt = sprintf("%0.1f", $irt);
            }
            
            # Modifications
            my @result=&promsMod::convertVarModStringSpectronaut($peptide);
            my (%modList, %varMods);
            if (scalar @result ==1) {
                my @resTab=split(/!/,$result[0]);

                if(!$massMods{$resTab[0]}) {
                    my $modID=&promsMod::getModificationIDfromString($dbh,$resTab[0],$resTab[1],undef);       #fetching modification ID
                    ($massMods{$resTab[0]})=$dbh->selectrow_array("SELECT MONO_MASS FROM MODIFICATION WHERE ID_MODIFICATION=$modID");
                }

                if ($resTab[1] eq "=") {
                    $modification=$resTab[0]." ("."N-term".")";
                    $varMods{"="} = $massMods{$resTab[0]};
                } else{
                    $modification = " @ ".$resTab[0]." (".$resTab[1].":".$resTab[2].")";
                    $varMods{$resTab[2]} = $massMods{$resTab[0]};
                }
            } else {
                foreach my $res (@result){
                    my @resTab=split(/!/,$res);
                    
                    if(!$massMods{$resTab[0]}) {
                        my $modID=&promsMod::getModificationIDfromString($dbh,$resTab[0],$resTab[1],undef);       #fetching modification ID
                        ($massMods{$resTab[0]})=$dbh->selectrow_array("SELECT MONO_MASS FROM MODIFICATION WHERE ID_MODIFICATION=$modID");
                    }
                    
                    if ($resTab[1] eq "=") {
                        $modList{$resTab[0]}{"N-term"}{""}=1;
                        $varMods{"="} = $massMods{$resTab[0]};
                    }
                    else{
                        $modList{$resTab[0]}{"$resTab[1]:"}{$resTab[2]}=1;
                        $varMods{$resTab[2]} = $massMods{$resTab[0]};
                    }
                }
                foreach my $varMod (keys(%modList)){
                    foreach my $res (keys %{$modList{$varMod}}){
                        if ($modification eq "") {
                            $modification=" @ ".$varMod." ($res".join(".",keys %{$modList{$varMod}{$res}}).")";
                        }
                        else{
                            $modification=$modification." @ ".$varMod." ($res".join(".",keys %{$modList{$varMod}{$res}}).")";
                        }
                    }
                }
            }
            
            next if($peptideInfo{$peptideInt}{$charge}{"$modification\_$irt"});
            
            $massParent = $lineContent[$colName2Index{'PrecursorMz'}];
            $massParent =~ s/,/\./;
            $massParent = sprintf("%0.4f", $massParent) if($massParent);
            
            $masseTot = &promsMod::mrCalc($sequencePeptide, \%varMods);
            $masseTot = sprintf("%0.4f", $masseTot) if($masseTot);
            
            my @protList = split(/;/, $lineContent[$colName2Index{'ProteinGroups'}]);
            $specificity = (1/scalar @protList)*100;
            $specificity = sprintf("%0.1f", $specificity) if($specificity);
            my $inversionLI = 0;
            
            
            ##> match peptide sequence on protein sequence and recover peptide positions
            if ($sequenceProtein=~m/$sequencePeptide/) {
                while ($sequenceProtein=~m/$sequencePeptide/g){
                    $firstPosition=length($`)+1;
                    $lastPosition=$firstPosition+length($sequencePeptide)-1;
                    push @begPeptide,$firstPosition;
                    push @endPeptide,$lastPosition;
                    $positionPeptide=$positionPeptide."$firstPosition-$lastPosition";
                }
            }
            else{ #for Leucine and Isoleucine inversion
                my $sequencePeptideMod=$sequencePeptide;
                $sequencePeptideMod=~s/[IL]/[IL]/g;
                if ($sequenceProtein=~m/$sequencePeptideMod/ ) {
                    $firstPosition=length($`)+1;
                    $lastPosition=$firstPosition+length($sequencePeptide)-1;
                    $positionPeptide=$positionPeptide."$firstPosition-$lastPosition";
                    $inversionLI=1;
                    push @begPeptide,$firstPosition;
                    push @endPeptide,$lastPosition;
                }
            }
            
            push @{$peptideInfo{$peptideInt}{$charge}{"$modification\_$irt"}},$inversionLI;
            push @{$peptideInfo{$peptideInt}{$charge}{"$modification\_$irt"}},$positionPeptide;
            push @{$peptideInfo{$peptideInt}{$charge}{"$modification\_$irt"}},$massParent;
            push @{$peptideInfo{$peptideInt}{$charge}{"$modification\_$irt"}},$masseTot;
            push @{$peptideInfo{$peptideInt}{$charge}{"$modification\_$irt"}},$specificity;
            $proteinList{"$peptideInt\_$charge"}=\@protList;
            $posBegPep{$peptideInt}=$firstPosition;
        }
    }
            
    close INFILE;
    my $trClass='row1';
    
    my $buttonName=$nbSearch."BUTTON";
    my $divName=$nbSearch."DIV";
    
    print "<CENTER><BR><FONT class=\"title2\">Results for \"$searchProt\"</FONT><BR>" unless $identLibType=~/UNIPROT|UniProt/;
    print qq
    |
    <TABLE width="70%" bgcolor="$darkColor">
    <TR>
        <BR><BR><TD><INPUT type="button" class="font11" id="peptide$buttonName" value="Hide" onclick="showHideFunction('peptide$nbSearch')" style="width:55px"/></TD>
        <TH bgcolor="$lightColor" align="left" width="100%" nowrap>&nbsp;Peptide list for $searchProt&nbsp;</TH>
    </TR>
    </TABLE>
    <DIV id="peptide$divName" ><BR>
    <SCRIPT LANGAGE="JavaScript">document.getElementById('waitDIV').style.display='none'</SCRIPT>
    <TABLE class="local" cellspacing=0>
    <TR bgcolor="$darkColor">
    <TH class="rbBorder" height=30px width=30px>&nbsp;#&nbsp;</TH>
    <TH class="rbBorder" height=30px>&nbsp;Sequence&nbsp;</TH>
    <TH class="rbBorder" height=30px>&nbsp;Modifications&nbsp;</TH>
    <TH class="rbBorder" height=30px>&nbsp;Position&nbsp;</TH>
    <TH class="rbBorder" height=30px>&nbsp;M/Z&nbsp;</TH>
    <TH class="rbBorder" height=30px>&nbsp;Charge&nbsp;</TH>
    <TH class="rbBorder" height=30px>&nbsp;iRT time&nbsp;</TH>
    <TH class="rbBorder" height=30px>&nbsp;Specificity (%)&nbsp;</TH>
    <TH class="bBorder"  height=30px>&nbsp;Found with&nbsp;</TH>
    </TR>
    |;
    my $nbPeptide=1;
    foreach my $peptide (sort {$posBegPep{$a} <=> $posBegPep{$b} || lc($a) cmp lc($b) || $a cmp $b} keys %peptideInfo){
        (my $pepseq=$peptide)=~s/n|\[[+-]?\d*\]|//g; # OpenSwath + Spectronaut Format
        foreach my $charge (sort keys %{$peptideInfo{$peptide}}){
            foreach my $peptideInfo (keys %{$peptideInfo{$peptide}{$charge}}){
                my ($modification,$irt)=split(/_/,$peptideInfo);
                my $inversionLI=$peptideInfo{$peptide}{$charge}{$peptideInfo}[0];
                my $positionPeptide=$peptideInfo{$peptide}{$charge}{$peptideInfo}[1];
                my $masseParent=$peptideInfo{$peptide}{$charge}{$peptideInfo}[2];
                my $masseTot=$peptideInfo{$peptide}{$charge}{$peptideInfo}[3];
                my $specificity=$peptideInfo{$peptide}{$charge}{$peptideInfo}[4];
                my $refProtList=$proteinList{"$peptide\_$charge"};
                
                my $mod=$modification unless $modification eq '-';
                $mod=($mod)? $mod : '';
                $modification=~s/^\s@//;
                $modification=~s/@/<BR>/g;
                print qq
                |<TR class=\"$trClass\" >
                <TD align="rigth" class="padtab" valign="center">&nbsp;$nbPeptide&nbsp;</TD>
                <TD align="left" class="padtab" valign="center">&nbsp;<span id="$nbPeptide" class="drawViewPept" onclick="pepView('$nbPeptide','$peptide\_$charge\_$mod','$irt','$masseParent','$masseTot','$libraryID')" onmouseover="popup('Click to display <B>Fragmentation Spectrum</B>.')" onmouseout="popout()">|; if ($inversionLI==1) {print "<FONT class=\"overlapRes2\"><B>$pepseq</B></FONT>*";$seqPB=1;}else{print "<B>$pepseq</B>";} print qq |</span>&nbsp;</TD>
                <TD align="center" class="padtab" valign="center">&nbsp;$modification&nbsp;</TD>
                <TD align="center" class="padtab" valign="center">&nbsp;|; if ($positionPeptide eq "") {print "No match on the protein."}else{print $positionPeptide} print qq | &nbsp;</TD>
                <TD align="center" class="padtab" valign="center">&nbsp;$masseParent&nbsp;</TD>
                <TD align="center" class="padtab" valign="center">&nbsp;$charge+&nbsp;</TD>
                <TD align="center" class="padtab" valign="center">&nbsp;$irt&nbsp;</TD>
                <TD align="center" class="padtab" valign="center">&nbsp;$specificity&nbsp;</TD>
                <TD align=\"center\" class=\"padtab\" valign=\"center\">&nbsp;|;
                if (scalar @{$refProtList}>5){
                    foreach my $protName (@{$refProtList}[0..4]){
                        if ($identLibType=~/UNIPROT|UniProt/){print "<A href=\"http://www.uniprot.org/uniprot/$protName\" target=\"_blank\">$protName</A>" ;}
                        else{print "<A href=\"http://www.uniprot.org/uniprot/?query=$protName\" target=\"_blank\">$protName</A>";}
                        print ", " unless $protName eq @{$refProtList}[4];
                    }
                    print "<A id=\"proteinList\" href=\"javascript:protListView($nbPeptide)\">...</A>";
                    print "<DIV id=\"DIVprotList$nbPeptide\" style=\"display:none\">";
                    for ($i=5;$i<=$#{$refProtList};$i++){
                        if ($identLibType=~/UNIPROT|UniProt/){print "<A href=\"http://www.uniprot.org/uniprot/@{$refProtList}[$i]\" target=\"_blank\">@{$refProtList}[$i]</A>";}
                        else{print "<A href=\"http://www.uniprot.org/uniprot/?query=@{$refProtList}[$i]\" target=\"_blank\">@{$refProtList}[$i]</A>";}
                        if ($i!=5 && $i%5==0 && $i!=$#{$refProtList}) {
                            print"<BR>";
                        }
                        elsif($i!=$#{$refProtList}){print ", ";}
                    }
                    print "</DIV>";
                }
                else{
                    foreach my $protName (@{$refProtList}){
                        if ($identLibType=~/UNIPROT|UniProt/){print "<A href=\"http://www.uniprot.org/uniprot/$protName\" target=\"_blank\">$protName</A>";}
                        else{print "<A href=\"http://www.uniprot.org/uniprot/?query=$protName\" target=\"_blank\">$protName</A>";}
                        print ", " unless $protName eq @{$refProtList}[-1];
                    } 
                }
                print "&nbsp;</TD></TR>";
                $nbPeptide++;
                $trClass=($trClass eq 'row1')? 'row2' : 'row1';
            }
        }
    }
    
    print "</TABLE>";
    if ($seqPB==1) {
        print "* Inversion of a Leucine and an Isoleucine in peptide sequence in comparison with the proteine sequence.";
    }
    print "</DIV><BR><BR>";


    print qq
    |
    <TABLE width="70%" bgcolor="$darkColor">
    <TR>
        <TD><INPUT type="button" class="font11" id="sequence$buttonName" value="Show" onclick="showHideFunction('sequence$nbSearch')" style="width:55px"/></TD>
        <TH bgcolor="$lightColor" align="left" width="100%" nowrap>&nbsp;Detailed sequence coverage for $searchProt&nbsp;</TH>
    </TR>
    </TABLE><BR><BR><BR>
    <DIV id="sequence$divName" style="display:none">
    |;
    
    my (%positionPeptide,$pepCoverage,@positionPeptideList);
    ###> Status of peptide positions
    $positionPeptide{0}=0; # first residue is set to 0
	foreach my $i (0..$#begPeptide) {
		$positionPeptide{$begPeptide[$i]-1}++; # -1 -> index not pos!
		$positionPeptide{$endPeptide[$i]-1}--;
    }
    
    ###> Coverage calcul
	my $prevPosition=0;
	my $matchBeg=0;
	foreach my $position (sort {$a<=>$b} keys %positionPeptide) {
        push @positionPeptideList,$position;
		next if $position==0;
		$positionPeptide{$position}+=$positionPeptide{$prevPosition}; 
		if ($matchBeg==0 && $positionPeptide{$position}>=1) { # match begins
			$matchBeg=$position;
		}
		elsif ($matchBeg>0 && $positionPeptide{$position}==0) { # match ends
			$pepCoverage+=($position-$matchBeg+1);
			$matchBeg=0;
		}
		$prevPosition=$position;
	}
    my $peptideCoverage=sprintf ("%.1f",$pepCoverage/length($sequenceProtein)*100);
    print "<B>Peptide coverage : ",$peptideCoverage,"%</B><BR>";
    
    ###>display protein sequence with peptides match
    my $aaClass='overlapRes0';
    my $overlap=0;
    my @protResidues=split(//,$sequenceProtein);
    my $prec=-1;
    my $lastPosIdx;
    print "<TABLE>";
    foreach my $posIdx (0..$#protResidues) {
        my $resPos=$posIdx+1;
        my $newline=0;
        if ($posIdx/60==int($posIdx/60)) {# new line
            $newline=1;
            if ($posIdx == 0) {print "<TR><TD align=\"right\" valign=\"top\" >$resPos&nbsp;</TD><TD class=\"seq\"><FONT class=\"$aaClass\">";}
            else{print "</FONT><BR><BR></TD></TR><TR><TD align=\"right\" valign=\"top\">$resPos&nbsp;</TD><TD class=\"seq\" ><FONT class=\"$aaClass\">";}
            $lastPosIdx=$posIdx+59;
            $lastPosIdx=$#protResidues if $lastPosIdx > $#protResidues;
        }
        elsif ($posIdx>0 && $posIdx/10==int($posIdx/10)) { # space
            print "&nbsp;";
        }
        ##>Match residue
        if (defined $positionPeptide{$posIdx}){
            if ($posIdx==0) {
                if ($positionPeptide{$posIdx}!=0) {$aaClass="overlapRes1"; print "</FONT><FONT class=\"$aaClass\">$protResidues[$posIdx]";}
                else{print $protResidues[$posIdx];}
                
            }
            elsif ($positionPeptide{$posIdx}>=1 ) {
                if ($positionPeptide{$positionPeptideList[$prec]}==0) {
                    $aaClass="overlapRes1";
                    print "</FONT><FONT class=\"$aaClass\">$protResidues[$posIdx]";
                }
                elsif($positionPeptide{$positionPeptideList[$prec]}>=1){
                    $aaClass="overlapRes1";
                    print "$protResidues[$posIdx]</FONT>";
                    print "<FONT class=\"$aaClass\">" unless $posIdx==$lastPosIdx;
                }
            }
            elsif($positionPeptide{$posIdx}==0){
                $aaClass="overlapRes0";
                print "$protResidues[$posIdx]</FONT>";
                print "<FONT class=\"$aaClass\">" unless $lastPosIdx == $#protResidues;
            }
            $prec ++;
        }
        ##>Normal residue
        else{
            print $protResidues[$posIdx];
        }
        ##>End
        if ($posIdx==$#protResidues) {
            print "<BR><BR></FONT></TD></TR>";
        }
        
    }
    print "<BR><BR></FONT></TD></TR></TABLE>";
    
    print "</DIV>";

    $dbh->disconnect;
}

####>Revision history<####
# 1.1.0 [ENHANCEMENT] Handles Spectronaut library format (VS 02/10/20)
# 1.0.10 Minor modif (MLP 19/12/17)
# 1.0.9 Add modification to allow the search of proteins for libraries with no uniprot ID (MLP 22/03/17)
# 1.0.8 Minor modif (MLP 31/01/17)
# 1.0.7 Add loading bar and Uniprot link on protein ACC (MLP 22/12/16)
# 1.0.6 Add &promsMod::cleanParameters verification (MLP 28/09/16)
# 1.0.5 Update 'searchLibrary' to sort peptides by position (MLP 22/07/16)
# 1.0.4 Minor modification of query names (MLP 22/06/16)
# 1.0.3 Add checkbox to focus uniprot research on a specie (MLP 18/05/16)
# 1.0.2 Minor modifications to list associated proteins for each peptide (MLP 02/05/16)
# 1.0.1 Change Swath_Lib path ($promsPath{data}/Swath_Lib to $promsPath{swath_lib}) (MLP 13/04/16)
# 1.0.0 Allow the search of peptide in a swath library (MLP 11/04/16)

