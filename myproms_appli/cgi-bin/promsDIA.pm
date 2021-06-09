#!/usr/local/bin/perl -w

################################################################################
# promsDIA.pm     1.1.4                                                      #
# Authors: M. Le Picard (Institut Curie)                                       #
# Contact: myproms@curie.fr                                                    #
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
package promsDIA;

require Exporter;

@ISA=qw(Exporter);
@EXPORT=qw();
@EXPORT_OK=qw();
$VERSION=1.00;

use strict;
use DBD::mysql;	# MySQL
use promsConfig;
use Data::Dumper;
use POSIX qw(strftime); # to get the time


sub exportLibrarySCRIPT {
    my ($dbh, $libraryID, $libraryName, $workDir, $refPromsPath, $startTime, $loadingDivID, $loadingSPAN, $w,$k,$p,$g,$f,$t,$y,$m,$labelling,$other,$lMax,$lMin,$s,$x,$n,$o,$pepMod,$processText,$protList)=@_;
    my $paramFile=$libraryName."_param";
    $processText=($processText)? $processText : '';
    
    ###################################
    ####<Processing submitted form>####
    ###################################
    my ($divID,$now,$waitTime,$status,$loading);
    if ($processText){
        $divID="document.getElementById('waitDIV')";
        $now=strftime("%s",localtime); # in sec
        $waitTime=strftime("%Hh %Mm %Ss",localtime($now-$startTime-3600));
        $status="Updated $waitTime ago";
        print "<BR><SCRIPT LANGUAGE=\"JavaScript\">$divID.innerHTML=\"\";$divID.innerHTML='$status';</SCRIPT>";
    }
    
    
    open(PARAMFILE,">","$workDir/$paramFile") or die ("open: $!");
    print PARAMFILE "Export for $k \nION_MASS_LIMITS:$lMin,$lMax\nION_TYPE:$s\nION_CHARGE:$x\nION_NUMBER_PER_PEPTIDE:$o,$n\n";
    print PARAMFILE "MAXIMUM_ERROR_ALLOWED:$p\nTIME_SCALE:$t\nUIS_ORDER:$y\nSWATH_FILE:$w\n";
    my $param='';
    if ($f){
        $param="-f $f ";
        print PARAMFILE "FASTA_FILE:$f\n";
    }
    if ($g){
        $param.="-g $g ";
        print PARAMFILE "FRAGMENT_MASS_MODIFICATIONS_ALLOWED:$g\n";
    }
    if ($m){
        $param.="-m $workDir/$m ";
        print PARAMFILE "MODIFICATION_FILE:$m\n";
    }
    if ($labelling){
        $param.="-i $labelling ";
        print PARAMFILE "LABELLING_FILE:$labelling\n";
    }
    
    my $libraryVersion=$dbh->selectrow_array("SELECT VERSION_NAME FROM SWATH_LIB WHERE ID_SWATH_LIB=$libraryID");
    print PARAMFILE "LIBRARY_VERSION:$libraryVersion\n";
    
    
    #################################################
    ###> compute excluded peptides modifications <###
    #################################################
    my %massAAave=&promsConfig::getMassAAave; 
    my @excludeMod;
    my $sthModInfo=$dbh->prepare("SELECT SLM.SPECIFICITY, MONO_MASS FROM MODIFICATION M, SWATH_LIB_MODIFICATION SLM WHERE SLM.ID_MODIFICATION=M.ID_MODIFICATION AND M.ID_MODIFICATION=? AND ID_SWATH_LIB=$libraryID");
    if($pepMod){
        foreach my $selectMod ($pepMod){
            $sthModInfo->execute($selectMod);
            my ($specificity,$monoMass)=$sthModInfo->fetchrow_array;
            foreach my $aa (split(//,$specificity)){
                if ($aa eq '='){
                    my $massExcluded=$monoMass+1.00794;
                    push @excludeMod, sprintf("%0.f",$massExcluded);
                }else{
                    my $massExcluded=$monoMass+$massAAave{$aa};
                    push @excludeMod, sprintf("%0.f",$massExcluded);
                }
            }   
        }
    }
    
    #####################################
    ###> recovering the protein list <###
    #####################################
    my @protList;
    if ($protList) {
        open(PROTLIST,"<$workDir/proteinList.txt") or die ("open: $!");
        while (<PROTLIST>) {
            foreach my $protein (split('[;,\s]',$_)){
                push @protList, $protein;
            }
        }
        close PROTLIST;
    }
    
    if ($processText){
       ###> Process DIV
        $loading="<progress value=\"0\" max=\"100\" style=\"width:400px\"></progress>";
        print "<SCRIPT LANGAGE=\"JavaScript\">$loadingDivID.innerHTML=\"\";$loadingDivID.innerHTML='<B>Step 1/2 :</B> Export library for PeakView ... <BR><BR> $loading';$loadingSPAN.innerHTML='0%';</SCRIPT>" if $k eq 'peakview';
        print "<SCRIPT LANGAGE=\"JavaScript\">$loadingDivID.innerHTML=\"\";$loadingDivID.innerHTML='$processText<BR><BR> $loading';$loadingSPAN.innerHTML='0%';</SCRIPT>" if $k eq 'openswath';
    }
    
    #######################################
    ###> Delete all N-ter acetylations <###
    #######################################
    my (%mzHash,@validProt);
    open(LIBFILE,"<","$workDir/$libraryName.sptxt") or die ("open: $!");
    open(OUTLIBFILE,">","$workDir/convert$k.sptxt") or die ("open: $!");
    my $saveLine;
    my ($match,$matchList)=(0,0);
    while (my $line=<LIBFILE>) {
        if ($line=~/^Name:(.*)/) {
            $match=0;
            if ($matchList) {
                print OUTLIBFILE $saveLine;
            }
            $matchList=0;
            $saveLine='';
            my $modOK=1;
            if(@excludeMod){
                foreach my $modif (@excludeMod){
                    if($1=~/\[$modif\]/){
                        $modOK=0;
                        last;
                    }
                }
                #if ($1=~/n\[43\]/) {
                #    $modOK=0;
                #}
            }
            if ($modOK) {
                $match=1;
                print OUTLIBFILE $line unless @protList;        ## export all the library
                $saveLine.=$line;
            }
        }
        else{
            if ($match==1) {
                $saveLine.=$line;
                print OUTLIBFILE $line unless @protList;        ## export all the library
                if ($line=~/^PrecursorMZ: (.+)/) {
                    $mzHash{"$1"}=1;
                }
                elsif ($line=~/^Comment: (.+)/ && @protList) {      ## export just the selected proteins
                    foreach my $protein (@protList){
                        if ($1=~/$protein/) {
                            push @validProt,$protein unless "@validProt"=~/$protein/;
                            $matchList=1;
                            last;
                        }
                    }
                }
            }
        }     
    }
    close LIBFILE;
    close OUTLIBFILE;
    
    ##############################
    ###> Run Spectrast2tsv.py <###
    ##############################
    my $command;
    my $finalFile;
    my $output="$workDir/sortieExport.txt";
    if ($k eq "peakview"){
        $command=$refPromsPath->{"msproteomicstools"}."/spectrast2tsv.py -l $lMin,$lMax -s $s -x $x -o $o -n $n -p $p $other $param -t $t -y $y -w $workDir/$w -k $k -a $workDir/$libraryName\"_peakview.tsv\" $workDir/convertpeakview.sptxt >$output 2>&1;";
        $command.="sed s/TRUE/FALSE/g $workDir/$libraryName\_peakview.tsv >$workDir/peakviewfile.txt;";
        $finalFile=$libraryName."_peakview.txt";
    }
    elsif ($k eq "openswath"){
        $finalFile=$libraryName."_openswath.tsv";
        $command=$refPromsPath->{"msproteomicstools"}."/spectrast2tsv.py -l $lMin,$lMax -s $s -x $x -o $o -n $n -p $p $other $param -t $t -y $y -w $workDir/$w -k $k -a $workDir/$finalFile $workDir/convertopenswath.sptxt >$output 2>&1;";
    }
    
    ###> Launch process
    my $childConvert = fork;
    unless ($childConvert) { # child here
        #>Disconnecting from server
        open STDOUT, '>/dev/null' or die "Can't open /dev/null: $!";
        open STDIN, '</dev/null' or die "Can't open /dev/null: $!";
        open STDERR, '>/dev/null' or die "Can't open /dev/null: $!";
        system $command;
        exit;
    }
    
    #############################################
    ###> Waiting for spectrast2tsv to finish <###
    #############################################
    my $massNumTot=scalar keys %mzHash;
    $massNumTot*=$n;           # ion number per peptide
    my $wait = 1;
    my $errorTxt = '';
    my $massNumber;
    my $prevNbLine=0;
    my $exportFile="$workDir/sortieExport.txt";
    my $waitNb=0;
    my $processFile=($k eq "peakview")? $libraryName."_peakview.tsv" : $libraryName."_openswath.tsv" ;
    while ($wait == 1) {
        sleep 30;
        if ($divID){
            $now=strftime("%s",localtime); # in sec
            $waitTime=strftime("%Hh %Mm %Ss",localtime($now-$startTime-3600));
            $status="Updated $waitTime ago";
            print "<SCRIPT LANGUAGE=\"JavaScript\">$divID.innerHTML=\"\";$divID.innerHTML='$status';</SCRIPT>";
        }
        
        ## loading process
        $massNumber=`cut -f1 $workDir/$processFile | uniq | wc -l` if (-s "$workDir/$processFile");
        $massNumber=($massNumber) ? $massNumber : 1;
        chomp $massNumber;
        my $percent=$massNumber/$massNumTot*100;
        $percent=sprintf("%.0f",$percent);
        $percent='100' if $percent>100;
        if ($processText){
            $loading="<progress value=\"$percent\" max=\"100\" style=\"width:400px\"></progress>";
            print "<SCRIPT LANGAGE=\"JavaScript\">$loadingDivID.innerHTML=\"\";$loadingDivID.innerHTML='<B>Step 1/2 :</B> Export library for PeakView ... <BR><BR> $loading';$loadingSPAN.innerHTML='$percent%';</SCRIPT>" if $k eq 'peakview';
            print "<SCRIPT LANGAGE=\"JavaScript\">$loadingDivID.innerHTML=\"\";$loadingDivID.innerHTML='$processText<BR><BR> $loading';$loadingSPAN.innerHTML='$percent%';</SCRIPT>" if $k eq 'openswath';
        }
        
        open(EXPORTFILE,">>$exportFile");
        print EXPORTFILE $percent;
        close EXPORTFILE;
        
        if(-s $exportFile){
            my $process = `grep -c 'Done.' $exportFile`;
            chomp $process;
            if ($process) { # finished normally
                $wait = 0;
            }
            else {
                $process = `grep -c 'Error' $exportFile`;
                chomp $process;
                if ($process) {
                    $errorTxt=`tail $exportFile`;
                    $wait = 0;
                }
                else{
                    $process=`grep -c 'This modification has not been recognized' $exportFile`;
                    chomp $process;
                    if ($process) {
                        $errorTxt='<BR>A modification has not been recognized, select the modification file in the export form (\'File with modifications delta mass\').';
                        $wait = 0;
                    }
                }
            }
        }
        my $errorExport = `grep -c 'Error' $exportFile`;
        chomp $errorExport;
        if ($errorExport) {
            $errorTxt=`tail $exportFile`;
            $wait = 0;
        }  
        
        if(-s "$workDir/$processFile" && $waitNb==2){
            my @infoLine=split(/\s/,`wc -l $workDir/$processFile`);
            my $nbLine=$infoLine[0];
            $wait=($nbLine>$prevNbLine)? 1 : 0;
            $prevNbLine=$nbLine;
            $waitNb=0;
        }
        $waitNb++;
    }
    
    if ($errorTxt) {
        if ($processText){
            print qq
|</FONT>
<BR><BR><BR><BR><BR>
<FONT class="title2" color="#DD0000">***ERROR: Data preprocessing failed : $errorTxt</FONT>
<BR>
|;
            return;
        }
        else{
            return ('error','error','error',$errorTxt);
            
        }
    }
    
    #############################
    ###> Add +50 to each iRT <###
    #############################
    if ($k eq "peakview") {
        my $numFileLine=`cat $workDir/peakviewfile.txt | wc -l`;
        open(INFILE,"<","$workDir/peakviewfile.txt");
        open(OUTFILE,">","$workDir/$libraryName\_peakview.txt");
        my $numLine=1;
        my $i=1;
        while (my $line=<INFILE>) {
            if ($line=~/^Q1/) {print OUTFILE $line;}
            else{
                my @lineInfo=split(/\t/,$line);
                my $RT=$lineInfo[2]+50;
                my $iRT=$lineInfo[12]+50;
                $lineInfo[2]=$RT;
                $lineInfo[12]=$iRT;
                print OUTFILE join("\t",@lineInfo);
            }
            if ($i==300) {
                $i=0;
                ## loading bar
                my $percent=$numLine/$numFileLine*100;
                $percent=sprintf("%.0f",$percent);
                $percent='100' if $percent>100;
                my $loading="<progress value=\"$percent\" max=\"100\" style=\"width:400px\"></progress>";
                print "<SCRIPT LANGAGE=\"JavaScript\">$loadingDivID.innerHTML=\"\";$loadingDivID.innerHTML='<B>Step 2/2 :</B> Conversion for PeakView ... <BR><BR> $loading';$loadingSPAN.innerHTML='$percent%';</SCRIPT>";
            }
            $i++;
            $numLine++; 
        }
        close OUTFILE;
        close INFILE;
    }
    
    
    print PARAMFILE "OTHER:Use theorical mass" if $other=~/-e/;
    print PARAMFILE "OTHER:Remove duplicate masses from labelling" if $other=~/-d/;
    
    
    my $missingProt;
    if (@validProt && @protList) {
        foreach my $protein (@protList){
            if ("@validProt"!~/$protein/) {
                $missingProt.=',' unless !$missingProt;
                $missingProt.=$protein;
            }
        }
    }
    print PARAMFILE "Desired proteins not found into the library : $missingProt" if $missingProt;
    close PARAMFILE;
    
    my $numValidProt= scalar @validProt if @validProt;
    $numValidProt=($numValidProt)? $numValidProt : 0;
    return ($finalFile,$paramFile,$numValidProt);
}
1;

##########################################################
################### TANDEM XML PARSING ###################
##########################################################
package XTandemXMLHandler; {    ## for .xml and .tandem
    my (%infoRank,%matches,%protListNoSeq,%query2matches,%matchSpectrum);
    my (%modificationList,%bestScoreTarget,%bestScoreDecoy);
    my ($minScore,$maxRank,$dataPath,$dbID,$modifNTer,$type,$dbh,$msFilename,$analysisID,%info,%rankSeq,%prot,$refProtSeq,$refQueryInfo,$refProtDes,$refProtOrg,$refProtLength,$refCharge,$refMassObs,$refMassExp,$refElutionTime,$refNumValid,$refMatchList,$refMaxProtMatch,$refMaxProtScore,$refMaxQueryScore,$refAnalysisMods,$refRankProtMatch,$refProtMW,$refProtDbRank);
    my $i=1;
    
    sub new {
        ($type,$msFilename,$dataPath,$modifNTer,$dbh,$dbID,$minScore,$analysisID,$maxRank,$refProtSeq,$refQueryInfo,$refProtDes,$refProtOrg,$refProtLength,$refCharge,$refMassObs,$refMassExp,$refElutionTime,$refNumValid,$refMatchList,$refMaxProtMatch,$refMaxProtScore,$refMaxQueryScore,$refAnalysisMods,$refRankProtMatch,$refProtMW,$refProtDbRank) = @_;
        my $self = bless ({}, $type);
        
        ###>Fetching parsing rules
        my ($parseRules)=$dbh->selectrow_array("SELECT PARSE_RULES FROM DATABANK,DATABANK_TYPE WHERE DATABANK.ID_DBTYPE=DATABANK_TYPE.ID_DBTYPE AND ID_DATABANK=$dbID");
        my @rules=split(',:,',$parseRules);
        my ($idRule)=($rules[0]=~/ID=(.+)/);
        my ($orgRule,$desRule);
        if ($rules[1]) {
            if ($rules[1]=~/DES=/) {($desRule=$rules[1])=~s/DES=//;}
            else {($orgRule=$rules[1])=~s/ORG=//;}
        }
        if ($rules[2]) {($orgRule=$rules[2])=~s/ORG=//;}
    
        $self->{'idRule'}=$idRule;
        $self->{'desRule'}=$desRule if $desRule;
        $self->{'orgRule'}=$orgRule if $orgRule;
        
        ###>Fetching modifications IDs and names
        my $modificationSth=$dbh->prepare("SELECT PSI_MS_NAME,ID_MODIFICATION,MONO_MASS FROM MODIFICATION ");
        $modificationSth->execute;
        while (my ($modName,$modID,$monoMass)=$modificationSth->fetchrow_array) {
            $modificationList{$modID}="$monoMass&$modName" if ($monoMass && $modName);
        }
        
        return $self;    
    }
    sub start_document {
        my ($self,$doc)=@_;
        $self->{'queryNum'}=0;
    }
    sub end_document {
        my $self=$_[0];
    }
    
    sub start_element {
        my ($self, $element) = @_;
        $self->{'Element'} = $element->{'Name'};
        
        if ($element->{'Name'} eq 'group' && defined $element->{'Attributes'}{'{}rt'}) {            ## one group correspond to one spectra
            $self->{'queryNum'}++;
            $self->{'onGroup'}=1;
            $self->{'charge'}=$element->{'Attributes'}{'{}z'}{'Value'};
            $self->{'massexp'}=$element->{'Attributes'}{'{}mh'}{'Value'};
            $self->{'rttime'}=$element->{'Attributes'}{'{}rt'}{'Value'};
            $self->{'scan'}=$element->{'Attributes'}{'{}id'}{'Value'};
            
        }
        elsif($element->{'Name'} eq 'protein'){	
            $self->{'isOnProtein'}=1;
        }
        elsif($element->{'Name'} eq 'note'){	
            $self->{'isOnNote'}=1;
        }
        elsif($element->{'Name'} eq 'peptide'){
            $self->{'isOnPeptide'}=1;
            $self->{'protLength'}=$element->{'Attributes'}{'{}end'}{'Value'};
        }
        elsif ($element->{'Name'} eq 'domain') {    ## peptide informations
            $self->{'isOnDomain'}=1;
            $self->{'massCalc'}=$element->{'Attributes'}{'{}mh'}{'Value'};
            $self->{'delta'}=$element->{'Attributes'}{'{}delta'}{'Value'};
            $self->{'MIS'}=$element->{'Attributes'}{'{}missed_cleavages'}{'Value'};
            $self->{'sequence'}=$element->{'Attributes'}{'{}seq'}{'Value'};
            my @idDomain=split(/\./,$element->{'Attributes'}{'{}id'}{'Value'});
            $self->{'rankdomain'}="$idDomain[0]"."$idDomain[1]"."$idDomain[-1]";
            my $prev=$element->{'Attributes'}{'{}pre'}{'Value'};
            my @aaPrev=split(//,$prev);
            my $post=$element->{'Attributes'}{'{}post'}{'Value'};
            my @aaPost=split(//,$post);
            my $pos=$element->{'Attributes'}{'{}start'}{'Value'};
            my ($aaPost,$aaPrev);
            if ($aaPrev[-1]=~/\[/) {$aaPrev="-";}
            else{$aaPrev=$aaPrev[-1];}
            if($aaPost[0]=~/\]/){$aaPost="+";}
            else{$aaPost=$aaPost[0];}
            $self->{'protContext'}="$pos,$aaPrev,$aaPost";
            $self->{'start'}=$pos;
            #$self->{'score'}=$element->{'Attributes'}{'{}hyperscore'}{'Value'};
            $self->{'score'}=$element->{'Attributes'}{'{}expect'}{'Value'};
        }
        elsif ($element->{'Name'} eq 'aa') {        ## peptide modifications
            $self->{'mod'}{$element->{'Attributes'}{'{}modified'}{'Value'}}{$element->{'Attributes'}{'{}at'}{'Value'}}=$element->{'Attributes'}{'{}type'}{'Value'};
        }	
    }
    
    sub end_element {
        my ($self, $element) = @_;
        if($element->{'Name'} eq 'note'){	
            $self->{'isOnNote'}=0;
            if ($self->{'isOnProtein'}) {
                my $decoyTag="";
                my $isDecoyProtein;
                if ($self->{'protDes'}=~/reverse_/) {
                    $isDecoyProtein=1;
                    $decoyTag="DECOY_";
                }
                else{$isDecoyProtein=0;}
                my ($identifier)=($self->{'protDes'}=~/$self->{'idRule'}/);
                $identifier="$decoyTag$identifier";
                $self->{'matchDecoyProtein'}=($isDecoyProtein)? 1 : 0;	
                $self->{'matchTargetProtein'}=($isDecoyProtein)? 0 : 1;	
                $self->{'identifier'}=$identifier;
                ($refProtDes->{$self->{'identifier'}})=($self->{'protDes'}=~/$self->{'desRule'}/) if $self->{'desRule'} && !$refProtDes->{$self->{'identifier'}};
                ($refProtOrg->{$self->{'identifier'}})=($self->{'protDes'}=~/$self->{'orgRule'}/) if $self->{'orgRule'} && !$refProtOrg->{$self->{'identifier'}};
                $refProtDbRank->{$self->{'identifier'}}=1;
            }
        }
        elsif($element->{'Name'} eq 'domain'){
            $self->{'isOnDomain'}=0;
            if ($self->{'mod'}) {
                ###> Recovering peptide modifications 
                my %modif;
                foreach my $modifMass (sort {$a <=> $b} keys %{$self->{'mod'}}){
                    my $varMod;
                    my $match=0;
                    my ($numModDB,$unimod);
                    foreach my $modID (keys %modificationList){
                        my ($mass,$modName)=split(/&/,$modificationList{$modID});
                        if ($mass > $modifMass-0.005 && $mass < $modifMass+0.005) {
                            $varMod=$modName;
                            $numModDB=$modID;
                            $match=1;
                            last;
                        }
                    }
                    if (! $match) {
                        my $line;
                        open(MODFILE,"<","$dataPath/Unified_Modification_PeakView.csv") or die ("open: $!"); #find modification in Unified_Modification_PeakView file
                        while ($line=<MODFILE>) {
                            my @lineStrg=split(/;/,$line);
                            if ($lineStrg[11] eq $modifMass || $lineStrg[11] eq ($modifMass=~s/,/./) || (($lineStrg[11] > $modifMass-0.005) && ($lineStrg[11] < $modifMass+0.005)) ) {
                                $varMod=$lineStrg[14];
                                $unimod=$lineStrg[7];
                                last;
                            }
                        }
                        close MODFILE;
                        $numModDB=$dbh->fetchrow_array("SELECT ID_MODIFICATION FROM MODIFICATION WHERE UNIMOD_ACC=$unimod") if $unimod;
                    }
                    $refAnalysisMods->{'V'}{$varMod}{'numModDB'}=$numModDB;
                    $modif{$modifMass}=$varMod;
                }
                
                foreach my $modifMass (sort {$modif{$a} cmp $modif{$b}} keys %{$self->{'mod'}}){
                    my $varMod=$modif{$modifMass};
                    my ($posPeptide,$aaMod);
                    foreach my $posProt (sort {$a <=> $b} keys %{$self->{'mod'}{$modifMass}}){
                        $aaMod=$self->{'mod'}{$modifMass}{$posProt};
                        my $pos;
                        $pos=$posProt-$self->{'start'}+1;   ## peptide's modification position
                        $posPeptide=($posPeptide)? $posPeptide.".".$pos : $pos;
                        
                        if ($modifNTer && $varMod eq 'Acetyl' && $pos==1 && $aaMod ne 'K' ) {
                            $aaMod='N-term';
                        }
                    }
                     if ($aaMod eq 'N-term') {
                        $self->{'varModString'}.=" + $varMod ($aaMod)";
                    }
                    else{$self->{'varModString'}.=" + $varMod ($aaMod:$posPeptide)";}
                    $refAnalysisMods->{'V'}{$varMod}{'specificity'}=$aaMod if not defined $refAnalysisMods->{'V'}{$varMod}{'specificity'};
                    $refAnalysisMods->{'V'}{$varMod}{'specificity'}.=$aaMod if (defined $refAnalysisMods->{'V'}{$varMod}{'specificity'} && $refAnalysisMods->{'V'}{$varMod}{'specificity'}!~/$aaMod/);
                    
                }
                $self->{'actualSequence'}="$self->{'sequence'}$self->{'varModString'}";
                $self->{'mod'}=();
            }
            else{
                $self->{'actualSequence'}=$self->{'sequence'};
                $self->{'varModString'}='';
            }
            $info{"$self->{'sequence'}_$self->{'varModString'}"}{"$self->{'score'},$self->{'massCalc'},$self->{'delta'},$self->{'MIS'},$self->{'matchDecoyProtein'},$self->{'matchTargetProtein'},$self->{'varModString'}"}=1;
            $rankSeq{"$self->{'sequence'}_$self->{'varModString'}"}=$self->{'rankdomain'} unless $rankSeq{"$self->{'sequence'}_$self->{'varModString'}"};
            push @{$prot{"$self->{'sequence'}_$self->{'score'}"}{'DECOY'}},"$self->{'identifier'}&$self->{'protContext'}" if $self->{'matchDecoyProtein'};
            push @{$prot{"$self->{'sequence'}_$self->{'score'}"}{'TARGET'}},"$self->{'identifier'}&$self->{'protContext'}" if $self->{'matchTargetProtein'};
            $self->{'varModString'}='';
        }
        elsif($element->{'Name'} eq 'protein'){
            $i++;
            if ($i==2000) {
                print ".";
                $i=0;
            }
            $self->{'protContext'}='';
            $self->{'protDes'}='';
            $self->{'isOnProtein'}=0;
        }
        elsif ($element->{'Name'} eq 'peptide') {
            if (!$self->{'protSeq'} && not exists $refProtLength->{$self->{'identifier'}}) {
                my $ident;
                if ($self->{'identifier'}=~/DECOY_/){($ident=$self->{'identifier'})=~s/DECOY_//;}
                else{$ident=$self->{'identifier'}}
                $protListNoSeq{$ident}{$analysisID}=1;
            }
            else{
                my $mass=&computeAAMass($self->{'protSeq'},'protein');
                $mass=sprintf "%.2f",$mass; 
                $refProtMW->{$self->{'identifier'}}=$mass if !$refProtMW->{$self->{'identifier'}};
                $refProtLength->{$self->{'identifier'}}=$self->{'protLength'} if !$refProtLength->{$self->{'identifier'}};
                $refProtSeq->{$self->{'identifier'}}=$self->{'protSeq'} unless $refProtSeq->{$self->{'identifier'}};
            }
            $self->{'isOnPeptide'}=0;
            $self->{'matchProtSeq'}=0;
            $self->{'protSeq'}='';
        }
        elsif($element->{'Name'} eq 'group' && $self->{'onGroup'}){
           foreach my $seqMod (sort {$rankSeq{$a} <=>$rankSeq{$b}} keys %info){
                foreach my $info (keys %{$info{$seqMod}}){
                    my ($score,$massCalc,$delta,$mis,$matchDecoyProtein,$matchTargetProtein,$varModString)=split(/,/,$info);
                    my ($sequence,$mod)=split(/_/,$seqMod);
                    $varModString='' unless $varModString;
                    if ($score >= $minScore ){
                        $self->{'hitNum'}++;
                        ###>Computing max query score 
                        if ($matchTargetProtein && ((not defined $refMaxQueryScore->{$self->{'queryNum'}}) || ($refMaxQueryScore->{$self->{'queryNum'}}<$score))) {
                            $refMaxQueryScore->{$self->{'queryNum'}}=$score;
                        }
                        elsif($matchDecoyProtein && ((not defined $refMaxQueryScore->{"-$self->{'queryNum'}"}) || ($refMaxQueryScore->{"-$self->{'queryNum'}"}<$score))){
                            $refMaxQueryScore->{"-$self->{'queryNum'}"}=$score;
                        }
                        $matches{$self->{'hitNum'}}{'massCalc'}=$massCalc-1.00794;              # mh = calculated peptide mass + a proton
                        $matches{$self->{'hitNum'}}{'delta'}=$delta;
                        $matches{$self->{'hitNum'}}{'MIS'}=$mis;
                        $matches{$self->{'hitNum'}}{'sequence'}=$sequence;
                        $matches{$self->{'hitNum'}}{'vmod'}=$varModString;
                        $matches{$self->{'hitNum'}}{'score'}=$score;
                        $matches{$self->{'hitNum'}}{'matchDecoyProtein'}=$matchDecoyProtein;	
                        $matches{$self->{'hitNum'}}{'matchTargetProtein'}=$matchTargetProtein;
                        my $actualSequence="$sequence$varModString";
                        $matches{$self->{'hitNum'}}{'actualSequence'}=$actualSequence;
                        if ($matchDecoyProtein) {
                            $self->{'decoy'}++;
                            $matchSpectrum{"-$self->{'queryNum'}"}{$self->{'hitNum'}}=$score;
                            foreach my $protContext (@{$prot{"$sequence\_$score"}{'DECOY'}}){
                                my ($identifier,$beg)=split(/&/,$protContext);
                                $matches{$self->{'hitNum'}}{'protein'}=$identifier;
                                $refMatchList->{$identifier}{$actualSequence}=1;
                                $infoRank{$self->{'hitNum'}}{$identifier}{$beg}=1;
                                $matches{$self->{'hitNum'}}{'numMatchProt'}++;
                            }
                            $refCharge->{"-$self->{'queryNum'}"}=$self->{'charge'};
                            $refMassExp->{"-$self->{'queryNum'}"}=sprintf('%.4f',$self->{'massexp'}-1.00794);       # mh = parent ion mass (plus a proton) from the spectrum
                            $refMassObs->{"-$self->{'queryNum'}"}=sprintf('%.4f',($matches{$self->{'hitNum'}}{'massCalc'}+($refCharge->{"-$self->{'queryNum'}"}*1.00794))/$refCharge->{"-$self->{'queryNum'}"});
                            (my $rt=$self->{'rttime'})=~s/PT|S//g;
                            $rt=sprintf('%.2f',$rt/60);
                            $refElutionTime->{"-$self->{'queryNum'}"}="sc$self->{'scan'};et$rt";
                        }
                        elsif($matchTargetProtein){
                            $self->{'target'}++;
                            $matchSpectrum{$self->{'queryNum'}}{$self->{'hitNum'}}=$score;
                            foreach my $protContext (@{$prot{"$sequence\_$score"}{'TARGET'}}){
                                my ($identifier,$beg)=split(/&/,$protContext);
                                $matches{$self->{'hitNum'}}{'protein'}=$identifier;
                                $refMatchList->{$identifier}{$actualSequence}=1;
                                $infoRank{$self->{'hitNum'}}{$identifier}{$beg}=1;
                                $matches{$self->{'hitNum'}}{'numMatchProt'}++;
                            }
                            $refCharge->{$self->{'queryNum'}}=$self->{'charge'};
                            $refMassExp->{$self->{'queryNum'}}=sprintf('%.4f',$self->{'massexp'}-1.00794);
                            $refMassObs->{$self->{'queryNum'}}=sprintf('%.4f',($matches{$self->{'hitNum'}}{'massCalc'}+($refCharge->{$self->{'queryNum'}}*1.00794))/$refCharge->{$self->{'queryNum'}});
                            (my $rt=$self->{'rttime'})=~s/PT|S//g;
                            $rt=sprintf('%.2f',$rt/60);
                            $refElutionTime->{$self->{'queryNum'}}="sc$self->{'scan'};et$rt";
                        }
                    }
                    else{$matchSpectrum{$self->{'queryNum'}}=undef;}
                    
                }
            }
            $self->{'onGroup'}=0;
            %info=();
            %prot=();
            %rankSeq=();
        }
        elsif ($element->{'Name'} eq 'bioml') {     ## end of the xml file
            ##Compute query rank
            foreach my $queryNum (keys %matchSpectrum){
                my $rank=0;
                foreach my $hitNum (sort {$matchSpectrum{$queryNum}{$b} cmp $matchSpectrum{$queryNum}{$a} || $a <=> $b} keys %{$matchSpectrum{$queryNum}}){
                    $rank++;
                    $query2matches{$queryNum}{$rank}=$hitNum;
                }
            }
            ##Compute best peptide score
            foreach my $queryNum (keys %query2matches){
                foreach my $rank (keys %{$query2matches{$queryNum}}){
                    next if $rank>$maxRank;
                    my $hitNum=$query2matches{$queryNum}{$rank};
                    if ($matches{$hitNum}{'matchTargetProtein'}) {
                        @{$bestScoreTarget{$matches{$hitNum}{'actualSequence'}}}=($matches{$hitNum}{'score'},$hitNum) if (!$bestScoreTarget{$matches{$hitNum}{'actualSequence'}} || $bestScoreTarget{$matches{$hitNum}{'actualSequence'}}[0] < $matches{$hitNum}{'score'});
                    }
                    if ($matches{$hitNum}{'matchDecoyProtein'}) {
                        @{$bestScoreDecoy{$matches{$hitNum}{'actualSequence'}}}=($matches{$hitNum}{'score'},$hitNum) if (!$bestScoreDecoy{$matches{$hitNum}{'actualSequence'}} || $bestScoreDecoy{$matches{$hitNum}{'actualSequence'}}[0] < $matches{$hitNum}{'score'});
                    }
                }
            }
            ## Storing into %maxProtMatch, %maxProtScore, %rankProtMatch and %queryInfo
            foreach my $queryNum (keys %query2matches){
                my ($rankTarget,$rankDecoy)=(0,0);
                foreach my $rank (sort keys %{$query2matches{$queryNum}}){
                    next if $rank>$maxRank;
                    my $select;
                    my $usedRankID="";
                    my $hitNum=$query2matches{$queryNum}{$rank};
                    if ($matches{$hitNum}{'matchTargetProtein'}) {
                        $select=($bestScoreTarget{$matches{$hitNum}{'actualSequence'}}[1]==$hitNum)? 0:-1;
                        $refNumValid->{$queryNum}++ if ($select==0);
                        $rankTarget++;
                        $usedRankID="$queryNum:$rankTarget";
                        my $dataStrg="MATCH=$matches{$hitNum}{'numMatchProt'},SRK=$rank,SEL=$select,MIS=$matches{$hitNum}{'MIS'},CALC=$matches{$hitNum}{'massCalc'},SC=$matchSpectrum{$queryNum}{$hitNum},DELT=$matches{$hitNum}{'delta'},SEQ=$matches{$hitNum}{'sequence'},";
                        $dataStrg.="VMOD=$matches{$hitNum}{'vmod'}," if $matches{$hitNum}{'vmod'};
                        push @{$refQueryInfo->{$queryNum}},$dataStrg;
                    }
                    elsif($matches{$hitNum}{'matchDecoyProtein'}){
                        $select=($bestScoreDecoy{$matches{$hitNum}{'actualSequence'}}[1]==$hitNum)? 0:-1;
                        $refNumValid->{$queryNum}++ if ($select==0);
                        $rankDecoy++;
                        $usedRankID="$queryNum:$rankDecoy";
                        my $dataStrg="MATCH=$matches{$hitNum}{'numMatchProt'},SRK=$rank,SEL=$select,MIS=$matches{$hitNum}{'MIS'},CALC=$matches{$hitNum}{'massCalc'},SC=$matchSpectrum{$queryNum}{$hitNum},DELT=$matches{$hitNum}{'delta'},SEQ=$matches{$hitNum}{'sequence'},";
                        $dataStrg.="VMOD=$matches{$hitNum}{'vmod'}," if $matches{$hitNum}{'vmod'};
                        push @{$refQueryInfo->{$queryNum}},$dataStrg;
                    }
                    foreach my $identifier (keys %{$infoRank{$hitNum}}) {
                        next unless $refMatchList->{$identifier};
                        $refRankProtMatch->{$usedRankID}{$identifier}=();
                        $refMaxProtScore->{$identifier}+=$matches{$hitNum}{'score'} if $select==0;
                        foreach my $beg (keys %{$infoRank{$hitNum}{$identifier}}){
                            push @{$refRankProtMatch->{$usedRankID}{$identifier}} , $beg;
                            $refMaxProtMatch->{$identifier}{$matches{$hitNum}{'actualSequence'}}++;
                        }
                    }
                }
            }
            foreach my $identifier (keys %{$refMatchList}) {
                delete($refMatchList->{$identifier}) unless $refMaxProtScore->{$identifier};
                delete($refProtSeq->{$identifier}) unless $refMaxProtScore->{$identifier};
            }
            if (scalar keys %protListNoSeq) {
                my @analysisID;
                push @analysisID,$analysisID;
                &promsMod::getProtInfo('silent',$dbh,$dbID,\@analysisID,$refProtDes,$refProtMW,$refProtOrg,$refProtLength,undef,\%protListNoSeq);
            }
        }
        delete ($self->{'Elements'}{$element->{'Name'}});
        $self->{'Element'} ='';			
    }
    
    sub characters {
        my ($self, $element) = @_;
        if ($self->{'isOnNote'} && $self->{'isOnProtein'}) {
            $self->{'protDes'}.=$element->{'Data'};
        }
        elsif($self->{'isOnPeptide'} && !$self->{'isOnDomain'}){
            ($self->{'protSeq'}.=$element->{'Data'})=~s/\s//g;
        }
    }
    
    ###> Function that computes the mass of a peptide or a protein given its sequence
    sub computeAAMass {
        my ($sequence,$type)=@_;
        my (%massAAave,%massATave);
        if($type eq 'protein'){# Average mass for protein mass computation
            %massAAave=&promsConfig::getMassAAave; # Average mass, needed for protein mass calculation
            %massATave=&promsConfig::getMassATave; # Average mass, needed for protein mass calculation
        }elsif($type eq 'peptide'){# Isotopic mass for Peptide mass computation
            %massAAave=&promsConfig::getMassAAmono; # Average mass, needed for peptide mass calculation
            %massATave=&promsConfig::getMassATMono; # Average mass, needed for peptide mass calculation
        }
    
        my %countAA;
        my $mass=0.000000000;
        foreach my $aa (split(//,uc($sequence)) ) {
            $countAA{$aa}++;
        }
        foreach my $aa (keys %countAA) {
            $mass+=($massAAave{$aa}*$countAA{$aa}) if $massAAave{$aa}; # some characters are not amino-acids
        }
        $mass+=($massATave{H}+$massATave{H}+$massATave{O}) if $mass; # H(Nter) + OH(Cter)
        return $mass;
    }
}
1;
##################################################################
################### TANDEM.PEP.XML XML PARSING ###################
##################################################################
package XTandemPEPXMLHandler; { ## only for .tandem.pep.xml
    my %infoRank;
    my (%matches,%scorePep,$protDBSeq,%query2matches,%matchSpectrum);
    my (%modificationList,%bestScoreTarget,%bestScoreDecoy);
    my ($refProtSeq,$type,$msFilename,$dataPath,$maxRank,$dbh,$dbID,$minScore,$analysisID,$refQueryInfo,$refProtDes,$refProtOrg,$refProtLength,$refCharge,$refMassObs,$refMassExp,$refElutionTime,$refNumValid,$refMatchList,$refMaxProtMatch,$refMaxProtScore,$refMaxQueryScore,$refAnalysisMods,$refRankProtMatch,$refProtMW,$refProtDbRank);
    my $i=1;

    sub new {
        ($type,$msFilename,$dataPath,$protDBSeq,$dbh,$dbID,$minScore,$maxRank,$analysisID,$refProtSeq,$refQueryInfo,$refProtDes,$refProtOrg,$refProtLength,$refCharge,$refMassObs,$refMassExp,$refElutionTime,$refNumValid,$refMatchList,$refMaxProtMatch,$refMaxProtScore,$refMaxQueryScore,$refAnalysisMods,$refRankProtMatch,$refProtMW,$refProtDbRank) = @_;
        my $self = bless ({}, $type);
        $self->{'protDBSeq'}=$protDBSeq;
        
        ###>Fetching parsing rules
        my ($parseRules)=$dbh->selectrow_array("SELECT PARSE_RULES FROM DATABANK,DATABANK_TYPE WHERE DATABANK.ID_DBTYPE=DATABANK_TYPE.ID_DBTYPE AND ID_DATABANK=$dbID");
        my @rules=split(',:,',$parseRules);
        my ($idRule)=($rules[0]=~/ID=(.+)/);
        my ($orgRule,$desRule);
        if ($rules[1]) {
            if ($rules[1]=~/DES=/) {($desRule=$rules[1])=~s/DES=//;}
            else {($orgRule=$rules[1])=~s/ORG=//;}
        }
        if ($rules[2]) {($orgRule=$rules[2])=~s/ORG=//;}
    
        $self->{'idRule'}=$idRule;
        $self->{'desRule'}=$desRule if $desRule;
        $self->{'orgRule'}=$orgRule if $orgRule;
        
        ###>Fetching modifications IDs and names
        my $modificationSth=$dbh->prepare("SELECT PSI_MS_NAME,ID_MODIFICATION,MONO_MASS FROM MODIFICATION ");
        $modificationSth->execute;
        while (my ($modName,$modID,$monoMass)=$modificationSth->fetchrow_array) {
            $modificationList{$modID}="$monoMass&$modName" if ($monoMass && $modName);
        }
        return $self;    
    }
    sub start_document {
        my ($self, $doc) = @_;
        $self->{'queryNum'}=0;
    }
    sub end_document {
        my $self = $_[0];
    }
    
    sub start_element {
        my ($self, $element) = @_;
        $self->{'Element'} = $element->{'Name'};
        
        if ($element->{'Name'} eq 'spectrum_query') {
            my @queryInf=split(/\./,$element->{'Attributes'}{'{}spectrum'}{'Value'});
            $self->{'queryNum'}++;
            $self->{'charge'}=$element->{'Attributes'}{'{}assumed_charge'}{'Value'};
            $self->{'massexp'}=$element->{'Attributes'}{'{}precursor_neutral_mass'}{'Value'};
            $self->{'rttime'}=sprintf('%.2f',$element->{'Attributes'}{'{}retention_time_sec'}{'Value'}/60);
            $self->{'scan'}=$element->{'Attributes'}{'{}start_scan'}{'Value'};
        }
        elsif ($element->{'Name'} eq 'search_hit') {
            $self->{'hitNum'}++;
            my $entry=$element->{'Attributes'}{'{}protein'}{'Value'};
            my $decoyTag="";
            my $isDecoyProtein;
            if ($entry=~/reverse_/) {
                $isDecoyProtein=1;
                $decoyTag="DECOY_";
            }
            else{$isDecoyProtein=0;}
            my $identifier="$decoyTag$entry";
            $self->{'numMatchProt'}=$element->{'Attributes'}{'{}num_tot_proteins'}{'Value'};
            
            $self->{'masscalc'}=$element->{'Attributes'}{'{}calc_neutral_pep_mass'}{'Value'};
            $self->{'delta'}=$element->{'Attributes'}{'{}massdiff'}{'Value'};
            $self->{'MIS'}=$element->{'Attributes'}{'{}num_missed_cleavages'}{'Value'};
            $self->{'sequence'}=$element->{'Attributes'}{'{}peptide'}{'Value'};
            $self->{'matchDecoyProtein'}=($isDecoyProtein)? 1 : 0;	
            $self->{'matchTargetProtein'}=($isDecoyProtein)? 0 : 1;	
            $self->{'identifier'}=$identifier;
            
            ($refProtDes->{$identifier})=($element->{'Attributes'}{'{}protein_descr'}{'Value'}=~/(.*)\sOS=/);
            ($refProtOrg->{$identifier})=($element->{'Attributes'}{'{}protein_descr'}{'Value'}=~/$self->{'orgRule'}/) if $self->{'orgRule'} && !$refProtOrg->{$identifier};
            $refProtDbRank->{$self->{'identifier'}}=1;
            $refProtLength->{$identifier}=length($self->{'protDBSeq'}{$identifier});
            my $mass=&computeAAMass($self->{'protDBSeq'}{$identifier},'protein');
            $mass=sprintf "%.2f",$mass; 
            $refProtMW->{$identifier}=$mass;
            
            $self->{'before'}=($element->{'Attributes'}{'{}peptide_prev_aa'}{'Value'})?  $element->{'Attributes'}{'{}peptide_prev_aa'}{'Value'} : '-';
            $self->{'after'}=($element->{'Attributes'}{'{}peptide_next_aa'}{'Value'})? $element->{'Attributes'}{'{}peptide_next_aa'}{'Value'} : '-';
            
            if ($self->{'protDBSeq'}{$identifier}=~/(\w?)$self->{'sequence'}(\w?)/){
                while ($self->{'protDBSeq'}{$identifier}=~/(\w?)$self->{'sequence'}(\w?)/g) {
                    my ($pepBeg,$pepEnd)=($1,$2);
                    my $startPos=$-[0]+1;
                    my $endPos=$+[0];
                    if ($pepBeg) {$startPos++;}
                    else{$pepBeg='-';}
                    if ($pepEnd) {$endPos--;}
                    else{$pepEnd='-';}
                    push @{$self->{'protContext'}}, "$startPos,$pepBeg,$pepEnd:$identifier";
                }
            }
            else{push @{$self->{'protContext'}}, "-,$self->{'before'},$self->{'after'}:$identifier";}
        }
        elsif($element->{'Name'} eq 'alternative_protein'){
            my $entry=$element->{'Attributes'}{'{}protein'}{'Value'};
            my $decoyTag="";
            my $isDecoyProtein;
            if ($entry=~/reverse_/) {
                $isDecoyProtein=1;
                $decoyTag="DECOY_";
            }
            else{$isDecoyProtein=0;}
            my $identifier="$decoyTag$entry";
            push @{$self->{'altidentifier'}},$identifier;
            ($refProtDes->{$identifier})=($element->{'Attributes'}{'{}protein_descr'}{'Value'}=~/(.*)\sOS=/);
            ($refProtOrg->{$identifier})=($element->{'Attributes'}{'{}protein_descr'}{'Value'}=~/$self->{'orgRule'}/) if $self->{'orgRule'} && !$refProtOrg->{$identifier};
            $refProtDbRank->{$identifier}=1;
            $refProtLength->{$identifier}=length($self->{'protDBSeq'}{$identifier});
            my $mass=&computeAAMass($self->{'protDBSeq'}{$identifier},'protein');
            $mass=sprintf "%.2f",$mass; 
            $refProtMW->{$identifier}=$mass;
            
            if ($self->{'protDBSeq'}{$identifier}=~/(\w?)$self->{'sequence'}(\w?)/){
                while ($self->{'protDBSeq'}{$identifier}=~/(\w?)$self->{'sequence'}(\w?)/g) {
                    my ($pepBeg,$pepEnd)=($1,$2);
                    my $startPos=$-[0]+1;
                    my $endPos=$+[0];
                    if ($pepBeg) {$startPos++;}
                    else{$pepBeg='-';}
                    if ($pepEnd) {$endPos--;}
                    else{$pepEnd='-';}
                    push @{$self->{'protContext'}}, "$startPos,$pepBeg,$pepEnd:$identifier";
                }
            }
            else{push @{$self->{'protContext'}}, "-,$self->{'before'},$self->{'after'}:$identifier";}
        }
        elsif($element->{'Name'} eq 'search_score'){
            if ($element->{'Attributes'}{'{}name'}{'Value'} eq 'hyperscore') {
                $self->{'score'}=$element->{'Attributes'}{'{}value'}{'Value'};
                if ($self->{'score'} >= $minScore){
                    $matchSpectrum{$self->{'queryNum'}}{$self->{'hitNum'}}=$self->{'score'};
                    $scorePep{$self->{'hitNum'}}=$self->{'score'};
                    if ($self->{'matchTargetProtein'} && ((not defined($refMaxQueryScore->{$self->{'queryNum'}})) || (($refMaxQueryScore->{$self->{'queryNum'}})<($element->{'Attributes'}{'{}value'}{'Value'})))) {
                        $refMaxQueryScore->{$self->{'queryNum'}}=$element->{'Attributes'}{'{}value'}{'Value'} ;
                    }
                    elsif($self->{'matchDecoyProtein'} && ((not defined($refMaxQueryScore->{"-$self->{'queryNum'}"})) || (($refMaxQueryScore->{"-$self->{'queryNum'}"})<($element->{'Attributes'}{'{}value'}{'Value'})))){
                        $refMaxQueryScore->{"-$self->{'queryNum'}"}=$element->{'Attributes'}{'{}value'}{'Value'} ;
                    }
                }
                else{
                    $matchSpectrum{$self->{'queryNum'}}=undef;
                }
            }
        }
        elsif($element->{'Name'} eq 'mod_aminoacid_mass'){
            $self->{'mod'}{$element->{'Attributes'}{'{}mass'}{'Value'}}{$element->{'Attributes'}{'{}position'}{'Value'}}=1;
        }
    }
    
    sub end_element {
        my ($self, $element) = @_;
        if ($element->{'Name'} eq 'modification_info' && exists $self->{'mod'}) {
            ###>Recovering peptides modifications
            my %modifList;
            my @sequence=split(//,$self->{'sequence'});
            my %massAAave=&promsConfig::getMassAAmono;
            foreach my $mod ( keys %{$self->{'mod'}}){
                foreach my $pos (keys %{$self->{'mod'}{$mod}}){
                    my $aaMod;
                    if ($pos eq '[') {  # N-term modif
                        $pos='=';
                        $aaMod=$sequence[0];
                    }		
                    elsif ($pos eq ']'){# C-term modif
                        $pos='-';
                        $aaMod=$sequence[-1];
                    }		
                    else{$aaMod=$sequence[$pos-1];}
                    
                    my $varMod;
                    my $aaMass=$massAAave{$aaMod};
                    my $modifMass=$mod-$aaMass;
                    my $match=0;
                    my ($numModDB,$unimod);
                    foreach my $modID (keys %modificationList){
                        my ($mass,$modName)=split(/&/,$modificationList{$modID});
                        if ($mass > $modifMass-0.005 && $mass < $modifMass+0.005) {
                            $varMod=$modName;
                            $numModDB=$modID;
                            $match=1;
                            last;
                        }
                    }
                    if (! $match) {
                        my $line;
                        open(MODFILE,"<","$dataPath/Unified_Modification_PeakView.csv") or die ("open: $!"); #find modification in Unified_Modification_PeakView file
                        while ($line=<MODFILE>) {
                            my @lineStrg=split(/;/,$line);
                            if ($lineStrg[11] eq $modifMass || $lineStrg[11] eq ($modifMass=~s/,/./) || (($lineStrg[11] > $modifMass-0.005) && ($lineStrg[11] < $modifMass+0.005)) ) {
                                $varMod=$lineStrg[14];
                                $unimod=$lineStrg[7];
                                last;
                            }
                        }
                        close MODFILE;
                        $numModDB=$dbh->fetchrow_array("SELECT ID_MODIFICATION FROM MODIFICATION WHERE UNIMOD_ACC=$unimod") if $unimod;
                    }
                    $modifList{$varMod}{$aaMod}{$pos}=1;
                    $refAnalysisMods->{'V'}{$varMod}{'numModDB'}=$numModDB;
                }
            }
            foreach my $varMod (sort {$a cmp $b} keys %modifList){
                my $posPeptide;
                foreach my $aaMod (sort {$a cmp $b} keys %{$modifList{$varMod}}){
                    foreach my $pos (sort {$a <=> $b} keys %{$modifList{$varMod}{$aaMod}}){
                        $posPeptide=($posPeptide)? $posPeptide.".".$pos : $pos;
                    }
                    $self->{'varModString'}.=" + $varMod ($aaMod:$posPeptide)";
                    $refAnalysisMods->{'V'}{$varMod}{'specificity'}=$aaMod if not defined $refAnalysisMods->{'V'}{$varMod}{'specificity'};
                    $refAnalysisMods->{'V'}{$varMod}{'specificity'}.=$aaMod if (defined $refAnalysisMods->{'V'}{$varMod}{'specificity'} && $refAnalysisMods->{'V'}{$varMod}{'specificity'}!~/$aaMod/);
                }
            }
            $self->{'actualSequence'}="$self->{'sequence'}$self->{'varModString'}";
            $self->{'mod'}=();
        }
        elsif($element->{'Name'} eq 'search_hit'){
            if (defined $matchSpectrum{$self->{'queryNum'}}) {
                my $actualSequence=($self->{'actualSequence'})? $self->{'actualSequence'} : $self->{'sequence'};
                $matches{$self->{'hitNum'}}{'masscalc'}=$self->{'masscalc'};
                $matches{$self->{'hitNum'}}{'delta'}=$self->{'delta'};
                $matches{$self->{'hitNum'}}{'MIS'}=$self->{'MIS'};
                $matches{$self->{'hitNum'}}{'sequence'}=$self->{'sequence'};
                $matches{$self->{'hitNum'}}{'matchDecoyProtein'}=$self->{'matchDecoyProtein'};	
                $matches{$self->{'hitNum'}}{'matchTargetProtein'}=$self->{'matchTargetProtein'};
                $matches{$self->{'hitNum'}}{'protein'}=$self->{'identifier'};
                $matches{$self->{'hitNum'}}{'actualSequence'}=$actualSequence;
                $matches{$self->{'hitNum'}}{'vmod'}=$self->{'varModString'} if $self->{'varModString'};
                $refMatchList->{$self->{'identifier'}}{$actualSequence}=1;
                
                if (@{$self->{'altidentifier'}}) {
                    foreach my $identifier (@{$self->{'altidentifier'}}){
                        $refMatchList->{$identifier}{$actualSequence}=1;
                    }
                }
                foreach my $protContext (@{$self->{'protContext'}}){
                    my ($beg,$identifier)=split(/:/,$protContext);
                    $infoRank{$self->{'hitNum'}}{$identifier}{$beg}=1;
                }
                if ($self->{'matchDecoyProtein'}) {
                    $refCharge->{"-$self->{'queryNum'}"}=$self->{'charge'};
                    $refMassExp->{"-$self->{'queryNum'}"}=$self->{'massexp'};
                    $refMassObs->{"-$self->{'queryNum'}"}=sprintf('%.4f',($matches{$self->{'hitNum'}}{'masscalc'}+($refCharge->{"-$self->{'queryNum'}"}*1.00794))/$refCharge->{"-$self->{'queryNum'}"});
                    $refElutionTime->{"-$self->{'queryNum'}"}="sc$self->{'scan'};et$self->{'rttime'}";
                }
                elsif($self->{'matchTargetProtein'}){
                    $refCharge->{$self->{'queryNum'}}=$self->{'charge'};
                    $refMassExp->{$self->{'queryNum'}}=$self->{'massexp'};
                    $refMassObs->{$self->{'queryNum'}}=sprintf('%.4f',($matches{$self->{'hitNum'}}{'masscalc'}+($refCharge->{$self->{'queryNum'}}*1.00794))/$refCharge->{$self->{'queryNum'}});
                    $refElutionTime->{$self->{'queryNum'}}="sc$self->{'scan'};et$self->{'rttime'}";
                }
            }	
            $self->{'protContext'}=();
            $self->{'actualSequence'}='';
            $self->{'varModString'}='';
            $self->{'altidentifier'}=();
            $i++;
            if ($i==2000) {
                print ".";
                $i=0;
            }
        }
        elsif ($element->{'Name'} eq 'msms_pipeline_analysis') {
            ##Compute query rank
            foreach my $queryNum (keys %matchSpectrum){
                my $rank=0;
                foreach my $hitNum (sort {$matchSpectrum{$queryNum}{$b} cmp $matchSpectrum{$queryNum}{$a}} keys %{$matchSpectrum{$queryNum}}){
                    $rank++;
                    $query2matches{$queryNum}{$rank}=$hitNum;
                }
            }
            ##>Compute max peptide score
            foreach my $queryNum (keys %query2matches){
                foreach my $rank (keys %{$query2matches{$queryNum}}){
                    next if $rank>$maxRank;
                    my $hitNum=$query2matches{$queryNum}{$rank};
                    if ($matches{$hitNum}{'matchTargetProtein'}) {
                        @{$bestScoreTarget{$matches{$hitNum}{'actualSequence'}}}=($scorePep{$hitNum},$hitNum) if (!$bestScoreTarget{$matches{$hitNum}{'actualSequence'}} || $bestScoreTarget{$matches{$hitNum}{'actualSequence'}}[0] < $scorePep{$hitNum});
                    }
                    if ($matches{$hitNum}{'matchDecoyProtein'}) {
                        @{$bestScoreDecoy{$matches{$hitNum}{'actualSequence'}}}=($scorePep{$hitNum},$hitNum) if (!$bestScoreDecoy{$matches{$hitNum}{'actualSequence'}} || $bestScoreDecoy{$matches{$hitNum}{'actualSequence'}}[0] < $scorePep{$hitNum});
                    }
                }
            }
            ###>Storing into %maxProtScore, %maxProtMatch, %rankProtMatch, and %queryInfo
            foreach my $queryNum (keys %query2matches){
                my ($rankTarget,$rankDecoy)=(0,0);
                foreach my $rank (sort keys %{$query2matches{$queryNum}}){
                    next if $rank>$maxRank;
                    my $usedRankID="";
                    my $hitNum=$query2matches{$queryNum}{$rank};
                    if ($matches{$hitNum}{'matchTargetProtein'}) {
                        my $select=($bestScoreTarget{$matches{$hitNum}{'actualSequence'}}[1]==$hitNum)? 0:-1;
                        $refNumValid->{$queryNum}=1 if (!$refNumValid->{$queryNum} && $select==0);
                        $rankTarget++;
                        $usedRankID="$queryNum:$rankTarget";
                        my $dataStrg="MATCH=$self->{'numMatchProt'},SRK=$rank,SEL=$select,MIS=$matches{$hitNum}{'MIS'},CALC=$matches{$hitNum}{'masscalc'},SC=$matchSpectrum{$queryNum}{$hitNum},DELT=$matches{$hitNum}{'delta'},SEQ=$matches{$hitNum}{'sequence'},";
                        $dataStrg.="VMOD=$matches{$hitNum}{'vmod'}," if $matches{$hitNum}{'vmod'};
                        push @{$refQueryInfo->{$queryNum}},$dataStrg;
                    }
                    elsif($matches{$hitNum}{'matchDecoyProtein'}){
                        my $select=($bestScoreDecoy{$matches{$hitNum}{'actualSequence'}}[1]==$hitNum)? 0:-1;
                        $refNumValid->{"-$queryNum"}=1 if (!$refNumValid->{"-$queryNum"} && $select==0);
                        $rankDecoy++;
                        $usedRankID="-$queryNum:$rankDecoy";
                        my $dataStrg="MATCH=$self->{'numMatchProt'},SRK=$rank,SEL=$select,MIS=$matches{$hitNum}{'MIS'},CALC=$matches{$hitNum}{'masscalc'},SC=$matchSpectrum{$queryNum}{$hitNum},DELT=$matches{$hitNum}{'delta'},SEQ=$matches{$hitNum}{'sequence'},";
                        $dataStrg.="VMOD=$matches{$hitNum}{'vmod'}," if $matches{$hitNum}{'vmod'};
                        push @{$refQueryInfo->{"-$queryNum"}},$dataStrg;
                    }
                    foreach my $identifier (keys %{$infoRank{$hitNum}}) {
                        next unless $refMatchList->{$identifier};
                        $refRankProtMatch->{$usedRankID}{$identifier}=();
                        $refMaxProtScore->{$identifier}+=$scorePep{$hitNum};
                        foreach my $beg (keys %{$infoRank{$hitNum}{$identifier}}){
                            push @{$refRankProtMatch->{$usedRankID}{$identifier}} , $beg;
                            $refMaxProtMatch->{$identifier}{$matches{$hitNum}{'actualSequence'}}++;
                        }
                        $refProtSeq->{$identifier}=$self->{'protDBSeq'}{$identifier} unless $refProtSeq->{$identifier};
                    }
                }
            }
        }
        delete ($self->{'Elements'}{$element->{'Name'}});
        $self->{'Element'} ='';			
    }
    
    ###> Function that computes the mass of a peptide or a protein given its sequence
    sub computeAAMass {
        my ($sequence,$type)=@_;
        my (%massAAave,%massATave);
        if($type eq 'protein'){# Average mass for protein mass computation
            %massAAave=&promsConfig::getMassAAave; # Average mass, needed for protein mass calculation
            %massATave=&promsConfig::getMassATave; # Average mass, needed for protein mass calculation
        }elsif($type eq 'peptide'){# Isotopic mass for Peptide mass computation
            %massAAave=&promsConfig::getMassAAmono; # Average mass, needed for peptide mass calculation
            %massATave=&promsConfig::getMassATMono; # Average mass, needed for peptide mass calculation
        }

        my %countAA;
        my $mass=0.000000000;
        foreach my $aa (split(//,uc($sequence)) ) {
            $countAA{$aa}++;
        }
        foreach my $aa (keys %countAA) {
            $mass+=($massAAave{$aa}*$countAA{$aa}) if $massAAave{$aa}; # some characters are not amino-acids
        }
        $mass+=($massATave{H}+$massATave{H}+$massATave{O}) if $mass; # H(Nter) + OH(Cter)

        return $mass;
    }
}
1;

####>Revision history<####
# 1.1.4 Encapsulated the promsPath within methods to avoid path retrieving issues when calling them by command line (VS 22/11/18)
# 1.1.3 Minor modification on the getProtInfo call (VS 16/11/2018)
# 1.1.2 Minor modif to handle library export errors (MLP 17/04/18)
# 1.1.1 Minor modif (MLP 10/01/18)
# 1.1.0 Add export function (MLP 06/12/17)
# 1.0.2 Minor correction (MLP 02/08/17)
# 1.0.1 Add modification to recover protein sequence (used to write analysis.fasta) (MLP 04/05/17)
# 1.0.0 Created (MLP 05/04/2017)