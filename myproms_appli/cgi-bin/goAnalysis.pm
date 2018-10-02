#############################################################################
# goAnalysis.pm               1.1.5                                         #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                   #
# Contact: myproms@curie.fr                                                 #
# Performs term enrichment analysis on a protein set using GO::TermFinder,  #
# Data provided by goAnalysis.cgi                                           #
# Drawing graph using GraphViz                                              #
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
package goAnalysis;

use GO::TermFinder;
use GraphViz;
use strict;
#use promsConfig;
#use promsMod;
use FileHandle;
#use LWP::UserAgent;
use File::Path qw(rmtree); # remove_tree
use File::Copy;

#use subs 'print';
my $verbose = 1;

#########################################
#### Fetching parameters (arguments) ####
#########################################
sub termFinder{
    my (%param) = @_;

    my $dbh = $param{dbh} or die 'Missing DBH';
    my $annotation = $param{annotation} or die 'Missing annotation';
    my %ontologies = %{$param{ontologies}} or die 'Missing ontology';
    my %protIDs = %{$param{data}} or die 'Missing data';
    my @aspects = @{$param{aspects}} or die 'Missing at least 1 aspect'; # (P,C,F)
    my $fileName = $param{fileName} or die 'Missing filename';
    my $drawGraph = (defined $param{drawGraph})?$param{drawGraph} : 2; # 0 -> no graph, 1 -> graph without unsignificant nodes, 2 -> complete graph
    my $threshold = ($param{threshold})?$param{threshold} : 0; # [0..1] if p-value criteria, [0..100] if FDR
    my $minPep = ($param{minPep})? $param{minPep}:1; # ignore proteins with less than $minPep peptides
    my @population = @{$param{population}}; # background population (list of protein NAMES that are annotated in GOA file, DB ids are useless)
    my $protNb = ($param{protNb})?$param{protNb}:undef;
    my $criteria = $param{criteria} or die 'Missing criteria'; #FDR or p-value
    my $method = ($param{method})?$param{method}:'none'; # FDR or p-value calculation method (FDRsimulation, bonferroni or none)
    my $slimOBO = ($param{slim})? $param{slim} : 0;
    $verbose = (defined $param{verbose})?$param{verbose}:1;

    my $pvalStrg; # to get the correct hash key name for p-values
    my $correction;
    if($method eq 'bonferroni'){
        $pvalStrg = 'CORRECTED_PVALUE';
        $correction = 'bonferroni';
    } else { # no multitest correction or FDR method
        $pvalStrg = 'PVALUE';
        $correction = 'none';
    }
    my $FDRsimulation = ($method eq 'FDRsimulation')? 1 : 0;
    #$minPep = 1 unless $minPep;

    my %promsPath=&promsConfig::getServerInfo('no_user');

    ####################################################
    #### Creating job file                          ####
    #### (readed by startGOAnalysis.cgi monitoring) ####
    ####################################################

    # Create job file (will be erased at the end)
    unless(-d "$promsPath{'tmp'}/GO"){
        mkdir "$promsPath{'tmp'}/GO";
    }
    mkdir "$promsPath{'tmp'}/GO/$fileName" or die $!;
    my $errorLogFile = "$promsPath{'tmp'}/GO/$fileName/errorLog.txt";

    ####################################
    #### SIG variable redefinitions ####
    ####################################
    local %SIG;

    # Script time out #
    $SIG{ALRM} = sub {die 'Timed out'};
    alarm 10000;

    # Printing warnings into errorLog (can be displayed within editGOAnalysis.cgi) #
    #$SIG{__WARN__} = sub {
    #    open WARNLOG, ">>$errorLogFile";
    #    print WARNLOG $_[0];
    #    close WARNLOG;
    #    warn $_[0];
    #};

    # Erasing files before any die #
    $SIG{__DIE__} = sub {
        rmtree("$promsPath{'tmp'}/GO/$fileName") if -d "$promsPath{'tmp'}/GO/$fileName";
        die $_[0];
    };

    #####################################
    #### Check if defined background ####
    #####################################

    my @additionalParameters;
    if(scalar @population){
        push @additionalParameters, (population => \@population);
    }

    ############################################################
    #### Check if number of proteins in organism is defined ####
    ############################################################

    if($protNb){
        push @additionalParameters, (totalNumGenes => $protNb);
    }



    ###########################################
    #### Annotation correction if slim OBO ####
    ###########################################
    if($slimOBO){
        subPrint("Annotation correction for slim...");
        &correctAnnotationForSlim($dbh,$annotation,\%ontologies,\@aspects);
        subPrint("OK <BR>\n");
    }


    ##############################################################
    #### Starting enrichment analysis for each aspect (P,C,F) ####
    ##############################################################

    my %drawnNodes;
    my @protList = keys %protIDs;
    foreach my $aspect (@aspects){

        my $ontology = $ontologies{$aspect};
        next unless &checkOntology($ontology);

        #-----------------------------------#
        # Checking for unannotated proteins #
        #-----------------------------------#
        my @unannotatedProteins;
        foreach my $protein (@protList){
            unless ($annotation->nameIsAnnotated( name => $protein , aspect => $aspect )){
                push @unannotatedProteins, $protIDs{$protein};
            }
        }

        #----------------------------#
        # Starting TermFinder module #
        #----------------------------#
        subPrint("Starting enrichment analysis for ".&getAspectName($aspect)."...");

        # Enrichment test
        my $termFinder;
        my @pvalues;
        {
            eval {
                # Redirecting STDOUT, GO::TermFinder prints some warnings in it
                local *STDOUT;
                open STDOUT, ">>$errorLogFile";
                local *STDERR = *STDOUT; # redirecting true warnings into STDOUT > errorLogFile
                $termFinder = GO::TermFinder->new(annotationProvider=> $annotation,
                                             ontologyProvider  => $ontology,
                                             aspect            => $aspect,
                                             @additionalParameters
                                            );
                @pvalues = $termFinder->findTerms(genes => \@protList, calculateFDR => $FDRsimulation, correction => $correction);
            } or do {
                subPrint("<FONT COLOR=\"red\">Error: ".$@."</FONT>\n");
                return 0;
            }
        }

        subPrint(" OK<BR>\n");

        #-----------------#
        # Preparing graph #
        #-----------------#
        subPrint("Drawing graph...");
        my $graph = GraphViz->new(
                        node => {shape => 'box', fontsize => 9, height => 0.2, width => 0.5}, # node default attributes
                        );

        # Adding root node
        $graph->add_node(getRootID($aspect) , label=>getAspectName($aspect) , shape => 'ellipse');
        $drawnNodes{getRootID($aspect)} = 1;

        #-----------------------#
        # Preparing result file #
        #-----------------------#
        open(RES,">$promsPath{'tmp'}/GO/$fileName/results_$aspect.txt"); RES->autoflush(1);
        print RES "### ",scalar @protList,"\t",$termFinder->totalNumGenes,"\n";

        #------------------------------#
        # Browsing each term (p-value) #
        #------------------------------#
        my $pvalNum = scalar @pvalues; # get total number of pvalues (i.e. terms) to calculate FDR by BH method
        my $i = 0; # p-value rank counter for FDR calculation (BH)
        my @significantPvalues;
        foreach my $pval (sort {$a->{$pvalStrg} <=> $b->{$pvalStrg}} @pvalues){
            $i++;
            my $node = $pval->{'NODE'};
            next if $node->term eq 'unannotated'; # Skip term if unannotated group

            #-----------------------------------#
            # Checking FDR or P-value threshold #
            #-----------------------------------#
            # P-value threshold
            if($criteria eq 'pval' and $threshold < 1){
                last if ($pval->{$pvalStrg} > $threshold);
            }
            # or FDR control
            elsif ($criteria eq 'FDR' and $threshold < 100){
                # Benjamini/Hochberg method
                if($method eq 'BH'){
                    last if $pval->{$pvalStrg} > $i * ($threshold/100) / $pvalNum;
                }
                # Simulation method (previously calculated by findTerms() method)
                else {
                    last if $pval->{'FDR_RATE'} > $threshold/100;
                }
            }
            push @significantPvalues, $pval;

            #----------------------#
            # Adding node to graph #
            #----------------------#
            my $goID = $node->goid();
            my $label = ($goID eq getRootID($aspect))? getAspectName($aspect): &cutLabel($node->term());
            $graph->add_node($goID,
                            label => $label,
                            style => 'filled',
                            fillcolor => &getPvalColor($pval->{$pvalStrg}),
                            URL => $goID );
            $drawnNodes{$goID} = 1; # indicates that this node is now drawn on the graph and has not to be drawn again

            #-------------------------------------------------------------#
            # Write in tab result file (GOID Term pval FDR proteinIDList) #
            #-------------------------------------------------------------#
            print RES $node->goid(),"\t",$node->term(),"\t",$pval->{$pvalStrg},"\t";
            print RES $pval->{'NUM_ANNOTATIONS'},"\t",$pval->{'TOTAL_NUM_ANNOTATIONS'},"\t",;
            my @protIDlist;
            foreach my $proteinName (values(%{$pval->{'ANNOTATED_GENES'}})){
                push @protIDlist, $protIDs{$proteinName};
            }
            print RES join(';',@protIDlist);
            if($criteria eq 'FDR'){
                if($method eq 'BH'){
                    my $correctedPvalue = ($pval->{$pvalStrg} * $pvalNum) / $i;
                    print RES "\t", $correctedPvalue;
                } else {
                    print RES "\t", $pval->{'FDR_RATE'};
                }
            }
            print RES "\n";
        }

        #----------------------------------------------------------------------#
        # Drawing node ancestors (of each node already drawn i.e. significant) #
        #----------------------------------------------------------------------#
        my $showAllNodes = ($drawGraph == 2)? 1:0;
        foreach my $pval (@significantPvalues){
            &drawAncestors($pval->{'NODE'},\$graph,\*RES,\%drawnNodes,$drawGraph);
        }

        #------------------------------------------#
        # Generating image and HTML image map data #
        #------------------------------------------#
        if($drawGraph >= 1){
            $graph->as_png("$promsPath{'tmp'}/GO/$fileName/graph_$aspect.png");
            $graph->as_imap("$promsPath{'tmp'}/GO/$fileName/imap_$aspect.tmp");
        }

        #-----------------------------------#
        # Checking for unannotated proteins #
        #-----------------------------------#
        if(scalar(@unannotatedProteins)){
            print RES "unannotated\t \t \t". scalar(@unannotatedProteins) ."\t \t", join(';',@unannotatedProteins),"\n";
        }

        #---------------------------------------------------------#
        # Checking for discarded proteins (if defined background) #
        # (Proteins not found in provided background)             #
        #---------------------------------------------------------#
        if(scalar(@population) && scalar $termFinder->discardedGenes){
            my @discardedProtIDs;
            foreach my $protName ($termFinder->discardedGenes){
                push @discardedProtIDs, $protIDs{$protName}
            }
            print RES "discarded\t \t \t". scalar(@discardedProtIDs) ."\t \t", join(';',@discardedProtIDs),"\n";
        }
        close(RES);
        subPrint(" OK<BR>\n");
    }
    #print JOB "#! FINISHED\n";
    #close(JOB);
    return 1;
}


sub drawAncestors{
    my ($node,$graphRef,$resFileHandle,$drawnNodes,$drawGraph) = @_;

     # fetching all parents of current node
    foreach my $parent ($node->parentNodes){

        next if $parent->goid() eq &getGoRootID;

        # parent node already drawn -> just add an edge between parent and current node
        if ($drawnNodes->{$parent->goid()}){
            $$graphRef->add_edge($parent->goid() => $node->goid())
        }
        # or not already drawn (= unsignificant nodes)
        else {
            # add a node with term name and grey background if user specified to show all nodes in graph
            if($drawGraph == 2){
                $$graphRef->add_node($parent->goid(),
                                     URL => $parent->goid(),
                                     label => &cutLabel($parent->term()),
                                     style => 'filled',
                                     fillcolor => '0, 0, 0.95'
                                     );
            }
            # else just add a small empty circle
            else {
                $$graphRef->add_node($parent->goid(),
                                     URL => $parent->goid(),
                                     label => '' ,
                                     shape => 'circle',
                                     height => 0.14,
                                     width => 0.14
                                     );
                print $resFileHandle $parent->goid(),"\t",$parent->term(),"\t",'2',"\t",'.',"\t",'.',"\t",'.',"\n"; #pval=2 means unsignificant but parent of a significant node, this info will be used in image map information fetching in displayGOAnalysis.cgi
            }
            # put that node drawn
            $drawnNodes->{$parent->goid()} = 1;
            # add an edge between this parent and its child
            $$graphRef->add_edge($parent->goid() => $node->goid());
            # re-applying function to this parent
            &drawAncestors($parent,$graphRef,$resFileHandle,$drawnNodes,$drawGraph);
        }
    }
}

sub getPvalColor{
    # HSV gradient (Hue,Saturation,Lightness)
    # Be sure that this function fits with the function in displayGOAnalysis.cgi
    my $pval = $_[0];
    $pval=1e-307 if $pval < 1e-307;
    my $logPval = -(log $pval)/(log 10);
    if($logPval >= 9){
        return "0, 1, 1";
    } elsif($logPval >= 5){
        my $hue = ($logPval*(-0.0375))+0.3375;
        return "$hue, 1, 1";
    } elsif($logPval >= 2) {
        my $saturation = ($logPval*(1/3))-(2/3);
        return "0.15, $saturation, 1";
    } else {
        return "white";
    }
}

sub cutLabel{
    if(length($_[0]) > 20){
        my @label = split(//,$_[0]);
        my $count = 0;
        foreach my $pos (0..scalar(@label)-1){
            $count++;
            if($count >= 20 && $label[$pos] eq ' '){
                $label[$pos] = "\n";
                $count = 0;
            }
        }
        return join '', @label;
    } else {
        return $_[0];
    }
}

# OBSOLETE
#sub giToUniprot{
#    # Function performing ID mapping at EBI, returns a hash table GI -> Uniprot_ID #
#    my @giList = @{$_[0]};
#    my %promsPath=&promsConfig::getServerInfo('no_user');
#
#    my %giToUniprot;
#    foreach my $gi (@giList){
#        $giToUniprot{$gi} = undef; # at the end, undef values will indicate unmapped gi identifiers
#    }
#
#    my $giString;
#    # Splitting request for each 50 proteins #
#    my $try = 0;
#    SPLICE:while(my @subGiList = splice(@giList, 0, 50)){
#        $try++;
#        $giString = join ' ',@subGiList;
#        my $params = {
#                from => 'P_GI',
#                to => 'ID',
#                format => 'tab',
#                query => $giString
#        };
#        my $agent = LWP::UserAgent->new;
#        $agent->timeout(360);
#        if (!$promsPath{'http_proxy'}) {$agent->env_proxy;}
#        elsif (lc($promsPath{'http_proxy'}) eq 'no') {$agent->no_proxy;}
#        else {$agent->proxy('http', $promsPath{'http_proxy'});}
#        push @{$agent->requests_redirectable}, 'POST';
#        my $response = $agent->post("http://www.uniprot.org/mapping/", $params);
#        while (my $wait = $response->header('Retry-After')) {
#                print STDERR "Waiting ($wait)...\n";
#                sleep $wait;
#                $response = $agent->get($response->base);
#        }
#        if ($response->is_success){
#            my @responseContent = split("\n",$response->content);
#            shift @responseContent;
#            foreach my $mapping (@responseContent){
#                my @IDs = split "\t", $mapping;
#                if($#IDs != 1){
#                    if($try <= 5){
#                        redo SPLICE;
#                    } else {
#                        &processFailure("Unexpected response content from Uniprot mapping service");
#                    }
#                }
#                #$giToUniprot{$IDs[0]} = $IDs[1] unless $giToUniprot{$IDs[0]};
#                push @{$giToUniprot{$IDs[0]}}, $IDs[1];
#            }
#            $try = 0;
#        } else {
#            if($try <= 5){
#                redo SPLICE;
#            } else {
#                &processFailure("got " . $response->status_line . " for " . $response->request->uri);
#            }
#        }
#        #print JOB "# NULL\n"; # just to alert monitoring that the process is still alive
#    }
#    return %giToUniprot;
#}

# OBSOLETE
#sub processFailure{
#    # Function used to process any failure of this script in order that startGOAnalysis.cgi monitoring gets this failure
#    # 1st argument should be the error message ($!)
#    #$dbh->disconnect;
#    if($_[0]){
#        print JOB "#! FAILED:<FONT color=red>Processing failed: $_[0]</FONT>";
#        close(JOB);
#        die $_[0];
#    } else {
#        print JOB "#! FAILED:<FONT color=red>Processing failed</FONT>";
#        close(JOB);
#        exit;
#    }
#}

sub correctAnnotationForSlim{
    my ($dbh,$annotation,$refHashSlimOntologies,$refAspects) = @_;

    my $completeOboID = $dbh->selectrow_array("SELECT ID_ONTOLOGY FROM ONTOLOGY WHERE STATUS=2 ORDER BY VERSION_DATE DESC LIMIT 0,1")
    or die 'No complete stored ontology in database as a reference in cases using slim ontologies';

    my %promsPath = &promsConfig::getServerInfo;
    my $completeOboFile = "$promsPath{obo}/$completeOboID.obo";

    # Storing complete ontologies to get path to valid terms in slim #
    my %completeOntology;
    foreach my $aspect (@$refAspects){
        $completeOntology{$aspect} = promsOboParser->new(
                                            ontologyFile => $completeOboFile ,
                                            aspect => $aspect);
    }

    # Browsing annotation structure #
    my $waitCount=0;
    while(my ($databaseID,$kAspect) = each %{$annotation->{'GO::AnnotationProvider::AnnotationParser::__goids'}}){
        while(my ($aspect,$kgoID) = each %{$kAspect}){
            next unless grep {$aspect eq $_ } @$refAspects; # skip unwanted aspect
            foreach my $goID (keys %{$kgoID}){
                next if $refHashSlimOntologies->{$aspect}->nodeFromId($goID); # skip if term is already found in slim
                next unless $completeOntology{$aspect}->nodeFromId($goID); # skip if term is not found in complete ontology
                my @validNodes;
                &getValidNodes($completeOntology{$aspect}->nodeFromId($goID),$refHashSlimOntologies->{$aspect},\@validNodes);
                # Updating annotation structure #
                foreach my $validNode (@validNodes){
                    $annotation->__storeGOID($databaseID,$validNode->goid,$aspect);
                }
                delete $annotation->{'GO::AnnotationProvider::AnnotationParser::__goids'}{$databaseID}{$aspect}{$goID};
                $waitCount++;
                if ($waitCount==5000) {
                    $waitCount=0;
                    &subPrint('.');
                }
            }
        }
    }
}

sub getValidNodes{
    my ($node,$slimOntology,$refNodeList,$refHashTerms) = @_;

    unless(defined $refHashTerms){ # this hash reference contains already mapped nodes and their parents in order to not map them again
        %{$refHashTerms} = ();
    }

    if($slimOntology->nodeFromId($node->goid)){# && !$refHashTerms->{$node->goid}){
        push @$refNodeList, $node;
        $refHashTerms->{$node->goid} = 1;
        foreach my $ancestor ($node->ancestors){
            $refHashTerms->{$ancestor->goid} = 1;
        }
    } else {
        my @parentNodes = $node->parentNodes;
        foreach my $parentNode (@parentNodes){
            &getValidNodes($parentNode,$slimOntology,$refNodeList,$refHashTerms);
        }
    }
}

sub checkOntology{
    my $ontology = shift;
    if(!$ontology->rootNode){
        my $aspect = getAspectName($ontology->__aspect);
        subPrint("<FONT color=\"red\">Root node in $aspect not found, it will not be processed.</FONT>\n");
        return 0;
    } else {
        return 1;
    }
}

sub selectID{
    # Return the first identifier from the input list that is annotated by the annotation provider
    my ($inputIdsRef, $annotation) = @_;

    my $chosenUniprot = $inputIdsRef->[0];

    if(scalar @{$inputIdsRef} > 1){
        foreach my $ID (@{$inputIdsRef}){
            if($annotation->nameIsAnnotated(name => $ID)){
                $chosenUniprot = $ID;
                last;
            }
        }
    }
    return $chosenUniprot;
}

sub getAspectName{
    my $aspectChar = shift;

    my %aspectName = ( 'P' => 'Biological Process',
                      'C' => 'Cellular Component',
                      'F' => 'Molecular Function');

    return $aspectName{$aspectChar} or die 'Invalid aspect';
}

sub getRootID{
    my $aspectChar = shift;

    my %rootID = ( 'P' => 'GO:0008150',
              'C' => 'GO:0005575',
              'F' => 'GO:0003674');

    return $rootID{$aspectChar};
}

sub getGoRootID{
    return 'GO:0003673'; # "Gene_Ontology" term
}

sub subPrint{

    if($verbose){
        while(my $text = shift){
            print $text; #,"\n"
        }
    }
    else {
        print '<!-- ',@_,' -->';
    }

}

sub getPval{
    my ($fileName, $thisGOID, $aspect) = @_;

    my %promsPath = &promsConfig::getServerInfo;
    open(RES,"$promsPath{'tmp'}/GO/$fileName/results_$aspect.txt");

    my $thisPval;
    while(<RES>){
        if(/^### (\d+)\t(\d+)/){ # 1st line in file
		next;
	}
        my ($goID,$goName,$pval,$numProt,$totNumProt,$protListStrg) = split /\t/;
        if($goID eq $thisGOID){
            $thisPval = $pval;
        }
    }
    close RES;

    return $thisPval;
}

sub storeData{
    my ($projectID, $goAnaID, $tmpName) = @_;

    my %promsPath = &promsConfig::getServerInfo;

    unlink "$promsPath{tmp}/GO/$tmpName/job.txt";
    unlink "$promsPath{tmp}/GO/$tmpName/stdout.txt";
    mkdir $promsPath{'go_unix'} unless -d $promsPath{'go_unix'};
    mkdir "$promsPath{go_unix}/project_$projectID" unless -d "$promsPath{go_unix}/project_$projectID";
    mkdir "$promsPath{go_unix}/project_$projectID/$goAnaID";
    foreach(glob("$promsPath{tmp}/GO/$tmpName/*")){
        move($_,"$promsPath{go_unix}/project_$projectID/$goAnaID/") || print STDERR "$_: $!";
    }
    rmtree("$promsPath{tmp}/GO/$tmpName") || print STDERR $!;
}

sub deleteGOAnalysis{
    # Removing data files
    #foreach my $dataFile (glob("$resultDirUnix/*")){
    #    unlink $dataFile or die "$dataFile : $!";
    #}
    #rmdir $resultDirUnix or die $!;

    my($dbh, $projectID, $goAnaID) = @_;

    # Firstly checking for go analysis children and delete them
    my $sthChildren = $dbh->prepare("SELECT ID_GOANALYSIS FROM GO_ANALYSIS WHERE ID_PARENT_GOANA=?");
    $sthChildren->execute($goAnaID);
    while(my ($childID) = $sthChildren->fetchrow_array){
        deleteGOAnalysis($dbh, $projectID, $childID);
    }
    $sthChildren->finish;

    my %promsPath=&promsConfig::getServerInfo();

    #remove_tree("$promsPath{go_unix}/project_$projectID/$goAnaID");
    rmtree("$promsPath{go_unix}/project_$projectID/$goAnaID");

    my $sthDelGOAna = $dbh->prepare("DELETE FROM GOANA_ANALYSIS WHERE ID_GOANALYSIS=?");
    my $sthDelGOQuanti = $dbh->prepare('DELETE FROM GOANA_QUANTIFICATION WHERE ID_GOANALYSIS=?');
    my $sthDelGO = $dbh->prepare("DELETE FROM GO_ANALYSIS WHERE ID_GOANALYSIS=?");
    $sthDelGOAna->execute($goAnaID) or die $!;
    $sthDelGOQuanti->execute($goAnaID) or die $!;
    $sthDelGO->execute($goAnaID) or die $!;
    $sthDelGOAna->finish;
    $sthDelGOQuanti->finish;
    $sthDelGO->finish;

}

sub formatGOName{
	my $goName = shift;

	# applying uppercase to term first letter, unless this letter must be lowercase (eg: mRNA, ncRNA...)
	unless($goName =~ /^\w{1,3}RNA/){
		$goName = ucfirst($goName);
	}

	return $goName;
}

#sub getProteinUniprotAc{
#    my ($dbh,$protID) = @_;
#
#    my @uniprotAcList;
#
#    my $sthProtId = $dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER MI, PROTEIN P
#                                  WHERE MI.ID_MASTER_PROTEIN=P.ID_MASTER_PROTEIN
#                                  AND P.ID_PROTEIN=?
#                                  AND ID_IDENTIFIER=(SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE=?)");
#    $sthProtId->execute($protID, "ID");
#    while(my ($uniprotAc) = $sthProtId->fetchrow_array){
#        push @uniprotAcList, $uniprotAc;
#    }
#    $sthProtId->finish;
#
#    return @uniprotAcList;
#}

sub getProteinIds{
    my ($dbh, $protID, $identifierID) = @_;

    my @idList;
    my $sthProtId = $dbh->prepare("SELECT VALUE FROM MASTERPROT_IDENTIFIER MI, PROTEIN P
                                  WHERE MI.ID_MASTER_PROTEIN=P.ID_MASTER_PROTEIN
                                  AND P.ID_PROTEIN=?
                                  AND ID_IDENTIFIER=?");
    $sthProtId->execute($protID, $identifierID);
    while(my ($id) = $sthProtId->fetchrow_array){
        push @idList, $id;
    }
    $sthProtId->finish;

    return @idList;
}

1;

####>Revision history<####
# 1.1.5 More print to keep server connection (PP 20/08/15)
# 1.1.4 Minor code cleaning (PP 12/11/14)
# 1.1.3 Check for p-value < 1e-307 before log (PP 01/07/14)
# 1.1.2 Uses rmtree instead of remove_tree (PP 10/03/14)
# 1.1.1 Change syntax $x=$y // $z to old ()? one (PP 10/03/14)
# 1.1.0 Catching 'die' from TermFinder process and display explicit error message (FY 24/10/13)
# 1.0.9 Better proxy declaration for LWP (PP 02/07/13)
# 1.0.8 New method to get protein identifiers according to the annotation format (FY 22/03/13)
# 1.0.7 Fixes revision history (PP 13/12/12)
# 1.0.6 Switching OboParser use to promsOboParser<BR>& Uniprot identifiers get by convertion are now selected based on their annotation quality in GOA file<BR>& Using Getopt module to parse script arguments<BR>& DIE signal redirection to print all fatal errors in job file (FY 19/07/12)
# 1.0.5 Corrected annotation for slim ontologies now stores only the most specific terms<BR>& Analysing uniprot mapping service response contents<BR>& Filtering category dataset with proteins only in current experiment (FY 23/04/12)
# 1.0.4 Redo (up to 5 times) if uniprot mapping fails<BR>& Skip terms in slim annotation correction if they are not found in complete ontology<BR>& Hide "Gene_Ontology" term appearing in graphs whene there is no threshold (FY 18/04/12)
# 1.0.3 Bug fix for backgroud population (FY 16/04/12)
# 1.0.2 Decrease size and increase time out of requests to uniprot mapping service (FY 13/04/12)
# 1.0.1 Add annotation correction when using slim ontologies (FY 03/04/12)
# 1.0.0 New script for performing GO term enrichment analysis with data provided by goAnalysis.cgi (FY 21/06/11)