#!/usr/local/bin/perl -w

################################################################################
# promsImportDataManager.pm         1.0.18                                     #
# Authors: M. Le Picard, V. Sabatet (Institut Curie)                           #
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

package promsImportDataManager;
require Exporter;

@ISA=qw(Exporter);
@EXPORT=qw();
@EXPORT_OK=qw();
$VERSION=1.00;

use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use strict;
use IO::Handle;
use promsConfig;
use promsMod;
use File::Copy;
use File::Basename;
use XML::Simple;
use POSIX qw(strftime); # to get the time
use File::Path qw(rmtree); # remove_tree
use String::Util qw(trim);
use File::Spec::Functions qw(splitpath); # Core module

$|=1; # buffer flush (allows instant printing)

# Project proteins information
my (%projectProt, %proteinLength, %proteinDescription, %incompleteProtein, %noDescriptionProtein, %noSequenceProtein, %noMWProtein);

# Protein information
my (%protLength, %protDes, %protSeq, %protMW, %protOrg, %numMatchProt);

sub new {
    my ($class, $userID, $projectID, $experimentID, $format, $quantiName, $quantiAnnot, $taxonomy, $instrument, $dataFile, $logFile, $dbIDRef) = @_;
    my $self = {
        _userID        => $userID,
        _projectID     => $projectID,
        _experimentID  => $experimentID,
        _format        => $format,
        _quantiName    => $quantiName,
        _quantiAnnot   => $quantiAnnot,
        _taxonomy      => $taxonomy,
        _instrument    => $instrument,
        _dataFile      => $dataFile,
        _logFile       => $logFile,
        _dbID          => $dbIDRef,
    };
    
    ###> Fetching project information
    my %promsPath=&promsConfig::getServerInfo('no_user');
    my $dbh = &promsConfig::dbConnect('no_user');
    (my $protVisibility, $self->{_identProject}) = $dbh->selectrow_array("SELECT PROT_VISIBILITY, ID_IDENTIFIER FROM PROJECT WHERE ID_PROJECT=".$self->{_projectID});
    $dbh->disconnect;
    
    $self->{_protVisibility} = ($protVisibility) ? $protVisibility : 0;
    $self->{analysisID} = {};
    $self->{peptidesID} = {};
    $self->{quantiID} = 0;
    
    bless $self, $class;
    return $self;
}

sub setPeptidesInformation {
    my ($self, $peptideInfoRef, $peptideProtRef, $peptideModifRef, $peptidePositionRef, $peptideListRef, $fragmentsInfosRef, $modifIDsRef, $assosProtPeptideRef, $dbRankProtRef) = @_;
    $self->{_peptideInfo} = $peptideInfoRef;
    $self->{_peptideProt} = $peptideProtRef;
    $self->{_peptideModif} = $peptideModifRef;
    $self->{_peptidePosition} = $peptidePositionRef;
    $self->{_peptideList} = $peptideListRef;
    $self->{_fragmentsInfos} = $fragmentsInfosRef;
    $self->{_modifIDs} = $modifIDsRef;
    $self->{_assosProtPeptide} = $assosProtPeptideRef;
    $self->{_dbRankProt} = $dbRankProtRef;
}

sub setDIAInfos {
    my ($self, $libIDs, $decoyNbRef, $targetNbRef, $fdr) = @_;
    $self->{_libIDs} = $libIDs;
    $self->{_decoyNb} = ($decoyNbRef) ? $decoyNbRef : '';
    $self->{_targetNb} = ($targetNbRef) ? $targetNbRef : '';
    $self->{_fdr} = ($fdr) ? $fdr : '';
}

sub setFilesInfos {
    my ($self, $analysisFileNamesRef, $sampleFileNamesRef) = @_;
    $self->{_analysisFileNames} = $analysisFileNamesRef;
    $self->{_sampleFileNames} = $sampleFileNamesRef;
}

sub setProtCoverage {
    my ($self, $protCoverageRef) = @_;
    $self->{_protCoverage} = $protCoverageRef;
}

sub setDesignExperiment {
    my ($self, $designExperiment) = @_;
    $self->{designExperiment} = $designExperiment;
}

sub getAnalysisID {
    my ($self) = @_;
    return $self->{analysisID};
}

sub getAnaObsID {
    my ($self) = @_;
    return $self->{anaObsID};
}

sub getPeptidesID {
    my ($self) = @_;
    return $self->{peptidesID};
}

sub getProteinsID {
    my ($self) = @_;
    return $self->{proteinsID};
}

sub getQuantiID {
    my ($self) = @_;
    return $self->{quantiID};
}

sub importData {
    my ($self, $peptideInfoRef, $peptideProtRef, $peptideModifRef, $peptidePositionRef, $peptideListRef, $fragmentsInfosRef, $modifIDsRef, $assosProtPeptideRef, $dbRankProtRef) = @_;

    $self->{_peptideInfo} = $peptideInfoRef;
    $self->{_peptideProt} = $peptideProtRef;
    $self->{_peptideModif} = $peptideModifRef;
    $self->{_peptidePosition} = $peptidePositionRef;
    $self->{_peptideList} = $peptideListRef;
    $self->{_fragmentsInfos} = $fragmentsInfosRef;
    $self->{_modifIDs} = $modifIDsRef;
    $self->{_assosProtPeptide} = $assosProtPeptideRef;
    $self->{_dbRankProt} = $dbRankProtRef;
    
    $self->_fetchProjectProteins;
    
    ##########################
    ### RT-ordered peptide ###
    ##########################
    my (%peptideRTOrder, %peptideRT);
    foreach my $analysis (keys %{$self->{_peptideList}}) {
        foreach my $peptide (keys %{$self->{_peptideList}{$analysis}}) {
            foreach my $rt (sort {$self->{_peptideList}{$analysis}{$peptide}{$a} <=> $self->{_peptideList}{$analysis}{$peptide}{$b}} keys %{$self->{_peptideList}{$analysis}{$peptide}}) {
                $peptideRT{$peptide} = $rt;
            }
        }
        
        my $RTorder=1;
        foreach my $peptide (sort {$peptideRT{$a} <=> $peptideRT{$b}} keys %peptideRT) {
            @{$peptideRTOrder{$analysis}{$peptide}}=($peptideRT{$peptide}, $RTorder);
            $RTorder++;
        }
    }
    
    ##########################################################################################
    ### Inserting data into tables SAMPLE,ANALYSIS,ANALYSIS_DATABANK,ANALYSIS_MODIFICATION ###
    ##########################################################################################
    my (%proteinID, %proteinList, %anaObsModifs, %anaSampleID, %displayPos);
    my %promsPath=&promsConfig::getServerInfo('no_user');
    my $dbh = &promsConfig::dbConnect('no_user');

    my $sthSample=$dbh->prepare("INSERT INTO SAMPLE (ID_SAMPLE,ID_EXPERIMENT,NAME,START_DATE,DISPLAY_POS,UPDATE_DATE,UPDATE_USER) VALUES (?,?,?,NOW(),?,NOW(),?)");
    my $sthAnalysis=$dbh->prepare("INSERT INTO ANALYSIS (ID_SAMPLE,NAME,START_DATE,DATA_FILE,VALID_STATUS,VALID_USER,LOWER_SCORES,FIRST_VALID_DATE,VALID_DATE,UPDATE_DATE,UPDATE_USER,WIFF_FILE,DECOY,FALSE_DISCOVERY,TAXONOMY,DISPLAY_POS,FILE_FORMAT,MS_TYPE,INSTRUMENT,VERIFIED_MG,MIN_SCORE,MAX_RANK) VALUES (?,?,NOW(),?,?,?,?,NOW(),NOW(),NOW(),?,?,?,?,?,?,?,?,?,?,?,?)");
    my $sthObs=$dbh->prepare("INSERT INTO OBSERVATION (ID_BIOSAMPLE, ID_ANALYSIS, TARGET_POS)  VALUES (NULL, ?, ?)");
    # my $sthAnaObsModif=$dbh->prepare("INSERT INTO OBS_MODIFICATION (ID_OBSERVATION,ID_MODIFICATION) VALUES (?,?)"); # TODO Add this insertion for PTMs Quantis
    my $sthAnaSwath=$dbh->prepare("INSERT INTO ANALYSIS_SWATH_LIB (ID_ANALYSIS,ID_SWATH_LIB,VERSION_NAME) VALUES (?,?,?)");
    my $sthAnaDB=$dbh->prepare("INSERT INTO ANALYSIS_DATABANK (ID_DATABANK,ID_ANALYSIS,DB_RANK) VALUES (?,?,?)");
    my $sthAnaMod=$dbh->prepare("INSERT INTO ANALYSIS_MODIFICATION (ID_ANALYSIS,ID_MODIFICATION,MODIF_TYPE,SPECIFICITY) VALUES (?,?,?,?)");
    my $sampleID=$dbh->selectrow_array("SELECT MAX(ID_SAMPLE) FROM SAMPLE");
    my $sampleDisplayPos=$dbh->selectrow_array("SELECT MAX(DISPLAY_POS) FROM SAMPLE WHERE ID_EXPERIMENT=".$self->{_experimentID});
    $sampleDisplayPos = 0 if(!$sampleDisplayPos);
    
    $self->appendLog("Storing analysis and samples data into database.");
    
    # Insertion into SAMPLE
    if($self->{designExperiment}) { # If a specific design should be followed, create it
        foreach my $sample (sort { lc($a) cmp lc($b) } keys %{$self->{designExperiment}{"conditions"}}) {
            $sthSample->execute(++$sampleID, $self->{_experimentID}, $sample, ++$sampleDisplayPos, $self->{_userID});

            foreach my $analysis (sort { $a <=> $b } keys %{$self->{designExperiment}{"conditions"}{$sample}}) {
                my %replicateInfos = %{$self->{designExperiment}{"conditions"}{$sample}{$analysis}};
                $anaSampleID{$replicateInfos{"anaName"}} = $sampleID;
            }
            @{$displayPos{$sampleID}} = ($sampleDisplayPos, 0);
        }
    } else {
        foreach my $analysis (sort { lc($a) cmp lc($b) } keys %{$self->{_analysisFileNames}}) {
            my $sampleName = $self->{_sampleFileNames}{$analysis};
            $sthSample->execute(++$sampleID, $self->{_experimentID}, $sampleName, ++$sampleDisplayPos, $self->{_userID});
            $anaSampleID{$analysis} = $sampleID;
            @{$displayPos{$sampleID}} = ($sampleDisplayPos, 0);
        }
    }
    $sthSample->finish;
    
    foreach my $analysis (sort { lc($a) cmp lc($b) } keys %{$self->{_analysisFileNames}}) {
        my $analysisFile = $self->{_analysisFileNames}{$analysis};
        my $analysisName = ($self->{_format} eq 'pkv') ? fileparse($analysisFile, qr/\.[^.]*/) : $analysis;
        $sampleID = $anaSampleID{$analysis};
        my $anaDisplayPos = ++$displayPos{$sampleID}[1];
        
        my $topPepNbDecoy = ($self->{_decoyNb} && $self->{_decoyNb}{$analysis}) ? $self->{_decoyNb}{$analysis} : 0;
        my $topPepNbTarget = ($self->{_targetNb} && $self->{_targetNb}{$analysis}) ? $self->{_targetNb}{$analysis} : 0;

        my $fileFormat=($self->{_format} eq 'pkv')? 'SWATH.PKV': ($self->{_format} eq 'prm')?  'SKYLINE.SKY' : ($self->{_format} eq 'openswath')? 'OPENSWATH.TSV' : 'SPECTRONAUT.XLS';
        my $msType=($self->{_format} eq 'pkv')? 'SWATH': ($self->{_format} eq 'prm')? 'TDA' : 'DIA';
        my $decoy=($self->{_format} eq 'pkv') ? "INT:SEARCH,FDR=".$self->{_fdr}.":precomputed" : '';

        # Insertion into ANALYSIS
        $sthAnalysis->execute($sampleID, $analysisName, $self->{_dataFile}, 2, $self->{_userID}, 0, $self->{_userID}, $analysisFile, $decoy, "$topPepNbDecoy:0:$topPepNbTarget:0", $self->{_taxonomy}, $anaDisplayPos, $fileFormat, $msType, $self->{_instrument}, 0, 0, 1);
        $sthAnalysis->finish;
        my $analysisID=$dbh->last_insert_id(undef, undef, 'ANALYSIS', 'ID_ANALYSIS');
        
        # Insertion into OBSERVATION
        $sthObs->execute($analysisID, 0);
        my $obsID=$dbh->last_insert_id(undef, undef, 'OBSERVATION', 'ID_OBSERVATION');

        $self->{analysisID}{$analysis} = $analysisID;
        $self->{anaObsID}{$analysis} = $obsID;
        
        ###Insertion into ANALYSIS_SWATH_LIB
        if($self->{_format} ne 'prm' && $self->{_libIDs}) {
            my @libID = split(/,/, $self->{_libIDs});
            foreach my $libID (@libID) {
                my ($versionName) = $dbh->selectrow_array("SELECT VERSION_NAME FROM SWATH_LIB WHERE ID_SWATH_LIB='$libID' ") or die "Couldn't prepare statement: ". $dbh->errstr;
                $sthAnaSwath->execute($analysisID, $libID, $versionName);		
            }
            $sthAnaSwath->finish;
        }
        
        # Insertion into ANALYSIS_DATABANK
        foreach my $dbInfo (@{$self->{_dbID}}) {
            my ($dbID,$dbRank) = split(/&/,$dbInfo);
            $sthAnaDB->execute($dbID, $analysisID, $dbRank);
        }

        foreach my $modifID (keys %{$self->{_modifIDs}}) {
            foreach my $modType (keys %{$self->{_modifIDs}{$modifID}}) {
                my @aaList = keys %{$self->{_modifIDs}{$modifID}{$modType}};
                my $aaSpecificity = join(",", @aaList);
                $sthAnaMod->execute($analysisID, $modifID, $modType, $aaSpecificity);
                $sthAnaMod->finish;
            }
            
            #$sthAnaObsModif->execute($obsID, $modifID); # TODO Add this in quantification sites
            #$sthAnaObsModif->finish;
        }

        foreach my $protein (keys %{$self->{_assosProtPeptide}{$analysis}}) {
            $proteinList{$protein} = 1;
        }
    }
    $dbh->commit;

    ## Recovering protein informations
    if($self->{_format} =~ /prm|openswath|spectronaut/) {
        $self->appendLog("Recovering protein informations.");
        $self->_fetchProteinFromTDAFile();
    }
    
    ###########################################
    #####Inserting data into table PROTEIN#####
    ###########################################
    $self->appendLog("Storing proteins data into database.");
    
    # Need %proteinList, %proteinID, %projectProt, %protLength, %protSeq, %protDes, %protMW, %protOrg, %noMWProtein, %noDescriptionProtein, %incompleteProtein, %noSequenceProtein)
    my $nbProtProcessed = 0;
    my $sthProtein=$dbh->prepare("INSERT INTO PROTEIN (ID_PROJECT,ID_MASTER_PROTEIN,ALIAS,IDENTIFIER,PROT_DES,PROT_LENGTH,PROT_SEQ,MW,COMMENTS,UPDATE_DATE,UPDATE_USER,ORGANISM) VALUES (?,NULL,?,?,?,?,?,?,NULL,NULL,NULL,?)");
    foreach my $protein (keys %proteinList) {
        ## recovering the alias of the protein in terms of the project identifier
        my $alias;
        if (!$self->{_identProject}) {
            $alias = $protein;
        } else {
            if ($protein=~/[sp|tr]\|(\w+)\|(\w+_\w+)/) { #uniprot ALL
                if ($self->{_identProject} == 1) {
                    $alias = $1;
                } elsif ($self->{_identProject} == 2) {
                    $alias = $2;
                }
            } else {
                $alias = $protein;
            }
        }

        ### if this protein exists in PROTEIN table
        if(exists $projectProt{$protein}) {
            my $sthProteinUpdate;
            $proteinID{$protein} = $projectProt{$protein};
            $self->{proteinsID}{$protein} = $projectProt{$protein};
            
            # protein length = null (PROTEIN table)
            if (exists $incompleteProtein{$protein}) {
                $sthProteinUpdate = $dbh->prepare("UPDATE PROTEIN SET PROT_LENGTH=? WHERE ID_PROTEIN=?");
                $sthProteinUpdate->execute($protLength{$protein}, $projectProt{$protein});
            }
            
            # no protein description (PROTEIN table)
            if (exists $noDescriptionProtein{$protein}) {
                $sthProteinUpdate = $dbh->prepare("UPDATE PROTEIN SET PROT_DES=? WHERE ID_PROTEIN=?");
                $sthProteinUpdate->execute($protDes{$protein}, $projectProt{$protein});
            }
            
            # no protein sequence (PROTEIN table)
            if (exists $noSequenceProtein{$protein}) {
                $sthProteinUpdate = $dbh->prepare("UPDATE PROTEIN SET PROT_SEQ=? WHERE ID_PROTEIN=?");
                $sthProteinUpdate->execute($protSeq{$protein}, $projectProt{$protein});
            }
            
            # no protein MW (PROTEIN table)
            if (exists $noMWProtein{$protein}) {
                $sthProteinUpdate = $dbh->prepare("UPDATE PROTEIN SET MW=? WHERE ID_PROTEIN=?");
                $sthProteinUpdate->execute($protMW{$protein}, $projectProt{$protein});
            }
        } else { # no duplicates !!
            my $species = ($protOrg{$protein}) ? $protOrg{$protein} : undef; #$self->{_taxonomy};
            $sthProtein->execute($self->{_projectID}, $alias, $protein, $protDes{$protein}, $protLength{$protein}, $protSeq{$protein}, $protMW{$protein}, $species) or die $dbh->errstr;
            my $protID = $dbh->last_insert_id(undef, undef, 'PROTEIN', 'ID_PROTEIN');
            $proteinID{$protein} = $protID;
            $self->{proteinsID}{$protein} = $protID;
        }
        
        $nbProtProcessed++;
        if($nbProtProcessed % 5000 == 0) {
            $self->appendLog("Storing proteins data into database ($nbProtProcessed processed)");
            $dbh->commit;
        }
    }
    $sthProtein->finish;
    $dbh->commit;



    ###################################################
    #### Inserting data into table ANALYSIS_PROTEIN ###
    ###################################################
    my $sthTopProt = $dbh->prepare("SELECT COUNT(*) FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=? AND VISIBILITY=2");
    my $sthBestVis = $dbh->prepare("SELECT MAX(VISIBILITY) FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=?");
    my $sthAnaProt = $dbh->prepare("INSERT INTO ANALYSIS_PROTEIN (ID_ANALYSIS,ID_PROTEIN,DB_RANK,CONF_LEVEL,SCORE,NUM_PEP,NUM_MATCH,PEP_COVERAGE,MATCH_GROUP,PEP_SPECIFICITY,VISIBILITY) VALUES (?,?,?,?,?,?,?,?,?,?,?)");
    
    my $nb = scalar keys %{$self->{analysisID}};
    my $currentAnaIndex = 1;
    
    ## Fetching protein score, pep number and protein specificity
    foreach my $analysis (sort { lc($a) cmp lc($b) } keys %{$self->{analysisID}}) {
        my (%matchList, %proteinScore, %bestPepSpecificity, @pepInfo);
        my $analysisID = $self->{analysisID}{$analysis};
        my $protScore;
        
        foreach my $protein (keys %{$self->{_assosProtPeptide}{$analysis}}) {
            my $pepSpecificity = 0;
            $protScore = 0;
            
            foreach my $peptide (keys %{$self->{_assosProtPeptide}{$analysis}{$protein}}) {
                $matchList{$protein}{$peptide} = 1;
                
                ### Protein score
                $protScore += $self->{_assosProtPeptide}{$analysis}{$protein}{$peptide}[1] if $self->{_assosProtPeptide}{$analysis}{$protein}{$peptide}[1];

                ### Max peptide specificity
                $pepSpecificity = ($self->{_peptideInfo}{$peptide}[2] && $pepSpecificity<$self->{_peptideInfo}{$peptide}[2]) ? $self->{_peptideInfo}{$peptide}[2] : $pepSpecificity;
            }
            $bestPepSpecificity{$protein} = $pepSpecificity;
            $proteinScore{$protein} = $protScore if $protScore;
        }

        
        ### Finding number of times proteins are at top of match group hierarchy
        $self->appendLog("Analysis $analysis : Computing proteins visibility ($currentAnaIndex/$nb)");
        my (%numProtTop, %bestProtVis, %proteinPepDens);
        
        foreach my $protein (keys %matchList) {
            my $protID = $proteinID{$protein};
            
            if ($projectProt{$protein}) { # proteins are already in project
                $sthTopProt->execute($protID);
                ($numProtTop{$protein}) = $sthTopProt->fetchrow_array;
                if ($self->{_protVisibility}) {
                    $sthBestVis->execute($protID);
                    $bestProtVis{$protein} = $sthBestVis->fetchrow_array;
                }
            } else { #### Proteins are new in project
                $numProtTop{$protein} = 0;
                $projectProt{$protein} = $protID; # add protein to project list
            }
            
            ### Peptide density
            $proteinPepDens{$protein} = $protLength{$protein}/(scalar (keys %{$self->{_assosProtPeptide}{$analysis}{$protein}}));
        }
        $sthTopProt->finish;
        $sthBestVis->finish;


        ### Finding match groups and protein visibility
        $self->appendLog("Analysis $analysis : Building match groups ($currentAnaIndex/$nb)");
        my (%visibility, %matchGroup, @sortedProt);
        if(%proteinScore) {
            @sortedProt = sort{scalar (keys %{$self->{_assosProtPeptide}{$analysis}{$b}})<=>scalar (keys %{$self->{_assosProtPeptide}{$analysis}{$a}}) || $numProtTop{$b}<=>$numProtTop{$a} || $proteinScore{$b}<=>$proteinScore{$a} || &deltaLength($protLength{$a},$protLength{$b},$proteinPepDens{$a},$proteinPepDens{$b})  || $protLength{$a}<=>$protLength{$b} } keys %matchList;
        } else {
            @sortedProt = sort{scalar (keys %{$self->{_assosProtPeptide}{$analysis}{$b}})<=>scalar (keys %{$self->{_assosProtPeptide}{$analysis}{$a}}) || $numProtTop{$b}<=>$numProtTop{$a} || &deltaLength($protLength{$a},$protLength{$b},$proteinPepDens{$a},$proteinPepDens{$b})  || $protLength{$a}<=>$protLength{$b} } keys %matchList;
        }
        &promsMod::createMatchGroups(\%matchList, \%matchGroup, \@sortedProt, \%visibility, \%bestProtVis, $self->{_protVisibility}, 0);
        
        ## Inserting data into table ANALYSIS_PROTEIN
        foreach my $protein (keys %{$self->{_assosProtPeptide}{$analysis}}) { # Replaced because when using %proteinID, it inserts, for each analysis, ALL proteins existing in all analysis, even though it does not contain these proteins initially !
            my $dbRank = $self->{_dbRankProt}{$protein};
            my $coverage = $self->{_protCoverage}{$analysis}{$protein};
            my $proteinID = $proteinID{$protein};
            my $protScore = ($proteinScore{$protein}) ? $proteinScore{$protein} : '';
            my $numMatch = $numMatchProt{$analysis}{$protein};
            my $matchGroup = $matchGroup{$protein};
            my $pepSpecificity = $bestPepSpecificity{$protein};
            my $visibility = $visibility{$protein};
            
            my $numPepNotImputed = 0;
            foreach my $peptide (keys %{$self->{_assosProtPeptide}{$analysis}{$protein}}) {
                $numPepNotImputed++ unless($self->{_assosProtPeptide}{$analysis}{$protein}{$peptide}[2]);
            }

            $sthAnaProt->execute($analysisID, $proteinID, $dbRank, 2, $protScore, $numPepNotImputed, $numMatch, $coverage, $matchGroup, $pepSpecificity, $visibility);
        }
        $sthAnaProt->finish;
        $currentAnaIndex++;
    }
    $dbh->commit;

    my $sthPeptide=$dbh->prepare("INSERT INTO PEPTIDE (ID_ANALYSIS,PEP_SEQ,PEP_LENGTH,QUERY_NUM,PEP_RANK,SEARCH_RANK,SCORE,MISS_CUT,MR_EXP,MR_CALC,MR_OBS,MR_DELTA,COMMENTS,SUBST,CHARGE,ELUTION_TIME,VALID_STATUS,DATA,SPEC_COUNT) VALUES (?,?,?,?,?,?,?,?,?,NULL,?,NULL,NULL,NULL,?,?,?,NULL,?)");
    my $sthModification=$dbh->prepare("INSERT INTO PEPTIDE_MODIFICATION (ID_PEPTIDE,ID_MODIFICATION,MODIF_TYPE,POS_STRING,REF_POS_STRING) VALUES (?,?,?,?,?)");
    my $sthAttrib=$dbh->prepare("INSERT INTO PEPTIDE_PROTEIN_ATTRIB (ID_PROTEIN,ID_PEPTIDE,ID_ANALYSIS,PEP_BEG,PEP_END,FLANKING_AA,IS_SPECIFIC) VALUES (?,?,?,?,?,?,?)");
    my $sthUpdatePeptide=$dbh->prepare("UPDATE PEPTIDE SET MR_CALC=?, MR_DELTA=? WHERE ID_PEPTIDE=?");
    
    ###############################################################
    ### Inserting data into tables PEPTIDE,PEPTIDE_MODIFICATION ###
    ###############################################################
    $self->appendLog("Storing peptides data into database.");
    my $nbPepProcessed = 0;
    my $nbPepAssocProcessed=0;

    foreach my $analysis (keys %{$self->{analysisID}}) {
        my $analysisID = $self->{analysisID}{$analysis};
        foreach my $peptide (keys %{$self->{_peptideList}{$analysis}}) {
            my $ghostPep = 0;
            foreach my $rt (sort {$self->{_peptideList}{$analysis}{$peptide}{$a} <=> $self->{_peptideList}{$analysis}{$peptide}{$b}} keys %{$self->{_peptideList}{$analysis}{$peptide}}) {
                my @pepSeq = split(/_/, $peptide);
                my $pepSequence = $pepSeq[0];
                $pepSequence=~s/\(UniMod:\d+\)//g; # OpenSwath
                $pepSequence=~s/(^\.|\.$)//g; # OpenSwath DIA
                $pepSequence=~s/\[([A-Za-z0-9 -]+)\s\(\w+\)\]//g; # TDA
                $pepSequence=~s/\(unimod:\d+\)//g; # TDA
                $pepSequence=~s/\[[+-]?\w+\.?\d*\]//g; # Peakview / Spectronaut
                $pepSequence=~s/\[[^\[]+\]//g; # Spectronaut             
                
                $ghostPep = ($self->{_peptideList}{$analysis}{$peptide}{$rt}[4]) ? 1 : 0; # Specific to spectronaut when using imputation
                my $pepLength = length $pepSequence;
                my $pepCharge = $pepSeq[1];
                my $mrObs = ($self->{_peptideInfo}{$peptide}[0]) ? $self->{_peptideInfo}{$peptide}[0] : $self->{_peptideList}{$analysis}{$peptide}{$rt}[0] ;
                my $mrExp = ($mrObs-1.007825032)*$pepCharge;
                my $missCut = ($self->{_peptideInfo}{$peptide}[1]) ? $self->{_peptideInfo}{$peptide}[1] : '0' ;
                my $pepScore = ($ghostPep) ? undef : ($self->{_peptideList}{$analysis}{$peptide}{$rt}[1] && $self->{_peptideList}{$analysis}{$peptide}{$rt}[1]=~/\d+/) ? $self->{_peptideList}{$analysis}{$peptide}{$rt}[1] : undef;
                
                $rt = sprintf("%0.2f", $rt);
                
                # Update peptide min retention time based on its fragments (to avoid undefined RT values due to #N/A value in data file)
                if($self->{_fragmentsInfos}{$analysis}{$peptide} && $rt eq "0.00") {
                    my $rtMean = 0;
                    my $nbFragments = 0;
                    foreach my $rt (sort {$self->{_fragmentsInfos}{$analysis}{$peptide}{$a} <=> $self->{_fragmentsInfos}{$analysis}{$peptide}{$b}} keys %{$self->{_fragmentsInfos}{$analysis}{$peptide}}) {
                        foreach my $fragment (keys %{$self->{_fragmentsInfos}{$analysis}{$peptide}{$rt}}) {
                            my $fragmentRT = ($self->{_fragmentsInfos}{$analysis}{$peptide}{$rt}{$fragment}[0]) ? $self->{_fragmentsInfos}{$analysis}{$peptide}{$rt}{$fragment}[0] : '0.00';
                            $rtMean += $fragmentRT;
                            $nbFragments++;
                        }
                    }
                    $rt = ($rtMean != 0 && $nbFragments != 0) ? $rtMean/$nbFragments : "0.00";
                }
                
                my $queryNum = ($ghostPep) ? undef : $peptideRTOrder{$analysis}{$peptide}[1];
                my $rank = ($ghostPep) ? undef : 1;
                my $validStatus = ($ghostPep) ? 0 : 1;
                my $elutionTime = ($pepScore) ? "sc".$pepScore.";et".$rt : "sc?;et".$rt;
                
                ###>Insertion into PEPTIDE
                $sthPeptide->execute($analysisID, $pepSequence, $pepLength, $queryNum, $rank, $rank, $pepScore, $missCut, $mrExp, $mrObs, $pepCharge, $elutionTime, $validStatus, 1);
                my $peptideID = $dbh->last_insert_id(undef, undef, 'PEPTIDE', 'ID_PEPTIDE');
                push(@{$self->{peptidesID}{$analysis}{$peptide}{$rt}}, $peptideID);
                
                ###>Insertion into PEPTIDE_MODIFICATION
                my %varMods;
                my $refPosStr; # Contains PTMs scores if any
                foreach my $modificationID (keys %{$self->{_peptideModif}{$analysis}{$peptide}}) {
                    my $position = "";
                    $refPosStr = undef;
                    foreach my $pos (sort {$a <=> $b} keys %{$self->{_peptideModif}{$analysis}{$peptide}{$modificationID}}) { # pos = position of modif OR "SCORE" -> represent PTMs score/probabilities
                        if($pos eq "-1") { # Correspond to PTM Probabilities
                            my $softwarePTMShort = ($self->{_format} eq 'spectronaut') ? "SPC" : ($self->{_format} eq 'openswath') ? 'OS' : '';
                            $refPosStr = "##PRB_$softwarePTMShort=".$self->{_peptideModif}{$analysis}{$peptide}{$modificationID}{$pos}; # <position1>:<proba1>;<position2>:<proba2>...
                        } else {
                            $pos = '=' if ($pos == 0);
                            $position = ($position) ? "$position.$pos" : $pos;
                            $varMods{$pos} = ($pos eq '=') ? $self->{_peptideModif}{$analysis}{$peptide}{$modificationID}{0} : $self->{_peptideModif}{$analysis}{$peptide}{$modificationID}{$pos};
                        }
                    }
                    $sthModification->execute($peptideID,$modificationID,'V',$position,$refPosStr);
                }
                
                my $mrCalc = &promsMod::mrCalc($pepSequence, \%varMods); ## $pepSequence = peptide sequence without modification ; %varMod => {modification position}=modification mass
                my $delta = $mrCalc-$mrExp;
                $sthUpdatePeptide->execute($mrCalc, $delta ,$peptideID);
                $nbPepProcessed++;
                
                if($nbPepProcessed % 10000 == 0) {
                    $self->appendLog("Storing peptides data into database ($nbPepProcessed processed)");
                    $dbh->commit;
                }
            }
        }
        $sthPeptide->finish;
        $sthUpdatePeptide->finish;
        $sthModification->finish;
        $dbh->commit;

        ########################################################
        ####Inserting data into tables PEPTIDE_PROTEIN_ATTRIB###
        ########################################################
        $self->appendLog("Storing peptides proteins association data into database.");
        foreach my $protein (keys %{$self->{_assosProtPeptide}{$analysis}}) {
            my $proteinID = $proteinID{$protein};
            foreach my $peptide (keys %{$self->{_assosProtPeptide}{$analysis}{$protein}}) {
                my $ghostPep = 0;
                foreach my $rt (keys %{$self->{peptidesID}{$analysis}{$peptide}}) {
                    $ghostPep = ($self->{_peptideList}{$analysis}{$peptide}{$rt}[4]) ? 1 : 0; # Specific to spectronaut when using imputation
                    foreach my $peptideID (@{$self->{peptidesID}{$analysis}{$peptide}{$rt}}) {
                        ### Max peptide specificity
                        my $specificity;
                        $specificity = ($self->{_peptideInfo}{$peptide}[2] && $self->{_peptideInfo}{$peptide}[2]==100) ? 1 : 0;
    
                        ### Peptide position on protein
                        foreach my $begPepInfo (keys %{$self->{_peptidePosition}{$analysis}{$protein}{$peptide}}) {
                            my @begPeptideInfo = split(/_/,$begPepInfo);
                            my $pepBeg = ($ghostPep) ? "-$begPeptideInfo[0]" : $begPeptideInfo[0];
                            my $endPepInfo = $self->{_peptidePosition}{$analysis}{$protein}{$peptide}{$begPepInfo};
                            my @endPeptideInfo = split(/_/,$self->{_peptidePosition}{$analysis}{$protein}{$peptide}{$begPepInfo});
                            my $pepEnd = ($ghostPep) ? "-$endPeptideInfo[0]" : $endPeptideInfo[0];
                            my $flanking = "$begPeptideInfo[1]$endPeptideInfo[1]";
                            $sthAttrib->execute($proteinID,$peptideID,$analysisID,$pepBeg,$pepEnd,$flanking,$specificity);
                        }
                        $nbPepAssocProcessed++;

                        if($nbPepAssocProcessed % 10000 == 0) {
                            $self->appendLog("Storing peptides proteins association data into database ($nbPepAssocProcessed processed)");
                            $dbh->commit;
                        }
                    }
                }
            }
        }
        $sthAttrib->finish;
        $dbh->commit;
    }

    ## Get all analysisID
    my @analysisIDs = ();
    foreach my $analysis (sort { lc($a) cmp lc($b) } keys %{$self->{analysisID}}) {
        push(@analysisIDs, $self->{analysisID}{$analysis});
    }
    
    #######################################################################
    ####Inserting data into tables QUANTIFICATION and ANA_QUANTIFICATION###
    #######################################################################
    if($self->{_quantiName} && $self->{_quantiAnnot}) {
        my $code = ($self->{_format} eq 'prm') ? 'XIC' : 'DIA';
        my $quantiMethod = $dbh->selectrow_array("SELECT ID_QUANTIFICATION_METHOD FROM QUANTIFICATION_METHOD WHERE CODE='$code'");
        
        $self->{quantiID} = $self->createQuantification($self->{_quantiName}, 'peptide', $self->{_quantiAnnot}, $quantiMethod, \@analysisIDs);
    }
    
    ############################
    ## UPDATE ANALYSIS STATUS ##
    ############################
    my $sthUpdateAnalysis = $dbh->prepare("UPDATE ANALYSIS SET VALID_STATUS=2 WHERE ID_ANALYSIS=?");
    foreach my $analysis (sort { lc($a) cmp lc($b) } keys %{$self->{analysisID}}) {
        my $analysisID = $self->{analysisID}{$analysis};
        $sthUpdateAnalysis->execute($analysisID);
    }
    
    ######################
    ## MAP NEW PROTEINS ##
    ######################
    system "./mapProteinIdentifiers.pl ".$self->{_userID}." ".join(',', @analysisIDs);
        
    $dbh->commit;
    $dbh->disconnect;
}

sub appendLog {
    my ($self, $content) = @_;
    
    if($self->{_logFile} && $content) {
        open(FILE,">>".$self->{_logFile}) || die "Error while opening ".$self->{_logFile}."\n";
        print FILE $content."\n";
        close FILE;
    }
}

##########################
####Check delta length#### Compares delta lengthes between 2 proteins with Peptide density of the smaller one
########################## (delta is considered not significant if smaller than pep density)
sub deltaLength {
	my ($l_a,$l_b,$d_a,$d_b)=@_;
	my $pepDensVal=($l_b > $l_a)? $d_a : $d_b;
	if (abs($l_b-$l_a) > $pepDensVal) {return $l_a<=>$l_b} else {return 0}
}

sub _fetchProjectProteins {
    my ($self) = @_;
    my %promsPath=&promsConfig::getServerInfo('no_user');
    my $dbh = &promsConfig::dbConnect('no_user');
    
    my $sthP=$dbh->prepare("SELECT ID_PROTEIN,IDENTIFIER,PROT_LENGTH,PROT_DES,PROT_SEQ,MW FROM PROTEIN WHERE ID_PROJECT=".$self->{_projectID});
    $sthP->execute;
    my $refProtData=$sthP->fetchall_arrayref; # reference to an array
    $sthP->finish;
    $dbh->disconnect;
    
    $self->appendLog("Scanning project's proteins.");
    
    my $nbLine=0;
    foreach my $refProt (@{$refProtData}) {
        next unless $refProt->[1]; # identifier is null (previous ID protection not removed)
        $projectProt{$refProt->[1]}=$refProt->[0]; # {identifier}=id
        $proteinLength{$refProt->[0]}=$refProt->[2]; #{id}=length
        $proteinDescription{$refProt->[0]}=$refProt->[3] unless !$refProt->[3]; #{id}=description

        if (!$refProt->[2]) { # no length
            $incompleteProtein{$refProt->[1]}=$refProt->[0]; # {identifier}=id
            if (!$refProt->[3] || $refProt->[3]=~/no\sdescription/i) { # no description
                $noDescriptionProtein{$refProt->[1]}=$refProt->[0]; # {identifier}=id
            }
        }
        if (!$refProt->[4]) { # no sequence
            $noSequenceProtein{$refProt->[1]}=$refProt->[0];    # {identifier}=id
        }
        if (!$refProt->[5]) { # no MW
            $noMWProtein{$refProt->[1]}=$refProt->[0];  # {identifier}=id
        }
    }
}

sub _fetchProteinFromTDAFile {
    my ($self) = @_;
    my @analysisList = values (%{$self->{analysisID}});
    my %protList;
    my $dbh = &promsConfig::dbConnect('no_user');
    
    foreach my $analysis (keys %{$self->{_assosProtPeptide}}){
        foreach my $protein (keys %{$self->{_assosProtPeptide}{$analysis}}){
            $protList{$protein}{$self->{analysisID}{$analysis}}=1;
        }
    }
    
    foreach my $dbInfo (@{$self->{_dbID}}) {
        my ($dbID,$dbRank)=split(/&/,$dbInfo);
        &promsMod::getProtInfo('silent', $dbh, $dbID, \@analysisList, \%protDes, \%protMW, \%protOrg, \%protLength, \%protSeq, \%protList);
    
        ## compute peptide specificity and protein coverage
        foreach my $anaName (keys %{$self->{_peptideProt}}){
            foreach my $pep (keys %{$self->{_peptideProt}{$anaName}}){
                next unless scalar keys %{$self->{_peptideProt}{$anaName}{$pep}};
                if (scalar keys %{$self->{_peptideProt}{$anaName}{$pep}} == 1){
                    $self->{_peptideInfo}{$pep}[2]=100;
                }else{
                    $self->{_peptideInfo}{$pep}[2]=(1/scalar keys %{$self->{_peptideProt}{$anaName}{$pep}})*100;
                }
            }
            
            ## protein coverage
            if(!$self->{_protCoverage}) {
                foreach my $protein (keys %{$self->{_peptidePosition}{$anaName}}){
                    next unless $protLength{$protein};
                    my %positionPeptide;
                    foreach my $peptide (keys %{$self->{_peptidePosition}{$anaName}{$protein}}){
                        foreach my $infoBeg (keys %{$self->{_peptidePosition}{$anaName}{$protein}{$peptide}}){
                            my ($beg,$begFlank)=split('_', $infoBeg);
                            my ($end,$endFlank)=split('_', $self->{_peptidePosition}{$anaName}{$protein}{$peptide}{$infoBeg});
                            $positionPeptide{$beg}++;
                            $positionPeptide{$end}--;
                            $numMatchProt{$anaName}{$protein}++;    #number of peptide match
                        }
                    }
                    my $hasPep=0;
                    my ($pepCoverage,$matchBeg);
                    foreach my $position (sort {$a <=> $b} keys %positionPeptide){
                        if($hasPep == 0){
                            $matchBeg=$position;
                        }
                        $hasPep+=$positionPeptide{$position};
                        if($hasPep==0){
                            $pepCoverage+=($position-$matchBeg+1);
                        }
                    }
                    my $coverage=sprintf ("%.1f",$pepCoverage/$protLength{$protein}*100);
                    $self->{_protCoverage}{$anaName}{$protein}=$coverage;
                }
            }
        }
    }
}

sub createQuantification {
    my ($self, $name, $focus, $annot, $method, $relatedAnalysisIDRef, $designID) = @_;
    my $dbh = &promsConfig::dbConnect('no_user');
    my @relatedAnalysisID = @{$relatedAnalysisIDRef};
    
    ## Inserting data into QUANTIFICATION
    my $sthQuanti = $dbh->prepare("INSERT INTO QUANTIFICATION (ID_MODIFICATION,ID_QUANTIFICATION_METHOD,ID_DESIGN,NAME,FOCUS,QUANTIF_ANNOT,STATUS,UPDATE_DATE,UPDATE_USER) VALUES (NULL,?,?,?,?,?,?,NOW(),?)");
    $sthQuanti->execute($method, $designID, $name, $focus, $annot, 1, ($self->{_userID}) ? $self->{_userID} : undef) or die $dbh->errstr;
    $sthQuanti->finish;
    my $quantiID = $dbh->last_insert_id(undef, undef, 'QUANTIFICATION', 'ID_QUANTIFICATION');
    
    ## Inserting data into ANA_QUANTIFICATION
    my $sthAnaQuanti = $dbh->prepare("INSERT INTO ANA_QUANTIFICATION (ID_QUANTIFICATION,ID_ANALYSIS,QUANTIF_FILE,IS_REFERENCE) VALUES (?,?,NULL,NULL)");
    foreach my $analysisID (@relatedAnalysisID) {
        $sthAnaQuanti->execute($quantiID, $analysisID);
    }
    $sthAnaQuanti->finish;
    
    $dbh->commit;
    $dbh->disconnect;
    return $quantiID;
}

# Need to have all IDs of analysis and observations computed with the importData function if not specified in $designExperiment{"conditions"}. expID comes from the one provided in object constructor
sub createExperimentDesign { # $design{"name"} + $design{"conditions"}{CONDITION_NAME}{REPLICATE_NAME} = {anaName => ANA_NAME, techRep => TECH_REP_NB, fraction => FRACTION_NB, anaID => ANA_ID, obsID => OBS_ID}
    my ($self, $designExperimentRef) = @_;
    my %designExperiment = ($designExperimentRef) ? %{$designExperimentRef} : %{$self->{designExperiment}};
    my %designExperimentID;
 
    if(%designExperiment && scalar keys %designExperiment > 1) {
        my $dbh = &promsConfig::dbConnect('no_user');
        my $sthDesign = $dbh->prepare("INSERT INTO DESIGN (ID_DESIGN,ID_EXPERIMENT,NAME,DES,UPDATE_DATE,UPDATE_USER) VALUES (?,?,?,?,NOW(),?)");
        my $sthState = $dbh->prepare("INSERT INTO EXPCONDITION (ID_EXPCONDITION,ID_DESIGN,NAME,DES) VALUES (?,?,?,?)");
        my $sthInsObsExp = $dbh->prepare("INSERT INTO OBS_EXPCONDITION (ID_EXPCONDITION,ID_OBSERVATION,FRACTION_GROUP,TECH_REP_GROUP) VALUES (?,?,?,?)");
        
        # Create empty Design
        my ($newDesignID) = $dbh->selectrow_array("SELECT MAX(ID_DESIGN) FROM DESIGN");
        $sthDesign->execute(++$newDesignID, $self->{_experimentID}, ($designExperiment{"name"}) ? $designExperiment{"name"} : "Quantification Design", "Automatically generated design", $self->{_userID}) || die $dbh->errstr();
        $designExperimentID{"ID"} = $newDesignID;
        my $stateIndex = 1;
        
        # Get max expConditionID (= State ID)
        my ($newExpConditionID) = $dbh->selectrow_array("SELECT MAX(ID_EXPCONDITION) FROM EXPCONDITION");
        
        foreach my $condition (sort { lc($a) cmp lc($b) } keys %{$designExperiment{"conditions"}}){
            
            # Add State
            if ($condition) {
                $sthState->execute(++$newExpConditionID, $newDesignID, $condition, "") || die $dbh->errstr();
                $designExperimentID{"STATE"}{$newExpConditionID}{"NAME"} = $condition; 
                $designExperimentID{"STATE"}{$newExpConditionID}{"INDEX"} = $stateIndex; 
                %{$designExperimentID{"STATE"}{$newExpConditionID}{"OBS"}} = ();
                
                foreach my $replicate (sort { $a <=> $b } keys %{$designExperiment{"conditions"}{$condition}}){
                    my %replicateInfos = %{$designExperiment{"conditions"}{$condition}{$replicate}};
                    my ($anaName, $techRep, $fraction) = ($replicateInfos{"anaName"}, $replicateInfos{"techRep"}, $replicateInfos{"fraction"});
                    my $obsID = ($replicateInfos{"obsID"}) ? $replicateInfos{"obsID"} : ($self->{anaObsID}) ? $self->{anaObsID}{$anaName} : undef;            
                    my $anaID = ($replicateInfos{"anaID"}) ? $replicateInfos{"anaID"} : ($self->{analysisID}) ? $self->{analysisID}{$anaName} : undef;            
    
                    # Link observations to experiment
                    if($obsID) {
                        $sthInsObsExp->execute($newExpConditionID, $obsID, ($fraction) ? $fraction : undef, ($techRep) ? $techRep : undef) || die $dbh->errstr();
                        $designExperimentID{"STATE"}{$newExpConditionID}{"OBS"}{$obsID} = $anaID;
                    }
                }
                
                $stateIndex++;
            }
        }
        
        $sthDesign->finish;
        $sthState->finish;
        $sthInsObsExp->finish;
        
        $dbh->commit;
        $dbh->disconnect;
        
        return \%designExperimentID;
    }
    
    return -1;
}

1;

####>Revision history<####
# 1.0.18 [BUGFIX] Count inserted peptides and proteins peptide association properly (VS 01/04/21) 
# 1.0.17 [BUGFIX] Fixed strict-related error (VS 02/02/21) 
# 1.0.16 [MINOR] Added log of proteins inserted in db (VS 29/01/21)
# 1.0.15 [MINOR] Avoid lock timeout by regularly committing changes during peptides and peptide_attribs insert (VS 11/12/2020)
# 1.0.14 [ENHANCEMENT] Handles imputed data from Spectronaut (VS 26/11/20)
# 1.0.13 [MINOR] Use &mrCalc from promsMod (VS 22/10/20)
# 1.0.12 [MINOR] Added spectronaut peptide sequence modifications pattern (VS 25/09/20)
# 1.0.11 [MINOR] Changed data sorting for sample/analysis insertion (VS 02/09/20)
# 1.0.10 [MODIF] Changed quantification design storage structure (VS 21/08/20)  
# 1.0.9 [BUGFIX] Do not insert all proteins for all analysis. Missing proteins where considered has ghost (VS 12/08/20)
# 1.0.8 [MODIF] Minor modif on skyline fileFormat to be consistent with handling of .sky files (VL 20/11/19)
# 1.0.7 [BUGFIX] Loop over proteinID keys instead of dbRankProt (dbRankProt misses keys in DIA) (VS 16/10/19) 
# 1.0.6 [MODIF] Fix unused name dbRankProtRef to dbRankProt and loop over its keys in importData to insert proteins in ANALYSIS_PROTEIN (VL 04/09/19) 
# 1.0.5 [MODIF] Changed createMatchGroup to fit with promsMod.pm prototype (VS 03/07/19)
# 1.0.4 [MODIF] Add RT to peptideID list returned after insertion (VS 02/07/19)
# 1.0.3 [MODIF] Add protein mapping after creation (VS 01/07/19)
# 1.0.2 [BUGFIX] Stabilize DIA data import (VS 11/06/19)
# 1.0.1 Added RT mean computation for TDA MS2 peptides (VS 27/11/18)
# 1.0.0 Creation (VS 12/11/18)
