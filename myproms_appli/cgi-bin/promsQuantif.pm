################################################################################
# promsQuantif.pm           1.3.4                                              #
# Authors: P. Poullet, G. Arras, S. Liva (Institut Curie)                      #
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
package promsQuantif;
require Exporter;

@ISA=qw(Exporter);
@EXPORT=qw();
@EXPORT_OK=qw();
$VERSION=1.00;

#use POSIX qw(strftime);
use strict;
use File::Path qw(rmtree); # remove_tree

my $MAX_INF_RATIO_DB=1000;
my $MIN_INF_RATIO_DB=0.001;
my $MAX_INF_RATIO=100; # 1000    ~~~Max ratio allowed before switching to infinite~~~
my $MIN_INF_RATIO=1/$MAX_INF_RATIO;

sub getExtremeRatios {
	return ($MIN_INF_RATIO,$MAX_INF_RATIO);
}

sub getProteinQuantifFamilies {
	my %proteinQuantifFamilies=('MEMBERS'=>{'RATIO'=>['PROT_RATIO_PEP','TNPQ','DIA'], # <- there should be more members soon
										'EMPAI'=>['EMPAI'],
										'SIN'=>['SIN'],
										'MQ'=>['MQ'],
									   },
								'NAME'=>{'RATIO'=>'Protein ratio',
										 'EMPAI'=>'emPAI',
										 'SIN'=>'Normalized Spectral Index',
										 'MQ'=>'MaxQuant intensities',
										},
								'MEASURES'=>{
									'MQ'=>[['MQ_INT','Intensity'],['MQ_IBAQ','iBAQ'],['MQ_LFQ','LFQ'],['MQ_SC','MS/MS count']],
									'EMPAI'=>[['EMPAI','emPAI'],['EMPAI_MOL','emPAI (Mol %)'],['EMPAI_MR','emPAI (Mr %)']],
									'SIN'=>[['SIN_SIN','SI<sub>N</sub>']]
								}
								);
	return %proteinQuantifFamilies;
}

sub getXicSoftware { # Only for peptide quantification!!!
	my ($dbh,$quantifID,$quantifAnnot)=@_; # quantifAnnot is optional
	my %code2Name=(PD=>'Proteome Discoverer',MCQ=>'MassChroQ',MAS=>'Mascot',PAR=>'Paragon',MQ=>'MaxQuant',
				   PKV=>'PeakView',OS=>'OpenSWATH',SKY=>'Skyline',
				   #''=>'Unknown'
				   );
	unless ($quantifAnnot) {
		($quantifAnnot)=$dbh->selectrow_array("SELECT QUANTIF_ANNOT FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$quantifID");
	}
	my ($xicSoftCode)=$quantifAnnot=~/::SOFTWARE=([^:]+)/; # remove software info for back compatibility
	my $version=0;
	if ($xicSoftCode) {
		($xicSoftCode,$version)=split(';',$xicSoftCode);
		$version=0 unless $version;
	}
	else {
		if ($quantifAnnot=~/EXTRACTION_ALGO=/) {$xicSoftCode='MCQ'; $version='2.0.1';}
		else {
			my ($fileFormat)=$dbh->selectrow_array("SELECT FILE_FORMAT FROM ANALYSIS A,ANA_QUANTIFICATION Q WHERE A.ID_ANALYSIS=Q.ID_ANALYSIS AND ID_QUANTIFICATION=$quantifID LIMIT 0,1");
			$xicSoftCode=($fileFormat=~/\.PDM\Z/)? 'PD' : ($fileFormat=~/^MASCOT/)? 'MAS' : ($fileFormat=~/PARAGON/)? 'PAR' : '?';
		}
	}
	return ($xicSoftCode,$code2Name{$xicSoftCode} || 'Unknown',$version);
}


sub getQuantifNormalizationName {
	my %normalizationNames=(
		#>Algo v3--->
		'none.none'=>		'None',
		'loess.none'=>		'Loess',
		'loess.scale'=>		'Loess & Scale',
		'median.none'=>		'Median',
		'median.scale'=>	'Median & Scale',
		'quantile'=>		'Quantile',
		# Options below are possible but not used
		#'none.scale'=>	'None & Scale',
		#'mean.none'=>	'Mean',
		#'mean.scale'=>	'Mean & Scale',

		#>Algo v2 (old)--->
		'loess.none.normalization'=>	'Loess',
		'loess.scale.normalization'=>	'Loess & Scale',
		'median.none.normalization'=>	'Median',
		'median.scale.normalization'=>	'Median & Scale',
		'quantile.normalization'=>		'Quantile',
		'global.mad.normalization'=>	'Global MAD',
		'global.normalization.mean'=>	'Global mean',
		'global.normalization.median'=>	'Global median',
		'global.mad.normalization'=>	'Global MAD'
	);
	return %normalizationNames;
}

###################################
####<<<Delete quantification>>>#### Use only if parent Analys(i/e)s will deleted too (Does not check for quantif-specific ghosts peptides)
###################################
sub deleteQuantification {

	my ($dbh,$projectID,$qID)=@_;
	my %promsPath=&promsConfig::getServerInfo; # WARNING: promsConfig must be declared in calling script!!!

	my ($focus)=$dbh->selectrow_array("SELECT FOCUS FROM QUANTIFICATION,QUANTIFICATION_METHOD WHERE QUANTIFICATION.ID_QUANTIFICATION_METHOD=QUANTIFICATION_METHOD.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$qID");

	if ($focus eq 'protein') {
		###<Delete links to parent quantifs...
		$dbh->do("DELETE FROM PARENT_QUANTIFICATION WHERE ID_QUANTIFICATION=$qID") || die $dbh->errstr();

		###<Delete PROTEIN_QUANTIFICATION & associated tables: PROTQUANTIF_MODRES & MODIFIED_RESIDUE (with JOIN!!!)
		$dbh->do("DELETE PQM FROM PROTQUANTIF_MODRES PQM INNER JOIN MODIFIED_RESIDUE R ON PQM.ID_MODIF_RES=R.ID_MODIF_RES WHERE R.ID_QUANTIFICATION=$qID");
		$dbh->do("DELETE FROM MODIFIED_RESIDUE WHERE ID_QUANTIFICATION=$qID");
		$dbh->do("DELETE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=$qID") || die $dbh->errstr();
	}
	else { # peptide
		#$dbh->do("DELETE FROM PEPTIDE_QUANTIFICATION WHERE ID_QUANTIFICATION=$qID") || die $dbh->errstr();
		###<Delete QUANTIF_REFRT information
		$dbh->do("DELETE FROM QUANTIF_REFRT WHERE ID_QUANTIFICATION=$qID") || die $dbh->errstr();
	}

	###<Delete ANA_QUANTIFICATION information
	$dbh->do("DELETE FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$qID") || die $dbh->errstr();

	###<Delete EXPCONDITION_QUANTIF information
	$dbh->do("DELETE FROM EXPCONDITION_QUANTIF WHERE ID_QUANTIFICATION=$qID") || die $dbh->errstr();

	###<Delete the directory and log file if one was created
	#remove_tree("$promsPath{quantification}/project_$projectID/quanti_$qID") if -e "$promsPath{quantification}/project_$projectID/quanti_$qID";
	rmtree("$promsPath{quantification}/project_$projectID/quanti_$qID") if -e ("$promsPath{quantification}/project_$projectID" && "$promsPath{quantification}/project_$projectID/quanti_$qID");
	unlink "$promsPath{logs}/quanti_$qID.log" if -e "$promsPath{logs}/quanti_$qID.log";

	###<Delete the QUANTIFICATION itself
	$dbh->do("DELETE FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$qID") || die $dbh->errstr();

}

#########################################################################
####<<<Get list of quantification parameters from quantif_info.txt>>>####
#########################################################################
sub getQuantificationParameters {
	my ($infoFile)=@_;
	my %params;
	open (INFO,$infoFile);
	my $section='';
	my $numConds=0;
	while (<INFO>) {
		if (/^PARAMETERS:/) {
			$section='parameters';
			next;
		}
		elsif (/^QUANTIFICATIONS:/) {
			$section='quantifications';
			next;
		}
		last if /^ANALYSES/;
		if ($section eq 'parameters' && /\w/) {
			chomp;
			my ($nameDB,$nameR,@paramValues)=split(/\t/,$_); # name-in-DB	name-in-R	value	extra-value
			@{$params{'DB'}{$nameDB}}=@paramValues if $nameDB;
			if ($nameR) {
				if ($nameR=~/name\.(grp|ratio)/) { # internal quantif only
				#	$paramValues[0]=".$paramValues[0]" if $paramValues[0]=~/^\d/; # problem with R adding X at begining of string if starts with digit (eg. 114)
				#	$paramValues[0]=&shortenName($paramValues[0],28); # max size allowed for figure annotation
					my ($statePos)=($nameDB=~/CONDITION_(\d+)/);
					push @{$params{'R'}{$nameR}},"State$statePos";
				}
				else {
					push @{$params{'R'}{$nameR}},$paramValues[0]; # same name can be found on multiple lines
				}
			}
			$numConds++ if $nameDB=~/CONDITION_\d+/; # conditions are listed in ascending order with no gaps (1,2,3,...,n)
		}
	}
	close INFO;
	$params{'DB'}{'FDR_ALPHA'}[0]*=100 if $params{'DB'}{'FDR_ALPHA'}; # converts % for DB (fraction for R)

if ($numConds) { # Ratios (n+1)/n --> for n conditions: num ratios = n*(n-1)/2
	foreach my $d (1..$numConds-1) {
		foreach my $n ($d+1..$numConds) {
			push @{$params{'DB'}{'RATIOS'}},"$n/$d";
		}
	}
}
	return %params;
}

######################################################################
####<<<Writes R parameters for Quantification.Function.R script>>>####
######################################################################
sub writeQuantifParameterFiles {
	my ($dataDir,$refParams)=@_;

	#<Character
	my @paramCharList=('normalization.method','pAdj.method','design','residual.variability'); # Super/SimpleRatio & LabelFree
	push @paramCharList,('quantification.method','bias.correction','name.grp','name.ratio','invRatio','alter','pAdj','typeTest','typeTab','metric','prot.ref'); # Ratio (old algo)
	push @paramCharList,('Design'); # for algo v2 + 'normalization.method','pAdj.method','pAdj'
	#push @paramCharList,('savegraph','displaygraph','Design','normalization.method','filepath'); # (Super/Simple)Ratio & LabelFree
	push @paramCharList,('contrasts.matrix','clusters'); # MSstats
	my %usedParams;
	my $matched=0;
	foreach my $param (@paramCharList) {
		next unless $refParams->{$param};
		next if $usedParams{$param}; # just to be safe (duplicate-proof!)
		$matched++;
		open(PARAM_CHAR,">$dataDir/param_char.txt") if $matched==1; # create file only if not empty
		foreach my $value (@{$refParams->{$param}}) {
			#if ($param=~/^name\./) { # name.grp, name.ratio
			#	$value=".$value" if $value=~/^\d/; # problem with R adding X at begining of string if starts with digit
			#	$value=&shortenName($value,28); # max size allowed for figure annotation
			#}
			print PARAM_CHAR "$param\t$value\n";
		}
		$usedParams{$param}=1;
	}
	close PARAM_CHAR if $matched;

	#<Numerical
	$matched=0;
	foreach my $param ('sup.var','grp','int.min','int.max','threshold.CV','lim.PERM','alpha','threshold.out','denom','conf.level') { # lim.PERM obsolete
		next unless $refParams->{$param};
		next if $usedParams{$param}; # just to be safe (duplicate-proof!)
		$matched++;
		open(PARAM_NUM,">$dataDir/param_numeric.txt") if $matched==1; # create file only if not empty
		foreach my $value (@{$refParams->{$param}}) {
			print PARAM_NUM "$param\t$value\n";
		}
		$usedParams{$param}=1;
	}
	close PARAM_NUM if $matched;
}


sub getDistinctQuantifNames { # TODO: test on mixed ratio & non-ratio quantifs!!!
	my ($dbh,$quantifListStrg,$customKey)=@_; # $quantifListStrg can be 'quantifID1_pos1:quantifID1_pos2:...' OR 'quantifID1_pos1_pos2:...'
	my @customKeys=($customKey)? split(':',$customKey) : (); # optional parameter
	my (%parentNames,%quantifNames,%labelingInfo,%stateInfo,%quantifRatioNames,%parents,%quantifs,%references,%ratios);
	my $maxTargetPosInQuantif=0;
	my %sthQuantif = (D=>$dbh->prepare("SELECT D.NAME,Q.NAME,Q.QUANTIF_ANNOT FROM DESIGN D,QUANTIFICATION Q WHERE D.ID_DESIGN=Q.ID_DESIGN AND Q.ID_QUANTIFICATION = ?"),
					  A=>$dbh->prepare("SELECT A.NAME,Q.NAME,Q.QUANTIF_ANNOT FROM ANALYSIS A,ANA_QUANTIFICATION AQ,QUANTIFICATION Q WHERE A.ID_ANALYSIS=AQ.ID_ANALYSIS AND AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND ID_DESIGN IS NULL AND Q.ID_QUANTIFICATION = ?")
					 );
	my $numQuantifTgPos=0;
	foreach my $quantif (split(':',$quantifListStrg)) {
		my ($quantifID,@targetPositions)=split('_',$quantif);
		unless ($quantifNames{$quantifID}) { # do not query twice
			$sthQuantif{D}->execute ($quantifID);
			my $quantifAnnot;
			($parentNames{$quantifID},$quantifNames{$quantifID},$quantifAnnot)=$sthQuantif{D}->fetchrow_array;
			unless ($parentNames{$quantifID}) { # Internal-analysis
				$sthQuantif{A}->execute ($quantifID);
				($parentNames{$quantifID},$quantifNames{$quantifID},$quantifAnnot)=$sthQuantif{A}->fetchrow_array;
			}
			if ($targetPositions[0] && $targetPositions[0] > 0 && !$labelingInfo{$quantifID}) { # skip non-target/channel quantifs (EMPAI,SIN)
				%{$labelingInfo{$quantifID}}=%{$stateInfo{$quantifID}}=();
				&extractQuantificationParameters($dbh,$quantifAnnot,$labelingInfo{$quantifID},$stateInfo{$quantifID}); # loads %{$labelingInfo{$quantifID}} & %{$stateInfo{$quantifID}}
			}
			$parents{$parentNames{$quantifID}}=1; # counts distinct parents names
			$quantifs{$quantifNames{$quantifID}}=1; # counts distinct quantifs names
			my $numTargetPos=($labelingInfo{$quantifID}{'RATIOS'})? scalar @{$labelingInfo{$quantifID}{'RATIOS'}} : ($labelingInfo{$quantifID}{'STATES'})? scalar @{$labelingInfo{$quantifID}{'STATES'}} : 0;
			$maxTargetPosInQuantif=$numTargetPos if $maxTargetPosInQuantif < $numTargetPos;
		}
		foreach my $targetPos (@targetPositions) {
			$numQuantifTgPos++;
			my $quantif=$quantifID.'_'.$targetPos;
			##<Ratio quantif
			my (@info,@junctions);
			if ($targetPos > 0) { # targetPos quantifs
				if ($labelingInfo{$quantifID}{'RATIOS'}) { # ratio quantifs
					my ($testStatePos,$refStatePos)=split(/\//,$labelingInfo{$quantifID}{'RATIOS'}[$targetPos-1]);
					my $normTag='';
					if ($testStatePos=~/%/) { # Super ratio
						$normTag='°';
						$testStatePos=~s/%\d+//;
						$refStatePos=~s/%\d+//;
					}
					if ($customKey) {
						@info=($parentNames{$quantifID},$quantifNames{$quantifID},$stateInfo{$quantifID}{$testStatePos}{NAME}.$normTag,$stateInfo{$quantifID}{$refStatePos}{NAME}.$normTag); # Design,quantif,ratioTest,ratioRef
						@junctions=(' > ',' : ','/');
					}
					$quantifRatioNames{'FULL'}{$quantif}="$parentNames{$quantifID} > $quantifNames{$quantifID} : $stateInfo{$quantifID}{$testStatePos}{NAME}$normTag/$stateInfo{$quantifID}{$refStatePos}{NAME}$normTag";
					$quantifRatioNames{'EXTENDED'}{$quantif}="$quantifNames{$quantifID} : $stateInfo{$quantifID}{$testStatePos}{NAME}$normTag/$stateInfo{$quantifID}{$refStatePos}{NAME}$normTag";
					$quantifRatioNames{'RATIO'}{$quantif}="$stateInfo{$quantifID}{$testStatePos}{NAME}$normTag/$stateInfo{$quantifID}{$refStatePos}{NAME}$normTag";
					$quantifRatioNames{'TEST'}{$quantif}="$stateInfo{$quantifID}{$testStatePos}{NAME}$normTag~";
					#$quantifRatioNames{'QUANTIF'}{$quantif}=$quantifNames{$quantifID};
					#$quantifRatioNames{'PARENT'}{$quantif}=$parentNames{$quantifID};
					$references{"$stateInfo{$quantifID}{$refStatePos}{NAME}$normTag"}=1; # counts distinct refs
					$ratios{"$stateInfo{$quantifID}{$testStatePos}{NAME}$normTag/$stateInfo{$quantifID}{$refStatePos}{NAME}$normTag"}=1;
				}
				else { # intensity quantifs
					if ($customKey) {
						@info=($parentNames{$quantifID},$quantifNames{$quantifID},$stateInfo{$quantifID}{$targetPos}{NAME}); # Design,quantif,ratioTest
						@junctions=(' > ',' : ');
					}
					$quantifRatioNames{'FULL'}{$quantif}="$parentNames{$quantifID} > $quantifNames{$quantifID} : $stateInfo{$quantifID}{$targetPos}{NAME}";
					$quantifRatioNames{'EXTENDED'}{$quantif}="$quantifNames{$quantifID} : $stateInfo{$quantifID}{$targetPos}{NAME}";
					$quantifRatioNames{'RATIO'}{$quantif}=$quantifRatioNames{'STATE'}{$quantif}=$quantifRatioNames{'TEST'}{$quantif}=$stateInfo{$quantifID}{$targetPos}{NAME};
				}
				$quantifRatioNames{'QUANTIF'}{$quantif}=$quantifNames{$quantifID};
				$quantifRatioNames{'PARENT'}{$quantif}=$parentNames{$quantifID};
			}
			##<Non-targetPos quantif
			else {
				if ($customKey) {
					@info=($quantifNames{$quantifID}); # no_parent,quantif
				}
				$quantifRatioNames{'FULL'}{$quantif}=$quantifRatioNames{'QUANTIF'}{$quantif}=$quantifRatioNames{'RATIO'}{$quantif}=$quantifRatioNames{'TEST'}{$quantif}=$quantifNames{$quantifID};
				#$quantifs{$quantifID}=1; # counts distinct quantifs
			}
			if ($customKey) {
				$quantifRatioNames{'CUSTOM'}{$quantif}='';
				my $j=-1; # previously used key index
				foreach my $i (0..$#customKeys) {
					next unless $customKeys[$i];
					$quantifRatioNames{'CUSTOM'}{$quantif}.=$junctions[$j] if ($j>=0 && $customKeys[$j]);
					$quantifRatioNames{'CUSTOM'}{$quantif}.=$info[$i].$junctions[$i];
					$j=$i;
				}
			}
		}
	}
	$sthQuantif{D}->finish;
	$sthQuantif{A}->finish;

	##<OPTIMAL name>## (based on item names NOT ids!)
	my $numQuantifNames=scalar keys %quantifs;
	if ($numQuantifNames == $numQuantifTgPos && $maxTargetPosInQuantif==1) { # All quantif ratios come from quantifs with different names && only 1 ratio each
		%{$quantifRatioNames{'OPTIMAL'}}=%{$quantifRatioNames{'QUANTIF'}};
	}
	elsif ($numQuantifTgPos && scalar keys %ratios == $numQuantifTgPos) { # All quantif ratios have different names
		if (scalar keys %references <= 1) { # same reference name for all or no ratio quantif -> skip
			%{$quantifRatioNames{'OPTIMAL'}}=%{$quantifRatioNames{'TEST'}};
		}
		else {
			%{$quantifRatioNames{'OPTIMAL'}}=%{$quantifRatioNames{'RATIO'}};
		}
	}
	elsif ($numQuantifTgPos) {
		%{$quantifRatioNames{'OPTIMAL'}}=%{$quantifRatioNames{'FULL'}};
		foreach my $quantif (split(':',$quantifListStrg)) {
			my ($quantifID,@targetPositions)=split('_',$quantif);
			foreach my $ratioPos (@targetPositions) {
				my $quantif=$quantifID.'_'.$ratioPos;
				if (scalar keys %parents == 1) {
					my $qParName=quotemeta($parentNames{$quantifID});
					$quantifRatioNames{'OPTIMAL'}{$quantif}=~s/$qParName > //;
				}
				if (scalar keys %quantifs == 1) {
					my $qQuantName=quotemeta($quantifNames{$quantifID});
					$quantifRatioNames{'OPTIMAL'}{$quantif}=~s/$qQuantName : //;
				}
				if (scalar keys %ratios == 1) {
					my $qRatioName=quotemeta($quantifRatioNames{'RATIO'}{$quantif});
					$quantifRatioNames{'OPTIMAL'}{$quantif}=~s/\s*[>:]*\s*$qRatioName//;
				}
				elsif (scalar keys %references == 1) {
					my $qTestName=quotemeta($quantifRatioNames{'TEST'}{$quantif});
					$quantifRatioNames{'OPTIMAL'}{$quantif}=~s/\/$qTestName//;
				}
			}
		}
	}
	if ($customKey) {
		foreach my $quantif (split(':',$quantifListStrg)) {
			$quantifRatioNames{'CUSTOM'}{$quantif}=$quantifRatioNames{'OPTIMAL'}{$quantif} unless $quantifRatioNames{'CUSTOM'}{$quantif}; # just to be safe
		}
	}

	return %quantifRatioNames;
}

sub extractQuantificationParameters {
	my ($dbh,$quantifAnnot,$refLabelingInfo,$refStateInfo)=@_;

	#my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
	my (@labelInfo)=split('::',$quantifAnnot);
	foreach my $infoStrg (@labelInfo) {
		my ($setting,$valueStrg)=split('=',$infoStrg);
		@{$refLabelingInfo->{$setting}}=split(';',$valueStrg);
	}
	my $numStates=0;
	if ($refLabelingInfo->{'STATES'}[0]=~/\d,#/) { # design quantif
		my $sthExpCondName=$dbh->prepare("SELECT NAME FROM EXPCONDITION WHERE ID_EXPCONDITION=?");
		if ($refLabelingInfo->{'RATIOS'}) { # ratio quantif
			for (my $i=0;$i<=$#{$refLabelingInfo->{'RATIOS'}}; $i++) {
				$refLabelingInfo->{'RATIOS'}[$i]=~s/#//g; # remove id flags
				my @ratioConds=split(/\//,$refLabelingInfo->{'RATIOS'}[$i]);
				foreach my $condID (@ratioConds) {
					last if $condID=~/%/; # Super ratio
					$sthExpCondName->execute($condID);
					($refStateInfo->{$condID}{'NAME'})=$sthExpCondName->fetchrow_array;
				}
			}
		}
		else { # state quantif (MaxQuant)
			foreach my $stateData (@{$refLabelingInfo->{STATES}}) { # 1,#1:#1,#1      1,#1:#2,#2
				$numStates++;
				$stateData=~s/#//g; # remove id flags
				my ($numRep,$replicates,$expCondID)=split(',',$stateData);
				$sthExpCondName->execute($expCondID);
				my ($expCondName)=$sthExpCondName->fetchrow_array;
				$refStateInfo->{$numStates}{'NAME'}=$expCondName;
			}
		}
		$sthExpCondName->finish;
	}
	else { # internal ratio quantif
		foreach my $stateData (@{$refLabelingInfo->{'STATES'}}) {
			$numStates++;
			(my $numRep,$refStateInfo->{$numStates}{'NAME'},my $repPosStrg)=split(',',$stateData);
			$refStateInfo->{$numStates}{'NAME'}=~s/\./\+/g;
		}
	}
}


############################################################################################
#######################<Quantification data collection routine>#############################
############################################################################################
sub fetchQuantificationData {
	my ($dbh,$refParams,$refQuantifInfo,$refQuantifValues,$refProteinInfo,$refSelectedElements,$refExcludedProteins)=@_; # $refExcludedProteins optional
	my $getProtInfo=1;
	unless ($refProteinInfo) {
		$refProteinInfo={}; # defined to prevent uninitialized value later on
		$getProtInfo=0;
	}
	my $verbose=$refParams->{VERBOSE} || 0;
	my @selectedQuantifications=@{$refParams->{QUANTIF_LIST}};
	my $selQuantifFamily;
	if ($refParams->{'QUANTIF_FAMILY'}) {$selQuantifFamily=$refParams->{'QUANTIF_FAMILY'};} # optional
	else {
		my ($quantifID)=(split('_',$selectedQuantifications[0]))[0];
		my ($qMethCode)=$dbh->selectrow_array("SELECT QM.CODE FROM QUANTIFICATION Q,QUANTIFICATION_METHOD QM WHERE Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=$quantifID");
		$selQuantifFamily=($qMethCode=~/PROT_RATIO_PEP|TNPQ/)? 'RATIO' : $qMethCode;
	}
	($refParams->{MIN_RATIO},$refParams->{MAX_RATIO},$refParams->{MIN_PVALUE})=(1000000,0.000001,1); # modified only when RATIO for view = log2,heatmap,explorAna and returned
	my %quantifParamInfo;
	my $view=$refParams->{'VIEW'} || 'log2';

	my ($protSelection,$isoSelection)=(0,0); # default
	my ($refSelectedProteins,$refSelectedIsoforms)=({},{});
	if ($refSelectedElements && keys %{$refSelectedElements}) {
		my $key=(keys %{$refSelectedElements})[0];
		if ($key && $key=~/\D/) { # isoform
			#<Check if a Protein quantif is used
			my %quantifs;
			foreach my $quantif (@selectedQuantifications) {
				my ($quantifID,$ratioPos)=split('_',$quantif);
				$quantifs{$quantifID}=1;
			}
			my ($existProtQuantif)=$dbh->selectrow_array("SELECT COUNT(*) FROM QUANTIFICATION WHERE ID_MODIFICATION IS NULL AND ID_QUANTIFICATION IN (".join(',',keys %quantifs).")");
			if ($existProtQuantif) {
				$protSelection=1;
				foreach my $modProtID (keys %{$refSelectedElements}) {
					my ($protID)=$modProtID=~/^(\d+)/;
					$refSelectedProteins->{$protID}=1;
				}
				#<Return info
				$refParams->{RETURN}='WARNING: Site restriction list was downgraded to standard protein list';
			}
			else {
				$isoSelection=1;
				$refSelectedIsoforms=$refSelectedElements;
			}
		}
		else {
			$protSelection=1;
			$refSelectedProteins=$refSelectedElements;
		}
	}
	my $restrictStrg=($protSelection && scalar keys %{$refSelectedProteins} <= 25)? 'AND PQ.ID_PROTEIN IN ('.join(',',keys %{$refSelectedProteins}).')' : '';
	my $exclusion=($refExcludedProteins && scalar keys %{$refExcludedProteins})? 1 : 0;
	#my ($quantifMethod)=$dbh->selectrow_array("SELECT CODE FROM QUANTIFICATION_METHOD WHERE ID_QUANTIFICATION_METHOD=$selQuantifMethodID"); # TO BE changed to multi-method label (multiple methods could generate protein ratios)

	####<RATIO|MQ>####
	if ($selQuantifFamily=~/RATIO|MQ/) {

		my $numPepCode=$refParams->{'NUM_PEP_CODE'};

		my $sthQinfo=$dbh->prepare("SELECT NAME,QUANTIF_ANNOT,ID_QUANTIFICATION_METHOD,ID_MODIFICATION FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
		#my $sthAna=$dbh->prepare("SELECT GROUP_CONCAT(ID_ANALYSIS SEPARATOR ',') FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=? GROUP BY ID_QUANTIFICATION"); <-- truncated if too many ana!!!
		my $sthAna=$dbh->prepare('SELECT ID_ANALYSIS FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=?');
		my $sthQMP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,NAME,CODE FROM QUANTIFICATION_PARAMETER WHERE ID_QUANTIFICATION_METHOD=?");
		#my $sthProtQ0=$dbh->prepare("SELECT ID_PROTEIN,QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=? AND ID_QUANTIF_PARAMETER=?");
		my $sthProtQ0=$dbh->prepare("SELECT PQ.ID_PROTEIN,GROUP_CONCAT(RESIDUE,POSITION ORDER BY POSITION SEPARATOR '.'),QUANTIF_VALUE
										FROM PROTEIN_QUANTIFICATION PQ
										LEFT JOIN PROTQUANTIF_MODRES PQMR ON PQ.ID_PROT_QUANTIF=PQMR.ID_PROT_QUANTIF
										LEFT JOIN MODIFIED_RESIDUE MR ON PQMR.ID_MODIF_RES=MR.ID_MODIF_RES AND MR.ID_QUANTIFICATION=PQ.ID_QUANTIFICATION
										WHERE PQ.ID_QUANTIFICATION=? AND ID_QUANTIF_PARAMETER=? $restrictStrg GROUP BY PQ.ID_PROT_QUANTIF");
		#my $sthProtQ=$dbh->prepare("SELECT ID_PROTEIN,QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION WHERE ID_QUANTIFICATION=? AND TARGET_POS=? AND ID_QUANTIF_PARAMETER=?");
		my $sthProtQ=$dbh->prepare("SELECT PQ.ID_PROTEIN,GROUP_CONCAT(RESIDUE,POSITION ORDER BY POSITION SEPARATOR '.'),QUANTIF_VALUE
										FROM PROTEIN_QUANTIFICATION PQ
										LEFT JOIN PROTQUANTIF_MODRES PQMR ON PQ.ID_PROT_QUANTIF=PQMR.ID_PROT_QUANTIF
										LEFT JOIN MODIFIED_RESIDUE MR ON PQMR.ID_MODIF_RES=MR.ID_MODIF_RES AND MR.ID_QUANTIFICATION=PQ.ID_QUANTIFICATION
										WHERE PQ.ID_QUANTIFICATION=? AND TARGET_POS=? AND ID_QUANTIF_PARAMETER=? $restrictStrg GROUP BY PQ.ID_PROT_QUANTIF");
		my $modifPrefixStrg='';	# adds [modifID] before modif site unless $refParams->{SEL_MOD_ID} is set
		foreach my $quantif (@selectedQuantifications) {
			print '.' if $verbose;
			$refParams->{MINUS_INF}{$quantif}=$refParams->{PLUS_INF}{$quantif}=0; # records existence of +/-inf values
			my ($quantifID,$ratioPos)=split('_',$quantif);
			unless ($refQuantifInfo->{$quantifID}) { # Needed only once for a given quantifID
				$sthQinfo->execute($quantifID);
				my ($quantifName,$quantifAnnot,$quantifMethID,$modifID)=$sthQinfo->fetchrow_array;
				$modifID=0 unless $modifID;
				$modifPrefixStrg=($modifID && (!$refParams->{SEL_MODIF_ID} || $refParams->{SEL_MODIF_ID} != $modifID))? "[$modifID]" : '';
				my (%labelingInfo,%stateInfo);
				&extractQuantificationParameters($dbh,$quantifAnnot,\%labelingInfo,\%stateInfo);
				$sthAna->execute($quantifID);
				my @anaIDList;
				while (my ($anaID)=$sthAna->fetchrow_array) {push @anaIDList,$anaID;}
				my @itemInfo=&promsMod::getItemInfo($dbh,'QUANTIFICATION',$quantifID);
				@{$refQuantifInfo->{$quantifID}}=($quantifName,\%labelingInfo,\%stateInfo,\@itemInfo,$modifID,join(',',@anaIDList)); # no Analysis in itemInfo if quantif is from Design
				$sthQMP->execute($quantifMethID);
				while (my ($paramID,$paramName,$paramCode)=$sthQMP->fetchrow_array) {
					@{$quantifParamInfo{$quantifID}{$paramCode}}=($paramID,$paramName);
				}
			}
			my $quantifSoftware=($refQuantifInfo->{$quantifID}[1]->{'SOFTWARE'})? $refQuantifInfo->{$quantifID}[1]->{'SOFTWARE'}[0] : 'myProMS';
			$quantifSoftware='MaxQuant' if $quantifSoftware eq 'MQ';
			my $ratioType=($refQuantifInfo->{$quantifID}[1]->{'RATIO_TYPE'})? $refQuantifInfo->{$quantifID}[1]->{'RATIO_TYPE'}[0] : 'Ratio';
			my @quantifParams;
			if ($ratioType eq 'None') {
				@quantifParams=($refParams->{MEASURE});
			}
			else {
				@quantifParams=('RATIO'); # Must start with RATIO!
				my $pvalueCode=($ratioType=~/S\w+Ratio/ || $refQuantifInfo->{$quantifID}[1]->{'FDR_CONTROL'}[0] eq 'TRUE')? 'PVAL_ADJ' : 'PVAL';
				push @quantifParams,$pvalueCode if $quantifSoftware ne 'MaxQuant';
			}
			#push @quantifParams,'DIST_PEP_USED' if ($numPepCode eq 'DIST_PEP_USED' && $ratioType=~/S\w+Ratio/ && $view ne 'heatmap');
			my (%filteredProtIDs,%allowedProtIDs);
			foreach my $paramCode (@quantifParams) {
				$sthProtQ->execute($quantifID,$ratioPos,$quantifParamInfo{$quantifID}{$paramCode}[0]);
				my $paramKey=($paramCode=~/PVAL/)? 'P_VALUE' : $paramCode; # don't care if PVAL_ADJ or PVAL
				if ($view=~/volcano|log2/) { # {quantif}{protID} !!!
					while (my ($protID,$modResStrg,$qValue)=$sthProtQ->fetchrow_array) {
						next if ($protSelection && !$refSelectedProteins->{$protID});
						next if ($exclusion && $refExcludedProteins->{$protID});
						#my $modProtID=($modResStrg)? $protID.'-'.&decodeModificationSite($modResStrg) : $protID; # quantif of modification
						my $modProtID=$protID;
						if ($modResStrg) {
							$modProtID.='-'.$modifPrefixStrg.&decodeModificationSite($modResStrg);
							next if ($isoSelection && !$refSelectedIsoforms->{$modProtID});
						}
						$allowedProtIDs{$modProtID}=1;
						$refQuantifValues->{$quantif}{$modProtID}{$paramKey}=$qValue;
						#@{$refProteinInfo->{$protID}}=();
						if ($paramCode eq 'RATIO') { # recording min/max (!infinite) fold changes
							if ($qValue < $refParams->{MIN_RATIO}) {
								if ($qValue==$MIN_INF_RATIO_DB) {$refParams->{MINUS_INF}{$quantif}=1;} else {$refParams->{MIN_RATIO}=$qValue;}
							}
							elsif ($qValue > $refParams->{MAX_RATIO}) {
								if ($qValue==$MAX_INF_RATIO_DB) {$refParams->{PLUS_INF}{$quantif}=1;} else {$refParams->{MAX_RATIO}=$qValue;}
							}
						}
						elsif ($paramKey eq 'P_VALUE') {
							$refParams->{MIN_PVALUE}=$qValue if ($qValue && $refParams->{MIN_PVALUE} > $qValue); # ignore pval=0
						}
					}
				}
				else { # heatmap,list,explorAna: {protID}{quantif} !!!
					while (my ($protID,$modResStrg,$qValue)=$sthProtQ->fetchrow_array) {
						next if ($protSelection && !$refSelectedProteins->{$protID});
						next if ($exclusion && $refExcludedProteins->{$protID});
						#my $modProtID=($modResStrg)? $protID.'-'.&decodeModificationSite($modResStrg) : $protID; # quantif of modification
						#if ($refParams->{MOD_SELECT}) { # ajax call on selected set of (modif) proteins
						#	if ($modResStrg) {next if !$refSelectedProteins->{$modProtID};} # mod vs mod
						#	else {next if (!$refSelectedProteins->{$protID}) && !map {/^$protID-/} keys %{$refSelectedProteins};} # std quantif called with modProtIDs in refSelect (<- log2 plot)
						#}
						#else {next if ($protSelection && !$refSelectedProteins->{$protID});} # classical
						my $modProtID=$protID;
						if ($modResStrg) {
							$modProtID.='-'.$modifPrefixStrg.&decodeModificationSite($modResStrg);
							next if ($isoSelection && !$refSelectedIsoforms->{$modProtID});
						}
						if ($refParams->{STRICT_FILTER} && $refParams->{STRICT_FILTER}{$paramCode}) {
							if ($paramCode eq 'RATIO' && $qValue > $refParams->{STRICT_FILTER}{INV_RATIO} && $qValue < $refParams->{STRICT_FILTER}{RATIO}) {$filteredProtIDs{$modProtID}=1; next;}
							elsif ($paramCode eq 'P_VALUE' && $qValue > $refParams->{STRICT_FILTER}{$paramCode}) {$filteredProtIDs{$modProtID}=1; next;}
							#elsif ($qValue < $refParams->{STRICT_FILTER}{$paramCode}) {$filteredProtIDs{$modProtID}=1; next;}
						}
						$allowedProtIDs{$modProtID}=1;
						$refQuantifValues->{$modProtID}{$quantif}{$paramKey}=$qValue;
						#@{$refProteinInfo->{$protID}}=();
						if ($paramCode eq 'RATIO') { # recording min/max (!infinite) fold changes
							#if ($qValue < $minRatio && $qValue != $MIN_INF_RATIO_DB) {$minRatio=$qValue;}
							#elsif ($qValue > $maxRatio && $qValue != $MAX_INF_RATIO_DB) {$maxRatio=$qValue;}
							if ($qValue < $refParams->{MIN_RATIO}) {
								if ($qValue==$MIN_INF_RATIO_DB) {$refParams->{MINUS_INF}{$quantif}=1;} else {$refParams->{MIN_RATIO}=$qValue;}
							}
							elsif ($qValue > $refParams->{MAX_RATIO}) {
								if ($qValue==$MAX_INF_RATIO_DB) {$refParams->{PLUS_INF}{$quantif}=1;} else {$refParams->{MAX_RATIO}=$qValue;}
							}
						}
						elsif ($paramKey eq 'P_VALUE') {
							$refParams->{MIN_PVALUE}=$qValue if ($qValue && $refParams->{MIN_PVALUE} > $qValue); # ignore pval=0
						}
					}
				}
			}

			$numPepCode='RAZ_UNI_PEP' if ($ratioType eq 'None' && $numPepCode eq 'NUM_PEP_USED'); # fall back for compatibility with MaxQuant Intensity quantifs
			if ($ratioType eq 'Ratio' || $numPepCode=~/^(PEPTIDES|RAZ_UNI_PEP|UNIQUE_PEP)$/) {
				$sthProtQ0->execute($quantifID,$quantifParamInfo{$quantifID}{$numPepCode}[0]);
				if ($view=~/volcano|log2/) {
					while (my ($protID,$modResStrg,$qValue)=$sthProtQ0->fetchrow_array) {
						my $modProtID=($modResStrg)? $protID.'-'.$modifPrefixStrg.&decodeModificationSite($modResStrg) : $protID; # quantif of modification
						#next if (!$refQuantifValues->{$quantif}{$modProtID} || !$refQuantifValues->{$quantif}{$modProtID}{'RATIO'}); # in some extreme cases: prevents undef values in volcano plot
						next unless $allowedProtIDs{$modProtID};
						$refQuantifValues->{$quantif}{$modProtID}{$numPepCode}=$qValue; # quantif => protID !!!
						@{$refProteinInfo->{$protID}}=();
					}
				}
				else {
					while (my ($protID,$modResStrg,$qValue)=$sthProtQ0->fetchrow_array) {
						my $modProtID=($modResStrg)? $protID.'-'.$modifPrefixStrg.&decodeModificationSite($modResStrg) : $protID; # quantif of modification
						next if $filteredProtIDs{$modProtID};
						if ($refParams->{STRICT_FILTER} && $refParams->{STRICT_FILTER}{$numPepCode} && $qValue < $refParams->{STRICT_FILTER}{$numPepCode}) {
							$filteredProtIDs{$modProtID}=1;
							next;
						}
						next unless $allowedProtIDs{$modProtID};
						$refQuantifValues->{$modProtID}{$quantif}{$numPepCode}=$qValue; # quantif => protID !!!
						@{$refProteinInfo->{$protID}}=();
					}
				}
			}
			elsif ($numPepCode=~/(NUM|DIST)_PEP_USED/) { # S(uper/imple)Ratio
				my @usedRatios=($ratioPos);
				if ($refQuantifInfo->{$quantifID}[1]->{'RATIOS'} && $refQuantifInfo->{$quantifID}[1]->{'RATIOS'}[$ratioPos-1]=~/%/) { # 2ndary ratio of SuperRatio
					my ($testStatePos,$refStatePos)=split(/\//,$refQuantifInfo->{$quantifID}[1]->{'RATIOS'}[$ratioPos-1]);
					$testStatePos=~s/\d+%//; # keep Pos only = TARGET_POS
					$refStatePos=~s/\d+%//; # keep Pos only = TARGET_POS
					push @usedRatios,($testStatePos,$refStatePos);
				}
				my %usedProtIDs; # Required if numPep is sum of both primary ratios. now min(prim ratios=> immmediate filtering possible)
				my $numPepInDB=0;
				R_POS:foreach my $rPos (@usedRatios) {
					$sthProtQ->execute($quantifID,$rPos,$quantifParamInfo{$quantifID}{$numPepCode}[0]);
					if ($view=~/volcano|log2/) {
						while (my ($protID,$modResStrg,$qValue)=$sthProtQ->fetchrow_array) {
							$numPepInDB=1 if $rPos==$ratioPos;
							my $modProtID=($modResStrg)? $protID.'-'.$modifPrefixStrg.&decodeModificationSite($modResStrg) : $protID; # quantif of modification
							next unless $allowedProtIDs{$modProtID};
							#$refQuantifValues->{$quantif}{$modProtID}{$numPepCode}+=$qValue; # SUM!!! used DIST_PEP_USED for more stringent estimation
							$refQuantifValues->{$quantif}{$modProtID}{$numPepCode}=$qValue if (!$refQuantifValues->{$quantif}{$modProtID}{$numPepCode} || $qValue < $refQuantifValues->{$quantif}{$modProtID}{$numPepCode});
							@{$refProteinInfo->{$protID}}=();
						}
					}
					else {
						while (my ($protID,$modResStrg,$qValue)=$sthProtQ->fetchrow_array) {
							$numPepInDB=1 if $rPos==$ratioPos;
							my $modProtID=($modResStrg)? $protID.'-'.$modifPrefixStrg.&decodeModificationSite($modResStrg) : $protID; # quantif of modification
							next if $filteredProtIDs{$modProtID};
							next unless $allowedProtIDs{$modProtID};
							#$refQuantifValues->{$modProtID}{$quantif}{$numPepCode}+=$qValue;
							$refQuantifValues->{$modProtID}{$quantif}{$numPepCode}=$qValue if (!$refQuantifValues->{$modProtID}{$quantif}{$numPepCode} || $qValue < $refQuantifValues->{$modProtID}{$quantif}{$numPepCode});
							$usedProtIDs{$modProtID}=$protID;
				#@{$refProteinInfo->{$protID}}=();
#print "$quantif:$rPos:$modProtID:$numPepCode=>ref=$refQuantifValues->{$modProtID}{$quantif}{$numPepCode} ($refQuantifValues->{$modProtID}{$quantif}{RATIO})<BR>\n";
						}
					}
					last R_POS if $numPepInDB;
				}
				# Required if numPep is based on both primary ratios.
				foreach my $modProtID (keys %usedProtIDs) { # apply filter after compiling both primary ratios num pep
					if ($refParams->{STRICT_FILTER} && $refParams->{STRICT_FILTER}{$numPepCode} && $refQuantifValues->{$modProtID}{$quantif}{$numPepCode} < $refParams->{STRICT_FILTER}{$numPepCode}) {
						$filteredProtIDs{$modProtID}=1;
					}
					else {@{$refProteinInfo->{$usedProtIDs{$modProtID}}}=();}
				}

###				}
###				else { # primary ratio
###					$sthProtQ->execute($quantifID,$ratioPos,$quantifParamInfo{$quantifID}{$numPepCode}[0]);
###					if ($view=~/volcano|log2/) {
###						while (my ($protID,$modResStrg,$qValue)=$sthProtQ->fetchrow_array) {
###							#next if ($protSelection && !$refSelectedProteins->{$protID});
###							#next if ($exclusion && $refExcludedProteins->{$protID});
###							my $modProtID=($modResStrg)? $protID.'-'.$modifPrefixStrg.&decodeModificationSite($modResStrg) : $protID; # quantif of modification
###							next unless $allowedProtIDs{$modProtID};
###							$refQuantifValues->{$quantif}{$modProtID}{$numPepCode}=$qValue;
###						}
###					}
###					else {
###						while (my ($protID,$modResStrg,$qValue)=$sthProtQ->fetchrow_array) {
###							#next if ($protSelection && !$refSelectedProteins->{$protID});
###							#next if ($exclusion && $refExcludedProteins->{$protID});
###							my $modProtID=($modResStrg)? $protID.'-'.$modifPrefixStrg.&decodeModificationSite($modResStrg) : $protID; # quantif of modification
###next if $filteredProtIDs{$modProtID};
###if ($refParams->{STRICT_FILTER} && $refParams->{STRICT_FILTER}{$numPepCode} && $qValue < $refParams->{STRICT_FILTER}{$numPepCode}) {
###	$filteredProtIDs{$modProtID}=1;
###	next;
###}
###							next unless $allowedProtIDs{$modProtID};
###							$refQuantifValues->{$modProtID}{$quantif}{$numPepCode}=$qValue;
###						}
###					}
###				}
			}
			foreach my $modProtID (keys %filteredProtIDs) {
				if ($refQuantifValues->{$modProtID} && $refQuantifValues->{$modProtID}{$quantif}) {
					delete $refQuantifValues->{$modProtID}{$quantif};
					if (scalar keys %{$refQuantifValues->{$modProtID}} == 0) {
						delete $refQuantifValues->{$modProtID};
						my ($protID)=($modProtID=~/^(\d+)/);
						if ($protID eq $modProtID || !map {/^$protID-/} keys %{$refQuantifValues}) {
							delete $refProteinInfo->{$protID};
						}
					}
				}
			}

		}
		$sthQinfo->finish;
		$sthAna->finish;
		$sthQMP->finish;
		$sthProtQ0->finish;
		$sthProtQ->finish;


		##>Clean 2ndary Ratios (for Quantifs before bug fix 23/02/15: ratios stored for noSuperRatioProteins)
		foreach my $quantif (@selectedQuantifications) {
			my ($quantifID,$ratioPos)=split('_',$quantif);
			if ($refQuantifInfo->{$quantifID}[1]->{'RATIOS'} && $refQuantifInfo->{$quantifID}[1]->{'RATIOS'}[$ratioPos-1]=~/%/) {
				if ($view=~/volcano|log2/) {
					foreach my $modProtID (keys %{$refQuantifValues->{$quantif}}) {
						unless ($refQuantifValues->{$quantif}{$modProtID}{$numPepCode}) {
							delete $refQuantifValues->{$quantif}{$modProtID};
							#my ($protID)=($modProtID=~/^(\d+)/);
							#delete $refProteinInfo->{$protID} if ($protID eq $modProtID || !map {/^$protID-/} keys %{$refQuantifValues->{$quantif}});
						}
					}
				}
				else {
					foreach my $modProtID (keys %{$refQuantifValues}) {
						next unless $refQuantifValues->{$modProtID}{$quantif};
						unless ($refQuantifValues->{$modProtID}{$quantif}{$numPepCode}) {
							delete $refQuantifValues->{$modProtID}{$quantif};
							if (scalar keys %{$refQuantifValues->{$modProtID}}==0) {
								delete $refQuantifValues->{$modProtID};
								#my ($protID)=($modProtID=~/^(\d+)/);
								#if ($protID eq $modProtID || !map {/^$protID-/} keys %{$refQuantifValues}) {
								#	delete $refProteinInfo->{$protID};
#print "DEL: $protID<BR>\n";
								#}
							}
						}
					}
				}
			}
		}
	}

	####<SIN or EMPAI>####
	elsif ($selQuantifFamily=~/SIN|EMPAI/) {
		my %paramCodes;
		my $sthQMP=$dbh->prepare("SELECT ID_QUANTIF_PARAMETER,QP.CODE FROM QUANTIFICATION_PARAMETER QP,QUANTIFICATION_METHOD QM WHERE QP.ID_QUANTIFICATION_METHOD=.QM.ID_QUANTIFICATION_METHOD AND QM.CODE='$selQuantifFamily'");
		$sthQMP->execute;
		while (my ($paramID,$paramCode)=$sthQMP->fetchrow_array) {
			$paramCodes{$paramID}=$paramCode;
		}
		$sthQMP->finish;

		my $sthQinfo=$dbh->prepare("SELECT NAME,ID_ANALYSIS FROM QUANTIFICATION Q,ANA_QUANTIFICATION AQ WHERE Q.ID_QUANTIFICATION=AQ.ID_QUANTIFICATION AND Q.ID_QUANTIFICATION=? LIMIT 0,1");
		my $sthProtQ0=$dbh->prepare("SELECT PQ.ID_PROTEIN,ID_QUANTIF_PARAMETER,QUANTIF_VALUE FROM PROTEIN_QUANTIFICATION PQ WHERE ID_QUANTIFICATION=? $restrictStrg");
		my $sthGetPepNum=$dbh->prepare("SELECT ID_PROTEIN,NUM_PEP FROM ANALYSIS_PROTEIN WHERE ID_ANALYSIS=?");
		foreach my $quantifID (@selectedQuantifications) {
			($quantifID)=(split('_',$quantifID))[0];($quantifID)=(split('_',$quantifID))[0]; # <ID_QUANTIFICATION>_0 when called for export option in startExporatoryAnalysis.cgi
			print '.' if $verbose;
			$sthQinfo->execute($quantifID);
			my ($quantifName,$anaID)=$sthQinfo->fetchrow_array;
			my @itemInfo=&promsMod::getItemInfo($dbh,'QUANTIFICATION',$quantifID);
			@{$refQuantifInfo->{$quantifID}}=($quantifName,undef,undef,\@itemInfo,$anaID);
			$sthProtQ0->execute($quantifID);
			while (my ($protID,$paramID,$qValue)=$sthProtQ0->fetchrow_array) {
				next if ($protSelection && !$refSelectedProteins->{$protID});
				next if ($exclusion && $refExcludedProteins->{$protID});
				@{$refProteinInfo->{$protID}}=();
				if ($view =~ /log2/) { # {quantif}{protID} !!!
					$refQuantifValues->{$quantifID}{$protID}{$paramCodes{$paramID}}=$qValue;
				}
				else { # heatmap,list,explorAna: {protID}{quantif} !!!
					$refQuantifValues->{$protID}{$quantifID}{$paramCodes{$paramID}}=$qValue;
				}
			}
			#<Peptides
			$sthGetPepNum->execute($anaID);
			if ($view =~ /log2/) { # {quantif}{protID} !!!
				while (my ($protID,$numPep)=$sthGetPepNum->fetchrow_array) {
					next if (!$refQuantifValues->{$quantifID} || !$refQuantifValues->{$quantifID}{$protID});
					$refQuantifValues->{$quantifID}{$protID}{PEPTIDES}=$numPep;
				}
			}
			else { # heatmap,list,explorAna: {protID}{quantif} !!!
				while (my ($protID,$numPep)=$sthGetPepNum->fetchrow_array) {
					next if (!$refQuantifValues->{$protID} || !$refQuantifValues->{$protID}{$quantifID});
					$refQuantifValues->{$protID}{$quantifID}{PEPTIDES}=$numPep;
				}
			}
		}
		$sthQinfo->finish;
		$sthProtQ0->finish;
		$sthGetPepNum->finish;
	}

	####<Fetching protein info>####
	if ($getProtInfo) { # not needed for exploratory analysis & GO quantif
		print '+' if $verbose;
		my ($extraFieldStrg,$sthGN);
		if ($view=~/list|export/) {
			$extraFieldStrg=($view eq 'list')? ',MW,PROT_DES,ORGANISM' : ',IDENTIFIER,PROT_DES,ORGANISM';
			my ($GNidentID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='GN'");
			$sthGN=$dbh->prepare("SELECT MI.VALUE FROM PROTEIN P,MASTERPROT_IDENTIFIER MI WHERE P.ID_MASTER_PROTEIN=MI.ID_MASTER_PROTEIN AND MI.ID_IDENTIFIER=$GNidentID AND P.ID_PROTEIN=? ORDER BY MI.RANK");
		}
		else {$extraFieldStrg='';}
		my $sthProtInfo=$dbh->prepare("SELECT ALIAS$extraFieldStrg,PROT_LENGTH FROM PROTEIN WHERE ID_PROTEIN=?");
		my $protCount=0;
		foreach my $protID (keys %{$refProteinInfo}) {
			$protCount++;
			if ($protCount==500) {
				print '.' if $verbose;
				$protCount=0;
			}
			$sthProtInfo->execute($protID);
			@{$refProteinInfo->{$protID}}=$sthProtInfo->fetchrow_array;
			$refProteinInfo->{$protID}[0]=~s/ .*//; # Clean protein alias from badly parsed characters in case identifier conversion
			$refProteinInfo->{$protID}[0]=~s/[,;']/\./g; # Clean MaxQuant crappy contaminant identifiers
			if ($view=~/list|export/) {
				$refProteinInfo->{$protID}[1]=sprintf "%.1f",$refProteinInfo->{$protID}[1]/1000; # MW
				my @geneList;
				$sthGN->execute($protID);
				while (my ($gene)=$sthGN->fetchrow_array) {push @geneList,$gene;}
				push @{$refProteinInfo->{$protID}},\@geneList;
			}
		}
		$sthProtInfo->finish;
		$sthGN->finish if $view=~/list|export/;
	}
}


###########################################################################################
######<Matching a bioSample property value to a list of quantifications & ratioPos>######## eg. used in Exploratory Analysis
###########################################################################################
sub getQuantificationsFromPropertyValue {
	my ($dbh,$propID,$propType,$propValue,$quantifList)=@_;

	my $queryMatchStrg;
	if ($propType eq 'T') {
		$queryMatchStrg='LIKE';
		$propValue=~s/^stepValue=\d:*//;
		$propValue='%'.$propValue;
	}
	else {$queryMatchStrg='=';}

	my %quantifRatios;
	foreach my $quantData (split(':',$quantifList)) {
		my ($quantID,$targetPos)=split('_',$quantData);
		push @{$quantifRatios{$quantID}},$targetPos;
	}

	#####<Find Obs matching prop value in experiment
	###my %observations;
	###my $sthObs=$dbh->prepare("SELECT O.ID_OBSERVATION,O.TARGET_POS FROM OBSERVATION O
	###							INNER JOIN BIOSAMPLE_PROPERTY BP ON O.ID_BIOSAMPLE=BP.ID_BIOSAMPLE
	###							INNER JOIN ANALYSIS A ON O.ID_ANALYSIS=A.ID_ANALYSIS
	###							INNER JOIN SAMPLE S ON A.ID_SAMPLE=S.ID_SAMPLE
	###							WHERE S.ID_EXPERIMENT=$experimentID AND BP.ID_PROPERTY=$propID AND BP.PROPERTY_VALUE $queryMatchStrg ?
	###						");
	###							# GROUP BY BP.ID_BIOSAMPLE # Dangerous to group in case returned Obs not used in explorAna quantifs!!!
	###
	###foreach my $propValue (@propertyValues) {
	###	$sthObs->execute($propValue);
	###	while (my ($obsID,$targetPos)=$sthObs->fetchrow_array) {
	###		$observations{$obsID}=$targetPos;
	###	}
	###}
	###$sthObs->finish;

	my $sthQType=$dbh->prepare("SELECT ID_DESIGN,QM.CODE,QUANTIF_ANNOT FROM QUANTIFICATION Q,QUANTIFICATION_METHOD QM WHERE Q.ID_QUANTIFICATION_METHOD=QM.ID_QUANTIFICATION_METHOD AND ID_QUANTIFICATION=?");
	my $sthNoRatioQ=$dbh->prepare("SELECT 1 FROM OBSERVATION O
										INNER JOIN BIOSAMPLE_PROPERTY BP ON O.ID_BIOSAMPLE=BP.ID_BIOSAMPLE
										INNER JOIN ANA_QUANTIFICATION AQ ON O.ID_ANALYSIS=AQ.ID_ANALYSIS
										WHERE AQ.ID_QUANTIFICATION=? AND BP.ID_PROPERTY=$propID AND BP.PROPERTY_VALUE $queryMatchStrg ? LIMIT 0,1");

	my $sthDesignQ=$dbh->prepare("SELECT 1 FROM OBSERVATION O
										INNER JOIN BIOSAMPLE_PROPERTY BP ON O.ID_BIOSAMPLE=BP.ID_BIOSAMPLE
										WHERE O.ID_OBSERVATION=? AND BP.ID_PROPERTY=$propID AND BP.PROPERTY_VALUE $queryMatchStrg ?");

	my $sthAnaQ=$dbh->prepare("SELECT ID_ANALYSIS FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=?");
	my $sthAnaObs=$dbh->prepare("SELECT TARGET_POS,GROUP_CONCAT(ID_MODIFICATION ORDER BY ID_MODIFICATION SEPARATOR ',')
										FROM OBSERVATION O
										LEFT JOIN OBS_MODIFICATION OM ON O.ID_OBSERVATION=OM.ID_OBSERVATION
										INNER JOIN BIOSAMPLE_PROPERTY BP ON O.ID_BIOSAMPLE=BP.ID_BIOSAMPLE
										WHERE O.ID_ANALYSIS=? AND BP.ID_PROPERTY=$propID AND BP.PROPERTY_VALUE $queryMatchStrg ? GROUP BY O.ID_OBSERVATION");
	my $sthParQ=$dbh->prepare("SELECT QUANTIF_ANNOT FROM QUANTIFICATION Q,PARENT_QUANTIFICATION PQ WHERE PQ.ID_PARENT_QUANTIFICATION=Q.ID_QUANTIFICATION AND PQ.ID_QUANTIFICATION=?");

	###<Quantif type
	my %matchedQuantifs=('TEST'=>{},'REF'=>{});
	QUANTIF:foreach my $quantifID (keys %quantifRatios) {
		$sthQType->execute($quantifID);
		my ($designID,$quantifMethod,$quantifAnnot)=$sthQType->fetchrow_array;

		##<No ratio quantif
		if ($quantifMethod=~/EMPAI|SIN|MQ/) {
			$sthNoRatioQ->execute($quantifID,$propValue);
			my ($matchProp)=$sthNoRatioQ->fetchrow_array;
			if ($matchProp) {
				$matchedQuantifs{TEST}{$quantifID.'_'.$quantifRatios{$quantifID}[0]}=1; # only 1 tgPos=0
			}
		}
		else { # PROT_RATIO_PEP or TNPQ
			my (%labelingInfo,%stateInfo);
			my ($labelStrg,@labelInfo)=split('::',$quantifAnnot);
			my $labelType=($labelStrg=~/LABEL=(.+)/)? uc($1) : 'FREE';
			foreach my $infoStrg (@labelInfo) {
				my ($setting,$valueStrg)=split('=',$infoStrg);
				@{$labelingInfo{$setting}}=split(';',$valueStrg);
			}
			TARGET_POS:foreach my $targetPos (@{$quantifRatios{$quantifID}}) {

				#<Find matching test & ref states
				my $ratioData=$labelingInfo{'RATIOS'}[$targetPos-1];
				#$ratioData=~s/#//g; # remove id flags
				my $stateType;
				foreach my $stateCode (split('/',$ratioData)) { # scans both test & ref states (test is return if at least 1 match found)
					$stateType=(!$stateType)? 'TEST' : 'REF';

					##<Design quantif
					if ($designID) {
						$stateCode=~s/%.+//; # super ratio
						my $currState;
						foreach my $state (@{$labelingInfo{'STATES'}}) {
							if ($state=~/,$stateCode\Z/) {
								$currState=(split(',',$state))[1]; # numRep,replic,<id/tgPos>
								last;
							}
						}
						$currState=~s/#//g; # remove id flags
						my ($firstObsID)=($currState=~/^(\d+):/); # Replicate bioSamples should have same properties!!! (no need to loop on all replics)
						$sthDesignQ->execute($firstObsID,$propValue);
						my ($matchProp)=$sthDesignQ->fetchrow_array;
						if ($matchProp) {
							$matchedQuantifs{$stateType}{$quantifID.'_'.$targetPos}=1;
							next TARGET_POS;
						}
					}
					##<Internal labeled quantif
					else {
						my $currState=$stateCode;

						#<Observations linked to analysis
						$sthAnaQ->execute($quantifID);
						my ($anaID)=$sthAnaQ->fetchrow_array;
						$sthAnaObs->execute($anaID);

						# iTRAQ
						if ($labelType eq 'ITRAQ') { # obs targetPos=quantif targetPos
							while (my ($tgPos,$modCode)=$sthAnaObs->fetchrow_array) {
								if ($tgPos==$targetPos) {
									$matchedQuantifs{$stateType}{$quantifID.'_'.$targetPos}=1;
									next TARGET_POS;
								}
							}
						}

						elsif ($labelType eq 'SILAC') {
							my %observations;
							while (my ($tgPos,$modCode)=$sthAnaObs->fetchrow_array) {
								$modCode=0 unless $modCode;
								$observations{$modCode}=1;
							}
							next QUANTIF unless scalar keys %observations; # No obs or no BioSample match for any label channel of this analysis

							# Find matching modCode for $targetPos
							$sthParQ->execute($quantifID);
							my ($pepQuantifAnnot)=$sthParQ->fetchrow_array;
							my ($pepLabelStrg,@pepLabelInfo)=split('::',$pepQuantifAnnot);
							foreach my $infoStrg (@pepLabelInfo) {
								my ($chanNum,$chanName,$labelStrg)=split(';',$infoStrg); # '1;Light;No label;'   '2;Heavy;Label:13C(6);K'
								if ($chanNum==$currState) {
									if ($labelStrg=~/^None#/) { # unlabeled channel
										if ($observations{0}) {
											$matchedQuantifs{$stateType}{$quantifID.'_'.$targetPos}=1;
											next TARGET_POS;
										}
									}
									else {
										my @mods;
										foreach my $label (split(/\@/,$labelStrg)) {
											my ($labelModifAlias,$labelModifName,$modifRes,$searchMod,$modID)=split(/#/,$label);
											unless ($modID) { #only defined for MassChroQ XIC
												$modID=&promsMod::getModificationIDfromString($dbh,$searchMod);
											}
											push @mods,$modID if $modID;
										}
										my $modCode=join(',',sort{$a<=>$b} @mods);
										if ($observations{$modCode}) {
											$matchedQuantifs{$stateType}{$quantifID.'_'.$targetPos}=1;
											next TARGET_POS;
										}
									}
									last;
								}
							}
						}
					}

				}

			}
		}
	}

	$sthQType->finish;
	$sthNoRatioQ->finish;
	$sthDesignQ->finish;
	$sthAnaQ->finish;
	$sthAnaObs->finish;
	$sthParQ->finish;

	return (scalar %{$matchedQuantifs{TEST}})? %{$matchedQuantifs{TEST}} : %{$matchedQuantifs{REF}}; # return matching test conds if any. If none return matching ref conds

}

sub convertTreatmentToText {
	my $treatmentValue=$_[0];
	return 'treated' if (!$treatmentValue || $treatmentValue eq '%' || $treatmentValue=~/^stepValue=\d\Z/); # only step data
	my %valueData;
	foreach my $attrVal (split('::',$treatmentValue)) { # stepValue=1::quantity=16::quantUnit=gray::duration=10::durUnit=min
		my ($attr,$val)=split('=',$attrVal);
		$valueData{$attr}=$val;
	}
	my $treatmentText="$valueData{quantity} $valueData{quantUnit}" if $valueData{'quantity'};
	$treatmentText.='-' if ($valueData{'quantity'} && $valueData{'duration'});
	$treatmentText.="$valueData{duration} $valueData{durUnit}" if $valueData{'duration'};
	return $treatmentText;
}

sub fetchSiteList {
	my ($dbh,$listID,$refList,$focusModID)=@_; # Optional $focusModID: remove prefix [modID] for this modif
	my %modifInfo=();
	my $sthSite=$dbh->prepare("SELECT ID_PROTEIN,GROUP_CONCAT('[',ID_MODIFICATION,']',RESIDUE,POSITION ORDER BY ID_MODIFICATION,POSITION SEPARATOR '.') FROM CATEGORY_PROTEIN CP,MODIFICATION_SITE MS WHERE CP.ID_CATEGORY_PROTEIN=MS.ID_CATEGORY_PROTEIN AND CP.ID_CATEGORY=? GROUP BY CP.ID_CATEGORY_PROTEIN");
# Formating improvement for ambiguous pos & compatible with multi-modifID sites BUT slower than above
# Ambiguous: '[9]-1227.+1234:1/6', normal: '[9]S190.S193:/' ... [&<other modif code>]
#SELECT RES1.P_ID,GROUP_CONCAT('[',RES1.M_ID,']',RES1.POS,':',COALESCE(MS.RESIDUE,''),'/',COALESCE(MS.POSITION,'') SEPARATOR '&') AS MODCODE FROM
#	(SELECT CP.ID_CATEGORY_PROTEIN AS CP_ID,ID_PROTEIN AS P_ID,ID_MODIFICATION AS M_ID,GROUP_CONCAT(RESIDUE,POSITION ORDER BY ID_MODIFICATION,POSITION SEPARATOR '.') AS POS
#		FROM CATEGORY_PROTEIN CP,MODIFICATION_SITE MS
#        WHERE CP.ID_CATEGORY_PROTEIN=MS.ID_CATEGORY_PROTEIN AND CP.ID_CATEGORY=171 AND RESIDUE NOT REGEXP '[0-9]'
#        GROUP BY CP.ID_CATEGORY_PROTEIN,ID_PROTEIN)
#	AS RES1
#LEFT JOIN MODIFICATION_SITE MS ON RES1.CP_ID=MS.ID_CATEGORY_PROTEIN AND RES1.M_ID=MS.ID_MODIFICATION AND MS.RESIDUE REGEXP '[0-9]'
#GROUP BY RES1.CP_ID
	$sthSite->execute($listID);
	while (my ($protID,$modCode)=$sthSite->fetchrow_array) {
		if ($modCode=~/\]\d/) { # (occurence) => ambiguous
			$modCode=~s/(\[\d+\]\d+)\.(.+)/$2\.$1/; # Occurence is likely first when ORDER BY POSITION => send at the end. WARNING: Won't work for multi-modif site
		}
		$modCode=~s/(.)\[\d+\]/$1/g; # remove all [modID] except first. WARNING: Won't be correct for multi-modif site
		$modCode=~s/^\[$focusModID\]// if $focusModID; # [\d+] will is kept if different from [$noPrefixModID]!!! wanted to distinguish from another modif. WARNING: Won't be correct for multi-modif site
		$refList->{"$protID-$modCode"}=1;
	}
	$sthSite->finish;
}

#--------- For insertion of quantified modification sites by run(XIC/SWATH)ProtQuantification.pl --------------
sub insertModifiedResidues {
	my ($protID,$refModResList,$refModResiduesID,$dbh,$sthInsModRes)=@_;
	my @newList;
	if ($refModResList->[0]=~/~/) {&encodeModifPosAmbiguity($refModResList,\@newList)} else {@newList=@{$refModResList};} # $refModResList must not be overwritten!
#print DEBUG "\tMODRES_ID='@newList'\n"; # unless $refModResiduesID->{$modRes};
	foreach my $modRes (@newList) {
		if (!$refModResiduesID->{$protID} || !$refModResiduesID->{$protID}{$modRes}) {
			my ($res,$pos)=($modRes=~/^(.)(.+)/); # Y180 -> (Y,180)
			$sthInsModRes->execute($res,$pos);
			$refModResiduesID->{$protID}{$modRes}=$dbh->last_insert_id(undef,undef,'MODIFIED_RESIDUE','ID_MODIF_RES');
		}
	}
}
sub insertProteinModResidues {
	my ($refModResList,$sthInsProtRes,$refModResiduesID,$protQuantifID)=@_;
	my @newList;
	if ($refModResList->[0]=~/~/) {&encodeModifPosAmbiguity($refModResList,\@newList)} else {@newList=@{$refModResList};} # $refModResList must not be overwritten!
#print DEBUG "\tLIST='@newList'\n"; # unless $refModResiduesID->{$modRes};
	foreach my $modRes (@newList) {
		$sthInsProtRes->execute($refModResiduesID->{$modRes},$protQuantifID);
	}
}
sub encodeModifPosAmbiguity { # (PP 2017/06/22)
	# ex1 (res~res): 123~135:2/3@<rt> --> (-123,+135,@<rt>,23) WARNING: in x/y x,y must be between [1..9]
	# ex2 (Nter~res): n0~7:2/3 --> (=0,+7,23)  "-" becomes "=" to indicate N-term (n0: 0 indicates Protein N-term)
	# ex3 (res~Cter): 251~c0:2/3 --> (-251,*0,23)  "+" becomes "*" to indicate C-term (c0: 0 indicates Protein C-term)
	# ex4 (Nter~Cter): n78~c85:2/3 --> (=78,*85,23)
	my ($refModResList,$refCodedModResList)=@_;
	my @modRes=split('[~:@/]',$refModResList->[0]);
	#@{$refCodedModResList}=('-'.$modRes[0],'+'.$modRes[1]);
	#push @{$refCodedModResList},'@'.$modRes[2] if $modRes[4]; # has RT data (optional)
	#push @{$refCodedModResList},$modRes[-2].$modRes[-1];
	# Starting pos
	if ($modRes[0]=~/^n(\d+)/) {@{$refCodedModResList}=('='.$1);} # N-ter
	else {@{$refCodedModResList}=('-'.$modRes[0]);} # res
	# Ending pos
	if ($modRes[1]=~/^c(\d+)/) {push @{$refCodedModResList},'*'.$1;} # C-ter
	else {push @{$refCodedModResList},'+'.$modRes[1];} # res
	# RT
	push @{$refCodedModResList},'@'.$modRes[2] if $modRes[4]; # has RT data (optional)
	# Matched & Available pos
	push @{$refCodedModResList},$modRes[-2].$modRes[-1];
}
#-----------------------------------------------------------------------------------

#sub convertSitesForDisplay { # optional: refParams->{SEL_MOD_ID=>[modifID] OR MODIF_FORMAT=>(name|code|code_html)}, refParams->{PROT_LIST} to also fetch protein list
#	my ($dbh,$refList,$selModifID)=@_;
#	my %modifInfo=();
#	my $sthSite=$dbh->prepare("SELECT ID_PROTEIN,GROUP_CONCAT('[',ID_MODIFICATION,']',RESIDUE,POSITION ORDER BY ID_MODIFICATION,POSITION SEPARATOR '.') FROM CATEGORY_PROTEIN CP,MODIFICATION_SITE MS WHERE CP.ID_CATEGORY_PROTEIN=MS.ID_CATEGORY_PROTEIN AND CP.ID_CATEGORY=? GROUP BY CP.ID_CATEGORY_PROTEIN");
#	$sthSite->execute($listID);
#	while (my ($protID,$modCode)=$sthSite->fetchrow_array) {
#		if ($modCode=~/\]\d/) { # (occurence) => ambiguous
#			$modCode=~s/(\[\d+\]\d+)\.(.+)/$2\.$1/; # Occurence is likely first when ORDER BY POSITION => send at the end. WARNING: Won't work for multi-modif site
#		}
#		if ($refParams->{SEL_MOD_ID} || !$refParams->{MODIF_FORMAT}) {
#			next if ($refParams->{SEL_MOD_ID} && $modCode !~ /^\[$refParams->{SEL_MOD_ID}\]/);
#			$modCode=~s/\[\d+\]//g;
#			$modCode=&decodeModificationSite($modCode);
#		}
#		else {
#			my $modifFormat=$refParams->{MODIF_FORMAT} || 'name';
#			my ($modID)=($modCode=~/\[(\d+)\]/);
#			$modCode=~s/\[\d+\]//g;
#			if ($modifFormat eq 'name') {
#				my ($psi,$interim)=$dbh->selectrow_array("SELECT PSI_MS_NAME,INTERIM_NAME FROM MODIFICATION WHERE ID_MODIFICATION=$modID");
#				$modifInfo{$modID}=($psi || $interim).':';
#			}
#			else {
#				unless ($modifInfo{$modID}) { # 1-letter code (with color if code_html)
#					my ($code,$color)=$dbh->selectrow_array("SELECT DISPLAY_CODE,DISPLAY_COLOR FROM MODIFICATION WHERE ID_MODIFICATION=$modID");
#					$modifInfo{$modID}=($modifFormat eq 'code_html')? '<FONT color="#'.$color.'">'.$code.'</FONT>:' : $code.':';
#				}
#			}
#			$modCode=$modifInfo{$modID}.&decodeModificationSite($modCode);
#		}
#		$refList->{"$protID-$modCode"}=1;
#		$refParams->{PROT_LIST}{$protID}=1 if $refParams->{PROT_LIST};
#	}
#	$sthSite->finish;
#}

sub decodeModificationSite { # modCode[,any text prefix eg; name of modification]
	# Former decodeModifPosAmbiguity
	my $modCode=$_[0];
	my $prefix=$_[1] || '';

	#<Handle [modID] if any
	my ($modID)=($modCode=~/\[(\d+)\]/); # in case [modID] is/are left
	$modCode=~s/\[\d+\]//g; # remove all to allow proper parsing...
	$prefix.='['.$modID.']' if $modID; # ... and add 1st again after parsing

	#<Parse cleaned modCode
	my ($begStrg,$endStrg,$rtStrg,$siteStrg)=('','','','');
	if ($modCode=~/[\-=\+\*]/) { # Ambiguity format WARNING: != classical terminal modif codes!!!
							  # - begins with normal res, + ends with normal res
							  # = begins with (prot/pep)N-term, * ends with (prot/pep)C-term
							  # n/c (lower case) non-ambiguous (prot/pep)(N/C)-term
		foreach my $resPos (split(/\./,$modCode)) {
			my ($res,$pos)=($resPos=~/^(.)(.+)/);
			if ($res eq '=') {$begStrg=($pos<=1)? 'ProtNt' : 'PepNt'.$pos;} # Nter beg
			elsif ($res eq '-') {$begStrg=$pos;} # res beg
			elsif ($res eq '*') {$endStrg=($pos==0)? 'ProtCt' : 'PepCt'.$pos;} # Cter end
			elsif ($res eq '+') {$endStrg=$pos;} # res end
			elsif ($res eq '@') {$rtStrg=$resPos;} # rt
			elsif ($res=~/^\d/) {$siteStrg=$res.'/'.$pos;} # sites
		}
		return $prefix.$begStrg.'~'.$endStrg.$rtStrg.':'.$siteStrg;
	}
	else {
		$modCode=~s/n(0|1)(\.|\Z)/ProtNt$2/;
		$modCode=~s/n(\d+)/PepNt$1/;
		$modCode=~s/c0(\.|\Z)/ProtCt$2/;
		$modCode=~s/c(\d+)/PepCt$1/;
		return $prefix.$modCode;
	}
}

#sub displayModificationSite { # convert starting [modifID] into (name|code|code_html):
#	my ($modProtID,$refModifInfo)=@_;
#	$modProtID=~s/\[(\d+)\]/$refModifInfo->{$1}:/;
#}


1;

####>Revision history
# 1.3.4 view=export handled in &fetchQuantificationData (PP 26/09/18)
# 1.3.3 Added subs for insertion of quantified modification sites by run(XIC/SWATH)ProtQuantification.pl (PP 19/07/18)
# 1.3.2 Clean protein alias from badly parsed characters in case identifier conversion (PP 09/07/18)
# 1.3.1 [Fix] bugs in &getQuantificationsFromPropertyValue for MaxQuant quanti (PP 20/06/18)
# 1.3.0 Deletes peptide quantification file $promsPath{quantification}/project_$projectID/quanti_$quantifID/peptide_quantification(_$targetPos).txt, not from DB (PP 30/05/18)
# 1.2.4 &getQuantificationsFromPropertyValue now scans both test and reference states. Returns test if 1+ matches; if not returns reference (PP 30/05/18)
# 1.2.3 Added 'clusters' parameter for MSstats in &writeQuantifParameterFiles (PP 24/01/18)
# 1.2.2 Fix &writeQuantifParameterFiles for log algo (PP 04/12/17)
# 1.2.1 New normalization option naming & changes in &writeQuantifParameterFiles (PP 27/11/17)
# 1.2.0 Renamed/updated &decodeModifPosAmbiguity to &decodeModificationSite, changes in fetchQuantificationData & added &fetchSiteList (PP 16/11/17)
# 1.1.9c Add &getXicSoftware (GA 08/11/17)
# 1.1.9b Minor modification to prevent error in export emPAI quantifications (GA 03/11/17)
# 1.1.9 Change query to recover analysis list from &fetchQuantificationData due to limitation in sql string (SL 10/07/17)
# 1.1.8 Extends &decodeModifPosAmbiguity for Prot/Pep N/C-term (PP 22/06/17)
# 1.1.7 Bug fix in fetchQuantificationData when isoform selection occurs (PP 22/02/17)
# 1.1.6 Minor modification on emPAI measures + add SIN in getProteinQuantifFamilies (GA 01/02/17)
# 1.1.5 Handles MaxQuant and SWATH quantifications (PP 23/01/17)
# 1.1.4 &fetchQuantificationData also returns min. p-value for RATIO quantif (PP 24/10/16)
# 1.1.3 Minor modification of &getProteinQuantifFamilies (GA 01/08/16)
# 1.1.2 Minor bug fix in &writeQuantifParameterFiles (PP 01/08/16)
# 1.1.1 Added MSstats parameters in &writeQuantifParameterFiles (PP 29/07/16)
# 1.1.0 Bug fix in OPTIMAL quantif name in &getDistinctQuantifNames (PP 24/03/16)
# 1.0.9 Add MaxQuant and emPAI quanti (GA 01/03/16)
# 1.0.8 Fix uninitialized $quantifAnnot in &getDistinctQuantifNames if internal quantif (PP 18/02/16)
# 1.0.7 Added &getProteinQuantifFamilies (PP 27/01/16)
# 1.0.6 Speed optimization & improvement of quantif optimal name in &getDistinctQuantifNames (PP 03/11/15)
# 1.0.5 &fetchQuantificationData minRatio & maxRatio no longer stay between ]$MIN_INF_RATIO & $MAX_INF_RATIO[ (PP 15/10/15)
# 1.0.4 modify split for quantifAnnot in &extractQuantificationParameters and add condition in &getDistinctQuantifNames (SL 07/09/15)
# 1.0.3 &fetchQuantificationData compatible with modification quantification (PP 04/09/15)
# 1.0.2 Added deleteQuantification, writeQuantifParameterFiles, getQuantificationParameters & getQuantifNormalizationName functions<BR>& SimpleRatio management (PP 12/11/14)
# 1.0.1 Added &getQuantificationsFromPropertyValue and &convertTreatmentToText used for exploratory analyses (PP 21/08/14)
# 1.0.0 New pm file for quantification management (PP 10/07/14)