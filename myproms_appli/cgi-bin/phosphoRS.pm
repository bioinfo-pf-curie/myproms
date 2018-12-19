###############################################################################
# phosphoRS.pm               1.3.1                                            #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                     #
# Contact: myproms@curie.fr                                                   #
# Class used to determinate phosphorylation sites on peptides using PhosphoRS #
###############################################################################
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

package phosphoRS;
use XML::Simple;
use POSIX qw(strftime);
use promsConfig;
use promsMod;
use strict;
use Storable 'dclone';
use File::Copy;

$| = 1;


sub new{
    # Class constructor
    # Usage:
    # my $phosphoRS = new phosphoRS(
    #						File => "<full path MS data file name >",
    #						ActivationType => ["CID", "ETD" or "HCD"],
    #						MassTolerance => x.x (default is 0.5),
    #						QueryFilter => \@queryNumList (optional, only these queries will be processed)
    #   				 )

    my ($class,%parameters) = @_;

    # Checking arguments #
    die "Missing argument: File" unless $parameters{'File'};
    unless ($parameters{'ActivationType'}){
		warn "Unspecified activation type (set to CID)";
		$parameters{'ActivationType'} = 'CID';
    }
    unless ($parameters{'MassTolerance'}){
		warn "Unspecified mass tolerance (set to 0.5)";
		$parameters{'MassTolerance'} = 0.5;
    }

    # Setting object attributes #
    my $this = {
		analysisID => $parameters{'AnaID'},
		fullJobDir => $parameters{'fullJobDir'}, # Full path to job directory
		dataFile => $parameters{'File'}, # name of the MS data input file
        inputRef => undef, # hash of input data, exportable by XMLSimple
        inputFile => undef, # file containing input data in XML format
        outputFile => undef, # file containing output data in XML format
        outputRef => undef, # hash converted data of XML output data
        isoforms => undef, # hash of output data in a format that can be efficiently queried by myProMS
		isoformMatch => undef # hash for retrieving isoform groups
    };

    $this->{dataFileFormat} = ($parameters{FileFormat})? $parameters{FileFormat} : 'MASCOT.DAT'; # could be PARAGON.XML or *.PDM

    bless($this,$class);

    if (uc($this->{dataFileFormat}) =~ /^MASCOT|\.PDM$/){
		$this->_addSpectraFromMascotFile(%parameters);
    }
	elsif (uc($this->{dataFileFormat}) eq 'PARAGON.XML') {
		$this->_addSpectraFromParagonFile(%parameters);
    }
	else {
		die "Unknown or unhandled file format";
    }

    return $this;
}



sub startAnalysis {
    # Method for starting the PhosphoRS analysis.
    # Use this method only if spectra were added.
    # Returns the exit code of the PhosphoRS software (0 = success, else be careful)

    my $this = shift;
    die "No input spectrum for analysis" unless defined $this->{'inputRef'};
    die "No phosphorylation in variable modification parameters for this file" unless $this->{'inputRef'}{'Phosphorylation'}{'Symbol'};
    my $inputFile = $this->_writeInputFile();
	my $fullJobDir=$this->{'fullJobDir'};
	my $outputFile = "$fullJobDir/output_data.xml";
    my %promsPath = &promsConfig::getServerInfo('no_user');
	my %cluster=&promsConfig::getClusterInfo;
	my $javaCommand=($cluster{'on'})? "$cluster{path}{java}/java" :  ($promsPath{java})? "$promsPath{java}/java" : 'java';
	my $outPRS = `$javaCommand -jar $promsPath{bin}/phosphoRS.jar $inputFile $outputFile 2> $fullJobDir/PRSerrors.txt`;
	my $prsError=`head -5 $fullJobDir/PRSerrors.txt` if -e "$fullJobDir/PRSerrors.txt";
	if ($prsError) {
		chomp($prsError);
		$prsError=~s/\n/<BR>\n/g;
		die "An unexpected error has occured:<BR>\n$prsError";
	}
	my ($code)=($outPRS =~ /Code:(\d+)/);

    if (-e $outputFile){
		$this->_storeResults($outputFile);
		
    }
    if ($code && ($code == 100 or $code == 200)) {
		my $xmlValidation = `xmllint --noout $inputFile --schema $promsPath{bin}/phosphoRS.xsd`;
		warn $xmlValidation;
    }

    return $code;
}

sub addSpectrumInfo{

    my ($this, %args) = @_;
    my $ref = $this->_getInputRef;

    my $peaks = $args{Peaks} or die "Missing peaks";
    my $queryNum = $args{Query} or die "Missing query number";
    my $charge = $args{Charge} or die "Missing charge value";
    my $activationType = ($args{Activation})? $args{Activation} : 'CID';

    if(defined $ref->{'Spectra'}{'Spectrum'}[$queryNum-1]){
		$ref->{'Spectra'}{'Spectrum'}[$queryNum-1]{'Peaks'}{'content'} = $peaks;
		$ref->{'Spectra'}{'Spectrum'}[$queryNum-1]{'ID'} = $queryNum;
		$ref->{'Spectra'}{'Spectrum'}[$queryNum-1]{'PrecursorCharge'} = $charge;
		$ref->{'Spectra'}{'Spectrum'}[$queryNum-1]{'ActivationTypes'} = $activationType;
		$this->_setInputRef($ref);
		return 1;
    }
	else {return 0;}
}

sub getIsoforms{
    # This method returns an array reference of all isoforms corresponding to the provided query number and rank.
    # Isoforms are array references where 1st element is the calculated probability and
    # the second element is an array reference of all phosphopositions.
    # Data is processed in a new hash structure first time this method is called, in order to be accessed efficiently the following times.

    my ($this, $queryNum, $peptideNum) = @_;
    my $outRef = $this->{'outputRef'} or die "Get isoforms before analysis";
    unless(defined $this->{'isoforms'}){ # result structure is built only first time this function is called
        my %isoforms;
        foreach my $spectrumRef (@{$outRef->{'Spectra'}{'Spectrum'}}){
			foreach my $peptideRef (@{$spectrumRef->{'Peptides'}{'Peptide'}}){
				foreach my $isoformRef ( @{$peptideRef->{'Isoforms'}{'Isoform'}}){
					my @positions;
					foreach my $phosphoSiteRef (@{$isoformRef->{'PhosphoSites'}{'PhosphoSite'}}){
						push @positions, $phosphoSiteRef->{'SeqPos'};
					}
					push @{$isoforms{$spectrumRef->{'ID'}}{$peptideRef->{'ID'}}}, [ $isoformRef->{'PepProb'}, \@positions ];
				}
			}
        }
        $this->{'isoforms'} = \%isoforms;
    }
    if(exists $this->{'isoforms'}{$queryNum}{$peptideNum}){
		return $this->{'isoforms'}{$queryNum}{$peptideNum};
#    } elsif ( my $altPepNum = $this->_getMatchedPhosphoPeptide($queryNum, $peptideNum)){
#	return $this->{'isoforms'}{$queryNum}{$altPepNum};
    }
	else {
		return undef;
    }
}


#------------------#
# Static functions #
#------------------#
sub getIsoformsFromFile{
    # Static function to parse XML outfile to get isoforms from 1 spectrum only
    # To get all isoforms from all spectrums, use &getIsoforms method with a phosphoRS ref object

    my ($fileName,$queryNum,$peptideNum) = @_;

    die "Missing query number" unless $queryNum;
    die "Missing rank" unless $peptideNum;

    open XML, $fileName or die "Cannot open $fileName: $!";
    my ($spectrumOK,$peptideOK) = (0,0);
    my $XMLstring; # will contain the XML part of interest
    while(<XML>){
	if(/<Spectrum ID="$queryNum">/){
	    $spectrumOK = 1;
	} elsif ($spectrumOK and /<Peptide ID="$peptideNum">/){
	    $peptideOK = 1;
	    $XMLstring = $_;
	} elsif ($spectrumOK and $peptideOK and $_ !~ /<\/Peptide>/){
	    $XMLstring .= $_;
	} elsif ($spectrumOK and $peptideOK and /<\/Peptide>/) {
	    $XMLstring .= $_;
	    last;
	}
    }
    close XML;

    my $xmlRef = XMLin($XMLstring, ForceArray => ['Peptide', 'Isoform', 'PhosphoSite']);

    my @isoforms;
    foreach my $isoformRef ( @{$xmlRef->{'Isoforms'}{'Isoform'}}){
	my @positions;
	foreach my $phosphoSiteRef (@{$isoformRef->{'PhosphoSites'}{'PhosphoSite'}}){
	    push @positions, $phosphoSiteRef->{'SeqPos'};
	}
	push @isoforms, [ $isoformRef->{'PepProb'}, $isoformRef->{'PepScore'}, \@positions ];
    }

    return \@isoforms;
}

sub cleanResults{
    my ($outputFile, $queryNumsRef) = @_;

    my $xmlRef = XMLin($outputFile, ForceArray => ['Spectrum', 'Peptide', 'Isoform', 'PhosphoSite']);
    my $newRef = dclone($xmlRef);
    $newRef->{'Spectra'}{'Spectrum'} = ();

    QUERY:foreach my $queryNum (@{$queryNumsRef}){
		next unless $queryNum;
		foreach my $spectrum (@{$xmlRef->{'Spectra'}{'Spectrum'}}){
			if ($spectrum->{'ID'} eq $queryNum) {
				push @{$newRef->{'Spectra'}{'Spectrum'}}, $spectrum;
				next QUERY;
			}
		}
    }

    my $xmlOut = XMLout($newRef, RootName => 'PhosphoRS_Results');

    open OUT, ">$outputFile" or die $!;
    print OUT $xmlOut;
    close OUT;
}

sub printIcon{
    # Static function to print a PRS notification icon for peptide lists
    # Argument: PRS database string -> PRS=status;score;positions

    my ($string,$refOptions) = @_;
	$refOptions={} unless $refOptions;
	$refOptions->{'format'}='print' unless $refOptions->{'format'};
	$refOptions->{'text'}='&nbsp;PRS&nbsp;' unless $refOptions->{'text'};
	$refOptions->{'posShift'}=0 unless $refOptions->{'posShift'};
	unless ($string) {
		if ($refOptions->{'format'}=~/^p/) {print ''; return;}
		else {return '';}
	}
	my @seq=($refOptions->{'pepSeq'})? split(//,$refOptions->{'pepSeq'}) : ();

    my ($status,$proba,$positions) = split /;/, $string;
    $proba = 1*(sprintf("%.1f", $proba));
    my ($fontColor,$backgroundColor,$popupText);
    if($status == 3){ # PRS has confirmed best Mascot isoform
	    $fontColor = 'black';
	    $backgroundColor = '#00ff2a';
	    $popupText = "Position(s) confirmed by PhosphoRS<BR>Probability=$proba%";
    }
	elsif ($status == 2) { # Isoform was changed for PRS
	    $fontColor = 'black';
	    $backgroundColor = 'yellow'; # 'red';
		$positions =~ s/[^\d\.]//g;
	    $positions =~ s/\./, /g;
	    #$positions =~ s/\[(\w+)\]([^\[,]+)/$2 \($1\) /g;
		$positions =~ s/(\d+)/$seq[$1-1].$1/eg if $refOptions->{'pepSeq'};
		$positions =~ s/(\d+)/$1+$refOptions->{'posShift'}/eg if $refOptions->{'posShift'};
		#$positions =~ s/(\d+)/$seq[$1-1].($1+$refOptions->{'posShift'})/eg; # both at once also possible
	    $popupText = "Position(s) changed by PhosphoRS (Previously: $positions)<BR>Probability=$proba%";
    }
	elsif ($status == 1) { # PRS != Mascot positions, unchanged because < threshold
	    $fontColor = 'black';
	    $backgroundColor = 'red'; # 'yellow';
		$positions =~ s/[^\d\.]//g;
	    $positions =~ s/\./, /g;
		$positions =~ s/(\d+)/$seq[$1-1].$1/eg if $refOptions->{'pepSeq'};
		$positions =~ s/(\d+)/($1+$refOptions->{'posShift'})/eg if $refOptions->{'posShift'};
	    $popupText = "Position(s) not matched by PhosphoRS ($positions)<BR>Probability=$proba% (below threshold)";
    }
	elsif ($status == 0) { # PRS = Mascot positions but < threshold
	    $fontColor = 'black';
	    $backgroundColor = 'orange'; # 'grey';
	    $popupText = "Position(s) matched by PhosphoRS<BR>Probability=$proba% (below threshold)";
    }
	elsif ($status == 4) {
	    $fontColor = 'grey';
	    $backgroundColor = 'black';
	    $popupText = "No isoform evidence found by PhosphoRS";
    }
	my $prsString="<FONT style=\"background-color:$backgroundColor;color:$fontColor;font-size:10px\" onmouseover=\"popup('$popupText');\" onmouseout=\"popout();\">$refOptions->{text}</FONT>";
	if ($refOptions->{'format'}=~/^p/) {print $prsString;}
	else {return $prsString;}
}

#---------------------#
# Private subroutines #
#---------------------#
sub _addSpectraFromMascotFile{
    # Fill inputRef with MS data file

    # arguments :
    # MassTolerance => x.x
    # File => ".dat filename"
    # QueryFilter => \@queryList (optional)

    my ($this,%parameters) = @_;

    my $ref = $this->_getInputRef();

    my $activationType = $parameters{'ActivationType'};
    $ref->{'MassTolerance'}{'Value'} = $parameters{'MassTolerance'};

    my %queryFilter;
    my $thereIsFilter = 0;
    if($parameters{'QueryFilter'}){
        $thereIsFilter = 1;
        foreach my $query (@{$parameters{'QueryFilter'}}){
            $queryFilter{$query} = 1;
        }
    }

    my $onPeptides = 0;
    my $onQuery = 0;
    my %modSymbol;
    my %symbolGroup;
    $symbolGroup{0} = 0; # 0 = no PTM
    open(FILE, $this->{'dataFile'}) or die "Unable to open file: $!";
    while(<FILE>){
        chomp;
        if(/^IT_MODS=(.+)$/){
            my @mods = split /\,/, $1;
            my $i =0;
            foreach my $modFullName (@mods){
                $i = _incrementSymbol($i);
                my ($modName,$aa) = ($modFullName =~ /^(.+)\s\((.+)\)\s?$/);
                $modName=~s/://g; # ':' is separator in modificationInfo value string
				if($modSymbol{$modName}){
					$symbolGroup{$i} = $modSymbol{$modName};
					my @modifValue = _getModifValue($ref,$modSymbol{$modName});
					$modifValue[-1] .= $aa;
					_setModifValue($ref,$modSymbol{$modName},@modifValue);
				}
				else {
					$modSymbol{$modName} = $i;
					$symbolGroup{$i} = $i;
					push @{$ref->{'ModificationInfos'}{'ModificationInfo'}}, {};
					$ref->{'ModificationInfos'}{'ModificationInfo'}[-1]{Symbol} = $i;
					$ref->{'ModificationInfos'}{'ModificationInfo'}[-1]{Value} = "$i:$modName:$modName:null:null:0:$aa";
					$ref->{'Phosphorylation'}{'Symbol'} = $i if $modName eq 'Phospho';
				}
            }
        }
		elsif (/^delta(\d+)=([^\,]+)\,.+$/){
            my $symbol = $symbolGroup{$1};
            my $deltaMass = $2;
            _setModifValue($ref,$symbol,undef,undef,undef,$deltaMass,undef,undef,undef);
        }
		elsif(/^NeutralLoss(\d+)=(.+)$/){
            my $symbol = $symbolGroup{$1};
            my $massLoss = $2;
			my @modifValue = _getModifValue($ref,$symbol);
            my $lossName = $modifValue[2] . 'Loss';
            _setModifValue($ref,$symbol,undef,undef,undef,undef,$lossName,$massLoss,undef);
        }
		elsif (/name=\"peptides\"/){
            $onPeptides = 1;
        }
		elsif ($onPeptides && /--gc0p4Jq0M2Yt08jU534c0p/){
            $onPeptides = 0;
        }
		elsif (/name=\"query(\d+)\"/){
            next if ($thereIsFilter and !$queryFilter{$1});
            $onQuery = $1;
        }
		elsif ($onQuery && /^charge=(\d+)\+/){
			next unless $ref->{'Spectra'}{'Spectrum'}[$onQuery-1]{'IdentifiedPhosphorPeptides'};
            $ref->{'Spectra'}{'Spectrum'}[$onQuery-1]{'ID'} = $onQuery;
            $ref->{'Spectra'}{'Spectrum'}[$onQuery-1]{'PrecursorCharge'} = $1;
            $ref->{'Spectra'}{'Spectrum'}[$onQuery-1]{'ActivationTypes'} = $activationType;
        }
		elsif($onQuery && /^Ions1=(.+)$/){
			next unless $ref->{'Spectra'}{'Spectrum'}[$onQuery-1]{'IdentifiedPhosphorPeptides'};
            $ref->{'Spectra'}{'Spectrum'}[$onQuery-1]{'Peaks'}{'content'} = $1;
        }
		elsif($onQuery && /--gc0p4Jq0M2Yt08jU534c0p/){
            $onQuery = 0;
        }
		elsif ($onPeptides && /^q(\d+)_p(\d+)=(.+);/){
            next if ($thereIsFilter and !$queryFilter{$1});
            my $query = $1;
            my $peptide = $2;
            my @values = split(/\,/, $3);
            my $sequence = $values[4];
            my $modifString = $values[6];
			my @modifList = split(//, $modifString);
			@modifList = map { $symbolGroup{$_} } @modifList; # replacing mod symbols by their group symbol
			$modifString = join('', @modifList) or die "$values[6]";
            substr($modifString, 1, 0) = '.';
            substr($modifString, -1, 0) = '.';
			push @{$ref->{'Spectra'}{'Spectrum'}[$query-1]{'IdentifiedPhosphorPeptides'}{'Peptide'}}, { ID => $peptide, Sequence => $sequence, ModificationInfo => $modifString };
        }
    }
    close(FILE);

    $this->_setInputRef($ref);
}

sub _addSpectraFromParagonFile{
    # Fill inputRef with MS data Paragon file
    # arguments :
    # MassTolerance => x.x
    # File => ".dat filename"

    my ($this,%parameters) = @_;

    my $activationType = $parameters{'ActivationType'};

    my $ref = $this->_getInputRef();

    $ref->{'MassTolerance'}{'Value'} = $parameters{'MassTolerance'};

    my $handler = ParagonSAXHandler->new($ref,$activationType);
    my $xmlparser = XML::SAX::ParserFactory->parser(Handler => $handler );
    $xmlparser->parse_uri($this->{dataFile});

    $this->_setInputRef($handler->{ref});

    package ParagonSAXHandler;
    use base qw(XML::SAX::Base);{

	our @el = (); # list of atoms with their symbol and isotope mass and probability
	our $currentMss = '';
	our $currentPry = '';
	our $currentSym = '';

	our @mod = ();
	our $currentTgt = '';
	our $currentNme = '';
	our $currentFma = '';
	our $currentRpF = '';
	our $currentNLF = '';

	our @phosphoSpectra = ();
	our $currentSpectrum;
	our $currentSpectrumID;
	our @currentPhosphoMatches = ();
	our $currentMatch;
	our $currentMatchID;
	our @currentMods;
	our @currentSubs;
	our $currentPeakData = '';

	our %interestVMods;

	sub new {
	    my ($class,$ref,$activationType) = @_;
	    return bless({ref => $ref, activationType => $activationType},$class);
	}

	sub start_document{
	    $currentSpectrumID = 0;
	}

	sub end_document{
	    # constructing the entire XML input-PRS-file reference
	    my $this = shift;

	    my %aaCode=('Alanine'=>'A',
			'Arginine'=>'R',
			'Asparagine'=>'N',
			'AsnOrAsp'=>'B',
			'Aspartic Acid'=>'D',
			'Cysteine'=>'C',
			'Glutamine'=>'Q',
			'GlnOrGlu'=>'Z',
			'Glutamic Acid'=>'E',
			'Glycine'=>'G',
			'Histidine'=>'H',
			'Isoleucine'=>'I',
			'Leucine'=>'L',
			'Lysine'=>'K',
			'Methionine'=>'M',
			'Phenylalanine'=>'F',
			'Proline'=>'P',
			'Serine'=>'S',
			'Threonine'=>'T',
			'Tryptophan'=>'W',
			'Tyrosine'=>'Y',
			'Valine'=>'V',
			'Selenocysteine'=>'U',
			);

	    # Atom masses
	    my %atomMass;
	    while (my $el = shift @el) {
		my $avgMass = 0;
		while (my $iso = shift @{$el->{Iso}}) {
		    $avgMass += $iso->{Mss} * $iso->{Pry};
		}
		$atomMass{$el->{Sym}} = $avgMass;
	    }
	    $this->{atomMass} = \%atomMass;

	    # ModificationInfos
	    my $modSymbol = 0;
	    my %modSymbol;
	    while (my $mod = shift @mod) {
		next unless (exists $interestVMods{$mod->{Nme}} && $interestVMods{$mod->{Nme}} == 1);
		$interestVMods{$mod->{Nme}} = 2; # avoid duplicates

		$modSymbol = _incrementSymbol($modSymbol);
		$modSymbol{$mod->{Nme}} = $modSymbol;

		$mod->{Nme} =~ s/://g; # ':' is separator in modificationInfo value string

		# Mass
		my $modMass = $this->_computeMolMass($mod->{Fma});
		$modMass -= $this->_computeMolMass($mod->{RpF}) if $mod->{RpF}; # atom(s) lost during bound

		# Neutral loss
		if ($mod->{Nme} eq 'Phospho') {
		    $mod->{NLF} = 'H3PO4';
		    # because the phospho neutralLoss is only indicated for Phospho(Ser,Thr) modification,
		    # and despite this specification, this is the modification called 'Phospho' which is used for phosphorylation located on Ser and Thr residues
		}

		my($neutralLossName,$neutralLossMass);
		if ($mod->{NLF}) {
		    $neutralLossName = $mod->{Nme} . 'Loss';
		    $neutralLossMass = $this->_computeMolMass($mod->{NLF});
		}
		else {
		    $neutralLossName = 'null';
		    $neutralLossMass = 0;
		}

		# Targets
		my $targets = '';
		while (my $tgt = shift @{$mod->{Tgt}}) {
		    $targets .= $aaCode{$tgt};
		}

		my $modValue = join ':', ($modSymbol,$mod->{Nme},$mod->{Nme},$modMass,$neutralLossName,$neutralLossMass,$targets);

		push @{$this->{ref}{ModificationInfos}{ModificationInfo}}, {Symbol => $modSymbol, Value => $modValue};

		if ($mod->{Nme} eq 'Phospho') {
		    $this->{ref}{Phosphorylation}{Symbol} = $modSymbol;
		    # Phospho (ST) is separated from other phosphos in Paragon files,
		    # but PhosphoRS does not seem to care about provided phosphorylable residues and always put 'STY'
		}

	    }

	    # Spectra
	    while (my $paragonSpectrum = shift @phosphoSpectra) {
		# Peaks
		my $peakString = '';
		my @peaks = split /\n/, $paragonSpectrum->{Peaks};
		shift @peaks; # 1st line is empty
		foreach my $peakData (@peaks){
		    my ($mz,$charge,$intensity) = split /\t/, $peakData;
		    $peakString .= "$mz:$intensity,";
		}
		$peakString =~ s/,$//;
		my $peaks = { content => $peakString };

		# Peptides
		my @peptides;
		while (my $match = shift @{$paragonSpectrum->{Matches}}) {
		    # Mods
		    my @modList = map {0} (1..length($match->{seq}));
		    while (my $modFeature = shift @{$match->{Mods}}) {
			# Be sure to have all phosphos in 1 entity
			my $modSymbol = ($modFeature->{mod} =~ /^Phospho\(?/)? $this->{ref}{Phosphorylation}{Symbol} : $modSymbol{$modFeature->{mod}};
			$modList[$modFeature->{pos} - 1] = $modSymbol;
		    }
		    my $modString = '0.'.join('',@modList).'.0';
		    # Subs
		    if ($match->{Subs}) {
			my @seqList = split //, $match->{seq};
			while (my $sub = shift @{$match->{Subs}}) {
			    $seqList[$sub->{pos}-1] = $sub->{'sub'};
			}
			$match->{seq} = join '',@seqList;
		    }


		    push @peptides, { ID => $match->{ID}, Sequence => $match->{seq}, ModificationInfo => $modString};
		}
		my $peptides = { Peptide => \@peptides};

		push @{$this->{ref}{'Spectra'}{'Spectrum'}}, { ID => $paragonSpectrum->{ID},
							      PrecursorCharge => $paragonSpectrum->{charge},
							      ActivationTypes => $this->{activationType},
							      Peaks => $peaks,
							      IdentifiedPhosphorPeptides => $peptides}
	    }

	}

	sub start_element{
	    my ($this,$element) = @_;

	    if ($this->_isOn('DataDictionary')) {
		if ($element->{Name} eq 'El') {
		    push @el, {};
		}
		elsif ($element->{'Name'} eq 'Mod'){
		    push @mod, {};
		}
	    }
	    elsif($element->{Name} eq 'SPECTRUM'){
		$currentSpectrumID++;
		$currentMatchID = 0;
		$currentSpectrum = { ID => $currentSpectrumID , charge => $element->{'Attributes'}->{'{}charge'}->{'Value'} };
	    }
	    elsif($this->_isOn('SPECTRUM')){
		if ($element->{Name} eq 'MATCH') {
		    $currentMatchID++;
		    $currentMatch = { ID => $currentMatchID, seq => $element->{'Attributes'}->{'{}seq'}->{'Value'}}
		}
		elsif($this->_isOn('MATCH') && $element->{Name} eq 'MOD_FEATURE'){
		    push @currentMods, { mod => $element->{'Attributes'}->{'{}mod'}->{'Value'}, pos => $element->{'Attributes'}->{'{}pos'}->{'Value'}};
		}
		elsif($this->_isOn('MATCH') && $element->{Name} eq 'SUBSTITUTION_FEATURE'){
		    push @currentSubs, { 'sub' => $element->{'Attributes'}->{'{}sub'}->{'Value'}, pos => $element->{'Attributes'}->{'{}pos'}->{'Value'}};
		}
	    }

	    $this->_toggleElement($element,1);
	}

	sub end_element{
	    my ($this,$element) = @_;

	    # Elements dictionary (used to compute modification masses)
	    if ($this->_isOn('DataDictionary')) {
		if ($this->_isOn('El')) {
		    if ($element->{Name} eq 'Iso') {
			push @{$el[-1]{Iso}}, {Mss => $currentMss, Pry => $currentPry};
			$currentMss = '';
			$currentPry = '';
		    }
		    elsif($element->{Name} eq 'Sym'){
			$el[-1]{Sym} = $currentSym;
			$currentSym = '';
		    }
		}
		elsif($this->_isOn('Mod')){
		    if ($element->{Name} eq 'Nme') {
			$mod[-1]{Nme} = $currentNme;
			$currentNme = '';
		    }
		    elsif($element->{Name} eq 'Fma'){
			$mod[-1]{Fma} = $currentFma;
			$currentFma = '';
		    }
		    elsif($element->{Name} eq 'RpF'){
			$mod[-1]{RpF} = $currentRpF;
			$currentRpF = '';
		    }
		    elsif($element->{Name} eq 'NLF'){
			$mod[-1]{NLF} = $currentNLF if $currentNLF =~ /\w/;
			$currentNLF = '';
		    }
		    elsif($element->{Name} eq 'Tgt'){
			push @{$mod[-1]{Tgt}}, $currentTgt;
			$currentTgt = '';
		    }
		}
	    }
	    elsif($this->_isOn('SPECTRUM')){
		if ($element->{Name} eq 'MSMSPEAKS') {
		    $currentSpectrum->{Peaks} = $currentPeakData;
		    $currentPeakData = '';
		}
		elsif($element->{Name} eq 'MATCH'){
		    # check if match contains phosphos
		    if (scalar grep {$_->{mod} eq 'Phospho'} @currentMods) {
			@{$currentMatch->{Mods}} = @currentMods;
			@{$currentMatch->{Subs}} = @currentSubs if scalar @currentSubs;
			push @currentPhosphoMatches, $currentMatch;
			foreach my $mod (@currentMods){ $interestVMods{$mod->{mod}} = 1 }; # keep only mods found in phosphopeptides
		    }
		    @currentMods = ();
		    @currentSubs = ();
		    $currentMatch = {};
		}
		elsif($element->{Name} eq 'SPECTRUM'){
		    if (scalar @currentPhosphoMatches) {
			@{$currentSpectrum->{Matches}} = @currentPhosphoMatches;
			push @phosphoSpectra, $currentSpectrum;
		    }
		    @currentPhosphoMatches = ();
		    $currentSpectrum = {};
		    print "<!-- -->\n" unless $currentSpectrumID % 1000; # avoids time out
		}
	    }

	    $this->_toggleElement($element,0);
	}

	sub characters{
	    my ($this,$element) = @_;

	    if ($this->_isOn('DataDictionary')) {

		# El
		if ($this->_isOn('El')) {
		    if ($this->_isOn('Sym')) {
			$currentSym .= $element->{Data};
		    }
		    elsif($this->_isOn('Iso')){
			if ($this->_isOn('Mss')) {
			    $currentMss .= $element->{Data};
			}
			elsif($this->_isOn('Pry')){
			    $currentPry .= $element->{Data};
			}
		    }
		}

		# Mod
		elsif ($this->_isOn('Mod')){
		    if ($this->_isOn('Nme')) {
			$currentNme .= $element->{'Data'};
		    }
		    elsif($this->_isOn('Fma')){
			$currentFma .= $element->{'Data'};
		    }
		    elsif($this->_isOn('RpF')){
			$currentRpF .= $element->{'Data'};
		    }
		    elsif($this->_isOn('NLF')){
			$currentNLF .= $element->{'Data'};
		    }
		    elsif($this->_isOn('Tgt')){
			$currentTgt .= $element->{'Data'};
		    }
		}
	    }
	    # Peaks
	    elsif($this->_isOn('SPECTRUM') && $this->_isOn('MSMSPEAKS')){
		$currentPeakData .= $element->{'Data'};
	    }

	}

	sub _toggleElement{
	    my ($this,$element,$boolean) = @_;

	    $this->{isOn}{$element->{Name}} = $boolean;
	}

	sub _isOn{
	    my ($this,$elementName) = @_;
	    my $elName=($this->{isOn}->{$elementName})? $this->{isOn}->{$elementName} : 0;
	    return $elName;
	}

	sub _computeMolMass{
	    my ($this,$fma) = @_;

	    my $mass = 0;

	    while ($fma =~ /([A-Z][a-z]?)(\d+)*/g) {
			my $sym = $1;
			my $n = (defined $2)? $2 : 1;

			$mass += $this->{atomMass}{$sym} * $n;
	    }

	    return $mass;
	}

	sub _incrementSymbol{
	    my $symbol = shift;

	    my @alphaNum = ((0..9),('A'..'Z'),('a'..'z'));
	    my $nextSymbol;

	    for(my $i=0;$i<=$#alphaNum;$i++){
			if ($alphaNum[$i] eq $symbol) {
				$nextSymbol = (defined $alphaNum[$i+1])? $alphaNum[$i+1] : undef;
				last;
			}
	    }

	    unless ($nextSymbol) {die 'Too many modifications on phosphopeptides. The process aborted.'}

	    return $nextSymbol;
	}
    }
}

#sub _alreadyThere{
#    my ($peptides, $sequence) = @_;
#
#    foreach my $phosphoPep ( @$peptides ){
#	if($phosphoPep->{Sequence} eq $sequence){
#	    return $phosphoPep->{ID};
#	}
#    }
#    return undef;
#}
#
#sub _matchPhosphoPeptide{
#    my ($this, $query,$isoPepID, $pepID) = @_;
#
#    $this->{isoformMatch}{$query}{$pepID} = $isoPepID;
#}
#
#sub _getMatchedPhosphoPeptide{
#    my ($this, $query, $pepID) = @_;
#
#    return $this->{isoformMatch}{$query}{$pepID};
#}

sub _writeInputFile{
    my $this = shift;
    #my %promsPath = &promsConfig::getServerInfo('no_user');
    #unless(-d "$promsPath{tmp}/phosphoRS"){
    #    mkdir "$promsPath{tmp}/phosphoRS";
    #}
    my $ref = $this->_getInputRef;

    # Cleaning undefined spectra #
    my $i=0;
    my $nbSpectra=0;
    while($i<=$#{$ref->{'Spectra'}{'Spectrum'}}){
		if($ref->{'Spectra'}{'Spectrum'}[$i]){
			$i++;
		}
		else{
			splice(@{$ref->{'Spectra'}{'Spectrum'}}, $i, 1);
		}
		$nbSpectra++;
		if ($nbSpectra > 100000) {
		    $nbSpectra=0;
		    print ".";
		}
    }
    # Cleaning undefined modification infos (when there are previously merges PTMs)
    my $j=0;
    while($j<=$#{$ref->{'ModificationInfos'}{'ModificationInfo'}}){
		if($ref->{'ModificationInfos'}{'ModificationInfo'}[$j]){
			$j++;
		}
		else {
			splice(@{$ref->{'ModificationInfos'}{'ModificationInfo'}}, $j , 1);
		}
    }
    my $xml = XMLout($ref, RootName => 'phosphoRSInput'); # This command can take a lot of time...
    #$xml =~ s/[\n\r][\s]*<Spectrum><\/Spectrum>//g; # cleaning empty elements
    my $xml2=""; # was added on July 9th 2015 to avoid Substitution loop
    foreach my $line (split(/\n/,$xml)) {
	    $line=~ s/[\s]*<Spectrum><\/Spectrum>//g; # cleaning empty elements
	    $xml2.="$line\n" unless $line eq '';
    }
    #my $inputFile = "$promsPath{'tmp'}/phosphoRS/".strftime("%Y%m%d%H%M%S",localtime).'.xml';
    my $inputFile = $this->{'fullJobDir'}.'/input_data.xml';
	open(XMLfile, ">$inputFile");
    print XMLfile $xml2;
    close(XMLfile);
    $this->{'inputFile'} = $inputFile;
    return $inputFile;
}

sub _setModifValue{ #static
    # used to modify modifications infos in input data
    # undef values are not changed
    my ($ref, $symbol, @newValues) = @_;

    # finding modification entity
    my ($modificationInfo,$modInfoIndex) = _getModificationInfo($ref,$symbol);
    my @currentValues = split(/:/, $modificationInfo->{Value});
    die "Processing error" unless scalar(@currentValues) == scalar(@newValues);
    for(my $i=0;$i<=$#currentValues;$i++){
        $currentValues[$i] = $newValues[$i] if $newValues[$i]; #skip undef
    }
    $modificationInfo->{'Value'} = join(':',@currentValues);

    $ref->{'ModificationInfos'}{'ModificationInfo'}[$modInfoIndex] = $modificationInfo;
}

sub _getModifValue{
    my ($ref, $symbol) = @_;

    my ($modificationInfo) = _getModificationInfo($ref,$symbol);

    return map {(defined $_ && length($_))? $_ : undef} split(/:/, $modificationInfo->{'Value'});
}

sub _getModificationInfo{
    my ($ref, $symbol) = @_;

    my $modificationInfo;
    my $modInfoIndex;
    for(my $i=0;$i<=$#{$ref->{'ModificationInfos'}{'ModificationInfo'}};$i++){
		if ($ref->{'ModificationInfos'}{'ModificationInfo'}[$i] && $ref->{'ModificationInfos'}{'ModificationInfo'}[$i]{Symbol} eq $symbol) {
			$modificationInfo = $ref->{'ModificationInfos'}{'ModificationInfo'}[$i];
			$modInfoIndex = $i;
			last;
		}
    }

    return ($modificationInfo,$modInfoIndex);
}

sub _incrementSymbol{
    my $symbol = shift;

    my @alphaNum = ((0..9),('A'..'Z'),('a'..'z'));
    my $nextSymbol;

    for(my $i=0;$i<=$#alphaNum;$i++){
		if ($alphaNum[$i] eq $symbol) {
			$nextSymbol = (defined $alphaNum[$i+1])? $alphaNum[$i+1] : undef;
			last;
		}
    }

    unless ($nextSymbol) {die 'Too many modifications on phosphopeptides. The process aborted.'}

    return $nextSymbol;
}

sub _storeResults{
    my ($this,$outputFile) = @_;
    $this->{'outputFile'} = $outputFile;
    $this->{outputRef} = XMLin($outputFile, ForceArray => ['Spectrum', 'Peptide', 'Isoform', 'PhosphoSite']);
    my @dataFilePath = split /\//, $this->{'dataFile'};
    my $dataFileName = pop @dataFilePath;
    $dataFileName =~ s/\.\w{3}$//;
    my $dir = join '/', @dataFilePath;
	my $anaID = $this->{'analysisID'};
    my $outputFileName = "$dir/PRS_ana_$anaID.xml";
	copy($outputFile,$outputFileName);
}

sub _getInputRef{
    my $this = shift;
    return $this->{inputRef};
}

sub _setInputRef{
    my ($this,$ref) = @_;
    $this->{inputRef} = $ref;
}
1;

####>Revision history<####
# 1.3.1 Add in _storeResults a copy of xml file (GA 12/12/18)
# 1.3.0 Major code update for improved background run (PP 09/11/18)
# 1.2.1 Changed $maxHours from 12 to 48 (PP 19/06/18)
# 1.2.0 Now uses &promsConfig::clusterInfo (PP 15/06/18)
# 1.1.9 Minor change in &cleanResults to handle unexpected undef $queryNum (PP 25/05/18)
# 1.1.8 Minor modif in _addSpectraFromMascotFile for modifications names containing final space (MLP 27/03/18)
# 1.1.7 Uses xxx_ana_&lt;anaID&gt;.xxx instead of data file name for PRS output file (PP 21/03/17)
# 1.1.6 Change in &printIcon for (PP 22/021/17)
# 1.1.5 Improved phospho.jar error management (21/10/16)
# 1.1.4 Uses $promsPath{java} if declared (PP 20/09/16)
# 1.1.3 Change maxHours (GA 14/03/16)
# 1.1.2 Changes to keep server connexion for big job on cluster  & minor change in &printIcon function (PP 09/09/15)
# 1.1.1 Move java call to cluster (GA 10/07/15)
# 1.1.0 Accepts MASCOT|SEQUEST.PDM file & bug fix for modification names containing ':' (PP 23/12/14)
# 1.0.9 Change syntax $x=$y // $z to old ()? one (PP 11/03/14)
# 1.0.8 Handling PARAGON format
# 1.0.7 Add filtering method to keep only validated queries and decrease result-file size (FY 08/04/13)
# 1.0.6 Deleting empty ModificationInfo markups to avoid code 107 errors (FY 20/12/12)
# 1.0.5 Removed global variable declarations (now local) (FY 07/11/12)
# 1.0.4 Add status 4 for no isoform case (FY 06/11/12)
# 1.0.3 PRS string parsing bug fix of in printIcon method (FY 04/10/12)
# 1.0.2 Testing file names definition before unlink commands<BR>& Management of any kind of phosphorylation (FY 12/06/12)
# 1.0.1 Add method to fill incomplete spectra data in Proteome Discoverer analyses +<BR> Add score extraction from XML output file+<BR> Unknown neutral loss value set to 0 instead of 'null' (FY 21/05/12)
# 1.0.0 New class used as an interface between myProMS and PhosphoRS java application (FY 01/03/12)