####################################################################################
# promsMod.pm           4.0.6                                                      #
# Authors: P. Poullet, G. Arras, S.Liva, M. Le Picard, V. Sabatet (Institut Curie) #
# Contact: myproms@curie.fr                                                    	   #
####################################################################################
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
package promsMod;
require Exporter;

@ISA=qw(Exporter);
@EXPORT=qw();
@EXPORT_OK=qw();
$VERSION=1.00;

use POSIX qw(strftime);
use strict;
use IO::File;  # Core module
use IO::Uncompress::Unzip qw(unzip $UnzipError); # needed by &extractSpectrumMSF
use File::Spec::Functions qw(splitpath); # Core module
use File::Path qw(mkpath rmtree); # remove_tree
use File::Copy qw(move);

###############################
####>>>myProMS Constants<<<####
###############################

##>Aliases for Status
sub getStatusAlias {
	my %statusAlias=('bioinfo'=>'Bioinformatician','mass'=>'Massist','manag'=>'Manager','bio'=>'Biologist');
	return $statusAlias{$_[0]};
}
##>Aliases for Profile
sub getProfileAlias {
	my %profileAlias=(
		'bioinfo'=>'Bioinformatician',
		'mass'=>'Massist',
		'manag'=>'Manager',
		'guest'=>'Guest',
		'user'=>'User',
		'power_user'=>'Power User',
		'super_user'=>'Super User',
		'administrator'=>'Administrator',
		'power_administrator'=>'Power Administrator',
		'super_administrator'=>'Super Administrator'
	);
	return $profileAlias{$_[0]};
}

##>Generic Names for tables
my %type=(
	'PROJECT'=>'Project','EXPERIMENT'=>'Experiment','SAMPLE'=>'Sample','ANALYSIS'=>'Analysis',
	'GEL2D'=>'2D-Gel','SPOT'=>'Spot',
	'GEL2D SAMPLE'=>'2D-Gel or free Sample','SAMPLE GEL2D'=>'free Sample or 2D-Gel',
	'SPOT SAMPLE'=>'Spot or free Sample','SAMPLE SPOT'=>'free Sample or Spot',
	'SPECTRUM'=>'Spectrum',
	'PEPTIDE'=>'Peptides','PROTEIN'=>'Proteins', #'PROT_CLUSTER'=>'Protein Cluster','GENE'=>'Gene',
	'CATEGORY'=>'List','CLASSIFICATION'=>'Custom Lists',
	'DATABANK'=>'Databank',
	'GO_ANALYSIS'=>'GO analysis',
	'DESIGN'=>'Design',
	'EXPCONDITION'=>'Condition',
	'QUANTIFICATION'=>'Quantification',
	'SPECIES'=>'Species',
	'EXPLORANA' => 'Exploratory Analyses',
	'METADATA' => 'Metadata',
	'PATHWAY_ANALYSIS' => 'Pathway analysis'
);

##> AA code
my %aaThree2OneLetter=(
	'Ala'=>'A','Arg'=>'R','Asp'=>'D','Asn'=>'N','Cys'=>'C','Glu'=>'E','Gln'=>'Q','Gly'=>'G',
	'His'=>'H','Ile'=>'I','Leu'=>'L','Lys'=>'K','Met'=>'M','Phe'=>'F','Pro'=>'P','Ser'=>'S',
	'Thr'=>'T','Trp'=>'W','Tyr'=>'Y','Val'=>'V'
);
sub getAAThree2OneLetter {
	return %aaThree2OneLetter;
}


##########################################
####<<< Fetching item simple infos >>>#### $_[0]: item type
##########################################
sub getItemParent {
	my %parentTable=(
		'PROJECT_ACCESS'=>['USER_LIST','USER_PROFILE'],
		'PROJECT'=>['PROJECT_ACCESS'],
		'EXPERIMENT'=>['PROJECT'],
		'SAMPLE'=>['EXPERIMENT','GEL2D'],
		'ANALYSIS'=>['SAMPLE'],
		'GEL2D'=>['EXPERIMENT'],
		'SPOT'=>['GEL2D'],
		'QUERY_VALIDATION'=>['ANALYSIS'],
		'SPECTRUM'=>['ANALYSIS'],
		'PEPTIDE'=>['ANALYSIS'],
		'DESIGN'=>['EXPERIMENT'],
		'EXPCONDITION'=>['DESIGN'],
		'QUANTIFICATION'=>['DESIGN'],
		'EXPLORANALYSIS' => ['EXPERIMENT'],
		'PATHWAY_ANALYSIS' => ['EXPERIMENT']
	);
	return @{$parentTable{uc($_[0])}}; # <= an array!!! (in case more than 1 parent)
}
sub getItemChild {
	my %childTable=('PROJECT'=>['EXPERIMENT'],
					'EXPERIMENT'=>['SAMPLE','GEL2D'],
					'GEL2D'=>['SPOT'],
					'SPOT'=>['ANALYSIS'], #['SAMPLE']
					'SAMPLE'=>['ANALYSIS'],
					'ANALYSIS'=>[''],
					'DESIGN'=>['EXPCONDITION','QUANTIFICATION'],
					'EXPCONDITION'=>[''],
					'QUANTIFICATION'=>[''],
					'EXPLORANALYSIS' => [''],
					'PATHWAY_ANALYSIS' => ['']
					);
	return @{$childTable{uc($_[0])}}; # <= an array!!! (in case more than 1 child)
}
sub getItemType {
	return $type{uc($_[0])};
}
# sub getItemPkName {
# 	return $pkColumn{uc($_[0])};
# }
#sub getItemFkName {
#	return $fkColumn{uc($_[0])};
#}
sub getItemPlurial {
 	my $item=uc($_[0]);
	#if ($item eq 'ANALYSIS'){
	#	return 'MS/MS Analyses';
	#}
	#elsif($item eq 'SAMPLE'){
	#	return 'Samples';
	#}
	#else{
	#	return 'Experiments';
	#}
	my $itemPlurial=($item eq 'ANALYSIS')? 'Analyses' : $type{$item}.'s';
	return $itemPlurial;
}

###########################################
####<<< Identifier Type description >>>####
###########################################
sub getIdentifierTypes {
	my %identifierTypes=('GI_ACCESSION'=>'NCBI gi number (gi|125987826)',
						 'NCBI_ALL'=>'NCBI all (gi|125987826|sp|P15311|EZRI_HUMAN)',
						 'UNIPROT_ACCESSION'=>'Uniprot accession (P15311)',
						 'UNIPROT_ID'=>'Uniprot id (EZRI_HUMAN)',
						 'UNIPROT_ALL'=>'Uniprot all (sp|P15311|EZRI_HUMAN)',
						 'REFSEQ_PROTEIN'=>'NCBI RefSeq Protein (NP_003370.2)',
						 'IPI_ACCESSION'=>'International Protein Index (IPI00843975)',
						 'FLYBASE_ID'=>'FlyBase Polypeptide (FBpp0085360)',
						 'UPS_ACCESSION'=>'UPS accession (P02768ups)',
						 'UPS_ID'=>'UPS id (ALBU_HUMAN_UPS)',
						 'UPS_ALL'=>'UPS all (P02768ups|ALBU_HUMAN_UPS)',
						 );
	return %identifierTypes;
}

#####################################################
####<<< Clean parameters used in system calls >>>#### (to prevent code injection)
#####################################################
sub cleanParameters {
	my @params;
	foreach my $param (@_) {
		$param=~s/;.*//g if $param; # remove ";" and everything after that
		push @params,$param;
	}
	return ($#params)? @params : $params[0];
}
sub cleanNumericalParameters {
	my @params;
	foreach my $param (@_) {
		$param=~s/\D+//g if $param; # remove everything after that is not a number
		push @params,$param;
	}
	return ($#params)? @params : $params[0];
}
####################################################################### $_[0]: string to be resize
####<<< Resize a string to a maximum size by truncating its end >>>#### $_[1]: maximum size wanted
#######################################################################
sub resize {
	my $string=$_[0] || '';
	my $maxSize=$_[1] || 100;
	if (length($string)>$maxSize) {
		my $newSize=$maxSize-3;
		$string=~s/^(.{$newSize}).*/$1/;
		$string.='...';
	}
	return $string;
}

########################################################################### $_[0]: string to be resize
####<<< Shorten a string to a maximum size by truncating its middle >>>#### $_[1]: maximum size wanted
###########################################################################
sub shortenName {
	my $text=($_[0])? $_[0] : '';
	my $maxLength=($_[1])? $_[1] : length($text);
	return $text if length($text) <= $maxLength;
	my $subTextSize=int((($maxLength-3)/2));
	my ($begText,$endText)=($text=~/^(.{$subTextSize}).*(.{$subTextSize})/);
	return "$begText...$endText";
}

###########################################
####<<< Make string HTML compatible >>>####
###########################################
sub HTMLcompatible {
	my $string=($_[0])? $_[0] : '';
# 	$string=~s/ /&nbsp;/g; # text cannot wrap
	$string=~s/</&lt;/g;
	$string=~s/>/&gt;/g;
	$string=~s/\n/<BR>/g;
	$string=~s/\r//g; # multi-line fields use '\r\n'!!!
	$string=~s/"/&quot;/g;
	$string=~s/'/&apos;/g;
	return $string;
}

#######################################
####<<< Make date more readable >>>####
#######################################
sub formatDate { # argument = date from DB
	return '--/--/----' unless $_[0];
	my ($date,$time)=split(/ /,$_[0]);
	my ($year,$month,$day)=split(/-/,$date);
	$year+=2000 if $year<2000; # oracle
	return "$day/$month/$year $time";
}

#############################################
####<<< Sort for alphanumeric strings >>>####
#############################################
sub sortSmart {
	my ($param1,$param2,$level,$ref_i,$lastIdx)=@_;
	$level=0 unless $level;
	if ($level==0) {
		$param1='?' unless $param1;
		$param2='?' unless $param2;
		my @parData1=split(/(\d+|[a-z]+|\W+)/,lc($param1));
		my @parData2=split(/(\d+|[a-z]+|\W+)/,lc($param2));
		#<Removing undef indexes (1st index & between number/text)
		foreach my $refArray (\@parData1,\@parData2) {
			my $i=0;
			while ($i <= $#{$refArray}) {
				if (length($refArray->[$i])) {$i++;}
				else {splice(@{$refArray},$i,1);}
			}
		}
		my $lastIdx=$#parData1;
		if ($#parData1 < $#parData2) {push @parData1,'';} # 1 more is enough
		elsif ($#parData1 > $#parData2) {push @parData2,'';} # 1 more is enough
		my $i=0;
		return &sortSmart(\@parData1,\@parData2,1,\$i,$lastIdx);
	}
	else { # level 1
		my $currentSort=0;
		if ($param1->[$$ref_i] && $param1->[$$ref_i] !~ /\D/ && $param2->[$$ref_i] && $param2->[$$ref_i] !~ /\D/) {
			$param1->[$$ref_i]=~s/^0+([^0])/$1/; #00100 -> 100 but 000 -> 000
			$param2->[$$ref_i]=~s/^0+([^0])/$1/;
			$currentSort=($param1->[$$ref_i]<=>$param2->[$$ref_i]);
		}
		else {$currentSort=($param1->[$$ref_i] cmp $param2->[$$ref_i]);}
		while ($currentSort==0 && $$ref_i<$lastIdx) {
			$$ref_i++;
			$currentSort=&sortSmart($param1,$param2,1,$ref_i,$lastIdx);
		}
		return $currentSort;
	}
}
sub preparePtmString { # S145.T234 -> 145.234. for better sort by sortSmart
	my ($ptmStrg)=@_;
	return '' unless $ptmStrg;
	$ptmStrg=~s/^-//;
	return $ptmStrg if $ptmStrg=~/\+/; # pos ambiguity format
	my $revStrg='';
	foreach my $modPos (split(/\./,$ptmStrg)) {
		#my ($res,$pos)=($modPos=~/^(.)(\d+)/);
		#$revStrg.=$2.$1;
		my ($pos)=($modPos=~/(\d+)$/);
		$revStrg.=$1.'.';
	}
	return $revStrg;
}

sub MaxQuantProbIcon {
	my ($prob,$refOpt)=@_;
	$refOpt={} unless $refOpt;
	$refOpt->{text}='&nbsp;MQ&nbsp;' unless defined $refOpt->{text};
	my $qStrg=($refOpt->{inPopup})? '\\' : '';
	my $iconStrg="<FONT style=${qStrg}'font-size:10px;background-color:rgb(";
	$iconStrg.=($prob >= 0.5)? (255-sprintf '%.0f',($prob-0.5)*510).',255,0' : sprintf '255,%.0f,0',$prob*510;
	$iconStrg.=")$qStrg'";
	$iconStrg.=" onmouseover=\"popup('$refOpt->{popupText}')\" onmouseout=\"popout()\"" if $refOpt->{popupText};
	$iconStrg.=">$refOpt->{text}</FONT>";
	return $iconStrg;
}


##################################
####<<< Fetching user info >>>####
##################################
sub getUserInfo {
	my $dbh=$_[0];
	my $qUserID=$dbh->quote($_[1]);
	my $projectID=$_[2]; # optional argument

	my %userProfile;

	###<Fetching user details>###
	my ($userName,$userStatus,$userLab,$userTel,$userEmail,$userInfo,$mascotIdStrg,$workgroup) = $dbh->selectrow_array("SELECT USER_NAME,USER_STATUS,LABORATORY,USER_TEL,USER_EMAIL,USER_INFO,MASCOT_IDS,WORK_GROUP FROM USER_LIST WHERE ID_USER = $qUserID");
	($userLab,$userTel,$userEmail,$userInfo,$mascotIdStrg,$workgroup)=&chkDef($userLab,$userTel,$userEmail,$userInfo,$mascotIdStrg,$workgroup);
	my @mascotIdList=split(',',$mascotIdStrg); # 'userID,groupID1,groupID2,...' userID can be empty

	###<All access status => return>###
	if ($userStatus=~/bioinfo|mass/) {
		###<Fetching user Project(s)>###
		if (defined $projectID) {
			##<Fetching data for a specific project>##
			$userProfile{$projectID}=$userStatus;
		}
		else {
			##<Fetching data for all projects>##
			my $sthAP=$dbh->prepare("SELECT ID_PROJECT FROM PROJECT");
			$sthAP->execute;
			while (my ($projectID)=$sthAP->fetchrow_array) {
				$userProfile{$projectID}=$userStatus;
			}
			$sthAP->finish;
		}
	}
	else { # user is a biologist or manager
		###<Fetching user accessible Project(s)>###
		if (defined $projectID) {
			##<Fetching data for a specific project>##
			if ($userStatus eq 'manag') {
				my ($projWorkgroup)=$dbh->selectrow_array("SELECT WORK_GROUP FROM PROJECT WHERE ID_PROJECT=$projectID");
				$projWorkgroup='' unless $projWorkgroup;
				$userProfile{$projectID}=$userStatus if ($workgroup && $workgroup eq $projWorkgroup);
			}
			if (!$userProfile{$projectID}) { # bio or manager with access to another Workgroup's project
				next if $userProfile{$projectID}; # in case manager still have undeleted project-based access (just to be safe)
				($userProfile{$projectID})=$dbh->selectrow_array("SELECT USER_PROFILE.NAME FROM USER_PROFILE,PROJECT_ACCESS WHERE PROJECT_ACCESS.ID_USER=$qUserID AND PROJECT_ACCESS.ID_PROJECT=$projectID AND USER_PROFILE.ID_PROFILE=PROJECT_ACCESS.ID_PROFILE");
			}
		}
		else {
			##<Fetching data for all projects associated with user>##
			if ($userStatus eq 'manag' && $workgroup) { # manager w/o workgroup has access to nothing!
				#my $sthWG=($workgroup)? $dbh->prepare("SELECT ID_PROJECT FROM PROJECT WHERE WORK_GROUP=?") : $dbh->prepare("SELECT ID_PROJECT FROM PROJECT WHERE WORK_GROUP=? OR WORK_GROUP IS NULL");
				my $sthWG=$dbh->prepare("SELECT ID_PROJECT FROM PROJECT WHERE WORK_GROUP=?");
				$sthWG->execute($workgroup);
				while (my ($projectID)=$sthWG->fetchrow_array) {
					$userProfile{$projectID}=$userStatus;
				}
				$sthWG->finish;
			}
			#<Granted access to other projects
			my $sthGA=$dbh->prepare ("SELECT PROJECT_ACCESS.ID_PROJECT,USER_PROFILE.NAME FROM USER_PROFILE,PROJECT_ACCESS WHERE PROJECT_ACCESS.ID_USER=$qUserID AND PROJECT_ACCESS.ID_PROFILE=USER_PROFILE.ID_PROFILE");
			$sthGA->execute;
			while (my ($projectID,$profileName)=$sthGA->fetchrow_array) {
				next if $userProfile{$projectID}; # in case manager still have undeleted project-based access (just to be safe)
				$userProfile{$projectID}=$profileName;
			}
			$sthGA->finish;
		}
	}
	return ($userName,$userStatus,\%userProfile,$userLab,$userTel,$userEmail,$userInfo,\@mascotIdList,$workgroup);
}

######################################################
####<<< Fetching parent data for provided item >>>####
######################################################
sub getItemInfo {
	my ($dbh,$startItem,$startID,$childID)=@_;
	$startItem=uc($startItem);
	my @itemInfo; # Contains item hierarchy  (project--->analysis) [i]{'ID'=>$itemID,'NAME'=>$itemName,'ITEM'=>$itemList[$i],'TYPE'=>$itemType,'POS'=>$displayPos,...} Optional: Analysis='VALID','DECOY'; project='STATUS'
	my @itemList=('ANNOTATION_ITEM','META_ANNOTATION','EXPLORANALYSIS','ANALYSIS','SAMPLE','QUANTIFICATION','DESIGN','SPOT','GEL2D','EXPERIMENT','PROJECT');

	my $itemOK=0;
	my $parentID=$startID;
	my $i=0;
	while ($i<=$#itemList) {
		$itemOK=1 if $itemList[$i] eq $startItem;
		unless ($itemOK) {
			$i++;
			next;
		}
		my $itemName='';
		my $itemID=$parentID;
		my $increment=1;
		if ($itemList[$i] eq 'ANALYSIS') {
			($itemName,$parentID,my $displayPos,my $validStatus,my $decoy)=$dbh->selectrow_array("SELECT NAME,ID_$itemList[$i+1],DISPLAY_POS,VALID_STATUS,DECOY FROM ANALYSIS WHERE ID_ANALYSIS=$itemID");
			unshift @itemInfo,{'ID'=>$itemID,'NAME'=>$itemName,'ITEM'=>$itemList[$i],'TYPE'=>$type{$itemList[$i]},'POS'=>$displayPos,'VALID'=>$validStatus,'DECOY'=>$decoy}; # add at begining of array
		}
		elsif ($itemList[$i] eq 'PROJECT') {
			($itemName,my $projectStatus)=$dbh->selectrow_array("SELECT NAME,STATUS FROM PROJECT WHERE ID_PROJECT=$itemID");
			$projectStatus=0 unless $projectStatus;
			unshift @itemInfo,{'ID'=>$itemID,'NAME'=>$itemName,'ITEM'=>$itemList[$i],'TYPE'=>$type{$itemList[$i]},'POS'=>0,'STATUS'=>$projectStatus}; # add at begining of array
		}
		else {
			my $displayPos=0;
			if ($itemList[$i] eq 'SAMPLE') {
				($itemName,my $expID,my $spotID,$displayPos)=$dbh->selectrow_array("SELECT NAME,ID_EXPERIMENT,ID_SPOT,DISPLAY_POS FROM SAMPLE WHERE ID_SAMPLE=$itemID");
				if ($spotID) {$parentID=$spotID;} # 2D Gel
				else { # no 2D Gel
					$parentID=$expID;
					$increment+=2; # skip SPOT and GEL2D tables
				}
				$increment+=2; # skip QUANTIFICATION and DESIGN tables
			}
			elsif ($itemList[$i] eq 'EXPERIMENT') {
				($itemName,$parentID,$displayPos)=$dbh->selectrow_array("SELECT NAME,ID_PROJECT,DISPLAY_POS FROM EXPERIMENT WHERE ID_EXPERIMENT=$itemID");
			}
			elsif ($itemList[$i] eq 'QUANTIFICATION') {
#print "QUANTI!<BR>";
				($itemName,$parentID,my $focus)=$dbh->selectrow_array("SELECT NAME,ID_DESIGN,FOCUS FROM QUANTIFICATION WHERE ID_QUANTIFICATION=$itemID");
				###> Change on 22/10/2012 to print the XIC extraction in selectOptionQuanti.cgi
				unless ($parentID) { # branch to analysis
					if ($focus eq 'protein') {
						($parentID)=$dbh->selectrow_array("SELECT ID_ANALYSIS FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$itemID LIMIT 0,1");
						$increment=-2; # back to analysis
						#branch to experiment? ($expID)=$dbh->selectrow_array("SELECT S.ID_EXPERIMENT FROM ANALYSIS A, SAMPLE S WHERE A.ID_SAMPLE=S.ID_SAMPLE AND A.ID_ANALYSIS=$anaID"); # as long as analyses are directly linked to sample!
						#branch to parent quanti? ($parentID)=$dbh->selectrow_array("SELECT ID_PARENT_QUANTIFICATION FROM PARENT_QUANTIFICATION WHERE ID_QUANTIFICATION=$itemID");
					}
					else { # peptide
						if ($childID) {$parentID=$childID;} # XIC quanti link to a prot quanti (selectOptionQuanti.cgi)
						else {
							($parentID)=$dbh->selectrow_array("SELECT ID_ANALYSIS FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$itemID ORDER BY IS_REFERENCE DESC LIMIT 0,1");
							$increment=-2; # back to analysis
						}
					}
				}
				#$parentID=$childID unless ($parentID);
			}
			elsif ($itemList[$i] eq 'DESIGN'){
#print "DESIGN!<BR>";
				($itemName,$parentID)=$dbh->selectrow_array("SELECT NAME,ID_EXPERIMENT FROM DESIGN WHERE ID_DESIGN=$itemID");
				$increment+=2; # skip SPOT and GEL2D tables
			}
			elsif ($itemList[$i] eq 'EXPLORANALYSIS') {
				($itemName,$parentID)=$dbh->selectrow_array("SELECT NAME,ID_EXPERIMENT FROM EXPLORANALYSIS WHERE ID_EXPLORANALYSIS=$itemID");
				$increment+=6;
			}
			elsif ($itemList[$i] eq 'META_ANNOTATION' || $itemList[$i] eq 'ANNOTATION_ITEM') {
				if($itemList[$i] eq 'ANNOTATION_ITEM') {
					$itemID = $dbh->selectrow_array("SELECT ID_META_ANNOTATION FROM ANNOTATION_ITEM WHERE ID_ANNOTATION_ITEM=$itemID");
				}

				($itemName, my $projectID, my $experimentID, my $sampleID)=$dbh->selectrow_array("SELECT NAME, ID_PROJECT, ID_EXPERIMENT, ID_SAMPLE FROM META_ANNOTATION WHERE ID_META_ANNOTATION=$itemID");
				
				if($sampleID) {
					$parentID = $sampleID;
					$increment += 2; # skip until Sample
				} elsif($experimentID) {
					$parentID = $experimentID;
					$increment += 8; # skip until Experiment
				} else {
					$parentID = $projectID;
					$increment += 9; # skip until Project
				}
			} else {
				($itemName,$parentID)=$dbh->selectrow_array("SELECT NAME,ID_$itemList[$i+1] FROM $itemList[$i] WHERE ID_$itemList[$i]=$itemID");
			}
			unshift @itemInfo,{'ID'=>$itemID,'NAME'=>$itemName,'ITEM'=>$itemList[$i],'TYPE'=>$type{$itemList[$i]},'POS'=>$displayPos}; # add at begining of array
		}
		$i+=$increment unless ($itemList[$i] eq 'QUANTIFICATION' && $childID && $parentID == $childID);
	}
	return @itemInfo;
}

###################################
####<<< Fetching project ID >>>####
###################################
sub getProjectID {
	my ($dbh,$id,$item)=@_;
	$id=&cleanNumericalParameters($id);
	$item=uc($item);
	my $projectID;
	if ($item eq 'PROJECT') {
		$projectID=$id;
	}
	elsif ($item eq 'EXPERIMENT') {
		($projectID)=$dbh->selectrow_array("SELECT ID_PROJECT FROM EXPERIMENT WHERE ID_EXPERIMENT=$id");
	}
	elsif ($item eq 'SAMPLE') {
		($projectID)=$dbh->selectrow_array("SELECT ID_PROJECT FROM EXPERIMENT,SAMPLE WHERE EXPERIMENT.ID_EXPERIMENT=SAMPLE.ID_EXPERIMENT AND ID_SAMPLE=$id");
	}
	elsif ($item eq 'ANALYSIS') {
		($projectID)=$dbh->selectrow_array("SELECT ID_PROJECT FROM EXPERIMENT,SAMPLE,ANALYSIS WHERE ANALYSIS.ID_SAMPLE=SAMPLE.ID_SAMPLE AND EXPERIMENT.ID_EXPERIMENT=SAMPLE.ID_EXPERIMENT AND ID_ANALYSIS=$id");
		#>Uncomment when direct connection from ANALYSIS to SPOT
		#unless ($projectID) {
		#	($projectID)=$dbh->selectrow_array("SELECT ID_PROJECT FROM EXPERIMENT,GEL2D,SPOT,SAMPLE,ANALYSIS WHERE ANALYSIS.ID_SPOT=SPOT.ID_SPOT AND SPOT.ID_GEL2D=GEL2D.ID_GEL2D AND EXPERIMENT.ID_EXPERIMENT=GEL2D.ID_EXPERIMENT AND ID_ANALYSIS=$id");
		#}
	}
	elsif ($item eq 'GEL2D') {
		($projectID)=$dbh->selectrow_array("SELECT ID_PROJECT FROM EXPERIMENT,GEL2D WHERE EXPERIMENT.ID_EXPERIMENT=GEL2D.ID_EXPERIMENT AND ID_GEL2D=$id");
	}
	elsif ($item eq 'SPOT') {
		($projectID)=$dbh->selectrow_array("SELECT ID_PROJECT FROM EXPERIMENT,GEL2D,SPOT WHERE SPOT.ID_GEL2D=GEL2D.ID_GEL2D AND EXPERIMENT.ID_EXPERIMENT=GEL2D.ID_EXPERIMENT AND ID_SPOT=$id");
	}
	elsif ($item eq 'DESIGN') {
		($projectID)=$dbh->selectrow_array("SELECT ID_PROJECT FROM DESIGN,EXPERIMENT WHERE DESIGN.ID_EXPERIMENT=EXPERIMENT.ID_EXPERIMENT AND ID_DESIGN=$id");
	}
	elsif ($item eq 'EXPCONDITION') {
		($projectID)=$dbh->selectrow_array("SELECT ID_PROJECT FROM EXPCONDITION,DESIGN,EXPERIMENT WHERE DESIGN.ID_EXPERIMENT=EXPERIMENT.ID_EXPERIMENT AND DESIGN.ID_DESIGN=EXPCONDITION.ID_DESIGN AND ID_EXPCONDITION=$id");
	}
	elsif ($item eq 'QUANTIFICATION') {
		my ($analysisID)=$dbh->selectrow_array("SELECT ID_ANALYSIS FROM ANA_QUANTIFICATION WHERE ID_QUANTIFICATION=$id LIMIT 0,1");
		if ($analysisID) {$projectID=&getProjectID($dbh,$analysisID,'ANALYSIS');}
		else {($projectID)=$dbh->selectrow_array("SELECT ID_PROJECT FROM QUANTIFICATION,DESIGN,EXPERIMENT WHERE EXPERIMENT.ID_EXPERIMENT=DESIGN.ID_EXPERIMENT AND DESIGN.ID_DESIGN=QUANTIFICATION.ID_DESIGN AND ID_QUANTIFICATION=$id");}
	}
	elsif ($item eq 'GO_ANALYSIS') {
		($projectID)=$dbh->selectrow_array("SELECT ID_PROJECT FROM EXPERIMENT,GO_ANALYSIS WHERE EXPERIMENT.ID_EXPERIMENT=GO_ANALYSIS.ID_EXPERIMENT AND ID_GOANALYSIS=$id");
	}
	elsif ($item eq 'EXPLORANA') {
		($projectID)=$dbh->selectrow_array("SELECT ID_PROJECT FROM EXPLORANALYSIS,EXPERIMENT WHERE EXPLORANALYSIS.ID_EXPERIMENT=EXPERIMENT.ID_EXPERIMENT AND ID_EXPLORANALYSIS=$id");
	}
	elsif ($item eq 'PATHWAY_ANALYSIS') {
		($projectID)=$dbh->selectrow_array("SELECT ID_PROJECT FROM PATHWAY_ANALYSIS,EXPERIMENT WHERE PATHWAY_ANALYSIS.ID_EXPERIMENT=EXPERIMENT.ID_EXPERIMENT AND ID_PATHWAY_ANALYSIS=$id");
	}
	return $projectID;
}

######################################
####<<< Fetching item children >>>#### !!! ONLY for items containing proteins !!!
######################################
sub getChildrenList {
	my ($dbh,$id,$item)=@_;
	$item=uc($item);
	my @childrenList;
	my (@children,%query);
	if ($item eq 'SAMPLE'){
		@children=('ANALYSIS');
		$query{'ANALYSIS'}="SELECT ID_ANALYSIS,NAME,DES,DISPLAY_POS,VALID_STATUS FROM ANALYSIS WHERE ID_SAMPLE=$id";
	}
	elsif ($item eq 'SPOT'){
		@children=('ANALYSIS'); # virtual (skipping SAMPLE)!!!!!!!!!!!
		#$query{'SAMPLE'}="SELECT ID_SAMPLE,NAME,DES,DISPLAY_POS,1 FROM SAMPLE WHERE ID_SPOT=$id";
		$query{'ANALYSIS'}="SELECT ID_ANALYSIS,ANALYSIS.NAME,ANALYSIS.DES,ANALYSIS.DISPLAY_POS,VALID_STATUS FROM SAMPLE,ANALYSIS WHERE SAMPLE.ID_SAMPLE=ANALYSIS.ID_SAMPLE AND ID_SPOT=$id";
	}
	elsif ($item eq 'GEL2D'){
		@children=('SPOT');
		$query{'SPOT'}="SELECT ID_SPOT,NAME,DES,0,1 FROM SPOT WHERE ID_GEL2D=$id ORDER BY NAME ASC";
	}
	elsif ($item eq 'EXPERIMENT'){
		@children=('GEL2D','SAMPLE'); #,'DESIGN'
		$query{'GEL2D'}="SELECT ID_GEL2D,NAME,DES,DISPLAY_POS,1 FROM GEL2D WHERE ID_EXPERIMENT=$id";
		$query{'SAMPLE'}="SELECT ID_SAMPLE,NAME,DES,DISPLAY_POS,1 FROM SAMPLE WHERE ID_EXPERIMENT=$id AND ID_SPOT IS NULL";
		#$query{'DESIGN'}="SELECT ID_DESIGN,NAME,DES,0,1 FROM DESIGN WHERE ID_EXPERIMENT=$id";
	}
	elsif ($item eq 'PROJECT'){
		@children=('EXPERIMENT');
		$query{'EXPERIMENT'}="SELECT ID_EXPERIMENT,NAME,DES,DISPLAY_POS,1 FROM EXPERIMENT WHERE ID_PROJECT=$id";
	}
	elsif ($item eq 'CLASSIFICATION') {
		@children=('CATEGORY');
		$query{'CATEGORY'}="SELECT ID_CATEGORY,NAME,DES,DISPLAY_POS,1 FROM CATEGORY WHERE ID_CLASSIFICATION=$id";
	}
	#elsif ($item eq 'DESIGN') {
	#	@children=('EXPCONDITION');
	#	$query{'EXPCONDITION'}="SELECT ID_EXPCONDITION,NAME,DES,2,1 FROM EXPCONDITION WHERE ID_DESIGN=$id";
	#}
	foreach my $child (@children) {
		my $sth=$dbh->prepare($query{$child});
		$sth->execute;
		my %list;
		my $virtualPos=0;
		while (my ($childID,$name,$des,$displayPos,$validStatus)=$sth->fetchrow_array) {
			($name,$des)=&chkDef($name,$des);
			$displayPos=++$virtualPos if $child eq 'SPOT';
			@{$list{$childID}}=($name,$des,$displayPos,$validStatus);
		}
		$sth->finish;
		push @childrenList,{'ITEM'=>$child,'TYPE'=>$type{$child},'LIST'=>\%list} if scalar keys %list;
	}
	return @childrenList;
}

#############################################################################
####<<< Update the multi_ana directory when an analysis is suppressed >>>####
#############################################################################
sub removeFromMultiAna {
	my ($anaID,$projectID,$anaFile)=@_;
	my %promsPath=&promsConfig::getServerInfo('no_user');
	return unless -e "$promsPath{valid}/multi_ana/proj_$projectID"; # just to be safe
	my $originalFile="$promsPath{valid}/multi_ana/proj_$projectID/$anaFile";
	my $copyFile="$promsPath{valid}/multi_ana/proj_$projectID/copy_$anaFile";
	my $copyLines=0;
	if (-e $originalFile) { # just to be safe
		open (FILE, $originalFile);
		open (FILE2, ">$copyFile");
		while (my $line=<FILE>) {
			chomp $line;
			next if $line==$anaID;
			print FILE2 "$line\n";
			$copyLines++;
		}
		close FILE;
		close FILE2;
		move($copyFile,$originalFile);
	}
	# If there is no more analysis associated to the MSF file, it has to be erased from the multi_ana directory
	if ($copyLines==0) { # msf file no longer used
		unlink $originalFile if -e $originalFile;
		(my $msfFile = $anaFile)=~s/\.ana/\.msf/;
		unlink "$promsPath{valid}/multi_ana/proj_$projectID/$msfFile" if -e "$promsPath{valid}/multi_ana/proj_$projectID/$msfFile"; # true only if called by deleteProjectItem.cgi
		# If there is no more MSF files associated to the analysis, the directory becomes useless and has to be deleted
		opendir my($dir), "$promsPath{valid}/multi_ana/proj_$projectID/";
		my @files = grep { -f "$promsPath{valid}/multi_ana/proj_$projectID/$_" } readdir $dir;
		if (@files==0){
			#remove_tree("$promsPath{valid}/multi_ana/proj_$projectID");
			rmtree("$promsPath{valid}/multi_ana/proj_$projectID");
		}
	}
}


###################################
####<<<Link class subroutine>>>#### DEPRECATED use &getProteinClass
###################################
#sub getLinkClass {
#	my ($confidence,$visibility)=@_;
#	my ($class,$titleString);
#	if ($confidence == 0) {
#		$titleString="This protein was identified based on ion profile comparison with other analyses.<BR>";
#		if ($visibility>0) {$class='bi';} # Bold+Italic
#		else {$class='li';} # Light+Italic
#	}
#	elsif ($confidence == 1) {
#		$titleString="Massists have assigned a <B>bad confidence level</B> to this protein.<BR>";
#		if ($visibility>0) {$class='bg';} # Bold+Gray
#		else {$class='lg';} # Light+Gray
#	}
#	else {
#		$titleString='';
#		if ($visibility>0) {$class='bb';} # Bold+Black
#		else {$class='lb';} # Light+Black
#	}
#	$titleString.='Click to view detailed information on protein.';
#	return ($class,$titleString);
#}

#############################################
####<<<Protein status class subroutine>>>####
#############################################
sub getProteinClass {
	my ($confidence,$visibility)=@_;
	my ($class,$titleString);
	if ($confidence == 0) {
		$titleString="<U>Recovered protein:</U> peptides retrieved by ion-profile comparison with other Analyses during quantification step.<BR>";
		$class='virtualProt'; #'noConf';
	}
	elsif ($confidence == 1) {
		$titleString="Massists have assigned a <B>bad confidence level</B> to this protein.<BR>";
		$class='lowConf';
	}
	else {
		$titleString='';
		$class='highConf';
	}
	$class.=($visibility)? ' visibleProt' : ' hiddenProt';
	$titleString.='Click to view detailed information on protein.';
	return ($class,$titleString);
}

#############################################
####<<<Printing the sequence alignment>>>####
#############################################
#sub printAlignFile{
#	my ($dbh,$projectID,$refAlignFile,$refClass)=@_;
#	my $offset=16;
#	my (%aliasAlign,%seqAlign,%listAlias,%listIdent);
#	my $count=0;
#	my $sth=$dbh->prepare("SELECT ALIAS FROM PROTEIN WHERE IDENTIFIER=? AND ID_PROJECT=$projectID");
#	foreach my $line (@{$refAlignFile}){
#		$count++;
## 		next unless $count>6;
#		next unless $count>2;
#		if ($line=~/^\S/){
#			my ($ident,$seq)=split(/\s+/,$line);
#			if ($listAlias{$ident}) {
#				$aliasAlign{$count}=$listAlias{$ident};
#			}
#			else {
#				$sth->execute($ident);
#				($aliasAlign{$count})=$sth->fetchrow_array;
#			}
#			$seqAlign{$count}=$seq;
#			$listIdent{$count}=$ident;
#		}
#		else {
#			if ($line=~/\S/){
#				#$offset+6;
#				$aliasAlign{$count}='&nbsp';
#				$line=~s/\s/_/gi;
#				$seqAlign{$count}=substr($line,$offset);
#				$seqAlign{$count}=~s/_/&nbsp;/gi;
#			}
#			else{
#				$aliasAlign{$count}='&nbsp';
#				$seqAlign{$count}='&nbsp';
#			}
#		}
#	}
#	$sth->finish;
#	print "<TABLE border=0>\n";
#	print "<TR><TH class=\"title3\" colspan=2>*=identity&nbsp&nbsp&nbsp;:=strong conservation&nbsp&nbsp&nbsp;.=weak conservation<BR><BR></TH></TR>\n";
#	foreach my $nbline(sort{$a<=>$b} keys %aliasAlign){
#		print '<TR>';
#		if ($listIdent{$nbline}) { # not a blank line
#			my $class=($$refClass{$listIdent{$nbline}})? "$$refClass{$listIdent{$nbline}}" : 'ib';
#			print "<TH align=right><FONT class=\"$class\">$aliasAlign{$nbline}</FONT></TH>";
#		}
#		else {print "<TD>$aliasAlign{$nbline}</TD>";} # blank line
#		print "<TD class=\"seq\">&nbsp;$seqAlign{$nbline}</TD>";
#		print "</TR>\n";
#	}
#	print "</TABLE>\n";
#}


##############################################
####<<<Fetching list of classifications>>>####
##############################################
sub getListClass {
	my ($dbh,$projectID)=@_;
	my %classList;
	my $sth=$dbh->prepare("SELECT ID_CLASSIFICATION,NAME,DES FROM CLASSIFICATION WHERE ID_PROJECT=$projectID"); # AND ID_CLASSIFICATION>0");
	$sth->execute();
	while (my ($id_class,$nameClass,$desClass)=$sth->fetchrow_array) {
		$desClass='' unless $desClass;
		@{$classList{$id_class}}=($nameClass,$desClass);
	}
	$sth->finish;
	return (%classList);
}

###########################################################
####<<<Applying Project-wide protein visibility rule>>>####
###########################################################
sub applyVisibilityRule {
	my ($dbh,$projectID,$protVisibility,$oldProtVisibility)=@_;
	($protVisibility)=$dbh->selectrow_array("SELECT PROT_VISIBILITY FROM PROJECT WHERE ID_PROJECT=$projectID") unless defined($protVisibility);
	$oldProtVisibility=$protVisibility unless defined($oldProtVisibility);

	##<List of analyses with quantif or GO data
	my %anaHasQuantifOrGO;
	my @sthAnaList=(
				$dbh->prepare("SELECT DISTINCT AQ.ID_ANALYSIS FROM ANA_QUANTIFICATION AQ,QUANTIFICATION Q,ANALYSIS_PROTEIN AP,PROTEIN P WHERE AQ.ID_QUANTIFICATION=Q.ID_QUANTIFICATION AND Q.FOCUS='protein' AND AQ.ID_ANALYSIS=AP.ID_ANALYSIS AND AP.ID_PROTEIN=P.ID_PROTEIN AND P.ID_PROJECT=$projectID"),
				$dbh->prepare("SELECT GA.ID_ANALYSIS FROM GOANA_ANALYSIS GA,GO_ANALYSIS G,EXPERIMENT E WHERE E.ID_EXPERIMENT=G.ID_EXPERIMENT AND G.ID_GOANALYSIS=GA.ID_GOANALYSIS AND ID_PROJECT=$projectID")
			   );
	foreach my $sth (@sthAnaList) {
		$sth->execute;
		while (my ($anaID)=$sth->fetchrow_array) {
			$anaHasQuantifOrGO{$anaID}=1;
		}
		$sth->finish;
	}
	my $anaListString=(scalar keys %anaHasQuantifOrGO)? 'AND ID_ANALYSIS NOT IN ('.join(',',keys %anaHasQuantifOrGO).')' : '';

	##<A protein is Visible everywhere if Alias or made Visible in at least 1 Match Group.
	if ($protVisibility==2) {
		my $sthVis=$dbh->prepare("SELECT DISTINCT(PROTEIN.ID_PROTEIN) FROM PROTEIN,ANALYSIS_PROTEIN WHERE VISIBILITY>0 AND PROTEIN.ID_PROTEIN=ANALYSIS_PROTEIN.ID_PROTEIN AND ID_PROJECT=$projectID");
		my $sthUp=$dbh->prepare("UPDATE ANALYSIS_PROTEIN SET VISIBILITY=1 WHERE ID_PROTEIN=? AND VISIBILITY=0 $anaListString");
		$sthVis->execute;
		while (my ($protID)=$sthVis->fetchrow_array) {
			$sthUp->execute($protID);
		}
		$sthVis->finish;
		$sthUp->finish;
	}
	##<A protein is Visible everywhere if Alias of at least 1 Match Group.
	elsif ($protVisibility==1) {
		$dbh->do("UPDATE ANALYSIS_PROTEIN SET VISIBILITY=0 WHERE VISIBILITY=1 AND ID_PROTEIN IN (SELECT DISTINCT ID_PROTEIN FROM PROTEIN WHERE ID_PROJECT=$projectID) $anaListString"); # reset all non-alias to hidden
		my $sthVis=$dbh->prepare("SELECT DISTINCT(PROTEIN.ID_PROTEIN) FROM PROTEIN,ANALYSIS_PROTEIN WHERE VISIBILITY=2 AND PROTEIN.ID_PROTEIN=ANALYSIS_PROTEIN.ID_PROTEIN AND ID_PROJECT=$projectID");
		my $sthUp=$dbh->prepare("UPDATE ANALYSIS_PROTEIN SET VISIBILITY=1 WHERE ID_PROTEIN=? AND VISIBILITY=0 $anaListString");
		$sthVis->execute;
		while (my ($protID)=$sthVis->fetchrow_array) {
			$sthUp->execute($protID);
		}
		$sthVis->finish;
		$sthUp->finish;
		&removeFromClassifications($dbh,$projectID);
	}
	##<A protein is Visible only when Alias of a Match Group.
	else { # ==0
		$dbh->do("UPDATE ANALYSIS_PROTEIN SET VISIBILITY=0 WHERE VISIBILITY=1 AND ID_PROTEIN IN (SELECT DISTINCT ID_PROTEIN FROM PROTEIN WHERE ID_PROJECT=$projectID) $anaListString");
		&removeFromClassifications($dbh,$projectID);
	}
}
#sub applyVisibilityRule_old {
#	my ($dbh,$projectID,$protVisibility,$oldProtVisibility)=@_;
#	($protVisibility)=$dbh->selectrow_array("SELECT PROT_VISIBILITY FROM PROJECT WHERE ID_PROJECT=$projectID") unless defined($protVisibility);
#	$oldProtVisibility=$protVisibility unless defined($oldProtVisibility);
#	if ($protVisibility==2) {
#		my $sthVis=$dbh->prepare("SELECT DISTINCT(PROTEIN.ID_PROTEIN) FROM PROTEIN,ANALYSIS_PROTEIN WHERE VISIBILITY>0 AND PROTEIN.ID_PROTEIN=ANALYSIS_PROTEIN.ID_PROTEIN AND ID_PROJECT=$projectID");
#		my $sthUp=$dbh->prepare("UPDATE ANALYSIS_PROTEIN SET VISIBILITY=1 WHERE ID_PROTEIN=? AND VISIBILITY=0");
#		$sthVis->execute;
#		while (my ($protID)=$sthVis->fetchrow_array) {
#			$sthUp->execute($protID);
#		}
#		$sthVis->finish;
#		$sthUp->finish;
#	}
#	elsif ($protVisibility==1) {
#		$dbh->do("UPDATE ANALYSIS_PROTEIN SET VISIBILITY=0 WHERE VISIBILITY=1 AND ID_PROTEIN IN (SELECT DISTINCT ID_PROTEIN FROM PROTEIN WHERE ID_PROJECT=$projectID)"); # rese
#		my $sthVis=$dbh->prepare("SELECT DISTINCT(PROTEIN.ID_PROTEIN) FROM PROTEIN,ANALYSIS_PROTEIN WHERE VISIBILITY=2 AND PROTEIN.ID_PROTEIN=ANALYSIS_PROTEIN.ID_PROTEIN AND ID_PROJECT=$projectID");
#		my $sthUp=$dbh->prepare("UPDATE ANALYSIS_PROTEIN SET VISIBILITY=1 WHERE ID_PROTEIN=? AND VISIBILITY=0");
#		$sthVis->execute;
#		while (my ($protID)=$sthVis->fetchrow_array) {
#			$sthUp->execute($protID);
#		}
#		$sthVis->finish;
#		$sthUp->finish;
#		&removeFromClassifications($dbh,$projectID);
#	}
#	else { # ==0
#		$dbh->do("UPDATE ANALYSIS_PROTEIN SET VISIBILITY=0 WHERE VISIBILITY=1 AND ID_PROTEIN IN (SELECT DISTINCT ID_PROTEIN FROM PROTEIN WHERE ID_PROJECT=$projectID)");
#		&removeFromClassifications($dbh,$projectID);
#	}
#}


################################################################
####<<<Removing non-visible proteins from classifications>>>####
################################################################
sub removeFromClassifications {
	my ($dbh,$projectID)=@_;
	my $sthCatVis=$dbh->prepare("SELECT P.ID_PROTEIN FROM CATEGORY_PROTEIN CP
									INNER JOIN PROTEIN P ON CP.ID_PROTEIN=P.ID_PROTEIN
									INNER JOIN ANALYSIS_PROTEIN AP ON CP.ID_PROTEIN=AP.ID_PROTEIN
									WHERE ID_PROJECT=$projectID
									GROUP BY AP.ID_PROTEIN HAVING MAX(AP.VISIBILITY)=0");
	my $sthDelProt1=$dbh->prepare("DELETE MS FROM MODIFICATION_SITE MS INNER JOIN CATEGORY_PROTEIN CP ON MS.ID_CATEGORY_PROTEIN=CP.ID_CATEGORY_PROTEIN AND CP.ID_PROTEIN=?");
	my $sthDelProt2=$dbh->prepare("DELETE FROM CATEGORY_PROTEIN WHERE ID_PROTEIN=?");
	$sthCatVis->execute;
	foreach (my ($protID)=$sthCatVis->fetchrow_array) {
		$sthDelProt1->execute($protID); # Sites (if any)
		$sthDelProt2->execute($protID); # Proteins
	}
	$sthCatVis->finish;
	$sthDelProt1->finish;
	$sthDelProt2->finish;
}
#sub removeFromClassifications_old {
#	my ($dbh,$projectID)=@_;
#	my $sthCat=$dbh->prepare("SELECT DISTINCT(PROTEIN.ID_PROTEIN) FROM CATEGORY_PROTEIN,PROTEIN WHERE PROTEIN.ID_PROTEIN=CATEGORY_PROTEIN.ID_PROTEIN AND ID_PROJECT=$projectID");
#	my $sthVis=$dbh->prepare("SELECT MAX(VISIBILITY) FROM ANALYSIS_PROTEIN WHERE ID_PROTEIN=?");
#	my $sthDel=$dbh->prepare("DELETE FROM CATEGORY_PROTEIN WHERE ID_PROTEIN=?");
#	$sthCat->execute;
#	foreach (my ($protID)=$sthCat->fetchrow_array) {
#		my ($bestVis)=$sthVis->execute($protID);
#		$sthDel->execute($protID) unless $bestVis;
#	}
#	$sthCat->finish;
#	$sthVis->finish;
#	$sthDel->finish;
#}


######################################
####<<<Compute peptide coverage>>>####
######################################
sub getCoverage {
	my ($protLength,$refStatus)=@_;
	my %status=%{$refStatus};
	my $coverage=0;
	my $prevPos;
	my $matchBeg=0;
	foreach my $pos (sort {$a<=>$b} keys %status) {
		unless ($prevPos) {$prevPos=$matchBeg=$pos; next;}
		$status{$pos}+=$status{$prevPos}; # summing up all status
		if ($matchBeg==0 && $status{$pos}>=1) { # match begins
			$matchBeg=$pos;
		}
		elsif ($matchBeg>0 && $status{$pos}==0) { # match ends
			$coverage+=($pos-$matchBeg+1);
			$matchBeg=0;
		}
		$prevPos=$pos;
	}
	$coverage=($coverage*100)/$protLength;
	return $coverage;
}

###############################################
####<<<Define variable as '' if !defined>>>####
###############################################
sub chkDef {
	my @values=@_;
	for (my $v=0;$v<=$#values;$v++) {
		$values[$v]='' if !defined $values[$v];
	}
	return @values;
}

##########################################################################
####<<<Print the page menu for protein/query list during validation>>>####
##########################################################################
sub printPageMenu {
	my ($page,$maxPage)=@_;
	return if $maxPage==1;
	my $prevDisab=($page==1)? 'disabled' : '';
	my $nextDisab=($page==$maxPage)? 'disabled' : '';
	print qq
|<TABLE border=0 cellspacing=0><TR><TH nowrap>Page
<SELECT onchange="changePage(this.value)" style="font-weight:bold">
|;
	foreach my $p (1..$maxPage) {
		print "<OPTION value='$p'";
		print ' selected' if $p==$page;
		print ">$p</OPTION>\n";
	}
	print qq
|</SELECT>/$maxPage
<INPUT type="button" style="width:25px" value="<" onclick="changePage($page-1)" $prevDisab/><INPUT type="button" style="width:25px" value=">" onclick="changePage($page+1)" $nextDisab/>
</TH></TR></TABLE>
|;
}

##################################################################
####<<<Extract zipped archive containing a (set of) file(s)>>>####
##################################################################
sub unzipArchive {
	my ($file, $dest, $refOpt) = @_; # archive full path, destination dir, options (see below)
	die 'Need a file argument' unless defined $file;
	$refOpt={} unless $refOpt;
	$refOpt->{mode}='silent' unless $refOpt->{mode}; # silent of verbose
	$refOpt->{txtBefore}='' unless $refOpt->{txtBefore}; # if verbose mode: this text is added before extracted file name
	$refOpt->{txtAfter}='' unless $refOpt->{txtAfter}; # if verbose mode: this text is added after extracted file name
	$refOpt->{move2Top}=(!$refOpt->{move2Top} || $refOpt->{move2Top} eq 'no')? 0 : 1; # (do not) move all files to parent ($dest) directory (no sub dir limit)
	$refOpt->{keepArchive}=(!$refOpt->{keepArchive} || $refOpt->{keepArchive} eq 'no')? 0 : 1; # (do not) delete archive file after extraction
	$dest = '.' unless defined $dest;
	my $u = IO::Uncompress::Unzip->new($file) or die "Cannot open $file: $UnzipError";
	my $status;
	my $prevName=' ';
	for ($status = 1; $status > 0; $status = $u->nextStream()) {
		my $header = $u->getHeaderInfo();
		my (undef, $path, $name) = splitpath($header->{Name});
		last if $name eq $prevName; # seems that last entry is looped
		$prevName=$name;
		$path=~s/\/$//; # remove trailing /
		my $destdir;
		if ($refOpt->{move2Top}) {
			$destdir = $dest;
		}
		else {
			$destdir = "$dest/$path";
				unless (-d $destdir) {
					mkpath($destdir) or die "Couldn't mkdir $destdir: $!";
				}
		}
		if (!$name) { # directory # $name =~ m!/?$! ||
			last if $status < 0;
			next;
		}
		print $refOpt->{txtBefore},"$path/$name..." if $refOpt->{mode} eq 'verbose';
		my $destfile = "$destdir/$name";
		my $buff;
		my $fh = IO::File->new($destfile, "w") or die "Couldn't write to $destfile: $!";
		my $count=0;
		while (($status = $u->read($buff)) > 0) {
			$fh->write($buff);
			$count++;
			if ($count==10000) {
				print '.' if $refOpt->{mode} eq 'verbose';
				$count=0;
			}
		}
		$fh->close();
		my $stored_time = $header->{'Time'};
		utime ($stored_time, $stored_time, $destfile) or die "Couldn't touch $destfile: $!";
		print " Done.",$refOpt->{txtAfter} if $refOpt->{mode} eq 'verbose';
	}
	die "Error processing $file: $!\n" if $status < 0;

	##<Move all files to parent directory if any subdirectory>##
	#if ($refOpt->{move2Top}) {
	#	my $existSubDir=1;
	#	while ($existSubDir) {
	#		$existSubDir=0; # reset
	#		opendir (DIR,$dest);
	#		while (my $item = readdir (DIR)) {
	#			next if $item=~/^\.+$/; #  skip '.' & '..' directories
	#			if (-d "$dest/$item") { # directory
	#				opendir (SUBDIR,"$dest/$item");
	#				while (my $subItem = readdir (SUBDIR)) {
	#					next if $subItem=~/^\.+$/; #  skip '.' & '..' directories
	#					move("$dest/$item/$subItem","$dest/$subItem");
	#					$existSubDir=1 if -d "$dest/$subItem"; # set to 1 to check for new subdir moved up
	#				}
	#				close SUBDIR;
	#				rmtree "$dest/$item"; # should be empty but rmtree safer than rmdir
	#			}
	#		}
	#		close DIR;
	#	}
	#}
	unlink $file unless $refOpt->{keepArchive};
	return;
}

########################################################
####<<<Extract valid/spectrum data from data file>>>####
########################################################
sub extractData {
	my ($refSecField,$dataFile,$subFile,$fileType,$msfFile)=@_; # msfFile only for Prot Discov

	###<Loading data file>###
	open (IN,$dataFile);
	my @fileData=<IN>;
	close IN;

	###<Parsing data>###
	open (OUT,">$subFile") || die "Cannot open $subFile";

	if ($fileType=~/\.(DAT|PDM)\Z/) {
		for my $i (0..3) {
			print OUT $fileData[$i];
		}
		my %indexCount;
		my $section='none';
		my $prevSection;
		my $index=4;
		my $i=4;
		while ($i <= $#fileData) {
			if ($fileData[$i]=~/name="(\w+)"/) {
				$prevSection=$section;
				$section=$1;

				##<Writing end of previous section if all lines were not written
				if ($refSecField->{$prevSection} && $refSecField->{$prevSection} ne 'all') {
					print OUT "--gc0p4Jq0M2Yt08jU534c0p\n";
					$index++;
				}
				if ($refSecField->{$section}) {
					print OUT $fileData[$i];		# current ligne
					print OUT $fileData[$i+1];		# following empty line
					$indexCount{$section}=$index;	# previous ligne: --gc0p4Jq0M2Yt08jU534c0p
					$index++;
					$index++;
					$i++;
				}
			}
			elsif ($refSecField->{$section} && ($refSecField->{$section} eq 'all' || ($fileData[$i]=~/^(q.+=)/ && $refSecField->{$section}{$1}))) {
				$fileData[$i]=~s/;\S*/;"AC":0:1:1:1/ if $section eq 'peptides';
				print OUT $fileData[$i];
				$index++;
			}
			$i++;
		}
		if ($fileType =~ /\.PDM/) { # Need to rewrite the spectrum
			my %promsPath=&promsConfig::getServerInfo;
			my ($parentFile,$refIons,$retTimeMin,$spectrumIDXML,$scanNumber,$charge)=&extractSpectrumMSF($msfFile,$refSecField->{extspectrumid});
			my $queryString="query";
			foreach my $section (keys %{$refSecField}) {
				$queryString=$section if $section =~ /query/;
			}
			my $retTimeSec=$retTimeMin*60;
			my ($i,$ions)=(0,"");
			my ($miny,$maxy,$minx,$maxx);
			foreach my $valeur (@{$refIons}){#copy all the ions for this specific query
				if (!$miny){
					$miny=$valeur->{'Y'};
					$maxy=$valeur->{'Y'};
					$minx=$valeur->{'X'};
					$maxx=$valeur->{'X'};
				}
				$ions=$ions.",$valeur->{'X'}:$valeur->{'Y'}";
				if($miny > $valeur->{'Y'}){
					$miny=$valeur->{'Y'};
				}
				if($maxy < $valeur->{'Y'}){
					$maxy=$valeur->{'Y'};
				}
				if($minx > $valeur->{'X'}){
					$minx=$valeur->{'X'};
				}
				if($maxx < $valeur->{'X'}){
					$maxx=$valeur->{'X'};
				}
				$i=$i+1;
			}
			$ions = substr $ions, 1;#remove the first comma
			print OUT "Content-Type: application/x-Mascot; name=\"$queryString\"\n\n";
			print OUT "title=Spectrum$spectrumIDXML%20scans%3a$scanNumber%2c\n";
			print OUT "rtinseconds=$retTimeSec\n";
			print OUT "charge=$charge+\n";
			print OUT "mass_min=$minx\n";
			print OUT "mass_max=$maxy\n";
			print OUT "int_min=$miny\n";
			print OUT "int_max=$maxy\n";
			print OUT "num_vals=$i\n";
			print OUT "num_used1=-1\n";
			print OUT "Ions1=$ions\n";
			print OUT "--gc0p4Jq0M2Yt08jU534c0p--\n";
		}
		print OUT "Content-Type: application/x-Mascot; name=\"index\"\n\n";
		foreach my $ind (sort {$indexCount{$a}<=>$indexCount{$b}} keys %indexCount) {
			print OUT "$ind=$indexCount{$ind}\n";
		}
		print OUT "--gc0p4Jq0M2Yt08jU534c0p--\n";
	}
	elsif ($fileType=~/XML\Z/) {
		if ($fileType eq 'PARAGON.XML') {
			###> 1- Copy all information of PGF old file where there are VARMODIFICATIONS info
			my $pgfFile=$dataFile;
			($pgfFile)=~s/\.xml/\.pgf/;
			my $onSearch=0;
			my %search;
			open (IN,$pgfFile);
			while (my $line=<IN>) {
				print OUT $line;
				chomp($line);
				if ($line =~ '--gc0p4Jq0M2Yt08jU534c0p' ){
					$onSearch=0;
				}
				elsif ($line =~ /Content-Type: Paragon; name=\"search\"/){
					$onSearch=1;
					next;
				}
				if ($onSearch && $line !~ 'Content' ){
					my ($searchID,$rawfileName)=split(/,/,$line);
					$search{$searchID}=$rawfileName;
				}
			}
			close IN;
			###> 2- Create spectrum
			my $onQuery=0;
			my $onMSMSpeaks=0;
			my $nbMatch=0;
			my ($extSpectrumID,$searchRank,$queryString)=($refSecField->{extspectrumid},$refSecField->{'searchRank'},"query$refSecField->{'queryNum'}");
			my ($i,$ions)=(0,"");
			my ($minInt,$maxInt,$minMass,$maxMass,$searchID,$charge);
			foreach my $line (@fileData) {
				chomp($line);
				$onQuery=1 if ($line =~ /<SPECTRUM .* xml:id="$extSpectrumID/ ) ;
				next unless $onQuery;
				if ($onMSMSpeaks) {
					my $peakInfo=$line;
					$peakInfo =~ s/<\!\[CDATA\[//g;
					$peakInfo =~ s/\]\]>//g;
					my ($mz,$chargeState,$intens)=split(/\t/,$peakInfo);
					$ions.=",$mz:$intens";
					($minInt,$maxInt,$minMass,$maxMass)=($intens,$intens,$mz,$mz)if (!$minInt);
					$minInt=$intens if $intens<$minInt;
					$maxInt=$intens if $intens>$maxInt;
					$minMass=$mz if $mz<$minMass;
					$maxMass=$mz if $mz>$maxMass;
					$i++;
				}
				$onMSMSpeaks=1 if ($line =~ /<MSMSPEAKS/);
				if ($line =~ /<MATCH/){
					$nbMatch++;
					print OUT "$line\n";
					if ($nbMatch == $searchRank) {
						($searchID)=($line =~ /<MATCH .* searches=\"(.+)\" seq=/);
						($charge)=($line =~ /<MATCH charge=\"(\d)\"/);
					}
				}
				last if ($line=~ /]]>/);
			}
			$ions = substr $ions, 1;#remove the first comma
			print OUT "Content-Type: application/x-ParagonPeakList; name=\"$queryString\"\n\n";
			print OUT "title=$search{$searchID} (spectrum number in ProteinPilot software: $extSpectrumID)\n";
			#print OUT "rtinseconds=$retTimeSec\n";
			print OUT "charge=$charge+\n";
			print OUT "mass_min=$minMass\n";
			print OUT "mass_max=$maxMass\n";
			print OUT "int_min=$minInt\n";
			print OUT "int_max=$maxInt\n";
			print OUT "num_vals=$i\n";
			print OUT "num_used1=-1\n";
			print OUT "Ions1=$ions\n";
			print OUT "--gc0p4Jq0M2Yt08jU534c0p--\n";
		}
		else {
			my $section='none';
			my $prevSection;
			foreach my $line (@fileData) {
				if ($line=~/name="(\w+)"/) {
					$prevSection=$section;
					$section=$1;
				}
				if ($refSecField->{$section} && $refSecField->{$section} eq 'all') {
					print OUT $line;
				}
			}
		}
	}
	close OUT;
}

######################################################################################
####<Extract a spectrum form a MSF file and returns a string representation of it>####
######################################################################################
sub extractSpectrumMSF { # requires IO::Uncompress::Unzip & XML::Simple
	my ($msfHandle,$spectrumID)=@_; #,$extractDir,$userID
	#unless ($extractDir) {
	#	my %promsPath=&promsConfig::getServerInfo;
	#	$extractDir=$promsPath{'tmp'};
	#}
	#$userID=$ENV{'REMOTE_USER'} unless $userID;
	my $dbsqlite;
	if (ref($msfHandle)){$dbsqlite=$msfHandle;} # $msfHandle is a dbi connection
	else {$dbsqlite = DBI->connect("dbi:SQLite:$msfHandle","","",{PrintError=>1,RaiseError=>1});} # $msfHandle is a path to an msf file

	#########################################################################################
	###<Extraction of the spectrum file contained in the db and conversion into a zipfile>###
	#########################################################################################
#print "'SELECT Spectrum FROM SpectrumHeaders, MassPeaks, Spectra WHERE Spectra.UniqueSpectrumID=SpectrumHeaders.UniqueSpectrumID AND SpectrumHeaders.MassPeakID=MassPeaks.MassPeakID AND SpectrumHeaders.SpectrumID=$spectrumID'<BR>\n";
	my ($proteomeDiscovererVersion)=$dbsqlite->selectrow_array("SELECT SoftwareVersion FROM SchemaInfo ORDER BY rowid ASC LIMIT 1"); # last record
	$proteomeDiscovererVersion=~s/^(\d\.\d+).*/$1/; # x.x.x.xx -> x.x (numerical value)
	my ($fileInfoTabName)=($proteomeDiscovererVersion >= 2.0)?"WorkflowInputFiles":"FileInfos";
	my ($queryString)=($proteomeDiscovererVersion >= 2.2)?"SELECT Spectrum, SpectrumFileName FROM MSnSpectrumInfo MS, MassSpectrumItems MSI WHERE MS.SpectrumID=MSI.ID AND MSI.ID=$spectrumID":"SELECT Spectrum, FileName FROM SpectrumHeaders, MassPeaks, Spectra, $fileInfoTabName WHERE Spectra.UniqueSpectrumID=SpectrumHeaders.UniqueSpectrumID AND SpectrumHeaders.MassPeakID=MassPeaks.MassPeakID AND $fileInfoTabName.FileID=MassPeaks.FileID AND SpectrumHeaders.SpectrumID=$spectrumID";
	my ($zipSpectrumData,$FileName)=$dbsqlite->selectrow_array($queryString);
	$dbsqlite->disconnect unless ref($msfHandle);
	##my $zipFile="$extractDir/spectrum_$userID.zip";
	##open (OUTPUT, ">$zipFile") || die "Error: Cannot write $zipFile file!\n";
	##binmode OUTPUT;
	##print OUTPUT $zipSpectrumData;
	##close(OUTPUT);
	##my $xmlString=`unzip -p $zipFile`;#copy in the stdout the xml file that is used directly in perl
	##unlink $zipFile;
	my $xmlString;
	my @splitName=split(/[\/\\]/,$FileName);
	unzip \$zipSpectrumData => \$xmlString or die "unzip failed: $UnzipError\n";
	#require XML::Simple;
	my $xml = new XML::Simple();
	my $xmlData = $xml->XMLin($xmlString);

	#######################################################################################
	###<Return the ions, the Retention Time, the spectrumidentifier and the scan number>###
	#######################################################################################
	return ($splitName[-1],$xmlData->{'PeakCentroids'}{'Peak'},$xmlData->{'Header'}{'SpectrumIdentifiers'}{'SpectrumIdentifier'}{'RetentionTime'},$xmlData->{'Header'}{'SpectrumID'},$xmlData->{'Header'}{'SpectrumIdentifiers'}{'SpectrumIdentifier'}{'ScanNumber'},$xmlData->{'PrecursorInfo'}{'Charge'});
}

############################################################
####<Fetching proteins MW and description from databank>####
############################################################
# WARNING: $refProtList does not have the same structure:
# $refProtList->{$protID}{$anaID}=1 or $refProtList->{$protID}{pepSeq}=1 (or $refProtList->{$protID}=1?)
sub getProtInfo { # called by scanDatank.pl, storeAnalysis.cgi (only for MSF) and importMaxQuant
	my ($action,$dbh,$databankID,$refAnaList,$refProtDes,$refProtMW,$refProtOrg,$refProtLength,$refProtSeq,$refProtList,$localDbFile,$prefixID)=@_;
	$prefixID='' unless $prefixID; # for MaxQuant import, dynamic addition of CON__ prefix to ids in contaminant db
	my @percent=(10,20,30,40,50,60,70,80,90,100); # needed to monitor process progression
	my %taxonIDs;
	my %promsPath=&promsConfig::getServerInfo('no_user'); # no_user because used by importAnalysis.pl
	my %mascotServers=&promsConfig::getMascotServers;
	my %massAAave=&promsConfig::getMassAAave; # Average mass, needed for protein mass calculation
	my %massATave=&promsConfig::getMassATave; # Average mass, needed for protein mass calculation
	my $userID=($ENV{'REMOTE_USER'})? $ENV{'REMOTE_USER'} : 'myproms'; # no REMOTE_USER when called from cron
	$refProtSeq = {} if !defined $refProtSeq;
	
	##>Check $refProtList structure (multi vs single analysis)
	my $multiAnalysis=(scalar @{$refAnaList} > 1)? 1 : 0;
	unless ($multiAnalysis) { # can still be a multi-analysis structure but with only 1 analysis!!!
		my $anyProtKey=(keys %{$refProtList})[0];
		if (ref $refProtList->{$anyProtKey}) {
			my $any2ndKey=(keys %{$refProtList->{$anyProtKey}})[0];
			$multiAnalysis=1 if $any2ndKey=~/^\d+$/; # $anaID vs $pepSeq
		}
	}

	my $analysisID=$refAnaList->[0]; # for single analysis (verbose mode OR PDM)

	print '<FONT class="title3">&nbsp;-Fetching protein annotations and sequences ' if $action eq 'verbose';
	my ($dbFile,$numEntry,$dbOrganism,$parseRules,$identType,$defIdentType,$isCrap)=$dbh->selectrow_array("SELECT FASTA_FILE,NUM_ENTRY,ORGANISM,PARSE_RULES,IDENTIFIER_TYPE,DEF_IDENT_TYPE,IS_CRAP FROM DATABANK,DATABANK_TYPE WHERE DATABANK.ID_DBTYPE=DATABANK_TYPE.ID_DBTYPE AND ID_DATABANK=$databankID");
	$defIdentType='UNKNOWN' unless $defIdentType;
	$identType=$defIdentType unless $identType;
	my ($mascotServer,$fileInfo,$dbankDir,$databankFile); # LWP connection to Mascot #,$proxyStrg1,$proxyStrg2);

	if ($localDbFile) { # eg. for MSF searches
		$databankFile=$localDbFile;
		$fileInfo='local search file';
		#$numEntry=`grep -c '>' $localDbFile`;
		#chomp ($numEntry);
		$numEntry=0;
		open(DB,$localDbFile);
		while(<DB>) {
			$numEntry++ if /^>/;
		}
		close DB;
		$numEntry=100 unless $numEntry;
	}
	elsif ($dbFile=~/:/) { # Mascot databank
		#require LWP::UserAgent;
		#$agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');

		($mascotServer,$dbankDir,my $dbFileName)=split(':',$dbFile);

		##<Proxy
		#if ($mascotServers{$mascotServer}[1]) {
		#	if ($mascotServers{$mascotServer}[1] eq 'no') {$agent->no_proxy($mascotServers{$mascotServer}[0]);}
		#	else {$agent->proxy('http', $mascotServers{$mascotServer}[1]);}
		#}
		#else {$agent->env_proxy;}
		#
		print "(Checking databank '$mascotServer > $dbankDir'..." if $action eq 'verbose';
		#
		#my $response = $agent->post("$mascotServers{$mascotServer}[0]/cgi/myproms4databanks.pl",
		#					  ['ACT'=>'dbFile',
		#					   'DB'=>$dbankDir,
		#					   'FILE'=>$dbFileName
		#					   ]
		#					  );
		#unless ($response->is_success) {
		#	return $!;
		#}
		#if ($response->content=~/^#Error/) {
		#	(my $errorText=$response->content)=~s/^#//;
		#	return "$errorText at $mascotServer databank '$dbankDir' update step.";
		#}
		##my ($mascotDbFile,$numEnt)=split(/\t/,$response->content);
		#my @responseLines=split(/\n/,$response->content); # myproms4databanks.pl can return multiple lines
		#my ($mascotDbFile,$numEnt)=split(/\t/,$responseLines[-1]);
		#my $newDbFileName=(split(/\//,$mascotDbFile))[-1];
		#if ($newDbFileName ne $dbFileName) {
		#	$dbh->do("UPDATE DATABANK SET VERSION_NAME=NULL,VERSION_DATE=NOW(),UPDATE_DATE=NOW(),FASTA_FILE='$mascotServer:$dbankDir:$newDbFileName',NUM_ENTRY=$numEnt WHERE ID_DATABANK=$databankID");
		#	$numEntry=$numEnt;
		#	$dbFileName=$newDbFileName;
		#}
		my $verbose=($action eq 'verbose')? $action : undef;
		my ($updated,$newDbFileName,$numEnt,$errorText)=&updateMascotDB($databankID,$verbose);
		if ($errorText) {
			return " $errorText at $mascotServer databank '$dbankDir' update step.";
		}
		elsif ($updated) {
			$numEntry=$numEnt;
			$dbFileName=$newDbFileName;
		}
		$fileInfo=$dbFileName;
		#($databankFile=$dbFileName)=~s/.*sequence/$mascotServers{$mascotServer}[2]\/sequence/ if $mascotServers{$mascotServer}[2]; # local access to Mascot banks
		if ($mascotServers{$mascotServer}[2]) { # local access to Mascot banks
			$databankFile="$mascotServers{$mascotServer}[2]/sequence/$dbankDir/current/$dbFileName"; # standard hierarchy
			unless (-e $databankFile) { # non-standard hierarchy or unexpected problem -> try retrieve path via LWP to Mascot
				require LWP::UserAgent;
				my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
				$agent->timeout(360);
				if ($mascotServers{$mascotServer}[1]) { # proxy settings
					if ($mascotServers{$mascotServer}[1] eq 'no') {$agent->no_proxy($mascotServers{$mascotServer}[0]);}
					else {$agent->proxy('http', $mascotServers{$mascotServer}[1]);}
				}
				else {$agent->env_proxy;}
				my $response = $agent->get("$mascotServers{$mascotServer}[0]/cgi/myproms4databanks.pl?ACT=list");
				while (my $wait = $response->header('Retry-After')) {
					sleep $wait;
					$response = $agent->get($response->base);
				}
				my @resultLines;
				if ($response->is_success) {
					@resultLines = split("\n",$response->content);
				}
				foreach (@resultLines) {
					next if /^#/;
					#chomp;
					if (/^(\S+)\s+(\S+)/) {
						if ($dbankDir eq $1) {
							my $fullDbankfile=$2;
							($databankFile=$fullDbankfile)=~s/.*sequence/$mascotServers{$mascotServer}[2]\/sequence/;
							last;
						}
					}
				}
			}
		}
		print ' Done.) ' if $action eq 'verbose';
	}
	else { # myProMS databank
		$databankFile="$promsPath{banks}/db_$databankID/$dbFile";
		$fileInfo=$dbFile;
	}

	print "from $fileInfo ($numEntry entries):</FONT><BR>\n" if $action eq 'verbose';
	#my $masterFastaFile=($multiAnalysis)? "$promsPath{valid}/analysis_$userID.fasta" : "$promsPath{valid}/ana_$analysisID/analysis.fasta";
	my $masterFastaFile="$promsPath{tmp}/analysis_$userID.fasta";
	my $maxProtID=scalar keys (%{$refProtList});

	##<Directory for analysis.fasta file(s)
	my ($anaValidStatus)=$dbh->selectrow_array("SELECT VALID_STATUS FROM ANALYSIS WHERE ID_ANALYSIS=$analysisID");
	my $anaParentDir;
	if (!$anaValidStatus || $anaValidStatus < 2) {$anaParentDir=$promsPath{'valid'};}
	else {
		my $projectID=&promsMod::getProjectID($dbh,$analysisID,'analysis');
		$anaParentDir="$promsPath{peptide}/proj_$projectID";
		mkdir $anaParentDir unless -e $anaParentDir;
	}

	###<Distant mascot server>###
	if ($mascotServer && !$mascotServers{$mascotServer}[2] && !$localDbFile) {
		if ($action eq 'verbose') {
			print qq|<TABLE cellpadding=0><TR><TH align=right class="title3">Mascot scan:</TH><TD><IFRAME name="watchScan_$analysisID" width=800 height=25 frameborder=0 src="$promsPath{cgi}/storeProjectItem.cgi?WATCH=1&ITEM_ID=$analysisID&MAXPROT=$maxProtID&DBSIZE=$numEntry&DIR=$anaParentDir"></IFRAME></TD><TH>|;
		}
		my $protString=join(':',(keys %{$refProtList}));
		#$parseRules=uri_escape($parseRules); # URL encoding unsafe parameters

#print "OUT=$masterFastaFile\nDB=>$dbankDir\nparseRules=>$parseRules\nprotList=>$protString\n";
		require LWP::UserAgent;
		my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
		$agent->timeout(360);
		if ($mascotServers{$mascotServer}[1]) { # proxy settings
			if ($mascotServers{$mascotServer}[1] eq 'no') {$agent->no_proxy($mascotServers{$mascotServer}[0]);}
			else {$agent->proxy('http', $mascotServers{$mascotServer}[1]);}
		}
		else {$agent->env_proxy;}
		my $response=$agent->post("$mascotServers{$mascotServer}[0]/cgi/myproms4databanks.pl",
									['ACT'=>'scan',
									'DB'=>$dbankDir,
									'parseRules'=>$parseRules,
									'identType'=>$identType,
									'protList'=>$protString
									],
									':content_file'=>$masterFastaFile
								   );
		unless ($response->is_success) {
			print "**ERROR: Bad answer from Mascot server!**</TH></TR></TABLE>\n" if $action eq 'verbose';
			return $!;
		}
		if ($response->content=~/^#Error/) {
			(my $errorText=$response->content)=~s/^#//;
			return "$errorText at $mascotServer databank '$dbankDir' scan step.</TH></TR></TABLE>\n";
		}

		if ($action eq 'verbose') {
			my $lwpEnd=0;
			my $count=0;
			my $count2=0;
			while ($lwpEnd==0 && $count2 < 10) { # 5 min
				if ($count > 30) {
					print '°';
					$count=0;
					$count2++;
				}
				sleep 1;
				my $lastLine=`tail -1 $masterFastaFile`;
				$lwpEnd=1 if ($lastLine && $lastLine=~/^##END/);
				$count++;
			}
			print "</TH></TR></TABLE>\n";
			#sleep 1;
			mkdir "$anaParentDir/ana_$analysisID" unless -e "$anaParentDir/ana_$analysisID";
			open(FLAG,">$anaParentDir/ana_$analysisID/end.txt");
			print FLAG '#';
			close FLAG;
		}

		##<Scanning fasta file for protein info
		if (-e $masterFastaFile) { # multi-analysis!!! $refProtList->{$anaID}{$protID}
			open(MFAS,$masterFastaFile);
			my %fileHandles; # file handles stored in variables
			#if ($multiAnalysis) {
				##<Opening all files
				foreach my $anaID (@{$refAnaList}) { # Dispatching protein entries in corresponding fasta files
					mkdir "$anaParentDir/ana_$anaID" unless -e "$anaParentDir/ana_$anaID";
					open ($fileHandles{$anaID},">>$anaParentDir/ana_$anaID/analysis.fasta"); # >> in case of multi-db search (&getProtInfo called for each db)
				}
			#}
			my ($protEntry, %matchedAnaList);
			my $previousAnnotString = '';
			my $sequence = '';
			while(<MFAS>) {
				if (/^>(.+)\n/) {
					my $annotStrg=$1;
					if($previousAnnotString) {
						foreach my $entryStrg (split('',$previousAnnotString)) {
							my ($identifier)=($entryStrg=~/^(\S+) ##/);
							$refProtSeq->{"$prefixID$identifier"}=$sequence;
						}
						$sequence = "";
					}
					$previousAnnotString = $annotStrg;
					
					#if ($multiAnalysis) { # Dispatching previous protein entry in corresponding fasta files
						foreach my $anaID (keys %matchedAnaList) { # if matched anaID => protEntry is set
							print {$fileHandles{$anaID}} $protEntry;
						}
						$protEntry=">$annotStrg\n"; # reset entry with current protein
						%matchedAnaList=();
					#}
					foreach my $entryStrg (split('',$annotStrg)) { # or '\001'
						my ($identifier,$des,$org,$mw,$length)=($entryStrg=~/^(\S+) ##DES=(.+) ##ORG=(.+) ##MW=(\S+) ##LEN=(\d+)/);
						$refProtDes->{"$prefixID$identifier"}=$des;
						$refProtOrg->{"$prefixID$identifier"}=$org;
						$refProtMW->{"$prefixID$identifier"}=$mw;
						$refProtLength->{"$prefixID$identifier"}=$length;
						if ($multiAnalysis) {
							foreach my $anaID (keys %{$refProtList->{"$prefixID$identifier"}}) {
								$matchedAnaList{$anaID}=1; # hash in case same analysis for 2 sub-entries ()
							}
						}
						else {$matchedAnaList{$analysisID}=1;}
					}
				} 
				#elsif ($multiAnalysis && $_ !~ /^#/) {$protEntry.=$_;} # add sequence line(s) except ##END
				# add sequence line(s) except ##END line
				elsif ($_ !~ /^#/) {
					$protEntry.=$_;
					$_=~s/\W+//g; # chomp not always enough? (Windows!) + cleans any non-AA
					$sequence .= $_;
				}
			}
			close MFAS;
			
			if($previousAnnotString && $sequence) {
				foreach my $entryStrg (split('',$previousAnnotString)) {
					my ($identifier)=($entryStrg=~/^(\S+) ##/);
					$refProtSeq->{"$prefixID$identifier"}=$sequence;
				}
			}
			
			#if ($multiAnalysis) { # Dispatching protein entries in corresponding fasta files
				foreach my $anaID (keys %matchedAnaList) { # last protEntry of master file
					print {$fileHandles{$anaID}} $protEntry;
				}
				##<Closing all files
				foreach my $anaID (@{$refAnaList}) {
					close $fileHandles{$anaID};
				}
				unlink $masterFastaFile; # delete master fasta file
			#}
		}
		else {print "**ERROR: No annotation retrieved**\n";}
	}

	###<myProMS databank OR local mascot server OR local DB file (MSF)>###
	elsif (-e $databankFile) {
		my @rules=split(',:,',$parseRules);
		my ($idRule)=($rules[0]=~/ID=(.+)/);
		my ($desRule,$orgRule);
		if ($rules[1]) {
			#if ($rules[1]=~/DES=/) {($desRule)=($rules[1]=~/DES=(.+)/);}
			#else {($orgRule)=($rules[1]=~/ORG=(.+)/);}
			if ($rules[1]=~/DES=/) {($desRule=$rules[1])=~s/DES=//;}
			else {($orgRule=$rules[1])=~s/ORG=//;}
		}
		#if ($rules[2]) {($orgRule)=($rules[2]=~/ORG=(.+)/);}
		if ($rules[2]) {($orgRule=$rules[2])=~s/ORG=//;}

		###<Scanning DBank file
		my $sequence='';
		my $matched;
		my (@matchedID,%matchedAnaList);
		my @limitValue;
		foreach my $pc (@percent) {push @limitValue,int(0.5+($maxProtID*$pc/100));}
		my $index=0;
		my $counter1=0;
		my $counter2=0;
		my $dotString='.';
		my $maxCounter2=int(0.5+($numEntry/100));
		if ($maxCounter2<2) {$maxCounter2=1; $dotString='...';}
		print "<FONT class=\"title3\">0%" if $action eq 'verbose';
		open (DB, $databankFile); # || die "can't open $databankFile\n";
		###open (FAS, ">$promsPath{valid}/ana_$analysisID/analysis.fasta");
		my %fileHandles; # file handles stored in variables
		##<Opening all files
		foreach my $anaID (@{$refAnaList}) { # Dispatching protein entries in corresponding fasta files
			mkdir "$anaParentDir/ana_$anaID" unless -e "$anaParentDir/ana_$anaID";
			open ($fileHandles{$anaID},">>$anaParentDir/ana_$anaID/analysis.fasta"); # >> in case of multi-db search (&getProtInfo called for each db)
		}
		while (<DB>) {
			if (/^>/) {
				if ($matched) {
					my %countAA;
					foreach my $aa (split(//,uc($sequence))) {
						$countAA{$aa}++;
					}
					my $mass=0;
					foreach my $aa (keys %countAA) {
						$mass+=($massAAave{$aa}*$countAA{$aa}) if $massAAave{$aa}; # some characters are not amino-acids
					}
					$mass+=($massATave{H}+$massATave{H}+$massATave{O}) if $mass; # H(Nter) + OH(Cter)
					$mass=sprintf "%.2f",$mass; # no need for more precision
					foreach my $identifier (@matchedID) { # already prefixed if any
						$refProtMW->{$identifier}=$mass;
						$refProtLength->{$identifier}=length($sequence);
						$refProtSeq->{$identifier}=$sequence;
					}
					###<Writing sequence to analysis.fasta
					###print FAS ">",join('',@matchedID),"\n$sequence\n";
					foreach my $anaID (keys %matchedAnaList) {
						print {$fileHandles{$anaID}} ">",join('',@matchedID),"\n$sequence\n";
					}
					@matchedID=();
					%matchedAnaList=();
					$sequence='';
					$matched=0;
					last if $counter1==$maxProtID; # all identifier were matched
				}
				chomp;
				(my $newLine=$_)=~s/^>\s*//;
				my @line=split(//,$newLine);
				foreach my $entry (@line) {
# 1					my ($identifier,$right)=($entry=~/\A(gi\|\d+)\S*\s*(.*)/);
# 2					my ($identifier,$right)=($entry=~/\A([^\|,\s]+\|?[^\,\s|]*)\S*\s*(.*)/); # (1st db|identifier or 1st identifier if no |) (des + org)
# 3					my ($identifier)=($entry=~/^(\w*\|?\w*\.?\d*)/);
					my ($identifier)=($entry=~/$idRule/);
					next unless $identifier; # just to be safe

					$identifier=$prefixID.$identifier; # prefixed from now on (if any)!

					if ($identType eq 'UNIPROT_ID') { # check for isoforms in 1st keyword before | & add to UNIPROT_ID
						$entry=~s/^sp\|//;
						if ($entry=~/^[^|]+-(\d+)\|/) {
							$identifier.="-$1";
						}
					}


					if (defined($refProtList->{$identifier})) {
# 						my ($right)=($entry=~/^\S+\s+(.*)/);
# 						next unless $right;
# 						if ($right=~/\[([^\[\]]+)\]\Z/ && $1 !~ /Contains/) { # no [ or ] allowed inside external []
# 							$protOrg{$identifier}=$1;
# 							(my $tempOrg=$protOrg{$identifier})=~s/([\(\)\|\\\*])/\\$1/g; # inactivate metacharacters
# 							($protDes{$identifier}=$right)=~s/ \[$tempOrg\]//;
# 						}
# 						else {$protDes{$identifier}=$right;}

						my $des; ($des)=($entry=~/$desRule/) if $desRule; # optional
						my $org; ($org)=($entry=~/$orgRule/) if $orgRule; # optional
						if ($org) {
							($des)=($entry=~/^\S+\s+(.+)$orgRule/) unless $desRule;
							$org=~s/\s+\Z//; # removing trailing spaces
							if ($org !~ /\D/) { # taxonID
								unless (defined $taxonIDs{$org}) {
									my $orgLine=`grep '^$org\[^0-9\]' $promsPath{banks}/taxonomy_codes.txt`;
									if ($orgLine) {($taxonIDs{$org})=($orgLine=~/^\d+\t(.+)\t/);}
									else {$taxonIDs{$org}=($dbOrganism)? $dbOrganism : '';}
								}
								$org=$taxonIDs{$org};
							}
							$refProtOrg->{$identifier}=$org;
						}
						else {
							($des)=($entry=~/^\S+\s+(.+)/) unless $des; # $desRule;
							#$refProtOrg->{$identifier}=$dbOrganism if $dbOrganism; # 1 single species in DB
						}
						$des='' unless $des;
						($refProtDes->{$identifier}=$des)=~s/\s+\Z//; # removing trailing spaces

						$matched=1;
						push @matchedID,$identifier;
						if ($multiAnalysis) { # local access to mascot server (called by scanDatank.pl)
							foreach my $anaID (keys %{$refProtList->{$identifier}}) {
								$matchedAnaList{$anaID}=1; # hash in case same analysis for 2 sub-entries ()
							}
						}
						else { # e.g. PDM (called by storeAnalysis.cgi)
							$matchedAnaList{$analysisID}=1;
						}

						##<Counter progress
						$counter1++;
						if ($action eq 'verbose' && $counter1>=$limitValue[$index]) {
							print "$percent[$index]%"; # keeping connection alive
							$index++;
						}
						$counter2++;
						if ($action eq 'verbose' && $counter2==$maxCounter2) {print $dotString; $counter2=0;}
					}
				}
				$counter2++;
				if ($action eq 'verbose' && $counter2==$maxCounter2) {print '.'; $counter2=0;}
			}
			elsif ($matched) {
				#chomp;
				$_=~s/\W+//g; # chomp not always enough? (Windows!) + cleans any non-AA
				$sequence.=$_;
			}
		}

		if ($matched) {
			my %countAA;
			foreach my $aa (split(//,uc($sequence))) {
				$countAA{$aa}++;
			}
			my $mass=0;
			foreach my $aa (keys %countAA) {
				$mass+=($massAAave{$aa}*$countAA{$aa}) if $massAAave{$aa}; # some characters are not amino-acids
			}
			$mass+=($massATave{H}+$massATave{H}+$massATave{O}) if $mass; # H(Nter) + OH(Cter)
			$mass=sprintf "%.2f",$mass; # no need for more precision
			foreach my $identifier (@matchedID) {
				$refProtMW->{"$prefixID$identifier"}=$mass;
				$refProtLength->{"$prefixID$identifier"}=length($sequence);
				$refProtSeq->{$identifier}=$sequence;
			}
			###<Writing sequence to analysis.fasta
			###print FAS '>',join('',@matchedID),"\n$sequence\n";
			foreach my $anaID (keys %matchedAnaList) {
				print {$fileHandles{$anaID}} ">",join('',@matchedID),"\n$sequence\n";
			}
		}
		# end of file
		##<Closing all files
		close DB;
		###close FAS;
		foreach my $anaID (@{$refAnaList}) {
			close $fileHandles{$anaID};
		}

		if ($action eq 'verbose') {
			print " (",$maxProtID-$counter1,"/$maxProtID unmatched proteins).</FONT>" if $counter1<$maxProtID;
			print "</FONT><BR>\n";
		}
	}

	##<Trying to annotate proteins from DB
	my ($sthMaster,$all2accRules,$sthProt);
	if ($identType=~/UNIPROT_/) {
		my $identCode=($identType eq 'UNIPROT_ID')? 'ID' : 'AC';
		my ($identID)=$dbh->selectrow_array("SELECT ID_IDENTIFIER FROM IDENTIFIER WHERE CODE='$identCode'");
		$sthMaster=$dbh->prepare("SELECT 1,PROT_DES,PROT_LENGTH,MW,PROT_SEQ,SCIENTIFIC_NAME FROM MASTERPROT_IDENTIFIER I,MASTER_PROTEIN MP,SPECIES S WHERE I.ID_MASTER_PROTEIN=MP.ID_MASTER_PROTEIN AND MP.ID_SPECIES=S.ID_SPECIES AND I.ID_IDENTIFIER=$identID AND I.VALUE=? ORDER BY UPDATE_DATE DESC LIMIT 1");
		if ($identType eq 'UNIPROT_ALL') {
			my ($accParseRules)=$dbh->selectrow_array("SELECT PARSE_RULES FROM DATABANK_TYPE WHERE DEF_IDENT_TYPE='UNIPROT_ACCESSION' ORDER BY ID_DBTYPE DESC LIMIT 1");
			my @rules=split(',:,',$accParseRules);
			($all2accRules)=($rules[0]=~/ID=(.+)/);
		}
		$sthProt=$dbh->prepare("SELECT IDENTIFIER,PROT_DES,PROT_LENGTH,MW,PROT_SEQ,ORGANISM FROM PROTEIN WHERE IDENTIFIER LIKE ? ORDER BY UPDATE_DATE DESC");
	}
	
	my ($protCount,$numMatched,$verboseAnnot)=(0,0,0);
	foreach my $identifier (keys %{$refProtList}) { # can be overwritten by next call of &getProtInfo if multi-db search
		next if $refProtDes->{$identifier};
		if ($protCount < 500 || $action ne 'verbose') {
			my $match=0;
			#<Check if protein already recorded in DB as master protein
			if ($sthMaster) {
				if ($protCount==0 && $action eq 'verbose') {
					print '<FONT class="title3">&nbsp;-Annotating unmatched proteins from database...';
					$verboseAnnot=1;
				}
				my ($usedIdentifier)=($all2accRules)? ($identifier=~/$all2accRules/) : ($identifier);
				$sthMaster->execute($usedIdentifier);
				($match,$refProtDes->{$identifier},$refProtLength->{$identifier},$refProtMW->{$identifier},$refProtSeq->{$identifier},$refProtOrg->{$identifier})=$sthMaster->fetchrow_array;
			}
			#<Check if protein already recorded in DB as simple protein
			if (!$match && $sthProt) {
				$sthProt->execute("%$identifier%");
				my $qIdent=quotemeta($identifier);
				while (my ($ident,$des,$ln,$mw,$seq,$org)=$sthProt->fetchrow_array) {
					if ($ident=~/(^|\|)$qIdent(\||$)/) {
						($refProtDes->{$identifier},$refProtLength->{$identifier},$refProtMW->{$identifier},$refProtSeq->{$identifier},$refProtOrg->{$identifier})=($des,$ln,$mw,$seq,$org);
						$match=1;
						last;
					}
				}
			}
			$numMatched++ if $match;
		}
		if ($protCount==500 && $verboseAnnot) {
			print '[Aborting: too many unmatched proteins]';
		}
		
		$refProtDes->{$identifier}='no description' unless $refProtDes->{$identifier};
		if ($dbOrganism) {
			$refProtOrg->{$identifier}=$dbOrganism if (!$refProtOrg->{$identifier} || $refProtOrg->{$identifier} eq 'unknown organism');
		}
		elsif (!$refProtOrg->{$identifier}) {
			$refProtOrg->{$identifier}='unknown organism';
		}
		#$refProtOrg->{$identifier}= 'unknown organism' unless $refProtOrg->{$identifier} ;
		$refProtMW->{$identifier}=0 unless $refProtMW->{$identifier};
		$refProtLength->{$identifier}=0 unless $refProtLength->{$identifier};
		
		$protCount++;
		print '.' if ($verboseAnnot && !($protCount % 50));
	}
	$sthMaster->finish if $sthMaster;
	$sthProt->finish if $sthProt;

	print " Done ($numMatched/$protCount additional proteins matched).</FONT><BR>\n" if $verboseAnnot;
	#print "<FONT class=\"title2\">&nbsp;Done.</FONT><BR><BR>\n" if $action eq 'verbose';
}

################################
####<Update Mascot databank>####
################################
sub updateMascotDB {
	my ($databankID,$verbose)=@_;
	my %mascotServers=&promsConfig::getMascotServers;

	my $dbh=&promsConfig::dbConnect('no_user');
	my ($dbFile,$numEntry)=$dbh->selectrow_array("SELECT FASTA_FILE,NUM_ENTRY FROM DATABANK WHERE ID_DATABANK=$databankID");
	$dbh->disconnect;

	my ($mascotServer,$dbankDir,$dbFileName)=split(/:/,$dbFile);
	my ($errorText,$updated,$newDbFileName);
	if ($mascotServers{$mascotServer}) {
		require LWP::UserAgent;
		my $agent = LWP::UserAgent->new(agent=>'libwww-perl myproms@curie.fr');
		$agent->timeout(360);
		#<Proxy
		if ($mascotServers{$mascotServer}[1]) {
			if ($mascotServers{$mascotServer}[1] eq 'no') {$agent->no_proxy($mascotServers{$mascotServer}[0]);}
			else {$agent->proxy('http', $mascotServers{$mascotServer}[1]);}
		}
		else {$agent->env_proxy;}
		my $responseData='';
		my $numChunks=0;
		$agent->request(HTTP::Request->new(GET => "$mascotServers{$mascotServer}[0]/cgi/myproms4databanks.pl?ACT=dbFile&DB=$dbankDir&FILE=$dbFileName"),
								sub {
									my ($chunk, $res) = @_;
									$responseData.=$chunk;
									print '.' if $verbose;
									$numChunks++;
									if ($numChunks > 500) {
										return(-1,'',0,'Update is taking too long');
									}
								});
		if ($responseData) {
			my @responseLines=split(/\n/,$responseData);
			my ($mascotDbFile,$numEnt)=split(/\t/,$responseLines[-1]);
			$newDbFileName=(split(/\//,$mascotDbFile))[-1];
			if ($newDbFileName=~/^#/) {
				$updated=-1;
				$responseData=~s/#Error: //;
				$errorText="Unexpected response from $mascotServer: \"$responseData\"";
			}
			elsif ($newDbFileName ne $dbFileName) {
				$updated=1;
				$numEntry=$numEnt;
				$dbh=&promsConfig::dbConnect('no_user');
				$dbh->do("UPDATE DATABANK SET VERSION_NAME=NULL,VERSION_DATE=NOW(),UPDATE_DATE=NOW(),FASTA_FILE='$mascotServer:$dbankDir:$newDbFileName',NUM_ENTRY=$numEntry WHERE ID_DATABANK=$databankID");
				$dbh->commit;
				$dbh->disconnect;
			}
			else {$updated=0;}
		}
		else {
			$updated=-1;
			$errorText="No response from $mascotServer";
		}
	}
	else {
		$updated=-1;
		$errorText="$mascotServer is not declared in myProMS configuration file.";
	}

	return($updated,$newDbFileName,$numEntry,$errorText);
}

########################################################
####<Deleting unused master proteins & dependencies>####
########################################################
sub deleteUnusedMasterProteins {
	my ($dbh,$refModifMasterProteins)=@_;
	my %modifSpecies;
	my $sthNumMP=$dbh->prepare("SELECT COUNT(ID_PROTEIN) FROM PROTEIN WHERE ID_MASTER_PROTEIN=?");
	my $sthSp=$dbh->prepare("SELECT ID_SPECIES FROM MASTER_PROTEIN WHERE ID_MASTER_PROTEIN=?");
	my $sthDelMI=$dbh->prepare("DELETE FROM MASTERPROT_IDENTIFIER WHERE ID_MASTER_PROTEIN=?");
	my $sthDelMP=$dbh->prepare("DELETE FROM MASTER_PROTEIN WHERE ID_MASTER_PROTEIN=?");
	my $checkProt=1;
	unless ($refModifMasterProteins) {
		$refModifMasterProteins={};
		my $sthAll=$dbh->prepare("SELECT MP.ID_MASTER_PROTEIN,ID_SPECIES FROM MASTER_PROTEIN MP LEFT OUTER JOIN PROTEIN P ON P.ID_MASTER_PROTEIN=MP.ID_MASTER_PROTEIN WHERE P.ID_MASTER_PROTEIN IS NULL");
		$sthAll->execute;
		while (my ($masterProtID,$speciesID)=$sthAll->fetchrow_array) {
			$refModifMasterProteins->{$masterProtID}=1;
			$modifSpecies{$speciesID}=1;
		}
		$sthAll->finish;
		$checkProt=0;
	}
	my $numMasterDel=0;
	foreach my $masterProtID (keys %{$refModifMasterProteins}) {
		$sthNumMP->execute($masterProtID);
		if ($checkProt) {
			my ($masterUsed)=$sthNumMP->fetchrow_array;
			next if $masterUsed; # still referenced by a protein
			$sthSp->execute($masterProtID);
			my ($speciesID)=$sthSp->fetchrow_array;
			$modifSpecies{$speciesID}=1 if $speciesID; # just to be safe
		}
		$sthDelMI->execute($masterProtID);
		$sthDelMP->execute($masterProtID);
		$numMasterDel++;
	}
	$sthNumMP->finish;
	$sthDelMI->finish;
	$sthDelMP->finish;

	#<Deleting unused species
	&deleteUnusedSpecies($dbh,\%modifSpecies);
	return $numMasterDel;
}

#################################
####<Deleting unused species>####
#################################
sub deleteUnusedSpecies {
	my ($dbh,$refModifSpecies)=@_;
	my $sthIsRef=$dbh->prepare("SELECT 1 FROM SPECIES WHERE ID_SPECIES=? AND IS_REFERENCE=1 LIMIT 0,1");
	my $sthNumSp1=$dbh->prepare("SELECT 1 FROM GOANNOTATION WHERE ID_SPECIES=? LIMIT 0,1");
	my $sthNumSp2=$dbh->prepare("SELECT 1 FROM IDENTIFIER WHERE ID_SPECIES=? LIMIT 0,1");
	my $sthNumSp3=$dbh->prepare("SELECT 1 FROM MASTER_PROTEIN WHERE ID_SPECIES=? LIMIT 0,1");
	my $sthDelSp=$dbh->prepare("DELETE FROM SPECIES WHERE ID_SPECIES=?");
	foreach my $speciesID (keys %{$refModifSpecies}) {
		$sthIsRef->execute($speciesID);
		my ($isRef)=$sthIsRef->fetchrow_array; # is reference
		unless ($isRef) {
			$sthNumSp1->execute($speciesID);
			my ($spUsed)=$sthNumSp1->fetchrow_array; # GO
			unless ($spUsed) {
				$sthNumSp2->execute($speciesID);
				($spUsed)=$sthNumSp2->fetchrow_array; # IDENTIFIER
				unless ($spUsed) {
					$sthNumSp3->execute($speciesID);
					($spUsed)=$sthNumSp3->fetchrow_array; # MASTER_PROTEIN
					$sthDelSp->execute($speciesID) unless $spUsed;
				}
			}
		}
	}
	$sthNumSp1->finish;
	$sthNumSp2->finish;
	$sthNumSp3->finish;
	$sthDelSp->finish;
}

######################################## creates directory if not exists
####<Clean directory or delete file>#### delete everything inside older than maxAge (directory itself is never deleted!)
######################################## returns the number of elements remaining
sub cleanDirectory {
	my ($dir,$maxAge)=@_;
		
	if (!-e $dir) {
		mkdir $dir;
		return 1;
	}
	return 1 if -l $dir; # skip links
	
	if ($maxAge) {
		if ($maxAge=~/^(\d+)s/i) {$maxAge=$1;} # seconds
		elsif ($maxAge=~/^(\d+)m/) {$maxAge=$1*60;} # minutes
		elsif ($maxAge=~/^(\d+)h/i) {$maxAge=$1*3600;} # hour
		elsif ($maxAge=~/^(\d+)d/i) {$maxAge=$1*3600*24;} # day
		elsif ($maxAge=~/^(\d+)M/) {$maxAge=$1*3600*24*30.5;} # Month
		elsif ($maxAge=~/^(\d+)y/i) {$maxAge=$1*3600*24*30.5*12;} # year
		else {$maxAge=~s/^(\d+).+/$1/;} # assumes seconds
	}
	$maxAge=3600*24*30 unless $maxAge; # Defaults to 1 month

	my $numChildren=0;
	my $now=time;
	if (-d $dir) { # is a directory
		return 1 if -e "$dir/no_delete.txt"; # safety net for dev data

		##>Clean previous job dirs & files
		opendir(my $DIR,$dir);
		while (defined (my $child = readdir($DIR))) {
			next if ($child eq '.' || $child eq '..');
			$numChildren++;
			next if -l "$dir/$child"; # skip links
			my $modTime=(stat("$dir/$child"))[9];
			if ($modTime && $now-$modTime > $maxAge) {
				if (-d "$dir/$child") { # old dir
					rmtree "$dir/$child";
				}
				else {unlink "$dir/$child";} # old file or link
				$numChildren--;
			}
		}
		close $DIR;
	}
	elsif (-e $dir) { # is a file
		my $modTime=(stat($dir))[9];
		if ($modTime && $now-$modTime > $maxAge) {unlink $dir;}
		else {$numChildren=1;}
	}

	return $numChildren;
}

######################################
####<Directory browsing functions>####
######################################
sub browseDirectory_JavaScript { # Not required if form element is 'SELECT' in &browseDirectory_getFiles
	print qq
|/**** Start of browseDirectory JavaScript functions ****/
function browseDirectoryExpandSelection(dirBox,dirIdx,chkBoxName) {
	var boxIdData=dirBox.id.split('_');
	var re=(boxIdData[1]=='0')? new RegExp(chkBoxName+':'+'browseDirectoryDIV_') : new RegExp(chkBoxName+':'+'browseDirectoryDIV_'+boxIdData[1]);
	var divList=document.getElementsByName(chkBoxName+':'+'browseDirectoryDIV');
	for (let i=(dirIdx+1); i<divList.length; i++) {
		if (divList[i].id.match(re)) {
			var divIdData=divList[i].id.split('_');
			document.getElementById(chkBoxName+':'+'browseDirectoryBOX_'+divIdData[1]).checked=dirBox.checked;
		}
		else {break;}
	}
	browseDirectoryUpdateNumSelectedFiles(chkBoxName);
}
function browseDirectoryUpdateNumSelectedFiles(chkBoxName) {
	var numFilesChecked=0;
	var filesChkBox=document.getElementsByName(chkBoxName);
	for (let i=0; i<filesChkBox.length; i++) {
		if (filesChkBox[i].checked) {numFilesChecked++;}
	}
	var fileStrg=(numFilesChecked <= 1)? 'file' : 'files';
	document.getElementById(chkBoxName+':'+'browseDirectorySPAN').innerHTML=numFilesChecked+' '+fileStrg+' selected';
}
function browseDirectoryExpandCollapsDirectory(dirId,dirIdx,myButton,chkBoxName) {
	var newVis,newValue;
	if (myButton.value=='-') { // -> hide
		newVis='none';
		newValue='+';
	}
	else { // -> show
		newVis='';
		newValue='-';
	}
	myButton.value=newValue;
	var re=new RegExp(chkBoxName+':'+'browseDirectoryDIV_'+dirId);
	var divList=document.getElementsByName(chkBoxName+':'+'browseDirectoryDIV');
	for (let i=(dirIdx+1); i<divList.length; i++) {
		if (divList[i].id.match(re)) {
			var divIdData=divList[i].id.split('_');
			var childButton=document.getElementById(chkBoxName+':'+'browseDirectoryBUTTON_'+divIdData[1]);
			if (childButton) {childButton.value=newValue;} // to make sure child button status matches parent
			divList[i].style.display=newVis;
		}
		else {break;}
	}
}
/**** End of browseDirectory JavaScript functions ****/
|;
}

sub browseDirectory_getFiles {
	my ($rootDir,$inputName,$refParam)=@_;
	##### $refParam keys: #####
	# -returnString: return result in a string. Do not print directly (default: print)
	# -maxDepth: max depth of sub directories (default: 5)
	# -fileMatch: full regex to be matched my files eg. "qr/\.xml/i" (default: undef, all files are matched)
	# -hideEmpty: hide directory with no (matching) files for checkbox only! (default: 1/true)
	# -inputType: checkbox or select (for single file selection) (optional) (default: checkbox)
	# -selectAttribStrg: string for "CSS class/id/style/JS event/disabled/..." for SELECT element only! (default: undef)
	$refParam={} unless $refParam;
	$refParam->{maxDepth}=5 unless $refParam->{maxDepth};
	$refParam->{hideEmpty}=1 unless defined($refParam->{hideEmpty});
	$refParam->{inputType}='checkbox' unless defined($refParam->{inputType}); # checkbox or select
	$refParam->{inputType}=lc($refParam->{inputType});
	$refParam->{selectAttribStrg}='' unless defined($refParam->{selectAttribStrg});

	###<Nested local recursive subroutine (must be defined before called)>###
	local *readDir = sub {
		my ($rootDir,$parentDir,$parentID,$refDirTree,$dirDepth,$refParam)=@_;
		my $currentDir=($parentDir)? "$rootDir/$parentDir" : $rootDir;
		opendir (my $DIR,$currentDir) || warn("Unable to read '$currentDir'!");
		my $matchPattern=$refParam->{'fileMatch'}; # can be undef
		my $rank=0;
		#while (defined (my $dirFile = readdir ($DIR))) {#}
		my @files = sort {&sortSmart($a,$b)} readdir($DIR);
		while (my $dirFile = shift @files) {
			next if $dirFile =~ /^\./; # (parent)dir or hidden file;
			my $fullDirFile="$currentDir/$dirFile";
			next if -l $fullDirFile; # do not follow symbolic links
			my $type=(-d $fullDirFile)? 'dir' : 'file';
			next if ($type eq 'file' && $matchPattern && $dirFile !~ $matchPattern);
			$rank++;
			my $id=($parentID)? "$parentID.$rank" : $rank;
			push @{$refDirTree},[$type,$id,$dirDepth,$parentDir,$dirFile,0]; # '0' disabled flag
			if ($type eq 'dir') { # directory
				my $refLeafIdx=$#{$refDirTree};
				my $newParentDir=($parentDir)? "$parentDir/$dirFile" : $dirFile;
				my ($numChildren)=($dirDepth==$refParam->{'maxDepth'})? 0 : &readDir($rootDir,$newParentDir,$id,$refDirTree,$dirDepth+1,$refParam);
				unless ($numChildren) {
					if ($refParam->{'hideEmpty'}) { # skip empty directory
						pop @{$refDirTree};
						$rank--;
					}
					else {$refDirTree->[$refLeafIdx][-1]=1;} # set disabled flag to 1
				}
			}
		}
		close $DIR;
		return $rank;
	};

	###<Read directories recursively>###
	my @dirTree=();
	&readDir($rootDir,'','',\@dirTree,0,$refParam);

	###<Display directories and files>###
	my $resultString='';
	if ($refParam->{inputType} eq 'select') {
		$resultString.="<SELECT name=\"$inputName\" $refParam->{selectAttribStrg}>\n";
		if (scalar @dirTree) {
			$resultString.="<OPTION value=\"\">-= Select =-</OPTION>\n";
			foreach my $f (0..$#dirTree) {
				my ($type,$id,$dirDepth,$parentDir,$file,$disabled)=@{$dirTree[$f]};
				next if $type eq 'dir';
				my $fileValue=($parentDir)? "$parentDir/$file" : $file;
				$resultString.="<OPTION value=\"$fileValue\">$fileValue</OPTION>\n";
			}
		}
		else {
			$resultString.="<OPTION value=\"\">***No matching files found.***</OPTION>\n";
		}
		$resultString.="<SELECT>\n";
	}
	else {
		$resultString.=qq
|<DIV class="darkBg" style="font-weight:bold">
<INPUT type="checkbox" id="$inputName:browseDirectoryBOX_0" value="/" onclick="browseDirectoryExpandSelection(this,-1,'$inputName')"/>/&nbsp;[<SPAN id=\"$inputName:browseDirectorySPAN\">0 file selected</SPAN>]
</DIV>
|;
		my @prevParents=(''); # index=depth
		my $prevDepth=0;
		foreach my $f (0..$#dirTree) {
			my ($type,$id,$dirDepth,$parentDir,$file,$disabled)=@{$dirTree[$f]};
			if ($dirDepth < $prevDepth) {
				splice(@prevParents,($dirDepth+1)); # splice works on length, not index
			}
			my ($directParentDir,$fileValue,$fileStrg)=('','','');
			if ($parentDir) {
				my $depth=0;
				foreach my $dir (split/\//,$parentDir) {
					$depth++;
					if ($prevParents[$depth] && $prevParents[$depth] eq $dir) {$fileStrg.='<FONT style="visibility:hidden">'.$dir.'/</FONT>';}
					else {$fileStrg.="$dir/";}
					$directParentDir=$dir;
				}
				$fileValue="$parentDir/$file";
			}
			else {$fileValue=$file;}
			$prevDepth=$dirDepth;
			my $divClass=($type eq 'dir')? 'darkBg' : 'lightBg';
			$resultString.="<DIV id=\"$inputName:browseDirectoryDIV_$id\" name=\"$inputName:browseDirectoryDIV\" class=\"$divClass\"><FONT style=\"visibility:hidden\">&nbsp;&nbsp;</FONT>";
			if ($type eq 'dir') {
				my $disabStrg=($disabled)? ' disabled' : '';
				$resultString.="$fileStrg<INPUT type=\"checkbox\" name=\"dir_$inputName\" id=\"$inputName:browseDirectoryBOX_$id\" value=\"$fileValue\" onclick=\"browseDirectoryExpandSelection(this,$f,'$inputName')\"$disabStrg/>$file/";
				unless ($disabled) {
					$resultString.="&nbsp;<INPUT type=\"button\" id=\"$inputName:browseDirectoryBUTTON_$id\" value=\"-\" style=\"height:15px;width:15px;padding:0;line-height:5px;\" onclick=\"browseDirectoryExpandCollapsDirectory('$id',$f,this,'$inputName')\"/>";
				}
				$prevParents[$dirDepth+1]=$file;
			}
			else {
				$resultString.="$fileStrg<INPUT type=\"checkbox\" name=\"$inputName\" id=\"$inputName:browseDirectoryBOX_$id\" value=\"$fileValue\" onclick=\"browseDirectoryUpdateNumSelectedFiles('$inputName')\"/>$file\n";
				$prevParents[$dirDepth]=$directParentDir;
			}
			$resultString.="</DIV>\n";
		}
		$resultString.="<DIV class=\"lightBg\" >No matching files found.</DIV>\n" unless scalar @dirTree;
	}

	###<Return or print>###
	if ($refParam->{returnString}) {return $resultString;}
	else {print $resultString;}
}

############################################################
####<Return search parameter fron analysis results File>####
############################################################
sub getSearchParam {
	my ($dbh,$anaID)=@_;
	my %promsPath=&promsConfig::getServerInfo;
	my ($fileFormat,$validStatus,$fileName,$rawName)=$dbh->selectrow_array("SELECT FILE_FORMAT,VALID_STATUS,DATA_FILE,NAME FROM ANALYSIS WHERE ID_ANALYSIS=$anaID");
	my %infoSubmit;
	my $projectID=&promsMod::getProjectID($dbh,$anaID,'analysis');
	$fileName=~s/\.xml/\.pgf/ if $fileFormat=~/(PHENYX|MASCOT)\.XML/;
	my $anaDir=($validStatus==2)? "$promsPath{peptide}/proj_$projectID/ana_$anaID" : "$promsPath{valid}/ana_$anaID";

	##<MaxQuant
	if ($fileFormat eq 'MAXQUANT.DIR') {
		$infoSubmit{'a:Search title'} = $rawName;
		if ($fileName eq 'parameters.txt') { # no mqpar.xml provided
			my %allParameters;
			open(PARAM,"$anaDir/parameters.txt");
			while(<PARAM>) {
				next if $.==1;
				s/\s+$//; # chomp is not enough <- Windows
				my ($param,$valueStrg)=split(/\t/,$_);
				$valueStrg='' unless $valueStrg;
				$allParameters{$param}=$valueStrg;
			}
			close PARAM;
			$infoSubmit{'c:MaxQuant version'}=$allParameters{'Version'} || 'Unknown';
			$infoSubmit{'f:Fixed modifications'}=join(', ',split(';',$allParameters{'Fixed modifications'})) || 'None';
			#$infoSubmit{'g:Variable modifications'}=join(', ',split(';',$allParameters{'Site tables'})) || 'None'; $infoSubmit{'g:Variable modifications'}=~s/Sites\.txt//g;
			my $lKey=97; # chr(97)='a'
			foreach my $param (keys %allParameters) {
				if ($param=~/^MS\/MS tol/) {
					push @{$infoSubmit{'h:MS/MS tolerance'}{chr($lKey).':'.$param}},$allParameters{$param};
					$lKey++;
				}
			}
			foreach my $file (split(';',$allParameters{'Fasta file'})) {
				my $fastaFile=(split(/[\\\/]/,$file))[-1];
				push @{$infoSubmit{'l:Databank'}{'a:File'}},$fastaFile;
			}
			$infoSubmit{'m:Include contaminants'}=$allParameters{'Include contaminants'};
			$infoSubmit{'n:Decoy mode'}=$allParameters{'Decoy mode'};
			
			my %summaryColNum;
			open(SUMM,"$anaDir/summary.txt");
			while(<SUMM>) {
				s/\s+$//; # chomp is not enough <- Windows
				my @parameters=split(/ *\t */,$_);
				if ($.==1) {
					my %msmsColNum;
					my $ncol=0;
					foreach my $colName (@parameters) {
						$summaryColNum{$colName}=$ncol;
						$ncol++;
					}
					next;
				}
				foreach my $i (0..$#parameters) {$parameters[$i]='' unless $parameters[$i];} # just to be safe
				# take only line 2 => assume same parameters for all raw files
				$infoSubmit{'d:Enzyme'}=$parameters[$summaryColNum{'Enzyme'}];
				$infoSubmit{'e:Max. missed cleavages'}=$parameters[$summaryColNum{'Max. missed cleavages'}] || 'Unknown';
				$infoSubmit{'g:Variable modifications'}=join(', ',split(';',$parameters[$summaryColNum{'Variable modifications'}])) || 'None';
				last;
			}
			close SUMM;
			$infoSubmit{'o:Parameter files'}='parameters.txt, summary.txt';
		}
		elsif ($fileName=~/\.xml$/) { # mqpar.xml
			require XML::Simple;
			my $xml = new XML::Simple();
			my $xmlParams = $xml->XMLin("$anaDir/$fileName",ForceArray=>['parameterGroup','string','short','int'],SuppressEmpty=>undef);	
			my $paramGrIdx=0; # Assumes only 1 parameter group!!!!!!!!!!!!!!!!!!!!!
			my $xmlGroupParams=$xmlParams->{parameterGroups}{parameterGroup}[$paramGrIdx];
			$infoSubmit{'c:MaxQuant version'} = $xmlParams->{maxQuantVersion} if $xmlParams->{maxQuantVersion};
			$infoSubmit{'d:Enzyme'}=join(',',@{$xmlGroupParams->{enzymes}{string}});
			$infoSubmit{'e:Max. missed cleavages'}=$xmlGroupParams->{maxMissedCleavages};
			$infoSubmit{'f:Fixed modifications'}=join(', ',@{$xmlParams->{fixedModifications}{string}}) if $xmlParams->{fixedModifications}{string}[0];
			$infoSubmit{'g:Variable modifications'}=join(', ',@{$xmlGroupParams->{variableModifications}{string}}) if $xmlGroupParams->{variableModifications}{string}[0];
			my $lKey=97; # chr(97)='a'
			foreach my $refMsmsParam (@{$xmlParams->{msmsParamsArray}{msmsParams}}) {
				my $name=$refMsmsParam->{Name};
				my $tol=$refMsmsParam->{MatchTolerance};
				my $unit=($refMsmsParam->{MatchToleranceInPpm} eq 'true')? 'ppm' : 'Da';
				push @{$infoSubmit{'h:MS/MS tolerance'}{chr($lKey).':'.$name}},"$tol $unit";
				$lKey++;
			}
			foreach my $file (@{$xmlParams->{fastaFiles}{string}}) {
				my $fastaFile=(split(/[\\\/]/,$file))[-1];
				push @{$infoSubmit{'l:Databank'}{'a:File'}},$fastaFile;
			}
			$infoSubmit{'m:Include contaminants'}=$xmlParams->{includeContaminants};
			$infoSubmit{'n:Decoy mode'}=$xmlParams->{decoyMode};
			$infoSubmit{'o:XML parameter file'}=$fileName;
		}
		return %infoSubmit;
	}
	##<All other search engines
	my $searchFile="$anaDir/$fileName";
	if ($validStatus==2 && $fileFormat=~/\.DAT\Z/) { # checking for minimal dat file
		(my $minDatFile=$fileName)=~s/\.dat/_min\.dat/;
		$searchFile=(-e "$anaDir/$minDatFile")? "$anaDir/$minDatFile" : "$anaDir/$fileName";
	}
	unless (-e $searchFile || $fileFormat =~ /SWATH|SKYLINE|PDM/) {
		$infoSubmit{'Error'}="Search file not found";
		return %infoSubmit;
	}

	##<Paragon
	if ($fileFormat eq 'PARAGON.XML') { # For PARAGON group files (ProteinPilot) converted into xml
		#delete ($infoSubmit{'Error'});
		#if ($validStatus<2) {$searchFile="$promsPath{'valid'}/ana_$anaID/$fileName";}
		#else {
		#	if(-e "$promsPath{peptide}/proj_$projectID/${anaID}_$fileName"){
		#		$searchFile="$promsPath{peptide}/proj_$projectID/${anaID}_$fileName";
		#	}
		#	my $pepFileName=sprintf "P%06d.xml",$anaID;
		#	$searchFile="$promsPath{peptide}/proj_$projectID/$pepFileName";
		#}
		open (FILE, $searchFile) || ($infoSubmit{'Error'}="Unable to open Search file");
		my $section="";
		while (my $line=<FILE>) {
			if($line =~ /<SEARCH_STAT num_proteins="(\d+)" num_residues="(\d+)"/) {
				push @{$infoSubmit{'l:Databank'}{'e:Sequences after taxonomy filter'}},$1;
				push @{$infoSubmit{'l:Databank'}{'f:Residues after taxonomy filter'}},$2;
				$section='SEARCH';
			}
			elsif ($line =~ /<PARAGON_METHOD_SETTINGS>/){
				$section="PARAGON_SETTINGS";
			}
			elsif ($line =~ /<\/PARAGON_METHOD_SETTINGS>/){
				$section="";
				last;
			}
			elsif ($line =~ /<PARAGON_VERSION>(.+)<\/PARAGON_VERSION>/){
				$infoSubmit{'c:Paragon version'}=$1;
			}
			#elsif ($line =~ /<QUANT_TYPE .* xml:id="(.+)" >/) {# only if a quantification was set
			#	$infoSubmit{'q:Quantification'}=$1;
			#}
			if ($section eq "SEARCH"){
				if ($line =~ /<MSTOLERANCE type="(\w+)" value="(.+)"/){
					my $value=sprintf("%.1f", $2);
					$infoSubmit{'h:Peptide tolerance'}=$value;
					$infoSubmit{'h:Peptide tolerance'}.=" $1";
				}
				elsif ($line =~ /<MSMSTOLERANCE type="(\w+)" value="(.+)"/){
					my $value=sprintf("%.1f", $2);
					$infoSubmit{'i:Fragment tolerance'}=$value;
					$infoSubmit{'i:Fragment tolerance'}.=" $1";
				}
				elsif ($line =~ /<FASTA filename="(.+)"/) {
					push @{$infoSubmit{'l:Databank'}{'b:File'}},$1;
					$section='';
				}
			}
			elsif($section eq "PARAGON_SETTINGS"){
				if ($line =~ /<UI_DIGESTION>(\w+)<\/UI_DIGESTION>/){
					$infoSubmit{'d:Enzyme'}=$1;
				}
				elsif ($line =~ /<UI_START_TIME>(.+)<\/UI_START_TIME>/) {
					$infoSubmit{'r:Search date'}=$1;
				}
				elsif ($line =~ /<UI_INSTRUMENT>(.+)<\/UI_INSTRUMENT>/) {
					$infoSubmit{'b:Instrument'}=$1;
				}
				elsif ($line =~ /<UI_SPECIES>(.+)<\/UI_SPECIES>/) {
					$infoSubmit{'m:Taxonomy'}=$1;
				}
				elsif ($line =~ /<UI_QUANT_TYPE>(.+)<\/UI_QUANT_TYPE>/) { # only if a quantification was set
					$infoSubmit{'q:Quantification'}=$1;
				}
			}
		}
		close FILE;

		$searchFile=~s/\.xml/\.pgf/;
		open (VMODFILE, $searchFile) || ($infoSubmit{'Error'}="Unable to open Search file");
		$infoSubmit{'g:Variable modifications'}="";
		my $onMasses=0;
		while (my $line=<VMODFILE>) {
			chomp($line);
			if ($line eq '--gc0p4Jq0M2Yt08jU534c0p' ){
				$onMasses=0;
				next;
			}
			elsif ($line eq 'Content-Type: Paragon; name="masses"'){
				$onMasses=1;
				next;
			}
			if ($onMasses) {
				my @vmod=split(/,/,$line);
				$infoSubmit{'g:Variable modifications'}.="$vmod[1],";
			}
		}
		close (VMODFILE);
		chop($infoSubmit{'g:Variable modifications'}) if $infoSubmit{'g:Variable modifications'};
	}
	elsif ($fileFormat eq 'SWATH.PKV') {			## SWATH quantification with PeakView
        my $searchParam=$dbh->selectrow_array("SELECT QUANTIF_ANNOT FROM QUANTIFICATION Q, ANA_QUANTIFICATION AQ WHERE Q.ID_QUANTIFICATION=AQ.ID_QUANTIFICATION AND ID_ANALYSIS=$anaID");
		foreach my $param (split(/::/,$searchParam)){
			my ($name,$value)=split(/=/,$param);
			next if $value=~/N\/A/;
			$infoSubmit{'g:Peptide tolerance'}=$value if $name=~/XIC_WIDTH/;
			$infoSubmit{'h:Extraction window'}=$value if $name=~/XIC_EXTRACTION_WINDOW/;
			if ($name=~/EXCLUDE_MODIFIED_PEPTIDES/) {
                $infoSubmit{'f:Exclude modified peptides'}=($value==1)? 'Selected': 'No selected';
            }
			$infoSubmit{'e:Transition per peptide'}=$value if $name=~/NB_TRANSITION_PER_PEPTIDE/;
			$infoSubmit{'d:Peptides per protein'}=$value if $name=~/NB_PEPTIDES_PER_PROTEIN/;
			push @{$infoSubmit{'b:Spectral library'}{'a:Name'}},$value if $name=~/LIBRARY_NAME/;
			push @{$infoSubmit{'b:Spectral library'}{'b:Version'}},$value if $name=~/LIBRARY_VERSION/;
			push @{$infoSubmit{'b:Spectral library'}{'c:Peptides'}},$value if $name=~/LIBRARY_PEPTIDES/;
			push @{$infoSubmit{'b:Spectral library'}{'d:Proteins'}},$value if $name=~/LIBRARY_PROTEINS/;
			push @{$infoSubmit{'b:Spectral library'}{'e:Specific proteins'}},$value if $name=~/LIBRARY_SPECIFICS_PROTEINS/;

			## Library export options
			push @{$infoSubmit{'c:Library export options'}{'a:Ion mass limits'}},$value if $name=~/ION_MASS_LIMITS/;
			push @{$infoSubmit{'c:Library export options'}{'b:Ion type'}},$value if $name=~/ION_TYPE/;
			push @{$infoSubmit{'c:Library export options'}{'c:Ion charge'}},$value if $name=~/ION_CHARGE/;
			push @{$infoSubmit{'c:Library export options'}{'d:Ion number per peptide'}},$value if $name=~/ION_NUMBER_PER_PEPTIDE/;
			push @{$infoSubmit{'c:Library export options'}{'e:Maximum permissible error'}},$value if $name=~/MAXIMUM_ERROR_ALLOWED/;
			push @{$infoSubmit{'c:Library export options'}{'f:Time scale'}},$value if $name=~/TIME_SCALE/;
			push @{$infoSubmit{'c:Library export options'}{'g:USI order'}},$value if $name=~/UIS_ORDER/;
			push @{$infoSubmit{'c:Library export options'}{'h:Swath windows file'}},$value if $name=~/SWATH_FILE/;
			push @{$infoSubmit{'c:Library export options'}{'i:Modification file'}},$value if $name=~/MODIFICATION_FILE/;
			push @{$infoSubmit{'c:Library export options'}{'j:Fragment mass modifications allowed'}},$value if $name=~/FRAGMENT_MASS_MODIFICATIONS_ALLOWED/;
			push @{$infoSubmit{'c:Library export options'}{'k:Labelling file'}},$value if $name=~/LABELLING_FILE/;
			push @{$infoSubmit{'c:Library export options'}{'l:Fasta file'}},$value if $name=~/FASTA_FILE/;
			push @{$infoSubmit{'c:Library export options'}{'m:Other'}},$value if $name=~/OTHER/;
		}
		
		my $sthDbID=$dbh->prepare("SELECT ID_DATABANK FROM ANALYSIS_DATABANK WHERE ID_ANALYSIS=$anaID");
		$sthDbID->execute;
		while (my $dbID=$sthDbID->fetchrow_array) {
            my ($dbName,$dbFile,$numEntry)=$dbh->selectrow_array("SELECT NAME,FASTA_FILE,NUM_ENTRY FROM DATABANK WHERE ID_DATABANK=$dbID");
			push @{$infoSubmit{'a:Databank'}{'a:Name'}},$dbName if $dbName;
			push @{$infoSubmit{'a:Databank'}{'b:File'}},$dbFile;
			push @{$infoSubmit{'a:Databank'}{'c:Sequences in databank'}},$numEntry;
        }
    }
	elsif($fileFormat eq 'OPENSWATH.TSV'){
		my $searchParam=$dbh->selectrow_array("SELECT QUANTIF_ANNOT FROM QUANTIFICATION Q, ANA_QUANTIFICATION AQ WHERE Q.ID_QUANTIFICATION=AQ.ID_QUANTIFICATION AND ID_ANALYSIS=$anaID");
		foreach my $param (split(/::/,$searchParam)){
			my ($name,$value)=split(/=/,$param);
			next if $value=~/N\/A/;
			push @{$infoSubmit{'b:Spectral library'}{'a:Name'}},$value if $name=~/LIBRARY_NAME/;
			push @{$infoSubmit{'b:Spectral library'}{'b:Version'}},$value if $name=~/LIBRARY_VERSION/;
			push @{$infoSubmit{'b:Spectral library'}{'c:Peptides'}},$value if $name=~/LIBRARY_PEPTIDES/;
			push @{$infoSubmit{'b:Spectral library'}{'d:Proteins'}},$value if $name=~/LIBRARY_PROTEINS/;
			push @{$infoSubmit{'b:Spectral library'}{'e:Specific proteins'}},$value if $name=~/LIBRARY_SPECIFICS_PROTEINS/;

			## Export options
			push @{$infoSubmit{'c:Library export options'}{'a:Ion mass limits'}},$value if $name=~/ION_MASS_LIMITS/;
			push @{$infoSubmit{'c:Library export options'}{'b:Ion type'}},$value if $name=~/ION_TYPE/;
			push @{$infoSubmit{'c:Library export options'}{'c:Ion charge'}},$value if $name=~/ION_CHARGE/;
			push @{$infoSubmit{'c:Library export options'}{'d:Ion number per peptide'}},$value if $name=~/ION_NUMBER_PER_PEPTIDE/;
			push @{$infoSubmit{'c:Library export options'}{'e:Maximum permissible error'}},$value if $name=~/MAXIMUM_ERROR_ALLOWED/;
			push @{$infoSubmit{'c:Library export options'}{'f:Time scale'}},$value if $name=~/TIME_SCALE/;
			push @{$infoSubmit{'c:Library export options'}{'g:USI order'}},$value if $name=~/UIS_ORDER/;
			push @{$infoSubmit{'c:Library export options'}{'h:Swath windows file'}},$value if $name=~/SWATH_FILE/;
			push @{$infoSubmit{'c:Library export options'}{'i:Modification file'}},$value if $name=~/MODIFICATION_FILE/;
			push @{$infoSubmit{'c:Library export options'}{'j:Fragment mass modifications allowed'}},$value if $name=~/FRAGMENT_MASS_MODIFICATIONS_ALLOWED/;
			push @{$infoSubmit{'c:Library export options'}{'k:Labelling file'}},$value if $name=~/LABELLING_FILE/;
			push @{$infoSubmit{'c:Library export options'}{'l:Fasta file'}},$value if $name=~/FASTA_FILE/;
			push @{$infoSubmit{'c:Library export options'}{'m:Other'}},$value if $name=~/OTHER/;
			
			## Library processing options
			push @{$infoSubmit{'d:Library processing options'}{'a:Precursor mz threshold'}},$value if $name=~/PRECURSOR_MZ_THRESHOLD/;
			push @{$infoSubmit{'d:Library processing options'}{'b:Product mz threshold'}},$value if $name=~/PRODUCT_MZ_THRESHOLD/;
			push @{$infoSubmit{'d:Library processing options'}{'c:Precursor lower mz limit'}},$value if $name=~/PRECURSOR_MIN_MZ/;
			push @{$infoSubmit{'d:Library processing options'}{'d:Precursor upper mz limit'}},$value if $name=~/PRECURSOR_MAX_MZ/;
			push @{$infoSubmit{'d:Library processing options'}{'f:Unimod modifications file'}},$value if $name=~/UNIMOD_MOD_FILE/;
			push @{$infoSubmit{'d:Library processing options'}{'f:Max nb alt localizations'}},$value if $name=~/MAX_NB_ALT_LOCALIZATIONS/;
			push @{$infoSubmit{'d:Library processing options'}{'e:IPF'}}, $value if $name=~/IPF/;

			## OpenSwath options
			push @{$infoSubmit{'e:OpenSwath Options'}{'a:IRT file'}}, $value if $name=~/IRT_FILE/;
			push @{$infoSubmit{'e:OpenSwath Options'}{'b:Swath windows file'}}, $value if $name=~/SWATH_WINDOWS/;
			push @{$infoSubmit{'e:OpenSwath Options'}{'c:FDR cutoff'}}, $value if $name=~/FDR_CUTOFF/;
			push @{$infoSubmit{'e:OpenSwath Options'}{'d:RT extraction windows'}}, $value if $name=~/RT_EXTRACTION_WINDOWS/;
			push @{$infoSubmit{'e:OpenSwath Options'}{'e:m/z extraction windows'}}, $value if $name=~/MZ_EXTRACTION_WINDOWS/;
			push @{$infoSubmit{'e:OpenSwath Options'}{'f:m/z extraction unit'}}, "ppm" if($name=~/MZ_EXTRACTION_WINDOWS_UNIT/);
			push @{$infoSubmit{'e:OpenSwath Options'}{'g:m/z correction function'}}, $value if $name=~/MZ_CORRECTION_FUNCTION/;
			push @{$infoSubmit{'e:OpenSwath Options'}{'h:Minimum peak width'}}, $value if $name=~/MIN_PEAK_WIDTH/;
			
			## PyProphet options
			push @{$infoSubmit{'f:PyProphet Options'}{'h:Method'}}, $value if $name=~/PYPROPHET_METHOD/;
			
			## TRIC options
			push @{$infoSubmit{'g:TRIC option'}{'m:Method'}},$value if $name=~/TRIC_METHOD/;
		}
		
		if(scalar @{$infoSubmit{'e:OpenSwath Options'}{'f:m/z extraction unit'}} == 0) {
			push @{$infoSubmit{'e:OpenSwath Options'}{'f:m/z extraction unit'}}, "Th";
		}

		my $sthDbID=$dbh->prepare("SELECT ID_DATABANK FROM ANALYSIS_DATABANK WHERE ID_ANALYSIS=$anaID");
		$sthDbID->execute;
		while (my $dbID=$sthDbID->fetchrow_array) {
            my ($dbName,$dbFile,$numEntry)=$dbh->selectrow_array("SELECT NAME,FASTA_FILE,NUM_ENTRY FROM DATABANK WHERE ID_DATABANK=$dbID");
			push @{$infoSubmit{'a:Databank'}{'a:Name'}},$dbName if $dbName;
			push @{$infoSubmit{'a:Databank'}{'b:File'}},$dbFile;
			push @{$infoSubmit{'a:Databank'}{'c:Sequences in databank'}},$numEntry;
        }
	}
	elsif($fileFormat eq 'SKYLINE.CSV') {
		## PARSE PARAMETERS
		my $searchParam=$dbh->selectrow_array("SELECT QUANTIF_ANNOT FROM QUANTIFICATION Q, ANA_QUANTIFICATION AQ WHERE Q.ID_QUANTIFICATION=AQ.ID_QUANTIFICATION AND ID_ANALYSIS=$anaID");

		foreach my $param (split(/::/,$searchParam)){
			my ($name,$value)=split(/=/,$param);
			next if $value=~/N\/A/;
			push @{$infoSubmit{'b:Library'}{'a:Cut-off score'}},$value if $name=~/LIB_CUTOFF/ && $value;
			push @{$infoSubmit{'b:Library'}{'b:Search files'}},$value if $name=~/LIB_SEARCH_FILES/ && $value;
			push @{$infoSubmit{'b:Library'}{'c:Include ambigous matches'}},$value if $name=~/LIB_AMBIGOUS/ && $value;
			
			if($name =~ /SOFTWARE/) {
				my ($soft, $softVersion) = split(/;/, $value);
				$infoSubmit{'c:Software'} = "Skyline ($softVersion)";
			}

			push @{$infoSubmit{'d:MS1 Filtering'}{'a:Precursor charges'}},$value if $name=~/MS1_PRECURSOR_CHARGE/ && $value;
			push @{$infoSubmit{'d:MS1 Filtering'}{'b:Isotope peaks included'}},$value if $name=~/MS1_ISOTOPE_PEAKS/ && $value;
			push @{$infoSubmit{'d:MS1 Filtering'}{'c:Precursor mass analyzer'}},$value if $name=~/MS1_MASS_ANALYZER/ && $value;
			push @{$infoSubmit{'d:MS1 Filtering'}{'d:Peaks'}},$value if $name=~/MS1_PEAKS/ && $value;
			push @{$infoSubmit{'d:MS1 Filtering'}{'e:Resolving power'}},$value if $name=~/MS1_POWER/ && $value;
			push @{$infoSubmit{'d:MS1 Filtering'}{'f:Resolving power charge'}},$value if $name=~/MS1_POWER_CHARGE/ && $value;
			
			push @{$infoSubmit{'e:Retention time'}{'a:Scan'}},$value if $name=~/RT_SCANS/ && $value;
			push @{$infoSubmit{'e:Retention time'}{'b:Filter'}},$value if $name=~/RT_SCANS_FILTER/ && $value;
			
			push @{$infoSubmit{'f:Fasta'}{'b:Comment'}},$value if $name=~/FASTA_COMMENT/ && $value;
			push @{$infoSubmit{'f:Fasta'}{'c:Enzyme'}},$value if $name=~/FASTA_ENZYME/ && $value;
			push @{$infoSubmit{'f:Fasta'}{'d:Max missed cleavage'}},$value if $name=~/FASTA_MAX_MC/ && $value;
			push @{$infoSubmit{'f:Fasta'}{'e:Filter proteins'}},"Having less than $value peptides" if $name=~/FASTA_FILTER/ && $value;
			
			if ($name =~ /MODIFS/ && $value) {
				$value =~ s/,/, /g;
				$infoSubmit{'g:Modifications'} = $value;
			}
		}
		my $sthDbID=$dbh->prepare("SELECT ID_DATABANK FROM ANALYSIS_DATABANK WHERE ID_ANALYSIS=$anaID");
		$sthDbID->execute;
		while (my $dbID=$sthDbID->fetchrow_array) {
            my ($dbName,$dbFile,$numEntry)=$dbh->selectrow_array("SELECT NAME,FASTA_FILE,NUM_ENTRY FROM DATABANK WHERE ID_DATABANK=$dbID");
			push @{$infoSubmit{'a:Databank'}{'a:Name'}},$dbName if $dbName;
			push @{$infoSubmit{'a:Databank'}{'b:File'}},$dbFile;
			push @{$infoSubmit{'a:Databank'}{'c:Sequences in databank'}},$numEntry;
        }
	}
	elsif ($fileFormat eq 'TDM.PEP.XML') {			## X! TANDEM
		open (FILE, $searchFile) || ($infoSubmit{'Error'}="Unable to open Search file");
		while (my $line=<FILE>) {
			if ($line=~/scoring, maximum missed cleavage sites/) {
				($infoSubmit{'e:Enzyme misscleavage'})=($line=~/value="(\d)"/);
			}
			elsif ($line=~/spectrum, parent monoisotopic mass error plus/) {
				($infoSubmit{'h:Peptide tolerance'})=($line=~/value="(.+)"/);
			}
			elsif ($line=~/spectrum, parent monoisotopic mass error units" value="(.+)"/) {
				$infoSubmit{'h:Peptide tolerance'}.=" $1";
			}
			elsif ($line=~/spectrum, fragment monoisotopic mass error" value="(.+)"/) {
				$infoSubmit{'i:Fragment tolerance'}="$1";
			}
			elsif ($line=~/spectrum, fragment monoisotopic mass error units" value="(\D+)"/) {
				$infoSubmit{'i:Fragment tolerance'}.=" $1";
			}
			elsif ($line=~/modelling, total spectra used/) {
				($infoSubmit{'o:Queries'})=($line=~/value=\"(\d+)\"/);
			}
			elsif ($line=~/process, version/) {
				($infoSubmit{'c:X! Tandem version'})=($line=~/value=\"X! Tandem Jackhammer TPP \((\d+\.\d+\.\d+\.\d) - LabKey, Insilicos, ISB\) \"/);
			}
			elsif ($line=~/process, start time" value="(.+)"/) {
				my @date=split(/:/,$1);
				$infoSubmit{'r:Search date'}="$date[2]/$date[1]/$date[0] $date[3]:$date[4]:$date[5]";
			}
            elsif ($line=~/list path, sequence source #1" value="(.+)"/) {
				push @{$infoSubmit{'d:Databank'}{'b:File'}},$1;
			}
			elsif ($line=~/aminoacid_modification aminoacid="(\D)" massdiff="(\d*\.?\d*)"/) {
				my $modif=$dbh->selectrow_array("SELECT PSI_MS_NAME FROM MODIFICATION WHERE MONO_MASS LIKE '%$2%'");
				$infoSubmit{'g:Variable modifications'}.=", " if $infoSubmit{'g:Variable modifications'};
				$infoSubmit{'g:Variable modifications'}.="$modif ($1)";
			}
			elsif ($line=~/terminal_modification terminus="n" massdiff="(\d*\.?\d*)"/) {
				my $modif=$dbh->selectrow_array("SELECT PSI_MS_NAME FROM MODIFICATION WHERE MONO_MASS LIKE '%$2%'");
				$infoSubmit{'g:Variable modifications'}.=", " if $infoSubmit{'g:Variable modifications'};
				$infoSubmit{'g:Variable modifications'}.="$modif (N-term)";
			}
			elsif ($line=~/enzymatic_search_constraint enzyme="(\D+)" max_num_internal_cleavages/) {
				$infoSubmit{'d:Enzyme'}=$1;
			}
		}
		my $sthDbID=$dbh->prepare("SELECT ID_DATABANK FROM ANALYSIS_DATABANK WHERE ID_ANALYSIS=$anaID");
		$sthDbID->execute;
		while (my $dbID=$sthDbID->fetchrow_array) {
            my ($dbName,$numEntry)=$dbh->selectrow_array("SELECT NAME,NUM_ENTRY FROM DATABANK WHERE ID_DATABANK=$dbID");
			push @{$infoSubmit{'d:Databank'}{'a:Name'}},$dbName if $dbName;
			push @{$infoSubmit{'d:Databank'}{'c:Sequences in databank'}},$numEntry;
        }
    }
	elsif($fileFormat=~/TDM.XML/){
		open (FILE, $searchFile) || ($infoSubmit{'Error'}="Unable to open Search file");
		while (my $line=<FILE>) {
			if ($line=~/scoring, maximum missed cleavage sites">(\d)<\/note>/) {
				$infoSubmit{'e:Enzyme misscleavage'}=$1;
			}
			elsif ($line=~/spectrum, parent monoisotopic mass error plus">(.+)<\/note>/) {
				$infoSubmit{'h:Peptide tolerance'}=$1;
			}
			elsif ($line=~/spectrum, parent monoisotopic mass error units">(.+)<\/note>/) {
				$infoSubmit{'h:Peptide tolerance'}.=" $1";
			}
			elsif ($line=~/spectrum, fragment monoisotopic mass error">(.+)<\/note>/) {
				$infoSubmit{'i:Fragment tolerance'}="$1";
			}
			elsif ($line=~/spectrum, fragment monoisotopic mass error units">(.+)<\/note>/) {
				$infoSubmit{'i:Fragment tolerance'}.=" $1";
			}
			elsif ($line=~/modelling, total spectra used">(\d+)<\/note>/) {
				$infoSubmit{'o:Queries'}=$1;
			}
			elsif ($line=~/process, start time">(\d+):(\d+):(\d+):(\d+):(\d+):(\d+)<\/note>/) {
				$infoSubmit{'r:Search date'}="$3/$2/$1 $4:$5:$6";
			}
			elsif ($line=~/process, version">X! Tandem Jackhammer TPP \((\d+\.\d+\.\d+\.\d) - LabKey, Insilicos, ISB\) <\/note>/) {
				$infoSubmit{'c:X! Tandem version'}=$1;
			}
            elsif ($line=~/list path, sequence source #1">(.+)<\/note>/) {
				push @{$infoSubmit{'d:Databank'}{'b:File'}},$1;
			}
			elsif($line=~/scoring, include reverse">(\D+)<\/note>/){
				$infoSubmit{'p:Decoy search'}=$1;
			}
			elsif ($line=~/residue, modification mass">(\d*\.?\d*)@(\D),?(\d*\.?\d*)@?(\D?)<\/note>/ || $line=~/refine, modification mass \d">(\d*\.?\d*)@(\D)<\/note>/) {
				my $mass=sprintf "%.4f",$1;
				my $modif=$dbh->selectrow_array("SELECT PSI_MS_NAME FROM MODIFICATION WHERE MONO_MASS LIKE '%$mass%'");
				$infoSubmit{'f:Fixed modifications'}.=", " if $infoSubmit{'f:Fixed modifications'};
				$infoSubmit{'f:Fixed modifications'}.="$modif ($2)";
				if ($3) {
                    my $mass=sprintf "%.4f",$3;
					my $modif=$dbh->selectrow_array("SELECT PSI_MS_NAME FROM MODIFICATION WHERE MONO_MASS LIKE '%$mass%'");
					$infoSubmit{'f:Fixed modifications'}.=", " if $infoSubmit{'f:Fixed modifications'};
					$infoSubmit{'f:Fixed modifications'}.="$modif ($4)";
                }
			}
			elsif($line=~/refine, potential N-terminus modifications">(\d+\.?\d+)<\/note>/){
				my $mass=sprintf "%.4f",$1;
				my $modif=$dbh->selectrow_array("SELECT PSI_MS_NAME FROM MODIFICATION WHERE MONO_MASS LIKE '%$mass%'");
				$infoSubmit{'g:Variable modifications'}.=", " if $infoSubmit{'g:Variable modifications'};
				$infoSubmit{'g:Variable modifications'}.="$modif (N-term)";
			}
			elsif ($line=~/residue, potential modification mass">(\d+\.?\d+)@(\D),?(\d*\.?\d*)@?(\D?)<\/note>/ ||$line=~/refine, potential modification mass \d">(\d+\.?\d+)@(\D)<\/note>/) {
				my $mass=sprintf "%.4f",$1;
				my $modif=$dbh->selectrow_array("SELECT PSI_MS_NAME FROM MODIFICATION WHERE MONO_MASS LIKE '%$mass%'");
				$infoSubmit{'g:Variable modifications'}.=", " if $infoSubmit{'g:Variable modifications'};
				$infoSubmit{'g:Variable modifications'}.="$modif ($2)";
				if ($3) {
                    my $mass=sprintf "%.4f",$3;
					my $modif=$dbh->selectrow_array("SELECT PSI_MS_NAME FROM MODIFICATION WHERE MONO_MASS LIKE '%$mass%'");
					$infoSubmit{'g:Variable modifications'}.=", " if $infoSubmit{'g:Variable modifications'};
					$infoSubmit{'g:Variable modifications'}.="$modif ($4)";
                }
			}
		}
		my $sthDbID=$dbh->prepare("SELECT ID_DATABANK FROM ANALYSIS_DATABANK WHERE ID_ANALYSIS=$anaID");
		$sthDbID->execute;
		while (my $dbID=$sthDbID->fetchrow_array) {
            my ($dbName,$numEntry)=$dbh->selectrow_array("SELECT NAME,NUM_ENTRY FROM DATABANK WHERE ID_DATABANK=$dbID");
			push @{$infoSubmit{'d:Databank'}{'a:Name'}},$dbName if $dbName;
			push @{$infoSubmit{'d:Databank'}{'c:Sequences in databank'}},$numEntry;
        }
	}
	else { # all other cases
		##<Proteome Discoverer
		if ($fileFormat =~ /\.PDM/) {# Get the Proteome-Discoverer's version
			(my $msfFileName = $fileName) =~ s/\_[0-9]*\.*[0-9]*\.pdm/\.msf/;
			my $msfFile=($validStatus==2)? "$promsPath{peptide}/proj_$projectID/ana_$anaID/$msfFileName" : "$promsPath{valid}/multi_ana/proj_$projectID/$msfFileName";
			if (-e $msfFile) {
				my $dbsqlite = DBI->connect( "dbi:SQLite:$msfFile", "", "", {PrintError => 1,RaiseError => 1});
				my ($softVersion)=$dbsqlite->selectrow_array("SELECT SoftwareVersion FROM SchemaInfo ORDER BY rowid ASC LIMIT 1");
				$dbsqlite->disconnect;
				$infoSubmit{'c:Proteome Discoverer version'} = $softVersion;
			}
		}
		
		open (FILE, $searchFile) || ($infoSubmit{'Error'}="Unable to open Search file");
		my $section='';
		my $modKey;
		while (my $line=<FILE>) {
			if ($line=~/[Mascot|peakListPhenyx]; name="(\w+)"/) {
				$section=$1;
				next;
			}
			#last if (($fileFormat eq 'MASCOT.DAT' && $line=~/Mascot; name="summary"/) || ($fileFormat eq 'PHENYX.XML' && $line=~/peakListPhenyx; name="query/));
			last if (($fileFormat eq 'MASCOT.DAT' && $section eq 'summary') || ($fileFormat eq 'PHENYX.XML' && $section=~/^query/));
			next if ($fileFormat eq 'MASCOT.DAT' && $section && $section !~/parameters|header|quantitation/);
			$line=~s/\s+\Z//; # delete \n and trailing spaces
			#>parameters section
			if ($section eq 'parameters') {
				if ($line=~ /^COM=(.+)/) {
					$infoSubmit{'a:Search title'} = $1;
				}
				elsif ($line=~ /^MODS=(.+)/) {
					$infoSubmit{'f:Fixed modifications'} = $1; # also updated in quantitation section (MASCOT.DAT)
				}
				elsif ($line=~ /^IT_MODS=(.+)/) {
					$infoSubmit{'g:Variable modifications'} = $1; # also updated in quantitation section (MASCOT.DAT)
				}
				elsif ($line=~ /^CLE=(.+)/) {
					$infoSubmit{'d:Enzyme'} = $1;
				}
				elsif ($line=~ /^PFA=(\d+)/) {
					$infoSubmit{'e:Enzyme misscleavage'} = $1;
				}
				elsif ($line=~ /^DB(\d*)=(.+)/) {
					push @{$infoSubmit{'l:Databank'}{'a:Name'}},$2;
				}
				elsif ($line=~ /^INSTRUMENT=(.+)/ ){
					$infoSubmit{'b:Instrument'} = $1;
				}
				elsif ($line=~ /^TOL=(.+)/ ){
					$infoSubmit{'h:Peptide tolerance'} = $1;
				}
				elsif ($line=~ /^TOLU=(.+)/) {
					$infoSubmit{'h:Peptide tolerance'} .= " $1";
				}
				elsif ($line=~ /^ITOL=(.+)/) {
					$infoSubmit{'i:Fragment tolerance'} = $1;
				}
				elsif ($line=~ /^ITOLU=(.+)/ ){
					$infoSubmit{'i:Fragment tolerance'} .= " $1";
				}
				elsif ($line=~ /^TAXONOMY=[\. ]*(.+)/ ){
					$infoSubmit{'m:Taxonomy'} = $1;
				}
				elsif ($line=~ /^QUANTITATION=(.*)/ ){
					$infoSubmit{'q:Quantification'} = $1; # also updated in quantification section (MASCOT.DAT)
				}
				elsif ($line=~ /^DECOY=(.*)/ ){
					$infoSubmit{'p:Decoy search'} = ($1)? 'Yes' : 'No';
				}
			}
			if ($section eq 'quantitation') {
				if ($line=~ /<modifications mode="([^"]+)"/) {
					$modKey=($1 eq 'fixed')? 'f:Fixed modifications' : 'g:Variable modifications';
				}
				elsif ($line=~ /<mod_file>(.+)<\/mod_file>/) {
					my $quantifMod=$1;
					my $quoteM=quotemeta($quantifMod);
					if (!$infoSubmit{$modKey} || $infoSubmit{$modKey} !~ /$quoteM/) {
						$infoSubmit{$modKey} .= ',' if $infoSubmit{$modKey};
						$infoSubmit{$modKey} .= $quantifMod;
					}
				}
				elsif ($line=~ /<component name="([^"]+)">/) {
					$infoSubmit{'q:Quantification'}.=":$1";
				}
			}
			#>"header" section
			else {
				if ($line=~ /^sequences(\d*)=(.+)/) { # multiple databanks
					my $dbIdx=($1)? $1-1 : 0;
					$infoSubmit{'l:Databank'}{'d:Sequences in databank'}[$dbIdx]=$2;
				}
				#elsif ($line=~ /^sequences_after_tax\d*=(.+)/ && (!$infoSubmit{'m:Taxonomy'} || ($infoSubmit{'m:Taxonomy'} && $infoSubmit{'m:Taxonomy'} ne 'All entries'))) {
				#	$infoSubmit{'n:Sequences (Taxonomy)'} = $1;
				#}
				elsif ($line=~ /^sequences_after_tax(\d*)=(.+)/) {
					my $dbIdx=($1)? $1-1 : 0;
					$infoSubmit{'l:Databank'}{'e:Sequences after taxonomy filter'}[$dbIdx]=$2;
				}
				elsif ($line=~ /^date=(.+)/) {
					$infoSubmit{'r:Search date'} = strftime("%d/%m/%Y %H:%M:%S",localtime($1));
				}
				elsif ($line=~ /^fastafile(\d*)=(.+)/) { # no \d is not indicative of single-db search!!!
					my $dbIdx=($1)? $1-1 : 0;
					$infoSubmit{'l:Databank'}{'b:File'}[$dbIdx]=$2;
				}
				elsif ($line=~ /^release(\d*)=(.+)/) { # no \d is not indicative of single-db search!!!
					my $dbIdx=($1)? $1-1 : 0;
					$infoSubmit{'l:Databank'}{'c:Release'}[$dbIdx]=$2;
				}
				elsif ($line=~ /^queries=(\d+)/) {
					$infoSubmit{'o:Queries'} = $1;
				}
				elsif ($line=~ /^version=(.+)/) {
					$infoSubmit{'c:Mascot version'} = $1;
				}
			}
		}
		close FILE;
	}

	return %infoSubmit;
}

############################################################
####<Get reporter info (should be in a specific file ?)>####
############################################################
sub getReporterInfo  {
	my $quantiMethodID = $_[0] ;
	my %quantiInfo ;
	my $dbh = $_[1];
	my $sth = $dbh->prepare("SELECT NAME, EXACT_MASS FROM QUANTIF_ELEMENT WHERE ID_QUANTIF_METHOD = $quantiMethodID") ;
	$sth-> execute () ;
	while (my ($elementName, $masse) = $sth->fetchrow_array ()) {
		$quantiInfo{$elementName}{"daValue"} =  $masse;
	}
	return %quantiInfo ;
}

############################################################
####<getDeltaReporter :return delta from reporter match>####
############################################################
sub getDeltaReporter {
	my $dbh = $_[0] ;
	my $instrument = $_[1] ;
	my $user = $_[2] ;
	my $deltaReporterFragment ;
	my ($usrParam,$rules, $deltFrag, $tolFrag, $nrLevel) = &getInstrParam ($dbh,$instrument,$user) ;
	if ($usrParam < 0){ #$usrParam = -1 : instrument not defined; 0 : defined for all user, 1 user defined
		$deltaReporterFragment = 0.3 ;
	}
	else {
		if ($deltFrag eq 'Da') {
			$deltaReporterFragment = $tolFrag ;
		}
		else {
			$deltaReporterFragment = 0.3 ;
		}
	}
	return $deltaReporterFragment ;
}
#############################################################
####<getInstrParam : return parameter from an instrument>####
#############################################################
sub getInstrParam { #  ($usrParam,$refFragRules,$deltaParent,$deltFrag,$tolFrag,$nrLevel) = getInstrParam ($dbh,$instrument,$userID)
	my $dbh = $_[0] ;
	my $localInstr= $dbh->quote($_[1]);
	my $localUser = $dbh->quote($_[2]);
	my ($rules,$deltaParent,$deltaFrag,$tolFrag,$nrLevel);
	my %fragmentationRules;
	my $usrParam = $dbh->selectrow_array ("SELECT COUNT(*) FROM INSTRUMENT where NAME=$localInstr and USE_STATUS='usr' and update_user=$localUser");
	if ($usrParam == 0) {
		my $defInstr = $dbh->selectrow_array ("SELECT COUNT(*) FROM INSTRUMENT where NAME=$localInstr and USE_STATUS='yes'");
		if ($defInstr>0) {
			($rules,$deltaParent,$deltaFrag,$tolFrag,$nrLevel) = $dbh->selectrow_array ("SELECT RULES,DELTA_PARENT,DELTA_FRAGMENT,TOL_FRAGMENT,NR_LEVEL FROM INSTRUMENT where NAME=$localInstr and USE_STATUS='yes'");
		}
		else {
			$usrParam=-1;
			%fragmentationRules=&promsConfig::getFragmentationRules($localInstr); # added by PP
			$deltaParent=&promsConfig::getDeltaParent($localInstr); # added by PP
			$deltaFrag=&promsConfig::getDeltaFragment($localInstr); # added by PP
			$tolFrag=&promsConfig::getFragmentTolerance($localInstr); # added by PP
			$nrLevel=&promsConfig::getNRLevel($localInstr); # added by PP
		}
	}
	else{
		($rules,$deltaParent,$deltaFrag,$tolFrag,$nrLevel) = $dbh->selectrow_array ("SELECT RULES,DELTA_PARENT,DELTA_FRAGMENT,TOL_FRAGMENT,NR_LEVEL FROM INSTRUMENT where NAME=$localInstr and USE_STATUS='usr' and update_user=$localUser");
	}
	if ($rules) { # $usrParam >= 0
		foreach my $line (split(',',$rules)) {
			my ($fragmentType,$value)=split('=',$line);
			$fragmentationRules{$fragmentType}=$value;
		}
	}

	return ($usrParam,\%fragmentationRules,$deltaParent,$deltaFrag,$tolFrag,$nrLevel);
}



############<<< Navigation tree functions & drawing  2.1.2 >>>#################
#############################################
####<<<Write JavaScript tree functions>>>####
#############################################
# $refTree->(depth,type,ID,labelClass,imageClass,selectable,name,popup,refChild1,...2,...3); array
# $refOptions->('STATUS'=>'depth:visibility:operture,...,...'; string   !!! Obsolete !!!
#				'SELECTION'=>(type,ID); array
#               'TITLE'=>'HTML-formatted text'; string
#               'BUTTON'=>'Exp/Collapse buttons style'; string
#               'MULTISELECT'=>0/1; # integer
#               'CHECKBOX'=>('item'=>0/1,...); hash
#               'CHECKED'=>('item:itemID'=>0/1,...); hash
#               'DISCHECK'=>('item:itemID'=>0/1,...); hash
#				'AUTOSAVE'=>'JS variable where tree status is saved'; string
sub writeJsTreeFunctions {
	my ($refTree,$refOptions)=@_;
	my %promsPath=&promsConfig::getServerInfo;

	###<Unvariable functions>###
	print qq
|//--->Creating JS array of tables<---//
var tableArray=new Array(); // array of pointers to tableInfo arrays
function createTableTree() {
	var selIndex;
	var tableInfo=new Array();
	// [0] visible/hidden
	// [1] open/closed
	// [2] depth
	// [3] item type
	// [4] id branch (type:id)
	// [5] labelClass ('' or null=>'item')
	// [6] selectable (0/1)
	// [7...n] list of children index
|;
	my $tableIndex=0;
	&createJsTableArray($refTree,\$tableIndex,$refOptions); # sets selIndex
	my $autoSaveCommand=($refOptions->{'AUTOSAVE'})? $refOptions->{'AUTOSAVE'}.'=writeTreeStatus();' : '/* no auto-save command */';
	print qq
|	//Defines which item is selected at beginning & if tree should be open to make it visible
	if (selIndex==null) selIndex=0; // in case starting selected item was not found
	selectItem(selIndex);
	if (tableArray[selIndex][0]=='none') {showSelectedItem(selIndex);} // make selected item visible
}
//--->Functions to manipulate tree<---//
function click2OpenClose(tabIndex,action,saveFlag) {
	if (tableArray[tabIndex].length==7) {return;} // item has no children
	if (!action) {action='auto';} // open,close,auto
	var iconImg=document.getElementsByName('iconImg')[tabIndex];
	var childVis;
	if (tableArray[tabIndex][1]==0 && action != 'close') { // closed -> open
		tableArray[tabIndex][1]=1;
		document.getElementsByName('expandImg')[tabIndex].src="$promsPath{images}/minus2.gif";
		childVis='block';
		if (document.getElementsByName('iconImg')[tabIndex].src.match('folder')) {
			document.getElementsByName('iconImg')[tabIndex].src="$promsPath{images}/open_folder.gif";
		}
	}
	else if (tableArray[tabIndex][1]==1 && action != 'open') { // open -> closed
		tableArray[tabIndex][1]=0;
		document.getElementsByName('expandImg')[tabIndex].src="$promsPath{images}/plus.gif";
		childVis='none';
		if (document.getElementsByName('iconImg')[tabIndex].src.match('folder')) {
			document.getElementsByName('iconImg')[tabIndex].src="$promsPath{images}/closed_folder.gif";
		}
	}
	// otherwise do nothing
	if (childVis) updateChildVisibility(tabIndex,childVis);
	if (saveFlag) {$autoSaveCommand}
}
function updateChildVisibility(tabIndex,childVis) {
	for (var i=7; i<tableArray[tabIndex].length; i++) { // looping through all children
		var childIndex=tableArray[tabIndex][i];
		tableArray[childIndex][0]=childVis; // update child visibility
		document.getElementById(childIndex).style.display=childVis;
		if (tableArray[childIndex][1]==1) { // update grandchildren if child is open
			updateChildVisibility(childIndex,childVis);
		}
	}
}
function updateTableTree(tabIndex,startIndex,childVis) {
	if (childVis=='block') { // expand
		if (tableArray[tabIndex].length > 7) { // item has children
			document.getElementsByName('expandImg')[tabIndex].src="$promsPath{images}/minus2.gif";
			//if (tableArray[tabIndex][3]=='folder') {
			if (document.getElementsByName('iconImg')[tabIndex].src.match('closed_folder')) {
				document.getElementsByName('iconImg')[tabIndex].src="$promsPath{images}/open_folder.gif";
			}
			tableArray[tabIndex][1]=1; // closed -> open
		}
		else {tableArray[tabIndex][1]=0;} // no child => closed
	}
	else { // collapse
		if (tableArray[tabIndex].length > 7) { // item has children
			document.getElementsByName('expandImg')[tabIndex].src="$promsPath{images}/plus.gif";
			//if (tableArray[tabIndex][3]=='folder') {
			if (document.getElementsByName('iconImg')[tabIndex].src.match('open_folder')) {
				document.getElementsByName('iconImg')[tabIndex].src="$promsPath{images}/closed_folder.gif";
			}
		}
		tableArray[tabIndex][1]=0; // open -> closed
	}
	if (tabIndex > startIndex) { // starting item visibility is not changed
		tableArray[tabIndex][0]=childVis;
	}
	document.getElementById(tabIndex).style.display=tableArray[tabIndex][0];

	for (var i=7; i<tableArray[tabIndex].length; i++) { // looping through all children
		updateTableTree(tableArray[tabIndex][i],startIndex,childVis);
	}
}
function showSelectedItem(selItem) { // ignores unmatched item
	var selIndex=getItemIndex(selItem);
	if (selIndex==null) return; // in case selItem is not found
	var path2Item=new Array();
	// Recording list of relevant tableArray indexes
	path2Item[0]=0; // root
	for (var i=1; i<selIndex; i++) {
		path2Item[tableArray[i][2]]=i; // depth
	}
	path2Item.length=tableArray[selIndex][2]; // truncate path if necessary (length = selIndex depth)
	// Opening tree	to show selected item
	for (var i=0; i<path2Item.length; i++) {
		click2OpenClose(path2Item[i],'open',0);
	}
}
function getItemIndex(selItem) { // branchID or index
	var selIndex;
	selItem+=''; // To string
	if (selItem.match(/:/)) { // branchID was passed
		for (var i=0; i<tableArray.length; i++) {
			if (tableArray[i][4]==selItem) {
				selIndex=i;
				break;
			}
		}
	}
	else {selIndex=selItem;} // index was passed
	return selIndex;
}
function getItemParent(tabIndex) {
	if (tabIndex==null) {tabIndex=selectedTab;}
	if (tabIndex==0) {return null;}
	var i=tabIndex-1;
	while (i>=0) {
		if (tableArray[i][2] < tableArray[tabIndex][2]) {break;} // compare depths
		i--;
	}
	return tableArray[i][4];
}

//--- Tree Status management ---//
function applyTreeStatus(treeStatus) { // type1:id1:id2:...,type2:id1:id2:...,typeX:id1:... (types can be repeated)
	if (!treeStatus) return;
	var dataBlocks=treeStatus.split(',');
	for (var d=0; d<dataBlocks.length; d++) {
		var typeData=dataBlocks[d].split(':');
		for (var i=1; i<typeData.length; i++) {
			showSelectedItem(typeData[0]+':'+typeData[i]); //branchID
		}
	}
}
function writeTreeStatus() {
	var statusObject=new Object();
	var deepestBranchID=tableArray[0][4];
	var prevDepth=0;
	var prevTendance=0;
	var newBranch=0;
	for (var i=1; i<tableArray.length; i++) { // looping through items in tree
		if (tableArray[i][0]=='none') continue; // skip all hidden items
		var tendance=tableArray[i][2]-prevDepth;
		if (tendance > 0) { // 1st child
			deepestBranchID=tableArray[i][4];
			newBranch=1; // scanning new branch
		}
		else if (tendance < 0 && prevTendance >= 0) { // record deepest 1st seebling
			var branchArray=deepestBranchID.split(':');
			if (!statusObject[branchArray[0]]) statusObject[branchArray[0]]=new Array();
			statusObject[branchArray[0]].push(branchArray[1]);
			newBranch=0;
		}
		prevTendance=tendance;
		prevDepth=tableArray[i][2];
	}
	if (newBranch) {
		var branchArray=deepestBranchID.split(':');
		if (!statusObject[branchArray[0]]) statusObject[branchArray[0]]=new Array();
		statusObject[branchArray[0]].push(branchArray[1]);
	}
	//To string
	var statusArray=new Array();
	for (var type in statusObject) {
		var tmpStrg=type+':'+statusObject[type].join(':');
		statusArray.push(tmpStrg);
	}
//alert('Status= '+statusArray.join(','));
	return statusArray.join(',');
}
function addItemToTree(newItems) { // string 'type:id1:id2,type2:id1:id2...'
	var treeStatus=writeTreeStatus();
	return treeStatus+','+newItems;
}
function deleteItemFromTree(selItem) {
	var selIndex;
	if (selItem==null) {selIndex=selectedTab;}
	else {
		selIndex=getItemIndex(selItem);
		if (selIndex==null) return; // in case selItem is not found
	}
	var brotherIndex;
	if (tableArray[selIndex-1][2]==tableArray[selIndex][2]) {brotherIndex=selIndex-1;} // upward brother
	else if (tableArray[selIndex+1] && tableArray[selIndex+1][2]==tableArray[selIndex][2]) {brotherIndex=selIndex+1;} // downward brother
	var parentIndex;
	for (var i=selIndex-1; i>=0; i--) {
		if (tableArray[i][2] < tableArray[selIndex][2]) { // parent was reached
			parentIndex=i;
			break;
		}
	}
	var treeStatus=writeTreeStatus();
	treeStatus+=',';
	treeStatus+=(brotherIndex)? tableArray[brotherIndex][4] : tableArray[parentIndex][4];

	return [treeStatus,tableArray[parentIndex][4]]; // new tree status + parent Item
}
|;

	###<Checkbox management functions>###
	if ($refOptions->{'CHECKBOX'}) {
		print qq
|function checkMyBox(tabIndex) {
	var itemBoxes=document.treeForm.checkIndex;
	var myBox=(itemBoxes.length)? itemBoxes[tabIndex] : itemBoxes;
		if (myBox.style.display=='none' \|\| myBox.disabled) return;
		if (myBox.checked) {myBox.checked=false;}
		else {
			myBox.checked=true;
			actionOnCheck(tabIndex);
		}

}
function updateCheckBoxes(tabIndex) {
	var lastTested=tabIndex;
	var itemBoxes=document.treeForm.checkIndex;
	var chkStatus=itemBoxes[tabIndex].checked;
	//Uncheck all parents
	if (chkStatus==false) uncheckMyParents(tabIndex);
	//Has children? => update all children
	if (tableArray[tabIndex].length>7) {
		for (var i=tabIndex+1;i<tableArray.length;i++) {
			if (tableArray[i][2]<=tableArray[tabIndex][2]) {break;} // brother or parent is reached
			if (itemBoxes[i].disabled) { // do not update a disabled childbox
				if (chkStatus) uncheckMyParents(i); // uncheck all parents
				continue;
			}
			itemBoxes[i].checked=chkStatus;
			lastTested=i;
		}
	}
	// Action performing function
	actionOnCheck(tabIndex);
	return lastTested; // needed in some cases
}
function uncheckMyParents(tabIndex) {
	var parentDepth=tableArray[tabIndex][2]-1;
	var itemBoxes=document.treeForm.checkIndex;
	for (var i=tabIndex-1;i>=0;i--) {
		if (tableArray[i][2]>parentDepth) {continue;}
		itemBoxes[i].checked=false;
		parentDepth--;
	}
}
function getCheckedItems() {
	//Processing tree checkboxes
	var refDepth=0;
	var skipChildren=0;
	var skipSubFrame=0;
	var checkedItemList=new Array();
|;
	foreach my $item (keys %{$refOptions->{'CHECKBOX'}}) {
		print "\tcheckedItemList['$item']=new Array();\n";
	}
	print qq
|	for (var i=0;i<tableArray.length;i++) {
		if (skipChildren==1 && tableArray[i][2]>refDepth) {
			continue;
		}
		skipChildren=0;
		if (document.treeForm.checkIndex[i].checked==true) {
			skipChildren=1;
//alert(tableArray[i][4]);
			var branch=tableArray[i][4].split(':');
			checkedItemList[branch[0]].push(branch[1]);
		}
		refDepth=tableArray[i][2];
	}
	// storing checked items in string
	var tempArray=new Array();
	for (var item in checkedItemList) {
		if (checkedItemList[item].length==0) {continue;}
		tempArray.push(item+':'+checkedItemList[item].join(','));
	}
	var itemString=tempArray.join('+');
	return itemString;
}
function applyCheckedStatus(checkedStrg) {
	var itemBoxes=document.treeForm.checkIndex;
	var categoryList=checkedStrg.split('+');
	for (var c=0; c<categoryList.length; c++) {
		var itemInfo=categoryList[c].split(':'); //item:id1,id2,...
		var idList=itemInfo[1].split(',');
		for (var d=0; d<idList.length; d++) {
			for (var i=0; i<itemBoxes.length; i++) {
				if (itemBoxes[i].disabled \|\| itemBoxes[i].style.display=='none' \|\| itemBoxes[i].checked) continue;
				if (tableArray[i][4]==itemInfo[0]+':'+idList[d]) itemBoxes[i].checked=true;
			}
		}
	}
}
|;
	}

	###<Item(s) selection functions>###
	if ($refOptions->{'MULTISELECT'}) {
		##<Multiple selectable items
		print qq
|function expandCollapseBranch(what,childVis) {
	if (what=='all') {updateTableTree(0,0,childVis);}
	else {
		for (var i=0; i<selectedItems.length; i++) {
			updateTableTree(selectedItems[i],selectedItems[i],childVis);
		}
	}
	$autoSaveCommand
}
//--->Function to select an item in the tree<---//
function selectItem(tabIndex) {
	if (tableArray[tabIndex][6]==0) {return;} // not selectable
	var selectedId=tableArray[tabIndex][4];
	if (tableArray[tabIndex][3]==selectedType) { // same item type selected
		// check if item is already selected
		if (document.getElementById(selectedId).className=='selItem') { // item is already selected -> unselect
			if (selectedItems.length>1) { // deselect only if more than 1 selected
				document.getElementById(selectedId).className=tableArray[tabIndex][5]; // default='item'
				var tempSelItems=new Array();
				for (var i=0; i<selectedItems.length; i++) {
					if (selectedItems[i] != tabIndex) {
						tempSelItems[tempSelItems.length]=selectedItems[i];
					}
				}
				selectedItems=tempSelItems;
			}
		}
		else { // select item
			document.getElementById(selectedId).className='selItem';
			selectedItems[selectedItems.length]=tabIndex;
			// Action performing function
			actionOnSelect(tabIndex);
		}
	}
	else { // new item type selected
		// unselect all previously selected items
		for (var i=0; i<selectedItems.length; i++) {
			document.getElementById(tableArray[selectedItems[i]][4]).className=tableArray[selectedItems[i]][5]; // default='item'
		}
		//select new item
		selectedItems.length=0;
		selectedItems[0]=tabIndex;
		document.getElementById(selectedId).className='selItem';
		selectedType=tableArray[tabIndex][3];
		// Action performing function
		actionOnSelect(tabIndex);
	}
}
function getSelectedBranchID(tabIndex) {
	var selectedArray=new Array();
	var selBranchType;
	if (tabIndex==null) { // not 0 ! return all selected branchIDs
		for (var i=0; i<selectedItems.length; i++) {
			var branchData=tableArray[selectedItems[i]][4].split(':');
			if (i==0) selBranchType=branchData[0];
			selectedArray.push(branchData[1]);
		}
		return selBranchType+':'+selectedArray.join(','); // a string !!!
	}
	else {return tableArray[selectedItems[tabIndex]][4];} // a string !!!
}

var selectedType='$refTree->[2]'; // top item type
var selectedItems=new Array();
selectedItems[0]=0;
|;
	}
	else {
		##<Single selectable item
		print qq
|function expandCollapseBranch(what,childVis) {
	if (what=='all') {updateTableTree(0,0,childVis);}
	else {updateTableTree(selectedTab,selectedTab,childVis);}
	$autoSaveCommand
}
//--->Function to select an item in the tree<---//
function selectItem(tabIndex) {
	if (tableArray[tabIndex][6]==0) {return;} // not selectable
	document.getElementById(tableArray[selectedTab][4]).className=tableArray[selectedTab][5]; // default='item'
	var selectedId=tableArray[tabIndex][4];
	document.getElementById(selectedId).className='selItem';
	selectedTab=tabIndex;
	// Action performing function
	actionOnSelect(tabIndex);
}
function getSelectedBranchID() {
	return tableArray[selectedTab][4];
}
var selectedTab=0; // default selection
|;
	}
}

##########################################
####<<<Create JavaScript tree array>>>#### Called by &writeJsTreeFunctions and by itself
##########################################
sub createJsTableArray {
# $refTree->(depth,type,ID,labelClass,imageClass,selectable,name,popup,refChild1,...2,...3); array
	my ($refTree,$refTabIndex,$refOptions)=@_;
	my $myIndex=${$refTabIndex};
	#my ($visCode,$operture)=($refOptions->{'STATUS'} && $refOptions->{'STATUS'}[$myIndex])? (split(/:/,$refOptions->{'STATUS'}[$myIndex]))[1,2] : (0,0);
	print "\tselIndex=$myIndex;\n" if ($refOptions->{'SELECTION'} && $refOptions->{'SELECTION'}[0] eq $refTree->[1] && $refOptions->{'SELECTION'}[1] eq $refTree->[2]); # eq not == because some ids are not always numerical
	my $visibility=($myIndex==0)? 'block' : 'none'; # Top item always visible
	my $depth=$refTree->[0];
	my $labelClass=($refTree->[3])? $refTree->[3] : 'item';
	print qq
|	tableArray[$myIndex]=new Array();
	tableInfo[$depth]=new Array();
	tableInfo[$depth][0]="$visibility"; // visibility
	tableInfo[$depth][1]=0; // operture
	tableInfo[$depth][2]=$depth; // depth
	tableInfo[$depth][3]="$refTree->[1]"; // item type
	tableInfo[$depth][4]="$refTree->[1]:$refTree->[2]"; // branchID
	tableInfo[$depth][5]="$labelClass" // labelClass
	tableInfo[$depth][6]=$refTree->[5]; // selectable
|;
	for (my $i=8;$i<=$#{$refTree};$i++) { # item's children
		${$refTabIndex}++;
		print "\ttableInfo[$depth][tableInfo[$depth].length]=${$refTabIndex};\n";
		&createJsTableArray($refTree->[$i],$refTabIndex,$refOptions);
	}
	print "\ttableArray[$myIndex]=tableInfo[$depth];\n";
}

###################################
####<<<Print navigation tree>>>#### Called to print tree & by itself
###################################
sub printTree {
# $refTree->(depth,type,ID,name,labelClass,imageClass,selectable,popup,refChild1,...2,...3); array
	my ($refTree,$refOptions,$refTabIndex,$lastBrother,@colList)=@_;
	my %iconImg=&promsConfig::getItemIcones;
	my %promsPath=&promsConfig::getServerInfo;
	my $myIndex=${$refTabIndex};
	if ($myIndex==0) {
		print "<FORM name=\"treeForm\">\n" if $refOptions->{'CHECKBOX'};
		print "$refOptions->{TITLE}" if $refOptions->{'TITLE'};
		if ($refOptions->{'BUTTON'}) {
			print qq
|<TABLE align=center cellpadding=1><TR><TH nowrap>
<INPUT type="button" value="Expand" style="$refOptions->{'BUTTON'}" onclick="expandCollapseBranch('item','block')" /><INPUT type="button" value="Collapse" style="$refOptions->{'BUTTON'}" onclick="expandCollapseBranch('item','none')" />
</TH></TR></TABLE>
|;
		}
	}
	#my ($visCode,$operture)=($refOptions->{'STATUS'}[$myIndex])? (split(/:/,$refOptions->{'STATUS'}[$myIndex]))[1,2] : (0,0);
	#my $visibility=($visCode==1)? 'block' : 'none';
	my $visibility=($myIndex==0)? 'block' : 'none'; # Top item always visible
	print "<TABLE id=\"$myIndex\" border=0 cellspacing=0 cellpadding=0 style=\"display:$visibility\"><TR>\n";
	foreach my $colImg (@colList) {
		print "<TD valign=middle><IMG src=\"$promsPath{images}/$colImg\" width=16 height=22 border=0\"></TD>";
	}
	unless ($refTree->[0]==0) {
		my $nodeImg;
		if ($lastBrother) {
			$nodeImg='lastnode.gif';
			push @colList,'space.gif';
		}
		else {
			$nodeImg='node.gif';
			push @colList,'vertline.gif';
		}
		print "<TD valign=middle><IMG src=\"$promsPath{images}/$nodeImg\" border=0></TD>\n";
	}
	my $expandImg;
	if ($#{$refTree}>=8) {$expandImg='plus.gif';}
	else {$expandImg='no_icone.gif';}
	my $chkboxStrg='';
	if ($refOptions->{'CHECKBOX'}) {
		$chkboxStrg="<INPUT type=\"checkbox\" name=\"checkIndex\" value=\"$myIndex\" onclick=\"updateCheckBoxes($myIndex)\"";
		$chkboxStrg.=' style="display:none"' unless $refOptions->{'CHECKBOX'}{$refTree->[1]}; # all checkboxes must be drawn to match table index
		$chkboxStrg.=' checked' if ($refOptions->{'CHECKED'} && $refOptions->{'CHECKED'}{"$refTree->[1]:$refTree->[2]"});
		$chkboxStrg.=' disabled' if ($refOptions->{'DISCHECK'} && $refOptions->{'DISCHECK'}{"$refTree->[1]:$refTree->[2]"});
		$chkboxStrg.='>';
	}
	my $labelClass=($refTree->[3])? $refTree->[3] : ($refTree->[5]==0)? 'noSelItem' : 'item';
	my $popupStrg='';
	if ($refTree->[7]) {
		if ($refTree->[7]=~/^#no processing#/) { # print string as is
			($popupStrg=$refTree->[7])=~s/^#no processing#//;
		}
		else {$popupStrg=" onmouseover=\"popup('$refTree->[7]')\" onmouseout=\"popout()\"";}
	}
	my $itemIcon=$refTree->[1]; $itemIcon.=":$refTree->[4]" if $refTree->[4];
	my $imageStrg="<IMG name=\"iconImg\" src=\"$promsPath{images}/";
	if ($iconImg{$itemIcon}) {$imageStrg.="$iconImg{$itemIcon}\" onclick=\"selectItem($myIndex)\">";}
	else {$imageStrg.="space.gif\" style=\"display:none\">";}
	my $itemID="$refTree->[1]:$refTree->[2]";
	print qq
|<TD valign=middle><IMG name="expandImg" src="$promsPath{images}/$expandImg" onclick="click2OpenClose($myIndex,'auto',1)"></TD>
<TD valign=middle nowrap>$chkboxStrg$imageStrg</TD>
<TH valign=middle align=left nowrap>&nbsp<A class="$labelClass" id="$itemID" href="javascript:selectItem($myIndex)" $popupStrg>$refTree->[6]</A></TH>
</TR></TABLE>
|;
	for (my $i=8;$i<=$#{$refTree};$i++) {
		${$refTabIndex}++;
		my $lastChild=($i==$#{$refTree})? 1 : 0;
		&printTree($refTree->[$i],$refOptions,$refTabIndex,$lastChild,@colList);
	}

	print "</FORM>\n" if ($refOptions->{'CHECKBOX'} && $myIndex==0);

	if ($myIndex==0) {print "<SCRIPT LANGUAGE=\"JavaScript\">createTableTree();</SCRIPT>\n";}

}

###################################
####<<< Prevent Right Click >>>####
###################################
sub noRightClick {}

########################################
####<<< Update validation history>>>####
########################################
# STATUS: 0=not reported, 1=reported, -1=reported BUT action was performed again in validation mode after report
sub updateAnalysisHistory {
	my ($dbh,$analysisID,$paramStrg,$valType,$userID) = @_; # $userID only for PhosphoRS
	$paramStrg = '' unless($paramStrg);
	my ($userName,$numVerifQueries,$numAllQueries,$numValidPeptides,$numDecoyPeptides,$numAllProteins,$numValidProteins);

	if ($valType ne 'endVal') { # $userName not needed & not defined if called from send2Biologist.cgi with $call='cmd'
		$userID = $ENV{'REMOTE_USER'} unless $userID;
		($userName) = $dbh->selectrow_array("SELECT USER_NAME FROM USER_LIST WHERE ID_USER='$userID'") if $userID;
		$userName='Unknown' unless $userName;
	}
	my ($maxStep) = $dbh->selectrow_array("SELECT MAX(STEP) FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID");
	$maxStep = 0 if (!$maxStep or $maxStep<0); # if NULL -> 0
	my ($lastValType,$lastValID,$lastStatus) = $dbh->selectrow_array("SELECT VAL_TYPE, ID_VAL_HISTORY, STATUS FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND STEP=$maxStep AND STATUS>=0");
	$lastValType='' unless $lastValType;

	# TODO: use the following statement for each INSERT
	my $sthNewEntry = $dbh->prepare("INSERT INTO VALIDATION_HISTORY(ID_ANALYSIS,STEP,STATUS,VAL_TYPE,PARAM_STRG,QUERY_VAL_STRG,PEP_VAL_STRG,PROT_VAL_STRG,START_DATE,VALID_USER)
					VALUES(?,?,?,?,?,?,?,?,NOW(),?)");
	my @newParams = ($analysisID, $maxStep+1, 0, $valType, undef, undef, undef, undef, $userName);

	if($valType eq "clear_all"){
		$dbh->do("DELETE FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND STATUS=0 AND VAL_TYPE!='filter' AND VAL_TYPE!='flag' AND VAL_TYPE!='prs'"); # prs (PP 10/08/15)
		$dbh->do("UPDATE VALIDATION_HISTORY SET STATUS=-1 WHERE ID_ANALYSIS=$analysisID AND VAL_TYPE!='filter' AND VAL_TYPE!='flag' AND VAL_TYPE!='prs'"); # prs (PP 10/08/15)
	}
	elsif ($valType eq "clear_auto"){
		$dbh->do("DELETE FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND STATUS=0 AND (VAL_TYPE='comp_s' OR VAL_TYPE='quali')");
		my $maxManual = $dbh->selectrow_array("SELECT MAX(STEP) FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND VAL_TYPE='manual'");
		if($maxManual){
			$dbh->do("UPDATE VALIDATION_HISTORY SET STATUS=-1 WHERE ID_ANALYSIS=$analysisID AND STEP<$maxManual");
			my ($numVerifQueries,$numAllQueries,$numValidPeptides,$numDecoyPeptides,$numAllProteins,$numValidProteins,$numFilteredProteins) = &countAnaStats($dbh,$analysisID);
			my ($idValHistory,$queryValStrg,$pepValStrg,$protValStrg) = $dbh->selectrow_array("SELECT ID_VAL_HISTORY, QUERY_VAL_STRG, PEP_VAL_STRG, PROT_VAL_STRG FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND STEP=$maxManual");
			$queryValStrg =~ s/V_Verif\=\d+;/V_Verif\=$numVerifQueries;/;
			$queryValStrg =~ s/V_All\=\d+;/V_All\=$numAllQueries;/;
			$pepValStrg =~ s/V_Valid\=\d+;/V_Valid\=$numValidPeptides;/;
			$pepValStrg =~ s/V_Decoy\=\d+;/V_Decoy\=$numDecoyPeptides;/;
			$protValStrg =~ s/V_Valid\=\d+;/V_Valid\=$numValidProteins;/;
			$protValStrg =~ s/V_All\=\d+;/V_All\=$numAllProteins;/;
			my $sthUpNbs = $dbh->prepare("UPDATE VALIDATION_HISTORY SET QUERY_VAL_STRG=?, PEP_VAL_STRG=?, PROT_VAL_STRG=? WHERE ID_VAL_HISTORY=$idValHistory");
			$sthUpNbs->execute($queryValStrg,$pepValStrg,$protValStrg);
			$sthUpNbs->finish;
		} else {
			$dbh->do("UPDATE VALIDATION_HISTORY SET STATUS=-1 WHERE ID_ANALYSIS=$analysisID");
		}
	}
	elsif($valType eq "report"){
		$dbh->do("DELETE FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND STATUS=-1");
		$dbh->do("UPDATE VALIDATION_HISTORY SET STATUS=1 WHERE ID_ANALYSIS=$analysisID AND STATUS=0");
		my $sthGetNb = $dbh->prepare("SELECT ID_VAL_HISTORY,VAL_TYPE, QUERY_VAL_STRG, PEP_VAL_STRG, PROT_VAL_STRG
					     FROM VALIDATION_HISTORY
					     WHERE ID_ANALYSIS=$analysisID AND STEP>0 ORDER BY STEP");
		$sthGetNb->execute;
		my $count = 0;
		while (my @lineElement = $sthGetNb->fetchrow_array){
			$count++;
			my $id = $lineElement[0];
			my $valType = $lineElement[1];
			if ($valType =~ /comp/ || $valType eq 'quali' || $valType eq 'manual' || $valType eq 'incexc' || $valType eq 'filter'){
				my @valStrgList;
				foreach my $valStrg ($lineElement[2],$lineElement[3],$lineElement[4]){
					# Copying Validate data (V_) to Reported data (R_)
					# V_ (validated) data are visible my massists, R_ (reported) are visible by biologists
					foreach my $e ( $valStrg =~ /V_(\w+)=\d+;/g){ # fetching data names
						my ($n) = ($valStrg =~ /V_$e=(\d+);/); # fetching corresponding validated values
						if($valStrg =~ /R_$e\=\d+;/){ # update reported values
							$valStrg =~ s/R_$e\=\d+;/R_$e\=$n;/;
						} else { # or append them if don't exist
							$valStrg .= "R_$e=$n;";
						}
					}
					push @valStrgList, $valStrg;
				}
				my $sthUpStrg = $dbh->prepare("UPDATE VALIDATION_HISTORY SET QUERY_VAL_STRG=?, PEP_VAL_STRG=?, PROT_VAL_STRG=? WHERE ID_VAL_HISTORY=$id");
				$sthUpStrg->execute(@valStrgList);
			}
			### Updating step number ###
			$dbh->do("UPDATE VALIDATION_HISTORY SET STEP=$count WHERE ID_VAL_HISTORY=$id");
		}
	}
	elsif ($valType eq "manual" || $valType eq "quali" || $valType eq "comp_s"){
		my ($numVerifQueries,$numAllQueries,$numValidPeptides,$numDecoyPeptides,$numAllProteins,$numValidProteins,$numFilteredProteins) = &countAnaStats($dbh,$analysisID);
		####Updating database###
		if ($valType eq "manual" && $lastValType eq "manual") {
			my $sthUpMan = $dbh->prepare("UPDATE VALIDATION_HISTORY SET QUERY_VAL_STRG=?, PEP_VAL_STRG=?, PROT_VAL_STRG=? WHERE STEP=$maxStep AND ID_ANALYSIS=$analysisID");
			my ($idValHistory,$queryValStrg,$pepValStrg,$protValStrg) = $dbh->selectrow_array("SELECT ID_VAL_HISTORY, QUERY_VAL_STRG, PEP_VAL_STRG, PROT_VAL_STRG FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND STEP=$maxStep");
			$queryValStrg =~ s/V_Verif\=\d+;/V_Verif\=$numVerifQueries;/;
			$queryValStrg =~ s/V_All\=\d+;/V_All\=$numAllQueries;/;
			$pepValStrg =~ s/V_Valid\=\d+;/V_Valid\=$numValidPeptides;/;
			$pepValStrg =~ s/V_Decoy\=\d+;/V_Decoy\=$numDecoyPeptides;/;
			$protValStrg =~ s/V_Valid\=\d+;/V_Valid\=$numValidProteins;/;
			$protValStrg =~ s/V_All\=\d+;/V_All\=$numAllProteins;/;
			$sthUpMan->execute($queryValStrg,$pepValStrg,$protValStrg);
			$sthUpMan->finish;
		} else {
			my $queryValStrg = "V_Verif=".$numVerifQueries.";V_All=".$numAllQueries.";"; # V_ -> validated but non yet reported
			my $pepValStrg = "V_Valid=".$numValidPeptides.";V_Decoy=".$numDecoyPeptides.";";
			my $protValStrg = "V_Valid=".$numValidProteins.";V_All=".$numAllProteins.";";
			#if ($valType eq "quali"){ $paramStrg = join ";", @{$paramList} ; }
			my $sthUp=$dbh->prepare("INSERT INTO VALIDATION_HISTORY (ID_ANALYSIS,STEP,VAL_TYPE,PARAM_STRG,QUERY_VAL_STRG,PEP_VAL_STRG,PROT_VAL_STRG,START_DATE,VALID_USER,STATUS) VALUES (?,?,?,?,?,?,?,NOW(),?,?)");
			$sthUp->execute($analysisID,$maxStep+1,$valType,$paramStrg,$queryValStrg,$pepValStrg,$protValStrg,$userName,0);
			$sthUp->finish;
		}
	}
	elsif ($valType eq "comp_f"){ # flagging comparison, only FlagUp and FlagDown peptides count will be displayed
		my $numFlagUpPeptides= 0; my $numFlagDownPeptides=0;
		my $sthGetIDQueries = $dbh->prepare("SELECT ID_QUERY FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID");
		$sthGetIDQueries->execute;
		my $info_pepStrg = '';
		foreach my $rank (1..10){
			$info_pepStrg .= "INFO_PEP$rank,";
		}
		$info_pepStrg =~ s/,\Z//;
		while(my $IDquery = $sthGetIDQueries->fetchrow_array){
			my @infoPep = $dbh->selectrow_array("SELECT $info_pepStrg FROM QUERY_VALIDATION WHERE ID_QUERY=$IDquery");
			RANK: foreach my $infoPep(@infoPep){
				last RANK unless($infoPep);
				my ($FLT) = ($infoPep =~ /FLT=(-?\d),?/);
				next RANK unless($FLT);
				if($FLT == 1){ $numFlagUpPeptides++ } elsif($FLT == -1){$numFlagDownPeptides++};
			}
		}
		$sthGetIDQueries->finish;
		my $sthUp=$dbh->prepare("INSERT INTO VALIDATION_HISTORY (ID_ANALYSIS,STEP,VAL_TYPE,PARAM_STRG,QUERY_VAL_STRG,PEP_VAL_STRG,PROT_VAL_STRG,START_DATE,VALID_USER,STATUS) VALUES (?,?,?,?,?,?,?,NOW(),?,?)");
		my $pepValStrg = "V_FlagUp=".$numFlagUpPeptides.";V_FlagDown=".$numFlagDownPeptides.";";
		$sthUp->execute($analysisID,$maxStep+1,$valType,$paramStrg,'',$pepValStrg,'',$userName,0);
		$sthUp->finish;
	}
	elsif ($valType eq "r_flag" && $lastValType ne "r_flags"){ # Removing flags
		my $sthUp=$dbh->prepare("INSERT INTO VALIDATION_HISTORY (ID_ANALYSIS,STEP,VAL_TYPE,PARAM_STRG,QUERY_VAL_STRG,PEP_VAL_STRG,PROT_VAL_STRG,START_DATE,VALID_USER,STATUS) VALUES (?,?,?,?,?,?,?,NOW(),?,?)");
		my $pepValStrg = "V_FlagUp=0;V_FlagDown=0;";
		$sthUp->execute($analysisID,$maxStep+1,$valType,$paramStrg,'',$pepValStrg,'',$userName,0);
		$sthUp->finish;
	}
	##>Modified by PP (12/11/04)
	elsif ($valType eq "filter") {
		my ($numVerifQueries,$numAllQueries,$numValidPeptides,$numDecoyPeptides,$numAllProteins,$numValidProteins,$numFilteredProteins) = &countAnaStats($dbh,$analysisID);
		my $protValStrg = "V_Filt=$numFilteredProteins;V_All=$numAllProteins;";
		$maxStep++;
		$dbh->do("INSERT INTO VALIDATION_HISTORY (ID_ANALYSIS,STEP,VAL_TYPE,PARAM_STRG,PROT_VAL_STRG,START_DATE,VALID_USER,STATUS) VALUES ($analysisID,$maxStep,'$valType','$paramStrg','$protValStrg',NOW(),'$userName',0)");
	}
	elsif ($valType eq "r_filt" && $lastValType ne "r_filt"){ # remove protein filter
		my ($numFilter)=$dbh->selectrow_array("SELECT COUNT(*) FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND VAL_TYPE='filter' AND STATUS>=0");
		if ($numFilter==1 && $lastValType eq 'filter') { # just delete filter entry
			my ($status) = $dbh->selectrow_array("SELECT STATUS FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND STEP=$maxStep");
			if($status == 0){
				$dbh->do("DELETE FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND STEP=$maxStep");
			} elsif($status == 1){
				$dbh->do("UPDATE VALIDATION_HISTORY SET STATUS=-1 WHERE ID_ANALYSIS=$analysisID AND STEP=$maxStep");
			}
		}
		else {
			my ($numVerifQueries,$numAllQueries,$numValidPeptides,$numDecoyPeptides,$numAllProteins,$numValidProteins,$numFilteredProteins) = &countAnaStats($dbh,$analysisID);
			#my $protValStrg = "V_Filt=".$numFilteredProteins.";V_All=".$numAllProteins.";";
			#my $sthUp=$dbh->prepare("INSERT INTO VALIDATION_HISTORY (ID_ANALYSIS,STEP,VAL_TYPE,PARAM_STRG,PROT_VAL_STRG,START_DATE,VALID_USER,STATUS) VALUES (?,?,?,?,?,NOW(),?,?)");
			#$sthUp->execute($analysisID,$maxStep+1,$valType,$paramStrg,$protValStrg,$userName,0);
			#$sthUp->finish;
			my $protValStrg = "V_Valid=$numValidProteins;V_All=$numAllProteins;";
			$maxStep++;
			$dbh->do("INSERT INTO VALIDATION_HISTORY (ID_ANALYSIS,STEP,VAL_TYPE,PARAM_STRG,PROT_VAL_STRG,START_DATE,VALID_USER,STATUS) VALUES ($analysisID,$maxStep,'$valType','$paramStrg','$protValStrg',NOW(),'$userName',0)");
		}
	}
	### Manual protein exclusion (or inclusion) ###
	elsif ($valType eq "incexc"){
		my @numList = &countAnaStats($dbh,$analysisID);
		my $numExcluProteins = $numList[7];
		my $protValStrg;
		my ($previousExclu) = $dbh->selectrow_array("SELECT PROT_VAL_STRG FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND VAL_TYPE='incexc' AND STEP<$maxStep AND STATUS>= 0 ORDER BY STEP DESC LIMIT 0,1");
		my $previousExcluProt = 0;
		if($previousExclu){
			($previousExcluProt) = ( $previousExclu =~ /V_Exclu\=(\d+);/);
		}
		if($lastValType eq 'incexc'){ # last entry was already a manual exclusion
			if($previousExcluProt == 0 && $numExcluProteins == 0){ # number of excluded proteins is 0, and was lastly registered as 0 in a previous entry
				my ($lastStatus) = $dbh->selectrow_array("SELECT STATUS FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND STEP=$maxStep");
				if($lastStatus == 0){ # last exclusion was not reported and can be deleted
					$dbh->do("DELETE FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND STEP=$maxStep");
				} elsif($lastStatus == 1){ # last exclusion was reported, status changes to -1
					my $sthUp = $dbh->prepare("UPDATE VALIDATION_HISTORY SET STATUS=? WHERE ID_ANALYSIS=$analysisID AND STEP=$maxStep");
					$sthUp->execute(-1);
					$sthUp->finish;
				}
			} else { # updating excluded proteins number in last entry
				($protValStrg) = $dbh->selectrow_array("SELECT PROT_VAL_STRG FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND STEP=$maxStep");
				my $sthUpIncExc = $dbh->prepare("UPDATE VALIDATION_HISTORY SET PROT_VAL_STRG=? WHERE STEP=$maxStep AND ID_ANALYSIS=$analysisID");
				$protValStrg =~ s/V_Exclu\=\d+;/V_Exclu\=$numExcluProteins;/;
				$sthUpIncExc->execute($protValStrg);
				$sthUpIncExc->finish;
			}
		}
		else{ # last entry is not a manual exclusion -> insert of a new line
			my $numAllProteins = $numList[4];
			my $protValStrg = "V_Exclu=$numExcluProteins;V_All=$numAllProteins;";
			my $sthUpIncExc=$dbh->prepare("INSERT INTO VALIDATION_HISTORY (ID_ANALYSIS,STEP,VAL_TYPE,PROT_VAL_STRG,START_DATE,VALID_USER,STATUS) VALUES (?,?,?,?,NOW(),?,?)");
			$sthUpIncExc->execute($analysisID,$maxStep+1,$valType,$protValStrg,$userName,0);
			$sthUpIncExc->finish;
		}

	}
	### End validation ###
	elsif($valType eq 'endVal'){
		$dbh->do("DELETE FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND STATUS=0");
		$dbh->do("UPDATE VALIDATION_HISTORY SET STATUS=1 WHERE ID_ANALYSIS=$analysisID AND STATUS=-1");
	}

	### PhosphoRS ###
	elsif($valType eq 'prs'){
		#my ($activationType) = ($paramStrg =~ /activationType:([^;]+);/);
		#my ($threshold) = ($paramStrg =~ /threshold:([^;]+);/);
		#my ($massDeviation) = ($paramStrg =~ /massDeviation:([^;]+);/);

		# Checking old PhosphoRS analysis #
		my ($previousStatus,$prsID) = $dbh->selectrow_array("SELECT STATUS,ID_VAL_HISTORY FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND VAL_TYPE=\"$valType\"");
		if(defined $previousStatus){
			if ($previousStatus == 0){
				$dbh->do("DELETE FROM VALIDATION_HISTORY WHERE ID_VAL_HISTORY=$prsID");
			} elsif ($previousStatus == 1) {
				$dbh->do("UPDATE VALIDATION_HISTORY SET STATUS=-1 WHERE ID_VAL_HISTORY=$prsID");
			}
		}

		# Count PRS data #
		my @infoPepList = map { "INFO_PEP$_" } (1..10);
		my $infoPepStrg = join(', ',@infoPepList);
		my $sthInfoPep = $dbh->prepare("SELECT $infoPepStrg FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID AND QUERY_NUM>0");
		my ($countConfirmed,$countChanged,$countFlagged,$countUnchanged,$countTotal) = (0,0,0,0,0);
		$sthInfoPep->execute;
		while(my @infoPep = $sthInfoPep->fetchrow_array){
			while (my $infoPep = shift @infoPep){
				if($infoPep =~ /VMOD=[^,]+Phospho[^,]+,/){
					$countTotal++;
					if($infoPep =~ /PRS=(\d);/){
						if($1 == 3){
							$countConfirmed++;
						} elsif ($1 == 2){
							$countChanged++;
						} elsif ($1 == 1){
							$countFlagged++;
						} elsif( $1 == 0){
							$countUnchanged++;
						}
					}
				}
			}
		}
		$sthInfoPep->finish;
		my $queryValStrg = "total:$countTotal;confirmed:$countConfirmed;changed:$countChanged;flagged:$countFlagged;unchanged:$countUnchanged;";

		my $sthNewPRS = $dbh->prepare("INSERT INTO VALIDATION_HISTORY (ID_ANALYSIS,STEP,VAL_TYPE,PARAM_STRG,QUERY_VAL_STRG,START_DATE,VALID_USER,STATUS) VALUES (?,?,?,?,?,NOW(),?,?)");
		$sthNewPRS->execute($analysisID,-1,$valType,$paramStrg,$queryValStrg,$userName,0);
		$sthNewPRS->finish;
	}

	### Activate lower-score peptides (manual) ###
	elsif($valType eq 'lowP_m'){
		unless ($lastValType eq 'lowP_m') {
			my $sthInsLP = $dbh->prepare("INSERT INTO VALIDATION_HISTORY (ID_ANALYSIS,STEP,VAL_TYPE,START_DATE,VALID_USER,STATUS) VALUES (?,?,?,NOW(),?,?)");
			$sthInsLP->execute($analysisID,$maxStep+1,$valType,$userName,0);
			$sthInsLP->finish;
		}
	}

	### Activate lower-score peptides (all) ###
	elsif($valType eq 'lowP_a'){

		# if last validation step was a low-score peptide manual activation (and unreported), merge it with the new one
		if ($lastValType eq 'lowP_m' && $lastStatus == 0) {
			my $sthUpVH = $dbh->prepare("UPDATE VALIDATION_HISTORY SET VAL_TYPE=? WHERE ID_VAL_HISTORY=?");
			$sthUpVH->execute($valType,$lastValID);
			$sthUpVH->finish;
		}
		elsif ($lastValType ne 'lowP_a'){
			$sthNewEntry->execute(@newParams);
		}
	}

	### Edit varMod ###
	#elsif($valType eq 'varMod'){
	#
	#	my ($vModInHistory) = $dbh->selectrow_array("SELECT COUNT(*) FROM VALIDATION_HISTORY WHERE ID_ANALYSIS=$analysisID AND VAL_TYPE=\"varMod\"");
	#	unless($vModInHistory){
	#		$newParams[1] = -1; # step
	#		$newParams[2] = 1 if $additionalParam; # if vMod change was made on validated data, consider this entry as reported
	#		$sthNewEntry->execute(@newParams);
	#	}
	#}

	### FDR during import ###
	#elsif($valType eq 'FDR'){
	#
	#	$newParams[1] = -1; #step
	#	$newParams[4] = $paramStrg;
	#
	#	$sthNewEntry->execute(@newParams);
	#}

	$sthNewEntry->finish;
	$dbh->commit;
}

sub countAnaStats{ # count queries/peptides/proteins in current state of an analysis
	my ($dbh,$analysisID) = @_;

	####Queries and peptides####
	my $sthQL = $dbh->prepare("SELECT QUERY_NUM,VALID_STATUS FROM QUERY_VALIDATION WHERE ID_ANALYSIS=$analysisID");
	$sthQL->execute;
	my $refQueryList = $sthQL->fetchall_arrayref;
	$sthQL->finish;


	my $numVerifQueries = 0;
	my $numAllQueries = 0;
	my $numValidPeptides = 0;
	my $numDecoyPeptides = 0;
	foreach my $refQuery (@{$refQueryList}){
		if ($refQuery->[0] >0){
			$numAllQueries++;
			$numVerifQueries++ if ($refQuery->[1] != -1);
			$numValidPeptides+=$refQuery->[1] if ($refQuery->[1]>0);
		} elsif ($refQuery->[0]<0){ #Decoy
			$numDecoyPeptides+=$refQuery->[1] if ($refQuery->[1]>0);
		}
	}

	####Proteins####
	my $sthCountProt=$dbh->prepare("SELECT COUNT(ID_PROT_VALID) FROM PROTEIN_VALIDATION WHERE ID_ANALYSIS=$analysisID AND IDENTIFIER NOT LIKE 'DECOY_%' AND SEL_STATUS>=? AND SEL_STATUS<=?");
	$sthCountProt->execute(-3,2);
	my ($numAllProteins) = $sthCountProt->fetchrow_array;
	$sthCountProt->execute(1,2);
	my ($numValidProteins) = $sthCountProt->fetchrow_array;
	$sthCountProt->execute(-3,-3);
	my ($numFilteredProteins) = $sthCountProt->fetchrow_array;
	$sthCountProt->execute(-2,-2);
	my ($numExcluProteins) = $sthCountProt->fetchrow_array;
	$sthCountProt->finish;

	return ($numVerifQueries,$numAllQueries,$numValidPeptides,$numDecoyPeptides,$numAllProteins,$numValidProteins,$numFilteredProteins,$numExcluProteins);
}

###########################################
####<<< AJAX save proteins to Lists >>>####
###########################################
# Ajax calls showProtQuantification.cgi
# Requires:
#	-JS function getXMLHTTP() to be defined
#	-An empty HTML DIV with id="saveProtDIV" to be defined
#   -A button with id="saveFormBUTTON" to display protein-version form
#   -A button with id="saveSiteFormBUTTON" to display site-version form 'optional)
sub printAjaxManageSaveProteins {
	my ($projectID,$refPromsPath,$path2ChkBoxes,$actionOnCreation)=@_;
	$path2ChkBoxes='document.protForm.chkProt' unless $path2ChkBoxes;
	$actionOnCreation='null' unless $actionOnCreation; # Optional callback function if a Theme or List is created (arguments passed are themeSel,listSel,XHR).
	print qq
|
//CHOOSE OR CREATE THEME, LIST and SAVE CHECKED PROTEIN --->
var saveListGlobals={
	themesFetched:false, // TODO: update all calling scripts from themesFetched to saveListGlobals.themesFetched
	listType:'PROT', // global JS variable because of possible swtich between PROT and SITE without new AJAX call for job=getThemes,save,cancel
	modifID:0
};
function ajaxManageSaveProteins(jobStrg,listType,modifID) { // Uses local XHR: No abort of global XHR if defined
	if (!listType) listType=saveListGlobals.listType \|\| 'PROT';
	saveListGlobals.listType=listType;
	var jobData=jobStrg.split(':');
	var job=jobData[0],themeID=jobData[1],listID=jobData[2];
	var formButtonID={PROT:'saveFormBUTTON',SITE:'saveSiteFormBUTTON'};
	var saveFormBut=document.getElementById(formButtonID[listType]);
	if (job=='cancel') {
		document.getElementById('saveProtDIV').style.display='none';
		if (saveFormBut.style.display=='none') {saveFormBut.style.display='block';} else {saveFormBut.style.visibility='visible';}
		return;
	}
	var saveProtButton=document.getElementById('saveProtBUTTON');
	var targetDiv;
	var paramStrg='ACT=ajaxSaveProt&job='+job+'&listType='+listType;
	if (job=='getThemes') { // classifications
		if (!modifID) {modifID=0;}
		saveListGlobals.modifID=modifID;
//console.log('getThemes',saveListGlobals.modifID);
		if (saveFormBut.style.display=='block') {saveFormBut.style.display='none';} else {saveFormBut.style.visibility='hidden';} // display:block indicates button should be removed not hidden
		/* Reset other list Type if exists */
		var otherType=(listType=='SITE')? 'PROT' :'SITE';
		var otherFormBut=document.getElementById(formButtonID[otherType]);
		if (otherFormBut) {
			if (otherFormBut.style.display=='none') {otherFormBut.style.display='block';} else {otherFormBut.style.visibility='visible';}
		}
		var entities={PROT:'proteins',SITE:'sites'};
		targetDiv=document.getElementById('saveProtDIV');
		targetDiv.style.display='block';
		if (saveListGlobals.themesFetched) { // already retrieved
			for (let i=1; i<=3; i++) {
				document.getElementById('listTypeEntity'+i).innerHTML=entities[listType];
			}
			if (document.getElementById('themeSELECT').selectedIndex) { // a theme was selected => reset lists
				//document.getElementById('themeSELECT').selectedIndex=0;
				ajaxManageSaveProteins(document.getElementById('themeSELECT').value,listType);
			}
			return;
		}
		/* Display wait bar */
		targetDiv.innerHTML='<BR><BR><IMG src="$refPromsPath->{images}/scrollbarGreen.gif"><BR><BR>'; // style="width:300px"
		paramStrg+='&projID=$projectID';
	}
	else if (job=='getLists') { // get lists
		targetDiv=document.getElementById('listSelectionDIV');
		var listSel=document.getElementById('listSelectionSELECT');
		var thNameInput=document.getElementById('themeNameIN');
		if (themeID==0) { // no selection
			targetDiv.innerHTML='';
			saveProtButton.disabled=true;
			return;
		}
		else if (themeID==-1) { // new theme
			thNameInput.style.display='block';
			saveProtButton.disabled=false;
		}
		else { // existing theme
			thNameInput.style.display='none';
			saveProtButton.disabled=false;
		}
		targetDiv.innerHTML='<BR><IMG src="$refPromsPath->{images}/scrollbarGreen.gif">';
		paramStrg+='&themeID='+themeID;
	}
	else if (job=='selectList') { // list selection
		document.getElementById('listNameIN').style.display=(listID==-1)? 'block' : 'none';
		saveProtButton.disabled=(listID==0)? true : false;
		return;
	}
	else { // saveProt
		paramStrg+='&modifID='+saveListGlobals.modifID;
		// Checking for Theme
		targetDiv=document.getElementById('saveProtDIV');
		var themeSel=document.getElementById('themeSELECT');
		var themeName=document.getElementById('themeNameIN').value;
		var listSel=document.getElementById('listSelectionSELECT');
		var listName=(listSel)? document.getElementById('listNameIN').value : null;
		if (themeSel.value=='getList:0') { // no Theme selected
			alert('ERROR: No Theme selected!');
			return;
		}
		else if (themeSel.value=='getLists:-1') { // new Theme & List
			if (!themeName \|\| themeName=='New Theme name') { // check Theme name
				alert('Provide a name for new Theme.');
				return;
			}
			else if (!listName \|\| listName=='New List name') { // check List name
				alert('Provide a name for new List.');
				return;
			}
			paramStrg+='&projID=$projectID&themeID=-1&themeName='+encodeURIComponent(themeName)+'&listID=-1&listName='+encodeURIComponent(listName);
		}
		else { // an existing Theme is selected
			var thData=themeSel.value.split(':');
			paramStrg+='&themeID='+thData[1];
			if (listSel.value.match(':0')) { // no List selected
				alert('ERROR: No List selected!');
				return;
			}
			else if (listSel.value.match(':-1')) { // check List name
				if (!listName \|\| listName=='New List name') { // check List name
					alert('Provide a name for new List.');
					return;
				}
				paramStrg+='&listID=-1&listName='+encodeURIComponent(listName);
			}
			else { // an existing List is selected
				var liData=listSel.value.split(':');
				paramStrg+='&listID='+liData[2];
			}
		}
		//List content management
		paramStrg+=(document.getElementById('replace0').checked)? '&replace=0' : '&replace=1';

		// Selected proteins
		var selProtList=[];
		var checkBoxList=$path2ChkBoxes;
		if (checkBoxList.length) {
			var firstHidden=true;
			for (var i=0; i < checkBoxList.length; i++) {
				if (checkBoxList[i].checked) {
					if (checkBoxList[i].className=='hiddenProt') {
						if (firstHidden) {
							firstHidden=false;
							alert('Hidden proteins cannot be stored in Lists');
						}
						continue;
					}
					/*Some values are anaID1:...anadIDx:protID*/
					var anaProtVal=checkBoxList[i].value.split(':');
					if (anaProtVal[anaProtVal.length-1].match('/')) { // modProtID with ambigutity => do not split!
						selProtList.push(checkBoxList[i].value);
					}
					else { // last/only Idx is protID
						selProtList.push(anaProtVal[anaProtVal.length-1]);
					}
				}
			}
		}
		else if (checkBoxList.checked) {
			if (checkBoxList.className=='hiddenProt') {
				alert('Hidden proteins cannot be stored in Lists');
				return;
			}
			/*Some values are anaID1:...anadIDx:protID*/
			var anaProtVal=checkBoxList.value.split(':');
			if (anaProtVal[anaProtVal.length-1].match('/')) { // modProtID with ambigutity => do not split!
				selProtList.push(checkBoxList.value);
			}
			else { // last/only Idx is protID
				selProtList.push(anaProtVal[anaProtVal.length-1]);
			}
		}
		if (selProtList.length==0) {
			alert('ERROR: No valid proteins selected!');
			return;
		}
		paramStrg+='&protID='+selProtList.join(','); // was ':' before (PP 21/01/15)
	}
	//Creation of the XMLHTTPRequest object
	var XHR = getXMLHTTP();
	if (!XHR) {
		return false;
	}
	XHR.open("POST","$refPromsPath->{cgi}/showProtQuantification.cgi",true);
	//Send the proper header information along with the request
	XHR.setRequestHeader("Content-type", "application/x-www-form-urlencoded; charset=UTF-8");
	XHR.setRequestHeader("Content-length", paramStrg.length);
	XHR.setRequestHeader("Connection", "close");
	// Processing response
	XHR.onreadystatechange=function() {
		if (XHR.readyState==4 && XHR.responseText) {
			targetDiv.innerHTML=XHR.responseText;
			if (job=='saveProt') {
				saveListGlobals.themesFetched=false;
				if ($actionOnCreation) {
					if (themeSel.value=='getLists:-1' \|\| listSel.value.match(':-1')) { // new Theme or List
						$actionOnCreation(themeSel,listSel,XHR);
					}
				}
				// Reset everything to normal after 5 sec
				setTimeout(function(){
								targetDiv.style.display='none';
								if (saveFormBut.style.display=='none') {saveFormBut.style.display='block';} else {saveFormBut.style.visibility='visible';}
							},5000);
			}
			else {saveListGlobals.themesFetched=true;}
		}
	}
	XHR.send(paramStrg);
}
// <--- END OF CHOOSE OR CREATE THEME...
|;
}



#####################################
####<<< Popup layer with info >>>####
#####################################
sub popupInfo {
	my $layerWidth=($_[0])? $_[0] : 200;
	my ($posTypeX,$fromX)=($_[1])? ($_[1],$_[2]) : ('rel',0); # 'rel' or 'abs'
	my ($posTypeY,$fromY)=($_[3])? ($_[3],$_[4]) : ('rel',20);
	print << "END_CODE";
/**********************************************************************************
PopupDescriptions
*   Copyright (C) 2001 Thomas Brattli
*   This script was released at DHTMLCentral.com
*   Visit for more great scripts!
*   This may be used and changed freely as long as this msg is intact!
*   We will also appreciate any links you could give us.
*
*   Made by Thomas Brattli
*
*   Script date: 09/04/2001 (keep this date to check versions)
*********************************************************************************/
function lib_bwcheck(){ //Browsercheck (needed)
	this.ver=navigator.appVersion
	this.agent=navigator.userAgent
	this.dom=document.getElementById?1:0
	this.opera5=(navigator.userAgent.indexOf("Opera")>-1 && document.getElementById)?1:0
	this.ie5=(this.ver.indexOf("MSIE 5")>-1 && this.dom && !this.opera5)?1:0;
	this.ie6=(this.ver.indexOf("MSIE 6")>-1 && this.dom && !this.opera5)?1:0;
	this.ie7=(this.ver.indexOf("MSIE 7")>-1 && this.dom && !this.opera5)?1:0;
	this.ie4=(document.all && !this.dom && !this.opera5)?1:0;
	this.ie=this.ie4||this.ie5||this.ie6||this.ie7
	this.mac=this.agent.indexOf("Mac")>-1
	this.ns6=(this.dom && parseInt(this.ver) >= 5) ?1:0;
	this.ns4=(document.layers && !this.dom)?1:0;
	this.bw=(this.ie7 || this.ie6 || this.ie5 || this.ie4 || this.ns4 || this.ns6 || this.opera5)
	return this
}
var bw=lib_bwcheck()

/***************************************************************************************
Variables to set:
***************************************************************************************/
layerWidth=$layerWidth; // width of popup layer
posTypeX = '$posTypeX'; // fromX is relative or absolute
fromX = $fromX; // How much from the actual mouse X should the description box appear?
posTypeY = '$posTypeY';  // fromY is relative or absolute
fromY = $fromY; // How much from the actual mouse Y should the description box appear?

//To set the font size, font type, border color or remove the border or whatever,
//change the clDescription class in the stylesheet.

//Makes crossbrowser object.
function makeObj(obj){
   	this.evnt=bw.dom? document.getElementById(obj):bw.ie4?document.all[obj]:bw.ns4?document.layers[obj]:0;
	if(!this.evnt) return false
	this.css=bw.dom||bw.ie4?this.evnt.style:bw.ns4?this.evnt:0;
   	this.wref=bw.dom||bw.ie4?this.evnt:bw.ns4?this.css.document:0;
	this.writeIt=b_writeIt;
	return this
}

// A unit of measure that will be added when setting the position of a layer.
var px = bw.ns4||window.opera?"":"px";

function b_writeIt(text){
	if (bw.ns4){this.wref.write(text);this.wref.close()}
	else this.wref.innerHTML = text
}

//Capturing mousemove
var descx = 0
var descy = 0
function popmousemove(e){descx=bw.ns4||bw.ns6?e.pageX:event.x; descy=bw.ns4||bw.ns6?e.pageY:event.y}

function getScrollPosition(){
        return Array((document.documentElement && document.documentElement.scrollLeft) || window.pageXOffset || self.pageXOffset || document.body.scrollLeft,(document.documentElement && document.documentElement.scrollTop) || window.pageYOffset || self.pageYOffset || document.body.scrollTop);
}

var oDesc;
//Shows the messages
function popup(text){
    if(oDesc){
		oDesc.writeIt('<div class="clDescription">'+text+'</div>')

		// Make sure layer stays inside window
		var scrollPosition = getScrollPosition();
		var winSize=(document.body)? document.body.clientWidth : window.innerWidth;
		var posX=(posTypeX=='rel')? descx+fromX : fromX;
		if(navigator.appName.match('Microsoft')) posX = posX + scrollPosition[0];
		oDesc.css.left=(posX-scrollPosition[0] > (winSize-layerWidth))? ((winSize-layerWidth)+scrollPosition[0])+px: (posX)+px;

		//if(bw.ie5||bw.ie6||bw.ie7) descy = descy+document.body.scrollTop
		if(navigator.appName.match('Microsoft')) descy = descy+document.body.scrollTop
		oDesc.css.top = (posTypeY=='rel')? (descy+fromY)+px : fromY+px;
		oDesc.css.visibility = "visible"
//alert(oDesc.css.left+" "+oDesc.css.top);
    }
}
//Hides it
function popout(){
	if(oDesc) oDesc.css.visibility = "hidden"
}
function setPopup(){
   	if(bw.ns4)document.captureEvents(Event.MOUSEMOVE)
    document.onmousemove = popmousemove;
	oDesc = new makeObj('divDescription')
}
END_CODE
}

####>VARIABLE MODIFICATIONS (PTMs)<####
sub getVariableModifications {
	my ($dbh)=@_;

	my $sthModInfo=$dbh->prepare("SELECT ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,DISPLAY_CODE,DISPLAY_COLOR,IS_LABEL,SPECIFICITY FROM MODIFICATION WHERE DISPLAY_CODE IS NOT NULL AND DISPLAY_CODE != '' AND DISPLAY_COLOR IS NOT NULL AND DISPLAY_COLOR != '' ORDER BY PSI_MS_NAME ASC, INTERIM_NAME ASC");
	$sthModInfo->execute;

	my %variableModifications=();
	while (my ($modificationID,$psiMsName,$interimName,$altNames,$dispCode,$dispColor,$isLabel,$specifity) = $sthModInfo->fetchrow_array) {
		if ($altNames) {
			$altNames=~ s/^##//;
			$altNames=(split('##',$altNames))[0]; # keep 1rst
		}
		my $name=$psiMsName || $interimName || $altNames || undef;
		next unless $name;

		@{$variableModifications{$modificationID}}=($name,$dispCode,"mod_$modificationID","color:#$dispColor;",$specifity);
		push @{$variableModifications{$modificationID}},1 if $isLabel;
	}
	$sthModInfo->finish;

	return %variableModifications;
}

##########################################################################################
####<<<Convert a varMod string into 4 values (modCode,Residues,Positions,occurence)>>>####
##########################################################################################
sub convertVarModString {
	my ($varModStrg)=@_;
	my ($residues,$positions)=('','');
	my ($occurence,$varModCode)=(0,'');

	$varModStrg=~s/^_+//; # deletes starting _ if any (Cochin)
	my ($varModText)=($varModStrg=~/^(.+)\s+\([^\(]/); # everything before the last ' (...'
	if (!$varModText) { # For MaxQuant, substitutions does not include the good regex. Ex : 'Asn->Asp' should become 'Asn->Asp (N)'
		if ($varModStrg =~ /->/) {
			my @res=split(/->/,$varModStrg);
			$varModStrg = "$varModStrg ($aaThree2OneLetter{$res[0]})" if $aaThree2OneLetter{$res[0]};
			($varModText)=($varModStrg=~/^(.+)\s+\([^\(]/);
		}
		return ($varModStrg,'','',0) unless $varModText;
	}
	# Change on the 15th of January 2018 because 'AzidoBiotinPEG4 Oxyde (K)' was recognized by this test
	#if ($varModText=~/(\d+) (.+)/) {
	if ($varModText=~/^(\d+) (.+)/) { # PMF with multiple occurences
		$occurence=$1;
		$varModCode=$2;
	}
	else {$varModCode=$varModText;}
	my ($modPosStrg)=($varModStrg=~/\(([^)]+)\)\Z/); # everything inside the last '(...)'
	if ($modPosStrg=~/[CN]\s*-\s*term(\s*\w*)/i) {
		my $aa=($1 && $1 ne 'inus')?"$1;":''; # be careful with terminus matching!!!
		if ($modPosStrg=~/N\s*-\s*term/i) {
			my $spec=($modPosStrg =~ /Protein/)? '-' : '=';
			$residues="$aa$spec";
			$positions=$spec; #1
		}
		else { # C-term
			my $spec=($modPosStrg =~ /Protein/)? '+' : '*';
			$residues="$aa$spec";
			$positions=$spec; #0
		}
		$occurence=1;
		$positions=-1 if $modPosStrg=~/:-1/; # Uncertainty for -term position (PP 25/01/13)
	}
	elsif ($modPosStrg!~/[a-z]/) { # !Unknown inside-text : ex: N-Acetyl (Protein)
		#($residues,$positions)=split(/:/,$modPosStrg); # MIS
		if ($modPosStrg=~/(.+):(.+)/) { # MIS
			$residues=$1;
			$positions=$2;
			$occurence=-1; # can be computed based on positions
		}
		else { # PMF
			$residues=$modPosStrg;
			$occurence=1 unless $occurence; # already defined if > 1
			$positions='' unless $positions;
		}
		if (length($residues)>1) {
			$residues=join(',',split(//,$residues));
		}
	}
	return ($varModCode,$residues,$positions,$occurence); # residues='STY'(= for N-term, * for C-term, - for Protein N-term, + for Protein C-term) // positions='2.5.9' (0 for N-term, 999 for C-term)

}

sub convertVarModStringSwath {
	my ($varModStrg)=@_;
	my @result;
	my ($residues,$positions)=('','');
	my ($occurence,$varModCode)=(0,'');
	my @pepMod=split("/",$varModStrg);
	foreach my $modif (@pepMod){
		next if (length($modif)==1 || length($modif)==2);
		my @mod=split(/,/,$modif);
		$residues=$mod[1];
		$varModCode=$mod[2];
		if ($mod[0]==-1) {$positions="=";} 		#(= for N-term)
		else{$positions=$mod[0]+1;}
		push(@result,"$varModCode!$residues!$positions");
	}
	return (@result);
}


##############################################################
####<<<String representation of a variable modification>>>####
##############################################################
sub toStringVariableModifications { # Too slow. Avoid if large dataset. Modify peptide query and use &decodeVarMod instead !!!!!!!!!!!!!!!!!!!!!!!!
	my ($dbh,$item,$itemID,$analysisID,$pepSeq,$extraParam,$posShift)=@_;
	$posShift=0 unless $posShift; # shift aa pos <=> pep beg-1 (from pos in peptide to pos in protein)
	my $extraQuery=($item eq 'QUERY' && $extraParam)? "AND PEP_RANK=$extraParam" : "";
	my $sthVmods=($item eq 'SPECTRUM')? $dbh->prepare("SELECT M.ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,SM.SPECIFICITY,POS_STRING FROM MODIFICATION M,SPECTRUM_MODIFICATION SM WHERE SM.ID_MODIFICATION=M.ID_MODIFICATION AND ID_SPECTRUM=$itemID AND SM.MODIF_TYPE='V' ORDER BY PSI_MS_NAME ASC, INTERIM_NAME ASC")
					: $dbh->prepare("SELECT M.ID_MODIFICATION,PSI_MS_NAME,INTERIM_NAME,SYNONYMES,AM.SPECIFICITY,POS_STRING,REF_POS_STRING FROM MODIFICATION M,${item}_MODIFICATION IM, ANALYSIS_MODIFICATION AM WHERE AM.ID_MODIFICATION=M.ID_MODIFICATION AND M.ID_MODIFICATION=IM.ID_MODIFICATION AND ID_$item=$itemID AND AM.ID_ANALYSIS=$analysisID $extraQuery ORDER BY PSI_MS_NAME ASC, INTERIM_NAME ASC"); # Only variable mods are stored in this table

	my @sequence=split(//,$pepSeq);
	my (%varMods,%excludedVarMods);
	my $refExcludedVM=($item eq 'PEPTIDE' && $extraParam)? $extraParam : \%excludedVarMods; # reference to a list of modIDs to skip (eg. labeled modifs)
	$sthVmods->execute;
	while (my ($modID,$psiMsName,$interimName,$altNames,$specificity,$posString,$refPosString) = $sthVmods->fetchrow_array) {
		next if $refExcludedVM->{$modID}; # skip this mod
		$altNames=~ s/##/,/g if $altNames;
		$altNames=~ s/^,//; $altNames=~ s/,\Z//; # remove starting & trailing "," if any
		my $name=($psiMsName)? $psiMsName : ($interimName)? $interimName : ($altNames)? $altNames : '';
		my %aa=();
		#my $oldPosTag=-1;
		my $vmodInString=scalar keys %varMods;
		foreach my $pos (split(/\./,$posString)) {
			#$oldPosTag++;
			#if ($pos eq '-1') { # eq because possible = - * +
			#	my @refPositions=split(/\./,$refPosString);
			#	$aa{$sequence[$refPositions[$oldPosTag]-1]}=1;
			#	next;
			#}
			if ($pos =~ /-\d+/) { # eq because possible = - * +
				foreach my $pos2 (split(/;/,$pos)) {
					if ($pos2 eq '--' || $pos2 eq '-=') {
						$aa{$sequence[0]}=1;
						next;
					}elsif ($pos2 eq '-+' || $pos2 eq '-*'){
						$aa{$sequence[$#sequence]}=1;
						next;
					}
					$aa{$sequence[abs($pos2)-1]}=1;
				}
				next;
			}
			# To distinguish Acetyl(K) and Acetyl(Protein N-term)
			if ($specificity && $pos =~ /=|-|\+|\*/ ) {
				my $toAdd='';
				if ($specificity =~ /=/ && $pos =~ /=/) {
					$toAdd='N-term';
				}
				elsif ($specificity =~ /-/ && $pos =~ /-/) {
					$toAdd='Protein N-term';
				}
				elsif ($specificity =~ /\*/ && $pos =~ /\*/) {
					$toAdd='C-term';
				}
				elsif ($specificity =~ /\+/ && $pos =~ /\+/) {
					$toAdd='Protein C-term';
				}
				$toAdd='?' if $pos=~ /-(\d+|=|-|\+|\*)/;
				if ($toAdd) {
					$pos=quotemeta($pos);
					#$posString=~s/\.*$pos\.*//;
					###> Modified the 29/10/13 because if $posString='6.-.10', then, it becomes '610' which is wrong!
					if ($posString=~/\.+$pos\.+/) {
						$posString=~s/\.*$pos\.*/\./;
					}else {
						$posString=~s/\.*$pos\.*//;
					}
					$varMods{"$name ($toAdd)"}=1;
					next;
				}
			}
			$aa{$sequence[$pos-1]}=1 unless $pos=~/\D/; # make sure it's a number, otherwise -> ?
		}
		#$posString=~s/-1/?/g if $posString; # -1 -> ? (position ambiguity)
		if ($posString) { # -1;-2;-3 means ambiguity in position 1, 2 and 3 and is replaced by ?
			$posString=~s/(-(\d+|=|-|\+|\*);{0,1})+/?/g;
			$posString=~s/\.\./\.\?\./g; # Add ambigous to match expected modification sites amount
			$posString=~s/\.$/\.\?/;
			$posString=~s/(\d+)/$1+$posShift/eg if $posShift; # e flag allows math inside regexp
		}else{
			$posString='?';
		}
		my $aaString=join('',sort{$a cmp $b} keys %aa);
		$aaString='?' unless $aaString;
		#if ($aaString) {
		next if $aaString eq '?' && $vmodInString - scalar keys %varMods != 0; # without this test, '+ Acetyl (N-term)' becomes '+ Acetyl (N-term) + Acetyl(?:?)'
		$varMods{"$name ($aaString:$posString)"}=1;
		#}
	}
	$sthVmods->finish;

	my $vmodString = (scalar keys %varMods)? ' + '.join(' + ',sort{lc($a) cmp lc($b)} keys %varMods) : ''; # prevents duplicates
	return $vmodString;
}

sub decodeVarMod {
	my ($dbh,$pepSeq,$modCode,$refModName,$posShift)=@_;
	# $modCode= '(&)<modID1>:<pos1.pos2.-pos3;-pos4;-pos5...>&<modID2>:<pos1.pos2...>...'
	# for modID1, there is an ambiguity on the third position between 3 potential sites
	# $refModName= reference of array {modID=>modName} empty at 1 call, updated as new modif are seen
	# $posShift= delta applied on modif pos (typically to display pos in protein instead of peptide)
	$modCode='' unless $modCode;
	$posShift=0 unless $posShift;
	my %modName;
	unless ($refModName) {
		$refModName=\%modName;
	}
	my @sequence=split(//,$pepSeq);
	my %varMods;
	my $sthVM=$dbh->prepare('SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES FROM MODIFICATION WHERE ID_MODIFICATION=?');
	$modCode=~s/^&//; # just in case
	foreach my $vCode (split(/&/,$modCode)) {
		my ($modID,$posString)=split(/:/,$vCode);
		unless ($refModName->{$modID}) {
			$sthVM->execute($modID);
			my ($psiMsName,$interimName,$altNames) = $sthVM->fetchrow_array;
			if ($altNames) {
				$altNames=~ s/##/,/g;
				$altNames=~ s/^,//; $altNames=~ s/,\Z//; # remove starting & trailing "," if any
			}
			$refModName->{$modID}=($psiMsName)? $psiMsName : ($interimName)? $interimName : ($altNames)? $altNames : '';
		}
		my %aa;
		my $vmodInString = scalar keys %varMods;
		my $modNtCt = '';
		my @newPos = ();
		my @ambigousSites = ();
		foreach my $pos (split(/\./,$posString)) {
			if ($pos =~ /;|-\d+/) { # in case of ambiguity of a position
				foreach my $posNeg (split(/;/,$pos)) {
					if ($posNeg eq '--' || $posNeg eq '-=') {
						$aa{$sequence[0]}=1;
						next;
					}
					elsif ($posNeg eq '-+' || $posNeg eq '-*') {
						$aa{$sequence[$#sequence]}=1;
						next;
					}
					$aa{$sequence[abs($posNeg)-1]}=1;
					push @ambigousSites, (-1*$posNeg)+$posShift;
				}
				push @newPos,'?';
			}
			elsif ($pos=~/\D/) { # possible = - * +
				$modNtCt=($pos eq '-')? 'Protein N-term' : ($pos eq '=')? 'N-term' : ($pos eq '+')? 'Protein C-term' : 'C-term'; # '*'
			}
			#else {$aa{$sequence[$pos-1]}=1; push @newPos,$pos+$posShift;}
			else {
				if ($sequence[$pos-1]) {
					$aa{$sequence[$pos-1]}=1;
				}
				else { # ERROR: pos is outside sequence's length!!!! Should never happen unless bug during data import
					#print "*$pepSeq*$modCode*$vCode*$posString*$pos*<BR>\n";
					$aa{'!'}=1;
				}
				push @newPos,$pos+$posShift;
			}
		}
		
		# Add ambigous to match expected modification sites amount
		my $remainingSites = () = $posString =~ /\.\./g; 
		$remainingSites += () = $posString =~ /\.$/g;
		for(my $i=0; $i<$remainingSites; $i++) {
			push @newPos,'?';
		}
		
		if(@ambigousSites) {
			$newPos[-1] .= "[".join(',', @ambigousSites)."]";
		}
		
		$varMods{"$refModName->{$modID} ($modNtCt)"}=1 if $modNtCt;
		$varMods{"$refModName->{$modID} (".join('',sort{$a cmp $b} keys %aa).':'.join('.',@newPos).")"}=1 if scalar keys %aa;
	}
	$sthVM->finish;
	return join(' + ',sort{lc($a) cmp lc($b)} keys %varMods);
}


#######################################################
####<<< Get the modification ID from myProMS DB >>>####
#######################################################
sub getModificationIDfromString { # requires XML::SAX::ParserFactory
	my ($dbh,$vmodString,$aaOneLetter,$vmods)=@_; # $aaOneLetter & $vmods are optional
	$vmodString =~ s/^\s+//; # remove first spaces if any ex: '   Acetyl' becomes 'Acetyl'
	$vmodString =~ s/\s+\Z//; # remove trailing spaces if any
	unless ($vmodString) {
		$dbh->rollback;
		$dbh->disconnect;
		die "Error in &promsMod::getModificationIDfromString: \$vmodString is empty or undef!\n";
	}
	$vmods={} unless $vmods; # anomymous hash reference. just to be safe
	if ($aaOneLetter) { # clean specif string
		$aaOneLetter=~s/[^\w\-\=\+\*]//g;
		$aaOneLetter=~s/[\d_]//g;
	}

	###> If this string was already referenced in myProMS MODIFICATION Table, get its number
	my $vmodStringUC=uc($vmodString);

	###> Known modifications problems !
	if ($vmodStringUC eq 'PROTEIN TERMINAL ACETYL') {
		$vmodStringUC='ACETYL';
	}
	elsif ($vmodStringUC eq 'GLNTHRGLYGLY') { # For this sumoylation, no other choice than changing the name of the VMOD (no synonyms in unimod_tables.xml version 05/04/2013)
		$vmodStringUC='QTGG';
	}
	$vmodStringUC=~s/-ADD//g; # Arg-add in Paragon becomes Arg in Unimod

	###> Look into myProMS database
	my ($modificationID,$specificityString)=$dbh->selectrow_array("SELECT ID_MODIFICATION,SPECIFICITY FROM MODIFICATION WHERE UPPER(PSI_MS_NAME)='$vmodStringUC' OR UPPER(INTERIM_NAME)='$vmodStringUC' LIMIT 1");
	unless ($modificationID) {
		($modificationID,$specificityString)=$dbh->selectrow_array("SELECT ID_MODIFICATION,SPECIFICITY FROM MODIFICATION WHERE UPPER(SYNONYMES) LIKE '%##$vmodStringUC##%' LIMIT 1");
	}
	$specificityString='' unless $specificityString;

	if ($modificationID) { # Modification was found in myProMS, check now the specificity.
		###> Has to check specificity
		if ($aaOneLetter) {
			foreach my $aa (split(//,$aaOneLetter)) { # prevents erroneous multi-letter strg to be inserted (eg STY instead of S,T,Y) PP 30/01/17
				my $qAA=quotemeta($aa); # * +
				if ($specificityString !~ /$qAA/) { # && $specificityString !~ /=|-|\+|\*/ !!!! Why "$specificityString !~ /=|-|\+|\*/"????? (PP 17/08/16)
					$specificityString.=($specificityString)? ",$aa" : $aa;
					#$dbh->do("UPDATE MODIFICATION SET SPECIFICITY='$specificityString,$qAA' WHERE ID_MODIFICATION=$modificationID");
					my $sthUpM=$dbh->prepare("UPDATE MODIFICATION SET SPECIFICITY=? WHERE ID_MODIFICATION=$modificationID");
					$sthUpM->execute($specificityString);
					$sthUpM->finish;
					$dbh->commit;
				}
			}
		}
	}
	###> Look into unimod_tables.xml
	else {
		my %promsPath = &promsConfig::getServerInfo('no_user');
		my $xmlFile="$promsPath{data}/unimod_tables.xml"; # unimod_tables.xml file
		if (!$vmods->{1}) { # if unimod_tables.xml was already read before, do not need to read it again!
			my $handler = UnimodXMLHandler->new($xmlFile,$vmods);
			require XML::SAX::ParserFactory; # to be safe
			my $xmlparser = XML::SAX::ParserFactory->parser(Handler => $handler );
			$xmlparser->parse_uri($xmlFile);
		}

		my $sthChkVMod=$dbh->prepare("SELECT ID_MODIFICATION FROM MODIFICATION WHERE UNIMOD_ACC=? LIMIT 1"); # LIMIT just to be safe (should match only 1 modif)
		my $sthAddVMod=$dbh->prepare("INSERT INTO MODIFICATION (ID_MODIFICATION,UNIMOD_ACC,PSI_MS_NAME,INTERIM_NAME,DES,SYNONYMES,COMPOSITION,MONO_MASS,AVGE_MASS,SPECIFICITY,IS_LABEL,VALID_STATUS) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)");
		my ($maxModID)=$dbh->selectrow_array("SELECT MAX(ID_MODIFICATION) FROM MODIFICATION");
		foreach my $unimodID (keys %{$vmods}) {
			my $shortNameUC=uc($vmods->{$unimodID}{'PSI_MS_NAME'});
			my $interimNameUC=uc($vmods->{$unimodID}{'INTERIM_NAME'});
			my $synonymesUC=($vmods->{$unimodID}{'SYNONYMES'})?uc($vmods->{$unimodID}{'SYNONYMES'}):undef;
			my $synonymMatch=0;
			if ($synonymesUC) {
				foreach my $syn (split(/##/,$synonymesUC)) {
					next unless $syn;
					$synonymMatch=1 if ($syn eq $vmodStringUC);
				}
			}

			if ( $vmodStringUC eq $shortNameUC || $vmodStringUC eq $interimNameUC || $synonymMatch ) { # Modification GlnGlnGlnThrGlyGly in Paragon is found as a synonym in Unimod
				($specificityString)=($vmods->{$unimodID}{'SPECIFICITY'})? join(',',keys %{$vmods->{$unimodID}{'SPECIFICITY'}}) : '';

				if ($aaOneLetter) {
					foreach my $aa (split(//,$aaOneLetter)) { # prevents erroneous multi-letter strg to be inserted (eg STY instead of S,T,Y) PP 30/01/17
						my $qAA=quotemeta($aa); # *+
						if ($specificityString !~ /$qAA/) { # && $specificityString !~ /=|-|\+|\*/
							$specificityString.=($specificityString)? ",$aa" : $aa;
						}
					}
				}
				$sthChkVMod->execute($unimodID);
				($modificationID)=$sthChkVMod->fetchrow_array;
				if (!$modificationID) { # Security : Modification not inserted before
					my $psiName=$vmods->{$unimodID}{'PSI_MS_NAME'}; # Acetyl
					my $interimName=$vmods->{$unimodID}{'INTERIM_NAME'}; # Acetyl
					my $des=$vmods->{$unimodID}{'DES'}; # Acetylation
					my $composition=$vmods->{$unimodID}{'COMPOSITION'}; # H(2) C(2) O
					my $monoMass=$vmods->{$unimodID}{'MONO_MASS'}; # 42.010565
					my $avgeMass=$vmods->{$unimodID}{'AVGE_MASS'}; # 42.010565
					my $synonymes=($vmods->{$unimodID}{'SYNONYMES'}) ? $vmods->{$unimodID}{'SYNONYMES'} : '';
					#my $isLabel=($psiName=~/label|(ITRAQ|TMT)\d+\s*plex/i || $interimName=~/label|(ITRAQ|TMT)\d+\s*plex/i || $des=~/label|(ITRAQ|TMT)\d+\s*plex/i)? 1 : 0;
					my $isLabel=$vmods->{$unimodID}{'IS_LABEL'} // 0;
					$modificationID=++$maxModID;
					$sthAddVMod->execute($modificationID,$unimodID,$psiName,$interimName,$des,$synonymes,$composition,$monoMass,$avgeMass,$specificityString,$isLabel,1);
					$dbh->commit;
				}
				#else { # Security : Modification inserted before -> code that would never be executed in theory
				#	($modificationID)=$dbh->selectrow_array("SELECT ID_MODIFICATION FROM MODIFICATION WHERE UNIMOD_ACC=$unimodID");
				#}
			}
		}
		$sthChkVMod->finish;

		###> The modification was not found in Unimod_tables.xml ; a modification with a 0 VALID_STATUS is inserted in myProMS
		if (!$modificationID) {
			my $isLabel=($vmodString=~/label/i)? 1 : 0;
			$modificationID=++$maxModID;
			$sthAddVMod->execute($modificationID,undef,undef,undef,undef,"##$vmodString##",undef,undef,undef,undef,$isLabel,0);
			$dbh->commit;
		}
		$sthAddVMod->finish;
	}

	return ($modificationID);
}

######################################
#######<Generates Match Groups>#######
######################################
sub updateMatchGroups {
	# builds match groups from scratch for a single analysis & updates table ANALYSIS_PROTEIN
	# entries must already be in table ANALYSIS_PROTEIN except for MATCH_GROUP & VISIBILITY that will be updated
	# parameters:
	#		-$dbh,
	#		-id of selected analysis,
	#		-empty hash ref for storing project-wide data to speed up next calls to &updateMatchGroups (optional, especially if only 1 call is performed)
	#		-hash ref of options:
	#				-VERBOSE: 0/1 (default 0), print progress status
	#				-COMMIT: 0/1 (default 0), commit after update
	#				-USE_GHOST: 0/1 (default 0), use ghost peptides to compute peptide/protein matches
	my ($dbh,$analysisID,$refProjectData,$refOptions)=@_;
	$refProjectData={BEST_VIS=>(),NUM_TOP=>(),VIS=>0} unless $refProjectData;
	$refProjectData->{BEST_VIS}=() unless $refProjectData->{BEST_VIS};
	$refProjectData->{NUM_TOP}=() unless $refProjectData->{NUM_TOP};
	my $verbose=$refOptions->{VERBOSE} || 0;
	my $commit=$refOptions->{COMMIT} || 0;
	my $useGhostPep=$refOptions->{USE_GHOST} || 0;

	print "<FONT class=\"title3\">&nbsp;-Fetching data to compute match groups..." if $verbose;
	my $projectID=&getProjectID($dbh,$analysisID,'analysis');
	unless (defined $refProjectData->{VIS}) {
		($refProjectData->{VIS})=$dbh->selectrow_array("SELECT PROT_VISIBILITY FROM PROJECT WHERE ID_PROJECT=$projectID");
		$refProjectData->{VIS}=0 unless $refProjectData->{VIS};
	}
	print '.' if $verbose;

	####<Fetching project proteins & number of times proteins are at top of match group
	unless (scalar keys %{$refProjectData->{BEST_VIS}}) {
		my $sthProjProt=$dbh->prepare("SELECT P.ID_PROTEIN,MAX(VISIBILITY) FROM ANALYSIS_PROTEIN AP,PROTEIN P WHERE AP.ID_PROTEIN=P.ID_PROTEIN AND ID_ANALYSIS != $analysisID AND ID_PROJECT=$projectID GROUP BY P.ID_PROTEIN");
		$sthProjProt->execute;
		while (my ($protID,$bestVis)=$sthProjProt->fetchrow_array) {
			$refProjectData->{BEST_VIS}{$protID}=$bestVis || 0;
		}
		$sthProjProt->finish;
		print '.' if $verbose;
		my $sthTop=$dbh->prepare("SELECT P.ID_PROTEIN,COUNT(*) FROM ANALYSIS_PROTEIN AP,PROTEIN P WHERE AP.ID_PROTEIN=P.ID_PROTEIN AND ID_ANALYSIS != $analysisID AND ID_PROJECT=$projectID AND VISIBILITY=2 GROUP BY P.ID_PROTEIN");
		$sthTop->execute;
		while (my ($protID,$numTop)=$sthTop->fetchrow_array) {
			$refProjectData->{NUM_TOP}{$protID}=$numTop;
		}
		$sthTop->finish;
		print '.' if $verbose;
	}

	####<Fetching Analysis protein data
	my ($preferredSpecies)=$dbh->selectrow_array("SELECT SCIENTIFIC_NAME FROM ANALYSIS A,SAMPLE S,EXPERIMENT E,SPECIES SP WHERE A.ID_SAMPLE=S.ID_SAMPLE AND S.ID_EXPERIMENT=E.ID_EXPERIMENT AND E.ID_SPECIES=SP.ID_SPECIES AND ID_ANALYSIS=$analysisID");

	my %matchList;
	my $ghostStrg=($useGhostPep)? '' : 'AND (VALID_STATUS > 0 OR PEP_BEG > 0)';
	my $sthProtPep=$dbh->prepare("SELECT ID_PROTEIN,PEP_SEQ FROM PEPTIDE P,PEPTIDE_PROTEIN_ATTRIB PPA WHERE P.ID_PEPTIDE=PPA.ID_PEPTIDE AND P.ID_ANALYSIS=$analysisID $ghostStrg");
	$sthProtPep->execute;
	while (my ($protID,$pepSeq)=$sthProtPep->fetchrow_array) {
		$matchList{$protID}{$pepSeq}=1;
	}
	$sthProtPep->finish;
	print '.' if $verbose;

	my (%numProtPeptides,%numProtTop,%bestProtVis,%proteinPepDens,%protSpeciesClass,%proteinScore,%proteinLength,%proteinAnnotQual);
	my $sthAnaProt=$dbh->prepare("SELECT AP.ID_PROTEIN,IDENTIFIER,SCORE,PROT_LENGTH,PROT_DES,ORGANISM FROM ANALYSIS_PROTEIN AP,PROTEIN P WHERE AP.ID_PROTEIN=P.ID_PROTEIN AND AP.ID_ANALYSIS=$analysisID");
	$sthAnaProt->execute;
	while (my ($protID,$identifier,$score,$length,$des,$organism)=$sthAnaProt->fetchrow_array) {
		$organism='Unknown' unless $organism;
		$numProtTop{$protID}=$refProjectData->{NUM_TOP}{$protID} || 0; # 0 if not already in project
		$numProtPeptides{$protID}=scalar (keys %{$matchList{$protID}});
		$proteinScore{$protID}=$score || 0;
		$proteinLength{$protID}=$length || 0;
		$proteinPepDens{$protID}=$proteinLength{$protID}/$numProtPeptides{$protID};
		$protSpeciesClass{$protID}=($preferredSpecies && $organism eq $preferredSpecies)? 1 : 0;
		$proteinAnnotQual{$protID}=($identifier=~/\b\w{1,5}_/)? 1 : ($identifier=~/\b\w{6}_/)? 2 : ($des=~/no\sdescription/i)? 6 : ($des=~/unnamed/i || $des=~/unknown/i)? 5 : ($des=~/hypothetical/i)? 4 : 3;  #PP
	}
	$sthAnaProt->finish;
	print '.' if $verbose;

	my @sortedProtIDs=sort{$numProtPeptides{$b}<=>$numProtPeptides{$a} || $protSpeciesClass{$b}<=>$protSpeciesClass{$a} || $numProtTop{$b}<=>$numProtTop{$a} || $proteinScore{$b}<=>$proteinScore{$a} || &deltaLength($proteinLength{$a},$proteinLength{$b},$proteinPepDens{$a},$proteinPepDens{$b}) || $proteinAnnotQual{$a}<=>$proteinAnnotQual{$b} || $proteinLength{$a}<=>$proteinLength{$b} || $a<=>$b} keys %matchList;
	print " Done</FONT><BR>\n" if $verbose;

	####<Create match groups
	my (%matchGroups,%visibility);
	&createMatchGroups(\%matchList,\%matchGroups,\%visibility,$refProjectData->{BEST_VIS},$refProjectData->{VIS},\@sortedProtIDs,$verbose);

	####<Update ANALYSIS_PROTEIN & project global data
	print "<FONT class=\"title3\">&nbsp;-Updating database..." if $verbose;
	my $sthUpAP=$dbh->prepare("UPDATE ANALYSIS_PROTEIN SET MATCH_GROUP=?,VISIBILITY=? WHERE ID_ANALYSIS=$analysisID AND ID_PROTEIN=?");
	my $count=0;
	foreach my $protID (keys %matchList) {
		$sthUpAP->execute($matchGroups{$protID},$visibility{$protID},$protID);
		$refProjectData->{BEST_VIS}{$protID}=$visibility{$protID} if (!defined $refProjectData->{BEST_VIS}{$protID} || $refProjectData->{BEST_VIS}{$protID} < $visibility{$protID});
		$refProjectData->{NUM_TOP}{$protID}++ if $visibility{$protID}==2;
		if ($verbose) {
			print '.' unless ++$count % 500;
		}
	}
	$dbh->commit if $commit;
	print " Done</FONT><BR>\n" if $verbose;

	####<Internal sub>####
	sub deltaLength {
		my ($l_a,$l_b,$d_a,$d_b)=@_;
		my $pepDensVal=($l_b > $l_a)? $d_a : $d_b;
		if (abs($l_b-$l_a) > $pepDensVal) {return $l_a<=>$l_b} else {return 0}
	}
}
sub createMatchGroups { # computes match groups from provided data & parameters (can also be called directly from a script)
	my ($refmatchList,$refmatchGroup,$refVisibility,$refBestProtVis,$protVisibility,$refsortedIdentifiers,$verbose)=@_;
	print "<FONT class=\"title3\">&nbsp;-Grouping proteins: 0%..." if $verbose;
	my $numGroup=0;
    %{$refmatchGroup}=(); # just to be safe
	%{$refVisibility}=(); # just to be safe
	my $count=0;
	my $numProt=scalar @{$refsortedIdentifiers};
	my $currPC=0;
	foreach my $i (0..$#{$refsortedIdentifiers}) {
		#if ($verbose) {
		#	$count++;
		#	if ($count==150 && $verbose) {print '.'; $count=0;}
		#}
		next if $refmatchGroup->{$refsortedIdentifiers->[$i]}; # already assigned to a match group
		$refmatchGroup->{$refsortedIdentifiers->[$i]}=++$numGroup;
		$refVisibility->{$refsortedIdentifiers->[$i]}=2; # Reference protein
		$refBestProtVis->{$refsortedIdentifiers->[$i]}=2; # update bestVis
#next; # SKIP grouping!!!
		foreach my $j ($i+1..$#{$refsortedIdentifiers}) {
			next if $refmatchGroup->{$refsortedIdentifiers->[$j]}; # already assigned to a match group
			if ($verbose) {
				$count++;
				if ($count >= 500000) {print '.'; $count=0;}
			}
			##<Comparing peptide contents of identifier#1 and identifier#2>## All peptides must match!
			my $matchOK=1;
			foreach my $seq (keys %{$refmatchList->{$refsortedIdentifiers->[$j]}}) {
				if (!$refmatchList->{$refsortedIdentifiers->[$i]}{$seq}) {
					delete $refmatchList->{$refsortedIdentifiers->[$i]}{$seq}; # to be safe
					$matchOK=0;
					last;
				}
			}
			if ($matchOK) {
				$refmatchGroup->{$refsortedIdentifiers->[$j]}=$refmatchGroup->{$refsortedIdentifiers->[$i]}; # $numGroup
				$refVisibility->{$refsortedIdentifiers->[$j]}=($protVisibility && defined($refBestProtVis->{$refsortedIdentifiers->[$j]}) && ($refBestProtVis->{$refsortedIdentifiers->[$j]}==2 || ($protVisibility==2 && $refBestProtVis->{$refsortedIdentifiers->[$j]})))? 1 : 0; # visible or hidden
				$refBestProtVis->{$refsortedIdentifiers->[$j]}=$refVisibility->{$refsortedIdentifiers->[$j]} if (!defined($refBestProtVis->{$refsortedIdentifiers->[$j]}) || $refVisibility->{$refsortedIdentifiers->[$j]} > $refBestProtVis->{$refsortedIdentifiers->[$j]}); # update bestVis
			}
		}
		if ($verbose) {
			my $newPC=10*int(10*(scalar keys %{$refmatchGroup})/$numProt);
			if ($newPC > $currPC) {
				for (my $pc=$currPC+10; $pc<=$newPC; $pc+=10) {
					print "$pc%";
					print "..." if $pc < 100;
				}
				$currPC=$newPC;
			}
		}
	}
	print " Done.</FONT><BR>\n" if $verbose;
	return $numGroup;
}


###############################################################################
###########################    UNIMOD PARSING   ###############################
###############################################################################
package UnimodXMLHandler; {

	my (%positions);

	sub new {#Constructor of the mzXMLHandler
		my ($type,$xmlFileName,$vmods)= @_;
		my $self = bless ({}, $type);
		$self->{'xmlFileName'}=$xmlFileName;
		$self->{'vmods'}=$vmods;
		return $self;
	}

	sub start_document {
		my ($self) = @_;
		%positions=();
		#print "Starting reading XML file to get VMODS !<BR>\n";
	}

	sub end_document {
		my ($self) = @_;
		#print "Finishing reading XML file<BR>\n";
	}

	sub start_element {# Read the elements of the XML
		my ($self, $element) = @_;

		if($element->{'Name'} eq "modifications"){
			$self->{'onmutation'}=1;
		}elsif ($element->{'Name'} eq "neutral_losses") {
			$self->{'onnl'}=1;
		}elsif ($element->{'Name'} eq "spec2nl") {
			$self->{'onspec2nl'}=1;
		}elsif ($element->{'Name'} eq "specificity") {
			$self->{'onspec'}=1;
		}elsif ($element->{'Name'} eq "alt_names") {
			$self->{'onaltNames'}=1;
		}elsif ($element->{'Name'} eq "positions") {
			$self->{'onposition'}=1;
		}
elsif ($element->{'Name'} eq "classifications") {
			$self->{'onclassification'}=1;
		}
#<List of modification types
if ($self->{'onclassification'} && $element->{'Name'} eq "classifications_row"){
	$self->{'label_record_id'}=$element->{'Attributes'}{'{}record_id'}{'Value'} if $element->{'Attributes'}{'{}classification'}{'Value'} eq 'Isotopic label'; # id=11
}
		#<modifications_row username_of_poster="JMcGouran" avge_mass="100.1191" ex_code_name="" mono_mass="100.063663" record_id="1257" full_name="Ub Bromide probe addition" code_name="Ub-Br2" approved="0" date_time_modified="2011-07-11 11:30:04" date_time_posted="2011-07-01 13:25:05" composition="H(8) C(4) N(2) O" group_of_poster="">
		elsif ($self->{'onmutation'} && $element->{'Name'} eq "modifications_row"){
			my $unimodID=$element->{'Attributes'}{'{}record_id'}{'Value'}; # 1
			$self->{'vmods'}->{$unimodID}{'PSI_MS_NAME'}=($element->{'Attributes'}{'{}ex_code_name'}{'Value'});
			$self->{'vmods'}->{$unimodID}{'PSI_MS_NAME'}=~s/&gt;/>/g;# In unimod_tables.xml, the modification  Pro->pyro-Glu becomes Pro-&gt;pyro-Glu
			$self->{'vmods'}->{$unimodID}{'INTERIM_NAME'}=$element->{'Attributes'}{'{}code_name'}{'Value'}; # Acetyl
			$self->{'vmods'}->{$unimodID}{'INTERIM_NAME'}=~s/&gt;/>/g;
			$self->{'vmods'}->{$unimodID}{'DES'}=$element->{'Attributes'}{'{}full_name'}{'Value'}; # Acetylation
			$self->{'vmods'}->{$unimodID}{'COMPOSITION'}=$element->{'Attributes'}{'{}composition'}{'Value'}; # H(2) C(2) O
			$self->{'vmods'}->{$unimodID}{'MONO_MASS'}=$element->{'Attributes'}{'{}mono_mass'}{'Value'}; # 42.010565
			$self->{'vmods'}->{$unimodID}{'AVGE_MASS'}=$element->{'Attributes'}{'{}avge_mass'}{'Value'}; # 42.010565
		}
		elsif ($self->{'onspec'} && $element->{'Name'} eq "specificity_row"){
			my $oneLetter=$element->{'Attributes'}{'{}one_letter'}{'Value'}; # H Y K ... or N-term / C-term
			if (length($oneLetter)>1) { # Code for terminal modifications.
				my $posString=$positions{$element->{'Attributes'}{'{}position_key'}{'Value'}};
				$oneLetter='=' if $posString eq 'Any N-term';
				$oneLetter='-' if $posString eq 'Protein N-term';
				$oneLetter='*' if $posString eq 'Any C-term';
				$oneLetter='+' if $posString eq 'Protein C-term';
			}
			else { # Gln->pyro-Glu on Q Any N-term -> it is important to keep this information
				my $posString=$positions{$element->{'Attributes'}{'{}position_key'}{'Value'}};
				$oneLetter.=';=' if $posString eq 'Any N-term';
				$oneLetter.=';-' if $posString eq 'Protein N-term';
				$oneLetter.=';*' if $posString eq 'Any C-term';
				$oneLetter.=';+' if $posString eq 'Protein C-term';
			}
			my $unimodID=$element->{'Attributes'}{'{}mod_key'}{'Value'};
			$self->{'vmods'}->{$unimodID}{'SPECIFICITY'}{$oneLetter}=1;
$self->{'vmods'}->{$unimodID}{'IS_LABEL'}=(!$self->{'vmods'}->{$unimodID}{'IS_LABEL'} && $element->{'Attributes'}{'{}classifications_key'}{'Value'}==$self->{'label_record_id'})? 1 : 0;
		}
		elsif ($self->{'onaltNames'} && $element->{'Name'} eq "alt_names_row") {
			my $unimodID=$element->{'Attributes'}{'{}mod_key'}{'Value'};
			if (!defined($self->{'vmods'}->{$unimodID}{'SYNONYMES'})) {
				$self->{'vmods'}->{$unimodID}{'SYNONYMES'}="##$element->{'Attributes'}{'{}alt_name'}{'Value'}##";
			}
			else{
				$self->{'vmods'}->{$unimodID}{'SYNONYMES'}.="$element->{'Attributes'}{'{}alt_name'}{'Value'}##";
			}
		}
		elsif ($self->{'onposition'} && $element->{'Name'} eq "positions_row") {
			my $posID=$element->{'Attributes'}{'{}record_id'}{'Value'};
			my $posString=$element->{'Attributes'}{'{}position'}{'Value'}; # Can take as value: '-', 'Anywhere', 'Any N-term', 'Any C-term', 'Protein N-term' & 'Protein C-term'
			$positions{$posID}=$posString; #
		}
	}

	sub end_element {# Read the elements of the XML
		my ($self, $element) = @_;

		if ($element->{'Name'} eq "modifications"){
			$self->{'onmutation'}=0;
		}elsif ($element->{'Name'} eq "neutral_losses") {
			$self->{'onnl'}=0;
		}elsif ($element->{'Name'} eq "spec2nl") {
			$self->{'onspec2nl'}=0;
		}elsif ($element->{'Name'} eq "specificity") {
			$self->{'onspec'}=0;
		}elsif ($element->{'Name'} eq "alt_names") {
			$self->{'onaltNames'}=0;
		}elsif ($element->{'Name'} eq "positions") {
			$self->{'onposition'}=0;
		}
		elsif ($element->{'Name'} eq "classifications") {
			$self->{'onclassification'}=0;
		}
	}
}

1;

####>Revision history
# 4.0.6 [Fix] bug in &getProtInfo when called with multiple analyses (PP 21/06/19)
# 4.0.5 &decodeVarMod now handles out of range PTM position with a '!' (PP 18/06/19)
# 4.0.4 Add Metadata handling in &getSearchParam (VS 25/05/19)
# 4.0.3 Change in SQL query in &getModificationIDfromString & improved verbose mode in &getProtInfo (PP 17/05/19)
# 4.0.2 &getProtInfo tries to fetch protein info from DB if no match found in databank file & update of &getSearchParam for MaxQuant and &convertVarModString (PP 21/02/19)
# 4.0.1 Added ambigous modification sites position in decodeVarMod (VS 14/12/18)
# 4.0.0 Handles multiple ambigous modification sites in decodeVarMod (VS 13/12/18)
# 3.9.9 Fixed PD search parameters displaying (VS 22/11/18)
# 3.9.8 Added TDA search parameters (VS 20/11/18)
# 3.9.7 Added protein sequence recovering in getProtInfo (VS 16/11/2018)
# 3.9.6 &updateAnalysisHistory takes userID as optional 5th argument (PP 07/11/18)
# 3.9.5 Improvement in &sortSmart (PP 02/11/18)		
# 3.9.4 Improvement in &cleanNumericalParameters (PP 02/11/18)
# 3.9.3 Minor modification in decodeVarMod (VS 23/10/2018)
# 3.9.2 Update in &deleteUnusedMasterProteins: checks all entries in MASTER_PROTEIN if no list is provided (PP 28/09/18)
# 3.9.1 Minor modification in updateAnalysisHistory (GA 27/06/18)
# 3.9.0 &updateMatchGroups: creates and updates match groups from scratch (PP 08/06/18)
# 3.8.6 Minor change in &updateAnalysisHistory to handle undefined $userID (PP 23/05/18)
# 3.8.5 Minor change in &resize (PP 16/04/18)
# 3.8.4 Modification of createMatchGroups arguments to include visibility computation (GA 27/03/18)
# 3.8.3 &deleteUnusedMasterProteins now returns number of master proteins deleted (PP 15/03/18)
# 3.8.2 Changed &cleanExportDirectory to more generic &cleanDirectory with optional $maxAge (PP 05/02/18)
# 3.8.1 Minor modification in convertVarModString for 'AzidoBiotinPEG4 Oxyde (K)' PTM (GA 15/01/18)
# 3.8.0 Minor modif to display openswath option (&getSearchParam) (MLP 06/12/17)
# 3.7.9 New options for browseDirectory_getFiles: including SELECT instead of checkbox (PP 04/12/17)
# 3.7.8 Move createMatchGroup from runMassChroQ.pl (GA 01/12/17)
# 3.7.7 Minor modif to display library export option (&getSearchParam) (MLP 30/11/17)
# 3.7.6 File smart sorting & multi-instance compatibility for browseDirectoryXXX functions (PP 09/11/17)
# 3.7.5 Add functions to display directory content (PP 03/11/17)
# 3.7.4 Add OPENSWATH.TSV to getSearchParam (MLP 21/09/17)
# 3.7.3 Minor change in &cleanParameters (PP 19/09/17)
# 3.7.2 Modification for PD 2.2 (GA 04/09/17)
# 3.7.1 Bug fix in JS ajaxManageSaveProteins function to allow 'Save protein'/Cancel/'Save protein' circle (PP 10/08/17)
# 3.7.0 Compatible with MODIFICATION_SITE & code cleaning (PP 02/08/17)
# 3.6.1 Added &cleanExportDirectory subroutine (PP 21/02/17)
# 3.6.0 Add option to shift aa position in &toStringVariableModifications and &decodeVarMod & added &MaxQuantProbIcon (PP 13/02/17)
# 3.5.9 Bug fix in &getProtInfo for multi-analysis data struture detection & update of &getModificationIDfromString and UnimodXMLHandler for detection of label modifications (PP 18/01/17)
# 3.5.8 Update in &getSearchParam for MaxQuant mqpar.xml file (PP 29/12/16)
# 3.5.7 Update &getSearchParam for X! Tandem and PeakView quantification (SWATH) (MLP 20/12/16)
# 3.5.6 added &unzipArchive to unzip any .zip archive & getProtInfo() now compatible with directly validated analyses & update in getSearchParam() (PP 17/11/16)
# 3.5.5 Added &cleanParameters function (PP 22/09/16)
# 3.5.4 Minor change in &getModificationIDfromString to allow SPECIFICITY update when field is empty (PP 17/08/16)
# 3.5.3 Minor change to prevent uninitialzed $defIdentType in &getProtInfo (PP 01/07/16)
# 3.5.2 Minor modification in &convertVarModStringSwath (MLP 20/05/2016)
# 3.5.1 Minor modification in sql query for &getSearchParam (GA 29/04/2016)
# 3.5.0 Add &convertVarModStringSwath function (MLP 28/04/2016)
# 3.4.9 Update &getSearchParam (GA 28/04/16)
# 3.4.8 Minor modification of &getProtInfo to deal with MaxQuant contaminant db and decodeVarMod to deal with ambiguities (GA 22/04/16)
# 3.4.7 Detects and moves UniProt isoform tag from ACC to ID if ID is selected in &getProInfo (PP 08/02/15)
# 3.4.6 Added &decodeVarMod for progressive replacement of &toStringVariableModifications (PP 20/01/16)
# 3.4.5 Adds display_pos info in &getItemInfo, default=0 (PP 10/12/15)
# 3.4.4 Minor change in &updateAnalysisHistory to prevent deletion of phosphoRS from validation history after "Clear all" (PP 10/08/15)
# 3.4.3 Merged with 3.4.2 from GA (PP 06/08/15)
# 3.4.2 Change tables used in extractSpectrumMSF for PD 2.0 version (GA 12/06/15)
# 3.4.1b Minor modif for better extraction of PARAGON.XML quantification parameter in &getSearchParam (PP 05/04/14)
# 3.4.1 Minor modif due to split-mode that changes $msfFileName (GA 27/02/15)
# 3.4.0 Fix bug due to new output of &extractSpectrumMSF (PP 18/02/15)
# 3.3.9 Get source rawfile from promsMod::extractSpectrumMSF (GA 28/01/15)
# 3.3.8 Improved detection of (anaID/mod)ProtID format in &printAjaxManageSaveProteins JS code (PP 23/01/15)
# 3.3.7 &extractSpectrumMSF unzips spectrum directly from string using IO::Uncompress::Unzip<BR>& accepts a filePath or a dbi connection<BR>&preparePtmString for better sortSmart of modified proteins (PP 26/12/14)
# 3.3.6 added PATHWAY_ANALYSIS to type, getProjectID, getItemChild and getItemParent functions (SL 21/10/14)
# 3.3.5 Moved &getQuantificationParameters & &writeQuantifParameterFiles to promsQunaitf.pm (PP 14/10/14)
# 3.3.4 Added new quantif parameters for Super/Simple & LabelFree ratios (PP 29/09/14)
# 3.3.3 Added printAjaxManageSaveProteins function (PP 26/08/14)
# 3.3.2 Bug fix in &getProtInfo for MSF created by 3.3.1 upgrade (PP 02/06/14)
# 3.3.1 Added new R parameters in &writeQuantifParameterFiles for SuperRatio quantification<BR>and bug fix for processing of local Mascot databank in &getProtInfo (PP 27/05/14)
# 3.3.0 add EXPLORANA to getItemInfo, getItemChild and getItemParent functions (SL 19/05/14)
# 3.2.9 PD version explicitly retrieved from last entry in SchemaInfo table (PP 12/05/14)
# 3.2.8 Bug fix in &getSearchParam for multi-db + moved &getIdentifierTypes from promsConfig.pm to promsMod.pm<BR>Bug fix in handling of Mascot DB when local access (PP 26/03/14)
# 3.2.7 Uses File::Path::rmtree instead of remove_tree (PP 10/03/14)
# 3.2.6 Change R name.grp values to stateX in &getQuantificationParameters (PP 26/02/14)
# 3.2.5 Removed most 'use package' (PP 13/11/13)
# 3.2.4 Moved Mascot databank update from &getProtInfo to &updateMascotDB (PP 07/11/13)
# 3.2.3 Minor modif of &getProtInfo due to update of myproms4databanks.pl (PP 04/11/13)
# 3.2.2 Minor modif in &toStringVariableModifications (GA 29/10/13)
# 3.2.1 Minor modif to get version of PD (PP 22/10/13)
# 3.2.0 Bug fix in &deleteUnusedSpecies (PP 10/10/13)
# 3.1.9 Minor modification for N-Term modifications in &toStringVariableModifications (GA 27/09/13)
# 3.1.8 Change in &toStringVariableModifications to handle more unusual cases (PP 26/09/13)
# 3.1.7 Change in &toStringVariableModifications for empty positions and &getSearchParam for Phenyx searches (GA 25/09/13)
# 3.1.6 Change in &convertVarModString for termini specificity (PP 24/09/13)
# 3.1.5 Prevents deletion of reference species in &deleteUnusedSpecies (PP 24/09/13)
# 3.1.4 New cases in updateValidationHistory (FY 17/09/13)
# 3.1.3 'no_user' access to &promsConfig::getServerInfo in &removeFromMultiAna for usage by .pl scripts (PP 13/09/13)
# 3.1.2 update &getItemInfo for XIC quantif (PP 09/09/13)
# 3.1.1 Minor bug correction in &getVariableModifications SQL query (PP 06/09/13)
# 3.1.0 Minor bug correction for ambiguous positions in toStringVariableModifications function (GA 02/09/13)
# 3.0.9 Modification of &applyVisibilityRule to prevent change in visibility in analyses with quantif or GO data (PP 20/08/13)
# 3.0.8 Merge of 3.0.6c & 3.0.7 (PP 24/07/13)
# 3.0.6c Change to &getModificationIDfromString for optimisation (PP 16/07/13)
# 3.0.6b Change in &toStringVariableModifications to allow modif skipping (PP 11/07/13)
# 3.0.6 Replacement of hard-coded mascot URL in &getProInfo (PP 04/07/13)
# 3.0.5 Corrections in &convertVarModString and &toStringVariableModifications (PP 18/06/13)
# 3.0.4 &getModificationIDfromString also extracts fixed/variable modifications from quantitation section (PP 04/06/13)
# 3.0.3 Minor correction in function &toStringVariableModifications specificity match for PTMs (GA 03/06/13)
# 3.0.2 Add IS_LABEL to &getModificationIDfromString function (PP 31/05/13)
# 3.0.1 Minor correction in &getModificationIDfromString to remove warning when uc($synonymesUC) is not defined (GA 28/05/13)
# 3.0.0 Add &toStringVariableModifications to create a string representation of a vmod (GA 23/05/13)
# 2.9.9 wget -> LWP code cleaning & change in &getModificationIDfromString (PP 21/05/13)
# 2.9.8 Minor modification in getModificationIDfromString for specificity matching (GA 21/05/13)
# 2.9.7 Bug correction in getModificationIDfromString for SYNONYMS matching<BR>Updated getVariableModifications from promsConfig version (GA 03/05/13)
# 2.9.6 Minor update in &formatDate (PP 29/04/13)
# 2.9.5 Minor bug correction for getModificationIDfromString (GA 23/04/13)
# 2.9.4 Minor modification for &convertVarModString, add * + - = distinction between Protein N/C-term and N/C-term (GA 17/04/13)
# 2.9.3 Additional parameters in &writeQuantifParameterFiles (PP 16/04/13)
# 2.9.2 Add UnimodXMLHandler and &getModificationIDfromString for new modification table in myProMS (GA 16/04/13)
# 2.9.1 Change in &extractData to include PARAGON spectrum for storeSpectrum.cgi -> refSpec (GA 22/03/13)
# 2.9.0 Check undef on $lastValType in &updateAnalysisHistory (PP 20/03/13)
# 2.8.9 Change in &removeFromMultiAna, &extractData and &extractSpectrumMSF to account for new send2Biologist v2.6.6 behavior (PP 07/03/13)
# 2.8.8 Update of &getSearchParam for new file path when validation is ended (PP 04/03/13)
# 2.8.7 &convertVarModString handles uncertainty for -term position & &getLinkClass handles ghost proteins & &getProteinClass (PP 06/02/13)
# 2.8.6 &deleteValidData moved to deleteProjectItem.cgi<BR>&getSearchParam and &getProtInfo compatible with multi-db search (PP 13/12/12)
# 2.8.5 Handles master proteins deletion. DAVID gene mapping is obsolete (PP 21/11/12)
# 2.8.4 Minor modification of getItemInfo to correct a bug for compareQuantification.cgi (GA 19/11/12)
# 2.8.3 Modification of extractData for PDM files (GA 13/11/12)
# 2.8.2 Merge 2.8.1GA & 2.8.1FY + modification of getSearchParam for new format handling (GA 12/11/12)
# 2.8.1fy Add functions to deal with categories (FY 09/11/12)
# 2.8.1ga Modification of getItemInfo for XIC extractions (GA 22/10/12)
# 2.8.1 Merge 2.8.0PP & 2.8.0GA (GA 15/10/12)
# 2.8.0PP Correction in &getUserInfo for manager (PP 18/09/12)
# 2.8.0ga Modification of functions (like getSearchParam) for PARAGON (GA 05/09/12)
# 2.8.0m Added project status & LWP call in &getItemInfo & GO_ANALYSIS in &getProjectID (PP 12/09/12)
# 2.8.0 New &getProtInfo subroutine for multi-analysis Databank scan (PP 13/08/12)
# 2.7.9 Improve access right definition for manager in &getUserInfo (PP 31/07/12)
# 2.7.8 Minor update: parseRules safe URL conversion (GA 26/07/12)
# 2.7.7 Minor modification of getSearchParam function (call by editProjectItem) (GA 08/06/12)
# 2.7.6 Minor improvement in navigation tree functions (PP 21/05/12)
# 2.7.5 Added Quantitification in &getSearchParam  (PP 07/05/12)
# 2.7.4 'no_user' in %promsPath declaration in &getProtInfo (PP 04/05/12)
# 2.7.3 Handles designless quantifs in &getItemInfo (PP 17/04/12)
# 2.7.2 Remove DESIGN & EXPCONDITION from & getCheckedItems (PP 06/04/12)
# 2.7.1 Add PhosphoRS analyses in validation history (FY 30/03/12)
# 2.7.0 Merge 2.5.8GA (Label-free) & 2.6.9PP (02/04/12)
# 2.6.9 &getUserInfo: no manager/project link if workgroups not set (PP 17/03/12)
# 2.6.8 Update popup() javascript function to correctly display the popup in scrolled windows (FY 19/01/12)
# 2.6.7 Manages peptide quantification data during analysis deletion (PP 24/11/11)
# 2.6.6 Added QUANTIFICATION to getProjectItem function & changes in proxy string (PP 14/10/11)
# 2.6.5 Labeled quantification parameters management functions (PP 23/09/11)
# 2.6.4 Manager, bug unit values in convertVarModString & quantif functions (PP 05/09/11)
# 2.6.3 add Mascot IDs to &getUserInfo returned data (PP 22/08/11)
# 2.6.2 Add GO Analysis type & bug in updateValidationHistory (FY 07/07/11)
# 2.6.1 minor change in convertVarModString to remove starting _ (Cochin) (PP 16/06/11)
# 2.6.0 full datafile conservation management (PP 01/06/11)
# 2.5.9 minor change in convertVarModString for parsing varMods in PMF format (no position info) (PP 10/05/11)
# 2.5.7 Add convertVarModString function (PP 29/04/11)
# 2.5.6 End validation management in updateAnalysisHistory function and minor bugs (FY 18/04/11)
# 2.5.5 Merges 2.5.4PP/2.5.4FY & changes in sub updateAnalysisHistory (PP 11/04/11)
# 2.5.4 Minor change in HTMLcompatible function (PP 04/03/11)<BR>Add updateAnalysisHistory function (FY 02/2011)
# 2.5.3 Minor change in getChildrenList sub: no child TYPE if no children (PP 21/02/11)
# 2.5.2 Minor improvement in getProtInfo sub (PP 05/01/11)
# 2.5.1 Minor bug in getProtInfo sub (PP)
