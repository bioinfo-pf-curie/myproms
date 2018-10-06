#!/usr/local/bin/perl -w

################################################################################
# editVarModification.cgi      1.3.1                                           #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Edits variable modifications position & peptide comments                     #
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
$|=1;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use promsConfig;
use promsMod;
use strict;
use phosphoRS;

#print header; warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};

####################
#    Parameters    #
####################
my $call=param('CALL');
my $varModID=param('ID');

#######################
#    Process Form     #
#######################
if (param('revert')) {
	#my $newComment=(param('editTextForm'))? param('editTextForm') : '';
	#$newComment =~ s/,/;/g;
	my $selModID=(param('revert') > 0)? param('revert') : 0;

	####>Connect to the database<####
	my $dbh=&promsConfig::dbConnect;

	if ($call eq 'listRanks') {
		my ($tag,$queryID,$qNum,$rank)=split('_',$varModID);
		if ($selModID) {
			my ($refPosString)= $dbh->selectrow_array("SELECT REF_POS_STRING FROM QUERY_MODIFICATION WHERE ID_MODIFICATION=$selModID AND ID_QUERY=$queryID AND PEP_RANK=$rank");
			$dbh->do("UPDATE QUERY_MODIFICATION SET POS_STRING='$refPosString',REF_POS_STRING=NULL WHERE ID_MODIFICATION=$selModID AND ID_QUERY=$queryID AND PEP_RANK=$rank");
		}
		else { # revert all
			my $sthPTMpos=$dbh->prepare("SELECT ID_MODIFICATION,REF_POS_STRING FROM QUERY_MODIFICATION WHERE ID_QUERY=$queryID AND PEP_RANK=$rank");
			$sthPTMpos->execute;
			while (my ($modID,$refPosString) = $sthPTMpos->fetchrow_array) {
				next unless $refPosString;# Position(s) of this modification was never changed so no need to update it !
				$dbh->do("UPDATE QUERY_MODIFICATION SET POS_STRING='$refPosString',REF_POS_STRING=NULL WHERE ID_MODIFICATION=$modID AND ID_QUERY=$queryID AND PEP_RANK=$rank");
			}
			$sthPTMpos->finish;
			#my ($rankInfo)=$dbh->selectrow_array("SELECT INFO_PEP$rank FROM QUERY_VALIDATION WHERE ID_QUERY=$queryID");
			#$rankInfo=~s/COM=[^,]*,?//; # remove old comments
			#$rankInfo.="COM=$newComment," if $newComment;
			#$dbh->do("UPDATE QUERY_VALIDATION SET INFO_PEP$rank='$rankInfo' WHERE ID_QUERY=$queryID");
		}
	}
	elsif ($call eq 'sequenceView') {
		my ($tag,$pepID)=split('_',$varModID);
		my $sthUpPep=$dbh->prepare("UPDATE PEPTIDE_MODIFICATION SET POS_STRING=?,REF_POS_STRING=? WHERE ID_MODIFICATION=? AND ID_PEPTIDE=$pepID");
		if ($selModID) {
			my ($refPosData)= $dbh->selectrow_array("SELECT REF_POS_STRING FROM PEPTIDE_MODIFICATION WHERE ID_MODIFICATION=$selModID AND ID_PEPTIDE=$pepID");
			my ($posString)=$refPosData=~/^([^#]+)/; # ignore MaxQuant ##PRB_MQ=xxxx if any
			(my $probData=$refPosData)=~s/^[^#]+//; $probData=undef unless $probData;
			$sthUpPep->execute($posString,$probData,$selModID);
		}
		else { # revert all
			#my $sthPTMpos=$dbh->prepare("SELECT ID_MODIFICATION,REF_POS_STRING FROM PEPTIDE_MODIFICATION WHERE ID_PEPTIDE=$pepID AND MODIF_TYPE='V'");
			my $sthPTMpos=$dbh->prepare("SELECT ID_MODIFICATION,REF_POS_STRING FROM PEPTIDE_MODIFICATION WHERE ID_PEPTIDE=$pepID"); # Change on 17/09
			$sthPTMpos->execute;
			while (my ($modID,$refPosData) = $sthPTMpos->fetchrow_array) {
				next unless $refPosData;# Position(s) of this modification was never changed so no need to update it !
				my ($posString)=$refPosData=~/^([^#]+)/; # ignore MaxQuant ##PRB_MQ=xxxx if any
				(my $probData=$refPosData)=~s/^[^#]+//; $probData=undef unless $probData;
				$sthUpPep->execute($posString,$probData,$modID);
			}
			$sthPTMpos->finish;
			#my $qComment=$dbh->quote($newComment);
			#$dbh->do("UPDATE PEPTIDE SET COMMENTS=$qComment WHERE ID_PEPTIDE=$pepID");
		}
		$sthUpPep->finish;
	}

	$dbh->commit;
	$dbh->disconnect;

	print header(-'charset' => 'UTF-8');
	print qq
|<HTML><HEAD>
<SCRIPT type="text/javascript">
if (opener) {opener.window.location.reload();}
window.location="./editVarModification.cgi?CALL=$call&ID=$varModID";
</SCRIPT>
</HEAD></HTML>
|;
	exit;
}
elsif (param('save')) {
	my $numVarMods=param('numVarMods');
	my @varModIDs=param('varModIDs');
	my @sequence = split(//,param('sequence'));
	my $seqLength=length(param('sequence'));
	my $newComment=(param('editTextForm'))? param('editTextForm') : '';
	$newComment =~ s/,/;/g;

	#Sequence composition
	my %seqComposition;
	foreach my $aa (@sequence) {$seqComposition{$aa}++;}

	####>Connect to the database<####
	my $dbh=&promsConfig::dbConnect;

	my ($tag,$queryID,$qNum,$rank,$pepID,$anaID,$pepLength);

	if ($call eq 'listRanks') {
		($tag,$queryID,$qNum,$rank)=split('_',$varModID);
		($anaID,my $rankInfo)=$dbh->selectrow_array("SELECT ID_ANALYSIS,INFO_PEP$rank FROM QUERY_VALIDATION WHERE ID_QUERY=$queryID");
		$rankInfo=~s/COM=[^,]*,?//; #remove old comments
		$rankInfo.="COM=$newComment," if $newComment;
		my $pepSequ=
		$dbh->do("UPDATE QUERY_VALIDATION SET INFO_PEP$rank='$rankInfo' WHERE ID_QUERY=$queryID");
	}
	else {
		($tag,$pepID)=split('_',$varModID);
		my $qComment=$dbh->quote($newComment);
		$dbh->do("UPDATE PEPTIDE SET COMMENTS=$qComment WHERE ID_PEPTIDE=$pepID");
		($anaID) = $dbh->selectrow_array("SELECT ID_ANALYSIS FROM PEPTIDE WHERE ID_PEPTIDE=$pepID");
	}

	foreach my $modID (@varModIDs) {
#print "$varModNames[$v] -> ";
		#next if $possibleAAs=~/-/; # N/C-term cannot be edited
		my $occurence=param("varModOccurences_$modID");
		my @posList=sort{$b<=>$a} param("chkVarMod_$modID");
		my ($anaSpecificity)=$dbh->selectrow_array("SELECT SPECIFICITY FROM ANALYSIS_MODIFICATION WHERE ID_ANALYSIS=$anaID AND ID_MODIFICATION=$modID AND MODIF_TYPE='V'");
		for (my $i=0; $i<= $#posList ; $i++) {
			if ($posList[$i] eq '-0') {
				$posList[$i]=($anaSpecificity =~ /-/)? '--':'-=';
			}elsif ($posList[$i]==0) {
				$posList[$i]=($anaSpecificity =~ /-/)? '-':'=';
			}elsif ($posList[$i]>$seqLength) {
				$posList[$i]=($anaSpecificity =~ /\+/)? '+':'*';
			}elsif ($posList[$i]< -$seqLength) {
				$posList[$i]=($anaSpecificity =~ /\+/)? '-+':'-*';
			}
		}

		#foreach my $unk ((1+scalar @posList)..$occurence) { # complete with unknown positions
		#	push @posList,-1;
		#}
		#
		#foreach my $pos (0..$#posList) { # prevent to put wrong positions for N/C-Term positions
		#	$posList[$pos]=1 if $posList[$pos] == 0;
		#	$posList[$pos]=$seqLength if $posList[$pos] == $seqLength+1;
		#}

		# Checking for residue type change for phosphorylations
		#if (($possibleAAs eq 'ST' || $possibleAAs eq 'Y') && $varModNames[$v] =~ /Phospho/i) {
		#	foreach my $pos (@posList) {
		#		if(($sequence[$pos-1] eq 'Y' && $possibleAAs eq 'ST')
		#		   || (($sequence[$pos-1] eq 'S' || $sequence[$pos-1] eq 'T') && $possibleAAs eq 'Y')){
		#			$possibleAAs = 'STY';
		#		}
		#	}
		#}
		#>Phosphorylation
		#if ($varModNames[$v] =~ /Phospho/i) {
		#	$possibleAAs='';
		#	$possibleAAs.='ST' if ($seqComposition{'S'} || $seqComposition{'T'});
		#	$possibleAAs.='Y' if $seqComposition{'Y'};
		#}
		##>Methylation
		#elsif ($varModNames[$v] =~ /^(Di|Tri)*Methyl/i) {
		#	$possibleAAs='';
		#	$possibleAAs.='K' if $seqComposition{'K'};
		#	$possibleAAs.='R' if $seqComposition{'R'};
		#}
		##>Acetylation
		#elsif ($varModNames[$v] =~ /^Acetyl/i) {
		#	$possibleAAs=($posList[0]<=0)? 'N-term' : 'K';
		#}

		my $posStrg='';
		my $nbMod=0;
		foreach my $pos (@posList) {
			if ($pos =~ /-(\d+|=|-|\+|\*)/ && $nbMod == $occurence-1) {
				$posStrg.=join(';',@posList[$nbMod..$#posList]);
				last;
			}
			$posStrg.=($nbMod == $occurence-1)?"$pos":"$pos.";
			$nbMod++;
		}

		#$varModNames[$v]=~s/\(.+\)/\($posStrg\)/; # converts PMF format into MIS
#print "$varModNames[$v]<BR>\n";

		if ($call eq 'listRanks') {
			my ($posString,$refPosString)=$dbh->selectrow_array("SELECT POS_STRING,REF_POS_STRING FROM QUERY_MODIFICATION WHERE ID_MODIFICATION=$modID AND ID_QUERY=$queryID AND PEP_RANK=$rank");
			if ($posString ne $posStrg) {
				$refPosString=$posString unless $refPosString;
				$dbh->do("UPDATE QUERY_MODIFICATION SET POS_STRING='$posStrg',REF_POS_STRING='$refPosString' WHERE ID_MODIFICATION=$modID AND ID_QUERY=$queryID AND PEP_RANK=$rank");
			}
		}
		elsif ($call eq 'sequenceView') {
			my ($posString,$refPosData)=$dbh->selectrow_array("SELECT POS_STRING,REF_POS_STRING FROM PEPTIDE_MODIFICATION WHERE ID_MODIFICATION=$modID AND ID_PEPTIDE=$pepID");
			if ($posString ne $posStrg) {
				my ($refPosString)=$refPosData=~/^([^#]+)/; # ignore MaxQuant ##PRB_MQ=xxxx if any
				$refPosString=$posString unless $refPosString;
				(my $probData=$refPosData)=~s/^[^#]+//;
				$refPosString.=$probData;
				$dbh->do("UPDATE PEPTIDE_MODIFICATION SET POS_STRING='$posStrg',REF_POS_STRING='$refPosString' WHERE ID_MODIFICATION=$modID AND ID_PEPTIDE=$pepID");
			}
		}
	}
	#my $varModString=' + '.join(' + ',@varModNames);
#print "$varModString<BR>\n";
#exit;

	$dbh->commit;
	$dbh->disconnect;

	print header(-charset=>'utf-8');
	print qq
|<HTML><HEAD>
<SCRIPT type="text/javascript">
if (opener) {opener.window.location.reload();}
window.close();
</SCRIPT>
</HEAD></HTML>
|;
	exit;
}

##################
#      Main      #
##################

####>Connect to the database<####
my $dbh=&promsConfig::dbConnect;

my ($pepSequence,$comments,$anaID,$prsStatus,$prsOutputFile,$fileFormat);
my $sthGetVM;
if ($call eq 'listRanks') {
	my ($tag,$queryID,$qNum,$rank)=split('_',$varModID);
	my ($rankInfo)=$dbh->selectrow_array("SELECT INFO_PEP$rank FROM QUERY_VALIDATION WHERE ID_QUERY=$queryID");
	($pepSequence)=($rankInfo=~/SEQ=(\w+)/);
	$sthGetVM=$dbh->prepare("SELECT ID_MODIFICATION,POS_STRING,REF_POS_STRING FROM QUERY_MODIFICATION WHERE ID_QUERY=$queryID AND PEP_RANK=$rank");
	#($varModString)=($rankInfo=~/VMOD=([^,]+)/);
	#($refModString)=($rankInfo=~/RMOD=([^,]+)/);
	($comments)=($rankInfo=~/COM=([^,]*),?/);
	if($rankInfo =~ /PRS=(\d);[^;,]*;[^;,]*,/){
		$prsStatus = $1;
	}
	$fileFormat='-'; # anything but MAXQUANT.DIR
	($anaID) = $dbh->selectrow_array("SELECT ID_ANALYSIS FROM QUERY_VALIDATION WHERE ID_QUERY=$queryID");
}
elsif ($call eq 'sequenceView') {
	my ($tag,$pepID)=split('_',$varModID);
	($pepSequence,$comments,my $data,$anaID,$fileFormat)=$dbh->selectrow_array("SELECT PEP_SEQ,P.COMMENTS,DATA,P.ID_ANALYSIS,FILE_FORMAT FROM PEPTIDE P,ANALYSIS A WHERE P.ID_ANALYSIS=A.ID_ANALYSIS AND ID_PEPTIDE=$pepID");
	$sthGetVM=$dbh->prepare("SELECT ID_MODIFICATION,POS_STRING,REF_POS_STRING FROM PEPTIDE_MODIFICATION WHERE ID_PEPTIDE=$pepID AND MODIF_TYPE='V'");
	if ($data) {
		if ($data =~ /PRS=(\d);/) {
			$prsStatus = $1;
		}
		#($refModString)=($data=~/REF_MOD=([^#]+)/);
	}
}

###> Retrieve MODIFICATIONS information such as name, specificity, position on the peptide and ref_pos is defined earlier
#my $sthPTMInfo=$dbh->prepare("SELECT PSI_MS_NAME,INTERIM_NAME,SPECIFICITY FROM MODIFICATION WHERE ID_MODIFICATION=?");
my $sthPTMInfo=$dbh->prepare("SELECT PSI_MS_NAME,INTERIM_NAME,SYNONYMES,AM.SPECIFICITY FROM MODIFICATION M,ANALYSIS_MODIFICATION AM WHERE M.ID_MODIFICATION=AM.ID_MODIFICATION AND AM.ID_ANALYSIS=$anaID AND M.ID_MODIFICATION=?");
$sthGetVM->execute;
my %ptmInfo;
while (my ($modID,$posStrg,$refPosData) = $sthGetVM->fetchrow_array) {
	$refPosData='' unless $refPosData;
	$sthPTMInfo->execute($modID);
	my ($psiName,$interName,$synonymes,$specificity)=$sthPTMInfo->fetchrow_array;
	my $modName=($psiName)? $psiName : ($interName)? $interName : ($synonymes)? $synonymes : "mod_$modID";
	$modName=~s/^##//; $modName=~s/##.*//; # in case synonymes
	@{$ptmInfo{$modID}}=($modName,$specificity,$posStrg,$refPosData);
}
$sthGetVM->finish;
$sthPTMInfo->finish;

$comments='' unless $comments;
if (defined $prsStatus) {
	my $projectID = promsMod::getProjectID($dbh, $anaID, 'ANALYSIS');
	($prsOutputFile, my $validStatus) = $dbh->selectrow_array("SELECT DATA_FILE, VALID_STATUS FROM ANALYSIS WHERE ID_ANALYSIS=$anaID");
	$prsOutputFile =~ s/\.\w{3}$/\.xml/;
	#my ($validStatus) = $dbh->selectrow_array("SELECT VALID_STATUS FROM ANALYSIS WHERE ID_ANALYSIS=$anaID");
	###> phosphoRS.pm changed the output filename for 1.7.1 modification
	if ($validStatus == 2) {
		$prsOutputFile = (-e "$promsPath{peptide}/proj_$projectID/ana_$anaID/PRS_$prsOutputFile")? "$promsPath{peptide}/proj_$projectID/ana_$anaID/PRS_$prsOutputFile" : "$promsPath{peptide}/proj_$projectID/ana_$anaID/PRS_ana_$anaID.xml";
	} else {
		$prsOutputFile = (-e "$promsPath{valid}/ana_$anaID/PRS_$prsOutputFile")?"$promsPath{valid}/ana_$anaID/PRS_$prsOutputFile" : "$promsPath{valid}/ana_$anaID/PRS_ana_$anaID.xml";
	}
}
$dbh->disconnect;

#$pepSequence.='XCMSTY'; # DEBUG

my @pepResidues=split('',$pepSequence);
my @pepAAs=('N-term');
push @pepAAs,@pepResidues;
#push @pepAAs,('C-term','Ambigous');
push @pepAAs,'C-term';
my $seqLength=length($pepSequence);
my $numVarMods=scalar keys %ptmInfo;
my (@varModNames,@varModIDs,%varModPositions,%varModOccurences);
foreach my $modID (keys %ptmInfo) {
	my ($varModName,$specificity,$posStrg,$refPosData)=@{$ptmInfo{$modID}};

	push @varModNames,$varModName;
	push @varModIDs,$modID;

	if ($posStrg) { # MIS
		foreach my $pos (split(/\./,$posStrg)) {
			#$pos=$seqLength+2 if $pos==-1; # -1 unknown position (can be repeated!!!)
			foreach my $pos2 (split(/;/,$pos)) { # For ambiguous positions...
				$pos2='-0' if ($pos2 =~/--|-=/);
				#$pos2=0 if ($pos2 =~/-\d!|=/); # Pattern matching is not working - Changed on 29/10/13
				$pos2=0 if ($pos2 =~/-(?!\d)|=/);
				$pos2=length($pepSequence)+1 if ($pos2 =~/\+|\*/);
				$pos2=-length($pepSequence)-1 if ($pos2 =~/-\+|-\*/);
				$varModPositions{$varModName}{$pos2}++;
			}
			$varModOccurences{$varModName}++;
		}
	}
	#if ($refPosStrg && $refPosStrg ne $posStrg){
	#	$refModString=1;
	#}elsif(!$refModString){ # Activate the revert button if one modification among several was changed !
	#	$refModString=0;
	#}
}


###########################
#      Starting HTML      #
###########################
my ($color1,$color2)=&promsConfig::getRowColors;
print header(-'content-encoding'=>'no');
warningsToBrowser(1);
print qq
|<HEAD>
<TITLE>Edit Modifications Position</TITLE>
<LINK rel="icon" type="image/png" href="$promsPath{images}/myProMS_icon.png">
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT type="text/javascript">
|;
&promsMod::popupInfo() if $fileFormat eq 'MAXQUANT.DIR';
print "var varModNameList=['",join("','",@varModNames),"'];\n";
print "var varModIDList=['",join("','",@varModIDs),"'];\n";
# Change the 1st of July 2013
#print qq
#|function checkForm(myForm) {
#	for (var v=0; v<varModNameList.length; v++) {
#		var chkVM=myForm['chkVarMod_'+varModIDList[v]];
#		var numChecked=0;
#		for (var i=0; i<chkVM.length; i++) {
#			if (chkVM[i].checked) {numChecked++;}
#		}
#		var chkUnk=myForm['chkUnknown_'+varModIDList[v]].checked;
#		if ((chkUnk && numChecked >= document.getElementById('varModOccurences_'+varModIDList[v]).value) \|\| (!chkUnk && numChecked != document.getElementById('varModOccurences_'+varModIDList[v]).value)) {
#			alert("ERROR: Number of modified positions for '"+varModNameList[v]+"' has changed.");
#			return false;
#		}
#	}
#	formWasSubmitted=1;
#//alert('OK'); return false;
#	return true;
#}
my $chkVMLength=$seqLength+1;# To avoid javascript function to fail when no ambiguity is possible (one oxidation site for example... Ambigous line does not exist)
print qq
|function checkForm(myForm) {
	for (var v=0; v<varModNameList.length; v++) {
		var chkVM=myForm['chkVarMod_'+varModIDList[v]];
		var numCheckedSure=0;
		var numCheckedUnsure=0;
		for (var i=0; i<$chkVMLength; i++) {
			if (chkVM[i].checked) {
				numCheckedSure+=1;
			}
			if (chkVM.length/2 == $chkVMLength+1) {
				if ((chkVM[i+chkVM.length/2].checked) {
					numCheckedUnsure+=1;
				}
				if (chkVM[i].checked && chkVM[i+$chkVMLength+1].checked) {
					var pos=chkVM[i].value;
					if (chkVM[i].value == 0) {
						pos='N-term';
					}
					if (chkVM[i].value == $chkVMLength) {
						pos='C-term';
					}
					alert("ERROR: "+varModNameList[v]+" cannot be sure and ambigous on the same residue "+pos+" !");
					return false;
				}
				
			}
		}
		var vmodOcc= document.getElementById('varModOccurences_'+varModIDList[v]).value;
		if (numCheckedSure > vmodOcc \|\| (numCheckedSure==vmodOcc && numCheckedUnsure>0)) {
			alert("ERROR: You have checked more '"+varModNameList[v]+"' than you should. (only "+vmodOcc+" found on this sequence)");
			return false;
		}
		if (numCheckedSure + numCheckedUnsure < vmodOcc) {
			alert("ERROR: You have checked less '"+varModNameList[v]+"' than you should. ("+vmodOcc+" found on this sequence)");
			return false;
		}
	}
	return false;
}
function closeWindow() {
	if (opener) {
		opener.document.getElementById(opener.selectedVarModId).style.color='#000000';
		opener.selectedVarModId=null;
	}
	window.close();
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Edit Modification Positions</FONT>
<BR>
<BR>
<FORM name="editVarMod" onsubmit="return(checkForm(this));" method="POST">
<INPUT type="hidden" name="CALL" value="$call">
<INPUT type="hidden" name="ID" value="$varModID">
<INPUT type="hidden" name="numVarMods" value="$numVarMods">
<INPUT type="hidden" name="sequence" value="$pepSequence">
|;
#foreach my $varModName (@varModNames) {
#	print "<INPUT type=\"hidden\" name=\"varModNames\" value=\"$varModName\">\n";
#}
foreach my $modID (keys %ptmInfo) {
	print "<INPUT type=\"hidden\" name=\"varModIDs\" value=\"$modID\">\n";
}
print qq
|<TABLE cellspacing=0 cellpadding=2 border=0>
<TR><TD colspan=3></TD>
|;
foreach my $i (1..($#pepAAs-1)) {print "<TH class=\"font11\">$i</TH>";}
print qq
|<TD colspan=2></TD></TR>
<TR bgcolor="$color2"><TH class="title3 rbBorder" nowrap>&nbsp;Modifications&nbsp;</TH><TH class="title3 rbBorder" nowrap>&nbsp;Confidence&nbsp;</TH>
|;
foreach my $i (0..$#pepAAs) {
	my $bClass=($i==$#pepAAs)? 'bBorder' : 'rbBorder';
	print "<TH class=\"title3 $bClass\">&nbsp;$pepAAs[$i]&nbsp;</TH>";
}
print "</TR>\n";
my $bgColor=$color1;
my $editableVarMods=0;
my $numRevertable=0;
foreach my $modID (sort{$ptmInfo{$a}[0] cmp $ptmInfo{$b}[0]} keys %ptmInfo) {
	my ($varModName,$specificity,$posStrg,$refPosData)=@{$ptmInfo{$modID}};

	my ($refPosStrg)=$refPosData=~/^([^#]+)/; # ignore MaxQuant ##PRB_MQ=xxxx if any
	(my $probStrg=$refPosData)=~s/^[^#]*##PRB_MQ=//;
	my %maxQuantProb;
	if ($probStrg) {
		while ($probStrg =~ /([^,:]+):([^,]+)/g) {
			$maxQuantProb{$1}=$2;
		}
	}

	# 1st pass
	my (%checkablePos,%checkedPos);
	foreach my $i (0..$#pepAAs) {
		#$checkablePos{$i}=1 unless ($possibleAAs !~ /$pepAAs[$i]/ || ($possibleAAs=~/-/ && $pepAAs[$i] !~ /-/)); # N-term vs N
		#my $posCode=($pepAAs[$i]=~/N-term/i)? '-' : ($pepAAs[$i]=~/C-term/i)? '\+' : $pepAAs[$i];
		#my $posCode=($i==0)? '-' : ($i==$#pepAAs)? '\+' : $pepAAs[$i]; # N/C-term
		my $posCode=($i==0)? '-|\=' : ($i==$#pepAAs)? '\+|\*' : $pepAAs[$i]; # N/C-term
		if ($specificity=~/$posCode/) {
			$checkablePos{$i}=1;
			if ($fileFormat eq 'MAXQUANT.DIR' && !$maxQuantProb{$i}) {
				$maxQuantProb{$i}=($probStrg || !$varModPositions{$varModName}{$i})? 0 : 1; # complete prob not recorded in REF_POS_STRING
			}
		}
		$checkedPos{$i}++ if $varModPositions{$varModName}{$i};
		$checkedPos{"-$i"}++ if $varModPositions{$varModName}{"-$i"};
	}
	my $numChecked=scalar keys %checkedPos;
	my $editable=($numChecked < scalar keys %checkablePos || $posStrg =~ /\-/)? 1 : 0;
	$editableVarMods+=$editable;
	my $revertable=($refPosStrg && $refPosStrg ne $posStrg)? 1 : 0;
	$numRevertable+=$revertable;
	my $rowSpan=($fileFormat eq 'MAXQUANT.DIR')? 2 : 1;
	print qq
|<TR bgcolor="$bgColor">
<TH valign=top align=right rowspan=$rowSpan nowrap>&nbsp;$varModName [x$varModOccurences{$varModName}]:
	<INPUT type="hidden" name="varModOccurences_$modID" id="varModOccurences_$modID" value="$varModOccurences{$varModName}"/></TH>
|;
	if ($fileFormat eq 'MAXQUANT.DIR') {
		print "<TH align=right>MaxQuant:</TH>"; # Confidence
		foreach my $i (0..($#pepAAs)) {
			print '<TH>';
			print &promsMod::MaxQuantProbIcon($maxQuantProb{$i},{text=>'&nbsp;&nbsp;&nbsp;&nbsp;',popupText=>'<B>Probability='.($maxQuantProb{$i}*100).'%</B>'}) if defined $maxQuantProb{$i};
			print '</TH>';
		}
		print "</TR>\n<TR bgcolor=\"$bgColor\">\n";
	}

	print "<TH align=right>Sure:</TH>";
	#>2nd pass (Sure)
	foreach my $i (0..($#pepAAs)) {
		my ($chkStatus,$visibilityStrg);
		if ($checkedPos{$i}) {
			$chkStatus='checked';
			$visibilityStrg='';
		}
		else {
			$chkStatus='';
			$visibilityStrg=($checkablePos{$i})? '' : 'style="visibility:hidden"';
		}
		print "<TH><INPUT type=\"checkbox\" name=\"chkVarMod_$modID\" value=\"$i\" $chkStatus $visibilityStrg";
		print ' onclick="this.checked=!this.checked;"' unless $editable; # forces box to stay checked
		print "/></TH>";
	}
	#>2nd pass (Ambigous)
	print "</TR>\n";
	if ($editable) {
		print "<TR bgcolor=\"$bgColor\">\n<TH>";
		print qq |&nbsp;<INPUT type="button" name="revert" value="Revert" onclick="if (confirm('Revert position(s) of $varModName to original state?')) {window.location='./editVarModification.cgi?CALL=$call&ID=$varModID&revert=$modID';}">&nbsp;| if $revertable;
		print "</TH>\n<TH align=right>Ambigous:</TH>\n";
		foreach my $i (0..($#pepAAs)) {
			my ($chkStatus,$visibilityStrg);
			if ($checkedPos{"-$i"} || ($i == -1 && $checkedPos{"-$i"} > 1 && !$checkedPos{$i})) { # For old ambigous storage...
				$chkStatus='checked';
				$visibilityStrg='';
			}
			else {
				$chkStatus='';
				$visibilityStrg=($checkablePos{$i})? '' : 'style="visibility:hidden"';
			}
			print "<TH><INPUT type=\"checkbox\" name=\"chkVarMod_$modID\" value=\"-$i\" $chkStatus $visibilityStrg/></TH>";
		}
		print "</TR>\n";
	}
	$bgColor=($bgColor eq $color2)? $color1 : $color2;
}

my $colSpan=2+scalar @pepAAs;
my $disSummit=($editableVarMods)? '' : 'disabled';
print qq
|<TR><TD colspan=$colSpan></TD></TR>
<TR bgcolor="$color2"><TD colspan=$colSpan><TABLE><TR>
<TH valign=top align=left nowrap>&nbsp;Comments :</TH><TD><TEXTAREA cols=54 name='editTextForm' rows=3 colspan=$colSpan>$comments</TEXTAREA></TD>
<TH width=100%>&nbsp;<INPUT type="submit" name="save" value="Save" style="width:90px" $disSummit>&nbsp;
<BR>&nbsp;<INPUT type="button" value="Cancel" style="width:90px" onclick="closeWindow();">&nbsp;|;
print qq |<BR>&nbsp;<INPUT type="button" name="revert" value="Revert All" style="width:90px" onclick="if (confirm('Revert position(s) of all modifications to original state?')) {window.location='./editVarModification.cgi?CALL=$call&ID=$varModID&revert=-1';}">&nbsp;| if $numRevertable > 1;
print qq |</TH>
</TR></TABLE></TD></TR>
</TABLE>
</FORM>
|;

#############
# phosphoRS #
#############
if (defined $prsStatus) {
	print qq
|<BR>
<FONT class="title">PhosphoRS Results</FONT>
<DIV id="waitDiv">
<IMG src="$promsPath{images}/scrollbarGreen.gif"/>
</DIV>
|;
	my $prsTable;
	if($prsStatus < 4){
		my $bgcolor = $color2;
		my ($qNum,$rank);
		if($call eq 'listRanks'){
			($qNum,$rank)=(split('_',$varModID))[2,3];
		} elsif ($call eq 'sequenceView'){
			my ($tag,$pepID)=split('_',$varModID);
			my $dbh = &promsConfig::dbConnect;
			($qNum,$rank) = $dbh->selectrow_array("SELECT QUERY_NUM, PEP_RANK FROM PEPTIDE WHERE ID_PEPTIDE=$pepID");
			$dbh->disconnect;
		}
		$prsTable = qq
|<TABLE border=0 cellspacing=0 align="center">
<TR bgcolor="$bgcolor"><TH class=\"rbBorder\">&nbsp;Isoform&nbsp;</TH><TH class=\"rbBorder\">&nbsp;Position(s)&nbsp;</TH><TH class=\"rbBorder\">&nbsp;Score&nbsp;</TH><TH class=\"bBorder\">&nbsp;Probability&nbsp;</TH></TR>
|;
		foreach my $isoformDataRef (@{phosphoRS::getIsoformsFromFile($prsOutputFile,$qNum,$rank)}){
			$bgcolor = ($bgcolor eq $color1)? $color2: $color1;
			my ($prob,$score, $positionRef) = @{$isoformDataRef};
			$prob = sprintf("%.2f", 100*$prob) . "%";
			$prob = "<0.01%" if $prob eq "0.00%";
			$score = sprintf("%.2f", $score);
			my @aaString = @pepAAs[1..$#pepAAs-2];
			foreach my $pos (@{$positionRef}){
				$aaString[$pos-1] = "<FONT style=\"font-weight:bold;color:red\">$aaString[$pos-1]</FONT>";
			}
			my $aaString = join('',@aaString);
			$prsTable .= "<TR bgcolor=\"$bgcolor\"><TD align=\"center\">&nbsp;$aaString&nbsp;</TD><TD align=\"center\">&nbsp;";
			$prsTable .= join(', ', @{$positionRef});
			$prsTable .= "&nbsp;</TD><TD align=\"center\">&nbsp;$score&nbsp;</TD><TD align=\"center\">&nbsp;$prob&nbsp;</TD></TR>";
		}
		$prsTable .= "</TABLE>";
	}
	elsif ($prsStatus == 4) {
		$prsTable = "<FONT class=\"title2\">No evidence for phospho-isoform</FONT>";
	}
	print qq
|<SCRIPT language="Javascript">
document.getElementById('waitDiv').style.visibility = 'hidden';
</SCRIPT>
$prsTable
|;
}

print qq
|</CENTER>
<BR><BR>
<DIV id="divDescription" class="clDescriptionCont">
<!--Empty div-->
</DIV>
<SCRIPT type="text/javascript">if ('$fileFormat'=='MAXQUANT.DIR') {setPopup();}</SCRIPT>
</BODY>
</HTML>
|;

####>Revision history<####
# 1.3.1 Minor modification (GA 31/01/18)
# 1.3.0 Minor modification to correct warning (GA 17/10/17)
# 1.2.9 Change PRS outputfile due to 1.1.7 phosphoRS.pm modification (GA 04/09/17)
# 1.2.8 Minor modification to make $sthPTMInfo work (GA 10/08/17)
# 1.2.7 Minor improvement & compatible with "##PRB_MQ=xxx" in QUERY/PEPTIDE_MODIFICATION.REF_POS_STRING for MaxQuant (PP 15/02/17)
# 1.2.6 Minor modif for N-Term modifications (GA 29/10/13)
# 1.2.5 Minor change in $sthPTMpos syntax to retrieve all PTMs (GA 17/09/13)
# 1.2.4 Minor bug correction for ambiguous positions (GA 02/09/13)
# 1.2.3 New ambigous handling (GA 28/06/13)
# 1.2.2 Modification because of bad strings in Sequenceview call and revertion not functional (GA 14/06/13)
# 1.2.1 Remove VMOD, VAR_MOD and RMOD from script (GA 03/06/13)
# 1.2.0 Minor modification following &convertVarModString update in promsMod (GA 17/04/13)
# 1.1.9 VMOD and RMOD tags don't need to be joined anymore in INFO_PEP string (FY 11/04/13)
# 1.1.8 PhosphoRS results are now correctly read for validaiton-ended analyses (FY 08/04/13)
# 1.1.7 Allows revertion to original PTMs position (PP 04/04/13)
# 1.1.6 Checks if comments is defined (PP 06/02/03)
# 1.1.5 Handles uncertainty on N/C-term position of Acetylation & reassignement to K (PP 25/01/13)
# 1.1.4 Added K,R as possible residues for (n)Methylation (PP 15/01/13)
# 1.1.3 Removed useless getSearchParam call (FY 12/11/12)
# 1.1.2 All S,T or Y residues are now selectable even for Phospho(ST) or Phospho(Y) var mods (FY 08/11/12)
# 1.1.1 Displaying message if no PRS isoform for current peptide (status=4) (FY 06/11/12)
# 1.1.0 Fixed missing comma and parsing rules for COM attribute in INFO_PEP (FY 20/07/12)
# 1.0.9 Improved display quality (PP 11/06/12)
# 1.0.8 Check if PEPTIDE.DATA is not null before PRS string parsing (PP 07/06/12)
# 1.0.7 Add score display in PRS results table (FY 25/05/12)
# 1.0.6 Add comments field for unknown position of modifications (GA 24/05/12)
# 1.0.5 PhosphoRS result view (FY 26/03/12)
# 1.0.4 changed '\s*\+\s*' to ' \+ ' when processing varModStrg (PP 28/07/11)
# 1.0.3 management of PMF format (PP 10/05/11)