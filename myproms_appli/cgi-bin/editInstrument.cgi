#!/usr/local/bin/perl -w

################################################################################
# editInstrument.cgi         1.1.8                                             #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
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
$| = 1;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;
use strict;

print header(-charset=>'UTF-8'); warningsToBrowser(1);
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my $userID=$ENV{'REMOTE_USER'};
my %fragmentsClassif=&promsConfig::getFragmentClassif;

###############################
####>Connecting to myProMS<####
###############################
my $dbh=&promsConfig::dbConnect;

#############################
####>Fetching parameters<####
#############################
my $action=param('ACT');
my $instrumentID=param('ID');
my $save=(param('save'))?1:0;


#############################
####>Special run in case of delete<####
#############################
if ($action eq 'delete') {
	$dbh->do("DELETE FROM INSTRUMENT WHERE ID_INSTRUMENT=$instrumentID");
	$dbh->commit();
	$dbh->disconnect();
	print header(-charset=>'UTF-8'); warningsToBrowser(1);
	print qq
|<HTML>
<BODY onLoad="window.location='./listInstruments.cgi'">
</BODY>
</HTML>
|;
	exit;
}

########generating fragment List
my %fragmentsList ;
foreach my $fragment (@{$fragmentsClassif{"N_term"}}) {
	$fragmentsList{$fragment}=0;
	foreach my $neutral (@{$fragmentsClassif{"neutral_loss"}}) {
		$fragmentsList{$fragment.$neutral}=0;
	}
}
foreach my $fragment (@{$fragmentsClassif{"C_term"}}) {
	$fragmentsList{$fragment}=0;
	foreach my $neutral (@{$fragmentsClassif{"neutral_loss"}}) {
		$fragmentsList{$fragment.$neutral}=0;
	}
}
foreach my $fragment (@{$fragmentsClassif{"intern"}}) {
	$fragmentsList{$fragment}=0;
	foreach my $neutral (@{$fragmentsClassif{"neutral_loss"}}) {
		$fragmentsList{$fragment.$neutral}=0;
	}
}
$fragmentsList{'+'}=0;
$fragmentsList{'2+'}=0;
$fragmentsList{'3+'} = 0;

#######################
####>Fetching data<#### from table INSTRUMENT
#######################
my ($name, $alias, $des, $comment, $usesStatus, $upDate, $upUser, $deltPar, $deltFrag, $tol, $nr, $rules) ;
my ($pageTitle, $alreadyInUse);
my %deltFragSelt = (
	'Da'=>''
);

my %deltParentSelt = (
	'Da'=>'',
	'ppm'=>''
);

my $errorRules= ' ';
my $saveAsNewStr =' ';



#####################################################
############ Saving results        #############################
#####################################################
if ($save==1) {
	$name = param ('name');
	$alias = param ('alias')?param ('alias'):'';
	$des = (param ('des'))?param ('des'):'';
	$comment= (param ('comments'))?param ('comments'):'';
	$usesStatus = param ('useStatus') ;
	$upDate = param ('date') ;
	#$upUser = param ('user') ;
	$deltPar = param ('deltPart') ;
	$deltFrag = param ('deltFrag') ;
	$tol = param ('tol') ;
	$nr = param ('nr') ;
	#$rules = param ('rules') ;
	$upUser=$userID;

	foreach my $fragment (keys %fragmentsList){
		my $fragmentsetting=param("$fragment");
		$fragmentsList{$fragment}=$fragmentsetting;
		$rules.=$fragment.'='.$fragmentsetting.',' if $fragmentsetting > 0;

	}
	chop $rules ;

	####Formatting value####
	$name=&promsMod::resize($name,50);
	$des=&promsMod::resize($des,50);
	$comment=&promsMod::resize($comment,50);
	$rules =~s/\s//g; #removing space, \n etc...

	####Fchecking for error####
	# my @rulesArray = split (/,/, $rules) ;
	# foreach my $line (@rulesArray) {
		# if (($line !~ /\w.?=\d/)&& ($line !~ /\d?[+-]=\d/)){
			# $errorRules.="$line,\n" ;
			# $errorData = 1;
		# }
	# }
	if ($usesStatus eq 'yes') {
		my $localName= $dbh->quote($name);
		if ($action eq 'add') {
			$alreadyInUse = $dbh->selectrow_array ("SELECT COUNT(*)FROM INSTRUMENT where NAME=$localName and USE_STATUS='yes'");
		}
		elsif 	($action eq 'edit'){
			$alreadyInUse = $dbh->selectrow_array ("SELECT COUNT(*)FROM INSTRUMENT where NAME=$localName and USE_STATUS='yes' and ID_INSTRUMENT!=$instrumentID");
		}
	}
	if ($usesStatus eq 'usr') {
		my $localName= $dbh->quote($name);
		my $localUser = $dbh->quote($userID);
		if ($action eq 'add') {
			$alreadyInUse = $dbh->selectrow_array ("SELECT COUNT(*)FROM INSTRUMENT where NAME=$localName and USE_STATUS='usr' and update_user=$localUser");
		}
		elsif 	($action eq 'edit'){
			$alreadyInUse = $dbh->selectrow_array ("SELECT COUNT(*)FROM INSTRUMENT where NAME=$localName and USE_STATUS='usr' and update_user=$localUser and ID_INSTRUMENT!=$instrumentID");
		}
	}

	if ($alreadyInUse==0){
		$name=$dbh->quote($name);
		$alias=$dbh->quote($alias);
		$des=$dbh->quote($des);
		$comment=$dbh->quote($comment);
		$rules=$dbh->quote($rules);
		$usesStatus=$dbh->quote($usesStatus);
		$upUser=$dbh->quote($upUser);
		$deltPar=$dbh->quote($deltPar);
		$deltFrag=$dbh->quote($deltFrag);
		$upDate=$dbh->quote($upDate);
		#print ("$ID, $name,$des,$comment,$usesStatus,$upDate,$upUser,$deltPar,$deltFrag,$tol,$nr,$rules");#DEBUG
		if ($action eq 'add') {
			my $ID = $dbh->selectrow_array ("SELECT MAX(ID_INSTRUMENT) FROM INSTRUMENT");
			$ID++ ;
			$dbh->do ("INSERT INTO INSTRUMENT (ID_INSTRUMENT,NAME,ALIAS,DES,COMMENTS,USE_STATUS,UPDATE_DATE,UPDATE_USER,DELTA_PARENT,DELTA_FRAGMENT,TOL_FRAGMENT,NR_LEVEL,RULES) VALUES($ID, $name,$alias,$des,$comment,$usesStatus,$upDate,$upUser,$deltPar,$deltFrag,$tol,$nr,$rules)") || die ();
		}
		elsif 	($action eq 'edit'){
			my $updateString = "NAME=$name,ALIAS=$alias,DES=$des,COMMENTS=$comment,USE_STATUS=$usesStatus,UPDATE_DATE=$upDate,UPDATE_USER=$upUser,";
			$updateString .= "DELTA_PARENT=$deltPar,DELTA_FRAGMENT=$deltFrag,TOL_FRAGMENT=$tol,NR_LEVEL=$nr,RULES=$rules";
			$dbh->do("UPDATE INSTRUMENT SET $updateString WHERE ID_INSTRUMENT=$instrumentID") || die $dbh->errstr;
		}
		$dbh->commit();
		$dbh->disconnect();
		#exit ;#DEBUG

		print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
</HEAD>
<BODY background='$promsPath{images}/bgProMS.gif'
<CENTER>
<BR><BR><IMG src='$promsPath{images}/engrenage.gif'>
</CENTER>
<FONT class=\"title2\">Saving Instrument definition ... </FONT>
|;
		sleep 3;
		print qq
|<SCRIPT LANGUAGE="JavaScript">
window.location ="./listInstruments.cgi";
</SCRIPT>
</BODY>
</HTML>
|;
		exit;
	}

	elsif ($alreadyInUse  > 0 ) {
		$pageTitle="Instrument already in Use";
	}
}

else { # save==0
	if ($action eq 'edit') {
		$saveAsNewStr ='&nbsp&nbsp<INPUT type="submit" name="save" value=" Save as new " onClick="saveAsNew();">';
		###>Fetching instrument info
		($name, $alias, $des, $comment, $usesStatus, $upDate, $upUser, $deltPar, $deltFrag, $tol, $nr, $rules)=
			$dbh->selectrow_array("SELECT NAME,ALIAS,DES, COMMENTS,USE_STATUS, UPDATE_DATE, UPDATE_USER,  DELTA_PARENT, DELTA_FRAGMENT,TOL_FRAGMENT, NR_LEVEL, RULES FROM INSTRUMENT WHERE ID_INSTRUMENT=$instrumentID");

		$pageTitle="Editing Instrument <FONT color='#DD0000'>$name</FONT>";
		$deltFragSelt{$deltFrag}='selected';
		$deltParentSelt{$deltPar}='selected';

		my @rulesArray = split (/,/, $rules) ;
		foreach my $line (@rulesArray) {
			my ($fragmentType,$value)=split (/=/,$line) ;
			$fragmentsList{$fragmentType}=$value;
		}
	}

	elsif ($action eq 'add') { # add

		$name= '' ;
		$alias= '';
		$des= '';
		$comment = '';
		$usesStatus = 'yes';
		$upDate = strftime("%Y-%m-%d %H:%M:%S", localtime);
		$upUser =$userID ;
		$deltPar = '' ;
		$deltFrag = '';
		$tol ='';
		$nr ='';
		$rules ='';

		####>Formatting HTML code<####
		$pageTitle='Adding new Instrument';
	}
}

$dbh->disconnect;

#######################
####>Starting HTML<####
#######################

print qq
|<HTML>
<HEAD>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
function saveAsNew() {
	document.InstForm.ACT.value='add';
}
function checkForm(myForm) {
	if (!myForm.name.value) {
		alert('Type a name for this instrument.');
		return false;
	}
	if (!myForm.tol.value) {
		alert('Type a tolerance for this instrument.');
		return false;
	}
	if (!myForm.nr.value) {
		alert('Type a noise reduction level for this instrument.');
		return false;
	}
	return true;
}
</SCRIPT>
</HEAD>
|;

print qq
|<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">$pageTitle</FONT>
<BR>
|;

####>FORM<####
print '<FORM name="InstForm" method="post" onsubmit="return checkForm(this);">';
print hidden(-name=>'ACT',-default=>"$action"),"\n";
print hidden(-name=>'ID',-default=>"$instrumentID"),"\n" if defined($instrumentID) ;

####>TABLE
my $yesSelected = ($usesStatus eq 'yes')?'selected':'';
my $noSelected= ($usesStatus eq 'no')?'selected':'';
my $usrSelected= ($usesStatus eq 'usr')?'selected':'';
my $resetString=($action eq 'add')? '  Clear  ' : ' Cancel changes ';
my $displayDate=&promsMod::formatDate($upDate);
my ($light,$dark)=&promsConfig::getRowColors;

print qq
|<TABLE border=0>
<TR><TD bgcolor=$dark><TABLE border=0 cellpadding=3>
	<TR>
		<TH align=right valign=top bgcolor=$dark >Name :</TH>
		<TD nowrap bgcolor=$light><INPUT type='text' name='name'size=50 value='$name'></TD>
	</TR>
	<TR>
		<TH align=right valign=top bgcolor=$dark >Alias :</TH>
		<TD nowrap bgcolor=$light><INPUT type='text' name='alias'size=50 value='$alias'> (for ProteinPilot software)</TD>
	</TR>
	<TR>
		<TH align=right valign=top bgcolor=$dark >Description :</TH>
		<TD nowrap bgcolor=$light><TEXTAREA name='des' rows=2 cols=78>$des</TEXTAREA></TD>
	</TR>
	<TR>
		<TH align=right valign=top bgcolor=$dark >Fragment tolerance :</TH>
		<TD nowrap bgcolor=$light ><INPUT type='text' name='tol' value='$tol' size=10> <SELECT name='deltFrag'>|;
foreach my $keys (keys %deltFragSelt) {
	print "<OPTION $deltFragSelt{$keys}>$keys</OPTION>\n";
}
print qq
|		</SELECT>
		</TD>
	</TR>
	<TR>
		<TH align=right valign=top bgcolor=$dark nowrap>Parent delta unit :</TH>
		<TD nowrap bgcolor=$light> <SELECT name='deltPart'>
|;
foreach my $key (keys %deltParentSelt) {
	print "<OPTION $deltParentSelt{$key}>$key</OPTION>\n";
}
print qq
|		</SELECT>
		</TD>
	</TR>
	<TR>
		<TH align=right valign=top bgcolor=$dark>Fragments detection :</TH>
		<TD nowrap bgcolor=$light>
		<TABLE><TR><TD align=right><B>Charges :</B></TD>
|;
my $i=0;
foreach my $fragment (sort{$a cmp $b} keys %fragmentsList) {
	$i++;
	my $selected0=($fragmentsList{$fragment} == 0)?'selected':'';
	my $selected1=($fragmentsList{$fragment} == 1)?'selected':'';
	my $selected2=($fragmentsList{$fragment} == 2)?'selected':'';

	my $displayedName = ($fragment eq '3+')? "&ge;3+" : $fragment;

	print qq
|<TD align=right>&nbsp;&nbsp;&nbsp;$displayedName:<SELECT name="$fragment">
		<OPTION value=0 $selected0>unused</OPTION>
		<OPTION value=1 $selected1>table only</OPTION>
		<OPTION value=2 $selected2>spectrum</OPTION>
	</SELECT>
	</TD>
|;
	if ($fragment eq '3+') {
		print "<TD></TD>";
		$i=5;
	}
	if ($i==5) {
		$i=0;
		print "</TR><TR>";
	}
}
print '</TR>' if $i>0;
#<TEXTAREA name='rules' rows='4' cols='65'>$rules</TEXTAREA>

print qq
|		</TABLE>
		</TD>
	</TR>
	<TR>
		<TH align=right valign=top bgcolor=$dark nowrap>Noise reduction level :</TH>
		<TD nowrap bgcolor=$light><INPUT type='text' name='nr' value='$nr' size=10> Minimal intensity (fraction of highest peak)</TD>
	</TR>
	<TR>
		<TH align=right valign=top bgcolor=$dark >Update user :</TH>
		<TD nowrap bgcolor=$light>$upUser<INPUT type=hidden name='user' value='$upUser'></TD>
	</TR>
	<TR>
		<TH align=right valign=top bgcolor=$dark >Date :</TH>
		<TD nowrap bgcolor=$light>$displayDate<INPUT type=hidden name='date' value='$upDate'></TD>
	</TR>
	<TR>
		<TH align=right valign=top bgcolor=$dark >Use status :</TH>
		<TD nowrap bgcolor=$light>
			<SELECT name='useStatus'>
				<OPTION value="yes" $yesSelected>yes</OPTION>
				<OPTION value="no" $noSelected>no</OPTION>
				<OPTION value="usr" $usrSelected>by '$upUser' only</OPTION>
			</SELECT>
		</TD>
	</TR>
	<TR>
		<TH align=right valign=top bgcolor=$dark >Comments :</TH>
		<TD nowrap bgcolor=$light><TEXTAREA name='comments' rows=5 cols=78>$comment</TEXTAREA></TD>
	</TR>

	<TR>
		<TH colspan=2>
		<INPUT type="submit" name="save" value=" Save ">
		$saveAsNewStr
		&nbsp&nbsp<INPUT type="reset" value="$resetString">
		&nbsp&nbsp<INPUT type="button" value="Cancel" onclick="window.location='./listInstruments.cgi'">
		</TH>
	</TR>
</TABLE>
</TD></TR>
</TABLE>
</FORM>
</CENTER>
</BODY>
</HTML>
|;

####>Revision history<####
# 1.1.8 charset=>'UTF-8' (PP 07/04/14)
# 1.1.7 GPL license (PP 23/09/13)
# 1.1.6 Managing >= 3+ charges (FY 18/03/13)
# 1.1.5 Minor modif to show-edit ALIAS for Paragon searches (GA 02/10/12)