#!/usr/local/bin/perl -w

################################################################################
# editUserInfo.cgi             2.2.0                                           #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Displays myProMS main entry page with links to different sections            #
# Called after login to server                                                 #
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
use promsConfig;
use promsMod;
use strict;

#print header; warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;
my ($color1,$color2)=&promsConfig::getRowColors;

#############################
####>Fetching parameters<####
#############################
my $userID=$ENV{'REMOTE_USER'};
my $action=(param('ACT'))? param('ACT') : 'view';
my $selUserID=(param('user'))? param('user') : ($action eq 'new')? '' : $userID;
my $pwdOK=param('pwdOK'); # previous password change was successful

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

####>Form submission<####
if (param('save')) { # submit button was clicked!!!!!
	&saveData;
	exit;
}
elsif ($action=~/close|open/) {
	my $isClosed=($action eq 'close')? 1 : undef;
	my $sthUA=$dbh->prepare("UPDATE USER_LIST SET IS_CLOSED=? WHERE ID_USER=?");
	$sthUA->execute($isClosed,$selUserID);
	$sthUA->finish;
	$dbh->commit;
	$action='view';
}

##############
####>MAIN<####
##############
my $qUserID=$dbh->quote($userID);
my ($userName,$userStatus,$userWorkgroup)=$dbh->selectrow_array("SELECT USER_NAME,USER_STATUS,WORK_GROUP FROM USER_LIST WHERE ID_USER=$qUserID");
$userWorkgroup='' unless $userWorkgroup;
my @statusList=('bio');
if ($userStatus ne 'bio') {
	push @statusList,'manag';
	if ($userStatus ne 'manag') {
		push @statusList,'mass';
		push @statusList,'bioinfo' if $userStatus eq 'bioinfo';
	}
}
my (%userData,%userProjects);
my ($statusAlias,$selUserName,$selUserStatus,$userLab,$userTel,$userEmail,$userInfo,$lastConnection,$userPref,$mascotIDs,$workgroup,$isClosed)=('','','','','','','','Never','','','',0);
if ($action eq 'new') {
	$workgroup=$userWorkgroup if $userStatus eq 'manag'; # users created by manager belong to his Workgroup
}
else {
	my $qSelUserID=$dbh->quote($selUserID);
	($selUserName,$selUserStatus,$userLab,$userTel,$userEmail,$userInfo,$lastConnection,$userPref,$mascotIDs,$workgroup,$isClosed)=$dbh->selectrow_array("SELECT USER_NAME,USER_STATUS,LABORATORY,USER_TEL,USER_EMAIL,USER_INFO,LAST_CONNECTION,USER_PREF,MASCOT_IDS,WORK_GROUP,IS_CLOSED FROM USER_LIST WHERE ID_USER=$qSelUserID");
	$statusAlias=&promsMod::getStatusAlias($selUserStatus);
	$userLab='' unless $userLab;
	$userTel='' unless $userTel;
	$userEmail='' unless $userEmail;
	$userInfo='' unless $userInfo;
	$isClosed=0 unless $isClosed;
	$workgroup='' if ($selUserStatus=~/bioinfo|mass/ || !$workgroup); # delete wGr if new status
	$mascotIDs='' if ($selUserStatus=~/bioinfo|mass/ || !defined $mascotIDs); # delete IDs if new status
	if ($lastConnection) {$lastConnection=&promsMod::formatDate($lastConnection);}
	else {$lastConnection='Unknown';}

	if ($userStatus ne 'bio') {

		####>Fetching list of accessible projects<####
		if ($action eq 'view' && $userID ne $selUserID && $selUserStatus=~/(bio|manag)\Z/) {
			my $wgFilter=($selUserStatus eq 'manag' && $workgroup)? " AND WORK_GROUP != '$workgroup'" : '';
			my $sthPA=$dbh->prepare("SELECT P.ID_PROJECT,P.NAME,DES,WORK_GROUP,UP.NAME FROM PROJECT P,PROJECT_ACCESS PA,USER_PROFILE UP WHERE P.ID_PROJECT=PA.ID_PROJECT AND PA.ID_PROFILE=UP.ID_PROFILE AND ID_USER=$qSelUserID$wgFilter");
			$sthPA->execute;
			while (my ($projID,$name,$des,$wkGr,$accessStatus)=$sthPA->fetchrow_array) {
				next if ($userStatus eq 'manag' && (!$workgroup || !$wkGr || $wkGr ne $workgroup)); # project is not in manager's workgroup => do not list
				@{$userProjects{$projID}}=($name,$wkGr,$accessStatus,$des);
			}
			$sthPA->finish;
		}

		####>Fetching list of users<####
		if ($userStatus eq 'manag' && !$userWorkgroup) { # special case: manager has no workgroup
			@{$userData{$userWorkgroup}{$userStatus}{$userID}}=($userName,$isClosed); # manager himself (there should be no other accessible users if no Wkg)
		}
		else {
			my $whereStrg=($userStatus eq 'manag')? "WHERE WORK_GROUP='$workgroup'" : '';
			my $sthU=$dbh->prepare("SELECT ID_USER,USER_NAME,USER_STATUS,WORK_GROUP,IS_CLOSED FROM USER_LIST $whereStrg");
			$sthU->execute;
			while (my ($uID,$name,$status,$wkGr,$isClos)=$sthU->fetchrow_array) {
				$wkGr='' unless $wkGr;
				$wkGr='ZZZ' if $status=~/mass|bioinfo/;
				@{$userData{$wkGr}{$status}{$uID}}=($name,$isClos);
			}
			$sthU->finish;
		}
	}
}
my %workgroupList;
if ($action=~/edit|new/ && $selUserID ne $userID && $userStatus ne 'manag') {
	#Users
	my $sthWG1=$dbh->prepare('SELECT DISTINCT(WORK_GROUP) FROM USER_LIST WHERE WORK_GROUP IS NOT NULL');
	$sthWG1->execute;
	while (my ($wkGr)=$sthWG1->fetchrow_array) {
		next unless $wkGr;
		$workgroupList{$wkGr}=1;
	}
	$sthWG1->finish;
	#Projects
	my $sthWG2=$dbh->prepare('SELECT DISTINCT(WORK_GROUP) FROM PROJECT WHERE WORK_GROUP IS NOT NULL');
	$sthWG2->execute;
	while (my ($wkGr)=$sthWG2->fetchrow_array) {
		next unless $wkGr;
		$workgroupList{$wkGr}=1;
	}
	$sthWG2->finish;
}

$dbh->disconnect;


####>UserPref Parsing<####
my $specAppCheck=($userPref=~/specApp=1/ || $action eq 'new')? 'checked' : '';
my $vertLabelCheck=($userPref=~/vertLab=1/)? 'checked' : '';
my $seqDbCheck=($userPref=~/seqDb=1/)? 'checked' : '';
my $spLibCheck=($userPref=~/spLib=1/)? 'checked' : '';
my $seqModCheck=($userPref=~/seqMod=1/)? 'checked' : '';
my $goCheck=($userPref=~/go=1/)? 'checked' : '';



####>HTML<####
print header(-charset=>'UTF-8');
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>User Info</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TD {background-color:$color1}
</STYLE>
<SCRIPT LANGUAGE="JavaScript">
function selectUser(action) {
	var selUser=(document.getElementById('user'))? document.getElementById('user').value : '$userID';
	window.location="./editUserInfo.cgi?ACT="+action+"&user="+selUser;
}
function updateWorkgroup(wkGr) {
	var visStatus=(wkGr=='#new#')? 'visible' : 'hidden';
	document.getElementById('newWorkgroup').style.visibility=visStatus;
}
|;
if ($userStatus eq 'bioinfo') {
	print qq
|function closeAccount() {
	if (confirm("Close '$selUserID' account?")) {
		window.location='./editUserInfo.cgi?ACT=close&user=$selUserID';
	}
}
function openAccount() {
	if (confirm("Open '$selUserID' account?")) {
		window.location='./editUserInfo.cgi?ACT=open&user=$selUserID';
	}
}
|;
}
print qq
|function checkForm(myForm) {
	// password
	if (action != 'edit' && (!myForm.pwd1.value \|\| !myForm.pwd2.value \|\| myForm.pwd1.value != myForm.pwd2.value)) {
		alert('Incorrect password. Please try again.');
		return false;
	}
	// login
	if (action == 'new' && !myForm.user.value) {
		alert('Enter a valid login.');
		return false;
	}
	// name
	if (action != 'pwd' && !myForm.name.value) {
		alert('Enter your full name.');
		return false;
	}
	// email
	if (action != 'pwd' && (!myForm.email.value \|\| !myForm.email.value.match(/.+@.+/))) {
		alert('Enter a valid e-mail adress.');
		return false;
	}
	//workgroup (only for bio & manag)
	if (myForm.workgroup && myForm.workgroup.value=='#new#' && !myForm.newWorkgroup.value) {
		alert('Enter a valid Workgroup name.');
		return false;
	}
	return true;
}
var action='$action';
function openProject(projID) {
	window.location="$promsPath{cgi}/openProject.cgi?&ACT=open&ID="+projID+"&branchID=project:"+projID;
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<DIV style="float:top">
<BR>
<TABLE><TR><TH bgcolor="$color2">
<FONT class="title2">&nbsp;Go to:</FONT><SELECT class="title action" onchange="if (this.value) {window.location='$promsPath{cgi}/'+this.value;}">
	<OPTION value="">-= Select =-</OPTION>
	<OPTION value="promsMain.cgi">Main Window</OPTION>
	<OPTION value="selectProject.cgi">Project Selection</OPTION>
</SELECT>
</TH></TR></TABLE>
</DIV>
<CENTER>
<FONT class="title1">User Information:</FONT><BR>
|;
if ($action ne 'view') {
	print "<FORM name=\"userForm\" onsubmit=\"return(checkForm(this));\">\n";
	print "<INPUT type=\"hidden\" name=\"ACT\" value=\"$action\" />\n";
	print "<INPUT type=\"hidden\" name=\"user\" value=\"$selUserID\" />\n" if $action ne 'new';
}
if ($action=~/view|pwd/) {
	print "<TABLE bgcolor=$color2>\n";
	print "<TR><TH align=right width=150><FONT class=\"title2\">Login :</FONT></TH><TD width=600><FONT class=\"title2\">$selUserID</FONT>";
	if ($action eq 'view' && ($selUserID eq $userID || $userStatus ne 'bio') && !$isClosed) {
		print "\n&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" value=\"Change Password\" onclick=\"window.location='./editUserInfo.cgi?ACT=pwd&user=$selUserID'\" />\n";
	}
	print "\n&nbsp;&nbsp;&nbsp;<FONT class=\"title3\">*** Account closed ***</FONT>\n" if $isClosed;
	if ($action eq 'view' && $userStatus eq 'bioinfo' && $selUserID ne $userID) {
		if ($isClosed) {print "\n&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" value=\"Open Account\" onclick=\"openAccount()\" />\n";}
		else {print "\n&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" value=\"Close Account\" onclick=\"closeAccount()\" />\n";}
	}
	print "</TD></TR>\n";
	if ($action eq 'pwd') {
		print "<TR><TH align=right><FONT class=\"title2\">New password :</FONT></TH><TD><INPUT type=\"password\" name=\"pwd1\" value=\"\" size=\"15\" />\n";
		print "&nbsp;&nbsp;<B>Confirm :</B><INPUT type=\"password\" name=\"pwd2\" value=\"\" size=\"15\" /></TD></TR>\n";
	}
	print "<TR><TH align=right><FONT class=\"title2\">Status :</FONT></TH><TD><FONT class=\"title2\">$statusAlias</FONT></TD></TR>\n";
	if ($selUserStatus=~/(manag|bio)\Z/) {
		# Workgroup
		print "<TR><TH align=right><FONT class=\"title3\">Workgroup :</FONT></TH><TD><FONT class=\"title3\">";
		if ($workgroup) {print $workgroup;} else {print 'None';}
		print "</FONT></TD></TR>\n";
		# Mascot ids
		print "<TR><TH align=right><FONT class=\"title3\">Mascot IDs :</FONT></TH><TD>";
		my ($mascotUserID,@mascotGroupIDs)=split(',',$mascotIDs);
		my $userIdStrg=(defined $mascotUserID && $mascotUserID ne '')? $mascotUserID : 'None';
		my $groupIdStrg=(scalar @mascotGroupIDs)? join(',',@mascotGroupIDs) : 'None';
		print "<B>User:</B> $userIdStrg&nbsp;&nbsp;&nbsp;&nbsp;<B>Group(s):</B> $groupIdStrg\n";
		print "</TD></TR>\n";
	}
	if ($action eq 'view' && $selUserStatus=~/manag|mass/) {
		print "<TR><TH align=right valign=top>Annotation access :</TH><TD>";
		my $someAccess=0;
		if ($seqDbCheck) {
			print '&bull; Sequence databanks & Species<BR>';
			$someAccess++;
		}
		if ($spLibCheck) {
			print '&bull; Spectral libraries & Species<BR>';
			$someAccess++;
		}
		if ($seqModCheck) {
			print '&bull; Sequence modifications & Species<BR>';
			$someAccess++;
		}
		if ($goCheck) {
			print '&bull; Gene Ontology & Species';
			$someAccess++;
		}
		print "None" unless $someAccess;
		print "</TD></TR>\n";
	}
	print qq
|<TR><TH align=right>User name :</TH><TD>$selUserName</TD></TR>
<TR><TH align=right>Laboratory :</TH><TD>$userLab</TD></TR>
<TR><TH align=right>Telephone :</TH><TD>$userTel</TD></TR>
<TR><TH align=right>E-mail :</TH><TD>$userEmail</TD></TR>
<TR><TH align=right>Other information :</TH><TD>$userInfo</TD></TR>
|;
	if ($action eq 'pwd') { # change password
		print "<TR><TH colspan=2><INPUT type=\"submit\" name=\"save\" value=\" Save \" />\n";
		print "&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" value=\" Cancel \" onclick=\"window.location='./editUserInfo.cgi?ACT=view'\" />\n";
		print "</TH></TR>\n";
	}
	else { # view
		print "<TR><TH align=right>Last connection :</TH><TD>$lastConnection</TD></TR>\n";
		if ($userStatus ne 'bio') {
			print "<TR><TH align=right valign=top>Other user :</TH><TD nowrap><SELECT id=\"user\" style=\"font-weight:bold;\">\n";
			#my @statusList=('bio','manag');
			#push @statusList,'mass' if $userStatus ne 'manag';
			#push @statusList,'bioinfo' if $userStatus eq 'bioinfo';
			foreach my $wkGr (sort{lc($a) cmp lc($b)} keys %userData) {
				my ($shiftStrg,$bulletStrg)=('','');
				if ($userStatus ne 'manag') {
					$shiftStrg='&nbsp;&nbsp;&nbsp;&nbsp;';
					$bulletStrg='-';
					my $wkGrLabel=($wkGr)? "Wg '$wkGr'" : 'No Workgroup';
					print "<OPTGROUP label=\"+$wkGrLabel:\"></OPTGROUP>\n" if $wkGr ne 'ZZZ';
				}
				foreach my $status (@statusList) {
					next unless $userData{$wkGr}{$status};
					my $usedStrg=($status=~/bioinfo|mass/)? '+' : $shiftStrg.$bulletStrg;
					print "<OPTGROUP label=\"$usedStrg",&promsMod::getStatusAlias($status),"s:\">\n";
					foreach my $user (sort{lc($userData{$wkGr}{$status}{$a}[0]) cmp lc($userData{$wkGr}{$status}{$b}[0])} keys %{$userData{$wkGr}{$status}}) {
						print '<OPTION ';
						print 'selected ' if $selUserID eq $user;
						print "value='$user'>$shiftStrg$userData{$wkGr}{$status}{$user}[0]";
						print ' [closed]' if $userData{$wkGr}{$status}{$user}[1];
						print "</OPTION>\n";
					}
					print "</OPTGROUP>\n";
				}
			}
			print "</SELECT>\n";
			print "<INPUT type=\"button\" value=\" View \" onclick=\"selectUser('view')\"/>\n";
			if ($userStatus ne 'bio') {
				print "<INPUT type=\"button\" value=\" Edit \" onclick=\"selectUser('edit')\"/>\n";
				print "&nbsp;&nbsp;&nbsp;&nbsp;<INPUT type=\"button\" value=\"Add New User\" onclick=\"window.location='./editUserInfo.cgi?ACT=new'\">\n";
			}
			print "</TD></TR>\n";
		}
# 		else {
# 			print "<TR><TH colspan=2><INPUT type=\"button\" value=\" Edit \" onclick=\"selectUser('edit')\"/></TH></TR>\n";
# 		}
	}
	print "</TABLE>\n";
}
else { # edit or new
	print "<TABLE bgcolor=$color2>\n";
	print "<TR><TH align=right width=150><FONT class=\"title2\">Login :</FONT></TH><TD width=500>";
	if ($action eq 'edit') {print "<FONT class=\"title2\">$selUserID</FONT></TD></TR>\n";}
	else { # new
		print "<INPUT type=\"text\" name=\"user\" value=\"\" size=\"20\" style=\"font-size:16px;font-weight:bold\" /></TD></TR>\n";
		print "<TR><TH align=right><FONT class=\"title2\">Password :</FONT></TH><TD><INPUT type=\"password\" name=\"pwd1\" value=\"\" size=\"15\" />\n";
		print "&nbsp;&nbsp;<B>Confirm :</B><INPUT type=\"password\" name=\"pwd2\" value=\"\" size=\"15\"/></TD></TR>\n";
	}
	print "<TR><TH align=right><FONT class=\"title2\">Status :</FONT></TH><TD><SELECT name=\"userStat\" style=\"font-size:16px;font-weight:bold\">\n";
	foreach my $status (@statusList) {
		# next if ($status eq 'bioinfo' && $userStatus ne 'bioinfo'); # only bioinfo can promote bioinfo
		# next if ($status ne 'bio' && $userStatus eq 'bio'); # bio cannot self-promote
		print "<OPTION value=\"$status\"";
		print ' selected' if $status eq $selUserStatus;
		print '>',&promsMod::getStatusAlias($status),"</OPTION>\n";
	}
	print "</SELECT></TD></TR>\n";
	if (!$selUserStatus || $selUserStatus=~/manag|bio\Z/) {
		# Workgroup
		print "<TR><TH align=right><FONT class=\"title3\">Workgroup :</FONT></TH><TD>";
		if ($userStatus=~/(manag|bio)\Z/) {
			print "<FONT class=\"title3\">";
			if ($workgroup) {print $workgroup;} else {print 'None';}
			print "</FONT><INPUT type=\"hidden\" name=\"workgroup\" value=\"$workgroup\"/>";
		}
		else {
			print "<SELECT name=\"workgroup\" class=\"title3\" onchange=\"updateWorkgroup(this.value)\"><OPTION value=''>None</OPTION>\n";
			foreach my $wkGr (sort{lc($a) cmp lc($b)} keys %workgroupList) {
				print "<OPTION value=\"$wkGr\"";
				print ' selected' if $wkGr eq $workgroup;
				print ">$wkGr</OPTION>\n";
			}
			print qq
|<OPTION value="#new#">-= New Workgroup =-</OPTION>
</SELECT>&nbsp;<INPUT type="text" name="newWorkgroup" id="newWorkgroup" value="" style="visibility:hidden"/>&nbsp;
|;
		}
		print "</TD></TR>\n";
		# Mascot ids
		print "<TR><TH align=right><FONT class=\"title3\">Mascot IDs :</FONT></TH><TD>";
		my ($mascotUserID,@mascotGroupIDs)=split(',',$mascotIDs);
		if ($userStatus=~/manag|bio\Z/) {
			my $userIdStrg=(defined $mascotUserID)? $mascotUserID : 'None';
			my $groupIdStrg=(scalar @mascotGroupIDs)? join(',',@mascotGroupIDs) : 'None';
			print "<B>User:</B> $userIdStrg&nbsp;&nbsp;&nbsp;&nbsp;<B>Group(s):</B> $groupIdStrg\n";
		}
		else {
			my $userIdStrg=(defined $mascotUserID)? $mascotUserID : '';
			my $groupIdStrg=(scalar @mascotGroupIDs)? join(',',@mascotGroupIDs) : '';
			foreach my $g (0..2) {
				$mascotGroupIDs[$g]='' unless $mascotGroupIDs[$g];
			}
			print qq
|<B>User:</B><INPUT type="text" name="mascotUserID" value="$userIdStrg" size="3">
&nbsp;&nbsp;&nbsp;&nbsp;<B>Group(s):</B>
<INPUT type="text" name="mascotGroupID" value="$groupIdStrg" size="12">&nbsp;<FONT class="font11"><I>(comma-separated list)</I></FONT>
|;
		}
		print "</TD></TR>\n";
	}

	if ($selUserStatus=~/manag|mass/) {
		print "<TR><TH align=right valign=top>Annotation access :</TH><TD>\n";
		if ($userStatus eq 'bioinfo') { # editable
			print qq
|<INPUT type="checkbox" name="seqDb" value="1" $seqDbCheck>Sequence databanks & Species<BR>
<INPUT type="checkbox" name="spLib" value="1" $spLibCheck>Spectral libraries & Species<BR>
<INPUT type="checkbox" name="seqMod" value="1" $seqModCheck>Sequence modifications<BR>
<INPUT type="checkbox" name="go" value="1" $goCheck>Gene Ontology & Species
|;
		}
		else { # not editable
			my $someAccess=0;
			if ($seqDbCheck) {
				print '&bull; Sequence databanks & Species<INPUT type="hidden" name="seqDb" value="1"><BR>';
				$someAccess++;
			}
			if ($spLibCheck) {
				print '&bull; Spectral libraries & Species<INPUT type="hidden" name="specLb" value="1"><BR>';
				$someAccess++;
			}
			if ($seqModCheck) {
				print '&bull; Sequence modifications<INPUT type="hidden" name="seqMod" value="1"><BR>';
				$someAccess++;
			}
			if ($goCheck) {
				print '&bull; Gene Ontology & Species<INPUT type="hidden" name="go" value="1">';
				$someAccess++;
			}
			print "None" unless $someAccess;
		}
		print "</TD></TR>\n";
	}

	print qq
|<TR><TH align=right>User name :</TH><TD><INPUT type="text" name="name" value="$selUserName" size="40"></TD></TR>
<TR><TH align=right>Laboratory :</TH><TD><INPUT type="text" name="lab" value="$userLab" size="65"></TD></TR>
<TR><TH align=right>Telephone :</TH><TD><INPUT type="text" name="tel" value="$userTel" size="20"></TD></TR>
<TR><TH align=right>E-mail :</TH><TD><INPUT type="text" name="email" value="$userEmail" size="40"></TD></TR>
<TR><TH align=right>Other information :</TH><TD><INPUT type="text" name="info" value="$userInfo" size="65"></TD></TR>
<TR><TH align=right valign=top>Preferences :</TH><TD>
	<INPUT type="checkbox" name="specApp" $specAppCheck>Use interactive spectrum<BR>
	<INPUT type="checkbox" name="vertLabel" $vertLabelCheck>Set label vertical
</TD></TR>
<TR><TH colspan=2>
	<INPUT type="submit" name="save" value=" Save " />
	&nbsp;&nbsp;&nbsp;<INPUT type="reset" value="Clear Changes" />
	&nbsp;&nbsp;&nbsp;<INPUT type="button" value=" Cancel " onclick="window.location='./editUserInfo.cgi?ACT=view'" />
</TH></TR>
</TABLE>
|;
}
print "</FORM>\n" if $action=~/new|edit|pwd/;

if ($action eq 'view' && $userID ne $selUserID && $selUserStatus=~/(bio|manag)\Z/) {
	print qq
|<BR><BR><FONT class="title1">Projects accessible by <FONT color="#DD0000">$selUserID</FONT>:</FONT><BR>
<TABLE cellspacing=0 cellpadding=4>
<TR bgcolor="$color2"><TH class="rbBorder title3">&nbsp;Name&nbsp;</TH><TH class="rbBorder title3">&nbsp;Workgroup&nbsp;</TH><TH class="rbBorder title3">&nbsp;Access right&nbsp;</TH><TH class="bBorder title3">&nbsp;Description&nbsp;</TH></TR>
|;
	my $bgColor=$color1;
	my $count=0;
	if ($selUserStatus eq 'manag' && $workgroup) {
		print "<TR><TD colspan=4>&nbsp;<FONT class=\"title3\">All projects in Workgroup <FONT color=\"#DD0000\">$workgroup</FONT>.<FONT></TD></TR>\n";
		$bgColor=$color2;
		$count++;
	}
	foreach my $projID (sort{lc($userProjects{$a}[0]) cmp lc($userProjects{$b}[0])} keys %userProjects) {
		$count++;
		my $wkGr=$userProjects{$projID}[1] || 'No Workgroup';
		my $accStatus=&promsMod::getProfileAlias($userProjects{$projID}[2]);
		my $des=$userProjects{$projID}[3] || '';
		print "<TR><TD style=\"background-color:$bgColor\">&nbsp;<INPUT type=\"button\" value=\"Open\" onclick=\"openProject($projID)\">&nbsp;<FONT class=\"title3\">$userProjects{$projID}[0]</FONT>&nbsp;</TD><TD style=\"background-color:$bgColor\">&nbsp;$wkGr&nbsp;</TD><TD style=\"background-color:$bgColor\">&nbsp;$accStatus&nbsp;</TD><TD style=\"background-color:$bgColor\">$des</TD></TR>\n";
		$bgColor=($bgColor eq $color1)? $color2 : $color1;
	}
	print "<TR><TD colspan=4>&nbsp;<FONT class=\"title3\">No projects found.</FONT></TD></TR>\n" unless $count;
}

if ($pwdOK) {
	print "<SCRIPT LANGUAGE=\"JavaScript\">alert('Password change was successful.');\n</SCRIPT>\n";
}

print "</CENTER></BODY>\n</HTML>\n";

#######################<<< SUBROUTINE >>>###########################

####################################
####< Storing data in DataBase >####
####################################
sub saveData {

	####<Password>####
	my $password;
	if ($action=~/pwd|new/) {
		$password=crypt(param('pwd1'),'ProMS');
		$password=$dbh->quote($password);
	}

	####<Change password>####
	if ($action eq 'pwd') {
		$dbh->do("UPDATE USER_LIST SET USER_PASSWORD=$password WHERE ID_USER='$selUserID'");
	}

	####<Edit or New>####
	else {
		my $selUserStatus=param('userStat');
		my $mascotIDs='';
		if ($selUserStatus=~/manag|bio/) {
			$mascotIDs=param('mascotUserID') if defined param('mascotUserID');
			$mascotIDs=~s/\D//g; # keep digits only;
			if (defined param('mascotGroupID')) {
				foreach my $grID (split(',',param('mascotGroupID'))) {
					$grID=~s/\D//g; # keep digits only;
					$mascotIDs.=",$grID" if $grID=~/\d/;
				}
			}
		}
		my $annotAccessStrg='';
		if ($selUserStatus=~/manag|mass/) {
			$annotAccessStrg.=(param('seqDb'))? ';seqDb=1' : ';seqDb=0';
			$annotAccessStrg.=(param('spLib'))? ';spLib=1' : ';spLib=0';
			$annotAccessStrg.=(param('seqMod'))? ';seqMod=1' : ';seqMod=0';
			$annotAccessStrg.=(param('go'))? ';go=1' : ';go=0';
		}
		$mascotIDs=$dbh->quote($mascotIDs);
		my $userName=param('name'); $userName=$dbh->quote($userName);
		my $userLab=param('lab'); $userLab=$dbh->quote($userLab);
		my $userTel=param('tel'); $userTel=$dbh->quote($userTel);
		my $userEmail=param('email'); $userEmail=$dbh->quote($userEmail);
		my $userInfo=param('info'); $userInfo=$dbh->quote($userInfo);
		my $specApp=(param('specApp'))?1:0;
		my $vertLabel=(param('vertLabel'))?1:0;
		my $userPref = "specApp=$specApp;vertLab=$vertLabel$annotAccessStrg"; $userPref=$dbh->quote($userPref);
		my $workgroup=(!param('workgroup'))? undef : (param('workgroup') eq '#new#')? param('newWorkgroup') : param('workgroup');
		$workgroup=$dbh->quote($workgroup);
		if ($action eq 'new') {
			##<Checking if user already exists
			my ($exist)=$dbh->selectrow_array("SELECT COUNT(*) FROM USER_LIST WHERE ID_USER='$selUserID'");
			if ($exist) {
				print header(-charset=>'UTF-8');
				print qq
|<HTML>
<HEAD>
<SCRIPT LANGUAGE="JavaScript">
alert("ERROR: User '$selUserID' already exists in database!");
window.location="./editUserInfo.cgi?ACT=view&user=$selUserID";
</SCRIPT></HTML>
|;
				exit;
			}
			else {
				$dbh->do("INSERT INTO USER_LIST (ID_USER,USER_PASSWORD,USER_STATUS,USER_NAME,LABORATORY,USER_TEL,USER_EMAIL,USER_INFO,USER_PREF,MASCOT_IDS,WORK_GROUP,START_DATE) VALUES ('$selUserID',$password,'$selUserStatus',$userName,$userLab,$userTel,$userEmail,$userInfo,$userPref,$mascotIDs,$workgroup,NOW())");
			}
		}
		else { # edit
			$dbh->do("UPDATE USER_LIST SET USER_NAME=$userName,USER_STATUS='$selUserStatus',LABORATORY=$userLab,USER_TEL=$userTel,USER_EMAIL=$userEmail,USER_INFO=$userInfo,USER_PREF=$userPref,MASCOT_IDS=$mascotIDs,WORK_GROUP=$workgroup WHERE ID_USER='$selUserID'");
			if ($selUserStatus=~/bioinfo|mass/) { # in case of promotion
				$dbh->do("DELETE FROM PROJECT_ACCESS WHERE ID_USER='$selUserID'"); # => delete all his project accesses
			}
			elsif ($selUserStatus eq 'manag') { # in case of promotion
				$dbh->do("DELETE FROM PROJECT_ACCESS WHERE ID_USER='$selUserID' AND ID_PROJECT IN (SELECT ID_PROJECT FROM PROJECT WHERE WORK_GROUP=$workgroup)");
			}
		}
	}
	#exit ; #debug
	$dbh->commit;
	$dbh->disconnect;

	###<Refresh window
	if ($action eq 'pwd') {
		print redirect("./editUserInfo.cgi?ACT=view&user=$selUserID&pwdOK=1");
	}
	else {
		print redirect("./editUserInfo.cgi?ACT=view&user=$selUserID");
	}
	exit;
}

####>Revision history<####
# 2.2.0 Added account close/open management (PP 02/10/18)
# 2.1.7 Added access control on Spectral libraries (PP 13/11/17)
# 2.1.6 Defaults to interactive spectrum for new user (PP 20/10/16)
# 2.1.5 Records account creation date in START_DATE field (PP 12/05/16)
# 2.1.4 Added list of projects accessible to user (PP 22/03/16)
# 2.1.3 Added Sequence modification access control (PP 24/07/13)
# 2.1.2 Called by promsMain.cgi (PP 13/12/12)
# 2.1.1 Minor bug fix (PP 18/09/12)
# 2.1.0 Adding restrictions to manager & massist on annotation files (PP 07/05/12)
# 2.0.9 Bug correction when upgrading to Manager (PP 17/10/11)
# 2.0.8 Allows Mascot groupID=0 (PP 10/10/11)
# 2.0.7 Manager management (PP 07/09/11)
# 2.0.6 Allows no Mascot userID but only group IDs (PP 23/08/11)
# 2.0.5 Management of MASCOT Ids (PP 14/03/2011)
