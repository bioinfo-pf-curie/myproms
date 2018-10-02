#!/usr/local/bin/perl -w

################################################################################
# projectAccess.cgi         2.0.6                                              #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Sets project accessibility to biologists                                     #
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

#######################
####>Configuration<####
#######################
my %promsPath=&promsConfig::getServerInfo;

##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;
my $userID=$ENV{'REMOTE_USER'};

###############################
####> Fetching parameters <####
###############################
my $projectID=param('ID');
&processForm if param('save');
my @userInfo=&promsMod::getUserInfo($dbh,$userID,$projectID);
my $projectAccess=$userInfo[2]->{$projectID};
my $editRight=($projectAccess=~/(bioinfo|mass|manag|administrator)/)? 1 : 0;
my ($projectName,$projWorkgroup) = $dbh->selectrow_array("SELECT NAME,WORK_GROUP FROM PROJECT WHERE ID_PROJECT=$projectID");
($projWorkgroup,my $projWorkGrStrg)=($projWorkgroup)? ($projWorkgroup,"Workgroup: $projWorkgroup") : ('','No Workgroup assigned');

##############################################################
####>Fetching all biologists and access rights to project<####
##############################################################
my %accessOptions; # @{...{key}}=(ID_PROFILE, display order, accessibility to current user)
my $sthAcID=$dbh->prepare("SELECT ID_PROFILE,NAME FROM USER_PROFILE");
$sthAcID->execute;
while (my ($acID,$acName)=$sthAcID->fetchrow_array) {
	@{$accessOptions{"$acName"}}=($acID); # [0] -> id
}
$sthAcID->finish;
push @{$accessOptions{'guest'}},(1,1); push @{$accessOptions{'user'}},(2,1); # [1] -> position [2] -> excludable from project
push @{$accessOptions{'power_user'}},(3,0); push @{$accessOptions{'super_user'}},(4,0);
push @{$accessOptions{'administrator'}},(5,1);
push @{$accessOptions{'power_administrator'}},(6,0); push @{$accessOptions{'super_administrator'}},(7,0);
if ($projectAccess ne 'bio') { # ecludable by manag, mass & bioinfo
	$accessOptions{'power_user'}[2]=$accessOptions{'power_administrator'}[2]=1;
	$accessOptions{'super_user'}[2]=$accessOptions{'super_administrator'}[2]=1;
}

####>Fetching list of managers & biologists<####
my (%nameList,%rightList,%userProjectAccess,%noAccessUsers);
my $sthUser=$dbh->prepare("SELECT ID_USER,USER_NAME,USER_STATUS,WORK_GROUP FROM USER_LIST WHERE USER_STATUS='bio' OR USER_STATUS='manag'");
$sthUser->execute;
while (my ($uID,$userName,$userStatus,$workgroup)=$sthUser->fetchrow_array) {
	my @uInfo = &promsMod::getUserInfo($dbh,$uID,$projectID);
	my $statusName=&promsMod::getStatusAlias($userStatus); # eq 'manag')? 'Manager' : 'Biologist';
	$workgroup='' unless $workgroup;
	if ($userStatus eq 'manag' && $workgroup && $workgroup eq $projWorkgroup) { # manager access
		$userProjectAccess{$uID}=$userStatus;
	}
	elsif ($uInfo[2]->{$projectID}) { # "bio"-type access
		$userProjectAccess{$uID}=$uInfo[2]->{$projectID};
	}
	else { # no access
		push @{$noAccessUsers{$workgroup}{$statusName.'s'}},$uID;
	}
	my $wkGrName=($workgroup)? $workgroup : 'None';
	push @{$nameList{$uID}},($userName,$statusName,$wkGrName);
}
$sthUser->finish;

$dbh->disconnect;

#######################
####>Starting HTML<####
#######################
my ($light,$dark)=&promsConfig::getRowColors;
my $bgColor=$light;

print header(-'content-encoding'=>'no',-charset=>'utf-8');
warningsToBrowser(1);

print qq
|<HTML>
<TITLE>Project Access</TITLE>
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<STYLE type="text/css">
TD {font-size:13px;}
</STYLE>
|;
if ($editRight) {
	print qq
|<SCRIPT LANGUAGE="JavaScript">
function showUserAccount(uID) {
	top.promsFrame.location="$promsPath{cgi}/editUserInfo.cgi?ACT=view&user="+uID;
}
function activateButtons() {
	document.accessForm.save.disabled=false;
	document.accessForm.cancel.disabled=false;
}
function cancelChanges() {
	window.location="./projectAccess.cgi?ID=$projectID";
}
</SCRIPT>
|;
}
print qq
|</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<FONT class="title">Accessibility to Project <FONT color="#DD0000">$projectName</FONT></FONT><FONT class="title2"><BR>$projWorkGrStrg.</FONT>
<BR><BR>
<FORM name="accessForm" method="post">
<TABLE border=0 cellspacing=0>
<TR><TH colspan=4><FONT style="font-size:18px">Users allowed to access this Project</FONT></TH></TR>
<TR bgcolor="$dark">
	<TH class="rbBorder">&nbsp;User&nbsp;</TH>
	<TH class="rbBorder">&nbsp;Status&nbsp;</TH>
	<TH class="rbBorder">&nbsp;Workgroup&nbsp;</TH>
	<TH class="bBorder">&nbsp;Access Right<SUP>*</SUP>&nbsp;</TH>
</TR>
|;
my $userString='';
my $numAccess=0;
foreach my $uID (sort{lc($nameList{$a}[0]) cmp lc($nameList{$b}[0])} keys %userProjectAccess) { #%nameList
	#next if $rightList{$uID} eq 'no access'; # do not list users with no access
	#next if $orphanManagers{$uID}; # skip orphan managers
	$numAccess++;
	print qq
|<TR bgcolor="$bgColor">
	<TH nowrap align=left>&nbsp;$nameList{$uID}[0] [$uID]&nbsp;</TH>
	<TH>&nbsp;$nameList{$uID}[1]&nbsp;</TH>
	<TH nowrap>&nbsp;$nameList{$uID}[2]&nbsp;</TH>
|;
	print "<TD nowrap>&nbsp;";
	if ($editRight && $userProjectAccess{$uID} ne 'manag') {
		print "<SELECT name=\"newAccess_$uID\" style=\"width:170px\" onchange=\"javascript:activateButtons()\">\n";
		if ($accessOptions{$userProjectAccess{$uID}} && $accessOptions{$userProjectAccess{$uID}}[2]==0) { # main user is not bioinfo nor mass nor manag
# 			print " disabled>\n";
			print "<OPTION value=''>",&promsMod::getProfileAlias($userProjectAccess{$uID}),"</OPTION>\n";
			print "<OPTION value='0'>No access</OPTION>\n"; # any administrator can reject anyone from project
		}
		else {
# 			print ">\n";
			print "<OPTION value='0'>No access</OPTION>\n" unless $uID eq $userID; # no self rejection
			foreach my $access (sort{$accessOptions{$a}[1]<=>$accessOptions{$b}[1]} keys %accessOptions) {
				next if $accessOptions{$access}[2]==0;
				print '<OPTION ';
				if ($userProjectAccess{$uID} eq $access) {print 'selected ';}
				print "value='$accessOptions{$access}[0]'>",&promsMod::getProfileAlias($access),"</OPTION>\n";
			}
			$userString.='::' if $userString;
			$userString.=$uID;
		}
		print "</SELECT>";
	}
	else {print &promsMod::getProfileAlias($userProjectAccess{$uID});}
	print "&nbsp;<INPUT type=\"button\" onclick=\"showUserAccount('$uID')\" value=\"Account\">" if $editRight;
	print "&nbsp;</TD></TR>\n";
	$bgColor=($bgColor eq $light)? $dark : $light;
}
unless ($numAccess) { # no listed user
	print "<TR bgcolor=\"$bgColor\"><TH colspan=4 align=left>&nbsp;No users</TH></TR>"
}
if ($editRight) {
	print qq
|<TR><TH colspan=4><BR><FONT style="font-size:18px;font-weight:bold;">Allow <SELECT name="newUser" style="font-size:15px;font-weight:bold;" onchange="javascript:activateButtons()">
<OPTION value="">-=Choose from List=-</OPTION>
|;
	foreach my $workgroup (sort{$a cmp $b} keys %noAccessUsers) {
		my $wkGrLabel=($workgroup)? "Workgroup '$workgroup'" : 'No Workgroup';
		print "<OPTGROUP label=\"+$wkGrLabel:\"></OPTGROUP>\n";
		foreach my $statusName (sort{$a cmp $b} keys %{$noAccessUsers{$workgroup}}) {
			print "<OPTGROUP label=\"&nbsp;&nbsp;&nbsp;&nbsp;-$statusName:\">";
			foreach my $uID (sort{lc($nameList{$a}[0]) cmp lc($nameList{$b}[0])} @{$noAccessUsers{$workgroup}{$statusName}}) {
				#next unless $rightList{$uID} eq 'no access';
				print "<OPTION value='$uID:$accessOptions{guest}[0]'>&nbsp;&nbsp;&nbsp;&nbsp;$nameList{$uID}[0] [$uID]</OPTION>\n";
			}
			print "</OPTGROUP>\n";
		}

	}
	print qq
|</SELECT> to access this Project.</FONT><BR><BR></TH></TR>
<TR><TD colspan=4 align=center>
<INPUT type="submit" name="save" value="Save" style="width:70px" disabled/>
&nbsp&nbsp&nbsp<INPUT type="button" name="cancel" value="Cancel Changes" style="width:150px" onclick="javascript:cancelChanges()" disabled/>
&nbsp&nbsp&nbsp<INPUT type="button" style="width:70px" value="End" onclick="top.promsFrame.optionFrame.selectOption()"/>
</TD></TR>
</TABLE>
<INPUT type="hidden" name="ID" value="$projectID"/>
<INPUT type="hidden" name="userList" value="$userString"/>
|;
	foreach my $uID (sort keys %userProjectAccess) {
		next if $userProjectAccess{$uID} eq 'manag';
		print "<INPUT type=\"hidden\" name=\"oldAccess_$uID\" value=\"$accessOptions{$userProjectAccess{$uID}}[0]\"/>\n";
	}
}
else {print "</TABLE>\n";}
print qq
|</FORM><BR>
<TABLE cellspacing=0 bgcolor=$dark>
<TR><TD><TABLE>
<TR><TD colspan=2 class="title2">&nbsp;Access rights description:</TD></TR>
<TR><TH align=right>Guest :</TH><TD bgcolor=$light>Read access to validated data.</TD></TR>
<TR><TH align=right>User :</TH><TD bgcolor=$light>Read/Write access to validated data.</TD></TR>
<TR><TH align=right>Administrator :</TH><TD bgcolor=$light><B>User</B> + Project access management.</TD></TR>
<TR><TH align=right>Power (User/Administrator) :</TH><TD bgcolor=$light><B>User/Admin.</B> with additional read/write access to non-validated protein data.</TD></TR>
<TR><TH align=right>Super (User/Administrator) :</TH><TD bgcolor=$light><B>User/Admin.</B> with full access rights on the current project.</TD></TR>
<TR><TH align=right>Manager :</TH><TD bgcolor=$light>Full access rights on all projects of a workgroup.</TD></TR>
</TABLE></TD></TR>
</TABLE>
</BODY>
</HTML>
|;


###################################
####>Processing Submitted form<####
###################################
sub processForm {
#print header,"\n";warningsToBrowser(1);

	####<Updating access to project>####
	my @userList=split('::',param('userList'));
	foreach my $uID (@userList) {
		my ($oldAccessID,$newAccessID)=(param("oldAccess_$uID"),param("newAccess_$uID"));
		##<No change
		next if $newAccessID == $oldAccessID;
		##<Access right is being deleted
		if ($newAccessID == 0) {
			$dbh->do("DELETE FROM PROJECT_ACCESS WHERE ID_USER='$uID' AND ID_PROJECT=$projectID");
		}
		##<Access right is being modified (admin<=>user)
		else {
			$dbh->do("UPDATE PROJECT_ACCESS SET ID_PROFILE=$newAccessID,UPDATE_DATE=NOW(),UPDATE_USER='$userID' WHERE ID_USER='$uID' AND ID_PROJECT=$projectID");
		}
	}

	####<Granting access to another biologist>####
	if (param('newUser')) {
		my ($uID,$accessID)=split(/:/,param('newUser'));
		$dbh->do("INSERT INTO PROJECT_ACCESS (ID_USER,ID_PROJECT,ID_PROFILE,UPDATE_DATE,UPDATE_USER) VALUES ('$uID',$projectID,$accessID,NOW(),'$userID')");
	}

	$dbh->commit;
}

####>Revision history<####
# 2.0.6 udpate UPDATE_DATE and UPDATE_USER fields & direct link to user account (PP 22/03/16)
# 2.0.5 GPL license (PP 23/09/13)
# 2.0.4 Manager are selectable if no workGroup assigned (PP 17/03/12)
# 2.0.3 Manager management (PP 18/09/11)
# 2.0.2 Now accessible to all users (but guest) but edition is restricted as before (PP 13/04/11)
# 2.0.1 Extension of Super user/admin access rights (PP 14/03/2011)
