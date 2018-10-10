#!/usr/local/bin/perl -w

################################################################################
# login.cgi      2.1.1                                                         #
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

####################################################
# Script allowing login to myProMS
# Written by Patrick Poullet
####################################################
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use DBI;
use DBD::mysql;
use POSIX qw(strftime); # to get the time
use promsConfig;
use promsMod;

#######################
####>Configuration<####
#######################
&skipLogin if $ENV{'REMOTE_USER'}; # already authentified by web server
my %promsPath=&promsConfig::getServerInfo('no_user');
my $badConnection=0;
my $oldUser='';
$oldUser=&connect2Server if param('user');
my $chkFullWin=(param('user') && !param('fullWin'))? '' : 'checked';

#######################
####>Starting HTML<####
#######################
print header;
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>Login to myProMS Server</TITLE>
<LINK rel="icon" type="image/png" href="$promsPath{images}/myProMS_icon.png">
<LINK rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
<SCRIPT LANGUAGE="JavaScript">
function checkWindow() {
	if (!opener) document.getElementById('fullWinBox').disabled=true; // in case redirection from bad cookie
}
function checkForm(myForm) {
	if (myForm.user.value && myForm.password.value) {
		return true;
	}
	alert('Please, enter your login and password.');
	return false;
}
</SCRIPT>
</HEAD>
<BODY onload="checkWindow()">
<CENTER>
<FORM method="post" onsubmit="return(checkForm(this));">
<FONT class="title">Login to myProMS Server</FONT><BR>
<TABLE cellspacing=0 cellpadding=3>
<TR><TH align=right class="title2">User:</TH><TD><INPUT type="text" name="user" value="$oldUser" size="15" /></TD></TR>
<TR><TH align=right class="title2">Password:</TH><TD><INPUT type="password" name="password" size="15" /></TD></TR>
<TR><TD class="title2" colspan=2><INPUT type="checkbox" name="fullWin" id="fullWinBox" value="1" $chkFullWin/> Run in dedicated window.</TD></TR>
<TR><TH colspan=2><INPUT type="submit" name="login" value=" Login " class="title2" />&nbsp&nbsp&nbsp<INPUT type="button" value=" Cancel " onclick="window.close()" /></TH></TR>
|;
if ($badConnection) {
	print '<TR><TH colspan=2><FONT style="font-weight:bold;color:#C50000">';
	if ($badConnection==2) {print 'Wrong password. Try again.';}
	elsif ($badConnection==3) {print 'Your account was closed. Contact your mass spec facility.';}
	else {print 'User unknown. Try again.';}
	print "</FONT></TH></TR>\n";
}
print qq
|</TABLE>
</FORM>
</CENTER>
</BODY>
</HTML>
|;


###################################
####<<<Connecting to myProMS>>>####
###################################
sub connect2Server {
	my $userID=param('user');
	my $password=crypt(param('password'),'ProMS');
	my $fullWin=(param('fullWin'))? 1 : 0;

	#####################################
	####<Checking login and password>####
	#####################################
	my $dbh=&promsConfig::dbConnect('no_user');
	my $sthL=$dbh->prepare("SELECT COUNT(*),ID_USER FROM USER_LIST WHERE ID_USER=?"); # retrieves true userID from db because for mySQL 'user_x'='user_x    '!!!!
	$sthL->execute($userID);
	my ($loginOK,$trueUserID)=$sthL->fetchrow_array;
	$sthL->finish;
	my ($userStatus,$closedAccount);
	if ($loginOK) { # user is OK
		$userID=$trueUserID; # replace to avoid keeping 'user_x    '
		my $sthP=$dbh->prepare("SELECT USER_STATUS,IS_CLOSED FROM USER_LIST WHERE ID_USER=? AND USER_PASSWORD=?");
		$sthP->execute($userID,$password);
		($userStatus,$closedAccount)=$sthP->fetchrow_array;
		$sthP->finish;
	}
	if ($userStatus && !$closedAccount) { # password & account are OK
		my $date=strftime("%Y-%m-%d %H:%M:%S",localtime);
		$dbh->do("UPDATE USER_LIST SET LAST_CONNECTION='$date' WHERE ID_USER='$userID'");
		$dbh->commit;
		$dbh->disconnect;

		####<Updating log file>####
		mkdir $promsPath{'logs'} unless -e $promsPath{'logs'};
		my ($year)=($date=~/^(\d+)/);
		open (LOG,">>$promsPath{logs}/connections_$year.log");
		print LOG "$date\t$userStatus\t$userID\n";
		close LOG;

		####<Starting HTML>####
		my $user_cookie=cookie(-name=>'myproms', -value=>"USER=$userID,KEY=".crypt($userID,"myproms_$ENV{REMOTE_ADDR}"), -path=>"$promsPath{cgi}/"); #, -expires=>'+72h'
		print header(-cookie=>$user_cookie);
		print qq
|<HTML>
<HEAD>
<TITLE>Login to myProMS Server</TITLE>
<SCRIPT LANGUAGE="JavaScript">
|;
		if ($fullWin) {
			print qq
|	window.moveTo(0,0);
	window.resizeTo(screen.availWidth,screen.availHeight);
	window.location="$promsPath{html}/proms.html";
|;
		}
		else {
			print qq
|if (opener) {
	opener.location="$promsPath{html}/proms.html";
	window.close();
}
else {window.location="$promsPath{html}/proms.html";}
|;
		}
		print qq
|</SCRIPT>
</HEAD>
</HTML>
|;
		exit;
	}
	else { # Bad login or account is closed
		$dbh->disconnect;
		$badConnection=($closedAccount)? 3 : ($loginOK)? 2 : 1;
		return $userID;
	}
}


################################
####<<<Skiping login form>>>####
################################
sub skipLogin {
	print header;
	print qq
|<HTML>
<HEAD>
<TITLE>Login to myProMS Server</TITLE>
<SCRIPT LANGUAGE="JavaScript">
	window.moveTo(0,0);
	window.resizeTo(screen.availWidth,screen.availHeight);
	window.location="./selectProject.cgi";
</SCRIPT>
</HEAD>
</HTML>
|;
	exit;
}

####>Revision history<####
# 2.1.1 Handles closed account & prevents full window option to be always checked (PP 11/06/18)
# 2.1.0 [Fix] bug due to mySQL considering 'user_x' equal to 'user_x   '!!! (PP 08/06/18)
# 2.0.9 Uses mkdir instead of make_path (PP 10/03/14)
# 2.0.8 GPL license (PP 23/09/13)
# 2.0.7 Records connection in file &lt;logs dir&gt;/connections_&lt;year&gt;.log (PP 21/02/13)
# 2.0.6 Finishing DBI statement handles before dbh disconnecting (FY 12/04/12)
# 2.0.5 used of prepare for critical SQL queries (PP 16/06/11)
