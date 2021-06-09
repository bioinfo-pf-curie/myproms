#!/usr/local/bin/perl -w

################################################################################
# promsMain.cgi         1.1.8                                                  #
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

#print header,"DEBUG<BR>\n"; warningsToBrowser(1); # DEBUG
#######################
####>Configuration<####
#######################
my ($lightColor,$darkColor)=&promsConfig::getRowColors;
my %promsPath=&promsConfig::getServerInfo;
my $start=(param('START'))? 1 : 0;


##########################
####>Connecting to DB<####
##########################
my $dbh=&promsConfig::dbConnect;

###############################
####>Fetching user details<####
###############################
my $userID=$ENV{'REMOTE_USER'};
#my ($userName,$userStatus,$refProfile,$userLab,$userTel,$userEmail,$userInfo,$refMascotIDs,$workgroup)=&promsMod::getUserInfo($dbh,$userID);
my ($userName,$userStatus,$lastConnection,$userPref)=$dbh->selectrow_array("SELECT USER_NAME,USER_STATUS,LAST_CONNECTION,USER_PREF FROM USER_LIST WHERE ID_USER='$userID'");
$dbh->disconnect;
my $fullStatus=&promsMod::getStatusAlias($userStatus);
$lastConnection=&promsMod::formatDate($lastConnection);
#my $superBio=($userStatus eq 'bio' && scalar(@{$refMascotIDs}))? 1 : 0;

if ($userStatus eq 'bio' && $start) {
	print header(-charset=>'utf-8');
	warningsToBrowser(1);
	print qq
|<HTML>
<HEAD>
<TITLE>myProMS Main Page</TITLE>
<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT LANGUAGE="JavaScript">
top.logoFrame.document.getElementById('userImg').src="$promsPath{images}/user_logged.gif";
top.logoFrame.document.getElementById('userInfo').innerHTML="$userName ($fullStatus), logged since $lastConnection.";
window.location="$promsPath{cgi}/selectProject.cgi";
</SCRIPT>
</HEAD>
</HTML>
|;
	exit;
}


my %colors=(0=>'#909090',1=>'#0000A5');
my (%hoverAction,%clickAction);

#>Users
my ($userString,$userDx)=($userStatus eq 'bio')? ('Account',200) : ('Users',160); #(1==1)?

#>Annotation data
my $okSeqDBs=($userStatus eq 'bioinfo' || $userPref=~/seqDb=1/)? 1 : 0;
($hoverAction{'seqDBs'},$clickAction{'seqDBs'})=($okSeqDBs==1)? ('on','selectDatabank.cgi') : ('silent','null');
my $okSpLib=($userStatus eq 'bioinfo' || $userPref=~/spLib=1/)? 1 : 0;
($hoverAction{'specLibs'},$clickAction{'specLibs'})=($okSpLib==1)? ('on','listSwathLibraries.cgi') : ('silent','null');
my $okSeqMod=($userStatus eq 'bioinfo' || $userPref=~/seqMod=1/)? 1 : 0;
($hoverAction{'seqMod'},$clickAction{'seqMod'})=($okSeqMod)? ('on','manageModifications.cgi') : ('silent','null');
my $okGoSp=($userStatus eq 'bioinfo' || $userPref=~/go=1/)? 1 : 0;
($hoverAction{'goAnnot'},$clickAction{'goAnnot'})=($okGoSp==1)? ('on','manageGOFiles.cgi') : ('silent','null');
my $okGSets=($userStatus eq 'bioinfo' || $userPref=~/gSet=1/)? 1 : 0;
($hoverAction{'geneSets'},$clickAction{'geneSets'})=($okGSets==1)? ('on','manageGeneSets.cgi') : ('silent','null');
my $okSpecies=($okSeqDBs || $okGoSp || $okSpLib || $okGSets)? 1 : 0;
($hoverAction{'species'},$clickAction{'species'})=($okSpecies)? ('on','manageSpecies.cgi') : ('silent','null');
my $okAnnotData=($okSeqDBs || $okSeqMod || $okGoSp || $okSpLib || $okGSets)? 1 : 0;
@{$hoverAction{'annotData'}}=($okAnnotData)? ('on','wait') : ('none','none');

#>Settings
my $okValidTpl=($userStatus eq 'bio')? 0 : 1;
($hoverAction{'validTpl'},$clickAction{'validTpl'})=($okValidTpl==1)? ('on','manageTemplates.cgi?ACT=list') : ('silent','null');
my $okRefPeptides=($userStatus eq 'bio')? 0 : 1;
($hoverAction{'refPeptides'},$clickAction{'refPeptides'})=($okRefPeptides==1)? ('on','manageReferencePeptides.cgi?ACT=list') : ('silent','null');

my $okLabelReagt=($userStatus eq 'bio')? 0 : 1;
($hoverAction{'labelReagents'},$clickAction{'labelReagents'})=($okLabelReagt==1)? ('on','manageLabelReagents.cgi?ACT=list') : ('silent','null');

my $okInstr=($userStatus eq 'bio')? 0 : 1;
($hoverAction{'instr'},$clickAction{'instr'})=($okInstr==1)? ('on','listInstruments.cgi') : ('silent','null');
my $okSettings=($okValidTpl || $okRefPeptides || $okLabelReagt || $okInstr)? 1 : 0;
@{$hoverAction{'settings'}}=($okSettings)? ('on','wait') : ('none','none');

#>Data mining
my $okDataMin=($userStatus=~/bioinfo|mass/)? 1 : 0;
($hoverAction{'proteins'},$clickAction{'proteins'})=($okDataMin==1)? ('on','proteinMining.cgi') : ('silent','null');
@{$hoverAction{'dataMin'}}=($okDataMin)? ('on','wait') : ('none','none');

#>Server
my $okLogStats=($userStatus=~/bioinfo|mass/)? 0 : 0;
($hoverAction{'logs'},$clickAction{'logs'})=($okLogStats==1)? ('on','????') : ('silent','null');
($hoverAction{'stats'},$clickAction{'stats'})=($okLogStats==1)? ('on','????') : ('silent','null');
my $okTest=($userStatus eq 'bioinfo')? 1 : 0;
($hoverAction{'test'},$clickAction{'test'})=($okTest==1)? ('on','testMyProMS.cgi?CALL=server') : ('silent','null');
my $okServer=($okLogStats || $okTest)? 1 : 0;
@{$hoverAction{'server'}}=($okServer)? ('on','wait') : ('none','none');

##############
####>HTML<####
##############
print header(-charset=>'utf-8');
warningsToBrowser(1);
print qq
|<HTML>
<HEAD>
<TITLE>myProMS Main Page</TITLE>
<SCRIPT src="$promsPath{html}/js/Raphael/raphael.js"></SCRIPT>
<SCRIPT LANGUAGE="JavaScript">
top.logoFrame.document.getElementById('userImg').src="$promsPath{images}/user_logged.gif";
top.logoFrame.document.getElementById('userInfo').innerHTML="$userName ($fullStatus), logged since $lastConnection.";
//Clear all tree status
//top.projectTreeStatus='';
//top.itemTreeStatus='';
window.onload=function() {
	var menu = Raphael("canvas",1000,700);
	//menu.image("$promsPath{images}/massSpectra2.png",50,100,900,561); // x0.75
	//menu.rect(0,0,1000,700,0).attr({stroke:'black'});
	menu.rect(0,0,1000,7,0).attr({stroke:'none',fill:'$colors{1}'});
	menu.rect(0,693,1000,7,0).attr({stroke:'none',fill:'$colors{1}'});
	//baseline
	var baseH=650;
	menu.path('M60 '+baseH+' l880 0').attr({stroke:'#f00','stroke-width':8,'stroke-linecap':'round'}); // baseline
	//Variables
	var defaultColor;

	// Users
	var uPeakX=200,uPeakH=200; // 140
	addPeak(uPeakX-40,20,1);
	addPeak(uPeakX-6,40,1);
	addPeak(uPeakX,uPeakH,4); // main
	addPeak(uPeakX+50,40,1);
	// Projects
	var pPeakX=280,pPeakH=500;
	addPeak(pPeakX-15,100,1);
	addPeak(pPeakX-7,20,1);
	addPeak(pPeakX,pPeakH,6); // main
	addPeak(pPeakX+7,180,1);
	addPeak(pPeakX+15,80,1);
	// Annotation data
	var adPeakX=370,adPeakH=320;
	addPeak(adPeakX-55,15,1);
	addPeak(adPeakX,adPeakH,3); // main
	addPeak(adPeakX+10,25,2);
	// Settings
	var setPeakX=480,setPeakH=400;
	addPeak(setPeakX-7,60,1);
	addPeak(setPeakX,setPeakH,2); // main
	addPeak(setPeakX+35,12,0);
	//Data mining
	var dmPeakX=580,dmPeakH=240;
	addPeak(dmPeakX-15,25,1);
	addPeak(dmPeakX,dmPeakH,1); // main
	//Server
	var serPeakX=760,serPeakH=180;
	addPeak(serPeakX,serPeakH,1); // main
	addPeak(serPeakX+25,10,0);


	  /********************************************************************/
	 /*                               TEXT                               */
	/********************************************************************/
	var textProperties={'text-anchor':'start','font-family':'verdana,arial,helvetica,sans-serif','font-size':38,'font-weight':'bold'};
	var popupTextProperties={'text-anchor':'start','font-family':'verdana,arial,helvetica,sans-serif','font-size':28,'font-weight':'bold'};
	var pathProperties={'stroke-width':3,'stroke-linecap':'round','stroke-linejoin':'round'};

	//Users-Account
	drawLink(uPeakX,uPeakH,uPeakX-$userDx,330,'$userString','$colors{1}','on','off','editUserInfo.cgi?ACT=view');

	//Projects
	drawLink(pPeakX,pPeakH,pPeakX-200,70,'Projects','$colors{1}','on','off','selectProject.cgi');

	//Server (From here: reverse order to keep popups visible -->)
	drawLink(serPeakX,serPeakH,serPeakX+20,400,'Server','$colors{$okServer}','$hoverAction{server}[0]','$hoverAction{server}[1]',null,
			 [
				['Log files','$colors{$okLogStats}','$hoverAction{logs}','$clickAction{logs}'],
				['Statistics','$colors{$okLogStats}','$hoverAction{stats}','$clickAction{stats}'],
				['Test','$colors{$okTest}','$hoverAction{test}','$clickAction{test}']
			 ]
			);

	//Data mining
	drawLink(dmPeakX,dmPeakH,dmPeakX+35,280,'Data mining','$colors{$okDataMin}','$hoverAction{dataMin}[0]','$hoverAction{dataMin}[1]',null,
			 [
				['Protein identification','$colors{$okDataMin}','$hoverAction{proteins}','$clickAction{proteins}']
			 ]
			);

	//Settings
	drawLink(setPeakX,setPeakH,setPeakX+20,170,'Settings','$colors{$okSettings}','$hoverAction{settings}[0]','$hoverAction{settings}[1]',null,
			 [
				['Validation templates','$colors{$okValidTpl}','$hoverAction{validTpl}','$clickAction{validTpl}'],
				['Reference peptides','$colors{$okRefPeptides}','$hoverAction{refPeptides}','$clickAction{refPeptides}'],
				['Label reagents','$colors{$okLabelReagt}','$hoverAction{labelReagents}','$clickAction{labelReagents}'],
				['Instruments','$colors{$okInstr}','$hoverAction{instr}','$clickAction{instr}']
			 ]
			);

	//Annotation data
	drawLink(adPeakX,adPeakH,adPeakX+75,70,'Annotation data','$colors{$okAnnotData}','$hoverAction{annotData}[0]','$hoverAction{annotData}[1]',null,
			 [
				['Sequence databanks','$colors{$okSeqDBs}','$hoverAction{seqDBs}','$clickAction{seqDBs}'],
				['Spectral libraries','$colors{$okSpLib}','$hoverAction{specLibs}','$clickAction{specLibs}'],
				['Sequence modifications','$colors{$okSeqMod}','$hoverAction{seqMod}','$clickAction{seqMod}'],
				['GO Annotations','$colors{$okGoSp}','$hoverAction{goAnnot}','$clickAction{goAnnot}'],
				['Gene Sets','$colors{$okGSets}','$hoverAction{geneSets}','$clickAction{geneSets}'],
				['Species','$colors{$okSpecies}','$hoverAction{species}','$clickAction{species}']
			 ]
			);

	//User's Guide
	drawSimpleLink(20,675,'start',"User's Guide",displayUsersGuide);

	//Disconnect
	drawSimpleLink(980,675,'end','Close session',top.closeSession);

	  /********************************************************************/
	 /*                            FUNCTIONS                             */
	/********************************************************************/
	function addPeak(posX,height,width) {
		var hw=Math.ceil(width/2);
		var peakPath='M'+(posX-hw)+' '+(baseH-3)+'l'+hw+' -'+height+'l'+hw+' '+height+'Z';
		menu.path(peakPath).attr({stroke:'#f00','stroke-width':3,'stroke-linecap':'round','stroke-linejoin':'round',fill:'#F00'});
	}

	function drawLink(peakX,peakH,startX,startY,label,color,hoverOn,hoverOff,clickAction,popupLinks) {
		/** Master link **/
		var masterLink=menu.text(startX,startY,label).attr(textProperties).attr({fill:color});
		var line=menu.path(drawPath(masterLink,peakX,peakH)).attr(pathProperties).attr({stroke:color}).data('text',masterLink);
		var mCursor=(hoverOn=='on')? 'pointer' : 'default';
		var mBox=masterLink.getBBox(); // for IE compatibility (a transparent box above text drives mouse events)
		var linkBox=menu.rect(mBox.x,mBox.y,mBox.width,mBox.height).attr({fill:'#FFF',opacity:0,stroke:'none',cursor:mCursor})
		.data({text:masterLink,color:color,hoverOn:hoverOn,hoverOff:hoverOff,line:line})
		.hover(function(){highlightLink(this,'hoverOn')},function(){highlightLink(this,'hoverOff')});
		if (clickAction) {
			linkBox.data('clickAction',clickAction)
			.click(function(){window.location='./'+this.data('clickAction');});
		}
		/** Popup box **/
		if (popupLinks) {
			/* Child links */
			var childLinks=[];
			childLinks[0]=null;
			linkBox.data('popup',childLinks);
			var childX=startX+15, childY=startY+50;
			for (var i=0; i<popupLinks.length; i++) {
				var j=i+1;
				//popupLinks[i]=[label,color,hoverOn,clickAction] hoverOff is allways 'wait'
				var subLink=menu.text(childX,childY,'· '+popupLinks[i][0]).attr(popupTextProperties).attr({fill:popupLinks[i][1]});
				var sCursor=(popupLinks[i][2]=='on')? 'pointer' : 'default';
				var sBox=subLink.getBBox();
				childLinks[j]=menu.rect(sBox.x,sBox.y,sBox.width,sBox.height).attr({fill:'#FFF',opacity:0,stroke:'none',cursor:sCursor})
				.data({text:subLink,color:popupLinks[i][1],hoverOn:popupLinks[i][2],parent:linkBox})
				.hover(function(){highlightLink(this,'hoverOn')},function(){highlightLink(this,'wait')});
				if (popupLinks[i][2]=='on') {
					childLinks[j].data('clickAction',popupLinks[i][3]);
					childLinks[j].click(function(){window.location='./'+this.data('clickAction');});
				}
				childY+=35;
			}
			/* Surrounding box */
			var childSet=menu.set(childLinks);
			var setBox=childSet.getBBox();
			var popupBox=menu.rect(setBox.x-10,setBox.y-5,setBox.width+20,setBox.height+10,10)
			.attr({'stroke-width':2,stroke:'#DD0000',fill:'#FFF'}).data('parent',linkBox)
			.hover(function(){highlightLink(this,'silent')},function(){highlightLink(this,'wait')});
			childSet.clear();
			childLinks[0]=popupBox; // 1st index is box
			for (var i=0; i<childLinks.length; i++) {
				if (i>0) childLinks[i].data('text').toFront().hide();
				childLinks[i].toFront().hide();
			}
		}
	}

	function drawSimpleLink(refX,refY,anchor,label,clickFunction) {
		var simpleLink=menu.text(refX,refY,label).attr(popupTextProperties).attr({fill:'$colors{1}','text-anchor':anchor});
		var lBox=simpleLink.getBBox(); // for IE compatibility (a transparent box above text drives mouse events)
		var simpleBox=menu.rect(lBox.x,lBox.y,lBox.width,lBox.height).attr({fill:'#FFF',opacity:0,stroke:'none',cursor:'pointer'})
		.data({text:simpleLink,color:'$colors{1}'})
		.hover(function(){highlightLink(this,'on')},function(){highlightLink(this,'off')})
		.click(function(){clickFunction();});
	}

	function drawPath(text,peakX,peakH) {
		var tBox=text.getBBox();
		var pathStrg;
		if (text.attr('x') < peakX) {
			pathStrg='M'+text.attr('x')+' '+(text.attr('y')+tBox.height/2)+' l'+tBox.width+' 0 L'+(peakX-2)+' '+(baseH-peakH-10);
		}
		else {
			pathStrg='M'+(text.attr('x')+tBox.width)+' '+(text.attr('y')+tBox.height/2)+' l-'+tBox.width+' 0 L'+(peakX+2)+' '+(baseH-peakH-10);
		}
		return pathStrg;
	}

	function highlightLink(linkBox,hoverType) {
		var action=(hoverType.match('hoverO'))? linkBox.data(hoverType) : hoverType;
//console.log(action);
	    if (action=='none') return;
	    var newColor, cursor='default';
	    var isPopup=(linkBox.data('parent'))? true : false;
	    if (action=='on' \|\| action=='silent') {
			if (action=='on') {
				newColor='#DD0000';
				cursor='pointer';
			}
			if (counter) clearInterval(counter);
			var newParent=(isPopup)? linkBox.data('parent') : linkBox;
			if (newParent != selParent) {
				if (selParent) highlightLink(selParent,'off');
				selParent=newParent;
			}
			if (linkBox.data('popup')) { // hover on parent link
				var childLinks=linkBox.data('popup');
				for (var i=0; i<childLinks.length; i++) {
					if (i>0) childLinks[i].data('text').show();
					childLinks[i].show();
				}
			}
	    }
	    else if (action=='wait') {
			visDuration=750; // millisec
			if (counter) clearInterval(counter);
			counter=setInterval(timer,250);
			if (isPopup) newColor=linkBox.data('color');
	    }
	    else { // action=off
			newColor=linkBox.data('color');
			if (linkBox.data('popup')) { // hover on parent link
				var childLinks=linkBox.data('popup');
				for (var i=0; i<childLinks.length; i++) {
					if (i>0) childLinks[i].data('text').hide();
					childLinks[i].hide();
				}
				visDuration=0;
			}
			selParent=null;
	    }
	    if (newColor) {
			//if (!text.data('isBox'))
			linkBox.data('text').attr({fill:newColor});
			linkBox.attr({cursor:cursor});
			if (linkBox.data('line')) {linkBox.data('line').attr({stroke:newColor});}
	    }
	}

	var visDuration=0;
	var counter;
	var selParent;

	function timer() {
	    visDuration=visDuration-250;
	    if (visDuration <= 0) {
			clearInterval(counter);
			//counter ended, do something here
			highlightLink(selParent,'off');
			//selParent.data('popup').hide();
			return;
	    }
	    //Do code for showing the number of seconds here
//console.log(visDuration);
	}

	function displayUsersGuide() {
		//window.open("http://myproms-demo.curie.fr/myProMSUsersGuide.pdf",'UsersGuideWindow');
		window.open("https://myproms-doc.readthedocs.io/en/latest",'UsersGuideWindow');
	}
}
</SCRIPT>
</HEAD>
<BODY background="$promsPath{images}/bgProMS.gif">
<CENTER>
<BR>
<IMG border=0 src="$promsPath{images}/myProMS48.gif">
<DIV id="canvas" style="text-align:center;"></DIV><!--position:absolute; left:0%; top:50px; width:100%; -->
</CENTER>
</BODY>
</HTML>
|;

####>Revision history<####
# 1.1.8 [FEATURE] Add Gene Sets management in Annotation data (VL 10/11/20)
# 1.1.7 Points to Readthedocs User's guide (PP 18/06/19)
# 1.1.6 Add 'Reference peptides' & 'Label reagents links' (PP 30/04/19)
# 1.1.5 Dedicated flag Spectral libraries access (PP 13/11/17)
# 1.1.4 Changed "SWATH libraries" to Spectral libraries" (PP 17/05/17)
# 1.1.3 Changed Data mining "Proteins" to "Protein identification" (PP 20/06/16)
# 1.1.2 Data mining enabled (PP 31/05/16)
# 1.1.1 Add SWATH libraries link. Data mining disabled (PP 18/03/16)
# 1.1.0 add proteins in data_mining part (SL 26/06/15)
# 1.0.9 Link to demo User's guide (PP 22/10/13)
# 1.0.8 Link to User's guide activated (PP 23/09/13)
# 1.0.7 Bug fix in annotation item activation (PP 24/07/13)
# 1.0.6 Minor change in Sequence Modification activation (PP 10/05/13)
# 1.0.5 Activate Sequence Modifications in Annotation Data (GA 23/04/13)
# 1.0.4 Disable User's guide link (PP 23/04/13)
# 1.0.3 Minor change in popup boxes display (PP 22/04/13)
# 1.0.2 Minor bug fix in settings for data mining link (PP 01/03/13)
# 1.0.1 Changed Links name & organisation (PP 21/02/13)
# 1.0.0 First stable version (PP 13/12/12)
# 0.0.1 created (PP 30/04/12)
