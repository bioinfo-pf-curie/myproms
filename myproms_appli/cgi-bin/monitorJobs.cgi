#!/usr/local/bin/perl -w

################################################################################
# monitorJobs.cgi       1.1.16                                                 #
# Authors: V. Sabatet (Institut Curie)                                         #
# Contact: myproms@curie.fr                                                    #
# Monitors all jobs (Quantifications, Libraries, PhosphoRS, DIA, TDA)          #
# called from processAnalyses.cgi                                              #
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
use strict;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use POSIX qw(strftime tzset floor); # to get the time
use File::stat;
use File::Basename;
use File::Path qw(rmtree); # remove_tree
use promsConfig;
use promsMod;
use promsQuantif;
use Time::Piece;
$ENV{'TZ'} = 'Europe/Paris';
tzset();

#######################
####>Configuration<####
#######################
my %promsPath = &promsConfig::getServerInfo;
my $userID = $ENV{'REMOTE_USER'};

my $LIMIT_STORAGE = "1 MONTH"; # Store for a selected amount of time (matching mySQL INTERVAL format)
my $REFRESH_TIME_INTERVAL = 5; # Refresh time interval (seconds)

my %quantifProcesses = (
    'XICMCQ'                     => 'Ext. ion chrom.',
    'EMPAI'                      => 'emPAI',
    'SIN'                        => 'SI<SUB>N</SUB>',
    'XICCORR'                    => 'XIC correction',
    'DESIGN'                     => 'Protein ratio',
    'DESIGN:SimpleRatio'         => 'Protein ratio',
    'DESIGN:SuperRatio'          => 'Super ratio',
    'PROT_RATIO_PEP'             => 'Protein ratio',
    'DESIGN:myProMS'             => 'Protein ratio',
    'DESIGN:PEP_INTENSITY'       => 'Protein ratio (Peptide intensity)',
    'DESIGN:PEP_INTENSITY:Abund' => 'Protein Abundance (DDA)',
    'DESIGN:PEP_RATIO'           => 'Protein ratio (Peptide ratio)',
    'DESIGN:MSstats'             => 'DIA-based Protein ratio (MSstats)',
    'DESIGN:DIA'                 => 'DIA-based Protein ratio (myProMS)',
    'DESIGN:DIA:Abund'           => 'DIA-based Protein Abundance',
    'DESIGN:TDA'                 => 'TDA (PRM/SRM/MRM)',
    'DESIGN:SSPA'                => 'SSP Analysis',
    'DESIGN:PROT_RULER'          => 'Proteomic Ruler',
    'HYDRO'                      => 'Peptides Hydrophobicity'
); #'SILAC'=>'SILAC','ITRAQ'=>'iTRAQ', 'DESIGN:SWATH/MSstats'=>'SWATH-based Protein ratio (MSstats)',

my %filters = (
    "STATUS" => (param('filterStatus')) ? [param('filterStatus')] : [],
    "TYPE" => (param("filterType")) ? [param("filterType")] : [],
    "USER" => (param("filterUser")) ? [param("filterUser")] : [],
    "PROJECT" => (param("filterProject")) ? [param("filterProject")] : [],
    "DATE" => (param("filterDateNumber") && param("filterDateNumber") > 0 && param("filterDateType")) ? param("filterDateNumber")." ".param("filterDateType") : '',
    "ID" => (param("filterID")) ? [split(',', promsMod::cleanParameters(param("filterID")))] : [],
);

# If grouped types exists, split them (i.e. Protein Ratio = DESIGN:myProMS;DESIGN;DESIGN:SimpleRatio;DESIGN:SuperRatio)
if($filters{"TYPE"}) {
    my @splittedTypes;
    foreach my $filterType (@{$filters{"TYPE"}}) {
        push(@splittedTypes, split(';', $filterType));
    }
    $filters{"TYPE"} = \@splittedTypes;
}
    
my $realTime = (param('realTime')) ? param('realTime') : 0;
my $dbh = &promsConfig::dbConnect;

# Get user status
my $userStatus = 'bio';
my @userInfo = &promsMod::getUserInfo($dbh, $userID);
$userStatus = $userInfo[1];

####################
####>Parameters<####
####################
if (param('AJAX')) {
    &ajaxUpdateJobsStatus if(param('AJAX') eq 'update');
    &ajaxDeleteJob        if(param('AJAX') eq 'delete');
    &ajaxStopJob          if(param('AJAX') eq 'stop');
    &ajaxInfosJob         if(param('AJAX') eq 'infoCluster');
    exit;
}

#######################
####>Starting HTML<####
#######################
my ($lightColor, $darkColor) = &promsConfig::getRowColors;
my $rowColor = $lightColor;

print header(-'content-encoding'=>'no',-charset=>'UTF-8'); #(-cache_control=>"no-cache, no-store, must-revalidate"); for AJAX in IE but works better in script called by AJAX function
warningsToBrowser(1);

print qq |
<html>
    <head>
        <title>Jobs Monitoring</title>
        <link rel="stylesheet" href="$promsPath{html}/promsStyle.css" type="text/css">
        <link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.8.2/css/all.css" integrity="sha384-oS3vJWv+0UjzBfQzYUhtDYW+Pj2yciDJxpsK1OYPAYjqT085Qq/1cq5FLXAZQ7Ay" crossorigin="anonymous">
        <style>
            div.job {
                margin: 12px 0 0 12px;
            }
            
            table.jobMeta {
                table-layout: fixed;
                overflow: scroll;
            }
            
            table.jobMeta tr td {
                padding: 2 15px 2 15px;
            }
            
            div.project {
                margin: 23px 0 32px 0;
            }
            
            .job .title3 {
                margin-bottom: 4px;
            }
            
            th.manageCol {
                width: 150px;
            }
            
            th.smallDescCol {
                width: 175px;
                text-align: center;
            }
            
            th.largeDescCol {
                width: 150px;
                text-align: center;
            }
            
            th.startedCol {
                width: 165px;
                text-align: center;
            }
            
            th.statusCol {
                max-width: 400px;
                min-width: 100px;
            }
            
            table#filtersTable {
                margin: auto;
            }
            
            table#filtersTable td {
                text-align: center;
                vertical-align: top;
            }
            
            table#filtersTable tr > th {
                padding-bottom: 7px;
            }
            
            select {
                resize: vertical;
                max-height: 160px;
            }
        </style>
        <script type="text/javascript">
            var noJobErrorUpdate={}; // list of jobs not to be updated because user is reading error message
            var timeOut; // Used to get real time job status
            var nbRefresh = 0;
            var currentQueries = 0; // Used to check if it is waiting for an ajax result
            
            function viewError(errorDivID) {
                var errorDiv = document.getElementById(errorDivID);
                if (errorDiv.style.display == 'none') {
                    errorDiv.style.display = '';
                    noJobErrorUpdate[errorDivID] = 1;
                    
                    // Stop job status monitoring
                    if(timeOut !== undefined) {
                        clearTimeout(timeOut);
                    }
                } else {
                    errorDiv.style.display = 'none';
                    delete noJobErrorUpdate[errorDivID];
                    
                    if(Object.keys(noJobErrorUpdate).length == 0) {
                        timeOut = setTimeout('ajaxMonitorJobs()', $REFRESH_TIME_INTERVAL*1000);
                    }
                }
            }
            
            function toggleProject(projectID) {
                var projectEl = document.getElementById('project:' + projectID);
                if(projectEl !== undefined) {
                    var jobsEl = projectEl.getElementsByClassName('job');
                    
                    //This will hide all matching elements except the first one
                    for(var i = 0; i < jobsEl.length; i++){
                       var jobEl = jobsEl[i];
                       jobEl.style.display = (jobEl.style.display === 'none') ? '' : 'none';
                    }
                    
                    var toggleIcon = document.getElementById('toggleIcon:' + projectID);
                    if(toggleIcon) {
                        toggleIcon.className = (toggleIcon.className == "far fa-plus-square") ? "far fa-minus-square" : "far fa-plus-square";
                    }
                }
            }
        
            // AJAX --->
            function viewJobRessources(jobID, realTime='') {
                var resourcesDiv = document.getElementById('resources:' + jobID);
                if (resourcesDiv.style.display == 'none' \|\| realTime) {
                    // Stop job status monitoring
                    if(timeOut !== undefined) {
                        clearTimeout(timeOut);
                    }
                
                    // Disabled job get ressources button
                    var getRessourcesButton = document.getElementById('realTimeResources:' + jobID);
                    if(getRessourcesButton) {
                        getRessourcesButton.disabled = true
                    }
                
                    //Creation of the XMLHTTPRequest object
                    var resourcesDivInner = document.getElementById('resourcesInner:' + jobID);
                    resourcesDivInner.innerHTML='<img src="$promsPath{images}/scrollbarGreen.gif" />';
                    resourcesDiv.style.display = '';
                    currentQueries++;
                    var XHR = getXMLHTTP();
                    XHR.open("GET", "./monitorJobs.cgi?AJAX=infoCluster&jobID=" + jobID + "&realTime=" + realTime, true);
                    XHR.onreadystatechange=function() {
                        if(XHR.readyState == 4) {
                            if(XHR.response) {
                                var newInnerHTML = XHR.response;
                                
                                if(realTime) { // Remove first line as it corresponds to the cluster job ID generated to retrieve real time information
                                    var lines = newInnerHTML.split('torque6.curie.fr');
                                    if(lines.length > 1) {
                                        newInnerHTML = lines.slice(1)[0];
                                    }
                                }
                                
                                if (resourcesDiv !== "undefined") {
                                    resourcesDivInner.innerHTML = newInnerHTML;
                                } else {
                                    resourcesDivInner.innerHTML = "Data not found";
                                }
                            }
                            
                            getRessourcesButton.disabled = false;
                            
                            currentQueries--;
                            if(currentQueries == 0) {
                                timeOut = setTimeout('ajaxMonitorJobs()', $REFRESH_TIME_INTERVAL*1000);
                            }
                        }
                    }
                    XHR.send(null);
                } else {
                    resourcesDiv.style.display = 'none';
                }
            }
            
            function deleteJob(jobID, act='delete', multi='') {
                if(timeOut !== undefined) {
                    clearTimeout(timeOut);
                }
                
                //Creation of the XMLHTTPRequest object
                var deleteButton = document.getElementById('delete:' + jobID);
                if(deleteButton) {
                    deleteButton.disabled = true;
                }
                
                //Creation of the XMLHTTPRequest object
                var XHR = getXMLHTTP();
                XHR.open("GET", "./monitorJobs.cgi?AJAX=" + act + "&jobID=" + jobID + "&multi=" + multi, true);
                XHR.onreadystatechange=function() {
                    if(XHR.readyState == 4) {
                        if(/RELOAD/.test(XHR.response)) {
                            window.location.reload();
                        } else {
                            alert(XHR.response);
                        }
                    }
                }
                XHR.send(null);
            }
            
            function getXMLHTTP() {
                var xhr = null;
                if(window.XMLHttpRequest) { // Firefox & others
                    xhr = new XMLHttpRequest();
                } else if(window.ActiveXObject) { // Internet Explorer
                    try {
                      xhr = new ActiveXObject("Msxml2.XMLHTTP");
                    } catch (e) {
                        try {
                            xhr = new ActiveXObject("Microsoft.XMLHTTP");
                        } catch (e1) {
                            xhr = null;
                        }
                    }
                } else { // XMLHttpRequest not supported by browser
                    alert("Your browser does not support XMLHTTPRequest objects...");
                }
                return xhr;
            }
    |;
    
    ###################################################
    #### Handles real time automatic status update ####
    ###################################################
    my $statusStr = "filterStatus=".join('&filterStatus=', @{$filters{STATUS}});
    my $typeStr = "filterType=".join('&filterType=', @{$filters{TYPE}});
    my $userStr = "filterUser=".join('&filterUser=', @{$filters{USER}});
    my $projectStr = "filterProject=".join('&filterProject=', @{$filters{PROJECT}});
    my $dateStr = $filters{DATE};
    my $idStr = join(',', @{$filters{ID}});
    
    print qq |
        var nTotJobs = 0;
        function updateJobsStatus(jobDataTxt) {
            if(timeOut !== undefined) {
                clearTimeout(timeOut);
            }
            
            if(jobDataTxt !== undefined) {
                var jobDataLine = jobDataTxt.split('\\n');
                let nJobs = 0;
                
                for (var i=0; i < jobDataLine.length; i++) { // last line is empty                
                    if (!jobDataLine[i].match('#'))
                        continue;
                    
                    var jobData = jobDataLine[i].split('##');
                    
                    // End of a job : Status changed to Error/Done/Stopped 
                    if(jobData[0] == 'RELOAD') {
                        timeOut = setTimeout('window.location.reload()', $REFRESH_TIME_INTERVAL*1000);
                        return;
                    }
                    
                    var realTimeUpdateButton = document.getElementById('realTimeResources:' + jobData[0]);
                    if(realTimeUpdateButton && jobData[1]) {
                        var regex = /Done\|Error\|Stopped/;
                        realTimeUpdateButton.style.display = (regex.test(jobData[1])) ? 'none' : 'inline-block';
                    }
                    
                    // Update job status    
                    var statusSpan = document.getElementById('status:' + jobData[0]);
                    if(statusSpan) {
                        statusSpan.innerHTML = jobData[1];
                        
                    } else { // New job has been detected
                        timeOut = setTimeout('window.location.reload()', $REFRESH_TIME_INTERVAL*1000);
                        return;
                    }
                    
                    // Update job elapsed time
                    var updatedSpan = document.getElementById('updated:' + jobData[0]);
                    if(updatedSpan) {
                        updatedSpan.innerHTML = jobData[2];
                    }
                    
                    // Print error status if any
                    var errorDivID = 'error:' + jobData[0];
                    if (!noJobErrorUpdate[errorDivID] && jobData[3]) {
                        document.getElementById(errorDivID).innerHTML = '<FIELDSET><LEGEND><B>Error message:</B></LEGEND>' + jobData[3] + '</FIELDSET>';
                    }
                    
                    nJobs++;
                }
                
                if(nJobs < nTotJobs) { // A job is over or has been deleted
                    timeOut = setTimeout('window.location.reload()', $REFRESH_TIME_INTERVAL*1000);
                    return;
                } else {
                    nTotJobs =  nJobs;
                }
            }
            
            timeOut = setTimeout('ajaxMonitorJobs()', $REFRESH_TIME_INTERVAL*1000);
        }
        
        function ajaxMonitorJobs() {
            let realTime = (nbRefresh != 0 && nbRefresh % 36 == 0) ? 'true' : ''; // Refresh with memory and time consumption

            //Creation of the XMLHTTPRequest object
            var XHR = getXMLHTTP();
            if (!XHR) { return false; }
            XHR.open("GET","./monitorJobs.cgi?AJAX=update&$statusStr&$typeStr&filterID=$idStr&filterDate=$dateStr&$userStr&$projectStr&realTime=" + realTime, true);
            XHR.onreadystatechange=function() {
                if (XHR.readyState==4) {
                    updateJobsStatus(XHR.responseText);
                    nbRefresh++;
                }
            }
            XHR.send(null);
        }
        // <--- AJAX
        </script>
    </head>
    <body background="$promsPath{images}/bgProMS.gif" >
        <span class="title">Jobs Monitoring</span>&nbsp;&nbsp;&nbsp;
        <input type="button" value="Close window" onclick="window.close()" />
        <br/><br/>
  |;
  
    # Filter by Date
    my $filterDateNumberSelected = (param("filterDateNumber")) ? param("filterDateNumber") : '';
    my $filterDateTypeSelected = (param("filterDateType")) ? param("filterDateType") : '';
    my @filterDateTypes = ('HOUR', 'DAY', 'WEEK', 'MONTH');
  
    print qq |
    <fieldset id="filtersFieldSet" style="margin: 20px auto 45px auto; width:70%;">
        <legend><b>Filter jobs:</b></legend>
        
        <form action='./monitorJobs.cgi' method='GET' style="margin: 7px 0 7px 0;">
            <table id='filtersTable'>
                <tr>
                    <th style='font-size:14px; min-width: 140px;'><label for='filterDateNumber'>From last</label></th>
                    <th style='font-size:14px; min-width: 110px;'><label for='filterStatus'>Status</label></th>
                    <th style='font-size:14px; min-width: 295px;'><label for='filterType'>Type</label></th>
    |; 
    
    if($userStatus =~ /bioinfo|mass|manag/) {
        print "<th style='font-size:14px; min-width: 180px;'><label for='filterUser'>User</label></th>";
        print "<th style='font-size:14px; min-width: 180px;'><label for='filterProject'>Project</label></th>";
    }
    
    print qq |
                    <th></th>
                </tr>
                <tr>
                    <td>
                        <input type="number" id='filterDateNumber' name='filterDateNumber' min="0" max="99" size="2" style="width: 45px;" onKeyUp="if(parseInt(this.value)>parseInt(this.max)){this.value=this.max;}else if(parseInt(this.value)<parseInt(this.min)){this.value=this.min;}" value='$filterDateNumberSelected'>
                        <select name='filterDateType'>
    |;
    
    
    foreach my $filterDateType (@filterDateTypes) {
        print("<option value='$filterDateType'");
        print(" selected") if($filterDateType eq $filterDateTypeSelected);
        print(">".ucfirst(lc($filterDateType))."(s)</option>\n");
    }
    
    print qq |
                        </select>
                    </td>
    |;

    my (@projectsID, @projectsName, @userIDs, @userNames);
    my (@status, @types, @projects, @users);
    my ($status, $types, $projects, $users) = $dbh->selectrow_array("SELECT GROUP_CONCAT(DISTINCT J.JOB_STATUS), GROUP_CONCAT(DISTINCT TYPE), GROUP_CONCAT(DISTINCT J.ID_PROJECT, ';', NAME), GROUP_CONCAT(DISTINCT J.ID_USER, ';', UL.USER_NAME) FROM JOB_HISTORY J LEFT JOIN PROJECT P ON P.ID_PROJECT=J.ID_PROJECT INNER JOIN USER_LIST UL ON UL.ID_USER=J.ID_USER") or die "Couldn't prepare statement: " . $dbh->errstr;
    @status = ('Queued', 'Running', 'Done', 'Stopped', 'Error');
    @types = split(',', $types) if($types);
    if ($projects) {
        @projects = split(',', $projects);
        @projectsID = map((split(';', $_))[0], @projects);
        @projectsName = map((split(';', $_))[1], @projects);
    }
    if ($users) {
        @users = split(',', $users) ;
        @userIDs = map((split(';', $_))[0], @users);
        @userNames = map((split(';', $_))[1], @users);
    }

    # Filter by Status
    print qq |
                    <td>
                        <select name='filterStatus' id='filterStatus' multiple>
    |;
    
    my $allSelected = (scalar @{$filters{STATUS}} == 0) ? 'selected' : '';
    print("<option value='' $allSelected>All</option>");
    foreach my $status (@status) {
        print("<option value='$status'");
        print(" selected") if(grep( /^$status$/, @{$filters{STATUS}}));
        print(">".$status."</option>\n");
    }
    
    print qq |
                        </select>
                    </td>
    |;
    
    # Filter by Type
    print qq |
                    <td>
                        <select name='filterType' id='filterType' multiple>
    |;

    push(@types, ('Import', 'Quantification'));
    @types = sort(@types);
    my %jobTypes;
    my $category = "";
    $allSelected = (scalar @{$filters{TYPE}} == 0) ? 'selected' : '';
    foreach my $type (sort @types) {
        if($type =~ /([^[]+)\[(.+?)\]/) {
            $category = $1;
            my $typeDesc = $2;
            $category =~ s/^\s+|\s+$//;
            $typeDesc = $quantifProcesses{$typeDesc} if($category eq "Quantification" && $quantifProcesses{$typeDesc});
            push(@{$jobTypes{($category) ? $category : "Others"}{$typeDesc}}, $type);
        } else {
            %{$jobTypes{$type}} = ();
        }
    }
    
    print("<option value='' $allSelected>All</option>");
    foreach my $category (sort keys %jobTypes) {
        print("<option value='$category'");
        print(" selected") if($category && grep( /^$category$/, @{$filters{TYPE}}));
        print(">".$category."</option>\n");
        foreach my $typeDesc (sort {lc $a cmp lc $b} keys %{$jobTypes{$category}}) {
            my $types = join(';', @{$jobTypes{$category}{$typeDesc}});
            print("<option value='$types'");
            my $regex = quotemeta($types);
            $regex =~ s/;/\$|^/;
            print(" selected") if($regex && grep( /^$regex$/, @{$filters{TYPE}}));
            print(">&nbsp;&nbsp;&nbsp;&nbsp;".$typeDesc."</option>\n");
        }
    }
        
    print qq |
                        </select>
                    </td>
    |;

    if($userStatus =~ /bioinfo|mass|manag/) {
        print qq |
                    <td>
                        <select id='filterUser' name='filterUser' multiple>
        |;
        
        my @indices = sort { $userNames[$a] cmp $userNames[$b] }  0 .. $#userNames;
        @userIDs = @userIDs[@indices];
        @userNames = @userNames[@indices];
        $allSelected = (scalar @{$filters{USER}} == 0) ? 'selected' : '';
        print("<option value='' $allSelected>All</option>");
        for(my $i=0; $i < scalar @userIDs; $i++) {
            my ($userID, $userName) = ($userIDs[$i], $userNames[$i]);
            print("<option value='$userID'");
            print(" selected") if(grep( /^$userID/, @{$filters{USER}}));
            print(">".$userName."</option>\n");
        }
                    
        print qq |
                        </select>
                    </td>
                    <td>
                        <select id='filterProject' name='filterProject' multiple>
        |;
    
        $allSelected = (scalar @{$filters{PROJECT}} == 0) ? 'selected' : '';
        print("<option value='' $allSelected>All</option>");
        @indices = sort { $projectsName[$a] cmp $projectsName[$b] }  0 .. $#projectsName;
        @projectsName = @projectsName[@indices];
        @projectsID = @projectsID[@indices];
        for(my $i=0; $i < scalar @projectsID; $i++) {
            my ($projectID, $projectName) = ($projectsID[$i], $projectsName[$i]);
            print("<option value='$projectID'");
            print(" selected") if(grep( /^$projectID/, @{$filters{PROJECT}}));
            print(">".$projectName."</option>\n");
        }
                    
        print qq |
                        </select>
                    </td>
        |;    
    }
    
    print qq |
                    <td style='vertical-align: middle;'><input type='submit' value='Apply filters'/></td>
                </tr>
            </table>
        </form>
    </fieldset>
    |;


#########################################
#### Check if some jobs are too old  ####
#########################################
checkForOldJobs();


##############################
#### Get jobs to display  ####
##############################
my @jobs = getJobs(\%filters);
if(!@jobs) { #### Check if there are some jobs to display ####
	print qq |
        <span class=\"title2\">No job found.</span>
    |;
}
else { ########### DISPLAY JOBS ############################
    
    my $currentProject = my $currentJob = my $projectDisplayName = '';
    my $currentlyInJob = 0;
    my $nbSubJobs = 0;
    my $hasRunningJobs = 0;
    
    foreach my $jobRef (@jobs) {
        my %job = %{$jobRef};
        $hasRunningJobs = 1 if($job{"status"} !~ /Done|Error|Stopped/ and !$hasRunningJobs);
        
        ## Check access rights
        if($job{projectID} > -1) {
            my @userInfo = &promsMod::getUserInfo($dbh, $userID, $job{projectID});
            my $projectAccess = ${$userInfo[2]}{$job{projectID}};
            $userStatus = $userInfo[1];
            $projectDisplayName = "Project '$job{projectName}' $projectAccess";
            next if ($projectAccess !~ /bioinfo|mass|manag/ && $userID ne $job{userID});
        } else {
            next if ($userStatus !~ /bioinfo|mass|manag/);
            $projectDisplayName = "$job{projectName}";
        }
        
        my ($masterJobCode, $jobRank) = split('\.', $job{dirName});
        if($currentlyInJob && $masterJobCode ne $currentJob) {
            print qq |
                        </tbody>
                    </table>
                </div>
            |;
            $currentlyInJob = 0;
            my $deleteDisplayStr = ($hasRunningJobs) ? "Stop and remove all" : "Remove all";
            $rowColor = $lightColor;
            print qq |
                <script>
                    var deleteMasterButton = document.getElementById('deleteMaster:$currentJob');
                    if(deleteMasterButton) {
                        deleteMasterButton.innerHTML = '<button id="delete:$currentJob" onclick="deleteJob(' + "'$currentJob', 'delete', 'multi')" + '"><!--<i class="fa fa-trash" aria-hidden="true"></i>-->$deleteDisplayStr</button>';
                    }
                </script>
             
            | if($nbSubJobs > 1);
            $nbSubJobs = 0;
            $hasRunningJobs = 0;
        }
        $currentJob = $masterJobCode;
        $nbSubJobs++;
        
        if($currentProject ne $job{projectID}) {
            print("</div>") if($currentProject);
            print qq |
                <div class='project' id='project:$job{projectID}'>
                    <div class="title2" style='cursor:pointer' onclick="toggleProject('$job{projectID}');"><i id='toggleIcon:$job{projectID}' class="far fa-minus-square"></i> $projectDisplayName :</div>
            |;
            $currentProject = $job{projectID};
        }
        
        printJob(\%job, $currentlyInJob, $userStatus);
        
        $rowColor = ($rowColor eq $lightColor) ? $darkColor : $lightColor;
        $currentlyInJob = 1;
    }
    
    my $deleteDisplayStr = ($hasRunningJobs) ? "Stop and remove all" : "Remove all";
    print("</table></div>") if($currentlyInJob); # End of job
    print("</div>"); # End of project
    print qq |
        <script>
            var deleteMasterButton = document.getElementById('deleteMaster:$currentJob');
            if(deleteMasterButton) {
                deleteMasterButton.innerHTML = '<button id="delete:$currentJob" onclick="deleteJob(' + "'$currentJob', 'delete', 'multi')" + '"><!--<i class="fa fa-trash" aria-hidden="true"></i>-->$deleteDisplayStr</button>';
            }
        </script>
     
    | if($nbSubJobs > 1);
}

print qq |
            <script type='text/javascript'>
                 timeOut = setTimeout('ajaxMonitorJobs()', 2000);
            </script>
        </body>
    </html>
|;

$dbh->disconnect;

##############################<<< SUBROUTINE >>>################################
sub ajaxUpdateJobsStatus {
	print header(-charset=>'UTF-8',-type=>'text/plain',-cache_control=>"no-cache, no-store, must-revalidate"); # for AJAX in IE
	warningsToBrowser(1);

	###<Fetching job list
	my @jobs = getJobs(\%filters);
	my $responseStrg = '';
    
    foreach my $jobRef (@jobs) {
        my %job = %{$jobRef};
        next if($job{status} !~ /Queued|Running/); # Do not check job status if it is : Stopped / Error / Done
        
        # Get status + error messages of current job
        my ($error, $status) = ('', '');
        my $startedTime = Time::Piece->strptime("$job{startedDate}", "%Y-%m-%d %H:%M:%S");
        my $now = Time::Piece->new()+3600; # Shift of +2 hours to correspond to GMT+1
        if($now > 1585447200 && $now < 1603584000) {
            $now+=3600;
        }
        my $timeDiff = ($now - $startedTime)->pretty;
        my $clusterError = ''; #(-s "$job{srcPath}/PBSerror.txt") ? $clusterInfo{'checkError'}->("$job{srcPath}/PBSerror.txt") : '';
        
        # Check cluster jobs for completion/errors
        if($job{"processID"}) {
            addJobInfos(\%job, $realTime);
        }
        
        ($status, $error) = getJobStatus(\%job);
        
        $error =~ s/\n/<br\/>/g if($error);
        $status =~ s/\n/<br\/>/g if($status);
        
        # Update job status if needed
        my $newStatus = '';
        
        if($job{"status"} ne 'Error' && (($error && $error !~ /Warning/) || $job{'infos'}{'clusterError'} || $clusterError || $job{'infos'}{'exceededMem'} || $job{'infos'}{'exceededTime'})) {
            $newStatus = 'Error';
        } elsif((($job{"processID"} && $job{'infos'}{'allExited'}) || (-s $job{logPath} && `tail -1 $job{logPath}` =~ /Ended/)) && $job{"status"} !~ /Error|Stopped|Done/) {
            $newStatus = 'Done';
        } elsif($job{'status'} eq 'Queued' && (-s $job{'logPath'} || $job{'infos'}{'minStart'})) {
            $newStatus = 'Running';
        }
        
        if($newStatus) {
            my $newInfos = "JOB_STATUS='$newStatus'";
            if($newStatus ne 'Running') {
                my $endedStr = 'NOW()'; # Put local current time as default ending of other options fail
                if($job{'infos'}{'maxEnd'} && $job{'infos'}{'minStart'}) { # Get cluster job end time if it was run on server (more precise)
                    my $timeDiffSec = $job{'infos'}{'maxEnd'}-$job{'infos'}{'minStart'};
                    $endedStr = "DATE_ADD(STARTED, INTERVAL $timeDiffSec SECOND)";
                } elsif(-s $job{logPath}) { # Get job end time by looking at log file last modification
                    my $endJob = stat($job{logPath})->[9]; # in sec
                    my $startJob = Time::Piece->strptime($jobRef->{startedDate}, '%Y-%m-%d %H:%M:%S')->strftime("%s");;
                    my $timeDiffSec = $endJob-$startJob;
                    $endedStr = "DATE_ADD(STARTED, INTERVAL $timeDiffSec SECOND)";
                }
                
                $newInfos .= ", ENDED=$endedStr";
            }
            
            $dbh->do("UPDATE JOB_HISTORY SET $newInfos WHERE ID_JOB='$job{ID}'");
            $dbh->commit;
            
            # Refresh if job end with either Done or Error status
            if($newStatus =~ /Done|Error/) {
                $responseStrg = "RELOAD".$responseStrg;
                
                if($newStatus eq 'Error') {
                    stopJobProcess($jobRef);
                }
            }
        }
        
        $responseStrg .= "$job{ID}##$status##($timeDiff)##$error\n";        
    }
    
	print $responseStrg;
}

sub getJobStatus {
    my ($jobRef) = @_;
    
    return ('', '') if(!$jobRef);
    
    my ($status, $error, $warning) = ('', '', '');
    my %job = %{$jobRef};
    my %clusterInfo;
    my $clusterError='';
    if ($job{'processID'} =~ /^C/) { # job is on cluster
        %clusterInfo = &promsConfig::getClusterInfo;
        $clusterError = (-s "$job{srcPath}/PBSerror.txt") ? $clusterInfo{'checkError'}->("$job{srcPath}/PBSerror.txt") : '';
    }
    
    if($job{status} eq 'Stopped') {
        return ('[ <b>Stopped</b> ]', '');
    } elsif($job{status} eq 'Done') {
        if(-s $job{logPath} && `tail -1 $job{logPath}` !~ /Ended/) {
            $warning = 'Job ended on cluster prematurely. Check server logs for potential errors.';
        }
    } elsif($job{status} ne 'Done' && -s $job{logPath} && `tail -2 $job{logPath}` =~ /Ended/) {
        return ('[ <b>Done</b> ]', '');
    } elsif($job{status} ne 'Error' && ($job{infos}{clusterError} || $clusterError || $job{infos}{"exceededMem"} || $job{infos}{"exceededTime"})) {
        open(ERROR, '>', $job{"errorPath"}) or die "Could not open file '$job{errorPath}' $!";
        
        if($job{'infos'}{'exceededTime'}) {
            print ERROR ("Cluster job \"$job{infos}{exceededTime}\" exceeded allowed run time !\n");
        } elsif($job{'infos'}{'exceededMem'}) {
            print ERROR ("Cluster job \"$job{infos}{exceededMem}\" exceeded allowed memory !\n");
        } elsif($job{'infos'}{'clusterError'} || $clusterError) {
            if($clusterError) {
                $status .= " Cluster returned error :";
                print ERROR "|$clusterError|";
            } else {
                $status .= "Job ended on cluster with error status";
                print ERROR " ($job{infos}{clusterError}), check server logs.";
            }
        }
        
        close ERROR;
    } elsif($job{status} eq 'Queued') {
        $status = "[ <b>Waiting for job to start...</b> ]";
    } elsif(-s $job{logPath}) {
        my $lastLine = `tail -1 $job{logPath}`;
        chomp $lastLine;
        
        # Special task for MassCroq XIC alignment (TODO : Put it in runXICQuantification)
        if($job{type} =~ /Quantification/ && $job{features}{TYPE} eq 'XICMCQ' && $lastLine =~ /^2\/3 / && -e "$job{srcPath}/masschroq_status2.txt") { # very long run not yet ended
            my $sampNum = 0;
            open(MCQ, "$job{srcPath}/masschroq_status2.txt");
            while (<MCQ>) {
                if (/^Parsing XML/) {
                    $status = 'Parsing : 2.1/3 Parsing XML input file...';
                    $sampNum=0;
                } elsif (/^Alignment method/) {
                    $status = 'Aligning : 2.2/3 Aligning MS run';
                    $sampNum=0;
                } elsif (/^Quantification method/) {
                    $status = 'Quantification : 2.3/3 Quantifying in MS run'; 
                    $sampNum=0;
                } elsif (/ quantification finished/) {
                    $status = 'Post matching : 2.4/3 Post matching in MS run';
                    $sampNum=0;
                } elsif (($status eq 'parsing' && /MS run 'samp\d+'/) || / MS run 'samp\d+'/) {
                    $sampNum++;
                }
            }
            close MCQ;
            $status .=" #$sampNum" if($sampNum);
        } else { # Get current status
            $status .= $lastLine;
        }
    }
    
    # Check for warning
    my $allLogFile = (-s $job{logPath}) ? `cat $job{logPath}` : '';
    
    if($allLogFile =~ /You had (\d+) "Peptide has unknown modification" errors/) {
            $warning = "You had $1 \"Peptide has unknown modification\" errors during SpectraST process.";
    } elsif($allLogFile =~ /warning/i) {
        $warning = `tail -1 $job{logPath}`;
        chomp $warning;
        my $stats = `tail -2 $job{logPath} | head -1`;
        chomp $stats;
        if($stats !~ /Ended/) {
            $status = $stats;
        }
    }
    
    # Check for errors
    if ($job{'processID'} =~ /^C/) { # job is on cluster
       $error = $clusterInfo{'checkError'}->($job{'errorPath'});
    }
    else {
        $error = `cat $job{errorPath}` if(-s $job{'errorPath'}); # Error file : $prsHomeDir/current/$anaID\_$job{dirName}\_error.txt
    }
    if($job{status} ne 'Queued') {
        $status = " -> ".$status if($job{status} =~ /Running|Error/ && $status);
        $status = ($error) ? "<b>Error</b>".$status : "<b>$job{status}</b>".$status;
        $status = "[ $status ]";
    }
    
    if($error) {
        chomp($error);
        $error =~ s/\n/<br\/>/g;
        $status .= " <INPUT type=\"button\" class='errorButton' value=\"Show/Hide message\" onclick=\"viewError('error:$job{ID}')\">";
    } elsif($warning) {
        $status.= " <FONT color=\"red\"><B>**WARNING**</B></FONT> " if($warning);
        $status .= " <INPUT type=\"button\" class='warningButton' value=\"Show/Hide message\" onclick=\"viewError('error:$job{ID}')\">";
        $error = "Warning : $warning" if($warning);
    }
    
    return ($status, $error);
}

sub getJobs {
    my @jobs;
    my ($filtersRef) = @_;

    my $getJobsQuery = "SELECT J.ID_JOB, J.ID_USER, J.ID_PROJECT, J.ID_JOB_CLUSTER, J.TYPE, J.JOB_STATUS, J.FEATURES, J.SRC_PATH, J.LOG_PATH, J.ERROR_PATH, J.STARTED, J.ENDED, IFNULL(UL.USER_NAME, 'Unknown'), IFNULL(P.NAME, 'Global') FROM JOB_HISTORY J LEFT JOIN USER_LIST UL ON UL.ID_USER=J.ID_USER LEFT JOIN PROJECT P ON P.ID_PROJECT=J.ID_PROJECT WHERE 1";
    
    if($filtersRef->{ID} && @{$filtersRef->{ID}}) {
        if($filtersRef->{MULTI}) {
            $getJobsQuery .= " AND J.ID_JOB REGEXP '".join("|", @{$filtersRef->{ID}})."'";
        } else {
            $getJobsQuery .= " AND J.ID_JOB IN('".join("','", @{$filtersRef->{ID}})."')";
        }
    } else {
        $getJobsQuery .= " AND J.JOB_STATUS IN ('".join("','", @{$filtersRef->{STATUS}})."') " if($filtersRef->{STATUS} && @{$filtersRef->{STATUS}});
        $getJobsQuery .= " AND J.ID_USER IN ('".join("','", @{$filtersRef->{USER}})."') " if($filtersRef->{USER} && @{$filtersRef->{USER}});
        $getJobsQuery .= " AND J.ID_PROJECT IN ('".join("','", @{$filtersRef->{PROJECT}})."') " if($filtersRef->{PROJECT} && @{$filtersRef->{PROJECT}});
        
        my $filtersTypeStr = (defined $filtersRef->{TYPE}) ? join("|", @{$filtersRef->{TYPE}}) : '';
        $filtersTypeStr =~ s/\[/\\\\[/g;
        $filtersTypeStr =~ s/\]/\\\\]/g;
        $getJobsQuery .= " AND TYPE REGEXP '$filtersTypeStr'" if($filtersRef->{TYPE} && @{$filtersRef->{TYPE}});
        
        if($filtersRef->{DATE}) {
            my $beforeAfterSign = '>';
            if(substr($filtersRef->{DATE}, 0, 1) eq '+') {
                $filtersRef->{DATE} = substr($filtersRef->{DATE}, 1);
                $beforeAfterSign = '<';
            }
            $getJobsQuery .= " AND STARTED $beforeAfterSign NOW() - INTERVAL ".$filtersRef->{DATE};
        }
    }
    $getJobsQuery .= " ORDER BY J.ID_PROJECT, Year(STARTED) DESC, Month(STARTED) DESC, Day(STARTED) DESC, Hour(STARTED) DESC, Minute(STARTED) DESC, J.TYPE, ID_JOB";
    my $sthJobs = $dbh->prepare($getJobsQuery);
    $sthJobs->execute();
    while (my ($ID, $userID, $projectID, $processID, $type, $status, $features, $srcPath, $logPath, $errorPath, $startedDate, $endDate, $userName, $projectName) = $sthJobs->fetchrow_array) {
        my %features = ();
        if($features) {
            foreach my $feature (split(';', $features)) {
                my ($name, $value) = split('=', $feature);
                $features{$name} = $value;
            }
        }
        
        my %job = (
            ID             => $ID,
            status         => $status,
            processID      => $processID,
            projectID      => ($projectID) ? $projectID : -1,
            projectName    => $projectName,
            type           => $type,
            features       => \%features,
            srcPath        => $srcPath,
            dirName        => basename($srcPath),
            parentPath     => dirname($srcPath),
            logPath        => $logPath,
            errorPath      => $errorPath,
            startedDate    => $startedDate,
            endDate        => $endDate,
            userID         => $userID,
            userName       => $userName,
        );
        push(@jobs, \%job);
    }
    $sthJobs->finish;
    
    return @jobs;
}

sub ajaxInfosJob {
    print header(-charset=>'UTF-8',-type=>'text/plain',-cache_control=>"no-cache, no-store, must-revalidate"); # for AJAX in IE
	warningsToBrowser(1);
    exit unless(param('jobID'));

    my @jobs = getJobs({ID => [split(',', param('jobID'))]});
    if(@jobs) {
        my %job = %{$jobs[0]};
        
        if($job{processID}) {
            printJobInfos(\%job, $realTime);
        } else {
            print("No job found on cluster.");
        }
    } else {
        print("Job(s) ".param('jobID')." not found");
    }
}

sub ajaxDeleteJob {
    print header(-charset=>'UTF-8',-type=>'text/plain',-cache_control=>"no-cache, no-store, must-revalidate"); # for AJAX in IE
	warningsToBrowser(1);
    exit unless(param('jobID'));
    my @jobs = getJobs({ID => [split(',', param('jobID'))], MULTI => param('multi')});
    
    if(@jobs) {
        my $allJobsID = '';
        foreach my $jobRef (@jobs) {
            deleteJob($jobRef);
        }
    }
    print("RELOAD");
}

sub ajaxStopJob {
    print header(-charset=>'UTF-8',-type=>'text/plain',-cache_control=>"no-cache, no-store, must-revalidate"); # for AJAX in IE
	warningsToBrowser(1);
    exit unless(param('jobID'));
    my @jobs = getJobs({ID => [split(',', param('jobID'))], MULTI => param('multi')});
    
    if(@jobs) {
        my $allJobsID = '';
        my $dbh=&promsConfig::dbConnect;
        
        $dbh->do("UPDATE JOB_HISTORY SET JOB_STATUS='Stopped' WHERE ID_JOB IN ('".join("', '", split(',', param('jobID')))."')"); # Set job status to stopped
        $dbh->commit;
        
        foreach my $jobRef (@jobs) {
            # Stop associated jobs if any
            stopJobProcess($jobRef);
        }
        
        print("RELOAD");
        $dbh->disconnect;
    } else {
        print("An error has occured, job(s) could not be deleted");
    }
}

# Delete a job + do job type-related actions
sub deleteJob {
    my ($jobRef) = @_;
    return if(!$jobRef);
    my %job = %{$jobRef};
    
    my $dbh=&promsConfig::dbConnect;
    $dbh->do("DELETE FROM JOB_HISTORY WHERE ID_JOB='$job{ID}'"); # Set job to deleted
    $dbh->commit;
    
    # Clean job data
    cleanJobData(\%job);
    
    # Stop associated jobs if required
    stopJobProcess(\%job) if($job{status} !~ /Error|Done|Stopped/);
    
    $dbh->disconnect;
}

sub cleanJobData {
    my ($jobRef) = @_;
    return if(!$jobRef);
    
    my %job = %{$jobRef};
    my $dbh=&promsConfig::dbConnect;
    
    if($job{type} =~ /Quantification/) {
        my $quantiID = $job{features}{"ID_QUANTIFICATION"};
        my $quantiType = $job{features}{"TYPE"};
        if ($quantiID && $quantiType ne 'XICCORR' && $job{status} !~ /Done/) { # quantif is recorded in DB BUT all data still in temp dir
            my $projID = ($job{projectID}) ? $job{projectID} : &promsMod::getProjectID($dbh, $quantiID, 'quantification');
            &promsQuantif::deleteQuantification($dbh, $projID, $quantiID, {KEEP_HISTORY=>1});
            $dbh->commit;
        }
        
        # Delete files
        system("rm -f $promsPath{tmp}/quantification/current/*$job{ID}*");
    } if($job{type} =~ /Import/) {
        if($job{type} =~ /MaxQuant/) {
            my $expID = $job{features}{"ID_EXPERIMENT"};
            if($expID) {
                my ($hasChild) = $dbh->selectrow_array("SELECT 1 FROM SAMPLE WHERE ID_EXPERIMENT=$expID");
                #if ($hasChild) {print "##DELETE_EXP\n";} # TODO Delete experiments corresponding to MaxQuant import
            }
        } elsif($job{type} =~ /DIA|TDA/) {
            if ($job{features}{"ACTION"} eq 'library') {
                $dbh->do("DELETE FROM SWATH_LIB WHERE ID_SWATH_LIB=$job{itemID}");
                $dbh->commit;
            }
        }
    } elsif($job{type} eq 'Phospho') {
        system("rm -f $promsPath{tmp}/phosphoRS/current/*$job{ID}*");
    } elsif($job{type} =~ /Functional Analysis/) {  # Only GSEA for now
        my $gseaID = $job{features}{"ID_GSEA"};
        if ($gseaID && $job{status} !~ /Done/) { # analysis recorded in DB but data still in temp dir
            my $projID = ($job{projectID}) ? $job{projectID} : &promsMod::getProjectID($dbh, $gseaID, 'PATHWAY_ANALYSIS');
            $dbh->do("DELETE FROM PATHWAYANA_QUANTIFICATION WHERE ID_PATHWAY_ANALYSIS=$gseaID") || die $dbh->errstr();
            $dbh->do("DELETE FROM PATHWAY_ANALYSIS WHERE ID_PATHWAY_ANALYSIS=$gseaID") || die $dbh->errstr();
            $dbh->commit;
        }
    }
    
    ## Delete job directory
    if (-e "$job{srcPath}") {
        system("rm -rf $job{srcPath}");
    }
    
    $dbh->disconnect;
}

sub checkForOldJobs {
    my @jobs = getJobs({DATE => "+$LIMIT_STORAGE"});
    
    if(@jobs) {
        my $allJobsID = '';
        foreach my $jobRef (@jobs) {
            deleteJob($jobRef);
        }
    }
}


## Stop Cluster jobs / Processes related to a specific job
sub stopJobProcess {
    my ($jobRef) = @_;
    return if (!$jobRef);
    
    my %job = %{$jobRef};
    if($job{processID}) {
        addJobInfos(\%job) if(!$job{infos});
        my (@clusterJobIDToKill, @processIDToKill);
        
        foreach my $jobID (split(';', $job{processID})) {
            my $jobPlatform = substr($jobID,0,1);
            my $jobID = substr($jobID,1);
            
            if($jobPlatform eq 'C' && !$job{infos}{$jobID}{end}) {
                push(@clusterJobIDToKill, $jobID);
            } elsif($jobPlatform eq 'L' && $job{infos}{$jobID}{status}) {
                push(@processIDToKill, $jobID);
            }
        }
        
        if(@clusterJobIDToKill) {
            my %clusterInfo = &promsConfig::getClusterInfo;
            $clusterInfo{'killJob'}->(\@clusterJobIDToKill);
        }
        
        if(@processIDToKill) {
            print("kill -9 ".join(' ', @processIDToKill)."&\n<br/>");
            system("kill -9 ".join(' ', @processIDToKill)."&");
        }
    }
    return;
}

sub addJobInfos {
    my ($jobRef, $realTime) = @_;
    return if(!$jobRef);

    my %processFlags = (
        'D'=>'Uninterruptible sleep (usually IO)',
        'R'=>'Running or runnable (on run queue)',
        'S'=>'Interruptible sleep (waiting for an event to complete)',
        'T'=>'Stopped, either by a job control signal or because it is being traced',
        'W'=>'paging',
        'X'=>'dead (should never be seen)',
        'Z'=>'Defunct ("zombie") process, terminated but not reaped by its parent',
    );
    
    if($jobRef->{"processID"}) {
        my %clusterInfo = &promsConfig::getClusterInfo;
        my $clusterError = 0;
        my $allExited = -1;
        my ($minStart, $maxEnd);
        my @clusterJobsID;
        
        foreach my $jobID (split(';', $jobRef->{"processID"})) {
            my $jobPlatform = substr($jobID,0,1);
            my $jobID = substr($jobID,1);
            
            if($jobPlatform eq 'C') {
                push(@clusterJobsID, $jobID);
            } else {
                my $localProcessInfo = `ps -p $jobID -o uid,pid,ppid,pri,status,pcpu,c,pmem,rss,sz,comm | tail -1 | sed 's/ \+/ /g' | tr -d '\n'`;
                my @infosHeader = ("uid", "pid", "ppid", "priority", "status", "pcpu", "cpu", "pmem", "mem", "vmem", "cmd");
                my %jobInfos = ();
                my $i = 0;
                
                if($localProcessInfo && $localProcessInfo !~ /\sUID/) {
                    foreach my $info (split(' ', $localProcessInfo)) {
                        $jobInfos{$infosHeader[$i]} = $info if($infosHeader[$i]);
                        $i++;
                    }
                    $jobInfos{"uid"} = $localProcessInfo;

                    $jobRef->{"infos"}{$jobID} = \%jobInfos;
                    $allExited = 0;
                } else {
                    $allExited = 1 if($allExited != 0);
                }
            }
        }
        
        if(@clusterJobsID) {
            my %jobsInfo = $clusterInfo{'jobInfos'}->(\@clusterJobsID, $realTime);
            %{$jobRef->{"infos"}} = ($jobRef->{"infos"}) ? (%{$jobRef->{"infos"}}, %jobsInfo) : %jobsInfo;
            
            if(scalar keys %jobsInfo > 0) {
                foreach my $jobID (keys %jobsInfo) {
                    my %jobInfo = %{$jobsInfo{$jobID}};

                    # Check for memory overuse
                    if($jobRef->{"infos"}{$jobID}{"resources_used.mem"} && $jobRef->{"infos"}{$jobID}{"resources_used.mem"} ne '?' && !$jobRef->{'infos'}{'exceededMem'}) {
                        my $memUse = sprintf("%.2f", (substr($jobRef->{"infos"}{$jobID}{"resources_used.mem"}, 0, -2)));
                        $memUse = ($jobRef->{"infos"}{$jobID}{"resources_used.mem"} =~ /gb/) ? $memUse*1024*1024 : ($jobRef->{"infos"}{$jobID}{"resources_used.mem"} =~ /mb/) ? $memUse*1024 : $memUse;
                        my $memAsked = sprintf("%.2f", (substr($jobRef->{"infos"}{$jobID}{"Resource_List.mem"}, 0, -2)));
                        $memAsked = ($jobRef->{"infos"}{$jobID}{"Resource_List.mem"} =~ /gb/) ? $memAsked*1024*1024 : ($jobRef->{"infos"}{$jobID}{"Resource_List.mem"} =~ /mb/) ? $memAsked*1024 : $memAsked;
                        $jobRef->{'infos'}{'exceededMem'} = $jobRef->{"infos"}{$jobID}{"jobname"} if($memUse+1024 >= $memAsked);
                    }
                                           
                    ## Check for time overspend
                    if($jobRef->{"infos"}{$jobID}{"Walltime.Remaining"} && $jobRef->{"infos"}{$jobID}{"Walltime.Remaining"} <= 0 && !$jobRef->{'infos'}{'exceededTime'}) {
                        $jobRef->{'infos'}{'exceededTime'} = $jobRef->{"infos"}{$jobID}{"jobname"};
                    }
                
                    $clusterError = $jobRef->{'infos'}{$jobID}{'Exit_status'} if(!$clusterError && $jobRef->{'infos'}{$jobID}{'Exit_status'} && int($jobRef->{'infos'}{$jobID}{'Exit_status'}) < 0);
                    $maxEnd = $jobRef->{"infos"}{$jobID}{"end"} if(!$maxEnd || ($jobRef->{"infos"}{$jobID}{"end"} && $jobRef->{"infos"}{$jobID}{"end"} > $maxEnd));
                    $minStart = $jobRef->{'infos'}{$jobID}{"start"} if(!$minStart || ($jobRef->{"infos"}{$jobID}{"start"} && $jobRef->{"infos"}{$jobID}{"start"} < $minStart));
                    
                    if(!defined $jobRef->{"infos"}{$jobID}{"Exit_status"}) {
                        $allExited = 0;
                    } elsif($allExited != 0) {
                        $allExited = 1;
                    }
                }
            } else {
                $allExited = 0;
            }
        }
        
        $jobRef->{'infos'}{'allExited'} = $allExited;
        $jobRef->{'infos'}{'clusterError'} = $clusterError;
        $jobRef->{'infos'}{'maxEnd'} = $maxEnd;
        $jobRef->{'infos'}{'minStart'} = $minStart;
    }
}

sub printJobInfos {
    my ($jobRef, $realTime) = @_;
    
    return if(!$jobRef);

    addJobInfos($jobRef, $realTime);
    my %job = %{$jobRef};
    my $jobResourcesStr = '';
    
    foreach my $jobID (split(';', $job{"processID"})) {
        my $jobPlatform = substr($jobID,0,1);
        $jobID = substr($jobID,1);
        
        $jobResourcesStr .= "<br/>" if($jobResourcesStr);
        if($jobPlatform eq 'C') {
            if(!$job{infos}{$jobID} || !$job{infos}{$jobID}{jobname}) {
                my $norClusterStr = ($realTime) ? ' nor in cluster job summary' : '';
                $jobResourcesStr .= qq |
                    <table>
                        <tr style='font-weight:bold'><td>Cluster job ID</td><td>$jobID</td></tr>
                        <tr><td style='font-weight:bold'>Status</td><td><span style='color:orange;font-weight:bold'>Not available in log files$norClusterStr yet</span><br/><i>(Get real time monitoring for more infos)</i></td></tr>
                    </table>
                |;
                next;
            }
            
            my %jobInfo = %{$job{infos}{$jobID}};
            $jobInfo{"resources_used.mem"} = (!$jobInfo{"resources_used.mem"}) ? '?' :
                                             ($jobInfo{"Resource_List.mem"} =~ /gb/) ?
                                                sprintf("%.2f", (substr($jobInfo{"resources_used.mem"}, 0, -2))/1024.0/1024.0) :
                                                sprintf('%.2f', (substr($jobInfo{"resources_used.mem"}, 0, -2))/1024.0);
            $jobInfo{"vmem"} = (!$jobInfo{"resources_used.vmem"}) ? '?' :
                                             (int(substr($jobInfo{"resources_used.vmem"}, 0, -2)) > 1024*1024) ?
                                                sprintf("%.2f", (substr($jobInfo{"resources_used.vmem"}, 0, -2))/(1024.0*1024.0))."gb" :
                                                sprintf('%.2f', (substr($jobInfo{"resources_used.vmem"}, 0, -2))/1024.0)."mb";
            $jobInfo{"resources_used.walltime"} = '?' if(!$jobInfo{"resources_used.walltime"});
            
            
            $jobInfo{"status"} = "<span style='color:orange'>Unknown</span>";
            if(defined $jobInfo{"Exit_status"}) {
                if($jobInfo{"Exit_status"} == 0) {
                    $jobInfo{"status"} = "<span style='color:green'>Ended (properly)</span>";
                } elsif($jobInfo{"Exit_status"} > 128) {
                    $jobInfo{"status"} = "<span style='color:red'>Killed</span>";
                } else {
                    $jobInfo{"status"} = "<span style='color:red'>Ended (Error : $jobInfo{Exit_status})</span>";
                }
            } elsif($jobInfo{"start"}) {
                $jobInfo{"status"} = "<span style='color:orange'>Running (Since ".strftime("%Y-%m-%d %H:%M:%S", localtime($jobInfo{"start"})).")</span>";
            } elsif($jobInfo{"qtime"}) {
                $jobInfo{"status"} = "<span style='color:grey'>Waiting to start</span>";
            }
            
            $jobResourcesStr .= qq |
                <table>
                    <tr style='font-weight:bold'><td>Job Name</td><td>$jobInfo{jobname}</td></tr>
                    <tr style='font-weight:bold'><td>Cluster job ID</td><td>$jobID</td></tr>
                    <tr style='font-weight:bold'><td>Status</td><td>$jobInfo{status}</td></tr>
                    <tr><td>Memory use</td><td><b>$jobInfo{'resources_used.mem'}/$jobInfo{'Resource_List.mem'}</b></td></tr>
                    <tr><td>VMemory use</td><td>$jobInfo{'vmem'}</td></tr>
            |;
            
            
            if(defined $jobInfo{"Exit_status"} && $jobInfo{"start"} && $jobInfo{"end"}) {
                my $startDate = strftime("%Y-%m-%d %H:%M:%S", localtime($jobInfo{"start"}));
                my $endDate = strftime("%Y-%m-%d %H:%M:%S", localtime($jobInfo{"end"}));
                $jobResourcesStr .= "<tr><td>Duration</td><td>$startDate -> $endDate</td></tr>";
            }
            
            $jobResourcesStr .= qq |
                    <tr><td>Run time</td><td><b>$jobInfo{'resources_used.walltime'}</b>/$jobInfo{'Resource_List.walltime'}</td></tr>
                </table>
            |;
        } elsif($jobPlatform eq 'L') {
            $jobResourcesStr .= qq |
                <table>
                    <tr style='font-weight:bold'><td>Processus PID</td><td>$jobID</td></tr>
            |;
            
            if(!defined $job{infos}{$jobID}) {
                $jobResourcesStr .= qq |
                        <tr style='font-weight:bold'><td>Status</td><td><span style='color:green'>Ended</span></td></tr>
                |;
            } else {
                my %jobInfo = %{$job{infos}{$jobID}};
                my $statusColor = ($jobInfo{status} eq '-') ? 'orange' : 'red';
                $jobInfo{"mem"} = (!$jobInfo{"mem"}) ? '?' :
                                    ($jobInfo{"mem"} > 1024) ?
                                    sprintf("%.2f", ($jobInfo{"mem"}/1024.0))."Mb" :
                                    sprintf("%.2f", $jobInfo{"mem"})."Kb";
                                    
                $jobInfo{"vmem"} = (!$jobInfo{"vmem"}) ? '?' :
                                    ($jobInfo{"vmem"} > 1024) ?
                                    sprintf("%.2f", ($jobInfo{"vmem"}/1024.0))."Mb" :
                                    sprintf("%.2f", $jobInfo{"vmem"})."Kb";
                                    
                $jobInfo{"status"} = 'Running' if($jobInfo{"status"} eq '-');
                
                $jobResourcesStr .= qq |
                        <tr style='font-weight:bold'><td>Name</td><td>$jobInfo{cmd}</td></tr>
                        <tr style='font-weight:bold'><td>Status</td><td><span style='color:$statusColor'>$jobInfo{status}</td></tr>
                        <tr><td>Parent Processus PID</td><td>$jobInfo{ppid}</td></tr>
                        <tr><td>Priority</td><td>$jobInfo{priority}</td></tr>
                        <tr><td>CPU usage</td><td>$jobInfo{cpu} ($jobInfo{pcpu}%)</td></tr>
                        <tr><td>Real + Virtual Memory usage</td><td>$jobInfo{mem} ($jobInfo{pmem}%) + $jobInfo{vmem}</td></tr>
                |;
            }
            
            $jobResourcesStr .= qq |
                </table>
            |;
        }
    }
    
    print($jobResourcesStr);
}

sub printJob {
    my ($jobRef, $currentlyInJob, $userStatus) = @_;
    my %job = %{$jobRef};
    my ($status, $error) = getJobStatus(\%job);
    my (@columns, @values);
    my ($masterJobCode, $jobRank) = split('\.', $job{dirName});
    
    # Set title depending on job type
    my $title = "&bull;&nbsp;";
    if($jobRank) {
        $title .= qq | <span id='deleteMaster:$masterJobCode'></span> |;
    }
    $title .= "&nbsp;Job: $masterJobCode&nbsp;&nbsp;-&nbsp;&nbsp;";
    
    if($jobRank) {
        push(@columns, "Subjob");
        push(@values, $jobRank);
    }
    
    # PhosphoRS
    if($job{type} eq 'Phospho') {
        
        # Build analysis hierarchy
        my $anaID = $job{features}{"ID_ANALYSIS"};
        my @anaInfo=&promsMod::getItemInfo($dbh, 'ANALYSIS', $anaID);
        my @anaHierarchy;
        foreach my $i (1..$#anaInfo) {
            push @anaHierarchy,$anaInfo[$i]{'NAME'};
        }
        my $anaName = join(' > ',@anaHierarchy);

        push(@columns, 'Analysis');
        push(@values, $anaName);
                             
        $title .= "Parameters: Threshold=$job{features}{PROB_THRESHOLD}%, Activation type=$job{features}{ACTIVATION_TYPE}, Mass tolerance=$job{features}{MASS_TOLERANCE} Da";
    }
    
    # All quantification types
    elsif($job{type} =~ /Quantification/) {
        my (@quantificationList, @designList, @experimentList, %designHierarchyStrg);
        my $section='';
        my $quantifType = '';
        $quantifType = $job{features}{TYPE} if($job{features}{TYPE});
        $quantifType .= ":".$job{features}{ALGO_TYPE} if($job{features}{ALGO_TYPE});
        my $quantifName = 'Not available';
        
        if(-e "$job{srcPath}/quantif_info.txt") {
            open (INFO,"$job{srcPath}/quantif_info.txt");
            while (<INFO>) {
                if (/^PARAMETERS:/) {
                    $section='parameters'; next;
                } elsif ($section eq 'parameters' && /^QUANTIF_NAME\t\t(.+)\n/) { # in case not quantifID yet
                    $quantifName = $1 if(!$quantifName);
                }
            }
            close INFO;
        }
    
        my $quantiItemType = (uc($quantifType) !~ /DESIGN|XICMCQ|HYDRO/) ? 'Analysis' : "Quantification";
        my $quantifTypeStr = ($quantifProcesses{$quantifType}) ? $quantifProcesses{$quantifType} : "Unknown ($quantifType)";

        push(@columns, "$quantiItemType name");
        
        if (uc($quantifType) !~ /DESIGN|XICMCQ|HYDRO/) { # MassChroQ=XIC
            #>Analyses quantif
            my $sthQN = $dbh->prepare("SELECT NAME FROM QUANTIFICATION WHERE ID_QUANTIFICATION=?");
            my $anaName = '-';
            if($job{features}{ID_ANALYSIS}) {
                my ($anaID, $parentQuantifID) = split(/\./, $job{features}{ID_ANALYSIS}); # anaID or anaID.parentQuantifID
                my @anaInfo = &promsMod::getItemInfo($dbh, 'ANALYSIS', $anaID);
                my @anaHierarchy;
                foreach my $i (1..$#anaInfo) {
                    push @anaHierarchy,$anaInfo[$i]{'NAME'};
                }
                $anaName = join(' > ',@anaHierarchy);
    
                #>Parent Quantification
                if ($parentQuantifID && uc($quantifType) eq 'XICCORR') { # =~ /SILAC|ITRAQ|TMT/
                    $sthQN->execute($parentQuantifID);
                    my ($parQuantifName) = $sthQN->fetchrow_array;
                    $anaName.="&nbsp;<br/>Using Quantification: $parQuantifName";
                }
                
            }
            
            push(@values, $anaName);
            $sthQN->finish;
        }
        else {
            my ($quantName, @quantHierarchy);
            my $quantID = abs($job{features}{"ID_QUANTIFICATION"});
            my @quantInfo = &promsMod::getItemInfo($dbh, ($job{features}{"ID_QUANTIFICATION"} > 0) ? 'QUANTIFICATION' : 'DESIGN', $quantID);
            foreach my $i (1..$#quantInfo) {
                push @quantHierarchy, $quantInfo[$i]{'NAME'};
            }
            $quantName = join(' > ', @quantHierarchy);
            $quantName .= ' > '.$quantifName if ($job{features}{"ID_QUANTIFICATION"} <= 0 && $quantifName); # no quantifID yet
            
            if($quantName eq 'Not available' && $quantifName) {
                $quantName = $quantifName;
            }
            
            $quantName = substr($quantName, 0, 20)."...".substr($quantName, length($quantName)-60) if(length($quantName) > 85); # TODO optimize quantification name displaying ?
            push(@values, $quantName);
        }
        
        $title .= "Quantification method: $quantifTypeStr";
        
    # All import types (DIA/TDA/TPP/MaxQuant)
    }
    elsif($job{type} =~ /Import|Spectral Library/) {
        my $processType = '';
        
        if($job{type} =~ /Import \[(.+)\]/) {
            $processType = "Import $1 data";
        }
        
        if ($job{type} =~ /MaxQuant/) { # MaxQuant
            my $expID = $job{features}{'ID_EXPERIMENT'};
            my $expName = 'Unknown';
            
            if($expID) {
                ($expName) = $dbh->selectrow_array("SELECT NAME FROM EXPERIMENT WHERE ID_EXPERIMENT=$expID");
            }
            
            push(@columns, "Experiment name");
            push(@values, $expName);
            
            $processType = "Import MaxQuant data";
        }
        elsif ($job{type} =~ /OpenSwath|TDA|TPP|Spectronaut/) { # DIA / TDA import
            my ($libName, $desc) = ('' , '');
            my $software = $job{features}{'SOFTWARE'};
            my $libID = $job{features}{"ID_LIBRARY"};
            my @libIDs = ($libID) ? split(/,/, $libID) : ();
            my $dbID = $job{features}{"ID_DATABANK"};
            my @dbIDs = ($dbID) ? split(/,/, $dbID) : ();
            my $useLib = 0;
            my @ids = (@libIDs) ? @libIDs : (@dbIDs) ? @dbIDs : ();
            if(@ids) {
                $useLib = ($software !~ /Skyline/ && @libIDs) ? 1 : 0; # GET LIBID (DIA) / DBID (TDA)
                foreach my $id (@ids) {
                    my $name;
                    if($useLib) {
                        ($name) = $dbh->selectrow_array("SELECT NAME FROM SWATH_LIB WHERE ID_SWATH_LIB=$id");
                    } else {
                        ($name) = $dbh->selectrow_array("SELECT NAME FROM DATABANK WHERE ID_DATABANK=$id");
                    }
                    
                    if($name) {
                        $libName .= "<br/>" if($libName);
                        $libName .= $name;
                    }
                }
                
                $desc = ($useLib) ? (scalar @ids == 1) ? "Library" : "Libraries" : (scalar @ids == 1) ? "Databank" : "Databanks";
            }
            
            push(@columns, 'Software');
            push(@values, ucfirst(lc($software)));
            
            if($libName) {
                push(@columns, "$desc Name");
                push(@values, $libName);
            }
            
            my $resultFile = $job{features}{"RESULT_FILE"};
            if($job{type} =~ /Spectronaut|OpenSwath/ && $resultFile) {
                push(@columns, 'Result file');
                push(@values, $resultFile);
            } elsif($job{type} =~ /TPP|Spectronaut/ && $job{type} =~ /Spectral Library/) {
                $processType = ($job{features}{"ACTION"} =~ /add/) ? 'Creating spectral library' : 'Updating spectral library';
            }
        }
        
        $title .= $processType;
    }

    # Functional Analysis (only GSEA for now)
    elsif($job{type} =~ /Functional Analysis/) {
        my ($funcAnaID, $parentQuantifID);
        my $funcAnaName = 'Not found';
        my $funcAnaType = '';
        $funcAnaType = $job{features}{TYPE} if ($job{features}{TYPE});

        if ($funcAnaType eq 'GSEA') {
            $funcAnaID = $job{features}{ID_GSEA};
            $parentQuantifID = $job{features}{ID_PARENT_QUANTIF};
            $funcAnaName = $job{features}{GSEA_NAME};
        }
        push(@columns, "$funcAnaType name");
        push(@values, $funcAnaName);

        if ($funcAnaType eq 'GSEA') {
            push(@columns, "Parent quantification");
            my ($quantifFullName, @quantifHierarchy);
            my @parQuantifInfo = &promsMod::getItemInfo($dbh, 'QUANTIFICATION', $parentQuantifID);
            foreach my $i (1..$#parQuantifInfo) {
                push @quantifHierarchy, $parQuantifInfo[$i]{'NAME'};
            }
            $quantifFullName = join(' > ', @quantifHierarchy);
            # TODO optimize name displaying ?
            $quantifFullName = substr($quantifFullName, 0, 20)."...".substr($quantifFullName, length($quantifFullName)-60) if(length($quantifFullName) > 85);
            push(@values, $quantifFullName);
        }

        $title .= "Functional Analysis method: $funcAnaType";
    }
    $title .= "&nbsp;&nbsp;-&nbsp;&nbsp;Launched by user : $job{userName}" if($job{userName});
    
    # Print table header if it is a new job
    if(!$currentlyInJob) {
        print qq |
            <div class='job'>
                <div class="title3">$title</div>
                <TABLE class='jobMeta' cellspacing=0>
                    <thead>
                        <TR bgcolor="$darkColor">
                            <TH class="rbBorder manageCol" style='min-width:176px'>Manage</TH>
        |;
        
        # Print job type-specific columns
        foreach my $column (@columns) {
            print qq |      <TH class="rbBorder">&nbsp;&nbsp;$column&nbsp;&nbsp;</TH> |;
        }
        
        print qq |
        
                            <TH class="rbBorder startedCol">&nbsp;&nbsp;Running time&nbsp;&nbsp;</TH>
                            <TH class="bBorder statusCol">&nbsp;&nbsp;Status&nbsp;&nbsp;</TH>
                        </TR>
                    </thead>
                    <tbody>
        |;
    }
    
    # Print common job information
    my $deleteButtonStr = ($job{"status"} !~ /Done|Error|Stopped/) ? 'Stop' : 'Remove';
    my $act = ($job{"status"} !~ /Done|Error|Stopped/) ? 'stop' : 'delete';
    print qq |
                    <TR bgcolor="$rowColor">
                        <TD style='text-align:center' valign='top'>
                            <button id='delete:$job{ID}' onclick="deleteJob('$job{ID}', '$act')"><!--<i class="fa fa-trash" aria-hidden="true"></i>-->$deleteButtonStr</button>
    |;
    
    if($userStatus =~ /bioinfo/) {
        print qq |
                            <button onclick="viewJobRessources('$job{ID}')"><!--<i class="fa fa-eye" aria-hidden="true"></i>-->Infos</button>
        |;
    }
    
    
    # Print job type-specific information
    foreach my $value (@values) {
        print qq |      <TD nowrap align='center' valign='top'>&nbsp;&nbsp;$value&nbsp;&nbsp;</TD> |;
    }
    
    print qq |          <TD nowrap align='center' valign='top'> |;
    
    # Print last common information
    if($job{endDate} && $job{status} !~ /Queued|Running/) {
        my $startedTime = Time::Piece->strptime("$job{startedDate}", "%Y-%m-%d %H:%M:%S");
        my $endedTime = Time::Piece->strptime("$job{endDate}", "%Y-%m-%d %H:%M:%S");
        my $timeDiff = ($endedTime - $startedTime)->pretty;
        print qq |
                                $job{startedDate} -> $job{endDate}<br/>
                                <b>Duration: $timeDiff</B><br/>
        |;
    } else {
        print qq |              Since <B>$job{startedDate}</B><br/> |;
    }
    
    print qq |          
                            <span id="updated:$job{ID}"></span>
                        </TD>
                        <TD nowrap style="padding:5px">
                            <SPAN id="status:$job{ID}" style='text-align:left; display:block;'>$status</SPAN>
                            <DIV id="error:$job{ID}" style="display:none"><FIELDSET><LEGEND><B>Message :</B></LEGEND>$error</FIELDSET></DIV>
                        </TD>    
                    </TR>
                    <TR id='resources:$job{ID}' style="display:none">
                        <td colspan="5">
                            <fieldset>
                                <legend><B>Job Infos:</B></legend><br/>
                                <button id='realTimeResources:$job{ID}' style="display:none;float: right;position: relative;right: 101px" onclick="viewJobRessources('$jobRef->{ID}', true)">
                                    <!--<i class="fa fa-sync-alt" aria-hidden="true"></i>--> Get real time infos
                                </button>
                                <div id="resourcesInner:$job{ID}" style="width: 100%; max-height: 260px !important; overflow-y: scroll; overflow-x: hidden;"></div>
                            </fieldset><br/>
                        </td>
                    </TR>
    |;
}

####>Revision history<####
# 1.1.16 [ENHANCEMENT] Use $clusterInfo{checkError}->() to filter out irrelevant cluster warnings from job error file (PP 04/06/21)
# 1.1.15 [FEATURE] Add Hydrophobicity jobs (VL 15/02/21)
# 1.1.14 [BUGFIX] Remove tracking of unstarted jobs on cluster that may caused jobs to fail (VS 08/03/21) 
# 1.1.13 [FEATURE] Add display and deletion of Functional Analysis [GSEA] jobs (VL 04/11/20)
# 1.1.12 [MINOR] Remove QUANTIFICATION_NAME from JOB_HISTORY table (VS 09/02/21)
# 1.1.11 [BUGFIX] Fix SQL query to match IDs and names independently from group_concat sort (VL 01/12/20)
# 1.1.10 [ENHANCEMENT] Extend categories to abundance quantifications for job monitoring (VL 19/11/20)
# 1.1.9 [ENHANCEMENT] Improved spectral library building displaying (VS 22/10/20)
# 1.1.8 [MINOR] Changed remaining JOB_HISTORY.STATUS to JOB_HISTORY.JOB_STATUS + Remove typo (VS 09/09/20)
# 1.1.7 [UPDATE] Changed JOB_HISTORY.STATUS to JOB_HISTORY.JOB_STATUS (PP 28/08/20)
# 1.1.6 [ENHANCEMENT] Added project filtering (VS 02/09/20)
# 1.1.5 [ENHANCEMENT] Improved filters displaying (VS 26/08/20)
# 1.1.4 [ENHANCEMENT] Added Spectronaut import job monitoring (VS 06/06/20)
# 1.1.3 [ENHANCEMENT] Uses hash reference for passing extra parameters to &promsQuantif::deleteQuantification (PP 03/06/20)
# 1.1.2 [ENHANCEMENT] Makes job type more readable in filtering section (VS 27/04/20)
# 1.1.1 [BUGFIX] Fix multiple job types filtering (VS 24/04/20)
# 1.1.0 [BUGFIX] Fix processing time computation (VS 07/04/20)
# 1.0.9 [CHANGE] Does not refresh page if any error details are displayed (VS 29/03/20)
# 1.0.8 [BUGFIX] Switch to GET method to submit filter form so that there are no refresh validation popup (VS 29/11/19)
# 1.0.7 [ENHANCEMENT] Add 'user' filter for non-guest users (VS 29/11/19)
# 1.0.6 [MINOR] Do not clean data on job stopping (VS 19/11/19)
# 1.0.5 [ENHANCEMENT] Handles multiple values for filtering parameters (VS 18/11/19)
# 1.0.4 [ENHANCEMENT] Improve reactivity on "stop" and "delete" actions (VS 17/11/19)
# 1.0.3 [ENHANCEMENT] Add "Stopped" job status (VS 14/11/19)
# 1.0.2 [ENHANCEMENT] Declare all JOB_HISTORY FIELDS in &getJobs SQL query (PP 08/11/19)
# 1.0.1 [FEATURE] Add current resources monitoring on the fly (VS 29/10/19)
# 1.0.0 [FEATURE] Merged all monitoring into this script (VS 09/10/19)
