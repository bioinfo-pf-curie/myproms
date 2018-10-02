/*
################################################################################
# linePlot.js   1.0.4                                                          #
# Authors: P. Poullet                                                          #
# Contact: patrick.poullet@curie.fr                                            #
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
*/


/*****************************************************************
                  Data point object
******************************************************************/
function dataPoint(set,data) {
    this.dataSet=set;
    this.label=data[0];
//    this.externalID=(data[1])? data[1] : data[0];
    this.x;
    this.y;
    this.valueList={};
    //this.size=(data[4])? data[4] : 3;
    this.point; // ref to the svg element
    this.pointList; // ref to list svg elements
    this.computeMean('x',data[1]);
    this.computeMean('y',data[2]);
	this.noLine=(data[3])? true : false; // not connected to flanking points
}
dataPoint.prototype = {
    size:10,
    computeMean: function(axis,valuesStrg) {
	var values=valuesStrg.split(':').sort(this.ascNumSort);
	var mean=0;
	for (var i=0; i<values.length; i++) {mean+=(1*values[i]);}
	this[axis]=Math.round(100*mean/values.length)/100;
	if (values.length > 1) {this.valueList[axis]=values;}
    },
    ascNumSort: function(a,b) {return a-b},
    getMinValue: function (axis) {
	return (this.valueList && this.valueList[axis])? this.valueList[axis][0] : this[axis];
    },
    getMaxValue: function (axis) {
	return (this.valueList && this.valueList[axis])? this.valueList[axis][this.valueList[axis].length-1] : this[axis];
    },
    info: function(type) {
		var infoStrg=this.label;
		if (type=='min') {
			return this.label;
		}
		else {
			if (this.dataSet.params.chart.dataSets.length > 1) infoStrg+='\nSet='+this.dataSet.params.name;
			infoStrg+='\nX='+(this.x*1).toFixed(3)+'\nY='+(this.y*1).toFixed(3);
		}
		return infoStrg;
    },
    getX: function() {
	    return this.x;
    },
    getY: function() {
	    return this.y;
    }
}


/*****************************************************************
                  Line Plot object
******************************************************************/
function linePlot(lineData) {
    var LP=this;
	this.chartID=registerInCuBioJS(this);
    this.divID=lineData.div;
    this.div=document.getElementById(this.divID);


	/******** Chart geometry ********/
    var TOP_SPACE=20;
    var BOTTOM_SPACE=20;
    var LEFT_SPACE=0;
    var RIGHT_SPACE=20;
    var XAXIS_SPACE=30;
    var YAXIS_SPACE=60;
    var chartX=LEFT_SPACE+YAXIS_SPACE+1;
    var chartY=TOP_SPACE+1;
    var chartW=(lineData.width)? lineData.width : 400;
    var chartH=(lineData.height)? lineData.height : 400;
	var canvasW=LEFT_SPACE+YAXIS_SPACE+chartW+RIGHT_SPACE;
    var canvasH=TOP_SPACE+chartH+XAXIS_SPACE+BOTTOM_SPACE;
	this.formPosition=(lineData.legendPosition)? lineData.legendPosition : 'right';


	/******** Chart variables & objects ********/
    this.dataSets=new Array();
	//this.dataSetLinks=new Object(); // List of links between dataPoint with same id in different dataSets
	this.selectAllDataSets=false;
	this.allowHighlight=false; // flag to allow or not highlighting
	//this.updateHighlight=lineData.updateHighlight; // object {callback:callback function in case delete/edit,editable:true/false (hl name edition flag)}
	//this.highlightedPoints=new Object(); // List of different user-defined labeled set of points
    //var colorList=['#000000','#4AA02C','#F660AB','#FBB917','#0000FF','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B'];
	var colorList=['#0000FF','#4AA02C','#000000','#F660AB','#FBB917','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B'];
	var colorIdx=0;
	this.axisXtext=lineData.axisX;
	this.axisYtext=lineData.axisY;
	this.axisXtitle=null; // SVG object
	this.axisYtitle=null; // SVG object
	this.forceXto0=0; this.forceYto0=0; // 0 or null: auto, 1: 0, 2: ~0
	this.sameScale=0;
	this.zoomable=false;
	this.showDataSets=true;
	this.searchable=false;
	this.selectable=false;
	this.customPointLabel=(lineData.pointLabel)? lineData.pointLabel : null;
	this.onPointExclusion=(lineData.onPointExclusion)? lineData.onPointExclusion : null;
	this.connectPoints=true;

    this.minValueX,this.maxValueX,this.minValueY,this.maxValueY;
	this.chartSettings={};
	this.chartSettings.reference={};
	this.chartSettings.current={};
	this.axisClosure=true; // draw 4 axes
	this.chartMarks=[];
	this.plotArea;
/*
	this.zoomText;
	this.dragArea;
	this.dragContext='area';
	this.panProcess={x:0,y:0,dx:0,dy:0}; //new Object();
	this.pointOnclick=lineData.pointOnclick;
	this.pointOnList=lineData.pointOnList;
*/

    /********************* Data import *********************/
    /***** Threshold lines *****/
/*
	this.editThreshold=false;
	this.thresholdLines=new Array();
	this.thresholdLines.push(new thresholdLine(LP,'X','',0,'#555555','- '));
	this.thresholdLines.push(new thresholdLine(LP,'Y','',0,'#555555','- '));
*/

    /***** Datasets *****/
    this.addDataSet=function(dsIdx,set) {
		this.dataSets[dsIdx]={};
		this.dataSets[dsIdx].params={};
		this.dataSets[dsIdx].params.chart=LP;
		this.dataSets[dsIdx].params.index=dsIdx;
		this.dataSets[dsIdx].params.name=(set.name)? set.name : 'dataSet '+dsIdx;
		this.dataSets[dsIdx].params.color=(set.color)? set.color : colorList[colorIdx++];
		this.dataSets[dsIdx].params.visible=true;
		this.dataSets[dsIdx].data=[];
		this.dataSets[dsIdx].line;
		if (colorIdx==colorList.length) colorIdx=0;
    }

    /***** Data points *****/
    this.addDataAsString=function(dsIdx,dataStrg) {
		var dataSet=dataStrg.split(';');
		for (var i=0; i < dataSet.length; i++) {
			var data=dataSet[i].split(',');
			this.dataSets[dsIdx].data.push(new dataPoint(this.dataSets[dsIdx],data));
			var lastPoint=this.dataSets[dsIdx].data[this.dataSets[dsIdx].data.length-1];
			if (this.minValueX==null) { // 1st point
				this.minValueX=lastPoint.getMinValue('x');
				this.maxValueX=lastPoint.getMaxValue('x');
				this.minValueY=lastPoint.getMinValue('y');
				this.maxValueY=lastPoint.getMaxValue('y');
			}
			else {
				this.minValueX=Math.min(this.minValueX,lastPoint.getMinValue('x'));
				this.maxValueX=Math.max(this.maxValueX,lastPoint.getMaxValue('x'));
				this.minValueY=Math.min(this.minValueY,lastPoint.getMinValue('y'));
				this.maxValueY=Math.max(this.maxValueY,lastPoint.getMaxValue('y'));
			}
//console.log(lastPoint.getX()+','+lastPoint.getY());
		}
    }

    /********************* LP display (Raphael) *********************/
    this.draw=function() {
		/* Cross-dataset links (links between dataPoint with same externalID) */
/*
		if (this.dataSets.length > 1) {
			for (var i=0; i < this.dataSets.length; i++) {
				for (var j=0; j < this.dataSets[i].data.length; j++) {
					var dataPointID=this.dataSets[i].data[j].externalID;
					if (!this.dataSetLinks[dataPointID]) {this.dataSetLinks[dataPointID]=new Array();}
					this.dataSetLinks[dataPointID].push(this.dataSets[i].data[j]);
				}
			}
		}
*/
		/* DIVs */
		var canvasDivID=this.divID+'_canvas';
		var formDivID=this.divID+'_form';
		//this.div.innerHTML='<TABLE><TR><TD style="min-width:'+canvasW+'px"><DIV id="'+canvasDivID+'"></DIV></TD><TD valign="top" style="min-width:200px"><DIV id="'+formDivID+'"></DIV></TD></TR></TABLE>';
		drawChartBorders(LP,canvasW,canvasDivID,formDivID);
		this.formDiv=document.getElementById(formDivID);
		/* Form menu */
		initializeForm(LP);
		/* Canvas */
		this.canvas=Raphael(canvasDivID,canvasW,canvasH);
		/* Back panel */
		this.backPanel=this.canvas.rect(0,0,canvasW,canvasH,0).attr({fill:'#fff',stroke:'#000'});
		/* Chart area */
		this.plotArea=this.canvas.rect(chartX,chartY,chartW,chartH,0)
		.attr({fill:'#F3F3F3',stroke:'none'}); //,cursor:"crosshair"
		/* Display chart */
		initializeChart(LP);
    }

    /***** Reset data values (x,y values) *****/
/*
    var tmpDataList=new Object();
    this.resetData=function(dataStrg) {
	    var dataSet=dataStrg.split(';');
	    for (var i=0; i < dataSet.length; i++) {
		    var data=dataSet[i].split(','); // [externalID,newX,newY]
		    tmpDataList[data[0]]=[data[1],data[2]];
	    }
    }
*/



}

/*
####>Revision history<####
# 1.0.4 Uses chart registering for dynamic html calls (PP 20/02/16)
# 1.0.3 GPL license (PP 23/09/13)
# 1.0.2 Commented LP.dataSetLinks declaration (PP 29/01/13)
# 1.0.1 Added noLine attribute to dataPoint for exclusion switch (PP 02/12/12)
# 1.0.0 Production version & form position option (PP 30/05/12)
# 0.0.1 starting script (PP 12/05/12)
*/
