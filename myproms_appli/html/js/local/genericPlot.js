/*
################################################################################
# genericPlot.js    1.0.8                                                      #
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
function gpDataPoint(set,data) {
//alert(data.length);
    this.dataSet=set;
    this.label=data[0];
    this.externalID=(data[1])? data[1] : data[0];
    this.x=null;
    this.y=null;
    this.valueList={};
    this.pointList=null; // ref to list svg elements
    this.computeMean('x',data[2]);
    this.computeMean('y',data[3]);
    this.size=(data[4])? data[4] : 10;
    this.point=null; // ref to the svg element
	this.highlightNames=[]; // new Array()
}
gpDataPoint.prototype = {
	computeMean: function(axis,valuesStrg) {
//console.log(axis+','+valuesStrg);
		var values=valuesStrg.split(':').sort(this.ascNumSort);
		var mean=0;
		for (var i=0; i<values.length; i++) {mean+=(1*values[i]);}
		this[axis]=Math.round(100*mean/values.length)/100;
		if (values.length > 1) {this.valueList[axis]=values;}
    },
	ascNumSort: function(a,b) {return a-b;},
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
                  Generic Plot object
******************************************************************/
function genericPlot(plotData) {
    var GP=this;
	this.chartID=registerInCuBioJS(this);

    this.divID=plotData.div;
    this.div=document.getElementById(this.divID);
	var canvasDivID=this.divID+'_canvas';
	var formDivID=this.divID+'_form';
	this.getDivID=function(){return canvasDivID;};
	this.exportAsImage=plotData.exportAsImage; // null or array [button name, image file name, action script path]

	var colorList=['#0000FF','#4AA02C','#000000','#F660AB','#FBB917','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B'];

    /******** Chart variables & objects ********/
    this.dataSets=[];
	this.dataSetLinks={}; // List of links between dataPoint with same id in different dataSets
	this.selectAllDataSets=false;
	this.thresholdLines=[];
	this.editThreshold=plotData.editThreshold; // if true: all th are editable
	this.existEditableThreshold=false; // default (al least 1 editable th)
	this.allowHighlight=plotData.allowHighlight; // flag to allow or not highlighting
	this.updateHighlight=plotData.updateHighlight; // object {callback:callback function in case delete/edit,editable:true/false (hl name edition flag)}
	this.highlightedPoints={}; // List of different user-defined labeled set of points
	this.pointOnclick=plotData.pointOnclick;
	this.pointOnList=plotData.pointOnList;
    //var pointDefaultSize=(plotData.pointDefaultSize)? plotData.pointDefaultSize : 10;
	this.pointOpacity=(plotData.pointOpacity)? plotData.pointOpacity : 1; // point opacity
	this.zoomable=plotData.zoomable;
	this.showDataSets=(plotData.hideDataSets)? false : undefined;
	this.searchable=(plotData.searchable)? plotData.searchable : (plotData.noSearch)? false : true;
	this.selectable=(plotData.noSelect)? false : true;
	this.customPointLabel=(plotData.pointLabel)? plotData.pointLabel : null;
	this.connectPoints=plotData.connectpoint; // undef or line or curve
	this.onPointExclusion=(plotData.onPointExclusion)? plotData.onPointExclusion : null;
	this.convertValue=(plotData.convertValue)? plotData.convertValue : function(axis,value) {return value;};
	this.sameScale=(plotData.sameScale)? true : false;
	this.chartSettings={reference:{},current:{}};
	this.chartMarks=[];
	this.plotArea=null;
	this.dragArea=null;
	this.dragContext='area';
	this.zoomText=null;
	this.panProcess={x:0,y:0,dx:0,dy:0}; //new Object();

	/***** Chart axes parameters *****/
	this.axisClosure=plotData.axisClosure; // draw 4 axes
	var axes=['X','Y','X2','Y2'];
	for (var i=0; i<axes.length; i++) {
		if (plotData['axis'+axes[i]]) {
			this['axis'+axes[i]+'text']=plotData['axis'+axes[i]].title;
			this['force'+axes[i]+'to0']=(plotData['axis'+axes[i]].forceTo0)? plotData['axis'+axes[i]].forceTo0 : 0; // 0 or null: auto, 1: 0, 2: ~0
			this['axis'+axes[i]+'title']=null; // SVG object
			this['minValue'+axes[i]]=null;
			this['maxValue'+axes[i]]=null;
			if (plotData['axis'+axes[i]].zoomable) this.zoomable=true; // overwritten
		}
	}
/*
	this.axisXtext=plotData.axisX;
	this.axisX2text=plotData.axisX2;
	this.axisYtext=plotData.axisY;
	this.axisY2text=plotData.axisY2;
	this.forceXto0=(plotData.forceXto0)? plotData.forceXto0 : 0; // 0 or null: auto, 1: 0, 2: ~0
	this.forceX2to0=(plotData.forceX2to0)? plotData.forceX2to0 : 0;
	this.forceYto0=(plotData.forceYto0)? plotData.forceYto0 : 0;
	this.forceY2to0=(plotData.forceY2to0)? plotData.forceY2to0 : 0;
	this.axisXtitle=null; // SVG object
	this.axisX2title=null; // SVG object
	this.axisYtitle=null; // SVG object
    this.minValueX=this.maxValueX=this.minValueY=this.maxValueY=null;
	this.minValueX2=this.maxValueX2=this.minValueY2=this.maxValueY2=null;
*/

	/******** Chart geometry ********/
    var BORDER_SPACE=10;
	var ZOOM_SPACE=10;
    var XAXIS_SPACE=40;
    var YAXIS_SPACE=50;
    var chartW=(plotData.width)? plotData.width : 400;
    var chartH=(plotData.height)? plotData.height : 400;
    var chartX=BORDER_SPACE+1;
	if (this.axisYtext) chartX+=YAXIS_SPACE;
    var chartY=BORDER_SPACE+1;
	if (this.zoomable) chartY+=ZOOM_SPACE;
	if (this.axisX2text) chartY+=XAXIS_SPACE;
	var canvasW=chartX+chartW+BORDER_SPACE-1;
	if (this.axisY2text) canvasW+=YAXIS_SPACE;
    var canvasH=chartY+chartH+BORDER_SPACE-1;
	if (this.axisXtext) canvasH+=XAXIS_SPACE;

    /********************* Data import *********************/
	/***** Threshold lines *****/
	this.addThreshold=function(th) { // [axis,name,value,color,dashString,keepAbove1,roundValue]
		var editable=false;
		if (th.editable) {
			this.existEditableThreshold=true;
			editable=true;
		}
		var TH=new thresholdLine(GP,th.axis,th.label,th.value,th.color,th.pattern,th.keepAbove1,th.roundValue,editable);
		this.thresholdLines.push(TH);
		if (this['minValue'+TH.axis]===null) {this['minValue'+th.axis]=this['maxValue'+TH.axis]=TH.value;}
		else {
			this['minValue'+TH.axis]=Math.min(this['minValue'+TH.axis],TH.value);
			this['maxValue'+TH.axis]=Math.max(this['maxValue'+TH.axis],TH.value);
		}
	};

    /***** Datasets *****/
    this.addDataSet=function(dsIdx,set) {
		this.dataSets[dsIdx]={};
		this.dataSets[dsIdx].params={};
		this.dataSets[dsIdx].params.chart=GP;
		this.dataSets[dsIdx].params.index=dsIdx;
		this.dataSets[dsIdx].params.axisX=(set.axisX && set.axisX.match('2'))? 'X2' : 'X'; // dual X axes
		this.dataSets[dsIdx].params.axisY=(set.axisY && set.axisY.match('2'))? 'Y2' : 'Y'; // dual Y axes
		this.dataSets[dsIdx].params.name=(set.name)? set.name : 'dataSet '+dsIdx;
		this.dataSets[dsIdx].params.color=(set.color)? set.color : colorList[dsIdx % colorList.length];
		this.dataSets[dsIdx].params.sizeName=(set.sizeName)? set.sizeName : null;
		this.dataSets[dsIdx].params.sizeRule=(set.sizeRule)? set.sizeRule : null;
		this.dataSets[dsIdx].params.line=(set.line)? set.line : null;
		this.dataSets[dsIdx].params.visible=true;
		this.dataSets[dsIdx].data=[];
		this.dataSets[dsIdx].line=null; // SVG path object
    };

    /***** Data points *****/
	this.addDataAsString=function(dsIdx,dataStrg) {
		var dataSet=dataStrg.split(';');
		for (var i=0; i < dataSet.length; i++) {
			var minValueXn='minValue'+this.dataSets[dsIdx].params.axisX;
			var maxValueXn='maxValue'+this.dataSets[dsIdx].params.axisX;
			var minValueYn='minValue'+this.dataSets[dsIdx].params.axisY;
			var maxValueYn='maxValue'+this.dataSets[dsIdx].params.axisY;
			var data=dataSet[i].split(',');
			var dp=new gpDataPoint(this.dataSets[dsIdx],data);
			this.dataSets[dsIdx].data.push(dp);
			var minX=dp.getMinValue('x');
			var maxX=dp.getMaxValue('x');
			var minY=dp.getMinValue('y');
			var maxY=dp.getMaxValue('y');
			if (this[minValueXn]==null) { // 1st point
				this[minValueXn]=minX;
				this[maxValueXn]=maxX;
			}
			else {
				this[minValueXn]=Math.min(this[minValueXn],minX);
				this[maxValueXn]=Math.max(this[maxValueXn],maxX);
			}
			if (this[minValueYn]==null) { // 1st point
				this[minValueYn]=minY;
				this[maxValueYn]=maxY;
			}
			else {
				this[minValueYn]=Math.min(this[minValueYn],minY);
				this[maxValueYn]=Math.max(this[maxValueYn],maxY);
			}
		}
    };

	/********************* Plot display (Raphael) *********************/
    this.draw=function() {
		/* DIVs */
		//this.div.innerHTML='<TABLE><TR><TD><DIV id="'+canvasDivID+'"></DIV></TD><TD valign="top"><DIV id="'+formDivID+'"></DIV></TD></TR></TABLE>';
		this.div.innerHTML='<TABLE><TR><TD style="min-width:'+canvasW+'px"><DIV id="'+canvasDivID+'"></DIV></TD><TD valign="top" style="min-width:200px"><DIV id="'+formDivID+'"></DIV></TD></TR></TABLE>';
		this.formDiv=document.getElementById(formDivID);
		/* Form menu */
initializeForm(GP);
		/* Canvas */
		this.canvas=Raphael(canvasDivID,canvasW,canvasH);
		/* Back panel */
		this.backPanel=this.canvas.rect(0,0,canvasW,canvasH,0).attr({fill:'#fff',stroke:'#000'});
		/* Chart area */
		this.plotArea=this.canvas.rect(chartX,chartY,chartW,chartH,0)
		.attr({fill:'#F3F3F3',stroke:'none'}) //,cursor:"crosshair"
		.drag(function(dx,dy,x,y,e) {extendDragging(GP,dx,dy,x,y,e);},
			  function(x,y,e) {setDragging(GP,x,y,e);},
			  function() {endDragging(GP);}
			  );
		if (this.zoomable) this.plotArea.dblclick(function(){zoomOut(GP)});


		/***** Cross-dataset links (links between dataPoint with same externalID) *****/
		if (this.dataSets.length > 1) {
			for (var i=0; i < this.dataSets.length; i++) {
				for (var j=0; j < this.dataSets[i].data.length; j++) {
					var dataPointID=this.dataSets[i].data[j].externalID;
					if (!this.dataSetLinks[dataPointID]) {this.dataSetLinks[dataPointID]=[];}
					this.dataSetLinks[dataPointID].push(this.dataSets[i].data[j]);
				}
			}
		}


		/***** Display chart(s) *****/
		initializeChart(GP);

    }

} // end of GP

/*
####>Revision history<####
# 1.0.8 [FEATURE] Optional external pre-search function (PP 29/10/19)
# 1.0.7 Minor bug fix in showDataSets initialization (PP 26/10/17)
# 1.0.6 Uses chart registering for dynamic html calls (PP 20/02/16)
# 1.0.5 Renamed dataPoint object to gpDataPoint for compatiblity with volcanoPlot library (PP 28/10/15)
# 1.0.4 Handle dataSet.params.sizeRule (PP 21/10/15)
# 1.0.3 PNG export option (PP 13/02/15)
# 1.0.2 Same scale on all axes & bug fix in threshold in chart range (PP 08/04/14)
# 1.0.1 GPL license (PP 23/09/13)
# 0.0.2 Added point opacity management (PP 15/07/13)
# 0.0.1 Started from volcanoPlot.js v.1.1.0 (PP 17/03/13)
*/
