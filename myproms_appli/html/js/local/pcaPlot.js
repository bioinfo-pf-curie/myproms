/*
################################################################################
# pcaPlot.js      1.1.1                                                        #
# Authors: P. Poullet (Institut Curie)                                         #
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
//function dataPoint(set,data) {
////alert(data.length);
//    this.dataSet=set;
//    this.label=data[0];
//    this.externalID=(data[1])? data[1] : data[0];
//    this.x=data[2];
//    this.y=data[3];
//    this.size=(data[4])? data[4] : 5;
//    this.point=null; // ref to the svg element
//	this.highlightNames=[]; // new Array()
//	this.isHidden=false;
//}
//dataPoint.prototype = {
//    info: function(type) {
//		var infoStrg=this.label;
//		if (type=='min') {
//			return this.label;
//		}
//		else {
//			if (this.dataSet.params.chart.dataSets.length > 1) infoStrg+='\nSet='+this.dataSet.params.name;
//			infoStrg+='\nX='+(this.x*1).toFixed(3)+'\nY='+(this.y*1).toFixed(3);
//		}
//		return infoStrg;
//    },
//    getX: function() {
//		return this.x;
//	},
//	getY: function() {
//		return this.y;
//    },
//	setX: function(newValue) {
//		this.x=newValue;
//	},
//	setY: function(newValue) {
//		this.y=newValue;
//	}
//}


/*****************************************************************
                  PCA Plot object
******************************************************************/
function pcaPlot(pca) {
    var PCA=this;
	this.chartID=registerInCuBioJS(this);
    this.divID=pca.div;
	this.exportAsImage=pca.exportAsImage;
	var canvasDivID; // initialized by this.draw()
	this.getDivID=function(){return canvasDivID;};
    this.div=document.getElementById(this.divID);

	/******** Chart geometry ********/
    var TOP_SPACE=20;
    var BOTTOM_SPACE=20;
    var LEFT_SPACE=10;
    var RIGHT_SPACE=20;
    var XAXIS_SPACE=30;
    var YAXIS_SPACE=50;
    var chartX=LEFT_SPACE+YAXIS_SPACE+1;
    var chartY=TOP_SPACE+1;
    var chartW=(pca.width)? pca.width : 400;
    var chartH=(pca.height)? pca.height : 400;
	var canvasW=LEFT_SPACE+YAXIS_SPACE+chartW+RIGHT_SPACE;
    var canvasH=TOP_SPACE+chartH+XAXIS_SPACE+BOTTOM_SPACE;


	/******** Chart variables & objects ********/
    this.dataSets=[];
	this.dataSetLinks={}; // List of links between dataPoint with same id in different dataSets
	this.selectAllDataSets=false;
	this.allowHighlight=pca.allowHighlight; // flag to allow or not highlighting
	this.updateHighlight=pca.updateHighlight; // object {callback:callback function in case delete/edit,editable:true/false (hl name edition flag)}
	this.highlightedPoints={}; // List of different user-defined labeled set of points
	this.hideUnmatchedPoints=false; // flag to hide/show points not matching active highlighting
    //var colorList=['#000000','#4AA02C','#F660AB','#FBB917','#0000FF','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B'];
	//var colorList=['#0000FF','#4AA02C','#000000','#F660AB','#FBB917','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B'];
	this.axisXtext=pca.axisX;
	this.axisYtext=pca.axisY;
	this.axisXtitle={}; // must be set to something for chartLibrary to reset with text SVG object
	this.axisYtitle={}; // must be set to something for chartLibrary to reset with text SVG object
	this.forceXto0=0; this.forceYto0=0; // 0 or null: auto, 1: 0, 2: ~0
	this.sameScale=1;
 	this.convertValue=function(axis,value) {
		return value;
	};
    this.minValueX=null;this.maxValueX=null;this.minValueY=null;this.maxValueY=null;
	this.chartSettings={};
	this.chartSettings.reference={};
	this.chartSettings.current={};
	this.axisClosure=true; // draw 4 axes
	this.chartMarks=[];
	this.zoomable=true;
	this.zoomText=null;
	this.plotArea=null;
	this.dragArea=null;
	this.dragContext='area';
	this.panProcess={x:0,y:0,dx:0,dy:0}; //new Object();
	this.pointOpacity=(pca.pointOpacity)? pca.pointOpacity : 0.7;
	this.pointOnclick=pca.pointOnclick;
	this.pointOnList=pca.pointOnList;
	this.selectable=true;
	this.searchable=pca.searchable || true;

	/****************** Data point object ******************************/
	var dataPoint=function(set,data) {
		this.dataSet=set;
		this.label=data[0];
		this.externalID=(data[1])? data[1] : data[0];
		this.x=data[2];
		this.y=data[3];
		this.z=(data[4])? data[4] : 0;
		this.size=15; // default
		this.point=null; // ref to the svg element
		this.highlightNames=[]; // new Array()
		this.isHidden=false;
	};
	dataPoint.prototype = {
		computeSize: function(refMinSize,minZ,scale) {
			//var r=(this.z-minZ)*scale;
			//this.size=3.14*r*r;
			this.size=refMinSize+((this.z-minZ)*(this.z-minZ)*3.14*scale);
		},
		info: function(type) {
			var infoStrg=this.label;
			if (type=='min') {
				return this.label;
			}
			else {
				if (this.dataSet.params.chart.dataSets.length > 1) infoStrg+='\nSet='+this.dataSet.params.name;
				infoStrg+='\nX='+(this.x*1).toFixed(3)+'\nY='+(this.y*1).toFixed(3);
				if (this.z) {
					infoStrg+='\nZ='+(this.z*1).toFixed(3);
				}
			}
			return infoStrg;
		},
		getX: function() { // needed by chartLibrary2.js
			return this.x;
		},
		getY: function() {
			return this.y;
		}
		//,
		//setX: function(newValue) {
		//	this.x=newValue;
		//},
		//setY: function(newValue) {
		//	this.y=newValue;
		//}
	};

    /********************* Data import *********************/
    /***** Threshold lines *****/
	this.editThreshold=false;
	this.thresholdLines=[];
	this.thresholdLines.push(new thresholdLine(PCA,'X','',0,'#555555','- '));
	this.thresholdLines.push(new thresholdLine(PCA,'Y','',0,'#555555','- '));


    /***** Datasets (only 1)*****/
	this.dataSets[0]={};
	this.dataSets[0].params={};
	this.dataSets[0].params.chart=PCA;
	this.dataSets[0].params.index=0;
	this.dataSets[0].params.name='PCA data';
	this.dataSets[0].params.color='#000000';
	this.dataSets[0].params.visible=true;
	this.dataSets[0].params.sizeRule={min:15,max:200,ratio:1};
	this.dataSets[0].data=[];

    /***** Data points *****/
	var minZ=1e200,maxZ=-1e-200;
    this.addDataAsString=function(dataStrg) {
		if (!this.dataSets[0]) {
			this.addDataSet(0,{name:'PCA data',color:'#000000'});
			this.dataSets[0].data=[];
		}
		var dataSet=dataStrg.split(';');
		for (var i=0; i < dataSet.length; i++) {
			var data=dataSet[i].split(',');
			var dp=new dataPoint(this.dataSets[0],data);
			this.dataSets[0].data.push(dp);
			if (this.minValueX===null) { // 1st point
				this.minValueX=dp.x;
				this.maxValueX=dp.x;
				this.minValueY=dp.y;
				this.maxValueY=dp.y;
			}
			else {
				this.minValueX=Math.min(this.minValueX,dp.x);
				this.maxValueX=Math.max(this.maxValueX,dp.x);
				this.minValueY=Math.min(this.minValueY,dp.y);
				this.maxValueY=Math.max(this.maxValueY,dp.y);
			}
			minZ=Math.min(minZ,dp.z);
			maxZ=Math.max(maxZ,dp.z);
		}
    };

    /********************* PCA display (Raphael) *********************/
    this.draw=function() {
		/* Compute point size based on z axis */
		var refMinSize=35,minZsize=0,scaleZ=1; // default (no z data)
		if (minZ < maxZ) { // z data
			var sizeRule=this.dataSets[0].params.sizeRule;
			refMinSize=sizeRule.min;
			minZsize=minZ;
			//var scaleZ=8/(maxZ-minZ); // max r=8 => max size ~200; min/max size overwritten by sizeRule
			scaleZ=(sizeRule.max-sizeRule.min)/((maxZ-minZ)*(maxZ-minZ)*3.14); // radius -> surface
		}
		for (let j=0; j < this.dataSets[0].data.length; j++) {
			//this.dataSets[0].data[j].computeSize(minZ,scaleZ);
			this.dataSets[0].data[j].computeSize(refMinSize,minZsize,scaleZ);
		}
		minZ=1e200; maxZ=-1e-200; // reset in case resetData() & redraw()
		/* Cross-dataset links (links between dataPoint with same externalID) */
		if (this.dataSets.length > 1) {
			for (let i=0; i < this.dataSets.length; i++) {
				for (let j=0; j < this.dataSets[i].data.length; j++) {
					var dataPointID=this.dataSets[i].data[j].externalID;
					if (!this.dataSetLinks[dataPointID]) {this.dataSetLinks[dataPointID]=[];}
					this.dataSetLinks[dataPointID].push(this.dataSets[i].data[j]);
				}
			}
		}
		/* DIVs */
		canvasDivID=this.divID+'_canvas';
		var formDivID=this.divID+'_form';
		this.div.innerHTML='<TABLE><TR><TD style="min-width:'+canvasW+'px"><DIV id="'+canvasDivID+'"></DIV></TD><TD valign="top" style="min-width:200px"><DIV id="'+formDivID+'"></DIV></TD></TR></TABLE>';
		this.formDiv=document.getElementById(formDivID);
		/* Form menu */
		initializeForm(PCA);
		/* Canvas */
		this.canvas=Raphael(canvasDivID,canvasW,canvasH);
		/* Back panel */
		this.canvas.rect(0,0,canvasW,canvasH,0).attr({fill:'#fff',stroke:'#000'});
		/* Chart area */
		this.plotArea=this.canvas.rect(chartX,chartY,chartW,chartH,0)
		.attr({fill:'#F3F3F3',stroke:'none'}) //,cursor:"crosshair"
		.drag(function(dx,dy,x,y,e) {extendDragging(PCA,dx,dy,x,y,e);},
			  function(x,y,e) {setDragging(PCA,x,y,e);},
			  function() {endDragging(PCA);}
			  )
		.dblclick(function(){zoomOut(PCA);});

		/* Display chart */
		initializeChart(PCA);
    };

	/***** Reset data values (x,y values) *****/
	var tmpDataList={};
	this.resetData=function(dataStrg) {
		var dataSet=dataStrg.split(';');
		for (var i=0; i < dataSet.length; i++) {
			var data=dataSet[i].split(','); // [externalID,newX,newY]
			tmpDataList[data[0]]=[data[1],data[2],data[3]];
			if (data[3]) {
				minZ=Math.min(minZ,data[3]);
				maxZ=Math.max(maxZ,data[3]);
			}
		}
	};

	/***** Redraw graph with new values *****/
	this.redraw=function(titleX,titleY) {
		//var chartX=this.plotArea.attr('x'); var chartX0=chartX-1;
		//var chartY=this.plotArea.attr('y');
		//var chartW=this.plotArea.attr('width');
		//var chartH=this.plotArea.attr('height');

		this.axisXtext=titleX;
		this.axisYtext=titleY;
		this.axisXtitle.attr({text:this.axisXtext});
		this.axisYtitle.attr({text:this.axisYtext});

		/*** Update x,y,z,size values ***/
		var sizeRule=PCA.dataSets[0].params.sizeRule;
		var scaleZ;
		if (minZ < maxZ) { // z data
			//scaleZ=8/(maxZ-minZ); // max r=8 => max size ~200; min size overwritten by sizeRule to 2
			scaleZ=(sizeRule.max-sizeRule.min)/((maxZ-minZ)*(maxZ-minZ)*3.14); // radius -> surface
		}
		var numMatches=0;
		for (var i=0; i < this.dataSets[0].data.length; i++) {
			var dp=this.dataSets[0].data[i];
			//dp.z=0; // default
			//dp.size=15; // default
			if (!tmpDataList[dp.externalID]) continue;
			numMatches++;
			dp.x=tmpDataList[dp.externalID][0];
			dp.y=tmpDataList[dp.externalID][1];
			dp.z=(tmpDataList[dp.externalID][2])? tmpDataList[dp.externalID][2] : 0;
			if (numMatches==1) { // normally <=> i=0 (just to be safe)
				this.minValueX=dp.x;
				this.maxValueX=dp.x;
				this.minValueY=dp.y;
				this.maxValueY=dp.y;
			}
			else {
				this.minValueX=Math.min(this.minValueX,dp.x);
				this.maxValueX=Math.max(this.maxValueX,dp.x);
				this.minValueY=Math.min(this.minValueY,dp.y);
				this.maxValueY=Math.max(this.maxValueY,dp.y);
			}
			if (minZ < maxZ) { // z data
				dp.computeSize(sizeRule.min,minZ,scaleZ);
			}
			else {
				dp.computeSize(35,0,0); // default size when no z data
			}
			//var usedSize=Math.max(Math.min(dp.size,sizeRule.max),sizeRule.min);
			//dp.point.attr({r:sizeRule.ratio*Math.sqrt(usedSize/3.14)});
			dp.point.attr({r:sizeRule.ratio*Math.sqrt(dp.size/3.14)});
		}
		tmpDataList={}; // clear list
		minZ=1e200; maxZ=-1e-200; // reset in case new resetData() & redraw()

		/*** Compute default range ***/
		computeDefaultRange(PCA);

		/*** Clear & plot chart ***/
		clearChart(PCA);
		plotChart(PCA);
	};

} // end of PCA

/*
####>Revision history<####
# 1.1.1 [FEATURE] Optional external pre-search function (PP 29/10/19)
# 1.1.0 Improved computation of point size based on 3rd dimension (PP 02/02/17)
# 1.0.9 Uses chart registering for dynamic html calls (PP 20/02/16)
# 1.0.8 Added 3rd dimension projection as point size (PP 08/11/15)
# 1.0.7 Added getDivID function && optional export image button (PP 08/08/14)
# 1.0.6 Minor improvements in dimension-change redraw (PP 02/07/14)
# 1.0.5 Compatibility with hide/show points not matching highlight (PP 14/06/14)
# 1.0.4 Minor compatibility & simplication changes (PP 17/03/13)
# 1.0.3 Handles properly dataset with no data (PP 26/07/12)
# 1.0.2 Compatibility with chartLibrary 1.0.2 (PP 13/05/12)
# 1.0.0 Production version (PP 04/04/12)
*/
