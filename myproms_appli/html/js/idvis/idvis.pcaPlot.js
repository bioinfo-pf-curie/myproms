/*
################################################################################
# idvis.pcaPlot.js      2.0.3                                                  #
# Authors: P. Poullet (Institut Curie)                                         #
# Contact: patrick.poullet@curie.fr                                            #
################################################################################
#----------------------------------CeCILL License-------------------------------
# This file is part of idvis (Interactive Data VISualization)
# idvis is a set of JavaScript libraries used as a layer above RaphaëlJS (http://raphaeljs.com) to draw interactive SVG images
# Copyright Institut Curie 2019
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
                  PCA Plot object
******************************************************************/
idvis.pcaPlot = function(pcaData) {
    var PCA=this;
    this.divID=pcaData.div;
    this.div=document.getElementById(this.divID);
    this.chartID=idvis.lib.registerChart(this);
	this.exportAsImage=pcaData.exportAsImage;
	var canvasDivID; // initialized by this.draw()
	this.getDivID=function(){return canvasDivID;};

	/******** Chart geometry ********/
    const TOP_SPACE=20,
		  BOTTOM_SPACE=20,
		  LEFT_SPACE=10,
		  RIGHT_SPACE=20,
		  XAXIS_SPACE=30,
		  YAXIS_SPACE=50,
		  chartX=LEFT_SPACE+YAXIS_SPACE+1,
		  chartY=TOP_SPACE+1,
		  chartW=pcaData.width || 400,
		  chartH=pcaData.height || 400,
		  canvasW=LEFT_SPACE+YAXIS_SPACE+chartW+RIGHT_SPACE,
		  canvasH=TOP_SPACE+chartH+XAXIS_SPACE+BOTTOM_SPACE;


	/******** Chart variables & objects ********/
    this.datasets=[];
	//this.datasetsLabel='Dataset'; // Not necessary
	this.datasetLinks={}; // List of links between dataPoint with same id in different datasets
	this.selectAllDatasets=false;
	this.allowHighlight=pcaData.allowHighlight; // flag to allow or not highlighting
	this.updateHighlight=pcaData.updateHighlight; // object {callback:callback function in case delete/edit,editable:true/false (hl name edition flag)}
	this.highlightedPoints={}; // List of different user-defined labeled set of points
	this.addHighlighting=function(hName,hColor,pointSet,matchPattern) { // for convenience: calls idvis.lib
		idvis.lib.addHighlighting(PCA,hName,hColor,pointSet,matchPattern);
	};
	this.hideUnmatchedPoints=false; // flag to hide/show points not matching active highlighting
	this.axisXtext=pcaData.axisX;
	this.axisYtext=pcaData.axisY;
	this.axisXtitle={}; // must be set to something for idvis to reset with text SVG object
	this.axisYtitle={}; // must be set to something for idvis to reset with text SVG object
	this.forceXto0=0; this.forceYto0=0; // 0 or null: auto, 1: 0, 2: ~0
	this.sameScale=1;
 	this.convertValue=function(axis,value) {return value;};
    this.minValueX=this.maxValueX=this.minValueY=this.maxValueY=null;
	this.chartSettings={reference:{},current:{}};
	this.axisClosure=true; // draw 4 axes
	this.chartMarks=[];
	this.zoomable=true;
	this.zoomText=null;
	this.plotArea=null;
	this.dragArea=null;
	this.dragContext='area';
	this.panProcess={x:0,y:0,dx:0,dy:0};
	this.pointOpacity=pcaData.pointOpacity || 0.7;
	this.pointOnclick=pcaData.pointOnclick;
	this.pointOnList=pcaData.pointOnList;
	this.selectable=true;
	this.searchable=pcaData.searchable || true;

    /********************* Data import *********************/
    /***** Threshold lines *****/
	this.features=[];
	this.features.push(new idvis.feature(PCA,{type:'line',axis:'X',label:'',value:0,color:'#555555',pattern:'- '}));
	this.features.push(new idvis.feature(PCA,{type:'line',axis:'Y',label:'',value:0,color:'#555555',pattern:'- '}));
	this.editableThresholds=[]; // no editable thresholds

    /***** Datasets (only 1)*****/
	idvis.lib.addDataset(PCA,0,{name:'PCA data',color:'#000000',pointSizeRule:{min:15,max:200,scale:1},pointLabels:{size:'Z2'}});

    /***** Data points *****/
	var minZ,maxZ;
    this.addDataAsString=function(dataStrg) {
		var dataset=dataStrg.split(';');
		for (let i=0; i < dataset.length; i++) {
			let data=dataset[i].split(',');
			let dp=new dataPoint(this.datasets[0],data);
			this.datasets[0].data.push(dp);
			if (this.minValueX===null) { // 1st point
				this.minValueX=this.maxValueX=dp.x;
				this.minValueY=this.maxValueY=dp.y;
				minZ=maxZ=dp.z;
			}
			else {
				this.minValueX=Math.min(this.minValueX,dp.x);
				this.maxValueX=Math.max(this.maxValueX,dp.x);
				this.minValueY=Math.min(this.minValueY,dp.y);
				this.maxValueY=Math.max(this.maxValueY,dp.y);
				minZ=Math.min(minZ,dp.z);
				maxZ=Math.max(maxZ,dp.z);
			}
		}
    };

    /********************* PCA display (Raphael) *********************/
    this.draw=function() {
		/* Compute point size based on z axis */
		let refMinSize=35,minZsize=0,scaleZ=1; // default (no z data)
		if (minZ < maxZ) { // z data
			let pointSizeRule=this.datasets[0].params.pointSizeRule;
			refMinSize=pointSizeRule.min;
			minZsize=minZ;
			//var scaleZ=8/(maxZ-minZ); // max r=8 => max size ~200; min/max size overwritten by pointSizeRule
			scaleZ=(pointSizeRule.max-pointSizeRule.min)/((maxZ-minZ)*(maxZ-minZ)*3.14); // radius -> surface
		}
		for (let j=0; j < this.datasets[0].data.length; j++) {
			//this.datasets[0].data[j].computeSize(minZ,scaleZ);
			this.datasets[0].data[j].computeSize(refMinSize,minZsize,scaleZ);
		}
		minZ=1e200; maxZ=-1e-200; // reset in case resetData() & redraw()
		/* Cross-dataset links (links between dataPoint with same externalID) */
		if (this.datasets.length > 1) {
			for (let i=0; i < this.datasets.length; i++) {
				for (let j=0; j < this.datasets[i].data.length; j++) {
					let dataPointID=this.datasets[i].data[j].externalID;
					if (!this.datasetLinks[dataPointID]) {this.datasetLinks[dataPointID]=[];}
					this.datasetLinks[dataPointID].push(this.datasets[i].data[j]);
				}
			}
		}
		/* DIVs */
		canvasDivID=this.divID+'_canvas';
		let formDivID=this.divID+'_form';
		this.div.innerHTML='<TABLE><TR><TD style="min-width:'+canvasW+'px"><DIV id="'+canvasDivID+'"></DIV></TD><TD valign="top" style="min-width:200px"><DIV id="'+formDivID+'"></DIV></TD></TR></TABLE>';
		idvis.lib.drawChartLayout(PCA,pcaData.formPosition,canvasW,canvasDivID,formDivID);
		/* Form menu */
		this.formDiv=document.getElementById(formDivID);
		idvis.lib.initializeForm(PCA);
		/* Canvas */
		this.canvas=Raphael(canvasDivID,canvasW,canvasH);
		/* Back panel */
		this.canvas.rect(0,0,canvasW,canvasH,0).attr({fill:'#fff',stroke:'#000'});
		/* Chart area */
		this.plotArea=this.canvas.rect(chartX,chartY,chartW,chartH,0)
		.attr({fill:'#F3F3F3',stroke:'none'}) //,cursor:"crosshair"
		.drag(function(dx,dy,x,y,e) {idvis.lib.extendDragging(PCA,dx,dy,x,y,e);},
			  function(x,y,e) {idvis.lib.setDragging(PCA,x,y,e);},
			  function() {idvis.lib.endDragging(PCA);}
			  )
		.dblclick(function(){idvis.lib.zoomOut(PCA);});

		/* Display chart */
		idvis.lib.initializeChart(PCA);
    };

	/***** Reset data values (x,y values) *****/
	var tmpDataList={};
	this.resetData=function(dataStrg) {
		let dataset=dataStrg.split(';');
		for (let i=0; i < dataset.length; i++) {
			let data=dataset[i].split(','); // [externalID,newX,newY]
			tmpDataList[data[0]]=[data[1],data[2],data[3]];
			if (data[3]) {
				minZ=Math.min(minZ,data[3]);
				maxZ=Math.max(maxZ,data[3]);
			}
		}
	};

	/***** Redraw graph with new values *****/
	this.redraw=function(titleX,titleY) {
		//let chartX=this.plotArea.attr('x'), chartX0=chartX-1,
		//	chartY=this.plotArea.attr('y'),
		//	chartW=this.plotArea.attr('width'),
		//	chartH=this.plotArea.attr('height');

		this.axisXtext=titleX;
		this.axisYtext=titleY;
		this.axisXtitle.attr({text:this.axisXtext});
		this.axisYtitle.attr({text:this.axisYtext});

		/*** Update x,y,z,size values ***/
		let pointSizeRule=PCA.datasets[0].params.pointSizeRule,
			scaleZ;
		if (minZ < maxZ) { // z data
			//scaleZ=8/(maxZ-minZ); // max r=8 => max size ~200; min size overwritten by pointSizeRule to 2
			scaleZ=(pointSizeRule.max-pointSizeRule.min)/((maxZ-minZ)*(maxZ-minZ)*3.14); // radius -> surface
		}
		let numMatches=0;
		for (let i=0; i < this.datasets[0].data.length; i++) {
			let dp=this.datasets[0].data[i];
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
				dp.computeSize(pointSizeRule.min,minZ,scaleZ);
			}
			else {
				dp.computeSize(35,0,0); // default size when no z data
			}
			//var usedSize=Math.max(Math.min(dp.size,pointSizeRule.max),pointSizeRule.min);
			//dp.point.attr({r:pointSizeRule.scale*Math.sqrt(usedSize/3.14)});
			dp.point.attr({r:pointSizeRule.scale*Math.sqrt(dp.size/3.14)});
		}
		tmpDataList={}; // clear list
		minZ=1e200; maxZ=-1e-200; // reset in case new resetData() & redraw()

		/*** Compute default range ***/
		idvis.lib.computeDefaultRange(PCA);

		/*** Clear & plot chart ***/
		idvis.lib.clearChart(PCA);
		idvis.lib.plotChart(PCA);
	};


/*================================ Nested objects =============================*/

	/****************** Data point object ******************************/
	var dataPoint = function(set,data) {
		this.dataset=set;
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
			let infoStrg=this.label;
			if (type=='min') {
				return this.label;
			}
			else {
				if (this.dataset.params.chart.datasets.length > 1) infoStrg+='\nSet='+this.dataset.params.name;
				infoStrg+='\n'+this.dataset.params.pointLabelX+'='+(this.x*1).toFixed(3)+'\n'+this.dataset.params.pointLabelY+'='+(this.y*1).toFixed(3);
				if (this.z) {
					infoStrg+='\n'+this.dataset.params.pointSizeName+'='+(this.z*1).toFixed(3);
				}
			}
			return infoStrg;
		},
		getX: function() { // needed by idvis.js
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


}; // end of PCA

/*
####>Revision history<####
# 2.0.3 [FEATURE] Accepts pcaData.formPosition attribute (PP 29/04/21)
# 2.0.2 [CHANGE] sizeRule to pointSizeRule (PP 18/06:20)
# 2.0.1 [FEATURE] Optional external pre-search function (PP 09/12/19)
# 2.0.0 Updated pcaPlot.js to code-encapsulated idvis.pcaPlot.js (PP 31/05/17)
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