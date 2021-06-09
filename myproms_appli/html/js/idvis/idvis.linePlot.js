/*
################################################################################
# idvis.linePlot.js     2.0.1                                                  #
# Authors: P. Poullet                                                          #
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
                  Line Plot object
******************************************************************/
idvis.linePlot = function(lineData) {
    var LP=this;
    this.chartType='linePlot';
	this.chartID=idvis.lib.registerChart(this);
    this.divID=lineData.div;
    this.div=document.getElementById(this.divID);


	/******** Chart geometry ********/
    const TOP_SPACE=20,
		  BOTTOM_SPACE=20,
		  LEFT_SPACE=0,
		  RIGHT_SPACE=20,
		  XAXIS_SPACE=30,
		  YAXIS_SPACE=60;
    var chartX=LEFT_SPACE+YAXIS_SPACE+1,
		chartY=TOP_SPACE+1,
		chartW=lineData.width || 400,
		chartH=lineData.height || 400;
	const canvasW=LEFT_SPACE+YAXIS_SPACE+chartW+RIGHT_SPACE,
		  canvasH=TOP_SPACE+chartH+XAXIS_SPACE+BOTTOM_SPACE;
	//this.formPosition=(lineData.legendPosition)? lineData.legendPosition : 'right';


	/******** Chart variables & objects ********/
    this.datasets=[];
	this.datasetsLabel=cpData.datasetsLabel || cpData.dataSetsLabel; // can be undef
	//this.datasetLinks={}; // List of links between dataPoint with same id in different datasets
	this.selectAllDatasets=false;
	this.allowHighlight=false; // flag to allow or not highlighting
	//this.updateHighlight=lineData.updateHighlight; // object {callback:callback function in case delete/edit,editable:true/false (hl name edition flag)}
	//this.highlightedPoints={}; // List of different user-defined labeled set of points
    //var colorList=['#000000','#4AA02C','#F660AB','#FBB917','#0000FF','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B'];
	//const colorList=['#0000FF','#4AA02C','#000000','#F660AB','#FBB917','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B'];
	this.axisXtext=lineData.axisX;
	this.axisYtext=lineData.axisY;
	this.axisXtitle=null; // SVG object
	this.axisYtitle=null; // SVG object
	this.forceXto0=0; this.forceYto0=0; // 0 or null: auto, 1: 0, 2: ~0
	this.sameScale=0;

	/* Zoom-related variables */
	this.zoomable=true; //false;
	this.zoomText=null;
	this.dragArea=null;
	this.dragContext='area';
	this.panProcess={x:0,y:0,dx:0,dy:0};

	this.showDatasets=true;
	this.searchable=false;
	this.selectable=false;
	this.customPointLabel=lineData.pointLabel || null;
	this.onPointExclusion=lineData.onPointExclusion || null;
	this.connectPoints=true;
    this.minValueX=this.maxValueX=this.minValueY=this.maxValueY=null;
	this.chartSettings={reference:{},current:{}};
	this.axisClosure=true; // draw 4 axes
	this.chartMarks=[];
	this.plotArea=null;

    /***** Datasets *****/
    this.addDataset=function(dsIdx,set) { // For convenience: idvis.lib.addDataset can be called directly
        if (!set.connection) set.connection=set.line || {};
		idvis.lib.addDataset(LP,dsIdx,set);
    };
	this.addDataSet=this.addDataset;

    /***** Data points *****/
    this.addDataAsString=function(dsIdx,dataStrg) {
		let dataset=dataStrg.split(';');
		for (let i=0; i < dataset.length; i++) {
			let data=dataset[i].split(',');
			this.datasets[dsIdx].data.push(new dataPoint(this.datasets[dsIdx],data));
			let lastPoint=this.datasets[dsIdx].data[this.datasets[dsIdx].data.length-1];
			if (this.minValueX===null) { // 1st point
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
		}
    };

    /********************* LP display (Raphael) *********************/
    this.draw=function() {
		/* DIVs */
		let canvasDivID=this.divID+'_canvas',
			formDivID=this.divID+'_form';
		//this.div.innerHTML='<TABLE><TR><TD style="min-width:'+canvasW+'px"><DIV id="'+canvasDivID+'"></DIV></TD><TD valign="top" style="min-width:200px"><DIV id="'+formDivID+'"></DIV></TD></TR></TABLE>';
		idvis.lib.drawChartLayout(LP,lineData.formPosition,canvasW,canvasDivID,formDivID);
		this.formDiv=document.getElementById(formDivID);
		/* Form menu */
		idvis.lib.initializeForm(LP);
		/* Canvas */
		this.canvas=Raphael(canvasDivID,canvasW,canvasH);
		/* Back panel */
		this.backPanel=this.canvas.rect(0,0,canvasW,canvasH,0).attr({fill:'#fff',stroke:'#000'});
		/* Chart area */
		this.plotArea=this.canvas.rect(chartX,chartY,chartW,chartH,0)
		.attr({fill:'#F3F3F3',stroke:'none'}) //; //,cursor:"crosshair"
		.drag(function(dx,dy,x,y,e) {idvis.lib.extendDragging(LP,dx,dy,x,y,e);},
			  function(x,y,e) {idvis.lib.setDragging(LP,x,y,e);},
			  function() {idvis.lib.endDragging(LP);}
			  )
		.dblclick(function(){idvis.lib.zoomOut(LP);});
		/* Display chart */
		idvis.lib.initializeChart(LP);
    };


/*====================================== Nested objects ===================================*/

	/*****************************************************************
                  Data point object
	******************************************************************/
	function dataPoint(set,data) {
		this.dataset=set;
		this.label=data[0];
//    	this.externalID=(data[1])? data[1] : data[0];
		this.x=this.y=null;
		this.valueList={};
		//this.size=(data[4])? data[4] : 3;
		this.point=null; // ref to the svg element
		this.pointList=null; // ref to list svg elements
		this.computeMean('x',data[1]);
		this.computeMean('y',data[2]);
		this.noLine=(data[3])? true : false; // not connected to flanking points
	}
	dataPoint.prototype = {
		size : 10,
		computeMean : function(axis,valuesStrg) {
			let values=valuesStrg.split(':').sort(idvis.sortNumber),
				mean=0;
			for (let i=0; i<values.length; i++) {mean+=(1*values[i]);}
			this[axis]=Math.round(100*mean/values.length)/100;
			if (values.length > 1) {this.valueList[axis]=values;}
		},
		getMinValue : function (axis) {
			return (this.valueList && this.valueList[axis])? this.valueList[axis][0] : this[axis];
		},
		getMaxValue : function (axis) {
			return (this.valueList && this.valueList[axis])? this.valueList[axis][this.valueList[axis].length-1] : this[axis];
		},
		info : function(type) {
			let infoStrg=this.label;
			if (type=='min') {
				return this.label;
			}
			else {
				if (this.dataset.params.chart.datasets.length > 1) infoStrg+='\nSet='+this.dataset.params.name;
				infoStrg+='\nX='+(this.x*1).toFixed(3)+'\nY='+(this.y*1).toFixed(3);
			}
			return infoStrg;
		},
		getX : function() {
			return this.x;
		},
		getY : function() {
			return this.y;
		}
	};

};

/*
####>Revision history<####
# 2.0.1 [UPDATE] Support for datasetsLabel (PP 06/05/21)
# 2.0.0 Updated linePlot.js to code-encapsulated idvis.linePlot.js (PP 13/06/17)
# 1.0.4 Uses chart registering for dynamic html calls (PP 20/02/16)
# 1.0.3 GPL license (PP 23/09/13)
# 1.0.2 Commented LP.datasetLinks declaration (PP 29/01/13)
# 1.0.1 Added noLine attribute to dataPoint for exclusion switch (PP 02/12/12)
# 1.0.0 Production version & form position option (PP 30/05/12)
# 0.0.1 starting script (PP 12/05/12)
*/