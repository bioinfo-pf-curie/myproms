/*
################################################################################
# volcanoPlot2.js      1.1.8                                                   #
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
function vpDataPoint(set,data) {
//console.log(data);
    this.dataSet=set;
    this.label=data[0];
    this.externalID=(data[1])? data[1] : data[0];
	if (data[2]=='-') {
		this.subChart=set.params.chart.subChart.minusInf;
		this.foldChange=data[2];
		this.pepFreq=data[3];
		this.x=Math.random();
	}
	else if (data[2]=='+') {
		this.subChart=set.params.chart.subChart.plusInf;
		this.foldChange=data[2];
		this.pepFreq=data[3];
		this.x=Math.random();
	}
	else {
		this.subChart=set.params.chart.subChart.volcano;
		this.foldChange=(data[2] >= 1)? (data[2]*1).toFixed(2) : '1/'+(1/data[2]).toFixed(2);
		this.pvalue=(data[3] < 1e-300)? 1e-300 : data[3];
		this.log2fc=set.params.chart.subChart.volcano.convertValue('X',data[2]);
		this.mlog10pv=set.params.chart.subChart.volcano.convertValue('Y',this.pvalue);
	}
    this.size=(data[4])? data[4] : 3;
    this.point=null; // ref to the svg element
	this.highlightNames=[];
}
vpDataPoint.prototype = {
    info: function(type) {
		var infoStrg=this.label;
		if (type=='min') {
			return this.label;
		}
		else {
			if (this.dataSet.params.chart.dataSets.length > 1) infoStrg+='\nSet='+this.dataSet.params.name;
			if (this.subChart.name=='volcano') {
				infoStrg+='\nFC='+this.foldChange+'\np-value';
				infoStrg+=(this.pvalue==1.e-300)? '<=1.e-300' : '='+this.pvalue;
			}
			else {
				infoStrg+=(this.foldChange=='-')? '\nFC=1/∞' : '\nFC=∞';// ∞
				infoStrg+='\npep. freq='+this.pepFreq;
			}
			if (this.dataSet.params.sizeName) {infoStrg+='\n'+this.dataSet.params.sizeName+'='+this.size;}
		}
		return infoStrg;
    },
    getX: function() {
		if (this.subChart.name=='volcano') {return this.log2fc;}
		else {return this.x;}
	},
	getY: function() {
		if (this.subChart.name=='volcano') {return this.mlog10pv;}
		else {return this.pepFreq;}
    }
}

/*****************************************************************
                  Volcano Plot object
******************************************************************/
function volcanoPlot(volcano) {
    var VP=this;
	this.chartID=registerInCuBioJS(this);
    this.divID=volcano.div;
    this.div=document.getElementById(this.divID);
	/* Graphic & form DIVs */
	var canvasDivID=this.divID+'_canvas';
	var formDivID=this.divID+'_form';
	this.getDivID=function(){return canvasDivID;};
	this.exportAsImage=volcano.exportAsImage; // null or array [button name, image file name, action script path]

    /******** Chart variables & objects ********/
	this.multiChart=false; // default
    this.dataSets=[];
	this.dataSetLinks={}; // List of links between dataPoint with same id in different dataSets
	this.selectAllDataSets=false;
	this.allowHighlight=volcano.allowHighlight; // flag to allow or not highlighting
	this.updateHighlight=volcano.updateHighlight; // object {callback:callback function in case delete/edit,editable:true/false (hl name edition flag)}
	this.highlightedPoints={}; // List of different user-defined labeled set of points
    var pointDefaultSize=(volcano.pointDefaultSize)? volcano.pointDefaultSize : 3;
	this.pointOpacity=(volcano.pointOpacity)? volcano.pointOpacity : 0.7; // point opacity
	//var colorList=['#000000','#4AA02C','#F660AB','#FBB917','#0000FF','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B'];
	var colorList=['#0000FF','#4AA02C','#000000','#F660AB','#FBB917','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B'];
	this.subChart={volcano:{},minusInf:{},plusInf:{}};
	this.activeCharts=[];
	var charts=this.subChart;
	var vp=charts.volcano;
	var mi=charts.minusInf;
	var pi=charts.plusInf;
	mi.name='minusInf';
	vp.name='volcano';
	pi.name='plusInf';

	/* Axes */
	vp.axisXtext='Log2(fold change B/A)'; //,minusInf:null,plusInf:null
	vp.axisYtext='-Log10(p-value)';
	mi.axisXtext=pi.axisXtext='';
	mi.axisYtext=pi.axisYtext='#Peptides/100 aa';
	vp.forceXto0=mi.forceXto0=pi.forceXto0=0;
	mi.forceYto0=pi.forceYto0=1; // 0 or null: auto, 1: 0, 2: ~0
	vp.forceYto0=2;
	mi.noTicksX='Only in A';
	pi.noTicksX='Only in B';
	vp.convertValue=function(axis,startValue) {
		var convertedValue=(axis=='X')? Math.log(startValue)/0.693147181 : -Math.log(startValue)/2.302585093; //LOG2 LOG10
		return convertedValue;
	};
    //this.minValueX,this.maxValueX,this.minValueY,this.maxValueY;
	for (var c in charts) {
		if (charts.hasOwnProperty(c)) {
			charts[c].mainChart=VP;
			charts[c].chartSettings={reference:{},current:{}};
			charts[c].axisClosure=false; // draw 2 axes only
			charts[c].chartMarks=[];
			charts[c].plotArea=null;
			charts[c].dragArea=null;
			charts[c].dragContext='area';
			if (c=='volcano') {
				charts[c].zoomable=true;
				charts[c].zoomText=null;
				charts[c].panProcess={x:0,y:0,dx:0,dy:0}; //new Object();
			}
		}
	}
	this.pointOnclick=volcano.pointOnclick;
	this.pointOnList=volcano.pointOnList;
	this.selectable=true;
	this.searchable=volcano.searchable || true; // always searchable

	/***** Threshold lines *****/
	this.editThreshold=true;
	this.thresholdLines=[];
	var minThresFC=(volcano.minFoldChange)? volcano.minFoldChange : (volcano.foldChange)? 1/volcano.foldChange : null;
	if (minThresFC) {
		this.thresholdLines.push(new thresholdLine(vp,'X','min. fold change',minThresFC,'#00A000','',1,1,true));
	}
	var maxThresFC=(volcano.maxFoldChange)? volcano.maxFoldChange : (volcano.foldChange)? volcano.foldChange : null;
	if (maxThresFC) {
		this.thresholdLines.push(new thresholdLine(vp,'X','max. fold change',maxThresFC,'#00A000','',1,1,true));
	}
	if (volcano.pValue) {
		this.thresholdLines.push(new thresholdLine(vp,'Y','max. p-value',volcano.pValue,'#FF0000','',0,0,true));
	}
vp.thresholdLines=this.thresholdLines;

    /********************* Data import *********************/
    /***** Datasets *****/
    this.addDataSet=function(dsIdx,set) {
		this.dataSets[dsIdx]={};
		this.dataSets[dsIdx].params={};
		this.dataSets[dsIdx].params.chart=VP;
		this.dataSets[dsIdx].params.index=dsIdx;
		this.dataSets[dsIdx].params.name=(set.name)? set.name : 'dataSet '+dsIdx;
		this.dataSets[dsIdx].params.color=(set.color)? set.color : colorList[dsIdx % colorList.length];
		this.dataSets[dsIdx].params.sizeName=(set.sizeName)? set.sizeName : null;
		this.dataSets[dsIdx].params.sizeRule=(set.sizeRule)? set.sizeRule : null;
		//this.dataSets[dsIdx].params.onclick=(set.onclick)? set.onclick : null;
		this.dataSets[dsIdx].params.visible=true;
		this.dataSets[dsIdx].data=[];
    }
    /***** Data points *****/
    this.addDataAsString=function(dsIdx,dataStrg) {
		var dataSet=dataStrg.split(';');
		for (var i=0; i < dataSet.length; i++) {
			var data=dataSet[i].split(',');
			if (isNaN(data[4]) && pointDefaultSize) {data[4]=pointDefaultSize;}
			var lastPoint=new vpDataPoint(this.dataSets[dsIdx],data);
			this.dataSets[dsIdx].data.push(lastPoint);
			var C=lastPoint.subChart;
			//if (chName != 'volcano') this.multiChart=true;
			if (C !== vp) this.multiChart=true;
			if (C.minValueX==null) { // 1st point
				C.minValueX=lastPoint.getX();
				C.maxValueX=lastPoint.getX();
				C.minValueY=lastPoint.getY();
				C.maxValueY=lastPoint.getY();
			}
			else {
				C.minValueX=Math.min(C.minValueX,lastPoint.getX());
				C.maxValueX=Math.max(C.maxValueX,lastPoint.getX());
				C.minValueY=Math.min(C.minValueY,lastPoint.getY());
				C.maxValueY=Math.max(C.maxValueY,lastPoint.getY());
			}
		}
		/* Overwrites limits for +/-infinity charts */
		if (this.multiChart) {
			mi.minValueX=pi.minValueX=0;
			mi.maxValueX=pi.maxValueX=1;
			mi.minValueY=pi.minValueY=0;
			if (!mi.maxValueY) mi.maxValueY=0;
			if (!pi.maxValueY) pi.maxValueY=0;
			if (mi.maxValueY < pi.maxValueY) {mi.maxValueY=pi.maxValueY;}
			else {pi.maxValueY=mi.maxValueY;}
		}
    }

    /******** Chart geometry ********/
    var TOP_SPACE=20;
    var BOTTOM_SPACE=20;
    var LEFT_SPACE=10;
	var MIDDLE_SPACE=20;
    var RIGHT_SPACE=20;
    var XAXIS_SPACE=30;
    var YAXIS_SPACE=50;
	//var chartX={},chartY={},chartW={},chartH={};
//var chartX=LEFT_SPACE+YAXIS_SPACE+1;
//var chartY=TOP_SPACE+1;
    vp.chartW=(volcano.width)? volcano.width : 400;
    vp.chartH=(volcano.height)? volcano.height : 400;
	vp.chartY=TOP_SPACE+1;
    var canvasH=TOP_SPACE+vp.chartH+XAXIS_SPACE+BOTTOM_SPACE;
	var canvasW;


	/********************* Volcano display (Raphael) *********************/
    this.draw=function() {

		/* Chart structure */
		if (this.multiChart) {
			this.activeCharts=['minusInf','volcano','plusInf'];
			mi.chartW=pi.chartW=50;
			mi.chartH=pi.chartH=vp.chartH;
			mi.chartY=pi.chartY=vp.chartY;
			mi.chartX=LEFT_SPACE+YAXIS_SPACE+1;
			vp.chartX=mi.chartX+MIDDLE_SPACE+mi.chartW+YAXIS_SPACE;
			pi.chartX=vp.chartX+MIDDLE_SPACE+vp.chartW+YAXIS_SPACE;
			canvasW=pi.chartX+pi.chartW+RIGHT_SPACE-1;
		}
		else {
			this.activeCharts=['volcano'];
			vp.chartX=LEFT_SPACE+YAXIS_SPACE+1;
			//canvasW=LEFT_SPACE+YAXIS_SPACE+chartW+RIGHT_SPACE;
			canvasW=vp.chartX+vp.chartW+RIGHT_SPACE-1;
		}

		/* Reset min/max ranges to make thresholds visible */
		if (minThresFC && vp.minValueX > vp.convertValue('X',minThresFC)) {vp.minValueX=vp.convertValue('X',minThresFC);}
		if (maxThresFC && vp.maxValueX < vp.convertValue('X',maxThresFC)) {vp.maxValueX=vp.convertValue('X',maxThresFC);}
		if (volcano.pValue && vp.maxValueY < vp.convertValue('Y',volcano.pValue)) {vp.maxValueY=vp.convertValue('Y',volcano.pValue);}
//console.log('minThX='+this.convertValue('X',minThresFC)+', maxThX='+this.convertValue('X',maxThresFC));
//console.log('minX='+this.minValueX+', maxX='+this.maxValueX);

		/* DIVs */
		//this.div.innerHTML='<TABLE><TR><TD><DIV id="'+canvasDivID+'"></DIV></TD><TD valign="top"><DIV id="'+formDivID+'"></DIV></TD></TR></TABLE>';
		this.div.innerHTML='<TABLE><TR><TD style="min-width:'+canvasW+'px"><DIV id="'+canvasDivID+'"></DIV></TD><TD valign="top" style="min-width:200px"><DIV id="'+formDivID+'"></DIV></TD></TR></TABLE>';
		this.formDiv=document.getElementById(formDivID);
		/* Form menu */
		initializeForm(VP);
		/* Canvas */
		this.canvas=Raphael(canvasDivID,canvasW,canvasH);
		/* Back panel */
		this.backPanel=this.canvas.rect(0,0,canvasW,canvasH,0).attr({fill:'#fff',stroke:'#000'});
		/* Chart area(s) */
		for (var i=0; i<this.activeCharts.length; i++) {
			var C=VP.subChart[this.activeCharts[i]];
			C.plotArea=this.canvas.rect(C.chartX,C.chartY,C.chartW,C.chartH,0)
			.attr({fill:'#F3F3F3',stroke:'none'}) //,cursor:"crosshair"
			.data({chart:C})
			.drag(function(dx,dy,x,y,e){extendDragging(this.data('chart'),dx,dy,x,y,e);}, //_extendDragging,
				  function(x,y,e){setDragging(this.data('chart'),x,y,e);}, //_setDragging,
				  function(){endDragging(this.data('chart'));} //_endDragging
				  );
		}
		vp.plotArea.dblclick(function(){zoomOut(this.data('chart'))});


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
		initializeChart(VP);

    }

} // end of VP

/*
####>Revision history<####
# 1.1.8 [FEATURE] Optional external pre-search function (PP 04/10/19)
# 1.1.7 "+/-∞" replaced by "Only in A/B" (PP 31/05/19)
# 1.1.6 Uses chart registering for dynamic html calls (PP 20/02/16)
# 1.1.5 Renamed dataPoint object to vpDataPoint for compatiblity with genericPlot library (PP 28/10/15)
# 1.1.4 Handle dataSet.params.sizeRule (PP 21/10/15)
# 1.1.3 Parameters for exporting SVG to PNG (PP 11/02/14)
# 1.1.2 GPL license (PP 23/09/13)
# 1.1.1 Added point opacity (PP 15/07/13)
# 1.1.0 Added sub graphes for +/- infinite ratios (PP 17/03/13)
# 1.0.7 Uncomment updateHighlight declaration (PP 29/01/13)
# 1.0.6 Handles number of datasets > colorList (PP 10/01/13)
# 1.0.5 Added definable pointDefaultSize, 3 otherwise (PP 21/09/12)
# 1.0.4 Min p-value set to 1e-300 to prevent JS from freezing (PP 30/07/12)
# 1.0.3 Handles properly dataset with no data (PP 26/07/12)
# 1.0.2 Skips size info in popup if no dataset sizeName (PP 06/07/12)
# 1.0.1 Compatibility with chartLibrary 1.0.2 (PP 13/05/12)
# 1.0.0 Production version (PP 01/03/12)
*/
