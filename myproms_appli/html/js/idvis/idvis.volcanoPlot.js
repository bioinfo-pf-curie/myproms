/*
################################################################################
# idvis.volcanoPlot.js      2.1.0                                              #
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
                  Volcano Plot object
******************************************************************/
idvis.volcanoPlot = function(vpData) {
   var VP=this;
   this.chartType=vpData.chartType || 'volcanoPlot'; // volcanoPlot or maPlot
   this.divID=vpData.div;
   this.div=document.getElementById(this.divID);
   this.chartID=idvis.lib.registerChart(this);
	/* Graphic & form DIVs */
	const canvasDivID=this.divID+'_canvas',
		  formDivID=this.divID+'_form';
	this.getDivID=function(){return canvasDivID;};
	this.exportAsImage=vpData.exportAsImage; // null or array [button name, image file name, action script path]

	/******** Chart variables & objects ********/
	this.multiChart=false; // default
	this.datasets=[];
	this.datasetsLabel=vpData.datasetsLabel || vpData.dataSetsLabel; // can be undef
	this.datasetLinks={}; // List of links between dataPoint with same id in different datasets
	this.selectAllDatasets=false;
	this.allowHighlight=vpData.allowHighlight; // flag to allow or not highlighting
	this.updateHighlight=vpData.updateHighlight; // object {callback:callback function in case delete/edit,editable:true/false (hl name edition flag)}
	this.multiPointDatasetRule=vpData.multiDatasetRule || vpData.multiPointDatasetRule; // string, updated by idvis if undefined
	this.highlightedPoints={}; // List of different user-defined labeled set of points
	this.addHighlighting=function(hName,hColor,pointSet,matchPattern) { // for convenience: calls idvis.lib
		idvis.lib.addHighlighting(VP,hName,hColor,pointSet,matchPattern);
	};
	if (this.chartType==='volcanoPlot') {
		this.pValueType=vpData.pValueType || 'p-value';
	}
	this.multiDatasetRule=vpData.multiDatasetRule; // can be undefined (updated by addDataset)
	var pointDefaultSize=vpData.pointDefaultSize || 3;
	this.pointOpacity=(vpData.pointOpacity)? vpData.pointOpacity : 0.7; // point opacity
	//var colorList=['#000000','#4AA02C','#F660AB','#FBB917','#0000FF','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B'];
	//const colorList=['#0000FF','#4AA02C','#000000','#F660AB','#FBB917','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B'];
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
    vp.axisYtext=(VP.chartType==='maPlot')? 'Mean Log2(Intensities A,B)' : '-Log10('+this.pValueType+')';
	mi.axisXtext=pi.axisXtext='';
    vp.forceXto0=mi.forceXto0=pi.forceXto0=0;
    if (VP.chartType==='volcanoPlot') {
        mi.axisYtext=pi.axisYtext='#Peptides/100 aa';
        mi.forceYto0=pi.forceYto0=1; // 0 or null: auto, 1: 0, 2: ~0
        vp.forceYto0=2;
    }
    else {
        mi.axisYtext='Log2(Intensity A)';
        pi.axisYtext='Log2(Intensity B)';
        mi.forceYto0=pi.forceYto0=vp.forceYto0=2;
    }
	
	mi.noTicksX='Ony in A'; //'-INF'; 
	pi.noTicksX='Only in B';  //'+INF';
    var LOG2=Math.log(2),
        LOG10=Math.log(10);
	vp.convertValue=function(axis,startValue,startValue2) {
		if (axis==='X') {return Math.log(startValue)/LOG2;}
        else { // Y
            if (VP.chartType==='volcanoPlot') {return -Math.log(startValue)/LOG10;}
            else { // maPlot
                if (startValue && startValue2) {return (Math.log(startValue)+Math.log(startValue2))/(2*LOG2);} // mean (log2A,log2B)
                else if (startValue) {return Math.log(startValue)/LOG2;}
                else {return Math.log(startValue2)/LOG2;}
            }
        }
	};
    if (this.chartType==='maPlot') {
        mi.convertValue=vp.convertValue;
        pi.convertValue=vp.convertValue;
    }
    //this.minValueX,this.maxValueX,this.minValueY,this.maxValueY;
	for (let c in charts) {
		if (charts.hasOwnProperty(c)) {
			charts[c].mainChart=VP;
			charts[c].chartSettings={reference:{},current:{}};
			charts[c].axisClosure=false; // draw 2 axes only
			charts[c].chartMarks=[];
			charts[c].plotArea=null;
			charts[c].dragArea=null;
			charts[c].dragContext='area';
			if (c==='volcano') {
				charts[c].zoomable=true;
				charts[c].zoomText=null;
				charts[c].panProcess={x:0,y:0,dx:0,dy:0};
			}
		}
	}
	this.pointOnClick=vpData.pointOnClick;
	this.pointOnList=vpData.pointOnList;
	this.selectable=true;
	this.searchable=vpData.searchable || true; // always searchable

	/***** Threshold lines *****/
	this.features=[];
	var minThresFC=(vpData.minFoldChange)? vpData.minFoldChange : (vpData.foldChange)? 1/vpData.foldChange : 0.5;
    this.features.push(new idvis.feature(vp,{type:'line',axis:'X',label:'min. fold change',value:minThresFC,color:'#00A000',keepAbove1:true,roundValue:true,editable:true}));
	var maxThresFC=(vpData.maxFoldChange)? vpData.maxFoldChange : (vpData.foldChange)? vpData.foldChange : 2;
	this.features.push(new idvis.feature(vp,{type:'line',axis:'X',label:'max. fold change',value:maxThresFC,color:'#00A000',keepAbove1:true,roundValue:true,editable:true}));
    var pValueThres;
   if (this.chartType==='volcanoPlot') {
        pValueThres=vpData.pValue || 0.05;
        this.features.push(new idvis.feature(vp,{type:'line',axis:'Y',label:'max. '+this.pValueType,value:pValueThres,color:'#FF0000',editable:true}));
    }
	vp.features=this.features;
    this.editableThresholds=this.features; // all thresholds are editable
    this.dimArea={label:'Dim non-significant proteins',active:true}; // Controls area where datapoints are dimmed
   
    /* Update points transparency */
    function _updatePointOpacity() {
        let minLog2fc=VP.features[0].value,
            maxLog2fc=VP.features[1].value,
            mlog10pv=(VP.chartType==='volcanoPlot')? VP.features[2].value : null;
        for (let dsIdx=0; dsIdx<VP.datasets.length; dsIdx++) {
            for (let i=0; i<VP.datasets[dsIdx].data.length; i++) {
				let dp=VP.datasets[dsIdx].data[i];
                if (dp.subChart.name !== 'volcano') continue;
                let pOpacity;
                if (VP.dimArea.active===true) {
                    if (VP.chartType==='volcanoPlot') {pOpacity=(dp.mlog10pv >= mlog10pv && (dp.log2fc <= minLog2fc || dp.log2fc >= maxLog2fc))? VP.pointOpacity : idvis.dimPointOpacity;}
                    else {pOpacity=(dp.log2fc <= minLog2fc || dp.log2fc >= maxLog2fc)? VP.pointOpacity : idvis.dimPointOpacity;} // maPlot
                }
                else {pOpacity=VP.pointOpacity;}
                //dp.point.attr({'fill-opacity':pOpacity});
                dp.point.animate({'fill-opacity':pOpacity},400,"easeInOut");
            }
        }
    }
    /* Callback after idvis.updateThreshold() => update points transparency */
    this.onThresholdUpdate=_updatePointOpacity;

    /********************* Data import *********************/
    /***** Datasets *****/
    this.addDataset=function(dsIdx,set) { // For convenience: idvis.lib.addDataset can be called directly
		idvis.lib.addDataset(VP,dsIdx,set);
    };
	this.addDataSet=this.addDataset; // just to be safe

    /***** Data points *****/
    //mi.minValueY=pi.minValueY=vp.minValueY=
    //mi.maxValueY=pi.maxValueY=vp.maxValueY=null; // null coherced to 0 if not changed (needed)
    this.addDataAsString=function(dsIdx,dataStrg) {
        var dataPointObject=(VP.chartType==='maPlot')? maDataPoint : vpDataPoint;
		let dataset=dataStrg.split(';');
		for (let i=0; i < dataset.length; i++) {
			let data=dataset[i].split(',');
			if (!idvis.isNumber(data[4]) && pointDefaultSize) {data[4]=pointDefaultSize;}
			//let lastPoint=new vpDataPoint(this.datasets[dsIdx],data);
			let lastPoint=new dataPointObject(this.datasets[dsIdx],data);
            this.datasets[dsIdx].data.push(lastPoint);
			let C=lastPoint.subChart;
			//if (chName != 'volcano') this.multiChart=true;
			if (C !== vp) this.multiChart=true;
			if (C.minValueX===undefined) { // 1st point
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
    };

    /******** Chart geometry ********/
    const TOP_SPACE=20;
    const BOTTOM_SPACE=20;
    const LEFT_SPACE=10;
	const MIDDLE_SPACE=20;
    const RIGHT_SPACE=20;
    const XAXIS_SPACE=30;
    const YAXIS_SPACE=50;
	//var chartX={},chartY={},chartW={},chartH={};
//var chartX=LEFT_SPACE+YAXIS_SPACE+1;
//var chartY=TOP_SPACE+1;
    vp.chartW=vpData.width || 400;
    vp.chartH=vpData.height || 400;
	vp.chartY=TOP_SPACE+1;
    var canvasH=TOP_SPACE+vp.chartH+XAXIS_SPACE+BOTTOM_SPACE;
	var canvasW;


	/********************* Volcano display (Raphael) *********************/
    this.draw = function() {

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
            
            /* Overwrites limits for +/-infinity charts */
			mi.minValueX=pi.minValueX=0;
			mi.maxValueX=pi.maxValueX=1;
            if  (this.chartType==='volcanoPlot') {
                mi.minValueY=pi.minValueY=0;
                if (!mi.maxValueY) mi.maxValueY=0;
                if (!pi.maxValueY) pi.maxValueY=0;
                if (mi.maxValueY < pi.maxValueY) {mi.maxValueY=pi.maxValueY;}
                else {pi.maxValueY=mi.maxValueY;}
            }
            else { // same Y range for maPlot
                let goodMin=[], goodMax=[];
                for (let C of [mi,vp,pi]) {
                    if (idvis.isNumber(C.minValueY)) goodMin.push(C.minValueY*1);
                    if (idvis.isNumber(C.maxValueY)) goodMax.push(C.maxValueY*1);
                }
                mi.minValueY=pi.minValueY=vp.minValueY=Math.min(...goodMin);
                mi.maxValueY=pi.maxValueY=vp.maxValueY=Math.max(...goodMax);
            }
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
		if (pValueThres && vp.maxValueY < vp.convertValue('Y',pValueThres)) {vp.maxValueY=vp.convertValue('Y',pValueThres);}
//console.log('minThX='+this.convertValue('X',minThresFC)+', maxThX='+this.convertValue('X',maxThresFC));
//console.log('minX='+this.minValueX+', maxX='+this.maxValueX);

		/* DIVs */
		//this.div.innerHTML='<TABLE><TR><TD style="min-width:'+canvasW+'px"><DIV id="'+canvasDivID+'"></DIV></TD><TD valign="top" style="min-width:200px"><DIV id="'+formDivID+'"></DIV></TD></TR></TABLE>';
		idvis.lib.drawChartLayout(VP,vpData.formPosition,canvasW,canvasDivID,formDivID);
		
		/* Form menu */
		this.formDiv=document.getElementById(formDivID);
		idvis.lib.initializeForm(VP);
		/* Canvas */
		this.canvas=Raphael(canvasDivID,canvasW,canvasH);
		/* Back panel */
		this.backPanel=this.canvas.rect(0,0,canvasW,canvasH,0).attr({fill:'#fff',stroke:'#000'});
		/* Chart area(s) */
        var _handleChartDrag=function(type) {
            return (type==='ext')? function(dx,dy,x,y,e){idvis.lib.extendDragging(this.data('chart'),dx,dy,x,y,e);} //_extendDragging,
                : (type==='set')? function(x,y,e){idvis.lib.setDragging(this.data('chart'),x,y,e);} //_setDragging,
                : function(){idvis.lib.endDragging(this.data('chart'));}; //_endDragging;
        };
		for (let i=0; i<this.activeCharts.length; i++) {
			let C=VP.subChart[this.activeCharts[i]];
			C.plotArea=this.canvas.rect(C.chartX,C.chartY,C.chartW,C.chartH,0)
			.attr({fill:'#F3F3F3',stroke:'none'}) //,cursor:"crosshair"
			.data({chart:C})
            /*
			.drag(function(dx,dy,x,y,e){idvis.lib.extendDragging(this.data('chart'),dx,dy,x,y,e);}, //_extendDragging,
				  function(x,y,e){idvis.lib.setDragging(this.data('chart'),x,y,e);}, //_setDragging,
				  function(){idvis.lib.endDragging(this.data('chart'));} //_endDragging
				  );
			*/
            .drag(_handleChartDrag('ext'),_handleChartDrag('set'),_handleChartDrag('end'));
		}
		vp.plotArea.dblclick(function(){idvis.lib.zoomOut(this.data('chart'));});


		/***** Cross-dataset links (links between dataPoint with same externalID) *****/
		if (this.datasets.length > 1) {
			for (let i=0; i < this.datasets.length; i++) {
				for (let j=0; j < this.datasets[i].data.length; j++) {
					let dataPointID=this.datasets[i].data[j].externalID;
					if (!this.datasetLinks[dataPointID]) {this.datasetLinks[dataPointID]=[];}
					this.datasetLinks[dataPointID].push(this.datasets[i].data[j]);
				}
			}
		}

		/***** Display chart(s) *****/
		idvis.lib.initializeChart(VP);
        
        /***** Change point opacity according to FC/PVAL thresholds *****/
        _updatePointOpacity();

    };
    

/*====================================== Nested objects ===================================*/

	/*****************************************************************
						Data point object
	******************************************************************/
     /**** Volcano Plot point ****/
	var vpDataPoint = function(set,data) {
		this.dataset=set;
		this.label=data[0];
		this.externalID=(data[1])? data[1] : data[0];
		if (data[2]==='-') {
			this.subChart=set.params.chart.subChart.minusInf;
			this.foldChange=data[2];
			this.pepFreq=data[3];
			this.x=Math.random();
		}
		else if (data[2]==='+') {
			this.subChart=set.params.chart.subChart.plusInf;
			this.foldChange=data[2];
			this.pepFreq=data[3];
			this.x=Math.random();
		}
		else {
			this.subChart=set.params.chart.subChart.volcano;
			this.foldChange=(data[2] >= 1)? (data[2]*1).toFixed(2) : '1/'+(1/data[2]).toFixed(2);
			this.log2fc=set.params.chart.subChart.volcano.convertValue('X',data[2]*1);
            this.pvalue=(data[3] < 1e-300)? 1e-300 : idvis.formatNumber(data[3]);
            this.mlog10pv=set.params.chart.subChart.volcano.convertValue('Y',this.pvalue);
		}
		this.size=(data[4])? data[4] : pointDefaultSize;
		this.point=null; // ref to the svg element
		this.highlightNames=[];
	};
	vpDataPoint.prototype = {
		info: function(type) {
			let infoStrg=this.label;
			if (type==='min') {
				return this.label;
			}
			else {
				if (this.dataset.params.chart.datasets.length > 1) infoStrg+='\nSet='+this.dataset.params.name;
				if (this.subChart.name==='volcano') {
                    infoStrg+='\nFC='+this.foldChange+'\n'+this.dataset.params.chart.pValueType;
					infoStrg+=(this.pvalue==1.e-300)? '<=1.e-300' : '='+this.pvalue;
				}
				else {
					infoStrg+=(this.foldChange==='-')? '\nFC=1/INF' : '\nFC=INF';
					infoStrg+='\npep. freq='+this.pepFreq;
				}
				if (this.dataset.params.pointSizeName) {infoStrg+='\n'+this.dataset.params.pointSizeName+'='+this.size;}
			}
			return infoStrg;
		},
		getX: function() {
			if (this.subChart.name==='volcano') {return this.log2fc;}
			else {return this.x;}
		},
		getY: function() {
			if (this.subChart.name==='volcano') {return this.mlog10pv;}
			else {return this.pepFreq;}
		}
	};
    
    /**** MA Plot point ****/
    var maDataPoint = function(set,data) {
//console.log(data);
		this.dataset=set;
		this.label=data[0];
		this.externalID=(data[1])? data[1] : data[0];
		if (data[2]==='-') {
			this.subChart=set.params.chart.subChart.minusInf;
			this.foldChange=data[2];
			this.x=Math.random();
		}
		else if (data[2]==='+') {
			this.subChart=set.params.chart.subChart.plusInf;
			this.foldChange=data[2];
			this.x=Math.random();
		}
		else {
			this.subChart=set.params.chart.subChart.volcano;
			this.foldChange=(data[2] >= 1)? (data[2]*1).toFixed(2) : '1/'+(1/data[2]).toFixed(2);
			this.log2fc=set.params.chart.subChart.volcano.convertValue('X',data[2]*1);
        }
        let intensities=data[3].split(':');
        this.intensityA=intensities[0]*1;
        this.intensityB=intensities[1]*1;
        this.meanLog=set.params.chart.subChart.volcano.convertValue('Y',this.intensityA,this.intensityB);
		this.size=(data[4])? data[4] : pointDefaultSize;
		this.point=null; // ref to the svg element
		this.highlightNames=[];
	};
	maDataPoint.prototype = {
		info: function(type) {
			let infoStrg=this.label;
			if (type==='min') {
				return this.label;
			}
			else {
				if (this.dataset.params.chart.datasets.length > 1) infoStrg+='\nSet='+this.dataset.params.name;
                infoStrg+=(this.foldChange==='-')? '\nFC=1/INF' : (this.foldChange==='+')? '\nFC=INF': '\nFC='+this.foldChange;
                infoStrg+='\nIntensity in A='+this.intensityA+'\nIntensity in B='+this.intensityB;
				if (this.dataset.params.pointSizeName) {infoStrg+='\n'+this.dataset.params.pointSizeName+'='+this.size;}
			}
			return infoStrg;
		},
		getX: function() {
			if (this.subChart.name==='volcano') {return this.log2fc;}
			else {return this.x;}
		},
		getY: function() {
			return this.meanLog;
		}
	};

}; // end of VP


/*
####>Revision history<####
# 2.1.0 [FEATURE] Support for MA-plot and form positioning & other minor changes (PP 05/05/21)
# 2.0.6 [BUGFIX] Minor bug fix in points opacity setting (PP 14/05/20)
# 2.0.5 [FEATURE] Handles p-value type (eg. adjusted) (PP 28/02/20)
# 2.0.4 [FEATURE] Accept an external search-text conversion function (PP 29/09/19)
# 2.0.3 "Not in S2/S1" replaced by "Only in A/B" (PP 03/06/19)
# 2.0.2 "+/-INF" replaced by "Not in S2/S1" (PP 30/05/19)
# 2.0.1 Non-significant points (according to thresholds) have higher transparency (PP 25/09/17) 
# 2.0.0 Updated volcanoPlot2.js to code-encapsulated idvis.volcanoPlot.js (PP 30/05/17)
# 1.1.6 Uses chart registering for dynamic html calls (PP 20/02/16)
# 1.1.5 Renamed dataPoint object to vpDataPoint for compatiblity with genericPlot library (PP 28/10/15)
# 1.1.4 Handle dataset.params.pointSizeRule (PP 21/10/15)
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