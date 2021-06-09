/*
################################################################################
# idvis.genericPlot.js    2.0.3                                                #
# requires idvis.js                                                            #
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
                  Generic Plot object
******************************************************************/
idvis.genericPlot = function(plotData) {
    var GP=this;
    this.chartType='genericPlot';
    this.divID=plotData.div;
    this.div=document.getElementById(this.divID);
    this.chartID=idvis.lib.registerChart(this);
	/* Graphic & form DIVs */
	const canvasDivID=this.divID+'_canvas',
		  formDivID=this.divID+'_form';
	this.getDivID=function(){return canvasDivID;};

	//var colorList=['#0000FF','#4AA02C','#000000','#F660AB','#FBB917','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B'];

    /******** Chart variables & objects ********/
    this.datasets=[];
	this.datasetsLabel=plotData.datasetsLabel || plotData.dataSetsLabel; // can be undef
	this.datasetLinks={}; // List of links between dataPoint with same id in different datasets
	this.selectAllDatasets=false;
	//this.thresholdLines=[];
	this.features=[];
	//this.editThreshold=plotData.editThreshold; // if true: all th are editable
	this.editableThresholds=[]; // list of editable features
    if (plotData.dimAreaData) { // {minX:<featureIdx>,maxX:<featureIdx>, minY:<featureIdx>, maxY:<featureIdx>,label:<chechbox label>}
        this.dimArea=plotData.dimAreaData;
        this.dimArea.rule=plotData.dimAreaData.rule || 'in';
        this.dimArea.active=true;
    }
	this.exportAsImage=plotData.exportAsImage;
	this.allowHighlight=plotData.allowHighlight; // flag to allow or not highlighting
	this.updateHighlight=plotData.updateHighlight; // object {callback:callback function in case delete/edit,editable:true/false (hl name edition flag)}
	this.highlightedPoints={}; // List of different user-defined labeled set of points
	this.addHighlighting=function(hName,hColor,pointSet,matchPattern) { // for convenience: calls idvis.lib
		idvis.lib.addHighlighting(GP,hName,hColor,pointSet,matchPattern);
	};
	this.multiPointDatasetRule=plotData.multiDatasetRule || plotData.multiPointDatasetRule; // string, updated by idvis if undefined
	this.pointOnClick=plotData.pointOnClick; // used by idvis _displayPointLabel()
	this.pointOnList=plotData.pointOnList;
    //var pointDefaultSize=(plotData.pointDefaultSize)? plotData.pointDefaultSize : 10;
	this.pointOpacity=plotData.pointOpacity || 1; // point opacity
	this.sameScale=(plotData.sameScale)? true : false;
	this.zoomable=plotData.zoomable;
	this.showDatasets=(plotData.hideDatasets)? false : undefined;
    this.searchable=(plotData.noSearch)? false : (plotData.searchable)? plotData.searchable : true;
	this.selectable=(plotData.noSelect)? false : true;
	this.customPointLabel=plotData.customPointLabel || null;
	this.connectPoints=plotData.connectpoint || null; // null or line or curve (global to all datasets or defined for each one)
//this.onPointExclusion=(plotData.onPointExclusion)? plotData.onPointExclusion : null;
	this.convertValue=plotData.convertValue || function(axis,value) {return value;};
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
	for (let i=0; i<axes.length; i++) {
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
    const BORDER_SPACE=10,
		  ZOOM_SPACE=10,
		  XAXIS_SPACE=40,
		  YAXIS_SPACE=50;
    var chartW=plotData.width || 400,
		chartH=(plotData.height)? plotData.height : 400,
		chartX=BORDER_SPACE+1;
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
	this.addThreshold=function() { // ...(th) Obsolete
		alert('WARNING: addThreshold() is obsolete. Use addFeature() instead !');
	};
	this.addFeature=function(fData) { // [axis,name,value,color,dashString,keepAbove1,roundValue]
		var feature=new idvis.feature(GP,fData); // [sh,minX*,maxX*,minY*,maxY*] *Only for custom & path
		if (feature.editable) {
			this.editableThresholds.push(feature);
		}
		//Add to list
		this.features.push(feature);
		//Update chart limit values
		if (feature.type=='line') {
			if (this['minValue'+feature.axis]===null) {this['minValue'+feature.axis]=this['maxValue'+feature.axis]=feature.startValue;}
			else {
				this['minValue'+feature.axis]=Math.min(this['minValue'+feature.axis],feature.startValue);
				this['maxValue'+feature.axis]=Math.max(this['maxValue'+feature.axis],feature.startValue);
			}
		}
		else if (feature.type=='text') {
			let axisY=(feature.axis=='X')? 'Y' : 'Y2';
			if (this['minValue'+feature.axis]===null) {
				this['minValue'+feature.axis]=this['maxValue'+feature.axis]=feature.x;
				this['minValue'+axisY]=this['maxValue'+faxisY]=feature.y;
			}
			else {
				this['minValue'+feature.axis]=Math.min(this['minValue'+feature.axis],feature.x);
				this['maxValue'+axisY]=Math.max(this['maxValue'+axisY],feature.y);
			}
		}
		else if (feature.minX !== null) { // custom, function
			let i0=(feature.axis=='X')? 0 : 2; // ['X','Y','X2','Y2'];
			for (let i=i0; i<=i0+1; i++) {
				let lAxis=(i==i0)? 'X' : 'Y';
				if (this['minValue'+axes[i]]===null) {this['minValue'+axes[i]]=this['maxValue'+axes[i]]=feature['min'+lAxis];}
				else {
					this['minValue'+axes[i]]=Math.min(this['minValue'+axes[i]],feature['min'+lAxis]);
					this['maxValue'+axes[i]]=Math.max(this['maxValue'+axes[i]],feature['max'+lAxis]);
				}
			}
		}
	};

    /***** Datasets *****/
    this.addDataset=function(dsIdx,set) {
		idvis.lib.addDataset(GP,dsIdx,set);
    };
	this.addDataSet=this.addDataset; // just to be safe

    /***** Data points *****/
	this.addDataAsString=function(dsIdx,dataStrg) {
        let minValueXn='minValue'+this.datasets[dsIdx].params.axisX,
            maxValueXn='maxValue'+this.datasets[dsIdx].params.axisX,
            minValueYn='minValue'+this.datasets[dsIdx].params.axisY,
            maxValueYn='maxValue'+this.datasets[dsIdx].params.axisY,
            dataset=dataStrg.split(';');
		for (let i=0; i < dataset.length; i++) {
			let data=dataset[i].split(','),
				dp=new dataPoint(this.datasets[dsIdx],data),
				xProp=(dp.valueList.x)? idvis.valuesProperties(dp.valueList.x) : {min:dp.x,max:dp.x},
				yProp=(dp.valueList.y)? idvis.valuesProperties(dp.valueList.y) : {min:dp.y,max:dp.y};
 			this.datasets[dsIdx].data.push(dp);
			if (this[minValueXn]===null) { // 1st point
				this[minValueXn]=xProp.min;
				this[maxValueXn]=xProp.max;
			}
			else {
				this[minValueXn]=Math.min(this[minValueXn],xProp.min);
				this[maxValueXn]=Math.max(this[maxValueXn],xProp.max);
			}
			if (this[minValueYn]===null) { // 1st point
				this[minValueYn]=yProp.min;
				this[maxValueYn]=yProp.max;
			}
			else {
				this[minValueYn]=Math.min(this[minValueYn],yProp.min);
				this[maxValueYn]=Math.max(this[maxValueYn],yProp.max);
			}
		}
    };

	/********************* Plot display (Raphael) *********************/
    this.draw=function() {
		/* DIVs */
		//this.div.innerHTML='<TABLE><TR><TD style="min-width:'+canvasW+'px"><DIV id="'+canvasDivID+'"></DIV></TD><TD valign="top" style="min-width:200px"><DIV id="'+formDivID+'"></DIV></TD></TR></TABLE>';
		idvis.lib.drawChartLayout(GP,plotData.formPosition,canvasW,canvasDivID,formDivID);
		this.formDiv=document.getElementById(formDivID);
		/* Form menu */
		idvis.lib.initializeForm(GP);
		/* Canvas */
		this.canvas=Raphael(canvasDivID,canvasW,canvasH);
		/* Back panel */
		this.backPanel=this.canvas.rect(0,0,canvasW,canvasH,0).attr({fill:'#fff',stroke:'#000'});
		/* Chart area */
		this.plotArea=this.canvas.rect(chartX,chartY,chartW,chartH,0)
		.attr({fill:'#F3F3F3',stroke:'none'}) //,cursor:"crosshair"
		.drag(function(dx,dy,x,y,e) {idvis.lib.extendDragging(GP,dx,dy,x,y,e);},
			  function(x,y,e) {idvis.lib.setDragging(GP,x,y,e);},
			  function() {idvis.lib.endDragging(GP);}
			  );
		if (this.zoomable) this.plotArea.dblclick(function(){idvis.lib.zoomOut(GP);});


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
		idvis.lib.initializeChart(GP);
        
         /***** Change point opacity according to thresholds *****/
        if (GP.dimArea) {_updatePointOpacity();}

    };

    /* Update points transparency */
    function _updatePointOpacity() {
        let minX=(idvis.isNumber(GP.dimArea.minX))? GP.features[GP.dimArea.minX].value : GP.chartSettings.reference.X.minRange, //GP.minValueX,
            maxX=(idvis.isNumber(GP.dimArea.maxX))? GP.features[GP.dimArea.maxX].value : GP.chartSettings.reference.X.maxRange, //GP.maxValueX,
            minY=(idvis.isNumber(GP.dimArea.minY))? GP.features[GP.dimArea.minY].value : GP.chartSettings.reference.Y.minRange, //GP.minValueY,
            maxY=(idvis.isNumber(GP.dimArea.maxY))? GP.features[GP.dimArea.maxY].value : GP.chartSettings.reference.Y.maxRange; //GP.maxValueY;
//console.log(minX,maxX,minY,maxY);
        let [inOpacity,outOpacity]=(GP.dimArea.rule==='in')? [idvis.dimPointOpacity,GP.pointOpacity] : [GP.pointOpacity,idvis.dimPointOpacity];
        for (let dsIdx=0; dsIdx<GP.datasets.length; dsIdx++) {
            if (GP.datasets[dsIdx].params.axisX==='X2' || GP.datasets[dsIdx].params.axisY==='Y2') continue; // only for primary axes
            for (let i=0; i<GP.datasets[dsIdx].data.length; i++) {
				let dp=GP.datasets[dsIdx].data[i];
                let pOpacity;
                if (GP.dimArea.active===true) {
//console.log(i,dp.x,dp.y);
                    pOpacity=(dp.x > minX && dp.x < maxX && dp.y > minY && dp.y < maxY)? inOpacity : outOpacity;
                }
                else {pOpacity=GP.pointOpacity;}
                //dp.point.attr({'fill-opacity':pOpacity});
                dp.point.animate({'fill-opacity':pOpacity},400,"easeInOut");
            }
        }
    }
    /* Callback after idvis.updateThreshold() => update points transparency */
    this.onThresholdUpdate=_updatePointOpacity;

/*====================================== Nested objects ===================================*/


	/*****************************************************************
					  Data point object
	******************************************************************/
	var dataPoint = function(set,data) { // data=[label,id,<x:values>,<y:values>,point size,]
//alert(data.length);
		this.dataset=set;
		this.label=data[0];
		this.externalID=(data[1])? data[1] : data[0];
		this.x=null;
		this.y=null;
		this.valueList={};
		this.pointList=null; // ref to list svg elements
        this.computePoints(data[2],data[3]);
		//this.computeMean('x',data[2]);
		//this.computeMean('y',data[3]);
		this.size=(data[4])? data[4] : 10;
		this.point=null; // ref to the svg element
		this.highlightNames=[]; // new Array()
	};
	dataPoint.prototype = {     
        computePoints: function(xValuesStrg,yValuesStrg) {
            let xValues=xValuesStrg.split(':'),
                yValues=yValuesStrg.split(':'),
                xProp=idvis.valuesProperties(xValues),
                yProp=idvis.valuesProperties(yValues);
            this.x=xProp.mean; //Math.round(100*xProp.mean)/100;
            this.y=yProp.mean; //Math.round(100*yProp.mean)/100;
            if (xValues.length > 1) {this.valueList.x=xValues;}
            if (yValues.length > 1) {this.valueList.y=yValues;}
        },
 /*       
		computeMean: function(axis,valuesStrg) {
	//console.log(axis+','+valuesStrg);
			let values=valuesStrg.split(':').sort(idvis.sortNumber),
				mean=0;
			for (let i=0; i<values.length; i++) {mean+=(1*values[i]);}
			this[axis]=Math.round(100*mean/values.length)/100;
			if (values.length > 1) {this.valueList[axis]=values;}
		},
		getMinValue: function (axis) {
			return (this.valueList && this.valueList[axis])? this.valueList[axis][0] : this[axis];
		},
		getMaxValue: function (axis) {
			return (this.valueList && this.valueList[axis])? this.valueList[axis][this.valueList[axis].length-1] : this[axis];
		},
*/
		info: function(type) {
			let infoStrg=this.label;
			if (type=='min') {
				return this.label;
			}
			else {
				if (this.dataset.params.chart.datasets.length > 1) infoStrg+='\nSet='+this.dataset.params.name;
				infoStrg+='\n'+this.dataset.params.pointLabelX+'='+(this.x*1).toFixed(3)+'\n'+this.dataset.params.pointLabelY+'='+(this.y*1).toFixed(3);
                if (this.dataset.params.pointSizeName) {infoStrg+='\n'+this.dataset.params.pointSizeName+'='+this.size;}
			}
			return infoStrg;
		},
		getX: function() {
			return this.x;
		},
		getY: function() {
			return this.y;
		}
	};

}; // end of GP



/*
####>Revision history<####
# 2.0.3 [FEATURE] Support for point shape (05/05/21)
# 2.0.2 [FEATURE] Optional external pre-search function (PP 06/12/19)
# 2.0.1 Multi-value points in X and Y axes (PP ../../18)
# 2.0.0 Moved from cubiojs to idvis namespace & updates (PP 14/06/17)
# 1.1.1 Change thresholds to features (../02/14)
# 1.1.0 Upgrade to cubiojs namespace (01/02/14)
# 1.0.1 GPL license (PP 23/09/13)
# 0.0.2 Added point opacity management (PP 15/07/13)
# 0.0.1 Started from volcanoPlot.js v.1.1.0 (PP 17/03/13)
*/